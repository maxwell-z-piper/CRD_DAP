"""Read the Script 05 handoff and build the approved BL spectral experiment."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

from astropy.io import fits
from astropy.table import Table
import numpy as np

from . import populations, templates
from .global_disk_io import file_hash, load_json

# Vacuum rest wavelengths from MaNGA DR17. The O III pair uses the DR17
# 4959:5007 flux tie. Each other listed line has a separate amplitude.
CORE_LINES = {
    'OII3726': (3727.092, 1.), 'OII3729': (3729.875, 1.),
    'Hdelta': (4102.899, 1.), 'Hgamma': (4341.691, 1.), 'Hbeta': (4862.691, 1.),
    'OIII4959': (4960.295, .35), 'OIII5007': (5008.240, 1.),
    'NI5198': (5199.349, 1.), 'NI5200': (5201.705, 1.),
}
OPTIONAL_LINES = {
    'H10': (3798.983, 1.), 'H9': (3836.479, 1.), 'NeIII3869': (3869.860, 1.),
    'H8': (3890.158, 1.), 'NeIII3967': (3968.590, 1.),
    'Hepsilon': (3971.202, 1.), 'HeII4686': (4687.015, 1.),
}
LINES = CORE_LINES | OPTIONAL_LINES
C_KMS = templates.C_KMS
STELLAR_COLUMNS = ('VA_GLOBAL_KMS', 'SIGMA_A_KMS', 'VB_GLOBAL_KMS', 'SIGMA_B_KMS')


def fingerprint(data):
    return hashlib.sha256(json.dumps(data, sort_keys=True, allow_nan=False).encode()).hexdigest()


def interval(values, name):
    values = tuple(map(float, values))
    if len(values) != 2 or not np.isfinite(values).all() or values[0] >= values[1]:
        raise ValueError(f'{name} must contain two increasing finite values.')
    return values


def gas_settings(cfg):
    enabled = bool(cfg.FIT_BL_GAS_EMISSION)
    if not enabled:
        return dict(enabled=False, seed=None, source='disabled'), None
    manual = [getattr(cfg, name, None)
              for name in ('MANGA_GAS_V_SEED_KMS', 'MANGA_GAS_SIGMA_SEED_KMS')]
    path = None
    source = dict(overrides=manual)
    if all(value is not None for value in manual):
        seed = list(map(float, manual))
    else:
        catalog = getattr(cfg, 'MANGA_DAPALL', None)
        if catalog is None:
            raise ValueError('Set MANGA_DAPALL or supply both manual gas seeds.')
        path = Path(catalog).expanduser().resolve()
        plateifu = getattr(cfg, 'MANGA_PLATEIFU', None) or str(cfg.TARGET_NAME)
        daptype = str(cfg.MANGA_DAPTYPE)
        with fits.open(path) as hdus:
            data = hdus[daptype].data
            rows = data[np.char.strip(data['PLATEIFU'].astype(str)) == plateifu]
            if len(rows) != 1 or not rows[0]['DAPDONE'] or int(rows[0]['DAPQUAL']) & (1 << 30):
                raise ValueError(f'No unique successful DAPall row for {plateifu}, {daptype}.')
            redshift, dispersion = map(float, (rows[0]['HA_Z'], rows[0]['HA_GSIGMA_1RE']))
        velocity = C_KMS*np.log((1+redshift)/(1+float(cfg.REDSHIFT))) if redshift > -1 else np.nan
        seed = [float(value) if value is not None else float(catalog_value)
                for value, catalog_value in zip(manual, (velocity, dispersion))]
        source.update(plateifu=plateifu, daptype=daptype,
                      ha_z=redshift if np.isfinite(redshift) else None,
                      ha_gsigma_1re_kms=dispersion if np.isfinite(dispersion) else None)
    if not np.isfinite(seed).all() or seed[1] <= 0:
        raise ValueError('Invalid DAPall gas seed. Supply a finite velocity and positive corrected dispersion.')
    half = float(getattr(cfg, 'BL_GAS_V_HALF_RANGE_KMS', 2000.))
    floor = float(getattr(cfg, 'BL_GAS_SIGMA_MIN_PIXELS', .01))
    upper = float(getattr(cfg, 'BL_GAS_SIGMA_MAX_KMS', 1000.))
    if not np.isfinite([half, floor, upper]).all() or min(half, floor, upper) <= 0:
        raise ValueError('BL gas bounds must be finite and positive.')
    return dict(enabled=True, seed=seed, source=source, half_range=half,
                sigma_floor_pixels=floor, sigma_max_kms=upper), path


def _upstream_paths(run5, run1, run2):
    run5 = Path(run5).expanduser().resolve()
    manifest_path = run5 / 'metadata/script05_manifest.json'
    manifest = load_json(manifest_path)
    saved = manifest['input_paths']
    run1 = Path(run1 or Path(saved['coordinates']).parent.parent).expanduser().resolve()
    run2 = Path(run2 or Path(saved['membership']).parent.parent).expanduser().resolve()
    return manifest, dict(
        handoff_manifest=manifest_path, handoff=run5/'products/RH3_global_kinematics.ecsv',
        membership=run2/'products/master_bin_membership.npz',
        bins=run2/'products/master_bins.fits', bin_table=run2/'products/master_bin_table.ecsv',
        bl_spectra=run2/'products/master_bin_spectra.fits',
        script2_manifest=run2/'metadata/script02_manifest.json',
        prepared=run1/'products/prepared_BL.fits',
        coordinates=run1/'products/BL_spatial_coordinates.npz',
        lsf=templates.discover_script1_lsf_product(run1, arm='BL'),
    )


def read_inputs(cfg, run5, run1, run2):
    manifest, paths = _upstream_paths(run5, run1, run2)
    old = manifest['contract']
    target, z = str(cfg.TARGET_NAME), float(cfg.REDSHIFT)
    if old['target'] != target or not np.isclose(old['redshift'], z, atol=1e-10, rtol=0):
        raise ValueError('Target or redshift differs from the Script 05 handoff.')
    table = Table.read(paths['handoff'], format='ascii.ecsv')
    bins = Table.read(paths['bin_table'], format='ascii.ecsv')
    ids = np.asarray(table['BIN_ID'])
    if (table.meta['fingerprint'] != manifest['fingerprint'] or table.meta['model'] != old['model']
            or load_json(paths['script2_manifest'])['target'] != target
            or not np.array_equal(ids, np.arange(len(table)))
            or not np.array_equal(ids, bins['BIN_ID'])
            or int(sum(table['SELECTED'])) != manifest['n_selected_bins']
            or int(sum(table['USABLE'])) != manifest['n_usable_bins']
            or np.any(table['USABLE'] & ~table['SELECTED'])):
        raise ValueError('Script 05 and Script 02 disagree on the model or bin domain.')
    if file_hash(paths['membership']) != old['inputs']['membership']:
        raise ValueError('Script 02 membership differs from the one used in Script 05.')
    eligible = np.asarray(table['SELECTED'] & table['USABLE'], dtype=bool)
    kin = np.column_stack([np.asarray(table[key], dtype=float) for key in STELLAR_COLUMNS])
    if (not eligible.any() or not np.isfinite(kin[eligible]).all()
            or np.any(kin[eligible][:, [1, 3]] <= 0)):
        raise ValueError('Script 05 supplies no usable stellar kinematics or contains invalid values.')

    with fits.open(paths['bl_spectra'], memmap=False) as h:
        wave, flux, noise = (np.asarray(h[f'BL_{key}'].data, dtype=float)
                             for key in ('WAVE', 'FLUX', 'UNCERT'))
        flags = np.asarray(h['BL_GOOD'].data)
        source_flux_unit = str(h['BL_FLUX'].header.get('SRCBUNIT', ''))
        aperture_scale = float(h[0].header.get('BLSCALE', 1.))
    if (wave.ndim != 1 or not np.isfinite(wave).all() or np.any(np.diff(wave) <= 0)
            or not (flux.shape == noise.shape == flags.shape == (len(table), len(wave)))
            or not np.isin(flags, (0, 1)).all()):
        raise ValueError('Malformed BL wave, flux, uncertainty or mask from Script 02.')
    good = flags.astype(bool)
    if not np.all(np.isfinite(flux[good]) & np.isfinite(noise[good]) & (noise[good] > 0)):
        raise ValueError('Invalid Script 02 BL flux or calibrated error on good pixels.')
    with fits.open(paths['prepared']) as h:
        header = h[0].header
        prepared_wave = h['WAVELENGTH'].data
        if not np.array_equal(prepared_wave, wave):
            raise ValueError('Script 01 and Script 02 use different observed BL grids.')
        science_medium = templates.infer_science_wavelength_medium(header)
    if cfg.SCIENCE_WAVELENGTH_MEDIUM not in ('auto', science_medium):
        raise ValueError('Configured science medium disagrees with prepared BL metadata.')
    template_medium = str(cfg.TEMPLATE_WAVELENGTH_MEDIUM).lower()
    rest_wave = templates.convert_wavelength_medium(wave/(1+z), science_medium, template_medium)
    velscale = C_KMS*float(np.median(np.diff(np.log(rest_wave))))
    lsf = templates.load_saved_lsf_model(paths['lsf'], polynomial_order=int(cfg.LSF_MODEL_WAVELENGTH_ORDER))
    lsf_range = templates.convert_wavelength_medium(
        np.array([lsf.empirical_min, lsf.empirical_max])/(1+z), science_medium, template_medium)

    with np.load(paths['coordinates'], allow_pickle=False) as spatial:
        east, north = spatial['x_arcsec'], spatial['y_arcsec']
    with np.load(paths['membership'], allow_pickle=False) as membership:
        y, x, bid = (membership[f'bl_{key}'] for key in ('y_pix', 'x_pix', 'bin_id'))
    bin_map = fits.getdata(paths['bins']).astype(int)
    if (east.shape != north.shape or east.shape != bin_map.shape
            or not np.array_equal(np.unique(bin_map[bin_map >= 0]), ids)
            or not np.array_equal(bin_map[y, x], bid)
            or np.count_nonzero(bin_map >= 0) != len(bid)):
        raise ValueError('BL coordinates, bin map and saved membership disagree.')

    gas, dapall = gas_settings(cfg)
    if dapall is not None:
        paths['dapall'] = dapall
    paths['xsl'] = Path(cfg.XSL_TEMPLATE_LIBRARY).expanduser().resolve()
    fit_range = interval(cfg.BL_FIT_REST_RANGE_ANGSTROM, 'BL fitting band')
    reference = interval(cfg.BL_FRACTION_REST_RANGE_ANGSTROM, 'BL fraction reference band')
    mdegree = int(cfg.BL_MDEGREE)
    if int(cfg.BL_DEGREE) != -1 or float(cfg.BL_REGUL_LIKELIHOOD) != 0 or mdegree < 0:
        raise ValueError('Script 06 requires no additive polynomial and an unregularized BL fit.')
    if gas['enabled']:
        v, sigma = gas['seed']
        half = gas['half_range']
        gas['bounds'] = [[v-half, v+half], [gas['sigma_floor_pixels']*velscale, gas['sigma_max_kms']]]
        if not gas['bounds'][1][0] <= sigma <= gas['bounds'][1][1]:
            raise ValueError('Gas dispersion seed falls outside its configured bounds.')
    rebin_fraction = float(cfg.BL_LOG_REBIN_MIN_VALID_FRACTION)
    if not 0 < rebin_fraction <= 1:
        raise ValueError('BL_LOG_REBIN_MIN_VALID_FRACTION must be in (0, 1].')
    settings = dict(fit_range=fit_range, reference=reference, mdegree=mdegree,
                    velscale=velscale, science_medium=science_medium, template_medium=template_medium,
                    gas=gas, rebin_fraction=rebin_fraction, ml_band=str(cfg.BL_MASS_TO_LIGHT_BAND),
                    lsf_order=int(cfg.LSF_MODEL_WAVELENGTH_ORDER),
                    fraction_tolerance=float(getattr(cfg, 'BL_FRACTION_CONSTRAINT_TOLERANCE', 1e-5)),
                    imf='Salpeter (configured XSL library)')
    if not np.isfinite(settings['fraction_tolerance']) or settings['fraction_tolerance'] <= 0:
        raise ValueError('BL fraction constraint tolerance must be finite and positive.')
    hashes = {key: file_hash(path) for key, path in paths.items()
              if key not in ('prepared', 'membership')}
    hashes['membership'] = old['inputs']['membership']
    hashes['prepared_bl_header'] = hashlib.sha256(header.tostring().encode()).hexdigest()
    contract = dict(target=target, redshift=z, model=old['model'], script05=manifest['fingerprint'],
                    hashes=hashes, settings=settings)
    return dict(paths=paths, contract=contract, table=table, eligible=eligible, kin=kin,
                wave=wave, rest_wave=rest_wave, flux=flux, noise=noise, good=good,
                source_flux_unit=source_flux_unit, aperture_scale=aperture_scale,
                lsf=lsf, lsf_range=lsf_range, east=east, north=north,
                bin_map=bin_map, bin_table=bins)


def default_review(cfg, inputs):
    return dict(
        lines=list(CORE_LINES) + list(getattr(cfg, 'BL_OPTIONAL_GAS_LINES', ())),
        masks=[
            dict(frame=frame, interval=list(interval(span, f'BL {frame} mask')), bins=None,
                 reason='Target config')
            for frame, key in (('observed', 'BL_MASK_OBSERVED_RANGES_ANGSTROM'),
                               ('rest', 'BL_MASK_REST_RANGES_ANGSTROM'))
            for span in getattr(cfg, key, ())
        ],
    )


def validate_review(review, n_bins):
    lines = review['lines']
    if (not isinstance(lines, list) or len(set(lines)) != len(lines)
            or not set(lines).issubset(LINES)):
        raise ValueError('Review lines must be unique names from the proposed BL line catalog.')
    if ('OIII4959' in lines) != ('OIII5007' in lines):
        raise ValueError('Approve or remove the tied OIII4959 and OIII5007 pair together.')
    for entry in review['masks']:
        if entry['frame'] not in ('observed', 'rest') or not entry['reason'].strip():
            raise ValueError('Each reviewed mask requires a wavelength frame and a reason.')
        interval(entry['interval'], 'Reviewed exclusion')
        bins = entry['bins']
        if bins is not None and (not bins or any(type(b) is not int or not 0 <= b < n_bins for b in bins)):
            raise ValueError('Mask bin IDs must refer to original Script 02 bins.')


def _target_fwhm(inputs, wave):
    settings = inputs['contract']['settings']
    return templates.observed_lsf_to_template_rest_fwhm(
        wave, lsf_model=inputs['lsf'], redshift=inputs['contract']['redshift'],
        science_medium=settings['science_medium'], template_medium=settings['template_medium'],
    )


def build_experiment(inputs, review):
    """Freeze the spectra, masks, line support and LSF after the reviewed exclusions."""
    validate_review(review, len(inputs['eligible']))
    settings = inputs['contract']['settings']
    z = inputs['contract']['redshift']
    vel = settings['velscale']
    kin = inputs['kin'][inputs['eligible']]
    # Allow six Gaussian sigmas of physical stellar support. This padding uses
    # the fixed RH3 widths, not the much wider gas search. Gas has its own grid padding.
    low_shift = float(np.min(kin[:, [0, 2]]-6*kin[:, [1, 3]]))
    high_shift = float(np.max(kin[:, [0, 2]]+6*kin[:, [1, 3]]))
    source_lo, source_hi = inputs['lsf_range'] + [6., -6.]
    # Preserve one extra template sample at either end for the log-grid edges.
    blue_padding = np.exp(max(0., high_shift+vel)/C_KMS)
    red_padding = np.exp(min(0., low_shift-vel)/C_KMS)
    safe_lo = source_lo*blue_padding
    safe_hi = source_hi*red_padding
    fit_lo = max(settings['fit_range'][0], safe_lo, inputs['rest_wave'][0])
    fit_hi = min(settings['fit_range'][1], safe_hi, inputs['rest_wave'][-1])
    ref_lo, ref_hi = settings['reference']
    if not fit_lo <= ref_lo < ref_hi <= fit_hi:
        raise ValueError('The fixed BL reference band is outside the supported stellar fitting band.')
    wave, edges = templates.make_log_wavelength_grid(fit_lo, fit_hi, velscale=vel)
    native_good = inputs['good'].copy()
    approved_log_good = np.ones((len(inputs['eligible']), len(wave)), dtype=bool)
    observed_log_wave = templates.convert_wavelength_medium(
        wave, settings['template_medium'], settings['science_medium'])*(1+z)
    for mask in review['masks']:
        observed = mask['frame'] == 'observed'
        native_grid = inputs['wave'] if observed else inputs['rest_wave']
        log_grid = observed_log_wave if observed else wave
        lo, hi = mask['interval']
        for grid, valid in ((native_grid, native_good), (log_grid, approved_log_good)):
            selected = (grid >= lo) & (grid <= hi)
            if mask['bins'] is None:
                valid[:, selected] = False
            else:
                valid[np.ix_(mask['bins'], selected)] = False
    native_good[:, ~np.isfinite(inputs['lsf'].evaluate(inputs['wave']))] = False
    shape = (len(inputs['eligible']), len(wave))
    galaxy = np.full(shape, np.nan)
    noise = np.full(shape, np.nan)
    good = np.zeros(shape, dtype=bool)
    scales = np.full(shape[0], np.nan)
    for bid in np.flatnonzero(inputs['eligible']):
        galaxy[bid], noise[bid], good[bid], _ = templates.rebin_spectrum_with_diagonal_noise(
            wavelength=inputs['rest_wave'], flux=inputs['flux'][bid],
            uncertainty=inputs['noise'][bid], good=native_good[bid],
            out_wavelength=wave, out_edges=edges, min_valid_fraction=settings['rebin_fraction'])
        good[bid] &= approved_log_good[bid]
        valid = good[bid]
        if valid.any():
            amplitude = np.abs(galaxy[bid, valid])
            scale = float(np.median(amplitude)) or float(np.percentile(amplitude, 75))
            if np.isfinite(scale) and scale > 0:
                scales[bid] = scale
                galaxy[bid] /= scale
                noise[bid] /= scale

    line_vac = np.array([val[0] for val in LINES.values()])
    line_wave = templates.convert_wavelength_medium(line_vac, 'vacuum', settings['template_medium'])
    line_obs = templates.convert_wavelength_medium(line_vac, 'vacuum', settings['science_medium'])*(1+z)
    gas = settings['gas']
    v_seed, s_seed = gas['seed'] if gas['enabled'] else (0., 0.)
    v_lo, v_hi = gas['bounds'][0] if gas['enabled'] else (0., 0.)
    fwhm = _target_fwhm(inputs, line_wave)
    center = line_wave*np.exp(v_seed/C_KMS)
    # Use the full FWHM at the nominal seed to classify coverage. Mark partial
    # if at least half this core and two pixels are available. Show the full
    # allowed center excursion separately for review.
    core_width = np.hypot(fwhm, 2.354820045*center*s_seed/C_KMS)
    core = np.column_stack((center-core_width/2, center+core_width/2))
    search = line_wave[:, None]*np.exp(np.array([v_lo, v_hi])/C_KMS)
    coverage = np.zeros((len(galaxy), len(LINES)))
    counts = np.zeros_like(coverage, dtype=int)
    for i, (lo, hi) in enumerate(core):
        if np.isfinite([lo, hi]).all() and hi > lo:
            widths = np.maximum(0., np.minimum(edges[1:], hi)-np.maximum(edges[:-1], lo))
            coverage[:, i] = np.clip((good @ widths)/(hi-lo), 0., 1.)
            counts[:, i] = (good & (widths > 0)).sum(axis=1)
    usable = ((coverage >= .5) & (counts >= 2) & np.isfinite(fwhm)
              & np.array([gas['enabled'] and name in review['lines'] for name in LINES]))
    source_range = (max(source_lo, fit_lo/blue_padding), min(source_hi, fit_hi/red_padding))
    data = dict(wavelength=wave, edges=edges, galaxy=galaxy, noise=noise, good=good,
                normalization=scales, line_wave=line_wave, line_obs=line_obs, core=core,
                search=search, coverage=coverage, usable=usable, fwhm=fwhm,
                effective_range=(fit_lo, fit_hi), template_range=source_range,
                native_good=native_good)
    return data


def build_templates(inputs, experiment):
    from ppxf.ppxf_util import gaussian

    settings = inputs['contract']['settings']
    gas = settings['gas']
    velscale = settings['velscale']
    wave = experiment['wavelength']
    basis = templates.prepare_xsl_bl_templates(
        xsl_path=inputs['paths']['xsl'], supported_range=experiment['template_range'],
        velscale=velscale, target_fwhm_rest=lambda w: _target_fwhm(inputs, w),
        template_medium=settings['template_medium'], normalization_range=settings['reference'],
    )
    meta = populations.load_ssp_metadata(inputs['paths']['xsl'], settings['template_medium'],
                                         settings['ml_band'])
    n_stellar = len(meta['ages'])*len(meta['metals'])
    if basis.templates.shape[1] != n_stellar:
        raise ValueError('BL stellar templates and population metadata have different grid sizes.')
    max_shift = max(map(abs, gas['bounds'][0])) + 6*gas['bounds'][1][1] if gas['enabled'] else 0.
    pad = int(np.ceil(max_shift/velscale))
    offsets = np.arange(-pad, len(basis.wavelength)+pad)
    template_wave = basis.wavelength[0]*np.exp(offsets*velscale/C_KMS)
    stellar = np.pad(basis.templates, ((pad, pad), (0, 0)))
    individual = experiment['usable']
    groups, columns, eligible = [], [], []
    for i, name in enumerate(LINES):
        if name == 'OIII4959':
            continue
        indices = [i]
        if name == 'OIII5007':
            indices.insert(0, list(LINES).index('OIII4959'))
        valid = [j for j in indices if individual[:, j].any()]
        if not valid:
            continue
        # Include both members of the tied doublet when their instrumental
        # profiles exist. The approved pixel mask still governs their leverage.
        model = [j for j in indices if np.isfinite(experiment['fwhm'][j])]
        centers = experiment['line_wave'][model]
        widths = experiment['fwhm'][model]
        ratios = np.array([list(LINES.values())[j][1] for j in model])
        profiles = gaussian(np.log(template_wave), centers, widths)
        profiles /= (np.gradient(template_wave) @ profiles)[None, :]
        columns.append(profiles @ ratios)
        eligible.append(individual[:, valid].any(axis=1))
        groups.append('OIII4959_5007' if name == 'OIII5007' else name)
    gas_matrix = np.column_stack(columns) if columns else np.empty((len(template_wave), 0))
    usable = np.column_stack(eligible) if eligible else np.empty((len(inputs['eligible']), 0), dtype=bool)
    setup = dict(
        matrix=np.column_stack((stellar, stellar, gas_matrix)),
        component=np.repeat([0, 1, 2], [n_stellar, n_stellar, len(groups)]),
        template_wavelength=template_wave, wavelength=wave,
        vsyst=C_KMS*float(np.log(template_wave[0]/wave[0])),
        velscale=velscale, n_stellar=n_stellar, gas_names=groups,
        gas_seed=gas['seed'] if gas['enabled'] else [],
        gas_bounds=gas.get('bounds'), mdegree=settings['mdegree'],
        fraction_tolerance=settings['fraction_tolerance'],
        stellar_normalization=basis.normalization, population_metadata=meta,
    )
    return setup, usable
