"""Load the verified nominal products used by the Script 07 ensemble."""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
from pathlib import Path

from astropy.table import Table
import numpy as np

from . import bl_population_io, global_disk, rh3_refit_io
from .config import load_config
from .global_disk_io import file_hash, load_json


@dataclass
class NominalRun:
    paths: dict
    provenance: dict
    problem: global_disk.DiskProblem
    fit: dict
    optimizer: dict
    grid_settings: dict
    spectra: dict
    tables: dict
    ids: dict
    setup: dict
    spatial: dict
    bl_settings: dict
    bl_display: dict


def _verified_models(path, table, ids, arm, fingerprint, spectra):
    """Load the winning models by saved bin ID and verify their index and scores."""
    models = np.full_like(spectra['noise'], np.nan, dtype=float)
    components = {name: np.full_like(models, np.nan)
                  for name in (('A', 'B', 'gas') if arm == 'BL' else ('A', 'B'))}
    weights = {}
    digest = hashlib.sha256()
    for bid in ids:
        profile_path = path / f'bin_{bid:04d}.npz'
        digest.update(bytes.fromhex(file_hash(profile_path)))
        with np.load(profile_path, allow_pickle=False) as saved:
            i = int(np.nanargmin(saved['chi2']))
            if (int(saved['bin_id']) != bid or str(saved['fingerprint']) != fingerprint
                    or saved['winning_attempt'][i] < 0
                    or not np.isclose(saved['fraction'][i], table[f'FA_{arm}'][bid], atol=1e-10, rtol=0)):
                raise ValueError(f'Winning {arm} profile disagrees with its nominal table for bin {bid}')
            models[bid] = saved['model_total'][i]
            weights[int(bid)] = saved['weights'][i]
            for key, values in components.items():
                values[bid] = saved[f'model_{key}'][i]
    good = spectra['good'][ids]
    errors = spectra['noise'][ids][good]
    if not np.all(np.isfinite(models[ids][good]) & np.isfinite(errors) & (errors > 0)):
        raise ValueError(f'Invalid nominal model or calibrated {arm} noise on fitted pixels')
    return models, components, weights, digest.hexdigest()


def _read_bl(path, signature, settings, nbin):
    with np.load(path/'products/BL_log_spectra.npz', allow_pickle=False) as data, np.load(
        path/'products/XSL_BL_templates.npz', allow_pickle=False
    ) as basis:
        if str(data['fingerprint']) != signature or str(basis['fingerprint']) != signature:
            raise ValueError('The nominal BL spectra or template fingerprint differs from Script 06')
        spectra = {k: data[k] for k in ('galaxy', 'noise', 'good', 'normalization_scale', 'wavelength')}
        templates = basis['templates']
        labels = basis['component']
        gas_names = basis['gas_names'].tolist()
        available = data['gas_template_available']
        gas = settings['gas']
        setup = dict(matrix=templates, component=labels, wavelength=spectra['wavelength'],
                     n_stellar=int(np.count_nonzero(labels == 0)), gas_names=gas_names,
                     gas_seed=gas['seed'] if gas['enabled'] else [], gas_bounds=gas.get('bounds'),
                     gas_by_bin=available, velscale=float(basis['velscale']), vsyst=float(basis['vsyst']),
                     mdegree=settings['mdegree'], fraction_tolerance=settings['fraction_tolerance'],
                     stellar_normalization=basis['stellar_normalization'],
                     population_metadata={k: basis[k] for k in ('ages', 'metals', 'masses', 'luminosity')})
        if (available.shape != (nbin, len(gas_names)) or templates.shape[1] != len(labels)
                or len(labels) != 2*setup['n_stellar'] + len(gas_names)
                or np.count_nonzero(labels == 1) != setup['n_stellar']):
            raise ValueError('Malformed BL stellar basis or reviewed gas masks')
        pixels = (nbin, len(spectra['wavelength']))
        if any(spectra[k].shape != pixels for k in ('galaxy', 'good', 'noise')):
            raise ValueError('Malformed BL log spectra')
        spectra['good'] = spectra['good'].astype(bool)
        calibrated = data['noise_before_residual_scale'] * data['residual_noise_scale'][:, None]
        if not np.allclose(spectra['noise'][spectra['good']], calibrated[spectra['good']],
                           equal_nan=True, rtol=3e-7, atol=0):
            raise ValueError('Saved BL errors disagree with the nominal residual multipliers')
        display = dict(source_flux_unit=str(data['source_flux_unit']),
                       applied_aperture_scale=float(data['applied_aperture_scale']))
    return setup, spectra, display


def read_nominal(cfg, runs, fit_options):
    r4, r5, r6 = (Path(runs[i]).expanduser().resolve() for i in (4, 5, 6))
    leaf5 = r5/'metadata/script05_manifest.json'
    leaf6 = r6/'metadata/script06_manifest.json'
    m5, m6 = load_json(leaf5), load_json(leaf6)
    c5, c6 = m5['contract'], m6['contract']
    if c6['script05'] != m5['fingerprint'] or c6['model'] != c5['model']:
        raise ValueError('Choose matching nominal Script 05 and 06 runs')
    if str(cfg.TARGET_NAME) != c5['target'] or not np.isclose(cfg.REDSHIFT, c5['redshift'], atol=1e-10, rtol=0):
        raise ValueError('The Script 07 target/redshift differs from the nominal analysis')
    saved5 = load_config(r5/'metadata/config_snapshot.py', validate=False)
    verified = rh3_refit_io.read_inputs(saved5, r4, runs)
    if any(c5['inputs'][key] != digest for key, digest in verified['contract']['inputs'].items()):
        raise ValueError('Supplied Script 01–04 products differ from the Script 05 handoff')
    fit = verified['fit']
    if c5['model'] != fit['name']:
        raise ValueError('Selected Script 04 model differs from the Script 05 handoff')
    paths = verified['paths'].copy()
    table_paths = dict(RH3=r5/'products/RH3_global_kinematics.ecsv',
                       BL=r6/'products/BL_stellar_populations.ecsv')
    tables = {arm: Table.read(p, format='ascii.ecsv') for arm, p in table_paths.items()}
    nbin = len(verified['coarse_min'])
    for arm, manifest in (('RH3', m5), ('BL', m6)):
        table = tables[arm]
        if (table.meta['fingerprint'] != manifest['fingerprint'] or table.meta['model'] != fit['name']
                or not np.array_equal(table['BIN_ID'], np.arange(nbin))
                or np.count_nonzero(table['USABLE']) != manifest['n_usable_bins']):
            raise ValueError(f'Nominal {arm} table does not agree with the saved run')
    if c6['hashes']['handoff'] != file_hash(table_paths['RH3']):
        raise ValueError('BL analysis used a different Script 05 table')
    rh_ids = np.flatnonzero(tables['RH3']['USABLE'])
    bl_ids = np.flatnonzero(tables['BL']['USABLE'] & tables['BL']['NOISE_CALIBRATED'] & tables['RH3']['USABLE'])
    if not len(rh_ids) or not len(bl_ids) or not np.isin(rh_ids, verified['ids']).all():
        raise ValueError('No usable paired nominal experiment or an invalid RH3 cohort')
    if not np.allclose(
        np.column_stack([tables['BL'][k][bl_ids] for k in bl_population_io.STELLAR_COLUMNS]),
        np.column_stack([tables['RH3'][k][bl_ids] for k in bl_population_io.STELLAR_COLUMNS]),
        atol=1e-8, rtol=0,
    ):
        raise ValueError('Nominal BL stellar kinematics disagree with the RH3 handoff')

    rh3 = verified['spectra']
    bl_setup, bl, bl_display = _read_bl(r6, m6['fingerprint'], c6['settings'], nbin)
    model_hash = {}
    for arm, path, spectrum, cohort, manifest in (
        ('RH3', r5, rh3, rh_ids, m5), ('BL', r6, bl, bl_ids, m6)
    ):
        models, components, weights, digest = _verified_models(
            path/'products'/f'{arm}_fraction_profiles', tables[arm], cohort, arm,
            manifest['fingerprint'], spectrum)
        if not np.isfinite(spectrum['normalization_scale'][cohort]).all() or np.any(
            spectrum['normalization_scale'][cohort] <= 0
        ):
            raise ValueError(f'Invalid nominal {arm} spectral normalization')
        spectrum.update(model=models, components=components, model_weights=weights)
        model_hash[arm] = digest

    paths.update(nominal_rh3=table_paths['RH3'], nominal_bl=table_paths['BL'],
                 bl_spectra=r6/'products/BL_log_spectra.npz', bl_templates=r6/'products/XSL_BL_templates.npz',
                 bl_review=r6/'metadata/BL_line_approval.json', manifest5=leaf5, manifest6=leaf6,
                 config4=r4/'metadata/config_snapshot.py', psf=r4/'products/adopted_RH3_PSF.npz')
    review = load_json(paths['bl_review'])
    if review['review_fingerprint'] != m6['fingerprint']:
        raise ValueError('BL line review differs from the approved nominal fit')
    disk_contract = load_json(paths['model_contract'])['contract']
    disk_source = disk_contract['provenance']
    disk_config = load_config(paths['config4'], validate=False)
    optimizer = fit_options(disk_config)
    with np.load(paths['psf'], allow_pickle=False) as psf:
        kernel, light, east, north = (psf[key] for key in ('kernel', 'light', 'east', 'north'))
    if (hashlib.sha256(kernel.tobytes()).hexdigest() != disk_contract['kernel_sha256']
            or hashlib.sha256(light.tobytes()).hexdigest() != disk_source['continuum_sha256']
            or not np.array_equal(east, verified['east']) or not np.array_equal(north, verified['north'])):
        raise ValueError('Selected Script 04 PSF and continuum products fail their provenance checks')
    s1, s2 = paths['coordinates'].parent.parent, paths['membership'].parent.parent
    paths.update(bl_coords=s1/'products/BL_spatial_coordinates.npz',
                 bin_table=s2/'products/master_bin_table.ecsv')
    if file_hash(paths['bl_coords']) != c6['hashes']['coordinates']:
        raise ValueError('The BL map coordinates differ from Script 06')
    with np.load(paths['bl_coords'], allow_pickle=False) as grids:
        bl_east, bl_north = (grids[key] for key in ('x_arcsec', 'y_arcsec'))
    with np.load(paths['membership'], allow_pickle=False) as membership:
        ry, rx, rb = (membership[f'rh3_{key}'] for key in ('y_pix', 'x_pix', 'bin_id'))
        by, bx, bb = (membership[f'bl_{key}'] for key in ('y_pix', 'x_pix', 'bin_id'))
    bl_map = np.full(bl_east.shape, -1, dtype=int)
    bl_map[by, bx] = bb
    nodes = np.asarray(fit['nodes'], dtype=float)
    geometry = np.asarray(disk_source['geometry_initial'], dtype=float)
    vmax = float(disk_config.SCRIPT04_VROT_MAX_KMS)
    bounds = np.vstack((np.asarray(disk_source['geometry_bounds']), np.tile([0., vmax], (2*len(nodes), 1))))
    psf_enabled = fit.get('fit_group', {}).get('psf_enabled', '_psf1_' in fit['name'])
    problem = global_disk.DiskProblem(
        east, north, light, ry, rx, rb, rh_ids, np.arange(nbin), None,
        *verified['axes'], nodes, bounds, geometry, min(float(fit['info']['initial_speed_kms']), vmax),
        int(fit['sign_a']), kernel if psf_enabled else None)
    c3 = load_config(paths['script3_config'], validate=False)
    provenance = dict(target=c5['target'], redshift=c5['redshift'], model=fit['name'],
                      source_hashes={key: file_hash(path) for key, path in paths.items()},
                      winning_profile_hashes=model_hash, rh3_eligible=rh_ids.tolist(),
                      bl_eligible=bl_ids.tolist(), coefficient=fit['coefficient'], optimizer=optimizer,
                      nominal_optimizer_has_stability_record='kinematic_refinement' in fit,
                      bl_review=review['review'])
    return NominalRun(
        dict(paths, run4=r4, run5=r5, run6=r6, run1=s1, run2=s2), provenance, problem, fit, optimizer,
        dict(sigma_start=float(c3.RH3_SIGMA_START_KMS)),
        dict(RH3=rh3, BL=bl), tables, dict(RH3=rh_ids, BL=bl_ids),
        dict(RH3=verified['setup'], BL=bl_setup),
        dict(RH3=(east, north, verified['bin_map']), BL=(bl_east, bl_north, bl_map)),
        c6['settings'], bl_display)
