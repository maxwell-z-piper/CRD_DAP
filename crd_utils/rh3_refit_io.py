"""Verify the saved Script 03 and Script 04 handoff for exact RH3 extraction."""
from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
from astropy.table import Table

from .config import load_config
from .global_disk_io import file_hash, load_json
from .rh3_refit import FRACTIONS


def read_inputs(cfg, script4_run, run_overrides):
    run4 = Path(script4_run).expanduser().resolve()
    paths = {
        'selection': run4 / 'products/selected_global_model.json',
        'model_contract': run4 / 'metadata/script04_inputs.json',
    }
    selection = load_json(paths['selection'])
    contract = load_json(paths['model_contract'])
    fingerprint = contract['fingerprint']
    if selection['fingerprint'] != fingerprint:
        raise ValueError('The selected model does not match the Script 04 input fingerprint')
    name = selection['model']
    if not name or not all(c.isalnum() or c == '_' for c in name):
        raise ValueError('Invalid selected Script 04 model name')
    # Use the run-local filenames so Windows products can be moved to a Mac.
    paths.update(fit=run4 / f'fits/{name}.json',
                 predictions=run4 / f'products/{name}_bin_predictions.ecsv')
    fit = load_json(paths['fit'])
    predictions = Table.read(paths['predictions'], format='ascii.ecsv')
    if (fit['fingerprint'] != fingerprint or predictions.meta.get('fingerprint') != fingerprint
            or fit['name'] != name or predictions.meta.get('model') != name):
        raise ValueError('Script 04 selection, fit and prediction table disagree')

    provenance = contract['contract']['provenance']
    runs = {
        n: Path(run_overrides.get(n) or provenance[f'source_script{n}_run']).expanduser().resolve()
        for n in (1, 2, 3)
    }
    paths.update(
        membership=runs[2] / 'products/master_bin_membership.npz',
        coordinates=runs[1] / 'products/RH3_spatial_coordinates.npz',
        cube=runs[3] / 'products/RH3_likelihood_cubes.npz',
        spectra=runs[3] / 'products/RH3_log_spectra_and_local_best_fits.npz',
        templates=runs[3] / 'products/XSL_RH3_templates.npz',
        script3_manifest=runs[3] / 'metadata/script03_manifest.json',
        script3_config=runs[3] / 'metadata/config_snapshot.py',
    )
    hashes = {name: file_hash(path) for name, path in paths.items()}
    for key, saved_key in (('cube', 'likelihood_sha256'), ('membership', 'membership_sha256')):
        if hashes[key] != provenance[saved_key]:
            raise ValueError(f'The supplied {key} differs from the product used by Script 04')

    manifest = load_json(paths['script3_manifest'])
    if str(manifest['target']) != str(cfg.TARGET_NAME) or not np.isclose(
        float(manifest['redshift']), float(cfg.REDSHIFT), atol=1e-10, rtol=0
    ):
        raise ValueError('Target or redshift differs from the Script 03 run')
    if not manifest.get('likelihood_width_scalar_calibrated'):
        raise ValueError('Script 05 requires the scalar-calibrated Script 03 likelihood and noise products')
    recorded_cfg = load_config(paths['script3_config'], validate=False)
    nbin = int(manifest['n_bins'])

    with np.load(paths['spectra'], allow_pickle=False) as z:
        spectra = {key: z[key] for key in ('wavelength', 'galaxy', 'noise', 'good', 'normalization_scale')}
        noise_hash = str(z['noise_model_hash'])
        noise_before = z['noise_before_residual_scale']
        residual_scale = z['residual_noise_scale']
    if noise_hash != manifest['noise_model']['residual_calibration_hash']:
        raise ValueError('Saved spectra and Script 03 noise manifest disagree')
    with np.load(paths['cube'], allow_pickle=False) as z:
        coarse_min = z['local_best_chi2_total']
        axes = [z[key] for key in ('VA_grid', 'VB_grid', 'fA_grid')]
        if str(z['noise_model_hash']) != noise_hash:
            raise ValueError('Script 03 likelihood and spectra use different noise experiments')
    if coarse_min.shape != (nbin,) or not np.all(np.isfinite(coarse_min)):
        raise ValueError('Missing finite Script 03 minima for the saved bins')
    for axis in axes:
        if axis.ndim != 1 or len(axis) < 2 or not np.all(np.isfinite(axis)) or not np.all(np.diff(axis) > 0):
            raise ValueError('Malformed Script 03 grid axis')

    wavelength = spectra['wavelength']
    with np.load(paths['templates'], allow_pickle=False) as z:
        basis = np.asarray(z['templates'], dtype=float)
        template_wave = z['wavelength']
        template_velscale = float(z['velscale_kms'])
        normalization_range = z['normalization_range'].tolist()
    for wave in (wavelength, template_wave):
        if wave.ndim != 1 or len(wave) < 5 or not np.all(np.isfinite(wave) & (wave > 0)) or not np.all(np.diff(wave) > 0):
            raise ValueError('Saved wavelength axes must be positive and strictly increasing')
    if basis.ndim != 2 or basis.shape[0] != len(template_wave) or basis.shape[1] < 1 or not np.all(np.isfinite(basis)):
        raise ValueError('Malformed saved RH3 template basis')
    velscale = float(manifest['velscale_kms'])
    if (not np.isfinite(velscale) or velscale <= 0
            or not np.isclose(template_velscale, velscale, rtol=1e-9, atol=0)
            or any(not np.allclose(np.diff(np.log(wave))*299792.458, velscale, rtol=1e-7, atol=1e-6)
                   for wave in (wavelength, template_wave))):
        raise ValueError('Saved spectra and templates have inconsistent velocity sampling')

    for key in ('galaxy', 'noise', 'good'):
        if spectra[key].shape != (nbin, len(wavelength)):
            raise ValueError(f'Malformed Script 03 {key} array')
    if not np.all(np.isin(spectra['good'], (0, 1))):
        raise ValueError('Saved good-pixel mask is not binary')
    good = spectra['good'].astype(bool)
    if np.any(good.sum(axis=1) < 5) or not np.all(
        np.isfinite(spectra['galaxy'][good]) & np.isfinite(spectra['noise'][good]) & (spectra['noise'][good] > 0)
    ):
        raise ValueError('Invalid fitted pixels or calibrated noise in saved spectra')
    scale = spectra['normalization_scale']
    if scale.shape != (nbin,) or not np.all(np.isfinite(scale) & (scale > 0)):
        raise ValueError('Invalid Script 03 spectral normalization')
    if (noise_before.shape != spectra['noise'].shape or residual_scale.shape != (nbin,)
            or not np.all(np.isfinite(residual_scale) & (residual_scale > 0))
            or not np.allclose(spectra['noise'][good], (noise_before*residual_scale[:, None])[good],
                               rtol=3e-7, atol=0)):
        raise ValueError('Saved calibrated errors disagree with the Script 03 residual multipliers')
    spectra['good'] = good

    ids = np.asarray(predictions['BIN_ID'])
    if (ids.ndim != 1 or ids.dtype.kind not in 'iu' or not len(ids)
            or np.any((ids < 0) | (ids >= nbin)) or len(np.unique(ids)) != len(ids)
            or ids.tolist() != fit['bin_ids'] or ids.tolist() != fit['info']['selected_bins']):
        raise ValueError('The selected Script 04 bin IDs disagree with the fit or Script 03 bin domain')
    velocities = np.column_stack((predictions['VA_GLOBAL_KMS'], predictions['VB_GLOBAL_KMS']))
    for column, key, axis in zip(velocities.T, ('va', 'vb'), axes[:2]):
        if (not np.all(np.isfinite(column) & (column >= axis[0]) & (column <= axis[-1]))
                or not np.allclose(column, fit[key], rtol=0, atol=1e-8)):
            raise ValueError('Script 04 prediction velocities disagree with the fit or exceed Script 03 support')
    nodes = np.asarray(fit['nodes'], dtype=float)
    theta = np.asarray(fit['theta'], dtype=float)
    if (nodes.ndim != 1 or len(nodes) < 1 or not np.all(np.isfinite(nodes) & (nodes > 0))
            or not np.all(np.diff(nodes) > 0) or theta.shape != (5 + 2*len(nodes),)
            or not np.all(np.isfinite(theta)) or not 0 < theta[3] < 90
            or np.any(theta[5:] < 0) or fit['sign_a'] not in (-1, 1)):
        raise ValueError('Malformed selected Script 04 geometry or radial curves')
    data_score = np.asarray(predictions['DATA_DELTA_CHI2'], dtype=float)
    if not np.all(np.isfinite(data_score)) or not np.allclose(data_score, fit['bin_score'], atol=1e-8, rtol=0):
        raise ValueError('Script 04 data scores disagree with the selected fit')

    with np.load(paths['coordinates'], allow_pickle=False) as z:
        east, north = z['x_arcsec'], z['y_arcsec']
    spatial_hash = hashlib.sha256(east.tobytes() + north.tobytes()).hexdigest()
    if spatial_hash != provenance['spatial_sha256'] or east.ndim != 2 or east.shape != north.shape:
        raise ValueError('RH3 coordinates differ from the Script 04 geometry')
    with np.load(paths['membership'], allow_pickle=False) as z:
        y, x, bins = (z[f'rh3_{key}'] for key in ('y_pix', 'x_pix', 'bin_id'))
    if (not (y.shape == x.shape == bins.shape) or y.ndim != 1
            or any(a.dtype.kind not in 'iu' for a in (y, x, bins))
            or not np.array_equal(np.unique(bins), np.arange(nbin))
            or np.any((y < 0) | (y >= east.shape[0]) | (x < 0) | (x >= east.shape[1]))
            or len(np.unique(np.ravel_multi_index((y, x), east.shape))) != len(y)
            or not np.all(np.isfinite(east[y, x]) & np.isfinite(north[y, x]))):
        raise ValueError('RH3 spatial membership is inconsistent with the saved bins and coordinates')
    bin_map = np.full(east.shape, -1, dtype=int)
    bin_map[y, x] = bins

    lo, hi = map(float, manifest['sigma_bounds_kms'])
    baseline = float(cfg.MANGA_SIGMA_BASE_KMS)
    if not 0 < lo < hi or not np.isfinite([lo, hi]).all() or not lo <= baseline <= hi:
        raise ValueError('Set MANGA_SIGMA_BASE_KMS inside the inherited Script 03 dispersion bounds')
    if float(getattr(cfg, 'RH3_FA_FINE_STEP', 0.01)) != 0.01:
        raise ValueError('The Script 05 protocol requires RH3_FA_FINE_STEP = 0.01')
    settings = dict(baseline=baseline, sigma_bounds=[lo, hi],
                    sigma_warning=float(recorded_cfg.RH3_SIGMA_BOUNDARY_WARNING_KMS),
                    fraction_tolerance=float(recorded_cfg.RH3_FRACTION_CONSTRAINT_TOLERANCE),
                    degree=int(manifest['degree']), mdegree=int(manifest['mdegree']),
                    regul=float(manifest['regul']), velscale=velscale)
    if (not np.isfinite([settings['sigma_warning'], settings['fraction_tolerance']]).all()
            or settings['sigma_warning'] < 0 or settings['fraction_tolerance'] <= 0
            or settings['degree'] < -1 or settings['mdegree'] < 0 or settings['regul'] != 0):
        raise ValueError('Invalid or unsupported inherited Script 03 fit controls')

    setup = dict(settings, wavelength=wavelength, template_wavelength=template_wave,
                 templates=np.concatenate((basis, basis), axis=1),
                 component=np.repeat([0, 1], basis.shape[1]))
    edge_cells = int(getattr(cfg, 'GRID_EDGE_WARNING_CELLS', 2))
    if edge_cells < 0:
        raise ValueError('GRID_EDGE_WARNING_CELLS must be nonnegative')
    contract05 = dict(target=str(cfg.TARGET_NAME), redshift=float(cfg.REDSHIFT),
                      inputs=hashes, script04_fingerprint=fingerprint, model=name,
                      settings=settings, fractions=FRACTIONS.tolist(),
                      fraction_definition=manifest['fraction_definition'],
                      normalization_range_angstrom=normalization_range,
                      noise_model_hash=noise_hash, grid_edge_cells=edge_cells)
    return dict(setup=setup, spectra=spectra, fit=fit, ids=ids, velocities=velocities,
                coarse_min=coarse_min, interpolated_chi2=coarse_min[ids]+data_score, axes=axes,
                east=east, north=north, bin_map=bin_map, paths=paths, contract=contract05)
