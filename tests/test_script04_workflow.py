"""Exercise the real CLI, spawn-safe workers, review gate and frozen outputs.

Small artificial likelihoods keep this an implementation check, not a test of
whether these observations can constrain a galaxy's rotation curve.
"""
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
from astropy.io import fits
from astropy.table import Table

from crd_utils.global_disk_io import save_json


def test_review_then_eight_models_and_resume(tmp_path):
    root = Path(__file__).resolve().parents[1]
    runs = [tmp_path / f'run{i}' for i in (1, 2, 3)]
    out = tmp_path / 'global'
    for run in runs:
        (run / 'products').mkdir(parents=True)
        (run / 'metadata').mkdir()
    yy, xx = np.indices((7, 7))
    east, north = xx - 3., yy - 3.
    light = np.exp(-(east**2 + north**2) / 12.)
    member = ((north >= 0).astype(int) * 2 + (east >= 0)).ravel()
    header = fits.Header(dict(NAXIS=3, CTYPE1='RA---TAN', CTYPE2='DEC--TAN',
                              CTYPE3='WAVE', CRVAL1=10., CRVAL2=20.,
                              CRVAL3=8000., CRPIX1=4., CRPIX2=4., CRPIX3=1.,
                              CDELT1=-1/3600, CDELT2=1/3600, CDELT3=1., IFUPA=0.))
    cube = np.broadcast_to(light, (5, 7, 7)).copy()
    fits.HDUList([
        fits.PrimaryHDU(cube, header),
        fits.ImageHDU(np.arange(8000., 8005.), name='WAVELENGTH'),
        fits.ImageHDU(np.ones_like(cube, dtype='uint8'), name='GOODMASK'),
        fits.ImageHDU(np.ones(5, dtype='uint8'), name='GOODWAVE'),
        fits.ImageHDU(np.ones((7, 7), dtype='uint8'), name='GOODSPAX'),
    ]).writeto(runs[0] / 'products/prepared_RH3.fits')
    np.savez(runs[0] / 'products/RH3_spatial_coordinates.npz', x_arcsec=east, y_arcsec=north)
    np.savez(runs[1] / 'products/master_bin_membership.npz',
             rh3_y_pix=yy.ravel(), rh3_x_pix=xx.ravel(), rh3_bin_id=member,
             rh3_x_arcsec=east.ravel(), rh3_y_arcsec=north.ravel())
    save_json(runs[1] / 'metadata/script02_manifest.json',
              dict(target='synthetic', rh3_window_observed_angstrom=[8000., 8004.]))
    save_json(runs[2] / 'metadata/script03_manifest.json', dict(target='synthetic', redshift=0.))
    axis = np.linspace(-200., 200., 9)
    aa, bb, ff = np.meshgrid(axis, axis, [.2, .5, .8], indexing='ij')
    chi = np.array([((aa-sign*50)/30)**2 + ((bb+sign*50)/30)**2 + (ff-.5)**2
                    for sign in (-1, -1, 1, 1)])
    np.savez(runs[2] / 'products/RH3_likelihood_cubes.npz', chi2_total=chi,
             VA_grid=axis, VB_grid=axis, fA_grid=[.2, .5, .8], one_velocity=[-20., -20., 20., 20.])
    settings = dict(
        TARGET_NAME='synthetic', REDSHIFT=0., PA_KIN_INITIAL_DEG=0., PA_KIN_ERR_DEG=1.,
        PA_BOUND_SIGMA=3., X0_INITIAL_ARCSEC=0., Y0_INITIAL_ARCSEC=0.,
        VSYS_INITIAL_KMS=0., CENTER_BOUND_ARCSEC=.1, VSYS_BOUND_KMS=10.,
        SCRIPT04_INCLINATION_BOUNDS_DEG=(25., 35.), SCRIPT04_INCLINATION_INITIAL_DEG=30.,
        MANGA_SIGMA_PEAK_POS_KMS=80., MANGA_SIGMA_PEAK_NEG_KMS=80., MANGA_SIGMA_BASE_KMS=40.,
        R_2SIGMA_POS_INITIAL_ARCSEC=1., R_2SIGMA_NEG_INITIAL_ARCSEC=1.,
        RING_SPACE_ARCSEC=2., R_START_FWHM_FACTOR=1., RING_DELTA_FACTOR=.5,
        SCRIPT04_DISK_A_ROTATION_SIGN=1, SCRIPT04_VROT_MAX_KMS=400.,
        SCRIPT04_VTYP_KMS=200., SCRIPT04_PSF_KERNEL_RADIUS_ARCSEC=3., SCRIPT04_PSF_QUADRATURE=3,
        SCRIPT04_RANDOM_SEED=7, SCRIPT04_POPULATION_FACTOR=1, SCRIPT04_GLOBAL_ITERATIONS=0,
        SCRIPT04_LOCAL_STARTS=1, SCRIPT04_LOCAL_ITERATIONS=1,
        SCRIPT04_OPTIMIZER_TOLERANCE=1e-3, SCRIPT04_BASELINE_COEFFICIENT=.002,
    )
    config = tmp_path / 'synthetic.py'
    config.write_text('\n'.join(f'{key} = {value!r}' for key, value in settings.items()) + '\n')
    # The PSF measurement function has separate native-pixel tests; this test
    # starts at its documented saved-product contract with a known profile.
    save_json(out / 'products/stellar_psf_measurements.json', dict(measurements=[
        dict(path='synthetic standard', candidates=dict(gaussian=dict(
            success=True, flags=[], native_pixel_arcsec=[1., 1.], parameters=[.6, .6, 0., 1.])))]))
    command = [sys.executable, str(root / 'scripts/04_fit_global_two_disk_model.py'),
               '--config', str(config), '--output-dir', str(out), '--psf-profile', 'gaussian',
               '--workers', '2']
    for index, run in enumerate(runs, 1):
        command += [f'--script{index}-run', str(run)]
    environment = dict(os.environ, MPLCONFIGDIR=str(tmp_path / 'matplotlib'), PYTHONDONTWRITEBYTECODE='1')

    def invoke(phase, *extra, success=True):
        result = subprocess.run(command + ['--phase', phase, *extra], env=environment,
                                capture_output=True, text=True, timeout=180)
        assert (result.returncode == 0) == success, result.stdout + result.stderr
        return result

    invoke('compare', success=False)  # Review is required; no automatic choice.
    invoke('calibrate')
    assert (out / 'figures/smoothing_calibration_2x6.pdf').stat().st_size > 1000
    assert not (out / 'products/eight_model_comparison.json').exists()
    invoke('review', '--coefficient', '.002', '--note', 'software test only')
    invoke('compare')
    comparison = json.loads((out / 'products/eight_model_comparison.json').read_text())
    assert len(comparison['models']) == 8
    checkpoint = out / 'fits/full_psf1_smooth1.json'
    original = checkpoint.read_bytes()
    invoke('compare', '--resume')
    assert checkpoint.read_bytes() == original
    invoke('choose-model', '--model-name', 'full_psf1_smooth1', '--note', 'software test only')
    selected = json.loads((out / 'products/selected_global_model.json').read_text())
    table = Table.read(selected['bin_predictions'])
    assert len(table) == 4
    assert np.all(np.isfinite(table['VA_GLOBAL_KMS']))
    fit = json.loads(checkpoint.read_text())
    assert np.isclose(fit['objective'], fit['data_score'] + fit['penalty'])
    assert fit['n_bins'] == 4
    adopted_kernel = (out / 'products/adopted_RH3_PSF.npz').read_bytes()
    # A changed ceiling must not silently reuse previous completed fits.
    config.write_text(config.read_text().replace('SCRIPT04_VROT_MAX_KMS = 400.0',
                                                 'SCRIPT04_VROT_MAX_KMS = 600.0'))
    failure = invoke('calibrate', '--resume', success=False)
    assert 'inputs/code changed' in failure.stderr
    assert (out / 'products/adopted_RH3_PSF.npz').read_bytes() == adopted_kernel
