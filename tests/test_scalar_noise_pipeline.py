"""Scientific checks for scalar calibration and its diagonal pPXF contract."""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

from astropy.io import fits
import numpy as np
import pytest

from crd_utils import binning, ppxf_utils, scalar_noise, spatial_noise


@dataclass
class Cube:
    arm: str
    flux: np.ndarray
    variance: np.ndarray
    good: np.ndarray
    good_spaxel: np.ndarray
    nwave: int
    wavelength: np.ndarray
    header: dict


def _hash(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_prepared_variance_allows_only_excluded_zero_padding():
    @dataclass
    class Prepared:
        arm: str
        variance: np.ndarray
        good: np.ndarray

    reference = np.array([[[2., 0., 3.]]])
    prepared = Prepared('BL', np.array([[[2., np.nan, 3.]]]),
                        np.array([[[True, False, True]]]))
    spatial_noise._check_prepared_variance(prepared, reference)

    prepared.variance[0, 0, 2] = 4.
    with pytest.raises(ValueError, match='1 samples .*1 marked good'):
        spatial_noise._check_prepared_variance(prepared, reference)

    prepared.variance[0, 0, 2] = 3.
    prepared.good[0, 0, 1] = True
    with pytest.raises(ValueError, match='1 samples .*1 marked good'):
        spatial_noise._check_prepared_variance(prepared, reference)


def test_correlated_spaxel_mock_replay_inflates_bin_noise_once(tmp_path):
    rng = np.random.default_rng(43)
    wave = np.arange(40.) + 5000.
    flux = np.zeros((1, 2, wave.size))
    variance = np.ones_like(flux)
    good = np.ones_like(flux, dtype=bool)
    good[:, :, 4] = False  # an inherited atmospheric/bad-channel mask
    cube = Cube('BL', flux, variance, good, np.ones((1, 2), bool),
                wave.size, wave, {'BUNIT': 'flux/arcsec2'})
    bins = np.zeros((1, 2), dtype=int)
    formal = binning.coadd_bin_spectra(cube, bins, n_bins=1,
                                      pixel_area_arcsec2=.09)
    science_i = tmp_path / 'science_icubes.fits'
    science_v = tmp_path / 'science_vcubes.fits'
    science_m = tmp_path / 'science_mcubes.fits'
    science_e = tmp_path / 'science_ecubes.fits'
    fits.PrimaryHDU(np.moveaxis(flux, -1, 0)).writeto(science_i)
    fits.PrimaryHDU(np.moveaxis(variance, -1, 0)).writeto(science_v)
    fits.PrimaryHDU(np.zeros((wave.size,1,2))).writeto(science_m)
    fits.PrimaryHDU(np.ones((wave.size,1,2))).writeto(science_e)
    mock_paths = []
    for i in range(80):
        # Both spaxels share each draw, so the sum has variance 4 rather than 2.
        draw = rng.normal(size=wave.size)
        mock = np.repeat(draw[None, None, :], 2, axis=1)
        path = tmp_path / f'mock_{i:03d}.fits'
        fits.PrimaryHDU(np.moveaxis(mock, -1, 0)).writeto(path)
        mock_paths.append(str(path))
    manifest = tmp_path / 'manifest.json'
    manifest.write_text(json.dumps(dict(
        schema='crd-noise-stack-ensemble-v1', complete=True, arm='BL',
        science_replay_passed=True, noise_icubes=mock_paths,
        noise_sha256={str(p):_hash(tmp_path / p.split('/')[-1]) for p in mock_paths},
        reference_products={'i':str(science_i), 'v':str(science_v),
                            'm':str(science_m), 'e':str(science_e)},
        reference_sha256={'i':_hash(science_i), 'v':_hash(science_v),
                          'm':_hash(science_m), 'e':_hash(science_e)})))
    corrected, scale = spatial_noise.calibrate(
        cube, bins, formal, manifest_path=manifest, n_bins=1,
        pixel_area_arcsec2=.09, min_member_fraction=1.,
        rest_range=(5000.,5039.), redshift=0.,
        output_path=tmp_path/'calibration.npz')
    assert 1.2 < scale[0] < 1.6
    np.testing.assert_allclose(corrected.uncertainty[formal.good],
                               formal.uncertainty[formal.good] * scale[0])
    assert not corrected.good[0,4]
    assert np.isnan(corrected.uncertainty[0,4])
    # pPXF receives the corrected sigma; chi-square must use it exactly once.
    total = ppxf_utils.total_chi2(np.ones(5), np.zeros(5),
                                   np.full(5, scale[0]), np.arange(5))
    assert np.isclose(total, 5 / scale[0]**2)


@pytest.mark.parametrize('residual_scale', [.8, 1.2])
def test_full_grid_best_residual_scales_noise_and_chi2_once(monkeypatch, tmp_path, residual_scale):
    script_path = Path(__file__).resolve().parents[1] / 'scripts' / '03_build_RH3_likelihood_cubes.py'
    spec = importlib.util.spec_from_file_location('crd_script03_full_grid_test', script_path)
    script = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(script)

    calls = {'grid': 0, 'best': 0}

    def single_fit(**kwargs):
        return SimpleNamespace(
            success=True, chi2_total=10., reduced_chi2=2.,
            velocity=np.array([0.]), sigma=np.array([40.]),
            bestfit=np.zeros(5), weights=np.array([1.]),
        )

    def full_grid(**kwargs):
        calls['grid'] += 1
        return SimpleNamespace(
            best_index=(0, 0, 0),
            chi2_total=np.array([[[3., 4.]]]),
            reduced_chi2=np.array([[[.6, .8]]]),
            sigma_a=np.array([[[50., 50.]]]),
            sigma_b=np.array([[[60., 60.]]]),
            fit_status=np.zeros((1, 1, 2), dtype=np.int8),
            sigma_boundary=np.zeros((1, 1, 2), dtype=np.uint8),
            n_failures=0, n_ppxf_failures=0,
            n_fixed_velocity_mismatch=0, n_sigma_boundary=0,
        )

    def best_fit(**kwargs):
        calls['best'] += 1
        galaxy = kwargs['galaxy']
        noise = kwargs['noise']
        z = (residual_scale / 1.4826) * np.array([-2., -1., 0., 1., 2.])
        model = galaxy - noise * z
        return SimpleNamespace(
            success=True, bestfit=model,
            chi2_total=ppxf_utils.total_chi2_diagonal(galaxy, model, noise,
                                                        kwargs['goodpixels']),
            reduced_chi2=1., sigma=np.array([50., 60.]),
            weights=np.array([.6, .4]), polyweights=None,
        )

    monkeypatch.setattr(script.ppxf_utils, 'fit_single_losvd', single_fit)
    monkeypatch.setattr(script.ppxf_grid, 'build_rh3_likelihood_cube', full_grid)
    monkeypatch.setattr(script.ppxf_utils, 'fit_fixed_two_component_state', best_fit)
    script._WORKER = dict(
        min_good_pixels=3, templates_single=np.ones((5, 1)),
        templates_two=np.ones((5, 2)), component=np.array([0, 1]), n_basis=1,
        velscale=1., lam=np.arange(5.) + 7000., lam_temp=np.arange(5.) + 7000.,
        va_grid=np.array([1.]), vb_grid=np.array([-1.]), fa_grid=np.array([.4, .6]),
        sigma_start=50., sigma_bounds=(20., 100.), sigma_boundary_tolerance=2.,
        single_start_velocity=0., single_velocity_bounds=(-2., 2.),
        degree=0, mdegree=0, regul=0., noise_model_hash='unit-test',
    )

    result = script._fit_one_bin_worker((0, np.arange(10., 15.), np.ones(5),
                                         np.ones(5, bool), np.ones(5)))
    assert calls == {'grid': 1, 'best': 1}
    assert result['residual_noise_scale'] == pytest.approx(residual_scale)
    np.testing.assert_allclose(result['noise'],
                               result['noise_before_residual_scale'] * residual_scale)
    np.testing.assert_allclose(result['chi2_total'],
                               np.array([[[3., 4.]]]) / residual_scale**2)
    assert result['one_chi2_total'] == pytest.approx(10. / residual_scale**2)
    assert result['best_chi2_total'] == pytest.approx(
        ppxf_utils.total_chi2_diagonal(result['galaxy'], result['best_bestfit'],
                                        result['noise'], np.arange(5)))

    checkpoint = tmp_path / 'bin_0000.npz'
    script._write_checkpoint(checkpoint, result)
    assert script._checkpoint_is_valid(checkpoint, (1, 1, 2), 'unit-test')
    assert not script._checkpoint_is_valid(checkpoint, (1, 1, 2), 'old-method')

    run = SimpleNamespace(products_dir=tmp_path, metadata_dir=tmp_path)
    product = scalar_noise.write_calibration(
        run, scales=np.array([result['residual_noise_scale']]),
        centers=np.array([result['residual_center']]),
        residuals=result['normalized_residual'][None, :], input_hash='unit-test')
    with np.load(product, allow_pickle=False) as saved:
        assert saved['raw_scale'][0] == pytest.approx(residual_scale)
        assert saved['adopted_scale'][0] == pytest.approx(residual_scale)
        assert not bool(saved['floor_applied'])
