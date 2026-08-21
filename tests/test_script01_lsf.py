import numpy as np

from crd_utils.psf_lsf import measure_arc_lsf


def test_measure_arc_lsf_recovers_gaussian_width():
    ny, nx = 180, 48
    wave_axis = 8400.0 + 0.20 * np.arange(ny)
    wavemap = np.repeat(wave_axis[:, None], nx, axis=1)

    slicemap = np.zeros((ny, nx), dtype=int)
    posmap = np.zeros((ny, nx), dtype=float)
    for x in range(nx):
        slicemap[:, x] = x // 16
        posmap[:, x] = x % 16

    sigma_a = 0.50
    lines = [8407.0, 8414.0, 8421.0, 8428.0]
    arc = np.zeros((ny, nx), dtype=float) + 10.0
    for center in lines:
        profile = 500.0 * np.exp(-0.5 * ((wave_axis - center) / sigma_a) ** 2)
        arc += profile[:, None]

    rng = np.random.default_rng(3)
    arc += rng.normal(0.0, 1.0, arc.shape)

    result = measure_arc_lsf(
        arc,
        wavemap,
        slicemap,
        posmap,
        wavegood0=8402.0,
        wavegood1=8433.0,
        polynomial_order=1,
        measure_spatial_variation=True,
        spatial_bins=2,
        peak_prominence_fraction=0.02,
        peak_height_percentile=60.0,
        min_peak_distance_pix=10,
        line_fit_half_width_pix=8,
        min_good_lines=6,
    )

    expected_fwhm = 2.354820045 * sigma_a
    assert result.n_lines_used >= 6
    assert np.isclose(np.median(result.fwhm_angstrom), expected_fwhm, atol=0.25)

    # The instrument-good interval can extend beyond the accepted line centers.
    # The LSF model must not silently extrapolate beyond its empirical support.
    assert result.measurement_wavelength_min >= result.wavelength_min
    assert result.measurement_wavelength_max <= result.wavelength_max
    below = result.measurement_wavelength_min - 1.0
    assert np.isnan(result.evaluate_fwhm(below))
    assert np.isfinite(result.evaluate_fwhm(below, allow_extrapolation=True))
