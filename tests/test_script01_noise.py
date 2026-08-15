import numpy as np

from crd_utils.noise import estimate_spectral_correlation, characterize_preliminary_noise


def test_spectral_correlation_white_noise_is_small():
    rng = np.random.default_rng(2)
    z = rng.normal(size=(100, 200))
    lags, corr = estimate_spectral_correlation(z, max_lag=5)
    assert lags[0] == 0
    assert np.isclose(corr[0], 1.0)
    assert np.nanmax(np.abs(corr[1:])) < 0.15


def test_preliminary_noise_diagnostic_runs():
    rng = np.random.default_rng(4)
    ny, nx, nw = 6, 7, 101
    wave = np.linspace(-1, 1, nw)
    smooth = 20.0 + 2.0 * wave + 0.5 * wave**2
    sigma = np.ones((ny, nx, nw))
    flux = smooth[None, None, :] + rng.normal(size=(ny, nx, nw))
    good = np.ones_like(flux, dtype=bool)
    image = np.median(flux, axis=-1)
    good_spaxel = np.ones((ny, nx), dtype=bool)

    result = characterize_preliminary_noise(
        flux,
        sigma,
        good,
        image,
        good_spaxel,
        max_spaxels=20,
        low_flux_percentile=50,
        savgol_window=21,
        max_lag=5,
    )
    assert result.n_spaxels_used > 0
    assert result.n_samples_used > 100
    assert result.variance_scale_factor > 0
