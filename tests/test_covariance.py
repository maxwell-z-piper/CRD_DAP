from __future__ import annotations

import numpy as np

from crd_utils import covariance


def test_build_whitener_matches_direct_quadratic_form_and_decouples_masked_pixels():
    n = 48
    sigma = np.linspace(0.8, 1.2, n)
    good = np.ones(n, dtype=bool)
    good[[7, 19, 20, 35]] = False
    rho = np.zeros((1, 6), dtype=float)
    rho[0, 0] = 1.0
    rho[0, 1] = 0.25
    rho[0, 2] = 0.08
    result = covariance.build_inverse_cholesky(
        sigma, good, scale=1.3, rho_by_block=rho, eigen_floor=1e-10
    )
    W = result.inv_cholesky
    assert W.shape == (n, n)
    assert np.max(np.abs(W[np.ix_(good, ~good)])) < 1e-10

    rng = np.random.default_rng(4)
    r = rng.normal(size=n)
    r[~good] = 123.0  # must not affect fitted rows
    q = W @ r
    chi_w = np.sum(q[np.flatnonzero(good)] ** 2)
    chi_fun = covariance.covariance_total_chi2(r, W, np.flatnonzero(good))
    np.testing.assert_allclose(chi_w, chi_fun, rtol=1e-12, atol=1e-12)


def test_bootstrap_band_uses_powerbins_and_detects_zero_consistency():
    rng = np.random.default_rng(10)
    nbin, nlag = 100, 8
    rho = rng.normal(0.0, 0.04, size=(nbin, nlag))
    rho[:, 0] = 1.0
    rho[:, 1] += 0.30
    band = covariance.bootstrap_simultaneous_band(
        rho, n_bootstrap=500, confidence=0.95, random_seed=9
    )
    assert band.lower[1] > 0.0
    assert band.consistent_with_zero[1] is np.False_ or not bool(band.consistent_with_zero[1])
    # At least one of the deliberately zero-centered longer lags should include 0.
    assert np.any(band.consistent_with_zero[2:])


def test_equal_wavelength_blocks_cover_every_pixel_once():
    b = covariance.equal_wavelength_blocks(101, 3)
    assert b.shape == (101,)
    assert set(np.unique(b)) == {0, 1, 2}
    assert np.max(np.bincount(b)) - np.min(np.bincount(b)) <= 1


def test_lag_correlation_recovers_short_range_signal():
    rng = np.random.default_rng(12)
    n = 5000
    x = np.zeros(n)
    for j in range(1, n):
        x[j] = 0.45 * x[j - 1] + rng.normal()
    rho, count = covariance.lag_correlation(
        x, np.ones(n, dtype=bool), max_lag=5, min_pairs=50
    )
    assert count[0, 1] > 4000
    assert 0.35 < rho[0, 1] < 0.55
    assert abs(rho[0, 5]) < abs(rho[0, 1])

def test_regularized_covariance_still_decouples_masked_pixels():
    """
    Regression test for the real Script-3 failure.

    Use a deliberately strong, oscillatory lag-correlation curve that does not
    define a positive-definite covariance matrix.  This forces the
    eigenvalue-regularization branch.

    The regularized fitted-pixel covariance must remain completely independent
    of excluded wavelength pixels.
    """

    n_pixel = 64

    sigma = np.ones(
        n_pixel,
        dtype=float,
    )

    good = np.ones(
        n_pixel,
        dtype=bool,
    )

    # Mimic internal wavelength-mask gaps.
    good[
        [7, 19, 20, 35]
    ] = False

    # Deliberately strong positive short-lag correlation followed by an
    # anticorrelated tail.  This is qualitatively similar to the structure
    # encountered in the real RH3 residuals and forces positive-definite repair.
    rho = np.array(
        [[
            1.00,
            0.95,
            0.84,
            0.70,
            0.57,
            0.44,
            0.32,
            0.19,
            0.08,
            0.00,
           -0.04,
           -0.08,
           -0.12,
           -0.14,
           -0.17,
           -0.22,
           -0.27,
           -0.35,
           -0.43,
           -0.47,
           -0.52,
        ]],
        dtype=float,
    )

    result = covariance.build_inverse_cholesky(
        sigma,
        good,
        scale=1.0,
        rho_by_block=rho,
        eigen_floor=1.0e-8,
    )

    W = result.inv_cholesky

    fitted_idx = np.flatnonzero(good)
    excluded_idx = np.flatnonzero(~good)

    # The deliberately pathological input should have required repair.
    assert result.regularized

    # Fitted whitened residuals must not depend on excluded pixels.
    coupling = W[
        np.ix_(
            fitted_idx,
            excluded_idx,
        )
    ]

    assert np.max(
        np.abs(coupling)
    ) < 1.0e-10

    # Stronger physical test:
    #
    # changing masked residuals by an absurd amount must leave fitted chi^2
    # completely unchanged.
    rng = np.random.default_rng(123)

    residual_1 = rng.normal(
        size=n_pixel,
    )

    residual_2 = residual_1.copy()

    residual_1[excluded_idx] = 0.0
    residual_2[excluded_idx] = 1.0e6

    chi2_1 = covariance.covariance_total_chi2(
        residual_1,
        W,
        fitted_idx,
    )

    chi2_2 = covariance.covariance_total_chi2(
        residual_2,
        W,
        fitted_idx,
    )

    np.testing.assert_allclose(
        chi2_1,
        chi2_2,
        rtol=1.0e-12,
        atol=1.0e-12,
    )
