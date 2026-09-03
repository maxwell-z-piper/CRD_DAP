"""Regression tests for the CRD_DAP pPXF-9.4.8 cached-whitener patch.

Run *inside the CRD_DAP conda environment after applying the patch*::

    python -m pytest -q ppxf_patch_9_4_8/tests/test_ppxf_cached_whitener.py

The central scientific test is direct equivalence between:

1. stock pPXF's full covariance input C; and
2. the patched cached path supplied with W=L^-1, where C=L L^T.

The two paths must produce the same objective and fit to numerical precision.
"""
from __future__ import annotations

import inspect
from unittest import mock

import numpy as np
import pytest
from scipy import linalg

import ppxf
import ppxf.ppxf as ppxf_module
from ppxf.ppxf import ppxf as PPXF

EXPECTED_VERSION = "9.4.8"


def _require_supported_installation():
    assert str(ppxf.__version__) == EXPECTED_VERSION, (
        f"Tests require pPXF {EXPECTED_VERSION}; found {ppxf.__version__}."
    )
    assert "noise_inv_cholesky" in inspect.signature(PPXF).parameters, (
        "CRD_DAP cached-whitener patch is not installed in this Python environment."
    )


def _problem(npix: int = 128, rho: float = 0.28, seed: int = 11):
    """Construct a deterministic small full-spectrum fitting problem."""
    rng = np.random.default_rng(seed)
    x = np.linspace(-1.0, 1.0, npix)

    # Two smooth, non-degenerate template shapes with absorption-like features.
    t1 = 1.0 - 0.24*np.exp(-0.5*((x + 0.23)/0.055)**2) - 0.15*np.exp(-0.5*((x - 0.34)/0.075)**2)
    t2 = 1.0 - 0.18*np.exp(-0.5*((x + 0.08)/0.085)**2) - 0.11*np.exp(-0.5*((x - 0.47)/0.045)**2)
    templates = np.column_stack([t1, t2])

    sigma = 0.018*(1.0 + 0.18*(x + 1.0)/2.0)
    lag = np.abs(np.subtract.outer(np.arange(npix), np.arange(npix)))
    corr = rho**lag
    cov = np.outer(sigma, sigma)*corr
    chol = linalg.cholesky(cov, lower=True, check_finite=True)
    inv_chol = linalg.solve_triangular(
        chol, np.eye(npix), lower=True, check_finite=True
    )

    noiseless = 0.62*t1 + 0.38*t2
    galaxy = noiseless + chol @ rng.standard_normal(npix)
    errvec = np.sqrt(np.diag(cov))
    goodpixels = np.arange(4, npix - 4, dtype=int)

    return templates, galaxy, errvec, cov, inv_chol, goodpixels


def _fit_one_component(*, templates, galaxy, noise, goodpixels, **kwargs):
    return PPXF(
        templates,
        galaxy,
        noise,
        55.0,
        start=[0.0, 165.0],
        moments=2,
        degree=2,
        mdegree=0,
        goodpixels=goodpixels,
        clean=False,
        quiet=True,
        **kwargs,
    )


def _fit_two_component_fixed_velocity(*, templates, galaxy, noise, goodpixels, **kwargs):
    templates2 = np.column_stack([templates, templates])
    component = np.array([0, 0, 1, 1], dtype=int)
    return PPXF(
        templates2,
        galaxy,
        noise,
        55.0,
        start=[[-110.0, 150.0], [125.0, 150.0]],
        moments=[2, 2],
        component=component,
        fraction=0.55,
        fixed=[[True, False], [True, False]],
        bounds=[
            [[-110.001, -109.999], [20.0, 320.0]],
            [[124.999, 125.001], [20.0, 320.0]],
        ],
        degree=2,
        mdegree=0,
        goodpixels=goodpixels,
        clean=False,
        quiet=True,
        **kwargs,
    )


def _assert_fit_equivalent(a, b):
    np.testing.assert_allclose(np.asarray(a.sol, dtype=float), np.asarray(b.sol, dtype=float), rtol=2e-9, atol=2e-9)
    np.testing.assert_allclose(np.asarray(a.bestfit, dtype=float), np.asarray(b.bestfit, dtype=float), rtol=2e-9, atol=2e-9)
    np.testing.assert_allclose(np.asarray(a.weights, dtype=float), np.asarray(b.weights, dtype=float), rtol=2e-9, atol=2e-9)
    np.testing.assert_allclose(float(a.chi2), float(b.chi2), rtol=2e-9, atol=2e-9)


def test_patch_signature_and_version():
    _require_supported_installation()


def test_cached_whitener_matches_stock_covariance_one_component():
    _require_supported_installation()
    templates, galaxy, errvec, cov, inv_chol, goodpixels = _problem()

    stock = _fit_one_component(
        templates=templates, galaxy=galaxy, noise=cov, goodpixels=goodpixels
    )
    cached = _fit_one_component(
        templates=templates,
        galaxy=galaxy,
        noise=errvec,  # valid placeholder; patched pPXF replaces internal weighting with W
        goodpixels=goodpixels,
        noise_inv_cholesky=inv_chol,
    )

    _assert_fit_equivalent(stock, cached)
    np.testing.assert_allclose(cached.noise, inv_chol, rtol=0.0, atol=0.0)
    assert cached.crd_dap_noise_inv_cholesky is True


def test_cached_whitener_matches_stock_covariance_fixed_two_component_state():
    """Exercise the pPXF mode closest to a CRD_DAP Script-3 grid state."""
    _require_supported_installation()
    templates, galaxy, errvec, cov, inv_chol, goodpixels = _problem(seed=19)

    stock = _fit_two_component_fixed_velocity(
        templates=templates, galaxy=galaxy, noise=cov, goodpixels=goodpixels
    )
    cached = _fit_two_component_fixed_velocity(
        templates=templates,
        galaxy=galaxy,
        noise=errvec,
        goodpixels=goodpixels,
        noise_inv_cholesky=inv_chol,
    )

    _assert_fit_equivalent(stock, cached)


def test_cached_path_skips_cholesky_inside_ppxf():
    _require_supported_installation()
    templates, galaxy, errvec, cov, inv_chol, goodpixels = _problem(npix=96)

    real_cholesky = ppxf_module.linalg.cholesky
    with mock.patch.object(ppxf_module.linalg, "cholesky", wraps=real_cholesky) as wrapped:
        _fit_one_component(
            templates=templates, galaxy=galaxy, noise=cov, goodpixels=goodpixels
        )
        assert wrapped.call_count >= 1

    with mock.patch.object(ppxf_module.linalg, "cholesky", wraps=real_cholesky) as wrapped:
        _fit_one_component(
            templates=templates,
            galaxy=galaxy,
            noise=errvec,
            goodpixels=goodpixels,
            noise_inv_cholesky=inv_chol,
        )
        assert wrapped.call_count == 0


def test_default_error_vector_path_is_unchanged():
    _require_supported_installation()
    templates, galaxy, errvec, _cov, _inv_chol, goodpixels = _problem(seed=23)
    fit = _fit_one_component(
        templates=templates, galaxy=galaxy, noise=errvec, goodpixels=goodpixels
    )
    assert fit.crd_dap_noise_inv_cholesky is False
    np.testing.assert_allclose(fit.noise, errvec, rtol=0.0, atol=0.0)


@pytest.mark.parametrize(
    "bad_builder,match",
    [
        (lambda w: w[:-1, :-1], "shape"),
        (lambda w: np.where(np.eye(w.shape[0], dtype=bool), np.nan, w), "finite"),
        (lambda w: w + np.triu(np.ones_like(w), 1)*1e-3, "lower-triangular"),
        (lambda w: w*np.where(np.eye(w.shape[0], dtype=bool), -1.0, 1.0), "positive diagonal"),
    ],
)
def test_invalid_cached_whitener_is_rejected(bad_builder, match):
    _require_supported_installation()
    templates, galaxy, errvec, _cov, inv_chol, goodpixels = _problem(npix=72)
    bad = bad_builder(inv_chol.copy())
    with pytest.raises(ValueError, match=match):
        _fit_one_component(
            templates=templates,
            galaxy=galaxy,
            noise=errvec,
            goodpixels=goodpixels,
            noise_inv_cholesky=bad,
        )
