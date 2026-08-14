r"""Noise, variance-rescaling, covariance, and correlated-realization helpers.

Because CRD_DAP converts :math:`\Delta\chi^2` into relative-likelihood weights,
the uncertainty scale is scientifically important. A multiplicative error in
the noise estimate changes the *width* of the likelihood basin even when the
best-fitting cell is unchanged.
"""

from __future__ import annotations

import numpy as np


def robust_variance_scale(residual: np.ndarray, sigma: np.ndarray) -> float:
    """Estimate a scalar uncertainty rescaling from normalized residuals.

    This simple helper is not the full covariance treatment. It is intended as
    one QC measure: if ``residual/sigma`` has a robust width far from unity, the
    formal variance scale should be investigated before likelihood inference.
    """
    r = np.asarray(residual, dtype=float)
    s = np.asarray(sigma, dtype=float)
    good = np.isfinite(r) & np.isfinite(s) & (s > 0)
    z = r[good] / s[good]
    if z.size < 10:
        raise ValueError("Too few valid residual samples to estimate variance scale.")
    median = np.median(z)
    mad = np.median(np.abs(z - median))
    robust_sigma = 1.4826 * mad
    return float(robust_sigma**2)


def estimate_spectral_covariance(*args, **kwargs):
    """Estimate wavelength-pixel covariance introduced by reduction/resampling.

    The exact estimator must be chosen after inspecting real stacked KCWI cubes.
    The result should support either a covariance matrix supplied directly to
    pPXF or a validated approximation when a full matrix is impractical.
    """
    raise NotImplementedError("Implemented after Script-1 residual/noise inspection.")


def generate_correlated_noise(covariance: np.ndarray, *, rng: np.random.Generator | None = None) -> np.ndarray:
    """Draw one zero-mean multivariate Gaussian noise realization."""
    cov = np.asarray(covariance, dtype=float)
    if cov.ndim != 2 or cov.shape[0] != cov.shape[1]:
        raise ValueError("covariance must be square.")
    rng = np.random.default_rng() if rng is None else rng
    return rng.multivariate_normal(np.zeros(cov.shape[0]), cov)
