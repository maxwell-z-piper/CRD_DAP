"""Spectral covariance utilities for CRD_DAP.

This module contains the low-level statistical machinery used by Script 3 to
turn empirically measured residual correlations into a positive-definite noise
model that pPXF can use.  It deliberately does *not* decide which covariance
model is scientifically preferred; model calibration/selection lives in
``crd_utils.covariance_calibration``.

Notation
--------
For PowerBin ``i`` and log-wavelength pixel ``j``:

``r_ij = F_ij - M_ij``
    Residual between the observed spectrum and the calibration pPXF model.
``sigma_ij``
    Formal one-sigma uncertainty propagated by Script 2 and the Script-3
    overlap rebinning.
``z_ij = r_ij / sigma_ij``
    Normalized residual before empirical variance rescaling.
``s_i``
    Per-bin multiplicative standard-deviation scale.
``R_i``
    Correlation matrix with unit diagonal.
``D_i = diag(sigma_i)``
    Diagonal matrix of formal errors.
``C_i = s_i^2 D_i R_i D_i``
    Adopted covariance matrix.

The cached-pPXF interface consumes ``W_i = L_i^-1`` where
``C_i = L_i L_i^T`` is a lower-triangular Cholesky factorization.  The pPXF
patch therefore avoids recomputing an invariant Cholesky factor for every
likelihood-grid state in the same PowerBin.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np
from scipy import linalg


MAD_TO_SIGMA = 1.482602218505602


@dataclass(frozen=True)
class BootstrapBand:
    """Simultaneous bootstrap confidence band for a lag-correlation curve."""

    lag: np.ndarray
    center: np.ndarray
    lower: np.ndarray
    upper: np.ndarray
    n_bootstrap: int
    confidence: float

    @property
    def consistent_with_zero(self) -> np.ndarray:
        """Boolean vector: confidence band contains zero at each lag."""
        return (self.lower <= 0.0) & (self.upper >= 0.0)


@dataclass(frozen=True)
class WhiteningResult:
    """Dense inverse-Cholesky operator plus numerical-QC metadata."""

    inv_cholesky: np.ndarray
    min_eigenvalue_before_regularization: float
    eigenvalue_floor_applied: float
    regularized: bool


def robust_standard_deviation(values: np.ndarray) -> float:
    """Return a robust Gaussian-equivalent standard deviation using MAD.

    The median is removed first.  MAD is preferred to an ordinary RMS for the
    covariance-calibration scale because a small number of unmodelled line or
    atmospheric residuals should not be allowed to set the noise amplitude for
    an entire PowerBin.
    """
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    if x.size < 5:
        return np.nan
    med = float(np.median(x))
    mad = float(np.median(np.abs(x - med)))
    if not np.isfinite(mad):
        return np.nan
    if mad > 0:
        return float(MAD_TO_SIGMA * mad)
    # Degenerate fallback.  It is intentionally not the primary estimator.
    std = float(np.std(x, ddof=1)) if x.size > 1 else np.nan
    return std if np.isfinite(std) and std > 0 else np.nan


def lag_correlation(
    values: np.ndarray,
    good: np.ndarray,
    *,
    max_lag: int,
    min_pairs: int = 25,
    block_index: np.ndarray | None = None,
    n_blocks: int | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Measure correlation coefficient versus integer spectral lag.

    Parameters
    ----------
    values
        One residual-like vector.  For covariance estimation this is normally
        ``z = residual / formal_sigma`` (optionally divided by the per-bin scale
        ``s_i``).  For whitening QC it is the whitened residual vector.
    good
        Boolean fitted-pixel mask on the same log-wavelength grid.
    max_lag
        Maximum pixel separation to inspect.
    min_pairs
        Minimum number of valid pixel pairs required to report a lag.
    block_index
        Optional integer block label for every wavelength pixel.  When supplied,
        a pair contributes to the block containing the pair midpoint.  This
        avoids forcing correlations to zero merely because a pair straddles a
        block boundary.
    n_blocks
        Number of blocks.  Required when ``block_index`` is supplied.

    Returns
    -------
    rho, pair_count
        Arrays with shape ``(n_blocks, max_lag + 1)``.  With no block indexing,
        ``n_blocks=1``.  Lag zero is fixed to one whenever the block has enough
        valid samples.
    """
    x = np.asarray(values, dtype=float)
    g = np.asarray(good, dtype=bool)
    if x.ndim != 1 or g.ndim != 1 or x.size != g.size:
        raise ValueError("values and good must be one-dimensional arrays of identical length")
    max_lag = int(max_lag)
    min_pairs = int(min_pairs)
    if max_lag < 1:
        raise ValueError("max_lag must be >= 1")
    if min_pairs < 3:
        raise ValueError("min_pairs must be >= 3")

    if block_index is None:
        blocks = np.zeros(x.size, dtype=int)
        nb = 1
    else:
        blocks = np.asarray(block_index, dtype=int)
        if blocks.shape != x.shape:
            raise ValueError("block_index must match values shape")
        if n_blocks is None:
            raise ValueError("n_blocks is required with block_index")
        nb = int(n_blocks)
        if nb < 1 or np.any((blocks < 0) | (blocks >= nb)):
            raise ValueError("block_index contains an invalid block label")

    rho = np.full((nb, max_lag + 1), np.nan, dtype=float)
    counts = np.zeros((nb, max_lag + 1), dtype=int)

    finite_good = g & np.isfinite(x)
    for b in range(nb):
        n0 = int(np.sum(finite_good & (blocks == b)))
        counts[b, 0] = n0
        if n0 >= min_pairs:
            rho[b, 0] = 1.0

    n = x.size
    for lag in range(1, max_lag + 1):
        left = np.arange(0, n - lag, dtype=int)
        right = left + lag
        pair_good = finite_good[left] & finite_good[right]
        if block_index is None:
            pair_block = np.zeros(left.size, dtype=int)
        else:
            midpoint = (left + right) // 2
            pair_block = blocks[midpoint]

        for b in range(nb):
            use = pair_good & (pair_block == b)
            count = int(np.sum(use))
            counts[b, lag] = count
            if count < min_pairs:
                continue
            a = np.asarray(x[left[use]], dtype=float)
            c = np.asarray(x[right[use]], dtype=float)
            a = a - np.mean(a)
            c = c - np.mean(c)
            den = float(np.sqrt(np.sum(a * a) * np.sum(c * c)))
            if den > 0 and np.isfinite(den):
                rho[b, lag] = float(np.sum(a * c) / den)

    return rho, counts


def pooled_median_correlation(per_bin_rho: np.ndarray) -> np.ndarray:
    """Robustly pool per-bin lag correlations with an equal-bin median."""
    arr = np.asarray(per_bin_rho, dtype=float)
    if arr.ndim < 2:
        raise ValueError("per_bin_rho must have bin as its first axis")
    with np.errstate(all="ignore"):
        pooled = np.nanmedian(arr, axis=0)
    return np.asarray(pooled, dtype=float)


def bootstrap_simultaneous_band(
    per_bin_rho: np.ndarray,
    *,
    n_bootstrap: int,
    confidence: float = 0.95,
    random_seed: int = 12345,
) -> BootstrapBand:
    """Bootstrap PowerBins and construct a simultaneous confidence band.

    Entire PowerBins are resampled with replacement.  Individual wavelength
    pixels are never independently resampled because doing so would destroy the
    spectral correlation being measured.

    A simultaneous band is formed from the bootstrap distribution of the
    maximum absolute deviation from the pooled central curve across *all*
    inspected non-zero lags and blocks.  This protects against declaring a
    long-lag correlation significant merely because many lags were inspected.
    """
    arr = np.asarray(per_bin_rho, dtype=float)
    if arr.ndim not in {2, 3}:
        raise ValueError("per_bin_rho must have shape (nbin,nlag) or (nbin,nblock,nlag)")
    nbin = arr.shape[0]
    if nbin < 3:
        raise ValueError("At least three PowerBins are required for bootstrap confidence bands")
    n_bootstrap = int(n_bootstrap)
    if n_bootstrap < 100:
        raise ValueError("n_bootstrap must be >= 100")
    confidence = float(confidence)
    if not 0.5 < confidence < 1.0:
        raise ValueError("confidence must lie in (0.5, 1)")

    center = pooled_median_correlation(arr)
    rng = np.random.default_rng(int(random_seed))
    boot = np.full((n_bootstrap,) + center.shape, np.nan, dtype=float)
    for b in range(n_bootstrap):
        ind = rng.integers(0, nbin, size=nbin)
        boot[b] = pooled_median_correlation(arr[ind])

    # Studentization is deliberately avoided because per-lag variance estimates
    # become unstable when the true correlation is near zero.  The maximum
    # absolute deviation gives a transparent simultaneous band in rho units.
    dev = np.abs(boot - center)
    if dev.shape[-1] > 1:
        dev[..., 0] = np.nan  # lag zero is fixed by definition and not tested
    flat = dev.reshape(n_bootstrap, -1)
    with np.errstate(all="ignore"):
        max_dev = np.nanmax(flat, axis=1)
    max_dev = max_dev[np.isfinite(max_dev)]
    if max_dev.size < max(20, n_bootstrap // 2):
        raise RuntimeError("Too few finite bootstrap realizations to construct a confidence band")
    half_width = float(np.quantile(max_dev, confidence))
    lower = np.clip(center - half_width, -1.0, 1.0)
    upper = np.clip(center + half_width, -1.0, 1.0)
    if center.shape[-1] > 0:
        lower[..., 0] = 1.0
        upper[..., 0] = 1.0

    lag = np.arange(center.shape[-1], dtype=int)
    return BootstrapBand(
        lag=lag,
        center=np.asarray(center, dtype=float),
        lower=np.asarray(lower, dtype=float),
        upper=np.asarray(upper, dtype=float),
        n_bootstrap=n_bootstrap,
        confidence=confidence,
    )


def zero_insignificant_lags(rho: np.ndarray, band: BootstrapBand) -> np.ndarray:
    """Set non-zero lags to zero wherever the simultaneous band contains zero."""
    out = np.asarray(rho, dtype=float).copy()
    if out.shape != band.center.shape:
        raise ValueError("rho and bootstrap band shapes disagree")
    consistent = band.consistent_with_zero
    out[consistent] = 0.0
    out[..., 0] = 1.0
    return out


def equal_wavelength_blocks(n_pixel: int, n_blocks: int) -> np.ndarray:
    """Return deterministic equal-width wavelength-block labels."""
    n_pixel = int(n_pixel)
    n_blocks = int(n_blocks)
    if n_pixel < n_blocks or n_blocks < 1:
        raise ValueError("Require n_pixel >= n_blocks >= 1")
    # floor(j * n_blocks / n_pixel) gives blocks differing by at most one pixel.
    j = np.arange(n_pixel, dtype=int)
    return np.minimum((j * n_blocks) // n_pixel, n_blocks - 1).astype(int)


def correlation_matrix_from_lags(
    n_pixel: int,
    good: np.ndarray,
    rho_by_block: np.ndarray,
    *,
    block_index: np.ndarray | None = None,
) -> np.ndarray:
    """Construct a full-grid correlation matrix from empirical lag curves.

    Rejected pixels are intentionally statistically decoupled from fitted pixels
    and retain unit variance.  This prevents a selected whitened residual from
    depending on a masked sample through the lower-triangular operator pPXF uses.
    For wavelength-block models, a pair receives the lag coefficient of the
    block containing its midpoint.
    """
    n_pixel = int(n_pixel)
    good = np.asarray(good, dtype=bool)
    if good.shape != (n_pixel,):
        raise ValueError("good mask shape disagrees with n_pixel")
    rho = np.asarray(rho_by_block, dtype=float)
    if rho.ndim == 1:
        rho = rho[None, :]
    if rho.ndim != 2 or rho.shape[1] < 1:
        raise ValueError("rho_by_block must be one- or two-dimensional")
    if not np.allclose(rho[:, 0], 1.0, atol=1e-12, rtol=0.0):
        raise ValueError("Every lag-correlation curve must have rho(0)=1")

    nb = rho.shape[0]
    if block_index is None:
        if nb != 1:
            raise ValueError("block_index is required for multiple rho blocks")
        blocks = np.zeros(n_pixel, dtype=int)
    else:
        blocks = np.asarray(block_index, dtype=int)
        if blocks.shape != (n_pixel,):
            raise ValueError("block_index shape disagrees with n_pixel")
        if np.any((blocks < 0) | (blocks >= nb)):
            raise ValueError("block_index contains invalid labels")

    R = np.eye(n_pixel, dtype=float)
    good_idx = np.flatnonzero(good)
    if good_idx.size == 0:
        return R
    max_lag = rho.shape[1] - 1

    # Work only on the good-good subspace.  Masked dimensions remain identity and
    # have zero covariance with every good dimension.
    for a_pos, a in enumerate(good_idx):
        upper = good_idx[a_pos + 1 :]
        if upper.size == 0:
            break
        lag = upper - a
        use = lag <= max_lag
        if not np.any(use):
            continue
        b = upper[use]
        k = lag[use]
        midpoint = (a + b) // 2
        block = blocks[midpoint]
        vals = rho[block, k]
        R[a, b] = vals
        R[b, a] = vals
    return R


def _nearest_positive_definite_correlation(R: np.ndarray, eigen_floor: float) -> tuple[np.ndarray, float, bool]:
    """Eigenvalue-floor a symmetric matrix and renormalize it to unit diagonal."""
    A = np.asarray(R, dtype=float)
    A = 0.5 * (A + A.T)
    eigval, eigvec = linalg.eigh(A, check_finite=True)
    min_before = float(np.min(eigval))
    floor = float(eigen_floor)
    clipped = np.maximum(eigval, floor)
    changed = bool(np.any(eigval < floor))
    if changed:
        A = (eigvec * clipped) @ eigvec.T
        diag = np.sqrt(np.clip(np.diag(A), floor, np.inf))
        A = A / np.outer(diag, diag)
        A = 0.5 * (A + A.T)
        np.fill_diagonal(A, 1.0)
    return A, min_before, changed


def build_inverse_cholesky(
    noise: np.ndarray,
    good: np.ndarray,
    *,
    scale: float,
    rho_by_block: np.ndarray,
    block_index: np.ndarray | None = None,
    eigen_floor: float = 1.0e-8,
) -> WhiteningResult:
    """Build ``W=L^-1`` for ``C=s^2 D R D`` on the full pPXF vector.

    The science covariance is defined on fitted pixels.  Excluded pixels are
    given a benign diagonal variance and are decorrelated from fitted pixels so
    pPXF's full-vector covariance API cannot reintroduce masked samples through
    whitening.
    """
    sigma = np.asarray(noise, dtype=float).copy()
    g = np.asarray(good, dtype=bool)
    if sigma.ndim != 1 or g.shape != sigma.shape:
        raise ValueError("noise and good must be one-dimensional arrays with identical shape")
    gp = np.flatnonzero(g & np.isfinite(sigma) & (sigma > 0))
    if gp.size < 5:
        raise ValueError("Too few finite positive good-pixel uncertainties for covariance construction")
    fill = float(np.median(sigma[gp]))
    bad = ~np.isfinite(sigma) | (sigma <= 0)
    sigma[bad] = fill
    scale = float(scale)
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError("scale must be finite and positive")

    R = correlation_matrix_from_lags(
        sigma.size, g, rho_by_block, block_index=block_index
    )
    sd = scale * sigma
    C = np.outer(sd, sd) * R
    # Fast path: empirical lag matrices are often already positive definite.
    # Avoid an O(N^3) eigendecomposition unless Cholesky actually fails.
    min_before = np.nan
    changed = False
    try:
        L = linalg.cholesky(C, lower=True, check_finite=True)
    except linalg.LinAlgError:
        R_pd, min_before, changed = _nearest_positive_definite_correlation(
            R, eigen_floor=float(eigen_floor)
        )
        C = np.outer(sd, sd) * R_pd
        L = linalg.cholesky(C, lower=True, check_finite=True)
    W = linalg.solve_triangular(
        L, np.eye(C.shape[0], dtype=float), lower=True, check_finite=True
    )

    # Fitted rows must not depend on excluded residuals.  This is a critical API
    # invariant when pPXF uses a full inverse-Cholesky matrix plus goodpixels.
    excluded = ~g
    if np.any(excluded) and np.any(g):
        coupling = W[np.ix_(g, excluded)]
        if coupling.size and float(np.max(np.abs(coupling))) > 1.0e-10:
            raise RuntimeError(
                "Covariance construction coupled fitted whitened residuals to excluded pixels; "
                "masked dimensions must remain statistically decoupled."
            )

    return WhiteningResult(
        inv_cholesky=np.asarray(W, dtype=float),
        min_eigenvalue_before_regularization=min_before,
        eigenvalue_floor_applied=float(eigen_floor),
        regularized=changed,
    )


def whiten_residuals(residual: np.ndarray, inv_cholesky: np.ndarray) -> np.ndarray:
    """Apply a cached inverse-Cholesky operator to a residual vector."""
    r = np.asarray(residual, dtype=float)
    W = np.asarray(inv_cholesky, dtype=float)
    if r.ndim != 1 or W.shape != (r.size, r.size):
        raise ValueError("residual and inverse-Cholesky dimensions disagree")
    return W @ r


def covariance_total_chi2(
    residual: np.ndarray,
    inv_cholesky: np.ndarray,
    goodpixels: np.ndarray,
) -> float:
    """Return ``r^T C^-1 r`` on the fixed fitted-pixel experiment."""
    q = whiten_residuals(residual, inv_cholesky)
    gp = np.asarray(goodpixels, dtype=int)
    if gp.ndim != 1 or gp.size == 0:
        raise ValueError("goodpixels must be a non-empty one-dimensional index array")
    return float(np.sum(q[gp] ** 2))


def correlation_change(a: np.ndarray, b: np.ndarray) -> float:
    """Maximum absolute change over finite lag-correlation entries."""
    aa = np.asarray(a, dtype=float)
    bb = np.asarray(b, dtype=float)
    if aa.shape != bb.shape:
        return np.inf
    finite = np.isfinite(aa) & np.isfinite(bb)
    if not np.any(finite):
        return np.inf
    return float(np.max(np.abs(aa[finite] - bb[finite])))


def relative_scale_change(a: np.ndarray, b: np.ndarray) -> float:
    """Maximum fractional change in finite positive per-bin scale factors."""
    aa = np.asarray(a, dtype=float)
    bb = np.asarray(b, dtype=float)
    if aa.shape != bb.shape:
        return np.inf
    finite = np.isfinite(aa) & np.isfinite(bb) & (aa > 0) & (bb > 0)
    if not np.any(finite):
        return np.inf
    return float(np.max(np.abs(bb[finite] - aa[finite]) / aa[finite]))
