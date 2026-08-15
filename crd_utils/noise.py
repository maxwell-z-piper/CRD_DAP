"""Noise-scale and spectral-covariance diagnostics.

Script 1 can only make a *preliminary* covariance assessment because no
high-quality stellar model exists yet.  The functions here therefore use
high-pass residuals from low-surface-brightness spatial samples to diagnose
whether the formal uncertainty scale and wavelength-pixel correlations are
obviously problematic.  The production covariance model should be revisited
once the first pPXF residuals are available.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.signal import savgol_filter


@dataclass(frozen=True)
class NoiseDiagnosticResult:
    normalized_residuals: np.ndarray
    variance_scale_factor: float
    lags: np.ndarray
    correlation: np.ndarray
    n_spaxels_used: int
    n_samples_used: int


def robust_variance_scale(residual: np.ndarray, sigma: np.ndarray) -> float:
    """Estimate a scalar variance rescaling from normalized residuals.

    If the robust width of ``residual/sigma`` is ``s``, multiplying the formal
    variance by ``s**2`` would bring the normalized width to approximately one.
    This is a diagnostic in Script 1, not an automatically applied correction
    unless the user explicitly enables that behavior.
    """
    r = np.asarray(residual, dtype=float)
    s = np.asarray(sigma, dtype=float)
    good = np.isfinite(r) & np.isfinite(s) & (s > 0)
    z = r[good] / s[good]
    if z.size < 10:
        raise ValueError("Too few valid residual samples to estimate variance scale")
    median = np.median(z)
    mad = np.median(np.abs(z - median))
    robust_sigma = 1.4826 * mad
    return float(robust_sigma**2)


def _choose_low_flux_spaxels(
    collapsed_image: np.ndarray,
    good_spaxel: np.ndarray,
    *,
    max_spaxels: int,
    low_flux_percentile: float,
) -> np.ndarray:
    image = np.asarray(collapsed_image, dtype=float)
    spatial_good = np.asarray(good_spaxel, dtype=bool) & np.isfinite(image)
    if not np.any(spatial_good):
        raise ValueError("No good spatial samples are available for noise diagnostics")

    values = image[spatial_good]
    threshold = np.nanpercentile(values, float(low_flux_percentile))
    candidate = np.argwhere(spatial_good & (image <= threshold))
    if candidate.size == 0:
        candidate = np.argwhere(spatial_good)

    # Deterministic selection keeps Script-1 diagnostics reproducible.  Spread
    # the selection through the ranked low-flux list rather than taking only a
    # compact patch of adjacent pixels.
    if candidate.shape[0] > int(max_spaxels):
        pick = np.linspace(0, candidate.shape[0] - 1, int(max_spaxels), dtype=int)
        candidate = candidate[pick]
    return candidate


def _safe_savgol_window(nwave: int, requested: int, polyorder: int) -> int:
    window = min(int(requested), nwave - (1 - nwave % 2))
    if window % 2 == 0:
        window -= 1
    minimum = polyorder + 3
    if minimum % 2 == 0:
        minimum += 1
    if window < minimum:
        window = minimum
    if window >= nwave:
        window = nwave - 1 if nwave % 2 == 0 else nwave
    if window % 2 == 0:
        window -= 1
    return int(window)


def preliminary_highpass_residuals(
    flux: np.ndarray,
    uncertainty: np.ndarray,
    good: np.ndarray,
    collapsed_image: np.ndarray,
    good_spaxel: np.ndarray,
    *,
    max_spaxels: int = 200,
    low_flux_percentile: float = 30.0,
    savgol_window: int = 31,
    savgol_polyorder: int = 2,
    clip_sigma: float = 5.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Build preliminary normalized high-pass residuals for selected spaxels.

    Returns
    -------
    normalized_matrix
        Shape ``(n_selected_spaxels, nwave)`` with NaN at unusable samples.
    selected_yx
        Integer ``(y, x)`` locations used for the diagnostic.
    """
    f = np.asarray(flux, dtype=float)
    s = np.asarray(uncertainty, dtype=float)
    g = np.asarray(good, dtype=bool)
    if f.shape != s.shape or f.shape != g.shape or f.ndim != 3:
        raise ValueError("flux, uncertainty, and good must be matching 3-D cubes")

    selected = _choose_low_flux_spaxels(
        collapsed_image,
        good_spaxel,
        max_spaxels=max_spaxels,
        low_flux_percentile=low_flux_percentile,
    )
    window = _safe_savgol_window(f.shape[-1], savgol_window, savgol_polyorder)
    rows: list[np.ndarray] = []
    used_coords: list[tuple[int, int]] = []

    for y, x in selected:
        spec = f[y, x].astype(float, copy=True)
        sig = s[y, x].astype(float, copy=True)
        use = g[y, x] & np.isfinite(spec) & np.isfinite(sig) & (sig > 0)
        if np.sum(use) < max(10, window):
            continue

        # Interpolate across masked points solely to obtain the smooth continuum
        # estimate.  Masked points remain NaN in the returned normalized residual.
        idx = np.arange(spec.size)
        interp = np.interp(idx, idx[use], spec[use])
        smooth = savgol_filter(interp, window_length=window, polyorder=savgol_polyorder, mode="interp")
        z = np.full(spec.size, np.nan, dtype=float)
        z[use] = (spec[use] - smooth[use]) / sig[use]

        finite = np.isfinite(z)
        if np.sum(finite) >= 10:
            med = np.nanmedian(z[finite])
            mad = np.nanmedian(np.abs(z[finite] - med))
            scale = 1.4826 * mad
            if np.isfinite(scale) and scale > 0:
                z[finite & (np.abs(z - med) > clip_sigma * scale)] = np.nan
        rows.append(z)
        used_coords.append((int(y), int(x)))

    if not rows:
        raise ValueError("No spaxels contained enough usable wavelengths for noise diagnostics")
    return np.vstack(rows), np.asarray(used_coords, dtype=int)


def estimate_spectral_correlation(normalized_matrix: np.ndarray, max_lag: int = 15) -> tuple[np.ndarray, np.ndarray]:
    """Estimate median wavelength-lag correlation from normalized residual rows."""
    z = np.asarray(normalized_matrix, dtype=float)
    if z.ndim != 2:
        raise ValueError("normalized_matrix must be 2-D: (spaxel, wavelength)")
    max_lag = int(max_lag)
    if max_lag < 0:
        raise ValueError("max_lag must be non-negative")

    correlations = []
    for lag in range(max_lag + 1):
        row_values = []
        for row in z:
            if lag == 0:
                a = row
                b = row
            else:
                a = row[:-lag]
                b = row[lag:]
            good = np.isfinite(a) & np.isfinite(b)
            if np.sum(good) < 10:
                continue
            aa = a[good] - np.mean(a[good])
            bb = b[good] - np.mean(b[good])
            denom = np.sqrt(np.sum(aa**2) * np.sum(bb**2))
            if denom > 0:
                row_values.append(float(np.sum(aa * bb) / denom))
        correlations.append(np.nan if not row_values else float(np.nanmedian(row_values)))

    lags = np.arange(max_lag + 1, dtype=int)
    corr = np.asarray(correlations, dtype=float)
    if corr.size and np.isfinite(corr[0]) and corr[0] != 0:
        corr = corr / corr[0]
    return lags, corr


def characterize_preliminary_noise(
    flux: np.ndarray,
    uncertainty: np.ndarray,
    good: np.ndarray,
    collapsed_image: np.ndarray,
    good_spaxel: np.ndarray,
    *,
    max_spaxels: int = 200,
    low_flux_percentile: float = 30.0,
    savgol_window: int = 31,
    savgol_polyorder: int = 2,
    max_lag: int = 15,
) -> NoiseDiagnosticResult:
    """Run the full Script-1 preliminary noise/covariance diagnostic."""
    z, selected = preliminary_highpass_residuals(
        flux,
        uncertainty,
        good,
        collapsed_image,
        good_spaxel,
        max_spaxels=max_spaxels,
        low_flux_percentile=low_flux_percentile,
        savgol_window=savgol_window,
        savgol_polyorder=savgol_polyorder,
    )
    finite = np.isfinite(z)
    vals = z[finite]
    if vals.size < 10:
        raise ValueError("Too few normalized residual samples for noise characterization")
    med = np.median(vals)
    mad = np.median(np.abs(vals - med))
    width = 1.4826 * mad
    variance_scale = float(width**2)
    lags, corr = estimate_spectral_correlation(z, max_lag=max_lag)
    return NoiseDiagnosticResult(
        normalized_residuals=z,
        variance_scale_factor=variance_scale,
        lags=lags,
        correlation=corr,
        n_spaxels_used=int(z.shape[0]),
        n_samples_used=int(np.sum(finite)),
    )


def save_noise_diagnostic(result: NoiseDiagnosticResult, path: str) -> None:
    """Save the numerical Script-1 noise diagnostic as a compressed NPZ file."""
    np.savez_compressed(
        path,
        normalized_residuals=result.normalized_residuals,
        variance_scale_factor=result.variance_scale_factor,
        lags=result.lags,
        correlation=result.correlation,
        n_spaxels_used=result.n_spaxels_used,
        n_samples_used=result.n_samples_used,
    )


def generate_correlated_noise(covariance: np.ndarray, *, rng: np.random.Generator | None = None) -> np.ndarray:
    """Draw one zero-mean multivariate Gaussian noise realization."""
    cov = np.asarray(covariance, dtype=float)
    if cov.ndim != 2 or cov.shape[0] != cov.shape[1]:
        raise ValueError("covariance must be square")
    rng = np.random.default_rng() if rng is None else rng
    return rng.multivariate_normal(np.zeros(cov.shape[0]), cov)
