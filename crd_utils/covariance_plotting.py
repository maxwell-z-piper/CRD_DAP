"""Diagnostic plots for Script-3 spectral-covariance calibration.

This module is intentionally separate from :mod:`crd_utils.plotting` so the
covariance patch can be added without replacing the repository's existing,
large plotting module.  Every figure generated here is documented in
``DIAGNOSTICS.md``.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def _finish(fig: plt.Figure, path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_covariance_validation_bins_map(
    bin_map: np.ndarray,
    values: np.ndarray,
    selection_table,
    path: str | Path,
    *,
    title: str,
    colorbar_label: str,
) -> Path:
    """Overlay the deterministic covariance-validation bins on a bin-valued map."""
    bmap = np.asarray(bin_map, dtype=int)
    vals = np.asarray(values, dtype=float)
    image = np.full(bmap.shape, np.nan, dtype=float)
    valid = (bmap >= 0) & (bmap < vals.size)
    image[valid] = vals[bmap[valid]]

    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    finite = image[np.isfinite(image)]
    if finite.size:
        lo, hi = np.nanpercentile(finite, [2, 98])
    else:
        lo, hi = None, None
    im = ax.imshow(image, origin="lower", aspect="equal", vmin=lo, vmax=hi)
    fig.colorbar(im, ax=ax, label=colorbar_label)

    selected_ids = np.asarray(selection_table["BIN_ID"], dtype=int)
    selected_mask = np.isin(bmap, selected_ids)
    if np.any(selected_mask):
        ax.contour(selected_mask.astype(float), levels=[0.5], linewidths=1.5, origin="lower")

    x = np.asarray(selection_table["X_CENTROID_PIXEL"], dtype=float)
    y = np.asarray(selection_table["Y_CENTROID_PIXEL"], dtype=float)
    is_peak = np.asarray(selection_table["IS_2SIGMA_PEAK"], dtype=bool)
    added_peak = np.asarray(selection_table["ADDED_FOR_2SIGMA"], dtype=bool)
    base_radial = np.asarray(selection_table["BASE_RADIAL_SELECTION"], dtype=bool)

    for xi, yi, bid, peak, added in zip(x, y, selected_ids, is_peak, added_peak):
        marker = "*" if peak else "o"
        markersize = 11 if added else 9
        ax.plot(xi, yi, marker=marker, linestyle="none", markersize=markersize)
        ax.text(xi + 0.5, yi + 0.5, str(int(bid)), fontsize=7)

    # Use only the 12 deterministic radial selections to draw the visual PA guide.
    radial = base_radial & np.isfinite(x) & np.isfinite(y)
    if np.sum(radial) >= 2:
        xx = x[radial]
        yy = y[radial]
        if np.ptp(xx) >= np.ptp(yy) and np.ptp(xx) > 0:
            m, b = np.polyfit(xx, yy, 1)
            grid = np.linspace(np.min(xx), np.max(xx), 100)
            ax.plot(grid, m * grid + b, linestyle="--", linewidth=1.0, label="PA_kin sampling axis")
        elif np.ptp(yy) > 0:
            m, b = np.polyfit(yy, xx, 1)
            grid = np.linspace(np.min(yy), np.max(yy), 100)
            ax.plot(m * grid + b, grid, linestyle="--", linewidth=1.0, label="PA_kin sampling axis")
        ax.legend(loc="best")

    ax.set_xlabel("Spatial x pixel")
    ax.set_ylabel("Spatial y pixel")
    ax.set_title(title)
    return _finish(fig, path)


def plot_covariance_lag_band(
    lag: np.ndarray,
    center: np.ndarray,
    lower: np.ndarray,
    upper: np.ndarray,
    path: str | Path,
    *,
    title: str,
    ylabel: str = "Residual correlation coefficient",
) -> Path:
    """Plot a lag-correlation curve with its simultaneous bootstrap band."""
    lag = np.asarray(lag, dtype=float)
    c = np.asarray(center, dtype=float)
    lo = np.asarray(lower, dtype=float)
    hi = np.asarray(upper, dtype=float)
    if c.ndim == 1:
        c = c[None, :]
        lo = lo[None, :]
        hi = hi[None, :]
    fig, ax = plt.subplots(figsize=(8.5, 5.0))
    for b in range(c.shape[0]):
        label = "full fit interval" if c.shape[0] == 1 else f"wavelength block {b + 1}"
        ax.plot(lag, c[b], marker="o", label=label)
        ax.fill_between(lag, lo[b], hi[b], alpha=0.2)
    ax.axhline(0.0, linestyle="--")
    ax.set_xlabel("Spectral lag (log-wavelength pixels)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(loc="best")
    return _finish(fig, path)


def plot_covariance_scale_vs_bin_properties(
    scale: np.ndarray,
    source_table,
    path: str | Path,
    *,
    title: str,
) -> Path:
    """Plot the per-bin noise scale against PowerBin size and achieved RH3 S/N."""
    scale = np.asarray(scale, dtype=float)
    npix = (
        np.asarray(source_table["NPIX_RH3"], dtype=float)
        if "NPIX_RH3" in source_table.colnames
        else np.arange(scale.size, dtype=float)
    )
    sn = (
        np.asarray(source_table["RH3_SN"], dtype=float)
        if "RH3_SN" in source_table.colnames
        else np.full(scale.size, np.nan)
    )
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
    axes[0].scatter(npix, scale, s=18, alpha=0.7)
    axes[0].set_xlabel("RH3 member spaxels per PowerBin")
    axes[0].set_ylabel("Empirical noise scale s_i")
    axes[1].scatter(sn, scale, s=18, alpha=0.7)
    axes[1].set_xlabel("RH3 achieved S/N")
    axes[1].set_ylabel("Empirical noise scale s_i")
    fig.suptitle(title)
    return _finish(fig, path)


def plot_covariance_residual_stack(
    wavelength: np.ndarray,
    residuals: np.ndarray,
    noise: np.ndarray,
    good: np.ndarray,
    path: str | Path,
    *,
    title: str,
) -> Path:
    """Show recurring normalized residual structure versus wavelength.

    The median normalized residual across PowerBins is a direct template/LSF
    mismatch diagnostic: stochastic noise averages toward zero, while a repeated
    residual tied to a stellar feature or wavelength-local reduction problem
    remains coherent across bins.  The middle 68% PowerBin envelope is shown as
    a robust visual scale; it is diagnostic only and is not used to construct C.
    """
    wave = np.asarray(wavelength, dtype=float)
    resid = np.asarray(residuals, dtype=float)
    sig = np.asarray(noise, dtype=float)
    mask = np.asarray(good, dtype=bool)
    if resid.shape != sig.shape or resid.shape != mask.shape or resid.ndim != 2:
        raise ValueError("residuals, noise, and good must share shape (nbin, npix)")
    if wave.ndim != 1 or wave.size != resid.shape[1]:
        raise ValueError("wavelength length must match the residual spectral axis")

    z = np.full_like(resid, np.nan, dtype=float)
    ok = mask & np.isfinite(resid) & np.isfinite(sig) & (sig > 0)
    z[ok] = resid[ok] / sig[ok]
    with np.errstate(all="ignore"):
        med = np.nanmedian(z, axis=0)
        lo = np.nanpercentile(z, 16.0, axis=0)
        hi = np.nanpercentile(z, 84.0, axis=0)
        n = np.sum(np.isfinite(z), axis=0)

    fig, ax = plt.subplots(figsize=(10.0, 5.0))
    ax.plot(wave, med, linewidth=1.0, label="median normalized residual")
    ax.fill_between(wave, lo, hi, alpha=0.2, label="16--84% PowerBin envelope")
    ax.axhline(0.0, linestyle="--", linewidth=1.0)
    ax.set_xlabel("Rest-frame wavelength (Angstrom)")
    ax.set_ylabel("(data - model) / formal sigma")
    ax.set_title(title)
    ax.legend(loc="best")

    ax2 = ax.twinx()
    ax2.plot(wave, n, linewidth=0.7, alpha=0.35)
    ax2.set_ylabel("Contributing PowerBins")
    return _finish(fig, path)


def plot_covariance_model_grid_comparison(
    va_grid: np.ndarray,
    vb_grid: np.ndarray,
    fa_grid: np.ndarray,
    cubes: dict[str, np.ndarray],
    path: str | Path,
    *,
    bin_id: int,
) -> Path:
    """Compare M1--M4 one-dimensional profile Delta-chi2 curves for one bin."""
    grids = [
        (np.asarray(va_grid, dtype=float), "V_A (km/s)"),
        (np.asarray(vb_grid, dtype=float), "V_B (km/s)"),
        (np.asarray(fa_grid, dtype=float), "f_A,RH3"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(14.0, 4.2), constrained_layout=True)
    for name, cube in cubes.items():
        arr = np.asarray(cube, dtype=float)
        finite = np.isfinite(arr)
        if not np.any(finite):
            continue
        delta = arr - float(np.nanmin(arr[finite]))
        for axis, (grid, xlabel) in enumerate(grids):
            reduce_axes = tuple(a for a in range(3) if a != axis)
            profile = np.nanmin(np.where(np.isfinite(delta), delta, np.nan), axis=reduce_axes)
            axes[axis].plot(grid, profile, marker="o", markersize=3, label=name)
            axes[axis].set_xlabel(xlabel)
            axes[axis].set_ylabel("profile Delta chi2")
            axes[axis].axhline(1.0, linestyle="--", linewidth=0.7)
            axes[axis].axhline(4.0, linestyle=":", linewidth=0.7)
    axes[0].legend(loc="best")
    fig.suptitle(f"Covariance-model full-grid comparison: PowerBin {int(bin_id)}")
    return _finish(fig, path)
