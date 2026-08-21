"""Diagnostic and publication plotting functions for CRD_DAP.

Every nontrivial function in this module must have a corresponding detailed
entry in ``DIAGNOSTICS.md``.  The plots are deliberately explicit about masks,
normalizations, and fitted quantities because a successful numerical exit code
is not sufficient evidence that the science result is trustworthy.
"""

from __future__ import annotations

from pathlib import Path
import json
from typing import Any

import numpy as np
import matplotlib.pyplot as plt

from .psf_lsf import ArcLSFResult, PSFEstimate
from .noise import NoiseDiagnosticResult


def save_plot_metadata(payload: dict[str, Any], path: str | Path) -> Path:
    """Save figure-generation metadata beside a diagnostic image."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n")
    return path


def _finish(fig: plt.Figure, path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_collapsed_continuum(
    image: np.ndarray,
    path: str | Path,
    *,
    title: str,
    peak_yx: tuple[float, float] | None = None,
    centroid_yx: tuple[float, float] | None = None,
) -> Path:
    """Plot one collapsed-continuum image with center diagnostics."""
    fig, ax = plt.subplots(figsize=(7, 5.5))
    finite = np.asarray(image, dtype=float)[np.isfinite(image)]
    if finite.size:
        vmin, vmax = np.nanpercentile(finite, [2, 98])
    else:
        vmin, vmax = None, None
    im = ax.imshow(image, origin="lower", aspect="equal", vmin=vmin, vmax=vmax)
    fig.colorbar(im, ax=ax, label="Collapsed continuum flux")
    if peak_yx is not None:
        y, x = peak_yx
        ax.plot(x, y, marker="x", linestyle="none", label="Smoothed peak")
    if centroid_yx is not None:
        y, x = centroid_yx
        ax.plot(x, y, marker="+", linestyle="none", label="Continuum centroid")
    if peak_yx is not None or centroid_yx is not None:
        ax.legend(loc="best")
    ax.set_xlabel("Spatial x pixel")
    ax.set_ylabel("Spatial y pixel")
    ax.set_title(title)
    return _finish(fig, path)



def plot_center_comparison(
    reference_image: np.ndarray,
    reference_wcs,
    *,
    bl_peak_sky,
    bl_centroid_sky,
    rh3_peak_sky,
    rh3_centroid_sky,
    path: str | Path,
    title: str = "BL / RH3 center comparison",
) -> Path:
    """Compare BL/RH3 peak and centroid estimates on the BL pixel grid.

    All four sky-coordinate estimates are transformed through the BL celestial
    WCS solely for a common visualization.  This plot does not alter either
    science cube and does not define the final kinematic center used by Script 4.
    """
    fig, ax = plt.subplots(figsize=(7, 5.5))
    image = np.asarray(reference_image, dtype=float)
    finite = image[np.isfinite(image)]
    if finite.size:
        vmin, vmax = np.nanpercentile(finite, [2, 98])
    else:
        vmin, vmax = None, None
    im = ax.imshow(image, origin="lower", aspect="equal", vmin=vmin, vmax=vmax)
    fig.colorbar(im, ax=ax, label="BL collapsed continuum flux")

    entries = [
        ("BL peak", bl_peak_sky, "x"),
        ("BL centroid", bl_centroid_sky, "+"),
        ("RH3 peak projected to BL", rh3_peak_sky, "o"),
        ("RH3 centroid projected to BL", rh3_centroid_sky, "s"),
    ]
    for label, sky, marker in entries:
        x, y = reference_wcs.world_to_pixel(sky)
        ax.plot(float(x), float(y), marker=marker, linestyle="none", label=label, markersize=8)

    ax.set_xlabel("BL spatial x pixel")
    ax.set_ylabel("BL spatial y pixel")
    ax.set_title(title)
    ax.legend(loc="best")
    return _finish(fig, path)

def plot_registration(
    reference_image: np.ndarray,
    moving_on_reference: np.ndarray,
    difference: np.ndarray,
    path: str | Path,
    *,
    reference_label: str = "BL",
    moving_label: str = "RH3",
    residual_shift_arcsec: tuple[float, float] | None = None,
    cross_correlation_valid: bool = True,
    status_reason: str | None = None,
    wavelength_label: str | None = None,
    reference_contrast_snr: float | None = None,
    moving_contrast_snr: float | None = None,
) -> Path:
    """Three-panel BL/RH3 WCS registration diagnostic.

    A failed morphology-contrast check is shown explicitly rather than plotting
    an unstable normalized-difference image with a misleading numerical shift.
    """
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), constrained_layout=True)
    arrays = [reference_image, moving_on_reference, difference]
    titles = [reference_label, f"{moving_label} on {reference_label} WCS", "Normalized difference"]
    for i, (ax, arr, title) in enumerate(zip(axes, arrays, titles)):
        finite = np.asarray(arr)[np.isfinite(arr)]
        if finite.size:
            lo, hi = np.nanpercentile(finite, [2, 98])
            im = ax.imshow(arr, origin="lower", aspect="equal", vmin=lo, vmax=hi)
            fig.colorbar(im, ax=ax, shrink=0.85)
        else:
            ax.imshow(np.zeros_like(reference_image, dtype=float), origin="lower", aspect="equal")
            if i == 2:
                ax.text(
                    0.5,
                    0.5,
                    "Cross-correlation\ninconclusive",
                    ha="center",
                    va="center",
                    transform=ax.transAxes,
                )
        ax.set_title(title)
        ax.set_xlabel("Spatial x pixel")
        ax.set_ylabel("Spatial y pixel")

    subtitle_parts = []
    if wavelength_label:
        subtitle_parts.append(wavelength_label)
    if reference_contrast_snr is not None and moving_contrast_snr is not None:
        subtitle_parts.append(
            f"contrast: {reference_label}={reference_contrast_snr:.2f}, "
            f"{moving_label}={moving_contrast_snr:.2f}"
        )

    if cross_correlation_valid and residual_shift_arcsec is not None:
        dy, dx = residual_shift_arcsec
        title = (
            f"Residual cross-correlation shift after WCS reprojection: "
            f"dx={dx:.3f}\", dy={dy:.3f}\""
        )
    else:
        title = status_reason or "Registration cross-correlation inconclusive"
    if subtitle_parts:
        title += "\n" + " | ".join(subtitle_parts)
    fig.suptitle(title)
    return _finish(fig, path)


def plot_valid_spaxels(
    good_spaxel: np.ndarray,
    good_fraction: np.ndarray,
    path: str | Path,
    *,
    title: str,
) -> Path:
    """Plot binary spatial validity and continuous good-wavelength fraction."""
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), constrained_layout=True)
    im0 = axes[0].imshow(good_spaxel.astype(float), origin="lower", vmin=0, vmax=1)
    fig.colorbar(im0, ax=axes[0], shrink=0.85, label="Usable (0/1)")
    axes[0].set_title("Final spatial usability")
    im1 = axes[1].imshow(good_fraction, origin="lower", vmin=0, vmax=1)
    fig.colorbar(im1, ax=axes[1], shrink=0.85, label="Good wavelength fraction")
    axes[1].set_title("Sample-level survival fraction")
    for ax in axes:
        ax.set_xlabel("Spatial x pixel")
        ax.set_ylabel("Spatial y pixel")
    fig.suptitle(title)
    return _finish(fig, path)


def plot_bad_wavelength_fraction(
    wavelength: np.ndarray,
    bad_fraction: np.ndarray,
    path: str | Path,
    *,
    title: str,
    threshold: float,
    wavegood0: float | None = None,
    wavegood1: float | None = None,
) -> Path:
    """Plot globally bad spatial-sample fraction as a function of wavelength."""
    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.plot(wavelength, bad_fraction)
    ax.axhline(threshold, linestyle="--", label=f"Global rejection threshold = {threshold:.2f}")
    if wavegood0 is not None:
        ax.axvline(float(wavegood0), linestyle=":", label="WAVGOOD0")
    if wavegood1 is not None:
        ax.axvline(float(wavegood1), linestyle=":", label="WAVGOOD1")
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel("Wavelength (Å)")
    ax.set_ylabel("Fraction of spatial samples bad")
    ax.set_title(title)
    ax.legend(loc="best")
    return _finish(fig, path)


def plot_lsf(result: ArcLSFResult, path: str | Path, *, title: str) -> Path:
    """Plot empirical arc-line FWHM measurements and supported LSF model.

    The polynomial is drawn only where accepted arc-line measurements directly
    constrain it.  Instrument-good edge regions outside that interval are
    shaded as *unconstrained extrapolation* rather than silently extending the
    fitted curve.
    """
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.scatter(
        result.wavelength,
        result.fwhm_angstrom,
        s=20,
        alpha=0.7,
        label="Accepted arc-line fits",
    )
    grid = np.linspace(
        result.measurement_wavelength_min,
        result.measurement_wavelength_max,
        500,
    )
    ax.plot(
        grid,
        result.evaluate_fwhm(grid),
        linewidth=2,
        label=f"Polynomial order {result.polynomial_order} (empirically supported)",
    )

    if result.measurement_wavelength_min > result.wavelength_min:
        ax.axvspan(
            result.wavelength_min,
            result.measurement_wavelength_min,
            alpha=0.12,
            label="Instrument-good but no accepted-line support",
        )
    if result.measurement_wavelength_max < result.wavelength_max:
        ax.axvspan(
            result.measurement_wavelength_max,
            result.wavelength_max,
            alpha=0.12,
        )
    ax.axvline(result.wavelength_min, linestyle=":", alpha=0.7)
    ax.axvline(result.wavelength_max, linestyle=":", alpha=0.7)

    ax.set_xlim(result.wavelength_min, result.wavelength_max)
    ax.set_xlabel("Wavelength (Å)")
    ax.set_ylabel("Instrumental FWHM (Å)")
    ax.set_title(title)
    ax.legend(loc="best")
    return _finish(fig, path)


def plot_lsf_spatial_variation(result: ArcLSFResult, path: str | Path, *, title: str) -> Path:
    """Plot residual line widths by slice relative to wavelength-only LSF model."""
    model = result.evaluate_fwhm(result.wavelength)
    frac = (result.fwhm_angstrom - model) / model
    fig, ax = plt.subplots(figsize=(9, 5))
    sc = ax.scatter(result.slice_id, 100.0 * frac, c=result.wavelength, s=25, alpha=0.75)
    fig.colorbar(sc, ax=ax, label="Wavelength (Å)")
    ax.axhline(0.0, linestyle="--")
    ax.set_xlabel("KCWI slice ID")
    ax.set_ylabel("FWHM residual from wavelength model (%)")
    ax.set_title(title)
    return _finish(fig, path)


def _median_and_interval_by_group(
    group: np.ndarray,
    values: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return group coordinate, median, and 16/84-percentile residuals."""
    g = np.asarray(group)
    v = np.asarray(values, dtype=float)
    coords = []
    medians = []
    low = []
    high = []
    for value in np.unique(g):
        m = (g == value) & np.isfinite(v)
        if not np.any(m):
            continue
        p16, p50, p84 = np.nanpercentile(v[m], [16, 50, 84])
        coords.append(float(value))
        medians.append(float(p50))
        low.append(float(p50 - p16))
        high.append(float(p84 - p50))
    return (
        np.asarray(coords),
        np.asarray(medians),
        np.asarray(low),
        np.asarray(high),
    )


def plot_lsf_spatial_summary(result: ArcLSFResult, path: str | Path, *, title: str) -> Path:
    """Summarize coherent LSF residuals after averaging over line measurements.

    The existing spatial-variation scatter plot intentionally shows every
    accepted line fit and therefore mixes measurement scatter with coherent
    spatial structure.  This companion figure shows group medians and 16--84%
    intervals, making the much smaller slice/position dependence visible.
    """
    model = result.evaluate_fwhm(result.wavelength)
    frac_percent = 100.0 * (result.fwhm_angstrom - model) / model

    slice_x, slice_med, slice_lo, slice_hi = _median_and_interval_by_group(
        result.slice_id, frac_percent
    )
    pos_x, pos_med, pos_lo, pos_hi = _median_and_interval_by_group(
        result.position_bin, frac_percent
    )

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
    axes[0].errorbar(
        slice_x,
        slice_med,
        yerr=np.vstack([slice_lo, slice_hi]),
        fmt="o",
        capsize=2,
    )
    axes[0].axhline(0.0, linestyle="--")
    axes[0].set_xlabel("KCWI slice ID")
    axes[0].set_ylabel("Median FWHM residual (%)")
    axes[0].set_title(f"Slice medians (RMS={100*result.slice_fractional_rms:.2f}%)")

    axes[1].errorbar(
        pos_x,
        pos_med,
        yerr=np.vstack([pos_lo, pos_hi]),
        fmt="o",
        capsize=2,
    )
    axes[1].axhline(0.0, linestyle="--")
    axes[1].set_xlabel("Position bin within slice")
    axes[1].set_ylabel("Median FWHM residual (%)")
    axes[1].set_title(
        f"Position-bin medians (RMS={100*result.position_fractional_rms:.2f}%)"
    )
    fig.suptitle(title)
    return _finish(fig, path)


def plot_psf_summary(estimate: PSFEstimate, path: str | Path, *, title: str) -> Path:
    """Record the adopted PSF estimate/provenance as a compact diagnostic figure."""
    fig, ax = plt.subplots(figsize=(7, 3.5))
    ax.axis("off")
    if np.isfinite(estimate.fwhm_arcsec):
        value = f"{estimate.fwhm_arcsec:.3f} arcsec"
    else:
        value = "Not determined"
    text = f"Adopted PSF FWHM: {value}\nSource: {estimate.source}\n{estimate.detail}"
    ax.text(0.5, 0.5, text, ha="center", va="center", transform=ax.transAxes)
    ax.set_title(title)
    return _finish(fig, path)


def plot_psf_comparison(bl: PSFEstimate, rh3: PSFEstimate, path: str | Path) -> Path:
    """Compare adopted BL and RH3 PSF values on one common angular scale."""
    fig, ax = plt.subplots(figsize=(7, 4.5))
    labels = ["BL", "RH3"]
    values = [bl.fwhm_arcsec, rh3.fwhm_arcsec]
    x = np.arange(2)
    finite = np.isfinite(values)
    if np.any(finite):
        ax.bar(x[finite], np.asarray(values)[finite])
    for xi, val in zip(x, values):
        if not np.isfinite(val):
            ax.text(xi, 0.05, "unknown", ha="center", va="bottom", transform=ax.get_xaxis_transform())
    ax.set_xticks(x, labels)
    ax.set_ylabel("PSF FWHM (arcsec)")
    ax.set_title("BL / RH3 PSF comparison")
    return _finish(fig, path)


def plot_normalized_residuals(result: NoiseDiagnosticResult, path: str | Path, *, title: str) -> Path:
    """Histogram preliminary normalized high-pass residuals for noise-scale QC."""
    z = result.normalized_residuals[np.isfinite(result.normalized_residuals)]
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.hist(z, bins=80, density=True, alpha=0.75)
    ax.set_xlim(np.nanpercentile(z, 0.5), np.nanpercentile(z, 99.5))
    ax.set_xlabel("(Flux − high-pass smooth model) / formal σ")
    ax.set_ylabel("Density")
    ax.set_title(
        f"{title}\nPreliminary variance scale factor = {result.variance_scale_factor:.3f}"
    )
    return _finish(fig, path)


def plot_spectral_covariance(result: NoiseDiagnosticResult, path: str | Path, *, title: str) -> Path:
    """Plot preliminary wavelength-lag correlation from high-pass residuals."""
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(result.lags, result.correlation, marker="o")
    ax.axhline(0.0, linestyle="--")
    ax.set_xlabel("Spectral lag (pixels)")
    ax.set_ylabel("Median residual correlation coefficient")
    ax.set_title(title)
    return _finish(fig, path)
