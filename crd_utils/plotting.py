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


def _registration_panel_limits(
    image: np.ndarray,
    mask: np.ndarray,
    *,
    lower_percentile: float = 5.0,
    upper_percentile: float = 99.5,
) -> tuple[float | None, float | None]:
    """Return robust display limits for one registration-image panel.

    Registration images can contain a small number of very large edge values
    from masked/reprojected IFU boundaries.  Those samples are important to the
    data-quality bookkeeping, but they should not control the color stretch of
    a morphology-comparison figure.  Limits are therefore measured only from
    the *common valid footprint* and use deliberately robust percentiles.

    This helper affects visualization only.  It does not change the images,
    masks, WCS, or cross-correlation calculation used for registration QC.
    """
    arr = np.asarray(image, dtype=float)
    use = np.asarray(mask, dtype=bool) & np.isfinite(arr)
    values = arr[use]
    if values.size == 0:
        return None, None

    lo, hi = np.nanpercentile(values, [lower_percentile, upper_percentile])
    if not np.isfinite(lo) or not np.isfinite(hi):
        return None, None
    if hi <= lo:
        lo = float(np.nanmin(values))
        hi = float(np.nanmax(values))
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        return None, None
    return float(lo), float(hi)


def plot_registration(
    reference_image: np.ndarray,
    moving_on_reference: np.ndarray,
    difference: np.ndarray,
    path: str | Path,
    *,
    reference_label: str = "BL",
    moving_label: str = "RH3",
    overlap: np.ndarray | None = None,
    residual_shift_arcsec: tuple[float, float] | None = None,
    cross_correlation_valid: bool = True,
    status_reason: str | None = None,
    wavelength_label: str | None = None,
    reference_contrast_snr: float | None = None,
    moving_contrast_snr: float | None = None,
) -> Path:
    """Three-panel BL/RH3 WCS-registration diagnostic.

    The first two panels are intentionally displayed only over the common valid
    footprint and use robust stretches measured from that same footprint.  This
    keeps thin IFU/reprojection edge artifacts from dominating the color scale
    and hiding the galaxy morphology that the diagnostic is meant to compare.

    The third panel is conditional:

    * when morphology cross-correlation is valid, show the normalized morphology
      difference with a symmetric zero-centered stretch;
    * when cross-correlation is inconclusive, do **not** show a numerically
      unstable residual image.  Instead show the common valid footprint and an
      explicit annotation that no morphology residual/shift was adopted.

    These choices affect the diagnostic visualization only.  They do not alter
    the science cubes or registration decision itself.
    """
    reference = np.asarray(reference_image, dtype=float)
    moving = np.asarray(moving_on_reference, dtype=float)
    diff = np.asarray(difference, dtype=float)

    if reference.shape != moving.shape:
        raise ValueError("reference_image and moving_on_reference must have the same shape")

    if overlap is None:
        common = np.isfinite(reference) & np.isfinite(moving)
    else:
        common = np.asarray(overlap, dtype=bool)
        if common.shape != reference.shape:
            raise ValueError("overlap must match the registration image shape")
        common &= np.isfinite(reference) & np.isfinite(moving)

    # Display the two morphology panels on exactly the same spatial support.
    # Outside the common valid footprint there is nothing meaningful to compare.
    reference_display = np.where(common, reference, np.nan)
    moving_display = np.where(common, moving, np.nan)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), constrained_layout=True)

    for ax, arr, title in (
        (axes[0], reference_display, f"{reference_label} (common valid footprint)"),
        (axes[1], moving_display, f"{moving_label} on {reference_label} WCS"),
    ):
        vmin, vmax = _registration_panel_limits(arr, common)
        im = ax.imshow(
            arr,
            origin="lower",
            aspect="equal",
            vmin=vmin,
            vmax=vmax,
        )
        fig.colorbar(im, ax=ax, shrink=0.85)
        ax.set_title(title)
        ax.set_xlabel("Spatial x pixel")
        ax.set_ylabel("Spatial y pixel")

    if cross_correlation_valid:
        # Only a trusted registration gets a morphology-residual panel.  Use a
        # symmetric color range around zero so positive/negative residuals have
        # equal visual weight and a single outlier cannot set the entire scale.
        diff_display = np.where(common & np.isfinite(diff), diff, np.nan)
        finite = diff_display[np.isfinite(diff_display)]
        if finite.size:
            vmax = float(np.nanpercentile(np.abs(finite), 98.0))
            if not np.isfinite(vmax) or vmax <= 0:
                vmax = None
                vmin = None
            else:
                vmin = -vmax
        else:
            vmin = vmax = None

        im = axes[2].imshow(
            diff_display,
            origin="lower",
            aspect="equal",
            vmin=vmin,
            vmax=vmax,
        )
        fig.colorbar(im, ax=axes[2], shrink=0.85, label="Normalized morphology residual")
        axes[2].set_title("Normalized morphology difference")
    else:
        # A failed/inconclusive registration should never produce a dramatic
        # residual plot.  Such a plot has no accepted astrometric interpretation
        # and can become pathological when one arm has near-zero continuum.
        footprint = common.astype(float)
        im = axes[2].imshow(
            footprint,
            origin="lower",
            aspect="equal",
            vmin=0.0,
            vmax=1.0,
        )
        fig.colorbar(im, ax=axes[2], shrink=0.85, label="Common valid footprint (0/1)")
        axes[2].set_title("Registration inconclusive")
        axes[2].text(
            0.5,
            0.5,
            "No morphology residual shown\nNo numerical shift adopted",
            ha="center",
            va="center",
            transform=axes[2].transAxes,
            bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.85},
        )

    axes[2].set_xlabel("Spatial x pixel")
    axes[2].set_ylabel("Spatial y pixel")

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
            "Residual cross-correlation shift after WCS reprojection: "
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



def plot_effective_exposure(
    median_exposure: np.ndarray,
    coverage_fraction: np.ndarray,
    path: str | Path,
    *,
    title: str,
) -> Path:
    """Plot median effective exposure and instrument-good coverage fraction.

    KcwiKit output grids may intentionally be larger than the illuminated IFU
    footprint.  This plot keeps zero-exposure padding visually distinct from
    genuinely observed low-S/N spaxels so downstream binning never confuses the
    two.
    """
    exp = np.asarray(median_exposure, dtype=float)
    cov = np.asarray(coverage_fraction, dtype=float)
    if exp.shape != cov.shape or exp.ndim != 2:
        raise ValueError(
            "median_exposure and coverage_fraction must be matching 2-D arrays"
        )

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), constrained_layout=True)

    finite_positive = exp[np.isfinite(exp) & (exp > 0)]
    vmax = (
        float(np.nanpercentile(finite_positive, 99))
        if finite_positive.size
        else None
    )
    im0 = axes[0].imshow(exp, origin="lower", vmin=0, vmax=vmax, aspect="equal")
    fig.colorbar(
        im0, ax=axes[0], shrink=0.85, label="Median effective exposure (s)"
    )
    axes[0].set_title("Effective exposure")

    im1 = axes[1].imshow(cov, origin="lower", vmin=0, vmax=1, aspect="equal")
    fig.colorbar(
        im1, ax=axes[1], shrink=0.85, label="Wavelength coverage fraction"
    )
    axes[1].set_title("Instrument-good wavelength coverage")

    for ax in axes:
        ax.set_xlabel("Spatial x pixel")
        ax.set_ylabel("Spatial y pixel")

    fig.suptitle(title)
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

# =============================================================================
# Script 2: BL master PowerBin diagnostics
# =============================================================================

def _finite_spatial_values(values: np.ndarray, mask: np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    use = np.asarray(mask, dtype=bool) & np.isfinite(arr)
    return arr[use]


def plot_binning_aperture(
    significance_proxy: np.ndarray,
    aperture_mask: np.ndarray,
    x_arcsec: np.ndarray,
    y_arcsec: np.ndarray,
    path: str | Path,
    *,
    threshold: float,
    center_xy_arcsec: tuple[float, float] = (0.0, 0.0),
    title: str = "BL stellar-body binning aperture",
) -> Path:
    """Show the smoothed detection proxy and final PowerBin input aperture."""
    proxy = np.asarray(significance_proxy, dtype=float)
    aperture = np.asarray(aperture_mask, dtype=bool)
    x = np.asarray(x_arcsec, dtype=float)
    y = np.asarray(y_arcsec, dtype=float)
    if proxy.shape != aperture.shape or proxy.shape != x.shape or proxy.shape != y.shape:
        raise ValueError("Script-2 aperture plot arrays must have matching 2-D shapes")

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.6), constrained_layout=True)
    finite = proxy[np.isfinite(proxy)]
    if finite.size:
        vmin, vmax = np.nanpercentile(finite, [2, 98])
    else:
        vmin = vmax = None
    im = axes[0].imshow(proxy, origin="lower", aspect="equal", vmin=vmin, vmax=vmax)
    fig.colorbar(im, ax=axes[0], shrink=0.85, label="Smoothed continuum significance proxy")
    axes[0].contour(aperture.astype(float), levels=[0.5], linewidths=1.0)
    axes[0].set_title(f"Detection proxy | threshold={float(threshold):.2f}")
    axes[0].set_xlabel("Spatial x pixel")
    axes[0].set_ylabel("Spatial y pixel")

    yy, xx = np.where(aperture)
    if yy.size:
        sc = axes[1].scatter(x[yy, xx], y[yy, xx], c=proxy[yy, xx], s=22, marker="s")
        fig.colorbar(sc, ax=axes[1], shrink=0.85, label="Detection proxy")
    axes[1].plot(center_xy_arcsec[0], center_xy_arcsec[1], marker="+", linestyle="none", markersize=10)
    axes[1].set_aspect("equal")
    axes[1].set_xlabel("Tangent-plane x (arcsec)")
    axes[1].set_ylabel("Tangent-plane y (arcsec)")
    axes[1].set_title(f"Final aperture | N={int(np.sum(aperture))} pixels")
    fig.suptitle(title)
    return _finish(fig, path)


def plot_master_bins(
    bin_map: np.ndarray,
    x_arcsec: np.ndarray,
    y_arcsec: np.ndarray,
    path: str | Path,
    *,
    center_xy_arcsec: tuple[float, float] = (0.0, 0.0),
    pa_kin_deg: float | None = None,
    title: str = "BL master PowerBins",
) -> Path:
    """Plot the physical BL-defined PowerBin membership in tangent-plane coordinates."""
    bin_map = np.asarray(bin_map, dtype=int)
    x = np.asarray(x_arcsec, dtype=float)
    y = np.asarray(y_arcsec, dtype=float)
    if bin_map.shape != x.shape or bin_map.shape != y.shape:
        raise ValueError("bin_map/x_arcsec/y_arcsec must match")
    use = bin_map >= 0
    fig, ax = plt.subplots(figsize=(7, 6))
    yy, xx = np.where(use)
    if yy.size:
        sc = ax.scatter(x[yy, xx], y[yy, xx], c=bin_map[yy, xx], s=24, marker="s")
        fig.colorbar(sc, ax=ax, label="PowerBin ID")
    ax.plot(center_xy_arcsec[0], center_xy_arcsec[1], marker="+", linestyle="none", markersize=12, label="Adopted center")
    if pa_kin_deg is not None and np.isfinite(pa_kin_deg):
        # Tangent-plane x/y are east/north-like offsets. Astronomical PA is
        # measured east of north, so dx=sin(PA), dy=cos(PA).
        pa = np.deg2rad(float(pa_kin_deg))
        finite_radius = np.hypot(x[use], y[use]) if np.any(use) else np.array([1.0])
        length = float(np.nanmax(finite_radius)) if finite_radius.size else 1.0
        dx = np.sin(pa) * length
        dy = np.cos(pa) * length
        ax.plot([-dx, dx], [-dy, dy], linestyle="--", label=r"$PA_{\rm kin}$ axis")
    ax.set_aspect("equal")
    ax.set_xlabel("Tangent-plane x (arcsec)")
    ax.set_ylabel("Tangent-plane y (arcsec)")
    ax.set_title(title)
    ax.legend(loc="best")
    return _finish(fig, path)


def _bin_values_on_map(bin_map: np.ndarray, values: np.ndarray) -> np.ndarray:
    bin_map = np.asarray(bin_map, dtype=int)
    vals = np.asarray(values, dtype=float)
    out = np.full(bin_map.shape, np.nan, dtype=float)
    use = (bin_map >= 0) & (bin_map < vals.size)
    out[use] = vals[bin_map[use]]
    return out


def plot_bin_value_map(
    bin_map: np.ndarray,
    values: np.ndarray,
    path: str | Path,
    *,
    title: str,
    colorbar_label: str,
    vmin: float | None = None,
    vmax: float | None = None,
) -> Path:
    """Plot one scalar per PowerBin on the BL membership grid."""
    image = _bin_values_on_map(bin_map, values)
    fig, ax = plt.subplots(figsize=(7, 5.5))
    im = ax.imshow(image, origin="lower", aspect="equal", vmin=vmin, vmax=vmax)
    fig.colorbar(im, ax=ax, label=colorbar_label)
    ax.set_xlabel("BL spatial x pixel")
    ax.set_ylabel("BL spatial y pixel")
    ax.set_title(title)
    return _finish(fig, path)


def plot_bl_rh3_sn_comparison(
    bin_map: np.ndarray,
    bl_sn: np.ndarray,
    rh3_sn: np.ndarray,
    path: str | Path,
    *,
    upper_percentile: float = 95.0,
) -> Path:
    """Compare BL and red/RH3 S/N on one identical PowerBin color scale."""
    bl = np.asarray(bl_sn, dtype=float)
    rh = np.asarray(rh3_sn, dtype=float)
    combined = np.concatenate([bl[np.isfinite(bl)], rh[np.isfinite(rh)]])
    if combined.size:
        vmax = float(np.nanpercentile(combined, float(upper_percentile)))
        if not np.isfinite(vmax) or vmax <= 0:
            vmax = None
    else:
        vmax = None
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.6), constrained_layout=True)
    images = [
        _bin_values_on_map(bin_map, bl),
        _bin_values_on_map(bin_map, rh),
    ]
    for ax, image, title in zip(axes, images, ("BL", "RH3 in BL-defined bins")):
        im = ax.imshow(image, origin="lower", aspect="equal", vmin=0.0, vmax=vmax)
        fig.colorbar(im, ax=ax, shrink=0.85, label="Robust continuum S/N per spectral pixel")
        ax.set_xlabel("BL spatial x pixel")
        ax.set_ylabel("BL spatial y pixel")
        ax.set_title(title)
    fig.suptitle(
        f"BL / RH3 achieved S/N | shared vmax=P{float(upper_percentile):g}(combined)"
    )
    return _finish(fig, path)


def plot_bin_transfer(
    bl_bin_map: np.ndarray,
    rh3_bin_map: np.ndarray,
    match_distance_arcsec: np.ndarray,
    path: str | Path,
    *,
    max_distance_arcsec: float,
) -> Path:
    """QC the sky-WCS transfer of BL membership onto the RH3 native grid."""
    bl = np.asarray(bl_bin_map, dtype=int)
    rh = np.asarray(rh3_bin_map, dtype=int)
    dist = np.asarray(match_distance_arcsec, dtype=float)
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.4), constrained_layout=True)
    for ax, image, title in (
        (axes[0], np.where(bl >= 0, bl, np.nan), "BL master bin map"),
        (axes[1], np.where(rh >= 0, rh, np.nan), "RH3 transferred membership"),
    ):
        im = ax.imshow(image, origin="lower", aspect="equal")
        fig.colorbar(im, ax=ax, shrink=0.82, label="PowerBin ID")
        ax.set_title(title)
        ax.set_xlabel("Spatial x pixel")
        ax.set_ylabel("Spatial y pixel")
    dshow = np.where(rh >= 0, dist, np.nan)
    im = axes[2].imshow(dshow, origin="lower", aspect="equal", vmin=0.0, vmax=float(max_distance_arcsec))
    fig.colorbar(im, ax=axes[2], shrink=0.82, label="Nearest BL pixel-center distance (arcsec)")
    axes[2].set_title("WCS membership-transfer distance")
    axes[2].set_xlabel("RH3 spatial x pixel")
    axes[2].set_ylabel("RH3 spatial y pixel")
    return _finish(fig, path)
