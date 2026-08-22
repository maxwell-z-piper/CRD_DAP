"""BL-defined master spatial binning utilities for CRD_DAP Script 2.

Script 2 creates exactly one adaptive spatial tessellation from the BL cube and
then transfers that *physical sky membership* to the red/RH3 stream.  This
module contains the scientific logic behind that stage so the driver remains
short and every assumption can be unit-tested independently.

Key design choices
------------------
* PowerBin is mandatory.  CRD_DAP never silently falls back to legacy VorBin.
* The PowerBin capacity is ``(S/N)^2``.  A callable capacity is used even for
  the current diagonal-noise baseline so a calibrated non-additive spatial
  covariance correction can be inserted later without changing the algorithm.
* Low native-spaxel S/N is *not* by itself a rejection criterion.  A smoothed
  continuum-detection aperture identifies the useful stellar body, after which
  PowerBin is allowed to combine faint spatial samples to the target S/N.
* Bin spectra are geometric aperture sums, not inverse-variance-weighted spatial
  averages.  This preserves the physical light weighting required by later
  bin-integrated disk models.
* Formal variances are propagated diagonally.  Script 1's spectral-correlation
  diagnostic is carried forward as provenance, but it is not converted into an
  unvalidated covariance matrix here.
"""

from __future__ import annotations

from dataclasses import dataclass
from importlib import metadata as importlib_metadata
from typing import Any

import numpy as np
from scipy.ndimage import binary_dilation, convolve, label


@dataclass(frozen=True)
class ContinuumWindowMaps:
    """Per-spaxel continuum signal/noise summaries for one wavelength window."""

    signal: np.ndarray
    noise: np.ndarray
    valid_fraction: np.ndarray
    wavelength_mask: np.ndarray
    observed_min: float
    observed_max: float
    n_channels: int


@dataclass(frozen=True)
class SNWindowCoverage:
    """Coverage diagnostics for one configured continuum-S/N window.

    This is a pure QC object.  It never changes the requested wavelength range
    or substitutes another range.  ``requested_*`` describe the full
    redshifted configuration interval; ``usable_*`` describe the portion that
    survives Script-1 ``GOODWAVE``.
    """

    requested_observed_min: float
    requested_observed_max: float
    usable_observed_min: float
    usable_observed_max: float
    requested_width_angstrom: float
    envelope_coverage_fraction: float
    usable_channel_fraction: float
    n_requested_channels: int
    n_usable_channels: int
    truncated_blue: bool
    truncated_red: bool


@dataclass(frozen=True)
class ApertureResult:
    """Smoothed galaxy-detection aperture used as the PowerBin input domain."""

    mask: np.ndarray
    significance_proxy: np.ndarray
    threshold: float
    n_pixels: int
    nearest_component_distance_pix: float


@dataclass(frozen=True)
class PowerBinResult:
    """Compact CRD_DAP representation of a PowerBin tessellation."""

    bin_map: np.ndarray
    aperture_y: np.ndarray
    aperture_x: np.ndarray
    bin_num: np.ndarray
    n_bins: int
    bin_capacity: np.ndarray
    generator_xy: np.ndarray
    generator_radius: np.ndarray
    npix_per_bin: np.ndarray
    single: np.ndarray
    rms_frac_percent: float
    powerbin_version: str


@dataclass(frozen=True)
class TransferResult:
    """Mapping of red/RH3 spaxels onto BL-defined physical PowerBins."""

    bin_map: np.ndarray
    match_distance_arcsec: np.ndarray
    n_candidate_spaxels: int
    n_assigned_spaxels: int
    assigned_fraction: float


@dataclass(frozen=True)
class CoaddedBinSpectra:
    """Aperture-summed spectra and formal diagonal uncertainties for all bins."""

    flux: np.ndarray
    uncertainty: np.ndarray
    good: np.ndarray
    contributing_spaxels: np.ndarray
    n_members: np.ndarray
    spatial_scale_factor: float
    spatial_scale_reason: str


@dataclass(frozen=True)
class BinSNDiagnostics:
    """Robust achieved-continuum S/N diagnostics for every spatial bin.

    ``sn`` is the production-facing achieved S/N value.  By default it is only
    defined when the median continuum signal is positive.  ``signed_sn`` keeps
    the signed ratio-of-medians for audit/debugging, while
    ``legacy_median_ratio`` records the older ``median(flux/uncertainty)``
    estimator.  Keeping both makes numerical pathologies visible without
    allowing a huge negative value to masquerade as a meaningful achieved S/N.
    """

    sn: np.ndarray
    signed_sn: np.ndarray
    legacy_median_ratio: np.ndarray
    median_flux: np.ndarray
    median_uncertainty: np.ndarray
    min_uncertainty: np.ndarray
    p05_uncertainty: np.ndarray
    negative_flux_fraction: np.ndarray
    n_good_channels: np.ndarray
    positive_continuum: np.ndarray


def observed_range_from_rest(rest_range: tuple[float, float], redshift: float) -> tuple[float, float]:
    """Convert a rest-frame wavelength interval to observed-frame Angstrom."""
    lo, hi = (float(rest_range[0]), float(rest_range[1]))
    if not np.isfinite(lo) or not np.isfinite(hi) or lo >= hi:
        raise ValueError("rest_range must contain finite increasing wavelengths")
    if redshift < 0:
        raise ValueError("redshift must be non-negative")
    factor = 1.0 + float(redshift)
    return lo * factor, hi * factor


def sn_window_coverage(
    wavelength: np.ndarray,
    good_wavelength: np.ndarray,
    *,
    rest_range: tuple[float, float],
    redshift: float,
) -> SNWindowCoverage:
    """Quantify whether a configured S/N window is fully supported by the cube.

    The diagnostic deliberately does **not** search for a better wavelength
    range and does not alter the requested interval.  Its sole purpose is to
    catch configuration/data mismatches such as using a production CaT window
    with an RL integration-test cube whose usable red edge truncates that
    interval.

    Two complementary fractions are reported:

    ``envelope_coverage_fraction``
        Fraction of the requested wavelength *width* lying inside the envelope
        spanned by Script-1 ``GOODWAVE`` channels.

    ``usable_channel_fraction``
        Fraction of native wavelength samples inside the requested interval
        that are actually marked good by Script 1.  Isolated masked channels
        therefore affect this number without being mistaken for an edge
        truncation.
    """
    wave = np.asarray(wavelength, dtype=float).ravel()
    good = np.asarray(good_wavelength, dtype=bool).ravel()
    if wave.size != good.size:
        raise ValueError("wavelength and good_wavelength must have the same length")
    finite_wave = np.isfinite(wave)
    good = good & finite_wave
    if not np.any(good):
        raise ValueError("No finite Script-1 GOODWAVE channels are available")

    req_lo, req_hi = observed_range_from_rest(rest_range, redshift)
    req_width = float(req_hi - req_lo)
    if req_width <= 0:
        raise ValueError("Requested observed S/N window has non-positive width")

    good_wave = wave[good]
    env_lo = float(np.nanmin(good_wave))
    env_hi = float(np.nanmax(good_wave))
    overlap_lo = max(req_lo, env_lo)
    overlap_hi = min(req_hi, env_hi)
    overlap_width = max(0.0, overlap_hi - overlap_lo)
    envelope_fraction = float(np.clip(overlap_width / req_width, 0.0, 1.0))

    requested_channels = finite_wave & (wave >= req_lo) & (wave <= req_hi)
    n_requested = int(np.sum(requested_channels))
    n_usable = int(np.sum(requested_channels & good))
    usable_channel_fraction = (
        float(n_usable / n_requested) if n_requested > 0 else 0.0
    )

    # Allow half a native wavelength sample when deciding whether an edge is
    # genuinely truncated, so floating-point WCS rounding does not create a
    # spurious warning for an otherwise exact match.
    finite_sorted = np.sort(wave[finite_wave])
    if finite_sorted.size > 1:
        dw = float(np.nanmedian(np.diff(finite_sorted)))
        tolerance = 0.5 * abs(dw) if np.isfinite(dw) else 0.0
    else:
        tolerance = 0.0

    return SNWindowCoverage(
        requested_observed_min=float(req_lo),
        requested_observed_max=float(req_hi),
        usable_observed_min=float(max(req_lo, env_lo)) if overlap_width > 0 else np.nan,
        usable_observed_max=float(min(req_hi, env_hi)) if overlap_width > 0 else np.nan,
        requested_width_angstrom=req_width,
        envelope_coverage_fraction=envelope_fraction,
        usable_channel_fraction=usable_channel_fraction,
        n_requested_channels=n_requested,
        n_usable_channels=n_usable,
        truncated_blue=bool(req_lo < env_lo - tolerance),
        truncated_red=bool(req_hi > env_hi + tolerance),
    )


def _safe_nanmedian(array: np.ndarray, axis: int) -> np.ndarray:
    """Nan-median without emitting all-NaN slice warnings."""
    arr = np.asarray(array, dtype=float)
    finite_count = np.sum(np.isfinite(arr), axis=axis)
    # np.nanmedian is still the most efficient implementation; suppress only the
    # expected RuntimeWarning for deliberately all-masked spatial elements.
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        out = np.nanmedian(arr, axis=axis)
    out = np.asarray(out, dtype=float)
    out[finite_count == 0] = np.nan
    return out


def continuum_window_maps(
    cube,
    *,
    rest_range: tuple[float, float],
    redshift: float,
    min_valid_fraction: float = 0.5,
) -> ContinuumWindowMaps:
    """Build robust per-spaxel signal/noise estimates in a science window.

    ``signal`` is the median flux-density sample and ``noise`` is the square root
    of the median formal variance over the usable channels in the requested
    rest-frame interval.  With equal spatial pixel area these scalars preserve
    the expected per-spectral-pixel aperture S/N scaling:

        S/N_bin = sum(signal_p) / sqrt(sum(noise_p^2)).

    This is intentionally a continuum-S/N proxy rather than a broad-band
    integrated S/N, so the configured target remains interpretable as an
    approximate per-spectral-pixel S/N in the fitting region.
    """
    min_valid_fraction = float(min_valid_fraction)
    if not 0.0 < min_valid_fraction <= 1.0:
        raise ValueError("min_valid_fraction must lie in (0, 1]")

    obs_lo, obs_hi = observed_range_from_rest(rest_range, redshift)
    wave = np.asarray(cube.wavelength, dtype=float)
    wave_ok = np.asarray(cube.good_wavelength, dtype=bool)
    wmask = wave_ok & (wave >= obs_lo) & (wave <= obs_hi)
    nchan = int(np.sum(wmask))
    if nchan < 2:
        raise ValueError(
            f"Requested rest window {rest_range} has only {nchan} usable observed channels "
            f"({obs_lo:.1f}--{obs_hi:.1f} A)"
        )

    good = np.asarray(cube.good[..., wmask], dtype=bool)
    flux = np.asarray(cube.flux[..., wmask], dtype=float)
    var = np.asarray(cube.variance[..., wmask], dtype=float)
    sample_good = good & np.isfinite(flux) & np.isfinite(var) & (var > 0)

    fraction = np.mean(sample_good, axis=-1)
    flux_masked = np.where(sample_good, flux, np.nan)
    var_masked = np.where(sample_good, var, np.nan)
    signal = _safe_nanmedian(flux_masked, axis=-1)
    noise = np.sqrt(_safe_nanmedian(var_masked, axis=-1))

    spatial_ok = (
        np.asarray(cube.good_spaxel, dtype=bool)
        & (fraction >= min_valid_fraction)
        & np.isfinite(signal)
        & np.isfinite(noise)
        & (noise > 0)
    )
    signal = np.where(spatial_ok, signal, np.nan)
    noise = np.where(spatial_ok, noise, np.nan)

    used_wave = wave[wmask]
    return ContinuumWindowMaps(
        signal=signal,
        noise=noise,
        valid_fraction=fraction,
        wavelength_mask=wmask,
        observed_min=float(np.nanmin(used_wave)),
        observed_max=float(np.nanmax(used_wave)),
        n_channels=nchan,
    )


def _smoothed_significance(
    signal: np.ndarray,
    noise: np.ndarray,
    valid: np.ndarray,
    sigma: float,
) -> np.ndarray:
    """Return a Gaussian-weighted local continuum significance map.

    This is used only to define the outer stellar-body aperture.  For a Gaussian
    kernel ``k`` and diagonal spatial variances, normalization of the kernel
    cancels in S/N, giving

        S/N = sum(k_p S_p) / sqrt(sum(k_p^2 N_p^2)).

    Missing pixels are excluded from both sums.  The calculation therefore gains
    significance for spatially extended low-surface-brightness light rather than
    merely smoothing the *native* per-pixel S/N, which would not reflect the
    information gained by combining neighboring pixels.  KcwiKit spatial
    covariance is not known yet, so this map is explicitly an aperture-detection
    proxy rather than a likelihood quantity.
    """
    signal = np.asarray(signal, dtype=float)
    noise = np.asarray(noise, dtype=float)
    valid = (
        np.asarray(valid, dtype=bool)
        & np.isfinite(signal)
        & np.isfinite(noise)
        & (noise > 0)
    )
    out = np.full(signal.shape, np.nan, dtype=float)
    if float(sigma) <= 0:
        out[valid] = signal[valid] / noise[valid]
        return out

    radius = max(1, int(np.ceil(4.0 * float(sigma))))
    yy, xx = np.mgrid[-radius : radius + 1, -radius : radius + 1]
    kernel = np.exp(-(xx**2 + yy**2) / (2.0 * float(sigma) ** 2))
    kernel /= np.sum(kernel)

    weighted_signal = convolve(
        np.where(valid, signal, 0.0),
        kernel,
        mode="constant",
        cval=0.0,
    )
    weighted_var = convolve(
        np.where(valid, noise**2, 0.0),
        kernel**2,
        mode="constant",
        cval=0.0,
    )
    kernel_coverage = convolve(
        valid.astype(float),
        kernel,
        mode="constant",
        cval=0.0,
    )

    ok = valid & (kernel_coverage > 1e-6) & (weighted_var > 0)
    out[ok] = weighted_signal[ok] / np.sqrt(weighted_var[ok])
    return out


def make_analysis_aperture(
    signal: np.ndarray,
    noise: np.ndarray,
    good_spaxel: np.ndarray,
    *,
    center_yx: tuple[float, float],
    mode: str = "auto_connected_sn",
    smooth_sigma_pix: float = 2.0,
    threshold: float = 1.0,
    dilate_pix: int = 2,
    center_max_distance_pix: float = 5.0,
    min_pixels: int = 25,
    max_radius_arcsec: float | None = None,
    x_arcsec: np.ndarray | None = None,
    y_arcsec: np.ndarray | None = None,
) -> ApertureResult:
    """Define the useful stellar-body aperture without rejecting faint bins.

    The automatic mode thresholds a *smoothed continuum significance proxy* and
    keeps the connected component nearest the independently measured galaxy
    center.  The selected component can be dilated before intersection with the
    Script-1 hard-good spaxel mask.  Once a spaxel lies inside this aperture,
    its native S/N is not used as a rejection criterion; PowerBin may combine it
    with neighbors to reach the target S/N.

    The significance image is only an aperture-definition diagnostic.  It is not
    propagated as a formal likelihood quantity.
    """
    signal = np.asarray(signal, dtype=float)
    noise = np.asarray(noise, dtype=float)
    good_spaxel = np.asarray(good_spaxel, dtype=bool)
    if signal.shape != noise.shape or signal.shape != good_spaxel.shape or signal.ndim != 2:
        raise ValueError("signal, noise, and good_spaxel must be matching 2-D arrays")

    valid = good_spaxel & np.isfinite(signal) & np.isfinite(noise) & (noise > 0)
    proxy = _smoothed_significance(
        signal,
        noise,
        valid,
        float(smooth_sigma_pix),
    )

    mode = str(mode).strip().lower()
    if mode == "all_good":
        aperture = valid.copy()
        nearest_distance = 0.0
    elif mode == "auto_connected_sn":
        detected = valid & np.isfinite(proxy) & (proxy >= float(threshold))
        labels, nlab = label(detected)
        if nlab == 0:
            raise RuntimeError(
                "Automatic binning aperture found no connected component above the configured "
                "smoothed continuum threshold. Lower BINNING_APERTURE_SN_THRESHOLD only after "
                "inspecting binning_aperture.png, or use an explicit aperture mode."
            )

        cy, cx = float(center_yx[0]), float(center_yx[1])
        best_label = None
        best_distance = np.inf
        for lab in range(1, int(nlab) + 1):
            yy, xx = np.where(labels == lab)
            if yy.size == 0:
                continue
            d = float(np.min(np.hypot(yy - cy, xx - cx)))
            if d < best_distance:
                best_distance = d
                best_label = lab
        if best_label is None or best_distance > float(center_max_distance_pix):
            raise RuntimeError(
                "Automatic binning aperture did not find a detected component sufficiently close "
                f"to the adopted center (nearest={best_distance:.2f} pix, allowed="
                f"{float(center_max_distance_pix):.2f} pix)."
            )
        aperture = labels == int(best_label)
        if int(dilate_pix) > 0:
            aperture = binary_dilation(aperture, iterations=int(dilate_pix))
        aperture &= valid
        nearest_distance = float(best_distance)
    else:
        raise ValueError("BINNING_APERTURE_MODE must be 'auto_connected_sn' or 'all_good'")

    if max_radius_arcsec is not None:
        if x_arcsec is None or y_arcsec is None:
            raise ValueError("x_arcsec/y_arcsec are required when max_radius_arcsec is set")
        radius = np.hypot(np.asarray(x_arcsec, dtype=float), np.asarray(y_arcsec, dtype=float))
        aperture &= np.isfinite(radius) & (radius <= float(max_radius_arcsec))

    n_pixels = int(np.sum(aperture))
    if n_pixels < int(min_pixels):
        raise RuntimeError(
            f"Binning aperture contains only {n_pixels} pixels; minimum is {int(min_pixels)}. "
            "Inspect the aperture diagnostic before changing thresholds."
        )

    return ApertureResult(
        mask=aperture,
        significance_proxy=proxy,
        threshold=float(threshold),
        n_pixels=n_pixels,
        nearest_component_distance_pix=nearest_distance,
    )


def spatial_noise_factor(npix: int, *, mode: str, alpha: float) -> float:
    """Return optional non-additive spatial-noise inflation for one candidate bin.

    ``mode='none'`` is the production-safe baseline until a spatial covariance
    law has been calibrated for the actual KcwiKit stack.  ``mode='log10'`` is
    provided only as an explicit hook for a future empirically measured law:

        factor = 1 + alpha log10(N_spaxel).

    CRD_DAP deliberately does not choose a non-zero ``alpha`` by default.
    """
    mode = str(mode).strip().lower()
    if mode == "none":
        return 1.0
    if mode == "log10":
        if alpha < 0:
            raise ValueError("POWERBIN_SPATIAL_COVARIANCE_ALPHA cannot be negative")
        return float(1.0 + float(alpha) * np.log10(max(int(npix), 1)))
    raise ValueError("POWERBIN_SPATIAL_COVARIANCE_MODE must be 'none' or 'log10'")


def bin_sn(
    indices: np.ndarray,
    signal: np.ndarray,
    noise: np.ndarray,
    *,
    covariance_mode: str = "none",
    covariance_alpha: float = 0.0,
) -> float:
    """Compute the continuum S/N proxy of one candidate PowerBin."""
    idx = np.asarray(indices, dtype=int)
    if idx.size == 0:
        return 0.0
    s = float(np.sum(np.asarray(signal, dtype=float)[idx]))
    n2 = float(np.sum(np.asarray(noise, dtype=float)[idx] ** 2))
    if not np.isfinite(s) or not np.isfinite(n2) or n2 <= 0:
        return 0.0
    factor = spatial_noise_factor(
        idx.size,
        mode=covariance_mode,
        alpha=float(covariance_alpha),
    )
    return float(s / (np.sqrt(n2) * factor))


def run_powerbin(
    x_arcsec: np.ndarray,
    y_arcsec: np.ndarray,
    signal_map: np.ndarray,
    noise_map: np.ndarray,
    aperture_mask: np.ndarray,
    *,
    target_sn: float,
    pixel_size_arcsec: float,
    covariance_mode: str = "none",
    covariance_alpha: float = 0.0,
    regul: bool = True,
    maxiter: int = 50,
    verbose: int = 1,
) -> PowerBinResult:
    """Run mandatory Cappellari PowerBin on the BL analysis aperture."""
    try:
        from powerbin import PowerBin
        try:
            pb_version = importlib_metadata.version("powerbin")
        except importlib_metadata.PackageNotFoundError:
            pb_version = "unknown"
    except Exception as exc:  # pragma: no cover - environment-specific failure
        raise ImportError(
            "CRD_DAP Script 2 requires the PowerBin package and never silently falls back "
            "to VorBin. Install it into the active crd_dap environment with:\n\n"
            "    python -m pip install powerbin\n"
        ) from exc

    aperture = np.asarray(aperture_mask, dtype=bool)
    y_idx, x_idx = np.where(aperture)
    if y_idx.size < 2:
        raise ValueError("PowerBin requires at least two aperture pixels")

    x = np.asarray(x_arcsec, dtype=float)[aperture]
    y = np.asarray(y_arcsec, dtype=float)[aperture]
    signal = np.asarray(signal_map, dtype=float)[aperture]
    noise = np.asarray(noise_map, dtype=float)[aperture]
    if not (
        np.all(np.isfinite(x))
        and np.all(np.isfinite(y))
        and np.all(np.isfinite(signal))
        and np.all(np.isfinite(noise))
        and np.all(noise > 0)
    ):
        raise ValueError("PowerBin input aperture contains non-finite coordinates/signal/noise")

    xy = np.column_stack([x, y])

    def capacity_spec(index: np.ndarray) -> float:
        sn = bin_sn(
            index,
            signal,
            noise,
            covariance_mode=covariance_mode,
            covariance_alpha=float(covariance_alpha),
        )
        return float(max(sn, 0.0) ** 2)

    pb = PowerBin(
        xy,
        capacity_spec,
        target_capacity=float(target_sn) ** 2,
        pixelsize=float(pixel_size_arcsec),
        regul=bool(regul),
        maxiter=int(maxiter),
        verbose=int(verbose),
    )

    bin_num = np.asarray(pb.bin_num, dtype=int)
    # PowerBin bin IDs are expected to be contiguous, but remap defensively so
    # every downstream FITS table uses a deterministic 0..N-1 convention.
    unique = np.unique(bin_num)
    remap = {int(old): new for new, old in enumerate(unique.tolist())}
    mapped = np.asarray([remap[int(v)] for v in bin_num], dtype=int)
    nbin = len(unique)

    bin_map = np.full(aperture.shape, -1, dtype=np.int32)
    bin_map[y_idx, x_idx] = mapped

    # Recalculate capacities from the final membership rather than assuming a
    # particular PowerBin attribute name/version.  This also guarantees that an
    # optional CRD_DAP covariance law is represented exactly in the saved value.
    capacities = np.empty(nbin, dtype=float)
    npix = np.empty(nbin, dtype=int)
    xybin = np.empty((nbin, 2), dtype=float)
    rbin = np.empty(nbin, dtype=float)
    single = np.empty(nbin, dtype=bool)
    for bid in range(nbin):
        idx = np.flatnonzero(mapped == bid)
        npix[bid] = idx.size
        capacities[bid] = capacity_spec(idx)
        single[bid] = idx.size == 1
        # Prefer PowerBin's generators/radii when accessible; otherwise use a
        # geometric centroid and equivalent radius as a stable fallback.
        try:
            old_id = unique[bid]
            xybin[bid] = np.asarray(pb.xybin[int(old_id)], dtype=float)
            rbin[bid] = float(np.asarray(pb.rbin)[int(old_id)])
        except Exception:
            xybin[bid] = [float(np.mean(x[idx])), float(np.mean(y[idx]))]
            rbin[bid] = float(np.sqrt(idx.size / np.pi) * pixel_size_arcsec)

    non_single = ~single
    if np.any(non_single):
        target_capacity = float(target_sn) ** 2
        rms_frac = float(
            100.0
            * np.sqrt(np.mean(((capacities[non_single] - target_capacity) / target_capacity) ** 2))
        )
    else:
        rms_frac = float("nan")

    return PowerBinResult(
        bin_map=bin_map,
        aperture_y=y_idx.astype(np.int32),
        aperture_x=x_idx.astype(np.int32),
        bin_num=mapped.astype(np.int32),
        n_bins=int(nbin),
        bin_capacity=capacities,
        generator_xy=xybin,
        generator_radius=rbin,
        npix_per_bin=npix,
        single=single,
        rms_frac_percent=rms_frac,
        powerbin_version=str(pb_version),
    )


def transfer_bin_map_by_wcs(
    bl_bin_map: np.ndarray,
    bl_wcs,
    rh3_wcs,
    rh3_good_spaxel: np.ndarray,
    *,
    bl_pixel_scale_arcsec: float,
    max_distance_arcsec: float,
) -> TransferResult:
    """Assign each red/RH3 spaxel to the nearest BL bin in sky coordinates.

    The current experiment uses the same slicer and KcwiKit spatial sampling in
    both arms.  Under that condition nearest BL-pixel-center membership is an
    accurate and transparent transfer of the physical aperture.  A strict
    maximum sky-plane center separation prevents accidental use when the two
    WCS grids are materially different.
    """
    bl_bin_map = np.asarray(bl_bin_map, dtype=int)
    rh_good = np.asarray(rh3_good_spaxel, dtype=bool)
    if bl_bin_map.ndim != 2 or rh_good.ndim != 2:
        raise ValueError("BL bin map and RH3 good-spaxel mask must be 2-D")

    yy, xx = np.indices(rh_good.shape, dtype=float)
    world = rh3_wcs.pixel_to_world(xx, yy)
    xb, yb = bl_wcs.world_to_pixel(world)
    xb = np.asarray(xb, dtype=float)
    yb = np.asarray(yb, dtype=float)
    xi = np.full(xb.shape, -999999, dtype=int)
    yi = np.full(yb.shape, -999999, dtype=int)
    finite_xy = np.isfinite(xb) & np.isfinite(yb)
    xi[finite_xy] = np.rint(xb[finite_xy]).astype(int)
    yi[finite_xy] = np.rint(yb[finite_xy]).astype(int)

    inside = (
        rh_good
        & np.isfinite(xb)
        & np.isfinite(yb)
        & (xi >= 0)
        & (xi < bl_bin_map.shape[1])
        & (yi >= 0)
        & (yi < bl_bin_map.shape[0])
    )
    distance = np.full(rh_good.shape, np.nan, dtype=float)
    distance[inside] = (
        np.hypot(xb[inside] - xi[inside], yb[inside] - yi[inside])
        * float(bl_pixel_scale_arcsec)
    )

    transfer = np.full(rh_good.shape, -1, dtype=np.int32)
    candidate = inside & (distance <= float(max_distance_arcsec))
    cy, cx = np.where(candidate)
    if cy.size:
        bids = bl_bin_map[yi[cy, cx], xi[cy, cx]]
        keep = bids >= 0
        transfer[cy[keep], cx[keep]] = bids[keep].astype(np.int32)

    # Denominator counts RH3 good spaxels whose WCS lands inside a BL science
    # bin, irrespective of the distance tolerance; this makes the completeness
    # fraction specifically diagnose the physical-membership transfer.
    possible = inside.copy()
    py, px = np.where(possible)
    if py.size:
        possible_bids = bl_bin_map[yi[py, px], xi[py, px]]
        possible[py, px] &= possible_bids >= 0
    n_candidate = int(np.sum(possible))
    n_assigned = int(np.sum(transfer >= 0))
    fraction = float(n_assigned / n_candidate) if n_candidate else 0.0

    return TransferResult(
        bin_map=transfer,
        match_distance_arcsec=distance,
        n_candidate_spaxels=n_candidate,
        n_assigned_spaxels=n_assigned,
        assigned_fraction=fraction,
    )


def _surface_brightness_spatial_scale(header: Any, pixel_area_arcsec2: float) -> tuple[float, str]:
    """Return multiplicative scale for geometric aperture summation.

    KcwiKit science stacks for this project use surface-brightness units per
    square arcsecond.  When that is explicit in ``BUNIT``, summing a physical
    aperture requires multiplication by the spatial pixel area.  For other
    units we conservatively preserve the native per-spaxel normalization.
    """
    bunit = str(header.get("BUNIT", "")).lower().replace(" ", "")
    per_arcsec2_tokens = ("/arcsec2", "/arcsec^2", "arcsec-2", "arcsec**-2")
    if any(tok in bunit for tok in per_arcsec2_tokens):
        return float(pixel_area_arcsec2), "BUNIT is surface brightness per arcsec^2"
    return 1.0, "BUNIT does not explicitly indicate per-arcsec^2 surface brightness"


def coadd_bin_spectra(
    cube,
    bin_map: np.ndarray,
    *,
    n_bins: int,
    pixel_area_arcsec2: float,
    min_member_fraction: float = 0.5,
) -> CoaddedBinSpectra:
    """Geometrically sum one cube into the supplied physical bin membership."""
    bin_map = np.asarray(bin_map, dtype=int)
    if bin_map.shape != cube.good_spaxel.shape:
        raise ValueError("bin_map must match the cube spatial shape")
    if not 0.0 < float(min_member_fraction) <= 1.0:
        raise ValueError("min_member_fraction must lie in (0, 1]")

    nwave = int(cube.nwave)
    flux_out = np.full((int(n_bins), nwave), np.nan, dtype=float)
    unc_out = np.full_like(flux_out, np.nan)
    good_out = np.zeros((int(n_bins), nwave), dtype=bool)
    ncontrib_out = np.zeros((int(n_bins), nwave), dtype=np.int16)
    n_members = np.zeros(int(n_bins), dtype=np.int32)

    scale, reason = _surface_brightness_spatial_scale(cube.header, pixel_area_arcsec2)
    for bid in range(int(n_bins)):
        members = bin_map == bid
        nmem = int(np.sum(members))
        n_members[bid] = nmem
        if nmem == 0:
            continue
        flux = np.asarray(cube.flux[members, :], dtype=float)
        var = np.asarray(cube.variance[members, :], dtype=float)
        good = np.asarray(cube.good[members, :], dtype=bool)
        sample_good = good & np.isfinite(flux) & np.isfinite(var) & (var > 0)
        ncontrib = np.sum(sample_good, axis=0)
        required = max(1, int(np.ceil(float(min_member_fraction) * nmem)))
        out_good = ncontrib >= required

        summed_flux = np.sum(np.where(sample_good, flux, 0.0), axis=0) * scale
        summed_var = np.sum(np.where(sample_good, var, 0.0), axis=0) * scale**2
        flux_out[bid, out_good] = summed_flux[out_good]
        unc_out[bid, out_good] = np.sqrt(summed_var[out_good])
        good_out[bid, out_good] = True
        ncontrib_out[bid] = np.asarray(ncontrib, dtype=np.int16)

    return CoaddedBinSpectra(
        flux=flux_out,
        uncertainty=unc_out,
        good=good_out,
        contributing_spaxels=ncontrib_out,
        n_members=n_members,
        spatial_scale_factor=float(scale),
        spatial_scale_reason=reason,
    )


def achieved_sn_diagnostics_per_bin(
    spectra: CoaddedBinSpectra,
    wavelength: np.ndarray,
    *,
    rest_range: tuple[float, float],
    redshift: float,
    min_good_channels: int = 10,
    require_positive_continuum: bool = True,
) -> BinSNDiagnostics:
    """Measure robust per-bin continuum S/N and retain audit diagnostics.

    Script 2 originally used ``median(flux/uncertainty)``.  That quantity is
    useful as a signed diagnostic, but it can become numerically enormous if a
    subset of formal uncertainties is pathologically small.  It can also become
    strongly negative when a sky-subtracted continuum fluctuates below zero.

    The production-facing estimator is therefore the *ratio of robust
    locations*

        S/N = median(flux) / median(uncertainty),

    evaluated over the configured science window.  This preserves the intended
    interpretation as a typical per-spectral-pixel continuum S/N while reducing
    sensitivity to individual tiny-variance samples.

    If ``require_positive_continuum`` is True, bins with non-positive median
    continuum are assigned ``NaN`` for ``sn`` rather than a huge negative number.
    Their signed value and detailed diagnostics are retained in the returned
    object so the condition is never hidden.  No flux or variance is modified.
    """
    min_good_channels = int(min_good_channels)
    if min_good_channels < 2:
        raise ValueError("min_good_channels must be at least 2")

    obs_lo, obs_hi = observed_range_from_rest(rest_range, redshift)
    wave = np.asarray(wavelength, dtype=float)
    w = (wave >= obs_lo) & (wave <= obs_hi)
    nbin = int(spectra.flux.shape[0])

    def full(fill=np.nan, dtype=float):
        return np.full(nbin, fill, dtype=dtype)

    sn = full()
    signed_sn = full()
    legacy = full()
    med_flux = full()
    med_unc = full()
    min_unc = full()
    p05_unc = full()
    neg_frac = full()
    n_good = full(0, dtype=np.int32)
    positive = full(False, dtype=bool)

    if not np.any(w):
        return BinSNDiagnostics(
            sn=sn, signed_sn=signed_sn, legacy_median_ratio=legacy,
            median_flux=med_flux, median_uncertainty=med_unc,
            min_uncertainty=min_unc, p05_uncertainty=p05_unc,
            negative_flux_fraction=neg_frac, n_good_channels=n_good,
            positive_continuum=positive,
        )

    for bid in range(nbin):
        good = (
            spectra.good[bid]
            & w
            & np.isfinite(spectra.flux[bid])
            & np.isfinite(spectra.uncertainty[bid])
            & (spectra.uncertainty[bid] > 0)
        )
        n = int(np.sum(good))
        n_good[bid] = n
        if n < min_good_channels:
            continue

        flux = np.asarray(spectra.flux[bid, good], dtype=float)
        unc = np.asarray(spectra.uncertainty[bid, good], dtype=float)
        ratios = flux / unc

        mf = float(np.nanmedian(flux))
        mu = float(np.nanmedian(unc))
        med_flux[bid] = mf
        med_unc[bid] = mu
        min_unc[bid] = float(np.nanmin(unc))
        p05_unc[bid] = float(np.nanpercentile(unc, 5.0))
        neg_frac[bid] = float(np.mean(flux < 0.0))
        legacy[bid] = float(np.nanmedian(ratios))

        if np.isfinite(mf) and np.isfinite(mu) and mu > 0.0:
            signed_sn[bid] = mf / mu
            positive[bid] = mf > 0.0
            if (not require_positive_continuum) or positive[bid]:
                sn[bid] = signed_sn[bid]

    return BinSNDiagnostics(
        sn=sn, signed_sn=signed_sn, legacy_median_ratio=legacy,
        median_flux=med_flux, median_uncertainty=med_unc,
        min_uncertainty=min_unc, p05_uncertainty=p05_unc,
        negative_flux_fraction=neg_frac, n_good_channels=n_good,
        positive_continuum=positive,
    )


def achieved_sn_per_bin(
    spectra: CoaddedBinSpectra,
    wavelength: np.ndarray,
    *,
    rest_range: tuple[float, float],
    redshift: float,
    min_good_channels: int = 10,
    require_positive_continuum: bool = True,
) -> np.ndarray:
    """Return the robust production-facing achieved continuum S/N per bin.

    This compatibility wrapper delegates to
    :func:`achieved_sn_diagnostics_per_bin`.
    """
    return achieved_sn_diagnostics_per_bin(
        spectra,
        wavelength,
        rest_range=rest_range,
        redshift=redshift,
        min_good_channels=min_good_channels,
        require_positive_continuum=require_positive_continuum,
    ).sn


def normalized_flux_weights(
    continuum_signal_map: np.ndarray,
    bin_map: np.ndarray,
    *,
    n_bins: int,
) -> np.ndarray:
    """Return positive continuum-light weights for every spatial pixel.

    The output has the spatial shape of ``bin_map``.  Within each bin, positive
    continuum estimates are normalized to sum to unity.  If noise fluctuations
    make every member non-positive, the bin falls back to equal geometric
    weights; the caller can record this circumstance from the returned map if
    desired.
    """
    signal = np.asarray(continuum_signal_map, dtype=float)
    bin_map = np.asarray(bin_map, dtype=int)
    if signal.shape != bin_map.shape:
        raise ValueError("continuum_signal_map and bin_map must match")
    weights = np.zeros(signal.shape, dtype=float)
    for bid in range(int(n_bins)):
        members = bin_map == bid
        if not np.any(members):
            continue
        vals = np.clip(np.where(np.isfinite(signal[members]), signal[members], 0.0), 0.0, None)
        total = float(np.sum(vals))
        if total > 0:
            weights[members] = vals / total
        else:
            weights[members] = 1.0 / int(np.sum(members))
    return weights


def bin_centroids(
    bin_map: np.ndarray,
    x_arcsec: np.ndarray,
    y_arcsec: np.ndarray,
    weights: np.ndarray,
    *,
    n_bins: int,
) -> dict[str, np.ndarray]:
    """Return geometric and continuum-light-weighted centroid coordinates."""
    bin_map = np.asarray(bin_map, dtype=int)
    x = np.asarray(x_arcsec, dtype=float)
    y = np.asarray(y_arcsec, dtype=float)
    weights = np.asarray(weights, dtype=float)
    out = {
        "x_geom": np.full(int(n_bins), np.nan),
        "y_geom": np.full(int(n_bins), np.nan),
        "x_flux": np.full(int(n_bins), np.nan),
        "y_flux": np.full(int(n_bins), np.nan),
    }
    for bid in range(int(n_bins)):
        m = bin_map == bid
        if not np.any(m):
            continue
        out["x_geom"][bid] = float(np.mean(x[m]))
        out["y_geom"][bid] = float(np.mean(y[m]))
        w = weights[m]
        sw = float(np.sum(w))
        if sw > 0:
            out["x_flux"][bid] = float(np.sum(x[m] * w) / sw)
            out["y_flux"][bid] = float(np.sum(y[m] * w) / sw)
    return out
