"""Cube preparation, continuum images, center estimates, and arm registration.

The functions here operate on the standardized ``(y, x, wavelength)`` arrays
returned by :mod:`crd_utils.io`.  They intentionally keep *hard data-quality*
operations separate from later science-fit masks.  In particular, low S/N is
not a reason to reject a spatial element in Script 1; Script 2 uses PowerBin to
combine low-S/N spaxels.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from scipy.ndimage import gaussian_filter, map_coordinates
from scipy.signal import fftconvolve


@dataclass(frozen=True)
class CenterEstimate:
    """Continuum-based center diagnostics for one cube."""

    peak_yx: tuple[float, float]
    centroid_yx: tuple[float, float]
    peak_sky: SkyCoord
    centroid_sky: SkyCoord


@dataclass(frozen=True)
class RegistrationResult:
    """BL/RH3 registration diagnostic on the BL spatial grid.

    ``cross_correlation_valid`` is deliberately explicit.  A numerical
    cross-correlation shift is only reported when *both* registration images
    contain enough spatial contrast to make morphology registration meaningful.
    Otherwise the returned shifts are NaN and Script 1 falls back to the
    independently measured sky-coordinate center comparison as the relevant QC.
    """

    moving_on_reference: np.ndarray
    difference: np.ndarray
    overlap: np.ndarray
    residual_shift_yx_pix: tuple[float, float]
    residual_shift_arcsec: tuple[float, float]
    residual_shift_radius_arcsec: float
    cross_correlation_valid: bool
    reference_contrast_snr: float
    moving_contrast_snr: float
    status_reason: str


def collapsed_continuum(
    flux: np.ndarray,
    good: np.ndarray | None = None,
    *,
    statistic: str = "median",
) -> np.ndarray:
    """Collapse a spectral cube into a robust 2-D continuum image.

    Parameters
    ----------
    flux
        Array with wavelength on the final axis.
    good
        Optional boolean mask with the same shape as ``flux``. True means the
        sample is usable.
    statistic
        ``"median"`` or ``"mean"``.  The median is the production default
        because it is less sensitive to isolated emission/sky residuals.
    """
    f = np.asarray(flux, dtype=float)
    if f.ndim != 3:
        raise ValueError("flux must be a 3-D cube with wavelength on the final axis")
    if good is not None:
        good = np.asarray(good, dtype=bool)
        if good.shape != f.shape:
            raise ValueError("good mask must match flux shape")
        f = np.where(good, f, np.nan)

    with np.errstate(all="ignore"):
        if statistic == "median":
            return np.nanmedian(f, axis=-1)
        if statistic == "mean":
            return np.nanmean(f, axis=-1)
    raise ValueError("statistic must be 'median' or 'mean'")


def smooth_image(image: np.ndarray, sigma_pix: float = 1.0) -> np.ndarray:
    """Gaussian-smooth an image without allowing NaNs to bias nearby pixels."""
    img = np.asarray(image, dtype=float)
    if sigma_pix <= 0:
        return img.copy()
    valid = np.isfinite(img)
    values = np.where(valid, img, 0.0)
    weights = gaussian_filter(valid.astype(float), sigma=sigma_pix, mode="nearest")
    smoothed = gaussian_filter(values, sigma=sigma_pix, mode="nearest")
    out = np.full_like(img, np.nan, dtype=float)
    good = weights > 0
    out[good] = smoothed[good] / weights[good]
    return out


def brightest_spaxel(image: np.ndarray) -> tuple[int, int]:
    """Return the brightest finite *positive* image pixel as ``(y, x)``."""
    image = np.asarray(image, dtype=float)
    valid = np.isfinite(image) & (image > 0)
    if not np.any(valid):
        raise ValueError("Image contains no finite positive pixels")
    masked = np.where(valid, image, -np.inf)
    return tuple(int(i) for i in np.unravel_index(np.argmax(masked), image.shape))


def flux_weighted_centroid(
    image: np.ndarray,
    *,
    min_percentile: float = 60.0,
) -> tuple[float, float]:
    """Return a robust positive-flux centroid ``(y, x)``.

    The percentile floor prevents large areas of low-level background noise from
    dominating the centroid of an extended galaxy.  This centroid is a *QC
    comparison*, not the final kinematic center used by Script 4.
    """
    img = np.asarray(image, dtype=float)
    finite_positive = np.isfinite(img) & (img > 0)
    if np.sum(finite_positive) < 3:
        raise ValueError("Too few positive pixels for a continuum centroid")

    threshold = np.nanpercentile(img[finite_positive], float(min_percentile))
    use = finite_positive & (img >= threshold)
    weights = np.where(use, img, 0.0)
    total = float(np.sum(weights))
    if not np.isfinite(total) or total <= 0:
        raise ValueError("Centroid weights are not positive/finite")

    yy, xx = np.indices(img.shape, dtype=float)
    y = float(np.sum(weights * yy) / total)
    x = float(np.sum(weights * xx) / total)
    return y, x


def estimate_continuum_center(
    image: np.ndarray,
    celestial_wcs: WCS,
    *,
    smooth_sigma_pix: float = 1.0,
    centroid_min_percentile: float = 60.0,
) -> CenterEstimate:
    """Measure continuum peak and centroid in both pixel and sky coordinates."""
    smoothed = smooth_image(image, smooth_sigma_pix)
    peak_y, peak_x = brightest_spaxel(smoothed)
    cen_y, cen_x = flux_weighted_centroid(
        smoothed,
        min_percentile=centroid_min_percentile,
    )
    peak_sky = celestial_wcs.pixel_to_world(float(peak_x), float(peak_y))
    centroid_sky = celestial_wcs.pixel_to_world(float(cen_x), float(cen_y))
    return CenterEstimate(
        peak_yx=(float(peak_y), float(peak_x)),
        centroid_yx=(float(cen_y), float(cen_x)),
        peak_sky=peak_sky,
        centroid_sky=centroid_sky,
    )


def spatial_offset_grids_arcsec(
    celestial_wcs: WCS,
    shape_yx: tuple[int, int],
    center_sky: SkyCoord,
) -> tuple[np.ndarray, np.ndarray]:
    """Return tangent-plane ``(x, y)`` offset grids in arcsec from ``center_sky``."""
    ny, nx = shape_yx
    yy, xx = np.indices((ny, nx), dtype=float)
    sky = celestial_wcs.pixel_to_world(xx, yy)
    dlon, dlat = center_sky.spherical_offsets_to(sky)
    return dlon.arcsec.astype(float), dlat.arcsec.astype(float)


def reproject_image_via_wcs(
    moving_image: np.ndarray,
    moving_wcs: WCS,
    reference_wcs: WCS,
    reference_shape: tuple[int, int],
    *,
    order: int = 1,
) -> tuple[np.ndarray, np.ndarray]:
    """Sample ``moving_image`` onto the reference celestial pixel grid.

    This lightweight diagnostic reprojection avoids a dependency on the
    ``reproject`` package.  It is used only to *measure/check registration*; the
    science cubes themselves are not resampled in Script 1, which avoids adding
    another interpolation/covariance operation before PowerBin.
    """
    moving = np.asarray(moving_image, dtype=float)
    ny, nx = reference_shape
    yy, xx = np.indices((ny, nx), dtype=float)
    sky = reference_wcs.pixel_to_world(xx, yy)
    mx, my = moving_wcs.world_to_pixel(sky)

    inside = (
        np.isfinite(mx)
        & np.isfinite(my)
        & (mx >= 0)
        & (my >= 0)
        & (mx <= moving.shape[1] - 1)
        & (my <= moving.shape[0] - 1)
    )

    coords = np.vstack([my.ravel(), mx.ravel()])
    filled = np.where(np.isfinite(moving), moving, 0.0)
    sampled = map_coordinates(filled, coords, order=order, mode="constant", cval=np.nan)
    sampled = sampled.reshape(reference_shape)
    sampled[~inside] = np.nan
    return sampled, inside


def _normalized_for_correlation(image: np.ndarray, mask: np.ndarray) -> np.ndarray:
    arr = np.asarray(image, dtype=float)
    use = mask & np.isfinite(arr)
    if np.sum(use) < 9:
        raise ValueError("Insufficient overlap for registration cross-correlation")
    values = arr[use]
    med = np.median(values)
    mad = np.median(np.abs(values - med))
    scale = 1.4826 * mad
    if not np.isfinite(scale) or scale <= 0:
        scale = np.std(values)
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError("Registration image has zero/invalid robust scale")
    out = np.zeros_like(arr, dtype=float)
    out[use] = (arr[use] - med) / scale
    return out




def registration_contrast_snr(
    image: np.ndarray,
    mask: np.ndarray,
    *,
    smooth_sigma_pix: float = 1.0,
) -> float:
    """Return a robust morphology-contrast statistic for registration QC.

    The statistic is intentionally simple and diagnostic rather than a formal
    photometric S/N.  After light NaN-aware smoothing, it compares the upper
    spatial tail of the image to the robust pixel-to-pixel background scale:

    ``contrast_snr = (P99 - median) / (1.4826 * MAD)``.

    A nearly flat collapsed-continuum image can otherwise produce a formally
    precise but scientifically meaningless cross-correlation shift after robust
    normalization.  Script 1 therefore requires this contrast statistic to
    exceed a configurable threshold in *both* arms before trusting the shift.
    """
    arr = smooth_image(np.asarray(image, dtype=float), sigma_pix=smooth_sigma_pix)
    use = np.asarray(mask, dtype=bool) & np.isfinite(arr)
    if np.sum(use) < 25:
        return float("nan")

    values = arr[use]
    med = float(np.nanmedian(values))
    mad = float(np.nanmedian(np.abs(values - med)))
    scale = 1.4826 * mad
    if not np.isfinite(scale) or scale <= 0:
        scale = float(np.nanstd(values))
    if not np.isfinite(scale) or scale <= 0:
        return 0.0

    high = float(np.nanpercentile(values, 99.0))
    return float(max(0.0, (high - med) / scale))


def _parabolic_offset(v_minus: float, v0: float, v_plus: float) -> float:
    """Subpixel peak offset for three equally spaced samples, clipped to ±1."""
    denom = v_minus - 2.0 * v0 + v_plus
    if not np.isfinite(denom) or abs(denom) < 1e-12:
        return 0.0
    offset = 0.5 * (v_minus - v_plus) / denom
    return float(np.clip(offset, -1.0, 1.0))


def residual_registration_shift(
    reference_image: np.ndarray,
    moving_on_reference: np.ndarray,
    overlap: np.ndarray,
    *,
    max_shift_pix: tuple[float, float] | None = None,
) -> tuple[float, float]:
    """Estimate residual ``(dy, dx)`` shift after WCS reprojection.

    The images are robustly normalized and cross-correlated.  A simple
    three-point parabola provides subpixel refinement near the correlation peak.
    The returned shift is the additional shift that would align the moving image
    to the reference image on the reference grid.

    ``max_shift_pix=(max_dy, max_dx)`` limits the search to a physically
    plausible residual around the WCS-predicted alignment.  Script 1 is testing
    *residual* registration, not trying to recover an arbitrarily large blind
    translation.  This prevents footprint edges or weak continuum structure
    from winning the correlation at tens of pixels from zero shift.
    """
    ref = np.asarray(reference_image, dtype=float)
    mov = np.asarray(moving_on_reference, dtype=float)
    use = np.asarray(overlap, dtype=bool) & np.isfinite(ref) & np.isfinite(mov)

    a = _normalized_for_correlation(ref, use)
    b = _normalized_for_correlation(mov, use)
    weight = use.astype(float)

    corr = fftconvolve(a, b[::-1, ::-1], mode="same")
    norm = fftconvolve(weight, weight[::-1, ::-1], mode="same")
    valid = norm > max(9.0, 0.1 * np.nanmax(norm))

    # Restrict the correlation search to a local window around zero residual
    # shift when requested.  Large offsets should be fixed in the WCS/stacking
    # stage rather than silently adopted by this diagnostic.
    cy, cx = np.array(corr.shape) // 2
    if max_shift_pix is not None:
        max_dy, max_dx = (float(max_shift_pix[0]), float(max_shift_pix[1]))
        if max_dy <= 0 or max_dx <= 0:
            raise ValueError("max_shift_pix values must be positive")
        yy, xx = np.indices(corr.shape)
        local = (
            (np.abs(yy - cy) <= max_dy)
            & (np.abs(xx - cx) <= max_dx)
        )
        valid &= local

    corr_norm = np.full_like(corr, -np.inf, dtype=float)
    corr_norm[valid] = corr[valid] / norm[valid]
    if not np.any(np.isfinite(corr_norm)):
        raise ValueError("No valid correlation samples inside registration search window")

    iy, ix = np.unravel_index(np.nanargmax(corr_norm), corr_norm.shape)
    dy = float(iy - cy)
    dx = float(ix - cx)

    if 0 < iy < corr_norm.shape[0] - 1:
        dy += _parabolic_offset(corr_norm[iy - 1, ix], corr_norm[iy, ix], corr_norm[iy + 1, ix])
    if 0 < ix < corr_norm.shape[1] - 1:
        dx += _parabolic_offset(corr_norm[iy, ix - 1], corr_norm[iy, ix], corr_norm[iy, ix + 1])
    return dy, dx


def register_cube_pair(
    reference_image: np.ndarray,
    reference_wcs: WCS,
    moving_image: np.ndarray,
    moving_wcs: WCS,
    *,
    reference_pixel_scale_arcsec: tuple[float, float],
    min_contrast_snr: float = 5.0,
    contrast_smooth_sigma_pix: float = 1.0,
    max_residual_shift_arcsec: float = 2.0,
) -> RegistrationResult:
    """Diagnose RH3-to-BL registration without resampling either science cube.

    The moving image is sampled onto the reference WCS only for this QC plot.
    Before cross-correlation, both images must pass a morphology-contrast test.
    This prevents a nearly featureless continuum image from producing a
    numerically well-defined but physically meaningless residual shift.
    """
    projected, overlap = reproject_image_via_wcs(
        moving_image,
        moving_wcs,
        reference_wcs,
        reference_image.shape,
    )

    ref = np.asarray(reference_image, dtype=float)
    use = overlap & np.isfinite(ref) & np.isfinite(projected)
    ref_contrast = registration_contrast_snr(
        ref, use, smooth_sigma_pix=contrast_smooth_sigma_pix
    )
    mov_contrast = registration_contrast_snr(
        projected, use, smooth_sigma_pix=contrast_smooth_sigma_pix
    )

    contrast_ok = (
        np.isfinite(ref_contrast)
        and np.isfinite(mov_contrast)
        and ref_contrast >= float(min_contrast_snr)
        and mov_contrast >= float(min_contrast_snr)
    )

    if contrast_ok:
        scale_x, scale_y = reference_pixel_scale_arcsec
        max_arc = float(max_residual_shift_arcsec)
        if not np.isfinite(max_arc) or max_arc <= 0:
            raise ValueError("max_residual_shift_arcsec must be positive")
        max_dy_pix = max_arc / float(scale_y)
        max_dx_pix = max_arc / float(scale_x)

        dy_candidate, dx_candidate = residual_registration_shift(
            ref,
            projected,
            use,
            max_shift_pix=(max_dy_pix, max_dx_pix),
        )

        # A best-fit peak at the edge of the allowed local search window means
        # the morphology correlation is trying to run away from the WCS solution.
        # Treat that as inconclusive rather than reporting the boundary as a
        # measured astrometric shift.
        boundary_margin_pix = 0.75
        hit_boundary = (
            abs(dy_candidate) >= max(0.0, max_dy_pix - boundary_margin_pix)
            or abs(dx_candidate) >= max(0.0, max_dx_pix - boundary_margin_pix)
        )

        # Normalize each image by a robust scale before differencing so the
        # display emphasizes morphology rather than passband flux normalization.
        ref_n = _normalized_for_correlation(ref, use)
        mov_n = _normalized_for_correlation(projected, use)
        diff = np.full_like(ref_n, np.nan)
        diff[use] = ref_n[use] - mov_n[use]

        if hit_boundary:
            dy = dx = dy_arc = dx_arc = float("nan")
            status = (
                "cross-correlation inconclusive: best peak reached the local "
                f"search boundary (max residual {max_arc:.3g} arcsec); "
                "inspect WCS/centers rather than adopting a large blind shift"
            )
            valid = False
        else:
            dy = float(dy_candidate)
            dx = float(dx_candidate)
            dx_arc = dx * scale_x
            dy_arc = dy * scale_y
            status = "cross-correlation valid within local residual-search window"
            valid = True
    else:
        dy = dx = dy_arc = dx_arc = float("nan")
        diff = np.full_like(ref, np.nan, dtype=float)
        status = (
            "cross-correlation inconclusive: insufficient morphology contrast "
            f"(reference={ref_contrast:.3g}, moving={mov_contrast:.3g}, "
            f"required>={float(min_contrast_snr):.3g})"
        )
        valid = False

    return RegistrationResult(
        moving_on_reference=projected,
        difference=diff,
        overlap=overlap,
        residual_shift_yx_pix=(dy, dx),
        residual_shift_arcsec=(dy_arc, dx_arc),
        residual_shift_radius_arcsec=float(np.hypot(dx_arc, dy_arc)),
        cross_correlation_valid=valid,
        reference_contrast_snr=float(ref_contrast),
        moving_contrast_snr=float(mov_contrast),
        status_reason=status,
    )
