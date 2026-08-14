"""XookSuut-style non-parametric two-disk rotation-model helpers.

The project adopts the concentric-ring philosophy and linear interpolation of
ring velocities rather than introducing a new analytic rotation-curve shape.
This module will eventually contain the complete custom objective that samples
per-bin RH3 profile-likelihood surfaces.
"""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from .geometry import deprojected_radius_azimuth


@dataclass(frozen=True)
class RingGrid:
    r_start: float
    r_final: float
    ring_space: float
    delta: float
    radii: np.ndarray


def make_ring_grid(
    fwhm_psf_arcsec: float,
    r_target_arcsec: float,
    *,
    ring_space_arcsec: float | None = None,
    delta_factor: float = 0.5,
    r_start_fwhm_factor: float = 1.0,
) -> RingGrid:
    """Construct the reproducible non-parametric radial-node grid.

    ``R_start`` is one delivered PSF FWHM from the center by default. The target
    outer radius (e.g. the adopted outer 2sigma radius) is snapped to the nearest
    allowed ring center. ``delta = delta_factor * ring_space``; the baseline
    choice of 0.5 produces non-overlapping nominal annuli whose full width equals
    the node spacing. PSF convolution still correlates neighboring radial
    measurements, so the annuli must not be called statistically independent.
    """
    if fwhm_psf_arcsec <= 0:
        raise ValueError("PSF FWHM must be positive.")
    ring_space = fwhm_psf_arcsec if ring_space_arcsec is None else float(ring_space_arcsec)
    if ring_space <= 0:
        raise ValueError("ring_space must be positive.")

    r_start = float(r_start_fwhm_factor) * float(fwhm_psf_arcsec)
    if r_target_arcsec < r_start:
        raise ValueError("Target final radius must be at or outside R_start.")

    n_final = int(np.round((r_target_arcsec - r_start) / ring_space))
    n_final = max(n_final, 0)
    r_final = r_start + n_final * ring_space
    radii = r_start + np.arange(n_final + 1, dtype=float) * ring_space
    delta = float(delta_factor) * ring_space
    return RingGrid(r_start, r_final, ring_space, delta, radii)


def interpolate_rotation_curve(radius: np.ndarray, ring_radii: np.ndarray, ring_velocities: np.ndarray) -> np.ndarray:
    """Linearly interpolate a non-parametric rotation curve between ring nodes.

    Interior to the first ring, the baseline model linearly connects
    ``V_rot(0)=0`` to the first free ring velocity. Values beyond the final ring
    are returned as NaN; the radial-extent model must explicitly decide which
    data are eligible rather than silently extrapolating the rotation curve.
    """
    r = np.asarray(radius, dtype=float)
    rr = np.asarray(ring_radii, dtype=float)
    vv = np.asarray(ring_velocities, dtype=float)
    if rr.ndim != 1 or vv.ndim != 1 or rr.size != vv.size:
        raise ValueError("ring_radii and ring_velocities must be equal-length 1-D arrays.")
    if np.any(np.diff(rr) <= 0):
        raise ValueError("ring_radii must be strictly increasing.")

    xp = np.concatenate(([0.0], rr))
    fp = np.concatenate(([0.0], vv))
    result = np.interp(r, xp, fp)
    result = np.asarray(result, dtype=float)
    result[r > rr[-1]] = np.nan
    return result


def projected_circular_velocity(
    x_arcsec: np.ndarray,
    y_arcsec: np.ndarray,
    *,
    pa_deg: float,
    inclination_deg: float,
    vsys_kms: float,
    ring_radii: np.ndarray,
    ring_velocities: np.ndarray,
) -> np.ndarray:
    r"""Evaluate the circular LOS model ``Vsys + Vrot(R) sin(i) cos(theta)``."""
    radius, azimuth = deprojected_radius_azimuth(x_arcsec, y_arcsec, pa_deg, inclination_deg)
    vrot = interpolate_rotation_curve(radius, ring_radii, ring_velocities)
    return vsys_kms + vrot * np.sin(np.deg2rad(inclination_deg)) * np.cos(azimuth)


def flux_weighted_bin_velocity(velocities: np.ndarray, weights: np.ndarray) -> tuple[float, float]:
    """Return flux-weighted mean LOS velocity and unresolved intra-bin shear.

    The shear term quantifies broadening introduced when a spatial PowerBin
    covers a velocity gradient. It is a QC quantity and may later enter the
    forward model if mocks show that centroid-only treatment biases dispersions.
    """
    v = np.asarray(velocities, dtype=float)
    w = np.asarray(weights, dtype=float)
    good = np.isfinite(v) & np.isfinite(w) & (w > 0)
    if not np.any(good):
        raise ValueError("No positive finite weights/velocities supplied.")
    v = v[good]
    w = w[good]
    mean = np.sum(w * v) / np.sum(w)
    shear2 = np.sum(w * (v - mean) ** 2) / np.sum(w)
    return float(mean), float(np.sqrt(max(shear2, 0.0)))
