"""Photometric/kinematic geometry helpers and the permanent Disk A/B convention."""

from __future__ import annotations

import numpy as np


def inclination_from_axis_ratio(q: np.ndarray | float, q0: np.ndarray | float) -> np.ndarray:
    r"""Convert observed axis ratio to inclination using a finite disk thickness.

    The adopted relation is

    .. math::

        \cos^2 i = \frac{q^2-q_0^2}{1-q_0^2},

    where ``q=b/a`` is the observed disk axis ratio and ``q0`` is the assumed
    intrinsic thickness. Inputs for which ``q < q0`` are clipped to edge-on
    (cos^2 i = 0) rather than producing an invalid square root.
    """
    q = np.asarray(q, dtype=float)
    q0 = np.asarray(q0, dtype=float)
    if np.any((q <= 0) | (q > 1)):
        raise ValueError("Observed axis ratio q must lie in (0, 1].")
    if np.any((q0 < 0) | (q0 >= 1)):
        raise ValueError("Intrinsic thickness q0 must lie in [0, 1).")

    cos2 = (q**2 - q0**2) / (1.0 - q0**2)
    cos2 = np.clip(cos2, 0.0, 1.0)
    return np.degrees(np.arccos(np.sqrt(cos2)))


def sample_inclination_prior(
    q: float,
    q_err: float,
    q0: float,
    q0_err: float,
    *,
    n_draws: int = 10000,
    seed: int = 12345,
) -> np.ndarray:
    """Propagate PyMorph axis-ratio and intrinsic-thickness uncertainty to i.

    Gaussian draws are truncated to physical ranges. The returned distribution
    can be summarized by median/percentiles and used to construct an informative
    inclination prior for the global disk model.
    """
    rng = np.random.default_rng(seed)
    q_draw = rng.normal(q, q_err, n_draws)
    q0_draw = rng.normal(q0, q0_err, n_draws)
    q_draw = np.clip(q_draw, 1e-4, 1.0)
    q0_draw = np.clip(q0_draw, 0.0, 0.999)
    return inclination_from_axis_ratio(q_draw, q0_draw)


def rotate_to_major_minor(
    x_arcsec: np.ndarray,
    y_arcsec: np.ndarray,
    pa_deg: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Rotate sky-plane offsets into major/minor-axis coordinates.

    This low-level helper assumes a documented internal angle convention. The
    final Script-1/4 implementation must explicitly translate FITS/WCS and
    ``fit_kinematic_pa`` conventions into this internal coordinate definition.
    """
    theta = np.deg2rad(pa_deg)
    x = np.asarray(x_arcsec, dtype=float)
    y = np.asarray(y_arcsec, dtype=float)
    major = x * np.cos(theta) + y * np.sin(theta)
    minor = -x * np.sin(theta) + y * np.cos(theta)
    return major, minor


def deprojected_radius_azimuth(
    x_arcsec: np.ndarray,
    y_arcsec: np.ndarray,
    pa_deg: float,
    inclination_deg: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Return disk-plane radius and azimuth from sky-plane offsets."""
    major, minor = rotate_to_major_minor(x_arcsec, y_arcsec, pa_deg)
    cos_i = np.cos(np.deg2rad(inclination_deg))
    if cos_i <= 0:
        raise ValueError("Inclination must be < 90 degrees for deprojection.")
    minor_disk = minor / cos_i
    radius = np.hypot(major, minor_disk)
    azimuth = np.arctan2(minor_disk, major)
    return radius, azimuth


def signed_pa_coordinate(x_arcsec: np.ndarray, y_arcsec: np.ndarray, pa_deg: float) -> np.ndarray:
    """Signed coordinate along the adopted PA axis used for deterministic labels."""
    major, _ = rotate_to_major_minor(x_arcsec, y_arcsec, pa_deg)
    return major
