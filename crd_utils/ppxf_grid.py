"""Explicit RH3 velocity/fraction-grid construction and profile-likelihood cubes."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np

from .ppxf_utils import PPXFStateResult, fit_fixed_two_component_state


FIT_STATUS_SUCCESS = np.int8(0)
FIT_STATUS_PPXF_FAILURE = np.int8(1)
FIT_STATUS_FIXED_VELOCITY_MISMATCH = np.int8(2)


@dataclass
class RH3LikelihoodCube:
    chi2_total: np.ndarray
    reduced_chi2: np.ndarray
    sigma_a: np.ndarray
    sigma_b: np.ndarray
    fit_status: np.ndarray
    sigma_boundary: np.ndarray
    va_grid: np.ndarray
    vb_grid: np.ndarray
    fa_grid: np.ndarray
    best_index: tuple[int, int, int] | None
    n_failures: int
    n_sigma_boundary: int


def uniform_grid(minimum: float, maximum: float, n: int) -> np.ndarray:
    """Inclusive uniformly sampled one-dimensional grid."""
    if n < 2:
        raise ValueError("n must be >= 2.")
    if maximum <= minimum:
        raise ValueError("maximum must exceed minimum.")
    return np.linspace(float(minimum), float(maximum), int(n))


def fraction_grid(minimum: float, maximum: float, step: float) -> np.ndarray:
    """Inclusive fraction grid robust to floating-point endpoint rounding."""
    if not (0 <= minimum < maximum <= 1):
        raise ValueError("Fraction bounds must satisfy 0 <= min < max <= 1.")
    if step <= 0:
        raise ValueError("step must be positive.")
    n = int(np.floor((maximum - minimum) / step + 0.5))
    values = minimum + np.arange(n + 1) * step
    values = values[values <= maximum + step * 1e-8]
    return np.round(values, 12)


def build_rh3_likelihood_cube(
    *,
    templates_two_component: np.ndarray,
    component: np.ndarray,
    galaxy: np.ndarray,
    noise: np.ndarray,
    velscale: float,
    lam: np.ndarray,
    lam_temp: np.ndarray,
    goodpixels: np.ndarray,
    va_grid: np.ndarray,
    vb_grid: np.ndarray,
    fa_grid: np.ndarray,
    sigma_start_a: float,
    sigma_start_b: float,
    sigma_bounds: tuple[float, float],
    degree: int,
    mdegree: int,
    regul: float = 0.0,
    sigma_boundary_tolerance_kms: float = 2.0,
    progress_callback: Callable[[int, int], None] | None = None,
) -> RH3LikelihoodCube:
    """Evaluate the complete exact ``(V_A,V_B,f_A)`` grid for one bin.

    The returned object is a profile-likelihood cube: pPXF minimizes over the
    two dispersions, template weights, and additive-polynomial coefficients at
    every exact velocity/fraction coordinate.
    """
    va_grid = np.asarray(va_grid, dtype=float)
    vb_grid = np.asarray(vb_grid, dtype=float)
    fa_grid = np.asarray(fa_grid, dtype=float)
    shape = (va_grid.size, vb_grid.size, fa_grid.size)

    chi2 = np.full(shape, np.inf, dtype=np.float64)
    reduced = np.full(shape, np.nan, dtype=np.float32)
    sig_a = np.full(shape, np.nan, dtype=np.float32)
    sig_b = np.full(shape, np.nan, dtype=np.float32)
    status = np.full(shape, FIT_STATUS_PPXF_FAILURE, dtype=np.int8)
    boundary = np.zeros(shape, dtype=np.uint8)

    smin, smax = map(float, sigma_bounds)
    tol = max(0.0, float(sigma_boundary_tolerance_kms))
    total_states = int(np.prod(shape))
    completed = 0

    for ia, va in enumerate(va_grid):
        for ib, vb in enumerate(vb_grid):
            for jf, fa in enumerate(fa_grid):
                result = fit_fixed_two_component_state(
                    templates_two_component=templates_two_component,
                    component=component,
                    galaxy=galaxy,
                    noise=noise,
                    velscale=velscale,
                    lam=lam,
                    lam_temp=lam_temp,
                    goodpixels=goodpixels,
                    velocity_a=float(va),
                    velocity_b=float(vb),
                    fraction_a=float(fa),
                    start_sigma_a=float(sigma_start_a),
                    start_sigma_b=float(sigma_start_b),
                    sigma_bounds=(smin, smax),
                    degree=int(degree),
                    mdegree=int(mdegree),
                    regul=float(regul),
                    keep_full=False,
                )
                idx = (ia, ib, jf)
                if result.success:
                    # The fixed-velocity coordinates are scientifically exact.
                    # Treat a surprising pPXF return value as a failed state
                    # instead of silently turning the grid into a half-cell fit.
                    if not (
                        np.isclose(result.velocity[0], va, atol=1.0e-6, rtol=0.0)
                        and np.isclose(result.velocity[1], vb, atol=1.0e-6, rtol=0.0)
                    ):
                        status[idx] = FIT_STATUS_FIXED_VELOCITY_MISMATCH
                    else:
                        chi2[idx] = result.chi2_total
                        reduced[idx] = result.reduced_chi2
                        sig_a[idx] = result.sigma[0]
                        sig_b[idx] = result.sigma[1]
                        status[idx] = FIT_STATUS_SUCCESS
                        near = (
                            (result.sigma[0] <= smin + tol)
                            or (result.sigma[0] >= smax - tol)
                            or (result.sigma[1] <= smin + tol)
                            or (result.sigma[1] >= smax - tol)
                        )
                        boundary[idx] = np.uint8(1 if near else 0)
                completed += 1
                if progress_callback is not None:
                    progress_callback(completed, total_states)

    finite = np.isfinite(chi2)
    best_index = None
    if np.any(finite):
        flat = int(np.nanargmin(np.where(finite, chi2, np.nan)))
        best_index = tuple(int(x) for x in np.unravel_index(flat, shape))

    return RH3LikelihoodCube(
        chi2_total=chi2,
        reduced_chi2=reduced,
        sigma_a=sig_a,
        sigma_b=sig_b,
        fit_status=status,
        sigma_boundary=boundary,
        va_grid=va_grid,
        vb_grid=vb_grid,
        fa_grid=fa_grid,
        best_index=best_index,
        n_failures=int(np.sum(status != FIT_STATUS_SUCCESS)),
        n_sigma_boundary=int(np.sum(boundary > 0)),
    )
