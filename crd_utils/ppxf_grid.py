"""Explicit RH3 velocity/fraction-grid construction and profile-likelihood cubes."""

from __future__ import annotations

import numpy as np


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


def build_rh3_likelihood_cube(*args, **kwargs):
    """Evaluate pPXF over explicit (V_A, V_B, f_A) coordinates for one bin."""
    raise NotImplementedError("Implemented with Script 3.")
