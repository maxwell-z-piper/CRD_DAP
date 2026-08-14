"""Scientifically important validation checks used throughout the pipeline."""

from __future__ import annotations

import numpy as np


def check_wavelength_convention(science_medium: str, template_medium: str) -> None:
    """Hard-fail if science and template wavelength media are inconsistent."""
    allowed = {"air", "vacuum"}
    s = science_medium.lower()
    t = template_medium.lower()
    if s not in allowed or t not in allowed:
        raise ValueError("Wavelength medium must be 'air' or 'vacuum'.")
    if s != t:
        raise ValueError(
            f"Science wavelengths are {s} while templates are {t}. Convert one dataset before pPXF fitting."
        )


def grid_edge_distance_cells(index: tuple[int, ...], shape: tuple[int, ...]) -> int:
    """Minimum number of cells between a selected index and any grid edge."""
    if len(index) != len(shape):
        raise ValueError("index and shape dimensionality differ.")
    distances = []
    for i, n in zip(index, shape):
        if not 0 <= i < n:
            raise IndexError("Grid index lies outside shape.")
        distances.extend([i, n - 1 - i])
    return int(min(distances))


def is_grid_edge_warning(index: tuple[int, ...], shape: tuple[int, ...], warning_cells: int = 2) -> bool:
    """True when selected state is within ``warning_cells`` of any edge."""
    return grid_edge_distance_cells(index, shape) <= warning_cells


def compare_center_estimates(*args, **kwargs):
    """Compare continuum peak/centroid, PyMorph, and final kinematic centers."""
    raise NotImplementedError("Implemented with Scripts 1 and 4.")
