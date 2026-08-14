"""Likelihood-surface utilities for RH3 and joint RH3+BL inference.

The RH3 grid stores a *profile likelihood*: at each explicit
``(V_A, V_B, f_A)`` coordinate, nuisance parameters such as dispersions and
stellar-template mixtures are optimized. Consequently the normalized weights
constructed here are called **relative-likelihood weights**, not formal Bayesian
posterior probabilities.
"""

from __future__ import annotations

from itertools import product
import numpy as np


def delta_chi2(chi2: np.ndarray, reference: float | None = None) -> np.ndarray:
    """Return ``chi2 - reference``, defaulting to the finite minimum."""
    arr = np.asarray(chi2, dtype=float)
    finite = np.isfinite(arr)
    if not np.any(finite):
        raise ValueError("chi2 array contains no finite values.")
    ref = np.nanmin(arr[finite]) if reference is None else float(reference)
    return arr - ref


def relative_likelihood_from_delta_chi2(delta: np.ndarray) -> np.ndarray:
    r"""Map :math:`\Delta\chi^2` to relative likelihood ``exp(-Delta/2)``.

    Values are shifted only through the supplied delta-chi2 definition; this
    function does not normalize them to sum to one.
    """
    delta = np.asarray(delta, dtype=float)
    out = np.zeros_like(delta, dtype=float)
    finite = np.isfinite(delta)
    out[finite] = np.exp(-0.5 * delta[finite])
    return out


def normalize_weights(weights: np.ndarray, mask: np.ndarray | None = None) -> np.ndarray:
    """Normalize non-negative weights, optionally within a boolean mask."""
    w = np.asarray(weights, dtype=float).copy()
    if mask is not None:
        mask = np.asarray(mask, dtype=bool)
        if mask.shape != w.shape:
            raise ValueError("mask and weights must have identical shapes.")
        w[~mask] = 0.0
    w[~np.isfinite(w)] = 0.0
    if np.any(w < 0):
        raise ValueError("Likelihood weights cannot be negative.")
    total = w.sum()
    if total <= 0:
        raise ValueError("Cannot normalize weights with zero total mass.")
    return w / total


def likelihood_mass_cell_counts(probability: np.ndarray, masses=(0.90, 0.95, 0.99, 0.995, 0.999)) -> dict[float, int]:
    """Number of highest-weight cells required to enclose requested masses.

    This is a diagnostic of likelihood concentration only. It is deliberately
    *not* used as a hard scientific cut in the baseline pipeline.
    """
    p = np.asarray(probability, dtype=float).ravel()
    p = p[np.isfinite(p) & (p > 0)]
    if p.size == 0:
        raise ValueError("No positive finite probability weights supplied.")
    p = p / p.sum()
    cumulative = np.cumsum(np.sort(p)[::-1])
    result = {}
    for mass in masses:
        if not 0 < mass <= 1:
            raise ValueError(f"Requested mass must be in (0, 1], got {mass}")
        result[float(mass)] = int(np.searchsorted(cumulative, mass, side="left") + 1)
    return result


def effective_number_of_states(probability: np.ndarray) -> float:
    """Return ``1/sum(p^2)`` for normalized discrete likelihood weights."""
    p = normalize_weights(probability)
    return float(1.0 / np.sum(p**2))


def neighbor_offsets(ndim: int, connectivity: str = "full") -> list[tuple[int, ...]]:
    """Return offsets for face-only or full immediate-neighbor connectivity.

    In three dimensions, ``connectivity='full'`` yields the 26 cells touching
    the current cell by a face, edge, or corner. ``'face'`` yields 6 neighbors.
    """
    if connectivity not in {"full", "face"}:
        raise ValueError("connectivity must be 'full' or 'face'.")
    offsets = []
    for offset in product((-1, 0, 1), repeat=ndim):
        if all(v == 0 for v in offset):
            continue
        if connectivity == "face" and sum(v != 0 for v in offset) != 1:
            continue
        offsets.append(offset)
    return offsets


def _valid_neighbor(index: tuple[int, ...], offset: tuple[int, ...], shape: tuple[int, ...]) -> tuple[int, ...] | None:
    candidate = tuple(i + d for i, d in zip(index, offset))
    if all(0 <= c < n for c, n in zip(candidate, shape)):
        return candidate
    return None


def descend_to_local_minimum(
    chi2: np.ndarray,
    start_index: tuple[int, ...],
    *,
    connectivity: str = "full",
) -> tuple[int, ...]:
    """Follow discrete steepest descent to the local minimum reached by a cell.

    This helper is intended for identifying the likelihood basin associated with
    the global-model anchor. It is a topological operation: cells in the same
    basin need not have similar chi2 values; they share the same downhill
    destination.

    Ties are resolved deterministically by lexicographic index ordering. The
    exact neighbor convention will be validated on injection/recovery cubes
    before publication use.
    """
    arr = np.asarray(chi2, dtype=float)
    if len(start_index) != arr.ndim:
        raise ValueError("start_index dimensionality does not match chi2 array.")
    if not np.isfinite(arr[start_index]):
        raise ValueError("start_index points to a non-finite chi2 cell.")

    offsets = neighbor_offsets(arr.ndim, connectivity)
    current = tuple(start_index)
    visited = set()

    while True:
        if current in visited:
            # A strict decrease should prevent cycles. Reaching this branch
            # indicates numerical/pathological input and should not be hidden.
            raise RuntimeError("Cycle encountered during chi2 steepest descent.")
        visited.add(current)

        candidates = [(arr[current], current)]
        for offset in offsets:
            neighbor = _valid_neighbor(current, offset, arr.shape)
            if neighbor is not None and np.isfinite(arr[neighbor]):
                candidates.append((arr[neighbor], neighbor))

        best_value = min(value for value, _ in candidates)
        best_indices = sorted(index for value, index in candidates if value == best_value)
        best = best_indices[0]

        if best_value < arr[current]:
            current = best
            continue
        return current


def basin_mask_for_anchor(
    chi2: np.ndarray,
    anchor_index: tuple[int, ...],
    *,
    connectivity: str = "full",
) -> tuple[np.ndarray, tuple[int, ...]]:
    """Return all finite cells descending to the same minimum as an anchor.

    This brute-force implementation is intentionally clear and appropriate for
    development-scale grids such as 17x17x9. If future grids become much larger,
    the mapping can be memoized/optimized without changing the scientific
    definition.
    """
    arr = np.asarray(chi2, dtype=float)
    anchor_minimum = descend_to_local_minimum(arr, anchor_index, connectivity=connectivity)
    mask = np.zeros(arr.shape, dtype=bool)

    for index in np.ndindex(arr.shape):
        if not np.isfinite(arr[index]):
            continue
        destination = descend_to_local_minimum(arr, index, connectivity=connectivity)
        if destination == anchor_minimum:
            mask[index] = True

    return mask, anchor_minimum
