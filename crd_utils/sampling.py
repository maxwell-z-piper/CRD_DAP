"""Likelihood-weighted RH3 state sampling and convergence diagnostics.

A crucial bookkeeping rule is implemented conceptually here: if RH3 grid cells
are *drawn according to their RH3 relative-likelihood weights*, the RH3
likelihood is already encoded in the sampling frequency. The downstream BL
weight must not multiply every draw by the RH3 likelihood again, which would
count RH3 twice.
"""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from .likelihood import normalize_weights


@dataclass(frozen=True)
class ConvergenceResult:
    converged: bool
    max_change: float
    tolerance: float
    per_quantile_changes: dict[str, float]


def sample_discrete_states(
    probability: np.ndarray,
    n_draws: int,
    *,
    rng: np.random.Generator | None = None,
) -> np.ndarray:
    """Draw flattened grid-cell indices according to normalized likelihood."""
    if n_draws < 1:
        raise ValueError("n_draws must be >= 1.")
    p = normalize_weights(probability).ravel()
    rng = np.random.default_rng() if rng is None else rng
    return rng.choice(p.size, size=n_draws, replace=True, p=p)


def effective_sample_size(weights: np.ndarray) -> float:
    r"""Return importance-weight effective sample size.

    .. math::
        N_\mathrm{eff} = \frac{(\sum w)^2}{\sum w^2}.
    """
    w = np.asarray(weights, dtype=float)
    good = np.isfinite(w) & (w >= 0)
    w = w[good]
    if w.size == 0 or np.sum(w) <= 0:
        return 0.0
    return float((np.sum(w) ** 2) / np.sum(w**2))


def posterior_quantiles(values: np.ndarray, weights: np.ndarray, quantiles=(0.16, 0.50, 0.84)) -> np.ndarray:
    """Weighted quantiles for a one-dimensional discrete solution family."""
    x = np.asarray(values, dtype=float)
    w = np.asarray(weights, dtype=float)
    good = np.isfinite(x) & np.isfinite(w) & (w >= 0)
    x = x[good]
    w = w[good]
    if x.size == 0 or np.sum(w) <= 0:
        raise ValueError("No finite positively weighted values supplied.")

    order = np.argsort(x)
    x = x[order]
    w = w[order]
    cumulative = np.cumsum(w) / np.sum(w)
    return np.interp(np.asarray(quantiles, dtype=float), cumulative, x)


def check_quantile_convergence(
    previous: np.ndarray,
    current: np.ndarray,
    *,
    fraction_of_ci: float = 0.05,
    absolute_floor: float = 0.005,
) -> ConvergenceResult:
    """Test stability of the 16th/50th/84th percentile summary.

    The tolerance is the larger of an absolute numerical floor and a fraction
    of the current 68% interval width. This is more stable and interpretable
    than demanding that printed values remain identical after decimal rounding.
    """
    previous = np.asarray(previous, dtype=float)
    current = np.asarray(current, dtype=float)
    if previous.shape != (3,) or current.shape != (3,):
        raise ValueError("previous and current must contain [P16, P50, P84].")
    width = max(float(current[2] - current[0]), 0.0)
    tolerance = max(float(absolute_floor), float(fraction_of_ci) * width)
    changes = np.abs(current - previous)
    labels = ("P16", "P50", "P84")
    per = {label: float(change) for label, change in zip(labels, changes)}
    max_change = float(np.max(changes))
    return ConvergenceResult(max_change < tolerance, max_change, tolerance, per)
