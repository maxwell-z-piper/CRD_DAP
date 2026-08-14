"""One- versus two-component controls and mock-calibrated detection metrics.

A simple textbook chi-square p-value is not assumed for the one/two-component
mixture comparison. Instead, single-component mocks calibrate the empirical
false-positive probability of the chosen detection statistic, while injected
true two-component mocks calibrate recovery/completeness as a function of the
observed problem parameters.
"""

from __future__ import annotations

import numpy as np


def two_component_statistic(chi2_one: float, chi2_two: float) -> float:
    """Return ``chi2_one - chi2_two``; positive values favor two components."""
    return float(chi2_one - chi2_two)


def empirical_false_positive_probability(observed_statistic: float, single_component_mock_statistics: np.ndarray) -> float:
    """Fraction of true-one-component mocks producing >= observed statistic."""
    mocks = np.asarray(single_component_mock_statistics, dtype=float)
    mocks = mocks[np.isfinite(mocks)]
    if mocks.size == 0:
        raise ValueError("No finite single-component mock statistics supplied.")
    # Add-one smoothing avoids reporting an exact zero from a finite simulation
    # ensemble while preserving a transparent frequentist interpretation.
    exceed = np.count_nonzero(mocks >= observed_statistic)
    return float((exceed + 1) / (mocks.size + 1))


def recovery_probability(success_flags: np.ndarray) -> float:
    """Empirical recovery fraction for matched two-component injection mocks."""
    flags = np.asarray(success_flags, dtype=bool)
    if flags.size == 0:
        raise ValueError("No recovery trials supplied.")
    return float(np.mean(flags))
