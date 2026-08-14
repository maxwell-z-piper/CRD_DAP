import numpy as np

from crd_utils.sampling import (
    sample_discrete_states,
    effective_sample_size,
    check_quantile_convergence,
)


def test_discrete_sampling_favors_high_probability_cell():
    rng = np.random.default_rng(1)
    draws = sample_discrete_states(np.array([0.9, 0.1]), 2000, rng=rng)
    assert np.mean(draws == 0) > 0.85


def test_effective_sample_size_uniform():
    assert np.isclose(effective_sample_size(np.ones(4)), 4.0)


def test_quantile_convergence_passes_small_change():
    previous = np.array([1.0, 2.0, 3.0])
    current = np.array([1.01, 2.01, 3.01])
    result = check_quantile_convergence(previous, current, fraction_of_ci=0.05, absolute_floor=0.005)
    assert result.converged
