import numpy as np

from crd_utils.likelihood import (
    delta_chi2,
    relative_likelihood_from_delta_chi2,
    normalize_weights,
    likelihood_mass_cell_counts,
    basin_mask_for_anchor,
)


def test_delta_chi2_and_relative_likelihood():
    chi2 = np.array([10.0, 12.0, 14.0])
    d = delta_chi2(chi2)
    assert np.allclose(d, [0.0, 2.0, 4.0])
    l = relative_likelihood_from_delta_chi2(d)
    assert np.allclose(l, [1.0, np.exp(-1.0), np.exp(-2.0)])


def test_normalize_weights():
    p = normalize_weights(np.array([1.0, 2.0, 1.0]))
    assert np.isclose(p.sum(), 1.0)
    assert np.allclose(p, [0.25, 0.5, 0.25])


def test_likelihood_mass_cell_counts():
    p = np.array([0.6, 0.3, 0.1])
    result = likelihood_mass_cell_counts(p, masses=(0.5, 0.9, 1.0))
    assert result[0.5] == 1
    assert result[0.9] == 2
    assert result[1.0] == 3


def test_basin_mask_for_anchor_separates_two_minima():
    # Simple 1-D landscape represented as an array. The two minima at indices
    # 1 and 5 are separated by a ridge. The anchor near index 0 should select
    # the left-hand basin.
    chi2 = np.array([2.0, 0.0, 1.0, 5.0, 2.0, 0.5, 2.0])
    mask, minimum = basin_mask_for_anchor(chi2, (0,), connectivity="full")
    assert minimum == (1,)
    assert mask[0]
    assert mask[1]
    assert not mask[5]
