from __future__ import annotations

import numpy as np
from astropy.table import Table

from crd_utils.covariance_calibration import (
    ModelComparison,
    choose_simplest_stable_model,
    compare_likelihood_cubes,
    select_representative_bins,
)


def _fake_bin_geometry(n_each_side=12):
    # 24 narrow bins along x/PA plus two off-axis distractors.
    x = np.r_[np.linspace(-12, -0.5, n_each_side), np.linspace(0.5, 12, n_each_side), -5, 5]
    y = np.r_[np.linspace(-0.15, 0.15, n_each_side), np.linspace(0.15, -0.15, n_each_side), 3.0, -3.0]
    n = x.size
    area = np.linspace(0.8, 4.0, n)
    tab = Table()
    tab["BIN_ID"] = np.arange(n, dtype=int)
    tab["AREA_ARCSEC2"] = area
    tab["X_GEOM_ARCSEC"] = x
    tab["Y_GEOM_ARCSEC"] = y
    tab["RH3_SN"] = np.linspace(50, 10, n)
    tab["NPIX_RH3"] = np.arange(1, n + 1)
    sigma = np.full(n, 80.0)
    # Strong off-center sigma peaks not guaranteed to coincide with six radial picks.
    sigma[4] = 170.0
    sigma[n_each_side + 7] = 180.0
    # simple stripe bin map so every bin has a pixel centroid
    bmap = np.tile(np.arange(n, dtype=int), (3, 1))
    return tab, bmap, sigma


def test_representative_bins_are_symmetric_radial_plus_optional_sigma_peaks():
    tab, bmap, sigma = _fake_bin_geometry()
    reps = select_representative_bins(
        tab, bmap, sigma, pa_deg=90.0, n_radial_total=12,
        corridor_diameter_factor=1.0, add_two_sigma_bins=True,
    )
    radial = np.asarray(reps["BASE_RADIAL_SELECTION"], dtype=bool)
    assert np.sum(radial) == 12
    sides = np.asarray(reps["PA_SIDE"], dtype=int)[radial]
    assert np.sum(sides < 0) == 6
    assert np.sum(sides > 0) == 6
    assert len(np.unique(np.asarray(reps["BIN_ID"], dtype=int))) == len(reps)
    assert 12 <= len(reps) <= 14


def _cube(center):
    shape = (5, 5, 3)
    idx = np.indices(shape)
    return ((idx[0]-center[0])**2 + (idx[1]-center[1])**2 + (idx[2]-center[2])**2).astype(float)


def test_requirement_b_agreement_is_grid_resolution_based():
    a = {1: _cube((2, 2, 1)), 2: _cube((2, 2, 1))}
    b = {1: _cube((3, 2, 1)), 2: _cube((2, 2, 1))}
    comp = compare_likelihood_cubes(a, b, max_cell_shift=1, delta_chi2_thresholds=(1.0, 4.0))
    assert comp.agree
    c = {1: _cube((4, 2, 1)), 2: _cube((2, 2, 1))}
    comp2 = compare_likelihood_cubes(a, c, max_cell_shift=1, delta_chi2_thresholds=(1.0, 4.0))
    assert not comp2.agree


def _cmp(agree):
    return ModelComparison(bool(agree), 0 if agree else 2, 0 if agree else 2, (), ())


def test_model_selection_uses_simplest_stable_candidate_and_does_not_default_to_m4():
    req = {"M1": False, "M2": True, "M3": True, "M4": True}
    pair = {
        ("M1", "M2"): _cmp(False), ("M1", "M3"): _cmp(False), ("M1", "M4"): _cmp(False),
        ("M2", "M3"): _cmp(True), ("M2", "M4"): _cmp(True), ("M3", "M4"): _cmp(True),
    }
    assert choose_simplest_stable_model(req, pair) == "M2"

    req_all = {m: True for m in ("M1", "M2", "M3", "M4")}
    pair_bad = {
        ("M1", "M2"): _cmp(False), ("M1", "M3"): _cmp(False), ("M1", "M4"): _cmp(False),
        ("M2", "M3"): _cmp(False), ("M2", "M4"): _cmp(False), ("M3", "M4"): _cmp(False),
    }
    import pytest
    with pytest.raises(RuntimeError, match="STABILITY_FAILURE"):
        choose_simplest_stable_model(req_all, pair_bad)
