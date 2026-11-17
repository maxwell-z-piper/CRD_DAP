"""Regression checks for Script 06 weight conventions and grouped recovery."""
import runpy
from pathlib import Path

import numpy as np
import pytest

from crd_utils import bl_population_products, populations, rh3_refit


def test_ssp_light_and_mass_conventions():
    meta = dict(ages=np.array([1., 9.]), metals=np.array([-.5, .1]),
                masses=np.array([[.8, .7], [.6, .5]]),
                luminosity=np.array([[2., 4.], [3., 5.]]))
    normalization = np.array([2., 4., 3., 5.])
    single, grid = populations.summarize_ssp_weights([0., 4., 0., 0.], normalization, meta)
    assert single['MASS_TO_LIGHT'] == pytest.approx(.7/4.)
    assert np.array_equal(grid['sfh_formed_mass'], [1., 0.])
    mixed, grid = populations.summarize_ssp_weights([2., 0., 0., 2.], normalization, meta)
    assert mixed['LIGHT_AGE_GYR'] == 5.
    assert 10**(mixed['LIGHT_LOGAGE_YR']-9.) == pytest.approx(3.)
    assert mixed['MASS_AGE_GYR'] == pytest.approx((.8 + .2*9)/1.)
    assert mixed['MASS_TO_LIGHT'] == pytest.approx(1./4.)
    assert grid['sfh_formed_mass'] == pytest.approx([1/1.4, .4/1.4])


def test_residual_scale_keeps_the_profile_minimum():
    setup = dict(n_stellar=1, stellar_normalization=np.array([2.]),
                 population_metadata=dict(ages=np.array([5.]), metals=np.array([0.]),
                                          masses=np.array([[.75]]), luminosity=np.array([[3.]])))
    residual = np.array([-.2, -.1, 0., .1, .2])
    profile = dict(chi2=np.array([3., 1., 2.]),
                   model_total=np.tile(np.ones(5), (3, 1)),
                   weights=np.tile([.5, .5], (3, 1)))
    spectrum = dict(galaxy=1.+residual, noise=np.ones(5), good=np.ones(5, dtype=bool))
    bl_population_products.finalize_profile(profile, setup, spectrum)
    assert profile['residual_scale'] == pytest.approx(.14826)
    assert profile['chi2_calibrated'] == pytest.approx(profile['chi2']/.14826**2)
    assert np.array_equal(profile['local_minimum_indices'], [1])
    assert profile['stellar_mass_fraction_A'] == .5


def test_group_move_on_saves_all_failed_fractions_in_one_decision(tmp_path, monkeypatch):
    script = Path(__file__).resolve().parents[1]/'scripts/06_fit_BL_two_component_populations.py'
    driver = runpy.run_path(str(script))
    f = rh3_refit.FRACTIONS
    pending = np.isin(f, [.01, .02, .99])
    profile = dict(bin_id=5, velocities=np.array([5., -5.]),
                   kinematics=np.array([5., 20., -5., 25.]),
                   fraction=f, pending=pending.copy(), moved_on=np.zeros(len(f), dtype=bool),
                   attempts=[[dict(reason='test failure' if p else '')] for p in pending])
    asked = []
    monkeypatch.setattr('builtins.input', lambda prompt: asked.append(prompt) or '3')
    import logging
    path = tmp_path/'bin_0005.npz'
    assert driver['recover'](profile, path, 'test', {}, {'gas': np.array([])},
                             True, logging.getLogger('test'))
    saved = rh3_refit.load_profile(path, 'test', 5, profile['velocities'])
    assert len(asked) == 1 and not saved['pending'].any()
    assert np.array_equal(saved['moved_on'], pending)
