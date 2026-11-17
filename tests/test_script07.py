"""Regression tests for Script 07 science aggregation and resumable draws."""
from contextlib import nullcontext
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from astropy.table import Table

from crd_utils import script07_fits as fits, script07_products as products
from crd_utils.global_disk_io import load_json
from tests.test_global_disk import problem_fixture


def driver():
    path = Path(__file__).absolute().parents[1]/'scripts/07_uncertainties_and_final_maps.py'
    spec = spec_from_file_location('script07_test_driver', path)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_seed_is_repeatable_and_independent_across_draws_and_arms():
    source = dict(good=np.array([[True, True, False]]), noise=np.ones((1, 3)), model=np.zeros((1, 3)))
    fits.initialize({'BL': {'gas_by_bin': [False]}}, {'RH3': source, 'BL': source}, {})
    red = fits.mock_spectrum('RH3', 0, 7, 4)['galaxy']
    np.testing.assert_array_equal(red, fits.mock_spectrum('RH3', 0, 7, 4)['galaxy'])
    assert not np.array_equal(red, fits.mock_spectrum('BL', 0, 7, 4)['galaxy'])
    assert not np.array_equal(red, fits.mock_spectrum('RH3', 0, 7, 5)['galaxy'])
    assert red[-1] == 0.


@pytest.mark.parametrize('mode, count, reason', [('smooth', 29, []), ('flat', 99, ['FLAT']),
                                               ('missing', 99, ['MISSING_FINE_SCORE'])])
def test_adaptive_scan_and_unattended_failed_point(monkeypatch, mode, count, reason):
    source = dict(good=np.ones((1, 5), bool), noise=np.ones((1, 5)), model=np.ones((1, 5)))
    fits.initialize({'RH3': {'baseline': 80.}}, {'RH3': source}, {})
    attempted = []

    def fit(setup, spectrum, vel, f, start):
        attempted.append(f)
        if mode == 'missing' and f == .52:
            return dict(chi2=None, flags=['NONFINITE_FIT'], status=None, reason='Test failure'), None
        score = 10. if mode == 'flat' else 100*(f-.53)**2
        return dict(chi2=score, flags=[], status=1, reason=''), dict(
            chi2=score, sigma=np.array([40., 50.]), achieved_fraction=f)

    monkeypatch.setattr(fits.rh3_refit, 'fit_fraction', fit)
    profile, pending = fits.fraction_worker(('RH3', 0, [110., -30.], 7, 1, True))
    assert profile['complete'] and not pending
    assert len(set(attempted)) == len(attempted) == count
    assert profile['fallback'] == reason
    assert profile['winner']['fraction'] == (.01 if mode == 'flat' else .53)
    if mode == 'missing':
        assert fits.audit_profile(profile)['missing_fractions'] == [.52]
        assert fits.audit_profile(profile)['moved_on'] == [.52]


def test_nominal_light_restoration_is_flux_weighted_and_mass_is_surviving():
    setup = dict(stellar_normalization=np.array([2., 1.]),
                 population_metadata=dict(ages=np.array([1., 9.]), metals=np.array([0.]),
                                          masses=np.array([[.5], [.8]]), luminosity=np.array([[2.], [1.]])))
    weights = np.asarray([[[8., 0.], [0., 1.]], [[0., 2.], [1., 0.]]])
    result, grid, n = products.aperture_statistics(weights, setup)
    assert n == 2
    assert grid.shape == (2, 2, 1)
    assert result['A_LIGHT_AGE_GYR'] == pytest.approx((8+18)/10)
    assert result['A_MASS_AGE_GYR'] == pytest.approx((2+1.6*9)/3.6)
    assert result['A_LIGHT_FRACTION'] == pytest.approx(10/12)
    assert result['A_MASS_TO_LIGHT'] == pytest.approx(3.6/10)


def test_resume_skips_saved_stages_and_cleans_only_committed_draws(tmp_path, monkeypatch):
    module = driver()
    p, truth = problem_fixture()
    va, vb, *_ = p.predict(truth)
    ids = dict(RH3=np.arange(6), BL=np.arange(6))
    setup = dict(n_stellar=1, stellar_normalization=np.array([1.]),
                 population_metadata=dict(ages=np.array([2.]), metals=np.array([-.1]),
                                          masses=np.array([[.8]]), luminosity=np.array([[1.]])))
    nominal = SimpleNamespace(ids=ids, setup={'RH3': {}, 'BL': setup}, problem=p, fit={'name': 'full_psf0_smooth0'},
                              spectra={'BL': {'normalization_scale': np.ones(6)}})
    (tmp_path/'draws').mkdir()
    calls = []

    def rh_grid(payload):
        bid, _, draw = payload
        calls.append(('GRID', draw, bid))
        return bid, p.chi2[bid], dict(bin_id=bid, minimum=0., elapsed_seconds=0.)

    def disk_fit(saved, cube):
        calls.append(('DISK',))
        return dict(theta=truth.tolist(), bin_ids=list(range(6)), va=va.tolist(), vb=vb.tolist(), flags=[])

    def fraction(payload):
        arm, bid, kin, seed, draw, unattended = payload
        calls.append((arm, draw, bid))
        if arm == 'BL' and draw == 2 and bid == 3 and calls.count(('BL', draw, bid)) == 1:
            raise KeyboardInterrupt
        profile = fits.empty_profile(arm, bid, kin)
        profile['complete'] = True
        if draw != 1 or bid != 0:
            profile['winner'] = dict(index=49, chi2=1., flags=[], status=1,
                                     fraction=.50, achieved=.50,
                                     **({'sigma': [40., 50.]} if arm == 'RH3' else
                                        {'weights': [.5, .5], 'gas': [0., 40.]}))
        return profile, []

    monkeypatch.setattr(fits, 'grid_worker', rh_grid)
    monkeypatch.setattr(fits, 'coherent_fit', disk_fit)
    monkeypatch.setattr(fits, 'fraction_worker', fraction)
    monkeypatch.setattr(fits, 'mock_spectrum', lambda *a: {})
    monkeypatch.setattr(module, 'bin_pool', lambda *a: nullcontext(None))
    (tmp_path/'metadata').mkdir()
    with pytest.raises(module.PendingReview, match='1 consecutive'):
        module.collect(tmp_path, nominal, 1, 7, 2, 1, True, {'RH3': None, 'BL': None})
    assert not list((tmp_path/'work').iterdir())
    assert load_json(tmp_path/'draws/draw_0001.json')['reasons'] == ['RH3_UNUSABLE_FRACTION']
    with pytest.raises(KeyboardInterrupt):
        module.collect(tmp_path, nominal, 1, 7, 2, 1, True, {'RH3': None, 'BL': None})
    assert (tmp_path/'work/draw_0002/04.json').exists()
    assert (tmp_path/'work/draw_0002/06/bin_0002.json').exists()
    accepted, attempts = module.collect(tmp_path, nominal, 1, 7, 2, 1, True, {'RH3': None, 'BL': None})
    assert [record['draw'] for record in accepted] == [2, 3]
    assert len(attempts) == 3
    assert calls.count(('BL', 2, 2)) == calls.count(('GRID', 2, 0)) == 1
    assert calls.count(('DISK',)) == 3
    assert not list((tmp_path/'work').iterdir())


def test_checkpointed_retry_replays_the_selected_start_once(tmp_path, monkeypatch):
    module = driver()
    source = dict(good=np.ones((1, 5), bool), noise=np.ones((1, 5)), model=np.ones((1, 5)))
    setup = {'baseline': 80.}
    fits.initialize({'RH3': setup}, {'RH3': source}, {})
    profile = fits.empty_profile('RH3', 0, [100., -100.])
    used = []

    def fit(setup, spectrum, vel, f, start):
        if f == .01 and list(start) == [80., 80.]:
            return dict(start=list(start), chi2=None, status=None, reason='Failed start', flags=[]), None
        if f == .01:
            used.append(list(start))
            if len(used) == 1:
                raise KeyboardInterrupt
        score = 10*(f-.4)**2
        return dict(start=list(start), chi2=score, status=1, reason='', flags=[]), dict(
            chi2=score, achieved_fraction=f, sigma=np.array(start))

    monkeypatch.setattr(fits.rh3_refit, 'fit_fraction', fit)
    profile, pending = fits.advance(profile, 42, 1, False)
    assert pending == [0]
    path = tmp_path/'profile.json'
    selected = SimpleNamespace(choose_retry=lambda *a: [93., 102.])
    monkeypatch.setattr(module.sys.stdin, 'isatty', lambda: True)

    def checkpoint(state):
        module.save_json(path, state)

    with pytest.raises(KeyboardInterrupt):
        module._resolve_recovery(profile, pending, selected, setup, 42, 1, checkpoint, False)
    restored = load_json(path)
    assert restored['retry_queue'] == [[0, [93., 102.]]]
    result = module._resolve_recovery(restored, [], None, setup, 42, 1, checkpoint, False)
    assert result['complete']
    assert used == [[93., 102.], [93., 102.]]
    assert len(result['attempts'][0]) == 2


def test_plot_revision_resumes_only_the_unchanged_experiment():
    from crd_utils.bl_population_io import fingerprint
    from crd_utils.global_disk_io import file_hash
    from crd_utils.script07_resume import _COMPATIBLE_REVISIONS, matches_checkpoint

    old_driver, old_products = _COMPATIBLE_REVISIONS[0]
    original = dict(driver_sha256=old_driver, code={'script07_products.py': old_products, 'fits.py': 'unchanged'},
                    seed=7, target_draws=100)
    saved = dict(contract=original, fingerprint=fingerprint(original))
    root = Path(__file__).absolute().parents[1]
    code = original['code'] | {'script07_products.py': file_hash(root/'crd_utils/script07_products.py')}
    revision = dict(original, driver_sha256=file_hash(root/'scripts/07_uncertainties_and_final_maps.py'), code=code)

    def state(contract):
        return dict(contract=contract, fingerprint=fingerprint(contract))

    assert matches_checkpoint(saved, state(revision))
    for prior_driver, prior_exporter in _COMPATIBLE_REVISIONS[1:]:
        prior = dict(original, driver_sha256=prior_driver, code=code | {'script07_products.py': prior_exporter})
        assert matches_checkpoint(state(prior), state(revision))
    assert matches_checkpoint(state(revision), state(revision))
    assert not matches_checkpoint(saved, state(dict(revision, seed=8)))
    assert not matches_checkpoint(saved, state(dict(revision, driver_sha256='changed fit driver')))
    assert not matches_checkpoint(saved, state(dict(revision, code=code | {'fits.py': 'changed'})))


def test_combined_rh3_dispersion_includes_the_separation_of_two_disks():
    moments = products.combined_rh3_moments(dict(
        VA_GLOBAL_KMS=np.array([100., 100., 100.]),
        VB_GLOBAL_KMS=np.array([-100., -100., -100.]),
        SIGMA_A_KMS=np.array([40., 40., 40.]),
        SIGMA_B_KMS=np.array([40., 40., 40.]),
        FA_RH3_ACHIEVED=np.array([.5, 1., np.nan]),
    ))
    np.testing.assert_allclose(moments['V_COMBINED_KMS'], [0., 100., np.nan])
    np.testing.assert_allclose(moments['SIGMA_COMBINED_KMS'],
                               [np.hypot(100., 40.), 40., np.nan])
