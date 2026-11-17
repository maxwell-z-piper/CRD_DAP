"""Regression checks for sharing search states without mixing scientific models."""
import importlib.util
import logging
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from crd_utils import global_disk
from crd_utils.global_disk_io import load_json, save_json
from tests.test_global_disk import problem_fixture


def driver():
    path=Path(__file__).resolve().parents[1]/'scripts/04_fit_global_two_disk_model.py'
    spec=importlib.util.spec_from_file_location('script04_shared_start_tests',path)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_shared_best_state_survives_bad_local_refinement(monkeypatch,caplog):
    problem,truth=problem_fixture()

    def explore(fun,bounds,*,init,**kwargs):
        return SimpleNamespace(x=init[0],population=init[:5],success=True,message='test search')

    def worsen(fun,state,**kwargs):
        bad=state.copy()
        bad[5:]=0.
        return SimpleNamespace(x=bad,fun=fun(bad),success=False,message='test local failure')

    monkeypatch.setattr(global_disk,'differential_evolution',explore)
    monkeypatch.setattr(global_disk,'minimize',worsen)
    with caplog.at_level(logging.INFO,logger='crd_utils.global_disk'):
        result=global_disk.fit_disk(problem,0.,global_iterations=0,local_starts=1,
                                    local_iterations=1,additional_starts=[truth])
    np.testing.assert_array_equal(result['theta'],truth)
    assert result['objective']==pytest.approx(problem.evaluate(truth,0.))
    assert result['additional_start_count']==1
    assert 'LOCAL_OPTIMIZER_NOT_CONVERGED' in result['flags']
    assert 'KINEMATIC_RESTART_STALLED' in result['flags']
    audit = result['kinematic_refinement']
    assert audit['stop_reason'] == 'no_improvement'
    assert len(audit['restarts']) == 1
    assert not audit['restarts'][0]['accepted']
    assert audit['restarts'][0]['node_span_kms'] is None
    assert 'geometry=[east0_arcsec,north0_arcsec' in caplog.text


@pytest.mark.parametrize(('los_scale', 'limit', 'expected_stop', 'expected_passes'), [
    (0.02, 10, 'stable', 6),
    (0.20, 6, 'restart_limit', 6),
])
def test_powell_restarts_use_full_kinematic_window(
        monkeypatch, los_scale, limit, expected_stop, expected_passes):
    reference, initial = problem_fixture()
    problem = SimpleNamespace(bounds=reference.bounds, nodes=reference.nodes,
                              bin_ids=reference.bin_ids, selected=reference.selected,
                              nnode=reference.nnode)

    def evaluate(state, coefficient, *, details=False):
        speed_change = state[5] - initial[5]
        score = 100. - speed_change
        if not details:
            return score
        nbin = len(problem.selected)
        return dict(objective=score, data_score=score, penalty=0.,
                    va=np.full(nbin, los_scale * speed_change), vb=np.zeros(nbin),
                    radial_padding_fraction=np.zeros(nbin))

    problem.evaluate = evaluate
    population = np.tile(initial, (5, 1))
    monkeypatch.setattr(global_disk, 'starting_models', lambda *args, **kwargs: population)
    monkeypatch.setattr(global_disk, 'differential_evolution',
                        lambda *args, **kwargs: SimpleNamespace(
                            x=initial, population=population, success=True, message='test search'))
    steps = iter([0., 2., 4., 6., 7.25, 8., 8.5])  # First call is the original local fit.

    def improve(fun, state, **kwargs):
        new_state = initial.copy()
        new_state[5] += next(steps)
        assert fun(new_state) <= fun(state)
        return SimpleNamespace(x=new_state, fun=fun(new_state), success=True, message='test local')

    monkeypatch.setattr(global_disk, 'minimize', improve)
    result = global_disk.fit_disk(problem, 0., local_starts=1, powell_max_restarts=limit)
    audit = result['kinematic_refinement']
    assert audit['stop_reason'] == expected_stop
    assert len(audit['restarts']) == expected_passes
    assert audit['restarts'][1]['node_span_kms'] is None
    assert audit['restarts'][2]['node_span_kms'] == pytest.approx(6.)
    assert audit['restarts'][-1]['node_span_kms'] == pytest.approx(2.5)
    assert audit['restarts'][-1]['los_span_kms'] == pytest.approx(2.5 * los_scale)
    assert result['objective'] == pytest.approx(91.5)
    assert result['theta'] == result['alternatives'][0]['theta']
    np.testing.assert_array_equal(result['alternatives'][1]['theta'], initial)
    assert ('KINEMATIC_RESTART_LIMIT_REACHED' in result['flags']) == (expected_stop == 'restart_limit')


def test_additional_starts_must_obey_the_target_problem():
    problem,truth=problem_fixture()
    with pytest.raises(ValueError,match='complete parameter vector'):
        global_disk.fit_disk(problem,0.,additional_starts=[truth[:-1]])
    outside=truth.copy();outside[5]=201.
    with pytest.raises(ValueError,match='bounds and likelihood support'):
        global_disk.fit_disk(problem,0.,additional_starts=[outside])


def test_dominated_fit_is_refined_with_compatible_states_and_resume_is_stable(tmp_path,monkeypatch):
    module=driver()
    problem,truth=problem_fixture()
    info=dict(extent='full',A=1.,label=dict(review_required=False))
    fingerprint='one-scientific-experiment'
    coefficient=.001
    bad=truth.copy();bad[5:]+=20.
    group=module.fit_group(problem,info)

    def record(name,state,c,**extra):
        detail=problem.evaluate(state,c,details=True)
        return dict(name=name,fingerprint=fingerprint,coefficient=c,theta=state.tolist(),
                    nodes=problem.nodes.tolist(),bin_ids=problem.bin_ids.tolist(),fit_group=group,
                    flags=[],alternatives=[dict(theta=state.tolist())],elapsed_seconds=0.,
                    **{k:v.tolist() if isinstance(v,np.ndarray) else v for k,v in detail.items()},**extra)

    saved=record('current',bad,coefficient)
    save_json(tmp_path/'fits/current.json',saved)
    save_json(tmp_path/'fits/other_coefficient.json',record('other_coefficient',truth,0.))
    incompatible=record('wrong_psf',truth,0.)
    incompatible['fit_group']=dict(group,psf_enabled=True)
    incompatible['theta']=[999.]  # Must be excluded before parameter validation.
    save_json(tmp_path/'fits/wrong_psf.json',incompatible)
    old=record('wrong_fingerprint',truth,0.)
    old['fingerprint']='different-experiment'
    old['theta']=[999.]
    save_json(tmp_path/'fits/wrong_fingerprint.json',old)
    states,sources=module.shared_candidates(problem,info,tmp_path,fingerprint)
    assert len(states)==2
    assert sources==['fits/current.json','fits/other_coefficient.json']

    def fit(payload):
        target,c,options,logpath,sources=payload
        assert any(np.array_equal(s,truth) for s in options['additional_starts'])
        result=record('current',truth,c)
        result['additional_start_sources']=sources
        return result

    monkeypatch.setattr(module,'_fit_job',fit)
    cfg=SimpleNamespace(SCRIPT04_RANDOM_SEED=1,SCRIPT04_POPULATION_FACTOR=1,
                        SCRIPT04_GLOBAL_ITERATIONS=0,SCRIPT04_LOCAL_STARTS=1,
                        SCRIPT04_LOCAL_ITERATIONS=1,SCRIPT04_OPTIMIZER_TOLERANCE=1e-5,
                        SCRIPT04_SHARED_REFINEMENT_PASSES=2)
    args=SimpleNamespace(workers=1,resume=False)
    jobs=[('current',problem,coefficient,info)]
    result=module.refine_shared_fits(jobs,{'current':saved},cfg,args,tmp_path,fingerprint,
                                     logging.getLogger('shared-start-test'))
    assert result['current']['objective']<saved['objective']
    assert result['current']['shared_start_check']['complete']
    assert list((tmp_path/'metadata/shared_start_history').glob('*_current.json'))
    # Replacing a fit must not erase an earlier useful candidate from the pool.
    retained,_=module.shared_candidates(problem,info,tmp_path,fingerprint)
    assert any(np.array_equal(s,bad) for s in retained)
    checkpoint=(tmp_path/'fits/current.json').read_bytes()
    args.resume=True
    module.refine_shared_fits(jobs,{'current':load_json(tmp_path/'fits/current.json')},cfg,args,
                              tmp_path,fingerprint,logging.getLogger('shared-start-test'))
    assert (tmp_path/'fits/current.json').read_bytes()==checkpoint
