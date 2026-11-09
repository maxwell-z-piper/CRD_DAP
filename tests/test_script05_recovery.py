"""Exercise retry accounting, persistence and selection of flagged finite fits."""
import importlib.util
import logging
from pathlib import Path

import numpy as np
import pytest

from crd_utils import rh3_refit


def driver():
    path = Path(__file__).resolve().parents[1] / 'scripts/05_extract_global_RH3_solution.py'
    spec = importlib.util.spec_from_file_location('script05_recovery_test', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def attempt(start=(102., 102.), reason='', chi2=1.):
    return dict(start=list(start), returned=True, status=0 if reason else 2,
                velocity=[-40., 60.], sigma=[50., 80.], chi2=chi2,
                achieved_fraction=.01, flags=['PPXF_NOT_CONVERGED'] if reason else [],
                reason=reason)


def checkpoint(tmp_path):
    setup = dict(wavelength=np.arange(5), component=np.array([0, 1]),
                 degree=-1, mdegree=0, baseline=102., sigma_bounds=[5., 250.])
    profile = rh3_refit.empty_profile(0, [-40., 60.], setup)
    for i in range(99):
        rh3_refit.record_attempt(profile, i, attempt(), dict(chi2=3.))
    profile['attempts'][0] = [attempt(reason='optimizer limit')]
    profile['pending'][:2] = True
    profile['chi2'][0] = 1.
    profile['chi2'][1] = np.nan
    profile['winning_attempt'][1] = -1
    path = tmp_path / 'bin_0000.npz'
    rh3_refit.save_profile(path, profile, 'test')
    inputs = dict(setup=setup, spectra={k: np.ones((1, 5)) for k in ('galaxy', 'noise', 'good')})
    return path, inputs


def load(path):
    return rh3_refit.load_profile(path, 'test', 0, [-40., 60.])


def test_unattended_then_interactive_review_preserves_best_flagged_result(tmp_path, monkeypatch, capsys):
    module = driver()
    path, inputs = checkpoint(tmp_path)
    profile = load(path)
    original = path.read_bytes()
    assert not module.recover(profile, path, 'test', inputs, False, logging.getLogger(__name__))
    assert path.read_bytes() == original

    entries = iter(['2', '87 117', '3'])
    monkeypatch.setattr('builtins.input', lambda prompt: next(entries))
    monkeypatch.setattr(rh3_refit, 'fit_fraction',
                        lambda setup, spectrum, velocity, fraction, start:
                        (attempt(start, chi2=2.), dict(chi2=2.)))
    assert module.recover(profile, path, 'test', inputs, True, logging.getLogger(__name__))
    result = load(path)
    assert not result['pending'].any()
    assert result['moved_on'][1]
    assert result['attempts'][0][1]['start'] == [87., 117.]
    assert result['winning_attempt'][0] == 0
    assert result['chi2'][0] == 1.
    assert np.isnan(result['chi2'][1])
    assert '87, 117' in capsys.readouterr().out


def test_interrupted_retry_counts_and_same_pair_is_latest(tmp_path, monkeypatch):
    module = driver()
    path, inputs = checkpoint(tmp_path)
    profile = load(path)
    entries = iter(['2', '88 116'])
    monkeypatch.setattr('builtins.input', lambda prompt: next(entries))

    def interrupted(*args):
        raise MemoryError('test interruption')

    monkeypatch.setattr(rh3_refit, 'fit_fraction', interrupted)
    with pytest.raises(MemoryError):
        module.recover(profile, path, 'test', inputs, True, logging.getLogger(__name__))
    profile = load(path)
    assert len(profile['attempts'][0]) == 2
    assert profile['attempts'][0][-1]['flags'] == ['INTERRUPTED_ATTEMPT']
    assert profile['pending'][0]

    starts = []

    def fail(setup, spectrum, velocity, fraction, start):
        starts.append(start)
        return attempt(start, 'optimizer limit'), None

    entries = iter(['1', '1', '3'])
    monkeypatch.setattr('builtins.input', lambda prompt: next(entries))
    monkeypatch.setattr(rh3_refit, 'fit_fraction', fail)
    module.recover(profile, path, 'test', inputs, True, logging.getLogger(__name__))
    result = load(path)
    assert starts == [[88., 116.], [88., 116.]]
    assert len(result['attempts'][0]) == 4
    assert not result['pending'].any()
