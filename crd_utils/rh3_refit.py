"""Exact RH3 fits and their resumable fraction profiles.

Keep spectral fitting independent of the terminal and run directory so later
uncertainty realizations can reuse the same fixed-velocity experiment.
"""
from __future__ import annotations

import json
from pathlib import Path
import time

import numpy as np

from .ppxf_utils import fit_fixed_two_component_state

FRACTIONS = np.arange(1, 100, dtype=float) / 100
MAX_ATTEMPTS = 4


def fit_fraction(setup, spectrum, velocities, fraction, start):
    """Return an audit record and, for a scorable result, compact spectral arrays."""
    record = dict(start=list(map(float, start)), returned=False, status=None,
                  velocity=None, sigma=None, chi2=None, achieved_fraction=None,
                  flags=[], reason='')
    result = fit_fixed_two_component_state(
        templates_two_component=setup['templates'], component=setup['component'],
        galaxy=spectrum['galaxy'], noise=spectrum['noise'],
        goodpixels=np.flatnonzero(spectrum['good']), lam=setup['wavelength'],
        lam_temp=setup['template_wavelength'], velscale=setup['velscale'],
        velocity_a=velocities[0], velocity_b=velocities[1], fraction_a=fraction,
        start_sigma_a=start[0], start_sigma_b=start[1],
        sigma_bounds=setup['sigma_bounds'], degree=setup['degree'],
        mdegree=setup['mdegree'], regul=setup['regul'], keep_full=True,
    )
    if not result.success:
        if result.error_message.startswith('MemoryError:'):
            raise MemoryError(result.error_message)
        record['reason'] = result.error_message
        record['flags'].append('PPXF_EXCEPTION')
        return record, None

    record['returned'] = True
    record['status'] = int(result.pp.status)
    for key, values in (('velocity', result.velocity), ('sigma', result.sigma)):
        record[key] = [float(v) if np.isfinite(v) else None for v in values]
    if np.isfinite(result.chi2_total):
        record['chi2'] = float(result.chi2_total)
    if record['status'] <= 0:
        record['flags'].append('PPXF_NOT_CONVERGED')
        record['reason'] = f'pPXF optimizer returned status {record["status"]}'

    weights = result.weights
    polyweights = result.polyweights if result.polyweights is not None else np.empty(0)
    if not all(np.all(np.isfinite(v)) for v in (
        result.velocity, result.sigma, result.chi2_total, result.bestfit, weights, polyweights
    )):
        record['flags'].append('NONFINITE_FIT')
        record['reason'] = 'Nonfinite fitted parameters, score or spectrum'
        return record, None
    if not np.allclose(result.velocity, velocities, rtol=0, atol=1e-6):
        record['flags'].append('FIXED_VELOCITY_MISMATCH')
        record['reason'] = 'Returned velocities differ from the fixed Script 04 pair'
        return record, None

    component = setup['component']
    total = weights.sum()
    achieved = float(weights[component == 0].sum() / total) if total > 0 else None
    record['achieved_fraction'] = achieved
    if achieved is None or abs(achieved - fraction) > setup['fraction_tolerance']:
        record['flags'].append('FIXED_FRACTION_MISMATCH')
        record['reason'] = 'Returned stellar weights do not satisfy the imposed light fraction'
        return record, None

    lo, hi = setup['sigma_bounds']
    if np.any(result.sigma < lo - 1e-6) or np.any(result.sigma > hi + 1e-6):
        record['flags'].append('SIGMA_OUT_OF_BOUNDS')
        record['reason'] = 'Returned dispersions violate the inherited bounds'
        return record, None
    tol = setup['sigma_warning']
    if np.any((result.sigma <= lo + tol) | (result.sigma >= hi - tol)):
        record['flags'].append('SIGMA_BOUNDARY')

    # pPXF's matrix starts with additive polynomial columns, followed by the
    # LOSVD-convolved templates in component order. The template columns
    # already carry any fitted multiplicative continuum.
    n_poly = len(polyweights)
    matrix = result.pp.matrix[:, n_poly:n_poly + len(weights)]
    model_a = matrix[:, component == 0] @ weights[component == 0]
    model_b = matrix[:, component == 1] @ weights[component == 1]
    additive = result.pp.apoly if result.pp.apoly is not None else 0.
    if not np.allclose(model_a + model_b + additive, result.bestfit, rtol=1e-8, atol=1e-10):
        raise RuntimeError('pPXF component spectra and additive continuum do not reconstruct the total model')
    return record, dict(
        chi2=result.chi2_total, sigma=result.sigma, achieved_fraction=achieved,
        model_A=model_a, model_B=model_b, model_total=result.bestfit, weights=weights,
        additive_coefficients=polyweights,
        multiplicative_coefficients=(result.pp.mpolyweights
                                     if result.pp.mpolyweights is not None else np.empty(0)),
    )


def empty_profile(bin_id, velocities, setup):
    n_fraction = len(FRACTIONS)
    n_wave = len(setup['wavelength'])
    profile = dict(bin_id=bin_id, velocities=np.asarray(velocities), fraction=FRACTIONS.copy(),
                   pending=np.ones(n_fraction, dtype=bool),
                   moved_on=np.zeros(n_fraction, dtype=bool),
                   winning_attempt=np.full(n_fraction, -1, dtype=int),
                   attempts=[[] for _ in FRACTIONS], elapsed_seconds=0.)
    shapes = dict(chi2=(), sigma=(2,), achieved_fraction=(), model_A=(n_wave,),
                  model_B=(n_wave,), model_total=(n_wave,),
                  weights=(len(setup['component']),),
                  additive_coefficients=(max(0, setup['degree'] + 1),),
                  multiplicative_coefficients=(setup['mdegree'],))
    profile.update({key: np.full((n_fraction, *shape), np.nan) for key, shape in shapes.items()})
    return profile


def record_attempt(profile, index, record, products):
    """Retain the best scorable result even if another start fails."""
    attempts = profile['attempts'][index]
    attempts.append(record)
    if products is not None:
        previous = profile['chi2'][index]
        if not np.isfinite(previous) or products['chi2'] < previous:
            for name, values in products.items():
                profile[name][index] = values
            profile['winning_attempt'][index] = len(attempts) - 1
    profile['pending'][index] = bool(record['reason']) and len(attempts) < MAX_ATTEMPTS


def fit_initial_profile(setup, bin_id, spectrum, velocities):
    started = time.perf_counter()
    profile = empty_profile(bin_id, velocities, setup)
    start = (setup['baseline'], setup['baseline'])
    for index, fraction in enumerate(FRACTIONS):
        record, products = fit_fraction(setup, spectrum, velocities, fraction, start)
        record_attempt(profile, index, record, products)
    profile['elapsed_seconds'] = time.perf_counter() - started
    return profile


def initialize_worker(setup):
    # The templates and fit settings are shared by all bins handled by this
    # process. Send them once at process startup.
    global _WORKER_SETUP
    _WORKER_SETUP = setup


def fit_worker_bin(bin_id, spectrum, velocities):
    return fit_initial_profile(_WORKER_SETUP, bin_id, spectrum, velocities)


def save_profile(path, profile, fingerprint):
    """An atomic per-bin product also serves as its retry checkpoint."""
    path = Path(path)
    payload = {key: value for key, value in profile.items() if key != 'attempts'}
    payload.update(fingerprint=fingerprint, attempts_json=json.dumps(profile['attempts'], allow_nan=False))
    temp = path.with_suffix('.npz.tmp')
    with temp.open('wb') as stream:
        np.savez_compressed(stream, **payload)
    temp.replace(path)


def load_profile(path, fingerprint, bin_id, velocities):
    with np.load(path, allow_pickle=False) as saved:
        profile = dict(saved)
    if (str(profile.pop('fingerprint')) != fingerprint or int(profile['bin_id']) != bin_id
            or not np.array_equal(profile['velocities'], velocities)
            or not np.array_equal(profile['fraction'], FRACTIONS)):
        raise ValueError(f'Incompatible Script 05 checkpoint: {path}')
    profile['attempts'] = json.loads(str(profile.pop('attempts_json')))
    if len(profile['attempts']) != len(FRACTIONS) or any(
        not 1 <= len(a) <= MAX_ATTEMPTS for a in profile['attempts']
    ):
        raise ValueError(f'Malformed Script 05 attempt history: {path}')
    return profile


def best_fraction(profile):
    """Return the best available finite profile index, with stable ties."""
    chi2 = profile['chi2']
    return int(np.nanargmin(chi2)) if np.any(np.isfinite(chi2)) else None
