"""Numerical workers and the compact fraction search used by Script 07."""
from __future__ import annotations

from dataclasses import replace
import signal
import time

import numpy as np

from . import global_disk, populations, ppxf_grid, rh3_refit

FRACTIONS = rh3_refit.FRACTIONS
COARSE = np.r_[np.arange(0, 99, 5), 98]
BL_STATS = ('LIGHT_AGE_GYR', 'LIGHT_LOGAGE_YR', 'LIGHT_METAL_DEX',
            'MASS_AGE_GYR', 'MASS_LOGAGE_YR', 'MASS_METAL_DEX', 'MASS_TO_LIGHT')
REJECT_FRACTION = 0.10


def initialize(setups, spectra, grid):
    global _SETUP, _SPECTRA, _GRID
    _SETUP, _SPECTRA, _GRID = setups, spectra, grid


def initialize_worker(*args):
    # The parent handles Ctrl-C and owns worker shutdown.
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    initialize(*args)


def mock_spectrum(arm, bid, seed, draw):
    saved = _SPECTRA[arm]
    good = saved['good'][bid]
    noise = saved['noise'][bid]
    galaxy = saved['model'][bid].copy()
    rng = np.random.default_rng(np.random.SeedSequence([seed, draw, int(arm == 'BL'), bid]))
    if draw > 0:
        galaxy[good] += rng.standard_normal(np.count_nonzero(good)) * noise[good]
    spectrum = dict(galaxy=galaxy, noise=noise, good=good)
    if arm == 'BL':
        spectrum['gas'] = _SETUP['BL']['gas_by_bin'][bid]
    return spectrum


def grid_worker(payload):
    bid, seed, draw = payload
    started = time.perf_counter()
    spectrum = mock_spectrum('RH3', bid, seed, draw)
    setup = _SETUP['RH3']
    cube = ppxf_grid.build_rh3_likelihood_cube(
        templates_two_component=setup['templates'], component=setup['component'],
        galaxy=spectrum['galaxy'], noise=spectrum['noise'], goodpixels=np.flatnonzero(spectrum['good']),
        lam=setup['wavelength'], lam_temp=setup['template_wavelength'], velscale=setup['velscale'],
        degree=setup['degree'], mdegree=setup['mdegree'], regul=setup['regul'],
        sigma_bounds=setup['sigma_bounds'], sigma_boundary_tolerance_kms=setup['sigma_warning'],
        **_GRID)
    chi = np.where(np.isfinite(cube.chi2_total), cube.chi2_total, np.inf)
    valid = cube.best_index is not None
    minimum = float(chi[cube.best_index]) if valid else None
    best = [float(axis[i]) for axis, i in zip((cube.va_grid, cube.vb_grid, cube.fa_grid),
                                             cube.best_index)] if valid else None
    audit = dict(bin_id=bid, minimum=minimum, best=best, failed_states=cube.n_failures,
                 ppxf_failures=cube.n_ppxf_failures,
                 velocity_mismatches=cube.n_fixed_velocity_mismatch,
                 sigma_boundaries=cube.n_sigma_boundary, failure_messages=cube.failure_message_counts,
                 elapsed_seconds=time.perf_counter()-started)
    return bid, (chi-minimum if valid else chi), audit


def coherent_fit(nominal, cubes):
    problem = replace(nominal.problem, chi2=cubes)
    try:
        return global_disk.fit_disk(problem, nominal.fit['coefficient'], **nominal.optimizer)
    except (RuntimeError, ValueError) as exc:
        # Let unexpected code/configuration problems stop the run and be fixed.
        rejected_messages = ('No supported initial disk model', 'No finite global disk solution')
        if not str(exc).startswith(rejected_messages):
            raise
        return dict(failure=str(exc), flags=['NO_VALID_COHERENT_MODEL'])


def select_window(scores):
    """Return the bracket about the coarse winner, or the full grid and reasons."""
    y = scores[COARSE]
    if not np.isfinite(y).all():
        return np.arange(99), ['MISSING_COARSE_SCORE']
    shifted = y-y.min()
    span = float(shifted.max())
    x = FRACTIONS[COARSE]
    quadratic = np.polynomial.polynomial.polyfit(x, shifted, 2)
    residual = np.max(np.abs(np.polynomial.polynomial.polyval(x, quadratic)-shifted))
    minima = (y[1:-1] < y[:-2]) & (y[1:-1] < y[2:])
    fallback = []
    if span <= 2:
        fallback.append('FLAT')
    if residual > max(0.5, 0.06*span):
        fallback.append('NONQUADRATIC')
    if np.count_nonzero(minima) >= 2:
        fallback.append('COMPETING_MINIMA')
    if fallback:
        return np.arange(99), fallback
    j = int(np.argmin(y))
    return np.arange(COARSE[max(0, j-1)], COARSE[min(len(COARSE)-1, j+1)]+1), []


def empty_profile(arm, bid, kin):
    return dict(arm=arm, bin_id=bid, kin=list(kin), attempts=[[] for _ in FRACTIONS],
                scores=[None]*len(FRACTIONS), moved_on=[], fallback=[], window=None,
                winner=None, retry_queue=[], complete=False, elapsed_seconds=0.)


def _fit(profile, index, start, setup, spectrum):
    started = time.perf_counter()
    arm = profile['arm']
    method = rh3_refit.fit_fraction if arm == 'RH3' else populations.fit_fraction
    record, result = method(setup, spectrum, profile['kin'], FRACTIONS[index], start)
    history = profile['attempts'][index]
    history.append(record)
    profile['elapsed_seconds'] += time.perf_counter()-started
    previous = profile['scores'][index]
    if result is not None and (previous is None or result['chi2'] < previous):
        score = float(result['chi2'])
        profile['scores'][index] = score
        winner = profile['winner']
        if winner is None or (score, index) < (winner['chi2'], winner['index']):
            product = dict(chi2=score, index=int(index), flags=record['flags'], status=record['status'],
                           fraction=float(FRACTIONS[index]), achieved=float(result['achieved_fraction']))
            if arm == 'RH3':
                product['sigma'] = result['sigma'].tolist()
            else:
                product['weights'] = result['weights'][:2*setup['n_stellar']].tolist()
                product['gas'] = result['gas_kinematics'].tolist()
            profile['winner'] = product


def pending_indices(profile, indices):
    moved = set(profile['moved_on'])
    return [i for i in indices if profile['attempts'][i][-1]['reason']
            and len(profile['attempts'][i]) < rh3_refit.MAX_ATTEMPTS and i not in moved]


def run_retry_queue(profile, setup, spectrum, checkpoint):
    """Persist a chosen retry before launching it, including grouped BL retries."""
    while profile['retry_queue']:
        index, start = profile['retry_queue'][0]
        _fit(profile, index, start, setup, spectrum)
        profile['retry_queue'].pop(0)
        checkpoint(profile)


def advance(profile, seed, draw, non_interactive):
    """Advance until finished or an initial trial requires a user decision."""
    arm, bid = profile['arm'], profile['bin_id']
    setup = _SETUP[arm]
    spectrum = mock_spectrum(arm, bid, seed, draw)
    initial = [setup['baseline']]*2 if arm == 'RH3' else setup['gas_seed']

    def scan(indices):
        for i in indices:
            if not profile['attempts'][i]:
                _fit(profile, i, initial, setup, spectrum)
        pending = pending_indices(profile, indices)
        if non_interactive:
            profile['moved_on'].extend(pending)
            return []
        return pending

    if profile['complete']:
        return profile, []
    if profile['window'] is None and setup.get('full_scan', False):
        profile['window'] = list(range(99))
        profile['fallback'] = ['REQUESTED_FULL_SCAN']
    if profile['window'] is None:
        pending = scan(COARSE)
        if pending:
            return profile, pending
        window, profile['fallback'] = select_window(np.asarray(profile['scores'], dtype=float))
        profile['window'] = window.tolist()
    pending = scan(profile['window'])
    if pending:
        return profile, pending
    if not profile['fallback'] and any(profile['scores'][i] is None for i in profile['window']):
        profile['fallback'] = ['MISSING_FINE_SCORE']
        profile['window'] = list(range(99))
        return advance(profile, seed, draw, non_interactive)
    profile['complete'] = True
    return profile, []


def fraction_worker(payload):
    arm, bid, kin, seed, draw, non_interactive = payload
    return advance(empty_profile(arm, bid, kin), seed, draw, non_interactive)


def audit_profile(profile):
    tried = [i for i, history in enumerate(profile['attempts']) if history]
    absent = [i for i in tried if profile['scores'][i] is None]
    issues = {str(i+1): [r for r in profile['attempts'][i] if r['reason'] or r['flags'] or
                        len(profile['attempts'][i]) > 1]
              for i in tried}
    # Retain failures and all manual recovery decisions, without every clean
    # optimizer trace or per-fraction spectral array.
    return dict(bin_id=profile['bin_id'], n_sampled=len(tried),
                n_attempts=sum(map(len, profile['attempts'])), fallback=profile['fallback'],
                elapsed_seconds=profile['elapsed_seconds'], flags=profile_flags(profile),
                missing_fractions=FRACTIONS[absent].tolist(),
                moved_on=FRACTIONS[profile['moved_on']].tolist(),
                attempts={i: records for i, records in issues.items() if records})


def bl_population(winner, setup, scale):
    """Compress the winning BL fit to population values and restored SSP weights."""
    basis = setup['population_metadata']
    weights = np.asarray(winner['weights']).reshape(2, setup['n_stellar'])
    result, masses = {}, []
    for label, luminosity in zip('AB', weights):
        stats, grids = populations.summarize_ssp_weights(luminosity, setup['stellar_normalization'], basis)
        masses.append(grids['total_surviving_mass'])
        result.update({f'{label}_{key}': stats[key] for key in BL_STATS})
    result.update(A_STELLAR_MASS_FRACTION=masses[0]/sum(masses), FA_BL=winner['fraction'], FA_BL_ACHIEVED=winner['achieved'],
                  CHI2_BL=winner['chi2'], GAS_V_KMS=winner['gas'][0], GAS_SIGMA_KMS=winner['gas'][1])
    return result, weights*scale


def profile_flags(profile):
    winner = profile['winner']
    flags = list(winner['flags']) if winner is not None else ['NO_USABLE_FIT']
    if winner is not None and winner['index'] in (0, len(FRACTIONS)-1):
        flags.append('FRACTION_EDGE')
    if any(history and score is None for history, score in zip(profile['attempts'], profile['scores'])):
        flags.append('INCOMPLETE_PROFILE')
    return flags
