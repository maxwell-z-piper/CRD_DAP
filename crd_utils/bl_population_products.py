"""Finalize BL noise calibration, summary tables and population products."""
from __future__ import annotations

import numpy as np

from . import populations, rh3_refit, scalar_noise


def write_npz(path, **arrays):
    """Publish an NPZ only after its complete contents have been written."""
    temporary = path.with_suffix(path.suffix+'.tmp')
    with temporary.open('wb') as stream:
        np.savez_compressed(stream, **arrays)
    temporary.replace(path)


def finalize_profile(profile, setup, spectrum):
    """Calibrate a bin once from its best unregularized fit."""
    best = rh3_refit.best_fraction(profile)
    profile.update(residual_scale=np.nan, residual_center=np.nan,
                   chi2_calibrated=np.full_like(profile['chi2'], np.nan),
                   standardized_residual=np.full_like(spectrum['galaxy'], np.nan))
    if best is None:
        return
    try:
        scale, center, residual = scalar_noise.from_best_grid_state(
            spectrum['galaxy'], profile['model_total'][best], spectrum['noise'],
            np.flatnonzero(spectrum['good']))
    except ValueError:
        # An undefined scale is explicit. The uncalibrated fit still has a
        # valid score and remains available to inspect.
        pass
    else:
        profile.update(residual_scale=scale, residual_center=center,
                       chi2_calibrated=scalar_noise.rescale_chi2(profile['chi2'], scale),
                       standardized_residual=residual)

    # Report all separated local minima as a descriptive shape diagnostic.
    # Interior points adjacent to failed fits do not establish a local minimum.
    score = profile['chi2']
    finite = np.isfinite(score)
    inner = (finite[:-2] & finite[1:-1] & finite[2:]
             & (score[1:-1] <= score[:-2]) & (score[1:-1] <= score[2:])
             & ((score[1:-1] < score[:-2]) | (score[1:-1] < score[2:])))
    endpoints = (finite[0] and finite[1] and score[0] < score[1],
                 finite[-1] and finite[-2] and score[-1] < score[-2])
    profile['local_minimum_indices'] = np.flatnonzero(np.r_[endpoints[0], inner, endpoints[1]])

    n = setup['n_stellar']
    pop = setup['population_metadata']
    norm = setup['stellar_normalization']
    stats_a, grids_a = populations.summarize_ssp_weights(profile['weights'][best, :n], norm, pop)
    stats_b, grids_b = populations.summarize_ssp_weights(profile['weights'][best, n:2*n], norm, pop)
    for label, grids, stats in (('A', grids_a, stats_a), ('B', grids_b, stats_b)):
        profile.update({f'{label}_{key}': value for key, value in grids.items()})
        profile.update({f'{label}_{key}': value for key, value in stats.items()})
    mass_a, mass_b = profile['A_total_surviving_mass'], profile['B_total_surviving_mass']
    profile['stellar_mass_fraction_A'] = mass_a/(mass_a+mass_b)


SUMMARY_STATS = ('LIGHT_AGE_GYR', 'LIGHT_LOGAGE_YR', 'LIGHT_METAL_DEX', 'LIGHT_OLDEST_SHARE',
                 'LIGHT_YOUNGEST_SHARE', 'LIGHT_METAL_LOW_SHARE', 'LIGHT_METAL_HIGH_SHARE',
                 'MASS_AGE_GYR', 'MASS_LOGAGE_YR', 'MASS_METAL_DEX', 'MASS_OLDEST_SHARE',
                 'MASS_YOUNGEST_SHARE', 'MASS_METAL_LOW_SHARE', 'MASS_METAL_HIGH_SHARE',
                 'MASS_TO_LIGHT')
CONTRASTS = tuple(key for key in SUMMARY_STATS if not key.endswith('_SHARE'))


def summary_row(inputs, experiment, setup, bid, profile=None, profile_path=''):
    upstream = inputs['table'][bid]
    flags = []
    row = dict(
        BIN_ID=bid, SELECTED=bool(upstream['SELECTED']),
        RH3_ELIGIBLE=bool(inputs['eligible'][bid]), USABLE=False, NOISE_CALIBRATED=False,
        VA_GLOBAL_KMS=float(inputs['kin'][bid, 0]), SIGMA_A_KMS=float(inputs['kin'][bid, 1]),
        VB_GLOBAL_KMS=float(inputs['kin'][bid, 2]), SIGMA_B_KMS=float(inputs['kin'][bid, 3]),
        RH3_FLAGS='' if np.ma.is_masked(upstream['FLAGS']) else str(upstream['FLAGS']),
        FA_RH3=float(upstream['FA_RH3']) if not np.ma.is_masked(upstream['FA_RH3']) else np.nan,
        FA_BL=np.nan, FA_BL_ACHIEVED=np.nan, A_STELLAR_MASS_FRACTION=np.nan,
        GAS_V_KMS=np.nan, GAS_SIGMA_KMS=np.nan, GAS_TEMPLATES='', PPXF_STATUS=-999,
        CHI2_BEFORE_RESIDUAL_SCALE=np.nan, CHI2_TOTAL=np.nan,
        RESIDUAL_NOISE_SCALE=np.nan, N_GOOD_PIXELS=int(experiment['good'][bid].sum()),
        N_MISSING_FRACTIONS=0, N_ATTEMPTS=0, N_LOCAL_MINIMA=0,
        PROFILE=str(profile_path),
    )
    row.update({f'{label}_{key}': np.nan for label in ('A', 'B') for key in SUMMARY_STATS})
    row.update({f'DELTA_{key}': np.nan for key in CONTRASTS})
    if profile is None:
        if not row['SELECTED']:
            flags.append('OUTSIDE_SELECTED_MODEL')
        elif not row['RH3_ELIGIBLE']:
            flags.append('NO_USABLE_RH3_KINEMATICS')
        else:
            flags.append('INSUFFICIENT_BL_PIXELS_OR_FLUX')
    else:
        row['N_ATTEMPTS'] = sum(map(len, profile['attempts']))
        row['N_MISSING_FRACTIONS'] = int(np.count_nonzero(~np.isfinite(profile['chi2'])))
        if row['N_MISSING_FRACTIONS']:
            flags.append('INCOMPLETE_PROFILE')
        best = rh3_refit.best_fraction(profile)
        if best is None:
            flags.append('NO_USABLE_BL_FIT')
        else:
            winner = profile['attempts'][best][profile['winning_attempt'][best]]
            flags.extend(winner['flags'])
            if best in (0, len(profile['fraction'])-1):
                flags.append('FRACTION_EDGE')
            if len(profile['local_minimum_indices']) > 1:
                flags.append('MULTIPLE_PROFILE_MINIMA')
            if not np.isfinite(profile['residual_scale']):
                flags.append('RESIDUAL_NOISE_SCALE_UNDEFINED')
            v, sigma = profile['gas_kinematics'][best]
            fitted_groups = setup['gas_by_bin'][bid]
            row.update(USABLE=True, FA_BL=profile['fraction'][best],
                       FA_BL_ACHIEVED=profile['achieved_fraction'][best],
                       A_STELLAR_MASS_FRACTION=profile['stellar_mass_fraction_A'],
                       GAS_V_KMS=v, GAS_SIGMA_KMS=sigma,
                       GAS_TEMPLATES='|'.join(np.asarray(setup['gas_names'], dtype=str)[fitted_groups]),
                       PPXF_STATUS=winner['status'],
                       CHI2_BEFORE_RESIDUAL_SCALE=profile['chi2'][best],
                       CHI2_TOTAL=profile['chi2_calibrated'][best],
                       RESIDUAL_NOISE_SCALE=profile['residual_scale'],
                       NOISE_CALIBRATED=bool(np.isfinite(profile['residual_scale'])),
                       N_LOCAL_MINIMA=len(profile['local_minimum_indices']))
            for label in ('A', 'B'):
                for key in SUMMARY_STATS:
                    row[f'{label}_{key}'] = float(profile[f'{label}_{key}'])
                # Flags denote weight on an edge, with only a numerical
                # tolerance. The continuous shares allow later science cuts.
                if profile[f'{label}_LIGHT_OLDEST_SHARE'] > 1e-6:
                    flags.append(f'{label}_OLDEST_AGE_WEIGHT')
                if profile[f'{label}_LIGHT_YOUNGEST_SHARE'] > 1e-6:
                    flags.append(f'{label}_YOUNGEST_AGE_WEIGHT')
                if max(profile[f'{label}_LIGHT_METAL_LOW_SHARE'],
                       profile[f'{label}_LIGHT_METAL_HIGH_SHARE']) > 1e-6:
                    flags.append(f'{label}_METALLICITY_EDGE_WEIGHT')
            for key in CONTRASTS:
                row[f'DELTA_{key}'] = row[f'A_{key}']-row[f'B_{key}']
    row['FLAGS'] = '|'.join(flags)
    return row
