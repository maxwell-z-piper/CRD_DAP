"""Fixed stellar kinematics, BL fraction profiles and XSL population summaries."""
from __future__ import annotations

from importlib.resources import files
from pathlib import Path
import time
import warnings

import numpy as np

from . import rh3_refit
from .ppxf_utils import _validate_spectrum_inputs, total_chi2_diagonal


def load_ssp_metadata(path, medium, band):
    """Return age/metal axes and luminosity and surviving mass per unit formed mass."""
    from ppxf import ppxf_util

    with np.load(path, allow_pickle=False) as data:
        ages, metals, masses = (np.asarray(data[key], dtype=float)
                                for key in ('ages', 'metals', 'masses'))
        spectra = np.asarray(data['templates'], dtype=float)
        wavelength = np.asarray(data['lam'], dtype=float)
    if (ages.ndim != 1 or metals.ndim != 1 or not ages.size or not metals.size
            or not np.all(np.isfinite(ages) & (ages > 0)) or not np.isfinite(metals).all()
            or np.any(np.diff(ages) <= 0) or np.any(np.diff(metals) <= 0)
            or masses.shape != (len(ages), len(metals))
            or spectra.shape != (len(wavelength), *masses.shape)
            or not np.all(np.isfinite(masses) & (masses > 0) & (masses <= 1))):
        raise ValueError('XSL population metadata must describe the complete age–metallicity grid.')
    # Use the unnormalized SSPs in Lsun/Angstrom per unit formed mass. Evaluate
    # the filter and the reference Sun on the same vacuum wavelength convention.
    from .templates import air_to_vacuum
    if medium == 'air':
        wavelength = air_to_vacuum(wavelength)
    solar_file = files('ppxf').joinpath('sps_models/spectra_sun_vega.npz')
    with solar_file.open('rb') as stream, np.load(stream, allow_pickle=False) as sun:
        sun_wave, sun_flux = sun['lam'], sun['flux_sun']
    phot = ppxf_util.synthetic_photometry(wavelength, spectra.reshape(len(wavelength), -1),
                                         [band], redshift=0., quiet=True)
    sun = ppxf_util.synthetic_photometry(sun_wave, sun_flux, [band], redshift=0., quiet=True)
    solar_luminosity_cgs = 3.828e33
    ten_pc_cm = 3.085677581491367e19
    luminosity = (np.asarray(phot.flux).reshape(masses.shape) * solar_luminosity_cgs
                  / (4*np.pi*ten_pc_cm**2 * float(np.asarray(sun.flux).item())))
    if not np.all(np.isfinite(luminosity) & (luminosity > 0)):
        raise ValueError(f'XSL spectra or the solar reference do not cover the {band} filter.')
    return dict(ages=ages, metals=metals, masses=masses, luminosity=luminosity)


def summarize_ssp_weights(weights, normalization, metadata):
    """Summarize one disk with explicit light, formed-mass and surviving-mass bases."""
    ages, metals = metadata['ages'], metadata['metals']
    shape = (len(ages), len(metals))
    light = np.asarray(weights).reshape(shape)
    formed = light / np.asarray(normalization).reshape(shape)
    surviving = formed * metadata['masses']
    luminosity = float(np.sum(formed * metadata['luminosity']))
    stats = {}
    for name, grid in (('LIGHT', light), ('MASS', surviving)):
        unit = grid / grid.sum()
        stats.update({
            f'{name}_AGE_GYR': float(np.sum(unit * ages[:, None])),
            f'{name}_LOGAGE_YR': float(np.sum(unit * (np.log10(ages)[:, None] + 9))),
            f'{name}_METAL_DEX': float(np.sum(unit * metals)),
            f'{name}_OLDEST_SHARE': float(unit[-1].sum()),
            f'{name}_YOUNGEST_SHARE': float(unit[0].sum()),
            f'{name}_METAL_LOW_SHARE': float(unit[:, 0].sum()),
            f'{name}_METAL_HIGH_SHARE': float(unit[:, -1].sum()),
        })
    stats['MASS_TO_LIGHT'] = float(surviving.sum() / luminosity)
    return stats, dict(
        light=light/light.sum(), formed_mass=formed/formed.sum(),
        surviving_mass=surviving/surviving.sum(), sfh_formed_mass=formed.sum(axis=1)/formed.sum(),
        total_surviving_mass=float(surviving.sum()),
    )


def fit_fraction(setup, spectrum, kinematics, fraction, start):
    """Evaluate one prescribed stellar fraction and retain an auditable result."""
    from ppxf.ppxf import ppxf

    galaxy, noise, wavelength, good = _validate_spectrum_inputs(
        spectrum['galaxy'], spectrum['noise'], setup['wavelength'], np.flatnonzero(spectrum['good']))
    gas = spectrum['gas']
    keep = np.r_[np.ones(setup['n_stellar']*2, dtype=bool), gas]
    labels = setup['component'][keep]
    has_gas = bool(gas.any())
    stellar = np.asarray(kinematics).reshape(2, 2)
    starts = stellar.tolist() + ([list(start)] if has_gas else [])
    bounds = [[[v-0.01, v+0.01], [max(0.001, s-0.01), s+0.01]] for v, s in stellar]
    if has_gas:
        bounds.append(setup['gas_bounds'])
    record = dict(start=list(start), status=None, flags=[], reason='', achieved_fraction=None,
                  gas_kinematics=None, warnings=[], chi2=None)
    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always', RuntimeWarning)
            pp = ppxf(
                setup['matrix'][:, keep], galaxy, noise, setup['velscale'], starts,
                component=labels, moments=[-2, -2] + ([2] if has_gas else []),
                bounds=bounds, fraction=float(fraction),
                gas_component=(labels == 2) if has_gas else None,
                gas_names=np.asarray(setup['gas_names'])[gas] if has_gas else None,
                goodpixels=good, lam=wavelength, vsyst=setup['vsyst'],
                degree=-1, mdegree=setup['mdegree'], regul=0, clean=False, quiet=True,
            )
        record['warnings'] = sorted({f'{Path(w.filename).name}:{w.lineno}: {w.message}' for w in caught})
    except MemoryError:
        raise
    except Exception as exc:
        record.update(flags=['PPXF_EXCEPTION'], reason=f'{type(exc).__name__}: {exc}')
        return record, None
    record['status'] = int(pp.status)
    if pp.status <= 0:
        record['flags'].append('PPXF_NOT_CONVERGED')
        record['reason'] = f'pPXF optimizer returned status {pp.status}'
    chi2 = total_chi2_diagonal(galaxy, pp.bestfit, noise, good)
    arrays = [pp.weights, pp.bestfit, chi2, *pp.sol]
    polynomial = pp.mpolyweights if pp.mpolyweights is not None else np.empty(0)
    if not all(np.isfinite(values).all() for values in (*arrays, polynomial)):
        record.update(reason='Nonfinite pPXF output')
        record['flags'].append('NONFINITE_FIT')
        return record, None
    record['chi2'] = chi2
    if not np.allclose(pp.sol[:2], stellar, rtol=0, atol=1e-6):
        raise RuntimeError('BL fit changed the fixed Script 05 stellar kinematics.')
    stellar_weight = pp.weights[labels < 2].sum()
    achieved = float(pp.weights[labels == 0].sum() / stellar_weight) if stellar_weight > 0 else None
    record['achieved_fraction'] = achieved
    if achieved is None or abs(achieved - fraction) > setup['fraction_tolerance']:
        record['flags'].append('FIXED_FRACTION_MISMATCH')
        record['reason'] = 'Stellar weights do not satisfy the imposed BL fraction'
        return record, None
    gas_solution = np.full(2, np.nan)
    if has_gas:
        gas_solution = np.asarray(pp.sol[2])
        record['gas_kinematics'] = gas_solution.tolist()
        lo, hi = np.asarray(setup['gas_bounds']).T
        if np.any((gas_solution < lo-1e-6) | (gas_solution > hi+1e-6)):
            record['flags'].append('GAS_OUT_OF_BOUNDS')
            record['reason'] = 'Returned gas kinematics violate the configured bounds'
            return record, None
        if np.any(np.minimum(gas_solution-lo, hi-gas_solution) <= 2):
            record['flags'].append('GAS_KINEMATIC_BOUNDARY')
        if not np.any(pp.weights[labels == 2] > 0):
            record['flags'].append('NO_FITTED_GAS_FLUX')
    else:
        record['flags'].append('NO_USABLE_GAS_LINES')
    a, b, g = (pp.matrix[:, labels == i] @ pp.weights[labels == i] for i in range(3))
    if not np.allclose(a+b+g, pp.bestfit, rtol=1e-8, atol=1e-10):
        raise RuntimeError('BL component spectra do not reconstruct the total model.')
    weights = np.zeros(len(keep))
    weights[keep] = pp.weights
    return record, dict(chi2=chi2, achieved_fraction=achieved, gas_kinematics=gas_solution,
                        model_A=a, model_B=b, model_gas=g, model_total=pp.bestfit,
                        weights=weights, multiplicative_coefficients=polynomial)


def fit_bl_fraction_grid(setup, bin_id, spectrum, kinematics):
    """Run one complete bin from the same initial gas seed at all 99 fractions."""
    started = time.perf_counter()
    fractions = rh3_refit.FRACTIONS
    n, nw = len(fractions), len(setup['wavelength'])
    profile = dict(bin_id=bin_id, velocities=np.asarray(kinematics)[[0, 2]],
                   kinematics=np.asarray(kinematics), fraction=fractions.copy(),
                   attempts=[[] for _ in fractions], pending=np.ones(n, dtype=bool),
                   moved_on=np.zeros(n, dtype=bool), winning_attempt=np.full(n, -1, dtype=int))
    shapes = dict(chi2=(), achieved_fraction=(), gas_kinematics=(2,),
                  weights=(len(setup['component']),), multiplicative_coefficients=(setup['mdegree'],))
    shapes.update({name: (nw,) for name in ('model_A', 'model_B', 'model_gas', 'model_total')})
    profile.update({key: np.full((n, *shape), np.nan) for key, shape in shapes.items()})
    for i, fraction in enumerate(fractions):
        record, result = fit_fraction(setup, spectrum, kinematics, fraction, setup['gas_seed'])
        rh3_refit.record_attempt(profile, i, record, result)
    profile['elapsed_seconds'] = time.perf_counter() - started
    return profile


def initialize_worker(setup):
    global _WORKER_SETUP
    _WORKER_SETUP = setup


def fit_worker_bin(bin_id, spectrum, kinematics):
    return fit_bl_fraction_grid(_WORKER_SETUP, bin_id, spectrum, kinematics)
