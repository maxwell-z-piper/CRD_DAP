"""Optional nominal BL sensitivity calculation with a matched XSL Kroupa library."""
from __future__ import annotations

from dataclasses import replace
from pathlib import Path
import shutil

from astropy.table import Table
import numpy as np

from . import bl_population_io
from .config import load_config
from .global_disk_io import file_hash, load_json, save_json
from .script07_fits import BL_STATS, audit_profile, bl_population
from .script07_products import _table, aperture_statistics


def run(directory, nominal, library, workers, unattended, driver, fraction_stage, bin_pool):
    library = Path(library).expanduser().resolve()
    digest = file_hash(library)
    output = directory/'products'
    state = directory/'metadata/kroupa_sensitivity.json'
    if state.exists() and load_json(state)['template_sha256'] != digest:
        raise ValueError('The Kroupa library changed. Use a separate Script 07 output directory.')
    if state.exists() and load_json(state).get('complete'):
        scratch = directory/'work'/'kroupa'
        if scratch.exists():
            shutil.rmtree(scratch)
        return
    manifest = load_json(nominal.paths['manifest6'])
    if digest == manifest['contract']['hashes']['xsl']:
        raise ValueError('Supply a Kroupa library that differs from the nominal SSP library')
    cfg = load_config(nominal.paths['run6']/'metadata/config_snapshot.py', validate=False)
    inputs = bl_population_io.read_inputs(cfg, nominal.paths['run5'],
                                          nominal.paths['run1'], nominal.paths['run2'])
    inputs['contract']['settings'] = nominal.bl_settings
    inputs['paths']['xsl'] = library
    experiment = bl_population_io.build_experiment(inputs, nominal.provenance['bl_review'])
    spectrum = nominal.spectra['BL']
    if not (np.array_equal(experiment['good'], spectrum['good']) and
            np.allclose(experiment['wavelength'], spectrum['wavelength'], atol=1e-8, rtol=0)):
        raise ValueError('Kroupa comparison must use the same approved BL pixels and wavelength sampling')
    setup, gas_masks = bl_population_io.build_templates(inputs, experiment)
    saved = nominal.setup['BL']
    if (setup['gas_names'] != saved['gas_names'] or not np.array_equal(gas_masks, saved['gas_by_bin'])
            or any(not np.array_equal(setup['population_metadata'][k], saved['population_metadata'][k])
                   for k in ('ages', 'metals'))):
        raise ValueError('Kroupa and nominal BL runs must have matched age/metallicity axes and gas masks')
    setup.update(gas_by_bin=gas_masks, full_scan=True)
    alternate = replace(nominal, setup=nominal.setup | {'BL': setup},
                        spectra=nominal.spectra | {'BL': spectrum | {'model': spectrum['galaxy']}})
    save_json(state, dict(template_sha256=digest, library=library, complete=False))
    scratch = directory/'work'/'kroupa'
    scratch.mkdir(exist_ok=True)
    ids = nominal.ids['BL']
    columns = bl_population_io.STELLAR_COLUMNS
    kinematics = {int(bid): [float(nominal.tables['BL'][column][bid]) for column in columns] for bid in ids}
    with bin_pool(workers, alternate) as pool:
        results, missing = fraction_stage(scratch, 'BL', ids, kinematics, alternate,
                                          pool, workers, 0, 0, unattended, driver)
    data = {'BIN_ID': ids}
    weights = np.full((len(ids), 2, setup['n_stellar']), np.nan)
    original_weights = np.full_like(weights, np.nan)
    measures = {}
    for i, bid in enumerate(ids):
        best = results[bid]['winner']
        if best is not None:
            scale = spectrum['normalization_scale'][bid]
            measures[bid], weights[i] = bl_population(best, setup, scale)
            original_weights[i] = spectrum['model_weights'][int(bid)][:2*setup['n_stellar']].reshape(2, -1)*scale
    for name in ('FA_BL', *(f'{disk}_{key}' for disk in 'AB' for key in BL_STATS)):
        values = np.array([measures.get(b, {}).get(name, np.nan) for b in ids])
        data[f'{name}_KROUPA'] = values
        data[f'{name}_KROUPA_MINUS_NOMINAL'] = values - np.asarray(nominal.tables['BL'][name][ids], dtype=float)
    _table(output/'Kroupa_BL_sensitivity.ecsv', Table(data))
    kroupa, kroupa_weights, count = aperture_statistics(weights, setup)
    baseline, salpeter_weights, _ = aperture_statistics(original_weights, saved)
    formed_fraction = {}
    for label, grid, basis in (('kroupa', kroupa_weights, setup), ('nominal', salpeter_weights, saved)):
        formed = grid / basis['stellar_normalization'].reshape(grid.shape[1:])
        age_mass = formed.sum(axis=2)
        formed_fraction[label] = (age_mass/age_mass.sum(axis=1, keepdims=True)).tolist()
    save_json(output/'Kroupa_aperture_sensitivity.json',
              dict(usable_bins=count, missing_bins=missing, nominal_on_common_bins=baseline,
                   kroupa=kroupa, kroupa_minus_nominal={k: value-baseline[k] for k, value in kroupa.items()},
                   age_gyr=setup['population_metadata']['ages'],
                   formed_mass_fraction_by_age=formed_fraction,
                   audit=[audit_profile(results[int(b)]) for b in ids]))
    save_json(state, dict(template_sha256=digest, library=library, complete=True))
    shutil.rmtree(scratch)
