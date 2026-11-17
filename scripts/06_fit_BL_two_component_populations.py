#!/usr/bin/env python3
"""Review BL line coverage, then fit populations at fixed Script 05 kinematics."""
from __future__ import annotations

import os
for _name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
              'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_name] = '1'

import argparse
from collections import Counter
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from importlib.metadata import version
import logging
from pathlib import Path
import sys
import time

from astropy.io import fits
from astropy.table import Table
import numpy as np

from crd_utils import bl_population_io as io
from crd_utils import bl_population_plotting as plots
from crd_utils import bl_population_products as products
from crd_utils import populations, ppxf_utils, rh3_refit, scalar_noise, stellar_psf, templates
from crd_utils.config import load_config, snapshot_config, write_config_manifest
from crd_utils.global_disk_io import file_hash, load_json, save_json
from crd_utils.logging_utils import RunContext, log_section, setup_pipeline_logger, setup_step_logger
from crd_utils.runtime import BinProgress, output_lock

SCRIPT = Path(__file__).stem

# The original Script 06 included the plotting source in its fit fingerprint.
LEGACY_CODE = {
    '06_fit_BL_two_component_populations.py': 'e904dd64bbe1ca51fb48f1d5e00ff1958e35c0bb0ea14dfe73df031e5323bcbd',
    'bl_population_plotting.py': 'a4901310e92b7307cc5870c56e887cbae403fbfecf0916cd1de2c7012d11c68d',
}


def arguments():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--phase', required=True, choices=('review', 'fit', 'psf'))
    p.add_argument('--config', type=Path, required=True)
    p.add_argument('--script5-run', type=Path, required=True)
    p.add_argument('--script1-run', type=Path, help='Override the saved BL preparation run path')
    p.add_argument('--script2-run', type=Path, help='Override the saved BL binning run path')
    p.add_argument('--output-dir', type=Path, required=True)
    p.add_argument('--workers', type=int, help='Independent BL bin profiles')
    p.add_argument('--approve-review', action='store_true',
                   help='Approve a previously generated, matching line-review PDF')
    p.add_argument('--resume', action='store_true', help='Reuse matching BL checkpoints in the fit phase')
    p.add_argument('--non-interactive', action='store_true', help='Save unresolved recovery decisions')
    args = p.parse_args()
    if args.approve_review and args.phase != 'review':
        p.error('--approve-review belongs to the review phase.')
    if args.resume and args.phase != 'fit':
        p.error('--resume belongs to the fit phase.')
    return args


def write_table(path, table):
    temporary = path.with_suffix(path.suffix+'.tmp')
    table.write(temporary, format='ascii.ecsv', overwrite=True)
    temporary.replace(path)


def write_review_products(run, inputs, experiment, review):
    rows = []
    ids = np.flatnonzero(inputs['eligible'])
    for index, line in enumerate(io.LINES):
        for bid in ids:
            rows.append(dict(
                BIN_ID=int(bid), LINE=line, REST_WAVELENGTH=experiment['line_wave'][index],
                OBSERVED_WAVELENGTH_AT_SYSTEMIC=experiment['line_obs'][index],
                CORE_LO=experiment['core'][index, 0], CORE_HI=experiment['core'][index, 1],
                CENTER_SEARCH_LO=experiment['search'][index, 0],
                CENTER_SEARCH_HI=experiment['search'][index, 1],
                CORE_FRACTION=experiment['coverage'][bid, index],
                PROPOSED=line in review['lines'] and inputs['contract']['settings']['gas']['enabled'],
                USABLE=bool(experiment['usable'][bid, index]),
            ))
    table = Table(rows=rows)
    table.meta.update(rest_medium=inputs['contract']['settings']['template_medium'],
                      observed_medium=inputs['contract']['settings']['science_medium'],
                      wavelength_unit='Angstrom',
                      core_definition='Nominal FWHM including BL LSF and gas seed dispersion',
                      threshold='At least half the core and at least two good log pixels',
                      partial_definition='Core fraction below one with nonzero coverage',
                      effective_fit_range=list(map(float, experiment['effective_range'])))
    write_table(run.products_dir/'BL_line_coverage.ecsv', table)
    plots.line_review(run.figures_dir/'BL_line_review.pdf', inputs, experiment, review)


def psf_phase(run, cfg, inputs, logger):
    """Measure blue standards for QA without using the result in the BL fitter."""
    paths = [Path(path).expanduser().resolve() for path in cfg.BL_PSF_STANDARD_CUBES]
    if not paths:
        raise ValueError('Set BL_PSF_STANDARD_CUBES to reduced native BL standard-star cubes.')
    settings = inputs['contract']['settings']
    redshift = inputs['contract']['redshift']
    observed = getattr(cfg, 'BL_PSF_OBSERVED_RANGE_ANGSTROM', None)
    measurements = []
    records = []
    for path in paths:
        medium = templates.infer_science_wavelength_medium(fits.getheader(path))
        if observed is None:
            band = templates.convert_wavelength_medium(settings['reference'],
                                                        settings['template_medium'], medium)*(1+redshift)
        else:
            band = templates.convert_wavelength_medium(io.interval(observed, 'BL PSF observed range'),
                                                        settings['science_medium'], medium)
        measurements.append(stellar_psf.measure_star(path, band))
        records.append(dict(path=path, sha256=file_hash(path), observed_band=band, medium=medium))
    widths = {
        kind: float(np.median([max(item['candidates'][kind]['fwhm_intrinsic_axes_arcsec'])
                              for item in measurements]))
        for kind in ('gaussian', 'moffat')
    }
    stellar_psf.plot_stars(measurements, run.figures_dir/'BL_stellar_PSF_QA.png')
    save_json(run.products_dir/'BL_stellar_psf.json',
              dict(inputs=records, target=inputs['contract']['target'], measurements=measurements,
                   median_larger_axis_fwhm_arcsec=widths, use='Standard-star QA only'))
    logger.info('Saved blue PSF QA | median larger-axis FWHM (arcsec): %s', widths)
    return 0


def review_phase(run, inputs, default, base, args, logger):
    proposed = run.metadata_dir/'BL_line_review.json'
    preview = run.metadata_dir/'BL_line_preview.json'
    approved = run.metadata_dir/'BL_line_approval.json'
    if (run.metadata_dir/'script06_resume_state.json').exists():
        raise ValueError('BL fitting has started. Use a new output directory to revise the line review.')
    if not proposed.exists():
        save_json(proposed, default)

    def approve(review, signature):
        files = load_json(preview)['products']
        for path in (run.figures_dir/'BL_line_review.pdf', run.products_dir/'BL_line_coverage.ecsv'):
            if file_hash(path) != files[path.name]:
                raise ValueError(f'Regenerate the line preview before approving: {path}')
        save_json(approved, dict(input_fingerprint=base, review_fingerprint=signature, review=review))
        logger.info('Approved BL lines and masks. Run the fit phase in this directory.')
        return 0

    while True:
        review = load_json(proposed)
        io.validate_review(review, len(inputs['eligible']))
        signature = io.fingerprint(dict(input_fingerprint=base, review=review))
        if args.approve_review:
            if not preview.is_file() or load_json(preview)['fingerprint'] != signature:
                raise ValueError('Preview the current line list and masks before approving them.')
            return approve(review, signature)
        experiment = io.build_experiment(inputs, review)
        write_review_products(run, inputs, experiment, review)
        files = (run.figures_dir/'BL_line_review.pdf', run.products_dir/'BL_line_coverage.ecsv')
        save_json(preview, dict(fingerprint=signature, products={p.name: file_hash(p) for p in files}))
        logger.info('Inspect %s and %s.', *files)
        logger.info('Edit %s to change the proposed lines or add documented wavelength masks.', proposed)
        if args.non_interactive or not sys.stdin.isatty():
            logger.info('Run the review phase again with --approve-review after inspecting the preview.')
            return 2
        print('1. Approve this preview\n2. Reload the edited review JSON and regenerate\n3. Pause')
        try:
            choice = input('Line review [1/2/3] > ').strip()
            while choice not in ('1', '2', '3'):
                choice = input('Choose 1, 2 or 3 > ').strip()
        except EOFError:
            return 2
        if choice == '1':
            if io.fingerprint(dict(input_fingerprint=base, review=load_json(proposed))) != signature:
                logger.info('The proposal changed. Regenerating the preview before approval.')
                continue
            return approve(review, signature)
        if choice == '3':
            return 2


def spectrum_for(experiment, setup, bid):
    return {key: experiment[key][bid] for key in ('galaxy', 'noise', 'good')} | {
        'gas': setup['gas_by_bin'][bid],
    }


def initial_profiles(inputs, experiment, setup, jobs, progress):
    if progress.workers == 1 and not progress.enabled:
        for bid in jobs:
            yield populations.fit_bl_fraction_grid(setup, bid, spectrum_for(experiment, setup, bid),
                                                    inputs['kin'][bid])
        return
    with ProcessPoolExecutor(max_workers=progress.workers, initializer=populations.initialize_worker,
                             initargs=(setup,)) as pool:
        indices = iter(jobs)

        def submit(bid):
            return pool.submit(populations.fit_worker_bin, bid, spectrum_for(experiment, setup, bid),
                               inputs['kin'][bid])

        running = {submit(next(indices)) for _ in range(progress.workers)}
        while running:
            done, running = wait(running, timeout=5. if progress.enabled else None,
                                 return_when=FIRST_COMPLETED)
            if not done:
                progress.refresh()
            for future in done:
                yield future.result()
                bid = next(indices, None)
                if bid is not None:
                    running.add(submit(bid))


def fraction_ranges(indices):
    groups = np.split(indices, np.flatnonzero(np.diff(indices) != 1)+1)
    f = rh3_refit.FRACTIONS
    return ', '.join(f'{f[g[0]]:.2f}–{f[g[-1]]:.2f}' if len(g) > 1 else f'{f[g[0]]:.2f}'
                     for g in groups)


def retry_choice(profile, pending, setup, has_gas):
    bid = int(profile['bin_id'])
    reasons = Counter(profile['attempts'][i][-1]['reason'] for i in pending)
    print(f'\nBin {bid} | {len(pending)} fractions pending | f_A {fraction_ranges(pending)}')
    for reason, count in reasons.items():
        print(f'  {count} × {reason}')
    print('1. Retry the pending group with its latest starts\n'
          '2. Enter one new gas start for the pending group\n'
          '3. Move on and retain the completed fits')
    while True:
        option = input('BL recovery [1/2/3] > ').strip()
        if option == '1':
            return [profile['attempts'][i][-1]['start'] for i in pending]
        if option == '3':
            return None
        if option == '2' and has_gas:
            try:
                start = list(map(float, input('Enter gas V and dispersion in km/s > ').split()))
            except ValueError:
                start = []
            if len(start) == 2 and all(np.isfinite(v) and lo <= v <= hi
                                       for v, (lo, hi) in zip(start, setup['gas_bounds'])):
                return [start]*len(pending)
            print(f'Enter finite starts within {setup["gas_bounds"]}.')
        elif option == '2':
            print('This bin has no usable gas templates. Choose 1 or 3.')
        else:
            print('Choose 1, 2 or 3.')


def recover(profile, path, signature, setup, spectrum, interactive, logger):
    while interactive and profile['pending'].any():
        pending = np.flatnonzero(profile['pending'])
        try:
            starts = retry_choice(profile, pending, setup, bool(spectrum['gas'].any()))
        except EOFError:
            return False
        bid = int(profile['bin_id'])
        if starts is None:
            profile['pending'][pending] = False
            profile['moved_on'][pending] = True
            rh3_refit.save_profile(path, profile, signature)
            logger.info('Bin %d moved on from %d fractions.', bid, len(pending))
            break
        logger.info('Bin %d retrying %d fractions.', bid, len(pending))
        for index, start in zip(pending, starts):
            reserved = dict(start=list(start), status=None, reason='Recovery attempt interrupted',
                            flags=['INTERRUPTED_ATTEMPT'])
            rh3_refit.record_attempt(profile, index, reserved, None)
            rh3_refit.save_profile(path, profile, signature)
            record, fit = populations.fit_fraction(setup, spectrum, profile['kinematics'],
                                                    profile['fraction'][index], start)
            profile['attempts'][index].pop()
            rh3_refit.record_attempt(profile, index, record, fit)
            rh3_refit.save_profile(path, profile, signature)
            logger.info('Bin %d f_A=%.2f recovery | status %s | %s', bid, profile['fraction'][index],
                        record['status'], record['reason'] or 'fit returned', extra={'file_only': True})
    return interactive


def fit_profiles(run, inputs, experiment, setup, signature, args, workers, logger):
    eligible = inputs['eligible'] & np.isfinite(experiment['normalization'])
    eligible &= experiment['good'].sum(axis=1) >= max(5, setup['mdegree']+3)
    ids = np.flatnonzero(eligible)
    if not len(ids):
        raise ValueError('No eligible bin has enough valid BL pixels and a nonzero normalization.')
    profile_dir = run.products_dir/'BL_fraction_profiles'
    profile_dir.mkdir(exist_ok=True)
    paths = {int(bid): profile_dir/f'bin_{bid:04d}.npz' for bid in ids}
    jobs = [bid for bid, path in paths.items() if not path.is_file()]
    interactive = sys.stdin.isatty() and not args.non_interactive
    rows = {}
    pending = 0

    def finish(profile, path):
        nonlocal interactive, pending
        bid = int(profile['bin_id'])
        spectrum = spectrum_for(experiment, setup, bid)
        interactive = recover(profile, path, signature, setup, spectrum, interactive, logger)
        if profile['pending'].any():
            pending += int(profile['pending'].sum())
            return
        if 'chi2_calibrated' not in profile:
            products.finalize_profile(profile, setup, spectrum)
            rh3_refit.save_profile(path, profile, signature)
        figure = run.figures_dir/'BL_fraction_bins'/f'bin_{bid:04d}.pdf'
        if not figure.exists():
            plots.bin_diagnostic(figure, bid, profile, setup['wavelength'], spectrum['galaxy'],
                                 spectrum['noise'], spectrum['good'], inputs['contract']['settings'])
        row = products.summary_row(inputs, experiment, setup, bid, profile, path.relative_to(run.run_dir))
        rows[bid] = row
        logger.info('Bin %d finalized | %d attempts | %d missing fractions | %s', bid,
                    row['N_ATTEMPTS'], row['N_MISSING_FRACTIONS'], row['FLAGS'] or 'no BL flags',
                    extra={'file_only': True})

    logger.info('%d eligible BL bins | %d saved profiles | %d new profiles | %d workers',
                len(ids), len(ids)-len(jobs), len(jobs), workers)
    for bid, path in paths.items():
        if path.is_file():
            kin = inputs['kin'][bid]
            profile = rh3_refit.load_profile(path, signature, bid, kin[[0, 2]])
            if not np.array_equal(profile['kinematics'], kin):
                raise ValueError(f'Script 06 checkpoint kinematics differ for bin {bid}.')
            finish(profile, path)
    if jobs:
        progress = BinProgress(len(ids), len(ids)-len(jobs), min(workers, len(jobs)))
        progress.refresh()
        try:
            for profile in initial_profiles(inputs, experiment, setup, jobs, progress):
                progress.clear()
                bid = int(profile['bin_id'])
                rh3_refit.save_profile(paths[bid], profile, signature)
                progress.checkpoint(profile['elapsed_seconds'])
                logger.info('Bin %d initial scan saved | %d/%d new scans | %.2f s/fit', bid,
                            progress.done, len(jobs), profile['elapsed_seconds']/99,
                            extra={'file_only': True})
                finish(profile, paths[bid])
                if progress.done < len(jobs):
                    progress.refresh()
        finally:
            progress.clear()
    if pending:
        logger.warning('%d fraction decisions remain. Resume the fit in a terminal.', pending)
        return 2, None
    table = Table(rows=[rows[bid] if bid in rows else products.summary_row(inputs, experiment, setup, bid)
                        for bid in range(len(inputs['eligible']))])
    return 0 if np.any(table['USABLE']) else 1, table


def export(run, inputs, experiment, setup, signature, table):
    settings = inputs['contract']['settings']
    table.meta.update(fingerprint=signature, script05=inputs['contract']['script05'],
                      target=inputs['contract']['target'], model=inputs['contract']['model'],
                      fraction_definition='Stellar template-weight share in the fixed BL reference band',
                      fraction_reference_angstrom=list(settings['reference']),
                      template_medium=settings['template_medium'], imf=settings['imf'],
                      mass_to_light_filter=settings['ml_band'], velocity_frame='Script 05 LOS frame',
                      chi2_definition='CHI2_TOTAL uses the BL residual scalar exactly once')
    for col in (*io.STELLAR_COLUMNS, 'GAS_V_KMS', 'GAS_SIGMA_KMS'):
        table[col].unit = 'km / s'
    table_path = run.products_dir/'BL_stellar_populations.ecsv'
    write_table(table_path, table)
    q = np.asarray(table['RESIDUAL_NOISE_SCALE'], dtype=float)
    products.write_npz(run.products_dir/'BL_log_spectra.npz',
                       fingerprint=signature, wavelength=setup['wavelength'],
                       galaxy=experiment['galaxy'], good=experiment['good'],
                       gas_template_available=setup['gas_by_bin'],
                       normalization_scale=experiment['normalization'],
                       noise_before_residual_scale=experiment['noise'], residual_noise_scale=q,
                       noise=experiment['noise']*q[:, None], source_flux_unit=inputs['source_flux_unit'],
                       applied_aperture_scale=inputs['aperture_scale'])
    products.write_npz(run.products_dir/'XSL_BL_templates.npz',
                       fingerprint=signature, templates=setup['matrix'],
                       wavelength=setup['template_wavelength'], component=setup['component'],
                       stellar_normalization=setup['stellar_normalization'],
                       gas_names=np.asarray(setup['gas_names'], dtype=str),
                       velscale=setup['velscale'], vsyst=setup['vsyst'],
                       template_medium=settings['template_medium'], imf=settings['imf'],
                       mass_to_light_filter=settings['ml_band'], **setup['population_metadata'])
    plots.population_maps(run.figures_dir/'BL_stellar_populations.pdf', table, inputs, settings['ml_band'])
    save_json(run.metadata_dir/'script06_manifest.json', dict(
        script=SCRIPT, fingerprint=signature, contract=inputs['contract'],
        input_paths=inputs['paths'], n_rh3_eligible=int(sum(table['RH3_ELIGIBLE'])),
        n_usable_bins=int(sum(table['USABLE'])), n_calibrated_bins=int(sum(table['NOISE_CALIBRATED'])),
        n_incomplete_profiles=int(sum(table['N_MISSING_FRACTIONS'] > 0)),
        fit_range_rest_angstrom=experiment['effective_range'],
        products=dict(population_table='products/BL_stellar_populations.ecsv',
                      fraction_profiles='products/BL_fraction_profiles',
                      spectra='products/BL_log_spectra.npz', templates='products/XSL_BL_templates.npz',
                      line_coverage='products/BL_line_coverage.ecsv',
                      line_review='figures/BL_line_review.pdf', approval='metadata/BL_line_approval.json'),
    ))


def main():
    args = arguments()
    cfg = load_config(args.config, validate=False)
    workers = args.workers if args.workers is not None else getattr(cfg, 'N_WORKERS', None)
    workers = 1 if workers is None else int(workers)
    if workers < 1:
        raise ValueError('--workers must be at least 1.')
    inputs = io.read_inputs(cfg, args.script5_run, args.script1_run, args.script2_run)
    default = io.default_review(cfg, inputs)
    contract = inputs['contract']
    contract['review_defaults'] = default
    contract['versions'] = {name: version(name) for name in ('numpy', 'scipy', 'ppxf', 'capfit', 'astropy')}
    contract['code'] = {Path(path).name: file_hash(path) for path in (
        __file__, io.__file__, products.__file__, populations.__file__,
        templates.__file__, ppxf_utils.__file__, rh3_refit.__file__, scalar_noise.__file__,
    )}
    base = io.fingerprint(contract)
    directory = args.output_dir.expanduser().resolve()
    directory.mkdir(parents=True, exist_ok=True)
    with output_lock(directory, '.script06.lock'):
        if any((directory/'metadata').glob('script0[1-5]*.json')):
            raise ValueError('Choose a Script 06 output directory separate from the upstream runs.')
        run = RunContext(directory, *(directory/name for name in ('logs', 'products', 'figures',
                                                                 'diagnostics', 'metadata')))
        for path in (run.logs_dir, run.products_dir, run.figures_dir, run.diagnostics_dir, run.metadata_dir):
            path.mkdir(exist_ok=True)
        master = setup_pipeline_logger(run)
        if sys.stderr.isatty():
            for handler in master.handlers:
                if not isinstance(handler, logging.FileHandler):
                    handler.addFilter(lambda record: not getattr(record, 'file_only', False))
        logger = setup_step_logger(run, SCRIPT)
        log_section(logger, 'Script 06 BL stellar populations')
        if args.phase == 'psf':
            return psf_phase(run, cfg, inputs, logger)
        if args.phase == 'review':
            return review_phase(run, inputs, default, base, args, logger)
        approved = load_json(run.metadata_dir/'BL_line_approval.json')
        if (approved['input_fingerprint'] != base
                and contract['code'][Path(__file__).name] == plots.LEGACY_COMPATIBLE_DRIVER_SHA256):
            # Keep the original provenance and checkpoint IDs for the prior
            # release only if the entire reconstructed fit contract matches.
            legacy = contract | {'code': contract['code'] | LEGACY_CODE}
            legacy_base = io.fingerprint(legacy)
            if approved['input_fingerprint'] == legacy_base:
                inputs['contract'] = contract = legacy
                base = legacy_base
        review = load_json(run.metadata_dir/'BL_line_review.json')
        signature = io.fingerprint(dict(input_fingerprint=base, review=review))
        if (approved['input_fingerprint'] != base or approved['review'] != review
                or approved['review_fingerprint'] != signature):
            raise ValueError('The BL inputs or reviewed masks changed. Regenerate and approve the line review.')
        state = run.metadata_dir/'script06_resume_state.json'
        if state.exists():
            if not args.resume or load_json(state)['fingerprint'] != signature:
                raise ValueError('Use --resume with matching inputs or start in a new output directory.')
        elif args.resume:
            raise ValueError('No Script 06 fit checkpoint exists to resume.')
        experiment = io.build_experiment(inputs, review)
        setup, gas_by_bin = io.build_templates(inputs, experiment)
        setup['gas_by_bin'] = gas_by_bin
        if not state.exists():
            snapshot_config(cfg, run.metadata_dir)
            write_config_manifest(cfg, run.metadata_dir/'config_manifest.json')
            save_json(state, dict(fingerprint=signature, contract=contract, approved_review=approved))
        logger.info('Fit range %.1f–%.1f Å | %d eligible RH3 bins | %d nominal fraction trials',
                    *experiment['effective_range'], inputs['eligible'].sum(), 99*inputs['eligible'].sum())
        size = 99*inputs['eligible'].sum()*(4*len(setup['wavelength'])+len(setup['component'])
                                           +setup['mdegree'])*8
        logger.info('Approximately %.2f GB of uncompressed fraction spectra and weights.', size/1e9)
        logger.info('Gas seed %s km/s | bounds %s | %s', setup['gas_seed'], setup['gas_bounds'],
                    contract['settings']['gas']['source'])
        try:
            status, table = fit_profiles(run, inputs, experiment, setup, signature, args, workers, logger)
            if table is not None:
                export(run, inputs, experiment, setup, signature, table)
                logger.info('Saved Script 06 handoff with %d usable BL bins.', sum(table['USABLE']))
            return status
        except (Exception, KeyboardInterrupt):
            logger.exception('Script 06 interrupted. Saved profiles remain available with --resume.')
            raise


if __name__ == '__main__':
    sys.exit(main())
