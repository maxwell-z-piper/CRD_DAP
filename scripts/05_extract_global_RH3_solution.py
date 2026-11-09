#!/usr/bin/env python3
"""Refit all 99 RH3 fractions at the selected coherent Script 04 velocities."""
from __future__ import annotations

import os
for _name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
              'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_name] = '1'

import argparse
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from contextlib import contextmanager
import hashlib
from importlib.metadata import version
import json
import logging
from pathlib import Path
from shutil import get_terminal_size
import sys
import time

import numpy as np
from astropy.table import Table

from crd_utils import ppxf_utils, rh3_refit as refit, rh3_refit_io, rh3_refit_plotting as plots
from crd_utils.config import load_config, snapshot_config, write_config_manifest
from crd_utils.global_disk_io import file_hash, load_json, save_json
from crd_utils.logging_utils import RunContext, log_section, setup_pipeline_logger, setup_step_logger

SCRIPT = Path(__file__).stem
PROGRESS_HEARTBEAT_SECONDS = 5.0


def format_duration(seconds):
    minutes, seconds = divmod(max(0, round(seconds)), 60)
    hours, minutes = divmod(minutes, 60)
    return f'{hours}:{minutes:02d}:{seconds:02d}' if hours else f'{minutes:02d}:{seconds:02d}'


class BinProgress:
    """Script 03 style terminal heartbeat for checkpointed initial scans."""

    def __init__(self, total, saved, workers):
        self.enabled = sys.stderr.isatty()
        self.total = total
        self.saved = saved
        self.workers = workers
        self.started = time.perf_counter()
        self.last_bin = self.started
        self.done = 0
        self.fit_seconds = 0.
        self.width = 0
        self.spinner_index = 0

    def checkpoint(self, fit_seconds):
        self.done += 1
        self.fit_seconds += fit_seconds
        self.last_bin = time.perf_counter()

    def eta(self):
        remaining = self.total - self.saved - self.done
        if not self.done:
            return None
        return self.fit_seconds / self.done * remaining / min(self.workers, max(1, remaining))

    def refresh(self):
        if not self.enabled:
            return
        now = time.perf_counter()
        done = self.saved + self.done
        fraction = done / self.total
        spinner = '|/-\\'[self.spinner_index % 4]
        self.spinner_index += 1
        elapsed = format_duration(now - self.started)
        estimate = self.eta()
        eta = 'ETA warming up' if estimate is None else f'ETA {format_duration(estimate)}'
        last = (f'last new bin {format_duration(now - self.last_bin)} ago'
                if self.done else 'last new bin: none yet')

        def line(bar_width):
            filled = int(fraction * bar_width)
            bar = '#' * filled + '-' * (bar_width - filled)
            return f'{spinner} RUNNING [{bar}] {done}/{self.total} bins ({fraction:5.1%})'

        columns = max(1, get_terminal_size().columns - 1)
        status = f'{line(24)} | workers={self.workers} | elapsed {elapsed} | {last} | {eta}'
        if len(status) > columns:
            status = f'{line(12)} | {elapsed} | {eta}'
        status = status[:columns]
        self.width = min(columns, max(self.width, len(status)))
        sys.stderr.write('\r' + status.ljust(self.width))
        sys.stderr.flush()

    def clear(self):
        if self.width:
            sys.stderr.write('\r' + ' ' * self.width + '\r')
            sys.stderr.flush()
            self.width = 0


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--script4-run', type=Path, required=True, help='Run with an explicitly selected model')
    for n in (1, 2, 3):
        parser.add_argument(f'--script{n}-run', type=Path, help='Override a relocated upstream run')
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--workers', type=int, help='Independent initial bin profiles')
    parser.add_argument('--resume', action='store_true', help='Reuse matching completed and pending bin profiles')
    parser.add_argument('--non-interactive', action='store_true', help='Save retry decisions for interactive resume')
    return parser.parse_args()


@contextmanager
def output_lock(directory):
    """Release the run lock even if the process is interrupted or killed."""
    with (directory / '.script05.lock').open('a+b') as stream:
        stream.seek(0)
        try:
            if os.name == 'nt':
                import msvcrt
                msvcrt.locking(stream.fileno(), msvcrt.LK_NBLCK, 1)
            else:
                import fcntl
                fcntl.flock(stream, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError as exc:
            raise RuntimeError(f'Another Script 05 process holds the run lock: {directory}') from exc
        yield


def spectrum_for(inputs, bin_id):
    return {key: inputs['spectra'][key][bin_id] for key in ('galaxy', 'noise', 'good')}


def initial_profiles(inputs, jobs, workers, progress):
    """Keep only one outstanding initial profile per worker."""
    setup = inputs['setup']
    velocities = inputs['velocities']
    ids = inputs['ids']
    if workers == 1 and not progress.enabled:
        for index in jobs:
            bid = int(ids[index])
            yield refit.fit_initial_profile(setup, bid, spectrum_for(inputs, bid), velocities[index])
        return
    with ProcessPoolExecutor(max_workers=workers, initializer=refit.initialize_worker,
                             initargs=(setup,)) as pool:
        indices = iter(jobs)
        running = set()

        def submit(index):
            bid = int(ids[index])
            running.add(pool.submit(refit.fit_worker_bin, bid,
                                    spectrum_for(inputs, bid), velocities[index]))

        for index in indices:
            submit(index)
            if len(running) == workers:
                break
        timeout = PROGRESS_HEARTBEAT_SECONDS if progress.enabled else None
        while running:
            done, running = wait(running, timeout=timeout, return_when=FIRST_COMPLETED)
            if not done:
                progress.refresh()
            for future in done:
                yield future.result()
                index = next(indices, None)
                if index is not None:
                    submit(index)


def choose_retry(profile, index, setup):
    attempts = profile['attempts'][index]
    last = attempts[-1]
    low, high = setup['sigma_bounds']
    suggestions = [round(setup['baseline']*scale) for scale in (.85, 1.15)]
    suggestions = ', '.join(str(s) for s in suggestions if low <= s <= high)
    print(f'\nBin {int(profile["bin_id"])} | f_A={profile["fraction"][index]:.2f} | '
          f'{len(attempts)-1} of 3 recovery attempts used\n'
          f'Last start {last["start"]} km/s | {last["reason"]}')
    print('1. Retry the same starting pair\n2. Enter new A/B starts\n3. Move on')
    while True:
        choice = input('Choice [1/2/3] > ').strip()
        if choice == '1':
            return last['start']
        if choice == '3':
            return None
        if choice == '2':
            if suggestions:
                print(f'Options near the original MaNGA baseline (km/s) — {suggestions}')
            values = input(f'Enter A and B in km/s, separated by space [{low:g} to {high:g}] > ')
            try:
                start = list(map(float, values.split()))
            except ValueError:
                start = []
            if len(start) == 2 and all(np.isfinite(v) and low <= v <= high for v in start):
                return start
            print('Enter two finite dispersions inside the fitted bounds.')
        else:
            print('Choose 1, 2 or 3.')


def recover(profile, path, fingerprint, inputs, interactive, logger):
    if not interactive:
        return False
    bid = int(profile['bin_id'])
    setup = inputs['setup']
    for index in np.flatnonzero(profile['pending']):
        while profile['pending'][index]:
            try:
                start = choose_retry(profile, index, setup)
            except EOFError:
                logger.warning('Interactive input ended. Pending retries remain saved.')
                return False
            if start is None:
                profile['pending'][index] = False
                profile['moved_on'][index] = True
                refit.save_profile(path, profile, fingerprint)
                logger.info('Bin %d f_A=%.2f | user moved on', bid, profile['fraction'][index])
                break

            # Reserve the attempt before running it so an interrupted call
            # cannot reset the retry budget or lose the chosen starting pair.
            interrupted = dict(start=start, returned=False, status=None, velocity=None, sigma=None,
                               chi2=None, achieved_fraction=None, flags=['INTERRUPTED_ATTEMPT'],
                               reason='Recovery attempt interrupted before a result was saved')
            refit.record_attempt(profile, index, interrupted, None)
            refit.save_profile(path, profile, fingerprint)
            started = time.perf_counter()
            record, products = refit.fit_fraction(setup, spectrum_for(inputs, bid),
                                                profile['velocities'], profile['fraction'][index], start)
            profile['attempts'][index].pop()
            refit.record_attempt(profile, index, record, products)
            profile['elapsed_seconds'] += time.perf_counter() - started
            refit.save_profile(path, profile, fingerprint)
            logger.info('Bin %d f_A=%.2f | start %s km/s | status %s | chi-square %s | %s',
                        bid, profile['fraction'][index], start, record['status'], record['chi2'],
                        ', '.join(record['flags']) or 'fit returned')
    return True


def summarize(profile, inputs, index, path):
    bid = int(profile['bin_id'])
    va, vb = profile['velocities']
    best = refit.best_fraction(profile)
    missing = int((~np.isfinite(profile['chi2'])).sum())
    row = dict(BIN_ID=bid, SELECTED=True, USABLE=best is not None,
               VA_GLOBAL_KMS=va, VB_GLOBAL_KMS=vb,
               SIGMA_A_KMS=np.nan, SIGMA_B_KMS=np.nan, FA_RH3=np.nan, FA_RH3_ACHIEVED=np.nan,
               V_COMBINED_KMS=np.nan, SIGMA_COMBINED_KMS=np.nan, CHI2_TOTAL=np.nan,
               CHI2_SCRIPT03_MIN=inputs['coarse_min'][bid],
               CHI2_SCRIPT04_INTERPOLATED=inputs['interpolated_chi2'][index],
               DELTA_CHI2_SCRIPT03_MIN=np.nan, N_MISSING_FRACTIONS=missing,
               N_ATTEMPTS=sum(map(len, profile['attempts'])), PPXF_STATUS=-999,
               PROFILE=str(path), FLAGS='')
    flags = []
    edge = inputs['contract']['grid_edge_cells']
    for v, axis in zip((va, vb), inputs['axes'][:2]):
        i = min(edge, len(axis)-1)
        if v <= axis[i] or v >= axis[-1-i]:
            flags.append('GRID_EDGE_WARNING')
            break
    if missing:
        flags.append('INCOMPLETE_PROFILE')
    if best is None:
        flags.append('NO_USABLE_FIT')
    else:
        f = profile['achieved_fraction'][best]
        sa, sb = profile['sigma'][best]
        chi2 = profile['chi2'][best]
        winner = profile['attempts'][best][profile['winning_attempt'][best]]
        row.update(SIGMA_A_KMS=sa, SIGMA_B_KMS=sb, FA_RH3=profile['fraction'][best], FA_RH3_ACHIEVED=f,
                   V_COMBINED_KMS=f*va+(1-f)*vb,
                   SIGMA_COMBINED_KMS=np.sqrt(f*sa**2+(1-f)*sb**2+f*(1-f)*(va-vb)**2),
                   CHI2_TOTAL=chi2, DELTA_CHI2_SCRIPT03_MIN=chi2-inputs['coarse_min'][bid],
                   PPXF_STATUS=winner['status'])
        flags.extend(winner['flags'])
        if best in (0, len(refit.FRACTIONS)-1):
            flags.append('FRACTION_EDGE')
    row['FLAGS'] = '|'.join(flags)
    return row


def run_profiles(run, inputs, fingerprint, workers, interactive, logger):
    bins_dir = run.products_dir / 'RH3_fraction_profiles'
    bins_dir.mkdir(exist_ok=True)
    ids = inputs['ids']
    paths = [bins_dir / f'bin_{bid:04d}.npz' for bid in ids]
    jobs = [i for i, path in enumerate(paths) if not path.exists()]
    rows = {}
    pending = 0
    index_for_bin = {int(bid): i for i, bid in enumerate(ids)}

    def finish_bin(profile, index, path):
        nonlocal interactive, pending
        bid = int(profile['bin_id'])
        interactive = recover(profile, path, fingerprint, inputs, interactive, logger)
        if profile['pending'].any():
            pending += int(profile['pending'].sum())
            return
        figure_path = run.figures_dir / 'RH3_fraction_bins' / path.with_suffix('.pdf').name
        if not figure_path.exists():
            plots.bin_diagnostic(figure_path, profile, spectrum_for(inputs, bid), inputs['setup']['wavelength'])
        rows[bid] = summarize(profile, inputs, index, path.relative_to(run.run_dir))
        logger.info('Bin %d finalized | %d attempts | %d missing fractions | %s', bid,
                    rows[bid]['N_ATTEMPTS'], rows[bid]['N_MISSING_FRACTIONS'],
                    rows[bid]['FLAGS'] or 'no best-fit flags', extra={'file_only': True})

    logger.info('%d selected bins | %d saved profiles | %d initial profiles to calculate | %d workers',
                len(ids), len(ids)-len(jobs), len(jobs), workers)
    for index, (bid, path) in enumerate(zip(ids, paths)):
        if path.exists():
            profile = refit.load_profile(path, fingerprint, bid, inputs['velocities'][index])
            finish_bin(profile, index, path)

    if jobs:
        progress = BinProgress(len(ids), len(ids)-len(jobs), min(workers, len(jobs)))
        logger.info('Progress counts saved initial bin scans. Retry decisions are handled separately.')
        progress.refresh()
        try:
            for profile in initial_profiles(inputs, jobs, progress.workers, progress):
                progress.clear()
                bid = int(profile['bin_id'])
                index = index_for_bin[bid]
                path = paths[index]
                refit.save_profile(path, profile, fingerprint)
                progress.checkpoint(profile['elapsed_seconds'])
                logger.info('Bin %d initial scan %d/%d | %.3f s/fit | elapsed %.1fs | initial-fit ETA %.1fs',
                            bid, progress.done, len(jobs), progress.fit_seconds/(progress.done*len(refit.FRACTIONS)),
                            time.perf_counter()-progress.started, progress.eta(), extra={'file_only': True})
                finish_bin(profile, index, path)
                if progress.done < len(jobs):
                    progress.refresh()
        finally:
            progress.clear()
        logger.info('Initial scans saved for all %d selected bins.', len(ids))
    if pending:
        logger.warning('%d fraction decisions are pending. Resume in an interactive terminal.', pending)
        return 2

    # Keep original bin IDs and make out-of-domain rows explicit for maps and BL.
    exemplar = next(iter(rows.values()))
    missing_row = {key: np.nan for key in exemplar}
    missing_row.update(SELECTED=False, USABLE=False, PPXF_STATUS=-999, PROFILE='',
                       FLAGS='OUTSIDE_SELECTED_MODEL', N_ATTEMPTS=0, N_MISSING_FRACTIONS=0)
    table = Table(rows=[rows.get(i, dict(missing_row, BIN_ID=i)) for i in range(len(inputs['coarse_min']))])
    table.meta.update(model=inputs['fit']['name'], fingerprint=fingerprint,
                      systemic_velocity_kms=inputs['fit']['theta'][4],
                      velocity_frame='Original Script 03 LOS frame',
                      fraction_definition=inputs['contract']['fraction_definition'])
    for col in ('VA_GLOBAL_KMS', 'VB_GLOBAL_KMS', 'SIGMA_A_KMS', 'SIGMA_B_KMS', 'V_COMBINED_KMS', 'SIGMA_COMBINED_KMS'):
        table[col].unit = 'km / s'
    table_path = run.products_dir / 'RH3_global_kinematics.ecsv'
    temp = table_path.with_suffix('.ecsv.tmp')
    table.write(temp, format='ascii.ecsv', overwrite=True)
    temp.replace(table_path)
    plots.kinematic_maps(run.figures_dir / 'RH3_global_kinematics.pdf', table, inputs)
    plots.major_axis_profile(run.figures_dir / 'RH3_major_axis_velocity_profile.pdf', inputs)
    save_json(run.metadata_dir / 'script05_manifest.json',
              dict(script=SCRIPT, fingerprint=fingerprint, contract=inputs['contract'],
                   input_paths=inputs['paths'], n_selected_bins=len(ids),
                   n_usable_bins=int(sum(table['USABLE'])),
                   n_incomplete_profiles=int(sum(table['N_MISSING_FRACTIONS'] > 0)),
                   products=dict(kinematics_table=table_path.relative_to(run.run_dir),
                                 fraction_profiles=bins_dir.relative_to(run.run_dir)),
                   uncertainty='Conditional exact RH3 fits at the selected coherent velocities'))
    logger.info('Saved Script 05 handoff: %s', table_path)
    return 0 if np.any(table['USABLE']) else 1


def main():
    args = arguments()
    cfg = load_config(args.config, validate=False)
    workers = args.workers if args.workers is not None else getattr(cfg, 'N_WORKERS', None)
    workers = 1 if workers is None else int(workers)
    if workers < 1:
        raise ValueError('--workers must be at least 1')
    overrides = {n: getattr(args, f'script{n}_run') for n in (1, 2, 3)}
    inputs = rh3_refit_io.read_inputs(cfg, args.script4_run, overrides)
    contract = inputs['contract']
    contract['code'] = {Path(path).name: file_hash(path) for path in (
        __file__, refit.__file__, rh3_refit_io.__file__, ppxf_utils.__file__
    )}
    contract['versions'] = {name: version(name) for name in ('numpy', 'scipy', 'ppxf', 'capfit')}
    fingerprint = hashlib.sha256(json.dumps(contract, sort_keys=True, allow_nan=False).encode()).hexdigest()

    directory = args.output_dir.expanduser().resolve()
    if not args.resume:
        directory.mkdir(parents=True, exist_ok=False)
    elif not directory.is_dir():
        raise FileNotFoundError(f'Script 05 run directory does not exist: {directory}')
    with output_lock(directory):
        run = RunContext(directory, *(directory/name for name in (
            'logs', 'products', 'figures', 'diagnostics', 'metadata'
        )))
        state = run.metadata_dir / 'script05_resume_state.json'
        if args.resume and load_json(state)['fingerprint'] != fingerprint:
            raise ValueError('Script 05 inputs or fit settings changed. Use a new output directory.')
        if not args.resume:
            for path in (run.logs_dir, run.products_dir, run.figures_dir, run.diagnostics_dir, run.metadata_dir):
                path.mkdir()
            snapshot_config(cfg, run.metadata_dir)
            write_config_manifest(cfg, run.metadata_dir / 'config_manifest.json')
            save_json(state, dict(fingerprint=fingerprint, contract=contract))
        master_logger = setup_pipeline_logger(run)
        if sys.stderr.isatty():
            for handler in master_logger.handlers:
                if not isinstance(handler, logging.FileHandler):
                    handler.addFilter(lambda record: not getattr(record, 'file_only', False))
        logger = setup_step_logger(run, SCRIPT)
        log_section(logger, 'Script 05 exact RH3 extraction')
        logger.info('Model %s | MaNGA dispersion seed %.3f km/s per disk | sigma bounds %s',
                    inputs['fit']['name'], inputs['setup']['baseline'], inputs['setup']['sigma_bounds'])
        nfit = len(inputs['ids'])*len(refit.FRACTIONS)
        setup = inputs['setup']
        size = 8*nfit*(3*len(setup['wavelength'])+len(setup['component'])
                       +max(0, setup['degree']+1)+setup['mdegree'])
        logger.info('%d nominal fits | %.2f GB uncompressed spectral products before plots',
                    nfit, size/1e9)
        try:
            return run_profiles(run, inputs, fingerprint, workers,
                                sys.stdin.isatty() and not args.non_interactive, logger)
        except (Exception, KeyboardInterrupt):
            logger.exception('Script 05 interrupted. Saved bin profiles remain available with --resume.')
            raise


if __name__ == '__main__':
    sys.exit(main())
