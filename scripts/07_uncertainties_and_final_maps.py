#!/usr/bin/env python3
"""Refit paired RH3 and BL noise realizations and export final uncertainties."""
from __future__ import annotations

import os
for _name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
              'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_name] = '1'

import argparse
from contextlib import contextmanager
import importlib.util
import json
from importlib.metadata import version
import logging
import multiprocessing
from pathlib import Path
import shutil
import sys
import time

import numpy as np

from crd_utils import global_disk, ppxf_grid, ppxf_utils, populations, rh3_refit, templates
from crd_utils import script07_fits as fits
from crd_utils import script07_io as io
from crd_utils import script07_products as products
from crd_utils.bl_population_io import fingerprint
from crd_utils.bl_population_products import write_npz
from crd_utils.config import load_config, snapshot_config
from crd_utils.global_disk_io import file_hash, load_json, save_json
from crd_utils.runtime import BinProgress, output_lock
from crd_utils.script07_resume import matches_checkpoint

LOG = logging.getLogger('CRD_DAP.Script07')
SCRIPTS = Path(__file__).absolute().parent


class PendingReview(Exception):
    """Exit with a resumable checkpoint after a requested review."""


def arguments():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config', type=Path, required=True)
    for n in range(1, 7):
        p.add_argument(f'--script{n}-run', type=Path, required=n >= 4,
                       help='Required nominal run' if n >= 4 else 'Override the recorded nominal run path')
    p.add_argument('--output-dir', type=Path, required=True)
    p.add_argument('--workers', type=int, help='Parallel bins within the active draw')
    p.add_argument('--draws', type=int, help='Accepted draws to collect (default 100 from the config)')
    p.add_argument('--resume', action='store_true', help='Continue the matching ensemble at its saved stage')
    p.add_argument('--non-interactive', action='store_true', help='Move on automatically from failed fractions')
    p.add_argument('--validate-inputs', action='store_true', help='Verify nominal handoffs without starting draws')
    p.add_argument('--kroupa-library', type=Path, help='Optional matched XSL library for a separate nominal BL sensitivity fit')
    return p.parse_args()


def load_driver(stage, name):
    spec = importlib.util.spec_from_file_location(f'script07_stage{stage}', SCRIPTS/f'{name}.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@contextmanager
def bin_pool(workers, nominal):
    p = nominal.problem
    grid = dict(va_grid=p.va_axis, vb_grid=p.vb_axis, fa_grid=p.fractions,
                sigma_start_a=nominal.grid_settings['sigma_start'],
                sigma_start_b=nominal.grid_settings['sigma_start'])
    spectra = {arm: {key: state[key] for key in ('noise', 'good', 'model')}
               for arm, state in nominal.spectra.items()}
    params = (nominal.setup, spectra, grid)
    fits.initialize(*params)
    if workers == 1:
        yield None
    else:
        with multiprocessing.get_context('spawn').Pool(workers, initializer=fits.initialize_worker,
                                                        initargs=params) as pool:
            yield pool


def completed_bins(pool, function, jobs, progress):
    if pool is None:
        yield from map(function, jobs)
        return
    outcomes = pool.imap_unordered(function, jobs, chunksize=1)
    for _ in jobs:
        while True:
            try:
                outcome = outcomes.next(timeout=5.)
            except multiprocessing.TimeoutError:
                progress.refresh()
            else:
                yield outcome
                break


def stage03(scratch, nominal, pool, workers, seed, draw):
    paths = {int(bid): scratch/'03'/f'bin_{bid:04d}.npz' for bid in nominal.ids['RH3']}
    (scratch/'03').mkdir(exist_ok=True)
    jobs = [(bid, seed, draw) for bid, path in paths.items() if not path.exists()]
    LOG.info('Draw %d | 03 RH3 grid | %d/%d saved bins', draw, len(paths)-len(jobs), len(paths))
    progress = BinProgress(len(paths), len(paths)-len(jobs), workers)
    try:
        progress.refresh()
        for bid, cube, audit in completed_bins(pool, fits.grid_worker, jobs, progress):
            write_npz(paths[bid], chi2=cube, audit=json.dumps(audit, allow_nan=False))
            progress.checkpoint(audit['elapsed_seconds'])
            progress.refresh()
    finally:
        progress.clear()

    audits, cubes = [], []
    disk_checkpoint = scratch/'04.json'
    for path in paths.values():
        with np.load(path, allow_pickle=False) as archive:
            audits.append(json.loads(str(archive['audit'])))
            if not disk_checkpoint.exists():
                cubes.append(archive['chi2'])
    return np.stack(cubes) if cubes else None, audits


def stage04(scratch, nominal, draw, grids, audits):
    path = scratch/'04.json'
    if path.exists():
        LOG.info('Draw %d | 04 coherent model | reused checkpoint', draw)
        return load_json(path)
    LOG.info('Draw %d | 04 coherent model | %s', draw, nominal.fit['name'])
    missing = [row['bin_id'] for row in audits if row['minimum'] is None]
    if missing:
        result = dict(failure='One or more fixed RH3 bins have no finite grid states',
                      no_grid_support=missing, flags=['NO_VALID_COHERENT_MODEL'])
    else:
        result = fits.coherent_fit(nominal, grids)
    save_json(path, result)
    return result


def _resolve_recovery(profile, pending, driver, setup, seed, draw, checkpoint, unattended):
    if not pending and not profile['retry_queue']:
        return profile
    arm, bid = profile['arm'], profile['bin_id']
    spectrum = fits.mock_spectrum(arm, bid, seed, draw)
    while pending or profile['retry_queue']:
        if profile['retry_queue']:
            fits.run_retry_queue(profile, setup, spectrum, checkpoint)
            profile, pending = fits.advance(profile, seed, draw, unattended)
            checkpoint(profile)
            continue
        checkpoint(profile)
        if unattended:
            profile['moved_on'].extend(pending)
        else:
            if not sys.stdin.isatty():
                raise PendingReview('Resume in a terminal for fit recovery, or pass --non-interactive to move on.')
            profile['fraction'] = fits.FRACTIONS.tolist()
            try:
                if arm == 'RH3':
                    pair = driver.choose_retry(profile, pending[0], setup)
                    indices, starts = [pending[0]], None if pair is None else [pair]
                else:
                    indices = pending
                    starts = driver.retry_choice(profile, pending, setup, bool(spectrum['gas'].any()))
            except EOFError as exc:
                raise PendingReview('Input ended. The pending fit decisions are saved.') from exc
            if starts is None:
                profile['moved_on'].extend(indices)
            else:
                profile['retry_queue'] = list(zip(indices, starts))
        checkpoint(profile)
        if not profile['retry_queue']:
            profile, pending = fits.advance(profile, seed, draw, unattended)
            checkpoint(profile)
    return profile


def fraction_stage(directory, arm, cohort, available, nominal, pool, workers,
                   seed, draw, unattended, driver):
    directory.mkdir(exist_ok=True)
    paths = {bid: directory/f'bin_{bid:04d}.json' for bid in available}
    jobs, profiles = [], {}
    for bid, path in paths.items():
        kin = available[bid]
        if not path.exists():
            jobs.append((arm, bid, kin, seed, draw, unattended))
            continue
        profile = load_json(path)
        if profile['complete'] and not profile['retry_queue']:
            profiles[bid] = profile
        else:
            # Recreate pending indices without running ahead of a chosen retry.
            if profile['retry_queue']:
                pending = []
            else:
                profile, pending = fits.advance(profile, seed, draw, unattended)
            profiles[bid] = _resolve_recovery(profile, pending, driver, nominal.setup[arm],
                                              seed, draw, lambda p: save_json(path, p), unattended)
            save_json(path, profiles[bid])
    LOG.info('Draw %d | %s fractions | %d/%d saved bins | %d without an upstream fit',
             draw, arm, len(profiles), len(paths), len(cohort)-len(paths))
    progress = BinProgress(len(paths), len(profiles), workers)
    try:
        if jobs:
            progress.refresh()
        for profile, pending in completed_bins(pool, fits.fraction_worker, jobs, progress):
            bid = profile['bin_id']
            path = paths[bid]
            progress.clear()
            save_json(path, profile)
            profiles[bid] = _resolve_recovery(profile, pending, driver, nominal.setup[arm],
                                              seed, draw, lambda p: save_json(path, p), unattended)
            save_json(path, profiles[bid])
            progress.checkpoint(profiles[bid]['elapsed_seconds'])
            progress.refresh()
    finally:
        progress.clear()
    missing = [int(bid) for bid in cohort if bid not in profiles or profiles[bid]['winner'] is None]
    return profiles, missing


def run_draw(directory, nominal, pool, workers, seed, draw, unattended, drivers):
    started = time.monotonic()
    scratch = directory/'work'/f'draw_{draw:04d}'
    scratch.mkdir(parents=True, exist_ok=True)
    grids, audit03 = stage03(scratch, nominal, pool, workers, seed, draw)
    disk = stage04(scratch, nominal, draw, grids, audit03)
    del grids
    ids = nominal.ids
    audit = dict(draw=draw, seed=seed, grid=audit03, disk_flags=disk.get('flags', []),
                 disk={k: v for k, v in disk.items() if k != 'history'},
                 cohorts={a: len(v) for a, v in ids.items()})
    rh3, bl, missing = {}, {}, {}
    reasons = []
    if 'failure' in disk:
        reasons.append('NO_VALID_COHERENT_MODEL')
    else:
        available = {int(bid): [va, vb] for bid, va, vb in zip(disk['bin_ids'], disk['va'], disk['vb'])}
        rh3, missing['RH3'] = fraction_stage(scratch/'05', 'RH3', ids['RH3'], available,
                                            nominal, pool, workers, seed, draw, unattended, drivers['RH3'])
        if len(missing['RH3']) / len(ids['RH3']) >= fits.REJECT_FRACTION:
            reasons.append('RH3_UNUSABLE_FRACTION')
        else:
            available_bl = {}
            for bid in ids['BL']:
                winner = rh3[bid]['winner']
                if winner is not None:
                    va, vb = available[bid]
                    sa, sb = winner['sigma']
                    available_bl[int(bid)] = [va, sa, vb, sb]
            bl, missing['BL'] = fraction_stage(scratch/'06', 'BL', ids['BL'], available_bl,
                                               nominal, pool, workers, seed, draw, unattended, drivers['BL'])
            if len(missing['BL']) / len(ids['BL']) >= fits.REJECT_FRACTION:
                reasons.append('BL_UNUSABLE_FRACTION')
    values, weights = products.measurements(nominal, disk, rh3, bl)
    aperture, aperture_ssp, usable = products.aperture_statistics(weights, nominal.setup['BL'])
    data_file = f'draw_{draw:04d}.npz'
    write_npz(directory/'draws'/data_file, **values, bl_weights=weights,
              aperture_ssp=aperture_ssp, theta=np.asarray(disk.get('theta', [])),
              bin_ids=nominal.problem.bin_ids)
    audit.update(accepted=not reasons, reasons=reasons, missing_bins=missing, data_file=data_file,
                 data_sha256=file_hash(directory/'draws'/data_file),
                 aperture=aperture, aperture_usable_bins=usable,
                 fit_audit={a: [fits.audit_profile(profiles[b]) for b in sorted(profiles)]
                            for a, profiles in (('RH3', rh3), ('BL', bl))},
                 winning_flags={a: {str(k): flags for k, p in profiles.items()
                                    if (flags := fits.profile_flags(p))}
                                for a, profiles in (('RH3', rh3), ('BL', bl))},
                 elapsed_this_session_seconds=time.monotonic()-started)
    # Publishing the audit is the draw commit. Resume counts only commits.
    save_json(directory/'draws'/f'draw_{draw:04d}.json', audit)
    shutil.rmtree(scratch)
    return audit


def committed_draws(directory):
    attempts = sorted((load_json(path) for path in (directory/'draws').glob('draw_*.json')),
                      key=lambda record: record['draw'])
    for index, record in enumerate(attempts, 1):
        name = f'draw_{index:04d}.npz'
        path = directory/'draws'/name
        if (record['draw'] != index or record['data_file'] != name or not path.is_file()
                or file_hash(path) != record['data_sha256']):
            raise ValueError('An ensemble commit is incomplete or has a missing draw product')
        scratch = directory/'work'/f'draw_{index:04d}'
        if scratch.exists():
            shutil.rmtree(scratch)
    return attempts


def collect(directory, nominal, workers, seed, target, rejection_limit, unattended, drivers):
    attempts = committed_draws(directory)
    accepted = [record for record in attempts if record['accepted']]
    if len(accepted) >= target:
        return accepted[:target], attempts
    review_path = directory/'metadata/rejection_review.json'
    reviewed = load_json(review_path)['reviewed_through'] if review_path.exists() else 0
    latest_accepted = accepted[-1]['draw'] if accepted else 0
    streak = len(attempts)-max(latest_accepted, reviewed)
    if streak >= rejection_limit:
        # Reaching this point via --resume acknowledges the preceding pause,
        # while all attempt statuses stay in the original audit.
        reviewed = len(attempts)
        save_json(review_path, dict(reviewed_through=reviewed, previous_streak=streak))
        LOG.info('Resuming after review through attempt %d', reviewed)
    with bin_pool(workers, nominal) as pool:
        while len(accepted) < target:
            draw = len(attempts)+1
            LOG.info('Attempt %d | %d/%d accepted', draw, len(accepted), target)
            record = run_draw(directory, nominal, pool, workers, seed, draw, unattended, drivers)
            attempts.append(record)
            if record['accepted']:
                accepted.append(record)
                latest_accepted = draw
            streak = draw-max(latest_accepted, reviewed)
            LOG.info('Attempt %d %s | %d/%d accepted | consecutive rejections %d',
                     draw, 'accepted' if record['accepted'] else f'rejected {record["reasons"]}',
                     len(accepted), target, streak)
            if streak >= rejection_limit:
                raise PendingReview(f'{streak} consecutive rejections. Inspect the draw audit, then use --resume to continue.')
    return accepted, attempts


def main():
    args = arguments()
    logging.basicConfig(level=logging.INFO, format='%(asctime)s | %(levelname)s | %(message)s')
    cfg = load_config(args.config, validate=False)
    defaults = (getattr(cfg, 'SCRIPT07_N_DRAWS', 100), getattr(cfg, 'SCRIPT07_MAX_CONSECUTIVE_REJECTIONS', 10),
                getattr(cfg, 'SCRIPT07_RANDOM_SEED', 20261114))
    target, reject_limit, seed = map(int, defaults)
    target = args.draws if args.draws is not None else target
    workers = args.workers if args.workers is not None else int(getattr(cfg, 'N_WORKERS', None) or 1)
    if min(target, reject_limit, workers) < 1 or seed < 0:
        raise ValueError('Draw target, rejection limit and workers must be positive, with a nonnegative seed')
    run_paths = {n: getattr(args, f'script{n}_run') for n in range(1, 7)}
    driver04 = load_driver(4, '04_fit_global_two_disk_model')
    nominal = io.read_nominal(cfg, run_paths, driver04.fit_options)
    directory = args.output_dir.expanduser().resolve()
    roots = {p.resolve() for key, p in nominal.paths.items() if key.startswith('run')}
    roots.add(nominal.paths['cube'].parent.parent.resolve())
    if any(directory == root or root in directory.parents or directory in root.parents for root in roots):
        raise ValueError('Choose a Script 07 output directory separate from nominal runs')
    LOG.info('Verified %s | %s | %d RH3 and %d BL eligible bins',
             nominal.provenance['target'], nominal.fit['name'], len(nominal.ids['RH3']), len(nominal.ids['BL']))
    if not nominal.provenance['nominal_optimizer_has_stability_record']:
        LOG.warning('The nominal Script 04 fit predates the Powell stability checkpoint. Mock fits use the current optimizer.')
    if args.validate_inputs:
        return 0

    drivers = {'RH3': load_driver(5, '05_extract_global_RH3_solution'),
               'BL': load_driver(6, '06_fit_BL_two_component_populations')}
    provenance = dict(nominal.provenance, seed=seed, target_draws=target,
                      versions={name: version(name) for name in ('numpy', 'scipy', 'ppxf', 'capfit', 'astropy')},
                      code={Path(module.__file__).name: file_hash(module.__file__)
                            for module in (fits, io, products, global_disk, ppxf_grid, ppxf_utils,
                                           populations, rh3_refit, templates, *drivers.values())},
                      driver_sha256=file_hash(__file__))
    state = dict(fingerprint=fingerprint(provenance), contract=provenance,
                 input_paths=nominal.paths)
    directory.mkdir(parents=True, exist_ok=True)
    with output_lock(directory, '.script07.lock'):
        metadata = directory/'metadata'
        path = metadata/'script07_resume_state.json'
        if path.exists():
            if not args.resume or not matches_checkpoint(load_json(path), state):
                raise ValueError('Use --resume with matching science inputs and settings, or choose a new output directory')
        elif args.resume or any(f.name != '.script07.lock' for f in directory.iterdir()):
            raise ValueError('No resumable Script 07 run here. Choose a new output directory')
        else:
            snapshot_config(cfg, metadata)
            save_json(path, state)
        for folder in ('draws', 'work'):
            (directory/folder).mkdir(exist_ok=True)
        handler = logging.FileHandler(directory/'pipeline.log', mode='a')
        handler.setFormatter(logging.Formatter('%(asctime)s | %(levelname)s | %(message)s'))
        logging.getLogger().addHandler(handler)
        try:
            accepted, attempts = collect(directory, nominal, workers, seed, target, reject_limit,
                                         args.non_interactive, drivers)
            products.export(directory, nominal, accepted, attempts)
            if args.kroupa_library is not None:
                from crd_utils import script07_systematics
                script07_systematics.run(directory, nominal, args.kroupa_library, workers,
                                         args.non_interactive, drivers['BL'], fraction_stage, bin_pool)
                products.export_index(directory, nominal)
        except PendingReview as exc:
            LOG.warning('%s', exc)
            return 2
        except KeyboardInterrupt:
            LOG.warning('Interrupted. Use --resume to continue from saved checkpoints.')
            return 130
        finally:
            logging.getLogger().removeHandler(handler)
            handler.close()
    LOG.info('Saved %d accepted draws and final science products in %s', target, directory)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
