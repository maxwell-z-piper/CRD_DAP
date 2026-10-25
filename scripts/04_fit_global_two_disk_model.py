#!/usr/bin/env python3
"""Script 04: measured stellar PSF, reviewed regularity, eight coherent fits.

Run --phase psf first; inspect both stellar candidates. Choose --psf-profile
explicitly for calibration. The calibration writes the agreed 2x6 review sheet.
The eight-model comparison never starts automatically: a reviewed coefficient
must first be selected. Script 05 will refit the selected continuous velocities
with pPXF; this stage only evaluates the saved scalar-diagonal likelihood cubes.
"""
from __future__ import annotations

import os
for _variable in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS',
                  'NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[_variable]='1'

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime,timezone
import hashlib
import json
import logging
from pathlib import Path
import shutil
import sys

import numpy as np
from astropy.table import Table

from crd_utils.config import load_config, snapshot_config, write_config_manifest
from crd_utils import global_disk, global_disk_io as data_io
from crd_utils import global_disk_plotting as plots
from crd_utils import stellar_psf


def arguments():
    parser=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--config',required=True,type=Path)
    parser.add_argument('--output-dir',required=True,type=Path)
    parser.add_argument('--phase',choices=('psf','calibrate','review','compare','choose-model'),default='calibrate')
    parser.add_argument('--script1-run',type=Path,help='Explicit local Script 01 directory; works after moving Windows products')
    parser.add_argument('--script2-run',type=Path)
    parser.add_argument('--script3-run',type=Path)
    parser.add_argument('--psf-profile',choices=('gaussian','moffat'),help='Explicit choice AFTER inspecting stellar QA')
    parser.add_argument('--workers',type=int,default=1,help='Parallel independent global fits, with one BLAS thread each')
    parser.add_argument('--resume',action='store_true',help='Reuse only fingerprint-matching completed fits')
    parser.add_argument('--baseline-coefficient',type=float,help='Center C=lambda*A for a new strength sweep')
    parser.add_argument('--coefficient',type=float,help='Reviewed C=lambda*A to adopt, including zero')
    parser.add_argument('--note',default='',help='Human rationale retained with a coefficient/model selection')
    parser.add_argument('--interactive',action='store_true',help='Three-choice terminal coefficient review')
    parser.add_argument('--model-name',help='One named comparison model to pass to Script 05')
    return parser.parse_args()


def setup_output(path):
    for name in ('metadata','products','figures','logs','fits'):
        (path/name).mkdir(parents=True,exist_ok=True)
    logger=logging.getLogger('CRD_DAP.04')
    logger.setLevel(logging.INFO)
    logger.propagate=False
    for handler in logger.handlers[:]:
        handler.close();logger.removeHandler(handler)
    formatter=logging.Formatter('%(asctime)s | %(levelname)-8s | %(name)s | %(message)s')
    for handler in (logging.StreamHandler(sys.stdout),logging.FileHandler(path/'logs/04_fit_global_two_disk_model.log')):
        handler.setFormatter(formatter);logger.addHandler(handler)
    return logger


def psf_phase(cfg,out,resume,logger):
    paths=[Path(path).expanduser().resolve() for path in cfg.SCRIPT04_PSF_STANDARD_CUBES]
    if not paths:
        raise ValueError('Set SCRIPT04_PSF_STANDARD_CUBES to reduced native standard-star cubes')
    observed=cfg.SCRIPT04_PSF_OBSERVED_RANGE_ANGSTROM
    if observed is None:
        observed=np.asarray(cfg.RH3_SN_REST_RANGE_ANGSTROM)*(1+cfg.REDSHIFT)
    request=dict(paths=[str(p) for p in paths],sha256=[data_io.file_hash(p) for p in paths],
                 observed_range=list(map(float,observed)),cutout=float(cfg.SCRIPT04_PSF_CUTOUT_ARCSEC),
                 quadrature=int(cfg.SCRIPT04_PSF_QUADRATURE),code=data_io.file_hash(Path(stellar_psf.__file__)))
    path=out/'products/stellar_psf_measurements.json'
    if path.exists():
        saved=data_io.load_json(path)
        if not resume or saved['request']!=request:
            raise ValueError('PSF output exists; use --resume with identical inputs or a new output directory')
        logger.info('Reusing measured standard-star candidates')
    else:
        measurements=[]
        for path_in in paths:
            logger.info('Measuring both PSF candidates: %s',path_in)
            measurements.append(stellar_psf.measure_star(path_in,observed,
                cutout_arcsec=cfg.SCRIPT04_PSF_CUTOUT_ARCSEC,quadrature=cfg.SCRIPT04_PSF_QUADRATURE))
        saved=dict(request=request,measurements=measurements,
                   flags=['STANDARD_NOT_CONTEMPORANEOUS_WITH_SCIENCE',
                          'STANDARD_TO_SCIENCE_PSF_TRANSFER_APPROXIMATION'])
        data_io.save_json(path,saved)
    stellar_psf.plot_stars(saved['measurements'],out/'figures/stellar_PSF_QA.png')
    logger.info('Inspect figures/stellar_PSF_QA.pdf; then explicitly choose --psf-profile gaussian or moffat')


def problem_setup(cfg,args,out,logger):
    if any(value is None for value in (args.script1_run,args.script2_run,args.script3_run)):
        raise ValueError('Supply --script1-run, --script2-run and --script3-run explicitly')
    data_io.required_measurements(cfg)
    if args.psf_profile is None:
        raise ValueError('Inspect stellar_PSF_QA.pdf and supply --psf-profile; no profile is assumed')
    psf_path=out/'products/stellar_psf_measurements.json'
    if not psf_path.exists():
        raise ValueError('First run --phase psf in this output directory')
    inputs=data_io.read_inputs(cfg,args.script1_run,args.script2_run,args.script3_run)
    measured=data_io.load_json(psf_path)
    kernel,psf=stellar_psf.effective_kernel(measured['measurements'],args.psf_profile,inputs['header'],
        radius_arcsec=cfg.SCRIPT04_PSF_KERNEL_RADIUS_ARCSEC,quadrature=cfg.SCRIPT04_PSF_QUADRATURE)
    cfg_values={key:value for key,value in cfg.as_dict().items() if key.startswith((
        'SCRIPT04_','MANGA_SIGMA_','R_2SIGMA_','PA_KIN_','Q0_','RING_','R_START_',
        'CENTER_BOUND_','VSYS_BOUND_','INCLINATION_PRIOR_','X0_INITIAL_','Y0_INITIAL_','VSYS_INITIAL_'))}
    # Scalar input values and stage-04 code matter; upstream path strings that
    # have moved between Windows/Mac do not replace verified product contents.
    contract=dict(provenance=inputs['provenance'],config=cfg_values,psf=psf,
                  psf_product_hash=data_io.file_hash(psf_path),
                  kernel_sha256=hashlib.sha256(kernel.tobytes()).hexdigest(),
                  code={p.name:data_io.file_hash(p) for p in (
                      Path(__file__),Path(global_disk.__file__),Path(data_io.__file__),Path(stellar_psf.__file__))})
    serialized=json.dumps(contract,sort_keys=True,default=str)
    fingerprint=hashlib.sha256(serialized.encode()).hexdigest()
    existing=out/'metadata/script04_inputs.json'
    if existing.exists() and data_io.load_json(existing)['fingerprint']!=fingerprint:
        raise ValueError('Stage-04 scientific inputs/code changed. Use a new output directory; completed fits are not interchangeable.')
    if not existing.exists():
        snapshot_config(cfg,out/'metadata')
        write_config_manifest(cfg,out/'metadata/config_manifest.json')
        data_io.save_json(existing,dict(fingerprint=fingerprint,contract=json.loads(serialized)))
    # Validate reuse BEFORE writing adopted products. A rejected changed run
    # must not replace the kernel belonging to existing completed fits.
    kernel_path=out/'products/adopted_RH3_PSF.npz'
    if not kernel_path.exists():
        temporary=kernel_path.with_suffix('.tmp.npz')
        np.savez_compressed(temporary,kernel=kernel,light=inputs['light'],
                            east=inputs['east'],north=inputs['north'])
        temporary.replace(kernel_path)
    if not (out/'metadata/adopted_RH3_PSF.json').exists():
        data_io.save_json(out/'metadata/adopted_RH3_PSF.json',psf)
    problems={};infos={}
    for extent in ('two_sigma','full'):
        for use_psf in (False,True):
            key=(extent,use_psf)
            problems[key],infos[key]=data_io.make_problem(inputs,cfg,kernel if use_psf else None,
                                                         psf['major_fwhm_arcsec'],extent)
            logger.info('%s PSF=%s: %d bins, %d nodes/disk, radial spacing %.4f arcsec, Rmax %.4f arcsec',
                        extent,use_psf,len(problems[key].selected),problems[key].nnode,
                        infos[key]['spacing_arcsec'],infos[key]['Rmax_arcsec'])
            if use_psf:
                plots.geometry_preview(problems[key],infos[key],out/f'figures/{extent}_radial_geometry.png')
    return inputs,problems,infos,fingerprint,psf


def _fit_job(payload):
    problem,coefficient,options,logpath=payload
    handler=logging.FileHandler(logpath)
    handler.setFormatter(logging.Formatter('%(asctime)s | %(message)s'))
    logger=logging.getLogger('crd_utils.global_disk')
    logger.handlers=[handler];logger.setLevel(logging.INFO)
    try:
        return global_disk.fit_disk(problem,coefficient,**options)
    finally:
        handler.close();logger.handlers=[]


def fit_options(cfg):
    return dict(seed=int(cfg.SCRIPT04_RANDOM_SEED),population_factor=int(cfg.SCRIPT04_POPULATION_FACTOR),
                global_iterations=int(cfg.SCRIPT04_GLOBAL_ITERATIONS),local_starts=int(cfg.SCRIPT04_LOCAL_STARTS),
                local_iterations=int(cfg.SCRIPT04_LOCAL_ITERATIONS),tolerance=float(cfg.SCRIPT04_OPTIMIZER_TOLERANCE))


def run_jobs(jobs,cfg,args,out,fingerprint,logger):
    """Atomic fit checkpoints, one parallel task per entire global fit."""
    completed={};pending=[]
    for name,problem,coefficient,info in jobs:
        path=out/f'fits/{name}.json'
        if path.exists():
            result=data_io.load_json(path)
            if not args.resume or result['fingerprint']!=fingerprint or result['coefficient']!=coefficient:
                raise ValueError(f'Existing fit {name}: use --resume with matching inputs or a new directory')
            completed[name]=result
            logger.info('Reusing completed fit %s',name)
        else:
            pending.append((name,problem,coefficient,info))
    def save(name,problem,coefficient,info,result):
        result.update(name=name,fingerprint=fingerprint,info=info,sign_a=problem.sign_a,
                      **{'lambda':float(coefficient/info['A'])})
        if info['label']['review_required']:
            result['flags'].append('DISK_LABEL_REFERENCE_REVIEW_REQUIRED')
        if problem.kernel is not None:
            result['flags'] += data_io.load_json(out/'metadata/adopted_RH3_PSF.json')['flags']
            from scipy.signal import fftconvolve
            coverage=fftconvolve((problem.light>0).astype(float),problem.kernel,mode='same')
            y,x,b=problem.member_y,problem.member_x,problem.member_bin
            affected=np.bincount(b,weights=(coverage[y,x]<1-1e-6),minlength=len(problem.bin_ids))>0
            result['psf_edge_affected_bins']=problem.bin_ids[problem.selected][affected[problem.selected]].tolist()
            if result['psf_edge_affected_bins']:
                result['flags'].append('PSF_SPATIAL_EDGE_OR_MASK_AFFECTED')
        else:
            result['psf_edge_affected_bins']=[]
        result['flags']=sorted(set(result['flags']))
        result['n_bins']=len(problem.selected)
        result['n_parameters']=len(problem.bounds)
        result['mean_bin_data_score']=result['data_score']/len(problem.selected)
        data_io.save_json(out/f'fits/{name}.json',result)
        completed[name]=result
        logger.info('Completed %s: Qdata %.6g, penalty %.6g, %.1fs; flags=%s',
                    name,result['data_score'],result['penalty'],result['elapsed_seconds'],result['flags'])
    if args.workers<1:
        raise ValueError('--workers must be at least 1')
    options=fit_options(cfg)
    if args.workers==1:
        for name,problem,coefficient,info in pending:
            logger.info('Starting %s; optimizer history: logs/%s.progress.log',name,name)
            result=_fit_job((problem,coefficient,options,out/f'logs/{name}.progress.log'))
            save(name,problem,coefficient,info,result)
    elif pending:
        with ProcessPoolExecutor(max_workers=min(args.workers,len(pending))) as pool:
            futures={pool.submit(_fit_job,(problem,c,options,out/f'logs/{name}.progress.log')):(name,problem,c,info)
                     for name,problem,c,info in pending}
            for future in as_completed(futures):
                save(*futures[future],future.result())
    return completed


def coefficient_name(value):
    return 'calibration_C'+format(float(value),'.12g')


def calibrate(cfg,args,out,problems,infos,fingerprint,logger):
    baseline=args.baseline_coefficient
    baseline=float(cfg.SCRIPT04_BASELINE_COEFFICIENT if baseline is None else baseline)
    if not np.isfinite(baseline) or baseline<=0:
        raise ValueError('Calibration baseline coefficient must be positive')
    values=np.r_[0.,baseline*2.**np.arange(-5,6)]
    problem,info=problems['full',True],infos['full',True]
    jobs=[(coefficient_name(c),problem,float(c),info) for c in values]
    fitted=run_jobs(jobs,cfg,args,out,fingerprint,logger)
    results=[fitted[job[0]] for job in jobs]
    control=results[0]
    if min(r['data_score'] for r in results)<control['data_score']-1e-7*max(1.,abs(control['data_score'])):
        if 'CONTROL_SEARCH_NEEDS_REFINEMENT' not in control['flags']:
            control['flags'].append('CONTROL_SEARCH_NEEDS_REFINEMENT')
            data_io.save_json(out/f'fits/{control["name"]}.json',control)
        logger.warning('A penalized run has lower data cost than the control: inspect alternative minima/convergence before selection')
    plots.calibration_summary(problem,results,out/'figures/smoothing_calibration_2x6.png',info)
    record=dict(fingerprint=fingerprint,baseline_coefficient=baseline,
                coefficients=values.tolist(),fits=[r['name'] for r in results],
                reference=info,figure='figures/smoothing_calibration_2x6.pdf')
    data_io.save_json(out/'products/smoothing_calibration.json',record)
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    data_io.save_json(out/f'metadata/smoothing_sweep_{stamp}.json',record)
    logger.info('Calibration finished. Review the ONE-PAGE 2x6 PDF, then run --phase review. The eight fits have NOT started.')


def review(cfg,args,out,problems,infos,fingerprint,logger):
    record=data_io.load_json(out/'products/smoothing_calibration.json')
    if record['fingerprint']!=fingerprint:
        raise ValueError('Calibration fingerprint mismatch')
    coefficient=args.coefficient
    if args.interactive:
        print('Inspect figures/smoothing_calibration_2x6.pdf before choosing.')
        print('1: weaker range; 2: stronger range; 3: adopt a tested value (0 = control)')
        choice=input('Choice [1/2/3]: ').strip()
        if choice in ('1','2'):
            args.baseline_coefficient=record['baseline_coefficient']*(.25 if choice=='1' else 4.)
            args.resume=True
            calibrate(cfg,args,out,problems,infos,fingerprint,logger)
            return
        if choice!='3':
            raise ValueError('No selection was made')
        print('Tested effective C values:',record['coefficients'])
        coefficient=float(input('Adopt C=lambda*A (not an extent-specific lambda): '))
        args.note=input('Reason for selection: ').strip()
    if coefficient is None or not np.isfinite(coefficient) or coefficient<0:
        raise ValueError('Supply --coefficient C, or --interactive for the three-choice review')
    path=out/f'fits/{coefficient_name(coefficient)}.json'
    if not path.exists():
        raise ValueError('This value has not been fitted/reviewed. Run --phase calibrate --baseline-coefficient VALUE --resume, inspect it, then select.')
    result=data_io.load_json(path)
    if result['fingerprint']!=fingerprint:
        raise ValueError('Selected fit belongs to different scientific inputs')
    selection=dict(fingerprint=fingerprint,coefficient=float(coefficient),reference_lambda=result['lambda'],
                   reference_A=result['info']['A'],fit=result['name'],note=args.note,
                   flags=result['flags'],selected_utc=datetime.now(timezone.utc).isoformat(),
                   interpretation='empirical human-reviewed regularity choice; not a posterior probability or optimality proof')
    data_io.save_json(out/'products/smoothing_selection.json',selection)
    logger.info('Saved reviewed C=%.6g. Run --phase compare explicitly to launch the eight configurations.',coefficient)


def export_fit(result,out):
    table=Table()
    table['BIN_ID']=result['bin_ids']
    for name,key in (('VA_GLOBAL_KMS','va'),('VB_GLOBAL_KMS','vb'),('FA_RH3_COARSE','fraction'),
                     ('DATA_DELTA_CHI2','bin_score'),('SHEAR_A_KMS','shear_a'),('SHEAR_B_KMS','shear_b'),
                     ('RADIAL_PADDING_LIGHT_FRACTION','radial_padding_fraction')):
        table[name]=result[key]
    table['PSF_EDGE_AFFECTED']=np.isin(result['bin_ids'],result.get('psf_edge_affected_bins',[]))
    table.meta.update(model=result['name'],fingerprint=result['fingerprint'],
                      velocity_convention='LOS km/s in the original Script 03 frame; no systemic subtraction applied',
                      fraction_status='coarse profile; exact refinement belongs to Script 05')
    path=out/f'products/{result["name"]}_bin_predictions.ecsv'
    table.write(path,format='ascii.ecsv',overwrite=True)
    return path


def compare(cfg,args,out,problems,infos,fingerprint,logger):
    selection=data_io.load_json(out/'products/smoothing_selection.json')
    if selection['fingerprint']!=fingerprint:
        raise ValueError('Reviewed coefficient belongs to different scientific inputs')
    coefficient=selection['coefficient']
    jobs=[]
    for extent in ('two_sigma','full'):
        for psf in (False,True):
            for regular in (False,True):
                name=f'{extent}_psf{int(psf)}_smooth{int(regular)}'
                value=coefficient if regular else 0.
                path=out/f'fits/{name}.json'
                # Reuse the two full-aperture PSF pilot solutions exactly when
                # their scientific fingerprint matches this comparison.
                pilot=out/f'fits/{coefficient_name(value)}.json'
                if extent=='full' and psf and not path.exists() and pilot.exists():
                    result=data_io.load_json(pilot)
                    if result['fingerprint']==fingerprint:
                        result['name']=name;data_io.save_json(path,result)
                jobs.append((name,problems[extent,psf],value,infos[extent,psf]))
    # Existing matching pilot copies are intentional reuse, not user overwrite.
    args.resume=True
    results=run_jobs(jobs,cfg,args,out,fingerprint,logger)
    for name,problem,_,_ in jobs:
        export_fit(results[name],out)
        plots.solution_diagnostics(problem,results[name],out/f'figures/{name}_diagnostics.png')
    comparison={}
    for psf in (False,True):
        for regular in (False,True):
            short=results[f'two_sigma_psf{int(psf)}_smooth{int(regular)}']
            full=results[f'full_psf{int(psf)}_smooth{int(regular)}']
            full_cost=dict(zip(full['bin_ids'],full['bin_score']))
            shared=np.array([full_cost[b] for b in short['bin_ids']])
            comparison[f'psf{int(psf)}_smooth{int(regular)}']=dict(
                common_bin_ids=short['bin_ids'],restricted_inner_data_score=short['data_score'],
                common_bin_count=len(shared),
                restricted_mean_inner_data_score=short['data_score']/len(shared),
                full_mean_inner_data_score=float(shared.mean()),
                restricted_parameter_count=short['n_parameters'],
                full_parameter_count=full['n_parameters'],
                full_inner_data_score=float(shared.sum()),
                full_outer_data_score=float(full['data_score']-shared.sum()),
                geometry_change=(np.asarray(full['theta'][:5])-short['theta'][:5]).tolist())
    plots.extent_comparison(results,out/'figures/radial_extent_rotationcurve_comparison.png')
    data_io.save_json(out/'products/eight_model_comparison.json',dict(fingerprint=fingerprint,
        models=list(results),smoothing_selection=selection,common_inner_comparisons=comparison,
        interpretation='compare common-domain data scores and geometry/curves; no automatic preferred model'))
    logger.info('Eight configurations saved. Choose a model explicitly with --phase choose-model --model-name NAME.')


def main():
    args=arguments()
    cfg=load_config(args.config,validate=False)
    out=args.output_dir.expanduser().resolve()
    logger=setup_output(out)
    try:
        if args.phase=='psf':
            psf_phase(cfg,out,args.resume,logger);return
        _,problems,infos,fingerprint,_=problem_setup(cfg,args,out,logger)
        if args.phase=='calibrate':
            calibrate(cfg,args,out,problems,infos,fingerprint,logger)
        elif args.phase=='review':
            review(cfg,args,out,problems,infos,fingerprint,logger)
        elif args.phase=='compare':
            compare(cfg,args,out,problems,infos,fingerprint,logger)
        else:
            comparison=data_io.load_json(out/'products/eight_model_comparison.json')
            if comparison['fingerprint']!=fingerprint or args.model_name not in comparison['models']:
                raise ValueError('Choose a named model from this run\'s eight-model comparison')
            result=data_io.load_json(out/f'fits/{args.model_name}.json')
            prediction_path=export_fit(result,out)
            data_io.save_json(out/'products/selected_global_model.json',dict(
                fingerprint=fingerprint,model=args.model_name,note=args.note,
                bin_predictions=str(prediction_path),fit=str(out/f'fits/{args.model_name}.json'),
                stage05_contract='Exact pPXF refits at VA_GLOBAL_KMS,VB_GLOBAL_KMS; fractions here remain coarse',
                noise='Use final saved Script 03 noise unchanged; do not apply scalar multipliers twice'))
            logger.info('Selected model recorded for Script 05: %s',args.model_name)
    except Exception:
        logger.exception('Script 04 stopped; completed fit checkpoints remain available')
        raise


if __name__=='__main__':
    main()
