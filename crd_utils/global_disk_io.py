"""Input/provenance contracts for Script 04, including portable run paths."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import warnings

import numpy as np
from astropy.io import fits

from .global_disk import DiskProblem, disk_coordinates, radial_nodes
from .geometry import sample_inclination_prior


def file_hash(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):
            digest.update(block)
    return digest.hexdigest()


def save_json(path,value):
    """Replace a completed small product atomically; never leave half a fit."""
    path=Path(path)
    path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_suffix(path.suffix+'.tmp')
    def clean(item):
        if isinstance(item,dict):
            return {key:clean(value) for key,value in item.items()}
        if isinstance(item,(list,tuple,np.ndarray)):
            return [clean(value) for value in item]
        if isinstance(item,(float,np.floating)) and not np.isfinite(item):
            return None
        if isinstance(item,np.generic):
            return item.item()
        if isinstance(item,Path):
            return str(item)
        return item
    temp.write_text(json.dumps(clean(value),indent=2,allow_nan=False)+'\n')
    temp.replace(path)


def load_json(path):
    return json.loads(Path(path).read_text())


def required_measurements(cfg):
    names=('MANGA_SIGMA_PEAK_POS_KMS','MANGA_SIGMA_PEAK_NEG_KMS','MANGA_SIGMA_BASE_KMS',
           'R_2SIGMA_POS_INITIAL_ARCSEC','R_2SIGMA_NEG_INITIAL_ARCSEC')
    missing=[name for name in names if getattr(cfg,name,None) is None]
    if missing:
        raise ValueError('Fill these reviewed MaNGA config inputs before calibration: '+', '.join(missing))
    values={name:float(getattr(cfg,name)) for name in names}
    if any(not np.isfinite(v) or v<=0 for v in values.values()):
        raise ValueError('MaNGA dispersions and peak distances must be finite and positive')
    base=values[names[2]]
    if min(values[names[0]],values[names[1]])<=base:
        raise ValueError('Both dispersion peaks must exceed the baseline for the equal-light seed estimate')
    return values


def inclination_constraints(cfg):
    """Convert existing PyMorph/q0 inputs into the approved hard interval.

    Explicit inclination inputs can override the catalog for another target.
    Otherwise reuse the existing project's q/q0 uncertainty propagation; save
    its catalog band/HDU and clipping convention instead of hiding them.
    """
    manual=getattr(cfg,'SCRIPT04_INCLINATION_BOUNDS_DEG',None)
    if manual is not None:
        lo,hi=map(float,manual)
        initial=getattr(cfg,'SCRIPT04_INCLINATION_INITIAL_DEG',None)
        initial=(lo+hi)/2 if initial is None else float(initial)
        source=dict(source='explicit hard inclination bounds')
    else:
        if cfg.SCRIPT04_PYMORPH_HDU is None:
            raise ValueError('Select SCRIPT04_PYMORPH_HDU after reviewing the catalog band, or supply explicit inclination bounds')
        hdu=int(cfg.SCRIPT04_PYMORPH_HDU)
        with fits.open(cfg.PYMORPH_VAC) as hdus:
            data=hdus[hdu].data
            matches=np.array([str(value).strip()==str(cfg.TARGET_NAME).strip() for value in data['PLATEIFU']])
            if np.sum(matches)!=1:
                raise ValueError('PyMorph PLATEIFU must match exactly one target row; set explicit inclination bounds otherwise')
            row=data[matches][0]
            if 'FLAG_FAILED_SE' in data.names and int(row['FLAG_FAILED_SE'])!=0:
                raise ValueError('PyMorph two-component photometry failed; supply reviewed inclination bounds')
            q,err=float(row['BA_SE_DISK']),float(row['BA_SE_DISK_ERR'])
        if not np.isfinite(q+err) or not 0<q<=1 or err<0:
            raise ValueError('Invalid PyMorph disk axis ratio/error')
        draws=sample_inclination_prior(q,err,cfg.Q0_INTRINSIC_DISK,cfg.Q0_INTRINSIC_DISK_ERR,
                                      n_draws=cfg.N_INCLINATION_PRIOR_DRAWS,seed=cfg.SCRIPT04_RANDOM_SEED)
        p16,initial,p84=np.percentile(draws,[16,50,84])
        factor=float(cfg.INCLINATION_PRIOR_SIGMA_FACTOR)
        lo=max(0.,initial-factor*(initial-p16))
        hi=min(np.nextafter(90.,0.),initial+factor*(p84-initial))
        source=dict(source='PyMorph disk axis ratio and existing q0 propagation',hdu=hdu,
                    q=q,q_error=err,q0=float(cfg.Q0_INTRINSIC_DISK),
                    q0_error=float(cfg.Q0_INTRINSIC_DISK_ERR),
                    sampling='existing Gaussian draws with physical clipping',
                    percentiles_16_50_84=[float(p16),float(initial),float(p84)])
    if not (0<=lo<initial<hi<90):
        raise ValueError('Inclination constraints must satisfy 0 <= lower < initial < upper < 90 degrees')
    return float(initial),(float(lo),float(hi)),source


def read_inputs(cfg,script1_run,script2_run,script3_run):
    """Read frozen upstream products; explicit paths support moved Windows runs.

    Only selected continuum channels are read from the prepared cube, avoiding
    a second full in-memory science cube. Its original GOODMASK is respected.
    Absolute paths in the original manifests remain historical provenance.
    """
    s1,s2,s3=map(lambda p:Path(p).expanduser().resolve(),(script1_run,script2_run,script3_run))
    m2=load_json(s2/'metadata/script02_manifest.json')
    m3=load_json(s3/'metadata/script03_manifest.json')
    if str(m3['target'])!=str(cfg.TARGET_NAME) or str(m2['target'])!=str(cfg.TARGET_NAME):
        raise ValueError('Target mismatch between config and upstream runs')
    if not np.isclose(float(m3['redshift']),float(cfg.REDSHIFT),rtol=0,atol=1e-10):
        raise ValueError('Current redshift differs from the selected Script 03 run')
    cube_path=s3/'products/RH3_likelihood_cubes.npz'
    with np.load(cube_path) as z:
        chi=np.asarray(z['chi2_total'],float)
        va,vb,fa=(np.asarray(z[name],float) for name in ('VA_grid','VB_grid','fA_grid'))
        one_velocity=np.asarray(z['one_velocity'],float)
    if chi.ndim!=4 or chi.shape[1:]!=(len(va),len(vb),len(fa)):
        raise ValueError('Malformed Script 03 likelihood axes')
    if any(len(a)<2 or not np.all(np.isfinite(a)) or np.any(np.diff(a)<=0) for a in (va,vb,fa)):
        raise ValueError('Likelihood coordinates must be finite and strictly increasing')
    finite=np.isfinite(chi)
    if not np.all(np.any(finite,axis=(1,2,3))):
        raise ValueError('At least one bin has no valid spectral states; no automatic bin dropping is allowed')
    reference=np.min(np.where(finite,chi,np.inf),axis=(1,2,3))
    chi=chi-reference[:,None,None,None]
    membership_path=s2/'products/master_bin_membership.npz'
    with np.load(membership_path) as z:
        membership={name:z[name].copy() for name in z.files}
    with np.load(s1/'products/RH3_spatial_coordinates.npz') as z:
        east,north=np.asarray(z['x_arcsec']),np.asarray(z['y_arcsec'])
    prepared_path=s1/'products/prepared_RH3.fits'
    with fits.open(prepared_path,memmap=True) as hdus:
        header=hdus[0].header.copy()
        wave=np.asarray(hdus['WAVELENGTH'].data)
        lo,hi=m2['rh3_window_observed_angstrom']
        use=(wave>=lo)&(wave<=hi)&np.asarray(hdus['GOODWAVE'].data,bool)
        if use.sum()<2:
            raise ValueError('Recorded RH3 continuum window has too few usable channels')
        data=np.asarray(hdus[0].data[use],float)
        good=np.asarray(hdus['GOODMASK'].data[use],bool)&np.isfinite(data)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore',RuntimeWarning)
            continuum=np.nanmedian(np.where(good,data,np.nan),axis=0)
        spatial=np.asarray(hdus['GOODSPAX'].data,bool)
    if continuum.shape!=east.shape or north.shape!=east.shape:
        raise ValueError('Prepared cube and spatial coordinates have different shapes')
    # Match Script 02's nonnegative light-weight convention, preserving GLOBAL
    # amplitudes. No per-bin renormalization may precede spatial convolution.
    light=np.where(spatial&np.isfinite(continuum),np.maximum(continuum,0.),0.)
    if light.max()<=0:
        raise ValueError('The RH3 continuum has no positive light')
    light/=light.max()
    y,x,b=(np.asarray(membership[name],int) for name in ('rh3_y_pix','rh3_x_pix','rh3_bin_id'))
    nbin=len(chi)
    if np.any(b<0) or np.any(b>=nbin) or not np.array_equal(np.unique(b),np.arange(nbin)):
        raise ValueError('Script 02 membership does not cover exactly the Script 03 bins')
    if not (np.allclose(east[y,x],membership['rh3_x_arcsec'],atol=1e-5) and
            np.allclose(north[y,x],membership['rh3_y_arcsec'],atol=1e-5)):
        raise ValueError('Member coordinates do not match the prepared RH3 geometry')
    flux=np.bincount(b,weights=light[y,x],minlength=nbin)
    if np.any(flux<=0):
        raise ValueError('Some bins have no positive continuum for PSF weighting; inspect the continuum map')
    initial_inc,inc_bounds,inc_source=inclination_constraints(cfg)
    vsys=getattr(cfg,'VSYS_INITIAL_KMS',None)
    if vsys is None:
        values=one_velocity[np.isfinite(one_velocity)]
        if not len(values):
            raise ValueError('Supply VSYS_INITIAL_KMS; no finite one-component velocities are available')
        vsys=float(np.median(values))
    geometry=np.array([cfg.X0_INITIAL_ARCSEC or 0.,cfg.Y0_INITIAL_ARCSEC or 0.,
                       cfg.PA_KIN_INITIAL_DEG,initial_inc,vsys],float)
    pa_half=float(cfg.PA_BOUND_SIGMA*cfg.PA_KIN_ERR_DEG)
    center=float(cfg.CENTER_BOUND_ARCSEC)
    bound=float(cfg.VSYS_BOUND_KMS)
    bounds=np.array([[geometry[0]-center,geometry[0]+center],[geometry[1]-center,geometry[1]+center],
                     [geometry[2]-pa_half,geometry[2]+pa_half],inc_bounds,[vsys-bound,vsys+bound]])
    if np.any(~np.isfinite(bounds)) or np.any(bounds[:,1]<=bounds[:,0]):
        raise ValueError('Geometry hard bounds must be finite and have nonzero width')
    provenance=dict(source_script1_run=str(s1),source_script2_run=str(s2),source_script3_run=str(s3),
                    historical_script1_path=m3.get('source_script1_run'),
                    historical_script2_path=m3.get('source_script2_run'),
                    likelihood_sha256=file_hash(cube_path),membership_sha256=file_hash(membership_path),
                    continuum_sha256=hashlib.sha256(light.tobytes()).hexdigest(),
                    spatial_sha256=hashlib.sha256(east.tobytes()+north.tobytes()).hexdigest(),
                    continuum_observed_angstrom=[float(lo),float(hi)],
                    inclination=inc_source,geometry_initial=geometry.tolist(),
                    geometry_bounds=bounds.tolist(),noise='saved final scalar-calibrated diagonal chi-square; no rescaling')
    return dict(east=east,north=north,light=light,member_y=y,member_x=x,member_bin=b,
                chi2=chi,va_axis=va,vb_axis=vb,fractions=fa,one_velocity=one_velocity,
                geometry=geometry,geometry_bounds=bounds,header=header,provenance=provenance,
                reference_chi2=reference)


def make_problem(inputs,cfg,kernel,psf_fwhm,extent):
    """Freeze whole-bin selection at the reference geometry for one extent."""
    measurements=required_measurements(cfg)
    initial=inputs['geometry']
    radius,cosine=disk_coordinates(inputs['east']-initial[0],inputs['north']-initial[1],initial[2],initial[3])
    y,x,b=inputs['member_y'],inputs['member_x'],inputs['member_bin']
    nbin=len(inputs['chi2'])
    cut=max(measurements['R_2SIGMA_POS_INITIAL_ARCSEC'],measurements['R_2SIGMA_NEG_INITIAL_ARCSEC'])
    spacing=float(cfg.RING_SPACE_ARCSEC or psf_fwhm)
    first=float(cfg.R_START_FWHM_FACTOR*psf_fwhm)
    if not np.isfinite(cfg.SCRIPT04_VTYP_KMS) or cfg.SCRIPT04_VTYP_KMS<=0:
        raise ValueError('SCRIPT04_VTYP_KMS must be finite and positive')
    if extent=='two_sigma':
        initial_nodes=radial_nodes(first,spacing,cut,outward=False)
        nominal_edge=initial_nodes[-1]+float(cfg.RING_DELTA_FACTOR)*spacing
        # Include complete bins touching the nominal domain; never cut spectra.
        touching=np.bincount(b,weights=(radius[y,x]<=nominal_edge),minlength=nbin)>0
        selected=np.flatnonzero(touching)
        if not len(selected):
            raise ValueError('No whole bins intersect the two-sigma domain')
        target=float(np.max(radius[y[np.isin(b,selected)],x[np.isin(b,selected)]]))
        nodes=radial_nodes(first,spacing,max(initial_nodes[-1],target),outward=True)
    elif extent=='full':
        selected=np.arange(nbin)
        nodes=radial_nodes(first,spacing,float(np.max(radius[y,x])),outward=True)
        nominal_edge=float(nodes[-1])
    else:
        raise ValueError('Unknown radial extent')
    speed_peak=np.mean([measurements['MANGA_SIGMA_PEAK_POS_KMS'],measurements['MANGA_SIGMA_PEAK_NEG_KMS']])
    base=measurements['MANGA_SIGMA_BASE_KMS']
    los_seed=np.sqrt(speed_peak**2-base**2)
    seed_speed=los_seed/np.sin(np.deg2rad(initial[3]))
    explicit=getattr(cfg,'SCRIPT04_DISK_A_ROTATION_SIGN',None)
    if explicit is None:
        # Reference side is negative signed major-axis coordinate. Identify
        # its outer single-component sense without relabeling by f_A.
        major=radius*cosine
        weight=np.bincount(b,weights=inputs['light'][y,x],minlength=nbin)
        location=np.bincount(b,weights=inputs['light'][y,x]*major[y,x],minlength=nbin)/weight
        reference=np.flatnonzero((location<0)&np.isfinite(inputs['one_velocity']))
        if not len(reference):
            raise ValueError('Cannot assign disk identity; supply SCRIPT04_DISK_A_ROTATION_SIGN after inspection')
        bid=int(reference[np.argmin(location[reference])])
        value=inputs['one_velocity'][bid]-initial[4]
        if value==0:
            raise ValueError('AMBIGUOUS_DISK_LABEL: set SCRIPT04_DISK_A_ROTATION_SIGN explicitly')
        sign=-int(np.sign(value))
        label=dict(source='outermost negative-major-axis bin single-component sign',reference_bin=bid,
                   reference_los_offset=float(value),sign_a=sign,review_required=True)
    else:
        sign=int(explicit)
        if sign not in (-1,1):
            raise ValueError('SCRIPT04_DISK_A_ROTATION_SIGN must be +1 or -1')
        label=dict(source='explicit reviewed sign',sign_a=sign,review_required=False)
    vmax=float(cfg.SCRIPT04_VROT_MAX_KMS)
    if not np.isfinite(vmax) or vmax<=0:
        raise ValueError('SCRIPT04_VROT_MAX_KMS must be positive')
    bounds=np.vstack((inputs['geometry_bounds'],np.tile([0.,vmax],(2*len(nodes),1))))
    problem=DiskProblem(**{key:inputs[key] for key in ('east','north','light','member_y','member_x','member_bin',
        'va_axis','vb_axis','fractions')},chi2=inputs['chi2'][selected],selected=selected,
        bin_ids=np.arange(nbin),nodes=nodes,bounds=bounds,initial_geometry=initial,
        initial_speed=min(seed_speed,vmax),sign_a=sign,kernel=kernel)
    info=dict(extent=extent,nominal_edge_arcsec=float(nominal_edge),nodes=nodes.tolist(),
              selected_bins=selected.tolist(),label=label,initial_speed_kms=float(seed_speed),
              spacing_arcsec=spacing,Rmax_arcsec=float(nodes[-1]),
              A=float((nodes[-1]/spacing)**4/cfg.SCRIPT04_VTYP_KMS**2),
              light_weighting='fixed globally normalized nonnegative observed RH3 continuum proxy',
              boundary_rule='whole bins intersecting reference-geometry domain; outward node coverage; constant outer continuation')
    return problem,info
