"""Measure a stellar PSF before the Script 04 kinematic fit.

Gaussian and Moffat candidates are both fitted, never selected by assumption.
Fits integrate over native rectangular pixels. Kernels retain that pixel
response once, then transfer its instrument orientation to the science cube.
This is an effective standard-star PSF transfer, not a replay of the science
stack's drizzle operation; that remaining approximation is explicitly flagged.
"""
from __future__ import annotations

from pathlib import Path
import warnings

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.optimize import least_squares

from .io import wavelength_axis_from_header


def profile_value(x,y,alpha_x,alpha_y,angle,beta,kind):
    """Unit-peak continuous stellar profile in native angular coordinates."""
    c,s = np.cos(angle),np.sin(angle)
    xx,yy = c*x+s*y,-s*x+c*y
    square = (xx/alpha_x)**2+(yy/alpha_y)**2
    return np.exp(-0.5*square) if kind=='gaussian' else (1+square)**(-beta)


def integrated_profile(x,y,pixel_vectors,parameters,kind,quadrature=5):
    """Average the continuous profile over actual native angular pixel cells."""
    ax,ay,angle,beta = parameters
    result = np.zeros_like(np.asarray(x),dtype=float)
    offsets = (np.arange(quadrature)+0.5)/quadrature-0.5
    for ox in offsets:
        for oy in offsets:
            dx,dy = pixel_vectors @ np.array([ox,oy])
            result += profile_value(x+dx,y+dy,ax,ay,angle,beta,kind)
    return result/quadrature**2


def measure_star(path, observed_range, *, cutout_arcsec=5., quadrature=5):
    """Fit both PSF candidates in a specified observed-Angstrom interval.

    Input is a reduced native standard cube, not a 1-D extracted standard.
    Invalid masks/flags are excluded. This returns QA arrays as well as fit
    parameters; an optimizer's success alone does not certify a reliable PSF.
    """
    path = Path(path)
    with fits.open(path,memmap=True) as hdus:
        header = hdus[0].header.copy()
        wave = wavelength_axis_from_header(header,hdus[0].shape[0])
        use = (wave>=observed_range[0]) & (wave<=observed_range[1])
        if np.sum(use)<3:
            raise ValueError(f"PSF interval {observed_range} is not covered by {path}")
        flux = np.asarray(hdus[0].data[use],dtype=float)
        good = np.isfinite(flux)
        for name in ('MASK','FLAGS'):
            if name in hdus:
                good &= hdus[name].data[use]==0
        with warnings.catch_warnings():
            warnings.simplefilter('ignore',RuntimeWarning)
            image = np.nanmedian(np.where(good,flux,np.nan),axis=0)
    if not np.any(np.isfinite(image)):
        raise ValueError(f"No usable stellar image in {path}")
    py,px = np.unravel_index(np.nanargmax(image),image.shape)
    wcs = WCS(header).celestial
    # Coordinates in the instrument's two spatial directions. Convert the
    # fitted angle to the science orientation only when constructing a kernel.
    scales = np.linalg.norm(wcs.pixel_scale_matrix,axis=0)*3600.
    yy,xx = np.indices(image.shape)
    x,y = (xx-px)*scales[0],(yy-py)*scales[1]
    select = np.isfinite(image) & (np.hypot(x,y)<=cutout_arcsec)
    values = image[select]
    if len(values)<20:
        raise ValueError(f"Too few usable native pixels for a PSF fit: {path}")
    background = float(np.nanmedian(image))
    amplitude = float(np.nanmax(values)-background)
    if amplitude<=0:
        raise ValueError(f"Standard star has no positive peak: {path}")
    scaled = (values-background)/amplitude
    xv,yv = x[select],y[select]
    pixel_vectors = np.diag(scales)
    candidates = {}
    for kind in ('gaussian','moffat'):
        # beta is irrelevant for Gaussian; keep it out of that optimizer.
        initial = [1.,0.,0.,0.,0.7,0.7,0.]
        lower = [0.,-2*scales[0],-2*scales[1],-2.,0.03,0.03,-np.pi/2]
        upper = [10.,2*scales[0],2*scales[1],2.,cutout_arcsec,cutout_arcsec,np.pi/2]
        if kind=='moffat':
            initial += [3.]
            lower += [1.01]
            upper += [30.]
        def model(p,full=False):
            amp,cx,cy,bg,ax,ay,angle = p[:7]
            beta = p[7] if kind=='moffat' else 1.
            xx,yy = (x,y) if full else (xv,yv)
            return bg+amp*integrated_profile(xx-cx,yy-cy,pixel_vectors,
                                            [ax,ay,angle,beta],kind,quadrature)
        result = least_squares(lambda p:model(p)-scaled,initial,bounds=(lower,upper),
                               max_nfev=2000,x_scale='jac',ftol=1e-9,xtol=1e-9,gtol=1e-9)
        ax,ay,angle = result.x[4:7]
        beta = result.x[7] if kind=='moffat' else 1.
        factor = np.sqrt(8*np.log(2)) if kind=='gaussian' else 2*np.sqrt(2**(1/beta)-1)
        flags = []
        if not result.success:
            flags.append('PSF_OPTIMIZER_NOT_CONVERGED')
        span = np.asarray(upper)-lower
        if np.any(np.minimum(result.x-lower,upper-result.x)<1e-3*span):
            flags.append('PSF_PARAMETER_BOUNDARY')
        candidates[kind] = dict(parameters=[float(ax),float(ay),float(angle),float(beta)],
            native_pixel_arcsec=scales.tolist(),fwhm_intrinsic_axes_arcsec=(factor*np.array([ax,ay])).tolist(),
            residual_rms=float(np.sqrt(np.mean(result.fun**2))),success=bool(result.success),
            flags=flags,model=model(result.x,True).tolist())
    return dict(path=str(path.resolve()),object=str(header.get('OBJECT','unknown')),
                date_obs=str(header.get('DATE-OBS','')),ut=str(header.get('UT','')),
                ifu_pa_deg=float(header.get('IFUPA',0.)),observed_range=list(observed_range),
                x=x.tolist(),y=y.tolist(),image=((image-background)/amplitude).tolist(),
                candidates=candidates)


def effective_kernel(measurements,kind,science_header,*,radius_arcsec=8.,quadrature=5):
    """Transfer each measured standard profile to the science orientation.

    Average unit-integral kernels, so an arbitrary standard exposure time or
    stellar count level never sets its weight. This retains native pixel
    integration once. Unknown science-vs-standard seeing/drizzle differences
    are provenance limitations, not silently fitted extra broadening.
    """
    if kind not in ('gaussian','moffat'):
        raise ValueError('Choose gaussian or moffat after inspecting stellar QA')
    if not measurements or not np.isfinite(radius_arcsec) or radius_arcsec<=0 or quadrature<1:
        raise ValueError('Measured standards, a positive kernel radius and quadrature are required')
    matrix = WCS(science_header).celestial.pixel_scale_matrix*3600.
    scale_min = np.linalg.svd(matrix,compute_uv=False).min()
    half = int(np.ceil(radius_arcsec/scale_min))
    yy,xx = np.mgrid[-half:half+1,-half:half+1]
    east,north = np.einsum('ij,j...->i...',matrix,np.array([xx,yy]))
    kernels = []
    measured_widths = []
    science_pa = np.deg2rad(float(science_header.get('IFUPA',0.)))
    # IFU native x is east-west at PA=0; the full WCS handles output pixels.
    rx = np.cos(science_pa)*east-np.sin(science_pa)*north
    ry = np.sin(science_pa)*east+np.cos(science_pa)*north
    for measurement in measurements:
        candidate = measurement['candidates'][kind]
        if not candidate['success']:
            raise ValueError(f"Review flagged PSF candidate before using it: {measurement['path']} {candidate['flags']}")
        pixel = np.diag(candidate['native_pixel_arcsec'])
        kernel = integrated_profile(rx,ry,pixel,candidate['parameters'],kind,quadrature)
        if kernel.sum()<=0 or not np.all(np.isfinite(kernel)):
            raise ValueError('Invalid measured PSF kernel')
        kernels.append(kernel/kernel.sum())
        # Numerical half-maximum widths include the native rectangular response.
        t = np.linspace(-radius_arcsec,radius_arcsec,4001)
        axis_widths = []
        covariance=np.array([[np.sum(kernel*rx*rx),np.sum(kernel*rx*ry)],
                             [np.sum(kernel*rx*ry),np.sum(kernel*ry*ry)]])/kernel.sum()
        _,vectors=np.linalg.eigh(covariance)
        direction=vectors[:,-1]
        angle=np.arctan2(direction[1],direction[0])
        for theta in (angle,angle+np.pi/2):
            profile = integrated_profile(t*np.cos(theta),t*np.sin(theta),pixel,
                                         candidate['parameters'],kind,quadrature)
            inside = np.flatnonzero(profile>=profile.max()/2)
            axis_widths.append(float(t[inside[-1]]-t[inside[0]]))
        measured_widths.append(max(axis_widths))
    kernel = np.mean(kernels,axis=0)
    kernel /= kernel.sum()
    return kernel,dict(profile=kind,major_fwhm_arcsec=float(np.median(measured_widths)),
                       individual_major_fwhm_arcsec=measured_widths,
                       kernel_radius_arcsec=float(radius_arcsec),
                       flags=sorted(set(['STANDARD_TO_SCIENCE_PSF_TRANSFER_APPROXIMATION',
                              'SCIENCE_DRIZZLE_PSF_NOT_REPLAYED']+
                              [flag for m in measurements for flag in m['candidates'][kind]['flags']])),
                       convention='pixel-integrated standard fit; native response retained once; equal unit-flux kernel mixture')


def plot_stars(measurements,path):
    """One stellar QA sheet: data, both residual maps and angular profiles."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes = plt.subplots(len(measurements),4,figsize=(16,4*len(measurements)),squeeze=False)
    for row,m in enumerate(measurements):
        data = np.asarray(m['image'],dtype=float)
        x,y=np.asarray(m['x']),np.asarray(m['y'])
        axes[row,0].pcolormesh(x,y,data,shading='nearest',vmin=0,vmax=1)
        axes[row,0].set_title(Path(m['path']).name+'\n'+m['ut'])
        for col,kind in enumerate(('gaussian','moffat'),1):
            fit=m['candidates'][kind]
            residual=data-np.asarray(fit['model'])
            axes[row,col].pcolormesh(x,y,residual,shading='nearest',cmap='RdBu_r',vmin=-0.1,vmax=0.1)
            axes[row,col].set_title(f"{kind}: RMS={fit['residual_rms']:.3g}\n{','.join(fit['flags']) or 'fit completed'}")
        for ax in axes[row,:3]:
            ax.set(xlim=(-5,5),ylim=(-5,5),aspect='equal',xlabel='Native x (arcsec)',ylabel='Native y (arcsec)')
        r=np.hypot(x,y)
        axes[row,3].scatter(r.ravel(),data.ravel(),s=3,c='k',alpha=.3,label='data')
        for kind in ('gaussian','moffat'):
            axes[row,3].scatter(r.ravel(),np.ravel(m['candidates'][kind]['model']),s=2,alpha=.3,label=kind)
        axes[row,3].set(xlabel='Angular radius (arcsec)',ylabel='Normalized stellar light',yscale='symlog',ylim=(-.03,2),xlim=(0,5))
        axes[row,3].legend()
    fig.tight_layout()
    fig.savefig(path,dpi=160)
    fig.savefig(Path(path).with_suffix('.pdf'))
    plt.close(fig)
