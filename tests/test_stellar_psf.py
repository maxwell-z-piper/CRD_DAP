"""Pixel integration and unit-flux transfer checks for non-Gaussian PSFs."""
import numpy as np
from astropy.io import fits
from crd_utils.stellar_psf import integrated_profile,effective_kernel


def test_native_pixel_integration_changes_core_and_preserves_symmetry():
    t=np.linspace(-3,3,101)
    p=[.5,.5,0.,3.]
    coarse=integrated_profile(t,np.zeros_like(t),np.diag([1.36,.29]),p,'moffat',7)
    point=integrated_profile(t,np.zeros_like(t),np.diag([.001,.001]),p,'moffat',3)
    assert coarse.max()<point.max()
    assert np.allclose(coarse,coarse[::-1])


def test_effective_kernel_is_normalized_and_keeps_sampling_provenance():
    header=fits.Header({'NAXIS':2,'CTYPE1':'RA---TAN','CTYPE2':'DEC--TAN',
                        'CRVAL1':10.,'CRVAL2':20.,'CRPIX1':1.,'CRPIX2':1.,
                        'CDELT1':-.3/3600,'CDELT2':.3/3600,'IFUPA':0.})
    measurement=dict(path='fixture',candidates={'moffat':dict(success=True,flags=[],
                      native_pixel_arcsec=[1.36,.29],parameters=[.6,.6,0.,3.])})
    kernel,summary=effective_kernel([measurement],'moffat',header,radius_arcsec=3.,quadrature=5)
    assert np.isclose(kernel.sum(),1)
    assert np.all(kernel>=0)
    assert summary['major_fwhm_arcsec']>1.
    assert 'SCIENCE_DRIZZLE_PSF_NOT_REPLAYED' in summary['flags']
