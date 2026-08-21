import numpy as np
from astropy.io import fits

from crd_utils.io import load_prepared_cube


def _header():
    h = fits.Header()
    h["NAXIS"] = 3
    h["CTYPE1"] = "RA---TAN"
    h["CTYPE2"] = "DEC--TAN"
    h["CTYPE3"] = "WAVE"
    h["CRVAL1"] = 150.0
    h["CRVAL2"] = 2.0
    h["CRVAL3"] = 5000.0
    h["CRPIX1"] = 2.0
    h["CRPIX2"] = 2.0
    h["CRPIX3"] = 1.0
    h["CDELT1"] = -0.3 / 3600.0
    h["CDELT2"] = 0.3 / 3600.0
    h["CDELT3"] = 1.0
    h["CUNIT1"] = "deg"
    h["CUNIT2"] = "deg"
    h["CUNIT3"] = "Angstrom"
    h["CRDDAP"] = True
    h["CRDSTEP"] = 1
    h["CRDARM"] = "BL"
    return h


def test_load_prepared_cube_preserves_script1_goodmask(tmp_path):
    shape = (5, 3, 4)  # FITS/native order: wavelength, y, x
    flux = np.ones(shape, dtype=np.float32)
    uncert = np.ones(shape, dtype=np.float32)
    good = np.ones(shape, dtype=np.uint8)
    good[2, 1, 1] = 0
    path = tmp_path / "prepared_BL.fits"
    fits.HDUList([
        fits.PrimaryHDU(flux, header=_header()),
        fits.ImageHDU(uncert, name="UNCERT"),
        fits.ImageHDU(np.zeros(shape, dtype=np.uint8), name="MASK"),
        fits.ImageHDU(good, name="GOODMASK"),
        fits.ImageHDU(np.ones((3, 4), dtype=np.uint8), name="GOODSPAX"),
        fits.ImageHDU(np.ones((3, 4), dtype=np.float32), name="GOODFRAC"),
        fits.ImageHDU(np.ones(5, dtype=np.uint8), name="GOODWAVE"),
        fits.ImageHDU(np.arange(5000.0, 5005.0), name="WAVELENGTH"),
    ]).writeto(path)

    cube = load_prepared_cube(path, expected_arm="BL")
    assert cube.shape == (3, 4, 5)
    assert not cube.good[1, 1, 2]
    assert cube.input_format == "prepared"
