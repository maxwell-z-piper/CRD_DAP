from pathlib import Path

import numpy as np
from astropy.io import fits

from crd_utils.io import load_kcwi_cube, save_prepared_cube, wavelength_axis_from_header


def _synthetic_cube(path: Path):
    nw, ny, nx = 20, 4, 5
    flux = np.ones((nw, ny, nx), dtype=np.float32)
    uncert = np.ones_like(flux) * 0.2
    mask = np.zeros_like(flux, dtype=np.uint8)
    flags = np.zeros_like(flux, dtype=np.uint16)

    # One bad sample and one channel that is bad everywhere.
    mask[3, 1, 1] = 1
    flags[7, :, :] = 128

    h = fits.Header()
    h["NAXIS"] = 3
    h["CTYPE1"] = "RA---TAN"
    h["CTYPE2"] = "DEC--TAN"
    h["CTYPE3"] = "WAVE"
    h["CUNIT1"] = "deg"
    h["CUNIT2"] = "deg"
    h["CUNIT3"] = "Angstrom"
    h["CRVAL1"] = 150.0
    h["CRVAL2"] = 2.0
    h["CRVAL3"] = 8400.0
    h["CRPIX1"] = 3.0
    h["CRPIX2"] = 2.5
    h["CRPIX3"] = 1.0
    h["CD1_1"] = -0.0001
    h["CD2_2"] = 0.0001
    h["CD3_3"] = 0.5
    h["WAVGOOD0"] = 8401.0
    h["WAVGOOD1"] = 8408.0
    h["SPECSYS"] = "BARYCENT"

    fits.HDUList(
        [
            fits.PrimaryHDU(flux, header=h),
            fits.ImageHDU(uncert, name="UNCERT"),
            fits.ImageHDU(mask, name="MASK"),
            fits.ImageHDU(flags, name="FLAGS"),
        ]
    ).writeto(path)


def test_wavelength_axis_linear():
    h = fits.Header()
    h["NAXIS"] = 3
    h["CTYPE3"] = "WAVE"
    h["CUNIT3"] = "Angstrom"
    h["CRVAL3"] = 5000.0
    h["CRPIX3"] = 1.0
    h["CD3_3"] = 2.0
    wave = wavelength_axis_from_header(h, 4, fits_axis=3)
    assert np.allclose(wave, [5000.0, 5002.0, 5004.0, 5006.0])


def test_load_and_save_prepared_cube(tmp_path):
    source = tmp_path / "cube.fits"
    _synthetic_cube(source)
    cube = load_kcwi_cube(
        source,
        arm="RH3",
        min_good_wavelength_fraction=0.5,
        bad_channel_fraction_threshold=0.5,
    )
    assert cube.shape == (4, 5, 20)
    assert cube.good_wavelength[7] == 0
    assert np.isclose(cube.wavelength[0], 8400.0)
    assert cube.pixel_scales_arcsec()[0] > 0

    dest = tmp_path / "prepared.fits"
    save_prepared_cube(cube, dest)
    with fits.open(dest) as hdul:
        assert "GOODMASK" in hdul
        assert "GOODSPAX" in hdul
        assert "GOODWAVE" in hdul
        assert hdul[0].data.shape == (20, 4, 5)
