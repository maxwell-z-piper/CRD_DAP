from pathlib import Path

import numpy as np
from astropy.io import fits

from crd_utils.io import load_kcwikit_stack, save_prepared_cube


def _header(nx=8, ny=6, nw=12):
    h = fits.Header()
    h["NAXIS"] = 3
    h["NAXIS1"] = nx
    h["NAXIS2"] = ny
    h["NAXIS3"] = nw
    h["CTYPE1"] = "RA---TAN"
    h["CTYPE2"] = "DEC--TAN"
    h["CTYPE3"] = "WAVE"
    h["CUNIT1"] = "deg"
    h["CUNIT2"] = "deg"
    h["CUNIT3"] = "Angstrom"
    h["CRVAL1"] = 119.64
    h["CRVAL2"] = 41.58
    h["CRVAL3"] = 5000.0
    h["CRPIX1"] = 4.5
    h["CRPIX2"] = 3.5
    h["CRPIX3"] = 1.0
    # 0.3 arcsec/pixel square output grid, north-up for the synthetic test.
    h["CD1_1"] = -0.3 / 3600.0
    h["CD1_2"] = 0.0
    h["CD2_1"] = 0.0
    h["CD2_2"] = 0.3 / 3600.0
    h["CD3_3"] = 1.0
    h["WAVGOOD0"] = 5001.0
    h["WAVGOOD1"] = 5010.0
    h["VCORR"] = -23.6
    h["VCORRTYP"] = "heliocentric"
    return h


def _make_stack(tmp_path: Path):
    nw, ny, nx = 12, 6, 8
    h = _header(nx=nx, ny=ny, nw=nw)

    flux = np.zeros((nw, ny, nx), dtype=np.float64)
    var = np.zeros_like(flux)
    mask = np.ones_like(flux, dtype=np.int16)
    exp = np.zeros_like(flux)

    # Real IFU footprint occupies only the central 4x4 region. The large output
    # padding should never cause all wavelength channels to be rejected.
    flux[:, 1:5, 2:6] = 2.0
    var[:, 1:5, 2:6] = 0.04
    mask[:, 1:5, 2:6] = 0
    exp[:, 1:5, 2:6] = 1200.0

    paths = {}
    for role, data in {
        "icube": flux,
        "vcube": var,
        "mcube": mask,
        "ecube": exp,
    }.items():
        path = tmp_path / f"test_{role}.fits"
        hh = h.copy()
        if role == "ecube":
            hh["BUNIT"] = "s"
        fits.PrimaryHDU(data, header=hh).writeto(path)
        paths[role] = path
    return paths


def test_load_kcwikit_stack_excludes_padding_from_channel_qc(tmp_path):
    paths = _make_stack(tmp_path)
    cube = load_kcwikit_stack(
        paths["icube"],
        paths["vcube"],
        paths["mcube"],
        paths["ecube"],
        arm="BL",
        min_good_wavelength_fraction=0.8,
        bad_channel_fraction_threshold=0.5,
    )

    assert cube.input_format == "kcwikit"
    assert cube.shape == (6, 8, 12)
    assert cube.flags is None
    assert cube.exposure is not None
    assert np.allclose(cube.uncertainty[cube.exposure > 0], 0.2)

    # Instrument-good channels survive even though most of the requested output
    # canvas is zero-exposure padding.
    assert np.sum(cube.good_wavelength) == 10
    assert np.sum(cube.good_spaxel) == 16
    assert np.allclose(cube.coverage_fraction_spaxel[1:5, 2:6], 1.0)
    assert np.allclose(cube.coverage_fraction_spaxel[0, :], 0.0)


def test_kcwikit_companion_wcs_mismatch_fails(tmp_path):
    paths = _make_stack(tmp_path)
    with fits.open(paths["vcube"], mode="update") as hdul:
        hdul[0].header["CRVAL1"] += 0.01

    try:
        load_kcwikit_stack(
            paths["icube"], paths["vcube"], paths["mcube"], paths["ecube"], arm="BL"
        )
    except ValueError as exc:
        assert "WCS does not match" in str(exc)
    else:
        raise AssertionError("Expected mismatched KcwiKit companion WCS to fail")


def test_prepared_kcwikit_cube_saves_exposure_and_coverage(tmp_path):
    paths = _make_stack(tmp_path)
    cube = load_kcwikit_stack(
        paths["icube"], paths["vcube"], paths["mcube"], paths["ecube"], arm="BL"
    )
    out = tmp_path / "prepared.fits"
    save_prepared_cube(cube, out)

    with fits.open(out) as hdul:
        assert "EXPOSURE" in hdul
        assert "COVFRAC" in hdul
        assert "FLAGS" not in hdul
        assert hdul[0].header["CRDINFMT"] == "kcwikit"
        assert hdul["EXPOSURE"].header["BUNIT"] == "s"
