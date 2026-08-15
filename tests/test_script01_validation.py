from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits

from crd_utils.io import discover_arc_sidecars
from crd_utils.validation import (
    infer_velocity_frame_from_header,
    infer_wavelength_medium_from_header,
    validate_arc_science_configuration,
    validate_script1_conventions,
)


def _kcwi_header(*, camera="BLUE", grating="BL", ofname="kb00001.fits"):
    h = fits.Header()
    h["NAXIS"] = 3
    h["CTYPE1"] = "RA---TAN"
    h["CTYPE2"] = "DEC--TAN"
    h["CTYPE3"] = ("WAVE", "Vacuum Wavelengths")
    h["VCORRTYP"] = "heliocentric"
    h["CAMERA"] = camera
    h["IFUNAM"] = "Large"
    h["BINNING"] = "2,2"
    h["OFNAME"] = ofname
    if camera == "BLUE":
        h["BGRATNAM"] = grating
        h["BCWAVE"] = 4600.0
    else:
        h["RGRATNAM"] = grating
        h["RCWAVE"] = 8500.0
    return h


def test_real_drp_style_convention_inference():
    h = _kcwi_header()
    assert infer_wavelength_medium_from_header(h) == "vacuum"
    assert infer_velocity_frame_from_header(h) == "heliocentric"

    result = validate_script1_conventions(
        h,
        science_medium="auto",
        template_medium="air",
        science_velocity_frame="auto",
    )
    assert result.science_medium == "vacuum"
    assert result.science_velocity_frame == "heliocentric"
    assert result.template_conversion_required


def test_arc_science_grating_mismatch_is_rejected():
    science = _kcwi_header(camera="RED", grating="RH3", ofname="kr00100.fits")
    arc = _kcwi_header(camera="RED", grating="RL", ofname="kr00025.fits")
    with pytest.raises(ValueError, match="grating"):
        validate_arc_science_configuration(science, arc, arm="RH3")


def test_sidecar_discovery_can_use_ofname_when_exposure_number_differs(tmp_path: Path):
    # This reproduces the real RED behavior found in the 2024-12-26 calibration:
    # the master arc filename and geometry-map filenames can carry different
    # exposure numbers while OFNAME still traces them to the same original arc.
    arc_header = _kcwi_header(camera="RED", grating="RL", ofname="kr241226_00025.fits")
    arc_header["NAXIS"] = 2
    arc_path = tmp_path / "kr241226_00025_marc.fits"
    fits.PrimaryHDU(np.ones((5, 6), dtype=np.float32), header=arc_header).writeto(arc_path)

    for kind in ("wavemap", "slicemap", "posmap"):
        map_header = arc_header.copy()
        fits.PrimaryHDU(np.ones((5, 6), dtype=np.float32), header=map_header).writeto(
            tmp_path / f"kr241226_00027_{kind}.fits"
        )

    found = discover_arc_sidecars(arc_path)
    assert found["wavemap"].name == "kr241226_00027_wavemap.fits"
    assert found["slicemap"].name == "kr241226_00027_slicemap.fits"
    assert found["posmap"].name == "kr241226_00027_posmap.fits"
