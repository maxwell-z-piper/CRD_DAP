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
        validate_arc_science_configuration(
            science,
            arc,
            arm="RH3",
            expected_grating="RH3",
        )


def test_configurable_red_grating_allows_matching_rl_science_and_arc():
    """RL is valid for Script-1 testing when it is explicitly requested."""
    science = _kcwi_header(camera="RED", grating="RL", ofname="kr00171.fits")
    arc = _kcwi_header(camera="RED", grating="RL", ofname="kr00025.fits")

    result = validate_arc_science_configuration(
        science,
        arc,
        arm="RH3",
        expected_grating="RL",
    )

    assert result.expected_grating == "RL"
    assert result.science_grating == "RL"
    assert result.arc_grating == "RL"
    assert result.science_camera == "RED"
    assert result.arc_camera == "RED"


def test_configurable_blue_grating_allows_matching_non_bl_science_and_arc():
    """The blue stream is configurable too; it is not hard-coded to BL."""
    science = _kcwi_header(camera="BLUE", grating="BM", ofname="kb00118.fits")
    arc = _kcwi_header(camera="BLUE", grating="BM", ofname="kb00021.fits")

    result = validate_arc_science_configuration(
        science,
        arc,
        arm="BL",
        expected_grating="BM",
    )

    assert result.expected_grating == "BM"
    assert result.science_grating == "BM"
    assert result.arc_grating == "BM"
    assert result.science_camera == "BLUE"
    assert result.arc_camera == "BLUE"


def test_configured_expected_grating_is_still_a_hard_requirement():
    """Configurability must not become an ignore-mismatch switch."""
    science = _kcwi_header(camera="RED", grating="RL", ofname="kr00171.fits")
    arc = _kcwi_header(camera="RED", grating="RL", ofname="kr00025.fits")

    with pytest.raises(ValueError, match="RH3_EXPECTED_GRATING"):
        validate_arc_science_configuration(
            science,
            arc,
            arm="RH3",
            expected_grating="RH3",
        )


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
