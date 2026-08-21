"""Scientifically important validation checks used throughout CRD_DAP.

Script 1 deliberately treats wavelength-frame and calibration-configuration
checks as hard scientific validation, not cosmetic metadata checks.  A small
silent mismatch in wavelength medium, radial-velocity correction, slicer,
grating, or central wavelength can propagate directly into the two-component
kinematic inference.

The functions in this module are intentionally conservative: when the metadata
are sufficiently explicit to infer a convention, the configured value is
checked against the FITS header.  When ``"auto"`` is requested, the header must
contain enough information to infer the convention unambiguously or the script
fails and asks the user to specify it explicitly.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord


_ALLOWED_WAVELENGTH_MEDIA = {"air", "vacuum"}
_ALLOWED_VELOCITY_FRAMES = {"heliocentric", "barycentric", "none", "topocentric"}


@dataclass(frozen=True)
class ConventionCheck:
    science_medium: str
    template_medium: str
    template_conversion_required: bool
    science_velocity_frame: str
    header_velocity_frame: str | None
    header_wavelength_medium: str | None


@dataclass(frozen=True)
class CalibrationMatch:
    """Summary of an arc/science configuration match.

    ``expected_grating`` is the grating requested by the target configuration.
    It is stored explicitly so the run products retain the distinction between
    the pipeline arm label (``BL``/``RH3``) and the actual instrument grating
    used for a particular dataset.
    """

    arm: str
    expected_grating: str
    science_camera: str
    arc_camera: str
    science_grating: str
    arc_grating: str
    science_ifu: str
    arc_ifu: str
    science_binning: str
    arc_binning: str
    science_cwave_angstrom: float | None
    arc_cwave_angstrom: float | None


def normalize_wavelength_medium(value: str) -> str:
    medium = str(value).strip().lower()
    if medium not in _ALLOWED_WAVELENGTH_MEDIA:
        raise ValueError("Wavelength medium must be 'air' or 'vacuum'")
    return medium


def normalize_velocity_frame(value: str) -> str:
    frame = str(value).strip().lower()
    aliases = {
        "helio": "heliocentric",
        "bary": "barycentric",
        "barycent": "barycentric",
        "heliocen": "heliocentric",
        "topocent": "topocentric",
        "topo": "topocentric",
    }
    frame = aliases.get(frame, frame)
    if frame not in _ALLOWED_VELOCITY_FRAMES:
        raise ValueError(
            "Velocity frame must be one of heliocentric, barycentric, none, or topocentric"
        )
    return frame


def _spectral_fits_axis(header: fits.Header) -> int | None:
    """Return the FITS axis carrying the wavelength coordinate when identifiable."""
    ndim = int(header.get("NAXIS", 0))
    for axis in range(1, ndim + 1):
        ctype = str(header.get(f"CTYPE{axis}", "")).strip().upper()
        if "AWAV" in ctype or "WAVE" in ctype:
            return axis
    return None


def infer_wavelength_medium_from_header(header: fits.Header) -> str | None:
    """Infer air/vacuum wavelength medium from KCWI/FITS spectral metadata.

    KCWI DRP cube construction writes ``CTYPE='AWAV'`` for air wavelengths.
    After the DRP's air-to-vacuum correction the uploaded real-world test cube
    uses ``CTYPE='WAVE'`` with the card comment ``Vacuum Wavelengths``.  We use
    those explicit metadata first and only fall back to a small set of common
    convention keywords if they are present.
    """
    axis = _spectral_fits_axis(header)
    if axis is not None:
        key = f"CTYPE{axis}"
        ctype = str(header.get(key, "")).strip().upper()
        if "AWAV" in ctype:
            return "air"

        # WAVE by itself is not universally sufficient to prove vacuum, so
        # inspect the FITS-card comment written by the KCWI DRP as well.
        try:
            comment = str(header.comments[key]).strip().lower()
        except Exception:
            comment = ""
        if "vacuum" in comment:
            return "vacuum"
        if "air wavelength" in comment:
            return "air"

    # Optional explicit convention keywords used by some FITS producers.
    for key in ("WAVEMED", "WAVEMEDM", "AIRORVAC", "VACUUM"):
        if key not in header:
            continue
        raw = header[key]
        if isinstance(raw, (bool, np.bool_)):
            if key == "VACUUM":
                return "vacuum" if bool(raw) else "air"
            continue
        text = str(raw).strip().lower()
        if "vac" in text:
            return "vacuum"
        if "air" in text:
            return "air"

    return None


def infer_velocity_frame_from_header(header: fits.Header) -> str | None:
    """Infer the spectral reference frame from standard/KCWI metadata.

    In addition to FITS ``SPECSYS``/``SSYSOBS``, KCWI DRP 1.2 products record
    the actually applied radial-velocity correction in ``VCORRTYP``.  The
    latter is therefore checked explicitly.
    """
    for key in ("VCORRTYP", "SPECSYS", "SSYSOBS"):
        if key not in header:
            continue
        text = str(header[key]).strip().lower()
        if "bary" in text:
            return "barycentric"
        if "helio" in text:
            return "heliocentric"
        if "topo" in text:
            return "topocentric"
        if text in {"none", "uncorrected", "native"}:
            return "none"
    return None


def validate_script1_conventions(
    header: fits.Header,
    *,
    science_medium: str,
    template_medium: str,
    science_velocity_frame: str,
) -> ConventionCheck:
    """Validate wavelength medium and velocity-frame metadata.

    ``science_medium`` and ``science_velocity_frame`` may be set to ``"auto"``
    (or the legacy value ``"UNKNOWN"`` for the velocity frame).  Auto mode does
    *not* guess: it requires the FITS header to identify the convention.  If an
    explicit config value is supplied and the header independently identifies a
    contradictory convention, Script 1 hard-fails.

    A science/template air-vacuum difference is allowed at Script 1 because the
    template-preparation stage will perform an explicit conversion before pPXF.
    The mismatch is recorded as ``template_conversion_required``.
    """
    header_medium = infer_wavelength_medium_from_header(header)
    requested_medium = str(science_medium).strip().lower()
    if requested_medium in {"auto", "unknown"}:
        if header_medium is None:
            raise ValueError(
                "SCIENCE_WAVELENGTH_MEDIUM is auto/unknown but the FITS header does not "
                "unambiguously identify air versus vacuum wavelengths. Set it explicitly."
            )
        sci_medium = header_medium
    else:
        sci_medium = normalize_wavelength_medium(requested_medium)
        if header_medium is not None and header_medium != sci_medium:
            raise ValueError(
                f"Config science wavelength medium is {sci_medium!r}, but FITS metadata "
                f"indicate {header_medium!r}. Resolve this before analysis."
            )

    temp_medium = normalize_wavelength_medium(template_medium)

    header_frame = infer_velocity_frame_from_header(header)
    requested_frame = str(science_velocity_frame).strip().lower()
    if requested_frame in {"auto", "unknown"}:
        if header_frame is None:
            raise ValueError(
                "SCIENCE_VELOCITY_FRAME is auto/unknown and the FITS header does not "
                "provide a recognized VCORRTYP/SPECSYS/SSYSOBS frame. Set it explicitly."
            )
        sci_frame = header_frame
    else:
        sci_frame = normalize_velocity_frame(requested_frame)
        if header_frame is not None and header_frame != sci_frame:
            raise ValueError(
                f"Config science velocity frame is {sci_frame!r}, but FITS metadata "
                f"indicate {header_frame!r}. Resolve this before analysis."
            )

    return ConventionCheck(
        science_medium=sci_medium,
        template_medium=temp_medium,
        template_conversion_required=(sci_medium != temp_medium),
        science_velocity_frame=sci_frame,
        header_velocity_frame=header_frame,
        header_wavelength_medium=header_medium,
    )


def check_wavelength_convention(science_medium: str, template_medium: str) -> None:
    """Hard-fail if two arrays about to be fit use different wavelength media."""
    s = normalize_wavelength_medium(science_medium)
    t = normalize_wavelength_medium(template_medium)
    if s != t:
        raise ValueError(
            f"Science wavelengths are {s} while templates are {t}. "
            "Convert one dataset before pPXF fitting."
        )


def _clean_header_text(value) -> str:
    return str(value).strip().upper()


def _grating_for_arm(header: fits.Header, arm: str) -> str:
    arm_u = str(arm).strip().upper()
    if arm_u == "BL":
        return _clean_header_text(header.get("BGRATNAM", ""))
    if arm_u == "RH3":
        return _clean_header_text(header.get("RGRATNAM", ""))
    raise ValueError("Calibration matching currently expects arm='BL' or arm='RH3'")


def _cwave_for_arm(header: fits.Header, arm: str) -> float | None:
    key = "BCWAVE" if str(arm).strip().upper() == "BL" else "RCWAVE"
    if key not in header:
        return None
    try:
        return float(header[key])
    except Exception:
        return None


def validate_arc_science_configuration(
    science_header: fits.Header,
    arc_header: fits.Header,
    *,
    arm: str,
    expected_grating: str | None = None,
    central_wavelength_tolerance_angstrom: float = 0.5,
) -> CalibrationMatch:
    """Hard-check that a master arc matches the science instrumental setup.

    Parameters
    ----------
    science_header, arc_header
        FITS headers for the science cube and candidate master arc.
    arm
        Pipeline arm label. ``BL`` selects blue-camera grating metadata and
        ``RH3`` selects red-camera grating metadata. The arm label intentionally
        remains stable even when Script 1 is exercised with a different grating.
    expected_grating
        Grating required for this run, e.g. ``"BL"``, ``"BM"``, ``"RH3"``, or
        ``"RL"``. If omitted, the legacy defaults are retained: ``BL`` for the
        blue arm and ``RH3`` for the red arm.
    central_wavelength_tolerance_angstrom
        Maximum allowed science/arc central-wavelength difference.

    Notes
    -----
    Making the expected grating configurable does *not* weaken calibration
    validation. The science grating must match the configured grating, the arc
    grating must match the configured grating, and the science and arc must
    still match one another. Camera, IFU slicer, detector binning, and central
    wavelength are checked as before.
    """
    arm_u = str(arm).strip().upper()
    if arm_u not in {"BL", "RH3"}:
        raise ValueError("Calibration matching currently expects arm='BL' or arm='RH3'")

    expected_camera = "BLUE" if arm_u == "BL" else "RED"

    if expected_grating is None:
        expected_grating_u = "BL" if arm_u == "BL" else "RH3"
    else:
        expected_grating_u = _clean_header_text(expected_grating)
        if not expected_grating_u:
            raise ValueError("expected_grating must be a non-empty grating name")

    sci_camera = _clean_header_text(science_header.get("CAMERA", ""))
    arc_camera = _clean_header_text(arc_header.get("CAMERA", ""))
    if sci_camera and sci_camera != expected_camera:
        raise ValueError(
            f"{arm_u} science cube reports CAMERA={sci_camera!r}, expected {expected_camera!r}."
        )
    if arc_camera and arc_camera != expected_camera:
        raise ValueError(
            f"{arm_u} master arc reports CAMERA={arc_camera!r}, expected {expected_camera!r}."
        )
    if sci_camera and arc_camera and sci_camera != arc_camera:
        raise ValueError(f"Science/arc CAMERA mismatch: {sci_camera!r} versus {arc_camera!r}")

    sci_grating = _grating_for_arm(science_header, arm_u)
    arc_grating = _grating_for_arm(arc_header, arm_u)
    if sci_grating and sci_grating != expected_grating_u:
        raise ValueError(
            f"Configured {arm_u} science cube reports grating {sci_grating!r}, "
            f"but {arm_u}_EXPECTED_GRATING requires {expected_grating_u!r}. "
            "Check the target configuration path and configured expected grating."
        )
    if arc_grating and arc_grating != expected_grating_u:
        raise ValueError(
            f"Configured {arm_u} master arc reports grating {arc_grating!r}, "
            f"but {arm_u}_EXPECTED_GRATING requires {expected_grating_u!r}. "
            "Use the matching grating calibration."
        )
    if sci_grating and arc_grating and sci_grating != arc_grating:
        raise ValueError(
            f"{arm_u} science/arc grating mismatch: {sci_grating!r} versus {arc_grating!r}. "
            "Use a master arc from the exact science grating configuration."
        )

    sci_ifu = _clean_header_text(science_header.get("IFUNAM", ""))
    arc_ifu = _clean_header_text(arc_header.get("IFUNAM", ""))
    if sci_ifu and arc_ifu and sci_ifu != arc_ifu:
        raise ValueError(f"Science/arc IFUNAM mismatch: {sci_ifu!r} versus {arc_ifu!r}")

    sci_bin = _clean_header_text(science_header.get("BINNING", ""))
    arc_bin = _clean_header_text(arc_header.get("BINNING", ""))
    if sci_bin and arc_bin and sci_bin != arc_bin:
        raise ValueError(f"Science/arc BINNING mismatch: {sci_bin!r} versus {arc_bin!r}")

    sci_cwave = _cwave_for_arm(science_header, arm_u)
    arc_cwave = _cwave_for_arm(arc_header, arm_u)
    if sci_cwave is not None and arc_cwave is not None:
        if abs(sci_cwave - arc_cwave) > float(central_wavelength_tolerance_angstrom):
            raise ValueError(
                f"Science/arc central wavelength mismatch for {arm_u}: "
                f"{sci_cwave:.4f} versus {arc_cwave:.4f} Angstrom."
            )

    return CalibrationMatch(
        arm=arm_u,
        expected_grating=expected_grating_u,
        science_camera=sci_camera,
        arc_camera=arc_camera,
        science_grating=sci_grating,
        arc_grating=arc_grating,
        science_ifu=sci_ifu,
        arc_ifu=arc_ifu,
        science_binning=sci_bin,
        arc_binning=arc_bin,
        science_cwave_angstrom=sci_cwave,
        arc_cwave_angstrom=arc_cwave,
    )


def grid_edge_distance_cells(index: tuple[int, ...], shape: tuple[int, ...]) -> int:
    """Minimum number of cells between a selected index and any grid edge."""
    if len(index) != len(shape):
        raise ValueError("index and shape dimensionality differ")
    distances = []
    for i, n in zip(index, shape):
        if not 0 <= i < n:
            raise IndexError("Grid index lies outside shape")
        distances.extend([i, n - 1 - i])
    return int(min(distances))


def is_grid_edge_warning(
    index: tuple[int, ...],
    shape: tuple[int, ...],
    warning_cells: int = 2,
) -> bool:
    """True when selected state is within ``warning_cells`` of any edge."""
    return grid_edge_distance_cells(index, shape) <= warning_cells


def sky_separation_arcsec(a: SkyCoord, b: SkyCoord) -> float:
    """Angular separation between two centers in arcsec."""
    return float(a.separation(b).arcsec)


def compare_center_estimates(
    estimates: Iterable[SkyCoord],
    *,
    warning_arcsec: float,
) -> tuple[float, bool]:
    """Return maximum pairwise center separation and warning state."""
    coords = list(estimates)
    if len(coords) < 2:
        return 0.0, False
    separations = []
    for i in range(len(coords)):
        for j in range(i + 1, len(coords)):
            separations.append(coords[i].separation(coords[j]).arcsec)
    maximum = float(np.max(separations))
    return maximum, maximum > float(warning_arcsec)
