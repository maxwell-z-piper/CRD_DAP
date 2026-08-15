"""Scientifically important validation checks used throughout CRD_DAP."""

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
    }
    frame = aliases.get(frame, frame)
    if frame not in _ALLOWED_VELOCITY_FRAMES:
        raise ValueError(
            "Velocity frame must be one of heliocentric, barycentric, none, or topocentric"
        )
    return frame


def infer_velocity_frame_from_header(header: fits.Header) -> str | None:
    """Infer a standard spectral reference frame when FITS metadata make it explicit."""
    for key in ("SPECSYS", "SSYSOBS"):
        if key not in header:
            continue
        text = str(header[key]).strip().lower()
        for token, normalized in (
            ("bary", "barycentric"),
            ("helio", "heliocentric"),
            ("topo", "topocentric"),
        ):
            if token in text:
                return normalized
    return None


def validate_script1_conventions(
    header: fits.Header,
    *,
    science_medium: str,
    template_medium: str,
    science_velocity_frame: str,
) -> ConventionCheck:
    """Validate explicit wavelength/frame metadata before Script 1 continues.

    A science/template air-vacuum *difference* is permitted at this stage only
    because template preparation will perform the explicit conversion before any
    pPXF fit.  The required conversion is returned and logged.  What is not
    permitted is an unknown science convention or a header/config velocity-frame
    contradiction.
    """
    sci_medium = normalize_wavelength_medium(science_medium)
    temp_medium = normalize_wavelength_medium(template_medium)

    if str(science_velocity_frame).strip().upper() == "UNKNOWN":
        inferred = infer_velocity_frame_from_header(header)
        if inferred is None:
            raise ValueError(
                "SCIENCE_VELOCITY_FRAME is UNKNOWN and the FITS header does not "
                "provide a recognized SPECSYS/SSYSOBS frame. Set it explicitly in the config."
            )
        sci_frame = inferred
    else:
        sci_frame = normalize_velocity_frame(science_velocity_frame)
        inferred = infer_velocity_frame_from_header(header)
        if inferred is not None and inferred != sci_frame:
            raise ValueError(
                f"Config science velocity frame is {sci_frame!r}, but FITS metadata "
                f"indicate {inferred!r}. Resolve this before analysis."
            )

    return ConventionCheck(
        science_medium=sci_medium,
        template_medium=temp_medium,
        template_conversion_required=(sci_medium != temp_medium),
        science_velocity_frame=sci_frame,
        header_velocity_frame=inferred,
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
