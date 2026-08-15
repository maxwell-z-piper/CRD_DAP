"""KCWI/KCRM FITS I/O and standardized prepared-cube products.

This module deliberately centralizes all assumptions about the KCWI DRP FITS
layout.  Downstream science code should work with :class:`KCWICube` rather than
repeating extension names, FITS-axis conventions, or mask semantics.

The current implementation targets KCWI DRP 1.2-style cube products, whose
science image is stored in ``PRIMARY`` and which may contain ``UNCERT``,
``MASK``, ``FLAGS``, and ``NOSKYSUB`` extensions.  ``UNCERT`` is treated as a
1-sigma standard-deviation image, matching the DRP's use of
``StdDevUncertainty``.

Internally CRD_DAP uses array order ``(y, x, wavelength)``.  The original FITS
axis order is retained in the object so that Script 1 can write a prepared FITS
product with the original WCS untouched.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import json
from typing import Any

import numpy as np
from astropy.io import fits
from astropy import units as u
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales


_SPECTRAL_CTYPE_TOKENS = ("WAVE", "AWAV", "FREQ", "VELO", "VRAD")


@dataclass
class KCWICube:
    """In-memory standardized representation of one KCWI/KCRM science cube.

    Parameters
    ----------
    path
        Source FITS file.
    arm
        Human-readable arm label, normally ``"BL"`` or ``"RH3"``.
    flux, uncertainty, variance
        Arrays in CRD_DAP order ``(y, x, wavelength)``.
    drp_mask
        Boolean DRP mask in standardized order.  ``True`` means unusable.
    flags
        Integer DRP bitmask in standardized order.  The production default is
        conservative: any non-zero flag is treated as unusable when building
        ``good``.  The raw values are preserved here for later inspection.
    noskysub
        Optional un-sky-subtracted cube in standardized order.
    wavelength
        One-dimensional wavelength axis in Angstrom.
    good
        Final hard-good sample mask after finite-value, DRP mask/flag,
        instrument-good wavelength, globally bad-channel, and spatial-spaxel
        requirements are combined.
    base_good
        Sample mask before channel/spaxel-level rejection.  Useful for QC.
    good_wavelength
        One-dimensional channel mask after ``WAVGOOD0/WAVGOOD1`` and bad-channel
        rejection.
    good_spaxel
        Two-dimensional spatial mask after the minimum-good-wavelength-fraction
        criterion.
    good_fraction_spaxel
        Fraction of instrument-good channels that survive ``base_good`` for each
        spatial element.
    bad_fraction_wavelength
        Fraction of spatial samples rejected in each wavelength channel before
        the global bad-channel cut is applied.
    header
        Copy of the primary FITS header.
    original_spectral_axis
        Numpy-axis index of the spectral dimension in the original FITS array.
    original_shape
        Shape of the original primary array.
    celestial_wcs
        Two-dimensional celestial WCS extracted from the FITS header.
    """

    path: Path
    arm: str
    flux: np.ndarray
    uncertainty: np.ndarray
    variance: np.ndarray
    drp_mask: np.ndarray
    flags: np.ndarray
    noskysub: np.ndarray | None
    wavelength: np.ndarray
    good: np.ndarray
    base_good: np.ndarray
    good_wavelength: np.ndarray
    good_spaxel: np.ndarray
    good_fraction_spaxel: np.ndarray
    bad_fraction_wavelength: np.ndarray
    header: fits.Header
    original_spectral_axis: int
    original_shape: tuple[int, ...]
    celestial_wcs: WCS

    @property
    def shape(self) -> tuple[int, int, int]:
        """Return standardized ``(ny, nx, nwave)`` shape."""
        return tuple(int(v) for v in self.flux.shape)

    @property
    def ny(self) -> int:
        return int(self.flux.shape[0])

    @property
    def nx(self) -> int:
        return int(self.flux.shape[1])

    @property
    def nwave(self) -> int:
        return int(self.flux.shape[2])

    def pixel_scales_arcsec(self) -> tuple[float, float]:
        """Return approximate celestial pixel scales ``(x, y)`` in arcsec/pixel."""
        scales = np.asarray(proj_plane_pixel_scales(self.celestial_wcs), dtype=float)
        if scales.size != 2 or not np.all(np.isfinite(scales)):
            raise ValueError(f"Could not determine two spatial pixel scales for {self.path}")
        return float(abs(scales[0]) * 3600.0), float(abs(scales[1]) * 3600.0)


def write_json(payload: dict[str, Any], path: str | Path) -> Path:
    """Write an indented JSON metadata file, creating parent directories."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n")
    return path


def read_json(path: str | Path) -> dict[str, Any]:
    """Read a JSON metadata file."""
    return json.loads(Path(path).read_text())


def inspect_fits_extensions(path: str | Path) -> list[dict[str, Any]]:
    """Return a compact machine-readable description of a FITS file's HDUs."""
    path = Path(path).expanduser().resolve()
    result: list[dict[str, Any]] = []
    with fits.open(path, memmap=True) as hdul:
        for idx, hdu in enumerate(hdul):
            data = getattr(hdu, "data", None)
            result.append(
                {
                    "index": idx,
                    "name": str(hdu.name),
                    "shape": None if data is None else list(data.shape),
                    "dtype": None if data is None else str(data.dtype),
                }
            )
    return result


def summarize_integer_flags(flags: np.ndarray) -> dict[int, int]:
    """Return exact counts of integer DRP flag values for provenance/QC logs.

    This summary intentionally reports *values* rather than attempting to assign
    undocumented scientific meanings to every bit.  The raw flag cube remains
    available, and CRD_DAP's default hard mask conservatively rejects non-zero
    values unless the target config explicitly chooses otherwise.
    """
    values, counts = np.unique(np.asarray(flags), return_counts=True)
    return {int(v): int(c) for v, c in zip(values, counts)}


def _find_spectral_fits_axis(header: fits.Header, ndim: int) -> int:
    """Return the one-based FITS axis number corresponding to wavelength."""
    for fits_axis in range(1, ndim + 1):
        ctype = str(header.get(f"CTYPE{fits_axis}", "")).upper()
        if any(token in ctype for token in _SPECTRAL_CTYPE_TOKENS):
            return fits_axis

    # KCWI DRP cubes are conventionally 3-D with wavelength on FITS axis 3.
    # We retain this as an explicit fallback only when the header lacks CTYPE.
    if ndim == 3 and ("CD3_3" in header or "CDELT3" in header):
        return 3

    raise ValueError(
        "Could not identify the spectral FITS axis from CTYPE keywords. "
        "CRD_DAP requires an explicit wavelength WCS."
    )


def _numpy_axis_from_fits_axis(fits_axis: int, ndim: int) -> int:
    """Convert a one-based FITS axis number to a zero-based numpy axis."""
    if not 1 <= fits_axis <= ndim:
        raise ValueError("FITS axis lies outside array dimensionality")
    return ndim - fits_axis


def _spectral_step_from_header(header: fits.Header, fits_axis: int) -> float:
    """Return the linear spectral increment in the header's native units."""
    cd_key = f"CD{fits_axis}_{fits_axis}"
    if cd_key in header:
        return float(header[cd_key])
    cdelt_key = f"CDELT{fits_axis}"
    if cdelt_key in header:
        return float(header[cdelt_key])
    raise ValueError(
        f"No {cd_key} or {cdelt_key} keyword found for the spectral axis; "
        "a linear wavelength solution is required by the current Script-1 loader."
    )


def wavelength_axis_from_header(
    header: fits.Header,
    n_pixels: int,
    *,
    fits_axis: int | None = None,
) -> np.ndarray:
    """Construct a linear wavelength axis in Angstrom from FITS WCS keywords.

    KCWI ``icubew``/``icubes`` products use a linear wavelength axis.  FITS pixel
    coordinates are one-based, so pixel index zero corresponds to FITS pixel 1.
    """
    if fits_axis is None:
        ndim = int(header.get("NAXIS", 3))
        fits_axis = _find_spectral_fits_axis(header, ndim)

    crval = float(header.get(f"CRVAL{fits_axis}", np.nan))
    crpix = float(header.get(f"CRPIX{fits_axis}", 1.0))
    if not np.isfinite(crval):
        raise ValueError(f"Missing/invalid CRVAL{fits_axis} for wavelength axis")
    step = _spectral_step_from_header(header, fits_axis)

    pix_fits = np.arange(n_pixels, dtype=float) + 1.0
    values = crval + (pix_fits - crpix) * step

    unit_text = str(header.get(f"CUNIT{fits_axis}", "Angstrom")).strip()
    try:
        unit = u.Unit(unit_text) if unit_text else u.AA
        values = (values * unit).to_value(u.AA)
    except Exception as exc:
        raise ValueError(
            f"Could not interpret spectral unit CUNIT{fits_axis}={unit_text!r}"
        ) from exc

    return np.asarray(values, dtype=float)


def _canonicalize(array: np.ndarray | None, spectral_numpy_axis: int) -> np.ndarray | None:
    """Move the original spectral numpy axis to the final position."""
    if array is None:
        return None
    arr = np.asarray(array)
    if arr.ndim != 3:
        raise ValueError(f"Expected a 3-D KCWI cube extension, received shape {arr.shape}")
    return np.moveaxis(arr, spectral_numpy_axis, -1)


def _uncanonicalize(array: np.ndarray, spectral_numpy_axis: int) -> np.ndarray:
    """Move standardized final wavelength axis back to the original numpy axis."""
    return np.moveaxis(np.asarray(array), -1, spectral_numpy_axis)


def _read_optional_hdu(hdul: fits.HDUList, name: str) -> np.ndarray | None:
    try:
        data = hdul[name].data
    except (KeyError, IndexError):
        return None
    if data is None:
        return None
    return np.asarray(data)


def build_quality_masks(
    flux: np.ndarray,
    uncertainty: np.ndarray,
    drp_mask: np.ndarray | None,
    flags: np.ndarray | None,
    wavelength: np.ndarray,
    *,
    wavegood0: float | None,
    wavegood1: float | None,
    reject_any_nonzero_flag: bool,
    min_good_wavelength_fraction: float,
    bad_channel_fraction_threshold: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Construct Script-1 hard-good masks and QC fractions.

    Returns
    -------
    good
        Final 3-D hard-good mask.
    base_good
        3-D mask before global channel and spatial-spaxel rejection.
    good_wavelength
        One-dimensional instrument/global channel mask.
    good_spaxel
        Two-dimensional spatial mask.
    bad_fraction_wavelength
        Fraction of spatial samples rejected at each wavelength before the
        global channel cut.
    """
    f = np.asarray(flux, dtype=float)
    s = np.asarray(uncertainty, dtype=float)
    if f.shape != s.shape or f.ndim != 3:
        raise ValueError("flux and uncertainty must be matching 3-D standardized cubes")

    base_good = np.isfinite(f) & np.isfinite(s) & (s > 0)

    if drp_mask is not None:
        m = np.asarray(drp_mask)
        if m.shape != f.shape:
            raise ValueError("DRP mask shape does not match flux cube")
        base_good &= ~(m.astype(bool))

    if flags is not None and reject_any_nonzero_flag:
        flg = np.asarray(flags)
        if flg.shape != f.shape:
            raise ValueError("FLAGS shape does not match flux cube")
        base_good &= flg == 0

    wave = np.asarray(wavelength, dtype=float)
    if wave.ndim != 1 or wave.size != f.shape[-1]:
        raise ValueError("wavelength axis does not match cube spectral dimension")

    instrument_good = np.isfinite(wave)
    if wavegood0 is not None and np.isfinite(wavegood0):
        instrument_good &= wave >= float(wavegood0)
    if wavegood1 is not None and np.isfinite(wavegood1):
        instrument_good &= wave <= float(wavegood1)

    # Bad-channel fraction is measured across all spatial elements using only
    # the sample-level data-quality mask.  Channels outside WAVGOOD are forced
    # to bad separately and therefore need not influence this diagnostic.
    n_spatial = f.shape[0] * f.shape[1]
    bad_fraction_wavelength = 1.0 - np.sum(base_good, axis=(0, 1)) / float(n_spatial)
    global_channel_good = bad_fraction_wavelength <= float(bad_channel_fraction_threshold)
    good_wavelength = instrument_good & global_channel_good

    if not np.any(instrument_good):
        raise ValueError("No wavelength samples survive WAVGOOD/instrument-good constraints")

    denom = int(np.sum(instrument_good))
    good_fraction_spaxel = np.sum(base_good[..., instrument_good], axis=-1) / float(denom)
    good_spaxel = good_fraction_spaxel >= float(min_good_wavelength_fraction)

    good = base_good & good_wavelength[None, None, :] & good_spaxel[..., None]
    return good, base_good, good_wavelength, good_spaxel, bad_fraction_wavelength


def load_kcwi_cube(
    path: str | Path,
    *,
    arm: str,
    reject_any_nonzero_flag: bool = True,
    min_good_wavelength_fraction: float = 0.80,
    bad_channel_fraction_threshold: float = 0.50,
) -> KCWICube:
    """Load a KCWI/KCRM cube and construct the standardized Script-1 data model.

    The loader requires ``UNCERT`` because the profile-likelihood analysis cannot
    proceed with an undefined noise model.  ``MASK``, ``FLAGS``, and ``NOSKYSUB``
    are read when present.  Missing ``MASK``/``FLAGS`` are represented by clean
    arrays, but their absence should be visible in Script-1 logs.
    """
    source = Path(path).expanduser().resolve()
    if not source.exists():
        raise FileNotFoundError(source)

    with fits.open(source, memmap=False) as hdul:
        if hdul[0].data is None:
            raise ValueError(f"PRIMARY contains no data: {source}")
        primary = np.asarray(hdul[0].data)
        header = hdul[0].header.copy()
        if primary.ndim != 3:
            raise ValueError(f"KCWI science cube must be 3-D; got {primary.shape} from {source}")

        uncert = _read_optional_hdu(hdul, "UNCERT")
        if uncert is None:
            raise ValueError(
                f"Required UNCERT extension is missing from {source}. "
                "CRD_DAP expects the KCWI DRP 1-sigma uncertainty product."
            )
        mask = _read_optional_hdu(hdul, "MASK")
        flags = _read_optional_hdu(hdul, "FLAGS")
        noskysub = _read_optional_hdu(hdul, "NOSKYSUB")

    fits_spec_axis = _find_spectral_fits_axis(header, primary.ndim)
    np_spec_axis = _numpy_axis_from_fits_axis(fits_spec_axis, primary.ndim)

    flux_c = np.asarray(_canonicalize(primary, np_spec_axis), dtype=float)
    uncert_c = np.asarray(_canonicalize(uncert, np_spec_axis), dtype=float)
    mask_c = _canonicalize(mask, np_spec_axis)
    flags_c = _canonicalize(flags, np_spec_axis)
    nosky_c = _canonicalize(noskysub, np_spec_axis)

    if flux_c.shape != uncert_c.shape:
        raise ValueError("PRIMARY and UNCERT shapes differ after axis standardization")

    if mask_c is None:
        mask_c = np.zeros(flux_c.shape, dtype=bool)
    else:
        mask_c = np.asarray(mask_c).astype(bool)

    if flags_c is None:
        flags_c = np.zeros(flux_c.shape, dtype=np.uint16)
    else:
        flags_c = np.asarray(flags_c)

    wave = wavelength_axis_from_header(header, flux_c.shape[-1], fits_axis=fits_spec_axis)
    wavegood0 = header.get("WAVGOOD0")
    wavegood1 = header.get("WAVGOOD1")

    good, base_good, good_wave, good_spaxel, bad_wave_frac = build_quality_masks(
        flux_c,
        uncert_c,
        mask_c,
        flags_c,
        wave,
        wavegood0=wavegood0,
        wavegood1=wavegood1,
        reject_any_nonzero_flag=reject_any_nonzero_flag,
        min_good_wavelength_fraction=min_good_wavelength_fraction,
        bad_channel_fraction_threshold=bad_channel_fraction_threshold,
    )

    try:
        celestial = WCS(header).celestial
        if celestial.pixel_n_dim != 2 or celestial.world_n_dim != 2:
            raise ValueError
    except Exception as exc:
        raise ValueError(
            f"A valid 2-D celestial WCS is required for BL/RH3 registration: {source}"
        ) from exc

    # Recompute good-fraction array here so it is directly available on the
    # object (build_quality_masks returns only the final boolean spatial mask).
    instrument_good = np.isfinite(wave)
    if wavegood0 is not None:
        instrument_good &= wave >= float(wavegood0)
    if wavegood1 is not None:
        instrument_good &= wave <= float(wavegood1)
    good_fraction_spaxel = np.sum(base_good[..., instrument_good], axis=-1) / float(
        np.sum(instrument_good)
    )

    return KCWICube(
        path=source,
        arm=str(arm),
        flux=flux_c,
        uncertainty=uncert_c,
        variance=uncert_c**2,
        drp_mask=mask_c,
        flags=flags_c,
        noskysub=None if nosky_c is None else np.asarray(nosky_c, dtype=float),
        wavelength=wave,
        good=good,
        base_good=base_good,
        good_wavelength=good_wave,
        good_spaxel=good_spaxel,
        good_fraction_spaxel=good_fraction_spaxel,
        bad_fraction_wavelength=bad_wave_frac,
        header=header,
        original_spectral_axis=np_spec_axis,
        original_shape=tuple(int(v) for v in primary.shape),
        celestial_wcs=celestial,
    )


def save_prepared_cube(cube: KCWICube, destination: str | Path, *, overwrite: bool = False) -> Path:
    """Write a standardized Script-1 FITS product while preserving input WCS.

    The primary data and native DRP extensions are written back in their original
    axis order.  CRD_DAP-specific quality products are appended as additional
    extensions.  Downstream scripts can therefore use the same WCS convention as
    the original DRP cube while gaining explicit hard-good masks.
    """
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)

    header = cube.header.copy()
    header["CRDDAP"] = (True, "Processed by CRD_DAP")
    header["CRDSTEP"] = (1, "CRD_DAP preparation step")
    header["CRDARM"] = (cube.arm, "CRD_DAP arm label")

    primary = fits.PrimaryHDU(
        _uncanonicalize(cube.flux, cube.original_spectral_axis).astype(np.float32),
        header=header,
    )
    hdus: list[fits.hdu.base.ExtensionHDU | fits.PrimaryHDU] = [primary]

    hdus.append(
        fits.ImageHDU(
            _uncanonicalize(cube.uncertainty, cube.original_spectral_axis).astype(np.float32),
            name="UNCERT",
        )
    )
    hdus.append(
        fits.ImageHDU(
            _uncanonicalize(cube.drp_mask.astype(np.uint8), cube.original_spectral_axis),
            name="MASK",
        )
    )
    hdus.append(
        fits.ImageHDU(
            _uncanonicalize(cube.flags, cube.original_spectral_axis),
            name="FLAGS",
        )
    )
    if cube.noskysub is not None:
        hdus.append(
            fits.ImageHDU(
                _uncanonicalize(cube.noskysub, cube.original_spectral_axis).astype(np.float32),
                name="NOSKYSUB",
            )
        )

    hdus.append(
        fits.ImageHDU(
            _uncanonicalize(cube.good.astype(np.uint8), cube.original_spectral_axis),
            name="GOODMASK",
        )
    )
    hdus.append(fits.ImageHDU(cube.good_spaxel.astype(np.uint8), name="GOODSPAX"))
    hdus.append(fits.ImageHDU(cube.good_fraction_spaxel.astype(np.float32), name="GOODFRAC"))
    hdus.append(fits.ImageHDU(cube.good_wavelength.astype(np.uint8), name="GOODWAVE"))
    hdus.append(fits.ImageHDU(cube.wavelength.astype(np.float64), name="WAVELENGTH"))

    fits.HDUList(hdus).writeto(destination, overwrite=overwrite)
    return destination


def _header_match_token(header: fits.Header, key: str) -> str:
    """Return a normalized string token used for calibration-file matching."""
    return str(header.get(key, "")).strip().upper()


def _sidecar_header_matches_arc(sidecar_header: fits.Header, arc_header: fits.Header) -> bool:
    """Return True when a geometry product plausibly belongs to ``arc_header``.

    Real KCWI reductions do not always preserve the master-arc exposure number
    in the *geometry-map filename*.  In the 2024-12-26 red-side test data used to
    validate CRD_DAP, for example, the master arc is named ``...00025_marc.fits``
    while the maps are named ``...00027_*map.fits``; their FITS ``OFNAME`` values
    still identify the same original ``...00025.fits`` calibration.  Header
    provenance is therefore a more reliable fallback than filename roots alone.
    """
    arc_ofname = _header_match_token(arc_header, "OFNAME")
    side_ofname = _header_match_token(sidecar_header, "OFNAME")
    if arc_ofname and side_ofname and arc_ofname != side_ofname:
        return False

    # These setup fields should agree for maps generated from the same geometry
    # solution.  Missing values are tolerated here; the science-vs-arc check is
    # stricter and occurs separately in validation.py.
    for key in ("CAMERA", "IFUNAM", "BINNING", "BGRATNAM", "RGRATNAM"):
        a = _header_match_token(arc_header, key)
        b = _header_match_token(sidecar_header, key)
        if a and b and a != b:
            return False
    return True


def discover_arc_sidecars(master_arc: str | Path) -> dict[str, Path]:
    """Discover KCWI DRP wavelength/slice/position maps for a master arc.

    Discovery proceeds in three increasingly permissive stages:

    1. the conventional identical filename root;
    2. a unique same-root glob match;
    3. a header-provenance search for maps whose ``OFNAME`` and instrumental
       configuration match the master arc.

    The third stage is required for real red-side reductions in which map
    filenames can carry a different exposure number from the associated master
    arc even though their FITS headers trace them to the same input arc.  Users
    can always bypass discovery with explicit paths in the target config.
    """
    arc = Path(master_arc).expanduser().resolve()
    stem = arc.stem
    lower = stem.lower()
    for suffix in ("_marcs", "_marc"):
        if lower.endswith(suffix):
            root = stem[: -len(suffix)]
            break
    else:
        root = stem

    arc_header = fits.getheader(arc, 0)
    result: dict[str, Path] = {}
    for kind in ("wavemap", "slicemap", "posmap"):
        exact = arc.with_name(f"{root}_{kind}.fits")
        if exact.exists():
            result[kind] = exact.resolve()
            continue

        same_root = sorted(arc.parent.glob(f"{root}*_{kind}.fits"))
        if len(same_root) == 1:
            result[kind] = same_root[0].resolve()
            continue
        if len(same_root) > 1:
            header_matches = []
            for candidate in same_root:
                try:
                    if _sidecar_header_matches_arc(fits.getheader(candidate, 0), arc_header):
                        header_matches.append(candidate.resolve())
                except Exception:
                    continue
            if len(header_matches) == 1:
                result[kind] = header_matches[0]
                continue

        # Filename-root discovery failed.  Search all nearby maps of this type
        # and use FITS provenance/setup metadata.  This is intentionally local to
        # the arc directory so it cannot wander through unrelated reductions.
        candidates = sorted(arc.parent.glob(f"*_{kind}.fits"))
        header_matches = []
        for candidate in candidates:
            try:
                chead = fits.getheader(candidate, 0)
            except Exception:
                continue
            if _sidecar_header_matches_arc(chead, arc_header):
                header_matches.append(candidate.resolve())

        if len(header_matches) == 1:
            result[kind] = header_matches[0]
        elif len(header_matches) == 0:
            raise FileNotFoundError(
                f"Could not discover {kind} sidecar for master arc {arc}. "
                "No filename-root or FITS-header provenance match was found. "
                "Supply an explicit config path."
            )
        else:
            raise RuntimeError(
                f"Multiple header-matched {kind} sidecars found for {arc}: "
                + ", ".join(str(m) for m in header_matches)
                + ". Supply an explicit config path to remove the ambiguity."
            )

    return result
