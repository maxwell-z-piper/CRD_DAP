"""KCWI/KCRM FITS I/O and standardized prepared-cube products.

This module deliberately centralizes all assumptions about the KCWI DRP FITS
layout.  Downstream science code should work with :class:`KCWICube` rather than
repeating extension names, FITS-axis conventions, or mask semantics.

The production implementation supports two science-input layouts:

1. a KCWI DRP multi-extension cube whose science image is stored in ``PRIMARY``
   and which contains ``UNCERT`` plus optional ``MASK``, ``FLAGS``, and
   ``NOSKYSUB`` extensions; and
2. a KcwiKit post-DRP stack represented by four matched files: ``*_icubes``
   (flux), ``*_vcubes`` (variance), ``*_mcubes`` (final binary stack mask), and
   ``*_ecubes`` (effective exposure time).

For KcwiKit stacks the final ``mcubes`` product is treated as a binary validity /
coverage mask: zero means that at least one unmasked input exposure contributed,
and non-zero means no valid contribution.  Original KCWI DRP bit-level FLAGS are
used by KcwiKit while stacking but are not recoverable from the final binary
stack mask, so CRD_DAP never pretends otherwise.

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
    exposure: np.ndarray | None
    drp_mask: np.ndarray
    flags: np.ndarray | None
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
    input_format: str = "drp"
    source_paths: dict[str, Path] | None = None
    coverage_fraction_spaxel: np.ndarray | None = None

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
    exposure: np.ndarray | None = None,
    wavegood0: float | None,
    wavegood1: float | None,
    reject_any_nonzero_flag: bool,
    min_good_wavelength_fraction: float,
    bad_channel_fraction_threshold: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Construct Script-1 hard-good masks and QC fractions.

    A subtle but important point for post-stacking cubes is that the requested
    KcwiKit output grid can be larger than the actually exposed IFU footprint.
    Those padded pixels must *not* count as bad spatial samples when deciding
    whether an entire wavelength channel is globally problematic.  Otherwise a
    deliberately generous output canvas could make every wavelength appear
    mostly bad.

    ``exposure`` therefore serves two roles when available:

    * a sample is geometrically covered only when its effective exposure is
      finite and strictly positive;
    * the global bad-channel fraction is evaluated only over the core spatial
      footprint that has coverage across the configured fraction of the
      instrument-good wavelength range.

    Returns
    -------
    good
        Final 3-D hard-good mask.
    base_good
        3-D sample mask before global channel and spatial-spaxel rejection.
    good_wavelength
        One-dimensional instrument/global channel mask.
    good_spaxel
        Two-dimensional spatial mask.
    bad_fraction_wavelength
        Fraction of covered/core spatial samples rejected at each wavelength.
    good_fraction_spaxel
        Fraction of instrument-good wavelengths that are fully usable per
        spatial element.
    coverage_fraction_spaxel
        Fraction of instrument-good wavelengths with positive geometric
        exposure/coverage per spatial element.
    """
    f = np.asarray(flux)
    s = np.asarray(uncertainty)
    if f.shape != s.shape or f.ndim != 3:
        raise ValueError("flux and uncertainty must be matching 3-D standardized cubes")

    wave = np.asarray(wavelength, dtype=float)
    if wave.ndim != 1 or wave.size != f.shape[-1]:
        raise ValueError("wavelength axis does not match cube spectral dimension")

    instrument_good = np.isfinite(wave)
    if wavegood0 is not None and np.isfinite(wavegood0):
        instrument_good &= wave >= float(wavegood0)
    if wavegood1 is not None and np.isfinite(wavegood1):
        instrument_good &= wave <= float(wavegood1)
    if not np.any(instrument_good):
        raise ValueError("No wavelength samples survive WAVGOOD/instrument-good constraints")

    # Geometric coverage is explicit for KcwiKit stacks.  For a native DRP cube
    # without an exposure cube, finite flux/noise samples are the best available
    # coverage proxy; quality masks are still applied separately below.
    if exposure is not None:
        exp = np.asarray(exposure)
        if exp.shape != f.shape:
            raise ValueError("Exposure cube shape does not match flux cube")
        coverage_sample = np.isfinite(exp) & (exp > 0)
    else:
        coverage_sample = np.isfinite(f) & np.isfinite(s) & (s > 0)

    base_good = coverage_sample & np.isfinite(f) & np.isfinite(s) & (s > 0)

    if drp_mask is not None:
        m = np.asarray(drp_mask)
        if m.shape != f.shape:
            raise ValueError("Mask shape does not match flux cube")
        base_good &= ~(m.astype(bool))

    if flags is not None and reject_any_nonzero_flag:
        flg = np.asarray(flags)
        if flg.shape != f.shape:
            raise ValueError("FLAGS shape does not match flux cube")
        base_good &= flg == 0

    denom = int(np.sum(instrument_good))
    coverage_fraction_spaxel = (
        np.sum(coverage_sample[..., instrument_good], axis=-1) / float(denom)
    )
    good_fraction_spaxel = (
        np.sum(base_good[..., instrument_good], axis=-1) / float(denom)
    )

    # Restrict wavelength-level QC to the real IFU footprint rather than the
    # zero-exposure output-grid padding.  The same minimum-fraction threshold is
    # intentionally used for the first implementation so this choice remains
    # transparent and configurable rather than introducing another hidden cut.
    coverage_core = coverage_fraction_spaxel >= float(min_good_wavelength_fraction)
    if not np.any(coverage_core):
        # A very short/test observation can have wavelength-dependent edge
        # coverage that fails the strict core definition. Fall back to any
        # spatial element with at least one covered instrument-good sample, but
        # keep the final good_spaxel criterion unchanged.
        coverage_core = coverage_fraction_spaxel > 0
    if not np.any(coverage_core):
        raise ValueError("No spatial samples have positive coverage in the instrument-good range")

    n_core = int(np.sum(coverage_core))
    bad_fraction_wavelength = 1.0 - (
        np.sum(base_good[coverage_core, :], axis=0) / float(n_core)
    )
    global_channel_good = bad_fraction_wavelength <= float(bad_channel_fraction_threshold)
    good_wavelength = instrument_good & global_channel_good

    good_spaxel = good_fraction_spaxel >= float(min_good_wavelength_fraction)
    good = base_good & good_wavelength[None, None, :] & good_spaxel[..., None]

    return (
        good,
        base_good,
        good_wavelength,
        good_spaxel,
        bad_fraction_wavelength,
        good_fraction_spaxel,
        coverage_fraction_spaxel,
    )


def _load_primary_array_and_header(
    path: str | Path,
    *,
    dtype: np.dtype | type | None = None,
) -> tuple[np.ndarray, fits.Header]:
    """Load one primary-HDU image cube and copy its header.

    KcwiKit final products are separate single-HDU files.  We make an explicit
    in-memory copy here so the returned array remains valid after the FITS file
    closes and so the configured float precision is under CRD_DAP control.
    """
    source = Path(path).expanduser().resolve()
    if not source.exists():
        raise FileNotFoundError(source)
    with fits.open(source, memmap=True) as hdul:
        if hdul[0].data is None:
            raise ValueError(f"PRIMARY contains no data: {source}")
        arr = np.array(hdul[0].data, dtype=dtype, copy=True)
        hdr = hdul[0].header.copy()
    if arr.ndim != 3:
        raise ValueError(f"Expected a 3-D primary cube in {source}; got shape {arr.shape}")
    return arr, hdr


def _wcs_signature(header: fits.Header) -> dict[str, Any]:
    """Return the spatial/spectral WCS values that must match companion cubes."""
    keys = [
        "NAXIS1", "NAXIS2", "NAXIS3",
        "CTYPE1", "CTYPE2", "CTYPE3",
        "CUNIT1", "CUNIT2", "CUNIT3",
        "CRVAL1", "CRVAL2", "CRVAL3",
        "CRPIX1", "CRPIX2", "CRPIX3",
        "CD1_1", "CD1_2", "CD1_3",
        "CD2_1", "CD2_2", "CD2_3",
        "CD3_1", "CD3_2", "CD3_3",
        "CDELT1", "CDELT2", "CDELT3",
    ]
    return {key: header[key] for key in keys if key in header}


def _headers_have_same_wcs(reference: fits.Header, other: fits.Header) -> tuple[bool, list[str]]:
    """Compare the WCS-bearing header values of two KcwiKit companion cubes."""
    ref = _wcs_signature(reference)
    oth = _wcs_signature(other)
    mismatches: list[str] = []
    for key in sorted(set(ref) | set(oth)):
        if key not in ref or key not in oth:
            mismatches.append(f"{key}: missing from one companion header")
            continue
        a, b = ref[key], oth[key]
        try:
            af, bf = float(a), float(b)
            if not np.isclose(af, bf, rtol=1e-10, atol=1e-12, equal_nan=True):
                mismatches.append(f"{key}: {a!r} != {b!r}")
        except (TypeError, ValueError):
            if str(a).strip().upper() != str(b).strip().upper():
                mismatches.append(f"{key}: {a!r} != {b!r}")
    return len(mismatches) == 0, mismatches


def load_kcwikit_stack(
    icube_path: str | Path,
    vcube_path: str | Path,
    mcube_path: str | Path,
    ecube_path: str | Path,
    *,
    arm: str,
    min_good_wavelength_fraction: float = 0.80,
    bad_channel_fraction_threshold: float = 0.50,
    float_dtype: str | np.dtype = "float32",
) -> KCWICube:
    """Load a four-file KcwiKit post-DRP stack into the CRD_DAP cube model.

    Parameters
    ----------
    icube_path, vcube_path, mcube_path, ecube_path
        KcwiKit stacked flux, variance, binary mask, and effective-exposure
        products, respectively.
    arm
        CRD_DAP arm label (normally ``BL`` or ``RH3``).
    float_dtype
        In-memory precision for the large stacked flux/variance/exposure cubes.
        ``float32`` is the default to keep a BL+RH3 Script-1 run practical on a
        workstation; extracted spectra may later be promoted to float64 for
        pPXF.  Set to ``float64`` if memory is plentiful and exact preservation
        of the KcwiKit storage precision is desired.

    Notes
    -----
    KcwiKit's current stacking implementation uses the input PyDRP ``FLAGS`` to
    reject bad samples, but its final ``mcubes`` product is set from whether the
    effective exposure is zero.  Consequently the original bit-level flag values
    are not present in the final stack and CRD_DAP stores a zero-valued placeholder
    ``flags`` array while explicitly recording ``input_format='kcwikit'``.
    """
    paths = {
        "icube": Path(icube_path).expanduser().resolve(),
        "vcube": Path(vcube_path).expanduser().resolve(),
        "mcube": Path(mcube_path).expanduser().resolve(),
        "ecube": Path(ecube_path).expanduser().resolve(),
    }
    for role, source in paths.items():
        if not source.exists():
            raise FileNotFoundError(f"KcwiKit {role} does not exist: {source}")

    dtype = np.dtype(float_dtype)
    primary, header = _load_primary_array_and_header(paths["icube"], dtype=dtype)
    variance_native, vheader = _load_primary_array_and_header(paths["vcube"], dtype=dtype)
    mask_native, mheader = _load_primary_array_and_header(paths["mcube"], dtype=np.int16)
    exposure_native, eheader = _load_primary_array_and_header(paths["ecube"], dtype=dtype)

    for role, arr in (
        ("vcube", variance_native),
        ("mcube", mask_native),
        ("ecube", exposure_native),
    ):
        if arr.shape != primary.shape:
            raise ValueError(
                f"KcwiKit {role} shape {arr.shape} does not match icube shape {primary.shape}"
            )

    for role, hdr in (("vcube", vheader), ("mcube", mheader), ("ecube", eheader)):
        same, mismatches = _headers_have_same_wcs(header, hdr)
        if not same:
            details = "; ".join(mismatches[:12])
            raise ValueError(
                f"KcwiKit {role} WCS does not match icube WCS. "
                f"First mismatch(es): {details}"
            )

    fits_spec_axis = _find_spectral_fits_axis(header, primary.ndim)
    np_spec_axis = _numpy_axis_from_fits_axis(fits_spec_axis, primary.ndim)

    flux_c = np.asarray(_canonicalize(primary, np_spec_axis), dtype=dtype)
    var_c = np.asarray(_canonicalize(variance_native, np_spec_axis), dtype=dtype)
    stack_mask_c = np.asarray(_canonicalize(mask_native, np_spec_axis), dtype=np.int16)
    exposure_c = np.asarray(_canonicalize(exposure_native, np_spec_axis), dtype=dtype)

    # The current KcwiKit final mask is binary.  Fail loudly if a future version
    # starts writing another convention so we can inspect it rather than silently
    # interpreting bit values incorrectly.
    unexpected_mask = (stack_mask_c != 0) & (stack_mask_c != 1)
    if np.any(unexpected_mask):
        sample_values = np.asarray(stack_mask_c[unexpected_mask][:20])
        raise ValueError(
            "KcwiKit final mcube contains values outside {0,1}; its mask semantics "
            f"may have changed. Example unexpected values: {sample_values!r}."
        )
    mask_c = stack_mask_c != 0

    # KcwiKit vcubes are variances, not 1-sigma uncertainties.  Only positive,
    # finite variance values are square-rooted. Invalid/padded samples remain NaN
    # here and are subsequently excluded by the hard-good mask.
    uncertainty_c = np.full(var_c.shape, np.nan, dtype=dtype)
    positive_var = np.isfinite(var_c) & (var_c > 0)
    uncertainty_c[positive_var] = np.sqrt(var_c[positive_var]).astype(dtype, copy=False)

    covered = np.isfinite(exposure_c) & (exposure_c > 0) & (~mask_c)
    bad_negative = covered & np.isfinite(var_c) & (var_c < 0)
    if np.any(bad_negative):
        raise ValueError(
            f"KcwiKit variance cube contains {int(np.sum(bad_negative))} negative "
            "variance samples inside valid exposure coverage."
        )

    # Original PyDRP flags are no longer available after KcwiKit's final binary
    # mask construction.  Do not infer or fabricate bit meanings.
    flags_c = None

    wave = wavelength_axis_from_header(header, flux_c.shape[-1], fits_axis=fits_spec_axis)
    wavegood0 = header.get("WAVGOOD0")
    wavegood1 = header.get("WAVGOOD1")

    (
        good,
        base_good,
        good_wave,
        good_spaxel,
        bad_wave_frac,
        good_fraction_spaxel,
        coverage_fraction_spaxel,
    ) = build_quality_masks(
        flux_c,
        uncertainty_c,
        mask_c,
        flags_c,
        wave,
        exposure=exposure_c,
        wavegood0=wavegood0,
        wavegood1=wavegood1,
        reject_any_nonzero_flag=False,
        min_good_wavelength_fraction=min_good_wavelength_fraction,
        bad_channel_fraction_threshold=bad_channel_fraction_threshold,
    )

    try:
        celestial = WCS(header).celestial
        if celestial.pixel_n_dim != 2 or celestial.world_n_dim != 2:
            raise ValueError
    except Exception as exc:
        raise ValueError(
            f"A valid 2-D celestial WCS is required for BL/RH3 registration: {paths['icube']}"
        ) from exc

    return KCWICube(
        path=paths["icube"],
        arm=str(arm),
        flux=flux_c,
        uncertainty=uncertainty_c,
        variance=var_c,
        exposure=exposure_c,
        drp_mask=mask_c,
        flags=flags_c,
        noskysub=None,
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
        input_format="kcwikit",
        source_paths=paths,
        coverage_fraction_spaxel=coverage_fraction_spaxel,
    )

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
        flags_c = None
    else:
        flags_c = np.asarray(flags_c)

    wave = wavelength_axis_from_header(header, flux_c.shape[-1], fits_axis=fits_spec_axis)
    wavegood0 = header.get("WAVGOOD0")
    wavegood1 = header.get("WAVGOOD1")

    (
        good,
        base_good,
        good_wave,
        good_spaxel,
        bad_wave_frac,
        good_fraction_spaxel,
        coverage_fraction_spaxel,
    ) = build_quality_masks(
        flux_c,
        uncert_c,
        mask_c,
        flags_c,
        wave,
        exposure=None,
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

    return KCWICube(
        path=source,
        arm=str(arm),
        flux=flux_c,
        uncertainty=uncert_c,
        variance=uncert_c**2,
        exposure=None,
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
        input_format="drp",
        source_paths={"cube": source},
        coverage_fraction_spaxel=coverage_fraction_spaxel,
    )



def load_prepared_cube(path: str | Path, *, expected_arm: str | None = None) -> KCWICube:
    """Load a standardized Script-1 ``prepared_*.fits`` product.

    Downstream stages must use the exact hard-good sample/spaxel/wavelength masks
    written by Script 1 rather than reconstructing them from the original DRP
    extensions.  This loader therefore requires ``GOODMASK``, ``GOODSPAX``,
    ``GOODWAVE``, and ``WAVELENGTH`` and preserves the original celestial WCS.

    Parameters
    ----------
    path
        Script-1 prepared FITS file.
    expected_arm
        Optional arm label (normally ``"BL"`` or ``"RH3"``) that is checked
        against the ``CRDARM`` provenance keyword.
    """
    source = Path(path).expanduser().resolve()
    if not source.exists():
        raise FileNotFoundError(f"Prepared cube does not exist: {source}")

    with fits.open(source, memmap=True) as hdul:
        primary = np.asarray(hdul[0].data)
        header = hdul[0].header.copy()
        if primary.ndim != 3:
            raise ValueError(f"Prepared KCWI cube must be 3-D; got {primary.shape} from {source}")
        if int(header.get("CRDSTEP", -1)) != 1:
            raise ValueError(
                f"{source} is not marked as a CRD_DAP Script-1 prepared cube (CRDSTEP=1 required)"
            )

        arm = str(header.get("CRDARM", "")).strip()
        if expected_arm is not None and arm.upper() != str(expected_arm).upper():
            raise ValueError(
                f"Prepared cube {source} reports CRDARM={arm!r}, expected {expected_arm!r}"
            )

        required = ("UNCERT", "GOODMASK", "GOODSPAX", "GOODWAVE", "WAVELENGTH")
        missing = [name for name in required if name not in hdul]
        if missing:
            raise ValueError(
                f"Prepared cube {source} is missing required extension(s): {', '.join(missing)}"
            )

        uncert_native = np.asarray(hdul["UNCERT"].data)
        good_native = np.asarray(hdul["GOODMASK"].data)
        good_spaxel = np.asarray(hdul["GOODSPAX"].data).astype(bool)
        good_fraction = np.asarray(hdul["GOODFRAC"].data, dtype=float) if "GOODFRAC" in hdul else good_spaxel.astype(float)
        coverage_fraction = (
            np.asarray(hdul["COVFRAC"].data, dtype=float) if "COVFRAC" in hdul else None
        )
        good_wave = np.asarray(hdul["GOODWAVE"].data).astype(bool).ravel()
        wavelength = np.asarray(hdul["WAVELENGTH"].data, dtype=float).ravel()
        mask_native = np.asarray(hdul["MASK"].data) if "MASK" in hdul else np.zeros_like(primary, dtype=np.uint8)
        flags_native = np.asarray(hdul["FLAGS"].data) if "FLAGS" in hdul else None
        exposure_native = np.asarray(hdul["EXPOSURE"].data) if "EXPOSURE" in hdul else None
        nosky_native = np.asarray(hdul["NOSKYSUB"].data) if "NOSKYSUB" in hdul else None

    fits_spec_axis = _find_spectral_fits_axis(header, primary.ndim)
    np_spec_axis = _numpy_axis_from_fits_axis(fits_spec_axis, primary.ndim)
    flux = np.asarray(_canonicalize(primary, np_spec_axis), dtype=float)
    uncertainty = np.asarray(_canonicalize(uncert_native, np_spec_axis), dtype=float)
    good = np.asarray(_canonicalize(good_native, np_spec_axis), dtype=bool)
    drp_mask = np.asarray(_canonicalize(mask_native, np_spec_axis), dtype=bool)
    flags = None if flags_native is None else _canonicalize(flags_native, np_spec_axis)
    exposure = None if exposure_native is None else np.asarray(_canonicalize(exposure_native, np_spec_axis), dtype=float)
    noskysub = None if nosky_native is None else np.asarray(_canonicalize(nosky_native, np_spec_axis), dtype=float)

    if flux.shape != uncertainty.shape or flux.shape != good.shape:
        raise ValueError("Prepared PRIMARY/UNCERT/GOODMASK shapes differ after axis standardization")
    if wavelength.size != flux.shape[-1] or good_wave.size != flux.shape[-1]:
        raise ValueError("Prepared WAVELENGTH/GOODWAVE length does not match spectral axis")
    if good_spaxel.shape != flux.shape[:2]:
        raise ValueError("Prepared GOODSPAX shape does not match spatial cube shape")

    variance = uncertainty**2
    # Script 1 does not save its pre-channel/pre-spaxel ``base_good`` mask because
    # downstream stages only need the final immutable hard-good mask.  Reuse the
    # final mask here rather than pretending to reconstruct an earlier state.
    base_good = good.copy()
    bad_fraction = np.full(wavelength.shape, np.nan, dtype=float)
    for j in range(wavelength.size):
        spatial = good_spaxel
        if np.any(spatial):
            bad_fraction[j] = 1.0 - float(np.mean(good[..., j][spatial]))

    try:
        celestial = WCS(header).celestial
        if celestial.pixel_n_dim != 2 or celestial.world_n_dim != 2:
            raise ValueError
    except Exception as exc:
        raise ValueError(f"Prepared cube lacks a valid 2-D celestial WCS: {source}") from exc

    return KCWICube(
        path=source,
        arm=arm or str(expected_arm or ""),
        flux=flux,
        uncertainty=uncertainty,
        variance=variance,
        exposure=exposure,
        drp_mask=drp_mask,
        flags=None if flags is None else np.asarray(flags),
        noskysub=noskysub,
        wavelength=wavelength,
        good=good,
        base_good=base_good,
        good_wavelength=good_wave,
        good_spaxel=good_spaxel,
        good_fraction_spaxel=good_fraction,
        bad_fraction_wavelength=bad_fraction,
        header=header,
        original_spectral_axis=np_spec_axis,
        original_shape=tuple(int(v) for v in primary.shape),
        celestial_wcs=celestial,
        input_format="prepared",
        source_paths={"prepared_cube": source},
        coverage_fraction_spaxel=coverage_fraction,
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
    header["CRDINFMT"] = (cube.input_format, "CRD_DAP science input format")

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
    if cube.exposure is not None:
        exp_hdu = fits.ImageHDU(
            _uncanonicalize(cube.exposure, cube.original_spectral_axis).astype(np.float32),
            name="EXPOSURE",
        )
        exp_hdu.header["BUNIT"] = "s"
        hdus.append(exp_hdu)
    if cube.flags is not None:
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
    if cube.coverage_fraction_spaxel is not None:
        hdus.append(
            fits.ImageHDU(cube.coverage_fraction_spaxel.astype(np.float32), name="COVFRAC")
        )
    hdus.append(fits.ImageHDU(cube.good_wavelength.astype(np.uint8), name="GOODWAVE"))
    hdus.append(fits.ImageHDU(cube.wavelength.astype(np.float64), name="WAVELENGTH"))

    fits.HDUList(hdus).writeto(destination, overwrite=overwrite)
    return destination



def summarize_binary_mask(mask: np.ndarray) -> dict[str, int]:
    """Return exact counts for a binary valid/invalid stack mask."""
    arr = np.asarray(mask)
    return {
        "valid_zero": int(np.sum(arr == 0)),
        "invalid_nonzero": int(np.sum(arr != 0)),
        "total": int(arr.size),
    }


def summarize_effective_exposure(exposure: np.ndarray | None) -> dict[str, float | int] | None:
    """Return compact coverage/exposure statistics for logs and manifests."""
    if exposure is None:
        return None
    arr = np.asarray(exposure, dtype=float)
    positive = arr[np.isfinite(arr) & (arr > 0)]
    if positive.size == 0:
        return {"positive_samples": 0, "total_samples": int(arr.size)}
    return {
        "positive_samples": int(positive.size),
        "total_samples": int(arr.size),
        "positive_fraction": float(positive.size / arr.size),
        "min_positive_s": float(np.min(positive)),
        "median_positive_s": float(np.median(positive)),
        "max_positive_s": float(np.max(positive)),
    }

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
