#!/usr/bin/env python3
"""CRD_DAP Script 03d: build an externally anchored RH3 atmospheric mask.

Purpose
-------
Script 03b showed that the weak RL development data contain recurrent, very
large negative residuals.  Those residuals are useful *diagnostics*, but they
should not be the primary definition of a production mask.  This script builds
an observed-frame mask from atmospheric information that is independent of the
stellar pPXF fit:

1. empirical sky emission measured from one or more pre-ZAP/sky-model cubes;
2. optional telluric transmission measured/modelled for the observation;
3. optional O2 A-band review interval (review only unless explicitly enabled).

The Script-03b recurrent-residual table is used only as a cross-check.  It does
not decide which wavelengths enter the atmospheric mask.

Recommended KSkyWizard inputs
-----------------------------
The strongest empirical reference is the *pre-ZAP* cropped/non-sky-subtracted
cube used by KSkyWizard (or a ZAP sky-model cube if one was saved).  The script
extracts a robust spatial median sky spectrum.  Supply the KSkyWizard source
mask when available; otherwise the spatial median is still robust when the
astronomical target occupies less than half of the field.

Example
-------
python scripts/03d_build_RH3_atmospheric_mask.py \\
    --script3-run runs/8143-1902_S03_20260828_122251 \\
    --config config/8143-1902.py \\
    --pre-zap-cube /path/to/kr241226_00172_cropped.fits \\
    --pre-zap-cube /path/to/kr241226_00173_cropped.fits \\
    --pre-zap-cube /path/to/kr241226_00174_cropped.fits \\
    --pre-zap-cube /path/to/kr241226_00175_cropped.fits \\
    --bins 231,328

If KSkyWizard/ZAP saved the subtracted sky itself, use --sky-model-cube instead.
A 1-D telluric model may be supplied with --telluric-spectrum.

Statistical contract
--------------------
This stage is read-only with respect to Script 3.  It performs no pPXF fits and
never edits the target configuration.  It writes an ECSV mask table that can be
passed to Script 03c for targeted pPXF refits.  Only after that validation should
RH3_ATMOSPHERIC_MASK_FILE point to the accepted table for a full Script-3 rerun.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import sys
from typing import Iterable

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.table import Table
from scipy import ndimage

import crd_utils as crd


SCRIPT_NAME = "03d_build_RH3_atmospheric_mask"


def _parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Build an RH3 fixed observed-frame mask from empirical sky/telluric references."
    )
    p.add_argument("--script3-run", required=True, help="Completed Script-3 run directory.")
    p.add_argument("--config", required=True, help="Current target config, e.g. config/8143-1902.py.")
    p.add_argument(
        "--pre-zap-cube", action="append", default=[],
        help=(
            "Pre-ZAP/non-sky-subtracted KSkyWizard cube. Repeat for multiple exposures. "
            "A robust spatial-median sky spectrum is extracted from each cube."
        ),
    )
    p.add_argument(
        "--sky-model-cube", action="append", default=[],
        help=(
            "Optional ZAP/KSkyWizard sky-model cube (the sky that was subtracted). Repeatable. "
            "This is treated as a direct empirical sky reference."
        ),
    )
    p.add_argument(
        "--post-zap-cube", action="append", default=[],
        help=(
            "Optional sky-subtracted cube(s), used only to display/measure remaining blank-sky residuals. "
            "They do not define the atmospheric mask."
        ),
    )
    p.add_argument(
        "--object-mask", action="append", default=[],
        help=(
            "Optional 2-D FITS object mask(s) for --pre-zap-cube; object pixels must be >0, sky=0. "
            "Give one mask to reuse for all references or one mask per pre-ZAP cube."
        ),
    )
    p.add_argument(
        "--telluric-spectrum", default=None,
        help=(
            "Optional 1-D telluric transmission reference (FITS image/table or ASCII/ECSV). "
            "If supplied, transmission below --telluric-threshold is independently masked."
        ),
    )
    p.add_argument(
        "--telluric-medium", choices=("air", "vacuum"), default="air",
        help="Wavelength medium of --telluric-spectrum when it cannot be inferred (default air).",
    )
    p.add_argument(
        "--telluric-threshold", type=float, default=0.97,
        help="Mask telluric reference samples with transmission below this value (default 0.97).",
    )
    p.add_argument(
        "--include-known-o2-a-in-mask", action="store_true",
        help=(
            "Also include a conservative O2 A-band review interval (7590--7705 A observed air). "
            "Default is review-only, not included, unless an actual telluric spectrum supports it."
        ),
    )
    p.add_argument(
        "--sky-line-sigma", type=float, default=8.0,
        help="Minimum high-pass empirical-sky significance for a skyline sample (default 8).",
    )
    p.add_argument(
        "--sky-reference-min-fraction", type=float, default=0.50,
        help="Minimum fraction of empirical references that must detect a skyline sample (default 0.50).",
    )
    p.add_argument(
        "--continuum-window-angstrom", type=float, default=75.0,
        help="Median-filter width used to remove broad continuum from sky references (default 75 A).",
    )
    p.add_argument(
        "--mask-padding-log-pixels", type=int, default=2,
        help="Dilate independently detected atmospheric features by this many Script-3 log pixels (default 2).",
    )
    p.add_argument(
        "--mask-bridge-log-pixels", type=int, default=1,
        help="Bridge false gaps of at most this many log pixels inside a detected feature (default 1).",
    )
    p.add_argument(
        "--recurrent-table", dest="recurrent_table", default=None,
        help=(
            "Optional Script-03b recurrent-wavelength ECSV. Default: "
            "<script3-run>/validation/03b/products/RH3_03b_recurrent_wavelengths.ecsv. "
            "Used for QC columns only; never for mask selection."
        ),
    )
    p.add_argument("--bins", default="231,328", help="Comma-separated Script-3 bins for overlay plots.")
    p.add_argument(
        "--output-dir", default=None,
        help="Default: <script3-run>/validation/03d_atmosphere",
    )
    return p


def _load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def _write_json(payload: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, sort_keys=True)
        f.write("\n")
    tmp.replace(path)


def _vacuum_to_air(vacuum_angstrom: np.ndarray) -> np.ndarray:
    vac = np.asarray(vacuum_angstrom, dtype=float)
    sigma2 = (1.0e4 / vac) ** 2
    n = 1.0 + 0.0000834254 + 0.02406147 / (130.0 - sigma2) + 0.00015998 / (38.9 - sigma2)
    return vac / n


def _air_to_vacuum(air_angstrom: np.ndarray) -> np.ndarray:
    air = np.asarray(air_angstrom, dtype=float)
    vac = air * 1.00028
    for _ in range(6):
        predicted_air = _vacuum_to_air(vac)
        vac *= air / predicted_air
    return vac


def _convert_medium(wave: np.ndarray, src: str, dst: str) -> np.ndarray:
    src = str(src).strip().lower()
    dst = str(dst).strip().lower()
    w = np.asarray(wave, dtype=float)
    if src == dst:
        return w.copy()
    if src == "air" and dst == "vacuum":
        return _air_to_vacuum(w)
    if src == "vacuum" and dst == "air":
        return _vacuum_to_air(w)
    raise ValueError(f"Unsupported wavelength-medium conversion {src!r} -> {dst!r}.")


def _rest_template_to_observed_science(
    rest_wave: np.ndarray, redshift: float, template_medium: str, science_medium: str
) -> np.ndarray:
    rest_science = _convert_medium(rest_wave, template_medium, science_medium)
    return rest_science * (1.0 + float(redshift))


def _log_center_edges(center: np.ndarray) -> np.ndarray:
    c = np.asarray(center, dtype=float)
    if c.ndim != 1 or c.size < 2 or np.any(c <= 0):
        raise ValueError("Expected a positive one-dimensional wavelength grid with >=2 samples.")
    logc = np.log(c)
    e = np.empty(c.size + 1, dtype=float)
    e[1:-1] = 0.5 * (logc[:-1] + logc[1:])
    e[0] = logc[0] - 0.5 * (logc[1] - logc[0])
    e[-1] = logc[-1] + 0.5 * (logc[-1] - logc[-2])
    return np.exp(e)


def _contiguous_true_groups(mask: np.ndarray) -> list[tuple[int, int]]:
    m = np.asarray(mask, dtype=bool)
    idx = np.flatnonzero(m)
    if idx.size == 0:
        return []
    breaks = np.where(np.diff(idx) > 1)[0]
    starts = np.r_[0, breaks + 1]
    ends = np.r_[breaks, idx.size - 1]
    return [(int(idx[s]), int(idx[e])) for s, e in zip(starts, ends)]


def _bridge_small_false_gaps(mask: np.ndarray, max_gap: int) -> np.ndarray:
    out = np.asarray(mask, dtype=bool).copy()
    if max_gap <= 0 or out.size < 3:
        return out
    inv = ~out
    for i0, i1 in _contiguous_true_groups(inv):
        if i0 == 0 or i1 == out.size - 1:
            continue
        if (i1 - i0 + 1) <= max_gap and out[i0 - 1] and out[i1 + 1]:
            out[i0:i1 + 1] = True
    return out


def _parse_bins(text: str, nbin: int) -> list[int]:
    out: list[int] = []
    for token in str(text).split(","):
        token = token.strip()
        if not token:
            continue
        bid = int(token)
        if bid < 0 or bid >= nbin:
            raise ValueError(f"Bin {bid} is outside valid range 0..{nbin - 1}.")
        if bid not in out:
            out.append(bid)
    return out


def _first_image_hdu(path: Path, min_ndim: int = 1):
    hdul = fits.open(path, memmap=True)
    for hdu in hdul:
        data = getattr(hdu, "data", None)
        if data is not None and np.ndim(data) >= min_ndim and np.issubdtype(np.asarray(data).dtype, np.number):
            return hdul, hdu
    hdul.close()
    raise ValueError(f"No numeric image HDU with ndim>={min_ndim} found in {path}.")


def _spectral_numpy_axis(header, ndim: int) -> tuple[int, int]:
    naxis = int(header.get("NAXIS", ndim))
    for fits_axis in range(1, naxis + 1):
        ctype = str(header.get(f"CTYPE{fits_axis}", "")).upper()
        if any(token in ctype for token in ("WAVE", "AWAV", "LAMB", "SPEC")):
            np_axis = ndim - fits_axis
            if 0 <= np_axis < ndim:
                return np_axis, fits_axis
    # KCWI 3-D cubes conventionally store wavelength on FITS axis 3 -> numpy axis 0.
    if ndim == 3:
        return 0, 3
    if ndim == 1:
        return 0, 1
    raise ValueError("Could not identify spectral axis from FITS CTYPE keywords.")


def _unit_to_angstrom_factor(unit: str) -> float:
    u = str(unit or "").strip().lower().replace("å", "angstrom")
    if not u or u in {"a", "aa", "angstrom", "angstroms"}:
        return 1.0
    if u in {"nm", "nanometer", "nanometers"}:
        return 10.0
    if u in {"um", "micron", "microns", "micrometer", "micrometers"}:
        return 1.0e4
    if u in {"m", "meter", "meters"}:
        return 1.0e10
    raise ValueError(f"Unsupported spectral CUNIT {unit!r}; expected Angstrom, nm, um, or m.")


def _linear_wavelength_from_header(header, n: int, fits_axis: int) -> np.ndarray:
    crval = header.get(f"CRVAL{fits_axis}")
    crpix = header.get(f"CRPIX{fits_axis}", 1.0)
    cdelt = header.get(f"CDELT{fits_axis}")
    if cdelt is None:
        cdelt = header.get(f"CD{fits_axis}_{fits_axis}")
    if crval is None or cdelt is None:
        raise ValueError(
            f"Cannot build wavelength axis: CRVAL{fits_axis}/CDELT{fits_axis} (or CD diagonal) missing."
        )
    pix = np.arange(n, dtype=float) + 1.0
    wave = float(crval) + (pix - float(crpix)) * float(cdelt)
    factor = _unit_to_angstrom_factor(header.get(f"CUNIT{fits_axis}", "Angstrom"))
    return wave * factor


def _load_object_mask(path: Path, spatial_shape: tuple[int, ...]) -> np.ndarray:
    hdul, hdu = _first_image_hdu(path, min_ndim=2)
    try:
        arr = np.array(hdu.data, copy=True)
    finally:
        hdul.close()
    if arr.ndim > 2:
        arr = np.squeeze(arr)
    if arr.shape != spatial_shape:
        raise ValueError(f"Object mask {path} shape {arr.shape} does not match cube spatial shape {spatial_shape}.")
    return np.asarray(arr > 0, dtype=bool)


def _extract_spatial_median_spectrum(path: Path, object_mask: Path | None = None) -> tuple[np.ndarray, np.ndarray, int]:
    hdul, hdu = _first_image_hdu(path, min_ndim=3)
    try:
        data = np.array(hdu.data, dtype=float, copy=True)
        header = hdu.header.copy()
    finally:
        hdul.close()
    spec_axis, fits_axis = _spectral_numpy_axis(header, data.ndim)
    cube = np.moveaxis(data, spec_axis, 0)
    wave = _linear_wavelength_from_header(header, cube.shape[0], fits_axis)
    spatial_shape = cube.shape[1:]
    if object_mask is not None:
        obj = _load_object_mask(object_mask, spatial_shape)
        sky_sel = ~obj
    else:
        sky_sel = np.ones(spatial_shape, dtype=bool)
    flat = cube.reshape(cube.shape[0], -1)
    sky_flat = sky_sel.reshape(-1)
    n_spax = int(np.count_nonzero(sky_flat))
    if n_spax < 5:
        raise ValueError(f"Only {n_spax} spatial sky samples remain in {path}; need at least 5.")
    # A spatial median is intentionally used: it is robust to faint target light and isolated CRs.
    spec = np.nanmedian(flat[:, sky_flat], axis=1)
    if wave[0] > wave[-1]:
        wave = wave[::-1]
        spec = spec[::-1]
    return np.asarray(wave, dtype=float), np.asarray(spec, dtype=float), n_spax


def _robust_highpass_z(wave: np.ndarray, spectrum: np.ndarray, window_angstrom: float) -> tuple[np.ndarray, np.ndarray, float]:
    w = np.asarray(wave, dtype=float)
    s = np.asarray(spectrum, dtype=float)
    finite = np.isfinite(w) & np.isfinite(s)
    if np.count_nonzero(finite) < 20:
        raise ValueError("Atmospheric reference contains too few finite samples.")
    dw = float(np.nanmedian(np.diff(w[finite])))
    if not np.isfinite(dw) or dw <= 0:
        raise ValueError("Atmospheric reference wavelength grid is not strictly increasing.")
    size = max(5, int(round(float(window_angstrom) / dw)))
    if size % 2 == 0:
        size += 1
    filled = s.copy()
    med = float(np.nanmedian(filled[finite]))
    filled[~finite] = med
    continuum = ndimage.median_filter(filled, size=size, mode="nearest")
    hp = filled - continuum
    # Estimate the local reference noise from the non-positive half so bright sky lines do not inflate it.
    center = float(np.nanmedian(hp[finite]))
    lower = hp[finite & (hp <= center)]
    if lower.size >= 10:
        sigma = 1.4826 * float(np.nanmedian(np.abs(lower - np.nanmedian(lower))))
    else:
        sigma = 1.4826 * float(np.nanmedian(np.abs(hp[finite] - center)))
    if not np.isfinite(sigma) or sigma <= 0:
        sigma = float(np.nanstd(hp[finite]))
    if not np.isfinite(sigma) or sigma <= 0:
        raise ValueError("Could not estimate finite positive scatter in atmospheric reference.")
    z = hp / sigma
    z[~finite] = np.nan
    hp[~finite] = np.nan
    return hp, z, sigma


def _load_telluric_spectrum(path: Path, fallback_medium: str) -> tuple[np.ndarray, np.ndarray, str]:
    suffix = path.suffix.lower()
    if suffix in {".fits", ".fit", ".fz"}:
        with fits.open(path, memmap=True) as hdul:
            # Prefer a table with obvious wavelength/transmission columns.
            for hdu in hdul:
                data = getattr(hdu, "data", None)
                names = list(getattr(data, "names", []) or []) if data is not None else []
                if names:
                    low = {str(n).lower(): n for n in names}
                    wname = next((low[k] for k in low if any(t in k for t in ("wave", "lambda", "lam"))), None)
                    tname = next((low[k] for k in low if any(t in k for t in ("trans", "telluric", "model"))), None)
                    if wname is not None and tname is not None:
                        wave = np.asarray(data[wname], dtype=float).ravel()
                        trans = np.asarray(data[tname], dtype=float).ravel()
                        medium = str(hdu.header.get("WAVEMED", hdu.header.get("AIRORVAC", fallback_medium))).strip().lower()
                        if medium not in {"air", "vacuum"}:
                            medium = fallback_medium
                        order = np.argsort(wave)
                        return wave[order], trans[order], medium
            # Otherwise accept a 1-D image with spectral WCS.
            for hdu in hdul:
                data = getattr(hdu, "data", None)
                if data is not None and np.ndim(data) == 1:
                    trans = np.asarray(data, dtype=float)
                    _, fits_axis = _spectral_numpy_axis(hdu.header, 1)
                    wave = _linear_wavelength_from_header(hdu.header, trans.size, fits_axis)
                    medium = str(hdu.header.get("WAVEMED", hdu.header.get("AIRORVAC", fallback_medium))).strip().lower()
                    if medium not in {"air", "vacuum"}:
                        medium = fallback_medium
                    order = np.argsort(wave)
                    return wave[order], trans[order], medium
        raise ValueError(f"Could not identify wavelength/transmission arrays in {path}.")

    tab = Table.read(path)
    names = {str(n).lower(): n for n in tab.colnames}
    wname = next((names[k] for k in names if any(t in k for t in ("wave", "lambda", "lam"))), None)
    tname = next((names[k] for k in names if any(t in k for t in ("trans", "telluric", "model"))), None)
    if wname is None or tname is None:
        if len(tab.colnames) < 2:
            raise ValueError(f"Telluric table {path} needs wavelength and transmission columns.")
        wname, tname = tab.colnames[:2]
    wave = np.asarray(tab[wname], dtype=float)
    trans = np.asarray(tab[tname], dtype=float)
    medium = str(tab.meta.get("WAVELENGTH_MEDIUM", fallback_medium)).strip().lower()
    if medium not in {"air", "vacuum"}:
        medium = fallback_medium
    order = np.argsort(wave)
    return wave[order], trans[order], medium


def _load_recurrent_qc(path: Path | None, nlog: int) -> dict[str, np.ndarray] | None:
    if path is None or not path.is_file():
        return None
    tab = Table.read(path, format="ascii.ecsv")
    if len(tab) != nlog:
        raise ValueError(f"03b recurrent table has {len(tab)} rows but Script-3 log grid has {nlog} pixels.")
    out: dict[str, np.ndarray] = {}
    for name in (
        "FRACTION_BINS_RAW_NEGATIVE_OUTLIER",
        "EXCESS_N_BINS_TWO_NONPATHOLOGY",
        "FRACTION_TWO_NONPATHOLOGY",
        "N_BINS_RAW_NEGATIVE_OUTLIER",
    ):
        if name in tab.colnames:
            out[name] = np.asarray(tab[name], dtype=float)
    return out


def _load_script3(run: Path) -> tuple[dict, dict, dict]:
    manifest_path = run / "metadata" / "script03_manifest.json"
    spectra_path = run / "products" / "RH3_log_spectra_and_local_best_fits.npz"
    if not manifest_path.is_file() or not spectra_path.is_file():
        raise FileNotFoundError("Script 03d requires a completed Script-3 run with manifest and saved spectra product.")
    manifest = _load_json(manifest_path)
    with np.load(spectra_path, allow_pickle=False) as z:
        spectra = {k: np.asarray(z[k]) for k in z.files}
    return manifest, spectra, {"manifest": str(manifest_path), "spectra": str(spectra_path)}


def _resolve_masks_for_references(paths: list[str], masks: list[str]) -> list[Path | None]:
    if not paths:
        return []
    if not masks:
        return [None] * len(paths)
    if len(masks) == 1:
        p = Path(masks[0]).expanduser().resolve()
        return [p] * len(paths)
    if len(masks) != len(paths):
        raise ValueError("Give zero object masks, one reusable mask, or exactly one mask per --pre-zap-cube.")
    return [Path(x).expanduser().resolve() for x in masks]


def _interpolate_reference_to_grid(
    wave: np.ndarray, spectrum: np.ndarray, grid: np.ndarray
) -> np.ndarray:
    finite = np.isfinite(wave) & np.isfinite(spectrum)
    if np.count_nonzero(finite) < 5:
        return np.full_like(grid, np.nan, dtype=float)
    w = np.asarray(wave[finite], dtype=float)
    s = np.asarray(spectrum[finite], dtype=float)
    order = np.argsort(w)
    w, s = w[order], s[order]
    out = np.interp(grid, w, s, left=np.nan, right=np.nan)
    out[(grid < w[0]) | (grid > w[-1])] = np.nan
    return out


def _build_interval_table(
    mask: np.ndarray,
    sky_mask: np.ndarray,
    tell_mask: np.ndarray,
    o2_mask: np.ndarray,
    wave_rest: np.ndarray,
    observed_science: np.ndarray,
    redshift: float,
    template_medium: str,
    science_medium: str,
    sky_z: np.ndarray,
    sky_fraction: np.ndarray,
    telluric_trans: np.ndarray | None,
    recurrence: dict[str, np.ndarray] | None,
) -> Table:
    rest_edges = _log_center_edges(wave_rest)
    obs_science_edges = _log_center_edges(observed_science)
    if science_medium == "vacuum":
        obs_vac_edges = obs_science_edges
        obs_air_edges = _vacuum_to_air(obs_vac_edges)
    else:
        obs_air_edges = obs_science_edges
        obs_vac_edges = _air_to_vacuum(obs_air_edges)

    rows = []
    for iid, (i0, i1) in enumerate(_contiguous_true_groups(mask)):
        sl = slice(i0, i1 + 1)
        has_sky = bool(np.any(sky_mask[sl]))
        has_tel = bool(np.any(tell_mask[sl]))
        has_o2 = bool(np.any(o2_mask[sl]))
        reasons = []
        if has_sky:
            reasons.append("EMPIRICAL_SKY_EMISSION")
        if has_tel:
            reasons.append("TELLURIC_REFERENCE")
        if has_o2:
            reasons.append("KNOWN_O2_A")
        reason = "+".join(reasons) if reasons else "ATMOSPHERIC_REFERENCE"
        peak_z = float(np.nanmax(sky_z[sl])) if np.any(np.isfinite(sky_z[sl])) else np.nan
        peak_frac = float(np.nanmax(sky_fraction[sl])) if np.any(np.isfinite(sky_fraction[sl])) else np.nan
        min_trans = (
            float(np.nanmin(telluric_trans[sl]))
            if telluric_trans is not None and np.any(np.isfinite(telluric_trans[sl]))
            else np.nan
        )
        raw_frac = np.nan
        raw_count = np.nan
        resid_excess = np.nan
        if recurrence is not None:
            if "FRACTION_BINS_RAW_NEGATIVE_OUTLIER" in recurrence:
                raw_frac = float(np.nanmax(recurrence["FRACTION_BINS_RAW_NEGATIVE_OUTLIER"][sl]))
            if "N_BINS_RAW_NEGATIVE_OUTLIER" in recurrence:
                raw_count = float(np.nanmax(recurrence["N_BINS_RAW_NEGATIVE_OUTLIER"][sl]))
            if "EXCESS_N_BINS_TWO_NONPATHOLOGY" in recurrence:
                resid_excess = float(np.nanmax(recurrence["EXCESS_N_BINS_TWO_NONPATHOLOGY"][sl]))
        rows.append((
            iid, True, reason, i0, i1, i1 - i0 + 1,
            float(rest_edges[i0]), float(rest_edges[i1 + 1]),
            float(obs_air_edges[i0]), float(obs_air_edges[i1 + 1]),
            float(obs_vac_edges[i0]), float(obs_vac_edges[i1 + 1]),
            float(obs_science_edges[i0]), float(obs_science_edges[i1 + 1]),
            peak_z, peak_frac, min_trans, raw_count, raw_frac, resid_excess,
        ))
    tab = Table(rows=rows, names=(
        "INTERVAL_ID", "INCLUDED_IN_ATMOSPHERIC_MASK", "REASON",
        "LOG_INDEX_LO", "LOG_INDEX_HI", "N_LOG_PIXELS",
        "REST_LO_ANGSTROM", "REST_HI_ANGSTROM",
        "OBSERVED_AIR_LO_ANGSTROM", "OBSERVED_AIR_HI_ANGSTROM",
        "OBSERVED_VACUUM_LO_ANGSTROM", "OBSERVED_VACUUM_HI_ANGSTROM",
        "OBSERVED_SCIENCE_LO_ANGSTROM", "OBSERVED_SCIENCE_HI_ANGSTROM",
        "PEAK_EMPIRICAL_SKY_Z", "PEAK_SKY_REFERENCE_DETECTION_FRACTION",
        "MIN_TELLURIC_TRANSMISSION", "QC_PEAK_RAW_NEGATIVE_BIN_COUNT",
        "QC_PEAK_RAW_NEGATIVE_BIN_FRACTION", "QC_PEAK_2C_NONPATHOLOGY_EXCESS",
    ))
    tab.meta["MASK_SELECTION_BASIS"] = "EMPIRICAL_ATMOSPHERIC_REFERENCE_ONLY"
    tab.meta["REDSHIFT"] = float(redshift)
    tab.meta["REST_WAVELENGTH_MEDIUM"] = str(template_medium)
    tab.meta["OBSERVED_SCIENCE_WAVELENGTH_MEDIUM"] = str(science_medium)
    tab.meta["RECURRENCE_USED_FOR_SELECTION"] = False
    return tab


def _write_config_snippet(path: Path, table: Table, science_medium: str) -> None:
    lines = [
        "# Accepted atmospheric mask generated by Script 03d.",
        "# Intervals are selected from empirical sky/telluric references, not from pPXF residuals.",
        f"# Wavelength medium: observed {science_medium} (native Script-3 science medium).",
        "RH3_MASK_OBSERVED_RANGES_ANGSTROM = [",
    ]
    for row in table:
        if not bool(row["INCLUDED_IN_ATMOSPHERIC_MASK"]):
            continue
        lo = float(row["OBSERVED_SCIENCE_LO_ANGSTROM"])
        hi = float(row["OBSERVED_SCIENCE_HI_ANGSTROM"])
        lines.append(f"    ({lo:.2f}, {hi:.2f}),  # {row['REASON']}")
    lines.append("]")
    lines.append("")
    lines.append("# Preferred production integration (avoids copying numbers manually):")
    lines.append(f"# RH3_ATMOSPHERIC_MASK_FILE = Path({str(path.with_name('RH3_03d_atmospheric_mask.ecsv'))!r})")
    lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")


def _plot_atmospheric_overlay(
    path: Path,
    obs: np.ndarray,
    sky_z: np.ndarray,
    sky_fraction: np.ndarray,
    sky_mask: np.ndarray,
    telluric_trans: np.ndarray | None,
    tell_mask: np.ndarray,
    final_mask: np.ndarray,
    recurrence: dict[str, np.ndarray] | None,
    post_resid_z: np.ndarray | None,
    science_medium: str,
) -> None:
    fig, axes = plt.subplots(4, 1, figsize=(13, 11), sharex=True)
    ax = axes[0]
    ax.plot(obs, sky_z, lw=0.9, label="median empirical pre-ZAP/sky-model high-pass significance")
    ax.axhline(0.0, lw=0.6)
    ax.set_ylabel("sky high-pass / robust sigma")
    ax.legend(fontsize=8)
    ax2 = ax.twinx()
    ax2.plot(obs, sky_fraction, lw=0.7, alpha=0.5, label="reference detection fraction")
    ax2.set_ylim(-0.02, 1.02)
    ax2.set_ylabel("reference fraction")

    ax = axes[1]
    if post_resid_z is not None:
        ax.plot(obs, post_resid_z, lw=0.8, label="post-ZAP blank-sky residual high-pass significance")
        ax.axhline(0.0, lw=0.6)
        ax.legend(fontsize=8)
    else:
        ax.text(0.5, 0.5, "No --post-zap-cube supplied", ha="center", va="center", transform=ax.transAxes)
    ax.set_ylabel("post-ZAP residual / sigma")

    ax = axes[2]
    if recurrence is not None:
        if "FRACTION_BINS_RAW_NEGATIVE_OUTLIER" in recurrence:
            ax.plot(obs, recurrence["FRACTION_BINS_RAW_NEGATIVE_OUTLIER"], lw=0.9, label="03b raw-negative bin fraction")
        if "EXCESS_N_BINS_TWO_NONPATHOLOGY" in recurrence:
            exc = recurrence["EXCESS_N_BINS_TWO_NONPATHOLOGY"]
            scale = np.nanmax(np.abs(exc))
            if np.isfinite(scale) and scale > 0:
                ax.plot(obs, exc / scale, lw=0.8, label="03b 2C localized excess (scaled)")
        ax.legend(fontsize=8)
    else:
        ax.text(0.5, 0.5, "No Script-03b recurrence table available", ha="center", va="center", transform=ax.transAxes)
    ax.set_ylabel("independent QC")

    ax = axes[3]
    if telluric_trans is not None:
        ax.plot(obs, telluric_trans, lw=0.9, label="telluric transmission")
        ax.set_ylim(-0.05, 1.05)
        ax.legend(fontsize=8)
    else:
        ax.text(0.5, 0.5, "No telluric spectrum supplied; O2 A band is review-only by default", ha="center", va="center", transform=ax.transAxes)
    ax.set_ylabel("transmission")
    ax.set_xlabel(f"Observed-frame wavelength ({science_medium} Angstrom)")

    for ax in axes:
        for i0, i1 in _contiguous_true_groups(final_mask):
            ax.axvspan(obs[i0], obs[i1], alpha=0.12)
    axes[0].set_title("RH3 atmospheric-reference mask: selection independent of pPXF residuals")
    fig.tight_layout()
    fig.savefig(path, dpi=170)
    plt.close(fig)


def _robust_ylim(galaxy: np.ndarray, good: np.ndarray) -> tuple[float, float]:
    vals = np.asarray(galaxy, dtype=float)[np.asarray(good, dtype=bool) & np.isfinite(galaxy)]
    if vals.size < 5:
        return -2.0, 2.0
    lo, hi = np.percentile(vals, [1.0, 99.0])
    pad = 0.15 * max(hi - lo, 1.0)
    return float(lo - pad), float(hi + pad)


def _plot_bin_overlay(path: Path, bid: int, wave_rest: np.ndarray, galaxy: np.ndarray, model: np.ndarray, good: np.ndarray, mask: np.ndarray) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(13, 7), sharex=True)
    retained = np.asarray(good, dtype=bool) & ~mask
    removed = np.asarray(good, dtype=bool) & mask
    axes[0].plot(wave_rest, galaxy, lw=0.8, label="saved RH3 spectrum")
    axes[0].plot(wave_rest, model, lw=0.9, label="saved local-best 2C model")
    if np.any(removed):
        axes[0].scatter(wave_rest[removed], galaxy[removed], marker="x", s=22, label="atmospheric-reference mask")
    ylim = _robust_ylim(galaxy, retained)
    axes[0].set_ylim(*ylim)
    axes[0].set_ylabel("Normalized flux")
    axes[0].set_title(f"Bin {bid}: externally anchored atmospheric mask (saved model; no refit)")
    axes[0].legend(fontsize=8)
    resid = np.full_like(galaxy, np.nan, dtype=float)
    valid = np.asarray(good, dtype=bool) & np.isfinite(galaxy) & np.isfinite(model)
    resid[valid] = galaxy[valid] - model[valid]
    axes[1].plot(wave_rest[retained], resid[retained], lw=0.8, label="retained residual")
    if np.any(removed):
        axes[1].scatter(wave_rest[removed], resid[removed], marker="x", s=22, label="removed residual")
    axes[1].axhline(0.0, lw=0.6)
    axes[1].set_ylabel("data-model")
    axes[1].set_xlabel("Rest-frame wavelength (Angstrom)")
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=170)
    plt.close(fig)


def main() -> int:
    args = _parser().parse_args()
    run = Path(args.script3_run).expanduser().resolve()
    cfg = crd.config.load_config(args.config)
    out = Path(args.output_dir).expanduser().resolve() if args.output_dir else run / "validation" / "03d_atmosphere"
    products = out / "products"
    figures = out / "figures"
    metadata = out / "metadata"
    for d in (products, figures, metadata):
        d.mkdir(parents=True, exist_ok=True)

    manifest, spectra, source_paths = _load_script3(run)
    wave_rest = np.asarray(spectra["wavelength"], dtype=float)
    galaxy_all = np.asarray(spectra["galaxy"], dtype=float)
    good_all = np.asarray(spectra["good"], dtype=bool)
    two_all = np.asarray(spectra["local_best_two_component_model"], dtype=float)
    nbin, nlog = galaxy_all.shape
    if wave_rest.size != nlog:
        raise ValueError("Saved Script-3 wavelength and galaxy arrays disagree in length.")

    redshift = float(manifest.get("redshift", cfg.REDSHIFT))
    template_medium = str(manifest.get("template_wavelength_medium", getattr(cfg, "TEMPLATE_WAVELENGTH_MEDIUM", "air"))).lower()
    science_medium = str(manifest.get("science_wavelength_medium", "")).lower()
    if science_medium not in {"air", "vacuum"}:
        raise ValueError("Script-3 manifest must record science_wavelength_medium as air or vacuum.")
    if template_medium not in {"air", "vacuum"}:
        raise ValueError("Script-3 manifest/config must record template wavelength medium as air or vacuum.")
    observed_science = _rest_template_to_observed_science(wave_rest, redshift, template_medium, science_medium)

    pre_paths = [Path(x).expanduser().resolve() for x in args.pre_zap_cube]
    sky_model_paths = [Path(x).expanduser().resolve() for x in args.sky_model_cube]
    post_paths = [Path(x).expanduser().resolve() for x in args.post_zap_cube]
    if not pre_paths and not sky_model_paths and args.telluric_spectrum is None and not args.include_known_o2_a_in_mask:
        raise ValueError(
            "No independent atmospheric reference was supplied. Give --pre-zap-cube, --sky-model-cube, "
            "--telluric-spectrum, or explicitly --include-known-o2-a-in-mask."
        )
    masks = _resolve_masks_for_references(args.pre_zap_cube, args.object_mask)

    reference_specs = []
    reference_meta = []
    for p, m in zip(pre_paths, masks):
        w, s, nsky = _extract_spatial_median_spectrum(p, m)
        reference_specs.append((w, s))
        reference_meta.append({"path": str(p), "kind": "pre-zap", "object_mask": None if m is None else str(m), "n_sky_spaxels": nsky})
    for p in sky_model_paths:
        w, s, nsky = _extract_spatial_median_spectrum(p, None)
        reference_specs.append((w, s))
        reference_meta.append({"path": str(p), "kind": "sky-model", "object_mask": None, "n_sky_spaxels": nsky})

    sky_z_stack = []
    sky_flux_stack = []
    sky_detect_stack = []
    for w, s in reference_specs:
        interp = _interpolate_reference_to_grid(w, s, observed_science)
        hp, z, _ = _robust_highpass_z(observed_science, interp, float(args.continuum_window_angstrom))
        sky_flux_stack.append(interp)
        sky_z_stack.append(z)
        sky_detect_stack.append(np.asarray(z >= float(args.sky_line_sigma), dtype=bool))

    if sky_z_stack:
        zarr = np.asarray(sky_z_stack, dtype=float)
        sky_z = np.nanmedian(zarr, axis=0)
        detect_fraction = np.mean(np.asarray(sky_detect_stack, dtype=float), axis=0)
        sky_reference_flux = np.nanmedian(np.asarray(sky_flux_stack, dtype=float), axis=0)
        sky_base = (
            np.isfinite(sky_z)
            & (sky_z >= float(args.sky_line_sigma))
            & (detect_fraction >= float(args.sky_reference_min_fraction))
        )
    else:
        sky_z = np.full(nlog, np.nan)
        detect_fraction = np.zeros(nlog, dtype=float)
        sky_reference_flux = np.full(nlog, np.nan)
        sky_base = np.zeros(nlog, dtype=bool)

    sky_mask = _bridge_small_false_gaps(sky_base, int(args.mask_bridge_log_pixels))
    if int(args.mask_padding_log_pixels) > 0 and np.any(sky_mask):
        sky_mask = ndimage.binary_dilation(sky_mask, iterations=int(args.mask_padding_log_pixels))

    # Post-ZAP spectra are QC only; they never select mask pixels.
    #
    # A well sky-subtracted *stacked* post-ZAP cube can have an almost perfectly
    # flat spatial-median spectrum after interpolation to the Script-3 grid.  In
    # that case the robust high-pass scatter is exactly (or numerically) zero.
    # That is not a science failure and must never abort construction of the
    # independent atmospheric mask, which has already been determined from the
    # pre-ZAP / sky-model references above.
    post_z_list = []
    post_flux_list = []
    post_qc_meta = []
    for p in post_paths:
        w, s, _ = _extract_spatial_median_spectrum(p, None)
        interp = _interpolate_reference_to_grid(w, s, observed_science)
        try:
            _, z, sigma_post = _robust_highpass_z(
                observed_science,
                interp,
                float(args.continuum_window_angstrom),
            )
        except ValueError as exc:
            post_qc_meta.append(
                {
                    "path": str(p),
                    "status": "skipped_degenerate_highpass",
                    "reason": str(exc),
                }
            )
            # Preserve the interpolated post-ZAP spectrum for the output table,
            # but omit a standardized-z curve because no finite positive scatter
            # exists with which to standardize it.
            post_flux_list.append(interp)
            print(
                f"WARNING: post-ZAP QC reference {p.name} has no finite positive "
                f"high-pass scatter and will be skipped for POST_ZAP_RESIDUAL_Z "
                f"diagnostics only. Atmospheric-mask selection is unaffected."
            )
            continue

        post_z_list.append(z)
        post_flux_list.append(interp)
        post_qc_meta.append(
            {
                "path": str(p),
                "status": "used",
                "highpass_sigma": float(sigma_post),
            }
        )

    post_resid_z = (
        np.nanmedian(np.asarray(post_z_list, dtype=float), axis=0)
        if post_z_list
        else None
    )

    telluric_trans = None
    tell_mask = np.zeros(nlog, dtype=bool)
    telluric_meta = None
    if args.telluric_spectrum is not None:
        tpath = Path(args.telluric_spectrum).expanduser().resolve()
        tw, tt, tmed = _load_telluric_spectrum(tpath, args.telluric_medium)
        tw_science = _convert_medium(tw, tmed, science_medium)
        telluric_trans = _interpolate_reference_to_grid(tw_science, tt, observed_science)
        tell_mask = np.isfinite(telluric_trans) & (telluric_trans < float(args.telluric_threshold))
        tell_mask = _bridge_small_false_gaps(tell_mask, int(args.mask_bridge_log_pixels))
        if int(args.mask_padding_log_pixels) > 0 and np.any(tell_mask):
            tell_mask = ndimage.binary_dilation(tell_mask, iterations=int(args.mask_padding_log_pixels))
        telluric_meta = {"path": str(tpath), "input_medium": tmed, "threshold": float(args.telluric_threshold)}

    # A-band is documented/reviewed independently, but is not auto-included by default.
    obs_air = _vacuum_to_air(observed_science) if science_medium == "vacuum" else observed_science
    o2_review = (obs_air >= 7590.0) & (obs_air <= 7705.0)
    o2_mask = o2_review.copy() if args.include_known_o2_a_in_mask else np.zeros(nlog, dtype=bool)

    final_mask = sky_mask | tell_mask | o2_mask
    if not np.any(final_mask):
        raise RuntimeError("Independent atmospheric references selected zero Script-3 log pixels.")

    recurrent_path = (
        Path(args.recurrent_table).expanduser().resolve()
        if args.recurrent_table
        else run / "validation" / "03b" / "products" / "RH3_03b_recurrent_wavelengths.ecsv"
    )
    recurrence = _load_recurrent_qc(recurrent_path if recurrent_path.is_file() else None, nlog)

    table = _build_interval_table(
        final_mask, sky_mask, tell_mask, o2_mask, wave_rest, observed_science,
        redshift, template_medium, science_medium, sky_z, detect_fraction,
        telluric_trans, recurrence,
    )
    table.meta["SKY_LINE_SIGMA_THRESHOLD"] = float(args.sky_line_sigma)
    table.meta["SKY_REFERENCE_MIN_FRACTION"] = float(args.sky_reference_min_fraction)
    table.meta["CONTINUUM_WINDOW_ANGSTROM"] = float(args.continuum_window_angstrom)
    table.meta["MASK_PADDING_LOG_PIXELS"] = int(args.mask_padding_log_pixels)
    table.meta["MASK_BRIDGE_LOG_PIXELS"] = int(args.mask_bridge_log_pixels)
    table.meta["O2_A_INCLUDED"] = bool(args.include_known_o2_a_in_mask)
    mask_path = products / "RH3_03d_atmospheric_mask.ecsv"
    table.write(mask_path, format="ascii.ecsv", overwrite=True)

    reference_table = Table()
    reference_table["REST_WAVELENGTH_ANGSTROM"] = wave_rest
    reference_table["OBSERVED_SCIENCE_WAVELENGTH_ANGSTROM"] = observed_science
    reference_table["EMPIRICAL_SKY_MEDIAN"] = sky_reference_flux
    reference_table["EMPIRICAL_SKY_Z"] = sky_z
    reference_table["SKY_REFERENCE_DETECTION_FRACTION"] = detect_fraction
    reference_table["ATMOSPHERIC_MASK"] = final_mask
    reference_table["SKY_EMISSION_MASK"] = sky_mask
    reference_table["TELLURIC_MASK"] = tell_mask
    reference_table["O2_A_REVIEW"] = o2_review
    if telluric_trans is not None:
        reference_table["TELLURIC_TRANSMISSION"] = telluric_trans
    if post_resid_z is not None:
        reference_table["POST_ZAP_RESIDUAL_Z"] = post_resid_z
    if recurrence is not None:
        for key, arr in recurrence.items():
            reference_table[f"QC_{key}"] = arr
    reference_table.meta["OBSERVED_SCIENCE_WAVELENGTH_MEDIUM"] = science_medium
    reference_table.meta["REDSHIFT"] = redshift
    reference_table.meta["MASK_SELECTION_BASIS"] = "EMPIRICAL_ATMOSPHERIC_REFERENCE_ONLY"
    reference_table.write(products / "RH3_03d_atmospheric_reference_spectrum.ecsv", format="ascii.ecsv", overwrite=True)

    snippet = products / "RH3_03d_atmospheric_mask_config_snippet.txt"
    _write_config_snippet(snippet, table, science_medium)

    _plot_atmospheric_overlay(
        figures / "RH3_03d_atmospheric_reference_overlay.png",
        observed_science, sky_z, detect_fraction, sky_mask,
        telluric_trans, tell_mask, final_mask, recurrence, post_resid_z, science_medium,
    )
    for bid in _parse_bins(args.bins, nbin):
        _plot_bin_overlay(
            figures / f"RH3_03d_bin_{bid:04d}_atmospheric_mask.png",
            bid, wave_rest, galaxy_all[bid], two_all[bid], good_all[bid], final_mask,
        )

    manifest_out = {
        "script": SCRIPT_NAME,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source_script3_run": str(run),
        "config": str(Path(args.config).expanduser()),
        "redshift": redshift,
        "science_wavelength_medium": science_medium,
        "template_wavelength_medium": template_medium,
        "mask_selection_basis": "empirical_atmospheric_reference_only",
        "recurrence_used_for_selection": False,
        "pre_zap_references": reference_meta,
        "sky_model_references": [str(p) for p in sky_model_paths],
        "post_zap_qc_references": [str(p) for p in post_paths],
        "post_zap_qc_status": post_qc_meta,
        "telluric_reference": telluric_meta,
        "known_o2_a_review_air_angstrom": [7590.0, 7705.0],
        "known_o2_a_included": bool(args.include_known_o2_a_in_mask),
        "n_log_pixels": int(nlog),
        "n_masked_log_pixels": int(np.count_nonzero(final_mask)),
        "masked_fit_fraction": float(np.mean(final_mask)),
        "n_intervals": int(len(table)),
        "sky_line_sigma": float(args.sky_line_sigma),
        "sky_reference_min_fraction": float(args.sky_reference_min_fraction),
        "continuum_window_angstrom": float(args.continuum_window_angstrom),
        "mask_padding_log_pixels": int(args.mask_padding_log_pixels),
        "mask_bridge_log_pixels": int(args.mask_bridge_log_pixels),
        "03b_recurrent_table_for_qc": str(recurrent_path) if recurrent_path.is_file() else None,
        "outputs": {
            "atmospheric_mask_ecsv": str(mask_path),
            "reference_spectrum_ecsv": str(products / "RH3_03d_atmospheric_reference_spectrum.ecsv"),
            "config_snippet_txt": str(snippet),
            "overlay_png": str(figures / "RH3_03d_atmospheric_reference_overlay.png"),
        },
        "source_products": source_paths,
        "warning": (
            "03d defines a candidate from atmospheric references independent of pPXF. "
            "Validate with 03c selected/full-grid tests before using it in a complete Script-3 rerun."
        ),
    }
    _write_json(manifest_out, metadata / "script03d_atmospheric_mask_manifest.json")

    print(f"Atmospheric mask: {len(table)} interval(s), {np.count_nonzero(final_mask)}/{nlog} log pixels "
          f"({100.0*np.mean(final_mask):.2f}% of fit grid).")
    print(f"Written: {mask_path}")
    print("Selection is based on empirical sky/telluric references; Script-03b recurrence is QC only.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
