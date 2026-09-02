#!/usr/bin/env python3
"""CRD_DAP Script 03e: trace surviving RH3 bad pixels back to individual exposures.

Purpose
-------
Scripts 03b--03d separate recurrent fitting pathology from fixed atmospheric
contamination.  Script 03e addresses the remaining question:

    Where do the still-catastrophic, *unmasked* spectral samples come from?

For explicitly selected PowerBins this script traces the strongest surviving
negative outliers in the saved Script-3 spectrum back through the individual
KSkyWizard products.  It compares, exposure by exposure:

1. the pre-ZAP cube (optional, but strongly recommended);
2. the corresponding post-ZAP ``*_zap_icubes.fits`` cube;
3. every unique native spatial pixel intersecting the BL-defined physical
   PowerBin, mapped by celestial WCS rather than by assuming identical arrays.

This allows a surviving bad sample to be distinguished qualitatively as, e.g.:

* strong positive pre-ZAP sky feature -> negative post-ZAP residual:
  sky-over-subtraction suspect;
* present in only one post-ZAP exposure:
  exposure-specific transient/reduction defect suspect;
* present in all exposures but only a small spatial subset:
  stable detector/spatial defect suspect;
* absent from individual post-ZAP cubes but present in Script 3:
  KcwiKit stack/reprojection/coaddition path suspect;
* recurrent across exposures and spatially coherent with no strong pre-ZAP sky:
  stable calibration/wavelength-dependent systematic suspect.

No automatic science mask is created by this script.  The classifications are
diagnostic labels only.

Ca II triplet audit
-------------------
The script also performs an explicit CaT information-preservation audit using
the standard air-vacuum CaT rest wavelengths (8498.02, 8542.09, 8662.14 A).
It reports whether each line is in the *current* Script-3 fit interval and
whether the accepted atmospheric mask overlaps its line center, +/-5 A core,
or +/-20 A kinematic-information window.  This is deliberately a warning/audit:
a genuinely contaminated CaT wavelength must not be silently retained just to
"protect" the line, but an atmospheric mask must also not erase CaT information
without making that loss explicit.

Typical Mac usage
-----------------
python scripts/03e_trace_RH3_bad_pixels.py \
    --script3-run /Users/maxpiper/CRD_DAP/runs/8143-1902_S03_20260828_122251 \
    --config config/8143-1902.py \
    --pre-zap-cube /Volumes/Mainframe/KCWI/KcwiKit_Redux/redux/kr241226_00172_icubed.fits \
    --pre-zap-cube /Volumes/Mainframe/KCWI/KcwiKit_Redux/redux/kr241226_00173_icubed.fits \
    --pre-zap-cube /Volumes/Mainframe/KCWI/KcwiKit_Redux/redux/kr241226_00174_icubed.fits \
    --pre-zap-cube /Volumes/Mainframe/KCWI/KcwiKit_Redux/redux/kr241226_00175_icubed.fits \
    --post-zap-cube /Volumes/Mainframe/KCWI/KcwiKit_Redux/redux/kskywizard/kr241226_00172_zap_icubes.fits \
    --post-zap-cube /Volumes/Mainframe/KCWI/KcwiKit_Redux/redux/kskywizard/kr241226_00173_zap_icubes.fits \
    --post-zap-cube /Volumes/Mainframe/KCWI/KcwiKit_Redux/redux/kskywizard/kr241226_00174_zap_icubes.fits \
    --post-zap-cube /Volumes/Mainframe/KCWI/KcwiKit_Redux/redux/kskywizard/kr241226_00175_zap_icubes.fits \
    --bins 231,328
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import math
from pathlib import Path
from pathlib import PureWindowsPath
import re

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales

import crd_utils as crd


SCRIPT_NAME = "03e_trace_RH3_bad_pixels"
CAT_REST_AIR_ANGSTROM = np.array([8498.02, 8542.09, 8662.14], dtype=float)
CAT_NAMES = ("CaT8498", "CaT8542", "CaT8662")


def _parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Trace surviving RH3 catastrophic pixels through individual pre/post-ZAP exposures."
    )
    p.add_argument("--script3-run", required=True, help="Completed Script-3 run directory.")
    p.add_argument(
        "--script2-run",
        default=None,
        help=(
            "Script-2 run containing products/master_bins.fits. If omitted, 03e first "
            "tries the Script-3 manifest path and then a same-basename run under RUNS_ROOT."
        ),
    )
    p.add_argument("--config", required=True, help="Current target config, e.g. config/8143-1902.py.")
    p.add_argument(
        "--post-zap-cube", action="append", default=[], required=True,
        help="Individual KSkyWizard post-ZAP *_zap_icubes.fits cube. Repeat once per exposure.",
    )
    p.add_argument(
        "--pre-zap-cube", action="append", default=[],
        help=(
            "Matching pre-ZAP *_icubed.fits cube. Repeat once per exposure. "
            "Strongly recommended because pre-vs-post sign changes are the cleanest "
            "test for sky over-subtraction."
        ),
    )
    p.add_argument(
        "--atmospheric-mask",
        default=None,
        help=(
            "Accepted Script-03d atmospheric mask ECSV. Default: "
            "<script3-run>/validation/03d_atmosphere/products/RH3_03d_atmospheric_mask.ecsv."
        ),
    )
    p.add_argument("--bins", required=True, help="Comma-separated PowerBin IDs, e.g. 231,328.")
    p.add_argument(
        "--top-n", type=int, default=8,
        help="Maximum surviving negative candidates traced per bin (default 8).",
    )
    p.add_argument(
        "--raw-negative-sigma", type=float, default=15.0,
        help=(
            "Minimum model-independent local negative significance for automatic candidate "
            "selection after atmospheric masking (default 15 robust-sigma)."
        ),
    )
    p.add_argument(
        "--candidate-separation-angstrom", type=float, default=2.5,
        help="Minimum rest-frame separation between automatically selected candidates (default 2.5 A).",
    )
    p.add_argument(
        "--local-sideband-inner-angstrom", type=float, default=4.0,
        help="Inner exclusion radius around a candidate when estimating local continuum (default 4 A).",
    )
    p.add_argument(
        "--local-sideband-outer-angstrom", type=float, default=22.0,
        help="Outer radius for local-continuum/scatter sidebands (default 22 A).",
    )
    p.add_argument(
        "--candidate-half-width-angstrom", type=float, default=1.5,
        help="Half-width used to summarize flux at one candidate in individual cubes (default 1.5 A).",
    )
    p.add_argument(
        "--window-half-width-angstrom", type=float, default=30.0,
        help="Half-width of per-candidate diagnostic plots in observed frame (default 30 A).",
    )
    p.add_argument(
        "--exposure-negative-z", type=float, default=8.0,
        help="Post-ZAP bin-median high-pass z threshold for an extreme negative exposure (default -8).",
    )
    p.add_argument(
        "--pre-sky-positive-z", type=float, default=8.0,
        help="Pre-ZAP bin-median high-pass z threshold for a strong positive sky feature (default +8).",
    )
    p.add_argument(
        "--spaxel-negative-z", type=float, default=5.0,
        help="Per-spaxel local negative threshold used for spatial-coherence fraction (default -5).",
    )
    p.add_argument(
        "--max-spaxel-traces", type=int, default=60,
        help="Maximum native spaxel traces drawn per exposure/panel (default 60).",
    )
    p.add_argument(
        "--cat-core-half-width-angstrom", type=float, default=5.0,
        help="CaT core half-width used for atmospheric-mask overlap audit (rest-frame air A; default 5).",
    )
    p.add_argument(
        "--cat-kinematic-half-width-angstrom", type=float, default=20.0,
        help="Broader CaT information-window half-width for mask-overlap audit (rest-frame air A; default 20).",
    )
    p.add_argument(
        "--output-dir", default=None,
        help="Default: <script3-run>/validation/03e_exposure_trace",
    )
    return p


def _load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def _write_json(payload: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, sort_keys=True)
        f.write("\n")


def _parse_bins(text: str, nbin: int) -> list[int]:
    out: list[int] = []
    for tok in str(text).split(","):
        tok = tok.strip()
        if not tok:
            continue
        bid = int(tok)
        if bid < 0 or bid >= nbin:
            raise ValueError(f"Bin {bid} is outside valid range 0..{nbin-1}.")
        if bid not in out:
            out.append(bid)
    if not out:
        raise ValueError("--bins did not contain any valid bin IDs.")
    return out


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
    raise ValueError(f"Unsupported wavelength medium conversion: {src!r} -> {dst!r}")


def _rest_template_to_observed_science(
    rest_wave: np.ndarray,
    redshift: float,
    template_medium: str,
    science_medium: str,
) -> np.ndarray:
    rest_science = _convert_medium(rest_wave, template_medium, science_medium)
    return rest_science * (1.0 + float(redshift))


def _resolve_script2_run(cfg, manifest: dict, explicit: str | None) -> Path:
    required_rel = Path("products") / "master_bins.fits"

    if explicit is not None:
        run = Path(explicit).expanduser().resolve()
        if not (run / required_rel).is_file():
            raise FileNotFoundError(f"Explicit Script-2 run lacks {required_rel}: {run}")
        return run

    raw = manifest.get("source_script2_run")
    if raw not in (None, ""):
        direct = Path(str(raw)).expanduser()
        if direct.is_dir() and (direct / required_rel).is_file():
            return direct.resolve()

        # A Script-3 run copied from Windows to macOS often retains a Windows
        # absolute path in its manifest.  Recover the run basename and look under
        # the current config RUNS_ROOT before forcing the user to repeat it.
        win_name = PureWindowsPath(str(raw)).name
        posix_name = Path(str(raw).replace("\\", "/")).name
        for name in (win_name, posix_name):
            if not name:
                continue
            candidate = Path(cfg.RUNS_ROOT).expanduser().resolve() / name
            if candidate.is_dir() and (candidate / required_rel).is_file():
                return candidate.resolve()

    raise FileNotFoundError(
        "Could not resolve the Script-2 run containing products/master_bins.fits. "
        "Pass --script2-run /Users/.../runs/<Script-2-run> explicitly."
    )


def _load_atmospheric_mask(
    path: Path | None,
    run3: Path,
    nlog: int,
    science_medium: str,
) -> tuple[np.ndarray, list[tuple[float, float]], Path | None]:
    if path is None:
        candidate = run3 / "validation" / "03d_atmosphere" / "products" / "RH3_03d_atmospheric_mask.ecsv"
        path = candidate if candidate.is_file() else None

    mask = np.zeros(nlog, dtype=bool)
    intervals: list[tuple[float, float]] = []
    if path is None:
        return mask, intervals, None

    path = Path(path).expanduser().resolve()
    tab = Table.read(path, format="ascii.ecsv")
    medium = str(tab.meta.get("OBSERVED_SCIENCE_WAVELENGTH_MEDIUM", "unknown")).strip().lower()
    if medium in {"air", "vacuum"} and science_medium in {"air", "vacuum"} and medium != science_medium:
        raise ValueError(
            f"Atmospheric mask is observed {medium}, but Script-3 science medium is {science_medium}."
        )

    for row in tab:
        if "INCLUDED_IN_ATMOSPHERIC_MASK" in tab.colnames and not bool(row["INCLUDED_IN_ATMOSPHERIC_MASK"]):
            continue
        i0 = int(row["LOG_INDEX_LO"])
        i1 = int(row["LOG_INDEX_HI"])
        if i0 < 0 or i1 >= nlog or i1 < i0:
            raise ValueError(f"Invalid mask log-index range {i0}..{i1} in {path}")
        mask[i0:i1 + 1] = True
        if "OBSERVED_SCIENCE_LO_ANGSTROM" in tab.colnames:
            intervals.append(
                (float(row["OBSERVED_SCIENCE_LO_ANGSTROM"]), float(row["OBSERVED_SCIENCE_HI_ANGSTROM"]))
            )
    return mask, intervals, path


def _local_raw_z(
    wave: np.ndarray,
    flux: np.ndarray,
    good: np.ndarray,
    inner: float,
    outer: float,
) -> np.ndarray:
    """Model-independent local deviation score at every wavelength sample."""
    w = np.asarray(wave, dtype=float)
    f = np.asarray(flux, dtype=float)
    g = np.asarray(good, dtype=bool) & np.isfinite(w) & np.isfinite(f)
    z = np.full(w.size, np.nan, dtype=float)

    for j in np.flatnonzero(g):
        d = np.abs(w - w[j])
        side = g & (d >= float(inner)) & (d <= float(outer))
        if np.count_nonzero(side) < 10:
            continue
        vals = f[side]
        med = float(np.median(vals))
        sig = float(1.4826 * np.median(np.abs(vals - med)))
        if not np.isfinite(sig) or sig <= 0:
            continue
        z[j] = (f[j] - med) / sig
    return z


def _select_candidates(
    wave: np.ndarray,
    flux: np.ndarray,
    good: np.ndarray,
    atmosphere_mask: np.ndarray,
    *,
    threshold: float,
    top_n: int,
    separation: float,
    inner: float,
    outer: float,
) -> tuple[np.ndarray, np.ndarray]:
    retained = np.asarray(good, dtype=bool) & ~np.asarray(atmosphere_mask, dtype=bool)
    z = _local_raw_z(wave, flux, retained, inner=inner, outer=outer)

    pool = np.flatnonzero(retained & np.isfinite(z) & (z <= -abs(float(threshold))))
    if pool.size == 0:
        # Never pretend there are no candidates merely because a threshold was
        # too conservative.  For diagnosis, fall back to the most negative local
        # deviations, while recording their actual significance.
        pool = np.flatnonzero(retained & np.isfinite(z))
        pool = pool[np.argsort(z[pool])[: max(1, int(top_n) * 3)]]
    else:
        pool = pool[np.argsort(z[pool])]

    chosen: list[int] = []
    for idx in pool:
        if any(abs(float(wave[idx]) - float(wave[k])) < float(separation) for k in chosen):
            continue
        chosen.append(int(idx))
        if len(chosen) >= int(top_n):
            break
    return np.asarray(chosen, dtype=int), z


def _frame_id(path: Path) -> str:
    m = re.search(r"_(\d{5})(?:_|\.|$)", path.name)
    return m.group(1) if m else path.stem


def _pair_cubes(pre_paths: list[Path], post_paths: list[Path]):
    if not post_paths:
        raise ValueError("At least one --post-zap-cube is required.")

    pre_by_id = {_frame_id(p): p for p in pre_paths}
    post_by_id = {_frame_id(p): p for p in post_paths}

    pairs = []
    for fid, post in sorted(post_by_id.items()):
        pairs.append((fid, pre_by_id.get(fid), post))

    unmatched_pre = sorted(set(pre_by_id) - set(post_by_id))
    if unmatched_pre:
        print(f"WARNING: pre-ZAP frames with no post-ZAP match will be ignored: {unmatched_pre}")

    if pre_paths and not any(pre is not None for _, pre, _ in pairs):
        if len(pre_paths) == len(post_paths):
            print("WARNING: frame IDs could not be matched; pairing pre/post cubes by supplied order.")
            pairs = [
                (_frame_id(post), pre, post)
                for pre, post in zip(pre_paths, post_paths)
            ]
        else:
            raise ValueError(
                "Pre/post frame IDs could not be matched and list lengths differ; "
                "supply matching exposure files."
            )
    return pairs


def _first_3d_hdu(hdul: fits.HDUList):
    for hdu in hdul:
        if hdu.data is not None and np.ndim(hdu.data) == 3:
            return hdu
    raise ValueError("No 3-D image HDU found.")


def _wavelength_axis(header: fits.Header, nwave: int) -> np.ndarray:
    if "CRVAL3" not in header or "CRPIX3" not in header:
        raise ValueError("Cube header lacks CRVAL3/CRPIX3 for a linear wavelength axis.")
    cdelt = header.get("CDELT3", header.get("CD3_3"))
    if cdelt is None:
        raise ValueError("Cube header lacks CDELT3/CD3_3.")
    pix = np.arange(nwave, dtype=float) + 1.0
    wave = float(header["CRVAL3"]) + (pix - float(header["CRPIX3"])) * float(cdelt)
    unit_text = str(header.get("CUNIT3", "Angstrom")).strip()
    try:
        unit = u.Unit(unit_text) if unit_text else u.AA
        wave = (wave * unit).to_value(u.AA)
    except Exception:
        # KCWI products conventionally use Angstrom; preserve values if a
        # non-standard unit string is present and make the assumption explicit.
        print(f"WARNING: could not parse CUNIT3={unit_text!r}; assuming wavelength values are Angstrom.")
    return np.asarray(wave, dtype=float)


class CubeView:
    def __init__(self, path: Path):
        self.path = Path(path).expanduser().resolve()
        with fits.open(self.path, memmap=True) as hdul:
            hdu = _first_3d_hdu(hdul)
            self.data = np.asarray(hdu.data, dtype=np.float32)
            # WCS may live in the image HDU or primary header.
            header = hdu.header.copy()
            if not any(k in header for k in ("CTYPE1", "CTYPE2")):
                header.extend(hdul[0].header, update=True)
            self.header = header

        if self.data.ndim != 3:
            raise ValueError(f"Expected a 3-D cube: {self.path}")
        self.wave = _wavelength_axis(self.header, self.data.shape[0])
        self.celestial = WCS(self.header).celestial
        scales = np.abs(proj_plane_pixel_scales(self.celestial)) * 3600.0
        self.pixel_scales_arcsec = np.asarray(scales, dtype=float)
        self.match_tolerance_arcsec = max(0.25, 0.80 * float(np.hypot(*self.pixel_scales_arcsec)))


def _bin_sky_coordinates(bin_map: np.ndarray, bin_wcs: WCS, bid: int) -> tuple[np.ndarray, np.ndarray, SkyCoord]:
    yy, xx = np.where(np.asarray(bin_map) == int(bid))
    if yy.size == 0:
        raise ValueError(f"Bin {bid} has no pixels in master BINMAP.")
    sky = bin_wcs.pixel_to_world(xx.astype(float), yy.astype(float))
    return yy, xx, sky


def _map_sky_to_cube(cube: CubeView, sky: SkyCoord):
    x, y = cube.celestial.world_to_pixel(sky)
    xi = np.rint(x).astype(int)
    yi = np.rint(y).astype(int)
    inside = (
        np.isfinite(x) & np.isfinite(y)
        & (xi >= 0) & (xi < cube.data.shape[2])
        & (yi >= 0) & (yi < cube.data.shape[1])
    )
    if not np.any(inside):
        return np.empty((0, 2), dtype=int), np.array([], dtype=float)

    xi = xi[inside]
    yi = yi[inside]
    sky_req = sky[inside]
    sky_pix = cube.celestial.pixel_to_world(xi.astype(float), yi.astype(float))
    dist = sky_req.separation(sky_pix).arcsec
    good = np.isfinite(dist) & (dist <= cube.match_tolerance_arcsec)

    coords = np.column_stack([yi[good], xi[good]])
    dist = np.asarray(dist[good], dtype=float)
    if coords.size == 0:
        return np.empty((0, 2), dtype=int), np.array([], dtype=float)

    # De-duplicate native spaxels; a coarse Large-slicer cube can map several
    # 0.3-arcsec BL output-grid samples to the same native spatial element.
    unique, first = np.unique(coords, axis=0, return_index=True)
    return unique, dist[first]


def _extract_spaxel_matrix(cube: CubeView, coords_yx: np.ndarray) -> np.ndarray:
    if coords_yx.size == 0:
        return np.empty((0, cube.data.shape[0]), dtype=float)
    rows = [cube.data[:, int(y), int(x)] for y, x in coords_yx]
    return np.asarray(rows, dtype=float)


def _local_measurements(
    wave: np.ndarray,
    spectra: np.ndarray,
    lam: float,
    *,
    center_half_width: float,
    side_inner: float,
    side_outer: float,
    spaxel_negative_z: float,
) -> dict:
    wave = np.asarray(wave, dtype=float)
    spectra = np.asarray(spectra, dtype=float)
    if spectra.ndim == 1:
        spectra = spectra[None, :]
    if spectra.size == 0:
        return {
            "n_spaxels": 0, "bin_delta": np.nan, "bin_z": np.nan,
            "bad_spaxel_fraction": np.nan, "median_spaxel_z": np.nan,
            "spaxel_z": np.array([], dtype=float), "spaxel_delta": np.array([], dtype=float),
        }

    center = np.abs(wave - float(lam)) <= float(center_half_width)
    d = np.abs(wave - float(lam))
    side = (d >= float(side_inner)) & (d <= float(side_outer))
    if np.count_nonzero(center) < 1 or np.count_nonzero(side) < 8:
        return {
            "n_spaxels": spectra.shape[0], "bin_delta": np.nan, "bin_z": np.nan,
            "bad_spaxel_fraction": np.nan, "median_spaxel_z": np.nan,
            "spaxel_z": np.full(spectra.shape[0], np.nan),
            "spaxel_delta": np.full(spectra.shape[0], np.nan),
        }

    zvals = np.full(spectra.shape[0], np.nan, dtype=float)
    delta = np.full(spectra.shape[0], np.nan, dtype=float)
    for i, spec in enumerate(spectra):
        sidevals = spec[side & np.isfinite(spec)]
        centervals = spec[center & np.isfinite(spec)]
        if sidevals.size < 8 or centervals.size == 0:
            continue
        med = float(np.median(sidevals))
        sig = float(1.4826 * np.median(np.abs(sidevals - med)))
        if not np.isfinite(sig) or sig <= 0:
            continue
        c = float(np.median(centervals))
        delta[i] = c - med
        zvals[i] = delta[i] / sig

    medspec = np.nanmedian(spectra, axis=0)
    sidevals = medspec[side & np.isfinite(medspec)]
    centervals = medspec[center & np.isfinite(medspec)]
    if sidevals.size >= 8 and centervals.size:
        bmed = float(np.median(sidevals))
        bsig = float(1.4826 * np.median(np.abs(sidevals - bmed)))
        bdelta = float(np.median(centervals) - bmed)
        bz = bdelta / bsig if np.isfinite(bsig) and bsig > 0 else np.nan
    else:
        bdelta = np.nan
        bz = np.nan

    finite_z = np.isfinite(zvals)
    badfrac = (
        float(np.mean(zvals[finite_z] <= -abs(float(spaxel_negative_z))))
        if np.any(finite_z) else np.nan
    )
    return {
        "n_spaxels": int(spectra.shape[0]),
        "bin_delta": bdelta,
        "bin_z": bz,
        "bad_spaxel_fraction": badfrac,
        "median_spaxel_z": float(np.nanmedian(zvals)) if np.any(finite_z) else np.nan,
        "spaxel_z": zvals,
        "spaxel_delta": delta,
    }


def _normalize_window(wave: np.ndarray, spectra: np.ndarray, lam: float, inner: float, outer: float):
    """Normalize each spaxel by its local sideband median for diagnostic plotting."""
    wave = np.asarray(wave, dtype=float)
    s = np.asarray(spectra, dtype=float)
    if s.ndim == 1:
        s = s[None, :]
    d = np.abs(wave - float(lam))
    side = (d >= float(inner)) & (d <= float(outer))
    out = s.copy()
    for i in range(s.shape[0]):
        vals = s[i, side & np.isfinite(s[i])]
        if vals.size:
            scale = float(np.median(vals))
            if np.isfinite(scale) and scale != 0:
                out[i] = s[i] / scale
    return out


def _interval_overlap_fraction(intervals: list[tuple[float, float]], lo: float, hi: float) -> float:
    if hi <= lo:
        return 0.0
    pieces = []
    for a, b in intervals:
        aa = max(float(lo), float(a))
        bb = min(float(hi), float(b))
        if bb > aa:
            pieces.append((aa, bb))
    if not pieces:
        return 0.0
    pieces.sort()
    merged = [list(pieces[0])]
    for a, b in pieces[1:]:
        if a <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], b)
        else:
            merged.append([a, b])
    covered = sum(b - a for a, b in merged)
    return float(covered / (hi - lo))


def _cat_audit(
    cfg,
    manifest: dict,
    intervals: list[tuple[float, float]],
    *,
    core_half_width: float,
    kin_half_width: float,
) -> Table:
    science_medium = str(manifest.get("science_wavelength_medium", "vacuum")).strip().lower()
    template_medium = str(manifest.get("template_wavelength_medium", "air")).strip().lower()
    z = float(manifest.get("redshift", cfg.REDSHIFT))
    fit_lo, fit_hi = [float(x) for x in manifest.get("fit_rest_range_angstrom", cfg.RH3_FIT_REST_RANGE_ANGSTROM)]

    rows = []
    for name, rest_air in zip(CAT_NAMES, CAT_REST_AIR_ANGSTROM):
        rest_template = float(_convert_medium(np.array([rest_air]), "air", template_medium)[0])
        obs_center = float(
            _rest_template_to_observed_science(
                np.array([rest_template]), z, template_medium, science_medium
            )[0]
        )
        center_masked = any(lo <= obs_center <= hi for lo, hi in intervals)

        def obs_window(half):
            lo_rest_air = rest_air - float(half)
            hi_rest_air = rest_air + float(half)
            rest_t = _convert_medium(np.array([lo_rest_air, hi_rest_air]), "air", template_medium)
            obs = _rest_template_to_observed_science(rest_t, z, template_medium, science_medium)
            return float(np.min(obs)), float(np.max(obs))

        core_lo, core_hi = obs_window(core_half_width)
        kin_lo, kin_hi = obs_window(kin_half_width)
        core_frac = _interval_overlap_fraction(intervals, core_lo, core_hi)
        kin_frac = _interval_overlap_fraction(intervals, kin_lo, kin_hi)

        in_fit = fit_lo <= rest_template <= fit_hi
        if not in_fit:
            status = "NOT_IN_CURRENT_SCRIPT3_FIT"
        elif center_masked:
            status = "LINE_CENTER_MASKED"
        elif kin_frac >= 0.25:
            status = "SUBSTANTIAL_KINEMATIC_WINDOW_MASK"
        elif kin_frac > 0:
            status = "PARTIAL_KINEMATIC_WINDOW_MASK"
        else:
            status = "CLEAR"

        rows.append((
            name, rest_air, rest_template, obs_center, in_fit, center_masked,
            core_lo, core_hi, core_frac, kin_lo, kin_hi, kin_frac, status,
        ))

    tab = Table(
        rows=rows,
        names=(
            "LINE", "REST_AIR_ANGSTROM", "REST_TEMPLATE_ANGSTROM",
            "OBSERVED_SCIENCE_CENTER_ANGSTROM", "IN_CURRENT_SCRIPT3_FIT",
            "LINE_CENTER_MASKED", "CORE_OBS_LO_ANGSTROM", "CORE_OBS_HI_ANGSTROM",
            "CORE_MASK_FRACTION", "KIN_OBS_LO_ANGSTROM", "KIN_OBS_HI_ANGSTROM",
            "KINEMATIC_WINDOW_MASK_FRACTION", "STATUS",
        ),
    )
    tab.meta["SCIENCE_WAVELENGTH_MEDIUM"] = science_medium
    tab.meta["TEMPLATE_WAVELENGTH_MEDIUM"] = template_medium
    tab.meta["REDSHIFT"] = z
    tab.meta["FIT_REST_RANGE_ANGSTROM"] = [fit_lo, fit_hi]
    tab.meta["CAT_CORE_HALF_WIDTH_REST_AIR_ANGSTROM"] = float(core_half_width)
    tab.meta["CAT_KINEMATIC_HALF_WIDTH_REST_AIR_ANGSTROM"] = float(kin_half_width)
    return tab


def _candidate_classification(
    exposure_rows: list[dict],
    *,
    exposure_negative_z: float,
    pre_sky_positive_z: float,
) -> tuple[str, str]:
    if not exposure_rows:
        return "NO_EXPOSURE_DATA", "UNRESOLVED"

    post_extreme = [
        np.isfinite(r["POST_BIN_Z"]) and r["POST_BIN_Z"] <= -abs(float(exposure_negative_z))
        for r in exposure_rows
    ]
    pre_positive = [
        np.isfinite(r["PRE_BIN_Z"]) and r["PRE_BIN_Z"] >= abs(float(pre_sky_positive_z))
        for r in exposure_rows
    ]
    n = len(exposure_rows)
    npost = int(np.sum(post_extreme))
    npre = int(np.sum(pre_positive))

    if npost == 0:
        recurrence = "NOT_REPRODUCED_IN_INDIVIDUAL_POST_ZAP"
    elif npost == 1:
        recurrence = "SINGLE_EXPOSURE"
    elif npost >= math.ceil(0.75 * n):
        recurrence = "RECURRENT_ACROSS_EXPOSURES"
    else:
        recurrence = "MULTI_EXPOSURE_PARTIAL"

    bad_fracs = np.array([r["POST_BAD_SPAXEL_FRACTION"] for r in exposure_rows], dtype=float)
    med_bad = float(np.nanmedian(bad_fracs)) if np.any(np.isfinite(bad_fracs)) else np.nan

    paired_over = sum(bool(pe and po) for pe, po in zip(pre_positive, post_extreme))
    if npost == 0:
        origin = "STACK_REPROJECTION_OR_COADD_SUSPECT"
    elif paired_over >= max(1, math.ceil(0.5 * npost)):
        origin = "SKY_OVER_SUBTRACTION_SUSPECT"
    elif npost == 1:
        origin = "EXPOSURE_SPECIFIC_ARTIFACT_SUSPECT"
    elif recurrence == "RECURRENT_ACROSS_EXPOSURES" and np.isfinite(med_bad) and med_bad >= 0.50:
        origin = "RECURRENT_SPATIALLY_COHERENT_SYSTEMATIC_SUSPECT"
    elif recurrence == "RECURRENT_ACROSS_EXPOSURES" and np.isfinite(med_bad) and med_bad <= 0.20:
        origin = "RECURRENT_SPATIALLY_LOCALIZED_DEFECT_SUSPECT"
    elif npre > 0:
        origin = "ATMOSPHERIC_OR_SKY_SUBTRACTION_RELATED_SUSPECT"
    else:
        origin = "UNRESOLVED"

    return recurrence, origin


def _plot_candidate(
    outpath: Path,
    *,
    bid: int,
    candidate_rank: int,
    rest_wave: float,
    observed_wave: float,
    script3_z: float,
    pairs,
    mapped,
    window_half_width: float,
    side_inner: float,
    side_outer: float,
    max_spaxel_traces: int,
):
    nrow = len(pairs)
    ncol = 2 if any(pre is not None for _, pre, _ in pairs) else 1
    fig, axes = plt.subplots(
        nrow, ncol,
        figsize=(13 if ncol == 2 else 8, max(3.0, 2.7 * nrow)),
        squeeze=False, sharex=False,
    )

    for i, (fid, pre_cube, post_cube) in enumerate(pairs):
        datasets = []
        if ncol == 2:
            datasets.append(("pre-ZAP", pre_cube, mapped[(fid, "pre")]))
        datasets.append(("post-ZAP", post_cube, mapped[(fid, "post")]))

        for j, (label, cube, coords) in enumerate(datasets):
            ax = axes[i, j]
            if cube is None or coords.size == 0:
                ax.text(0.5, 0.5, "no usable mapped spaxels", ha="center", va="center", transform=ax.transAxes)
                ax.set_title(f"{fid} {label}")
                continue

            matrix = _extract_spaxel_matrix(cube, coords)
            win = np.abs(cube.wave - float(observed_wave)) <= float(window_half_width)
            if not np.any(win):
                ax.text(0.5, 0.5, "candidate outside cube", ha="center", va="center", transform=ax.transAxes)
                ax.set_title(f"{fid} {label}")
                continue

            norm = _normalize_window(cube.wave, matrix, observed_wave, side_inner, side_outer)
            draw = norm[: int(max_spaxel_traces)]
            for row in draw:
                ax.plot(cube.wave[win], row[win], lw=0.45, alpha=0.18)
            med = np.nanmedian(norm, axis=0)
            ax.plot(cube.wave[win], med[win], lw=1.4, label="native-spaxel median")
            ax.axvline(observed_wave, ls="--", lw=0.8)
            ax.set_title(f"frame {fid} | {label} | mapped native spaxels={matrix.shape[0]}")
            ax.set_ylabel("local-normalized flux")
            ax.legend(fontsize=7, loc="best")

    for ax in axes[-1]:
        ax.set_xlabel(f"Observed wavelength (A) | candidate={observed_wave:.2f} A")

    fig.suptitle(
        f"Script 03e bin {bid} candidate {candidate_rank} | "
        f"rest={rest_wave:.2f} A | Script-3 local raw z={script3_z:.1f}",
        fontsize=11,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(outpath, dpi=160)
    plt.close(fig)


def _plot_bin_overview(
    outpath: Path,
    *,
    bid: int,
    wave_rest: np.ndarray,
    galaxy: np.ndarray,
    good: np.ndarray,
    atmospheric_mask: np.ndarray,
    candidates: np.ndarray,
    candidate_z: np.ndarray,
    classifications: list[dict],
    exposure_ids: list[str],
):
    fig, axes = plt.subplots(2, 1, figsize=(13, 7), sharex=False)

    ax = axes[0]
    ax.plot(wave_rest, galaxy, lw=0.8, label="saved Script-3 normalized spectrum")
    masked = good & atmospheric_mask & np.isfinite(galaxy)
    if np.any(masked):
        ax.scatter(wave_rest[masked], galaxy[masked], marker="x", s=15, label="03d atmospheric mask")
    if candidates.size:
        ax.scatter(
            wave_rest[candidates], galaxy[candidates], marker="o", s=34,
            facecolors="none", label="03e surviving candidate"
        )
    ax.set_xlabel("Rest-frame wavelength (A)")
    ax.set_ylabel("Normalized flux")
    ax.set_title(f"Bin {bid}: surviving model-independent negative candidates")
    ax.legend(fontsize=8)

    ax = axes[1]
    if classifications:
        matrix = np.full((len(exposure_ids), len(classifications)), np.nan, dtype=float)
        labels = []
        for j, row in enumerate(classifications):
            labels.append(f"{row['REST_WAVE_ANGSTROM']:.1f}")
            by_exp = row["EXPOSURE_Z_BY_ID"]
            for i, fid in enumerate(exposure_ids):
                matrix[i, j] = by_exp.get(fid, np.nan)
        im = ax.imshow(matrix, origin="upper", aspect="auto")
        ax.set_yticks(np.arange(len(exposure_ids)), labels=exposure_ids)
        ax.set_xticks(np.arange(len(labels)), labels=labels, rotation=45, ha="right")
        ax.set_xlabel("Candidate rest wavelength (A)")
        ax.set_ylabel("post-ZAP exposure")
        fig.colorbar(im, ax=ax, label="post-ZAP bin-median local z")
    else:
        ax.text(0.5, 0.5, "No candidates", ha="center", va="center", transform=ax.transAxes)

    fig.tight_layout()
    fig.savefig(outpath, dpi=160)
    plt.close(fig)


def main() -> int:
    args = _parser().parse_args()
    run3 = Path(args.script3_run).expanduser().resolve()
    cfg = crd.load_config(args.config, validate=True, strict_paths=False)
    manifest = _load_json(run3 / "metadata" / "script03_manifest.json")

    out = (
        Path(args.output_dir).expanduser().resolve()
        if args.output_dir is not None
        else run3 / "validation" / "03e_exposure_trace"
    )
    products = out / "products"
    figures = out / "figures"
    metadata = out / "metadata"
    for p in (products, figures, metadata):
        p.mkdir(parents=True, exist_ok=True)

    with np.load(run3 / "products" / "RH3_log_spectra_and_local_best_fits.npz", allow_pickle=False) as z:
        spectra = {k: np.asarray(z[k]) for k in z.files}

    wave_rest = np.asarray(spectra["wavelength"], dtype=float)
    galaxy_all = np.asarray(spectra["galaxy"], dtype=float)
    noise_all = np.asarray(spectra["noise"], dtype=float)
    good_all = np.asarray(spectra["good"], dtype=bool)

    nbin = galaxy_all.shape[0]
    bins = _parse_bins(args.bins, nbin)

    science_medium = str(manifest.get("science_wavelength_medium", "vacuum")).strip().lower()
    template_medium = str(manifest.get("template_wavelength_medium", cfg.TEMPLATE_WAVELENGTH_MEDIUM)).strip().lower()
    redshift = float(manifest.get("redshift", cfg.REDSHIFT))
    observed_grid = _rest_template_to_observed_science(
        wave_rest, redshift, template_medium, science_medium
    )

    mask_arg = Path(args.atmospheric_mask).expanduser().resolve() if args.atmospheric_mask else None
    atmospheric_mask, atmospheric_intervals, atmospheric_path = _load_atmospheric_mask(
        mask_arg, run3, wave_rest.size, science_medium
    )

    run2 = _resolve_script2_run(cfg, manifest, args.script2_run)
    with fits.open(run2 / "products" / "master_bins.fits") as hdul:
        bin_map = np.asarray(hdul[0].data, dtype=int)
        bin_wcs = WCS(hdul[0].header).celestial

    pre_paths = [Path(p).expanduser().resolve() for p in args.pre_zap_cube]
    post_paths = [Path(p).expanduser().resolve() for p in args.post_zap_cube]
    pairs_paths = _pair_cubes(pre_paths, post_paths)

    pairs = []
    for fid, pre_path, post_path in pairs_paths:
        pre_cube = CubeView(pre_path) if pre_path is not None else None
        post_cube = CubeView(post_path)
        pairs.append((fid, pre_cube, post_cube))

    # CaT audit first because preserving kinematic information is a core science
    # requirement, not merely an after-the-fact plotting concern.
    cat = _cat_audit(
        cfg, manifest, atmospheric_intervals,
        core_half_width=float(args.cat_core_half_width_angstrom),
        kin_half_width=float(args.cat_kinematic_half_width_angstrom),
    )
    cat.meta["ATMOSPHERIC_MASK_SOURCE"] = str(atmospheric_path) if atmospheric_path else "none"
    cat_path = products / "RH3_03e_CaT_mask_audit.ecsv"
    cat.write(cat_path, format="ascii.ecsv", overwrite=True)

    print("\nCa II triplet atmospheric-mask audit:")
    for row in cat:
        print(
            f"  {row['LINE']}: rest_air={row['REST_AIR_ANGSTROM']:.2f} A, "
            f"observed_{science_medium}={row['OBSERVED_SCIENCE_CENTER_ANGSTROM']:.2f} A, "
            f"in_current_fit={bool(row['IN_CURRENT_SCRIPT3_FIT'])}, "
            f"center_masked={bool(row['LINE_CENTER_MASKED'])}, "
            f"+/-20A mask fraction={row['KINEMATIC_WINDOW_MASK_FRACTION']:.3f}, "
            f"status={row['STATUS']}"
        )

    exposure_rows: list[dict] = []
    candidate_rows: list[dict] = []

    for bid in bins:
        _, _, bin_sky = _bin_sky_coordinates(bin_map, bin_wcs, bid)

        # Map this physical PowerBin independently to every pre/post-ZAP cube.
        mapped = {}
        for fid, pre_cube, post_cube in pairs:
            if pre_cube is not None:
                coords, dist = _map_sky_to_cube(pre_cube, bin_sky)
                mapped[(fid, "pre")] = coords
                mapped[(fid, "pre_dist")] = dist
            else:
                mapped[(fid, "pre")] = np.empty((0, 2), dtype=int)
                mapped[(fid, "pre_dist")] = np.array([], dtype=float)

            coords, dist = _map_sky_to_cube(post_cube, bin_sky)
            mapped[(fid, "post")] = coords
            mapped[(fid, "post_dist")] = dist

        candidate_idx, raw_z = _select_candidates(
            wave_rest,
            galaxy_all[bid],
            good_all[bid],
            atmospheric_mask,
            threshold=float(args.raw_negative_sigma),
            top_n=int(args.top_n),
            separation=float(args.candidate_separation_angstrom),
            inner=float(args.local_sideband_inner_angstrom),
            outer=float(args.local_sideband_outer_angstrom),
        )

        print(f"\nBin {bid}: selected {candidate_idx.size} surviving candidates after 03d atmospheric masking.")

        bin_classifications = []
        for rank, idx in enumerate(candidate_idx, start=1):
            rest_lam = float(wave_rest[idx])
            obs_lam = float(observed_grid[idx])
            rows_this_candidate = []

            print(
                f"  candidate {rank}: rest={rest_lam:.2f} A -> observed {science_medium}={obs_lam:.2f} A "
                f"| local raw z={raw_z[idx]:.1f}"
            )

            for fid, pre_cube, post_cube in pairs:
                pre_coords = mapped[(fid, "pre")]
                post_coords = mapped[(fid, "post")]

                pre_matrix = _extract_spaxel_matrix(pre_cube, pre_coords) if pre_cube is not None else np.empty((0, 0))
                post_matrix = _extract_spaxel_matrix(post_cube, post_coords)

                pre_meas = (
                    _local_measurements(
                        pre_cube.wave, pre_matrix, obs_lam,
                        center_half_width=float(args.candidate_half_width_angstrom),
                        side_inner=float(args.local_sideband_inner_angstrom),
                        side_outer=float(args.local_sideband_outer_angstrom),
                        spaxel_negative_z=float(args.spaxel_negative_z),
                    )
                    if pre_cube is not None else
                    {
                        "n_spaxels": 0, "bin_delta": np.nan, "bin_z": np.nan,
                        "bad_spaxel_fraction": np.nan, "median_spaxel_z": np.nan,
                    }
                )
                post_meas = _local_measurements(
                    post_cube.wave, post_matrix, obs_lam,
                    center_half_width=float(args.candidate_half_width_angstrom),
                    side_inner=float(args.local_sideband_inner_angstrom),
                    side_outer=float(args.local_sideband_outer_angstrom),
                    spaxel_negative_z=float(args.spaxel_negative_z),
                )

                pre_dist = mapped[(fid, "pre_dist")]
                post_dist = mapped[(fid, "post_dist")]

                row = {
                    "BIN_ID": int(bid),
                    "CANDIDATE_RANK": int(rank),
                    "LOG_INDEX": int(idx),
                    "REST_WAVE_ANGSTROM": rest_lam,
                    "OBSERVED_SCIENCE_WAVE_ANGSTROM": obs_lam,
                    "SCRIPT3_LOCAL_RAW_Z": float(raw_z[idx]),
                    "EXPOSURE_ID": str(fid),
                    "PRE_CUBE": str(pre_cube.path) if pre_cube is not None else "",
                    "POST_CUBE": str(post_cube.path),
                    "N_PRE_MAPPED_SPAXELS": int(pre_meas["n_spaxels"]),
                    "N_POST_MAPPED_SPAXELS": int(post_meas["n_spaxels"]),
                    "PRE_MEDIAN_MAP_DISTANCE_ARCSEC": float(np.nanmedian(pre_dist)) if pre_dist.size else np.nan,
                    "POST_MEDIAN_MAP_DISTANCE_ARCSEC": float(np.nanmedian(post_dist)) if post_dist.size else np.nan,
                    "PRE_BIN_DELTA": float(pre_meas["bin_delta"]),
                    "PRE_BIN_Z": float(pre_meas["bin_z"]),
                    "POST_BIN_DELTA": float(post_meas["bin_delta"]),
                    "POST_BIN_Z": float(post_meas["bin_z"]),
                    "POST_BAD_SPAXEL_FRACTION": float(post_meas["bad_spaxel_fraction"]),
                    "POST_MEDIAN_SPAXEL_Z": float(post_meas["median_spaxel_z"]),
                }
                exposure_rows.append(row)
                rows_this_candidate.append(row)

            recurrence, origin = _candidate_classification(
                rows_this_candidate,
                exposure_negative_z=float(args.exposure_negative_z),
                pre_sky_positive_z=float(args.pre_sky_positive_z),
            )

            z_by_id = {r["EXPOSURE_ID"]: r["POST_BIN_Z"] for r in rows_this_candidate}
            bad_fracs = np.array([r["POST_BAD_SPAXEL_FRACTION"] for r in rows_this_candidate], dtype=float)
            med_bad = float(np.nanmedian(bad_fracs)) if np.any(np.isfinite(bad_fracs)) else np.nan

            c_row = {
                "BIN_ID": int(bid),
                "CANDIDATE_RANK": int(rank),
                "LOG_INDEX": int(idx),
                "REST_WAVE_ANGSTROM": rest_lam,
                "OBSERVED_SCIENCE_WAVE_ANGSTROM": obs_lam,
                "SCRIPT3_LOCAL_RAW_Z": float(raw_z[idx]),
                "N_EXPOSURES": len(rows_this_candidate),
                "N_POST_EXTREME_NEGATIVE": int(sum(
                    np.isfinite(r["POST_BIN_Z"]) and r["POST_BIN_Z"] <= -abs(float(args.exposure_negative_z))
                    for r in rows_this_candidate
                )),
                "N_PRE_STRONG_POSITIVE": int(sum(
                    np.isfinite(r["PRE_BIN_Z"]) and r["PRE_BIN_Z"] >= abs(float(args.pre_sky_positive_z))
                    for r in rows_this_candidate
                )),
                "MEDIAN_POST_BAD_SPAXEL_FRACTION": med_bad,
                "EXPOSURE_RECURRENCE_CLASS": recurrence,
                "ORIGIN_SUSPECT_CLASS": origin,
                "EXPOSURE_Z_BY_ID": z_by_id,  # plotting-only; not written directly to ECSV
            }
            candidate_rows.append(c_row)
            bin_classifications.append(c_row)

            _plot_candidate(
                figures / f"bin_{bid:04d}_candidate_{rank:02d}_{obs_lam:.1f}A.png",
                bid=bid,
                candidate_rank=rank,
                rest_wave=rest_lam,
                observed_wave=obs_lam,
                script3_z=float(raw_z[idx]),
                pairs=pairs,
                mapped=mapped,
                window_half_width=float(args.window_half_width_angstrom),
                side_inner=float(args.local_sideband_inner_angstrom),
                side_outer=float(args.local_sideband_outer_angstrom),
                max_spaxel_traces=int(args.max_spaxel_traces),
            )

        _plot_bin_overview(
            figures / f"bin_{bid:04d}_exposure_trace_overview.png",
            bid=bid,
            wave_rest=wave_rest,
            galaxy=galaxy_all[bid],
            good=good_all[bid],
            atmospheric_mask=atmospheric_mask,
            candidates=candidate_idx,
            candidate_z=raw_z,
            classifications=bin_classifications,
            exposure_ids=[fid for fid, _, _ in pairs],
        )

    # ECSV products
    exposure_names = [
        "BIN_ID", "CANDIDATE_RANK", "LOG_INDEX", "REST_WAVE_ANGSTROM",
        "OBSERVED_SCIENCE_WAVE_ANGSTROM", "SCRIPT3_LOCAL_RAW_Z", "EXPOSURE_ID",
        "PRE_CUBE", "POST_CUBE", "N_PRE_MAPPED_SPAXELS", "N_POST_MAPPED_SPAXELS",
        "PRE_MEDIAN_MAP_DISTANCE_ARCSEC", "POST_MEDIAN_MAP_DISTANCE_ARCSEC",
        "PRE_BIN_DELTA", "PRE_BIN_Z", "POST_BIN_DELTA", "POST_BIN_Z",
        "POST_BAD_SPAXEL_FRACTION", "POST_MEDIAN_SPAXEL_Z",
    ]
    exp_table = Table(rows=[[r[n] for n in exposure_names] for r in exposure_rows], names=exposure_names)
    exp_table.meta["SOURCE_SCRIPT3_RUN"] = str(run3)
    exp_table.meta["SOURCE_SCRIPT2_RUN"] = str(run2)
    exp_table.meta["ATMOSPHERIC_MASK"] = str(atmospheric_path) if atmospheric_path else "none"
    exp_table.meta["SCIENCE_WAVELENGTH_MEDIUM"] = science_medium
    exp_table.write(products / "RH3_03e_exposure_trace.ecsv", format="ascii.ecsv", overwrite=True)

    cand_names = [
        "BIN_ID", "CANDIDATE_RANK", "LOG_INDEX", "REST_WAVE_ANGSTROM",
        "OBSERVED_SCIENCE_WAVE_ANGSTROM", "SCRIPT3_LOCAL_RAW_Z",
        "N_EXPOSURES", "N_POST_EXTREME_NEGATIVE", "N_PRE_STRONG_POSITIVE",
        "MEDIAN_POST_BAD_SPAXEL_FRACTION", "EXPOSURE_RECURRENCE_CLASS",
        "ORIGIN_SUSPECT_CLASS",
    ]
    cand_table = Table(rows=[[r[n] for n in cand_names] for r in candidate_rows], names=cand_names)
    cand_table.meta["DIAGNOSTIC_ONLY"] = True
    cand_table.meta["NOTE"] = (
        "Origin classes are hypotheses from pre/post-ZAP recurrence and spatial coherence. "
        "Script 03e never turns these classes into a science mask automatically."
    )
    cand_table.write(products / "RH3_03e_candidate_classification.ecsv", format="ascii.ecsv", overwrite=True)

    manifest_out = {
        "script": SCRIPT_NAME,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source_script3_run": str(run3),
        "source_script2_run": str(run2),
        "config": str(Path(args.config).expanduser().resolve()),
        "atmospheric_mask": str(atmospheric_path) if atmospheric_path else None,
        "science_wavelength_medium": science_medium,
        "template_wavelength_medium": template_medium,
        "redshift": redshift,
        "bins": bins,
        "pre_post_pairs": [
            {
                "exposure_id": fid,
                "pre_zap_cube": str(pre.path) if pre is not None else None,
                "post_zap_cube": str(post.path),
            }
            for fid, pre, post in pairs
        ],
        "candidate_selection": {
            "raw_negative_sigma": float(args.raw_negative_sigma),
            "top_n": int(args.top_n),
            "minimum_separation_rest_angstrom": float(args.candidate_separation_angstrom),
            "selection_is_after_atmospheric_mask": True,
            "selection_is_model_independent": True,
        },
        "classification_thresholds": {
            "post_exposure_negative_z": -abs(float(args.exposure_negative_z)),
            "pre_sky_positive_z": abs(float(args.pre_sky_positive_z)),
            "post_spaxel_negative_z": -abs(float(args.spaxel_negative_z)),
        },
        "products": {
            "exposure_trace_ecsv": str(products / "RH3_03e_exposure_trace.ecsv"),
            "candidate_classification_ecsv": str(products / "RH3_03e_candidate_classification.ecsv"),
            "cat_mask_audit_ecsv": str(cat_path),
            "figures_dir": str(figures),
        },
    }
    _write_json(manifest_out, metadata / "script03e_manifest.json")

    print(f"\nScript 03e complete: {out}")
    print(f"  Exposure trace: {products / 'RH3_03e_exposure_trace.ecsv'}")
    print(f"  Candidate classes: {products / 'RH3_03e_candidate_classification.ecsv'}")
    print(f"  CaT audit: {cat_path}")
    print("  IMPORTANT: 03e classifications are diagnostic only and do not modify Script-3 goodpixels.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
