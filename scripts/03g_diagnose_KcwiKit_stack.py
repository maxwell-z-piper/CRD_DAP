#!/usr/bin/env python3
"""CRD_DAP Script 03g: forensic diagnosis of KcwiKit RH3 stack artifacts.

03f established that the catastrophic RH3 samples already exist in the final
KcwiKit stack.  This diagnostic therefore asks whether a single KSkyWizard
exposure, a nearby input voxel, or the KcwiKit reprojection/stack is responsible.

DIAGNOSTIC ONLY: no masks or science products are modified.

Typical Mac usage
-----------------
python scripts/03g_diagnose_KcwiKit_stack.py \
    --script3-run /Users/maxpiper/CRD_DAP/runs/8143-1902_S03_20260828_122251 \
    --config config/8143-1902.py \
    --post-zap-cube /Volumes/Mainframe/KCWI/KcwiKit_Redux/redux/kskywizard/kr241226_00172_zap_icubes.fits \
    --post-zap-cube /Volumes/Mainframe/KCWI/KcwiKit_Redux/redux/kskywizard/kr241226_00173_zap_icubes.fits \
    --post-zap-cube /Volumes/Mainframe/KCWI/KcwiKit_Redux/redux/kskywizard/kr241226_00174_zap_icubes.fits \
    --post-zap-cube /Volumes/Mainframe/KCWI/KcwiKit_Redux/redux/kskywizard/kr241226_00175_zap_icubes.fits

Alternate final-stack experiment
--------------------------------
Add all four explicit final-stack paths plus a label, for example:

    --final-icube /path/red_bilinear_icubes.fits \
    --final-vcube /path/red_bilinear_vcubes.fits \
    --final-mcube /path/red_bilinear_mcubes.fits \
    --final-ecube /path/red_bilinear_ecubes.fits \
    --stack-label bilinear_all4

The tested stack must preserve the original Script-2 output grid. 03g validates
shape, spectral centers, and celestial WCS before applying the saved PowerBin
membership raster.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from itertools import combinations
import json
import math
from pathlib import Path, PureWindowsPath
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS

import crd_utils as crd
from crd_utils import io


def parser():
    p = argparse.ArgumentParser(description="Trace bad final KcwiKit RH3 voxels into individual post-ZAP exposures.")
    p.add_argument("--script3-run", required=True)
    p.add_argument("--config", required=True)
    p.add_argument("--script2-run", default=None)
    p.add_argument("--script1-run", default=None)
    p.add_argument("--trace-table", default=None, help="03f stage-transition ECSV; auto-detected by default.")
    p.add_argument("--post-zap-cube", action="append", default=[], required=True)

    # By default 03g uses the four KcwiKit stack files recorded by Script 01.
    # These explicit overrides make controlled alternate-stack experiments
    # possible without rerunning Scripts 01--03 or editing provenance manifests.
    p.add_argument("--final-icube", default=None,
                   help="Explicit alternate KcwiKit science stack (*_icubes.fits).")
    p.add_argument("--final-vcube", default=None,
                   help="Explicit alternate KcwiKit variance stack (*_vcubes.fits).")
    p.add_argument("--final-mcube", default=None,
                   help="Explicit alternate KcwiKit mask stack (*_mcubes.fits).")
    p.add_argument("--final-ecube", default=None,
                   help="Explicit alternate KcwiKit effective-exposure stack (*_ecubes.fits).")
    p.add_argument("--stack-label", default=None,
                   help=("Short label for this tested stack, e.g. bilinear_all4, omit173, "
                         "single173. Used in output paths/metadata."))
    p.add_argument("--grid-match-tolerance-arcsec", type=float, default=0.05,
                   help=("Maximum celestial-WCS mismatch allowed relative to the original "
                         "Script-1 KcwiKit stack at test grid points (default 0.05 arcsec). "
                         "Alternate stacks must use the same final output grid so the saved "
                         "Script-2 membership raster remains valid."))

    p.add_argument("--bins", default=None)
    p.add_argument("--max-output-spaxels", type=int, default=4)
    p.add_argument("--spatial-radius", type=int, default=2, help="Native input-spaxel radius; 2 gives a 5x5 search.")
    p.add_argument("--spectral-radius", type=int, default=4, help="Search +/- this many native input channels.")
    p.add_argument("--input-extreme-z", type=float, default=8.0)
    p.add_argument("--final-pathology-z", type=float, default=20.0)
    p.add_argument("--sideband-inner", type=float, default=4.0)
    p.add_argument("--sideband-outer", type=float, default=18.0)
    p.add_argument("--plot-half-width", type=float, default=22.0)
    p.add_argument("--exposure-time-tolerance", type=float, default=2.0)
    p.add_argument("--thum0", default=None)
    p.add_argument("--thum", default=None)
    p.add_argument("--output-dir", default=None)
    return p


def read_json(path):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def write_json(obj, path):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2, sort_keys=True, default=str)
        f.write("\n")


def any_basename(value):
    s = str(value)
    return PureWindowsPath(s).name or Path(s.replace("\\", "/")).name


def resolve_run(explicit, manifest_value, runs_root, required_rel, flag):
    cand = []
    if explicit:
        cand.append(Path(explicit).expanduser())
    if manifest_value:
        cand.append(Path(str(manifest_value)).expanduser())
        name = any_basename(manifest_value)
        if name:
            cand.append(Path(runs_root) / name)
    for c in cand:
        try:
            c = c.resolve()
        except Exception:
            pass
        if c.is_dir() and (c / required_rel).is_file():
            return c
    raise FileNotFoundError(f"Could not resolve run containing {required_rel}; pass {flag} explicitly.")


def frame_id(path):
    m = re.search(r"_(\d{5})(?:_|\.|$)", Path(path).name)
    return m.group(1) if m else Path(path).stem


def safe_label(value):
    if value is None:
        return None
    s = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value).strip())
    s = s.strip("._-")
    if not s:
        raise ValueError("--stack-label becomes empty after sanitization.")
    return s


def choose_final_stack(args, manifest_source):
    """Return tested final-stack paths and whether they are explicit overrides."""
    explicit = {
        "icube": args.final_icube,
        "vcube": args.final_vcube,
        "mcube": args.final_mcube,
        "ecube": args.final_ecube,
    }
    supplied = [v is not None for v in explicit.values()]
    if any(supplied) and not all(supplied):
        missing = [k for k, v in explicit.items() if v is None]
        raise ValueError(
            "Explicit alternate-stack mode requires all four of "
            "--final-icube/--final-vcube/--final-mcube/--final-ecube. "
            f"Missing: {missing}"
        )

    if all(supplied):
        source = {k: Path(v).expanduser().resolve() for k, v in explicit.items()}
        mode = "explicit_override"
    else:
        source = dict(manifest_source)
        mode = "script01_manifest"

    for k, p in source.items():
        if not Path(p).is_file():
            raise FileNotFoundError(f"Tested final {k} does not exist: {p}")

    label = safe_label(args.stack_label)
    if label is None:
        if mode == "explicit_override":
            # Derive something stable from the science-stack basename.
            name = Path(source["icube"]).stem
            name = re.sub(r"_icubes?$", "", name, flags=re.IGNORECASE)
            label = safe_label(name) or "alternate_stack"
        else:
            label = "script01_control"

    return source, mode, label


def _celestial_grid_points(shape_wyx):
    _, ny, nx = shape_wyx
    return np.asarray([
        [0.0, 0.0],
        [nx - 1.0, 0.0],
        [0.0, ny - 1.0],
        [nx - 1.0, ny - 1.0],
        [(nx - 1.0) / 2.0, (ny - 1.0) / 2.0],
    ], dtype=float)


def validate_alternate_grid(reference_icube, tested_icube, tolerance_arcsec):
    """Require the alternate stack to preserve the original final output grid.

    Script-2 RH3 membership is an integer raster on the original KcwiKit grid.
    Direct alternate-stack comparisons are only exact if that grid is retained.
    """
    with fits.open(reference_icube, memmap=False) as href, fits.open(tested_icube, memmap=False) as htest:
        aref = np.asarray(href[0].data)
        atest = np.asarray(htest[0].data)
        href0 = href[0].header.copy()
        htest0 = htest[0].header.copy()

    report = {
        "reference_icube": str(reference_icube),
        "tested_icube": str(tested_icube),
        "reference_shape": list(aref.shape),
        "tested_shape": list(atest.shape),
        "shape_match": bool(aref.shape == atest.shape),
        "wcs_tolerance_arcsec": float(tolerance_arcsec),
    }
    if aref.shape != atest.shape:
        report["grid_match"] = False
        report["reason"] = "shape_mismatch"
        return report

    # Spectral axis: require effectively identical native wavelength centers.
    wref = wavelength_axis(href0, aref.shape[0])
    wtest = wavelength_axis(htest0, atest.shape[0])
    dw = np.asarray(wtest - wref, dtype=float)
    report["max_abs_spectral_grid_difference_angstrom"] = float(np.nanmax(np.abs(dw)))
    spectral_match = bool(np.allclose(wref, wtest, rtol=0.0, atol=1.0e-6))
    report["spectral_grid_match"] = spectral_match

    # Celestial output grid: compare the same output pixel coordinates.
    wr = WCS(href0).celestial
    wt = WCS(htest0).celestial
    pts = _celestial_grid_points(aref.shape)
    seps = []
    for x, y in pts:
        sr = wr.pixel_to_world(float(x), float(y))
        st = wt.pixel_to_world(float(x), float(y))
        seps.append(float(sr.separation(st).arcsec))
    report["celestial_test_point_separations_arcsec"] = seps
    report["max_celestial_grid_difference_arcsec"] = float(np.nanmax(seps))
    celestial_match = bool(
        np.all(np.isfinite(seps))
        and np.nanmax(seps) <= float(tolerance_arcsec)
    )
    report["celestial_grid_match"] = celestial_match
    report["grid_match"] = bool(spectral_match and celestial_match)
    if not spectral_match:
        report["reason"] = "spectral_grid_mismatch"
    elif not celestial_match:
        report["reason"] = "celestial_grid_mismatch"
    else:
        report["reason"] = "ok"
    return report


def first_3d_hdu(hdul):
    for h in hdul:
        if h.data is not None and np.ndim(h.data) == 3:
            return h
    raise ValueError("No 3-D image HDU found.")


def named_3d(hdul, names):
    lookup = {str(h.name).upper(): h for h in hdul}
    for name in names:
        h = lookup.get(name.upper())
        if h is not None and h.data is not None and np.ndim(h.data) == 3:
            return np.asarray(h.data)
    return None


def wavelength_axis(header, n):
    step = header.get("CDELT3", header.get("CD3_3"))
    if step is None or "CRVAL3" not in header or "CRPIX3" not in header:
        raise ValueError("Cube lacks a simple linear spectral WCS (CRVAL3/CRPIX3/CDELT3).")
    pix = np.arange(n, dtype=float) + 1.0
    wave = float(header["CRVAL3"]) + (pix - float(header["CRPIX3"])) * float(step)
    unit_text = str(header.get("CUNIT3", "Angstrom")).strip()
    try:
        unit = u.Unit(unit_text) if unit_text else u.AA
        wave = (wave * unit).to_value(u.AA)
    except Exception:
        print(f"WARNING: cannot parse CUNIT3={unit_text!r}; assuming Angstrom.")
    return np.asarray(wave, dtype=float)


class Exposure:
    def __init__(self, path):
        self.path = Path(path).expanduser().resolve()
        self.frame = frame_id(self.path)
        with fits.open(self.path, memmap=False) as h:
            sh = first_3d_hdu(h)
            self.data = np.asarray(sh.data, dtype=float)
            hdr = sh.header.copy()
            if not any(k in hdr for k in ("CTYPE1", "CTYPE2")):
                hdr.extend(h[0].header, update=True)
            self.wave = wavelength_axis(hdr, self.data.shape[0])
            self.celestial = WCS(hdr).celestial
            self.exptime = float(h[0].header.get("EXPTIME", hdr.get("EXPTIME", np.nan)))
            self.mask = named_3d(h, ("MASK",))
            self.flags = named_3d(h, ("FLAGS",))
            self.uncert = named_3d(h, ("UNCERT", "UNCERTAINTY"))
            self.inventory = [
                {"index": i, "name": str(hh.name), "shape": None if hh.data is None else list(np.shape(hh.data)),
                 "dtype": None if hh.data is None else str(np.asarray(hh.data).dtype)}
                for i, hh in enumerate(h)
            ]
        if self.data.ndim != 3:
            raise ValueError(f"Expected (wave,y,x) data in {self.path}")


def raw_primary(path):
    with fits.open(path, memmap=False) as h:
        a = np.asarray(h[0].data)
    if a.ndim != 3:
        raise ValueError(f"Expected (wave,y,x) primary cube: {path}")
    return a


def local_z(wave, flux, center, inner, outer):
    w = np.asarray(wave, float)
    f = np.asarray(flux, float)
    finite = np.isfinite(w) & np.isfinite(f)
    d = np.abs(w - float(center))
    side = finite & (d >= inner) & (d <= outer)
    if np.sum(side) < 8:
        return np.nan
    med = float(np.median(f[side]))
    sig = float(1.4826 * np.median(np.abs(f[side] - med)))
    if not np.isfinite(sig) or sig <= 0:
        return np.nan
    jgood = np.flatnonzero(finite)
    j = int(jgood[np.argmin(np.abs(w[jgood] - float(center)))])
    return float((f[j] - med) / sig)


def local_norm(wave, flux, center, inner, outer):
    w = np.asarray(wave, float)
    f = np.asarray(flux, float)
    d = np.abs(w - center)
    side = np.isfinite(f) & (d >= inner) & (d <= outer)
    out = f.copy()
    if np.sum(side) >= 8:
        med = float(np.median(f[side]))
        if np.isfinite(med) and med != 0:
            out = f / med
    return out


def bilinear_spectrum(cube, x, y):
    nw, ny, nx = cube.shape
    if not np.isfinite(x) or not np.isfinite(y):
        return np.full(nw, np.nan), []
    x0, y0 = int(math.floor(x)), int(math.floor(y))
    dx, dy = x - x0, y - y0
    terms = [
        (y0, x0, (1-dx)*(1-dy)), (y0, x0+1, dx*(1-dy)),
        (y0+1, x0, (1-dx)*dy), (y0+1, x0+1, dx*dy),
    ]
    terms = [(yy, xx, ww) for yy, xx, ww in terms if 0 <= yy < ny and 0 <= xx < nx and ww > 0]
    if not terms:
        return np.full(nw, np.nan), []
    arr = np.asarray([cube[:, yy, xx] for yy, xx, _ in terms], float)
    wt = np.asarray([ww for _, _, ww in terms], float)[:, None]
    finite = np.isfinite(arr)
    num = np.nansum(np.where(finite, arr * wt, 0.0), axis=0)
    den = np.sum(np.where(finite, wt, 0.0), axis=0)
    out = np.full(nw, np.nan)
    ok = den > 0
    out[ok] = num[ok] / den[ok]
    return out, terms


def inspect_footprint(exp, final_wcs, xout, yout, lam, spatial_radius, spectral_radius, inner, outer, extreme_z):
    sky = final_wcs.pixel_to_world(float(xout), float(yout))
    xin, yin = exp.celestial.world_to_pixel(sky)
    xin, yin = float(xin), float(yin)
    inside = (np.isfinite(xin) and np.isfinite(yin) and -0.5 <= xin <= exp.data.shape[2]-0.5 and -0.5 <= yin <= exp.data.shape[1]-0.5)
    bil, terms = bilinear_spectrum(exp.data, xin, yin)
    j = int(np.argmin(np.abs(exp.wave - lam)))
    j0, j1 = max(0, j-spectral_radius), min(exp.data.shape[0], j+spectral_radius+1)
    coords = []
    if inside:
        xc, yc = int(round(xin)), int(round(yin))
        for yy in range(yc-spatial_radius, yc+spatial_radius+1):
            for xx in range(xc-spatial_radius, xc+spatial_radius+1):
                if 0 <= yy < exp.data.shape[1] and 0 <= xx < exp.data.shape[2]:
                    coords.append((yy, xx))

    z_exact = []
    flux_exact = []
    worst_z, worst = np.nan, None
    for yy, xx in coords:
        spec = np.asarray(exp.data[:, yy, xx], float)
        z_exact.append(local_z(exp.wave, spec, lam, inner, outer))
        flux_exact.append(float(spec[j]) if np.isfinite(spec[j]) else np.nan)
        for jj in range(j0, j1):
            z = local_z(exp.wave, spec, float(exp.wave[jj]), inner, outer)
            if np.isfinite(z) and (not np.isfinite(worst_z) or z < worst_z):
                worst_z = float(z)
                worst = (yy, xx, jj, float(exp.wave[jj]), float(spec[jj]))

    z_exact = np.asarray(z_exact, float)
    flux_exact = np.asarray(flux_exact, float)
    bil_z = local_z(exp.wave, bil, lam, inner, outer)

    nearest_mask = nearest_flags = nearest_uncert = np.nan
    if inside:
        yy = int(np.clip(round(yin), 0, exp.data.shape[1]-1))
        xx = int(np.clip(round(xin), 0, exp.data.shape[2]-1))
        if exp.mask is not None and exp.mask.shape == exp.data.shape:
            nearest_mask = float(exp.mask[j, yy, xx])
        if exp.flags is not None and exp.flags.shape == exp.data.shape:
            nearest_flags = float(exp.flags[j, yy, xx])
        if exp.uncert is not None and exp.uncert.shape == exp.data.shape:
            nearest_uncert = float(exp.uncert[j, yy, xx])

    return {
        "frame": exp.frame, "path": str(exp.path), "exptime": exp.exptime,
        "mapped_x": xin, "mapped_y": yin, "inside": inside,
        "nearest_wave_index": j, "nearest_wave": float(exp.wave[j]),
        "bilinear_flux": float(bil[j]) if np.isfinite(bil[j]) else np.nan,
        "bilinear_z": float(bil_z) if np.isfinite(bil_z) else np.nan,
        "n_neighbors": len(coords),
        "median_neighbor_flux": float(np.nanmedian(flux_exact)) if np.any(np.isfinite(flux_exact)) else np.nan,
        "min_neighbor_flux": float(np.nanmin(flux_exact)) if np.any(np.isfinite(flux_exact)) else np.nan,
        "median_neighbor_z": float(np.nanmedian(z_exact)) if np.any(np.isfinite(z_exact)) else np.nan,
        "min_neighbor_z": float(np.nanmin(z_exact)) if np.any(np.isfinite(z_exact)) else np.nan,
        "nearby_worst_z": float(worst_z) if np.isfinite(worst_z) else np.nan,
        "nearby_worst_record": worst,
        "nearby_extreme": bool(np.isfinite(worst_z) and worst_z <= -abs(extreme_z)),
        "nearest_mask": nearest_mask, "nearest_flags": nearest_flags, "nearest_uncert": nearest_uncert,
        "bilinear_spectrum": bil, "bilinear_terms": terms,
    }


def possible_subsets(total, exposures, tol):
    if not np.isfinite(total):
        return []
    finite = [(e.frame, e.exptime) for e in exposures if np.isfinite(e.exptime)]
    ans = []
    for r in range(1, len(finite)+1):
        for combo in combinations(finite, r):
            if abs(sum(v for _, v in combo) - total) <= tol:
                ans.append(tuple(fid for fid, _ in combo))
    return ans


def classify(final_z, possible, inputs, final_thresh, input_thresh):
    if not np.isfinite(final_z) or final_z > -abs(final_thresh):
        return "FINAL_OUTPUT_NOT_CATASTROPHIC"
    nearby = [r["frame"] for r in inputs if r["nearby_extreme"]]
    bilbad = [r["frame"] for r in inputs if np.isfinite(r["bilinear_z"]) and r["bilinear_z"] <= -abs(input_thresh)]
    if len(possible) == 1 and len(possible[0]) == 1:
        f = possible[0][0]
        if f in nearby or f in bilbad:
            return f"SINGLE_CONTRIBUTOR_{f}_HAS_NEGATIVE_INPUT"
        return f"SINGLE_CONTRIBUTOR_{f}_LOCAL_INPUT_NORMAL"
    if len(nearby) == 1:
        return f"ONE_FRAME_NEARBY_NEGATIVE_INPUT_{nearby[0]}"
    if len(nearby) > 1:
        return "MULTIPLE_FRAMES_HAVE_NEARBY_NEGATIVE_INPUT"
    if bilbad:
        return "FRACTIONAL_INPUT_SAMPLE_ALREADY_NEGATIVE"
    return "ALL_LOCAL_INPUT_NEIGHBORHOODS_NORMAL_STACK_REPROJECTION_SUSPECT"


def plot_candidate(path, bid, rank, rest, obs, final_wave, output_rows, exposures, input_by_xy, half, inner, outer, stack_label):
    nrows = 2 + len(exposures)
    fig, axes = plt.subplots(nrows, 1, figsize=(13, 2.8*nrows), squeeze=False)
    axes = axes[:, 0]

    ax = axes[0]
    for r in output_rows:
        m = np.abs(final_wave - obs) <= half
        norm = local_norm(final_wave, r["final_spectrum"], obs, inner, outer)
        ax.plot(final_wave[m], norm[m], lw=1.0,
                label=f"out({r['x']},{r['y']}) E={r['ecube']:.0f}s z={r['final_z']:.1f}")
    ax.axvline(obs, ls="--", lw=0.8)
    ax.set_ylabel("local-normalized flux")
    ax.set_title("Final KcwiKit member spaxels")
    ax.legend(fontsize=7)

    ax = axes[1]
    labels = [f"({r['x']},{r['y']})" for r in output_rows]
    ax.bar(np.arange(len(labels)), [r["ecube"] for r in output_rows])
    ax.set_xticks(np.arange(len(labels)), labels=labels)
    ax.set_ylabel("e-cube [s]")
    ax.set_title("Effective exposure at exact candidate channel")

    for ii, exp in enumerate(exposures):
        ax = axes[2+ii]
        for r in output_rows:
            ir = next(q for q in input_by_xy[(r["y"], r["x"])] if q["frame"] == exp.frame)
            m = np.abs(exp.wave - obs) <= half
            norm = local_norm(exp.wave, ir["bilinear_spectrum"], obs, inner, outer)
            ax.plot(exp.wave[m], norm[m], lw=1.0,
                    label=(f"out({r['x']},{r['y']}) -> in({ir['mapped_x']:.2f},{ir['mapped_y']:.2f}) "
                           f"z={ir['bilinear_z']:.1f}, near-min-z={ir['nearby_worst_z']:.1f}"))
        ax.axvline(obs, ls="--", lw=0.8)
        ax.set_ylabel("local-normalized flux")
        ax.set_title(f"Input frame {exp.frame}: fractional sample + 5x5 neighborhood search")
        ax.legend(fontsize=7)
    axes[-1].set_xlabel("Observed wavelength [A]")
    fig.suptitle(
        f"03g [{stack_label}] | bin {bid} candidate {rank} | "
        f"rest={rest:.2f} A | observed={obs:.2f} A"
    )
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_thum(path, xy, outpath, title):
    with fits.open(path, memmap=False) as h:
        imgs = [np.asarray(hh.data, float) for hh in h if hh.data is not None and np.ndim(hh.data) == 2]
    if not imgs:
        return
    ncol = min(2, len(imgs)); nrow = int(math.ceil(len(imgs)/ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(6*ncol, 5*nrow), squeeze=False)
    for i, ax in enumerate(axes.ravel()):
        if i >= len(imgs):
            ax.axis("off"); continue
        im = imgs[i]
        vals = im[np.isfinite(im)]
        lo, hi = np.percentile(vals, [5, 99]) if vals.size else (0, 1)
        ax.imshow(im, origin="lower", vmin=lo, vmax=hi)
        if xy:
            ax.scatter([x for x, y in xy], [y for x, y in xy], marker="x", s=60)
        ax.set_title(f"extension {i}")
    fig.suptitle(title); fig.tight_layout(rect=(0,0,1,0.97)); fig.savefig(outpath, dpi=160); plt.close(fig)


def main():
    args = parser().parse_args()
    run3 = Path(args.script3_run).expanduser().resolve()
    cfg = crd.load_config(args.config, validate=True, strict_paths=False)
    runs_root = Path(cfg.RUNS_ROOT).expanduser().resolve()

    s3m = read_json(run3 / "metadata" / "script03_manifest.json")
    run2 = resolve_run(args.script2_run, s3m.get("source_script2_run"), runs_root,
                       Path("products/master_bin_membership.npz"), "--script2-run")
    s2m = read_json(run2 / "metadata" / "script02_manifest.json")
    run1 = resolve_run(args.script1_run, s2m.get("source_script1_run") or s3m.get("source_script1_run"),
                       runs_root, Path("metadata/script01_manifest.json"), "--script1-run")

    trace_path = Path(args.trace_table).expanduser().resolve() if args.trace_table else run3/"validation"/"03f_stage_trace"/"products"/"RH3_03f_stage_transition_trace.ecsv"
    trace = Table.read(trace_path, format="ascii.ecsv")
    if args.bins:
        wanted = {int(x.strip()) for x in args.bins.split(",") if x.strip()}
        trace = trace[[int(v) in wanted for v in trace["BIN_ID"]]]
    if len(trace) == 0: raise ValueError("No 03f candidates selected.")

    m1 = read_json(run1 / "metadata" / "script01_manifest.json")
    manifest_source = {
        k: Path(v).expanduser().resolve()
        for k, v in m1["RH3_source_paths"].items()
        if k in ("icube", "vcube", "mcube", "ecube")
    }
    if set(manifest_source) != {"icube", "vcube", "mcube", "ecube"}:
        raise KeyError("Script-1 manifest does not contain all four RH3 KcwiKit source paths.")

    source, source_mode, stack_label = choose_final_stack(args, manifest_source)

    # Each alternate stack gets its own output subtree by default. This prevents
    # bilinear / leave-one-out / single-frame experiments from overwriting one another.
    out = (
        Path(args.output_dir).expanduser().resolve()
        if args.output_dir
        else run3 / "validation" / "03g_kcwikit_forensics" / stack_label
    )
    products, figures, metadata = out/"products", out/"figures", out/"metadata"
    for d in (products, figures, metadata):
        d.mkdir(parents=True, exist_ok=True)

    grid_report = validate_alternate_grid(
        manifest_source["icube"],
        source["icube"],
        float(args.grid_match_tolerance_arcsec),
    )
    write_json(grid_report, products / "RH3_03g_stack_grid_match.json")
    if not grid_report["grid_match"]:
        raise ValueError(
            "The tested alternate KcwiKit stack does not preserve the original "
            "Script-2 output grid, so the saved RH3 membership raster cannot be "
            "applied exactly. Rebuild the alternate stack with the same final "
            "dimension/orientation/xpix/ypix/WCS geometry. Grid report: "
            f"{products / 'RH3_03g_stack_grid_match.json'}"
        )

    print("\n03g tested stack")
    print("-----------------")
    print(f"label: {stack_label}")
    print(f"source mode: {source_mode}")
    for k in ("icube", "vcube", "mcube", "ecube"):
        print(f"{k}: {source[k]}")
    print(
        "grid match to Script-1 control: PASS | "
        f"max celestial difference={grid_report['max_celestial_grid_difference_arcsec']:.6f} arcsec | "
        f"max spectral difference={grid_report['max_abs_spectral_grid_difference_angstrom']:.6g} A"
    )

    stack = io.load_kcwikit_stack(source["icube"], source["vcube"], source["mcube"], source["ecube"],
                                  arm="RH3", min_good_wavelength_fraction=float(cfg.MIN_GOOD_WAVELENGTH_FRACTION),
                                  bad_channel_fraction_threshold=float(cfg.BAD_CHANNEL_FRACTION_THRESHOLD),
                                  float_dtype=str(cfg.STACK_FLOAT_DTYPE))
    ri, rv, rm, rexp = [raw_primary(source[k]) for k in ("icube","vcube","mcube","ecube")]
    with fits.open(source["icube"], memmap=False) as h: final_wcs = WCS(h[0].header).celestial
    final_wave = np.asarray(stack.wavelength, float)

    exposures = sorted([Exposure(p) for p in args.post_zap_cube], key=lambda e: e.frame)
    write_json({e.frame:{"path":str(e.path),"exptime":e.exptime,"shape":list(e.data.shape),"extensions":e.inventory} for e in exposures},
               products/"RH3_03g_input_exposure_inventory.json")

    with np.load(run2/"products"/"master_bin_membership.npz", allow_pickle=False) as z:
        ymem = np.asarray(z["rh3_y_pix"], int); xmem = np.asarray(z["rh3_x_pix"], int); bmem = np.asarray(z["rh3_bin_id"], int)

    output_rows, input_rows, cand_rows = [], [], []
    bad_xy = set()

    for c in trace:
        bid, rank, j = int(c["BIN_ID"]), int(c["CANDIDATE_RANK"]), int(c["NATIVE_INDEX"])
        rest, obs = float(c["REST_WAVE_ANGSTROM"]), float(c["OBSERVED_WAVE_ANGSTROM"])
        member = bmem == bid
        yy, xx = ymem[member], xmem[member]
        valid = np.asarray([stack.good[int(y), int(x), j] for y,x in zip(yy,xx)], bool)
        yy, xx = yy[valid], xx[valid]
        vals = np.asarray([stack.flux[int(y), int(x), j] for y,x in zip(yy,xx)], float)
        if vals.size == 0: continue
        choose = np.argsort(vals)[:min(args.max_output_spaxels, vals.size)]
        per_candidate = []; by_xy = {}

        for ii in choose:
            y, x = int(yy[ii]), int(xx[ii]); bad_xy.add((x,y))
            spec = np.asarray(stack.flux[y,x,:], float)
            fz = local_z(final_wave, spec, obs, args.sideband_inner, args.sideband_outer)
            evalue = float(rexp[j,y,x]); poss = possible_subsets(evalue, exposures, args.exposure_time_tolerance)
            inputs = [inspect_footprint(e, final_wcs, x, y, obs, args.spatial_radius, args.spectral_radius,
                                       args.sideband_inner, args.sideband_outer, args.input_extreme_z) for e in exposures]
            by_xy[(y,x)] = inputs
            cls = classify(fz, poss, inputs, args.final_pathology_z, args.input_extreme_z)
            subset_text = ";".join("+".join(s) for s in poss)
            r = {"x":x,"y":y,"ecube":evalue,"final_z":fz,"final_spectrum":spec,"classification":cls}
            per_candidate.append(r)
            output_rows.append((bid,rank,rest,obs,j,x,y,float(ri[j,y,x]),float(rv[j,y,x]),float(rm[j,y,x]),evalue,float(fz),subset_text,len(poss),cls))
            for ir in inputs:
                input_rows.append((bid,rank,rest,obs,x,y,ir["frame"],ir["path"],ir["mapped_x"],ir["mapped_y"],ir["inside"],ir["exptime"],
                                   ir["nearest_wave_index"],ir["nearest_wave"],ir["bilinear_flux"],ir["bilinear_z"],ir["n_neighbors"],
                                   ir["median_neighbor_flux"],ir["min_neighbor_flux"],ir["median_neighbor_z"],ir["min_neighbor_z"],
                                   ir["nearby_worst_z"],ir["nearby_extreme"],ir["nearest_mask"],ir["nearest_flags"],ir["nearest_uncert"]))

        classes = [r["classification"] for r in per_candidate]
        if any("ONE_FRAME_NEARBY_NEGATIVE_INPUT" in q or "SINGLE_CONTRIBUTOR" in q and "HAS_NEGATIVE_INPUT" in q for q in classes):
            cclass = "ONE_INPUT_FRAME_IS_PLAUSIBLE_CULPRIT"
        elif any("ALL_LOCAL_INPUT_NEIGHBORHOODS_NORMAL" in q for q in classes):
            cclass = "BAD_FINAL_OUTPUT_WITH_NORMAL_LOCAL_INPUTS"
        elif any("MULTIPLE_FRAMES" in q for q in classes):
            cclass = "MULTIPLE_INPUT_FRAMES_HAVE_NEARBY_NEGATIVES"
        else:
            cclass = "MIXED_OR_UNRESOLVED"
        ev = np.asarray([r["ecube"] for r in per_candidate], float)
        cand_rows.append((bid,rank,rest,obs,j,len(vals),len(per_candidate),float(np.nanmin(vals)),float(np.nansum(vals)),
                          float(np.nanmin(ev)),float(np.nanmax(ev)),cclass))

        plot_candidate(figures/f"bin_{bid:04d}_candidate_{rank:02d}_kcwikit_forensics.png", bid, rank, rest, obs,
                       final_wave, per_candidate, exposures, by_xy, args.plot_half_width,
                       args.sideband_inner, args.sideband_outer, stack_label)

    ot = Table(rows=output_rows, names=("BIN_ID","CANDIDATE_RANK","REST_WAVE_ANGSTROM","OBSERVED_WAVE_ANGSTROM","NATIVE_INDEX","X","Y",
                                        "FINAL_FLUX","FINAL_VARIANCE","FINAL_MASK_RAW","ECUBE_SECONDS","FINAL_LOCAL_Z",
                                        "POSSIBLE_EXPOSURE_SUBSETS","N_POSSIBLE_SUBSETS","CLASSIFICATION"))
    ot.meta["DIAGNOSTIC_ONLY"] = True
    ot.write(products/"RH3_03g_final_output_voxels.ecsv", format="ascii.ecsv", overwrite=True)

    it = Table(rows=input_rows, names=("BIN_ID","CANDIDATE_RANK","REST_WAVE_ANGSTROM","OBSERVED_WAVE_ANGSTROM","FINAL_X","FINAL_Y","EXPOSURE_ID","POST_ZAP_CUBE",
                                       "MAPPED_INPUT_X","MAPPED_INPUT_Y","INSIDE_SPATIAL_FOOTPRINT","INPUT_EXPTIME_SECONDS","INPUT_NEAREST_WAVE_INDEX",
                                       "INPUT_NEAREST_WAVE_ANGSTROM","BILINEAR_FLUX_EXACT","BILINEAR_LOCAL_Z","N_SPATIAL_NEIGHBORS",
                                       "MEDIAN_NEIGHBOR_FLUX_EXACT","MIN_NEIGHBOR_FLUX_EXACT","MEDIAN_NEIGHBOR_LOCAL_Z_EXACT","MIN_NEIGHBOR_LOCAL_Z_EXACT",
                                       "NEARBY_MOST_NEGATIVE_LOCAL_Z","NEARBY_EXTREME_NEGATIVE_PRESENT","NEAREST_MASK_VALUE","NEAREST_FLAGS_VALUE","NEAREST_UNCERT_VALUE"))
    it.meta["SPATIAL_RADIUS"] = args.spatial_radius; it.meta["SPECTRAL_RADIUS_CHANNELS"] = args.spectral_radius
    it.write(products/"RH3_03g_input_footprint_trace.ecsv", format="ascii.ecsv", overwrite=True)

    ct = Table(rows=cand_rows, names=("BIN_ID","CANDIDATE_RANK","REST_WAVE_ANGSTROM","OBSERVED_WAVE_ANGSTROM","NATIVE_INDEX",
                                      "N_VALID_MEMBER_OUTPUT_SPAXELS","N_TRACED_OUTPUT_SPAXELS","MOST_NEGATIVE_MEMBER_FLUX","SUM_VALID_MEMBER_FLUX",
                                      "MIN_TRACED_ECUBE_SECONDS","MAX_TRACED_ECUBE_SECONDS","CANDIDATE_CLASSIFICATION"))
    ct.write(products/"RH3_03g_candidate_summary.ecsv", format="ascii.ecsv", overwrite=True)

    if args.thum0: plot_thum(Path(args.thum0).expanduser().resolve(), sorted(bad_xy), figures/"RH3_03g_bad_output_spaxels_on_thum0.png", "03g bad output spaxels on thum0")
    if args.thum: plot_thum(Path(args.thum).expanduser().resolve(), sorted(bad_xy), figures/"RH3_03g_bad_output_spaxels_on_thum.png", "03g bad output spaxels on aligned thum")

    write_json({
                "script":"03g_diagnose_KcwiKit_stack.py",
                "created_utc":datetime.now(timezone.utc).isoformat(),
                "script1_run":str(run1),"script2_run":str(run2),"script3_run":str(run3),
                "trace_table":str(trace_path),
                "stack_label":stack_label,
                "final_stack_source_mode":source_mode,
                "script01_control_kcwikit_sources":{k:str(v) for k,v in manifest_source.items()},
                "tested_final_kcwikit_sources":{k:str(v) for k,v in source.items()},
                "stack_grid_match_report":grid_report,
                "post_zap_exposures":[{"frame":e.frame,"path":str(e.path),"exptime":e.exptime} for e in exposures],
                "limitation":"Bilinear local sampling is diagnostic and does not exactly reproduce Montage drizzle.",
                "parameters":vars(args)
               }, metadata/"script03g_manifest.json")

    print(f"\nSCRIPT 03g CANDIDATE SUMMARY [{stack_label}]")
    print("=" * (31 + len(stack_label)))
    for r in ct:
        print(f"bin={int(r['BIN_ID'])} cand={int(r['CANDIDATE_RANK'])} rest={float(r['REST_WAVE_ANGSTROM']):.2f} A | "
              f"members={int(r['N_VALID_MEMBER_OUTPUT_SPAXELS'])} | E={float(r['MIN_TRACED_ECUBE_SECONDS']):.0f}--{float(r['MAX_TRACED_ECUBE_SECONDS']):.0f}s | "
              f"{r['CANDIDATE_CLASSIFICATION']}")
    print(f"\nProducts: {products}")
    print(f"Figures:  {figures}")
    print("No science products were modified.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
