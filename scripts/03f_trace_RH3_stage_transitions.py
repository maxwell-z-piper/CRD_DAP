#!/usr/bin/env python3
"""CRD_DAP Script 03f: locate the first pipeline stage that creates RH3 artifacts.

This script is intentionally diagnostic-only.  It does not edit masks, cubes,
or any production products.

It traces the exact 03e surviving bad-pixel candidates through:

    final KcwiKit stacked four-cube product
        -> Script 01 prepared_RH3.fits
        -> Script 02 RH3 native binned spectrum
        -> Script 03 log-rebinning
        -> Script 03 saved normalized spectrum

The goal is to identify the *first* transition at which a candidate becomes
pathological, instead of inferring the origin from the final pPXF residual.

A second, full-cube audit compares the KcwiKit final ICUBE directly against
Script 01 prepared_RH3.fits.  Script 01 is expected to preserve the science
flux apart from the explicitly documented float32 storage conversion.

Typical Mac usage
-----------------
python scripts/03f_trace_RH3_stage_transitions.py \
    --script3-run /Users/maxpiper/CRD_DAP/runs/8143-1902_S03_20260828_122251 \
    --config config/8143-1902.py

If copied Windows manifests cannot be resolved automatically, add:
    --script1-run /Users/maxpiper/CRD_DAP/runs/<S01>
    --script2-run /Users/maxpiper/CRD_DAP/runs/<S02>
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path, PureWindowsPath

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.table import Table

import crd_utils as crd
from crd_utils import binning, io, templates


SCRIPT_NAME = "03f_trace_RH3_stage_transitions"


def _parser():
    p = argparse.ArgumentParser(
        description="Trace RH3 catastrophic pixels through KcwiKit -> Script1 -> Script2 -> Script3."
    )
    p.add_argument("--script3-run", required=True)
    p.add_argument("--config", required=True)
    p.add_argument("--script2-run", default=None)
    p.add_argument("--script1-run", default=None)
    p.add_argument(
        "--candidate-table",
        default=None,
        help=(
            "03e RH3_03e_candidate_classification.ecsv. Default: "
            "<script3-run>/validation/03e_exposure_trace/products/..."
        ),
    )
    p.add_argument(
        "--bins",
        default=None,
        help="Optional comma-separated subset of candidate-table BIN_ID values.",
    )
    p.add_argument(
        "--window-half-width-rest-angstrom",
        type=float, default=25.0,
        help="Half-width of stage-ladder diagnostic windows in rest-template Angstrom.",
    )
    p.add_argument(
        "--local-sideband-inner-angstrom", type=float, default=4.0
    )
    p.add_argument(
        "--local-sideband-outer-angstrom", type=float, default=22.0
    )
    p.add_argument(
        "--pathology-z", type=float, default=20.0,
        help="Absolute local robust-z threshold used only to label a stage as pathological.",
    )
    p.add_argument(
        "--output-dir", default=None,
        help="Default: <script3-run>/validation/03f_stage_trace",
    )
    return p


def _read_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def _write_json(obj, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, sort_keys=True, default=str) + "\n")


def _basename_from_any_path(value: str) -> str:
    s = str(value)
    return PureWindowsPath(s).name or Path(s.replace("\\", "/")).name


def _resolve_run(explicit, manifest_value, runs_root: Path, required_rel: Path, label: str) -> Path:
    candidates = []
    if explicit:
        candidates.append(Path(explicit).expanduser())
    if manifest_value not in (None, ""):
        candidates.append(Path(str(manifest_value)).expanduser())
        name = _basename_from_any_path(str(manifest_value))
        if name:
            candidates.append(runs_root / name)

    for p in candidates:
        try:
            r = p.resolve()
        except Exception:
            r = p
        if r.is_dir() and (r / required_rel).is_file():
            return r

    raise FileNotFoundError(
        f"Could not resolve {label} run containing {required_rel}. "
        f"Pass the corresponding --{label.lower().replace(' ', '')}-run explicitly."
    )


def _load_script2_arrays(run2: Path):
    path = run2 / "products" / "master_bin_spectra.fits"
    with fits.open(path, memmap=False) as h:
        wave = np.asarray(h["RH3_WAVE"].data, dtype=float).ravel()
        flux = np.asarray(h["RH3_FLUX"].data, dtype=float)
        unc = np.asarray(h["RH3_UNCERT"].data, dtype=float)
        good = np.asarray(h["RH3_GOOD"].data, dtype=bool)
        try:
            ncontrib = np.asarray(h["RH3_NCONTRIB"].data)
        except Exception:
            ncontrib = None
    return path, wave, flux, unc, good, ncontrib


def _load_script3_selected(run3: Path):
    path = run3 / "products" / "RH3_log_spectra_and_local_best_fits.npz"
    with np.load(path, allow_pickle=False) as z:
        out = {k: np.asarray(z[k]) for k in z.files}
    return path, out


def _candidate_path(run3: Path, explicit: str | None) -> Path:
    if explicit:
        p = Path(explicit).expanduser().resolve()
    else:
        p = (
            run3 / "validation" / "03e_exposure_trace" / "products"
            / "RH3_03e_candidate_classification.ecsv"
        )
    if not p.is_file():
        raise FileNotFoundError(
            f"03e candidate table not found: {p}. Pass --candidate-table explicitly."
        )
    return p


def _filter_candidate_table(tab: Table, bins_text: str | None) -> Table:
    if bins_text is None:
        return tab
    wanted = {int(x.strip()) for x in str(bins_text).split(",") if x.strip()}
    keep = np.array([int(x) in wanted for x in tab["BIN_ID"]], dtype=bool)
    return tab[keep]


def _script1_source_paths(run1: Path) -> dict[str, Path]:
    manifest = _read_json(run1 / "metadata" / "script01_manifest.json")
    raw = manifest.get("RH3_source_paths", {})
    needed = ("icube", "vcube", "mcube", "ecube")
    missing = [k for k in needed if k not in raw]
    if missing:
        raise KeyError(
            "Script-1 manifest does not contain the KcwiKit RH3 source paths "
            f"needed for an identity audit: missing {missing}"
        )
    return {k: Path(str(raw[k])).expanduser().resolve() for k in needed}


def _full_primary_identity_audit(icube: Path, prepared: Path) -> dict:
    """Compare raw final KcwiKit ICUBE against Script-1 prepared PRIMARY.

    The prepared product is deliberately float32, so the strongest test is
    prepared == raw.astype(float32), not prepared == raw float64 bit-for-bit.
    """
    with fits.open(icube, memmap=True) as hin, fits.open(prepared, memmap=True) as hout:
        a = np.asarray(hin[0].data)
        b = np.asarray(hout[0].data)

        result = {
            "kcwikit_path": str(icube),
            "prepared_path": str(prepared),
            "kcwikit_shape": list(a.shape),
            "prepared_shape": list(b.shape),
            "kcwikit_dtype": str(a.dtype),
            "prepared_dtype": str(b.dtype),
        }
        if a.shape != b.shape:
            result["shape_match"] = False
            return result
        result["shape_match"] = True

        finite = np.isfinite(a) & np.isfinite(b)
        result["n_total"] = int(a.size)
        result["n_both_finite"] = int(np.sum(finite))
        result["n_finite_mismatch"] = int(np.sum(np.isfinite(a) != np.isfinite(b)))

        if np.any(finite):
            af = np.asarray(a[finite], dtype=np.float64)
            bf = np.asarray(b[finite], dtype=np.float64)
            d = bf - af
            result["max_abs_difference_vs_float64_input"] = float(np.max(np.abs(d)))
            result["median_abs_difference_vs_float64_input"] = float(np.median(np.abs(d)))
            denom = np.maximum(np.abs(af), np.finfo(float).tiny)
            rel = np.abs(d) / denom
            result["p99_relative_difference_vs_float64_input"] = float(np.percentile(rel, 99))
            result["n_sign_flips_nonzero"] = int(
                np.sum((af != 0) & (bf != 0) & (np.signbit(af) != np.signbit(bf)))
            )

            expected = np.asarray(a[finite], dtype=np.float32)
            actual = np.asarray(b[finite], dtype=np.float32)
            result["n_mismatch_vs_expected_float32_cast"] = int(np.sum(actual != expected))
            result["fraction_mismatch_vs_expected_float32_cast"] = float(np.mean(actual != expected))
            result["prepared_matches_documented_float32_cast_exactly"] = bool(np.array_equal(actual, expected))
        return result


def _make_transfer_map(membership_path: Path, spatial_shape: tuple[int, int]) -> np.ndarray:
    with np.load(membership_path, allow_pickle=False) as z:
        y = np.asarray(z["rh3_y_pix"], dtype=int)
        x = np.asarray(z["rh3_x_pix"], dtype=int)
        bid = np.asarray(z["rh3_bin_id"], dtype=int)
    m = np.full(spatial_shape, -1, dtype=int)
    inside = (y >= 0) & (x >= 0) & (y < spatial_shape[0]) & (x < spatial_shape[1])
    if not np.all(inside):
        raise ValueError("Some Script-2 RH3 membership pixels fall outside the prepared RH3 grid.")
    m[y, x] = bid
    return m


def _pixel_area_arcsec2(cube) -> float:
    sx, sy = cube.pixel_scales_arcsec()
    return float(sx * sy)


def _local_z(wave, flux, center, inner, outer):
    w = np.asarray(wave, float)
    f = np.asarray(flux, float)
    finite = np.isfinite(w) & np.isfinite(f)
    d = np.abs(w - float(center))
    side = finite & (d >= float(inner)) & (d <= float(outer))
    if np.sum(side) < 8:
        return np.nan
    vals = f[side]
    med = float(np.median(vals))
    sig = float(1.4826 * np.median(np.abs(vals - med)))
    if not np.isfinite(sig) or sig <= 0:
        return np.nan
    j = np.flatnonzero(finite)[np.argmin(np.abs(w[finite] - float(center)))]
    return float((f[j] - med) / sig)


def _normalize_local(wave, flux, center, inner, outer):
    w = np.asarray(wave, float)
    f = np.asarray(flux, float)
    out = f.copy()
    d = np.abs(w - float(center))
    side = np.isfinite(f) & (d >= float(inner)) & (d <= float(outer))
    if np.sum(side) >= 8:
        med = float(np.median(f[side]))
        if np.isfinite(med) and med != 0:
            out = f / med
    return out


def _log_edges_from_centers(log_wave):
    w = np.asarray(log_wave, float)
    lw = np.log(w)
    mid = 0.5 * (lw[:-1] + lw[1:])
    edges = np.empty(w.size + 1, dtype=float)
    edges[1:-1] = np.exp(mid)
    edges[0] = np.exp(lw[0] - 0.5 * (lw[1] - lw[0]))
    edges[-1] = np.exp(lw[-1] + 0.5 * (lw[-1] - lw[-2]))
    return edges


def _rest_template_wave(wave_obs, redshift, science_medium, template_medium):
    rest_science = np.asarray(wave_obs, float) / (1.0 + float(redshift))
    return templates.convert_wavelength_medium(rest_science, science_medium, template_medium)


def _fresh_rebin(
    wave_obs, flux, unc, good, saved_log_wave,
    redshift, science_medium, template_medium, min_valid_fraction
):
    rest = _rest_template_wave(wave_obs, redshift, science_medium, template_medium)
    edges = _log_edges_from_centers(saved_log_wave)
    return templates.rebin_spectrum_with_diagonal_noise(
        wavelength=rest,
        flux=flux,
        uncertainty=unc,
        good=good,
        out_wavelength=saved_log_wave,
        out_edges=edges,
        min_valid_fraction=float(min_valid_fraction),
    )


def _compare_arrays(a, b, good=None):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    use = np.isfinite(a) & np.isfinite(b)
    if good is not None:
        use &= np.asarray(good, bool)
    if not np.any(use):
        return {"n": 0, "max_abs": np.nan, "median_abs": np.nan, "p99_abs": np.nan}
    d = np.abs(a[use] - b[use])
    return {
        "n": int(d.size),
        "max_abs": float(np.max(d)),
        "median_abs": float(np.median(d)),
        "p99_abs": float(np.percentile(d, 99)),
    }


def _first_pathological_stage(zs: dict[str, float], threshold: float) -> str:
    order = [
        "KCWIKIT_FINAL_STACK_COADD",
        "SCRIPT1_PREPARED_COADD",
        "SCRIPT2_SAVED_NATIVE",
        "SCRIPT3_FRESH_LOG_REBIN",
        "SCRIPT3_SAVED_LOG",
    ]
    for key in order:
        z = zs.get(key, np.nan)
        if np.isfinite(z) and z <= -abs(float(threshold)):
            return key
    return "NO_STAGE_CROSSES_THRESHOLD"


def _plot_ladder(
    path: Path, *, bid: int, rank: int,
    center_rest: float, center_obs: float,
    wave_native, stack_flux, prep_flux, s2_flux,
    wave_log, fresh_log_norm, saved_log,
    half_width_rest, redshift, science_medium, template_medium,
    inner, outer
):
    rest_native = _rest_template_wave(
        wave_native, redshift, science_medium, template_medium
    )
    fig, axes = plt.subplots(2, 1, figsize=(12, 7), sharex=False)

    mnat = np.abs(rest_native - center_rest) <= half_width_rest
    ax = axes[0]
    ax.plot(rest_native[mnat],
            _normalize_local(rest_native, stack_flux, center_rest, inner, outer)[mnat],
            lw=1.1, label="final KcwiKit stack: exact Script-2-bin coadd")
    ax.plot(rest_native[mnat],
            _normalize_local(rest_native, prep_flux, center_rest, inner, outer)[mnat],
            lw=1.0, label="Script-1 prepared: exact same coadd")
    ax.plot(rest_native[mnat],
            _normalize_local(rest_native, s2_flux, center_rest, inner, outer)[mnat],
            lw=0.9, label="Script-2 saved native bin spectrum")
    ax.axvline(center_rest, ls="--", lw=0.8)
    ax.set_ylabel("local-normalized flux")
    ax.set_title("Native-wavelength stage comparison")
    ax.legend(fontsize=8)

    mlog = np.abs(wave_log - center_rest) <= half_width_rest
    ax = axes[1]
    ax.plot(wave_log[mlog],
            _normalize_local(wave_log, fresh_log_norm, center_rest, inner, outer)[mlog],
            lw=1.1, label="fresh rebin of Script-2 native spectrum")
    ax.plot(wave_log[mlog],
            _normalize_local(wave_log, saved_log, center_rest, inner, outer)[mlog],
            lw=0.9, label="Script-3 saved normalized spectrum")
    ax.axvline(center_rest, ls="--", lw=0.8)
    ax.set_xlabel("Rest-frame template wavelength (A)")
    ax.set_ylabel("local-normalized flux")
    ax.set_title("Script-3 rebin/checkpoint comparison")
    ax.legend(fontsize=8)

    fig.suptitle(
        f"Script 03f | bin {bid} candidate {rank} | "
        f"rest={center_rest:.2f} A | observed={center_obs:.2f} A"
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(path, dpi=160)
    plt.close(fig)


def main():
    args = _parser().parse_args()
    run3 = Path(args.script3_run).expanduser().resolve()
    cfg = crd.load_config(args.config, validate=True, strict_paths=False)
    runs_root = Path(cfg.RUNS_ROOT).expanduser().resolve()

    s3manifest = _read_json(run3 / "metadata" / "script03_manifest.json")
    run2 = _resolve_run(
        args.script2_run,
        s3manifest.get("source_script2_run"),
        runs_root,
        Path("products/master_bin_spectra.fits"),
        "Script2",
    )
    s2manifest = _read_json(run2 / "metadata" / "script02_manifest.json")
    run1 = _resolve_run(
        args.script1_run,
        s2manifest.get("source_script1_run") or s3manifest.get("source_script1_run"),
        runs_root,
        Path("products/prepared_RH3.fits"),
        "Script1",
    )

    out = (
        Path(args.output_dir).expanduser().resolve()
        if args.output_dir
        else run3 / "validation" / "03f_stage_trace"
    )
    products = out/"products"
    figures = out/"figures"
    metadata = out/"metadata"
    for p in (products, figures, metadata):
        p.mkdir(parents=True, exist_ok=True)

    cand_path = _candidate_path(run3, args.candidate_table)
    candidates = _filter_candidate_table(
        Table.read(cand_path, format="ascii.ecsv"), args.bins
    )
    if len(candidates) == 0:
        raise ValueError("No 03e candidates remain after --bins selection.")

    source_paths = _script1_source_paths(run1)
    prepared_path = run1/"products"/"prepared_RH3.fits"

    # ------------------------------------------------------------------
    # A. Full-cube identity audit: this is the direct Script-1 test.
    # ------------------------------------------------------------------
    identity = _full_primary_identity_audit(source_paths["icube"], prepared_path)
    _write_json(identity, products/"RH3_03f_script1_primary_identity.json")

    print("\nScript-1 full-cube flux identity audit")
    print("--------------------------------------")
    for k, v in identity.items():
        print(f"{k}: {v}")

    # Load both through the same CRD_DAP data model.
    stack = io.load_kcwikit_stack(
        source_paths["icube"],
        source_paths["vcube"],
        source_paths["mcube"],
        source_paths["ecube"],
        arm="RH3",
        min_good_wavelength_fraction=float(cfg.MIN_GOOD_WAVELENGTH_FRACTION),
        bad_channel_fraction_threshold=float(cfg.BAD_CHANNEL_FRACTION_THRESHOLD),
        float_dtype=getattr(cfg, "STACK_FLOAT_DTYPE", "float32"),
    )
    prepared = io.load_prepared_cube(prepared_path, expected_arm="RH3")

    internal_flux_cmp = _compare_arrays(stack.flux, prepared.flux)
    internal_unc_cmp = _compare_arrays(stack.uncertainty, prepared.uncertainty)
    internal_good_equal = bool(np.array_equal(stack.good, prepared.good))

    # ------------------------------------------------------------------
    # B. Reconstruct the exact Script-2 RH3 bin membership and coadd.
    # ------------------------------------------------------------------
    membership = run2/"products"/"master_bin_membership.npz"
    transfer_map = _make_transfer_map(membership, prepared.good_spaxel.shape)

    _, wave2, flux2, unc2, good2, ncontrib2 = _load_script2_arrays(run2)
    nbin = flux2.shape[0]
    pixel_area = _pixel_area_arcsec2(prepared)
    min_member_fraction = float(cfg.BIN_SPECTRUM_MIN_MEMBER_FRACTION)

    stack_coadd = binning.coadd_bin_spectra(
        stack, transfer_map, n_bins=nbin,
        pixel_area_arcsec2=pixel_area,
        min_member_fraction=min_member_fraction,
    )
    prepared_coadd = binning.coadd_bin_spectra(
        prepared, transfer_map, n_bins=nbin,
        pixel_area_arcsec2=pixel_area,
        min_member_fraction=min_member_fraction,
    )

    s2_flux_cmp_from_prepared = _compare_arrays(
        prepared_coadd.flux, flux2, good=good2
    )
    s2_unc_cmp_from_prepared = _compare_arrays(
        prepared_coadd.uncertainty, unc2, good=good2
    )

    # ------------------------------------------------------------------
    # C. Reproduce Script-3 log rebin from the saved Script-2 spectrum.
    # ------------------------------------------------------------------
    _, s3 = _load_script3_selected(run3)
    wave3 = np.asarray(s3["wavelength"], float)
    gal3 = np.asarray(s3["galaxy"], float)
    good3 = np.asarray(s3["good"], bool)
    scale3 = np.asarray(s3["normalization_scale"], float)

    science_medium = str(
        s3manifest.get("science_wavelength_medium", "vacuum")
    ).strip().lower()
    template_medium = str(
        s3manifest.get("template_wavelength_medium", cfg.TEMPLATE_WAVELENGTH_MEDIUM)
    ).strip().lower()
    redshift = float(s3manifest.get("redshift", cfg.REDSHIFT))

    # Fresh rebin all requested bins only.
    needed_bins = sorted({int(x) for x in candidates["BIN_ID"]})
    fresh = {}
    for bid in needed_bins:
        f, n, g, vf = _fresh_rebin(
            wave2, flux2[bid], unc2[bid], good2[bid], wave3,
            redshift, science_medium, template_medium,
            float(cfg.RH3_LOG_REBIN_MIN_VALID_FRACTION),
        )
        fresh[bid] = (f, n, g, vf)

    # ------------------------------------------------------------------
    # D. Candidate-by-candidate stage ladder.
    # ------------------------------------------------------------------
    rows = []
    for row in candidates:
        bid = int(row["BIN_ID"])
        rank = int(row["CANDIDATE_RANK"])
        center_rest = float(row["REST_WAVE_ANGSTROM"])
        center_obs = float(row["OBSERVED_SCIENCE_WAVE_ANGSTROM"])

        fresh_flux = np.asarray(fresh[bid][0], float)
        scale = float(scale3[bid])
        fresh_norm = fresh_flux / scale if np.isfinite(scale) and scale > 0 else fresh_flux*np.nan

        zs = {
            "KCWIKIT_FINAL_STACK_COADD": _local_z(
                _rest_template_wave(wave2, redshift, science_medium, template_medium),
                stack_coadd.flux[bid], center_rest,
                args.local_sideband_inner_angstrom, args.local_sideband_outer_angstrom
            ),
            "SCRIPT1_PREPARED_COADD": _local_z(
                _rest_template_wave(wave2, redshift, science_medium, template_medium),
                prepared_coadd.flux[bid], center_rest,
                args.local_sideband_inner_angstrom, args.local_sideband_outer_angstrom
            ),
            "SCRIPT2_SAVED_NATIVE": _local_z(
                _rest_template_wave(wave2, redshift, science_medium, template_medium),
                flux2[bid], center_rest,
                args.local_sideband_inner_angstrom, args.local_sideband_outer_angstrom
            ),
            "SCRIPT3_FRESH_LOG_REBIN": _local_z(
                wave3, fresh_norm, center_rest,
                args.local_sideband_inner_angstrom, args.local_sideband_outer_angstrom
            ),
            "SCRIPT3_SAVED_LOG": _local_z(
                wave3, gal3[bid], center_rest,
                args.local_sideband_inner_angstrom, args.local_sideband_outer_angstrom
            ),
        }

        first = _first_pathological_stage(zs, args.pathology_z)

        j2 = int(np.nanargmin(np.abs(
            _rest_template_wave(wave2, redshift, science_medium, template_medium) - center_rest
        )))
        j3 = int(np.nanargmin(np.abs(wave3 - center_rest)))

        rows.append((
            bid, rank, center_rest, center_obs,
            int(j2), float(wave2[j2]),
            int(j3), float(wave3[j3]),
            float(zs["KCWIKIT_FINAL_STACK_COADD"]),
            float(zs["SCRIPT1_PREPARED_COADD"]),
            float(zs["SCRIPT2_SAVED_NATIVE"]),
            float(zs["SCRIPT3_FRESH_LOG_REBIN"]),
            float(zs["SCRIPT3_SAVED_LOG"]),
            first,
            float(stack_coadd.flux[bid, j2]) if np.isfinite(stack_coadd.flux[bid, j2]) else np.nan,
            float(prepared_coadd.flux[bid, j2]) if np.isfinite(prepared_coadd.flux[bid, j2]) else np.nan,
            float(flux2[bid, j2]) if np.isfinite(flux2[bid, j2]) else np.nan,
            float(fresh_norm[j3]) if np.isfinite(fresh_norm[j3]) else np.nan,
            float(gal3[bid, j3]) if np.isfinite(gal3[bid, j3]) else np.nan,
            bool(good2[bid, j2]),
            bool(fresh[bid][2][j3]),
            bool(good3[bid, j3]),
        ))

        _plot_ladder(
            figures/f"bin_{bid:04d}_candidate_{rank:02d}_stage_ladder.png",
            bid=bid, rank=rank,
            center_rest=center_rest, center_obs=center_obs,
            wave_native=wave2,
            stack_flux=stack_coadd.flux[bid],
            prep_flux=prepared_coadd.flux[bid],
            s2_flux=flux2[bid],
            wave_log=wave3,
            fresh_log_norm=fresh_norm,
            saved_log=gal3[bid],
            half_width_rest=float(args.window_half_width_rest_angstrom),
            redshift=redshift,
            science_medium=science_medium,
            template_medium=template_medium,
            inner=float(args.local_sideband_inner_angstrom),
            outer=float(args.local_sideband_outer_angstrom),
        )

    names = (
        "BIN_ID","CANDIDATE_RANK","REST_WAVE_ANGSTROM","OBSERVED_WAVE_ANGSTROM",
        "NATIVE_INDEX","NATIVE_OBS_WAVE_ANGSTROM","LOG_INDEX_NEAREST","LOG_REST_WAVE_ANGSTROM",
        "Z_KCWIKIT_FINAL_STACK_COADD","Z_SCRIPT1_PREPARED_COADD","Z_SCRIPT2_SAVED_NATIVE",
        "Z_SCRIPT3_FRESH_LOG_REBIN","Z_SCRIPT3_SAVED_LOG","FIRST_PATHOLOGICAL_STAGE",
        "FLUX_KCWIKIT_FINAL_STACK_COADD","FLUX_SCRIPT1_PREPARED_COADD","FLUX_SCRIPT2_SAVED_NATIVE",
        "FLUX_SCRIPT3_FRESH_LOG_NORM","FLUX_SCRIPT3_SAVED_LOG",
        "GOOD_SCRIPT2_NATIVE","GOOD_FRESH_LOG","GOOD_SCRIPT3_SAVED_LOG"
    )
    tab = Table(rows=rows, names=names)
    tab.meta["DIAGNOSTIC_ONLY"] = True
    tab.meta["PATHOLOGY_LOCAL_Z_THRESHOLD"] = float(args.pathology_z)
    tab.meta["SOURCE_SCRIPT1_RUN"] = str(run1)
    tab.meta["SOURCE_SCRIPT2_RUN"] = str(run2)
    tab.meta["SOURCE_SCRIPT3_RUN"] = str(run3)
    tab.meta["SOURCE_KCWIKIT_ICUBE"] = str(source_paths["icube"])
    tab.write(
        products/"RH3_03f_stage_transition_trace.ecsv",
        format="ascii.ecsv", overwrite=True
    )

    summary = {
        "script": SCRIPT_NAME,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script1_run": str(run1),
        "script2_run": str(run2),
        "script3_run": str(run3),
        "candidate_table": str(cand_path),
        "kcwikit_source_paths": {k: str(v) for k, v in source_paths.items()},
        "script1_full_primary_identity": identity,
        "kcwikit_vs_prepared_internal_flux": internal_flux_cmp,
        "kcwikit_vs_prepared_internal_uncertainty": internal_unc_cmp,
        "kcwikit_vs_prepared_goodmask_exact_equal": internal_good_equal,
        "prepared_manual_coadd_vs_script2_saved_flux": s2_flux_cmp_from_prepared,
        "prepared_manual_coadd_vs_script2_saved_uncertainty": s2_unc_cmp_from_prepared,
        "pixel_area_arcsec2": pixel_area,
        "bin_spectrum_min_member_fraction": min_member_fraction,
        "science_wavelength_medium": science_medium,
        "template_wavelength_medium": template_medium,
        "redshift": redshift,
        "products": {
            "trace_table": str(products/"RH3_03f_stage_transition_trace.ecsv"),
            "script1_identity": str(products/"RH3_03f_script1_primary_identity.json"),
            "figures": str(figures),
        },
    }
    _write_json(summary, metadata/"script03f_manifest.json")

    print("\nCross-stage exactness summary")
    print("-----------------------------")
    print("KcwiKit internal flux vs prepared:", internal_flux_cmp)
    print("KcwiKit internal uncertainty vs prepared:", internal_unc_cmp)
    print("KcwiKit GOOD vs prepared GOOD exact:", internal_good_equal)
    print("Prepared exact coadd vs Script2 flux:", s2_flux_cmp_from_prepared)
    print("Prepared exact coadd vs Script2 uncertainty:", s2_unc_cmp_from_prepared)

    print("\nCandidate first-pathology stages")
    print("--------------------------------")
    for r in tab:
        print(
            f"bin={int(r['BIN_ID'])} candidate={int(r['CANDIDATE_RANK'])} "
            f"rest={float(r['REST_WAVE_ANGSTROM']):.2f} A -> "
            f"{r['FIRST_PATHOLOGICAL_STAGE']} | "
            f"z(stack,S1,S2,S3fresh,S3saved)="
            f"({r['Z_KCWIKIT_FINAL_STACK_COADD']:.1f},"
            f"{r['Z_SCRIPT1_PREPARED_COADD']:.1f},"
            f"{r['Z_SCRIPT2_SAVED_NATIVE']:.1f},"
            f"{r['Z_SCRIPT3_FRESH_LOG_REBIN']:.1f},"
            f"{r['Z_SCRIPT3_SAVED_LOG']:.1f})"
        )

    print(f"\nScript 03f complete: {out}")
    print("No production files were modified.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
