#!/usr/bin/env python3
"""CRD_DAP Script 03c: targeted pPXF test of RH3 fixed-wavelength masks.

This is a CHEAP/SELECTIVE follow-up to Script 03d (preferred) or the legacy
Script-03b residual-derived mask experiment. It reads a completed Script-3 run
and reruns pPXF only for explicitly selected PowerBins. It never modifies the
completed Script-3 run and it does not require Scripts 1/2 to be repeated.

Typical workflow
----------------
1. Run 03b for residual/recurrence QC.
2. Run 03d with empirical pre-ZAP/sky-model data to build
   ``RH3_03d_atmospheric_mask.ecsv`` independently of pPXF residuals.
3. Run this script on pathological examples, e.g. ``--bins 231,328``.
4. Start with ``--mode selected`` (one 1C fit + the old local-best 2C state per
   bin). If the spectra/residuals improve, use ``--mode local-grid`` or
   ``--mode full-grid`` on one or two bins before pointing
   ``RH3_ATMOSPHERIC_MASK_FILE`` at the accepted table and rebuilding Script 3.

Important statistical point
---------------------------
The candidate wavelength mask is fixed *before* every pPXF state in a tested
bin.  No state-dependent sigma clipping is used.  Thus any local/full grid
comparison remains a valid comparison on one common set of spectral samples.
"""

from __future__ import annotations

import os
for _name in (
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
):
    os.environ.setdefault(_name, "1")

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astropy.table import Table

import crd_utils as crd
from crd_utils import ppxf_utils


def _parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Targeted pPXF refit test for externally anchored or diagnostic RH3 fixed-wavelength masks."
    )
    p.add_argument("--script3-run", required=True, help="Completed Script-3 run directory.")
    p.add_argument("--config", required=True, help="Current target config, e.g. config/8143-1902.py.")
    p.add_argument(
        "--mask-table",
        default=None,
        help=(
            "Fixed observed-frame mask ECSV. Preferred source is Script 03d's "
            "RH3_03d_atmospheric_mask.ecsv. If omitted, 03c tries RH3_ATMOSPHERIC_MASK_FILE, "
            "then the default 03d product, then the legacy 03b residual-derived candidate table, "
            "and finally RH3_MASK_OBSERVED_RANGES_ANGSTROM from the config."
        ),
    )
    p.add_argument(
        "--candidate-mask-table",
        default=None,
        help="Deprecated alias for --mask-table; retained for backward compatibility.",
    )
    p.add_argument("--bins", required=True, help="Comma-separated PowerBin IDs, e.g. 231,328.")
    p.add_argument(
        "--mode",
        choices=("selected", "local-grid", "full-grid"),
        default="selected",
        help=(
            "selected: refit 1C and the previously selected 2C coordinate only; "
            "local-grid: test a small neighborhood around the old 2C minimum; "
            "full-grid: recompute all 2601 states for the requested bins."
        ),
    )
    p.add_argument(
        "--velocity-radius-cells", type=int, default=1,
        help="For local-grid, +/- this many VA/VB grid cells around the old minimum (default 1).",
    )
    p.add_argument(
        "--fraction-radius-cells", type=int, default=1,
        help="For local-grid, +/- this many fA cells around the old minimum (default 1).",
    )
    p.add_argument(
        "--output-dir", default=None,
        help="Default: <script3-run>/validation/03c_candidate_mask_test",
    )
    return p


def _load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def _parse_bins(text: str, nbin: int) -> list[int]:
    out: list[int] = []
    for token in str(text).split(","):
        token = token.strip()
        if not token:
            continue
        bid = int(token)
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
    # Fixed-point iteration is robust at optical/NIR wavelengths.
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


def _saved_rest_to_observed_science(
    rest_wave: np.ndarray, redshift: float, template_medium: str, science_medium: str
) -> np.ndarray:
    rest_science = _convert_medium(rest_wave, template_medium, science_medium)
    return rest_science * (1.0 + float(redshift))


def _candidate_mask_from_table(table_path: Path, nlog: int) -> tuple[np.ndarray, list[tuple[float, float]], str]:
    tab = Table.read(table_path, format="ascii.ecsv")
    mask = np.zeros(nlog, dtype=bool)
    intervals: list[tuple[float, float]] = []
    medium = str(tab.meta.get("OBSERVED_SCIENCE_WAVELENGTH_MEDIUM", "unknown")).strip().lower()
    for row in tab:
        if "INCLUDED_IN_ATMOSPHERIC_MASK" in tab.colnames and not bool(row["INCLUDED_IN_ATMOSPHERIC_MASK"]):
            continue
        if "INCLUDED_IN_CANDIDATE" in tab.colnames and not bool(row["INCLUDED_IN_CANDIDATE"]):
            continue
        i0 = int(row["LOG_INDEX_LO"])
        i1 = int(row["LOG_INDEX_HI"])
        if i0 < 0 or i1 >= nlog or i1 < i0:
            raise ValueError(f"Invalid log-index interval in {table_path}: {i0}..{i1}")
        mask[i0:i1 + 1] = True
        if "OBSERVED_SCIENCE_LO_ANGSTROM" in tab.colnames:
            intervals.append((float(row["OBSERVED_SCIENCE_LO_ANGSTROM"]), float(row["OBSERVED_SCIENCE_HI_ANGSTROM"])))
    return mask, intervals, medium


def _candidate_mask_from_config(
    cfg, wave_rest: np.ndarray, template_medium: str, science_medium: str
) -> tuple[np.ndarray, list[tuple[float, float]], str]:
    intervals = [(float(lo), float(hi)) for lo, hi in getattr(cfg, "RH3_MASK_OBSERVED_RANGES_ANGSTROM", [])]
    if not intervals:
        raise FileNotFoundError(
            "No 03b candidate-mask table was found and RH3_MASK_OBSERVED_RANGES_ANGSTROM is empty in the config."
        )
    obs = _saved_rest_to_observed_science(wave_rest, float(cfg.REDSHIFT), template_medium, science_medium)
    mask = np.zeros(wave_rest.size, dtype=bool)
    for lo, hi in intervals:
        mask |= (obs >= lo) & (obs <= hi)
    return mask, intervals, science_medium


def _chi2(galaxy, model, noise, good) -> float:
    g = np.asarray(good, dtype=bool)
    valid = g & np.isfinite(galaxy) & np.isfinite(model) & np.isfinite(noise) & (noise > 0)
    if np.count_nonzero(valid) < 5:
        return np.inf
    r = (np.asarray(galaxy)[valid] - np.asarray(model)[valid]) / np.asarray(noise)[valid]
    return float(np.sum(r * r))


def _robust_stats(galaxy, model, noise, good) -> dict:
    g = np.asarray(good, dtype=bool)
    valid = g & np.isfinite(galaxy) & np.isfinite(model) & np.isfinite(noise) & (noise > 0)
    if np.count_nonzero(valid) < 5:
        return {"median_abs_r": np.nan, "max_abs_r": np.nan, "mad_r": np.nan, "top1_frac": np.nan}
    r = (np.asarray(galaxy)[valid] - np.asarray(model)[valid]) / np.asarray(noise)[valid]
    a = np.abs(r)
    med = float(np.median(r))
    mad = float(1.4826 * np.median(np.abs(r - med)))
    chi = r * r
    total = float(np.sum(chi))
    return {
        "median_abs_r": float(np.median(a)),
        "max_abs_r": float(np.max(a)),
        "mad_r": mad,
        "top1_frac": float(np.max(chi) / total) if total > 0 else np.nan,
    }


def _achieved_fraction(weights: np.ndarray | None, n_templates: int) -> float:
    if weights is None:
        return np.nan
    w = np.asarray(weights, dtype=float)
    if w.size != 2 * n_templates:
        return np.nan
    denom = float(np.sum(w))
    return float(np.sum(w[:n_templates]) / denom) if np.isfinite(denom) and denom > 0 else np.nan


def _grid_indices(values: np.ndarray, center: float, radius: int) -> np.ndarray:
    i = int(np.argmin(np.abs(np.asarray(values, dtype=float) - float(center))))
    lo = max(0, i - int(radius))
    hi = min(len(values), i + int(radius) + 1)
    return np.arange(lo, hi, dtype=int)


def _refit_one(
    templates, galaxy, noise, wave, wave_temp, good, old_v, old_sigma, manifest, cfg
):
    support = manifest.get("template_support", {})
    bounds = support.get("single_velocity_bounds_kms")
    if bounds is None:
        va = np.asarray(manifest.get("grid", {}).get("VA_kms", []), dtype=float)
        if va.size:
            bounds = [float(np.min(va) - cfg.RH3_SINGLE_VELOCITY_MARGIN_KMS), float(np.max(va) + cfg.RH3_SINGLE_VELOCITY_MARGIN_KMS)]
        else:
            bounds = [-600.0, 600.0]
    return ppxf_utils.fit_single_losvd(
        templates=templates,
        galaxy=galaxy,
        noise=noise,
        velscale=float(manifest["velscale_kms"]),
        lam=wave,
        lam_temp=wave_temp,
        goodpixels=np.flatnonzero(good),
        # Match production Script 3 rather than warm-starting from the old solution.
        start_velocity=float(cfg.RH3_SINGLE_VELOCITY_START_KMS),
        start_sigma=float(cfg.RH3_SIGMA_START_KMS),
        velocity_bounds=(float(bounds[0]), float(bounds[1])),
        sigma_bounds=(float(cfg.RH3_SIGMA_MIN_KMS), float(cfg.RH3_SIGMA_MAX_KMS)),
        degree=int(cfg.RH3_DEGREE),
        mdegree=int(cfg.RH3_MDEGREE),
        regul=float(cfg.RH3_REGUL),
        keep_full=True,
    )


def _refit_state(
    templates_two, component, galaxy, noise, wave, wave_temp, good,
    va, vb, fa, sigma_a, sigma_b, manifest, cfg,
):
    return ppxf_utils.fit_fixed_two_component_state(
        templates_two_component=templates_two,
        component=component,
        galaxy=galaxy,
        noise=noise,
        velscale=float(manifest["velscale_kms"]),
        lam=wave,
        lam_temp=wave_temp,
        goodpixels=np.flatnonzero(good),
        velocity_a=float(va),
        velocity_b=float(vb),
        fraction_a=float(fa),
        # Match production Script 3: every exact state starts from the same sigma initializer.
        start_sigma_a=float(cfg.RH3_SIGMA_START_KMS),
        start_sigma_b=float(cfg.RH3_SIGMA_START_KMS),
        sigma_bounds=(float(cfg.RH3_SIGMA_MIN_KMS), float(cfg.RH3_SIGMA_MAX_KMS)),
        degree=int(cfg.RH3_DEGREE),
        mdegree=int(cfg.RH3_MDEGREE),
        regul=float(cfg.RH3_REGUL),
        keep_full=True,
    )


def _plot_bin(path: Path, bid: int, wave, galaxy, noise, good_old, candidate_mask, old_two, new_one, new_two):
    retained = good_old & ~candidate_mask
    fig, axes = plt.subplots(2, 1, figsize=(14, 8), sharex=True)

    # Robust display range from retained data/models so rejected catastrophic points do not flatten the plot.
    pool = []
    for arr in (galaxy[retained], new_one[retained] if new_one is not None else [], new_two[retained] if new_two is not None else []):
        a = np.asarray(arr, dtype=float)
        a = a[np.isfinite(a)]
        if a.size:
            pool.append(a)
    if pool:
        allv = np.concatenate(pool)
        lo, hi = np.percentile(allv, [1.0, 99.0])
        pad = 0.15 * max(hi - lo, 1e-6)
        ylim = (lo - pad, hi + pad)
    else:
        ylim = None

    axes[0].plot(wave, galaxy, lw=0.7, label="saved normalized RH3")
    axes[0].plot(wave, old_two, lw=0.8, label="old local-best 2C")
    if new_one is not None:
        axes[0].plot(wave, new_one, lw=1.0, label="masked 1C refit")
    if new_two is not None:
        axes[0].plot(wave, new_two, lw=1.0, label="masked 2C refit")
    masked = good_old & candidate_mask
    if np.any(masked):
        if ylim is not None:
            ymark = np.clip(galaxy[masked], ylim[0], ylim[1])
        else:
            ymark = galaxy[masked]
        axes[0].scatter(wave[masked], ymark, marker="x", s=24, label="candidate-masked pixel")
    if ylim is not None:
        axes[0].set_ylim(*ylim)
    axes[0].set_ylabel("Normalized flux")
    axes[0].set_title(f"Bin {bid}: targeted candidate-mask pPXF test")
    axes[0].legend(fontsize=8, ncol=2)

    for model, label in ((new_one, "masked 1C"), (new_two, "masked 2C")):
        if model is None:
            continue
        valid = retained & np.isfinite(galaxy) & np.isfinite(model) & np.isfinite(noise) & (noise > 0)
        r = np.full_like(galaxy, np.nan, dtype=float)
        r[valid] = (galaxy[valid] - model[valid]) / noise[valid]
        axes[1].plot(wave[valid], r[valid], lw=0.7, label=label)
    axes[1].axhline(0.0, lw=0.8)
    axes[1].set_xlabel("Rest-frame wavelength (Angstrom)")
    axes[1].set_ylabel("(data-model)/noise")
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=170)
    plt.close(fig)


def main() -> int:
    args = _parser().parse_args()
    cfg = crd.load_config(args.config, validate=True, strict_paths=True)
    run = Path(args.script3_run).expanduser().resolve()
    if not run.is_dir():
        raise FileNotFoundError(run)

    out = Path(args.output_dir).expanduser().resolve() if args.output_dir else run / "validation" / "03c_candidate_mask_test"
    (out / "figures").mkdir(parents=True, exist_ok=True)
    (out / "products").mkdir(parents=True, exist_ok=True)
    (out / "metadata").mkdir(parents=True, exist_ok=True)

    manifest = _load_json(run / "metadata" / "script03_manifest.json")
    with np.load(run / "products" / "RH3_log_spectra_and_local_best_fits.npz", allow_pickle=False) as z:
        spectra = {k: np.asarray(z[k]) for k in z.files}
    with np.load(run / "products" / "RH3_likelihood_cubes.npz", allow_pickle=False) as z:
        cube = {k: np.asarray(z[k]) for k in z.files}
    with np.load(run / "products" / "XSL_RH3_templates.npz", allow_pickle=False) as z:
        template_product = {k: np.asarray(z[k]) for k in z.files}

    wave = np.asarray(spectra["wavelength"], dtype=float)
    galaxy_all = np.asarray(spectra["galaxy"], dtype=float)
    noise_all = np.asarray(spectra["noise"], dtype=float)
    good_all = np.asarray(spectra["good"], dtype=bool)
    old_two_all = np.asarray(spectra["local_best_two_component_model"], dtype=float)
    templates = np.asarray(template_product["templates"], dtype=float)
    wave_temp = np.asarray(template_product["wavelength"], dtype=float)
    n_templates = templates.shape[1]
    templates_two = np.column_stack([templates, templates])
    component = np.concatenate([np.zeros(n_templates, dtype=int), np.ones(n_templates, dtype=int)])

    nbin = galaxy_all.shape[0]
    bins = _parse_bins(args.bins, nbin)

    science_medium = str(manifest.get("science_wavelength_medium", cfg.SCIENCE_WAVELENGTH_MEDIUM)).strip().lower()
    template_medium = str(manifest.get("template_wavelength_medium", cfg.TEMPLATE_WAVELENGTH_MEDIUM)).strip().lower()

    explicit_mask_arg = args.mask_table if args.mask_table is not None else args.candidate_mask_table
    candidate_paths = []
    if explicit_mask_arg is not None:
        candidate_paths.append(Path(explicit_mask_arg).expanduser().resolve())
    else:
        configured_file = getattr(cfg, "RH3_ATMOSPHERIC_MASK_FILE", None)
        if configured_file not in (None, ""):
            cp = Path(configured_file).expanduser()
            if not cp.is_absolute():
                cp = Path(cfg.PROJECT_ROOT) / cp
            candidate_paths.append(cp.resolve())
        candidate_paths.extend([
            run / "validation" / "03d_atmosphere" / "products" / "RH3_03d_atmospheric_mask.ecsv",
            run / "validation" / "03b" / "products" / "RH3_03b_candidate_observed_masks.ecsv",
        ])

    table_path = next((p for p in candidate_paths if p.is_file()), None)
    if table_path is not None:
        candidate_mask, intervals, mask_medium = _candidate_mask_from_table(table_path, wave.size)
        mask_source = str(table_path)
        if mask_medium not in {"air", "vacuum", "unknown"}:
            raise ValueError(f"Unrecognized observed mask medium in mask table: {mask_medium}")
        if mask_medium in {"air", "vacuum"} and science_medium in {"air", "vacuum"} and mask_medium != science_medium:
            raise ValueError(
                f"Mask table is in observed {mask_medium} but Script-3 science medium is {science_medium}. "
                "Regenerate the atmospheric mask from this Script-3 run; do not silently mix wavelength media."
            )
    else:
        candidate_mask, intervals, mask_medium = _candidate_mask_from_config(cfg, wave, template_medium, science_medium)
        mask_source = "target config RH3_MASK_OBSERVED_RANGES_ANGSTROM"

    if not np.any(candidate_mask):
        raise RuntimeError("Candidate mask contains zero saved log-grid pixels.")

    va_grid = np.asarray(cube["VA_grid"], dtype=float)
    vb_grid = np.asarray(cube["VB_grid"], dtype=float)
    fa_grid = np.asarray(cube["fA_grid"], dtype=float)

    rows = []
    full_payload = {}
    for bid in bins:
        galaxy = galaxy_all[bid]
        noise = noise_all[bid]
        good_old = good_all[bid]
        good = good_old & ~candidate_mask
        if np.count_nonzero(good) < int(cfg.RH3_MIN_GOOD_LOG_PIXELS):
            raise RuntimeError(
                f"Bin {bid}: candidate mask leaves only {np.count_nonzero(good)} good log pixels, "
                f"below RH3_MIN_GOOD_LOG_PIXELS={cfg.RH3_MIN_GOOD_LOG_PIXELS}."
            )

        old_one_v = float(cube["one_velocity"][bid])
        old_one_s = float(cube["one_sigma"][bid])
        old_va = float(cube["local_best_VA"][bid])
        old_vb = float(cube["local_best_VB"][bid])
        old_fa = float(cube["local_best_fA"][bid])
        old_sa = float(cube["local_best_sigma_A"][bid])
        old_sb = float(cube["local_best_sigma_B"][bid])

        one = _refit_one(
            templates, galaxy, noise, wave, wave_temp, good,
            old_one_v, old_one_s, manifest, cfg,
        )
        if not one.success:
            raise RuntimeError(f"Bin {bid}: masked 1C refit failed: {one.error_message}")

        selected = _refit_state(
            templates_two, component, galaxy, noise, wave, wave_temp, good,
            old_va, old_vb, old_fa, old_sa, old_sb, manifest, cfg,
        )
        if not selected.success:
            raise RuntimeError(f"Bin {bid}: masked selected 2C refit failed: {selected.error_message}")

        best_state = (old_va, old_vb, old_fa, selected)
        grid_result = None
        if args.mode in {"local-grid", "full-grid"}:
            if args.mode == "full-grid":
                iva = np.arange(va_grid.size)
                ivb = np.arange(vb_grid.size)
                ifa = np.arange(fa_grid.size)
            else:
                iva = _grid_indices(va_grid, old_va, args.velocity_radius_cells)
                ivb = _grid_indices(vb_grid, old_vb, args.velocity_radius_cells)
                ifa = _grid_indices(fa_grid, old_fa, args.fraction_radius_cells)

            chi = np.full((iva.size, ivb.size, ifa.size), np.inf, dtype=float)
            siga = np.full_like(chi, np.nan)
            sigb = np.full_like(chi, np.nan)
            best = None
            for ai, ia in enumerate(iva):
                for bi, ib in enumerate(ivb):
                    for fi, iff in enumerate(ifa):
                        res = _refit_state(
                            templates_two, component, galaxy, noise, wave, wave_temp, good,
                            va_grid[ia], vb_grid[ib], fa_grid[iff], old_sa, old_sb, manifest, cfg,
                        )
                        if not res.success:
                            continue
                        chi[ai, bi, fi] = res.chi2_total
                        siga[ai, bi, fi] = res.sigma[0]
                        sigb[ai, bi, fi] = res.sigma[1]
                        if best is None or res.chi2_total < best[0]:
                            best = (res.chi2_total, float(va_grid[ia]), float(vb_grid[ib]), float(fa_grid[iff]), res)
            if best is None:
                raise RuntimeError(f"Bin {bid}: all {args.mode} 2C states failed.")
            best_state = (best[1], best[2], best[3], best[4])
            grid_result = {
                "VA": va_grid[iva], "VB": vb_grid[ivb], "fA": fa_grid[ifa],
                "chi2_total": chi, "sigma_A": siga, "sigma_B": sigb,
            }
            full_payload[f"bin_{bid:04d}_VA"] = grid_result["VA"]
            full_payload[f"bin_{bid:04d}_VB"] = grid_result["VB"]
            full_payload[f"bin_{bid:04d}_fA"] = grid_result["fA"]
            full_payload[f"bin_{bid:04d}_chi2_total"] = grid_result["chi2_total"]
            full_payload[f"bin_{bid:04d}_sigma_A"] = grid_result["sigma_A"]
            full_payload[f"bin_{bid:04d}_sigma_B"] = grid_result["sigma_B"]

        new_va, new_vb, new_fa, two = best_state
        old_two = old_two_all[bid]
        old_chi_retained = _chi2(galaxy, old_two, noise, good)
        old_stats = _robust_stats(galaxy, old_two, noise, good)
        new_stats = _robust_stats(galaxy, two.bestfit, noise, good)
        achieved = _achieved_fraction(two.weights, n_templates)

        rows.append((
            bid, args.mode, int(np.count_nonzero(good_old)), int(np.count_nonzero(good)), int(np.count_nonzero(good_old & candidate_mask)),
            float(cube["one_chi2_total"][bid]), float(one.chi2_total), float(one.velocity[0]), float(one.sigma[0]),
            old_va, old_vb, old_fa, old_sa, old_sb, float(cube["local_best_chi2_total"][bid]), old_chi_retained,
            float(new_va), float(new_vb), float(new_fa), achieved, float(two.sigma[0]), float(two.sigma[1]), float(two.chi2_total),
            old_stats["median_abs_r"], new_stats["median_abs_r"], old_stats["max_abs_r"], new_stats["max_abs_r"],
            old_stats["mad_r"], new_stats["mad_r"], old_stats["top1_frac"], new_stats["top1_frac"],
        ))

        _plot_bin(
            out / "figures" / f"bin_{bid:04d}_candidate_mask_refit.png",
            bid, wave, galaxy, noise, good_old, candidate_mask, old_two,
            one.bestfit, two.bestfit,
        )

        print(
            f"Bin {bid}: kept {np.count_nonzero(good)}/{np.count_nonzero(good_old)} pixels | "
            f"2C old coordinate=({old_va:.0f},{old_vb:.0f},f={old_fa:.2f}) -> "
            f"tested best=({new_va:.0f},{new_vb:.0f},f={new_fa:.2f}) | "
            f"chi2 retained old model={old_chi_retained:.3g}, masked refit={two.chi2_total:.3g}"
        )

    names = (
        "BIN_ID", "MODE", "N_GOOD_BEFORE", "N_GOOD_AFTER", "N_MASKED_CANDIDATE",
        "OLD_1C_CHI2_FULL", "NEW_1C_CHI2_MASKED", "NEW_1C_V_KMS", "NEW_1C_SIGMA_KMS",
        "OLD_VA_KMS", "OLD_VB_KMS", "OLD_FA", "OLD_SIGMA_A_KMS", "OLD_SIGMA_B_KMS", "OLD_2C_CHI2_FULL", "OLD_2C_MODEL_CHI2_ON_MASKED_PIXSET",
        "NEW_BEST_VA_KMS", "NEW_BEST_VB_KMS", "NEW_BEST_FA_GRID", "NEW_BEST_FA_ACHIEVED", "NEW_SIGMA_A_KMS", "NEW_SIGMA_B_KMS", "NEW_2C_CHI2_MASKED",
        "OLD_MODEL_MEDIAN_ABS_R_ON_MASKED_PIXSET", "NEW_MODEL_MEDIAN_ABS_R", "OLD_MODEL_MAX_ABS_R_ON_MASKED_PIXSET", "NEW_MODEL_MAX_ABS_R",
        "OLD_MODEL_ROBUST_R_MAD_ON_MASKED_PIXSET", "NEW_MODEL_ROBUST_R_MAD", "OLD_MODEL_TOP1_CHI2_FRAC_ON_MASKED_PIXSET", "NEW_MODEL_TOP1_CHI2_FRAC",
    )
    result_table = Table(rows=rows, names=names)
    result_table.meta["SOURCE_SCRIPT3_RUN"] = str(run)
    result_table.meta["CONFIG"] = str(Path(args.config).expanduser().resolve())
    result_table.meta["MASK_SOURCE"] = mask_source
    result_table.meta["OBSERVED_MASK_WAVELENGTH_MEDIUM"] = mask_medium
    result_table.meta["N_MASK_INTERVALS"] = len(intervals)
    result_table.meta["N_LOG_PIXELS_MASKED"] = int(np.count_nonzero(candidate_mask))
    result_table.meta["MODE"] = args.mode
    result_table.write(out / "products" / "RH3_03c_candidate_mask_refit_summary.ecsv", format="ascii.ecsv", overwrite=True)

    if full_payload:
        np.savez_compressed(out / "products" / "RH3_03c_test_grids.npz", **full_payload)

    manifest_out = {
        "script": "03c_test_RH3_candidate_masks",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source_script3_run": str(run),
        "config": str(Path(args.config).expanduser().resolve()),
        "mask_source": mask_source,
        "observed_mask_wavelength_medium": mask_medium,
        "mask_intervals": [[float(lo), float(hi)] for lo, hi in intervals],
        "n_log_pixels_masked": int(np.count_nonzero(candidate_mask)),
        "bins": bins,
        "mode": args.mode,
        "scientific_note": (
            "This is a targeted mask test only. Do not substitute 03c products for a complete Script-3 likelihood run. "
            "If the externally anchored 03d mask is accepted, point RH3_ATMOSPHERIC_MASK_FILE to that table (preferred) "
            "or copy the validated intervals into RH3_MASK_OBSERVED_RANGES_ANGSTROM, then rerun Script 3."
        ),
    }
    with (out / "metadata" / "script03c_manifest.json").open("w", encoding="utf-8") as f:
        json.dump(manifest_out, f, indent=2)
        f.write("\n")

    print(f"\n03c complete. Results: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
