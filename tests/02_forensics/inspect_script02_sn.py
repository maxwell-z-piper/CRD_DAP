#!/usr/bin/env python3
"""Inspect suspicious achieved-S/N values from a completed CRD_DAP Script-2 run.

This utility is intentionally diagnostic-only.  It never modifies the saved
spectra, variance, bin membership, or the original Script-2 products.  It reads
``master_bin_spectra.fits``, recomputes the legacy and robust S/N diagnostics,
prints the most suspicious bins, and writes a compact ECSV report plus three
single-panel figures for the worst RH3 bin.

Typical use
-----------
python scripts/inspect_script02_sn.py \
    --config config/8143-1902.py \
    --script2-run runs/8143-1902_S02_20260821_075709

If ``--script2-run`` is omitted, the newest complete ``*_S02_*`` run for the
configured target is used.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import json

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.table import Table

import crd_utils as crd
from crd_utils import binning


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Inspect Script-2 bin S/N pathologies.")
    parser.add_argument("--config", required=True, help="Target configuration file")
    parser.add_argument("--script2-run", default=None, help="Completed Script-2 run directory")
    parser.add_argument("--top", type=int, default=20, help="Number of suspicious RH3 bins to print")
    return parser


def _find_run(cfg, explicit: str | None) -> Path:
    if explicit is not None:
        run = Path(explicit).expanduser().resolve()
    else:
        root = Path(cfg.RUNS_ROOT).expanduser().resolve()
        safe = str(cfg.TARGET_NAME).replace(" ", "_").replace("/", "-")
        candidates = [
            p for p in root.glob(f"{safe}_S02_*")
            if (p / "products" / "master_bin_spectra.fits").is_file()
        ]
        if not candidates:
            raise FileNotFoundError(f"No complete Script-2 run found for {cfg.TARGET_NAME!r} under {root}")
        candidates.sort(key=lambda p: (p.stat().st_mtime, p.name), reverse=True)
        run = candidates[0]
    spectra = run / "products" / "master_bin_spectra.fits"
    if not spectra.is_file():
        raise FileNotFoundError(f"Missing Script-2 spectra product: {spectra}")
    return run


def _load_arm(hdul, arm: str) -> tuple[np.ndarray, binning.CoaddedBinSpectra]:
    wave = np.asarray(hdul[f"{arm}_WAVE"].data, dtype=float)
    flux = np.asarray(hdul[f"{arm}_FLUX"].data, dtype=float)
    unc = np.asarray(hdul[f"{arm}_UNCERT"].data, dtype=float)
    good = np.asarray(hdul[f"{arm}_GOOD"].data, dtype=bool)
    ncontrib = np.asarray(hdul[f"{arm}_NCONTRIB"].data, dtype=np.int16)
    nbin = flux.shape[0]
    spec = binning.CoaddedBinSpectra(
        flux=flux,
        uncertainty=unc,
        good=good,
        contributing_spaxels=ncontrib,
        n_members=np.zeros(nbin, dtype=np.int32),
        spatial_scale_factor=1.0,
        spatial_scale_reason="diagnostic reload of already-coadded spectra",
    )
    return wave, spec


def _diagnostic_table(diag: binning.BinSNDiagnostics, old_sn: np.ndarray | None) -> Table:
    tab = Table()
    nbin = diag.sn.size
    tab["BIN_ID"] = np.arange(nbin, dtype=int)
    if old_sn is not None and len(old_sn) == nbin:
        tab["OLD_TABLE_SN"] = np.asarray(old_sn, dtype=float)
    tab["ROBUST_SN"] = diag.sn
    tab["SIGNED_SN"] = diag.signed_sn
    tab["LEGACY_MEDIAN_RATIO"] = diag.legacy_median_ratio
    tab["MEDIAN_FLUX"] = diag.median_flux
    tab["MEDIAN_UNCERT"] = diag.median_uncertainty
    tab["MIN_UNCERT"] = diag.min_uncertainty
    tab["P05_UNCERT"] = diag.p05_uncertainty
    tab["NEGATIVE_FLUX_FRACTION"] = diag.negative_flux_fraction
    tab["N_GOOD_CHANNELS"] = diag.n_good_channels
    tab["POSITIVE_CONTINUUM"] = diag.positive_continuum
    return tab


def _print_worst(tab: Table, top: int) -> int:
    # Prioritize bins that reproduce the original symptom: the most negative
    # legacy median(flux/uncertainty).  If those are unavailable, use signed S/N.
    key = np.asarray(tab["LEGACY_MEDIAN_RATIO"], dtype=float)
    finite = np.isfinite(key)
    if not np.any(finite):
        raise RuntimeError("No finite RH3 S/N diagnostics were found")
    order = np.flatnonzero(finite)[np.argsort(key[finite])]
    worst = int(order[0])

    print("\nMost suspicious RH3 bins (sorted by legacy median flux/uncertainty):")
    print(
        " BIN   old/table      robust       signed       legacy      med_flux      med_unc     min_unc   negfrac  nchan"
    )
    print("-" * 118)
    for bid in order[: max(int(top), 1)]:
        old = float(tab["OLD_TABLE_SN"][bid]) if "OLD_TABLE_SN" in tab.colnames else np.nan
        print(
            f"{bid:4d}  {old:11.3g}  {float(tab['ROBUST_SN'][bid]):11.3g}  "
            f"{float(tab['SIGNED_SN'][bid]):11.3g}  {float(tab['LEGACY_MEDIAN_RATIO'][bid]):11.3g}  "
            f"{float(tab['MEDIAN_FLUX'][bid]):11.3g}  {float(tab['MEDIAN_UNCERT'][bid]):11.3g}  "
            f"{float(tab['MIN_UNCERT'][bid]):9.3g}  {float(tab['NEGATIVE_FLUX_FRACTION'][bid]):7.3f}  "
            f"{int(tab['N_GOOD_CHANNELS'][bid]):5d}"
        )
    return worst


def _plot_worst(run: Path, wave: np.ndarray, spec: binning.CoaddedBinSpectra, diag, bid: int, cfg) -> list[Path]:
    outdir = run / "figures"
    outdir.mkdir(parents=True, exist_ok=True)
    lo, hi = binning.observed_range_from_rest(tuple(cfg.RH3_SN_REST_RANGE_ANGSTROM), float(cfg.REDSHIFT))
    window = (wave >= lo) & (wave <= hi)
    good = spec.good[bid] & np.isfinite(spec.flux[bid]) & np.isfinite(spec.uncertainty[bid]) & (spec.uncertainty[bid] > 0)

    paths = []

    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.plot(wave[good], spec.flux[bid, good], lw=0.8)
    ax.axhline(0.0, ls="--", lw=1.0)
    ax.axvspan(lo, hi, alpha=0.12)
    ax.set_xlabel("Observed wavelength (Å)")
    ax.set_ylabel("Coadded flux")
    ax.set_title(f"RH3 bin {bid}: flux | median window flux={diag.median_flux[bid]:.4g}")
    p = outdir / f"RH3_SN_pathology_bin{bid:04d}_flux.png"
    fig.tight_layout(); fig.savefig(p, dpi=160); plt.close(fig); paths.append(p)

    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.plot(wave[good], spec.uncertainty[bid, good], lw=0.8)
    ax.axvspan(lo, hi, alpha=0.12)
    ax.set_xlabel("Observed wavelength (Å)")
    ax.set_ylabel("Formal uncertainty")
    ax.set_title(
        f"RH3 bin {bid}: uncertainty | median={diag.median_uncertainty[bid]:.4g}, "
        f"min={diag.min_uncertainty[bid]:.4g}"
    )
    p = outdir / f"RH3_SN_pathology_bin{bid:04d}_uncertainty.png"
    fig.tight_layout(); fig.savefig(p, dpi=160); plt.close(fig); paths.append(p)

    fig, ax = plt.subplots(figsize=(9, 4.5))
    g = good & window
    ratio = np.full(wave.shape, np.nan, dtype=float)
    ratio[g] = spec.flux[bid, g] / spec.uncertainty[bid, g]
    ax.plot(wave[g], ratio[g], marker=".", ms=3, lw=0.7)
    ax.axhline(0.0, ls="--", lw=1.0)
    ax.set_xlabel("Observed wavelength (Å)")
    ax.set_ylabel("Flux / formal uncertainty")
    ax.set_title(
        f"RH3 bin {bid}: S/N samples | robust={diag.signed_sn[bid]:.4g}, "
        f"legacy median={diag.legacy_median_ratio[bid]:.4g}"
    )
    p = outdir / f"RH3_SN_pathology_bin{bid:04d}_ratio.png"
    fig.tight_layout(); fig.savefig(p, dpi=160); plt.close(fig); paths.append(p)
    return paths


def main() -> int:
    args = _parser().parse_args()
    cfg = crd.load_config(args.config, validate=True, strict_paths=False)
    run = _find_run(cfg, args.script2_run)
    product = run / "products" / "master_bin_spectra.fits"
    table_path = run / "products" / "master_bin_table.ecsv"

    print(f"Script-2 run: {run}")
    print(f"Spectra:      {product}")

    with fits.open(product, memmap=True) as hdul:
        rh_wave, rh_spec = _load_arm(hdul, "RH3")
        diag = binning.achieved_sn_diagnostics_per_bin(
            rh_spec,
            rh_wave,
            rest_range=tuple(cfg.RH3_SN_REST_RANGE_ANGSTROM),
            redshift=float(cfg.REDSHIFT),
            min_good_channels=int(getattr(cfg, "BIN_SN_MIN_GOOD_CHANNELS", 10)),
            require_positive_continuum=bool(getattr(cfg, "BIN_SN_REQUIRE_POSITIVE_CONTINUUM", True)),
        )

        old_sn = None
        if table_path.is_file():
            old = Table.read(table_path, format="ascii.ecsv")
            if "RH3_SN" in old.colnames:
                old_sn = np.asarray(old["RH3_SN"], dtype=float)

        report = _diagnostic_table(diag, old_sn)
        worst = _print_worst(report, args.top)

        report_path = run / "products" / "RH3_SN_diagnostic_report.ecsv"
        report.write(report_path, format="ascii.ecsv", overwrite=True)
        figures = _plot_worst(run, rh_wave, rh_spec, diag, worst, cfg)

    nonpositive = np.isfinite(diag.signed_sn) & (diag.signed_sn <= 0)
    extreme = (
        (np.isfinite(diag.signed_sn) & (np.abs(diag.signed_sn) > float(getattr(cfg, "BIN_SN_EXTREME_ABS_WARNING", 1000.0))))
        | (np.isfinite(diag.legacy_median_ratio) & (np.abs(diag.legacy_median_ratio) > float(getattr(cfg, "BIN_SN_EXTREME_ABS_WARNING", 1000.0))))
    )
    summary = {
        "script2_run": str(run),
        "n_bins": int(diag.sn.size),
        "n_valid_positive_sn": int(np.sum(np.isfinite(diag.sn))),
        "n_nonpositive_signed_sn": int(np.sum(nonpositive)),
        "n_extreme_diagnostic_bins": int(np.sum(extreme)),
        "worst_bin": int(worst),
        "worst_signed_sn": float(diag.signed_sn[worst]),
        "worst_legacy_median_ratio": float(diag.legacy_median_ratio[worst]),
        "worst_median_flux": float(diag.median_flux[worst]),
        "worst_median_uncertainty": float(diag.median_uncertainty[worst]),
    }
    metadata = run / "metadata"
    metadata.mkdir(parents=True, exist_ok=True)
    summary_path = metadata / "RH3_SN_diagnostic_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")

    print("\nInterpretation guide:")
    print("  * ROBUST_SN is the new production-facing ratio median(flux)/median(uncertainty).")
    print("  * SIGNED_SN preserves that value even when the median continuum is negative.")
    print("  * LEGACY_MEDIAN_RATIO reproduces the old median(flux/uncertainty) diagnostic.")
    print("  * A negative median continuum is treated as 'no positive achieved S/N', not clipped data.")
    print("  * If legacy is enormous but signed is ordinary, tiny-uncertainty samples were dominating the old estimator.")
    print("  * If both are enormous, inspect the formal uncertainty figure/window placement before trusting the variance.")
    print(f"\nSaved report:  {report_path}")
    print(f"Saved summary: {summary_path}")
    for p in figures:
        print(f"Saved figure:  {p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
