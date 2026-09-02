#!/usr/bin/env python3
"""Optional development utility: scan continuum-S/N windows in Script-2 spectra.

This script is deliberately *not* imported by the production CRD_DAP pipeline.
It exists only to diagnose integration-test data whose configured production
window is a poor match to the available grating setup (for example, an RL cube
used to exercise code intended for RH3 CaT observations).

The utility never edits the target config and never chooses a replacement
science window automatically.  It reports how candidate observed-frame windows
behave so a developer can understand whether a pathology is localized to a bad
edge/sky region or is present across the spectrum.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table
import matplotlib.pyplot as plt


def _parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--script2-run", required=True, help="Completed Script-2 run directory")
    p.add_argument("--arm", default="RH3", choices=("BL", "RH3"))
    p.add_argument("--width", type=float, default=200.0, help="Candidate window width in observed Angstrom")
    p.add_argument("--step", type=float, default=50.0, help="Spacing between candidate window centers in Angstrom")
    p.add_argument("--min-wave", type=float, default=None, help="Optional observed-frame scan lower bound")
    p.add_argument("--max-wave", type=float, default=None, help="Optional observed-frame scan upper bound")
    p.add_argument("--min-good-channels", type=int, default=20)
    p.add_argument("--extreme-abs-sn", type=float, default=1000.0)
    p.add_argument("--print-top", type=int, default=10)
    return p


def _load_spectra(run_dir: Path, arm: str):
    path = run_dir / "products" / "master_bin_spectra.fits"
    if not path.is_file():
        raise FileNotFoundError(f"Missing Script-2 spectra product: {path}")
    prefix = arm.upper()
    with fits.open(path, memmap=True) as hdul:
        wave = np.asarray(hdul[f"{prefix}_WAVE"].data, dtype=float).ravel()
        flux = np.asarray(hdul[f"{prefix}_FLUX"].data, dtype=float)
        unc = np.asarray(hdul[f"{prefix}_UNCERT"].data, dtype=float)
        good = np.asarray(hdul[f"{prefix}_GOOD"].data, dtype=bool)
    return path, wave, flux, unc, good


def _evaluate_window(wave, flux, unc, good, lo, hi, min_good_channels, extreme_abs_sn):
    w = np.isfinite(wave) & (wave >= lo) & (wave <= hi)
    n_native = int(np.sum(w))
    nbin = flux.shape[0]
    signed = np.full(nbin, np.nan, dtype=float)
    positive = np.zeros(nbin, dtype=bool)
    enough = np.zeros(nbin, dtype=bool)

    for bid in range(nbin):
        use = (
            w
            & good[bid]
            & np.isfinite(flux[bid])
            & np.isfinite(unc[bid])
            & (unc[bid] > 0)
        )
        n = int(np.sum(use))
        if n < int(min_good_channels):
            continue
        enough[bid] = True
        med_flux = float(np.nanmedian(flux[bid, use]))
        med_unc = float(np.nanmedian(unc[bid, use]))
        if not np.isfinite(med_flux) or not np.isfinite(med_unc) or med_unc <= 0:
            continue
        signed[bid] = med_flux / med_unc
        positive[bid] = med_flux > 0

    valid_positive = enough & positive & np.isfinite(signed)
    evaluable = enough & np.isfinite(signed)
    extreme = evaluable & (np.abs(signed) > float(extreme_abs_sn))

    return {
        "n_native_channels": n_native,
        "n_bins_evaluable": int(np.sum(evaluable)),
        "n_bins_positive": int(np.sum(valid_positive)),
        "positive_fraction": float(np.mean(positive[enough])) if np.any(enough) else np.nan,
        "median_positive_sn": float(np.nanmedian(signed[valid_positive])) if np.any(valid_positive) else np.nan,
        "p10_positive_sn": float(np.nanpercentile(signed[valid_positive], 10.0)) if np.any(valid_positive) else np.nan,
        "extreme_fraction": float(np.mean(extreme[evaluable])) if np.any(evaluable) else np.nan,
    }


def _plot(table: Table, path: Path, arm: str) -> None:
    center = np.asarray(table["CENTER_A"], dtype=float)
    positive = np.asarray(table["POSITIVE_FRACTION"], dtype=float)
    median_sn = np.asarray(table["MEDIAN_POSITIVE_SN"], dtype=float)
    extreme = np.asarray(table["EXTREME_FRACTION"], dtype=float)

    fig, ax1 = plt.subplots(figsize=(9.0, 5.5))
    ax1.plot(center, positive, marker="o", label="Positive-continuum bin fraction")
    ax1.plot(center, extreme, marker="o", label="Extreme-|S/N| bin fraction")
    ax1.set_xlabel("Observed window center (Angstrom)")
    ax1.set_ylabel("Bin fraction")
    ax1.set_ylim(-0.02, 1.02)

    ax2 = ax1.twinx()
    ax2.plot(center, median_sn, marker="s", label="Median positive-bin S/N")
    ax2.set_ylabel("Median robust continuum S/N")

    handles1, labels1 = ax1.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(handles1 + handles2, labels1 + labels2, loc="best")
    ax1.set_title(f"{arm} Script-2 observed-frame S/N-window scan")
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=170, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    args = _parser().parse_args()
    run_dir = Path(args.script2_run).expanduser().resolve()
    source, wave, flux, unc, good = _load_spectra(run_dir, args.arm)

    if args.width <= 0 or args.step <= 0:
        raise ValueError("--width and --step must be positive")
    if args.min_good_channels < 2:
        raise ValueError("--min-good-channels must be at least 2")

    finite_wave = wave[np.isfinite(wave)]
    lo_bound = float(np.nanmin(finite_wave)) if args.min_wave is None else float(args.min_wave)
    hi_bound = float(np.nanmax(finite_wave)) if args.max_wave is None else float(args.max_wave)
    if hi_bound <= lo_bound:
        raise ValueError("Scan upper bound must exceed lower bound")

    half = 0.5 * float(args.width)
    first_center = lo_bound + half
    last_center = hi_bound - half
    if last_center < first_center:
        raise ValueError("Requested scan range is narrower than --width")
    centers = np.arange(first_center, last_center + 0.5 * args.step, float(args.step))

    rows = []
    for center in centers:
        lo = float(center - half)
        hi = float(center + half)
        stats = _evaluate_window(
            wave,
            flux,
            unc,
            good,
            lo,
            hi,
            args.min_good_channels,
            args.extreme_abs_sn,
        )
        rows.append(
            (
                float(center),
                lo,
                hi,
                stats["n_native_channels"],
                stats["n_bins_evaluable"],
                stats["n_bins_positive"],
                stats["positive_fraction"],
                stats["median_positive_sn"],
                stats["p10_positive_sn"],
                stats["extreme_fraction"],
            )
        )

    table = Table(
        rows=rows,
        names=(
            "CENTER_A",
            "LO_A",
            "HI_A",
            "N_NATIVE_CHANNELS",
            "N_BINS_EVALUABLE",
            "N_BINS_POSITIVE",
            "POSITIVE_FRACTION",
            "MEDIAN_POSITIVE_SN",
            "P10_POSITIVE_SN",
            "EXTREME_FRACTION",
        ),
    )

    out_table = run_dir / "products" / f"{args.arm}_SN_window_scan.ecsv"
    table.write(out_table, format="ascii.ecsv", overwrite=True)
    out_plot = run_dir / "figures" / f"{args.arm}_SN_window_scan.png"
    _plot(table, out_plot, args.arm)

    # This ranking is deliberately diagnostic only.  It does not update the
    # target config or claim that the highest-ranked interval is scientifically
    # appropriate for a production RH3/CaT analysis.
    positive_fraction = np.asarray(table["POSITIVE_FRACTION"], dtype=float)
    median_sn = np.asarray(table["MEDIAN_POSITIVE_SN"], dtype=float)
    finite = np.isfinite(positive_fraction) & np.isfinite(median_sn)
    ids = np.flatnonzero(finite)
    ids = ids[np.lexsort((-median_sn[ids], -positive_fraction[ids]))]

    print(f"Loaded: {source}")
    print(f"Scanned {len(table)} observed-frame windows for {args.arm}.")
    print(f"Saved table: {out_table}")
    print(f"Saved figure: {out_plot}")
    print("\nBest-behaved windows by positive-continuum fraction, then median S/N")
    print("(development diagnostic only; the production pipeline does not auto-select a window)\n")
    print(" center_A    lo_A    hi_A   positive_frac   median_SN   extreme_frac")
    for idx in ids[: max(0, int(args.print_top))]:
        print(
            f"{table['CENTER_A'][idx]:9.1f} {table['LO_A'][idx]:7.1f} {table['HI_A'][idx]:7.1f} "
            f"{table['POSITIVE_FRACTION'][idx]:13.3f} {table['MEDIAN_POSITIVE_SN'][idx]:11.2f} "
            f"{table['EXTREME_FRACTION'][idx]:12.3f}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
