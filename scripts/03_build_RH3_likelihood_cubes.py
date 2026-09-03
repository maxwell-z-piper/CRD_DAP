#!/usr/bin/env python3
"""CRD_DAP Script 3: build independent RH3 profile-likelihood cubes.

Before evaluating the full grid, this stage calibrates the spectral noise
covariance on the exact RH3 log-wavelength experiment using high-quality,
multi-start two-component pPXF residuals.  Four increasingly flexible covariance
models (M1--M4) are iterated to convergence, tested with simultaneous PowerBin-
bootstrap whitening diagnostics, and compared through complete production grids
for a deterministic set of 12--14 representative PowerBins.  The least-complex
residual-adequate and scientifically stable covariance model is then frozen.

For every BL-defined master PowerBin the production stage evaluates the complete
explicit ``(V_A, V_B, f_A,RH3)`` grid with that frozen covariance.  The three
grid coordinates are held fixed at exact values.  pPXF profiles over
``sigma_A``, ``sigma_B``, the full duplicated XSL SSP basis, and additive-
continuum nuisance terms.  The validated pPXF-9.4.8 cached-whitener interface
factorizes each PowerBin covariance once and reuses the same inverse-Cholesky
operator for all states in that bin.

The fundamental saved quantity is TOTAL covariance-aware chi-square, not pPXF's
reduced chi2.  The resulting surface is a profile likelihood; normalized
exp(-Delta chi2/2) values used downstream are relative-likelihood weights, not
Bayesian posterior probabilities.

Restartability
--------------
One completed PowerBin is the atomic checkpoint.  ``--resume`` reuses those
checkpoint files after verifying the configuration fingerprint.  Checkpoints
are intentionally deleted only after all final consolidated products and the
Script-3 manifest have been written successfully.  A failed/interrupted run
therefore keeps everything needed to continue.

Parallelism
-----------
``--workers N`` means N Python worker processes.  BLAS/OpenMP thread pools are
capped at one thread per process before NumPy/pPXF are imported.  Thus on a
four-core machine ``--workers 3`` launches three compute workers and avoids
four-way nested BLAS threading.  CPU affinity is not pinned; final scheduling is
still controlled by the operating system.

Interactive progress
--------------------
During the expensive per-bin grid stage, an interactive terminal gets a single
carriage-return status line with a spinner, bin-level progress bar, elapsed time,
time since the last completed bin, and ETA once at least one new bin has
finished.  The heartbeat is terminal-only and is not written into the science
log file.  It therefore shows that the parent process is alive and has not yet
received a worker exception without cluttering the permanent log.
"""

from __future__ import annotations

# These limits MUST be set before NumPy/pPXF (or crd_utils, which imports NumPy)
# are imported.  One --workers process therefore maps approximately to one CPU
# core instead of each worker spawning its own BLAS/OpenMP thread pool.
import os
for _name in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
):
    os.environ.setdefault(_name, "1")

import argparse
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
import warnings

# CAPFIT occasionally evaluates actred/prered with a zero predicted reduction
# during an otherwise valid pPXF optimization.  Those two RuntimeWarnings are
# harmless terminal noise for this brute-force stage.  Suppress ONLY those
# exact CAPFIT messages; pipeline logger warnings and all other Python warnings
# remain visible.  This top-level filter is also executed in spawned workers.
warnings.filterwarnings(
    "ignore",
    message=r"(divide by zero|invalid value) encountered in scalar divide",
    category=RuntimeWarning,
    module=r"capfit\.capfit",
)
from datetime import datetime, timezone
import hashlib
import json
import logging
from pathlib import Path
import platform
import shutil
import sys
import time
from types import SimpleNamespace

import numpy as np
from astropy.io import fits
from astropy.table import Table

import crd_utils as crd
from crd_utils import covariance_calibration, io, plotting, ppxf_grid, ppxf_utils, templates
import crd_utils.covariance_plotting as covariance_plotting


CHECKPOINT_SCHEMA_VERSION = 3
# Extra log-grid pixels used only as template-support edge safety.
TEMPLATE_EDGE_SAFETY_PIXELS = 4
# Seconds between terminal-only heartbeat refreshes while waiting for a bin.
PROGRESS_HEARTBEAT_SECONDS = 5.0
PROGRESS_BAR_WIDTH = 24
_WORKER = {}


def _format_duration(seconds: float) -> str:
    """Format a positive duration compactly for the interactive status line."""
    seconds = max(0, int(round(float(seconds))))
    hours, rem = divmod(seconds, 3600)
    minutes, secs = divmod(rem, 60)
    if hours:
        return f"{hours:d}:{minutes:02d}:{secs:02d}"
    return f"{minutes:02d}:{secs:02d}"


def _progress_bar(done: int, total: int, width: int = PROGRESS_BAR_WIDTH) -> str:
    """Return a fixed-width ASCII bin-level progress bar."""
    total = max(1, int(total))
    done = min(max(0, int(done)), total)
    frac = done / total
    filled = min(width, int(frac * width))
    return "[" + "#" * filled + "-" * (width - filled) + "]"


def _write_terminal_status(text: str, previous_width: int = 0) -> int:
    """Write one carriage-return status line without touching the logfile."""
    if not sys.stderr.isatty():
        return 0
    width = max(previous_width, len(text))
    sys.stderr.write("\r" + text.ljust(width))
    sys.stderr.flush()
    return width


def _clear_terminal_status(previous_width: int) -> None:
    """Erase the dynamic status line before a normal logger message."""
    if previous_width > 0 and sys.stderr.isatty():
        sys.stderr.write("\r" + " " * previous_width + "\r")
        sys.stderr.flush()


def _fit_status_line(
    *,
    spinner: str,
    completed_total: int,
    nbin: int,
    workers: int,
    elapsed: float,
    since_last_completion: float,
    eta_seconds: float | None,
) -> str:
    """Build the terminal-only Script-3 heartbeat/progress line."""
    pct = 100.0 * completed_total / max(1, nbin)
    if eta_seconds is None or not np.isfinite(eta_seconds):
        eta_text = "ETA warming up"
        last_text = "last new bin: none yet"
    else:
        eta_text = f"ETA {_format_duration(eta_seconds)}"
        last_text = f"last new bin {_format_duration(since_last_completion)} ago"
    return (
        f"{spinner} RUNNING {_progress_bar(completed_total, nbin)} "
        f"{completed_total}/{nbin} bins ({pct:5.1f}%) | workers={workers} | "
        f"elapsed {_format_duration(elapsed)} | {last_text} | {eta_text}"
    )


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build exact RH3 (V_A,V_B,f_A) pPXF profile-likelihood cubes."
    )
    parser.add_argument(
        "--config",
        required=True,
        help="Path to a target-specific config derived from config/target_config_template.py",
    )
    parser.add_argument(
        "--script2-run",
        default=None,
        help=(
            "Explicit Script-2 run directory containing products/master_bin_spectra.fits. "
            "If omitted, the newest complete Script-2 run for the target is used."
        ),
    )
    parser.add_argument(
        "--script1-run",
        default=None,
        help=(
            "Explicit Script-1 run directory containing products/prepared_RH3.fits and "
            "the saved RH3 LSF products. If omitted, Script 3 uses the source_script1_run "
            "recorded in the Script-2 manifest. This option is useful when run directories "
            "have been renamed or moved after Script 2 was created."
        ),
    )
    parser.add_argument(
        "--run-name",
        default=None,
        help=(
            "Explicit Script-3 run-directory name. With --resume this identifies the "
            "existing run to continue; without --resume it names a new run."
        ),
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=None,
        help=(
            "Number of Python worker processes / approximate CPU cores. "
            "If omitted, config N_WORKERS is used when set, otherwise cpu_count-1."
        ),
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Resume an interrupted Script-3 run from per-bin checkpoints.",
    )
    return parser


def _safe_target(name: str) -> str:
    return str(name).replace(" ", "_").replace("/", "-")


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def _required_script2_products(run_dir: Path) -> dict[str, Path]:
    products = run_dir / "products"
    return {
        "spectra": products / "master_bin_spectra.fits",
        "bins": products / "master_bins.fits",
        "table": products / "master_bin_table.ecsv",
        "membership": products / "master_bin_membership.npz",
        "manifest": run_dir / "metadata" / "script02_manifest.json",
    }


def _is_complete_script2_run(run_dir: Path) -> bool:
    req = _required_script2_products(run_dir)
    return all(req[key].is_file() for key in ("spectra", "bins", "table", "membership", "manifest"))


def _is_usable_script1_run(run_dir: Path) -> bool:
    """Return True when the minimum Script-1 products required by Script 3 exist."""
    run_dir = Path(run_dir)
    return (run_dir / "products" / "prepared_RH3.fits").is_file()


def _resolve_script1_run(s02_manifest: dict, explicit: str | None) -> tuple[Path, Path | None]:
    """Resolve the Script-1 provenance path, allowing an explicit path override.

    Returns
    -------
    resolved_run
        Script-1 run actually used by Script 3.
    manifest_run
        Script-1 path recorded by Script 2, or None if the manifest lacks it.
    """
    manifest_value = s02_manifest.get("source_script1_run")
    manifest_run = (
        Path(manifest_value).expanduser().resolve()
        if manifest_value not in (None, "")
        else None
    )

    if explicit is not None:
        run = Path(explicit).expanduser().resolve()
        if not run.is_dir():
            raise FileNotFoundError(f"Explicit Script-1 run does not exist: {run}")
        if not _is_usable_script1_run(run):
            raise FileNotFoundError(
                "Explicit Script-1 run is missing the prepared RH3 product required by "
                f"Script 3: {run / 'products' / 'prepared_RH3.fits'}"
            )
        return run, manifest_run

    if manifest_run is None:
        raise FileNotFoundError(
            "The Script-2 manifest does not contain source_script1_run. "
            "Pass --script1-run explicitly."
        )
    if not _is_usable_script1_run(manifest_run):
        raise FileNotFoundError(
            "The Script-1 run recorded in the Script-2 manifest is unavailable or has "
            "been moved/renamed. Recorded path: "
            f"{manifest_run}. Pass --script1-run <current-script1-run> explicitly."
        )
    return manifest_run, manifest_run


def _find_script2_run(cfg, explicit: str | None) -> Path:
    if explicit is not None:
        run = Path(explicit).expanduser().resolve()
        if not _is_complete_script2_run(run):
            missing = [str(p) for p in _required_script2_products(run).values() if not p.is_file()]
            raise FileNotFoundError(
                "Explicit Script-2 run is incomplete. Missing:\n  - " + "\n  - ".join(missing)
            )
        return run
    root = Path(cfg.RUNS_ROOT).expanduser().resolve()
    safe = _safe_target(cfg.TARGET_NAME)
    candidates = [
        p for p in root.iterdir()
        if p.is_dir() and p.name.startswith(safe) and "_S02_" in p.name and _is_complete_script2_run(p)
    ] if root.exists() else []
    if not candidates:
        raise FileNotFoundError(
            f"No complete Script-2 run found for target {cfg.TARGET_NAME!r} under {root}. "
            "Pass --script2-run explicitly."
        )
    candidates.sort(key=lambda p: (p.stat().st_mtime, p.name), reverse=True)
    return candidates[0]


def _new_run(cfg, run_name: str | None):
    if run_name is None:
        run_name = f"{_safe_target(cfg.TARGET_NAME)}_S03_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    return crd.create_run_context(cfg, run_name=run_name)


def _existing_run_context(run_dir: Path):
    """Minimal context object matching the public path attributes used by drivers."""
    run_dir = Path(run_dir).expanduser().resolve()
    ctx = SimpleNamespace(
        run_dir=run_dir,
        products_dir=run_dir / "products",
        figures_dir=run_dir / "figures",
        metadata_dir=run_dir / "metadata",
        logs_dir=run_dir / "logs",
    )
    for path in (ctx.products_dir, ctx.figures_dir, ctx.metadata_dir, ctx.logs_dir):
        path.mkdir(parents=True, exist_ok=True)
    return ctx


def _find_resume_run(cfg, run_name: str | None) -> Path:
    root = Path(cfg.RUNS_ROOT).expanduser().resolve()
    if run_name is not None:
        run = root / run_name
        if not run.is_dir():
            raise FileNotFoundError(f"Requested resume run does not exist: {run}")
        return run
    safe = _safe_target(cfg.TARGET_NAME)
    candidates = []
    if root.exists():
        for p in root.iterdir():
            if not (p.is_dir() and p.name.startswith(safe) and "_S03_" in p.name):
                continue
            if (p / "checkpoints").is_dir() and (p / "metadata" / "script03_resume_state.json").is_file():
                candidates.append(p)
    if not candidates:
        raise FileNotFoundError(
            f"No resumable Script-3 run found for target {cfg.TARGET_NAME!r} under {root}."
        )
    candidates.sort(key=lambda p: (p.stat().st_mtime, p.name), reverse=True)
    return candidates[0]


def _setup_logger(run) -> logging.Logger:
    """Append Script-3 messages to terminal, master pipeline log, and step log."""
    logger = logging.getLogger("CRD_DAP.03_build_RH3_likelihood_cubes")
    logger.setLevel(logging.INFO)
    logger.propagate = False
    for handler in list(logger.handlers):
        logger.removeHandler(handler)
        try:
            handler.close()
        except Exception:
            pass
    fmt = logging.Formatter(
        "%(asctime)s | %(levelname)-8s | CRD_DAP.03_build_RH3_likelihood_cubes | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    stream = logging.StreamHandler(sys.stdout)
    stream.setFormatter(fmt)
    logger.addHandler(stream)
    for path in (run.run_dir / "pipeline.log", run.logs_dir / "03_build_RH3_likelihood_cubes.log"):
        fh = logging.FileHandler(path, mode="a")
        fh.setFormatter(fmt)
        logger.addHandler(fh)
    return logger


def _section(logger, title: str, char: str = "-") -> None:
    line = char * max(18, len(title))
    logger.info(line)
    logger.info(title)
    logger.info(line)


def _read_script2_manifest(run_dir: Path) -> dict:
    return dict(io.read_json(run_dir / "metadata" / "script02_manifest.json"))


def _load_script2_arrays(products: dict[str, Path]):
    with fits.open(products["spectra"], memmap=False) as hdul:
        wave = np.asarray(hdul["RH3_WAVE"].data, dtype=float)
        flux = np.asarray(hdul["RH3_FLUX"].data, dtype=float)
        unc = np.asarray(hdul["RH3_UNCERT"].data, dtype=float)
        good = np.asarray(hdul["RH3_GOOD"].data, dtype=bool)
    with fits.open(products["bins"], memmap=False) as hdul:
        bin_map = np.asarray(hdul[0].data, dtype=int)
    table = Table.read(products["table"], format="ascii.ecsv")
    if not (flux.shape == unc.shape == good.shape):
        raise ValueError("RH3 FLUX/UNCERT/GOOD arrays in master_bin_spectra.fits disagree in shape.")
    if flux.ndim != 2 or flux.shape[1] != wave.size:
        raise ValueError("Unexpected RH3 master-bin spectra shape.")
    if len(table) != flux.shape[0]:
        raise ValueError("master_bin_table row count does not match RH3 spectra bin count.")
    return wave, flux, unc, good, bin_map, table


def _load_prepared_rh3_header(script1_run: Path):
    path = script1_run / "products" / "prepared_RH3.fits"
    if not path.is_file():
        raise FileNotFoundError(
            f"Script-1 prepared RH3 product is required to recover wavelength convention: {path}"
        )
    return fits.getheader(path, 0)


def _resolve_science_medium(cfg, prepared_header) -> str:
    configured = str(getattr(cfg, "SCIENCE_WAVELENGTH_MEDIUM", "auto")).strip().lower()
    inferred = templates.infer_science_wavelength_medium(prepared_header)
    if configured == "auto":
        return inferred
    if configured not in {"air", "vacuum"}:
        raise ValueError("SCIENCE_WAVELENGTH_MEDIUM must be auto, air, or vacuum.")
    if configured != inferred:
        raise ValueError(
            f"Configured science wavelength medium {configured!r} disagrees with prepared RH3 metadata {inferred!r}."
        )
    return configured


def _merge_mask_intervals(intervals: list[tuple[float, float]]) -> list[tuple[float, float]]:
    clean = []
    for lo, hi in intervals:
        lo = float(lo)
        hi = float(hi)
        if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
            raise ValueError(f"Invalid observed-frame mask interval ({lo}, {hi}).")
        clean.append((lo, hi))
    if not clean:
        return []
    clean.sort(key=lambda x: x[0])
    merged = [list(clean[0])]
    for lo, hi in clean[1:]:
        if lo <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], hi)
        else:
            merged.append([lo, hi])
    return [(float(lo), float(hi)) for lo, hi in merged]


def _resolve_observed_masks(cfg, science_medium: str) -> tuple[list[tuple[float, float]], dict]:
    """Resolve only *additional/manual* observed-frame masks for Script 3.

    Production atmospheric masking now happens upstream:

        CRD_DRP atmospheric mask
            -> Script 1 GOODWAVE/GOODMASK
            -> Script 2 RH3_GOOD
            -> Script 3 good_native

    ``RH3_ATMOSPHERIC_MASK_FILE`` is therefore a retired 03d-era pathway.  A
    non-None value is treated as a configuration error so the same atmospheric
    wavelengths cannot be applied through two independent mechanisms.
    """
    legacy = getattr(cfg, "RH3_ATMOSPHERIC_MASK_FILE", None)
    if legacy not in (None, ""):
        raise ValueError(
            "RH3_ATMOSPHERIC_MASK_FILE is a retired CRD_DAP 03d-era input. "
            "The production atmospheric mask must come from the validated "
            "CRD_DRP reduction manifest and is already inherited through "
            "Script-1/Script-2 GOOD masks. Set RH3_ATMOSPHERIC_MASK_FILE=None. "
            "Use RH3_MASK_OBSERVED_RANGES_ANGSTROM only for intentional "
            "additional/manual science exclusions."
        )

    manual = [
        (float(lo), float(hi))
        for lo, hi in getattr(cfg, "RH3_MASK_OBSERVED_RANGES_ANGSTROM", [])
    ]
    combined = _merge_mask_intervals(manual)
    provenance = {
        "manual_config_ranges": [[float(lo), float(hi)] for lo, hi in manual],
        "legacy_atmospheric_mask_file_disabled": True,
        "effective_additional_ranges": [
            [float(lo), float(hi)] for lo, hi in combined
        ],
        "science_wavelength_medium": str(science_medium),
    }
    return combined, provenance

def _apply_native_masks(
    native_wave_obs: np.ndarray,
    native_wave_rest_template: np.ndarray,
    native_good: np.ndarray,
    observed_masks: list[tuple[float, float]],
    rest_masks: list[tuple[float, float]],
    lsf_model: templates.SavedLSFModel,
) -> np.ndarray:
    good = np.asarray(native_good, dtype=bool).copy()
    good &= np.isfinite(native_wave_obs)
    good &= (
        (native_wave_obs >= float(lsf_model.empirical_min))
        & (native_wave_obs <= float(lsf_model.empirical_max))
    )
    for lo, hi in observed_masks:
        good &= ~((native_wave_obs >= float(lo)) & (native_wave_obs <= float(hi)))
    for lo, hi in rest_masks:
        good &= ~(
            (native_wave_rest_template >= float(lo))
            & (native_wave_rest_template <= float(hi))
        )
    return good


def _determine_velscale(native_rest_template: np.ndarray, cfg) -> float:
    explicit = getattr(cfg, "RH3_VELSCALE_KMS", None)
    if explicit is not None:
        value = float(explicit)
        if value <= 0:
            raise ValueError("RH3_VELSCALE_KMS must be positive or None.")
        return value
    dlog = np.diff(np.log(np.asarray(native_rest_template, dtype=float)))
    dlog = dlog[np.isfinite(dlog) & (dlog > 0)]
    if dlog.size < 5:
        raise ValueError("Cannot determine RH3 velocity scale from wavelength array.")
    return float(templates.C_KMS * np.nanmedian(dlog))


def _prepare_log_binned_spectra(
    *,
    wave_obs: np.ndarray,
    flux: np.ndarray,
    uncertainty: np.ndarray,
    good_native: np.ndarray,
    redshift: float,
    science_medium: str,
    template_medium: str,
    fit_range: tuple[float, float],
    velscale: float,
    cfg,
    observed_masks: list[tuple[float, float]],
    rest_masks: list[tuple[float, float]],
    lsf_model: templates.SavedLSFModel,
):
    rest_science = np.asarray(wave_obs, dtype=float) / (1.0 + float(redshift))
    rest_template = templates.convert_wavelength_medium(
        rest_science, science_medium, template_medium
    )
    log_wave, log_edges = templates.make_log_wavelength_grid(
        float(fit_range[0]), float(fit_range[1]), velscale=float(velscale)
    )
    nbin = flux.shape[0]
    nlog = log_wave.size
    gal = np.full((nbin, nlog), np.nan, dtype=np.float64)
    noise = np.full((nbin, nlog), np.nan, dtype=np.float64)
    good = np.zeros((nbin, nlog), dtype=bool)
    valid_fraction = np.zeros((nbin, nlog), dtype=np.float32)

    for bid in range(nbin):
        native = _apply_native_masks(
            wave_obs,
            rest_template,
            good_native[bid],
            observed_masks,
            rest_masks,
            lsf_model,
        )
        f, n, g, frac = templates.rebin_spectrum_with_diagonal_noise(
            wavelength=rest_template,
            flux=flux[bid],
            uncertainty=uncertainty[bid],
            good=native,
            out_wavelength=log_wave,
            out_edges=log_edges,
            min_valid_fraction=float(cfg.RH3_LOG_REBIN_MIN_VALID_FRACTION),
        )
        gal[bid] = f
        noise[bid] = n
        good[bid] = g
        valid_fraction[bid] = frac
    return log_wave, log_edges, gal, noise, good, valid_fraction


def _normalize_galaxy_for_ppxf(galaxy, noise, good):
    gp = np.flatnonzero(good & np.isfinite(galaxy) & np.isfinite(noise) & (noise > 0))
    if gp.size == 0:
        return galaxy, noise, np.nan
    scale = float(np.nanmedian(np.abs(np.asarray(galaxy)[gp])))
    if not np.isfinite(scale) or scale <= 0:
        scale = float(np.nanpercentile(np.abs(np.asarray(galaxy)[gp]), 75.0))
    if not np.isfinite(scale) or scale <= 0:
        return galaxy, noise, np.nan
    return np.asarray(galaxy, dtype=float) / scale, np.asarray(noise, dtype=float) / scale, scale


def _script03_velocity_support(va, vb, cfg, velscale):
    """Compute template support from the actual one- and two-component domains."""
    va = np.asarray(va, dtype=float)
    vb = np.asarray(vb, dtype=float)
    two_component_max_abs = float(max(np.max(np.abs(va)), np.max(np.abs(vb))))
    single_bounds = (
        float(np.min(va) - cfg.RH3_SINGLE_VELOCITY_MARGIN_KMS),
        float(np.max(va) + cfg.RH3_SINGLE_VELOCITY_MARGIN_KMS),
    )
    single_component_max_abs = float(max(abs(single_bounds[0]), abs(single_bounds[1])))
    max_kinematic_abs = max(two_component_max_abs, single_component_max_abs)
    dispersion_margin = float(cfg.RH3_TEMPLATE_PADDING_SIGMA) * float(cfg.RH3_SIGMA_MAX_KMS)
    edge_safety = float(TEMPLATE_EDGE_SAFETY_PIXELS) * float(velscale)
    total_padding = max_kinematic_abs + dispersion_margin + edge_safety
    return {
        "two_component_max_abs_kms": two_component_max_abs,
        "single_velocity_bounds_kms": single_bounds,
        "single_component_max_abs_kms": single_component_max_abs,
        "dispersion_margin_kms": dispersion_margin,
        "edge_safety_kms": edge_safety,
        "total_padding_kms": float(total_padding),
    }


def _coverage_preflight_bin_ids(good_all):
    """Select bins spanning the earliest, latest, and widest good-pixel support."""
    good = np.asarray(good_all, dtype=bool)
    first = np.full(good.shape[0], good.shape[1], dtype=int)
    last = np.full(good.shape[0], -1, dtype=int)
    span = np.full(good.shape[0], -1, dtype=int)
    for bid in range(good.shape[0]):
        gp = np.flatnonzero(good[bid])
        if gp.size:
            first[bid] = int(gp[0])
            last[bid] = int(gp[-1])
            span[bid] = int(gp[-1] - gp[0])
    valid = np.flatnonzero(last >= 0)
    if valid.size == 0:
        raise RuntimeError("No Script-3 bin has valid log-wavelength pixels for pPXF preflight.")
    return sorted({
        int(valid[np.argmin(first[valid])]),
        int(valid[np.argmax(last[valid])]),
        int(valid[np.argmax(span[valid])]),
    })


def _run_ppxf_template_coverage_preflight(
    *, logger, prepared, templates_two, component, gal_all, noise_all, good_all,
    log_wave, va, vb, fa, support, cfg, velscale
):
    """Run real pPXF coverage checks before any worker processes are launched."""
    bids = _coverage_preflight_bin_ids(good_all)
    fraction_probe = float(fa[len(fa) // 2])
    corner_states = [
        (float(va[0]), float(vb[0])),
        (float(va[0]), float(vb[-1])),
        (float(va[-1]), float(vb[0])),
        (float(va[-1]), float(vb[-1])),
    ]
    single_bounds = tuple(support["single_velocity_bounds_kms"])

    logger.info(
        "Template-support budget: max|V| two-component=%.1f km/s; one-component bounds=%.1f--%.1f km/s; "
        "dispersion margin=%.1f km/s; %d-pixel edge safety=%.1f km/s; total padding=%.1f km/s",
        support["two_component_max_abs_kms"], single_bounds[0], single_bounds[1],
        support["dispersion_margin_kms"], TEMPLATE_EDGE_SAFETY_PIXELS,
        support["edge_safety_kms"], support["total_padding_kms"],
    )
    logger.info(
        "pPXF wavelength-coverage preflight: galaxy log grid=%.2f--%.2f A; prepared templates=%.2f--%.2f A; test bins=%s",
        float(log_wave[0]), float(log_wave[-1]), float(prepared.wavelength[0]),
        float(prepared.wavelength[-1]), bids,
    )

    for bid in bids:
        gp = np.flatnonzero(np.asarray(good_all[bid], dtype=bool))
        galaxy_n, noise_n, scale = _normalize_galaxy_for_ppxf(gal_all[bid], noise_all[bid], good_all[bid])
        if not np.isfinite(scale):
            raise RuntimeError(f"pPXF preflight bin {bid} could not be normalized.")
        logger.info(
            "pPXF preflight bin %d good-pixel wavelength range: %.2f--%.2f A (%d pixels)",
            bid, float(log_wave[gp[0]]), float(log_wave[gp[-1]]), gp.size,
        )

        single = ppxf_utils.fit_single_losvd(
            templates=prepared.templates, galaxy=galaxy_n, noise=noise_n,
            velscale=velscale, lam=log_wave, lam_temp=prepared.wavelength,
            goodpixels=gp, start_velocity=float(cfg.RH3_SINGLE_VELOCITY_START_KMS),
            start_sigma=float(cfg.RH3_SIGMA_START_KMS), velocity_bounds=single_bounds,
            sigma_bounds=(float(cfg.RH3_SIGMA_MIN_KMS), float(cfg.RH3_SIGMA_MAX_KMS)),
            degree=int(cfg.RH3_DEGREE), mdegree=int(cfg.RH3_MDEGREE),
            regul=float(cfg.RH3_REGUL), keep_full=False,
        )
        if not single.success:
            raise RuntimeError(
                "SCRIPT03_PPXF_PREFLIGHT_FAILED | One-component control failed before worker launch. "
                f"Bin={bid}; galaxy_good_range={log_wave[gp[0]]:.2f}--{log_wave[gp[-1]]:.2f} A; "
                f"template_range={prepared.wavelength[0]:.2f}--{prepared.wavelength[-1]:.2f} A; "
                f"velocity_bounds={single_bounds}; error={single.error_message}"
            )

        for va0, vb0 in corner_states:
            state = ppxf_utils.fit_fixed_two_component_state(
                templates_two_component=templates_two, component=component, galaxy=galaxy_n,
                noise=noise_n, velscale=velscale, lam=log_wave, lam_temp=prepared.wavelength,
                goodpixels=gp, velocity_a=va0, velocity_b=vb0, fraction_a=fraction_probe,
                start_sigma_a=float(cfg.RH3_SIGMA_START_KMS),
                start_sigma_b=float(cfg.RH3_SIGMA_START_KMS),
                sigma_bounds=(float(cfg.RH3_SIGMA_MIN_KMS), float(cfg.RH3_SIGMA_MAX_KMS)),
                degree=int(cfg.RH3_DEGREE), mdegree=int(cfg.RH3_MDEGREE),
                regul=float(cfg.RH3_REGUL), keep_full=False,
            )
            if not state.success:
                raise RuntimeError(
                    "SCRIPT03_PPXF_PREFLIGHT_FAILED | Exact two-component corner state failed before worker launch. "
                    f"Bin={bid}; VA={va0:.1f}; VB={vb0:.1f}; fA={fraction_probe:.3f}; "
                    f"galaxy_good_range={log_wave[gp[0]]:.2f}--{log_wave[gp[-1]]:.2f} A; "
                    f"template_range={prepared.wavelength[0]:.2f}--{prepared.wavelength[-1]:.2f} A; "
                    f"error={state.error_message}"
                )

    logger.info(
        "pPXF template-coverage preflight: PASS (%d representative bin(s); one-component control + four two-component grid corners each)",
        len(bids),
    )
    return bids


def _worker_init(constants: dict) -> None:
    global _WORKER
    _WORKER = constants


def _fit_one_bin_worker(payload):
    """Fit one PowerBin with the frozen, selected covariance model.

    The dense inverse-Cholesky operator is constructed exactly once at worker
    entry for this PowerBin and then reused by the one-component control, all
    2601 fixed likelihood states, and the exact refit of the local minimum.
    """
    bid, galaxy, noise, good, valid_fraction = payload
    c = _WORKER
    min_good = int(c["min_good_pixels"])
    gp = np.flatnonzero(np.asarray(good, dtype=bool))
    if gp.size < min_good:
        raise RuntimeError(
            f"Bin {bid} has only {gp.size} good log pixels; configured minimum is {min_good}."
        )
    galaxy_n, noise_n, scale = _normalize_galaxy_for_ppxf(galaxy, noise, good)
    if not np.isfinite(scale):
        raise RuntimeError(f"Bin {bid} could not be normalized for pPXF.")

    cov_model = c["covariance_model"]
    effective_noise, inv_chol, whitening = covariance_calibration.effective_noise_and_whitener(
        cov_model,
        bin_id=int(bid),
        noise=noise_n,
        good=np.asarray(good, dtype=bool),
        eigen_floor=float(c["covariance_eigen_floor"]),
    )

    single = ppxf_utils.fit_single_losvd(
        templates=c["templates_single"],
        galaxy=galaxy_n,
        noise=effective_noise,
        velscale=c["velscale"],
        lam=c["lam"],
        lam_temp=c["lam_temp"],
        goodpixels=gp,
        start_velocity=float(c["single_start_velocity"]),
        start_sigma=float(c["sigma_start"]),
        velocity_bounds=tuple(c["single_velocity_bounds"]),
        sigma_bounds=tuple(c["sigma_bounds"]),
        degree=int(c["degree"]),
        mdegree=int(c["mdegree"]),
        regul=float(c["regul"]),
        noise_inv_cholesky=inv_chol,
        keep_full=True,
    )

    single_diagnostic = (
        "SUCCESS "
        f"(V={float(single.velocity[0]):.3f} km/s, "
        f"sigma={float(single.sigma[0]):.3f} km/s, "
        f"chi2_total={float(single.chi2_total):.6g})"
        if single.success
        else f"FAILED: {single.error_message or 'unknown one-component pPXF error'}"
    )

    cube = ppxf_grid.build_rh3_likelihood_cube(
        templates_two_component=c["templates_two"],
        component=c["component"],
        galaxy=galaxy_n,
        noise=effective_noise,
        velscale=c["velscale"],
        lam=c["lam"],
        lam_temp=c["lam_temp"],
        goodpixels=gp,
        va_grid=c["va_grid"],
        vb_grid=c["vb_grid"],
        fa_grid=c["fa_grid"],
        sigma_start_a=float(c["sigma_start"]),
        sigma_start_b=float(c["sigma_start"]),
        sigma_bounds=tuple(c["sigma_bounds"]),
        degree=int(c["degree"]),
        mdegree=int(c["mdegree"]),
        regul=float(c["regul"]),
        noise_inv_cholesky=inv_chol,
        sigma_boundary_tolerance_kms=float(c["sigma_boundary_tolerance"]),
    )
    if cube.best_index is None:
        lines = [
            f"Bin {bid}: every two-component grid state failed or was rejected.",
            f"One-component control: {single_diagnostic}",
            (
                "Two-component grid summary: "
                f"ppxf_failures={cube.n_ppxf_failures}, "
                f"fixed_velocity_mismatches={cube.n_fixed_velocity_mismatch}, "
                f"total_states={int(np.prod(cube.chi2_total.shape))}."
            ),
        ]
        if cube.first_failure_state is not None:
            va0, vb0, fa0 = cube.first_failure_state
            lines.append(
                "First failed state: "
                f"VA={va0:.3f} km/s, VB={vb0:.3f} km/s, fA={fa0:.3f}."
            )
            lines.append(
                "First state error: "
                f"{cube.first_failure_message or 'unknown failure'}"
            )
        if cube.failure_message_counts:
            lines.append("Most common two-component failure messages:")
            for msg, count in cube.failure_message_counts[:5]:
                lines.append(f"  {count}/{int(np.prod(cube.chi2_total.shape))}: {msg}")
        raise RuntimeError("\n".join(lines))

    ia, ib, jf = cube.best_index
    best_va = float(c["va_grid"][ia])
    best_vb = float(c["vb_grid"][ib])
    best_fa = float(c["fa_grid"][jf])
    best = ppxf_utils.fit_fixed_two_component_state(
        templates_two_component=c["templates_two"],
        component=c["component"],
        galaxy=galaxy_n,
        noise=effective_noise,
        velscale=c["velscale"],
        lam=c["lam"],
        lam_temp=c["lam_temp"],
        goodpixels=gp,
        velocity_a=best_va,
        velocity_b=best_vb,
        fraction_a=best_fa,
        start_sigma_a=float(c["sigma_start"]),
        start_sigma_b=float(c["sigma_start"]),
        sigma_bounds=tuple(c["sigma_bounds"]),
        degree=int(c["degree"]),
        mdegree=int(c["mdegree"]),
        regul=float(c["regul"]),
        noise_inv_cholesky=inv_chol,
        keep_full=True,
    )
    if not best.success:
        raise RuntimeError(f"Bin {bid}: exact refit of local best state failed: {best.error_message}")

    n_basis = int(c["n_basis"])
    w = np.asarray(best.weights, dtype=float)
    wa = float(np.sum(w[:n_basis]))
    wb = float(np.sum(w[n_basis:2 * n_basis]))
    achieved_fa = wa / (wa + wb) if (wa + wb) > 0 else np.nan

    return {
        "bin_id": int(bid),
        "schema_version": CHECKPOINT_SCHEMA_VERSION,
        "covariance_model": np.asarray(str(cov_model.name)),
        "covariance_hash": np.asarray(str(c["covariance_hash"])),
        "covariance_scale": np.asarray(float(cov_model.scale[int(bid)])),
        "covariance_regularized": np.asarray(bool(whitening.regularized) if whitening is not None else False),
        "normalization_scale": float(scale),
        "n_good_pixels": int(gp.size),
        "valid_fraction": np.asarray(valid_fraction, dtype=np.float32),
        "galaxy": np.asarray(galaxy_n, dtype=np.float32),
        # Preserve the normalized *formal* Script-2 uncertainty as the canonical
        # science noise vector.  For M1 only, pPXF receives s_i*noise; keeping the
        # formal vector here prevents downstream Script 5 from accidentally
        # applying s_i twice when it reconstructs the frozen covariance model.
        "noise": np.asarray(noise_n, dtype=np.float32),
        "ppxf_noise": np.asarray(effective_noise, dtype=np.float32),
        "good": np.asarray(good, dtype=np.uint8),
        "chi2_total": cube.chi2_total,
        "reduced_chi2": cube.reduced_chi2,
        "sigma_a": cube.sigma_a,
        "sigma_b": cube.sigma_b,
        "fit_status": cube.fit_status,
        "sigma_boundary": cube.sigma_boundary,
        "n_failures": int(cube.n_failures),
        "n_ppxf_failures": int(cube.n_ppxf_failures),
        "n_fixed_velocity_mismatch": int(cube.n_fixed_velocity_mismatch),
        "n_sigma_boundary": int(cube.n_sigma_boundary),
        "one_success": bool(single.success),
        "one_chi2_total": float(single.chi2_total),
        "one_reduced_chi2": float(single.reduced_chi2),
        "one_velocity": float(single.velocity[0]) if single.success else np.nan,
        "one_sigma": float(single.sigma[0]) if single.success else np.nan,
        "one_bestfit": (
            np.asarray(single.bestfit, dtype=np.float32)
            if single.bestfit is not None else np.full_like(galaxy_n, np.nan, dtype=np.float32)
        ),
        "one_weights": (
            np.asarray(single.weights, dtype=np.float32)
            if single.weights is not None else np.full(n_basis, np.nan, dtype=np.float32)
        ),
        "best_index": np.asarray(cube.best_index, dtype=np.int16),
        "best_va": best_va,
        "best_vb": best_vb,
        "best_fa_grid": best_fa,
        "best_fa_achieved": float(achieved_fa),
        "best_sigma_a": float(best.sigma[0]),
        "best_sigma_b": float(best.sigma[1]),
        "best_chi2_total": float(best.chi2_total),
        "best_reduced_chi2": float(best.reduced_chi2),
        "best_bestfit": np.asarray(best.bestfit, dtype=np.float32),
        "best_weights": np.asarray(best.weights, dtype=np.float32),
        "best_polyweights": (
            np.asarray(best.polyweights, dtype=np.float32)
            if best.polyweights is not None else np.empty(0, dtype=np.float32)
        ),
    }

def _checkpoint_path(checkpoints: Path, bid: int) -> Path:
    return checkpoints / f"bin_{int(bid):04d}.npz"


def _write_checkpoint(path: Path, result: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("wb") as handle:
        np.savez_compressed(handle, **result)
    os.replace(tmp, path)


def _checkpoint_is_valid(
    path: Path,
    expected_shape: tuple[int, int, int],
    covariance_hash: str | None = None,
) -> bool:
    if not path.is_file():
        return False
    try:
        with np.load(path, allow_pickle=False) as data:
            version = int(np.asarray(data["schema_version"]).reshape(-1)[0])
            if version != CHECKPOINT_SCHEMA_VERSION or tuple(data["chi2_total"].shape) != expected_shape:
                return False
            if covariance_hash is not None:
                saved = str(np.asarray(data["covariance_hash"]).reshape(-1)[0])
                if saved != str(covariance_hash):
                    return False
            return True
    except Exception:
        return False


def _read_checkpoint(path: Path) -> dict:
    with np.load(path, allow_pickle=False) as data:
        return {key: np.asarray(data[key]) for key in data.files}


def _save_template_product(run, prepared: templates.PreparedTemplates, lsf_model, cfg) -> Path:
    path = run.products_dir / "XSL_RH3_templates.npz"
    payload = {
        "templates": np.asarray(prepared.templates, dtype=np.float32),
        "wavelength": np.asarray(prepared.wavelength, dtype=np.float64),
        "velscale_kms": np.asarray(prepared.velscale),
        "native_template_fwhm": np.asarray(prepared.native_template_fwhm, dtype=np.float32),
        "target_fwhm": np.asarray(prepared.target_fwhm, dtype=np.float32),
        "convolution_fwhm": np.asarray(prepared.convolution_fwhm, dtype=np.float32),
        "normalization": np.asarray(prepared.normalization, dtype=np.float64),
        "normalization_range": np.asarray(prepared.normalization_range, dtype=float),
        "template_medium": np.asarray(prepared.template_medium),
        "lsf_source": np.asarray(str(lsf_model.source_path)),
        "lsf_empirical_observed_range": np.asarray([lsf_model.empirical_min, lsf_model.empirical_max]),
    }
    if prepared.age is not None:
        payload["age"] = np.asarray(prepared.age, dtype=float)
    if prepared.metallicity is not None:
        payload["metallicity"] = np.asarray(prepared.metallicity, dtype=float)
    np.savez_compressed(path, **payload)
    return path


def _consolidate(
    *,
    run,
    checkpoints: Path,
    nbin: int,
    grids: tuple[np.ndarray, np.ndarray, np.ndarray],
    lam: np.ndarray,
    source_table: Table,
    bin_map: np.ndarray,
    cfg,
    logger,
    covariance_model_name: str,
    covariance_hash: str,
):
    va, vb, fa = grids
    shape = (nbin, va.size, vb.size, fa.size)
    chi2 = np.full(shape, np.inf, dtype=np.float64)
    rchi2 = np.full(shape, np.nan, dtype=np.float32)
    siga = np.full(shape, np.nan, dtype=np.float32)
    sigb = np.full(shape, np.nan, dtype=np.float32)
    status = np.full(shape, ppxf_grid.FIT_STATUS_PPXF_FAILURE, dtype=np.int8)
    boundary = np.zeros(shape, dtype=np.uint8)

    nlog = lam.size
    galaxy = np.full((nbin, nlog), np.nan, dtype=np.float32)
    noise = np.full((nbin, nlog), np.nan, dtype=np.float32)
    ppxf_noise = np.full((nbin, nlog), np.nan, dtype=np.float32)
    good = np.zeros((nbin, nlog), dtype=np.uint8)
    valid_fraction = np.zeros((nbin, nlog), dtype=np.float32)
    one_model = np.full((nbin, nlog), np.nan, dtype=np.float32)
    two_model = np.full((nbin, nlog), np.nan, dtype=np.float32)

    one_success = np.zeros(nbin, dtype=bool)
    one_v = np.full(nbin, np.nan)
    one_s = np.full(nbin, np.nan)
    one_chi = np.full(nbin, np.nan)
    one_rchi = np.full(nbin, np.nan)
    best_va = np.full(nbin, np.nan)
    best_vb = np.full(nbin, np.nan)
    best_fa = np.full(nbin, np.nan)
    achieved_fa = np.full(nbin, np.nan)
    best_sa = np.full(nbin, np.nan)
    best_sb = np.full(nbin, np.nan)
    best_chi = np.full(nbin, np.nan)
    best_rchi = np.full(nbin, np.nan)
    n_good = np.zeros(nbin, dtype=int)
    n_fail = np.zeros(nbin, dtype=int)
    n_bound = np.zeros(nbin, dtype=int)

    one_weights = []
    two_weights = []
    polyweights = []
    scales = np.full(nbin, np.nan)
    covariance_scale = np.full(nbin, np.nan)
    covariance_regularized = np.zeros(nbin, dtype=bool)

    for bid in range(nbin):
        data = _read_checkpoint(_checkpoint_path(checkpoints, bid))
        chi2[bid] = data["chi2_total"]
        rchi2[bid] = data["reduced_chi2"]
        siga[bid] = data["sigma_a"]
        sigb[bid] = data["sigma_b"]
        status[bid] = data["fit_status"]
        boundary[bid] = data["sigma_boundary"]
        galaxy[bid] = data["galaxy"]
        noise[bid] = data["noise"]
        ppxf_noise[bid] = data["ppxf_noise"] if "ppxf_noise" in data else data["noise"]
        good[bid] = data["good"]
        valid_fraction[bid] = data["valid_fraction"]
        one_model[bid] = data["one_bestfit"]
        two_model[bid] = data["best_bestfit"]
        one_success[bid] = bool(np.asarray(data["one_success"]).reshape(-1)[0])
        one_v[bid] = float(data["one_velocity"])
        one_s[bid] = float(data["one_sigma"])
        one_chi[bid] = float(data["one_chi2_total"])
        one_rchi[bid] = float(data["one_reduced_chi2"])
        best_va[bid] = float(data["best_va"])
        best_vb[bid] = float(data["best_vb"])
        best_fa[bid] = float(data["best_fa_grid"])
        achieved_fa[bid] = float(data["best_fa_achieved"])
        best_sa[bid] = float(data["best_sigma_a"])
        best_sb[bid] = float(data["best_sigma_b"])
        best_chi[bid] = float(data["best_chi2_total"])
        best_rchi[bid] = float(data["best_reduced_chi2"])
        n_good[bid] = int(data["n_good_pixels"])
        n_fail[bid] = int(data["n_failures"])
        n_bound[bid] = int(data["n_sigma_boundary"])
        scales[bid] = float(data["normalization_scale"])
        covariance_scale[bid] = float(data["covariance_scale"])
        covariance_regularized[bid] = bool(np.asarray(data["covariance_regularized"]).reshape(-1)[0])
        one_weights.append(np.asarray(data["one_weights"], dtype=np.float32))
        two_weights.append(np.asarray(data["best_weights"], dtype=np.float32))
        polyweights.append(np.asarray(data["best_polyweights"], dtype=np.float32))

    cube_path = run.products_dir / "RH3_likelihood_cubes.npz"
    np.savez_compressed(
        cube_path,
        chi2_total=chi2,
        reduced_chi2=rchi2,
        sigma_A=siga,
        sigma_B=sigb,
        fit_status=status,
        sigma_boundary=boundary,
        VA_grid=va,
        VB_grid=vb,
        fA_grid=fa,
        one_velocity=one_v,
        one_sigma=one_s,
        one_chi2_total=one_chi,
        one_reduced_chi2=one_rchi,
        local_best_VA=best_va,
        local_best_VB=best_vb,
        local_best_fA=best_fa,
        local_best_fA_achieved=achieved_fa,
        local_best_sigma_A=best_sa,
        local_best_sigma_B=best_sb,
        local_best_chi2_total=best_chi,
        local_best_reduced_chi2=best_rchi,
        n_good_pixels=n_good,
        n_failed_states=n_fail,
        n_sigma_boundary_states=n_bound,
        covariance_scale=covariance_scale,
        covariance_regularized=covariance_regularized.astype(np.uint8),
        covariance_model_name=np.asarray(str(covariance_model_name)),
        covariance_model_hash=np.asarray(str(covariance_hash)),
    )

    selected_path = run.products_dir / "RH3_log_spectra_and_local_best_fits.npz"
    np.savez_compressed(
        selected_path,
        wavelength=np.asarray(lam, dtype=np.float64),
        galaxy=galaxy,
        # ``noise`` remains the normalized formal uncertainty. ``ppxf_noise`` is
        # included only for audit; it differs from ``noise`` for M1.
        noise=noise,
        ppxf_noise=ppxf_noise,
        good=good,
        valid_fraction=valid_fraction,
        normalization_scale=scales,
        covariance_scale=covariance_scale,
        covariance_regularized=covariance_regularized.astype(np.uint8),
        covariance_model_name=np.asarray(str(covariance_model_name)),
        covariance_model_hash=np.asarray(str(covariance_hash)),
        one_component_model=one_model,
        local_best_two_component_model=two_model,
        one_component_weights=np.stack(one_weights),
        local_best_two_component_weights=np.stack(two_weights),
        local_best_additive_polyweights=np.stack(polyweights) if polyweights and polyweights[0].size else np.empty((nbin, 0), dtype=np.float32),
    )

    summary = source_table.copy()
    summary["RH3_1C_SUCCESS"] = one_success
    summary["RH3_1C_V_KMS"] = one_v
    summary["RH3_1C_SIGMA_KMS"] = one_s
    summary["RH3_1C_CHI2"] = one_chi
    summary["RH3_LOCAL_VA_KMS"] = best_va
    summary["RH3_LOCAL_VB_KMS"] = best_vb
    summary["RH3_LOCAL_FA"] = best_fa
    summary["RH3_LOCAL_FA_ACHIEVED"] = achieved_fa
    summary["RH3_LOCAL_SIGMAA_KMS"] = best_sa
    summary["RH3_LOCAL_SIGMAB_KMS"] = best_sb
    summary["RH3_LOCAL_2C_CHI2"] = best_chi
    summary["RH3_LOCAL_DELTA_CHI2_1C_2C"] = one_chi - best_chi
    summary["RH3_GRID_FAILED_STATES"] = n_fail
    summary["RH3_GRID_SIGMA_BOUNDARY_STATES"] = n_bound
    summary["RH3_GRID_GOOD_LOG_PIXELS"] = n_good
    summary["RH3_COV_SCALE"] = covariance_scale
    summary["RH3_COV_PD_REGULARIZED"] = covariance_regularized
    table_path = run.products_dir / "RH3_local_likelihood_summary.ecsv"
    summary.write(table_path, format="ascii.ecsv", overwrite=True)

    delta_local = one_chi - best_chi
    plotting.plot_bin_value_map(
        bin_map, one_v, run.figures_dir / "RH3_single_component_velocity.png",
        title="RH3 one-component stellar velocity", colorbar_label="V (km/s)",
    )
    plotting.plot_bin_value_map(
        bin_map, one_s, run.figures_dir / "RH3_single_component_sigma.png",
        title="RH3 one-component stellar dispersion", colorbar_label="sigma (km/s)",
    )
    for values, filename, title, label in (
        (best_va, "RH3_local_best_VA.png", "Independent RH3 cube minimum: V_A", "V_A (km/s)"),
        (best_vb, "RH3_local_best_VB.png", "Independent RH3 cube minimum: V_B", "V_B (km/s)"),
        (best_fa, "RH3_local_best_fA.png", "Independent RH3 cube minimum: f_A", "f_A,RH3"),
        (best_sa, "RH3_local_best_sigmaA.png", "Independent RH3 cube minimum: sigma_A", "sigma_A (km/s)"),
        (best_sb, "RH3_local_best_sigmaB.png", "Independent RH3 cube minimum: sigma_B", "sigma_B (km/s)"),
        (delta_local, "RH3_local_one_vs_two_delta_chi2.png", "Local one- vs two-component improvement", "chi2_1comp - chi2_2comp"),
        (n_fail / float(va.size * vb.size * fa.size), "RH3_grid_failure_fraction.png", "RH3 pPXF grid-state failure fraction", "Failed state fraction"),
        (n_bound / float(va.size * vb.size * fa.size), "RH3_sigma_boundary_fraction.png", "RH3 grid states near a sigma bound", "Boundary state fraction"),
    ):
        plotting.plot_bin_value_map(bin_map, values, run.figures_dir / filename, title=title, colorbar_label=label)

    if bool(getattr(cfg, "SAVE_PER_BIN_DIAGNOSTICS", True)):
        per_bin = run.figures_dir / "RH3_likelihood_bins"
        per_bin.mkdir(parents=True, exist_ok=True)
        logger.info("Writing %d per-bin RH3 likelihood diagnostics", nbin)
        for bid in range(nbin):
            plotting.plot_rh3_likelihood_bin(
                va_grid=va,
                vb_grid=vb,
                fa_grid=fa,
                chi2_total=chi2[bid],
                wavelength=lam,
                galaxy=galaxy[bid],
                noise=noise[bid],
                good=good[bid].astype(bool),
                one_model=one_model[bid],
                two_model=two_model[bid],
                output_path=per_bin / f"bin_{bid:04d}.png",
                bin_id=bid,
                one_chi2=one_chi[bid],
                two_chi2=best_chi[bid],
                best_state=(best_va[bid], best_vb[bid], best_fa[bid], best_sa[bid], best_sb[bid]),
            )

    return cube_path, selected_path, table_path


def main() -> int:
    args = _parser().parse_args()
    cfg = crd.load_config(args.config, validate=True, strict_paths=True)

    if args.workers is not None and args.workers < 1:
        raise ValueError("--workers must be >= 1.")
    physical = os.cpu_count() or 1
    if args.workers is not None:
        workers = int(args.workers)
    elif getattr(cfg, "N_WORKERS", None) is not None:
        workers = int(cfg.N_WORKERS)
    else:
        workers = max(1, physical - 1)

    if args.resume:
        run_dir = _find_resume_run(cfg, args.run_name)
        run = _existing_run_context(run_dir)
        resume_state_path = run.metadata_dir / "script03_resume_state.json"
        state = dict(io.read_json(resume_state_path))
        current_hash = _sha256(Path(args.config).expanduser().resolve())
        if current_hash != state.get("config_sha256"):
            raise RuntimeError(
                "The current config file differs from the config used to create this "
                "Script-3 run. Resume is refused because changing a wavelength range, "
                "grid, polynomial, or sigma bound would mix incompatible checkpoints."
            )
        source_run = Path(state["source_script2_run"]).expanduser().resolve()
        if args.script2_run is not None and source_run != Path(args.script2_run).expanduser().resolve():
            raise RuntimeError("--script2-run disagrees with the source run stored in the resume state.")
        state_script1 = state.get("source_script1_run")
        script1_run = (
            Path(state_script1).expanduser().resolve()
            if state_script1 not in (None, "")
            else None
        )
        resume_script1_override = None
        if args.script1_run is not None:
            explicit_script1 = Path(args.script1_run).expanduser().resolve()
            if script1_run is not None and script1_run != explicit_script1:
                existing_checkpoints = list((run.run_dir / "checkpoints").glob("bin_*.npz"))
                if existing_checkpoints:
                    raise RuntimeError(
                        "--script1-run disagrees with the source Script-1 run stored in the "
                        "resume state, and completed per-bin checkpoints already exist. "
                        "Changing Script-1 provenance after fitting has begun could mix "
                        "different LSF/wavelength assumptions. Start a new Script-3 run instead."
                    )
                resume_script1_override = script1_run
            script1_run = explicit_script1
        if script1_run is not None and not _is_usable_script1_run(script1_run):
            raise FileNotFoundError(
                "The Script-1 run stored for this resume is unavailable or incomplete: "
                f"{script1_run}."
            )
    else:
        source_run = _find_script2_run(cfg, args.script2_run)
        script1_run = None
        resume_script1_override = None
        run = _new_run(cfg, args.run_name)

    logger = _setup_logger(run)
    start_time = time.perf_counter()
    quality_flags: list[str] = []
    checkpoints = run.run_dir / "checkpoints"
    checkpoints.mkdir(parents=True, exist_ok=True)

    try:
        _section(logger, "CRD_DAP SCRIPT 3: BUILD RH3 PROFILE-LIKELIHOOD CUBES", "=")
        logger.info("Run directory: %s", run.run_dir)
        logger.info("Source Script-2 run: %s", source_run)
        logger.info("Python: %s", sys.version.replace("\n", " "))
        logger.info("Platform: %s", platform.platform())
        logger.info("Target: %s | redshift=%.8f", cfg.TARGET_NAME, float(cfg.REDSHIFT))
        logger.info(
            "Parallelism: workers=%d process(es), os.cpu_count()=%d; BLAS/OpenMP threads capped at 1 per worker",
            workers, physical,
        )
        if workers > physical:
            logger.warning(
                "WORKER_OVERSUBSCRIPTION | --workers=%d exceeds os.cpu_count()=%d; this may reduce performance.",
                workers, physical,
            )

        products = _required_script2_products(source_run)
        if not _is_complete_script2_run(source_run):
            raise FileNotFoundError(f"Incomplete Script-2 run: {source_run}")
        s02_manifest = _read_script2_manifest(source_run)
        if args.resume and script1_run is not None:
            manifest_value = s02_manifest.get("source_script1_run")
            manifest_script1_run = (
                Path(manifest_value).expanduser().resolve()
                if manifest_value not in (None, "")
                else None
            )
        else:
            script1_run, manifest_script1_run = _resolve_script1_run(
                s02_manifest, args.script1_run
            )

        logger.info("Source Script-1 run: %s", script1_run)
        if (
            args.script1_run is not None
            and manifest_script1_run is not None
            and script1_run != manifest_script1_run
        ):
            logger.warning(
                "SCRIPT1_RUN_PATH_OVERRIDE | Script-2 manifest records %s, but Script 3 "
                "will use explicit --script1-run %s. This is appropriate when the same "
                "Script-1 run directory was moved or renamed; verify that the products "
                "really correspond to the Script-2 provenance.",
                manifest_script1_run, script1_run,
            )
        if args.resume and resume_script1_override is not None:
            logger.warning(
                "RESUME_SCRIPT1_PATH_OVERRIDE | Resume state recorded %s, but no completed "
                "per-bin checkpoints exist, so the explicit --script1-run path %s is safe "
                "to adopt. The resume metadata will be updated before fitting begins.",
                resume_script1_override, script1_run,
            )
            state["source_script1_run"] = str(script1_run)
            io.write_json(state, run.metadata_dir / "script03_resume_state.json")

        if not args.resume:
            state = {
                "script": "03_build_RH3_likelihood_cubes",
                "target": str(cfg.TARGET_NAME),
                "source_script2_run": str(source_run),
                "source_script1_run": str(script1_run),
                "source_reduction_manifest": s02_manifest.get("source_reduction_manifest"),
                "source_reduction_manifest_sha256": s02_manifest.get(
                    "source_reduction_manifest_sha256"
                ),
                "config_path": str(Path(args.config).expanduser().resolve()),
                "config_sha256": _sha256(Path(args.config).expanduser().resolve()),
                "checkpoint_schema_version": CHECKPOINT_SCHEMA_VERSION,
                "created_utc": datetime.now(timezone.utc).isoformat(),
            }
            io.write_json(state, run.metadata_dir / "script03_resume_state.json")

        _section(logger, "1. Load Script-2 master spectra and Script-1 LSF")
        wave_obs, flux, unc, good_native, bin_map, source_table = _load_script2_arrays(products)
        nbin = flux.shape[0]
        header = _load_prepared_rh3_header(script1_run)
        science_medium = _resolve_science_medium(cfg, header)
        template_medium = str(cfg.TEMPLATE_WAVELENGTH_MEDIUM).strip().lower()
        if template_medium not in {"air", "vacuum"}:
            raise ValueError("TEMPLATE_WAVELENGTH_MEDIUM must be 'air' or 'vacuum'.")
        lsf_path = templates.discover_script1_lsf_product(
            script1_run,
            arm="RH3",
            explicit=getattr(cfg, "RH3_LSF_PRODUCT", None),
        )
        lsf_model = templates.load_saved_lsf_model(
            lsf_path, polynomial_order=int(cfg.LSF_MODEL_WAVELENGTH_ORDER)
        )
        logger.info(
            "RH3 master spectra: nbin=%d, native channels=%d, observed wavelength=%.1f--%.1f A",
            nbin, wave_obs.size, float(np.nanmin(wave_obs)), float(np.nanmax(wave_obs)),
        )
        logger.info(
            "Wavelength media: science=%s, XSL/templates=%s | Script-1 LSF=%s | empirical support=%.1f--%.1f A observed",
            science_medium, template_medium, lsf_path, lsf_model.empirical_min, lsf_model.empirical_max,
        )

        fit_range = tuple(float(x) for x in cfg.RH3_FIT_REST_RANGE_ANGSTROM)
        rest_template_native = templates.convert_wavelength_medium(
            wave_obs / (1.0 + float(cfg.REDSHIFT)), science_medium, template_medium
        )
        velscale = _determine_velscale(rest_template_native, cfg)
        logger.info(
            "RH3 fitting interval: %.1f--%.1f A rest (%s) | log-grid velocity scale=%.3f km/s/pixel",
            fit_range[0], fit_range[1], template_medium, velscale,
        )
        observed_masks, observed_mask_provenance = _resolve_observed_masks(cfg, science_medium)
        rest_masks = [(float(lo), float(hi)) for lo, hi in getattr(cfg, "RH3_MASK_REST_RANGES_ANGSTROM", [])]

        inherited_atmosphere = (
            (s02_manifest.get("source_script1_atmospheric_masks") or {}).get("RH3")
        )
        if inherited_atmosphere:
            logger.info(
                "RH3 CRD_DRP atmospheric mask is already inherited through Script-2 "
                "RH3_GOOD: source=%s | masked=%s/%s native channels",
                inherited_atmosphere.get("source"),
                inherited_atmosphere.get("n_masked_pixels"),
                inherited_atmosphere.get("n_native_pixels"),
            )
        else:
            logger.warning(
                "NO_CRD_DRP_ATMOSPHERIC_PROVENANCE | Script-2 manifest does not record "
                "an inherited RH3 CRD_DRP atmospheric mask. This is acceptable only for "
                "legacy/testing runs."
            )

        if observed_masks:
            logger.info(
                "RH3 additional fixed observed-frame masks (%s science wavelength medium): %s",
                science_medium, observed_masks,
            )
        if rest_masks:
            logger.info("RH3 fixed rest-frame masks (%s template medium): %s", template_medium, rest_masks)
        if float(cfg.RH3_SIGMA_START_KMS) < 3.0 * velscale:
            logger.warning(
                "SIGMA_START_UNDERSAMPLED | RH3_SIGMA_START_KMS=%.2f is <3 log pixels (%.2f km/s); pPXF may be less stable.",
                float(cfg.RH3_SIGMA_START_KMS), 3.0 * velscale,
            )

        _section(logger, "2. Construct the fixed log-wavelength experiment")
        log_wave, log_edges, gal_all, noise_all, good_all, valid_frac_all = _prepare_log_binned_spectra(
            wave_obs=wave_obs,
            flux=flux,
            uncertainty=unc,
            good_native=good_native,
            redshift=float(cfg.REDSHIFT),
            science_medium=science_medium,
            template_medium=template_medium,
            fit_range=fit_range,
            velscale=velscale,
            cfg=cfg,
            observed_masks=observed_masks,
            rest_masks=rest_masks,
            lsf_model=lsf_model,
        )
        good_counts = np.sum(good_all, axis=1)
        logger.info(
            "Log spectra: %d pixels/bin; good-pixel count median=%d, min=%d, max=%d",
            log_wave.size, int(np.median(good_counts)), int(np.min(good_counts)), int(np.max(good_counts)),
        )

        # The science representation intentionally retains NaNs in rejected
        # log-grid samples. pPXF, however, requires its complete galaxy/noise
        # vectors to be finite and its complete noise vector to be positive even
        # when those samples are absent from goodpixels. Verify here that every
        # *fitted* sample is already valid; ppxf_utils will replace only invalid
        # excluded samples in private copies passed to pPXF.
        good_bool = np.asarray(good_all, dtype=bool)
        invalid_good_galaxy = good_bool & ~np.isfinite(gal_all)
        invalid_good_noise = good_bool & (~np.isfinite(noise_all) | (noise_all <= 0))
        if np.any(invalid_good_galaxy) or np.any(invalid_good_noise):
            bad_bins = np.flatnonzero(
                np.any(invalid_good_galaxy | invalid_good_noise, axis=1)
            )
            raise RuntimeError(
                "SCRIPT03_INVALID_GOODPIXELS | The fixed good-pixel mask includes "
                "non-finite galaxy flux or non-finite/non-positive noise. This is a "
                "science-mask inconsistency and must not be repaired with pPXF placeholders. "
                f"Affected bins={bad_bins.size}/{nbin}; examples={bad_bins[:12].tolist()}"
            )

        excluded = ~good_bool
        galaxy_placeholders = excluded & ~np.isfinite(gal_all)
        noise_placeholders = excluded & (~np.isfinite(noise_all) | (noise_all <= 0))
        placeholder_union = galaxy_placeholders | noise_placeholders
        placeholder_counts = np.sum(placeholder_union, axis=1)
        placeholder_bins = np.flatnonzero(placeholder_counts > 0)
        if placeholder_bins.size:
            logger.info(
                "pPXF excluded-pixel API sanitization required for %d/%d bins: "
                "galaxy placeholders=%d, noise placeholders=%d; median/max affected "
                "samples per affected bin=%d/%d. Only private pPXF input copies are "
                "filled; these samples remain outside goodpixels and the saved science "
                "arrays retain their original NaN/mask representation.",
                int(placeholder_bins.size), nbin,
                int(np.sum(galaxy_placeholders)), int(np.sum(noise_placeholders)),
                int(np.median(placeholder_counts[placeholder_bins])),
                int(np.max(placeholder_counts[placeholder_bins])),
            )
        else:
            logger.info(
                "pPXF excluded-pixel API sanitization: no placeholder fills are required."
            )

        too_short = np.flatnonzero(good_counts < int(cfg.RH3_MIN_GOOD_LOG_PIXELS))
        if too_short.size:
            raise RuntimeError(
                f"{too_short.size}/{nbin} bins have fewer than RH3_MIN_GOOD_LOG_PIXELS="
                f"{int(cfg.RH3_MIN_GOOD_LOG_PIXELS)} usable log samples. Example bins: {too_short[:12].tolist()}"
            )

        _section(logger, "3. Prepare full XSL SSP basis at the empirical RH3 LSF")
        va = ppxf_grid.uniform_grid(cfg.RH3_VA_MIN_KMS, cfg.RH3_VA_MAX_KMS, cfg.RH3_VA_N)
        vb = ppxf_grid.uniform_grid(cfg.RH3_VB_MIN_KMS, cfg.RH3_VB_MAX_KMS, cfg.RH3_VB_N)
        fa = ppxf_grid.fraction_grid(cfg.RH3_FA_MIN, cfg.RH3_FA_MAX, cfg.RH3_FA_STEP)
        support = _script03_velocity_support(va, vb, cfg, velscale)
        padding_kms = float(support["total_padding_kms"])

        def target_fwhm(rest_template_wave):
            return templates.observed_lsf_to_template_rest_fwhm(
                rest_template_wave,
                lsf_model=lsf_model,
                redshift=float(cfg.REDSHIFT),
                science_medium=science_medium,
                template_medium=template_medium,
            )

        prepared = templates.prepare_xsl_rh3_templates(
            xsl_path=cfg.XSL_TEMPLATE_LIBRARY,
            fit_range=fit_range,
            velscale=velscale,
            target_fwhm_rest=target_fwhm,
            template_medium=template_medium,
            velocity_padding_kms=padding_kms,
        )
        template_path = _save_template_product(run, prepared, lsf_model, cfg)
        n_basis = prepared.templates.shape[1]
        templates_two = np.column_stack([prepared.templates, prepared.templates])
        component = np.concatenate([
            np.zeros(n_basis, dtype=int),
            np.ones(n_basis, dtype=int),
        ])
        logger.info(
            "Prepared XSL RH3 basis: %d SSPs -> %d duplicated two-component templates; template pixels=%d",
            n_basis, templates_two.shape[1], prepared.templates.shape[0],
        )
        logger.info(
            "f_A,RH3 convention: each SSP independently normalized to unit mean stellar flux over %.1f--%.1f A rest (%s); additive polynomial is excluded from f_A",
            fit_range[0], fit_range[1], template_medium,
        )
        logger.info("Prepared template product: %s", template_path)

        # pPXF only requires equal velscale, not equal starting wavelength.  Check
        # the actual logarithmic spacing explicitly before launching a long run.
        temp_velscale = templates.C_KMS * np.nanmedian(np.diff(np.log(prepared.wavelength)))
        if not np.isclose(temp_velscale, velscale, rtol=1.0e-7, atol=1.0e-7):
            raise RuntimeError(
                f"Template/galaxy velocity-scale mismatch: templates={temp_velscale}, galaxy={velscale} km/s."
            )

        preflight_bins = _run_ppxf_template_coverage_preflight(
            logger=logger, prepared=prepared, templates_two=templates_two, component=component,
            gal_all=gal_all, noise_all=noise_all, good_all=good_all, log_wave=log_wave,
            va=va, vb=vb, fa=fa, support=support, cfg=cfg, velscale=velscale,
        )

        if not bool(getattr(cfg, "RH3_COVARIANCE_ENABLE", True)):
            raise RuntimeError(
                "RH3_COVARIANCE_DISABLED | Production Script 3 now requires empirical covariance "
                "calibration before likelihood cubes are generated. Set RH3_COVARIANCE_ENABLE=True."
            )

        _section(logger, "4. Calibrate and select the RH3 spectral covariance model")
        if args.resume and covariance_calibration.calibration_products_complete(run):
            calibration = covariance_calibration.load_calibration_run(run)
            logger.info(
                "Resume covariance calibration: reusing selected %s model from %s",
                calibration.selected_model_name,
                calibration.selection_json,
            )
        else:
            existing_grid_checkpoints = list(checkpoints.glob("bin_*.npz"))
            if args.resume and existing_grid_checkpoints:
                raise RuntimeError(
                    "RESUME_COVARIANCE_PRODUCTS_MISSING | Per-bin likelihood checkpoints exist but the "
                    "saved covariance-calibration decision is incomplete. Refusing to recalibrate around "
                    "already-computed likelihood states because that could mix incompatible chi-square metrics."
                )
            calibration = covariance_calibration.run_script03_covariance_calibration(
                cfg=cfg,
                logger=logger,
                run=run,
                source_table=source_table,
                bin_map=bin_map,
                templates_single=prepared.templates,
                templates_two=templates_two,
                component=component,
                galaxy_all=gal_all,
                noise_all=noise_all,
                good_all=good_all,
                wavelength=log_wave,
                lam_temp=prepared.wavelength,
                velscale=velscale,
                va_grid=va,
                vb_grid=vb,
                fa_grid=fa,
                sigma_bounds=(float(cfg.RH3_SIGMA_MIN_KMS), float(cfg.RH3_SIGMA_MAX_KMS)),
                single_velocity_bounds=tuple(support["single_velocity_bounds_kms"]),
                degree=int(cfg.RH3_DEGREE),
                mdegree=int(cfg.RH3_MDEGREE),
                regul=float(cfg.RH3_REGUL),
                sigma_start=float(cfg.RH3_SIGMA_START_KMS),
                sigma_boundary_tolerance=float(cfg.RH3_SIGMA_BOUNDARY_WARNING_KMS),
                workers=workers,
                plotting_module=covariance_plotting,
            )

        selected_covariance = calibration.selected_model
        covariance_hash = calibration.selected_model_hash
        logger.info(
            "Frozen production covariance: %s | hash=%s | iterations=%d | Requirement A=%s",
            calibration.selected_model_name,
            covariance_hash[:16],
            selected_covariance.n_iterations,
            "PASS" if selected_covariance.requirement_a_pass else "FAIL",
        )
        if not selected_covariance.requirement_a_pass:
            raise RuntimeError(
                "SELECTED_COVARIANCE_REQUIREMENT_A_FAILED | Internal model-selection inconsistency; "
                "production likelihood fitting is blocked."
            )
        if int(selected_covariance.numerical_regularization_count) > 0:
            quality_flags.append("RH3_COVARIANCE_PD_REGULARIZATION")
            logger.warning(
                "RH3_COVARIANCE_PD_REGULARIZATION | The selected covariance model required "
                "positive-definiteness eigenvalue flooring for %d calibration whiteners. "
                "Inspect covariance diagnostics before publication.",
                int(selected_covariance.numerical_regularization_count),
            )

        state["covariance_model"] = str(calibration.selected_model_name)
        state["covariance_model_hash"] = str(covariance_hash)
        state["covariance_selection_json"] = str(calibration.selection_json)
        io.write_json(state, run.metadata_dir / "script03_resume_state.json")

        _section(logger, "5. Evaluate covariance-calibrated 3-D profile-likelihood cubes")
        expected_shape = (va.size, vb.size, fa.size)
        completed = [
            bid for bid in range(nbin)
            if _checkpoint_is_valid(_checkpoint_path(checkpoints, bid), expected_shape, covariance_hash)
        ]
        pending = [bid for bid in range(nbin) if bid not in set(completed)]
        if args.resume:
            logger.info("Resume scan: %d/%d covariance-compatible bin checkpoints found; %d bins remain", len(completed), nbin, len(pending))
        else:
            logger.info("Production grid per bin: %d x %d x %d = %d exact covariance-aware states", va.size, vb.size, fa.size, int(np.prod(expected_shape)))

        constants = {
            "templates_single": prepared.templates,
            "templates_two": templates_two,
            "component": component,
            "n_basis": n_basis,
            "velscale": velscale,
            "lam": log_wave,
            "lam_temp": prepared.wavelength,
            "va_grid": va,
            "vb_grid": vb,
            "fa_grid": fa,
            "sigma_start": float(cfg.RH3_SIGMA_START_KMS),
            "sigma_bounds": (float(cfg.RH3_SIGMA_MIN_KMS), float(cfg.RH3_SIGMA_MAX_KMS)),
            "sigma_boundary_tolerance": float(cfg.RH3_SIGMA_BOUNDARY_WARNING_KMS),
            "single_start_velocity": float(getattr(cfg, "RH3_SINGLE_VELOCITY_START_KMS", 0.0)),
            "single_velocity_bounds": tuple(support["single_velocity_bounds_kms"]),
            "degree": int(cfg.RH3_DEGREE),
            "mdegree": int(cfg.RH3_MDEGREE),
            "regul": float(cfg.RH3_REGUL),
            "min_good_pixels": int(cfg.RH3_MIN_GOOD_LOG_PIXELS),
            "covariance_model": selected_covariance,
            "covariance_hash": str(covariance_hash),
            "covariance_eigen_floor": float(getattr(cfg, "RH3_COVARIANCE_EIGEN_FLOOR", 1.0e-8)),
        }

        if pending:
            fit_start = time.perf_counter()
            last_completion = fit_start
            done_now = 0
            spinner_chars = "|/-\\"
            spinner_index = 0
            status_width = 0

            logger.info(
                "Interactive heartbeat: terminal status refresh every %.0f s while waiting for completed bins; progress is measured at the PowerBin checkpoint level.",
                PROGRESS_HEARTBEAT_SECONDS,
            )

            with ProcessPoolExecutor(
                max_workers=workers,
                initializer=_worker_init,
                initargs=(constants,),
            ) as pool:
                future_to_bid = {
                    pool.submit(
                        _fit_one_bin_worker,
                        (bid, gal_all[bid], noise_all[bid], good_all[bid], valid_frac_all[bid]),
                    ): bid
                    for bid in pending
                }
                outstanding = set(future_to_bid)

                # Draw immediately so the user sees activity even before the
                # first expensive 2601-state PowerBin has completed.
                now = time.perf_counter()
                status_width = _write_terminal_status(
                    _fit_status_line(
                        spinner=spinner_chars[spinner_index % len(spinner_chars)],
                        completed_total=len(completed) + done_now,
                        nbin=nbin,
                        workers=workers,
                        elapsed=now - fit_start,
                        since_last_completion=now - last_completion,
                        eta_seconds=None,
                    ),
                    status_width,
                )

                while outstanding:
                    done_set, outstanding = wait(
                        outstanding,
                        timeout=PROGRESS_HEARTBEAT_SECONDS,
                        return_when=FIRST_COMPLETED,
                    )
                    now = time.perf_counter()

                    if not done_set:
                        spinner_index += 1
                        if done_now > 0:
                            seconds_per_bin = (now - fit_start) / done_now
                            eta_seconds = (len(pending) - done_now) * seconds_per_bin
                        else:
                            eta_seconds = None
                        status_width = _write_terminal_status(
                            _fit_status_line(
                                spinner=spinner_chars[spinner_index % len(spinner_chars)],
                                completed_total=len(completed) + done_now,
                                nbin=nbin,
                                workers=workers,
                                elapsed=now - fit_start,
                                since_last_completion=now - last_completion,
                                eta_seconds=eta_seconds,
                            ),
                            status_width,
                        )
                        continue

                    for future in done_set:
                        bid = future_to_bid[future]
                        _clear_terminal_status(status_width)
                        status_width = 0
                        try:
                            result = future.result()
                        except Exception:
                            logger.exception(
                                "Bin %d failed before checkpoint completion; aborting Script 3 so the run can be resumed.",
                                bid,
                            )
                            for f in outstanding:
                                f.cancel()
                            raise

                        _write_checkpoint(_checkpoint_path(checkpoints, bid), result)
                        done_now += 1
                        last_completion = time.perf_counter()
                        elapsed = last_completion - fit_start
                        rate = elapsed / done_now
                        remain = (len(pending) - done_now) * rate
                        logger.info(
                            "Bin %d complete | %d/%d new bins | failed states=%d/%d (pPXF=%d, fixed-V mismatch=%d) | sigma-boundary states=%d | elapsed=%.1f min | ETA=%.1f min",
                            bid,
                            done_now,
                            len(pending),
                            int(result["n_failures"]),
                            int(np.prod(expected_shape)),
                            int(result["n_ppxf_failures"]),
                            int(result["n_fixed_velocity_mismatch"]),
                            int(result["n_sigma_boundary"]),
                            elapsed / 60.0,
                            remain / 60.0,
                        )

                    if outstanding:
                        spinner_index += 1
                        now = time.perf_counter()
                        seconds_per_bin = (now - fit_start) / done_now if done_now else np.nan
                        eta_seconds = (len(pending) - done_now) * seconds_per_bin if done_now else None
                        status_width = _write_terminal_status(
                            _fit_status_line(
                                spinner=spinner_chars[spinner_index % len(spinner_chars)],
                                completed_total=len(completed) + done_now,
                                nbin=nbin,
                                workers=workers,
                                elapsed=now - fit_start,
                                since_last_completion=now - last_completion,
                                eta_seconds=eta_seconds,
                            ),
                            status_width,
                        )

                _clear_terminal_status(status_width)
        else:
            logger.info("All %d bin checkpoints are already complete; proceeding directly to consolidation.", nbin)

        # Final safety check: every bin checkpoint must be valid before we create
        # products that downstream scripts may interpret as complete.
        missing = [bid for bid in range(nbin) if not _checkpoint_is_valid(_checkpoint_path(checkpoints, bid), expected_shape, covariance_hash)]
        if missing:
            raise RuntimeError(f"Missing/invalid checkpoints after fitting: {missing[:20]}")

        _section(logger, "6. Consolidate likelihood products and diagnostics")
        cube_path, selected_path, table_path = _consolidate(
            run=run,
            checkpoints=checkpoints,
            nbin=nbin,
            grids=(va, vb, fa),
            lam=log_wave,
            source_table=source_table,
            bin_map=bin_map,
            cfg=cfg,
            logger=logger,
            covariance_model_name=str(calibration.selected_model_name),
            covariance_hash=str(covariance_hash),
        )

        # QC summary from the consolidated table is intentionally local-minimum
        # only. Script 4, not Script 3, decides which component labeling/basin is
        # physically coherent across the galaxy.
        summary = Table.read(table_path, format="ascii.ecsv")
        frac_diff = np.abs(np.asarray(summary["RH3_LOCAL_FA_ACHIEVED"], dtype=float) - np.asarray(summary["RH3_LOCAL_FA"], dtype=float))
        if np.nanmax(frac_diff) > float(cfg.RH3_FRACTION_CONSTRAINT_TOLERANCE):
            quality_flags.append("RH3_FRACTION_CONSTRAINT_WARNING")
            logger.warning(
                "RH3_FRACTION_CONSTRAINT_WARNING | max |achieved-grid f_A|=%.3g exceeds tolerance %.3g",
                float(np.nanmax(frac_diff)), float(cfg.RH3_FRACTION_CONSTRAINT_TOLERANCE),
            )
        if np.any(~np.asarray(summary["RH3_1C_SUCCESS"], dtype=bool)):
            quality_flags.append("RH3_ONE_COMPONENT_FAILURES")
        if np.any(np.asarray(summary["RH3_GRID_FAILED_STATES"], dtype=int) > 0):
            quality_flags.append("RH3_GRID_FIT_FAILURES")
        if np.any(np.asarray(summary["RH3_GRID_SIGMA_BOUNDARY_STATES"], dtype=int) > 0):
            quality_flags.append("RH3_SIGMA_BOUNDARY_STATES")

        manifest = {
            "script": "03_build_RH3_likelihood_cubes",
            "target": str(cfg.TARGET_NAME),
            "redshift": float(cfg.REDSHIFT),
            "source_script2_run": str(source_run),
            "source_script1_run": str(script1_run),
            "source_script2_quality_flags": list(s02_manifest.get("quality_flags", [])),
            "n_bins": int(nbin),
            "workers": int(workers),
            "thread_caps": {
                name: os.environ.get(name) for name in (
                    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                    "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
                )
            },
            "fit_rest_range_angstrom": [float(fit_range[0]), float(fit_range[1])],
            "science_wavelength_medium": science_medium,
            "template_wavelength_medium": template_medium,
            "observed_mask_wavelength_medium": science_medium,
            "source_reduction_manifest": s02_manifest.get("source_reduction_manifest"),
            "source_reduction_manifest_sha256": s02_manifest.get(
                "source_reduction_manifest_sha256"
            ),
            "inherited_rh3_atmospheric_mask": inherited_atmosphere,
            "atmospheric_mask_application_stage": (
                "CRD_DRP -> Script 1 GOODMASK/GOODWAVE -> Script 2 RH3_GOOD"
            ),
            "configured_observed_mask_ranges_angstrom": [
                [float(lo), float(hi)] for lo, hi in getattr(cfg, "RH3_MASK_OBSERVED_RANGES_ANGSTROM", [])
            ],
            # Retain these legacy manifest keys as empty values so older readers
            # fail softly while production no longer uses the 03d file pathway.
            "atmospheric_mask_file": None,
            "atmospheric_mask_file_ranges_angstrom": [],
            "effective_observed_mask_ranges_angstrom": observed_mask_provenance[
                "effective_additional_ranges"
            ],
            "configured_rest_mask_ranges_angstrom": [
                [float(lo), float(hi)] for lo, hi in getattr(cfg, "RH3_MASK_REST_RANGES_ANGSTROM", [])
            ],
            "velscale_kms": float(velscale),
            "template_support": {
                "two_component_max_abs_velocity_kms": float(support["two_component_max_abs_kms"]),
                "single_velocity_bounds_kms": [float(x) for x in support["single_velocity_bounds_kms"]],
                "dispersion_padding_sigma_multiple": float(cfg.RH3_TEMPLATE_PADDING_SIGMA),
                "dispersion_margin_kms": float(support["dispersion_margin_kms"]),
                "edge_safety_pixels": int(TEMPLATE_EDGE_SAFETY_PIXELS),
                "edge_safety_kms": float(support["edge_safety_kms"]),
                "total_velocity_padding_kms": float(support["total_padding_kms"]),
                "prepared_template_range_angstrom": [float(prepared.wavelength[0]), float(prepared.wavelength[-1])],
                "ppxf_preflight_bins": [int(x) for x in preflight_bins],
                "ppxf_preflight_passed": True,
            },
            "grid": {
                "VA_kms": va.tolist(),
                "VB_kms": vb.tolist(),
                "fA": fa.tolist(),
            },
            "sigma_start_kms": float(cfg.RH3_SIGMA_START_KMS),
            "sigma_bounds_kms": [float(cfg.RH3_SIGMA_MIN_KMS), float(cfg.RH3_SIGMA_MAX_KMS)],
            "degree": int(cfg.RH3_DEGREE),
            "mdegree": int(cfg.RH3_MDEGREE),
            "regul": float(cfg.RH3_REGUL),
            "profile_likelihood": True,
            "fixed_coordinates": ["V_A", "V_B", "f_A_RH3"],
            "profiled_nuisance": ["sigma_A", "sigma_B", "SSP weights", "additive polynomial"],
            "chi2_statistic": (
                "total covariance-aware chi-square r^T C^-1 r on one fixed good-pixel set per bin; "
                "M1 reduces to the diagonal special case"
            ),
            "noise_model": {
                "selected_model": str(calibration.selected_model_name),
                "selected_model_hash": str(covariance_hash),
                "per_bin_scale": "s_i calibrated from robust normalized pPXF residual scatter",
                "correlation": "empirical lag correlation selected by M1--M4 residual adequacy + representative full-grid stability",
                "cached_ppxf_keyword": "noise_inv_cholesky",
                "required_ppxf_version": str(getattr(cfg, "RH3_COVARIANCE_REQUIRED_PPXF_VERSION", "9.4.8")),
                "bootstrap_n": int(getattr(cfg, "RH3_COVARIANCE_BOOTSTRAP_N", 2000)),
                "bootstrap_confidence": float(getattr(cfg, "RH3_COVARIANCE_BOOTSTRAP_CONFIDENCE", 0.95)),
                "max_lag_pixels": int(getattr(cfg, "RH3_COVARIANCE_MAX_LAG", 20)),
                "wavelength_blocks_M3_M4": int(getattr(cfg, "RH3_COVARIANCE_WAVELENGTH_BLOCKS", 3)),
                "convergence_tolerance": float(getattr(cfg, "RH3_COVARIANCE_CONVERGENCE_TOL", 0.01)),
                "max_iterations": int(getattr(cfg, "RH3_COVARIANCE_MAX_ITER", 5)),
                "representative_radial_bins": int(getattr(cfg, "RH3_COVARIANCE_VALIDATION_RADIAL_BINS", 12)),
                "representative_total_bins": int(len(calibration.representative_bins)),
                "model_selection_json": str(calibration.selection_json),
                "candidate_model_product": str(calibration.covariance_product),
            },
            "likelihood_width_covariance_calibrated": True,
            "likelihood_width_publication_final": bool(
                str(getattr(cfg, "RH3_EXPECTED_GRATING", "RH3")).strip().upper() == "RH3"
            ),
            "likelihood_width_publication_final_note": (
                "Covariance calibration passed. Publication-final RH3 interpretation additionally "
                "requires the production RH3 grating/data; RL or other surrogate integration runs "
                "remain development tests even when their covariance is internally calibrated."
            ),
            "fraction_definition": (
                "stellar-template light fraction of component A over RH3_FIT_REST_RANGE_ANGSTROM; "
                "each XSL SSP independently normalized to unit mean flux over that band; additive polynomial excluded"
            ),
            "template_product": str(template_path),
            "lsf_product": str(lsf_path),
            "lsf_empirical_observed_range_angstrom": [lsf_model.empirical_min, lsf_model.empirical_max],
            "checkpoints_deleted_on_success": bool(cfg.SCRIPT03_DELETE_CHECKPOINTS_ON_SUCCESS),
            "quality_flags": quality_flags,
            "products": {
                "likelihood_cubes_npz": str(cube_path),
                "log_spectra_and_local_best_fits_npz": str(selected_path),
                "local_likelihood_summary_ecsv": str(table_path),
                "prepared_xsl_templates_npz": str(template_path),
                "covariance_candidates_npz": str(calibration.covariance_product),
                "covariance_calibration_fits_npz": str(calibration.calibration_fit_product),
                "covariance_validation_grids_npz": str(calibration.validation_grid_product),
                "covariance_validation_bins_ecsv": str(run.products_dir / "covariance_validation_bins.ecsv"),
                "covariance_model_selection_json": str(calibration.selection_json),
                "covariance_iteration_history_ecsv": str(run.products_dir / "RH3_covariance_iteration_history.ecsv"),
                "covariance_model_comparison_ecsv": str(run.products_dir / "RH3_covariance_model_comparison.ecsv"),
            },
        }
        io.write_json(manifest, run.metadata_dir / "script03_manifest.json")

        if bool(cfg.SCRIPT03_DELETE_CHECKPOINTS_ON_SUCCESS):
            shutil.rmtree(checkpoints)
            logger.info("Successful consolidation complete; deleted intermediary checkpoint directory: %s", checkpoints)
        else:
            logger.info("Successful consolidation complete; checkpoint retention requested by config: %s", checkpoints)

        _section(logger, "SCRIPT 3 COMPLETE", "=")
        logger.info("Likelihood cube product: %s", cube_path)
        logger.info("Selected-fit spectra product: %s", selected_path)
        logger.info("Per-bin summary: %s", table_path)
        logger.info("Quality flags: %s", quality_flags)
        logger.info("Elapsed time: %.2f min", (time.perf_counter() - start_time) / 60.0)
        logger.info(
            "Covariance calibration and representative M1--M4 full-grid validation passed; "
            "the saved Delta-chi2 widths use the frozen selected RH3 covariance model."
        )
        return 0

    except Exception:
        logger.exception(
            "Script 3 failed/interrupted. Completed per-bin checkpoints remain in %s and may be reused with --resume.",
            checkpoints,
        )
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
