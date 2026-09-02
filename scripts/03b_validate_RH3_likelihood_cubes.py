#!/usr/bin/env python3
"""CRD_DAP Script 03b: validate completed Script-3 RH3 likelihood cubes.

This is a CHEAP, READ-ONLY post-processing stage.  It never calls pPXF and it
never changes the Script-3 science products.  Its purpose is to extract enough
numerical evidence from an expensive Script-3 run that most diagnostics can be
changed/repeated without rebuilding the likelihood cubes.

Validation contract
-------------------
Script 03b separates four questions that must not be conflated:

1. Mechanical integrity
   * product shapes and grids agree;
   * selected-fit chi-square recomputes from the saved spectrum/model/noise;
   * the duplicated-template light-fraction bookkeeping is internally
     consistent when the saved XSL product is available.

2. Component-swap symmetry
   Before Script 4 assigns physical disk labels, the model must be invariant
   under

       (V_A, V_B, f_A) <-> (V_B, V_A, 1-f_A).

   The script measures this over the full 3-D cube and near the supported part
   of the cube.  Large deviations are a numerical/implementation warning.

3. Information/topology/boundaries
   * effective number of likelihood-supported states;
   * entropy/effective support and 90--99.9% mass cell counts;
   * velocity/fraction grid-edge likelihood mass;
   * sigma-boundary-state fractions;
   * strict and epsilon-thick local-minimum plateau counts.

   The epsilon thickness is derived from the computed A/B swap mismatch near
   the low-chi-square part of each cube.  It is a NUMERICAL topology scale, not
   a confidence threshold and not the final Script-5 watershed definition.

4. Spectral leverage / bad-pixel / residual-scale diagnosis
   For the saved one-component and local-best two-component models, compute

       r_j = (F_j - M_j) / sigma_j
       chi2_j = r_j**2

   and quantify how much total chi-square comes from the worst 1, 5, 10, and
   top-1-percent pixels.  Recurrent high-|r| wavelengths across many bins are
   saved explicitly in both rest-frame and (when redshift is available)
   observed-frame wavelength.

   The validator also distinguishes a *localized wavelength problem* from a
   *bin-wide residual-scale problem*.  For every bin it measures the formal
   noise distribution, robust residual scatter, median |r|, fraction of good
   pixels above the recurrent-residual threshold, galaxy/model excursions, and
   saved additive-polynomial amplitude.  A diagnostic ``GLOBAL_RESIDUAL_PATHOLOGY``
   flag identifies bins for which most fitted pixels are already above the
   chosen standardized-residual threshold; this is a development diagnostic,
   not an automatic masking/rescaling instruction.

Outputs are written below ``<script3-run>/validation/03b`` by default.  The
source Script-3 products are never modified.  In addition to validation, 03b can
propose a conservative fixed observed-frame mask from recurrent *raw negative*
outliers and quantify its impact on the already-saved selected models.  This
proposal remains diagnostic: 03b never edits the target config and never reruns
pPXF.
"""

from __future__ import annotations

import argparse
from collections import OrderedDict
from datetime import datetime, timezone
import json
import logging
import os
import re
from pathlib import Path
import sys
from typing import Iterable

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, SymLogNorm

from astropy.io import fits
from astropy.table import Table

try:
    from scipy import ndimage
except Exception as exc:  # pragma: no cover - environment diagnostic
    raise ImportError(
        "Script 03b requires scipy.ndimage for plateau-aware topology diagnostics."
    ) from exc


SCRIPT_NAME = "CRD_DAP.03b_validate_RH3_likelihood_cubes"
MASS_LEVELS = (0.90, 0.95, 0.99, 0.995, 0.999)
RESIDUAL_THRESHOLDS = (5.0, 10.0, 20.0, 50.0, 100.0)
TOPOLOGY_LOW_CHI_FRACTION = 0.25
SWAP_NEAR_MIN_DELTA_CHI2 = 25.0  # numerical diagnostic only; not a confidence cut
SWAP_REL_PASS = 1.0e-3
SWAP_REL_WARN = 1.0e-2
CHI2_RECOMPUTE_RTOL = 5.0e-4
CHI2_RECOMPUTE_ATOL = 1.0e-3
LEVERAGE_TOP1_WARN = 0.25
LEVERAGE_TOP10_WARN = 0.50
LEVERAGE_MAX_R_WARN = 20.0
GRID_EDGE_WEIGHT_WARN = 0.25
SIGMA_BOUNDARY_STATE_WARN = 0.50
GLOBAL_RESIDUAL_ABS_R_DEFAULT = 10.0
GLOBAL_RESIDUAL_PIXEL_FRACTION_DEFAULT = 0.80
ROBUST_RESIDUAL_TO_NOISE_WARN = 5.0
MODEL_MEDIAN_OFFSET_NOISE_WARN = 10.0
RAW_NEGATIVE_OUTLIER_SIGMA_DEFAULT = 20.0


def _parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Cheap validation/post-processing of a completed Script-3 RH3 likelihood run."
    )
    p.add_argument(
        "--script3-run",
        required=True,
        help="Completed Script-3 run directory containing products/RH3_likelihood_cubes.npz.",
    )
    p.add_argument(
        "--script2-run",
        default=None,
        help=(
            "Optional current Script-2 run path for spatial bin maps/source table. "
            "If omitted, Script 03b uses source_script2_run in the Script-3 manifest."
        ),
    )
    p.add_argument(
        "--output-dir",
        default=None,
        help="Output directory. Default: <script3-run>/validation/03b",
    )
    p.add_argument(
        "--topology-epsilon",
        default="auto",
        help=(
            "Plateau thickness in total chi-square units, or 'auto' (default). "
            "Auto derives one epsilon per bin from robust low-chi A/B swap mismatch."
        ),
    )
    p.add_argument(
        "--representative-count",
        type=int,
        default=12,
        help="Maximum number of automatically selected representative-bin diagnostic figures.",
    )
    p.add_argument(
        "--recurrent-residual-threshold",
        type=float,
        default=10.0,
        help="|standardized residual| threshold used for recurrent-wavelength counts (default 10).",
    )
    p.add_argument(
        "--redshift",
        type=float,
        default=None,
        help=(
            "Optional target redshift used to add observed-frame wavelength columns/plots. "
            "If omitted, Script 03b first looks for 'redshift' in the Script-3 manifest. "
            "Older Script-3 manifests may not contain it."
        ),
    )
    p.add_argument(
        "--config",
        default=None,
        help=(
            "Optional target configuration file. It is used only as a fallback source of "
            "REDSHIFT when an older Script-3 manifest does not store redshift. The config is "
            "not modified and no science paths are opened by Script 03b."
        ),
    )
    p.add_argument(
        "--candidate-mask-mode",
        choices=("off", "raw-negative", "extended"),
        default="off",
        help=(
            "Legacy residual-derived candidate-mask experiment. Default is OFF because production "
            "mask selection should be externally anchored with Script 03d. 'raw-negative' uses only "
            "model-independent recurrent negative outliers as a diagnostic candidate; 'extended' "
            "also includes strong recurrent 2C residual excesses. Neither mode is the production "
            "source of truth unless explicitly chosen for a controlled test."
        ),
    )
    p.add_argument(
        "--candidate-raw-negative-fraction",
        type=float,
        default=0.20,
        help=(
            "Minimum fraction of all bins with a model-independent raw negative outlier at a "
            "wavelength for it to enter the conservative candidate mask (default 0.20)."
        ),
    )
    p.add_argument(
        "--candidate-residual-excess-fraction",
        type=float,
        default=0.20,
        help=(
            "For --candidate-mask-mode extended, minimum localized 2C residual excess divided by "
            "the number of non-pathology bins (default 0.20)."
        ),
    )
    p.add_argument(
        "--candidate-mask-padding-pixels",
        type=int,
        default=1,
        help="Dilate each selected log-wavelength feature by this many pixels (default 1).",
    )
    p.add_argument(
        "--candidate-mask-bridge-pixels",
        type=int,
        default=1,
        help=(
            "Merge selected wavelength groups separated by at most this many unselected log pixels "
            "before padding (default 1)."
        ),
    )
    p.add_argument(
        "--candidate-mask-max-fit-fraction",
        type=float,
        default=0.15,
        help=(
            "Warn if the proposed fixed mask removes more than this fraction of the saved Script-3 "
            "log grid (default 0.15)."
        ),
    )
    p.add_argument(
        "--mask-impact-bins",
        default=None,
        help=(
            "Optional comma-separated bin IDs for additional before/after candidate-mask figures, "
            "for example '231,328'. No pPXF refitting is performed."
        ),
    )

    p.add_argument(
        "--global-residual-abs-r-threshold",
        type=float,
        default=GLOBAL_RESIDUAL_ABS_R_DEFAULT,
        help=(
            "Absolute standardized-residual threshold used for the bin-wide residual-pathology "
            "diagnostic (default 10). This does not mask or rescale data."
        ),
    )
    p.add_argument(
        "--global-residual-pixel-fraction",
        type=float,
        default=GLOBAL_RESIDUAL_PIXEL_FRACTION_DEFAULT,
        help=(
            "Flag a bin as GLOBAL_RESIDUAL_PATHOLOGY when at least this fraction of its valid "
            "pixels exceeds --global-residual-abs-r-threshold in the local-best 2C residual "
            "(default 0.80). Diagnostic only."
        ),
    )
    p.add_argument(
        "--raw-negative-outlier-sigma",
        type=float,
        default=RAW_NEGATIVE_OUTLIER_SIGMA_DEFAULT,
        help=(
            "Model-independent raw-spectrum negative-outlier threshold in robust-MAD sigma units "
            "below each bin's median flux (default 20). Diagnostic only."
        ),
    )
    p.add_argument(
        "--no-spatial-maps",
        action="store_true",
        help="Skip Script-2 bin-map lookup and spatial diagnostic maps.",
    )
    return p


def _setup_logger(output_dir: Path) -> logging.Logger:
    output_dir.mkdir(parents=True, exist_ok=True)
    logs = output_dir / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger(SCRIPT_NAME)
    logger.setLevel(logging.INFO)
    logger.propagate = False
    for h in list(logger.handlers):
        logger.removeHandler(h)
        try:
            h.close()
        except Exception:
            pass
    fmt = logging.Formatter(
        "%(asctime)s | %(levelname)-8s | CRD_DAP.03b_validate_RH3_likelihood_cubes | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(fmt)
    logger.addHandler(sh)
    fh = logging.FileHandler(logs / "03b_validate_RH3_likelihood_cubes.log", mode="w")
    fh.setFormatter(fmt)
    logger.addHandler(fh)
    return logger


def _section(logger: logging.Logger, title: str, char: str = "-") -> None:
    line = char * max(18, len(title))
    logger.info(line)
    logger.info(title)
    logger.info(line)


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


def _require(path: Path, label: str) -> Path:
    if not path.is_file():
        raise FileNotFoundError(f"Missing {label}: {path}")
    return path


def _safe_float(value, default=np.nan) -> float:
    try:
        out = float(value)
        return out if np.isfinite(out) else float(default)
    except Exception:
        return float(default)


def _morton2000_refractive_index(vacuum_angstrom: np.ndarray) -> np.ndarray:
    """Dry-air refractive index; same Morton (2000) expression as crd_utils.templates."""
    wave = np.asarray(vacuum_angstrom, dtype=float)
    if np.any(wave < 2000.0):
        raise ValueError("Morton2000 air/vacuum conversion requires wavelength >= 2000 A.")
    sigma2 = (1.0e4 / wave) ** 2
    return (
        1.0
        + 8.34254e-5
        + 2.406147e-2 / (130.0 - sigma2)
        + 1.5998e-4 / (38.9 - sigma2)
    )


def _vacuum_to_air(wavelength_angstrom: np.ndarray) -> np.ndarray:
    wave = np.asarray(wavelength_angstrom, dtype=float)
    return wave / _morton2000_refractive_index(wave)


def _air_to_vacuum(wavelength_angstrom: np.ndarray, *, maxiter: int = 12) -> np.ndarray:
    air = np.asarray(wavelength_angstrom, dtype=float)
    vac = air * 1.00028
    for _ in range(int(maxiter)):
        new = air * _morton2000_refractive_index(vac)
        if np.all(np.abs(new - vac) <= 1.0e-10 * np.maximum(vac, 1.0)):
            vac = new
            break
        vac = new
    return vac


def _convert_wavelength_medium(wavelength_angstrom: np.ndarray, from_medium: str, to_medium: str) -> np.ndarray:
    src = str(from_medium).strip().lower()
    dst = str(to_medium).strip().lower()
    if src not in {"air", "vacuum"} or dst not in {"air", "vacuum"}:
        raise ValueError("Wavelength medium must be 'air' or 'vacuum'.")
    wave = np.asarray(wavelength_angstrom, dtype=float)
    if src == dst:
        return wave.copy()
    if src == "vacuum" and dst == "air":
        return _vacuum_to_air(wave)
    return _air_to_vacuum(wave)


def _observed_wavelength_grids(rest_wave: np.ndarray, redshift: float, template_medium: str, logger) -> tuple[np.ndarray, np.ndarray]:
    """Return observed-air and observed-vacuum wavelength arrays from the saved rest grid."""
    rest = np.asarray(rest_wave, dtype=float)
    if not np.isfinite(redshift):
        nan = np.full_like(rest, np.nan, dtype=float)
        return nan, nan
    medium = str(template_medium).strip().lower()
    if medium not in {"air", "vacuum"}:
        logger.warning(
            "Unknown saved template wavelength medium %r; exact observed air/vacuum conversion unavailable.",
            template_medium,
        )
        nan = np.full_like(rest, np.nan, dtype=float)
        return nan, nan
    rest_vac = _convert_wavelength_medium(rest, medium, "vacuum")
    observed_vac = rest_vac * (1.0 + float(redshift))
    observed_air = _vacuum_to_air(observed_vac)
    return observed_air, observed_vac


def _mad_sigma(values: np.ndarray) -> float:
    """Robust Gaussian-equivalent scatter 1.4826*MAD, ignoring non-finite values."""
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return np.nan
    med = float(np.median(v))
    return float(1.4826 * np.median(np.abs(v - med)))


def _read_redshift_from_config(config_path: str | None) -> float:
    """Read REDSHIFT from a target config without importing CRD_DAP or validating paths."""
    if config_path in (None, ""):
        return np.nan
    import importlib.util

    path = Path(config_path).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"--config does not exist: {path}")
    spec = importlib.util.spec_from_file_location("crd_dap_03b_target_config", path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load target config: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if not hasattr(module, "REDSHIFT"):
        raise AttributeError(f"Target config does not define REDSHIFT: {path}")
    z = float(module.REDSHIFT)
    if not np.isfinite(z) or z < 0:
        raise ValueError(f"Target config REDSHIFT must be finite and non-negative; got {z!r}")
    return z


def _resolve_redshift(args_redshift, manifest: dict, logger: logging.Logger, config_path: str | None = None) -> float:
    """Resolve target redshift with explicit > manifest > config precedence."""
    if args_redshift is not None:
        z = float(args_redshift)
        if not np.isfinite(z) or z < 0:
            raise ValueError("--redshift must be finite and non-negative.")
        logger.info("Observed-frame wavelength conversion: using explicit --redshift %.8f", z)
        return z
    for key in ("redshift", "target_redshift", "z"):
        if key in manifest:
            z = _safe_float(manifest.get(key))
            if np.isfinite(z) and z >= 0:
                logger.info("Observed-frame wavelength conversion: using Script-3 manifest %s=%.8f", key, z)
                return z
    zcfg = _read_redshift_from_config(config_path)
    if np.isfinite(zcfg):
        logger.info("Observed-frame wavelength conversion: using target config REDSHIFT=%.8f", zcfg)
        return zcfg
    logger.warning(
        "No target redshift is stored in this Script-3 manifest and no usable --config/--redshift "
        "was supplied. Observed-frame wavelength columns will be NaN."
    )
    return np.nan


def _log_center_edges(wave: np.ndarray) -> np.ndarray:
    """Return bin edges for a strictly increasing, approximately log-spaced wavelength grid."""
    w = np.asarray(wave, dtype=float)
    if w.ndim != 1 or w.size < 2 or np.any(~np.isfinite(w)) or np.any(np.diff(w) <= 0):
        raise ValueError("Saved wavelength grid must be finite and strictly increasing.")
    edges = np.empty(w.size + 1, dtype=float)
    edges[1:-1] = np.sqrt(w[:-1] * w[1:])
    edges[0] = w[0] * w[0] / edges[1]
    edges[-1] = w[-1] * w[-1] / edges[-2]
    return edges


def _bridge_small_false_gaps(mask: np.ndarray, max_gap: int) -> np.ndarray:
    """Fill False runs of length <= max_gap that are bracketed by selected pixels."""
    out = np.asarray(mask, dtype=bool).copy()
    max_gap = int(max_gap)
    if max_gap <= 0 or out.size == 0:
        return out
    true_idx = np.flatnonzero(out)
    if true_idx.size < 2:
        return out
    for left, right in zip(true_idx[:-1], true_idx[1:]):
        gap = int(right - left - 1)
        if 0 < gap <= max_gap:
            out[left + 1:right] = True
    return out


def _contiguous_true_groups(mask: np.ndarray) -> list[tuple[int, int]]:
    idx = np.flatnonzero(np.asarray(mask, dtype=bool))
    if idx.size == 0:
        return []
    groups: list[tuple[int, int]] = []
    start = prev = int(idx[0])
    for value in idx[1:]:
        value = int(value)
        if value == prev + 1:
            prev = value
            continue
        groups.append((start, prev))
        start = prev = value
    groups.append((start, prev))
    return groups


def _parse_bin_id_list(text: str | None, nbin: int) -> list[int]:
    if text in (None, ""):
        return []
    out: list[int] = []
    for token in str(text).split(","):
        token = token.strip()
        if not token:
            continue
        bid = int(token)
        if bid < 0 or bid >= nbin:
            raise ValueError(f"--mask-impact-bins contains out-of-range bin {bid}; valid 0..{nbin-1}")
        if bid not in out:
            out.append(bid)
    return out


def _convert_rest_edges_to_observed(
    rest_edges: np.ndarray, redshift: float, template_medium: str, science_medium: str
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Convert rest-grid edges to observed air, vacuum, and science-medium wavelengths."""
    rest = np.asarray(rest_edges, dtype=float)
    if not np.isfinite(redshift):
        nan = np.full_like(rest, np.nan)
        return nan, nan, nan
    tmed = str(template_medium).strip().lower()
    smed = str(science_medium).strip().lower()
    if tmed not in {"air", "vacuum"} or smed not in {"air", "vacuum"}:
        nan = np.full_like(rest, np.nan)
        return nan, nan, nan
    rest_vac = _convert_wavelength_medium(rest, tmed, "vacuum")
    obs_vac = rest_vac * (1.0 + float(redshift))
    obs_air = _vacuum_to_air(obs_vac)
    obs_science = obs_vac if smed == "vacuum" else obs_air
    return obs_air, obs_vac, obs_science


def _candidate_mask_products(
    *,
    wave: np.ndarray,
    raw_negative_counts: np.ndarray,
    rank_score: np.ndarray,
    nbin: int,
    n_nonpathology: int,
    redshift: float,
    template_medium: str,
    science_medium: str,
    mode: str,
    raw_fraction: float,
    residual_excess_fraction: float,
    padding_pixels: int,
    bridge_pixels: int,
) -> tuple[np.ndarray, Table, Table]:
    """Build conservative candidate fixed-wavelength masks and review-only intervals."""
    nlog = len(wave)
    raw_fraction = float(raw_fraction)
    residual_excess_fraction = float(residual_excess_fraction)
    if not (0 < raw_fraction <= 1):
        raise ValueError("--candidate-raw-negative-fraction must lie in (0, 1].")
    if not (0 < residual_excess_fraction <= 1):
        raise ValueError("--candidate-residual-excess-fraction must lie in (0, 1].")
    if padding_pixels < 0 or bridge_pixels < 0:
        raise ValueError("candidate mask padding/bridge pixels must be >= 0.")

    raw_base = np.asarray(raw_negative_counts, dtype=float) >= raw_fraction * float(max(1, nbin))
    residual_base = np.asarray(rank_score, dtype=float) >= residual_excess_fraction * float(max(1, n_nonpathology))

    selected_base = np.zeros(nlog, dtype=bool)
    if mode == "raw-negative":
        selected_base = raw_base.copy()
    elif mode == "extended":
        selected_base = raw_base | residual_base
    elif mode != "off":
        raise ValueError(f"Unknown candidate mask mode: {mode}")

    selected = _bridge_small_false_gaps(selected_base, bridge_pixels)
    if padding_pixels > 0 and np.any(selected):
        selected = ndimage.binary_dilation(selected, iterations=int(padding_pixels))

    edges = _log_center_edges(wave)
    obs_air_edges, obs_vac_edges, obs_science_edges = _convert_rest_edges_to_observed(
        edges, redshift, template_medium, science_medium
    )

    def build_table(mask: np.ndarray, included: bool) -> Table:
        rows = []
        for interval_id, (i0, i1) in enumerate(_contiguous_true_groups(mask)):
            sl = slice(i0, i1 + 1)
            raw_hit = bool(np.any(raw_base[sl]))
            resid_hit = bool(np.any(residual_base[sl]))
            if raw_hit and resid_hit:
                reason = "RAW_NEGATIVE+RECURRENT_RESIDUAL"
            elif raw_hit:
                reason = "RAW_NEGATIVE"
            else:
                reason = "RECURRENT_RESIDUAL_REVIEW"
            rows.append((
                interval_id, included, reason, i0, i1, i1 - i0 + 1,
                float(edges[i0]), float(edges[i1 + 1]),
                float(obs_air_edges[i0]) if np.isfinite(obs_air_edges[i0]) else np.nan,
                float(obs_air_edges[i1 + 1]) if np.isfinite(obs_air_edges[i1 + 1]) else np.nan,
                float(obs_vac_edges[i0]) if np.isfinite(obs_vac_edges[i0]) else np.nan,
                float(obs_vac_edges[i1 + 1]) if np.isfinite(obs_vac_edges[i1 + 1]) else np.nan,
                float(obs_science_edges[i0]) if np.isfinite(obs_science_edges[i0]) else np.nan,
                float(obs_science_edges[i1 + 1]) if np.isfinite(obs_science_edges[i1 + 1]) else np.nan,
                int(np.max(np.asarray(raw_negative_counts)[sl])) if i1 >= i0 else 0,
                float(np.max(np.asarray(raw_negative_counts, dtype=float)[sl]) / float(max(1, nbin))),
                float(np.max(np.asarray(rank_score, dtype=float)[sl])) if i1 >= i0 else np.nan,
            ))
        names = (
            "INTERVAL_ID", "INCLUDED_IN_CANDIDATE", "REASON", "LOG_INDEX_LO", "LOG_INDEX_HI", "N_LOG_PIXELS",
            "REST_LO_ANGSTROM", "REST_HI_ANGSTROM",
            "OBSERVED_AIR_LO_ANGSTROM", "OBSERVED_AIR_HI_ANGSTROM",
            "OBSERVED_VACUUM_LO_ANGSTROM", "OBSERVED_VACUUM_HI_ANGSTROM",
            "OBSERVED_SCIENCE_LO_ANGSTROM", "OBSERVED_SCIENCE_HI_ANGSTROM",
            "PEAK_RAW_NEGATIVE_BIN_COUNT", "PEAK_RAW_NEGATIVE_BIN_FRACTION", "PEAK_LOCALIZED_RESIDUAL_EXCESS",
        )
        return Table(rows=rows, names=names)

    candidate = build_table(selected, True)
    review_mask = _bridge_small_false_gaps(residual_base & ~raw_base, bridge_pixels)
    if padding_pixels > 0 and np.any(review_mask):
        review_mask = ndimage.binary_dilation(review_mask, iterations=int(padding_pixels))
    review = build_table(review_mask, False)

    for table in (candidate, review):
        table.meta["CANDIDATE_MASK_MODE"] = mode
        table.meta["RAW_NEGATIVE_BIN_FRACTION_THRESHOLD"] = raw_fraction
        table.meta["RESIDUAL_EXCESS_NONPATHOLOGY_FRACTION_THRESHOLD"] = residual_excess_fraction
        table.meta["PADDING_LOG_PIXELS"] = int(padding_pixels)
        table.meta["BRIDGE_LOG_PIXELS"] = int(bridge_pixels)
        table.meta["REDSHIFT"] = None if not np.isfinite(redshift) else float(redshift)
        table.meta["REST_WAVELENGTH_MEDIUM"] = str(template_medium)
        table.meta["OBSERVED_SCIENCE_WAVELENGTH_MEDIUM"] = str(science_medium)
    return np.asarray(selected, dtype=bool), candidate, review


def _candidate_mask_impact_table(
    spectra: dict, candidate_mask: np.ndarray, global_r_thr: float, global_r_frac: float, raw_negative_sigma: float
) -> Table:
    """Recompute selected-model diagnostics after masking, without refitting pPXF."""
    galaxy_all = np.asarray(spectra["galaxy"], dtype=float)
    noise_all = np.asarray(spectra["noise"], dtype=float)
    good_all = np.asarray(spectra["good"], dtype=bool)
    one_all = np.asarray(spectra["one_component_model"], dtype=float)
    two_all = np.asarray(spectra["local_best_two_component_model"], dtype=float)
    rows = []
    for bid in range(galaxy_all.shape[0]):
        galaxy = galaxy_all[bid]
        noise = noise_all[bid]
        good0 = good_all[bid]
        good1 = good0 & ~candidate_mask
        one0, _, _ = _residual_metrics(galaxy, one_all[bid], noise, good0)
        one1, _, _ = _residual_metrics(galaxy, one_all[bid], noise, good1)
        two0, r20, _ = _residual_metrics(galaxy, two_all[bid], noise, good0)
        two1, r21, _ = _residual_metrics(galaxy, two_all[bid], noise, good1)
        d0 = _noise_residual_diagnostics(galaxy, two_all[bid], noise, good0)
        d1 = _noise_residual_diagnostics(galaxy, two_all[bid], noise, good1)
        frac0 = float(np.mean(np.abs(r20[np.isfinite(r20)]) > global_r_thr)) if np.any(np.isfinite(r20)) else np.nan
        frac1 = float(np.mean(np.abs(r21[np.isfinite(r21)]) > global_r_thr)) if np.any(np.isfinite(r21)) else np.nan

        def n_raw(good):
            valid = good & np.isfinite(galaxy)
            vals = galaxy[valid]
            if vals.size < 3:
                return 0
            med = float(np.median(vals))
            scatter = _mad_sigma(vals)
            if not np.isfinite(scatter) or scatter <= 0:
                return 0
            return int(np.count_nonzero((vals - med) / scatter < -raw_negative_sigma))

        rows.append((
            bid, int(np.count_nonzero(good0)), int(np.count_nonzero(good1)), int(np.count_nonzero(good0 & candidate_mask)),
            one0["chi2_recomputed"], one1["chi2_recomputed"],
            two0["chi2_recomputed"], two1["chi2_recomputed"],
            two0["top1_frac"], two1["top1_frac"],
            two0["top10_frac"], two1["top10_frac"],
            two0["max_abs_r"], two1["max_abs_r"],
            d0["robust_residual_to_noise"], d1["robust_residual_to_noise"],
            d0["median_abs_r"], d1["median_abs_r"],
            frac0, frac1, bool(np.isfinite(frac0) and frac0 >= global_r_frac), bool(np.isfinite(frac1) and frac1 >= global_r_frac),
            n_raw(good0), n_raw(good1),
        ))
    return Table(rows=rows, names=(
        "BIN_ID", "N_GOOD_BEFORE", "N_GOOD_AFTER", "N_MASKED_CANDIDATE",
        "ONE_CHI2_SAME_MODEL_BEFORE", "ONE_CHI2_SAME_MODEL_AFTER",
        "TWO_CHI2_SAME_MODEL_BEFORE", "TWO_CHI2_SAME_MODEL_AFTER",
        "TWO_TOP1_CHI2_FRAC_BEFORE", "TWO_TOP1_CHI2_FRAC_AFTER",
        "TWO_TOP10_CHI2_FRAC_BEFORE", "TWO_TOP10_CHI2_FRAC_AFTER",
        "TWO_MAX_ABS_R_BEFORE", "TWO_MAX_ABS_R_AFTER",
        "TWO_ROBUST_RESIDUAL_TO_NOISE_BEFORE", "TWO_ROBUST_RESIDUAL_TO_NOISE_AFTER",
        "TWO_MEDIAN_ABS_R_BEFORE", "TWO_MEDIAN_ABS_R_AFTER",
        "TWO_GLOBAL_ABS_R_FRACTION_BEFORE", "TWO_GLOBAL_ABS_R_FRACTION_AFTER",
        "GLOBAL_RESIDUAL_PATHOLOGY_BEFORE", "GLOBAL_RESIDUAL_PATHOLOGY_AFTER",
        "GALAXY_N_RAW_NEGATIVE_OUTLIERS_BEFORE", "GALAXY_N_RAW_NEGATIVE_OUTLIERS_AFTER",
    ))


def _write_candidate_config_snippet(path: Path, table: Table, science_medium: str) -> None:
    lines = [
        "# Candidate observed-frame fixed mask proposed by Script 03b.",
        "# IMPORTANT: values are in the native SCIENCE wavelength medium used by Script 3",
        f"# (science medium for this run: {science_medium}).",
        "# Validate with 03c before copying into the production target config.",
        "RH3_MASK_OBSERVED_RANGES_ANGSTROM = [",
    ]
    for row in table:
        lo = float(row["OBSERVED_SCIENCE_LO_ANGSTROM"])
        hi = float(row["OBSERVED_SCIENCE_HI_ANGSTROM"])
        reason = str(row["REASON"])
        if np.isfinite(lo) and np.isfinite(hi):
            lines.append(f"    ({lo:.2f}, {hi:.2f}),  # {reason}")
    lines.append("]")
    lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")


def _plot_candidate_mask_impact(impact: Table, output_path: Path) -> None:
    before = np.asarray(impact["TWO_TOP1_CHI2_FRAC_BEFORE"], dtype=float)
    after = np.asarray(impact["TWO_TOP1_CHI2_FRAC_AFTER"], dtype=float)
    rb = np.asarray(impact["TWO_ROBUST_RESIDUAL_TO_NOISE_BEFORE"], dtype=float)
    ra = np.asarray(impact["TWO_ROBUST_RESIDUAL_TO_NOISE_AFTER"], dtype=float)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    use = np.isfinite(before) & np.isfinite(after)
    axes[0].scatter(before[use], after[use], s=14, alpha=0.6)
    if np.any(use):
        lo = min(np.min(before[use]), np.min(after[use]))
        hi = max(np.max(before[use]), np.max(after[use]))
        axes[0].plot([lo, hi], [lo, hi], ls="--", lw=0.8)
    axes[0].set_xlabel("Worst-pixel chi2 fraction before")
    axes[0].set_ylabel("Worst-pixel chi2 fraction after")
    axes[0].set_title("Candidate mask: spectral leverage")
    use = np.isfinite(rb) & np.isfinite(ra) & (rb > 0) & (ra > 0)
    axes[1].scatter(rb[use], ra[use], s=14, alpha=0.6)
    if np.any(use):
        lo = min(np.min(rb[use]), np.min(ra[use]))
        hi = max(np.max(rb[use]), np.max(ra[use]))
        axes[1].plot([lo, hi], [lo, hi], ls="--", lw=0.8)
        axes[1].set_xscale("log")
        axes[1].set_yscale("log")
    axes[1].set_xlabel("Robust residual/noise before")
    axes[1].set_ylabel("Robust residual/noise after")
    axes[1].set_title("Same saved 2C model; no refit")
    fig.tight_layout()
    fig.savefig(output_path, dpi=160)
    plt.close(fig)


def _plot_mask_impact_bin(
    bid: int, spectra: dict, wave: np.ndarray, candidate_mask: np.ndarray, path: Path
) -> None:
    galaxy = np.asarray(spectra["galaxy"][bid], dtype=float)
    noise = np.asarray(spectra["noise"][bid], dtype=float)
    good = np.asarray(spectra["good"][bid], dtype=bool)
    two = np.asarray(spectra["local_best_two_component_model"][bid], dtype=float)
    after = good & ~candidate_mask
    fig, axes = plt.subplots(2, 1, figsize=(13, 8), sharex=True)
    axes[0].plot(wave, galaxy, lw=0.7, label="saved RH3 spectrum")
    axes[0].plot(wave, two, lw=1.0, label="saved local-best 2C model")
    if np.any(good & candidate_mask):
        axes[0].scatter(wave[good & candidate_mask], galaxy[good & candidate_mask], s=18, marker="x", label="candidate-masked")
    axes[0].legend(fontsize=8)
    axes[0].set_ylabel("Normalized flux")
    axes[0].set_title(f"Bin {bid}: proposed fixed mask (saved model; no refit)")
    r = np.full_like(galaxy, np.nan)
    valid = good & np.isfinite(galaxy) & np.isfinite(two) & np.isfinite(noise) & (noise > 0)
    r[valid] = (galaxy[valid] - two[valid]) / noise[valid]
    axes[1].plot(wave[after], r[after], lw=0.7, label="retained residual")
    if np.any(good & candidate_mask):
        axes[1].scatter(wave[good & candidate_mask], r[good & candidate_mask], s=18, marker="x", label="removed residual")
    axes[1].axhline(0.0, lw=0.8)
    axes[1].set_xlabel("Rest-frame wavelength (Angstrom)")
    axes[1].set_ylabel("(data-model)/noise")
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def _noise_residual_diagnostics(galaxy, model, noise, good) -> dict:
    """Summarize formal noise, residual scale, model excursions, and standardized residuals."""
    galaxy = np.asarray(galaxy, dtype=float)
    model = np.asarray(model, dtype=float)
    noise = np.asarray(noise, dtype=float)
    good = np.asarray(good, dtype=bool)
    base = good & np.isfinite(galaxy) & np.isfinite(model) & np.isfinite(noise)
    valid = base & (noise > 0)
    out = {
        "n_good": int(np.count_nonzero(good)),
        "n_valid": int(np.count_nonzero(valid)),
        "noise_nonpositive_frac": np.nan,
        "noise_min": np.nan, "noise_p01": np.nan, "noise_median": np.nan,
        "noise_p99": np.nan, "noise_max": np.nan,
        "residual_median": np.nan, "residual_mad_sigma": np.nan,
        "median_abs_residual": np.nan, "median_abs_r": np.nan,
        "robust_residual_to_noise": np.nan, "model_median_offset_noise": np.nan,
        "model_min": np.nan, "model_p01": np.nan, "model_median": np.nan,
        "model_p99": np.nan, "model_max": np.nan,
        "galaxy_negative_frac": np.nan,
    }
    if np.any(base):
        out["noise_nonpositive_frac"] = float(np.mean(noise[base] <= 0))
    if not np.any(valid):
        return out
    n = noise[valid]
    g = galaxy[valid]
    m = model[valid]
    residual = g - m
    r = residual / n
    out.update({
        "noise_min": float(np.min(n)),
        "noise_p01": float(np.percentile(n, 1.0)),
        "noise_median": float(np.median(n)),
        "noise_p99": float(np.percentile(n, 99.0)),
        "noise_max": float(np.max(n)),
        "residual_median": float(np.median(residual)),
        "residual_mad_sigma": _mad_sigma(residual),
        "median_abs_residual": float(np.median(np.abs(residual))),
        "median_abs_r": float(np.median(np.abs(r))),
        "model_min": float(np.min(m)),
        "model_p01": float(np.percentile(m, 1.0)),
        "model_median": float(np.median(m)),
        "model_p99": float(np.percentile(m, 99.0)),
        "model_max": float(np.max(m)),
        "galaxy_negative_frac": float(np.mean(g < 0)),
    })
    if out["noise_median"] > 0:
        out["robust_residual_to_noise"] = out["residual_mad_sigma"] / out["noise_median"]
        out["model_median_offset_noise"] = abs(out["residual_median"]) / out["noise_median"]
    return out


def _path_is_foreign_platform(value: str) -> bool:
    """Detect manifest paths written on another OS before Path.resolve mangles them."""
    s = str(value)
    windows_abs = bool(re.match(r"^[A-Za-z]:[\\/]", s)) or s.startswith("\\\\")
    posix_abs = s.startswith("/")
    if os.name == "posix" and windows_abs:
        return True
    if os.name == "nt" and posix_abs:
        return True
    return False


def _load_products(script3_run: Path):
    products = script3_run / "products"
    cube_path = _require(products / "RH3_likelihood_cubes.npz", "Script-3 likelihood cube")
    spectra_path = _require(
        products / "RH3_log_spectra_and_local_best_fits.npz",
        "Script-3 selected-fit spectra",
    )
    manifest_path = _require(script3_run / "metadata" / "script03_manifest.json", "Script-3 manifest")
    template_path = products / "XSL_RH3_templates.npz"
    summary_path = products / "RH3_local_likelihood_summary.ecsv"

    with np.load(cube_path, allow_pickle=False) as d:
        cube = {k: np.asarray(d[k]) for k in d.files}
    with np.load(spectra_path, allow_pickle=False) as d:
        spectra = {k: np.asarray(d[k]) for k in d.files}
    templates = None
    if template_path.is_file():
        with np.load(template_path, allow_pickle=False) as d:
            templates = {k: np.asarray(d[k]) for k in d.files}
    manifest = _load_json(manifest_path)
    summary = Table.read(summary_path, format="ascii.ecsv") if summary_path.is_file() else None
    return cube_path, spectra_path, template_path, cube, spectra, templates, manifest, summary


def _validate_shapes(cube: dict, spectra: dict, manifest: dict) -> tuple[int, int, int, int, int]:
    required_cube = (
        "chi2_total", "reduced_chi2", "sigma_A", "sigma_B", "fit_status",
        "sigma_boundary", "VA_grid", "VB_grid", "fA_grid", "one_chi2_total",
        "local_best_chi2_total", "n_failed_states", "n_sigma_boundary_states",
    )
    required_spec = (
        "wavelength", "galaxy", "noise", "good", "one_component_model",
        "local_best_two_component_model", "local_best_two_component_weights",
    )
    missing = [k for k in required_cube if k not in cube] + [k for k in required_spec if k not in spectra]
    if missing:
        raise KeyError("Script-3 products are missing required arrays: " + ", ".join(missing))

    chi = cube["chi2_total"]
    if chi.ndim != 4:
        raise ValueError(f"chi2_total must be 4-D (bin,VA,VB,fA); got {chi.shape}")
    nbin, nva, nvb, nfa = chi.shape
    if nva != cube["VA_grid"].size or nvb != cube["VB_grid"].size or nfa != cube["fA_grid"].size:
        raise ValueError("Likelihood cube axes disagree with saved grid arrays.")
    for key in ("reduced_chi2", "sigma_A", "sigma_B", "fit_status", "sigma_boundary"):
        if cube[key].shape != chi.shape:
            raise ValueError(f"{key} shape {cube[key].shape} disagrees with chi2_total {chi.shape}")
    wave = spectra["wavelength"]
    nlog = wave.size
    for key in ("galaxy", "noise", "good", "one_component_model", "local_best_two_component_model"):
        if spectra[key].shape != (nbin, nlog):
            raise ValueError(f"{key} shape {spectra[key].shape} != expected {(nbin, nlog)}")
    if int(manifest.get("n_bins", nbin)) != nbin:
        raise ValueError("Script-3 manifest n_bins disagrees with likelihood cube.")
    return nbin, nva, nvb, nfa, nlog


def _build_swap_index(values: np.ndarray, target_values: np.ndarray, tol: float = 1e-10) -> np.ndarray:
    """Map each value onto an equal target value; fail instead of silently approximating."""
    idx = []
    for val in np.asarray(values, dtype=float):
        hits = np.where(np.isclose(np.asarray(target_values, dtype=float), val, rtol=0.0, atol=tol))[0]
        if hits.size != 1:
            raise ValueError(f"Could not uniquely map swap grid coordinate {val}")
        idx.append(int(hits[0]))
    return np.asarray(idx, dtype=int)


def _component_swap_cube(arr: np.ndarray, va: np.ndarray, vb: np.ndarray, fa: np.ndarray) -> np.ndarray:
    """Return arr evaluated at (V_B,V_A,1-f_A), aligned to original coordinates."""
    if arr.shape[-3:] != (va.size, vb.size, fa.size):
        raise ValueError("Array shape is incompatible with supplied swap grids.")
    # For every original VA value find the matching coordinate on the VB axis,
    # and vice versa. This keeps the implementation correct even if the two
    # arrays are stored as distinct objects.
    va_to_vb = _build_swap_index(va, vb)
    vb_to_va = _build_swap_index(vb, va)
    fa_comp = _build_swap_index(1.0 - fa, fa, tol=1e-8)
    # Desired swapped cell at [ia,ib,j] is original [vb_to_va[ib], va_to_vb[ia], fa_comp[j]].
    out = np.empty_like(arr)
    for ia in range(va.size):
        for ib in range(vb.size):
            out[..., ia, ib, :] = arr[..., vb_to_va[ib], va_to_vb[ia], fa_comp]
    return out


def _finite_delta(chi: np.ndarray) -> tuple[np.ndarray, float]:
    finite = np.isfinite(chi)
    if not np.any(finite):
        return np.full_like(chi, np.inf, dtype=float), np.nan
    cmin = float(np.nanmin(chi[finite]))
    d = np.full_like(chi, np.inf, dtype=float)
    d[finite] = np.asarray(chi[finite], dtype=float) - cmin
    d[d < 0] = 0.0
    return d, cmin


def _likelihood_weights(delta: np.ndarray) -> np.ndarray:
    w = np.zeros(delta.shape, dtype=np.float64)
    finite = np.isfinite(delta)
    if not np.any(finite):
        return w
    exponent = -0.5 * np.clip(np.asarray(delta[finite], dtype=np.float64), 0.0, 1490.0)
    vals = np.exp(exponent)
    s = float(np.sum(vals))
    if s > 0 and np.isfinite(s):
        w[finite] = vals / s
    return w


def _mass_cell_counts(weights: np.ndarray, levels: Iterable[float] = MASS_LEVELS) -> list[int]:
    vals = np.asarray(weights[np.isfinite(weights) & (weights > 0)], dtype=float)
    if vals.size == 0:
        return [0 for _ in levels]
    vals = np.sort(vals)[::-1]
    c = np.cumsum(vals)
    return [int(np.searchsorted(c, level, side="left") + 1) for level in levels]


def _entropy_metrics(weights: np.ndarray) -> tuple[float, float, float]:
    p = np.asarray(weights[weights > 0], dtype=float)
    if p.size == 0:
        return np.nan, np.nan, np.nan
    neff_ipr = 1.0 / float(np.sum(p * p))
    entropy = -float(np.sum(p * np.log(p)))
    neff_entropy = float(np.exp(entropy))
    return neff_ipr, entropy, neff_entropy


def _percentile_or_nan(values: np.ndarray, q: float) -> float:
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    return float(np.percentile(v, q)) if v.size else np.nan


def _robust_swap_epsilon(delta: np.ndarray, swap_abs: np.ndarray) -> float:
    """Estimate a numerical topology thickness from low-chi swap mismatch.

    Use the lowest-rank quartile of finite likelihood cells rather than a
    statistical Delta-chi-square confidence cut. This asks: how asymmetric is
    the computation in the portion of the cube closest to its minima?
    """
    finite = np.isfinite(delta) & np.isfinite(swap_abs)
    if not np.any(finite):
        return 0.0
    dvals = delta[finite]
    cutoff = np.quantile(dvals, TOPOLOGY_LOW_CHI_FRACTION)
    vals = swap_abs[finite & (delta <= cutoff)]
    if vals.size == 0:
        vals = swap_abs[finite]
    med = float(np.median(vals))
    mad = float(np.median(np.abs(vals - med)))
    robust = med + 3.0 * 1.4826 * mad
    p95 = float(np.percentile(vals, 95.0))
    eps = max(1.0e-8, min(max(robust, med), p95 if p95 > 0 else max(robust, med)))
    return float(eps)


def _topology_metrics(chi: np.ndarray, epsilon: float) -> tuple[int, int, int, float]:
    """Count strict minima and epsilon-thick minimum plateaus with 26-connectivity."""
    arr = np.asarray(chi, dtype=float)
    finite = np.isfinite(arr)
    if not np.any(finite):
        return 0, 0, 0, np.nan
    work = np.where(finite, arr, np.inf)
    footprint = np.ones((3, 3, 3), dtype=bool)
    neigh_min = ndimage.minimum_filter(work, footprint=footprint, mode="constant", cval=np.inf)
    strict = finite & np.isclose(work, neigh_min, rtol=0.0, atol=1e-12)
    thick = finite & (work <= neigh_min + float(max(0.0, epsilon)))
    structure = np.ones((3, 3, 3), dtype=int)
    _, n_strict = ndimage.label(strict, structure=structure)
    labels, n_thick = ndimage.label(thick, structure=structure)
    if n_thick:
        sizes = np.bincount(labels.ravel())[1:]
        largest = int(np.max(sizes)) if sizes.size else 0
    else:
        largest = 0
    return int(n_strict), int(n_thick), largest, largest / float(np.count_nonzero(finite))


def _fraction_on_grid_edges(weights: np.ndarray) -> tuple[float, float, float]:
    vel_edge = np.zeros(weights.shape, dtype=bool)
    vel_edge[0, :, :] = True
    vel_edge[-1, :, :] = True
    vel_edge[:, 0, :] = True
    vel_edge[:, -1, :] = True
    fa_edge = np.zeros(weights.shape, dtype=bool)
    fa_edge[:, :, 0] = True
    fa_edge[:, :, -1] = True
    return (
        float(np.sum(weights[vel_edge])),
        float(np.sum(weights[fa_edge])),
        float(np.sum(weights[vel_edge | fa_edge])),
    )


def _residual_metrics(galaxy, model, noise, good) -> tuple[dict, np.ndarray, np.ndarray]:
    galaxy = np.asarray(galaxy, dtype=float)
    model = np.asarray(model, dtype=float)
    noise = np.asarray(noise, dtype=float)
    good = np.asarray(good, dtype=bool)
    valid = good & np.isfinite(galaxy) & np.isfinite(model) & np.isfinite(noise) & (noise > 0)
    r = np.full(galaxy.shape, np.nan, dtype=np.float64)
    chi_pix = np.full(galaxy.shape, np.nan, dtype=np.float64)
    if not np.any(valid):
        metrics = {k: np.nan for k in (
            "chi2_recomputed", "max_abs_r", "top1_frac", "top5_frac", "top10_frac", "top1pct_frac"
        )}
        for t in RESIDUAL_THRESHOLDS:
            metrics[f"n_abs_r_gt_{int(t)}"] = 0
        return metrics, r, chi_pix
    r[valid] = (galaxy[valid] - model[valid]) / noise[valid]
    chi_pix[valid] = r[valid] ** 2
    vals = chi_pix[valid]
    order = np.sort(vals)[::-1]
    total = float(np.sum(order))

    def frac_top(n: int) -> float:
        if total <= 0 or not np.isfinite(total):
            return np.nan
        return float(np.sum(order[: min(n, order.size)]) / total)

    n1pct = max(1, int(np.ceil(0.01 * order.size)))
    metrics = {
        "chi2_recomputed": total,
        "max_abs_r": float(np.nanmax(np.abs(r[valid]))),
        "top1_frac": frac_top(1),
        "top5_frac": frac_top(5),
        "top10_frac": frac_top(10),
        "top1pct_frac": frac_top(n1pct),
    }
    for t in RESIDUAL_THRESHOLDS:
        metrics[f"n_abs_r_gt_{int(t)}"] = int(np.count_nonzero(np.abs(r[valid]) > t))
    return metrics, r, chi_pix


def _find_bin_map_and_source_table(manifest: dict, explicit_script2: str | None, nbin: int, logger):
    run_value = explicit_script2 if explicit_script2 is not None else manifest.get("source_script2_run")
    if run_value in (None, ""):
        logger.warning("No Script-2 path available; spatial maps and source-table augmentation will be skipped.")
        return None, None, None
    if explicit_script2 is None and _path_is_foreign_platform(str(run_value)):
        logger.warning(
            "Script-3 manifest records a Script-2 path from a different operating system (%s). "
            "It will not be resolved on this machine; pass --script2-run <current-path> to enable spatial maps.",
            run_value,
        )
        return None, None, None
    run = Path(run_value).expanduser().resolve()
    bin_path = run / "products" / "master_bins.fits"
    table_path = run / "products" / "master_bin_table.ecsv"
    if not bin_path.is_file():
        logger.warning("Script-2 master bin map unavailable at %s; spatial maps will be skipped.", bin_path)
        return None, None, run
    with fits.open(bin_path, memmap=False) as hdul:
        bin_map = np.asarray(hdul[0].data, dtype=int)
    source_table = Table.read(table_path, format="ascii.ecsv") if table_path.is_file() else None
    if source_table is not None and len(source_table) != nbin:
        logger.warning("Script-2 source table has %d rows but Script 3 has %d bins; table augmentation skipped.", len(source_table), nbin)
        source_table = None
    return bin_map, source_table, run


def _plot_bin_map(bin_map: np.ndarray, values: np.ndarray, path: Path, title: str, label: str, *, norm=None):
    image = np.full(bin_map.shape, np.nan, dtype=float)
    for bid, value in enumerate(np.asarray(values, dtype=float)):
        image[bin_map == bid] = value
    fig, ax = plt.subplots(figsize=(8, 6))
    im = ax.imshow(image, origin="lower", interpolation="nearest", norm=norm)
    ax.set_title(title)
    ax.set_xlabel("BL spatial x pixel")
    ax.set_ylabel("BL spatial y pixel")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(label)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def _plot_hist(values, path: Path, title: str, xlabel: str, logx: bool = False):
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    fig, ax = plt.subplots(figsize=(7, 5))
    if v.size:
        if logx:
            v = v[v > 0]
            if v.size:
                lo, hi = np.nanmin(v), np.nanmax(v)
                if hi > lo:
                    bins = np.geomspace(lo, hi, 40)
                    ax.hist(v, bins=bins)
                    ax.set_xscale("log")
                else:
                    ax.hist(v, bins=20)
            else:
                ax.text(0.5, 0.5, "No positive finite values", transform=ax.transAxes, ha="center")
        else:
            ax.hist(v, bins=40)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Number of bins")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def _robust_ylim(galaxy: np.ndarray, models: list[np.ndarray], good: np.ndarray) -> tuple[float, float]:
    vals = [np.asarray(galaxy, dtype=float)[good]]
    vals.extend(np.asarray(m, dtype=float)[good] for m in models)
    x = np.concatenate([v[np.isfinite(v)] for v in vals if v.size])
    if x.size == 0:
        return -1.0, 1.0
    lo, hi = np.percentile(x, [1.0, 99.0])
    med = np.median(x)
    mad = 1.4826 * np.median(np.abs(x - med))
    span = max(hi - lo, 8.0 * mad, 1e-3)
    return med - 0.65 * span, med + 0.65 * span


def _plot_representative_bin(
    bid: int,
    cube: dict,
    spectra: dict,
    metrics: Table,
    path: Path,
    residual_threshold: float,
):
    chi = np.asarray(cube["chi2_total"][bid], dtype=float)
    delta, _ = _finite_delta(chi)
    prof = np.min(delta, axis=2)
    safe_chi = np.where(np.isfinite(chi), chi, np.inf)
    best_f_idx = np.argmin(safe_chi, axis=2)
    fa_grid = np.asarray(cube["fA_grid"], dtype=float)
    best_f_map = np.asarray(fa_grid[best_f_idx], dtype=float)
    best_f_map[~np.any(np.isfinite(chi), axis=2)] = np.nan
    va = np.asarray(cube["VA_grid"], dtype=float)
    vb = np.asarray(cube["VB_grid"], dtype=float)
    wave = np.asarray(spectra["wavelength"], dtype=float)
    galaxy = np.asarray(spectra["galaxy"][bid], dtype=float)
    noise = np.asarray(spectra["noise"][bid], dtype=float)
    good = np.asarray(spectra["good"][bid], dtype=bool)
    one = np.asarray(spectra["one_component_model"][bid], dtype=float)
    two = np.asarray(spectra["local_best_two_component_model"][bid], dtype=float)
    _, r1, _ = _residual_metrics(galaxy, one, noise, good)
    _, r2, _ = _residual_metrics(galaxy, two, noise, good)

    fig = plt.figure(figsize=(15, 11))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.0, 1.1, 0.85])
    ax0 = fig.add_subplot(gs[0, 0])
    # log1p keeps topology visible even when pathological pixels make Delta-chi2 enormous.
    im0 = ax0.imshow(
        np.log10(1.0 + prof), origin="lower", aspect="auto",
        extent=[va[0], va[-1], vb[0], vb[-1]], interpolation="nearest",
    )
    ax0.set_title("Profiled velocity surface: log10(1 + Delta chi2)")
    ax0.set_xlabel("V_A (km/s)")
    ax0.set_ylabel("V_B (km/s)")
    fig.colorbar(im0, ax=ax0, label="log10(1 + profiled Delta chi2)")

    ax1 = fig.add_subplot(gs[0, 1])
    im1 = ax1.imshow(
        best_f_map, origin="lower", aspect="auto",
        extent=[va[0], va[-1], vb[0], vb[-1]], interpolation="nearest",
        vmin=float(np.min(fa_grid)), vmax=float(np.max(fa_grid)),
    )
    ax1.set_title("Best discrete f_A at each velocity pair")
    ax1.set_xlabel("V_A (km/s)")
    ax1.set_ylabel("V_B (km/s)")
    fig.colorbar(im1, ax=ax1, label="f_A,RH3")

    ax2 = fig.add_subplot(gs[1, :])
    ax2.plot(wave, galaxy, lw=0.7, label="RH3 spectrum")
    ax2.plot(wave, one, lw=1.1, label="1-component")
    ax2.plot(wave, two, lw=1.1, label="local-best 2-component")
    ylim = _robust_ylim(galaxy, [one, two], good)
    ax2.set_ylim(*ylim)
    clipped = good & np.isfinite(galaxy) & ((galaxy < ylim[0]) | (galaxy > ylim[1]))
    if np.any(clipped):
        ymark = np.where(galaxy[clipped] < ylim[0], ylim[0], ylim[1])
        ax2.scatter(wave[clipped], ymark, marker="x", s=22, label="display-clipped spectrum pixel")
    ax2.set_title("Spectrum and selected models (robust display limits; fitting data unchanged)")
    ax2.set_xlabel("Rest-frame wavelength (Angstrom)")
    ax2.set_ylabel("Normalized flux density")
    ax2.legend(loc="best", fontsize=8)

    ax3 = fig.add_subplot(gs[2, :])
    ax3.plot(wave, r1, lw=0.7, label="1-component standardized residual")
    ax3.plot(wave, r2, lw=0.7, label="2-component standardized residual")
    ax3.axhline(residual_threshold, ls="--", lw=0.8)
    ax3.axhline(-residual_threshold, ls="--", lw=0.8)
    # Residual display is also robust so one 1e5-sigma point does not hide the rest.
    rr = np.concatenate([np.abs(r1[np.isfinite(r1)]), np.abs(r2[np.isfinite(r2)])])
    rlim = max(6.0, float(np.percentile(rr, 99.0)) if rr.size else 6.0)
    rlim = min(50.0, 1.25 * rlim)
    ax3.set_ylim(-rlim, rlim)
    ax3.set_xlabel("Rest-frame wavelength (Angstrom)")
    ax3.set_ylabel("(data-model)/noise")
    ax3.legend(loc="best", fontsize=8)

    row = metrics[bid]
    fig.suptitle(
        f"Script 03b validation bin {bid} | N_eff={row['NEFF_IPR']:.2f} | "
        f"swap p95/scale={row['SWAP_REL_P95']:.3g} | "
        f"2C top1 chi2 leverage={row['TWO_TOP1_CHI2_FRAC']:.3f} | "
        f"global |r| frac={row['TWO_GLOBAL_ABS_R_FRACTION']:.3f} | "
        f"resid/noise={row['TWO_ROBUST_RESIDUAL_TO_NOISE']:.2g} | "
        f"sigma-boundary={row['SIGMA_BOUNDARY_STATE_FRAC']:.3f}",
        fontsize=11,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(path, dpi=160)
    plt.close(fig)


def _select_representative_bins(table: Table, limit: int) -> list[tuple[str, int]]:
    """Choose bins by diagnostic category, deduplicating while preserving priority."""
    n = len(table)
    candidates: list[tuple[str, int]] = []

    def argmax(col):
        a = np.asarray(table[col], dtype=float)
        return int(np.nanargmax(a)) if np.any(np.isfinite(a)) else None

    def argmin(col):
        a = np.asarray(table[col], dtype=float)
        return int(np.nanargmin(a)) if np.any(np.isfinite(a)) else None

    for label, bid in (
        ("worst_two_component_pixel_leverage", argmax("TWO_TOP1_CHI2_FRAC")),
        ("largest_two_component_standardized_residual", argmax("TWO_MAX_ABS_R")),
        ("largest_binwide_standardized_residual_fraction", argmax("TWO_GLOBAL_ABS_R_FRACTION")),
        ("largest_robust_residual_to_noise", argmax("TWO_ROBUST_RESIDUAL_TO_NOISE")),
        ("largest_additive_polynomial_amplitude", argmax("ADDITIVE_POLY_MAX_ABS_COEFF")),
        ("largest_raw_negative_outlier_count", argmax("GALAXY_N_NEGATIVE_OUTLIERS")),
        ("worst_swap_symmetry", argmax("SWAP_REL_P95")),
        ("least_informative_likelihood", argmax("NEFF_IPR")),
        ("most_concentrated_likelihood", argmin("NEFF_IPR")),
        ("largest_sigma_boundary_fraction", argmax("SIGMA_BOUNDARY_STATE_FRAC")),
        ("largest_velocity_edge_likelihood_mass", argmax("VELOCITY_EDGE_WEIGHT_MASS")),
        ("largest_fA_edge_likelihood_mass", argmax("FA_EDGE_WEIGHT_MASS")),
        ("most_epsilon_minimum_plateaus", argmax("N_EPSILON_MIN_PLATEAUS")),
        ("largest_epsilon_minimum_plateau", argmax("LARGEST_EPSILON_PLATEAU_FRAC")),
        ("largest_local_1C_minus_2C_delta_chi2", argmax("LOCAL_DELTA_CHI2_1C_2C")),
        ("smallest_local_1C_minus_2C_delta_chi2", argmin("LOCAL_DELTA_CHI2_1C_2C")),
    ):
        if bid is not None:
            candidates.append((label, bid))

    # Add two useful development candidates when possible: a clean-ish
    # informative bin and a clean-ish broad negative-control candidate.
    top1 = np.asarray(table["TWO_TOP1_CHI2_FRAC"], dtype=float)
    maxr = np.asarray(table["TWO_MAX_ABS_R"], dtype=float)
    neff = np.asarray(table["NEFF_IPR"], dtype=float)
    global_bad = np.asarray(table["GLOBAL_RESIDUAL_PATHOLOGY"], dtype=bool) if "GLOBAL_RESIDUAL_PATHOLOGY" in table.colnames else np.zeros(n, dtype=bool)
    robust_ratio = np.asarray(table["TWO_ROBUST_RESIDUAL_TO_NOISE"], dtype=float) if "TWO_ROBUST_RESIDUAL_TO_NOISE" in table.colnames else np.full(n, np.nan)
    clean = (
        np.isfinite(top1) & np.isfinite(maxr) & np.isfinite(neff) & np.isfinite(robust_ratio)
        & (~global_bad) & (top1 < 0.10) & (maxr < 20.0) & (robust_ratio < ROBUST_RESIDUAL_TO_NOISE_WARN)
    )
    if np.any(clean):
        idx = np.where(clean)[0]
        candidates.append(("clean_informative_candidate", int(idx[np.argmin(neff[idx])])))
        candidates.append(("clean_broad_negative_control_candidate", int(idx[np.argmax(neff[idx])])))

    out: list[tuple[str, int]] = []
    seen = set()
    for label, bid in candidates:
        if bid in seen or bid < 0 or bid >= n:
            continue
        out.append((label, bid))
        seen.add(bid)
        if len(out) >= max(1, int(limit)):
            break
    return out


def _template_validation(templates: dict | None, spectra: dict, logger) -> dict:
    if templates is None:
        logger.warning("XSL_RH3_templates.npz not found; template-normalization cross-check skipped.")
        return {"available": False}
    if "templates" not in templates or "wavelength" not in templates:
        logger.warning("Saved XSL product lacks templates/wavelength arrays; normalization cross-check skipped.")
        return {"available": True, "normalization_checked": False}
    t = np.asarray(templates["templates"], dtype=float)
    tw = np.asarray(templates["wavelength"], dtype=float)
    if t.ndim != 2:
        return {"available": True, "normalization_checked": False, "shape": list(t.shape)}
    if t.shape[0] == tw.size:
        matrix = t
    elif t.shape[1] == tw.size:
        matrix = t.T
    else:
        return {"available": True, "normalization_checked": False, "shape": list(t.shape)}
    norm_range = np.asarray(templates.get("normalization_range", [tw[0], tw[-1]]), dtype=float).ravel()
    use = (tw >= norm_range[0]) & (tw <= norm_range[-1])
    means = np.nanmean(matrix[use], axis=0)
    dev = np.abs(means - 1.0)
    result = {
        "available": True,
        "normalization_checked": True,
        "shape": list(t.shape),
        "n_templates": int(matrix.shape[1]),
        "normalization_range_angstrom": [float(norm_range[0]), float(norm_range[-1])],
        "median_abs_mean_minus_one": float(np.nanmedian(dev)),
        "max_abs_mean_minus_one": float(np.nanmax(dev)),
    }
    logger.info(
        "XSL normalization check: %d templates | median |mean-1|=%.3g | max=%.3g",
        result["n_templates"], result["median_abs_mean_minus_one"], result["max_abs_mean_minus_one"],
    )
    return result


def main() -> int:
    args = _parser().parse_args()
    script3_run = Path(args.script3_run).expanduser().resolve()
    if not script3_run.is_dir():
        raise FileNotFoundError(f"Script-3 run does not exist: {script3_run}")
    output_dir = (
        Path(args.output_dir).expanduser().resolve()
        if args.output_dir is not None
        else script3_run / "validation" / "03b"
    )
    products_dir = output_dir / "products"
    figures_dir = output_dir / "figures"
    rep_dir = figures_dir / "representative_bins"
    for p in (products_dir, figures_dir, rep_dir, output_dir / "metadata"):
        p.mkdir(parents=True, exist_ok=True)
    logger = _setup_logger(output_dir)

    _section(logger, "CRD_DAP SCRIPT 03b: VALIDATE RH3 LIKELIHOOD CUBES", "=")
    logger.info("Source Script-3 run: %s", script3_run)
    logger.info("Validation output: %s", output_dir)
    logger.info("This stage is read-only and performs no pPXF fits.")
    logger.info("Production atmospheric-mask selection is delegated to Script 03d; 03b recurrence products are diagnostic/QC only unless a legacy --candidate-mask-mode is explicitly requested.")

    _section(logger, "1. Load and mechanically validate Script-3 products")
    (
        cube_path, spectra_path, template_path, cube, spectra, templates, manifest, source_summary
    ) = _load_products(script3_run)
    nbin, nva, nvb, nfa, nlog = _validate_shapes(cube, spectra, manifest)
    va = np.asarray(cube["VA_grid"], dtype=float)
    vb = np.asarray(cube["VB_grid"], dtype=float)
    fa = np.asarray(cube["fA_grid"], dtype=float)
    wave = np.asarray(spectra["wavelength"], dtype=float)
    redshift = _resolve_redshift(args.redshift, manifest, logger, args.config)
    template_medium = str(manifest.get("template_wavelength_medium", "unknown")).strip().lower()
    science_medium = str(manifest.get("science_wavelength_medium", "unknown")).strip().lower()
    observed_air_wave, observed_vacuum_wave = _observed_wavelength_grids(wave, redshift, template_medium, logger)
    observed_wave = observed_air_wave
    global_r_thr = float(args.global_residual_abs_r_threshold)
    global_r_frac = float(args.global_residual_pixel_fraction)
    raw_negative_sigma = float(args.raw_negative_outlier_sigma)
    if not np.isfinite(global_r_thr) or global_r_thr <= 0:
        raise ValueError("--global-residual-abs-r-threshold must be finite and > 0.")
    if not np.isfinite(global_r_frac) or not (0 < global_r_frac <= 1):
        raise ValueError("--global-residual-pixel-fraction must be in (0, 1].")
    if not np.isfinite(raw_negative_sigma) or raw_negative_sigma <= 0:
        raise ValueError("--raw-negative-outlier-sigma must be finite and > 0.")
    logger.info("Likelihood shape: bins=%d, VA=%d, VB=%d, fA=%d (%d states/bin)", nbin, nva, nvb, nfa, nva * nvb * nfa)
    logger.info("Saved log spectrum: %d pixels, %.2f--%.2f A rest", nlog, wave[0], wave[-1])
    if np.isfinite(redshift) and np.all(np.isfinite(observed_air_wave)):
        logger.info(
            "Observed-frame diagnostic grid: air %.2f--%.2f A; vacuum %.2f--%.2f A (z=%.8f)",
            observed_air_wave[0], observed_air_wave[-1], observed_vacuum_wave[0], observed_vacuum_wave[-1], redshift,
        )
    if np.any(np.asarray(cube["n_failed_states"], dtype=int) > 0):
        logger.warning("Some bins contain failed grid states; see per-bin table.")
    else:
        logger.info("Grid-state execution check: PASS | zero failed states in all %d bins", nbin)

    template_result = _template_validation(templates, spectra, logger)

    _section(logger, "2. A/B component-swap symmetry and likelihood information content")
    chi_all = np.asarray(cube["chi2_total"], dtype=np.float64)
    chi_swap = _component_swap_cube(chi_all, va, vb, fa)
    swap_abs_all = np.abs(chi_all - chi_swap)
    # Under A<->B exchange, sigma_A at one coordinate should equal sigma_B
    # at the swapped coordinate (and vice versa). Fit-status and sigma-boundary
    # masks should likewise be symmetric. These are additional bookkeeping
    # checks beyond chi-square symmetry itself.
    sigma_a_all = np.asarray(cube["sigma_A"], dtype=float)
    sigma_b_all = np.asarray(cube["sigma_B"], dtype=float)
    sigma_b_swap = _component_swap_cube(sigma_b_all, va, vb, fa)
    sigma_a_swap = _component_swap_cube(sigma_a_all, va, vb, fa)
    status_all = np.asarray(cube["fit_status"])
    status_swap = _component_swap_cube(status_all, va, vb, fa)
    boundary_all = np.asarray(cube["sigma_boundary"])
    boundary_swap = _component_swap_cube(boundary_all, va, vb, fa)

    # Per-bin arrays.
    cols: OrderedDict[str, np.ndarray] = OrderedDict()
    cols["BIN_ID"] = np.arange(nbin, dtype=int)
    float_cols = [
        "CHI2_MIN", "LOCAL_DELTA_CHI2_1C_2C",
        "SWAP_ABS_MEDIAN", "SWAP_ABS_P95", "SWAP_ABS_P99", "SWAP_ABS_MAX",
        "SWAP_NEAR_MIN_ABS_P95", "SWAP_REL_P95", "SWAP_LOCAL_MIN_GAP",
        "SIGMA_SWAP_ABS_MEDIAN", "SIGMA_SWAP_ABS_P95", "SIGMA_SWAP_ABS_MAX",
        "SWAP_STATUS_MISMATCH_FRAC", "SWAP_SIGMA_BOUNDARY_MISMATCH_FRAC",
        "TOPOLOGY_EPSILON_CHI2",
        "NEFF_IPR", "LIKELIHOOD_ENTROPY", "NEFF_ENTROPY",
        "VELOCITY_EDGE_WEIGHT_MASS", "FA_EDGE_WEIGHT_MASS", "ANY_GRID_EDGE_WEIGHT_MASS",
        "SIGMA_BOUNDARY_STATE_FRAC", "LARGEST_EPSILON_PLATEAU_FRAC",
        "ONE_CHI2_RECOMPUTED", "ONE_CHI2_SAVED", "ONE_CHI2_RELERR", "ONE_MAX_ABS_R",
        "ONE_TOP1_CHI2_FRAC", "ONE_TOP5_CHI2_FRAC", "ONE_TOP10_CHI2_FRAC", "ONE_TOP1PCT_CHI2_FRAC",
        "TWO_CHI2_RECOMPUTED", "TWO_CHI2_SAVED", "TWO_CHI2_RELERR", "TWO_MAX_ABS_R",
        "TWO_TOP1_CHI2_FRAC", "TWO_TOP5_CHI2_FRAC", "TWO_TOP10_CHI2_FRAC", "TWO_TOP1PCT_CHI2_FRAC",
        "GALAXY_MIN", "GALAXY_P01", "GALAXY_MEDIAN", "GALAXY_P99", "GALAXY_MAX",
        "GALAXY_NEGATIVE_FRAC", "GALAXY_ROBUST_SCATTER", "GALAXY_MIN_ROBUST_Z",
        "GALAXY_NEGATIVE_OUTLIER_FRAC",
        "NOISE_MIN", "NOISE_P01", "NOISE_MEDIAN", "NOISE_P99", "NOISE_MAX", "NOISE_NONPOSITIVE_FRAC",
        "ONE_RESIDUAL_MEDIAN", "ONE_RESIDUAL_MAD_SIGMA", "ONE_MEDIAN_ABS_RESIDUAL",
        "ONE_MEDIAN_ABS_R", "ONE_ROBUST_RESIDUAL_TO_NOISE", "ONE_MODEL_MEDIAN_OFFSET_NOISE",
        "TWO_RESIDUAL_MEDIAN", "TWO_RESIDUAL_MAD_SIGMA", "TWO_MEDIAN_ABS_RESIDUAL",
        "TWO_MEDIAN_ABS_R", "TWO_ROBUST_RESIDUAL_TO_NOISE", "TWO_MODEL_MEDIAN_OFFSET_NOISE",
        "ONE_MODEL_MIN", "ONE_MODEL_P01", "ONE_MODEL_MEDIAN", "ONE_MODEL_P99", "ONE_MODEL_MAX",
        "TWO_MODEL_MIN", "TWO_MODEL_P01", "TWO_MODEL_MEDIAN", "TWO_MODEL_P99", "TWO_MODEL_MAX",
        "TWO_GLOBAL_ABS_R_FRACTION", "ADDITIVE_POLY_MAX_ABS_COEFF", "ADDITIVE_POLY_L2",
    ]
    int_cols = [
        "N_FINITE_STATES", "N_MASS_90", "N_MASS_95", "N_MASS_99", "N_MASS_995", "N_MASS_999",
        "N_STRICT_MIN_PLATEAUS", "N_EPSILON_MIN_PLATEAUS", "LARGEST_EPSILON_PLATEAU_CELLS",
        "N_FAILED_STATES", "N_SIGMA_BOUNDARY_STATES", "N_GOOD_PIXELS", "N_VALID_RESIDUAL_PIXELS",
        "GALAXY_N_NEGATIVE_OUTLIERS",
    ]
    bool_cols = [
        "LOCAL_BEST_VELOCITY_EDGE", "LOCAL_BEST_FA_EDGE", "LOCAL_BEST_SIGMA_EDGE",
        "LEVERAGE_WARNING", "GLOBAL_RESIDUAL_PATHOLOGY", "RESIDUAL_SCALE_WARNING",
        "MODEL_MEDIAN_OFFSET_WARNING",
    ]
    for name in float_cols:
        cols[name] = np.full(nbin, np.nan, dtype=float)
    for name in int_cols:
        cols[name] = np.zeros(nbin, dtype=int)
    for name in bool_cols:
        cols[name] = np.zeros(nbin, dtype=bool)
    for prefix in ("ONE", "TWO"):
        for t in RESIDUAL_THRESHOLDS:
            cols[f"{prefix}_N_ABS_R_GT_{int(t)}"] = np.zeros(nbin, dtype=int)

    one_r_all = np.full((nbin, nlog), np.nan, dtype=np.float32)
    two_r_all = np.full((nbin, nlog), np.nan, dtype=np.float32)
    one_chi_pix_all = np.full((nbin, nlog), np.nan, dtype=np.float32)
    two_chi_pix_all = np.full((nbin, nlog), np.nan, dtype=np.float32)

    explicit_topo = None
    if str(args.topology_epsilon).strip().lower() != "auto":
        explicit_topo = float(args.topology_epsilon)
        if explicit_topo < 0:
            raise ValueError("--topology-epsilon must be non-negative or 'auto'.")

    sigma_bounds = manifest.get("sigma_bounds_kms", [np.nan, np.nan])
    sigma_lo, sigma_hi = float(sigma_bounds[0]), float(sigma_bounds[1])
    local_sa = np.asarray(cube.get("local_best_sigma_A", np.full(nbin, np.nan)), dtype=float)
    local_sb = np.asarray(cube.get("local_best_sigma_B", np.full(nbin, np.nan)), dtype=float)
    best_va = np.asarray(cube.get("local_best_VA", np.full(nbin, np.nan)), dtype=float)
    best_vb = np.asarray(cube.get("local_best_VB", np.full(nbin, np.nan)), dtype=float)
    best_fa = np.asarray(cube.get("local_best_fA", np.full(nbin, np.nan)), dtype=float)

    recurrent_thr = float(args.recurrent_residual_threshold)
    one_recurrent_counts = np.zeros(nlog, dtype=int)
    two_recurrent_counts = np.zeros(nlog, dtype=int)
    raw_negative_recurrent_counts = np.zeros(nlog, dtype=int)
    raw_negative_outlier_all = np.zeros((nbin, nlog), dtype=bool)

    for bid in range(nbin):
        chi = chi_all[bid]
        delta, cmin = _finite_delta(chi)
        weights = _likelihood_weights(delta)
        finite = np.isfinite(chi)
        cols["CHI2_MIN"][bid] = cmin
        cols["N_FINITE_STATES"][bid] = int(np.count_nonzero(finite))
        cols["LOCAL_DELTA_CHI2_1C_2C"][bid] = float(cube["one_chi2_total"][bid] - cube["local_best_chi2_total"][bid])

        swap_abs = swap_abs_all[bid]
        vals = swap_abs[finite & np.isfinite(swap_abs)]
        cols["SWAP_ABS_MEDIAN"][bid] = _percentile_or_nan(vals, 50)
        cols["SWAP_ABS_P95"][bid] = _percentile_or_nan(vals, 95)
        cols["SWAP_ABS_P99"][bid] = _percentile_or_nan(vals, 99)
        cols["SWAP_ABS_MAX"][bid] = _percentile_or_nan(vals, 100)
        near = finite & np.isfinite(swap_abs) & (np.minimum(delta, _finite_delta(chi_swap[bid])[0]) <= SWAP_NEAR_MIN_DELTA_CHI2)
        cols["SWAP_NEAR_MIN_ABS_P95"][bid] = _percentile_or_nan(swap_abs[near], 95)
        dyn_scale = max(1.0, _percentile_or_nan(delta[finite], 95))
        cols["SWAP_REL_P95"][bid] = cols["SWAP_ABS_P95"][bid] / dyn_scale if np.isfinite(dyn_scale) else np.nan
        try:
            idx = np.unravel_index(np.nanargmin(chi), chi.shape)
            swapped_min_chi = chi_swap[bid][idx]
            cols["SWAP_LOCAL_MIN_GAP"][bid] = abs(float(chi[idx]) - float(swapped_min_chi))
        except Exception:
            pass
        sigma_swap_diff = np.concatenate([
            np.abs(sigma_a_all[bid] - sigma_b_swap[bid]).ravel(),
            np.abs(sigma_b_all[bid] - sigma_a_swap[bid]).ravel(),
        ])
        cols["SIGMA_SWAP_ABS_MEDIAN"][bid] = _percentile_or_nan(sigma_swap_diff, 50)
        cols["SIGMA_SWAP_ABS_P95"][bid] = _percentile_or_nan(sigma_swap_diff, 95)
        cols["SIGMA_SWAP_ABS_MAX"][bid] = _percentile_or_nan(sigma_swap_diff, 100)
        valid_swap = finite & np.isfinite(chi_swap[bid])
        if np.any(valid_swap):
            cols["SWAP_STATUS_MISMATCH_FRAC"][bid] = float(np.mean(status_all[bid][valid_swap] != status_swap[bid][valid_swap]))
            cols["SWAP_SIGMA_BOUNDARY_MISMATCH_FRAC"][bid] = float(np.mean(boundary_all[bid][valid_swap] != boundary_swap[bid][valid_swap]))
        eps = explicit_topo if explicit_topo is not None else _robust_swap_epsilon(delta, swap_abs)
        cols["TOPOLOGY_EPSILON_CHI2"][bid] = eps

        neff_ipr, entropy, neff_entropy = _entropy_metrics(weights)
        cols["NEFF_IPR"][bid] = neff_ipr
        cols["LIKELIHOOD_ENTROPY"][bid] = entropy
        cols["NEFF_ENTROPY"][bid] = neff_entropy
        counts = _mass_cell_counts(weights)
        for name, value in zip(("N_MASS_90", "N_MASS_95", "N_MASS_99", "N_MASS_995", "N_MASS_999"), counts):
            cols[name][bid] = value
        velmass, famass, anymass = _fraction_on_grid_edges(weights)
        cols["VELOCITY_EDGE_WEIGHT_MASS"][bid] = velmass
        cols["FA_EDGE_WEIGHT_MASS"][bid] = famass
        cols["ANY_GRID_EDGE_WEIGHT_MASS"][bid] = anymass

        nstrict, nthick, largest, largest_frac = _topology_metrics(chi, eps)
        cols["N_STRICT_MIN_PLATEAUS"][bid] = nstrict
        cols["N_EPSILON_MIN_PLATEAUS"][bid] = nthick
        cols["LARGEST_EPSILON_PLATEAU_CELLS"][bid] = largest
        cols["LARGEST_EPSILON_PLATEAU_FRAC"][bid] = largest_frac

        nfailed = int(np.asarray(cube["n_failed_states"])[bid])
        nbound = int(np.asarray(cube["n_sigma_boundary_states"])[bid])
        cols["N_FAILED_STATES"][bid] = nfailed
        cols["N_SIGMA_BOUNDARY_STATES"][bid] = nbound
        cols["SIGMA_BOUNDARY_STATE_FRAC"][bid] = nbound / float(nva * nvb * nfa)
        cols["LOCAL_BEST_VELOCITY_EDGE"][bid] = (
            np.isclose(best_va[bid], va[0]) or np.isclose(best_va[bid], va[-1]) or
            np.isclose(best_vb[bid], vb[0]) or np.isclose(best_vb[bid], vb[-1])
        )
        cols["LOCAL_BEST_FA_EDGE"][bid] = np.isclose(best_fa[bid], fa[0]) or np.isclose(best_fa[bid], fa[-1])
        cols["LOCAL_BEST_SIGMA_EDGE"][bid] = (
            (np.isfinite(sigma_lo) and (np.isclose(local_sa[bid], sigma_lo, atol=1e-3) or np.isclose(local_sb[bid], sigma_lo, atol=1e-3))) or
            (np.isfinite(sigma_hi) and (np.isclose(local_sa[bid], sigma_hi, atol=1e-3) or np.isclose(local_sb[bid], sigma_hi, atol=1e-3)))
        )

        galaxy = np.asarray(spectra["galaxy"][bid], dtype=float)
        noise = np.asarray(spectra["noise"][bid], dtype=float)
        good = np.asarray(spectra["good"][bid], dtype=bool)
        one_model = np.asarray(spectra["one_component_model"][bid], dtype=float)
        two_model = np.asarray(spectra["local_best_two_component_model"][bid], dtype=float)
        one_m, r1, c1 = _residual_metrics(galaxy, one_model, noise, good)
        two_m, r2, c2 = _residual_metrics(galaxy, two_model, noise, good)
        one_diag = _noise_residual_diagnostics(galaxy, one_model, noise, good)
        two_diag = _noise_residual_diagnostics(galaxy, two_model, noise, good)
        cols["N_GOOD_PIXELS"][bid] = int(two_diag["n_good"])
        cols["N_VALID_RESIDUAL_PIXELS"][bid] = int(two_diag["n_valid"])
        for name in ("noise_min", "noise_p01", "noise_median", "noise_p99", "noise_max", "noise_nonpositive_frac"):
            cols[name.upper()][bid] = two_diag[name]
        cols["GALAXY_NEGATIVE_FRAC"][bid] = two_diag["galaxy_negative_frac"]
        for prefix, diag in (("ONE", one_diag), ("TWO", two_diag)):
            cols[f"{prefix}_RESIDUAL_MEDIAN"][bid] = diag["residual_median"]
            cols[f"{prefix}_RESIDUAL_MAD_SIGMA"][bid] = diag["residual_mad_sigma"]
            cols[f"{prefix}_MEDIAN_ABS_RESIDUAL"][bid] = diag["median_abs_residual"]
            cols[f"{prefix}_MEDIAN_ABS_R"][bid] = diag["median_abs_r"]
            cols[f"{prefix}_ROBUST_RESIDUAL_TO_NOISE"][bid] = diag["robust_residual_to_noise"]
            cols[f"{prefix}_MODEL_MEDIAN_OFFSET_NOISE"][bid] = diag["model_median_offset_noise"]
            cols[f"{prefix}_MODEL_MIN"][bid] = diag["model_min"]
            cols[f"{prefix}_MODEL_P01"][bid] = diag["model_p01"]
            cols[f"{prefix}_MODEL_MEDIAN"][bid] = diag["model_median"]
            cols[f"{prefix}_MODEL_P99"][bid] = diag["model_p99"]
            cols[f"{prefix}_MODEL_MAX"][bid] = diag["model_max"]

        valid_r2 = np.isfinite(r2)
        if np.any(valid_r2):
            cols["TWO_GLOBAL_ABS_R_FRACTION"][bid] = float(np.mean(np.abs(r2[valid_r2]) > global_r_thr))
        cols["GLOBAL_RESIDUAL_PATHOLOGY"][bid] = (
            np.isfinite(cols["TWO_GLOBAL_ABS_R_FRACTION"][bid])
            and cols["TWO_GLOBAL_ABS_R_FRACTION"][bid] >= global_r_frac
        )
        cols["RESIDUAL_SCALE_WARNING"][bid] = (
            np.isfinite(cols["TWO_ROBUST_RESIDUAL_TO_NOISE"][bid])
            and cols["TWO_ROBUST_RESIDUAL_TO_NOISE"][bid] > ROBUST_RESIDUAL_TO_NOISE_WARN
        )
        cols["MODEL_MEDIAN_OFFSET_WARNING"][bid] = (
            np.isfinite(cols["TWO_MODEL_MEDIAN_OFFSET_NOISE"][bid])
            and cols["TWO_MODEL_MEDIAN_OFFSET_NOISE"][bid] > MODEL_MEDIAN_OFFSET_NOISE_WARN
        )
        if "local_best_additive_polyweights" in spectra:
            pw = np.asarray(spectra["local_best_additive_polyweights"][bid], dtype=float)
            pw = pw[np.isfinite(pw)]
            if pw.size:
                cols["ADDITIVE_POLY_MAX_ABS_COEFF"][bid] = float(np.max(np.abs(pw)))
                cols["ADDITIVE_POLY_L2"][bid] = float(np.sqrt(np.sum(pw * pw)))

        one_r_all[bid] = r1.astype(np.float32)
        two_r_all[bid] = r2.astype(np.float32)
        one_chi_pix_all[bid] = c1.astype(np.float32)
        two_chi_pix_all[bid] = c2.astype(np.float32)
        one_recurrent_counts += np.asarray(np.abs(r1) > recurrent_thr, dtype=int)
        two_recurrent_counts += np.asarray(np.abs(r2) > recurrent_thr, dtype=int)

        for prefix, m, saved in (
            ("ONE", one_m, float(cube["one_chi2_total"][bid])),
            ("TWO", two_m, float(cube["local_best_chi2_total"][bid])),
        ):
            cols[f"{prefix}_CHI2_RECOMPUTED"][bid] = m["chi2_recomputed"]
            cols[f"{prefix}_CHI2_SAVED"][bid] = saved
            denom = max(1.0, abs(saved))
            cols[f"{prefix}_CHI2_RELERR"][bid] = abs(m["chi2_recomputed"] - saved) / denom
            cols[f"{prefix}_MAX_ABS_R"][bid] = m["max_abs_r"]
            cols[f"{prefix}_TOP1_CHI2_FRAC"][bid] = m["top1_frac"]
            cols[f"{prefix}_TOP5_CHI2_FRAC"][bid] = m["top5_frac"]
            cols[f"{prefix}_TOP10_CHI2_FRAC"][bid] = m["top10_frac"]
            cols[f"{prefix}_TOP1PCT_CHI2_FRAC"][bid] = m["top1pct_frac"]
            for t in RESIDUAL_THRESHOLDS:
                cols[f"{prefix}_N_ABS_R_GT_{int(t)}"][bid] = m[f"n_abs_r_gt_{int(t)}"]

        gvalid = good & np.isfinite(galaxy)
        if np.any(gvalid):
            gidx = np.where(gvalid)[0]
            gv = galaxy[gvalid]
            gmed = float(np.median(gv))
            gscatter = _mad_sigma(gv)
            cols["GALAXY_MIN"][bid] = float(np.min(gv))
            cols["GALAXY_P01"][bid] = float(np.percentile(gv, 1.0))
            cols["GALAXY_MEDIAN"][bid] = gmed
            cols["GALAXY_P99"][bid] = float(np.percentile(gv, 99.0))
            cols["GALAXY_MAX"][bid] = float(np.max(gv))
            cols["GALAXY_ROBUST_SCATTER"][bid] = gscatter
            if np.isfinite(gscatter) and gscatter > 0:
                gz = (gv - gmed) / gscatter
                neg_out = gz < -raw_negative_sigma
                n_neg = int(np.count_nonzero(neg_out))
                cols["GALAXY_N_NEGATIVE_OUTLIERS"][bid] = n_neg
                cols["GALAXY_NEGATIVE_OUTLIER_FRAC"][bid] = n_neg / float(gv.size)
                cols["GALAXY_MIN_ROBUST_Z"][bid] = float(np.min(gz))
                if n_neg:
                    raw_negative_outlier_all[bid, gidx[neg_out]] = True
                    raw_negative_recurrent_counts[gidx[neg_out]] += 1

        cols["LEVERAGE_WARNING"][bid] = (
            (np.isfinite(cols["TWO_TOP1_CHI2_FRAC"][bid]) and cols["TWO_TOP1_CHI2_FRAC"][bid] > LEVERAGE_TOP1_WARN) or
            (np.isfinite(cols["TWO_TOP10_CHI2_FRAC"][bid]) and cols["TWO_TOP10_CHI2_FRAC"][bid] > LEVERAGE_TOP10_WARN) or
            (np.isfinite(cols["TWO_MAX_ABS_R"][bid]) and cols["TWO_MAX_ABS_R"][bid] > LEVERAGE_MAX_R_WARN)
        )

    metrics = Table(cols)

    # Augment with useful Script-2 S/N columns when available, without making
    # their names an API requirement.
    bin_map = source_table = resolved_script2 = None
    if not args.no_spatial_maps:
        bin_map, source_table, resolved_script2 = _find_bin_map_and_source_table(
            manifest, args.script2_run, nbin, logger
        )
    if source_table is not None:
        for name in source_table.colnames:
            upper = name.upper()
            if ("SN" in upper or "SNR" in upper) and name not in metrics.colnames:
                try:
                    metrics[name] = source_table[name]
                except Exception:
                    pass

    _section(logger, "3. Mechanical chi-square and f_A bookkeeping cross-checks")
    max_one_rel = _percentile_or_nan(np.asarray(metrics["ONE_CHI2_RELERR"], dtype=float), 100)
    max_two_rel = _percentile_or_nan(np.asarray(metrics["TWO_CHI2_RELERR"], dtype=float), 100)
    chi2_recompute_pass = (
        np.isfinite(max_one_rel) and np.isfinite(max_two_rel) and
        max_one_rel <= CHI2_RECOMPUTE_RTOL and max_two_rel <= CHI2_RECOMPUTE_RTOL
    )
    logger.info(
        "Selected-fit chi2 recomputation: %s | max relative error 1C=%.3g, 2C=%.3g",
        "PASS" if chi2_recompute_pass else "WARN", max_one_rel, max_two_rel,
    )

    fraction_check = {"available": False}
    weights2 = np.asarray(spectra["local_best_two_component_weights"], dtype=float)
    if templates is not None and "templates" in templates and weights2.ndim == 2:
        t = np.asarray(templates["templates"])
        tw = np.asarray(templates.get("wavelength", []))
        n_basis = t.shape[1] if t.ndim == 2 and t.shape[0] == tw.size else (t.shape[0] if t.ndim == 2 and t.shape[1] == tw.size else 0)
        if n_basis > 0 and weights2.shape[1] >= 2 * n_basis:
            wa = np.sum(weights2[:, :n_basis], axis=1)
            wb = np.sum(weights2[:, n_basis:2*n_basis], axis=1)
            with np.errstate(divide="ignore", invalid="ignore"):
                faw = wa / (wa + wb)
            saved_achieved = np.asarray(cube.get("local_best_fA_achieved", np.full(nbin, np.nan)), dtype=float)
            diff = np.abs(faw - saved_achieved)
            fraction_check = {
                "available": True,
                "n_basis": int(n_basis),
                "median_abs_difference": _percentile_or_nan(diff, 50),
                "max_abs_difference": _percentile_or_nan(diff, 100),
            }
            logger.info(
                "Duplicated-template f_A bookkeeping: max |f_A(weights)-saved achieved|=%.3g",
                fraction_check["max_abs_difference"],
            )

    _section(logger, "4. Summarize likelihood topology, boundaries, and spectral leverage")
    swap_rel_p95_global = _percentile_or_nan(np.asarray(metrics["SWAP_REL_P95"], dtype=float), 95)
    if swap_rel_p95_global <= SWAP_REL_PASS:
        swap_status = "PASS"
    elif swap_rel_p95_global <= SWAP_REL_WARN:
        swap_status = "WARN"
    else:
        swap_status = "FAIL"
    leverage_frac = float(np.mean(np.asarray(metrics["LEVERAGE_WARNING"], dtype=bool)))
    sigma_bound_warn_frac = float(np.mean(np.asarray(metrics["SIGMA_BOUNDARY_STATE_FRAC"], dtype=float) > SIGMA_BOUNDARY_STATE_WARN))
    velocity_edge_warn_frac = float(np.mean(np.asarray(metrics["VELOCITY_EDGE_WEIGHT_MASS"], dtype=float) > GRID_EDGE_WEIGHT_WARN))
    fa_edge_warn_frac = float(np.mean(np.asarray(metrics["FA_EDGE_WEIGHT_MASS"], dtype=float) > GRID_EDGE_WEIGHT_WARN))
    logger.info("Swap-symmetry status: %s | 95th percentile across-bin SWAP_REL_P95=%.3g", swap_status, swap_rel_p95_global)
    logger.info("Likelihood information: median N_eff(IPR)=%.2f states; median N_95=%d states",
                _percentile_or_nan(np.asarray(metrics["NEFF_IPR"], dtype=float), 50),
                int(np.nanmedian(np.asarray(metrics["N_MASS_95"], dtype=float))))
    logger.info("Spectral leverage warning fraction: %.1f%% of bins", 100.0 * leverage_frac)
    logger.info("Sigma-boundary >50%% state fraction: %.1f%% of bins", 100.0 * sigma_bound_warn_frac)
    logger.info("Velocity-edge likelihood mass >25%%: %.1f%% of bins", 100.0 * velocity_edge_warn_frac)
    logger.info("f_A-edge likelihood mass >25%%: %.1f%% of bins", 100.0 * fa_edge_warn_frac)

    global_bad = np.asarray(metrics["GLOBAL_RESIDUAL_PATHOLOGY"], dtype=bool)
    raw_negative_bin = np.asarray(metrics["GALAXY_N_NEGATIVE_OUTLIERS"], dtype=int) > 0
    residual_scale_bad = np.asarray(metrics["RESIDUAL_SCALE_WARNING"], dtype=bool)
    model_offset_bad = np.asarray(metrics["MODEL_MEDIAN_OFFSET_WARNING"], dtype=bool)
    n_global_bad = int(np.count_nonzero(global_bad))
    n_global_clean = int(nbin - n_global_bad)
    global_bad_frac = n_global_bad / float(nbin)
    logger.info(
        "Bin-wide residual pathology: %d/%d bins (%.1f%%) have >=%.0f%% of valid pixels with |r_2C| > %.1f",
        n_global_bad, nbin, 100.0 * global_bad_frac, 100.0 * global_r_frac, global_r_thr,
    )
    logger.info(
        "Model-independent raw negative outliers: %d/%d bins contain at least one pixel below median - %.1f robust-MAD sigma",
        int(np.count_nonzero(raw_negative_bin)), nbin, raw_negative_sigma,
    )
    n_overlap_raw_global = int(np.count_nonzero(raw_negative_bin & global_bad))
    frac_global_with_raw = n_overlap_raw_global / float(max(1, n_global_bad))
    frac_raw_with_global = n_overlap_raw_global / float(max(1, np.count_nonzero(raw_negative_bin)))
    logger.info(
        "Raw-outlier/global-residual overlap: %d bins | %.1f%% of GLOBAL_RESIDUAL_PATHOLOGY bins have a raw negative outlier; "
        "%.1f%% of raw-outlier bins are globally pathological",
        n_overlap_raw_global, 100.0 * frac_global_with_raw, 100.0 * frac_raw_with_global,
    )
    if np.any(raw_negative_recurrent_counts > 0):
        iraw = int(np.argmax(raw_negative_recurrent_counts))
        if np.isfinite(redshift) and np.isfinite(observed_air_wave[iraw]):
            logger.info(
                "Most recurrent raw negative-outlier wavelength: rest %.2f A -> observed air %.2f A | %d bins",
                wave[iraw], observed_air_wave[iraw], int(raw_negative_recurrent_counts[iraw]),
            )
        else:
            logger.info(
                "Most recurrent raw negative-outlier wavelength: rest %.2f A | %d bins",
                wave[iraw], int(raw_negative_recurrent_counts[iraw]),
            )
    logger.info(
        "Residual-scale diagnostics: %.1f%% of bins have robust residual scatter / median formal noise > %.1f; "
        "%.1f%% have |median(data-model)| / median noise > %.1f",
        100.0 * float(np.mean(residual_scale_bad)), ROBUST_RESIDUAL_TO_NOISE_WARN,
        100.0 * float(np.mean(model_offset_bad)), MODEL_MEDIAN_OFFSET_NOISE_WARN,
    )
    for label, mask in (("global-pathology", global_bad), ("non-pathology", ~global_bad)):
        if np.any(mask):
            logger.info(
                "%s bins: n=%d | median noise=%.4g | median residual MAD=%.4g | "
                "median robust residual/noise=%.3g | median |r|=%.3g | median galaxy min=%.4g | median model=%.4g",
                label, int(np.count_nonzero(mask)),
                _percentile_or_nan(np.asarray(metrics["NOISE_MEDIAN"], dtype=float)[mask], 50),
                _percentile_or_nan(np.asarray(metrics["TWO_RESIDUAL_MAD_SIGMA"], dtype=float)[mask], 50),
                _percentile_or_nan(np.asarray(metrics["TWO_ROBUST_RESIDUAL_TO_NOISE"], dtype=float)[mask], 50),
                _percentile_or_nan(np.asarray(metrics["TWO_MEDIAN_ABS_R"], dtype=float)[mask], 50),
                _percentile_or_nan(np.asarray(metrics["GALAXY_MIN"], dtype=float)[mask], 50),
                _percentile_or_nan(np.asarray(metrics["TWO_MODEL_MEDIAN"], dtype=float)[mask], 50),
            )

    # Recurrent-wavelength counts split by the bin-wide pathology flag. This
    # separates a broad baseline problem from wavelength-localized features.
    one_hit = np.isfinite(one_r_all) & (np.abs(one_r_all) > recurrent_thr)
    two_hit = np.isfinite(two_r_all) & (np.abs(two_r_all) > recurrent_thr)
    one_recurrent_clean = np.sum(one_hit[~global_bad], axis=0, dtype=int) if n_global_clean else np.zeros(nlog, dtype=int)
    two_recurrent_clean = np.sum(two_hit[~global_bad], axis=0, dtype=int) if n_global_clean else np.zeros(nlog, dtype=int)
    one_recurrent_bad = np.sum(one_hit[global_bad], axis=0, dtype=int) if n_global_bad else np.zeros(nlog, dtype=int)
    two_recurrent_bad = np.sum(two_hit[global_bad], axis=0, dtype=int) if n_global_bad else np.zeros(nlog, dtype=int)
    baseline_one = float(np.median(one_recurrent_counts))
    baseline_two = float(np.median(two_recurrent_counts))
    baseline_one_clean = float(np.median(one_recurrent_clean)) if n_global_clean else 0.0
    baseline_two_clean = float(np.median(two_recurrent_clean)) if n_global_clean else 0.0
    one_excess = one_recurrent_counts.astype(float) - baseline_one
    two_excess = two_recurrent_counts.astype(float) - baseline_two
    one_excess_clean = one_recurrent_clean.astype(float) - baseline_one_clean
    two_excess_clean = two_recurrent_clean.astype(float) - baseline_two_clean
    logger.info(
        "Recurrent-wavelength baseline (median across wavelength): all bins 1C=%.1f, 2C=%.1f; "
        "non-pathology bins 1C=%.1f, 2C=%.1f",
        baseline_one, baseline_two, baseline_one_clean, baseline_two_clean,
    )
    rank_score = np.maximum(one_excess_clean, two_excess_clean) if n_global_clean else np.maximum(one_excess, two_excess)
    if np.any(np.isfinite(rank_score)):
        imax = int(np.nanargmax(rank_score))
        if np.isfinite(redshift):
            logger.info(
                "Strongest localized recurrent feature: rest %.2f A -> observed air %.2f A | excess non-pathology bins=%.0f",
                wave[imax], observed_air_wave[imax], rank_score[imax],
            )
        else:
            logger.info(
                "Strongest localized recurrent feature: rest %.2f A | excess non-pathology bins=%.0f",
                wave[imax], rank_score[imax],
            )

    metrics_path = products_dir / "RH3_03b_validation_summary.ecsv"
    metrics.write(metrics_path, format="ascii.ecsv", overwrite=True)

    np.savez_compressed(
        products_dir / "RH3_03b_raw_diagnostics.npz",
        wavelength=wave,
        wavelength_rest=wave,
        wavelength_observed_air=observed_air_wave,
        wavelength_observed_vacuum=observed_vacuum_wave,
        redshift=np.asarray(redshift),
        standardized_residual_one=one_r_all,
        standardized_residual_two=two_r_all,
        chi2_pixel_one=one_chi_pix_all,
        chi2_pixel_two=two_chi_pix_all,
        swap_abs_delta_chi2=np.asarray(swap_abs_all, dtype=np.float32),
        recurrent_abs_residual_threshold=np.asarray(recurrent_thr),
        recurrent_bin_count_one=one_recurrent_counts,
        recurrent_bin_count_two=two_recurrent_counts,
        recurrent_bin_count_one_nonpathology=one_recurrent_clean,
        recurrent_bin_count_two_nonpathology=two_recurrent_clean,
        recurrent_excess_one=one_excess,
        recurrent_excess_two=two_excess,
        recurrent_excess_one_nonpathology=one_excess_clean,
        recurrent_excess_two_nonpathology=two_excess_clean,
        global_residual_pathology=np.asarray(global_bad, dtype=np.uint8),
        raw_negative_outlier=np.asarray(raw_negative_outlier_all, dtype=np.uint8),
        raw_negative_outlier_sigma=np.asarray(raw_negative_sigma),
        raw_negative_recurrent_count=raw_negative_recurrent_counts,
    )

    recurrence = Table()
    # Preserve the original generic column for backwards compatibility while
    # making the wavelength convention explicit in new columns.
    recurrence["WAVELENGTH_ANGSTROM"] = wave
    recurrence["REST_FRAME_WAVELENGTH_ANGSTROM"] = wave
    recurrence["OBSERVED_FRAME_WAVELENGTH_ANGSTROM"] = observed_air_wave
    recurrence["OBSERVED_AIR_WAVELENGTH_ANGSTROM"] = observed_air_wave
    recurrence["OBSERVED_VACUUM_WAVELENGTH_ANGSTROM"] = observed_vacuum_wave
    recurrence["OBSERVED_SCIENCE_WAVELENGTH_ANGSTROM"] = (observed_vacuum_wave if science_medium == "vacuum" else observed_air_wave if science_medium == "air" else np.full_like(wave, np.nan))
    recurrence["N_BINS_ONE_ABS_R_GT_THRESHOLD"] = one_recurrent_counts
    recurrence["N_BINS_TWO_ABS_R_GT_THRESHOLD"] = two_recurrent_counts
    recurrence["FRACTION_BINS_ONE"] = one_recurrent_counts / float(nbin)
    recurrence["FRACTION_BINS_TWO"] = two_recurrent_counts / float(nbin)
    recurrence["BASELINE_N_BINS_ONE"] = np.full(nlog, baseline_one)
    recurrence["BASELINE_N_BINS_TWO"] = np.full(nlog, baseline_two)
    recurrence["EXCESS_N_BINS_ONE"] = one_excess
    recurrence["EXCESS_N_BINS_TWO"] = two_excess
    recurrence["N_BINS_NONPATHOLOGY"] = np.full(nlog, n_global_clean, dtype=int)
    recurrence["N_BINS_PATHOLOGY"] = np.full(nlog, n_global_bad, dtype=int)
    recurrence["N_BINS_ONE_NONPATHOLOGY"] = one_recurrent_clean
    recurrence["N_BINS_TWO_NONPATHOLOGY"] = two_recurrent_clean
    recurrence["FRACTION_ONE_NONPATHOLOGY"] = one_recurrent_clean / float(max(1, n_global_clean))
    recurrence["FRACTION_TWO_NONPATHOLOGY"] = two_recurrent_clean / float(max(1, n_global_clean))
    recurrence["BASELINE_N_BINS_ONE_NONPATHOLOGY"] = np.full(nlog, baseline_one_clean)
    recurrence["BASELINE_N_BINS_TWO_NONPATHOLOGY"] = np.full(nlog, baseline_two_clean)
    recurrence["EXCESS_N_BINS_ONE_NONPATHOLOGY"] = one_excess_clean
    recurrence["EXCESS_N_BINS_TWO_NONPATHOLOGY"] = two_excess_clean
    recurrence["N_BINS_ONE_PATHOLOGY"] = one_recurrent_bad
    recurrence["N_BINS_TWO_PATHOLOGY"] = two_recurrent_bad
    recurrence["N_BINS_RAW_NEGATIVE_OUTLIER"] = raw_negative_recurrent_counts
    recurrence["FRACTION_BINS_RAW_NEGATIVE_OUTLIER"] = raw_negative_recurrent_counts / float(nbin)
    recurrence["RANK_SCORE_LOCALIZED_EXCESS"] = rank_score
    recurrence.meta["ABS_STANDARDIZED_RESIDUAL_THRESHOLD"] = recurrent_thr
    recurrence.meta["GLOBAL_RESIDUAL_ABS_R_THRESHOLD"] = global_r_thr
    recurrence.meta["GLOBAL_RESIDUAL_PIXEL_FRACTION"] = global_r_frac
    recurrence.meta["RAW_NEGATIVE_OUTLIER_SIGMA"] = raw_negative_sigma
    recurrence.meta["N_GLOBAL_RESIDUAL_PATHOLOGY"] = n_global_bad
    recurrence.meta["N_GLOBAL_RESIDUAL_NONPATHOLOGY"] = n_global_clean
    recurrence.meta["REDSHIFT"] = None if not np.isfinite(redshift) else float(redshift)
    recurrence.meta["REST_WAVELENGTH_CONVENTION"] = str(manifest.get("template_wavelength_medium", "unknown"))
    recurrence.meta["OBSERVED_SCIENCE_WAVELENGTH_CONVENTION"] = science_medium
    recurrence.write(products_dir / "RH3_03b_recurrent_wavelengths.ecsv", format="ascii.ecsv", overwrite=True)
    order = np.argsort(np.nan_to_num(rank_score, nan=-np.inf))[::-1]
    topn = min(100, wave.size)
    recurrence_ranked = recurrence[order[:topn]]
    recurrence_ranked.write(
        products_dir / "RH3_03b_top_recurrent_wavelengths.ecsv",
        format="ascii.ecsv", overwrite=True,
    )

    candidate_mask = np.zeros(nlog, dtype=bool)
    candidate_table = Table()
    review_table = Table()
    impact_table = None
    candidate_table_path = products_dir / "RH3_03b_candidate_observed_masks.ecsv"
    review_table_path = products_dir / "RH3_03b_review_recurrent_intervals.ecsv"
    impact_table_path = products_dir / "RH3_03b_candidate_mask_impact.ecsv"
    config_snippet_path = products_dir / "RH3_03b_candidate_mask_config_snippet.txt"

    if args.candidate_mask_mode != "off":
        if not np.isfinite(redshift):
            logger.warning(
                "Candidate observed-frame mask generation skipped because redshift is unavailable. "
                "Supply --config or --redshift and rerun 03b."
            )
        elif science_medium not in {"air", "vacuum"}:
            logger.warning(
                "Candidate observed-frame mask generation skipped because science wavelength medium is %r.",
                science_medium,
            )
        else:
            candidate_mask, candidate_table, review_table = _candidate_mask_products(
                wave=wave,
                raw_negative_counts=raw_negative_recurrent_counts,
                rank_score=rank_score,
                nbin=nbin,
                n_nonpathology=n_global_clean,
                redshift=redshift,
                template_medium=template_medium,
                science_medium=science_medium,
                mode=args.candidate_mask_mode,
                raw_fraction=float(args.candidate_raw_negative_fraction),
                residual_excess_fraction=float(args.candidate_residual_excess_fraction),
                padding_pixels=int(args.candidate_mask_padding_pixels),
                bridge_pixels=int(args.candidate_mask_bridge_pixels),
            )
            candidate_table.write(candidate_table_path, format="ascii.ecsv", overwrite=True)
            review_table.write(review_table_path, format="ascii.ecsv", overwrite=True)
            _write_candidate_config_snippet(config_snippet_path, candidate_table, science_medium)
            frac_masked_grid = float(np.mean(candidate_mask))
            logger.info(
                "Candidate fixed-wavelength mask (%s): %d interval(s), %d/%d log pixels (%.2f%% of fit grid).",
                args.candidate_mask_mode, len(candidate_table), int(np.count_nonzero(candidate_mask)), nlog, 100.0 * frac_masked_grid,
            )
            logger.info(
                "Candidate config values are written in observed %s wavelengths because Script 3 applies "
                "RH3_MASK_OBSERVED_RANGES_ANGSTROM directly to the native science wavelength array.",
                science_medium,
            )
            if frac_masked_grid > float(args.candidate_mask_max_fit_fraction):
                logger.warning(
                    "CANDIDATE_MASK_TOO_BROAD | proposed mask removes %.1f%% of the fit grid (>%.1f%% warning threshold).",
                    100.0 * frac_masked_grid, 100.0 * float(args.candidate_mask_max_fit_fraction),
                )
            if len(review_table):
                logger.info(
                    "Review-only recurrent-residual intervals: %d. These are NOT part of the conservative "
                    "raw-negative candidate unless --candidate-mask-mode extended is requested.",
                    len(review_table),
                )

            impact_table = _candidate_mask_impact_table(
                spectra, candidate_mask, global_r_thr, global_r_frac, raw_negative_sigma
            )
            impact_table.meta["NO_PPXF_REFIT"] = True
            impact_table.meta["CANDIDATE_MASK_MODE"] = args.candidate_mask_mode
            impact_table.write(impact_table_path, format="ascii.ecsv", overwrite=True)
            _plot_candidate_mask_impact(impact_table, figures_dir / "RH3_03b_candidate_mask_impact.png")

            gb0 = np.asarray(impact_table["GLOBAL_RESIDUAL_PATHOLOGY_BEFORE"], dtype=bool)
            gb1 = np.asarray(impact_table["GLOBAL_RESIDUAL_PATHOLOGY_AFTER"], dtype=bool)
            lev0 = np.asarray(impact_table["TWO_TOP1_CHI2_FRAC_BEFORE"], dtype=float)
            lev1 = np.asarray(impact_table["TWO_TOP1_CHI2_FRAC_AFTER"], dtype=float)
            rr0 = np.asarray(impact_table["TWO_ROBUST_RESIDUAL_TO_NOISE_BEFORE"], dtype=float)
            rr1 = np.asarray(impact_table["TWO_ROBUST_RESIDUAL_TO_NOISE_AFTER"], dtype=float)
            logger.info(
                "Candidate-mask SAME-MODEL impact (no pPXF refit): global-pathology bins %d -> %d; "
                "median top1 chi2 fraction %.3g -> %.3g; median robust residual/noise %.3g -> %.3g.",
                int(np.count_nonzero(gb0)), int(np.count_nonzero(gb1)),
                _percentile_or_nan(lev0, 50), _percentile_or_nan(lev1, 50),
                _percentile_or_nan(rr0, 50), _percentile_or_nan(rr1, 50),
            )

            impact_bins = _parse_bin_id_list(args.mask_impact_bins, nbin)
            for bid in impact_bins:
                _plot_mask_impact_bin(
                    bid, spectra, wave, candidate_mask,
                    figures_dir / f"RH3_03b_candidate_mask_bin_{bid:04d}.png",
                )

    # Compact residual/noise table and two-row group summary for quick upload/review.
    residual_cols = [
        "BIN_ID", "GLOBAL_RESIDUAL_PATHOLOGY", "TWO_GLOBAL_ABS_R_FRACTION",
        "NOISE_MIN", "NOISE_P01", "NOISE_MEDIAN", "NOISE_P99", "NOISE_MAX", "NOISE_NONPOSITIVE_FRAC",
        "GALAXY_MIN", "GALAXY_P01", "GALAXY_MEDIAN", "GALAXY_P99", "GALAXY_MAX", "GALAXY_NEGATIVE_FRAC",
        "GALAXY_ROBUST_SCATTER", "GALAXY_MIN_ROBUST_Z", "GALAXY_N_NEGATIVE_OUTLIERS", "GALAXY_NEGATIVE_OUTLIER_FRAC",
        "TWO_RESIDUAL_MEDIAN", "TWO_RESIDUAL_MAD_SIGMA", "TWO_MEDIAN_ABS_RESIDUAL", "TWO_MEDIAN_ABS_R",
        "TWO_ROBUST_RESIDUAL_TO_NOISE", "TWO_MODEL_MEDIAN_OFFSET_NOISE",
        "TWO_MODEL_MIN", "TWO_MODEL_P01", "TWO_MODEL_MEDIAN", "TWO_MODEL_P99", "TWO_MODEL_MAX",
        "ADDITIVE_POLY_MAX_ABS_COEFF", "ADDITIVE_POLY_L2",
        "RESIDUAL_SCALE_WARNING", "MODEL_MEDIAN_OFFSET_WARNING",
    ]
    Table(metrics[residual_cols]).write(
        products_dir / "RH3_03b_residual_scale_diagnostics.ecsv", format="ascii.ecsv", overwrite=True
    )
    group_rows = []
    for label, mask in (("GLOBAL_RESIDUAL_PATHOLOGY", global_bad), ("NONPATHOLOGY", ~global_bad)):
        if not np.any(mask):
            continue
        group_rows.append((
            label, int(np.count_nonzero(mask)),
            _percentile_or_nan(np.asarray(metrics["NOISE_MEDIAN"], dtype=float)[mask], 50),
            _percentile_or_nan(np.asarray(metrics["TWO_RESIDUAL_MAD_SIGMA"], dtype=float)[mask], 50),
            _percentile_or_nan(np.asarray(metrics["TWO_ROBUST_RESIDUAL_TO_NOISE"], dtype=float)[mask], 50),
            _percentile_or_nan(np.asarray(metrics["TWO_MEDIAN_ABS_R"], dtype=float)[mask], 50),
            _percentile_or_nan(np.asarray(metrics["GALAXY_MIN"], dtype=float)[mask], 50),
            _percentile_or_nan(np.asarray(metrics["TWO_MODEL_MEDIAN"], dtype=float)[mask], 50),
            _percentile_or_nan(np.asarray(metrics["ADDITIVE_POLY_MAX_ABS_COEFF"], dtype=float)[mask], 50),
        ))
    group_table = Table(
        rows=group_rows,
        names=(
            "GROUP", "N_BINS", "MEDIAN_NOISE", "MEDIAN_RESIDUAL_MAD_SIGMA",
            "MEDIAN_ROBUST_RESIDUAL_TO_NOISE", "MEDIAN_ABS_R", "MEDIAN_GALAXY_MIN",
            "MEDIAN_TWO_MODEL", "MEDIAN_ADDITIVE_POLY_MAX_ABS_COEFF",
        ),
    )
    group_table.write(products_dir / "RH3_03b_residual_group_summary.ecsv", format="ascii.ecsv", overwrite=True)

    representatives = _select_representative_bins(metrics, args.representative_count)
    rep_table = Table(
        rows=[(label, bid) for label, bid in representatives],
        names=("CATEGORY", "BIN_ID"),
        dtype=("U64", int),
    )
    rep_table.write(products_dir / "RH3_03b_representative_bins.ecsv", format="ascii.ecsv", overwrite=True)
    logger.info("Automatically selected representative bins: %s", representatives)

    _section(logger, "5. Write cheap validation figures")
    _plot_hist(metrics["NEFF_IPR"], figures_dir / "RH3_03b_neff_hist.png", "RH3 likelihood information content", "N_eff (IPR)", logx=True)
    _plot_hist(metrics["SWAP_REL_P95"], figures_dir / "RH3_03b_swap_relative_error_hist.png", "A/B swap-symmetry relative error", "Per-bin swap p95 / robust Delta-chi2 scale", logx=True)
    _plot_hist(metrics["TWO_TOP1_CHI2_FRAC"], figures_dir / "RH3_03b_top1_leverage_hist.png", "Two-component worst-pixel chi2 leverage", "Fraction of chi2 from worst pixel")
    _plot_hist(metrics["SIGMA_BOUNDARY_STATE_FRAC"], figures_dir / "RH3_03b_sigma_boundary_hist.png", "Sigma-boundary state fraction", "Fraction of grid states near sigma bound")
    _plot_hist(metrics["N_EPSILON_MIN_PLATEAUS"], figures_dir / "RH3_03b_topology_plateau_count_hist.png", "Epsilon-thick local-minimum plateaus", "Number of plateaus")
    _plot_hist(metrics["NOISE_MEDIAN"], figures_dir / "RH3_03b_formal_noise_median_hist.png", "Median formal noise per RH3 bin", "median normalized formal noise", logx=True)
    _plot_hist(metrics["TWO_ROBUST_RESIDUAL_TO_NOISE"], figures_dir / "RH3_03b_robust_residual_to_noise_hist.png", "Robust residual scatter / formal noise", "1.4826 MAD(data-model) / median(noise)", logx=True)
    _plot_hist(metrics["TWO_GLOBAL_ABS_R_FRACTION"], figures_dir / "RH3_03b_binwide_residual_fraction_hist.png", "Bin-wide standardized-residual pathology", f"fraction of valid pixels with |r| > {global_r_thr:g}")
    _plot_hist(metrics["GALAXY_N_NEGATIVE_OUTLIERS"], figures_dir / "RH3_03b_raw_negative_outlier_count_hist.png", "Model-independent raw negative outliers per bin", f"N pixels below median - {raw_negative_sigma:g} robust-MAD sigma")

    fig, axes = plt.subplots(2, 1, figsize=(11, 8), sharex=True)
    axes[0].plot(wave, one_recurrent_counts, lw=1.0, label="1-component")
    axes[0].plot(wave, two_recurrent_counts, lw=1.0, label="local-best 2-component")
    axes[0].axhline(baseline_one, ls="--", lw=0.8, label="1C wavelength median")
    axes[0].axhline(baseline_two, ls=":", lw=0.8, label="2C wavelength median")
    axes[0].set_ylabel(f"N bins with |r| > {recurrent_thr:g}")
    axes[0].set_title("Recurrent standardized residuals: absolute count and localized excess")
    axes[0].legend(fontsize=8, ncol=2)
    axes[1].plot(wave, one_excess_clean, lw=1.0, label="1C non-pathology excess")
    axes[1].plot(wave, two_excess_clean, lw=1.0, label="2C non-pathology excess")
    axes[1].axhline(0.0, lw=0.8)
    axes[1].set_xlabel("Rest-frame wavelength (Angstrom)")
    axes[1].set_ylabel("Excess affected bins above wavelength median")
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(figures_dir / "RH3_03b_recurrent_residual_wavelengths.png", dpi=160)
    plt.close(fig)

    if np.isfinite(redshift):
        fig, axes = plt.subplots(2, 1, figsize=(11, 8), sharex=True)
        axes[0].plot(observed_wave, one_recurrent_counts, lw=1.0, label="1-component")
        axes[0].plot(observed_wave, two_recurrent_counts, lw=1.0, label="local-best 2-component")
        axes[0].axhline(baseline_one, ls="--", lw=0.8)
        axes[0].axhline(baseline_two, ls=":", lw=0.8)
        axes[0].set_ylabel(f"N bins with |r| > {recurrent_thr:g}")
        axes[0].set_title(f"Observed-frame recurrent residuals (z={redshift:.6f})")
        axes[0].legend(fontsize=8)
        axes[1].plot(observed_wave, one_excess_clean, lw=1.0, label="1C non-pathology excess")
        axes[1].plot(observed_wave, two_excess_clean, lw=1.0, label="2C non-pathology excess")
        axes[1].axhline(0.0, lw=0.8)
        axes[1].set_xlabel("Observed-frame air wavelength (Angstrom)")
        axes[1].set_ylabel("Excess affected bins")
        axes[1].legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(figures_dir / "RH3_03b_recurrent_residual_wavelengths_observed.png", dpi=160)
        plt.close(fig)

    fig, ax = plt.subplots(figsize=(11, 5))
    ax.plot(wave, raw_negative_recurrent_counts, lw=1.0)
    ax.set_xlabel("Rest-frame wavelength (Angstrom)")
    ax.set_ylabel("N bins with raw negative outlier")
    ax.set_title(f"Model-independent raw negative outliers (< median - {raw_negative_sigma:g} robust-MAD sigma)")
    fig.tight_layout()
    fig.savefig(figures_dir / "RH3_03b_raw_negative_outlier_wavelengths.png", dpi=160)
    plt.close(fig)
    if np.isfinite(redshift) and np.all(np.isfinite(observed_air_wave)):
        fig, ax = plt.subplots(figsize=(11, 5))
        ax.plot(observed_air_wave, raw_negative_recurrent_counts, lw=1.0)
        ax.set_xlabel("Observed-frame air wavelength (Angstrom)")
        ax.set_ylabel("N bins with raw negative outlier")
        ax.set_title(f"Observed-frame model-independent raw negative outliers (z={redshift:.6f})")
        fig.tight_layout()
        fig.savefig(figures_dir / "RH3_03b_raw_negative_outlier_wavelengths_observed.png", dpi=160)
        plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 6))
    x = np.asarray(metrics["NOISE_MEDIAN"], dtype=float)
    y = np.asarray(metrics["TWO_RESIDUAL_MAD_SIGMA"], dtype=float)
    use = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
    ax.scatter(x[use & (~global_bad)], y[use & (~global_bad)], s=18, alpha=0.7, label="non-pathology")
    ax.scatter(x[use & global_bad], y[use & global_bad], s=18, alpha=0.7, label="global residual pathology")
    if np.any(use):
        lo = min(np.min(x[use]), np.min(y[use]))
        hi = max(np.max(x[use]), np.max(y[use]))
        ax.plot([lo, hi], [lo, hi], ls="--", lw=0.8, label="residual scatter = formal noise")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Median formal noise")
    ax.set_ylabel("Robust 2C residual scatter (1.4826 MAD)")
    ax.set_title("Formal noise versus realized residual scale")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(figures_dir / "RH3_03b_noise_vs_residual_scatter.png", dpi=160)
    plt.close(fig)

    if bin_map is not None:
        spatial_specs = (
            ("NEFF_IPR", "RH3_03b_neff_map.png", "RH3 likelihood effective number of states", "N_eff (IPR)"),
            ("SWAP_REL_P95", "RH3_03b_swap_relative_error_map.png", "A/B swap-symmetry relative error", "swap p95 / Delta-chi2 scale"),
            ("TWO_TOP1_CHI2_FRAC", "RH3_03b_top1_leverage_map.png", "Worst-pixel leverage in local-best 2C fit", "fraction of chi2 from worst pixel"),
            ("TWO_MAX_ABS_R", "RH3_03b_max_standardized_residual_map.png", "Largest local-best 2C standardized residual", "max |(data-model)/noise|"),
            ("VELOCITY_EDGE_WEIGHT_MASS", "RH3_03b_velocity_edge_mass_map.png", "Likelihood mass on V-grid boundary", "weight mass"),
            ("FA_EDGE_WEIGHT_MASS", "RH3_03b_fA_edge_mass_map.png", "Likelihood mass on f_A-grid boundary", "weight mass"),
            ("N_EPSILON_MIN_PLATEAUS", "RH3_03b_topology_plateau_count_map.png", "Epsilon-thick local-minimum plateau count", "number of plateaus"),
            ("SIGMA_BOUNDARY_STATE_FRAC", "RH3_03b_sigma_boundary_fraction_map.png", "Sigma-boundary grid-state fraction", "state fraction"),
            ("TWO_GLOBAL_ABS_R_FRACTION", "RH3_03b_binwide_residual_fraction_map.png", "Bin-wide standardized-residual pathology", f"fraction pixels |r|>{global_r_thr:g}"),
            ("TWO_ROBUST_RESIDUAL_TO_NOISE", "RH3_03b_robust_residual_to_noise_map.png", "Robust residual scatter / formal noise", "MAD residual / median noise"),
            ("NOISE_MEDIAN", "RH3_03b_formal_noise_median_map.png", "Median formal noise", "normalized formal noise"),
            ("ADDITIVE_POLY_MAX_ABS_COEFF", "RH3_03b_additive_poly_amplitude_map.png", "Maximum absolute additive-polynomial coefficient", "max |poly coefficient|"),
            ("GALAXY_N_NEGATIVE_OUTLIERS", "RH3_03b_raw_negative_outlier_count_map.png", "Model-independent raw negative-outlier count", "N raw negative outliers"),
        )
        for col, filename, title, label in spatial_specs:
            _plot_bin_map(bin_map, np.asarray(metrics[col], dtype=float), figures_dir / filename, title, label)
    else:
        logger.warning("Spatial maps skipped because a usable Script-2 master bin map was not available.")

    for category, bid in representatives:
        _plot_representative_bin(
            bid, cube, spectra, metrics,
            rep_dir / f"{category}__bin_{bid:04d}.png",
            residual_threshold=recurrent_thr,
        )

    _section(logger, "6. Validation contract summary")
    mechanical_status = "PASS" if chi2_recompute_pass and np.all(np.asarray(metrics["N_FAILED_STATES"]) == 0) else "WARN"
    leverage_status = "WARN" if leverage_frac > 0 else "PASS"
    residual_scale_status = "WARN" if global_bad_frac > 0 else "PASS"
    boundary_status = "WARN" if (sigma_bound_warn_frac > 0 or velocity_edge_warn_frac > 0 or fa_edge_warn_frac > 0) else "PASS"
    overall = "PASS"
    if mechanical_status != "PASS" or swap_status == "FAIL":
        overall = "FAIL"
    elif swap_status == "WARN" or leverage_status == "WARN" or residual_scale_status == "WARN" or boundary_status == "WARN":
        overall = "WARN"

    logger.info("Mechanical integrity: %s", mechanical_status)
    logger.info("A/B swap symmetry: %s", swap_status)
    logger.info("Spectral leverage: %s", leverage_status)
    logger.info("Bin-wide residual scale/pathology: %s", residual_scale_status)
    logger.info("Grid/sigma boundaries: %s", boundary_status)
    logger.info("SCRIPT 03b OVERALL: %s", overall)
    logger.info("Important: WARN/FAIL here diagnoses the DEVELOPMENT RUN; it does not force the weak RL test data to resemble a successful CRD decomposition.")

    manifest_out = {
        "script": "03b_validate_RH3_likelihood_cubes",
        "schema_version": 2,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source_script3_run": str(script3_run),
        "source_script2_run": str(resolved_script2) if resolved_script2 is not None else manifest.get("source_script2_run"),
        "source_products": {
            "likelihood_cubes": str(cube_path),
            "selected_fits": str(spectra_path),
            "xsl_templates": str(template_path) if template_path.is_file() else None,
        },
        "redshift": None if not np.isfinite(redshift) else float(redshift),
        "shape": {"n_bins": nbin, "n_VA": nva, "n_VB": nvb, "n_fA": nfa, "states_per_bin": nva*nvb*nfa, "n_log_pixels": nlog},
        "validation_contract": {
            "mechanical_integrity": mechanical_status,
            "swap_symmetry": swap_status,
            "spectral_leverage": leverage_status,
            "residual_scale": residual_scale_status,
            "boundaries": boundary_status,
            "overall": overall,
        },
        "thresholds": {
            "swap_near_min_delta_chi2_diagnostic": SWAP_NEAR_MIN_DELTA_CHI2,
            "swap_relative_pass": SWAP_REL_PASS,
            "swap_relative_warn": SWAP_REL_WARN,
            "topology_epsilon_mode": "explicit" if explicit_topo is not None else "auto_from_low_chi_swap_mismatch",
            "topology_epsilon_explicit": explicit_topo,
            "topology_low_chi_rank_fraction": TOPOLOGY_LOW_CHI_FRACTION,
            "residual_thresholds": list(RESIDUAL_THRESHOLDS),
            "recurrent_abs_residual_threshold": recurrent_thr,
            "global_residual_abs_r_threshold": global_r_thr,
            "global_residual_pixel_fraction": global_r_frac,
            "robust_residual_to_noise_warn": ROBUST_RESIDUAL_TO_NOISE_WARN,
            "model_median_offset_noise_warn": MODEL_MEDIAN_OFFSET_NOISE_WARN,
            "raw_negative_outlier_sigma": raw_negative_sigma,
            "candidate_mask_mode": args.candidate_mask_mode,
            "candidate_raw_negative_fraction": float(args.candidate_raw_negative_fraction),
            "candidate_residual_excess_fraction": float(args.candidate_residual_excess_fraction),
            "candidate_mask_padding_pixels": int(args.candidate_mask_padding_pixels),
            "candidate_mask_bridge_pixels": int(args.candidate_mask_bridge_pixels),
            "leverage_top1_warn": LEVERAGE_TOP1_WARN,
            "leverage_top10_warn": LEVERAGE_TOP10_WARN,
            "leverage_max_abs_r_warn": LEVERAGE_MAX_R_WARN,
            "grid_edge_weight_warn": GRID_EDGE_WEIGHT_WARN,
            "sigma_boundary_state_warn": SIGMA_BOUNDARY_STATE_WARN,
        },
        "global_metrics": {
            "swap_rel_p95_across_bins_p95": swap_rel_p95_global,
            "median_neff_ipr": _percentile_or_nan(np.asarray(metrics["NEFF_IPR"], dtype=float), 50),
            "median_n_mass_95": float(np.nanmedian(np.asarray(metrics["N_MASS_95"], dtype=float))),
            "fraction_bins_leverage_warning": leverage_frac,
            "n_bins_global_residual_pathology": n_global_bad,
            "fraction_bins_global_residual_pathology": global_bad_frac,
            "n_bins_with_raw_negative_outliers": int(np.count_nonzero(raw_negative_bin)),
            "fraction_bins_with_raw_negative_outliers": float(np.mean(raw_negative_bin)),
            "n_bins_raw_negative_and_global_residual_pathology": n_overlap_raw_global,
            "fraction_global_residual_pathology_with_raw_negative_outlier": frac_global_with_raw,
            "fraction_raw_negative_outlier_bins_with_global_residual_pathology": frac_raw_with_global,
            "median_recurrent_count_one_all_wavelengths": baseline_one,
            "median_recurrent_count_two_all_wavelengths": baseline_two,
            "median_recurrent_count_one_nonpathology": baseline_one_clean,
            "median_recurrent_count_two_nonpathology": baseline_two_clean,
            "median_two_robust_residual_to_noise": _percentile_or_nan(np.asarray(metrics["TWO_ROBUST_RESIDUAL_TO_NOISE"], dtype=float), 50),
            "median_two_median_abs_r": _percentile_or_nan(np.asarray(metrics["TWO_MEDIAN_ABS_R"], dtype=float), 50),
            "fraction_bins_sigma_boundary_state_fraction_gt_0p5": sigma_bound_warn_frac,
            "fraction_bins_velocity_edge_weight_mass_gt_0p25": velocity_edge_warn_frac,
            "fraction_bins_fA_edge_weight_mass_gt_0p25": fa_edge_warn_frac,
            "max_one_component_chi2_recompute_relative_error": max_one_rel,
            "max_two_component_chi2_recompute_relative_error": max_two_rel,
        },
        "template_validation": template_result,
        "fraction_bookkeeping": fraction_check,
        "representative_bins": [{"category": cat, "bin_id": int(bid)} for cat, bid in representatives],
        "products": {
            "validation_summary_ecsv": str(metrics_path),
            "raw_diagnostics_npz": str(products_dir / "RH3_03b_raw_diagnostics.npz"),
            "recurrent_wavelengths_ecsv": str(products_dir / "RH3_03b_recurrent_wavelengths.ecsv"),
            "top_recurrent_wavelengths_ecsv": str(products_dir / "RH3_03b_top_recurrent_wavelengths.ecsv"),
            "candidate_observed_masks_ecsv": str(candidate_table_path) if candidate_table_path.is_file() else None,
            "review_recurrent_intervals_ecsv": str(review_table_path) if review_table_path.is_file() else None,
            "candidate_mask_impact_ecsv": str(impact_table_path) if impact_table_path.is_file() else None,
            "candidate_mask_config_snippet_txt": str(config_snippet_path) if config_snippet_path.is_file() else None,
            "residual_scale_diagnostics_ecsv": str(products_dir / "RH3_03b_residual_scale_diagnostics.ecsv"),
            "residual_group_summary_ecsv": str(products_dir / "RH3_03b_residual_group_summary.ecsv"),
            "representative_bins_ecsv": str(products_dir / "RH3_03b_representative_bins.ecsv"),
            "figures_dir": str(figures_dir),
        },
        "scientific_caveat": (
            "Script-3 likelihood widths remain non-publication-final while the spectral variance/covariance "
            "model is uncalibrated. Script 03b information-content/topology metrics are development diagnostics."
        ),
    }
    _write_json(manifest_out, output_dir / "metadata" / "script03b_validation_manifest.json")
    logger.info("Per-bin validation table: %s", metrics_path)
    logger.info("Raw residual/swap diagnostics: %s", products_dir / "RH3_03b_raw_diagnostics.npz")
    logger.info("Residual/noise diagnostics: %s", products_dir / "RH3_03b_residual_scale_diagnostics.ecsv")
    logger.info("Residual group summary: %s", products_dir / "RH3_03b_residual_group_summary.ecsv")
    logger.info("Recurrent wavelengths (rest + observed when z available): %s", products_dir / "RH3_03b_recurrent_wavelengths.ecsv")
    if candidate_table_path.is_file():
        logger.info("Candidate observed-frame masks: %s", candidate_table_path)
        logger.info("Candidate config snippet: %s", config_snippet_path)
        logger.info("Candidate-mask same-model impact: %s", impact_table_path)
    logger.info("Representative-bin figures: %s", rep_dir)
    logger.info("Validation manifest: %s", output_dir / "metadata" / "script03b_validation_manifest.json")
    _section(logger, "SCRIPT 03b COMPLETE", "=")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
