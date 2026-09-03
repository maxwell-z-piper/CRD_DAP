"""Script-3 covariance calibration and model-selection helpers.

The production philosophy is intentionally conservative:

1. obtain high-quality free two-component pPXF residuals for every PowerBin;
2. estimate the residual variance scale and spectral correlation anew for the
   current run;
3. calibrate four increasingly flexible covariance candidates;
4. require each candidate to converge and test whether it whitens the residuals;
5. evaluate the *full production likelihood grid* for a deterministic set of
   representative PowerBins under all four candidates; and
6. adopt the least-complex model that passes residual QC and is scientifically
   stable against the more complex passing candidates.

The covariance is then frozen for the full Script-3 grid.  It is never learned
separately for individual (V_A,V_B,f_A) states.
"""

from __future__ import annotations

from dataclasses import dataclass
import inspect
from pathlib import Path
from typing import Any, Iterable

import numpy as np
from astropy.table import Table

from . import covariance
from . import ppxf_utils


MODEL_ORDER = ("M1", "M2", "M3", "M4")
MODEL_DESCRIPTION = {
    "M1": "diagonal formal errors with per-bin residual scale s_i",
    "M2": "one common wavelength-stationary R plus per-bin scale s_i",
    "M3": "common wavelength-block-dependent R plus per-bin scale s_i",
    "M4": "per-bin wavelength-block-dependent R_i plus per-bin scale s_i",
}


@dataclass
class FreeTwoComponentFit:
    """Best multi-start free two-component calibration fit for one PowerBin."""

    success: bool
    chi2_total: float
    reduced_chi2: float
    velocity: np.ndarray
    sigma: np.ndarray
    fraction_a: float
    bestfit: np.ndarray
    weights: np.ndarray
    start_index: int
    n_successful_starts: int
    error_message: str = ""


@dataclass
class CandidateModel:
    """Compact covariance candidate representation.

    ``rho`` has shape ``(1,nlag)`` for M1/M2, ``(nblock,nlag)`` for M3, and
    ``(nbin,nblock,nlag)`` for M4.  The dense covariance/whitener is deliberately
    *not* stored: it is reconstructed and factorized once per PowerBin when
    needed, which keeps products small while preserving exact reproducibility.
    """

    name: str
    scale: np.ndarray
    rho: np.ndarray
    block_index: np.ndarray
    converged: bool
    n_iterations: int
    final_scale_change: float
    final_rho_change: float
    requirement_a_pass: bool = False
    whitened_scale_median: float = np.nan
    whitened_band_center: np.ndarray | None = None
    whitened_band_lower: np.ndarray | None = None
    whitened_band_upper: np.ndarray | None = None
    max_lag: int = 0
    numerical_regularization_count: int = 0

    def rho_for_bin(self, bin_id: int) -> np.ndarray:
        if self.name in {"M1", "M2", "M3"}:
            return np.asarray(self.rho, dtype=float)
        if self.name == "M4":
            return np.asarray(self.rho[int(bin_id)], dtype=float)
        raise ValueError(f"Unknown covariance model {self.name!r}")


@dataclass(frozen=True)
class ModelComparison:
    """Requirement-B comparison of two full representative likelihood cubes."""

    agree: bool
    max_minimum_cell_shift: int
    max_interval_edge_cell_shift: int
    failing_bins: tuple[int, ...]
    details: tuple[dict[str, Any], ...]


def require_cached_ppxf_patch(expected_version: str = "9.4.8") -> None:
    """Hard-fail unless the validated CRD_DAP cached-whitener pPXF patch exists."""
    import ppxf
    from ppxf.ppxf import ppxf as PPXF

    version = str(getattr(ppxf, "__version__", "unknown"))
    if version != str(expected_version):
        raise RuntimeError(
            f"CRD_DAP covariance-aware Script 3 requires validated pPXF {expected_version}; found {version}."
        )
    if "noise_inv_cholesky" not in inspect.signature(PPXF).parameters:
        raise RuntimeError(
            "CRD_DAP_PPXF_CACHED_WHITENER_PATCH_MISSING | The active pPXF installation "
            "does not expose noise_inv_cholesky. Apply the repository's ppxf_patch_9_4_8 "
            "inside the CRD_DAP conda environment and run its regression tests before Script 3."
        )


def validate_covariance_config(cfg) -> None:
    """Validate Script-3 covariance settings close to the machinery that uses them.

    Keeping these checks in this module avoids requiring a wholesale replacement
    of ``crd_utils/config.py`` when the covariance feature is added to an existing
    CRD_DAP checkout.  The target config still records every numerical choice and
    is snapshotted in the normal run provenance.
    """
    if not isinstance(getattr(cfg, "RH3_COVARIANCE_ENABLE", True), (bool, np.bool_)):
        raise ValueError("RH3_COVARIANCE_ENABLE must be boolean")
    if not bool(getattr(cfg, "RH3_COVARIANCE_ENABLE", True)):
        raise RuntimeError(
            "RH3_COVARIANCE_DISABLED | Production Script 3 now requires the residual-based "
            "covariance calibration before likelihood cubes are evaluated."
        )
    if int(getattr(cfg, "RH3_COVARIANCE_MAX_LAG", 20)) < 1:
        raise ValueError("RH3_COVARIANCE_MAX_LAG must be >= 1")
    if int(getattr(cfg, "RH3_COVARIANCE_MIN_PAIRS", 25)) < 3:
        raise ValueError("RH3_COVARIANCE_MIN_PAIRS must be >= 3")
    if int(getattr(cfg, "RH3_COVARIANCE_BOOTSTRAP_N", 2000)) < 100:
        raise ValueError("RH3_COVARIANCE_BOOTSTRAP_N must be >= 100")
    ci = float(getattr(cfg, "RH3_COVARIANCE_BOOTSTRAP_CONFIDENCE", 0.95))
    if not 0.5 < ci < 1.0:
        raise ValueError("RH3_COVARIANCE_BOOTSTRAP_CONFIDENCE must lie in (0.5, 1)")
    if int(getattr(cfg, "RH3_COVARIANCE_WAVELENGTH_BLOCKS", 3)) < 2:
        raise ValueError("RH3_COVARIANCE_WAVELENGTH_BLOCKS must be >= 2")
    tol = float(getattr(cfg, "RH3_COVARIANCE_CONVERGENCE_TOL", 0.01))
    if not np.isfinite(tol) or tol <= 0:
        raise ValueError("RH3_COVARIANCE_CONVERGENCE_TOL must be finite and positive")
    if int(getattr(cfg, "RH3_COVARIANCE_MAX_ITER", 5)) < 1:
        raise ValueError("RH3_COVARIANCE_MAX_ITER must be >= 1")
    white_tol = float(getattr(cfg, "RH3_COVARIANCE_WHITENED_SCALE_TOL", 0.05))
    if not 0 < white_tol < 1:
        raise ValueError("RH3_COVARIANCE_WHITENED_SCALE_TOL must lie in (0, 1)")
    sep = tuple(float(x) for x in getattr(
        cfg, "RH3_COVARIANCE_CALIBRATION_SEPARATION_FRACTIONS", (0.15, 0.40, 0.75)
    ))
    if len(sep) < 2 or any((not np.isfinite(x)) or x <= 0 or x > 1 for x in sep):
        raise ValueError(
            "RH3_COVARIANCE_CALIBRATION_SEPARATION_FRACTIONS requires at least two finite values in (0, 1]"
        )
    nrep = int(getattr(cfg, "RH3_COVARIANCE_VALIDATION_RADIAL_BINS", 12))
    if nrep < 4 or nrep % 2:
        raise ValueError("RH3_COVARIANCE_VALIDATION_RADIAL_BINS must be an even integer >= 4")
    if float(getattr(cfg, "RH3_COVARIANCE_PA_CORRIDOR_DIAMETER_FACTOR", 1.0)) <= 0:
        raise ValueError("RH3_COVARIANCE_PA_CORRIDOR_DIAMETER_FACTOR must be positive")
    inner = float(getattr(cfg, "RH3_COVARIANCE_2SIGMA_INNER_RADIUS_FRACTION", 0.10))
    outer = float(getattr(cfg, "RH3_COVARIANCE_2SIGMA_OUTER_RADIUS_FRACTION", 0.95))
    if not 0 <= inner < outer <= 1:
        raise ValueError("RH3_COVARIANCE_2SIGMA radius fractions must satisfy 0 <= inner < outer <= 1")
    if float(getattr(cfg, "RH3_COVARIANCE_EIGEN_FLOOR", 1.0e-8)) <= 0:
        raise ValueError("RH3_COVARIANCE_EIGEN_FLOOR must be positive")
    if int(getattr(cfg, "RH3_COVARIANCE_MODEL_AGREEMENT_MAX_CELL_SHIFT", 1)) < 0:
        raise ValueError("RH3_COVARIANCE_MODEL_AGREEMENT_MAX_CELL_SHIFT cannot be negative")
    dchi = tuple(float(x) for x in getattr(
        cfg, "RH3_COVARIANCE_MODEL_AGREEMENT_DELTA_CHI2", (1.0, 4.0)
    ))
    if not dchi or any((not np.isfinite(x)) or x <= 0 for x in dchi):
        raise ValueError("RH3_COVARIANCE_MODEL_AGREEMENT_DELTA_CHI2 must contain positive finite thresholds")


def calibration_velocity_seeds(
    single_velocity: float,
    *,
    va_bounds: tuple[float, float],
    vb_bounds: tuple[float, float],
    separation_fractions: Iterable[float],
) -> list[tuple[float, float]]:
    """Construct deterministic multi-start velocity seeds around a 1-C solution.

    Both component orderings are included for every separation.  This is not a
    prior on the physical disk velocities; it is only a robust optimization
    strategy for obtaining a good model spectrum whose residuals can be used to
    calibrate the noise covariance.
    """
    va_lo, va_hi = map(float, va_bounds)
    vb_lo, vb_hi = map(float, vb_bounds)
    lo = max(va_lo, vb_lo)
    hi = min(va_hi, vb_hi)
    if hi <= lo:
        raise ValueError("VA/VB calibration velocity bounds do not overlap")
    center = float(single_velocity) if np.isfinite(single_velocity) else 0.0
    center = float(np.clip(center, lo, hi))
    span = float(hi - lo)
    seeds: list[tuple[float, float]] = []
    for frac in separation_fractions:
        frac = float(frac)
        if not 0 < frac <= 1:
            raise ValueError("Calibration separation fractions must lie in (0,1]")
        sep = frac * span
        a = float(np.clip(center - 0.5 * sep, va_lo, va_hi))
        b = float(np.clip(center + 0.5 * sep, vb_lo, vb_hi))
        seeds.extend([(a, b), (b, a)])
    # Stable deduplication after clipping near bounds.
    unique: list[tuple[float, float]] = []
    for seed in seeds:
        if not any(np.allclose(seed, old, atol=1e-9, rtol=0.0) for old in unique):
            unique.append(seed)
    if len(unique) < 2:
        raise RuntimeError("Could not construct at least two distinct two-component calibration starts")
    return unique


def fit_free_two_component_multistart(
    *,
    templates_two_component: np.ndarray,
    component: np.ndarray,
    galaxy: np.ndarray,
    noise: np.ndarray,
    velscale: float,
    lam: np.ndarray,
    lam_temp: np.ndarray,
    goodpixels: np.ndarray,
    single_velocity: float,
    start_sigma: float,
    va_bounds: tuple[float, float],
    vb_bounds: tuple[float, float],
    sigma_bounds: tuple[float, float],
    degree: int,
    mdegree: int,
    regul: float,
    separation_fractions: Iterable[float],
    noise_inv_cholesky: np.ndarray | None = None,
) -> FreeTwoComponentFit:
    """Run robust multi-start free two-component pPXF and keep the best model."""
    seeds = calibration_velocity_seeds(
        single_velocity,
        va_bounds=va_bounds,
        vb_bounds=vb_bounds,
        separation_fractions=separation_fractions,
    )
    best = None
    successes = 0
    errors: list[str] = []
    for j, (va0, vb0) in enumerate(seeds):
        result = ppxf_utils.fit_free_two_component_losvd(
            templates_two_component=templates_two_component,
            component=component,
            galaxy=galaxy,
            noise=noise,
            velscale=velscale,
            lam=lam,
            lam_temp=lam_temp,
            goodpixels=goodpixels,
            start_velocity_a=va0,
            start_velocity_b=vb0,
            start_sigma_a=float(start_sigma),
            start_sigma_b=float(start_sigma),
            velocity_bounds_a=va_bounds,
            velocity_bounds_b=vb_bounds,
            sigma_bounds=sigma_bounds,
            degree=int(degree),
            mdegree=int(mdegree),
            regul=float(regul),
            noise_inv_cholesky=noise_inv_cholesky,
            keep_full=True,
        )
        if not result.success:
            errors.append(result.error_message)
            continue
        successes += 1
        if best is None or float(result.chi2_total) < float(best[1].chi2_total):
            best = (j, result)

    if best is None:
        n = np.asarray(galaxy).size
        return FreeTwoComponentFit(
            success=False,
            chi2_total=np.inf,
            reduced_chi2=np.inf,
            velocity=np.asarray([np.nan, np.nan]),
            sigma=np.asarray([np.nan, np.nan]),
            fraction_a=np.nan,
            bestfit=np.full(n, np.nan),
            weights=np.empty(0),
            start_index=-1,
            n_successful_starts=0,
            error_message="; ".join(errors[:4]) or "All multi-start calibration fits failed",
        )

    start_index, result = best
    weights = np.asarray(result.weights, dtype=float)
    n_basis = weights.size // 2
    wa = float(np.sum(weights[:n_basis]))
    wb = float(np.sum(weights[n_basis : 2 * n_basis]))
    fraction = wa / (wa + wb) if (wa + wb) > 0 else np.nan
    return FreeTwoComponentFit(
        success=True,
        chi2_total=float(result.chi2_total),
        reduced_chi2=float(result.reduced_chi2),
        velocity=np.asarray(result.velocity, dtype=float),
        sigma=np.asarray(result.sigma, dtype=float),
        fraction_a=float(fraction),
        bestfit=np.asarray(result.bestfit, dtype=float),
        weights=weights,
        start_index=int(start_index),
        n_successful_starts=int(successes),
    )


def _residual_scale_and_rho(
    residuals: np.ndarray,
    noise: np.ndarray,
    good: np.ndarray,
    *,
    max_lag: int,
    min_pairs: int,
    block_index: np.ndarray,
    n_blocks: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Measure per-bin s_i and blockwise rho_i(k) from calibration residuals."""
    residuals = np.asarray(residuals, dtype=float)
    noise = np.asarray(noise, dtype=float)
    good = np.asarray(good, dtype=bool)
    if not (residuals.shape == noise.shape == good.shape):
        raise ValueError("residuals/noise/good arrays must have identical shape")
    nbin = residuals.shape[0]
    scale = np.full(nbin, np.nan, dtype=float)
    rho = np.full((nbin, n_blocks, max_lag + 1), np.nan, dtype=float)
    pairs = np.zeros((nbin, n_blocks, max_lag + 1), dtype=int)

    for bid in range(nbin):
        gp = good[bid] & np.isfinite(residuals[bid]) & np.isfinite(noise[bid]) & (noise[bid] > 0)
        z = np.full(residuals.shape[1], np.nan, dtype=float)
        z[gp] = residuals[bid, gp] / noise[bid, gp]
        s = covariance.robust_standard_deviation(z[gp])
        if not np.isfinite(s) or s <= 0:
            continue
        scale[bid] = s
        z[gp] /= s
        rr, nn = covariance.lag_correlation(
            z,
            gp,
            max_lag=max_lag,
            min_pairs=min_pairs,
            block_index=block_index if n_blocks > 1 else None,
            n_blocks=n_blocks if n_blocks > 1 else None,
        )
        rho[bid] = rr
        pairs[bid] = nn
    return scale, rho, pairs


def estimate_candidate_model(
    name: str,
    *,
    residuals: np.ndarray,
    noise: np.ndarray,
    good: np.ndarray,
    max_lag: int,
    min_pairs: int,
    n_wavelength_blocks: int,
    n_bootstrap: int,
    bootstrap_confidence: float,
    random_seed: int,
) -> tuple[CandidateModel, dict[str, Any]]:
    """Estimate one M1--M4 covariance candidate from current fit residuals."""
    name = str(name).upper()
    if name not in MODEL_ORDER:
        raise ValueError(f"Unknown covariance candidate {name!r}")
    nbin, nlog = np.asarray(residuals).shape
    if name in {"M1", "M2"}:
        n_blocks = 1
        block_index = np.zeros(nlog, dtype=int)
    else:
        n_blocks = int(n_wavelength_blocks)
        block_index = covariance.equal_wavelength_blocks(nlog, n_blocks)

    scale, per_bin_rho, pair_count = _residual_scale_and_rho(
        residuals,
        noise,
        good,
        max_lag=int(max_lag),
        min_pairs=int(min_pairs),
        block_index=block_index,
        n_blocks=n_blocks,
    )
    if np.any(~np.isfinite(scale)):
        bad = np.flatnonzero(~np.isfinite(scale))
        raise RuntimeError(f"{name}: could not estimate a finite residual scale for bins {bad[:12].tolist()}")

    raw_band = covariance.bootstrap_simultaneous_band(
        per_bin_rho,
        n_bootstrap=int(n_bootstrap),
        confidence=float(bootstrap_confidence),
        random_seed=int(random_seed),
    )

    if name == "M1":
        adopted = np.zeros((1, max_lag + 1), dtype=float)
        adopted[:, 0] = 1.0
    elif name in {"M2", "M3"}:
        adopted = covariance.zero_insignificant_lags(raw_band.center, raw_band)
    else:  # M4: preserve bin-specific coefficients only where the global band requires them.
        support = ~raw_band.consistent_with_zero
        support[..., 0] = True
        adopted = np.asarray(per_bin_rho, dtype=float).copy()
        pooled = np.asarray(raw_band.center, dtype=float)
        for bid in range(nbin):
            missing = ~np.isfinite(adopted[bid])
            if np.any(missing):
                adopted[bid][missing] = np.broadcast_to(pooled, adopted[bid].shape)[missing]
            adopted[bid][~support] = 0.0
            adopted[bid, :, 0] = 1.0

    model = CandidateModel(
        name=name,
        scale=np.asarray(scale, dtype=float),
        rho=np.asarray(adopted, dtype=float),
        block_index=np.asarray(block_index, dtype=int),
        converged=False,
        n_iterations=0,
        final_scale_change=np.inf,
        final_rho_change=np.inf,
        max_lag=int(max_lag),
    )
    diagnostics = {
        "per_bin_raw_rho": per_bin_rho,
        "pair_count": pair_count,
        "raw_band": raw_band,
    }
    return model, diagnostics


def build_bin_whitener(
    model: CandidateModel,
    *,
    bin_id: int,
    noise: np.ndarray,
    good: np.ndarray,
    eigen_floor: float,
) -> covariance.WhiteningResult | None:
    """Build the frozen full-vector whitener for one bin/candidate.

    M1 remains on pPXF's ordinary diagonal-noise path and therefore returns
    ``None``; callers should multiply the one-dimensional noise vector by s_i.
    """
    bid = int(bin_id)
    if model.name == "M1":
        return None
    return covariance.build_inverse_cholesky(
        np.asarray(noise, dtype=float),
        np.asarray(good, dtype=bool),
        scale=float(model.scale[bid]),
        rho_by_block=model.rho_for_bin(bid),
        block_index=model.block_index if model.rho_for_bin(bid).shape[0] > 1 else None,
        eigen_floor=float(eigen_floor),
    )


def effective_noise_and_whitener(
    model: CandidateModel,
    *,
    bin_id: int,
    noise: np.ndarray,
    good: np.ndarray,
    eigen_floor: float,
) -> tuple[np.ndarray, np.ndarray | None, covariance.WhiteningResult | None]:
    """Return pPXF noise vector and optional cached whitener for one PowerBin."""
    bid = int(bin_id)
    noise = np.asarray(noise, dtype=float)
    if model.name == "M1":
        return noise * float(model.scale[bid]), None, None
    whitening = build_bin_whitener(
        model, bin_id=bid, noise=noise, good=good, eigen_floor=eigen_floor
    )
    return noise, whitening.inv_cholesky, whitening


def whitened_requirement_a(
    model: CandidateModel,
    *,
    residuals: np.ndarray,
    noise: np.ndarray,
    good: np.ndarray,
    max_lag: int,
    min_pairs: int,
    n_bootstrap: int,
    bootstrap_confidence: float,
    random_seed: int,
    eigen_floor: float,
    whitened_scale_tolerance: float,
    diagnostic_wavelength_blocks: int,
) -> CandidateModel:
    """Evaluate Requirement A on whitened residuals.

    Requirement A deliberately inspects the whitened residuals in wavelength
    blocks even when the candidate itself is stationary (M1/M2).  Otherwise a
    positive residual correlation around one spectral region could cancel a
    negative correlation elsewhere and make a poor stationary model appear to
    pass.  The simultaneous bootstrap band spans *all* diagnostic blocks and
    non-zero lags in one family-wise test.
    """
    residuals = np.asarray(residuals, dtype=float)
    noise = np.asarray(noise, dtype=float)
    good = np.asarray(good, dtype=bool)
    nbin, nlog = residuals.shape
    n_blocks = int(diagnostic_wavelength_blocks)
    block_index = covariance.equal_wavelength_blocks(nlog, n_blocks)
    per_bin_rho = np.full((nbin, n_blocks, max_lag + 1), np.nan, dtype=float)
    white_scale = np.full(nbin, np.nan, dtype=float)
    regularized = 0

    for bid in range(nbin):
        gp = np.flatnonzero(good[bid])
        if model.name == "M1":
            q = np.full(nlog, np.nan, dtype=float)
            q[gp] = residuals[bid, gp] / (noise[bid, gp] * float(model.scale[bid]))
        else:
            whitening = build_bin_whitener(
                model, bin_id=bid, noise=noise[bid], good=good[bid], eigen_floor=eigen_floor
            )
            regularized += int(whitening.regularized)
            # Excluded residuals are filled with zero only for matrix multiplication;
            # fitted rows are guaranteed not to couple to excluded dimensions.
            r = np.where(np.isfinite(residuals[bid]), residuals[bid], 0.0)
            q = covariance.whiten_residuals(r, whitening.inv_cholesky)
        white_scale[bid] = covariance.robust_standard_deviation(q[gp])
        rr, _ = covariance.lag_correlation(
            q,
            good[bid],
            max_lag=max_lag,
            min_pairs=min_pairs,
            block_index=block_index,
            n_blocks=n_blocks,
        )
        per_bin_rho[bid] = rr

    band = covariance.bootstrap_simultaneous_band(
        per_bin_rho,
        n_bootstrap=int(n_bootstrap),
        confidence=float(bootstrap_confidence),
        random_seed=int(random_seed) + 1000,
    )
    zero_ok = bool(np.all(band.consistent_with_zero[..., 1:]))
    median_scale = float(np.nanmedian(white_scale))
    scale_ok = bool(
        np.isfinite(median_scale)
        and abs(median_scale - 1.0) <= float(whitened_scale_tolerance)
    )

    model.requirement_a_pass = bool(zero_ok and scale_ok)
    model.whitened_scale_median = median_scale
    model.whitened_band_center = np.asarray(band.center, dtype=float)
    model.whitened_band_lower = np.asarray(band.lower, dtype=float)
    model.whitened_band_upper = np.asarray(band.upper, dtype=float)
    model.numerical_regularization_count = int(regularized)
    return model

def _bin_pixel_centroids(bin_map: np.ndarray, nbin: int) -> tuple[np.ndarray, np.ndarray]:
    y = np.full(nbin, np.nan, dtype=float)
    x = np.full(nbin, np.nan, dtype=float)
    for bid in range(nbin):
        yy, xx = np.where(np.asarray(bin_map) == bid)
        if yy.size:
            y[bid] = float(np.mean(yy))
            x[bid] = float(np.mean(xx))
    return y, x


def _pa_coordinates(x: np.ndarray, y: np.ndarray, pa_deg: float) -> tuple[np.ndarray, np.ndarray]:
    """Rotate east/north tangent-plane coordinates onto PA and cross-PA axes.

    Astronomical PA is measured from north through east, so a unit vector along
    PA is (sin PA, cos PA) in (east, north) coordinates.
    """
    th = np.deg2rad(float(pa_deg))
    xpar = np.asarray(x, dtype=float) * np.sin(th) + np.asarray(y, dtype=float) * np.cos(th)
    xperp = np.asarray(x, dtype=float) * np.cos(th) - np.asarray(y, dtype=float) * np.sin(th)
    return xpar, xperp


def select_representative_bins(
    source_table: Table,
    bin_map: np.ndarray,
    one_component_sigma: np.ndarray,
    *,
    pa_deg: float,
    n_radial_total: int = 12,
    corridor_diameter_factor: float = 1.0,
    add_two_sigma_bins: bool = True,
    two_sigma_inner_fraction: float = 0.10,
    two_sigma_outer_fraction: float = 0.95,
) -> Table:
    """Select 12 deterministic PA samples plus up to two 2-sigma candidates.

    With the default ``n_radial_total=12`` the selector places six target radii
    on each signed side of ``PA_kin`` at normalized radii 1/7,...,6/7.  The
    physical radial step therefore scales automatically with the usable size of
    the galaxy.  A PowerBin is normally eligible when its centroid lies within
    one median equivalent PowerBin diameter of the PA axis.  If no unused bin
    satisfies that corridor for one target, the nearest unused same-side bin is
    used and explicitly flagged so the validation set is never silently short.

    After the symmetric 12-bin radial set is frozen, the strongest off-center
    local one-component-dispersion maximum on each side is identified inside
    the same PA corridor.  If that candidate is already one of the radial bins,
    the existing row is marked as containing the 2-sigma candidate; otherwise
    it is *added*.  The default validation set therefore contains 12--14 unique
    PowerBins while preserving the original symmetric radial sampling.
    """
    if int(n_radial_total) < 4 or int(n_radial_total) % 2:
        raise ValueError("n_radial_total must be an even integer >= 4")
    n_side = int(n_radial_total) // 2
    required = ("BIN_ID", "AREA_ARCSEC2", "X_GEOM_ARCSEC", "Y_GEOM_ARCSEC", "RH3_SN")
    missing = [k for k in required if k not in source_table.colnames]
    if missing:
        raise ValueError("Representative-bin selection requires columns: " + ", ".join(missing))

    bid = np.asarray(source_table["BIN_ID"], dtype=int)
    x = np.asarray(source_table["X_GEOM_ARCSEC"], dtype=float)
    y = np.asarray(source_table["Y_GEOM_ARCSEC"], dtype=float)
    area = np.asarray(source_table["AREA_ARCSEC2"], dtype=float)
    sn = np.asarray(source_table["RH3_SN"], dtype=float)
    sigma = np.asarray(one_component_sigma, dtype=float)
    if sigma.shape != bid.shape:
        raise ValueError("one_component_sigma length must equal source_table length")
    xpar, xperp = _pa_coordinates(x, y, float(pa_deg))
    med_area = float(np.nanmedian(area[np.isfinite(area) & (area > 0)]))
    if not np.isfinite(med_area) or med_area <= 0:
        raise ValueError("Could not determine a positive median PowerBin area")
    median_diameter = 2.0 * np.sqrt(med_area / np.pi)
    corridor = float(corridor_diameter_factor) * median_diameter
    ypix, xpix = _bin_pixel_centroids(np.asarray(bin_map, dtype=int), len(source_table))

    selected: list[dict[str, Any]] = []
    used: set[int] = set()
    side_extent: dict[int, float] = {}
    radial_positions = np.arange(1, n_side + 1, dtype=float) / float(n_side + 1)

    def add_row(
        index: int,
        *,
        reason: str,
        side: int,
        target_norm: float | None,
        fallback: bool,
        base_radial: bool,
        added_for_sigma: bool,
        two_sigma_corridor_fallback: bool = False,
    ):
        b = int(bid[index])
        if b in used:
            return
        used.add(b)
        selected.append({
            "BIN_ID": b,
            "SELECTION_REASON": reason,
            "PA_SIDE": int(side),
            "BASE_RADIAL_SELECTION": bool(base_radial),
            "IS_2SIGMA_PEAK": bool(added_for_sigma),
            "ADDED_FOR_2SIGMA": bool(added_for_sigma),
            "TWO_SIGMA_CORRIDOR_FALLBACK": bool(two_sigma_corridor_fallback),
            "TARGET_NORMALIZED_RADIUS": np.nan if target_norm is None else float(target_norm),
            "SIGNED_PA_RADIUS_ARCSEC": float(xpar[index]),
            "ABS_PA_RADIUS_ARCSEC": float(abs(xpar[index])),
            "X_PERP_ARCSEC": float(xperp[index]),
            "CORRIDOR_HALF_WIDTH_ARCSEC": float(corridor),
            "MEDIAN_EQUIVALENT_BIN_DIAMETER_ARCSEC": float(median_diameter),
            "CORRIDOR_FALLBACK": bool(fallback),
            "ONECOMP_SIGMA_KMS": float(sigma[index]),
            "RH3_SN": float(sn[index]),
            "AREA_ARCSEC2": float(area[index]),
            "NPIX_RH3": int(source_table["NPIX_RH3"][index]) if "NPIX_RH3" in source_table.colnames else -1,
            "X_GEOM_ARCSEC": float(x[index]),
            "Y_GEOM_ARCSEC": float(y[index]),
            "X_CENTROID_PIXEL": float(xpix[b]) if 0 <= b < xpix.size else np.nan,
            "Y_CENTROID_PIXEL": float(ypix[b]) if 0 <= b < ypix.size else np.nan,
        })

    for side in (-1, +1):
        same = np.flatnonzero(np.isfinite(xpar) & (side * xpar > 0))
        if same.size < n_side:
            raise RuntimeError(f"Only {same.size} PowerBins lie on PA side {side:+d}; need {n_side}")
        on_axis = same[np.isfinite(xperp[same]) & (np.abs(xperp[same]) <= corridor)]
        extent_pool = on_axis if on_axis.size else same
        extent = float(np.max(np.abs(xpar[extent_pool])))
        if not np.isfinite(extent) or extent <= 0:
            raise RuntimeError(f"Could not determine a positive PA extent on side {side:+d}")
        side_extent[side] = extent
        for j, u in enumerate(radial_positions, start=1):
            target = side * float(u) * extent
            available = np.asarray([idx for idx in same if int(bid[idx]) not in used], dtype=int)
            in_corridor = available[np.isfinite(xperp[available]) & (np.abs(xperp[available]) <= corridor)]
            fallback = False
            pool = in_corridor
            if pool.size == 0:
                pool = available
                fallback = True
            if pool.size == 0:
                raise RuntimeError(f"No unused PowerBin remains on PA side {side:+d}")
            if fallback:
                # Locked fallback rule: when the median-diameter corridor contains
                # no unused bin, preserve the requested side and choose the centroid
                # closest to PA_kin; radial proximity breaks ties.
                order = np.lexsort((np.abs(xpar[pool] - target), np.abs(xperp[pool])))
            else:
                # Inside the corridor, select the PowerBin centroid closest to the
                # requested point on the PA axis; perpendicular distance breaks ties.
                dist2 = (xpar[pool] - target) ** 2 + xperp[pool] ** 2
                order = np.lexsort((np.abs(xperp[pool]), dist2))
            idx = int(pool[order[0]])
            add_row(
                idx,
                reason=f"radial_{'positive' if side > 0 else 'negative'}_{j}",
                side=side,
                target_norm=float(u),
                fallback=fallback,
                base_radial=True,
                added_for_sigma=False,
                two_sigma_corridor_fallback=False,
            )

    # Mark or add one off-centre dispersion candidate per side.  The 10--95%
    # default radial window prevents the central dispersion peak and extreme
    # aperture edge from being mislabeled as the classical off-centre 2-sigma
    # region.  This is a deterministic validation-bin locator, not a claim that
    # the one-component peak alone proves a physical two-sigma galaxy.
    two_sigma_candidates: dict[int, int | None] = {-1: None, +1: None}
    if add_two_sigma_bins:
        for side in (-1, +1):
            extent = side_extent[side]
            radial_pool = np.flatnonzero(
                np.isfinite(xpar)
                & np.isfinite(xperp)
                & np.isfinite(sigma)
                & (side * xpar > 0)
                & (np.abs(xpar) >= float(two_sigma_inner_fraction) * extent)
                & (np.abs(xpar) <= float(two_sigma_outer_fraction) * extent)
            )
            if radial_pool.size == 0:
                continue
            same = radial_pool[np.abs(xperp[radial_pool]) <= corridor]
            sigma_corridor_fallback = False
            if same.size == 0:
                # Guarantee a same-side off-centre candidate whenever finite sigma
                # measurements exist.  The fallback is explicit in the ECSV/QC.
                min_perp = float(np.nanmin(np.abs(xperp[radial_pool])))
                tol_perp = max(1.0e-9, 1.0e-6 * max(1.0, min_perp))
                same = radial_pool[np.abs(np.abs(xperp[radial_pool]) - min_perp) <= tol_perp]
                sigma_corridor_fallback = True
            order_r = same[np.argsort(np.abs(xpar[same]))]
            local: list[int] = []
            for j, idx in enumerate(order_r):
                left = sigma[order_r[j - 1]] if j > 0 else -np.inf
                right = sigma[order_r[j + 1]] if j + 1 < order_r.size else -np.inf
                if sigma[idx] >= left and sigma[idx] >= right:
                    local.append(int(idx))
            pool = np.asarray(local if local else same.tolist(), dtype=int)
            idx = int(pool[np.nanargmax(sigma[pool])])
            two_sigma_candidates[side] = int(bid[idx])
            if int(bid[idx]) in used:
                for row in selected:
                    if int(row["BIN_ID"]) == int(bid[idx]):
                        row["IS_2SIGMA_PEAK"] = True
                        row["TWO_SIGMA_CORRIDOR_FALLBACK"] = bool(sigma_corridor_fallback)
                        break
            else:
                add_row(
                    idx,
                    reason=f"two_sigma_peak_{'positive' if side > 0 else 'negative'}",
                    side=side,
                    target_norm=None,
                    fallback=False,
                    base_radial=False,
                    added_for_sigma=True,
                    two_sigma_corridor_fallback=bool(sigma_corridor_fallback),
                )

    out = Table(rows=selected)
    out.meta["PA_KIN_DEG"] = float(pa_deg)
    out.meta["N_RADIAL_REQUESTED"] = int(n_radial_total)
    out.meta["N_SELECTED"] = len(out)
    out.meta["N_ADDED_FOR_2SIGMA"] = int(np.sum(np.asarray(out["ADDED_FOR_2SIGMA"], dtype=bool)))
    out.meta["N_2SIGMA_CORRIDOR_FALLBACK"] = int(np.sum(np.asarray(out["TWO_SIGMA_CORRIDOR_FALLBACK"], dtype=bool)))
    out.meta["N_2SIGMA_CANDIDATES_REPRESENTED"] = int(np.sum([v is not None for v in two_sigma_candidates.values()]))
    out.meta["TWO_SIGMA_CANDIDATE_NEG_BIN_ID"] = -1 if two_sigma_candidates[-1] is None else int(two_sigma_candidates[-1])
    out.meta["TWO_SIGMA_CANDIDATE_POS_BIN_ID"] = -1 if two_sigma_candidates[+1] is None else int(two_sigma_candidates[+1])
    out.meta["MEDIAN_EQUIVALENT_BIN_DIAMETER_ARCSEC"] = float(median_diameter)
    out.meta["CORRIDOR_HALF_WIDTH_ARCSEC"] = float(corridor)
    out.meta["RADIAL_TARGETS_NORMALIZED"] = [float(x) for x in radial_positions]
    return out

def _profile_interval(profile_delta: np.ndarray, threshold: float) -> tuple[int, int] | None:
    use = np.flatnonzero(np.isfinite(profile_delta) & (profile_delta <= float(threshold)))
    if use.size == 0:
        return None
    return int(use[0]), int(use[-1])


def _cube_summary_indices(chi2: np.ndarray, thresholds: tuple[float, ...]) -> dict[str, Any]:
    arr = np.asarray(chi2, dtype=float)
    finite = np.isfinite(arr)
    if not np.any(finite):
        raise ValueError("Likelihood cube has no finite state")
    best_flat = int(np.nanargmin(np.where(finite, arr, np.nan)))
    best = tuple(int(x) for x in np.unravel_index(best_flat, arr.shape))
    delta = arr - float(arr[best])
    result: dict[str, Any] = {"best": best, "intervals": {}}
    for axis in range(3):
        reduce_axes = tuple(a for a in range(3) if a != axis)
        prof = np.nanmin(np.where(np.isfinite(delta), delta, np.nan), axis=reduce_axes)
        result["intervals"][axis] = {
            float(t): _profile_interval(prof, float(t)) for t in thresholds
        }
    return result


def compare_likelihood_cubes(
    cubes_a: dict[int, np.ndarray],
    cubes_b: dict[int, np.ndarray],
    *,
    max_cell_shift: int = 1,
    delta_chi2_thresholds: tuple[float, ...] = (1.0, 4.0),
) -> ModelComparison:
    """Compare Requirement-B scientific outputs at the production-grid resolution.

    Agreement is defined relative to what the production grid can actually
    resolve: best-state indices and the edges of one-dimensional profile-
    likelihood regions at the configured Delta-chi2 levels may move by at most
    ``max_cell_shift`` grid cells on any axis.
    """
    common = sorted(set(cubes_a) & set(cubes_b))
    if set(cubes_a) != set(cubes_b) or not common:
        raise ValueError("Candidate cube dictionaries must contain the same representative bins")
    worst_min = 0
    worst_edge = 0
    failures: list[int] = []
    details: list[dict[str, Any]] = []
    for bid in common:
        sa = _cube_summary_indices(cubes_a[bid], delta_chi2_thresholds)
        sb = _cube_summary_indices(cubes_b[bid], delta_chi2_thresholds)
        min_shift = int(max(abs(a - b) for a, b in zip(sa["best"], sb["best"])))
        edge_shift = 0
        missing_interval = False
        for axis in range(3):
            for threshold in delta_chi2_thresholds:
                ia = sa["intervals"][axis][float(threshold)]
                ib = sb["intervals"][axis][float(threshold)]
                if ia is None or ib is None:
                    missing_interval = True
                    edge_shift = max(edge_shift, max_cell_shift + 1)
                    continue
                edge_shift = max(edge_shift, abs(ia[0] - ib[0]), abs(ia[1] - ib[1]))
        agree = (not missing_interval) and min_shift <= int(max_cell_shift) and edge_shift <= int(max_cell_shift)
        if not agree:
            failures.append(int(bid))
        worst_min = max(worst_min, min_shift)
        worst_edge = max(worst_edge, edge_shift)
        details.append({
            "bin_id": int(bid),
            "best_a": sa["best"],
            "best_b": sb["best"],
            "minimum_cell_shift": min_shift,
            "interval_edge_cell_shift": edge_shift,
            "agree": bool(agree),
        })
    return ModelComparison(
        agree=not failures,
        max_minimum_cell_shift=int(worst_min),
        max_interval_edge_cell_shift=int(worst_edge),
        failing_bins=tuple(failures),
        details=tuple(details),
    )


def choose_simplest_stable_model(
    requirement_a: dict[str, bool],
    pairwise: dict[tuple[str, str], ModelComparison],
) -> str:
    """Apply the locked M1--M4 selection rule.

    The least-complex model must pass Requirement A and agree with every more
    complex model that also passes Requirement A.  M4 is *not* selected merely
    because it is the final model in the hierarchy: if simpler models pass A but
    all passing models materially disagree, selection fails and the production
    grid must not start.  M4 is accepted directly only when it is the sole model
    that passes residual adequacy.
    """
    passing = [m for m in MODEL_ORDER if bool(requirement_a.get(m, False))]
    if not passing:
        raise RuntimeError("COVARIANCE_MODEL_NO_REQUIREMENT_A_PASS | No covariance candidate whitened the residuals adequately")

    for i, model in enumerate(MODEL_ORDER[:-1]):
        if model not in passing:
            continue
        more = [m for m in MODEL_ORDER[i + 1 :] if m in passing]
        if not more:
            return model
        if all(pairwise[(model, other)].agree for other in more):
            return model

    if passing == ["M4"]:
        return "M4"

    raise RuntimeError(
        "COVARIANCE_MODEL_STABILITY_FAILURE | Residual-adequate covariance models "
        "do not yield scientifically stable representative likelihood surfaces. "
        "Do not launch the production grid; inspect template mismatch, wavelength "
        "non-stationarity, bin-dependent covariance, masks, and LSF treatment."
    )


def save_candidate_models(models: dict[str, CandidateModel], path: str | Path) -> Path:
    """Save compact calibrated covariance candidates without dense Nlambda^2 matrices."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload: dict[str, Any] = {"model_order": np.asarray(MODEL_ORDER)}
    for name, model in models.items():
        prefix = f"{name}_"
        payload[prefix + "scale"] = np.asarray(model.scale, dtype=float)
        payload[prefix + "rho"] = np.asarray(model.rho, dtype=float)
        payload[prefix + "block_index"] = np.asarray(model.block_index, dtype=np.int16)
        payload[prefix + "converged"] = np.asarray(bool(model.converged))
        payload[prefix + "n_iterations"] = np.asarray(int(model.n_iterations))
        payload[prefix + "final_scale_change"] = np.asarray(float(model.final_scale_change))
        payload[prefix + "final_rho_change"] = np.asarray(float(model.final_rho_change))
        payload[prefix + "requirement_a_pass"] = np.asarray(bool(model.requirement_a_pass))
        payload[prefix + "whitened_scale_median"] = np.asarray(float(model.whitened_scale_median))
        payload[prefix + "whitened_band_center"] = np.asarray(model.whitened_band_center, dtype=float)
        payload[prefix + "whitened_band_lower"] = np.asarray(model.whitened_band_lower, dtype=float)
        payload[prefix + "whitened_band_upper"] = np.asarray(model.whitened_band_upper, dtype=float)
        payload[prefix + "max_lag"] = np.asarray(int(model.max_lag))
        payload[prefix + "numerical_regularization_count"] = np.asarray(int(model.numerical_regularization_count))
    np.savez_compressed(path, **payload)
    return path


def load_candidate_models(path: str | Path) -> dict[str, CandidateModel]:
    """Reload compact covariance candidates saved by :func:`save_candidate_models`."""
    out: dict[str, CandidateModel] = {}
    with np.load(Path(path), allow_pickle=False) as data:
        order = [str(x) for x in np.asarray(data["model_order"]).tolist()]
        for name in order:
            p = f"{name}_"
            out[name] = CandidateModel(
                name=name,
                scale=np.asarray(data[p + "scale"], dtype=float),
                rho=np.asarray(data[p + "rho"], dtype=float),
                block_index=np.asarray(data[p + "block_index"], dtype=int),
                converged=bool(np.asarray(data[p + "converged"]).reshape(-1)[0]),
                n_iterations=int(np.asarray(data[p + "n_iterations"]).reshape(-1)[0]),
                final_scale_change=float(np.asarray(data[p + "final_scale_change"]).reshape(-1)[0]),
                final_rho_change=float(np.asarray(data[p + "final_rho_change"]).reshape(-1)[0]),
                requirement_a_pass=bool(np.asarray(data[p + "requirement_a_pass"]).reshape(-1)[0]),
                whitened_scale_median=float(np.asarray(data[p + "whitened_scale_median"]).reshape(-1)[0]),
                whitened_band_center=np.asarray(data[p + "whitened_band_center"], dtype=float),
                whitened_band_lower=np.asarray(data[p + "whitened_band_lower"], dtype=float),
                whitened_band_upper=np.asarray(data[p + "whitened_band_upper"], dtype=float),
                max_lag=int(np.asarray(data[p + "max_lag"]).reshape(-1)[0]),
                numerical_regularization_count=int(np.asarray(data[p + "numerical_regularization_count"]).reshape(-1)[0]),
            )
    return out

# -----------------------------------------------------------------------------
# High-level Script-3 orchestration
# -----------------------------------------------------------------------------

from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
import hashlib
import json
import time

from . import ppxf_grid

_CAL_WORKER: dict[str, Any] = {}


@dataclass
class CalibrationRunResult:
    selected_model_name: str
    selected_model: CandidateModel
    models: dict[str, CandidateModel]
    representative_bins: Table
    initial_one_velocity: np.ndarray
    initial_one_sigma: np.ndarray
    covariance_product: Path
    calibration_fit_product: Path
    validation_grid_product: Path
    selection_json: Path
    selected_model_hash: str


def _normalize_spectrum(galaxy: np.ndarray, noise: np.ndarray, good: np.ndarray):
    gp = np.flatnonzero(np.asarray(good, dtype=bool) & np.isfinite(galaxy) & np.isfinite(noise) & (np.asarray(noise) > 0))
    if gp.size == 0:
        return np.asarray(galaxy, dtype=float), np.asarray(noise, dtype=float), np.nan
    scale = float(np.nanmedian(np.abs(np.asarray(galaxy, dtype=float)[gp])))
    if not np.isfinite(scale) or scale <= 0:
        scale = float(np.nanpercentile(np.abs(np.asarray(galaxy, dtype=float)[gp]), 75.0))
    if not np.isfinite(scale) or scale <= 0:
        return np.asarray(galaxy, dtype=float), np.asarray(noise, dtype=float), np.nan
    return np.asarray(galaxy, dtype=float) / scale, np.asarray(noise, dtype=float) / scale, scale


def _normalize_all(gal_all, noise_all, good_all):
    nbin = np.asarray(gal_all).shape[0]
    gn = np.full_like(np.asarray(gal_all, dtype=float), np.nan, dtype=float)
    nn = np.full_like(np.asarray(noise_all, dtype=float), np.nan, dtype=float)
    scales = np.full(nbin, np.nan, dtype=float)
    for bid in range(nbin):
        gn[bid], nn[bid], scales[bid] = _normalize_spectrum(gal_all[bid], noise_all[bid], good_all[bid])
    if np.any(~np.isfinite(scales)):
        bad = np.flatnonzero(~np.isfinite(scales))
        raise RuntimeError(f"Covariance calibration could not normalize bins {bad[:12].tolist()}")
    return gn, nn, scales


def _cal_worker_init(constants: dict[str, Any]) -> None:
    global _CAL_WORKER
    _CAL_WORKER = constants


def _single_control_worker(payload):
    """Run the preliminary 1-C control without making it a hard science gate.

    The two-component residual calibration is the scientifically relevant fit. A
    failed one-component control therefore returns a flagged row and lets the
    multi-start calibration fall back to the configured systemic start velocity.
    This preserves the pre-existing Script-3 philosophy that a 1-C failure alone
    must not invalidate a spectrum that can still be fit successfully with two
    stellar components.
    """
    bid, galaxy, noise, good = payload
    c = _CAL_WORKER
    gp = np.flatnonzero(np.asarray(good, dtype=bool))
    result = ppxf_utils.fit_single_losvd(
        templates=c["templates_single"], galaxy=galaxy, noise=noise,
        velscale=c["velscale"], lam=c["lam"], lam_temp=c["lam_temp"],
        goodpixels=gp, start_velocity=c["single_start_velocity"],
        start_sigma=c["sigma_start"], velocity_bounds=c["single_velocity_bounds"],
        sigma_bounds=c["sigma_bounds"], degree=c["degree"], mdegree=c["mdegree"],
        regul=c["regul"], keep_full=True,
    )
    if not result.success:
        return {
            "bin_id": int(bid),
            "success": False,
            "velocity": np.nan,
            "sigma": np.nan,
            "chi2": np.inf,
            "bestfit": np.full(np.asarray(galaxy).shape, np.nan, dtype=np.float32),
            "error_message": str(result.error_message or "unknown one-component pPXF failure"),
        }
    return {
        "bin_id": int(bid),
        "success": True,
        "velocity": float(result.velocity[0]),
        "sigma": float(result.sigma[0]),
        "chi2": float(result.chi2_total),
        "bestfit": np.asarray(result.bestfit, dtype=np.float32),
        "error_message": "",
    }


def _free_calibration_worker(payload):
    bid, galaxy, noise, good, single_velocity, model_payload = payload
    c = _CAL_WORKER
    gp = np.flatnonzero(np.asarray(good, dtype=bool))
    if model_payload is None:
        effective_noise = noise
        W = None
        regularized = False
    else:
        model_name = model_payload["name"]
        scale = float(model_payload["scale"])
        if model_name == "M1":
            effective_noise = np.asarray(noise, dtype=float) * scale
            W = None
            regularized = False
        else:
            rho = np.asarray(model_payload["rho"], dtype=float)
            block_index = np.asarray(model_payload["block_index"], dtype=int)
            whitening = covariance.build_inverse_cholesky(
                noise, good, scale=scale, rho_by_block=rho,
                block_index=block_index if rho.shape[0] > 1 else None,
                eigen_floor=float(c["eigen_floor"]),
            )
            effective_noise = noise
            W = whitening.inv_cholesky
            regularized = whitening.regularized

    fit = fit_free_two_component_multistart(
        templates_two_component=c["templates_two"], component=c["component"],
        galaxy=galaxy, noise=effective_noise, velscale=c["velscale"],
        lam=c["lam"], lam_temp=c["lam_temp"], goodpixels=gp,
        single_velocity=float(single_velocity), start_sigma=c["sigma_start"],
        va_bounds=c["va_bounds"], vb_bounds=c["vb_bounds"], sigma_bounds=c["sigma_bounds"],
        degree=c["degree"], mdegree=c["mdegree"], regul=c["regul"],
        separation_fractions=c["separation_fractions"], noise_inv_cholesky=W,
    )
    if not fit.success:
        raise RuntimeError(f"Covariance free-two-component calibration failed for bin {bid}: {fit.error_message}")
    return {
        "bin_id": int(bid),
        "velocity": np.asarray(fit.velocity, dtype=np.float32),
        "sigma": np.asarray(fit.sigma, dtype=np.float32),
        "fraction_a": float(fit.fraction_a),
        "chi2": float(fit.chi2_total),
        "reduced_chi2": float(fit.reduced_chi2),
        "bestfit": np.asarray(fit.bestfit, dtype=np.float32),
        "start_index": int(fit.start_index),
        "n_successful_starts": int(fit.n_successful_starts),
        "regularized": bool(regularized),
    }


def _validation_grid_worker(payload):
    model_name, bid, galaxy, noise, good, model_payload = payload
    c = _CAL_WORKER
    gp = np.flatnonzero(np.asarray(good, dtype=bool))
    scale = float(model_payload["scale"])
    if model_name == "M1":
        effective_noise = np.asarray(noise, dtype=float) * scale
        W = None
    else:
        rho = np.asarray(model_payload["rho"], dtype=float)
        block_index = np.asarray(model_payload["block_index"], dtype=int)
        whitening = covariance.build_inverse_cholesky(
            noise, good, scale=scale, rho_by_block=rho,
            block_index=block_index if rho.shape[0] > 1 else None,
            eigen_floor=float(c["eigen_floor"]),
        )
        effective_noise = noise
        W = whitening.inv_cholesky
    cube = ppxf_grid.build_rh3_likelihood_cube(
        templates_two_component=c["templates_two"], component=c["component"],
        galaxy=galaxy, noise=effective_noise, velscale=c["velscale"],
        lam=c["lam"], lam_temp=c["lam_temp"], goodpixels=gp,
        va_grid=c["va_grid"], vb_grid=c["vb_grid"], fa_grid=c["fa_grid"],
        sigma_start_a=c["sigma_start"], sigma_start_b=c["sigma_start"],
        sigma_bounds=c["sigma_bounds"], degree=c["degree"], mdegree=c["mdegree"],
        regul=c["regul"], noise_inv_cholesky=W,
        sigma_boundary_tolerance_kms=c["sigma_boundary_tolerance"],
    )
    if cube.best_index is None:
        raise RuntimeError(f"Representative covariance grid failed completely for {model_name} bin {bid}")
    return {
        "model": str(model_name), "bin_id": int(bid),
        "chi2": np.asarray(cube.chi2_total, dtype=np.float64),
        "fit_status": np.asarray(cube.fit_status, dtype=np.int8),
    }


def _parallel_map(worker, payloads, constants, *, workers: int, logger, label: str):
    """Small crash-safe ProcessPool helper for calibration phases."""
    payloads = list(payloads)
    results = []
    start = time.perf_counter()
    with ProcessPoolExecutor(max_workers=int(workers), initializer=_cal_worker_init, initargs=(constants,)) as pool:
        futures = {pool.submit(worker, payload): payload for payload in payloads}
        outstanding = set(futures)
        done_count = 0
        while outstanding:
            done, outstanding = wait(outstanding, return_when=FIRST_COMPLETED)
            for future in done:
                try:
                    results.append(future.result())
                except Exception:
                    for other in outstanding:
                        other.cancel()
                    raise
                done_count += 1
                if done_count == 1 or done_count == len(payloads) or done_count % max(1, len(payloads)//10) == 0:
                    logger.info("%s: %d/%d complete (%.1f min)", label, done_count, len(payloads), (time.perf_counter()-start)/60.0)
    return results


def _fit_result_arrays(results, nbin: int, nlog: int):
    velocity = np.full((nbin, 2), np.nan, dtype=float)
    sigma = np.full((nbin, 2), np.nan, dtype=float)
    fraction = np.full(nbin, np.nan, dtype=float)
    chi2 = np.full(nbin, np.nan, dtype=float)
    rchi2 = np.full(nbin, np.nan, dtype=float)
    bestfit = np.full((nbin, nlog), np.nan, dtype=np.float32)
    start_index = np.full(nbin, -1, dtype=int)
    nstarts = np.zeros(nbin, dtype=int)
    regularized = np.zeros(nbin, dtype=bool)
    for row in results:
        bid = int(row["bin_id"])
        velocity[bid] = row["velocity"]
        sigma[bid] = row["sigma"]
        fraction[bid] = row["fraction_a"]
        chi2[bid] = row["chi2"]
        rchi2[bid] = row["reduced_chi2"]
        bestfit[bid] = row["bestfit"]
        start_index[bid] = row["start_index"]
        nstarts[bid] = row["n_successful_starts"]
        regularized[bid] = row["regularized"]
    return {
        "velocity": velocity, "sigma": sigma, "fraction": fraction,
        "chi2": chi2, "reduced_chi2": rchi2, "bestfit": bestfit,
        "start_index": start_index, "n_successful_starts": nstarts,
        "regularized": regularized,
    }


def _model_payload(model: CandidateModel, bid: int) -> dict[str, Any]:
    return {
        "name": model.name,
        "scale": float(model.scale[int(bid)]),
        "rho": np.asarray(model.rho_for_bin(int(bid)), dtype=float),
        "block_index": np.asarray(model.block_index, dtype=np.int16),
    }


def _candidate_hash(model: CandidateModel, wavelength: np.ndarray, good: np.ndarray) -> str:
    h = hashlib.sha256()
    h.update(str(model.name).encode())
    for arr in (model.scale, model.rho, model.block_index, np.asarray(wavelength, dtype=float), np.asarray(good, dtype=np.uint8)):
        a = np.ascontiguousarray(arr)
        h.update(str(a.dtype).encode()); h.update(str(a.shape).encode()); h.update(a.tobytes())
    return h.hexdigest()


def run_script03_covariance_calibration(
    *,
    cfg,
    logger,
    run,
    source_table: Table,
    bin_map: np.ndarray,
    templates_single: np.ndarray,
    templates_two: np.ndarray,
    component: np.ndarray,
    galaxy_all: np.ndarray,
    noise_all: np.ndarray,
    good_all: np.ndarray,
    wavelength: np.ndarray,
    lam_temp: np.ndarray,
    velscale: float,
    va_grid: np.ndarray,
    vb_grid: np.ndarray,
    fa_grid: np.ndarray,
    sigma_bounds: tuple[float, float],
    single_velocity_bounds: tuple[float, float],
    degree: int,
    mdegree: int,
    regul: float,
    sigma_start: float,
    sigma_boundary_tolerance: float,
    workers: int,
    plotting_module,
) -> CalibrationRunResult:
    """Run the complete pre-grid RH3 covariance calibration and M1--M4 selection."""
    validate_covariance_config(cfg)
    require_cached_ppxf_patch(str(getattr(cfg, "RH3_COVARIANCE_REQUIRED_PPXF_VERSION", "9.4.8")))

    nbin, nlog = np.asarray(galaxy_all).shape
    gal_n, noise_n, normalization = _normalize_all(galaxy_all, noise_all, good_all)
    max_lag = int(getattr(cfg, "RH3_COVARIANCE_MAX_LAG", 20))
    min_pairs = int(getattr(cfg, "RH3_COVARIANCE_MIN_PAIRS", 25))
    n_boot = int(getattr(cfg, "RH3_COVARIANCE_BOOTSTRAP_N", 2000))
    boot_ci = float(getattr(cfg, "RH3_COVARIANCE_BOOTSTRAP_CONFIDENCE", 0.95))
    random_seed = int(getattr(cfg, "RH3_COVARIANCE_RANDOM_SEED", 12345))
    n_blocks = int(getattr(cfg, "RH3_COVARIANCE_WAVELENGTH_BLOCKS", 3))
    conv_tol = float(getattr(cfg, "RH3_COVARIANCE_CONVERGENCE_TOL", 0.01))
    max_iter = int(getattr(cfg, "RH3_COVARIANCE_MAX_ITER", 5))
    eigen_floor = float(getattr(cfg, "RH3_COVARIANCE_EIGEN_FLOOR", 1.0e-8))
    white_tol = float(getattr(cfg, "RH3_COVARIANCE_WHITENED_SCALE_TOL", 0.05))
    sep_frac = tuple(float(x) for x in getattr(cfg, "RH3_COVARIANCE_CALIBRATION_SEPARATION_FRACTIONS", (0.15, 0.40, 0.75)))

    constants = {
        "templates_single": np.asarray(templates_single, dtype=float),
        "templates_two": np.asarray(templates_two, dtype=float),
        "component": np.asarray(component, dtype=int),
        "velscale": float(velscale), "lam": np.asarray(wavelength, dtype=float),
        "lam_temp": np.asarray(lam_temp, dtype=float),
        "single_start_velocity": float(getattr(cfg, "RH3_SINGLE_VELOCITY_START_KMS", 0.0)),
        "single_velocity_bounds": tuple(single_velocity_bounds),
        "sigma_start": float(sigma_start), "sigma_bounds": tuple(sigma_bounds),
        "degree": int(degree), "mdegree": int(mdegree), "regul": float(regul),
        "va_bounds": (float(np.min(va_grid)), float(np.max(va_grid))),
        "vb_bounds": (float(np.min(vb_grid)), float(np.max(vb_grid))),
        "separation_fractions": sep_frac, "eigen_floor": eigen_floor,
        "va_grid": np.asarray(va_grid, dtype=float), "vb_grid": np.asarray(vb_grid, dtype=float),
        "fa_grid": np.asarray(fa_grid, dtype=float),
        "sigma_boundary_tolerance": float(sigma_boundary_tolerance),
    }

    logger.info("Covariance calibration: initial one-component controls for %d PowerBins", nbin)
    single_rows = _parallel_map(
        _single_control_worker,
        [(bid, gal_n[bid], noise_n[bid], good_all[bid]) for bid in range(nbin)],
        constants, workers=workers, logger=logger, label="Covariance 1-C controls",
    )
    one_v = np.full(nbin, np.nan)
    one_s = np.full(nbin, np.nan)
    one_chi = np.full(nbin, np.nan)
    one_success = np.zeros(nbin, dtype=bool)
    one_error = np.full(nbin, "", dtype="U512")
    one_model = np.full((nbin, nlog), np.nan, dtype=np.float32)
    for row in single_rows:
        bid = int(row["bin_id"])
        one_success[bid] = bool(row.get("success", True))
        one_v[bid] = row["velocity"]
        one_s[bid] = row["sigma"]
        one_chi[bid] = row["chi2"]
        one_model[bid] = row["bestfit"]
        one_error[bid] = str(row.get("error_message", ""))[:512]
    failed_1c = np.flatnonzero(~one_success)
    if failed_1c.size:
        logger.warning(
            "COVARIANCE_ONE_COMPONENT_CONTROL_FAILURES | %d/%d preliminary 1-C controls failed. "
            "These bins remain eligible for the free two-component calibration, using the configured "
            "single-velocity start as the seed center. Example bins=%s",
            int(failed_1c.size), nbin, failed_1c[:12].tolist(),
        )
        one_v[failed_1c] = float(getattr(cfg, "RH3_SINGLE_VELOCITY_START_KMS", 0.0))

    # Deterministic representative selection is based on the preliminary 1-C
    # sigma map before covariance model choice, so all M1--M4 candidates are
    # tested on exactly the same spatial sample.
    pa = getattr(cfg, "PA_KIN_INITIAL_DEG", None)
    if pa is None or not np.isfinite(float(pa)):
        raise RuntimeError("RH3 covariance representative-bin selection requires finite PA_KIN_INITIAL_DEG")
    reps = select_representative_bins(
        source_table, bin_map, one_s, pa_deg=float(pa),
        n_radial_total=int(getattr(cfg, "RH3_COVARIANCE_VALIDATION_RADIAL_BINS", 12)),
        corridor_diameter_factor=float(getattr(cfg, "RH3_COVARIANCE_PA_CORRIDOR_DIAMETER_FACTOR", 1.0)),
        add_two_sigma_bins=bool(getattr(cfg, "RH3_COVARIANCE_ADD_2SIGMA_BINS", True)),
        two_sigma_inner_fraction=float(getattr(cfg, "RH3_COVARIANCE_2SIGMA_INNER_RADIUS_FRACTION", 0.10)),
        two_sigma_outer_fraction=float(getattr(cfg, "RH3_COVARIANCE_2SIGMA_OUTER_RADIUS_FRACTION", 0.95)),
    )
    rep_path = run.products_dir / "covariance_validation_bins.ecsv"
    reps.write(rep_path, format="ascii.ecsv", overwrite=True)
    n_fallback = int(np.sum(np.asarray(reps["CORRIDOR_FALLBACK"], dtype=bool)))
    if n_fallback:
        logger.warning(
            "COVARIANCE_REPRESENTATIVE_BIN_CORRIDOR_FALLBACK | %d/%d baseline radial selections "
            "had no unused centroid inside the median-diameter PA corridor; the nearest-to-PA "
            "same-side PowerBin was used and flagged in covariance_validation_bins.ecsv.",
            n_fallback, len(reps),
        )
    n_sigma_fallback = int(np.sum(np.asarray(reps["TWO_SIGMA_CORRIDOR_FALLBACK"], dtype=bool)))
    if n_sigma_fallback:
        logger.warning(
            "COVARIANCE_2SIGMA_CORRIDOR_FALLBACK | %d off-centre sigma candidate(s) had no finite "
            "PowerBin inside the nominal PA corridor; the closest-to-PA radial candidate was used "
            "and explicitly flagged.",
            n_sigma_fallback,
        )
    logger.info(
        "Covariance representative set: %d bins (%d fixed radial + %d added 2-sigma; %d 2-sigma candidate(s) represented)",
        len(reps),
        int(np.sum(np.asarray(reps["BASE_RADIAL_SELECTION"], dtype=bool))),
        int(np.sum(np.asarray(reps["ADDED_FOR_2SIGMA"], dtype=bool))),
        int(np.sum(np.asarray(reps["IS_2SIGMA_PEAK"], dtype=bool))),
    )

    plotting_module.plot_covariance_validation_bins_map(
        bin_map, one_s, reps, run.figures_dir / "covariance_validation_bins_on_sigma.png",
        title="Script-3 covariance validation bins on preliminary RH3 stellar dispersion",
        colorbar_label="One-component sigma (km/s)",
    )
    plotting_module.plot_covariance_validation_bins_map(
        bin_map, np.asarray(source_table["RH3_SN"], dtype=float), reps,
        run.figures_dir / "covariance_validation_bins_on_RH3_SN.png",
        title="Script-3 covariance validation bins on RH3 achieved S/N",
        colorbar_label="RH3 achieved S/N",
    )

    logger.info("Covariance calibration: baseline free two-component multi-start fits with formal diagonal noise")
    base_rows = _parallel_map(
        _free_calibration_worker,
        [(bid, gal_n[bid], noise_n[bid], good_all[bid], one_v[bid], None) for bid in range(nbin)],
        constants, workers=workers, logger=logger, label="Baseline covariance residual fits",
    )
    base = _fit_result_arrays(base_rows, nbin, nlog)
    base_resid = gal_n - np.asarray(base["bestfit"], dtype=float)
    plotting_module.plot_covariance_residual_stack(
        wavelength, base_resid, noise_n, good_all,
        run.figures_dir / "RH3_covariance_baseline_residual_stack.png",
        title="Baseline free two-component residuals before covariance calibration",
    )

    models: dict[str, CandidateModel] = {}
    final_fits: dict[str, dict[str, np.ndarray]] = {}
    final_raw_bands: dict[str, covariance.BootstrapBand] = {}
    iteration_history: list[dict[str, Any]] = []
    history_path = run.products_dir / "RH3_covariance_iteration_history.ecsv"

    for model_index, name in enumerate(MODEL_ORDER):
        logger.info("Covariance candidate %s: %s", name, MODEL_DESCRIPTION[name])
        current, diag = estimate_candidate_model(
            name, residuals=base_resid, noise=noise_n, good=good_all,
            max_lag=max_lag, min_pairs=min_pairs, n_wavelength_blocks=n_blocks,
            n_bootstrap=n_boot, bootstrap_confidence=boot_ci,
            random_seed=random_seed + 100 * model_index,
        )
        final_raw_bands[name] = diag["raw_band"]
        converged = False
        scale_change = np.inf; rho_change = np.inf
        fit_arrays = None
        for iteration in range(1, max_iter + 1):
            payloads = [
                (bid, gal_n[bid], noise_n[bid], good_all[bid], one_v[bid], _model_payload(current, bid))
                for bid in range(nbin)
            ]
            rows = _parallel_map(
                _free_calibration_worker, payloads, constants, workers=workers, logger=logger,
                label=f"{name} covariance iteration {iteration}",
            )
            fit_arrays = _fit_result_arrays(rows, nbin, nlog)
            residuals = gal_n - np.asarray(fit_arrays["bestfit"], dtype=float)
            updated, diag = estimate_candidate_model(
                name, residuals=residuals, noise=noise_n, good=good_all,
                max_lag=max_lag, min_pairs=min_pairs, n_wavelength_blocks=n_blocks,
                n_bootstrap=n_boot, bootstrap_confidence=boot_ci,
                random_seed=random_seed + 100 * model_index + iteration,
            )
            scale_change = covariance.relative_scale_change(current.scale, updated.scale)
            rho_change = covariance.correlation_change(current.rho, updated.rho)
            logger.info("%s iteration %d: max fractional Delta s=%.5f; max |Delta rho|=%.5f; convergence threshold=%.5f", name, iteration, scale_change, rho_change, conv_tol)
            iteration_history.append({
                "MODEL": name,
                "ITERATION": int(iteration),
                "MAX_FRACTIONAL_SCALE_CHANGE": float(scale_change),
                "MAX_ABS_RHO_CHANGE": float(rho_change),
                "CONVERGED": bool(scale_change < conv_tol and rho_change < conv_tol),
            })
            Table(rows=iteration_history).write(history_path, format="ascii.ecsv", overwrite=True)
            current = updated
            final_raw_bands[name] = diag["raw_band"]
            if scale_change < conv_tol and rho_change < conv_tol:
                converged = True
                current.converged = True; current.n_iterations = iteration
                current.final_scale_change = scale_change; current.final_rho_change = rho_change
                break

        if not converged:
            current.converged = False
            current.n_iterations = max_iter
            current.final_scale_change = scale_change
            current.final_rho_change = rho_change
            Table(rows=iteration_history).write(history_path, format="ascii.ecsv", overwrite=True)
            logger.warning(
                "COVARIANCE_CALIBRATION_MAXITER_REACHED | %s did not converge after %d iterations: "
                "max fractional Delta s=%.5f, max |Delta rho|=%.5f. Inspect calibration-fit residuals, "
                "template/LSF mismatch, atmospheric residuals, wavelength non-stationarity, low-information "
                "bins, and multi-start consistency. Production likelihood fitting is blocked.",
                name, max_iter, scale_change, rho_change,
            )
            failure_path = run.metadata_dir / "RH3_covariance_model_selection.json"
            failure_path.write_text(json.dumps({
                "selected_model": None,
                "model_order": list(MODEL_ORDER),
                "failed_model": name,
                "production_grid_allowed": False,
                "failure": "COVARIANCE_CALIBRATION_MAXITER_REACHED",
                "final_scale_change": float(scale_change),
                "final_rho_change": float(rho_change),
                "convergence_tolerance": float(conv_tol),
                "max_iterations": int(max_iter),
                "iteration_history_ecsv": str(history_path),
                "guidance": [
                    "Inspect free two-component calibration-fit residuals for coherent stellar-feature mismatch.",
                    "Check whether residual correlation changes strongly across wavelength blocks.",
                    "Check CRD_DRP atmospheric masks for surviving sky/telluric residuals.",
                    "Verify the RH3 LSF and template-resolution match.",
                    "Inspect dependence of s_i and rho_i on PowerBin size and S/N.",
                    "Check whether multi-start fits land in substantially different fit-quality basins.",
                    "Do not loosen the 0.01 convergence tolerance merely to force a production run.",
                ],
            }, indent=2, sort_keys=True) + "\n")
            raise RuntimeError(
                f"COVARIANCE_CALIBRATION_MAXITER_REACHED | {name} failed the locked {conv_tol:.3f} "
                f"convergence tolerance after {max_iter} iterations. See {history_path} and {failure_path}."
            )

        # One final fit under the converged model provides the residuals used by
        # Requirement A, rather than evaluating whitening on residuals produced
        # under the immediately preceding (merely close) covariance estimate.
        rows = _parallel_map(
            _free_calibration_worker,
            [(bid, gal_n[bid], noise_n[bid], good_all[bid], one_v[bid], _model_payload(current, bid)) for bid in range(nbin)],
            constants, workers=workers, logger=logger, label=f"{name} final converged residual fits",
        )
        fit_arrays = _fit_result_arrays(rows, nbin, nlog)
        residuals = gal_n - np.asarray(fit_arrays["bestfit"], dtype=float)
        current = whitened_requirement_a(
            current, residuals=residuals, noise=noise_n, good=good_all,
            max_lag=max_lag, min_pairs=min_pairs, n_bootstrap=n_boot,
            bootstrap_confidence=boot_ci, random_seed=random_seed + 5000 + model_index,
            eigen_floor=eigen_floor, whitened_scale_tolerance=white_tol,
            diagnostic_wavelength_blocks=n_blocks,
        )
        models[name] = current
        final_fits[name] = fit_arrays
        logger.info("%s Requirement A: %s | whitened median robust sigma=%.4f | simultaneous %.1f%% band contains zero at all non-zero lags=%s", name, "PASS" if current.requirement_a_pass else "FAIL", current.whitened_scale_median, 100*boot_ci, bool(np.all((current.whitened_band_lower[..., 1:] <= 0) & (current.whitened_band_upper[..., 1:] >= 0))))

        plotting_module.plot_covariance_lag_band(
            np.arange(max_lag + 1), final_raw_bands[name].center,
            final_raw_bands[name].lower, final_raw_bands[name].upper,
            run.figures_dir / f"RH3_covariance_{name}_raw_lag_band.png",
            title=f"{name} calibration residual lag correlation before whitening",
        )
        plotting_module.plot_covariance_lag_band(
            np.arange(max_lag + 1), current.whitened_band_center,
            current.whitened_band_lower, current.whitened_band_upper,
            run.figures_dir / f"RH3_covariance_{name}_whitened_lag_band.png",
            title=f"{name} whitened residual lag correlation",
        )
        plotting_module.plot_covariance_scale_vs_bin_properties(
            current.scale, source_table, run.figures_dir / f"RH3_covariance_{name}_scale_vs_bin_properties.png",
            title=f"{name} per-PowerBin empirical noise scale",
        )

    # Persist the complete convergence audit before Requirement B begins.
    Table(rows=iteration_history).write(history_path, format="ascii.ecsv", overwrite=True)
    cov_path = save_candidate_models(models, run.products_dir / "RH3_covariance_candidates.npz")
    fit_path = run.products_dir / "RH3_covariance_calibration_fits.npz"
    fit_payload: dict[str, Any] = {
        "wavelength": np.asarray(wavelength, dtype=float),
        "normalization_scale": normalization,
        "initial_one_success": one_success.astype(np.uint8),
        "initial_one_velocity": one_v,
        "initial_one_sigma": one_s,
        "initial_one_chi2": one_chi,
        "initial_one_error": one_error,
        "initial_one_model": one_model,
    }
    for name, fit in final_fits.items():
        for key, value in fit.items():
            fit_payload[f"{name}_{key}"] = value
    np.savez_compressed(fit_path, **fit_payload)

    # Requirement B: full production grid for every selected representative bin
    # under every covariance candidate.  No narrowed Delta-V stencil is used.
    rep_ids = [int(x) for x in np.asarray(reps["BIN_ID"], dtype=int)]
    tasks = []
    for name in MODEL_ORDER:
        model = models[name]
        for bid in rep_ids:
            tasks.append((name, bid, gal_n[bid], noise_n[bid], good_all[bid], _model_payload(model, bid)))
    logger.info("Requirement B: evaluating %d representative model/bin cubes = %d exact pPXF states", len(tasks), len(tasks)*int(np.size(va_grid)*np.size(vb_grid)*np.size(fa_grid)))
    grid_rows = _parallel_map(
        _validation_grid_worker, tasks, constants, workers=workers, logger=logger,
        label="Covariance Requirement-B full grids",
    )
    candidate_cubes: dict[str, dict[int, np.ndarray]] = {m: {} for m in MODEL_ORDER}
    grid_payload: dict[str, Any] = {
        "model_order": np.asarray(MODEL_ORDER), "representative_bin_ids": np.asarray(rep_ids, dtype=int),
        "VA_grid": np.asarray(va_grid, dtype=float), "VB_grid": np.asarray(vb_grid, dtype=float), "fA_grid": np.asarray(fa_grid, dtype=float),
    }
    for row in grid_rows:
        name = row["model"]; bid = int(row["bin_id"])
        candidate_cubes[name][bid] = row["chi2"]
        grid_payload[f"{name}_bin_{bid:04d}_chi2"] = row["chi2"]
        grid_payload[f"{name}_bin_{bid:04d}_status"] = row["fit_status"]
    grid_path = run.products_dir / "RH3_covariance_model_validation_grids.npz"
    np.savez_compressed(grid_path, **grid_payload)

    comparison_dir = run.figures_dir / "RH3_covariance_model_comparison_bins"
    comparison_dir.mkdir(parents=True, exist_ok=True)
    for bid in rep_ids:
        plotting_module.plot_covariance_model_grid_comparison(
            va_grid, vb_grid, fa_grid,
            {name: candidate_cubes[name][bid] for name in MODEL_ORDER},
            comparison_dir / f"bin_{bid:04d}.png",
            bin_id=bid,
        )

    pairwise: dict[tuple[str, str], ModelComparison] = {}
    pairwise_json: dict[str, Any] = {}
    comparison_rows: list[dict[str, Any]] = []
    max_shift = int(getattr(cfg, "RH3_COVARIANCE_MODEL_AGREEMENT_MAX_CELL_SHIFT", 1))
    thresholds = tuple(float(x) for x in getattr(cfg, "RH3_COVARIANCE_MODEL_AGREEMENT_DELTA_CHI2", (1.0, 4.0)))
    for i, left in enumerate(MODEL_ORDER):
        for right in MODEL_ORDER[i+1:]:
            comp = compare_likelihood_cubes(candidate_cubes[left], candidate_cubes[right], max_cell_shift=max_shift, delta_chi2_thresholds=thresholds)
            pairwise[(left, right)] = comp
            pairwise_json[f"{left}_vs_{right}"] = {
                "agree": comp.agree,
                "max_minimum_cell_shift": comp.max_minimum_cell_shift,
                "max_interval_edge_cell_shift": comp.max_interval_edge_cell_shift,
                "failing_bins": list(comp.failing_bins),
                "details": list(comp.details),
            }
            for detail in comp.details:
                comparison_rows.append({
                    "MODEL_SIMPLE": left,
                    "MODEL_COMPLEX": right,
                    "BIN_ID": int(detail["bin_id"]),
                    "BEST_SIMPLE": str(tuple(detail["best_a"])),
                    "BEST_COMPLEX": str(tuple(detail["best_b"])),
                    "MINIMUM_CELL_SHIFT": int(detail["minimum_cell_shift"]),
                    "PROFILE_INTERVAL_EDGE_CELL_SHIFT": int(detail["interval_edge_cell_shift"]),
                    "AGREE": bool(detail["agree"]),
                })
            logger.info("Requirement B %s vs %s: %s | max minimum shift=%d cell(s), max profile-interval edge shift=%d cell(s)", left, right, "AGREE" if comp.agree else "DIFFER", comp.max_minimum_cell_shift, comp.max_interval_edge_cell_shift)

    comparison_table_path = run.products_dir / "RH3_covariance_model_comparison.ecsv"
    Table(rows=comparison_rows).write(comparison_table_path, format="ascii.ecsv", overwrite=True)

    requirement_a = {name: bool(models[name].requirement_a_pass) for name in MODEL_ORDER}
    try:
        selected_name = choose_simplest_stable_model(requirement_a, pairwise)
    except RuntimeError as exc:
        failure_path = run.metadata_dir / "RH3_covariance_model_selection.json"
        failure_payload = {
            "selected_model": None,
            "model_order": list(MODEL_ORDER),
            "requirement_a": requirement_a,
            "requirement_b_pairwise": pairwise_json,
            "representative_bins_ecsv": str(rep_path),
            "model_comparison_ecsv": str(comparison_table_path),
            "n_representative_bins": len(rep_ids),
            "representative_bin_ids": rep_ids,
            "agreement_max_cell_shift": max_shift,
            "agreement_delta_chi2_thresholds": list(thresholds),
            "bootstrap_confidence": boot_ci,
            "bootstrap_n": n_boot,
            "convergence_tolerance": conv_tol,
            "max_iterations": max_iter,
            "production_grid_allowed": False,
            "failure": str(exc),
            "guidance": [
                "Inspect the multi-start calibration residual spectra for coherent stellar-feature mismatch.",
                "Compare raw and whitened lag bands for wavelength non-stationarity.",
                "Check CRD_DRP atmospheric masks for residual sky/telluric structure.",
                "Check the RH3 LSF model and template-resolution match.",
                "Inspect covariance scale versus PowerBin size/SN for bin-dependent behavior.",
                "Do not relax convergence/agreement thresholds solely to force a production run.",
            ],
        }
        failure_path.write_text(json.dumps(failure_payload, indent=2, sort_keys=True, default=str) + "\n")
        raise
    selected = models[selected_name]
    selected_hash = _candidate_hash(selected, wavelength, good_all)
    logger.info("Covariance model selection: %s | %s", selected_name, MODEL_DESCRIPTION[selected_name])
    selected_residuals = gal_n - np.asarray(final_fits[selected_name]["bestfit"], dtype=float)
    plotting_module.plot_covariance_residual_stack(
        wavelength, selected_residuals, noise_n, good_all,
        run.figures_dir / "RH3_covariance_selected_residual_stack.png",
        title=f"Selected {selected_name} covariance-calibration residuals",
    )

    selection_path = run.metadata_dir / "RH3_covariance_model_selection.json"
    selection_payload = {
        "selected_model": selected_name,
        "selected_model_description": MODEL_DESCRIPTION[selected_name],
        "selected_model_hash": selected_hash,
        "model_order": list(MODEL_ORDER),
        "requirement_a": requirement_a,
        "requirement_b_pairwise": pairwise_json,
        "representative_bins_ecsv": str(rep_path),
        "model_comparison_ecsv": str(comparison_table_path),
        "iteration_history_ecsv": str(history_path),
        "n_representative_bins": len(rep_ids),
        "representative_bin_ids": rep_ids,
        "agreement_max_cell_shift": max_shift,
        "agreement_delta_chi2_thresholds": list(thresholds),
        "bootstrap_confidence": boot_ci,
        "bootstrap_n": n_boot,
        "convergence_tolerance": conv_tol,
        "max_iterations": max_iter,
        "production_grid_allowed": True,
    }
    selection_path.write_text(json.dumps(selection_payload, indent=2, sort_keys=True, default=str) + "\n")

    return CalibrationRunResult(
        selected_model_name=selected_name,
        selected_model=selected,
        models=models,
        representative_bins=reps,
        initial_one_velocity=one_v,
        initial_one_sigma=one_s,
        covariance_product=cov_path,
        calibration_fit_product=fit_path,
        validation_grid_product=grid_path,
        selection_json=selection_path,
        selected_model_hash=selected_hash,
    )


def calibration_products_complete(run) -> bool:
    """Return True when the saved pre-production covariance decision is reusable."""
    required = [
        run.products_dir / "RH3_covariance_candidates.npz",
        run.products_dir / "RH3_covariance_calibration_fits.npz",
        run.products_dir / "RH3_covariance_model_validation_grids.npz",
        run.products_dir / "covariance_validation_bins.ecsv",
        run.metadata_dir / "RH3_covariance_model_selection.json",
    ]
    return all(Path(p).is_file() for p in required)


def load_calibration_run(run) -> CalibrationRunResult:
    """Reload a previously completed covariance-calibration decision for resume."""
    if not calibration_products_complete(run):
        raise FileNotFoundError("Saved Script-3 covariance calibration is incomplete")
    cov_path = run.products_dir / "RH3_covariance_candidates.npz"
    fit_path = run.products_dir / "RH3_covariance_calibration_fits.npz"
    grid_path = run.products_dir / "RH3_covariance_model_validation_grids.npz"
    rep_path = run.products_dir / "covariance_validation_bins.ecsv"
    selection_path = run.metadata_dir / "RH3_covariance_model_selection.json"
    models = load_candidate_models(cov_path)
    selection = json.loads(selection_path.read_text())
    selected_name = str(selection["selected_model"])
    if selected_name not in models:
        raise RuntimeError("Saved covariance selection names a model absent from RH3_covariance_candidates.npz")
    with np.load(fit_path, allow_pickle=False) as data:
        one_v = np.asarray(data["initial_one_velocity"], dtype=float)
        one_s = np.asarray(data["initial_one_sigma"], dtype=float)
    reps = Table.read(rep_path, format="ascii.ecsv")
    return CalibrationRunResult(
        selected_model_name=selected_name,
        selected_model=models[selected_name],
        models=models,
        representative_bins=reps,
        initial_one_velocity=one_v,
        initial_one_sigma=one_s,
        covariance_product=cov_path,
        calibration_fit_product=fit_path,
        validation_grid_product=grid_path,
        selection_json=selection_path,
        selected_model_hash=str(selection["selected_model_hash"]),
    )
