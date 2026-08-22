"""Standardized pPXF wrappers for CRD_DAP spectral likelihood work.

Script 3 uses these functions to ensure that every likelihood-grid state is
compared on exactly the same log-wavelength samples, with the same uncertainty
model and the same continuum treatment.  The wrappers intentionally import
pPXF lazily so lightweight CRD_DAP utilities and unit tests do not require the
optional science dependency merely to import the package.

The statistically important quantity returned here is the *total* chi-square
on the fixed good-pixel set.  ``pp.chi2`` is also retained as a reduced-chi2 QC
quantity, but relative likelihoods must be formed from differences in total
chi-square, not from differences in ``pp.chi2``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np


@dataclass
class PPXFStateResult:
    """Compact scalar/result bundle from one pPXF fit."""

    success: bool
    chi2_total: float
    reduced_chi2: float
    velocity: np.ndarray
    sigma: np.ndarray
    bestfit: np.ndarray | None = None
    weights: np.ndarray | None = None
    polyweights: np.ndarray | None = None
    error_message: str = ""
    pp: Any | None = None


def _import_ppxf():
    try:
        from ppxf.ppxf import ppxf
    except Exception as exc:  # pragma: no cover - depends on optional package
        raise ImportError(
            "Script 3 requires the optional pPXF dependency. Install the "
            "repository science extras before running this stage."
        ) from exc
    return ppxf


def _validate_spectrum_inputs(
    galaxy: np.ndarray,
    noise: np.ndarray,
    lam: np.ndarray,
    goodpixels: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    galaxy = np.asarray(galaxy, dtype=float)
    noise = np.asarray(noise, dtype=float)
    lam = np.asarray(lam, dtype=float)
    goodpixels = np.asarray(goodpixels, dtype=int)

    if galaxy.ndim != 1 or noise.ndim != 1 or lam.ndim != 1:
        raise ValueError("galaxy, noise, and lam must be one-dimensional.")
    if not (galaxy.size == noise.size == lam.size):
        raise ValueError("galaxy, noise, and lam must have identical lengths.")
    if goodpixels.ndim != 1 or goodpixels.size < 5:
        raise ValueError("goodpixels must contain at least five spectral samples.")
    if np.any(goodpixels < 0) or np.any(goodpixels >= galaxy.size):
        raise ValueError("goodpixels contains an out-of-range index.")
    if np.any(~np.isfinite(noise[goodpixels])) or np.any(noise[goodpixels] <= 0):
        raise ValueError("noise must be finite and positive on every good pixel.")
    if np.any(~np.isfinite(galaxy[goodpixels])):
        raise ValueError("galaxy must be finite on every good pixel.")
    if np.any(~np.isfinite(lam)) or np.any(np.diff(lam) <= 0):
        raise ValueError("lam must be finite and strictly increasing.")
    return galaxy, noise, lam, goodpixels


def total_chi2_diagonal(
    galaxy: np.ndarray,
    bestfit: np.ndarray,
    noise: np.ndarray,
    goodpixels: np.ndarray,
) -> float:
    """Return total chi-square for a diagonal spectral-noise model."""
    gp = np.asarray(goodpixels, dtype=int)
    resid = (np.asarray(galaxy, dtype=float)[gp] - np.asarray(bestfit, dtype=float)[gp])
    sig = np.asarray(noise, dtype=float)[gp]
    return float(np.sum((resid / sig) ** 2))


def _extract_solution(pp, n_components: int) -> tuple[np.ndarray, np.ndarray]:
    sol = pp.sol
    if n_components == 1:
        arr = np.asarray(sol, dtype=float)
        return np.asarray([arr[0]], dtype=float), np.asarray([arr[1]], dtype=float)

    if len(sol) != n_components:
        raise ValueError(
            f"Expected {n_components} pPXF component solutions but received {len(sol)}."
        )
    vel = np.asarray([np.asarray(x, dtype=float)[0] for x in sol], dtype=float)
    sig = np.asarray([np.asarray(x, dtype=float)[1] for x in sol], dtype=float)
    return vel, sig


def fit_single_losvd(
    *,
    templates: np.ndarray,
    galaxy: np.ndarray,
    noise: np.ndarray,
    velscale: float,
    lam: np.ndarray,
    lam_temp: np.ndarray,
    goodpixels: np.ndarray,
    start_velocity: float,
    start_sigma: float,
    velocity_bounds: tuple[float, float],
    sigma_bounds: tuple[float, float],
    degree: int = 4,
    mdegree: int = 0,
    regul: float = 0.0,
    keep_full: bool = True,
) -> PPXFStateResult:
    """Run the one-component RH3 control fit on a fixed good-pixel set."""
    ppxf = _import_ppxf()
    galaxy, noise, lam, goodpixels = _validate_spectrum_inputs(
        galaxy, noise, lam, goodpixels
    )
    templates = np.asarray(templates, dtype=float)
    lam_temp = np.asarray(lam_temp, dtype=float)
    if templates.ndim != 2 or templates.shape[0] != lam_temp.size:
        raise ValueError("templates must have shape (n_template_pixels, n_templates).")

    try:
        pp = ppxf(
            templates,
            galaxy,
            noise,
            float(velscale),
            start=[float(start_velocity), float(start_sigma)],
            moments=2,
            bounds=[
                [float(velocity_bounds[0]), float(velocity_bounds[1])],
                [float(sigma_bounds[0]), float(sigma_bounds[1])],
            ],
            goodpixels=goodpixels,
            degree=int(degree),
            mdegree=int(mdegree),
            regul=float(regul),
            lam=lam,
            lam_temp=lam_temp,
            clean=False,
            quiet=True,
        )
        vel, sig = _extract_solution(pp, 1)
        total = total_chi2_diagonal(galaxy, pp.bestfit, noise, goodpixels)
        return PPXFStateResult(
            success=True,
            chi2_total=total,
            reduced_chi2=float(pp.chi2),
            velocity=vel,
            sigma=sig,
            bestfit=np.asarray(pp.bestfit, dtype=float).copy() if keep_full else None,
            weights=np.asarray(pp.weights, dtype=float).copy() if keep_full else None,
            polyweights=(
                None
                if getattr(pp, "polyweights", None) is None or not keep_full
                else np.asarray(pp.polyweights, dtype=float).copy()
            ),
            pp=pp if keep_full else None,
        )
    except Exception as exc:
        return PPXFStateResult(
            success=False,
            chi2_total=np.inf,
            reduced_chi2=np.inf,
            velocity=np.asarray([np.nan]),
            sigma=np.asarray([np.nan]),
            error_message=f"{type(exc).__name__}: {exc}",
        )


def fit_fixed_two_component_state(
    *,
    templates_two_component: np.ndarray,
    component: np.ndarray,
    galaxy: np.ndarray,
    noise: np.ndarray,
    velscale: float,
    lam: np.ndarray,
    lam_temp: np.ndarray,
    goodpixels: np.ndarray,
    velocity_a: float,
    velocity_b: float,
    fraction_a: float,
    start_sigma_a: float,
    start_sigma_b: float,
    sigma_bounds: tuple[float, float],
    degree: int = 4,
    mdegree: int = 0,
    regul: float = 0.0,
    keep_full: bool = False,
) -> PPXFStateResult:
    """Profile one exact ``(V_A, V_B, f_A)`` state.

    ``V_A`` and ``V_B`` are held fixed with pPXF's ``fixed`` keyword.  Only
    the two dispersions and linear nuisance quantities are optimized.  The
    component fraction is imposed with pPXF's ``fraction`` constraint.
    """
    ppxf = _import_ppxf()
    galaxy, noise, lam, goodpixels = _validate_spectrum_inputs(
        galaxy, noise, lam, goodpixels
    )
    templates_two_component = np.asarray(templates_two_component, dtype=float)
    component = np.asarray(component, dtype=int)
    lam_temp = np.asarray(lam_temp, dtype=float)

    if not 0.0 < float(fraction_a) < 1.0:
        raise ValueError("fraction_a must lie strictly between 0 and 1.")
    if templates_two_component.ndim != 2:
        raise ValueError("templates_two_component must be two-dimensional.")
    if templates_two_component.shape[0] != lam_temp.size:
        raise ValueError("Template wavelength axis does not match lam_temp.")
    if component.size != templates_two_component.shape[1]:
        raise ValueError("component must contain one label per template column.")
    if set(np.unique(component)) != {0, 1}:
        raise ValueError("Script-3 two-component fits require component labels {0, 1}.")

    # V_A and V_B are exact likelihood coordinates, not free search parameters.
    # pPXF still requires finite bounds when bounds are supplied for the free
    # dispersions, so give each fixed velocity only a tiny bookkeeping interval
    # around its requested value.  The ``fixed`` keyword is what enforces the
    # exact coordinate.  Keeping these bounds local also prevents pPXF's
    # lam/lam_temp template-coverage logic from interpreting a fixed state as a
    # spurious +/-2000 km/s velocity search.
    velocity_half_width = max(1.0e-3, 1.0e-3 * abs(float(velscale)))
    va = float(velocity_a)
    vb = float(velocity_b)
    bounds = [
        [
            [va - velocity_half_width, va + velocity_half_width],
            [float(sigma_bounds[0]), float(sigma_bounds[1])],
        ],
        [
            [vb - velocity_half_width, vb + velocity_half_width],
            [float(sigma_bounds[0]), float(sigma_bounds[1])],
        ],
    ]

    try:
        pp = ppxf(
            templates_two_component,
            galaxy,
            noise,
            float(velscale),
            start=[
                [float(velocity_a), float(start_sigma_a)],
                [float(velocity_b), float(start_sigma_b)],
            ],
            moments=[2, 2],
            component=component,
            fraction=float(fraction_a),
            fixed=[[True, False], [True, False]],
            bounds=bounds,
            goodpixels=goodpixels,
            degree=int(degree),
            mdegree=int(mdegree),
            regul=float(regul),
            lam=lam,
            lam_temp=lam_temp,
            clean=False,
            quiet=True,
        )
        vel, sig = _extract_solution(pp, 2)
        total = total_chi2_diagonal(galaxy, pp.bestfit, noise, goodpixels)
        return PPXFStateResult(
            success=True,
            chi2_total=total,
            reduced_chi2=float(pp.chi2),
            velocity=vel,
            sigma=sig,
            bestfit=np.asarray(pp.bestfit, dtype=float).copy() if keep_full else None,
            weights=np.asarray(pp.weights, dtype=float).copy() if keep_full else None,
            polyweights=(
                None
                if getattr(pp, "polyweights", None) is None or not keep_full
                else np.asarray(pp.polyweights, dtype=float).copy()
            ),
            pp=pp if keep_full else None,
        )
    except Exception as exc:
        return PPXFStateResult(
            success=False,
            chi2_total=np.inf,
            reduced_chi2=np.inf,
            velocity=np.asarray([np.nan, np.nan]),
            sigma=np.asarray([np.nan, np.nan]),
            error_message=f"{type(exc).__name__}: {exc}",
        )
