"""Standardized pPXF wrappers for CRD_DAP spectral likelihood work.

All likelihood states for one PowerBin are compared on exactly the same
log-wavelength samples, with the same continuum treatment and the same *frozen*
noise model.  For covariance-aware fits CRD_DAP uses the validated pPXF-9.4.8
``noise_inv_cholesky`` patch: a precomputed ``W=L^-1`` is reused for every
state rather than repeatedly factorizing the same covariance matrix.

The statistically important quantity returned here is total chi-square.  For a
diagonal noise vector this is ``sum((r/sigma)^2)``.  For a cached covariance
whitener it is ``sum((W r)[goodpixels]^2)``.  pPXF's reduced ``pp.chi2`` is kept
only as a QC quantity.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from .covariance import covariance_total_chi2


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
    noise_inv_cholesky: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray | None]:
    """Validate a spectrum and build pPXF-safe *private* input vectors.

    CRD_DAP keeps rejected log-grid samples as NaN. pPXF validates complete input
    vectors before applying ``goodpixels``, so only already-excluded invalid
    samples receive benign local placeholders. Fitted samples are never repaired.

    When a cached inverse-Cholesky matrix is supplied it must match the complete
    pPXF vector. Its fitted rows must not depend on excluded residuals; the
    covariance constructor enforces this invariant when the matrix is built.
    """
    galaxy = np.array(galaxy, dtype=float, copy=True)
    noise = np.array(noise, dtype=float, copy=True)
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

    W = None
    if noise_inv_cholesky is not None:
        W = np.asarray(noise_inv_cholesky, dtype=float)
        if W.shape != (galaxy.size, galaxy.size):
            raise ValueError(
                f"noise_inv_cholesky must have shape {(galaxy.size, galaxy.size)}"
            )
        if not np.all(np.isfinite(W)):
            raise ValueError("noise_inv_cholesky must contain only finite values")
        if not np.allclose(W, np.tril(W), rtol=0.0, atol=1e-12):
            raise ValueError("noise_inv_cholesky must be lower triangular")
        if np.any(np.diag(W) <= 0):
            raise ValueError("noise_inv_cholesky must have a positive diagonal")

    excluded = np.ones(galaxy.size, dtype=bool)
    excluded[goodpixels] = False
    bad_galaxy = excluded & ~np.isfinite(galaxy)
    bad_noise = excluded & (~np.isfinite(noise) | (noise <= 0))

    if np.any(bad_galaxy):
        galaxy_fill = float(np.median(galaxy[goodpixels]))
        if not np.isfinite(galaxy_fill):
            galaxy_fill = 0.0
        galaxy[bad_galaxy] = galaxy_fill

    if np.any(bad_noise):
        noise_fill = float(np.median(noise[goodpixels]))
        if not np.isfinite(noise_fill) or noise_fill <= 0:
            raise ValueError("Could not construct a finite positive pPXF noise placeholder.")
        noise[bad_noise] = noise_fill

    if np.any(~np.isfinite(galaxy)):
        raise ValueError(
            "galaxy contains non-finite samples outside the validated excluded-pixel placeholder case."
        )
    if np.any(~np.isfinite(noise)) or np.any(noise <= 0):
        raise ValueError(
            "noise contains non-finite/non-positive samples outside the validated excluded-pixel placeholder case."
        )

    return galaxy, noise, lam, goodpixels, W


def total_chi2_diagonal(
    galaxy: np.ndarray,
    bestfit: np.ndarray,
    noise: np.ndarray,
    goodpixels: np.ndarray,
) -> float:
    """Return total chi-square for a diagonal spectral-noise model."""
    gp = np.asarray(goodpixels, dtype=int)
    resid = np.asarray(galaxy, dtype=float)[gp] - np.asarray(bestfit, dtype=float)[gp]
    sig = np.asarray(noise, dtype=float)[gp]
    return float(np.sum((resid / sig) ** 2))


def total_chi2(
    galaxy: np.ndarray,
    bestfit: np.ndarray,
    noise: np.ndarray,
    goodpixels: np.ndarray,
    *,
    noise_inv_cholesky: np.ndarray | None = None,
) -> float:
    """Return the exact total chi-square for the adopted diagonal/covariance model."""
    if noise_inv_cholesky is None:
        return total_chi2_diagonal(galaxy, bestfit, noise, goodpixels)
    residual = np.asarray(galaxy, dtype=float) - np.asarray(bestfit, dtype=float)
    return covariance_total_chi2(residual, noise_inv_cholesky, goodpixels)


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


def _ppxf_noise_kwargs(noise_inv_cholesky: np.ndarray | None) -> dict[str, Any]:
    return {} if noise_inv_cholesky is None else {"noise_inv_cholesky": noise_inv_cholesky}


def _success_result(pp, n_components: int, galaxy, noise, goodpixels, W, keep_full) -> PPXFStateResult:
    vel, sig = _extract_solution(pp, n_components)
    total = total_chi2(
        galaxy, pp.bestfit, noise, goodpixels, noise_inv_cholesky=W
    )
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
    noise_inv_cholesky: np.ndarray | None = None,
    keep_full: bool = True,
) -> PPXFStateResult:
    """Run the one-component RH3 control fit on a fixed good-pixel set."""
    ppxf = _import_ppxf()
    galaxy, noise, lam, goodpixels, W = _validate_spectrum_inputs(
        galaxy, noise, lam, goodpixels, noise_inv_cholesky
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
            **_ppxf_noise_kwargs(W),
        )
        return _success_result(pp, 1, galaxy, noise, goodpixels, W, keep_full)
    except Exception as exc:
        return PPXFStateResult(
            success=False,
            chi2_total=np.inf,
            reduced_chi2=np.inf,
            velocity=np.asarray([np.nan]),
            sigma=np.asarray([np.nan]),
            error_message=f"{type(exc).__name__}: {exc}",
        )


def fit_free_two_component_losvd(
    *,
    templates_two_component: np.ndarray,
    component: np.ndarray,
    galaxy: np.ndarray,
    noise: np.ndarray,
    velscale: float,
    lam: np.ndarray,
    lam_temp: np.ndarray,
    goodpixels: np.ndarray,
    start_velocity_a: float,
    start_velocity_b: float,
    start_sigma_a: float,
    start_sigma_b: float,
    velocity_bounds_a: tuple[float, float],
    velocity_bounds_b: tuple[float, float],
    sigma_bounds: tuple[float, float],
    degree: int = 4,
    mdegree: int = 0,
    regul: float = 0.0,
    noise_inv_cholesky: np.ndarray | None = None,
    keep_full: bool = True,
) -> PPXFStateResult:
    """Free two-component fit used only to obtain high-quality calibration residuals.

    Both component velocities and dispersions are nonlinear free parameters and
    the component light ratio is determined by the duplicated template weights.
    The returned component labels have no physical meaning at this stage; only
    the quality of the summed model spectrum matters for covariance calibration.
    """
    ppxf = _import_ppxf()
    galaxy, noise, lam, goodpixels, W = _validate_spectrum_inputs(
        galaxy, noise, lam, goodpixels, noise_inv_cholesky
    )
    templates_two_component = np.asarray(templates_two_component, dtype=float)
    component = np.asarray(component, dtype=int)
    lam_temp = np.asarray(lam_temp, dtype=float)
    if templates_two_component.ndim != 2 or templates_two_component.shape[0] != lam_temp.size:
        raise ValueError("Two-component templates have incompatible dimensions.")
    if component.size != templates_two_component.shape[1] or set(np.unique(component)) != {0, 1}:
        raise ValueError("Two-component calibration requires one component label per template and labels {0,1}.")

    sb = (float(sigma_bounds[0]), float(sigma_bounds[1]))
    bounds = [
        [[float(velocity_bounds_a[0]), float(velocity_bounds_a[1])], [sb[0], sb[1]]],
        [[float(velocity_bounds_b[0]), float(velocity_bounds_b[1])], [sb[0], sb[1]]],
    ]
    try:
        pp = ppxf(
            templates_two_component,
            galaxy,
            noise,
            float(velscale),
            start=[
                [float(start_velocity_a), float(start_sigma_a)],
                [float(start_velocity_b), float(start_sigma_b)],
            ],
            moments=[2, 2],
            component=component,
            bounds=bounds,
            goodpixels=goodpixels,
            degree=int(degree),
            mdegree=int(mdegree),
            regul=float(regul),
            lam=lam,
            lam_temp=lam_temp,
            clean=False,
            quiet=True,
            **_ppxf_noise_kwargs(W),
        )
        return _success_result(pp, 2, galaxy, noise, goodpixels, W, keep_full)
    except Exception as exc:
        return PPXFStateResult(
            success=False,
            chi2_total=np.inf,
            reduced_chi2=np.inf,
            velocity=np.asarray([np.nan, np.nan]),
            sigma=np.asarray([np.nan, np.nan]),
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
    noise_inv_cholesky: np.ndarray | None = None,
    keep_full: bool = False,
) -> PPXFStateResult:
    """Profile one exact ``(V_A,V_B,f_A)`` state under a frozen noise model."""
    ppxf = _import_ppxf()
    galaxy, noise, lam, goodpixels, W = _validate_spectrum_inputs(
        galaxy, noise, lam, goodpixels, noise_inv_cholesky
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

    velocity_half_width = max(1.0e-3, 1.0e-3 * abs(float(velscale)))
    va = float(velocity_a)
    vb = float(velocity_b)
    bounds = [
        [[va - velocity_half_width, va + velocity_half_width], [float(sigma_bounds[0]), float(sigma_bounds[1])]],
        [[vb - velocity_half_width, vb + velocity_half_width], [float(sigma_bounds[0]), float(sigma_bounds[1])]],
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
            **_ppxf_noise_kwargs(W),
        )
        return _success_result(pp, 2, galaxy, noise, goodpixels, W, keep_full)
    except Exception as exc:
        return PPXFStateResult(
            success=False,
            chi2_total=np.inf,
            reduced_chi2=np.inf,
            velocity=np.asarray([np.nan, np.nan]),
            sigma=np.asarray([np.nan, np.nan]),
            error_message=f"{type(exc).__name__}: {exc}",
        )
