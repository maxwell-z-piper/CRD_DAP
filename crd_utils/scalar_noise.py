"""Residual-noise calibration from the exact best state of each RH3 grid.

Script 2 supplies the spatially calibrated noise vector. Script 3 fits its
complete unregularized grid once with that vector, then measures a residual
multiplier from the best grid state. A uniform multiplier changes chi-square
amplitudes, but not the location of the minimum of an unregularized grid.
"""
from __future__ import annotations

import hashlib
import json

import numpy as np


SCHEMA = "crd-rh3-full-grid-residual-noise-v2"


def fingerprint(*arrays, settings):
    """Bind resumed bin checkpoints to every input of the noise experiment."""
    payload = dict(settings, residual_method=SCHEMA, residual_floor=None)
    digest = hashlib.sha256(json.dumps(payload, sort_keys=True).encode())
    for value in arrays:
        arr = np.ascontiguousarray(value)
        digest.update(str(arr.shape).encode())
        digest.update(str(arr.dtype).encode())
        digest.update(arr.tobytes())
    return digest.hexdigest()


def from_best_grid_state(galaxy, bestfit, spatial_noise, goodpixels):
    """Return raw MAD scale, residual center, and normalized residual vector.

    The scale is deliberately allowed below one. The same fixed good pixels
    must have been used for the full grid and its exact best-state refit.
    """
    galaxy = np.asarray(galaxy, dtype=float)
    bestfit = np.asarray(bestfit, dtype=float)
    spatial_noise = np.asarray(spatial_noise, dtype=float)
    gp = np.asarray(goodpixels, dtype=int)
    if galaxy.ndim != 1 or bestfit.shape != galaxy.shape or spatial_noise.shape != galaxy.shape:
        raise ValueError("Best-grid residual arrays must be matching one-dimensional spectra")
    if gp.size == 0 or np.any(gp < 0) or np.any(gp >= galaxy.size):
        raise ValueError("Best-grid residuals require nonempty valid goodpixels")
    if (np.any(~np.isfinite(galaxy[gp])) or np.any(~np.isfinite(bestfit[gp]))
            or np.any(~np.isfinite(spatial_noise[gp])) or np.any(spatial_noise[gp] <= 0)):
        raise ValueError("Best-grid residual inputs are invalid on goodpixels")
    z = (galaxy[gp] - bestfit[gp]) / spatial_noise[gp]
    center = float(np.median(z))
    scale = float(1.4826 * np.median(np.abs(z - center)))
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError("Best-grid residual scatter is non-finite or non-positive")
    residual = np.full(galaxy.shape, np.nan, dtype=np.float32)
    residual[gp] = z
    return scale, center, residual


def rescale_chi2(chi2, scale):
    """Apply one bin-wide sigma multiplier to a diagonal chi-square result."""
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError("Residual scale must be finite and positive")
    return np.asarray(chi2) / float(scale) ** 2


def write_calibration(run, *, scales, centers, residuals, input_hash):
    """Save the all-bin audit product after checkpoint consolidation."""
    scales = np.asarray(scales, dtype=float)
    centers = np.asarray(centers, dtype=float)
    residuals = np.asarray(residuals, dtype=np.float32)
    if scales.ndim != 1 or centers.shape != scales.shape or residuals.shape[0] != scales.size:
        raise ValueError("Residual calibration arrays have incompatible shapes")
    if np.any(~np.isfinite(scales)) or np.any(scales <= 0):
        raise ValueError("Residual calibration scales must be finite and positive")
    product = run.products_dir / "RH3_scalar_noise_calibration.npz"
    metadata = run.metadata_dir / "RH3_scalar_noise_calibration.json"
    np.savez_compressed(product,
        raw_scale=scales,
        adopted_scale=scales,
        residual_center=centers,
        normalized_residual=residuals,
        input_hash=np.asarray(str(input_hash)),
        source=np.asarray("exact best state of full unregularized likelihood grid"),
        floor_applied=np.asarray(False))
    metadata.write_text(json.dumps(dict(
        schema=SCHEMA,
        input_hash=str(input_hash),
        product=str(product),
        floor_applied=False,
        explanation="MAD of exact full-grid best-state residuals; all bins; diagonal pPXF",
    ), indent=2) + "\n")
    return product
