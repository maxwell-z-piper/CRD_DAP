"""BL-defined master PowerBins and exact transfer of spatial membership to RH3."""

from __future__ import annotations

import numpy as np


def shared_sn_color_limits(bl_sn: np.ndarray, rh3_sn: np.ndarray, *, upper_percentile: float = 95.0, lower: float | None = None) -> tuple[float, float]:
    """Return one common color scale for side-by-side BL/RH3 S/N maps."""
    values = np.concatenate([np.asarray(bl_sn, dtype=float).ravel(), np.asarray(rh3_sn, dtype=float).ravel()])
    values = values[np.isfinite(values)]
    if values.size == 0:
        raise ValueError("No finite S/N values supplied.")
    vmax = float(np.percentile(values, upper_percentile))
    vmin = float(np.min(values) if lower is None else lower)
    if vmax <= vmin:
        vmax = float(np.max(values))
    return vmin, vmax


def apply_bin_membership(*args, **kwargs):
    """Coadd another cube using a previously defined physical bin membership."""
    raise NotImplementedError("Implemented with Script 2 after registration data model is fixed.")
