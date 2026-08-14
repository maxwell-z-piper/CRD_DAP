"""Cube preparation, collapsed-continuum images, registration, and masks.

The final implementation should distinguish hard data-quality masks from
analysis-specific wavelength masks. Low S/N alone is not a reason to reject a
spaxel; PowerBin is responsible for combining low-S/N spatial samples.
"""

from __future__ import annotations

import numpy as np


def collapsed_continuum(flux: np.ndarray, good: np.ndarray | None = None, *, statistic: str = "median") -> np.ndarray:
    """Collapse a spectral cube into a robust 2-D continuum image.

    Parameters
    ----------
    flux
        Array with wavelength on the final axis.
    good
        Optional boolean mask with the same shape as ``flux``. True means the
        pixel is usable.
    statistic
        ``'median'`` or ``'mean'``. The production implementation may also add
        inverse-variance weighting once Script 1's uncertainty model is known.
    """
    f = np.asarray(flux, dtype=float)
    if f.ndim != 3:
        raise ValueError("flux must be a 3-D cube with wavelength on the final axis.")
    if good is not None:
        good = np.asarray(good, dtype=bool)
        if good.shape != f.shape:
            raise ValueError("good mask must match flux shape.")
        f = np.where(good, f, np.nan)
    if statistic == "median":
        return np.nanmedian(f, axis=-1)
    if statistic == "mean":
        return np.nanmean(f, axis=-1)
    raise ValueError("statistic must be 'median' or 'mean'.")


def brightest_spaxel(image: np.ndarray) -> tuple[int, int]:
    """Return the brightest finite *positive* image pixel as an initializer."""
    image = np.asarray(image, dtype=float)
    valid = np.isfinite(image) & (image > 0)
    if not np.any(valid):
        raise ValueError("Image contains no finite positive pixels.")
    masked = np.where(valid, image, -np.inf)
    return tuple(int(i) for i in np.unravel_index(np.argmax(masked), image.shape))


def register_cube_pair(*args, **kwargs):
    """Register BL and RH3 cubes onto a common sky-coordinate system.

    The production implementation will use WCS/sky coordinates and collapsed
    continuum structure, not naive array-index equality when the two cubes have
    different WCS sampling.
    """
    raise NotImplementedError("Implemented with Script 1.")
