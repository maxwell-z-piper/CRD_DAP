"""Empirical PSF and instrumental LSF characterization.

The production pipeline requires the KCWI DRP master arc calibration products
and measures the instrumental line-spread function from unresolved arc lines.
The preferred model can depend on wavelength and, if supported by the data,
slice/spatial position. A nominal resolving power is a fallback diagnostic, not
the primary production LSF.
"""

from __future__ import annotations

import numpy as np


def gaussian_fwhm_to_sigma(fwhm: np.ndarray | float) -> np.ndarray:
    """Convert Gaussian FWHM to sigma in the same units."""
    return np.asarray(fwhm, dtype=float) / np.sqrt(8.0 * np.log(2.0))


def required_template_convolution_sigma(data_sigma: np.ndarray, template_sigma: np.ndarray) -> np.ndarray:
    r"""Return :math:`\sqrt{\sigma_\mathrm{data}^2-\sigma_\mathrm{temp}^2}`.

    Negative values indicate that the templates are broader than the data and
    cannot be matched by further smoothing the templates. The production code
    should treat this as an explicit compatibility warning/error rather than
    silently taking an absolute value.
    """
    data = np.asarray(data_sigma, dtype=float)
    template = np.asarray(template_sigma, dtype=float)
    diff2 = data**2 - template**2
    result = np.full(np.broadcast(data, template).shape, np.nan, dtype=float)
    good = diff2 >= 0
    result[good] = np.sqrt(diff2[good])
    return result


def measure_arc_lsf(*args, **kwargs):
    """Fit unresolved master-arc lines and infer FWHM(lambda, slice).

    Implemented with Script 1 after the exact master-arc product structure and
    usable arc-line list have been inspected.
    """
    raise NotImplementedError("Implemented with Script 1 using real master arcs.")


def measure_cube_psf(*args, **kwargs):
    """Estimate the effective delivered PSF of a stacked cube when possible."""
    raise NotImplementedError("Implemented with Script 1 after data inspection.")
