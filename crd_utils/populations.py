"""Two-component BL stellar-population fitting and population summaries.

Disk A and Disk B receive independent full SSP template libraries. Their RH3
kinematics are either fixed to the exact global RH3 solution (primary anchored
fit) or propagated through the RH3-supported state family (likelihood-propagated
fit). The BL stellar light fraction is wavelength-band specific and must not be
forced to equal the RH3 fraction.
"""

from __future__ import annotations


def fit_bl_fraction_grid(*args, **kwargs):
    """Evaluate BL population fits on an explicit stellar f_A grid."""
    raise NotImplementedError("Implemented with Script 6.")


def summarize_ssp_weights(*args, **kwargs):
    """Calculate luminosity/mass-weighted age, metallicity, M/L, and SFH summaries."""
    raise NotImplementedError("Implemented with Script 6 after template metadata is fixed.")
