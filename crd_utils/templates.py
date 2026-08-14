"""XSL SSP template preparation for RH3 kinematics and BL populations.

The baseline plan uses the full selected XSL SSP grid in RH3 rather than
reducing the library solely to save computation time. Template processing must
preserve age/metallicity metadata, enforce a common normalization convention,
match the measured instrumental LSF, and maintain explicit air/vacuum wavelength
metadata.
"""

from __future__ import annotations


def prepare_xsl_rh3_templates(*args, **kwargs):
    """Crop, resolution-match, normalize, and log-rebin XSL SSPs for RH3."""
    raise NotImplementedError("Implemented with Script 3 after Script-1 LSF products exist.")


def prepare_xsl_bl_templates(*args, **kwargs):
    """Prepare two identical SSP template libraries for BL Disk A/B fitting."""
    raise NotImplementedError("Implemented with Script 6.")
