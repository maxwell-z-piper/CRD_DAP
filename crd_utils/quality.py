"""Central definitions for quality-control flags.

Science scripts should not invent ad-hoc strings for warnings. Keeping all flag
names and descriptions here ensures that maps, tables, logs, and the final MRT
use the same vocabulary.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class QualityFlagDefinition:
    name: str
    description: str


FLAG_DEFINITIONS = {
    "GRID_EDGE_WARNING": QualityFlagDefinition(
        "GRID_EDGE_WARNING",
        "Selected RH3 solution lies within the configured number of grid cells of a parameter-grid edge.",
    ),
    "LOW_TWO_COMPONENT_SUPPORT": QualityFlagDefinition(
        "LOW_TWO_COMPONENT_SUPPORT",
        "Mock-calibrated one-vs-two-component comparison indicates weak or unreliable two-component support.",
    ),
    "OUTSIDE_TWO_SIGMA_REGION": QualityFlagDefinition(
        "OUTSIDE_TWO_SIGMA_REGION",
        "PowerBin lies outside the primary 2sigma-limited radial model but may still be used in the full-aperture robustness run.",
    ),
    "HIGH_BIN_SHEAR": QualityFlagDefinition(
        "HIGH_BIN_SHEAR",
        "Model-predicted intra-bin velocity shear is large enough that centroid-only treatment may bias the LOSVD.",
    ),
    "POOR_RING_COVERAGE": QualityFlagDefinition(
        "POOR_RING_COVERAGE",
        "A non-parametric rotation-curve ring has inadequate spatial/azimuthal data coverage.",
    ),
    "BL_JOINT_NOT_CONVERGED": QualityFlagDefinition(
        "BL_JOINT_NOT_CONVERGED",
        "Likelihood-propagated BL population sampling reached the maximum draw count without satisfying convergence criteria.",
    ),
    "MC_NOT_CONVERGED": QualityFlagDefinition(
        "MC_NOT_CONVERGED",
        "End-to-end Monte Carlo uncertainty propagation reached its maximum realization count without convergence.",
    ),
    "AGE_GRID_EDGE_A": QualityFlagDefinition(
        "AGE_GRID_EDGE_A", "Disk A population weights accumulate at the age boundary of the SSP template grid."
    ),
    "AGE_GRID_EDGE_B": QualityFlagDefinition(
        "AGE_GRID_EDGE_B", "Disk B population weights accumulate at the age boundary of the SSP template grid."
    ),
    "METALLICITY_GRID_EDGE_A": QualityFlagDefinition(
        "METALLICITY_GRID_EDGE_A", "Disk A population weights accumulate at the metallicity boundary of the SSP grid."
    ),
    "METALLICITY_GRID_EDGE_B": QualityFlagDefinition(
        "METALLICITY_GRID_EDGE_B", "Disk B population weights accumulate at the metallicity boundary of the SSP grid."
    ),
    "RADIAL_EXTENT_DISAGREEMENT": QualityFlagDefinition(
        "RADIAL_EXTENT_DISAGREEMENT",
        "2sigma-limited and full-aperture global models give materially different solutions over their shared radial domain.",
    ),
    "AMBIGUOUS_DISK_LABEL": QualityFlagDefinition(
        "AMBIGUOUS_DISK_LABEL",
        "The signed-PA outer-branch convention could not assign Disk A/B unambiguously and requires inspection.",
    ),
    "WAVELENGTH_CONVENTION_ERROR": QualityFlagDefinition(
        "WAVELENGTH_CONVENTION_ERROR",
        "Science and template wavelength conventions or velocity frames are incompatible.",
    ),
    "CENTER_DISAGREEMENT": QualityFlagDefinition(
        "CENTER_DISAGREEMENT",
        "Independent continuum center estimates disagree by more than the configured tolerance.",
    ),
    "REGISTRATION_OFFSET_WARNING": QualityFlagDefinition(
        "REGISTRATION_OFFSET_WARNING",
        "Residual BL/RH3 cross-correlation shift after WCS reprojection exceeds the configured tolerance.",
    ),
    "PSF_NOT_MEASURED": QualityFlagDefinition(
        "PSF_NOT_MEASURED",
        "No explicit or header-based PSF FWHM was available for this arm; a valid PSF must be supplied before PSF-dependent modeling.",
    ),
    "LSF_SPATIAL_VARIATION": QualityFlagDefinition(
        "LSF_SPATIAL_VARIATION",
        "Master-arc line widths show spatial/slice variation larger than the configured fractional tolerance.",
    ),
    "NOISE_SCALE_WARNING": QualityFlagDefinition(
        "NOISE_SCALE_WARNING",
        "Preliminary normalized residuals imply a variance scale materially different from unity; final likelihood work requires follow-up noise calibration.",
    ),
    "SPECTRAL_COVARIANCE_WARNING": QualityFlagDefinition(
        "SPECTRAL_COVARIANCE_WARNING",
        "Preliminary high-pass residuals show non-negligible wavelength-pixel correlation that must be revisited with pPXF residuals.",
    ),
}


def describe_flag(name: str) -> str:
    """Return the human-readable meaning of a standardized quality flag."""
    try:
        return FLAG_DEFINITIONS[name].description
    except KeyError as exc:
        raise KeyError(f"Unknown CRD_DAP quality flag: {name}") from exc
