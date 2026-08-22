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
    "BINNING_APERTURE_WARNING": QualityFlagDefinition(
        "BINNING_APERTURE_WARNING",
        "The automatically detected stellar-body aperture is unusually small or requires manual scientific inspection before accepting the PowerBin tessellation.",
    ),
    "BIN_TRANSFER_INCOMPLETE": QualityFlagDefinition(
        "BIN_TRANSFER_INCOMPLETE",
        "The red/RH3 native grid does not recover the configured fraction of BL-defined physical bin membership through the celestial WCS transfer.",
    ),
    "LOW_BL_BIN_SN": QualityFlagDefinition(
        "LOW_BL_BIN_SN",
        "One or more final BL PowerBins have measured continuum S/N substantially below the configured target despite the PowerBin capacity target.",
    ),
    "SPATIAL_COVARIANCE_UNCALIBRATED": QualityFlagDefinition(
        "SPATIAL_COVARIANCE_UNCALIBRATED",
        "PowerBin capacities and formal coadded uncertainties currently use diagonal spatial variance because no validated KcwiKit spatial covariance law has been supplied.",
    ),
    "SN_WINDOW_COVERAGE_WARNING": QualityFlagDefinition(
        "SN_WINDOW_COVERAGE_WARNING",
        "A configured Script-2 continuum-S/N window extends beyond the Script-1 GOODWAVE envelope; the requested window is retained and the truncation is recorded as QC rather than automatically replaced.",
    ),
    "NONPOSITIVE_BIN_CONTINUUM": QualityFlagDefinition(
        "NONPOSITIVE_BIN_CONTINUUM",
        "One or more bin S/N windows have non-positive median continuum, so a positive achieved-S/N value is undefined; signed diagnostics are retained for inspection.",
    ),
    "EXTREME_BIN_SN_DIAGNOSTIC": QualityFlagDefinition(
        "EXTREME_BIN_SN_DIAGNOSTIC",
        "One or more signed or legacy bin S/N diagnostics exceed the configured extreme-value warning threshold and should be inspected for tiny formal uncertainties or edge/systematic residuals.",
    ),
    "BIN_SN_ESTIMATOR_DISAGREEMENT": QualityFlagDefinition(
        "BIN_SN_ESTIMATOR_DISAGREEMENT",
        "The robust ratio-of-medians S/N and legacy median(flux/uncertainty) estimator disagree by more than the configured factor for one or more bins.",
    ),
}


def describe_flag(name: str) -> str:
    """Return the human-readable meaning of a standardized quality flag."""
    try:
        return FLAG_DEFINITIONS[name].description
    except KeyError as exc:
        raise KeyError(f"Unknown CRD_DAP quality flag: {name}") from exc
