"""Public convenience namespace for CRD_DAP helper functions.

Science-driver scripts may use ``import crd_utils as crd`` for common public
functions while implementation details remain split into focused modules.
"""

from .config import PipelineConfig, load_config, validate_config, snapshot_config, write_config_manifest
from .logging_utils import RunContext, create_run_context, setup_pipeline_logger, setup_step_logger, log_section
from .likelihood import (
    delta_chi2,
    relative_likelihood_from_delta_chi2,
    normalize_weights,
    likelihood_mass_cell_counts,
    effective_number_of_states,
    descend_to_local_minimum,
    basin_mask_for_anchor,
)
from .geometry import inclination_from_axis_ratio, sample_inclination_prior, deprojected_radius_azimuth, signed_pa_coordinate
from .disk_model import RingGrid, make_ring_grid, interpolate_rotation_curve, projected_circular_velocity, flux_weighted_bin_velocity
from .sampling import sample_discrete_states, effective_sample_size, posterior_quantiles, check_quantile_convergence
from .runtime import TimingResult, timed_block, estimate_parallel_walltime
from .quality import FLAG_DEFINITIONS, describe_flag

__all__ = [
    "PipelineConfig",
    "load_config",
    "validate_config",
    "snapshot_config",
    "write_config_manifest",
    "RunContext",
    "create_run_context",
    "setup_pipeline_logger",
    "setup_step_logger",
    "log_section",
    "delta_chi2",
    "relative_likelihood_from_delta_chi2",
    "normalize_weights",
    "likelihood_mass_cell_counts",
    "effective_number_of_states",
    "descend_to_local_minimum",
    "basin_mask_for_anchor",
    "inclination_from_axis_ratio",
    "sample_inclination_prior",
    "deprojected_radius_azimuth",
    "signed_pa_coordinate",
    "RingGrid",
    "make_ring_grid",
    "interpolate_rotation_curve",
    "projected_circular_velocity",
    "flux_weighted_bin_velocity",
    "sample_discrete_states",
    "effective_sample_size",
    "posterior_quantiles",
    "check_quantile_convergence",
    "TimingResult",
    "timed_block",
    "estimate_parallel_walltime",
    "FLAG_DEFINITIONS",
    "describe_flag",
]
