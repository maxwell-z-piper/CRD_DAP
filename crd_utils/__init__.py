"""Public convenience namespace for CRD_DAP helper functions.

Science-driver scripts may use ``import crd_utils as crd`` for common public
functions while implementation details remain split into focused modules.
"""

from .config import PipelineConfig, load_config, validate_config, validate_input_paths, snapshot_config, write_config_manifest
from .logging_utils import RunContext, create_run_context, setup_pipeline_logger, setup_step_logger, log_section
from .io import (
    KCWICube,
    load_kcwi_cube,
    save_prepared_cube,
    inspect_fits_extensions,
    summarize_integer_flags,
    wavelength_axis_from_header,
    discover_arc_sidecars,
    write_json,
    read_json,
)
from .cube_utils import (
    CenterEstimate,
    RegistrationResult,
    collapsed_continuum,
    smooth_image,
    brightest_spaxel,
    flux_weighted_centroid,
    estimate_continuum_center,
    spatial_offset_grids_arcsec,
    register_cube_pair,
)
from .psf_lsf import (
    ArcLSFResult,
    PSFEstimate,
    gaussian_fwhm_to_sigma,
    gaussian_sigma_to_fwhm,
    required_template_convolution_sigma,
    measure_arc_lsf,
    measure_arc_lsf_from_files,
    save_arc_lsf_result,
    estimate_psf,
)
from .noise import (
    NoiseDiagnosticResult,
    robust_variance_scale,
    characterize_preliminary_noise,
    estimate_spectral_correlation,
    generate_correlated_noise,
)
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
    "KCWICube",
    "load_kcwi_cube",
    "save_prepared_cube",
    "inspect_fits_extensions",
    "wavelength_axis_from_header",
    "discover_arc_sidecars",
    "write_json",
    "read_json",
    "CenterEstimate",
    "RegistrationResult",
    "collapsed_continuum",
    "smooth_image",
    "brightest_spaxel",
    "flux_weighted_centroid",
    "estimate_continuum_center",
    "spatial_offset_grids_arcsec",
    "register_cube_pair",
    "ArcLSFResult",
    "PSFEstimate",
    "gaussian_fwhm_to_sigma",
    "gaussian_sigma_to_fwhm",
    "required_template_convolution_sigma",
    "measure_arc_lsf",
    "measure_arc_lsf_from_files",
    "save_arc_lsf_result",
    "estimate_psf",
    "NoiseDiagnosticResult",
    "robust_variance_scale",
    "characterize_preliminary_noise",
    "estimate_spectral_correlation",
    "generate_correlated_noise",
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
