"""API-contract tests for CRD_DAP Script 1.

These tests exist specifically to prevent partial file updates from leaving the
Script-1 driver and helper modules out of sync.  They are intentionally simple:
if the driver calls a public helper, the helper must exist with the keyword
arguments the driver expects.
"""

from __future__ import annotations

import inspect

from crd_utils import cube_utils, io, noise, plotting, psf_lsf, validation


def test_script1_plotting_api_is_complete():
    required = {
        "plot_bad_wavelength_fraction",
        "plot_center_comparison",
        "plot_collapsed_continuum",
        "plot_effective_exposure",
        "plot_lsf",
        "plot_lsf_spatial_summary",
        "plot_lsf_spatial_variation",
        "plot_normalized_residuals",
        "plot_psf_comparison",
        "plot_psf_summary",
        "plot_registration",
        "plot_spectral_covariance",
        "plot_valid_spaxels",
    }
    missing = sorted(name for name in required if not hasattr(plotting, name))
    assert not missing, f"Missing Script-1 plotting helpers: {missing}"


def test_registration_api_accepts_local_search_radius():
    sig = inspect.signature(cube_utils.register_cube_pair)
    assert "max_residual_shift_arcsec" in sig.parameters
    assert "min_contrast_snr" in sig.parameters
    assert "contrast_smooth_sigma_pix" in sig.parameters


def test_script1_io_api_is_complete():
    required = {
        "inspect_fits_extensions",
        "load_kcwi_cube",
        "load_kcwikit_stack",
        "save_prepared_cube",
        "summarize_binary_mask",
        "summarize_effective_exposure",
        "summarize_integer_flags",
        "write_json",
    }
    missing = sorted(name for name in required if not hasattr(io, name))
    assert not missing, f"Missing Script-1 I/O helpers: {missing}"


def test_script1_noise_api_is_complete():
    required = {"characterize_preliminary_noise", "save_noise_diagnostic"}
    missing = sorted(name for name in required if not hasattr(noise, name))
    assert not missing, f"Missing Script-1 noise helpers: {missing}"


def test_script1_lsf_api_is_complete():
    required = {
        "estimate_psf",
        "measure_arc_lsf_from_files",
        "read_primary_header",
        "save_arc_lsf_result",
    }
    missing = sorted(name for name in required if not hasattr(psf_lsf, name))
    assert not missing, f"Missing Script-1 PSF/LSF helpers: {missing}"


def test_script1_validation_api_is_complete():
    required = {
        "sky_separation_arcsec",
        "validate_arc_science_configuration",
        "validate_script1_conventions",
    }
    missing = sorted(name for name in required if not hasattr(validation, name))
    assert not missing, f"Missing Script-1 validation helpers: {missing}"
