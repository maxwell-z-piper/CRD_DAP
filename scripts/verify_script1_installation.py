#!/usr/bin/env python3
"""Fast preflight check for a synchronized CRD_DAP Script-1 installation.

Run this before the full real-data Script 1 when files have been replaced or
merged.  It does not open the large KCWI cubes.  It verifies that the installed
helper modules expose the API expected by the current driver and that the key
configuration controls added during real-data QC are present.
"""

from __future__ import annotations

import argparse
import inspect
from pathlib import Path
import sys

import crd_utils as crd
from crd_utils import cube_utils, io, noise, plotting, psf_lsf, validation


def _require(module, names):
    missing = [name for name in names if not hasattr(module, name)]
    if missing:
        raise RuntimeError(
            f"{module.__name__} is missing required Script-1 API: {missing}. "
            "This usually means files from different CRD_DAP revisions were mixed."
        )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, help="Target-specific config file")
    args = parser.parse_args()

    print("CRD_DAP Script-1 preflight")
    print(f"Python: {sys.version.split()[0]}")
    print(f"crd_utils: {Path(crd.__file__).resolve()}")
    print(f"cube_utils: {Path(cube_utils.__file__).resolve()}")
    print(f"plotting: {Path(plotting.__file__).resolve()}")
    print(f"psf_lsf: {Path(psf_lsf.__file__).resolve()}")

    _require(plotting, [
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
    ])
    _require(io, [
        "inspect_fits_extensions", "load_kcwi_cube", "load_kcwikit_stack",
        "save_prepared_cube", "summarize_binary_mask",
        "summarize_effective_exposure", "summarize_integer_flags", "write_json",
    ])
    _require(noise, ["characterize_preliminary_noise", "save_noise_diagnostic"])
    _require(psf_lsf, [
        "estimate_psf", "measure_arc_lsf_from_files", "read_primary_header",
        "save_arc_lsf_result",
    ])
    _require(validation, [
        "sky_separation_arcsec", "validate_arc_science_configuration",
        "validate_script1_conventions",
    ])

    reg_sig = inspect.signature(cube_utils.register_cube_pair)
    for name in (
        "min_contrast_snr",
        "contrast_smooth_sigma_pix",
        "max_residual_shift_arcsec",
    ):
        if name not in reg_sig.parameters:
            raise RuntimeError(
                f"register_cube_pair is missing keyword {name!r}; signature is {reg_sig}"
            )

    cfg = crd.load_config(args.config)
    required_cfg = [
        "BL_EXPECTED_GRATING",
        "RH3_EXPECTED_GRATING",
        "LSF_EDGE_EXTRAPOLATION_WARNING_ANGSTROM",
        "REGISTRATION_USE_COMMON_WAVELENGTH_IF_AVAILABLE",
        "REGISTRATION_MIN_COMMON_RANGE_ANGSTROM",
        "REGISTRATION_MIN_COMMON_CHANNELS",
        "REGISTRATION_MIN_CONTRAST_SNR",
        "REGISTRATION_CONTRAST_SMOOTH_SIGMA_PIX",
        "REGISTRATION_MAX_RESIDUAL_SHIFT_ARCSEC",
    ]
    missing_cfg = [name for name in required_cfg if not hasattr(cfg, name)]
    if missing_cfg:
        raise RuntimeError(f"Target config is missing current Script-1 settings: {missing_cfg}")

    print(f"Config: {Path(args.config).resolve()}")
    print(f"Expected gratings: BL={cfg.BL_EXPECTED_GRATING}, RH3-stream={cfg.RH3_EXPECTED_GRATING}")
    print(f"register_cube_pair signature: {reg_sig}")
    print("PREFLIGHT PASS: Script-1 driver/helper API is synchronized.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
