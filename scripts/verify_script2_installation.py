#!/usr/bin/env python3
"""Fast API/dependency preflight for CRD_DAP Script 2.

This check intentionally does not open the large KCWI FITS cubes.  It verifies
that the active Python environment imports the expected repository files, that
PowerBin is available, and that the Script-2 helper APIs/configuration keys are
present before a real run is attempted.
"""

from __future__ import annotations

import argparse
import inspect
from importlib import metadata as importlib_metadata
from pathlib import Path
import sys


def _parser():
    p = argparse.ArgumentParser(description="Verify CRD_DAP Script-2 installation/API.")
    p.add_argument("--config", required=True)
    return p


def main() -> int:
    args = _parser().parse_args()
    errors: list[str] = []

    try:
        from powerbin import PowerBin
        try:
            powerbin_version = importlib_metadata.version("powerbin")
        except importlib_metadata.PackageNotFoundError:
            powerbin_version = "unknown"
        print(f"PowerBin: PASS | version={powerbin_version} | {PowerBin}")
        pb_params = inspect.signature(PowerBin).parameters
        for required in ("target_capacity", "pixelsize", "regul", "maxiter", "verbose"):
            if required not in pb_params:
                errors.append(f"Installed PowerBin API lacks constructor parameter {required!r}")
    except Exception as exc:
        errors.append(
            "PowerBin import failed. Install into this exact environment with "
            "`python -m pip install powerbin`. Original error: " + repr(exc)
        )

    import crd_utils
    from crd_utils import binning, io, plotting

    print(f"Python:      {sys.executable}")
    print(f"crd_utils:   {Path(crd_utils.__file__).resolve()}")
    print(f"binning.py:  {Path(binning.__file__).resolve()}")
    print(f"io.py:       {Path(io.__file__).resolve()}")
    print(f"plotting.py: {Path(plotting.__file__).resolve()}")

    required_functions = {
        binning: [
            "continuum_window_maps",
            "make_analysis_aperture",
            "run_powerbin",
            "transfer_bin_map_by_wcs",
            "coadd_bin_spectra",
            "achieved_sn_per_bin",
            "normalized_flux_weights",
            "bin_centroids",
        ],
        io: ["load_prepared_cube"],
        plotting: [
            "plot_binning_aperture",
            "plot_master_bins",
            "plot_bin_value_map",
            "plot_bl_rh3_sn_comparison",
            "plot_bin_transfer",
            # Script-1 functions that must remain after Script-2 additions.
            "plot_effective_exposure",
            "plot_registration",
            "plot_lsf",
        ],
    }
    for module, names in required_functions.items():
        for name in names:
            if not hasattr(module, name):
                errors.append(f"{module.__name__}.{name} is missing")

    # API-contract checks guard against the mixed-revision problem encountered
    # during Script-1 development.
    if "max_distance_arcsec" not in inspect.signature(binning.transfer_bin_map_by_wcs).parameters:
        errors.append("transfer_bin_map_by_wcs lacks max_distance_arcsec")
    if "min_member_fraction" not in inspect.signature(binning.coadd_bin_spectra).parameters:
        errors.append("coadd_bin_spectra lacks min_member_fraction")

    cfg = crd_utils.load_config(args.config, validate=True, strict_paths=False)
    required_cfg = [
        "BL_TARGET_SN",
        "BL_BINNING_REST_RANGE_ANGSTROM",
        "RH3_SN_REST_RANGE_ANGSTROM",
        "BINNING_APERTURE_MODE",
        "BINNING_APERTURE_SN_THRESHOLD",
        "POWERBIN_SPATIAL_COVARIANCE_MODE",
        "BIN_TRANSFER_MAX_DISTANCE_ARCSEC",
        "BIN_SPECTRUM_MIN_MEMBER_FRACTION",
    ]
    for key in required_cfg:
        if not hasattr(cfg, key):
            errors.append(f"Config is missing {key}")

    if errors:
        print("\nPREFLIGHT FAIL")
        for error in errors:
            print(f"  - {error}")
        return 1

    print("\nPREFLIGHT PASS: Script-2 driver/helper API is synchronized.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
