import inspect
from pathlib import Path

from crd_utils import binning, io, plotting


def test_script02_support_api_is_complete():
    for name in (
        "continuum_window_maps",
        "make_analysis_aperture",
        "run_powerbin",
        "transfer_bin_map_by_wcs",
        "coadd_bin_spectra",
        "achieved_sn_per_bin",
        "normalized_flux_weights",
        "bin_centroids",
    ):
        assert hasattr(binning, name), name
    assert hasattr(io, "load_prepared_cube")
    for name in (
        "plot_binning_aperture",
        "plot_master_bins",
        "plot_bin_value_map",
        "plot_bl_rh3_sn_comparison",
        "plot_bin_transfer",
        # Guard against dropping existing Script-1 plotting helpers.
        "plot_effective_exposure",
        "plot_registration",
        "plot_lsf",
    ):
        assert hasattr(plotting, name), name


def test_script02_driver_helper_keywords_match():
    assert "max_distance_arcsec" in inspect.signature(binning.transfer_bin_map_by_wcs).parameters
    assert "min_member_fraction" in inspect.signature(binning.coadd_bin_spectra).parameters


def test_script02_driver_exists():
    repo_root = Path(__file__).resolve().parents[1]
    assert (repo_root / "scripts" / "02_make_master_BL_bins.py").is_file()
