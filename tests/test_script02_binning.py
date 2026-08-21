import sys
import types

import numpy as np

from crd_utils import binning


def test_observed_range_from_rest():
    lo, hi = binning.observed_range_from_rest((4800.0, 5400.0), 0.04)
    assert np.isclose(lo, 4992.0)
    assert np.isclose(hi, 5616.0)


def test_bin_sn_diagonal_and_explicit_log10_correction():
    signal = np.array([4.0, 4.0, 4.0, 4.0])
    noise = np.ones(4)
    idx = np.arange(4)
    assert np.isclose(binning.bin_sn(idx, signal, noise), 8.0)
    corrected = binning.bin_sn(
        idx,
        signal,
        noise,
        covariance_mode="log10",
        covariance_alpha=1.0,
    )
    assert corrected < 8.0


def test_auto_aperture_keeps_center_component_and_not_remote_island():
    signal = np.zeros((21, 21), dtype=float)
    noise = np.ones_like(signal)
    good = np.ones_like(signal, dtype=bool)
    signal[8:13, 8:13] = 5.0
    signal[1:4, 1:4] = 9.0

    result = binning.make_analysis_aperture(
        signal,
        noise,
        good,
        center_yx=(10.0, 10.0),
        smooth_sigma_pix=0.5,
        threshold=1.0,
        dilate_pix=0,
        center_max_distance_pix=3.0,
        min_pixels=5,
    )
    assert result.mask[10, 10]
    assert not result.mask[2, 2]



def test_smoothed_aperture_can_detect_extended_light_below_native_pixel_threshold():
    signal = np.zeros((31, 31), dtype=float)
    noise = np.ones_like(signal)
    good = np.ones_like(signal, dtype=bool)
    # Every native source pixel has S/N=0.8 < the aperture threshold of 2, but
    # the spatially extended source becomes significant when neighboring pixels
    # are combined by the Gaussian detection kernel.
    signal[12:19, 12:19] = 0.8

    result = binning.make_analysis_aperture(
        signal,
        noise,
        good,
        center_yx=(15.0, 15.0),
        smooth_sigma_pix=1.5,
        threshold=2.0,
        dilate_pix=0,
        center_max_distance_pix=2.0,
        min_pixels=5,
    )
    assert result.mask[15, 15]
    assert result.significance_proxy[15, 15] > 2.0

def test_run_powerbin_wrapper_uses_callable_capacity(monkeypatch):
    class FakePowerBin:
        def __init__(self, xy, capacity_spec, **kwargs):
            assert np.isfinite(capacity_spec(np.array([0, 1])))
            self.bin_num = np.array([0, 0, 1, 1])
            self.xybin = np.array([[0.5, 0.0], [2.5, 0.0]])
            self.rbin = np.array([0.5, 0.5])

    fake = types.ModuleType("powerbin")
    fake.PowerBin = FakePowerBin
    monkeypatch.setitem(sys.modules, "powerbin", fake)

    x = np.arange(4, dtype=float)[None, :]
    y = np.zeros_like(x)
    signal = np.full_like(x, 5.0)
    noise = np.ones_like(x)
    aperture = np.ones_like(x, dtype=bool)
    result = binning.run_powerbin(
        x,
        y,
        signal,
        noise,
        aperture,
        target_sn=7.0,
        pixel_size_arcsec=1.0,
        verbose=0,
    )
    assert result.n_bins == 2
    assert result.bin_map.tolist() == [[0, 0, 1, 1]]
