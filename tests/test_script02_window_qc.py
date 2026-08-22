import numpy as np

from crd_utils import binning, quality


def test_sn_window_coverage_detects_red_edge_truncation():
    wave = np.arange(8800.0, 8931.0, 1.0)
    good = np.ones(wave.size, dtype=bool)
    result = binning.sn_window_coverage(
        wave,
        good,
        rest_range=(8820.0, 9060.0),
        redshift=0.0,
    )
    assert not result.truncated_blue
    assert result.truncated_red
    assert result.usable_observed_max == 8930.0
    assert 0.0 < result.envelope_coverage_fraction < 1.0


def test_sn_window_coverage_does_not_treat_internal_mask_as_edge_truncation():
    wave = np.arange(5000.0, 5101.0, 1.0)
    good = np.ones(wave.size, dtype=bool)
    good[45:56] = False
    result = binning.sn_window_coverage(
        wave,
        good,
        rest_range=(5010.0, 5090.0),
        redshift=0.0,
    )
    assert not result.truncated_blue
    assert not result.truncated_red
    assert np.isclose(result.envelope_coverage_fraction, 1.0)
    assert result.usable_channel_fraction < 1.0


def test_sn_window_coverage_flag_is_standardized():
    text = quality.describe_flag("SN_WINDOW_COVERAGE_WARNING")
    assert "GOODWAVE" in text
