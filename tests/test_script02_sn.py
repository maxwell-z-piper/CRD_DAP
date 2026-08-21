import numpy as np

from crd_utils import binning


def _spec(flux, unc, good=None):
    flux = np.asarray(flux, dtype=float)
    unc = np.asarray(unc, dtype=float)
    if good is None:
        good = np.ones_like(flux, dtype=bool)
    return binning.CoaddedBinSpectra(
        flux=flux,
        uncertainty=unc,
        good=np.asarray(good, dtype=bool),
        contributing_spaxels=np.ones_like(flux, dtype=np.int16),
        n_members=np.ones(flux.shape[0], dtype=np.int32),
        spatial_scale_factor=1.0,
        spatial_scale_reason="test",
    )


def test_achieved_sn_uses_ratio_of_medians_for_positive_continuum():
    spec = _spec(
        [[10.0, 11.0, 9.0, 10.0, 10.5]],
        [[2.0, 2.0, 2.0, 2.0, 2.0]],
    )
    wave = np.arange(5, dtype=float) + 5000.0
    d = binning.achieved_sn_diagnostics_per_bin(
        spec,
        wave,
        rest_range=(4999.0, 5005.0),
        redshift=0.0,
        min_good_channels=2,
    )
    assert np.isclose(d.sn[0], 5.0)
    assert np.isclose(d.signed_sn[0], 5.0)
    assert d.positive_continuum[0]


def test_nonpositive_continuum_is_nan_but_signed_diagnostic_is_preserved():
    # Mimics the symptom from the real RL test: a negative continuum combined
    # with a very small formal uncertainty.  The production-facing achieved S/N
    # must be undefined, not a giant negative number.
    spec = _spec(
        [[-1.0, -1.1, -0.9, -1.0, -1.05]],
        [[1e-5, 1e-5, 1e-5, 1e-5, 1e-5]],
    )
    wave = np.arange(5, dtype=float) + 8800.0
    d = binning.achieved_sn_diagnostics_per_bin(
        spec,
        wave,
        rest_range=(8799.0, 8805.0),
        redshift=0.0,
        min_good_channels=2,
        require_positive_continuum=True,
    )
    assert np.isnan(d.sn[0])
    assert d.signed_sn[0] < -1e4
    assert d.legacy_median_ratio[0] < -1e4
    assert not d.positive_continuum[0]


def test_too_few_good_channels_returns_nan():
    spec = _spec([[10.0, 10.0, 10.0]], [[2.0, 2.0, 2.0]], [[True, False, False]])
    wave = np.array([5000.0, 5001.0, 5002.0])
    d = binning.achieved_sn_diagnostics_per_bin(
        spec,
        wave,
        rest_range=(4999.0, 5003.0),
        redshift=0.0,
        min_good_channels=2,
    )
    assert np.isnan(d.sn[0])
    assert d.n_good_channels[0] == 1


def test_compatibility_wrapper_matches_diagnostics_sn():
    spec = _spec([[6.0, 6.0, 6.0]], [[2.0, 2.0, 2.0]])
    wave = np.array([5000.0, 5001.0, 5002.0])
    kwargs = dict(rest_range=(4999.0, 5003.0), redshift=0.0, min_good_channels=2)
    d = binning.achieved_sn_diagnostics_per_bin(spec, wave, **kwargs)
    sn = binning.achieved_sn_per_bin(spec, wave, **kwargs)
    assert np.allclose(sn, d.sn, equal_nan=True)
