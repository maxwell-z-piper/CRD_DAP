import numpy as np

from crd_utils.cube_utils import registration_contrast_snr, residual_registration_shift
from crd_utils import plotting


def test_registration_contrast_rejects_flat_image():
    image = np.ones((40, 50), dtype=float)
    mask = np.ones_like(image, dtype=bool)
    assert registration_contrast_snr(image, mask) == 0.0


def test_registration_contrast_detects_compact_source():
    yy, xx = np.indices((60, 70), dtype=float)
    source = 20.0 * np.exp(-0.5 * (((xx - 35.0) / 4.0) ** 2 + ((yy - 30.0) / 4.0) ** 2))
    rng = np.random.default_rng(12)
    image = source + rng.normal(0.0, 0.25, source.shape)
    mask = np.ones_like(image, dtype=bool)
    assert registration_contrast_snr(image, mask, smooth_sigma_pix=1.0) > 5.0


def test_registration_search_window_bounds_runaway_peak():
    yy, xx = np.indices((80, 90), dtype=float)
    reference = np.exp(-0.5 * (((xx - 45.0) / 3.0) ** 2 + ((yy - 40.0) / 3.0) ** 2))
    # Deliberately put the same compact source far outside a 4-pixel residual
    # search window. The returned local candidate must remain inside the window
    # rather than jumping tens of pixels.
    moving = np.exp(-0.5 * (((xx - 45.0) / 3.0) ** 2 + ((yy - 65.0) / 3.0) ** 2))
    mask = np.ones_like(reference, dtype=bool)
    dy, dx = residual_registration_shift(
        reference, moving, mask, max_shift_pix=(4.0, 4.0)
    )
    assert abs(dy) <= 4.0
    assert abs(dx) <= 4.0


def test_effective_exposure_plotter_is_available():
    # Regression test for the Script-1 QC merge that accidentally dropped this
    # pre-existing plotting helper.
    assert callable(plotting.plot_effective_exposure)
