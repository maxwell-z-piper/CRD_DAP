import numpy as np

from crd_utils.cube_utils import registration_contrast_snr


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
