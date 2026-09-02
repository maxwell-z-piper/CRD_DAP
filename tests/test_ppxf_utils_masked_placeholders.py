"""Regression tests for pPXF-safe handling of already-masked samples."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import sys

import numpy as np
import pytest


MODULE_PATH = Path(__file__).resolve().parents[1] / "crd_utils" / "ppxf_utils.py"
spec = importlib.util.spec_from_file_location("crd_dap_patch_ppxf_utils", MODULE_PATH)
ppxf_utils = importlib.util.module_from_spec(spec)
assert spec.loader is not None
sys.modules[spec.name] = ppxf_utils
spec.loader.exec_module(ppxf_utils)


def _base_vectors():
    galaxy = np.array([1.0, np.nan, 3.0, 4.0, 5.0, np.inf, 7.0, 8.0])
    noise = np.array([1.0, np.nan, 2.0, 2.5, 3.0, 0.0, 4.0, 5.0])
    lam = np.arange(8.0) + 7000.0
    goodpixels = np.array([0, 2, 3, 4, 6, 7])
    return galaxy, noise, lam, goodpixels


def test_excluded_invalid_samples_are_sanitized_without_mutating_inputs():
    galaxy, noise, lam, goodpixels = _base_vectors()
    galaxy_original = galaxy.copy()
    noise_original = noise.copy()

    g2, n2, l2, gp2 = ppxf_utils._validate_spectrum_inputs(
        galaxy, noise, lam, goodpixels
    )

    # pPXF's full-vector contract is satisfied.
    assert np.all(np.isfinite(g2))
    assert np.all(np.isfinite(n2))
    assert np.all(n2 > 0)

    # Fitted samples are scientifically untouched.
    np.testing.assert_array_equal(g2[goodpixels], galaxy[goodpixels])
    np.testing.assert_array_equal(n2[goodpixels], noise[goodpixels])
    np.testing.assert_array_equal(l2, lam)
    np.testing.assert_array_equal(gp2, goodpixels)

    # Caller-owned science arrays preserve their NaN/Inf/zero mask provenance.
    np.testing.assert_equal(galaxy, galaxy_original)
    np.testing.assert_equal(noise, noise_original)


def test_nonpositive_noise_inside_goodpixels_is_rejected():
    galaxy, noise, lam, goodpixels = _base_vectors()
    noise[3] = 0.0
    with pytest.raises(ValueError, match="finite and positive on every good pixel"):
        ppxf_utils._validate_spectrum_inputs(galaxy, noise, lam, goodpixels)


def test_nonfinite_noise_inside_goodpixels_is_rejected():
    galaxy, noise, lam, goodpixels = _base_vectors()
    noise[3] = np.nan
    with pytest.raises(ValueError, match="finite and positive on every good pixel"):
        ppxf_utils._validate_spectrum_inputs(galaxy, noise, lam, goodpixels)


def test_nonfinite_galaxy_inside_goodpixels_is_rejected():
    galaxy, noise, lam, goodpixels = _base_vectors()
    galaxy[3] = np.nan
    with pytest.raises(ValueError, match="galaxy must be finite on every good pixel"):
        ppxf_utils._validate_spectrum_inputs(galaxy, noise, lam, goodpixels)
