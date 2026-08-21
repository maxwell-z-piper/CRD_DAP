from types import SimpleNamespace

import numpy as np

from crd_utils import binning


def test_geometric_coadd_scales_surface_brightness_by_pixel_area():
    flux = np.ones((2, 2, 4), dtype=float)
    variance = np.ones_like(flux)
    good = np.ones_like(flux, dtype=bool)
    cube = SimpleNamespace(
        flux=flux,
        variance=variance,
        good=good,
        good_spaxel=np.ones((2, 2), dtype=bool),
        nwave=4,
        header={"BUNIT": "1e-16 erg/s/cm2/arcsec2/Angstrom"},
    )
    bin_map = np.array([[0, 0], [1, 1]], dtype=int)
    result = binning.coadd_bin_spectra(
        cube,
        bin_map,
        n_bins=2,
        pixel_area_arcsec2=0.09,
        min_member_fraction=1.0,
    )
    # Two unit-surface-brightness pixels per bin, each covering 0.09 arcsec^2.
    assert np.allclose(result.flux, 0.18)
    assert np.allclose(result.uncertainty, np.sqrt(2.0) * 0.09)
    assert np.all(result.good)
    assert np.array_equal(result.n_members, [2, 2])


def test_coadd_masks_wavelength_when_too_few_members_are_valid():
    flux = np.ones((1, 2, 3), dtype=float)
    variance = np.ones_like(flux)
    good = np.ones_like(flux, dtype=bool)
    good[0, 1, 1] = False
    cube = SimpleNamespace(
        flux=flux,
        variance=variance,
        good=good,
        good_spaxel=np.ones((1, 2), dtype=bool),
        nwave=3,
        header={"BUNIT": "arbitrary"},
    )
    bin_map = np.array([[0, 0]], dtype=int)
    result = binning.coadd_bin_spectra(
        cube,
        bin_map,
        n_bins=1,
        pixel_area_arcsec2=1.0,
        min_member_fraction=1.0,
    )
    assert result.good[0, 0]
    assert not result.good[0, 1]
    assert result.good[0, 2]
    assert np.isnan(result.flux[0, 1])
