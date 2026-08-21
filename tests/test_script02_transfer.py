import numpy as np
from astropy.wcs import WCS

from crd_utils.binning import transfer_bin_map_by_wcs


def _wcs():
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [150.0, 2.0]
    w.wcs.crpix = [3.0, 3.0]
    w.wcs.cdelt = np.array([-0.3 / 3600.0, 0.3 / 3600.0])
    return w


def test_identity_wcs_transfers_membership_exactly():
    bin_map = np.array([
        [-1, -1, -1, -1, -1],
        [-1,  0,  0,  1, -1],
        [-1,  0,  0,  1, -1],
        [-1,  2,  2,  2, -1],
        [-1, -1, -1, -1, -1],
    ])
    good = np.ones_like(bin_map, dtype=bool)
    result = transfer_bin_map_by_wcs(
        bin_map,
        _wcs(),
        _wcs(),
        good,
        bl_pixel_scale_arcsec=0.3,
        max_distance_arcsec=0.05,
    )
    assert np.array_equal(result.bin_map[bin_map >= 0], bin_map[bin_map >= 0])
    assert np.isclose(result.assigned_fraction, 1.0)
