from __future__ import annotations

import numpy as np

from crd_utils import covariance
from crd_utils.ppxf_utils import total_chi2


def test_total_chi2_uses_cached_whitener_exactly():
    n = 30
    noise = np.linspace(1.0, 1.5, n)
    good = np.ones(n, dtype=bool)
    rho = np.zeros((1, 4)); rho[0, 0] = 1; rho[0, 1] = 0.2
    W = covariance.build_inverse_cholesky(noise, good, scale=1.1, rho_by_block=rho).inv_cholesky
    galaxy = np.linspace(0, 1, n)
    bestfit = galaxy + 0.03*np.sin(np.arange(n))
    r = galaxy - bestfit
    expected = np.sum((W @ r)**2)
    got = total_chi2(galaxy, bestfit, noise, np.arange(n), noise_inv_cholesky=W)
    np.testing.assert_allclose(got, expected, rtol=1e-12, atol=1e-12)
