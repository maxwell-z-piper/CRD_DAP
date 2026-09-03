#!/usr/bin/env python3
"""Small runtime benchmark: stock covariance pPXF versus cached inverse-Cholesky.

This is intentionally a benchmark, not a pass/fail unit test.  Run after the
patch is installed, e.g.::

    python ppxf_patch_9_4_8/tests/benchmark_ppxf_cached_whitener.py --npix 1169 --repeats 10

The synthetic template problem is much smaller than CRD_DAP's 416-template XSL
basis, so the absolute times are not a predictor of the final run.  The useful
quantity is the repeated covariance-setup overhead that disappears in cached
mode on the same machine/environment.
"""
from __future__ import annotations

import argparse
import inspect
import time

import numpy as np
from scipy import linalg
import ppxf
from ppxf.ppxf import ppxf as PPXF


def make_problem(npix: int, rho: float = 0.3):
    x = np.linspace(-1.0, 1.0, npix)
    t1 = 1.0 - 0.25*np.exp(-0.5*((x + 0.2)/0.05)**2) - 0.12*np.exp(-0.5*((x - 0.37)/0.07)**2)
    t2 = 1.0 - 0.16*np.exp(-0.5*((x + 0.05)/0.08)**2) - 0.10*np.exp(-0.5*((x - 0.48)/0.05)**2)
    templates = np.column_stack([t1, t2])
    galaxy = 0.6*t1 + 0.4*t2
    sigma = np.full(npix, 0.02)
    lag = np.abs(np.subtract.outer(np.arange(npix), np.arange(npix)))
    cov = np.outer(sigma, sigma)*(rho**lag)
    L = linalg.cholesky(cov, lower=True)
    W = linalg.solve_triangular(L, np.eye(npix), lower=True)
    good = np.arange(4, npix - 4, dtype=int)
    return templates, galaxy, sigma, cov, W, good


def fit(templates, galaxy, noise, good, **kwargs):
    return PPXF(
        templates, galaxy, noise, 55.0, [0.0, 165.0],
        moments=2, degree=2, mdegree=0, goodpixels=good,
        clean=False, quiet=True, **kwargs
    )


def timed(repeats, fn):
    samples = []
    for _ in range(repeats):
        t0 = time.perf_counter()
        fn()
        samples.append(time.perf_counter() - t0)
    return np.asarray(samples)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--npix", type=int, default=512)
    ap.add_argument("--repeats", type=int, default=5)
    args = ap.parse_args()

    if str(ppxf.__version__) != "9.4.8":
        raise RuntimeError(f"Expected pPXF 9.4.8; found {ppxf.__version__}")
    if "noise_inv_cholesky" not in inspect.signature(PPXF).parameters:
        raise RuntimeError("Cached-whitener patch is not installed.")

    templates, galaxy, errvec, cov, W, good = make_problem(args.npix)

    # Warm-up imports/FFT machinery.
    fit(templates, galaxy, errvec, good, noise_inv_cholesky=W)

    stock = timed(args.repeats, lambda: fit(templates, galaxy, cov, good))
    cached = timed(args.repeats, lambda: fit(templates, galaxy, errvec, good, noise_inv_cholesky=W))

    print(f"pPXF version: {ppxf.__version__}")
    print(f"npix={args.npix}, repeats={args.repeats}")
    print(f"stock covariance median : {np.median(stock):.4f} s/fit")
    print(f"cached whitener median  : {np.median(cached):.4f} s/fit")
    print(f"ratio stock/cached      : {np.median(stock)/np.median(cached):.3f}x")
    print("NOTE: benchmark the real CRD_DAP 416-template fit separately before forecasting total runtime.")


if __name__ == "__main__":
    main()
