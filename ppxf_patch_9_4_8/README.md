# CRD_DAP pPXF 9.4.8 cached-whitener patch

## Purpose

CRD_DAP will reuse one fixed spectral covariance model for many pPXF calls within a PowerBin. Stock pPXF accepts a full covariance matrix, but it converts that matrix into its internal whitening operator during every pPXF construction. This patch adds an optional keyword:

```python
noise_inv_cholesky=W
```

where

$$
C = L L^T
$$
$$
W = L^{-1}
$$

$C$ is the spectral covariance matrix, $L$ is its lower-triangular Cholesky factor, and $W$ is the inverse lower-Cholesky factor. Supplying $W$ allows CRD_DAP to factor a frozen covariance matrix once and reuse the exact same whitening operator for all subsequent pPXF fits for that spectrum.

The normal positional `noise` argument is still required by pPXF. In cached mode CRD_DAP should pass a valid positive 1-D uncertainty vector there; after standard pPXF input validation the patched code replaces pPXF's internal weighting operator with the supplied $W$.

## Scope

The patch is deliberately locked to **pPXF 9.4.8**.

Only one installed pPXF source file is modified:

```text
<active environment>/site-packages/ppxf/ppxf.py
```

No `capfit.py`, utility modules, examples, or other pPXF files are changed.

A backup is created beside the original file before any modification.

## Why this package does not include a replacement `ppxf.py`

The current pPXF license permits non-commercial use and personal/internal modification but states that the code may not be redistributed. For that reason this CRD_DAP patch package does **not** contain a copied or modified pPXF source file. Instead, `apply_ppxf_9_4_8_cached_whitener_patch.py` modifies the user's own installed pPXF 9.4.8 file in place after validating the active environment and source structure.

This also makes the environment boundary explicit: the patcher prints the exact `ppxf.py` path before editing it and, by default, refuses to modify a pPXF installation outside the active Python `sys.prefix`.

## Conda environment behavior

A normal conda environment has its own Python executable and its own environment-specific `site-packages`. Therefore, if `crd_dap` is activated and pPXF was installed normally into that environment, patching its `ppxf.py` affects that environment only.

There are exceptions: `pip install --user`, a shared/editable install, or `PYTHONPATH` can cause Python to import a package from outside the active conda environment. The patcher checks for this and refuses by default if the imported `ppxf.py` is outside `sys.prefix`.

## Installation

Activate the CRD_DAP environment first:

```bash
conda activate crd_dap
```

Confirm pPXF 9.4.8 is the version imported by that environment:

```bash
python -c "import ppxf, sys; print(sys.executable); print(ppxf.__version__)"
```

Dry-run the patcher:

```bash
python ppxf_patch_9_4_8/apply_ppxf_9_4_8_cached_whitener_patch.py --check
```

The output reports:

- Python executable;
- active `sys.prefix`;
- pPXF version;
- exact installed `ppxf.py` path;
- whether that path is inside the active environment;
- whether the patch is already installed;
- whether the pPXF 9.4.8 source structure is patchable.

Apply the patch:

```bash
python ppxf_patch_9_4_8/apply_ppxf_9_4_8_cached_whitener_patch.py --apply
```

The patcher:

1. requires pPXF exactly `9.4.8`;
2. requires the imported `ppxf.py` to live under the active `sys.prefix` unless explicitly overridden;
3. finds the pPXF covariance preprocessing block structurally with Python's AST rather than depending on fragile line numbers;
4. adds the `noise_inv_cholesky` keyword;
5. installs the supplied inverse-Cholesky matrix directly as pPXF's internal whitening operator after the ordinary `noise` input has passed standard validation;
6. compiles the modified source before writing it;
7. saves the original `ppxf.py` as a backup;
8. writes the change atomically;
9. launches a fresh Python subprocess to verify that the patched keyword imports correctly;
10. automatically restores the backup if that post-write smoke test fails.

## Regression tests

After applying the patch:

```bash
python -m pytest -q ppxf_patch_9_4_8/tests/test_ppxf_cached_whitener.py
```

The tests require direct numerical equivalence between stock covariance pPXF and cached-whitener pPXF for both:

- a one-component nonlinear fit;
- a fixed-velocity two-component fit resembling one CRD_DAP Script-3 grid state.

The tests compare the pPXF solution, best-fitting spectrum, template weights, and reduced chi-square to tight numerical tolerances. They also verify that the cached path does not call SciPy's Cholesky routine inside pPXF and that malformed cached matrices are rejected.

## Benchmark

After the tests pass:

```bash
python ppxf_patch_9_4_8/tests/benchmark_ppxf_cached_whitener.py --npix 1169 --repeats 10
```

This compares repeated stock-covariance fits against repeated cached-whitener fits on the same environment. The included synthetic problem has only a small template basis, so do not use the absolute runtime to forecast the final 416-template CRD_DAP run. The final Script-3 implementation should additionally benchmark a small number of real XSL fits.

## Restore the original pPXF installation

```bash
python ppxf_patch_9_4_8/apply_ppxf_9_4_8_cached_whitener_patch.py --restore
```

The patcher restores the backup made at installation time.

Alternatively, reinstalling pPXF 9.4.8 inside the environment will replace the local modification.

## CRD_DAP integration contract

For PowerBin `i`, once covariance calibration has converged:

$$
C_i = L_i L_i^T
$$
$$
W_i = L_i^{-1}
$$

CRD_DAP should save and freeze $W_i$. Every likelihood-grid state for that PowerBin then uses the same $W_i$ through `noise_inv_cholesky`.

The covariance/whitener must be recalculated if the actual spectral experiment changes, including bin membership, wavelength grid, masks, reduction, or other operations that change the noise covariance. Merely refining $(V_A, V_B, f_A)$ later does not require recalibration; the same frozen RH3 whitener can be reused by later RH3 refinement stages.
