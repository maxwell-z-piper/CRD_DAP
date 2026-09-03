# CRD_DAP Script-3 covariance-calibration update

This bundle updates Script 3 so the RH3 profile-likelihood cubes are evaluated
with an empirically calibrated spectral covariance model rather than the old
development-only diagonal-noise approximation.

## Prerequisite

The `crd_dap` conda environment must contain the validated pPXF 9.4.8
cached-whitener patch in the repository's `ppxf_patch_9_4_8/` directory.  The
active pPXF object must expose the `noise_inv_cholesky` keyword. Script 3 checks
this and hard-fails before calibration if the patch is absent.

## Replace these repository files

```text
scripts/03_build_RH3_likelihood_cubes.py
crd_utils/ppxf_utils.py
crd_utils/ppxf_grid.py
METHODS.md
DIAGNOSTICS.md
```

## Add these new repository files

```text
crd_utils/covariance.py
crd_utils/covariance_calibration.py
crd_utils/covariance_plotting.py
tests/test_covariance.py
tests/test_covariance_calibration.py
tests/test_ppxf_utils_covariance.py
```

The existing `crd_utils/plotting.py`, `crd_utils/config.py`, and
`crd_utils/__init__.py` are deliberately **not** replaced by this bundle. The
new covariance plotting code lives in `covariance_plotting.py`, and the new
covariance config validation lives close to the machinery in
`covariance_calibration.py`, minimizing regression risk in central utilities.

## Add the covariance settings to configuration

Copy the contents of

```text
config/SCRIPT03_COVARIANCE_CONFIG_BLOCK.txt
```

into the Script-3 portion of both:

1. `config/target_config_template.py`, and
2. the target-specific config used for the run (for example `config/8143-1902.py`).

The agreed production defaults include:

- simultaneous 95% PowerBin bootstrap bands;
- 2000 bootstrap realizations;
- lag search through 20 log-wavelength pixels;
- M1--M4 covariance hierarchy;
- convergence only when both `max fractional Delta s < 0.01` and
  `max |Delta rho| < 0.01`;
- maximum five covariance iterations;
- 12 deterministic radial validation bins (6 per PA side) plus up to two added
  off-center high-dispersion candidates;
- full production `(V_A,V_B,f_A)` grids for Requirement B;
- production grid blocked when covariance convergence or model stability fails.

## Validation after installing the files

Run in the same `crd_dap` environment where the pPXF patch already passed:

```bash
cd /Users/maxpiper/CRD_DAP
conda activate crd_dap

python -m compileall -q scripts crd_utils tests

python -m pytest -q \
    ppxf_patch_9_4_8/tests/test_ppxf_cached_whitener.py \
    tests/test_covariance.py \
    tests/test_covariance_calibration.py \
    tests/test_ppxf_utils_covariance.py
```

Then run Script 3 normally, for example:

```bash
python scripts/03_build_RH3_likelihood_cubes.py \
    --config config/8143-1902.py \
    --script1-run /Users/maxpiper/CRD_DAP/runs/8143-1902_script01 \
    --script2-run /Users/maxpiper/CRD_DAP/runs/8143-1902_script02 \
    --workers 3
```

Do not use `--resume` to continue an old diagonal-noise Script-3 run. The
checkpoint schema is now version 3 and each production checkpoint is tied to a
specific selected covariance-model hash. Start a new Script-3 run for the new
statistical experiment.

## Important products before the full production grid

The calibration gate writes, among others:

```text
products/covariance_validation_bins.ecsv
products/RH3_covariance_iteration_history.ecsv
products/RH3_covariance_candidates.npz
products/RH3_covariance_calibration_fits.npz
products/RH3_covariance_model_validation_grids.npz
products/RH3_covariance_model_comparison.ecsv
metadata/RH3_covariance_model_selection.json
```

The all-bin 2601-state likelihood calculation begins only after the saved model
selection has `production_grid_allowed=true`.

## Validation performed while generating this bundle

- Python `compileall`: PASS for the full bundle.
- Low-level covariance, representative-bin/model-selection, and explicit
  covariance-chi-square tests: **8 passed** in the generation environment using
  a minimal Astropy Table compatibility stub (Astropy itself is not installed
  in the generation container).
- The true pPXF cached-whitener regression suite could not be rerun in this
  generation container because pPXF is not installed here. The user's actual
  `crd_dap` environment already ran the independent pPXF patch suite with
  **9 passed**, and benchmarked the cached path at 2.438x faster than stock
  covariance input for the 1169-pixel test.

The first real Script-3 run remains an integration test of the complete
Astropy+pPXF+XSL+real-data path. Inspect the covariance diagnostics before
allowing a long production run to proceed unattended.
