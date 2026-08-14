# `crd_utils` package map

`crd_utils` is a Python **package** (a directory of modules), not one monolithic class. The seven science-driver scripts should remain readable and call reusable, heavily documented functions from this package.

Current module responsibilities:

- `config.py` — target-configuration loading, validation, snapshots, machine-readable manifests.
- `logging_utils.py` — timestamped run directories, terminal + master log + step logs.
- `io.py` — FITS/table/metadata I/O abstraction.
- `cube_utils.py` — collapsed continuum, masks, cube registration helpers.
- `noise.py` — variance scaling, spectral covariance, correlated noise realizations.
- `psf_lsf.py` — empirical master-arc LSF and delivered PSF characterization.
- `binning.py` — BL-defined PowerBin membership and transfer to RH3.
- `templates.py` — XSL SSP wavelength/LSF/normalization preparation.
- `ppxf_utils.py` — standardized low-level pPXF wrappers.
- `ppxf_grid.py` — explicit RH3 `(V_A, V_B, f_A)` grid construction and evaluation.
- `likelihood.py` — profile-likelihood transforms, basin topology, likelihood-concentration diagnostics.
- `model_selection.py` — one/two-component statistic and mock-calibrated false-positive/recovery metrics.
- `geometry.py` — PyMorph-derived inclination, coordinate transforms, signed PA convention.
- `disk_model.py` — XookSuut-style ring grid, interpolation, projected velocities, bin-integrated shear.
- `populations.py` — BL two-component stellar-population fits and SSP summaries.
- `sampling.py` — RH3-likelihood-weighted sampling, ESS, convergence, eventual joint-minimum refinement.
- `validation.py` — hard scientific checks and grid/boundary validation.
- `quality.py` — centralized quality-flag definitions.
- `runtime.py` — timings, wall-time forecasting, later multicore benchmark support.
- `plotting.py` — all diagnostic/publication plot implementations; every nontrivial plot must be documented in `../DIAGNOSTICS.md`.

Several modules currently contain deliberately explicit `NotImplementedError` placeholders. This is intentional: the interfaces and scientific responsibilities are defined now, while implementation that depends on the exact real KCWI DRP file structure will be added alongside the corresponding science script.
