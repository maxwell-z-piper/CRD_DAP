# `crd_utils` package map

`crd_utils` is a Python **package** (a directory of modules), not one monolithic class. The seven science-driver scripts should remain readable and call reusable, heavily documented functions from this package.

Current module responsibilities:

- `config.py` — target-configuration loading, validation, snapshots, machine-readable manifests.
- `logging_utils.py` — timestamped run directories, terminal + master log + step logs.
- `io.py` — FITS/table/metadata I/O abstraction.
- `cube_utils.py` — collapsed continuum, masks, cube registration helpers.
- `noise.py` — preliminary Script-1 noise and spectral-correlation diagnostics.
- `spatial_noise.py` — empirical per-bin scaling from fixed-stack Gaussian noise ensembles.
- `scalar_noise.py` — per-bin residual scaling for Script-3 diagonal pPXF fits.
- `psf_lsf.py` — empirical master-arc LSF and delivered PSF characterization.
- `binning.py` — BL-defined PowerBin membership and transfer to RH3.
- `templates.py` — XSL SSP wavelength/LSF/normalization preparation.
- `ppxf_utils.py` — standardized low-level pPXF wrappers.
- `ppxf_grid.py` — explicit RH3 `(V_A, V_B, f_A)` grid construction and evaluation.
- `likelihood.py` — profile-likelihood transforms, basin topology, likelihood-concentration diagnostics.
- `model_selection.py` — one/two-component statistic and mock-calibrated false-positive/recovery metrics.
- `geometry.py` — PyMorph-derived inclination, coordinate transforms, signed PA convention.
- `disk_model.py` — XookSuut-style ring grid, interpolation, projected velocities, bin-integrated shear.
- `global_disk.py` — Script-04 fixed-sense two-disk projection, per-fraction likelihood interpolation, curvature and global search.
- `global_disk_io.py` — portable upstream contracts, whole-bin domains, reviewed geometry inputs and atomic products.
- `stellar_psf.py` — native-pixel-integrated Gaussian/Moffat standard fits, effective PSF transfer and stellar QA.
- `global_disk_plotting.py` — the 2×6 strength-review sheet, geometry and eight-model comparison figures.
- `populations.py` — BL two-component stellar-population fits and SSP summaries.
- `sampling.py` — RH3-likelihood-weighted sampling, ESS, convergence, eventual joint-minimum refinement.
- `validation.py` — hard scientific checks and grid/boundary validation.
- `quality.py` — centralized quality-flag definitions.
- `runtime.py` — timings, wall-time forecasting, later multicore benchmark support.
- `plotting.py` — all diagnostic/publication plot implementations; every nontrivial plot must be documented in `../DIAGNOSTICS.md`.

Script-1 functionality is implemented in `io.py`, `cube_utils.py`, `noise.py`, `psf_lsf.py`, `validation.py`, and `plotting.py`. Script 04 uses the `global_disk*` modules; the older `disk_model.py` primitives retain their original contracts for existing callers. Other later-stage modules may still contain explicit `NotImplementedError` placeholders, filled only when their driver is implemented and tested.
