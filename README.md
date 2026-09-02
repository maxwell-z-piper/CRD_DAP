# CRD_DAP

**CRD_DAP** is a research pipeline for decomposing counter-rotating stellar disks (CRDs) observed with Keck/KCWI/KCRM. The project is being developed around paired BL and RH3 observations and is designed to preserve as much of the spectral likelihood information as possible from the high-resolution RH3 kinematic fits before propagating those constraints into two-component stellar-population fits in BL.

The repository contains the configuration, logging, validation, diagnostic-documentation architecture, the `crd_utils` helper package, and implemented science-driver stages 01--03: cube preparation/registration, master BL PowerBins, and independent RH3 profile-likelihood cubes. Later science-driver scripts described in `METHODS.md` will be implemented sequentially.

**Production CRD_DAP runs begin from a validated `CRD_DRP` reduction package, not from arbitrary cube paths.** The upstream reduction pipeline must be completed first, and its `CRD_DRP_reduction_manifest.json` is the formal handoff contract between reduction and science analysis.

## Scientific design in one sentence

Rather than reducing every spatial bin to one independently selected two-component LOSVD, the pipeline constructs per-bin RH3 profile-likelihood surfaces, uses a globally coherent XookSuut-style two-disk rotation model to select the physically consistent solution family, and propagates the RH3-supported kinematic uncertainty into separate BL stellar-population measurements.


## Required upstream reduction: CRD_DRP

`CRD_DAP` is the **science-analysis pipeline**. The production BL and red/RH3 cubes must first be processed and validated by the separate **CRD_DRP** reduction repository.

The intended handoff is:

```text
KCWI/KCRM raw data
        ↓
KCWI DRP / KSkyWizard / KcwiKit
        ↓
CRD_DRP stack-integrity validation
        ↓
CRD_DRP BLUE/RED atmospheric-mask construction
        ↓
CRD_DRP reduction-package validation
        ↓
CRD_DRP_reduction_manifest.json
        ↓
CRD_DAP Script 01
        ↓
Scripts 02--07
```

The CRD_DRP manifest records the validated BLUE and RED KcwiKit `i/v/m/e` stacks, the finalized atmospheric masks, file hashes, wavelength-grid provenance, and the reduction-package PASS/FAIL state. For production analysis, Script 01 should therefore be launched with:

```bash
python scripts/01_prepare_and_register_cubes.py \
    --config config/<target>.py \
    --reduction-manifest /path/to/CRD_DRP_reduction_manifest.json
```

Script 01 requires the CRD_DRP package and both arms to be `PASS`, maps `BLUE -> BL` and `RED -> RH3`, verifies the referenced files and SHA256 hashes, checks that each atmospheric mask matches the native science wavelength grid/medium, and then folds those wavelength exclusions into the prepared `GOODWAVE`/`GOODMASK` products.

Direct science-cube paths in the target configuration remain available for legacy/testing workflows, but they are **not the preferred production handoff**. Master arcs, PyMorph/XSL inputs, target redshift, and analysis settings remain CRD_DAP configuration inputs.

The optional CaT mask audit lives upstream in CRD_DRP. It is a reduction-QC tool for checking whether enough of the expected Ca II triplet survives the finalized RED atmospheric mask; it is not a required CRD_DAP stage.

## Repository layout

```text
CRD_DAP/
├── README.md
├── METHODS.md
├── DIAGNOSTICS.md
├── pyproject.toml
├── .gitignore
├── config/
│   ├── __init__.py
│   └── target_config_template.py
├── crd_utils/
│   ├── __init__.py
│   ├── config.py
│   ├── logging_utils.py
│   ├── io.py
│   ├── cube_utils.py
│   ├── noise.py
│   ├── psf_lsf.py
│   ├── binning.py
│   ├── templates.py
│   ├── ppxf_utils.py
│   ├── ppxf_grid.py
│   ├── likelihood.py
│   ├── model_selection.py
│   ├── geometry.py
│   ├── disk_model.py
│   ├── populations.py
│   ├── sampling.py
│   ├── validation.py
│   ├── quality.py
│   ├── runtime.py
│   └── plotting.py
├── scripts/
│   ├── README.md
│   ├── 01_prepare_and_register_cubes.py
│   ├── 02_make_master_BL_bins.py
│   └── 03_build_RH3_likelihood_cubes.py
├── archive/
│   └── script03_forensics/
│       └── 03a--03g historical development/forensic utilities
├── tests/
├── examples/
├── data/
│   ├── README.md
│   ├── 01_prepare_and_register_cubes.py
│   ├── 02_make_master_BL_bins.py
│   └── 03_build_RH3_likelihood_cubes.py
└── runs/
```


## Current implementation status

- **Script 01:** implemented with the production `CRD_DRP_reduction_manifest.json` handoff. It verifies the validated reduction package, resolves BLUE/RED science products, applies the finalized CRD_DRP atmospheric masks, and preserves both native and combined mask provenance in the prepared cubes.
- **Script 02:** implemented and integration-tested; BL defines one master PowerBin tessellation and the same physical memberships are transferred to the red/RH3 stream. CRD_DRP wavelength exclusions are inherited automatically through Script-1 `GOODWAVE`/`GOODMASK`.
- **Script 03:** implemented for exact independent `(V_A,V_B,f_A,RH3)` profile-likelihood cubes, one-component controls, empirical-LSF XSL preparation, restart checkpoints, and multiprocessing. The production atmospheric mask is inherited from Scripts 01--02 rather than rebuilt inside Script 03. The first development pass uses formal diagonal spectral uncertainties and explicitly marks likelihood widths as uncalibrated until residual variance/covariance is revisited.
  Script 3 can take both `--script1-run` and `--script2-run` explicitly, which is useful when completed upstream run directories have been moved or renamed and older manifest paths are stale.
- **Historical 03a--03g utilities:** no longer part of the production execution path. They should be retained under `archive/script03_forensics/` (or Git history) as provenance for the reduction/masking investigation.
- **Scripts 04--07:** design documented; implementation follows after Script-3 integration testing.

## Intended science scripts

The executable analysis will be organized into seven readable driver scripts:

1. `01_prepare_and_register_cubes.py`
2. `02_make_master_BL_bins.py`
3. `03_build_RH3_likelihood_cubes.py`
4. `04_fit_global_two_disk_model.py`
5. `05_extract_global_RH3_solution.py`
6. `06_fit_BL_two_component_populations.py`
7. `07_uncertainties_and_final_maps.py`

The drivers should remain short and readable. Most scientific logic belongs in `crd_utils`, where functions are heavily documented and testable in isolation.

## Required external data for a target

A production target requires, at minimum:

- a **validated CRD_DRP reduction manifest** (`CRD_DRP_reduction_manifest.json`), which points to the final BLUE/RED `i/v/m/e` stacks and their atmospheric masks;
- the corresponding KCWI/KCRM master-arc calibration product(s);
- the MaNGA PyMorph VAC or the relevant local extract from it;
- the XSL SSP template library used for pPXF fitting;
- initial target information such as redshift and a stellar kinematic PA estimate.

The science-stack filenames and atmospheric-mask filenames are resolved from the CRD_DRP manifest rather than independently re-entered for a normal production run. The master arcs remain CRD_DAP inputs because the adopted plan measures the instrumental LSF empirically from calibration data rather than relying only on a nominal resolving power.

## Configuration

Copy `config/target_config_template.py` to a target-specific file and edit the calibration, ancillary-data, target, and analysis settings there. The helper `crd_utils.config.load_config()` validates required entries, and each run saves an immutable snapshot of the configuration alongside the outputs.

For production science, the reduced BL/RH3 stack paths are supplied by `--reduction-manifest`; they should not be treated as an independent second source of truth in the target config. The old direct KcwiKit/DRP science-input fields remain only for legacy/testing mode.

After Script 01, Scripts 02 and 03 continue to use explicit upstream run directories in the normal way:

```bash
python scripts/02_make_master_BL_bins.py     --config config/<target>.py     --script1-run /path/to/script01_run

python scripts/03_build_RH3_likelihood_cubes.py     --config config/<target>.py     --script1-run /path/to/script01_run     --script2-run /path/to/script02_run     --workers <N>
```

`RH3_ATMOSPHERIC_MASK_FILE` is a retired 03d-era production pathway and should remain `None`. Any deliberate *additional* science masks belong in `RH3_MASK_OBSERVED_RANGES_ANGSTROM` / `RH3_MASK_REST_RANGES_ANGSTROM`; the routine sky/telluric mask already comes from CRD_DRP.

## Logging and reproducibility

Each run is intended to create a timestamped directory under `runs/` containing:

- a snapshot of the target configuration;
- `pipeline.log`;
- a step-specific log for every science script;
- numerical products;
- diagnostic figures;
- metadata needed to reproduce the run, including the source CRD_DRP manifest and propagated reduction/mask provenance when production inputs are used.

The logging utilities are designed to write to both the terminal and log files. Long pPXF calculations can therefore be left unattended without losing the progress and timing history. During Script 3's expensive per-bin likelihood stage, an interactive terminal additionally shows one carriage-return heartbeat/progress line (spinner, completed PowerBins, elapsed time, time since the last completed bin, and ETA after the first completion). That dynamic line is deliberately **not** copied into the permanent log. It indicates that the parent process remains alive and has not yet received a worker exception; it is not a substitute for the fit/QC diagnostics.

## Documentation

- `METHODS.md` is the detailed scientific and statistical design document.
- `DIAGNOSTICS.md` is the detailed catalogue of every planned diagnostic figure, what it shows, how it is calculated, how to interpret it, and what action to take if it looks wrong.

Both files should evolve together with the code. A new science diagnostic should not be considered complete until it is documented in `DIAGNOSTICS.md`.

## Development philosophy

The project deliberately favors scientific completeness, traceability, and explicit diagnostics over minimizing runtime or disk usage. Expensive calculations may be cached and parallelized, but approximations should be introduced only when their effect has been tested.

The numerical settings in the template configuration are **development defaults**, not universal scientific truths. Values such as grid resolution, polynomial degree, Monte Carlo count, and geometry-prior widths are expected to be validated on mocks and real data before publication.


Script 3 now derives template wavelength padding from both the exact two-component grid and the one-component control velocity domain, then runs a real pPXF wavelength-coverage preflight before launching worker processes.
