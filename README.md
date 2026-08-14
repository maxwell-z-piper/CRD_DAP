# CRD_DAP

**CRD_DAP** is a research pipeline for decomposing counter-rotating stellar disks (CRDs) observed with Keck/KCWI/KCRM. The project is being developed around paired BL and RH3 observations and is designed to preserve as much of the spectral likelihood information as possible from the high-resolution RH3 kinematic fits before propagating those constraints into two-component stellar-population fits in BL.

The current repository contains the **pre-Script-1 infrastructure**: configuration, logging, validation, diagnostic-documentation architecture, and the `crd_utils` helper package. The seven science-driver scripts described in `METHODS.md` will be implemented next.

## Scientific design in one sentence

Rather than reducing every spatial bin to one independently selected two-component LOSVD, the pipeline constructs per-bin RH3 profile-likelihood surfaces, uses a globally coherent XookSuut-style two-disk rotation model to select the physically consistent solution family, and propagates the RH3-supported kinematic uncertainty into separate BL stellar-population measurements.

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
│   └── README.md
├── tests/
├── examples/
├── data/
│   └── README.md
└── runs/
```

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

The configuration architecture assumes the analysis will eventually require, at minimum:

- a fully reduced and stacked BL cube;
- a fully reduced and stacked RH3 cube;
- the corresponding KCWI DRP master arc calibration product(s);
- the MaNGA PyMorph VAC or the relevant local extract from it;
- the XSL SSP template library used for pPXF fitting;
- initial target information such as redshift and a stellar kinematic PA estimate.

The master arcs are treated as required inputs because the adopted plan measures the instrumental LSF empirically from the calibration data rather than relying only on a nominal resolving power.

## Configuration

Copy `config/target_config_template.py` to a target-specific file and edit the paths and analysis settings there. The helper `crd_utils.config.load_config()` validates required entries, and each run should save an immutable snapshot of the configuration alongside the outputs.

## Logging and reproducibility

Each run is intended to create a timestamped directory under `runs/` containing:

- a snapshot of the target configuration;
- `pipeline.log`;
- a step-specific log for every science script;
- numerical products;
- diagnostic figures;
- metadata needed to reproduce the run.

The logging utilities are designed to write to both the terminal and log files. Long pPXF calculations can therefore be left unattended without losing the progress and timing history.

## Documentation

- `METHODS.md` is the detailed scientific and statistical design document.
- `DIAGNOSTICS.md` is the detailed catalogue of every planned diagnostic figure, what it shows, how it is calculated, how to interpret it, and what action to take if it looks wrong.

Both files should evolve together with the code. A new science diagnostic should not be considered complete until it is documented in `DIAGNOSTICS.md`.

### Markdown mathematics

Repository documentation follows GitHub's native MathJax Markdown convention. Use `$...$` for inline mathematics and standalone `$$...$$` blocks for displayed equations. Do not use `\(...\)` or `\[...\]` delimiters in repository Markdown files. This keeps mathematical notation rendered consistently in GitHub's **Preview** view.

## Development philosophy

The project deliberately favors scientific completeness, traceability, and explicit diagnostics over minimizing runtime or disk usage. Expensive calculations may be cached and parallelized, but approximations should be introduced only when their effect has been tested.

The numerical settings in the template configuration are **development defaults**, not universal scientific truths. Values such as grid resolution, polynomial degree, Monte Carlo count, and geometry-prior widths are expected to be validated on mocks and real data before publication.
