# CRD_DAP pipeline calls

Quick reference for **DAP Scripts 01–05**. Mac / Bash calls are the default;
Windows / PowerShell equivalents are in their own section.

Scripts 01–03 follow the DAP handoff calls at the bottom of
[INSTRUCTIONS_COMPLETE_NOISE_REDUCTION.md](/Users/maxpiper/CRD_DRP/INSTRUCTIONS_COMPLETE_NOISE_REDUCTION.md).
Script 04 requires a reviewed PSF, smoothing strength and model selection before Script 05 can run.
Update paths/run names for another target or reduction. Use fresh run names
when starting a new run.

## Mac: paths

Set once per terminal. The reduction manifest below is the existing Mac
`_red_zap` handoff; substitute a newly generated manifest when appropriate.
The target config must point to the matching science/noise products.

```bash
export DAP=/Users/maxpiper/CRD_DAP
export DAP_PY=/Users/maxpiper/miniforge3/envs/crd_dap/bin/python
export CONFIG="$DAP/config/8143-1902.py"
export REDUCTION_MANIFEST=/Users/maxpiper/CRD_DRP/validation_products/03_reduction_package_red_zap/CRD_DRP_reduction_manifest.json

export S01_NAME=8143-1902_scalar_redzap_S01
export S02_NAME=8143-1902_scalar_redzap_S02
export S03_NAME=8143-1902_scalar_redzap_S03
export SCRIPT1_RUN="$DAP/runs/$S01_NAME"
export SCRIPT2_RUN="$DAP/runs/$S02_NAME"
export SCRIPT3_RUN="$DAP/runs/$S03_NAME"
export SCRIPT04_DIR="$DAP/runs/8143-1902_global04"
export SCRIPT05_DIR="$DAP/runs/8143-1902_exact05"

cd "$DAP"
```

### 01 — Prepare and register cubes

Use the completed DRP reduction-package manifest:

```bash
"$DAP_PY" "$DAP/scripts/01_prepare_and_register_cubes.py" \
  --config "$CONFIG" --reduction-manifest "$REDUCTION_MANIFEST" --run-name "$S01_NAME"
```

### 02 — Master BL bins and spatial noise scalars

```bash
"$DAP_PY" "$DAP/scripts/02_make_master_BL_bins.py" \
  --config "$CONFIG" --script1-run "$SCRIPT1_RUN" --run-name "$S02_NAME"
```

### 03 — RH3 likelihood cubes and residual noise scalars

```bash
"$DAP_PY" "$DAP/scripts/03_build_RH3_likelihood_cubes.py" \
  --config "$CONFIG" --script1-run "$SCRIPT1_RUN" --script2-run "$SCRIPT2_RUN" --run-name "$S03_NAME" \
  --workers 3
```

To resume the same interrupted Script-03 run:

```bash
"$DAP_PY" "$DAP/scripts/03_build_RH3_likelihood_cubes.py" \
  --config "$CONFIG" --script1-run "$SCRIPT1_RUN" --script2-run "$SCRIPT2_RUN" --run-name "$S03_NAME" \
  --workers 3 --resume
```

### 04 — Coherent two-disk models

Fill the MaNGA inputs, select the PyMorph HDU (or explicit inclination inputs),
and check the standard-star paths in the config. With manual inclination fields
set to `None`, the code calculates inclination/bounds from the selected table.
Adding stage-04 inputs alone does not require rerunning Scripts 01–03.

**A. Measure both PSF candidates.**

```bash
"$DAP_PY" "$DAP/scripts/04_fit_global_two_disk_model.py" \
  --config "$CONFIG" --output-dir "$SCRIPT04_DIR" --phase psf
```

Open `$SCRIPT04_DIR/figures/stellar_PSF_QA.pdf`. This 4x4 panel plots the normalized stellar light of the actual exposure, the data - a `gaussian` fit, the data - a `moffat` fit, and then those residuals against their radius. We are most interested in the residual plots and look for four things: 
One   - Is there unrecognizable structure in the residuals? If so, this indicates the model has the wrong shape, width, or centering. 
Two   - Does the model reduce the residuals both in the center and the surrounding light?
Three - Are the results consistent across exposures?
Four  - Are the fitted parameters being sensible and converging?
Which ever model answers 'yes' to the above questions is the one you should choose. Note, it has been shown that the `moffat` model slightly out-performs the `gaussian` model for KCWI. Below is an example of choosing `moffat`.

```bash
export PSF_PROFILE=moffat
```

**B. Run the twelve-fit smoothing calibration.**

```bash
"$DAP_PY" "$DAP/scripts/04_fit_global_two_disk_model.py" \
  --config "$CONFIG" --output-dir "$SCRIPT04_DIR" --script1-run "$SCRIPT1_RUN" \
  --script2-run "$SCRIPT2_RUN" --script3-run "$SCRIPT3_RUN" --psf-profile "$PSF_PROFILE" \
  --phase calibrate --workers 4
```

Note - you can resume this fitting with the `--resume` arguemnt. If you pause this stage and would like to start it anew, you must delete all entries in "$SCRIPT04_DIR/runs/fits/", then rerun the above command.

**C. Inspect `figures/smoothing_calibration_2x6.pdf`, then review/select strength.**

```bash
"$DAP_PY" "$DAP/scripts/04_fit_global_two_disk_model.py" \
  --config "$CONFIG" --output-dir "$SCRIPT04_DIR" --script1-run "$SCRIPT1_RUN" \
  --script2-run "$SCRIPT2_RUN" --script3-run "$SCRIPT3_RUN" --psf-profile "$PSF_PROFILE" --phase review --interactive --workers 4 --resume
```

This step will give you three options to choose a smoothing model coefficient:
1 - this will redo the 12 fits with a range of `weaker` coefficients
2 - this will redo the 12 fits with a range of `stronger` coefficients
3 - this will allow you to choose which of the tested coefficients to use for this stage of the pipeline.

Note that you can select `0` to use the control fit that does not have any smoothing coefficient, if desired.

**D. After selecting strength, run the eight-model comparison.**

```bash
"$DAP_PY" "$DAP/scripts/04_fit_global_two_disk_model.py" \
  --config "$CONFIG" --output-dir "$SCRIPT04_DIR" --script1-run "$SCRIPT1_RUN" \
  --script2-run "$SCRIPT2_RUN" --script3-run "$SCRIPT3_RUN" --psf-profile "$PSF_PROFILE" \
  --phase compare --workers 4 --resume
```

This step runs the eight different models given the smoothing coefficient chosen earlier. Those eight models are: smoothing on OR smoothing off: fit on full aperture / restricted up to 2-sigma region, PSF off/on. This creates `$SCRIPT04_DIR/figures/radial_extent_rotationcurve_comparison.pdf`, where you can compare all eight models to choose for the next step.

**E. Inspect the comparisons and choose a model.** The name below is an example.

```bash
"$DAP_PY" "$DAP/scripts/04_fit_global_two_disk_model.py" \
  --config "$CONFIG" --output-dir "$SCRIPT04_DIR" --script1-run "$SCRIPT1_RUN" \
  --script2-run "$SCRIPT2_RUN" --script3-run "$SCRIPT3_RUN" --psf-profile "$PSF_PROFILE" \
  --phase choose-model --model-name full_psf1_smooth1 \
  --note "Replace with your reason for choosing this model"
```

The choices are `{two_sigma,full}_psf{0,1}_smooth{0,1}`. Selection writes
`products/selected_global_model.json`. The PSF, strength and model examples here
are not automatic scientific choices.

### 05 — Exact RH3 refits and kinematic maps

Use the selected model saved by Script 04. Keep the same terminal available if you want to review failed fractions as they occur.

```bash
"$DAP_PY" "$DAP/scripts/05_extract_global_RH3_solution.py" \
  --config "$CONFIG" --script4-run "$SCRIPT04_DIR" --output-dir "$SCRIPT05_DIR" \
  --script1-run "$SCRIPT1_RUN" --script2-run "$SCRIPT2_RUN" --script3-run "$SCRIPT3_RUN" \
  --workers 4
```

Append `--resume` to this call after interruption. Use `--non-interactive` for unattended initial scans. When a run exits with status 2, remove `--non-interactive` and resume in a terminal to review the saved decisions.

## Windows: PowerShell

These use `C:/Users/maxhi/CRD_DAP` and the run names recorded in your Windows
metadata. Use Conda-initialized PowerShell. `conda run` avoids assuming a Windows
Miniforge installation path. Update the Windows config's data, templates/arcs,
PyMorph and standard-star paths before running.

These are Windows command equivalents; they have not been executed on the
Windows host from this Mac.

```powershell
$Dap = 'C:/Users/maxhi/CRD_DAP'
$Config = "$Dap/config/8143-1902.py"
$ReductionManifest = 'C:/Users/maxhi/CRD_DRP/validation_products/03_reduction_package_red_zap/CRD_DRP_reduction_manifest.json'

$S01Name = '8143-1902_scalar_redzap_S01'
$S02Name = '8143-1902_scalar_redzap_S02'
$S03Name = '8143-1902_scalar_redzap_S03'
$Script1Run = "$Dap/runs/$S01Name"
$Script2Run = "$Dap/runs/$S02Name"
$Script3Run = "$Dap/runs/$S03Name"
$Script04Dir = "$Dap/runs/8143-1902_global04"
$Script05Dir = "$Dap/runs/8143-1902_exact05"

Set-Location $Dap
```

### 01 — Prepare and register cubes

```powershell
conda run --no-capture-output -n crd_dap python "$Dap/scripts/01_prepare_and_register_cubes.py" `
  --config "$Config" --reduction-manifest "$ReductionManifest" --run-name "$S01Name"
```

### 02 — Master BL bins and spatial noise scalars

```powershell
conda run --no-capture-output -n crd_dap python "$Dap/scripts/02_make_master_BL_bins.py" `
  --config "$Config" --script1-run "$Script1Run" --run-name "$S02Name"
```

### 03 — RH3 likelihood cubes and residual noise scalars

```powershell
conda run --no-capture-output -n crd_dap python "$Dap/scripts/03_build_RH3_likelihood_cubes.py" `
  --config "$Config" --script1-run "$Script1Run" --script2-run "$Script2Run" --run-name "$S03Name" `
  --workers 7
```

To resume the same interrupted run:

```powershell
conda run --no-capture-output -n crd_dap python "$Dap/scripts/03_build_RH3_likelihood_cubes.py" `
  --config "$Config" --script1-run "$Script1Run" --script2-run "$Script2Run" --run-name "$S03Name" `
  --workers 7 --resume
```

Seven workers is an example for the eight-core Windows machine; adjust for
available RAM and responsiveness.

### 04 — Coherent two-disk models

Use the same config/review steps as above, with Windows input paths.

**A. Measure both PSF candidates.**

```powershell
conda run --no-capture-output -n crd_dap python "$Dap/scripts/04_fit_global_two_disk_model.py" `
  --config "$Config" --output-dir "$Script04Dir" --phase psf
```

Inspect PSF QA, then set the reviewed profile; this is an example:

```powershell
$PsfProfile = 'moffat'
```

**B. Run the twelve-fit smoothing calibration.**

```powershell
conda run --no-capture-output -n crd_dap python "$Dap/scripts/04_fit_global_two_disk_model.py" `
  --config "$Config" --output-dir "$Script04Dir" --script1-run "$Script1Run" `
  --script2-run "$Script2Run" --script3-run "$Script3Run" --psf-profile "$PsfProfile" --phase calibrate `
  --workers 4
```

**C. Inspect `figures/smoothing_calibration_2x6.pdf`, then review/select strength.**

```powershell
conda run --no-capture-output -n crd_dap python "$Dap/scripts/04_fit_global_two_disk_model.py" `
  --config "$Config" --output-dir "$Script04Dir" --script1-run "$Script1Run" `
  --script2-run "$Script2Run" --script3-run "$Script3Run" --psf-profile "$PsfProfile" --phase review `
  --interactive --workers 4 --resume
```

The interactive menu offers **1: weaker range**, **2: stronger range**, or **3: adopt a tested coefficient**. Weaker/stronger divides/multiplies the current baseline by four and launches another twelve-entry comparison; inspect its new sheet and run this review command again. Overlapping coefficients are reused, subject to the shared-candidate refinement check. Selecting a coefficient does not automatically launch the eight-model comparison.

**D. After selecting strength, run the eight-model comparison.**

```powershell
conda run --no-capture-output -n crd_dap python "$Dap/scripts/04_fit_global_two_disk_model.py" `
  --config "$Config" --output-dir "$Script04Dir" --script1-run "$Script1Run" `
  --script2-run "$Script2Run" --script3-run "$Script3Run" --psf-profile "$PsfProfile" --phase compare `
  --workers 4 --resume
```

**E. Inspect the comparisons and choose a model.** The name below is an example.

```powershell
conda run --no-capture-output -n crd_dap python "$Dap/scripts/04_fit_global_two_disk_model.py" `
  --config "$Config" --output-dir "$Script04Dir" --script1-run "$Script1Run" `
  --script2-run "$Script2Run" --script3-run "$Script3Run" --psf-profile "$PsfProfile" `
  --phase choose-model --model-name full_psf1_smooth1 `
  --note "Replace with your reason for choosing this model"
```

### 05 — Exact RH3 refits and kinematic maps

Run with the activated environment in a terminal for interactive retries.

```powershell
conda activate crd_dap
python "$Dap/scripts/05_extract_global_RH3_solution.py" `
  --config "$Config" --script4-run "$Script04Dir" --output-dir "$Script05Dir" `
  --script1-run "$Script1Run" --script2-run "$Script2Run" --script3-run "$Script3Run" `
  --workers 4
```

Use the same output directory with `--resume` to continue saved work.

## Updated Script 04 search

The shared-start refinement and log-header update changes the search fingerprint. Use a fresh Script 04 output directory, run its PSF phase and repeat the Script 04 commands above; the existing Scripts 01-03 runs remain the inputs. Keep old runs intact and do not edit their fingerprints. Within the new run, `--resume` reuses matching completed fits while allowing additional refinement when a newly available compatible state scores better. `SCRIPT04_SHARED_REFINEMENT_PASSES` is a config numerical budget (default 3); unfinished checks are flagged. Candidate sharing never crosses radial-domain/PSF/node-grid boundaries.

## Available arguments

| Script | Arguments |
|---|---|
| 01 | Required: `--config`. Handoff: `--reduction-manifest`. Optional: `--run-name`, `--skip-reduction-hash-check`. |
| 02 | Required: `--config`. Optional: `--script1-run`, `--run-name`. Explicit upstream paths are recommended. |
| 03 | Required: `--config`. Optional: `--script1-run`, `--script2-run`, `--run-name`, `--workers`, `--resume`. |
| 04 | Required: `--config`, `--output-dir`. Options: `--phase`, `--script1-run`, `--script2-run`, `--script3-run`, `--psf-profile`, `--workers`, `--resume`, `--baseline-coefficient`, `--coefficient`, `--interactive`, `--model-name`, `--note`. |
| 05 | Required `--config`, `--script4-run`, `--output-dir`. Optional `--script1-run`, `--script2-run`, `--script3-run`, `--workers`, `--resume`, `--non-interactive`. Upstream overrides are useful after relocating runs. |

- **Script 04 phases:** `psf`, `calibrate`, `review`, `compare`, `choose-model`.
  Except for `psf`, all three upstream-run paths and `--psf-profile` are currently
  required by the workflow. `review` needs `--interactive` or `--coefficient`;
  `compare` needs a saved strength selection; `choose-model` needs `--model-name`.
- **Change the calibration range:** use the calibration call with
  `--baseline-coefficient VALUE --resume`.
- **Noninteractive selection:** use the review call with
  `--coefficient VALUE --note "reason"` instead of `--interactive`.
  The selected effective coefficient `C=lambda*A` must already have a matching
  reviewed calibration fit. Zero selects the unpenalized control.
- **Resume:** use the original inputs and run name/output directory. Script 03
  checkpoints completed bins; Script 04 checkpoints completed whole fits.
  Scripts 01/02 have no `--resume` or `--workers` switch.
- **Fresh run after changes:** use new output/run names for changed inputs or
  science settings. Changing stage-04 bounds, PSF or search settings requires
  a new stage-04 directory.
- **Integrity option:** `--skip-reduction-hash-check` is available only in
  Script 01. It does not bypass Script 02's noise-manifest checks and is omitted
  from the default calls.
- **Config settings:** the intrinsic speed ceiling (400 km/s), MaNGA inputs,
  inclination inputs, standard-star paths and optimizer budgets are in the
  target config, not additional command-line flags.
- **Full help:** append `--help` to the relevant script call.

The Script 05 handoff is `products/RH3_global_kinematics.ecsv` and the associated per-bin fraction profiles. Scripts 06–07 are planned.
