# CRD_DAP Diagnostic Figure Catalogue

## Purpose

This file documents every planned pipeline diagnostic and publication-candidate figure. It is intentionally verbose. The diagnostic plots are a major part of the scientific quality control: a successful terminal exit code is not sufficient evidence that the analysis is trustworthy.

For every implemented plot, this document should record:

- exact filename;
- producing science script;
- helper function in `crd_utils.plotting`;
- scientific purpose;
- mathematical quantity shown;
- input products;
- panel layout;
- axes and units;
- color normalization and clipping;
- symbols/overlays;
- what a healthy result looks like;
- warning signs;
- recommended action when a warning appears;
- associated numerical files and log messages;
- whether the plot is QC-only, an analysis diagnostic, or a publication candidate.

If a plotting implementation changes, update this file at the same time.

---

# Conventions used throughout diagnostic maps

Unless a plot explicitly states otherwise:

- spatial coordinates should be in arcsec relative to the adopted KCWI center;
- the same image orientation should be used across all BL/RH3 maps;
- the signed kinematic-major-axis direction should be indicated when component identity matters;
- PowerBin maps should render the physical bin polygons/memberships rather than misleadingly interpolating them into smooth continuous fields;
- missing/flagged bins should have a visually distinct masked appearance rather than being assigned a numerical color;
- important color limits and percentiles should be saved in a JSON metadata sidecar when `SAVE_PLOT_METADATA_JSON=True`.

---

# Script 1 — cube preparation, calibration, registration, PSF/LSF, and noise

## `BL_collapsed_continuum.png`

**Script:** 01  
**Helper:** planned `plot_collapsed_continuum()`  
**Class:** QC + potential methods figure panel

### Purpose

Shows the robust two-dimensional BL continuum surface-brightness distribution obtained by collapsing clean wavelength samples. It verifies that the stacked cube contains the expected galaxy morphology and provides the primary visual context for center estimation and registration.

### Quantity

Conceptually,

\[
I_{\mathrm{BL}}(x,y)
=
\operatorname{median}_{\lambda\in\mathrm{good}}
F_{\mathrm{BL}}(x,y,\lambda),
\]

or the final adopted weighted equivalent.

### Overlays

- brightest-positive-spaxel initializer;
- continuum centroid;
- PyMorph center if transformed reliably;
- IFU footprint/boundary;
- masked spatial pixels.

### Healthy result

Galaxy light is spatially coherent, the peak/centroid are plausible, and there are no isolated extreme pixels dominating the morphology.

### Warning signs

- single-spaxel spike much brighter than surrounding continuum;
- strong negative/positive striping;
- galaxy unexpectedly truncated or offset;
- center estimates differing by more than a meaningful fraction of the PSF.

### Action

Inspect DQ/masks, wavelength-collapse range, stacking, WCS, and center estimation before continuing.

### Numerical/log companions

- collapsed image array;
- center comparison table;
- count/fraction of valid spatial pixels.

---

## `RH3_collapsed_continuum.png`

Same philosophy as the BL continuum image, but for RH3. It is especially important for checking whether the CaT-region continuum morphology and center are spatially consistent with BL after registration.

Healthy BL and RH3 images need not have identical surface-brightness distributions because the passbands differ, but their astrometric structure should be mutually plausible.

---

## `BL_RH3_registration.png`

**Script:** 01  
**Class:** critical QC

### Purpose

Demonstrates that BL and RH3 have been placed onto the same physical sky coordinate system before a single BL-defined PowerBin map is transferred to RH3.

### Suggested panels

1. BL collapsed continuum;
2. RH3 collapsed continuum transformed into BL coordinates;
3. overlay or normalized difference image.

### Overlays

- adopted common center;
- measured registration shift \((\Delta x,\Delta y)\);
- common IFU overlap region.

### Healthy result

Galaxy morphology aligns to within a small fraction of the delivered PSF or the documented registration uncertainty.

### Warning signs

- coherent dipole residual in the difference image;
- center offset comparable to a PowerBin dimension;
- WCS rotation/scale mismatch;
- non-overlapping regions later treated as shared apertures.

### Action

Do not proceed to Script 2 until registration is resolved.

---

## `geometry_center_comparison.png`

**Script:** 01 / updated in 04  
**Class:** QC

### Purpose

Compares independent center estimates:

- brightest continuum spaxel;
- continuum centroid/photometric fit;
- transformed PyMorph center, if available;
- final Script-4 kinematic center once known.

### Interpretation

Agreement supports a stable geometry initializer. A significant mismatch may indicate dust, asymmetric morphology, registration problems, or a genuine offset between photometric and kinematic centers.

### Action

Large disagreements are logged and should inform center bounds in Script 4 rather than being silently averaged.

---

## `BL_valid_spaxels.png` / `RH3_valid_spaxels.png`

**Script:** 01  
**Class:** QC

### Purpose

Maps which spatial elements survive hard data-quality requirements and the fraction of usable wavelengths per spaxel.

### Recommended display

Either two panels per arm:

- binary usable/unusable map;
- continuous good-wavelength fraction.

### Warning signs

Large contiguous missing regions, slice-like patterns, or a surprising number of rejected central spaxels.

### Action

Inspect DRP mask interpretation and `MIN_GOOD_WAVELENGTH_FRACTION` before PowerBin.

---

## `BL_bad_wavelength_fraction.png` / `RH3_bad_wavelength_fraction.png`

**Script:** 01  
**Class:** QC

### Purpose

Shows, as a function of wavelength, the fraction of spatial samples flagged unusable. This distinguishes isolated bad spaxels from wavelength channels that are globally problematic.

### Axes

- x: wavelength, Å;
- y: fraction of spatial samples bad, 0–1.

### Overlays

- hard channel-rejection threshold;
- instrument-good wavelength boundaries;
- major known emission/sky regions for context if useful.

### Warning signs

Broad sections with unexpectedly high rejection fraction or sharp features coincident with reduction artifacts.

---

## `BL_LSF.png` / `RH3_LSF.png`

**Script:** 01  
**Class:** critical QC + methods/publication candidate

### Purpose

Documents the empirical instrumental LSF measured from the required master arc.

### Quantity

Measured arc-line FWHM or Gaussian \(\sigma\) as a function of wavelength, optionally separated by slice/spatial position.

### Recommended display

- measured individual arc-line widths;
- smooth adopted LSF model;
- nominal instrumental-resolution expectation as a secondary reference;
- XSL template LSF overlaid in compatible units.

### Warning signs

- large unexplained scatter between neighboring lines;
- strong spatial dependence ignored by the adopted model;
- template width exceeding the data width over important wavelengths;
- poor smooth model residuals.

### Action

Do not trust derived stellar dispersions until the cause is understood. If the template is intrinsically broader than RH3 over critical regions, reconsider template/data matching strategy.

---

## `LSF_spatial_variation.png`

**Script:** 01  
**Class:** QC

### Purpose

Tests whether one wavelength-only LSF is adequate or whether the master arc shows meaningful slice/spatial variation.

### Display

FWHM residual relative to the wavelength-only mean/model versus slice or spatial coordinate for representative wavelengths.

### Healthy result

Variation is small relative to the precision required for \(\sigma_A,\sigma_B\).

### Warning sign

Coherent slice-dependent resolution shifts comparable to the astrophysical dispersion uncertainty.

---

## `BL_PSF.png` / `RH3_PSF.png`

**Script:** 01  
**Class:** QC

### Purpose

Documents the effective delivered spatial resolution separately for both arms.

### Important comparison

The difference

\[
\Delta\mathrm{FWHM}_{\mathrm{PSF}}
=
\mathrm{FWHM}_{\mathrm{BL}}
-
\mathrm{FWHM}_{\mathrm{RH3}}
\]

informs whether Script 6 needs explicit PSF-aware transfer of the RH3 velocity model to BL.

---

## `BL_RH3_PSF_comparison.png`

**Script:** 01  
**Class:** QC

### Purpose

Side-by-side or overlaid comparison of the two arm PSFs on the same angular scale.

### Warning sign

A large difference relative to typical PowerBin dimensions, which means “same aperture” does not imply the same effective spatial weighting.

---

## `noise_normalized_residuals_BL.png` / `noise_normalized_residuals_RH3.png`

**Script:** 01, possibly updated after initial pPXF fits  
**Class:** critical QC

### Purpose

Tests whether formal uncertainties place normalized residuals on the expected scale.

### Quantity

\[
z(\lambda)=\frac{F(\lambda)-M(\lambda)}{\sigma_F(\lambda)}.
\]

### Display

Histogram and/or quantile comparison to a unit-width Gaussian after obvious model-mismatch regions are handled appropriately.

### Warning signs

- width far from unity;
- strong heavy tails;
- asymmetric distribution;
- wavelength-correlated residuals.

### Action

Revisit variance scaling and covariance characterization before using \(e^{-\Delta\chi^2/2}\) as a relative-likelihood mapping.

---

## `spectral_covariance_BL.png` / `spectral_covariance_RH3.png`

**Script:** 01  
**Class:** critical QC

### Purpose

Shows preliminary wavelength-pixel correlation introduced by cube reconstruction, interpolation, and stacking.

### Suggested display

Correlation coefficient as a function of spectral lag, for diagnostic inspection only.

### Healthy result

The measured correlation is reported clearly. The production likelihood uses scalar noise calibration and does not fit spectral lag coefficients; its confidence widths require separate injection/recovery validation.

---

# Script 2 — master PowerBins

## `master_bins.png`

**Script:** 02
**Class:** major QC + publication candidate

### Purpose

Shows the physical BL-defined PowerBin tessellation used by every later stage.

### Overlays

- bin boundaries;
- adopted center;
- signed \(PA_{\mathrm{kin}}\) axis;
- optional effective radius;
- optional initial/updated \(2\sigma\) locations;
- IFU footprint.

### Healthy result

Bins are compact, follow the usable footprint, and become larger only where lower surface brightness requires it.

### Warning signs

Very elongated/nonphysical bins, tiny disconnected islands, or bins extending across masked gaps in a way that would make one spectrum spatially incoherent.

---

## `BL_SN_per_bin.png` / `RH3_SN_per_bin.png`

**Script:** 02
**Class:** QC

Separate convenience versions of achieved S/N maps.

---

## `BL_RH3_SN_comparison.png`

**Script:** 02  
**Class:** major QC + possible paper figure

### Purpose

Directly compares achieved S/N in the *same physical PowerBins* between the BL and RH3 arms.

### Panels

- left: BL S/N;
- right: RH3 S/N.

### Color normalization

Both panels share exactly the same `vmin` and `vmax`.

Default upper limit:

\[
v_{\max}
=P_{95}\left(\{S/N_{\mathrm{BL}}\}\cup\{S/N_{\mathrm{RH3}}\}\right).
\]

Values above the 95th percentile saturate but remain numerically preserved in data tables.

### Interpretation

Equivalent colors mean equivalent S/N. This makes it immediately obvious whether RH3 has substantially higher information content in the BL-defined bins.

---

## `bin_area_map.png`

**Script:** 02  
**Class:** QC

Maps physical bin area or number of member spaxels. This helps interpret later population/kinematic resolution and identify bins where intra-bin velocity shear may become important.

---

# Script 3 — RH3 scalar noise calibration and profile likelihood

After Script 3 completes, inspect `products/RH3_scalar_noise_calibration.npz`: `raw_scale`, `adopted_scale`, and `normalized_residual` cover every PowerBin. The raw and adopted scales are identical; values below one are retained. Review unusually large or small scales and structured best-fit residuals before interpreting likelihood widths. The full grid was fitted with the Script-2 spatially corrected uncertainty, and every saved chi-square was divided by the square of the residual scale. The saved `noise` in `RH3_log_spectra_and_local_best_fits.npz` is the final likelihood denominator; `noise_before_residual_scale` records what pPXF received during the grid fit.

## `RH3_single_component_velocity.png` and `RH3_single_component_sigma.png`

**Script:** 03  
**Class:** QC

These maps reveal large velocity/dispersion discontinuities, edge values and low-information regions in the independent one-component control fits. They are not disk-A/B estimates.

## `RH3_likelihood_bins/bin_XXXX.png`

**Script:** 03  
**Class:** per-bin fit diagnostic

Inspect the local two-component likelihood surface and best spectral model for grid-edge minima, failed states, broad or disjoint basins, and coherent stellar-feature residuals. The physical disk labels are assigned later by the global model.

## `RH3_local_one_vs_two_delta_chi2.png`

**Script:** 03  
**Class:** QC

The difference in total scalar-calibrated diagonal chi-square between one- and two-component fits is an internal comparison, not a formal detection significance. Reliability and interval coverage require injection/recovery tests. Spectral lag correlation measured in Script 1 remains a diagnostic; it is not included in the scalar likelihood.

---

# Script 4 — radial geometry and global two-disk model

The implemented products below replace the earlier proposed Script-4 figure
list. Commands and review gates are documented in [SCRIPT04.md](SCRIPT04.md).
All figures have PNG previews and single-page PDF companions.

## stellar_PSF_QA.pdf

Each Feige 24 exposure has a continuum image, Gaussian and Moffat residual maps,
and angular light profiles. Map axes are native angular coordinates, respecting
rectangular pixels. Check centering, residual structure, flags and differences
between exposures. RMS is descriptive, not a formal Gaussian/Moffat test.
Profile selection remains explicit. The transferred PSF does not replay the
exact science-stack drizzle response, and standards were taken hours earlier.

## full_radial_geometry.pdf and two_sigma_radial_geometry.pdf

The continuum map shows reference ring ellipses and center. The companion
histogram shows selected member-spaxel radii and records bin count, extent and
the provisional disk-label source. Check node density relative to the measured
larger-axis FWHM and inspect boundary-crossing bins. Membership remains fixed
during optimization.

## smoothing_calibration_2x6.pdf

One page, exactly two rows of six entries: control plus eleven ordered strengths.
Each entry contains:

- signed A/B intrinsic rotation curves, with faint control curves;
- a per-bin data-cost-change map relative to that same control;
- actual geometry, changes from control, C, lambda, Q_data, penalty, RMS curve
  distances and fit flags.

All maps share a zero-centered diverging color scale: blue improves the data
fit, red worsens it. The color scale never includes the smoothing penalty.
Adjacent-strength RMS changes appear in a strip above the entries. Curve
distances exclude the imposed central anchor and are descriptive, not errors,
significance levels or quality ranks. Open the PDF and zoom to inspect text.

The control should have the lowest unpenalized data cost if the searches have
found their appropriate minima. CONTROL_SEARCH_NEEDS_REFINEMENT means a
penalized search found a better data score; compare optimizer alternatives and
budgets before interpreting the sweep. No best strength is automatically chosen.

## model_diagnostics.pdf for each of the eight named models

The actual filenames are, for example, full_psf1_smooth1_diagnostics.pdf.
Six panels show radial curves, predicted A/B LOS maps, coarse profiled fraction,
per-bin data Delta-chi-square from independent minima, and search history.
Search time versus objective is numerical progress, not posterior sampling.
The title retains convergence, boundary, PSF-transfer and coverage flags.

## extent_comparison.pdf

Four panels pair restricted/full curves for each PSF/regularity setting. Compare
common radii before assessing outer behavior. The accompanying
products/eight_model_comparison.json reports common-inner-bin data sums and
means, separate outer contributions, parameter counts and geometry changes.
Different bin sets' total scores cannot rank model evidence. No uncertainty
bands or automatic preferred model are implied.

## Machine-readable products and logs

Each fits/model_name.json retains the joint parameter vector, radial grid,
per-bin data scores, coarse fractions, shear, radial-padding light fractions,
alternative refined starts, objective history and convergence flags. Prediction
ECSV tables are saved in products/. A selected_global_model.json is written only
after explicit model selection, identifying the Script-05 input velocities.

The main log summarizes inputs, dimensions and completed fits. Separate
logs/model_name.progress.log files record live search data cost, penalty,
geometry and elapsed time. The scientific input fingerprint, config snapshot,
PSF choice, strength sweeps and selection rationale are retained. Resume reuses
only matching completed whole fits.

ROTATION_SPEED_CEILING identifies an intrinsic speed near the configurable
400 km/s maximum; ROTATION_SPEED_BOUNDARY also includes zero. GEOMETRY_BOUNDARY
identifies a geometry limit. Non-convergence flags do not invalidate the software
workflow, but they prevent claiming a settled optimum. PSF edge and outer
constant-speed padding diagnostics describe where the spatial approximation
matters; they do not create additional data or remove those bins.

---

# Script 5 — exact RH3 extraction, basins, grid validation, one/two-component reliability

## `RH3_selected_basin.png`

**Script:** 05  
**Class:** major methods diagnostic

### Purpose

Explains which 3-D likelihood basin is associated with the global XookSuut-style anchor.

### Recommended panels

1. profiled \(\Delta\chi^2(V_A,V_B)\);
2. basin membership projected into velocity space;
3. best \(f_A\) surface;
4. markers for independent cube global minimum, global-model anchor, and selected basin minimum.

### Interpretation

The selected basin need not contain the absolute individual-cube minimum. A slightly worse local basin may be selected because it is globally coherent across the galaxy.

---

## `RH3_basin_likelihood_concentration.png`

**Script:** 05  
**Class:** QC/methods

### Purpose

Separates *basin topology* from *statistical concentration*.

### Quantity

For basin cell \(c\),

\[
p_c
\propto
\exp(-\Delta\chi_c^2/2).
\]

### Suggested display

Cumulative normalized likelihood mass versus number/fraction of highest-weight basin cells, with reference levels at 90%, 95%, 99%, 99.5%, and 99.9%.

### Important interpretation

A basin may contain the whole cube while 99% of its likelihood is concentrated in a small subset. A large basin alone does not mean RH3 is uninformative.

### Scientific role

The listed mass levels are diagnostics only. They are not hard candidate-selection cutoffs.

---

## `RH3_cube_global_vs_selected_minimum.png`

**Script:** 05  
**Class:** critical QC

### Quantity

\[
\Delta\chi^2_{\mathrm{cube},i}
=
\chi_i^2(\mathrm{selected})
-
\min_{V_A,V_B,f_A}\chi_i^2.
\]

### Purpose

Tests whether the global smooth-disk constraint is choosing a low-likelihood solution that the individual spectrum strongly rejects.

### Healthy result

Many bins near zero; some modestly positive bins are acceptable when the global model chooses a nearby coherent basin.

### Warning sign

Large positive values in numerous bins, especially coherent spatial regions.

---

## `RH3_selected_delta_chi2_map.png`

**Script:** 05  
**Class:** QC

Maps the selected solution's \(\Delta\chi^2\) relative to the minimum of its own selected basin. Distinguishes “wrong basin relative to global cube” from “far up the wall of the chosen basin.”

---

## `RH3_grid_edge_distance.png`

**Script:** 05  
**Class:** critical QC

### Quantity

Minimum number of grid cells between the selected discrete state and any \(V_A\), \(V_B\), or \(f_A\) boundary.

### Flag

`GRID_EDGE_WARNING` when the distance is ≤2 cells.

### Action

Expand and rerun the relevant grid rather than extrapolating beyond it.

---

## `RH3_exact_vs_grid_solution.png`

**Script:** 05  
**Class:** method-validation QC

### Purpose

Compares coarse-grid/interpolated quantities to the exact pPXF refit at the continuous global velocities.

### Suggested maps

- \(V_A^{\mathrm{exact}}-V_A^{\mathrm{grid}}\);
- \(V_B^{\mathrm{exact}}-V_B^{\mathrm{grid}}\);
- \(\sigma_A^{\mathrm{exact}}-\sigma_A^{\mathrm{grid}}\);
- \(\sigma_B^{\mathrm{exact}}-\sigma_B^{\mathrm{grid}}\);
- \(f_A^{\mathrm{fine}}-f_A^{\mathrm{coarse}}\).

### Interpretation

Small differences validate that the cube is adequate for likelihood topology while the final exact refit removes quantization.

---

## `RH3_fA_refinement.png`

**Script:** 05  
**Class:** per-bin QC

Plots \(\chi^2(f_{A,\mathrm{RH3}})\) at the exact global velocities across the fine 0.01-step local fraction scan.

### Warning signs

Minimum at the refinement interval edge or a very broad/flat profile.

### Action

Automatically widen the fine interval if the edge is hit.

---

## `two_component_mock_calibration.png`

**Script:** 05 / model calibration  
**Class:** critical methods diagnostic

### Purpose

Shows the mock-data calibration connecting the one-vs-two-component statistic and recovery behavior to physical/observational conditions.

### Suggested content

Separate panels or clearly documented variants for:

- single-component false-positive distribution of \(T\);
- two-component recovery fraction versus \(\Delta V\);
- recovery versus \(f_A\);
- recovery versus S/N;
- interpolation region in \((S/N,\Delta V,f_A,\sigma_A,\sigma_B)\) space.

### Important

Do not label the empirical recovery fraction as the posterior probability that an observed spectrum “is two-component.”

---

## `two_component_false_positive_map.png`

**Script:** 05  
**Class:** science-quality map

Maps

\[
p_{\mathrm{false},i}
=
P(T\ge T_i\mid\mathrm{true\ one\ component})
\]

from matched single-component simulations.

Low values mean a true single-component spectrum rarely creates as strong a two-component preference.

---

## `two_component_recovery_probability_map.png`

**Script:** 05  
**Class:** science-quality map

Maps empirical two-component recovery/completeness under conditions comparable to the selected bin solution.

Low recovery does not necessarily mean the observed bin is one-component; it means the decomposition lies in a regime where the method is not reliably complete.

---

# Script 6 — BL stellar populations and RH3-likelihood propagation

## `BL_primary_spectral_fit.png`

**Script:** 06  
**Class:** major per-bin science diagnostic

### Purpose

Shows the RH3-anchored BL population fit at exact best RH3 kinematics.

### Panels

- observed BL spectrum;
- total model;
- Disk A stellar contribution;
- Disk B stellar contribution;
- gas contribution where present;
- residuals.

### Interpretation

The two stellar model components should combine smoothly without one component being used to fit isolated noise spikes. Residuals should not show broad template mismatch.

---

## `BL_joint_profile_spectral_fit.png`

Same display as the primary fit, but using the deterministically refined minimum of the joint RH3+BL profile statistic.

Comparing this plot with `BL_primary_spectral_fit.png` helps determine whether allowing RH3-supported kinematic movement changes the population decomposition visibly.

---

## `BL_emission_model.png`

**Script:** 06  
**Class:** QC

### Purpose

Shows fitted gas emission separately from the stellar components. Particularly important around Hβ and other blue emission features that could bias age-sensitive absorption features.

### Important

Gas-template weights are excluded from \(f_{A,\mathrm{BL}}\).

---

## `BL_fit_residuals.png`

**Script:** 06  
**Class:** critical QC

Residual spectrum, ideally with normalized residuals and masked regions indicated. Used to diagnose template mismatch, insufficient polynomial treatment, emission residuals, and noise-scale problems.

---

## `BL_fA_profile.png`

**Script:** 06  
**Class:** major inference diagnostic

### Purpose

Shows the profile joint statistic or BL statistic versus explicit stellar \(f_{A,\mathrm{BL}}\) for the anchored and/or joint-profile solution.

### Coarse/fine display

Plot coarse 0.1-step evaluations and fine 0.01-step refinement distinctly so the viewer can see that the final fraction is not a coarse-grid artifact.

---

## `BL_fA_marginalized_distribution.png`

**Script:** 06  
**Class:** major inference diagnostic

### Purpose

Shows the RH3-likelihood-propagated support for \(f_{A,\mathrm{BL}}\) after integrating/reweighting over the sampled RH3-supported kinematic family.

### Report

- mode;
- median;
- 16th/84th percentiles;
- anchored best fraction;
- joint-profile best fraction.

### Terminology

Use “likelihood-propagated distribution” or equivalent unless the implementation has become a formal Bayesian posterior.

---

## `RH3_BL_joint_likelihood.png`

**Script:** 06  
**Class:** methods/science diagnostic

### Purpose

Visualizes how BL information modifies the RH3-supported solution family.

### Suggested view

Profiled \((V_A,V_B)\) RH3 contours with markers for:

- RH3 global anchor;
- RH3 basin minimum;
- best sampled joint state;
- deterministically refined joint profile state.

Optional joint-likelihood contours may be shown only with precise documentation of how the RH3 sampling/reweighting is represented.

---

## `RH3_to_joint_kinematic_shift.png`

**Script:** 06  
**Class:** robustness QC

Maps the differences between exact RH3-anchored kinematics and the joint-profile kinematics selected after BL information is allowed to act within the RH3-supported family.

Large coherent shifts deserve inspection.

---

## `joint_profile_refinement_path.png`

**Script:** 06  
**Class:** optimizer QC

### Purpose

Shows the deterministic local descent path used to find the exact discrete joint-profile minimum starting from the best sampled state.

### Display

On a suitable 2-D projection of the RH3 basin, mark:

- starting sampled cell;
- successive lower-\(J\) neighboring cells;
- final local joint minimum.

### Interpretation

A short path means the sampling already found the neighborhood well. A very long path may indicate a broad shallow valley or insufficient sampling coverage.

---

## `BL_sampling_convergence.png`

**Script:** 06  
**Class:** critical QC

### Purpose

Demonstrates that likelihood-propagated population summaries stabilize as more RH3-weighted draws are added.

### x-axis

Cumulative number of RH3 draws.

### y-series

Separate panels or figures for:

- \(\log Age_A\), \(\log Age_B\);
- \([M/H]_A\), \([M/H]_B\);
- \(f_{A,\mathrm{BL}}\);
- \(\Delta\log Age\);
- \(\Delta[M/H]\).

Show median and optionally 16th/84th percentiles.

### Convergence markers

Mark each batch check as pass/fail and the point where the required number of consecutive passes is achieved.

### Critical interpretation

Each plotted point uses **all draws accumulated up to that point**, not just the newest batch.

---

## `BL_sampling_effective_sample_size.png`

**Script:** 06  
**Class:** sampling QC

### Purpose

Tracks effective sample size as BL likelihood reweights the RH3-drawn state ensemble.

A very low effective sample size means the final joint result is dominated by a small number of sampled states and may require more draws or broader direct evaluation.

---

## `BL_unique_cells_vs_draws.png`

**Script:** 06  
**Class:** computational + likelihood-concentration diagnostic

### x-axis

Cumulative RH3 draws.

### y-axis

Number of unique RH3 basin cells encountered/evaluated.

### Interpretation

Rapid plateau at a small number of unique cells indicates strongly concentrated RH3 support; nearly linear growth indicates a broad distribution.

---

## `BL_population_stability_age.png`

**Script:** 06  
**Class:** robustness diagnostic

Shows how inferred ages vary across the RH3-supported kinematic states and BL fractions, weighted or marked by joint support.

Purpose: demonstrate whether age conclusions are stable against the remaining RH3 kinematic uncertainty.

---

## `BL_population_stability_metallicity.png`

Same philosophy for metallicity.

---

## `BL_template_weight_A.png` / `BL_template_weight_B.png`

**Script:** 06  
**Class:** population QC/science diagnostic

### Purpose

Shows SSP weight distribution over age-metallicity space for each disk for selected fits.

### Warning signs

Weights pile against the oldest/youngest or highest/lowest metallicity template-grid boundary.

---

## `BL_template_grid_edge_flags.png`

**Script:** 06  
**Class:** critical QC map

Spatially maps bins with `AGE_GRID_EDGE_*` or `METALLICITY_GRID_EDGE_*` flags.

A coherent region of boundary hits may indicate the SSP library does not span the required populations.

---

# Script 7 — final science maps, uncertainties, robustness, and runtime

## `final_VA_map.png` / `final_VB_map.png`

**Script:** 07  
**Class:** publication science maps

Final globally coherent stellar velocity maps with uncertainties/quality masks incorporated.

---

## `final_sigmaA_map.png` / `final_sigmaB_map.png`

Final RH3 stellar-dispersion maps. High-shear/LSF-sensitive bins should be clearly flagged or masked according to the adopted quality policy.

---

## `final_deltaV_map.png`

Maps

\[
\Delta V=|V_A-V_B|.
\]

Useful both scientifically and for interpreting decomposition reliability.

---

## `final_fA_RH3_map.png`

Maps the refined RH3 stellar fraction. This fraction refers to the RH3 fitting band and should not be directly equated with the BL fraction.

---

## `final_age_A_map.png` / `final_age_B_map.png`

**Script:** 07  
**Class:** major publication science maps

Each PowerBin is assigned the inferred stellar age of that disk. The spatial resolution is the BL PowerBin, not the native spaxel scale.

---

## `final_metallicity_A_map.png` / `final_metallicity_B_map.png`

Same structure for \([M/H]\).

---

## `final_fA_BL_map.png`

Maps the stellar Disk A fraction over the adopted usable BL fitting range.

Comparison with `final_fA_RH3_map.png` may itself be scientifically informative because different stellar populations can have wavelength-dependent component light fractions.

---

## `disk_contrast_delta_age.png`

**Script:** 07  
**Class:** major publication science map

### Quantity

Prefer logarithmic age contrast when appropriate:

\[
\Delta\log Age
=
\log Age_A-\log Age_B.
\]

### Interpretation

Positive and negative values identify which disk is older. The sign convention must be written directly in the figure caption/metadata.

---

## `disk_contrast_delta_metallicity.png`

**Quantity**

\[
\Delta[M/H]
=
[M/H]_A-[M/H]_B.
\]

This is a direct physical contrast between the two stellar disks.

---

## `method_delta_age_A.png` / `method_delta_age_B.png`

**Script:** 07  
**Class:** robustness maps

### Purpose

These are *not* disk A-minus-B contrasts. They compare inference methods:

\[
\Delta Age_{A,\mathrm{method}}
=
Age_{A,\mathrm{likelihood-propagated}}
-
Age_{A,\mathrm{anchored}}.
\]

Same for Disk B.

### Interpretation

Near-zero maps show that propagating RH3 kinematic uncertainty does not materially change the age inference.

---

## `method_delta_metallicity_A.png` / `method_delta_metallicity_B.png`

Same method-robustness comparison for metallicity.

---

## `age_uncertainty_A.png` / `age_uncertainty_B.png`

**Script:** 07  
**Class:** uncertainty maps

Maps the selected statistical uncertainty summary (e.g. half-width or asymmetric interval information encoded in companion tables) derived from end-to-end MC realizations.

---

## `metallicity_uncertainty_A.png` / `metallicity_uncertainty_B.png`

Same for metallicity.

---

## `radial_age_profiles.png`

**Script:** 07  
**Class:** publication science figure

Radial profiles of Disk A and Disk B age using deprojected radius. Must account for spatial correlations when interpreting gradients; individual neighboring PowerBins are not independent measurements.

---

## `radial_metallicity_profiles.png`

Same for metallicity.

---

## `MC_convergence.png`

**Script:** 07  
**Class:** critical uncertainty QC

### Purpose

Demonstrates that final statistical uncertainty summaries stabilize with increasing number of end-to-end realizations.

### x-axis

Cumulative \(N_{\mathrm{MC}}\).

### y-series

Representative or aggregate 16th/50th/84th percentile summaries for the principal science quantities.

### Convergence

Use a documented batch/interval stability rule rather than selecting \(N_{\mathrm{MC}}\) solely because it is a round number.

---

## `systematic_robustness_summary.png`

**Script:** 07  
**Class:** major methods diagnostic

### Purpose

Summarizes how major science conclusions change under reasonable alternative analysis assumptions, such as:

- ring spacing;
- radial extent;
- inclination prior;
- PSF/LSF perturbation;
- polynomial degree;
- RH3/fraction grid refinement;
- mildly relaxed geometry;
- alternative SSP/template assumptions.

### Important

These are systematic/method-robustness shifts, not automatically random errors to be added in quadrature.

---

## `disk_age_difference_significance.png` / `disk_metallicity_difference_significance.png`

**Script:** 07  
**Class:** optional future science diagnostic

### Status

Do not implement as a simplistic binary significance map until the treatment of spatial correlation and multiple comparisons is explicitly decided.

Potential quantities include the fraction of MC/likelihood-propagated realizations satisfying a signed disk contrast, but interpretation must be documented carefully.

---

# Runtime / computation diagnostics

## `MC_runtime_scaling.png`

**Script:** benchmark utility / 07  
**Class:** computational planning

### Purpose

Shows measured or carefully labeled estimated wall time as a function of worker/core count.

### Inputs

Measured timings from one or more full MC realizations and, ideally, direct benchmark runs at multiple worker counts.

### Important

Do not present ideal linear scaling as measured performance. Distinguish measured points from forecasts.

---

# Standard plot metadata sidecars

For important plots, save an adjacent JSON file when enabled. Example fields:

```json
{
  "filename": "BL_RH3_SN_comparison.png",
  "target": "example",
  "run_id": "example_YYYYMMDD_HHMMSS",
  "vmin": 35.0,
  "vmax": 96.4,
  "vmax_definition": "95th percentile of combined BL+RH3 bin S/N",
  "n_bins": 152,
  "config_snapshot": "metadata/config_snapshot.py"
}
```

The metadata should store plot-defining numerical choices that cannot be recovered reliably by inspecting the PNG alone.

---

# Diagnostic review philosophy

Before a final science run is accepted, the user should be able to answer all of the following from the diagnostic suite:

1. Are BL and RH3 spatially registered?
2. Are the instrumental LSF and template LSF compatible and measured accurately enough for two-component dispersions?
3. Does the formal noise scale match observed residuals, and is spectral covariance handled?
4. Does the BL PowerBin tessellation make physical sense, and what RH3 S/N is achieved in the same bins?
5. Do representative RH3 likelihood surfaces genuinely contain meaningful two-component structure?
6. Is the RH3 velocity/fraction grid wide and fine enough?
7. Are the measured RH3 \(2\sigma\) peak radii and XookSuut-style radial nodes sensible?
8. Are non-parametric ring velocities constrained by adequate spatial coverage?
9. Does the global disk model fit the selected spectral likelihood basins without strongly overriding individual spectra?
10. Are centroid-only and bin-integrated velocities consistent, or is intra-bin shear important?
11. Do the \(2\sigma\)-limited and full-aperture global solutions agree over their common radial range?
12. Do exact continuous RH3 refits agree with the coarse-grid likelihood solution?
13. Are the one-vs-two-component results reliable under matched mock calibrations?
14. Do BL spectral residuals, gas fits, and template weights look physically reasonable?
15. Are BL fractions and stellar-population results stable when RH3 kinematic uncertainty is propagated?
16. Did likelihood-weighted Script-6 sampling converge, and is effective sample size adequate?
17. Are age/metallicity results limited by SSP-template boundaries?
18. Did the end-to-end Monte Carlo uncertainty calculation converge?
19. Are the principal science conclusions robust to important analysis choices and SPS assumptions?
20. Can every final figure/table be traced back to the configuration, run logs, and numerical products that produced it?

If the answer to an important item is “no” or “unknown,” the pipeline should make that visible rather than allowing a visually polished final map to hide the unresolved issue.
