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

$$
I_{\mathrm{BL}}(x,y)=\mathrm{median}_{\lambda\in\mathrm{good}}F_{\mathrm{BL}}(x,y,\lambda),
$$

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

1. BL registration continuum, displayed only over the common valid BL/RH3 footprint with a robust common-footprint-derived stretch;
2. RH3 registration continuum transformed into BL coordinates, on the same common valid footprint and with its own robust stretch;
3. normalized morphology difference image **only when the cross-correlation passes the registration QC**. If registration is inconclusive, this panel instead shows the common-valid-footprint mask and explicitly states that no residual image or numerical shift was adopted.

When a sufficiently wide common instrument-good wavelength interval exists, both registration images are collapsed over that same observed-frame interval and the interval is written into the figure title/metadata. Otherwise the diagnostic falls back to the full-arm continuum images.

### Overlays / annotations

- wavelength interval or full-arm mode used to construct the registration images;
- robust morphology-contrast statistic for each arm;
- measured residual shift $(\Delta x,\Delta y)$ when the cross-correlation is valid;
- an explicit `inconclusive` annotation rather than a numerical shift when the registration QC fails;
- the common valid footprint in the third panel for inconclusive cases, instead of an unstable normalized residual;
- robust image stretches measured only from common-valid pixels so narrow IFU/reprojection edge artifacts do not dominate the morphology display.

### Healthy result

Either (1) both images have sufficient contrast and the residual morphology shift is small relative to the delivered PSF / registration requirement, or (2) morphology cross-correlation is explicitly declared inconclusive while the independent WCS peak/centroid comparison remains mutually consistent. A low-contrast image is **not** evidence for a real astrometric offset.

The cross-correlation is also restricted to a configurable local residual-search radius around the WCS-predicted alignment. If the best peak reaches that boundary, the diagnostic is **inconclusive** rather than adopting the boundary value as an astrometric measurement. This guards against IFU-footprint edges or weak passband structure producing a spurious large shift.

### Warning signs

- coherent dipole residual in a valid normalized-difference image;
- center offset comparable to a PowerBin dimension;
- WCS rotation/scale mismatch;
- non-overlapping regions later treated as shared apertures;
- `REGISTRATION_INCONCLUSIVE`, which means the numerical cross-correlation shift must not be interpreted.

### Action

Resolve a genuine center/WCS disagreement before Script 2. If the morphology cross-correlation alone is inconclusive because one passband has insufficient contrast, use the independent sky-coordinate center agreement and appropriate external astrometric information rather than forcing the cross-correlation.

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

Documents the empirical instrumental LSF measured from the required master arc and its DRP wavelength/slice/position geometry maps. Before this figure is produced, Script 1 verifies that the calibration matches the science cube in camera, grating, slicer/IFU, binning, and central wavelength. Geometry-map discovery is also checked through FITS provenance rather than relying only on filename exposure numbers.

### Quantity

Measured arc-line FWHM or Gaussian $\sigma$ as a function of wavelength, optionally separated by slice/spatial position.

### Recommended display

- measured individual arc-line widths;
- smooth adopted LSF model drawn **only across the empirically supported accepted-line interval**;
- instrument-good `WAVGOOD0/1` boundaries;
- shaded edge regions that are instrument-good but lack direct accepted-line support;
- nominal instrumental-resolution expectation as a secondary reference, if useful;
- XSL template LSF overlaid in compatible units when available.

### Warning signs

- large unexplained scatter between neighboring lines;
- strong spatial dependence ignored by the adopted model;
- template width exceeding the data width over important wavelengths;
- poor smooth model residuals;
- a large gap between `WAVGOOD0/1` and the bluest/reddest accepted arc-line measurements.

### Action

Do not trust derived stellar dispersions until the cause is understood. If the template is intrinsically broader than RH3 over critical regions, reconsider template/data matching strategy. If `LSF_EMPIRICAL_COVERAGE_GAP` is raised, later fitting must mask wavelengths outside the empirical LSF support or introduce an independently validated LSF model; the polynomial must not be silently extrapolated.

---

## `BL_LSF_spatial_variation.png` / `RH3_LSF_spatial_variation.png`

**Script:** 01  
**Class:** QC

### Purpose

Tests whether one wavelength-only LSF is adequate or whether the master arc shows meaningful slice/spatial variation.

### Display

FWHM residual relative to the wavelength-only mean/model versus slice or spatial coordinate for representative wavelengths.

### Healthy result

Variation is small relative to the precision required for $\sigma_A,\sigma_B$. The plotted spatial metric should describe coherent slice/position offsets from the global wavelength-only LSF model, not simply the raw RMS of every individual arc-line residual. The latter is saved separately because it also contains line-fitting scatter.

### Warning sign

Coherent slice-dependent resolution shifts comparable to the astrophysical dispersion uncertainty.

---

## `BL_LSF_spatial_summary.png` / `RH3_LSF_spatial_summary.png`

**Script:** 01  
**Class:** critical QC + methods/publication candidate

### Purpose

Separates coherent field dependence of the LSF from the much larger scatter of individual arc-line fits. This is the preferred visual diagnostic for deciding whether a wavelength-only LSF is adequate.

### Panels

1. median fractional FWHM residual for each KCWI slice, with the 16th--84th percentile interval of accepted line measurements in that slice;
2. median fractional FWHM residual for each within-slice position bin, again with the 16th--84th percentile interval.

The title reports the RMS of the group medians used by the numerical QC metrics.

### Healthy result

Group medians remain close to zero and their RMS is small compared with the precision required for the stellar-dispersion measurements, even if the line-by-line scatter plot appears visually broad.

### Warning signs

- coherent monotonic or step-like offsets with slice ID;
- position-bin medians systematically displaced from zero;
- group-median RMS approaching the astrophysical dispersion precision requirement.

### Action

If coherent spatial variation is non-negligible, carry a slice/position-dependent LSF into later bin-level template matching rather than using one global wavelength-only curve.

---

## `BL_PSF.png` / `RH3_PSF.png`

**Script:** 01  
**Class:** QC

### Purpose

Documents the effective delivered spatial resolution separately for both arms.

### Important comparison

The difference

$$
\Delta\mathrm{FWHM}_{\mathrm{PSF}}=\mathrm{FWHM}_{\mathrm{BL}}-\mathrm{FWHM}_{\mathrm{RH3}}
$$

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

$$
z(\lambda)=\frac{F(\lambda)-M(\lambda)}{\sigma_F(\lambda)}.
$$

### Display

Histogram and/or quantile comparison to a unit-width Gaussian after obvious model-mismatch regions are handled appropriately.

### Warning signs

- width far from unity;
- strong heavy tails;
- asymmetric distribution;
- wavelength-correlated residuals.

### Action

Revisit variance scaling and covariance characterization before using $e^{-\Delta\chi^2/2}$ as a relative-likelihood mapping.

---

## `spectral_covariance_BL.png` / `spectral_covariance_RH3.png`

**Script:** 01  
**Class:** critical QC

### Purpose

Visualizes wavelength-pixel covariance introduced by cube reconstruction, interpolation, and stacking.

### Suggested display

Correlation coefficient as a function of spectral lag, and optionally a covariance/correlation matrix for representative spectra.

### Healthy result

Either covariance is negligible beyond a small lag, or a stable structure is measured and incorporated/approximated in pPXF fits and MC noise generation.

---

# Script 2 — master PowerBins

## `binning_aperture.png`

**Script:** 02  
**Class:** critical QC

### Purpose

Shows the distinction between the Script-1 hard-good data footprint and the smaller region of the BL cube that Script 2 considers to contain useful stellar light for adaptive binning. This prevents exposed blank sky from being accreted into enormous low-information PowerBins.

### Panels

- smoothed BL continuum-detection/significance proxy with the configured threshold and final aperture boundary;
- accepted aperture in tangent-plane coordinates with the adopted galaxy center.

### Important interpretation

The aperture threshold defines the *outer analysis domain*. It is not a native-spaxel S/N cut. Low-S/N spaxels inside the connected stellar-body aperture remain available to PowerBin.

### Healthy result

One connected region contains the galaxy center, follows the visible stellar body, excludes the large blank-sky footprint, and does not carve holes through plausible low-surface-brightness galaxy light.

### Warning signs

- aperture locks onto a remote artifact rather than the central galaxy;
- only the bright nucleus survives;
- aperture reaches the full exposed rectangle because sky structure is being interpreted as galaxy light;
- disconnected islands or detector/stack edges dominate the detection proxy.

### Action

Inspect the BL continuum image and Script-1 validity maps before changing the threshold. If the automatic aperture is not scientifically defensible, do not proceed merely because PowerBin can technically run.

---

## `master_bins.png`

**Script:** 02  
**Class:** major QC + publication candidate

### Purpose

Shows the physical BL-defined PowerBin tessellation used by every later stage. This is the authoritative spatial partition; RH3 is not independently binned.

### Overlays

- member pixels / bin identities;
- adopted center;
- signed $PA_{\mathrm{kin}}$ axis when available;
- IFU/analysis footprint through the absence of bins outside the accepted aperture.

### Healthy result

Bins are compact and connected, small in high-S/N regions, and grow gradually toward lower surface brightness. The tessellation should cover the useful stellar body without extending across excluded sky or masked gaps.

### Warning signs

Very elongated/nonphysical bins, tiny disconnected islands, extreme bin-size jumps driven by an incorrect aperture, or a bin layout that does not trace the visible galaxy.

### Numerical companions

- `master_bins.fits`;
- `master_bin_table.ecsv`;
- `master_bin_membership.npz`;
- PowerBin capacity and RMS fractional scatter in the log/manifest.

---

## `BL_RH3_bin_transfer.png`

**Script:** 02  
**Class:** critical cross-arm QC

### Purpose

Verifies that the single BL-defined physical tessellation has been transferred onto the red/RH3 native grid through the celestial WCS rather than recomputed independently.

### Panels

- BL master bin membership;
- RH3 native pixels labeled by transferred BL bin ID;
- distance between each assigned RH3 pixel center and the nearest BL member-pixel center in arcseconds.

### Healthy result

The red/RH3 membership reproduces the same sky apertures, bin IDs remain spatially coherent, and transfer distances are comfortably below the configured limit over nearly all candidate pixels.

### Warning signs

- large unassigned holes inside the BL-defined galaxy aperture;
- coherent edge offsets between arms;
- transfer distances clustering near the maximum allowed value;
- mismatched spatial scale/rotation making nearest-pixel membership inappropriate.

### Action

Do not independently re-PowerBin RH3 to make the plot look cleaner. Resolve the WCS/spatial-grid problem or implement an area-overlap aperture transfer if the native grids are genuinely different.

---

## `BL_SN_per_bin.png` / `RH3_SN_per_bin.png`

**Script:** 02  
**Class:** QC

### Purpose

Maps the robust achieved continuum S/N per spectral pixel for each final extracted spectrum. The plotted value is `median(flux) / median(uncertainty)` over the configured window and is shown only when the median continuum is positive. BL uses the configured BL binning/fitting window; RH3 uses its independently configured S/N-evaluation window but the **same BL-defined physical bins**.

### Healthy result

Most non-single BL bins lie near the target S/N, while single high-surface-brightness pixels may exceed it. RH3 S/N may differ substantially because it does not drive the tessellation. Bins with undefined achieved S/N are rendered in neutral gray inside the existing BL master-bin footprint, while sky outside the tessellation remains white. Gray therefore means "the physical bin exists, but this S/N measurement is undefined," not failed BL-to-RH3 membership transfer.

### Warning signs

A large population of BL bins far below the target can indicate a mismatch between the PowerBin S/N proxy and the final extracted-spectrum measurement, poor wavelength coverage, or unmodeled spatial covariance. A bin with an enormous signed or legacy S/N should be inspected with `scripts/inspect_script02_sn.py`; the production-facing map does not clip the spectrum but treats non-positive continuum as undefined positive S/N.

### Numerical companions

`master_bin_table.ecsv` preserves `BL_SN_SIGNED`, `RH3_SN_SIGNED`, `BL_SN_LEGACY`, `RH3_SN_LEGACY`, median continuum flux, median/minimum uncertainty, negative-flux fraction, and valid-channel counts. These columns are the first place to look when an achieved-S/N value appears pathological.

Script 2 also logs the requested and usable observed-frame wavelength interval for each arm. `SN_WINDOW_COVERAGE_WARNING` means the configured window is truncated by Script-1 `GOODWAVE`; the pipeline does not search for or adopt a substitute window automatically.

---

## `BL_RH3_SN_comparison.png`

**Script:** 02  
**Class:** major QC + possible paper figure

### Purpose

Directly compares achieved S/N in the *same physical PowerBins* between the BL and RH3 arms.

### Panels

- left: BL S/N;
- right: RH3 S/N assigned onto the BL master bin map.

### Color normalization

Both panels share exactly the same `vmin` and `vmax`.

Default upper limit:

$$
v_{\max}=P_{95}\left(\{S/N_{\mathrm{BL}}\}\cup\{S/N_{\mathrm{RH3}}\}\right).
$$

Values above the 95th percentile saturate visually but remain numerically preserved in the table and FITS products.

### Interpretation

Equivalent colors mean equivalent S/N. This makes it immediately obvious whether the red/RH3 spectra have enough information in the BL-defined apertures and prevents independent autoscaling from hiding cross-arm differences. Undefined S/N bins are gray in-place rather than disappearing, so the viewer can distinguish an S/N pathology from a missing spatial bin.

---

## `RH3_SN_window_scan.png` / `BL_SN_window_scan.png`

**Script:** development utility `scripts/scan_script02_sn_windows.py`  
**Class:** optional integration-test QC; not part of the production pipeline

### Purpose

Scans fixed-width **observed-frame** continuum windows through the already extracted Script-2 bin spectra. This is useful when non-production test data use a grating whose useful wavelength range is a poor match to the production CaT S/N window.

### Quantities

- fraction of evaluable bins with positive median continuum;
- fraction of evaluable bins with extreme signed S/N;
- median robust S/N among positive-continuum bins.

The companion `*_SN_window_scan.ecsv` additionally records window bounds, native channel count, number of evaluable/positive bins, and the 10th-percentile positive-bin S/N.

### Interpretation

Use this plot only to diagnose whether pathologies are localized in wavelength. It does not define the production RH3 science window, never edits a target config, and should not be used to optimize away real sky/telluric or reduction problems.

---

## `RH3_SN_pathology_bin####_flux.png` / `_uncertainty.png` / `_ratio.png`

**Script:** `inspect_script02_sn.py`
**Class:** targeted QC/debugging

### Purpose

These files are generated only when the standalone Script-2 S/N inspection utility is run. They display the coadded flux, formal uncertainty, and sample-by-sample $F/\sigma$ values for the most negative legacy RH3 S/N bin in the configured window.

### Interpretation

If the continuum itself is negative while the uncertainty is ordinary, the correct conclusion is that the bin has no positive continuum detection in that window. If the legacy median-of-ratios is enormous while the robust ratio-of-medians is modest, unusually small formal-uncertainty samples were dominating the old estimator. If both are enormous, inspect the formal variance and wavelength-window placement before trusting the S/N diagnostic. The utility never edits the science products.

### Numerical companion

`products/RH3_SN_diagnostic_report.ecsv` contains the old table value, robust achieved S/N, signed robust S/N, legacy median-of-ratios, median flux, median/minimum/5th-percentile uncertainty, negative-flux fraction, and valid-channel count for every bin.

---

## `bin_area_map.png`

**Script:** 02  
**Class:** QC

### Purpose

Maps the physical area of every BL master PowerBin in square arcseconds.

### Interpretation

Large outer bins are expected as surface brightness falls, but their size directly sets the spatial resolution of later kinematic/population maps and the scale over which intra-bin velocity shear must eventually be evaluated.

### Warning signs

A few bins becoming dramatically larger than neighboring bins can indicate an overly permissive analysis aperture, poor S/N, or masked geometry that deserves inspection.

### Numerical companions

Area and BL/RH3 member counts are saved per bin in `master_bin_table.ecsv`, while exact native member pixels are saved in `master_bin_membership.npz`.

---

# Script 3 — RH3 profile likelihood and preliminary kinematics

Script 3 is intentionally verbose diagnostically because its independent 3-D spectral likelihood cubes are the information product used by every later decomposition stage. A clean-looking local minimum is not, by itself, a disk assignment.

## `RH3_single_component_velocity.png`

**Script:** 03  
**Helper:** `plot_bin_value_map()`  
**Class:** QC + science comparison

### Purpose

Maps the single-LOSVD stellar velocity fit in every BL-defined master PowerBin, using the same XSL basis, LSF, fixed good pixels, and continuum convention as the two-component likelihood calculation.

### Healthy result

A coherent large-scale rotation pattern is present where the data contain useful stellar information. Individual outer bins may be noisy; this map is not spatially regularized.

### Warning signs

Large isolated velocity excursions, widespread bound-sensitive fits, or structure inconsistent with the raw spectra should trigger inspection of wavelength zero points, the LSF, and the fixed fit mask before Script 4.

---

## `RH3_single_component_sigma.png`

**Script:** 03  
**Helper:** `plot_bin_value_map()`  
**Class:** critical geometry input + possible publication figure

### Purpose

Maps the independent one-component stellar dispersion measurements. For real RH3 data this is the higher-resolution starting point for locating the $2\sigma$ morphology used by the later radial-domain analysis.

### Important limitation

The present Script-3 run uses a development diagonal noise model; formal likelihood widths/errors are not publication-final until residual variance/covariance is calibrated. The map itself is still useful for first-pass morphology and integration testing.

---

## `RH3_local_best_VA.png` / `RH3_local_best_VB.png`

**Script:** 03  
**Helper:** `plot_bin_value_map()`  
**Class:** analysis diagnostic only

### Purpose

Displays the exact velocity-grid coordinates of the **independent minimum of each individual PowerBin's 3-D cube**.

### Interpretation

These maps are deliberately allowed to look ragged, swap components, or jump between secondary minima. Script 3 has imposed no spatial coherence. Their scientific role is to show what each spectrum independently prefers before Script 4 asks whether a globally coherent two-disk trajectory can pass through the likelihood surfaces.

Do not interpret the labels A/B here as final physical disk identities.

---

## `RH3_local_best_fA.png`

**Script:** 03  
**Helper:** `plot_bin_value_map()`  
**Class:** analysis diagnostic only

### Quantity

The explicit fraction-grid coordinate of the independent local cube minimum.

### Definition

Every XSL SSP is independently normalized to unit mean stellar flux density over `RH3_FIT_REST_RANGE_ANGSTROM` before the basis is duplicated into components A and B. Therefore the grid variable

$$
f_{A,\mathrm{RH3}}
=
\frac{\sum_jw_{A,j}}{\sum_jw_{A,j}+\sum_jw_{B,j}}
$$

is the fraction of fitted **stellar-template light** assigned to component A over the RH3 fitting passband. It is not a mass fraction, and the additive polynomial is excluded from the component-light definition.

### Warning signs

A large fraction of local minima at 0.1 or 0.9 can mean the data frequently prefer an effectively one-component solution, the fraction grid is truncating a profile, or the decomposition is weak. Component swaps naturally map $f_A\leftrightarrow1-f_A$ and are not resolved until Script 4.

---

## `RH3_local_best_sigmaA.png` / `RH3_local_best_sigmaB.png`

**Script:** 03  
**Helper:** `plot_bin_value_map()`  
**Class:** QC + analysis diagnostic

### Purpose

Shows the nuisance dispersions profiled at each bin's independent local 3-D minimum.

### Warning signs

Values repeatedly close to the configured 5 or 250 km/s bounds indicate that the nuisance profile may be truncated or that a two-component state is using an extreme width to compensate for an inappropriate velocity/fraction hypothesis.

Use `RH3_sigma_boundary_fraction.png` to see whether bound contact is isolated or widespread through the full cube.

---

## `RH3_local_one_vs_two_delta_chi2.png`

**Script:** 03  
**Helper:** `plot_bin_value_map()`  
**Class:** QC; **not** a formal significance map

### Quantity

$$
T_i
=
\chi^2_{1C,i}
-
\min_{V_A,V_B,f_A}\chi^2_{2C,i}.
$$

Both terms use the same fixed fitted pixels and formal diagonal uncertainty vector.

### Interpretation

Positive values mean that the best local two-component grid state fits better than the one-component control. Do **not** label this as a p-value or a sigma detection. The mapping from $T_i$ to reliable two-component recovery must be calibrated with injection/recovery mocks.

---

## `RH3_grid_failure_fraction.png`

**Script:** 03  
**Helper:** `plot_bin_value_map()`  
**Class:** critical numerical QC

### Quantity

Fraction of the configured $V_A\times V_B\times f_A$ states in each PowerBin that did not produce a valid exact pPXF state.

A healthy production calculation should normally be zero or extremely small. Coherent regions of failures indicate a numerical, masking, template-coverage, or bound problem that should be resolved before Script 4.

---

## `RH3_sigma_boundary_fraction.png`

**Script:** 03  
**Helper:** `plot_bin_value_map()`  
**Class:** numerical + scientific QC

### Quantity

Fraction of exact grid states whose profiled $\sigma_A$ or $\sigma_B$ lies within `RH3_SIGMA_BOUNDARY_WARNING_KMS` of the configured 5--250 km/s nuisance bounds.

### Interpretation

Some poor high-$\chi^2$ hypotheses can naturally profile to an extreme sigma and are not automatically a problem. A large boundary fraction **near the supported likelihood basin**, however, means the nuisance parameter range may be affecting the likelihood geometry and must be investigated.

---

## `RH3_likelihood_bins/bin_XXXX.png`

**Script:** 03  
**Helper:** `plot_rh3_likelihood_bin()`  
**Class:** per-bin analysis diagnostic

### Purpose

Provides a compact inspection of the independent 3-D profile-likelihood cube and its corresponding spectrum.

### Implemented content

The figure profiles the 3-D total-$\chi^2$ cube over $f_A$ to show the $V_A/V_B$ likelihood structure, displays the best fraction across velocity space, marks the independent local cube minimum, and compares the one-component and local-best two-component spectral models/residuals on the same fixed good pixels.

### Healthy result

High-information bins can show a localized or clearly structured minimum/basin. Lower-information bins may legitimately show broad valleys or near-symmetries. That breadth is information and should not be "cleaned" by propagating neighboring starting guesses.

### Warning signs

- minimum on a velocity/fraction-grid edge;
- widespread failed cells;
- local best fit visibly poor despite a numerical minimum;
- a nearly featureless surface in a bin being described as an independent two-component detection;
- a sharp surface whose width is trusted before the noise/covariance model is calibrated.

---

## Numerical companions for Script 3

### `products/RH3_likelihood_cubes.npz`

The fundamental per-state scalar product. It contains the full total-$\chi^2$, reduced-$\chi^2$, $\sigma_A$, $\sigma_B$, fit-status, and sigma-boundary arrays plus all three grid axes.

### `products/RH3_log_spectra_and_local_best_fits.npz`

Stores the fixed log-rebinned spectra/noise/masks together with the complete one-component model and the independent local-best two-component model/template weights for every bin.

### `products/RH3_local_likelihood_summary.ecsv`

Human-readable per-bin table containing the one-component control, local cube minimum, achieved-vs-requested $f_A$ check, local $\Delta\chi^2$, failed-state counts, and sigma-boundary-state counts.

### `products/XSL_RH3_templates.npz`

Records the actual LSF-matched, log-rebinned, passband-normalized XSL basis used by the run, including the template normalization factors and target/template/convolution FWHM arrays.

### Restart checkpoints

During an incomplete run, `checkpoints/bin_XXXX.npz` files are the atomic restart units used by `--resume`. They are deleted only after a complete successful consolidation and manifest write when `SCRIPT03_DELETE_CHECKPOINTS_ON_SUCCESS=True`.

---

## `RH3_grid_resolution_comparison.png`

**Script:** later Script-3 validation utility/run  
**Class:** critical method-validation diagnostic; **not yet implemented by the baseline driver**

### Purpose

Compare representative coarse and finer profile-likelihood grids to demonstrate that the baseline 17×17×9 grid does not materially alter the selected basin topology. If the topology or global solution moves materially, the production grid must be refined.

---

# Script 4 — radial geometry and global two-disk model

## `xooksuut_radial_model.png`

**Script:** 04  
**Class:** major methods/publication figure

### Purpose

Displays the exact non-parametric radial-ring geometry used by the global two-disk model.

### Must show

- galaxy center;
- signed $PA_{\mathrm{kin}}$ axis;
- projected ring ellipses;
- $R_{\mathrm{start}}$;
- each ring center;
- $R_{2\sigma,+}$ and $R_{2\sigma,-}$;
- adopted $R_{2\sigma,\mathrm{max}}$;
- snapped $R_{\mathrm{final}}$;
- final annulus $R_{\mathrm{final}}\pm\delta$;
- PSF scale;
- PowerBin footprint.

### Mathematical definitions

$$
R_{\mathrm{start}}=\mathrm{FWHM}_{\mathrm{PSF}},
$$

$$
\delta=0.5\,\mathrm{ring\_space},
$$

and the snapped outer-radius relation documented in `METHODS.md`.

### Healthy result

The radial nodes resolve the galaxy at approximately the delivered PSF scale and the final $2\sigma$-limited ring encloses the relevant two-component region without an obvious excessive gap.

---

## `disk_label_convention.png`

**Script:** 04  
**Class:** critical bookkeeping QC + methods

### Purpose

Makes Disk A/B identity reproducible and visually auditable.

### Shows

- oriented signed major-axis line;
- negative and positive direction labels;
- reference outer spatial bin/region;
- preliminary outer single-component velocity sign;
- branch assigned Disk A;
- branch assigned Disk B.

### Warning

If component identity is ambiguous, the plot must display the `AMBIGUOUS_DISK_LABEL` state rather than drawing a confident assignment.

---

## `ring_coverage.png`

**Script:** 04  
**Class:** QC

### Purpose

Quantifies how strongly each non-parametric ring is actually constrained by usable spatial data.

### Possible panels

- number of PowerBins contributing by radius;
- azimuthal coverage fraction;
- weighted coverage metric.

### Warning signs

A fitted ring controlled by only one tiny spatial region or strongly one-sided azimuthal sampling.

### Action

Flag `POOR_RING_COVERAGE`; consider coarser ring spacing or report the affected ring as weakly constrained rather than smoothing it silently.

---

## `global_model_VA.png` / `global_model_VB.png`

**Script:** 04  
**Class:** science diagnostic

Maps the bin-integrated model-predicted Disk A and B velocities at the best global $\Theta$.

---

## `global_model_residuals_A.png` / `global_model_residuals_B.png`

**Script:** 04  
**Class:** critical QC

### Quantity

Difference between the representative selected spectral solution and the global model prediction, in km/s, using a precisely documented definition of the selected per-bin velocity.

### Warning signs

Coherent residual structures, bars/noncircular motions, or one side of the galaxy consistently offset.

---

## `global_normalized_residuals.png`

**Script:** 04  
**Class:** QC

Normalizes residuals by an appropriate velocity uncertainty/likelihood scale where available. This highlights statistically meaningful deviations that may be visually small in km/s.

---

## `best_theta_rotationcurves.png`

**Script:** 04  
**Class:** major science/publication figure

### Purpose

Shows the final non-parametric Disk A and Disk B rotation curves.

### Display

- ring-node velocities as explicit markers;
- linear interpolation between nodes;
- optional uncertainty bands once Script 7 is complete;
- $R_{\mathrm{start}}$, $R_{2\sigma}$, and radial-limit markers.

### Important

Do not draw an analytic smooth fit that was not used by the inference.

---

## `Q_optimization_history.png`

**Script:** 04  
**Class:** optimizer QC

### Purpose

Shows the history of

$$
Q(\Theta)=\sum_iq_i(\Theta)
$$

across optimizer evaluations/iterations.

### Recommended curves

- raw $Q$ for every evaluated state;
- best-so-far $Q$.

### Healthy result

Best-so-far curve approaches a plateau and later evaluations fail to find substantial improvements.

### Warning signs

Continual late improvements, repeated large jumps between basins, or termination immediately after a new best solution.

---

## `geometry_initial_vs_final.png`

**Script:** 04  
**Class:** QC

Compares initial and fitted values of center, PA, inclination, and systemic velocity, including external uncertainty ranges/bounds.

Large movement to a hard prior/boundary should be obvious.

---

## `intra_bin_velocity_shear.png`

**Script:** 04  
**Class:** critical QC

### Quantity

$$
\sigma_{\mathrm{shear},i}^2=\frac{\sum_p I_p(V_p-\bar V_i)^2}{\sum_pI_p}.
$$

### Purpose

Identifies PowerBins whose internal velocity gradient can broaden the integrated spectrum and mimic intrinsic stellar dispersion.

### Warning signs

Shear comparable to the precision of $\sigma_A$ or $\sigma_B$.

### Action

Flag `HIGH_BIN_SHEAR` and consider explicit forward-model treatment or rebinning validation.

---

## `centroid_vs_bin_integrated_velocity.png`

**Script:** 04  
**Class:** approximation-validation QC

### Quantity

$$
\Delta V_i=V_{i,\mathrm{bin-integrated}}-V_{i,\mathrm{centroid}}.
$$

### Purpose

Directly demonstrates whether the simpler centroid approximation would have been safe.

### Healthy result

Differences are much smaller than the kinematic uncertainties for nearly all bins.

---

# Radial-extent comparison diagnostics

## `radial_extent_bin_usage.png`

**Script:** 04  
**Class:** methods/QC

Shows which PowerBins enter the $2\sigma$-limited model and which additional bins enter the full-aperture model.

---

## `radial_extent_rotationcurve_comparison.png`

**Script:** 04 / finalized 07  
**Class:** major robustness/publication candidate

### Purpose

Compares

$$
V_A^{2\sigma}(R),\quad V_A^{\mathrm{full}}(R)
$$

and the same for Disk B.

### Critical comparison region

The shared inner radial domain. The outer extension naturally exists only for the full-aperture model.

### Interpretation

- close overlap: outer low-minority-fraction data are consistent and may safely help constrain $\Theta$;
- large inner shifts: outer bins are influencing the decomposition enough that the $2\sigma$-limited solution should likely remain fiducial.

---

## `radial_extent_geometry_comparison.png`

**Script:** 04  
**Class:** robustness QC

Compares $x_0,y_0,PA,i,V_{\mathrm{sys}}$ between the two radial-extent fits, with uncertainties when available.

---

## `radial_extent_velocity_difference.png`

**Script:** 04  
**Class:** robustness QC

PowerBin map of the difference between velocities predicted by the $2\sigma$-limited and full-aperture best global solutions over their common domain.

---

# Script 5 — exact RH3 extraction, basins, grid validation, one/two-component reliability

## `RH3_selected_basin.png`

**Script:** 05  
**Class:** major methods diagnostic

### Purpose

Explains which 3-D likelihood basin is associated with the global XookSuut-style anchor.

### Recommended panels

1. profiled $\Delta\chi^2(V_A,V_B)$;
2. basin membership projected into velocity space;
3. best $f_A$ surface;
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

For basin cell $c$,

$$
p_c\propto\exp(-\Delta\chi_c^2/2).
$$

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

$$
\Delta\chi^2_{\mathrm{cube},i}=\chi_i^2(\mathrm{selected})-\min_{V_A,V_B,f_A}\chi_i^2.
$$

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

Maps the selected solution's $\Delta\chi^2$ relative to the minimum of its own selected basin. Distinguishes “wrong basin relative to global cube” from “far up the wall of the chosen basin.”

---

## `RH3_grid_edge_distance.png`

**Script:** 05  
**Class:** critical QC

### Quantity

Minimum number of grid cells between the selected discrete state and any $V_A$, $V_B$, or $f_A$ boundary.

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

- $V_A^{\mathrm{exact}}-V_A^{\mathrm{grid}}$;
- $V_B^{\mathrm{exact}}-V_B^{\mathrm{grid}}$;
- $\sigma_A^{\mathrm{exact}}-\sigma_A^{\mathrm{grid}}$;
- $\sigma_B^{\mathrm{exact}}-\sigma_B^{\mathrm{grid}}$;
- $f_A^{\mathrm{fine}}-f_A^{\mathrm{coarse}}$.

### Interpretation

Small differences validate that the cube is adequate for likelihood topology while the final exact refit removes quantization.

---

## `RH3_fA_refinement.png`

**Script:** 05  
**Class:** per-bin QC

Plots $\chi^2(f_{A,\mathrm{RH3}})$ at the exact global velocities across the fine 0.01-step local fraction scan.

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

- single-component false-positive distribution of $T$;
- two-component recovery fraction versus $\Delta V$;
- recovery versus $f_A$;
- recovery versus S/N;
- interpolation region in $(S/N,\Delta V,f_A,\sigma_A,\sigma_B)$ space.

### Important

Do not label the empirical recovery fraction as the posterior probability that an observed spectrum “is two-component.”

---

## `two_component_false_positive_map.png`

**Script:** 05  
**Class:** science-quality map

Maps

$$
p_{\mathrm{false},i}=P(T\ge T_i\mid\mathrm{true\ one\ component})
$$

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

Gas-template weights are excluded from $f_{A,\mathrm{BL}}$.

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

Shows the profile joint statistic or BL statistic versus explicit stellar $f_{A,\mathrm{BL}}$ for the anchored and/or joint-profile solution.

### Coarse/fine display

Plot coarse 0.1-step evaluations and fine 0.01-step refinement distinctly so the viewer can see that the final fraction is not a coarse-grid artifact.

---

## `BL_fA_marginalized_distribution.png`

**Script:** 06  
**Class:** major inference diagnostic

### Purpose

Shows the RH3-likelihood-propagated support for $f_{A,\mathrm{BL}}$ after integrating/reweighting over the sampled RH3-supported kinematic family.

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

Profiled $(V_A,V_B)$ RH3 contours with markers for:

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
- successive lower-$J$ neighboring cells;
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

- $\log Age_A$, $\log Age_B$;
- $[M/H]_A$, $[M/H]_B$;
- $f_{A,\mathrm{BL}}$;
- $\Delta\log Age$;
- $\Delta[M/H]$.

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

$$
\Delta V=|V_A-V_B|.
$$

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

Same structure for $[M/H]$.

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

$$
\Delta\log Age=\log Age_A-\log Age_B.
$$

### Interpretation

Positive and negative values identify which disk is older. The sign convention must be written directly in the figure caption/metadata.

---

## `disk_contrast_delta_metallicity.png`

**Quantity**

$$
\Delta[M/H]=[M/H]_A-[M/H]_B.
$$

This is a direct physical contrast between the two stellar disks.

---

## `method_delta_age_A.png` / `method_delta_age_B.png`

**Script:** 07  
**Class:** robustness maps

### Purpose

These are *not* disk A-minus-B contrasts. They compare inference methods:

$$
\Delta Age_{A,\mathrm{method}}=Age_{A,\mathrm{likelihood-propagated}}-Age_{A,\mathrm{anchored}}.
$$

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

Cumulative $N_{\mathrm{MC}}$.

### y-series

Representative or aggregate 16th/50th/84th percentile summaries for the principal science quantities.

### Convergence

Use a documented batch/interval stability rule rather than selecting $N_{\mathrm{MC}}$ solely because it is a round number.

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
7. Are the measured RH3 $2\sigma$ peak radii and XookSuut-style radial nodes sensible?
8. Are non-parametric ring velocities constrained by adequate spatial coverage?
9. Does the global disk model fit the selected spectral likelihood basins without strongly overriding individual spectra?
10. Are centroid-only and bin-integrated velocities consistent, or is intra-bin shear important?
11. Do the $2\sigma$-limited and full-aperture global solutions agree over their common radial range?
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
