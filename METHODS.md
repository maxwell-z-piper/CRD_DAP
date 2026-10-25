# CRD_DAP Scientific and Statistical Methods

## Purpose of this document

This file is the authoritative design document for the CRD_DAP analysis pipeline. It records not only **what** each stage does, but **why** the analysis is structured this way, what assumptions are being made, which quantities are profile-likelihood quantities rather than formal posteriors, how uncertainty is propagated, and which numerical choices must still be validated on simulations or real data.

The intent is that a future reader should be able to reconstruct the scientific logic without relying on memory or old chat logs. If the implementation changes, this document should be updated at the same time.

---

# 1. Scientific goal

The project targets galaxies containing two counter-rotating stellar disks. The primary goals are to recover spatially coherent maps of the two stellar LOSVDs and then measure the stellar populations of the two components separately.

The two KCWI/KCRM arms play complementary roles:

- **RH3** provides high spectral resolution around the Ca II triplet region and is the primary source of two-component stellar kinematics.
- **BL** provides the broad blue wavelength information needed for age, metallicity, and other stellar-population measurements.

The core methodological idea is to avoid reducing each spatial bin to one independently selected two-component velocity solution. Instead, the RH3 spectrum in each spatial bin is compressed into a profile-likelihood surface

$$
\chi_i^2(V_A,V_B,f_{A,\mathrm{RH3}}),
$$

and a global two-disk rotation model selects a spatially coherent trajectory through all of those surfaces. The family of RH3-supported kinematic states is then propagated into the BL stellar-population decomposition.

This preserves more spectral information than constructing an independent best-fit velocity map first and smoothing or fitting that map afterward.

---

# 2. Pipeline overview

The seven science-driver scripts are intended to be:

1. `01_prepare_and_register_cubes.py`
2. `02_make_master_BL_bins.py`
3. `03_build_RH3_likelihood_cubes.py`
4. `04_fit_global_two_disk_model.py`
5. `05_extract_global_RH3_solution.py`
6. `06_fit_BL_two_component_populations.py`
7. `07_uncertainties_and_final_maps.py`

The driver scripts should stay compact. Detailed calculations belong in the `crd_utils` package.

A schematic flow is

$$
\mathrm{BL+RH3\ cubes}
\rightarrow
\mathrm{prepared/registered\ cubes}
\rightarrow
\mathrm{BL\ master\ PowerBins}
\rightarrow
\chi_i^2(V_A,V_B,f_A)
\rightarrow
\mathrm{global\ two-disk\ model}
\rightarrow
\mathrm{exact\ RH3\ solution}
\rightarrow
\mathrm{BL\ populations}
\rightarrow
\mathrm{uncertainties/final\ maps}.
$$

---

# 3. Reproducibility architecture

## 3.1 Target configuration

All target-specific values and analysis hyperparameters are stored in a target configuration derived from `config/target_config_template.py`. Science functions should not contain hidden target-specific constants.

Required paths include:

- stacked BL cube;
- stacked RH3 cube;
- BL master arc;
- RH3 master arc;
- PyMorph VAC or local target extract;
- XSL SSP library.

The master arcs are required because the adopted baseline measures the instrumental line-spread function empirically from calibration lines.

Each run snapshots the exact configuration into its run directory.

## 3.2 Logging

Every run has a master `pipeline.log` and one log per science stage. Messages are written to both the terminal and file logs.

Long-running stages must log:

- start time;
- input files;
- relevant configuration parameters;
- progress;
- per-bin or per-batch timing;
- estimated remaining time where meaningful;
- warnings and quality flags;
- completion summary.

A long pPXF calculation should therefore be able to run unattended.

## 3.3 Diagnostic documentation

Every scientifically meaningful plot is documented in `DIAGNOSTICS.md` with:

- filename;
- producing script/helper;
- scientific purpose;
- mathematical definition;
- axes and units;
- plotting normalization;
- overlays;
- healthy interpretation;
- warning signs;
- recommended response to a warning;
- associated numerical products and log fields.

A new diagnostic is not considered finished until its `DIAGNOSTICS.md` entry exists.

---

# 4. Script 1 — prepare and register cubes

## 4.1 Inputs

The production science input for each arm is the four-file KcwiKit post-DRP stack:

- `*_icubes.fits` — stacked science flux;
- `*_vcubes.fits` — stacked variance;
- `*_mcubes.fits` — final binary stack mask;
- `*_ecubes.fits` — effective exposure time in seconds.

Script 1 therefore requires the four BL stack files and the four RH3 stack files, plus:

- BL master arc;
- RH3 master arc;
- configuration metadata.

The science cubes are assumed to have already passed through the KCWI DRP and KcwiKit stacking. Script 1 is science preparation, not detector-level reduction or cube stacking. A native DRP multi-extension cube remains supported for legacy/testing use, but is not the production stacked-data path.

KcwiKit's current stacking implementation uses the original PyDRP `FLAGS` cube to reject invalid contributing samples. The final stacked `mcubes` product is then binary: zero means valid contributing exposure exists and one means no valid contribution. The original bit-level DRP flag values therefore cannot be reconstructed from the final stack and CRD_DAP does not attempt to invent them.

## 4.2 Internal cube representation

The pipeline should expose, for each arm:

$$
F(x,y,\lambda),
$$

variance and 1-sigma uncertainty,

$$
\mathrm{Var}_F(x,y,\lambda),\qquad \sigma_F(x,y,\lambda)=\sqrt{\mathrm{Var}_F(x,y,\lambda)},
$$

effective exposure time,

$$
E(x,y,\lambda),
$$

quality/mask information, wavelength coordinates, spatial WCS, spectral-resolution metadata, and relevant reduction/stacking metadata.

Downstream scripts should not need to know KcwiKit filenames, FITS-axis ordering, or the distinction between variance and standard deviation. Script 1 validates that the four companion cubes in each arm have identical array shapes and identical spatial/spectral WCS before combining them into the internal data model.

The KcwiKit stacks used during development are stored as 64-bit floating-point arrays and can be several GB across both arms. The production default is therefore to load/write prepared large cubes as float32 for practical memory use, while allowing extracted spectra to be promoted to float64 for pPXF. This choice is configurable rather than hidden.

## 4.3 Hard bad-pixel / bad-spaxel treatment

A KcwiKit spectral sample is unusable when, for example:

- the final stacked mask is non-zero;
- the effective exposure is zero or non-finite;
- the flux is non-finite;
- the variance is non-finite or non-positive;
- it lies outside the instrumentally valid wavelength range.

The requested KcwiKit output grid can intentionally be larger than the illuminated IFU footprint. Zero-exposure padding must not count as a population of "bad spaxels" when deciding whether an entire wavelength channel is globally unusable. Script 1 therefore computes wavelength-level bad-sample fractions only over the spatial footprint with real wavelength coverage, while the zero-exposure padding remains excluded from every downstream science mask.

A spatial spaxel may be excluded if the fraction of usable wavelength samples falls below the configurable threshold. **Low S/N by itself is not a reason to reject a spaxel**; Script 2 uses PowerBin to combine low-S/N spatial information.

## 4.4 Wavelength masks

Three concepts must remain distinct:

1. **instrument-good wavelengths** — valid according to calibration/data quality;
2. **RH3 fit-good wavelengths** — wavelengths used for RH3 stellar fitting;
3. **BL fit-good wavelengths** — wavelengths used for BL stellar/gas fitting.

Scientifically valid wavelengths should not be permanently discarded in Script 1 merely because a later fit masks them.

## 4.5 Wavelength conventions

Air/vacuum wavelength conventions and the adopted heliocentric/barycentric frame must be known explicitly. The pipeline must hard-fail rather than proceed with science and template wavelengths in inconsistent conventions.

A convention mismatch can mimic a velocity zero-point error and would contaminate both the RH3 likelihood surfaces and the inferred systemic velocity.

## 4.6 Collapsed-continuum images

For each arm, construct a robust 2-D continuum image by collapsing clean continuum wavelengths:

$$
I(x,y)=\rm{median}_{\lambda\in\mathrm{good}}F(x,y,\lambda)
$$

or an appropriately weighted equivalent.

These images are used for:

- visual QC;
- approximate center estimation;
- BL/RH3 spatial-registration checks;
- identifying pathological spatial regions.

They should be called `*_collapsed_continuum.png`, not “white-light” images, to avoid implying a standard broadband filter response.

## 4.7 Initial galaxy center

The primary KCWI-centered initializer should come from the collapsed continuum image, not from manual eyeballing.

Useful center estimates include:

- brightest positive continuum spaxel after light smoothing;
- continuum centroid/2-D photometric center;
- transformed PyMorph center if an accurate SDSS-to-KCWI coordinate transform is available.

These estimates are compared. The final global kinematic center remains a free parameter in Script 4 within sensible bounds.

## 4.8 BL/RH3 registration

BL and RH3 must be placed in a common sky-coordinate system. Identical array indices should not be assumed to represent identical physical positions unless WCS verification proves that they do.

When the two arms share a sufficiently wide instrument-good observed-frame wavelength interval, the registration diagnostic should collapse the *same wavelength interval* in both cubes before morphology cross-correlation. This minimizes wavelength-dependent stellar-population/color structure as a false source of apparent astrometric shift. The minimum overlap width and number of channels are configurable. If no useful overlap exists, the normal full-arm continuum images may still be compared.

A numerical morphology cross-correlation shift is trusted only when both registration images contain sufficient spatial contrast. Script 1 records a robust contrast statistic for each image and returns the cross-correlation as **inconclusive** rather than forcing a shift when either image is nearly featureless. In that case, the independent sky-coordinate peak/centroid comparison remains the relevant astrometric QC and no spurious residual shift is adopted.

Because the cubes have already been placed on celestial WCS, this step measures only a **local residual** registration. The morphology cross-correlation is therefore restricted to a configurable maximum residual-search radius around zero shift. A correlation peak that runs to this boundary is classified as inconclusive rather than being reported as a large blind translation; a genuinely large astrometric disagreement should be repaired in the stacking/WCS stage.

The science cubes themselves are not resampled in Script 1. The output must allow Script 2 to apply one physical BL-defined PowerBin membership to the RH3 cube using saved sky/tangent-plane coordinates.

## 4.9 PSF characterization

Measure or adopt separate effective delivered PSFs:

$$
\mathrm{PSF}_{\mathrm{BL}},\qquad \mathrm{PSF}_{\mathrm{RH3}}.
$$

The arms use the same slicer geometry but may not have identical effective PSFs. This matters because using the same geometric PowerBin does not guarantee identical spatial weighting of stellar light in both arms.

The RH3 PSF will also define the default XookSuut-style radial node spacing.

## 4.10 Empirical LSF from master arcs

The primary LSF comes from unresolved lines in the required master-arc products.

The expected grating in each pipeline data stream is explicit configuration metadata: `BL_EXPECTED_GRATING` for the blue stream and `RH3_EXPECTED_GRATING` for the red stream. Their production defaults are `BL` and `RH3`, respectively. This configurability is intended to support controlled Script-1 validation on other KCWI/KCRM setups without weakening calibration safeguards. For every run, the science cube and master arc must both report the configured expected grating and must still match one another in camera, slicer/IFU, detector binning, and central wavelength. A noncanonical grating is logged prominently and must not be interpreted as production BL/RH3 science in later grating-specific stages unless those stages have been explicitly adapted.

Aim to estimate

$$
\mathrm{FWHM}_{\mathrm{inst}}(\lambda,s),
$$

where $s$ is slice/spatial position if the data support measuring a meaningful spatial dependence. If spatial variation is negligible, collapse to a wavelength-only model

$$
\mathrm{FWHM}_{\mathrm{inst}}(\lambda).
$$

The measured LSF is compared against the XSL template LSF. The required Gaussian template broadening is conceptually

$$
\sigma_{\mathrm{conv}}^2 = \sigma_{\mathrm{data}}^2-\sigma_{\mathrm{template}}^2.
$$

If the right-hand side is negative, the templates are broader than the data at that wavelength and the pipeline must report this explicitly rather than taking an absolute value or silently forcing a convolution.

The instrument-good interval from `WAVGOOD0/1` and the **empirical accepted-line support** of the LSF model are distinct quantities. The wavelength polynomial may be fit using line measurements that occupy only a subset of the instrument-good range. Script 1 therefore stores both intervals. By default, evaluating the saved LSF model outside the bluest-to-reddest accepted line centers returns `NaN` rather than silently extrapolating the polynomial. A configurable edge-gap threshold raises `LSF_EMPIRICAL_COVERAGE_GAP` when a substantial instrument-good region lacks direct arc-line support. Later stellar fitting must mask such wavelengths or supply an independently validated LSF treatment.

Individual arc-line widths can show substantially more scatter than the coherent resolution variation across the field. Script 1 therefore preserves both: the full line-by-line residual diagnostic and a spatially averaged summary of slice and within-slice position-bin medians. The latter is the more relevant test of whether one wavelength-only LSF is adequate for the science.

## 4.11 Noise and covariance characterization

The absolute $\chi^2$ scale matters because later stages use

$$
L_{\mathrm{rel}}\propto e^{-\Delta\chi^2/2}.
$$

If uncertainties are underestimated by a factor $a$, $\chi^2$ differences become artificially too large by approximately $a^{-2}$, producing likelihood basins that are too narrow.

The adopted spectral-fit likelihood uses a wavelength-dependent, scalar-calibrated
uncertainty vector and a **diagonal covariance between wavelength samples**. It
therefore treats errors at different wavelengths as uncorrelated for the fit;
this is an analysis assumption, not a claim that the reduced spectra are
physically free of spectral correlation. No spectral lag coefficients or dense
covariance matrices are fitted. Script 1's high-pass correlation measurement is
a preliminary diagnostic and does not change the likelihood or set its noise
scale. Scripts 2 and 3 calibrate the uncertainty amplitude as described below.
Any profile-likelihood interval is conditional on this diagonal approximation.

---

# 5. Script 2 — BL-defined master PowerBins

## 5.1 Stage input and provenance

Script 2 consumes the immutable products written by one completed Script-1 run:

- `prepared_BL.fits`;
- `prepared_RH3.fits`;
- `BL_spatial_coordinates.npz`;
- `RH3_spatial_coordinates.npz`;
- the preliminary Script-1 spectral-correlation products when available.

The selected Script-1 run directory is written into the Script-2 FITS headers and JSON manifest. Script 2 must not reconstruct the original DRP/KcwiKit masks; it uses the final `GOODMASK`, `GOODSPAX`, and `GOODWAVE` state saved by Script 1.

## 5.2 BL alone defines the spatial tessellation

PowerBin is run **once** on BL. The baseline target is

$$
(S/N)_{\mathrm{BL}}\approx40.
$$

The configured target is interpreted as an approximate continuum S/N per spectral pixel in a user-selected BL rest-frame fitting window. For each BL spatial sample $p$, Script 2 calculates a robust continuum signal $S_p$ and formal noise $N_p$ from the median flux density and median variance over usable wavelengths in that window.

For a candidate bin with member set $\mathcal B_i$, the baseline diagonal-noise capacity is

$$
C_i=\left[\frac{\sum_{p\in\mathcal B_i}S_p}{\sqrt{\sum_{p\in\mathcal B_i}N_p^2}}\right]^2.
$$

PowerBin therefore receives the callable capacity $C_i$ and the target capacity

$$
C_{\mathrm{target}}=(S/N)_{\mathrm{BL,target}}^2.
$$

The callable form is retained even in the diagonal-noise baseline because PowerBin permits non-additive capacities. A future empirically calibrated KcwiKit **spatial** covariance correction can therefore be inserted without changing the tessellation architecture. Script 1 measured spectral correlation along wavelength; that result is not by itself a valid spatial covariance law, so Script 2 does not invent one.

PowerBin is mandatory for the production pipeline. Script 2 must never silently fall back to legacy VorBin.

## 5.3 Useful stellar-body aperture

The Script-1 hard-good footprint identifies where data exist, but it includes blank sky around the galaxy. Giving every exposed sky pixel to PowerBin would allow very large low-information outer bins to form.

Script 2 therefore defines a separate **analysis aperture** before PowerBin. The baseline automatic mode uses a smoothed BL continuum-detection proxy, selects the connected component nearest the independently adopted galaxy center, optionally dilates that component, and finally intersects it with the Script-1 good-spaxel mask.

This aperture threshold is a boundary-definition diagnostic, **not** a per-spaxel science S/N cut. Once a pixel is inside the accepted stellar-body aperture, low native S/N alone does not remove it; PowerBin is allowed to combine faint neighboring pixels until the target bin S/N is reached.

The aperture threshold, smoothing scale, dilation, minimum size, and optional maximum physical radius are configurable and must be inspected on real data rather than treated as fundamental constants.

## 5.4 Apply identical physical memberships to RH3

RH3 is **not** independently PowerBinned.

The BL membership map is transferred to the red/RH3 native grid through the celestial WCS. For the baseline paired observations the two arms use the same slicer and KcwiKit spatial sampling, so each valid red/RH3 pixel center is transformed to BL pixel coordinates and assigned to the nearest BL member pixel only when the sky-plane separation is within the configured tolerance.

Script 2 first checks that the BL and RH3 spatial pixel scales agree within tolerance. If the grids are materially different, nearest-pixel transfer is rejected rather than silently pretending that the apertures are identical; a future area-overlap/polygon transfer would then be required.

For bin $i$, the resulting spectra represent the same geometric sky aperture as closely as the native sampling permits:

$$
F_i^{\mathrm{BL}}(\lambda)
$$

and

$$
F_i^{\mathrm{RH3}}(\lambda).
$$

The achieved RH3 S/N is measured but never drives the bin boundaries.

### Achieved-S/N estimator and pathology handling

The final extracted-spectrum S/N is a QC/reporting statistic, not the PowerBin capacity itself. For each bin and configured continuum window, Script 2 now defines the robust signed diagnostic

$$
(S/N)_{i,\mathrm{signed}} = \frac{\mathrm{median}[F_i(\lambda)]}{\mathrm{median}[\sigma_i(\lambda)]}.
$$

This ratio-of-medians is less sensitive than the former $\mathrm{median}[F_i(\lambda)/\sigma_i(\lambda)]$ estimator to isolated wavelength samples with spuriously tiny formal uncertainties. The older estimator is still saved in the per-bin table as an audit diagnostic.

A non-positive median continuum does not represent a meaningful positive achieved S/N. When `BIN_SN_REQUIRE_POSITIVE_CONTINUUM=True`, the production-facing `BL_SN`/`RH3_SN` value is therefore saved as `NaN` for such a bin, while the signed value, legacy estimator, median flux, median uncertainty, minimum uncertainty, negative-flux fraction, and number of good channels remain explicitly saved. No flux or variance is clipped or rescaled by this rule.

Extreme signed values and large disagreement between the robust and legacy estimators generate quality flags and should trigger inspection with `scripts/inspect_script02_sn.py`.

### S/N-window coverage QC

Before measuring achieved S/N, Script 2 compares each configured rest-frame window with the Script-1 `GOODWAVE` envelope after transforming the requested interval to observed wavelength. This is a configuration/data-consistency check only. The pipeline logs the requested observed interval, the usable interval, the fraction of the requested wavelength width covered, and the fraction of native requested channels that remain usable.

If the requested interval extends beyond the Script-1 usable wavelength envelope, Script 2 raises `SN_WINDOW_COVERAGE_WARNING`. The configured window is **not** moved, optimized, or replaced automatically. This matters for development data such as an RL integration fixture being passed through a pipeline whose production red-arm window is designed for RH3 CaT observations. Production RH3 data should normally place the CaT window well inside the intended grating coverage; the warning exists to catch configuration mistakes and non-production test cases rather than to make poor data look artificially good.

An optional development-only utility, `scripts/scan_script02_sn_windows.py`, can scan fixed-width observed-frame windows in an already completed Script-2 spectral product. It reports the positive-continuum bin fraction, median positive-bin S/N, lower-tail S/N, and extreme-|S/N| fraction. This utility never edits the target config and is not imported by the production pipeline.

## 5.5 Geometric aperture spectra and scalar-calibrated uncertainties

Spatial extraction is a **geometric sum**, not an inverse-variance-weighted spatial average. This preserves the actual light mixture that later bin-integrated disk models must predict.

For one wavelength sample,

$$
F_i(\lambda)=\sum_{p\in\mathcal B_i}F_p(\lambda).
$$

If the input cube is explicitly in surface-brightness units per square arcsecond, the sum is multiplied by the spatial pixel area so the extracted spectrum represents integrated aperture flux. Formal diagonal variance is propagated as

$$
\mathrm{Var}_i(\lambda)=\sum_{p\in\mathcal B_i}\mathrm{Var}_p(\lambda).
$$

A wavelength sample is accepted for the bin only when at least the configured fraction of its spatial members contributes a valid sample. The formal error is multiplied by an empirical spatial scale measured from mock stack sums for that same bin. Both arrays are saved.

The Script-1 preliminary spectral-correlation sequence is copied into the Script-2 spectral product as provenance/QC information. It is **not** converted into a final covariance matrix because that high-pass diagnostic was intentionally preliminary.

## 5.6 Continuum-light weights and bin geometry

For later PSF/bin-integrated modeling, Script 2 saves the complete native membership and a normalized positive continuum-light weight for every member pixel. Within one bin,

$$
\sum_{p\in\mathcal B_i}w_{i,p}=1.
$$

Both geometric and continuum-light-weighted centroids are saved. If all measured continuum values in a bin are non-positive because of noise, equal geometric weights are used as a transparent fallback rather than producing undefined weights.

## 5.7 Per-bin products

The authoritative Script-2 products include:

- `master_bins.fits`: BL master membership, analysis aperture, per-bin BL/RH3 S/N maps, area map, transferred RH3 membership, and WCS-transfer distance;
- `master_bin_spectra.fits`: wavelength, aperture-summed flux, formal uncertainty, valid-sample mask, and number of contributing spaxels for both arms, plus preliminary Script-1 spectral-correlation sequences when available;
- `master_bin_table.ecsv`: one row per bin with member counts, area, geometric/flux centroids, sky position, robust positive-continuum BL/RH3 S/N, signed/legacy S/N audit diagnostics, continuum/uncertainty summaries, PowerBin capacity, RH3 transfer coverage, and single-pixel status;
- `master_bin_membership.npz`: complete BL and RH3 native pixel membership, tangent-plane coordinates, continuum-light weights, transfer distances, and PowerBin generator information;
- `script02_manifest.json`: source Script-1 run, configuration/provenance, covariance assumptions, transfer completeness, quality flags, and product paths.

These products are intentionally generous because later stages need exact membership rather than only a plotted bin map.

## 5.8 Shared S/N plot normalization

BL and RH3 S/N maps use the same color normalization. The default upper limit is the 95th percentile of the **combined** BL+RH3 bin-S/N distribution,

$$
v_{\max}=P_{95}\left(\{S/N_{\mathrm{BL}}\}\cup\{S/N_{\mathrm{RH3}}\}\right),
$$

so a given color means the same S/N in both panels while a few very high-S/N bins do not destroy contrast elsewhere.

## 5.9 Script-2 quality flags

Important non-fatal flags include:

- `BINNING_TARGET_SN_NOMINAL`: PowerBin membership was defined using nominal BL noise before the per-bin mock calibration;
- `BIN_TRANSFER_INCOMPLETE`: fewer red/RH3 pixels than expected could be assigned to the BL-defined physical apertures;
- `LOW_BL_BIN_SN`: the measured S/N of one or more final BL spectra falls substantially below the configured target despite the PowerBin proxy;
- `SN_WINDOW_COVERAGE_WARNING`: a configured achieved-S/N window extends beyond the Script-1 `GOODWAVE` envelope; the pipeline records the truncation but does not change the window automatically;
- `NONPOSITIVE_BIN_CONTINUUM`: one or more configured S/N windows have non-positive median continuum, so a positive achieved-S/N value is undefined for those bins;
- `EXTREME_BIN_SN_DIAGNOSTIC`: a signed robust or legacy S/N diagnostic exceeds the configured numerical warning threshold;
- `BIN_SN_ESTIMATOR_DISAGREEMENT`: the robust ratio-of-medians and legacy median-of-ratios disagree by more than the configured factor;
- `BINNING_APERTURE_WARNING`: the useful-stellar-body aperture requires manual inspection.

A flag is not an instruction to tune parameters until the warning disappears. It records where a scientific assumption or data limitation must be inspected.

---

# 6. Script 3 — RH3 profile-likelihood cubes

## 6.1 Main product

For every BL-defined master PowerBin, construct the complete independent spectral surface

$$
\boxed{\chi_i^2(V_A,V_B,f_{A,\mathrm{RH3}})}.
$$

The baseline grid is

$$
17\times17\times9,
$$

with the velocity and fraction axes defined by the target configuration. The final grid extent/resolution must still pass later convergence tests; Script 3 preserves the full coarse cube so that those tests do not require changing the global-disk logic.

Script 3 does **not** impose spatial coherence. Every PowerBin is interrogated independently. Component identity and the physically coherent trajectory through these independent cubes enter only in Script 4.

## 6.2 Exact meaning of a grid coordinate

At a grid coordinate

$$
(V_A,V_B,f_A),
$$

all three quantities are fixed exactly. In particular, this is intentionally different from the earlier Mitzkus-style development grid in which pPXF could move each velocity by half a grid cell.

For Script 3, pPXF is called with the component velocities fixed and profiles only over nuisance quantities including

$$
\sigma_A,\qquad \sigma_B,
$$

the stellar-template weights within each component, and the additive-continuum coefficients.

Because the velocities are fixed, their pPXF bounds are only tiny bookkeeping intervals around the requested coordinate; Script 3 does **not** give a fixed grid state an artificial broad velocity-search range. This is important both conceptually and for pPXF's automatic `lam`/`lam_temp` template-coverage logic.

The complete $V_A/V_B$ plane is retained, including $V_A=V_B$ and the component-swapped half of the cube. The symmetry

$$
(V_A,V_B,f_A)\longleftrightarrow(V_B,V_A,1-f_A)
$$

is therefore preserved rather than broken locally. Script 4 assigns physical Disk A/B identities from the global spatial solution.

## 6.3 Statistical interpretation: profile likelihood

Because nuisance parameters are optimized at every explicit $(V_A,V_B,f_A)$ coordinate, the Script-3 surface is a **profile likelihood**, not a fully marginalized Bayesian posterior.

The fundamental statistic saved by the cube is the **total** chi-square on the fixed fitted pixels,

$$
\chi^2_{\rm total}
=
\sum_{p\in\mathrm{good}}
\left[\frac{F_p-M_p}{\sigma_p}\right]^2
$$

for the current diagonal-noise development model. pPXF's reduced $\chi^2$ is retained as QC, but downstream relative-likelihood calculations must use differences in total chi-square.

When later writing

$$
w_c\propto\exp(-\Delta\chi_c^2/2),
$$

these are called **relative-likelihood weights**, not posterior probabilities.

## 6.4 RH3 fitting interval and fixed fitted pixels

The production template configuration uses

```python
RH3_FIT_REST_RANGE_ANGSTROM = (8470.0, 8900.0)
```

for the RH3 kinematic fit. This is deliberately separate from `RH3_SN_REST_RANGE_ANGSTROM`, which is a Script-2 continuum/S/N and spatial-light-weight reference window rather than the Script-3 science fit interval.

The current `8143-1902.py` integration fixture contains RL rather than RH3 data, so it temporarily uses a broad interior RL Script-3 fit interval. This is an integration-test override only; the production template keeps the CaT-region RH3 interval above.

For each PowerBin, Script 3 constructs one fixed `goodpixels` set from:

- Script-2 spectral validity;
- the configured Script-3 rest-frame fitting interval;
- the empirically supported Script-1 LSF interval;
- finite positive formal uncertainty;
- optional configured observed-frame sky/telluric masks;
- optional configured rest-frame masks.

That identical pixel set is used by the one-component control and by **every** two-component grid state in the bin. `clean=True` or state-dependent clipping is not allowed in the likelihood cube because allowing different states to discard different pixels would make their total chi-square values non-comparable.

### pPXF interface placeholders for already-excluded samples

The rebinned CRD_DAP science arrays deliberately retain `NaN` at rejected log-grid samples so the masking provenance remains explicit. pPXF imposes a stricter API requirement: it validates the complete input vectors before applying `goodpixels`, requiring the full galaxy vector to be finite and the full noise vector to be finite and positive. Script 3 therefore performs a narrow interface-only sanitization in `crd_utils.ppxf_utils`. Every index in `goodpixels` is first required to contain finite flux and finite positive uncertainty. Only invalid samples that are **already excluded** from `goodpixels` are replaced in private copies passed to pPXF, using benign finite placeholders. The original science/checkpoint arrays are not modified, and the placeholder samples contribute neither to the pPXF optimization nor to CRD_DAP's explicit total $\chi^2$, which is evaluated only on `goodpixels`. If an invalid value ever appears *inside* `goodpixels`, Script 3 hard-fails instead of repairing it.

## 6.5 Wavelength / velocity convention

The Script-2 observed wavelength vector is brought close to the galaxy rest frame using the configured target redshift. The science wavelength medium is validated from the prepared Script-1 product and converted to the configured XSL/template medium when required.

The Script-3 velocity coordinates therefore represent residual LOS velocities relative to the adopted target redshift. Any small error in that initial systemic redshift is handled later by the free global $V_{\rm sys}$ parameter in Script 4.

The galaxy and templates are supplied to pPXF with their explicit wavelength arrays. No additional `vsyst` offset is introduced in this rest-frame construction.

## 6.6 Full XSL SSP basis and empirical LSF matching

Use the **full XSL SSP grid** for RH3 rather than reducing the basis purely for speed.

The RH3-prepared template matrix is produced by:

1. loading the complete configured XSL SSP library;
2. computing template padding from the **actual fitting domains**: the most extreme exact two-component grid velocity, the full one-component control velocity bounds, the configured multiple of the maximum allowed dispersion, and a small log-grid edge-safety margin;
3. cropping to that padded wavelength range;
4. evaluating the Script-1 empirical wavelength-dependent RH3 LSF only inside its measured arc-line support;
5. broadening the XSL templates in quadrature where the data LSF is broader;
6. applying the fixed light-fraction normalization convention below;
7. log-rebinning at the same velocity scale as the galaxy spectrum.

The pipeline must **not** take an absolute value when

$$
\mathrm{FWHM}_{\rm data}^2-\mathrm{FWHM}_{\rm template}^2<0.
$$

If the XSL SSP basis is intrinsically broader than the real RH3 data over part of the required interval, Script 3 currently hard-fails and reports the incompatibility rather than silently degrading or misrepresenting the spectral resolution. A future production treatment may broaden the galaxy to a validated common target LSF, but that operation must also propagate the resulting noise covariance and should be adopted explicitly rather than hidden inside template preparation.

Before any `ProcessPoolExecutor` workers are launched, Script 3 performs a real **pPXF wavelength-coverage preflight** using the installed pPXF version. It selects the small set of bins that collectively span the earliest, latest, and widest good-pixel wavelength support; runs the one-component control; and evaluates the four extreme two-component velocity-grid corners at a representative fraction. The log records the galaxy and prepared-template wavelength endpoints, the padding budget, the test bins, and PASS/FAIL. A failure aborts before the expensive cube calculation begins.

## 6.7 Definition of $f_{A,\mathrm{RH3}}$

This convention is scientifically important and must remain explicit.

Before the XSL basis is duplicated into components A and B, **each SSP spectrum is independently normalized to unit mean stellar flux density over the fixed Script-3 RH3 fitting interval**:

$$
\left\langle T_j(\lambda)\right\rangle_{\lambda\in\mathrm{RH3\ fit\ band}}=1.
$$

Both components then receive identical copies of this normalized SSP basis. pPXF's two-component `fraction` constraint is

$$
f_{A,\mathrm{RH3}}
=
\frac{\sum_j w_{A,j}}
{\sum_j w_{A,j}+\sum_j w_{B,j}}.
$$

Because all SSPs share the same passband normalization, this has the intended interpretation:

> **$f_{A,\mathrm{RH3}}$ is the fraction of fitted stellar-template light assigned to component A over the RH3 fitting passband.**

It is **not** a stellar-mass fraction. The additive polynomial is a nuisance continuum term and is not counted as either disk's stellar light. Script 6 defines an analogous BL light fraction over the usable BL passband, and therefore

$$
f_{A,\mathrm{RH3}}\neq f_{A,\mathrm{BL}}
$$

is physically allowed and potentially informative about different disk stellar populations.

The achieved fraction reconstructed from the returned component template weights is saved as QC and compared with the requested grid fraction.

## 6.8 Dispersion and continuum nuisance parameters

The baseline Script-3 dispersion settings are

```python
RH3_SIGMA_START_KMS = 60.0
RH3_SIGMA_MIN_KMS = 5.0
RH3_SIGMA_MAX_KMS = 250.0
```

The extended 250 km/s upper bound is intentionally generous. If a state wants a very broad component, the pipeline should expose that behavior rather than skewing the profile likelihood by forcing it against the older 180 km/s simulation bound. States within the configured warning distance of either sigma bound are flagged and counted.

The baseline kinematic continuum treatment is

```python
RH3_DEGREE = 4
RH3_MDEGREE = 0
```

and

```python
RH3_REGUL = 0.0
```

for every fit whose chi-square enters the likelihood cube.

## 6.9 Scalar uncertainty calibration

The reduced-exposure `UNCERT` arrays specify independent Gaussian input draws for a fixed KcwiKit stacking experiment. The science shifts, masks, drizzle setting, exposure selection, and stack geometry are replayed exactly. In Script 2 each mock stack is summed through **every actual PowerBin**, using the same Script-1 hard-good samples and the same spatial membership. For bin $i$, the spatial multiplier is the square root of the median across usable science-window channels of the empirical mock variance divided by the formal binned variance. It multiplies the entire native-grid error spectrum for that bin without changing its wavelength dependence. The formal error and multiplier are saved separately. The BL PowerBin tessellation is made from nominal input S/N; achieved S/N uses the calibrated errors.

This experiment measures the spatial coaddition penalty of the fixed reduction and stacking process. It does not repeat or assign uncertainty to detector calibration, the fitted ZAP sky model, telluric correction, or measured shifts. The red ZAP exposures were cropped by KSkyWizard to `WAVGOOD0`–`WAVGOOD1` before sky subtraction; the new stack’s shorter wavelength range is inherited from that crop.

Script 3 rebins these spatially corrected RH3 spectra to its fixed log-wavelength experiment. For every PowerBin it evaluates the complete, unregularized $V_A\times V_B\times f_A$ grid once using the Script-2 spatially calibrated error spectrum and a fixed good-pixel set. It refits the exact best grid state to retain the model spectrum, then measures normalized residuals $z=(g-m)/\sigma_{spatial}$. The residual multiplier is $q_i=1.4826\,\mathrm{median}|z-\mathrm{median}(z)|$. Values above and below one are both used; no floor is applied. This estimator is conditional on the selected best-fitting model and can itself be affected by fit flexibility or model mismatch. Every bin is fitted; none is chosen as representative.

The final error spectrum is $\sigma_{fit,i}=q_i\sigma_{spatial,i}$. Since $q_i$ is uniform across the fitted wavelengths of one bin, and the grid fits are unregularized, each state's total and reduced chi-square is divided by $q_i^2$ without rerunning the 2601 grid states. The best grid location is unchanged. The saved total chi-square is $\sum_{j\in good}[(g_j-m_j)/\sigma_{fit,j}]^2$; the grid profiles over the two dispersions, template weights and additive continuum as before. The saved `noise` vector is this final likelihood denominator, while `noise_before_residual_scale` is the vector supplied to pPXF during the grid fit. Neither scalar should be applied again. The calibration product records the unfloored scales, residuals and an input hash. Per-bin checkpoints include those diagnostics and use a new schema version; older checkpoints cannot be resumed into this method.

No spectral lag coefficients or dense covariance matrices enter the likelihood. Script-1 high-pass autocorrelation remains an observational diagnostic, not a correction or a condition that must be driven to zero. The spatial and residual scalars calibrate uncertainty amplitude; they do not remove correlation between wavelength samples. Profile-likelihood intervals are therefore conditional on the adopted diagonal approximation. The likelihood-width publication flag remains false while the final data and inference checks are incomplete; it does not require a future spectral-covariance model.

## 6.10 Restartability and CPU workers

Script 3 writes one atomic checkpoint per completed bin and verifies configuration and noise-experiment hashes on resume. A fresh Script-2 run is required when changing the stacked science cube, atmospheric mask or mock ensemble. `--workers N` controls the expensive likelihood-grid process pool; the best-state residual calibration occurs within each completed bin. Checkpoints are deleted only after successful consolidation when the configured cleanup switch is enabled.

# 7. Script 4 — coherent two-disk model

The implemented driver is scripts/04_fit_global_two_disk_model.py. See
[SCRIPT04.md](SCRIPT04.md) for commands and the required human review steps.
The design below supersedes earlier plans for this stage. It fits the saved
Script-3 spectral likelihoods; it does not call XookSuut, refit pPXF, or produce
posterior samples.

## 7.1 Coordinates, geometry and permanent disk labels

Both disks share a fitted center, astronomical PA (east of north), inclination,
and systemic velocity. Each has its own nonnegative intrinsic radial speed
magnitudes. Their rotation signs remain opposite and fixed at every radius.
The configurable numerical ceiling is initially 400 km/s per disk; it is not a
physical maximum or the same quantity as the saved LOS velocity-grid bound.

The permanent A convention refers to the outer branch on the negative signed
major-axis side. Its provisional sign follows the outermost reference bin's
single-component velocity and is flagged for review. An explicit reviewed sign
can replace it. Local light dominance never reassigns A/B labels.

Existing center, PA and systemic-velocity ranges are hard bounds. Inclination
bounds come from a deliberately selected PyMorph disk-axis-ratio table and the
existing q/q0 propagation, or explicit reviewed inclination inputs. The catalog
table/band is not silently selected. The inherited q/q0 sampler uses Gaussian
draws with physical clipping, which is recorded in provenance. These bounds are
not new posterior measurements.

## 7.2 Fixed radial grids and whole-bin domains

The first free node is one adopted RH3 PSF FWHM from the center. Default node
spacing is that FWHM, with nominal annulus half-width half the spacing. Set
U(0)=0 and interpolate linearly through the free nodes. Hold the last-node
speed constant over outer padding; do not create an unmeasured extra node.

The restricted domain starts from the larger of the two user-supplied MaNGA
dispersion-peak distances, snapped to the nearest node. Include complete bins
touching the outer annulus, then extend nodes outward to cover those complete
bins. The full domain includes every upstream bin and rounds node coverage
outward. Membership, node dimension and reference geometry are fixed before
optimization. Changing trial geometry cannot improve a score by dropping bins.
The actual restricted extent may exceed its nominal edge because of a large
crossing bin.

MaNGA peak and baseline dispersions give the equal-light initial LOS speed
sqrt(sigma_peak^2 - sigma_base^2), or half the inferred disk separation.
Correct by sin(i) for an intrinsic-speed seed. Vary starting shapes, amplitudes,
A/B asymmetries and geometry; these seeds add no likelihood term.

## 7.3 Measured stellar PSF and fixed light proxy

Fit both Gaussian and Moffat profiles to Feige 24 continuum images in the RH3
continuum interval, integrating over native rectangular KCWI pixels. The user
inspects residuals and profiles before choosing a candidate. Rotate the
pixel-integrated profiles into the science orientation and average unit-flux
kernels. The median larger-axis effective FWHM sets the radial-node scale.

This is an effective standard-star PSF transfer approximation. The standards
are not contemporaneous with the science exposures, and the exact science
stack/drizzle response is not replayed. Preserve those limitations as flags.
Do not interpret this as a precise recovered science-stack PSF.

Use the observed RH3 continuum as the fixed common light proxy I for both
disks. Collapse the recorded Script-2 continuum window respecting the prepared
masks; set negative continuum weights to zero as in Script 2. Preserve global
relative brightness before convolution, rather than convolving independently
normalized bin weights. No intrinsic or component-specific light model is fitted.

For each projected velocity field V, predict first moments as:

    I_model = P * I
    V_model = [P * (I V)] / I_model

Sum numerator and denominator through each entire bin aperture. For the
no-PSF mode, use direct continuum-weighted aperture means. Save the corresponding
second-moment velocity shear as a diagnostic; do not use it to replace pPXF's
component dispersions. Report received light originating beyond the final node
and flag spatial edges or missing light. No flux is invented beyond the image.

The observed-light substitution follows a precedent in
[DiGiorgio Zanger et al. (2024), section III, equations 4–6](https://arxiv.org/html/2407.11908v1#S3).
It remains an approximation for these two-component spectra.

## 7.4 Data objective and fraction profiling

For each trial global model, evaluate every saved RH3 fraction plane by bilinear
interpolation in V_A and V_B; only then minimize over fraction. Profiling corner
minima first could mix incompatible fractions and artificially lower a score.
Linear interpolation in fraction cannot create a lower minimum between the
saved planes, so exact finer fraction fitting remains Script 05's task.

Subtract each bin's independent minimum once, then sum the bin contributions:

    Q_data = sum_i min_f Delta_chi2_i[V_A,i, V_B,i, f]

No extrapolation or boundary clipping is allowed. A trial lacking a supported
state for any selected bin is rejected in full; no bins disappear from its
objective. Failed zero-weight interpolation corners do not invalidate an exact
valid grid state.

Read the final Script-3 scalar-calibrated chi-square surfaces unchanged. Do not
apply either noise multiplier again. Both wavelength independence within bins
and independence of bin contributions remain declared working approximations.
The fitting objective is not a covariance-corrected likelihood.

## 7.5 Experimental curvature penalty

Following [Sellwood & Zánmar Sánchez (2010), Appendix A](https://arxiv.org/html/0912.5493v2#A1),
apply squared second differences to the two independent radial speed curves:

    P = lambda A sum_disks sum_j (U_j - 2 U_(j+1) + U_(j+2))^2
    A = (Rmax / DeltaR)^4 / Vtyp^2
    C = lambda A

The first retained triple is (U1,U2,U3): omit the triple involving the imposed
zero-speed origin. This central omission is our explicit choice, not a
prescription attributed to the paper. Remaining terms can still influence U1.
There are no artificial zero-slope endpoint conditions.

C is fixed across extent comparisons; lambda_extent=C/A_extent. Keep common
spacing and fixed Vtyp=200 km/s for reporting. Vtyp does not constrain the fitted
speeds. C has units (km/s)^-2 for the implemented velocity units. Outer models
can add penalty terms, so their total penalty need not equal the inner penalty.

Optimize Q_data + P, retaining both parts separately. The unpenalized control
uses C=0. These are regularized optimization solutions, not evidence values.

## 7.6 Manual strength calibration and eight comparisons

Before the final comparison, fit a full-aperture PSF-weighted control and eleven
positive strengths C=0.002*2^j for j=-5 through +5. The literature-motivated
starting value is a trial, not a validated strength for our data objective.

The single-page two-row by six-column sheet shows curves, per-bin data-cost
changes, geometry and changes, separate scores, flags, and descriptive RMS
curve distances. The user can shift the range weaker/stronger or select a
tested coefficient, including zero. Custom values must be fitted and reviewed.
Record the selection rationale and freeze C before the final comparison.

Run all combinations of restricted/full extent, PSF off/on, and curvature
off/on: eight configurations. Matching pilot fits can be reused. Compare
extent pairs on common inner bins and radii, report outer-bin costs separately,
and retain parameter counts and mean per-bin costs. These means are descriptive;
do not call them reduced chi-square or normalize the optimization objective by
bin count. No automatic model preference is assigned.

## 7.7 Search, provenance and downstream contract

Use differential-evolution population exploration followed by multiple bounded
Powell refinements of promising starts, evaluating the same objective throughout.
Save alternative local solutions, search histories, convergence/boundary flags,
geometry, radial curves, per-bin LOS velocities, coarse fractions and scores.
Numerical search budgets are configurable and their exhaustion remains visible.
Search trajectories are not posterior draws or rotation-curve error bars.

Parallelize independent entire fits, with one numerical-library thread per
worker. A completed fit is an atomic checkpoint. Verify code/input/config
fingerprints before reuse. Changed science inputs or search settings require a
new output directory; explicit upstream paths support transferred Windows runs
without rewriting original provenance.

The user selects one of the eight models explicitly. The selected-model product
points to its joint geometry/curves and bin predictions for exact Script-05 pPXF
refits. Final uncertainties are a later inference task, conditional on the
declared noise and spatial-model assumptions. Successful software tests do not
establish that the low-quality development data constrain a physical solution.

References and reusable BibTeX are retained in references/SCRIPT04_REFERENCES.md
and references/script04.bib. The ring framework is inspired by
[XookSuut](https://arxiv.org/abs/2110.05095); our spectral-cube objective and
optimizer are a custom implementation, not that package's Bayesian machinery.

---

# 8. Script 5 — exact RH3 solution and quality assessment

## 8.1 Final exact pPXF refits

The global disk model predicts continuous velocities that generally do not land exactly on the Script-3 grid.

For every bin, rerun pPXF at

$$
V_{A,i}^{\mathrm{global}},
\qquad
V_{B,i}^{\mathrm{global}}
$$

rather than reporting only interpolated nuisance parameters from the coarse grid.

This removes grid quantization from final RH3 measurements.

## 8.2 Fine RH3 light-fraction refinement

The coarse fraction grid maps the global likelihood topology, but the published RH3 light fraction should not be quantized at 0.1.

Around the coarse preferred fraction, rerun an explicit fine grid with default step

$$
\Delta f_{A,\mathrm{RH3}}=0.01.
$$

A local interval is tried first. If the fine minimum reaches an interval edge, expand the interval automatically; the full physical fraction range can be scanned if needed.

The script should print the number of additional pPXF fits, average fit time, and estimated remaining runtime.

## 8.3 Grid-edge warnings

Any selected RH3 state within two grid cells of a velocity or fraction-grid edge receives a warning. This identifies truncation risk: the true likelihood valley may extend beyond the computed cube.

No extrapolation beyond a profile-likelihood cube is allowed.

## 8.4 Basin identification

The global model provides the continuous $(V_A,V_B)$ location. For discrete basin bookkeeping, locate the nearest velocity-grid cell and choose the RH3 fraction that minimizes $\chi^2$ at that velocity coordinate. This defines the discrete **anchor cell**.

The anchor need not be the global minimum of the individual RH3 cube.

Starting from the anchor, follow decreasing $\chi^2$ through neighboring 3-D cells until reaching a local minimum. Every finite cell whose downhill path terminates at that same minimum belongs to the selected basin.

The baseline implementation can use full immediate 3-D connectivity (26 neighbors in a 3-D grid), but this must be validated on injection/recovery likelihood cubes before publication.

A basin is a topological object. It may contain the whole cube if only one local minimum exists.

## 8.5 Relative-likelihood weights inside the basin

For selected-basin cell $c$, define

$$
\Delta\chi_c^2=\chi_c^2-\chi_{\mathrm{basin,min}}^2
$$

and

$$
w_c=\exp\left(-\frac{\Delta\chi_c^2}{2}\right).
$$

Normalize within the selected basin:

$$
p_c=\frac{w_c}{\sum_{j\in\mathcal B}w_j}.
$$

These are **normalized relative-likelihood weights conditional on the selected basin**, not formal Bayesian posterior cell probabilities.

The diagnostic counts of cells needed to enclose 90%, 95%, 99%, 99.5%, and 99.9% of the normalized likelihood are retained as measures of likelihood concentration only. They are not hard science cuts.

## 8.6 Why a large basin is not necessarily a problem

A basin can occupy the whole cube while nearly all likelihood mass remains near its minimum. Basin extent describes topology; likelihood concentration describes statistical precision.

This distinction is central to Script 6.

## 8.7 One-vs-two-component mock calibration

For the final selected/refined two-component solution, compute the one-vs-two-component statistic

$$
T=\chi^2_{1\mathrm{comp}}-\chi^2_{2\mathrm{comp}}.
$$

Use mocks to derive two distinct quantities:

### Empirical false-positive probability

From true one-component mock spectra,

$$
p_{\mathrm{false}}=P(T\ge T_{\mathrm{obs}}\mid\mathrm{true\ one\ component}).
$$

### Two-component recovery probability

From injected true two-component mocks matched/interpolated in

$$
S/N,\quad\Delta V,\quad f_A,\quad\sigma_A,\quad\sigma_B,
$$

estimate

$$
p_{\mathrm{recover}}=P(\mathrm{successful\ recovery}\mid\mathrm{true\ two\ component,\ conditions}).
$$

Do not interpret either quantity as a Bayesian probability that the observed galaxy “is two-component.”

---

# 9. Script 6 — BL two-component stellar populations

## 9.1 BL stellar-light fraction definition

Define the BL fraction over the globally adopted usable BL fitting band after consistent template normalization:

$$
f_{A,\mathrm{BL}}=\frac{\sum_jw_{A,j}}{\sum_jw_{A,j}+\sum_jw_{B,j}}.
$$

The fraction is purely stellar. Gas-template weights are excluded.

Because the two disks may have different stellar populations,

$$
f_{A,\mathrm{BL}}\ne f_{A,\mathrm{RH3}}
$$

is physically allowed and potentially informative.

## 9.2 Gas emission

Fit gas emission simultaneously in BL where relevant:

$$
F_{\mathrm{model}}=F_{\star,A}+F_{\star,B}+F_{\mathrm{gas}}.
$$

Gas has its own kinematics and does not contribute to the stellar disk light fraction.

## 9.3 Continuum treatment

Keep additive and multiplicative polynomial choices explicit configuration parameters.

Baseline likelihood fit:

- no additive polynomial (`degree=-1`);
- modest multiplicative polynomial, degree validated on mocks.

The purpose is to absorb continuum/flux-calibration mismatch without allowing additive terms to arbitrarily alter absorption-line depths used for stellar-population inference.

## 9.4 BL fraction search

Use a two-stage explicit fraction search:

1. coarse grid with approximately $\Delta f_A=0.1$;
2. local refinement near the preferred region with approximately $\Delta f_A=0.01$.

If the fine solution reaches the local refinement edge, expand automatically.

## 9.5 Three population products

### A. RH3-anchored solution

Fix the exact global RH3 kinematics and fit BL populations while profiling the BL stellar fraction.

This answers: **What populations are inferred if the best RH3 kinematics are adopted exactly?**

### B. RH3-likelihood-propagated solution

Propagate the selected RH3 basin's normalized relative-likelihood distribution into the BL analysis.

This answers: **How do the BL population inferences change across the family of kinematic states supported by RH3?**

### C. Joint profile solution

Search for the lowest joint statistic

$$
J(c,f_{A,\mathrm{BL}})=\Delta\chi^2_{\mathrm{RH3},c}+\Delta\chi^2_{\mathrm{BL}}(c,f_{A,\mathrm{BL}})
$$

within the selected physical basin, with deterministic local refinement around the best sampled state.

This gives the best discrete joint RH3+BL profile solution without relying on random sampling to happen to hit the exact local joint minimum.

## 9.6 Direct evaluation versus likelihood-weighted sampling

If the selected RH3 basin contains at most `N_DIRECT` cells (development default 200), evaluate every state exactly.

If the basin is larger, draw RH3 states according to

$$
p_c\propto e^{-\Delta\chi_c^2/2}.
$$

Conceptually, each basin cell has “tickets” in proportion to its RH3 relative-likelihood weight.

Every draw is independent from the same fixed RH3 distribution. The **draws themselves do not build on one another**. What accumulates is the population inference calculated from the cumulative draw ensemble.

After each new sampling batch, recompute population summaries using **all draws obtained so far**.

## 9.7 Avoid double-counting RH3

If states are drawn according to $p_{\mathrm{RH3}}(c)$, RH3 has already entered through the sampling frequency. Do **not** multiply each sampled draw by $L_{\mathrm{RH3}}$ again.

A repeated cell may require only one expensive BL pPXF calculation; its multiplicity is retained when constructing the cumulative sample.

## 9.8 BL weighting within sampled states

For a sampled RH3 state, BL supplies

$$
L_{\mathrm{BL}}(c,f_A)
\propto
\exp\left[-\frac{\Delta\chi^2_{\mathrm{BL}}}{2}\right].
$$

Because the proposal draws already follow the RH3 distribution, the BL likelihood reweights the sampled ensemble toward the joint RH3+BL target.

## 9.9 Sampling convergence

After a minimum number of RH3 draws, add batches and recalculate the weighted 16th, 50th, and 84th percentiles of scientifically important quantities, including

$$
\log Age_A,\quad\log Age_B,
$$

$$
[M/H]_A,\quad[M/H]_B,
$$

$$
f_{A,\mathrm{BL}},
$$

$$
\Delta\log Age,\quad\Delta[M/H].
$$

For each quantity, define the current 68% interval width

$$
W=P_{84}-P_{16}.
$$

The tolerance is

$$
T=\max(T_{\mathrm{floor}},\epsilon W),
$$

with development default $\epsilon=0.05$.

Require changes in $P_{16},P_{50},P_{84}$ to remain below tolerance for several consecutive cumulative-batch comparisons (development default three).

A human-readable rounded-value stability summary may also be logged, but literal decimal rounding is not the formal stopping rule.

## 9.10 Maximum draw safety limit

Sampling is bounded by `RH3_MAX_DRAWS`. If the convergence criteria are not met:

- save all partial state/results/caches;
- raise/catch a convergence exception;
- mark the bin `BL_JOINT_NOT_CONVERGED`;
- continue the larger run when practical;
- return a non-zero overall status or explicit summary requiring inspection.

If every basin cell has been evaluated, the result becomes exact over the discrete basin and sampling stops regardless of the batch criterion.

## 9.11 Deterministic local refinement of the joint profile minimum

Sampling identifies a promising region. Starting from the best sampled joint state:

1. evaluate unevaluated neighboring RH3 cells;
2. profile/refine $f_{A,\mathrm{BL}}$;
3. move to a neighbor if it decreases $J$;
4. repeat until no neighbor improves $J$.

All BL fits are cached, so previously evaluated states are never repeated unnecessarily.

## 9.12 Template-grid edge warnings

Flag population solutions that accumulate against SSP-library boundaries in age or metallicity. A boundary-limited solution should not be reported as if its mean age/metallicity were fully constrained by the data.

## 9.13 Population-difference distributions

For every weighted solution, calculate directly

$$
\Delta\log Age=\log Age_A-\log Age_B
$$

and

$$
\Delta[M/H]=[M/H]_A-[M/H]_B.
$$

This preserves correlations between the two disk population measurements better than subtracting two separately summarized medians after the fact.

---

# 10. Script 7 — uncertainties, systematics, and final maps

## 10.1 Statistical uncertainty should propagate through the full chain

The gold-standard realization is end-to-end:

$$
\mathrm{perturb\ RH3}
\rightarrow
\mathrm{RH3\ likelihood}
\rightarrow
\mathrm{global\ disk\ model}
\rightarrow
\mathrm{exact\ RH3\ extraction}
\rightarrow
\mathrm{perturb\ BL}
\rightarrow
\mathrm{population\ fit}.
$$

This captures the fact that BL populations are conditional on uncertain RH3 kinematics.

## 10.2 Noise realizations

Noise realizations for the final science uncertainty stage must preserve the measured spatial and residual scalar calibrations and be checked by injection/recovery. Script 1 spectral autocorrelation is diagnostic only; it does not define a production spectral covariance model.

## 10.3 Monte Carlo count is convergence-driven

Do not choose an enormous fixed value such as $10^5$ by default. Monte Carlo sampling error falls approximately as

$$
N_{\mathrm{MC}}^{-1/2}.
$$

A development strategy is:

- minimum: about 200 realizations;
- check every 50–100 new realizations;
- expected final range: a few hundred to roughly 1000 if convergence supports it;
- safety maximum: about 2000 until real timings and convergence behavior are known.

The final number is determined from both statistical convergence and measured runtime.

## 10.4 Parallelism

Nominal bin-based stages may parallelize across PowerBins.

The end-to-end Monte Carlo should preferentially parallelize independent realizations across workers and avoid nested uncontrolled parallelism inside each realization.

After the first complete nominal pipeline run, benchmark one or more full realizations and forecast wall time for the user's available core count (e.g. 3 versus 20–40 cores).

## 10.5 Three uncertainty layers

Keep distinct:

### Statistical uncertainty

Noise/covariance realizations through the full inference chain.

### Analysis-choice robustness

Reasonable alternatives such as:

- ring spacing;
- $R_{\mathrm{start}}$;
- $R_{\mathrm{final}}$;
- 2σ-limited versus full-aperture model;
- inclination prior;
- PSF uncertainty;
- LSF uncertainty;
- polynomial degree;
- fraction-grid resolution;
- velocity-grid resolution;
- shared versus mildly relaxed geometry.

### SPS/template systematic

Repeat representative fits with a suitable alternative population model/template assumption when possible.

Do not automatically collapse these three categories into one quadrature error bar. They answer different scientific questions.

## 10.6 Final spatial maps

Every BL PowerBin receives age/metallicity measurements for Disk A and Disk B. Therefore these are valid 2-D maps at PowerBin resolution:

$$
Age_A(x,y),\quad Age_B(x,y),
$$

$$
[M/H]_A(x,y),\quad [M/H]_B(x,y).
$$

Also map disk contrasts

$$
\Delta\log Age(x,y),\qquad\Delta[M/H] (x,y)
$$



and method-robustness differences between anchored and likelihood-propagated solutions.

## 10.7 Spatial correlation / multiple comparisons

Neighboring bins are not fully independent because of the PSF and the shared global RH3 model. Scientific interpretation should emphasize coherent spatial/radial patterns and galaxy-wide component differences rather than treating every bin as an independent hypothesis test.

If binary significance maps are eventually used, false-discovery/multiple-comparison treatment should be considered explicitly.

## 10.8 Final machine-readable table

The final per-PowerBin table should contain at least:

- bin ID;
- sky/bin coordinates;
- deprojected radius;
- BL and RH3 S/N;
- $V_A,V_B$ and uncertainties;
- $\sigma_A,\sigma_B$ and uncertainties;
- $f_{A,\mathrm{RH3}}$ and uncertainty/likelihood summaries;
- $f_{A,\mathrm{BL}}$ and uncertainty/likelihood summaries;
- age and metallicity for both disks;
- direct disk contrasts;
- one/two-component mock-calibration metrics;
- analysis mode (direct/sampled) for Script 6;
- convergence metrics;
- standardized quality flags.

An Astropy ECSV/FITS table is preferred internally because it preserves metadata and units. A journal machine-readable-table format can be generated at publication time.

---

# 11. Profile likelihood, relative likelihood, and terminology

The RH3 grid is a profile likelihood because nuisance parameters are optimized at each explicit coordinate. Therefore

$$
w_c=e^{-\Delta\chi_c^2/2}
$$

is a **relative-likelihood weight**.

It is useful for:

- ranking states;
- describing likelihood concentration;
- likelihood-weighted RH3 state sampling;
- constructing a practical RH3-supported solution family.

It should not automatically be called a Bayesian posterior probability unless a full posterior model, including priors and nuisance-parameter marginalization, has actually been constructed.

Formal final uncertainties rely heavily on end-to-end noise realizations and robustness tests.

---

# 12. Grid-resolution and edge validation

The nominal 17×17×9 RH3 grid is a development choice. Representative bins must be repeated on finer grids to test whether:

- selected global velocities change;
- basin topology changes;
- likelihood concentration changes;
- Script-6 population summaries change.

A selected solution within two cells of an edge in $V_A$, $V_B$, or $f_A$ is automatically flagged.

---

# 13. Data-storage philosophy

Do not intentionally discard scientifically useful products merely to minimize disk usage. Before a large production run, estimate storage requirements so sufficient capacity is available.

Scalar 17×17×9 cubes are cheap. Full template-weight arrays and model spectra can become large, but may be worth preserving for reproducibility and diagnosis. Storage policy should be informed by actual measured output size, not imposed in advance as a science-limiting approximation.

---

# 14. Runtime philosophy

Runtime is a practical constraint, not a scientific objective.

The pipeline should:

- cache repeated expensive fits;
- parallelize naturally independent work;
- print measured per-fit/per-bin timings;
- estimate remaining runtime;
- benchmark end-to-end MC scaling;
- fail safely rather than run indefinitely.

Approximations made solely for speed must be validated against a more complete calculation.

---

# 15. Current choices that remain calibration parameters

The following are intentionally **not frozen scientific truths**:

- intrinsic disk thickness $q_0$ prior;
- exact geometry-prior widths;
- BL multiplicative-polynomial degree;
- final RH3 velocity-grid spacing/extent;
- exact one/two-component mock calibration and reliability thresholds;
- final Monte Carlo realization count;
- final worker/core count;
- whether full spatial covariance or a validated approximation is computationally practical;
- whether the 26-neighbor basin definition is optimal for real likelihood surfaces.

These should be resolved using simulations, real diagnostic plots, and runtime benchmarks rather than guessed in advance.

---

# 16. Definition of “ready for science”

A target should not be treated as publication-ready merely because all scripts execute. The following should be true:

- wavelength conventions are verified;
- noise-scale calibration has converged, and the diagonal spectral-noise assumption is stated explicitly;
- empirical LSF/template compatibility is verified;
- grid-edge warnings are resolved or understood;
- velocity-grid convergence has been demonstrated;
- 2σ-limited and full-aperture global solutions have been compared;
- one/two-component support has been mock-calibrated;
- BL population fits do not unknowingly sit on template-grid boundaries;
- Script-6 propagation converges or exact basin enumeration is achieved;
- end-to-end MC uncertainties converge;
- major analysis-choice/systematic tests have been performed;
- all important diagnostics have been inspected and documented.

Only then should final maps, radial trends, and machine-readable measurements be treated as the scientific result.

### Script-3 state-failure bookkeeping

The profile-likelihood surface is constructed only from states for which pPXF
returns a valid solution at the requested exact velocities. State-level pPXF
exceptions and fixed-velocity mismatches are not assigned finite chi-square;
they are recorded as distinct fit-status classes. If an entire bin contains no
valid state, Script 3 aborts rather than constructing a fictitious likelihood
surface, and reports the underlying pPXF exception frequencies for diagnosis.
This bookkeeping is numerical quality control and does not alter the scientific
likelihood definition.
