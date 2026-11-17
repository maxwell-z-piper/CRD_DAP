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

**Implementation status, reviewed 2026-11-07.** Scripts 01–05 are executable. Sections 4–8 describe the implemented stages. Sections 9–10 document the planned population and uncertainty stages. Walkthrough numbers describe execution order.

Use [PIPELINE_CALLS.md](PIPELINE_CALLS.md) for commands, arguments and Script 04 review gates. Section 7.8 defines the calibration-log fields.

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

Scripts 01-03 use a master `pipeline.log` and a stage-specific log. Script 04 currently uses `logs/04_fit_global_two_disk_model.log` for its main workflow and `logs/<fit_name>.progress.log` for each independent optimizer fit. Main workflow messages go to the terminal and file; the per-fit search trace is written to its own file. These traces are optimization histories, not posterior samples. Section 7.8 defines their fields; Section 7.9 records the implemented log header and the new-search restart requirement.

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

## 4.0 Execution walkthrough - implemented

Driver: `scripts/01_prepare_and_register_cubes.py`. This stage prepares already reduced cubes; it does not perform detector reduction or stacking.

1. **Establish the run and inputs.** Load the target config, create the run directory, snapshot configuration/provenance and initialize logs. When `--reduction-manifest` is supplied, validate the CRD_DRP handoff and resolve its science stacks, atmospheric masks and noise-experiment provenance; otherwise use the explicit legacy/config input path.
2. **Load both arms and validate their data.** Read each KcwiKit flux/variance/mask/exposure quartet, or the supported native-DRP input, and check array shapes, WCS, instrument metadata and wavelength/velocity conventions. Construct hard-good spectral, spatial and wavelength masks from valid flux, positive variance and exposure/quality information.
3. **Apply inherited atmospheric exclusions.** For a reduction-manifest run, combine the supplied atmospheric exclusions with the native-good masks while retaining their separate provenance. Later scripts inherit this mask state rather than reconstructing it.
4. **Build continuum images and center diagnostics.** Collapse usable continuum channels, estimate the smoothed peak and centroid in each arm, transform them to sky positions and compare the two arms.
5. **Check registration and save common coordinates.** Compare the arms through celestial WCS, using matched continuum wavelengths where possible and a bounded morphology cross-correlation only when trustworthy. Save the chosen common sky origin, registration assessment and east/north coordinate grids; do not resample the science cubes or silently apply an inconclusive shift.
6. **Write data-quality diagnostics.** Plot spatial validity, wavelength coverage and exposure/mask information so the accepted footprint and exclusions can be inspected.
7. **Record the initial PSF provenance.** Adopt configured/header seeing values when available and save their source or missing-value flags. This is not a fit to the extended galaxy and is distinct from Script 04's later standard-star PSF measurement.
8. **Measure the empirical spectral LSF.** Validate the master arcs against the science setup, fit usable calibration lines, summarize wavelength/spatial behavior and preserve the accepted-line support and coverage warnings.
9. **Measure preliminary noise diagnostics.** Record high-pass residual/scatter and wavelength-correlation diagnostics. These do not fit spectral lag coefficients or set the later scalar likelihood calibration.
10. **Save the prepared handoff.** Write `prepared_BL.fits`, `prepared_RH3.fits`, the two spatial-coordinate products, LSF/noise products and metadata/quality summaries. Inspect the diagnostics before using this run as Script 02's source.

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

The adopted spectral-fit likelihood uses a wavelength-dependent, scalar-calibrated uncertainty vector and a **diagonal covariance between wavelength samples**. It therefore treats errors at different wavelengths as uncorrelated for the fit; this is an analysis assumption, not a claim that the reduced spectra are physically free of spectral correlation. No spectral lag coefficients or dense covariance matrices are fitted. Script 1's high-pass correlation measurement is a preliminary diagnostic and does not change the likelihood or set its noise scale. Scripts 2 and 3 calibrate the uncertainty amplitude as described below. Any profile-likelihood interval is conditional on this diagonal approximation.

---

# 5. Script 2 — BL-defined master PowerBins

## 5.0 Execution walkthrough - implemented

Driver: `scripts/02_make_master_BL_bins.py`. BL defines one spatial tessellation that is then applied to both arms.

1. **Resolve and record the upstream run.** Load config, select the completed Script 01 run, snapshot provenance and report inherited quality/atmospheric-mask information.
2. **Read the frozen preparation products.** Load both prepared cubes and coordinate grids. Check the science/noise-reference compatibility and spatial sampling; retain Script 01's accepted masks.
3. **Construct continuum S/N maps and the stellar-body aperture.** Measure the configured BL/RH3 reference windows, calculate nominal signal/noise and identify the configured useful galaxy footprint. This separates the analysis aperture from the larger exposed field without discarding every faint member pixel.
4. **Run PowerBin once on BL.** Use the nominal BL continuum capacity and configured target S/N to define bins. The target is nominal; it is not a guarantee of the later empirically calibrated achieved S/N.
5. **Transfer the membership to RH3.** Map native RH3 pixel centers into the BL tessellation through WCS, enforce the permitted matching distance/sampling checks and report coverage. Do not independently tessellate RH3.
6. **Coadd spectra and calibrate spatial noise.** Sum the accepted aperture members and formal variances in each arm. Pass the saved Gaussian-noise stacks through the same masks, memberships and coaddition; compare their empirical variance with the formal variance to obtain each bin's spatial multiplier. Save the formal uncertainty and multiplier separately and apply the multiplier across that bin's native error spectrum.
7. **Assess achieved spectral quality.** Calculate continuum S/N from the calibrated errors, compare the nominal PowerBin and achieved values, and flag poor window coverage, nonpositive continuum and other documented pathologies.
8. **Record light weights and geometry.** Compute continuum-light weights, bin centroids, sky positions, areas and arm-to-arm coverage. Save the actual per-spaxel memberships/coordinates/weights used by later aperture integrations.
9. **Save the binning handoff.** Write `master_bins.fits`, `master_bin_spectra.fits`, `master_bin_table.ecsv`, `master_bin_membership.npz`, both spatial-noise calibration products, diagnostics and `metadata/script02_manifest.json`. Review the aperture, membership transfer and S/N diagnostics before Script 03.

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

## 6.0 Execution walkthrough - implemented

Driver: `scripts/03_build_RH3_likelihood_cubes.py`. Every spatial bin is fitted independently; spatial coherence is not imposed here.

1. **Initialize or validate a resumed run.** Load config and the explicit upstream runs, establish worker/thread settings and verify the checkpoint schema, config and source-run consistency before reusing completed bins.
2. **Read RH3 spectra and empirical resolution.** Load Script 02's binned RH3 spectra, calibrated errors and masks together with the Script 01 LSF/header. Resolve wavelength media, the adopted redshift, fit interval and fixed atmospheric/rest-frame exclusions.
3. **Freeze the log-wavelength experiment.** Rebin spectra and propagate formal errors onto the selected velocity-scale grid. Define a fixed good-pixel set per bin and reject invalid fitted samples. Finite interface placeholders are permitted only for already excluded samples, as described in Section 6.4.
4. **Prepare the full XSL basis and check coverage.** Match the empirical LSF within measured support, normalize each SSP over the fixed fitting band, prepare padded templates at the science velocity scale and run the installed pPXF coverage preflight before starting the expensive bin fits.
5. **Freeze the noise/fitting fingerprint.** Record the spectra, noise, masks, templates, grid and fit settings defining the full-grid residual-calibration experiment. Scan for compatible completed-bin checkpoints.
6. **Evaluate each pending bin.** Fit the one-component control and every configured `V_A x V_B x f_A` state using the Script 02 error spectrum. At a two-component state, keep both velocities and the light fraction fixed while profiling dispersions, SSP weights and the allowed continuum terms. Retain failed-state classes and boundary diagnostics rather than assigning fictitious finite scores.
7. **Calibrate that bin from its best grid state.** Refit the best exact state to retain its model spectrum, measure `q = 1.4826 * MAD[(galaxy-model)/sigma_spatial]` on the fixed good pixels and retain values above or below one without a floor. Set `sigma_fit = q * sigma_spatial` and divide the grid/control chi-square values by `q^2`; the uniform unregularized rescaling does not require repeating the full grid.
8. **Checkpoint and consolidate.** Save each completed bin atomically. Once all bins finish, assemble `RH3_likelihood_cubes.npz`, `RH3_log_spectra_and_local_best_fits.npz`, `XSL_RH3_templates.npz`, `RH3_scalar_noise_calibration.npz`, the local-summary table, diagnostics and manifest. Delete intermediate checkpoints only after successful consolidation when configured to do so. Later stages consume the final calibrated chi-square/noise directly without multiplying either scalar again.

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

The implemented driver is scripts/04_fit_global_two_disk_model.py. See [PIPELINE_CALLS.md](PIPELINE_CALLS.md) for commands and the required human review steps. The design below supersedes earlier plans for this stage. It fits the saved Script-3 spectral likelihoods; it does not call XookSuut, refit pPXF, or produce posterior samples.

## 7.0 Execution walkthrough - implemented

Driver: `scripts/04_fit_global_two_disk_model.py`. Its phases are explicit commands, not one automatic uninterrupted workflow. The sequence is PSF measurement, profile review, smoothing calibration, strength review, eight-model comparison and final model selection.

1. **Supply the reviewed initialization inputs.** Set the MaNGA dispersion peaks/baseline and peak distances, choose the PyMorph table or explicit inclination bounds, check the PA/geometry bounds and standard-star paths, and pass the three actual upstream run directories. Original Windows paths in manifests remain historical provenance when products have moved.
2. **Measure the standard stars (`--phase psf`).** Fit pixel-integrated Gaussian and Moffat profiles to each reduced standard-star continuum image and save `stellar_psf_measurements.json` and the PSF QA sheet. Inspect residual structure, parameters and flags; explicitly choose `--psf-profile gaussian` or `moffat` for subsequent phases.
3. **Load and freeze the model inputs.** Read the final Script 03 likelihood cubes, Script 02 memberships/continuum-window record and Script 01 prepared RH3 image/coordinates. Subtract each bin's independent cube minimum, build the fixed nonnegative RH3 continuum proxy, derive inclination constraints, transfer the selected PSF profile and verify the scientific fingerprint before reusing outputs.
4. **Establish both radial domains.** Build the restricted two-sigma and full-aperture whole-bin selections, radial nodes, geometry/speed bounds and provisional disk-label convention. Save the reference-geometry plots. These plots show the reference geometry, not the optimizer's final answer.
5. **Run the strength pilot (`--phase calibrate`).** Use the full-aperture, PSF-enabled problem for an unpenalized control and eleven positive `C = lambda*A` strengths around the trial baseline. For each strength, explore a population with differential evolution and refine several promising starts with bounded Powell optimization. Select the accepted alternative with the lowest `Qdata + P`; parallel workers handle independent entire fits.
6. **Save and inspect the twelve-fit sheet.** Write fit checkpoints and the two-row by six-column comparison of rotation curves, per-bin data-cost changes, geometry, flags and descriptive curve differences. This does not automatically choose a strength or launch the eight final fits.
7. **Review the strength (`--phase review`).** Inspect the sheet, request a weaker/stronger sweep or explicitly adopt an already fitted coefficient, including zero, and record the reason. Retain the effective `C` rather than one extent-specific `lambda`.
8. **Compare eight configurations (`--phase compare`).** Fit restricted/full extent, PSF off/on and smoothing off/on, reusing compatible full-aperture PSF pilot fits where appropriate. Save each solution's curves, geometry, bin predictions, data/penalty costs and flags. Compare extents on their common inner domain and report outer contributions separately.
9. **Select a model (`--phase choose-model`).** After reviewing the comparisons, name the selected model and record the rationale in `selected_global_model.json`. Its predicted continuous velocities are the handoff for the exact Script 05 refits; coarse fractions and optimizer trajectories are not final error bars.
10. **Resume only compatible work.** `--resume` reuses completed whole-fit checkpoints and reruns unfinished fits from their starts. Without `--resume`, calibration refuses to overwrite an existing completed fit. Matching reruns regenerate figures; changed scientific inputs/search settings require a new output directory.

## 7.1 Coordinates, geometry and permanent disk labels

Both disks share a fitted center, astronomical PA (east of north), inclination, and systemic velocity. Each has its own nonnegative intrinsic radial speed magnitudes. Their rotation signs remain opposite and fixed at every radius. The configurable numerical ceiling is initially 400 km/s per disk; it is not a physical maximum or the same quantity as the saved LOS velocity-grid bound.

The permanent A convention refers to the outer branch on the negative signed major-axis side. Its provisional sign follows the outermost reference bin's single-component velocity and is flagged for review. An explicit reviewed sign can replace it. Local light dominance never reassigns A/B labels.

Existing center, PA and systemic-velocity ranges are hard bounds. Inclination bounds come from a deliberately selected PyMorph disk-axis-ratio table and the existing q/q0 propagation, or explicit reviewed inclination inputs. The catalog table/band is not silently selected. The inherited q/q0 sampler uses Gaussian draws with physical clipping, which is recorded in provenance. These bounds are not new posterior measurements.

The coordinate calculation uses east-positive and north-positive sky offsets: `major = east*sin(PA) + north*cos(PA)` and `minor = east*cos(PA) - north*sin(PA)`, followed by deprojection of the minor coordinate by `cos(inclination)`. PA=0 degrees therefore aligns the major axis north-south, and PA=90 degrees aligns it east-west. Script 04 sky maps now display north up and east left by reversing the horizontal display axis only. This does not negate the model coordinates, rotate the physical PA or alter saved fit results. The IFU footprint orientation is distinct from the galaxy major axis. The geometry preview hides zero-weight continuum cells, including unused rectangular cube padding, rather than drawing that padding as a black footprint. Positive continuum samples and rings retain the same sky coordinates; this display mask does not change model light weights or rotate either layer.

## 7.2 Fixed radial grids and whole-bin domains

The first free node is one adopted RH3 PSF FWHM from the center. Default node spacing is that FWHM, with nominal annulus half-width half the spacing. Set U(0)=0 and interpolate linearly through the free nodes. Hold the last-node speed constant over outer padding; do not create an unmeasured extra node.

The restricted domain starts from the larger of the two user-supplied MaNGA dispersion-peak distances, snapped to the nearest node. Include complete bins touching the outer annulus, then extend nodes outward to cover those complete bins. The full domain includes every upstream bin and rounds node coverage outward. Membership, node dimension and reference geometry are fixed before optimization. Changing trial geometry cannot improve a score by dropping bins. The actual restricted extent may exceed its nominal edge because of a large crossing bin.

MaNGA peak and baseline dispersions give the equal-light initial LOS speed sqrt(sigma_peak^2 - sigma_base^2), or half the inferred disk separation. Correct by sin(i) for an intrinsic-speed seed. Vary starting shapes, amplitudes, A/B asymmetries and geometry; these seeds add no likelihood term.

## 7.3 Measured stellar PSF and fixed light proxy

Fit both Gaussian and Moffat profiles to Feige 24 continuum images in the RH3 continuum interval, integrating over native rectangular KCWI pixels. The user inspects residuals and profiles before choosing a candidate. Rotate the pixel-integrated profiles into the science orientation and average unit-flux kernels. The median larger-axis effective FWHM sets the radial-node scale.

This is an effective standard-star PSF transfer approximation. The standards are not contemporaneous with the science exposures, and the exact science stack/drizzle response is not replayed. Preserve those limitations as flags. Do not interpret this as a precise recovered science-stack PSF.

Use the observed RH3 continuum as the fixed common light proxy I for both disks. Collapse the recorded Script-2 continuum window respecting the prepared masks; set negative continuum weights to zero as in Script 2. Preserve global relative brightness before convolution, rather than convolving independently normalized bin weights. No intrinsic or component-specific light model is fitted.

For each projected velocity field V, predict first moments as:

    I_model = P * I
    V_model = [P * (I V)] / I_model

Sum numerator and denominator through each entire bin aperture. For the no-PSF mode, use direct continuum-weighted aperture means. Save the corresponding second-moment velocity shear as a diagnostic; do not use it to replace pPXF's component dispersions. Report received light originating beyond the final node and flag spatial edges or missing light. No flux is invented beyond the image.

The observed-light substitution follows a precedent in [DiGiorgio Zanger et al. (2024), section III, equations 4–6](https://arxiv.org/html/2407.11908v1#S3). It remains an approximation for these two-component spectra.

## 7.4 Data objective and fraction profiling

For each trial global model, evaluate every saved RH3 fraction plane by bilinear interpolation in V_A and V_B; only then minimize over fraction. Profiling corner minima first could mix incompatible fractions and artificially lower a score. Linear interpolation in fraction cannot create a lower minimum between the saved planes, so exact finer fraction fitting remains Script 05's task.

Subtract each bin's independent minimum once, then sum the bin contributions:

    Q_data = sum_i min_f Delta_chi2_i[V_A,i, V_B,i, f]

No extrapolation or boundary clipping is allowed. A trial lacking a supported state for any selected bin is rejected in full; no bins disappear from its objective. Failed zero-weight interpolation corners do not invalidate an exact valid grid state.

Read the final Script-3 scalar-calibrated chi-square surfaces unchanged. Do not apply either noise multiplier again. Both wavelength independence within bins and independence of bin contributions remain declared working approximations. The fitting objective is not a covariance-corrected likelihood.

## 7.5 Experimental curvature penalty

Following [Sellwood & Zánmar Sánchez (2010), Appendix A](https://arxiv.org/html/0912.5493v2#A1), apply squared second differences to the two independent radial speed curves:

    P = lambda A sum_disks sum_j (U_j - 2 U_(j+1) + U_(j+2))^2
    A = (Rmax / DeltaR)^4 / Vtyp^2
    C = lambda A

The first retained triple is (U1,U2,U3): omit the triple involving the imposed zero-speed origin. This central omission is our explicit choice, not a prescription attributed to the paper. Remaining terms can still influence U1. There are no artificial zero-slope endpoint conditions.

C is fixed across extent comparisons; lambda_extent=C/A_extent. Keep common spacing and fixed Vtyp=200 km/s for reporting. Vtyp does not constrain the fitted speeds. C has units (km/s)^-2 for the implemented velocity units. Outer models can add penalty terms, so their total penalty need not equal the inner penalty.

Optimize Q_data + P, retaining both parts separately. The unpenalized control uses C=0. These are regularized optimization solutions, not evidence values.

## 7.6 Manual strength calibration and eight comparisons

Before the final comparison, fit a full-aperture PSF-weighted control and eleven positive strengths C=0.002*2^j for j=-5 through +5. The literature-motivated starting value is a trial, not a validated strength for our data objective.

The single-page two-row by six-column sheet shows curves, per-bin data-cost changes, geometry and changes, separate scores, flags, and descriptive RMS curve distances. The user can shift the range weaker/stronger or select a tested coefficient, including zero. Custom values must be fitted and reviewed. Record the selection rationale and freeze C before the final comparison.

Run all combinations of restricted/full extent, PSF off/on, and curvature off/on: eight configurations. Matching pilot fits can be reused. Compare extent pairs on common inner bins and radii, report outer-bin costs separately, and retain parameter counts and mean per-bin costs. These means are descriptive; do not call them reduced chi-square or normalize the optimization objective by bin count. No automatic model preference is assigned.

## 7.7 Search, provenance and downstream contract

Use differential-evolution population exploration followed by multiple bounded Powell refinements of promising starts, evaluating the same objective throughout. Completed selected states and accepted local alternatives are shared only between fits with the same scientific fingerprint, radial extent, PSF setting, radial nodes and bin membership. Each candidate is rescored at the receiving fit's own coefficient; its old penalized score is never copied across coefficients. Shared candidates enter the population/local-start pool, and a local refinement cannot discard a better valid starting state. Save alternative local solutions, search histories, convergence/boundary flags, geometry, radial curves, per-bin LOS velocities, coarse fractions and scores. Numerical search budgets are configurable and their exhaustion remains visible. Search trajectories are not posterior draws or rotation-curve error bars.

Restart Powell from the lowest-objective state after the original multi-start search and retain only finite improvements. For each node speed and each modeled bin's line-of-sight velocity, compute its maximum minus minimum over the last `SCRIPT04_POWELL_STABILITY_WINDOW` retained states (default four, counting the original winner initially). Check both disks and stop when the largest node-speed span is at most `SCRIPT04_POWELL_NODE_TOL_KMS` (3 km/s) and the largest projected span is at most `SCRIPT04_POWELL_LOS_TOL_KMS` (0.25 km/s). Begin testing as soon as the window is full. These tolerances were tested by replaying the same full, PSF-convolved smoothed development model from three fresh search seeds and a search initialized with the saved candidates. They measure changes within a search trajectory.

Keep the best finite state if the search reaches `SCRIPT04_POWELL_MAX_RESTARTS` (default 100) or a restart fails to improve the objective. Record `KINEMATIC_RESTART_LIMIT_REACHED` or `KINEMATIC_RESTART_STALLED` respectively when the velocity criterion remains unmet. The original local alternatives remain in the fit record for comparison and candidate sharing.

After each calibration or eight-model batch, cross-check each result against the saved compatible candidate pool, including earlier sweep ranges and archived pre-refinement states. Refit cases beaten by a known candidate, up to `SCRIPT04_SHARED_REFINEMENT_PASSES` passes (default 3). This is a numerical budget, not a physical prior or proof of the global optimum. Save candidate sources, per-fit checks and previous versions under `metadata/shared_start_history/`. If a better known state remains when the pass budget is exhausted, retain `SHARED_START_REFINEMENT_INCOMPLETE` and allow further work on resume. Candidate sets are frozen before dispatch so worker completion order does not change another job's starting pool.

Parallelize independent entire fits, with one numerical-library thread per worker. A completed fit is an atomic checkpoint. Compatible completed fits are reused but still cross-checked when new shared candidates appear. Verify code/input/config fingerprints before reuse. Changed science inputs or search settings require a new output directory; explicit upstream paths support transferred Windows runs without rewriting original provenance.

The user selects one of the eight models explicitly. The selected-model product points to its joint geometry/curves and bin predictions for exact Script-05 pPXF refits. Final uncertainties are a later inference task, conditional on the declared noise and spatial-model assumptions. Successful software tests do not establish that the low-quality development data constrain a physical solution.

The [Script 04 reference notes](references/SCRIPT04_REFERENCES.md) describe the design precedents. Reusable BibTeX is in [crd_dap.bib](references/crd_dap.bib). The ring framework is inspired by [XookSuut](https://arxiv.org/abs/2110.05095); our spectral-cube objective and optimizer are a custom implementation, not that package's Bayesian machinery.

---

## 7.8 Reading Script 04 calibration logs

The main workflow log is `logs/04_fit_global_two_disk_model.log`. Each independent fit has `logs/<fit_name>.progress.log`; for example, `calibration_C0.001.progress.log` follows the full-aperture, PSF-enabled pilot with effective smoothing coefficient `C=0.001`. The coefficient is fixed during that fit; it is not a noise scale, a convergence tolerance or the reported extent-specific `lambda`.

| Field | Definition and units |
|---|---|
| Timestamp | Local wall-clock time when this progress record was written. |
| `step` | Combined progress-history index across differential-evolution callbacks, successive Powell refinements and local-result bookkeeping. It is not a bin number, number of spectral fits, posterior sample, evaluation count or percentage complete. |
| `Qdata` | Sum over selected bins of the profiled/interpolated chi-square excess above each bin's independent Script 03 minimum; see Section 7.4. This is neither raw spectral chi-square nor reduced chi-square and is not expected to approach one. |
| `P` | Current curvature penalty `C * sum(second_difference(U)^2)` for both radial speed curves, excluding the imposed-origin triplet; see Section 7.5. |
| `geometry` | Ordered array `[east0_arcsec, north0_arcsec, PA_deg, inclination_deg, Vsys_kms]`. East/north offsets refer to the prepared coordinate origin, PA is east of north, and systemic velocity is in the original Script 03 residual-velocity frame. |
| `elapsed` | Wall-clock seconds since this entire independent fit began, including setup/search time; it is not exposure time, summed CPU time or elapsed time for the whole twelve-fit calibration. |

For each candidate, the optimizer varies shared geometry and both disks' radial speeds, projects their velocities, applies the adopted continuum/PSF aperture weighting, profiles the saved fraction planes and minimizes `objective = Qdata + P`. It does not rerun pPXF or recalibrate uncertainty amplitudes. The progress line prints only the five geometry parameters; all node speeds still participate in the calculation.

A steep linear rise has zero second difference even if its slope is large. For example, a retained triple `(100,200,300) km/s` adds zero curvature cost, while `(100,200,100) km/s` adds `C*200^2`, or 40 when `C=0.001`. An initial `P=0` can therefore describe a zero-curvature starting curve even during a penalized fit. Conversely, a slightly larger `Qdata` can be accepted if the penalty decreases enough to lower their sum.

The global phase uses `SCRIPT04_GLOBAL_ITERATIONS` as a differential-evolution generation budget. Local refinement uses up to `SCRIPT04_LOCAL_STARTS` distinct starting states, each with its own `SCRIPT04_LOCAL_ITERATIONS` Powell budget. These optimizer operations are not all the same unit. Defaults of 80 global generations and four local starts do not imply that the combined `step` counter stops at 80.

The implementation prints the first valid callback and then every fifth combined history entry. It also appends accepted local-result records to the saved history without necessarily printing a matching progress line. Moving to another local starting state can raise the displayed cost. The fit saves the lowest objective reached.

The Powell restart lines report the restart number, whether it improved the objective, the retained and trial objectives, and the rolling node and projected velocity spans in km/s. Spans are populated once the window is full. The final line reports the stopping reason. Read the saved `fits/<fit_name>.json` for the selected parameters, alternatives, histories, convergence flags, per-bin predictions, separate objective terms, and the `kinematic_refinement` audit.

`GLOBAL_SEARCH_BUDGET_REACHED` means differential evolution did not declare success within its search budget; it is not a program crash or a proof of a global optimum. `LOCAL_OPTIMIZER_NOT_CONVERGED` marks an unsuccessful selected local refinement. `GEOMETRY_BOUNDARY` identifies a fitted geometry parameter near a hard bound; `ROTATION_SPEED_BOUNDARY` and `ROTATION_SPEED_CEILING` do the same for speed limits. Inspect these alongside the saved solutions rather than changing bounds solely to remove flags.

Fit progress logs and the main log append across attempts. Use timestamps and restarted step/elapsed counters to distinguish attempts; old log lines are not evidence that their checkpoint still exists in the current directory. A completed whole-fit JSON checkpoint, with matching fingerprint, is the reusable result.

## 7.9 Log header and search-version restart requirement

**Implemented 2026-10-26.** Every new independent fit attempt starts with field definitions, including the geometry order/units, `objective = Qdata + P` and the meaning of the combined progress counter. It also records the initial population size, shared-candidate count and best initial objective. Existing append-mode log records are not rewritten.

The shared-start refinement, log header, and Powell stability check change hashed optimizer/driver code. Use a new Script 04 output directory to run the updated search on the original Scripts 01–03 products. Run the PSF phase in the new directory, then repeat calibration, coefficient review, comparison and model selection. Subsequent `--resume` calls within that new, matching run reuse completed fits and continue any necessary shared-start refinement.

---

# 8. Script 5 — exact RH3 solution and quality assessment

## 8.0 Execution walkthrough

`scripts/05_extract_global_RH3_solution.py` extracts conditional RH3 solutions at the chosen coherent Script 04 velocities. Its inputs are the selected fit, the verified spatial membership and coordinates, and the prepared Script 03 spectra, template basis and calibrated errors.

1. **Verify inputs.** Require the selected-model record, fit and prediction table to agree. Verify the Script 03 cube, spatial membership and coordinates against the saved Script 04 provenance. Load the saved spectral settings, verify the calibrated noise on fitted pixels and record the input and numerical-environment fingerprints. Use explicit upstream directory overrides after moving a run between machines.
2. **Fit each selected bin.** Evaluate every fraction in the fixed 0.01–0.99 grid. Both component velocities stay at the continuous Script 04 predictions. Each fraction starts from the MaNGA baseline for both dispersions and profiles the two widths independently with the established SSP basis and continuum terms.
3. **Checkpoint each initial scan and review retries.** Each worker completes one bin profile. The main process writes the bin product before displaying pending recovery decisions. Save each retry choice and its outcome atomically.
4. **Select the best available state and generate diagnostics.** Use the minimum total chi-square among scorable results and retain its flags. Save the full fraction profile, per-bin figure, galaxy-wide kinematic maps and signed major-axis velocity profile.
5. **Write the Script 06 handoff.** The ECSV table provides the adopted kinematics and identifies each bin's spectral product. Its original bin IDs include explicit rows for bins outside the selected model. Formal uncertainty propagation is described in Section 10.

## 8.1 Fixed spectral experiment

The saved Script 03 inputs supply the normalized spectra, wavelength grid, good-pixel mask, final scalar-calibrated diagonal errors, LSF-matched templates, velocity sampling, dispersion bounds and continuum settings. Evaluate total chi-square on those same good pixels. The scalar calibrations are already included in the saved errors. Use the prescribed RH3 stellar-template light fraction in its original normalization passband and keep the additive continuum separately recorded.

The 99 fractions are fixed trial coordinates. At every fraction the starting widths are both `MANGA_SIGMA_BASE_KMS` and are then fitted independently. With $V_A,V_B$ fixed, the selected fraction provides its own fitted dispersions and component spectra. The original Script 03 cube remains an upstream comparison product. Exact refitting here refers to direct pPXF evaluation at the continuous coherent velocities.

## 8.2 Numerical recovery and resume

The parent process displays the bin ID, fraction, most recent dispersion pair and reason for recovery. An exception, unusable numerical result or reported non-convergence triggers review. Boundary proximity is recorded as a warning. You can retry the most recent starting pair, enter a new A/B pair, or move on. The new-input prompt suggests the original MaNGA baseline minus and plus 15%, rounded to km/s when within the inherited bounds. The budget is one initial fit plus up to three additional calls at each fraction.

Retain all attempt records. A returned fit is scorable if its outputs are finite and satisfy the fixed coordinates and dispersion bounds. A finite non-converged fit remains available with its optimizer flag. Keep the lowest scorable chi-square at each fraction and then across the bin. An unscorable fraction remains missing. Report a best available minimum from an incomplete profile and leave the solution missing for a bin with no scorable fraction.

Use the same output directory and `--resume` to reuse matching profiles. An unattended run retains unresolved retry decisions and exits with status 2. Resume from an interactive terminal to review those decisions. Interrupted recovery calls retain the reserved attempt and chosen start. Worker or memory interruptions preserve saved profiles, and the worker count can be reduced on resume.

## 8.3 Handoff and diagnostic products

Each `products/RH3_fraction_profiles/bin_NNNN.npz` holds all 99 fractions, best available results per fraction, component and total spectra, template and polynomial coefficients, achieved fractions and attempt history. This file also serves as the restart checkpoint. The saved Script 03 spectra and errors provide the common observed data for residual calculation. Store bin products independently to bound memory usage.

`products/RH3_global_kinematics.ecsv` provides the selected velocities, refined dispersions, fractions, fit scores and flags by original bin ID. It includes the independent Script 03 minimum and the interpolated Script 04 score on the same chi-square scale. A direct finer fit can improve on the coarse-grid minimum. The manifest records the selected model, verified input hashes, fitting settings and relative product paths.

Review the 99-point profile, residuals, optimizer status, dispersion-boundary flags and the fine fraction endpoint flag. The coarse velocity-grid warning uses the configured distance from the original velocity edges. The per-bin PDFs and best-state galaxy-wide figures are described in DIAGNOSTICS. Their velocities are displayed relative to the fitted systemic velocity. The table preserves the original Script 03 frame.

---

# 9. Script 6 — BL two-component stellar populations

## 9.0 Execution walkthrough

The driver `scripts/06_fit_BL_two_component_populations.py` fits the BL spectra conditional on the exact Script 05 stellar kinematics. Each bin retains its original ID and disk labels. The same two velocities and dispersions are fixed at every fraction. The BL fraction and stellar mixtures are determined by the blue data.

1. **Verify the handoff.** Read the completed Script 05 table, the matched Script 02 BL spectra and membership, and the prepared Script 01 BL wavelength grid and empirical LSF. Preserve Script 05 flags and fit only bins in the selected model with usable RH3 kinematics and valid BL inputs.
2. **Prepare the line review.** Sum the eligible aperture spectra on the native observed grid, retaining the original flux scales and showing the number of contributors per channel. Shift the line list to the target redshift and check nominal core coverage per bin. The review PDF displays the proposed gas lines, fitting band, reference band, available pixels and velocity search span. Approve the proposal or edit the saved JSON and regenerate the preview.
3. **Prepare the BL fit.** Rebin the BL spectrum and spatially calibrated error to a common rest-frame log grid. Use the inherited quality and atmospheric masks plus approved exclusions. Match the XSL templates to the measured BL instrumental resolution and prepare gas templates. Preserve the empirical LSF domain and the wavelength support needed by the fixed stellar kinematics. The actual fit range may be narrower than the requested range.
4. **Evaluate the full BL fraction grid.** Fit all 99 fractions from 0.01 to 0.99. Each initial trial uses the same approved pixels, fixed stellar kinematics and MaNGA gas seed. Save each bin atomically as soon as its scan returns. In a terminal, the progress bar refreshes while workers run and clears for recovery prompts.
5. **Review grouped recovery.** A failed or non-converged fraction is recorded with its reason. Retry the pending group using the latest starts, enter a new gas start for that group, or move on while retaining usable states. Each fraction has up to three recovery attempts after its initial fit. An unattended fit saves unresolved decisions and exits with status 2. Resume in a terminal with the same reviewed output directory.
6. **Calibrate and summarize.** Select the best available unregularized fit, apply the BL residual scalar once, and summarize both disks' stellar populations and flags. Write a table, the spectra and noise arrays, the BL template basis, per-bin fit profiles, and diagnostic figures for review and Script 07.

## 9.1 Fitting band, line coverage and approved masks

The requested BL fitting interval defaults to 3700–5350 Å in the XSL rest-frame medium. The fixed 4800–5300 Å reference band defines the light-fraction weights for every bin. It is independent of the wider interval used to fit the populations. Templates are matched to the empirical BL LSF. Stellar wavelength padding follows the fixed Script 05 velocities and dispersions, and separate numerical padding accommodates the wider gas bounds. The optional `--phase psf` measures the BL standard cubes with the same Gaussian/Moffat helper as Script 04 and saves QA for the final science package. The nominal BL population fit uses the Script 05 kinematics directly.

The inherited BL mask comes from Script 02. Additional intervals in the review JSON have a wavelength frame, endpoints, an optional list of original bin IDs and a reason. Observed intervals use the observed science wavelength medium. Rest intervals use the XSL template medium. A proposed gas line is usable when at least half of the nominal core and at least two fitted log pixels are available. The core is one FWHM wide, including the measured BL LSF and the MaNGA seed dispersion. The wider gas center excursion is shown separately. Partial coverage is recorded per bin. The approved pixel mask is fixed across all fractions.

## 9.2 Gas and continuum

The core list is [O II] 3726/3729, Hδ, Hγ, Hβ, [O III] 4959/5007 and [N I] 5198/5200. Higher Balmer lines, [Ne III] and He II are available for review. The [O III] doublet uses the MaNGA DR17 flux ratio and a shared template. Other gas amplitudes are independent, with one gas velocity and dispersion per bin. Templates with usable cores enter the fit and the inherited mask still controls which pixels are scored. Gas templates are matched to the BL LSF.

One target-wide gas seed is taken from the matching local MaNGA DR17 DAPall row. Convert central `HA_Z` to the configured pPXF rest frame and use the already instrument-corrected `HA_GSIGMA_1RE`. The velocity bounds are ±2000 km/s around that seed by default. The dispersion runs from 0.01 BL log pixels to 1000 km/s by default. The bounds and manual gas seeds are configurable. Each trial fits one shared eighth-order multiplicative polynomial for the two stellar components by default through `BL_MDEGREE = 8`. The nominal fit has no additive polynomial and uses `BL_REGUL_LIKELIHOOD = 0`.

## 9.3 BL fraction and noise

Each XSL SSP is normalized to unit mean flux in the fixed reference band before the template library is copied to the two disks. The adopted BL fraction is the sum of Disk A stellar weights divided by the sum of both disks' stellar weights. It excludes gas and is independent of the RH3 light fraction. Save the imposed and returned BL fractions and verify their agreement before accepting a score.

Start with Script 02's BL errors, which already include the spatial scalar. From the minimum unregularized fit measure the BL residual scalar as 1.4826 times the median absolute deviation of error-normalized residuals on the approved pixels. Multiply the spatially calibrated errors by this positive scalar, including values below one, and divide all of the bin's scores by its square. The uniform rescaling preserves the minimum. Save both error vectors and the scalar. An undefined scalar is flagged and leaves the calibrated errors and scores missing. The diagonal noise approximation leaves wavelength correlations outside the adopted model.

## 9.4 Population summaries

The nominal fit uses the configured XSL Salpeter templates for both disks. Retain the full age and metallicity grid. Save both a weighted arithmetic age in Gyr and a weighted mean log10 age in years, plus weighted metallicity on the template library axis. Report these statistics using reference-band light weights and surviving mass weights. Report Disk A minus Disk B for the corresponding summaries.

Undo the SSP light normalization to obtain formed-mass weights. Multiply by the SSP surviving mass fraction, which includes remnants, for the surviving-mass statistics. Compute each disk's stellar mass-to-light ratio in the configurable rest-frame filter, SDSS g by default, using synthetic filter luminosities of the original template mixture in solar units. Save the per-disk age–metallicity light, formed-mass and surviving-mass grids and the age-binned formed-mass history. These grids and histories describe the nominal unregularized solution.

## 9.5 Outputs and quality flags

`products/BL_stellar_populations.ecsv` gives the adopted kinematics, BL fraction and gas fit, calibrated score, population summaries and flags for every original bin. `products/BL_fraction_profiles/bin_NNNN.npz` holds all attempted fractions, fit scores, template and continuum weights, component and total spectra for each accepted state, recovery history and best-state population grids. `products/BL_log_spectra.npz` preserves the approved pixel masks, flux normalization, and both error vectors. `products/XSL_BL_templates.npz` records the matched basis and population metadata. The manifest records inputs, settings and review provenance. See DIAGNOSTICS for the figures.

Preserve Script 05 flags and distinguish missing or incomplete BL profiles, optimizer failures, gas bounds, fraction endpoints, multiple sampled local minima and SSP-grid edge weight. The local minima count describes the profile shape. Store the continuous light and mass shares on the template boundaries so science cuts can be chosen after reviewing the data. Formal uncertainty propagation belongs to Script 07.

---

# 10. Script 7 — uncertainties, systematics, and final maps

## 10.0 Execution walkthrough

Driver `scripts/07_uncertainties_and_final_maps.py` estimates sampling uncertainty under the adopted two-disk model and scalar-calibrated diagonal noise. Use one paired ensemble to measure the joint variation of the rotation curves, per-bin kinematics, and stellar populations.

1. **Load the nominal analysis.** Gather the selected global disk model, the best-fit RH3 and BL spectra, the calibrated errors, and the saved masks and spatial memberships. Define fixed eligible bin sets before drawing noise.
2. **Generate one mock observation.** Add a Gaussian draw to each saved best-fit RH3 and BL spectrum. Use separate noise draws for the two arms and retain the identity of the realization throughout the chain.
3. **Refit the coherent model.** Evaluate the full nominal Script 03 velocity–fraction grid on the noisy RH3 spectra. Use the same grid bounds, spacing, and fraction planes as the real-data fit. The current grid has $17\times17\times9$ states per bin. Refit the selected Script 04 model on these scores.
4. **Refit dispersions and populations.** Pass the recovered coherent velocities and the same noisy RH3 spectra to Script 05. Refit dispersions and RH3 fractions, then pass the recovered velocities and dispersions to Script 06 and fit the noisy BL spectra. Use the fraction search in §10.3 for these repeated fits.
5. **Summarize 100 accepted realizations.** Save the recovered quantities for every draw, the run and bin failure records, and the statistical summaries on a common set of reporting radii and bin IDs.
6. **Evaluate analysis choices and export the science package.** Report the sensitivity checks alongside the nominal results and statistical errors. Include machine-readable measurements, figures, and the provenance of the input products and mock ensemble.

With `--kroupa-library`, Script 07 refits the BL populations using a supplied matched XSL Kroupa library. This separate comparison uses the nominal observed BL spectra, the calibrated errors, and all 99 fractions. Hold the kinematics, masks and other fit settings fixed to isolate the IMF effect. Report the changes in each disk and in their contrasts, including age, metallicity, stellar mass-to-light ratio, disk mass ratio and star formation history.

The science-package index links the original calibration and QA products alongside the exported results. The final science package brings together the bin maps, coherent model, kinematic and population products, PSF and noise QA, numerical tables, figures and provenance. Include an annotated integrated BL plus RH3 spectrum on a common rest-wavelength axis with the gap visible. Restore per-bin flux scales, sum over a matched modeled aperture and preserve the observed, total-fit, Disk A, Disk B and gas flux arrays.

## 10.1 Coherent-model refits

Refit every parameter that was free in the selected nominal Script 04 fit, including the center, PA, inclination, systemic velocity, and radial speeds where applicable, using the same bounds. Hold the bin set, disk labels and rotation signs, radial extent and node locations, PSF and LSF treatment, continuum light proxy, and smoothing strength at their adopted nominal settings. This measures how the fitted quantities vary with noise for a fixed analysis procedure.

Use the Script 04 search protocol and numerical settings adopted for the nominal fit, including the broad starting-model pool and the Powell velocity stability check in §7.7. Score candidates against each mock's own Script 03 cubes. Retain the restart outcome and convergence flags with each realization.

Script 07 fixes the RH3 bin set before drawing noise, using only bins with a usable best-fit coadded model spectrum in the original Script 05 handoff. Refit the global model on this same set for every mock. Use the intersection with the original usable Script 06 bins for BL population refits. Keep the original nominal failures outside the ensemble and count mock failures against the eligible set for each stage. Retain the original Script 04 fit as the nominal curve. If it used additional bins, record the difference in coverage when reporting the ensemble curve errors. On the 8143-1902 development run, omitting the three bins without usable Script 05 fits changed intrinsic node speeds by at most 0.13 km/s in matched Script 04 refits.

Replace a realization if it has no valid coherent fit or if at least 10% of the fixed eligible bins in either downstream stage lack a usable result after recovery. A finite, supported fit that stalled or reached the Powell cap remains usable and flagged. Apply the fraction within the corresponding stage and retain the draw ID, failure reason, and affected bins for every attempt. Below that threshold, retain successful measurements and record missing results for the affected bins. Boundary hits and weak constraints are recorded separately from unusable fits. Log the overall rejection fraction because preferential failure of difficult draws can narrow the reported spread.

## 10.2 Noise realizations and reported errors

Use the saved scalar-calibrated error vectors to generate and score the mocks, with the original spatial and residual factors applied once. Hold the calibration fixed across realizations. Draw independent Gaussian errors under the adopted diagonal-noise approximation. The uncertainty results are conditional on that calibration and noise model.

Keep the nominal measurement as the point estimate and report the sample standard deviation across usable realizations as the statistical error. Save the 16th and 84th percentiles to show asymmetric distributions, along with the number of usable draws and the frequency of boundary hits for each quantity. Compute Disk A minus Disk B and radial summaries within each draw before calculating their spread. The default target is 100 accepted realizations (`SCRIPT07_N_DRAWS`). Each quantity records its own usable-draw count. Rejections do not count toward the target.

## 10.3 Coarse-to-fine light-fraction search for Monte Carlo refits

Evaluate 21 coarse fractions at $f_A=0.01,0.06,\ldots,0.96,0.99$. Use the total chi-square scores under the fixed, calibrated noise. Let $S$ be the maximum coarse chi-square minus the minimum. Reuse the coarse fits and evaluate the remaining points of the full 0.01-spaced grid if any of the following conditions holds.

- **Flat profile.** $S\leq 2$. A change this small across the sampled fraction range gives little guidance about which interval to refine.
- **Departure from a quadratic.** Fit an unweighted least-squares quadratic to the 21 coarse scores as a function of fraction. Trigger a full scan if the largest absolute residual exceeds $\max(0.5,0.06S)$ in total chi-square units. The absolute term filters out small departures and the relative term allows proportionate variation in profiles with a large span. The maximum residual makes the check sensitive to individual abrupt changes.
- **Competing minima.** At least two interior coarse points each have chi-square strictly lower than both immediate coarse neighbors. These minima can compete even when the overall profile is close to quadratic.

If none of the conditions holds, fill in the 0.01-spaced points between the best coarse point and its immediate coarse neighbors, refining only inward at an endpoint. Select the lowest chi-square from the evaluated pPXF fits. A typical interior refinement uses 29 fits in total. A coarse or fine point with no usable score after recovery also triggers a full scan. If some states still have no usable score after a full scan, flag the incomplete profile and retain the lowest finite score. Nominal Scripts 05 and 06 use all 99 fractions.

The cutoffs are empirical search heuristics chosen by replaying the shortcut on the existing test profiles. The 0.05 scan alone missed 5 of 575 complete RH3 profiles and 4 of 578 complete BL profiles. The three triggers together recovered all nine saved full-scan winners and averaged 62.3 RH3 fits and 38.0 BL fits per profile. A 7% relative residual threshold missed BL bin 564 in that replay. The numbers describe agreement with the saved scans on the development data. A one-draw pilot with input kinematics held fixed tested 40 noisy bins per arm against identical 99-state scans. All winners, RH3 dispersions and BL population summaries agreed. The shortcut averaged 81.1 fraction fits per RH3 bin and 48.0 per BL bin.

## 10.4 Parallelism

Process draws in attempt order and parallelize independent bins within the active stage. This gives the terminal one owner for fit-recovery choices and keeps scratch storage bounded by one active draw. Workers each use one numerical-library thread. The user can change `--workers` on resume. Use the Script 05 individual and Script 06 grouped recovery prompts with up to three additional attempts per fraction. With `--non-interactive`, record Move on for unresolved prompts and continue within the usual quality rules.

Restart with `--resume` and the same output directory. The saved draw ID, bin ID and arm regenerate identical noise. Verify the nominal inputs and the run fingerprint, skip committed draws, reuse completed bins or the completed whole-disk fit, and continue from the first unfinished part of the active draw. Pause at the configurable consecutive rejection limit (`SCRIPT07_MAX_CONSECUTIVE_REJECTIONS`, default 10), including in noninteractive mode. Resuming after review acknowledges the pause and preserves the full rejection history.

Store accepted and rejected attempts as compact numerical samples and audit records. After committing a draw, remove its run-owned grids, spectral intermediates and bin checkpoints. Create publication figures only during the final export. Final figures and tables can be recreated with `--resume` after all draws are committed.

## 10.5 Three uncertainty layers

Keep distinct:

### Statistical uncertainty

Gaussian spectral-noise realizations through the full inference chain under the adopted scalar-calibrated diagonal model.

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



and differences from any approved robustness comparisons.

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
- two-component mock recovery, fit failure, and boundary-hit metrics;
- approved BL gas line list, wavelength exclusions and calibrated errors;
- usable-draw counts and sample standard deviations, with 16th and 84th percentiles;
- standardized quality flags.

An Astropy ECSV/FITS table is preferred internally because it preserves metadata and units. A journal machine-readable-table format can be generated at publication time.

## 10.9 Aperture-wide stellar populations and light fractions

Report one row per disk for the fitted KCWI aperture alongside the resolved measurements. Use the adopted nominal footprint and eligible BL bins for both disks, and retain the aperture definition and usable coverage for each draw. The spatial weighting follows [Coccato et al. (2013), §4.3 and Table 2](https://arxiv.org/pdf/1210.7807#page=7), who reported component luminosity-weighted averages and radial profiles from spatially resolved IFU measurements. Script 07 applies that spatial weighting to the full-spectrum pPXF population fits. See [the reference index](references/README.md#crd-decomposition-and-detectability-precedents) and BibTeX key `Coccato2013`.

For each disk and bin, use the disk's stellar template weights in the rest-frame 4800–5300 Å BL reference band and restore the saved per-bin flux scale before summing across bins. The input spectra already coadd the spaxels within each bin, so each bin contributes its recovered total light once. Compute the aperture mean of each disk's age and metallicity by weighting the bin measurements by that disk's restored reference-band light. Report the arithmetic mean age in Gyr and the mean of log age in log years separately, and report metallicity as the weighted mean of the SSP [M/H] values in dex. For the mass-weighted counterparts, use that disk's inferred surviving stellar mass in each bin, restoring the bin scale and applying the Script 06 XSL mass convention.

For the BL aperture light fractions, divide each disk's summed stellar reference-band light by the sum for both disks over the same bins. The fraction uses stellar light only. Label the passband and fitted footprint in the table. Derive the aperture mass-to-light ratio from summed surviving mass divided by summed luminosity in the configured M/L filter, which defaults to SDSS `g`.

Calculate the aperture statistics and Disk A minus Disk B contrasts separately for each realization. Report the nominal result as the point estimate and the across-draw sample standard deviation, 16th and 84th percentiles, and number of usable realizations. Track missing-bin coverage in the draw records so the aperture summary can be interpreted alongside the resolved maps.

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
- reliability thresholds for interpreting weak two-component measurements;
- whether a particular target warrants additional Monte Carlo realizations beyond the initial 100;
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
- two-component recovery and its limitations have been assessed on mocks;
- BL population fits do not unknowingly sit on template-grid boundaries;
- the approved Script 06 BL fraction profiles and population flags have been reviewed;
- end-to-end Monte Carlo error bars and draw failure rates have been reviewed;
- major analysis-choice/systematic tests have been performed;
- all important diagnostics have been inspected and documented.

Only then should final maps, radial trends, and machine-readable measurements be treated as the scientific result.

### Script-3 state-failure bookkeeping

The profile-likelihood surface is constructed only from states for which pPXF returns a valid solution at the requested exact velocities. State-level pPXF exceptions and fixed-velocity mismatches are not assigned finite chi-square; they are recorded as distinct fit-status classes. If an entire bin contains no valid state, Script 3 aborts rather than constructing a fictitious likelihood surface, and reports the underlying pPXF exception frequencies for diagnosis. This bookkeeping is numerical quality control and does not alter the scientific likelihood definition.
