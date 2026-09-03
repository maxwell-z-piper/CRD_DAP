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

Script 1 must therefore characterize:

- residual normalized-width / variance rescaling;
- spectral covariance introduced by interpolation/resampling/stacking;
- whether a full pPXF covariance matrix or a validated approximation is appropriate.

This is a required precondition for trusting profile-likelihood widths.

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

## 5.5 Geometric aperture spectra and formal uncertainties

Spatial extraction is a **geometric sum**, not an inverse-variance-weighted spatial average. This preserves the actual light mixture that later bin-integrated disk models must predict.

For one wavelength sample,

$$
F_i(\lambda)=\sum_{p\in\mathcal B_i}F_p(\lambda).
$$

If the input cube is explicitly in surface-brightness units per square arcsecond, the sum is multiplied by the spatial pixel area so the extracted spectrum represents integrated aperture flux. Formal diagonal variance is propagated as

$$
\mathrm{Var}_i(\lambda)=\sum_{p\in\mathcal B_i}\mathrm{Var}_p(\lambda).
$$

A wavelength sample is accepted for the bin only when at least the configured fraction of its spatial members contributes a valid sample.

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

- `SPATIAL_COVARIANCE_UNCALIBRATED`: the current PowerBin capacity uses formal diagonal spatial variance because no empirical spatial covariance law has yet been validated;
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

## 6.9 Why covariance is calibrated inside Script 3

Script 3 is the first stage that has both ingredients needed for a production spectral-noise model: the **exact log-wavelength spectra actually passed to pPXF** and a flexible stellar model capable of removing the two-component galaxy spectrum. The preliminary Script-1 autocorrelation diagnostic remains useful as an early warning, but it is not used as the final covariance matrix because its residuals are produced by smooth high-pass subtraction rather than by a stellar fit, and because Script 2 spatially coadds spaxels while Script 3 performs the final overlap rebinning onto the pPXF grid.

The production ordering is therefore

```text
Script-2 RH3 PowerBin spectra
        ↓
fixed Script-3 log-wavelength experiment
        ↓
high-quality multi-start two-component pPXF calibration fits
        ↓
residual variance/correlation calibration
        ↓
M1--M4 covariance adequacy + representative full-grid validation
        ↓
freeze one selected covariance model
        ↓
full 17×17×9 profile-likelihood cubes
```

The covariance is calibrated **before** the expensive all-bin likelihood calculation. The complete 2601-state grid is not required to estimate the covariance: covariance calibration needs a sufficiently faithful model spectrum and its residuals, not the full profile-likelihood surface. The complete production grid is evaluated under all four candidate covariance structures only for a small deterministic validation set (Section 6.15), after which the selected covariance model is frozen and used for every production PowerBin.

The covariance calibration is repeated for every new RH3 Script-3 run because it depends on the actual reduction, stacking, PowerBin coaddition, masking, fitting interval, and log-rebinning. Once selected for that fixed experiment, however, it is a **one-and-done statistical calibration**: later RH3 refinement stages should reconstruct and reuse the same frozen model rather than relearn covariance from each refined kinematic state.

## 6.10 Covariance-calibration pPXF residuals

For every PowerBin, Script 3 first performs a one-component control and then a deterministic **multi-start free two-component pPXF fit**. The calibration fit uses the same scientifically relevant ingredients as the later likelihood calculation: the same XSL basis, wavelength grid, good-pixel mask, LSF treatment, polynomial convention, and dispersion bounds. Both stellar components have free velocities and dispersions. Multiple velocity-separation seeds and both A/B orderings are used so that one poor initial guess is unlikely to leave a false two-component residual pattern.

The free two-component calibration fit is not itself interpreted as the physical decomposition. Its component labels may swap, and its optimized velocities need not coincide with the globally coherent solution eventually selected by Scripts 4--5. Its purpose is narrower: obtain a model spectrum that reproduces the observed stellar spectrum well enough that the remaining high-frequency residuals can diagnose the measurement/reduction noise.

For PowerBin $i$ and log-wavelength pixel $j$, define

$$
F_{ij}
$$

as the observed normalized flux density,

$$
M_{ij}
$$

as the best multi-start two-component pPXF calibration model,

$$
r_{ij}=F_{ij}-M_{ij}
$$

as the residual, and

$$
\sigma_{ij}
$$

as the normalized formal one-sigma uncertainty propagated from Script 2 through Script 3's overlap rebinning. The formal normalized residual is

$$
z_{ij}=\frac{r_{ij}}{\sigma_{ij}}.
$$

The per-bin empirical noise-amplitude factor is

$$
s_i = \operatorname{robust\,std}(z_{ij}),
$$

where Script 3 uses the Gaussian-equivalent median absolute deviation (MAD) rather than an ordinary RMS so that a small number of unmodelled pixels cannot set the uncertainty scale for the entire PowerBin. Thus $s_i>1$ means the residual scatter is larger than predicted by the formal errors, while $s_i<1$ means it is smaller. Importantly, $s_i$ is an **overall standard-deviation scale**; it is not a multiplier on the correlation coefficient itself.

A failed one-component control does not automatically prevent covariance calibration if the multi-start two-component calibration fit succeeds. The failure is retained as QC. A PowerBin for which the multi-start two-component calibration fit cannot obtain a valid model is a hard covariance-calibration failure, because its residuals are not trustworthy enough to define the noise model.

## 6.11 Lag correlations, bootstrap uncertainty, and template-mismatch protection

After dividing each residual spectrum by both its formal uncertainty and its fitted scale $s_i$, Script 3 measures the correlation between samples separated by $k$ log-wavelength pixels. For PowerBin $i$,

$$
\rho_i(k)
=
\operatorname{Corr}\!\left(
\frac{r_{ij}}{s_i\sigma_{ij}},
\frac{r_{i,j+k}}{s_i\sigma_{i,j+k}}
\right).
$$

Lag zero is fixed by definition,

$$
\rho_i(0)=1.
$$

Nonzero lags are measured out to the configured generous maximum (`RH3_COVARIANCE_MAX_LAG`; default 20 pixels). The pipeline does **not** assume beforehand that covariance vanishes after a particular number of pixels. Instead, it measures the lag curve and lets the data determine which nonzero lags have statistically supported correlation.

Uncertainty on the pooled lag curve is obtained with a **PowerBin bootstrap**. Entire PowerBins are resampled with replacement; wavelength pixels are never independently resampled, because doing so would destroy the correlation being measured. For each bootstrap realization, Script 3 recomputes the robust pooled lag curve. A simultaneous confidence band is then formed from the bootstrap distribution of the largest absolute departure from the central curve over all inspected nonzero lags (and, where relevant, all wavelength blocks). The production default is a simultaneous 95% band from 2000 bootstrap realizations.

A measured lag is called statistically consistent with zero when

$$
0\in[\rho_{\rm lower}(k),\rho_{\rm upper}(k)]
$$

for that simultaneous band. This is different from saying that the measured coefficient is numerically small: a small coefficient with a very small uncertainty can still be significantly nonzero. Lags whose simultaneous band contains zero are set to zero in the compact adopted correlation model.

The same blockwise calculation is also a safeguard against **template mismatch**. In reality the residual can be written schematically as

$$
r(\lambda)=n(\lambda)+\delta(\lambda),
$$

where $n$ is stochastic measurement/reduction noise and $\delta$ is deterministic spectral-model mismatch. An imperfect line depth, LSF, abundance pattern, or kinematic line shape can produce neighboring residuals with the same sign and therefore mimic statistical covariance. Genuine reconstruction/rebinning covariance is expected to depend primarily on pixel separation, whereas template mismatch is often tied to particular wavelength regions or absorption features. Script 3 therefore tests whitened residual correlations in several broad wavelength blocks even when the candidate covariance itself is wavelength-stationary, and saves residual-stack diagnostics so repeated feature-locked residual structure can be recognized. A model is not accepted merely because positive and negative correlations from different parts of the spectrum cancel in a full-band average.

## 6.12 Definitions of $D_i$, $R_i$, $C_i$, and the four candidate models

For PowerBin $i$, let

$$
D_i=\operatorname{diag}(\sigma_{i1},\sigma_{i2},\ldots,\sigma_{iN})
$$

be the diagonal matrix containing the formal one-sigma uncertainties of the fixed Script-3 wavelength experiment, and let $R_i$ be a correlation matrix with unit diagonal. The covariance model is

$$
\boxed{C_i=s_i^2D_iR_iD_i}.
$$

The four candidates are deliberately nested from least to most flexible:

1. **M1 — diagonal:** $R_i=I$. Each PowerBin has its own $s_i$, but wavelength pixels are otherwise treated as independent.
2. **M2 — common stationary correlation:** all bins share one empirically measured lag correlation $R_{\rm common}$, while every bin retains its own $s_i$.
3. **M3 — common wavelength-block correlation:** all bins share a set of empirically measured lag curves, but the curve is allowed to differ between broad wavelength blocks. Each bin still has its own $s_i$.
4. **M4 — per-bin wavelength-block correlation:** each PowerBin is allowed its own blockwise lag coefficients $R_i$, but only at lags for which the ensemble bootstrap shows globally supported nonzero correlation. This restriction prevents a single noisy residual spectrum from inventing arbitrary long-range covariance.

M2--M4 do not discard the Script-2 variances. They add empirical off-diagonal structure to those formal per-pixel errors. The hierarchy also separates two physically distinct questions: $s_i$ describes **how much** noise remains in a bin, whereas $R_i$ describes **which wavelength samples move together**.

For blockwise models, Script 3 divides the fixed wavelength vector into deterministic equal-width blocks (default three). A pair of pixels is assigned the lag coefficient of the block containing the pair midpoint. This block model is not claimed to be a physical kernel; it is a deliberately simple diagnostic extension used to test whether the stationary approximation is adequate. If M2 passes the residual and likelihood-stability tests, the more complex block structure is not retained.

## 6.13 Exact covariance-aware pPXF likelihood and cached whitening

For a model residual vector $\mathbf r_i$, the Gaussian correlated-noise statistic is

$$
\boxed{\chi_i^2=\mathbf r_i^{\mathsf T}C_i^{-1}\mathbf r_i}.
$$

Factor the covariance as

$$
C_i=L_iL_i^{\mathsf T},
$$

where $L_i$ is the lower-triangular Cholesky factor, and define the whitening operator

$$
W_i=L_i^{-1}.
$$

The whitened residual is

$$
\mathbf q_i=W_i\mathbf r_i,
$$

so the same statistic can be written

$$
\boxed{\chi_i^2=\mathbf q_i^{\mathsf T}\mathbf q_i}.
$$

CRD_DAP uses the validated pPXF 9.4.8 cached-whitener patch in `ppxf_patch_9_4_8/`. The patch adds the `noise_inv_cholesky` interface and has been regression-tested against stock full-covariance pPXF. For a fixed PowerBin, Script 3 constructs $W_i$ **once** and passes the same operator to every pPXF state in that bin. This avoids repeating an identical Cholesky factorization for all 2601 states without changing pPXF's covariance-aware objective function or nuisance-parameter optimization.

The active pPXF version and the presence of the cached-whitener keyword are hard-checked before covariance calibration begins. Script 3 refuses a production covariance run if the validated interface is absent.

Rejected wavelength samples are decorrelated from fitted samples in the full-vector covariance representation. They remain excluded by `goodpixels`, and the whitening operator is explicitly checked so a fitted whitened residual cannot depend on a masked residual through matrix multiplication.

If an empirical lag matrix is not numerically positive definite, Script 3 applies a small eigenvalue floor and renormalizes the correlation matrix to unit diagonal before Cholesky factorization. Any such use is counted and flagged (`RH3_COVARIANCE_PD_REGULARIZATION`) so it cannot occur silently.

## 6.14 Iterative calibration, convergence, and Requirement A

There is a small circular dependence: the initial pPXF residuals are obtained with a provisional covariance model, while changing the covariance can slightly change the best-fitting stellar model and therefore its residuals. Each candidate M1--M4 is therefore calibrated iteratively:

```text
initial multi-start 2-C fit
        ↓
estimate s_i and rho(k)
        ↓
covariance-aware multi-start 2-C refit
        ↓
re-estimate s_i and rho(k)
        ↓
repeat until stable
```

Between iterations $n$ and $n+1$, Script 3 monitors

$$
\Delta_s
=
\max_i
\left|
\frac{s_i^{(n+1)}-s_i^{(n)}}{s_i^{(n)}}
\right|
$$

and

$$
\Delta_\rho
=
\max_k
\left|
\rho^{(n+1)}(k)-\rho^{(n)}(k)
\right|,
$$

with the maximum extended over wavelength blocks and PowerBins where applicable. The locked convergence rule is

$$
\boxed{\Delta_s<0.01\quad\text{and}\quad\Delta_\rho<0.01.}
$$

These are documented numerical development tolerances, not universal physical constants. Their adequacy should be checked in sensitivity tests before the final publication run.

The hard maximum is

```python
RH3_COVARIANCE_MAX_ITER = 5
```

iterations. If a candidate reaches the hard stop without convergence, Script 3 writes the iteration history and a failure record, emits `COVARIANCE_CALIBRATION_MAXITER_REACHED`, and **does not launch the production likelihood grid**. The diagnostic guidance directs the user to inspect calibration-fit quality, coherent stellar-feature residuals/template mismatch, atmospheric residuals, the LSF, wavelength non-stationarity, bin-size/S/N dependence, and multi-start stability. The intended response is investigation, not merely loosening the 0.01 criterion until the run passes.

After convergence, **Requirement A** tests whether the adopted candidate actually explains the residual statistics. Residuals are whitened with the candidate covariance. A passing candidate requires:

- a median robust whitened-residual standard deviation within the configured tolerance of unity (development default ±5%); and
- the simultaneous 95% PowerBin-bootstrap confidence band to contain zero at every tested nonzero lag in every diagnostic wavelength block.

Thus Requirement A asks directly whether the supposedly whitened residuals behave like approximately independent unit-scale noise.

## 6.15 Representative bins and Requirement B: scientific stability of the likelihood surface

Passing residual whiteness is necessary but not sufficient. An unnecessarily flexible covariance model can fit noise in the covariance estimate itself and alter the scientific likelihood surface. Script 3 therefore evaluates **the complete production grid under all four covariance models** for a deterministic representative subset of PowerBins before choosing the production model.

The baseline validation set contains 12 PowerBins: six on each signed side of $PA_{\rm kin}$. On each side, the requested normalized radii are

$$
\frac{1}{7},\frac{2}{7},\ldots,\frac{6}{7}
$$

of that side's usable major-axis radial extent. This keeps the number of validation bins fixed while automatically scaling the physical radial step to galaxy size and sampling inner, intermediate, and outer regions.

The normal major-axis corridor half-width is based on the median equivalent circular PowerBin diameter. If $\widetilde A$ is the median PowerBin area,

$$
\boxed{D_{\rm med}=2\sqrt{\frac{\widetilde A}{\pi}}.}
$$

A bin is normally eligible when its centroid satisfies

$$
|x_\perp|\le D_{\rm med}
$$

for the default corridor factor. Within the corridor, the unused PowerBin centroid closest to the requested point on the PA axis is selected. If no unused bin remains in the corridor, Script 3 still guarantees a representative sample by choosing the unused same-side centroid closest to $PA_{\rm kin}$, using radial proximity as a tie-breaker; the fallback is explicitly flagged in the ECSV and log.

After these symmetric 12 radial bins are fixed, Script 3 searches the preliminary one-component dispersion profile for the strongest off-center local maximum on each side inside a configurable 10--95% radial window. If the corresponding PowerBin is already in the 12, that row is marked as containing the candidate $2\sigma$ region. If it is not, the bin is **added rather than substituted**, so the final validation set contains 12--14 unique PowerBins. This retains symmetric radial coverage while guaranteeing representation of the region expected to provide particularly strong two-component information. The automatic peak is a validation locator, not by itself a claim that the galaxy is physically a classical $2\sigma$ system.

For every representative bin and every M1--M4 candidate, Script 3 evaluates the exact production

$$
17\times17\times9=2601
$$

state grid. No narrowed $\Delta V$ stencil is substituted for this validation.

**Requirement B** compares the scientific outputs of candidate covariance models at the resolution of the production grid. It compares the location of the minimum and the edges of one-dimensional profile-$\Delta\chi^2$ intervals at configured levels (development defaults $\Delta\chi^2=1$ and 4). The default numerical stability allowance is one production-grid cell on any axis. These are grid-resolution/sensitivity settings rather than physical priors and should be checked before publication.

The final decision rule is:

> Choose the **least complex** M1--M4 model that passes Requirement A and whose representative full-grid likelihood surfaces are stable against every more complex model that also passes Requirement A.

Examples:

```text
M1 fails A
M2 passes A and agrees with M3/M4
M3 passes A
M4 passes A
→ adopt M2
```

or

```text
M1/M2 fail A
M3 passes A and agrees with M4
M4 passes A
→ adopt M3
```

If residual-adequate models materially disagree, the pipeline does **not** automatically choose M4. It writes a model-selection failure with diagnostic guidance and blocks the production grid, because disagreement can indicate covariance overfitting, wavelength non-stationarity, template/LSF mismatch, or insufficient information.

Once selected, the covariance model is frozen. PowerBin $i$ may have its own $s_i$ and, if required by the selected model, its own $R_i$, but that same $C_i$ is used for every $(V_A,V_B,f_A)$ state within that bin. This fixed metric is required for the state's $\chi^2$ values to be directly comparable.

## 6.16 One-component control and interpretation of the calibrated likelihood cube

Each production PowerBin still receives a single stellar LOSVD control using the same final wavelength experiment and the **selected frozen covariance metric**. This saves

$$
V_{\star,1C},\qquad\sigma_{\star,1C},\qquad\chi^2_{1C}.
$$

The local raw comparison

$$
T_i
=
\chi^2_{1C,i}
-
\min_{V_A,V_B,f_A}\chi^2_{2C,i}
$$

is useful QC but is not converted directly into a textbook p-value. Secure two-component recovery is calibrated later with mocks because mixture-model regularity assumptions are not guaranteed here.

For each candidate state $c$ in one PowerBin, Script 3 saves the covariance-aware total $\chi^2_i(c)$. The local relative surface is

$$
\Delta\chi_i^2(c)
=
\chi_i^2(c)-\min_{c'}\chi_i^2(c'),
$$

and the corresponding relative-likelihood weight is

$$
w_i(c)\propto\exp[-\Delta\chi_i^2(c)/2].
$$

These are **relative-likelihood weights**, not Bayesian posterior probabilities. After the covariance-calibration requirements pass, the width of this Script-3 surface uses the empirically selected frozen spectral covariance model and is no longer marked as the old development-only diagonal-noise likelihood width.

## 6.17 Saved covariance and likelihood products

In addition to the original likelihood products, Script 3 saves the complete covariance decision trail. Important products are:

```text
products/RH3_covariance_candidates.npz
products/RH3_covariance_calibration_fits.npz
products/RH3_covariance_iteration_history.ecsv
products/covariance_validation_bins.ecsv
products/RH3_covariance_model_validation_grids.npz
products/RH3_covariance_model_comparison.ecsv
metadata/RH3_covariance_model_selection.json
```

`RH3_covariance_candidates.npz` stores the compact $s_i$, lag-correlation coefficients, wavelength-block labels, convergence/QC state, and selected-model metadata. Dense $N_\lambda\times N_\lambda$ inverse-Cholesky matrices are intentionally **not** written for every bin because they are large and can be reconstructed exactly from the compact model plus the saved formal uncertainty vector and good-pixel mask. Downstream RH3 refinement should use this compact saved model to reconstruct/factor each $C_i$ once and reuse it.

The main products remain:

```text
products/RH3_likelihood_cubes.npz
products/RH3_log_spectra_and_local_best_fits.npz
products/RH3_local_likelihood_summary.ecsv
products/XSL_RH3_templates.npz
metadata/script03_manifest.json
```

The likelihood product records the selected covariance-model name/hash and per-bin scale. `RH3_log_spectra_and_local_best_fits.npz` deliberately preserves `noise` as the normalized **formal** Script-2 uncertainty. An additional `ppxf_noise` vector is saved for audit. This distinction prevents a downstream stage from applying $s_i$ twice when reconstructing the frozen covariance.

## 6.18 Restartability and CPU workers

One completed production PowerBin is the atomic likelihood-grid checkpoint:

```text
checkpoints/bin_0000.npz
checkpoints/bin_0001.npz
...
```

Checkpoint schema version and the selected covariance-model hash are checked on resume. A checkpoint generated under a different covariance decision is therefore not silently reused.

The pre-production covariance products are also reusable on resume, but only when the complete saved model-selection decision exists. If per-bin likelihood checkpoints exist while covariance-calibration products are missing or incomplete, Script 3 refuses to recalibrate around those existing states because doing so could mix incompatible $\chi^2$ metrics in one final cube.

The config file is SHA-256 checked before resuming. Script 3 normally obtains the Script-1 provenance path from the Script-2 manifest, but also accepts `--script1-run <path>` for a moved/renamed run. Once completed likelihood checkpoints exist, changing Script-1 provenance is refused because it could mix different wavelength/LSF assumptions.

After every bin has completed and **all consolidated numerical products plus the Script-3 manifest have been written successfully**, the checkpoint directory is deleted by default:

```python
SCRIPT03_DELETE_CHECKPOINTS_ON_SUCCESS = True
```

`--workers N` launches $N$ Python worker processes and caps BLAS/OpenMP thread pools to one thread per worker before NumPy/pPXF import. During both calibration and the expensive production grid, permanent logs record bin-level completion and failures. The production grid retains the terminal-only heartbeat behavior described previously: a spinner, PowerBin checkpoint progress, elapsed time, time since the newest completed bin, and an ETA after timing information is available.

CAPFIT can emit repeated `RuntimeWarning`s when its internal predicted reduction is zero (`actred/prered`) even when pPXF ultimately returns a valid fit. Script 3 suppresses only the two exact CAPFIT scalar-division warnings by warning category, message, and module. Pipeline warnings and unrelated Python warnings remain enabled.

---

# 7. Script 4 — global two-disk model

## 7.1 Radial-node philosophy

Adopt an XookSuut-style non-parametric concentric-ring parameterization rather than an in-house analytic rotation curve.

Default radial choices:

$$
R_{\mathrm{start}}=1\times\mathrm{FWHM}_{\mathrm{PSF,RH3}},
$$

$$
\mathrm{ring\_space}\approx\mathrm{FWHM}_{\mathrm{PSF,RH3}},
$$

$$
\delta=0.5\,\mathrm{ring\_space}.
$$

The first free node is therefore one PSF FWHM from the center. Interior to the first node,

$$
V_{\mathrm{rot}}(0)=0
$$

and the rotation curve is linearly connected to

$$
V_{\mathrm{rot}}(R_{\mathrm{start}})=V_1.
$$

The nominal annuli with full width $2\delta$ do not overlap when $\delta=0.5\,\mathrm{ring\_space}$, but they should not be described as statistically independent because the PSF correlates neighboring spatial information.

## 7.2 Two outer-radius analyses

Always support two radial extents.

### A. $2\sigma$-limited model

Measure the positive- and negative-side $2\sigma$ peak radii from the preliminary RH3 single-component dispersion map, using the MaNGA values only as initial guidance.

Adopt

$$
R_{2\sigma, \mathrm{adopted}} = \max(R_{2\sigma,+}, R_{2\sigma,-}).
$$

Snap it to the nearest allowed ring center:

$$
R_{\mathrm{final}} = R_{\mathrm{start}} + \mathrm{round} \left[ \frac{R_{2\sigma,\mathrm{adopted}}-R_{\mathrm{start}}} {\mathrm{ring\_space}} \right] \mathrm{ring\_space}.
$$

### B. Full-aperture model

Extend the radial grid sufficiently far to include the maximum deprojected radius reached by valid KCWI PowerBin member spaxels, snapping consistently to the ring grid.

The purpose is to test whether the lower-minority-fraction outer data remain consistent with the inner two-disk solution.

### Interpretation

Compare both fits over their common inner radial domain. If the inner rotation curves and global geometry are stable, the full-aperture solution can become the preferred result because it uses all available data. If the solutions disagree materially, the $2\sigma$-limited model becomes the safer fiducial solution and the full-aperture run is treated as a robustness test.

The code should not silently decide based on an undocumented threshold. It should save both solutions and quantify their differences.

## 7.3 Global parameter vector

Baseline:

$$
\Theta=
\{x_0,y_0,PA,i,V_{\mathrm{sys}},
V_{A,1}\ldots V_{A,K},
V_{B,1}\ldots V_{B,K}\}.
$$

Both disks share center, PA, inclination, and systemic velocity initially. They have independent non-parametric ring velocities.

No explicit smoothness penalty is imposed on the ring velocities in the baseline model.

## 7.4 Disk A/B identity convention

Avoid “dominant disk” language because dominance can change with radius and wavelength.

Adopt a reproducible signed major-axis coordinate. Orient the 180°-degenerate PA axis so that increasing signed coordinate follows a documented image-coordinate convention. Walking across the IFU from the negative to positive side:

- Disk A is the outer-branch velocity component encountered first on the reference side;
- Disk B is its counter-rotating counterpart, commonly associated with the inner branch.

If both branches are present at the outermost reference position, use the branch whose velocity sign matches the preliminary outer single-component RH3 velocity. If the assignment remains ambiguous, flag it for human inspection rather than silently guessing.

This convention must remain fixed across RH3, BL, maps, tables, and reruns.

## 7.5 Initial inclination and geometry constraints

Use PyMorph disk axis ratio

$$
q=b/a
$$

and its uncertainty to derive an inclination distribution through

$$
\cos^2 i=\frac{q^2-q_0^2}{1-q_0^2},
$$

where $q_0$ is an assumed intrinsic disk thickness with its own uncertainty.

PyMorph does not directly provide the final inclination uncertainty; that uncertainty is propagated from $q$, $q_0$, and their adopted distributions.

Useful PyMorph quantities include disk axis ratio, disk PA, fitted disk center, and their fit errors. The photometric PA is mainly a cross-check because the stellar kinematic PA from `fit_kinematic_pa` is the primary orientation initializer.

The global parameters have sensible bounds/priors:

- $PA$ informed by the `fit_kinematic_pa` 1σ error;
- $i$ informed by the PyMorph-derived distribution;
- $x_0,y_0$ constrained around the KCWI continuum/photometric center;
- $V_{\mathrm{sys}}$ free within a generous range around the preliminary RH3 estimate.

`fit_kinematic_pa` need not provide a formal $V_{\mathrm{sys}}$ error; Script 7 will determine the uncertainty through the global-model realization ensemble.

## 7.6 Model projection

For a sky position, transform to disk-plane radius $R$ and azimuth $\theta$, linearly interpolate the non-parametric rotation curve, and project

$$
V_{\mathrm{LOS}}=V_{\mathrm{sys}}+V_{\mathrm{rot}}(R)\sin i\cos\theta.
$$

Disk A and B use their own $V_{\mathrm{rot}}(R)$ curves with opposite rotation senses.

## 7.7 Bin-integrated model velocities

Do not blindly evaluate only at the PowerBin centroid.

For member spaxels $p$ in bin $i$, compute a flux-weighted mean model velocity

$$
\bar V_i=\frac{\sum_p I_pV_p}{\sum_p I_p},
$$

and model-predicted unresolved shear

$$
\sigma_{\mathrm{shear},i}^2=\frac{\sum_p I_p(V_p-\bar V_i)^2}{\sum_p I_p}.
$$

PSF effects should be included in the spatial weighting when required. The centroid-only approximation may still be used if diagnostics demonstrate that its difference from the bin-integrated prediction is negligible.

## 7.8 Global objective

For trial $\Theta$, predict bin-integrated velocities

$$
V_{A,i}(\Theta),\qquad V_{B,i}(\Theta).
$$

Evaluate/interpolate the RH3 profile-likelihood cube at those velocities and profile over the local RH3 fraction:

$$
q_i(\Theta)=\min_{f_{A,i}}\Delta\chi_i^2\left[V_{A,i}(\Theta),V_{B,i}(\Theta),f_{A,i}\right].
$$

Then

$$
Q(\Theta)=\sum_i q_i(\Theta).
$$

Minimize or sample the global objective while preserving the parameter history.

## 7.9 Geometry robustness

The baseline shared geometry is intentionally conservative. Later robustness tests should allow limited alternatives such as

$$
PA_A\ne PA_B
$$

if scientifically justified, without making the baseline model unnecessarily underconstrained.

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

Use the covariance-aware noise model characterized in Script 1 rather than independent per-pixel perturbations if the data show significant spectral covariance.

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
- noise scale/covariance is validated;
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
