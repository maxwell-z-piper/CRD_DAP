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

- stacked BL science cube;
- stacked RH3 science cube;
- BL master arc;
- RH3 master arc;
- configuration metadata.

The science cubes are assumed to have already passed through the KCWI DRP. Script 1 is science preparation, not detector-level reduction.

## 4.2 Internal cube representation

The pipeline should expose, for each arm:

$$
F(x,y,\lambda),
$$

an uncertainty or variance cube,

$$
\sigma_F(x,y,\lambda),
$$

quality/mask information, wavelength coordinates, spatial WCS, spectral resolution metadata, and relevant reduction-frame metadata.

Downstream scripts should not need to know raw FITS extension names.

For production analysis, a stacked science cube must preserve a propagated `UNCERT` product. `MASK`, `FLAGS`, and `NOSKYSUB` should also be retained whenever available. A hand-stacked flux-only cube may still be useful for visual or structural testing, but it is not a valid input to the profile-likelihood analysis because the absolute $\chi^2$ scale would be undefined.

## 4.3 Hard bad-pixel / bad-spaxel treatment

A spectral sample is unusable when, for example:

- the DRP quality information marks it invalid;
- the flux is non-finite;
- the uncertainty is non-finite or non-positive;
- it lies outside the instrumentally valid wavelength range.

The exact DRP flag interpretation will be implemented after inspecting real products.

A spatial spaxel may be excluded if the fraction of usable wavelength samples falls below the configurable threshold. **Low S/N by itself is not a reason to reject a spaxel**; Script 2 uses PowerBin to combine low-S/N spatial information.

## 4.4 Wavelength masks

Three concepts must remain distinct:

1. **instrument-good wavelengths** — valid according to calibration/data quality;
2. **RH3 fit-good wavelengths** — wavelengths used for RH3 stellar fitting;
3. **BL fit-good wavelengths** — wavelengths used for BL stellar/gas fitting.

Scientifically valid wavelengths should not be permanently discarded in Script 1 merely because a later fit masks them.

## 4.5 Wavelength conventions

Air/vacuum wavelength conventions and the adopted heliocentric/barycentric frame must be known explicitly. The preferred Script-1 configuration is `auto`, meaning that CRD_DAP reads the convention from explicit KCWI DRP metadata and hard-fails if the header is ambiguous. In the real DRP products used to validate Script 1, the wavelength medium is encoded by the spectral `CTYPE`/comment and the applied radial-velocity correction is recorded by `VCORRTYP`. An explicit config value is still allowed, but it is cross-checked against the header rather than blindly trusted.

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

The output must allow Script 2 to apply one physical BL-defined PowerBin membership to the RH3 cube.

## 4.9 PSF characterization

Measure or adopt separate effective delivered PSFs:

$$
\mathrm{PSF}_{\mathrm{BL}},\qquad \mathrm{PSF}_{\mathrm{RH3}}.
$$

The arms use the same slicer geometry but may not have identical effective PSFs. This matters because using the same geometric PowerBin does not guarantee identical spatial weighting of stellar light in both arms.

The RH3 PSF will also define the default XookSuut-style radial node spacing.

## 4.10 Empirical LSF from master arcs

The primary LSF comes from unresolved lines in the required master-arc products. The master arc is used together with its DRP `wavemap`, `slicemap`, and `posmap` geometry products. Script 1 first attempts conventional filename-root discovery, but also supports FITS-header provenance matching because real RED reductions can assign the geometry maps a different exposure number from the `*_marc.fits` filename. Explicit sidecar paths remain available in the target configuration.

Before any LSF is accepted, the master arc must match the science cube in camera, arm-specific grating, slicer/IFU, detector binning, and central wavelength. This prevents a different RED grating calibration (for example RL) from being used accidentally for RH3 merely because both are red-side products.

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

## 5.1 BL alone defines the spatial tessellation

PowerBin is run once on BL at the target S/N, initially expected to be around

$$
(S/N)_{\mathrm{BL}}\approx40.
$$

The resulting membership map is the master spatial grid for both arms.

RH3 is **not** independently PowerBinned in the production analysis.

## 5.2 Apply identical physical memberships to RH3

For bin $i$ with member-spaxel set $\mathcal B_i$, extract

$$
F_i^{\mathrm{BL}}(\lambda)
$$

and

$$
F_i^{\mathrm{RH3}}(\lambda)
$$

from the same physical sky region after registration.

The achieved RH3 S/N is measured but does not drive the bin boundaries.

## 5.3 Per-bin products

Store at least:

- bin ID;
- member spaxels;
- centroid;
- physical area;
- BL S/N;
- RH3 S/N;
- coadded spectra and uncertainties/covariances;
- flux weights needed for bin-integrated disk-model calculations.

## 5.4 Shared S/N plot normalization

BL and RH3 S/N maps use the same color normalization. The default upper limit is the 95th percentile of the **combined** BL+RH3 bin-S/N distribution, so a given color means the same S/N in both panels while a few extremely high-S/N bins do not destroy contrast elsewhere.

---

# 6. Script 3 — RH3 profile-likelihood cubes

## 6.1 Main product

For every PowerBin construct

$$
\chi_i^2(V_A,V_B,f_{A,\mathrm{RH3}}).
$$

The explicit grid is initially planned at approximately

$$
17\times17\times9,
$$

but the final grid resolution/extent must pass convergence tests.

## 6.2 Meaning of a grid cell

At a fixed grid coordinate

$$
(V_A,V_B,f_A),
$$

those three parameters define the profile coordinate. pPXF optimizes nuisance quantities such as

$$
\sigma_A,\quad\sigma_B,
$$

stellar-template mixtures, and continuum terms.

The result is therefore a **profile likelihood**, not a fully marginalized Bayesian posterior.

## 6.3 Full XSL SSP library

Use the full selected XSL SSP grid for RH3 rather than reducing the template basis purely to save computation time. Scientific completeness is prioritized over runtime.

The RH3-prepared template matrix is produced by:

1. loading the selected XSL SSP library;
2. retaining the needed age/metallicity metadata;
3. ensuring wavelength coverage safely brackets the observed/rest-frame fitting region and allowed velocity shifts;
4. converting wavelength medium if needed;
5. matching the measured RH3 LSF;
6. applying a common normalization convention;
7. log-rebinning consistently for pPXF.

## 6.4 No regularization in likelihood fits

Any fit whose $\chi^2$ enters the RH3 likelihood cube must use

$$
\mathrm{regul}=0.
$$

Regularization may later be used to visualize smooth star-formation histories, but not to define the likelihood surfaces.

## 6.5 Saved nuisance information

Storage minimization is not a scientific goal. Save all scalar products that materially help later interpretation, such as:

- $\chi^2$;
- reduced $\chi^2$ where meaningful;
- $\sigma_A$;
- $\sigma_B$;
- fit-status flags;
- relevant continuum/quality summaries;
- enough template-weight information to reconstruct scientifically important selected states.

Full model spectra/template weights may be saved broadly if storage permits; at minimum they must be saved for scientifically selected and diagnostic solutions.

## 6.6 One-component control

Also fit a single stellar LOSVD per PowerBin to obtain

$$
\chi^2_{1\mathrm{comp},i}.
$$

The one-vs-two-component statistic can be defined as

$$
T_i = \chi^2_{1\mathrm{comp}, i} - \chi^2_(2\mathrm{comp},i}.
$$


A textbook $\chi^2$ p-value is **not** assumed because mixture-model regularity conditions may fail. Significance/reliability is calibrated from mocks in Script 5/model-selection utilities.

## 6.7 Preliminary single-component RH3 maps

Before fixing the final Script-4 radial domain, construct simple single-component RH3 maps

$$
V_{\star,\mathrm{single}}(x,y),
\qquad
\sigma_{\star,\mathrm{single}}(x,y).
$$

These provide a higher-resolution KCWI measurement of the $2\sigma$ structure than the initial MaNGA map.

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
