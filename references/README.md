# CRD pipeline references

Reviewed 2026-10-24 against the five exported CRD Decomposition chats, the current pipeline documentation, and subsequent design discussions. This is a curated reference library for the pipeline, its scientific precedents and alternatives discussed along the way.

BibTeX entries are in [script04.bib](script04.bib). The filename and three original citation keys are retained for compatibility; the bibliography now covers the pipeline as a whole. [SCRIPT04_REFERENCES.md](SCRIPT04_REFERENCES.md) keeps the more detailed Script 04 design notes.

The notes below distinguish implemented methods, scientific context, future-stage precedents and methods we set aside. A citation does not mean that we adopted every part of that paper's analysis. Sources with verified journal metadata use their published records; other entries identify the arXiv version or software repository actually checked. Record software versions, commits and local patches with the eventual science run.

## Spectral fitting, binning and photometric inputs

- **Cappellari & Emsellem (2004)** — [Parametric Recovery of Line-of-Sight Velocity Distributions from Absorption-Line Spectra of Galaxies via Penalized Likelihood](https://doi.org/10.1086/381875). Publications of the Astronomical Society of the Pacific, 116, 138–147. BibTeX: `CappellariEmsellem2004`. Foundational pPXF reference for the pipeline's spectral fitting.

- **Cappellari (2017)** — [Improving the full spectrum fitting method: accurate convolution with Gauss-Hermite functions](https://doi.org/10.1093/mnras/stw3020). Monthly Notices of the Royal Astronomical Society, 466, 798–811. BibTeX: `Cappellari2017`. Primary pPXF methodology, convolution and multiple-component reference. Our current fits retain the declared scalar-diagonal noise approximation.

- **Cappellari (2023)** — [Full spectrum fitting with photometry in pPXF: stellar population versus dynamical masses, non-parametric star formation history and metallicity for 3200 LEGA-C galaxies at redshift z ∼ 0.8](https://doi.org/10.1093/mnras/stad2597). Monthly Notices of the Royal Astronomical Society, 526, 3273–3300. BibTeX: `Cappellari2023`. Modern pPXF optimization/population reference; discussed for later template-systematic checks. The complete LEGA-C/photometric analysis is not adopted.

- **Cappellari (2025)** — [PowerBin: fast adaptive data binning with Centroidal Power Diagrams](https://doi.org/10.1093/mnras/staf1726). Monthly Notices of the Royal Astronomical Society, 544, 1432–1446. BibTeX: `Cappellari2025PowerBin`. Actual Script-02 adaptive binning implementation; the production binning citation.

- **Verro et al. (2022)** — [The X-shooter Spectral Library (XSL): Data Release 3](https://doi.org/10.1051/0004-6361/202142388). Astronomy & Astrophysics, 660, A34. BibTeX: `Verro2022XSLDR3`. Empirical stellar-library provenance underlying XSL; distinct from its SSP-model paper.

- **Verro et al. (2022)** — [Modelling simple stellar populations in the near-ultraviolet to near-infrared with the X-shooter Spectral Library (XSL)](https://doi.org/10.1051/0004-6361/202142387). Astronomy & Astrophysics, 661, A50. BibTeX: `Verro2022XSLSSP`. SSP-model reference for the configured XSL basis; cite alongside DR3.

- **Fischer et al. (2019)** — [SDSS-IV MaNGA PyMorph Photometric and Deep Learning Morphological Catalogues and implications for bulge properties and stellar angular momentum](https://doi.org/10.1093/mnras/sty3135). Monthly Notices of the Royal Astronomical Society, 483, 2057–2077. BibTeX: `Fischer2019PyMorph`. PyMorph methodology and model-preference flags; axis ratios are inputs, while our code derives inclination.

- **Domínguez Sánchez et al. (2022)** — [SDSS-IV DR17: final release of MaNGA PyMorph photometric and deep-learning morphological catalogues](https://doi.org/10.1093/mnras/stab3089). Monthly Notices of the Royal Astronomical Society, 509, 4024–4036. BibTeX: `DominguezSanchez2022PyMorph`. Actual DR17 catalog provenance, band/HDU ordering and flags. It does not specify our q0 assumptions or inclination hard bounds.

- **Meert et al. (2016)** — [A catalogue of two-dimensional photometric decompositions in the SDSS-DR7 spectroscopic main galaxy sample: extension to g and i bands](https://doi.org/10.1093/mnras/stv2475). Monthly Notices of the Royal Astronomical Society, 455, 2440–2452. BibTeX: `Meert2016Multiband`. Recent r-band discussion and wavelength-dependent bulge/disc structure; a diagnostic reference, not an additional pipeline fit.


## Coherent disk geometry, PSF and smoothing

- **López-Cobá et al. (2024)** — [XookSuut: A Bayesian Tool for Modeling Circular and Non-Circular Flows on 2D Velocity Maps](https://doi.org/10.22201/ia.01851101p.2024.60.01.03). Revista Mexicana de Astronomía y Astrofísica, 60, 19–39. BibTeX: `LopezCoba2023XookSuut`. Published radial-ring framework precedent. Script 04 uses a custom spectral-likelihood objective/optimizer, not the XookSuut sampler. The old citation key is retained; publication year is 2024.

- **Sellwood & Zánmar Sánchez (2010)** — [Quantifying non-circular streaming motions in disc galaxies](https://doi.org/10.1111/j.1365-2966.2010.16430.x). Monthly Notices of the Royal Astronomical Society, 404, 1733–1744. BibTeX: `SellwoodZanmarSanchez2010`. Adopted second-difference penalty and A normalization, Appendix A. Section 4.4 footnote 6 supplies Vtyp=200 km/s and lambda=0.002/A. Our central omission and manual fixed-C selection are project choices.

- **DiGiorgio Zanger et al. (2024)** — [The Strength of Bisymmetric Modes in SDSS-IV/MaNGA Barred Galaxy Kinematics](https://arxiv.org/abs/2407.11908). arXiv:2407.11908. BibTeX: `DiGiorgioZanger2024`. Nirvana: section III equations 4–6 and following paragraph motivate intensity-weighted PSF convolution and the observed-brightness proxy. Cited arXiv version; its complete barred-disk analysis is not adopted.

- **Leclercq et al. (2024)** — [Linking Mg II and [O II] spatial distribution to ionizing photon escape in confirmed LyC leakers and non-leakers](https://doi.org/10.1051/0004-6361/202449362). Astronomy & Astrophysics, 687, A73. BibTeX: `Leclercq2024PSF`. KCWI Gaussian/Moffat standard-star PSF precedent; Appendix B in publication, Appendix A in preprint. Does not validate our time-separated standard-to-science drizzle transfer.

- **Krajnović et al. (2006)** — [Kinemetry: a generalization of photometry to the higher moments of the line-of-sight velocity distribution](https://doi.org/10.1111/j.1365-2966.2005.09902.x). Monthly Notices of the Royal Astronomical Society, 366, 787–802. BibTeX: `Krajnovic2006Kinemetry`. Appendix-C / fit_kinematic_pa reference previously discussed for estimating the supplied PA. Script 04 consumes a reviewed PA/error; full kinemetry is not implemented.


## CRD decomposition and detectability precedents

- **Bevacqua et al. (2022)** — [SDSS-IV MaNGA: integral-field kinematics and stellar population of a sample of galaxies with counter-rotating stellar discs selected from about 4000 galaxies](https://doi.org/10.1093/mnras/stab3732). Monthly Notices of the Royal Astronomical Society, 511, 139–157. BibTeX: `Bevacqua2022`. Central CRD identification and two-component velocity chi-square-map precedent; our exact velocity/fraction grid and global selection are extensions.

- **Mitzkus et al. (2017)** — [Dominant dark matter and a counter-rotating disc: MUSE view of the low-luminosity S0 galaxy NGC 5102](https://doi.org/10.1093/mnras/stw2677). Monthly Notices of the Royal Astronomical Society, 464, 4789–4806. BibTeX: `Mitzkus2017`. Two-component pPXF, CaT and velocity-grid precedent. Our production grid fixes velocities exactly rather than permitting within-cell motion.

- **Coccato et al. (2011)** — [Dating the formation of the counter-rotating stellar disc in the spiral galaxy NGC 5719 by disentangling its stellar populations](https://doi.org/10.1111/j.1745-3933.2011.01016.x). Monthly Notices of the Royal Astronomical Society: Letters, 412, L113–L117. BibTeX: `Coccato2011`. IFU component-spectrum/population decomposition and synthetic-recovery precedent. Does not establish our noise model or imply later uncertainty tests are complete.

- **Coccato et al. (2013)** — [Spectroscopic evidence of distinct stellar populations in the counter-rotating stellar disks of NGC 3593 and NGC 4550](https://doi.org/10.1051/0004-6361/201220460). Astronomy & Astrophysics, 549, A3. BibTeX: `Coccato2013`. Two-dimensional IFU decomposition and separate component-population maps.

- **Johnston et al. (2013)** — [Disentangling the stellar populations in the counter-rotating disc galaxy NGC 4550](https://doi.org/10.1093/mnras/sts121). Monthly Notices of the Royal Astronomical Society, 428, 1296–1302. BibTeX: `Johnston2013`. Long-slit decomposition and Monte Carlo uncertainty precedent; distinguish this dataset from the IFU studies.

- **Pizzella et al. (2018)** — [Evidence for the formation of the young counter-rotating stellar disk from gas acquired by IC 719](https://doi.org/10.1051/0004-6361/201731712). Astronomy & Astrophysics, 616, A22. BibTeX: `Pizzella2018`. MUSE/CaT component kinematics and populations; residual-based S/N discussion. Relevant precedent, not a verbatim scalar-calibration recipe.

- **Katkov et al. (2016)** — [Stellar counter-rotation in lenticular galaxy NGC 448](https://doi.org/10.1093/mnras/stw1452). Monthly Notices of the Royal Astronomical Society, 461, 2068–2079. BibTeX: `Katkov2016`. Long-slit decomposition and noisy-best-fit-spectrum/refitting uncertainty precedent; current Script 04 does not generate final error bars.

- **Bao et al. (2024)** — [Uncovering the formation of the counter-rotating stellar discs in SDSS J074834.64+444117.8](https://doi.org/10.1093/mnras/stae243). Monthly Notices of the Royal Astronomical Society, 528, 2643–2652. BibTeX: `Bao2024`. MaNGA two-component pPXF, CaT and population comparison; its orbit-based dynamical method is not adopted.

- **Rubino et al. (2021)** — [Detectability of large-scale counter-rotating stellar disks in galaxies with integral-field spectroscopy](https://doi.org/10.1051/0004-6361/202140702). Astronomy & Astrophysics, 654, A30. BibTeX: `Rubino2021`. Detectability/LOSVD-blending context for target selection and recovery tests; not a universal threshold for our data.


## Noise-calibration and statistical reading

- **Krajnović et al. (2015)** — [Unveiling the counter-rotating nature of the kinematically distinct core in NGC 5813 with MUSE](https://doi.org/10.1093/mnras/stv958). Monthly Notices of the Royal Astronomical Society, 452, 2–18. BibTeX: `Krajnovic2015NGC5813`. Chat-5 residual-noise/Monte Carlo pPXF precedent. Our best-grid-state MAD scalar is our specific implementation.

- **Zhuang et al. (2024)** — [Metals in Star-forming Galaxies with KCWI. I. Methodology and First Results on the Abundances of Iron, Magnesium, and Oxygen](https://doi.org/10.3847/1538-4357/ad5ff8). The Astrophysical Journal, 972, 182. BibTeX: `Zhuang2024KCWI`. KCWI reduction and noise/SNR methodology comparison. Our fixed-stack replay measures its own correction instead of importing their coefficient.

- **Hogg et al. (2010)** — [Data analysis recipes: Fitting a model to data](https://arxiv.org/abs/1008.4686). arXiv:1008.4686. BibTeX: `HoggBovyLang2010`. Chat-4 reading for Gaussian-noise likelihoods and model-fitting interpretation.

- **Andrae (2010)** — [Error estimation in astronomy: A guide](https://arxiv.org/abs/1009.2755). arXiv:1009.2755. BibTeX: `Andrae2010`. Chat-4 reading on grid/likelihood errors and Monte Carlo; no implication that final pipeline uncertainties are already calibrated.

- **Andrae et al. (2010)** — [Dos and don'ts of reduced chi-squared](https://arxiv.org/abs/1012.3754). arXiv:1012.3754. BibTeX: `AndraeEtAl2010`. Reading behind the distinction between total chi-square differences and reduced-chi-square diagnostics.


## Instrument, atmosphere and upstream software

- **Morrissey et al. (2018)** — [The Keck Cosmic Web Imager Integral Field Spectrograph](https://doi.org/10.3847/1538-4357/aad597). The Astrophysical Journal, 864, 93. BibTeX: `Morrissey2018KCWI`. KCWI instrument reference; red-channel configuration and actual software versions need separate provenance.

- **Soto et al. (2016)** — [ZAP -- Enhanced PCA Sky Subtraction for Integral Field Spectroscopy](https://doi.org/10.1093/mnras/stw474). Monthly Notices of the Royal Astronomical Society, 458, 3210–3220. BibTeX: `Soto2016ZAP`. Upstream PCA sky subtraction. The Gaussian-noise experiment does not refit the ZAP model for every realization.

- **Hanuschik (2003)** — [A flux-calibrated, high-resolution atlas of optical sky emission from UVES](https://doi.org/10.1051/0004-6361:20030885). Astronomy & Astrophysics, 407, 1157–1164. BibTeX: `Hanuschik2003`. External atmospheric-line identification/cross-check in Chat 4. Actual masks use observation-specific sky/telluric information; the full atlas is not automatically loaded.

- **Chen et al. (2021)** — [The KBSS–KCWI survey: the connection between extended Ly α haloes and galaxy azimuthal angle at z ∼ 2–3](https://doi.org/10.1093/mnras/stab2383). Monthly Notices of the Royal Astronomical Society, 508, 19–43. BibTeX: `Chen2021KCWI`. Explicitly requested by KcwiKit's CITATION.cff; also a KCWI reduction/noise comparison. Its empirical bin correction is not copied into our replay calibration.

- **Prochaska et al. (2020)** — [PypeIt: The Python Spectroscopic Data Reduction Pipeline](https://doi.org/10.21105/joss.02308). Journal of Open Source Software, 5, 2308. BibTeX: `Prochaska2020PypeIt`. Upstream PypeIt/telluric context for KSkyWizard; this project uses KCWI_DRP for detector reduction.

- **Prusinski & Chen (2024)** — [KCWIKit: KCWI Post-Processing and Improvements](https://ascl.net/2404.003). Astrophysics Source Code Library, ascl:2404.003. BibTeX: `PrusinskiChen2024KCWIKit`. Actual stacking/post-processing software. Authors/ASCL record verified in official and local CITATION.cff; record the release/commit used.

- **Neill et al. (2023)** — [KCWI_DRP: Keck Cosmic Web Imager Data Reduction Pipeline](https://ascl.net/2301.019). Astrophysics Source Code Library, ascl:2301.019. BibTeX: `Neill2023KCWIDRP`. Upstream detector-reduction software record; keep installed release and patches in actual-run provenance.

- **KSkyWizard contributors (software repository)** — [KSkyWizard: sky subtraction and telluric correction for KCWI/KCRM](https://github.com/zhuyunz/KSkyWizard). BibTeX: `KSkyWizardSoftware`. Upstream software link. No dedicated peer-reviewed KSkyWizard paper was established in this audit; do not invent publication metadata.


## MaNGA data and CRD sample context

- **Bundy et al. (2015)** — [Overview of the SDSS-IV MaNGA Survey: Mapping Nearby Galaxies at Apache Point Observatory](https://doi.org/10.1088/0004-637X/798/1/7). The Astrophysical Journal, 798, 7. BibTeX: `Bundy2015MaNGA`. Survey provenance for MaNGA targets and ancillary maps.

- **Westfall et al. (2019)** — [The Data Analysis Pipeline for the SDSS-IV MaNGA IFU Galaxy Survey: Overview](https://doi.org/10.3847/1538-3881/ab44a2). The Astronomical Journal, 158, 231. BibTeX: `Westfall2019MaNGADAP`. Provenance/interpretation of supplied MaNGA kinematic maps; distinct from our custom CRD_DAP.

- **Gasymov et al. (2025)** — [Stellar CoRGI in MaNGA: Stellar Counterrotating Galaxies Identified in the MaNGA Survey](https://doi.org/10.3847/1538-4365/ae0590). The Astrophysical Journal Supplement Series, 281, 19. BibTeX: `Gasymov2025CoRGI`. CRD sample/LOSVD and selective spectral-decomposition comparison; NBURSTS is not our fitting engine.

- **Piper & Crocker (2026)** — [A New Identification Method and Sample of Counter-Rotating Disk Galaxies in SDSS-IV MaNGA DR17](https://arxiv.org/abs/2603.01925). arXiv:2603.01925. BibTeX: `PiperCrocker2026`. User project and target-identification context from the earlier chats; cited verified arXiv version rather than inferred journal metadata.

- **Bao et al. (2026)** — [Properties of Galaxies with Counter-rotating Stellar Disks in the MaNGA Survey](https://arxiv.org/abs/2603.04044). arXiv:2603.04044. BibTeX: `Bao2026`. Comparison sample and Chat-5 target literature checks. Its classifications do not add optimizer constraints; cited arXiv version.


## Historical alternatives — not current pipeline methods

- **Cappellari & Copin (2003)** — [Adaptive spatial binning of integral-field spectroscopic data using Voronoi tessellations](https://doi.org/10.1046/j.1365-8711.2003.06541.x). Monthly Notices of the Royal Astronomical Society, 342, 345–354. BibTeX: `CappellariCopin2003Voronoi`. Legacy binning background. Current Script 02 uses PowerBin with no silent VorBin fallback.

- **Greco & Brandt (2016)** — [The Measurement, Treatment, and Impact of Spectral Covariance and Bayesian Priors in Integral-Field Spectroscopy of Exoplanets](https://doi.org/10.3847/1538-4357/833/2/134). The Astrophysical Journal, 833, 134. BibTeX: `GrecoBrandt2016`. Chat-4 covariance discussion in exoplanet IFS; not a KCWI/CRD calibration prescription or our current likelihood model.

- **Czekala et al. (2015)** — [Constructing a Flexible Likelihood Function for Spectroscopic Inference](https://doi.org/10.1088/0004-637X/812/2/128). The Astrophysical Journal, 812, 128. BibTeX: `Czekala2015`. Correlated-residual kernels and noise-versus-model-mismatch discussion; those kernels are not retained in the scalar-diagonal baseline.

- **Liu et al. (2020)** — [Covariance-regularized Reconstruction of Data Cubes in Integral Field Spectroscopy and Application to MaNGA Data](https://doi.org/10.3847/1538-3881/ab5aeb). The Astronomical Journal, 159, 22. BibTeX: `Liu2020CubeCovariance`. Cube-resampling/covariance discussion. We do not use this alternative reconstruction; its preprint appeared in 2019.

- **Cappellari (2002)** — [Efficient Multi-Gaussian Expansion of galaxies](https://doi.org/10.1046/j.1365-8711.2002.05412.x). Monthly Notices of the Royal Astronomical Society, 333, 400–410. BibTeX: `Cappellari2002MGE`. Alternative intrinsic-light modeling discussed before selecting the observed RH3 continuum proxy; not adopted.

- **Ogle et al. (2019)** — [A Break in Spiral Galaxy Scaling Relations at the Upper Limit of Galaxy Mass](https://doi.org/10.3847/2041-8213/ab459e). The Astrophysical Journal Letters. BibTeX: `Ogle2019SuperSpirals`. Extreme rotation-speed discussion only; not a universal physical bound and not the basis of our final 400 km/s ceiling.

- **Corsini (2014)** — [Counter-Rotation in Disk Galaxies](https://arxiv.org/abs/1403.1263). arXiv:1403.1263. BibTeX: `Corsini2014Review`. General formation/observational review linked in Chat 5; cited arXiv version, not an optimizer prescription.

## Citation corrections and scope

- **XookSuut:** the published paper is from **2024**, volume 60, pages 19–39. The earlier `LopezCoba2023XookSuut` key remains valid, but its metadata now describe the published article.

- **Coccato versus Johnston (2013):** Coccato's two-galaxy IFU study is A&A 549, A3. MNRAS 428, 1296–1302 is Johnston's NGC 4550 paper. Both are included with their correct links.

- **Pizzella et al. (2018):** the IC 719 paper is A&A 616, A22. The MNRAS 475, 232 link in an earlier discussion belongs to a different paper and is not used for this citation.

- **XSL:** DR3 and the SSP models are separate Verro et al. (2022) papers, so both entries are retained.

- **KcwiKit:** its citation instructions request Chen et al. (2021). That paper and the ASCL software record are included; template version/DOI placeholders from software metadata are not copied.

- The earlier shorthand “Pizzella 2004 review” could not be matched confidently to a specific publication, so no guessed entry is added. The verified Corsini review is linked above.

- The scalar-diagonal noise baseline does not adopt the historical covariance kernels, covariance-aware cube reconstruction, VorBin or MGE alternatives. Later recovery/error-bar stages remain separate work; listing their precedents does not assert that they have already been implemented or validated.

- This index excludes unrelated observing logistics and coding/AI discussions, and does not import entire bibliographies from papers pasted into the old chats.
