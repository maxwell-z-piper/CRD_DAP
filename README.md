# CRD_DAP

CRD_DAP prepares the KCWI blue and red stacks, defines BL PowerBins, and fits RH3 profile-likelihood cubes. The target-specific inputs and science controls are in `config/8143-1902.py`; the full sequence is described in `METHODS.md` and `scripts/README.md`.

The current uncertainty path has two measured scalars for every actual bin. Script 2 coadds independently drawn reduced-exposure noise cubes through the fixed KcwiKit stack and compares each bin’s mock variance with its formal binned variance. Script 3 fits the complete RH3 likelihood grid once with that spatially corrected noise, estimates a robust residual multiplier from the exact best grid state, and rescales the saved noise and chi-square surface by that multiplier, including when it is below one. The fit assumes diagonal covariance between wavelength samples. Script-1 spectral correlation remains a diagnostic; likelihood intervals are conditional on the diagonal assumption.

The red science stack should be built from the four KSkyWizard `_zap_icubes.fits` exposures. The earlier red stack used the pre-ZAP `_icubes.fits` exposures. The exact stacking and mock replay commands are in `/Users/maxpiper/CRD_DRP/noise/README.md`. KSkyWizard cropped the ZAP inputs to its good wavelength range; the new red stack therefore has less spectral coverage. Build a fresh blue stack from the recorded exposure as well; the saved blue flux cannot be reproduced exactly. Rebuild the CRD_DRP atmospheric masks and reduction manifest, then run DAP Scripts 1–3 anew. Prior Script-3 checkpoints use another grid and noise model and cannot be resumed.

From this repository, validation can be run in the `crd_dap` environment:

```bash
/Users/maxpiper/miniforge3/envs/crd_dap/bin/python -m pytest -q tests/test_scalar_noise_pipeline.py tests/test_ppxf_utils_masked_placeholders.py tests/test_script02_api_contract.py
```

Script 3's `RH3_log_spectra_and_local_best_fits.npz` stores `noise` as the final scalar-corrected likelihood denominator. pPXF evaluated the grid using `noise_before_residual_scale`; its saved chi-square surface was then divided by the square of the best-state residual multiplier. Do not multiply the saved `noise` or chi-square by either scalar again. Old Script-3 checkpoints use a different calibration method and cannot be resumed with this version. The Script-3 manifest marks publication-final likelihood widths false while final data and inference checks remain incomplete.
