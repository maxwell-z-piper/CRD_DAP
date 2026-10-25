# Target configuration files

Copy `target_config_template.py` to one file per target and edit that copy.

Example:

```bash
cp config/target_config_template.py config/9028-1901.py
```

The target file contains both required input paths and configurable numerical-analysis settings. Development defaults are labeled clearly and should not be treated as universal scientific values.

Every pipeline run saves a verbatim `config_snapshot.py` and JSON configuration manifest inside the run metadata directory.

Script 04 accepts reviewed MaNGA measurements directly in section
11 of the target config; supplying map files is optional:

- `MANGA_SIGMA_PEAK_POS_KMS` and `MANGA_SIGMA_PEAK_NEG_KMS`: instrumentally
  corrected dispersions at the two off-center peaks, in km/s.
- `MANGA_SIGMA_BASE_KMS`: a representative underlying component dispersion,
  using the same correction convention, in km/s.
- `R_2SIGMA_POS_INITIAL_ARCSEC` and `R_2SIGMA_NEG_INITIAL_ARCSEC`: positive
  projected major-axis distances of those peaks from the adopted center, in
  arcsec. POS/NEG labels spatial sides; it does not make the input values signed.

These values initialize the two-disk search; they do not impose narrow velocity
bounds. Leave placeholders unset until the measurements have been entered.
Scripts 01-03 do not consume these fields, so adding them does not change their
scientific calculations. Whole-config hashes will nevertheless change; do not
edit completed runs' saved config snapshots to match the new file.

Script 04's commands and review sequence are in [../SCRIPT04.md](../SCRIPT04.md).
`SCRIPT04_VROT_MAX_KMS = 400.0` is an intrinsic speed ceiling; a ceiling flag can
motivate changing it to 600 and rerunning this stage in a new output directory.
The saved Script-3 LOS grids still restrict supported predictions.

Select `SCRIPT04_PYMORPH_HDU` after reviewing the catalog photometric band, or
supply explicit inclination bounds and initializer. No catalog band is silently
selected. The PSF profile is also an explicit command-line choice after viewing
the standard-star QA. Numerical search budgets are exposed separately from the
scientific bounds and selected regularity strength.
