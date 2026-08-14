# Science-driver scripts

This directory will contain the seven executable pipeline stages described in `../METHODS.md`.

The drivers should intentionally stay compact. They should orchestrate work by calling thoroughly documented functions in `crd_utils` rather than containing large blocks of analysis logic themselves.

Planned scripts:

```text
01_prepare_and_register_cubes.py
02_make_master_BL_bins.py
03_build_RH3_likelihood_cubes.py
04_fit_global_two_disk_model.py
05_extract_global_RH3_solution.py
06_fit_BL_two_component_populations.py
07_uncertainties_and_final_maps.py
```

Each script should:

1. load and validate a target configuration;
2. create or attach to a run directory;
3. initialize terminal + file logging;
4. print a concise run header describing inputs and major settings;
5. call the relevant `crd_utils` functions;
6. save numerical products before making optional presentation plots;
7. write a concise completion summary with timing and warnings;
8. return a non-zero exit status if a scientifically important stage fails, while preserving completed products whenever practical.
