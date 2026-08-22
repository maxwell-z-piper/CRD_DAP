# Science-driver scripts

This directory contains the executable pipeline stages described in `../METHODS.md`. Scripts 1--3 are implemented; later stages are added sequentially as their helper functions are completed and tested.

The drivers should intentionally stay compact. They should orchestrate work by calling thoroughly documented functions in `crd_utils` rather than containing large blocks of analysis logic themselves.

Science-driver sequence:

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

## Script 3 command-line behavior

`03_build_RH3_likelihood_cubes.py` supports:

```bash
python scripts/03_build_RH3_likelihood_cubes.py \
    --config config/<target>.py \
    --script1-run runs/<script1-run> \
    --script2-run runs/<script2-run> \
    --workers 3
```


`--script1-run` is optional when the Script-2 manifest still points to the correct Script-1 run. Use it explicitly when upstream run directories have been moved or renamed. If the explicit Script-1 path differs from the path recorded by Script 2, Script 3 logs a provenance warning but uses the explicit path.

`--workers N` means N Python worker processes (approximately N compute cores); BLAS/OpenMP thread pools are capped at one thread per worker. CPU affinity is not pinned, so the operating system still controls scheduling.

An interrupted run can be continued with `--resume`. One PowerBin is one atomic checkpoint. The configuration hash is verified before any checkpoint is reused. By default the intermediary checkpoint directory is deleted only after every final Script-3 product and the manifest have been written successfully.

The baseline driver intentionally does **not** expose `--bins` or `--max-bins`; add a targeted diagnostic interface later only if troubleshooting demonstrates a real need.
