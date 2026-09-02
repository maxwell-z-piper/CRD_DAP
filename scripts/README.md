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

While the likelihood workers are active, an interactive terminal shows a single-line heartbeat refreshed every 5 seconds, for example:

```text
/ RUNNING [###---------------------] 63/484 bins ( 13.0%) | workers=3 | elapsed 42:18 | last bin 00:11 ago | ETA 4:42:07
```

The progress unit is one completed PowerBin checkpoint, not an individual pPXF grid state. Before the first bin finishes the line shows `ETA warming up`. The heartbeat is written directly to the terminal with a carriage return, so it does not flood the permanent log file. It means the parent process is alive and has not received a worker exception; it cannot by itself prove that a numerically slow worker is scientifically healthy.

The two repetitive CAPFIT scalar-division `RuntimeWarning`s generated during otherwise valid pPXF optimization are suppressed narrowly by message and module. Pipeline logger warnings, pPXF exceptions, and unrelated Python warnings remain visible.

An interrupted run can be continued with `--resume`. One PowerBin is one atomic checkpoint. The configuration hash is verified before any checkpoint is reused. By default the intermediary checkpoint directory is deleted only after every final Script-3 product and the manifest have been written successfully.

The baseline driver intentionally does **not** expose `--bins` or `--max-bins`; add a targeted diagnostic interface later only if troubleshooting demonstrates a real need.

Before launching workers, Script 3 computes template padding from both the exact two-component velocity grid and the wider one-component control velocity bounds, adds the configured dispersion-kernel support and a small log-grid edge margin, and then performs a real pPXF wavelength-coverage preflight. The preflight must pass for the one-component control and the four extreme two-component velocity-grid corners on representative widest-coverage bins before multiprocessing begins.

Masked log-grid samples may remain `NaN` in the saved CRD_DAP science arrays. Because pPXF validates the complete galaxy/noise vectors before applying `goodpixels`, the common pPXF wrapper makes private copies and supplies finite positive placeholders **only at already-excluded indices**. Values on `goodpixels` are never repaired: an invalid fitted sample is a hard failure. The placeholders therefore have no likelihood weight and do not alter the saved science arrays.

### Script-3 pPXF failure diagnostics

A failure raised through `concurrent.futures.process._RemoteTraceback` does not
by itself mean that `--workers` caused the error. pPXF state exceptions occur in
worker processes, so `ProcessPoolExecutor` necessarily transports them back to
the parent process this way.

If every two-component state in a PowerBin fails, Script 3 reports the
one-component-control status, first failed `(V_A,V_B,f_A)` coordinate, exact
first pPXF exception, and counts of the most common state-level errors. For a
partially successful cube, the per-bin completion line separately reports the
number of pPXF failures and fixed-velocity mismatches.
