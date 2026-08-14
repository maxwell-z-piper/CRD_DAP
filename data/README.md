# Data directory

Large science inputs should not be committed to Git.

The target configuration should point to the actual locations of:

- stacked BL cube;
- stacked RH3 cube;
- BL/RH3 master arc calibration products;
- PyMorph VAC or a local target extract;
- XSL SSP templates;
- optional MaNGA maps or ancillary target data used for initialization/validation.

The pipeline should never assume that these data are located inside the repository. Absolute paths or paths resolved relative to a user-defined data root are acceptable.
