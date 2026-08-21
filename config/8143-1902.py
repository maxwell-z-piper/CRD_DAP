"""Target-specific configuration template for CRD_DAP.

Copy this file to a new target-specific configuration, for example::

    config/9028-1901.py

and edit only that copy.

The science-driver scripts should import configuration values through
``crd_utils.config.load_config`` rather than importing this module directly.
That loader performs validation and writes a configuration snapshot into each
run directory so that every numerical result can be traced to the exact
settings that produced it.

IMPORTANT
---------
The values below are development defaults and placeholders. They are not all
scientifically validated choices. Settings that still require mock-data or
real-data calibration are explicitly marked as such.
"""

from pathlib import Path

# =============================================================================
# 1. TARGET IDENTITY
# =============================================================================

TARGET_NAME = "8143-1902"
MANGA_ID = "1-44047"
REDSHIFT = 0.04138  # dimensionless; required

# =============================================================================
# 2. REQUIRED INPUT PATHS
# =============================================================================
# Production CRD_DAP science cubes are expected to come from KcwiKit post-DRP
# stacking. KcwiKit writes four matched single-HDU cubes per arm:
#
#   *_icubes.fits  science flux
#   *_vcubes.fits  variance (NOT 1-sigma uncertainty)
#   *_mcubes.fits  final binary stack mask; 0=valid contribution, 1=no valid data
#   *_ecubes.fits  effective exposure time in seconds
#
# Keep these products together and never substitute one arm's companion cube
# for another. Script 1 verifies that all four files have identical shape and
# spatial/spectral WCS before using them.
SCIENCE_INPUT_FORMAT = "kcwikit"  # supported: "kcwikit" (production), "drp" (legacy/testing)

BL_ICUBE = Path("/Users/maxpiper/Desktop/CRD_Decomposition/CRD_DAP_TestData/test_blue_icubes.fits")
BL_VCUBE = Path("/Users/maxpiper/Desktop/CRD_Decomposition/CRD_DAP_TestData/test_blue_vcubes.fits")
BL_MCUBE = Path("/Users/maxpiper/Desktop/CRD_Decomposition/CRD_DAP_TestData/test_blue_mcubes.fits")
BL_ECUBE = Path("/Users/maxpiper/Desktop/CRD_Decomposition/CRD_DAP_TestData/test_blue_ecubes.fits")

RH3_ICUBE = Path("/Users/maxpiper/Desktop/CRD_Decomposition/CRD_DAP_TestData/test_red_icubes.fits")
RH3_VCUBE = Path("/Users/maxpiper/Desktop/CRD_Decomposition/CRD_DAP_TestData/test_red_vcubes.fits")
RH3_MCUBE = Path("/Users/maxpiper/Desktop/CRD_Decomposition/CRD_DAP_TestData/test_red_mcubes.fits")
RH3_ECUBE = Path("/Users/maxpiper/Desktop/CRD_Decomposition/CRD_DAP_TestData/test_red_ecubes.fits")

# Optional legacy/native-DRP input mode. These are ignored when
# SCIENCE_INPUT_FORMAT="kcwikit". A DRP cube must contain an UNCERT extension.
BL_CUBE = Path("/path/to/stacked_or_single_BL_icubew.fits")
RH3_CUBE = Path("/path/to/stacked_or_single_RH3_icubew.fits")

BL_MASTER_ARC = Path("/Users/maxpiper/Desktop/CRD_Decomposition/KCWI_TestData/kb241226_00021_marc.fits")
RH3_MASTER_ARC = Path("/Users/maxpiper/Desktop/CRD_Decomposition/KCWI_TestData/kr241226_00025_marc.fits")

PYMORPH_VAC = Path("/Users/maxpiper/Desktop/CRD_Decomposition/CRD_DAP_TestData/manga-pymorph-DR17.fits")
XSL_TEMPLATE_LIBRARY = Path("/Users/maxpiper/Desktop/CRD_Decomposition/CRD_DAP_TestData/spectra_xsl_9.0.npz")

# Optional ancillary products used for initialization or comparison. These are
# not required for the core KCWI analysis.
MANGA_SIGMA_MAP = None
MANGA_VELOCITY_MAP = None

# =============================================================================
# 3. OUTPUT / RUN DIRECTORY
# =============================================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]
RUNS_ROOT = PROJECT_ROOT / "runs"

# =============================================================================
# 4. WAVELENGTH AND VELOCITY CONVENTIONS
# =============================================================================
# Script 1 must verify these against FITS metadata and template metadata. A
# science/template medium difference is allowed only when it is explicit and
# recorded for conversion during template preparation; two arrays may never be
# passed to pPXF in inconsistent wavelength media because that can mimic a
# velocity zero-point offset.

# "auto" is the recommended production setting: Script 1 infers the science
# convention from the reduced cube metadata and hard-fails if it cannot do so.
# An explicit "air"/"vacuum" value is also allowed and is checked against the
# header rather than blindly trusted.
SCIENCE_WAVELENGTH_MEDIUM = "auto"
TEMPLATE_WAVELENGTH_MEDIUM = "air"    # set to actual XSL convention used

# Likewise, "auto" reads the applied radial-velocity correction from KCWI DRP
# metadata (notably VCORRTYP). Explicit values are allowed and cross-checked.
SCIENCE_VELOCITY_FRAME = "auto"

# =============================================================================
# 5. INITIAL GEOMETRY
# =============================================================================

# Stellar kinematic PA from fit_kinematic_pa or equivalent, in degrees.
PA_KIN_INITIAL_DEG = None
PA_KIN_ERR_DEG = None

# Systemic velocity initializer in km/s. This may come from the preliminary RH3
# single-component fit / fit_kinematic_pa correction. It is NOT treated as an
# exact measurement; V_sys remains a free global parameter.
VSYS_INITIAL_KMS = None

# Optional explicit initial center in prepared-cube coordinates (arcsec). If
# None, Script 1 derives a continuum-peak / centroid initializer.
X0_INITIAL_ARCSEC = None
Y0_INITIAL_ARCSEC = None

# Disk identity convention. The PA axis is assigned a reproducible sign so that
# the analysis can walk from negative to positive signed major-axis coordinate.
# Disk A is the outer-branch component encountered first on the reference side;
# Disk B is its counter-rotating counterpart.
DISK_LABEL_CONVENTION = "signed_pa_outer_branch_first"

# =============================================================================
# 6. PHOTOMETRIC INCLINATION PRIOR
# =============================================================================
# The inclination itself is derived from the PyMorph disk axis ratio rather
# than assumed to be provided directly by the VAC. The intrinsic disk thickness
# q0 treatment must be calibrated/justified before publication.

PYMORPH_USE_DISK_COMPONENT = True
Q0_INTRINSIC_DISK = 0.2       # DEVELOPMENT DEFAULT; validate scientifically
Q0_INTRINSIC_DISK_ERR = 0.05  # DEVELOPMENT DEFAULT
N_INCLINATION_PRIOR_DRAWS = 10000

# =============================================================================
# 7. PSF / LSF SETTINGS
# =============================================================================

# If None, derive empirically. Manual values may be supplied for testing only.
BL_PSF_FWHM_ARCSEC = None
RH3_PSF_FWHM_ARCSEC = None

LSF_MODE = "master_arc"  # baseline production mode
LSF_MODEL_WAVELENGTH_ORDER = 2
LSF_MEASURE_SPATIAL_VARIATION = True

# Expected grating in each pipeline data stream. The stream labels remain BL and
# RH3 because downstream science products are organized around the intended
# BL+RH3 experiment, but Script 1 can validate other setups without disabling
# calibration checks. For normal production observations leave these at BL/RH3.
# For an integration test using, e.g., an RL red cube and matching RL master arc,
# set RH3_EXPECTED_GRATING = "RL" in that target's config only.
BL_EXPECTED_GRATING = "BL"
RH3_EXPECTED_GRATING = "RL"

# KCWI DRP writes *_wavemap.fits, *_slicemap.fits, and *_posmap.fits geometry
# products associated with the master-arc wavelength solution. Leave these as
# None to let Script 1 discover them automatically. Discovery first tries the
# conventional filename root and then falls back to FITS-header provenance
# (OFNAME + instrumental setup), because real RED reductions can give the maps a
# different exposure number from the *_marc.fits filename. Explicit paths always
# remain available when a reduction contains multiple ambiguous candidates.
BL_ARC_WAVEMAP = None
BL_ARC_SLICEMAP = None
BL_ARC_POSMAP = None
RH3_ARC_WAVEMAP = None
RH3_ARC_SLICEMAP = None
RH3_ARC_POSMAP = None

# Master-arc line detection/fitting defaults. These are deliberately exposed
# because unusual lamp spectra or configurations may require tuning after the
# Script-1 LSF diagnostics are inspected.
LSF_SPATIAL_BINS = 3
ARC_PEAK_PROMINENCE_FRACTION = 0.05
ARC_PEAK_HEIGHT_PERCENTILE = 70.0
ARC_MIN_PEAK_DISTANCE_PIX = 4
ARC_LINE_FIT_HALF_WIDTH_PIX = 5
ARC_MIN_GOOD_LINES = 6
ARC_LSF_SIGMA_CLIP = 3.0
LSF_SPATIAL_VARIATION_WARNING_FRACTION = 0.10

# A polynomial LSF is only empirically supported between the bluest and reddest
# accepted arc-line measurements. If an instrument-good edge extends farther
# than this many Angstroms beyond that support, Script 1 raises a QC flag and
# downstream fitting must not silently extrapolate the LSF.
LSF_EDGE_EXTRAPOLATION_WARNING_ANGSTROM = 100.0

# Master arc and science cube must share camera, grating, slicer, binning, and
# central wavelength. This tolerance applies only to the central-wavelength
# header comparison; it is not an LSF tolerance.
ARC_SCIENCE_CWAVE_TOLERANCE_ANGSTROM = 0.5

# An extended galaxy cube cannot provide an unbiased PSF estimate by itself.
# Explicit values above take precedence. If they are None, Script 1 may use a
# recognized seeing/FWHM value from the science-cube header. The header-key
# search is recorded in the output provenance.
PSF_HEADER_KEYS = ("SEEING", "FWHM", "GUIDFWHM")

# =============================================================================
# 8. SCRIPT 1: DATA-QUALITY / NOISE SETTINGS
# =============================================================================

# KcwiKit stacks are stored as float64 by default and can be several GB per
# BL+RH3 dataset. Script 1 uses float32 in memory/prepared products by default to
# keep the full two-arm preparation practical on a workstation. Extracted 1-D
# spectra may later be promoted to float64 for pPXF. Set to "float64" if desired.
STACK_FLOAT_DTYPE = "float32"

MIN_GOOD_WAVELENGTH_FRACTION = 0.80
BAD_CHANNEL_FRACTION_THRESHOLD = 0.50

# This switch applies only to native DRP multi-extension input. KcwiKit already
# used the original PyDRP FLAGS while constructing its stack and the final
# mcubes file is a binary validity/coverage mask rather than the original bitmask.
REJECT_ANY_NONZERO_DRP_FLAG = True

# Continuum-center and registration diagnostics. The BL continuum peak is the
# default common coordinate origin because BL defines the master PowerBins;
# Script 4 later refits the kinematic center within explicit bounds.
CENTER_SMOOTH_SIGMA_PIX = 1.0
CENTER_CENTROID_MIN_PERCENTILE = 60.0
COMMON_CENTER_SOURCE = "BL_peak"  # currently supported: "BL_peak", "RH3_peak"
CENTER_WARNING_ARCSEC = 0.5
REGISTRATION_WARNING_ARCSEC = 0.25

# Registration-image construction and trust criteria. If the two arms share a
# sufficiently wide instrument-good wavelength interval, Script 1 collapses the
# same observed-frame wavelengths in both arms before cross-correlation. This
# suppresses false offsets from wavelength-dependent morphology. If there is no
# useful overlap, full-arm continuum images are used instead. In either case the
# numerical shift is trusted only when both images pass the morphology-contrast
# screen; otherwise Script 1 records REGISTRATION_INCONCLUSIVE and relies on the
# independent WCS/center comparison rather than reporting a spurious shift.
REGISTRATION_USE_COMMON_WAVELENGTH_IF_AVAILABLE = True
REGISTRATION_MIN_COMMON_RANGE_ANGSTROM = 50.0
REGISTRATION_MIN_COMMON_CHANNELS = 20
REGISTRATION_MIN_CONTRAST_SNR = 5.0
REGISTRATION_CONTRAST_SMOOTH_SIGMA_PIX = 1.0

# Script 1 measures only the *residual* shift after the two cubes have already
# been placed on sky WCS. Limit the morphology cross-correlation to this local
# radius so IFU edges or weak continuum structure cannot drive a spurious large
# translation. A correlation peak at this boundary is treated as inconclusive.
REGISTRATION_MAX_RESIDUAL_SHIFT_ARCSEC = 2.0

# Noise/covariance characterization. Exact implementation is developed in
# crd_utils.noise; these switches record the intended production behavior.
ESTIMATE_VARIANCE_RESCALING = True
ESTIMATE_SPECTRAL_COVARIANCE = True
USE_COVARIANCE_IN_PPXF_WHEN_SUPPORTED = True

# Script 1's covariance estimate is intentionally preliminary because no pPXF
# residual model exists yet. It uses high-pass residuals from low-continuum
# spaxels as a QC screen and is revisited after the first stellar fits.
NOISE_DIAGNOSTIC_MAX_SPAXELS = 200
NOISE_LOW_FLUX_PERCENTILE = 30.0
NOISE_SAVGOL_WINDOW = 31
NOISE_SAVGOL_POLYORDER = 2
NOISE_MAX_SPECTRAL_LAG = 15
NOISE_VARIANCE_SCALE_WARNING_FRACTION = 0.20
NOISE_CORRELATION_WARNING_ABS = 0.10
APPLY_PRELIMINARY_VARIANCE_RESCALING = False

# =============================================================================
# 9. SCRIPT 2: MASTER BL POWERBINS
# =============================================================================

BL_TARGET_SN = 7.0

# Rest-frame continuum windows used to define the BL binning S/N proxy and to
# report achieved S/N in the transferred red/RH3 spectra.  Script 2 converts
# these to observed frame with REDSHIFT and intersects them with Script-1's
# GOODWAVE mask.  These are S/N-measurement windows, not hard science-fit masks.
BL_BINNING_REST_RANGE_ANGSTROM = (4800.0, 5400.0)
RH3_SN_REST_RANGE_ANGSTROM = (8470.0, 8700.0)
BINNING_MIN_VALID_WINDOW_FRACTION = 0.50

# The PowerBin input domain should represent the useful stellar body rather than
# every exposed sky pixel in the KcwiKit canvas.  The automatic mode thresholds
# a spatially smoothed continuum-significance proxy, keeps the connected
# component nearest the Script-1 galaxy center, dilates it slightly, and then
# intersects it with Script-1's hard GOODSPAX mask.  Native low-S/N spaxels
# inside this aperture are retained for PowerBin to combine.
BINNING_APERTURE_MODE = "auto_connected_sn"  # supported: auto_connected_sn, all_good
BINNING_APERTURE_SMOOTH_SIGMA_PIX = 2.0
BINNING_APERTURE_SN_THRESHOLD = 2.0           # DEVELOPMENT DEFAULT; inspect diagnostic
BINNING_APERTURE_DILATE_PIX = 2
BINNING_APERTURE_CENTER_MAX_DISTANCE_PIX = 5.0
BINNING_APERTURE_MIN_PIXELS = 25
BINNING_APERTURE_MAX_RADIUS_ARCSEC = None      # optional explicit outer limit

# PowerBin is mandatory.  A callable capacity is used so a calibrated
# non-additive spatial-covariance law can be inserted later.  The current safe
# baseline does NOT invent such a correction: Script 1 measured spectral
# correlation, not a validated spatial covariance law for coadded spaxels.
POWERBIN_SPATIAL_COVARIANCE_MODE = "none"      # supported: none, log10
POWERBIN_SPATIAL_COVARIANCE_ALPHA = 0.0        # only used for mode="log10"
POWERBIN_REGUL = True
POWERBIN_MAXITER = 50
POWERBIN_VERBOSE = 1

# Transfer the BL-defined membership to the red/RH3 native grid through the two
# celestial WCS solutions.  The current experiment uses matched spatial
# sampling; Script 2 fails rather than silently pretending that strongly
# different pixel scales represent identical physical apertures.
BIN_TRANSFER_MAX_DISTANCE_ARCSEC = 0.25
BIN_TRANSFER_MAX_PIXEL_SCALE_FRACTIONAL_DIFFERENCE = 0.05
BIN_TRANSFER_MIN_ASSIGNED_FRACTION = 0.95

# A wavelength sample in a coadded bin spectrum is retained only when at least
# this fraction of that bin's member spaxels have a Script-1 GOODMASK sample.
BIN_SPECTRUM_MIN_MEMBER_FRACTION = 0.50

# Achieved-S/N QC uses a robust ratio-of-medians estimator rather than the
# older median(flux/uncertainty).  The latter is retained in the output table
# for audit/debugging because it can expose tiny-variance pathologies.  A bin
# with non-positive median continuum has no meaningful positive "achieved S/N",
# so its production-facing S/N is NaN while its signed diagnostic is preserved.
BIN_SN_MIN_GOOD_CHANNELS = 10
BIN_SN_REQUIRE_POSITIVE_CONTINUUM = True
BIN_SN_EXTREME_ABS_WARNING = 1000.0       # QC warning only; never clips data
BIN_SN_ESTIMATOR_DISAGREEMENT_FACTOR = 10.0  # compare robust vs legacy estimator

# Flag bins whose measured BL S/N is substantially below the requested target.
# PowerBin equalizes a continuum proxy, so a modest tolerance is expected.
BINNING_LOW_SN_WARNING_FRACTION = 0.80

# Shared upper color limit for BL/RH3 S/N comparison plots. The scale should be
# derived from the combined distribution from both arms.
SN_PLOT_UPPER_PERCENTILE = 95.0

# =============================================================================
# 10. SCRIPT 3: RH3 LIKELIHOOD CUBES
# =============================================================================

# Velocity-grid arrays will be constructed inclusively from min/max/number of
# points. The exact extent/resolution must later pass a grid-convergence test.
RH3_VA_MIN_KMS = -400.0
RH3_VA_MAX_KMS = 400.0
RH3_VA_N = 17

RH3_VB_MIN_KMS = -400.0
RH3_VB_MAX_KMS = 400.0
RH3_VB_N = 17

RH3_FA_MIN = 0.10
RH3_FA_MAX = 0.90
RH3_FA_STEP = 0.10

# Use the full selected XSL SSP grid for RH3 rather than reducing the template
# library purely for speed. Scientific completeness takes priority over runtime.
RH3_USE_FULL_XSL_SSP_GRID = True

# Fits whose chi^2 values enter likelihood calculations must be unregularized.
RH3_REGUL = 0.0

# =============================================================================
# 11. SCRIPT 4: GLOBAL XOOKSUUT-STYLE DISK MODEL
# =============================================================================

# Ring spacing defaults to one delivered RH3 PSF FWHM when None.
RING_SPACE_ARCSEC = None
RING_DELTA_FACTOR = 0.5  # delta = 0.5 * ring_space
R_START_FWHM_FACTOR = 1.0

# Initial 2sigma radii may be supplied from MaNGA but the production value is
# intended to be updated from the preliminary single-component RH3 sigma map.
R_2SIGMA_POS_INITIAL_ARCSEC = None
R_2SIGMA_NEG_INITIAL_ARCSEC = None

# Always support two radial-extent fits:
#   1. snapped outer 2sigma radius;
#   2. full usable KCWI aperture.
RUN_TWO_SIGMA_EXTENT_MODEL = True
RUN_FULL_APERTURE_MODEL = True

# No explicit smoothness regularization on ring velocities.
ROTATION_CURVE_SMOOTHNESS_PENALTY = 0.0

# Geometry bounds are expressed as configurable multiples of external/initial
# uncertainties where available. Exact production values should be calibrated.
PA_BOUND_SIGMA = 3.0
CENTER_BOUND_ARCSEC = 1.0       # DEVELOPMENT DEFAULT
VSYS_BOUND_KMS = 100.0          # generous development default
INCLINATION_PRIOR_SIGMA_FACTOR = 2.0

# =============================================================================
# 12. SCRIPT 5: EXACT RH3 SOLUTION / QUALITY CONTROL
# =============================================================================

GRID_EDGE_WARNING_CELLS = 2

# Fine local fraction refinement around the coarse RH3 solution.
RH3_FA_FINE_STEP = 0.01
RH3_FA_FINE_HALF_WIDTH = 0.10

# Runtime estimate printing for expensive refinement stages.
PRINT_RUNTIME_ESTIMATES = True

# =============================================================================
# 13. SCRIPT 6: BL TWO-COMPONENT POPULATIONS
# =============================================================================

# Stellar-population likelihood fits remain unregularized. A separate
# regularized SFH visualization may be produced only after the preferred
# kinematic/fraction solution has been established.
BL_REGUL_LIKELIHOOD = 0.0

# Keep polynomial treatment explicit and configurable. Baseline has no additive
# polynomial; multiplicative degree should be validated on mocks.
BL_DEGREE = -1
BL_MDEGREE = 8  # DEVELOPMENT DEFAULT; validate on mocks

FIT_BL_GAS_EMISSION = True

# BL stellar-light fraction grid: coarse mapping followed by local refinement.
BL_FA_COARSE_MIN = 0.10
BL_FA_COARSE_MAX = 0.90
BL_FA_COARSE_STEP = 0.10
BL_FA_FINE_STEP = 0.01
BL_FA_FINE_HALF_WIDTH = 0.10

# If the selected RH3 basin is small enough, evaluate every state exactly.
N_DIRECT = 200

# Otherwise draw RH3 states according to normalized RH3 relative-likelihood
# weights and increase the cumulative sample until population summaries are
# stable. These are numerical-analysis defaults, not physical assumptions.
RH3_SAMPLE_BATCH_SIZE = 500
RH3_MIN_DRAWS = 1000
RH3_MAX_DRAWS = 20000
N_CONSECUTIVE_CONVERGENCE = 3
CONVERGENCE_FRACTION_OF_CI = 0.05
CONVERGENCE_FLOOR_LOGAGE = 0.005
CONVERGENCE_FLOOR_METALLICITY = 0.005
CONVERGENCE_FLOOR_FA = 0.005
RANDOM_SEED = 12345

# =============================================================================
# 14. SCRIPT 7: MONTE CARLO / SYSTEMATIC UNCERTAINTIES
# =============================================================================

# End-to-end MC count is deliberately convergence-controlled. These are starting
# values only; a runtime benchmark should be performed after one full nominal
# pipeline run before the final publication calculation.
MC_MIN_REALIZATIONS = 200
MC_BATCH_SIZE = 50
MC_MAX_REALIZATIONS = 2000
MC_CONSECUTIVE_CONVERGENCE = 3

# Number of parallel workers. None means choose at runtime / command line.
N_WORKERS = None

# =============================================================================
# 15. DIAGNOSTIC SETTINGS
# =============================================================================

BASIN_DIAGNOSTIC_MASSES = [0.90, 0.95, 0.99, 0.995, 0.999]
SAVE_PLOT_METADATA_JSON = True
SAVE_PER_BIN_DIAGNOSTICS = True

# =============================================================================
# 16. DEVELOPMENT / DEBUGGING
# =============================================================================

OVERWRITE_EXISTING_PRODUCTS = False
FAIL_ON_WARNING = False
