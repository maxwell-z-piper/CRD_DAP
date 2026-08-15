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

TARGET_NAME = "CHANGE_ME"
MANGA_ID = "CHANGE_ME"
REDSHIFT = None  # dimensionless; required

# =============================================================================
# 2. REQUIRED INPUT PATHS
# =============================================================================
# The pipeline is intentionally explicit about calibration inputs. The master
# arcs are required because the adopted analysis measures the instrumental LSF
# empirically rather than relying only on nominal resolving power.

BL_CUBE = Path("/path/to/stacked_BL_icubew.fits")
RH3_CUBE = Path("/path/to/stacked_RH3_icubew.fits")

BL_MASTER_ARC = Path("/path/to/BL_master_arc.fits")
RH3_MASTER_ARC = Path("/path/to/RH3_master_arc.fits")

PYMORPH_VAC = Path("/path/to/pymorph_vac.fits")
XSL_TEMPLATE_LIBRARY = Path("/path/to/xsl_ssp_library")

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

SCIENCE_WAVELENGTH_MEDIUM = "vacuum"  # "air" or "vacuum"
TEMPLATE_WAVELENGTH_MEDIUM = "air"    # set to actual XSL convention used

# Record whether the reduced cube wavelength solution has already received a
# heliocentric/barycentric correction. Exact allowed values will be validated
# by crd_utils.validation once Script 1 is implemented.
SCIENCE_VELOCITY_FRAME = "UNKNOWN"

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

# KCWI DRP normally writes *_wavemap.fits, *_slicemap.fits, and *_posmap.fits
# products derived from the same master-arc root. Leave these as None to let
# Script 1 discover those files automatically beside BL_MASTER_ARC/RH3_MASTER_ARC.
# Set explicit paths only for non-standard reductions or renamed files.
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

# An extended galaxy cube cannot provide an unbiased PSF estimate by itself.
# Explicit values above take precedence. If they are None, Script 1 may use a
# recognized seeing/FWHM value from the science-cube header. The header-key
# search is recorded in the output provenance.
PSF_HEADER_KEYS = ("SEEING", "FWHM", "GUIDFWHM")

# =============================================================================
# 8. SCRIPT 1: DATA-QUALITY / NOISE SETTINGS
# =============================================================================

MIN_GOOD_WAVELENGTH_FRACTION = 0.80
BAD_CHANNEL_FRACTION_THRESHOLD = 0.50
REJECT_ANY_NONZERO_DRP_FLAG = True

# Continuum-center and registration diagnostics. The BL continuum peak is the
# default common coordinate origin because BL defines the master PowerBins;
# Script 4 later refits the kinematic center within explicit bounds.
CENTER_SMOOTH_SIGMA_PIX = 1.0
CENTER_CENTROID_MIN_PERCENTILE = 60.0
COMMON_CENTER_SOURCE = "BL_peak"  # currently supported: "BL_peak", "RH3_peak"
CENTER_WARNING_ARCSEC = 0.5
REGISTRATION_WARNING_ARCSEC = 0.25

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

BL_TARGET_SN = 40.0

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
