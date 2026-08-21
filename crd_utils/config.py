"""Configuration loading, validation, and run-time snapshot utilities.

The pipeline intentionally keeps target-specific values out of the science
functions. A target configuration is an ordinary Python file so that units,
comments, and derived path definitions remain transparent to the researcher.
This module loads that file into a read-only-ish namespace, validates the
entries needed by the pipeline, and writes a verbatim snapshot into each run
folder.

Why snapshot the configuration?
-------------------------------
A result such as a population map is only scientifically reproducible if we can
recover the exact velocity-grid resolution, fraction-grid spacing, polynomial
degree, PSF/LSF settings, convergence tolerances, and input paths used to make
it. Git history alone is insufficient because users can rerun old code with a
new configuration. Every run therefore keeps its own copy.
"""

from __future__ import annotations

from dataclasses import dataclass
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import json
import shutil
from types import ModuleType
from typing import Any, Iterable


REQUIRED_CONFIG_KEYS = (
    "TARGET_NAME",
    "REDSHIFT",
    "SCIENCE_INPUT_FORMAT",
    "BL_MASTER_ARC",
    "RH3_MASTER_ARC",
    "PYMORPH_VAC",
    "XSL_TEMPLATE_LIBRARY",
    "RUNS_ROOT",
)

_KCWIKIT_STACK_KEYS = (
    "BL_ICUBE", "BL_VCUBE", "BL_MCUBE", "BL_ECUBE",
    "RH3_ICUBE", "RH3_VCUBE", "RH3_MCUBE", "RH3_ECUBE",
)

_DRP_CUBE_KEYS = ("BL_CUBE", "RH3_CUBE")


@dataclass(frozen=True)
class PipelineConfig:
    """Thin wrapper around the loaded Python configuration module.

    Attribute access is forwarded to the underlying module. The wrapper keeps
    the source path so it can be copied into a run directory.
    """

    module: ModuleType
    source_path: Path

    def __getattr__(self, name: str) -> Any:
        return getattr(self.module, name)

    def as_dict(self) -> dict[str, Any]:
        """Return public uppercase configuration values as a plain dictionary."""
        return {
            key: value
            for key, value in vars(self.module).items()
            if key.isupper() and not key.startswith("_")
        }


def load_config(path: str | Path, *, validate: bool = True, strict_paths: bool = True) -> PipelineConfig:
    """Load a target configuration from an arbitrary Python file.

    Parameters
    ----------
    path
        Path to the target-specific Python configuration file.
    validate
        Run structural validation immediately after loading.
    strict_paths
        If True, required science/calibration paths must already exist. Set to
        False for dry runs, documentation builds, or unit tests that are only
        checking configuration syntax.
    """
    source = Path(path).expanduser().resolve()
    if not source.exists():
        raise FileNotFoundError(f"Configuration file does not exist: {source}")

    spec = spec_from_file_location("crd_target_config", source)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not construct import specification for {source}")

    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    cfg = PipelineConfig(module=module, source_path=source)

    if validate:
        validate_config(cfg, strict_paths=strict_paths)
    return cfg


def validate_config(cfg: PipelineConfig, *, strict_paths: bool = True) -> None:
    """Validate the entries required before a science run can begin.

    This function deliberately performs only configuration-level checks. More
    specialized scientific validation (wavelength conventions, FITS extensions,
    template/data LSF compatibility, etc.) belongs in ``crd_utils.validation``
    and is executed by the relevant pipeline stage.
    """
    missing = [key for key in REQUIRED_CONFIG_KEYS if not hasattr(cfg.module, key)]
    if missing:
        raise ValueError("Missing required configuration keys: " + ", ".join(missing))

    if cfg.REDSHIFT is None:
        raise ValueError("REDSHIFT must be specified in the target configuration.")
    if float(cfg.REDSHIFT) < 0:
        raise ValueError("REDSHIFT must be non-negative.")

    input_format = str(cfg.SCIENCE_INPUT_FORMAT).strip().lower()
    if input_format not in {"kcwikit", "drp"}:
        raise ValueError("SCIENCE_INPUT_FORMAT must be either 'kcwikit' or 'drp'.")

    format_keys = _KCWIKIT_STACK_KEYS if input_format == "kcwikit" else _DRP_CUBE_KEYS
    missing_format = [key for key in format_keys if not hasattr(cfg.module, key)]
    if missing_format:
        raise ValueError(
            f"SCIENCE_INPUT_FORMAT={input_format!r} requires configuration keys: "
            + ", ".join(missing_format)
        )

    if strict_paths:
        path_keys = (
            *format_keys,
            "BL_MASTER_ARC",
            "RH3_MASTER_ARC",
            "PYMORPH_VAC",
            "XSL_TEMPLATE_LIBRARY",
        )
        absent = []
        for key in path_keys:
            value = getattr(cfg, key)
            if value is None:
                absent.append(f"{key}=None")
                continue
            path = Path(value).expanduser()
            if not path.exists():
                absent.append(f"{key}={path}")
        if absent:
            raise FileNotFoundError(
                "Required input path(s) do not exist:\n  - " + "\n  - ".join(absent)
            )

    if getattr(cfg, "RING_DELTA_FACTOR", 0.5) <= 0:
        raise ValueError("RING_DELTA_FACTOR must be positive.")
    if getattr(cfg, "GRID_EDGE_WARNING_CELLS", 2) < 0:
        raise ValueError("GRID_EDGE_WARNING_CELLS cannot be negative.")
    if getattr(cfg, "N_DIRECT", 200) < 1:
        raise ValueError("N_DIRECT must be >= 1.")

    min_good = float(getattr(cfg, "MIN_GOOD_WAVELENGTH_FRACTION", 0.80))
    if not 0.0 <= min_good <= 1.0:
        raise ValueError("MIN_GOOD_WAVELENGTH_FRACTION must lie between 0 and 1.")
    bad_channel = float(getattr(cfg, "BAD_CHANNEL_FRACTION_THRESHOLD", 0.50))
    if not 0.0 <= bad_channel <= 1.0:
        raise ValueError("BAD_CHANNEL_FRACTION_THRESHOLD must lie between 0 and 1.")
    stack_dtype = str(getattr(cfg, "STACK_FLOAT_DTYPE", "float32")).lower()
    if stack_dtype not in {"float32", "float64"}:
        raise ValueError("STACK_FLOAT_DTYPE must be 'float32' or 'float64'.")

    # The pipeline stream names are BL/RH3, but Script 1 may be used to prepare
    # and validate other blue/red grating setups. Keep the historical BL/RH3
    # expectations as safe defaults while allowing target configs to override
    # either grating explicitly. The exact name is then cross-checked against
    # both the science and master-arc FITS headers.
    for key, default in (("BL_EXPECTED_GRATING", "BL"), ("RH3_EXPECTED_GRATING", "RH3")):
        grating = str(getattr(cfg, key, default)).strip()
        if not grating:
            raise ValueError(f"{key} must be a non-empty grating name.")
    if int(getattr(cfg, "LSF_MODEL_WAVELENGTH_ORDER", 2)) < 0:
        raise ValueError("LSF_MODEL_WAVELENGTH_ORDER cannot be negative.")
    if int(getattr(cfg, "ARC_MIN_GOOD_LINES", 6)) < 3:
        raise ValueError("ARC_MIN_GOOD_LINES must be >= 3.")
    if float(getattr(cfg, "REGISTRATION_WARNING_ARCSEC", 0.25)) <= 0:
        raise ValueError("REGISTRATION_WARNING_ARCSEC must be positive.")
    if float(getattr(cfg, "REGISTRATION_MIN_COMMON_RANGE_ANGSTROM", 50.0)) < 0:
        raise ValueError("REGISTRATION_MIN_COMMON_RANGE_ANGSTROM cannot be negative.")
    if int(getattr(cfg, "REGISTRATION_MIN_COMMON_CHANNELS", 20)) < 1:
        raise ValueError("REGISTRATION_MIN_COMMON_CHANNELS must be >= 1.")
    if float(getattr(cfg, "REGISTRATION_MIN_CONTRAST_SNR", 5.0)) <= 0:
        raise ValueError("REGISTRATION_MIN_CONTRAST_SNR must be positive.")
    if float(getattr(cfg, "REGISTRATION_CONTRAST_SMOOTH_SIGMA_PIX", 1.0)) < 0:
        raise ValueError("REGISTRATION_CONTRAST_SMOOTH_SIGMA_PIX cannot be negative.")
    max_reg = float(getattr(cfg, "REGISTRATION_MAX_RESIDUAL_SHIFT_ARCSEC", 2.0))
    if max_reg <= 0:
        raise ValueError("REGISTRATION_MAX_RESIDUAL_SHIFT_ARCSEC must be positive.")
    warn_reg = float(getattr(cfg, "REGISTRATION_WARNING_ARCSEC", 0.25))
    if max_reg <= warn_reg:
        raise ValueError(
            "REGISTRATION_MAX_RESIDUAL_SHIFT_ARCSEC must exceed "
            "REGISTRATION_WARNING_ARCSEC."
        )
    if float(getattr(cfg, "LSF_EDGE_EXTRAPOLATION_WARNING_ANGSTROM", 100.0)) < 0:
        raise ValueError("LSF_EDGE_EXTRAPOLATION_WARNING_ANGSTROM cannot be negative.")


def validate_input_paths(cfg: PipelineConfig, path_keys: Iterable[str]) -> dict[str, Path]:
    """Validate only the input paths required by a particular pipeline stage.

    The repository-wide configuration contains inputs that are not needed by
    every script.  Stage-specific validation avoids forcing, for example, the
    XSL library or an RH3 science cube to exist merely to exercise BL-only
    Script-1 development/validation mode.  Production science runs still request
    all paths required by that stage explicitly.
    """
    resolved: dict[str, Path] = {}
    absent: list[str] = []
    for key in path_keys:
        if not hasattr(cfg, key):
            absent.append(f"{key}=<missing config key>")
            continue
        value = getattr(cfg, key)
        if value is None:
            absent.append(f"{key}=None")
            continue
        path = Path(value).expanduser().resolve()
        if not path.exists():
            absent.append(f"{key}={path}")
        else:
            resolved[key] = path
    if absent:
        raise FileNotFoundError(
            "Required input path(s) for this stage do not exist:\n  - " + "\n  - ".join(absent)
        )
    return resolved


def snapshot_config(cfg: PipelineConfig, destination_dir: str | Path) -> Path:
    """Copy the exact configuration file into a run directory."""
    destination_dir = Path(destination_dir)
    destination_dir.mkdir(parents=True, exist_ok=True)
    destination = destination_dir / "config_snapshot.py"
    shutil.copy2(cfg.source_path, destination)
    return destination


def _json_safe(value: Any) -> Any:
    """Convert common configuration values to JSON-serializable objects."""
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    if isinstance(value, (list, tuple)):
        return [_json_safe(v) for v in value]
    if isinstance(value, dict):
        return {str(k): _json_safe(v) for k, v in value.items()}
    return repr(value)


def write_config_manifest(cfg: PipelineConfig, destination: str | Path) -> Path:
    """Write public configuration values to JSON for machine-readable provenance."""
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    payload = {key: _json_safe(value) for key, value in cfg.as_dict().items()}
    destination.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return destination
