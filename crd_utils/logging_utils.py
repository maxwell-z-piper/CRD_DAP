"""Run-directory creation and pipeline logging.

Long RH3/BL pPXF calculations are intentionally allowed to run unattended.
Every science script therefore writes the same human-readable messages to the
terminal and to persistent log files. The master ``pipeline.log`` records the
complete run; step-specific logs make it easier to inspect one stage.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
import logging
from pathlib import Path
import sys
from typing import Optional

from .config import PipelineConfig, snapshot_config, write_config_manifest


@dataclass(frozen=True)
class RunContext:
    """Filesystem locations associated with one immutable pipeline run."""

    run_dir: Path
    logs_dir: Path
    products_dir: Path
    figures_dir: Path
    diagnostics_dir: Path
    metadata_dir: Path


def create_run_context(cfg: PipelineConfig, *, run_name: str | None = None) -> RunContext:
    """Create a timestamped output tree and snapshot the configuration."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    safe_target = str(cfg.TARGET_NAME).replace(" ", "_").replace("/", "-")
    name = run_name or f"{safe_target}_{timestamp}"

    run_dir = Path(cfg.RUNS_ROOT).expanduser().resolve() / name
    if run_dir.exists():
        raise FileExistsError(f"Run directory already exists: {run_dir}")

    logs = run_dir / "logs"
    products = run_dir / "products"
    figures = run_dir / "figures"
    diagnostics = run_dir / "diagnostics"
    metadata = run_dir / "metadata"

    for directory in (logs, products, figures, diagnostics, metadata):
        directory.mkdir(parents=True, exist_ok=False)

    snapshot_config(cfg, metadata)
    write_config_manifest(cfg, metadata / "config_manifest.json")

    return RunContext(run_dir, logs, products, figures, diagnostics, metadata)


def _formatter() -> logging.Formatter:
    return logging.Formatter(
        fmt="%(asctime)s | %(levelname)-8s | %(name)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def setup_pipeline_logger(
    run: RunContext,
    *,
    name: str = "CRD_DAP",
    level: int = logging.INFO,
) -> logging.Logger:
    """Create the master terminal + ``pipeline.log`` logger."""
    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.propagate = False

    # Remove stale handlers when the same Python interpreter starts a new run.
    for handler in list(logger.handlers):
        logger.removeHandler(handler)
        handler.close()

    console = logging.StreamHandler(sys.stdout)
    console.setLevel(level)
    console.setFormatter(_formatter())

    file_handler = logging.FileHandler(run.logs_dir / "pipeline.log")
    file_handler.setLevel(level)
    file_handler.setFormatter(_formatter())

    logger.addHandler(console)
    logger.addHandler(file_handler)
    return logger


def setup_step_logger(
    run: RunContext,
    step_name: str,
    *,
    master_logger_name: str = "CRD_DAP",
    level: int = logging.INFO,
) -> logging.Logger:
    """Create a step logger that writes to its own file and the master handlers.

    The child logger propagates to the already configured master logger, so a
    message appears in the terminal and ``pipeline.log``. A dedicated file
    handler additionally records only this stage.
    """
    logger = logging.getLogger(f"{master_logger_name}.{step_name}")
    logger.setLevel(level)
    logger.propagate = True

    for handler in list(logger.handlers):
        logger.removeHandler(handler)
        handler.close()

    step_file = logging.FileHandler(run.logs_dir / f"{step_name}.log")
    step_file.setLevel(level)
    step_file.setFormatter(_formatter())
    logger.addHandler(step_file)
    return logger


def log_section(logger: logging.Logger, title: str, char: str = "-") -> None:
    """Write a visually obvious section divider to terminal and log files."""
    line = char * max(12, len(title))
    logger.info(line)
    logger.info(title)
    logger.info(line)
