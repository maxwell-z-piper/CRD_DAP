"""I/O helpers for FITS cubes, tables, arrays, and provenance metadata.

This module will centralize file-format assumptions so science functions do not
scatter FITS-extension names and table-column conventions throughout the code.
Script 1 will be the first major consumer.
"""

from __future__ import annotations

from pathlib import Path
import json
from typing import Any


def write_json(payload: dict[str, Any], path: str | Path) -> Path:
    """Write an indented JSON metadata file, creating parent directories."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n")
    return path


def read_json(path: str | Path) -> dict[str, Any]:
    """Read a JSON metadata file."""
    return json.loads(Path(path).read_text())


def load_kcwi_cube(*args, **kwargs):
    """Load a KCWI/KCRM cube into the pipeline's internal data model.

    Implementation is deferred to Script 1 because it must be matched to the
    exact KCWI DRP product structure in the user's reduced data. The eventual
    return object should expose flux, variance/uncertainty, masks/flags,
    wavelength coordinates, WCS, and relevant calibration metadata without
    requiring downstream code to know FITS extension names.
    """
    raise NotImplementedError("Implemented with Script 1 after inspecting real KCWI DRP products.")
