"""All diagnostic and publication plotting functions.

Every nontrivial plotting function added here must have a matching, detailed
entry in ``DIAGNOSTICS.md``. The docstring should explain the scientific purpose,
not merely the Matplotlib commands. Important numerical plot settings should be
written to optional JSON sidecar metadata for reproducibility.
"""

from __future__ import annotations


def save_plot_metadata(*args, **kwargs):
    """Save figure-generation metadata beside a diagnostic image."""
    raise NotImplementedError("Implemented with the first plotting functions in Script 1.")
