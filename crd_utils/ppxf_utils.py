"""Low-level pPXF wrappers shared by kinematic and population stages.

All fits whose chi2 values enter a likelihood calculation must be unregularized.
A separate regularized fit may later be used only to visualize a smooth SFH
after the kinematic/fraction solution has been selected.
"""

from __future__ import annotations


def fit_single_losvd(*args, **kwargs):
    """Run one-component stellar pPXF with standardized masks/noise handling."""
    raise NotImplementedError("Implemented with Script 3 preliminary/control fits.")


def fit_fixed_two_component_state(*args, **kwargs):
    """Fit one explicit (V_A, V_B, f_A) profile-likelihood grid state."""
    raise NotImplementedError("Implemented with Script 3.")
