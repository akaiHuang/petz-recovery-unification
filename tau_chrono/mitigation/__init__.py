"""Error mitigation powered by tau-chrono noise characterization.

Provides tau-informed variants of standard error mitigation techniques:
- ZNE (Zero Noise Extrapolation) with selective gate amplification
- PEC (Probabilistic Error Cancellation) using Kraus-derived quasi-probabilities

These techniques leverage tau-chrono's per-gate noise characterization
to focus mitigation effort where it matters most.
"""

from .zne import tau_informed_zne, selective_amplify

__all__ = ["tau_informed_zne", "selective_amplify"]
