"""
Quantum Inspire bridge for tau-chrono.

.. deprecated::
    This module is maintained for backward compatibility only.
    Use ``tau_chrono.backends.quantum_inspire`` instead.

All functionality has been moved to ``tau_chrono.backends.quantum_inspire``.
This module re-exports every public name from there so that existing code
like ``from tau_chrono.qi_bridge import get_qi_backend`` continues to work.
"""

from __future__ import annotations

import warnings as _warnings

_warnings.warn(
    "tau_chrono.qi_bridge is deprecated. "
    "Use tau_chrono.backends.quantum_inspire instead.",
    DeprecationWarning,
    stacklevel=2,
)

# Re-export everything from the new canonical location
from .backends.quantum_inspire import (  # noqa: F401
    CircuitAnalysis,
    GateCharacterization,
    QuantumInspireBackend,
    compare_gates,
    extract_gate_kraus,
    get_qi_backend,
    is_available,
    list_qi_backends,
    run_tau_analysis,
)

__all__ = [
    "GateCharacterization",
    "CircuitAnalysis",
    "QuantumInspireBackend",
    "get_qi_backend",
    "list_qi_backends",
    "extract_gate_kraus",
    "run_tau_analysis",
    "compare_gates",
    "is_available",
]
