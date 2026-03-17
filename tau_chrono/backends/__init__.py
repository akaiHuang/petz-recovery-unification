"""
tau-chrono backend adapters.

Provides a unified interface for running quantum circuits on different
hardware and simulator backends.  All backends satisfy the
``TauChronoBackend`` protocol defined in ``base.py``.

Available backends:
    - ``SimulatorBackend``: local numpy-based simulator (always available)
    - ``IBMBackend``: IBM Quantum via qiskit-ibm-runtime (optional)
    - ``QuantumInspireBackend``: Quantum Inspire via qiskit-quantuminspire (optional)

Auto-detection:
    - ``available_backends()``: returns dict of available backend classes
    - ``SIMULATOR_AVAILABLE``, ``IBM_AVAILABLE``, ``QI_AVAILABLE``: flags

Usage:
    >>> from tau_chrono.backends import SimulatorBackend, available_backends
    >>> sim = SimulatorBackend(num_qubits=1)
    >>> print(available_backends())
    {'simulator': <class 'SimulatorBackend'>, ...}
"""

from __future__ import annotations

from typing import Dict, List, Type

# Base protocol (always available)
from .base import BackendInfo, TauChronoBackend

# Simulator backend (always available -- numpy only)
from .simulator import NoiseModel, SimCircuit, SimulatorBackend

SIMULATOR_AVAILABLE = True

# IBM Quantum backend (optional)
IBM_AVAILABLE = False
IBMBackend = None
try:
    from .ibm import IBMBackend as _IBMBackend
    from .ibm import get_ibm_backend, is_available as _ibm_check, list_ibm_backends

    IBM_AVAILABLE = _ibm_check()
    if IBM_AVAILABLE:
        IBMBackend = _IBMBackend
except ImportError:
    pass

# Quantum Inspire backend (optional)
QI_AVAILABLE = False
QuantumInspireBackend = None
try:
    from .quantum_inspire import (
        CircuitAnalysis,
        GateCharacterization,
        QuantumInspireBackend as _QIBackend,
        compare_gates,
        extract_gate_kraus,
        get_qi_backend,
        is_available as _qi_check,
        list_qi_backends,
        run_tau_analysis,
    )

    QI_AVAILABLE = _qi_check()
    if QI_AVAILABLE:
        QuantumInspireBackend = _QIBackend
except ImportError:
    pass


def available_backends() -> Dict[str, Type]:
    """Return a dict of available backend classes.

    Returns
    -------
    dict
        Mapping from backend name to backend class.
        Always includes ``"simulator"``.  Other backends
        are included only if their dependencies are installed.

    Examples
    --------
    >>> from tau_chrono.backends import available_backends
    >>> backends = available_backends()
    >>> print(list(backends.keys()))
    ['simulator']  # or ['simulator', 'ibm', 'quantum_inspire'] if deps installed
    """
    result: Dict[str, Type] = {
        "simulator": SimulatorBackend,
    }

    if IBM_AVAILABLE and IBMBackend is not None:
        result["ibm"] = IBMBackend

    if QI_AVAILABLE and QuantumInspireBackend is not None:
        result["quantum_inspire"] = QuantumInspireBackend

    return result


def available_backend_names() -> List[str]:
    """Return a list of available backend names.

    Returns
    -------
    list of str
        Available backend names.
    """
    return list(available_backends().keys())


__all__ = [
    # Base
    "TauChronoBackend",
    "BackendInfo",
    # Simulator (always available)
    "SimulatorBackend",
    "SimCircuit",
    "NoiseModel",
    "SIMULATOR_AVAILABLE",
    # IBM (optional)
    "IBMBackend",
    "IBM_AVAILABLE",
    "get_ibm_backend",
    "list_ibm_backends",
    # Quantum Inspire (optional)
    "QuantumInspireBackend",
    "QI_AVAILABLE",
    "get_qi_backend",
    "list_qi_backends",
    "extract_gate_kraus",
    "run_tau_analysis",
    "compare_gates",
    "GateCharacterization",
    "CircuitAnalysis",
    # Discovery
    "available_backends",
    "available_backend_names",
]
