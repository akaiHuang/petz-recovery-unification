"""
Abstract backend protocol for tau-chrono.

Defines the interface that all quantum backends must implement.
Any object satisfying the ``TauChronoBackend`` protocol can be used
with tau-chrono's tomography and analysis pipelines.

Public API:
    - ``TauChronoBackend``: runtime-checkable Protocol for backends
    - ``BackendInfo``: metadata dataclass returned by ``backend.info()``
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import (
    Dict,
    List,
    Optional,
    Protocol,
    Sequence,
    runtime_checkable,
)

import numpy as np
from numpy.typing import NDArray


@dataclass
class BackendInfo:
    """Metadata for a quantum backend."""

    name: str
    num_qubits: int
    provider: str = "unknown"
    is_simulator: bool = False
    max_shots: int = 100_000
    basis_gates: List[str] = field(default_factory=list)
    extra: Dict[str, object] = field(default_factory=dict)


@runtime_checkable
class TauChronoBackend(Protocol):
    """Protocol that all tau-chrono backends must satisfy.

    Any object with the following attributes and methods can be used
    as a backend, regardless of its concrete class:

    Attributes
    ----------
    name : str
        Human-readable backend name.
    num_qubits : int
        Number of qubits available on this backend.

    Methods
    -------
    run_circuits(circuits, shots=1024)
        Execute a batch of circuits and return measurement counts.
    info()
        Return backend metadata.
    """

    @property
    def name(self) -> str:
        """Human-readable backend name."""
        ...

    @property
    def num_qubits(self) -> int:
        """Number of qubits available."""
        ...

    def run_circuits(
        self,
        circuits: Sequence,
        shots: int = 1024,
    ) -> List[Dict[str, int]]:
        """Execute circuits and return measurement counts.

        Parameters
        ----------
        circuits : sequence
            Circuits to execute.  The concrete type depends on the
            backend (Qiskit ``QuantumCircuit``, OpenQASM string, etc.).
        shots : int
            Number of measurement shots per circuit.

        Returns
        -------
        list of dict
            One counts dictionary per circuit.
            Keys are bitstrings (e.g. ``"01"``), values are counts.
        """
        ...

    def info(self) -> BackendInfo:
        """Return metadata about this backend."""
        ...
