"""
IBM Quantum backend adapter for tau-chrono.

Connects to IBM Quantum hardware and simulators via ``qiskit-ibm-runtime``.
Gracefully falls back if the package is not installed.

Public API:
    - ``IBMBackend``: adapter wrapping an IBM Quantum backend
    - ``list_ibm_backends``: list available IBM backends
    - ``get_ibm_backend``: connect to a specific backend by name

Requires:
    pip install qiskit qiskit-ibm-runtime
"""

from __future__ import annotations

from typing import Dict, List, Optional, Sequence

from .base import BackendInfo, TauChronoBackend

# --------------------------------------------------------------------------
# Availability check
# --------------------------------------------------------------------------

_IBM_AVAILABLE = False
_IBM_IMPORT_ERROR: Optional[str] = None

try:
    from qiskit_ibm_runtime import (
        QiskitRuntimeService,
        SamplerV2,
        EstimatorV2,
    )
    from qiskit.compiler import transpile
    from qiskit import QuantumCircuit

    _IBM_AVAILABLE = True
except ImportError as exc:
    _IBM_IMPORT_ERROR = (
        f"qiskit-ibm-runtime is required for IBM Quantum backends.\n"
        f"Install with: pip install qiskit qiskit-ibm-runtime\n"
        f"Original error: {exc}"
    )


def is_available() -> bool:
    """Return True if IBM Quantum dependencies are installed."""
    return _IBM_AVAILABLE


def _require_ibm() -> None:
    """Raise ImportError if IBM deps are missing."""
    if not _IBM_AVAILABLE:
        raise ImportError(_IBM_IMPORT_ERROR)


# --------------------------------------------------------------------------
# Backend discovery
# --------------------------------------------------------------------------

def list_ibm_backends(
    channel: str = "ibm_quantum",
    instance: Optional[str] = None,
    min_qubits: int = 0,
    simulator: bool = True,
    operational: bool = True,
) -> List[str]:
    """List available IBM Quantum backends.

    Parameters
    ----------
    channel : str
        IBM Quantum channel (``"ibm_quantum"`` or ``"ibm_cloud"``).
    instance : str, optional
        IBM Quantum instance (e.g. ``"hub/group/project"``).
    min_qubits : int
        Minimum number of qubits required.
    simulator : bool
        Include simulator backends.
    operational : bool
        Only include operational (non-retired) backends.

    Returns
    -------
    list of str
        Backend names.
    """
    _require_ibm()

    kwargs = {"channel": channel}
    if instance is not None:
        kwargs["instance"] = instance

    service = QiskitRuntimeService(**kwargs)
    backends = service.backends(
        min_num_qubits=min_qubits,
        simulator=simulator,
        operational=operational,
    )
    return [b.name for b in backends]


def get_ibm_backend(
    backend_name: str = "ibm_brisbane",
    channel: str = "ibm_quantum",
    instance: Optional[str] = None,
    optimization_level: int = 0,
) -> "IBMBackend":
    """Connect to an IBM Quantum backend.

    Parameters
    ----------
    backend_name : str
        Backend name (e.g. ``"ibm_brisbane"``, ``"ibm_sherbrooke"``).
    channel : str
        IBM Quantum channel.
    instance : str, optional
        IBM Quantum instance.
    optimization_level : int
        Transpiler optimization level (default 0 for tau-chrono:
        no gate rewriting).

    Returns
    -------
    IBMBackend
        A tau-chrono backend adapter.
    """
    _require_ibm()

    kwargs = {"channel": channel}
    if instance is not None:
        kwargs["instance"] = instance

    service = QiskitRuntimeService(**kwargs)
    raw_backend = service.backend(backend_name)

    return IBMBackend(
        raw_backend=raw_backend,
        service=service,
        optimization_level=optimization_level,
    )


# --------------------------------------------------------------------------
# IBMBackend adapter
# --------------------------------------------------------------------------

class IBMBackend:
    """Adapter wrapping an IBM Quantum backend for tau-chrono.

    Satisfies the ``TauChronoBackend`` protocol.

    Parameters
    ----------
    raw_backend
        The underlying ``qiskit_ibm_runtime`` backend object.
    service
        The ``QiskitRuntimeService`` instance.
    optimization_level : int
        Transpiler optimization level.  Default 0 (no gate rewriting)
        to preserve the exact circuit structure for tomography.
    """

    def __init__(
        self,
        raw_backend,
        service=None,
        optimization_level: int = 0,
    ):
        _require_ibm()
        self._backend = raw_backend
        self._service = service
        self._optimization_level = optimization_level

    @property
    def name(self) -> str:
        return self._backend.name

    @property
    def num_qubits(self) -> int:
        return self._backend.num_qubits

    @property
    def raw_backend(self):
        """Access the underlying qiskit-ibm-runtime backend object."""
        return self._backend

    def info(self) -> BackendInfo:
        """Return backend metadata."""
        config = getattr(self._backend, "configuration", None)
        basis_gates = []
        if config is not None and callable(config):
            cfg = config()
            basis_gates = getattr(cfg, "basis_gates", [])
        elif hasattr(self._backend, "basis_gates"):
            basis_gates = list(self._backend.basis_gates or [])

        is_sim = False
        name_lower = self._backend.name.lower()
        if "simulator" in name_lower or "aer" in name_lower:
            is_sim = True

        return BackendInfo(
            name=self._backend.name,
            num_qubits=self._backend.num_qubits,
            provider="IBM Quantum",
            is_simulator=is_sim,
            max_shots=100_000,
            basis_gates=basis_gates,
        )

    def run_circuits(
        self,
        circuits: Sequence,
        shots: int = 1024,
    ) -> List[Dict[str, int]]:
        """Execute circuits on the IBM backend using SamplerV2.

        Parameters
        ----------
        circuits : sequence of QuantumCircuit
            Qiskit circuits to execute.
        shots : int
            Shots per circuit.

        Returns
        -------
        list of dict
            Measurement counts for each circuit.
        """
        _require_ibm()

        # Transpile with the configured optimization level
        transpiled = transpile(
            list(circuits),
            self._backend,
            optimization_level=self._optimization_level,
        )

        # Use SamplerV2 for execution
        sampler = SamplerV2(backend=self._backend)

        # Submit all circuits
        results_list = []
        for circuit in transpiled:
            job = sampler.run([circuit], shots=shots)
            result = job.result()

            # Extract counts from SamplerV2 result
            counts = self._extract_counts(result, 0)
            results_list.append(counts)

        return results_list

    def run_circuits_batch(
        self,
        circuits: Sequence,
        shots: int = 1024,
    ) -> List[Dict[str, int]]:
        """Execute circuits as a single batch job (more efficient).

        Parameters
        ----------
        circuits : sequence of QuantumCircuit
            Qiskit circuits to execute.
        shots : int
            Shots per circuit.

        Returns
        -------
        list of dict
            Measurement counts for each circuit.
        """
        _require_ibm()

        transpiled = transpile(
            list(circuits),
            self._backend,
            optimization_level=self._optimization_level,
        )

        sampler = SamplerV2(backend=self._backend)
        job = sampler.run(transpiled, shots=shots)
        result = job.result()

        results_list = []
        for i in range(len(transpiled)):
            counts = self._extract_counts(result, i)
            results_list.append(counts)

        return results_list

    @staticmethod
    def _extract_counts(result, circuit_index: int) -> Dict[str, int]:
        """Extract counts dict from a SamplerV2 result.

        Handles different result formats across qiskit-ibm-runtime versions.
        """
        try:
            # SamplerV2 result format (qiskit-ibm-runtime >= 0.20)
            pub_result = result[circuit_index]
            data = pub_result.data

            # Try to get counts from the classical register
            if hasattr(data, "meas"):
                bitarray = data.meas
            elif hasattr(data, "c"):
                bitarray = data.c
            else:
                # Find first classical register
                for attr_name in dir(data):
                    if not attr_name.startswith("_"):
                        attr = getattr(data, attr_name)
                        if hasattr(attr, "get_counts"):
                            bitarray = attr
                            break
                else:
                    raise AttributeError(
                        "Cannot find classical register in result data"
                    )

            return dict(bitarray.get_counts())

        except (IndexError, AttributeError, TypeError):
            # Fallback: try legacy result format
            try:
                return result.get_counts(circuit_index)
            except Exception:
                raise RuntimeError(
                    f"Cannot extract counts from result for circuit "
                    f"{circuit_index}. Result type: {type(result)}"
                )

    def run_estimator(
        self,
        circuits: Sequence,
        observables: Sequence,
        shots: Optional[int] = None,
    ) -> List[float]:
        """Run circuits with EstimatorV2 for expectation values.

        Parameters
        ----------
        circuits : sequence of QuantumCircuit
            Circuits to run.
        observables : sequence
            Observable operators (SparsePauliOp or similar).
        shots : int, optional
            Shots per circuit.  None uses the default.

        Returns
        -------
        list of float
            Expectation values for each (circuit, observable) pair.
        """
        _require_ibm()

        transpiled = transpile(
            list(circuits),
            self._backend,
            optimization_level=self._optimization_level,
        )

        estimator = EstimatorV2(backend=self._backend)

        pubs = list(zip(transpiled, observables))
        if shots is not None:
            job = estimator.run(pubs, shots=shots)
        else:
            job = estimator.run(pubs)

        result = job.result()

        values = []
        for i in range(len(pubs)):
            pub_result = result[i]
            values.append(float(pub_result.data.evs))

        return values

    def __repr__(self) -> str:
        return (
            f"IBMBackend(name='{self.name}', "
            f"qubits={self.num_qubits}, "
            f"opt_level={self._optimization_level})"
        )
