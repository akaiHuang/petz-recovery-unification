"""
Quantum Inspire backend adapter for tau-chrono.

Connects to Quantum Inspire (quantum-inspire.com) via the Qiskit provider.
This module is the canonical location for QI integration; the legacy
``tau_chrono.qi_bridge`` module re-exports from here.

Public API:
    - ``QuantumInspireBackend``: adapter wrapping a QI backend
    - ``list_qi_backends``: list available QI backends
    - ``get_qi_backend``: connect to a specific QI backend
    - ``extract_gate_kraus``: extract Kraus ops via process tomography
    - ``run_tau_analysis``: full pipeline (tomography -> Kraus -> tau)
    - ``compare_gates``: compare tau across multiple gates
    - ``GateCharacterization``: result dataclass
    - ``CircuitAnalysis``: result dataclass

Requires:
    pip install quantuminspire qiskit-quantuminspire
    qi login  (opens browser for authentication)
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Dict, List, Optional, Sequence

import numpy as np
from numpy.typing import NDArray

from .base import BackendInfo, TauChronoBackend

# --------------------------------------------------------------------------
# Availability check
# --------------------------------------------------------------------------

_QI_AVAILABLE = False
_QI_IMPORT_ERROR: Optional[str] = None

try:
    from qiskit_quantuminspire.qi_provider import QIProvider
    from qiskit.compiler import transpile

    _QI_AVAILABLE = True
except ImportError as exc:
    _QI_IMPORT_ERROR = (
        f"qiskit-quantuminspire is required for Quantum Inspire backends.\n"
        f"Install: pip install quantuminspire qiskit-quantuminspire\n"
        f"Then run: qi login\n"
        f"Original error: {exc}"
    )


def is_available() -> bool:
    """Return True if Quantum Inspire dependencies are installed."""
    return _QI_AVAILABLE


def _require_qi() -> None:
    """Raise ImportError if QI deps are missing."""
    if not _QI_AVAILABLE:
        raise ImportError(_QI_IMPORT_ERROR)


# --------------------------------------------------------------------------
# Backend discovery
# --------------------------------------------------------------------------

def list_qi_backends() -> List[str]:
    """List available Quantum Inspire backends.

    Returns
    -------
    list of str
        Backend names (e.g. ``"QX emulator"``, ``"Starmon-7"``).
    """
    _require_qi()
    provider = QIProvider()
    backends = provider.backends()
    return [b.name for b in backends]


def get_qi_backend(
    backend_name: str = "QX emulator",
    optimization_level: int = 0,
) -> "QuantumInspireBackend":
    """Connect to a Quantum Inspire backend.

    Parameters
    ----------
    backend_name : str
        Backend name (e.g. ``"QX emulator"``, ``"Starmon-7"``, ``"Tuna-5"``).
    optimization_level : int
        Transpiler optimization level (default 0: no gate rewriting).

    Returns
    -------
    QuantumInspireBackend
        A tau-chrono backend adapter.
    """
    _require_qi()
    provider = QIProvider()
    raw_backend = provider.get_backend(backend_name)
    return QuantumInspireBackend(
        raw_backend=raw_backend,
        optimization_level=optimization_level,
    )


# --------------------------------------------------------------------------
# QuantumInspireBackend adapter
# --------------------------------------------------------------------------

class QuantumInspireBackend:
    """Adapter wrapping a Quantum Inspire backend for tau-chrono.

    Satisfies the ``TauChronoBackend`` protocol.

    Parameters
    ----------
    raw_backend
        The underlying Qiskit backend from QIProvider.
    optimization_level : int
        Transpiler optimization level (default 0).
    """

    def __init__(
        self,
        raw_backend,
        optimization_level: int = 0,
    ):
        _require_qi()
        self._backend = raw_backend
        self._optimization_level = optimization_level

    @property
    def name(self) -> str:
        return self._backend.name

    @property
    def num_qubits(self) -> int:
        try:
            config = self._backend.configuration()
            return config.n_qubits
        except AttributeError:
            return getattr(self._backend, 'num_qubits', 0)

    @property
    def raw_backend(self):
        """Access the underlying QI Qiskit backend."""
        return self._backend

    def info(self) -> BackendInfo:
        """Return backend metadata."""
        try:
            config = self._backend.configuration()
            is_sim = getattr(config, "simulator", False)
            basis_gates = getattr(config, "basis_gates", [])
            max_shots = getattr(config, "max_shots", 100_000)
            n_qubits = config.n_qubits
        except AttributeError:
            is_sim = False
            basis_gates = []
            max_shots = 100_000
            n_qubits = self.num_qubits

        return BackendInfo(
            name=self._backend.name,
            num_qubits=n_qubits,
            provider="Quantum Inspire",
            is_simulator=is_sim,
            max_shots=max_shots,
            basis_gates=list(basis_gates),
        )

    def run_circuits(
        self,
        circuits: Sequence,
        shots: int = 1024,
    ) -> List[Dict[str, int]]:
        """Execute circuits on the Quantum Inspire backend.

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
        _require_qi()

        transpiled = transpile(
            list(circuits),
            self._backend,
            optimization_level=self._optimization_level,
        )

        job = self._backend.run(transpiled, shots=shots)
        job.wait_for_final_state(timeout=600)
        result = job.result()

        return [result.get_counts(i) for i in range(len(transpiled))]

    def __repr__(self) -> str:
        return (
            f"QuantumInspireBackend(name='{self.name}', "
            f"qubits={self.num_qubits}, "
            f"opt_level={self._optimization_level})"
        )


# --------------------------------------------------------------------------
# Gate characterization (migrated from qi_bridge.py)
# --------------------------------------------------------------------------

@dataclass
class GateCharacterization:
    """Result of characterizing a single gate."""

    gate_name: str
    kraus_ops: List[NDArray]
    n_kraus: int
    tau_worst: float
    tau_avg: float
    cptp_error: float
    process_fidelity: float


@dataclass
class CircuitAnalysis:
    """Result of analyzing a circuit with real hardware noise."""

    circuit_name: str
    gate_names: List[str]
    gate_chars: List[GateCharacterization]
    composition: object  # CompositionResult
    tau_naive_total: float
    tau_bayesian_total: float
    improvement_percent: float


def extract_gate_kraus(
    gate_circuit_fn: Callable,
    backend,
    gate_name: str = "gate",
    shots: int = 4096,
    verbose: bool = True,
    timeout: int = 600,
) -> GateCharacterization:
    """Extract Kraus operators for a gate via process tomography.

    Parameters
    ----------
    gate_circuit_fn : callable
        Function(qc, qubit) that appends the gate to a circuit.
    backend
        Qiskit-compatible backend or tau-chrono backend adapter.
    gate_name : str
        Human-readable name.
    shots : int
        Shots per tomography circuit.
    verbose : bool
        Print progress.
    timeout : int
        Max wait time in seconds for job completion.

    Returns
    -------
    GateCharacterization
        Full characterization including Kraus operators.
    """
    from ..tomography import run_process_tomography_1q
    from ..petz import tau_parameter, apply_channel, fidelity

    # If it's a tau-chrono backend adapter, unwrap to get the raw backend
    raw = getattr(backend, "raw_backend", backend)

    kraus_ops = run_process_tomography_1q(
        gate_circuit_fn, raw, shots=shots, verbose=verbose, timeout=timeout
    )

    # CPTP error
    total = sum(K.conj().T @ K for K in kraus_ops)
    cptp_err = float(np.linalg.norm(total - np.eye(2)))

    # Compute tau for standard input states
    sigma = np.eye(2, dtype=complex) / 2
    test_states = [
        np.array([[1, 0], [0, 0]], dtype=complex),
        np.array([[0, 0], [0, 1]], dtype=complex),
        np.array([[0.5, 0.5], [0.5, 0.5]], dtype=complex),
    ]

    taus = [tau_parameter(rho, kraus_ops, sigma) for rho in test_states]
    tau_worst = max(taus)
    tau_avg = sum(taus) / len(taus)

    fids = []
    for rho in test_states:
        rho_out = apply_channel(rho, kraus_ops)
        fids.append(fidelity(rho, rho_out))
    f_avg = sum(fids) / len(fids)
    process_fidelity = (2 * f_avg + 1) / 3

    if verbose:
        print(f"\n--- {gate_name} characterization ---")
        print(f"  Kraus rank: {len(kraus_ops)}")
        print(f"  CPTP error: {cptp_err:.2e}")
        print(f"  tau (worst): {tau_worst:.4f}")
        print(f"  tau (avg):   {tau_avg:.4f}")
        print(f"  Process fidelity: {process_fidelity:.4f}")

    return GateCharacterization(
        gate_name=gate_name,
        kraus_ops=kraus_ops,
        n_kraus=len(kraus_ops),
        tau_worst=tau_worst,
        tau_avg=tau_avg,
        cptp_error=cptp_err,
        process_fidelity=process_fidelity,
    )


def run_tau_analysis(
    gate_fns: List[Callable],
    gate_names: List[str],
    backend,
    rho: Optional[NDArray] = None,
    sigma: Optional[NDArray] = None,
    shots: int = 4096,
    verbose: bool = True,
) -> CircuitAnalysis:
    """Full pipeline: characterize gates, then run Bayesian composition.

    Parameters
    ----------
    gate_fns : list of callable
        Gate functions, each with signature (qc, qubit).
    gate_names : list of str
        Names for each gate.
    backend
        Qiskit-compatible backend or tau-chrono backend adapter.
    rho : NDArray, optional
        Input state (default: |0>).
    sigma : NDArray, optional
        Reference state (default: I/2).
    shots : int
        Shots per tomography circuit.
    verbose : bool
        Print progress.

    Returns
    -------
    CircuitAnalysis
        Complete analysis with Bayesian vs naive comparison.
    """
    from ..bayesian import bayesian_compose

    if rho is None:
        rho = np.array([[1, 0], [0, 0]], dtype=complex)
    if sigma is None:
        sigma = np.eye(2, dtype=complex) / 2

    gate_chars = []
    for fn, name in zip(gate_fns, gate_names):
        if verbose:
            print(f"\n{'=' * 50}")
            print(f"Characterizing gate: {name}")
            print(f"{'=' * 50}")
        char = extract_gate_kraus(fn, backend, name, shots, verbose)
        gate_chars.append(char)

    channels = [gc.kraus_ops for gc in gate_chars]

    if verbose:
        print(f"\n{'=' * 50}")
        print("Running Bayesian composition analysis...")
        print(f"{'=' * 50}")

    composition = bayesian_compose(
        channels=channels,
        sigma_0=sigma,
        rho=rho,
        channel_names=gate_names,
    )

    if verbose:
        print(f"\n--- Circuit Analysis Results ---")
        print(f"  Gates: {' -> '.join(gate_names)}")
        print(f"  tau (naive/multiplicative): "
              f"{composition.tau_multiplicative_total:.4f}")
        print(f"  tau (Bayesian):             "
              f"{composition.tau_bayesian_total:.4f}")
        print(f"  Improvement:                "
              f"{composition.improvement_percent:.1f}%")
        print(f"  Composition inequality:     "
              f"{'HOLDS' if composition.composition_holds else 'VIOLATED'}")
        print(f"\n  Per-gate breakdown:")
        for gr in composition.gate_results:
            print(f"    {gr.channel_name}: tau_naive={gr.tau_naive:.4f}, "
                  f"tau_eff={gr.tau_eff:.4f}, class={gr.classification}")

    return CircuitAnalysis(
        circuit_name=" -> ".join(gate_names),
        gate_names=gate_names,
        gate_chars=gate_chars,
        composition=composition,
        tau_naive_total=composition.tau_multiplicative_total,
        tau_bayesian_total=composition.tau_bayesian_total,
        improvement_percent=composition.improvement_percent,
    )


def compare_gates(
    gate_fns: Dict[str, Callable],
    backend,
    shots: int = 4096,
    verbose: bool = True,
) -> Dict[str, GateCharacterization]:
    """Compare multiple gates on the same backend.

    Parameters
    ----------
    gate_fns : dict
        Mapping from gate name to gate function.
    backend
        Qiskit-compatible backend or tau-chrono backend adapter.
    shots : int
        Shots per tomography circuit.
    verbose : bool
        Print progress.

    Returns
    -------
    dict
        Mapping from gate name to GateCharacterization.
    """
    results = {}
    for name, fn in gate_fns.items():
        results[name] = extract_gate_kraus(fn, backend, name, shots, verbose)

    if verbose:
        print(f"\n{'=' * 60}")
        print("GATE COMPARISON SUMMARY")
        print(f"{'=' * 60}")
        print(f"{'Gate':<15} {'Kraus':>6} {'tau_worst':>10} "
              f"{'tau_avg':>10} {'F_proc':>8}")
        print("-" * 60)
        for name, gc in results.items():
            print(f"{name:<15} {gc.n_kraus:>6} {gc.tau_worst:>10.4f} "
                  f"{gc.tau_avg:>10.4f} {gc.process_fidelity:>8.4f}")

    return results
