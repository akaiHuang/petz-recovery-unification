"""
Local numpy-based simulator backend for tau-chrono.

Provides a zero-dependency (beyond numpy) quantum simulator that can
run 1-qubit and 2-qubit process tomography circuits.  Supports
optional noise channels (depolarizing, amplitude damping, dephasing)
applied after each gate.

This backend works entirely without Qiskit -- it interprets a simple
gate-list representation and simulates measurement sampling directly.

Public API:
    - ``SimulatorBackend``: main backend class
    - ``NoiseModel``: configurable noise specification
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import (
    Callable,
    Dict,
    List,
    Optional,
    Sequence,
    Tuple,
    Union,
)

import numpy as np
from numpy.typing import NDArray

from .base import BackendInfo, TauChronoBackend


# -----------------------------------------------------------------------
# Standard 1-qubit gates (matrix definitions)
# -----------------------------------------------------------------------

_I = np.eye(2, dtype=complex)
_X = np.array([[0, 1], [1, 0]], dtype=complex)
_Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
_Z = np.array([[1, 0], [0, -1]], dtype=complex)
_H = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
_S = np.array([[1, 0], [0, 1j]], dtype=complex)
_SDG = np.array([[1, 0], [0, -1j]], dtype=complex)
_T = np.array([[1, 0], [0, np.exp(1j * np.pi / 4)]], dtype=complex)
_TDG = np.array([[1, 0], [0, np.exp(-1j * np.pi / 4)]], dtype=complex)

GATE_MATRICES: Dict[str, NDArray] = {
    "id": _I,
    "i": _I,
    "x": _X,
    "y": _Y,
    "z": _Z,
    "h": _H,
    "s": _S,
    "sdg": _SDG,
    "t": _T,
    "tdg": _TDG,
}


def _rx(theta: float) -> NDArray:
    """Rotation around X axis."""
    return np.array([
        [np.cos(theta / 2), -1j * np.sin(theta / 2)],
        [-1j * np.sin(theta / 2), np.cos(theta / 2)],
    ], dtype=complex)


def _ry(theta: float) -> NDArray:
    """Rotation around Y axis."""
    return np.array([
        [np.cos(theta / 2), -np.sin(theta / 2)],
        [np.sin(theta / 2), np.cos(theta / 2)],
    ], dtype=complex)


def _rz(theta: float) -> NDArray:
    """Rotation around Z axis."""
    return np.array([
        [np.exp(-1j * theta / 2), 0],
        [0, np.exp(1j * theta / 2)],
    ], dtype=complex)


# -----------------------------------------------------------------------
# Noise model
# -----------------------------------------------------------------------

@dataclass
class NoiseModel:
    """Simple noise model for the simulator.

    Parameters
    ----------
    depolarizing_rate : float
        Single-qubit depolarizing error rate per gate (0 to 1).
        Applies the channel: rho -> (1-p)*rho + p/3*(X rho X + Y rho Y + Z rho Z)
    amplitude_damping_rate : float
        T1-like amplitude damping parameter gamma per gate (0 to 1).
    dephasing_rate : float
        T2-like dephasing parameter per gate (0 to 1).
    readout_error : float
        Probability of a bit-flip on measurement (0 to 0.5).
    """

    depolarizing_rate: float = 0.0
    amplitude_damping_rate: float = 0.0
    dephasing_rate: float = 0.0
    readout_error: float = 0.0

    def apply(self, rho: NDArray) -> NDArray:
        """Apply all noise channels to a density matrix."""
        if self.depolarizing_rate > 0:
            p = self.depolarizing_rate
            rho = (1 - p) * rho + (p / 3) * (
                _X @ rho @ _X + _Y @ rho @ _Y + _Z @ rho @ _Z
            )
        if self.amplitude_damping_rate > 0:
            gamma = self.amplitude_damping_rate
            K0 = np.array([[1, 0], [0, np.sqrt(1 - gamma)]], dtype=complex)
            K1 = np.array([[0, np.sqrt(gamma)], [0, 0]], dtype=complex)
            rho = K0 @ rho @ K0.conj().T + K1 @ rho @ K1.conj().T
        if self.dephasing_rate > 0:
            p = self.dephasing_rate
            rho = (1 - p) * rho + p * _Z @ rho @ _Z
        return rho


# -----------------------------------------------------------------------
# Circuit representation (Qiskit-free)
# -----------------------------------------------------------------------

@dataclass
class SimCircuit:
    """Lightweight circuit representation for the simulator.

    Each instruction is a tuple of (gate_name, qubit_indices, params).
    """

    num_qubits: int
    num_clbits: int
    instructions: List[Tuple[str, List[int], List[float]]] = field(
        default_factory=list
    )
    name: str = ""

    def add_gate(
        self,
        gate: str,
        qubits: List[int],
        params: Optional[List[float]] = None,
    ) -> None:
        """Append a gate instruction."""
        self.instructions.append((gate.lower(), qubits, params or []))

    def add_measurement(self, qubit: int, clbit: int) -> None:
        """Append a measurement instruction."""
        self.instructions.append(("measure", [qubit, clbit], []))


# -----------------------------------------------------------------------
# Qiskit circuit parser
# -----------------------------------------------------------------------

def _parse_qiskit_circuit(qc) -> SimCircuit:
    """Convert a Qiskit QuantumCircuit to our SimCircuit representation.

    Parameters
    ----------
    qc : qiskit.circuit.QuantumCircuit
        The Qiskit circuit to convert.

    Returns
    -------
    SimCircuit
        Equivalent lightweight circuit.
    """
    sc = SimCircuit(
        num_qubits=qc.num_qubits,
        num_clbits=qc.num_clbits,
        name=getattr(qc, "name", ""),
    )

    for instruction in qc.data:
        # Handle both old-style (op, qargs, cargs) and new-style
        # CircuitInstruction
        if hasattr(instruction, "operation"):
            op = instruction.operation
            qargs = instruction.qubits
            cargs = instruction.clbits
        else:
            op, qargs, cargs = instruction

        gate_name = op.name.lower()

        if gate_name == "barrier":
            continue

        # Get qubit indices
        qubit_indices = []
        for q in qargs:
            if hasattr(q, "_index"):
                qubit_indices.append(q._index)
            elif hasattr(qc, "find_bit"):
                qubit_indices.append(qc.find_bit(q).index)
            else:
                qubit_indices.append(int(q))

        if gate_name == "measure":
            clbit_indices = []
            for c in cargs:
                if hasattr(c, "_index"):
                    clbit_indices.append(c._index)
                elif hasattr(qc, "find_bit"):
                    clbit_indices.append(qc.find_bit(c).index)
                else:
                    clbit_indices.append(int(c))
            for qi, ci in zip(qubit_indices, clbit_indices):
                sc.add_measurement(qi, ci)
        else:
            params = [float(p) for p in getattr(op, "params", [])]
            sc.add_gate(gate_name, qubit_indices, params)

    return sc


# -----------------------------------------------------------------------
# Density matrix simulation engine
# -----------------------------------------------------------------------

def _simulate_circuit(
    circuit: SimCircuit,
    noise: Optional[NoiseModel] = None,
    shots: int = 1024,
    rng: Optional[np.random.Generator] = None,
) -> Dict[str, int]:
    """Simulate a circuit using density matrix evolution.

    Parameters
    ----------
    circuit : SimCircuit
        The circuit to simulate.
    noise : NoiseModel, optional
        Noise model to apply after each gate.
    shots : int
        Number of measurement shots.
    rng : np.random.Generator, optional
        Random number generator for reproducibility.

    Returns
    -------
    dict
        Measurement counts (bitstring -> count).
    """
    if rng is None:
        rng = np.random.default_rng()

    n = circuit.num_qubits
    dim = 2 ** n

    # Initial state: |0...0>
    rho = np.zeros((dim, dim), dtype=complex)
    rho[0, 0] = 1.0

    # Measurement map: clbit -> qubit
    meas_map: Dict[int, int] = {}

    for gate_name, qubits, params in circuit.instructions:
        if gate_name == "measure":
            qubit_idx, clbit_idx = qubits
            meas_map[clbit_idx] = qubit_idx
            continue

        # Get gate matrix
        if gate_name in GATE_MATRICES:
            U = GATE_MATRICES[gate_name]
        elif gate_name == "rx" and len(params) >= 1:
            U = _rx(params[0])
        elif gate_name == "ry" and len(params) >= 1:
            U = _ry(params[0])
        elif gate_name == "rz" and len(params) >= 1:
            U = _rz(params[0])
        elif gate_name == "cx" or gate_name == "cnot":
            # 2-qubit CNOT gate handled below
            U = np.array([
                [1, 0, 0, 0],
                [0, 1, 0, 0],
                [0, 0, 0, 1],
                [0, 0, 1, 0],
            ], dtype=complex)
        elif gate_name == "cz":
            U = np.array([
                [1, 0, 0, 0],
                [0, 1, 0, 0],
                [0, 0, 1, 0],
                [0, 0, 0, -1],
            ], dtype=complex)
        elif gate_name == "swap":
            U = np.array([
                [1, 0, 0, 0],
                [0, 0, 1, 0],
                [0, 1, 0, 0],
                [0, 0, 0, 1],
            ], dtype=complex)
        else:
            raise ValueError(f"Unknown gate: {gate_name}")

        # Build full unitary for multi-qubit system
        if len(qubits) == 1 and n == 1:
            full_U = U
        elif len(qubits) == 1 and n > 1:
            # Embed single-qubit gate into n-qubit space
            full_U = _embed_1q_gate(U, qubits[0], n)
        elif len(qubits) == 2 and n == 2:
            if qubits == [0, 1]:
                full_U = U
            else:
                # Swap qubits for reversed order
                SWAP = np.array([
                    [1, 0, 0, 0],
                    [0, 0, 1, 0],
                    [0, 1, 0, 0],
                    [0, 0, 0, 1],
                ], dtype=complex)
                full_U = SWAP @ U @ SWAP
        elif len(qubits) == 2 and n > 2:
            full_U = _embed_2q_gate(U, qubits[0], qubits[1], n)
        else:
            raise ValueError(
                f"Gate {gate_name} on {len(qubits)} qubits not supported"
            )

        # Apply unitary
        rho = full_U @ rho @ full_U.conj().T

        # Apply noise (single-qubit noise after each gate)
        if noise is not None and len(qubits) == 1:
            if n == 1:
                rho = noise.apply(rho)
            else:
                rho = _apply_1q_noise(noise, rho, qubits[0], n)

    # Measurement sampling
    # Determine which qubits to measure and in which order
    if not meas_map:
        # No explicit measurements, measure all qubits
        meas_map = {i: i for i in range(n)}

    num_clbits = circuit.num_clbits or n
    probs = np.real(np.diag(rho))
    probs = np.maximum(probs, 0)  # numerical safety
    probs /= probs.sum()

    # Sample outcomes
    outcomes = rng.choice(dim, size=shots, p=probs)

    counts: Dict[str, int] = {}
    for outcome in outcomes:
        # Build bitstring from classical bit mapping
        bits = format(outcome, f"0{n}b")
        # bits[0] is qubit n-1 (MSB), bits[-1] is qubit 0 (LSB)

        clbits = ["0"] * num_clbits
        for clbit_idx, qubit_idx in meas_map.items():
            # Extract bit for this qubit
            bit_val = bits[n - 1 - qubit_idx]

            # Apply readout error
            if noise is not None and noise.readout_error > 0:
                if rng.random() < noise.readout_error:
                    bit_val = "1" if bit_val == "0" else "0"

            clbits[num_clbits - 1 - clbit_idx] = bit_val

        bitstring = "".join(clbits)
        counts[bitstring] = counts.get(bitstring, 0) + 1

    return counts


def _embed_1q_gate(U: NDArray, qubit: int, n: int) -> NDArray:
    """Embed a 1-qubit gate into an n-qubit Hilbert space."""
    ops = [_I] * n
    ops[qubit] = U
    result = ops[0]
    for op in ops[1:]:
        result = np.kron(result, op)
    return result


def _embed_2q_gate(
    U: NDArray, q0: int, q1: int, n: int
) -> NDArray:
    """Embed a 2-qubit gate into an n-qubit Hilbert space.

    The gate acts on qubits q0 (control) and q1 (target), using
    permutation matrices to move them into position.
    """
    dim = 2 ** n
    # Build permutation that maps (q0, q1) to (0, 1)
    qubit_order = list(range(n))
    # Move q0 to position 0 and q1 to position 1
    qubit_order.remove(q0)
    qubit_order.remove(q1)
    perm = [q0, q1] + qubit_order

    # Build permutation matrix
    P = np.zeros((dim, dim), dtype=complex)
    for i in range(dim):
        bits = [(i >> (n - 1 - k)) & 1 for k in range(n)]
        new_bits = [bits[perm[k]] for k in range(n)]
        j = sum(b << (n - 1 - k) for k, b in enumerate(new_bits))
        P[j, i] = 1.0

    # Gate on first two qubits, identity on rest
    gate_full = U
    for _ in range(n - 2):
        gate_full = np.kron(gate_full, _I)

    return P.conj().T @ gate_full @ P


def _apply_1q_noise(
    noise: NoiseModel, rho: NDArray, qubit: int, n: int
) -> NDArray:
    """Apply single-qubit noise channel to a specific qubit in n-qubit state."""
    dim = 2 ** n
    # Partial trace / apply noise using Kraus-style approach
    # For efficiency, we build the 1-qubit reduced state, apply noise,
    # then embed back.  However, this loses correlations.  For a proper
    # simulation we apply the noise Kraus operators directly.

    if noise.depolarizing_rate > 0:
        p = noise.depolarizing_rate
        for pauli in [_X, _Y, _Z]:
            K = np.sqrt(p / 3) * pauli
            K_full = _embed_1q_gate(K, qubit, n)
            rho = (1 - p) * rho + K_full @ rho @ K_full.conj().T
            # Reset for next Pauli -- we need the original rho
        # Correct: apply all at once
        rho_new = (1 - noise.depolarizing_rate) * rho
        for pauli in [_X, _Y, _Z]:
            K_full = _embed_1q_gate(pauli, qubit, n)
            rho_new = rho_new + (noise.depolarizing_rate / 3) * (
                K_full @ rho @ K_full.conj().T
            )
        rho = rho_new

    if noise.amplitude_damping_rate > 0:
        gamma = noise.amplitude_damping_rate
        K0 = np.array([[1, 0], [0, np.sqrt(1 - gamma)]], dtype=complex)
        K1 = np.array([[0, np.sqrt(gamma)], [0, 0]], dtype=complex)
        K0_full = _embed_1q_gate(K0, qubit, n)
        K1_full = _embed_1q_gate(K1, qubit, n)
        rho = K0_full @ rho @ K0_full.conj().T + K1_full @ rho @ K1_full.conj().T

    if noise.dephasing_rate > 0:
        p = noise.dephasing_rate
        Z_full = _embed_1q_gate(_Z, qubit, n)
        rho = (1 - p) * rho + p * Z_full @ rho @ Z_full.conj().T

    return rho


# -----------------------------------------------------------------------
# SimulatorBackend class
# -----------------------------------------------------------------------

class SimulatorBackend:
    """Local numpy-based quantum simulator backend.

    Satisfies the ``TauChronoBackend`` protocol.  No external
    dependencies beyond numpy.

    Parameters
    ----------
    num_qubits : int
        Number of qubits to simulate (default 1).
    noise_model : NoiseModel, optional
        Noise model to apply during simulation.
    seed : int, optional
        Random seed for reproducibility.

    Examples
    --------
    >>> from tau_chrono.backends.simulator import SimulatorBackend, NoiseModel
    >>> backend = SimulatorBackend(num_qubits=1)
    >>> # Use with tomography pipeline:
    >>> from tau_chrono.tomography import simulate_process_tomography_1q
    >>> # Or run raw circuits:
    >>> from tau_chrono.backends.simulator import SimCircuit
    >>> c = SimCircuit(num_qubits=1, num_clbits=1)
    >>> c.add_gate("h", [0])
    >>> c.add_measurement(0, 0)
    >>> counts = backend.run_circuits([c], shots=1000)
    """

    def __init__(
        self,
        num_qubits: int = 1,
        noise_model: Optional[NoiseModel] = None,
        seed: Optional[int] = None,
    ):
        self._name = "tau-chrono-simulator"
        self._num_qubits = num_qubits
        self._noise_model = noise_model
        self._rng = np.random.default_rng(seed)

    @property
    def name(self) -> str:
        return self._name

    @property
    def num_qubits(self) -> int:
        return self._num_qubits

    @property
    def noise_model(self) -> Optional[NoiseModel]:
        return self._noise_model

    @noise_model.setter
    def noise_model(self, model: Optional[NoiseModel]) -> None:
        self._noise_model = model

    def info(self) -> BackendInfo:
        """Return backend metadata."""
        return BackendInfo(
            name=self._name,
            num_qubits=self._num_qubits,
            provider="tau-chrono (local)",
            is_simulator=True,
            max_shots=10_000_000,
            basis_gates=["id", "x", "y", "z", "h", "s", "sdg",
                         "t", "tdg", "rx", "ry", "rz", "cx", "cz"],
        )

    def run_circuits(
        self,
        circuits: Sequence,
        shots: int = 1024,
    ) -> List[Dict[str, int]]:
        """Execute circuits and return measurement counts.

        Accepts either ``SimCircuit`` objects or Qiskit
        ``QuantumCircuit`` objects (auto-detected and converted).

        Parameters
        ----------
        circuits : sequence
            Circuits to execute.
        shots : int
            Shots per circuit.

        Returns
        -------
        list of dict
            Measurement counts for each circuit.
        """
        results = []
        for circuit in circuits:
            # Auto-detect Qiskit circuits and convert
            if not isinstance(circuit, SimCircuit):
                circuit = _parse_qiskit_circuit(circuit)

            counts = _simulate_circuit(
                circuit,
                noise=self._noise_model,
                shots=shots,
                rng=self._rng,
            )
            results.append(counts)
        return results

    def run_tomography_1q(
        self,
        kraus_ops: List[NDArray],
        shots: int = 100_000,
        verbose: bool = True,
    ) -> List[NDArray]:
        """Run simulated 1-qubit process tomography using Kraus operators.

        This is a convenience method that wraps
        ``simulate_process_tomography_1q`` from the tomography module.

        Parameters
        ----------
        kraus_ops : list of NDArray
            The "true" Kraus operators to simulate.
        shots : int
            Simulated shots per circuit.
        verbose : bool
            Print progress.

        Returns
        -------
        list of NDArray
            Reconstructed Kraus operators.
        """
        from ..tomography import simulate_process_tomography_1q

        return simulate_process_tomography_1q(
            kraus_ops, shots=shots, verbose=verbose
        )

    def __repr__(self) -> str:
        noise_str = ""
        if self._noise_model is not None:
            parts = []
            if self._noise_model.depolarizing_rate > 0:
                parts.append(f"depol={self._noise_model.depolarizing_rate:.3f}")
            if self._noise_model.amplitude_damping_rate > 0:
                parts.append(f"ad={self._noise_model.amplitude_damping_rate:.3f}")
            if self._noise_model.dephasing_rate > 0:
                parts.append(f"deph={self._noise_model.dephasing_rate:.3f}")
            if self._noise_model.readout_error > 0:
                parts.append(f"ro={self._noise_model.readout_error:.3f}")
            if parts:
                noise_str = f", noise=[{', '.join(parts)}]"
        return f"SimulatorBackend(qubits={self._num_qubits}{noise_str})"
