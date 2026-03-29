#!/usr/bin/env python3
"""
QEC Experiment: 3-Qubit Repetition Code on QuTech Tuna-9
=========================================================

Runs a 3-qubit bit-flip repetition code (5 qubits: 3 data + 2 ancilla)
on the QuTech Tuna-9 superconducting processor (9 qubits) via
Quantum Inspire.

Three decoding strategies are compared:

  Strategy A -- No error correction:
    Encode |psi> -> |psi psi psi>, wait, measure q0 directly.
    This is the raw (uncorrected) logical error rate.

  Strategy B -- Standard syndrome correction:
    Encode, extract parity syndromes from two ancilla qubits,
    apply majority-vote correction assuming uniform noise.

  Strategy C -- tau-chrono noise-informed correction:
    Same circuit as B, but the majority vote is WEIGHTED by
    per-qubit tau values obtained from tau-chrono process tomography.
    A qubit with higher tau (more noise) is trusted less.

The experiment measures logical error rate for each strategy over
4096 shots and compares them.

Usage:
  python qec_tuna9.py                   # Run on Tuna-9 hardware
  python qec_tuna9.py --local           # Local noisy simulation
  python qec_tuna9.py --shots 8192      # More shots for better stats

Output:
  Prints summary table + saves JSON to results/qec_tuna9_results.json

Requirements:
  pip install qiskit quantuminspire qiskit-quantuminspire
  qi login

Author: Sheng-Kai Huang / QDA
Date: 2026-03-19
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from dataclasses import dataclass, asdict
from typing import Dict, List, Optional, Tuple

import numpy as np

# ---------------------------------------------------------------------------
# tau-chrono imports
# ---------------------------------------------------------------------------
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import tau_chrono as tc
from tau_chrono.petz import tau_parameter, apply_channel, fidelity
from tau_chrono.tomography import simulate_process_tomography_1q

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

BACKEND_NAME = "Tuna-9"
DEFAULT_SHOTS = 4096

# Local simulation noise parameters (5% depolarizing per CNOT, as specified)
LOCAL_1Q_ERROR = 0.01      # 1% depolarizing per 1Q gate
LOCAL_2Q_ERROR = 0.05      # 5% depolarizing per 2Q gate (CNOT)
LOCAL_IDLE_ERROR = 0.005   # idle depolarizing during barrier/wait

# Data qubit and ancilla qubit indices.
# These will be remapped to physical qubits after checking connectivity.
DATA_QUBITS = [0, 1, 2]
ANCILLA_QUBITS = [3, 4]


# ---------------------------------------------------------------------------
# Result dataclasses
# ---------------------------------------------------------------------------

@dataclass
class QubitTauCharacterization:
    """tau-chrono characterization for a single data qubit."""
    qubit_index: int
    physical_qubit: int
    tau_worst: float
    tau_avg: float
    process_fidelity: float


@dataclass
class StrategyResult:
    """Result of running one decoding strategy."""
    name: str
    description: str
    logical_error_rate: float
    total_shots: int
    num_errors: int
    raw_counts: Optional[Dict[str, int]] = None


@dataclass
class ExperimentResult:
    """Full experiment result."""
    backend: str
    timestamp: str
    shots: int
    qubit_mapping: Dict[str, int]
    tau_characterizations: List[dict]
    strategy_results: List[dict]
    improvement_B_over_A: float
    improvement_C_over_B: float
    improvement_C_over_A: float


# ---------------------------------------------------------------------------
# Backend status check
# ---------------------------------------------------------------------------

def check_backend_status(backend_name: str) -> Tuple[str, Optional[object]]:
    """Check whether the QI backend is available and idle.

    Returns
    -------
    (status, backend_obj)
        status is one of: "idle", "busy", "offline", "not_found", "error"
        backend_obj is the raw QI backend if found, else None.
    """
    try:
        from qiskit_quantuminspire.qi_provider import QIProvider
        provider = QIProvider()
        for b in provider.backends():
            if b.name == backend_name:
                status_val = getattr(b, 'status', None)
                if status_val is not None:
                    status_str = (status_val.value
                                  if hasattr(status_val, 'value')
                                  else str(status_val))
                else:
                    # Fallback: try operational flag
                    try:
                        bstatus = b.status()
                        if bstatus.operational:
                            status_str = "idle"
                        else:
                            status_str = "offline"
                    except Exception:
                        status_str = "unknown"
                return status_str, b
        return "not_found", None
    except ImportError:
        return "error", None
    except Exception as e:
        return f"error: {e}", None


# ---------------------------------------------------------------------------
# Connectivity check
# ---------------------------------------------------------------------------

def get_qubit_mapping(backend) -> Dict[str, int]:
    """Determine optimal physical qubit mapping for the repetition code.

    The circuit needs:
      q0 -- CNOT -- q1 -- CNOT -- q2   (data qubits, linear chain)
      q0 -- CNOT -- q3                  (ancilla for parity(q0,q1))
      q1 -- CNOT -- q3
      q1 -- CNOT -- q4                  (ancilla for parity(q1,q2))
      q2 -- CNOT -- q4

    Required connectivity edges:
      (q0,q1), (q0,q3), (q1,q2), (q1,q3), (q1,q4), (q2,q4)

    Returns a mapping: {"d0": phys, "d1": phys, "d2": phys, "a0": phys, "a1": phys}
    """
    # Try to get coupling map from backend
    coupling_map = None
    try:
        config = backend.configuration()
        coupling_map = config.coupling_map
    except Exception:
        pass

    if coupling_map is None:
        try:
            coupling_map = backend.coupling_map
            if hasattr(coupling_map, 'get_edges'):
                coupling_map = list(coupling_map.get_edges())
        except Exception:
            pass

    if coupling_map is None:
        # Fallback: assume linear connectivity for Tuna-9
        # 0-1-2-3-4-5-6-7-8
        print("  WARNING: Could not read coupling map. Assuming linear 0-1-2-3-4-5-6-7-8.")
        coupling_map = [[i, i+1] for i in range(8)] + [[i+1, i] for i in range(8)]

    # Build adjacency set (undirected)
    edges = set()
    for pair in coupling_map:
        a, b = int(pair[0]), int(pair[1])
        edges.add((a, b))
        edges.add((b, a))

    def connected(a, b):
        return (a, b) in edges

    # Get number of qubits
    n_qubits = 9
    try:
        config = backend.configuration()
        n_qubits = config.n_qubits
    except Exception:
        pass

    # Brute-force search over all 5-qubit assignments
    # (for 9 qubits this is C(9,5)*5! = 126*120 = 15120, very fast)
    best_mapping = None

    from itertools import permutations, combinations
    for combo in combinations(range(n_qubits), 5):
        for perm in permutations(combo):
            d0, d1, d2, a0, a1 = perm
            # Check all required edges
            if (connected(d0, d1) and connected(d1, d2) and
                connected(d0, a0) and connected(d1, a0) and
                connected(d1, a1) and connected(d2, a1)):
                best_mapping = {
                    "d0": d0, "d1": d1, "d2": d2,
                    "a0": a0, "a1": a1,
                }
                break
        if best_mapping is not None:
            break

    if best_mapping is None:
        # Relax: allow transpiler to handle routing.
        # Use a default linear assignment and rely on SWAP insertion.
        print("  WARNING: No 5-qubit subset satisfies all connectivity constraints.")
        print("           Using default assignment [0,1,2,3,4]; transpiler will add SWAPs.")
        best_mapping = {"d0": 0, "d1": 1, "d2": 2, "a0": 3, "a1": 4}

    return best_mapping


# ---------------------------------------------------------------------------
# Circuit construction
# ---------------------------------------------------------------------------

def build_circuit_strategy_A(qubit_map: Dict[str, int]) -> "QuantumCircuit":
    """Strategy A: Encode, wait, measure q0 only (no correction).

    This measures the raw logical error rate without any syndrome extraction
    or correction. We prepare |0> -> |000>, let noise accumulate, and read
    out just the first data qubit.
    """
    from qiskit import QuantumCircuit

    n_phys = max(qubit_map.values()) + 1
    qc = QuantumCircuit(n_phys, 1, name="Strategy_A_no_correction")

    d0, d1, d2 = qubit_map["d0"], qubit_map["d1"], qubit_map["d2"]

    # Encoding: |0> -> |000>
    qc.cx(d0, d1)
    qc.cx(d0, d2)

    # Noise accumulation barrier
    qc.barrier()

    # Add identity-equivalent gates to increase circuit depth (more noise)
    # Two pairs of X gates = identity, but the gate execution adds noise
    for q in [d0, d1, d2]:
        qc.id(q)
        qc.id(q)

    qc.barrier()

    # Measure d0 only (no syndrome, no correction)
    qc.measure(d0, 0)

    return qc


def build_circuit_strategy_B(qubit_map: Dict[str, int]) -> "QuantumCircuit":
    """Strategy B: Encode, wait, syndrome extraction, measure all 5 qubits.

    Syndrome is extracted into ancilla qubits:
      a0 = parity(d0, d1)
      a1 = parity(d1, d2)

    All 5 qubits are measured. Classical post-processing applies
    majority-vote correction.
    """
    from qiskit import QuantumCircuit

    n_phys = max(qubit_map.values()) + 1
    qc = QuantumCircuit(n_phys, 5, name="Strategy_B_standard_correction")

    d0, d1, d2 = qubit_map["d0"], qubit_map["d1"], qubit_map["d2"]
    a0, a1 = qubit_map["a0"], qubit_map["a1"]

    # Encoding: |0> -> |000>
    qc.cx(d0, d1)
    qc.cx(d0, d2)

    # Noise accumulation barrier
    qc.barrier()

    # Add identity gates for noise accumulation
    for q in [d0, d1, d2]:
        qc.id(q)
        qc.id(q)

    qc.barrier()

    # Syndrome extraction
    # a0 measures parity of (d0, d1)
    qc.cx(d0, a0)
    qc.cx(d1, a0)

    # a1 measures parity of (d1, d2)
    qc.cx(d1, a1)
    qc.cx(d2, a1)

    qc.barrier()

    # Measure all: data qubits on classical bits 0-2, ancilla on 3-4
    qc.measure(d0, 0)
    qc.measure(d1, 1)
    qc.measure(d2, 2)
    qc.measure(a0, 3)
    qc.measure(a1, 4)

    return qc


# Strategy C uses the same circuit as B; the difference is in post-processing.
# So build_circuit_strategy_C = build_circuit_strategy_B.


# ---------------------------------------------------------------------------
# Classical decoding
# ---------------------------------------------------------------------------

def decode_strategy_A(counts: Dict[str, int]) -> Tuple[float, int, int]:
    """Decode Strategy A: raw readout of d0.

    Returns (logical_error_rate, num_errors, total_shots).
    We encoded |0>, so any measurement of '1' on d0 is an error.
    """
    total = sum(counts.values())
    errors = 0
    for bitstring, count in counts.items():
        # Qiskit little-endian: rightmost bit = classical bit 0 = d0
        bit_d0 = int(bitstring[-1])
        if bit_d0 == 1:
            errors += count
    return errors / total, errors, total


def decode_strategy_B(counts: Dict[str, int]) -> Tuple[float, int, int]:
    """Decode Strategy B: standard majority-vote correction.

    Syndrome table for 3-qubit repetition code:
      s0 s1 | error location | correction
      0  0  | no error       | none
      1  0  | d0             | flip d0
      1  1  | d1             | flip d1
      0  1  | d2             | flip d2

    After correction, the logical value is the corrected d0.
    We encoded |0>, so corrected d0 = 1 is a logical error.
    """
    total = sum(counts.values())
    errors = 0

    for bitstring, count in counts.items():
        # Qiskit little-endian: bit positions from right
        # bit 0 = d0, bit 1 = d1, bit 2 = d2, bit 3 = a0, bit 4 = a1
        bits = bitstring.zfill(5)  # Ensure 5 characters
        # Reverse because Qiskit is little-endian
        d0 = int(bits[-(0+1)])  # classical bit 0
        d1 = int(bits[-(1+1)])  # classical bit 1
        d2 = int(bits[-(2+1)])  # classical bit 2
        s0 = int(bits[-(3+1)])  # classical bit 3 = a0 = parity(d0,d1)
        s1 = int(bits[-(4+1)])  # classical bit 4 = a1 = parity(d1,d2)

        # Apply syndrome-based correction
        corrected_d0 = d0
        corrected_d1 = d1
        corrected_d2 = d2

        if s0 == 1 and s1 == 0:
            # Error on d0
            corrected_d0 ^= 1
        elif s0 == 1 and s1 == 1:
            # Error on d1
            corrected_d1 ^= 1
        elif s0 == 0 and s1 == 1:
            # Error on d2
            corrected_d2 ^= 1
        # s0==0, s1==0: no single-qubit error detected

        # Logical value = corrected d0 (we could also use majority vote
        # of corrected data, but for the repetition code they should agree)
        if corrected_d0 != 0:
            errors += count

    return errors / total, errors, total


def decode_strategy_C(
    counts: Dict[str, int],
    tau_values: Dict[int, float],
    qubit_map: Dict[str, int],
) -> Tuple[float, int, int]:
    """Decode Strategy C: tau-weighted correction.

    Instead of the fixed syndrome lookup, we use per-qubit tau values
    to weight the correction decision. A qubit with higher tau is
    deemed noisier and trusted less.

    The decoding logic:
      1. Extract raw data bits and syndrome bits.
      2. If syndrome indicates an error, decide WHICH qubit to correct
         using tau-weighted trust: the qubit with the HIGHEST tau
         (most noise) among the syndrome-implicated qubits is corrected.
      3. For ambiguous syndromes (e.g., multi-qubit errors), use
         weighted majority vote on all 3 data qubits where the vote
         weight for qubit i is (1 - tau_i).

    Parameters
    ----------
    counts : dict
        Measurement counts from the 5-qubit circuit.
    tau_values : dict
        Mapping from physical qubit index to tau_worst value.
    qubit_map : dict
        Logical-to-physical qubit mapping.
    """
    total = sum(counts.values())
    errors = 0

    # Get tau for each data qubit
    tau_d0 = tau_values.get(qubit_map["d0"], 0.1)
    tau_d1 = tau_values.get(qubit_map["d1"], 0.1)
    tau_d2 = tau_values.get(qubit_map["d2"], 0.1)

    # Trust weights: lower tau = more trustworthy
    # w_i = 1 - tau_i (higher weight = more trusted)
    w0 = 1.0 - tau_d0
    w1 = 1.0 - tau_d1
    w2 = 1.0 - tau_d2

    for bitstring, count in counts.items():
        bits = bitstring.zfill(5)
        d0 = int(bits[-(0+1)])
        d1 = int(bits[-(1+1)])
        d2 = int(bits[-(2+1)])
        s0 = int(bits[-(3+1)])
        s1 = int(bits[-(4+1)])

        corrected_d0 = d0
        corrected_d1 = d1
        corrected_d2 = d2

        if s0 == 0 and s1 == 0:
            # No error detected by syndrome.
            # But we can still do a weighted majority vote in case
            # the syndrome itself had errors.
            pass

        elif s0 == 1 and s1 == 0:
            # Syndrome implicates d0 or the ancilla.
            # Standard: flip d0. Tau-informed: flip the NOISIER of d0/d1.
            # (The syndrome says d0 and d1 disagree; flip the one we trust less)
            if tau_d0 >= tau_d1:
                corrected_d0 ^= 1  # d0 is noisier, likely flipped
            else:
                corrected_d1 ^= 1  # d1 is noisier

        elif s0 == 0 and s1 == 1:
            # Syndrome implicates d2 or d1.
            if tau_d2 >= tau_d1:
                corrected_d2 ^= 1
            else:
                corrected_d1 ^= 1

        elif s0 == 1 and s1 == 1:
            # Syndrome implicates d1 (standard interpretation).
            # But could also be two errors on d0 and d2.
            # Tau-informed: if d1 is the noisiest, flip d1 (standard).
            # If d1 is actually the cleanest, consider the two-error case.
            if tau_d1 >= tau_d0 and tau_d1 >= tau_d2:
                corrected_d1 ^= 1
            else:
                # Two-error scenario: flip d0 and d2 if they are both noisier
                # than d1. Otherwise fall back to standard (flip d1).
                if tau_d0 + tau_d2 > 2 * tau_d1:
                    corrected_d0 ^= 1
                    corrected_d2 ^= 1
                else:
                    corrected_d1 ^= 1

        # After correction, determine logical value via weighted majority vote
        # Vote: each qubit votes for its corrected value, weighted by trust
        vote_0 = w0 * (1 - corrected_d0) + w1 * (1 - corrected_d1) + w2 * (1 - corrected_d2)
        vote_1 = w0 * corrected_d0 + w1 * corrected_d1 + w2 * corrected_d2

        logical = 0 if vote_0 >= vote_1 else 1

        # We encoded |0>, so logical=1 is an error
        if logical != 0:
            errors += count

    return errors / total, errors, total


# ---------------------------------------------------------------------------
# tau-chrono qubit characterization
# ---------------------------------------------------------------------------

def characterize_qubits_hardware(
    backend,
    qubit_map: Dict[str, int],
    shots: int = 4096,
) -> Dict[int, QubitTauCharacterization]:
    """Characterize each data qubit using tau-chrono process tomography.

    Runs an identity gate (do-nothing) on each data qubit and measures
    tau, which captures the idle noise on that qubit.

    Parameters
    ----------
    backend
        Quantum Inspire backend (QuantumInspireBackend wrapper).
    qubit_map
        Logical-to-physical qubit mapping.
    shots
        Shots per tomography circuit.

    Returns
    -------
    dict mapping physical qubit index -> QubitTauCharacterization
    """
    from tau_chrono.backends.quantum_inspire import extract_gate_kraus

    results = {}
    sigma = np.eye(2, dtype=complex) / 2
    test_states = [
        np.array([[1, 0], [0, 0]], dtype=complex),
        np.array([[0, 0], [0, 1]], dtype=complex),
        np.array([[0.5, 0.5], [0.5, 0.5]], dtype=complex),
    ]

    for logical_name in ["d0", "d1", "d2"]:
        phys_q = qubit_map[logical_name]
        logical_idx = int(logical_name[1])

        print(f"  Characterizing {logical_name} (physical qubit {phys_q})...")

        # Identity gate = idle noise characterization
        def gate_id(qc, q):
            qc.id(q)

        char = extract_gate_kraus(
            gate_id, backend, f"id_q{phys_q}",
            shots=shots, verbose=False,
        )

        results[phys_q] = QubitTauCharacterization(
            qubit_index=logical_idx,
            physical_qubit=phys_q,
            tau_worst=char.tau_worst,
            tau_avg=char.tau_avg,
            process_fidelity=char.process_fidelity,
        )
        print(f"    tau_worst = {char.tau_worst:.4f}, "
              f"tau_avg = {char.tau_avg:.4f}, "
              f"F_proc = {char.process_fidelity:.4f}")

    return results


def characterize_qubits_local(
    qubit_map: Dict[str, int],
    noise_profile: Optional[Dict[str, float]] = None,
) -> Dict[int, QubitTauCharacterization]:
    """Simulate per-qubit characterization with non-uniform noise.

    Creates intentionally non-uniform noise to demonstrate the advantage
    of noise-informed correction.

    Parameters
    ----------
    qubit_map
        Logical-to-physical qubit mapping.
    noise_profile
        Optional dict with keys "d0", "d1", "d2" mapping to depolarizing
        error rates. If None, uses a default non-uniform profile.
    """
    if noise_profile is None:
        # Non-uniform noise: d1 is the noisiest qubit
        noise_profile = {
            "d0": 0.02,   # 2% depolarizing
            "d1": 0.08,   # 8% depolarizing (noisy qubit!)
            "d2": 0.03,   # 3% depolarizing
        }

    sigma = np.eye(2, dtype=complex) / 2
    test_states = [
        np.array([[1, 0], [0, 0]], dtype=complex),
        np.array([[0, 0], [0, 1]], dtype=complex),
        np.array([[0.5, 0.5], [0.5, 0.5]], dtype=complex),
    ]

    results = {}
    for logical_name in ["d0", "d1", "d2"]:
        phys_q = qubit_map[logical_name]
        logical_idx = int(logical_name[1])
        p = noise_profile[logical_name]

        kraus = tc.depolarizing(p)

        taus = [tau_parameter(rho, kraus, sigma) for rho in test_states]
        tau_worst = max(taus)
        tau_avg = sum(taus) / len(taus)

        fids = [fidelity(rho, apply_channel(rho, kraus)) for rho in test_states]
        f_proc = (2 * (sum(fids) / len(fids)) + 1) / 3

        results[phys_q] = QubitTauCharacterization(
            qubit_index=logical_idx,
            physical_qubit=phys_q,
            tau_worst=tau_worst,
            tau_avg=tau_avg,
            process_fidelity=f_proc,
        )
        print(f"  {logical_name} (phys {phys_q}): p_depol = {p:.3f}, "
              f"tau_worst = {tau_worst:.4f}, tau_avg = {tau_avg:.4f}")

    return results


# ---------------------------------------------------------------------------
# Local noisy simulation
# ---------------------------------------------------------------------------

def simulate_circuit_local(
    circuit_name: str,
    qubit_map: Dict[str, int],
    shots: int,
    noise_profile: Optional[Dict[str, float]] = None,
) -> Dict[str, int]:
    """Simulate a QEC circuit with realistic depolarizing noise.

    Instead of using Qiskit Aer (which may not be installed), this
    performs a direct density-matrix simulation with noise injection
    after each gate layer.

    Parameters
    ----------
    circuit_name
        "A" or "B" (C uses the same circuit as B).
    qubit_map
        Logical-to-physical qubit mapping.
    shots
        Number of measurement shots.
    noise_profile
        Per-qubit depolarizing rates. If None, uses default non-uniform.

    Returns
    -------
    dict
        Measurement counts (bitstrings -> counts).
    """
    if noise_profile is None:
        noise_profile = {
            "d0": 0.02,
            "d1": 0.08,
            "d2": 0.03,
        }

    # We simulate by tracking the probability of each computational basis state.
    # For 5 qubits, 2^5 = 32 states -- tractable.
    if circuit_name == "A":
        n_bits = 1
        n_qubits = 3  # Only data qubits matter
    else:
        n_bits = 5
        n_qubits = 5  # 3 data + 2 ancilla

    # Use a simplified stochastic simulation:
    # For each shot, track the 3 data qubits + 2 ancilla as classical bits,
    # applying probabilistic bit-flip errors.
    rng = np.random.default_rng(seed=None)

    # Error rates per qubit
    p_d0 = noise_profile.get("d0", 0.02)
    p_d1 = noise_profile.get("d1", 0.08)
    p_d2 = noise_profile.get("d2", 0.03)
    p_cnot = LOCAL_2Q_ERROR   # 5% per CNOT
    p_ancilla = 0.01          # ancilla measurement noise

    counts = {}

    for _ in range(shots):
        # Initialize |00000>
        d = [0, 0, 0, 0, 0]  # d0, d1, d2, a0, a1

        # ENCODING: CNOT(d0, d1), CNOT(d0, d2)
        # CNOT flips target if control is 1
        # Since we start in |0>, these are no-ops on the state.
        # But they inject CNOT noise.
        # CNOT error model: with probability p_cnot, apply random Pauli
        if rng.random() < p_cnot:
            error_type = rng.integers(0, 3)
            if error_type == 0:    # X on control
                d[0] ^= 1
            elif error_type == 1:  # X on target
                d[1] ^= 1
            else:                  # X on both
                d[0] ^= 1
                d[1] ^= 1

        if rng.random() < p_cnot:
            error_type = rng.integers(0, 3)
            if error_type == 0:
                d[0] ^= 1
            elif error_type == 1:
                d[2] ^= 1
            else:
                d[0] ^= 1
                d[2] ^= 1

        # IDLE NOISE (barrier + id gates)
        # Each data qubit gets independent depolarizing noise
        if rng.random() < p_d0:
            d[0] ^= 1
        if rng.random() < p_d1:
            d[1] ^= 1
        if rng.random() < p_d2:
            d[2] ^= 1

        if circuit_name == "A":
            # Strategy A: just measure d0
            # Measurement noise
            meas_d0 = d[0]
            if rng.random() < 0.005:
                meas_d0 ^= 1
            result = str(meas_d0)
        else:
            # SYNDROME EXTRACTION
            # CNOT(d0, a0): a0 ^= d0
            d[3] ^= d[0]
            if rng.random() < p_cnot:
                error_type = rng.integers(0, 3)
                if error_type == 0:
                    d[0] ^= 1
                elif error_type == 1:
                    d[3] ^= 1
                else:
                    d[0] ^= 1
                    d[3] ^= 1

            # CNOT(d1, a0): a0 ^= d1
            d[3] ^= d[1]
            if rng.random() < p_cnot:
                error_type = rng.integers(0, 3)
                if error_type == 0:
                    d[1] ^= 1
                elif error_type == 1:
                    d[3] ^= 1
                else:
                    d[1] ^= 1
                    d[3] ^= 1

            # CNOT(d1, a1): a1 ^= d1
            d[4] ^= d[1]
            if rng.random() < p_cnot:
                error_type = rng.integers(0, 3)
                if error_type == 0:
                    d[1] ^= 1
                elif error_type == 1:
                    d[4] ^= 1
                else:
                    d[1] ^= 1
                    d[4] ^= 1

            # CNOT(d2, a1): a1 ^= d2
            d[4] ^= d[2]
            if rng.random() < p_cnot:
                error_type = rng.integers(0, 3)
                if error_type == 0:
                    d[2] ^= 1
                elif error_type == 1:
                    d[4] ^= 1
                else:
                    d[2] ^= 1
                    d[4] ^= 1

            # Measurement noise on all 5 qubits
            measured = list(d)
            for i in range(5):
                if rng.random() < 0.005:
                    measured[i] ^= 1

            # Build bitstring: Qiskit little-endian convention
            # classical bit 0 = d0 (rightmost), ..., bit 4 = a1 (leftmost)
            result = ""
            for i in reversed(range(5)):
                result += str(measured[i])

        counts[result] = counts.get(result, 0) + 1

    return counts


# ---------------------------------------------------------------------------
# Main experiment
# ---------------------------------------------------------------------------

def run_experiment(
    local: bool = False,
    shots: int = DEFAULT_SHOTS,
    backend_name: str = BACKEND_NAME,
) -> ExperimentResult:
    """Run the full 3-strategy QEC experiment.

    Parameters
    ----------
    local : bool
        If True, simulate locally with depolarizing noise.
    shots : int
        Shots per circuit.
    backend_name : str
        QI backend name (ignored if local=True).
    """
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")

    print("=" * 72)
    print("QEC EXPERIMENT: 3-Qubit Repetition Code")
    print("=" * 72)
    print(f"  Mode:      {'LOCAL SIMULATION' if local else f'HARDWARE ({backend_name})'}")
    print(f"  Shots:     {shots:,}")
    print(f"  Timestamp: {timestamp}")
    print("=" * 72)

    backend = None
    raw_backend = None

    if not local:
        # Check backend status
        print(f"\nChecking {backend_name} status...")
        status, raw_be = check_backend_status(backend_name)
        print(f"  Status: {status}")

        if "idle" not in str(status).lower() and "unknown" not in str(status).lower():
            print(f"\n  {backend_name} is not idle (status: {status}).")
            print(f"  Use --local to run a local simulation instead.")
            print(f"  Or wait and try again later.")
            sys.exit(1)

        # Connect via tau-chrono backend
        from tau_chrono.backends.quantum_inspire import get_qi_backend
        backend = get_qi_backend(backend_name, optimization_level=0)
        raw_backend = backend.raw_backend
        print(f"  Connected: {backend.name} ({backend.num_qubits} qubits)")

    # -----------------------------------------------------------------------
    # Step 1: Determine qubit mapping
    # -----------------------------------------------------------------------
    print(f"\n--- Step 1: Qubit Mapping ---")
    if local:
        qubit_map = {"d0": 0, "d1": 1, "d2": 2, "a0": 3, "a1": 4}
        print(f"  Using default mapping: {qubit_map}")
    else:
        qubit_map = get_qubit_mapping(raw_backend)
        print(f"  Optimal mapping: {qubit_map}")

    # -----------------------------------------------------------------------
    # Step 2: Characterize data qubits with tau-chrono
    # -----------------------------------------------------------------------
    print(f"\n--- Step 2: tau-chrono Qubit Characterization ---")
    if local:
        tau_chars = characterize_qubits_local(qubit_map)
    else:
        tau_chars = characterize_qubits_hardware(backend, qubit_map, shots)

    # Build tau lookup for decoder
    tau_values = {phys: char.tau_worst for phys, char in tau_chars.items()}
    print(f"\n  tau values: { {k: f'{v:.4f}' for k, v in tau_values.items()} }")

    # Identify noisiest qubit
    noisiest_phys = max(tau_values, key=tau_values.get)
    noisiest_logical = [k for k, v in qubit_map.items()
                        if v == noisiest_phys and k.startswith("d")][0]
    print(f"  Noisiest data qubit: {noisiest_logical} (phys {noisiest_phys}, "
          f"tau = {tau_values[noisiest_phys]:.4f})")

    # -----------------------------------------------------------------------
    # Step 3: Run Strategy A (no correction)
    # -----------------------------------------------------------------------
    print(f"\n--- Step 3: Strategy A (No Error Correction) ---")
    t0 = time.time()

    if local:
        counts_A = simulate_circuit_local("A", qubit_map, shots)
    else:
        from qiskit.compiler import transpile
        circ_A = build_circuit_strategy_A(qubit_map)
        transpiled_A = transpile(circ_A, raw_backend, optimization_level=0)
        job_A = raw_backend.run(transpiled_A, shots=shots)
        job_A.wait_for_final_state(timeout=600)
        counts_A = job_A.result().get_counts(0)

    ler_A, err_A, tot_A = decode_strategy_A(counts_A)
    dt_A = time.time() - t0
    print(f"  Logical error rate:  {ler_A:.4f}  ({err_A}/{tot_A} errors)")
    print(f"  Time: {dt_A:.1f}s")

    result_A = StrategyResult(
        name="A",
        description="No error correction (raw readout)",
        logical_error_rate=ler_A,
        total_shots=tot_A,
        num_errors=err_A,
    )

    # -----------------------------------------------------------------------
    # Step 4: Run Strategy B (standard syndrome correction)
    # -----------------------------------------------------------------------
    print(f"\n--- Step 4: Strategy B (Standard Syndrome Correction) ---")
    t0 = time.time()

    if local:
        counts_B = simulate_circuit_local("B", qubit_map, shots)
    else:
        circ_B = build_circuit_strategy_B(qubit_map)
        transpiled_B = transpile(circ_B, raw_backend, optimization_level=0)
        job_B = raw_backend.run(transpiled_B, shots=shots)
        job_B.wait_for_final_state(timeout=600)
        counts_B = job_B.result().get_counts(0)

    ler_B, err_B, tot_B = decode_strategy_B(counts_B)
    dt_B = time.time() - t0
    print(f"  Logical error rate:  {ler_B:.4f}  ({err_B}/{tot_B} errors)")
    print(f"  Time: {dt_B:.1f}s")

    result_B = StrategyResult(
        name="B",
        description="Standard syndrome correction (uniform majority vote)",
        logical_error_rate=ler_B,
        total_shots=tot_B,
        num_errors=err_B,
    )

    # -----------------------------------------------------------------------
    # Step 5: Run Strategy C (tau-informed correction)
    # -----------------------------------------------------------------------
    print(f"\n--- Step 5: Strategy C (tau-chrono Noise-Informed Correction) ---")
    t0 = time.time()

    # Strategy C uses the SAME circuit as B, but different decoding
    if local:
        counts_C = simulate_circuit_local("B", qubit_map, shots)
    else:
        # Re-use the same circuit as B (or run again for independent statistics)
        circ_C = build_circuit_strategy_B(qubit_map)
        circ_C.name = "Strategy_C_tau_informed"
        transpiled_C = transpile(circ_C, raw_backend, optimization_level=0)
        job_C = raw_backend.run(transpiled_C, shots=shots)
        job_C.wait_for_final_state(timeout=600)
        counts_C = job_C.result().get_counts(0)

    ler_C, err_C, tot_C = decode_strategy_C(counts_C, tau_values, qubit_map)
    dt_C = time.time() - t0
    print(f"  Logical error rate:  {ler_C:.4f}  ({err_C}/{tot_C} errors)")
    print(f"  Time: {dt_C:.1f}s")

    result_C = StrategyResult(
        name="C",
        description="tau-chrono noise-informed correction (weighted majority vote)",
        logical_error_rate=ler_C,
        total_shots=tot_C,
        num_errors=err_C,
    )

    # -----------------------------------------------------------------------
    # Step 6: Compare results
    # -----------------------------------------------------------------------
    print(f"\n{'=' * 72}")
    print("RESULTS COMPARISON")
    print(f"{'=' * 72}")

    # Improvement metrics (handle zero denominators)
    def safe_improvement(baseline, improved):
        if baseline <= 0:
            return 0.0
        return (baseline - improved) / baseline * 100

    imp_B_A = safe_improvement(ler_A, ler_B)
    imp_C_B = safe_improvement(ler_B, ler_C)
    imp_C_A = safe_improvement(ler_A, ler_C)

    print(f"\n  {'Strategy':<45} {'LER':>10} {'Errors':>10}")
    print(f"  {'-'*45} {'-'*10} {'-'*10}")
    print(f"  {'A: No correction (raw readout)':<45} {ler_A:>10.4f} {err_A:>10d}")
    print(f"  {'B: Standard syndrome (uniform vote)':<45} {ler_B:>10.4f} {err_B:>10d}")
    print(f"  {'C: tau-chrono informed (weighted vote)':<45} {ler_C:>10.4f} {err_C:>10d}")

    print(f"\n  Improvements:")
    print(f"    B over A (syndrome helps):          {imp_B_A:+.1f}%")
    print(f"    C over B (tau-chrono helps):         {imp_C_B:+.1f}%")
    print(f"    C over A (total improvement):        {imp_C_A:+.1f}%")

    print(f"\n  Per-qubit tau values used for weighting:")
    for logical_name in ["d0", "d1", "d2"]:
        phys = qubit_map[logical_name]
        tau_w = tau_values[phys]
        trust = 1.0 - tau_w
        print(f"    {logical_name} (phys {phys}): tau = {tau_w:.4f}, "
              f"trust weight = {trust:.4f}")

    print(f"\n  INTERPRETATION:")
    if imp_C_B > 0:
        print(f"    tau-chrono noise-informed decoding IMPROVES over standard "
              f"correction by {imp_C_B:.1f}%.")
        print(f"    The noisiest qubit ({noisiest_logical}) was correctly "
              f"down-weighted in the majority vote.")
    elif imp_C_B == 0:
        print(f"    tau-chrono decoding matches standard correction.")
        print(f"    With uniform noise, both strategies are equivalent (expected).")
    else:
        print(f"    In this run, tau-chrono did not improve over standard correction.")
        print(f"    This can happen with small shot counts or nearly uniform noise.")

    # -----------------------------------------------------------------------
    # Build and save result
    # -----------------------------------------------------------------------
    experiment_result = ExperimentResult(
        backend="local_simulation" if local else backend_name,
        timestamp=timestamp,
        shots=shots,
        qubit_mapping=qubit_map,
        tau_characterizations=[asdict(c) for c in tau_chars.values()],
        strategy_results=[asdict(result_A), asdict(result_B), asdict(result_C)],
        improvement_B_over_A=imp_B_A,
        improvement_C_over_B=imp_C_B,
        improvement_C_over_A=imp_C_A,
    )

    # Save JSON
    results_dir = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "results",
    )
    os.makedirs(results_dir, exist_ok=True)
    json_path = os.path.join(results_dir, "qec_tuna9_results.json")

    with open(json_path, "w") as f:
        json.dump(asdict(experiment_result), f, indent=2, default=str)
    print(f"\n  Results saved to: {json_path}")

    return experiment_result


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="QEC 3-Qubit Repetition Code on QuTech Tuna-9",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python qec_tuna9.py                   # Run on Tuna-9 hardware
  python qec_tuna9.py --local           # Local noisy simulation
  python qec_tuna9.py --local --shots 8192
  python qec_tuna9.py --backend "Tuna-5"
        """,
    )
    parser.add_argument(
        "--local", action="store_true",
        help="Run local simulation with 5%% depolarizing noise (no hardware needed)",
    )
    parser.add_argument(
        "--backend", default=BACKEND_NAME,
        help=f"Quantum Inspire backend name (default: {BACKEND_NAME})",
    )
    parser.add_argument(
        "--shots", type=int, default=DEFAULT_SHOTS,
        help=f"Shots per circuit (default: {DEFAULT_SHOTS})",
    )
    parser.add_argument(
        "--status", action="store_true",
        help="Check backend status and exit",
    )

    args = parser.parse_args()

    if args.status:
        status, _ = check_backend_status(args.backend)
        print(f"{args.backend}: {status}")
        if "idle" in str(status).lower():
            print("Backend is ready. You can run the experiment.")
        else:
            print("Backend is not idle. Use --local for simulation.")
        return

    run_experiment(
        local=args.local,
        shots=args.shots,
        backend_name=args.backend,
    )
    print("\nDone.")


if __name__ == "__main__":
    main()
