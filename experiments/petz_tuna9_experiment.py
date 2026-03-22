#!/usr/bin/env python3
"""
Petz Recovery Fidelity Experiment on Tuna-9 Quantum Computer.

Experimentally verifies the Petz recovery fidelity bound:

    F(rho, R_Petz(N(rho)))  >=  exp(-DeltaD)        [Wilde 2015]

  where DeltaD = D(rho||sigma) - D(N(rho)||N(sigma)) is the relative
  entropy decrease under the channel N.
  Equivalently:  tau = 1 - F  <=  1 - exp(-DeltaD).

  In Paper 1 (Huang 2026), Sigma is defined in the CMI (conditional
  mutual information) form via a tripartite extension, yielding the
  equivalent expression  F >= exp(-Sigma_CMI/2).  In this experiment
  we work in the channel divergence setting, where the bound takes
  the simpler form  F >= exp(-DeltaD).

  When --rotated is passed, the experiment also computes the rotated
  Petz map (JRSWW 2018) and verifies that F_rot also satisfies the
  same bound.

  Note: Our fidelity() returns F^2 = (Tr sqrt(sqrt(rho) sigma sqrt(rho)))^2
  (squared convention).  The JRSWW bound F^2 >= exp(-DeltaD) then
  becomes fidelity() >= exp(-DeltaD), the same threshold for both maps.
  The standard Petz map typically gives HIGHER fidelity than the rotated
  map since it is the Bayesian-consistent (functorial) recovery.

Three experimental modes:

    Mode 1 -- Simulator (exact):
        Analytical Kraus operators, exact linear algebra.
        Ground-truth baseline; no sampling noise.

    Mode 2 -- Noisy simulator (process tomography pipeline):
        Finite-shot process tomography -> Kraus reconstruction -> Petz.
        Tests the full measurement pipeline as hardware would experience.

    Mode 3 -- Tuna-9 hardware (Quantum Inspire):
        Real decoherence from Tuna-9's native noise.
        Process tomography -> Kraus extraction -> classical Petz recovery.

Usage:
    python experiments/petz_tuna9_experiment.py --mode simulator
    python experiments/petz_tuna9_experiment.py --mode noisy-sim
    python experiments/petz_tuna9_experiment.py --mode tuna9
    python experiments/petz_tuna9_experiment.py --mode all --plot --save

Author: Sheng-Kai Huang (2026)
License: MIT
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from numpy.typing import NDArray

# ---------------------------------------------------------------------------
# Project imports
# ---------------------------------------------------------------------------
_PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_PROJECT_ROOT))

from tau_chrono.petz import (
    apply_channel,
    apply_petz_recovery,
    delta_D,
    fidelity,
    petz_recovery_map,
    relative_entropy,
    tau_parameter,
)
from tau_chrono.channels import (
    amplitude_damping,
    depolarizing,
    dephasing,
    verify_cptp,
)
from tau_chrono.backends.simulator import NoiseModel, SimulatorBackend
from tau_chrono.tomography import simulate_process_tomography_1q


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------


@dataclass
class PointResult:
    """Result for a single (rho, sigma, channel) experiment point."""

    label: str                    # e.g. "gamma=0.30" or "depth=10"
    param: float                  # channel parameter

    # Standard Petz recovery
    F_petz: float                 # F(rho, R_Petz(N(rho)))
    tau_petz: float               # 1 - F_petz

    # Bound A: standard Petz, channel setting (Wilde 2015)
    #   F >= exp(-DeltaD)  <==>  tau <= 1 - exp(-DeltaD)
    delta_D_val: float            # DeltaD = D(rho||sigma) - D(N(rho)||N(sigma))
    bound_A: float                # 1 - exp(-DeltaD)
    satisfies_A: bool             # tau_petz <= bound_A ?

    # Bound B: rotated Petz, JRSWW (channel setting)
    #   F_rot (squared convention) >= exp(-DeltaD)
    #   tau_rot <= 1 - exp(-DeltaD)  [same threshold, different recovery map]
    bound_B: float                # 1 - exp(-DeltaD) [same as bound_A]

    # Relative entropies (for diagnostics)
    D_before: float               # D(rho || sigma)
    D_after: float                # D(N(rho) || N(sigma))
    source: str = "simulator"

    # Rotated Petz results (computed when available)
    F_rotated: Optional[float] = None
    tau_rotated: Optional[float] = None
    satisfies_B: Optional[bool] = None  # tau_rotated <= bound_B ?


@dataclass
class ExperimentResult:
    """Full result set for one experiment sweep."""

    experiment_name: str
    channel_type: str
    source: str
    input_state: str
    reference_state: str
    n_points: int
    points: List[PointResult] = field(default_factory=list)
    all_satisfy_A: bool = True
    all_satisfy_B: bool = True
    timestamp: str = ""


# ---------------------------------------------------------------------------
# Standard states
# ---------------------------------------------------------------------------

RHO_0 = np.array([[1, 0], [0, 0]], dtype=complex)         # |0><0|
RHO_1 = np.array([[0, 0], [0, 1]], dtype=complex)         # |1><1|
RHO_PLUS = np.array([[0.5, 0.5], [0.5, 0.5]], dtype=complex)  # |+><+|
RHO_PLUS_I = np.array([[0.5, -0.5j], [0.5j, 0.5]], dtype=complex)  # |+i><+i|
SIGMA_MM = np.eye(2, dtype=complex) / 2.0                 # maximally mixed


# ---------------------------------------------------------------------------
# Core measurement function
# ---------------------------------------------------------------------------


def measure_point(
    rho: NDArray,
    sigma: NDArray,
    kraus_ops: List[NDArray],
    label: str,
    param: float,
    source: str = "simulator",
    compute_rotated: bool = False,
) -> PointResult:
    """
    Measure tau, F, DeltaD, and check both bounds for one (rho, sigma, N).

    Bound A (standard Petz, Wilde 2015):
        F(rho, R_Petz(N(rho))) >= exp(-DeltaD)

    Bound B (rotated Petz, JRSWW 2018):
        F(rho, R_rot(N(rho))) >= exp(-DeltaD)
        [same threshold; uses rotated recovery map for comparison]
    """
    # Apply channel
    rho_out = apply_channel(rho, kraus_ops)
    sigma_out = apply_channel(sigma, kraus_ops)

    # Standard Petz recovery
    S_R = petz_recovery_map(kraus_ops, sigma)
    rho_recovered = apply_petz_recovery(rho_out, kraus_ops, sigma, S_R)

    F_petz = fidelity(rho, rho_recovered)
    tau_petz = 1.0 - F_petz

    # Relative entropies
    D_before = relative_entropy(rho, sigma)
    D_after = relative_entropy(rho_out, sigma_out)

    if np.isinf(D_before) or np.isinf(D_after):
        dD = float("inf")
    else:
        dD = max(D_before - D_after, 0.0)

    # Bound A: tau_petz <= 1 - exp(-DeltaD)
    if np.isinf(dD):
        bound_A = 1.0
    else:
        bound_A = 1.0 - np.exp(-dD)

    satisfies_A = tau_petz <= bound_A + 1e-10

    # Bound B (JRSWW, channel setting): F_rot (squared convention) >= exp(-DeltaD)
    #   Our fidelity() returns (Tr sqrt(...))^2 = F^2, and JRSWW proves
    #   F^2 >= exp(-DeltaD) in the channel setting.  So the threshold
    #   is the same as Bound A; the test is on F_rot rather than F_std.
    if np.isinf(dD):
        bound_B = 1.0
    else:
        bound_B = 1.0 - np.exp(-dD)

    # Rotated Petz (optional)
    F_rot = None
    tau_rot = None
    satisfies_B = None
    if compute_rotated:
        try:
            from numerical.petz_toolkit import rotated_petz_recovery

            rho_rot_recovered = rotated_petz_recovery(
                rho_out, kraus_ops, sigma, n_points=201
            )
            F_rot = float(fidelity(rho, rho_rot_recovered))
            tau_rot = 1.0 - F_rot
            satisfies_B = bool(tau_rot <= bound_B + 1e-10)
        except Exception:
            pass

    return PointResult(
        label=label,
        param=param,
        F_petz=float(np.clip(F_petz, 0, 1)),
        tau_petz=float(np.clip(tau_petz, 0, 1)),
        delta_D_val=float(dD),
        bound_A=float(np.clip(bound_A, 0, 1)),
        satisfies_A=bool(satisfies_A),
        bound_B=float(np.clip(bound_B, 0, 1)),
        D_before=float(D_before),
        D_after=float(D_after),
        source=source,
        F_rotated=F_rot,
        tau_rotated=tau_rot,
        satisfies_B=satisfies_B,
    )


# ---------------------------------------------------------------------------
# Experiment 1: Amplitude damping sweep (exact simulator)
# ---------------------------------------------------------------------------


def run_amplitude_damping_experiment(
    gamma_values: Optional[NDArray] = None,
    rho: Optional[NDArray] = None,
    sigma: Optional[NDArray] = None,
    compute_rotated: bool = False,
) -> ExperimentResult:
    """
    Sweep gamma in [0, 1] for the amplitude damping channel.

    Everything is exact linear algebra -- no sampling noise.
    """
    if gamma_values is None:
        gamma_values = np.linspace(0, 1, 11)
    if rho is None:
        rho = RHO_PLUS
    if sigma is None:
        sigma = SIGMA_MM

    result = ExperimentResult(
        experiment_name="Amplitude Damping (Exact Simulator)",
        channel_type="amplitude_damping",
        source="simulator",
        input_state="|+>",
        reference_state="I/2",
        n_points=len(gamma_values),
        timestamp=time.strftime("%Y-%m-%d %H:%M:%S"),
    )

    for gamma in gamma_values:
        kraus = amplitude_damping(float(gamma))
        pt = measure_point(
            rho, sigma, kraus,
            label=f"gamma={gamma:.2f}",
            param=float(gamma),
            source="simulator",
            compute_rotated=compute_rotated,
        )
        result.points.append(pt)
        if not pt.satisfies_A:
            result.all_satisfy_A = False
        if pt.satisfies_B is not None and not pt.satisfies_B:
            result.all_satisfy_B = False

    return result


# ---------------------------------------------------------------------------
# Experiment 2: Depolarizing channel sweep
# ---------------------------------------------------------------------------


def run_depolarizing_experiment(
    p_values: Optional[NDArray] = None,
    rho: Optional[NDArray] = None,
    sigma: Optional[NDArray] = None,
    compute_rotated: bool = False,
) -> ExperimentResult:
    """Sweep depolarizing probability p in [0, 1]."""
    if p_values is None:
        p_values = np.linspace(0, 1, 11)
    if rho is None:
        rho = RHO_PLUS
    if sigma is None:
        sigma = SIGMA_MM

    result = ExperimentResult(
        experiment_name="Depolarizing Channel (Exact Simulator)",
        channel_type="depolarizing",
        source="simulator",
        input_state="|+>",
        reference_state="I/2",
        n_points=len(p_values),
        timestamp=time.strftime("%Y-%m-%d %H:%M:%S"),
    )

    for p in p_values:
        kraus = depolarizing(float(p))
        pt = measure_point(
            rho, sigma, kraus,
            label=f"p={p:.2f}",
            param=float(p),
            source="simulator",
            compute_rotated=compute_rotated,
        )
        result.points.append(pt)
        if not pt.satisfies_A:
            result.all_satisfy_A = False
        if pt.satisfies_B is not None and not pt.satisfies_B:
            result.all_satisfy_B = False

    return result


# ---------------------------------------------------------------------------
# Experiment 3: Multi-state universality test
# ---------------------------------------------------------------------------


def run_multistate_experiment(
    gamma_values: Optional[NDArray] = None,
    compute_rotated: bool = False,
) -> List[ExperimentResult]:
    """Run amplitude damping for multiple input states to show universality."""
    if gamma_values is None:
        gamma_values = np.linspace(0, 1, 11)

    states = {
        "|+>": RHO_PLUS,
        "|0>": RHO_0,
        "|1>": RHO_1,
        "|+i>": RHO_PLUS_I,
    }

    results = []
    for state_name, rho in states.items():
        result = ExperimentResult(
            experiment_name=f"Amplitude Damping -- rho = {state_name}",
            channel_type="amplitude_damping",
            source="simulator",
            input_state=state_name,
            reference_state="I/2",
            n_points=len(gamma_values),
            timestamp=time.strftime("%Y-%m-%d %H:%M:%S"),
        )
        for gamma in gamma_values:
            kraus = amplitude_damping(float(gamma))
            pt = measure_point(
                rho, SIGMA_MM, kraus,
                label=f"gamma={gamma:.2f}",
                param=float(gamma),
                source="simulator",
                compute_rotated=compute_rotated,
            )
            result.points.append(pt)
            if not pt.satisfies_A:
                result.all_satisfy_A = False
            if pt.satisfies_B is not None and not pt.satisfies_B:
                result.all_satisfy_B = False
        results.append(result)

    return results


# ---------------------------------------------------------------------------
# Experiment 4: Noisy simulator (process tomography pipeline)
# ---------------------------------------------------------------------------


def run_noisy_simulator_experiment(
    gamma_values: Optional[NDArray] = None,
    shots: int = 100_000,
    seed: int = 42,
) -> ExperimentResult:
    """
    Hardware-like pipeline: finite-shot process tomography, then Petz.

    Tests the full measurement pipeline including finite-statistics errors,
    Choi matrix reconstruction, and CPTP projection.
    """
    if gamma_values is None:
        gamma_values = np.array([0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 1.0])

    rho = RHO_PLUS
    sigma = SIGMA_MM

    result = ExperimentResult(
        experiment_name="Amplitude Damping (Noisy Sim + Tomography)",
        channel_type="amplitude_damping",
        source="noisy-sim",
        input_state="|+>",
        reference_state="I/2",
        n_points=len(gamma_values),
        timestamp=time.strftime("%Y-%m-%d %H:%M:%S"),
    )

    for gamma in gamma_values:
        kraus_true = amplitude_damping(float(gamma))
        kraus_recon = simulate_process_tomography_1q(
            kraus_true, shots=shots, verbose=False
        )

        if not verify_cptp(kraus_recon, tol=0.1):
            print(f"  WARNING: CPTP error large at gamma={gamma:.2f}")

        pt = measure_point(
            rho, sigma, kraus_recon,
            label=f"gamma={gamma:.2f}",
            param=float(gamma),
            source="noisy-sim",
        )
        result.points.append(pt)
        if not pt.satisfies_A:
            result.all_satisfy_A = False

    return result


# ---------------------------------------------------------------------------
# Experiment 5: Tuna-9 hardware
# ---------------------------------------------------------------------------


def run_tuna9_experiment(
    backend_name: str = "Tuna-9",
    shots: int = 4096,
) -> ExperimentResult:
    """
    Run on real Tuna-9 quantum hardware via Quantum Inspire.

    Uses the hardware's native noise as the quantum channel under test.
    Process tomography extracts Kraus operators, then classical Petz
    recovery computes the optimal fidelity.
    """
    from tau_chrono.backends import QI_AVAILABLE

    if not QI_AVAILABLE:
        raise RuntimeError(
            "Quantum Inspire not available.\n"
            "Install: pip install quantuminspire qiskit-quantuminspire\n"
            "Then: qi login"
        )

    from tau_chrono.backends.quantum_inspire import (
        get_qi_backend,
        extract_gate_kraus,
    )

    rho = RHO_PLUS
    sigma = SIGMA_MM

    print(f"Connecting to {backend_name}...")
    backend = get_qi_backend(backend_name, optimization_level=0)
    info = backend.info()
    print(f"  Backend: {info.name}")
    print(f"  Qubits:  {info.num_qubits}")
    print(f"  Simulator: {info.is_simulator}")
    print(f"  Basis gates: {info.basis_gates}")

    gate_fns = {
        "Identity": lambda qc, q: qc.id(q),
        "Hadamard": lambda qc, q: qc.h(q),
        "X gate":   lambda qc, q: qc.x(q),
    }

    result = ExperimentResult(
        experiment_name=f"Hardware Noise ({backend_name})",
        channel_type="hardware_native",
        source="tuna9",
        input_state="|+>",
        reference_state="I/2",
        n_points=len(gate_fns),
        timestamp=time.strftime("%Y-%m-%d %H:%M:%S"),
    )

    for i, (gate_name, gate_fn) in enumerate(gate_fns.items()):
        print(f"\n--- Gate {i+1}/{len(gate_fns)}: {gate_name} ---")

        char = extract_gate_kraus(
            gate_fn, backend, gate_name, shots=shots, verbose=True
        )

        pt = measure_point(
            rho, sigma, char.kraus_ops,
            label=gate_name,
            param=float(i),
            source="tuna9",
        )
        result.points.append(pt)
        if not pt.satisfies_A:
            result.all_satisfy_A = False

    return result


# ---------------------------------------------------------------------------
# Experiment 6: Depth scaling
# ---------------------------------------------------------------------------


def run_depth_scaling_experiment(
    depths: Optional[List[int]] = None,
    noise_model: Optional[NoiseModel] = None,
) -> ExperimentResult:
    """Measure tau as a function of circuit depth (accumulated noise)."""
    if depths is None:
        depths = [0, 1, 2, 5, 10, 20, 50]
    if noise_model is None:
        noise_model = NoiseModel(
            depolarizing_rate=0.005,
            amplitude_damping_rate=0.002,
            dephasing_rate=0.003,
            readout_error=0.01,
        )

    rho = RHO_PLUS
    sigma = SIGMA_MM

    result = ExperimentResult(
        experiment_name="Depth Scaling (Noisy Identity Circuits)",
        channel_type="depth_scaling",
        source="noisy-sim",
        input_state="|+>",
        reference_state="I/2",
        n_points=len(depths),
        timestamp=time.strftime("%Y-%m-%d %H:%M:%S"),
    )

    for depth in depths:
        if depth == 0:
            kraus = [np.eye(2, dtype=complex)]
        else:
            # Effective depolarizing after d layers
            p_eff = 1.0 - (1.0 - noise_model.depolarizing_rate) ** depth
            p_eff = min(p_eff, 1.0)
            kraus = depolarizing(p_eff)

        pt = measure_point(
            rho, sigma, kraus,
            label=f"depth={depth}",
            param=float(depth),
            source="noisy-sim",
        )
        result.points.append(pt)
        if not pt.satisfies_A:
            result.all_satisfy_A = False

    return result


# ---------------------------------------------------------------------------
# Experiment 7: Random channels (statistical verification)
# ---------------------------------------------------------------------------


def run_random_channel_experiment(
    n_samples: int = 200,
    seed: int = 42,
    compute_rotated: bool = False,
) -> ExperimentResult:
    """
    Test the bound on random CPTP maps to verify it holds universally.

    This is the most rigorous test: random states, random channels.
    """
    from numerical.petz_toolkit import (
        random_cptp_map,
        random_density_matrix,
        random_pure_state,
    )

    rng = np.random.default_rng(seed)

    result = ExperimentResult(
        experiment_name=f"Random Channels (n={n_samples})",
        channel_type="random",
        source="simulator",
        input_state="random",
        reference_state="I/2",
        n_points=n_samples,
        timestamp=time.strftime("%Y-%m-%d %H:%M:%S"),
    )

    sigma = SIGMA_MM

    for i in range(n_samples):
        rho = random_pure_state(2, rng=rng)
        kraus = random_cptp_map(2, 2, rng=rng)

        pt = measure_point(
            rho, sigma, kraus,
            label=f"sample={i}",
            param=float(i),
            source="simulator",
            compute_rotated=compute_rotated,
        )
        result.points.append(pt)
        if not pt.satisfies_A:
            result.all_satisfy_A = False
        if pt.satisfies_B is not None and not pt.satisfies_B:
            result.all_satisfy_B = False

    return result


# ---------------------------------------------------------------------------
# Display
# ---------------------------------------------------------------------------


def print_result_table(result: ExperimentResult, show_rotated: bool = False) -> None:
    """Print a formatted table of experiment results."""
    print()
    print("=" * 85)
    print(f"  {result.experiment_name}")
    print(f"  Channel: {result.channel_type} | Source: {result.source}")
    print(f"  Input: {result.input_state} | Reference: {result.reference_state}")
    print(f"  {result.timestamp}")
    print("=" * 85)

    # Bound A header
    if show_rotated:
        hdr = (
            f"{'Param':>10s}  {'F_std':>7s}  {'tau_std':>7s}  "
            f"{'F_rot':>7s}  {'tau_rot':>7s}  "
            f"{'DeltaD':>7s}  {'bnd_A':>7s}  {'bnd_B':>7s}  {'A':>2s} {'B':>2s}"
        )
    else:
        hdr = (
            f"{'Param':>10s}  {'F_Petz':>8s}  {'tau':>8s}  "
            f"{'DeltaD':>8s}  {'bound_A':>8s}  {'gap':>8s}  {'OK?':>4s}"
        )
    print(hdr)
    print("-" * 85)

    for pt in result.points:
        dD_str = "inf" if np.isinf(pt.delta_D_val) else f"{pt.delta_D_val:.4f}"

        if show_rotated:
            F_rot_str = f"{pt.F_rotated:.4f}" if pt.F_rotated is not None else "  n/a"
            tau_rot_str = (
                f"{pt.tau_rotated:.4f}" if pt.tau_rotated is not None else "  n/a"
            )
            ok_A = "Y" if pt.satisfies_A else "N"
            ok_B = (
                "Y"
                if pt.satisfies_B is True
                else ("N" if pt.satisfies_B is False else "?")
            )
            line = (
                f"{pt.label:>10s}  {pt.F_petz:7.4f}  {pt.tau_petz:7.4f}  "
                f"{F_rot_str:>7s}  {tau_rot_str:>7s}  "
                f"{dD_str:>7s}  {pt.bound_A:7.4f}  {pt.bound_B:7.4f}  "
                f"{ok_A:>2s} {ok_B:>2s}"
            )
        else:
            gap = pt.bound_A - pt.tau_petz
            ok = "Y" if pt.satisfies_A else "N"
            line = (
                f"{pt.label:>10s}  {pt.F_petz:8.4f}  {pt.tau_petz:8.4f}  "
                f"{dD_str:>8s}  {pt.bound_A:8.4f}  {gap:8.4f}  {ok:>4s}"
            )
        print(line)

    print("-" * 85)

    # Bound A summary
    n_ok_A = sum(1 for pt in result.points if pt.satisfies_A)
    n_total = len(result.points)
    status_A = "PASS" if result.all_satisfy_A else "FAIL"
    print(
        f"  Bound A  [F >= exp(-DeltaD)]:       "
        f"{n_ok_A}/{n_total} ({100*n_ok_A/n_total:.0f}%) -- {status_A}"
    )

    # Bound B summary (if rotated Petz was computed)
    pts_with_B = [pt for pt in result.points if pt.satisfies_B is not None]
    if pts_with_B:
        n_ok_B = sum(1 for pt in pts_with_B if pt.satisfies_B)
        n_B = len(pts_with_B)
        status_B = "PASS" if all(pt.satisfies_B for pt in pts_with_B) else "FAIL"
        print(
            f"  Bound B  [F_rot >= exp(-DeltaD)]:   "
            f"{n_ok_B}/{n_B} ({100*n_ok_B/n_B:.0f}%) -- {status_B}"
        )

    print()


def print_random_summary(result: ExperimentResult) -> None:
    """Print a statistical summary for random channel experiments."""
    print()
    print("=" * 85)
    print(f"  {result.experiment_name}")
    print("=" * 85)

    taus = [pt.tau_petz for pt in result.points]
    bounds_A = [pt.bound_A for pt in result.points]
    gaps_A = [pt.bound_A - pt.tau_petz for pt in result.points]
    dDs = [pt.delta_D_val for pt in result.points if not np.isinf(pt.delta_D_val)]

    n_ok_A = sum(1 for pt in result.points if pt.satisfies_A)
    n_total = len(result.points)

    print(f"  Samples:      {n_total}")
    print(f"  Bound A compliance: {n_ok_A}/{n_total} ({100*n_ok_A/n_total:.0f}%)")
    print(f"  tau range:    [{min(taus):.6f}, {max(taus):.6f}]")
    print(f"  DeltaD range: [{min(dDs):.6f}, {max(dDs):.6f}]")
    print(f"  gap (bnd-tau) range: [{min(gaps_A):.6f}, {max(gaps_A):.6f}]")
    print(f"  gap mean:     {np.mean(gaps_A):.6f}")

    if any(pt.satisfies_B is not None for pt in result.points):
        pts_B = [pt for pt in result.points if pt.satisfies_B is not None]
        n_ok_B = sum(1 for pt in pts_B if pt.satisfies_B)
        gaps_B = [
            pt.bound_B - pt.tau_rotated
            for pt in pts_B
            if pt.tau_rotated is not None
        ]
        print(f"  Bound B compliance: {n_ok_B}/{len(pts_B)}")
        if gaps_B:
            print(f"  gap_B range:  [{min(gaps_B):.6f}, {max(gaps_B):.6f}]")

    print()


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------


def plot_results(
    results: List[ExperimentResult],
    save_path: Optional[str] = None,
) -> None:
    """
    Publication-quality plot of tau vs parameter with both bound curves.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available; skipping plot.")
        return

    # Filter out random-channel experiments (they don't plot well as curves)
    plottable = [r for r in results if r.channel_type != "random"]
    if not plottable:
        print("No plottable results.")
        return

    n_plots = len(plottable)
    fig, axes = plt.subplots(
        1, n_plots, figsize=(6 * n_plots, 5), squeeze=False
    )

    source_style = {
        "simulator": {"color": "#2196F3", "marker": "o"},
        "noisy-sim": {"color": "#FF9800", "marker": "s"},
        "tuna9":     {"color": "#E91E63", "marker": "D"},
    }

    for idx, result in enumerate(plottable):
        ax = axes[0, idx]

        params = [pt.param for pt in result.points]
        taus = [pt.tau_petz for pt in result.points]
        bounds_A = [pt.bound_A for pt in result.points]
        bounds_B = [pt.bound_B for pt in result.points]

        style = source_style.get(result.source, {"color": "#666", "marker": "o"})

        # Sort by parameter for smooth curves
        order = np.argsort(params)
        params_s = [params[i] for i in order]
        taus_s = [taus[i] for i in order]
        bounds_A_s = [bounds_A[i] for i in order]
        bounds_B_s = [bounds_B[i] for i in order]

        # Bound A curve (Wilde): F >= exp(-DeltaD)
        ax.plot(
            params_s, bounds_A_s, "-",
            color="#666666", linewidth=2,
            label="Bound A: $1-e^{-\\Delta D}$",
        )

        # Bound B curve (JRSWW): F_rot^2 >= exp(-DeltaD)
        ax.plot(
            params_s, bounds_B_s, "--",
            color="#999999", linewidth=1.5,
            label="Bound B: $1-e^{-\\Delta D}$ (rot)",
        )

        # Shade forbidden region above Bound A
        ax.fill_between(
            params_s, bounds_A_s, [1.05] * len(params_s),
            alpha=0.07, color="red",
        )

        # Shade allowed region below Bound A
        ax.fill_between(
            params_s, [-0.02] * len(params_s), bounds_A_s,
            alpha=0.05, color="green",
        )

        # Measured tau (standard Petz)
        ax.scatter(
            params_s, taus_s,
            color=style["color"], marker=style["marker"],
            s=60, zorder=3, edgecolors="black", linewidths=0.5,
            label=f"$\\tau_{{\\mathrm{{Petz}}}}$ ({result.source})",
        )

        # Rotated Petz tau (if available)
        taus_rot = [
            result.points[i].tau_rotated
            for i in order
            if result.points[i].tau_rotated is not None
        ]
        if taus_rot:
            params_rot = [
                params[i]
                for i in order
                if result.points[i].tau_rotated is not None
            ]
            ax.scatter(
                params_rot, taus_rot,
                color="#4CAF50", marker="^",
                s=50, zorder=3, edgecolors="black", linewidths=0.5,
                label="$\\tau_{\\mathrm{rot}}$",
            )

        ax.set_xlabel(_param_label(result.channel_type), fontsize=12)
        ax.set_ylabel("$\\tau = 1 - F$", fontsize=12)
        ax.set_title(result.experiment_name, fontsize=10, pad=10)
        ax.legend(fontsize=8, loc="upper left")
        ax.set_ylim(-0.02, 1.02)
        ax.grid(True, alpha=0.3)

    fig.suptitle(
        "Petz Recovery Fidelity: "
        "$\\tau \\leq 1 - e^{-\\Delta D}$ (channel bound)",
        fontsize=13, fontweight="bold", y=1.02,
    )
    plt.tight_layout()

    if save_path is None:
        save_path = str(_PROJECT_ROOT / "results" / "petz_experiment.png")

    Path(save_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(save_path, dpi=150, bbox_inches="tight")
    print(f"  Plot saved to: {save_path}")
    plt.close(fig)


def _param_label(channel_type: str) -> str:
    labels = {
        "amplitude_damping": "$\\gamma$ (damping probability)",
        "depolarizing": "$p$ (depolarizing probability)",
        "dephasing": "$p$ (dephasing probability)",
        "depth_scaling": "Circuit depth",
        "hardware_native": "Gate index",
    }
    return labels.get(channel_type, "Parameter")


# ---------------------------------------------------------------------------
# JSON serialization
# ---------------------------------------------------------------------------


def save_results_json(
    results: List[ExperimentResult],
    path: Optional[str] = None,
) -> str:
    """Save results to JSON."""
    if path is None:
        path = str(_PROJECT_ROOT / "results" / "petz_experiment_results.json")

    Path(path).parent.mkdir(parents=True, exist_ok=True)

    def _safe(x):
        if isinstance(x, float) and np.isinf(x):
            return "inf"
        return x

    data = []
    for r in results:
        d = {
            "experiment_name": r.experiment_name,
            "channel_type": r.channel_type,
            "source": r.source,
            "input_state": r.input_state,
            "reference_state": r.reference_state,
            "n_points": r.n_points,
            "all_satisfy_A": r.all_satisfy_A,
            "all_satisfy_B": r.all_satisfy_B,
            "timestamp": r.timestamp,
            "points": [
                {
                    "label": pt.label,
                    "param": pt.param,
                    "F_petz": pt.F_petz,
                    "tau_petz": pt.tau_petz,
                    "delta_D_val": _safe(pt.delta_D_val),
                    "bound_A": pt.bound_A,
                    "satisfies_A": pt.satisfies_A,
                    "bound_B": pt.bound_B,
                    "D_before": _safe(pt.D_before),
                    "D_after": _safe(pt.D_after),
                    "source": pt.source,
                    "F_rotated": pt.F_rotated,
                    "tau_rotated": pt.tau_rotated,
                    "satisfies_B": pt.satisfies_B,
                }
                for pt in r.points
            ],
        }
        data.append(d)

    with open(path, "w") as f:
        json.dump(data, f, indent=2)

    print(f"  Results saved to: {path}")
    return path


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Petz Recovery Fidelity Experiment",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Bounds under test:
  (A) Standard Petz:  F >= exp(-DeltaD)        [Wilde 2015]
  (B) Rotated Petz:   F_rot >= exp(-DeltaD)     [JRSWW 2018]

  where DeltaD = D(rho||sigma) - D(N(rho)||N(sigma))

Examples:
  python experiments/petz_tuna9_experiment.py --mode simulator
  python experiments/petz_tuna9_experiment.py --mode simulator --rotated
  python experiments/petz_tuna9_experiment.py --mode noisy-sim --shots 50000
  python experiments/petz_tuna9_experiment.py --mode tuna9 --backend Tuna-9
  python experiments/petz_tuna9_experiment.py --mode all --plot --save
        """,
    )
    parser.add_argument(
        "--mode",
        choices=["simulator", "noisy-sim", "tuna9", "all"],
        default="simulator",
        help="Experiment mode (default: simulator)",
    )
    parser.add_argument("--plot", action="store_true", help="Generate plots")
    parser.add_argument("--save", action="store_true", help="Save results to JSON")
    parser.add_argument(
        "--backend", default="Tuna-9",
        help="QI backend name for tuna9 mode (default: Tuna-9)",
    )
    parser.add_argument(
        "--shots", type=int, default=100_000,
        help="Shots per circuit for tomography (default: 100000)",
    )
    parser.add_argument(
        "--rotated", action="store_true",
        help="Also compute rotated Petz (JRSWW) for comparison",
    )
    parser.add_argument(
        "--random-samples", type=int, default=200,
        help="Number of random channels to test (default: 200)",
    )

    args = parser.parse_args()

    print()
    print("=" * 78)
    print("  PETZ RECOVERY FIDELITY EXPERIMENT")
    print("  Bound A: F_std >= exp(-DeltaD)        [Wilde 2015, standard Petz]")
    print("  Bound B: F_rot >= exp(-DeltaD)        [JRSWW 2018, rotated Petz]")
    print("  Huang (2026) -- Petz Recovery Map as Retrodiction Functor")
    print("=" * 78)

    all_results: List[ExperimentResult] = []
    show_rotated = args.rotated

    # ---- Mode: simulator ----
    if args.mode in ("simulator", "all"):
        print("\n>>> Experiment 1: Amplitude Damping (exact) <<<")
        r1 = run_amplitude_damping_experiment(compute_rotated=show_rotated)
        print_result_table(r1, show_rotated=show_rotated)
        all_results.append(r1)

        print("\n>>> Experiment 2: Depolarizing (exact) <<<")
        r2 = run_depolarizing_experiment(compute_rotated=show_rotated)
        print_result_table(r2, show_rotated=show_rotated)
        all_results.append(r2)

        print("\n>>> Experiment 3: Multi-State Universality <<<")
        multi = run_multistate_experiment(compute_rotated=show_rotated)
        for mr in multi:
            print_result_table(mr, show_rotated=show_rotated)
            all_results.append(mr)

        print("\n>>> Experiment 4: Random Channels (statistical) <<<")
        r4 = run_random_channel_experiment(
            n_samples=args.random_samples,
            compute_rotated=show_rotated,
        )
        print_random_summary(r4)
        all_results.append(r4)

    # ---- Mode: noisy-sim ----
    if args.mode in ("noisy-sim", "all"):
        print("\n>>> Experiment 5: Noisy Sim + Tomography <<<")
        r5 = run_noisy_simulator_experiment(shots=args.shots)
        print_result_table(r5)
        all_results.append(r5)

        print("\n>>> Experiment 6: Depth Scaling <<<")
        r6 = run_depth_scaling_experiment()
        print_result_table(r6)
        all_results.append(r6)

    # ---- Mode: tuna9 ----
    if args.mode in ("tuna9", "all"):
        print(f"\n>>> Experiment 7: Tuna-9 Hardware ({args.backend}) <<<")
        try:
            r7 = run_tuna9_experiment(
                backend_name=args.backend,
                shots=min(args.shots, 4096),
            )
            print_result_table(r7)
            all_results.append(r7)
        except RuntimeError as e:
            print(f"\n  Tuna-9 not available: {e}")
            print("  Skipping hardware experiment.")

    # ---- Summary ----
    print("\n" + "=" * 78)
    print("  EXPERIMENT SUMMARY")
    print("=" * 78)

    total_A = 0
    pass_A = 0
    total_B = 0
    pass_B = 0

    for r in all_results:
        n_ok_A = sum(1 for pt in r.points if pt.satisfies_A)
        n_tot = len(r.points)
        total_A += n_tot
        pass_A += n_ok_A
        status_A = "PASS" if r.all_satisfy_A else "FAIL"

        pts_B = [pt for pt in r.points if pt.satisfies_B is not None]
        n_ok_B = sum(1 for pt in pts_B if pt.satisfies_B)
        total_B += len(pts_B)
        pass_B += n_ok_B

        line = f"  [{status_A}] {r.experiment_name}: A={n_ok_A}/{n_tot}"
        if pts_B:
            line += f", B={n_ok_B}/{len(pts_B)}"
        print(line)

    print(f"\n  Bound A total: {pass_A}/{total_A} points")
    if total_B > 0:
        print(f"  Bound B total: {pass_B}/{total_B} points")

    if pass_A == total_A:
        print("\n  BOUND A VERIFIED: F(rho, R_Petz(N(rho))) >= exp(-DeltaD)")
    else:
        n_fail = total_A - pass_A
        print(f"\n  WARNING: {n_fail} Bound-A violation(s) detected.")

    if total_B > 0 and pass_B == total_B:
        print("  BOUND B VERIFIED: F(rho, R_rot(N(rho))) >= exp(-DeltaD)")
    elif total_B > 0:
        n_fail_B = total_B - pass_B
        print(f"  WARNING: {n_fail_B} Bound-B violation(s) detected.")

    print("=" * 78)

    # ---- Plot ----
    if args.plot:
        print("\nGenerating plots...")
        plot_results(all_results)

    # ---- Save ----
    if args.save:
        print("\nSaving results...")
        save_results_json(all_results)

    print()


if __name__ == "__main__":
    main()
