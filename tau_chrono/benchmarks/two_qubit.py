"""
2-qubit benchmark for tau-chrono Bayesian Noise Tracker.

Four circuits comparing Bayesian sigma-tracking vs naive multiplicative:
    A: Bell state + idle noise (4 rounds)
    B: SWAP via 3 CNOTs with noise between each
    C: Random 2-qubit circuit (4 layers)
    D: VQE variational ansatz (4 layers)

Run as CLI::

    python -m tau_chrono.benchmarks.two_qubit
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray
from typing import List
from dataclasses import dataclass

from tau_chrono.channels import (
    amplitude_damping,
    depolarizing,
    verify_cptp,
    two_qubit_depolarizing,
    correlated_dephasing,
    local_noise,
    crosstalk_dephasing,
    amplitude_damping_2q,
    unitary_channel,
    CNOT,
    CNOT_REV,
    I4,
)
from tau_chrono.petz import apply_channel, tau_parameter
from tau_chrono.bayesian import compose_kraus_compressed
from tau_chrono.utils import commutator_norm

KrausList = List[NDArray[np.complexfloating]]
DensityMatrix = NDArray[np.complexfloating]

ket0 = np.array([1, 0], dtype=complex)
ket1 = np.array([0, 1], dtype=complex)


def Ry(theta: float) -> NDArray:
    """Single-qubit Y-rotation gate."""
    return np.array([
        [np.cos(theta / 2), -np.sin(theta / 2)],
        [np.sin(theta / 2),  np.cos(theta / 2)],
    ], dtype=complex)


def Rz(theta: float) -> NDArray:
    """Single-qubit Z-rotation gate."""
    return np.array([
        [np.exp(-1j * theta / 2), 0],
        [0, np.exp(1j * theta / 2)],
    ], dtype=complex)


# ---------------------------------------------------------------------------
# Result dataclasses
# ---------------------------------------------------------------------------


@dataclass
class GateResult2Q:
    """Result for a single gate in a 2-qubit circuit."""
    gate_index: int
    channel_name: str
    tau_naive: float
    tau_eff: float
    improvement_pct: float
    classification: str  # UNITARY, NORMAL, NEAR_SATURATED, SATURATED


@dataclass
class CircuitResult:
    """Full result of a 2-qubit circuit benchmark."""
    circuit_name: str
    n_gates: int
    n_noise_gates: int
    noise_description: str
    gate_results: List[GateResult2Q]
    tau_bayesian_total: float
    tau_multiplicative_total: float
    improvement_pct: float
    sum_tau_naive: float
    sum_tau_eff: float
    improvement_pergate_pct: float
    best_per_gate_improvement: float
    composition_lhs: float
    composition_rhs: float
    composition_holds: bool


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def make_2q_sigma(eigenvalues: list) -> DensityMatrix:
    """Create a diagonal 2-qubit reference state from eigenvalue list (normalised)."""
    ev = np.array(eigenvalues, dtype=float)
    ev = ev / ev.sum()
    return np.diag(ev).astype(complex)


def run_circuit_benchmark(
    circuit_name: str,
    channels: List[KrausList],
    channel_names: List[str],
    rho_0: DensityMatrix,
    sigma_0: DensityMatrix,
    noise_description: str,
) -> CircuitResult:
    """Run a full circuit benchmark comparing naive vs Bayesian tau."""
    n_gates = len(channels)
    gate_results: List[GateResult2Q] = []
    rho_curr = rho_0.copy()
    sigma_curr = sigma_0.copy()
    n_noise = 0

    for i in range(n_gates):
        kraus = channels[i]
        is_unitary = (len(kraus) == 1 and
                      np.linalg.norm(kraus[0] @ kraus[0].conj().T - I4) < 1e-10)

        if is_unitary:
            tau_n = 0.0
            tau_e = 0.0
            imp = 0.0
            classification = "UNITARY"
        else:
            tau_n = tau_parameter(rho_0, kraus, sigma_0)
            tau_e = tau_parameter(rho_curr, kraus, sigma_curr)
            imp = (tau_n - tau_e) / tau_n * 100.0 if tau_n > 1e-12 else 0.0
            n_noise += 1
            if tau_e < 1e-10:
                classification = "SATURATED"
            elif commutator_norm(
                apply_channel(rho_curr, kraus),
                apply_channel(sigma_curr, kraus),
            ) < 0.01:
                classification = "NEAR_SATURATED"
            else:
                classification = "NORMAL"

        gate_results.append(GateResult2Q(
            gate_index=i,
            channel_name=channel_names[i],
            tau_naive=tau_n,
            tau_eff=tau_e,
            improvement_pct=imp,
            classification=classification,
        ))
        rho_curr = apply_channel(rho_curr, kraus)
        sigma_curr = apply_channel(sigma_curr, kraus)

    # Total composed channel
    composed = compose_kraus_compressed(channels)
    tau_bayesian_total = tau_parameter(rho_0, composed, sigma_0)

    prod_fid = 1.0
    for gr in gate_results:
        prod_fid *= (1.0 - gr.tau_naive)
    tau_mult = 1.0 - prod_fid

    improvement = (tau_mult - tau_bayesian_total) / tau_mult * 100.0 if tau_mult > 1e-15 else 0.0

    noise_gates = [gr for gr in gate_results if gr.classification != "UNITARY"]
    sum_naive = sum(gr.tau_naive for gr in noise_gates)
    sum_eff = sum(gr.tau_eff for gr in noise_gates)
    imp_pergate = (sum_naive - sum_eff) / sum_naive * 100.0 if sum_naive > 1e-15 else 0.0
    best_gate = max((gr.improvement_pct for gr in noise_gates), default=0.0)

    sum_sqrt_eff = sum(np.sqrt(max(gr.tau_eff, 0.0)) for gr in gate_results)
    sqrt_total = np.sqrt(max(tau_bayesian_total, 0.0))
    comp_holds = sqrt_total <= sum_sqrt_eff + 1e-8

    return CircuitResult(
        circuit_name=circuit_name,
        n_gates=n_gates,
        n_noise_gates=n_noise,
        noise_description=noise_description,
        gate_results=gate_results,
        tau_bayesian_total=tau_bayesian_total,
        tau_multiplicative_total=tau_mult,
        improvement_pct=improvement,
        sum_tau_naive=sum_naive,
        sum_tau_eff=sum_eff,
        improvement_pergate_pct=imp_pergate,
        best_per_gate_improvement=best_gate,
        composition_lhs=sqrt_total,
        composition_rhs=sum_sqrt_eff,
        composition_holds=comp_holds,
    )


# ---------------------------------------------------------------------------
# Circuit definitions
# ---------------------------------------------------------------------------


def circuit_A_bell_state(rng: np.random.Generator,
                         sigma_0: DensityMatrix) -> CircuitResult:
    """Circuit A: Bell state + idle noise (4 rounds)."""
    psi_plus0 = np.kron((ket0 + ket1) / np.sqrt(2), ket0)
    rho_0 = np.outer(psi_plus0, psi_plus0.conj())

    channels: List[KrausList] = [unitary_channel(CNOT)]
    names = ["CNOT"]

    for k in range(4):
        gamma = 0.03 + 0.015 * k
        dp = 0.03 + 0.01 * k
        channels.append(local_noise(amplitude_damping(gamma), amplitude_damping(gamma)))
        names.append(f"LocalAD({gamma:.3f})_{k}")
        channels.append(local_noise(depolarizing(dp), depolarizing(dp)))
        names.append(f"LocalDep({dp:.2f})_{k}")
        channels.append(correlated_dephasing(0.03))
        names.append(f"CorrDeph(0.03)_{k}")
        theta = 0.15 * (k + 1)
        U_rz = np.kron(Rz(theta), Rz(-theta))
        channels.append(unitary_channel(U_rz))
        names.append(f"Rz({theta:.2f})_{k}")

    return run_circuit_benchmark(
        "A: Bell + Idle Noise", channels, names, rho_0, sigma_0,
        "4x[AD(0.03-0.075) + Dep(0.03-0.06) + CorrDeph(0.03) + Rz]",
    )


def circuit_B_swap(rng: np.random.Generator,
                   sigma_0: DensityMatrix) -> CircuitResult:
    """Circuit B: SWAP via 3 CNOTs with noise between each."""
    psi_plus = (ket0 + ket1) / np.sqrt(2)
    psi_init = np.kron(psi_plus, ket0)
    rho_0 = np.outer(psi_init, psi_init.conj())

    channels: List[KrausList] = []
    names: List[str] = []

    for step, (cnot_gate, cnot_name) in enumerate([
        (CNOT, "CNOT(0,1)"), (CNOT_REV, "CNOT(1,0)"), (CNOT, "CNOT(0,1)")
    ]):
        channels.append(unitary_channel(cnot_gate))
        names.append(cnot_name)
        channels.append(local_noise(depolarizing(0.05), depolarizing(0.05)))
        names.append(f"LocalDep(0.05)_{step}")
        channels.append(local_noise(amplitude_damping(0.04), amplitude_damping(0.04)))
        names.append(f"LocalAD(0.04)_{step}")
        channels.append(correlated_dephasing(0.04))
        names.append(f"CorrDeph(0.04)_{step}")
        channels.append(crosstalk_dephasing(0.03, 0.2))
        names.append(f"Xtalk(0.03)_{step}")

    return run_circuit_benchmark(
        "B: SWAP (3 CNOTs)", channels, names, rho_0, sigma_0,
        "3x[CNOT + Dep(0.05) + AD(0.04) + CorrDeph(0.04) + Xtalk(0.03)]",
    )


def circuit_C_random_2q(rng: np.random.Generator,
                        sigma_0: DensityMatrix) -> CircuitResult:
    """Circuit C: Random 2-qubit circuit with 4 layers."""
    psi = rng.standard_normal(4) + 1j * rng.standard_normal(4)
    psi /= np.linalg.norm(psi)
    rho_0 = np.outer(psi, psi.conj())

    channels: List[KrausList] = []
    names: List[str] = []
    noise_params = [
        (0.06, 0.04, 0.03),
        (0.07, 0.05, 0.04),
        (0.08, 0.06, 0.03),
        (0.09, 0.07, 0.05),
    ]

    for layer_idx in range(4):
        theta1 = rng.uniform(0, 2 * np.pi)
        theta2 = rng.uniform(0, 2 * np.pi)
        phi1 = rng.uniform(0, 2 * np.pi)
        phi2 = rng.uniform(0, 2 * np.pi)
        U_rot = np.kron(Rz(phi1) @ Ry(theta1), Rz(phi2) @ Ry(theta2))
        channels.append(unitary_channel(U_rot))
        names.append(f"Rot_L{layer_idx}")
        dp, ag, dph = noise_params[layer_idx]
        channels.append(local_noise(depolarizing(dp), amplitude_damping(ag)))
        names.append(f"Noise_L{layer_idx}")
        if layer_idx % 2 == 0:
            channels.append(unitary_channel(CNOT))
            names.append(f"CNOT_L{layer_idx}")
        else:
            channels.append(unitary_channel(CNOT_REV))
            names.append(f"CNOTrev_L{layer_idx}")
        channels.append(correlated_dephasing(dph))
        names.append(f"CorrDeph_L{layer_idx}")

    return run_circuit_benchmark(
        "C: Random 2Q (4 layers)", channels, names, rho_0, sigma_0,
        "4x[Rot + LocalNoise + CNOT + CorrDeph], escalating noise",
    )


def circuit_D_variational(rng: np.random.Generator,
                          sigma_0: DensityMatrix) -> CircuitResult:
    """Circuit D: VQE variational ansatz -- 4 layers."""
    psi_00 = np.kron(ket0, ket0)
    rho_0 = np.outer(psi_00, psi_00.conj())
    thetas = rng.uniform(0, np.pi, size=8)

    channels: List[KrausList] = []
    names: List[str] = []
    noise_schedule = [
        (0.03, 0.02, 0.04, 0.03),
        (0.04, 0.03, 0.05, 0.04),
        (0.04, 0.03, 0.05, 0.04),
        (0.05, 0.04, 0.06, 0.05),
    ]

    for layer in range(4):
        U = np.kron(Ry(thetas[2 * layer]), Ry(thetas[2 * layer + 1]))
        channels.append(unitary_channel(U))
        names.append(f"Ry_L{layer}")
        dp, ag, dep2q, dph = noise_schedule[layer]
        channels.append(local_noise(depolarizing(dp), depolarizing(dp)))
        names.append(f"LocDep({dp})_L{layer}")
        channels.append(amplitude_damping_2q(ag))
        names.append(f"LocAD({ag})_L{layer}")
        if layer % 2 == 0:
            channels.append(unitary_channel(CNOT))
            names.append(f"CNOT_L{layer}")
        else:
            channels.append(unitary_channel(CNOT_REV))
            names.append(f"CNOTrev_L{layer}")
        channels.append(two_qubit_depolarizing(dep2q))
        names.append(f"2QDep({dep2q})_L{layer}")
        channels.append(correlated_dephasing(dph))
        names.append(f"CorrDeph({dph})_L{layer}")

    return run_circuit_benchmark(
        "D: VQE Ansatz (4 layers)", channels, names, rho_0, sigma_0,
        "4x[Ry + LocDep + LocAD + CNOT + 2QDep + CorrDeph]",
    )


# ---------------------------------------------------------------------------
# Pretty printing
# ---------------------------------------------------------------------------


def print_circuit_result(result: CircuitResult, verbose: bool = True) -> None:
    """Print a formatted circuit result."""
    print(f"\n{'=' * 82}")
    print(f"  Circuit: {result.circuit_name}")
    print(f"  Layers: {result.n_gates}  |  Noise layers: {result.n_noise_gates}"
          f"  |  {result.noise_description}")
    print(f"{'=' * 82}")

    if verbose:
        print(f"\n  {'Gate':<26} {'tau_naive':>10} {'tau_eff':>10} "
              f"{'Improve%':>10} {'Class':<16}")
        print(f"  {'-'*26} {'-'*10} {'-'*10} {'-'*10} {'-'*16}")
        for gr in result.gate_results:
            if gr.classification == "UNITARY":
                print(f"  {gr.channel_name:<26} {'--':>10} "
                      f"{'--':>10} {'--':>10} {'UNITARY':<16}")
            else:
                print(f"  {gr.channel_name:<26} {gr.tau_naive:>10.6f} "
                      f"{gr.tau_eff:>10.6f} {gr.improvement_pct:>+10.2f}% "
                      f"{gr.classification:<16}")

    print(f"\n  --- Total circuit ---")
    print(f"    tau_multiplicative:    {result.tau_multiplicative_total:.6f}")
    print(f"    tau_bayesian:          {result.tau_bayesian_total:.6f}")
    print(f"    IMPROVEMENT:           {result.improvement_pct:+.2f}%")
    print(f"    Composition holds:     {result.composition_holds}")


def print_summary_table(results: List[CircuitResult]) -> None:
    """Print a comparison table of all circuits."""
    print(f"\n{'=' * 96}")
    print(f"  TAU-CHRONO 2-QUBIT BENCHMARK -- SUMMARY")
    print(f"{'=' * 96}")
    print(f"\n  {'Circuit':<28} {'#Noise':>6} {'tau_mult':>10} {'tau_bayes':>10} "
          f"{'Overall%':>10} {'BestGate%':>10} {'Comp':>5}")
    print(f"  {'-'*28} {'-'*6} {'-'*10} {'-'*10} {'-'*10} {'-'*10} {'-'*5}")

    for r in results:
        comp = "OK" if r.composition_holds else "FAIL"
        print(f"  {r.circuit_name:<28} {r.n_noise_gates:>6} "
              f"{r.tau_multiplicative_total:>10.6f} {r.tau_bayesian_total:>10.6f} "
              f"{r.improvement_pct:>+10.2f}% "
              f"{r.best_per_gate_improvement:>+10.2f}% {comp:>5}")

    avg_imp = np.mean([r.improvement_pct for r in results])
    max_imp = max(r.improvement_pct for r in results)
    all_comp = all(r.composition_holds for r in results)
    print(f"\n  Avg improvement: {avg_imp:+.2f}%  |  Max: {max_imp:+.2f}%"
          f"  |  All compositions hold: {all_comp}")
    print(f"{'=' * 96}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def run() -> List[CircuitResult]:
    """Run all 4 circuits and return results."""
    rng = np.random.default_rng(42)
    sigma_0 = make_2q_sigma([0.4, 0.3, 0.2, 0.1])
    results = [
        circuit_A_bell_state(rng, sigma_0),
        circuit_B_swap(rng, sigma_0),
        circuit_C_random_2q(rng, sigma_0),
        circuit_D_variational(rng, sigma_0),
    ]
    return results


def main() -> None:
    """Pretty-print all 2-qubit benchmarks."""
    np.random.seed(42)
    rng = np.random.default_rng(42)

    print("=" * 96)
    print("  TAU-CHRONO BAYESIAN NOISE TRACKER -- 2-Qubit Benchmark")
    print("=" * 96)

    # CPTP verification
    print("\n--- CPTP Verification (2-qubit channels) ---")
    from tau_chrono.channels import dephasing as deph_fn
    test_channels = [
        ("2Q Depolarizing(0.1)", two_qubit_depolarizing(0.1)),
        ("Corr Dephasing(0.1)", correlated_dephasing(0.1)),
        ("Local AD(0.1)xAD(0.1)", local_noise(amplitude_damping(0.1), amplitude_damping(0.1))),
        ("Crosstalk(0.1,0.3)", crosstalk_dephasing(0.1, 0.3)),
    ]
    for name, kraus in test_channels:
        ok = verify_cptp(kraus)
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {len(kraus)} Kraus ops")

    sigma_0 = make_2q_sigma([0.4, 0.3, 0.2, 0.1])

    results = [
        circuit_A_bell_state(rng, sigma_0),
        circuit_B_swap(rng, sigma_0),
        circuit_C_random_2q(rng, sigma_0),
        circuit_D_variational(rng, sigma_0),
    ]

    for r in results:
        print_circuit_result(r)

    print_summary_table(results)

    # Verification
    print(f"\n{'=' * 96}")
    print("  VERIFICATION")
    print(f"{'=' * 96}")
    t1 = all(r.composition_holds for r in results)
    print(f"  [{'PASS' if t1 else 'FAIL'}] All composition inequalities hold")
    t2 = all(r.improvement_pct > 0 for r in results)
    print(f"  [{'PASS' if t2 else 'FAIL'}] All circuits show positive improvement")
    max_imp = max(r.improvement_pct for r in results)
    t3 = max_imp >= 15
    print(f"  [{'PASS' if t3 else 'FAIL'}] >=15% improvement achieved (max={max_imp:.2f}%)")
    print(f"{'=' * 96}")


if __name__ == "__main__":
    main()
