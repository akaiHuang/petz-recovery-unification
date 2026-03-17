"""
Canonical 5-gate single-qubit benchmark for tau-chrono.

Circuit: AD(0.05) -> Dep(0.08) -> AD(0.12) -> Deph(0.03) -> AD(0.10)

    rho   = |+><+| = [[0.5, 0.5], [0.5, 0.5]]
    sigma = diag(0.8, 0.2)

Run as CLI::

    python -m tau_chrono.benchmarks.single_qubit
"""

from __future__ import annotations

import numpy as np

from tau_chrono.channels import (
    amplitude_damping,
    depolarizing,
    dephasing,
    verify_cptp,
)
from tau_chrono.bayesian import bayesian_compose
from tau_chrono.petz import apply_channel


def run() -> dict:
    """Run the canonical 5-gate benchmark and return results dict."""
    circuit = [
        ("G0: AD(0.05)",   amplitude_damping(0.05)),
        ("G1: Dep(0.08)",  depolarizing(0.08)),
        ("G2: AD(0.12)",   amplitude_damping(0.12)),
        ("G3: Deph(0.03)", dephasing(0.03)),
        ("G4: AD(0.10)",   amplitude_damping(0.10)),
    ]
    channel_names = [name for name, _ in circuit]
    channels = [ops for _, ops in circuit]

    rho = 0.5 * np.array([[1, 1], [1, 1]], dtype=complex)
    sigma_0 = np.diag([0.8, 0.2]).astype(complex)

    # Verify CPTP
    for name, ops in circuit:
        assert verify_cptp(ops), f"{name} failed CPTP check"

    result = bayesian_compose(channels, sigma_0, rho, channel_names)

    return {
        "gate_results": result.gate_results,
        "tau_bayesian": result.tau_bayesian_total,
        "tau_multiplicative": result.tau_multiplicative_total,
        "improvement_pct": result.improvement_percent,
        "composition_holds": result.composition_holds,
        "composition_slack": result.composition_slack,
    }


def main() -> None:
    """Pretty-print the 5-gate benchmark."""
    print("=" * 78)
    print("  TAU-CHRONO -- Canonical 5-Gate Single-Qubit Benchmark")
    print("=" * 78)

    circuit = [
        ("G0: AD(0.05)",   amplitude_damping(0.05)),
        ("G1: Dep(0.08)",  depolarizing(0.08)),
        ("G2: AD(0.12)",   amplitude_damping(0.12)),
        ("G3: Deph(0.03)", dephasing(0.03)),
        ("G4: AD(0.10)",   amplitude_damping(0.10)),
    ]
    channel_names = [name for name, _ in circuit]
    channels = [ops for _, ops in circuit]

    rho = 0.5 * np.array([[1, 1], [1, 1]], dtype=complex)
    sigma_0 = np.diag([0.8, 0.2]).astype(complex)

    # CPTP
    print("\n--- CPTP Verification ---")
    for name, ops in circuit:
        ok = verify_cptp(ops)
        print(f"  {name}: CPTP = {ok}")

    # Bayesian composition
    result = bayesian_compose(channels, sigma_0, rho, channel_names)

    # Per-gate table
    print(f"\n{'=' * 78}")
    print("  PER-GATE ANALYSIS")
    print("=" * 78)
    header = (f"{'Gate':<18} {'tau_naive':>10} {'tau_eff':>10} {'Ratio':>8} "
              f"{'pre_[.,.]':>10} {'post_[.,.]':>10} {'Class':<16}")
    print(header)
    print("-" * 78)

    for gr in result.gate_results:
        ratio = gr.tau_eff / gr.tau_naive if gr.tau_naive > 1e-15 else 0.0
        print(f"  {gr.channel_name:<16} {gr.tau_naive:>10.6f} {gr.tau_eff:>10.6f} "
              f"{ratio:>8.4f} {gr.pre_commutator_norm:>10.6f} "
              f"{gr.commutator_norm:>10.6f} {gr.classification:<16}")

    # Summary
    print(f"\n{'=' * 78}")
    print("  COMPOSITION ANALYSIS")
    print("=" * 78)
    print(f"\n  tau_bayesian (exact composed):            {result.tau_bayesian_total:.6f}")
    print(f"  tau_multiplicative (1-prod(1-tau_naive)): {result.tau_multiplicative_total:.6f}")
    print(f"  Improvement:                             {result.improvement_percent:.2f}%")
    print(f"\n  Composition inequality: sqrt(tau_total) <= sum sqrt(tau_i^eff)")
    print(f"    LHS = {result.composition_lhs:.6f}  <=  RHS = {result.composition_rhs:.6f}")
    print(f"    Slack = {result.composition_slack:.6f}   HOLDS: {result.composition_holds}")

    # Sigma evolution
    print(f"\n{'=' * 78}")
    print("  BAYESIAN SIGMA EVOLUTION")
    print("=" * 78)
    sigma_current = sigma_0.copy()
    rho_current = rho.copy()
    for name, ops in circuit:
        sigma_current = apply_channel(sigma_current, ops)
        rho_current = apply_channel(rho_current, ops)
        s_diag = np.diag(sigma_current).real
        r_diag = np.diag(rho_current).real
        r_off = abs(rho_current[0, 1])
        print(f"  After {name}:")
        print(f"    sigma diag = ({s_diag[0]:.6f}, {s_diag[1]:.6f})")
        print(f"    rho   diag = ({r_diag[0]:.6f}, {r_diag[1]:.6f}), |off-diag| = {r_off:.6f}")

    # Verification
    print(f"\n{'=' * 78}")
    print("  VERIFICATION SUMMARY")
    print("=" * 78)
    tests = 0
    total = 3

    t1 = result.composition_holds
    print(f"  [{'PASS' if t1 else 'FAIL'}] Composition inequality holds")
    if t1: tests += 1

    t2 = result.improvement_percent > 0
    print(f"  [{'PASS' if t2 else 'FAIL'}] Bayesian improvement > 0% "
          f"(got {result.improvement_percent:.2f}%)")
    if t2: tests += 1

    g3 = result.gate_results[3]
    t3 = g3.tau_eff <= g3.tau_naive + 1e-10
    print(f"  [{'PASS' if t3 else 'FAIL'}] G3 tau_eff <= tau_naive")
    if t3: tests += 1

    print(f"\n  Tests passed: {tests}/{total}")
    print("=" * 78)


if __name__ == "__main__":
    main()
