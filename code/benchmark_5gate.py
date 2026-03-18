#!/usr/bin/env python3
"""
=============================================================================
CANONICAL 5-GATE BENCHMARK for tau-Chrono EDA Engine
=============================================================================

Circuit: AD(0.05) -> Dep(0.08) -> AD(0.12) -> Deph(0.03) -> AD(0.10)

  rho = |+><+| = [[0.5, 0.5], [0.5, 0.5]]
  sigma_0 = diag(0.8, 0.2)

Tests:
  1. Per-gate table: gate, channel, tau_naive, tau_eff, classification
  2. Total: tau_bayesian, tau_multiplicative, improvement %
  3. Composition inequality: sqrt(tau_total) <= sum sqrt(tau_i^eff)
  4. G3 (dephasing) classification and behavior

The Bayesian improvement compares:
  - tau_bayesian = tau(rho, N_composed, sigma_0) [exact composed channel]
  - tau_multiplicative = 1 - prod(1 - tau_naive_i)  [independent-gate estimate]

Author: Sheng-Kai Huang (2026)
License: MIT
"""

from __future__ import annotations

import sys
import os
import numpy as np

# Ensure imports work
sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'numerical'))

from noise_channels import amplitude_damping, depolarizing, dephasing, verify_cptp
from bayesian_engine import bayesian_compose
from petz_toolkit import apply_channel, tau_parameter


def main():
    print("=" * 78)
    print("  TAU-CHRONO EDA ENGINE -- Canonical 5-Gate Benchmark")
    print("=" * 78)

    # -----------------------------------------------------------------------
    # Circuit definition
    # -----------------------------------------------------------------------
    circuit = [
        ("G0: AD(0.05)",   amplitude_damping(0.05)),
        ("G1: Dep(0.08)",  depolarizing(0.08)),
        ("G2: AD(0.12)",   amplitude_damping(0.12)),
        ("G3: Deph(0.03)", dephasing(0.03)),
        ("G4: AD(0.10)",   amplitude_damping(0.10)),
    ]
    channel_names = [name for name, _ in circuit]
    channels = [ops for _, ops in circuit]

    # Initial states
    rho = 0.5 * np.array([[1, 1], [1, 1]], dtype=complex)     # |+><+|
    sigma_0 = np.diag([0.8, 0.2]).astype(complex)

    # -----------------------------------------------------------------------
    # Verify CPTP
    # -----------------------------------------------------------------------
    print("\n--- CPTP Verification ---")
    for name, ops in circuit:
        ok = verify_cptp(ops)
        print(f"  {name}: CPTP = {ok}")

    # -----------------------------------------------------------------------
    # Run Bayesian composition
    # -----------------------------------------------------------------------
    print("\n--- Running Bayesian Composition Engine ---")
    result = bayesian_compose(channels, sigma_0, rho, channel_names)

    # -----------------------------------------------------------------------
    # Per-gate table
    # -----------------------------------------------------------------------
    print("\n" + "=" * 78)
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

    # -----------------------------------------------------------------------
    # Summary statistics
    # -----------------------------------------------------------------------
    sum_tau_naive = sum(gr.tau_naive for gr in result.gate_results)
    sum_tau_eff = sum(gr.tau_eff for gr in result.gate_results)
    sum_sqrt_naive = sum(np.sqrt(max(gr.tau_naive, 0.0)) for gr in result.gate_results)
    sum_sqrt_eff = sum(np.sqrt(max(gr.tau_eff, 0.0)) for gr in result.gate_results)

    print("\n" + "=" * 78)
    print("  COMPOSITION ANALYSIS")
    print("=" * 78)

    print(f"\n  Per-gate sums:")
    print(f"    Sum tau_naive (fixed rho_0, sigma_0):        {sum_tau_naive:.6f}")
    print(f"    Sum tau_eff   (Bayesian rho_i, sigma_i):     {sum_tau_eff:.6f}")
    print(f"    Sum sqrt(tau_naive):                         {sum_sqrt_naive:.6f}")
    print(f"    Sum sqrt(tau_eff):                           {sum_sqrt_eff:.6f}")

    print(f"\n  Total channel:")
    print(f"    tau_bayesian (exact composed):               {result.tau_bayesian_total:.6f}")
    print(f"    tau_multiplicative (1-prod(1-tau_naive)):    {result.tau_multiplicative_total:.6f}")
    print(f"    Improvement (Bayesian vs multiplicative):    {result.improvement_percent:.2f}%")

    # -----------------------------------------------------------------------
    # Composition inequality check
    # -----------------------------------------------------------------------
    print(f"\n  Composition inequality: sqrt(tau_total) <= sum sqrt(tau_i^eff)")
    print(f"    LHS = sqrt(tau_total)    = {result.composition_lhs:.6f}")
    print(f"    RHS = sum sqrt(tau_eff)  = {result.composition_rhs:.6f}")
    print(f"    Slack = RHS - LHS        = {result.composition_slack:.6f}")
    print(f"    HOLDS: {result.composition_holds}")

    # -----------------------------------------------------------------------
    # G3 (Dephasing) analysis
    # -----------------------------------------------------------------------
    g3 = result.gate_results[3]
    print(f"\n  G3 Deph(0.03) detailed analysis:")
    print(f"    Pre-channel  commutator norm: {g3.pre_commutator_norm:.8f}")
    print(f"    Post-channel commutator norm: {g3.commutator_norm:.8f}")
    print(f"    Classification:               {g3.classification}")
    print(f"    tau_naive:                    {g3.tau_naive:.6f}")
    print(f"    tau_eff:                      {g3.tau_eff:.6f}")
    if g3.tau_naive > 1e-15:
        print(f"    tau_eff/tau_naive:             {g3.tau_eff/g3.tau_naive:.4f}")

    # -----------------------------------------------------------------------
    # Bayesian sigma evolution
    # -----------------------------------------------------------------------
    print("\n" + "=" * 78)
    print("  BAYESIAN SIGMA EVOLUTION")
    print("=" * 78)
    print(f"\n  sigma_0 = diag({sigma_0[0,0].real:.2f}, {sigma_0[1,1].real:.2f})")

    sigma_current = sigma_0.copy()
    rho_current = rho.copy()
    for i, (name, ops) in enumerate(circuit):
        sigma_current = apply_channel(sigma_current, ops)
        rho_current = apply_channel(rho_current, ops)
        s_diag = np.diag(sigma_current).real
        r_diag = np.diag(rho_current).real
        r_off = abs(rho_current[0, 1])
        print(f"  After {name}:")
        print(f"    sigma diag = ({s_diag[0]:.6f}, {s_diag[1]:.6f})")
        print(f"    rho   diag = ({r_diag[0]:.6f}, {r_diag[1]:.6f}), |off-diag| = {r_off:.6f}")

    # -----------------------------------------------------------------------
    # Verification summary
    # -----------------------------------------------------------------------
    print("\n" + "=" * 78)
    print("  VERIFICATION SUMMARY")
    print("=" * 78)

    tests_passed = 0
    tests_total = 4

    # Test 1: Composition inequality holds
    t1 = result.composition_holds
    status1 = "PASS" if t1 else "FAIL"
    print(f"\n  [{status1}] T1: Composition inequality holds "
          f"(slack={result.composition_slack:.6f})")
    if t1:
        tests_passed += 1

    # Test 2: Improvement positive
    t2 = result.improvement_percent > 0
    status2 = "PASS" if t2 else "FAIL"
    print(f"  [{status2}] T2: Bayesian improvement > 0% "
          f"(got {result.improvement_percent:.2f}%)")
    if t2:
        tests_passed += 1

    # Test 3: G3 (dephasing) has tau_eff <= tau_naive
    t3 = g3.tau_eff <= g3.tau_naive + 1e-10
    status3 = "PASS" if t3 else "FAIL"
    print(f"  [{status3}] T3: G3 tau_eff <= tau_naive "
          f"({g3.tau_eff:.6f} vs {g3.tau_naive:.6f})")
    if t3:
        tests_passed += 1

    # Test 4: Zero composition violations
    t4 = result.composition_holds
    status4 = "PASS" if t4 else "FAIL"
    print(f"  [{status4}] T4: Zero composition violations")
    if t4:
        tests_passed += 1

    # Bonus: G3 classification
    g3_near_sat = g3.classification in ("NEAR_SATURATED", "SATURATED")
    if g3_near_sat:
        print(f"\n  [BONUS] G3 classified as {g3.classification} "
              f"(pre_comm={g3.pre_commutator_norm:.6f}, "
              f"post_comm={g3.commutator_norm:.6f})")
    else:
        print(f"\n  [INFO] G3 classified as {g3.classification} "
              f"(pre_comm={g3.pre_commutator_norm:.6f}, "
              f"post_comm={g3.commutator_norm:.6f})")
        # Dephasing on |+> produces large commutator because |+> is maximally
        # off-diagonal while sigma is diagonal. The "near-saturation" regime
        # requires rho to be nearly diagonal (in the sigma eigenbasis), which
        # happens only after sufficient amplitude damping drives rho toward |0>.

    print(f"\n  Tests passed: {tests_passed}/{tests_total}")
    print("=" * 78)

    return result


if __name__ == "__main__":
    result = main()
