#!/usr/bin/env python3
"""
tau-chrono + Quantum Inspire Demo

End-to-end demonstration:
  1. Connect to Quantum Inspire (simulator or real hardware)
  2. Run process tomography to extract real Kraus operators
  3. Feed into tau-chrono Bayesian composition engine
  4. Compare Bayesian vs naive noise tracking

Usage:
  # Step 1: Install dependencies
  pip install quantuminspire qiskit-quantuminspire

  # Step 2: Login (one-time, opens browser)
  qi login

  # Step 3: Run this script
  python qi_demo.py                    # QX emulator (simulator)
  python qi_demo.py --backend Starmon-7  # Real hardware!
  python qi_demo.py --local            # Local simulation (no QI needed)

Author: Sheng-Kai Huang
"""

import argparse
import sys
import os
import numpy as np

# Add parent dir to path for tau_chrono
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import tau_chrono as tc
from tau_chrono.tomography import simulate_process_tomography_1q


def run_local_demo():
    """Run demo with simulated tomography (no hardware needed)."""
    print("=" * 60)
    print("tau-chrono QDA Demo (Local Simulation)")
    print("=" * 60)

    # Define "true" noise channels (simulating what hardware would give)
    noise_models = {
        "Amplitude Damping (γ=0.05)": tc.amplitude_damping(0.05),
        "Depolarizing (p=0.03)": tc.depolarizing(0.03),
        "Dephasing (p=0.04)": tc.dephasing(0.04),
    }

    print("\n--- Step 1: Simulated Process Tomography ---")
    print("Extracting Kraus operators via simulated measurements...\n")

    reconstructed = {}
    for name, true_kraus in noise_models.items():
        print(f"\nChannel: {name}")
        kraus_recon = simulate_process_tomography_1q(
            true_kraus, shots=100000, verbose=True
        )
        reconstructed[name] = kraus_recon

        # Compare reconstruction quality
        sigma = np.eye(2, dtype=complex) / 2
        rho = np.array([[1, 0], [0, 0]], dtype=complex)

        tau_true = tc.tau_parameter(rho, true_kraus, sigma)
        tau_recon = tc.tau_parameter(rho, kraus_recon, sigma)
        print(f"  tau (true):          {tau_true:.6f}")
        print(f"  tau (reconstructed): {tau_recon:.6f}")
        print(f"  Relative error:      {abs(tau_true - tau_recon) / max(tau_true, 1e-10) * 100:.1f}%")

    # Step 2: Bayesian composition with reconstructed channels
    print("\n" + "=" * 60)
    print("--- Step 2: Bayesian Composition Analysis ---")
    print("=" * 60)

    # Simulate a 5-gate circuit using reconstructed noise
    channel_sequence = [
        reconstructed["Amplitude Damping (γ=0.05)"],
        reconstructed["Depolarizing (p=0.03)"],
        reconstructed["Dephasing (p=0.04)"],
        reconstructed["Amplitude Damping (γ=0.05)"],
        reconstructed["Depolarizing (p=0.03)"],
    ]
    channel_names = [
        "H + AD", "CNOT + Dep", "Rz + Deph", "X + AD", "CNOT + Dep"
    ]

    rho = np.array([[1, 0], [0, 0]], dtype=complex)
    sigma = np.eye(2, dtype=complex) / 2

    result = tc.bayesian_compose(
        channels=channel_sequence,
        sigma_0=sigma,
        rho=rho,
        channel_names=channel_names,
    )

    print(f"\nCircuit: {' → '.join(channel_names)}")
    print(f"\nPer-gate analysis:")
    print(f"{'Gate':<15} {'tau_naive':>10} {'tau_eff':>10} {'Class':>16}")
    print("-" * 55)
    for gr in result.gate_results:
        print(f"{gr.channel_name:<15} {gr.tau_naive:>10.6f} "
              f"{gr.tau_eff:>10.6f} {gr.classification:>16}")

    print(f"\nTotal results:")
    print(f"  tau (multiplicative/naive): {result.tau_multiplicative_total:.6f}")
    print(f"  tau (Bayesian):             {result.tau_bayesian_total:.6f}")
    print(f"  Improvement:                {result.improvement_percent:.1f}%")
    print(f"  Composition inequality:     {'HOLDS ✓' if result.composition_holds else 'VIOLATED ✗'}")
    print(f"  Slack:                       {result.composition_slack:.6f}")

    return result


def run_qi_demo(backend_name: str = "QX emulator", shots: int = 4096):
    """Run demo on Quantum Inspire backend."""
    from tau_chrono.qi_bridge import (
        get_qi_backend,
        list_qi_backends,
        extract_gate_kraus,
        run_tau_analysis,
    )

    print("=" * 60)
    print(f"tau-chrono QDA Demo (Quantum Inspire: {backend_name})")
    print("=" * 60)

    # Connect to backend
    print(f"\nConnecting to {backend_name}...")
    try:
        backend = get_qi_backend(backend_name)
        print(f"Connected: {backend.name}")
    except Exception as e:
        print(f"Connection failed: {e}")
        print("\nAvailable backends:")
        try:
            for name in list_qi_backends():
                print(f"  - {name}")
        except Exception:
            pass
        print("\nMake sure you have run: qi login")
        return None

    # Define gates to test
    def gate_h(qc, q):
        qc.h(q)

    def gate_x(qc, q):
        qc.x(q)

    def gate_rz(qc, q):
        qc.rz(np.pi / 4, q)  # T-gate equivalent

    def gate_sx(qc, q):
        qc.sx(q)  # sqrt(X)

    # Step 1: Characterize individual gates
    print(f"\n--- Step 1: Gate Characterization (shots={shots}) ---\n")

    gate_fns = [gate_h, gate_x, gate_rz, gate_sx]
    gate_names = ["H", "X", "Rz(π/4)", "SX"]

    gate_chars = []
    for fn, name in zip(gate_fns, gate_names):
        char = extract_gate_kraus(fn, backend, name, shots, verbose=True)
        gate_chars.append(char)

    # Step 2: Analyze a circuit
    print(f"\n--- Step 2: Circuit Analysis ---\n")

    circuit_fns = [gate_h, gate_rz, gate_x, gate_h, gate_sx]
    circuit_names = ["H", "Rz(π/4)", "X", "H", "SX"]

    analysis = run_tau_analysis(
        gate_fns=circuit_fns,
        gate_names=circuit_names,
        backend=backend,
        shots=shots,
        verbose=True,
    )

    # Step 3: Summary
    print(f"\n{'='*60}")
    print("FINAL SUMMARY")
    print(f"{'='*60}")
    print(f"Backend:     {backend_name}")
    print(f"Circuit:     {' → '.join(circuit_names)}")
    print(f"tau (naive): {analysis.tau_naive_total:.6f}")
    print(f"tau (Bayes): {analysis.tau_bayesian_total:.6f}")
    print(f"Improvement: {analysis.improvement_percent:.1f}%")
    print(f"\nThis means Bayesian noise tracking reduces estimated error by "
          f"{analysis.improvement_percent:.1f}% compared to naive independent "
          f"gate models.")

    if backend_name != "QX emulator":
        print(f"\n*** These are REAL HARDWARE results from {backend_name}! ***")

    return analysis


def main():
    parser = argparse.ArgumentParser(
        description="tau-chrono QDA Demo with Quantum Inspire"
    )
    parser.add_argument(
        "--backend", default="QX emulator",
        help="QI backend name (default: 'QX emulator')"
    )
    parser.add_argument(
        "--shots", type=int, default=4096,
        help="Shots per tomography circuit (default: 4096)"
    )
    parser.add_argument(
        "--local", action="store_true",
        help="Run local simulation only (no QI connection needed)"
    )
    parser.add_argument(
        "--list-backends", action="store_true",
        help="List available QI backends and exit"
    )

    args = parser.parse_args()

    if args.list_backends:
        from tau_chrono.qi_bridge import list_qi_backends
        print("Available Quantum Inspire backends:")
        for name in list_qi_backends():
            print(f"  - {name}")
        return

    if args.local:
        run_local_demo()
    else:
        run_qi_demo(args.backend, args.shots)


if __name__ == "__main__":
    main()
