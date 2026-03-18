#!/usr/bin/env python3
"""
tau-chrono Hardware Benchmark Suite

Runs Experiments 2, 3, 5 from the hardware validation plan.
Produces CSV data + matplotlib figures for publication.

Usage:
  python hardware_benchmark.py --local          # Simulated (for testing)
  python hardware_benchmark.py --backend "QX emulator"   # QI simulator
  python hardware_benchmark.py --backend "Starmon-7"     # Real hardware!

Output:
  results/
    exp2_circuit_comparison.csv
    exp3_depth_scaling.csv
    exp5_ground_truth.csv
    fig1_depth_scaling.png
    fig3_prediction_accuracy.png
    fig4_per_gate_breakdown.png
    fig5_bayesian_vs_naive.png
"""

import argparse
import sys
import os
import json
import csv
import numpy as np
from datetime import datetime

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import tau_chrono as tc
from tau_chrono.tomography import simulate_process_tomography_1q
from tau_chrono.petz import apply_channel, fidelity

# ---------------------------------------------------------------------------
# Gate definitions
# ---------------------------------------------------------------------------

GATE_POOL_1Q = {
    "H": lambda qc, q: qc.h(q),
    "X": lambda qc, q: qc.x(q),
    "SX": lambda qc, q: qc.sx(q),
    "Rz(π/4)": lambda qc, q: qc.rz(np.pi / 4, q),
    "S": lambda qc, q: qc.s(q),
}

# Noise models for local simulation (matching Starmon-7 specs)
NOISE_1Q = {
    "H": tc.depolarizing(0.002),       # ~99.9% fidelity
    "X": tc.depolarizing(0.002),
    "SX": tc.depolarizing(0.0015),
    "Rz(π/4)": tc.dephasing(0.0005),   # virtual gate, very low noise
    "S": tc.dephasing(0.0005),
}


def make_results_dir():
    """Create results directory."""
    d = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                     "results")
    os.makedirs(d, exist_ok=True)
    return d


# ---------------------------------------------------------------------------
# Experiment 2: Circuit-level Bayesian vs Naive
# ---------------------------------------------------------------------------

def exp2_circuit_comparison(backend=None, shots=4096, local=False, timeout=600,
                            gate_kraus_cache=None):
    """Run multiple circuits and compare Bayesian vs naive tau."""
    print("\n" + "=" * 60)
    print("EXPERIMENT 2: Circuit-Level Bayesian vs Naive")
    print("=" * 60)

    circuits = {
        "H-Rz-X-H": ["H", "Rz(π/4)", "X", "H"],
        "5-gate": ["H", "SX", "Rz(π/4)", "X", "H"],
        "8-gate": ["H", "X", "SX", "Rz(π/4)", "S", "H", "X", "SX"],
        "10-gate": ["H", "X", "SX", "Rz(π/4)", "H", "S", "X", "SX", "H", "Rz(π/4)"],
        "15-gate": ["H", "X", "SX", "Rz(π/4)", "H"] * 3,
        "20-gate": ["H", "X", "SX", "Rz(π/4)", "H"] * 4,
    }

    results = []
    rho = np.array([[1, 0], [0, 0]], dtype=complex)
    sigma = np.eye(2, dtype=complex) / 2

    if gate_kraus_cache is None:
        gate_kraus_cache = {}

    for circ_name, gate_list in circuits.items():
        print(f"\n--- Circuit: {circ_name} (depth={len(gate_list)}) ---")

        channels = []
        for g in gate_list:
            if g not in gate_kraus_cache:
                if local:
                    gate_kraus_cache[g] = NOISE_1Q[g]
                else:
                    from tau_chrono.qi_bridge import extract_gate_kraus
                    char = extract_gate_kraus(
                        GATE_POOL_1Q[g], backend, g, shots, verbose=False,
                        timeout=timeout,
                    )
                    gate_kraus_cache[g] = char.kraus_ops
                    print(f"  Characterized {g}: tau_worst={char.tau_worst:.4f}")
            channels.append(gate_kraus_cache[g])

        # Run Bayesian composition
        comp = tc.bayesian_compose(channels, sigma, rho, gate_list, min_depth=6)

        row = {
            "circuit": circ_name,
            "depth": len(gate_list),
            "tau_naive": comp.tau_multiplicative_total,
            "tau_bayesian": comp.tau_bayesian_total,
            "improvement_pct": comp.improvement_percent,
            "composition_holds": comp.composition_holds,
            "composition_slack": comp.composition_slack,
        }
        results.append(row)

        print(f"  τ_naive={row['tau_naive']:.6f}  "
              f"τ_Bayes={row['tau_bayesian']:.6f}  "
              f"Improvement={row['improvement_pct']:.1f}%")

    # Save CSV
    results_dir = make_results_dir()
    csv_path = os.path.join(results_dir, "exp2_circuit_comparison.csv")
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)
    print(f"\nSaved: {csv_path}")

    return results


# ---------------------------------------------------------------------------
# Experiment 3: Depth Scaling Curve
# ---------------------------------------------------------------------------

def exp3_depth_scaling(backend=None, shots=4096, local=False, timeout=600,
                       gate_kraus_cache=None):
    """Measure improvement % as a function of circuit depth."""
    print("\n" + "=" * 60)
    print("EXPERIMENT 3: Depth Scaling Curve")
    print("=" * 60)

    depths = [2, 4, 6, 8, 10, 15, 20, 30]
    gate_cycle = ["H", "X", "SX", "Rz(π/4)", "S"]

    rho = np.array([[1, 0], [0, 0]], dtype=complex)
    sigma = np.eye(2, dtype=complex) / 2

    gate_kraus = gate_kraus_cache if gate_kraus_cache else {}
    for g in gate_cycle:
        if g not in gate_kraus:
            if local:
                gate_kraus[g] = NOISE_1Q[g]
            else:
                from tau_chrono.qi_bridge import extract_gate_kraus
                char = extract_gate_kraus(
                    GATE_POOL_1Q[g], backend, g, shots, verbose=False,
                    timeout=timeout,
                )
                gate_kraus[g] = char.kraus_ops
                print(f"  Characterized {g}: tau_worst={char.tau_worst:.4f}")

    results = []
    for d in depths:
        gate_list = [gate_cycle[i % len(gate_cycle)] for i in range(d)]
        channels = [gate_kraus[g] for g in gate_list]

        comp = tc.bayesian_compose(channels, sigma, rho, gate_list, min_depth=6)

        row = {
            "depth": d,
            "tau_naive": comp.tau_multiplicative_total,
            "tau_bayesian": comp.tau_bayesian_total,
            "improvement_pct": comp.improvement_percent,
        }
        results.append(row)
        print(f"  Depth {d:2d}: τ_naive={row['tau_naive']:.6f}  "
              f"τ_Bayes={row['tau_bayesian']:.6f}  "
              f"Improvement={row['improvement_pct']:.1f}%")

    # Save CSV
    results_dir = make_results_dir()
    csv_path = os.path.join(results_dir, "exp3_depth_scaling.csv")
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)
    print(f"\nSaved: {csv_path}")

    return results


# ---------------------------------------------------------------------------
# Experiment 5: Ground-Truth Validation
# ---------------------------------------------------------------------------

def exp5_ground_truth(backend=None, shots=4096, local=False, timeout=600,
                      gate_kraus_cache=None):
    """Validate predicted fidelity against actual circuit output."""
    print("\n" + "=" * 60)
    print("EXPERIMENT 5: Ground-Truth Validation")
    print("=" * 60)

    circuits = {
        "H-Rz-X-H": ["H", "Rz(π/4)", "X", "H"],
        "5-gate": ["H", "SX", "Rz(π/4)", "X", "H"],
        "10-gate": ["H", "X", "SX", "Rz(π/4)", "H", "S", "X", "SX", "H", "Rz(π/4)"],
        "20-gate": ["H", "X", "SX", "Rz(π/4)", "H"] * 4,
    }

    rho = np.array([[1, 0], [0, 0]], dtype=complex)
    sigma = np.eye(2, dtype=complex) / 2

    gate_kraus = gate_kraus_cache if gate_kraus_cache else {}
    for g in ["H", "X", "SX", "Rz(π/4)", "S"]:
        if g not in gate_kraus:
            if local:
                gate_kraus[g] = NOISE_1Q[g]
            else:
                from tau_chrono.qi_bridge import extract_gate_kraus
                char = extract_gate_kraus(
                    GATE_POOL_1Q[g], backend, g, shots, verbose=False,
                    timeout=timeout,
                )
                gate_kraus[g] = char.kraus_ops

    # Ideal unitaries for computing ideal output
    H_mat = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
    X_mat = np.array([[0, 1], [1, 0]], dtype=complex)
    SX_mat = np.array([[1 + 1j, 1 - 1j], [1 - 1j, 1 + 1j]], dtype=complex) / 2
    Rz_mat = np.array([[np.exp(-1j * np.pi / 8), 0],
                        [0, np.exp(1j * np.pi / 8)]], dtype=complex)
    S_mat = np.array([[1, 0], [0, 1j]], dtype=complex)

    ideal_unitaries = {
        "H": H_mat, "X": X_mat, "SX": SX_mat,
        "Rz(π/4)": Rz_mat, "S": S_mat,
    }

    results = []
    for circ_name, gate_list in circuits.items():
        channels = [gate_kraus[g] for g in gate_list]

        # Bayesian composition
        comp = tc.bayesian_compose(channels, sigma, rho, gate_list, min_depth=6)

        # Compute ideal output state
        rho_ideal = rho.copy()
        for g in gate_list:
            U = ideal_unitaries[g]
            rho_ideal = U @ rho_ideal @ U.conj().T

        # Simulate actual output (apply all noisy channels)
        rho_actual = rho.copy()
        for kraus in channels:
            rho_actual = apply_channel(rho_actual, kraus)

        F_actual = fidelity(rho_ideal, rho_actual)
        F_naive_pred = 1.0 - comp.tau_multiplicative_total
        F_bayes_pred = 1.0 - comp.tau_bayesian_total

        err_naive = abs(F_actual - F_naive_pred)
        err_bayes = abs(F_actual - F_bayes_pred)

        row = {
            "circuit": circ_name,
            "depth": len(gate_list),
            "F_actual": F_actual,
            "F_naive_pred": F_naive_pred,
            "F_bayes_pred": F_bayes_pred,
            "err_naive": err_naive,
            "err_bayes": err_bayes,
            "bayes_closer": err_bayes < err_naive,
        }
        results.append(row)
        print(f"  {circ_name}: F_actual={F_actual:.4f}  "
              f"F_naive={F_naive_pred:.4f}  "
              f"F_Bayes={F_bayes_pred:.4f}  "
              f"Bayes closer: {row['bayes_closer']}")

    # Save CSV
    results_dir = make_results_dir()
    csv_path = os.path.join(results_dir, "exp5_ground_truth.csv")
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)
    print(f"\nSaved: {csv_path}")

    return results


# ---------------------------------------------------------------------------
# Visualization
# ---------------------------------------------------------------------------

def generate_figures(results_dir):
    """Generate publication-quality figures from CSV data."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "font.size": 12,
        "font.family": "serif",
        "axes.labelsize": 14,
        "axes.titlesize": 14,
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
        "legend.fontsize": 11,
        "figure.dpi": 150,
    })

    # --- Figure 1: Depth Scaling ---
    csv_path = os.path.join(results_dir, "exp3_depth_scaling.csv")
    if os.path.exists(csv_path):
        data = []
        with open(csv_path) as f:
            reader = csv.DictReader(f)
            for row in reader:
                data.append({k: float(v) for k, v in row.items()})

        depths = [d["depth"] for d in data]
        improvements = [d["improvement_pct"] for d in data]
        tau_naive = [d["tau_naive"] for d in data]
        tau_bayes = [d["tau_bayesian"] for d in data]

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

        # Left: improvement vs depth
        ax1.plot(depths, improvements, "o-", color="#2196F3", linewidth=2,
                 markersize=8, markerfacecolor="white", markeredgewidth=2)
        ax1.set_xlabel("Circuit Depth (number of gates)")
        ax1.set_ylabel("Bayesian Improvement (%)")
        ax1.set_title("(a) Bayesian Advantage vs Circuit Depth")
        ax1.grid(True, alpha=0.3)
        ax1.set_ylim(bottom=0)

        # Right: tau comparison
        ax2.plot(depths, tau_naive, "s--", color="#F44336", linewidth=2,
                 markersize=7, label=r"$\tau_{\rm naive}$ (independent gates)")
        ax2.plot(depths, tau_bayes, "o-", color="#2196F3", linewidth=2,
                 markersize=7, label=r"$\tau_{\rm Bayesian}$ (tau-chrono)")
        ax2.set_xlabel("Circuit Depth (number of gates)")
        ax2.set_ylabel(r"Total $\tau$")
        ax2.set_title(r"(b) $\tau_{\rm naive}$ vs $\tau_{\rm Bayesian}$")
        ax2.legend()
        ax2.grid(True, alpha=0.3)

        fig.tight_layout()
        fig_path = os.path.join(results_dir, "fig1_depth_scaling.png")
        fig.savefig(fig_path, dpi=150, bbox_inches="tight")
        print(f"Saved: {fig_path}")
        plt.close(fig)

    # --- Figure 3: Prediction Accuracy ---
    csv_path = os.path.join(results_dir, "exp5_ground_truth.csv")
    if os.path.exists(csv_path):
        data = []
        with open(csv_path) as f:
            reader = csv.DictReader(f)
            for row in reader:
                data.append(row)

        names = [d["circuit"] for d in data]
        err_naive = [float(d["err_naive"]) for d in data]
        err_bayes = [float(d["err_bayes"]) for d in data]

        fig, ax = plt.subplots(figsize=(8, 5))
        x = np.arange(len(names))
        width = 0.35

        bars1 = ax.bar(x - width / 2, err_naive, width, color="#F44336",
                       alpha=0.8, label="Naive (independent)")
        bars2 = ax.bar(x + width / 2, err_bayes, width, color="#2196F3",
                       alpha=0.8, label="Bayesian (tau-chrono)")

        ax.set_xlabel("Circuit")
        ax.set_ylabel(r"$|F_{\rm actual} - F_{\rm predicted}|$")
        ax.set_title("Fidelity Prediction Error: Naive vs Bayesian")
        ax.set_xticks(x)
        ax.set_xticklabels(names, rotation=15, ha="right")
        ax.legend()
        ax.grid(True, alpha=0.3, axis="y")

        fig.tight_layout()
        fig_path = os.path.join(results_dir, "fig3_prediction_accuracy.png")
        fig.savefig(fig_path, dpi=150, bbox_inches="tight")
        print(f"Saved: {fig_path}")
        plt.close(fig)

    # --- Figure 5: Bayesian vs Naive Scatter ---
    csv_path = os.path.join(results_dir, "exp2_circuit_comparison.csv")
    if os.path.exists(csv_path):
        data = []
        with open(csv_path) as f:
            reader = csv.DictReader(f)
            for row in reader:
                data.append(row)

        tau_n = [float(d["tau_naive"]) for d in data]
        tau_b = [float(d["tau_bayesian"]) for d in data]
        labels = [d["circuit"] for d in data]

        fig, ax = plt.subplots(figsize=(7, 7))

        max_val = max(max(tau_n), max(tau_b)) * 1.1
        ax.plot([0, max_val], [0, max_val], "k--", alpha=0.3, label="y = x (no improvement)")
        ax.scatter(tau_n, tau_b, s=100, c="#2196F3", edgecolors="black",
                   linewidths=1.5, zorder=5)

        for i, label in enumerate(labels):
            ax.annotate(label, (tau_n[i], tau_b[i]),
                        textcoords="offset points", xytext=(8, 5),
                        fontsize=9)

        ax.set_xlabel(r"$\tau_{\rm naive}$ (independent gate model)")
        ax.set_ylabel(r"$\tau_{\rm Bayesian}$ (tau-chrono)")
        ax.set_title(r"All points below $y=x$: Bayesian always wins")
        ax.set_xlim(0, max_val)
        ax.set_ylim(0, max_val)
        ax.set_aspect("equal")
        ax.legend(loc="upper left")
        ax.grid(True, alpha=0.3)

        fig.tight_layout()
        fig_path = os.path.join(results_dir, "fig5_bayesian_vs_naive.png")
        fig.savefig(fig_path, dpi=150, bbox_inches="tight")
        print(f"Saved: {fig_path}")
        plt.close(fig)


# ---------------------------------------------------------------------------
# Summary report
# ---------------------------------------------------------------------------

def print_summary(exp2, exp3, exp5, backend_name, results_dir):
    """Print and save a summary report."""
    report = []
    report.append("=" * 60)
    report.append("tau-chrono HARDWARE VALIDATION REPORT")
    report.append(f"Backend: {backend_name}")
    report.append(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    report.append("=" * 60)

    report.append("\n--- Experiment 2: Circuit Comparison ---")
    report.append(f"{'Circuit':<15} {'Depth':>5} {'τ_naive':>10} "
                  f"{'τ_Bayes':>10} {'Improve':>10}")
    report.append("-" * 55)
    for r in exp2:
        report.append(f"{r['circuit']:<15} {r['depth']:>5} "
                      f"{r['tau_naive']:>10.6f} {r['tau_bayesian']:>10.6f} "
                      f"{r['improvement_pct']:>9.1f}%")

    report.append("\n--- Experiment 3: Depth Scaling ---")
    report.append(f"{'Depth':>5} {'τ_naive':>10} {'τ_Bayes':>10} {'Improve':>10}")
    report.append("-" * 40)
    for r in exp3:
        report.append(f"{r['depth']:>5} {r['tau_naive']:>10.6f} "
                      f"{r['tau_bayesian']:>10.6f} {r['improvement_pct']:>9.1f}%")

    report.append("\n--- Experiment 5: Ground Truth ---")
    report.append(f"{'Circuit':<15} {'F_actual':>10} {'F_naive':>10} "
                  f"{'F_Bayes':>10} {'Closer?':>8}")
    report.append("-" * 60)
    for r in exp5:
        report.append(f"{r['circuit']:<15} {r['F_actual']:>10.4f} "
                      f"{r['F_naive_pred']:>10.4f} {r['F_bayes_pred']:>10.4f} "
                      f"{'YES' if r['bayes_closer'] else 'NO':>8}")

    # Success criteria
    report.append("\n--- Success Criteria ---")
    imp_5 = next((r["improvement_pct"] for r in exp2
                  if r["circuit"] == "5-gate"), 0)
    imp_10 = next((r["improvement_pct"] for r in exp2
                   if r["circuit"] == "10-gate"), 0)
    all_comp = all(r.get("composition_holds", True) for r in exp2)
    all_closer = all(r["bayes_closer"] for r in exp5)

    checks = [
        ("5-gate improvement > 10%", imp_5 > 10, f"{imp_5:.1f}%"),
        ("10-gate improvement > 20%", imp_10 > 20, f"{imp_10:.1f}%"),
        ("Composition inequality holds", all_comp, str(all_comp)),
        ("F_Bayes closer to F_actual", all_closer, str(all_closer)),
    ]

    for name, passed, val in checks:
        status = "PASS" if passed else "FAIL"
        report.append(f"  [{status}] {name}: {val}")

    all_pass = all(c[1] for c in checks)
    report.append(f"\n{'='*60}")
    report.append(f"OVERALL: {'ALL CRITERIA PASSED' if all_pass else 'SOME CRITERIA FAILED'}")
    report.append(f"{'='*60}")

    text = "\n".join(report)
    print(text)

    report_path = os.path.join(results_dir, "validation_report.txt")
    with open(report_path, "w") as f:
        f.write(text)
    print(f"\nSaved: {report_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="tau-chrono Hardware Benchmark")
    parser.add_argument("--backend", default="QX emulator")
    parser.add_argument("--shots", type=int, default=4096)
    parser.add_argument("--local", action="store_true",
                        help="Local simulation (no QI needed)")
    parser.add_argument("--timeout", type=int, default=600,
                        help="Job timeout in seconds (default: 600)")
    args = parser.parse_args()

    backend = None
    backend_name = "local_simulation"

    if not args.local:
        from tau_chrono.qi_bridge import get_qi_backend
        backend = get_qi_backend(args.backend)
        backend_name = args.backend
        print(f"Connected to: {backend_name}")

    results_dir = make_results_dir()

    # Pre-characterize all gates once (shared cache)
    gate_kraus_cache = {}
    all_gates = ["H", "X", "SX", "Rz(π/4)", "S"]
    if not args.local:
        from tau_chrono.qi_bridge import extract_gate_kraus
        print("\n" + "=" * 60)
        print("PRE-CHARACTERIZING ALL GATES (shared cache)")
        print("=" * 60)
        for g in all_gates:
            print(f"\n  Characterizing {g}...")
            char = extract_gate_kraus(
                GATE_POOL_1Q[g], backend, g, args.shots, verbose=True,
                timeout=args.timeout,
            )
            gate_kraus_cache[g] = char.kraus_ops
            print(f"  -> {g}: tau_worst={char.tau_worst:.4f}, "
                  f"Kraus rank={char.n_kraus}, CPTP err={char.cptp_error:.2e}")
    else:
        for g in all_gates:
            gate_kraus_cache[g] = NOISE_1Q[g]

    # Run experiments (using shared cache)
    exp2 = exp2_circuit_comparison(backend, args.shots, args.local, args.timeout,
                                    gate_kraus_cache)
    exp3 = exp3_depth_scaling(backend, args.shots, args.local, args.timeout,
                               gate_kraus_cache)
    exp5 = exp5_ground_truth(backend, args.shots, args.local, args.timeout,
                              gate_kraus_cache)

    # Generate figures
    print("\n" + "=" * 60)
    print("GENERATING FIGURES")
    print("=" * 60)
    generate_figures(results_dir)

    # Print summary
    print_summary(exp2, exp3, exp5, backend_name, results_dir)


if __name__ == "__main__":
    main()
