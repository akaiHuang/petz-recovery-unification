#!/usr/bin/env python3
"""
Analyze the Petz recovery bound τ ≤ 1 - exp(-Σ/2) across consciousness states.

This script reads the pre-computed Σ and τ values from compute_sigma_tau.py
and performs statistical tests:

1. Does the bound hold for all individual epochs?
2. How tight is the bound per stage?
3. Is there a monotonic relationship between consciousness level and Σ?
4. What is the "gap" Δ = (1-exp(-Σ/2)) - τ per stage?

Usage:
    python scripts/analyze_bound.py
"""

import os
import sys
import csv
from pathlib import Path

import numpy as np

STAGE_ORDER = ["Wake", "N1", "N2", "N3", "REM"]
STAGE_COLORS = {
    "Wake": "#e74c3c",
    "N1": "#f39c12",
    "N2": "#2ecc71",
    "N3": "#3498db",
    "REM": "#9b59b6",
}

# Consciousness level ordering (expected Σ ordering)
CONSCIOUSNESS_LEVELS = {
    "N3": 0,    # Deepest sleep = lowest consciousness
    "N2": 1,
    "N1": 2,
    "REM": 3,   # Dream state
    "Wake": 4,  # Full consciousness
}


def load_data(csv_path: Path) -> dict:
    """Load sigma-tau data from CSV."""
    data = {s: {"sigma": [], "tau": []} for s in STAGE_ORDER}

    with open(csv_path, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            stage = row["stage"]
            if stage in STAGE_ORDER:
                data[stage]["sigma"].append(float(row["sigma"]))
                data[stage]["tau"].append(float(row["tau"]))

    return data


def test_bound_violation(data: dict):
    """Test if τ ≤ 1 - exp(-Σ/2) holds for all data points."""
    print("\n" + "=" * 70)
    print("TEST 1: Petz Bound Violations")
    print("  Testing τ ≤ 1 - exp(-Σ/2) for each epoch")
    print("=" * 70)

    total_points = 0
    total_violations = 0

    for stage in STAGE_ORDER:
        sigmas = np.array(data[stage]["sigma"])
        taus = np.array(data[stage]["tau"])
        if len(sigmas) == 0:
            continue

        bounds = 1.0 - np.exp(-sigmas / 2.0)
        violations = np.sum(taus > bounds)
        n = len(sigmas)
        total_points += n
        total_violations += violations

        # Gap: how far below the bound
        gaps = bounds - taus
        mean_gap = np.mean(gaps)
        min_gap = np.min(gaps)

        status = "PASS" if violations == 0 else f"FAIL ({violations}/{n})"
        print(f"  {stage:6s}: {n:5d} points, violations: {violations:5d}, "
              f"gap_mean: {mean_gap:.4f}, gap_min: {min_gap:.4f}  [{status}]")

    viol_pct = 100 * total_violations / total_points if total_points > 0 else 0
    print(f"\n  TOTAL: {total_points} points, {total_violations} violations "
          f"({viol_pct:.2f}%)")

    if total_violations == 0:
        print("  >>> RESULT: Petz bound holds perfectly for all data points.")
    elif viol_pct < 1:
        print("  >>> RESULT: Petz bound holds for >99% of data points.")
        print("      Small violations may be due to finite-sample estimation noise.")
    else:
        print("  >>> RESULT: Significant violations detected.")
        print("      This may indicate issues with the Σ/τ estimation method.")


def test_consciousness_ordering(data: dict):
    """Test if Σ increases monotonically with consciousness level."""
    print("\n" + "=" * 70)
    print("TEST 2: Consciousness Level vs Σ Ordering")
    print("  Expected: N3 < N2 < N1 < REM < Wake")
    print("=" * 70)

    means = {}
    for stage in STAGE_ORDER:
        sigmas = data[stage]["sigma"]
        if sigmas:
            means[stage] = np.mean(sigmas)

    if len(means) < 2:
        print("  Insufficient data for ordering test.")
        return

    # Sort by mean Σ
    sorted_stages = sorted(means.keys(), key=lambda s: means[s])
    expected_order = sorted(means.keys(), key=lambda s: CONSCIOUSNESS_LEVELS[s])

    print(f"\n  Expected order (by consciousness): {' < '.join(expected_order)}")
    print(f"  Observed order (by mean Σ):        {' < '.join(sorted_stages)}")

    # Compute Spearman rank correlation
    observed_ranks = [sorted_stages.index(s) for s in STAGE_ORDER if s in means]
    expected_ranks = [expected_order.index(s) for s in STAGE_ORDER if s in means]

    n = len(observed_ranks)
    if n > 2:
        d_sq = sum((o - e) ** 2 for o, e in zip(observed_ranks, expected_ranks))
        rho = 1 - 6 * d_sq / (n * (n ** 2 - 1))
        print(f"\n  Spearman rank correlation: ρ = {rho:.3f}")
        if rho > 0.8:
            print("  >>> RESULT: Strong positive correlation with consciousness level.")
        elif rho > 0.4:
            print("  >>> RESULT: Moderate correlation with consciousness level.")
        else:
            print("  >>> RESULT: Weak/no correlation. Further investigation needed.")

    # Print mean values
    print(f"\n  Mean Σ by stage:")
    for stage in sorted_stages:
        level = CONSCIOUSNESS_LEVELS[stage]
        print(f"    {stage:6s} (level {level}): Σ = {means[stage]:.4f}")


def test_bound_tightness(data: dict):
    """Analyze how tight the Petz bound is per stage."""
    print("\n" + "=" * 70)
    print("TEST 3: Bound Tightness (Gap Analysis)")
    print("  Gap Δ = [1-exp(-Σ/2)] - τ  (smaller = tighter)")
    print("=" * 70)

    for stage in STAGE_ORDER:
        sigmas = np.array(data[stage]["sigma"])
        taus = np.array(data[stage]["tau"])
        if len(sigmas) == 0:
            continue

        bounds = 1.0 - np.exp(-sigmas / 2.0)
        gaps = bounds - taus
        saturation = taus / bounds  # How close τ is to the bound

        print(f"\n  {stage}:")
        print(f"    N = {len(sigmas)}")
        print(f"    Gap: mean = {np.mean(gaps):.4f}, "
              f"std = {np.std(gaps):.4f}, "
              f"median = {np.median(gaps):.4f}")
        print(f"    Saturation τ/bound: mean = {np.mean(saturation):.3f}, "
              f"max = {np.max(saturation):.3f}")

    print(f"\n  Interpretation:")
    print(f"    Saturation near 1.0 = bound is tight (Petz map nearly optimal)")
    print(f"    Saturation near 0.0 = bound is loose (room for better retrodiction)")
    print(f"    If N3 has higher saturation than Wake:")
    print(f"      → deep sleep is closer to perfect retrodiction")
    print(f"      → supports 'closed system' interpretation")


def generate_summary_plot(data: dict, results_dir: Path):
    """Generate the final analysis figure."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    fig.suptitle(
        r"Petz Recovery Bound $\tau \leq 1 - e^{-\Sigma/2}$"
        "\nacross Consciousness States (EEG Sleep-EDF)",
        fontsize=14, fontweight="bold"
    )

    # --- Panel A: τ vs Σ scatter ---
    ax = axes[0, 0]
    sigma_range = np.linspace(0, 2.5, 300)
    bound_curve = 1.0 - np.exp(-sigma_range / 2.0)
    ax.plot(sigma_range, bound_curve, "k-", linewidth=2, label="Petz bound", zorder=10)
    ax.fill_between(sigma_range, bound_curve, 1.0, alpha=0.08, color="red")

    for stage in STAGE_ORDER:
        s = data[stage]["sigma"]
        t = data[stage]["tau"]
        if s:
            ax.scatter(s, t, c=STAGE_COLORS[stage], label=stage,
                       alpha=0.2, s=8, edgecolors="none")
    ax.set_xlabel(r"$\Sigma$")
    ax.set_ylabel(r"$\tau$")
    ax.set_title("(A) All epochs")
    ax.legend(fontsize=8)
    ax.set_xlim(0, None)
    ax.set_ylim(0, 0.5)
    ax.grid(True, alpha=0.2)

    # --- Panel B: Mean ± std per stage ---
    ax = axes[0, 1]
    ax.plot(sigma_range, bound_curve, "k-", linewidth=2, zorder=10)
    ax.fill_between(sigma_range, bound_curve, 1.0, alpha=0.08, color="red")

    for stage in STAGE_ORDER:
        s = np.array(data[stage]["sigma"])
        t = np.array(data[stage]["tau"])
        if len(s) == 0:
            continue
        ax.errorbar(np.mean(s), np.mean(t), xerr=np.std(s), yerr=np.std(t),
                     fmt="o", markersize=14, color=STAGE_COLORS[stage],
                     label=stage, capsize=5, linewidth=2, zorder=5)
    ax.set_xlabel(r"$\langle\Sigma\rangle$")
    ax.set_ylabel(r"$\langle\tau\rangle$")
    ax.set_title("(B) Stage averages")
    ax.legend(fontsize=10)
    ax.set_xlim(0, None)
    ax.set_ylim(0, 0.5)
    ax.grid(True, alpha=0.2)

    # --- Panel C: Σ distribution by stage ---
    ax = axes[1, 0]
    stage_sigmas = [data[s]["sigma"] for s in STAGE_ORDER if data[s]["sigma"]]
    stage_labels = [s for s in STAGE_ORDER if data[s]["sigma"]]
    if stage_sigmas:
        bp = ax.boxplot(stage_sigmas, labels=stage_labels, patch_artist=True)
        for patch, stage in zip(bp["boxes"], stage_labels):
            patch.set_facecolor(STAGE_COLORS[stage])
            patch.set_alpha(0.6)
    ax.set_ylabel(r"$\Sigma$")
    ax.set_title(r"(C) $\Sigma$ distribution")
    ax.grid(True, alpha=0.2, axis="y")

    # --- Panel D: Gap distribution by stage ---
    ax = axes[1, 1]
    gaps_by_stage = []
    gap_labels = []
    for stage in STAGE_ORDER:
        s = np.array(data[stage]["sigma"])
        t = np.array(data[stage]["tau"])
        if len(s) == 0:
            continue
        gaps = (1.0 - np.exp(-s / 2.0)) - t
        gaps_by_stage.append(gaps)
        gap_labels.append(stage)

    if gaps_by_stage:
        bp = ax.boxplot(gaps_by_stage, labels=gap_labels, patch_artist=True)
        for patch, stage in zip(bp["boxes"], gap_labels):
            patch.set_facecolor(STAGE_COLORS[stage])
            patch.set_alpha(0.6)
    ax.axhline(y=0, color="red", linestyle="--", linewidth=1, alpha=0.5)
    ax.set_ylabel(r"Gap $\Delta = (1-e^{-\Sigma/2}) - \tau$")
    ax.set_title("(D) Bound tightness (Δ > 0 means bound holds)")
    ax.grid(True, alpha=0.2, axis="y")

    plt.tight_layout()
    fig_path = results_dir / "figures" / "petz_bound_analysis.png"
    fig.savefig(str(fig_path), dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"\n  Saved: {fig_path}")


def main():
    script_dir = Path(__file__).resolve().parent
    project_dir = script_dir.parent
    results_dir = project_dir / "results"
    csv_path = results_dir / "statistics" / "sigma_tau_by_stage.csv"

    print("=" * 70)
    print("Petz Recovery Bound Analysis")
    print(r"Testing: τ ≤ 1 − exp(−Σ/2) across consciousness states")
    print("=" * 70)

    if not csv_path.exists():
        print(f"\nData file not found: {csv_path}")
        print("Run compute_sigma_tau.py first.")
        sys.exit(1)

    data = load_data(csv_path)

    # Print data summary
    total = sum(len(data[s]["sigma"]) for s in STAGE_ORDER)
    print(f"\nLoaded {total} data points from {csv_path}")
    for stage in STAGE_ORDER:
        print(f"  {stage:6s}: {len(data[stage]['sigma']):5d} points")

    # Run tests
    test_bound_violation(data)
    test_consciousness_ordering(data)
    test_bound_tightness(data)

    # Generate plots
    print("\nGenerating analysis plots...")
    generate_summary_plot(data, results_dir)

    # Final summary
    print(f"\n{'=' * 70}")
    print("SUMMARY")
    print(f"{'=' * 70}")
    print("This analysis tests the Petz recovery bound from Huang (2025):")
    print("  τ ≤ 1 − exp(−Σ/2)")
    print()
    print("If the bound holds and Σ correlates with consciousness level:")
    print("  → Retrodiction quality (1-τ) is fundamentally limited by")
    print("    entropy production (Σ) in the brain")
    print("  → Deep sleep (low Σ) → good retrodiction (low τ)")
    print("  → Waking (high Σ) → poor retrodiction (high τ)")
    print("  → The 'flow of time' (irreversibility) varies with consciousness")
    print()
    print("Connection to core theory:")
    print("  Σ = D(ρ_spacetime || ρ_matter) in different limits")
    print("  Zero-entropy environment → Σ = 0 → τ = 0 → perfect retrodiction")
    print("  This EEG experiment is a biological proxy for the same physics.")


if __name__ == "__main__":
    main()
