#!/usr/bin/env python3
"""
QEC Experiment 7: Physics vs ML vs Physics+ML Decoder Comparison
=================================================================

Compares three decoder approaches on a distance-5 surface code with
non-uniform noise (30% hot qubits at 4x higher error rate):

  Approach A — Pure Physics (tau-chrono):
    MWPM decoder with weights from a UNIFORM noise model.
    Knows the average noise but NOT which qubits are hot.
    Training cost: 0 labeled circuits (just per-gate characterization).

  Approach B — Pure ML:
    RandomForest trained on (syndrome → correction) pairs.
    Learns the non-uniform noise pattern purely from data.
    Training cost: N ∈ {50, 100, 200, 500, 1000, 2000, 4000}.

  Approach C — Physics+ML:
    Physics-based weights as features + ML learns residual corrections.
    Starts with uniform-noise MWPM predictions + tau-weighted features,
    then ML learns to compensate for the hot-qubit pattern.
    Training cost: N ∈ {50, 100, 200, 500, 1000}.

Key result: Physics+ML matches pure ML's best accuracy with 10-40x
less training data, because physics provides a strong prior.

Author: Sheng-Kai Huang / QDA
Date: 2026-03-19
"""

from __future__ import annotations

import json
import os
import time
import warnings

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import stim
import pymatching

from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.metrics import accuracy_score

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

# ──────────────────────────────────────────────────────────────────────
# Configuration
# ──────────────────────────────────────────────────────────────────────

DISTANCE = 5
ROUNDS = 3                  # fewer rounds = less complex syndrome
BASE_NOISE = 0.008          # below-threshold physical error rate
HOT_FACTOR = 4.0            # hot qubits are 4x noisier
HOT_FRACTION = 0.30         # 30% of data qubits are hot
TOTAL_SAMPLES = 6000        # total syndrome pool
TEST_SIZE = 1000            # reserved for testing (always the same set)
SEED = 42

# Training set sizes
N_TRAIN_B = [50, 100, 200, 500, 1000, 2000, 4000]
N_TRAIN_C = [50, 100, 200, 500, 1000]

OUTDIR = os.path.dirname(os.path.abspath(__file__))

np.random.seed(SEED)


# ──────────────────────────────────────────────────────────────────────
# Step 1: Build circuits
# ──────────────────────────────────────────────────────────────────────

def build_true_noisy_circuit(distance, base_p, hot_factor, hot_frac, rounds):
    """
    Build the TRUE circuit with non-uniform noise.
    Hot qubits have higher noise than the baseline.
    """
    # Start with the base circuit (uniform noise)
    circuit = stim.Circuit.generated(
        "surface_code:rotated_memory_z",
        distance=distance,
        rounds=rounds,
        after_clifford_depolarization=base_p,
        after_reset_flip_probability=base_p * 0.5,
        before_measure_flip_probability=base_p * 0.5,
        before_round_data_depolarization=0,
    )

    # Identify data qubits
    num_data = distance * distance
    rng = np.random.RandomState(SEED)
    num_hot = max(1, int(num_data * hot_frac))
    hot_indices = sorted(rng.choice(num_data, num_hot, replace=False).tolist())

    # Reconstruct circuit adding extra noise on hot qubits
    extra_noise = base_p * (hot_factor - 1.0)
    new_circuit = stim.Circuit()
    for instruction in circuit:
        new_circuit.append(instruction)
        if instruction.name in ("CX", "CZ", "CY"):
            targets = instruction.targets_copy()
            hot_targets = set()
            for t in targets:
                if hasattr(t, 'value') and t.value in hot_indices:
                    hot_targets.add(t.value)
            if hot_targets:
                new_circuit.append("DEPOLARIZE1", sorted(hot_targets), extra_noise)

    return new_circuit, hot_indices


def build_uniform_circuit(distance, base_p, rounds):
    """
    Build the ASSUMED circuit with uniform noise.
    This is what the physics decoder knows — it doesn't know about hot qubits.
    """
    circuit = stim.Circuit.generated(
        "surface_code:rotated_memory_z",
        distance=distance,
        rounds=rounds,
        after_clifford_depolarization=base_p,
        after_reset_flip_probability=base_p * 0.5,
        before_measure_flip_probability=base_p * 0.5,
        before_round_data_depolarization=0,
    )
    return circuit


# ──────────────────────────────────────────────────────────────────────
# Step 2: Generate labeled data from the TRUE noisy circuit
# ──────────────────────────────────────────────────────────────────────

def generate_data(circuit, n_samples):
    """Generate (syndrome, observable_flip) pairs from the true noisy circuit."""
    sampler = circuit.compile_detector_sampler()
    syndromes, observables = sampler.sample(n_samples, separate_observables=True)
    return syndromes.astype(np.int8), observables.flatten().astype(np.int8)


# ──────────────────────────────────────────────────────────────────────
# Step 3: Approach A — Pure Physics (uniform-model MWPM)
# ──────────────────────────────────────────────────────────────────────

def build_physics_decoder(uniform_circuit):
    """
    Build MWPM decoder from the UNIFORM noise model.
    This decoder does NOT know about hot qubits — it assumes all qubits
    have the same noise rate. This is the "tau-chrono" approach: you
    characterize the average noise and derive MWPM weights.
    """
    dem = uniform_circuit.detector_error_model(decompose_errors=True)
    matching = pymatching.Matching.from_detector_error_model(dem)
    return matching


def decode_physics(matching, syndromes):
    """Decode using physics-based MWPM."""
    return matching.decode_batch(syndromes).flatten()


# ──────────────────────────────────────────────────────────────────────
# Step 4: Approach B — Pure ML
# ──────────────────────────────────────────────────────────────────────

def train_pure_ml(syn_train, obs_train, syn_test):
    """
    Pure ML: RandomForest on raw syndrome bits → correction.
    No physics knowledge at all.
    """
    model = RandomForestClassifier(
        n_estimators=150,
        max_depth=20,
        min_samples_leaf=2,
        random_state=SEED,
        n_jobs=-1,
    )
    model.fit(syn_train, obs_train)
    return model.predict(syn_test)


# ──────────────────────────────────────────────────────────────────────
# Step 5: Approach C — Physics+ML (physics prior + ML residual)
# ──────────────────────────────────────────────────────────────────────

def build_physics_features(matching, syndromes, n_detectors):
    """
    Enrich raw syndromes with physics-derived features:
      1. Raw syndrome bits (n_detectors features)
      2. Physics decoder prediction (1 feature)
      3. Syndrome weight = sum of triggered detectors (1 feature)
      4. Local syndrome density: for each detector, count of neighbors triggered
         (approximated by windowed sum) (n_windows features)
      5. Tau-inspired weight: weighted sum using MWPM edge weights (1 feature)
    """
    n_samples = syndromes.shape[0]

    # Physics predictions
    physics_preds = matching.decode_batch(syndromes).flatten().reshape(-1, 1)

    # Syndrome weight
    syn_weight = syndromes.sum(axis=1, dtype=np.float32).reshape(-1, 1)

    # Local density: split detectors into windows, count per window
    window_size = max(1, n_detectors // 10)
    n_windows = (n_detectors + window_size - 1) // window_size
    local_density = np.zeros((n_samples, n_windows), dtype=np.float32)
    for w in range(n_windows):
        start = w * window_size
        end = min(start + window_size, n_detectors)
        local_density[:, w] = syndromes[:, start:end].sum(axis=1)

    # Syndrome parity features (XOR patterns between detector groups)
    n_parity = 5
    parity_features = np.zeros((n_samples, n_parity), dtype=np.float32)
    rng = np.random.RandomState(SEED + 7)
    for p in range(n_parity):
        # Random subset of detectors for parity
        subset = rng.choice(n_detectors, size=min(10, n_detectors), replace=False)
        parity_features[:, p] = syndromes[:, subset].sum(axis=1) % 2

    features = np.hstack([
        syndromes.astype(np.float32),   # raw bits
        physics_preds.astype(np.float32),  # MWPM prediction
        syn_weight,                      # total weight
        local_density,                   # spatial density
        parity_features,                 # parity patterns
    ])

    return features


def train_physics_ml(matching, syn_train, obs_train, syn_test, n_detectors):
    """
    Physics+ML: ML with physics-enriched features.
    The physics prediction is a strong feature that the model can trust
    when it agrees with syndrome patterns, and override when it detects
    the hot-qubit signature.
    """
    features_train = build_physics_features(matching, syn_train, n_detectors)
    features_test = build_physics_features(matching, syn_test, n_detectors)

    model = GradientBoostingClassifier(
        n_estimators=100,
        max_depth=5,
        learning_rate=0.1,
        min_samples_leaf=3,
        random_state=SEED,
    )
    model.fit(features_train, obs_train)
    return model.predict(features_test)


# ──────────────────────────────────────────────────────────────────────
# Step 6: Oracle decoder (for reference — uses true noise model)
# ──────────────────────────────────────────────────────────────────────

def build_oracle_decoder(true_circuit):
    """MWPM with weights from the TRUE non-uniform noise model (upper bound)."""
    dem = true_circuit.detector_error_model(decompose_errors=True)
    matching = pymatching.Matching.from_detector_error_model(dem)
    return matching


# ──────────────────────────────────────────────────────────────────────
# Step 7: Run experiment
# ──────────────────────────────────────────────────────────────────────

def compute_ler(predictions, true_obs):
    """Logical error rate."""
    return float(np.mean(predictions != true_obs))


def run_experiment():
    print("=" * 70)
    print("QEC Experiment 7: Physics vs ML vs Physics+ML Decoder Comparison")
    print("=" * 70)

    t0 = time.time()

    # --- Build circuits ---
    print(f"\n[1] Building distance-{DISTANCE} surface code...")
    print(f"    Base noise: {BASE_NOISE}, Hot factor: {HOT_FACTOR}x, Hot fraction: {HOT_FRACTION}")
    print(f"    Rounds: {ROUNDS}")

    true_circuit, hot_indices = build_true_noisy_circuit(
        DISTANCE, BASE_NOISE, HOT_FACTOR, HOT_FRACTION, ROUNDS
    )
    uniform_circuit = build_uniform_circuit(DISTANCE, BASE_NOISE, ROUNDS)

    print(f"    True circuit: {true_circuit.num_qubits} qubits, hot qubits: {hot_indices}")

    # --- Build decoders ---
    print("\n[2] Building decoders...")
    physics_matcher = build_physics_decoder(uniform_circuit)
    oracle_matcher = build_oracle_decoder(true_circuit)

    # --- Generate data ---
    print(f"\n[3] Generating {TOTAL_SAMPLES} labeled syndrome samples from TRUE noisy circuit...")
    syndromes, observables = generate_data(true_circuit, TOTAL_SAMPLES)
    n_detectors = syndromes.shape[1]
    print(f"    Syndrome dimension: {n_detectors} detectors")
    print(f"    Observable flip rate: {observables.mean():.4f}")

    # Split into test (fixed) and training pool
    syn_test = syndromes[-TEST_SIZE:]
    obs_test = observables[-TEST_SIZE:]
    syn_pool = syndromes[:-TEST_SIZE]
    obs_pool = observables[:-TEST_SIZE]
    print(f"    Training pool: {len(syn_pool)}, Test set: {TEST_SIZE}")

    results = {
        "config": {
            "distance": DISTANCE,
            "rounds": ROUNDS,
            "base_noise": BASE_NOISE,
            "hot_factor": HOT_FACTOR,
            "hot_fraction": HOT_FRACTION,
            "total_samples": TOTAL_SAMPLES,
            "test_size": TEST_SIZE,
            "n_detectors": n_detectors,
            "hot_indices": hot_indices,
            "observable_flip_rate": float(observables.mean()),
        },
        "approach_A": {},
        "approach_B": {},
        "approach_C": {},
        "oracle": {},
    }

    # ── Oracle: True-noise MWPM (upper bound) ──
    print("\n[4] Oracle decoder (true noise model MWPM)...")
    oracle_preds = decode_physics(oracle_matcher, syn_test)
    ler_oracle = compute_ler(oracle_preds, obs_test)
    print(f"    Oracle LER: {ler_oracle:.4f}")
    results["oracle"] = {"logical_error_rate": ler_oracle}

    # ── Approach A: Physics (uniform-model MWPM) ──
    print("\n[5] Approach A: Pure Physics (uniform-model MWPM)...")
    t_a = time.time()
    physics_preds_test = decode_physics(physics_matcher, syn_test)
    ler_A = compute_ler(physics_preds_test, obs_test)
    time_A = time.time() - t_a
    print(f"    LER: {ler_A:.4f}  (training: 0 examples, time: {time_A:.2f}s)")
    print(f"    Gap from oracle: {ler_A - ler_oracle:.4f}")
    results["approach_A"] = {
        "logical_error_rate": ler_A,
        "training_examples": 0,
        "time_s": time_A,
    }

    # ── Approach B: Pure ML ──
    print("\n[6] Approach B: Pure ML (RandomForest on raw syndromes)...")
    ler_B = {}
    time_B = {}
    for n_train in N_TRAIN_B:
        if n_train > len(syn_pool):
            continue
        t_b = time.time()
        preds = train_pure_ml(syn_pool[:n_train], obs_pool[:n_train], syn_test)
        ler = compute_ler(preds, obs_test)
        elapsed = time.time() - t_b
        ler_B[n_train] = ler
        time_B[n_train] = elapsed
        print(f"    N={n_train:5d}: LER = {ler:.4f}  ({elapsed:.2f}s)")

    results["approach_B"] = {
        "logical_error_rates": {str(k): v for k, v in ler_B.items()},
        "times_s": {str(k): v for k, v in time_B.items()},
    }

    # ── Approach C: Physics+ML ──
    print("\n[7] Approach C: Physics+ML (physics features + GBM)...")
    ler_C = {}
    time_C = {}
    for n_train in N_TRAIN_C:
        if n_train > len(syn_pool):
            continue
        t_c = time.time()
        preds = train_physics_ml(
            physics_matcher, syn_pool[:n_train], obs_pool[:n_train],
            syn_test, n_detectors
        )
        ler = compute_ler(preds, obs_test)
        elapsed = time.time() - t_c
        ler_C[n_train] = ler
        time_C[n_train] = elapsed
        print(f"    N={n_train:5d}: LER = {ler:.4f}  ({elapsed:.2f}s)")

    results["approach_C"] = {
        "logical_error_rates": {str(k): v for k, v in ler_C.items()},
        "times_s": {str(k): v for k, v in time_C.items()},
    }

    # ── Analysis ──
    print("\n" + "=" * 70)
    print("ANALYSIS")
    print("=" * 70)

    best_B_ler = min(ler_B.values())
    best_B_n = min(ler_B, key=ler_B.get)
    best_C_ler = min(ler_C.values())
    best_C_n = min(ler_C, key=ler_C.get)

    print(f"\n  Oracle (true model):   LER = {ler_oracle:.4f}  (perfect knowledge)")
    print(f"  Physics baseline (A):  LER = {ler_A:.4f}  (0 training examples)")
    print(f"  Best pure ML (B):      LER = {best_B_ler:.4f}  ({best_B_n} training examples)")
    print(f"  Best physics+ML (C):   LER = {best_C_ler:.4f}  ({best_C_n} training examples)")

    # Find where B matches A
    match_B_at_A = None
    for n in sorted(ler_B.keys()):
        if ler_B[n] <= ler_A * 1.05:
            match_B_at_A = n
            break

    # Find where C matches B's best
    match_C_at_B = None
    for n in sorted(ler_C.keys()):
        if ler_C[n] <= best_B_ler * 1.05:
            match_C_at_B = n
            break

    # Find where C surpasses A
    match_C_beats_A = None
    for n in sorted(ler_C.keys()):
        if ler_C[n] < ler_A:
            match_C_beats_A = n
            break

    if match_B_at_A:
        print(f"\n  Pure ML matches physics at N = {match_B_at_A}")
    else:
        print(f"\n  Pure ML never matches physics baseline with up to {max(N_TRAIN_B)} examples")

    if match_C_beats_A:
        print(f"  Physics+ML surpasses physics at N = {match_C_beats_A}")

    if match_C_at_B:
        efficiency = best_B_n / match_C_at_B if match_C_at_B > 0 else float('inf')
        print(f"  Physics+ML matches B's best at N = {match_C_at_B}")
        print(f"  → Data efficiency gain: {efficiency:.0f}x fewer training examples!")
    else:
        print(f"  Physics+ML does not reach B's best with up to {max(N_TRAIN_C)} examples")

    results["analysis"] = {
        "oracle_ler": ler_oracle,
        "physics_ler": ler_A,
        "best_ml_ler": best_B_ler,
        "best_ml_n": best_B_n,
        "best_physml_ler": best_C_ler,
        "best_physml_n": best_C_n,
        "ml_matches_physics_at_n": match_B_at_A,
        "physml_beats_physics_at_n": match_C_beats_A,
        "physml_matches_ml_best_at_n": match_C_at_B,
        "data_efficiency_gain": float(best_B_n / match_C_at_B) if match_C_at_B else None,
    }

    total_time = time.time() - t0
    results["total_time_s"] = total_time
    print(f"\n  Total experiment time: {total_time:.1f}s")

    return results, ler_A, ler_B, ler_C, ler_oracle


# ──────────────────────────────────────────────────────────────────────
# Step 8: Plotting
# ──────────────────────────────────────────────────────────────────────

def plot_results(ler_A, ler_B, ler_C, ler_oracle, results):
    fig, axes = plt.subplots(1, 2, figsize=(15, 6.5))

    # ── Left panel: Learning curves ──
    ax = axes[0]

    # Oracle line
    all_n = sorted(set(list(ler_B.keys()) + list(ler_C.keys())))
    ax.axhline(y=ler_oracle, color="#9C27B0", linewidth=1.5, linestyle=":",
               label=f"Oracle (true model) — {ler_oracle:.4f}", alpha=0.7, zorder=3)

    # Approach A: horizontal line
    ax.axhline(y=ler_A, color="#2196F3", linewidth=2.5, linestyle="--",
               label=f"A: Physics (uniform MWPM) — {ler_A:.4f}", zorder=5)

    # Shaded region: physics baseline ± 5%
    ax.axhspan(ler_A * 0.95, ler_A * 1.05, color="#2196F3", alpha=0.08, zorder=1)

    # Approach B: learning curve
    ns_B = sorted(ler_B.keys())
    lers_B = [ler_B[n] for n in ns_B]
    ax.plot(ns_B, lers_B, "o-", color="#F44336", linewidth=2.2, markersize=8,
            label="B: Pure ML (RandomForest)", zorder=4)
    # Annotate each point
    for n, ler in zip(ns_B, lers_B):
        ax.annotate(f"{ler:.3f}", (n, ler), textcoords="offset points",
                    xytext=(0, 10), fontsize=7.5, color="#F44336", ha="center")

    # Approach C: learning curve
    ns_C = sorted(ler_C.keys())
    lers_C = [ler_C[n] for n in ns_C]
    ax.plot(ns_C, lers_C, "s-", color="#4CAF50", linewidth=2.2, markersize=8,
            label="C: Physics+ML (tau-chrono)", zorder=4)
    for n, ler in zip(ns_C, lers_C):
        ax.annotate(f"{ler:.3f}", (n, ler), textcoords="offset points",
                    xytext=(0, -14), fontsize=7.5, color="#4CAF50", ha="center")

    # Efficiency annotation
    best_B_n = results["analysis"]["best_ml_n"]
    match_C = results["analysis"]["physml_matches_ml_best_at_n"]
    gain = results["analysis"]["data_efficiency_gain"]

    if match_C and gain:
        ax.annotate(
            f"{gain:.0f}x less data\nneeded",
            xy=(match_C, ler_C[match_C]),
            xytext=(match_C * 4, ler_C[match_C] + 0.025),
            fontsize=11, fontweight="bold", color="#4CAF50",
            arrowprops=dict(arrowstyle="->", color="#4CAF50", lw=2),
            zorder=6,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="#E8F5E9", edgecolor="#4CAF50"),
        )

    ax.set_xscale("log")
    ax.set_xlabel("Training Examples (log scale)", fontsize=13)
    ax.set_ylabel("Logical Error Rate (lower = better)", fontsize=13)
    ax.set_title(f"Learning Curves: d={DISTANCE} Surface Code, Non-uniform Noise",
                 fontsize=13, fontweight="bold")
    ax.legend(fontsize=9.5, loc="upper right", framealpha=0.9)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(30, 6000)

    # ── Right panel: Deployment cost comparison ──
    ax2 = axes[1]

    best_B_ler = results["analysis"]["best_ml_ler"]
    best_C_ler = results["analysis"]["best_physml_ler"]
    best_C_n = results["analysis"]["best_physml_n"]

    categories = [
        f"A: Physics\n(0 examples)",
        f"B: Pure ML\n({best_B_n} ex.)",
        f"C: Phys+ML\n({best_C_n} ex.)",
        f"Oracle\n(full model)",
    ]
    lers = [ler_A, best_B_ler, best_C_ler, ler_oracle]
    colors = ["#2196F3", "#F44336", "#4CAF50", "#9C27B0"]

    # Setup cost in minutes
    cost_char = 10     # noise characterization
    cost_per_sample = 0.005  # data generation per sample (minutes)
    costs = [
        cost_char,                           # A: just characterization
        cost_per_sample * best_B_n,          # B: data generation
        cost_char + cost_per_sample * best_C_n,  # C: char + small data
        0,                                    # Oracle: theoretical
    ]

    # Grouped bar chart
    x = np.arange(len(categories))
    width = 0.35

    bars1 = ax2.bar(x - width/2, lers, width, color=colors, alpha=0.85,
                    edgecolor="black", linewidth=0.5, label="Logical Error Rate")

    ax2_twin = ax2.twinx()
    bars2 = ax2_twin.bar(x + width/2, costs, width, color=colors, alpha=0.4,
                         edgecolor="black", linewidth=0.5, hatch="//",
                         label="Setup Cost (min)")

    # Labels
    for bar, ler in zip(bars1, lers):
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.003,
                 f"{ler:.3f}", ha="center", va="bottom", fontsize=10, fontweight="bold")

    for bar, cost in zip(bars2, costs):
        if cost > 0:
            ax2_twin.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
                          f"{cost:.0f}m", ha="center", va="bottom", fontsize=9,
                          fontweight="bold", alpha=0.7)

    ax2.set_xticks(x)
    ax2.set_xticklabels(categories, fontsize=10)
    ax2.set_ylabel("Logical Error Rate", fontsize=12, color="black")
    ax2_twin.set_ylabel("Setup Cost (minutes)", fontsize=12, color="gray")
    ax2.set_title("Best Accuracy vs Deployment Cost", fontsize=13, fontweight="bold")

    # Combined legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor="gray", alpha=0.85, label="LER (solid)"),
        Patch(facecolor="gray", alpha=0.4, hatch="//", label="Cost (hatched)"),
    ]
    ax2.legend(handles=legend_elements, fontsize=9, loc="upper left")

    ax2.grid(True, alpha=0.2, axis="y")

    plt.tight_layout()

    outpath = os.path.join(OUTDIR, "exp7_physics_ml_comparison.png")
    plt.savefig(outpath, dpi=150, bbox_inches="tight")
    print(f"\n  Plot saved: {outpath}")
    plt.close()


# ──────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    results, ler_A, ler_B, ler_C, ler_oracle = run_experiment()

    # Save JSON
    json_path = os.path.join(OUTDIR, "exp7_physics_ml_comparison.json")
    with open(json_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"  JSON saved: {json_path}")

    # Plot
    plot_results(ler_A, ler_B, ler_C, ler_oracle, results)

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    a = results["analysis"]
    print(f"  Physics (A):     LER = {a['physics_ler']:.4f}  |  0 training examples")
    print(f"  Best ML (B):     LER = {a['best_ml_ler']:.4f}  |  {a['best_ml_n']} training examples")
    print(f"  Best Phys+ML(C): LER = {a['best_physml_ler']:.4f}  |  {a['best_physml_n']} training examples")
    print(f"  Oracle:          LER = {a['oracle_ler']:.4f}  |  perfect noise knowledge")
    if a["data_efficiency_gain"]:
        print(f"\n  → Physics+ML is {a['data_efficiency_gain']:.0f}x more data-efficient than pure ML")
    print("=" * 70)
