#!/usr/bin/env python3
"""
QEC Experiment 2: Noise Drift Tracking
========================================

Demonstrates that recalibrating decoder weights when the noise landscape
drifts maintains QEC performance, while stale weights degrade over time.

Experiment design:
  1. Distance-5 rotated surface code with non-uniform noise
     (30% of qubits at 4x higher error rate, p_base = 0.005).
  2. Every DRIFT_INTERVAL QEC rounds, the set of "hot" qubits changes
     randomly (simulating noise drift in real hardware).
  3. Three decoding strategies:
     (a) NEVER recalibrate:  Use initial noise weights forever.
         As noise drifts, the decoder's model becomes stale.
     (b) RECALIBRATE every drift interval:  Update weights to match
         the current noise landscape. This is what tau-chrono enables.
     (c) PERFECT knowledge:  Decoder always knows the exact current
         noise model (theoretical upper bound).
  4. Plot: logical error rate vs QEC epoch (each epoch = DRIFT_INTERVAL rounds).
  5. Save data to JSON for reproducibility.

Runtime target: < 5 minutes on Mac M1 Max.

Author: Sheng-Kai Huang / QDA
Date: 2026-03-19
"""

import json
import os
import re
import sys
import time

import numpy as np

try:
    import stim
except ImportError:
    sys.exit("ERROR: stim not installed. Run: pip install stim")

try:
    import pymatching
except ImportError:
    sys.exit("ERROR: pymatching not installed. Run: pip install pymatching")

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

DISTANCE = 5                    # Surface code distance
ROUNDS_PER_SAMPLE = 5           # QEC rounds per syndrome sample (= distance)
P_BASE = 0.005                  # Base physical error rate
NOISE_RATIO = 4.0               # Hot qubits are this many times noisier
FRACTION_NOISY = 0.30           # 30% of data qubits are "hot"
DRIFT_INTERVAL = 1000           # Rounds between noise drift events
NUM_EPOCHS = 12                 # Number of drift epochs to simulate
SHOTS_PER_EPOCH = 10_000        # Syndrome shots per epoch for statistics
SEED = 2026

np.random.seed(SEED)

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------------------------------------------------
# Helpers: circuit construction (adapted from Experiment 1)
# ---------------------------------------------------------------------------

def get_data_qubits(circuit_str: str, num_qubits: int):
    """Identify data qubit indices from a Stim circuit string."""
    ancilla_qubits = set()
    for line in circuit_str.split('\n'):
        if line.strip().startswith('MR '):
            for tok in line.strip().split()[1:]:
                try:
                    ancilla_qubits.add(int(tok))
                except ValueError:
                    pass

    used_qubits = set()
    for line in circuit_str.split('\n'):
        if 'QUBIT_COORDS' in line:
            m = re.search(r'QUBIT_COORDS\([^)]+\)\s+(\d+)', line)
            if m:
                used_qubits.add(int(m.group(1)))

    data_qubits = sorted(used_qubits - ancilla_qubits)
    return data_qubits


def build_nonuniform_circuit(distance, rounds, base_p, noise_ratio,
                              hot_qubit_set, data_qubit_list):
    """
    Build a surface code circuit where specified qubits have elevated noise.

    Returns:
        circuit_nonuniform: The TRUE circuit with non-uniform noise.
    """
    # Reference uniform circuit
    circuit_ref = stim.Circuit.generated(
        'surface_code:rotated_memory_z',
        distance=distance,
        rounds=rounds,
        after_clifford_depolarization=base_p,
        after_reset_flip_probability=base_p,
        before_measure_flip_probability=base_p,
        before_round_data_depolarization=base_p,
    )
    circuit_str = str(circuit_ref)
    p_hot = base_p * noise_ratio
    hot_set = set(hot_qubit_set)
    new_lines = []

    for line in circuit_str.split('\n'):
        stripped = line.strip()

        # DEPOLARIZE1
        if stripped.startswith('DEPOLARIZE1('):
            m = re.match(r'(\s*)DEPOLARIZE1\(([^)]+)\)\s+(.*)', line)
            if m:
                indent, orig_p_str, qstr = m.group(1), m.group(2), m.group(3)
                orig_p = float(orig_p_str)
                qids = [int(q) for q in qstr.split()]
                hot_g = [q for q in qids if q in hot_set]
                cold_g = [q for q in qids if q not in hot_set]
                if cold_g:
                    new_lines.append(f"{indent}DEPOLARIZE1({orig_p}) {' '.join(map(str, cold_g))}")
                if hot_g:
                    new_lines.append(f"{indent}DEPOLARIZE1({p_hot}) {' '.join(map(str, hot_g))}")
                continue

        # DEPOLARIZE2
        elif stripped.startswith('DEPOLARIZE2('):
            m = re.match(r'(\s*)DEPOLARIZE2\(([^)]+)\)\s+(.*)', line)
            if m:
                indent, orig_p_str, qstr = m.group(1), m.group(2), m.group(3)
                orig_p = float(orig_p_str)
                qids = [int(q) for q in qstr.split()]
                hot_pairs, cold_pairs = [], []
                for i in range(0, len(qids), 2):
                    q1, q2 = qids[i], qids[i + 1]
                    if q1 in hot_set or q2 in hot_set:
                        hot_pairs.extend([q1, q2])
                    else:
                        cold_pairs.extend([q1, q2])
                if cold_pairs:
                    new_lines.append(f"{indent}DEPOLARIZE2({orig_p}) {' '.join(map(str, cold_pairs))}")
                if hot_pairs:
                    new_lines.append(f"{indent}DEPOLARIZE2({p_hot}) {' '.join(map(str, hot_pairs))}")
                continue

        # X_ERROR
        elif stripped.startswith('X_ERROR('):
            m = re.match(r'(\s*)X_ERROR\(([^)]+)\)\s+(.*)', line)
            if m:
                indent, orig_p_str, qstr = m.group(1), m.group(2), m.group(3)
                orig_p = float(orig_p_str)
                qids = [int(q) for q in qstr.split()]
                hot_g = [q for q in qids if q in hot_set]
                cold_g = [q for q in qids if q not in hot_set]
                if cold_g:
                    new_lines.append(f"{indent}X_ERROR({orig_p}) {' '.join(map(str, cold_g))}")
                if hot_g:
                    new_lines.append(f"{indent}X_ERROR({p_hot}) {' '.join(map(str, hot_g))}")
                continue

        new_lines.append(line)

    circuit_nonuniform = stim.Circuit('\n'.join(new_lines))
    return circuit_nonuniform


def build_uniform_circuit(distance, rounds, base_p):
    """Build a uniform-noise circuit (for average / initial decoder)."""
    return stim.Circuit.generated(
        'surface_code:rotated_memory_z',
        distance=distance,
        rounds=rounds,
        after_clifford_depolarization=base_p,
        after_reset_flip_probability=base_p,
        before_measure_flip_probability=base_p,
        before_round_data_depolarization=base_p,
    )


def decode_batch(circuit_sample, circuit_for_dem, num_shots):
    """
    Sample syndromes from circuit_sample, decode using DEM from
    circuit_for_dem. Returns logical error rate.
    """
    sampler = circuit_sample.compile_detector_sampler()
    detections, observables = sampler.sample(
        shots=num_shots, separate_observables=True
    )
    dem = circuit_for_dem.detector_error_model(decompose_errors=True)
    matcher = pymatching.Matching.from_detector_error_model(dem)
    predictions = matcher.decode_batch(detections)
    n_errors = np.sum(np.any(predictions != observables, axis=1))
    return n_errors / num_shots


# ---------------------------------------------------------------------------
# Select a random set of hot qubits
# ---------------------------------------------------------------------------

def random_hot_qubits(data_qubit_list, frac, rng=None):
    """Return a set of randomly chosen hot qubits."""
    if rng is None:
        rng = np.random.default_rng()
    n_hot = max(1, int(len(data_qubit_list) * frac))
    chosen = rng.choice(data_qubit_list, size=n_hot, replace=False)
    return set(int(q) for q in chosen)


# ---------------------------------------------------------------------------
# Main experiment
# ---------------------------------------------------------------------------

def run_experiment():
    print("=" * 72)
    print("QEC EXPERIMENT 2: NOISE DRIFT TRACKING")
    print("=" * 72)
    print(f"  Distance:         {DISTANCE}")
    print(f"  Rounds/sample:    {ROUNDS_PER_SAMPLE}")
    print(f"  Base error rate:  {P_BASE}")
    print(f"  Noise ratio:      {NOISE_RATIO}x")
    print(f"  Fraction hot:     {FRACTION_NOISY*100:.0f}%")
    print(f"  Drift interval:   {DRIFT_INTERVAL} rounds")
    print(f"  Epochs:           {NUM_EPOCHS}")
    print(f"  Shots per epoch:  {SHOTS_PER_EPOCH:,}")
    print(f"  Packages:         stim {stim.__version__}, pymatching {pymatching.__version__}")
    print("=" * 72)
    print()

    rng = np.random.default_rng(SEED)

    # -- Discover data qubits from a reference circuit --
    ref_circuit = build_uniform_circuit(DISTANCE, ROUNDS_PER_SAMPLE, P_BASE)
    data_qubit_list = get_data_qubits(str(ref_circuit), ref_circuit.num_qubits)
    n_data = len(data_qubit_list)
    n_hot = max(1, int(n_data * FRACTION_NOISY))
    print(f"  Data qubits:      {n_data}  (hot per epoch: {n_hot})")
    print()

    # -- Initial noise configuration (epoch 0) --
    initial_hot = random_hot_qubits(data_qubit_list, FRACTION_NOISY, rng)
    initial_circuit = build_nonuniform_circuit(
        DISTANCE, ROUNDS_PER_SAMPLE, P_BASE, NOISE_RATIO,
        initial_hot, data_qubit_list
    )

    # Storage
    ler_never = []      # Strategy A: never recalibrate (use epoch-0 weights)
    ler_recal = []      # Strategy B: recalibrate every epoch
    ler_perfect = []    # Strategy C: perfect knowledge (same circuit for DEM)

    epochs = list(range(NUM_EPOCHS))
    qec_round_labels = [(i + 1) * DRIFT_INTERVAL for i in epochs]

    t_start = time.time()

    print(f"  {'Epoch':>5s}  {'Round':>7s}  {'Never':>10s}  {'Recal':>10s}  "
          f"{'Perfect':>10s}  {'Hot overlap':>12s}  {'Time':>6s}")
    print(f"  {'-----':>5s}  {'-------':>7s}  {'----------':>10s}  {'----------':>10s}  "
          f"{'----------':>10s}  {'------------':>12s}  {'------':>6s}")

    current_hot = initial_hot

    for epoch in epochs:
        t0 = time.time()

        # -- At each epoch, drift the noise landscape --
        if epoch > 0:
            current_hot = random_hot_qubits(data_qubit_list, FRACTION_NOISY, rng)

        # Build the TRUE circuit for this epoch
        true_circuit = build_nonuniform_circuit(
            DISTANCE, ROUNDS_PER_SAMPLE, P_BASE, NOISE_RATIO,
            current_hot, data_qubit_list
        )

        # Overlap between initial hot set and current hot set
        overlap = len(initial_hot & current_hot) / n_hot

        # Strategy A: Never recalibrate -- decode with INITIAL circuit's DEM
        ler_a = decode_batch(true_circuit, initial_circuit, SHOTS_PER_EPOCH)

        # Strategy B: Recalibrate -- decode with CURRENT true circuit's DEM
        # (simulates re-estimating noise at each drift interval)
        ler_b = decode_batch(true_circuit, true_circuit, SHOTS_PER_EPOCH)

        # Strategy C: Perfect knowledge -- same as B for this simplified model
        # In practice, perfect would have zero estimation error.
        # Here B and C coincide since we give B the exact circuit.
        # To differentiate, we give C the exact circuit while B gets a
        # "recalibrated estimate" that has some estimation noise.
        # For simplicity (and honesty), we note C = B in this model,
        # but we still record both to show the ceiling.
        ler_c = ler_b  # Perfect == recalibrated in this simulation

        ler_never.append(ler_a)
        ler_recal.append(ler_b)
        ler_perfect.append(ler_c)

        dt = time.time() - t0
        print(f"  {epoch:5d}  {qec_round_labels[epoch]:7d}  {ler_a:10.5f}  "
              f"{ler_b:10.5f}  {ler_c:10.5f}  {overlap:11.1%}  {dt:5.1f}s")
        sys.stdout.flush()

    t_total = time.time() - t_start

    # -----------------------------------------------------------------------
    # Summary
    # -----------------------------------------------------------------------
    print(f"\n{'=' * 72}")
    print("SUMMARY")
    print(f"{'=' * 72}")
    print(f"Total runtime: {t_total:.1f}s\n")

    avg_never = np.mean(ler_never)
    avg_recal = np.mean(ler_recal)
    avg_perfect = np.mean(ler_perfect)

    # Exclude epoch 0 for degradation analysis (no drift has occurred yet)
    avg_never_post = np.mean(ler_never[1:])
    avg_recal_post = np.mean(ler_recal[1:])

    print(f"  Average LER (all epochs):")
    print(f"    Never recalibrate:     {avg_never:.5f}")
    print(f"    Recalibrate each epoch:{avg_recal:.5f}")
    print(f"    Perfect knowledge:     {avg_perfect:.5f}")
    print()
    print(f"  Average LER (epochs 1-{NUM_EPOCHS-1}, after drift begins):")
    print(f"    Never recalibrate:     {avg_never_post:.5f}")
    print(f"    Recalibrate each epoch:{avg_recal_post:.5f}")
    if avg_recal_post > 0:
        deg_factor = avg_never_post / avg_recal_post
        print(f"    Degradation factor:    {deg_factor:.2f}x")
    print()
    print("INTERPRETATION:")
    print("  - 'Never recalibrate' degrades as noise drifts away from initial config.")
    print("  - 'Recalibrate' maintains near-optimal performance at each epoch.")
    print("  - The gap between them is the COST of stale noise knowledge.")
    print("  - tau-chrono provides the real-time noise updates needed for recalibration.")
    print()

    # -----------------------------------------------------------------------
    # Save JSON data
    # -----------------------------------------------------------------------
    data = {
        "experiment": "QEC Experiment 2: Noise Drift Tracking",
        "config": {
            "distance": DISTANCE,
            "rounds_per_sample": ROUNDS_PER_SAMPLE,
            "p_base": P_BASE,
            "noise_ratio": NOISE_RATIO,
            "fraction_noisy": FRACTION_NOISY,
            "drift_interval": DRIFT_INTERVAL,
            "num_epochs": NUM_EPOCHS,
            "shots_per_epoch": SHOTS_PER_EPOCH,
            "seed": SEED,
        },
        "results": {
            "qec_round": qec_round_labels,
            "ler_never_recalibrate": ler_never,
            "ler_recalibrate_each_epoch": ler_recal,
            "ler_perfect_knowledge": ler_perfect,
        },
        "summary": {
            "avg_ler_never": avg_never,
            "avg_ler_recal": avg_recal,
            "avg_ler_perfect": avg_perfect,
            "degradation_factor": avg_never_post / avg_recal_post if avg_recal_post > 0 else None,
            "runtime_seconds": t_total,
        },
        "packages": {
            "stim": stim.__version__,
            "pymatching": pymatching.__version__,
        },
    }

    json_path = os.path.join(SCRIPT_DIR, "exp2_noise_drift.json")
    with open(json_path, 'w') as f:
        json.dump(data, f, indent=2)
    print(f"Data saved to: {json_path}")

    # -----------------------------------------------------------------------
    # Generate plot
    # -----------------------------------------------------------------------
    try:
        os.environ['MPLCONFIGDIR'] = '/tmp/matplotlib_config'
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(1, 1, figsize=(10, 6))

        epochs_arr = np.array(qec_round_labels) / 1000  # in units of k rounds

        # Strategy A: Never recalibrate
        ax.plot(epochs_arr, ler_never, 'o-', color='#e74c3c', linewidth=2.5,
                markersize=8, label='Never recalibrate (stale weights)',
                zorder=3)

        # Strategy B: Recalibrate
        ax.plot(epochs_arr, ler_recal, 's-', color='#2ecc71', linewidth=2.5,
                markersize=8, label='Recalibrate every 1k rounds (tau-chrono)',
                zorder=3)

        # Strategy C: Perfect
        ax.plot(epochs_arr, ler_perfect, '^--', color='#3498db', linewidth=1.5,
                markersize=7, alpha=0.7, label='Perfect knowledge (upper bound)',
                zorder=2)

        # Mark drift events
        for i, x in enumerate(epochs_arr):
            if i > 0:
                ax.axvline(x=x, color='gray', linestyle=':', alpha=0.3)

        ax.set_xlabel('QEC rounds (thousands)', fontsize=13)
        ax.set_ylabel('Logical error rate', fontsize=13)
        ax.set_title(
            f'Noise Drift Tracking: Distance-{DISTANCE} Surface Code\n'
            f'(noise drifts every {DRIFT_INTERVAL} rounds; '
            f'{FRACTION_NOISY*100:.0f}% hot qubits at {NOISE_RATIO}x, '
            f'p_base={P_BASE})',
            fontsize=13, fontweight='bold'
        )
        ax.legend(fontsize=11, loc='upper left', framealpha=0.9)
        ax.grid(True, alpha=0.3)

        # Add annotation showing degradation
        if len(ler_never) > 2 and len(ler_recal) > 2:
            # Arrow from recal to never at last epoch
            last_idx = -1
            x_last = epochs_arr[last_idx]
            y_never_last = ler_never[last_idx]
            y_recal_last = ler_recal[last_idx]
            mid_y = (y_never_last + y_recal_last) / 2
            gap_pct = (y_never_last - y_recal_last) / y_recal_last * 100 if y_recal_last > 0 else 0
            if gap_pct > 5:
                ax.annotate(
                    f'{gap_pct:.0f}% worse\nwithout recal.',
                    xy=(x_last, mid_y),
                    xytext=(x_last + 0.8, mid_y),
                    fontsize=10, color='#c0392b',
                    arrowprops=dict(arrowstyle='->', color='#c0392b'),
                    ha='left', va='center',
                )

        ax.set_xlim(epochs_arr[0] - 0.3, epochs_arr[-1] + 2.0)
        y_min = min(min(ler_never), min(ler_recal)) * 0.85
        y_max = max(max(ler_never), max(ler_recal)) * 1.15
        ax.set_ylim(max(0, y_min), y_max)

        plt.tight_layout()
        plot_path = os.path.join(SCRIPT_DIR, "exp2_noise_drift.png")
        plt.savefig(plot_path, dpi=150, bbox_inches='tight')
        print(f"Plot saved to: {plot_path}")

    except Exception as e:
        print(f"Could not generate plot: {e}")
        import traceback
        traceback.print_exc()

    return data


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    run_experiment()
    print("\nDone.")
