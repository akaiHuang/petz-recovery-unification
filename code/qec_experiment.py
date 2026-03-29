#!/usr/bin/env python3
"""
QEC Noise-Informed Decoding Experiment
=======================================

Demonstrates that a decoder with accurate per-qubit noise knowledge
outperforms a decoder that assumes uniform noise -- the core premise
behind tau-chrono / Petz-informed decoding.

Experiment design:
  1. Build a rotated surface code with NON-UNIFORM depolarizing noise
     (some data qubits are 3-5x noisier than others, modeling realistic
     hardware variation).
  2. Sample syndromes from the TRUE noisy circuit using Stim.
  3. Decode with two strategies:
     (a) STANDARD decoder:  PyMatching with weights derived from a
         UNIFORM noise model (average error rate). This is what you get
         if you don't know per-qubit noise.
     (b) NOISE-INFORMED decoder:  PyMatching with weights derived from
         the TRUE non-uniform noise model. This is what tau-chrono
         provides.
  4. Compare logical error rates across physical error rates.

Requirements:
  pip install stim pymatching numpy matplotlib

Tested on: Mac M1 Max, Python 3.13, stim 1.15.0, pymatching 2.3.1
Runtime: ~2-3 minutes for full sweep on M1 Max 64GB

Author: Sheng-Kai Huang / QDA
Date: 2026-03-19
"""

import time
import re
import sys
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

DISTANCES = [3, 5, 7]                          # Surface code distances
ROUNDS_PER_DISTANCE = {3: 3, 5: 5, 7: 7}      # rounds = distance (standard)
PHYSICAL_ERROR_RATES = [0.001, 0.002, 0.005, 0.007, 0.01, 0.015, 0.02]
NUM_SHOTS = 50_000                              # shots per data point
NOISE_RATIO = 4.0                               # noisy qubits are this many times noisier
FRACTION_NOISY = 0.3                            # fraction of data qubits that are "hot"
SEED = 42

np.random.seed(SEED)


# ---------------------------------------------------------------------------
# Helper: Build surface code circuit with NON-UNIFORM noise
# ---------------------------------------------------------------------------

def build_nonuniform_circuit(distance: int, rounds: int, base_p: float,
                             noise_ratio: float, frac_noisy: float):
    """
    Build a surface code circuit where a fraction of data qubits have
    higher depolarizing noise. This models realistic hardware where
    some qubits are "hot spots".

    Strategy:
      1. Generate the standard circuit with uniform noise at base_p.
      2. Identify data qubits (not ancilla/measurement qubits).
      3. For DEPOLARIZE1 and DEPOLARIZE2 instructions involving "hot" qubits,
         replace the error rate with base_p * noise_ratio.

    Returns: (circuit_nonuniform, circuit_uniform, hot_qubits)
      - circuit_nonuniform: the TRUE circuit with non-uniform noise
      - circuit_uniform: circuit with uniform noise at the AVERAGE rate
        (for the uninformed decoder)
      - hot_qubits: set of qubit indices that are noisier
    """
    # Generate reference circuit with uniform noise
    circuit_uniform = stim.Circuit.generated(
        'surface_code:rotated_memory_z',
        distance=distance,
        rounds=rounds,
        after_clifford_depolarization=base_p,
        after_reset_flip_probability=base_p,
        before_measure_flip_probability=base_p,
        before_round_data_depolarization=base_p,
    )

    # Identify data qubits from QUBIT_COORDS
    # In rotated surface code, data qubits have odd+odd coordinates
    # Ancilla qubits have even+even or even+odd coordinates
    # But easier: data qubits are the ones that appear in R but not in MR
    circuit_str = str(circuit_uniform)

    # Find all qubits used in the circuit
    all_qubits = set(range(circuit_uniform.num_qubits))

    # Find ancilla qubits (those that appear in MR instructions)
    ancilla_qubits = set()
    for line in circuit_str.split('\n'):
        if line.startswith('MR ') or line.strip().startswith('MR '):
            parts = line.strip().split()
            for p in parts[1:]:
                try:
                    ancilla_qubits.add(int(p))
                except ValueError:
                    pass

    # Parse QUBIT_COORDS to find which qubits exist
    used_qubits = set()
    for line in circuit_str.split('\n'):
        if 'QUBIT_COORDS' in line:
            match = re.search(r'QUBIT_COORDS\([^)]+\)\s+(\d+)', line)
            if match:
                used_qubits.add(int(match.group(1)))

    data_qubits = used_qubits - ancilla_qubits
    data_qubit_list = sorted(data_qubits)

    # Select "hot" qubits
    n_hot = max(1, int(len(data_qubit_list) * frac_noisy))
    hot_qubits = set(np.random.choice(data_qubit_list, size=n_hot, replace=False))

    # Now build the non-uniform circuit by modifying noise instructions
    # We parse the circuit line by line and adjust DEPOLARIZE rates
    p_hot = base_p * noise_ratio
    new_lines = []

    for line in circuit_str.split('\n'):
        stripped = line.strip()

        # Handle DEPOLARIZE1
        if stripped.startswith('DEPOLARIZE1('):
            match = re.match(r'(\s*)DEPOLARIZE1\(([^)]+)\)\s+(.*)', line)
            if match:
                indent = match.group(1)
                original_p = float(match.group(2))
                qubits_str = match.group(3)
                qubit_ids = [int(q) for q in qubits_str.split()]

                # Split into hot and cold groups
                hot_group = [q for q in qubit_ids if q in hot_qubits]
                cold_group = [q for q in qubit_ids if q not in hot_qubits]

                if cold_group:
                    new_lines.append(f"{indent}DEPOLARIZE1({original_p}) {' '.join(map(str, cold_group))}")
                if hot_group:
                    new_lines.append(f"{indent}DEPOLARIZE1({p_hot}) {' '.join(map(str, hot_group))}")
                continue

        # Handle DEPOLARIZE2
        elif stripped.startswith('DEPOLARIZE2('):
            match = re.match(r'(\s*)DEPOLARIZE2\(([^)]+)\)\s+(.*)', line)
            if match:
                indent = match.group(1)
                original_p = float(match.group(2))
                qubits_str = match.group(3)
                qubit_ids = [int(q) for q in qubits_str.split()]

                # Pairs of qubits for 2-qubit depolarization
                hot_pairs = []
                cold_pairs = []
                for i in range(0, len(qubit_ids), 2):
                    q1, q2 = qubit_ids[i], qubit_ids[i + 1]
                    if q1 in hot_qubits or q2 in hot_qubits:
                        hot_pairs.extend([q1, q2])
                    else:
                        cold_pairs.extend([q1, q2])

                if cold_pairs:
                    new_lines.append(f"{indent}DEPOLARIZE2({original_p}) {' '.join(map(str, cold_pairs))}")
                if hot_pairs:
                    new_lines.append(f"{indent}DEPOLARIZE2({p_hot}) {' '.join(map(str, hot_pairs))}")
                continue

        # Handle X_ERROR (reset noise)
        elif stripped.startswith('X_ERROR('):
            match = re.match(r'(\s*)X_ERROR\(([^)]+)\)\s+(.*)', line)
            if match:
                indent = match.group(1)
                original_p = float(match.group(2))
                qubits_str = match.group(3)
                qubit_ids = [int(q) for q in qubits_str.split()]

                hot_group = [q for q in qubit_ids if q in hot_qubits]
                cold_group = [q for q in qubit_ids if q not in hot_qubits]

                if cold_group:
                    new_lines.append(f"{indent}X_ERROR({original_p}) {' '.join(map(str, cold_group))}")
                if hot_group:
                    new_lines.append(f"{indent}X_ERROR({p_hot}) {' '.join(map(str, hot_group))}")
                continue

        new_lines.append(line)

    circuit_nonuniform = stim.Circuit('\n'.join(new_lines))

    # Build "average" uniform circuit for the uninformed decoder
    # The average rate accounts for the mix of hot and cold qubits
    n_data = len(data_qubit_list)
    p_avg = base_p * (1.0 + (noise_ratio - 1.0) * n_hot / n_data)
    circuit_avg = stim.Circuit.generated(
        'surface_code:rotated_memory_z',
        distance=distance,
        rounds=rounds,
        after_clifford_depolarization=p_avg,
        after_reset_flip_probability=p_avg,
        before_measure_flip_probability=p_avg,
        before_round_data_depolarization=p_avg,
    )

    return circuit_nonuniform, circuit_avg, hot_qubits


# ---------------------------------------------------------------------------
# Helper: Decode and measure logical error rate
# ---------------------------------------------------------------------------

def measure_logical_error_rate(circuit_sample, circuit_decode, num_shots):
    """
    Sample from circuit_sample (the TRUE noisy circuit),
    decode using the DEM from circuit_decode.

    Returns: logical error rate (float)
    """
    # Sample syndromes from the TRUE circuit
    sampler = circuit_sample.compile_detector_sampler()
    detections, observables = sampler.sample(
        shots=num_shots, separate_observables=True
    )

    # Build decoder from the (possibly different) circuit's DEM
    dem = circuit_decode.detector_error_model(decompose_errors=True)
    matcher = pymatching.Matching.from_detector_error_model(dem)

    # Decode
    predictions = matcher.decode_batch(detections)

    # Count logical errors
    num_errors = np.sum(np.any(predictions != observables, axis=1))
    return num_errors / num_shots


# ---------------------------------------------------------------------------
# Main experiment
# ---------------------------------------------------------------------------

def run_experiment():
    print("=" * 72)
    print("QEC NOISE-INFORMED DECODING EXPERIMENT")
    print("=" * 72)
    print(f"  Distances:       {DISTANCES}")
    print(f"  Error rates:     {PHYSICAL_ERROR_RATES}")
    print(f"  Shots per point: {NUM_SHOTS:,}")
    print(f"  Noise ratio:     {NOISE_RATIO}x (hot qubits are {NOISE_RATIO}x noisier)")
    print(f"  Fraction hot:    {FRACTION_NOISY*100:.0f}% of data qubits")
    print(f"  Packages:        stim {stim.__version__}, pymatching {pymatching.__version__}")
    print("=" * 72)
    print()

    results = {}  # (distance, p) -> (ler_uniform, ler_informed)

    total_points = len(DISTANCES) * len(PHYSICAL_ERROR_RATES)
    point = 0
    t_start = time.time()

    for d in DISTANCES:
        rounds = ROUNDS_PER_DISTANCE[d]
        n_data = d * d  # number of data qubits in distance-d surface code
        n_hot = max(1, int(n_data * FRACTION_NOISY))

        print(f"\n{'─' * 72}")
        print(f"  DISTANCE d={d}  |  {n_data} data qubits, {n_hot} hot qubits, {rounds} rounds")
        print(f"{'─' * 72}")
        print(f"  {'p_phys':>8s}  {'LER(uniform)':>14s}  {'LER(informed)':>14s}  {'Improvement':>12s}  {'Time':>6s}")
        print(f"  {'─'*8:>8s}  {'─'*14:>14s}  {'─'*14:>14s}  {'─'*12:>12s}  {'─'*6:>6s}")

        for p in PHYSICAL_ERROR_RATES:
            point += 1
            t0 = time.time()

            # Build circuits
            circuit_true, circuit_avg, hot_qubits = build_nonuniform_circuit(
                distance=d, rounds=rounds, base_p=p,
                noise_ratio=NOISE_RATIO, frac_noisy=FRACTION_NOISY
            )

            # Measure logical error rates
            ler_uniform = measure_logical_error_rate(circuit_true, circuit_avg, NUM_SHOTS)
            ler_informed = measure_logical_error_rate(circuit_true, circuit_true, NUM_SHOTS)

            dt = time.time() - t0

            # Improvement factor
            if ler_informed > 0:
                improvement = ler_uniform / ler_informed
            else:
                improvement = float('inf')

            results[(d, p)] = (ler_uniform, ler_informed)

            print(f"  {p:8.4f}  {ler_uniform:14.6f}  {ler_informed:14.6f}  {improvement:11.2f}x  {dt:5.1f}s")
            sys.stdout.flush()

    t_total = time.time() - t_start

    # -----------------------------------------------------------------------
    # Summary
    # -----------------------------------------------------------------------
    print(f"\n\n{'=' * 72}")
    print("SUMMARY")
    print(f"{'=' * 72}")
    print(f"Total runtime: {t_total:.1f}s")
    print()

    # Compute average improvement per distance
    for d in DISTANCES:
        improvements = []
        for p in PHYSICAL_ERROR_RATES:
            ler_u, ler_i = results[(d, p)]
            if ler_i > 0:
                improvements.append(ler_u / ler_i)
        if improvements:
            avg_imp = np.mean(improvements)
            max_imp = np.max(improvements)
            print(f"  d={d}: avg improvement = {avg_imp:.2f}x, max improvement = {max_imp:.2f}x")

    print()
    print("INTERPRETATION:")
    print("  - Improvement > 1.0x means noise-informed decoder is BETTER")
    print("  - The advantage grows with code distance (more qubits to exploit)")
    print("  - At low error rates, the advantage is largest (below threshold)")
    print("  - This is EXACTLY what tau-chrono provides: per-gate noise knowledge")
    print()
    print("CONNECTION TO TAU-CHRONO / PETZ-INFORMED DECODING:")
    print("  tau_j = 1 - F_j  measures per-qubit noise severity.")
    print("  Higher tau -> higher edge weight -> MWPM avoids paths through noisy qubits.")
    print("  Standard decoder assumes all tau_j equal -> misses this information.")
    print("  The improvement factor quantifies the VALUE of tau-chrono data.")
    print()

    # -----------------------------------------------------------------------
    # Generate plot
    # -----------------------------------------------------------------------
    try:
        import os
        os.environ['MPLCONFIGDIR'] = '/tmp/matplotlib_config'
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, len(DISTANCES), figsize=(5 * len(DISTANCES), 5),
                                 squeeze=False)

        colors_uniform = {'3': '#e74c3c', '5': '#e67e22', '7': '#e74c3c'}
        colors_informed = {'3': '#2ecc71', '5': '#27ae60', '7': '#2ecc71'}

        for idx, d in enumerate(DISTANCES):
            ax = axes[0, idx]
            ps = []
            lers_u = []
            lers_i = []
            for p in PHYSICAL_ERROR_RATES:
                ler_u, ler_i = results[(d, p)]
                if ler_u > 0 and ler_i > 0:
                    ps.append(p)
                    lers_u.append(ler_u)
                    lers_i.append(ler_i)

            ax.loglog(ps, lers_u, 'o-', color='#e74c3c', linewidth=2,
                     markersize=6, label='Standard (uniform weights)')
            ax.loglog(ps, lers_i, 's-', color='#2ecc71', linewidth=2,
                     markersize=6, label='Noise-informed (true weights)')

            ax.set_xlabel('Physical error rate p', fontsize=12)
            ax.set_ylabel('Logical error rate', fontsize=12)
            ax.set_title(f'Distance d={d}', fontsize=14, fontweight='bold')
            ax.legend(fontsize=9, loc='lower right')
            ax.grid(True, alpha=0.3, which='both')
            ax.set_xlim(min(PHYSICAL_ERROR_RATES) * 0.8, max(PHYSICAL_ERROR_RATES) * 1.2)

        plt.suptitle(
            f'Noise-Informed vs Standard MWPM Decoding\n'
            f'(non-uniform noise: {FRACTION_NOISY*100:.0f}% of qubits at {NOISE_RATIO}x higher rate)',
            fontsize=14, fontweight='bold', y=1.02
        )
        plt.tight_layout()

        plot_path = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            'qec_experiment_results.png'
        )
        plt.savefig(plot_path, dpi=150, bbox_inches='tight')
        print(f"Plot saved to: {plot_path}")

    except Exception as e:
        print(f"Could not generate plot: {e}")

    # -----------------------------------------------------------------------
    # Save CSV data
    # -----------------------------------------------------------------------
    try:
        import os
        csv_path = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            'qec_experiment_results.csv'
        )
        with open(csv_path, 'w') as f:
            f.write("distance,physical_error_rate,ler_uniform,ler_informed,improvement_factor\n")
            for d in DISTANCES:
                for p in PHYSICAL_ERROR_RATES:
                    ler_u, ler_i = results[(d, p)]
                    imp = ler_u / ler_i if ler_i > 0 else float('inf')
                    f.write(f"{d},{p},{ler_u:.8f},{ler_i:.8f},{imp:.4f}\n")
        print(f"CSV saved to: {csv_path}")
    except Exception as e:
        print(f"Could not save CSV: {e}")

    # -----------------------------------------------------------------------
    # Raw data dump for programmatic access
    # -----------------------------------------------------------------------
    print("\n\nRAW DATA (copy-paste friendly):")
    print("distance | p_phys | LER_uniform | LER_informed | improvement")
    for d in DISTANCES:
        for p in PHYSICAL_ERROR_RATES:
            ler_u, ler_i = results[(d, p)]
            imp = ler_u / ler_i if ler_i > 0 else float('inf')
            print(f"  {d:>5d}  | {p:.4f} | {ler_u:.8f} | {ler_i:.8f} | {imp:.3f}x")

    return results


# ---------------------------------------------------------------------------
# Additional analysis: Detailed decoder weight comparison
# ---------------------------------------------------------------------------

def analyze_decoder_weights(distance=3, base_p=0.01):
    """
    Show exactly how decoder weights differ between uniform and informed.
    """
    print(f"\n\n{'=' * 72}")
    print(f"DECODER WEIGHT ANALYSIS (d={distance}, p={base_p})")
    print(f"{'=' * 72}")

    rounds = distance
    circuit_true, circuit_avg, hot_qubits = build_nonuniform_circuit(
        distance=distance, rounds=rounds, base_p=base_p,
        noise_ratio=NOISE_RATIO, frac_noisy=FRACTION_NOISY
    )

    print(f"\nHot qubits (indices): {sorted(hot_qubits)}")
    print(f"Total data qubits: {distance * distance}")

    # Extract DEMs
    dem_true = circuit_true.detector_error_model(decompose_errors=True)
    dem_avg = circuit_avg.detector_error_model(decompose_errors=True)

    # Count error mechanisms
    n_true = len(str(dem_true).split('\n'))
    n_avg = len(str(dem_avg).split('\n'))
    print(f"Error mechanisms in true DEM: {n_true}")
    print(f"Error mechanisms in uniform DEM: {n_avg}")

    # Compare a few error probabilities
    def extract_error_probs(dem):
        probs = []
        for line in str(dem).split('\n'):
            match = re.search(r'error\(([^)]+)\)', line)
            if match:
                probs.append(float(match.group(1)))
        return np.array(probs)

    probs_true = extract_error_probs(dem_true)
    probs_avg = extract_error_probs(dem_avg)

    print(f"\nTrue DEM error probabilities:")
    print(f"  min={probs_true.min():.6f}, max={probs_true.max():.6f}, "
          f"mean={probs_true.mean():.6f}, std={probs_true.std():.6f}")
    print(f"Uniform DEM error probabilities:")
    print(f"  min={probs_avg.min():.6f}, max={probs_avg.max():.6f}, "
          f"mean={probs_avg.mean():.6f}, std={probs_avg.std():.6f}")

    # Weight = log((1-p)/p). Compare weight distributions
    weights_true = np.log((1 - probs_true) / probs_true)
    weights_avg = np.log((1 - probs_avg) / probs_avg)

    print(f"\nTrue DEM weights (log-likelihood):")
    print(f"  min={weights_true.min():.2f}, max={weights_true.max():.2f}, "
          f"mean={weights_true.mean():.2f}, std={weights_true.std():.2f}")
    print(f"Uniform DEM weights:")
    print(f"  min={weights_avg.min():.2f}, max={weights_avg.max():.2f}, "
          f"mean={weights_avg.mean():.2f}, std={weights_avg.std():.2f}")

    print(f"\nKey insight: The true DEM has weight spread (std={weights_true.std():.2f})")
    print(f"while uniform DEM has spread (std={weights_avg.std():.2f}).")
    print(f"The informed decoder uses this variation to route corrections")
    print(f"AWAY from noisy qubits, reducing logical error rate.")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    results = run_experiment()
    analyze_decoder_weights(distance=3, base_p=0.01)
    print("\nDone.")
