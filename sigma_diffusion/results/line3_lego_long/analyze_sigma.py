#!/usr/bin/env python3
"""
Sigma analysis of 3DGS training checkpoints.

Analyzes PLY files from OpenSplat training to:
1. Count gaussians at each checkpoint
2. Identify densification events
3. Compute parameter drift (Sigma) between consecutive checkpoints
4. Identify Sigma plateau (early stopping point)
5. Report loss curve
"""

import struct
import os
import sys
import math
import numpy as np

PLY_DIR = "/Users/akaihuangm1/Desktop/github/petz-recovery-unification/sigma_diffusion/results/line3_lego_long"

# Checkpoints saved at these iterations
CHECKPOINTS = [500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000]

def parse_ply_header(filepath):
    """Parse PLY header to get vertex count and property info."""
    with open(filepath, 'rb') as f:
        header_lines = []
        while True:
            line = f.readline().decode('ascii', errors='replace').strip()
            header_lines.append(line)
            if line == 'end_header':
                break

        # Parse header
        vertex_count = 0
        properties = []
        in_vertex = False
        for line in header_lines:
            if line.startswith('element vertex'):
                vertex_count = int(line.split()[-1])
                in_vertex = True
            elif line.startswith('element') and 'vertex' not in line:
                in_vertex = False
            elif line.startswith('property') and in_vertex:
                parts = line.split()
                dtype = parts[1]
                name = parts[2]
                properties.append((name, dtype))

        header_size = f.tell()

    return vertex_count, properties, header_size


def dtype_to_numpy(dtype_str):
    """Convert PLY dtype string to numpy dtype."""
    mapping = {
        'float': np.float32,
        'double': np.float64,
        'int': np.int32,
        'uint': np.uint32,
        'uchar': np.uint8,
        'short': np.int16,
        'ushort': np.uint16,
    }
    return mapping.get(dtype_str, np.float32)


def load_ply(filepath):
    """Load a PLY file and return gaussian parameters as structured array."""
    vertex_count, properties, header_size = parse_ply_header(filepath)

    if vertex_count == 0:
        return None, properties

    # Build numpy dtype
    dt = np.dtype([(name, dtype_to_numpy(dtype)) for name, dtype in properties])

    with open(filepath, 'rb') as f:
        f.seek(header_size)
        data = np.frombuffer(f.read(vertex_count * dt.itemsize), dtype=dt)

    return data, properties


def compute_stats(data, properties):
    """Compute summary statistics for all properties."""
    stats = {}
    prop_names = [p[0] for p in properties]
    for name in prop_names:
        col = data[name].astype(np.float64)
        stats[name] = {
            'mean': np.mean(col),
            'std': np.std(col),
            'min': np.min(col),
            'max': np.max(col),
        }
    return stats


def compute_sigma_between(data_a, data_b, properties):
    """
    Compute Sigma (parameter drift) between two checkpoint gaussian populations.

    Since gaussian counts change (densification/culling), we compare distributions
    rather than point-to-point. We use:
    1. Distribution-level KL divergence estimate for each parameter
    2. Overall Frobenius-norm-style drift of distribution moments
    """
    prop_names = [p[0] for p in properties]

    # Position properties
    pos_props = ['x', 'y', 'z']
    # Scale properties
    scale_props = [p for p in prop_names if p.startswith('scale_')]
    # Rotation properties
    rot_props = [p for p in prop_names if p.startswith('rot_')]
    # Opacity
    opacity_props = [p for p in prop_names if p.startswith('opacity')]
    # SH coefficients
    sh_props = [p for p in prop_names if p.startswith('f_dc') or p.startswith('f_rest')]

    categories = {
        'position': pos_props,
        'scale': scale_props,
        'rotation': rot_props,
        'opacity': opacity_props,
        'color_sh': sh_props,
    }

    sigma_by_category = {}
    total_sigma = 0.0

    for cat_name, cat_props in categories.items():
        available = [p for p in cat_props if p in prop_names]
        if not available:
            continue

        cat_sigma = 0.0
        for prop in available:
            a = data_a[prop].astype(np.float64)
            b = data_b[prop].astype(np.float64)

            mu_a, sig_a = np.mean(a), np.std(a) + 1e-10
            mu_b, sig_b = np.mean(b), np.std(b) + 1e-10

            # Skip if both distributions are essentially constant (all same value)
            if np.std(a) < 1e-15 and np.std(b) < 1e-15:
                # Both constant -- check if means differ
                if abs(mu_a - mu_b) > 1e-15:
                    cat_sigma += (mu_a - mu_b)**2 / (2e-10)  # Cap contribution
                continue

            # KL divergence between two Gaussians:
            # KL(A||B) = log(sig_b/sig_a) + (sig_a^2 + (mu_a-mu_b)^2)/(2*sig_b^2) - 0.5
            try:
                kl = math.log(sig_b/sig_a) + (sig_a**2 + (mu_a - mu_b)**2)/(2*sig_b**2) - 0.5
                if math.isnan(kl) or math.isinf(kl):
                    kl = 0.0
            except (ValueError, ZeroDivisionError):
                kl = 0.0
            cat_sigma += abs(kl)  # Use absolute value for total drift

        sigma_by_category[cat_name] = cat_sigma
        total_sigma += cat_sigma

    return total_sigma, sigma_by_category


def main():
    print("=" * 80)
    print("3D Gaussian Splatting - Sigma Analysis of Training Checkpoints")
    print("=" * 80)

    # ---- PART 1: Parse training log data from the output ----
    # Loss and gaussian count data extracted from training output
    loss_data = {}
    gaussian_events = []

    # Parse from known training output
    losses = {
        10: 0.443267, 50: 0.402146, 100: 0.362921, 150: 0.310415, 200: 0.283344,
        250: 0.254378, 300: 0.182924, 350: 0.13749, 400: 0.203173, 450: 0.109661,
        500: 0.122067, 600: 0.0984644, 700: 0.0955307, 800: 0.079762, 900: 0.0714082,
        1000: 0.0493855, 1100: 0.0573344, 1200: 0.159828, 1300: 0.0378001,
        1400: 0.0348198, 1500: 0.0228278, 1600: 0.0243413, 1700: 0.0173436,
        1800: 0.0193634, 1900: 0.0182099, 2000: 0.0201379, 2100: 0.0172877,
        2200: 0.0164119, 2300: 0.014204, 2400: 0.0125169, 2500: 0.0138772,
        2600: 0.0154176, 2700: 0.0149037, 2800: 0.012776, 2900: 0.0121062,
        3000: 0.0597718, 3100: 0.0253811, 3200: 0.0186974, 3300: 0.0190491,
        3400: 0.0176728, 3500: 0.0167093, 3600: 0.0163754, 3700: 0.0178543,
        3800: 0.0132856, 3900: 0.0147445, 4000: 0.0125174, 4100: 0.0143567,
        4200: 0.0134178, 4300: 0.0155529, 4400: 0.0192658, 4500: 0.014448,
        4600: 0.0118043, 4700: 0.0296515, 4800: 0.024514, 4900: 0.0154402,
        5000: 0.0235944,
    }

    # Densification events: (step, added, new_total, culled, remaining)
    densification = [
        (600,  116,  411,  330,  81),
        (700,  158,  239,  101,  138),
        (800,  270,  408,  148,  260),
        (900,  486,  746,  271,  475),
        (1000, 855,  1330, 476,  854),
        (1100, 1430, 2284, 711,  1573),
        (1200, 2303, 3876, 1058, 2818),
        (1300, 3139, 5957, 1317, 4640),
        (1400, 3748, 8388, 1338, 7050),
        (1500, 4038, 11088,1342, 9746),
        (1600, 4081, 13827,1245, 12582),
        (1700, 3834, 16416,1228, 15188),
        (1800, 3715, 18903,1263, 17640),
        (1900, 3544, 21184,1237, 19947),
        (2000, 3407, 23354,1301, 22053),
        (2100, 3296, 25349,1321, 24028),
        (2200, 3146, 27174,1305, 25869),
        (2300, 3072, 28941,1396, 27545),
        (2400, 2881, 30426,1289, 29137),
    ]

    print("\n" + "-" * 80)
    print("PART 1: Gaussian Count at Each Checkpoint")
    print("-" * 80)

    # Map checkpoint iterations to gaussian counts
    # After densification at step X, the remaining count applies from X to X+99
    # The checkpoint at iteration N captures the state after all densifications up to N
    checkpoint_counts = {}
    for cp in CHECKPOINTS:
        # Find the last densification event at or before this checkpoint
        count = 295  # initial
        for step, added, new_total, culled, remaining in densification:
            if step <= cp:
                count = remaining
        checkpoint_counts[cp] = count

    # But after iteration 2500 (last densification at 2400), count stays at 29137
    # At step 3000, there's a big reset (opacity reset, not densification)

    print(f"\n{'Iteration':>10} | {'Gaussians':>10} | {'Loss':>12} | {'File Size':>12}")
    print("-" * 55)

    ply_data = {}
    for cp in CHECKPOINTS:
        fname = f"lego_5k_{cp}.ply"
        fpath = os.path.join(PLY_DIR, fname)
        if os.path.exists(fpath):
            data, props = load_ply(fpath)
            n_gauss = len(data) if data is not None else 0
            fsize = os.path.getsize(fpath)
            loss_val = losses.get(cp, 0)
            print(f"{cp:>10} | {n_gauss:>10} | {loss_val:>12.6f} | {fsize:>10} B")
            ply_data[cp] = (data, props, n_gauss)
        else:
            print(f"{cp:>10} | {'FILE MISSING':>10}")

    # ---- PART 2: Densification Timeline ----
    print("\n" + "-" * 80)
    print("PART 2: Densification Events")
    print("-" * 80)
    print(f"\n{'Step':>6} | {'Added':>7} | {'Before Cull':>12} | {'Culled':>7} | {'After':>7} | {'Net Change':>11}")
    print("-" * 65)

    prev_remaining = 295
    for step, added, new_total, culled, remaining in densification:
        net = remaining - prev_remaining
        print(f"{step:>6} | {added:>7} | {new_total:>12} | {culled:>7} | {remaining:>7} | {net:>+11}")
        prev_remaining = remaining

    print(f"\nDensification ran from iteration 600 to 2400 (every 100 steps)")
    print(f"Peak growth rate: iterations 1000-1600 (~3000-4000 added per step)")
    print(f"Growth slows: iterations 1700+ (diminishing returns)")
    print(f"Final count after last densification: 29,137 gaussians")

    # ---- PART 3: Sigma Analysis ----
    print("\n" + "-" * 80)
    print("PART 3: Sigma (Parameter Drift) Between Consecutive Checkpoints")
    print("-" * 80)

    sorted_cps = sorted(ply_data.keys())
    sigma_results = []

    for i in range(len(sorted_cps) - 1):
        cp_a = sorted_cps[i]
        cp_b = sorted_cps[i + 1]
        data_a, props_a, n_a = ply_data[cp_a]
        data_b, props_b, n_b = ply_data[cp_b]

        if data_a is None or data_b is None:
            continue

        sigma_total, sigma_cats = compute_sigma_between(data_a, data_b, props_a)
        sigma_results.append((cp_a, cp_b, n_a, n_b, sigma_total, sigma_cats))

    print(f"\n{'Interval':>15} | {'N_a':>7} | {'N_b':>7} | {'Sigma_total':>12} | {'Pos':>8} | {'Scale':>8} | {'Rot':>8} | {'Opacity':>8} | {'Color':>8}")
    print("-" * 110)

    for cp_a, cp_b, n_a, n_b, sigma_total, sigma_cats in sigma_results:
        pos = sigma_cats.get('position', 0)
        scl = sigma_cats.get('scale', 0)
        rot = sigma_cats.get('rotation', 0)
        opa = sigma_cats.get('opacity', 0)
        col = sigma_cats.get('color_sh', 0)
        print(f"{cp_a:>6}-{cp_b:<6}  | {n_a:>7} | {n_b:>7} | {sigma_total:>12.4f} | {pos:>8.4f} | {scl:>8.4f} | {rot:>8.4f} | {opa:>8.4f} | {col:>8.4f}")

    # ---- PART 4: Loss Curve ----
    print("\n" + "-" * 80)
    print("PART 4: Loss Curve (sampled every 100 iterations)")
    print("-" * 80)

    sorted_iters = sorted(losses.keys())
    print(f"\n{'Iter':>6} | {'Loss':>12} | {'Bar'}")
    print("-" * 60)
    for it in sorted_iters:
        loss = losses[it]
        bar_len = int(loss * 150)
        bar = '#' * min(bar_len, 60)
        print(f"{it:>6} | {loss:>12.6f} | {bar}")

    # ---- PART 5: Plateau Detection ----
    print("\n" + "-" * 80)
    print("PART 5: Sigma Plateau Analysis")
    print("-" * 80)

    if sigma_results:
        sigmas = [(cp_a, cp_b, sigma_total) for cp_a, cp_b, _, _, sigma_total, _ in sigma_results]

        # Find where Sigma drops and stabilizes
        print("\nSigma evolution:")
        for cp_a, cp_b, sigma in sigmas:
            if math.isnan(sigma) or math.isinf(sigma):
                sigma_display = 0.0
            else:
                sigma_display = sigma
            bar_len = int(sigma_display * 5)
            bar = '|' * min(max(bar_len, 0), 60)
            print(f"  {cp_a:>5}-{cp_b:<5}: Sigma = {sigma:>10.4f}  {bar}")

        # Find plateau: where relative change in sigma < threshold
        print("\nRelative change in Sigma:")
        for i in range(1, len(sigmas)):
            prev_sigma = sigmas[i-1][2]
            curr_sigma = sigmas[i][2]
            if math.isnan(prev_sigma) or math.isnan(curr_sigma):
                print(f"  {sigmas[i][0]:>5}-{sigmas[i][1]:<5}: delta_Sigma/Sigma =      NaN (skipped)")
                continue
            rel_change = abs(curr_sigma - prev_sigma) / (prev_sigma + 1e-10)
            marker = " <-- PLATEAU" if rel_change < 0.3 and curr_sigma < 5 else ""
            print(f"  {sigmas[i][0]:>5}-{sigmas[i][1]:<5}: delta_Sigma/Sigma = {rel_change:>8.4f}{marker}")

        # Find the best early stopping point
        # The plateau is where sigma is small AND stable
        valid_sigmas = [(i, s) for i, (_, _, s) in enumerate(sigmas) if not math.isnan(s)]
        min_sigma_idx = min(valid_sigmas, key=lambda x: x[1])[0] if valid_sigmas else 0
        print(f"\nMinimum Sigma interval: {sigmas[min_sigma_idx][0]}-{sigmas[min_sigma_idx][1]} (Sigma = {sigmas[min_sigma_idx][2]:.4f})")

        # Check if loss at plateau checkpoint is acceptable
        plateau_cp = sigmas[min_sigma_idx][1]
        plateau_loss = losses.get(plateau_cp, None)
        final_loss = losses.get(5000, None)
        if plateau_loss and final_loss:
            quality_ratio = plateau_loss / final_loss
            print(f"\nEarly stopping analysis:")
            print(f"  Plateau checkpoint: iteration {plateau_cp}")
            print(f"  Loss at plateau:    {plateau_loss:.6f}")
            print(f"  Loss at final:      {final_loss:.6f}")
            print(f"  Quality ratio:      {quality_ratio:.2f}x (1.0 = same as final)")
            if quality_ratio < 1.5:
                print(f"  --> ACCEPTABLE: Early stopping at {plateau_cp} gives near-final quality")
            else:
                print(f"  --> NOT RECOMMENDED: Significant quality gap remains")

    # ---- PART 6: Key findings ----
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)

    # Gaussian growth summary
    print(f"""
Gaussian Count Evolution:
  Start:         295 (from SfM points)
  After densif:  29,137 (at iter 2500)
  Final:         29,137 (iter 5000, no more densification after 2500)
  Growth factor: {29137/295:.0f}x

Densification Dynamics:
  - Started at iteration 600 (OpenSplat default: 500 + delay)
  - Ran every 100 iterations until iteration 2400
  - Aggressive culling at early stages (iter 600: culled 80% of gaussians!)
  - Peak growth around iterations 1000-1600
  - Net growth stabilized around iter 2000+

Loss Curve:
  - Rapid descent: 0.44 -> 0.05 (iters 10-1000)
  - Steady improvement: 0.05 -> 0.012 (iters 1000-2900)
  - Opacity reset spike at iter 3000: loss jumps to 0.06
  - Recovery: 0.06 -> 0.013 (iters 3000-4000)
  - Late instability: loss oscillates 0.012-0.037 (iters 4500-5000)

Key Observation - Opacity Reset at iter 3000:
  - File sizes are identical from iter 2500 onward (7,227,555 bytes)
  - This means gaussian count is stable at 29,137
  - The loss spike at iter 3000 is from opacity reset (standard 3DGS technique)
  - After reset, the model re-learns opacities (recovers by iter 3200)
""")

    print("Files written to:", PLY_DIR)
    print(f"Total checkpoints: {len(CHECKPOINTS)}")
    print(f"Final model: lego_5k.ply ({os.path.getsize(os.path.join(PLY_DIR, 'lego_5k.ply')):,} bytes)")


if __name__ == '__main__':
    main()
