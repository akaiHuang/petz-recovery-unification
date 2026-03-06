#!/usr/bin/env python3
"""
Extended Decoder Alpha-Independence Test (Large d)
===================================================

Extends the alpha-independence test from d=3,5,7 to d=3,5,7,9,11.
Tests with 3 decoders: MWPM, Unwt-MWPM, Union-Find.
Skips slow decoders (BP+OSD, Greedy-NN) to keep runtime manageable.

Prediction 1: ln R(eps,d) = alpha(p)*d + beta(eps),
where alpha depends only on physical error rate p (NOT on decoder choice).

Previous test at d=3,5,7 with 6 decoders gave CV=0.46 (inconclusive),
but among 4 functional decoders pairwise differences were <=1.4 sigma.

Author: Sheng-Kai Huang
Date: 2026-03-06
"""

import numpy as np
import re
import time
import os
import sys
from collections import defaultdict
from typing import Dict, Tuple, List, Optional

import stim
import pymatching

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D


def flush_print(*args, **kwargs):
    """Print with immediate flush to avoid output buffering."""
    print(*args, **kwargs, flush=True)


# ============================================================================
# Configuration
# ============================================================================

DISTANCES = [3, 5, 7, 9, 11]
P_VALUES = [0.003, 0.005, 0.007, 0.01]

# Post-selection parameter
PS_EPSILON = 0.20  # discard top 20% syndrome weight

# Shots per (d, p, decoder) -- scaled for speed vs statistics
SHOTS_FAST = {3: 100_000, 5: 100_000, 7: 100_000, 9: 50_000, 11: 30_000}
SHOTS_UF   = {3: 10_000,  5: 10_000,  7: 5_000,   9: 3_000,  11: 2_000}

DECODER_NAMES_ORDERED = ['MWPM', 'Unwt-MWPM', 'Union-Find']

# Output paths
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FIG_PATH_PNG = os.path.join(SCRIPT_DIR, 'fig_alpha_independence_large_d.png')
FIG_PATH_PDF = os.path.join(SCRIPT_DIR, 'fig_alpha_independence_large_d.pdf')


# ============================================================================
# Helper: DEM -> Check Matrix extraction
# ============================================================================

def dem_to_check_matrix(dem: stim.DetectorErrorModel, n_det: int, n_obs: int = 1):
    """
    Extract parity check matrix, observable matrix, and error probabilities
    from a stim DetectorErrorModel.
    """
    errors = []
    for inst in dem.flattened():
        if inst.type == 'error':
            targets = inst.targets_copy()
            p = inst.args_copy()[0]
            dets = [t.val for t in targets if t.is_relative_detector_id()]
            obs = [t.val for t in targets if t.is_logical_observable_id()]
            errors.append((p, dets, obs))

    n_errors = len(errors)
    H = np.zeros((n_det, n_errors), dtype=np.uint8)
    obs_matrix = np.zeros((n_obs, n_errors), dtype=np.uint8)
    probs = np.zeros(n_errors)

    for j, (p, dets, obs_list) in enumerate(errors):
        probs[j] = max(p, 1e-15)
        for d_idx in dets:
            H[d_idx, j] = 1
        for o in obs_list:
            obs_matrix[o, j] = 1

    llrs = np.log((1 - probs) / probs)
    return H, obs_matrix, probs, llrs


# ============================================================================
# Decoder Implementations
# ============================================================================

class MWPMDecoder:
    """Minimum-Weight Perfect Matching decoder via pymatching."""
    def __init__(self, dem, circuit=None, **kwargs):
        self.matching = pymatching.Matching.from_detector_error_model(dem)
        self.name = "MWPM"

    def decode_batch(self, detections, obs_flat):
        pred = self.matching.decode_batch(detections)
        return pred.flatten() != obs_flat


class UnweightedMWPMDecoder:
    """MWPM with uniform edge weights."""
    UNIFORM_P = 0.01

    def __init__(self, dem, circuit=None, **kwargs):
        self.name = "Unwt-MWPM"
        dem_str = str(dem)
        uniform_dem_str = re.sub(
            r'error\([^)]+\)',
            f'error({self.UNIFORM_P})',
            dem_str
        )
        uniform_dem = stim.DetectorErrorModel(uniform_dem_str)
        self.matching = pymatching.Matching.from_detector_error_model(uniform_dem)

    def decode_batch(self, detections, obs_flat):
        pred = self.matching.decode_batch(detections)
        return pred.flatten() != obs_flat


class UnionFindDecoder:
    """Union-Find decoder via ldpc package."""
    def __init__(self, dem, circuit=None, **kwargs):
        from ldpc import UnionFindDecoder as LdpcUF
        self.name = "Union-Find"
        n_det = kwargs.get('n_det', circuit.num_detectors if circuit else None)
        n_obs = kwargs.get('n_obs', 1)
        H, obs_matrix, probs, llrs = dem_to_check_matrix(dem, n_det, n_obs)
        self.H = H
        self.obs_matrix = obs_matrix
        self.llrs = llrs
        self.uf = LdpcUF(H, uf_method='peeling')

    def decode_batch(self, detections, obs_flat):
        num_shots = detections.shape[0]
        errors = np.zeros(num_shots, dtype=bool)
        for i in range(num_shots):
            syndrome = detections[i].astype(np.uint8)
            estimate = self.uf.decode(syndrome, self.llrs)
            pred_obs = (self.obs_matrix @ estimate) % 2
            errors[i] = (pred_obs[0] != obs_flat[i])
        return errors


# ============================================================================
# Simulation Engine
# ============================================================================

def generate_circuit(distance: int, p: float, rounds: Optional[int] = None) -> stim.Circuit:
    """Generate a rotated surface code memory-Z experiment circuit."""
    if rounds is None:
        rounds = distance
    return stim.Circuit.generated(
        'surface_code:rotated_memory_z',
        rounds=rounds,
        distance=distance,
        after_clifford_depolarization=p,
        after_reset_flip_probability=p,
        before_measure_flip_probability=p,
        before_round_data_depolarization=p,
    )


def compute_wilson_ci(errors: np.ndarray) -> Tuple[float, float]:
    """Compute error rate with 95% Wilson score confidence interval half-width."""
    n = len(errors)
    if n == 0:
        return np.nan, np.nan
    rate = errors.mean()
    z = 1.96
    denom = 1 + z**2 / n
    hw = z * np.sqrt(rate * (1 - rate) / n + z**2 / (4 * n**2)) / denom
    return rate, hw


def compute_R_with_error(base_rate, base_n, ps_rate, ps_n):
    """
    Compute improvement ratio R = base_rate / ps_rate with propagated uncertainty.
    Uses delta method for ratio of two binomial proportions.
    """
    if ps_rate <= 0 or base_rate <= 0:
        return np.nan, np.nan
    R = base_rate / ps_rate

    se_base = np.sqrt(base_rate * (1 - base_rate) / base_n) if base_n > 0 else 0
    se_ps = np.sqrt(ps_rate * (1 - ps_rate) / ps_n) if ps_n > 0 else 0

    rel_var = (se_base / base_rate)**2 + (se_ps / ps_rate)**2
    R_se = R * np.sqrt(rel_var)
    return R, 1.96 * R_se  # 95% CI


# ============================================================================
# Main Simulation
# ============================================================================

def run_simulation():
    """
    Run the full alpha-independence test with extended distances.

    For each (d, p):
      1. Generate circuit and sample detection events
      2. Decode with each decoder (separate shot counts for fast vs slow)
      3. Compute logical error rate with and without post-selection
      4. Compute ln R for each (decoder, d, p) combination
    """
    flush_print()
    flush_print("=" * 75)
    flush_print("  EXTENDED DECODER ALPHA-INDEPENDENCE TEST (Large d)")
    flush_print("  Prediction 1: ln R(eps,d) = alpha(p)*d + beta(eps)")
    flush_print("  alpha should be INDEPENDENT of decoder choice")
    flush_print("=" * 75)
    flush_print()

    # results[decoder_name][d][p] = {
    #   'base_rate', 'base_ci', 'base_n',
    #   'ps_rate', 'ps_ci', 'ps_n', 'R', 'R_ci', 'lnR', 'lnR_err'
    # }
    results = {}
    for name in DECODER_NAMES_ORDERED:
        results[name] = {}

    total_configs = len(DISTANCES) * len(P_VALUES)
    config_idx = 0

    for d in DISTANCES:
        for p in P_VALUES:
            config_idx += 1
            flush_print(f"  [{config_idx}/{total_configs}] d={d}, p={p:.4f}")

            # Generate circuit
            circuit = generate_circuit(d, p)
            dem = circuit.detector_error_model(decompose_errors=True)
            sampler = circuit.compile_detector_sampler()
            n_det = circuit.num_detectors
            n_obs = circuit.num_observables
            kwargs = {'n_det': n_det, 'n_obs': n_obs}

            # Sample max shots needed (fast decoders need more)
            max_shots = SHOTS_FAST[d]
            t0 = time.time()
            det_all, obs_all = sampler.sample(max_shots, separate_observables=True)
            obs_flat_all = obs_all.flatten()
            weights_all = det_all.sum(axis=1).astype(float)
            t_sample = time.time() - t0
            flush_print(f"    Sampled {max_shots:,} shots in {t_sample:.1f}s")

            # Build decoders
            t_build = time.time()
            decoders = [
                MWPMDecoder(dem, circuit, **kwargs),
                UnweightedMWPMDecoder(dem, circuit, **kwargs),
                UnionFindDecoder(dem, circuit, **kwargs),
            ]
            t_build = time.time() - t_build

            for decoder in decoders:
                name = decoder.name
                if d not in results[name]:
                    results[name][d] = {}

                # Determine shot count for this decoder
                if name == 'Union-Find':
                    n_use = min(SHOTS_UF[d], max_shots)
                else:
                    n_use = max_shots

                det_use = det_all[:n_use]
                obs_use = obs_flat_all[:n_use]
                w_use = weights_all[:n_use]

                t_dec = time.time()
                is_error = decoder.decode_batch(det_use, obs_use)
                t_dec = time.time() - t_dec

                base_rate, base_ci = compute_wilson_ci(is_error)

                # Post-selection: keep bottom (1 - eps) fraction by syndrome weight
                cutoff = np.quantile(w_use, 1.0 - PS_EPSILON)
                mask = w_use <= cutoff
                n_kept = mask.sum()

                if n_kept > 0 and is_error[mask].size > 0:
                    ps_rate, ps_ci = compute_wilson_ci(is_error[mask])
                    # Floor: if ps_rate is 0, set to 1/n_kept
                    if ps_rate == 0:
                        ps_rate = 1.0 / n_kept
                        ps_ci = 0.0
                    R, R_ci = compute_R_with_error(base_rate, n_use, ps_rate, n_kept)
                else:
                    ps_rate, ps_ci = np.nan, np.nan
                    R, R_ci = np.nan, np.nan

                # Compute ln R
                if not np.isnan(R) and R > 0:
                    lnR = np.log(R)
                    if R_ci > 0 and R > 0:
                        lnR_err = R_ci / R  # delta method
                    else:
                        lnR_err = 0.1
                else:
                    lnR = np.nan
                    lnR_err = np.nan

                results[name][d][p] = {
                    'base_rate': base_rate, 'base_ci': base_ci, 'base_n': n_use,
                    'ps_rate': ps_rate, 'ps_ci': ps_ci, 'ps_n': n_kept,
                    'R': R, 'R_ci': R_ci,
                    'lnR': lnR, 'lnR_err': lnR_err,
                }

                flush_print(f"    {name:14s}: p_L={base_rate:.6f}, "
                            f"p_L_ps={ps_rate:.6f}, "
                            f"R={R:.3f}, ln R={lnR:.4f}  "
                            f"({t_dec:.1f}s, {n_use:,} shots)")

            flush_print()

    return results


# ============================================================================
# Alpha Extraction via Linear Fit
# ============================================================================

def extract_alphas(results, p_fixed=None):
    """
    For each decoder, fit ln R = alpha * d + beta across distances.

    If p_fixed is given, fit at that single p value.
    Returns dict: decoder_name -> (alpha, alpha_err, beta, beta_err, d_values, lnR_values, lnR_errors)
    """
    alphas = {}

    for name in DECODER_NAMES_ORDERED:
        d_values = []
        lnR_values = []
        lnR_errors = []

        p_list = [p_fixed] if p_fixed is not None else P_VALUES

        for d in DISTANCES:
            if d not in results[name]:
                continue
            for p in p_list:
                if p not in results[name][d]:
                    continue
                data = results[name][d][p]
                lnR = data['lnR']
                lnR_err = data['lnR_err']
                if np.isnan(lnR) or np.isnan(lnR_err):
                    continue

                d_values.append(d)
                lnR_values.append(lnR)
                lnR_errors.append(lnR_err)

        if len(d_values) < 2:
            continue

        d_arr = np.array(d_values, dtype=float)
        lnR_arr = np.array(lnR_values)
        lnR_err_arr = np.array(lnR_errors)

        # Weighted linear fit: lnR = alpha * d + beta
        w = 1.0 / (lnR_err_arr**2 + 1e-10)
        W = w.sum()
        Wd = (w * d_arr).sum()
        WlnR = (w * lnR_arr).sum()
        Wdd = (w * d_arr**2).sum()
        WdlnR = (w * d_arr * lnR_arr).sum()

        det = W * Wdd - Wd**2
        if abs(det) > 1e-20:
            alpha = (W * WdlnR - Wd * WlnR) / det
            beta = (Wdd * WlnR - Wd * WdlnR) / det
            alpha_err = np.sqrt(W / det)
            beta_err = np.sqrt(Wdd / det)
        else:
            coeffs = np.polyfit(d_arr, lnR_arr, 1)
            alpha = coeffs[0]
            beta = coeffs[1]
            alpha_err = 0.0
            beta_err = 0.0

        alphas[name] = {
            'alpha': alpha, 'alpha_err': alpha_err,
            'beta': beta, 'beta_err': beta_err,
            'd_values': d_arr, 'lnR_values': lnR_arr,
            'lnR_errors': lnR_err_arr,
        }

    return alphas


def compute_cv_and_pairwise(alphas):
    """Compute CV and pairwise sigma-differences from alpha dict."""
    names = [n for n in DECODER_NAMES_ORDERED if n in alphas]
    if len(names) < 2:
        return None, None, names

    alpha_vals = np.array([alphas[n]['alpha'] for n in names])
    alpha_errs = np.array([alphas[n]['alpha_err'] for n in names])

    mean_a = alpha_vals.mean()
    std_a = alpha_vals.std()
    cv = std_a / abs(mean_a) if abs(mean_a) > 1e-10 else float('inf')

    pairwise = []
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            diff = abs(alpha_vals[i] - alpha_vals[j])
            combined_err = np.sqrt(alpha_errs[i]**2 + alpha_errs[j]**2)
            sig = diff / combined_err if combined_err > 0 else float('inf')
            pairwise.append((names[i], names[j], alpha_vals[i], alpha_vals[j],
                             diff, combined_err, sig))

    return cv, pairwise, names


# ============================================================================
# Analysis and Printing
# ============================================================================

def print_analysis(results, alphas_by_p):
    """Print detailed analysis of alpha-independence."""
    flush_print()
    flush_print("=" * 75)
    flush_print("  ALPHA-INDEPENDENCE ANALYSIS")
    flush_print("=" * 75)
    flush_print()

    for p in P_VALUES:
        alphas = alphas_by_p[p]
        cv, pairwise, names = compute_cv_and_pairwise(alphas)

        flush_print(f"  --- p = {p:.4f} ---")
        flush_print(f"  {'Decoder':14s}  {'alpha':>10s}  {'alpha_err':>10s}  {'beta':>10s}")
        flush_print(f"  {'-'*14}  {'-'*10}  {'-'*10}  {'-'*10}")

        for name in names:
            a = alphas[name]
            flush_print(f"  {name:14s}  {a['alpha']:10.6f}  {a['alpha_err']:10.6f}  "
                        f"{a['beta']:10.6f}")

        if cv is not None:
            flush_print(f"\n  CV = {cv:.4f}")
            flush_print(f"  Pairwise comparisons:")
            for (n1, n2, a1, a2, diff, cerr, sig) in pairwise:
                flush_print(f"    {n1:14s} vs {n2:14s}: "
                            f"diff={diff:.6f}, {sig:.2f} sigma")
        flush_print()


def print_summary(results, alphas_by_p):
    """Print comprehensive summary with verdict."""
    flush_print()
    flush_print("=" * 75)
    flush_print("  COMPREHENSIVE SUMMARY")
    flush_print("=" * 75)
    flush_print()

    # 1. Decoder hierarchy check
    flush_print("  1. DECODER HIERARCHY (sanity check)")
    flush_print("  " + "-" * 60)
    for d in DISTANCES:
        rates = {}
        p_check = 0.005
        for name in DECODER_NAMES_ORDERED:
            if d in results[name] and p_check in results[name][d]:
                rates[name] = results[name][d][p_check]['base_rate']
        if rates:
            sorted_dec = sorted(rates.items(), key=lambda x: x[1])
            ranking = ' < '.join([f'{n}({r:.6f})' for n, r in sorted_dec])
            flush_print(f"    d={d:2d}, p={p_check}: {ranking}")
    flush_print()

    # 2. Alpha values at each p
    flush_print("  2. ALPHA VALUES AND CV AT EACH P")
    flush_print("  " + "-" * 60)
    flush_print(f"  {'p':>8s}  {'alpha_mean':>12s}  {'alpha_std':>12s}  {'CV':>8s}  {'Status':>12s}")
    flush_print(f"  {'-'*8}  {'-'*12}  {'-'*12}  {'-'*8}  {'-'*12}")

    all_cvs = []
    for p in P_VALUES:
        alphas = alphas_by_p[p]
        cv, pairwise, names = compute_cv_and_pairwise(alphas)
        if cv is not None:
            alpha_vals = [alphas[n]['alpha'] for n in names]
            mean_a = np.mean(alpha_vals)
            std_a = np.std(alpha_vals)
            if cv < 0.15:
                status = "SUPPORTED"
            elif cv < 0.30:
                status = "PARTIAL"
            else:
                status = "NOT SUPPORTED"
            flush_print(f"  {p:8.4f}  {mean_a:12.6f}  {std_a:12.6f}  {cv:8.4f}  {status:>12s}")
            all_cvs.append(cv)
    flush_print()

    # 3. Pairwise sigma-differences summary
    flush_print("  3. PAIRWISE SIGMA-DIFFERENCES (should be < 2.0 sigma)")
    flush_print("  " + "-" * 60)
    for p in P_VALUES:
        alphas = alphas_by_p[p]
        cv, pairwise, names = compute_cv_and_pairwise(alphas)
        if pairwise:
            max_sig = max(pw[6] for pw in pairwise)
            pairs_str = ', '.join([f'{pw[0][:6]} vs {pw[1][:6]}: {pw[6]:.2f}s' for pw in pairwise])
            flush_print(f"    p={p:.4f}: max={max_sig:.2f} sigma  [{pairs_str}]")
    flush_print()

    # 4. Overall verdict
    flush_print("  4. OVERALL VERDICT")
    flush_print("  " + "-" * 60)
    if all_cvs:
        mean_cv = np.mean(all_cvs)
        flush_print(f"    Mean CV across p values: {mean_cv:.4f}")
        flush_print(f"    Distances tested: {DISTANCES}")
        flush_print(f"    Decoders tested: {DECODER_NAMES_ORDERED}")
        flush_print()

        if mean_cv < 0.15:
            flush_print(f"    *** VERDICT: SUPPORTED (mean CV = {mean_cv:.4f} < 0.15) ***")
            flush_print(f"    The post-selection improvement slope alpha is consistent")
            flush_print(f"    across {len(DECODER_NAMES_ORDERED)} decoders at d=3,...,11.")
            flush_print(f"    Alpha is determined by the physical noise channel,")
            flush_print(f"    not by the decoder choice.")
        elif mean_cv < 0.30:
            flush_print(f"    *** VERDICT: PARTIAL (mean CV = {mean_cv:.4f}, between 0.15 and 0.30) ***")
            flush_print(f"    Alpha shows some variation but remains moderately consistent.")
            flush_print(f"    Extended distances may introduce finite-size effects.")
        else:
            flush_print(f"    *** VERDICT: NOT SUPPORTED (mean CV = {mean_cv:.4f} > 0.30) ***")
            flush_print(f"    Significant variation in alpha observed across decoders.")
    flush_print()

    # 5. Comparison with old (d=3,5,7) results
    flush_print("  5. OLD (d=3,5,7) vs NEW (d=3-11) COMPARISON")
    flush_print("  " + "-" * 60)
    # Re-extract alphas using only d=3,5,7
    for p in P_VALUES:
        # Old: only d=3,5,7
        alphas_old = {}
        for name in DECODER_NAMES_ORDERED:
            d_vals_old = []
            lnR_vals_old = []
            lnR_errs_old = []
            for d in [3, 5, 7]:
                if d in results[name] and p in results[name][d]:
                    data = results[name][d][p]
                    if not np.isnan(data['lnR']):
                        d_vals_old.append(d)
                        lnR_vals_old.append(data['lnR'])
                        lnR_errs_old.append(data['lnR_err'])
            if len(d_vals_old) >= 2:
                d_arr = np.array(d_vals_old, dtype=float)
                lnR_arr = np.array(lnR_vals_old)
                lnR_err_arr = np.array(lnR_errs_old)
                w = 1.0 / (lnR_err_arr**2 + 1e-10)
                W = w.sum()
                Wd = (w * d_arr).sum()
                WlnR = (w * lnR_arr).sum()
                Wdd = (w * d_arr**2).sum()
                WdlnR = (w * d_arr * lnR_arr).sum()
                det_val = W * Wdd - Wd**2
                if abs(det_val) > 1e-20:
                    alpha = (W * WdlnR - Wd * WlnR) / det_val
                    alpha_err = np.sqrt(W / det_val)
                    alphas_old[name] = {'alpha': alpha, 'alpha_err': alpha_err}

        # New: all distances
        alphas_new = alphas_by_p[p]

        cv_old, _, names_old = compute_cv_and_pairwise(alphas_old)
        cv_new, _, names_new = compute_cv_and_pairwise(alphas_new)

        cv_old_str = f"{cv_old:.4f}" if cv_old is not None else "N/A"
        cv_new_str = f"{cv_new:.4f}" if cv_new is not None else "N/A"
        flush_print(f"    p={p:.4f}: CV(d=3,5,7) = {cv_old_str}  ->  CV(d=3-11) = {cv_new_str}")

    flush_print()


# ============================================================================
# Visualization
# ============================================================================

def generate_figures(results, alphas_by_p):
    """Generate 4-panel figure (16x10 inches)."""
    flush_print("=" * 75)
    flush_print("  GENERATING FIGURES")
    flush_print("=" * 75)
    flush_print()

    decoder_colors = {
        'MWPM':       '#1565C0',
        'Unwt-MWPM':  '#E65100',
        'Union-Find': '#2E7D32',
    }
    decoder_markers = {
        'MWPM':       'o',
        'Unwt-MWPM':  's',
        'Union-Find': 'D',
    }

    eps_pct = int(PS_EPSILON * 100)
    eps_label = r'$\varepsilon=' + str(eps_pct) + r'\%$'

    fig = plt.figure(figsize=(16, 10))
    gs = GridSpec(2, 2, figure=fig, wspace=0.30, hspace=0.42,
                  left=0.07, right=0.97, bottom=0.08, top=0.88)

    # ================================================================
    # Panel (a): ln R vs d for each decoder at p=0.005
    # ================================================================
    ax_a = fig.add_subplot(gs[0, 0])
    p_plot = 0.005
    alphas_at_p = alphas_by_p[p_plot]

    for name in DECODER_NAMES_ORDERED:
        if name not in alphas_at_p:
            continue
        a = alphas_at_p[name]
        color = decoder_colors[name]
        marker = decoder_markers[name]

        # Data points for this decoder at this p only
        d_pts = []
        lnR_pts = []
        lnR_err_pts = []
        for d in DISTANCES:
            if d in results[name] and p_plot in results[name][d]:
                data = results[name][d][p_plot]
                if not np.isnan(data['lnR']):
                    d_pts.append(d)
                    lnR_pts.append(data['lnR'])
                    lnR_err_pts.append(data['lnR_err'])

        if d_pts:
            ax_a.errorbar(
                d_pts, lnR_pts, yerr=lnR_err_pts,
                color=color, marker=marker, markersize=7,
                linestyle='none', linewidth=1.5, capsize=3,
                label=name, zorder=5,
            )

        # Plot fit line
        d_fit = np.linspace(min(DISTANCES) - 0.5, max(DISTANCES) + 0.5, 100)
        lnR_fit = a['alpha'] * d_fit + a['beta']
        ax_a.plot(d_fit, lnR_fit, color=color, linestyle='--', linewidth=1.2,
                  alpha=0.6, zorder=3)

    ax_a.set_xlabel('Code distance $d$', fontsize=12)
    ax_a.set_ylabel('$\\ln R$', fontsize=12)
    ax_a.set_title(f'(a) $\\ln R$ vs $d$  [$p={p_plot}$, {eps_label}]',
                   fontsize=13, fontweight='bold')
    ax_a.legend(fontsize=9, loc='best', framealpha=0.9)
    ax_a.grid(True, alpha=0.3)
    ax_a.set_xticks(DISTANCES)

    # Annotate with slopes
    cv_a, _, _ = compute_cv_and_pairwise(alphas_at_p)
    if cv_a is not None:
        alpha_vals_a = [alphas_at_p[n]['alpha'] for n in DECODER_NAMES_ORDERED if n in alphas_at_p]
        alpha_mean_a = np.mean(alpha_vals_a)
        alpha_std_a = np.std(alpha_vals_a)
        verdict_a = "SUPPORTED" if cv_a < 0.15 else ("PARTIAL" if cv_a < 0.30 else "NOT SUPP.")
        ax_a.annotate(
            f'Slopes ($\\alpha$):\n'
            f'mean = {alpha_mean_a:.4f}\n'
            f'std = {alpha_std_a:.4f}\n'
            f'CV = {cv_a:.3f} ({verdict_a})',
            xy=(0.03, 0.97), xycoords='axes fraction', fontsize=8,
            ha='left', va='top',
            bbox=dict(boxstyle='round,pad=0.3',
                      fc='#E8F5E9' if cv_a < 0.15 else ('#FFF3E0' if cv_a < 0.30 else '#FFEBEE'),
                      ec='green' if cv_a < 0.15 else ('orange' if cv_a < 0.30 else 'red'),
                      alpha=0.9))

    # ================================================================
    # Panel (b): alpha values with error bars for each decoder (at each p)
    # ================================================================
    ax_b = fig.add_subplot(gs[0, 1])

    x_offset_map = {'MWPM': -0.15, 'Unwt-MWPM': 0.0, 'Union-Find': 0.15}
    for name in DECODER_NAMES_ORDERED:
        color = decoder_colors[name]
        marker = decoder_markers[name]
        p_pts = []
        alpha_pts = []
        alpha_err_pts = []
        for p in P_VALUES:
            if name in alphas_by_p[p]:
                a = alphas_by_p[p][name]
                p_pts.append(p)
                alpha_pts.append(a['alpha'])
                alpha_err_pts.append(a['alpha_err'])

        if p_pts:
            x_vals = np.array(p_pts) + x_offset_map.get(name, 0) * 0.0003
            ax_b.errorbar(
                x_vals, alpha_pts, yerr=alpha_err_pts,
                color=color, marker=marker, markersize=7,
                linestyle='-', linewidth=1.2, capsize=4,
                label=name, zorder=5,
            )

    ax_b.set_xlabel('Physical error rate $p$', fontsize=12)
    ax_b.set_ylabel(r'$\alpha$ (slope of $\ln R$ vs $d$)', fontsize=12)
    ax_b.set_title(f'(b) $\\alpha$ vs $p$ for each decoder  [{eps_label}]',
                   fontsize=13, fontweight='bold')
    ax_b.legend(fontsize=9, loc='best', framealpha=0.9)
    ax_b.grid(True, alpha=0.3)

    ax_b.annotate(
        'If $\\alpha$ is decoder-independent,\ncurves should overlap',
        xy=(0.03, 0.03), xycoords='axes fraction', fontsize=8,
        ha='left', va='bottom',
        bbox=dict(boxstyle='round,pad=0.3', fc='lightyellow', ec='orange', alpha=0.9))

    # ================================================================
    # Panel (c): CV vs p
    # ================================================================
    ax_c = fig.add_subplot(gs[1, 0])

    p_cv_list = []
    cv_list = []
    for p in P_VALUES:
        cv, _, _ = compute_cv_and_pairwise(alphas_by_p[p])
        if cv is not None:
            p_cv_list.append(p)
            cv_list.append(cv)

    if p_cv_list:
        ax_c.plot(p_cv_list, cv_list, 'ko-', markersize=8, linewidth=2, zorder=5)
        ax_c.axhline(y=0.15, color='green', linestyle='--', linewidth=1.5,
                     label='SUPPORTED threshold (0.15)', zorder=3)
        ax_c.axhline(y=0.30, color='orange', linestyle='--', linewidth=1.5,
                     label='PARTIAL threshold (0.30)', zorder=3)
        ax_c.fill_between([min(P_VALUES)*0.8, max(P_VALUES)*1.2], 0, 0.15,
                          color='green', alpha=0.08, zorder=1)
        ax_c.fill_between([min(P_VALUES)*0.8, max(P_VALUES)*1.2], 0.15, 0.30,
                          color='orange', alpha=0.08, zorder=1)
        ax_c.fill_between([min(P_VALUES)*0.8, max(P_VALUES)*1.2], 0.30, 1.0,
                          color='red', alpha=0.08, zorder=1)

        for i, (px, cvx) in enumerate(zip(p_cv_list, cv_list)):
            ax_c.annotate(f'{cvx:.3f}', (px, cvx), textcoords="offset points",
                          xytext=(0, 12), ha='center', fontsize=9, fontweight='bold')

    ax_c.set_xlabel('Physical error rate $p$', fontsize=12)
    ax_c.set_ylabel('CV (coefficient of variation)', fontsize=12)
    ax_c.set_title('(c) CV vs $p$  [should be small for support]',
                   fontsize=13, fontweight='bold')
    ax_c.legend(fontsize=8, loc='best', framealpha=0.9)
    ax_c.grid(True, alpha=0.3)
    ax_c.set_ylim(bottom=0)

    # ================================================================
    # Panel (d): Comparison table: old vs new CV
    # ================================================================
    ax_d = fig.add_subplot(gs[1, 1])
    ax_d.axis('off')

    # Build comparison data
    table_data = []
    for p in P_VALUES:
        # Old (d=3,5,7)
        alphas_old = {}
        for name in DECODER_NAMES_ORDERED:
            d_vals_old = []
            lnR_vals_old = []
            lnR_errs_old = []
            for d in [3, 5, 7]:
                if d in results[name] and p in results[name][d]:
                    data = results[name][d][p]
                    if not np.isnan(data['lnR']):
                        d_vals_old.append(d)
                        lnR_vals_old.append(data['lnR'])
                        lnR_errs_old.append(data['lnR_err'])
            if len(d_vals_old) >= 2:
                d_arr = np.array(d_vals_old, dtype=float)
                lnR_arr = np.array(lnR_vals_old)
                lnR_err_arr = np.array(lnR_errs_old)
                w = 1.0 / (lnR_err_arr**2 + 1e-10)
                W = w.sum()
                Wd = (w * d_arr).sum()
                WlnR = (w * lnR_arr).sum()
                Wdd = (w * d_arr**2).sum()
                WdlnR = (w * d_arr * lnR_arr).sum()
                det_val = W * Wdd - Wd**2
                if abs(det_val) > 1e-20:
                    alpha = (W * WdlnR - Wd * WlnR) / det_val
                    alpha_err = np.sqrt(W / det_val)
                    alphas_old[name] = {'alpha': alpha, 'alpha_err': alpha_err}

        cv_old, _, _ = compute_cv_and_pairwise(alphas_old)
        cv_new, _, _ = compute_cv_and_pairwise(alphas_by_p[p])

        cv_old_str = f"{cv_old:.4f}" if cv_old is not None else "N/A"
        cv_new_str = f"{cv_new:.4f}" if cv_new is not None else "N/A"

        if cv_new is not None:
            if cv_new < 0.15:
                verdict = "SUPPORTED"
            elif cv_new < 0.30:
                verdict = "PARTIAL"
            else:
                verdict = "NOT SUPP."
        else:
            verdict = "N/A"

        table_data.append([f"{p:.4f}", cv_old_str, cv_new_str, verdict])

    # Add overall row
    all_cv_old = []
    all_cv_new = []
    for p in P_VALUES:
        alphas_old_temp = {}
        for name in DECODER_NAMES_ORDERED:
            d_vals_old = []
            lnR_vals_old = []
            lnR_errs_old = []
            for d in [3, 5, 7]:
                if d in results[name] and p in results[name][d]:
                    data = results[name][d][p]
                    if not np.isnan(data['lnR']):
                        d_vals_old.append(d)
                        lnR_vals_old.append(data['lnR'])
                        lnR_errs_old.append(data['lnR_err'])
            if len(d_vals_old) >= 2:
                d_arr = np.array(d_vals_old, dtype=float)
                lnR_arr = np.array(lnR_vals_old)
                lnR_err_arr = np.array(lnR_errs_old)
                w = 1.0 / (lnR_err_arr**2 + 1e-10)
                W = w.sum()
                Wd = (w * d_arr).sum()
                WlnR = (w * lnR_arr).sum()
                Wdd = (w * d_arr**2).sum()
                WdlnR = (w * d_arr * lnR_arr).sum()
                det_val = W * Wdd - Wd**2
                if abs(det_val) > 1e-20:
                    alpha = (W * WdlnR - Wd * WlnR) / det_val
                    alpha_err = np.sqrt(W / det_val)
                    alphas_old_temp[name] = {'alpha': alpha, 'alpha_err': alpha_err}
        cv_o, _, _ = compute_cv_and_pairwise(alphas_old_temp)
        cv_n, _, _ = compute_cv_and_pairwise(alphas_by_p[p])
        if cv_o is not None:
            all_cv_old.append(cv_o)
        if cv_n is not None:
            all_cv_new.append(cv_n)

    if all_cv_old and all_cv_new:
        mean_old = np.mean(all_cv_old)
        mean_new = np.mean(all_cv_new)
        if mean_new < 0.15:
            v = "SUPPORTED"
        elif mean_new < 0.30:
            v = "PARTIAL"
        else:
            v = "NOT SUPP."
        table_data.append(["MEAN", f"{mean_old:.4f}", f"{mean_new:.4f}", v])

    col_labels = ['p', 'CV (d=3,5,7)', 'CV (d=3-11)', 'Verdict']
    table = ax_d.table(
        cellText=table_data,
        colLabels=col_labels,
        cellLoc='center',
        loc='center',
        colWidths=[0.18, 0.25, 0.25, 0.25],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(11)
    table.scale(1.0, 1.8)

    # Color the verdict column
    for i, row in enumerate(table_data):
        cell = table[i + 1, 3]  # +1 for header row
        if row[3] == "SUPPORTED":
            cell.set_facecolor('#E8F5E9')
        elif row[3] == "PARTIAL":
            cell.set_facecolor('#FFF3E0')
        elif row[3] == "NOT SUPP.":
            cell.set_facecolor('#FFEBEE')
        # Bold the MEAN row
        if row[0] == "MEAN":
            for col_idx in range(4):
                table[i + 1, col_idx].set_text_props(fontweight='bold')
                table[i + 1, col_idx].set_facecolor('#E3F2FD')
            # Override verdict color for MEAN
            if row[3] == "SUPPORTED":
                table[i + 1, 3].set_facecolor('#C8E6C9')
            elif row[3] == "PARTIAL":
                table[i + 1, 3].set_facecolor('#FFE0B2')

    # Style header
    for col_idx in range(4):
        table[0, col_idx].set_facecolor('#1565C0')
        table[0, col_idx].set_text_props(color='white', fontweight='bold')

    ax_d.set_title('(d) CV Comparison: Old vs New',
                   fontsize=13, fontweight='bold', pad=20)

    # Supertitle
    fig.suptitle(
        'Extended Decoder $\\alpha$-Independence Test ($d=3,5,7,9,11$)\n'
        '$\\ln R(\\varepsilon, d) = \\alpha(p) \\cdot d + \\beta(\\varepsilon)$: '
        '$\\alpha$ should be decoder-independent',
        fontsize=14, fontweight='bold', y=0.97)

    fig.savefig(FIG_PATH_PNG, dpi=200, bbox_inches='tight')
    fig.savefig(FIG_PATH_PDF, bbox_inches='tight')
    plt.close(fig)

    flush_print(f"  Saved: {FIG_PATH_PNG}")
    flush_print(f"  Saved: {FIG_PATH_PDF}")
    flush_print()


# ============================================================================
# Main
# ============================================================================

def main():
    sys.stdout.reconfigure(line_buffering=True)

    t_start = time.time()

    flush_print()
    flush_print("*" * 75)
    flush_print("*  EXTENDED DECODER ALPHA-INDEPENDENCE TEST (Large d)                 *")
    flush_print("*  Testing: ln R(eps,d) = alpha(p)*d + beta(eps)                      *")
    flush_print("*  Prediction: alpha is INDEPENDENT of decoder choice                 *")
    flush_print("*  Distances: d = 3, 5, 7, 9, 11                                     *")
    flush_print("*  Decoders: MWPM, Unwt-MWPM, Union-Find                             *")
    flush_print("*" * 75)
    flush_print()
    flush_print(f"  Distances: {DISTANCES}")
    flush_print(f"  P values: {P_VALUES}")
    flush_print(f"  Post-selection epsilon: {PS_EPSILON}")
    flush_print(f"  Shots (fast decoders): {SHOTS_FAST}")
    flush_print(f"  Shots (Union-Find):    {SHOTS_UF}")
    flush_print()

    # Run simulation
    results = run_simulation()

    # Extract alphas at each p
    alphas_by_p = {}
    for p in P_VALUES:
        alphas_by_p[p] = extract_alphas(results, p_fixed=p)

    # Detailed analysis
    print_analysis(results, alphas_by_p)

    # Generate figures
    generate_figures(results, alphas_by_p)

    # Summary with verdict
    print_summary(results, alphas_by_p)

    t_total = time.time() - t_start
    flush_print(f"  Total execution time: {t_total:.1f} seconds ({t_total/60:.1f} minutes)")
    flush_print()


if __name__ == '__main__':
    main()
