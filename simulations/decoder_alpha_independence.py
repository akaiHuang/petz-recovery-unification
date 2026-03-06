#!/usr/bin/env python3
"""
Decoder Alpha-Independence Test for Prediction 1
=================================================

Tests the prediction that post-selection improvement ln R(eps,d) = alpha(p)*d + beta(eps),
where alpha is INDEPENDENT of decoder choice.

This is the key test: if alpha is decoder-independent, it means the post-selection
improvement scales with code distance in a universal way determined by the physical
noise channel (entropy production Sigma), not by the recovery strategy (decoder).

Decoders tested:
  1. MWPM (weighted)      -- pymatching, industry standard
  2. Unweighted MWPM      -- pymatching with uniform weights
  3. Greedy NN             -- nearest-neighbor greedy matching
  4. Union-Find            -- ldpc package, cluster-growth decoder
  5. BP+OSD                -- ldpc package, belief propagation + ordered statistics
  6. Lookup Table          -- exhaustive single-error enumeration (d=3 only)

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
from ldpc import UnionFindDecoder as LdpcUF, BpOsdDecoder as LdpcBpOsd

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

DISTANCES = [3, 5, 7]
P_VALUES = [0.003, 0.005, 0.007, 0.01, 0.015]

# Post-selection parameters
PS_FRACTIONS = [0.10, 0.20, 0.30]  # epsilon values
PS_EPSILON_DEFAULT = 0.20           # default for alpha extraction

# Monte Carlo shots (for fast decoders: MWPM, Unwt-MWPM)
SHOTS_PER_POINT = 100_000  # per (d, p) combination

# Output paths
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FIG_PATH_PNG = os.path.join(SCRIPT_DIR, 'fig_decoder_alpha_independence.png')
FIG_PATH_PDF = os.path.join(SCRIPT_DIR, 'fig_decoder_alpha_independence.pdf')


# ============================================================================
# Helper: DEM -> Check Matrix extraction
# ============================================================================

def dem_to_check_matrix(dem: stim.DetectorErrorModel, n_det: int, n_obs: int = 1):
    """
    Extract parity check matrix, observable matrix, and error probabilities
    from a stim DetectorErrorModel.

    Returns:
        H: (n_det, n_errors) check matrix
        obs_matrix: (n_obs, n_errors) observable matrix
        probs: (n_errors,) error probabilities
        llrs: (n_errors,) log-likelihood ratios
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
        probs[j] = max(p, 1e-15)  # avoid log(0)
        for d in dets:
            H[d, j] = 1
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


class GreedyNNDecoder:
    """Greedy Nearest-Neighbor decoder using detector coordinates."""
    def __init__(self, dem, circuit=None, **kwargs):
        self.name = "Greedy-NN"
        assert circuit is not None

        det_coords = circuit.get_detector_coordinates()
        self.num_detectors = circuit.num_detectors
        self.coords = np.zeros((self.num_detectors, 3))
        for det_id, coord in det_coords.items():
            c = list(coord)
            self.coords[det_id] = c[:3] if len(c) >= 3 else c + [0] * (3 - len(c))

        self.boundary_detectors = set()
        for inst in dem.flattened():
            if inst.type == 'error':
                targets = inst.targets_copy()
                has_obs = any(t.is_logical_observable_id() for t in targets)
                if has_obs:
                    for t in targets:
                        if t.is_relative_detector_id():
                            self.boundary_detectors.add(t.val)

    def decode_batch(self, detections, obs_flat):
        num_shots = detections.shape[0]
        errors = np.zeros(num_shots, dtype=bool)
        for i in range(num_shots):
            pred = self._decode_single(detections[i])
            errors[i] = (pred != obs_flat[i])
        return errors

    def _decode_single(self, detection_events):
        triggered = np.where(detection_events)[0]
        if len(triggered) == 0:
            return 0

        observable_flips = 0
        active_coords = self.coords[triggered].copy()
        active_ids = list(triggered)
        is_boundary = [d in self.boundary_detectors for d in active_ids]
        paired = [False] * len(active_ids)

        while True:
            unpaired = [j for j in range(len(active_ids)) if not paired[j]]
            if len(unpaired) == 0:
                break
            if len(unpaired) == 1:
                j = unpaired[0]
                paired[j] = True
                if is_boundary[j]:
                    observable_flips += 1
                break

            best_dist = float('inf')
            best_pair = (unpaired[0], unpaired[1])
            for ai in range(len(unpaired)):
                for bi in range(ai + 1, len(unpaired)):
                    ja, jb = unpaired[ai], unpaired[bi]
                    dist = np.linalg.norm(active_coords[ja] - active_coords[jb])
                    if dist < best_dist:
                        best_dist = dist
                        best_pair = (ja, jb)

            paired[best_pair[0]] = True
            paired[best_pair[1]] = True

        return observable_flips % 2


class UnionFindDecoder:
    """Union-Find decoder via ldpc package."""
    def __init__(self, dem, circuit=None, **kwargs):
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


class BpOsdDecoder:
    """Belief Propagation + Ordered Statistics Decoding via ldpc package."""
    def __init__(self, dem, circuit=None, **kwargs):
        self.name = "BP+OSD"
        n_det = kwargs.get('n_det', circuit.num_detectors if circuit else None)
        n_obs = kwargs.get('n_obs', 1)
        H, obs_matrix, probs, llrs = dem_to_check_matrix(dem, n_det, n_obs)
        self.H = H
        self.obs_matrix = obs_matrix
        self.bposd = LdpcBpOsd(
            H,
            channel_probs=list(probs),
            max_iter=50,
            bp_method='product_sum',
            osd_method='osd_cs',
            osd_order=4,
        )

    def decode_batch(self, detections, obs_flat):
        num_shots = detections.shape[0]
        errors = np.zeros(num_shots, dtype=bool)
        for i in range(num_shots):
            syndrome = detections[i].astype(np.uint8)
            estimate = self.bposd.decode(syndrome)
            pred_obs = (self.obs_matrix @ estimate) % 2
            errors[i] = (pred_obs[0] != obs_flat[i])
        return errors


class LookupTableDecoder:
    """
    Lookup table decoder using single-error enumeration from DEM.
    Only practical for small codes (d=3).

    For each possible syndrome pattern arising from a single error mechanism,
    stores the most likely observable flip. For unseen syndromes, defaults
    to no flip (which is a poor approximation for high-weight syndromes).
    """
    def __init__(self, dem, circuit=None, **kwargs):
        self.name = "Lookup"
        n_det = kwargs.get('n_det', circuit.num_detectors if circuit else None)

        # Build lookup: syndrome -> (prob_no_flip, prob_flip)
        self.lookup = {}
        for inst in dem.flattened():
            if inst.type == 'error':
                targets = inst.targets_copy()
                p = inst.args_copy()[0]
                dets = tuple(sorted([t.val for t in targets if t.is_relative_detector_id()]))
                obs = [t.val for t in targets if t.is_logical_observable_id()]
                obs_flip = 1 if len(obs) > 0 else 0

                if dets not in self.lookup:
                    self.lookup[dets] = [0.0, 0.0]
                self.lookup[dets][obs_flip] += p

    def decode_batch(self, detections, obs_flat):
        num_shots = detections.shape[0]
        errors = np.zeros(num_shots, dtype=bool)
        for i in range(num_shots):
            triggered = tuple(sorted(np.where(detections[i])[0]))
            probs_obs = self.lookup.get(triggered, [1.0, 0.0])
            prediction = 1 if probs_obs[1] > probs_obs[0] else 0
            errors[i] = (prediction != obs_flat[i])
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

    # Relative uncertainties (standard errors)
    se_base = np.sqrt(base_rate * (1 - base_rate) / base_n) if base_n > 0 else 0
    se_ps = np.sqrt(ps_rate * (1 - ps_rate) / ps_n) if ps_n > 0 else 0

    # Delta method for ratio
    rel_var = (se_base / base_rate)**2 + (se_ps / ps_rate)**2
    R_se = R * np.sqrt(rel_var)
    return R, 1.96 * R_se  # 95% CI


# ============================================================================
# Main Simulation
# ============================================================================

def get_decoders_for_distance(d, dem, circuit):
    """
    Return list of decoder instances appropriate for this distance.
    Lookup table is only used for d=3.
    """
    n_det = circuit.num_detectors
    n_obs = circuit.num_observables
    kwargs = {'n_det': n_det, 'n_obs': n_obs}

    decoders = [
        MWPMDecoder(dem, circuit, **kwargs),
        UnweightedMWPMDecoder(dem, circuit, **kwargs),
        UnionFindDecoder(dem, circuit, **kwargs),
        BpOsdDecoder(dem, circuit, **kwargs),
        GreedyNNDecoder(dem, circuit, **kwargs),
    ]
    if d == 3:
        decoders.append(LookupTableDecoder(dem, circuit, **kwargs))

    return decoders


def run_simulation():
    """
    Run the full alpha-independence test.

    For each (d, p):
      1. Generate circuit and sample detection events (same samples for all decoders)
      2. Decode with each decoder
      3. Compute logical error rate with and without post-selection
      4. Compute ln R for each (decoder, d, p, epsilon) combination
    """
    print()
    print("=" * 75)
    print("  DECODER ALPHA-INDEPENDENCE TEST")
    print("  Prediction 1: ln R(eps,d) = alpha(p)*d + beta(eps)")
    print("  alpha should be INDEPENDENT of decoder choice")
    print("=" * 75)
    print()

    # results[decoder_name][d][p] = {
    #   'base_rate': float, 'base_ci': float, 'base_n': int,
    #   'ps': {eps: {'rate': float, 'ci': float, 'n': int, 'R': float, 'R_ci': float}}
    # }
    results = {}
    all_decoder_names = set()

    total_configs = len(DISTANCES) * len(P_VALUES)
    config_idx = 0

    for d in DISTANCES:
        for p in P_VALUES:
            config_idx += 1
            shots = SHOTS_PER_POINT

            print(f"  [{config_idx}/{total_configs}] d={d}, p={p:.4f}, shots={shots:,}")

            # Generate circuit and sample ONCE
            circuit = generate_circuit(d, p)
            dem = circuit.detector_error_model(decompose_errors=True)
            sampler = circuit.compile_detector_sampler()

            t0 = time.time()
            det, obs = sampler.sample(shots, separate_observables=True)
            obs_flat = obs.flatten()
            weights = det.sum(axis=1).astype(float)

            # Get decoders
            decoders = get_decoders_for_distance(d, dem, circuit)

            # For slow decoders (per-shot decoding), limit shots significantly
            # BP+OSD and Union-Find decode one shot at a time and are ~100x slower
            SLOW_DECODERS = {'Greedy-NN', 'Union-Find', 'BP+OSD', 'Lookup'}
            # Scale slow decoder shots by distance: d=3 is fast, d=7 is very slow
            slow_scale = {3: 10_000, 5: 3_000, 7: 1_000}
            SLOW_MAX_SHOTS = slow_scale.get(d, 1_000)

            for decoder in decoders:
                name = decoder.name
                all_decoder_names.add(name)

                if name not in results:
                    results[name] = {}
                if d not in results[name]:
                    results[name][d] = {}

                # Determine shot count for this decoder
                if name in SLOW_DECODERS:
                    n_use = SLOW_MAX_SHOTS
                else:
                    n_use = shots

                det_use = det[:n_use]
                obs_use = obs_flat[:n_use]
                w_use = weights[:n_use]

                t_dec = time.time()
                is_error = decoder.decode_batch(det_use, obs_use)
                t_dec = time.time() - t_dec

                base_rate, base_ci = compute_wilson_ci(is_error)

                # Post-selection
                ps_data = {}
                for eps in PS_FRACTIONS:
                    cutoff = np.quantile(w_use, 1.0 - eps)
                    mask = w_use <= cutoff
                    n_kept = mask.sum()

                    if n_kept > 0 and is_error[mask].size > 0:
                        ps_rate, ps_ci = compute_wilson_ci(is_error[mask])
                        R, R_ci = compute_R_with_error(base_rate, n_use, ps_rate, n_kept)
                    else:
                        ps_rate, ps_ci = np.nan, np.nan
                        R, R_ci = np.nan, np.nan

                    ps_data[eps] = {
                        'rate': ps_rate, 'ci': ps_ci, 'n': n_kept,
                        'R': R, 'R_ci': R_ci,
                    }

                results[name][d][p] = {
                    'base_rate': base_rate,
                    'base_ci': base_ci,
                    'base_n': n_use,
                    'ps': ps_data,
                }

                lnR_20 = np.log(ps_data[0.20]['R']) if ps_data[0.20]['R'] > 0 and not np.isnan(ps_data[0.20]['R']) else np.nan
                print(f"    {name:14s}: p_L={base_rate:.6f}, "
                      f"R(20%)={ps_data[0.20]['R']:.3f}, "
                      f"ln R={lnR_20:.4f}, "
                      f"({t_dec:.1f}s, {n_use:,} shots)")

            t_total = time.time() - t0
            print(f"    Total: {t_total:.1f}s")
            print()

    return results, sorted(all_decoder_names)


# ============================================================================
# Alpha Extraction
# ============================================================================

def extract_alphas(results, decoder_names, p_fixed=0.005, eps_fixed=0.20):
    """
    For each decoder, compute ln R vs d at fixed (p, eps) and do linear fit.
    ln R(d) = alpha * d + const
    Returns dict: decoder_name -> (alpha, alpha_err, beta, beta_err, d_values, lnR_values)
    """
    alphas = {}

    for name in decoder_names:
        d_values = []
        lnR_values = []
        lnR_errors = []

        for d in DISTANCES:
            if d not in results[name] or p_fixed not in results[name][d]:
                continue
            data = results[name][d][p_fixed]
            ps = data['ps'].get(eps_fixed, None)
            if ps is None or np.isnan(ps['R']) or ps['R'] <= 0:
                continue

            lnR = np.log(ps['R'])
            # Propagate R uncertainty to ln R
            if ps['R_ci'] > 0 and ps['R'] > 0:
                lnR_err = ps['R_ci'] / ps['R']  # delta method
            else:
                lnR_err = 0.1  # fallback

            d_values.append(d)
            lnR_values.append(lnR)
            lnR_errors.append(lnR_err)

        if len(d_values) >= 2:
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
                # Unweighted fit
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


def extract_alphas_all_p(results, decoder_names, eps_fixed=0.20):
    """Extract alpha for each decoder at each p value."""
    alphas_by_p = {}
    for p in P_VALUES:
        alphas_by_p[p] = extract_alphas(results, decoder_names, p_fixed=p, eps_fixed=eps_fixed)
    return alphas_by_p


# ============================================================================
# Analysis and Printing
# ============================================================================

def print_analysis(results, decoder_names, alphas):
    """Print detailed analysis of alpha-independence."""
    print()
    print("=" * 75)
    print("  ALPHA-INDEPENDENCE ANALYSIS")
    print("=" * 75)
    print()

    # Print alpha values for each decoder
    print(f"  {'Decoder':14s}  {'alpha':>10s}  {'alpha_err':>10s}  {'beta':>10s}  {'beta_err':>10s}")
    print(f"  {'-'*14}  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*10}")

    alpha_values = []
    alpha_errors = []
    decoder_names_with_alpha = []

    for name in decoder_names:
        if name in alphas:
            a = alphas[name]
            print(f"  {name:14s}  {a['alpha']:10.6f}  {a['alpha_err']:10.6f}  "
                  f"{a['beta']:10.6f}  {a['beta_err']:10.6f}")
            alpha_values.append(a['alpha'])
            alpha_errors.append(a['alpha_err'])
            decoder_names_with_alpha.append(name)

    print()

    if len(alpha_values) >= 2:
        alpha_arr = np.array(alpha_values)
        alpha_mean = alpha_arr.mean()
        alpha_std = alpha_arr.std()
        CV = alpha_std / abs(alpha_mean) if abs(alpha_mean) > 1e-10 else float('inf')

        print(f"  Alpha statistics:")
        print(f"    Mean:   {alpha_mean:.6f}")
        print(f"    Std:    {alpha_std:.6f}")
        print(f"    CV:     {CV:.4f}  (coefficient of variation = std/|mean|)")
        print()

        if CV < 0.10:
            print(f"    RESULT: CV = {CV:.4f} < 0.10")
            print(f"    ==> ALPHA-INDEPENDENCE SUPPORTED")
            print(f"    The slope alpha is consistent across all {len(alpha_values)} decoders.")
        elif CV < 0.20:
            print(f"    RESULT: CV = {CV:.4f} (between 0.10 and 0.20)")
            print(f"    ==> WEAK SUPPORT for alpha-independence")
            print(f"    Some variation observed but still relatively small.")
        else:
            print(f"    RESULT: CV = {CV:.4f} > 0.20")
            print(f"    ==> ALPHA-INDEPENDENCE NOT CLEARLY SUPPORTED")
            print(f"    Significant variation in alpha across decoders.")

        print()

        # Pairwise alpha ratios
        print(f"  Pairwise alpha comparisons:")
        for i in range(len(decoder_names_with_alpha)):
            for j in range(i+1, len(decoder_names_with_alpha)):
                name_i = decoder_names_with_alpha[i]
                name_j = decoder_names_with_alpha[j]
                a_i = alpha_values[i]
                a_j = alpha_values[j]
                ratio = a_i / a_j if abs(a_j) > 1e-10 else float('inf')
                diff = abs(a_i - a_j)
                combined_err = np.sqrt(alpha_errors[i]**2 + alpha_errors[j]**2)
                sig = diff / combined_err if combined_err > 0 else float('inf')
                print(f"    {name_i:14s} vs {name_j:14s}: "
                      f"ratio={ratio:.3f}, diff={diff:.6f}, {sig:.1f}sigma")
        print()

    # Print beta values
    print(f"  Beta (intercept) values -- expected to differ between decoders:")
    for name in decoder_names_with_alpha:
        a = alphas[name]
        print(f"    {name:14s}: beta = {a['beta']:.6f} +/- {a['beta_err']:.6f}")
    print()

    if len(alpha_values) >= 2:
        beta_values = [alphas[n]['beta'] for n in decoder_names_with_alpha]
        beta_arr = np.array(beta_values)
        beta_range = beta_arr.max() - beta_arr.min()
        print(f"    Beta range: {beta_range:.6f}")
        print(f"    (Beta SHOULD vary across decoders -- it captures decoder quality)")
        print()


# ============================================================================
# Visualization
# ============================================================================

def generate_figures(results, decoder_names, alphas, alphas_by_p,
                     p_fixed=0.005, eps_fixed=0.20):
    """Generate 4-panel figure for alpha-independence test."""
    print("=" * 75)
    print("  GENERATING FIGURES")
    print("=" * 75)
    print()

    # Color scheme for decoders
    decoder_colors = {
        'MWPM':       '#1565C0',
        'Unwt-MWPM':  '#E65100',
        'Union-Find': '#2E7D32',
        'BP+OSD':     '#7B1FA2',
        'Greedy-NN':  '#C62828',
        'Lookup':     '#00838F',
    }
    decoder_markers = {
        'MWPM':       'o',
        'Unwt-MWPM':  's',
        'Union-Find': 'D',
        'BP+OSD':     '^',
        'Greedy-NN':  'v',
        'Lookup':     'P',
    }

    # Pre-build the epsilon label to avoid f-string backslash issues
    eps_pct = int(eps_fixed * 100)
    eps_label = r'$\varepsilon=' + str(eps_pct) + r'\%$'

    fig = plt.figure(figsize=(16, 10))
    gs = GridSpec(2, 2, figure=fig, wspace=0.30, hspace=0.38,
                  left=0.07, right=0.97, bottom=0.08, top=0.90)

    # ================================================================
    # Panel (a): ln R vs d for all decoders at fixed p=0.005, eps=20%
    # ================================================================
    ax_a = fig.add_subplot(gs[0, 0])

    for name in decoder_names:
        if name not in alphas:
            continue
        a = alphas[name]
        color = decoder_colors.get(name, '#333333')
        marker = decoder_markers.get(name, 'x')

        ax_a.errorbar(
            a['d_values'], a['lnR_values'],
            yerr=a['lnR_errors'],
            color=color, marker=marker, markersize=7,
            linestyle='none', linewidth=1.5, capsize=3,
            label=name, zorder=5,
        )

        # Plot fit line
        d_fit = np.linspace(min(DISTANCES) - 0.5, max(DISTANCES) + 0.5, 50)
        lnR_fit = a['alpha'] * d_fit + a['beta']
        ax_a.plot(d_fit, lnR_fit, color=color, linestyle='--', linewidth=1.2,
                  alpha=0.6, zorder=3)

    ax_a.set_xlabel('Code distance $d$', fontsize=12)
    ax_a.set_ylabel('$\\ln R$', fontsize=12)
    ax_a.set_title(f'(a) $\\ln R$ vs $d$  [$p={p_fixed}$, {eps_label}]',
                   fontsize=13, fontweight='bold')
    ax_a.legend(fontsize=8, loc='best', framealpha=0.9)
    ax_a.grid(True, alpha=0.3)
    ax_a.set_xticks(DISTANCES)

    # Annotate with slopes
    alpha_vals = [alphas[n]['alpha'] for n in decoder_names if n in alphas]
    if alpha_vals:
        alpha_mean = np.mean(alpha_vals)
        alpha_std = np.std(alpha_vals)
        ax_a.annotate(
            f'Slopes (alpha):\n'
            f'mean = {alpha_mean:.4f}\n'
            f'std = {alpha_std:.4f}\n'
            f'CV = {alpha_std/abs(alpha_mean):.3f}' if abs(alpha_mean) > 1e-10 else 'N/A',
            xy=(0.03, 0.97), xycoords='axes fraction', fontsize=8,
            ha='left', va='top',
            bbox=dict(boxstyle='round,pad=0.3', fc='#E8F5E9', ec='green', alpha=0.9))

    # ================================================================
    # Panel (b): alpha vs decoder type (bar plot with error bars)
    # ================================================================
    ax_b = fig.add_subplot(gs[0, 1])

    names_with_alpha = [n for n in decoder_names if n in alphas]
    x_pos = np.arange(len(names_with_alpha))
    alpha_vals_plot = [alphas[n]['alpha'] for n in names_with_alpha]
    alpha_errs_plot = [alphas[n]['alpha_err'] for n in names_with_alpha]
    bar_colors = [decoder_colors.get(n, '#333333') for n in names_with_alpha]

    bars = ax_b.bar(x_pos, alpha_vals_plot, yerr=alpha_errs_plot,
                    color=bar_colors, edgecolor='white', linewidth=0.5,
                    capsize=5, alpha=0.85, zorder=4)

    # Draw mean line
    if len(alpha_vals_plot) >= 2:
        alpha_mean = np.mean(alpha_vals_plot)
        ax_b.axhline(y=alpha_mean, color='red', linestyle='--', linewidth=1.5,
                     label=f'Mean = {alpha_mean:.4f}', zorder=3)
        # Shade +/- 1 std
        alpha_std = np.std(alpha_vals_plot)
        ax_b.axhspan(alpha_mean - alpha_std, alpha_mean + alpha_std,
                     color='red', alpha=0.1, zorder=2)
        ax_b.legend(fontsize=9, loc='best')

    ax_b.set_xticks(x_pos)
    ax_b.set_xticklabels(names_with_alpha, rotation=30, ha='right', fontsize=9)
    ax_b.set_ylabel('$\\alpha$ (slope of $\\ln R$ vs $d$)', fontsize=12)
    ax_b.set_title(f'(b) $\\alpha$ vs Decoder Type  [$p={p_fixed}$, {eps_label}]',
                   fontsize=13, fontweight='bold')
    ax_b.grid(True, alpha=0.3, axis='y')

    if len(alpha_vals_plot) >= 2:
        CV = np.std(alpha_vals_plot) / abs(np.mean(alpha_vals_plot)) if abs(np.mean(alpha_vals_plot)) > 1e-10 else float('inf')
        verdict = "SUPPORTED" if CV < 0.10 else ("WEAK" if CV < 0.20 else "NOT SUPPORTED")
        ax_b.annotate(
            f'CV = {CV:.3f}\nalpha-indep: {verdict}',
            xy=(0.97, 0.97), xycoords='axes fraction', fontsize=9,
            ha='right', va='top',
            bbox=dict(boxstyle='round,pad=0.3',
                      fc='#E8F5E9' if CV < 0.10 else ('#FFF3E0' if CV < 0.20 else '#FFEBEE'),
                      ec='green' if CV < 0.10 else ('orange' if CV < 0.20 else 'red'),
                      alpha=0.9))

    # ================================================================
    # Panel (c): beta vs decoder type
    # ================================================================
    ax_c = fig.add_subplot(gs[1, 0])

    beta_vals_plot = [alphas[n]['beta'] for n in names_with_alpha]
    beta_errs_plot = [alphas[n]['beta_err'] for n in names_with_alpha]

    bars_c = ax_c.bar(x_pos, beta_vals_plot, yerr=beta_errs_plot,
                      color=bar_colors, edgecolor='white', linewidth=0.5,
                      capsize=5, alpha=0.85, zorder=4)

    ax_c.set_xticks(x_pos)
    ax_c.set_xticklabels(names_with_alpha, rotation=30, ha='right', fontsize=9)
    ax_c.set_ylabel('$\\beta$ (intercept)', fontsize=12)
    ax_c.set_title(f'(c) $\\beta$ vs Decoder Type  [$p={p_fixed}$, {eps_label}]',
                   fontsize=13, fontweight='bold')
    ax_c.grid(True, alpha=0.3, axis='y')

    if len(beta_vals_plot) >= 2:
        beta_range = max(beta_vals_plot) - min(beta_vals_plot)
        ax_c.annotate(
            f'$\\beta$ range = {beta_range:.4f}\n(SHOULD vary across decoders)',
            xy=(0.97, 0.97), xycoords='axes fraction', fontsize=9,
            ha='right', va='top',
            bbox=dict(boxstyle='round,pad=0.3', fc='#FFF3E0', ec='orange', alpha=0.9))

    # ================================================================
    # Panel (d): R_decoder / R_MWPM vs d -- should be constant if alpha same
    # ================================================================
    ax_d = fig.add_subplot(gs[1, 1])

    ref_decoder = 'MWPM'
    if ref_decoder in results:
        for name in decoder_names:
            if name == ref_decoder or name not in results:
                continue
            color = decoder_colors.get(name, '#333333')
            marker = decoder_markers.get(name, 'x')

            d_vals = []
            ratio_vals = []
            ratio_errs = []

            for d in DISTANCES:
                if (d not in results[name] or p_fixed not in results[name][d] or
                    d not in results[ref_decoder] or p_fixed not in results[ref_decoder][d]):
                    continue

                R_name = results[name][d][p_fixed]['ps'].get(eps_fixed, {}).get('R', np.nan)
                R_ref = results[ref_decoder][d][p_fixed]['ps'].get(eps_fixed, {}).get('R', np.nan)
                R_name_ci = results[name][d][p_fixed]['ps'].get(eps_fixed, {}).get('R_ci', np.nan)
                R_ref_ci = results[ref_decoder][d][p_fixed]['ps'].get(eps_fixed, {}).get('R_ci', np.nan)

                if np.isnan(R_name) or np.isnan(R_ref) or R_ref <= 0:
                    continue

                ratio = R_name / R_ref
                # Propagate errors for ratio
                if R_name > 0 and R_ref > 0 and not np.isnan(R_name_ci) and not np.isnan(R_ref_ci):
                    rel_var = (R_name_ci / R_name)**2 + (R_ref_ci / R_ref)**2
                    ratio_err = ratio * np.sqrt(rel_var)
                else:
                    ratio_err = 0.0

                d_vals.append(d)
                ratio_vals.append(ratio)
                ratio_errs.append(ratio_err)

            if d_vals:
                ax_d.errorbar(
                    d_vals, ratio_vals,
                    yerr=ratio_errs,
                    color=color, marker=marker, markersize=7,
                    linestyle='-', linewidth=1.5, capsize=3,
                    label=f'{name}/MWPM', zorder=5,
                )

    ax_d.axhline(y=1.0, color='gray', linestyle='--', linewidth=1, alpha=0.5, zorder=2)
    ax_d.set_xlabel('Code distance $d$', fontsize=12)
    ax_d.set_ylabel(f'$R_{{\\mathrm{{decoder}}}} / R_{{\\mathrm{{{ref_decoder}}}}}$', fontsize=12)
    ax_d.set_title(f'(d) $R$ Ratio vs $d$  [$p={p_fixed}$, {eps_label}]',
                   fontsize=13, fontweight='bold')
    ax_d.legend(fontsize=8, loc='best', framealpha=0.9)
    ax_d.grid(True, alpha=0.3)
    ax_d.set_xticks(DISTANCES)

    ax_d.annotate(
        'If $\\alpha$ is decoder-independent,\nratios should be constant in $d$',
        xy=(0.03, 0.03), xycoords='axes fraction', fontsize=8,
        ha='left', va='bottom',
        bbox=dict(boxstyle='round,pad=0.3', fc='lightyellow', ec='orange', alpha=0.9))

    # Supertitle
    fig.suptitle(
        'Decoder $\\alpha$-Independence Test for Post-Selection Improvement\n'
        '$\\ln R(\\varepsilon, d) = \\alpha(p) \\cdot d + \\beta(\\varepsilon)$: '
        '$\\alpha$ should be decoder-independent',
        fontsize=14, fontweight='bold', y=0.98)

    fig.savefig(FIG_PATH_PNG, dpi=200, bbox_inches='tight')
    fig.savefig(FIG_PATH_PDF, bbox_inches='tight')
    plt.close(fig)

    print(f"  Saved: {FIG_PATH_PNG}")
    print(f"  Saved: {FIG_PATH_PDF}")
    print()


# ============================================================================
# Summary
# ============================================================================

def print_summary(results, decoder_names, alphas, alphas_by_p):
    """Print comprehensive summary."""
    print()
    print("=" * 75)
    print("  COMPREHENSIVE SUMMARY")
    print("=" * 75)
    print()

    # 1. Decoder hierarchy check
    print("  1. DECODER HIERARCHY (sanity check)")
    print("  " + "-" * 50)
    for d in DISTANCES:
        print(f"    d={d}:")
        for p in [0.005]:
            rates = {}
            for name in decoder_names:
                if d in results.get(name, {}) and p in results[name][d]:
                    rates[name] = results[name][d][p]['base_rate']
            if rates:
                sorted_dec = sorted(rates.items(), key=lambda x: x[1])
                ranking = ' < '.join([f'{n}({r:.5f})' for n, r in sorted_dec])
                print(f"      p={p}: {ranking}")
    print()

    # 2. Alpha independence across all p values
    print("  2. ALPHA VALUES ACROSS ALL P VALUES")
    print("  " + "-" * 50)
    for p in P_VALUES:
        a_p = alphas_by_p.get(p, {})
        if not a_p:
            continue
        vals = [a_p[n]['alpha'] for n in decoder_names if n in a_p]
        if len(vals) >= 2:
            mean_a = np.mean(vals)
            std_a = np.std(vals)
            cv = std_a / abs(mean_a) if abs(mean_a) > 1e-10 else float('inf')
            status = "OK" if cv < 0.10 else ("WEAK" if cv < 0.20 else "FAIL")
            print(f"    p={p:.4f}: alpha_mean={mean_a:.5f}, "
                  f"alpha_std={std_a:.5f}, CV={cv:.3f} [{status}]")
    print()

    # 3. Overall verdict
    all_cvs = []
    for p in P_VALUES:
        a_p = alphas_by_p.get(p, {})
        vals = [a_p[n]['alpha'] for n in decoder_names if n in a_p]
        if len(vals) >= 2:
            mean_a = np.mean(vals)
            std_a = np.std(vals)
            if abs(mean_a) > 1e-10:
                all_cvs.append(std_a / abs(mean_a))

    if all_cvs:
        mean_cv = np.mean(all_cvs)
        print(f"  3. OVERALL VERDICT")
        print(f"  " + "-" * 50)
        print(f"    Mean CV across all p values: {mean_cv:.4f}")
        print(f"    Number of decoders tested: {len(decoder_names)}")
        print(f"    Number of decoder pairs: {len(decoder_names) * (len(decoder_names)-1) // 2}")
        print()
        if mean_cv < 0.10:
            print(f"    CONCLUSION: ALPHA-INDEPENDENCE IS SUPPORTED (CV < 0.10)")
            print(f"    The post-selection improvement slope alpha is consistent")
            print(f"    across {len(decoder_names)} different decoders, confirming that")
            print(f"    it is determined by the physical noise channel (Sigma),")
            print(f"    not by the decoder choice.")
        elif mean_cv < 0.20:
            print(f"    CONCLUSION: WEAK SUPPORT FOR ALPHA-INDEPENDENCE (0.10 < CV < 0.20)")
            print(f"    Alpha shows some variation but remains relatively consistent.")
        else:
            print(f"    CONCLUSION: ALPHA-INDEPENDENCE NOT CLEARLY SUPPORTED (CV > 0.20)")
            print(f"    Significant variation in alpha observed across decoders.")
            print(f"    This may indicate decoder-dependent corrections to the")
            print(f"    leading-order prediction, or insufficient statistics.")
    print()

    # 4. Physical interpretation
    print("  4. PHYSICAL INTERPRETATION")
    print("  " + "-" * 50)
    print("""
    Prediction 1 states: ln R(eps, d) = alpha(p) * d + beta(eps)
    where alpha depends only on the physical error rate p (through
    the entropy production Sigma of the noise channel), while beta
    captures the decoder-dependent baseline improvement.

    If alpha is decoder-independent:
      - The d-scaling of post-selection improvement is UNIVERSAL
      - It is determined by the physical noise channel alone
      - This supports the retrodiction-theoretic interpretation:
        syndrome weight ~ Sigma (entropy production), and
        filtering on low Sigma improves fidelity by exp(-Sigma/2)
        regardless of how the recovery is performed.

    The beta term captures decoder quality:
      - Better decoders (MWPM > Unwt > Greedy) have different baselines
      - But the RATE at which post-selection helps (per unit d) is the same
""")


# ============================================================================
# Main
# ============================================================================

def main():
    # Force line-buffered output so progress is visible in real time
    sys.stdout.reconfigure(line_buffering=True)

    t_start = time.time()

    print()
    print("*" * 75)
    print("*  DECODER ALPHA-INDEPENDENCE TEST                                    *")
    print("*  Testing: ln R(eps,d) = alpha(p)*d + beta(eps)                      *")
    print("*  Prediction: alpha is INDEPENDENT of decoder choice                 *")
    print("*" * 75)
    print()
    print(f"  Distances: {DISTANCES}")
    print(f"  P values: {P_VALUES}")
    print(f"  Post-selection fractions: {PS_FRACTIONS}")
    print(f"  Shots per point: {SHOTS_PER_POINT:,}")
    print(f"  Decoders: MWPM, Unwt-MWPM, Union-Find, BP+OSD, Greedy-NN, Lookup(d=3 only)")
    print()

    # Run simulation
    results, decoder_names = run_simulation()

    # Extract alphas
    p_fixed = 0.005
    eps_fixed = 0.20
    alphas = extract_alphas(results, decoder_names, p_fixed=p_fixed, eps_fixed=eps_fixed)
    alphas_by_p = extract_alphas_all_p(results, decoder_names, eps_fixed=eps_fixed)

    # Analysis
    print_analysis(results, decoder_names, alphas)

    # Generate figures
    generate_figures(results, decoder_names, alphas, alphas_by_p)

    # Summary
    print_summary(results, decoder_names, alphas, alphas_by_p)

    t_total = time.time() - t_start
    print(f"  Total execution time: {t_total:.1f} seconds")
    print()


if __name__ == '__main__':
    main()
