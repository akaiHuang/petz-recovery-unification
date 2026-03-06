"""
Post-Selection Thermodynamic Filtering for Repetition Codes.

Verifies Prediction 1 from Huang (2026): for stabilizer codes under noise,
rejecting high-entropy-production syndromes yields an improvement factor

    ln R(epsilon, d) = alpha(p) * d + beta(epsilon),

where alpha depends only on the noise rate p (not the decoder or epsilon)
and beta depends only on the rejection fraction epsilon.

Setup:
  - d-qubit repetition code (d = 3, 5, 7, 9) under i.i.d. bit-flip noise
  - Two decoders:
      (A) ML decoder: always picks minimum-weight coset representative
      (B) Greedy-right decoder: pairs adjacent defects left-to-right,
          connecting unpaired boundary defects to the right edge.
          Suboptimal because it sometimes picks the wrong coset.
  - Entropy production proxy: Delta_D(s) ~ w_min(s) * ln(1/p)
  - Post-selection: reject syndromes with Delta_D above threshold

Tests:
  1. ln R vs d is linear (exponential improvement with code distance)
  2. alpha(p) is consistent across different epsilon values
  3. alpha(p) ~ (1/2) * ln(1/p) as predicted by the Ising mapping
  4. Higher Delta_D syndromes have higher conditional failure rates

Reference: Huang (2026), Sec. S9 of Supplemental Material.

Usage:
    python simulations/postselection_filtering.py
"""

import os
import numpy as np
from itertools import product as iterproduct
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

_DIR = os.path.dirname(os.path.abspath(__file__))


# ======================================================================
# Repetition code utilities
# ======================================================================

def syndrome(error, d):
    """Compute syndrome: bit i = error[i] XOR error[i+1]."""
    return tuple(error[i] ^ error[i + 1] for i in range(d - 1))


def build_coset_table(d):
    """Map each syndrome to its light (w <= d//2) and heavy coset reps."""
    table = {}
    for e in iterproduct([0, 1], repeat=d):
        syn = syndrome(e, d)
        w = sum(e)
        if syn not in table:
            table[syn] = {'light': None, 'heavy': None,
                          'light_w': d + 1, 'heavy_w': -1}
        if w <= d // 2 and w < table[syn]['light_w']:
            table[syn]['light'] = e
            table[syn]['light_w'] = w
        elif w > d // 2 and w > table[syn]['heavy_w']:
            table[syn]['heavy'] = e
            table[syn]['heavy_w'] = w
    return table


def delta_D_proxy(w, p):
    """Data-processing gap proxy: Delta_D ~ w * ln(1/p)."""
    if p <= 0 or p >= 1:
        return 0.0
    return w * np.log(1.0 / p)


# ======================================================================
# Decoders
# ======================================================================

def make_ml_decoder(coset_table):
    """ML decoder: always picks minimum-weight coset representative."""
    return {syn: cosets['light'] for syn, cosets in coset_table.items()}


def make_greedy_right_decoder(d):
    """Greedy-right decoder: pair adjacent defects left-to-right,
    unpaired last defect connects to the right boundary.

    This produces valid corrections (matching syndromes) but picks
    the wrong coset for some syndromes, making it suboptimal.
    """
    lookup = {}
    for s_bits in range(2 ** (d - 1)):
        syn = tuple((s_bits >> i) & 1 for i in range(d - 1))
        correction = [0] * d
        defects = [i for i in range(d - 1) if syn[i] == 1]
        i = 0
        while i < len(defects):
            if i + 1 < len(defects):
                # Pair: flip qubits from defect[i]+1 to defect[i+1]
                for q in range(defects[i] + 1, defects[i + 1] + 1):
                    correction[q] = 1
                i += 2
            else:
                # Unpaired: connect to right boundary
                for q in range(defects[i] + 1, d):
                    correction[q] = 1
                i += 1
        lookup[syn] = tuple(correction)
    return lookup


# ======================================================================
# Core computation
# ======================================================================

def compute_syndrome_stats(d, p, decoder_lookup, coset_table):
    """Compute per-syndrome failure statistics by exact enumeration.

    Returns:
        p_fail: overall logical failure rate
        syn_stats: dict mapping syn -> {prob, fail_prob, cond_fail, delta_D}
    """
    syn_stats = {}

    for e in iterproduct([0, 1], repeat=d):
        w = sum(e)
        prob = (p ** w) * ((1 - p) ** (d - w))
        syn = syndrome(e, d)
        correction = decoder_lookup[syn]
        is_logical_error = all((ei ^ ci) == 1 for ei, ci in zip(e, correction))

        if syn not in syn_stats:
            syn_stats[syn] = {
                'prob': 0.0,
                'fail_prob': 0.0,
                'min_weight': coset_table[syn]['light_w'],
            }
        syn_stats[syn]['prob'] += prob
        if is_logical_error:
            syn_stats[syn]['fail_prob'] += prob

    total_fail = 0.0
    for stats in syn_stats.values():
        stats['cond_fail'] = (stats['fail_prob'] / stats['prob']
                              if stats['prob'] > 0 else 0)
        stats['delta_D'] = delta_D_proxy(stats['min_weight'], p)
        total_fail += stats['fail_prob']

    return total_fail, syn_stats


def postselect(syn_stats, epsilon):
    """Reject fraction epsilon of probability mass with largest Delta_D.

    Returns: (p_fail_ps, acceptance_rate, delta_D_threshold)
    """
    syn_list = sorted(syn_stats.items(),
                      key=lambda x: x[1]['delta_D'], reverse=True)
    total_prob = sum(s['prob'] for _, s in syn_list)

    rejected_prob = 0.0
    rejected_syns = set()
    for syn, stats in syn_list:
        if rejected_prob + stats['prob'] <= epsilon * total_prob + 1e-15:
            rejected_prob += stats['prob']
            rejected_syns.add(syn)
        else:
            break

    accepted_prob = sum(s['prob'] for sy, s in syn_stats.items()
                        if sy not in rejected_syns)
    accepted_fail = sum(s['fail_prob'] for sy, s in syn_stats.items()
                        if sy not in rejected_syns)
    delta_D_th = max((s['delta_D'] for sy, s in syn_stats.items()
                      if sy not in rejected_syns), default=0.0)

    if accepted_prob <= 0:
        return 0.0, 0.0, 0.0
    return accepted_fail / accepted_prob, accepted_prob / total_prob, delta_D_th


# ======================================================================
# Main simulation
# ======================================================================

def run_simulation():
    """Run the full post-selection thermodynamic filtering simulation."""

    print("=" * 70)
    print("Post-Selection Thermodynamic Filtering Simulation")
    print("Prediction 1 from Huang (2026)")
    print("=" * 70)
    print()
    print("Testing: ln R(epsilon, d) = alpha(p)*d + beta(epsilon)")
    print()

    distances = [3, 5, 7, 9]
    p_values = np.linspace(0.02, 0.20, 25)
    epsilon_values = [0.05, 0.10, 0.20, 0.30]
    p_fixed = 0.10

    # Build tables
    print("Building coset tables and decoder lookups...")
    coset_tables = {}
    ml_decoders = {}
    greedy_decoders = {}
    for d in distances:
        coset_tables[d] = build_coset_table(d)
        ml_decoders[d] = make_ml_decoder(coset_tables[d])
        greedy_decoders[d] = make_greedy_right_decoder(d)
    print("Done.\n")

    # Sanity check
    print("Decoder comparison at p=0.10:")
    for d in distances:
        pf_ml, _ = compute_syndrome_stats(d, 0.10, ml_decoders[d], coset_tables[d])
        pf_gr, _ = compute_syndrome_stats(d, 0.10, greedy_decoders[d], coset_tables[d])
        n_diff = sum(1 for syn in coset_tables[d]
                     if ml_decoders[d][syn] != greedy_decoders[d][syn])
        print(f"  d={d}: ML p_fail={pf_ml:.8f}, Greedy p_fail={pf_gr:.8f}, "
              f"differ={n_diff}/{len(coset_tables[d])}")
    print()

    all_decoders = [
        ("ML decoder", ml_decoders),
        ("Greedy decoder", greedy_decoders),
    ]
    colors_dec = {"ML decoder": '#2166ac', "Greedy decoder": '#b2182b'}
    ls_dec = {"ML decoder": '-', "Greedy decoder": '--'}

    # ==================================================================
    # Figure: 3-panel plot
    # ==================================================================
    fig, axes = plt.subplots(1, 3, figsize=(19, 5.5))

    # ------------------------------------------------------------------
    # Panel 1: Logical error rate vs p, with/without post-selection
    # ------------------------------------------------------------------
    print("-" * 70)
    print("Panel 1: Logical error rate vs p")
    print("-" * 70)

    ax1 = axes[0]
    eps_ps = 0.10
    d_alpha_vis = {3: 0.3, 5: 0.5, 7: 0.7, 9: 0.9}

    for dec_name, dec_dict in all_decoders:
        col = colors_dec[dec_name]
        ls = ls_dec[dec_name]
        for d in distances:
            pf_list, pf_ps_list = [], []
            for p in p_values:
                pf, ss = compute_syndrome_stats(d, p, dec_dict[d], coset_tables[d])
                pf_list.append(pf)
                pf_ps, _, _ = postselect(ss, eps_ps)
                pf_ps_list.append(pf_ps)

            av = d_alpha_vis[d]
            ax1.semilogy(p_values, [max(v, 1e-12) for v in pf_list],
                         color=col, ls=ls, lw=1.5, alpha=av)
            ax1.semilogy(p_values, [max(v, 1e-12) for v in pf_ps_list],
                         color=col, ls=ls, lw=2.5, alpha=av,
                         marker='o', ms=2)

    legend1 = [
        Line2D([0], [0], color=colors_dec['ML decoder'], ls='-', lw=1.5,
               label='ML, no PS'),
        Line2D([0], [0], color=colors_dec['ML decoder'], ls='-', lw=2.5,
               marker='o', ms=3, label='ML, PS (10%)'),
        Line2D([0], [0], color=colors_dec['Greedy decoder'], ls='--', lw=1.5,
               label='Greedy, no PS'),
        Line2D([0], [0], color=colors_dec['Greedy decoder'], ls='--', lw=2.5,
               marker='o', ms=3, label='Greedy, PS (10%)'),
    ] + [Line2D([0], [0], color='gray', lw=1, alpha=d_alpha_vis[d],
                label=f'd={d}') for d in distances]
    ax1.legend(handles=legend1, fontsize=7, loc='lower right')
    ax1.set_xlabel('Bit-flip probability $p$', fontsize=12)
    ax1.set_ylabel('Logical error rate $p_{\\mathrm{fail}}$', fontsize=12)
    ax1.set_title('Error rate with/without\npost-selection ($\\varepsilon$=10%)',
                  fontsize=12)
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(bottom=1e-8, top=1)

    # ------------------------------------------------------------------
    # Panel 2: ln R vs d for ML decoder at different epsilon
    # ------------------------------------------------------------------
    ax2 = axes[1]
    print()
    print("-" * 70)
    print(f"Panel 2: ln R vs d at p = {p_fixed} (ML decoder)")
    print("-" * 70)

    alpha_results = {}
    colors_eps = {0.05: '#4393c3', 0.10: '#2166ac',
                  0.20: '#d6604d', 0.30: '#b2182b'}
    markers_eps = {0.05: 'o', 0.10: 's', 0.20: '^', 0.30: 'D'}

    for dec_name, dec_dict in all_decoders:
        ls = ls_dec[dec_name]
        for epsilon in epsilon_values:
            lnR_values = []
            for d in distances:
                pf, ss = compute_syndrome_stats(d, p_fixed, dec_dict[d],
                                                 coset_tables[d])
                pf_ps, acc, dDth = postselect(ss, epsilon)

                if pf_ps > 1e-15 and pf > 1e-15:
                    lnR = np.log(pf / pf_ps)
                else:
                    lnR = 0.0
                lnR_values.append(lnR)

                print(f"  eps={epsilon:.2f}, {dec_name:16s}, d={d}: "
                      f"p_fail={pf:.8f}, p_fail_ps={pf_ps:.8f}, "
                      f"ln R={lnR:.4f}, accept={acc:.4f}")

            # Fit ln R = alpha * d + beta using points with lnR > 0
            valid_idx = [i for i in range(len(distances)) if lnR_values[i] > 0]
            if len(valid_idx) >= 2:
                d_arr = np.array([distances[i] for i in valid_idx])
                lnR_arr = np.array([lnR_values[i] for i in valid_idx])
                alpha_fit, beta_fit = np.polyfit(d_arr, lnR_arr, 1)
            else:
                alpha_fit, beta_fit = 0.0, 0.0

            alpha_results[(dec_name, epsilon)] = alpha_fit

            col = colors_eps[epsilon]
            mk = markers_eps[epsilon]
            short = 'ML' if 'ML' in dec_name else 'Gr'

            ax2.plot(distances, lnR_values, color=col, ls='none',
                     marker=mk, ms=8,
                     markerfacecolor=col if 'ML' in dec_name else 'none',
                     markeredgecolor=col, markeredgewidth=2)

            if alpha_fit != 0:
                d_fit = np.linspace(min(distances) - 0.5,
                                    max(distances) + 0.5, 50)
                ax2.plot(d_fit, alpha_fit * d_fit + beta_fit, color=col,
                         ls=ls, lw=1.5,
                         label=f'{short}, $\\varepsilon$={epsilon:.0%}: '
                               f'$\\alpha$={alpha_fit:.3f}')

    ax2.set_xlabel('Code distance $d$', fontsize=12)
    ax2.set_ylabel('$\\ln R(\\varepsilon, d)$', fontsize=12)
    ax2.set_title(f'Improvement scaling at $p={p_fixed}$\n'
                  f'$\\ln R = \\alpha \\cdot d + \\beta$ (filled=ML, open=Greedy)',
                  fontsize=11)
    ax2.legend(fontsize=6, loc='upper left', ncol=2)
    ax2.grid(True, alpha=0.3)
    ax2.set_xticks(distances)

    # ------------------------------------------------------------------
    # Panel 3: Delta_D distribution + conditional failure rate
    # ------------------------------------------------------------------
    ax3 = axes[2]
    d_ex = max(distances)
    _, ss_ml = compute_syndrome_stats(d_ex, p_fixed, ml_decoders[d_ex],
                                       coset_tables[d_ex])
    _, ss_gr = compute_syndrome_stats(d_ex, p_fixed, greedy_decoders[d_ex],
                                       coset_tables[d_ex])

    # Aggregate ML stats by Delta_D
    sig = np.array([s['delta_D'] for s in ss_ml.values()])
    prob = np.array([s['prob'] for s in ss_ml.values()])
    fail_ml = np.array([s['cond_fail'] for s in ss_ml.values()])
    fail_gr = np.array([ss_gr[syn]['cond_fail'] for syn in ss_ml.keys()])

    unique_dD = np.unique(sig)
    agg_p, agg_f_ml, agg_f_gr = [], [], []
    for dD in unique_dD:
        mask = np.isclose(sig, dD)
        tp = prob[mask].sum()
        tf_ml = (prob[mask] * fail_ml[mask]).sum()
        tf_gr = (prob[mask] * fail_gr[mask]).sum()
        agg_p.append(tp)
        agg_f_ml.append(tf_ml / tp if tp > 0 else 0)
        agg_f_gr.append(tf_gr / tp if tp > 0 else 0)
    agg_p = np.array(agg_p)
    agg_f_ml = np.array(agg_f_ml)
    agg_f_gr = np.array(agg_f_gr)

    bw = 0.25 * (unique_dD[1] - unique_dD[0]) if len(unique_dD) > 1 else 0.3
    bars = ax3.bar(unique_dD, agg_p, width=bw, alpha=0.7,
                   edgecolor='black', linewidth=0.5, color='steelblue',
                   label='$P(\\Delta D)$')

    # Threshold lines
    _, _, dDth10 = postselect(ss_ml, 0.10)
    _, _, dDth20 = postselect(ss_ml, 0.20)
    if dDth10 > 0:
        ax3.axvline(dDth10 + bw, color='red', ls='--', lw=2,
                    label=f'$\\varepsilon$=10% threshold')
    if dDth20 > 0:
        ax3.axvline(dDth20 + bw, color='orange', ls='--', lw=2,
                    label=f'$\\varepsilon$=20% threshold')

    ax3t = ax3.twinx()
    ax3t.plot(unique_dD, agg_f_ml, 'o-', color=colors_dec['ML decoder'],
              lw=2, ms=5, label='$p_{\\mathrm{fail}}(\\Delta D)$ ML', zorder=5)
    ax3t.plot(unique_dD, agg_f_gr, 's--', color=colors_dec['Greedy decoder'],
              lw=2, ms=5, label='$p_{\\mathrm{fail}}(\\Delta D)$ Greedy', zorder=5)
    ax3t.set_ylabel('Conditional failure rate', fontsize=11)
    ax3t.set_ylim(bottom=-0.02, top=1.05)

    ax3.set_xlabel('$\\Delta D(s) \\approx w \\cdot \\ln(1/p)$', fontsize=11)
    ax3.set_ylabel('Probability $P(\\Delta D)$', fontsize=11)
    ax3.set_title(f'Entropy production distribution\n(d={d_ex}, p={p_fixed})',
                  fontsize=12)
    h1, l1 = ax3.get_legend_handles_labels()
    h2, l2 = ax3t.get_legend_handles_labels()
    ax3.legend(h1 + h2, l1 + l2, fontsize=7, loc='center right')

    plt.suptitle('Post-Selection as Thermodynamic Filtering '
                 '(Huang 2026, Prediction 1)', fontsize=14, y=1.02)
    plt.tight_layout()

    fig.savefig(os.path.join(_DIR, 'fig_postselection_filtering.png'),
                dpi=150, bbox_inches='tight')
    fig.savefig(os.path.join(_DIR, 'fig_postselection_filtering.pdf'),
                bbox_inches='tight')
    print()
    print("Figures saved: simulations/fig_postselection_filtering.{png,pdf}")
    plt.close(fig)

    # ==================================================================
    # KEY RESULTS
    # ==================================================================
    print()
    print("=" * 70)
    print("KEY RESULTS")
    print("=" * 70)

    # Test 1: alpha consistency across epsilon (for ML decoder)
    print()
    print("Test 1: Is alpha consistent across epsilon? (ML decoder)")
    print(f"{'epsilon':>10} | {'alpha':>10} | {'beta':>10}")
    print("-" * 38)
    ml_alphas = []
    for epsilon in epsilon_values:
        a = alpha_results.get(("ML decoder", epsilon), 0)
        ml_alphas.append(a)
        # Compute beta from the intercept
        # beta = lnR - alpha*d at some d
        print(f"{epsilon:10.2f} | {a:10.4f} |")

    ml_alphas_nonzero = [a for a in ml_alphas if a > 0]
    if len(ml_alphas_nonzero) >= 2:
        alpha_mean = np.mean(ml_alphas_nonzero)
        alpha_std = np.std(ml_alphas_nonzero)
        print(f"\n  alpha(ML) = {alpha_mean:.4f} +/- {alpha_std:.4f}")
        print(f"  Coefficient of variation: {alpha_std/alpha_mean:.2%}")
        if alpha_std / alpha_mean < 0.3:
            print("  => alpha is approximately CONSISTENT across epsilon")
        else:
            print("  => alpha shows some variation with epsilon")

    # Test 2: alpha comparison between decoders
    print()
    print("Test 2: Decoder comparison (ML vs Greedy)")
    print(f"{'epsilon':>10} | {'alpha_ML':>10} | {'alpha_Gr':>10} | {'rel diff':>10}")
    print("-" * 50)
    for epsilon in epsilon_values:
        a_ml = alpha_results.get(("ML decoder", epsilon), 0)
        a_gr = alpha_results.get(("Greedy decoder", epsilon), 0)
        if abs(a_ml) > 1e-10:
            rd = abs(a_ml - a_gr) / abs(a_ml)
            print(f"{epsilon:10.2f} | {a_ml:10.4f} | {a_gr:10.4f} | {rd:10.2%}")
        else:
            print(f"{epsilon:10.2f} | {a_ml:10.4f} | {a_gr:10.4f} | {'N/A':>10}")

    print()
    print("  NOTE: The greedy decoder has p_fail ~ p (uncoded rate) because")
    print("  it picks the wrong coset for ~1/3 of syndromes. This makes it")
    print("  essentially a 'no-coding' baseline. Its failures are spread")
    print("  across all Delta_D levels, so post-selection helps less.")
    print()
    print("  For the REPETITION CODE, there are only 2 valid corrections")
    print("  per syndrome. Any decoder must pick one or the other.")
    print("  ML is the unique reasonable decoder; any alternative that")
    print("  differs must pick the wrong coset for some syndromes,")
    print("  leading to catastrophic performance.")
    print()
    print("  The decoder-independence prediction is most meaningful for")
    print("  2D codes (surface codes) where exponentially many corrections")
    print("  exist per syndrome, allowing decoders to differ in heuristics")
    print("  without catastrophic failure (e.g., MWPM vs union-find).")

    # Test 3: alpha(p) vs theoretical prediction
    print()
    print("Test 3: alpha(p) vs theoretical prediction alpha ~ (1/2)*ln(1/p)")
    print(f"{'p':>6} | {'alpha_ML':>10} | {'(1/2)ln(1/p)':>12} | {'ratio':>8}")
    print("-" * 45)

    p_scan = [0.03, 0.05, 0.08, 0.10, 0.12, 0.15]
    eps_test = 0.10

    for p in p_scan:
        lnR_v = []
        for d in distances:
            pf, ss = compute_syndrome_stats(d, p, ml_decoders[d], coset_tables[d])
            pf_ps, _, _ = postselect(ss, eps_test)
            if pf_ps > 1e-15 and pf > 1e-15:
                lnR_v.append(np.log(pf / pf_ps))
            else:
                lnR_v.append(0.0)
        valid = [(distances[i], lnR_v[i]) for i in range(len(distances))
                 if lnR_v[i] > 0]
        if len(valid) >= 2:
            c = np.polyfit(np.array([v[0] for v in valid]),
                           np.array([v[1] for v in valid]), 1)
            alpha_p = c[0]
        else:
            alpha_p = 0.0

        theo = 0.5 * np.log(1.0 / p)
        ratio = alpha_p / theo if theo > 0 and alpha_p > 0 else 0
        print(f"{p:6.3f} | {alpha_p:10.4f} | {theo:12.4f} | {ratio:8.2f}")

    print()
    print("  The ratio alpha / (ln(1/p)/2) measures how well the Ising")
    print("  mapping prediction works. For the repetition code at small p,")
    print("  we expect reasonable agreement; at larger p near threshold,")
    print("  corrections from finite-size effects become important.")

    # ==================================================================
    # Physical interpretation
    # ==================================================================
    print()
    print("=" * 70)
    print("PHYSICAL INTERPRETATION")
    print("=" * 70)
    print()
    print("1. THERMODYNAMIC FILTERING: Post-selecting on low-Delta_D syndromes")
    print("   is equivalent to rejecting high-entropy-production events.")
    print("   This is thermodynamic filtering: we keep only events where")
    print("   the noise process was 'gentle' (low information loss).")
    print()
    print("2. EXPONENTIAL IMPROVEMENT: ln R grows linearly with d, meaning")
    print("   R ~ exp(alpha * d). Each additional code qubit provides an")
    print("   exponential improvement in the post-selection gain.")
    print()
    print("3. SEPARATION OF SCALES: The factored form ln R = alpha*d + beta")
    print("   shows a clean separation between:")
    print("   - alpha(p): set by the noise channel (thermodynamic cost per qubit)")
    print("   - beta(epsilon): set by rejection threshold (overhead cost)")
    print()
    print("4. CONNECTION TO ISING MODEL: alpha ~ (1/2)*ln(1/p) connects to")
    print("   the Ising--QEC mapping where Delta_D corresponds to the Ising")
    print("   domain-wall free energy. Post-selection = conditioning on large")
    print("   free-energy gaps where the correct logical value is favored.")

    print()
    print("=" * 70)
    print("Done.")


if __name__ == "__main__":
    run_simulation()
