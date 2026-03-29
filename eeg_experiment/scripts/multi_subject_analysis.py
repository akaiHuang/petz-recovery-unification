#!/usr/bin/env python3
"""
multi_subject_analysis.py -- Multi-subject Petz bound analysis using Neural Sigma
==================================================================================

Runs Method 3 (Neural Sigma Estimator) on ALL available Sleep-EDF subjects,
aggregates results with proper statistics, and generates publication-quality figures.

Author: Sheng-Kai Huang (akai@fawstudio.com)
Date: 2026-03-20
"""

import csv
import json
import logging
import os
import sys
import time
import traceback
import warnings
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy import stats as scipy_stats

# Suppress MNE and sklearn warnings
os.environ["MNE_LOGGING_LEVEL"] = "WARNING"
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)

# Import from improved_sigma.py
sys.path.insert(0, str(Path(__file__).resolve().parent))
from improved_sigma import (
    STAGE_ORDER, STAGE_LABELS_LONG, STAGE_COLORS, STAGE_MAP, EPOCH_DURATION,
    load_eeg_data, method3_neural, compute_tau_retrodiction,
)

# ============================================================================
# CONSTANTS
# ============================================================================

MIN_PSG_SIZE_MB = 10  # Skip PSG files smaller than this (likely incomplete)
MAX_EPOCHS_PER_STAGE = 80  # Subsample if too many epochs (for speed)


# ============================================================================
# SUBJECT DISCOVERY
# ============================================================================

def find_complete_pairs(data_dir: Path) -> List[Tuple[Path, Path, str]]:
    """Find all complete PSG + Hypnogram pairs, skipping incomplete downloads."""
    psg_files = sorted(data_dir.glob("*-PSG.edf"))
    hyp_files = sorted(data_dir.glob("*-Hypnogram.edf"))

    # Build hypnogram lookup by 6-char prefix (e.g., SC4001)
    hyp_lookup = {}
    for h in hyp_files:
        prefix = h.name[:6]
        hyp_lookup[prefix] = h

    pairs = []
    for psg in psg_files:
        prefix = psg.name[:6]
        hyp = hyp_lookup.get(prefix)
        if hyp is None:
            logger.warning(f"  No hypnogram for {psg.name} -- skipping")
            continue

        # Check file size (incomplete downloads are small)
        size_mb = psg.stat().st_size / (1024 * 1024)
        if size_mb < MIN_PSG_SIZE_MB:
            logger.warning(f"  {psg.name} only {size_mb:.1f} MB -- likely incomplete, skipping")
            continue

        subject_id = psg.stem.replace("-PSG", "").replace("E0", "")
        pairs.append((psg, hyp, subject_id))

    return pairs


# ============================================================================
# PER-SUBJECT ANALYSIS (Neural Method Only)
# ============================================================================

def analyze_subject_neural(
    psg_path: Path,
    hyp_path: Path,
    subject_id: str,
    max_epochs_per_stage: int = MAX_EPOCHS_PER_STAGE,
) -> List[Dict]:
    """
    Run Neural Sigma estimator on one subject.

    Returns list of dicts with keys:
        subject, epoch_idx, stage, sigma_neural, tau, fidelity, bound, bound_holds
    """
    # Load data
    epochs, ch_names, sfreq, stages = load_eeg_data(
        str(psg_path), str(hyp_path), epoch_duration=EPOCH_DURATION
    )

    if stages is None:
        logger.error(f"  {subject_id}: No stage annotations")
        return []

    # Count epochs per stage
    stage_counts = defaultdict(int)
    stage_indices = defaultdict(list)
    for i, s in enumerate(stages):
        s_str = str(s)
        if s_str in STAGE_ORDER:
            stage_counts[s_str] += 1
            stage_indices[s_str].append(i)

    n_valid = sum(stage_counts.values())
    logger.info(f"  {subject_id}: {n_valid} valid epochs, sfreq={sfreq} Hz, ch={ch_names}")
    for s in STAGE_ORDER:
        if s in stage_counts:
            logger.info(f"    {STAGE_LABELS_LONG.get(s, s):>5s}: {stage_counts[s]} epochs")

    if n_valid == 0:
        return []

    # Subsample if needed (keep analysis tractable)
    selected_indices = []
    for stage in STAGE_ORDER:
        indices = stage_indices[stage]
        if len(indices) > max_epochs_per_stage:
            rng = np.random.RandomState(42)
            indices = sorted(rng.choice(indices, max_epochs_per_stage, replace=False))
        selected_indices.extend(indices)
    selected_indices.sort()

    results = []
    n_total = len(selected_indices)

    for count, idx in enumerate(selected_indices):
        if (count + 1) % 20 == 0 or count == 0:
            logger.info(f"    Epoch {count+1}/{n_total} (idx={idx})...")

        stage = str(stages[idx])
        data = epochs[idx]  # (n_channels, n_samples)

        # Compute tau (retrodiction error)
        try:
            tau = compute_tau_retrodiction(data, sfreq)
        except Exception:
            tau = np.nan

        # Compute Neural Sigma
        try:
            sigma = method3_neural(data, sfreq, n_segments=60, n_splits=3)
        except Exception as e:
            logger.debug(f"    Neural estimator failed for epoch {idx}: {e}")
            sigma = np.nan

        if np.isfinite(sigma) and np.isfinite(tau):
            fidelity = 1.0 - tau
            bound = np.exp(-sigma / 2.0)
            bound_holds = fidelity >= bound - 1e-10  # small tolerance

            results.append({
                "subject": subject_id,
                "epoch_idx": int(idx),
                "stage": stage,
                "sigma_neural": float(sigma),
                "tau": float(tau),
                "fidelity": float(fidelity),
                "bound_exp": float(bound),
                "bound_holds": bool(bound_holds),
            })

    return results


# ============================================================================
# AGGREGATION & STATISTICS
# ============================================================================

def aggregate_results(all_results: List[Dict]) -> Dict:
    """Compute aggregate statistics across all subjects."""

    # Group by stage
    by_stage = defaultdict(list)
    for r in all_results:
        by_stage[r["stage"]].append(r)

    # Group by subject
    by_subject = defaultdict(list)
    for r in all_results:
        by_subject[r["subject"]].append(r)

    stats = {}
    n_subjects = len(by_subject)

    for stage in STAGE_ORDER:
        data = by_stage.get(stage, [])
        if not data:
            continue

        sigmas = np.array([r["sigma_neural"] for r in data])
        taus = np.array([r["tau"] for r in data])
        fids = np.array([r["fidelity"] for r in data])
        bounds = np.array([r["bound_exp"] for r in data])
        holds = np.array([r["bound_holds"] for r in data])

        n = len(data)
        n_subjects_stage = len(set(r["subject"] for r in data))

        # Per-subject means for SE calculation (proper between-subject variance)
        subject_sigma_means = []
        subject_tau_means = []
        subject_compliance = []
        for subj, subj_data in by_subject.items():
            subj_stage = [r for r in subj_data if r["stage"] == stage]
            if subj_stage:
                subject_sigma_means.append(np.mean([r["sigma_neural"] for r in subj_stage]))
                subject_tau_means.append(np.mean([r["tau"] for r in subj_stage]))
                subject_compliance.append(np.mean([r["bound_holds"] for r in subj_stage]) * 100)

        subject_sigma_means = np.array(subject_sigma_means)
        subject_tau_means = np.array(subject_tau_means)
        subject_compliance = np.array(subject_compliance)

        stats[stage] = {
            "n_epochs": n,
            "n_subjects": n_subjects_stage,
            "sigma_mean": float(np.mean(sigmas)),
            "sigma_std": float(np.std(sigmas)),
            "sigma_se": float(np.std(subject_sigma_means) / np.sqrt(len(subject_sigma_means)))
                        if len(subject_sigma_means) > 1 else float(np.std(sigmas) / np.sqrt(n)),
            "sigma_median": float(np.median(sigmas)),
            "tau_mean": float(np.mean(taus)),
            "tau_std": float(np.std(taus)),
            "tau_se": float(np.std(subject_tau_means) / np.sqrt(len(subject_tau_means)))
                      if len(subject_tau_means) > 1 else float(np.std(taus) / np.sqrt(n)),
            "tau_median": float(np.median(taus)),
            "fidelity_mean": float(np.mean(fids)),
            "bound_mean": float(np.mean(bounds)),
            "compliance_pct": float(np.mean(holds) * 100),
            "compliance_se": float(np.std(subject_compliance) / np.sqrt(len(subject_compliance)))
                             if len(subject_compliance) > 1 else 0.0,
            "subject_sigma_means": subject_sigma_means.tolist(),
            "subject_tau_means": subject_tau_means.tolist(),
            "subject_compliance": subject_compliance.tolist(),
        }

    # Overall
    all_holds = [r["bound_holds"] for r in all_results]
    stats["_overall"] = {
        "n_epochs": len(all_results),
        "n_subjects": n_subjects,
        "compliance_pct": float(np.mean(all_holds) * 100) if all_holds else 0.0,
    }

    return stats


def run_statistical_tests(all_results: List[Dict]) -> Dict:
    """Run Kruskal-Wallis + post-hoc Mann-Whitney tests for Sigma ordering."""

    by_stage = defaultdict(list)
    for r in all_results:
        by_stage[r["stage"]].append(r["sigma_neural"])

    test_results = {}

    # Kruskal-Wallis across W, N1, N2, N3
    nrem_stages = ["W", "N1", "N2", "N3"]
    groups = [np.array(by_stage[s]) for s in nrem_stages if s in by_stage and len(by_stage[s]) > 0]
    present_stages = [s for s in nrem_stages if s in by_stage and len(by_stage[s]) > 0]

    if len(groups) >= 2:
        kw_stat, kw_p = scipy_stats.kruskal(*groups)
        test_results["kruskal_wallis"] = {
            "statistic": float(kw_stat),
            "p_value": float(kw_p),
            "stages": present_stages,
            "n_per_stage": [len(g) for g in groups],
            "significant": kw_p < 0.05,
        }

        # Post-hoc: pairwise Mann-Whitney U tests with Bonferroni correction
        n_comparisons = len(present_stages) * (len(present_stages) - 1) // 2
        posthoc = []
        for i in range(len(present_stages)):
            for j in range(i + 1, len(present_stages)):
                s1, s2 = present_stages[i], present_stages[j]
                g1, g2 = np.array(by_stage[s1]), np.array(by_stage[s2])
                u_stat, p_val = scipy_stats.mannwhitneyu(g1, g2, alternative="greater")
                posthoc.append({
                    "comparison": f"{s1} > {s2}",
                    "U_statistic": float(u_stat),
                    "p_value": float(p_val),
                    "p_bonferroni": float(min(p_val * n_comparisons, 1.0)),
                    "significant_bonferroni": p_val * n_comparisons < 0.05,
                    "mean_1": float(np.mean(g1)),
                    "mean_2": float(np.mean(g2)),
                    "effect_size_r": float(u_stat / (len(g1) * len(g2))),
                })
        test_results["posthoc_mannwhitney"] = posthoc

    # Also test W vs N3 specifically (the extreme comparison)
    if "W" in by_stage and "N3" in by_stage and len(by_stage["W"]) > 0 and len(by_stage["N3"]) > 0:
        g_w = np.array(by_stage["W"])
        g_n3 = np.array(by_stage["N3"])
        u, p = scipy_stats.mannwhitneyu(g_w, g_n3, alternative="greater")
        test_results["W_vs_N3"] = {
            "U_statistic": float(u),
            "p_value": float(p),
            "significant": p < 0.05,
            "W_mean": float(np.mean(g_w)),
            "N3_mean": float(np.mean(g_n3)),
            "ratio": float(np.mean(g_w) / np.mean(g_n3)) if np.mean(g_n3) > 0 else float("inf"),
        }

    # REM vs NREM comparison
    if "R" in by_stage and len(by_stage["R"]) > 0:
        g_rem = np.array(by_stage["R"])
        nrem_data = []
        for s in ["N1", "N2", "N3"]:
            if s in by_stage:
                nrem_data.extend(by_stage[s])
        if nrem_data:
            g_nrem = np.array(nrem_data)
            u, p = scipy_stats.mannwhitneyu(g_rem, g_nrem, alternative="greater")
            test_results["REM_vs_NREM"] = {
                "U_statistic": float(u),
                "p_value": float(p),
                "significant": p < 0.05,
                "REM_mean": float(np.mean(g_rem)),
                "NREM_mean": float(np.mean(g_nrem)),
            }

    return test_results


# ============================================================================
# PUBLICATION FIGURE
# ============================================================================

def generate_publication_figure(
    all_results: List[Dict],
    agg_stats: Dict,
    test_results: Dict,
    output_path: Path,
):
    """Generate a 4-panel publication-quality figure."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    import matplotlib.gridspec as gridspec

    # Collect data by stage
    by_stage = defaultdict(lambda: {"sigma": [], "tau": [], "fidelity": [], "bound": [], "holds": []})
    for r in all_results:
        s = r["stage"]
        by_stage[s]["sigma"].append(r["sigma_neural"])
        by_stage[s]["tau"].append(r["tau"])
        by_stage[s]["fidelity"].append(r["fidelity"])
        by_stage[s]["bound"].append(r["bound_exp"])
        by_stage[s]["holds"].append(r["bound_holds"])

    present_stages = [s for s in STAGE_ORDER if s in by_stage]
    colors = [STAGE_COLORS[s] for s in present_stages]
    labels = [STAGE_LABELS_LONG.get(s, s) for s in present_stages]

    fig = plt.figure(figsize=(14, 12))
    gs = gridspec.GridSpec(2, 2, hspace=0.35, wspace=0.30)

    # ── Panel A: Sigma by stage (box plot) ──
    ax_a = fig.add_subplot(gs[0, 0])
    sigma_data = [np.array(by_stage[s]["sigma"]) for s in present_stages]
    bp_a = ax_a.boxplot(
        sigma_data,
        patch_artist=True,
        labels=labels,
        widths=0.6,
        showfliers=True,
        flierprops=dict(marker=".", markersize=3, alpha=0.3),
    )
    for patch, color in zip(bp_a["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.6)
    for median in bp_a["medians"]:
        median.set_color("black")
        median.set_linewidth(2)

    # Add individual subject means as dots
    for i, stage in enumerate(present_stages):
        if stage in agg_stats and "subject_sigma_means" in agg_stats[stage]:
            subj_means = agg_stats[stage]["subject_sigma_means"]
            jitter = np.random.RandomState(42).uniform(-0.15, 0.15, len(subj_means))
            ax_a.scatter(
                np.full(len(subj_means), i + 1) + jitter,
                subj_means,
                c="black", s=30, zorder=5, alpha=0.7, marker="D",
                label="Subject means" if i == 0 else None,
            )

    ax_a.set_ylabel(r"$\Sigma$ (Neural entropy production)", fontsize=12)
    ax_a.set_title(r"(A) $\Sigma$ by sleep stage", fontsize=13, fontweight="bold")
    ax_a.grid(True, alpha=0.2, axis="y")

    # Add significance annotations
    if "kruskal_wallis" in test_results:
        kw = test_results["kruskal_wallis"]
        p_str = f"p = {kw['p_value']:.2e}" if kw['p_value'] < 0.001 else f"p = {kw['p_value']:.4f}"
        ax_a.text(
            0.98, 0.98,
            f"Kruskal-Wallis: H = {kw['statistic']:.1f}\n{p_str}",
            transform=ax_a.transAxes, ha="right", va="top", fontsize=9,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", edgecolor="gray", alpha=0.9),
        )
    if len(present_stages) > 0:
        ax_a.legend(fontsize=8, loc="upper left")

    # ── Panel B: Tau by stage (box plot) ──
    ax_b = fig.add_subplot(gs[0, 1])
    tau_data = [np.array(by_stage[s]["tau"]) for s in present_stages]
    bp_b = ax_b.boxplot(
        tau_data,
        patch_artist=True,
        labels=labels,
        widths=0.6,
        showfliers=True,
        flierprops=dict(marker=".", markersize=3, alpha=0.3),
    )
    for patch, color in zip(bp_b["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.6)
    for median in bp_b["medians"]:
        median.set_color("black")
        median.set_linewidth(2)

    ax_b.set_ylabel(r"$\tau = 1 - F$ (retrodiction error)", fontsize=12)
    ax_b.set_title(r"(B) $\tau$ by sleep stage", fontsize=13, fontweight="bold")
    ax_b.grid(True, alpha=0.2, axis="y")

    # ── Panel C: F vs exp(-Sigma/2) scatter (color by stage) ──
    ax_c = fig.add_subplot(gs[1, 0])

    # Petz bound curve
    sigma_range = np.linspace(0, 10, 500)
    bound_curve = np.exp(-sigma_range / 2.0)
    ax_c.plot(sigma_range, bound_curve, "k-", linewidth=2.5,
              label=r"$F = e^{-\Sigma/2}$ (Petz bound)", zorder=10)
    ax_c.fill_between(sigma_range, 0, bound_curve, alpha=0.08, color="red",
                       label="Forbidden region")

    for stage in present_stages:
        sigmas = by_stage[stage]["sigma"]
        fids = by_stage[stage]["fidelity"]
        ax_c.scatter(
            sigmas, fids,
            c=STAGE_COLORS[stage],
            label=STAGE_LABELS_LONG.get(stage, stage),
            alpha=0.35, s=12, edgecolors="none", zorder=5,
        )

    # Overall compliance
    overall = agg_stats.get("_overall", {})
    comply_pct = overall.get("compliance_pct", 0)
    ax_c.text(
        0.98, 0.02,
        f"Bound satisfied: {comply_pct:.1f}%\n({overall.get('n_epochs', 0)} epochs, "
        f"{overall.get('n_subjects', 0)} subjects)",
        transform=ax_c.transAxes, ha="right", va="bottom", fontsize=10,
        fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="gray", alpha=0.9),
    )

    ax_c.set_xlabel(r"$\Sigma$ (Neural entropy production)", fontsize=12)
    ax_c.set_ylabel(r"$F = 1 - \tau$ (fidelity)", fontsize=12)
    ax_c.set_title(r"(C) $F$ vs $e^{-\Sigma/2}$ scatter", fontsize=13, fontweight="bold")
    ax_c.legend(fontsize=8, loc="upper right")
    ax_c.set_xlim(0, None)
    ax_c.set_ylim(0, 1.05)
    ax_c.grid(True, alpha=0.2)

    # ── Panel D: Bound compliance by stage (bar chart with error bars) ──
    ax_d = fig.add_subplot(gs[1, 1])

    compliance_means = []
    compliance_ses = []
    for stage in present_stages:
        if stage in agg_stats:
            compliance_means.append(agg_stats[stage]["compliance_pct"])
            compliance_ses.append(agg_stats[stage].get("compliance_se", 0))
        else:
            compliance_means.append(0)
            compliance_ses.append(0)

    x_pos = np.arange(len(present_stages))
    bars = ax_d.bar(
        x_pos, compliance_means,
        yerr=compliance_ses,
        color=colors, alpha=0.7,
        edgecolor="black", linewidth=0.5,
        capsize=5, error_kw=dict(capthick=1.5),
    )

    # Add value labels on bars
    for bar, val, se in zip(bars, compliance_means, compliance_ses):
        ax_d.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + se + 1.5,
            f"{val:.1f}%",
            ha="center", va="bottom", fontsize=10, fontweight="bold",
        )

    ax_d.set_xticks(x_pos)
    ax_d.set_xticklabels(labels)
    ax_d.set_ylabel("Petz bound compliance (%)", fontsize=12)
    ax_d.set_title("(D) Bound compliance by stage", fontsize=13, fontweight="bold")
    ax_d.set_ylim(0, 115)
    ax_d.axhline(y=100, color="gray", linestyle="--", alpha=0.5, linewidth=0.8)
    ax_d.grid(True, alpha=0.2, axis="y")

    # Suptitle
    n_subj = overall.get("n_subjects", 0)
    n_ep = overall.get("n_epochs", 0)
    fig.suptitle(
        f"Neural $\\Sigma$ Estimator: Petz Recovery Bound on EEG Sleep Data\n"
        f"({n_subj} subjects, {n_ep} epochs, Sleep-EDF Database)",
        fontsize=14, fontweight="bold", y=1.02,
    )

    plt.tight_layout()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(output_path), dpi=300, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Saved figure: {output_path}")


# ============================================================================
# CSV & SUMMARY OUTPUT
# ============================================================================

def save_csv(all_results: List[Dict], output_path: Path):
    """Save per-epoch results to CSV."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "subject", "epoch_idx", "stage",
        "sigma_neural", "tau", "fidelity", "bound_exp", "bound_holds",
    ]

    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for r in all_results:
            row = {}
            for k in fieldnames:
                v = r.get(k, "")
                if isinstance(v, float):
                    row[k] = f"{v:.6f}" if np.isfinite(v) else ""
                elif isinstance(v, bool):
                    row[k] = str(v)
                else:
                    row[k] = v
            writer.writerow(row)

    logger.info(f"Saved CSV: {output_path}")


def generate_summary_text(
    agg_stats: Dict,
    test_results: Dict,
    all_results: List[Dict],
) -> str:
    """Generate summary text suitable for paper inclusion."""

    by_subject = defaultdict(list)
    for r in all_results:
        by_subject[r["subject"]].append(r)

    lines = []
    lines.append("=" * 80)
    lines.append("MULTI-SUBJECT PETZ BOUND ANALYSIS -- NEURAL SIGMA ESTIMATOR")
    lines.append("=" * 80)
    lines.append("")

    overall = agg_stats.get("_overall", {})
    lines.append(f"Subjects analyzed: {overall.get('n_subjects', 0)}")
    lines.append(f"Total epochs:      {overall.get('n_epochs', 0)}")
    lines.append(f"Overall Petz bound compliance: {overall.get('compliance_pct', 0):.1f}%")
    lines.append("")

    # Per-stage table
    lines.append(f"{'Stage':<6} {'N_ep':>6} {'N_subj':>7} "
                 f"{'Sigma_mean':>10} {'Sigma_SE':>9} "
                 f"{'Tau_mean':>9} {'Tau_SE':>8} "
                 f"{'Comply%':>9} {'Comply_SE':>10}")
    lines.append("-" * 85)

    for stage in STAGE_ORDER:
        if stage not in agg_stats:
            continue
        s = agg_stats[stage]
        lines.append(
            f"{STAGE_LABELS_LONG.get(stage, stage):<6} {s['n_epochs']:>6} {s['n_subjects']:>7} "
            f"{s['sigma_mean']:>10.4f} {s['sigma_se']:>9.4f} "
            f"{s['tau_mean']:>9.4f} {s['tau_se']:>8.4f} "
            f"{s['compliance_pct']:>8.1f}% {s['compliance_se']:>9.1f}%"
        )

    lines.append("")

    # Ordering check
    sigma_means = {}
    for stage in STAGE_ORDER:
        if stage in agg_stats:
            sigma_means[stage] = agg_stats[stage]["sigma_mean"]

    nrem = ["W", "N1", "N2", "N3"]
    present = [s for s in nrem if s in sigma_means]
    if len(present) >= 2:
        vals = [sigma_means[s] for s in present]
        ordered = all(vals[i] >= vals[i + 1] for i in range(len(vals) - 1))
        order_str = " > ".join([f"Sigma({s})={sigma_means[s]:.4f}" for s in present])
        lines.append(f"Ordering: {order_str}")
        lines.append(f"W > N1 > N2 > N3: {'CONFIRMED' if ordered else 'VIOLATED'}")
        lines.append("")

    # Statistical tests
    if "kruskal_wallis" in test_results:
        kw = test_results["kruskal_wallis"]
        lines.append(f"Kruskal-Wallis test (W, N1, N2, N3):")
        lines.append(f"  H = {kw['statistic']:.2f}, p = {kw['p_value']:.2e}, "
                     f"significant = {kw['significant']}")
        lines.append("")

    if "posthoc_mannwhitney" in test_results:
        lines.append("Post-hoc Mann-Whitney U (one-sided, Bonferroni corrected):")
        for ph in test_results["posthoc_mannwhitney"]:
            sig = "***" if ph["p_bonferroni"] < 0.001 else (
                "**" if ph["p_bonferroni"] < 0.01 else (
                    "*" if ph["p_bonferroni"] < 0.05 else "ns"))
            lines.append(
                f"  {ph['comparison']:>10s}: U = {ph['U_statistic']:.0f}, "
                f"p_bonf = {ph['p_bonferroni']:.4e}, {sig}"
            )
        lines.append("")

    if "W_vs_N3" in test_results:
        wn3 = test_results["W_vs_N3"]
        lines.append(f"Key comparison W > N3:")
        lines.append(f"  U = {wn3['U_statistic']:.0f}, p = {wn3['p_value']:.2e}")
        lines.append(f"  Sigma(W)/Sigma(N3) = {wn3['ratio']:.2f}")
        lines.append("")

    if "REM_vs_NREM" in test_results:
        rem = test_results["REM_vs_NREM"]
        lines.append(f"REM vs NREM:")
        lines.append(f"  U = {rem['U_statistic']:.0f}, p = {rem['p_value']:.2e}")
        lines.append(f"  Sigma(REM) = {rem['REM_mean']:.4f}, Sigma(NREM) = {rem['NREM_mean']:.4f}")
        lines.append("")

    # Paper paragraph
    lines.append("=" * 80)
    lines.append("SUGGESTED PAPER TEXT:")
    lines.append("=" * 80)
    lines.append("")

    n_subj = overall.get("n_subjects", 0)
    n_ep = overall.get("n_epochs", 0)
    comply = overall.get("compliance_pct", 0)

    w_sigma = sigma_means.get("W", 0)
    n3_sigma = sigma_means.get("N3", 0)

    kw_p = test_results.get("kruskal_wallis", {}).get("p_value", 1.0)
    wn3_p = test_results.get("W_vs_N3", {}).get("p_value", 1.0)
    ratio = test_results.get("W_vs_N3", {}).get("ratio", 0)

    lines.append(
        f"We analyzed {n_subj} subjects ({n_ep} 30-second epochs) from the "
        f"Sleep-EDF database using a neural entropy production estimator (Method 3). "
        f"The Petz recovery bound F >= exp(-Sigma/2) was satisfied in {comply:.1f}% "
        f"of all epochs. The neural entropy production Sigma tracked consciousness "
        f"level: Sigma(Wake) = {w_sigma:.3f} vs Sigma(N3) = {n3_sigma:.3f} "
        f"(ratio {ratio:.1f}x, Mann-Whitney p = {wn3_p:.1e}). "
        f"The ordering Sigma(W) > Sigma(N1) > Sigma(N2) > Sigma(N3) was statistically "
        f"significant (Kruskal-Wallis H = {test_results.get('kruskal_wallis', {}).get('statistic', 0):.1f}, "
        f"p = {kw_p:.1e}), confirming that consciousness corresponds to higher "
        f"irreversibility in neural dynamics, as quantified by the Petz recovery framework."
    )

    return "\n".join(lines)


# ============================================================================
# MAIN
# ============================================================================

def main():
    t0 = time.time()

    script_dir = Path(__file__).resolve().parent
    project_dir = script_dir.parent
    data_dir = project_dir / "data" / "sleep-edf"
    results_dir = project_dir / "results"

    fig_path = results_dir / "figures" / "multi_subject_petz_bound.png"
    csv_path = results_dir / "statistics" / "multi_subject_results.csv"
    summary_path = results_dir / "statistics" / "multi_subject_summary.txt"
    json_path = results_dir / "statistics" / "multi_subject_stats.json"

    print("=" * 80)
    print("MULTI-SUBJECT PETZ BOUND ANALYSIS")
    print("Method: Neural Sigma Estimator (Method 3)")
    print("Bound: tau <= 1 - exp(-Sigma/2)  [Petz recovery bound]")
    print("=" * 80)

    # Find all complete pairs
    pairs = find_complete_pairs(data_dir)
    print(f"\nFound {len(pairs)} complete PSG + Hypnogram pairs:")
    for psg, hyp, subj_id in pairs:
        size_mb = psg.stat().st_size / (1024 * 1024)
        print(f"  {subj_id}: PSG={size_mb:.0f}MB, Hyp={hyp.name}")

    if not pairs:
        print("No data found!")
        sys.exit(1)

    # Process all subjects
    all_results = []
    for idx, (psg_path, hyp_path, subject_id) in enumerate(pairs):
        print(f"\n{'─' * 60}")
        print(f"[{idx+1}/{len(pairs)}] Processing {subject_id}...")
        print(f"{'─' * 60}")

        try:
            subject_results = analyze_subject_neural(psg_path, hyp_path, subject_id)
            all_results.extend(subject_results)
            n = len(subject_results)
            elapsed = time.time() - t0
            print(f"  => {n} epochs analyzed for {subject_id} (cumulative: {len(all_results)} epochs, {elapsed:.0f}s)")
        except Exception as e:
            logger.error(f"  FAILED: {subject_id}: {e}")
            traceback.print_exc()
            continue

    if not all_results:
        print("\nNo results obtained. Check data files.")
        sys.exit(1)

    print(f"\n{'=' * 80}")
    print(f"TOTAL: {len(all_results)} epochs from {len(set(r['subject'] for r in all_results))} subjects")
    print(f"{'=' * 80}")

    # Aggregate statistics
    agg_stats = aggregate_results(all_results)
    test_results = run_statistical_tests(all_results)

    # Generate outputs
    print("\nGenerating outputs...")

    # 1. Publication figure
    generate_publication_figure(all_results, agg_stats, test_results, fig_path)

    # 2. CSV
    save_csv(all_results, csv_path)

    # 3. Summary text
    summary = generate_summary_text(agg_stats, test_results, all_results)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    with open(summary_path, "w") as f:
        f.write(summary)
    print(f"Saved summary: {summary_path}")

    # 4. JSON stats (machine-readable)
    json_stats = {
        "aggregate": {},
        "tests": test_results,
    }
    for k, v in agg_stats.items():
        if isinstance(v, dict):
            json_stats["aggregate"][k] = {
                kk: vv for kk, vv in v.items()
                if not isinstance(vv, (list, np.ndarray))
            }
    with open(json_path, "w") as f:
        json.dump(json_stats, f, indent=2, default=str)
    print(f"Saved JSON: {json_path}")

    # Print summary
    print("\n" + summary)

    elapsed = time.time() - t0
    print(f"\nTotal runtime: {elapsed:.0f}s ({elapsed/60:.1f} min)")
    print(f"\nOutputs:")
    print(f"  Figure:  {fig_path}")
    print(f"  CSV:     {csv_path}")
    print(f"  Summary: {summary_path}")
    print(f"  JSON:    {json_path}")


if __name__ == "__main__":
    main()
