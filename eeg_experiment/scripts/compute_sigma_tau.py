#!/usr/bin/env python3
"""
Compute retrodiction metrics Σ (relative entropy production) and τ (retrodiction error)
from EEG data, grouped by consciousness state.

Theory (Huang 2025, Paper 1):
    τ = 1 - F          (retrodiction error = 1 - fidelity)
    F ≥ exp(-Σ/2)      (Petz recovery bound)
    τ ≤ 1 - exp(-Σ/2)  (equivalent form)

Operational definitions for EEG:
    Σ: KL divergence between forward and backward transition matrices
       of the discretized EEG state (binned power spectrum)
    τ: Retrodiction error — how well can we reconstruct epoch t-1 from epoch t
       using the Petz recovery map (optimal Bayesian retrodiction)

Usage:
    python scripts/compute_sigma_tau.py
"""

import os
import sys
import warnings
from pathlib import Path

import numpy as np

os.environ["MNE_LOGGING_LEVEL"] = "WARNING"
import mne

# Reuse constants from load_and_preview
STAGE_MAP = {
    "Sleep stage W": "Wake",
    "Sleep stage 1": "N1",
    "Sleep stage 2": "N2",
    "Sleep stage 3": "N3",
    "Sleep stage 4": "N3",
    "Sleep stage R": "REM",
    "Sleep stage ?": "Unknown",
    "Movement time": "Movement",
}
STAGE_ORDER = ["Wake", "N1", "N2", "N3", "REM"]
STAGE_COLORS = {
    "Wake": "#e74c3c",
    "N1": "#f39c12",
    "N2": "#2ecc71",
    "N3": "#3498db",
    "REM": "#9b59b6",
}
EPOCH_DURATION = 30.0


def load_eeg_and_stages(psg_path: Path, hyp_path: Path):
    """Load EEG and parse into list of (epoch_data, stage_label) tuples."""
    raw = mne.io.read_raw_edf(str(psg_path), preload=True, verbose=False)

    # Pick EEG channels
    eeg_picks = [ch for ch in raw.ch_names if "EEG" in ch]
    if eeg_picks:
        raw.pick(eeg_picks)
    raw.filter(0.5, 45.0, verbose=False)

    # Parse annotations
    annot = mne.read_annotations(str(hyp_path))
    raw.set_annotations(annot)

    sfreq = raw.info["sfreq"]
    epochs = []

    for ann in raw.annotations:
        stage = STAGE_MAP.get(ann["description"], None)
        if stage not in STAGE_ORDER:
            continue
        onset = ann["onset"]
        duration = ann["duration"]
        n_epochs = int(duration // EPOCH_DURATION)
        for i in range(n_epochs):
            t_start = onset + i * EPOCH_DURATION
            start_sample = int(t_start * sfreq)
            n_samples = int(EPOCH_DURATION * sfreq)
            if start_sample < 0 or start_sample + n_samples > raw.n_times:
                continue
            data = raw.get_data(start=start_sample, stop=start_sample + n_samples)
            signal = data.mean(axis=0)  # Average across channels
            epochs.append((signal, stage, t_start))

    return epochs, sfreq


def signal_to_state_vector(signal: np.ndarray, sfreq: float,
                            n_bins: int = 20) -> np.ndarray:
    """
    Convert a 30-second EEG signal to a discrete probability distribution
    over frequency bins (normalized power spectrum).

    This is the "state" for retrodiction purposes.
    """
    from scipy.signal import welch

    freqs, psd = welch(signal, fs=sfreq, nperseg=min(256, len(signal)))

    # Keep frequencies 0.5-45 Hz
    mask = (freqs >= 0.5) & (freqs <= 45.0)
    psd = psd[mask]
    freqs = freqs[mask]

    # Bin into n_bins frequency bins
    bin_edges = np.linspace(0.5, 45.0, n_bins + 1)
    binned = np.zeros(n_bins)
    for i in range(n_bins):
        idx = (freqs >= bin_edges[i]) & (freqs < bin_edges[i + 1])
        if idx.any():
            binned[i] = np.trapezoid(psd[idx], freqs[idx])

    # Normalize to probability distribution
    total = binned.sum()
    if total > 0:
        binned /= total
    else:
        binned = np.ones(n_bins) / n_bins

    # Add small epsilon to avoid log(0)
    eps = 1e-12
    binned = binned + eps
    binned /= binned.sum()

    return binned


def compute_kl_divergence(p: np.ndarray, q: np.ndarray) -> float:
    """KL divergence D(p || q)."""
    mask = (p > 0) & (q > 0)
    return np.sum(p[mask] * np.log(p[mask] / q[mask]))


def compute_sigma_forward_backward(states: list[np.ndarray]) -> list[float]:
    """
    Compute Σ for each consecutive pair of epochs.

    Σ_t = D(p_t || p_{t-1}) + D(p_{t-1} || p_t)
        = symmetrized KL divergence (Jeffreys divergence)

    This measures the irreversibility of the t-1 -> t transition:
    how distinguishable is the forward process from the backward process.
    """
    sigmas = []
    for i in range(1, len(states)):
        p_prev = states[i - 1]
        p_curr = states[i]
        # Symmetrized KL = forward + backward
        sigma = compute_kl_divergence(p_curr, p_prev) + \
                compute_kl_divergence(p_prev, p_curr)
        sigmas.append(sigma)
    return sigmas


def compute_tau_retrodiction(states: list[np.ndarray]) -> list[float]:
    """
    Compute τ (retrodiction error) for each consecutive pair.

    The Petz recovery map for classical distributions is the Bayesian retrodiction:
        R(p_curr) = p_prior * (T^dagger p_curr / p_prior)
    For our discrete PSD states, we implement:
        - Forward: observe p_t (current state)
        - Retrodiction: estimate p_{t-1} from p_t using Bayesian update
        - τ = 1 - F(p_{t-1}, R(p_t))

    Fidelity for classical distributions (Bhattacharyya coefficient):
        F(p, q) = (Σ sqrt(p_i * q_i))^2
    """
    taus = []
    for i in range(1, len(states)):
        p_prev = states[i - 1]  # True previous state
        p_curr = states[i]      # Current state

        # Petz retrodiction (Bayesian): given uniform prior, the retrodicted
        # previous state is proportional to p_curr (trivial case).
        # With empirical prior (average of all previous states up to t-1):
        if i == 1:
            prior = p_prev.copy()
        else:
            prior = np.mean(states[:i], axis=0)
            prior = prior / prior.sum()

        # Bayesian retrodiction: R(p_curr) ∝ prior * likelihood
        # Here likelihood ∝ p_curr (transition kernel is approximately identity
        # for consecutive epochs). The retrodicted state:
        retrodicted = prior * p_curr
        retrodicted = retrodicted / retrodicted.sum()

        # Classical fidelity (Bhattacharyya coefficient squared)
        bc = np.sum(np.sqrt(p_prev * retrodicted))
        fidelity = bc ** 2

        tau = 1.0 - fidelity
        taus.append(max(0.0, tau))

    return taus


def main():
    script_dir = Path(__file__).resolve().parent
    project_dir = script_dir.parent
    data_dir = project_dir / "data" / "sleep-edf"
    results_dir = project_dir / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    (results_dir / "figures").mkdir(exist_ok=True)
    (results_dir / "statistics").mkdir(exist_ok=True)

    print("=" * 70)
    print("Computing Σ (entropy production) and τ (retrodiction error)")
    print("Testing: τ ≤ 1 − exp(−Σ/2)  [Petz recovery bound]")
    print("=" * 70)

    # Find data files — match PSG and Hypnogram by subject prefix (first 6 chars)
    psg_files = sorted(data_dir.glob("*-PSG.edf"))
    hyp_files = sorted(data_dir.glob("*-Hypnogram.edf"))
    if not psg_files:
        print(f"No data in {data_dir}. Run download_data.py first.")
        sys.exit(1)

    # Build hypnogram lookup by subject prefix
    hyp_lookup = {}
    for h in hyp_files:
        prefix = h.name[:6]  # e.g., "SC4001"
        hyp_lookup[prefix] = h

    # Collect results: (sigma, tau, stage) tuples
    results_by_stage = {s: {"sigma": [], "tau": []} for s in STAGE_ORDER}

    for psg_path in psg_files:
        prefix = psg_path.name[:6]
        hyp_path = hyp_lookup.get(prefix)
        if hyp_path is None:
            continue

        subject_id = psg_path.stem.replace("-PSG", "")
        print(f"\nProcessing {subject_id}...")

        try:
            epochs, sfreq = load_eeg_and_stages(psg_path, hyp_path)
        except Exception as e:
            print(f"  ERROR: {e}")
            continue

        print(f"  Total epochs: {len(epochs)}")

        # Group consecutive epochs by stage
        # We need consecutive same-stage epochs to compute meaningful Σ and τ
        current_stage = None
        current_run = []

        def process_run(run, stage):
            """Process a run of consecutive same-stage epochs."""
            if len(run) < 3:
                return  # Need at least 3 for meaningful statistics

            states = [signal_to_state_vector(sig, sfreq) for sig, _, _ in run]
            sigmas = compute_sigma_forward_backward(states)
            taus = compute_tau_retrodiction(states)

            for s, t in zip(sigmas, taus):
                results_by_stage[stage]["sigma"].append(s)
                results_by_stage[stage]["tau"].append(t)

        for signal, stage, t_start in epochs:
            if stage == current_stage:
                current_run.append((signal, stage, t_start))
            else:
                if current_run and current_stage:
                    process_run(current_run, current_stage)
                current_stage = stage
                current_run = [(signal, stage, t_start)]

        # Process last run
        if current_run and current_stage:
            process_run(current_run, current_stage)

        for stage in STAGE_ORDER:
            n = len(results_by_stage[stage]["sigma"])
            print(f"  {stage:6s}: {n} transition pairs")

    # --- Summary ---
    print(f"\n{'=' * 70}")
    print("RESULTS: Σ and τ by consciousness state")
    print(f"{'=' * 70}")
    print(f"{'Stage':>8s}  {'N':>6s}  {'Σ_mean':>10s}  {'Σ_std':>10s}  "
          f"{'τ_mean':>10s}  {'τ_std':>10s}  {'bound':>10s}  {'holds?':>8s}")
    print("-" * 80)

    for stage in STAGE_ORDER:
        sigmas = results_by_stage[stage]["sigma"]
        taus = results_by_stage[stage]["tau"]
        if not sigmas:
            print(f"{stage:>8s}  {'N/A':>6s}")
            continue

        s_mean = np.mean(sigmas)
        s_std = np.std(sigmas)
        t_mean = np.mean(taus)
        t_std = np.std(taus)
        bound = 1.0 - np.exp(-s_mean / 2.0)

        # Check if bound holds (on average)
        holds = "YES" if t_mean <= bound + 0.01 else "NEAR" if t_mean <= bound + 0.05 else "NO"

        print(f"{stage:>8s}  {len(sigmas):6d}  {s_mean:10.4f}  {s_std:10.4f}  "
              f"{t_mean:10.4f}  {t_std:10.4f}  {bound:10.4f}  {holds:>8s}")

    # --- Save results ---
    import csv
    csv_path = results_dir / "statistics" / "sigma_tau_by_stage.csv"
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["stage", "sigma", "tau"])
        for stage in STAGE_ORDER:
            for s, t in zip(results_by_stage[stage]["sigma"],
                           results_by_stage[stage]["tau"]):
                writer.writerow([stage, f"{s:.6f}", f"{t:.6f}"])
    print(f"\nSaved raw data: {csv_path}")

    # --- Plot ---
    plot_sigma_tau(results_by_stage, results_dir)

    print("\nDone.")


def plot_sigma_tau(results_by_stage: dict, results_dir: Path):
    """Plot τ vs Σ with the Petz bound overlay."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(16, 7))

    # --- Panel 1: Scatter plot τ vs Σ ---
    ax = axes[0]

    sigma_range = np.linspace(0, 2.0, 200)
    bound_curve = 1.0 - np.exp(-sigma_range / 2.0)
    ax.plot(sigma_range, bound_curve, "k--", linewidth=2,
            label=r"$\tau = 1 - e^{-\Sigma/2}$ (Petz bound)", zorder=10)
    ax.fill_between(sigma_range, bound_curve, 1.0, alpha=0.1, color="red",
                     label="Forbidden region")

    for stage in STAGE_ORDER:
        sigmas = results_by_stage[stage]["sigma"]
        taus = results_by_stage[stage]["tau"]
        if not sigmas:
            continue
        ax.scatter(sigmas, taus, c=STAGE_COLORS[stage], label=stage,
                   alpha=0.3, s=10, edgecolors="none")

    ax.set_xlabel(r"$\Sigma$ (relative entropy production)", fontsize=12)
    ax.set_ylabel(r"$\tau$ (retrodiction error)", fontsize=12)
    ax.set_title(r"Testing $\tau \leq 1 - e^{-\Sigma/2}$ across consciousness states",
                 fontsize=13)
    ax.legend(fontsize=9, loc="lower right")
    ax.set_xlim(0, None)
    ax.set_ylim(0, 1)
    ax.grid(True, alpha=0.3)

    # --- Panel 2: Mean ± std per stage ---
    ax = axes[1]

    for stage in STAGE_ORDER:
        sigmas = results_by_stage[stage]["sigma"]
        taus = results_by_stage[stage]["tau"]
        if not sigmas:
            continue
        s_mean = np.mean(sigmas)
        s_std = np.std(sigmas)
        t_mean = np.mean(taus)
        t_std = np.std(taus)
        ax.errorbar(s_mean, t_mean, xerr=s_std, yerr=t_std,
                     fmt="o", markersize=12, color=STAGE_COLORS[stage],
                     label=stage, capsize=5, linewidth=2, zorder=5)

    ax.plot(sigma_range, bound_curve, "k--", linewidth=2, zorder=10)
    ax.fill_between(sigma_range, bound_curve, 1.0, alpha=0.1, color="red")

    ax.set_xlabel(r"$\langle\Sigma\rangle$", fontsize=12)
    ax.set_ylabel(r"$\langle\tau\rangle$", fontsize=12)
    ax.set_title("Mean retrodiction metrics by consciousness state", fontsize=13)
    ax.legend(fontsize=10)
    ax.set_xlim(0, None)
    ax.set_ylim(0, 1)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    fig_path = results_dir / "figures" / "sigma_tau_petz_bound.png"
    fig.savefig(str(fig_path), dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {fig_path}")


if __name__ == "__main__":
    main()
