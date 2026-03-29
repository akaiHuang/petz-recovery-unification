#!/usr/bin/env python3
"""
Load Sleep-EDF EEG data and generate preview plots.

This script:
1. Loads .edf files using MNE-Python
2. Parses hypnogram annotations (sleep stage labels)
3. Extracts EEG channels (Fpz-Cz, Pz-Oz)
4. Plots 30-second samples for each sleep stage
5. Generates summary statistics

Usage:
    python scripts/load_and_preview.py
"""

import os
import sys
import warnings
from pathlib import Path

import numpy as np

# Suppress MNE info messages for cleaner output
os.environ["MNE_LOGGING_LEVEL"] = "WARNING"

import mne

# -------------------------------------------------------------------
# Constants
# -------------------------------------------------------------------
STAGE_MAP = {
    "Sleep stage W": "Wake",
    "Sleep stage 1": "N1",
    "Sleep stage 2": "N2",
    "Sleep stage 3": "N3",
    "Sleep stage 4": "N3",   # N3 + N4 merged per AASM rules
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

EPOCH_DURATION = 30.0  # seconds (standard for sleep scoring)
EEG_CHANNELS = ["EEG Fpz-Cz", "EEG Pz-Oz"]


def find_data_files(data_dir: Path) -> list[tuple[Path, Path]]:
    """
    Find paired PSG + Hypnogram files.

    Sleep-EDF naming convention:
      PSG:       SC4SSNE0-PSG.edf       (SS=subject, N=night)
      Hypnogram: SC4SSNXX-Hypnogram.edf (XX varies: EC, EH, EJ, EP, etc.)

    We match by the first 6 characters (e.g., SC4001).
    """
    psg_files = sorted(data_dir.glob("*-PSG.edf"))
    hyp_files = sorted(data_dir.glob("*-Hypnogram.edf"))

    # Build lookup: subject prefix (first 6 chars) -> hypnogram path
    hyp_lookup = {}
    for hyp in hyp_files:
        prefix = hyp.name[:6]  # e.g., "SC4001"
        hyp_lookup[prefix] = hyp

    pairs = []
    for psg in psg_files:
        prefix = psg.name[:6]
        if prefix in hyp_lookup:
            pairs.append((psg, hyp_lookup[prefix]))
        else:
            print(f"  Warning: no hypnogram for {psg.name} (prefix={prefix})")
    return pairs


def load_eeg_and_hypnogram(psg_path: Path, hyp_path: Path):
    """
    Load EEG data and sleep stage annotations.

    Returns:
        raw: MNE Raw object (EEG channels only)
        epochs_by_stage: dict mapping stage name -> list of (start_sec, duration_sec)
    """
    # Load PSG
    raw = mne.io.read_raw_edf(str(psg_path), preload=True, verbose=False)

    # Pick EEG channels
    available_channels = raw.ch_names
    eeg_picks = [ch for ch in EEG_CHANNELS if ch in available_channels]
    if not eeg_picks:
        # Try partial matching
        eeg_picks = [ch for ch in available_channels if "EEG" in ch]
    if eeg_picks:
        raw.pick(eeg_picks)

    # Band-pass filter: 0.5-45 Hz
    raw.filter(0.5, 45.0, verbose=False)

    # Load hypnogram annotations
    annot = mne.read_annotations(str(hyp_path))
    raw.set_annotations(annot)

    # Parse annotations into epochs by stage
    epochs_by_stage = {stage: [] for stage in STAGE_ORDER}

    for ann in raw.annotations:
        stage_name = STAGE_MAP.get(ann["description"], None)
        if stage_name and stage_name in STAGE_ORDER:
            onset = ann["onset"]  # in seconds from start
            duration = ann["duration"]
            # Split long annotations into 30-second epochs
            n_epochs = int(duration // EPOCH_DURATION)
            for i in range(n_epochs):
                t_start = onset + i * EPOCH_DURATION
                epochs_by_stage[stage_name].append((t_start, EPOCH_DURATION))

    return raw, epochs_by_stage


def compute_epoch_features(raw, t_start: float, duration: float) -> dict:
    """
    Compute features for a single 30-second epoch.

    Returns dict with:
      - power_delta, power_theta, power_alpha, power_beta, power_gamma
      - spectral_entropy
      - rms_amplitude
    """
    sfreq = raw.info["sfreq"]
    start_sample = int(t_start * sfreq)
    n_samples = int(duration * sfreq)

    # Bounds check
    if start_sample < 0 or start_sample + n_samples > raw.n_times:
        return None

    data = raw.get_data(start=start_sample, stop=start_sample + n_samples)

    # Average across channels
    signal = data.mean(axis=0)

    # RMS amplitude
    rms = np.sqrt(np.mean(signal ** 2))

    # Power spectral density via Welch
    from scipy.signal import welch
    freqs, psd = welch(signal, fs=sfreq, nperseg=min(256, len(signal)))

    # Band powers
    def band_power(fmin, fmax):
        idx = np.where((freqs >= fmin) & (freqs <= fmax))[0]
        if len(idx) == 0:
            return 0.0
        return np.trapezoid(psd[idx], freqs[idx])

    total_power = np.trapezoid(psd, freqs)
    if total_power == 0:
        total_power = 1e-30

    power_delta = band_power(0.5, 4.0)
    power_theta = band_power(4.0, 8.0)
    power_alpha = band_power(8.0, 13.0)
    power_beta = band_power(13.0, 30.0)
    power_gamma = band_power(30.0, 45.0)

    # Normalized PSD for spectral entropy
    psd_norm = psd / psd.sum() if psd.sum() > 0 else psd
    psd_norm = psd_norm[psd_norm > 0]
    spectral_entropy = -np.sum(psd_norm * np.log2(psd_norm))

    return {
        "rms_uv": rms * 1e6,  # Convert to microvolts
        "power_delta": power_delta / total_power,
        "power_theta": power_theta / total_power,
        "power_alpha": power_alpha / total_power,
        "power_beta": power_beta / total_power,
        "power_gamma": power_gamma / total_power,
        "spectral_entropy": spectral_entropy,
        "total_power": total_power,
    }


def plot_stage_samples(raw, epochs_by_stage, subject_id: str, results_dir: Path):
    """Plot 30-second EEG samples for each sleep stage."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(len(STAGE_ORDER), 1, figsize=(14, 12), sharex=False)
    fig.suptitle(
        f"EEG Samples by Sleep Stage — {subject_id}\n"
        f"Testing τ ≤ 1−exp(−Σ/2) across consciousness states",
        fontsize=13, fontweight="bold"
    )

    sfreq = raw.info["sfreq"]

    for i, stage in enumerate(STAGE_ORDER):
        ax = axes[i]
        epochs = epochs_by_stage[stage]

        if not epochs:
            ax.text(0.5, 0.5, f"{stage}: no data", transform=ax.transAxes,
                    ha="center", va="center", fontsize=12)
            ax.set_ylabel(stage)
            continue

        # Pick a representative epoch (middle of the recording for that stage)
        mid_idx = len(epochs) // 2
        t_start, duration = epochs[mid_idx]

        start_sample = int(t_start * sfreq)
        n_samples = int(duration * sfreq)
        if start_sample + n_samples > raw.n_times:
            n_samples = raw.n_times - start_sample

        data = raw.get_data(start=start_sample, stop=start_sample + n_samples)
        times = np.arange(n_samples) / sfreq

        # Plot first channel
        signal_uv = data[0] * 1e6  # Convert to microvolts
        ax.plot(times, signal_uv, color=STAGE_COLORS[stage], linewidth=0.5, alpha=0.8)
        ax.set_ylabel(f"{stage}\n(μV)", fontsize=10)
        ax.set_xlim(0, duration)

        # Add epoch count
        ax.text(0.98, 0.95, f"n={len(epochs)} epochs",
                transform=ax.transAxes, ha="right", va="top",
                fontsize=9, bbox=dict(boxstyle="round,pad=0.3",
                                       facecolor=STAGE_COLORS[stage],
                                       alpha=0.2))

        if i < len(STAGE_ORDER) - 1:
            ax.set_xticklabels([])

    axes[-1].set_xlabel("Time (seconds)", fontsize=11)
    plt.tight_layout()

    fig_path = results_dir / "figures" / f"stage_samples_{subject_id}.png"
    fig_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(fig_path), dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {fig_path}")
    return fig_path


def plot_hypnogram(epochs_by_stage, subject_id: str, results_dir: Path):
    """Plot hypnogram (sleep stage vs time)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    stage_to_num = {"Wake": 4, "REM": 3, "N1": 2, "N2": 1, "N3": 0}

    all_epochs = []
    for stage, epochs in epochs_by_stage.items():
        for t_start, dur in epochs:
            all_epochs.append((t_start, stage))

    if not all_epochs:
        return None

    all_epochs.sort(key=lambda x: x[0])
    times_h = [e[0] / 3600.0 for e in all_epochs]
    stages_num = [stage_to_num.get(e[1], -1) for e in all_epochs]

    fig, ax = plt.subplots(figsize=(14, 3))
    ax.step(times_h, stages_num, where="post", color="#2c3e50", linewidth=1)
    ax.fill_between(times_h, stages_num, step="post", alpha=0.3, color="#3498db")
    ax.set_yticks(list(stage_to_num.values()))
    ax.set_yticklabels(list(stage_to_num.keys()))
    ax.set_xlabel("Time (hours)")
    ax.set_title(f"Hypnogram — {subject_id}")
    ax.invert_yaxis()
    plt.tight_layout()

    fig_path = results_dir / "figures" / f"hypnogram_{subject_id}.png"
    fig_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(fig_path), dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {fig_path}")
    return fig_path


def plot_spectral_summary(all_features: dict, results_dir: Path):
    """
    Plot average spectral features per stage across all subjects.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    fig.suptitle(
        "EEG Features by Sleep Stage (all subjects)\n"
        "Consciousness-Retrodiction Experiment",
        fontsize=13, fontweight="bold"
    )

    bands = ["power_delta", "power_theta", "power_alpha", "power_beta", "power_gamma"]
    band_labels = ["Delta\n(0.5-4 Hz)", "Theta\n(4-8 Hz)", "Alpha\n(8-13 Hz)",
                    "Beta\n(13-30 Hz)", "Gamma\n(30-45 Hz)"]

    # --- Panel 1: Band power by stage ---
    ax = axes[0]
    x = np.arange(len(bands))
    width = 0.15
    for j, stage in enumerate(STAGE_ORDER):
        if stage not in all_features or not all_features[stage]:
            continue
        means = [np.mean([f[b] for f in all_features[stage]]) for b in bands]
        ax.bar(x + j * width, means, width, label=stage,
               color=STAGE_COLORS[stage], alpha=0.8)
    ax.set_xticks(x + width * 2)
    ax.set_xticklabels(band_labels, fontsize=8)
    ax.set_ylabel("Relative Power")
    ax.set_title("Band Power Distribution")
    ax.legend(fontsize=8)

    # --- Panel 2: Spectral entropy by stage ---
    ax = axes[1]
    entropies = {}
    for stage in STAGE_ORDER:
        if stage in all_features and all_features[stage]:
            vals = [f["spectral_entropy"] for f in all_features[stage]]
            entropies[stage] = vals

    if entropies:
        positions = range(len(entropies))
        bp = ax.boxplot(
            [entropies[s] for s in STAGE_ORDER if s in entropies],
            labels=[s for s in STAGE_ORDER if s in entropies],
            patch_artist=True
        )
        for patch, stage in zip(bp["boxes"], [s for s in STAGE_ORDER if s in entropies]):
            patch.set_facecolor(STAGE_COLORS[stage])
            patch.set_alpha(0.6)
    ax.set_ylabel("Spectral Entropy (bits)")
    ax.set_title("Spectral Entropy\n(proxy for Σ)")

    # --- Panel 3: RMS amplitude by stage ---
    ax = axes[2]
    rms_vals = {}
    for stage in STAGE_ORDER:
        if stage in all_features and all_features[stage]:
            vals = [f["rms_uv"] for f in all_features[stage]]
            rms_vals[stage] = vals

    if rms_vals:
        bp = ax.boxplot(
            [rms_vals[s] for s in STAGE_ORDER if s in rms_vals],
            labels=[s for s in STAGE_ORDER if s in rms_vals],
            patch_artist=True
        )
        for patch, stage in zip(bp["boxes"], [s for s in STAGE_ORDER if s in rms_vals]):
            patch.set_facecolor(STAGE_COLORS[stage])
            patch.set_alpha(0.6)
    ax.set_ylabel("RMS Amplitude (μV)")
    ax.set_title("Signal Amplitude")

    plt.tight_layout()
    fig_path = results_dir / "figures" / "spectral_summary_all_subjects.png"
    fig_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(fig_path), dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {fig_path}")


def save_statistics(all_features: dict, results_dir: Path):
    """Save summary statistics as CSV."""
    import csv

    stats_dir = results_dir / "statistics"
    stats_dir.mkdir(parents=True, exist_ok=True)

    csv_path = stats_dir / "stage_summary_statistics.csv"
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "Stage", "N_epochs",
            "Spectral_Entropy_mean", "Spectral_Entropy_std",
            "RMS_uV_mean", "RMS_uV_std",
            "Delta_mean", "Theta_mean", "Alpha_mean", "Beta_mean", "Gamma_mean",
        ])

        for stage in STAGE_ORDER:
            if stage not in all_features or not all_features[stage]:
                continue
            feats = all_features[stage]
            n = len(feats)
            se_mean = np.mean([f["spectral_entropy"] for f in feats])
            se_std = np.std([f["spectral_entropy"] for f in feats])
            rms_mean = np.mean([f["rms_uv"] for f in feats])
            rms_std = np.std([f["rms_uv"] for f in feats])
            d_mean = np.mean([f["power_delta"] for f in feats])
            t_mean = np.mean([f["power_theta"] for f in feats])
            a_mean = np.mean([f["power_alpha"] for f in feats])
            b_mean = np.mean([f["power_beta"] for f in feats])
            g_mean = np.mean([f["power_gamma"] for f in feats])

            writer.writerow([
                stage, n,
                f"{se_mean:.4f}", f"{se_std:.4f}",
                f"{rms_mean:.2f}", f"{rms_std:.2f}",
                f"{d_mean:.4f}", f"{t_mean:.4f}", f"{a_mean:.4f}",
                f"{b_mean:.4f}", f"{g_mean:.4f}",
            ])

    print(f"  Saved: {csv_path}")


def main():
    # Paths
    script_dir = Path(__file__).resolve().parent
    project_dir = script_dir.parent
    data_dir = project_dir / "data" / "sleep-edf"
    results_dir = project_dir / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    (results_dir / "figures").mkdir(exist_ok=True)
    (results_dir / "statistics").mkdir(exist_ok=True)

    print("=" * 70)
    print("EEG Consciousness-Retrodiction Experiment")
    print("Load and Preview Sleep-EDF Data")
    print("=" * 70)

    # Find data files
    pairs = find_data_files(data_dir)
    if not pairs:
        print(f"\nNo data files found in {data_dir}")
        print("Run: python scripts/download_data.py")
        sys.exit(1)

    print(f"\nFound {len(pairs)} subject recordings.")

    # Aggregate features across all subjects
    all_features = {stage: [] for stage in STAGE_ORDER}

    for psg_path, hyp_path in pairs:
        subject_id = psg_path.stem.replace("-PSG", "")
        print(f"\n{'─' * 50}")
        print(f"Loading {subject_id}...")

        try:
            raw, epochs_by_stage = load_eeg_and_hypnogram(psg_path, hyp_path)
        except Exception as e:
            print(f"  ERROR loading: {e}")
            continue

        # Print summary
        print(f"  Channels: {raw.ch_names}")
        print(f"  Sampling rate: {raw.info['sfreq']} Hz")
        print(f"  Duration: {raw.n_times / raw.info['sfreq'] / 3600:.1f} hours")
        print(f"  Epochs per stage:")
        for stage in STAGE_ORDER:
            n = len(epochs_by_stage[stage])
            print(f"    {stage:6s}: {n:5d} epochs ({n * 30 / 60:.0f} min)")

        # Plot hypnogram
        plot_hypnogram(epochs_by_stage, subject_id, results_dir)

        # Plot stage samples
        plot_stage_samples(raw, epochs_by_stage, subject_id, results_dir)

        # Compute features (sample up to 100 epochs per stage for speed)
        print(f"  Computing spectral features...")
        max_epochs_per_stage = 100
        for stage in STAGE_ORDER:
            epochs = epochs_by_stage[stage]
            if not epochs:
                continue
            # Sample evenly
            indices = np.linspace(0, len(epochs) - 1,
                                   min(max_epochs_per_stage, len(epochs)),
                                   dtype=int)
            for idx in indices:
                t_start, dur = epochs[idx]
                feats = compute_epoch_features(raw, t_start, dur)
                if feats is not None:
                    all_features[stage].append(feats)

    # Summary across all subjects
    print(f"\n{'=' * 70}")
    print("AGGREGATE SUMMARY")
    print(f"{'=' * 70}")
    for stage in STAGE_ORDER:
        n = len(all_features[stage])
        if n == 0:
            print(f"  {stage:6s}: no data")
            continue
        se = np.mean([f["spectral_entropy"] for f in all_features[stage]])
        se_std = np.std([f["spectral_entropy"] for f in all_features[stage]])
        print(f"  {stage:6s}: {n:4d} epochs, "
              f"spectral_entropy = {se:.3f} +/- {se_std:.3f} bits")

    print()
    print("Interpretation for consciousness-retrodiction hypothesis:")
    print("  Higher spectral entropy ~ higher Sigma (entropy production)")
    print("  Expected ordering: N3 < N2 < N1 < REM < Wake")
    print("  This will feed into compute_sigma_tau.py for Petz bound testing.")

    # Plot and save
    print(f"\nGenerating summary plots...")
    plot_spectral_summary(all_features, results_dir)
    save_statistics(all_features, results_dir)

    print(f"\nDone. Results in: {results_dir}")


if __name__ == "__main__":
    main()
