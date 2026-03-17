#!/usr/bin/env python3
"""
Line 3 Experiment: Sigma-Guided 3D Gaussian Splatting.

Uses the Sigma analysis from sigma_analysis.py to determine an optimal
early stopping point, then re-runs OpenSplat with reduced iterations.

The key insight from Petz recovery theory: when Sigma (the relative entropy
between consecutive states) drops to near-zero, the optimization has reached
a retrodictable regime -- further steps are informationally redundant.

Usage:
    python run_sigma_guided.py                              # Auto-detect from sigma_profile.json
    python run_sigma_guided.py --stop-at 15000              # Manual early stop
    python run_sigma_guided.py --threshold 0.0005           # Custom Sigma threshold
    python run_sigma_guided.py --compare                    # Run and compare with baseline
"""

import argparse
import json
import shutil
import sys
import time
from pathlib import Path

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[3]
RESULTS_DIR = PROJECT_ROOT / "sigma_diffusion" / "results" / "line3"
CHECKPOINT_DIR = RESULTS_DIR / "checkpoints"
GUIDED_DIR = RESULTS_DIR / "guided"


def load_sigma_profile(results_dir: Path) -> dict | None:
    """Load sigma_profile.json from previous analysis."""
    profile_path = results_dir / "sigma_profile.json"
    if not profile_path.exists():
        return None
    with open(profile_path) as f:
        return json.load(f)


def determine_early_stop(
    profile: dict | None,
    manual_stop: int | None = None,
    threshold: float | None = None,
) -> int:
    """
    Determine the optimal early stopping iteration.

    Priority:
    1. Manual override (--stop-at)
    2. Sigma threshold search (--threshold)
    3. Suggested stop from sigma_profile.json
    4. Default: 50% of total iterations
    """
    if manual_stop is not None:
        print(f"[INFO] Using manual early stop: iteration {manual_stop}")
        return manual_stop

    if profile is None:
        print("[WARN] No sigma_profile.json found. Run sigma_analysis.py first.")
        print("[WARN] Falling back to default 50% early stop.")
        return 15000  # default for 30k training

    # Try threshold-based detection
    if threshold is not None:
        records = profile.get("sigma_records", [])
        for record in records:
            if record["sigma_mean"] < threshold:
                stop = record["iteration_from"]
                print(f"[INFO] Sigma drops below {threshold} at iteration {stop}")
                return stop
        print(f"[WARN] Sigma never drops below {threshold}. Using suggested stop.")

    # Use suggested stop from analysis
    summary = profile.get("summary", {})
    suggested = summary.get("suggested_early_stop")
    if suggested is not None:
        print(f"[INFO] Using suggested early stop from Sigma analysis: iteration {suggested}")
        return suggested

    # Fallback: find the iteration where Sigma is at 10% of its peak
    records = profile.get("sigma_records", [])
    if records:
        sigma_values = [r["sigma_mean"] for r in records]
        peak_sigma = max(sigma_values)
        target = peak_sigma * 0.1
        for record in records:
            if record["sigma_mean"] < target and record["iteration_to"] > records[0]["iteration_to"]:
                stop = record["iteration_from"]
                print(f"[INFO] Sigma reaches 10% of peak ({target:.6f}) at iteration {stop}")
                return stop

    # Last resort
    print("[WARN] Could not determine early stop. Using 50% default.")
    return 15000


def run_guided_training(
    data_dir: Path,
    stop_iteration: int,
    output_dir: Path,
    opensplat_bin: Path | None = None,
) -> dict:
    """Run OpenSplat with the determined early stop iteration."""
    # Import from run_baseline
    sys.path.insert(0, str(SCRIPT_DIR))
    from run_baseline import find_opensplat, print_install_instructions, run_opensplat

    if opensplat_bin is None:
        opensplat_bin = find_opensplat()
        if opensplat_bin is None:
            print_install_instructions()
            sys.exit(1)

    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n{'='*72}")
    print(f"  Line 3: Sigma-Guided 3DGS Training")
    print(f"{'='*72}")
    print(f"  Early stop iteration: {stop_iteration}")
    print(f"  Data directory:       {data_dir}")
    print(f"  Output directory:     {output_dir}")
    print(f"{'='*72}\n")

    # Run with reduced iterations
    # Use fewer checkpoints for the guided run (just start, middle, end)
    checkpoint_every = max(500, stop_iteration // 10)

    metadata = run_opensplat(
        opensplat_bin=opensplat_bin,
        data_dir=data_dir,
        iterations=stop_iteration,
        checkpoint_every=checkpoint_every,
        output_dir=output_dir / "checkpoints",
    )

    # Save guided training info
    guided_info = {
        "type": "sigma_guided",
        "stop_iteration": stop_iteration,
        "training_metadata": metadata,
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
    }

    info_path = output_dir / "guided_info.json"
    with open(info_path, "w") as f:
        json.dump(guided_info, f, indent=2)

    return guided_info


def compare_results(
    baseline_dir: Path,
    guided_dir: Path,
    output_dir: Path,
):
    """Quick comparison of baseline vs guided training."""
    print(f"\n{'='*72}")
    print(f"  Comparison: Baseline vs Sigma-Guided")
    print(f"{'='*72}")

    # Load metadata
    baseline_meta_path = baseline_dir / "training_metadata.json"
    guided_meta_path = guided_dir / "guided_info.json"

    baseline_meta = None
    guided_meta = None

    if baseline_meta_path.exists():
        with open(baseline_meta_path) as f:
            baseline_meta = json.load(f)
    else:
        print("[WARN] Baseline metadata not found.")

    if guided_meta_path.exists():
        with open(guided_meta_path) as f:
            guided_meta = json.load(f)
    else:
        print("[WARN] Guided metadata not found.")

    comparison = {
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "baseline": {},
        "guided": {},
        "savings": {},
    }

    if baseline_meta:
        baseline_time = baseline_meta.get("elapsed_seconds", 0)
        baseline_iters = baseline_meta.get("iterations", 30000)
        baseline_losses = baseline_meta.get("losses", [])
        baseline_final_loss = baseline_losses[-1]["loss"] if baseline_losses else None

        comparison["baseline"] = {
            "iterations": baseline_iters,
            "elapsed_seconds": baseline_time,
            "final_loss": baseline_final_loss,
            "num_checkpoints": len(baseline_meta.get("checkpoints", [])),
        }

        print(f"\n  Baseline:")
        print(f"    Iterations:   {baseline_iters}")
        print(f"    Time:         {baseline_time:.1f}s")
        if baseline_final_loss:
            print(f"    Final loss:   {baseline_final_loss:.6f}")

    if guided_meta:
        gm = guided_meta.get("training_metadata", {})
        guided_time = gm.get("elapsed_seconds", 0)
        guided_iters = guided_meta.get("stop_iteration", 0)
        guided_losses = gm.get("losses", [])
        guided_final_loss = guided_losses[-1]["loss"] if guided_losses else None

        comparison["guided"] = {
            "iterations": guided_iters,
            "elapsed_seconds": guided_time,
            "final_loss": guided_final_loss,
            "stop_iteration": guided_iters,
        }

        print(f"\n  Sigma-Guided:")
        print(f"    Iterations:   {guided_iters}")
        print(f"    Time:         {guided_time:.1f}s")
        if guided_final_loss:
            print(f"    Final loss:   {guided_final_loss:.6f}")

    if baseline_meta and guided_meta:
        baseline_time = comparison["baseline"]["elapsed_seconds"]
        guided_time = comparison["guided"]["elapsed_seconds"]
        baseline_iters = comparison["baseline"]["iterations"]
        guided_iters = comparison["guided"]["iterations"]

        time_savings = (1 - guided_time / baseline_time) * 100 if baseline_time > 0 else 0
        iter_savings = (1 - guided_iters / baseline_iters) * 100 if baseline_iters > 0 else 0

        comparison["savings"] = {
            "iteration_reduction_pct": iter_savings,
            "time_reduction_pct": time_savings,
            "iterations_saved": baseline_iters - guided_iters,
            "time_saved_seconds": baseline_time - guided_time,
        }

        print(f"\n  Savings:")
        print(f"    Iteration reduction: {iter_savings:.1f}%")
        print(f"    Time reduction:      {time_savings:.1f}%")
        print(f"    Iterations saved:    {baseline_iters - guided_iters}")
        print(f"    Time saved:          {baseline_time - guided_time:.1f}s")

    # Note about quality comparison
    print(f"\n  [NOTE] For quality comparison (PSNR/SSIM), run evaluate.py")

    print(f"{'='*72}")

    # Save comparison
    comp_path = output_dir / "comparison_summary.json"
    with open(comp_path, "w") as f:
        json.dump(comparison, f, indent=2)
    print(f"\n[INFO] Comparison saved to: {comp_path}")

    return comparison


def find_data_dir() -> Path | None:
    """Find the default data directory (same logic as run_baseline)."""
    candidates = [
        SCRIPT_DIR / "data" / "minimal_test",
        SCRIPT_DIR / "data" / "garden",
        SCRIPT_DIR / "data" / "lego",
        SCRIPT_DIR / "data" / "custom",
    ]
    for c in candidates:
        if c.is_dir() and (c / "images").is_dir():
            return c
    return None


def main():
    parser = argparse.ArgumentParser(
        description="Sigma-guided 3DGS training with early stopping"
    )
    parser.add_argument(
        "--data-dir", type=str, default=None,
        help="Path to COLMAP-formatted scene directory"
    )
    parser.add_argument(
        "--stop-at", type=int, default=None,
        help="Manual early stop iteration (overrides Sigma analysis)"
    )
    parser.add_argument(
        "--threshold", type=float, default=None,
        help="Sigma threshold for early stop detection"
    )
    parser.add_argument(
        "--output-dir", type=str, default=str(GUIDED_DIR),
        help="Output directory for guided training results"
    )
    parser.add_argument(
        "--compare", action="store_true",
        help="Also compare with baseline after training"
    )
    parser.add_argument(
        "--opensplat-bin", type=str, default=None,
        help="Path to OpenSplat binary"
    )
    args = parser.parse_args()

    # Load Sigma profile
    profile = load_sigma_profile(RESULTS_DIR)

    # Determine early stop
    stop_iteration = determine_early_stop(profile, args.stop_at, args.threshold)

    # Find data directory
    if args.data_dir:
        data_dir = Path(args.data_dir)
    else:
        data_dir = find_data_dir()

    if data_dir is None:
        # If we don't have data, we can still do the analysis part
        if profile is not None:
            print(f"\n[INFO] No data directory found, but Sigma profile is available.")
            print(f"[INFO] Suggested early stop: iteration {stop_iteration}")
            print(f"[INFO] To run guided training, provide --data-dir")
            sys.exit(0)
        else:
            print("[ERROR] No data directory found. Run prepare_data.sh first.")
            sys.exit(1)

    opensplat_bin = Path(args.opensplat_bin) if args.opensplat_bin else None
    output_dir = Path(args.output_dir)

    # Run guided training
    guided_info = run_guided_training(
        data_dir=data_dir,
        stop_iteration=stop_iteration,
        output_dir=output_dir,
        opensplat_bin=opensplat_bin,
    )

    # Compare if requested
    if args.compare:
        compare_results(CHECKPOINT_DIR, output_dir, RESULTS_DIR)


if __name__ == "__main__":
    main()
