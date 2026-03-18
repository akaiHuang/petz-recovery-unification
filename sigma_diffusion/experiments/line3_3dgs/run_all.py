#!/usr/bin/env python3
"""
Line 3 Experiment: Full Orchestrator

Runs the complete Line 3 experiment pipeline:
  1. prepare_data  -- Generate or download test scene
  2. baseline      -- Run OpenSplat for 30,000 iterations with checkpoints
  3. analysis      -- Compute Sigma profile from checkpoint PLY files
  4. guided        -- Re-run with Sigma-determined early stopping
  5. evaluate      -- Compare baseline vs guided (PSNR, SSIM, Sigma)

Usage:
    python run_all.py                           # Run everything
    python run_all.py --start-from analysis     # Skip data prep and baseline
    python run_all.py --data-dir data/garden    # Use specific scene
    python run_all.py --iterations 10000        # Quick test with fewer iterations
    python run_all.py --dry-run                 # Show what would be run

Each step can also be run independently:
    bash prepare_data.sh
    python run_baseline.py
    python sigma_analysis.py
    python run_sigma_guided.py
    python evaluate.py
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[3]
RESULTS_DIR = PROJECT_ROOT / "sigma_diffusion" / "results" / "line3"

STEPS = ["prepare", "baseline", "analysis", "guided", "evaluate"]


def banner(text: str, width: int = 72):
    """Print a banner."""
    print(f"\n{'#' * width}")
    print(f"#  {text}")
    print(f"{'#' * width}\n")


def check_prerequisites() -> dict:
    """Check what's available on this system."""
    status = {
        "python": sys.executable,
        "python_version": sys.version.split()[0],
        "opensplat": None,
        "colmap": None,
        "numpy": False,
        "matplotlib": False,
        "PIL": False,
        "scipy": False,
    }

    # Check OpenSplat
    sys.path.insert(0, str(SCRIPT_DIR))
    try:
        from run_baseline import find_opensplat
        opensplat = find_opensplat()
        if opensplat:
            status["opensplat"] = str(opensplat)
    except ImportError:
        pass

    # Check COLMAP
    if shutil.which("colmap"):
        status["colmap"] = shutil.which("colmap")

    # Check Python packages
    for pkg in ["numpy", "matplotlib", "PIL", "scipy"]:
        try:
            __import__(pkg)
            status[pkg] = True
        except ImportError:
            status[pkg] = False

    return status


def print_prerequisites(status: dict):
    """Print prerequisite check results."""
    banner("Prerequisite Check")

    def icon(ok):
        return "[OK]" if ok else "[--]"

    print(f"  {icon(True)} Python:      {status['python_version']} ({status['python']})")
    print(f"  {icon(status['opensplat'])} OpenSplat:   {status['opensplat'] or 'NOT FOUND'}")
    print(f"  {icon(status['colmap'])} COLMAP:      {status['colmap'] or 'not found (optional)'}")
    print(f"  {icon(status['numpy'])} NumPy:       {'installed' if status['numpy'] else 'NOT FOUND'}")
    print(f"  {icon(status['matplotlib'])} Matplotlib:  {'installed' if status['matplotlib'] else 'NOT FOUND (optional)'}")
    print(f"  {icon(status['PIL'])} Pillow:      {'installed' if status['PIL'] else 'NOT FOUND (optional)'}")
    print(f"  {icon(status['scipy'])} SciPy:       {'installed' if status['scipy'] else 'NOT FOUND (optional)'}")
    print()


def print_opensplat_instructions():
    """Print instructions for building OpenSplat."""
    print("""
================================================================================
  OpenSplat is required for the baseline and guided training steps.

  To build OpenSplat on Mac (M1/M2/M3):

    cd /Users/akaihuangm1/Desktop/github/
    git clone https://github.com/pierotofy/OpenSplat.git
    cd OpenSplat

    # Install dependencies
    brew install cmake
    pip install torch torchvision  # or use conda

    # Build
    mkdir build && cd build
    cmake .. -DCMAKE_BUILD_TYPE=Release -DGPU_RUNTIME=MPS
    make -j$(sysctl -n hw.ncpu)

    # Verify
    ./opensplat --help

  Then re-run: python run_all.py
================================================================================
""")


def install_missing_packages(status: dict):
    """Install missing Python packages."""
    missing = []
    if not status["numpy"]:
        missing.append("numpy")
    if not status["matplotlib"]:
        missing.append("matplotlib")
    if not status["PIL"]:
        missing.append("Pillow")
    if not status["scipy"]:
        missing.append("scipy")

    if missing:
        print(f"[INFO] Installing missing packages: {', '.join(missing)}")
        subprocess.run(
            [sys.executable, "-m", "pip", "install"] + missing,
            check=True,
        )
        print("[INFO] Packages installed successfully.")


# ---------------------------------------------------------------------------
# Step Runners
# ---------------------------------------------------------------------------
def step_prepare(args) -> bool:
    """Step 1: Prepare data."""
    banner("Step 1/5: Prepare Data")

    data_dir = SCRIPT_DIR / "data"
    if args.data_dir:
        data_dir = Path(args.data_dir)
        if data_dir.is_dir() and (data_dir / "images").is_dir():
            print(f"[INFO] Using existing data directory: {data_dir}")
            return True

    # Check if data already exists
    for scene in ["minimal_test", "garden", "lego", "custom"]:
        scene_dir = SCRIPT_DIR / "data" / scene
        if scene_dir.is_dir() and (scene_dir / "images").is_dir():
            print(f"[INFO] Found existing scene: {scene_dir}")
            return True

    # Run prepare_data.sh
    script = SCRIPT_DIR / "prepare_data.sh"
    if not script.exists():
        print("[ERROR] prepare_data.sh not found!")
        return False

    mode = args.scene_mode if hasattr(args, "scene_mode") else "--minimal"

    result = subprocess.run(
        ["bash", str(script), mode],
        cwd=str(SCRIPT_DIR),
    )
    return result.returncode == 0


def step_baseline(args) -> bool:
    """Step 2: Run baseline training."""
    banner("Step 2/5: Baseline Training (OpenSplat)")

    cmd = [sys.executable, str(SCRIPT_DIR / "run_baseline.py")]

    if args.data_dir:
        cmd.extend(["--data-dir", args.data_dir])
    if args.iterations:
        cmd.extend(["--iterations", str(args.iterations)])
    if args.checkpoint_every:
        cmd.extend(["--checkpoint-every", str(args.checkpoint_every)])
    if args.opensplat_bin:
        cmd.extend(["--opensplat-bin", args.opensplat_bin])

    result = subprocess.run(cmd, cwd=str(SCRIPT_DIR))
    return result.returncode == 0


def step_analysis(args) -> bool:
    """Step 3: Sigma analysis."""
    banner("Step 3/5: Sigma Analysis")

    cmd = [
        sys.executable, str(SCRIPT_DIR / "sigma_analysis.py"),
        "--threshold", str(args.sigma_threshold),
    ]

    result = subprocess.run(cmd, cwd=str(SCRIPT_DIR))
    return result.returncode == 0


def step_guided(args) -> bool:
    """Step 4: Sigma-guided training."""
    banner("Step 4/5: Sigma-Guided Training")

    cmd = [sys.executable, str(SCRIPT_DIR / "run_sigma_guided.py")]

    if args.data_dir:
        cmd.extend(["--data-dir", args.data_dir])
    if args.opensplat_bin:
        cmd.extend(["--opensplat-bin", args.opensplat_bin])

    result = subprocess.run(cmd, cwd=str(SCRIPT_DIR))
    return result.returncode == 0


def step_evaluate(args) -> bool:
    """Step 5: Evaluation."""
    banner("Step 5/5: Evaluation")

    cmd = [sys.executable, str(SCRIPT_DIR / "evaluate.py")]

    if args.data_dir:
        cmd.extend(["--data-dir", args.data_dir])

    result = subprocess.run(cmd, cwd=str(SCRIPT_DIR))
    return result.returncode == 0


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Line 3 Experiment: Full 3DGS Sigma Optimization Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python run_all.py                           # Full pipeline
  python run_all.py --start-from analysis     # Skip to analysis (data/baseline done)
  python run_all.py --iterations 5000         # Quick test
  python run_all.py --dry-run                 # Preview steps
  python run_all.py --skip baseline guided    # Only prepare + analyze + evaluate
"""
    )
    parser.add_argument(
        "--start-from", type=str, choices=STEPS, default="prepare",
        help="Start from this step (skip earlier steps)"
    )
    parser.add_argument(
        "--skip", nargs="+", choices=STEPS, default=[],
        help="Skip specific steps"
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Show what would be run without executing"
    )
    parser.add_argument(
        "--data-dir", type=str, default=None,
        help="Path to scene data directory"
    )
    parser.add_argument(
        "--scene-mode", type=str, default="--minimal",
        choices=["--minimal", "--lego", "--garden"],
        help="Scene to prepare (default: --minimal)"
    )
    parser.add_argument(
        "--iterations", type=int, default=30000,
        help="Total training iterations (default: 30000)"
    )
    parser.add_argument(
        "--checkpoint-every", type=int, default=500,
        help="Checkpoint interval (default: 500)"
    )
    parser.add_argument(
        "--sigma-threshold", type=float, default=0.001,
        help="Sigma threshold for plateau detection (default: 0.001)"
    )
    parser.add_argument(
        "--opensplat-bin", type=str, default=None,
        help="Path to OpenSplat binary"
    )
    parser.add_argument(
        "--install-deps", action="store_true",
        help="Auto-install missing Python dependencies"
    )
    args = parser.parse_args()

    # Header
    banner("Line 3 Experiment: 3D Gaussian Splatting + Sigma Optimization")
    print(f"  Project root:  {PROJECT_ROOT}")
    print(f"  Script dir:    {SCRIPT_DIR}")
    print(f"  Results dir:   {RESULTS_DIR}")
    print(f"  Iterations:    {args.iterations}")
    print(f"  Checkpoints:   every {args.checkpoint_every}")
    print(f"  Sigma threshold: {args.sigma_threshold}")
    print()

    # Check prerequisites
    status = check_prerequisites()
    print_prerequisites(status)

    # Install missing packages if requested
    if args.install_deps:
        install_missing_packages(status)

    # Determine which steps to run
    start_idx = STEPS.index(args.start_from)
    steps_to_run = [s for i, s in enumerate(STEPS) if i >= start_idx and s not in args.skip]

    # Check if OpenSplat is needed
    needs_opensplat = "baseline" in steps_to_run or "guided" in steps_to_run
    if needs_opensplat and not status["opensplat"]:
        print("[WARN] OpenSplat not found. Training steps will be skipped.")
        if "baseline" in steps_to_run:
            steps_to_run.remove("baseline")
        if "guided" in steps_to_run:
            steps_to_run.remove("guided")

        if not steps_to_run or steps_to_run == ["prepare"]:
            print_opensplat_instructions()
            sys.exit(1)

    # Dry run
    if args.dry_run:
        banner("Dry Run: Steps to Execute")
        for i, step in enumerate(steps_to_run):
            print(f"  [{i+1}] {step}")
        print("\n  (Use without --dry-run to execute)")
        sys.exit(0)

    # Execute steps
    step_functions = {
        "prepare": step_prepare,
        "baseline": step_baseline,
        "analysis": step_analysis,
        "guided": step_guided,
        "evaluate": step_evaluate,
    }

    results = {}
    overall_start = time.time()

    for step_name in steps_to_run:
        step_start = time.time()
        func = step_functions[step_name]

        try:
            success = func(args)
        except KeyboardInterrupt:
            print(f"\n[INTERRUPTED] Step '{step_name}' was interrupted.")
            success = False
        except Exception as e:
            print(f"\n[ERROR] Step '{step_name}' failed with exception: {e}")
            success = False

        elapsed = time.time() - step_start
        results[step_name] = {
            "success": success,
            "elapsed_seconds": elapsed,
        }

        if success:
            print(f"\n[OK] Step '{step_name}' completed in {elapsed:.1f}s")
        else:
            print(f"\n[FAILED] Step '{step_name}' failed after {elapsed:.1f}s")

            # Decide whether to continue
            if step_name in ("prepare", "baseline") and step_name != steps_to_run[-1]:
                print("[INFO] Cannot continue without this step. Stopping.")
                break
            else:
                print("[INFO] Continuing with remaining steps...")

    # Final summary
    overall_elapsed = time.time() - overall_start

    banner("Pipeline Summary")
    print(f"  Total time: {overall_elapsed:.1f}s\n")

    for step_name, result in results.items():
        status_str = "[OK]  " if result["success"] else "[FAIL]"
        print(f"  {status_str} {step_name:12s}  ({result['elapsed_seconds']:.1f}s)")

    # Save pipeline results
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    pipeline_path = RESULTS_DIR / "pipeline_results.json"
    with open(pipeline_path, "w") as f:
        json.dump({
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
            "total_elapsed_seconds": overall_elapsed,
            "steps": results,
            "args": vars(args),
        }, f, indent=2)
    print(f"\n  Pipeline results: {pipeline_path}")

    # Check overall success
    all_success = all(r["success"] for r in results.values())
    if all_success:
        print(f"\n  All steps completed successfully!")
        print(f"  Results directory: {RESULTS_DIR}")
        print(f"  Key outputs:")
        for f_name in ["sigma_profile.json", "sigma_heatmap.png", "metrics.json", "comparison.png"]:
            f_path = RESULTS_DIR / f_name
            if f_path.exists():
                print(f"    - {f_path}")
    else:
        print(f"\n  Some steps failed. Check output above for details.")
        sys.exit(1)


if __name__ == "__main__":
    main()
