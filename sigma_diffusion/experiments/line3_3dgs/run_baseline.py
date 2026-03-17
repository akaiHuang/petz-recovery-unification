#!/usr/bin/env python3
"""
Line 3 Experiment: Run OpenSplat baseline training with checkpoint capture.

Runs OpenSplat for 30,000 iterations, saving PLY checkpoints every 500 iterations.
Captures training loss from stdout for later Sigma analysis.

Usage:
    python run_baseline.py                          # Use minimal_test scene
    python run_baseline.py --data-dir data/garden   # Use specific scene
    python run_baseline.py --iterations 10000       # Fewer iterations
    python run_baseline.py --checkpoint-every 1000  # Less frequent checkpoints
"""

import argparse
import json
import os
import re
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
CHECKPOINT_DIR = RESULTS_DIR / "checkpoints"
OPENSPLAT_DEFAULT = Path("/Users/akaihuangm1/Desktop/github/OpenSplat/build/opensplat")

# Alternative locations to search for OpenSplat
OPENSPLAT_SEARCH_PATHS = [
    OPENSPLAT_DEFAULT,
    Path.home() / "OpenSplat" / "build" / "opensplat",
    Path("/usr/local/bin/opensplat"),
    Path("/opt/homebrew/bin/opensplat"),
]


def find_opensplat() -> Path | None:
    """Search for the OpenSplat binary in known locations."""
    # Check environment variable first
    env_path = os.environ.get("OPENSPLAT_BIN")
    if env_path:
        p = Path(env_path)
        if p.is_file() and os.access(p, os.X_OK):
            return p

    # Check PATH
    which_result = shutil.which("opensplat")
    if which_result:
        return Path(which_result)

    # Check known locations
    for path in OPENSPLAT_SEARCH_PATHS:
        if path.is_file() and os.access(path, os.X_OK):
            return path

    return None


def print_install_instructions():
    """Print clear instructions for building OpenSplat on Mac."""
    print("""
================================================================================
  OpenSplat not found!
================================================================================

OpenSplat binary not found at any of the expected locations:
  - /Users/akaihuangm1/Desktop/github/OpenSplat/build/opensplat
  - $OPENSPLAT_BIN environment variable
  - System PATH

To build OpenSplat on Mac (M1/M2/M3):

  1. Clone the repository:
     cd /Users/akaihuangm1/Desktop/github/
     git clone https://github.com/pierotofy/OpenSplat.git
     cd OpenSplat

  2. Install dependencies:
     brew install cmake libtorch
     # Or with conda:
     # conda install pytorch::pytorch torchvision -c pytorch

  3. Build:
     mkdir build && cd build
     cmake .. -DCMAKE_BUILD_TYPE=Release -DGPU_RUNTIME=MPS
     make -j$(sysctl -n hw.ncpu)

  4. Verify:
     ./opensplat --help

  5. Re-run this script:
     python run_baseline.py

Alternative: set OPENSPLAT_BIN environment variable to point to your binary:
     export OPENSPLAT_BIN=/path/to/opensplat
================================================================================
""")


def find_data_dir() -> Path | None:
    """Find the default data directory."""
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


def run_opensplat(
    opensplat_bin: Path,
    data_dir: Path,
    iterations: int = 30000,
    checkpoint_every: int = 500,
    output_dir: Path | None = None,
) -> dict:
    """
    Run OpenSplat with checkpoint saving.

    OpenSplat saves the final model as a .ply file. To get intermediate checkpoints,
    we run OpenSplat multiple times with increasing iteration counts, or use the
    --save-every flag if available.

    Returns a dict with training metadata.
    """
    if output_dir is None:
        output_dir = CHECKPOINT_DIR

    output_dir.mkdir(parents=True, exist_ok=True)

    metadata = {
        "data_dir": str(data_dir),
        "iterations": iterations,
        "checkpoint_every": checkpoint_every,
        "opensplat_bin": str(opensplat_bin),
        "start_time": time.time(),
        "checkpoints": [],
        "losses": [],
    }

    # Strategy: Run OpenSplat in segments to capture checkpoints.
    # OpenSplat supports --save-every for intermediate saves in some versions.
    # We try that first; if it fails, fall back to segmented runs.

    print(f"\n{'='*72}")
    print(f"  Line 3 Baseline: OpenSplat Training")
    print(f"{'='*72}")
    print(f"  Data directory:    {data_dir}")
    print(f"  Total iterations:  {iterations}")
    print(f"  Checkpoint every:  {checkpoint_every}")
    print(f"  Output directory:  {output_dir}")
    print(f"  OpenSplat binary:  {opensplat_bin}")
    print(f"{'='*72}\n")

    # Check OpenSplat capabilities
    help_output = subprocess.run(
        [str(opensplat_bin), "--help"],
        capture_output=True, text=True, timeout=10
    )
    help_text = help_output.stdout + help_output.stderr
    has_save_every = "--save-every" in help_text

    if has_save_every:
        print("[INFO] OpenSplat supports --save-every. Using native checkpointing.")
        success = _run_with_native_checkpoints(
            opensplat_bin, data_dir, iterations, checkpoint_every, output_dir, metadata
        )
    else:
        print("[INFO] OpenSplat does not support --save-every.")
        print("[INFO] Using segmented training to capture checkpoints.")
        success = _run_segmented(
            opensplat_bin, data_dir, iterations, checkpoint_every, output_dir, metadata
        )

    metadata["end_time"] = time.time()
    metadata["elapsed_seconds"] = metadata["end_time"] - metadata["start_time"]
    metadata["success"] = success

    # Save metadata
    meta_path = output_dir / "training_metadata.json"
    with open(meta_path, "w") as f:
        json.dump(metadata, f, indent=2)
    print(f"\n[INFO] Training metadata saved to: {meta_path}")

    return metadata


def _run_with_native_checkpoints(
    opensplat_bin, data_dir, iterations, checkpoint_every, output_dir, metadata
):
    """Run OpenSplat with --save-every flag for native checkpointing."""
    output_ply = output_dir / "final.ply"

    cmd = [
        str(opensplat_bin),
        str(data_dir),
        "--num-iters", str(iterations),
        "--save-every", str(checkpoint_every),
        "-o", str(output_ply),
    ]

    print(f"[CMD] {' '.join(cmd)}\n")

    return _execute_and_capture(cmd, output_dir, metadata)


def _run_segmented(
    opensplat_bin, data_dir, iterations, checkpoint_every, output_dir, metadata
):
    """Run OpenSplat multiple times with increasing iterations."""
    checkpoints = list(range(checkpoint_every, iterations + 1, checkpoint_every))
    if checkpoints[-1] != iterations:
        checkpoints.append(iterations)

    all_success = True

    for i, num_iters in enumerate(checkpoints):
        ply_name = f"checkpoint_{num_iters:06d}.ply"
        output_ply = output_dir / ply_name

        if output_ply.exists():
            print(f"[SKIP] Checkpoint already exists: {ply_name}")
            metadata["checkpoints"].append({
                "iteration": num_iters,
                "file": str(output_ply),
                "skipped": True,
            })
            continue

        print(f"\n[RUN {i+1}/{len(checkpoints)}] Training to iteration {num_iters}...")

        cmd = [
            str(opensplat_bin),
            str(data_dir),
            "--num-iters", str(num_iters),
            "-o", str(output_ply),
        ]

        print(f"[CMD] {' '.join(cmd)}")

        success = _execute_and_capture(cmd, output_dir, metadata, iteration=num_iters)
        if not success:
            print(f"[ERROR] Training failed at iteration {num_iters}")
            all_success = False
            break

        if output_ply.exists():
            size_mb = output_ply.stat().st_size / (1024 * 1024)
            print(f"[OK] Checkpoint saved: {ply_name} ({size_mb:.1f} MB)")
            metadata["checkpoints"].append({
                "iteration": num_iters,
                "file": str(output_ply),
                "size_mb": size_mb,
            })
        else:
            # OpenSplat may save with a different name
            print(f"[WARN] Expected output not found: {output_ply}")
            # Look for any new .ply files
            ply_files = list(output_dir.glob("*.ply"))
            if ply_files:
                latest = max(ply_files, key=lambda p: p.stat().st_mtime)
                target = output_dir / ply_name
                if latest != target:
                    shutil.copy2(latest, target)
                    print(f"[INFO] Renamed {latest.name} -> {ply_name}")

    return all_success


def _execute_and_capture(cmd, output_dir, metadata, iteration=None):
    """Execute a command, capture stdout, and parse loss values."""
    log_path = output_dir / "training.log"
    loss_pattern = re.compile(
        r"(?:iter|step|iteration)\s*[:#]?\s*(\d+).*?(?:loss|Loss|LOSS)\s*[=:]\s*([\d.eE+-]+)",
        re.IGNORECASE,
    )
    # Also try a simpler pattern for different log formats
    loss_pattern_simple = re.compile(
        r"(\d+)\s+.*?([\d.]+(?:e[+-]?\d+)?)\s*$",
        re.IGNORECASE,
    )

    try:
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )

        with open(log_path, "a") as log_file:
            for line in process.stdout:
                # Write to log
                log_file.write(line)
                # Print to terminal
                sys.stdout.write(line)
                sys.stdout.flush()

                # Try to parse loss
                match = loss_pattern.search(line)
                if not match:
                    match = loss_pattern_simple.search(line.strip())

                if match:
                    try:
                        iter_num = int(match.group(1))
                        loss_val = float(match.group(2))
                        metadata["losses"].append({
                            "iteration": iter_num,
                            "loss": loss_val,
                        })
                    except (ValueError, IndexError):
                        pass

        process.wait()
        return process.returncode == 0

    except subprocess.TimeoutExpired:
        process.kill()
        print("[ERROR] Process timed out")
        return False
    except Exception as e:
        print(f"[ERROR] Process failed: {e}")
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Run OpenSplat baseline for Line 3 Sigma experiment"
    )
    parser.add_argument(
        "--data-dir", type=str, default=None,
        help="Path to COLMAP-formatted scene directory"
    )
    parser.add_argument(
        "--iterations", type=int, default=30000,
        help="Total training iterations (default: 30000)"
    )
    parser.add_argument(
        "--checkpoint-every", type=int, default=500,
        help="Save checkpoint every N iterations (default: 500)"
    )
    parser.add_argument(
        "--output-dir", type=str, default=None,
        help="Output directory for checkpoints"
    )
    parser.add_argument(
        "--opensplat-bin", type=str, default=None,
        help="Path to OpenSplat binary"
    )
    args = parser.parse_args()

    # Find OpenSplat
    if args.opensplat_bin:
        opensplat = Path(args.opensplat_bin)
        if not opensplat.is_file():
            print(f"[ERROR] OpenSplat binary not found at: {opensplat}")
            sys.exit(1)
    else:
        opensplat = find_opensplat()
        if opensplat is None:
            print_install_instructions()
            sys.exit(1)

    print(f"[INFO] Using OpenSplat: {opensplat}")

    # Find data directory
    if args.data_dir:
        data_dir = Path(args.data_dir)
    else:
        data_dir = find_data_dir()
        if data_dir is None:
            print("[ERROR] No data directory found. Run prepare_data.sh first:")
            print("  cd", SCRIPT_DIR)
            print("  bash prepare_data.sh")
            sys.exit(1)

    if not data_dir.is_dir():
        print(f"[ERROR] Data directory not found: {data_dir}")
        sys.exit(1)

    if not (data_dir / "images").is_dir():
        print(f"[ERROR] No 'images' subdirectory in: {data_dir}")
        sys.exit(1)

    print(f"[INFO] Using data: {data_dir}")

    # Output directory
    output_dir = Path(args.output_dir) if args.output_dir else CHECKPOINT_DIR

    # Run training
    metadata = run_opensplat(
        opensplat_bin=opensplat,
        data_dir=data_dir,
        iterations=args.iterations,
        checkpoint_every=args.checkpoint_every,
        output_dir=output_dir,
    )

    # Summary
    print(f"\n{'='*72}")
    print(f"  Training Complete")
    print(f"{'='*72}")
    print(f"  Elapsed time:  {metadata['elapsed_seconds']:.1f}s")
    print(f"  Checkpoints:   {len(metadata['checkpoints'])}")
    print(f"  Loss records:  {len(metadata['losses'])}")
    print(f"  Success:       {metadata['success']}")
    if metadata['losses']:
        final_loss = metadata['losses'][-1]['loss']
        print(f"  Final loss:    {final_loss:.6f}")
    print(f"{'='*72}")


if __name__ == "__main__":
    main()
