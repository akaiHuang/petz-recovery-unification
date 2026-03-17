#!/usr/bin/env python3
"""
Line 1 Orchestrator: Run the full baseline -> sigma analysis -> guided -> evaluate pipeline.

This is the main entry point for the Line 1 experiment (Z-Image / FLUX
image diffusion with Sigma optimization).

Stages:
    1. Baseline generation -- full-step diffusion with per-step Sigma tracking
    2. Sigma analysis -- identify low-information steps from the Sigma profile
    3. Guided generation -- re-run with adaptive step-skipping
    4. Evaluation -- compare baseline vs guided (PSNR, SSIM, LPIPS, visual grid)

Usage:
    # Quick test (10 images)
    python run_all.py --num-images 10

    # Full run (50 images)
    python run_all.py --num-images 50

    # Custom settings
    python run_all.py --num-images 20 --num-steps 28 --device mps --dtype bfloat16

    # Resume from guided step (skip baseline)
    python run_all.py --num-images 10 --skip-baseline

    # Only evaluate existing results
    python run_all.py --evaluate-only
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import sys
import time
from pathlib import Path
from typing import Any, Dict, Optional

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
logger = logging.getLogger(__name__)


def check_dependencies() -> Dict[str, bool]:
    """Check which dependencies are available."""
    deps = {}

    # Core
    try:
        import torch
        deps["torch"] = True
        deps["torch_version"] = torch.__version__
        deps["mps_available"] = torch.backends.mps.is_available()
        deps["cuda_available"] = torch.cuda.is_available()
    except ImportError:
        deps["torch"] = False

    # Diffusers
    try:
        import diffusers
        deps["diffusers"] = True
        deps["diffusers_version"] = diffusers.__version__
    except ImportError:
        deps["diffusers"] = False

    # Image processing
    try:
        from PIL import Image
        deps["pillow"] = True
    except ImportError:
        deps["pillow"] = False

    try:
        import matplotlib
        deps["matplotlib"] = True
    except ImportError:
        deps["matplotlib"] = False

    try:
        import scipy
        deps["scipy"] = True
    except ImportError:
        deps["scipy"] = False

    try:
        import lpips
        deps["lpips"] = True
    except ImportError:
        deps["lpips"] = False

    # Sigma diffusion module
    try:
        import sigma_diffusion
        deps["sigma_diffusion"] = True
    except ImportError:
        deps["sigma_diffusion"] = False

    return deps


def print_dependency_report(deps: Dict[str, bool]) -> None:
    """Print a colored dependency report."""
    logger.info("=" * 60)
    logger.info("DEPENDENCY CHECK")
    logger.info("=" * 60)

    required = ["torch", "diffusers", "pillow"]
    optional = ["matplotlib", "scipy", "lpips", "sigma_diffusion"]

    for name in required:
        status = "OK" if deps.get(name, False) else "MISSING (required)"
        version = deps.get(f"{name}_version", "")
        extra = f" ({version})" if version else ""
        logger.info(f"  [{status}] {name}{extra}")

    if deps.get("torch"):
        mps = "available" if deps.get("mps_available") else "not available"
        cuda = "available" if deps.get("cuda_available") else "not available"
        logger.info(f"    MPS: {mps}")
        logger.info(f"    CUDA: {cuda}")

    for name in optional:
        status = "OK" if deps.get(name, False) else "missing (optional)"
        logger.info(f"  [{status}] {name}")

    logger.info("=" * 60)

    # Check required
    missing = [n for n in required if not deps.get(n, False)]
    if missing:
        logger.error(
            f"Missing required dependencies: {', '.join(missing)}\n"
            "Install with:\n"
            "  pip install torch diffusers transformers accelerate Pillow"
        )
        return False
    return True


def run_all(
    num_images: int = 10,
    num_steps: int = 28,
    device: str = "mps",
    dtype: str = "bfloat16",
    seed: int = 42,
    threshold_percentile: float = 25.0,
    min_keep: int = 8,
    skip_baseline: bool = False,
    evaluate_only: bool = False,
    no_lpips: bool = False,
    output_root: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Run the full Line 1 experiment pipeline.

    Returns a summary dict with all results.
    """
    if output_root is None:
        output_root = str(
            Path(__file__).parent.parent.parent / "results" / "line1"
        )

    baseline_dir = os.path.join(output_root, "baseline")
    guided_dir = os.path.join(output_root, "guided")
    eval_dir = os.path.join(output_root, "evaluation")

    os.makedirs(output_root, exist_ok=True)

    summary: Dict[str, Any] = {
        "experiment": "line1_image_diffusion",
        "config": {
            "num_images": num_images,
            "num_steps": num_steps,
            "device": device,
            "dtype": dtype,
            "seed": seed,
            "threshold_percentile": threshold_percentile,
            "min_keep": min_keep,
        },
        "timestamp_start": time.strftime("%Y-%m-%dT%H:%M:%S"),
    }

    # -----------------------------------------------------------------------
    # Dependency check
    # -----------------------------------------------------------------------
    deps = check_dependencies()
    deps_ok = print_dependency_report(deps)
    summary["dependencies"] = deps

    if not deps_ok and not evaluate_only:
        logger.error("Cannot proceed without required dependencies.")
        return summary

    # -----------------------------------------------------------------------
    # Stage 1: Baseline generation
    # -----------------------------------------------------------------------
    if not evaluate_only and not skip_baseline:
        logger.info("")
        logger.info("=" * 60)
        logger.info("STAGE 1: BASELINE GENERATION")
        logger.info("=" * 60)

        from run_baseline import run_baseline

        t0 = time.time()
        try:
            baseline_results = run_baseline(
                num_images=num_images,
                num_steps=num_steps,
                output_dir=baseline_dir,
                device=device,
                dtype=dtype,
                seed=seed,
            )
            summary["baseline"] = {
                "status": "success",
                "num_images": len(baseline_results),
                "avg_total_sigma": float(
                    sum(r.total_sigma for r in baseline_results)
                    / max(len(baseline_results), 1)
                ),
                "avg_generation_time_s": float(
                    sum(r.generation_time_s for r in baseline_results)
                    / max(len(baseline_results), 1)
                ),
                "wall_time_s": time.time() - t0,
            }
            logger.info(
                f"Baseline complete in {summary['baseline']['wall_time_s']:.1f}s"
            )
        except Exception as e:
            logger.error(f"Baseline generation failed: {e}")
            summary["baseline"] = {"status": "failed", "error": str(e)}
            # Try to continue if baseline files exist
            if not os.path.exists(os.path.join(baseline_dir, "sigma_profile.json")):
                return summary
    else:
        if skip_baseline:
            logger.info("Skipping baseline (--skip-baseline flag).")
        summary["baseline"] = {"status": "skipped"}

    # Verify baseline outputs exist
    baseline_sigma_path = os.path.join(baseline_dir, "sigma_profile.json")
    if not evaluate_only and not os.path.exists(baseline_sigma_path):
        logger.error(
            f"Baseline sigma profile not found at {baseline_sigma_path}.\n"
            "Cannot proceed with guided generation."
        )
        return summary

    # -----------------------------------------------------------------------
    # Stage 2: Sigma-guided generation
    # -----------------------------------------------------------------------
    if not evaluate_only:
        logger.info("")
        logger.info("=" * 60)
        logger.info("STAGE 2: SIGMA-GUIDED GENERATION")
        logger.info("=" * 60)

        from run_sigma_guided import run_sigma_guided

        t0 = time.time()
        try:
            speedup_metrics = run_sigma_guided(
                baseline_dir=baseline_dir,
                output_dir=guided_dir,
                threshold_percentile=threshold_percentile,
                min_keep=min_keep,
                device=device,
                dtype=dtype,
                seed=seed,
            )
            summary["guided"] = {
                "status": "success",
                "step_reduction_pct": speedup_metrics.get("step_reduction_pct"),
                "wall_clock_speedup": speedup_metrics.get("wall_clock_speedup"),
                "theoretical_speedup": speedup_metrics.get(
                    "speedup_ratio_theoretical"
                ),
                "wall_time_s": time.time() - t0,
            }
            logger.info(
                f"Guided generation complete in "
                f"{summary['guided']['wall_time_s']:.1f}s"
            )
        except Exception as e:
            logger.error(f"Guided generation failed: {e}")
            summary["guided"] = {"status": "failed", "error": str(e)}
            # Try to continue with evaluation if guided images exist
            if not os.path.isdir(os.path.join(guided_dir, "images")):
                return summary
    else:
        summary["guided"] = {"status": "skipped (evaluate-only)"}

    # -----------------------------------------------------------------------
    # Stage 3: Evaluation
    # -----------------------------------------------------------------------
    logger.info("")
    logger.info("=" * 60)
    logger.info("STAGE 3: EVALUATION")
    logger.info("=" * 60)

    # Check if we have the minimum needed for evaluation
    baseline_imgs = os.path.join(baseline_dir, "images")
    guided_imgs = os.path.join(guided_dir, "images")

    if not os.path.isdir(baseline_imgs) or not os.path.isdir(guided_imgs):
        logger.error(
            "Cannot evaluate: missing image directories.\n"
            f"  Baseline: {baseline_imgs} "
            f"({'exists' if os.path.isdir(baseline_imgs) else 'MISSING'})\n"
            f"  Guided: {guided_imgs} "
            f"({'exists' if os.path.isdir(guided_imgs) else 'MISSING'})"
        )
        summary["evaluation"] = {"status": "failed", "error": "Missing image dirs"}
        return summary

    from evaluate import evaluate

    t0 = time.time()
    try:
        eval_results = evaluate(
            baseline_dir=baseline_dir,
            guided_dir=guided_dir,
            output_dir=eval_dir,
            compute_lpips_flag=not no_lpips,
        )
        summary["evaluation"] = {
            "status": "success",
            "num_pairs": eval_results.get("num_pairs"),
            "psnr_mean": (
                eval_results.get("psnr", {}).get("mean")
            ),
            "ssim_mean": (
                eval_results.get("ssim", {}).get("mean")
            ),
            "wall_time_s": time.time() - t0,
        }
        if eval_results.get("lpips"):
            summary["evaluation"]["lpips_mean"] = eval_results["lpips"].get("mean")
        if eval_results.get("speedup"):
            summary["evaluation"]["speedup"] = eval_results["speedup"]
    except Exception as e:
        logger.error(f"Evaluation failed: {e}")
        summary["evaluation"] = {"status": "failed", "error": str(e)}

    # -----------------------------------------------------------------------
    # Final summary
    # -----------------------------------------------------------------------
    summary["timestamp_end"] = time.strftime("%Y-%m-%dT%H:%M:%S")
    total_time = sum(
        v.get("wall_time_s", 0)
        for v in [
            summary.get("baseline", {}),
            summary.get("guided", {}),
            summary.get("evaluation", {}),
        ]
        if isinstance(v, dict)
    )
    summary["total_wall_time_s"] = total_time

    # Save summary
    summary_path = os.path.join(output_root, "experiment_summary.json")
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2, default=str)
    logger.info(f"Saved experiment summary to {summary_path}")

    # Print final report
    logger.info("")
    logger.info("=" * 60)
    logger.info("LINE 1 EXPERIMENT COMPLETE")
    logger.info("=" * 60)
    logger.info(f"  Total wall time: {total_time:.1f}s")
    logger.info(f"  Output directory: {output_root}")
    logger.info("")
    logger.info("  Results:")
    for stage in ["baseline", "guided", "evaluation"]:
        s = summary.get(stage, {})
        status = s.get("status", "unknown") if isinstance(s, dict) else "unknown"
        logger.info(f"    {stage}: {status}")
    logger.info("")

    if summary.get("evaluation", {}).get("status") == "success":
        ev = summary["evaluation"]
        logger.info("  Quality metrics:")
        if ev.get("psnr_mean") is not None:
            logger.info(f"    PSNR: {ev['psnr_mean']:.2f} dB")
        if ev.get("ssim_mean") is not None:
            logger.info(f"    SSIM: {ev['ssim_mean']:.4f}")
        if ev.get("lpips_mean") is not None:
            logger.info(f"    LPIPS: {ev['lpips_mean']:.4f}")
        if ev.get("speedup"):
            sp = ev["speedup"]
            logger.info(
                f"    Speedup: {sp.get('wall_clock_speedup', 'N/A')}x "
                f"(step reduction: {sp.get('step_reduction_pct', 'N/A')}%)"
            )

    logger.info("=" * 60)

    return summary


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Line 1 Orchestrator: Full pipeline for Z-Image/FLUX diffusion "
            "with Sigma optimization."
        ),
    )
    parser.add_argument(
        "--num-images",
        type=int,
        default=10,
        help="Number of images to generate (default: 10 for quick test, 50 for full).",
    )
    parser.add_argument(
        "--num-steps",
        type=int,
        default=28,
        help="Number of denoising steps (default: 28).",
    )
    parser.add_argument(
        "--device",
        type=str,
        default="mps",
        choices=["mps", "cuda", "cpu"],
        help="PyTorch device (default: mps for Mac M1).",
    )
    parser.add_argument(
        "--dtype",
        type=str,
        default="bfloat16",
        choices=["bfloat16", "float16", "float32"],
        help="Float precision (default: bfloat16).",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed (default: 42).",
    )
    parser.add_argument(
        "--threshold-percentile",
        type=float,
        default=25.0,
        help="Sigma percentile threshold for step skipping (default: 25).",
    )
    parser.add_argument(
        "--min-keep",
        type=int,
        default=8,
        help="Minimum number of denoising steps to keep (default: 8).",
    )
    parser.add_argument(
        "--skip-baseline",
        action="store_true",
        help="Skip baseline generation (use existing results).",
    )
    parser.add_argument(
        "--evaluate-only",
        action="store_true",
        help="Only run evaluation on existing results.",
    )
    parser.add_argument(
        "--no-lpips",
        action="store_true",
        help="Skip LPIPS computation.",
    )
    parser.add_argument(
        "--output-root",
        type=str,
        default=None,
        help="Root output directory (default: sigma_diffusion/results/line1/).",
    )
    args = parser.parse_args()

    run_all(
        num_images=args.num_images,
        num_steps=args.num_steps,
        device=args.device,
        dtype=args.dtype,
        seed=args.seed,
        threshold_percentile=args.threshold_percentile,
        min_keep=args.min_keep,
        skip_baseline=args.skip_baseline,
        evaluate_only=args.evaluate_only,
        no_lpips=args.no_lpips,
        output_root=args.output_root,
    )


if __name__ == "__main__":
    main()
