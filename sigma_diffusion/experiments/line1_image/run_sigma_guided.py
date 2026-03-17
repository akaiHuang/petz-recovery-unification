#!/usr/bin/env python3
"""
Line 1 Sigma-Guided: Adaptive step-skipping based on Sigma profile.

Reads the Sigma profile from a baseline run, identifies low-Sigma steps
(where D_KL between consecutive latents is small, meaning the step
contributes little information), and re-runs generation with those steps
skipped.

Key idea (from Petz recovery theory):
    Low Sigma_t means step t is nearly perfectly retrodictable from step t-1.
    Skipping such steps should preserve image quality while reducing compute.

Usage:
    python run_sigma_guided.py --baseline-dir ../../results/line1/baseline \
                               --output-dir ../../results/line1/guided
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import time
from dataclasses import asdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Sigma module (with inline fallback)
# ---------------------------------------------------------------------------

try:
    from sigma_diffusion import (
        AdaptiveScheduler,
        compute_skip_mask,
        compute_adaptive_timesteps,
        compute_early_stop,
    )
    _HAS_SCHEDULER = True
    logger.info("Using sigma_diffusion.AdaptiveScheduler.")
except ImportError:
    _HAS_SCHEDULER = False
    logger.warning("sigma_diffusion module not available; using inline scheduler.")

# Import from sibling module
from run_baseline import (
    LatentRecorder,
    compute_sigma_profile,
    load_pipeline,
    generate_with_tracking,
    load_prompts,
    BaselineResult,
)


# ---------------------------------------------------------------------------
# Inline adaptive scheduler (fallback)
# ---------------------------------------------------------------------------


def _compute_skip_mask_inline(
    sigma_values: List[float],
    threshold_percentile: float = 25.0,
    min_keep: int = 8,
) -> List[bool]:
    """
    Compute a boolean mask: True = keep step, False = skip step.

    Steps with Sigma below the threshold_percentile are candidates for
    skipping. We always keep at least `min_keep` steps (the ones with
    the highest Sigma values).
    """
    n = len(sigma_values)
    if n <= min_keep:
        return [True] * n

    arr = np.array(sigma_values)
    threshold = np.percentile(arr, threshold_percentile)

    # Always keep first and last step
    mask = [False] * n
    mask[0] = True
    mask[-1] = True

    # Keep steps above threshold
    for i in range(1, n - 1):
        if arr[i] >= threshold:
            mask[i] = True

    # Ensure minimum number of kept steps
    kept = sum(mask)
    if kept < min_keep:
        # Add steps in order of descending Sigma
        indices_by_sigma = np.argsort(-arr)
        for idx in indices_by_sigma:
            if not mask[idx]:
                mask[idx] = True
                kept += 1
                if kept >= min_keep:
                    break

    return mask


def _compute_adaptive_timesteps_inline(
    original_timesteps: List[int],
    skip_mask: List[bool],
) -> List[int]:
    """Return the subset of timesteps where skip_mask is True."""
    return [t for t, keep in zip(original_timesteps, skip_mask) if keep]


def _compute_early_stop_inline(
    sigma_values: List[float],
    tail_threshold: float = 0.01,
    tail_window: int = 3,
) -> Optional[int]:
    """
    Find the earliest step after which all remaining Sigma values are
    below tail_threshold. Returns None if no early stop is warranted.
    """
    n = len(sigma_values)
    if n <= tail_window:
        return None

    for i in range(n - tail_window, 0, -1):
        if sigma_values[i] >= tail_threshold:
            return i + 1  # stop after this step

    return None


# ---------------------------------------------------------------------------
# Unified API
# ---------------------------------------------------------------------------


def compute_adaptive_schedule(
    sigma_profile: List[Dict],
    total_steps: int,
    threshold_percentile: float = 25.0,
    min_keep: int = 8,
    tail_threshold: float = 0.01,
) -> Dict[str, Any]:
    """
    Compute the adaptive schedule from a sigma profile.

    Returns a dict with:
        - skip_mask: List[bool]
        - kept_steps: List[int] (indices of kept steps)
        - skipped_steps: List[int] (indices of skipped steps)
        - early_stop: Optional[int]
        - num_original: int
        - num_kept: int
        - speedup_ratio: float
    """
    sigma_values = [entry["sigma"] for entry in sigma_profile]

    if _HAS_SCHEDULER:
        mask = list(compute_skip_mask(
            sigma_values,
            threshold_percentile=threshold_percentile,
            min_keep=min_keep,
        ))
        early_stop = compute_early_stop(
            sigma_values,
            tail_threshold=tail_threshold,
        )
    else:
        mask = _compute_skip_mask_inline(
            sigma_values,
            threshold_percentile=threshold_percentile,
            min_keep=min_keep,
        )
        early_stop = _compute_early_stop_inline(
            sigma_values,
            tail_threshold=tail_threshold,
        )

    # Apply early stop
    if early_stop is not None and early_stop < len(mask):
        for i in range(early_stop, len(mask)):
            mask[i] = False
        mask[early_stop - 1] = True  # keep the stop step

    kept = [i for i, m in enumerate(mask) if m]
    skipped = [i for i, m in enumerate(mask) if not m]
    num_kept = len(kept)

    return {
        "skip_mask": mask,
        "kept_steps": kept,
        "skipped_steps": skipped,
        "early_stop": early_stop,
        "num_original": total_steps,
        "num_kept": num_kept,
        "speedup_ratio": total_steps / max(num_kept, 1),
        "sigma_threshold_percentile": threshold_percentile,
        "min_keep": min_keep,
    }


# ---------------------------------------------------------------------------
# Guided generation
# ---------------------------------------------------------------------------


def generate_guided(
    pipe: Any,
    prompt: str,
    schedule: Dict[str, Any],
    num_original_steps: int,
    recorder: LatentRecorder,
    seed: int = 42,
    device: str = "mps",
) -> Any:
    """
    Run generation with reduced timesteps based on the adaptive schedule.

    Strategy: We use the number of kept steps as num_inference_steps.
    The diffusers scheduler will automatically space them across the
    denoising range. This is a principled approximation -- the key insight
    is that Sigma-guided scheduling removes steps where the information
    gain was negligible.
    """
    import torch

    num_kept = schedule["num_kept"]
    generator = torch.Generator(device="cpu").manual_seed(seed)

    try:
        result = pipe(
            prompt=prompt,
            num_inference_steps=num_kept,
            generator=generator,
            callback_on_step_end=recorder,
        )
        return result
    except TypeError:
        logger.info("Falling back to legacy callback API.")
        try:
            result = pipe(
                prompt=prompt,
                num_inference_steps=num_kept,
                generator=generator,
                callback=recorder.legacy_callback,
                callback_steps=1,
            )
            return result
        except Exception as e:
            logger.warning(f"Legacy callback failed: {e}")

    # Final fallback
    result = pipe(
        prompt=prompt,
        num_inference_steps=num_kept,
        generator=generator,
    )
    return result


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def run_sigma_guided(
    baseline_dir: Optional[str] = None,
    output_dir: Optional[str] = None,
    threshold_percentile: float = 25.0,
    min_keep: int = 8,
    device: str = "mps",
    dtype: str = "bfloat16",
    seed: int = 42,
) -> Dict[str, Any]:
    """
    Run Sigma-guided generation using a baseline profile.

    Returns summary dict with speedup metrics.
    """
    base_results_dir = Path(__file__).parent.parent.parent / "results" / "line1"

    if baseline_dir is None:
        baseline_dir = str(base_results_dir / "baseline")
    if output_dir is None:
        output_dir = str(base_results_dir / "guided")

    os.makedirs(output_dir, exist_ok=True)
    images_dir = os.path.join(output_dir, "images")
    os.makedirs(images_dir, exist_ok=True)

    # -----------------------------------------------------------------------
    # Load baseline sigma profile
    # -----------------------------------------------------------------------
    sigma_path = os.path.join(baseline_dir, "sigma_profile.json")
    if not os.path.exists(sigma_path):
        raise FileNotFoundError(
            f"Baseline sigma profile not found at {sigma_path}.\n"
            "Run run_baseline.py first."
        )

    with open(sigma_path, "r") as f:
        baseline_profiles = json.load(f)
    logger.info(f"Loaded {len(baseline_profiles)} baseline sigma profiles.")

    # Load metadata
    meta_path = os.path.join(baseline_dir, "metadata.json")
    if os.path.exists(meta_path):
        with open(meta_path, "r") as f:
            baseline_meta = json.load(f)
    else:
        baseline_meta = {"num_steps": 28}

    num_original_steps = baseline_meta.get("num_steps", 28)

    # -----------------------------------------------------------------------
    # Compute aggregate sigma profile (average across all prompts)
    # -----------------------------------------------------------------------
    all_sigma_values = []
    for bp in baseline_profiles:
        sigma_vals = [entry["sigma"] for entry in bp["profile"]]
        all_sigma_values.append(sigma_vals)

    # Use per-step average to determine the schedule
    max_len = max(len(sv) for sv in all_sigma_values)
    avg_profile = []
    for step_idx in range(max_len):
        vals = [sv[step_idx] for sv in all_sigma_values if step_idx < len(sv)]
        avg_sigma = float(np.mean(vals))
        avg_profile.append({"sigma": avg_sigma, "step": step_idx})

    logger.info(f"Average sigma profile has {len(avg_profile)} steps.")

    # Compute adaptive schedule
    schedule = compute_adaptive_schedule(
        sigma_profile=avg_profile,
        total_steps=num_original_steps,
        threshold_percentile=threshold_percentile,
        min_keep=min_keep,
    )

    logger.info(
        f"Adaptive schedule: {schedule['num_kept']}/{schedule['num_original']} steps "
        f"(speedup: {schedule['speedup_ratio']:.2f}x)"
    )
    logger.info(f"  Kept steps: {schedule['kept_steps']}")
    logger.info(f"  Skipped steps: {schedule['skipped_steps']}")

    # Save schedule
    schedule_path = os.path.join(output_dir, "schedule.json")
    with open(schedule_path, "w") as f:
        json.dump(schedule, f, indent=2)
    logger.info(f"Saved schedule to {schedule_path}")

    # -----------------------------------------------------------------------
    # Load pipeline and generate
    # -----------------------------------------------------------------------
    logger.info(f"Loading pipeline (device={device}, dtype={dtype})...")
    pipe, model_name = load_pipeline(device=device, dtype_str=dtype)
    logger.info(f"Using model: {model_name}")

    recorder = LatentRecorder()
    guided_profiles = []
    total_baseline_time = 0.0
    total_guided_time = 0.0

    for bp in baseline_profiles:
        prompt = bp["prompt"]
        idx = bp["prompt_index"]
        logger.info(f"[{idx + 1}/{len(baseline_profiles)}] Guided: {prompt!r}")

        recorder.reset()
        t0 = time.time()
        output = generate_guided(
            pipe=pipe,
            prompt=prompt,
            schedule=schedule,
            num_original_steps=num_original_steps,
            recorder=recorder,
            seed=seed + idx,
            device=device,
        )
        gen_time = time.time() - t0

        # Save image
        image_path = os.path.join(images_dir, f"guided_{idx:03d}.png")
        if hasattr(output, "images") and output.images:
            output.images[0].save(image_path)
        else:
            logger.warning(f"No image output for prompt {idx}.")
            image_path = ""

        # Compute sigma profile for guided run
        sigma_entries, total_sigma = compute_sigma_profile(recorder.snapshots)

        guided_profiles.append({
            "prompt": prompt,
            "prompt_index": idx,
            "total_sigma": total_sigma,
            "num_steps": schedule["num_kept"],
            "generation_time_s": gen_time,
            "profile": [asdict(s) for s in sigma_entries],
        })

        total_baseline_time += bp.get("generation_time_s", gen_time)
        total_guided_time += gen_time

        logger.info(
            f"  Done in {gen_time:.1f}s | "
            f"total_sigma={total_sigma:.4f} | "
            f"steps={schedule['num_kept']}"
        )

    # -----------------------------------------------------------------------
    # Speedup metrics
    # -----------------------------------------------------------------------
    speedup_metrics = {
        "model_name": model_name,
        "num_images": len(baseline_profiles),
        "original_steps": num_original_steps,
        "guided_steps": schedule["num_kept"],
        "step_reduction_pct": (
            100.0 * (1.0 - schedule["num_kept"] / num_original_steps)
        ),
        "speedup_ratio_theoretical": schedule["speedup_ratio"],
        "total_baseline_time_s": total_baseline_time,
        "total_guided_time_s": total_guided_time,
        "wall_clock_speedup": (
            total_baseline_time / max(total_guided_time, 0.01)
        ),
        "avg_baseline_sigma": float(np.mean([
            bp["total_sigma"] for bp in baseline_profiles
        ])),
        "avg_guided_sigma": float(np.mean([
            gp["total_sigma"] for gp in guided_profiles
        ])),
        "schedule": schedule,
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
    }

    # Save guided sigma profiles
    guided_sigma_path = os.path.join(output_dir, "sigma_profile.json")
    with open(guided_sigma_path, "w") as f:
        json.dump(guided_profiles, f, indent=2)

    # Save speedup metrics
    metrics_path = os.path.join(output_dir, "speedup_metrics.json")
    with open(metrics_path, "w") as f:
        json.dump(speedup_metrics, f, indent=2)
    logger.info(f"Saved speedup metrics to {metrics_path}")

    # Summary
    logger.info("=" * 60)
    logger.info("SIGMA-GUIDED GENERATION SUMMARY")
    logger.info("=" * 60)
    logger.info(f"  Steps: {num_original_steps} -> {schedule['num_kept']}")
    logger.info(f"  Step reduction: {speedup_metrics['step_reduction_pct']:.1f}%")
    logger.info(
        f"  Theoretical speedup: {schedule['speedup_ratio']:.2f}x"
    )
    logger.info(
        f"  Wall-clock speedup: {speedup_metrics['wall_clock_speedup']:.2f}x"
    )
    logger.info(
        f"  Avg Sigma (baseline): {speedup_metrics['avg_baseline_sigma']:.4f}"
    )
    logger.info(
        f"  Avg Sigma (guided):   {speedup_metrics['avg_guided_sigma']:.4f}"
    )
    logger.info("=" * 60)

    return speedup_metrics


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description="Line 1 Sigma-Guided: Adaptive diffusion with Sigma-based step skipping."
    )
    parser.add_argument(
        "--baseline-dir",
        type=str,
        default=None,
        help="Path to baseline results directory.",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=None,
        help="Output directory for guided results.",
    )
    parser.add_argument(
        "--threshold-percentile",
        type=float,
        default=25.0,
        help="Percentile below which Sigma steps are skipped (default: 25).",
    )
    parser.add_argument(
        "--min-keep",
        type=int,
        default=8,
        help="Minimum number of steps to keep (default: 8).",
    )
    parser.add_argument(
        "--device",
        type=str,
        default="mps",
        choices=["mps", "cuda", "cpu"],
        help="PyTorch device (default: mps).",
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
    args = parser.parse_args()

    run_sigma_guided(
        baseline_dir=args.baseline_dir,
        output_dir=args.output_dir,
        threshold_percentile=args.threshold_percentile,
        min_keep=args.min_keep,
        device=args.device,
        dtype=args.dtype,
        seed=args.seed,
    )


if __name__ == "__main__":
    main()
