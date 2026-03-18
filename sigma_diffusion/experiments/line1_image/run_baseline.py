#!/usr/bin/env python3
"""
Line 1 Baseline: Z-Image / FLUX diffusion with per-step Sigma tracking.

Runs the full diffusion process, recording latent statistics at every step.
Computes Sigma_t = D_KL(p_{t-1} || p_t) using Gaussian approximation of
the latent distribution at each denoising step.

Theoretical basis:
    Sigma_t = D_KL(p_{t-1} || p_t)
    Petz fidelity bound:  F >= exp(-Sigma / 2)

See: Huang (2026), "Petz Recovery Map as a Universal Retrodiction Functor."

Usage:
    python run_baseline.py --num-images 10 --output-dir ../../results/line1/baseline
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import sys
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Sigma computation (inline fallback if sigma_diffusion module not ready)
# ---------------------------------------------------------------------------

try:
    from sigma_diffusion import kl_divergence_gaussian, SigmaTracker
    logger.info("Using sigma_diffusion module for Sigma computation.")
    _HAS_SIGMA_MODULE = True
except ImportError:
    logger.warning(
        "sigma_diffusion module not available; using inline Gaussian KL."
    )
    _HAS_SIGMA_MODULE = False


def _kl_divergence_gaussian_inline(
    mu_p: np.ndarray,
    var_p: np.ndarray,
    mu_q: np.ndarray,
    var_q: np.ndarray,
) -> float:
    """
    KL(p || q) for diagonal Gaussians.

    D_KL = 0.5 * sum[ log(var_q/var_p) + (var_p + (mu_p-mu_q)^2)/var_q - 1 ]
    """
    eps = 1e-8
    var_p = np.maximum(var_p, eps)
    var_q = np.maximum(var_q, eps)
    kl = 0.5 * np.sum(
        np.log(var_q / var_p)
        + (var_p + (mu_p - mu_q) ** 2) / var_q
        - 1.0
    )
    return float(kl)


def compute_kl(
    mu_p: np.ndarray,
    var_p: np.ndarray,
    mu_q: np.ndarray,
    var_q: np.ndarray,
) -> float:
    """Compute KL divergence, using the module if available."""
    if _HAS_SIGMA_MODULE:
        return float(kl_divergence_gaussian(mu_p, var_p, mu_q, var_q))
    return _kl_divergence_gaussian_inline(mu_p, var_p, mu_q, var_q)


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------


@dataclass
class LatentSnapshot:
    """Statistics of the latent tensor at a single denoising step."""
    step: int
    timestep: float
    mean: float
    std: float
    var: float
    min_val: float
    max_val: float
    shape: List[int]
    # Per-channel stats (for richer analysis)
    channel_means: List[float] = field(default_factory=list)
    channel_vars: List[float] = field(default_factory=list)


@dataclass
class SigmaEntry:
    """Sigma value between consecutive steps."""
    step: int
    sigma: float
    petz_fidelity_bound: float  # exp(-sigma/2)


@dataclass
class BaselineResult:
    """Full result for one prompt."""
    prompt: str
    prompt_index: int
    num_steps: int
    latent_stats: List[LatentSnapshot]
    sigma_profile: List[SigmaEntry]
    total_sigma: float
    generation_time_s: float
    image_path: str


# ---------------------------------------------------------------------------
# Latent callback
# ---------------------------------------------------------------------------


class LatentRecorder:
    """
    Callback for diffusers pipelines to record latent statistics at each step.

    Compatible with both `callback_on_step_end` (newer API) and the legacy
    `callback` parameter.
    """

    def __init__(self) -> None:
        self.snapshots: List[LatentSnapshot] = []
        self._step_counter: int = 0

    def reset(self) -> None:
        self.snapshots = []
        self._step_counter = 0

    def __call__(
        self,
        pipe: Any,
        step_index: int,
        timestep: Any,
        callback_kwargs: Dict[str, Any],
    ) -> Dict[str, Any]:
        """callback_on_step_end signature for diffusers >= 0.25."""
        # Extract latents from callback kwargs
        latents = callback_kwargs.get("latents", None)
        if latents is None:
            # Some pipelines use "denoised" or the pipe's internal state
            latents = getattr(pipe, "_current_latents", None)

        if latents is not None:
            self._record(latents, step_index, float(timestep))

        return callback_kwargs

    def legacy_callback(
        self,
        step: int,
        timestep: int,
        latents: Any,
    ) -> None:
        """Legacy callback(step, timestep, latents) for older diffusers."""
        self._record(latents, step, float(timestep))

    def _record(self, latents: Any, step: int, timestep: float) -> None:
        """Record statistics from a latent tensor."""
        import torch

        with torch.no_grad():
            # Move to CPU to avoid MPS issues with some stat operations
            lat_cpu = latents.detach().cpu().float()

            snapshot = LatentSnapshot(
                step=step,
                timestep=timestep,
                mean=float(lat_cpu.mean()),
                std=float(lat_cpu.std()),
                var=float(lat_cpu.var()),
                min_val=float(lat_cpu.min()),
                max_val=float(lat_cpu.max()),
                shape=list(lat_cpu.shape),
                channel_means=[],
                channel_vars=[],
            )

            # Per-channel stats (dim=0 is batch, dim=1 is channel)
            if lat_cpu.ndim >= 4:
                # Shape: [B, C, H, W]
                for c in range(lat_cpu.shape[1]):
                    ch = lat_cpu[:, c, :, :]
                    snapshot.channel_means.append(float(ch.mean()))
                    snapshot.channel_vars.append(float(ch.var()))
            elif lat_cpu.ndim == 3:
                # Shape: [B, seq_len, hidden_dim] (FLUX-like)
                snapshot.channel_means = [float(lat_cpu.mean())]
                snapshot.channel_vars = [float(lat_cpu.var())]

            self.snapshots.append(snapshot)
            self._step_counter += 1


# ---------------------------------------------------------------------------
# Sigma profile computation
# ---------------------------------------------------------------------------


def compute_sigma_profile(
    snapshots: List[LatentSnapshot],
) -> Tuple[List[SigmaEntry], float]:
    """
    Compute Sigma_t = D_KL(p_{t-1} || p_t) between consecutive latent
    distributions, using the Gaussian approximation (per-channel).

    Returns (list of SigmaEntry, total_sigma).
    """
    entries: List[SigmaEntry] = []
    total_sigma = 0.0

    for i in range(1, len(snapshots)):
        prev = snapshots[i - 1]
        curr = snapshots[i]

        if prev.channel_means and curr.channel_means:
            mu_p = np.array(prev.channel_means)
            var_p = np.array(prev.channel_vars)
            mu_q = np.array(curr.channel_means)
            var_q = np.array(curr.channel_vars)
        else:
            mu_p = np.array([prev.mean])
            var_p = np.array([prev.var])
            mu_q = np.array([curr.mean])
            var_q = np.array([curr.var])

        sigma = compute_kl(mu_p, var_p, mu_q, var_q)
        petz_bound = float(np.exp(-sigma / 2.0))

        entries.append(SigmaEntry(
            step=curr.step,
            sigma=sigma,
            petz_fidelity_bound=petz_bound,
        ))
        total_sigma += sigma

    return entries, total_sigma


# ---------------------------------------------------------------------------
# Pipeline loader
# ---------------------------------------------------------------------------


def load_pipeline(device: str = "mps", dtype_str: str = "bfloat16"):
    """
    Load a diffusion pipeline. Tries Z-Image first, then FLUX.1-schnell.

    Returns (pipeline, model_name).
    """
    import torch

    dtype_map = {
        "bfloat16": torch.bfloat16,
        "float16": torch.float16,
        "float32": torch.float32,
    }
    dtype = dtype_map.get(dtype_str, torch.bfloat16)

    # --- Attempt 1: Z-Image ---
    try:
        from diffusers import ZImagePipeline
        logger.info("Loading Z-Image pipeline...")
        pipe = ZImagePipeline.from_pretrained(
            "Tongyi-MAI/Z-Image",
            torch_dtype=dtype,
        )
        pipe = pipe.to(device)
        return pipe, "z-image-base"
    except (ImportError, Exception) as e:
        logger.warning(f"Z-Image not available: {e}")

    # --- Attempt 2: FLUX.1-schnell ---
    try:
        from diffusers import FluxPipeline
        logger.info("Loading FLUX.1-schnell pipeline...")
        pipe = FluxPipeline.from_pretrained(
            "black-forest-labs/FLUX.1-schnell",
            torch_dtype=dtype,
        )
        pipe = pipe.to(device)
        return pipe, "FLUX.1-schnell"
    except (ImportError, Exception) as e:
        logger.warning(f"FLUX.1-schnell not available: {e}")

    # --- Attempt 3: Stable Diffusion XL (widely available) ---
    try:
        from diffusers import StableDiffusionXLPipeline
        logger.info("Loading Stable Diffusion XL (fallback)...")
        pipe = StableDiffusionXLPipeline.from_pretrained(
            "stabilityai/stable-diffusion-xl-base-1.0",
            torch_dtype=dtype,
            variant="fp16" if dtype == torch.float16 else None,
        )
        pipe = pipe.to(device)
        return pipe, "SDXL-base-1.0"
    except (ImportError, Exception) as e:
        logger.warning(f"SDXL not available: {e}")

    # --- Attempt 4: Stable Diffusion 2.1 ---
    try:
        from diffusers import StableDiffusionPipeline
        logger.info("Loading Stable Diffusion 2.1 (fallback)...")
        pipe = StableDiffusionPipeline.from_pretrained(
            "stabilityai/stable-diffusion-2-1",
            torch_dtype=dtype,
        )
        pipe = pipe.to(device)
        return pipe, "SD-2.1"
    except (ImportError, Exception) as e:
        logger.warning(f"SD 2.1 not available: {e}")

    raise RuntimeError(
        "No diffusion pipeline available. Install diffusers and download a model:\n"
        "  pip install diffusers transformers accelerate\n"
        "Then try: huggingface-cli download stabilityai/stable-diffusion-xl-base-1.0"
    )


# ---------------------------------------------------------------------------
# MPS-safe generation
# ---------------------------------------------------------------------------


def generate_with_tracking(
    pipe: Any,
    prompt: str,
    num_steps: int,
    recorder: LatentRecorder,
    seed: int = 42,
    device: str = "mps",
) -> Any:
    """
    Run generation with latent recording.
    Handles MPS-specific issues with CPU fallback for problematic operations.
    """
    import torch

    generator = torch.Generator(device="cpu").manual_seed(seed)

    try:
        # Try with callback_on_step_end (diffusers >= 0.25)
        result = pipe(
            prompt=prompt,
            num_inference_steps=num_steps,
            generator=generator,
            callback_on_step_end=recorder,
        )
        return result
    except TypeError:
        # Older diffusers: use legacy callback
        logger.info("Falling back to legacy callback API.")
        try:
            result = pipe(
                prompt=prompt,
                num_inference_steps=num_steps,
                generator=generator,
                callback=recorder.legacy_callback,
                callback_steps=1,
            )
            return result
        except Exception as e:
            logger.warning(f"Legacy callback failed: {e}")

    # Final fallback: no callback, just generate
    logger.warning("Running without step-level tracking (no callback support).")
    result = pipe(
        prompt=prompt,
        num_inference_steps=num_steps,
        generator=generator,
    )
    return result


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def load_prompts(path: Optional[str] = None, num_images: int = 10) -> List[str]:
    """Load prompts from file, return first num_images."""
    if path is None:
        path = str(Path(__file__).parent / "prompts.txt")

    with open(path, "r") as f:
        prompts = [line.strip() for line in f if line.strip()]

    return prompts[:num_images]


def run_baseline(
    num_images: int = 10,
    num_steps: int = 28,
    output_dir: Optional[str] = None,
    device: str = "mps",
    dtype: str = "bfloat16",
    seed: int = 42,
) -> List[BaselineResult]:
    """
    Run baseline diffusion generation with Sigma tracking.

    Args:
        num_images: Number of prompts to process.
        num_steps: Number of denoising steps.
        output_dir: Where to save results.
        device: PyTorch device ("mps", "cuda", "cpu").
        dtype: Float precision ("bfloat16", "float16", "float32").
        seed: Random seed for reproducibility.

    Returns:
        List of BaselineResult for each generated image.
    """
    if output_dir is None:
        output_dir = str(
            Path(__file__).parent.parent.parent / "results" / "line1" / "baseline"
        )

    os.makedirs(output_dir, exist_ok=True)
    images_dir = os.path.join(output_dir, "images")
    os.makedirs(images_dir, exist_ok=True)

    # Load prompts
    prompts = load_prompts(num_images=num_images)
    logger.info(f"Loaded {len(prompts)} prompts.")

    # Load pipeline
    logger.info(f"Loading pipeline (device={device}, dtype={dtype})...")
    pipe, model_name = load_pipeline(device=device, dtype_str=dtype)
    logger.info(f"Using model: {model_name}")

    # Run generation for each prompt
    recorder = LatentRecorder()
    results: List[BaselineResult] = []

    for idx, prompt in enumerate(prompts):
        logger.info(f"[{idx + 1}/{len(prompts)}] Generating: {prompt!r}")
        recorder.reset()

        t0 = time.time()
        output = generate_with_tracking(
            pipe=pipe,
            prompt=prompt,
            num_steps=num_steps,
            recorder=recorder,
            seed=seed + idx,
            device=device,
        )
        gen_time = time.time() - t0

        # Save image
        image_path = os.path.join(images_dir, f"baseline_{idx:03d}.png")
        if hasattr(output, "images") and output.images:
            output.images[0].save(image_path)
        else:
            logger.warning(f"No image output for prompt {idx}.")
            image_path = ""

        # Compute Sigma profile
        sigma_profile, total_sigma = compute_sigma_profile(recorder.snapshots)

        result = BaselineResult(
            prompt=prompt,
            prompt_index=idx,
            num_steps=num_steps,
            latent_stats=recorder.snapshots,
            sigma_profile=sigma_profile,
            total_sigma=total_sigma,
            generation_time_s=gen_time,
            image_path=image_path,
        )
        results.append(result)

        logger.info(
            f"  Done in {gen_time:.1f}s | "
            f"total_sigma={total_sigma:.4f} | "
            f"steps_recorded={len(recorder.snapshots)}"
        )

    # -----------------------------------------------------------------------
    # Save outputs
    # -----------------------------------------------------------------------

    # 1. sigma_profile.json
    sigma_data = []
    for r in results:
        sigma_data.append({
            "prompt": r.prompt,
            "prompt_index": r.prompt_index,
            "total_sigma": r.total_sigma,
            "num_steps": r.num_steps,
            "generation_time_s": r.generation_time_s,
            "profile": [asdict(s) for s in r.sigma_profile],
        })

    sigma_path = os.path.join(output_dir, "sigma_profile.json")
    with open(sigma_path, "w") as f:
        json.dump(sigma_data, f, indent=2)
    logger.info(f"Saved sigma profile to {sigma_path}")

    # 2. latent_stats.json
    latent_data = []
    for r in results:
        latent_data.append({
            "prompt": r.prompt,
            "prompt_index": r.prompt_index,
            "stats": [asdict(s) for s in r.latent_stats],
        })

    latent_path = os.path.join(output_dir, "latent_stats.json")
    with open(latent_path, "w") as f:
        json.dump(latent_data, f, indent=2)
    logger.info(f"Saved latent stats to {latent_path}")

    # 3. Run metadata
    meta = {
        "model_name": model_name,
        "device": device,
        "dtype": dtype,
        "num_steps": num_steps,
        "num_images": len(results),
        "seed": seed,
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
    }
    meta_path = os.path.join(output_dir, "metadata.json")
    with open(meta_path, "w") as f:
        json.dump(meta, f, indent=2)
    logger.info(f"Saved metadata to {meta_path}")

    logger.info(
        f"Baseline complete: {len(results)} images, "
        f"avg_sigma={np.mean([r.total_sigma for r in results]):.4f}"
    )

    return results


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description="Line 1 Baseline: Diffusion generation with per-step Sigma tracking."
    )
    parser.add_argument(
        "--num-images",
        type=int,
        default=10,
        help="Number of images to generate (default: 10).",
    )
    parser.add_argument(
        "--num-steps",
        type=int,
        default=28,
        help="Number of denoising steps (default: 28).",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=None,
        help="Output directory (default: sigma_diffusion/results/line1/baseline).",
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

    run_baseline(
        num_images=args.num_images,
        num_steps=args.num_steps,
        output_dir=args.output_dir,
        device=args.device,
        dtype=args.dtype,
        seed=args.seed,
    )


if __name__ == "__main__":
    main()
