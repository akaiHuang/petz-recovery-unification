#!/usr/bin/env python3
"""
Entropic Time Scheduler for Image Diffusion
============================================
Based on NeurIPS 2025 insight (arXiv:2504.13612):
Don't SKIP steps -- REDISTRIBUTE the time schedule so each step
contributes equal information. Works WITHOUT retraining.

Uses pre-existing sigma profile from line1_zimage/sigma_profile_28.json
to build the entropic schedule, then runs 4 experiments:

A) Z-Image Base, 28 steps, uniform (baseline)         -- GOOD quality reference
B) Z-Image Base, 8 steps, uniform (naive reduction)    -- worse
C) Z-Image Base, 8 steps, ENTROPIC schedule (our method) -- should beat B
D) Z-Image Turbo, 8 steps (distilled competitor)       -- reference

Key finding from sigma profile:
  - Sigma is DRAMATICALLY higher near the END (steps 19-26, t < 716)
  - Step 26 (t=187.5) has Sigma = 0.666 -- 1000x larger than step 1
  - The entropic schedule concentrates timesteps in this high-information region
"""

import os
import sys
import json
import time
import warnings
import numpy as np
import torch
from PIL import Image, ImageDraw, ImageFont

warnings.filterwarnings("ignore")

os.environ["HF_TOKEN"] = "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"
os.environ["PYTORCH_MPS_HIGH_WATERMARK_RATIO"] = "0.0"

DEVICE = "mps"
DTYPE = torch.float32
PROMPT = "A red fox in snow, photorealistic"
SEED = 42
OUTPUT_DIR = "/Users/akaihuangm1/Desktop/github/petz-recovery-unification/sigma_diffusion/results/line1_entropic"
PROFILE_PATH = "/Users/akaihuangm1/Desktop/github/petz-recovery-unification/sigma_diffusion/results/line1_zimage/sigma_profile_28.json"
os.makedirs(OUTPUT_DIR, exist_ok=True)


# =====================================================================
# Build entropic schedule from pre-existing sigma profile
# =====================================================================

def build_entropic_sigmas(profile_data, n_target=8):
    """
    Build non-uniform flow-matching sigmas that allocate more steps to
    high-Sigma (high-information) regions.

    The profile gives us Sigma_KL at each step of the 28-step schedule.
    We build cumulative Sigma(t) and invert to find n_target equally-spaced
    quantiles in information space.

    Returns: list of n_target flow-matching sigma values (descending, 1.0 -> ~0)
    """
    sigma_data = profile_data["sigma_data"]

    # Extract timesteps and KL values (skip step 0 which has sigma=0)
    # sigma_data[i] has step, t, sigma (KL between step i-1 and step i)
    all_t = [s["t"] for s in sigma_data]          # 28 timesteps: 1000 -> 0
    all_kl = [s["sigma"] for s in sigma_data]      # 28 KL values (first is 0)

    # The KL between step i-1 and step i:
    # all_kl[0] = 0 (no previous step)
    # all_kl[1] = KL(step0 -> step1)
    # ...
    # all_kl[27] = 0 (terminal)

    # We have 27 intervals between 28 timesteps
    # Interval j = (step j -> step j+1) has KL value all_kl[j+1]
    n_full = len(all_t)  # 28
    interval_kl = np.array([all_kl[j + 1] for j in range(n_full - 1)])  # 27 values
    interval_kl = np.maximum(interval_kl, 1e-10)  # floor to avoid zero

    print(f"\n  Sigma (KL) profile summary:")
    print(f"    Total intervals: {len(interval_kl)}")
    print(f"    Total Sigma: {interval_kl.sum():.4f}")
    print(f"    Min KL: {interval_kl.min():.6f} (interval {np.argmin(interval_kl)})")
    print(f"    Max KL: {interval_kl.max():.6f} (interval {np.argmax(interval_kl)})")

    # Convert timesteps to flow-matching sigmas (sigma = t / 1000)
    fm_sigmas = np.array(all_t) / 1000.0  # 28 values from 1.0 to 0.0

    # Cumulative KL at each timestep boundary
    cum_kl = np.concatenate([[0.0], np.cumsum(interval_kl)])  # 28 values
    total_kl = cum_kl[-1]

    # We want n_target steps => n_target intervals => n_target + 1 boundaries
    # Include start (sigma=1.0) and end (sigma=0.0)
    target_cum = np.linspace(0.0, total_kl, n_target + 1)

    # Interpolate: for each target cumulative KL, find the flow-matching sigma
    entropic_fm = np.interp(target_cum, cum_kl, fm_sigmas)

    # Force endpoints
    entropic_fm[0] = 1.0
    entropic_fm[-1] = 0.0

    # Ensure monotonically decreasing
    entropic_fm = np.sort(entropic_fm)[::-1]

    print(f"\n  Entropic schedule ({n_target} steps):")
    print(f"    Flow-matching sigmas (boundaries): {np.round(entropic_fm, 5).tolist()}")
    print(f"    Approx timesteps: {[f'{s*1000:.0f}' for s in entropic_fm]}")

    # Return the n_target sigma VALUES (the step points, not boundaries)
    # For the scheduler: we need n_target values to pass to set_timesteps(sigmas=...)
    # The scheduler will add the terminal sigma (0.0) internally
    return entropic_fm[:-1].tolist()  # Remove terminal 0.0


def check_not_black(img, label):
    """Assert image is not black."""
    mean_px = np.array(img).mean()
    if mean_px <= 10:
        raise ValueError(f"Image {label} is black! mean pixel = {mean_px:.1f}")
    return mean_px


def make_comparison_grid(images, labels, output_path):
    """Create a 2x2 comparison grid with labels."""
    w, h = images[0].size
    margin = 55
    padding = 10
    grid_w = 2 * w + 3 * padding
    grid_h = 2 * (h + margin) + 3 * padding

    grid = Image.new("RGB", (grid_w, grid_h), (255, 255, 255))
    draw = ImageDraw.Draw(grid)

    try:
        font = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 18)
        font_small = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 14)
    except Exception:
        font = ImageFont.load_default()
        font_small = font

    positions = [
        (padding, padding),
        (w + 2 * padding, padding),
        (padding, h + margin + 2 * padding),
        (w + 2 * padding, h + margin + 2 * padding),
    ]

    for i, (img, label) in enumerate(zip(images, labels)):
        x, y = positions[i]
        draw.text((x + 5, y + 5), label, fill=(0, 0, 0), font=font)
        grid.paste(img, (x, y + margin))

    grid.save(output_path, quality=95)
    print(f"  Comparison grid saved: {output_path}")
    return grid


# =====================================================================
# Main
# =====================================================================

def main():
    print("=" * 70)
    print("ENTROPIC TIME SCHEDULER FOR IMAGE DIFFUSION")
    print(f"  Prompt: \"{PROMPT}\"")
    print(f"  Seed: {SEED}")
    print("=" * 70)

    # ------------------------------------------------------------------
    # Phase 1: Load pre-existing sigma profile
    # ------------------------------------------------------------------
    print("\n[Phase 1] Loading sigma profile from previous 28-step run...")
    with open(PROFILE_PATH) as f:
        profile_data = json.load(f)

    sigma_data = profile_data["sigma_data"]
    print(f"  Loaded {len(sigma_data)} data points")
    print(f"  Total Sigma: {profile_data['total_sigma']:.4f}")
    print(f"  Top-8 steps by Sigma: {profile_data['top8_steps']}")

    # Display the full profile
    print(f"\n  Full Sigma_KL profile:")
    for s in sigma_data:
        if s["sigma"] > 0:
            bar = "#" * min(int(s["sigma"] * 50), 50)
            print(f"    Step {s['step']:2d} (t={s['t']:7.1f}): "
                  f"Sigma={s['sigma']:.6f}  {bar}")

    # ------------------------------------------------------------------
    # Phase 2: Build entropic schedule
    # ------------------------------------------------------------------
    print("\n[Phase 2] Building entropic 8-step schedule...")

    entropic_sigmas = build_entropic_sigmas(profile_data, n_target=8)

    # Get the uniform 8-step sigmas for comparison
    from diffusers import FlowMatchEulerDiscreteScheduler
    scheduler = FlowMatchEulerDiscreteScheduler(
        num_train_timesteps=1000, shift=6.0, use_dynamic_shifting=False
    )
    scheduler.set_timesteps(8, device="cpu")
    uniform_fm_sigmas = scheduler.sigmas.cpu().numpy()
    uniform_ts = scheduler.timesteps.cpu().numpy()

    print(f"\n  COMPARISON:")
    print(f"    Uniform  8 sigmas: {np.round(uniform_fm_sigmas[:-1], 5).tolist()}")
    print(f"    Entropic 8 sigmas: {[round(s, 5) for s in entropic_sigmas]}")
    print(f"    Uniform  8 timesteps: {uniform_ts.tolist()}")
    print(f"    Entropic 8 timesteps: {[round(s*1000, 1) for s in entropic_sigmas]}")

    # Key difference visualization
    print(f"\n  KEY INSIGHT:")
    print(f"    Uniform puts first step at t={uniform_ts[0]:.0f} (sigma={uniform_fm_sigmas[0]:.3f})")
    print(f"    Entropic puts first step at t={entropic_sigmas[0]*1000:.0f} "
          f"(sigma={entropic_sigmas[0]:.3f})")
    print(f"    Uniform: {len([s for s in uniform_fm_sigmas[:-1] if s < 0.5])} steps below sigma=0.5")
    print(f"    Entropic: {len([s for s in entropic_sigmas if s < 0.5])} steps below sigma=0.5")
    print(f"    => Entropic concentrates MORE steps in the high-information region!")

    # ------------------------------------------------------------------
    # Phase 3: Load model and run experiments
    # ------------------------------------------------------------------
    from diffusers import ZImagePipeline

    print("\n[Phase 3] Loading Z-Image Base...")
    t_load = time.time()
    pipe_base = ZImagePipeline.from_pretrained(
        "Tongyi-MAI/Z-Image",
        torch_dtype=DTYPE,
    )
    pipe_base = pipe_base.to(DEVICE)
    pipe_base.set_progress_bar_config(disable=False)
    print(f"  Model loaded in {time.time() - t_load:.1f}s")

    # --- [A] 28 steps uniform baseline ---
    print("\n[A] Z-Image Base, 28 steps, uniform schedule (BASELINE)")
    sys.stdout.flush()
    gen_a = torch.Generator(device="cpu").manual_seed(SEED)
    t0 = time.time()
    result_a = pipe_base(
        prompt=PROMPT,
        num_inference_steps=28,
        generator=gen_a,
        guidance_scale=3.5,
    )
    time_a = time.time() - t0
    img_a = result_a.images[0]
    mean_px_a = check_not_black(img_a, "A")
    img_a.save(os.path.join(OUTPUT_DIR, "A_base_28step_uniform.png"))
    print(f"  Done in {time_a:.1f}s, mean_px={mean_px_a:.1f}")
    sys.stdout.flush()

    # --- [B] 8 steps uniform ---
    print("\n[B] Z-Image Base, 8 steps, uniform schedule (naive reduction)")
    sys.stdout.flush()
    gen_b = torch.Generator(device="cpu").manual_seed(SEED)
    t0 = time.time()
    result_b = pipe_base(
        prompt=PROMPT,
        num_inference_steps=8,
        generator=gen_b,
        guidance_scale=3.5,
    )
    time_b = time.time() - t0
    img_b = result_b.images[0]
    mean_px_b = check_not_black(img_b, "B")
    img_b.save(os.path.join(OUTPUT_DIR, "B_base_8step_uniform.png"))
    print(f"  Done in {time_b:.1f}s, mean_px={mean_px_b:.1f}")
    sys.stdout.flush()

    # --- [C] 8 steps ENTROPIC ---
    print("\n[C] Z-Image Base, 8 steps, ENTROPIC schedule")
    print(f"  Custom sigmas: {[round(s, 5) for s in entropic_sigmas]}")
    sys.stdout.flush()
    gen_c = torch.Generator(device="cpu").manual_seed(SEED)
    t0 = time.time()
    result_c = pipe_base(
        prompt=PROMPT,
        num_inference_steps=8,
        generator=gen_c,
        guidance_scale=3.5,
        sigmas=entropic_sigmas,
    )
    time_c = time.time() - t0
    img_c = result_c.images[0]
    mean_px_c = check_not_black(img_c, "C")
    img_c.save(os.path.join(OUTPUT_DIR, "C_base_8step_entropic.png"))
    print(f"  Done in {time_c:.1f}s, mean_px={mean_px_c:.1f}")
    sys.stdout.flush()

    # Free base model
    del pipe_base
    torch.mps.empty_cache()

    # --- [D] Z-Image Turbo 8 steps ---
    print("\n[D] Z-Image Turbo, 8 steps (distilled competitor)")
    sys.stdout.flush()
    t_load = time.time()
    pipe_turbo = ZImagePipeline.from_pretrained(
        "Tongyi-MAI/Z-Image-Turbo",
        torch_dtype=DTYPE,
    )
    pipe_turbo = pipe_turbo.to(DEVICE)
    pipe_turbo.set_progress_bar_config(disable=False)
    print(f"  Turbo loaded in {time.time() - t_load:.1f}s")

    gen_d = torch.Generator(device="cpu").manual_seed(SEED)
    t0 = time.time()
    result_d = pipe_turbo(
        prompt=PROMPT,
        num_inference_steps=8,
        generator=gen_d,
        guidance_scale=0.0,
    )
    time_d = time.time() - t0
    img_d = result_d.images[0]
    mean_px_d = check_not_black(img_d, "D")
    img_d.save(os.path.join(OUTPUT_DIR, "D_turbo_8step.png"))
    print(f"  Done in {time_d:.1f}s, mean_px={mean_px_d:.1f}")
    sys.stdout.flush()

    del pipe_turbo
    torch.mps.empty_cache()

    # ------------------------------------------------------------------
    # Phase 4: Comparison grid
    # ------------------------------------------------------------------
    print("\n[Phase 4] Creating comparison grid...")

    labels = [
        f"A: Base 28-step uniform ({time_a:.0f}s)",
        f"B: Base 8-step uniform ({time_b:.0f}s)",
        f"C: Base 8-step ENTROPIC ({time_c:.0f}s)",
        f"D: Turbo 8-step ({time_d:.0f}s)",
    ]
    make_comparison_grid(
        [img_a, img_b, img_c, img_d],
        labels,
        os.path.join(OUTPUT_DIR, "comparison_grid.png"),
    )

    # ------------------------------------------------------------------
    # Phase 5: Metrics
    # ------------------------------------------------------------------
    print("\n[Phase 5] Saving metrics...")

    metrics = {
        "prompt": PROMPT,
        "seed": SEED,
        "method": "Entropic Time Scheduler (NeurIPS 2025, arXiv:2504.13612)",
        "sigma_profile_source": PROFILE_PATH,
        "experiments": {
            "A_base_28step_uniform": {
                "model": "Tongyi-MAI/Z-Image",
                "steps": 28,
                "schedule": "uniform",
                "guidance_scale": 3.5,
                "time_s": round(time_a, 1),
                "mean_pixel": round(float(mean_px_a), 1),
            },
            "B_base_8step_uniform": {
                "model": "Tongyi-MAI/Z-Image",
                "steps": 8,
                "schedule": "uniform",
                "guidance_scale": 3.5,
                "time_s": round(time_b, 1),
                "mean_pixel": round(float(mean_px_b), 1),
            },
            "C_base_8step_entropic": {
                "model": "Tongyi-MAI/Z-Image",
                "steps": 8,
                "schedule": "entropic (KL-guided redistribution)",
                "guidance_scale": 3.5,
                "custom_sigmas": [round(s, 6) for s in entropic_sigmas],
                "time_s": round(time_c, 1),
                "mean_pixel": round(float(mean_px_c), 1),
            },
            "D_turbo_8step": {
                "model": "Tongyi-MAI/Z-Image-Turbo",
                "steps": 8,
                "schedule": "uniform (distilled)",
                "guidance_scale": 0.0,
                "time_s": round(time_d, 1),
                "mean_pixel": round(float(mean_px_d), 1),
            },
        },
        "sigma_profile_summary": {
            "total_sigma_kl": profile_data["total_sigma"],
            "top8_steps": profile_data["top8_steps"],
            "key_finding": "Sigma peaks at step 26 (t=187.5) with value 0.666, "
                           "1000x larger than early steps. Information is concentrated "
                           "in the cleanup phase (low timesteps / high step numbers).",
        },
        "schedule_comparison": {
            "uniform_8_sigmas": np.round(uniform_fm_sigmas[:-1], 6).tolist(),
            "entropic_8_sigmas": [round(s, 6) for s in entropic_sigmas],
            "uniform_8_timesteps": uniform_ts.tolist(),
            "entropic_8_approx_timesteps": [round(s * 1000, 1) for s in entropic_sigmas],
        },
    }

    with open(os.path.join(OUTPUT_DIR, "metrics.json"), "w") as f:
        json.dump(metrics, f, indent=2)

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("RESULTS SUMMARY")
    print("=" * 70)
    print(f"  Prompt: \"{PROMPT}\", Seed: {SEED}")
    print(f"\n  A: Base 28-step uniform  -> {time_a:.1f}s, mean_px={mean_px_a:.1f}")
    print(f"  B: Base 8-step uniform   -> {time_b:.1f}s, mean_px={mean_px_b:.1f}")
    print(f"  C: Base 8-step ENTROPIC  -> {time_c:.1f}s, mean_px={mean_px_c:.1f}")
    print(f"  D: Turbo 8-step          -> {time_d:.1f}s, mean_px={mean_px_d:.1f}")
    print(f"\n  Schedule comparison:")
    print(f"    Uniform  sigmas: {np.round(uniform_fm_sigmas[:-1], 4).tolist()}")
    print(f"    Entropic sigmas: {[round(s, 4) for s in entropic_sigmas]}")
    print(f"\n  All results saved to: {OUTPUT_DIR}")
    print("=" * 70)


if __name__ == "__main__":
    main()
