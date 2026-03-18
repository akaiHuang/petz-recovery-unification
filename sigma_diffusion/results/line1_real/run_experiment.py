#!/usr/bin/env python3
"""
Sigma-Guided Diffusion Experiment  (v2 — improved selection)
=============================================================
Compare: Baseline (50 steps) vs DDIM-uniform (8 steps) vs Sigma-guided (8 steps)

Key insight from v1: naive "pick top-K by Sigma" clusters all selected steps
in the mid-range, missing critical early (coarse structure) and late (detail)
steps.  v2 uses a STRATIFIED selection: divide the full schedule into K bins
and pick the highest-Sigma step within each bin.  This ensures coverage while
still prioritizing informative transitions.
"""
import os
os.environ["HF_TOKEN"] = "YOUR_HF_TOKEN_HERE"
os.environ["PYTORCH_MPS_HIGH_WATERMARK_RATIO"] = "0.0"

import json
import time
import warnings
import numpy as np
import torch
from pathlib import Path
from PIL import Image

warnings.filterwarnings("ignore")

# ─── Config ───
RESULTS_DIR = Path(__file__).parent
SEED = 42
REDUCED_STEPS = 8
BASELINE_STEPS = 50
MODEL_ID = "stable-diffusion-v1-5/stable-diffusion-v1-5"
DTYPE = torch.float16
DEVICE = "mps"
GUIDANCE_SCALE = 7.5
IMG_SIZE = 512

PROMPTS = [
    "A photorealistic red fox sitting in a snowy forest, golden hour light",
    "An astronaut riding a horse on Mars, cinematic lighting, 4k",
    "A futuristic city skyline at sunset, reflections on water, concept art",
]

print(f"Device: {DEVICE}")
print(f"Model: {MODEL_ID}")
print(f"Baseline steps: {BASELINE_STEPS}, Reduced steps: {REDUCED_STEPS}")
print(f"Prompts: {len(PROMPTS)}")
print()

# ─── Load Model ───
from diffusers import StableDiffusionPipeline, DDIMScheduler

print("Loading model...")
t0 = time.time()
pipe = StableDiffusionPipeline.from_pretrained(
    MODEL_ID, torch_dtype=DTYPE,
    safety_checker=None, requires_safety_checker=False,
)
pipe = pipe.to(DEVICE)
pipe.set_progress_bar_config(disable=True)
print(f"Model loaded in {time.time()-t0:.1f}s")

# ─── Helpers ───
def compute_psnr(img1, img2):
    mse = np.mean((img1.astype(np.float64) - img2.astype(np.float64)) ** 2)
    if mse == 0: return float('inf')
    return 10 * np.log10(255.0**2 / mse)

def compute_ssim(img1, img2):
    g1 = np.mean(img1.astype(np.float64), axis=2)
    g2 = np.mean(img2.astype(np.float64), axis=2)
    mu1, mu2 = g1.mean(), g2.mean()
    sig1_sq, sig2_sq = g1.var(), g2.var()
    sig12 = np.mean((g1 - mu1) * (g2 - mu2))
    C1, C2 = (0.01*255)**2, (0.03*255)**2
    return float(((2*mu1*mu2+C1)*(2*sig12+C2))/((mu1**2+mu2**2+C1)*(sig1_sq+sig2_sq+C2)))

def gaussian_kl(mu1, sig1, mu2, sig2):
    sig1 = np.maximum(sig1, 1e-8)
    sig2 = np.maximum(sig2, 1e-8)
    kl = np.log(sig2/sig1) + (sig1**2 + (mu1-mu2)**2)/(2*sig2**2) - 0.5
    return float(np.sum(kl))

def run_with_custom_timesteps(pipe, prompt, timesteps_list, seed, guidance_scale=7.5):
    """Run diffusion with specific timesteps."""
    scheduler = DDIMScheduler.from_config(pipe.scheduler.config)
    scheduler.set_timesteps(BASELINE_STEPS)

    with torch.no_grad():
        text_input = pipe.tokenizer(prompt, padding="max_length",
            max_length=pipe.tokenizer.model_max_length, truncation=True, return_tensors="pt")
        text_emb = pipe.text_encoder(text_input.input_ids.to(DEVICE))[0]
        uncond_input = pipe.tokenizer("", padding="max_length",
            max_length=pipe.tokenizer.model_max_length, return_tensors="pt")
        uncond_emb = pipe.text_encoder(uncond_input.input_ids.to(DEVICE))[0]
        text_embeddings = torch.cat([uncond_emb, text_emb])

    gen = torch.Generator(device="cpu").manual_seed(seed)
    latent_shape = (1, pipe.unet.config.in_channels, IMG_SIZE//8, IMG_SIZE//8)
    latents = torch.randn(latent_shape, generator=gen, dtype=DTYPE).to(DEVICE)
    latents = latents * scheduler.init_noise_sigma

    ts_tensor = torch.tensor(timesteps_list, dtype=torch.long)
    for t in ts_tensor:
        latent_model_input = torch.cat([latents]*2)
        latent_model_input = scheduler.scale_model_input(latent_model_input, t)
        with torch.no_grad():
            noise_pred = pipe.unet(latent_model_input, t.to(DEVICE),
                encoder_hidden_states=text_embeddings).sample
        u, c = noise_pred.chunk(2)
        noise_pred = u + guidance_scale * (c - u)
        latents = scheduler.step(noise_pred, t, latents).prev_sample

    with torch.no_grad():
        image = pipe.vae.decode(latents / pipe.vae.config.scaling_factor).sample
    image = (image/2+0.5).clamp(0,1)
    image_np = image.detach().cpu().permute(0,2,3,1).float().numpy()
    image_np = (image_np*255).round().astype(np.uint8)
    return Image.fromarray(image_np[0])

# ═════════════════════════════════════════════════════════
# PHASE 1: Sigma Profiling
# ═════════════════════════════════════════════════════════
print("=" * 60)
print("PHASE 1: Collecting Sigma profile from baseline run")
print("=" * 60)

pipe.scheduler = DDIMScheduler.from_config(pipe.scheduler.config)
pipe.scheduler.set_timesteps(BASELINE_STEPS)
baseline_timesteps = pipe.scheduler.timesteps.tolist()

latent_history = []

def record_callback(pipe_obj, step_idx, timestep, callback_kwargs):
    lat = callback_kwargs["latents"].detach().cpu().float().numpy()
    latent_history.append({
        "step": step_idx,
        "timestep": int(timestep),
        "means": lat.mean(axis=(0,2,3)).tolist(),
        "stds": lat.std(axis=(0,2,3)).tolist(),
    })
    return callback_kwargs

print("Running baseline (50 steps) to collect latent statistics...")
t0 = time.time()
gen = torch.Generator(device="cpu").manual_seed(SEED)
_ = pipe(PROMPTS[0], num_inference_steps=BASELINE_STEPS, generator=gen,
         guidance_scale=GUIDANCE_SCALE, callback_on_step_end=record_callback)
print(f"  Done in {time.time()-t0:.1f}s, {len(latent_history)} snapshots")

# Compute per-step Sigma
sigma_profile = []
for i in range(1, len(latent_history)):
    prev, curr = latent_history[i-1], latent_history[i]
    kl = gaussian_kl(np.array(prev["means"]), np.array(prev["stds"]),
                     np.array(curr["means"]), np.array(curr["stds"]))
    sigma_profile.append({
        "step": curr["step"],
        "timestep": curr["timestep"],
        "sigma_kl": kl,
    })

# ─── Selection Strategies ───
# Strategy A: Naive top-K (for reference)
naive_sorted = sorted(sigma_profile, key=lambda x: abs(x["sigma_kl"]), reverse=True)
naive_timesteps = sorted([s["timestep"] for s in naive_sorted[:REDUCED_STEPS]], reverse=True)

# Strategy B: Stratified — divide into K bins, pick max-Sigma per bin
n = len(sigma_profile)
bin_size = n // REDUCED_STEPS
stratified_timesteps = []
for b in range(REDUCED_STEPS):
    start = b * bin_size
    end = start + bin_size if b < REDUCED_STEPS-1 else n
    bin_entries = sigma_profile[start:end]
    best = max(bin_entries, key=lambda x: abs(x["sigma_kl"]))
    stratified_timesteps.append(best["timestep"])
stratified_timesteps = sorted(stratified_timesteps, reverse=True)

# DDIM uniform
pipe.scheduler.set_timesteps(REDUCED_STEPS)
ddim_timesteps = pipe.scheduler.timesteps.tolist()

print(f"\nDDIM-uniform timesteps:     {ddim_timesteps}")
print(f"Sigma-naive timesteps:      {naive_timesteps}")
print(f"Sigma-stratified timesteps: {stratified_timesteps}")

print(f"\nFull Sigma profile:")
for s in sorted(sigma_profile, key=lambda x: x["step"]):
    tags = []
    if s["timestep"] in naive_timesteps: tags.append("NAIVE")
    if s["timestep"] in stratified_timesteps: tags.append("STRAT")
    if s["timestep"] in ddim_timesteps: tags.append("DDIM")
    tag = f"  <-- {', '.join(tags)}" if tags else ""
    print(f"  step={s['step']:3d}  t={s['timestep']:4d}  Sigma={s['sigma_kl']:.6f}{tag}")

# Save sigma profile
with open(RESULTS_DIR / "sigma_profile.json", "w") as f:
    json.dump({
        "prompt": PROMPTS[0],
        "baseline_steps": BASELINE_STEPS,
        "reduced_steps": REDUCED_STEPS,
        "all_sigma_values": sorted(sigma_profile, key=lambda x: x["step"]),
        "ddim_uniform_timesteps": ddim_timesteps,
        "sigma_naive_timesteps": naive_timesteps,
        "sigma_stratified_timesteps": stratified_timesteps,
    }, f, indent=2)
print(f"Saved: sigma_profile.json")

# ═════════════════════════════════════════════════════════
# PHASE 2: Four-way comparison on all prompts
# ═════════════════════════════════════════════════════════
print("\n" + "=" * 60)
print("PHASE 2: Four-way comparison on 3 prompts")
print("  [A] Baseline (50 steps)")
print("  [B] DDIM-uniform (8 steps)")
print("  [C] Sigma-naive (8 steps, top-K)")
print("  [D] Sigma-stratified (8 steps, binned)")
print("=" * 60)

methods_info = {
    "baseline": {"steps": BASELINE_STEPS, "ts": None},
    "ddim": {"steps": REDUCED_STEPS, "ts": ddim_timesteps},
    "sigma_naive": {"steps": REDUCED_STEPS, "ts": naive_timesteps},
    "sigma_strat": {"steps": REDUCED_STEPS, "ts": stratified_timesteps},
}

all_results = {}
all_images = {}

for pi, prompt in enumerate(PROMPTS):
    pk = f"prompt_{pi+1}"
    print(f"\n{'─'*55}")
    print(f"Prompt {pi+1}/{len(PROMPTS)}: {prompt[:55]}...")
    print(f"{'─'*55}")

    all_results[pk] = {"prompt": prompt}
    all_images[pk] = {}

    # Baseline
    print(f"  [A] Baseline ({BASELINE_STEPS} steps)...", end="", flush=True)
    t0 = time.time()
    pipe.scheduler = DDIMScheduler.from_config(pipe.scheduler.config)
    gen = torch.Generator(device="cpu").manual_seed(SEED)
    img_base = pipe(prompt, num_inference_steps=BASELINE_STEPS, generator=gen,
                    guidance_scale=GUIDANCE_SCALE).images[0]
    print(f" {time.time()-t0:.1f}s")
    img_base.save(RESULTS_DIR / f"prompt{pi+1}_baseline.png")
    all_images[pk]["baseline"] = img_base

    # DDIM uniform
    print(f"  [B] DDIM-uniform ({REDUCED_STEPS} steps)...", end="", flush=True)
    t0 = time.time()
    img_ddim = run_with_custom_timesteps(pipe, prompt, ddim_timesteps, SEED, GUIDANCE_SCALE)
    print(f" {time.time()-t0:.1f}s")
    img_ddim.save(RESULTS_DIR / f"prompt{pi+1}_ddim8.png")
    all_images[pk]["ddim"] = img_ddim

    # Sigma naive
    print(f"  [C] Sigma-naive ({REDUCED_STEPS} steps)...", end="", flush=True)
    t0 = time.time()
    img_naive = run_with_custom_timesteps(pipe, prompt, naive_timesteps, SEED, GUIDANCE_SCALE)
    print(f" {time.time()-t0:.1f}s")
    img_naive.save(RESULTS_DIR / f"prompt{pi+1}_sigma_naive8.png")
    all_images[pk]["sigma_naive"] = img_naive

    # Sigma stratified
    print(f"  [D] Sigma-stratified ({REDUCED_STEPS} steps)...", end="", flush=True)
    t0 = time.time()
    img_strat = run_with_custom_timesteps(pipe, prompt, stratified_timesteps, SEED, GUIDANCE_SCALE)
    print(f" {time.time()-t0:.1f}s")
    img_strat.save(RESULTS_DIR / f"prompt{pi+1}_sigma_strat8.png")
    all_images[pk]["sigma_strat"] = img_strat

    # Metrics
    base_np = np.array(img_base)
    for method, img in [("ddim", img_ddim), ("sigma_naive", img_naive), ("sigma_strat", img_strat)]:
        arr = np.array(img)
        psnr = compute_psnr(base_np, arr)
        ssim = compute_ssim(base_np, arr)
        all_results[pk][f"psnr_{method}"] = round(float(psnr), 3)
        all_results[pk][f"ssim_{method}"] = round(float(ssim), 4)

    # Report
    print(f"\n  Results vs baseline:")
    for m in ["ddim", "sigma_naive", "sigma_strat"]:
        p = all_results[pk][f"psnr_{m}"]
        s = all_results[pk][f"ssim_{m}"]
        print(f"    {m:16s}  PSNR={p:6.2f} dB   SSIM={s:.4f}")

    # Winner
    best_method = max(["ddim", "sigma_naive", "sigma_strat"],
                      key=lambda m: all_results[pk][f"psnr_{m}"])
    all_results[pk]["winner_psnr"] = best_method
    print(f"    WINNER (PSNR): {best_method}")

# ═════════════════════════════════════════════════════════
# PHASE 3: Summary
# ═════════════════════════════════════════════════════════
print("\n" + "=" * 60)
print("PHASE 3: Summary")
print("=" * 60)

# Count wins
wins = {"ddim": 0, "sigma_naive": 0, "sigma_strat": 0}
for pk in [f"prompt_{i+1}" for i in range(len(PROMPTS))]:
    wins[all_results[pk]["winner_psnr"]] += 1

print(f"\nWins (PSNR):  DDIM={wins['ddim']}  Naive={wins['sigma_naive']}  Strat={wins['sigma_strat']}  (out of {len(PROMPTS)})")

# Averages
for m in ["ddim", "sigma_naive", "sigma_strat"]:
    avg_p = np.mean([all_results[f"prompt_{i+1}"][f"psnr_{m}"] for i in range(len(PROMPTS))])
    avg_s = np.mean([all_results[f"prompt_{i+1}"][f"ssim_{m}"] for i in range(len(PROMPTS))])
    print(f"  {m:16s}  avg PSNR={avg_p:.2f} dB  avg SSIM={avg_s:.4f}")

# Pairwise deltas: strat vs ddim
deltas_psnr = []
deltas_ssim = []
for i in range(len(PROMPTS)):
    pk = f"prompt_{i+1}"
    dp = all_results[pk]["psnr_sigma_strat"] - all_results[pk]["psnr_ddim"]
    ds = all_results[pk]["ssim_sigma_strat"] - all_results[pk]["ssim_ddim"]
    deltas_psnr.append(dp)
    deltas_ssim.append(ds)
    print(f"  Prompt {i+1}: Sigma-strat vs DDIM: PSNR delta={dp:+.2f} dB, SSIM delta={ds:+.4f}")

avg_dp = np.mean(deltas_psnr)
avg_ds = np.mean(deltas_ssim)
print(f"\n  Average Sigma-strat vs DDIM: PSNR={avg_dp:+.2f} dB, SSIM={avg_ds:+.4f}")

all_results["summary"] = {
    "wins": {k: int(v) for k, v in wins.items()},
    "avg_delta_strat_vs_ddim_psnr": round(float(avg_dp), 3),
    "avg_delta_strat_vs_ddim_ssim": round(float(avg_ds), 4),
    "ddim_timesteps": ddim_timesteps,
    "sigma_naive_timesteps": naive_timesteps,
    "sigma_stratified_timesteps": stratified_timesteps,
    "model": MODEL_ID,
    "baseline_steps": BASELINE_STEPS,
    "reduced_steps": REDUCED_STEPS,
    "seed": SEED,
}

with open(RESULTS_DIR / "comparison_metrics.json", "w") as f:
    json.dump(all_results, f, indent=2, default=str)
print(f"\nSaved: comparison_metrics.json")

# ─── Visualization: 3x4 grid ───
print("\nCreating comparison grid...")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, axes = plt.subplots(3, 4, figsize=(24, 18))
fig.suptitle(
    "Sigma-Guided Diffusion Step Selection\n"
    f"Model: SD v1.5 | Seed: {SEED} | Baseline: {BASELINE_STEPS} steps | Reduced: {REDUCED_STEPS} steps",
    fontsize=15, fontweight="bold", y=0.98
)

method_keys = ["baseline", "ddim", "sigma_naive", "sigma_strat"]
method_labels = [
    f"Baseline\n({BASELINE_STEPS} steps)",
    f"DDIM-uniform\n({REDUCED_STEPS} steps)",
    f"Sigma-naive\n({REDUCED_STEPS} steps, top-K)",
    f"Sigma-stratified\n({REDUCED_STEPS} steps, binned)"
]

for pi in range(3):
    pk = f"prompt_{pi+1}"
    r = all_results[pk]
    for mi, (mk, ml) in enumerate(zip(method_keys, method_labels)):
        ax = axes[pi][mi]
        ax.imshow(np.array(all_images[pk][mk]))
        ax.axis("off")
        if pi == 0:
            ax.set_title(ml, fontsize=12, fontweight="bold")
        if mk != "baseline":
            psnr = r[f"psnr_{mk}"]
            ssim = r[f"ssim_{mk}"]
            is_winner = (r["winner_psnr"] == mk)
            color = "forestgreen" if is_winner else "crimson"
            tag = " BEST" if is_winner else ""
            ax.text(5, 25, f"PSNR: {psnr:.1f} dB{tag}\nSSIM: {ssim:.3f}",
                   fontsize=10, color="white", fontfamily="monospace",
                   bbox=dict(boxstyle="round,pad=0.3", facecolor=color, alpha=0.8))

    axes[pi][0].set_ylabel(f"P{pi+1}", fontsize=14, fontweight="bold", rotation=0, labelpad=25)

fig.text(0.5, 0.01,
         f"Wins: DDIM={wins['ddim']} | Sigma-naive={wins['sigma_naive']} | Sigma-strat={wins['sigma_strat']}  |  "
         f"Avg strat-vs-DDIM: PSNR={avg_dp:+.2f} dB, SSIM={avg_ds:+.4f}",
         fontsize=13, ha="center", fontweight="bold",
         bbox=dict(boxstyle="round,pad=0.5", facecolor="lightyellow", alpha=0.9))

plt.tight_layout(rect=[0.02, 0.03, 1, 0.94])
fig.savefig(RESULTS_DIR / "comparison_grid.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"Saved: comparison_grid.png")

# ─── Sigma profile plot ───
print("Creating Sigma profile plot...")
fig2, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 8), gridspec_kw={"height_ratios": [3, 1]})

steps = [s["step"] for s in sorted(sigma_profile, key=lambda x: x["step"])]
sigmas = [s["sigma_kl"] for s in sorted(sigma_profile, key=lambda x: x["step"])]

# Color bars
bar_colors = []
for s in sorted(sigma_profile, key=lambda x: x["step"]):
    t = s["timestep"]
    if t in stratified_timesteps and t in ddim_timesteps:
        bar_colors.append("gold")
    elif t in stratified_timesteps:
        bar_colors.append("forestgreen")
    elif t in ddim_timesteps:
        bar_colors.append("crimson")
    else:
        bar_colors.append("lightgray")

ax1.bar(range(len(steps)), sigmas, color=bar_colors, edgecolor="gray", linewidth=0.5)
ax1.set_ylabel("Sigma (KL divergence)", fontsize=12)
ax1.set_title("Per-Step Sigma Profile with Selected Timesteps", fontsize=13, fontweight="bold")
ax1.legend(handles=[
    plt.Rectangle((0,0),1,1, fc="forestgreen", label="Sigma-stratified"),
    plt.Rectangle((0,0),1,1, fc="crimson", label="DDIM-uniform"),
    plt.Rectangle((0,0),1,1, fc="gold", label="Both"),
    plt.Rectangle((0,0),1,1, fc="lightgray", label="Not selected"),
], loc="upper left", fontsize=9)

for i, s in enumerate(sorted(sigma_profile, key=lambda x: x["step"])):
    if s["timestep"] in stratified_timesteps:
        ax1.annotate(f't={s["timestep"]}', (i, s["sigma_kl"]),
                    textcoords="offset points", xytext=(0, 5),
                    fontsize=7, ha="center", fontweight="bold", color="darkgreen")

# Bottom: timestep number line showing selections
ts_all = [s["timestep"] for s in sorted(sigma_profile, key=lambda x: x["step"])]
ax2.scatter(range(len(ts_all)), [0]*len(ts_all), c="lightgray", s=20, zorder=1)
for i, t in enumerate(ts_all):
    if t in stratified_timesteps:
        ax2.scatter(i, 0, c="forestgreen", s=80, zorder=3, marker="^")
    if t in ddim_timesteps:
        ax2.scatter(i, -0.1, c="crimson", s=80, zorder=3, marker="v")
ax2.set_xlabel("Denoising Step", fontsize=12)
ax2.set_yticks([-0.1, 0])
ax2.set_yticklabels(["DDIM", "Sigma-strat"], fontsize=10)
ax2.set_ylim(-0.25, 0.15)

fig2.tight_layout()
fig2.savefig(RESULTS_DIR / "sigma_profile_plot.png", dpi=150)
plt.close()
print(f"Saved: sigma_profile_plot.png")

print("\n" + "=" * 60)
print("EXPERIMENT COMPLETE")
print("=" * 60)
print(f"Results directory: {RESULTS_DIR}")
for f in sorted(RESULTS_DIR.glob("*.png")) + sorted(RESULTS_DIR.glob("*.json")):
    print(f"  {f.name}")
