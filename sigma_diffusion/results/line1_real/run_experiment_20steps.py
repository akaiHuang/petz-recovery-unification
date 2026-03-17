#!/usr/bin/env python3
"""
Sigma-Guided Diffusion Experiment  (v3 — 20 steps for visual quality)
=====================================================================
Same structure as v2, but with 20 reduced steps (40% of baseline 50).
This should produce visually meaningful images while still showing
the Sigma advantage.
"""
import os
os.environ["HF_TOKEN"] = "YOUR_HF_TOKEN_HERE"
os.environ["PYTORCH_MPS_HIGH_WATERMARK_RATIO"] = "0.0"
os.environ["MPLCONFIGDIR"] = "/tmp/mpl_config"

import json
import time
import warnings
import numpy as np
import torch
from pathlib import Path
from PIL import Image

warnings.filterwarnings("ignore")

RESULTS_DIR = Path(__file__).parent
SEED = 42
REDUCED_STEPS = 20
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

print(f"Config: {BASELINE_STEPS} baseline -> {REDUCED_STEPS} reduced steps ({REDUCED_STEPS/BASELINE_STEPS*100:.0f}%)")

from diffusers import StableDiffusionPipeline, DDIMScheduler

print("Loading model...")
t0 = time.time()
pipe = StableDiffusionPipeline.from_pretrained(
    MODEL_ID, torch_dtype=DTYPE,
    safety_checker=None, requires_safety_checker=False,
)
pipe = pipe.to(DEVICE)
pipe.set_progress_bar_config(disable=True)
print(f"Loaded in {time.time()-t0:.1f}s")

# ─── Helpers ───
def compute_psnr(a, b):
    mse = np.mean((a.astype(np.float64) - b.astype(np.float64))**2)
    return float('inf') if mse == 0 else 10*np.log10(255**2/mse)

def compute_ssim(a, b):
    g1 = np.mean(a.astype(np.float64), axis=2)
    g2 = np.mean(b.astype(np.float64), axis=2)
    mu1, mu2 = g1.mean(), g2.mean()
    C1, C2 = (0.01*255)**2, (0.03*255)**2
    return float(((2*mu1*mu2+C1)*(2*np.mean((g1-mu1)*(g2-mu2))+C2))/((mu1**2+mu2**2+C1)*(g1.var()+g2.var()+C2)))

def gaussian_kl(m1, s1, m2, s2):
    s1, s2 = np.maximum(s1,1e-8), np.maximum(s2,1e-8)
    return float(np.sum(np.log(s2/s1) + (s1**2+(m1-m2)**2)/(2*s2**2) - 0.5))

def run_custom(pipe, prompt, ts_list, seed):
    sched = DDIMScheduler.from_config(pipe.scheduler.config)
    sched.set_timesteps(BASELINE_STEPS)
    with torch.no_grad():
        ti = pipe.tokenizer(prompt, padding="max_length", max_length=pipe.tokenizer.model_max_length, truncation=True, return_tensors="pt")
        te = pipe.text_encoder(ti.input_ids.to(DEVICE))[0]
        ui = pipe.tokenizer("", padding="max_length", max_length=pipe.tokenizer.model_max_length, return_tensors="pt")
        ue = pipe.text_encoder(ui.input_ids.to(DEVICE))[0]
        emb = torch.cat([ue, te])
    gen = torch.Generator(device="cpu").manual_seed(seed)
    lat = torch.randn((1,pipe.unet.config.in_channels,IMG_SIZE//8,IMG_SIZE//8), generator=gen, dtype=DTYPE).to(DEVICE)
    lat = lat * sched.init_noise_sigma
    for t in torch.tensor(ts_list, dtype=torch.long):
        inp = torch.cat([lat]*2)
        inp = sched.scale_model_input(inp, t)
        with torch.no_grad():
            np_ = pipe.unet(inp, t.to(DEVICE), encoder_hidden_states=emb).sample
        u, c = np_.chunk(2)
        np_ = u + GUIDANCE_SCALE*(c-u)
        lat = sched.step(np_, t, lat).prev_sample
    with torch.no_grad():
        img = pipe.vae.decode(lat/pipe.vae.config.scaling_factor).sample
    img = (img/2+0.5).clamp(0,1).detach().cpu().permute(0,2,3,1).float().numpy()
    return Image.fromarray((img[0]*255).round().astype(np.uint8))

# ═══ PHASE 1: Sigma profile ═══
print("\n--- Phase 1: Sigma profiling ---")
pipe.scheduler = DDIMScheduler.from_config(pipe.scheduler.config)
pipe.scheduler.set_timesteps(BASELINE_STEPS)

latent_history = []
def cb(p, si, ts, kw):
    l = kw["latents"].detach().cpu().float().numpy()
    latent_history.append({"step":si, "timestep":int(ts),
        "means":l.mean(axis=(0,2,3)).tolist(), "stds":l.std(axis=(0,2,3)).tolist()})
    return kw

gen = torch.Generator(device="cpu").manual_seed(SEED)
_ = pipe(PROMPTS[0], num_inference_steps=BASELINE_STEPS, generator=gen,
         guidance_scale=GUIDANCE_SCALE, callback_on_step_end=cb)

sigma_profile = []
for i in range(1, len(latent_history)):
    p, c = latent_history[i-1], latent_history[i]
    kl = gaussian_kl(np.array(p["means"]), np.array(p["stds"]),
                     np.array(c["means"]), np.array(c["stds"]))
    sigma_profile.append({"step":c["step"], "timestep":c["timestep"], "sigma_kl":kl})

# Stratified selection: K bins
n = len(sigma_profile)
K = REDUCED_STEPS
bin_size = n // K
strat_ts = []
for b in range(K):
    s, e = b*bin_size, (b+1)*bin_size if b < K-1 else n
    best = max(sigma_profile[s:e], key=lambda x: abs(x["sigma_kl"]))
    strat_ts.append(best["timestep"])
strat_ts = sorted(strat_ts, reverse=True)

# DDIM uniform 20 steps
pipe.scheduler.set_timesteps(REDUCED_STEPS)
ddim_ts = pipe.scheduler.timesteps.tolist()

print(f"DDIM:  {ddim_ts}")
print(f"Strat: {strat_ts}")

# ═══ PHASE 2: Comparisons ═══
print("\n--- Phase 2: Comparisons ---")
all_results = {}
all_images = {}

for pi, prompt in enumerate(PROMPTS):
    pk = f"prompt_{pi+1}"
    print(f"\nPrompt {pi+1}: {prompt[:50]}...")
    all_results[pk] = {"prompt": prompt}
    all_images[pk] = {}

    # Baseline
    print("  Baseline...", end="", flush=True)
    t0 = time.time()
    pipe.scheduler = DDIMScheduler.from_config(pipe.scheduler.config)
    gen = torch.Generator(device="cpu").manual_seed(SEED)
    ib = pipe(prompt, num_inference_steps=BASELINE_STEPS, generator=gen, guidance_scale=GUIDANCE_SCALE).images[0]
    print(f" {time.time()-t0:.1f}s")
    ib.save(RESULTS_DIR / f"p{pi+1}_base.png")
    all_images[pk]["baseline"] = ib

    # DDIM
    print("  DDIM...", end="", flush=True)
    t0 = time.time()
    id_ = run_custom(pipe, prompt, ddim_ts, SEED)
    print(f" {time.time()-t0:.1f}s")
    id_.save(RESULTS_DIR / f"p{pi+1}_ddim20.png")
    all_images[pk]["ddim"] = id_

    # Sigma-strat
    print("  Sigma-strat...", end="", flush=True)
    t0 = time.time()
    is_ = run_custom(pipe, prompt, strat_ts, SEED)
    print(f" {time.time()-t0:.1f}s")
    is_.save(RESULTS_DIR / f"p{pi+1}_sigma20.png")
    all_images[pk]["sigma_strat"] = is_

    # Metrics
    bn = np.array(ib)
    for m, im in [("ddim", id_), ("sigma_strat", is_)]:
        an = np.array(im)
        all_results[pk][f"psnr_{m}"] = round(compute_psnr(bn, an), 3)
        all_results[pk][f"ssim_{m}"] = round(compute_ssim(bn, an), 4)

    best = max(["ddim","sigma_strat"], key=lambda m: all_results[pk][f"psnr_{m}"])
    all_results[pk]["winner"] = best
    print(f"  DDIM:  PSNR={all_results[pk]['psnr_ddim']:.2f} SSIM={all_results[pk]['ssim_ddim']:.4f}")
    print(f"  Sigma: PSNR={all_results[pk]['psnr_sigma_strat']:.2f} SSIM={all_results[pk]['ssim_sigma_strat']:.4f}")
    print(f"  Winner: {best}")

# ═══ Summary ═══
print("\n" + "="*60)
print("SUMMARY (20 steps / 50 baseline)")
print("="*60)

wins = {"ddim":0, "sigma_strat":0}
for pk in [f"prompt_{i+1}" for i in range(len(PROMPTS))]:
    wins[all_results[pk]["winner"]] += 1
    dp = all_results[pk]["psnr_sigma_strat"] - all_results[pk]["psnr_ddim"]
    ds = all_results[pk]["ssim_sigma_strat"] - all_results[pk]["ssim_ddim"]
    print(f"  {pk}: delta PSNR={dp:+.2f} dB, delta SSIM={ds:+.4f} -> {all_results[pk]['winner']}")

avg_dp = np.mean([all_results[f"prompt_{i+1}"]["psnr_sigma_strat"] - all_results[f"prompt_{i+1}"]["psnr_ddim"] for i in range(len(PROMPTS))])
avg_ds = np.mean([all_results[f"prompt_{i+1}"]["ssim_sigma_strat"] - all_results[f"prompt_{i+1}"]["ssim_ddim"] for i in range(len(PROMPTS))])

print(f"\nWins: DDIM={wins['ddim']} Sigma-strat={wins['sigma_strat']}")
print(f"Avg delta: PSNR={avg_dp:+.3f} dB, SSIM={avg_ds:+.4f}")

all_results["summary"] = {
    "wins": {k:int(v) for k,v in wins.items()},
    "avg_delta_psnr": round(float(avg_dp), 3),
    "avg_delta_ssim": round(float(avg_ds), 4),
    "ddim_timesteps": ddim_ts,
    "sigma_timesteps": strat_ts,
    "baseline_steps": BASELINE_STEPS,
    "reduced_steps": REDUCED_STEPS,
}

with open(RESULTS_DIR / "comparison_metrics_20.json", "w") as f:
    json.dump(all_results, f, indent=2, default=str)

# ─── Grid ───
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, axes = plt.subplots(3, 3, figsize=(18, 18))
fig.suptitle(
    f"Sigma-Guided Step Selection: 20/{BASELINE_STEPS} steps\n"
    f"Model: SD v1.5 | Seed: {SEED}",
    fontsize=15, fontweight="bold", y=0.98
)

for pi in range(3):
    pk = f"prompt_{pi+1}"
    r = all_results[pk]
    for mi, (mk, ml) in enumerate(zip(
        ["baseline", "ddim", "sigma_strat"],
        [f"Baseline ({BASELINE_STEPS})", f"DDIM-uniform ({REDUCED_STEPS})", f"Sigma-strat ({REDUCED_STEPS})"]
    )):
        ax = axes[pi][mi]
        ax.imshow(np.array(all_images[pk][mk]))
        ax.axis("off")
        if pi == 0:
            ax.set_title(ml, fontsize=13, fontweight="bold")
        if mk != "baseline":
            p = r[f"psnr_{mk}"]
            s = r[f"ssim_{mk}"]
            w = r["winner"] == mk
            c = "forestgreen" if w else "crimson"
            t = " BEST" if w else ""
            ax.text(5, 25, f"PSNR: {p:.1f}{t}\nSSIM: {s:.3f}",
                   fontsize=11, color="white", fontfamily="monospace",
                   bbox=dict(boxstyle="round,pad=0.3", facecolor=c, alpha=0.8))

fig.text(0.5, 0.01,
         f"Wins: DDIM={wins['ddim']} | Sigma={wins['sigma_strat']}  |  "
         f"Avg delta: PSNR={avg_dp:+.3f} dB, SSIM={avg_ds:+.4f}",
         fontsize=13, ha="center", fontweight="bold",
         bbox=dict(boxstyle="round,pad=0.5", facecolor="lightyellow", alpha=0.9))

plt.tight_layout(rect=[0, 0.03, 1, 0.94])
fig.savefig(RESULTS_DIR / "comparison_grid_20.png", dpi=150, bbox_inches="tight")
plt.close()

print(f"\nSaved: comparison_grid_20.png, comparison_metrics_20.json")
print("DONE")
