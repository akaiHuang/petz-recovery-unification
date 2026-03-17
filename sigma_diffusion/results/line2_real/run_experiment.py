#!/usr/bin/env python3
"""
Run the MDLM Sigma-Guided experiment with practical parameters for M1 Max.
This is the actual experiment runner -- tuned for feasible wall-clock time.
"""

import os
import sys
import time
import json
import math
from pathlib import Path

os.environ["HF_TOKEN"] = "YOUR_HF_TOKEN_HERE"
os.environ["TOKENIZERS_PARALLELISM"] = "false"

import numpy as np
import torch
import transformers

# Add the experiment module
sys.path.insert(0, str(Path(__file__).parent))
from mdlm_mps_experiment import (
    MDLMModel, load_mdlm_weights, MASK_INDEX, DEVICE,
    generate_baseline, generate_uniform_skip, generate_sigma_guided,
    profile_sigma_all_steps, compute_perplexity_gpt2,
    subs_parameterization, loglinear_noise, ddpm_cache_update,
    sample_categorical, compute_step_sigma
)

RESULTS_DIR = Path(__file__).parent


def main():
    print("=" * 70)
    print("REAL MDLM Sigma-Guided Diffusion Experiment")
    print("Model: kuleshov-group/mdlm-owt (NeurIPS 2024)")
    print("Device: Apple MPS (M1 Max)")
    print("=" * 70)

    # ── Load model ──
    print("\n[1/6] Loading pretrained MDLM model...")
    model = MDLMModel(
        vocab_size=50258, hidden_dim=768, cond_dim=128,
        n_blocks=12, n_heads=12, dropout=0.1)
    load_mdlm_weights(model)
    model = model.to(DEVICE).eval()
    n_params = sum(p.numel() for p in model.parameters())
    print(f"  Parameters: {n_params/1e6:.1f}M on {DEVICE}")

    tokenizer = transformers.AutoTokenizer.from_pretrained("gpt2")
    if tokenizer.pad_token is None:
        tokenizer.pad_token = tokenizer.eos_token

    # ── Experiment parameters ──
    # Chosen to be feasible on M1 Max in ~10 minutes total
    FULL_STEPS = 128       # baseline step count
    REDUCED_STEPS = 64     # 50% step reduction
    SEQ_LEN = 256          # moderate sequence length
    BATCH_SIZE = 4         # enough samples for statistics
    N_RUNS = 3             # repeat for stability

    print(f"\n  Parameters:")
    print(f"    Full steps:    {FULL_STEPS}")
    print(f"    Reduced steps: {REDUCED_STEPS} (50% reduction)")
    print(f"    Seq length:    {SEQ_LEN}")
    print(f"    Batch size:    {BATCH_SIZE}")
    print(f"    Repeats:       {N_RUNS}")

    # ── Warm-up ──
    print("\n[2/6] Warm-up...")
    dummy = MASK_INDEX * torch.ones(1, 64, dtype=torch.int64, device=DEVICE)
    _ = model(dummy, torch.tensor([0.5], device=DEVICE))
    print("  Done.")

    # ── Phase 1: Profile Sigma ──
    print(f"\n[3/6] Profiling Sigma over {FULL_STEPS} steps...")
    sigma_profile = profile_sigma_all_steps(
        model, FULL_STEPS, seq_len=SEQ_LEN, batch_size=1, device=DEVICE)

    sigmas = [s['sigma'] for s in sigma_profile]
    print(f"\n  Sigma statistics over {FULL_STEPS} steps:")
    print(f"    Min:    {min(sigmas):.8f}")
    print(f"    Max:    {max(sigmas):.8f}")
    print(f"    Mean:   {np.mean(sigmas):.8f}")
    print(f"    Std:    {np.std(sigmas):.8f}")
    print(f"    Median: {np.median(sigmas):.8f}")

    # Which steps have highest Sigma?
    sorted_steps = sorted(sigma_profile, key=lambda s: s['sigma'], reverse=True)
    print(f"\n  Top-10 highest Sigma steps:")
    for s in sorted_steps[:10]:
        print(f"    Step {s['step']:3d} (t={s['timestep']:.4f}): Sigma={s['sigma']:.8f}")
    print(f"\n  Bottom-10 lowest Sigma steps:")
    for s in sorted_steps[-10:]:
        print(f"    Step {s['step']:3d} (t={s['timestep']:.4f}): Sigma={s['sigma']:.8f}")

    # ── Phase 2: Generate and compare ──
    all_runs = {'baseline': [], 'uniform_skip': [], 'sigma_guided': []}

    for run_idx in range(N_RUNS):
        print(f"\n[4/6] Generation run {run_idx+1}/{N_RUNS}...")

        # A: Baseline
        print(f"  A. Baseline ({FULL_STEPS} steps)...")
        tokens_a, time_a, _ = generate_baseline(
            model, FULL_STEPS, seq_len=SEQ_LEN,
            batch_size=BATCH_SIZE, device=DEVICE)
        text_a = tokenizer.batch_decode(tokens_a.cpu(), skip_special_tokens=True)
        ppl_a, ppl_a_per = compute_perplexity_gpt2(text_a, device=DEVICE)
        print(f"    Time: {time_a:.2f}s, PPL: {ppl_a:.2f}")
        all_runs['baseline'].append({
            'time': time_a, 'ppl': ppl_a, 'ppl_per': ppl_a_per, 'text': text_a})

        # B: Uniform skip
        print(f"  B. Uniform skip ({REDUCED_STEPS} steps)...")
        tokens_b, time_b = generate_uniform_skip(
            model, FULL_STEPS, REDUCED_STEPS,
            seq_len=SEQ_LEN, batch_size=BATCH_SIZE, device=DEVICE)
        text_b = tokenizer.batch_decode(tokens_b.cpu(), skip_special_tokens=True)
        ppl_b, ppl_b_per = compute_perplexity_gpt2(text_b, device=DEVICE)
        print(f"    Time: {time_b:.2f}s, PPL: {ppl_b:.2f}")
        all_runs['uniform_skip'].append({
            'time': time_b, 'ppl': ppl_b, 'ppl_per': ppl_b_per, 'text': text_b})

        # C: Sigma-guided
        print(f"  C. Sigma-guided ({REDUCED_STEPS} of {FULL_STEPS} steps)...")
        tokens_c, time_c = generate_sigma_guided(
            model, FULL_STEPS, REDUCED_STEPS,
            sigma_profile, seq_len=SEQ_LEN,
            batch_size=BATCH_SIZE, device=DEVICE)
        text_c = tokenizer.batch_decode(tokens_c.cpu(), skip_special_tokens=True)
        ppl_c, ppl_c_per = compute_perplexity_gpt2(text_c, device=DEVICE)
        print(f"    Time: {time_c:.2f}s, PPL: {ppl_c:.2f}")
        all_runs['sigma_guided'].append({
            'time': time_c, 'ppl': ppl_c, 'ppl_per': ppl_c_per, 'text': text_c})

    # ── Aggregate results ──
    print(f"\n[5/6] Aggregating results over {N_RUNS} runs...")
    summary = {}
    for method, runs in all_runs.items():
        ppls = [r['ppl'] for r in runs]
        times = [r['time'] for r in runs]
        summary[method] = {
            'ppl_mean': np.mean(ppls),
            'ppl_std': np.std(ppls),
            'time_mean': np.mean(times),
            'time_std': np.std(times),
            'ppls': ppls,
            'times': times,
        }

    # ── Print summary ──
    print("\n" + "=" * 70)
    print("FINAL RESULTS (averaged over {} runs)".format(N_RUNS))
    print("=" * 70)
    print(f"{'Method':<20} {'Steps':<8} {'Time(s)':<12} {'PPL':<15} {'Speedup':<10}")
    print("-" * 65)
    base_time = summary['baseline']['time_mean']
    for method, s in summary.items():
        steps = FULL_STEPS if method == 'baseline' else REDUCED_STEPS
        speedup = base_time / s['time_mean'] if s['time_mean'] > 0 else 0
        print(f"{method:<20} {steps:<8} "
              f"{s['time_mean']:.2f}+/-{s['time_std']:.2f}  "
              f"{s['ppl_mean']:.2f}+/-{s['ppl_std']:.2f}  "
              f"{speedup:.2f}x")

    # ── Verdict ──
    ppl_base = summary['baseline']['ppl_mean']
    ppl_uni = summary['uniform_skip']['ppl_mean']
    ppl_sig = summary['sigma_guided']['ppl_mean']

    print("\n" + "-" * 65)
    print("ANALYSIS:")
    print(f"  Baseline PPL:      {ppl_base:.2f}")
    print(f"  Uniform skip PPL:  {ppl_uni:.2f} (delta: {ppl_uni - ppl_base:+.2f})")
    print(f"  Sigma-guided PPL:  {ppl_sig:.2f} (delta: {ppl_sig - ppl_base:+.2f})")

    if ppl_sig < ppl_uni:
        delta = ppl_uni - ppl_sig
        pct = delta / ppl_uni * 100
        print(f"\n  >>> Sigma-guided WINS over uniform skip by {delta:.2f} PPL ({pct:.1f}%)")
    elif ppl_sig > ppl_uni:
        delta = ppl_sig - ppl_uni
        pct = delta / ppl_uni * 100
        print(f"\n  >>> Sigma-guided LOSES to uniform skip by {delta:.2f} PPL ({pct:.1f}%)")
    else:
        print(f"\n  >>> TIE")

    # Is the result within noise?
    combined_std = np.sqrt(summary['uniform_skip']['ppl_std']**2 +
                           summary['sigma_guided']['ppl_std']**2)
    if combined_std > 0:
        z_score = abs(ppl_sig - ppl_uni) / combined_std
        print(f"  Effect size: {z_score:.2f} sigma (combined std={combined_std:.2f})")
        if z_score < 1:
            print(f"  NOTE: Difference is within 1-sigma noise -- NOT statistically significant")
        elif z_score < 2:
            print(f"  NOTE: Difference is marginal (1-2 sigma)")
        else:
            print(f"  NOTE: Difference appears statistically significant (>{z_score:.1f} sigma)")

    # ── Show sample text ──
    print("\n" + "-" * 65)
    print("SAMPLE TEXT (first sample from last run):")
    for method in ['baseline', 'uniform_skip', 'sigma_guided']:
        text = all_runs[method][-1]['text'][0][:300]
        print(f"\n  [{method}]:")
        print(f"  {text}")

    # ── Save everything ──
    print(f"\n[6/6] Saving results...")

    # Sigma profile
    with open(RESULTS_DIR / "sigma_profile.json", 'w') as f:
        json.dump(sigma_profile, f, indent=2)

    # Summary
    with open(RESULTS_DIR / "experiment_summary.json", 'w') as f:
        json.dump({
            'config': {
                'full_steps': FULL_STEPS,
                'reduced_steps': REDUCED_STEPS,
                'seq_len': SEQ_LEN,
                'batch_size': BATCH_SIZE,
                'n_runs': N_RUNS,
                'model': 'kuleshov-group/mdlm-owt',
                'device': DEVICE,
                'model_params_M': n_params / 1e6,
            },
            'sigma_stats': {
                'min': float(min(sigmas)),
                'max': float(max(sigmas)),
                'mean': float(np.mean(sigmas)),
                'std': float(np.std(sigmas)),
                'median': float(np.median(sigmas)),
            },
            'results': {
                method: {
                    'ppl_mean': float(s['ppl_mean']),
                    'ppl_std': float(s['ppl_std']),
                    'time_mean': float(s['time_mean']),
                    'time_std': float(s['time_std']),
                }
                for method, s in summary.items()
            },
            'sigma_guided_wins': bool(ppl_sig < ppl_uni),
            'verdict': 'sigma_guided_wins' if ppl_sig < ppl_uni else
                       'uniform_wins' if ppl_sig > ppl_uni else 'tie',
        }, f, indent=2)

    # All run data
    serializable_runs = {}
    for method, runs in all_runs.items():
        serializable_runs[method] = [
            {'time': r['time'], 'ppl': r['ppl'], 'ppl_per': r['ppl_per'],
             'text': r['text']}
            for r in runs
        ]
    with open(RESULTS_DIR / "all_runs.json", 'w') as f:
        json.dump(serializable_runs, f, indent=2, default=str)

    print(f"  Saved to {RESULTS_DIR}/")
    print("  Files: sigma_profile.json, experiment_summary.json, all_runs.json")
    print("\nDone.")


if __name__ == "__main__":
    main()
