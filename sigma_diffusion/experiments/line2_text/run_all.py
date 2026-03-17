#!/usr/bin/env python3
"""
Line 2: Diffusion LLM + Sigma Optimization
Goal: 1000 steps → 100 steps (10X) for text generation

Uses a minimal discrete diffusion approach over GPT2 tokenizer.
Falls back gracefully if MDLM is not available.
"""

import sys
import os
import json
import time
import numpy as np
from pathlib import Path


class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super().default(obj)

# Add project root
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

RESULTS_DIR = ROOT / "sigma_diffusion" / "results" / "line2"


def check_deps():
    """Check available dependencies."""
    status = {}
    try:
        import torch
        status["torch"] = torch.__version__
        status["mps"] = torch.backends.mps.is_available()
    except ImportError:
        status["torch"] = None
        return status

    try:
        from transformers import GPT2Tokenizer, GPT2LMHeadModel
        status["transformers"] = True
    except ImportError:
        status["transformers"] = False

    return status


def run_discrete_diffusion_baseline(
    n_samples: int = 10,
    seq_len: int = 64,
    T: int = 200,
    device: str = "mps",
):
    """
    Minimal discrete diffusion text generation with Sigma tracking.

    Uses absorbing-state diffusion: tokens are progressively masked,
    then unmasked by a pretrained GPT2 model predicting masked positions.
    """
    import torch
    from transformers import GPT2Tokenizer, GPT2LMHeadModel

    print(f"\n{'='*60}")
    print(f"  Line 2: Discrete Diffusion LLM — Baseline ({T} steps)")
    print(f"{'='*60}\n")

    tokenizer = GPT2Tokenizer.from_pretrained("gpt2")
    model = GPT2LMHeadModel.from_pretrained("gpt2").to(device).eval()
    vocab_size = tokenizer.vocab_size
    mask_token_id = tokenizer.eos_token_id  # use EOS as mask

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    all_sigma_profiles = []
    all_texts_baseline = []
    rng = np.random.default_rng(42)

    for sample_idx in range(n_samples):
        print(f"  Sample {sample_idx+1}/{n_samples}...")

        # Start from fully masked sequence
        x = torch.full((1, seq_len), mask_token_id, dtype=torch.long, device=device)

        sigma_profile = []
        prev_probs = None

        # Linear unmasking schedule: at step t, unmask (T-t)/T fraction
        for step in range(T):
            t_frac = 1.0 - step / T  # goes from 1.0 (all masked) to ~0 (all unmasked)

            # Get model predictions for masked positions
            with torch.no_grad():
                logits = model(x).logits[0]  # (seq_len, vocab_size)
                probs = torch.softmax(logits, dim=-1)  # (seq_len, vocab_size)

            # Compute Sigma (KL divergence from previous step)
            if prev_probs is not None:
                p = probs.float().cpu().numpy()
                q = prev_probs
                # Per-position KL, averaged
                eps = 1e-10
                kl_per_pos = np.sum(p * np.log((p + eps) / (q + eps)), axis=-1)
                sigma_t = float(np.mean(kl_per_pos))
            else:
                sigma_t = 0.0

            prev_probs = probs.float().cpu().numpy()

            cumulative = sum(s["sigma_t"] for s in sigma_profile) + sigma_t
            petz_bound = float(np.exp(-cumulative / 2.0))

            sigma_profile.append({
                "step": step,
                "sigma_t": sigma_t,
                "cumulative_sigma": cumulative,
                "petz_fidelity_bound": petz_bound,
                "masked_fraction": float(t_frac),
            })

            # Unmask some positions
            n_to_unmask = max(1, int(seq_len / T))
            masked_positions = (x[0] == mask_token_id).nonzero(as_tuple=True)[0]

            if len(masked_positions) > 0:
                n_unmask = min(n_to_unmask, len(masked_positions))
                # Pick positions with highest confidence
                mask_probs = probs[masked_positions]
                confidence = mask_probs.max(dim=-1).values
                _, top_idx = confidence.topk(min(n_unmask, len(confidence)))
                positions_to_unmask = masked_positions[top_idx]

                # Sample tokens for those positions
                for pos in positions_to_unmask:
                    sampled = torch.multinomial(probs[pos], 1)
                    x[0, pos] = sampled

        # Decode final text
        text = tokenizer.decode(x[0].cpu().tolist(), skip_special_tokens=True)
        all_texts_baseline.append(text)
        all_sigma_profiles.append(sigma_profile)

        print(f"    Text: {text[:80]}...")
        print(f"    Total Σ: {sigma_profile[-1]['cumulative_sigma']:.4f}")
        print(f"    Petz F: {sigma_profile[-1]['petz_fidelity_bound']:.4f}")

    # Save results
    with open(RESULTS_DIR / "baseline_texts.json", "w") as f:
        json.dump(all_texts_baseline, f, indent=2, cls=NumpyEncoder)

    with open(RESULTS_DIR / "baseline_sigma_profiles.json", "w") as f:
        json.dump(all_sigma_profiles, f, indent=2, cls=NumpyEncoder)

    # Compute average sigma profile
    avg_sigma = np.mean([
        [s["sigma_t"] for s in prof] for prof in all_sigma_profiles
    ], axis=0)

    np.save(RESULTS_DIR / "avg_sigma_profile.npy", avg_sigma)

    print(f"\n  Baseline complete. Avg total Σ: {np.sum(avg_sigma):.4f}")
    return avg_sigma, all_sigma_profiles, all_texts_baseline


def run_sigma_guided(
    avg_sigma: np.ndarray,
    n_samples: int = 10,
    seq_len: int = 64,
    T_original: int = 200,
    target_speedup: float = 10.0,
    device: str = "mps",
):
    """Re-run with Sigma-guided step skipping."""
    import torch
    from transformers import GPT2Tokenizer, GPT2LMHeadModel

    # Determine which steps to keep
    budget = max(int(T_original / target_speedup), 10)

    # Keep steps with highest Sigma
    top_indices = np.argsort(np.abs(avg_sigma))[::-1][:budget]
    active_steps = sorted(top_indices)

    actual_speedup = T_original / len(active_steps)

    print(f"\n{'='*60}")
    print(f"  Line 2: Sigma-Guided ({len(active_steps)} steps, {actual_speedup:.1f}X)")
    print(f"{'='*60}\n")

    # Quality loss estimate
    kept_sigma = np.sum(np.abs(avg_sigma[active_steps]))
    total_sigma = np.sum(np.abs(avg_sigma))
    skipped_sigma = total_sigma - kept_sigma
    quality_loss = 1.0 - np.exp(-skipped_sigma / 2.0)
    print(f"  Kept Σ: {kept_sigma:.4f} / {total_sigma:.4f} ({kept_sigma/total_sigma*100:.1f}%)")
    print(f"  Estimated quality loss (Petz): {quality_loss:.4f}")

    tokenizer = GPT2Tokenizer.from_pretrained("gpt2")
    model = GPT2LMHeadModel.from_pretrained("gpt2").to(device).eval()
    mask_token_id = tokenizer.eos_token_id

    all_texts_guided = []
    all_sigma_profiles = []

    for sample_idx in range(n_samples):
        print(f"  Sample {sample_idx+1}/{n_samples}...")

        x = torch.full((1, seq_len), mask_token_id, dtype=torch.long, device=device)
        sigma_profile = []
        prev_probs = None

        for step_idx, step in enumerate(active_steps):
            t_frac = 1.0 - step / T_original

            with torch.no_grad():
                logits = model(x).logits[0]
                probs = torch.softmax(logits, dim=-1)

            if prev_probs is not None:
                p = probs.float().cpu().numpy()
                q = prev_probs
                eps = 1e-10
                kl = np.mean(np.sum(p * np.log((p + eps) / (q + eps)), axis=-1))
                sigma_t = float(kl)
            else:
                sigma_t = 0.0

            prev_probs = probs.float().cpu().numpy()

            cumulative = sum(s["sigma_t"] for s in sigma_profile) + sigma_t
            sigma_profile.append({
                "step": step,
                "sigma_t": sigma_t,
                "cumulative_sigma": cumulative,
                "petz_fidelity_bound": float(np.exp(-cumulative / 2.0)),
            })

            # Unmask more positions per step (to compensate for fewer steps)
            step_size = active_steps[step_idx + 1] - step if step_idx + 1 < len(active_steps) else T_original - step
            n_to_unmask = max(1, int(seq_len * step_size / T_original))

            masked_positions = (x[0] == mask_token_id).nonzero(as_tuple=True)[0]
            if len(masked_positions) > 0:
                n_unmask = min(n_to_unmask, len(masked_positions))
                mask_probs = probs[masked_positions]
                confidence = mask_probs.max(dim=-1).values
                _, top_idx = confidence.topk(min(n_unmask, len(confidence)))
                positions_to_unmask = masked_positions[top_idx]

                for pos in positions_to_unmask:
                    sampled = torch.multinomial(probs[pos], 1)
                    x[0, pos] = sampled

        text = tokenizer.decode(x[0].cpu().tolist(), skip_special_tokens=True)
        all_texts_guided.append(text)
        all_sigma_profiles.append(sigma_profile)
        print(f"    Text: {text[:80]}...")

    with open(RESULTS_DIR / "guided_texts.json", "w") as f:
        json.dump(all_texts_guided, f, indent=2, cls=NumpyEncoder)

    with open(RESULTS_DIR / "guided_sigma_profiles.json", "w") as f:
        json.dump(all_sigma_profiles, f, indent=2, cls=NumpyEncoder)

    schedule = {
        "original_steps": T_original,
        "guided_steps": len(active_steps),
        "speedup": actual_speedup,
        "quality_loss_estimate": quality_loss,
        "active_steps": [int(s) for s in active_steps],
    }
    with open(RESULTS_DIR / "schedule.json", "w") as f:
        json.dump(schedule, f, indent=2, cls=NumpyEncoder)

    print(f"\n  Guided complete. {actual_speedup:.1f}X speedup.")
    return all_texts_guided, schedule


def evaluate(baseline_texts, guided_texts):
    """Compare text quality between baseline and guided."""
    import torch
    from transformers import GPT2Tokenizer, GPT2LMHeadModel

    print(f"\n{'='*60}")
    print(f"  Line 2: Evaluation")
    print(f"{'='*60}\n")

    tokenizer = GPT2Tokenizer.from_pretrained("gpt2")
    model = GPT2LMHeadModel.from_pretrained("gpt2").eval()

    def compute_perplexity(texts, model, tokenizer, device="cpu"):
        total_loss = 0.0
        total_tokens = 0
        model = model.to(device)
        for text in texts:
            if not text.strip():
                continue
            inputs = tokenizer(text, return_tensors="pt", truncation=True, max_length=64).to(device)
            with torch.no_grad():
                outputs = model(**inputs, labels=inputs["input_ids"])
                total_loss += outputs.loss.item() * inputs["input_ids"].shape[1]
                total_tokens += inputs["input_ids"].shape[1]
        return float(np.exp(total_loss / max(total_tokens, 1)))

    ppl_baseline = compute_perplexity(baseline_texts, model, tokenizer)
    ppl_guided = compute_perplexity(guided_texts, model, tokenizer)

    # Distinct-n metrics
    def distinct_n(texts, n=2):
        all_ngrams = []
        for text in texts:
            tokens = text.split()
            ngrams = [tuple(tokens[i:i+n]) for i in range(len(tokens)-n+1)]
            all_ngrams.extend(ngrams)
        if not all_ngrams:
            return 0.0
        return len(set(all_ngrams)) / len(all_ngrams)

    d2_baseline = distinct_n(baseline_texts, 2)
    d2_guided = distinct_n(guided_texts, 2)

    metrics = {
        "perplexity_baseline": ppl_baseline,
        "perplexity_guided": ppl_guided,
        "perplexity_ratio": ppl_guided / max(ppl_baseline, 1e-10),
        "distinct2_baseline": d2_baseline,
        "distinct2_guided": d2_guided,
    }

    with open(RESULTS_DIR / "metrics.json", "w") as f:
        json.dump(metrics, f, indent=2, cls=NumpyEncoder)

    print(f"  Perplexity (baseline): {ppl_baseline:.2f}")
    print(f"  Perplexity (guided):   {ppl_guided:.2f}")
    print(f"  Distinct-2 (baseline): {d2_baseline:.4f}")
    print(f"  Distinct-2 (guided):   {d2_guided:.4f}")

    # Save example texts side by side
    examples = []
    for i, (b, g) in enumerate(zip(baseline_texts[:5], guided_texts[:5])):
        examples.append({"baseline": b, "guided": g})

    with open(RESULTS_DIR / "example_comparisons.json", "w") as f:
        json.dump(examples, f, indent=2, ensure_ascii=False)

    return metrics


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--num-samples", type=int, default=5)
    parser.add_argument("--seq-len", type=int, default=64)
    parser.add_argument("--steps", type=int, default=200)
    parser.add_argument("--target-speedup", type=float, default=10.0)
    args = parser.parse_args()

    status = check_deps()
    print("Dependencies:", json.dumps(status, indent=2))

    if not status.get("transformers"):
        print("ERROR: transformers not installed")
        return

    device = "mps" if status.get("mps") else "cpu"
    print(f"Using device: {device}")

    t0 = time.time()

    # Phase 1: Baseline
    avg_sigma, profiles, baseline_texts = run_discrete_diffusion_baseline(
        n_samples=args.num_samples,
        seq_len=args.seq_len,
        T=args.steps,
        device=device,
    )

    # Phase 2: Sigma-guided
    guided_texts, schedule = run_sigma_guided(
        avg_sigma,
        n_samples=args.num_samples,
        seq_len=args.seq_len,
        T_original=args.steps,
        target_speedup=args.target_speedup,
        device=device,
    )

    # Phase 3: Evaluate
    metrics = evaluate(baseline_texts, guided_texts)

    elapsed = time.time() - t0
    summary = {
        "experiment": "Line 2: Diffusion LLM",
        "model": "GPT2 + discrete diffusion",
        "device": device,
        "baseline_steps": args.steps,
        "guided_steps": schedule["guided_steps"],
        "speedup": schedule["speedup"],
        "perplexity_baseline": metrics["perplexity_baseline"],
        "perplexity_guided": metrics["perplexity_guided"],
        "elapsed_seconds": elapsed,
    }

    with open(RESULTS_DIR / "experiment_summary.json", "w") as f:
        json.dump(summary, f, indent=2, cls=NumpyEncoder)

    print(f"\n{'='*60}")
    print(f"  SUMMARY: {schedule['speedup']:.1f}X speedup")
    print(f"  Perplexity: {metrics['perplexity_baseline']:.1f} → {metrics['perplexity_guided']:.1f}")
    print(f"  Time: {elapsed:.1f}s")
    print(f"  Results: {RESULTS_DIR}")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    main()
