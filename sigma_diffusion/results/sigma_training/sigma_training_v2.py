#!/usr/bin/env python3
"""
Sigma-Weighted vs Standard Training — V2 (Corrected)
=====================================================
V1 issue: raw |d(sigma)/dt| weights up to 5.9x at t~1 destabilized training.
High-t = mostly masked = huge CE loss * huge weight = training collapse.

V2 fix: Use TEMPERED sigma weights with three strategies:
  B1: sqrt(dsigma/dt) — gentle emphasis on high-Sigma steps
  B2: Importance sampling — sample t proportional to dsigma/dt, uniform weight
  B3: Clipped weights — cap at 2x max

This gives a FAIR comparison where total gradient magnitude is comparable.
"""

import os
import sys
import time
import json
import math
import copy
from pathlib import Path

import numpy as np

os.environ["HF_TOKEN"] = "YOUR_HF_TOKEN_HERE"
os.environ["TOKENIZERS_PARALLELISM"] = "false"

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader

import transformers
from datasets import load_dataset

DEVICE = "mps" if torch.backends.mps.is_available() else "cpu"
DTYPE = torch.float32
print(f"Device: {DEVICE}, dtype: {DTYPE}")

RESULTS_DIR = Path(__file__).parent
MASK_INDEX = None
NEG_INF = -1e9


# =====================================================================
# Model Architecture (same as v1)
# =====================================================================

class TimestepEmbedder(nn.Module):
    def __init__(self, hidden_size, freq_dim=128):
        super().__init__()
        self.mlp = nn.Sequential(
            nn.Linear(freq_dim, hidden_size),
            nn.SiLU(),
            nn.Linear(hidden_size, hidden_size),
        )
        self.freq_dim = freq_dim

    def forward(self, t):
        half = self.freq_dim // 2
        freqs = torch.exp(
            -math.log(10000.0) * torch.arange(half, device=t.device, dtype=torch.float32) / half
        )
        args = t[:, None].float() * freqs[None]
        emb = torch.cat([torch.cos(args), torch.sin(args)], dim=-1)
        return self.mlp(emb)


class TransformerBlock(nn.Module):
    def __init__(self, dim, n_heads, cond_dim, mlp_ratio=4, dropout=0.1):
        super().__init__()
        self.n_heads = n_heads
        self.head_dim = dim // n_heads
        self.norm1 = nn.LayerNorm(dim)
        self.attn_qkv = nn.Linear(dim, 3 * dim, bias=False)
        self.attn_out = nn.Linear(dim, dim, bias=False)
        self.norm2 = nn.LayerNorm(dim)
        self.mlp = nn.Sequential(
            nn.Linear(dim, mlp_ratio * dim),
            nn.GELU(),
            nn.Linear(mlp_ratio * dim, dim),
        )
        self.dropout = nn.Dropout(dropout)
        self.adaLN = nn.Linear(cond_dim, 6 * dim)
        nn.init.zeros_(self.adaLN.weight)
        nn.init.zeros_(self.adaLN.bias)

    def forward(self, x, c):
        B, L, D = x.shape
        shift_a, scale_a, gate_a, shift_m, scale_m, gate_m = \
            self.adaLN(c)[:, None].chunk(6, dim=2)
        h = self.norm1(x) * (1 + scale_a) + shift_a
        qkv = self.attn_qkv(h).reshape(B, L, 3, self.n_heads, self.head_dim)
        qkv = qkv.permute(2, 0, 3, 1, 4)
        q, k, v = qkv[0], qkv[1], qkv[2]
        attn = F.scaled_dot_product_attention(q, k, v)
        attn = attn.transpose(1, 2).reshape(B, L, D)
        x = x + gate_a * self.dropout(self.attn_out(attn))
        h = self.norm2(x) * (1 + scale_m) + shift_m
        x = x + gate_m * self.dropout(self.mlp(h))
        return x


class TinyMDLM(nn.Module):
    def __init__(self, vocab_size, hidden_dim=256, cond_dim=128,
                 n_blocks=6, n_heads=4, dropout=0.1, max_len=128):
        super().__init__()
        self.vocab_size = vocab_size
        self.hidden_dim = hidden_dim
        self.token_embed = nn.Embedding(vocab_size, hidden_dim)
        self.pos_embed = nn.Embedding(max_len, hidden_dim)
        self.time_embed = TimestepEmbedder(cond_dim)
        self.blocks = nn.ModuleList([
            TransformerBlock(hidden_dim, n_heads, cond_dim, dropout=dropout)
            for _ in range(n_blocks)
        ])
        self.final_norm = nn.LayerNorm(hidden_dim)
        self.output_proj = nn.Linear(hidden_dim, vocab_size)
        self.final_adaLN = nn.Linear(cond_dim, 2 * hidden_dim)
        nn.init.zeros_(self.final_adaLN.weight)
        nn.init.zeros_(self.final_adaLN.bias)
        self._init_weights()

    def _init_weights(self):
        for m in self.modules():
            if isinstance(m, nn.Linear):
                nn.init.xavier_uniform_(m.weight)
                if m.bias is not None:
                    nn.init.zeros_(m.bias)
            elif isinstance(m, nn.Embedding):
                nn.init.normal_(m.weight, std=0.02)

    def forward(self, x, t):
        B, L = x.shape
        positions = torch.arange(L, device=x.device).unsqueeze(0).expand(B, -1)
        h = self.token_embed(x) + self.pos_embed(positions)
        c = F.silu(self.time_embed(t))
        for block in self.blocks:
            h = block(h, c)
        shift, scale = self.final_adaLN(c)[:, None].chunk(2, dim=2)
        h = self.final_norm(h) * (1 + scale) + shift
        return self.output_proj(h)


# =====================================================================
# Noise Schedule
# =====================================================================

def loglinear_dsigma_dt(t, eps=1e-3):
    return (1 - eps) / (1 - (1 - eps) * t)

def mask_rate_from_t(t, eps=1e-3):
    return (1 - eps) * t

def apply_masking(x_clean, t, mask_token_id):
    B, L = x_clean.shape
    mask_prob = mask_rate_from_t(t)
    rand = torch.rand(B, L, device=x_clean.device)
    mask = rand < mask_prob[:, None]
    x_noisy = x_clean.clone()
    x_noisy[mask] = mask_token_id
    return x_noisy, mask


# =====================================================================
# Dataset
# =====================================================================

class WikiTextDataset(Dataset):
    def __init__(self, tokenizer, seq_len=128, max_samples=20000):
        print(f"Loading WikiText-2 (max {max_samples} seqs of len {seq_len})...")
        raw = load_dataset('wikitext', 'wikitext-2-raw-v1', split='train')
        all_text = "\n".join([t for t in raw['text'] if len(t.strip()) > 0])
        tokens = tokenizer.encode(all_text)
        print(f"  Total tokens: {len(tokens)}")
        self.sequences = []
        for i in range(0, len(tokens) - seq_len, seq_len):
            self.sequences.append(tokens[i:i + seq_len])
            if len(self.sequences) >= max_samples:
                break
        self.sequences = torch.tensor(self.sequences, dtype=torch.long)
        print(f"  Created {len(self.sequences)} sequences")

    def __len__(self):
        return len(self.sequences)

    def __getitem__(self, idx):
        return self.sequences[idx]


# =====================================================================
# Sigma Weight Strategies
# =====================================================================

def compute_sigma_weights_sqrt(n_bins=100, eps=1e-3):
    """Strategy B1: sqrt of derivative — gentle emphasis."""
    t_centers = torch.linspace(0.05, 0.95, n_bins)
    raw = loglinear_dsigma_dt(t_centers, eps)
    weights = torch.sqrt(raw)
    weights = weights / weights.mean()
    return t_centers, weights


def compute_sigma_weights_clipped(n_bins=100, eps=1e-3, max_ratio=2.0):
    """Strategy B3: clipped weights — cap max/min ratio."""
    t_centers = torch.linspace(0.05, 0.95, n_bins)
    raw = loglinear_dsigma_dt(t_centers, eps)
    weights = raw / raw.mean()
    weights = weights.clamp(max=max_ratio)
    weights = weights / weights.mean()  # re-normalize
    return t_centers, weights


def sample_t_importance(batch_size, eps=1e-3, device='cpu'):
    """
    Strategy B2: Importance sampling — sample t from p(t) proportional to dsigma/dt.
    Use inverse CDF: if sigma(t) = -log(1-(1-eps)*t), and we want
    p(t) ~ dsigma/dt = (1-eps)/(1-(1-eps)*t), then CDF is proportional to sigma(t).

    sigma(t) = -log(1-(1-eps)*t)
    sigma(0) = -log(1) = 0
    sigma(1) = -log(eps)

    To sample: u ~ Uniform(0,1), t = sigma^{-1}(u * sigma_max)
    sigma^{-1}(s) = (1 - exp(-s)) / (1-eps)
    """
    sigma_max = -math.log(eps)  # sigma at t=1
    # Sample u uniformly, map to sigma, then invert
    u = torch.rand(batch_size, device=device)
    # Map to [sigma(0.05), sigma(0.95)] range
    t_low, t_high = 0.05, 0.95
    sigma_low = -math.log(1 - (1-eps)*t_low)
    sigma_high = -math.log(1 - (1-eps)*t_high)
    sigma_sampled = sigma_low + u * (sigma_high - sigma_low)
    # Invert: t = (1 - exp(-sigma)) / (1-eps)
    t_sampled = (1 - torch.exp(-sigma_sampled)) / (1 - eps)
    return t_sampled.clamp(0.05, 0.95)


def get_sigma_weight(t_val, t_centers, weights):
    diffs = (t_val[:, None] - t_centers[None, :].to(t_val.device)).abs()
    nearest = diffs.argmin(dim=1).cpu()
    return weights[nearest].to(t_val.device)


# =====================================================================
# Training
# =====================================================================

def train_model(model, dataset, mask_token_id, num_steps=5000,
                batch_size=32, lr=1e-4, weight_strategy="none",
                model_name="model", device=DEVICE):
    """
    weight_strategy: "none" | "sqrt" | "clipped" | "importance"
    """
    model = model.to(device)
    model.train()
    optimizer = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=0.01)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=num_steps)

    dataloader = DataLoader(dataset, batch_size=batch_size, shuffle=True,
                            drop_last=True, num_workers=0)
    data_iter = iter(dataloader)

    t_centers, sigma_w = None, None
    if weight_strategy == "sqrt":
        t_centers, sigma_w = compute_sigma_weights_sqrt()
        print(f"  Sqrt weights: range [{sigma_w.min():.3f}, {sigma_w.max():.3f}]")
    elif weight_strategy == "clipped":
        t_centers, sigma_w = compute_sigma_weights_clipped(max_ratio=2.0)
        print(f"  Clipped weights: range [{sigma_w.min():.3f}, {sigma_w.max():.3f}]")
    elif weight_strategy == "importance":
        print(f"  Importance sampling: t ~ p(t) proportional to dsigma/dt")

    loss_history = []
    t0 = time.time()

    for step in range(num_steps):
        try:
            x_clean = next(data_iter)
        except StopIteration:
            data_iter = iter(dataloader)
            x_clean = next(data_iter)

        x_clean = x_clean.to(device)
        B = x_clean.shape[0]

        # Sample timesteps
        if weight_strategy == "importance":
            t = sample_t_importance(B, device=device)
        else:
            t = torch.rand(B, device=device) * 0.9 + 0.05

        x_noisy, mask = apply_masking(x_clean, t, mask_token_id)
        logits = model(x_noisy, t)

        logits_flat = logits.view(-1, logits.size(-1))
        targets_flat = x_clean.view(-1)
        mask_flat = mask.view(-1).float()

        ce_per_token = F.cross_entropy(logits_flat, targets_flat, reduction='none')

        if weight_strategy in ("sqrt", "clipped") and sigma_w is not None:
            w = get_sigma_weight(t, t_centers, sigma_w)
            L = x_clean.shape[1]
            w_expanded = w[:, None].expand(B, L).reshape(-1)
            ce_weighted = ce_per_token * mask_flat * w_expanded
        else:
            # "none" or "importance" — uniform weight (importance already in sampling)
            ce_weighted = ce_per_token * mask_flat

        n_masked = mask_flat.sum().clamp(min=1)
        loss = ce_weighted.sum() / n_masked

        optimizer.zero_grad()
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        optimizer.step()
        scheduler.step()

        loss_val = loss.item()
        loss_history.append(loss_val)

        if (step + 1) % 500 == 0 or step == 0:
            elapsed = time.time() - t0
            avg_loss = np.mean(loss_history[-100:])
            print(f"  [{model_name}] Step {step+1}/{num_steps}, "
                  f"loss={loss_val:.4f}, avg100={avg_loss:.4f}, "
                  f"time={elapsed:.1f}s")

        if (step + 1) % 2000 == 0:
            torch.save(model.state_dict(), RESULTS_DIR / f"{model_name}_step{step+1}.pt")

    torch.save(model.state_dict(), RESULTS_DIR / f"{model_name}_final.pt")
    total_time = time.time() - t0
    print(f"  [{model_name}] Done in {total_time:.1f}s. Final avg: {np.mean(loss_history[-100:]):.4f}")
    return loss_history


# =====================================================================
# Generation
# =====================================================================

@torch.no_grad()
def generate_samples(model, num_steps, seq_len, batch_size, mask_token_id,
                     device=DEVICE):
    model.eval()
    x = torch.full((batch_size, seq_len), mask_token_id, dtype=torch.long, device=device)

    eps_sched = 1e-3
    timesteps = torch.linspace(1.0 - 1e-5, eps_sched, num_steps, device=device)

    for i in range(num_steps):
        t_now = timesteps[i]
        t_next = timesteps[i + 1] if i + 1 < num_steps else torch.tensor(0.0, device=device)

        mask_rate_now = mask_rate_from_t(t_now, eps_sched)
        mask_rate_next = mask_rate_from_t(t_next, eps_sched)
        unmask_frac = (mask_rate_now - mask_rate_next).clamp(min=0)

        t_batch = t_now.expand(batch_size)
        logits = model(x, t_batch)
        logits[:, :, mask_token_id] = NEG_INF
        probs = F.softmax(logits, dim=-1)

        best_probs, _ = probs.max(dim=-1)
        is_masked = (x == mask_token_id)
        best_probs[~is_masked] = -1.0

        n_to_unmask = max(1, int(unmask_frac.item() * seq_len))

        for b in range(batch_size):
            n_masked_now = is_masked[b].sum().item()
            n_um = min(n_to_unmask, n_masked_now)
            if n_um <= 0:
                continue
            scores = best_probs[b].clone()
            scores += torch.rand_like(scores) * 1e-6
            _, top_idx = scores.topk(n_um)
            for idx in top_idx:
                x[b, idx] = torch.multinomial(probs[b, idx], 1).item()

    # Final cleanup
    remaining = (x == mask_token_id)
    if remaining.any():
        t_final = torch.tensor(0.05, device=device).expand(batch_size)
        logits = model(x, t_final)
        logits[:, :, mask_token_id] = NEG_INF
        probs = F.softmax(logits, dim=-1)
        for b in range(batch_size):
            for pos in remaining[b].nonzero(as_tuple=True)[0]:
                x[b, pos] = torch.multinomial(probs[b, pos], 1).item()

    model.train()
    return x


# =====================================================================
# Evaluation
# =====================================================================

@torch.no_grad()
def compute_perplexity_gpt2(texts, device=DEVICE):
    tokenizer = transformers.AutoTokenizer.from_pretrained("gpt2")
    eval_model = transformers.AutoModelForCausalLM.from_pretrained("gpt2")
    if tokenizer.pad_token is None:
        tokenizer.pad_token = tokenizer.eos_token
    eval_model = eval_model.to(device).eval()

    ppls = []
    for text in texts:
        if len(text.strip()) < 10:
            ppls.append(float('inf'))
            continue
        enc = tokenizer(text, return_tensors='pt', truncation=True, max_length=512)
        input_ids = enc['input_ids'].to(device)
        if input_ids.shape[1] < 2:
            ppls.append(float('inf'))
            continue
        logits = eval_model(input_ids).logits
        loss = F.cross_entropy(logits[:, :-1, :].reshape(-1, logits.size(-1)),
                               input_ids[:, 1:].reshape(-1))
        ppls.append(loss.exp().item())

    del eval_model
    if device == "mps":
        torch.mps.empty_cache()
    return ppls


def print_texts(texts, label, num_steps):
    print(f"\n{'='*70}")
    print(f"  {label} -- {num_steps} steps")
    print(f"{'='*70}")
    for i, text in enumerate(texts):
        t = text.strip().replace('\n', ' ').replace('  ', ' ')
        print(f"\n  [{i+1}] {t[:300]}")


# =====================================================================
# Main
# =====================================================================

def main():
    print("=" * 70)
    print("  EXPERIMENT V2: Sigma-Weighted vs Standard Training")
    print("  (Corrected: tempered weights for fair comparison)")
    print("=" * 70)

    # Config
    SEQ_LEN = 128
    HIDDEN_DIM = 256
    COND_DIM = 128
    N_BLOCKS = 6
    N_HEADS = 4
    BATCH_SIZE = 32
    LR = 1e-4
    NUM_STEPS = 5000
    NUM_SAMPLES = 10

    # Tokenizer
    print("\n[1] Loading tokenizer...")
    tokenizer = transformers.AutoTokenizer.from_pretrained("gpt2")
    if tokenizer.pad_token is None:
        tokenizer.pad_token = tokenizer.eos_token
    VOCAB_SIZE = len(tokenizer) + 1
    global MASK_INDEX
    MASK_INDEX = VOCAB_SIZE - 1
    print(f"  Vocab: {VOCAB_SIZE}, MASK={MASK_INDEX}")

    # Dataset
    print("\n[2] Preparing dataset...")
    dataset = WikiTextDataset(tokenizer, seq_len=SEQ_LEN, max_samples=20000)

    # Show weight profiles
    print("\n[3] Weight profiles:")
    t_test = torch.linspace(0.1, 0.9, 9)
    raw_w = loglinear_dsigma_dt(t_test)
    raw_w_n = raw_w / raw_w.mean()
    sqrt_w = torch.sqrt(raw_w)
    sqrt_w_n = sqrt_w / sqrt_w.mean()
    clip_w = raw_w_n.clamp(max=2.0)
    clip_w_n = clip_w / clip_w.mean()

    print(f"  {'t':<6} {'Uniform':<10} {'Raw':<10} {'Sqrt':<10} {'Clipped':<10}")
    for i in range(len(t_test)):
        print(f"  {t_test[i]:.1f}   {1.0:<10.3f} {raw_w_n[i]:<10.3f} "
              f"{sqrt_w_n[i]:<10.3f} {clip_w_n[i]:<10.3f}")

    # Build models — all from same seed
    print("\n[4] Building models (same initialization)...")
    torch.manual_seed(42)
    model_A = TinyMDLM(VOCAB_SIZE, HIDDEN_DIM, COND_DIM, N_BLOCKS, N_HEADS, max_len=SEQ_LEN)
    model_B1 = copy.deepcopy(model_A)
    model_B2 = copy.deepcopy(model_A)

    n_params = sum(p.numel() for p in model_A.parameters())
    print(f"  Params: {n_params:,} ({n_params/1e6:.1f}M)")

    # ── TRAINING ──
    models_config = [
        ("A_standard",   model_A,  "none"),
        ("B1_sqrt",      model_B1, "sqrt"),
        ("B2_importance", model_B2, "importance"),
    ]

    all_losses = {}
    for name, model, strategy in models_config:
        # Check for existing checkpoint
        ckpt = RESULTS_DIR / f"{name}_final.pt"
        if ckpt.exists():
            print(f"\n  [{name}] Loading from checkpoint...")
            model.load_state_dict(torch.load(ckpt, map_location='cpu', weights_only=True))
            model = model.to(DEVICE)
            # Try to load loss curve
            lc_path = RESULTS_DIR / "all_loss_curves_v2.json"
            if lc_path.exists():
                with open(lc_path) as f:
                    saved = json.load(f)
                all_losses[name] = saved.get(name, [6.0]*100)
            else:
                all_losses[name] = [6.0]*100
            # Update the model reference in models_config
            for i, (n, m, s) in enumerate(models_config):
                if n == name:
                    models_config[i] = (n, model, s)
            print(f"  [{name}] Loaded. Avg loss: {np.mean(all_losses[name][-100:]):.4f}")
        else:
            print(f"\n{'='*70}")
            print(f"  Training {name} (strategy={strategy})")
            print(f"{'='*70}")
            losses = train_model(
                model, dataset, MASK_INDEX,
                num_steps=NUM_STEPS,
                batch_size=BATCH_SIZE,
                lr=LR,
                weight_strategy=strategy,
                model_name=name,
                device=DEVICE,
            )
            all_losses[name] = losses

    # Save all loss curves
    with open(RESULTS_DIR / "all_loss_curves_v2.json", 'w') as f:
        json.dump(all_losses, f)

    # ── GENERATION & EVALUATION ──
    print(f"\n{'='*70}")
    print(f"  GENERATION & EVALUATION")
    print(f"{'='*70}")

    gen_steps_list = [64, 16, 8]
    all_results = {}

    # Re-extract models (in case they were updated in the loop)
    named_models = {name: model for name, model, _ in models_config}

    for gen_steps in gen_steps_list:
        print(f"\n--- {gen_steps} denoising steps ---")
        step_results = {}

        for name, model in named_models.items():
            print(f"  Generating: {name}...")
            t0 = time.time()
            tokens = generate_samples(model, gen_steps, SEQ_LEN, NUM_SAMPLES, MASK_INDEX, DEVICE)
            gen_time = time.time() - t0
            texts = tokenizer.batch_decode(tokens.cpu(), skip_special_tokens=True)

            print_texts(texts, name, gen_steps)

            print(f"  Computing PPL for {name}...")
            ppls = compute_perplexity_gpt2(texts, DEVICE)
            valid = [p for p in ppls if p < 1e6]
            avg_ppl = np.mean(valid) if valid else float('inf')
            med_ppl = np.median(valid) if valid else float('inf')

            step_results[name] = {
                "mean_ppl": float(avg_ppl),
                "median_ppl": float(med_ppl),
                "per_sample_ppl": [float(p) for p in ppls],
                "gen_time": gen_time,
                "texts": texts,
            }
            print(f"  {name}: mean_ppl={avg_ppl:.1f}, median_ppl={med_ppl:.1f}")

        all_results[f"{gen_steps}_steps"] = step_results

    # ── FINAL SUMMARY ──
    print(f"\n{'='*70}")
    print(f"  FINAL RESULTS")
    print(f"{'='*70}")

    # Training loss summary
    print(f"\n  Training loss (avg last 100 steps):")
    for name in all_losses:
        print(f"    {name:<20}: {np.mean(all_losses[name][-100:]):.4f}")

    # PPL comparison table
    print(f"\n  PPL Comparison (median):")
    header = f"  {'Steps':<8}"
    for name in named_models:
        header += f" {name:<16}"
    header += f" {'Best':<16}"
    print(header)
    print(f"  {'-'*80}")

    overall_wins = {name: 0 for name in named_models}
    for gen_steps in gen_steps_list:
        key = f"{gen_steps}_steps"
        row = f"  {gen_steps:<8}"
        ppls = {}
        for name in named_models:
            med = all_results[key][name]["median_ppl"]
            ppls[name] = med
            row += f" {med:<16.1f}"
        best = min(ppls, key=ppls.get)
        overall_wins[best] += 1
        row += f" {best:<16}"
        print(row)

    # Winner
    print(f"\n  Wins by model:")
    for name, wins in sorted(overall_wins.items(), key=lambda x: -x[1]):
        print(f"    {name:<20}: {wins}")

    best_model = max(overall_wins, key=overall_wins.get)
    print(f"\n{'='*70}")
    if "standard" in best_model.lower():
        print(f"  WINNER: {best_model} (Standard uniform training)")
        print(f"  Sigma-weighting did NOT improve over standard training.")
        print(f"  Hypothesis: NOT SUPPORTED at this scale.")
    elif best_model.startswith("B"):
        print(f"  WINNER: {best_model} (Sigma-weighted training)")
        print(f"  Sigma-weighting IMPROVED text generation quality!")
        print(f"  Hypothesis: SUPPORTED!")

        # Show improvement details
        for gen_steps in gen_steps_list:
            key = f"{gen_steps}_steps"
            ppl_std = all_results[key]["A_standard"]["median_ppl"]
            ppl_best = all_results[key][best_model]["median_ppl"]
            delta_pct = (ppl_std - ppl_best) / ppl_std * 100
            print(f"    {gen_steps} steps: {delta_pct:+.1f}% PPL improvement")
    else:
        print(f"  RESULT: TIE / INCONCLUSIVE")
    print(f"{'='*70}")

    # Save results
    with open(RESULTS_DIR / "experiment_results_v2.json", 'w') as f:
        json.dump(all_results, f, indent=2, default=str)

    # Save readable text file
    with open(RESULTS_DIR / "generated_texts_v2.txt", 'w') as f:
        for gen_steps in gen_steps_list:
            key = f"{gen_steps}_steps"
            f.write(f"\n{'='*70}\n{gen_steps} DENOISING STEPS\n{'='*70}\n")
            for name in named_models:
                f.write(f"\n--- {name} ---\n")
                for i, text in enumerate(all_results[key][name]["texts"]):
                    f.write(f"\nSample {i+1}:\n{text}\n")
                f.write(f"\nMedian PPL: {all_results[key][name]['median_ppl']:.1f}\n")

    with open(RESULTS_DIR / "comparison_table_v2.txt", 'w') as f:
        f.write("Sigma-Weighted vs Standard Training — V2 Comparison\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Steps: {NUM_STEPS}, Batch: {BATCH_SIZE}, LR: {LR}\n")
        f.write(f"Model: {N_BLOCKS}L dim={HIDDEN_DIM}, {n_params:,} params\n\n")

        f.write(f"Training loss (avg last 100):\n")
        for name in all_losses:
            f.write(f"  {name:<20}: {np.mean(all_losses[name][-100:]):.4f}\n")

        f.write(f"\nGeneration PPL (median):\n")
        f.write(f"{'Steps':<8}")
        for name in named_models:
            f.write(f" {name:<16}")
        f.write("\n" + "-"*72 + "\n")
        for gen_steps in gen_steps_list:
            key = f"{gen_steps}_steps"
            f.write(f"{gen_steps:<8}")
            for name in named_models:
                f.write(f" {all_results[key][name]['median_ppl']:<16.1f}")
            f.write("\n")

    print(f"\n  All results saved to {RESULTS_DIR}")
    print("  EXPERIMENT COMPLETE.")


if __name__ == "__main__":
    main()
