#!/usr/bin/env python3
"""
Sigma-Weighted vs Standard Training for Discrete Diffusion LM
==============================================================
Self-contained, memory-efficient version.
Trains models ONE AT A TIME to avoid OOM.
"""

import os, sys, time, json, math
from pathlib import Path
import numpy as np

os.environ["HF_TOKEN"] = "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"
os.environ["TOKENIZERS_PARALLELISM"] = "false"

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader
import transformers
from datasets import load_dataset

DEVICE = "mps" if torch.backends.mps.is_available() else "cpu"
RESULTS_DIR = Path("/Users/akaihuangm1/Desktop/github/petz-recovery-unification/sigma_diffusion/results/sigma_training")

# ── Config ──
SEQ_LEN = 128
HIDDEN_DIM = 256
COND_DIM = 128
N_BLOCKS = 6
N_HEADS = 4
BATCH_SIZE = 32
LR = 1e-4
NUM_STEPS = 5000
VOCAB_SIZE = 50258  # GPT-2 vocab + 1 mask token
MASK_INDEX = 50257
NUM_SAMPLES = 10

print(f"Device: {DEVICE}")


# =====================================================================
# Model
# =====================================================================

class TimestepEmbedder(nn.Module):
    def __init__(self, hidden_size, freq_dim=128):
        super().__init__()
        self.mlp = nn.Sequential(
            nn.Linear(freq_dim, hidden_size), nn.SiLU(),
            nn.Linear(hidden_size, hidden_size))
        self.freq_dim = freq_dim

    def forward(self, t):
        half = self.freq_dim // 2
        freqs = torch.exp(-math.log(10000.0) *
                          torch.arange(half, device=t.device, dtype=torch.float32) / half)
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
            nn.Linear(dim, mlp_ratio * dim), nn.GELU(),
            nn.Linear(mlp_ratio * dim, dim))
        self.dropout = nn.Dropout(dropout)
        self.adaLN = nn.Linear(cond_dim, 6 * dim)
        nn.init.zeros_(self.adaLN.weight)
        nn.init.zeros_(self.adaLN.bias)

    def forward(self, x, c):
        B, L, D = x.shape
        mods = self.adaLN(c)[:, None].chunk(6, dim=2)
        shift_a, scale_a, gate_a, shift_m, scale_m, gate_m = mods
        h = self.norm1(x) * (1 + scale_a) + shift_a
        qkv = self.attn_qkv(h).reshape(B, L, 3, self.n_heads, self.head_dim)
        qkv = qkv.permute(2, 0, 3, 1, 4)
        attn = F.scaled_dot_product_attention(qkv[0], qkv[1], qkv[2])
        attn = attn.transpose(1, 2).reshape(B, L, D)
        x = x + gate_a * self.dropout(self.attn_out(attn))
        h = self.norm2(x) * (1 + scale_m) + shift_m
        x = x + gate_m * self.dropout(self.mlp(h))
        return x


class TinyMDLM(nn.Module):
    def __init__(self):
        super().__init__()
        self.token_embed = nn.Embedding(VOCAB_SIZE, HIDDEN_DIM)
        self.pos_embed = nn.Embedding(SEQ_LEN, HIDDEN_DIM)
        self.time_embed = TimestepEmbedder(COND_DIM)
        self.blocks = nn.ModuleList([
            TransformerBlock(HIDDEN_DIM, N_HEADS, COND_DIM) for _ in range(N_BLOCKS)])
        self.final_norm = nn.LayerNorm(HIDDEN_DIM)
        self.output_proj = nn.Linear(HIDDEN_DIM, VOCAB_SIZE)
        self.final_adaLN = nn.Linear(COND_DIM, 2 * HIDDEN_DIM)
        nn.init.zeros_(self.final_adaLN.weight)
        nn.init.zeros_(self.final_adaLN.bias)
        # Init
        for m in self.modules():
            if isinstance(m, nn.Linear):
                nn.init.xavier_uniform_(m.weight)
                if m.bias is not None:
                    nn.init.zeros_(m.bias)
            elif isinstance(m, nn.Embedding):
                nn.init.normal_(m.weight, std=0.02)

    def forward(self, x, t):
        B, L = x.shape
        pos = torch.arange(L, device=x.device).unsqueeze(0).expand(B, -1)
        h = self.token_embed(x) + self.pos_embed(pos)
        c = F.silu(self.time_embed(t))
        for block in self.blocks:
            h = block(h, c)
        shift, scale = self.final_adaLN(c)[:, None].chunk(2, dim=2)
        h = self.final_norm(h) * (1 + scale) + shift
        return self.output_proj(h)


# =====================================================================
# Dataset
# =====================================================================

def load_data(tokenizer):
    print("Loading WikiText-2...")
    raw = load_dataset("wikitext", "wikitext-2-raw-v1", split="train")
    all_text = "\n".join([t for t in raw["text"] if len(t.strip()) > 0])
    tokens = tokenizer.encode(all_text)
    seqs = []
    for i in range(0, len(tokens) - SEQ_LEN, SEQ_LEN):
        seqs.append(tokens[i:i + SEQ_LEN])
        if len(seqs) >= 20000:
            break
    data = torch.tensor(seqs, dtype=torch.long)
    print(f"  {len(data)} sequences of length {SEQ_LEN}")
    return data


# =====================================================================
# Noise + Masking
# =====================================================================

def apply_masking(x_clean, t):
    B, L = x_clean.shape
    mask_prob = 0.999 * t  # mask rate = (1-eps)*t
    rand = torch.rand(B, L, device=x_clean.device)
    mask = rand < mask_prob[:, None]
    x_noisy = x_clean.clone()
    x_noisy[mask] = MASK_INDEX
    return x_noisy, mask


def loglinear_dsigma_dt(t, eps=1e-3):
    return (1 - eps) / (1 - (1 - eps) * t)


def sample_t_importance(B, device):
    """Sample t proportional to dsigma/dt using inverse CDF."""
    eps = 1e-3
    sigma_low = -math.log(1 - (1 - eps) * 0.05)
    sigma_high = -math.log(1 - (1 - eps) * 0.95)
    u = torch.rand(B, device=device)
    sigma = sigma_low + u * (sigma_high - sigma_low)
    t = (1 - torch.exp(-sigma)) / (1 - eps)
    return t.clamp(0.05, 0.95)


# =====================================================================
# Training
# =====================================================================

def train_one_model(data, strategy, name, init_state_dict):
    """Train a single model with given strategy. Returns (model, loss_history)."""
    print(f"\n{'='*60}")
    print(f"  Training: {name} (strategy={strategy})")
    print(f"{'='*60}")

    model = TinyMDLM()
    model.load_state_dict(init_state_dict)
    model = model.to(DEVICE)
    model.train()

    optimizer = torch.optim.AdamW(model.parameters(), lr=LR, weight_decay=0.01)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=NUM_STEPS)

    # Precompute sigma weights for sqrt strategy
    sigma_w = None
    t_bins = None
    if strategy == "sqrt":
        t_bins = torch.linspace(0.05, 0.95, 100)
        raw_w = loglinear_dsigma_dt(t_bins)
        sigma_w = torch.sqrt(raw_w)
        sigma_w = sigma_w / sigma_w.mean()
        print(f"  Sqrt weights: [{sigma_w.min():.3f}, {sigma_w.max():.3f}]")

    # DataLoader
    dataset_t = torch.utils.data.TensorDataset(data)
    loader = DataLoader(dataset_t, batch_size=BATCH_SIZE, shuffle=True, drop_last=True)
    data_iter = iter(loader)

    loss_history = []
    t0 = time.time()

    for step in range(NUM_STEPS):
        try:
            (x_clean,) = next(data_iter)
        except StopIteration:
            data_iter = iter(loader)
            (x_clean,) = next(data_iter)

        x_clean = x_clean.to(DEVICE)
        B = x_clean.shape[0]

        # Sample timesteps
        if strategy == "importance":
            t = sample_t_importance(B, DEVICE)
        else:
            t = torch.rand(B, device=DEVICE) * 0.9 + 0.05

        x_noisy, mask = apply_masking(x_clean, t)
        logits = model(x_noisy, t)

        ce = F.cross_entropy(logits.view(-1, VOCAB_SIZE), x_clean.view(-1), reduction="none")
        mask_flat = mask.view(-1).float()

        if strategy == "sqrt" and sigma_w is not None:
            diffs = (t[:, None] - t_bins[None, :].to(DEVICE)).abs()
            nearest = diffs.argmin(dim=1).cpu()
            w = sigma_w[nearest].to(DEVICE)
            w_exp = w[:, None].expand(B, SEQ_LEN).reshape(-1)
            ce_weighted = ce * mask_flat * w_exp
        else:
            ce_weighted = ce * mask_flat

        loss = ce_weighted.sum() / mask_flat.sum().clamp(min=1)

        optimizer.zero_grad()
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        optimizer.step()
        scheduler.step()

        loss_history.append(loss.item())

        if (step + 1) % 500 == 0 or step == 0:
            elapsed = time.time() - t0
            avg = np.mean(loss_history[-100:])
            print(f"  [{name}] Step {step+1}/{NUM_STEPS}, "
                  f"loss={loss.item():.4f}, avg100={avg:.4f}, "
                  f"time={elapsed:.0f}s")

    torch.save(model.state_dict(), RESULTS_DIR / f"{name}_final.pt")
    total = time.time() - t0
    final_avg = np.mean(loss_history[-100:])
    print(f"  [{name}] Done in {total:.0f}s. Final avg loss: {final_avg:.4f}")

    return model, loss_history


# =====================================================================
# Generation
# =====================================================================

@torch.no_grad()
def generate(model, num_steps, batch_size=NUM_SAMPLES):
    model.eval()
    x = torch.full((batch_size, SEQ_LEN), MASK_INDEX, dtype=torch.long, device=DEVICE)
    timesteps = torch.linspace(1.0 - 1e-5, 1e-3, num_steps, device=DEVICE)

    for i in range(num_steps):
        t_now = timesteps[i]
        t_next = timesteps[i + 1] if i + 1 < num_steps else torch.tensor(0.0, device=DEVICE)
        unmask_frac = (0.999 * t_now - 0.999 * t_next).clamp(min=0)

        logits = model(x, t_now.expand(batch_size))
        logits[:, :, MASK_INDEX] = -1e9
        probs = F.softmax(logits, dim=-1)

        best_probs, _ = probs.max(dim=-1)
        is_masked = (x == MASK_INDEX)
        best_probs[~is_masked] = -1.0

        n_to_unmask = max(1, int(unmask_frac.item() * SEQ_LEN))

        for b in range(batch_size):
            n_m = is_masked[b].sum().item()
            n_um = min(n_to_unmask, n_m)
            if n_um <= 0:
                continue
            scores = best_probs[b] + torch.rand_like(best_probs[b]) * 1e-6
            _, top_idx = scores.topk(n_um)
            for idx in top_idx:
                x[b, idx] = torch.multinomial(probs[b, idx], 1).item()

    # Final cleanup
    remaining = (x == MASK_INDEX)
    if remaining.any():
        logits = model(x, torch.tensor(0.05, device=DEVICE).expand(batch_size))
        logits[:, :, MASK_INDEX] = -1e9
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
def compute_ppl(texts):
    tokenizer = transformers.AutoTokenizer.from_pretrained("gpt2")
    judge = transformers.AutoModelForCausalLM.from_pretrained("gpt2").to(DEVICE).eval()
    if tokenizer.pad_token is None:
        tokenizer.pad_token = tokenizer.eos_token

    ppls = []
    for text in texts:
        if len(text.strip()) < 10:
            ppls.append(float("inf"))
            continue
        enc = tokenizer(text, return_tensors="pt", truncation=True, max_length=512)
        ids = enc["input_ids"].to(DEVICE)
        if ids.shape[1] < 2:
            ppls.append(float("inf"))
            continue
        logits = judge(ids).logits
        loss = F.cross_entropy(logits[:, :-1].reshape(-1, logits.size(-1)),
                               ids[:, 1:].reshape(-1))
        ppls.append(loss.exp().item())

    del judge
    torch.mps.empty_cache()
    return ppls


# =====================================================================
# Main
# =====================================================================

def main():
    print("=" * 60)
    print("  Sigma-Weighted vs Standard Training Experiment")
    print("=" * 60)

    tokenizer = transformers.AutoTokenizer.from_pretrained("gpt2")
    data = load_data(tokenizer)

    # Create initial weights (shared seed)
    torch.manual_seed(42)
    init_model = TinyMDLM()
    init_state = init_model.state_dict()
    n_params = sum(p.numel() for p in init_model.parameters())
    print(f"Model: {N_BLOCKS}L, dim={HIDDEN_DIM}, heads={N_HEADS}, params={n_params:,}")
    del init_model

    # Check for existing checkpoints
    strategies = [
        ("A_standard", "none"),
        ("B_sqrt", "sqrt"),
        ("C_importance", "importance"),
    ]

    models = {}
    all_losses = {}

    for name, strategy in strategies:
        ckpt = RESULTS_DIR / f"{name}_final.pt"
        if ckpt.exists():
            print(f"\n  {name}: loading from checkpoint")
            m = TinyMDLM()
            m.load_state_dict(torch.load(ckpt, map_location="cpu", weights_only=True))
            m = m.to(DEVICE)
            models[name] = m
            all_losses[name] = []  # No loss history available
        else:
            m, losses = train_one_model(data, strategy, name, init_state)
            models[name] = m
            all_losses[name] = losses

        # Free GPU memory between models
        torch.mps.empty_cache()

    # Save loss curves
    with open(RESULTS_DIR / "loss_curves_final.json", "w") as f:
        json.dump(all_losses, f)

    # ── GENERATION & EVALUATION ──
    print(f"\n{'='*60}")
    print(f"  GENERATION & EVALUATION")
    print(f"{'='*60}")

    gen_steps_list = [64, 16, 8]
    all_results = {}

    for gs in gen_steps_list:
        print(f"\n--- {gs} denoising steps ---")
        step_results = {}

        for name, model in models.items():
            print(f"\n  Generating: {name} ({gs} steps)...")
            t0 = time.time()
            tokens = generate(model, gs)
            gen_time = time.time() - t0
            texts = tokenizer.batch_decode(tokens.cpu(), skip_special_tokens=True)

            # Print all texts
            print(f"\n  === {name} -- {gs} steps ===")
            for i, text in enumerate(texts):
                t_clean = text.strip().replace("\n", " ").replace("  ", " ")
                print(f"  [{i+1}] {t_clean[:250]}")

            # Perplexity
            print(f"  Computing PPL...")
            ppls = compute_ppl(texts)
            valid = [p for p in ppls if p < 1e6]
            mean_ppl = np.mean(valid) if valid else float("inf")
            med_ppl = np.median(valid) if valid else float("inf")
            print(f"  {name}: mean_ppl={mean_ppl:.1f}, median_ppl={med_ppl:.1f}, time={gen_time:.1f}s")

            step_results[name] = {
                "mean_ppl": float(mean_ppl),
                "median_ppl": float(med_ppl),
                "ppls": [float(p) for p in ppls],
                "gen_time": gen_time,
                "texts": texts,
            }

        all_results[f"{gs}_steps"] = step_results

    # ── FINAL SUMMARY ──
    print(f"\n{'='*60}")
    print(f"  FINAL RESULTS")
    print(f"{'='*60}")

    # Training losses
    print(f"\n  Training loss (avg last 100):")
    for name in all_losses:
        if all_losses[name]:
            print(f"    {name:<20}: {np.mean(all_losses[name][-100:]):.4f}")
        else:
            print(f"    {name:<20}: (loaded from checkpoint)")

    # PPL table
    model_names = list(models.keys())
    print(f"\n  {'Steps':<8}", end="")
    for name in model_names:
        print(f" {name:<18}", end="")
    print(f" {'Best':<18}")
    print(f"  {'-'*76}")

    wins = {n: 0 for n in model_names}
    for gs in gen_steps_list:
        key = f"{gs}_steps"
        ppls = {}
        row = f"  {gs:<8}"
        for name in model_names:
            med = all_results[key][name]["median_ppl"]
            ppls[name] = med
            row += f" {med:<18.1f}"
        best = min(ppls, key=ppls.get)
        wins[best] += 1
        row += f" {best:<18}"
        print(row)

    # Winner declaration
    print(f"\n  Wins: {wins}")
    best_model = max(wins, key=wins.get)

    print(f"\n{'='*60}")
    if best_model == "A_standard":
        print(f"  WINNER: A_standard (Standard uniform training)")
        print(f"  Sigma-weighting did NOT help at this scale.")
    else:
        print(f"  WINNER: {best_model} (Sigma-weighted training)")
        print(f"  Sigma-weighting IMPROVED text generation!")
        for gs in gen_steps_list:
            key = f"{gs}_steps"
            p_std = all_results[key]["A_standard"]["median_ppl"]
            p_best = all_results[key][best_model]["median_ppl"]
            pct = (p_std - p_best) / p_std * 100
            print(f"    {gs} steps: {pct:+.1f}% PPL improvement")
    print(f"{'='*60}")

    # Save
    with open(RESULTS_DIR / "final_results.json", "w") as f:
        json.dump(all_results, f, indent=2, default=str)

    with open(RESULTS_DIR / "final_generated_texts.txt", "w") as f:
        for gs in gen_steps_list:
            key = f"{gs}_steps"
            f.write(f"\n{'='*60}\n  {gs} DENOISING STEPS\n{'='*60}\n")
            for name in model_names:
                f.write(f"\n--- {name} ---\n")
                f.write(f"Median PPL: {all_results[key][name]['median_ppl']:.1f}\n")
                for i, text in enumerate(all_results[key][name]["texts"]):
                    f.write(f"\nSample {i+1}:\n{text}\n")

    with open(RESULTS_DIR / "final_comparison.txt", "w") as f:
        f.write("Sigma-Weighted vs Standard Training\n")
        f.write("=" * 60 + "\n")
        f.write(f"Steps: {NUM_STEPS}, Batch: {BATCH_SIZE}, LR: {LR}\n")
        f.write(f"Model: {N_BLOCKS}L dim={HIDDEN_DIM}, {n_params:,} params\n\n")
        f.write(f"{'Steps':<8}")
        for n in model_names:
            f.write(f" {n:<18}")
        f.write("\n" + "-"*74 + "\n")
        for gs in gen_steps_list:
            key = f"{gs}_steps"
            f.write(f"{gs:<8}")
            for n in model_names:
                f.write(f" {all_results[key][n]['median_ppl']:<18.1f}")
            f.write("\n")
        f.write(f"\nWINNER: {best_model}\n")

    print(f"\n  Results saved to {RESULTS_DIR}")
    print("  DONE.")


if __name__ == "__main__":
    main()
