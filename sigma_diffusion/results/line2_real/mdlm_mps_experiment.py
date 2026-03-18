#!/usr/bin/env python3
"""
MDLM (Masked Diffusion Language Model) on Apple MPS -- Sigma-Guided Experiment
===============================================================================
Loads the pretrained kuleshov-group/mdlm-owt checkpoint (NeurIPS 2024),
replaces flash_attn with standard PyTorch SDPA so it runs on Apple Silicon,
then compares:
  1. Baseline (full 256 steps)
  2. Uniform step reduction (128 steps)
  3. Sigma-guided step selection (128 steps, keeping only high-Sigma steps)

Measures: generation time, perplexity (GPT-2 as judge), Sigma profile.
"""

import os
import sys
import time
import json
import math
import typing
from pathlib import Path
from dataclasses import dataclass, field

import numpy as np

# ── HF auth ──────────────────────────────────────────────────────────
os.environ["HF_TOKEN"] = "YOUR_HF_TOKEN_HERE"
os.environ["TOKENIZERS_PARALLELISM"] = "false"

import torch
import torch.nn as nn
import torch.nn.functional as F
from einops import rearrange
from safetensors import safe_open
from huggingface_hub import hf_hub_download
import transformers

# =====================================================================
# 1.  MPS-compatible MDLM model  (no flash_attn)
# =====================================================================

DEVICE = "mps" if torch.backends.mps.is_available() else "cpu"
DTYPE = torch.float32   # MPS doesn't support bf16 well; use fp32
print(f"Device: {DEVICE}, dtype: {DTYPE}")


class LayerNorm(nn.Module):
    def __init__(self, dim):
        super().__init__()
        self.weight = nn.Parameter(torch.ones([dim]))
        self.dim = dim

    def forward(self, x):
        x = F.layer_norm(x.float(), [self.dim])
        return x * self.weight[None, None, :]


class TimestepEmbedder(nn.Module):
    def __init__(self, hidden_size, frequency_embedding_size=256):
        super().__init__()
        self.mlp = nn.Sequential(
            nn.Linear(frequency_embedding_size, hidden_size, bias=True),
            nn.SiLU(),
            nn.Linear(hidden_size, hidden_size, bias=True))
        self.frequency_embedding_size = frequency_embedding_size

    @staticmethod
    def timestep_embedding(t, dim, max_period=10000):
        half = dim // 2
        freqs = torch.exp(
            -math.log(max_period)
            * torch.arange(start=0, end=half, dtype=torch.float32)
            / half).to(device=t.device)
        args = t[:, None].float() * freqs[None]
        embedding = torch.cat([torch.cos(args), torch.sin(args)], dim=-1)
        if dim % 2:
            embedding = torch.cat(
                [embedding, torch.zeros_like(embedding[:, :1])], dim=-1)
        return embedding

    def forward(self, t):
        t_freq = self.timestep_embedding(t, self.frequency_embedding_size)
        t_emb = self.mlp(t_freq)
        return t_emb


class EmbeddingLayer(nn.Module):
    def __init__(self, dim, vocab_dim):
        super().__init__()
        self.embedding = nn.Parameter(torch.empty((vocab_dim, dim)))
        torch.nn.init.kaiming_uniform_(self.embedding, a=math.sqrt(5))

    def forward(self, x):
        return self.embedding[x]


class Rotary(nn.Module):
    def __init__(self, dim, base=10_000):
        super().__init__()
        inv_freq = 1.0 / (base ** (torch.arange(0, dim, 2).float() / dim))
        self.register_buffer('inv_freq', inv_freq)
        self.seq_len_cached = None
        self.cos_cached = None
        self.sin_cached = None

    def forward(self, x, seq_dim=1):
        seq_len = x.shape[seq_dim]
        if seq_len != self.seq_len_cached:
            self.seq_len_cached = seq_len
            t = torch.arange(seq_len, device=x.device).type_as(self.inv_freq)
            freqs = torch.einsum("i,j->ij", t, self.inv_freq.clone())
            emb = torch.cat((freqs, freqs), dim=-1).to(x.device)
            self.cos_cached = emb.cos()
            self.sin_cached = emb.sin()
        return self.cos_cached, self.sin_cached


def apply_rotary_pos_emb_manual(q, k, cos, sin):
    """Apply RoPE to q and k.
    q, k: (batch, heads, seq, dim)
    cos, sin: (seq, dim)
    """
    cos = cos[None, None, :, :]  # (1, 1, seq, dim)
    sin = sin[None, None, :, :]

    def rotate_half(x):
        x1, x2 = x[..., :x.shape[-1]//2], x[..., x.shape[-1]//2:]
        return torch.cat((-x2, x1), dim=-1)

    q_rot = q * cos + rotate_half(q) * sin
    k_rot = k * cos + rotate_half(k) * sin
    return q_rot, k_rot


class DDiTBlock(nn.Module):
    """DiT block with standard PyTorch attention (no flash_attn)."""
    def __init__(self, dim, n_heads, cond_dim, mlp_ratio=4, dropout=0.1):
        super().__init__()
        self.n_heads = n_heads
        self.head_dim = dim // n_heads

        self.norm1 = LayerNorm(dim)
        self.attn_qkv = nn.Linear(dim, 3 * dim, bias=False)
        self.attn_out = nn.Linear(dim, dim, bias=False)
        self.dropout1 = nn.Dropout(dropout)

        self.norm2 = LayerNorm(dim)
        self.mlp = nn.Sequential(
            nn.Linear(dim, mlp_ratio * dim, bias=True),
            nn.GELU(approximate='tanh'),
            nn.Linear(mlp_ratio * dim, dim, bias=True))
        self.dropout2 = nn.Dropout(dropout)
        self.dropout = dropout

        self.adaLN_modulation = nn.Linear(cond_dim, 6 * dim, bias=True)
        self.adaLN_modulation.weight.data.zero_()
        self.adaLN_modulation.bias.data.zero_()

    def forward(self, x, rotary_cos_sin, c):
        batch_size, seq_len, dim = x.shape

        # AdaLN modulation
        (shift_msa, scale_msa, gate_msa,
         shift_mlp, scale_mlp, gate_mlp) = self.adaLN_modulation(c)[:, None].chunk(6, dim=2)

        # Self-attention
        x_skip = x
        x_norm = self.norm1(x) * (1 + scale_msa) + shift_msa

        qkv = self.attn_qkv(x_norm)
        qkv = rearrange(qkv, 'b s (three h d) -> three b h s d',
                         three=3, h=self.n_heads)
        q, k, v = qkv[0], qkv[1], qkv[2]

        # Apply rotary embeddings
        cos, sin = rotary_cos_sin
        q, k = apply_rotary_pos_emb_manual(q, k, cos, sin)

        # Standard scaled dot-product attention (MPS-compatible)
        attn_out = F.scaled_dot_product_attention(q, k, v, dropout_p=0.0)
        attn_out = rearrange(attn_out, 'b h s d -> b s (h d)')

        x = x_skip + gate_msa * F.dropout(
            self.attn_out(attn_out), p=self.dropout, training=self.training)

        # MLP
        x_norm2 = self.norm2(x) * (1 + scale_mlp) + shift_mlp
        x = x + gate_mlp * F.dropout(
            self.mlp(x_norm2), p=self.dropout, training=self.training)

        return x


class DDitFinalLayer(nn.Module):
    def __init__(self, hidden_size, out_channels, cond_dim):
        super().__init__()
        self.norm_final = LayerNorm(hidden_size)
        self.linear = nn.Linear(hidden_size, out_channels)
        self.linear.weight.data.zero_()
        self.linear.bias.data.zero_()
        self.adaLN_modulation = nn.Linear(cond_dim, 2 * hidden_size, bias=True)
        self.adaLN_modulation.weight.data.zero_()
        self.adaLN_modulation.bias.data.zero_()

    def forward(self, x, c):
        shift, scale = self.adaLN_modulation(c)[:, None].chunk(2, dim=2)
        x = self.norm_final(x) * (1 + scale) + shift
        x = self.linear(x)
        return x


class MDLMModel(nn.Module):
    """MPS-compatible MDLM backbone. Same architecture, no flash_attn."""
    def __init__(self, vocab_size=50258, hidden_dim=768, cond_dim=128,
                 n_blocks=12, n_heads=12, dropout=0.1):
        super().__init__()
        self.vocab_size = vocab_size
        self.vocab_embed = EmbeddingLayer(hidden_dim, vocab_size)
        self.sigma_map = TimestepEmbedder(cond_dim)
        self.rotary_emb = Rotary(hidden_dim // n_heads)

        self.blocks = nn.ModuleList([
            DDiTBlock(hidden_dim, n_heads, cond_dim, dropout=dropout)
            for _ in range(n_blocks)
        ])
        self.output_layer = DDitFinalLayer(hidden_dim, vocab_size, cond_dim)

    def forward(self, indices, sigma):
        """
        indices: (B, L) int64 token ids
        sigma: (B,) float32 noise level
        Returns: (B, L, V) logits
        """
        x = self.vocab_embed(indices)
        c = F.silu(self.sigma_map(sigma))
        cos, sin = self.rotary_emb(x)

        for block in self.blocks:
            x = block(x, (cos, sin), c)
        logits = self.output_layer(x, c)
        return logits


def load_mdlm_weights(model: MDLMModel):
    """Load pretrained weights from kuleshov-group/mdlm-owt."""
    print("Downloading MDLM-OWT weights from HuggingFace...")
    path = hf_hub_download('kuleshov-group/mdlm-owt', 'model.safetensors')
    print(f"  Weights file: {path}")

    state_dict = {}
    with safe_open(path, framework='pt', device='cpu') as f:
        for key in f.keys():
            state_dict[key] = f.get_tensor(key)

    # The HF checkpoint has keys like "backbone.blocks.0.norm1.weight"
    # Our model has keys like "blocks.0.norm1.weight"
    # Strip the "backbone." prefix
    cleaned = {}
    for k, v in state_dict.items():
        new_key = k.replace("backbone.", "")
        cleaned[new_key] = v

    missing, unexpected = model.load_state_dict(cleaned, strict=False)
    print(f"  Loaded weights: {len(cleaned)} tensors")
    if missing:
        print(f"  Missing keys: {missing}")
    if unexpected:
        print(f"  Unexpected keys: {unexpected}")
    return model


# =====================================================================
# 2.  Diffusion sampling logic (from diffusion.py, adapted for MPS)
# =====================================================================

MASK_INDEX = 50257   # vocab_size = 50258, mask is last token
NEG_INF = -1000000.0


def subs_parameterization(logits, xt, mask_index=MASK_INDEX):
    """Apply SUBS parameterization to logits."""
    logits[:, :, mask_index] += NEG_INF
    logits = logits - torch.logsumexp(logits, dim=-1, keepdim=True)
    unmasked = (xt != mask_index)
    logits[unmasked] = NEG_INF
    logits[unmasked, xt[unmasked]] = 0
    return logits


def sample_categorical(categorical_probs):
    gumbel_norm = 1e-10 - (torch.rand_like(categorical_probs) + 1e-10).log()
    return (categorical_probs / gumbel_norm).argmax(dim=-1)


def loglinear_noise(t, eps=1e-3):
    """LogLinear noise schedule: sigma(t) = -log(1 - (1-eps)*t)."""
    sigma = -torch.log1p(-(1 - eps) * t)
    dsigma = (1 - eps) / (1 - (1 - eps) * t)
    return sigma, dsigma


@torch.no_grad()
def ddpm_cache_update(model, x, t, dt, p_x0=None):
    """One DDPM caching update step. Returns (p_x0_cache, x_next)."""
    sigma, _ = loglinear_noise(t)
    if t.ndim > 1:
        t_flat = t.squeeze(-1)
    else:
        t_flat = t

    move_chance_t = t_flat[:, None, None]
    move_chance_s = (t_flat - dt)[:, None, None]

    if p_x0 is None:
        logits = model(x, sigma.squeeze(-1) if sigma.ndim > 1 else sigma)
        log_p_x0 = subs_parameterization(logits, x)
        p_x0 = log_p_x0.exp()

    q_xs = p_x0 * (move_chance_t - move_chance_s)
    q_xs[:, :, MASK_INDEX] = move_chance_s[:, :, 0]
    _x = sample_categorical(q_xs)

    copy_flag = (x != MASK_INDEX).to(x.dtype)
    x_next = copy_flag * x + (1 - copy_flag) * _x
    return p_x0, x_next


@torch.no_grad()
def compute_step_sigma(model, x, t, dt):
    """Compute Sigma (KL divergence) for a single diffusion step.

    Sigma = KL(p(x0|xt, t) || p(x0|xt, t-dt))
    Approximated by comparing the model's output distributions at t vs t-dt.
    """
    sigma_t, _ = loglinear_noise(t)
    sigma_s, _ = loglinear_noise(t - dt)

    logits_t = model(x, sigma_t.squeeze(-1) if sigma_t.ndim > 1 else sigma_t)
    log_p_t = subs_parameterization(logits_t.clone(), x)

    logits_s = model(x, sigma_s.squeeze(-1) if sigma_s.ndim > 1 else sigma_s)
    log_p_s = subs_parameterization(logits_s.clone(), x)

    # Only compute KL over masked positions
    masked = (x == MASK_INDEX)

    p_t = log_p_t.exp()

    # KL(p_t || p_s) = sum p_t * (log_p_t - log_p_s)
    kl = (p_t * (log_p_t - log_p_s)).sum(dim=-1)  # (B, L)
    kl = kl * masked.float()

    # Average over batch and masked positions
    n_masked = masked.float().sum(dim=-1).clamp(min=1)
    kl_per_sample = kl.sum(dim=-1) / n_masked  # (B,)
    return kl_per_sample.mean().item()


@torch.no_grad()
def generate_baseline(model, num_steps, seq_len=1024, batch_size=1,
                      device=DEVICE, record_sigma=False):
    """Standard MDLM generation with all steps."""
    x = MASK_INDEX * torch.ones(batch_size, seq_len, dtype=torch.int64, device=device)
    eps = 1e-5
    timesteps = torch.linspace(1, eps, num_steps + 1, device=device)
    dt = (1 - eps) / num_steps
    p_x0_cache = None
    sigma_values = []

    t0 = time.time()
    for i in range(num_steps):
        t = timesteps[i] * torch.ones(batch_size, 1, device=device)

        if record_sigma and i % max(1, num_steps // 64) == 0:
            # Record Sigma at sampled timesteps
            sigma_val = compute_step_sigma(model, x, t, dt)
            sigma_values.append((i, timesteps[i].item(), sigma_val))

        p_x0_cache, x_next = ddpm_cache_update(model, x, t, dt, p_x0=p_x0_cache)
        if not torch.allclose(x_next, x):
            p_x0_cache = None
        x = x_next

    # Final denoising
    t = timesteps[-1] * torch.ones(batch_size, 1, device=device)
    sigma_final, _ = loglinear_noise(t)
    logits = model(x, sigma_final.squeeze(-1) if sigma_final.ndim > 1 else sigma_final)
    x = logits.argmax(dim=-1)

    elapsed = time.time() - t0
    return x, elapsed, sigma_values


@torch.no_grad()
def profile_sigma_all_steps(model, num_steps, seq_len=1024, batch_size=1,
                            device=DEVICE):
    """Run full generation and record Sigma at EVERY step."""
    x = MASK_INDEX * torch.ones(batch_size, seq_len, dtype=torch.int64, device=device)
    eps = 1e-5
    timesteps = torch.linspace(1, eps, num_steps + 1, device=device)
    dt = (1 - eps) / num_steps
    sigma_values = []

    print(f"  Profiling Sigma over {num_steps} steps...")
    for i in range(num_steps):
        t = timesteps[i] * torch.ones(batch_size, 1, device=device)
        sigma_val = compute_step_sigma(model, x, t, dt)
        sigma_values.append({
            'step': i,
            'timestep': timesteps[i].item(),
            'sigma': sigma_val
        })

        # Also advance the state
        p_x0, x_next = ddpm_cache_update(model, x, t, dt)
        x = x_next

        if (i + 1) % 16 == 0:
            print(f"    Step {i+1}/{num_steps}, t={timesteps[i].item():.4f}, "
                  f"Sigma={sigma_val:.6f}")

    return sigma_values


@torch.no_grad()
def generate_sigma_guided(model, num_total_steps, num_keep_steps,
                          sigma_profile, seq_len=1024, batch_size=1,
                          device=DEVICE):
    """Sigma-guided generation: only execute the highest-Sigma steps.

    sigma_profile: list of {'step': int, 'sigma': float} from profiling
    """
    # Select which steps to keep based on highest Sigma
    sorted_by_sigma = sorted(sigma_profile, key=lambda s: s['sigma'], reverse=True)
    keep_steps = set(s['step'] for s in sorted_by_sigma[:num_keep_steps])

    # Statistics
    kept_sigmas = [s['sigma'] for s in sorted_by_sigma[:num_keep_steps]]
    skipped_sigmas = [s['sigma'] for s in sorted_by_sigma[num_keep_steps:]]
    print(f"  Sigma-guided: keeping {num_keep_steps}/{num_total_steps} steps")
    print(f"    Kept Sigma range: [{min(kept_sigmas):.6f}, {max(kept_sigmas):.6f}]")
    if skipped_sigmas:
        print(f"    Skipped Sigma range: [{min(skipped_sigmas):.6f}, {max(skipped_sigmas):.6f}]")

    x = MASK_INDEX * torch.ones(batch_size, seq_len, dtype=torch.int64, device=device)
    eps = 1e-5
    timesteps = torch.linspace(1, eps, num_total_steps + 1, device=device)
    dt = (1 - eps) / num_total_steps

    t0 = time.time()
    steps_executed = 0
    for i in range(num_total_steps):
        if i not in keep_steps:
            continue
        t = timesteps[i] * torch.ones(batch_size, 1, device=device)
        _, x = ddpm_cache_update(model, x, t, dt)
        steps_executed += 1

    # Final denoising
    t = timesteps[-1] * torch.ones(batch_size, 1, device=device)
    sigma_final, _ = loglinear_noise(t)
    logits = model(x, sigma_final.squeeze(-1) if sigma_final.ndim > 1 else sigma_final)
    x = logits.argmax(dim=-1)

    elapsed = time.time() - t0
    print(f"    Steps executed: {steps_executed}")
    return x, elapsed


@torch.no_grad()
def generate_uniform_skip(model, num_total_steps, num_keep_steps,
                          seq_len=1024, batch_size=1, device=DEVICE):
    """Uniform step reduction: evenly space the kept steps."""
    x = MASK_INDEX * torch.ones(batch_size, seq_len, dtype=torch.int64, device=device)
    eps = 1e-5
    # Use fewer timesteps directly
    timesteps = torch.linspace(1, eps, num_keep_steps + 1, device=device)
    dt = (1 - eps) / num_keep_steps

    t0 = time.time()
    p_x0_cache = None
    for i in range(num_keep_steps):
        t = timesteps[i] * torch.ones(batch_size, 1, device=device)
        p_x0_cache, x_next = ddpm_cache_update(model, x, t, dt, p_x0=p_x0_cache)
        if not torch.allclose(x_next, x):
            p_x0_cache = None
        x = x_next

    # Final denoising
    t = timesteps[-1] * torch.ones(batch_size, 1, device=device)
    sigma_final, _ = loglinear_noise(t)
    logits = model(x, sigma_final.squeeze(-1) if sigma_final.ndim > 1 else sigma_final)
    x = logits.argmax(dim=-1)

    elapsed = time.time() - t0
    return x, elapsed


# =====================================================================
# 3.  Perplexity evaluation (GPT-2 as judge)
# =====================================================================

@torch.no_grad()
def compute_perplexity_gpt2(text_samples, device=DEVICE, model_name="gpt2"):
    """Compute perplexity of generated text using GPT-2."""
    print(f"  Loading {model_name} for perplexity evaluation...")
    tokenizer = transformers.AutoTokenizer.from_pretrained(model_name)
    eval_model = transformers.AutoModelForCausalLM.from_pretrained(model_name)

    if tokenizer.pad_token is None:
        tokenizer.pad_token = tokenizer.eos_token

    eval_model = eval_model.to(device).eval()

    encoded = tokenizer(
        text_samples,
        return_tensors='pt',
        truncation=True,
        padding=True,
        max_length=1024,
        return_attention_mask=True)

    input_ids = encoded['input_ids'].to(device)
    attn_mask = encoded['attention_mask'].to(device)

    logits = eval_model(input_ids, attention_mask=attn_mask).logits
    shift_logits = logits[:, :-1, :].contiguous()
    shift_labels = input_ids[:, 1:].contiguous()
    shift_mask = attn_mask[:, 1:].contiguous()

    loss = F.cross_entropy(
        shift_logits.view(-1, shift_logits.size(-1)),
        shift_labels.view(-1),
        reduction='none')
    loss = loss.view(shift_labels.shape)

    # Per-sample perplexity
    per_sample_nll = (loss * shift_mask).sum(dim=-1) / shift_mask.sum(dim=-1).clamp(min=1)
    per_sample_ppl = per_sample_nll.exp()

    avg_ppl = per_sample_ppl.mean().item()
    del eval_model
    torch.mps.empty_cache() if device == "mps" else None
    return avg_ppl, per_sample_ppl.tolist()


# =====================================================================
# 4.  Main experiment
# =====================================================================

def main():
    results_dir = Path(__file__).parent

    print("=" * 70)
    print("MDLM Sigma-Guided Diffusion Experiment on Apple MPS")
    print("=" * 70)

    # ── Load model ──
    print("\n[1/5] Building MPS-compatible MDLM model...")
    model = MDLMModel(
        vocab_size=50258,
        hidden_dim=768,
        cond_dim=128,
        n_blocks=12,
        n_heads=12,
        dropout=0.1)

    load_mdlm_weights(model)
    model = model.to(DEVICE).eval()

    n_params = sum(p.numel() for p in model.parameters())
    print(f"  Model parameters: {n_params / 1e6:.1f}M")
    print(f"  Device: {DEVICE}")

    # ── Load tokenizer ──
    print("\n[2/5] Loading GPT-2 tokenizer...")
    tokenizer = transformers.AutoTokenizer.from_pretrained("gpt2")
    if tokenizer.pad_token is None:
        tokenizer.pad_token = tokenizer.eos_token

    # ── Configuration ──
    NUM_FULL_STEPS = 256
    NUM_REDUCED_STEPS = 128
    SEQ_LEN = 256       # shorter for speed on MPS
    BATCH_SIZE = 2       # small batch for M1 memory
    NUM_PROFILE_STEPS = 128  # profile with moderate steps

    print(f"\n  Config: full={NUM_FULL_STEPS} steps, reduced={NUM_REDUCED_STEPS}, "
          f"seq_len={SEQ_LEN}, batch={BATCH_SIZE}")

    # ── Warm-up ──
    print("\n[3/5] Warm-up forward pass...")
    dummy_x = MASK_INDEX * torch.ones(1, 64, dtype=torch.int64, device=DEVICE)
    dummy_sigma = torch.tensor([0.5], device=DEVICE)
    _ = model(dummy_x, dummy_sigma)
    print("  Warm-up done.")

    # ── Profile Sigma ──
    print(f"\n[4/5] Profiling Sigma over {NUM_PROFILE_STEPS} steps...")
    sigma_profile = profile_sigma_all_steps(
        model, NUM_PROFILE_STEPS, seq_len=SEQ_LEN,
        batch_size=1, device=DEVICE)

    # Save sigma profile
    profile_path = results_dir / "sigma_profile.json"
    with open(profile_path, 'w') as f:
        json.dump(sigma_profile, f, indent=2)
    print(f"  Saved Sigma profile to {profile_path}")

    # Analyze the profile
    sigmas = [s['sigma'] for s in sigma_profile]
    print(f"\n  Sigma statistics:")
    print(f"    Min:    {min(sigmas):.6f}")
    print(f"    Max:    {max(sigmas):.6f}")
    print(f"    Mean:   {np.mean(sigmas):.6f}")
    print(f"    Std:    {np.std(sigmas):.6f}")
    print(f"    Median: {np.median(sigmas):.6f}")

    # ── Generate samples ──
    print(f"\n[5/5] Generating samples...")
    all_results = {}

    # --- A: Baseline (full steps) ---
    print(f"\n  A. Baseline ({NUM_FULL_STEPS} steps)...")
    tokens_full, time_full, _ = generate_baseline(
        model, NUM_FULL_STEPS, seq_len=SEQ_LEN,
        batch_size=BATCH_SIZE, device=DEVICE)
    text_full = tokenizer.batch_decode(tokens_full.cpu(), skip_special_tokens=True)
    for i, t in enumerate(text_full):
        print(f"    Sample {i}: {t[:200]}...")

    ppl_full, ppl_full_per = compute_perplexity_gpt2(text_full, device=DEVICE)
    print(f"    Time: {time_full:.2f}s, PPL: {ppl_full:.2f}")
    all_results['baseline'] = {
        'steps': NUM_FULL_STEPS,
        'time': time_full,
        'perplexity': ppl_full,
        'perplexity_per_sample': ppl_full_per,
        'samples': text_full
    }

    # --- B: Uniform skip (reduced steps) ---
    print(f"\n  B. Uniform skip ({NUM_REDUCED_STEPS} steps)...")
    tokens_uniform, time_uniform = generate_uniform_skip(
        model, NUM_FULL_STEPS, NUM_REDUCED_STEPS,
        seq_len=SEQ_LEN, batch_size=BATCH_SIZE, device=DEVICE)
    text_uniform = tokenizer.batch_decode(tokens_uniform.cpu(), skip_special_tokens=True)
    for i, t in enumerate(text_uniform):
        print(f"    Sample {i}: {t[:200]}...")

    ppl_uniform, ppl_uniform_per = compute_perplexity_gpt2(text_uniform, device=DEVICE)
    print(f"    Time: {time_uniform:.2f}s, PPL: {ppl_uniform:.2f}")
    all_results['uniform_skip'] = {
        'steps': NUM_REDUCED_STEPS,
        'time': time_uniform,
        'perplexity': ppl_uniform,
        'perplexity_per_sample': ppl_uniform_per,
        'samples': text_uniform
    }

    # --- C: Sigma-guided (keep high-Sigma steps) ---
    # We need to map the profile steps to the full-step schedule
    # Profile was done with NUM_PROFILE_STEPS steps; we select top NUM_REDUCED_STEPS/2
    # from those, but we need to be careful about step mapping.
    #
    # Since sigma-guided uses the SAME total step count as baseline but only executes
    # a subset, we re-profile at the baseline step count.
    print(f"\n  Profiling Sigma at {NUM_FULL_STEPS} steps for sigma-guided selection...")
    sigma_profile_full = profile_sigma_all_steps(
        model, NUM_FULL_STEPS, seq_len=SEQ_LEN,
        batch_size=1, device=DEVICE)

    print(f"\n  C. Sigma-guided ({NUM_REDUCED_STEPS} of {NUM_FULL_STEPS} steps)...")
    tokens_sigma, time_sigma = generate_sigma_guided(
        model, NUM_FULL_STEPS, NUM_REDUCED_STEPS,
        sigma_profile_full, seq_len=SEQ_LEN,
        batch_size=BATCH_SIZE, device=DEVICE)
    text_sigma = tokenizer.batch_decode(tokens_sigma.cpu(), skip_special_tokens=True)
    for i, t in enumerate(text_sigma):
        print(f"    Sample {i}: {t[:200]}...")

    ppl_sigma, ppl_sigma_per = compute_perplexity_gpt2(text_sigma, device=DEVICE)
    print(f"    Time: {time_sigma:.2f}s, PPL: {ppl_sigma:.2f}")
    all_results['sigma_guided'] = {
        'steps': NUM_REDUCED_STEPS,
        'time': time_sigma,
        'perplexity': ppl_sigma,
        'perplexity_per_sample': ppl_sigma_per,
        'samples': text_sigma
    }

    # ── Summary ──
    print("\n" + "=" * 70)
    print("RESULTS SUMMARY")
    print("=" * 70)
    print(f"{'Method':<25} {'Steps':<8} {'Time(s)':<10} {'PPL':<12} {'Speedup':<10}")
    print("-" * 65)
    for name, r in all_results.items():
        speedup = time_full / r['time'] if r['time'] > 0 else 0
        print(f"{name:<25} {r['steps']:<8} {r['time']:<10.2f} {r['perplexity']:<12.2f} {speedup:<10.2f}x")

    # ── Verdict ──
    print("\n" + "-" * 65)
    ppl_b = all_results['baseline']['perplexity']
    ppl_u = all_results['uniform_skip']['perplexity']
    ppl_s = all_results['sigma_guided']['perplexity']

    if ppl_s < ppl_u:
        print(f"VERDICT: Sigma-guided WINS over uniform skip")
        print(f"  PPL improvement: {ppl_u - ppl_s:.2f} "
              f"({(ppl_u - ppl_s)/ppl_u * 100:.1f}% better)")
    elif ppl_s > ppl_u:
        print(f"VERDICT: Sigma-guided LOSES to uniform skip")
        print(f"  PPL degradation: {ppl_s - ppl_u:.2f} "
              f"({(ppl_s - ppl_u)/ppl_u * 100:.1f}% worse)")
    else:
        print(f"VERDICT: Sigma-guided TIES with uniform skip")

    print(f"\n  Baseline PPL:      {ppl_b:.2f}")
    print(f"  Uniform skip PPL:  {ppl_u:.2f} (delta from baseline: {ppl_u - ppl_b:+.2f})")
    print(f"  Sigma-guided PPL:  {ppl_s:.2f} (delta from baseline: {ppl_s - ppl_b:+.2f})")

    # ── Save results ──
    # Remove non-serializable items
    save_results = {}
    for k, v in all_results.items():
        save_results[k] = {kk: vv for kk, vv in v.items()}

    results_path = results_dir / "experiment_results.json"
    with open(results_path, 'w') as f:
        json.dump(save_results, f, indent=2, default=str)
    print(f"\nResults saved to {results_path}")

    # Save sigma profile for the full run
    profile_full_path = results_dir / "sigma_profile_full.json"
    with open(profile_full_path, 'w') as f:
        json.dump(sigma_profile_full, f, indent=2)
    print(f"Full sigma profile saved to {profile_full_path}")


if __name__ == "__main__":
    main()
