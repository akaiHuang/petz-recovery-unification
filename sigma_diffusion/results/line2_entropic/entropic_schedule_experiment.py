#!/usr/bin/env python3
"""
Entropic Time Scheduler for MDLM Diffusion Language Model
==========================================================
Implements the key insight from NeurIPS 2025 (arXiv:2504.13612):

    Don't skip steps. REDISTRIBUTE the time schedule so each step
    contributes EQUAL information (equal Sigma / D_KL).

Approach:
1. Profile the full 128-step generation to measure D_KL at each step
2. Build the CDF of cumulative Sigma
3. Invert the CDF to create a non-uniform time schedule where 64 points
   are placed so that each interval carries equal information
4. Compare:
   A) Baseline: 128 uniform steps
   B) Uniform skip: 64 uniform steps
   C) Entropic: 64 NON-UNIFORM steps (redistributed by Sigma profile)

Device: Apple M1 Max (MPS)
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

# -- HF auth --
os.environ["HF_TOKEN"] = "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"
os.environ["TOKENIZERS_PARALLELISM"] = "false"

import torch
import torch.nn as nn
import torch.nn.functional as F
from einops import rearrange
from safetensors import safe_open
from huggingface_hub import hf_hub_download
import transformers
from transformers import GPT2TokenizerFast

# =====================================================================
# 1.  MPS-compatible MDLM model  (no flash_attn)
# =====================================================================

DEVICE = "mps" if torch.backends.mps.is_available() else "cpu"
DTYPE = torch.float32
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
    cos = cos[None, None, :, :]
    sin = sin[None, None, :, :]

    def rotate_half(x):
        x1, x2 = x[..., :x.shape[-1]//2], x[..., x.shape[-1]//2:]
        return torch.cat((-x2, x1), dim=-1)

    q_rot = q * cos + rotate_half(q) * sin
    k_rot = k * cos + rotate_half(k) * sin
    return q_rot, k_rot


class DDiTBlock(nn.Module):
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
        (shift_msa, scale_msa, gate_msa,
         shift_mlp, scale_mlp, gate_mlp) = self.adaLN_modulation(c)[:, None].chunk(6, dim=2)

        x_skip = x
        x_norm = self.norm1(x) * (1 + scale_msa) + shift_msa

        qkv = self.attn_qkv(x_norm)
        qkv = rearrange(qkv, 'b s (three h d) -> three b h s d',
                         three=3, h=self.n_heads)
        q, k, v = qkv[0], qkv[1], qkv[2]

        cos, sin = rotary_cos_sin
        q, k = apply_rotary_pos_emb_manual(q, k, cos, sin)

        attn_out = F.scaled_dot_product_attention(q, k, v, dropout_p=0.0)
        attn_out = rearrange(attn_out, 'b h s d -> b s (h d)')

        x = x_skip + gate_msa * F.dropout(
            self.attn_out(attn_out), p=self.dropout, training=self.training)

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
        x = self.vocab_embed(indices)
        c = F.silu(self.sigma_map(sigma))
        cos, sin = self.rotary_emb(x)
        for block in self.blocks:
            x = block(x, (cos, sin), c)
        logits = self.output_layer(x, c)
        return logits


def load_mdlm_weights(model: MDLMModel):
    print("Downloading MDLM-OWT weights from HuggingFace...")
    path = hf_hub_download('kuleshov-group/mdlm-owt', 'model.safetensors')
    print(f"  Weights file: {path}")

    state_dict = {}
    with safe_open(path, framework='pt', device='cpu') as f:
        for key in f.keys():
            state_dict[key] = f.get_tensor(key)

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
# 2.  Diffusion sampling logic
# =====================================================================

MASK_INDEX = 50257
NEG_INF = -1000000.0


def subs_parameterization(logits, xt, mask_index=MASK_INDEX):
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
    sigma = -torch.log1p(-(1 - eps) * t)
    dsigma = (1 - eps) / (1 - (1 - eps) * t)
    return sigma, dsigma


@torch.no_grad()
def ddpm_cache_update(model, x, t, dt, p_x0=None):
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
    """Compute D_KL for a single diffusion step."""
    sigma_t, _ = loglinear_noise(t)
    sigma_s, _ = loglinear_noise(t - dt)

    logits_t = model(x, sigma_t.squeeze(-1) if sigma_t.ndim > 1 else sigma_t)
    log_p_t = subs_parameterization(logits_t.clone(), x)

    logits_s = model(x, sigma_s.squeeze(-1) if sigma_s.ndim > 1 else sigma_s)
    log_p_s = subs_parameterization(logits_s.clone(), x)

    masked = (x == MASK_INDEX)
    p_t = log_p_t.exp()

    kl = (p_t * (log_p_t - log_p_s)).sum(dim=-1)
    kl = kl * masked.float()

    n_masked = masked.float().sum(dim=-1).clamp(min=1)
    kl_per_sample = kl.sum(dim=-1) / n_masked
    return kl_per_sample.mean().item()


# =====================================================================
# 3.  Generation functions
# =====================================================================

@torch.no_grad()
def generate_baseline(model, num_steps, seq_len=256, batch_size=1,
                      device=DEVICE):
    """Standard MDLM generation with uniform time steps."""
    x = MASK_INDEX * torch.ones(batch_size, seq_len, dtype=torch.int64, device=device)
    eps = 1e-5
    timesteps = torch.linspace(1, eps, num_steps + 1, device=device)
    dt = (1 - eps) / num_steps
    p_x0_cache = None

    t0 = time.time()
    for i in range(num_steps):
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


@torch.no_grad()
def generate_uniform_skip(model, num_keep_steps, seq_len=256, batch_size=1,
                           device=DEVICE):
    """Uniform step reduction: evenly space the kept steps in [1, eps]."""
    x = MASK_INDEX * torch.ones(batch_size, seq_len, dtype=torch.int64, device=device)
    eps = 1e-5
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


@torch.no_grad()
def generate_entropic(model, entropic_timesteps, seq_len=256, batch_size=1,
                      device=DEVICE):
    """
    Entropic Time Scheduler generation.

    entropic_timesteps: 1D tensor of NON-UNIFORM time points, sorted descending
                        from ~1.0 to ~eps. Length = num_steps + 1.
    Each consecutive pair (t_i, t_{i+1}) defines one step.
    The dt varies per step: dt_i = t_i - t_{i+1}.
    """
    num_steps = len(entropic_timesteps) - 1
    x = MASK_INDEX * torch.ones(batch_size, seq_len, dtype=torch.int64, device=device)

    t0 = time.time()
    p_x0_cache = None
    for i in range(num_steps):
        t_val = entropic_timesteps[i]
        dt_val = entropic_timesteps[i] - entropic_timesteps[i + 1]

        t = t_val * torch.ones(batch_size, 1, device=device)
        # dt must be a scalar for the update
        p_x0_cache, x_next = ddpm_cache_update(model, x, t, dt_val.item(),
                                                 p_x0=p_x0_cache)
        if not torch.allclose(x_next, x):
            p_x0_cache = None
        x = x_next

    # Final denoising at the last timestep
    t_final = entropic_timesteps[-1]
    t = t_final * torch.ones(batch_size, 1, device=device)
    sigma_final, _ = loglinear_noise(t)
    logits = model(x, sigma_final.squeeze(-1) if sigma_final.ndim > 1 else sigma_final)
    x = logits.argmax(dim=-1)

    elapsed = time.time() - t0
    return x, elapsed


@torch.no_grad()
def profile_sigma_all_steps(model, num_steps, seq_len=256, batch_size=1,
                             device=DEVICE):
    """Run full generation and record D_KL at EVERY step."""
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

        # Advance the state
        p_x0, x_next = ddpm_cache_update(model, x, t, dt)
        x = x_next

        if (i + 1) % 16 == 0:
            print(f"    Step {i+1}/{num_steps}, t={timesteps[i].item():.4f}, "
                  f"Sigma={sigma_val:.6f}")

    return sigma_values


def build_entropic_schedule(sigma_profile, num_entropic_steps, eps=1e-5):
    """
    Build the Entropic Time Schedule from a Sigma profile.

    Key idea (arXiv:2504.13612):
        Place time points so that each interval carries EQUAL cumulative
        information (equal Sigma). This is done by inverting the CDF of Sigma.

    Args:
        sigma_profile: list of {'step': int, 'timestep': float, 'sigma': float}
        num_entropic_steps: number of steps in the new schedule (e.g. 64)
        eps: minimum time value

    Returns:
        entropic_timesteps: tensor of shape (num_entropic_steps + 1,),
                            sorted descending from ~1.0 to ~eps
    """
    # Extract timesteps and sigma values (sorted by timestep, descending)
    profile_sorted = sorted(sigma_profile, key=lambda s: s['timestep'], reverse=True)
    t_values = np.array([s['timestep'] for s in profile_sorted])
    sigma_values = np.array([s['sigma'] for s in profile_sorted])

    # Ensure all sigma values are non-negative
    sigma_values = np.maximum(sigma_values, 0.0)

    # Build cumulative distribution of Sigma
    # cumulative_sigma[i] = sum of sigma from step 0 to step i
    cumulative_sigma = np.cumsum(sigma_values)
    total_sigma = cumulative_sigma[-1]

    if total_sigma < 1e-12:
        print("  WARNING: Total Sigma is near zero. Falling back to uniform schedule.")
        return torch.linspace(1, eps, num_entropic_steps + 1)

    # Normalize to [0, 1]
    cdf = cumulative_sigma / total_sigma

    # We want num_entropic_steps intervals, so num_entropic_steps + 1 points.
    # Target CDF values: evenly spaced from 0 to 1
    target_cdf = np.linspace(0, 1, num_entropic_steps + 1)

    # Invert the CDF: for each target CDF value, find the corresponding timestep
    # Prepend (0, t_values[0]) to handle CDF=0 -> t=1.0
    cdf_with_start = np.concatenate([[0.0], cdf])
    t_with_start = np.concatenate([[t_values[0]], t_values])

    # Interpolate: target_cdf -> timestep
    # np.interp requires xp to be increasing, cdf is increasing
    entropic_t = np.interp(target_cdf, cdf_with_start, t_with_start)

    # Ensure boundaries
    entropic_t[0] = max(entropic_t[0], t_values[0])  # start at ~1.0
    entropic_t[-1] = max(entropic_t[-1], eps)  # end at eps

    # Ensure strictly decreasing
    for i in range(1, len(entropic_t)):
        if entropic_t[i] >= entropic_t[i-1]:
            entropic_t[i] = entropic_t[i-1] - 1e-7

    print(f"\n  Entropic Schedule built:")
    print(f"    Total Sigma: {total_sigma:.4f}")
    print(f"    Sigma per entropic step: {total_sigma / num_entropic_steps:.6f}")
    print(f"    Time range: [{entropic_t[0]:.6f}, {entropic_t[-1]:.6f}]")
    print(f"    First 5 dt: {[f'{entropic_t[i]-entropic_t[i+1]:.6f}' for i in range(5)]}")
    print(f"    Last  5 dt: {[f'{entropic_t[i]-entropic_t[i+1]:.6f}' for i in range(len(entropic_t)-6, len(entropic_t)-1)]}")

    # Compare with uniform schedule
    uniform_dt = (t_values[0] - eps) / num_entropic_steps
    print(f"    Uniform dt would be: {uniform_dt:.6f}")

    return torch.tensor(entropic_t, dtype=torch.float32)


# =====================================================================
# 4.  Perplexity evaluation (GPT-2 as judge)
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

    per_sample_nll = (loss * shift_mask).sum(dim=-1) / shift_mask.sum(dim=-1).clamp(min=1)
    per_sample_ppl = per_sample_nll.exp()

    avg_ppl = per_sample_ppl.mean().item()
    del eval_model
    if device == "mps":
        torch.mps.empty_cache()
    return avg_ppl, per_sample_ppl.tolist()


# =====================================================================
# 5.  Main experiment
# =====================================================================

def main():
    results_dir = Path(__file__).parent

    print("=" * 70)
    print("ENTROPIC TIME SCHEDULER FOR MDLM")
    print("NeurIPS 2025 (arXiv:2504.13612) approach applied to discrete diffusion")
    print("=" * 70)

    # -- Load model --
    print("\n[1/6] Building MPS-compatible MDLM model...")
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

    # -- Load tokenizer --
    print("\n[2/6] Loading GPT-2 tokenizer...")
    tokenizer = GPT2TokenizerFast.from_pretrained('gpt2')
    if tokenizer.pad_token is None:
        tokenizer.pad_token = tokenizer.eos_token

    # -- Configuration --
    NUM_BASELINE_STEPS = 128
    NUM_REDUCED_STEPS = 64
    SEQ_LEN = 256
    NUM_SAMPLES = 3

    print(f"\n  Config:")
    print(f"    Baseline steps:   {NUM_BASELINE_STEPS}")
    print(f"    Reduced steps:    {NUM_REDUCED_STEPS}")
    print(f"    Sequence length:  {SEQ_LEN}")
    print(f"    Samples per method: {NUM_SAMPLES}")

    # -- Warm-up --
    print("\n[3/6] Warm-up forward pass...")
    dummy_x = MASK_INDEX * torch.ones(1, 64, dtype=torch.int64, device=DEVICE)
    dummy_sigma = torch.tensor([0.5], device=DEVICE)
    _ = model(dummy_x, dummy_sigma)
    print("  Warm-up done.")

    # -- Profile Sigma --
    print(f"\n[4/6] Profiling Sigma over {NUM_BASELINE_STEPS} steps...")
    sigma_profile = profile_sigma_all_steps(
        model, NUM_BASELINE_STEPS, seq_len=SEQ_LEN,
        batch_size=1, device=DEVICE)

    # Save profile
    profile_path = results_dir / "sigma_profile_128steps.json"
    with open(profile_path, 'w') as f:
        json.dump(sigma_profile, f, indent=2)
    print(f"  Saved Sigma profile to {profile_path}")

    # Analyze the profile
    sigmas = [s['sigma'] for s in sigma_profile]
    print(f"\n  Sigma statistics ({NUM_BASELINE_STEPS} steps):")
    print(f"    Min:    {min(sigmas):.6f}")
    print(f"    Max:    {max(sigmas):.6f}")
    print(f"    Mean:   {np.mean(sigmas):.6f}")
    print(f"    Std:    {np.std(sigmas):.6f}")
    print(f"    Sum:    {np.sum(sigmas):.6f}")

    # Show where Sigma is concentrated
    sigmas_arr = np.array(sigmas)
    top_quarter = sigmas_arr[:len(sigmas)//4].sum()
    second_quarter = sigmas_arr[len(sigmas)//4:len(sigmas)//2].sum()
    third_quarter = sigmas_arr[len(sigmas)//2:3*len(sigmas)//4].sum()
    bottom_quarter = sigmas_arr[3*len(sigmas)//4:].sum()
    total = sigmas_arr.sum()
    print(f"\n  Sigma distribution by quarter (t=1 -> t=0):")
    print(f"    Q1 (t~1.0 -> t~0.75): {top_quarter:.4f} ({top_quarter/total*100:.1f}%)")
    print(f"    Q2 (t~0.75 -> t~0.50): {second_quarter:.4f} ({second_quarter/total*100:.1f}%)")
    print(f"    Q3 (t~0.50 -> t~0.25): {third_quarter:.4f} ({third_quarter/total*100:.1f}%)")
    print(f"    Q4 (t~0.25 -> t~0.00): {bottom_quarter:.4f} ({bottom_quarter/total*100:.1f}%)")

    # -- Build Entropic Schedule --
    print(f"\n[5/6] Building entropic schedule ({NUM_REDUCED_STEPS} steps)...")
    entropic_timesteps = build_entropic_schedule(
        sigma_profile, NUM_REDUCED_STEPS)
    entropic_timesteps = entropic_timesteps.to(DEVICE)

    # Save entropic schedule
    schedule_path = results_dir / "entropic_schedule.json"
    with open(schedule_path, 'w') as f:
        json.dump({
            'timesteps': entropic_timesteps.cpu().tolist(),
            'num_steps': NUM_REDUCED_STEPS,
            'description': 'Non-uniform time schedule where each interval carries equal D_KL'
        }, f, indent=2)
    print(f"  Saved entropic schedule to {schedule_path}")

    # Visualize schedule density
    entropic_t_np = entropic_timesteps.cpu().numpy()
    dts = np.diff(entropic_t_np)  # negative, since descending
    print(f"\n  Entropic schedule density:")
    print(f"    Smallest dt (densest region):  {np.min(np.abs(dts)):.6f}")
    print(f"    Largest  dt (sparsest region): {np.max(np.abs(dts)):.6f}")
    print(f"    Ratio (largest/smallest):      {np.max(np.abs(dts))/np.min(np.abs(dts)):.1f}x")

    # -- Generate samples --
    print(f"\n[6/6] Generating {NUM_SAMPLES} samples per method...")
    print("=" * 70)

    all_results = {
        'baseline': {'texts': [], 'times': [], 'steps': NUM_BASELINE_STEPS},
        'uniform_skip': {'texts': [], 'times': [], 'steps': NUM_REDUCED_STEPS},
        'entropic': {'texts': [], 'times': [], 'steps': NUM_REDUCED_STEPS},
    }

    for sample_idx in range(NUM_SAMPLES):
        print(f"\n{'='*70}")
        print(f"SAMPLE {sample_idx + 1}/{NUM_SAMPLES}")
        print(f"{'='*70}")

        # Set a different random seed for each sample, but same seed across methods
        seed = 42 + sample_idx * 1000

        # --- A: Baseline (128 steps) ---
        torch.manual_seed(seed)
        if DEVICE == "mps":
            torch.mps.manual_seed(seed)
        print(f"\n  A. Baseline ({NUM_BASELINE_STEPS} steps)...")
        tokens_base, time_base = generate_baseline(
            model, NUM_BASELINE_STEPS, seq_len=SEQ_LEN,
            batch_size=1, device=DEVICE)
        text_base = tokenizer.decode(tokens_base[0].cpu(), skip_special_tokens=True)
        all_results['baseline']['texts'].append(text_base)
        all_results['baseline']['times'].append(time_base)
        print(f"  Time: {time_base:.2f}s")
        print(f"  TEXT: {text_base[:500]}")
        print(f"  [... total {len(text_base)} chars]")

        # --- B: Uniform skip (64 steps) ---
        torch.manual_seed(seed)
        if DEVICE == "mps":
            torch.mps.manual_seed(seed)
        print(f"\n  B. Uniform skip ({NUM_REDUCED_STEPS} steps)...")
        tokens_uniform, time_uniform = generate_uniform_skip(
            model, NUM_REDUCED_STEPS, seq_len=SEQ_LEN,
            batch_size=1, device=DEVICE)
        text_uniform = tokenizer.decode(tokens_uniform[0].cpu(), skip_special_tokens=True)
        all_results['uniform_skip']['texts'].append(text_uniform)
        all_results['uniform_skip']['times'].append(time_uniform)
        print(f"  Time: {time_uniform:.2f}s")
        print(f"  TEXT: {text_uniform[:500]}")
        print(f"  [... total {len(text_uniform)} chars]")

        # --- C: Entropic (64 non-uniform steps) ---
        torch.manual_seed(seed)
        if DEVICE == "mps":
            torch.mps.manual_seed(seed)
        print(f"\n  C. Entropic ({NUM_REDUCED_STEPS} non-uniform steps)...")
        tokens_entropic, time_entropic = generate_entropic(
            model, entropic_timesteps, seq_len=SEQ_LEN,
            batch_size=1, device=DEVICE)
        text_entropic = tokenizer.decode(tokens_entropic[0].cpu(), skip_special_tokens=True)
        all_results['entropic']['texts'].append(text_entropic)
        all_results['entropic']['times'].append(time_entropic)
        print(f"  Time: {time_entropic:.2f}s")
        print(f"  TEXT: {text_entropic[:500]}")
        print(f"  [... total {len(text_entropic)} chars]")

    # -- Compute perplexity --
    print(f"\n{'='*70}")
    print("PERPLEXITY EVALUATION (GPT-2 as judge)")
    print(f"{'='*70}")

    ppls = {}
    for method_name in ['baseline', 'uniform_skip', 'entropic']:
        texts = all_results[method_name]['texts']
        print(f"\n  Evaluating {method_name} ({len(texts)} samples)...")
        avg_ppl, per_sample_ppl = compute_perplexity_gpt2(texts, device=DEVICE)
        ppls[method_name] = {
            'avg': avg_ppl,
            'per_sample': per_sample_ppl
        }
        all_results[method_name]['ppl_avg'] = avg_ppl
        all_results[method_name]['ppl_per_sample'] = per_sample_ppl
        print(f"    Average PPL: {avg_ppl:.2f}")
        for i, p in enumerate(per_sample_ppl):
            print(f"    Sample {i+1} PPL: {p:.2f}")

    # -- Summary --
    print(f"\n{'='*70}")
    print("RESULTS SUMMARY")
    print(f"{'='*70}")
    print(f"{'Method':<20} {'Steps':<8} {'Avg Time':<12} {'PPL':<12} {'Speedup':<10}")
    print("-" * 62)

    baseline_time = np.mean(all_results['baseline']['times'])
    for method_name in ['baseline', 'uniform_skip', 'entropic']:
        r = all_results[method_name]
        avg_time = np.mean(r['times'])
        speedup = baseline_time / avg_time if avg_time > 0 else 0
        ppl = r.get('ppl_avg', float('nan'))
        print(f"{method_name:<20} {r['steps']:<8} {avg_time:<12.2f} {ppl:<12.2f} {speedup:<10.2f}x")

    # -- Verdict --
    print(f"\n{'='*70}")
    print("VERDICT")
    print(f"{'='*70}")

    ppl_b = ppls['baseline']['avg']
    ppl_u = ppls['uniform_skip']['avg']
    ppl_e = ppls['entropic']['avg']

    print(f"  Baseline ({NUM_BASELINE_STEPS} steps):  PPL = {ppl_b:.2f}")
    print(f"  Uniform  ({NUM_REDUCED_STEPS} steps):   PPL = {ppl_u:.2f} (delta: {ppl_u - ppl_b:+.2f})")
    print(f"  Entropic ({NUM_REDUCED_STEPS} steps):   PPL = {ppl_e:.2f} (delta: {ppl_e - ppl_b:+.2f})")

    if ppl_e < ppl_u:
        improvement = (ppl_u - ppl_e) / ppl_u * 100
        print(f"\n  >>> ENTROPIC WINS over uniform skip by {improvement:.1f}% PPL <<<")
        print(f"  The information-theoretic redistribution of time steps works!")
    elif ppl_e > ppl_u:
        degradation = (ppl_e - ppl_u) / ppl_u * 100
        print(f"\n  >>> Entropic LOSES to uniform skip by {degradation:.1f}% PPL <<<")
    else:
        print(f"\n  >>> Entropic TIES with uniform skip <<<")

    speedup_e = baseline_time / np.mean(all_results['entropic']['times'])
    print(f"\n  Speedup vs baseline: {speedup_e:.2f}x (with {NUM_REDUCED_STEPS}/{NUM_BASELINE_STEPS} steps)")

    # -- Print ALL generated text for visual inspection --
    print(f"\n{'='*70}")
    print("FULL GENERATED TEXT (for visual quality inspection)")
    print(f"{'='*70}")

    for method_name in ['baseline', 'uniform_skip', 'entropic']:
        print(f"\n--- {method_name.upper()} ---")
        for i, text in enumerate(all_results[method_name]['texts']):
            print(f"\n  [Sample {i+1}]:")
            # Print first 800 chars for readability
            print(f"  {text[:800]}")
            if len(text) > 800:
                print(f"  [...truncated, total {len(text)} chars]")
            print()

    # -- Save all results --
    save_data = {}
    for method_name, r in all_results.items():
        save_data[method_name] = {
            'steps': r['steps'],
            'avg_time': float(np.mean(r['times'])),
            'times': r['times'],
            'ppl_avg': r.get('ppl_avg', None),
            'ppl_per_sample': r.get('ppl_per_sample', None),
            'texts': r['texts'],
        }

    results_path = results_dir / "entropic_experiment_results.json"
    with open(results_path, 'w') as f:
        json.dump(save_data, f, indent=2, default=str)
    print(f"\nResults saved to {results_path}")

    print(f"\nDone!")


if __name__ == "__main__":
    main()
