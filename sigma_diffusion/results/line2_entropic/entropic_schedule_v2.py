#!/usr/bin/env python3
"""
Entropic Time Scheduler v2 for MDLM Diffusion Language Model
==============================================================
Key fix from v1: The D_KL between consecutive steps is near-zero in MDLM
because masked diffusion has a different information flow than continuous
diffusion. In masked diffusion:
  - t = 1.0: everything is masked
  - t = eps: everything is unmasked
  - The probability of unmasking = t (the masking probability at time t)

The CORRECT information measure for the entropic scheduler is:
  I(t) = |d/dt [entropy of the model's predicted distribution]| * unmask_rate(t)

In practice, we measure for each step:
  1. H(t): average entropy of model's softmax distribution over masked tokens
  2. unmask_count: number of tokens that got unmasked in this step
  3. confidence: 1 - H(t)/log(V) (normalized inverse entropy)
  4. info_content = unmask_count * confidence (= tokens of useful information revealed)

The entropic schedule places MORE steps where info_content is HIGH:
  - Near t=1: many tokens are masked but model is uncertain -> moderate info
  - Near t~0.3-0.5: many tokens unmask AND model is confident -> HIGH info
  - Near t=0: few tokens left to unmask -> low info

This creates a non-uniform schedule that concentrates steps in the
"sweet spot" where the model is both active and confident.
"""

import os
import sys
import time
import json
import math
from pathlib import Path

import numpy as np

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
# 1.  MPS-compatible MDLM model (identical to v1)
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
        freqs = torch.exp(-math.log(max_period) * torch.arange(start=0, end=half, dtype=torch.float32) / half).to(device=t.device)
        args = t[:, None].float() * freqs[None]
        embedding = torch.cat([torch.cos(args), torch.sin(args)], dim=-1)
        if dim % 2:
            embedding = torch.cat([embedding, torch.zeros_like(embedding[:, :1])], dim=-1)
        return embedding
    def forward(self, t):
        t_freq = self.timestep_embedding(t, self.frequency_embedding_size)
        return self.mlp(t_freq)


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
    return q * cos + rotate_half(q) * sin, k * cos + rotate_half(k) * sin


class DDiTBlock(nn.Module):
    def __init__(self, dim, n_heads, cond_dim, mlp_ratio=4, dropout=0.1):
        super().__init__()
        self.n_heads = n_heads
        self.head_dim = dim // n_heads
        self.norm1 = LayerNorm(dim)
        self.attn_qkv = nn.Linear(dim, 3 * dim, bias=False)
        self.attn_out = nn.Linear(dim, dim, bias=False)
        self.norm2 = LayerNorm(dim)
        self.mlp = nn.Sequential(
            nn.Linear(dim, mlp_ratio * dim, bias=True),
            nn.GELU(approximate='tanh'),
            nn.Linear(mlp_ratio * dim, dim, bias=True))
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
        qkv = rearrange(qkv, 'b s (three h d) -> three b h s d', three=3, h=self.n_heads)
        q, k, v = qkv[0], qkv[1], qkv[2]
        cos, sin = rotary_cos_sin
        q, k = apply_rotary_pos_emb_manual(q, k, cos, sin)
        attn_out = F.scaled_dot_product_attention(q, k, v, dropout_p=0.0)
        attn_out = rearrange(attn_out, 'b h s d -> b s (h d)')
        x = x_skip + gate_msa * F.dropout(self.attn_out(attn_out), p=self.dropout, training=self.training)
        x_norm2 = self.norm2(x) * (1 + scale_mlp) + shift_mlp
        x = x + gate_mlp * F.dropout(self.mlp(x_norm2), p=self.dropout, training=self.training)
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
        return self.linear(x)


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
            for _ in range(n_blocks)])
        self.output_layer = DDitFinalLayer(hidden_dim, vocab_size, cond_dim)
    def forward(self, indices, sigma):
        x = self.vocab_embed(indices)
        c = F.silu(self.sigma_map(sigma))
        cos, sin = self.rotary_emb(x)
        for block in self.blocks:
            x = block(x, (cos, sin), c)
        return self.output_layer(x, c)


def load_mdlm_weights(model):
    print("Downloading MDLM-OWT weights from HuggingFace...")
    path = hf_hub_download('kuleshov-group/mdlm-owt', 'model.safetensors')
    state_dict = {}
    with safe_open(path, framework='pt', device='cpu') as f:
        for key in f.keys():
            state_dict[key] = f.get_tensor(key)
    cleaned = {k.replace("backbone.", ""): v for k, v in state_dict.items()}
    missing, unexpected = model.load_state_dict(cleaned, strict=False)
    print(f"  Loaded {len(cleaned)} tensors. Missing: {len(missing)}, Unexpected: {len(unexpected)}")
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
    t_flat = t.squeeze(-1) if t.ndim > 1 else t
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


# =====================================================================
# 3.  Information-theoretic profiling for masked diffusion
# =====================================================================

@torch.no_grad()
def profile_information_flow(model, num_steps, seq_len=256, batch_size=1,
                              device=DEVICE):
    """
    Profile the ACTUAL information flow in MDLM generation.

    For each step, we measure:
    1. n_unmasked: how many tokens transitioned from masked -> unmasked
    2. avg_entropy: average entropy of model's prediction over masked tokens
    3. max_prob: average max probability (confidence) over masked tokens
    4. info_content: composite measure = n_unmasked * (1 - H/log(V))

    The entropic schedule should equalize info_content across steps.
    """
    x = MASK_INDEX * torch.ones(batch_size, seq_len, dtype=torch.int64, device=device)
    eps = 1e-5
    timesteps = torch.linspace(1, eps, num_steps + 1, device=device)
    dt = (1 - eps) / num_steps
    log_vocab = math.log(50257)  # exclude mask token

    profile = []
    print(f"  Profiling information flow over {num_steps} steps...")

    for i in range(num_steps):
        t = timesteps[i] * torch.ones(batch_size, 1, device=device)
        sigma, _ = loglinear_noise(t)

        # Get model predictions
        logits = model(x, sigma.squeeze(-1) if sigma.ndim > 1 else sigma)
        log_probs = subs_parameterization(logits.clone(), x)

        # Compute entropy over masked positions
        masked = (x == MASK_INDEX)
        n_masked_before = masked.float().sum().item()

        if n_masked_before > 0:
            probs = log_probs.exp()
            # Entropy: H = -sum(p * log(p)), only over masked positions
            # Avoid log(0) by using log_probs directly
            entropy_per_pos = -(probs * log_probs).sum(dim=-1)  # (B, L)
            entropy_masked = (entropy_per_pos * masked.float()).sum() / max(n_masked_before, 1)
            avg_entropy = entropy_masked.item()

            # Max probability (confidence)
            max_probs = probs.max(dim=-1).values  # (B, L)
            avg_max_prob = (max_probs * masked.float()).sum().item() / max(n_masked_before, 1)
        else:
            avg_entropy = 0.0
            avg_max_prob = 1.0

        # Execute the step
        p_x0 = log_probs.exp()
        _, x_next = ddpm_cache_update(model, x, t, dt, p_x0=p_x0)

        # Count newly unmasked tokens
        was_masked = (x == MASK_INDEX)
        now_unmasked = (x_next != MASK_INDEX) & was_masked
        n_unmasked = now_unmasked.float().sum().item()

        # Confidence score: 1 - normalized entropy
        confidence = max(0, 1.0 - avg_entropy / log_vocab)

        # Information content: tokens revealed * confidence
        info_content = n_unmasked * confidence

        # Also: a pure "transition density" measure
        # = how much the mask probability changes in this interval
        # In MDLM: mask_prob = t, so d(mask_prob)/d(step) = dt
        # But the ACTUAL unmask rate depends on the model's predictions
        unmask_rate = n_unmasked / max(n_masked_before, 1)

        profile.append({
            'step': i,
            'timestep': timesteps[i].item(),
            'n_masked_before': n_masked_before,
            'n_unmasked': n_unmasked,
            'avg_entropy': avg_entropy,
            'avg_max_prob': avg_max_prob,
            'confidence': confidence,
            'info_content': info_content,
            'unmask_rate': unmask_rate,
        })

        x = x_next

        if (i + 1) % 16 == 0:
            print(f"    Step {i+1:3d}/{num_steps}, t={timesteps[i].item():.4f}, "
                  f"masked={n_masked_before:.0f}, unmasked={n_unmasked:.0f}, "
                  f"H={avg_entropy:.3f}, conf={confidence:.3f}, "
                  f"info={info_content:.2f}")

    return profile


def build_entropic_schedule_v2(info_profile, num_entropic_steps, eps=1e-5):
    """
    Build entropic schedule from information flow profile.

    Uses info_content = n_unmasked * confidence as the quantity to equalize.
    Places MORE time points where info_content is HIGH.

    Since masked diffusion has a clear structure:
    - Early (t~1): many masked, low confidence -> moderate info
    - Middle (t~0.3-0.5): active unmasking with growing confidence -> HIGH info
    - Late (t~0): few masked left -> low info

    The schedule will have:
    - Moderate density near t=1
    - HIGH density in the "sweet spot" (t~0.2-0.5)
    - Low density near t=0
    """
    profile_sorted = sorted(info_profile, key=lambda s: s['timestep'], reverse=True)
    t_values = np.array([s['timestep'] for s in profile_sorted])
    info_values = np.array([s['info_content'] for s in profile_sorted])

    # Ensure non-negative
    info_values = np.maximum(info_values, 0.0)

    # Add a small floor to avoid zero-info regions getting no coverage at all
    # This ensures we don't completely skip any time region
    floor = np.mean(info_values) * 0.01
    info_values_smoothed = info_values + floor

    # Build CDF
    cumulative = np.cumsum(info_values_smoothed)
    total = cumulative[-1]

    if total < 1e-12:
        print("  WARNING: Total info near zero. Falling back to uniform schedule.")
        return torch.linspace(1, eps, num_entropic_steps + 1)

    cdf = cumulative / total

    # Target: num_entropic_steps + 1 evenly spaced CDF values
    target_cdf = np.linspace(0, 1, num_entropic_steps + 1)

    # Prepend start point
    cdf_full = np.concatenate([[0.0], cdf])
    t_full = np.concatenate([[t_values[0]], t_values])

    # Invert CDF
    entropic_t = np.interp(target_cdf, cdf_full, t_full)

    # Fix boundaries
    entropic_t[0] = t_values[0]
    entropic_t[-1] = max(entropic_t[-1], eps)

    # Ensure strictly decreasing
    for i in range(1, len(entropic_t)):
        if entropic_t[i] >= entropic_t[i-1]:
            entropic_t[i] = entropic_t[i-1] - 1e-7

    # Analysis
    dts = np.abs(np.diff(entropic_t))
    print(f"\n  Entropic Schedule v2 built ({num_entropic_steps} steps):")
    print(f"    Total info_content: {total:.2f}")
    print(f"    Info per entropic step: {total / num_entropic_steps:.4f}")
    print(f"    Time range: [{entropic_t[0]:.6f}, {entropic_t[-1]:.6f}]")
    print(f"    Smallest dt: {dts.min():.6f} (densest)")
    print(f"    Largest  dt: {dts.max():.6f} (sparsest)")
    print(f"    Ratio: {dts.max()/max(dts.min(), 1e-10):.1f}x")

    # Show density by region
    n = len(entropic_t) - 1
    q1 = sum(1 for t in entropic_t if t > 0.75) - 1
    q2 = sum(1 for t in entropic_t if 0.5 < t <= 0.75)
    q3 = sum(1 for t in entropic_t if 0.25 < t <= 0.5)
    q4 = sum(1 for t in entropic_t if t <= 0.25) - 1  # -1 for endpoint
    print(f"    Steps in t=[0.75,1.0]: {q1}")
    print(f"    Steps in t=[0.50,0.75]: {q2}")
    print(f"    Steps in t=[0.25,0.50]: {q3}")
    print(f"    Steps in t=[0.00,0.25]: {q4}")
    uniform_each = num_entropic_steps // 4
    print(f"    (Uniform would have ~{uniform_each} per quarter)")

    return torch.tensor(entropic_t, dtype=torch.float32)


# =====================================================================
# 4.  Generation functions
# =====================================================================

@torch.no_grad()
def generate_baseline(model, num_steps, seq_len=256, batch_size=1, device=DEVICE):
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
    t = timesteps[-1] * torch.ones(batch_size, 1, device=device)
    sigma_final, _ = loglinear_noise(t)
    logits = model(x, sigma_final.squeeze(-1) if sigma_final.ndim > 1 else sigma_final)
    x = logits.argmax(dim=-1)
    return x, time.time() - t0


@torch.no_grad()
def generate_uniform_skip(model, num_keep_steps, seq_len=256, batch_size=1, device=DEVICE):
    x = MASK_INDEX * torch.ones(batch_size, seq_len, dtype=torch.int64, device=device)
    eps = 1e-5
    timesteps = torch.linspace(1, eps, num_keep_steps + 1, device=device)
    dt = (1 - eps) / num_keep_steps
    p_x0_cache = None
    t0 = time.time()
    for i in range(num_keep_steps):
        t = timesteps[i] * torch.ones(batch_size, 1, device=device)
        p_x0_cache, x_next = ddpm_cache_update(model, x, t, dt, p_x0=p_x0_cache)
        if not torch.allclose(x_next, x):
            p_x0_cache = None
        x = x_next
    t = timesteps[-1] * torch.ones(batch_size, 1, device=device)
    sigma_final, _ = loglinear_noise(t)
    logits = model(x, sigma_final.squeeze(-1) if sigma_final.ndim > 1 else sigma_final)
    x = logits.argmax(dim=-1)
    return x, time.time() - t0


@torch.no_grad()
def generate_entropic(model, entropic_timesteps, seq_len=256, batch_size=1, device=DEVICE):
    """Generate with non-uniform time schedule. Each (t_i, t_{i+1}) is one step."""
    num_steps = len(entropic_timesteps) - 1
    x = MASK_INDEX * torch.ones(batch_size, seq_len, dtype=torch.int64, device=device)
    p_x0_cache = None
    t0 = time.time()
    for i in range(num_steps):
        t_val = entropic_timesteps[i]
        dt_val = entropic_timesteps[i] - entropic_timesteps[i + 1]
        t = t_val * torch.ones(batch_size, 1, device=device)
        p_x0_cache, x_next = ddpm_cache_update(model, x, t, dt_val.item(), p_x0=p_x0_cache)
        if not torch.allclose(x_next, x):
            p_x0_cache = None
        x = x_next
    t_final = entropic_timesteps[-1]
    t = t_final * torch.ones(batch_size, 1, device=device)
    sigma_final, _ = loglinear_noise(t)
    logits = model(x, sigma_final.squeeze(-1) if sigma_final.ndim > 1 else sigma_final)
    x = logits.argmax(dim=-1)
    return x, time.time() - t0


# =====================================================================
# 5.  Perplexity evaluation
# =====================================================================

@torch.no_grad()
def compute_perplexity_gpt2(text_samples, device=DEVICE):
    print(f"  Computing PPL with GPT-2...")
    tokenizer = transformers.AutoTokenizer.from_pretrained("gpt2")
    eval_model = transformers.AutoModelForCausalLM.from_pretrained("gpt2")
    if tokenizer.pad_token is None:
        tokenizer.pad_token = tokenizer.eos_token
    eval_model = eval_model.to(device).eval()

    encoded = tokenizer(text_samples, return_tensors='pt', truncation=True,
                        padding=True, max_length=1024, return_attention_mask=True)
    input_ids = encoded['input_ids'].to(device)
    attn_mask = encoded['attention_mask'].to(device)

    logits = eval_model(input_ids, attention_mask=attn_mask).logits
    shift_logits = logits[:, :-1, :].contiguous()
    shift_labels = input_ids[:, 1:].contiguous()
    shift_mask = attn_mask[:, 1:].contiguous()

    loss = F.cross_entropy(shift_logits.view(-1, shift_logits.size(-1)),
                           shift_labels.view(-1), reduction='none')
    loss = loss.view(shift_labels.shape)
    per_sample_nll = (loss * shift_mask).sum(dim=-1) / shift_mask.sum(dim=-1).clamp(min=1)
    per_sample_ppl = per_sample_nll.exp()

    avg_ppl = per_sample_ppl.mean().item()
    del eval_model
    if device == "mps":
        torch.mps.empty_cache()
    return avg_ppl, per_sample_ppl.tolist()


# =====================================================================
# 6.  Main experiment
# =====================================================================

def main():
    results_dir = Path(__file__).parent

    print("=" * 70)
    print("ENTROPIC TIME SCHEDULER v2 FOR MDLM")
    print("Information-flow based schedule for masked diffusion")
    print("=" * 70)

    # -- Load model --
    print("\n[1/7] Building MDLM model...")
    model = MDLMModel(vocab_size=50258, hidden_dim=768, cond_dim=128,
                      n_blocks=12, n_heads=12, dropout=0.1)
    load_mdlm_weights(model)
    model = model.to(DEVICE).eval()
    print(f"  {sum(p.numel() for p in model.parameters()) / 1e6:.1f}M parameters on {DEVICE}")

    # -- Tokenizer --
    print("\n[2/7] Loading GPT-2 tokenizer...")
    tokenizer = GPT2TokenizerFast.from_pretrained('gpt2')
    if tokenizer.pad_token is None:
        tokenizer.pad_token = tokenizer.eos_token

    # -- Config --
    NUM_BASELINE = 128
    NUM_REDUCED = 64
    SEQ_LEN = 256
    NUM_SAMPLES = 3
    print(f"\n  Baseline: {NUM_BASELINE} steps, Reduced: {NUM_REDUCED} steps")
    print(f"  Seq len: {SEQ_LEN}, Samples: {NUM_SAMPLES}")

    # -- Warmup --
    print("\n[3/7] Warmup...")
    _ = model(MASK_INDEX * torch.ones(1, 64, dtype=torch.int64, device=DEVICE),
              torch.tensor([0.5], device=DEVICE))
    print("  Done.")

    # -- Profile information flow --
    print(f"\n[4/7] Profiling information flow ({NUM_BASELINE} steps)...")
    info_profile = profile_information_flow(
        model, NUM_BASELINE, seq_len=SEQ_LEN, batch_size=1, device=DEVICE)

    # Save profile
    profile_path = results_dir / "info_profile_128steps.json"
    with open(profile_path, 'w') as f:
        json.dump(info_profile, f, indent=2)
    print(f"  Saved to {profile_path}")

    # Analyze
    infos = [s['info_content'] for s in info_profile]
    unmasks = [s['n_unmasked'] for s in info_profile]
    confs = [s['confidence'] for s in info_profile]
    print(f"\n  Information flow statistics:")
    print(f"    Total info_content: {sum(infos):.2f}")
    print(f"    Total tokens unmasked: {sum(unmasks):.0f} / {SEQ_LEN}")
    print(f"    Info range: [{min(infos):.4f}, {max(infos):.4f}]")
    print(f"    Confidence range: [{min(confs):.4f}, {max(confs):.4f}]")

    # By quarter
    n = len(infos)
    for qi, (lo, hi, label) in enumerate([
        (0, n//4, "Q1 t=[1.0,0.75]"),
        (n//4, n//2, "Q2 t=[0.75,0.50]"),
        (n//2, 3*n//4, "Q3 t=[0.50,0.25]"),
        (3*n//4, n, "Q4 t=[0.25,0.00]"),
    ]):
        q_info = sum(infos[lo:hi])
        q_unmask = sum(unmasks[lo:hi])
        print(f"    {label}: info={q_info:.2f}, unmask={q_unmask:.0f}")

    # -- Build entropic schedule --
    print(f"\n[5/7] Building entropic schedule ({NUM_REDUCED} steps)...")
    entropic_timesteps = build_entropic_schedule_v2(info_profile, NUM_REDUCED)
    entropic_timesteps = entropic_timesteps.to(DEVICE)

    schedule_path = results_dir / "entropic_schedule_v2.json"
    with open(schedule_path, 'w') as f:
        json.dump({
            'timesteps': entropic_timesteps.cpu().tolist(),
            'num_steps': NUM_REDUCED,
            'method': 'info_content = n_unmasked * confidence'
        }, f, indent=2)

    # -- Generate samples --
    print(f"\n[6/7] Generating {NUM_SAMPLES} samples per method...")

    all_results = {
        'baseline': {'texts': [], 'times': [], 'steps': NUM_BASELINE},
        'uniform_skip': {'texts': [], 'times': [], 'steps': NUM_REDUCED},
        'entropic': {'texts': [], 'times': [], 'steps': NUM_REDUCED},
    }

    for si in range(NUM_SAMPLES):
        seed = 42 + si * 1000
        print(f"\n{'='*60}")
        print(f"SAMPLE {si+1}/{NUM_SAMPLES} (seed={seed})")
        print(f"{'='*60}")

        # A: Baseline
        torch.manual_seed(seed)
        if DEVICE == "mps": torch.mps.manual_seed(seed)
        tokens, elapsed = generate_baseline(model, NUM_BASELINE, seq_len=SEQ_LEN,
                                            batch_size=1, device=DEVICE)
        text = tokenizer.decode(tokens[0].cpu(), skip_special_tokens=True)
        all_results['baseline']['texts'].append(text)
        all_results['baseline']['times'].append(elapsed)
        print(f"\n  A. BASELINE ({NUM_BASELINE} steps, {elapsed:.2f}s):")
        print(f"  {text[:600]}")

        # B: Uniform skip
        torch.manual_seed(seed)
        if DEVICE == "mps": torch.mps.manual_seed(seed)
        tokens, elapsed = generate_uniform_skip(model, NUM_REDUCED, seq_len=SEQ_LEN,
                                                 batch_size=1, device=DEVICE)
        text = tokenizer.decode(tokens[0].cpu(), skip_special_tokens=True)
        all_results['uniform_skip']['texts'].append(text)
        all_results['uniform_skip']['times'].append(elapsed)
        print(f"\n  B. UNIFORM ({NUM_REDUCED} steps, {elapsed:.2f}s):")
        print(f"  {text[:600]}")

        # C: Entropic
        torch.manual_seed(seed)
        if DEVICE == "mps": torch.mps.manual_seed(seed)
        tokens, elapsed = generate_entropic(model, entropic_timesteps, seq_len=SEQ_LEN,
                                             batch_size=1, device=DEVICE)
        text = tokenizer.decode(tokens[0].cpu(), skip_special_tokens=True)
        all_results['entropic']['texts'].append(text)
        all_results['entropic']['times'].append(elapsed)
        print(f"\n  C. ENTROPIC ({NUM_REDUCED} non-uniform steps, {elapsed:.2f}s):")
        print(f"  {text[:600]}")

    # -- Perplexity --
    print(f"\n{'='*60}")
    print(f"[7/7] PERPLEXITY EVALUATION")
    print(f"{'='*60}")

    ppls = {}
    for method in ['baseline', 'uniform_skip', 'entropic']:
        texts = all_results[method]['texts']
        avg_ppl, per_ppl = compute_perplexity_gpt2(texts, device=DEVICE)
        ppls[method] = avg_ppl
        all_results[method]['ppl_avg'] = avg_ppl
        all_results[method]['ppl_per_sample'] = per_ppl
        print(f"  {method:20s}: avg PPL = {avg_ppl:8.2f}  per-sample = {[f'{p:.1f}' for p in per_ppl]}")

    # -- Summary --
    print(f"\n{'='*70}")
    print("RESULTS SUMMARY")
    print(f"{'='*70}")
    base_time = np.mean(all_results['baseline']['times'])
    print(f"{'Method':<20} {'Steps':<8} {'Time(s)':<10} {'PPL':<12} {'Speedup':<10}")
    print("-" * 60)
    for method in ['baseline', 'uniform_skip', 'entropic']:
        r = all_results[method]
        t_avg = np.mean(r['times'])
        speedup = base_time / t_avg
        print(f"{method:<20} {r['steps']:<8} {t_avg:<10.2f} {r['ppl_avg']:<12.2f} {speedup:<10.2f}x")

    # -- Verdict --
    ppl_b = ppls['baseline']
    ppl_u = ppls['uniform_skip']
    ppl_e = ppls['entropic']

    print(f"\n{'='*70}")
    print("VERDICT")
    print(f"{'='*70}")
    print(f"  Baseline  ({NUM_BASELINE} steps): PPL = {ppl_b:.2f}")
    print(f"  Uniform   ({NUM_REDUCED} steps):  PPL = {ppl_u:.2f}  (delta from baseline: {ppl_u - ppl_b:+.2f})")
    print(f"  Entropic  ({NUM_REDUCED} steps):  PPL = {ppl_e:.2f}  (delta from baseline: {ppl_e - ppl_b:+.2f})")

    if ppl_e < ppl_u:
        imp = (ppl_u - ppl_e) / ppl_u * 100
        print(f"\n  >>> ENTROPIC WINS over uniform skip by {imp:.1f}% PPL <<<")
        if ppl_e < ppl_b:
            print(f"  >>> ENTROPIC even BEATS the {NUM_BASELINE}-step baseline! <<<")
    elif ppl_e > ppl_u:
        deg = (ppl_e - ppl_u) / ppl_u * 100
        print(f"\n  >>> Entropic LOSES to uniform skip by {deg:.1f}% PPL <<<")
    else:
        print(f"\n  >>> Entropic TIES with uniform skip <<<")

    # Quality assessment
    print(f"\n{'='*70}")
    print("TEXT QUALITY ASSESSMENT")
    print(f"{'='*70}")
    for method in ['baseline', 'uniform_skip', 'entropic']:
        texts = all_results[method]['texts']
        ppls_list = all_results[method]['ppl_per_sample']
        print(f"\n--- {method.upper()} (PPL={all_results[method]['ppl_avg']:.1f}) ---")
        for i, (text, ppl) in enumerate(zip(texts, ppls_list)):
            # Simple quality metrics
            words = text.split()
            unique_ratio = len(set(words)) / max(len(words), 1)
            avg_word_len = np.mean([len(w) for w in words]) if words else 0
            print(f"\n  Sample {i+1} (PPL={ppl:.1f}, {len(words)} words, "
                  f"unique ratio={unique_ratio:.2f}):")
            print(f"  {text[:500]}")
            if len(text) > 500:
                print(f"  [...{len(text)} chars total]")

    # -- Save --
    save_data = {}
    for method, r in all_results.items():
        save_data[method] = {
            'steps': r['steps'],
            'avg_time': float(np.mean(r['times'])),
            'times': r['times'],
            'ppl_avg': r.get('ppl_avg'),
            'ppl_per_sample': r.get('ppl_per_sample'),
            'texts': r['texts'],
        }
    results_path = results_dir / "entropic_v2_results.json"
    with open(results_path, 'w') as f:
        json.dump(save_data, f, indent=2, default=str)
    print(f"\nResults saved to {results_path}")
    print("Done!")


if __name__ == "__main__":
    main()
