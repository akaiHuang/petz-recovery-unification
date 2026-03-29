#!/usr/bin/env python3
"""
Entropic Time Scheduler v3 for MDLM
====================================
Key improvements over v2:
1. Use ANALYTICAL expected unmask rate instead of stochastic single-trajectory count
2. Average the confidence profile over MULTIPLE profiling runs
3. Use smoothed (Gaussian-filtered) confidence profile to reduce noise
4. Compare multiple schedule strategies:
   - Uniform: t_i = 1 - i/N (baseline for 64 steps)
   - Entropic-analytical: equal dSigma per step using analytical MDLM rate
   - Entropic-empirical: equal info using measured confidence + analytical rate
   - Cosine: t_i = cos(pi*i/(2N))^2 (common in image diffusion)

MDLM Theory:
- Masking probability at time t: mask_prob(t) = t
- Expected number of masked tokens: n_mask(t) = L * t  (where L = seq_len)
- Expected unmask rate: dn_unmask/dt = L  (constant! for linear schedule)
- But the NOISE schedule is loglinear: sigma(t) = -log(1 - (1-eps)*t)
  So dsigma/dt = (1-eps)/(1-(1-eps)*t), which increases as t->0
- The model's confidence (1 - H/log(V)) increases as t decreases
- Information rate = dsigma/dt * confidence(t)
"""

import os
import sys
import time
import json
import math
from pathlib import Path

import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.interpolate import interp1d

os.environ["HF_TOKEN"] = "YOUR_HF_TOKEN_HERE"
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
# Model definition (same as v1/v2)
# =====================================================================

DEVICE = "mps" if torch.backends.mps.is_available() else "cpu"
DTYPE = torch.float32
print(f"Device: {DEVICE}")


class LayerNorm(nn.Module):
    def __init__(self, dim):
        super().__init__()
        self.weight = nn.Parameter(torch.ones([dim]))
        self.dim = dim
    def forward(self, x):
        return F.layer_norm(x.float(), [self.dim]) * self.weight[None, None, :]

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
        freqs = torch.exp(-math.log(max_period) * torch.arange(0, half, dtype=torch.float32) / half).to(t.device)
        args = t[:, None].float() * freqs[None]
        emb = torch.cat([torch.cos(args), torch.sin(args)], dim=-1)
        if dim % 2: emb = torch.cat([emb, torch.zeros_like(emb[:, :1])], dim=-1)
        return emb
    def forward(self, t):
        return self.mlp(self.timestep_embedding(t, self.frequency_embedding_size))

class EmbeddingLayer(nn.Module):
    def __init__(self, dim, vocab_dim):
        super().__init__()
        self.embedding = nn.Parameter(torch.empty((vocab_dim, dim)))
        torch.nn.init.kaiming_uniform_(self.embedding, a=math.sqrt(5))
    def forward(self, x): return self.embedding[x]

class Rotary(nn.Module):
    def __init__(self, dim, base=10_000):
        super().__init__()
        self.register_buffer('inv_freq', 1.0 / (base ** (torch.arange(0, dim, 2).float() / dim)))
        self.seq_len_cached = self.cos_cached = self.sin_cached = None
    def forward(self, x, seq_dim=1):
        seq_len = x.shape[seq_dim]
        if seq_len != self.seq_len_cached:
            self.seq_len_cached = seq_len
            t = torch.arange(seq_len, device=x.device).type_as(self.inv_freq)
            freqs = torch.einsum("i,j->ij", t, self.inv_freq.clone())
            emb = torch.cat((freqs, freqs), dim=-1).to(x.device)
            self.cos_cached, self.sin_cached = emb.cos(), emb.sin()
        return self.cos_cached, self.sin_cached

def apply_rotary(q, k, cos, sin):
    cos, sin = cos[None, None, :, :], sin[None, None, :, :]
    def rot(x):
        x1, x2 = x[..., :x.shape[-1]//2], x[..., x.shape[-1]//2:]
        return torch.cat((-x2, x1), dim=-1)
    return q*cos + rot(q)*sin, k*cos + rot(k)*sin

class DDiTBlock(nn.Module):
    def __init__(self, dim, n_heads, cond_dim, mlp_ratio=4, dropout=0.1):
        super().__init__()
        self.n_heads, self.dropout = n_heads, dropout
        self.norm1, self.norm2 = LayerNorm(dim), LayerNorm(dim)
        self.attn_qkv = nn.Linear(dim, 3*dim, bias=False)
        self.attn_out = nn.Linear(dim, dim, bias=False)
        self.mlp = nn.Sequential(nn.Linear(dim, mlp_ratio*dim), nn.GELU(approximate='tanh'), nn.Linear(mlp_ratio*dim, dim))
        self.adaLN_modulation = nn.Linear(cond_dim, 6*dim)
        self.adaLN_modulation.weight.data.zero_(); self.adaLN_modulation.bias.data.zero_()
    def forward(self, x, rcs, c):
        B, S, D = x.shape
        sh1,sc1,g1,sh2,sc2,g2 = self.adaLN_modulation(c)[:,None].chunk(6, dim=2)
        xn = self.norm1(x)*(1+sc1)+sh1
        qkv = rearrange(self.attn_qkv(xn), 'b s (t h d)->t b h s d', t=3, h=self.n_heads)
        q, k = apply_rotary(qkv[0], qkv[1], *rcs)
        a = F.scaled_dot_product_attention(q, k, qkv[2], dropout_p=0.0)
        x = x + g1*F.dropout(self.attn_out(rearrange(a,'b h s d->b s (h d)')), p=self.dropout, training=self.training)
        xn2 = self.norm2(x)*(1+sc2)+sh2
        x = x + g2*F.dropout(self.mlp(xn2), p=self.dropout, training=self.training)
        return x

class DDitFinalLayer(nn.Module):
    def __init__(self, hs, oc, cd):
        super().__init__()
        self.norm_final = LayerNorm(hs)
        self.linear = nn.Linear(hs, oc); self.linear.weight.data.zero_(); self.linear.bias.data.zero_()
        self.adaLN_modulation = nn.Linear(cd, 2*hs); self.adaLN_modulation.weight.data.zero_(); self.adaLN_modulation.bias.data.zero_()
    def forward(self, x, c):
        sh, sc = self.adaLN_modulation(c)[:,None].chunk(2, dim=2)
        return self.linear(self.norm_final(x)*(1+sc)+sh)

class MDLMModel(nn.Module):
    def __init__(self, vocab_size=50258, hidden_dim=768, cond_dim=128, n_blocks=12, n_heads=12, dropout=0.1):
        super().__init__()
        self.vocab_size = vocab_size
        self.vocab_embed = EmbeddingLayer(hidden_dim, vocab_size)
        self.sigma_map = TimestepEmbedder(cond_dim)
        self.rotary_emb = Rotary(hidden_dim // n_heads)
        self.blocks = nn.ModuleList([DDiTBlock(hidden_dim, n_heads, cond_dim, dropout=dropout) for _ in range(n_blocks)])
        self.output_layer = DDitFinalLayer(hidden_dim, vocab_size, cond_dim)
    def forward(self, idx, sigma):
        x = self.vocab_embed(idx)
        c = F.silu(self.sigma_map(sigma))
        rcs = self.rotary_emb(x)
        for b in self.blocks: x = b(x, rcs, c)
        return self.output_layer(x, c)

def load_mdlm_weights(model):
    path = hf_hub_download('kuleshov-group/mdlm-owt', 'model.safetensors')
    sd = {}
    with safe_open(path, framework='pt', device='cpu') as f:
        for k in f.keys(): sd[k] = f.get_tensor(k)
    cleaned = {k.replace("backbone.",""): v for k,v in sd.items()}
    model.load_state_dict(cleaned, strict=False)
    print(f"  Loaded {len(cleaned)} weight tensors")
    return model

# =====================================================================
# Diffusion primitives
# =====================================================================

MASK_INDEX = 50257
NEG_INF = -1e6

def subs_param(logits, xt):
    logits[:,:,MASK_INDEX] += NEG_INF
    logits = logits - torch.logsumexp(logits, dim=-1, keepdim=True)
    u = (xt != MASK_INDEX)
    logits[u] = NEG_INF; logits[u, xt[u]] = 0
    return logits

def sample_cat(probs):
    g = 1e-10 - (torch.rand_like(probs)+1e-10).log()
    return (probs / g).argmax(dim=-1)

def loglinear_noise(t, eps=1e-3):
    return -torch.log1p(-(1-eps)*t), (1-eps)/(1-(1-eps)*t)

@torch.no_grad()
def ddpm_update(model, x, t, dt, p_x0=None):
    sigma, _ = loglinear_noise(t)
    tf = t.squeeze(-1) if t.ndim > 1 else t
    mt, ms = tf[:,None,None], (tf-dt)[:,None,None]
    if p_x0 is None:
        p_x0 = subs_param(model(x, sigma.squeeze(-1) if sigma.ndim>1 else sigma), x).exp()
    q = p_x0 * (mt - ms)
    q[:,:,MASK_INDEX] = ms[:,:,0]
    _x = sample_cat(q)
    cf = (x != MASK_INDEX).to(x.dtype)
    return p_x0, cf*x + (1-cf)*_x

# =====================================================================
# Profiling: measure confidence as function of t
# =====================================================================

@torch.no_grad()
def profile_confidence(model, num_steps, seq_len=256, batch_size=1, device=DEVICE):
    """Profile model confidence (1 - H/logV) at each timestep.
    Returns: list of (timestep, confidence, n_masked)"""
    x = MASK_INDEX * torch.ones(batch_size, seq_len, dtype=torch.int64, device=device)
    eps = 1e-5
    ts = torch.linspace(1, eps, num_steps+1, device=device)
    dt = (1-eps)/num_steps
    logV = math.log(50257)
    results = []

    for i in range(num_steps):
        t = ts[i] * torch.ones(batch_size, 1, device=device)
        sigma, _ = loglinear_noise(t)
        logits = model(x, sigma.squeeze(-1) if sigma.ndim>1 else sigma)
        lp = subs_param(logits.clone(), x)
        masked = (x == MASK_INDEX)
        nm = masked.float().sum().item()
        if nm > 0:
            p = lp.exp()
            H = -(p * lp).sum(dim=-1)
            H_avg = (H * masked.float()).sum().item() / nm
            conf = max(0, 1 - H_avg/logV)
        else:
            conf = 1.0
        results.append({'step': i, 't': ts[i].item(), 'confidence': conf, 'n_masked': nm})
        # advance
        _, x = ddpm_update(model, x, t, dt, p_x0=lp.exp())
        if (i+1) % 32 == 0:
            print(f"    Step {i+1}/{num_steps}: t={ts[i].item():.4f}, conf={conf:.4f}, masked={nm:.0f}")

    return results


# =====================================================================
# Schedule builders
# =====================================================================

def build_uniform_schedule(N, eps=1e-5):
    """Standard uniform schedule."""
    return torch.linspace(1, eps, N+1)


def build_cosine_schedule(N, eps=1e-5):
    """Cosine schedule: t_i = cos^2(pi*i/(2N))."""
    i = torch.arange(N+1, dtype=torch.float64)
    t = torch.cos(math.pi * i / (2*N))**2
    t = t * (1 - eps) + eps  # scale to [eps, 1]
    t[0] = 1.0; t[-1] = eps
    return t.float()


def build_entropic_analytical(N, eps=1e-5):
    """
    Analytical entropic schedule for MDLM loglinear noise.

    The loglinear noise schedule has dsigma/dt = (1-eps)/(1-(1-eps)*t).
    To equalize dsigma per step, we need equal increments in sigma space.

    sigma(t) = -log(1 - (1-eps)*t)
    sigma(1) = -log(eps) ~ 6.9  (for eps=1e-3)
    sigma(eps) = -log(1 - (1-eps)*eps) ~ eps

    Uniform in sigma: sigma_i = sigma_max * i/N
    Invert: t_i = (1 - exp(-sigma_i)) / (1-eps)
    """
    sigma_max = -math.log(1 - (1-1e-3)*1.0)  # = -log(eps_noise) for t=1
    sigma_min = -math.log(1 - (1-1e-3)*eps)   # for t=eps

    sigma_vals = np.linspace(sigma_max, sigma_min, N+1)
    t_vals = (1 - np.exp(-sigma_vals)) / (1 - 1e-3)

    t_vals[0] = 1.0; t_vals[-1] = eps
    # ensure strictly decreasing
    for i in range(1, len(t_vals)):
        if t_vals[i] >= t_vals[i-1]:
            t_vals[i] = t_vals[i-1] - 1e-8
    return torch.tensor(t_vals, dtype=torch.float32)


def build_entropic_empirical(confidence_profile, N, eps=1e-5):
    """
    Entropic schedule using empirical confidence profile.

    Information rate at time t:
      I(t) = dsigma/dt * confidence(t)
           = [(1-eps_noise)/(1-(1-eps_noise)*t)] * confidence(t)

    We want each step to carry equal integral of I(t)dt.
    """
    # Extract and smooth confidence
    prof = sorted(confidence_profile, key=lambda s: s['t'], reverse=True)
    t_arr = np.array([s['t'] for s in prof])
    c_arr = np.array([s['confidence'] for s in prof])

    # Smooth with Gaussian filter (sigma=3 steps) to reduce noise
    c_smooth = gaussian_filter1d(c_arr, sigma=3)
    c_smooth = np.maximum(c_smooth, 0.01)

    # Compute I(t) = dsigma/dt * confidence(t) at each profiled point
    eps_noise = 1e-3
    dsigma_dt = (1 - eps_noise) / (1 - (1 - eps_noise) * t_arr)
    info_rate = dsigma_dt * c_smooth

    # Integrate: cumulative info from t=1 downward
    # Since t_arr is decreasing, dt is negative. We use |dt|.
    dt_arr = np.abs(np.diff(t_arr))
    # Trapezoidal integration
    info_increments = 0.5 * (info_rate[:-1] + info_rate[1:]) * dt_arr
    cum_info = np.concatenate([[0], np.cumsum(info_increments)])
    total_info = cum_info[-1]

    if total_info < 1e-12:
        print("  WARNING: Total info near zero, using uniform.")
        return build_uniform_schedule(N, eps)

    # Normalize CDF
    cdf = cum_info / total_info

    # Target: equally spaced CDF values
    target = np.linspace(0, 1, N+1)

    # Invert
    entropic_t = np.interp(target, cdf, t_arr)
    entropic_t[0] = 1.0
    entropic_t[-1] = max(eps, entropic_t[-1])

    # Ensure strictly decreasing
    for i in range(1, len(entropic_t)):
        if entropic_t[i] >= entropic_t[i-1]:
            entropic_t[i] = entropic_t[i-1] - 1e-8

    return torch.tensor(entropic_t, dtype=torch.float32)


def analyze_schedule(name, ts, N):
    """Print schedule statistics."""
    ts_np = ts.numpy() if isinstance(ts, torch.Tensor) else np.array(ts)
    dts = np.abs(np.diff(ts_np))
    # Count steps per quarter
    q1 = sum(1 for t in ts_np if t > 0.75) - 1
    q2 = sum(1 for t in ts_np if 0.5 < t <= 0.75)
    q3 = sum(1 for t in ts_np if 0.25 < t <= 0.5)
    q4 = sum(1 for t in ts_np if t <= 0.25) - 1
    print(f"  {name}:")
    print(f"    dt range: [{dts.min():.6f}, {dts.max():.6f}], ratio={dts.max()/max(dts.min(),1e-10):.1f}x")
    print(f"    Steps by quarter: Q1={q1}, Q2={q2}, Q3={q3}, Q4={q4}")
    return q1, q2, q3, q4


# =====================================================================
# Generation
# =====================================================================

@torch.no_grad()
def generate_with_schedule(model, timesteps_tensor, seq_len=256, batch_size=1, device=DEVICE):
    """Generate using arbitrary time schedule (tensor of decreasing t values)."""
    N = len(timesteps_tensor) - 1
    x = MASK_INDEX * torch.ones(batch_size, seq_len, dtype=torch.int64, device=device)
    p_x0_cache = None
    t0 = time.time()
    for i in range(N):
        t_val = timesteps_tensor[i]
        dt_val = (timesteps_tensor[i] - timesteps_tensor[i+1]).item()
        t = t_val * torch.ones(batch_size, 1, device=device)
        p_x0_cache, x_next = ddpm_update(model, x, t, dt_val, p_x0=p_x0_cache)
        if not torch.allclose(x_next, x):
            p_x0_cache = None
        x = x_next
    # Final
    t = timesteps_tensor[-1] * torch.ones(batch_size, 1, device=device)
    sf, _ = loglinear_noise(t)
    logits = model(x, sf.squeeze(-1) if sf.ndim>1 else sf)
    x = logits.argmax(dim=-1)
    return x, time.time() - t0


# =====================================================================
# PPL
# =====================================================================

@torch.no_grad()
def compute_ppl(texts, device=DEVICE):
    tok = transformers.AutoTokenizer.from_pretrained("gpt2")
    mdl = transformers.AutoModelForCausalLM.from_pretrained("gpt2").to(device).eval()
    if tok.pad_token is None: tok.pad_token = tok.eos_token
    enc = tok(texts, return_tensors='pt', truncation=True, padding=True, max_length=1024, return_attention_mask=True)
    ids, mask = enc['input_ids'].to(device), enc['attention_mask'].to(device)
    logits = mdl(ids, attention_mask=mask).logits
    sl, lab, sm = logits[:,:-1,:].contiguous(), ids[:,1:].contiguous(), mask[:,1:].contiguous()
    loss = F.cross_entropy(sl.view(-1,sl.size(-1)), lab.view(-1), reduction='none').view(lab.shape)
    nll = (loss*sm).sum(-1) / sm.sum(-1).clamp(min=1)
    ppl = nll.exp()
    avg = ppl.mean().item()
    del mdl
    if device=="mps": torch.mps.empty_cache()
    return avg, ppl.tolist()


# =====================================================================
# Main
# =====================================================================

def main():
    rd = Path(__file__).parent

    print("="*70)
    print("ENTROPIC TIME SCHEDULER v3 FOR MDLM")
    print("Analytical + empirical information-rate schedule")
    print("="*70)

    # Model
    print("\n[1] Loading model...")
    model = MDLMModel()
    load_mdlm_weights(model)
    model = model.to(DEVICE).eval()
    print(f"  {sum(p.numel() for p in model.parameters())/1e6:.1f}M params")

    # Tokenizer
    tokenizer = GPT2TokenizerFast.from_pretrained('gpt2')
    if tokenizer.pad_token is None: tokenizer.pad_token = tokenizer.eos_token

    # Config
    N_BASE = 128
    N_RED = 64
    SEQ = 256
    N_SAMP = 3

    # Warmup
    _ = model(MASK_INDEX*torch.ones(1,64,dtype=torch.int64,device=DEVICE), torch.tensor([0.5],device=DEVICE))

    # Profile confidence (average over 3 runs)
    print(f"\n[2] Profiling confidence ({N_BASE} steps, 3 runs)...")
    all_confs = []
    for run in range(3):
        print(f"  Run {run+1}/3:")
        torch.manual_seed(100 + run*777)
        if DEVICE=="mps": torch.mps.manual_seed(100 + run*777)
        prof = profile_confidence(model, N_BASE, seq_len=SEQ, device=DEVICE)
        all_confs.append(prof)

    # Average confidence profiles
    avg_profile = []
    for i in range(N_BASE):
        t_val = all_confs[0][i]['t']
        c_avg = np.mean([all_confs[r][i]['confidence'] for r in range(3)])
        nm_avg = np.mean([all_confs[r][i]['n_masked'] for r in range(3)])
        avg_profile.append({'step': i, 't': t_val, 'confidence': c_avg, 'n_masked': nm_avg})

    # Save
    with open(rd / "confidence_profile_avg.json", 'w') as f:
        json.dump(avg_profile, f, indent=2)

    # Print confidence by quarter
    print("\n  Averaged confidence by quarter:")
    for qi, (lo, hi, label) in enumerate([
        (0, 32, "Q1 t=[1.0,0.75]"), (32, 64, "Q2 t=[0.75,0.50]"),
        (64, 96, "Q3 t=[0.50,0.25]"), (96, 128, "Q4 t=[0.25,0.00]")]):
        cs = [avg_profile[i]['confidence'] for i in range(lo, hi)]
        print(f"    {label}: avg_conf={np.mean(cs):.4f} (range [{min(cs):.4f}, {max(cs):.4f}])")

    # Build schedules
    print(f"\n[3] Building schedules ({N_RED} steps)...")
    schedules = {
        'baseline_128': build_uniform_schedule(N_BASE),
        'uniform_64': build_uniform_schedule(N_RED),
        'cosine_64': build_cosine_schedule(N_RED),
        'entropic_analytical_64': build_entropic_analytical(N_RED),
        'entropic_empirical_64': build_entropic_empirical(avg_profile, N_RED),
    }

    for name, ts in schedules.items():
        analyze_schedule(name, ts, N_RED if '64' in name else N_BASE)

    # Save schedules
    for name, ts in schedules.items():
        with open(rd / f"schedule_{name}.json", 'w') as f:
            json.dump({'timesteps': ts.tolist(), 'name': name}, f, indent=2)

    # Generate
    print(f"\n[4] Generating {N_SAMP} samples x {len(schedules)} methods...")

    results = {name: {'texts': [], 'times': [], 'steps': len(ts)-1} for name, ts in schedules.items()}

    for si in range(N_SAMP):
        seed = 42 + si * 1000
        print(f"\n{'='*60}")
        print(f"SAMPLE {si+1}/{N_SAMP} (seed={seed})")
        print(f"{'='*60}")

        for name, ts in schedules.items():
            torch.manual_seed(seed)
            if DEVICE=="mps": torch.mps.manual_seed(seed)
            ts_dev = ts.to(DEVICE)
            tokens, elapsed = generate_with_schedule(model, ts_dev, seq_len=SEQ, device=DEVICE)
            text = tokenizer.decode(tokens[0].cpu(), skip_special_tokens=True)
            results[name]['texts'].append(text)
            results[name]['times'].append(elapsed)
            n_steps = len(ts) - 1
            print(f"\n  {name} ({n_steps} steps, {elapsed:.2f}s):")
            print(f"  {text[:400]}")

    # PPL
    print(f"\n{'='*60}")
    print("[5] PERPLEXITY EVALUATION")
    print(f"{'='*60}")

    for name in results:
        texts = results[name]['texts']
        avg_ppl, per_ppl = compute_ppl(texts, device=DEVICE)
        results[name]['ppl_avg'] = avg_ppl
        results[name]['ppl_per_sample'] = per_ppl
        print(f"  {name:30s}: PPL={avg_ppl:8.2f}  {[f'{p:.0f}' for p in per_ppl]}")

    # Summary
    print(f"\n{'='*70}")
    print("RESULTS SUMMARY")
    print(f"{'='*70}")
    base_time = np.mean(results['baseline_128']['times'])
    print(f"{'Method':<32} {'Steps':<7} {'Time':<8} {'PPL':<10} {'Speed':<8}")
    print("-"*65)
    for name in ['baseline_128', 'uniform_64', 'cosine_64',
                  'entropic_analytical_64', 'entropic_empirical_64']:
        r = results[name]
        t_avg = np.mean(r['times'])
        speedup = base_time / t_avg
        print(f"{name:<32} {r['steps']:<7} {t_avg:<8.2f} {r['ppl_avg']:<10.2f} {speedup:<8.2f}x")

    # Verdict
    ppl_base = results['baseline_128']['ppl_avg']
    ppl_uni = results['uniform_64']['ppl_avg']

    print(f"\n{'='*70}")
    print("DETAILED COMPARISON vs UNIFORM-64")
    print(f"{'='*70}")
    for name in ['cosine_64', 'entropic_analytical_64', 'entropic_empirical_64']:
        ppl = results[name]['ppl_avg']
        delta = ppl - ppl_uni
        pct = delta / ppl_uni * 100
        sign = "BETTER" if delta < 0 else "WORSE"
        print(f"  {name}: PPL={ppl:.2f} vs uniform {ppl_uni:.2f} -> {sign} by {abs(pct):.1f}%")

    # Best method
    methods_64 = ['uniform_64', 'cosine_64', 'entropic_analytical_64', 'entropic_empirical_64']
    best = min(methods_64, key=lambda n: results[n]['ppl_avg'])
    print(f"\n  BEST 64-step method: {best} (PPL={results[best]['ppl_avg']:.2f})")
    print(f"  vs baseline 128-step: delta = {results[best]['ppl_avg'] - ppl_base:+.2f}")

    # Full text output
    print(f"\n{'='*70}")
    print("FULL TEXT SAMPLES")
    print(f"{'='*70}")
    for name in ['baseline_128', 'uniform_64', 'entropic_empirical_64']:
        r = results[name]
        print(f"\n--- {name.upper()} (PPL={r['ppl_avg']:.1f}) ---")
        for i, (text, ppl) in enumerate(zip(r['texts'], r['ppl_per_sample'])):
            words = text.split()
            print(f"\n  [Sample {i+1}] PPL={ppl:.1f}, {len(words)} words:")
            print(f"  {text[:600]}")
            if len(text)>600: print(f"  [...{len(text)} chars total]")

    # Save
    save = {}
    for name, r in results.items():
        save[name] = {
            'steps': r['steps'],
            'avg_time': float(np.mean(r['times'])),
            'ppl_avg': r.get('ppl_avg'),
            'ppl_per_sample': r.get('ppl_per_sample'),
            'texts': r['texts'],
        }
    with open(rd / "entropic_v3_results.json", 'w') as f:
        json.dump(save, f, indent=2)
    print(f"\nSaved to {rd / 'entropic_v3_results.json'}")
    print("Done!")


if __name__ == "__main__":
    main()
