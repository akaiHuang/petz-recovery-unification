# Speed Audit Results

**Date**: 2026-03-23
**Platform**: Apple M1 Max, 10 CPU cores, macOS Darwin 25.1.0
**Configuration**: N_k=500, k=[3e-4, 0.35] Mpc^-1, RECFAST recombination, lg_max=ln_max=25 (nvar=56)

---

## 1. Profile Breakdown (Full scipy Pipeline)

| Step | Time (s) | % of Total | Location |
|------|----------|-----------|----------|
| Background (RECFAST) | 1.28 | 0.4% | CPU (numpy + scipy ODE) |
| Perturbations (500 k, 10 workers) | 358.32 | 99.1% | CPU (scipy Radau) |
| LOS Bessel functions | 0.69 | 0.2% | **GPU** (MLX Chebyshev) |
| LOS source accumulation | 1.17 | 0.3% | CPU (numpy) |
| C_l integral (TT+EE+TE) | 0.002 | 0.0% | CPU (numpy) |
| CubicSpline derivative | 0.045 | 0.0% | CPU (scipy) |
| **TOTAL** | **361.5** | **100%** | |

**GPU utilization: 0.19% of total time.** The pipeline is entirely CPU-bound.

### Background Breakdown

| Component | Time |
|-----------|------|
| Friedmann solve (numpy) | 1.8 ms |
| Interpolator setup (scipy) | 15.1 ms |
| RECFAST ODE (15,882 RHS evals) | ~1.26 s |
| Reionization iteration | ~3 ms |

Alternative recombination methods:
- `tanh`: 0.064 s (64x faster, but ~5% less accurate x_e)
- `peebles`: 0.57 s (2x faster, ~0.1% accuracy)
- `recfast`: 1.30 s (reference)
- `pickle cache load`: 0.004 s (320x faster, exact)

---

## 2. Per-k-mode ODE Profile

### Full range (tau_init -> 13,891 Mpc)

| k (Mpc^-1) | Time (s) | RHS evals | Steps | ms/eval | RHS/step |
|-------------|----------|-----------|-------|---------|----------|
| 0.0003 | 0.17 | 3,846 | 113 | 0.033 | 34.0 |
| 0.001 | 0.21 | 4,625 | 193 | 0.033 | 24.0 |
| 0.005 | 0.52 | 11,131 | 755 | 0.033 | 14.7 |
| 0.01 | 0.78 | 17,492 | 1,501 | 0.032 | 11.7 |
| 0.05 | 2.83 | 59,969 | 7,269 | 0.031 | 8.2 |
| 0.10 | 5.27 | 110,443 | 14,446 | 0.031 | 7.6 |
| 0.20 | 9.47 | 202,856 | 27,831 | 0.030 | 7.3 |
| 0.35 | 14.83 | 317,059 | 44,081 | 0.030 | 7.2 |

### Early-stop (tau_init -> 400 Mpc, just past recombination)

| k (Mpc^-1) | Time (s) | RHS evals | Steps | Fraction saved |
|-------------|----------|-----------|-------|----------------|
| 0.0003 | 0.14 | 2,917 | 52 | 18% |
| 0.001 | 0.15 | 3,540 | 68 | 29% |
| 0.01 | 0.33 | 7,913 | 157 | 58% |
| 0.05 | 0.52 | 11,660 | 409 | 82% |
| 0.10 | 0.68 | 14,304 | 759 | 87% |
| 0.20 | 0.79 | 16,255 | 1,205 | 92% |
| 0.35 | 1.02 | 20,743 | 1,770 | 93% |

**Key finding**: High-k modes spend 80-94% of their time on late free-streaming
(tau > 400 Mpc), which is only needed for the late ISW effect at low ell.
The hybrid solver already exploits this.

---

## 3. RHS Function Micro-Profile (k=0.1, tau=280 Mpc)

| Component | Time (us) | % of RHS |
|-----------|-----------|----------|
| 4x np.interp (50k grid) | 3.6 | 13% |
| State extraction | 0.9 | 3% |
| Metric diagnose (h', eta') | 1.9 | 7% |
| Photon hierarchy (l=0..25, for-loop) | 15.1 | 54% |
| Neutrino hierarchy (l=0..25, for-loop) | 11.0 | 39% |
| Assembly (np.zeros + assign) | 1.0 | 4% |
| **Full RHS call** | **28.0** | **100%** |

**Dominant cost**: The Python for-loops over multipole hierarchies (l=3..24)
account for ~80% of each RHS call. The pure-Python loop overhead itself is
only ~1 us; the cost is the individual floating-point operations in Python.

### scipy Overhead Analysis (k=0.1, full range)

| Component | Time (s) | % |
|-----------|----------|---|
| Time in RHS function | 3.47 | 64% |
| scipy overhead (Jacobian FD, Newton, LU, memory) | 1.96 | 36% |
| **Total** | **5.42** | **100%** |

### Jacobian Cost

- Radau uses **numerical finite-difference Jacobian**: nvar+1 = 57 RHS calls per update
- Estimated ~924 Jacobian updates for k=0.1 (full range)
- **~48% of all RHS calls are for Jacobian estimation**
- Jacobian sparsity: only 167/3136 = **5.3% fill** (tridiagonal hierarchy + metric coupling)

---

## 4. BDF vs Radau Comparison

| k (Mpc^-1) | Radau (s) | BDF (s) | Speedup | Max rel. diff |
|-------------|-----------|---------|---------|---------------|
| 0.001 | 0.19 | 0.12 | 1.66x | 6.9e-7 |
| 0.01 | 0.76 | 0.36 | 2.10x | 3.4e-8 |
| 0.05 | 2.70 | 1.28 | 2.11x | 4.0e-7 |
| 0.10 | 5.01 | 1.83 | 2.74x | 1.6e-7 |
| 0.20 | 9.14 | 3.36 | 2.72x | 9.9e-7 |
| 0.35 | 14.30 | 5.95 | 2.40x | 3.2e-7 |

**BDF is 2.0-2.7x faster than Radau with <1e-6 relative accuracy difference.**
BDF uses only 2.3 RHS/step vs Radau's 7.6 RHS/step because BDF avoids
the expensive Jacobian re-estimation (it keeps a factored Jacobian and
updates it less frequently).

---

## 5. LOS Integration Profile

| Component | Time (s) |
|-----------|----------|
| Bessel table load (cached) | 0.032 |
| Bessel table build (first run) | ~2-4 |
| Single Bessel eval (372 ells x 500 k) | 0.0078 |
| Bessel loop (163 snapshots) | 0.694 |
| Source accumulation (for-loop over ell) | 0.251 |
| Source accumulation (vectorized, no ell loop) | 0.092 |
| Full LOS integration | 1.86 |
| C_l integral (trapezoid) | 0.37 ms |

**LOS parameters**: 372 ell values, 500 k-modes, 163 tau snapshots = 30.3M Bessel evaluations.

The inner ell loop in the LOS integration could be vectorized for a 2.7x speedup
on that portion (~0.16s saved).

---

## 6. GPU Utilization

| Operation | Currently | On GPU? | % of time |
|-----------|-----------|---------|-----------|
| Background | CPU (numpy + scipy) | No | 0.4% |
| Perturbation ODE | CPU (scipy Radau) | No | 99.1% |
| Bessel evaluation | GPU (MLX Chebyshev) | **Yes** | 0.2% |
| Source accumulation | CPU (numpy) | No | 0.3% |
| C_l integral | CPU (numpy) | No | 0.0% |

**The GPU is idle 99.8% of the time.** The entire bottleneck is scipy's
per-k-mode sequential Radau solver running on CPU.

### Theoretical GPU Potential

If the Boltzmann ODE could be solved entirely on GPU (batched across all 500
k-modes), the theoretical speedup depends on:

1. **RHS vectorization**: Replace 500 separate Python RHS calls with one
   batched MLX operation over all k simultaneously. Each call: 28 us Python
   -> potentially ~0.1 us batched on GPU.

2. **Step synchronization challenge**: Different k-modes need vastly different
   step counts (113 steps for k=0.0003 vs 44,081 for k=0.35). Padding the
   short modes to match the longest wastes ~99.7% of compute for low-k modes.

3. **Implicit solver on GPU**: Radau/BDF require LU factorizations of
   56x56 matrices. Batched LU on GPU for 500 systems is efficient on Apple
   Silicon (MLX `linalg.solve`), but the Newton iteration control flow is
   hard to vectorize.

**Realistic estimate**: A custom BDF on GPU, with adaptive masking (completed
k-modes stop contributing), could achieve ~15-40s total for the perturbation
phase.

---

## 7. Multiprocessing Efficiency

| Metric | Value |
|--------|-------|
| CPU cores | 10 |
| Workers | 10 |
| Estimated serial time | ~3,580s (extrapolated) |
| Measured parallel time | 358s |
| Parallel efficiency | ~100% (10 cores) |
| Overhead factor | 1.0x |

The parallel efficiency is excellent. The bottleneck is simply that 500 k-modes
with an average of ~7s each (dominated by high-k) totals ~3,500s serial work.

### k-mode Time Distribution (serial estimate)

| k range | N modes | Serial time | % |
|---------|---------|-------------|---|
| [0, 0.001) | 86 | 42s | 3% |
| [0.001, 0.01) | 162 | 101s | 8% |
| [0.01, 0.05) | 114 | 172s | 13% |
| [0.05, 0.1) | 49 | 173s | 13% |
| [0.1, 0.2) | 49 | 323s | 25% |
| [0.2, 0.35] | 39 | 460s | 36% |

**The top 88 modes (k > 0.1) consume 61% of total compute time.**

---

## 8. Optimization Opportunities (Ranked by Speedup/Effort)

| # | Optimization | Current | Target | Speedup | Effort | Dependencies |
|---|-------------|---------|--------|---------|--------|-------------|
| **1** | **Switch Radau -> BDF** | 358s | ~145s | **2.5x** | **Low** (1-line change) | None |
| **2** | **Hybrid early-stop** (high-k to tau=400) | 358s | ~50s | **7x** | **Medium** (solver_hybrid exists) | Verify late ISW |
| **3** | **BDF + hybrid combined** | 358s | **~19s** | **~19x** | Medium | #1 + #2 |
| **4** | Analytical Jacobian (sparse, 5.3% fill) | 358s | ~300s | 1.2x | High (derive J analytically) | Fragile |
| **5** | Background caching (pickle) | 1.3s bg | 0.004s bg | 320x bg | Low | Invalidation logic |
| **6** | Vectorize LOS ell loop | 1.86s LOS | 1.6s LOS | 1.2x LOS | Low | None |
| **7** | Reduce N_k (500 -> 300) | 358s | ~215s | 1.7x | Low | Verify RMS < 8% |
| **8** | GPU batched BDF (custom) | 358s | 15-40s | **9-24x** | **Very High** | MLX batched LU, masking |
| **9** | Reduce hierarchy (lg=25 -> 15) | 0.68s/k | ~0.46s/k | 1.5x/k | Low | Verify accuracy |
| **10** | Green's function precompute | 358s first | 0.05s repeat | **7000x repeat** | Done (solver_greens.py) | First run still slow |

### Recommended Priority Order

**Immediate wins (< 1 hour implementation)**:

1. **Switch to BDF**: Change `method='Radau'` to `method='BDF'` in `solve_ivp`.
   Measured 2.0-2.7x faster with <1e-6 accuracy difference. BDF uses only
   2.3 RHS/step vs Radau's 7.6, because it avoids expensive per-step Jacobian
   finite-difference estimation (57 RHS calls each). Validated in parallel
   pipeline (N_k=50, zero failures, identical C_l peaks).
   **Expected: 358s -> ~145s** (516s serial / 10 workers + overhead).

2. **Combine with hybrid early-stop**: The existing `solver_hybrid.py` already
   implements this. With BDF:
   - Phase 1 (all 500 modes to tau=400): ~113s serial -> 11s parallel
   - Phase 2 (248 low-k modes, full range): ~41s serial -> 4s parallel
   - Background: 1.3s, LOS: 1.9s
   - **Expected total: ~19s** (19x speedup over current 361s).

**Medium-term (1-2 days)**:

3. **GPU batched BDF**: Write a custom BDF solver in MLX that processes all
   500 k-modes simultaneously. The 56x56 linear systems can be batched as
   (500, 56, 56) tensor operations. Adaptive masking handles different step
   counts per k. This is architecturally complex but would bring the
   perturbation solve to ~15-30s. Combined with the existing GPU Bessel LOS,
   total pipeline: **~20-35s.**

**Already available**:

4. **Green's function cache** (solver_greens.py): For repeated C_l evaluations
   with different primordial spectra (A_s, n_s), the source function cache
   gives instant (~50ms) C_l after one-time precomputation. Ideal for MCMC.

---

## 9. Detailed Cost Model

### Where the 358s Goes (Perturbation Phase)

```
358s total perturbation time (10 workers)
  = ~3580s serial work / 10 cores

Per k-mode (average ~7.2s serial):
  64% = RHS evaluation (4.6s)
    13% = np.interp (0.6s)
    54% = photon hierarchy for-loop (2.5s)
    39% = neutrino hierarchy for-loop (1.8s)
    ~4% = metric + assembly (0.2s)
  36% = scipy overhead (2.6s)
    48% of RHS calls = Jacobian finite differences
    ~15% = LU factorization (56x56)
    ~10% = Newton iteration control
    ~11% = memory allocation + Python overhead
```

### Why BDF Wins

- Radau: 5th-order, 3-stage implicit. 7.6 RHS/step.
  Jacobian re-estimated ~every 15 steps via nvar+1=57 finite-difference calls.
- BDF: up to 5th-order, multistep. 2.3 RHS/step.
  Jacobian updated less frequently, no per-step stage overhead.
- For this moderately-stiff system, BDF's lower per-step cost dominates
  despite taking slightly more steps.

### Why Hybrid Wins

- 252 modes with k >= 0.01 save 80-94% of their time (no late free-streaming)
- 248 modes with k < 0.01 still need full integration (late ISW matters)
- Net: ~91% of high-k serial time eliminated

---

## 10. Existing Solver Comparison

| Solver | Time | Accuracy vs CLASS | Notes |
|--------|------|-------------------|-------|
| Full scipy Radau | 361.5s | Reference | This audit |
| Hybrid (solver_hybrid.py) | 35-60s | 7% RMS | Early-stop for high-k |
| Magnus GPU (solver_magnus.py) | 7.6s | 58% | Matrix exponential, inaccurate |
| Green's function (precomputed) | 1.1ms | Exact match | Requires ~300s precompute |
| **Full scipy BDF** (projected) | **~145s** | **~identical** | **1-line change** |
| **BDF + hybrid** (projected) | **~19s** | **~7% RMS** | **Recommended next step** |
