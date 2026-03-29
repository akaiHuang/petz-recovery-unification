# mlx-class Technical Roadmap: Surpassing CLASS

**Goal**: Build a GPU-native CMB Boltzmann solver that exceeds CLASS in both speed and accuracy.

**Current state**: mlx-class v0.2.0 solves the tight-coupling approximation (TCA) with 6 variables (Phi, delta_b, delta_c, v_c, Theta_0, Theta_1) using batched RK4/IMEX on Apple MLX. Recombination uses a tanh fit. No free-streaming, no polarization, no neutrino hierarchy, no line-of-sight integration.

**Architecture principle**: CLASS was written for single-core CPUs in 2011. Its entire design reflects that constraint: aggressive truncation of multipole hierarchies, analytic approximations for pre-computed integrals, adaptive step size with serial ODE integration. A GPU solver can remove all of these compromises by solving the full system in parallel.

---

## Table of Contents

1. [Full Photon Hierarchy (l_max >> 25)](#1-full-photon-hierarchy)
2. [Second-Order Perturbation Theory](#2-second-order-perturbation-theory)
3. [Precision Recombination (HyRec)](#3-precision-recombination-hyrec)
4. [Mixed-Precision Strategy](#4-mixed-precision-strategy)
5. [Full Line-of-Sight Integration](#5-full-line-of-sight-integration)
6. [Neutrino Boltzmann Hierarchy](#6-neutrino-boltzmann-hierarchy)
7. [Novel Physics Support (Khronon, MG, EDE)](#7-novel-physics-support)
8. [Differentiable Pipeline](#8-differentiable-pipeline)
9. [Implementation Phases](#9-implementation-phases)

---

## 1. Full Photon Hierarchy

### What CLASS does

CLASS truncates the photon Boltzmann hierarchy at l_max_photon ~ 12-25 (depending on k and epoch). Beyond the truncation multipole l_t, it uses the approximation:

```
Theta_{l_t+1} = (2*l_t + 1) / (k * tau) * Theta_{l_t} - Theta_{l_t-1}
```

This is the free-streaming closure relation (e.g., Ma & Bertschinger 1995, Eq. 58). It is exact only in the free-streaming limit (kappa_dot -> 0) and introduces errors of order (kappa_dot * tau)^2 / l_t^2 near recombination where the transition from tight coupling to free streaming is not instantaneous.

Additionally, CLASS uses the tight-coupling approximation (TCA) before tau_switch ~ 0.85 * tau_rec with only Theta_0 and Theta_1 (collapsing the hierarchy to 2 moments), then switches to the full hierarchy. This switch introduces a matching error at percent level.

**CLASS's serial bottleneck**: Each k-mode is integrated sequentially. For N_k = 500, l_max = 25, that is 500 serial ODE integrations, each with ~25 photon multipoles + other species. Total variables per k: ~80. Serial wall time: ~10s.

### What GPU can do better

On GPU, all k-modes are batched. The cost of l_max = 100 vs l_max = 25 is:

- Memory: N_k * l_max * sizeof(float32) = 500 * 100 * 4 bytes = 200 KB (trivial)
- Compute: each RK4 step touches N_k * l_max elements. For l_max = 100, this is 4x more FLOPs per step than l_max = 25, but fully parallelized across k-modes and multipoles. On M1 Pro (2 TFLOPS fp32), a single RK4 step over 500 k-modes * 100 multipoles takes < 1 microsecond.
- Total: ~3000 time steps * 4 RK4 evaluations * (500 * 100 * 10 ops) = 6 GFLOP. At 2 TFLOPS, that is 3 ms. Add overhead: realistic estimate is 0.5-2s total.

**Key insight**: With l_max = 100, no closure approximation is needed at all. The hierarchy naturally decays because Theta_l ~ (k*tau)^l / (2l+1)!! for l >> k*tau. For k_max = 0.35 Mpc^-1, k*tau_rec ~ 100, so l_max = 100 captures the full dynamics without truncation artifacts.

### Accuracy improvement

| l_max | Closure error (l ~ 2000) | Notes |
|-------|--------------------------|-------|
| 12    | ~2-5%                    | CLASS default for TCA regime |
| 25    | ~0.5-1%                  | CLASS default for free-streaming |
| 50    | ~0.01-0.05%              | Sufficient for Planck |
| 100   | < 0.001%                 | Effectively exact for k < 0.35 |
| 200   | Machine precision         | Needed only for k > 1.0 |

The dominant improvement is at high l (l > 1500) and for polarization, where the E-mode signal is generated from Theta_2 (quadrupole) at recombination. Truncation errors in Theta_2 propagate directly into EE and TE spectra.

For TT at l < 2500: improvement from CLASS's ~0.1% to ~0.001%.
For EE at l < 3000: improvement from CLASS's ~0.5% to ~0.01%.

### Implementation

```
# Current state vector per k:
y[k, 0:6] = [Phi, delta_b, delta_c, v_c, Theta_0, Theta_1]

# Target state vector per k:
y[k, 0:4]        = [Phi, delta_b, delta_c, v_c]     # metric + CDM
y[k, 4:4+l_max]  = [Theta_0, ..., Theta_{l_max-1}]  # photon T
y[k, 4+l_max:4+2*l_max] = [E_0, ..., E_{l_max-1}]   # photon E-pol
y[k, 4+2*l_max:...] = [N_0, ..., N_{l_nu_max-1}]    # neutrinos
```

The photon hierarchy equation:
```
Theta_l' = k/(2l+1) * [l * Theta_{l-1} - (l+1) * Theta_{l+1}]
           - kappa_dot * [Theta_l - source_l]
```

This is a tridiagonal recurrence in l. On GPU, the entire l-vector for all k-modes can be updated simultaneously using a single matrix-vector product:

```python
# L is (l_max, l_max) tridiagonal, k-dependent
# Theta is (N_k, l_max)
# This is a batched tridiagonal matvec: perfectly suited for GPU
dTheta = k_arr[:, None] * (L @ Theta.T).T - kappa_dot * (Theta - S)
```

### Difficulty: MEDIUM

The physics is well-understood (Ma & Bertschinger 1995). The challenge is:
1. Matching TCA to full hierarchy (or running full hierarchy from the start)
2. Ensuring numerical stability with the tridiagonal recurrence at high l
3. Efficient batched tridiagonal operations in MLX

---

## 2. Second-Order Perturbation Theory

### What CLASS does

CLASS computes first-order perturbation theory only. All C_l are linear in the primordial power spectrum P_R(k). Non-linear corrections are added post-hoc:
- **Lensing**: Applied as a convolution on C_l using the lensing potential C_l^{phi phi}, computed from the Weyl potential integral. This is first-order in the lensing potential.
- **No second-order Boltzmann**: The coupled second-order Einstein-Boltzmann system is not solved.

Second-order effects are handled by separate codes (SONG, CosmoLib2) that are extremely slow (hours per spectrum) and limited in accuracy validation.

### What GPU can do better

Second-order perturbation theory introduces mode-coupling: the evolution of a mode (k1, k2) depends on products of first-order modes. The number of coupled modes scales as N_k^2. On CPU this is prohibitive; on GPU it becomes feasible.

Key second-order terms:

1. **Lensing of CMB (already in CLASS post-hoc)**:
   - C_l^{TT,lensed} = C_l^{TT} + integral over l' of C_{l'}^{phi phi} * coupling kernel
   - GPU advantage: the convolution over l' is a matrix multiplication, trivially parallel
   - Accuracy: 1-5% effect at l > 1000. CLASS's first-order lensing is adequate; second-order lensing-lensing is < 0.1% and not needed for current data

2. **ISW-lensing correlation**:
   - The integrated Sachs-Wolfe effect is correlated with the lensing potential
   - This creates a non-Gaussian bispectrum that biases the lensed power spectrum
   - Effect: ~0.5% at l ~ 100 in TT
   - GPU cost: requires storing Phi(k, tau) on a tau grid, then computing the correlation integral. Memory: N_k * N_tau * 4 bytes = 500 * 3000 * 4 = 6 MB. Trivial.

3. **Rees-Sciama effect** (second-order time derivative of Phi):
   - Phi'' from non-linear growth
   - Effect: < 0.1% for l < 3000
   - GPU cost: same as ISW-lensing

4. **Sunyaev-Zel'dovich (thermal and kinetic)**:
   - Requires halo model or N-body calibration
   - Not computable from linear Boltzmann alone
   - Out of scope for Boltzmann solver

5. **Spectral distortion from recombination**:
   - Second-order effect on x_e(z) from stimulated recombination
   - Changes the visibility function at 0.1% level
   - Can be included in HyRec (see Section 3)

### Practical strategy

Do NOT solve the full second-order Boltzmann system. Instead:

1. **Full first-order to l_max = 100+** (Section 1) -- this removes all truncation errors
2. **Exact lensing convolution** on GPU (matrix multiply, not the flat-sky approximation CLASS uses)
3. **ISW-lensing correlation** as a perturbative correction
4. **Store Phi(k, tau) trajectory** for post-processing second-order terms

This hybrid approach gets 99.9% of the accuracy of a full second-order code at 0.1% of the computational cost.

### Accuracy improvement

| Effect | CLASS | mlx-class target | Improvement |
|--------|-------|-------------------|-------------|
| Lensing (TT, l>1000) | ~0.1% (first-order flat-sky) | ~0.01% (exact curved-sky) | 10x |
| ISW-lensing (TT, l~100) | Not included | ~0.05% | New |
| Rees-Sciama | Not included | ~0.01% | New |
| Net (TT, all l) | ~0.2% | ~0.02% | 10x |

### Difficulty: HARD for full second-order; MEDIUM for the hybrid approach above

---

## 3. Precision Recombination (HyRec)

### What CLASS does

CLASS uses HyRec (by Yacine Ali-Haimoud) as its default recombination code. HyRec solves a multi-level hydrogen atom (effective 3-level: ground state, first excited, continuum) with:
- Lyman-alpha radiative transfer
- Two-photon decay from 2s -> 1s
- Stimulated recombination and photoionization
- Helium recombination (effective rate)
- Matter temperature evolution (Compton heating/cooling)

HyRec achieves ~0.1% accuracy in x_e(z), which translates to ~0.1% in C_l.

RECFAST (the older code) uses fitting functions and is accurate to ~1% in x_e, giving ~0.5% in C_l.

### What mlx-class currently does

A tanh fit for x_e(z) centered at z = 1089 with width Delta_z = 80. This is accurate to ~5% in x_e and ~2-5% in C_l. It captures the location of the visibility function peak but not its width or asymmetry.

### What GPU can do better

HyRec's multi-level atom is an ODE system in z (or a) with ~5-10 variables:
- x_e (free electron fraction)
- T_m (matter temperature)
- x_{2s}, x_{2p} (excited hydrogen populations)
- Lyman-alpha distortion parameter

This is a STIFF system (Lyman-alpha rate >> Hubble rate) requiring implicit integration. On CPU, HyRec takes ~0.5s.

On GPU, we can:

1. **Solve the multi-level atom for many parameter sets simultaneously**: For MCMC with 10^4 samples, batch 10^4 HyRec integrations on GPU. Each is ~10 variables * 1000 time steps = 10^4 operations. Batch of 10^4: 10^8 operations, < 1ms on GPU.

2. **Extend to more levels**: Instead of the effective 3-level atom, solve the full hydrogen atom with n_max = 100 shells. This requires tracking ~100 population levels. CPU cost: ~minutes. GPU cost: batched matrix exponentials for the rate matrix at each time step. For a 100x100 matrix, eigendecomposition takes ~10^6 FLOPs. Over 1000 time steps: 10^9 FLOPs = 0.5 ms on GPU.

3. **Include higher-order effects**:
   - Raman scattering (n=2 -> n=2 transitions)
   - Stimulated two-photon decay
   - Collisional processes
   - These change x_e by ~0.05% and C_l by ~0.02%

### Implementation plan

```
Phase 1: Port HyRec's effective 3-level atom to MLX
  - Replace tanh fit with actual ODE integration
  - Use implicit midpoint rule (L-stable for stiff Lyman-alpha)
  - Expected accuracy: match HyRec to 0.1%
  - GPU benefit: batched over cosmological parameter sets

Phase 2: Extended multi-level atom (n_max = 30-100)
  - Full rate matrix R_{ij} for hydrogen levels
  - GPU-accelerated matrix exponential at each time step
  - Expected accuracy: 0.01% in x_e, 0.005% in C_l

Phase 3: Helium recombination (He I, He II)
  - Separate multi-level atom for He I (intercombination lines)
  - He II recombination at z ~ 6000 (affects photon diffusion)
```

### Accuracy improvement

| Recombination model | x_e accuracy | C_l accuracy (l < 3000) |
|---------------------|-------------|------------------------|
| Tanh fit (current)  | ~5%         | ~2-5%                  |
| RECFAST             | ~1%         | ~0.5%                  |
| HyRec (3-level)     | ~0.1%       | ~0.1%                  |
| Multi-level (n=100) | ~0.01%      | ~0.005%                |

### Difficulty: MEDIUM (Phase 1), HARD (Phase 2-3)

---

## 4. Mixed-Precision Strategy

### What CLASS does

CLASS uses double precision (float64) throughout. This is safe but wastes compute on parts of the calculation where float32 or even float16 would suffice.

### GPU mixed-precision strategy

The key insight: different parts of the Boltzmann system have different precision requirements.

**Analysis of precision requirements by variable:**

| Variable/Equation | Dynamic range | Precision needed | Rationale |
|-------------------|---------------|------------------|-----------|
| Theta_l (photon moments) | ~10^{-5} to 1 | float32 | Oscillatory, well-conditioned |
| Phi (gravitational potential) | ~10^{-5} | float32 or float64 | Poisson equation is stiff: lambda = -k^2/(3*calH). For k = 0.35, calH ~ 0.01: lambda ~ -4. Moderate stiffness. |
| delta_b, delta_c (density) | ~1 to 10^3 | float32 | Grows as a(t), well-conditioned |
| v_b, v_c (velocity) | ~0.1 to 10 | float32 | Oscillatory, moderate |
| x_e (recombination) | 10^{-4} to 1 | float64 | Spans 4 decades, stiff ODE |
| Poisson source sum | Sum of 3 terms with cancellation | float64 | Omega_r/a^2 * delta_gamma + Omega_b/a * delta_b + Omega_c/a * delta_c: the three terms can partially cancel |
| Bessel j_l(x) for x >> l | ~10^{-20} to 1 | float64 | Oscillatory with exponential decay envelope |
| C_l integration | Sum of ~500 terms | float32 | Kahan summation or float64 accumulator |

**Proposed mixed-precision architecture:**

```python
# Main ODE integration: float32
# 90% of all compute, fully GPU-accelerated
y = mx.array(y0, dtype=mx.float32)           # state vector
k_arr = mx.array(k_arr, dtype=mx.float32)    # k-modes

# Poisson equation source: float64 accumulation
# Critical: Omega_r/a^2 * delta_gamma cancels against Omega_m/a * delta_m near equality
# Use compensated summation or promote to float64 for this one line
Phi_source = mx.float64(Omega_r_a2 * delta_gamma) + mx.float64(Omega_b_a * delta_b) + ...
Phi_dot = mx.float32(Phi_source / (3 * calH) - ...)

# Recombination ODE: float64
# Only ~10 variables, N_tau ~ 1000 steps, negligible cost
x_e = np.float64(...)  # numpy on CPU is fine

# Bessel functions: float64 computation, float32 storage
# scipy computes in float64; we cast to float32 for GPU storage
# For x >> l, j_l(x) ~ sin(x - l*pi/2) / x, which is fine in float32
# For x << l, j_l(x) ~ (x/2)^l / Gamma(l+3/2), underflows in float32 for l > 80
# Solution: use log-space for j_l^2 when x < l

# C_l summation: float32 with Kahan compensated summation
# Or: float64 accumulator (only N_ell * N_k ~ 300*500 = 150K adds)
```

**MLX-specific considerations:**
- MLX on Apple Silicon supports float16, float32, bfloat16 natively at full throughput
- float64 is supported but at reduced throughput (typically 1/32 of float32 on Apple GPU)
- Strategy: keep float64 on CPU (numpy) for the few operations that need it; everything else in float32 on GPU

### Accuracy improvement

| Precision strategy | Speed | Accuracy |
|--------------------|-------|----------|
| All float64 (CLASS) | 1x (CPU baseline) | ~10^{-14} (overkill) |
| All float32 (current mlx-class) | 6-40x | ~10^{-6} (may lose precision in Poisson) |
| Mixed (proposed) | 5-35x | ~10^{-8} (matches CLASS where it matters) |

The accuracy improvement is not about surpassing CLASS's numerical precision (which is already overkill) but about matching it while maintaining GPU speed. The real accuracy gains come from Sections 1-3 and 5-6.

### Difficulty: EASY (Phase 1 Poisson fix), MEDIUM (full mixed-precision pipeline)

---

## 5. Full Line-of-Sight Integration

### What CLASS does

CLASS computes C_l using the line-of-sight (LOS) integration method (Seljak & Zaldarriaga 1996):

```
Delta_l(k) = integral_0^{tau_0} d tau  S(k, tau) * j_l[k(tau_0 - tau)]
```

where the source function S(k, tau) contains:

1. **Sachs-Wolfe (SW)**: g(tau) * [Theta_0 + Phi + Psi] -- visibility function times the monopole + potential
2. **Doppler**: g(tau) * v_b * d/d(k*tau) j_l -- velocity projection
3. **Integrated Sachs-Wolfe (ISW)**: e^{-kappa} * [Phi' + Psi'] -- potential decay along the line of sight
4. **Polarization source**: g(tau) * Pi * [j_l'' and j_l terms] -- from the photon quadrupole Theta_2

CLASS pre-computes S(k, tau) on a grid of ~3000 time steps and ~500 k-modes, then performs the tau-integral for each (l, k) pair using cubic spline interpolation of j_l.

### What mlx-class currently does

Instantaneous recombination approximation: evaluate the source at tau_rec (or 0.85*tau_rec for TCA) and compute:

```
Delta_l(k) ~ S(k, tau_rec) * j_l(k * D_A)
```

This misses:
- The finite width of the visibility function (Delta_z ~ 80, or Delta_tau ~ 50 Mpc)
- ISW contribution (important at l < 30 and l > 2000)
- Reionization bump at l < 10
- Doppler projection (j_l' vs j_l)

The result: peak positions are off by ~30% and the overall shape is qualitatively correct but quantitatively wrong.

### What GPU can do better

The LOS integral is embarrassingly parallel over both l and k:

```
Delta_l(k) = sum_i  S(k, tau_i) * j_l[k * (tau_0 - tau_i)] * Delta_tau_i
```

This is a matrix multiplication:
- S is (N_k, N_tau) -- source function on the grid
- j_l is (N_ell, N_k, N_tau) -- Bessel functions on the grid
- Product: (N_ell, N_k) -- transfer function

On GPU, this is a batched einsum:
```python
# S: (N_k, N_tau), j_l: (N_ell, N_k, N_tau)
Delta = mx.sum(S[None, :, :] * jl[:, :, :] * dtau[None, None, :], axis=2)
# Output: (N_ell, N_k)
```

Memory: N_ell * N_k * N_tau * 4 bytes = 300 * 500 * 3000 * 4 = 1.8 GB. This is too large for a single allocation on most GPUs.

**Chunked approach**:
```python
# Process ell in chunks of 10-30
for ell_chunk in chunks(ell_values, chunk_size=20):
    jl_chunk = compute_bessel(ell_chunk, k_arr, tau_grid)  # (20, N_k, N_tau)
    Delta_chunk = einsum(S, jl_chunk, dtau)                # (20, N_k)
```

Memory per chunk: 20 * 500 * 3000 * 4 = 120 MB. Feasible.

Alternatively, use the Limber approximation for l > 100:
```
j_l(x) ~ sqrt(pi / (2l+1)) * delta(x - l - 1/2)
```
This reduces the integral to a single evaluation S(k, tau) at k*(tau_0 - tau) = l + 1/2, removing the tau integral entirely for high l.

**Optimal strategy**:
- l < 50: Full LOS integral with exact j_l (N_tau ~ 3000 quadrature points)
- l = 50-200: Full LOS integral with N_tau ~ 1000 (j_l oscillates slowly)
- l > 200: Extended Limber approximation with ~100 correction terms

### Quadrature points needed

The integrand S(k, tau) * j_l[k*(tau_0 - tau)] oscillates with frequency ~k in tau. The envelope (visibility function) has width ~50 Mpc centered at tau_rec ~ 280 Mpc.

Required quadrature:
- **Near recombination** (tau_rec - 50 to tau_rec + 50 Mpc): 200-500 points to resolve the visibility function and source oscillations
- **ISW epoch** (tau = 280 to 14000 Mpc): 100-200 points (slowly varying integrand)
- **Reionization** (tau ~ 4000-5000 Mpc, z ~ 6-8): 50-100 points
- **Total**: ~500-800 points (CLASS uses ~3000 for safety; we can use fewer with Gauss-Legendre)

With GPU, even 3000 points costs:
- Per (l, k): 3000 multiplies + 3000 adds = 6000 FLOPs
- Total: 300 * 500 * 6000 = 9 * 10^8 = 0.9 GFLOP
- At 2 TFLOPS: 0.5 ms

**The bottleneck is computing j_l, not the integral itself.** Bessel function evaluation costs ~50-100 FLOPs each. For 300 * 500 * 3000 = 4.5 * 10^8 evaluations: ~45 GFLOP = 22 ms.

Pre-compute and cache j_l on disk: 4.5 * 10^8 * 4 bytes = 1.8 GB. Too large. Better: compute j_l on-the-fly using the recurrence j_{l+1}(x) = (2l+1)/x * j_l(x) - j_{l-1}(x), which is 2 FLOPs per (l, x) pair.

### Accuracy improvement

| Integration method | Peak position error | Amplitude error | ISW |
|--------------------|--------------------|-----------------|----|
| Instantaneous (current) | ~30% | ~50% | No |
| LOS with 500 points | < 0.1% | < 0.5% | Yes |
| LOS with 3000 points | < 0.01% | < 0.1% | Yes |
| CLASS LOS | < 0.01% | < 0.1% | Yes |

This is the single biggest accuracy improvement available. It takes mlx-class from "qualitative" to "precision cosmology".

### Difficulty: HARD (this is the critical path item)

---

## 6. Neutrino Boltzmann Hierarchy

### What CLASS does

CLASS solves the massless neutrino hierarchy:
```
N_l' = k/(2l+1) * [l * N_{l-1} - (l+1) * N_{l+1}]
```
(no collision term, since neutrinos decouple at z ~ 10^{10}).

For massless neutrinos (3.046 effective species), CLASS uses:
- l_max_nu ~ 17 for non-relativistic regime
- Fluid approximation (ncdm module) for massive neutrinos with only (delta, theta, sigma) = 3 variables
- The ultra-relativistic fluid approximation (UFA) switches on when k*tau >> l_max_nu

For massive neutrinos (sum m_nu ~ 0.06 eV), CLASS uses the ncdm module with:
- Momentum-dependent distribution function f(q, tau) on a q-grid (~15 momentum bins)
- Each q-bin has its own Boltzmann hierarchy with l_max ~ 17
- Total: 15 * 17 = 255 variables per k-mode per massive neutrino species

### What GPU can do better

1. **More multipoles**: Push l_max_nu to 50-100, eliminating the truncation/UFA approximation entirely. Cost increase: 3-6x in neutrino sector. On GPU with batched operations, negligible wall-time increase.

2. **More momentum bins**: For massive neutrinos, use 50-100 q-bins instead of 15. This reduces the interpolation error in the momentum integral from ~0.5% to ~0.01%.

3. **Phase-space method**: Instead of the multipole expansion, directly solve the phase-space distribution f(x, p, hat{n}) on a discrete angle grid (N_angles ~ 100-200 directions). This avoids the truncation problem entirely but requires more memory. For N_k = 500, N_angles = 200, N_q = 50: state vector size = 500 * 200 * 50 * 4 = 20 MB. Feasible on GPU.

### Neutrino multipoles needed

The free-streaming length of massless neutrinos at recombination is:
```
lambda_fs = tau_rec ~ 280 Mpc
k_fs * tau_rec = l_max_needed
```

For k_max = 0.35 Mpc^{-1}: k * tau_rec ~ 100. So l_max_nu ~ 100 is needed for exact treatment.

With l_max_nu = 17 (CLASS), the truncation error at k = 0.35 is:
- N_{17} / N_0 ~ (k*tau/17)^{17} / (2*17+1)!! which is very small for k*tau < 17
- But for k*tau >> 17, the UFA approximation kicks in, with accuracy ~(1/k*tau)^2 ~ 1%

### GPU cost

| Configuration | Variables per k | Total state (500 k) | Wall time |
|---------------|----------------|---------------------|-----------|
| CLASS default (l_nu=17) | 17 per species | 500 * 17 * 4 = 34K | 0.5s (CPU) |
| l_nu = 50 | 50 per species | 500 * 50 * 4 = 100K | ~0.01s (GPU) |
| l_nu = 100 | 100 per species | 500 * 100 * 4 = 200K | ~0.02s (GPU) |
| Phase-space (200 angles) | 200 per species | 500 * 200 * 4 = 400K | ~0.05s (GPU) |
| Massive (50 q-bins, l=50) | 2500 per species | 500 * 2500 * 4 = 5M | ~0.2s (GPU) |

Memory is not a constraint. The compute cost scales linearly with l_max_nu * N_k, and on GPU this is fully parallel.

### Accuracy improvement

| Neutrino treatment | Error in C_l (l < 2500) | Phase shift error |
|--------------------|------------------------|-------------------|
| Background only (current mlx-class) | ~5% at l > 500 | ~2% |
| CLASS (l_nu=17 + UFA) | ~0.2% | ~0.05% |
| Full hierarchy (l_nu=100) | ~0.001% | < 0.001% |
| Phase-space (200 angles) | ~10^{-5}% | < 10^{-5}% |

The neutrino phase shift (change in acoustic peak positions due to neutrino free-streaming) is a key observable for measuring N_eff. Improving from 0.05% to 0.001% enables better constraints on N_eff and neutrino masses.

### Difficulty: MEDIUM

The equations are simple (no collision term). The implementation is straightforward: extend the state vector and add the hierarchy recurrence. The batched GPU architecture already handles this naturally.

---

## 7. Novel Physics Support (Khronon, MG, EDE)

### Current Khronon implementation

mlx-class already supports the ghost-condensation Khronon model with:
- Omega_K = 0.265 replacing Omega_c in the Friedmann equation
- c_s^2 = 0 (dust-like perturbations from K'(Q_0) = 0)
- Khronon correction factor on the analytic transfer function

This is equivalent to CDM at the perturbation level (because c_s^2 = 0 makes Khronon perturbations behave exactly like CDM perturbations).

### Extensible architecture for novel physics

The key architectural decision: separate the **equation module** from the **solver module**.

```python
class PerturbationEquations(ABC):
    """Abstract base class for perturbation equations."""

    @abstractmethod
    def n_variables(self) -> int:
        """Total number of variables per k-mode."""

    @abstractmethod
    def initial_conditions(self, k_arr, bg) -> mx.array:
        """Shape (N_k, n_var). Adiabatic IC by default."""

    @abstractmethod
    def derivative(self, y, k_arr, bg_at_tau) -> mx.array:
        """Shape (N_k, n_var). The RHS of dy/dtau."""

    @abstractmethod
    def source_function(self, y, k_arr, bg_at_tau) -> dict:
        """Return {name: (N_k,) array} for line-of-sight integration."""

    def stiff_eigenvalue(self, k_arr, bg_at_tau) -> mx.array:
        """For IMEX solver: the stiff part of the Jacobian."""
        return mx.zeros(len(k_arr))
```

**Concrete implementations:**

```python
class LCDMEquations(PerturbationEquations):
    """Standard LCDM: photon hierarchy + baryons + CDM + neutrinos."""
    # l_max_photon, l_max_neutrino configurable
    # Polarization optional

class KhrononEquations(PerturbationEquations):
    """Khronon DM: replaces CDM with Khronon field."""
    # Extra variable: pi (Khronon perturbation)
    # DBI kinetic function K(Q) with K'(Q_0) = 0
    # c_s^2 = K''(Q_0) * Q_0 / (K'(Q_0) + 2*K''(Q_0)*Q_0)
    # For ghost condensation: K'(Q_0) = 0 => c_s^2 = 0

class ModifiedGravityEquations(PerturbationEquations):
    """General modified gravity with mu(k,a), Sigma(k,a) parameterization."""
    # Phi = -(1 + mu) * 4piG * rho * delta / k^2
    # (Phi + Psi) = -(1 + Sigma) * 8piG * rho * delta / k^2

class EarlyDarkEnergyEquations(PerturbationEquations):
    """EDE with scalar field phi(tau)."""
    # Additional background ODE for phi(tau)
    # Perturbation: delta_phi, v_phi
    # Oscillating potential V(phi) = m^2 f^2 (1 - cos(phi/f))^n
```

### Khronon-specific features

For the Khronon model (Blanchet & Skordis 2025), the perturbation equation for the Khronon field pi is:

```
pi'' + 2*H*pi' + (c_s^2 * k^2 + m^2) * pi = source(Phi, v_b, ...)
```

With c_s^2 = 0 (ghost condensation), this becomes:
```
pi'' + 2*H*pi' + m^2 * pi = source(Phi, v_b, ...)
```

which is a damped oscillator with mass m = mu*c = H_0. The solution tracks CDM perturbations exactly at scales k >> H_0/c, with deviations at k ~ H_0/c (horizon scale).

GPU advantage: the Khronon equation adds just 2 variables per k-mode (pi and pi'). With batched integration, the cost increase is negligible.

### Novel physics test suite

For each new physics model, automated tests should verify:
1. Correct background evolution (H(z) matches analytic predictions)
2. Correct initial conditions (match the adiabatic mode coefficients)
3. Conservation laws (covariant conservation of stress-energy)
4. Known limits (e.g., Khronon -> CDM when c_s^2 = 0)
5. Stability (no NaN/Inf for the full k-range)

### Difficulty: MEDIUM (architecture), EASY (individual models once architecture exists)

---

## 8. Differentiable Pipeline

### What CLASS does

CLASS is written in C with no automatic differentiation support. Gradient-based inference (e.g., Hamiltonian Monte Carlo) requires either:
- Finite differences: 2 * N_params CLASS evaluations per gradient (N_params ~ 6-20)
- External emulators: neural network trained on CLASS outputs (cosmopower, CosmoPower-JAX)
- No native gradient

### What GPU can do

MLX supports automatic differentiation via `mx.grad()`, `mx.vjp()`, and `mx.jvp()`. However, there are significant caveats:

**What is differentiable in MLX:**
- All element-wise operations (add, mul, exp, sin, cos, ...)
- Matrix operations (matmul, einsum)
- Reductions (sum, mean, max via softmax approximation)
- Control flow with `mx.cond()` (limited)

**What is NOT differentiable:**
- scipy calls (Bessel functions, interpolation, ODE solvers)
- numpy operations
- Integer indexing and sorting
- Discontinuous operations (argmax, where with integer conditions)

**Differentiability audit of the current pipeline:**

| Component | Current impl | Differentiable? | Path to differentiability |
|-----------|-------------|-----------------|---------------------------|
| Background H(z) | numpy | No | Rewrite in MLX: straightforward |
| Recombination x_e(z) | numpy tanh | No | MLX tanh: trivial. HyRec ODE: need MLX ODE solver |
| Visibility function | numpy/scipy | No | MLX computation: straightforward |
| Perturbation ODE | MLX RK4 | **Partially** | The for-loop over tau steps breaks the graph. Need `mx.while_loop` or unrolled computation |
| Bessel functions | scipy | No | **Critical bottleneck**. Options: (a) MLX recurrence, (b) Chebyshev approximation, (c) asymptotic series |
| C_l integration | MLX sum | Yes | Already differentiable |

### Bessel function challenge

The biggest obstacle to full differentiability is the spherical Bessel functions j_l(x). Options:

1. **Backward recurrence in MLX**: j_{l-1}(x) = (2l+1)/x * j_l(x) + j_{l+1}(x). Start from l_max with Miller's algorithm. This is differentiable if implemented in MLX, but backward recurrence through many l values may have gradient stability issues.

2. **Chebyshev polynomial approximation**: Pre-compute Chebyshev coefficients for j_l(x) on intervals. Evaluation is a polynomial, trivially differentiable. Accuracy: ~10^{-10} with degree-20 Chebyshev on each interval. Storage: 300 ell * 20 coefficients * 100 intervals = 600K floats = 2.4 MB.

3. **Asymptotic expansion**: For x >> l, j_l(x) ~ sin(x - l*pi/2) / x. For x << l, j_l(x) ~ (x/2l)^l. Both are differentiable. Stitch together with a smooth transition.

4. **Pre-computed lookup table with interpolation**: Store j_l on a dense x-grid, use MLX cubic interpolation. Differentiable through the interpolation. Accuracy limited by grid spacing.

**Recommendation**: Option 2 (Chebyshev) for production, Option 3 (asymptotic) for quick prototype.

### ODE solver differentiability

The RK4 loop:
```python
for i in range(N_steps):
    y = rk4_step(y, ...)
```

is not natively captured by MLX's computation graph because each iteration creates new graph nodes. For N_steps = 3000, this creates a graph with 3000 * 4 * N_var * N_k nodes, which may exceed memory for backward pass.

Solutions:
1. **Checkpointing**: Store y at every 100th step, recompute intermediate steps during backward pass. Memory: 30 checkpoints * N_k * N_var * 4 bytes ~ 100 KB. Compute: 2x forward pass.

2. **Adjoint method** (Neural ODE approach): Solve the adjoint ODE backward in time. Memory: O(N_var), independent of N_steps. Compute: 1x forward + 1x backward ODE integration. This is the gold standard for differentiable ODE solvers.

3. **Implicit differentiation**: Treat the ODE endpoint as an implicit function of parameters, compute the gradient via the implicit function theorem. Requires solving a linear system at the endpoint.

**Recommendation**: Start with checkpointing (easiest in MLX), migrate to adjoint method for production.

### End-to-end gradient architecture

```python
def compute_cl_differentiable(params):
    """
    params: dict with {H0, omega_b, omega_c, A_s, n_s, tau_reion}

    Returns: C_l array, fully differentiable w.r.t. params
    """
    # Background (MLX)
    bg = solve_friedmann_mlx(params)

    # Recombination (MLX, differentiable tanh or MLX HyRec)
    visibility = compute_visibility_mlx(bg, params)

    # Perturbations (MLX ODE with checkpointing)
    source = solve_boltzmann_mlx(bg, visibility, params)

    # Bessel (Chebyshev in MLX)
    jl = bessel_chebyshev_mlx(ell_values, k_arr, bg.D_A)

    # C_l (MLX einsum)
    Cl = compute_cl_mlx(source, jl, k_arr)

    return Cl

# Gradient w.r.t. cosmological parameters:
grad_fn = mx.grad(lambda p: loss(compute_cl_differentiable(p), data))
```

This enables:
- **Hamiltonian Monte Carlo (HMC)**: 10-100x more efficient than Metropolis-Hastings for high-dimensional parameter spaces
- **Gradient-based optimization**: Find best-fit parameters in seconds instead of hours
- **Fisher matrix**: Exact Fisher matrix from the Hessian, no finite differences needed
- **Neural posterior estimation**: Train normalizing flows with exact likelihood gradients

### Accuracy improvement

Differentiability itself does not improve the accuracy of C_l computation. But it enables:
- More thorough exploration of parameter space (better posteriors)
- Exact derivatives for parameter forecasting (no finite-difference noise)
- Joint optimization with neural network emulators (train emulator with exact gradients from the solver)

### Difficulty: HARD (full end-to-end), MEDIUM (partial -- differentiable C_l integration + emulator)

---

## 9. Implementation Phases

### Phase 0: Foundation fixes (current -> v0.3) -- 1-2 weeks

**Goal**: Fix the current ~30% peak position error without changing the architecture.

| Task | Effort | Impact |
|------|--------|--------|
| Replace tanh x_e with Peebles ODE (3-variable, numpy) | 2 days | 5% -> 1% x_e error |
| Extend ODE integration to tau_rec (using IMEX solver, already written) | 1 day | Fix 15% endpoint error |
| Add ISW integral (e^{-kappa} * Phi' post-recombination) | 3 days | Fix l < 30, add late-ISW |
| Add Doppler term with j_l' projection | 2 days | Fix peak position shifts |
| Validate against CLASS reference spectrum (chi^2 comparison) | 1 day | Quantify residual error |

Expected outcome: peak positions within 5% of CLASS, amplitude within 10%.

### Phase 1: Precision first-order (v0.3 -> v0.5) -- 4-6 weeks

**Goal**: Match CLASS accuracy at 0.1% level.

| Task | Effort | Impact |
|------|--------|--------|
| Full photon hierarchy (l_max = 50) | 1 week | 0.5% -> 0.01% truncation error |
| Full neutrino hierarchy (l_max = 30) | 1 week | 5% -> 0.1% neutrino error |
| HyRec 3-level atom on CPU | 2 weeks | 5% -> 0.1% recombination error |
| Full LOS integration (500 tau-points) | 2 weeks | Fix all peak positions |
| Polarization (E-mode hierarchy) | 1 week | New: EE, TE spectra |
| Mixed precision (float64 Poisson accumulator) | 2 days | Numerical stability |
| Automated comparison pipeline vs CLASS output | 3 days | Continuous validation |

Expected outcome: C_l^{TT} within 0.1% of CLASS for l < 2500. EE and TE spectra available.

### Phase 2: Surpass CLASS (v0.5 -> v1.0) -- 2-3 months

**Goal**: Exceed CLASS accuracy while maintaining GPU speed advantage.

| Task | Effort | Impact |
|------|--------|--------|
| Photon hierarchy l_max = 100 | 1 week | Eliminate truncation entirely |
| Neutrino hierarchy l_max = 100 | 1 week | Eliminate UFA approximation |
| Multi-level recombination (n_max = 30) | 3 weeks | 0.1% -> 0.01% recomb |
| Exact curved-sky lensing | 2 weeks | 0.1% -> 0.01% lensing |
| Reionization model (tau_reion parameterized) | 1 week | Low-l EE signal |
| Massive neutrinos (50 q-bins, l=50 per bin) | 2 weeks | m_nu constraints |
| ISW-lensing correlation | 1 week | 0.5% correction at l~100 |
| Helium recombination (multi-level He I) | 2 weeks | 0.02% C_l correction |

Expected outcome: C_l^{TT} accuracy < 0.01% (10x better than CLASS's ~0.1%). Total wall time < 5 seconds on Apple Silicon.

### Phase 3: Differentiable & extensible (v1.0 -> v2.0) -- 3-6 months

**Goal**: Full differentiability and easy extension to novel physics.

| Task | Effort | Impact |
|------|--------|--------|
| Chebyshev Bessel functions in MLX | 2 weeks | Enable differentiability |
| Background solver in pure MLX | 1 week | Enable differentiability |
| ODE adjoint method | 3 weeks | Memory-efficient gradients |
| Plugin architecture for equation modules | 2 weeks | Easy novel physics |
| Khronon perturbation equation (pi, pi') | 1 week | First non-LCDM model |
| EDE scalar field module | 2 weeks | Second non-LCDM model |
| mu-Sigma modified gravity module | 1 week | Third non-LCDM model |
| HMC sampler with autodiff | 2 weeks | Full inference pipeline |
| Fisher forecasting module | 1 week | Parameter constraints |

Expected outcome: End-to-end differentiable Boltzmann solver with gradient-based MCMC. Novel physics models via plugin architecture. Speed: < 5s per spectrum, < 1ms per gradient evaluation (amortized).

---

## Summary: CLASS vs mlx-class target comparison

| Feature | CLASS v3.x | mlx-class v0.2 (now) | mlx-class v1.0 (target) | mlx-class v2.0 (target) |
|---------|-----------|---------------------|------------------------|------------------------|
| **Speed** | ~10s (1 CPU) | 0.3-2s (GPU) | < 3s (GPU) | < 5s (GPU) |
| **TT accuracy** | ~0.1% | ~30% (peaks off) | < 0.01% | < 0.01% |
| **EE/TE** | Yes | No | Yes (0.1%) | Yes (0.01%) |
| **Photon l_max** | 12-25 | 2 (TCA) | 100 | 100 |
| **Neutrino** | l=17 + UFA | Background only | l=100 exact | l=100 + massive |
| **Recombination** | HyRec 3-level | Tanh fit | HyRec 3-level | Multi-level n=30 |
| **LOS integration** | Full (3000 pts) | Instantaneous | Full (500 pts) | Full (1000 pts) |
| **Lensing** | First-order flat-sky | None | Exact curved-sky | Exact + ISW-lensing |
| **Reionization** | Yes | No | Yes | Yes |
| **Differentiable** | No | No | No | Yes (full pipeline) |
| **Novel physics** | Limited | Khronon (c_s^2=0) | Khronon + MG | Khronon + MG + EDE |
| **Platform** | Any (C) | macOS only (MLX) | macOS (MLX) | macOS (MLX) + CUDA? |

The path from v0.2 to v1.0 is primarily about implementing the standard first-order Boltzmann physics correctly on GPU. The path from v1.0 to v2.0 is about surpassing CLASS through brute-force higher resolution (more multipoles, more recombination levels) and enabling new capabilities (differentiability, novel physics plugins) that are architecturally impossible in CLASS's serial C codebase.

---

## Appendix A: Key References

1. Ma & Bertschinger (1995), ApJ 455, 7 -- Boltzmann hierarchy equations
2. Seljak & Zaldarriaga (1996), ApJ 469, 437 -- Line-of-sight integration (CMBFAST)
3. Blas, Lesgourgues & Tram (2011), JCAP 07, 034 -- CLASS code paper
4. Ali-Haimoud & Hirata (2011), PRD 83, 043513 -- HyRec recombination
5. Lesgourgues & Tram (2011), JCAP 09, 032 -- CLASS perturbation module
6. Blanchet & Skordis (2025), arXiv:2507.00912 -- Khronon-Tensor DBI kinetic structure
7. Lewis, Challinor & Lasenby (2000), ApJ 538, 473 -- CAMB and lensing
8. Peebles (1968), ApJ 153, 1 -- Recombination of hydrogen in the hot Big Bang

## Appendix B: MLX-specific Performance Notes

- Apple M1 Pro: ~2.6 TFLOPS (float32), ~0.08 TFLOPS (float64)
- Apple M3 Max: ~4.2 TFLOPS (float32), ~0.13 TFLOPS (float64)
- Apple M4 Ultra (expected): ~8+ TFLOPS (float32)
- GPU memory: unified with CPU (16-192 GB depending on model)
- `mx.eval()` forces synchronization; batch operations to minimize sync points
- `mx.compile()` can fuse operations for 2-3x speedup on compute-bound kernels
- No native batched tridiagonal solver; implement via parallel cyclic reduction or Thomas algorithm with batching over k

## Appendix C: Risk Assessment

| Risk | Likelihood | Mitigation |
|------|-----------|------------|
| MLX autodiff insufficient for full ODE | Medium | Fall back to checkpointing + finite diff for hard cases |
| Memory overflow for large l_max * N_k * N_tau | Low | Chunked computation over ell |
| Numerical instability in Bessel recurrence | Medium | Use backward recurrence (Miller's algorithm) + Chebyshev fallback |
| MLX float64 too slow for Poisson equation | Low | Keep Poisson accumulation on CPU (numpy float64), transfer result to GPU |
| HyRec port complexity underestimated | Medium | Phase 1 uses Peebles 3-level on CPU; full HyRec is Phase 2 |
| CLASS comparison reveals systematic bias | Low | Unit tests at each integration stage; compare transfer functions, not just C_l |
