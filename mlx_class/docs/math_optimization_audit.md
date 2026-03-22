# Mathematical Optimization Audit: mlx_class

**Date**: 2026-03-23
**Scope**: Full Boltzmann solver pipeline from perturbations to C_l
**Method**: Eigenvalue analysis, sparsity profiling, algorithmic cost modeling

---

## 1. A-Matrix Structure and Sparsity

The Boltzmann system y'(tau) = A(tau) y has dimension 56 x 56:
- 4 metric/matter variables: eta, delta_c, delta_b, theta_b
- 26 photon multipoles: F_g,0 ... F_g,25
- 26 neutrino multipoles: F_n,0 ... F_n,25

### Sparsity

At k=0.1, tau=200: **167 / 3136 nonzero (5.3% fill)**

| Block | Size | Fill | max\|A\| | Structure |
|-------|------|------|----------|-----------|
| Metric internal | 4x4 | 56% | 2.79 | Dense |
| Metric -> Photon | 26x4 | 7.7% | 14.88 | h' coupling |
| Metric -> Neutrino | 26x4 | 6.7% | 2.10 | h' coupling |
| Photon -> Metric | 4x26 | 3.8% | 0.21 | eta' feedback |
| **Photon internal** | **26x26** | **11.4%** | **1.12** | **Tridiagonal + diagonal** |
| Photon -> Neutrino | 26x26 | 0.4% | 0.0025 | Negligible |
| Neutrino -> Metric | 4x26 | 2.9% | 0.0013 | Negligible |
| Neutrino -> Photon | 26x26 | 0.4% | 0.0017 | Negligible |
| **Neutrino internal** | **26x26** | **7.8%** | **0.10** | **Purely tridiagonal** |

**Key structural finding**: The photon and neutrino hierarchies for l >= 3 are **purely tridiagonal** (verified: off-tridiagonal elements = 0 to machine precision). The coupling between photons and neutrinos is negligible (max 0.005). The metric sector couples to hierarchies only through h' and eta' constraints, creating a **rank-7 dense correction** to an otherwise nearly block-diagonal system.

### Decomposition: A = A_diag + A_tridiag + A_metric

| Component | max\|eigenvalue\| | Structure | Stiffness source |
|-----------|------------------|-----------|-----------------|
| A_diag (Thomson) | 2.79 (at tau=200) | Diagonal | -\|kappa_dot\| on photon l >= 1 |
| A_tridiag (streaming) | 0.10 | Tridiagonal | +/- i*k oscillatory |
| A_metric (constraints) | rank 7 | Low-rank dense | k^2/calH Poisson |

---

## 2. Eigenvalue Analysis

### Eigenvalue phase diagram (k = 0.1 Mpc^-1)

| tau (Mpc) | \|kappa_dot\| | max\|eig\| | Stiffness | Regime | n_negative | n_oscillatory |
|-----------|--------------|-----------|-----------|--------|------------|---------------|
| 10 | 928 | 62881 | >10^15 | Tight coupling | 54 | 50 |
| 50 | 32.5 | 439 | >10^14 | Tight coupling | 54 | 50 |
| 100 | 6.4 | 44 | >10^13 | Tight coupling | 54 | 50 |
| 200 | 1.1 | 3.9 | 660 | Transition | 54 | 50 |
| 250 | 0.31 | 0.89 | 150 | Transition | 54 | 52 |
| 270 | 0.11 | 0.30 | 50 | Decoupling | 54 | 52 |
| 280 | 0.06 | 0.15 | 26 | Decoupling | 54 | 52 |
| 300 | 0.003 | 0.10 | 17 | Free streaming | 54 | 52 |
| 500 | ~0 | 0.10 | 17 | Free streaming | 54 | 52 |

### Eigenvalue structure by regime

**Tight coupling (tau < 170 Mpc)**: 25 large negative eigenvalues from Thomson scattering (photon l >= 1 damped at rate -\|kd\|). The remaining 31 eigenvalues are small (acoustic oscillation + CDM growth). Stiffness ratio > 10^6.

**Transition (170 < tau < 273 Mpc)**: Thomson eigenvalues fade from O(10) to O(0.1). Stiffness drops from 10^3 to 50. Still needs implicit treatment.

**Free streaming (tau > 295 Mpc)**: All eigenvalues have \|Re\| < 0.01. 52 oscillatory modes (neutrino + photon streaming). Stiffness ~ 17. Purely explicit-friendly.

### Critical finding: The stiff-to-non-stiff transition is SHARP

The stiffness drops from 10^6 to 17 in a ~60 Mpc window around recombination (tau ~ 235-295). This means:
- 85% of the integration range (tau < 235) is dominated by Thomson stiffness
- 10% (235-295) is the transition region
- 5% (295-330) matters for the visibility function, and is non-stiff

---

## 3. Optimal Time-Stepping Strategy

| Phase | tau range | Method | Cost/step | Why | Steps needed |
|-------|-----------|--------|-----------|-----|--------------|
| **Early** | 0.01 - 235 Mpc | TCA-0 Magnus | O(31^3) | Eliminates Thomson stiffness by not evolving photon hierarchy. Log-spaced grid natural. | ~200 |
| **Transition** | 235 - 295 Mpc | IMEX or Implicit | O(56) or O(56^3) | Thomson fading but still O(1). Need photon quadrupole for Silk damping. | ~30 |
| **Free streaming** | 295 - 330 Mpc | Explicit RK4 | O(56 * N_k) | dt_max = 2.8/0.10 = 28 Mpc. Only ~2-3 steps needed. Visibility peak in this range. | ~5 |
| **Late ISW** | 330 - 3000 Mpc | Explicit RK4 | O(56 * N_k) | Sparse snapshots sufficient. Source varies slowly. | ~15 |

**Current approach**: Two-phase (TCA-0 Magnus + full Magnus or NDF15).

**Proposed improvement**: Three-phase with IMEX for transition.

**Impact**: The transition region currently forces either (a) high-order implicit solver for the full system (expensive) or (b) TCA-0 past its validity (inaccurate). An IMEX splitting that treats Thomson damping implicitly (diagonal solve = trivial) and streaming explicitly (tridiagonal = cheap) would handle this region at O(n) cost per step per mode instead of O(n^3).

### IMEX splitting details

The A-matrix decomposes cleanly:
```
A_implicit = diag(-|kd|, -|kd|, ..., 0, 0, ...)  [Thomson on photon sector]
           + rank-7 correction [metric coupling]
A_explicit = tridiagonal [streaming k*l/(2l+1)]
```

Implicit solve: `(I - dt * A_implicit)^{-1}` is diagonal + rank-7, solvable via Sherman-Morrison in O(n + 7^2) per mode.

---

## 4. Analytical Source Function Feasibility

The Hu-Sugiyama (1996) analytical approximation gives:

```
(Theta_0 + Psi)(k, tau_rec) = A_SW(k) * cos(k*r_s + phi(k)) * D_silk(k) + zero_shift(k)
```

| Component | Formula | Accuracy |
|-----------|---------|----------|
| Acoustic oscillation | cos(k * r_s) | Peak position < 0.5% |
| Phase shift phi(k) | phi_inf * (1 - exp(-alpha * k/k_eq)) | 1% positions, 5% heights |
| Silk damping | exp(-(k/k_D)^2) | 10% at k = k_D, worse at higher k |
| Driving amplitude | D(k) from Eisenstein-Hu transfer | 5-10% |
| Baryon zero-shift | R/(3(1+R)) * T(k) | 5% |

**Overall**: Peak positions to 1-2%, peak heights to 5-10%, damping tail to 10-20%.

**Verdict**: NOT sufficient for sub-percent precision. ODE through recombination is required for the full pipeline. However, the analytical form is valuable for:
1. **Validation**: Quick check of ODE output
2. **Initial guess for iterative solvers**
3. **Gradient preconditioning** in MCMC

---

## 5. Visibility-Weighted Integration Range

The visibility function g(tau) = |kd| * exp(-kappa) peaks at tau = 280.4 Mpc with FWHM = 22.5 Mpc.

| Threshold | tau range (Mpc) | Width (Mpc) | Contribution |
|-----------|----------------|-------------|--------------|
| g > 1% of peak | 245 - 346 | 101 | SW + Doppler (99%) |
| g > 0.1% of peak | 239 - 518 | 279 | + early ISW wing |
| g > 0.01% of peak | 235 - 8978 | 8743 | + reionization bump |

**Optimal quadrature** (for the SW + Doppler integral):
- 25 Gauss-Legendre points over [230, 330] Mpc: captures 99%+ of g(tau)
- 15 points over [137, 230] for early ISW
- 20 points over [330, 3000] for late ISW
- **Total: 60 snapshots** (vs current ~110 in solver_sync_gpu.py)

The early ISW window (tau_eq to tau_rec) contributes to l ~ 100-300 (the "early ISW plateau"). The late ISW window (tau > tau_rec) contributes to l < 20 (the late ISW rise).

---

## 6. k-Mode Importance Sampling

The integrand for C_l peaks at k ~ l / D_A with width ~ sqrt(l) / D_A.

| ell | k_peak (Mpc^-1) | dk_width | k_range needed |
|-----|-----------------|----------|----------------|
| 2 | 0.00014 | 0.00045 | [0.00001, 0.002] |
| 100 | 0.0072 | 0.00072 | [0.004, 0.011] |
| 220 | 0.016 | 0.001 | [0.011, 0.021] |
| 1000 | 0.072 | 0.002 | [0.061, 0.083] |
| 2500 | 0.180 | 0.004 | [0.162, 0.198] |

**Current**: 500 k-modes from 3e-4 to 0.35, geomspaced. For any given ell, only ~50-100 k-modes are in the relevant range. The others contribute negligibly.

**Optimization options**:
1. **Per-ell k-selection**: For each ell, only sum over k-modes with |k - l/D_A| < 5*dk. Reduces effective N_k per ell from 500 to ~80.
2. **Adaptive k-grid**: Place k-points densely at the acoustic peaks of the source function (period 2*pi/r_s ~ 0.044 Mpc^-1). Need ~10 points per oscillation cycle.
3. **k-ell joint sparse sampling**: Compute Delta_l(k) only where the Bessel function is non-negligible (j_l(k*D_A) > epsilon).

**Estimated savings**: 3-5x reduction in the C_l integration cost.

---

## 7. Bessel Function Optimization

### Current: Piecewise Chebyshev (bessel_gpu.py)
- Table build: ~2s (forward recurrence + scipy transition fix + DCT)
- GPU eval: ~0.026s for 372 ells x 5000 k-points (570x vs scipy)
- Memory: ~6.5 MB

### Optimization 1: Extended Limber approximation

For l > 30, the Bessel integral can be replaced by the Limber approximation:
```
Delta_l(k) ~ sqrt(pi/(2l+1)) * S(k, tau_l) / k
```
where tau_l satisfies k * chi(tau_l) = l + 1/2.

| ell range | Error (standard Limber) | Error (extended 2nd order) |
|-----------|------------------------|---------------------------|
| l = 10 | 1.2% | 0.007% |
| l = 20 | 0.6% | 0.002% |
| l = 50 | 0.2% | ~0% |
| l > 100 | < 0.1% | ~0% |

**Strategy**: Full Bessel for l = 2..30 (29 values), Limber for l = 31..2500.
**Impact**: Bessel evaluations reduced by **92%** (from 372 to 29 ells).

### Optimization 2: Precomputed j_l tensor

Precompute `j_l(k * chi_i)` for all (l, k, tau) at once: shape (372, 500, 60) float32 = 45 MB.
Then the LOS integral becomes a single einsum/matmul operation.

### Optimization 3: GPU-native forward recurrence

For the 29 low-ell values still needing full Bessel, use forward recurrence on GPU:
```python
j_{l+1}(x) = (2l+1)/x * j_l(x) - j_{l-1}(x)
```
This is stable for x > 1.5*l and fully vectorized over all x-points.

---

## 8. Autodiff Feasibility

| Pipeline Step | Differentiable? | Blocker | Fix |
|--------------|----------------|---------|-----|
| Background (Friedmann) | YES | numpy (can rewrite in MLX) | Rewrite ~100 lines |
| Recombination (RECFAST) | NO | scipy ODE | Precompute table, interpolate in MLX |
| Adiabatic ICs | YES | Pure algebra | Trivial |
| TCA-0 A-matrix build | YES | numpy arrays | Rewrite in MLX |
| **Magnus matrix exp** | **YES** | **MLX-native** | **Already differentiable** |
| NDF15 batched | NO | Newton iteration | Replace with Magnus |
| **RK4 steps** | **YES** | **MLX-native** | **Already differentiable** |
| scipy solve_ivp | NO | Python callback | Replace with Magnus/RK4 |
| Gauge transform | YES | Pure algebra | Trivial |
| Phi' + Psi' derivative | PARTIAL | CubicSpline | Replace with finite diff on MLX |
| **Bessel Chebyshev eval** | **YES*** | **MLX gather + arccos** | Need to verify gather grad |
| LOS integral | YES | MLX sum + multiply | Trivial |
| C_l integration | YES | Simple sum | Trivial |

### Path to full autodiff

The **Magnus expansion solver** is the key enabler:
```
params -> A(tau; params) -> Omega = integral A dtau -> expm(Omega) -> y(tau_i)
```

Every step is differentiable:
- `A(tau; params)`: linear in cosmological parameters
- `Omega`: quadrature (differentiable)
- `expm(Omega)`: matrix exponential has well-defined derivative
- Propagation: `y_{n+1} = expm(Omega_n) @ y_n` is a chain of matmuls

This enables: **dC_l / d(Omega_b, Omega_c, h, A_s, n_s, tau_reio)** in one backward pass.

**Impact**: Gradient-based MCMC (Hamiltonian Monte Carlo) for 6-parameter cosmological inference. HMC with gradients is typically **10-100x faster** than Metropolis-Hastings, especially in high dimensions.

**Effort**: Medium (2-3 weeks). Requires:
1. Rewrite background in MLX (~1 day)
2. Precompute recombination table (~1 day)
3. Ensure all intermediate operations are MLX-native (~1 week)
4. Implement backward pass and validate gradients (~1 week)

---

## 9. GPU Utilization Bottlenecks

Current pipeline timing (solver_sync_gpu.py, 500 k-modes):

| Step | Time | GPU/CPU | Bottleneck |
|------|------|---------|-----------|
| Background | 0.8s | CPU | numpy + scipy ODE |
| ODE (NDF15) | 2-5s | GPU | Matrix inverse on CPU stream |
| Gauge transform | 1.2s | Mixed | **Python loop over tau snapshots** |
| Phi' + Psi' | 0.5s | CPU | **CubicSpline loop over N_k modes** |
| LOS integration | 3.0s | GPU | **Loop over tau snapshots** |
| **Total** | **7-10s** | | |

### Fixes (no algorithmic changes, just engineering):

| Fix | Current | Proposed | Expected speedup |
|-----|---------|----------|-----------------|
| Batch gauge transform | Loop + mx.eval per snapshot | Single batched mx.eval | 10x -> 0.12s |
| GPU finite difference | scipy CubicSpline per mode | mx.diff on GPU tensor | 100x -> 0.005s |
| Precompute Bessel tensor | Loop over tau, eval per snapshot | One batch eval, then matmul | 3-5x -> 0.6s |
| Eliminate CPU<->GPU transfers | Multiple np<->mx round-trips | Keep everything on MLX graph | 2x overall |

**Combined engineering speedup**: 7-10s -> **2-3s** (3-4x faster without any algorithmic change).

---

## 10. Hierarchy Truncation

Current: lg_max = ln_max = 25 (fixed for all k-modes).

The source function only depends on F_g,0, F_g,1, F_g,2 (monopole, dipole, quadrupole). Higher multipoles matter only through:
1. Back-reaction on F_g,2 via tridiagonal streaming
2. Silk damping (information leaking to high l)

| k (Mpc^-1) | k * tau_dec | Optimal lg_max | Savings vs 25 |
|------------|------------|----------------|---------------|
| 0.001 | 0.3 | 6 | 76% |
| 0.01 | 3.0 | 6 | 76% |
| 0.05 | 15 | 24 | 4% |
| 0.1 | 30 | 25 | 0% |
| 0.35 | 105 | 25 | 0% |

**Strategy**: For k < 0.05 (low-ell modes, ~50% of k-range in geomspace), use lg_max = 8. For k > 0.05, use lg_max = 25.

**Impact**: ~40% reduction in system size for half the k-modes -> ~20% overall ODE speedup.

**Caveat**: The TCA-0 approach already eliminates the photon hierarchy for the early phase. This optimization mainly helps the full-hierarchy phase.

---

## 11. Recommended Optimizations (Ranked by Impact)

| # | Optimization | Speed gain | Accuracy gain | Effort | Risk |
|---|-------------|-----------|--------------|--------|------|
| **1** | **IMEX splitting for transition region** | 3-5x on ODE | Removes TCA error in transition | 1-2 weeks | Low (well-studied method) |
| **2** | **Limber approximation for l > 30** | 5-10x on LOS | < 0.01% error for l > 30 | 2-3 days | Very low |
| **3** | **Batch GPU operations (engineering)** | 3-4x overall | None | 3-5 days | Very low |
| **4** | **Precomputed Bessel tensor + matmul LOS** | 3-5x on LOS | Identical | 2-3 days | Very low |
| **5** | **k-mode importance sampling** | 2-3x on C_l | Same or better | 1 week | Low |
| **6** | **Autodiff via Magnus pipeline** | 10-100x for MCMC | Enables HMC | 2-3 weeks | Medium |
| **7** | **Adaptive hierarchy truncation** | 20% on ODE | Negligible loss | 2-3 days | Low |
| **8** | **Gauss-Legendre quadrature for visibility** | 2x on snapshot count | Better (optimal quadrature) | 1 day | Very low |
| **9** | **Matrix exp: n_terms 16->12** | 25% on expm | ~float32 limited | 1 hour | None |

### Combined impact estimate

Applying optimizations 1-5 and 7-9:
- ODE phase: 5x faster (IMEX + adaptive truncation + engineering)
- LOS phase: 10x faster (Limber + precomputed tensor + k-sampling)
- Overall: **current 7-10s -> 1-2s** for 500 k-modes, l_max = 2500

Adding optimization 6 (autodiff):
- MCMC parameter estimation: **current ~hours -> minutes** per chain

---

## Appendix A: Eigenvalue Data Tables

### Full eigenvalue spectrum at k=0.1, tau=50 (tight coupling)

- 25 Thomson-damped eigenvalues: Re(lambda) in [-439, -29.3]
  - Largest (most negative): lambda = -439.3 (h' constraint amplification)
  - Cluster at Re = -32.5 = -|kd| (pure Thomson)
- 31 small eigenvalues: |Re| < 0.04, |Im| < 0.10
  - Acoustic oscillation pair: Im ~ +/- 0.058 (sound waves)
  - CDM growing mode: Re ~ -0.04 (gravitational instability)
  - Neutrino free-streaming: |Im| ~ 0.01-0.10

### Full eigenvalue spectrum at k=0.1, tau=300 (free streaming)

- 0 Thomson-damped eigenvalues (decoupled)
- 56 small eigenvalues: |Re| < 0.01, |Im| < 0.10
  - 52 oscillatory (photon + neutrino streaming)
  - 4 near-zero (metric perturbations decaying)
  - Stiffness ratio: 17 (benign)

## Appendix B: Block-Diagonalizability

At tau=50 (tight coupling), the photon l >= 2 sector has:
- Internal eigenvalue range: max|eig| = 32.5 (responsible for 7% of total stiffness)
- Coupling to rest of system: max|A| = 0.25 (outgoing), 0.07 (incoming)

The dominant stiffness source is NOT the photon internal block, but the **h' constraint coupling** which amplifies the Thomson rate by a factor of calH / k^2, producing the largest eigenvalue (-439 at tau=50).

This confirms: the TCA-0 approach (eliminating the photon hierarchy) is the correct strategy for the early phase. The residual stiffness from h' constraint requires either implicit treatment or Magnus expansion.

## Appendix C: Tridiagonal Structure Verification

Post-decoupling (tau=300), the photon and neutrino streaming hierarchies for l = 3..25 are **exactly tridiagonal**:
```
A[l, l-1] = k * l / (2l+1)      [upward streaming]
A[l, l+1] = -k * (l+1) / (2l+1)  [downward streaming]
```
Off-tridiagonal elements: 0 to machine precision.

This confirms that tridiagonal solvers (Thomas algorithm) can be used for the streaming part in an IMEX splitting, at O(n) cost instead of O(n^3).
