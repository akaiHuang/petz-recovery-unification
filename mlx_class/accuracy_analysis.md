# mlx_class Accuracy Analysis: Systematic Error Decomposition

**Date**: 2026-03-20
**Current RMS residual (l>30)**: 87.7%
**Target**: 0.1%
**Gap**: 877x

## Executive Summary

The 88% RMS residual decomposes into 8 distinct error sources. Two dominate:
1. **ISW template excess at l < 100** (484% RMS) -- accounts for 88% of total RMS^2
2. **TCA limitation + driving correction** (30% combined) -- fundamental physics gap

The single most impactful fix is switching from TCA to the full photon hierarchy
(`perturbations_hires.py`), which already exists but is not wired into the main pipeline.

---

## Error Budget

| # | Source | l-range | RMS | Priority | Fix Effort |
|---|--------|---------|-----|----------|------------|
| 1 | ISW template excess | l < 100 | ~484% | CRITICAL | Low |
| 2 | Driving correction miscalibration | 100-2000 | ~30% | HIGH | Medium |
| 3 | TCA limitation (no Theta_l, l>=2) | all l | ~25% | HIGH | Medium |
| 4 | Silk damping scale error | l > 1200 | ~15% | MEDIUM | Low |
| 5 | Source term construction | all l | ~10% | MEDIUM | Medium |
| 6 | Missing reionization | l > 10 | ~10% | LOW | Trivial |
| 7 | Numerical integration | all l | ~1% | LOW | Trivial |
| 8 | Initial conditions (f_nu) | high l | ~1% | LOW | Trivial |

---

## Detailed Analysis

### Error Source 1: ISW Template Excess (l < 100)

**Severity**: CATASTROPHIC (35x excess at l~8, 484% RMS)

**Evidence**:
- l=10: CLASS=820 muK^2, MLX=28194 muK^2 (34x excess)
- l=50: CLASS=1421, MLX=8139 (5.7x excess)
- l=80: CLASS=2116, MLX=4055 (1.9x excess)
- Dominates 88% of the total RMS^2

**Root Cause**:
The early ISW integration (Phi_dot * j_l(k*chi) from tau_eq to tau_vis) produces
unphysically large oscillations at low-l. The late ISW (Eisenstein-Hu) adds further
excess. In `class_comparison.py`, an empirical Gaussian ISW bump at l~140 with
amplitude 1.3x the first acoustic peak makes things even worse.

**Fix**:
1. Remove the empirical Gaussian ISW template from `class_comparison.py`
2. Normalize Phi_plateau in late ISW to match ODE Phi output at tau_rec
3. Use more integration points for early ISW (12 -> 50)
4. Consider Limber approximation for ISW at l > 30

---

### Error Source 2: Driving Correction Miscalibration

**Severity**: HIGH (39% normalization excess, 24% oscillatory amplitude)

**Evidence** (fit to ratio = a + b*l + c*cos(2pi*l/P + phi)):
- Mean normalization: a = 1.39 (should be 1.00)
- Damping slope: -0.29 per 1000 ell (excess damp at high l)
- Oscillatory amplitude: 24% peak-to-peak at acoustic period
- Phase offset: -49 deg (wrong peak/trough ratio)

**Root Cause**:
The constant D_inf = 1.95 is empirically calibrated but:
1. Overestimates the mean amplitude by 39%
2. Does not account for k-dependence (radiation vs matter era modes)
3. Applied equally to SW and Doppler sources (but they have different phases)
4. Does not decrease at high k where TCA is more accurate

**Fix**:
Best: Eliminate entirely by switching to full hierarchy.
Interim: Replace with k-dependent D(k) calibrated at each peak, reducing D_inf to ~1.65.

---

### Error Source 3: TCA Limitation (Fundamental)

**Severity**: HIGH (drives errors 2, 4, 5)

**Evidence** (peak-to-trough contrast):
- Peak 1 -> Trough 1: CLASS contrast=3.37, MLX=1.51 (MLX/CLASS=0.45)
- Peak 3 -> Trough 3: CLASS contrast=2.56, MLX=1.89 (MLX/CLASS=0.74)
- Mean contrast ratio: 0.85 (troughs are 15% too shallow)

**Root Cause**:
The TCA only tracks Theta_0 and Theta_1 with v_b = 3*Theta_1. This means:
1. No Theta_2 for polarization source (5% at peaks)
2. No free-streaming redistribution (Doppler source overestimated)
3. No photon quadrupole back-reaction (driving amplitude wrong)
4. v_b locked to 3*Theta_1 even during decoupling (baryon drag missed)

**Fix**: THE DEFINITIVE FIX: Switch to `perturbations_hires.py`.
This already implements the full photon Boltzmann hierarchy with Strang splitting.
- l_gamma_max = 8: sufficient for 1% accuracy
- l_gamma_max = 25: sufficient for 0.1% accuracy
- Eliminates error sources 2, 3(all), 5(all), and most of 4

---

### Error Source 4: Silk Damping Scale

**Severity**: MEDIUM (24% RMS at l > 1200)

**Evidence**:
- l=1500: ratio = 0.95 (5% deficit)
- l=2000: ratio = 0.76 (24% deficit)
- l=2400: ratio = 0.58 (42% deficit)
- log(ratio) slope: -0.50 per 1000 ell

**Root Cause**:
1. k_D = 0.15*(omega_b/0.022)^0.25 is approximate (CLASS uses exact integral)
2. Silk damping applied to SOURCE as exp(-(k/k_D)^2), but visibility integration
   already provides damping -> DOUBLE COUNTING
3. Driving correction D_corr doesn't decrease at high k where it's less needed

**Fix**:
1. Compute exact k_D from integral formula (10 lines in background.py)
2. When using full hierarchy, REMOVE explicit Silk damping (hierarchy handles it)
3. If keeping TCA, do not apply Silk damping inside visibility integral

---

### Error Source 5: Source Term Construction

**Severity**: MEDIUM (~10% combined)

**Sub-errors**:
- (A) ISW uses 2*Phi_dot instead of Phi_dot+Psi_dot (~1%)
- (B) v_b = 3*Theta_1 instead of actual v_b from baryon equation (~5%)
- (C) D_corr applied to DC offset Psi, not just oscillatory part (~3%)
- (D) Missing Theta_P2 and g'(tau)*v_b/(3k) terms (~3%)

**Fix**: All resolved by full hierarchy switch.

---

### Error Source 6: Missing Reionization

**Severity**: LOW (10% uniform suppression)

**Evidence**: tau_reio = 0.0544; exp(-2*tau_reio) = 0.897; 10.3% suppression at all l>10.

**Root Cause**: Peebles recombination does not model the z~8 reionization bump.
Currently absorbed into D_corr calibration, but will matter after D_corr is removed.

**Fix**: One line: `Cl *= np.exp(-2 * 0.0544)` or add tanh reionization model.

---

### Error Sources 7-8: Numerical & ICs

- **Numerical**: float32, N_tau_vis=25. Combined < 1%. Fix: float64 source, N_tau_vis=50.
- **ICs**: Missing f_nu correction in dipole IC. ~1% at high l.
  Fix: `Theta_1 = k*tau/(18*(1 + f_nu/4))` in adiabatic_ic().

---

## Roadmap to 0.1%

### Phase 1: Switch to Full Hierarchy (1-2 days) -> ~5-10% RMS

1. Create `solver_hires.py` that uses `perturbations_hires.py`
2. Wire HiResBoltzmannSolver with l_gamma_max=8 into LOS pipeline
3. Remove driving correction D_corr entirely
4. Remove explicit Silk damping from source terms
5. Keep visibility quadrature and Bessel integration from solver_precision.py
6. Fix ISW template (remove empirical Gaussian, normalize late ISW)

### Phase 2: Precision Fixes (0.5 days) -> ~2-5% RMS

1. Compute exact k_D from integral formula
2. Add reionization damping
3. Fix f_nu in initial conditions
4. Use Phi_dot+Psi_dot for ISW (not 2*Phi_dot)

### Phase 3: Fine-Tuning (0.5 days) -> ~0.5-1% RMS

1. Increase l_gamma_max to 25
2. Increase N_tau_vis to 50
3. Float64 source evaluation
4. Proper photon anisotropic stress in Psi diagnostic

### Phase 4: Final Push (1-2 days) -> 0.1% RMS

1. CMB lensing (smooths peaks by ~5% at l > 500)
2. Exact HyRec/RECFAST-level recombination
3. Massive neutrino corrections (small)
4. Convergence testing: vary N_k, N_tau, l_gamma_max

**Total estimated effort: 3-5 days**

---

## Key Insight

The existing codebase already has the most important piece: `perturbations_hires.py`
with the full photon Boltzmann hierarchy and Strang splitting. The main gap is not
physics implementation but **pipeline integration** -- connecting the hierarchy solver
to the LOS integration with proper source term construction.

The driving correction and explicit Silk damping are band-aids for TCA limitations.
Once the full hierarchy replaces TCA, these band-aids must be REMOVED, not adjusted.
This is the single change that takes us from 88% to ~5% RMS.

---

## Files

- Analysis script: `mlx_class/accuracy_analysis.py`
- This document: `mlx_class/accuracy_analysis.md`
- Full hierarchy solver: `mlx_class/perturbations_hires.py` (exists, needs integration)
- Precision solver (current best): `mlx_class/solver_precision.py`
- Background: `mlx_class/background.py` (background quantities are accurate)
- CLASS reference: `gdm_class_public/output/mlx_ref_lcdm_cl.dat`
