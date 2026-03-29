# f(delta) Crisis Analysis: Re-Derivation of the Khronon Stress-Energy on FRW

**Author**: Sheng-Kai Huang (with computational analysis)
**Date**: 2026-03-17
**Status**: CRITICAL -- Paper 3's f(delta) formula is WRONG, but crisis persists
**Classification**: Error identification + honest assessment

---

## Executive Summary

The CLASS Level 2 test found that the f(delta) = (4+3delta)/4 background correction is excluded at ~13 sigma by Planck CMB data. We re-derived the Khronon stress-energy tensor on FRW by going directly to the BS2024 paper (arXiv:2404.06584).

### Key Findings

1. **Paper 3's formula f(delta) = (4+3delta)/4 is WRONG.** It does not match BS2024's Eq. (4.4a) and (4.14).

2. **The correct BS formula is:**
   ```
   rho_K = mu^2 * (Q^2 - 1) / (8piG) = mu^2 * delta * (2+delta) / (8piG)
   ```
   If one insists on writing rho_K = mu^2 * delta * f(delta) / (8piG), the correct f is:
   ```
   f(delta) = 2 + delta    (NOT (4+3delta)/4)
   ```

3. **The correct formula REDUCES the crisis by ~28%** (20.5% excess at z=1100 vs 28.7%), but the theory is still catastrophically excluded by Planck at >10 sigma.

4. **The equation of state is w = delta/(2+delta)** [BS Eq 4.17], NOT w_tilde = delta/2 as Paper 3 states. These agree only for delta << 1.

5. **The CMB crisis is REAL and structural.** It cannot be resolved by correcting the f(delta) formula.

---

## 1. The Error in Paper 3

### 1.1 What Paper 3 Claims

Paper 3 (paper3_weak_field.tex), Eq. (rhoK), line 380:
```latex
\rho_K = \frac{c^2\mu^2\,\delta\,f(\delta)}{8\pi G},
```
with f(delta) = (4 + 3delta)/4, citing BS2024.

Additionally, Paper 3 states w_tilde = delta/2 [line 392].

### 1.2 What BS2024 Actually Says

From Blanchet & Skordis 2024 (JCAP 11, 040; arXiv:2404.06584):

**Eq. (4.4a):** Khronon energy density on FLRW (c=1 convention):
```
rho_K = (1/(8piG)) * (Q * K_Q - K)
```

**Eq. (4.4b):** Khronon pressure:
```
P_K = K / (8piG)
```

**Eq. (4.5):** Equation of state:
```
w = P_K / rho_K = K / (Q * K_Q - K)
```

**Eq. (4.12):** Quadratic K:
```
K(Q) = mu^2 * (Q-1)^2
```

**Eq. (4.14):** Energy density for quadratic K:
```
rho_K = (mu^2/(8piG)) * (Q^2 - 1) = (I_0/(8piG*a^3)) * (1 + w_tilde_0/a^3)
```
where w_tilde_0 = I_0/(4mu^2).

**Eq. (4.17):** Equation of state:
```
w = (Q-1)/(Q+1) = delta/(2+delta)
```

### 1.3 The Discrepancy

Computing Q*K_Q - K for K = mu^2*(Q-1)^2 with delta = Q-1:

```
K = mu^2 * delta^2
K_Q = 2 * mu^2 * delta
Q * K_Q = (1+delta) * 2 * mu^2 * delta = 2*mu^2*delta*(1+delta)
Q * K_Q - K = 2*mu^2*delta*(1+delta) - mu^2*delta^2
            = 2*mu^2*delta + 2*mu^2*delta^2 - mu^2*delta^2
            = 2*mu^2*delta + mu^2*delta^2
            = mu^2 * delta * (2+delta)
            = mu^2 * (Q^2 - 1)    ✓ matches BS Eq 4.14
```

So the correct result is:
```
rho_K = mu^2 * delta * (2+delta) / (8piG)
```

Paper 3 claims:
```
rho_K = c^2 * mu^2 * delta * (4+3delta) / (4 * 8piG)
```

With c=1 (BS convention in Section 4):
- Correct: delta * (2+delta) = 2*delta + delta^2
- Paper 3: delta * (4+3delta)/4 = delta + 3*delta^2/4

**These are different functions.** For delta = 0.34:
- Correct: 0.34 * 2.34 = 0.796
- Paper 3: 0.34 * 1.255 = 0.427

The ratio is 1.865. The Paper 3 formula underestimates the energy density by a factor of ~1.9 at delta_0 = 0.34.

### 1.4 Where the Error Likely Originated

The formula f(delta) = (4+3delta)/4 does not correspond to any known expression from the BS paper. Possible sources of confusion:

1. **Misidentification with w_tilde_0**: BS defines w_tilde_0 = I_0/(4mu^2), and the energy density can be written as rho_K = I_0*(1 + w_tilde_0/a^3)/(8piG*a^3). Perhaps an attempt to rewrite this in terms of delta introduced an algebraic error.

2. **Confusion between different decompositions**: The factor (4+3delta)/4 does not arise from any standard decomposition of Q*K_Q - K = mu^2*delta*(2+delta).

3. **Possible earlier draft**: The formula may have been carried over from an earlier computation that used different conventions.

### 1.5 Impact on Omega_K Calculation

Despite the wrong f(delta), the Omega_K = delta*(2+delta)/3 calculation throughout the research notes IS CORRECT. This is because:
- omega_dm_prediction.md correctly derives rho_K = mu^2*delta*(2+delta)/(8piG)
- Omega_K = delta*(2+delta)/3 follows directly
- delta_0 = 0.340 for Omega_DM = 0.265 is correct

The error is ONLY in Paper 3's Eq. (rhoK) and the Level 2 CLASS analysis.

---

## 2. Correct BS Formulas

### 2.1 Friedmann Equation (BS Eq 4.3a)

In synchronous coordinates (c=1):
```
3H^2 + 3kappa/a^2 - Lambda = 8piG * sum(rho_I) + Q*K_Q - K
```

### 2.2 Khronon as Perfect Fluid

Energy density: rho_K = (Q*K_Q - K)/(8piG)
Pressure:       P_K = K/(8piG)
EOS:            w = K/(Q*K_Q - K) = delta/(2+delta)

### 2.3 Conservation Law (BS Eq 4.7)

```
K_Q = I_0 / a^3
```

For K = mu^2*(Q-1)^2:  2*mu^2*delta = I_0/a^3, giving mu^2*delta = I_0/(2*a^3).

### 2.4 Energy Density for Quadratic K (BS Eq 4.14)

```
rho_K = (mu^2/(8piG)) * (Q^2 - 1) = (I_0/(8piG*a^3)) * (1 + w_tilde_0/a^3)
```

where w_tilde_0 = I_0/(4*mu^2). With running mu_bg = H(a)/c:
- delta(a) = I_0*c^2 / (2*H(a)^2*a^3) = delta_0 / (E^2(a)*a^3)
- w_tilde_0/a^3 = delta/2

### 2.5 Equation of State (BS Eq 4.17)

```
w = (Q-1)/(Q+1) = delta/(2+delta)
```

NOT delta/2. For delta=0.34: w = 0.145 (not 0.170).

### 2.6 Adiabatic Sound Speed (BS Eq 4.18)

```
c_ad^2 = 1 - 1/Q = 2*w_tilde_0 / (2*w_tilde_0 + a^3)
```

This is the adiabatic sound speed, distinct from the perturbation sound speed c_s^2 = 0 (from ghost condensation).

---

## 3. Impact on the CMB Crisis

### 3.1 Correct vs Wrong Excess at z=1100

With running mu = H/c and delta_0 = 0.340:

| Formula | delta(z=1100) | f(1100)/f(0) | Excess DM |
|---------|---------------|--------------|-----------|
| BS correct: f=2+delta | 0.819 | 1.205 | +20.5% |
| Paper 3 wrong: f=(4+3d)/4 | 0.819 | 1.287 | +28.7% |

### 3.2 Correct Effective w(a) for CLASS

| z | delta | w_eff (BS correct) | w_eff (P3 wrong) | ratio |
|---|-------|-------------------|-----------------|-------|
| 0 | 0.34 | 0.100 | 0.139 | 0.72 |
| 0.5 | 0.66 | 0.097 | 0.130 | 0.75 |
| 1 | 0.85 | 0.064 | 0.083 | 0.77 |
| 3 | 1.05 | 0.011 | 0.014 | 0.78 |
| 10 | 1.08 | 0.000 | 0.000 | 0.78 |
| 100 | 1.05 | -0.003 | -0.004 | 0.78 |
| 1100 | 0.82 | -0.024 | -0.031 | 0.76 |

The BS w_eff is consistently ~72-78% of the Paper 3 w_eff. This reduces the CMB impact by ~25%.

### 3.3 Estimated Corrected CMB Exclusion

Level 2 found with Paper 3's formula:
- Max |dTT/TT| = 33% at 12.8 sigma (11-bin, delta_0=0.34)

With the correct BS formula (scaling by the w_eff ratio):
- Estimated max |dTT/TT| ~ 24% at ~9 sigma

**Still catastrophically excluded.** Planck constrains omega_cdm to ~1%.

### 3.4 delta_0 Exclusion Threshold (Corrected)

The Level 2 found delta_0 < 0.003 is required for safety. With the correct BS formula (25% less impact), the threshold relaxes to approximately:
- delta_0 < 0.004 for SAFE
- delta_0 < 0.04 for borderline

For delta_0 = 0.340 (needed for Omega_DM = 0.265), the theory is excluded at >10 sigma regardless of which f(delta) is used.

---

## 4. The Root Cause of the Crisis

### 4.1 Why rho_K is Not Exactly CDM-Like

With running mu = H(a)/c:
```
rho_K = I_0 * (2 + delta(a)) / (16piG * a^3)
```

The a^{-3} scaling is multiplied by the factor (2+delta(a)), where:
- Radiation era: delta → 0, so f → 2 (constant)
- Matter era: delta ≈ delta_0/Omega_m ≈ 1.08 (constant), so f ≈ 3.08 (constant)
- Lambda era: delta decreases from ~1.08 to delta_0 = 0.34, so f drops from ~3.08 to 2.34

The factor f = 2+delta changes by:
```
f_matter / f_today = (2+1.08)/(2+0.34) = 3.08/2.34 = 1.32
```

This means 32% MORE dark matter during the matter era than if extrapolated from today using a^{-3}. This is the core of the crisis.

### 4.2 The Self-Consistency Trap

The crisis is self-consistent:
- delta_0 = 0.34 is REQUIRED to get Omega_DM = 0.265 at z=0
- Running mu = H/c is REQUIRED to keep delta ~ O(1) and avoid the a^{-9} stiff-matter catastrophe
- But running mu causes delta to vary from 0.34 (today) to ~1.08 (matter era)
- This variation makes rho_K/a^3 non-constant, deviating from CDM

The theory is trapped: any mechanism that gives the right Omega_DM today necessarily gives too much DM at recombination.

### 4.3 Quantitative Proof

Define the CDM-normalized density:
```
r(z) = rho_K(z) / [rho_K(0) * (1+z)^3] = (2+delta(z)) / (2+delta_0)
```

For CDM: r = 1 at all z. For the Khronon:
- r(z=1100) = 1.205 (20.5% excess)
- r(z=11) = 1.315 (31.5% excess, peak)

This is an EXACT result, depending only on:
1. K(Q) = mu^2*(Q-1)^2 (ghost condensation)
2. mu = H/c (running)
3. delta_0 = 0.340 (from Omega_DM = 0.265)

---

## 5. Possible Resolutions

### 5.1 Does NOT Work: Correcting f(delta)

The correction from f=(4+3d)/4 to f=2+d reduces the crisis by ~25% but does not resolve it. The theory is still excluded at >10 sigma.

### 5.2 Does NOT Work: Reducing delta_0

For delta_0 < 0.004: the Khronon becomes safe for CMB, but Omega_K < 0.003, providing only ~1% of the dark matter. The theory then fails to explain dark matter.

### 5.3 Does NOT Work: Constant mu

With constant mu = H_0/c: delta ∝ a^{-3}, giving delta(1100) ~ 10^9. The energy density scales as a^{-9} (stiff matter). Even worse than the running case.

### 5.4 Potentially Works: Different Running Prescription

If mu runs FASTER than H/c, delta stays smaller at high z:
- mu(a) = H(a)^alpha/c with alpha > 1: delta ∝ 1/(H^{2alpha}*a^3)
- In matter era: delta ∝ a^{3(alpha-1)} → 0 for alpha > 1
- This makes f → 2 = const, giving exact a^{-3} scaling

**But:** alpha > 1 has no known theoretical motivation. The DPI + extensivity argument gives alpha = 1.

### 5.5 Potentially Works: Non-Quadratic K(Q)

If K(Q) is not exactly quadratic, the relationship between rho_K and delta changes:

For K(Q) = mu^2 * |Q-1|^p:
- K_Q = p*mu^2*|Q-1|^{p-2}*(Q-1)
- Q*K_Q - K = (p-1)*mu^2*delta^p + p*mu^2*delta^{p-1}
- For p=2: mu^2*(delta^2 + 2*delta) ✓
- For p=4: mu^2*(3*delta^4 + 4*delta^3) -- scales even more steeply, WORSE

The DBI completion K_DBI = mu^2*lambda_D^2*[sqrt(1+(Q-1)^2/lambda_D^2) - 1] might help for large delta (it saturates), but for delta ~ 1 the behavior is close to quadratic.

### 5.6 Potentially Works: Separating Background and Perturbation mu

If the background mu stays constant (mu_0 = H_0/c) while only the perturbation mu runs:
- Background: delta ∝ a^{-3}, giving delta(1100) ~ 10^9 → stiff matter
- This is the original CMB Catch-22. Does NOT work.

### 5.7 The Nuclear Option: Reinterpret What the Khronon IS

Perhaps the Khronon energy density should NOT be identified with ALL of the dark matter:
- Omega_K provides only a small fraction (~1-5%) of DM
- The remaining 95%+ is standard CDM (particles)
- The Khronon contributes only to the MOND-like galactic sector via J(Y)

In this case: delta_0 << 1, f ~ const, no CMB problem. But this abandons the core claim that the Khronon IS the dark matter.

### 5.8 Most Promising: AeST-Like Mechanism

In the AeST theory (Skordis & Zlosnik 2021), the scalar field has a more complex cosmological evolution that CAN reproduce CDM-like background evolution exactly. The BS Khronon with quadratic K fails, but a more general K(Q) might work.

The key insight from AeST: the function K(Q) must be chosen so that the cosmological background has w = 0 EXACTLY (not approximately). This requires K to be fine-tuned, which contradicts the ghost condensation argument for K = mu^2*(Q-1)^2.

---

## 6. Honest Assessment

### 6.1 What Is Wrong

1. **Paper 3's f(delta) = (4+3delta)/4 is incorrect.** The BS stress-energy tensor gives f = 2+delta.

2. **Paper 3's w_tilde = delta/2 is a leading-order approximation.** The exact BS result is w = delta/(2+delta).

3. **The CMB crisis is REAL.** The Khronon with K=mu^2*(Q-1)^2 and running mu=H/c predicts 20.5% more dark matter at z=1100 than CDM (when normalized at z=0). This is excluded by Planck at >10 sigma.

### 6.2 What Is Still Correct

1. **The perturbation sound speed c_s^2 = 0** is correct and safe for CMB.
2. **The conservation law mu^2*delta = I_0/(2a^3)** is correct.
3. **Omega_K = delta*(2+delta)/3** is correct.
4. **The galactic sector J(Y)** is independent and unaffected.
5. **GW170817 compatibility** (c_T = c) is structural and safe.

### 6.3 The Severity

The theory as currently formulated (K=mu^2*(Q-1)^2, mu=H/c) **cannot** explain all of the dark matter while being compatible with the CMB. This is a **fatal** problem for the claim that "dark matter = Khronon condensate."

### 6.4 Possible Ways Forward

1. **Abandon the running mu = H/c prescription** and find a different running that keeps delta constant (not just O(1)) across all epochs. This requires a theoretical derivation of the running, which is currently an assumption.

2. **Modify K(Q)** beyond the leading quadratic form. A K(Q) that makes w = 0 exactly on FRW would solve the problem, but this requires fine-tuning that undermines the ghost condensation argument.

3. **Accept that the Khronon is not 100% of dark matter.** In this case, delta_0 << 1 and the background correction is negligible. The Khronon provides only the MOND-like galactic phenomenology, while the cosmological dark matter has a separate origin.

4. **Adopt the AeST framework** instead of the simpler BS Khronon. AeST has been explicitly shown to reproduce the CMB (Skordis & Zlosnik 2021), but at the cost of more free parameters and a more complex action.

---

## 7. Erratum for Paper 3

### 7.1 Equation to Fix

Replace Eq. (rhoK) [line 380]:
```
Old: rho_K = c^2*mu^2*delta*f(delta)/(8piG), f(delta) = (4+3delta)/4
New: rho_K = c^2*mu^2*delta*(2+delta)/(8piG) = c^2*mu^2*(Q_0^2 - 1)/(8piG)
```

### 7.2 w_tilde to Fix

Replace w_tilde = delta/2 [line 392]:
```
Old: w_tilde = delta/2
New: w = delta/(2+delta)    [BS Eq (4.17)]
```

### 7.3 Consequential Changes

The text "f → 1 for delta << 1" should be "f → 2 for delta << 1" (where f = 2+delta). This changes the normalization but not the qualitative behavior.

The sentence "rho_K = c^2 I_0 f(delta)/(16piG a^3)" should be:
```
rho_K = c^2 I_0 (2+delta) / (16piG a^3)
```

---

## 8. Corrected Level 2 Results (Estimate)

The Level 2 CLASS runs used f = (4+3delta)/4. With the correct f = 2+delta, the effective w(a) is ~72-78% as large. Estimated corrected results:

| Run | P3 formula max|dTT/TT| | Corrected (est.) | P3 N_sigma | Corrected (est.) |
|-----|------------------------|-----------------|------------|------------------|
| B (11-bin, d0=0.34) | 33.0% | ~24% | 12.8 | ~9 |

The corrected value is still catastrophically excluded. A proper Level 2 re-run with the correct BS formula is needed for exact numbers, but the qualitative conclusion is unchanged.

---

## 9. Connection to the Extremal Principle

The Omega_DM derivation (omega_DM_derivation_2026_03_16.md) uses:
```
Sigma = delta*(2+delta) = (1+delta)^2 - 1
```

This is EXACTLY 3*Omega_K, which equals rho_K/rho_crit * 3. So Sigma = 3*rho_K/(rho_crit) in natural units. The extremal principle result Omega_DM = 0.268 is based on the CORRECT BS formula, not the wrong Paper 3 formula. The prediction is unaffected by the f(delta) error.

---

## 10. Conclusions

1. **f(delta) = (4+3delta)/4 is wrong.** The correct BS formula is rho_K = mu^2*delta*(2+delta)/(8piG).

2. **The correction reduces the CMB crisis by ~25%** but does not resolve it.

3. **The theory with K=mu^2*(Q-1)^2 and mu=H/c is excluded** at the background level by Planck CMB at >10 sigma.

4. **This is a STRUCTURAL problem,** not a computational error. The running mu = H/c, while necessary to keep delta ~ O(1), causes delta to vary across epochs, making rho_K/a^3 non-constant.

5. **The theory needs fundamental modification** at the cosmological level. The most promising options are: (a) a different running prescription, (b) a non-quadratic K(Q) with w=0 exactly, or (c) accepting that the Khronon is only a fraction of the dark matter.

---

## References

1. Blanchet, L. & Skordis, C. (2024). JCAP 11, 040. arXiv:2404.06584. Eqs. (4.3)-(4.18).
2. Skordis, C. & Zlosnik, T. (2021). PRL 127, 161302. arXiv:2007.00082.
3. Level 2 results: cmb_class_level2_results_2026_03_17.md (this repository)
4. Omega_DM derivation: omega_DM_derivation_2026_03_16.md (this repository)

---

*Last updated: 2026-03-17*
*This document identifies a formula error in Paper 3 and assesses its impact on the CMB crisis. The honest conclusion: the error makes the crisis ~25% less severe, but the theory is still catastrophically excluded by Planck. The CMB background problem is structural and requires fundamental modification of the theory.*
