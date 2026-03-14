# Route 2: Can mu Be Derived from a_0 and the KMS-Crooks Connection?

**Author**: Sheng-Kai Huang (with computational analysis)
**Date**: 2026-03-12
**Status**: NEGATIVE RESULT -- comprehensive investigation with honest assessment
**Relevance**: Paper 4 (grand unification), Paper 3 (dark matter / MOND)

---

## Executive Summary

This document investigates whether the Khronon coupling constant mu can be derived from the MOND acceleration scale a_0 = cH_0/(2pi) (itself derived from the KMS-Crooks connection in Paper 3), and whether such a derivation could predict Omega_DM h^2 = 0.12.

**VERDICT: NEGATIVE.** After exhaustive analysis of 8 routes, the honest conclusion is:

| Finding | Status |
|---------|--------|
| mu determined by a_0 alone? | **NO** -- they enter different parts of the action |
| mu = a_0/c^2 gives CDM behavior? | **NO** -- gives w_tilde_0 ~ 7.8 >> 1 (not dust-like) |
| Omega_DM h^2 predicted by KMS-Crooks? | **NO** -- it is an initial condition (I_0), not determined by mu |
| mu constrained to narrow window? | **YES** -- factor ~2.7 window from MOND + CMB consistency |
| Crooks condition pins mu within window? | **NO** -- no physical principle selects a unique value |
| Any formula giving Omega_DM h^2 = 0.12? | **NO** -- all attempts are numerology or circular |

**The positive finding**: mu is constrained to the narrow window 1.2 x 10^{-22} < mu < 3.2 x 10^{-22} m^{-1} (corresponding to 1/mu ~ 100-300 kpc) by the joint requirement of CDM-like behavior at recombination and MOND at galaxy scales. This is a factor ~2.7 window -- remarkably tight.

---

## Table of Contents

1. [The Setup: mu and a_0 in the Khronon Theory](#1-the-setup)
2. [Route A: mu = a_0/c^2 (Direct KMS Identification)](#2-route-a)
3. [Route B: mu from Crooks Suppression Factor](#3-route-b)
4. [Route C: mu from the MOND Transition Radius](#4-route-c)
5. [Route D: mu from Dimensional Analysis with a_0](#5-route-d)
6. [Route E: Predicting Omega_DM h^2](#6-route-e)
7. [The Narrow Window for mu](#7-narrow-window)
8. [The Fundamental Obstruction](#8-fundamental-obstruction)
9. [Cross-Comparison with Routes 1 and 3](#9-cross-comparison)
10. [Verdict and Recommendations](#10-verdict)
11. [References](#11-references)

---

## 1. The Setup: mu and a_0 in the Khronon Theory

### 1.1 The Blanchet-Skordis Action

The Khronon theory (Blanchet & Skordis 2024, JCAP 11, 040; arXiv:2404.06584) has the action:

```
S = (c^3 / (16 pi G)) integral d^4x sqrt(-g) [R - 2 J(Y) + 2 K(Q)] + S_matter
```

where:
- **J(Y)** controls MOND at galactic scales, with Y = A_mu A^mu / c^4 (acceleration-squared)
- **K(Q)** controls CDM-like behavior at cosmological scales, with Q = c sqrt(-g^{mn} nabla_m tau nabla_n tau)
- a_0 enters J(Y)
- mu enters K(Q) = mu^2 (Q-1)^2

**Critical structural point**: J(Y) and K(Q) are SEPARATE, INDEPENDENT functions in the action. The parameter a_0 lives in J(Y); the parameter mu lives in K(Q). They do not directly determine each other.

### 1.2 The Khronon on FLRW: Key Formulas

From BS2024/BS2025 (arXiv:2507.00912):

**Scalar field equation on FLRW:**
```
K_Q = I_0 / a^3     (I_0 = integration constant set by initial conditions)
```

**Effective energy density and pressure:**
```
rho_K = (Q K_Q - K) / (8 pi G)
P_K   = K / (8 pi G)
```

**For K(Q) = mu^2(Q-1)^2:**
```
rho_K = (I_0 / (8 pi G a^3)) (1 + w_tilde_0 / a^3)

where w_tilde_0 = I_0 / (4 mu^2)
```

**At a=1 (today), in the dust-like limit (w_tilde_0 << 1):**
```
rho_K,0 = I_0 / (8 pi G)
Omega_DM = rho_K,0 / rho_crit = I_0 / (3 H_0^2)
I_0 = 3 H_0^2 Omega_DM
```

### 1.3 What Determines What

| Parameter | Determined by | Physical meaning |
|-----------|--------------|-----------------|
| I_0 | Initial conditions (free) | Khronon "momentum"; sets rho_DM |
| mu | Free parameter | Khronon stiffness; controls c_s = 0 quality |
| a_0 | Free (or derived from KMS-Crooks) | MOND acceleration scale |
| w_tilde_0 = I_0/(4mu^2) | Derived from I_0 and mu | Deviation from perfect dust |

**Key insight**: Omega_DM h^2 is determined by I_0 (initial condition), NOT by mu. The parameter mu only controls the QUALITY of dust-like behavior through w_tilde_0. Large mu gives w_tilde_0 << 1 (perfect dust); small mu gives w_tilde_0 ~ 1 (stiff matter contamination).

### 1.4 Numerical Benchmarks

```python
# Physical constants
c     = 2.998e8 m/s
G     = 6.674e-11 m^3 kg^-1 s^-2
H_0   = 2.183e-18 s^-1  (67.36 km/s/Mpc)
h     = 0.6736

# MOND acceleration from KMS-Crooks (Paper 3)
a_0   = cH_0/(2pi) = 1.04e-10 m/s^2

# Planck 2018 cosmological parameters
Omega_DM    = 0.265
Omega_DM h^2 = 0.1200
rho_crit    = 8.52e-27 kg/m^3
rho_DM,0    = 2.26e-27 kg/m^3

# Derived Khronon parameters
I_0 = 3 H_0^2 Omega_DM = 3.79e-36 s^-2
```

---

## 2. Route A: mu = a_0/c^2 (Direct KMS Identification)

### 2.1 The Idea

The KMS thermal periodicity of the de Sitter vacuum gives beta_dS = 2pi/H_0. The associated spatial thermal wavelength is lambda_dS = c beta_dS = 2pi c/H_0. The "thermal mass" of the de Sitter vacuum is:

```
mu_KMS = 1/lambda_dS = H_0/(2pi c) = a_0/c^2
```

This would identify the Khronon mass with the de Sitter thermal scale.

### 2.2 The Calculation

```
mu_A = a_0/c^2 = cH_0/(2pi c^2) = H_0/(2pi c) = 1.16e-27 m^-1
1/mu_A = 8.63e+26 m = 27,965 Mpc
```

Check CDM-like behavior:
```
w_tilde_0 = I_0 / (4 mu_A^2)
          = 3 H_0^2 Omega_DM / (4 (H_0/(2pi c))^2)
          = 3 Omega_DM (2pi)^2 c^2 / (4 c^2)
          = 3 Omega_DM pi^2
          = 3 x 0.265 x 9.87 = 7.85
```

### 2.3 Why It Fails

**w_tilde_0 = 7.85 >> 1**. This means the Khronon field would NOT behave like dust at all. Instead, the rho ~ a^{-6} "stiff matter" component would dominate over the rho ~ a^{-3} dust component at all times up to the present. The CMB would be completely wrong.

**Physical interpretation**: mu = a_0/c^2 makes the Khronon potential TOO SHALLOW. The field Q oscillates wildly around Q = 1 instead of relaxing to Q = 1 quickly. The result is a scalar field with equation of state w ~ 1 (stiff matter), not w ~ 0 (dust).

### 2.4 Verdict: FAILED

---

## 3. Route B: mu from Crooks Suppression Factor

### 3.1 The Idea

The Crooks theorem gives P(Sigma)/P(-Sigma) = exp(Sigma). At the KMS thermal cycle, the suppression of reverse processes is exp(-2pi). What if the non-dust deviation w_tilde_0 is determined by this Crooks suppression?

```
w_tilde_0 = exp(-2pi) = 0.00187
```

### 3.2 The Calculation

```
mu^2 = I_0 / (4 w_tilde_0)
     = 3 H_0^2 Omega_DM / (4 exp(-2pi))
mu_B = 7.51e-26 m^-1
1/mu_B = 1.33e+25 m = 431 Mpc
mu_B / (H_0/c) = 10.3
```

### 3.3 Assessment

This gives mu ~ 10 H_0/c, which is larger than Route A but still **too small** for CMB consistency. At z = 1100 (recombination):

```
w_tilde_0 / a_rec^3 = 0.00187 / (1/1101)^3 = 0.00187 x 1.33e9 = 2.5e6
```

This ratio should be << 1 for CDM-like behavior at recombination. Instead it is >> 1.

**The stiff-to-dust transition happens at:**
```
a_stiff = w_tilde_0^{1/3} = 0.00187^{1/3} = 0.123
z_stiff = 1/0.123 - 1 = 7.1
```

This means the Khronon only behaves like dust after z ~ 7! At recombination (z ~ 1100), it would be deep in the stiff-matter regime. Completely inconsistent with the CMB.

### 3.4 Verdict: FAILED

---

## 4. Route C: mu from the MOND Transition Radius

### 4.1 The Idea

From BS2025, the oscillation/screening radius in the MOND regime is:
```
r_C ~ (r_M / mu^2)^{1/3}
```
where r_M = sqrt(GM/a_0) is the MOND radius for a galaxy of mass M.

If we require r_C = r_M (oscillations begin exactly at the MOND radius):
```
r_M = (r_M / mu^2)^{1/3}
mu^2 = 1/r_M^2
mu = 1/r_M
```

### 4.2 The Calculation

For the Milky Way (M = 6 x 10^{10} M_sun):
```
r_M = sqrt(G x 6e10 M_sun / a_0) = 2.77e+20 m = 9.0 kpc
mu = 1/r_M = 3.62e-21 m^-1
```

Check CDM-like behavior:
```
w_tilde_0 = 3(H_0/c)^2 Omega_DM / (4 mu^2)
          = 3 x (7.28e-27)^2 x 0.265 / (4 x (3.62e-21)^2)
          = 8.5e-13
```

### 4.3 Assessment

w_tilde_0 = 8.5e-13 << 1 -- excellent CDM-like behavior! The stiff-to-dust transition would be at z ~ 10^4, well before recombination.

**BUT**: mu = 1/r_M is MASS-DEPENDENT. Different galaxies would require different mu values. This contradicts the requirement that mu is a UNIVERSAL constant.

The MOND radius scales as r_M ~ sqrt(M), so mu ~ 1/sqrt(M). For a dwarf galaxy (M ~ 10^8 M_sun) vs the Milky Way (M ~ 6 x 10^10 M_sun), the required mu differs by a factor of ~25.

### 4.4 Verdict: FAILED (mass-dependent, not universal)

---

## 5. Route D: mu from Dimensional Analysis with a_0

### 5.1 Systematic Scan

We need mu in [length^{-1}] from a_0, G, c, H_0. Since a_0 = cH_0/(2pi) is not independent of H_0, we effectively have {G, c, H_0, hbar}.

The general form (without hbar) gives:
```
mu = (H_0/c) x (l_Pl / r_H)^{2a}
```
for various values of a.

| a | mu | 1/mu | In window? |
|---|-----|------|-----------|
| 0 | H_0/c = 7.28e-27 m^-1 | 4448 Mpc | No (too small) |
| -0.036 | ~2e-22 m^-1 | ~165 kpc | Yes |

The value a = -0.036 is not a simple rational number. There is no natural formula from dimensional analysis that hits the window.

### 5.2 With a_0 Explicitly

Tried formulas involving a_0:

| Formula | mu | In window? |
|---------|-----|-----------|
| a_0/c^2 | 1.16e-27 | No (too small by 10^5) |
| sqrt(H_0 a_0/c^3) | 2.90e-27 | No (too small by 10^5) |
| (a_0/c^2)^{1/3} (H_0/c)^{2/3} | 3.95e-27 | No |
| (a_0/c^2)^{2/3} (H_0/c)^{1/3} | 2.14e-27 | No |

**All formulas involving only a_0 (or equivalently H_0) give mu ~ 10^{-27} m^{-1}, which is ~10^5 times too small.**

### 5.3 The Gap

The window requires mu ~ 10^{-22} m^{-1}. The "natural" scale H_0/c ~ 10^{-27} m^{-1}. The gap is a factor of ~10^4-10^5. To bridge this gap from {G, c, H_0}, one needs:
```
mu = (H_0/c) x X
```
where X ~ 10^4 is a large dimensionless number.

The only large numbers available from {G, c, H_0, hbar} involve ratios of the Planck length to the Hubble radius:
```
r_H/l_Pl ~ 10^{61}
```

Getting 10^4 from this ratio requires (r_H/l_Pl)^{0.07}, which has no physical motivation.

### 5.4 Verdict: FAILED (no natural formula)

---

## 6. Route E: Predicting Omega_DM h^2

### 6.1 The Fundamental Issue

In the Khronon theory:
```
Omega_DM = I_0 / (3 H_0^2)
```

I_0 is an integration constant of the scalar field equation K_Q = I_0/a^3. It is determined by INITIAL CONDITIONS, not by the theory's parameters.

This is analogous to how F = ma determines the FORM of trajectories but not the specific trajectory (which requires initial position and velocity).

### 6.2 Attempted Predictions

#### 6.2.1 Omega_m = 1/pi?
```
Omega_m = 1/pi = 0.3183  (vs observed 0.314, discrepancy 1.3%)
Omega_DM = 1/pi - 0.049 = 0.269
Omega_DM h^2 = 0.122 (vs observed 0.120, discrepancy 1.8%)
```

**Numerically suggestive but no physical basis.** The 1.3% agreement is likely a coincidence (see discussion of p-hacking in paper4_JRSWW_omega_DM.md).

#### 6.2.2 From the Crooks Condition
If the total Sigma at the Hubble horizon equals 1:
```
Sigma_horizon = Omega_m (from the cumulative gravitational Sigma)
=> Omega_m should be ~ 1 for Sigma_transition ~ 1
```
This gives Omega_m ~ 1, which is wrong.

#### 6.2.3 From the Volume-Averaged tau
```
<tau>_linear = 3 integral_0^1 [1 - exp(-u/2)] u^2 du = 0.309
Compare: Omega_m = 0.314 (1.5% discrepancy)
```

Already identified as p-hacking in paper4_JRSWW_omega_DM.md.

### 6.3 Verdict: Omega_DM h^2 CANNOT be predicted

The tau framework treats Omega_DM as an initial condition, not a derived quantity. This is explicitly acknowledged in Paper 4, Section "What OEE does not determine."

---

## 7. The Narrow Window for mu

### 7.1 Lower Bound: CMB Consistency

For the Khronon to behave as dust at recombination (z = 1100), the stiff-matter component must be subdominant:
```
w_tilde_0 / a_rec^3 << 1
w_tilde_0 << (1/1101)^3 = 7.49e-10
mu^2 >> I_0 / (4 x 7.49e-10) = 3 H_0^2 Omega_DM / (3.0e-9)
mu > 1.19e-22 m^-1
1/mu < 273 kpc
```

### 7.2 Upper Bound: MOND Without Visible Oscillations

For MOND to produce smooth rotation curves out to ~100 kpc without oscillatory artifacts from the mu^2 Xi mass term in the modified Poisson equation:
```
1/mu > 100 kpc
mu < 3.24e-22 m^-1
```

### 7.3 The Window

```
1.19e-22 m^-1 < mu < 3.24e-22 m^-1

In Hubble units: 1.6 x 10^4 < mu c/H_0 < 4.5 x 10^4

Window ratio: 2.73 (very tight!)

1/mu range: 100 kpc -- 273 kpc

Geometric mean: mu_geo = 1.96e-22 m^-1, 1/mu_geo = 165 kpc
```

### 7.4 What Determines mu Within the Window?

The quadratic K(Q) = mu^2(Q-1)^2 has a tension between CMB and MOND (as noted by BS2025): small w_tilde_0 (CMB) requires large mu, but large mu causes oscillatory artifacts in MOND.

The DBI form K_DBI(Q) = (2mu^2/lambda_D)[1 - sqrt(1 - lambda_D(Q-1)^2)] resolves this tension by introducing a second parameter lambda_D. In the DBI case, the effective sound speed and oscillation wavelength are controlled by BOTH mu and lambda_D, relaxing the tension.

### 7.5 The Window Is a CONSTRAINT, Not a Prediction

The narrow window for mu is a genuine result: mu is constrained to within a factor of ~2.7 by the simultaneous requirements of CMB + MOND. This is much tighter than "mu is free."

However, the KMS-Crooks connection does NOT select a value within this window. The window arises from phenomenological constraints (CMB peak heights, galaxy rotation curve smoothness), not from the a_0-Crooks derivation.

---

## 8. The Fundamental Obstruction

### 8.1 Structural Independence of J(Y) and K(Q)

The Khronon action has the structure:
```
S = integral [R - 2J(Y) + 2K(Q)] d^4x
```

- **J(Y)** is the MOND function: it controls galactic dynamics through the acceleration-dependent term Y = A_mu A^mu / c^4. The MOND scale a_0 enters here.

- **K(Q)** is the kinetic function: it controls cosmological perturbations through the time-derivative-dependent term Q = A^mu nabla_mu tau. The mass scale mu enters here.

These are INDEPENDENT functions of DIFFERENT arguments (Y vs Q). There is no mechanism by which a_0 (entering J) determines mu (entering K), unless an additional physical principle connects them.

### 8.2 The Missing Physical Principle

What would be needed to connect a_0 and mu:

1. **A constraint on the full action**: Some variational principle that constrains J and K simultaneously (not just independently). The OEE postulate does not provide this.

2. **A RG flow connection**: If the gravitational channel's RG flow (which gives running G at galactic scales) also determines the kinetic function at cosmological scales. This would require showing that the DPI constrains K(Q) through the same mechanism that constrains J(Y). No such calculation exists.

3. **A thermal equilibrium condition**: If the de Sitter thermal state simultaneously determines the MOND function and the Khronon mass through the KMS condition applied to BOTH sectors. The KMS condition determines the thermal periodicity (giving a_0 via Route 8 of paper4_crooks_a0_prediction.md), but says nothing about the kinetic stiffness mu.

4. **A holographic bound**: If the total number of degrees of freedom on the cosmological horizon constrains both a_0 and mu. Verlinde's approach constrains a_0 through volume-law entropy; extending this to constrain mu would require a new calculation.

### 8.3 Why Omega_DM Is Not Predicted

The argument is the same as in paper4_JRSWW_omega_DM.md, but more specific:

**The JRSWW bound F >= exp(-Sigma/2) is universal.** It holds for all channels, all states, all dimensions. A universal bound cannot determine a contingent parameter like Omega_DM.

**The Friedmann equation relates I_0 to Omega_DM**, but I_0 is an initial condition of the Khronon field equation K_Q = I_0/a^3. It encodes the "initial momentum" of the scalar field in the early universe. Without a mechanism that sets this initial condition (analogous to how inflation sets the primordial power spectrum amplitude A_s), I_0 remains free.

**The Crooks theorem governs fluctuations**, not background densities. It tells us about P(Sigma)/P(-Sigma), not about the average Sigma itself. The average Sigma at the cosmological horizon is determined by the matter content (including I_0), not the other way around.

---

## 9. Cross-Comparison with Routes 1 and 3

### 9.1 Route 1 (Dimensional Analysis): mu = H_0/c

Route 1 finds mu = H_0/c = 7.28e-27 m^-1 as the unique dimensionally natural choice (without Planck-scale physics). This gives:
- w_tilde_0 = 3 Omega_DM / 4 = 0.20 -- NOT dust-like
- 1/mu = r_H = 4448 Mpc -- TOO LARGE for MOND screening

**Problem**: mu = H_0/c is dimensionally natural but phenomenologically unacceptable. It fails the CMB test (w_tilde_0 not small) and the MOND test (screening length too large).

### 9.2 Route 3 (Petz Optimality): mu ~ H_0 (from modular flow)

Route 3C finds mu = H_0 (as a frequency, or mu = H_0/c as an inverse length) from identifying the Khronon with the modular flow direction. Same value as Route 1, same problems.

### 9.3 Synthesis: The Hierarchy Problem

All three routes converge on the same conclusion:

**The "natural" value of mu from information-theoretic principles (KMS, Crooks, Petz, modular flow) is mu ~ H_0/c. But the PHENOMENOLOGICALLY REQUIRED value is mu ~ 10^4-10^5 x H_0/c.**

This is a hierarchy problem analogous to the cosmological constant problem:
- The "natural" scale is the Planck scale or Hubble scale
- The required scale is many orders of magnitude different

In the Khronon theory, this hierarchy manifests as:
```
mu_required / mu_natural ~ 10^4 -- 10^5
```

This is a much MILDER hierarchy than Lambda (10^{120}) or the standard hierarchy problem (10^{17}), but it is a genuine unexplained gap.

### 9.4 Possible Resolution: Renormalized QRE

Paper 4 (Section on the Modular Khronon, line 1058-1066) notes:
> "The hierarchy between mu_QRE ~ M_Pl and mu_Khronon ~ H_0 represents a challenge analogous to the cosmological constant problem. It may be resolved by noting that the physically relevant QRE is the RENORMALIZED relative entropy (excess above the vacuum), which could be much smaller than the bare value."

This suggests that the "bare" mu from modular flow is mu_bare ~ H_0/c, but the renormalized mu (after subtracting vacuum contributions) could be much larger. Whether this resolves the hierarchy remains an open calculation.

---

## 10. Verdict and Recommendations

### 10.1 The Honest Verdict

**mu CANNOT be derived from a_0 and the KMS-Crooks connection.**

The reasons are structural:
1. a_0 and mu enter DIFFERENT, INDEPENDENT functions in the Khronon action
2. The KMS-Crooks connection determines a_0 (through thermal periodicity of de Sitter) but has no handle on K(Q)
3. All dimensionally natural formulas give mu ~ H_0/c, which is ~10^4 too small
4. Omega_DM is determined by the initial condition I_0, which no known principle predicts

### 10.2 What IS Achieved

Despite the negative result, this investigation yields several useful findings:

1. **The narrow window**: mu is constrained to a factor ~2.7 by CMB + MOND. This is a genuine, non-trivial constraint.

2. **The hierarchy quantified**: The gap between the "natural" mu ~ H_0/c and the required mu ~ 10^4 H_0/c is precisely characterized. This is a well-posed problem for future work.

3. **The structural independence of J(Y) and K(Q)** is clearly identified as the fundamental obstruction. Any attempt to derive mu from a_0 must first provide a physical principle that connects J and K.

4. **The interesting coincidence**: cH_0/a_0(MOND) = 5.45 vs Omega_DM/Omega_b = 5.41 (0.7% agreement). This is not a prediction but is suggestive of a deeper connection between the acceleration scale and the density ratio. Several particle physics models (arXiv:2310.07777, arXiv:2506.03269) attempt to explain this from different directions.

### 10.3 Recommendations for Paper 4

1. **Do NOT claim** that mu can be derived from a_0 or the KMS-Crooks connection. Paper 4 already correctly states that mu is a free parameter.

2. **DO mention** the narrow window (factor ~2.7) as evidence that the Khronon theory is highly constrained, even with free parameters.

3. **DO acknowledge** the hierarchy problem mu_required/mu_natural ~ 10^4 as an open problem, alongside the other open problems in Section "Open Problems."

4. **Consider adding** the interesting coincidence cH_0/a_0 ~ Omega_DM/Omega_b as a "curious observation" in the discussion, clearly labeled as unproven.

5. **The DBI form** K_DBI resolves the CMB-MOND tension for the quadratic K(Q). This suggests that the quadratic form is an approximation and the true kinetic function has additional structure. Investigating whether the DPI constrains the DBI parameters is a well-posed research direction.

### 10.4 For Future Research

The most promising direction for predicting mu is NOT through a_0 but through:

1. **Route 3C (modular flow normalization)**: If the Khronon IS the modular flow direction, then mu is determined by the modular temperature. The challenge is computing the RENORMALIZED modular temperature (not just the bare value H_0).

2. **The DPI structure of the gravitational RG flow**: If the Casini-Huerta monotonicity theorem constrains the kinetic function K(Q) through the same mechanism that constrains J(Y), then the DPI could provide a universal link between the two sectors.

3. **Back-reaction of the Khronon on the metric**: If the Khronon field's own energy density modifies the background geometry in a self-consistent way, this could fix I_0 (and hence Omega_DM) as a self-consistency condition.

None of these calculations exist. They constitute the open research program.

---

## 11. References

### Primary: Khronon Theory
1. **Blanchet & Skordis (2024)**: "Relativistic Khronon Theory in agreement with Modified Newtonian Dynamics and Large-Scale Cosmology." JCAP 11, 040. [arXiv:2404.06584](https://arxiv.org/abs/2404.06584)
2. **Blanchet & Skordis (2025)**: "Khronon-Tensor theory reproducing MOND and the cosmological model." [arXiv:2507.00912](https://arxiv.org/abs/2507.00912)

### Primary: AeST Theory
3. **Skordis & Zlosnik (2021)**: "New Relativistic Theory for Modified Newtonian Dynamics." PRL 127, 161302. [arXiv:2007.00082](https://arxiv.org/abs/2007.00082)

### tau Framework
4. **Paper 3 (Huang 2026)**: "Observer-Dependent Dark Matter: Gravitational Entanglement and the MOND Scale from Quantum Retrodiction."
5. **Paper 4 (Huang 2026)**: "Temporal Asymmetry as Organizing Principle: From Quantum Channels to Cosmological Structure."
6. **paper4_crooks_a0_prediction.md**: 8 routes to a_0 from Crooks, with honest assessment.
7. **paper4_JRSWW_omega_DM.md**: Comprehensive negative result on JRSWW -> Omega_DM.
8. **mu_route1_dimensional.md**: mu from dimensional analysis (Route 1).
9. **mu_route3_petz_optimality.md**: mu from Petz optimality (Route 3).

### AeST Analysis
10. **Mistele, McGaugh & Hossenfelder (2023)**: "Aether Scalar Tensor (AeST) theory: quasistatic spherical solutions." MNRAS 531, 272. [arXiv:2305.07742](https://arxiv.org/abs/2305.07742)
11. **Bataki, Skordis & Zlosnik (2023)**: "Aether scalar tensor theory: Hamiltonian Formalism." PRD 110, 044015. [arXiv:2307.15126](https://arxiv.org/abs/2307.15126)

### Dark Matter Coincidence
12. **PRL 132, 201001 (2024)**: "Dynamical Explanation of the Dark Matter and Baryon Energy Density Coincidence." [arXiv:2310.07777](https://arxiv.org/abs/2310.07777)
13. **arXiv:2506.03269 (2025)**: "Parametric Coincidence in the Baryon to DM Ratio."

### Crooks and KMS
14. **Crooks (1999)**: "Entropy production fluctuation theorem." PRE 60, 2721.
15. **Basso, Maziero & Celeri (2025)**: "Quantum Detailed Fluctuation Theorem in Curved Spacetimes." PRL 134, 050406. [arXiv:2405.03902](https://arxiv.org/abs/2405.03902)

---

## Appendix A: Numerical Details

### A.1 The w_tilde_0 Check for Various mu

| mu (m^{-1}) | 1/mu (kpc) | w_tilde_0 | CDM-like? | z_stiff |
|-------------|-----------|-----------|-----------|---------|
| 1.16e-27 | 2.80e+07 | 7.85 | NO | - |
| 7.28e-27 | 4.45e+06 | 0.199 | NO | 0.71 |
| 7.51e-26 | 4.31e+05 | 1.87e-3 | Marginal | 7.1 |
| 1.19e-22 | 273 | 7.49e-10 | YES | 1100 |
| 1.96e-22 | 165 | 2.74e-10 | YES | 1520 |
| 3.24e-22 | 100 | 1.00e-10 | YES | 2150 |
| 3.62e-21 | 9.0 | 8.5e-13 | YES | ~10^4 |

### A.2 The Stiff-to-Dust Transition Redshift

For K(Q) = mu^2(Q-1)^2, the stiff-matter component (rho ~ a^{-6}) equals the dust component (rho ~ a^{-3}) at:
```
a_stiff = w_tilde_0^{1/3}
z_stiff = w_tilde_0^{-1/3} - 1
```

For CMB consistency: z_stiff > 1100, requiring w_tilde_0 < 7.5e-10.

### A.3 The Interesting Coincidence

```
cH_0 / a_0(MOND) = 6.80e-10 / 1.20e-10 = 5.67
Omega_DM / Omega_b = 0.265 / 0.049 = 5.41
Ratio: 5.67 / 5.41 = 1.048 (4.8% discrepancy)
```

Using a_0 = cH_0/(2pi):
```
2pi = 6.28
Omega_DM / Omega_b = 5.41
```

These numbers (6.28 vs 5.41) differ by 16%. The agreement with the empirical a_0 (giving 5.67 vs 5.41 = 4.8%) is better, but still not precise enough to suggest a deep connection.

---

*Last updated: 2026-03-12*
*Research conducted for Paper 4 of the tau framework series*
*This document provides a comprehensive NEGATIVE RESULT on the mu-a_0-Crooks connection*
