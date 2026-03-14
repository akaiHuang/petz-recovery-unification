# Computing Omega_DM h^2 from mu_0 = H_0/c: Complete Algebra

**Author**: Sheng-Kai Huang (with computational analysis)
**Date**: 2026-03-12
**Status**: Complete calculation with honest assessment
**Relevance**: Paper 4 (grand unification), links to mu_route1_dimensional.md and mu_route2_a0_crooks.md

---

## Executive Summary

Starting from mu_0 = H_0/c (the unique dimensionally natural Khronon coupling), we compute Omega_DM h^2 step by step through the Blanchet-Skordis Khronon FRW cosmology. The result:

| Finding | Value | Status |
|---------|-------|--------|
| Omega_K as function of delta = Q_0 - 1 | Omega_K = (delta^2 + 2*delta)/3 | **DERIVED** (algebra) |
| delta needed for Omega_DM = 0.265 | delta = 0.340 | Computed (quadratic solved exactly) |
| delta needed for Omega_m = 0.315 | delta = 0.395 | Computed (quadratic solved exactly) |
| rho_K scaling with redshift | rho_K ~ a^{-3} ONLY IF w_tilde_0 << 1 | **PROBLEM: w_tilde_0 ~ O(1) for mu = H_0/c** |
| w_tilde_0 for mu = H_0/c, Omega_DM = 0.265 | w_tilde_0 = 0.170 | NOT dust-like |
| Omega_DM h^2 = 0.12 predicted? | **NO** | One free parameter (I_0 or equivalently delta) remains |

**CRITICAL FINDING**: mu = H_0/c gives a **cosmological coincidence** (rho_K ~ rho_crit) but **fails the dust-like requirement** (w_tilde_0 ~ 0.2, not << 1). The Khronon with mu = H_0/c does NOT behave like CDM at recombination. This is consistent with Route 2's finding that the phenomenologically required mu is ~ 10^4 times larger than H_0/c.

**THERE IS NO MAJOR RESULT**: Omega_DM h^2 = 0.12 does NOT fall out naturally. It remains a free parameter (initial condition I_0).

---

## 1. Background Cosmology: Khronon on FRW

### 1.1 Setup and Conventions

The Blanchet-Skordis Khronon action (arXiv:2404.06584):

```
S = (c^3 / (16 pi G)) integral d^4x sqrt(-g) [R - 2 J(Y) + 2 K(Q)] + S_m
```

On FRW background ds^2 = -c^2 dt^2 + a(t)^2 d_vec{x}^2:
- The Khronon background: T = t (comoving time)
- Q_0 = c * sqrt(-g^{00} (nabla_0 T)^2) = c * sqrt(1/c^2 * 1) = 1 for comoving choice

**Important**: Blanchet & Skordis (2024, Sec. 2.4) show that the FLRW solution allows Q_0 != 1 if T = alpha * t + f(t) for some function f(t). The background value Q_0 is determined by initial conditions.

### 1.2 Khronon Stress-Energy on FRW

From BS2024, Eq. (2.12)-(2.14), the Khronon contribution to the Friedmann equations:

**Energy density** (appears in H^2 equation):
```
rho_K = (c^4 / (16 pi G)) * 2 * (Q_0 K' - K)       ... (*)
```

where K' = dK/dQ evaluated at Q_0.

**Dimension check on (*):**
- [c^4/(16 pi G)] = [m^4 s^{-4}] / [m^3 kg^{-1} s^{-2}] = [kg m s^{-2}] = [J/m] ???

Let me be more careful. The standard Friedmann equation is:

```
H^2 = (8 pi G / (3 c^2)) * rho
```

where rho has dimensions [energy/volume] = [kg m^{-1} s^{-2}].

From the action with prefactor c^3/(16 pi G), the effective energy density from K(Q) is:

```
rho_K = (c^4 / (8 pi G)) * (Q_0 K' - K)
```

**Dimension check:**
```
[c^4/(8 pi G)] = [m^4 s^{-4}] / [m^3 kg^{-1} s^{-2}] = [kg m s^{-2}]
```

That has dimensions [kg m s^{-2}] = [N]. We need [kg m^{-1} s^{-2}]. Something is off.

**Correction**: K(Q) has dimensions [m^{-2}] (same as R), and Q_0 is dimensionless, so K' = dK/dQ also has dimensions [m^{-2}]. Therefore:

```
[c^4/(8 pi G)] * [K] = [m^4 s^{-4}] / [m^3 kg^{-1} s^{-2}] * [m^{-2}]
                      = [kg m s^{-2}] * [m^{-2}]
                      = [kg m^{-1} s^{-2}]   CHECK (energy density)
```

Good. So:

```
rho_K = (c^4 / (8 pi G)) * (Q_0 K_Q - K)
```

where K_Q = dK/dQ.

### 1.3 For K(Q) = mu^2 (Q - 1)^2

```
K = mu^2 (Q_0 - 1)^2
K_Q = 2 mu^2 (Q_0 - 1)
```

Define delta := Q_0 - 1. Then:

```
Q_0 K_Q - K = (1 + delta) * 2 mu^2 delta - mu^2 delta^2
             = 2 mu^2 delta (1 + delta) - mu^2 delta^2
             = 2 mu^2 delta + 2 mu^2 delta^2 - mu^2 delta^2
             = 2 mu^2 delta + mu^2 delta^2
             = mu^2 delta (2 + delta)
             = mu^2 (Q_0^2 - 1)
```

**Check**: Q_0^2 - 1 = (1+delta)^2 - 1 = 1 + 2*delta + delta^2 - 1 = 2*delta + delta^2 = delta(2 + delta). Consistent.

Therefore:

```
rho_K = (c^4 / (8 pi G)) * mu^2 * delta * (2 + delta)
      = (c^4 mu^2 / (8 pi G)) * (delta^2 + 2 delta)
```

### 1.4 Substituting mu = H_0/c

```
mu^2 = H_0^2 / c^2
```

```
rho_K = (c^4 / (8 pi G)) * (H_0^2 / c^2) * (delta^2 + 2 delta)
      = (c^2 H_0^2 / (8 pi G)) * (delta^2 + 2 delta)
```

**Dimension check:**
```
[c^2 H_0^2 / G] = [m^2 s^{-2}] * [s^{-2}] / [m^3 kg^{-1} s^{-2}]
                 = [m^2 s^{-4}] / [m^3 kg^{-1} s^{-2}]
                 = [m^2 s^{-4} * kg s^2 / m^3]
                 = [kg m^{-1} s^{-2}]   CHECK (energy density)
```

Now recall:
```
rho_crit = 3 c^2 H_0^2 / (8 pi G)
```

**Dimension check on rho_crit:**
```
[c^2 H_0^2 / G] = [kg m^{-1} s^{-2}]   (same as above)   CHECK
```

Wait -- the standard Friedmann equation is H^2 = 8 pi G rho / (3 c^2), so:
```
rho_crit = 3 c^2 H_0^2 / (8 pi G)
```

**Dimension check:**
```
[3 c^2 H_0^2 / (8 pi G)] = [m^2 s^{-2}] * [s^{-2}] / [m^3 kg^{-1} s^{-2}]
                           = [m^2 s^{-4}] * [kg s^2 m^{-3}]
                           = [kg m^{-1} s^{-2}]   CHECK
```

Therefore:

```
Omega_K = rho_K / rho_crit
        = [(c^2 H_0^2 / (8 pi G)) * (delta^2 + 2 delta)] / [3 c^2 H_0^2 / (8 pi G)]
        = (delta^2 + 2 delta) / 3
```

---

## 2. The Master Formula

```
===============================================================
  Omega_K = (delta^2 + 2 delta) / 3     where delta = Q_0 - 1
===============================================================
```

This is EXACT for K(Q) = mu^2(Q-1)^2 with mu = H_0/c on FRW background.

**Special cases:**
- delta = 0: Omega_K = 0 (comoving Khronon, no dark matter)
- delta = 1: Omega_K = (1 + 2)/3 = 1.00
- delta = -1: Omega_K = (1 - 2)/3 = -1/3 (negative energy! unphysical)
- delta << 1: Omega_K ~ 2 delta/3 (linear approximation)
- delta >> 1: Omega_K ~ delta^2/3 (quadratic approximation)

### 2.1 What delta Gives Omega_DM = 0.265?

Solve (delta^2 + 2 delta)/3 = 0.265:
```
delta^2 + 2 delta - 0.795 = 0
delta = [-2 + sqrt(4 + 3.18)] / 2 = [-2 + sqrt(7.18)] / 2 = [-2 + 2.680] / 2 = 0.340
```

Wait, let me redo this carefully:
```
delta^2 + 2 delta = 3 * 0.265 = 0.795
delta^2 + 2 delta - 0.795 = 0
delta = [-2 +/- sqrt(4 + 4 * 0.795)] / 2
      = [-2 +/- sqrt(4 + 3.180)] / 2
      = [-2 +/- sqrt(7.180)] / 2
      = [-2 +/- 2.6796] / 2
```

Taking the positive root:
```
delta = (-2 + 2.6796) / 2 = 0.6796 / 2 = 0.3398
```

**Hmm**, this differs from Route 1's estimate of delta = 0.892. Let me check.

Route 1 (mu_route1_dimensional.md, Section 5.4) states:
```
Omega_K = delta^2/3 = Omega_DM = 0.265  =>  delta = sqrt(3 * 0.265) = 0.892
```

But that uses the QUADRATIC approximation delta^2/3, dropping the 2*delta term. The full formula is (delta^2 + 2*delta)/3.

With the full formula:
```
delta = 0.340  =>  Omega_K = (0.340^2 + 2 * 0.340)/3 = (0.1156 + 0.680)/3 = 0.7956/3 = 0.265   CHECK
```

So **delta = 0.340** gives Omega_DM = 0.265 exactly.

**IMPORTANT CORRECTION**: Route 1's formula rho_K = c^2 mu^2 (Q_0^2 - 1)/(8 pi G) is correct but Route 1's Section 5.4 erroneously dropped the linear term. For the quadratic K(Q) = mu^2(Q-1)^2, the correct expression is:

```
rho_K = (c^4 mu^2 / (8 pi G)) * (Q_0^2 - 1) = (c^4 mu^2 / (8 pi G)) * (delta^2 + 2 delta)
```

The delta^2/3 approximation only works when 2*delta >> delta^2, i.e., delta << 2. Since delta = 0.34 << 2, the linear term 2*delta DOMINATES over delta^2. Route 1's quadratic approximation was inappropriate.

### 2.2 What delta Gives Omega_m = 0.315?

```
delta^2 + 2 delta = 3 * 0.315 = 0.945
delta = [-2 + sqrt(4 + 3.780)] / 2 = [-2 + sqrt(7.780)] / 2 = [-2 + 2.789] / 2 = 0.395
```

Check: (0.395^2 + 2*0.395)/3 = (0.1560 + 0.790)/3 = 0.946/3 = 0.315.  CHECK

### 2.3 Summary Table

| Target | Omega | Required delta | Q_0 | delta << 1? |
|--------|-------|---------------|-----|-------------|
| Omega_DM | 0.265 | 0.340 | 1.340 | Marginal (linear term 0.680 vs quadratic 0.116) |
| Omega_m | 0.315 | 0.395 | 1.395 | Marginal |
| Omega_m = 1/3 | 0.333 | 0.414 | 1.414 = sqrt(2) | Interesting! |

**Note**: delta = sqrt(2) - 1 = 0.414 gives Omega_K = 1/3 exactly. This is because:
```
(delta^2 + 2 delta)/3 = ((sqrt(2)-1)^2 + 2(sqrt(2)-1))/3
= (2 - 2 sqrt(2) + 1 + 2 sqrt(2) - 2)/3
= 1/3
```

So Q_0 = sqrt(2) gives Omega_K = 1/3. This is algebraically clean but there is no known physical principle selecting Q_0 = sqrt(2).

---

## 3. Friedmann Equation: Complete Structure

### 3.1 Modified Friedmann Equation

With the Khronon contribution:

```
H^2 = (8 pi G / (3 c^2)) * (rho_m + rho_r + rho_Lambda + rho_K)
```

Dividing by H_0^2:

```
E^2(z) = Omega_m (1+z)^3 + Omega_r (1+z)^4 + Omega_Lambda + Omega_K(z)
```

where E(z) = H(z)/H_0 and the question is: **how does Omega_K scale with redshift?**

### 3.2 The Key Question: Does rho_K Scale as a^{-3}?

From BS2024 (Sec. 2.4), the Khronon field equation on FRW gives:

```
K_Q = I_0 / a^3
```

where I_0 is an integration constant (determined by initial conditions).

For K(Q) = mu^2 (Q - 1)^2:
```
K_Q = 2 mu^2 (Q_0 - 1) = 2 mu^2 delta
```

Therefore:
```
delta(a) = I_0 / (2 mu^2 a^3)
```

So delta EVOLVES with the scale factor: delta ~ a^{-3}.

### 3.3 rho_K(a) Computation

Substituting delta(a) = I_0/(2 mu^2 a^3) into rho_K:

```
rho_K(a) = (c^4 mu^2 / (8 pi G)) * [delta(a)^2 + 2 delta(a)]

         = (c^4 mu^2 / (8 pi G)) * [(I_0/(2 mu^2 a^3))^2 + 2 I_0/(2 mu^2 a^3)]

         = (c^4 mu^2 / (8 pi G)) * [I_0^2/(4 mu^4 a^6) + I_0/(mu^2 a^3)]

         = (c^4 / (8 pi G)) * [I_0^2/(4 mu^2 a^6) + I_0/a^3]
```

Define w_tilde_0 = I_0/(4 mu^2) (following BS2024 notation). Then:

```
rho_K(a) = (c^4 I_0 / (8 pi G)) * [w_tilde_0/a^6 + 1/a^3]

         = (c^4 I_0 / (8 pi G a^3)) * [1 + w_tilde_0/a^3]
```

**Dimension check on I_0:**
From K_Q = I_0/a^3 and [K_Q] = [m^{-2}] (K_Q = dK/dQ, K has dimensions m^{-2}, Q dimensionless):
```
[I_0] = [m^{-2}]    (since a is dimensionless)
```

Check rho_K dimensions:
```
[c^4 I_0 / G] = [m^4 s^{-4}] * [m^{-2}] / [m^3 kg^{-1} s^{-2}]
              = [m^2 s^{-4}] * [kg s^2 / m^3]
              = [kg m^{-1} s^{-2}]   CHECK
```

### 3.4 The Two Components

rho_K has two pieces:
1. **Dust-like** (a^{-3}): rho_dust = c^4 I_0 / (8 pi G a^3)
2. **Stiff-like** (a^{-6}): rho_stiff = c^4 I_0 w_tilde_0 / (8 pi G a^6)

The ratio at scale factor a:
```
rho_stiff / rho_dust = w_tilde_0 / a^3
```

At a = 1 (today): ratio = w_tilde_0
At recombination (a = 1/1101): ratio = w_tilde_0 * 1101^3 = 1.33 x 10^9 * w_tilde_0

**For CDM-like behavior at recombination, we need:**
```
w_tilde_0 * 1101^3 << 1
=> w_tilde_0 << 7.5 x 10^{-10}
```

### 3.5 What is w_tilde_0 for mu = H_0/c?

At a = 1 (today):
```
delta_0 = I_0 / (2 mu^2)
```

And:
```
w_tilde_0 = I_0 / (4 mu^2) = delta_0 / 2
```

For delta_0 = 0.340 (giving Omega_DM = 0.265):
```
w_tilde_0 = 0.340 / 2 = 0.170
```

For delta_0 = 0.395 (giving Omega_m = 0.315):
```
w_tilde_0 = 0.395 / 2 = 0.198
```

**PROBLEM: w_tilde_0 ~ 0.17 - 0.20 is NOT << 1.**

At recombination:
```
rho_stiff/rho_dust = 0.17 * (1101)^3 = 0.17 * 1.33 x 10^9 = 2.3 x 10^8
```

**This is catastrophic.** The stiff-matter component (a^{-6}) would COMPLETELY dominate over the dust component (a^{-3}) at recombination. The Khronon field would NOT behave like CDM.

---

## 4. The O(1) Prefactor Problem

### 4.1 Route 1's Claim Revisited

Route 1 claimed rho_eff = mu^2 c^2 / (8 pi G) = rho_crit/3, then said the ratio 0.265/0.333 = 0.795 gives (Q_0-1)^2 = 0.795, hence Q_0-1 = 0.892.

This used the **incorrect** quadratic-only approximation. With the full formula:

```
Omega_K = (delta^2 + 2 delta)/3
```

Setting Omega_K = 0.265: delta = 0.340 (not 0.892).

**The correct picture**: The dominant contribution at small delta is the LINEAR term 2*delta/3, not the quadratic delta^2/3. Route 1 was confused about which term dominates.

### 4.2 Can (Q_0 - 1) Be Derived?

The question is whether delta_0 = Q_0 - 1 = 0.340 can be derived from a physical principle.

**Answer: Almost certainly not.**

Reason: delta(a) = I_0/(2 mu^2 a^3) evolves with time. Its value today (delta_0) is determined by the initial condition I_0, which is set in the early universe. This is analogous to how the baryon density Omega_b is determined by baryogenesis -- it is a contingent fact about our universe's history, not derivable from fundamental constants.

### 4.3 Counting Free Parameters

With mu = H_0/c (fixed by dimensional analysis):

| Parameter | Status | Determines |
|-----------|--------|-----------|
| mu = H_0/c | FIXED (dimensional analysis) | Scale of Khronon kinetic term |
| I_0 | FREE (initial condition) | Omega_DM (or Omega_m) |
| delta_0 = I_0/(2 mu^2) | FREE (equivalent to I_0) | Same as above |

**One free parameter remains: I_0 (or equivalently delta_0 or Omega_DM).**

This is exactly the same number of free parameters as LCDM (which has Omega_DM h^2 as a free parameter fit to data).

---

## 5. Scaling with Redshift: Detailed Evolution

### 5.1 Complete rho_K(z) with mu = H_0/c

From Section 3.3, with mu = H_0/c:

```
rho_K(z) = rho_K,0 * (1+z)^3 * [1 + w_tilde_0 (1+z)^3] / [1 + w_tilde_0]
```

where rho_K,0 = rho_K at z=0.

**Effective equation of state:**

At late times (w_tilde_0 << a^3, i.e., z << z_stiff): w_eff ~ 0 (dust)
At early times (w_tilde_0 >> a^3, i.e., z >> z_stiff): w_eff ~ 1 (stiff matter)

The stiff-to-dust transition redshift:
```
z_stiff = w_tilde_0^{-1/3} - 1
```

For w_tilde_0 = 0.17: z_stiff = (0.17)^{-1/3} - 1 = 1.80 - 1 = 0.80

**The transition happens at z ~ 0.8, which is AFTER recombination (z ~ 1100)!**

This means for ALL of the pre-recombination era, the Khronon field behaves as STIFF MATTER (w = 1, rho ~ a^{-6}), not as dust.

### 5.2 Implications for the CMB

For the CMB acoustic peaks to be reproduced, the dark matter must:
1. Be non-relativistic (w ~ 0) before recombination
2. Have perturbations with c_s = 0 (no pressure support)
3. Have density ~ Omega_DM * rho_crit at z ~ 1100

With mu = H_0/c, condition (1) is CATASTROPHICALLY VIOLATED. The Khronon energy density at recombination would be dominated by the a^{-6} component:

```
rho_K(z=1100) ~ rho_K,0 * w_tilde_0 * (1101)^6 / (1 + w_tilde_0)
              ~ 0.17 * rho_K,0 * 1.76 x 10^{18}
              ~ 3 x 10^{17} * rho_K,0
```

This is ~10^{17} times larger than the expected CDM density at recombination (which should be ~ 1101^3 * rho_DM,0 ~ 10^9 * rho_DM,0). Off by 8 orders of magnitude.

### 5.3 The Dust-Like Requirement Sets a LOWER Bound on mu

For CDM-like behavior at z = 1100:
```
w_tilde_0 < 1/(1101)^3 = 7.5 x 10^{-10}
```

Since w_tilde_0 = I_0/(4 mu^2) and I_0 = 2 mu^2 delta_0:
```
w_tilde_0 = delta_0/2
```

So we need delta_0 < 1.5 x 10^{-9}. But then:
```
Omega_K = (delta_0^2 + 2 delta_0)/3 ~ 2 delta_0/3 < 10^{-9}
```

**This is 9 orders of magnitude below the observed Omega_DM = 0.265!**

### 5.4 The Fundamental Tension

With mu = H_0/c:
- **To get Omega_DM ~ 0.26**: need delta_0 ~ 0.34, giving w_tilde_0 ~ 0.17 (NOT dust-like)
- **To be dust-like at CMB**: need delta_0 < 10^{-9}, giving Omega_K < 10^{-9} (negligible)

**These two requirements are INCOMPATIBLE for mu = H_0/c.**

This is the same conclusion reached in mu_route2_a0_crooks.md (Section 9.1) and confirms that mu = H_0/c is phenomenologically UNACCEPTABLE despite being dimensionally natural.

---

## 6. Comparison with Skordis-Zlosnik (2021)

### 6.1 The Skordis-Zlosnik Theory

Skordis & Zlosnik (2021, PRL 127, 161302; arXiv:2007.00082) constructed the first relativistic MOND theory that fits the CMB. Their theory (AeST) uses:

- An Einstein-aether vector field A^mu (unit timelike)
- A scalar field phi (with shift symmetry)
- The metric g_{mu nu}

The Khronon limit (where A^mu is hypersurface-orthogonal) gives:

```
S_AeST = (c^3/16piG) integral sqrt(-g) [R + F(Y,Q) + lambda(g_{mu nu} A^mu A^nu + c^2)] d^4x + S_phi + S_m
```

where F(Y,Q) is a free function of:
- Y = A^mu;nu A_{mu;nu} (acceleration squared)
- Q = c * A^mu nabla_mu phi (related to Khronon Q in the hypersurface-orthogonal limit)

### 6.2 What Determines Omega_DM in AeST?

In AeST/Khronon:
1. **mu** (the mass scale in K(Q)) is a FREE PARAMETER
2. **I_0** (the Khronon "momentum" integration constant) is a FREE PARAMETER
3. **Omega_DM is determined by I_0**: Omega_DM = I_0 c^2 / (3 H_0^2)

**Skordis & Zlosnik do not predict Omega_DM.** They fit it to the CMB data.

Their fitted parameters (from arXiv:2007.00082, Table I):
- They work in units where 8piG = c = 1
- The scalar field coupling and mass scale are adjusted to fit Planck CMB TT, TE, EE spectra
- The effective dark matter density is an input (fit parameter)

### 6.3 Their Value of mu

Skordis & Zlosnik use a kinetic function of the form:
```
F(Q) = -2 alpha_K Q^2 + 2 mu_K^2 (Q - beta_K)^2
```

The mass scale mu_K is analogous to mu in Blanchet-Skordis. From their CMB fits:
- mu_K is chosen such that w_tilde_0 << 1 at recombination
- This requires mu_K >> H_0/c (by factors of ~10^4)
- The exact value is degenerate with other parameters in the fit

**Their mu is NOT H_0/c.** It is much larger, chosen phenomenologically.

### 6.4 Key Difference from tau Framework Claim

| Aspect | Skordis-Zlosnik | tau Framework (mu = H_0/c) |
|--------|----------------|---------------------------|
| mu | Free, fit to CMB (~10^4 H_0/c) | Claimed unique (H_0/c) |
| Omega_DM | Free (fit parameter) | Free (I_0 initial condition) |
| CDM-like at CMB? | YES (mu large enough) | NO (w_tilde_0 ~ 0.2) |
| MOND at galaxies? | YES (J(Y) separately tuned) | YES (claimed from a_0 = cH_0/2pi) |
| Predictive? | No (multiple free parameters) | No (I_0 free; mu fails CMB) |

---

## 7. The Key Question: Can Omega_DM h^2 Be Predicted?

### 7.1 Assessment

With mu_0 = H_0/c (derived) and delta_0 = Q_0 - 1 as the ONLY remaining free parameter:

**If delta_0 could be derived -> Omega_DM h^2 predicted (zero free parameters).**

But:
1. delta_0 is an initial condition from the scalar field equation K_Q = I_0/a^3
2. No known physical principle fixes I_0
3. Even if we could fix delta_0, mu = H_0/c gives w_tilde_0 ~ delta_0/2, which must be << 10^{-9} for CMB consistency -- making Omega_K negligibly small
4. **There is no value of delta_0 that simultaneously gives Omega_DM ~ 0.26 AND w_tilde_0 << 1 for mu = H_0/c**

### 7.2 The Catch-22

The problem is a fundamental tension:
```
Omega_K = (delta_0^2 + 2 delta_0)/3    (want ~ 0.26)
w_tilde_0 = delta_0/2                  (want << 10^{-9})
```

These two requirements demand:
```
delta_0 ~ 0.34    (from Omega_K)
delta_0 << 10^{-9} (from w_tilde_0)
```

**INCOMPATIBLE.** The only escape is to make mu >> H_0/c, so that:
```
w_tilde_0 = I_0/(4 mu^2) << 1
```
while keeping I_0 = 3 H_0^2 Omega_DM / c^2 fixed by observations.

This gives:
```
mu^2 >> I_0/4 = 3 H_0^2 Omega_DM / (4 c^2)
mu >> H_0/c * sqrt(3 Omega_DM / 4) ~ 0.45 H_0/c
```

For w_tilde_0 < 10^{-9}: mu > H_0/c * sqrt(Omega_DM/(4 * 10^{-9}/3)) ~ 1.4 x 10^4 * H_0/c

This is precisely the narrow window found in mu_route2_a0_crooks.md: mu ~ 10^4 H_0/c.

### 7.3 Final Verdict

```
=================================================================
  Omega_DM h^2 = 0.12 CANNOT be predicted from mu_0 = H_0/c.

  Reasons:
  1. Omega_DM is controlled by I_0 (initial condition), not mu
  2. mu = H_0/c gives w_tilde_0 ~ O(1), violating CMB constraints
  3. There is no value of delta_0 that works for BOTH Omega and CMB
  4. The phenomenologically required mu is ~10^4 * H_0/c, not H_0/c

  THIS IS NOT A MAJOR RESULT.
  It is a NEGATIVE RESULT confirming Route 2's findings.
=================================================================
```

---

## 8. What IS Achieved

Despite the negative result on the prediction, this calculation establishes several useful results:

### 8.1 Corrected Formula

The correct expression for the Khronon energy density fraction is:
```
Omega_K = (delta^2 + 2 delta)/3
```
where the LINEAR term (2 delta/3) dominates for delta < 2. Route 1's quadratic approximation (delta^2/3) was incorrect for the relevant range delta ~ 0.3-1.

### 8.2 The Catch-22 Quantified

The fundamental tension between Omega_K and w_tilde_0 for mu = H_0/c is now precisely quantified:
- To get Omega_K = 0.265: delta_0 = 0.340, w_tilde_0 = 0.170
- To get w_tilde_0 < 7.5 x 10^{-10}: delta_0 < 1.5 x 10^{-9}, Omega_K < 10^{-9}
- Gap: 8 orders of magnitude in delta_0, 17 orders of magnitude in Omega_K

### 8.3 The Q_0 = sqrt(2) Curiosity

The algebraically clean fact that Q_0 = sqrt(2) gives Omega_K = 1/3 exactly is noted but has no known physical basis.

### 8.4 mu = H_0/c as Dimensionally Unique but Phenomenologically Excluded

Route 1's conclusion that mu = H_0/c is dimensionally unique remains valid. But the present calculation shows it is phenomenologically EXCLUDED by the CMB. This is the "Khronon hierarchy problem":

```
mu_natural = H_0/c                    (dimensionally unique)
mu_required ~ 10^4 H_0/c              (phenomenologically required)
Hierarchy ratio ~ 10^4                 (mild but real)
```

---

## 9. Possible Resolutions

### 9.1 Running mu(k) = k

If mu runs with scale (anomalous dimension eta = 1 from Paper 3):
```
mu(k) = k
```
then at the Hubble scale k ~ H_0/c, mu ~ H_0/c (dimensional analysis value). But at smaller scales (galaxy/CMB scales), mu would be much larger:
```
At CMB scales (k ~ 0.1 Mpc^{-1} ~ 3 x 10^{-23} m^{-1}): mu ~ 3 x 10^{-23} >> 7 x 10^{-27}
```

This would give mu ~ 4000 * H_0/c at CMB scales, close to the required ~10^4 factor!

**However**: This interpretation requires that mu literally runs with wavenumber k, which is NOT established. The eta = 1 running applies to Newton's constant G(r) in Paper 3, not directly to the Khronon mass mu. Extending it to mu is SPECULATIVE.

### 9.2 DBI Completion

The DBI form K_DBI = (2mu^2/lambda_D)[1 - sqrt(1 - lambda_D(Q-1)^2)] introduces a second parameter lambda_D that modifies the stiff/dust ratio. In the DBI case, the effective w_tilde_0 can be much smaller even with mu = H_0/c, IF lambda_D is appropriately tuned.

### 9.3 Two-Regime mu

Perhaps mu is not a constant but has two regimes:
- mu_cosmo ~ 10^4 H_0/c (sets CDM behavior at perturbation level)
- mu_MOND ~ H_0/c (sets the overall energy scale at background level)

This would require a mechanism by which mu effectively changes between the background and perturbation levels -- a kind of "environmental" mu.

---

## 10. Numerical Summary

### 10.1 Physical Constants Used

```
H_0   = 67.4 km/s/Mpc = 2.184 x 10^{-18} s^{-1}
c     = 2.998 x 10^8 m/s
G     = 6.674 x 10^{-11} m^3 kg^{-1} s^{-2}
h     = 0.674
```

### 10.2 Derived Quantities

```
mu_0  = H_0/c = 7.28 x 10^{-27} m^{-1}
rho_crit = 3 c^2 H_0^2 / (8 pi G) = 8.53 x 10^{-27} kg/m^3
Omega_DM = 0.265   (Planck 2018)
Omega_DM h^2 = 0.1200 +/- 0.0012
Omega_b = 0.049
Omega_m = 0.314
```

### 10.3 Key Results for mu = H_0/c

| Quantity | Formula | Value |
|----------|---------|-------|
| delta_0 for Omega_DM = 0.265 | Solve quadratic | 0.340 |
| delta_0 for Omega_m = 0.315 | Solve quadratic | 0.395 |
| w_tilde_0 for Omega_DM = 0.265 | delta_0/2 | 0.170 |
| w_tilde_0 for Omega_m = 0.315 | delta_0/2 | 0.198 |
| z_stiff for Omega_DM case | (0.170)^{-1/3} - 1 | 0.80 |
| rho_stiff/rho_dust at z=1100 | w_tilde_0 * 1101^3 | 2.3 x 10^8 |
| CDM-like at recombination? | w_tilde_0 * 1101^3 << 1? | **NO** |
| mu required for CDM at z=1100 | > 1.19 x 10^{-22} m^{-1} | ~10^4 * H_0/c |

---

## 11. Relationship to Paper 4 Text

### 11.1 Paper 4, Section "What OEE does not determine" (line 924-958)

Paper 4 correctly states:
> "Omega_DM h^2 ~ 0.12 is not predicted. The Khronon background energy density vanishes for K(1) = 0; the observed dark matter density is an input parameter (initial condition), exactly as Earth's mass is an input for F = ma."

This is confirmed and strengthened by the present calculation. Not only is Omega_DM an initial condition, but the mu = H_0/c value is INCOMPATIBLE with CDM behavior at the CMB.

### 11.2 Paper 4, Open Problem 5 (line 1417-1421)

Paper 4 states:
> "Background dark matter energy. On FRW, Sigma_grav = 0 and K(1) = 0, giving zero background Khronon energy. The observed Omega_DM h^2 ~ 0.12 requires a non-trivial background profile -- an input, not a prediction."

This is also confirmed. The non-trivial profile (Q_0 != 1, i.e., delta != 0) is required, and its magnitude is a free parameter.

---

## 12. Conclusions

### 12.1 The Computation

Starting from mu_0 = H_0/c:

```
rho_K = (c^2 H_0^2 / (8 pi G)) * (delta^2 + 2 delta)

Omega_K = (delta^2 + 2 delta) / 3

w_tilde_0 = delta_0 / 2
```

For Omega_DM = 0.265: delta_0 = 0.340, w_tilde_0 = 0.170.

### 12.2 The Problem

w_tilde_0 = 0.17 means the Khronon is NOT dust-like. The stiff-matter component dominates at z > 0.8, making the theory INCOMPATIBLE with the CMB.

### 12.3 The Bottom Line

**Omega_DM h^2 = 0.12 does NOT fall out naturally from mu = H_0/c.** It requires:
1. A value of mu much larger than H_0/c (~10^4 times)
2. An initial condition I_0 tuned to give the observed Omega_DM
3. Neither of these is predicted by the tau framework

The tau framework provides a compelling SCALE (rho ~ H_0^2/G ~ rho_crit) and a compelling INTERPRETATION (dark matter = observer's inaccessible information) but does NOT predict the precise numerical value Omega_DM h^2 = 0.12.

This is consistent with the framework's self-described status as an "organizational principle" (like thermodynamics) rather than a "predictive theory" (like QCD).

---

## References

1. Blanchet, L. & Skordis, C. (2024). JCAP 11, 040. arXiv:2404.06584
2. Blanchet, L. & Skordis, C. (2025). arXiv:2507.00912
3. Skordis, C. & Zlosnik, T. (2021). PRL 127, 161302. arXiv:2007.00082
4. Planck Collaboration (2018). arXiv:1807.06209
5. mu_route1_dimensional.md (this repository)
6. mu_route2_a0_crooks.md (this repository)
7. paper4_JRSWW_omega_DM.md (this repository)

---

*Last updated: 2026-03-12*
*This document provides a complete, dimension-checked computation of Omega_DM from mu = H_0/c, yielding a NEGATIVE RESULT on the prediction of Omega_DM h^2 = 0.12.*
