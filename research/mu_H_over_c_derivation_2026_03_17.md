# Can mu_0 = H_0/c Be Derived from First Principles?

**Author**: Sheng-Kai Huang (with computational analysis)
**Date**: 2026-03-17
**Status**: COMPREHENSIVE ANALYSIS -- five approaches evaluated, verdict delivered
**Purpose**: Determine the epistemic status of mu_0 = H_0/c: derivable, or irreducible assumption?

---

## Executive Summary

We investigate five approaches to deriving the Khronon coupling constant mu_0 = H_0/c from first principles within the tau/Petz framework. The result is a **stratified picture**:

| Approach | Method | Result | Classification |
|----------|--------|--------|----------------|
| 1 | Petz bound at cosmological horizon | mu ~ H/c from Sigma ~ O(1) inside horizon | **HEURISTIC** |
| 2 | Running consistency (RG boundary condition) | mu_0 = H_0/c is the IR boundary condition of mu(k) = k | **TAUTOLOGICAL** |
| 3 | Extremal retrodictability (joint optimization) | Does NOT independently fix mu_0; requires mu_0 = H_0/c as input | **NEGATIVE** |
| 4 | Scale invariance of on-shell Khronon action | Identifies mu = H/c as the unique scale-invariant choice | **SUGGESTIVE** |
| 5 | DPI maximal rate and eta_mu = 1 | DPI constrains eta >= 0; extensivity selects eta = 1; mu_0 = H_0/c is boundary condition, not output | **PARTIAL** |

**Overall Verdict**: mu_0 = H_0/c **cannot be derived** from the Petz framework alone. It is a **naturality selection** -- the unique choice without Planck-scale physics -- elevated by multiple convergent arguments to the level of a well-motivated postulate. The deepest available argument is the modular flow identification (Route 3C of the Petz optimality analysis): if the Khronon IS the modular flow direction, then mu = H at the Hubble scale follows from the KMS periodicity. But this identification itself is structural, not derived.

**The minimal assumptions** to get mu_0 = H_0/c are:
1. The Khronon has dimensions [length^{-1}] (from the action)
2. No Planck-scale physics enters (d = 0 in the dimensional analysis)
3. The only cosmological scale is H_0/c

These are the same assumptions needed for a_0 ~ cH_0 in MOND -- and like that relation, mu_0 = H_0/c may be a fundamental coincidence that requires a deeper theory to explain, or it may be an irreducible boundary condition of our universe.

---

## 1. Approach 1: From the Petz Bound Applied to the Cosmological Horizon

### 1.1 Setup

The de Sitter horizon has:
- Temperature: T_dS = hbar H/(2pi k_B c)
- Entropy: S_dS = pi c^3/(G hbar H^2) (Gibbons-Hawking)
- Radius: L_H = c/H

The Petz bound states: F >= exp(-Sigma/2).

For the cosmological Khronon channel crossing the de Sitter horizon, we parameterize the entropy production as:
```
Sigma_dS = mu * L_H = mu * (c/H)
```
The idea is that Sigma measures the "information cost" of crossing the horizon, which depends on the Khronon coupling mu times the horizon size.

### 1.2 The Argument

At the horizon, the channel is maximally lossy: tau -> 1 (complete information loss), so formally Sigma -> infinity and F -> 0.

But JUST INSIDE the horizon (at proper distance epsilon from the horizon), the entropy production is finite. The natural condition for the "onset" of horizon behavior is Sigma ~ O(1), because:
- Sigma << 1: the channel is nearly reversible (far from the horizon)
- Sigma >> 1: the channel is nearly irreversible (deep inside the horizon)
- Sigma ~ O(1): the transition regime

Setting Sigma_transition ~ 1:
```
mu * (c/H) ~ 1
mu ~ H/c
```

### 1.3 A More Precise Version

From Paper 2, the gravitational entropy production is Sigma_grav = -ln(-g_00). For the de Sitter metric in static coordinates:
```
g_00 = -(1 - H^2 r^2/c^2)
Sigma_grav(r) = -ln(1 - H^2 r^2/c^2)
```

At the horizon r = c/H: Sigma -> infinity (as expected).

The Khronon contribution to Sigma depends on Q - 1. On the de Sitter background with the running mu = H/c:
```
Sigma_Khronon = (Q-1)^2 * delta(z) * (2 + delta(z))
```

The total Sigma at the Hubble radius should satisfy the Petz bound: F >= exp(-Sigma_total/2). The condition that the Khronon maintains O(1) fidelity at the horizon (i.e., F not exponentially suppressed) requires Sigma_Khronon ~ O(1) at the Hubble scale.

From the FRW framework: delta_0 = 0.34, giving Sigma_Khronon = delta_0(2+delta_0) = 0.80. This IS O(1), consistent with mu_0 = H_0/c.

### 1.4 What This Proves and What It Doesn't

**What it shows**: If Sigma_Khronon ~ O(1) at the Hubble scale, then mu * (c/H) ~ O(1), hence mu ~ H/c. The O(1) Sigma condition is physically natural -- it means the observer can still partially recover information at the cosmological horizon, neither trivially (Sigma << 1) nor hopelessly (Sigma >> 1).

**What it doesn't show**: WHY Sigma should be O(1). This is an aesthetic argument, not a mathematical derivation. One could equally argue for Sigma = 2pi (the "thermal" value) or Sigma = ln 2 (the "bit" value). Each gives a different O(1) coefficient.

**The deeper problem**: The argument uses Sigma ~ O(1) as a PRINCIPLE, but does not derive this principle from the Petz framework. The Petz bound F >= exp(-Sigma/2) is satisfied for ANY Sigma >= 0; it does not prefer O(1) values.

### 1.5 Classification

```
STATUS: HEURISTIC
Strength: Gives the right answer to order of magnitude
Weakness: The O(1) Sigma condition is not derived from any principle
What's needed: A principle that selects Sigma ~ O(1) at the horizon
```

---

## 2. Approach 2: From the Running Consistency (IR Boundary Condition)

### 2.1 The Argument

If mu runs as mu(k) = k (from eta_mu = 1, supported by 5 convergent arguments -- see mu_running_proof.md), then evaluating at the IR cutoff:
```
mu(k = H_0/c) = H_0/c
```

This gives mu_0 = H_0/c "automatically" as the boundary condition of the running at the largest scale.

### 2.2 Why This Is Tautological

This argument is circular in a subtle way:

1. The running is mu(k) = mu_0 * (k/k_H)^{eta_mu} with k_H = H_0/c
2. For eta_mu = 1: mu(k) = mu_0 * (k c/H_0)
3. At k = H_0/c: mu(H_0/c) = mu_0 * 1 = mu_0

This is just evaluating the running at its own boundary condition. It tells us that the running IS CONSISTENT with mu_0 = H_0/c, but it does not DERIVE mu_0 = H_0/c.

The ACTUAL content of this approach is: **given** that mu runs as mu(k) = k (from eta_mu = 1), the boundary condition at the Hubble scale is mu(k_H) = k_H = H_0/c. But the running mu(k) = k implicitly USES the boundary condition mu_0 = H_0/c to eliminate the arbitrary normalization constant.

To see this explicitly: with eta_mu = 1, the most general running is:
```
mu(k) = A * k
```
where A is a dimensionless constant set by the boundary condition. The choice A = 1 (i.e., mu = k in natural units) is equivalent to choosing the boundary condition mu(k_H) = k_H = H_0/c.

### 2.3 What Would Make This Non-Tautological

If we could derive A = 1 from an independent principle, this approach would be non-trivial. The dimensional transmutation argument (Approach C in mu_running_proof.md) states that A = 1 is the unique value making mu(k) independent of the boundary condition scale. This is the closest thing to a non-trivial input.

But "independent of the boundary condition scale" is itself a CHOICE, not a derivation. Other values of A are dimensionally consistent (they just make the running explicitly depend on H_0).

### 2.4 The Fixed-Point Argument

The strongest version of this approach uses the fixed-point structure of the dimensionless coupling tilde_mu = mu/k:

1. If tilde_mu has an IR fixed point tilde_mu* = 1, then mu(k) -> k at large scales
2. At the Hubble scale: mu(k_H) = tilde_mu* * k_H = H_0/c

The question reduces to: why is tilde_mu* = 1 and not some other O(1) number?

From Approach D of mu_running_proof.md, the beta function is:
```
beta_{tilde_mu} = (eta_mu - 1) * tilde_mu + higher order
```

At the fixed point: beta = 0, so either tilde_mu* = 0 (trivial) or eta_mu = 1 (the non-trivial fixed point exists for any tilde_mu*). The specific value tilde_mu* is determined by the higher-order terms, which are unknown.

**Net result**: The fixed-point analysis constrains eta_mu = 1 but NOT the value of tilde_mu*. The boundary condition A = tilde_mu* = 1 remains an input.

### 2.5 Classification

```
STATUS: TAUTOLOGICAL (as a derivation of mu_0)
         STRONG EVIDENCE (for eta_mu = 1 and the running mu(k) = k)
The running is well-motivated; the boundary condition is not derived.
```

---

## 3. Approach 3: From the Extremal Retrodictability Principle (Joint Optimization)

### 3.1 The Idea

The Omega_DM derivation (omega_DM_derivation_2026_03_16.md) found delta_0 = 0.343 from maximizing the retrodictability functional R[delta_0]. This used mu_0 = H_0/c as a fixed input. Can we instead treat BOTH mu_0 and delta_0 as free and let the extremal principle fix both?

### 3.2 The Functional with Two Free Parameters

The retrodictability functional with general mu_0 (not necessarily H_0/c):
```
R[mu_0, delta_0] = integral_0^infinity F(z) * Sigma(z) * (dt/dz) dz
```

where:
- F(z) = 1/(1 + delta(z))
- Sigma(z) = delta(z) * (2 + delta(z))
- delta(z) = delta_0 / (E^2(z) * a^3)
- E^2(a) = Omega_r/a^4 + Omega_m/a^3 + Omega_Lambda
- Omega_K = mu_0^2 * delta_0 * (2 + delta_0) / (3 H_0^2/c^2)

Wait -- this is where the structure matters. Let me be more careful.

### 3.3 The Dependence on mu_0

The Khronon density parameter is:
```
Omega_K = c^2 mu_0^2 * delta_0 * (2 + delta_0) / (3 H_0^2)
```

For mu_0 = H_0/c, this reduces to Omega_K = delta_0(2+delta_0)/3 (the formula used in the Omega_DM derivation).

For general mu_0, define the dimensionless ratio:
```
x = c mu_0 / H_0
```

Then:
```
Omega_K = x^2 * delta_0 * (2 + delta_0) / 3
```

The conservation law gives:
```
delta(a) = delta_0 / (E^2(a) * a^3)
```
where we've assumed the running mu_bg(a) = H(a)/c * x (generalizing the running with the same eta_mu = 1 but different normalization).

Actually, this needs more care. If mu runs as mu(k) = x * k (with general normalization x), then:
```
mu_bg(a) = x * H(a)/c
delta(a) = I_0 / (2 mu^2(a) * a^3) = I_0 c^2 / (2 x^2 H^2(a) * a^3)
         = delta_0 * H_0^2 / (H^2(a) * a^3)
         = delta_0 / (E^2(a) * a^3)
```

The crucial point: delta(a) depends on delta_0 but NOT on x! This is because mu^2 * delta = I_0/(2a^3) depends only on I_0 and a, not on how mu is partitioned into "normalization x" and "running H/c."

Therefore:
```
Omega_K = delta_0(2 + delta_0) * x^2 / 3
Omega_m = Omega_b + delta_0(2 + delta_0) * x^2 / 3
```

And the Friedmann equation involves Omega_m, which depends on x.

### 3.4 The Two-Parameter Optimization

The retrodictability functional R[x, delta_0] depends on both parameters:

1. Through the integrand: F * Sigma depends on delta(z), which depends on delta_0 and E(z), and E(z) depends on Omega_m = Omega_b + delta_0(2+delta_0) x^2/3
2. Through the proper time weighting: dt/dz = 1/((1+z) H_0 E(z)), and E(z) depends on Omega_m(x, delta_0)

Setting dR/d(delta_0) = 0 and dR/dx = 0 simultaneously would fix both.

### 3.5 The Critical Analysis

**Here is the problem.** Let us examine how R depends on x at fixed delta_0:

The product F * Sigma = delta(2+delta)/(1+delta) depends on delta(z) = delta_0/(E^2 a^3). The expansion rate E(z) depends on Omega_m(x, delta_0). For fixed delta_0, increasing x increases Omega_m, which:
1. Makes the universe younger (less proper time for the integral)
2. Changes E(z), which changes delta(z)

But here is the key: Omega_K = (delta_0(2+delta_0)/3) * x^2, so x and delta_0 are DEGENERATE in their effect on the Friedmann equation. The only combination that enters the Friedmann equation is:
```
Omega_K = delta_0(2+delta_0)/3 * x^2
```

This means R is a function of (delta_0, Omega_K) or equivalently (delta_0, x^2 * delta_0(2+delta_0)), NOT independently of x and delta_0.

The extremal conditions are:
```
dR/d(delta_0) = 0  at fixed Omega_K: determines how delta_0 is partitioned (given Omega_K)
dR/d(Omega_K) = 0: determines the total dark matter density
```

But the FIRST condition is trivial! At fixed Omega_K, the integrand F * Sigma = delta(z)(2+delta(z))/(1+delta(z)) depends only on delta(z), which depends only on delta_0/(E^2 a^3). Since E^2 depends on Omega_K (fixed), F*Sigma depends only on delta_0. And the proper-time weighting dt/dz depends on E(z), which depends only on Omega_K. So R depends on delta_0 ONLY through the integrand, not through the weighting.

At fixed Omega_K, we need:
```
dR/d(delta_0)|_{Omega_K} = integral [d(F*Sigma)/d(delta_0)] * dt/dz dz = 0
```

But F * Sigma = (1+delta) - 1/(1+delta) is monotonically increasing in delta, and delta propto delta_0, so d(F*Sigma)/d(delta_0) > 0 everywhere.

**This means dR/d(delta_0) > 0 at fixed Omega_K -- the functional ALWAYS wants larger delta_0!**

The only thing preventing delta_0 -> infinity is that Omega_K = delta_0(2+delta_0)x^2/3 is fixed, which forces x -> 0 as delta_0 -> infinity.

**There is NO extremum in the (x, delta_0) plane at fixed Omega_K.** The retrodictability functional has a RIDGE in the (x, delta_0) plane along the curve x^2 delta_0(2+delta_0) = 3 Omega_K*, where Omega_K* is the optimal dark matter density. Along this ridge, R is constant (because only Omega_K matters for R, and it's fixed). The functional is flat in the direction along the ridge.

### 3.6 Physical Interpretation

The retrodictability extremal principle determines Omega_DM (the total dark matter energy density), NOT how it is partitioned into mu_0 and delta_0 separately. This makes physical sense: the recoverable entropy depends on how much dark matter there is (Omega_K), not on the internal parameterization (mu_0 vs delta_0).

Therefore:
```
R fixes: Omega_K = 0.268  (total dark matter density)
R does NOT fix: the partition into mu_0 and delta_0 separately
```

To fix mu_0, we need ADDITIONAL information beyond the retrodictability principle.

### 3.7 Could a Modified Functional Break the Degeneracy?

If R depended on mu_0 and delta_0 independently (not just through Omega_K), the degeneracy would be broken. This would happen if:

1. The running changes with the normalization x (but it doesn't -- the running is mu(k) = x*k, and delta(z) is x-independent)
2. The fidelity F depends on mu_0 directly (but F = 1/(1+delta), and delta doesn't depend on x)
3. There are perturbation effects beyond the background (e.g., the Jeans scale depends on mu_0 directly)

Option 3 is interesting: the Khronon sound speed c_s = 0, but the Jeans-like scale 1/mu sets the scale below which the Khronon behaves differently from CDM. If the retrodictability functional included a scale-dependent weighting (not just proper-time weighting), it might distinguish between different mu_0 values.

But this would require going beyond the FRW background, into the perturbation theory -- a significant extension of the current framework.

### 3.8 Classification

```
STATUS: NEGATIVE
The extremal retrodictability principle fixes Omega_DM but NOT mu_0.
The degeneracy is fundamental: Omega_K depends on mu_0^2 * delta_0(2+delta_0),
and R depends only on Omega_K.
mu_0 and delta_0 are individually undetermined by this principle.
```

---

## 4. Approach 4: Scale Invariance of the On-Shell Khronon Action

### 4.1 The Idea

On a de Sitter background (H = const, a(t) = e^{Ht}), compute the Khronon on-shell action with the running mu_bg(a) = H/c. If S_K is independent of H (scale-invariant), this would identify mu = H/c as special.

### 4.2 The Calculation

The Khronon action density is:
```
L_K = (c^4/(16 pi G)) * mu^2 * (Q_0 - 1)^2
```

On the de Sitter background, Q_0 is determined by the conservation law. At the de Sitter fixed point (constant H), the Khronon expansion scalar is:
```
Q_0 = 1/(N) * d(a^3)/dt / (3 a^3) = ...
```

Actually, let me be precise. Q is defined as Q = c sqrt(-g^{ab} d_a tau d_b tau) for the Khronon scalar tau. On the FRW background with tau = t:
```
Q_0 = c / sqrt(-g_00) = c / c = 1
```
(in coordinates where g_00 = -c^2).

Wait -- this gives Q_0 = 1 identically on ANY FRW background. So Q_0 - 1 = 0 and the on-shell action vanishes!

### 4.3 The Non-Trivial Configuration

The non-trivial Khronon configuration has Q_0 = 1 + delta, where delta != 0. This arises from a Khronon that is NOT aligned with the cosmic time: tau != t. The conservation law K'(Q_0) = I_0/(a^3 Q_0) gives:
```
delta = I_0 / (2 mu^2 a^3)
```

The on-shell action per unit comoving volume is:
```
S_K = integral dt a^3 * (c^4/(16 pi G)) * mu^2 * delta^2
    = integral dt a^3 * (c^4/(16 pi G)) * mu^2 * I_0^2 / (4 mu^4 a^6)
    = integral dt (c^4 I_0^2) / (64 pi G mu^2 a^3)
```

With mu = H/c and H = const (de Sitter), a = e^{Ht}:
```
S_K = (c^4 I_0^2) / (64 pi G (H/c)^2) * integral dt e^{-3Ht}
    = (c^6 I_0^2) / (64 pi G H^2) * [1/(3H)]
    = (c^6 I_0^2) / (192 pi G H^3)
```

This depends on H! Specifically, S_K propto 1/H^3 propto Volume of de Sitter space (which scales as (c/H)^3).

### 4.4 The Energy Density Route

The Khronon energy density is:
```
rho_K = (c^2 mu^2 / (8 pi G)) * delta * (2 + delta)
```

With mu = H/c and delta = delta_0 (constant in the matter era, where H^2 propto a^{-3}):
```
rho_K = (H^2 / (8 pi G)) * delta_0 * (2 + delta_0)
      = rho_crit * delta_0(2+delta_0) / 3
```

This is proportional to rho_crit at ALL epochs (because mu = H/c tracks the Hubble rate). This is a form of "tracking behavior" -- the Khronon energy density maintains a constant fraction of the critical density.

The interesting observation: with mu = H/c, the ratio rho_K/rho_crit = delta_0(2+delta_0)/3 is a CONSTANT (independent of time, of H, of the expansion history). Any other choice of mu_0 would give a rho_K/rho_crit that evolves with time.

### 4.5 The Tracking Condition as a Selection Principle

**Principle**: The Khronon coupling mu_0 should be chosen such that the Khronon energy density tracks the critical density at all epochs.

This requires:
```
rho_K / rho_crit = const  =>  mu^2 * delta / rho_crit = const
```

Since rho_crit = 3H^2/(8 pi G) and mu^2 delta = I_0 c^2/(2 a^3) (from the conservation law), this gives:
```
I_0 c^2 / (2 a^3) * (2+delta) / (3 H^2/(8 pi G)) = const
```

In the matter era (delta ~ const, H^2 propto a^{-3}):
```
I_0 c^2 * (2+delta) * 8 pi G / (2 * 3 * H_0^2 * Omega_m) = const
```

This is automatically satisfied (both numerator and denominator scale as a^{-3}). The tracking condition holds for ANY mu_0 when the running is mu_bg = H/c, as long as delta is approximately constant.

Wait -- this is because the conservation law mu^2 delta = I_0/(2 a^3) and rho_crit propto H^2 propto a^{-3} (matter era) both scale as a^{-3}. The tracking is a CONSEQUENCE of the running, not a selection principle.

### 4.6 A Better Version: Scale Invariance of the Dimensionless Action

Define the dimensionless Khronon action:
```
tilde_S_K = (8 pi G / c^4) * S_K
```

On de Sitter with mu = x * H/c (general normalization):
```
tilde_S_K propto x^{-2} * (I_0/a^3)^2 / H^3 * integral
```

For mu = H/c specifically (x = 1), using the on-shell delta:
```
tilde_S_K propto I_0^2 / H^3 * integral dt e^{-3Ht}
            propto I_0^2 / H^4
```

The dimensionless action scales as 1/H^4. For de Sitter, H is the only scale, so this is just S_dS^{-1} ~ G H^2 ~ 1/S_Gibbons-Hawking. This is NOT scale-invariant.

### 4.7 Alternative: The Action Ratio

Consider the ratio of the Khronon action to the Einstein-Hilbert action:
```
S_K / S_EH = mu^2 * delta^2 / R
```

where R = 12 H^2/c^2 (de Sitter Ricci scalar). With mu = H/c:
```
S_K / S_EH = (H/c)^2 * delta^2 / (12 H^2/c^2) = delta^2/12
```

This ratio IS scale-invariant (independent of H)! For delta_0 = 0.343: S_K/S_EH = 0.010.

For general mu_0 = x H/c:
```
S_K / S_EH = x^2 delta^2 / 12
```

Since delta propto 1/(mu^2 a^3) propto 1/(x^2 H^2 a^3), and delta is evaluated at a = 1:
```
delta_0 = I_0 c^2 / (2 x^2 H_0^2)
S_K / S_EH = x^2 * I_0^2 c^4 / (4 x^4 H_0^4 * 12) = I_0^2 c^4 / (48 x^2 H_0^4)
```

This DOES depend on x! The ratio is minimized as x -> infinity (maximum coupling = minimum delta). There is no natural selection of x = 1 from this ratio.

### 4.8 Classification

```
STATUS: SUGGESTIVE (but does not uniquely derive mu_0 = H_0/c)

The tracking behavior (rho_K/rho_crit = const) is a nice feature of
mu = H/c, but it holds automatically for any normalization when the
running has eta_mu = 1. The scale invariance of S_K/S_EH requires x = 1
only if we additionally impose delta = const (which itself requires the
specific initial condition I_0 that gives Omega_K = 0.268).

No independent principle has been found that uniquely selects x = 1
from the on-shell action analysis.
```

---

## 5. Approach 5: From DPI (Data Processing Inequality) and eta_mu = 1

### 5.1 Review: What DPI Actually Constrains

The DPI states that for any CPTP map N: D(N(rho) || N(sigma)) <= D(rho || sigma).

Applied to the gravitational RG flow (coarse-graining from scale k to k' < k):
- D(rho_k || sigma_k) >= D(rho_{k'} || sigma_{k'})
- This implies Sigma(k) is non-increasing under coarse-graining
- For the Khronon: mu^2(k)/k^2 (the dimensionless coupling) must be non-decreasing toward the IR
- This gives eta_mu >= 0 (mu grows toward larger scales)

**DPI alone does not fix eta_mu = 1.** It only gives a lower bound eta_mu >= 0.

### 5.2 The Extensivity Argument

The specific value eta_mu = 1 comes from the extensivity assumption:

1. At scales r >> c/H_0, entanglement entropy transitions from area-law to volume-law (Verlinde 2016, Casini-Huerta 2017)
2. The Khronon-mediated force correction in position space scales as r^{-(2+eta_mu)}
3. The ratio to Newtonian force: delta_F/F_N ~ r^{-eta_mu}
4. For eta_mu > 1: force correction grows without bound (IR divergent, violates de Sitter bound)
5. For eta_mu < 1: force correction decays (IR irrelevant, contradicts extensivity)
6. For eta_mu = 1: marginal (logarithmic correction to potential)

This gives eta_mu = d - 3 = 1 in d = 4.

### 5.3 Does the DPI Maximal Rate Argument Work?

The idea: perhaps the DPI selects not just eta >= 0 but specifically eta = 1 as the "maximal DPI rate" -- the fastest running consistent with the data processing inequality.

**Analysis**: The DPI constrains:
```
d Sigma(k) / d ln k <= 0
```
(entropy production cannot increase under coarse-graining).

For the Khronon: Sigma(k) ~ mu^2(k)/k^2 ~ k^{2(eta_mu - 1)}. The DPI rate is:
```
|d Sigma / d ln k| = 2|eta_mu - 1| * Sigma(k)
```

For eta_mu > 1: Sigma increases toward UV, which is allowed by DPI (coarse-graining reduces it).
For eta_mu < 1: Sigma increases toward IR, which VIOLATES the DPI unless the magnitude is bounded.
For eta_mu = 1: Sigma is scale-independent (logarithmic), and |d Sigma/d ln k| = 0 -- the SLOWEST allowed rate!

**Wait -- this gives eta_mu = 1 as the SLOWEST rate, not the fastest!** The DPI is most easily satisfied when Sigma doesn't change much (eta_mu = 1), and is hardest to satisfy when eta_mu << 1 (where Sigma grows rapidly toward the IR).

So the DPI selects eta_mu >= 1 (not = 1) for strict satisfaction at all scales. The extensivity argument then provides the upper bound eta_mu <= 1 (from the de Sitter entropy bound). Together:

```
DPI: eta_mu >= 1  (for the Khronon channel to be a valid CPTP map)
Extensivity: eta_mu <= 1  (for the entropy not to exceed the de Sitter bound)
Together: eta_mu = 1  (uniquely!)
```

### 5.4 Can This Fix mu_0?

Even with eta_mu = 1 established, the DPI does not fix the boundary condition. The DPI constrains the RATE of running but not the NORMALIZATION. The boundary condition mu_0 = mu(k_H) remains a free parameter.

However, if we combine the DPI with the modular flow identification (Route 3C of mu_route3_petz_optimality.md), there is a path:

1. The modular flow at the Hubble scale has KMS periodicity beta = 2pi/H
2. The modular temperature is T_mod = H/(2pi)
3. If the Khronon coupling is the modular temperature: mu = 2pi T_mod = H
4. At the Hubble scale: mu(k_H) = H_0 (in natural units where c = 1), i.e., mu_0 = H_0/c

This gives mu_0 = H_0/c from the modular flow normalization. But the identification "Khronon coupling = modular temperature" is itself a POSTULATE (the OEE postulate), not derived from the DPI.

### 5.5 The Hierarchy of Assumptions

To go from the DPI to mu_0 = H_0/c, one needs:

**Level 1 (rigorous)**: DPI gives eta_mu >= 0
**Level 2 (semi-rigorous)**: DPI + extensivity gives eta_mu = 1
**Level 3 (speculative)**: Modular flow identification gives mu(k_H) = H_0/c

Each level adds assumptions:
- Level 1: CPTP (standard)
- Level 2: de Sitter extensivity (well-motivated but not proved)
- Level 3: Khronon = modular flow direction (structural identification, not derived)

### 5.6 Classification

```
STATUS: PARTIAL
DPI constrains the running (eta_mu >= 0, and = 1 with extensivity).
DPI does NOT fix the normalization (mu_0).
The normalization requires the modular flow identification (speculative).
```

---

## 6. Synthesis: The Derivation Landscape

### 6.1 What Has Been Achieved

```
LEVEL 1: PROVED (mathematical facts)
  - [mu] = [length^{-1}] from the Khronon action
  - mu_0 = H_0/c is the UNIQUE choice without Planck-scale physics
  - mu_0 = 2pi k_B T_dS/(hbar c) (algebraic identity)
  - eta_mu = 1 at any non-trivial fixed point of tilde_mu = mu/k

LEVEL 2: STRONGLY MOTIVATED (convergent evidence)
  - eta_mu = 1 from DPI + extensivity (3 independent approaches)
  - mu(k) = k bridges the hierarchy (10^4 gap between H_0/c and kpc^{-1})
  - Modular flow gives mu ~ H at the Hubble scale
  - Tracking behavior: rho_K/rho_crit = const with mu = H/c

LEVEL 3: NOT ACHIEVED
  - Derivation of mu_0 = H_0/c from the Petz bound or DPI alone
  - Independent principle that selects the O(1) coefficient (why mu_0 = H_0/c
    and not, say, 2pi H_0/c or H_0/(2pi c))
  - Joint determination of mu_0 and delta_0 from a single variational principle
```

### 6.2 The Honest Answer

**mu_0 = H_0/c is NOT derivable from the Petz framework with current tools.**

It is the unique naturality choice, supported by multiple convergent arguments, and elevated to a well-motivated postulate. The closest thing to a "derivation" is the chain:

```
Gibbons-Hawking temperature T_dS
  -> modular flow periodicity beta = 2pi/(H/c) [in length units]
  -> Khronon coupling = inverse periodicity: mu = H/c
```

But this requires identifying the Khronon with the modular flow direction, which is the OEE postulate -- a structural identification, not a mathematical theorem.

### 6.3 Comparison with Similar Situations in Physics

| Situation | Status | Analogy |
|-----------|--------|---------|
| a_0 ~ cH_0 in MOND | Observed coincidence, not derived | Closest analogy to mu_0 = H_0/c |
| Lambda ~ H_0^2 (cosmological constant) | The cosmological constant problem | mu_0 = H_0/c avoids this by NOT involving G or hbar |
| m_e/m_p ~ 1/1836 | Not derived from QED | mu_0/k_Pl ~ 10^{-61}: WOULD be such a problem if Planck physics entered |
| alpha_EM ~ 1/137 | Not derived from QED | But derived from GUT unification conditions |

The most relevant comparison is the MOND a_0 ~ cH_0 coincidence:
- Both relate a "local" parameter (mu_0, a_0) to a cosmological scale (H_0)
- Both are "naturality" arguments, not derivations
- Both are 40+ years old (MOND since 1983) with no deeper explanation
- The tau framework CONNECTS them: a_0 = c^2 mu_0/(2pi) = cH_0/(2pi)

If mu_0 = H_0/c could be derived, it would simultaneously explain the MOND coincidence -- a 40-year-old puzzle.

### 6.4 The Minimal Assumption Set

To use mu_0 = H_0/c in the theory, one needs EXACTLY ONE of:

**(A) Dimensional analysis**: [mu] = m^{-1}, no Planck-scale physics, unique choice -> mu_0 = H_0/c
**(B) De Sitter thermodynamics**: Khronon Compton wavelength = Hubble radius -> mu_0 = H_0/c
**(C) Modular flow**: Khronon = modular flow, KMS periodicity = 2pi c/H -> mu_0 = H_0/c
**(D) Running boundary condition**: mu(k) = k, evaluated at k_H = H_0/c -> mu_0 = H_0/c

All four are equivalent statements of the same physical content. None is more "fundamental" than the others -- they are four faces of the same diamond.

### 6.5 What Would Count as a Genuine Derivation

A genuine derivation would require showing that the Petz recovery framework + the Khronon action + some well-established principle REQUIRES mu = H/c, in the sense that any other value leads to a mathematical inconsistency or physical pathology.

**The closest attempt**: The DPI + extensivity argument shows that eta_mu != 1 leads to pathologies (IR divergence or extensivity violation). But this constrains the RUNNING, not the NORMALIZATION. The normalization mu_0 = H_0/c is consistent but not required.

**What's missing**: A principle that fixes the normalization. Possibilities include:
1. A microscopic calculation of the Khronon self-energy that gives A = 1 in mu(k) = A*k
2. A consistency condition from the full quantum gravity theory that requires tilde_mu* = 1 at the IR fixed point
3. A holographic argument relating the bulk Khronon coupling to the boundary conformal dimension
4. A no-go theorem showing that any other normalization leads to inconsistency with observations (this would be phenomenological, not first-principles)

None of these exist currently.

---

## 7. Connection to the Omega_DM Extremal Principle

### 7.1 What the Extremal Principle Assumes

The Omega_DM = 0.268 derivation (omega_DM_derivation_2026_03_16.md) uses:
1. F = 1/(1+delta) -- canonical tau framework choice
2. Sigma = delta(2+delta) -- exact Khronon stress-energy
3. Proper-time weighting
4. **mu_0 = H_0/c** -- assumed, not derived
5. **Running mu_bg = H(a)/c** -- assumed, not derived

### 7.2 What It Derives

Given assumptions 1-5, the extremal principle derives:
- delta_0 = 0.343
- Omega_DM = 0.268 (1.1% from observation)

### 7.3 What It Cannot Derive

As shown in Approach 3 above, the extremal principle CANNOT simultaneously derive mu_0 and delta_0. The functional R depends on Omega_K = mu_0^2 delta_0(2+delta_0)/(3H_0^2/c^2), not on mu_0 and delta_0 separately. The degeneracy is exact.

### 7.4 The Updated Chain of Logic

The correct derivation chain is:
```
Step 1: mu_0 = H_0/c           [ASSUMED -- naturality / modular flow]
Step 2: eta_mu = 1               [DERIVED from DPI + extensivity]
Step 3: mu_bg(a) = H(a)/c        [FOLLOWS from Steps 1 + 2]
Step 4: Omega_K = delta_0(2+delta_0)/3  [FOLLOWS from Step 1]
Step 5: delta_0 = 0.343          [DERIVED from retrodictability extremum]
Step 6: Omega_DM = 0.268         [FOLLOWS from Steps 4 + 5]
```

The input parameters are:
- mu_0 = H_0/c (1 assumption -- naturality)
- eta_mu = 1 (1 derivation -- DPI + extensivity)
- Omega_b = 0.049 (1 measured input)

The output is: Omega_DM = 0.268 (1.1% from observation).

### 7.5 Could a Future Derivation of mu_0 Complete the Chain?

If mu_0 = H_0/c could be derived from first principles, the derivation chain would have:
- 0 free dark matter parameters
- 1 measured input (Omega_b)
- 1 output (Omega_DM = 0.268)

This would be a genuine prediction of the dark matter density from information-theoretic principles plus the baryon density. It would be comparable to (but more powerful than) Verlinde's (2016) order-of-magnitude prediction Omega_DM ~ rho_crit/3.

---

## 8. Honest Classification Table

| Claim | Status | What's needed to upgrade |
|-------|--------|-------------------------|
| [mu] = m^{-1} from the action | **PROVED** | Nothing |
| mu_0 = H_0/c is unique without Planck physics | **PROVED** (dimensional analysis) | Nothing |
| mu_0 = 2pi k_B T_dS/(hbar c) | **PROVED** (algebraic identity) | Nothing |
| eta_mu = 1 from DPI + extensivity | **DERIVED** (semi-rigorous) | Microscopic calculation of Khronon self-energy |
| mu(k) = k bridges the hierarchy | **DERIVED** (numerical check) | Full Boltzmann code verification |
| mu_0 = H_0/c from Petz bound at horizon | **HEURISTIC** (Sigma ~ O(1) argument) | Principle selecting O(1) Sigma |
| mu_0 = H_0/c from running boundary condition | **TAUTOLOGICAL** | Derivation of A = 1 in mu(k) = Ak |
| mu_0 = H_0/c from extremal retrodictability | **NEGATIVE** (degeneracy in mu_0 vs delta_0) | New functional that breaks degeneracy |
| mu_0 = H_0/c from on-shell action | **SUGGESTIVE** (tracking behavior) | Principle selecting tracking |
| mu_0 = H_0/c from DPI maximal rate | **PARTIAL** (constrains running, not normalization) | Modular flow derivation |
| mu_0 = H_0/c from modular flow | **SPECULATIVE** (structural identification) | Proof that Khronon = modular flow |
| mu_0 = H_0/c is derivable (overall) | **NOT ACHIEVED** | Any of the above upgrades |

---

## 9. Recommendations

### 9.1 For the Papers

**Paper 3** (weak field): Continue to state mu_0 = H_0/c as "the unique naturality selection without Planck-scale physics" with the label [ASSUMED]. The 5 convergent arguments for eta_mu = 1 elevate this from "arbitrary" to "well-motivated," but it remains an assumption.

**Paper 4** (unification): The honest framing is:
> "The Khronon coupling mu_0 = H_0/c is the unique selection from dimensional analysis without Planck-scale physics. This is a naturality argument -- the same epistemic status as a_0 ~ cH_0 in MOND. Unlike MOND, the tau framework provides a structural explanation (de Sitter thermodynamics + modular flow), but not a mathematical derivation. A rigorous derivation from first principles remains the most important open problem of the framework."

### 9.2 For Future Research

**Priority 1**: Derive A = 1 in mu(k) = Ak from a microscopic calculation. This is equivalent to deriving mu_0 = H_0/c. The most promising route is the Khronon self-energy at non-perturbative level, using the Casini-Huerta modular Hamiltonian.

**Priority 2**: Find a variational principle that breaks the (mu_0, delta_0) degeneracy. This requires going beyond the FRW background retrodictability functional.

**Priority 3**: Test the prediction observationally. The running mu(k) = k predicts a specific scale dependence of the dark matter phenomenology. If the running is confirmed, the boundary condition mu_0 = H_0/c becomes an observational fact (from fitting the running at the Hubble scale), even if it's not derived theoretically.

### 9.3 The Philosophical Point

Sheng-Kai's philosophy -- "providing a new argument that connects known things, making things 'right,' is more meaningful than finding new equations" -- is directly relevant here. The insight is NOT that mu_0 = H_0/c is derivable (it isn't, currently). The insight is that this single parameter connects:

```
de Sitter temperature --> Khronon coupling --> dark matter density
                      --> MOND acceleration --> rotation curves
                      --> CMB compatibility (via running)
                      --> cosmic coincidence (via extremal retrodictability)
```

This unification of seemingly disparate phenomena through a single parameter is the physical content. Whether that parameter is "derived" or "assumed" is important for mathematical completeness but does not diminish the physical content.

The precedent is Newton's gravitational constant G. Newton did not derive G from first principles -- he measured it and showed that the SAME G explains falling apples, planetary orbits, and tidal forces. The value of G became derivable only in the context of string theory (where it relates to the string coupling and compactification volume), 300 years later. The tau framework may be in a similar position: mu_0 = H_0/c unifies many phenomena, and its derivation may require a deeper theory not yet available.

---

## 10. Summary of All Approaches

### Approach 1: Petz Bound at Cosmological Horizon
- **Method**: Sigma ~ O(1) inside the de Sitter horizon
- **Result**: mu ~ H/c (order of magnitude)
- **Classification**: HEURISTIC
- **Strength**: Gives the right scale; physically intuitive
- **Weakness**: O(1) Sigma condition not derived from any principle
- **Minimal assumptions**: Petz bound + "transition regime" criterion

### Approach 2: Running Consistency (IR Boundary Condition)
- **Method**: mu(k) = k evaluated at k = H_0/c
- **Result**: mu_0 = H_0/c (exact)
- **Classification**: TAUTOLOGICAL (as derivation of mu_0)
- **Strength**: Exact; connects to well-motivated running
- **Weakness**: Uses mu_0 = H_0/c as input (boundary condition)
- **Minimal assumptions**: eta_mu = 1 + boundary condition

### Approach 3: Extremal Retrodictability (Joint Optimization)
- **Method**: Maximize R[mu_0, delta_0] simultaneously
- **Result**: Cannot fix mu_0 (degeneracy with delta_0)
- **Classification**: NEGATIVE
- **Strength**: Clearly identifies the obstruction (Omega_K degeneracy)
- **Weakness**: The functional depends only on Omega_K, not on mu_0 and delta_0 separately
- **Minimal assumptions**: F = 1/(1+delta), Sigma = delta(2+delta), proper-time weighting

### Approach 4: Scale Invariance of On-Shell Khronon Action
- **Method**: Require S_K/S_EH to be H-independent
- **Result**: Tracking behavior with mu = H/c (suggestive but not unique)
- **Classification**: SUGGESTIVE
- **Strength**: Identifies nice property of mu = H/c (tracking)
- **Weakness**: Tracking holds for any normalization with eta_mu = 1
- **Minimal assumptions**: On-shell Khronon action + scale invariance

### Approach 5: DPI Maximal Rate
- **Method**: DPI constrains eta_mu >= 1; extensivity gives eta_mu <= 1; together: eta_mu = 1
- **Result**: Fixes the running but NOT the normalization
- **Classification**: PARTIAL
- **Strength**: Rigorous constraint on the running
- **Weakness**: Does not fix mu_0 (needs modular flow for normalization)
- **Minimal assumptions**: CPTP + de Sitter extensivity

---

## 11. The Bottom Line

**mu_0 = H_0/c is the most well-motivated assumption in the tau framework, but it IS an assumption.**

The five approaches analyzed here converge on a consistent picture but fail to close the gap between "well-motivated" and "derived." The deepest available argument chain is:

```
[STANDARD] Gibbons-Hawking 1977: de Sitter has temperature T_dS = hbar H/(2pi k_B c)
[PROVED]   Dimensional analysis: mu_0 = H_0/c is unique without Planck physics
[DERIVED]  DPI + extensivity: eta_mu = 1, so mu runs as mu(k) = A*k
[SPECULATIVE] Modular flow: Khronon = modular flow direction, so A = 1
[ASSUMED]  Combining: mu_0 = H_0/c
```

The weakest link is the modular flow identification (going from [DERIVED] to [SPECULATIVE]). If this could be proved, the entire chain would be on solid ground.

Until then, mu_0 = H_0/c should be treated as the framework's **foundational postulate** -- analogous to the speed of light c in special relativity or Newton's constant G in general relativity: a parameter whose value is set by nature and whose derivation awaits a deeper theory.

---

## References

### Dimensional Analysis and Naturality
1. Blanchet, L. & Skordis, C. (2024). JCAP 11, 040. arXiv:2404.06584 [Khronon action]
2. Milgrom, M. (1983). ApJ 270, 365 [a_0 ~ cH_0 coincidence]
3. Verlinde, E. (2016). SciPost Phys. 2, 016. arXiv:1611.02269 [Volume-law entropy -> DM]

### de Sitter Thermodynamics
4. Gibbons, G.W. & Hawking, S.W. (1977). PRD 15, 2738
5. Chandrasekaran, Longo, Penington, Witten (2023). JHEP 2023, 82

### DPI, Extensivity, and Running
6. Jacobson, T. (2016). PRL 116, 201101 [Entanglement equilibrium]
7. Casini, H., Huerta, M. & Myers, R.C. (2011). JHEP 1105, 036 [Modular Hamiltonian]
8. Kumar, S. (2025). arXiv:2509.05246 [Marginal IR running of gravity]

### Petz Recovery
9. Petz, D. (1988). QJM 39(1), 97-108
10. Junge, Renner, Sutter, Wilde, Winter (2018). Ann. Henri Poincare 19, 2955

### tau Framework
11. Huang, S.-K. (2026). Paper 1: Petz recovery unification
12. Huang, S.-K. (2026). Paper 3: Weak-field phenomenology
13. Huang, S.-K. (2026). Paper 4: Temporal asymmetry as organizing principle

### Previous Research Notes (this repository)
14. mu_synthesis.md (2026-03-12) [Master reference for mu story]
15. mu_running_proof.md (2026-03-12) [Five approaches to eta_mu = 1]
16. mu_route3_petz_optimality.md (2026-03-12) [Can Petz fix mu?]
17. omega_DM_derivation_2026_03_16.md (2026-03-16) [Retrodictability extremal principle]
18. mathematical_audit_2026_03_16.md (2026-03-16) [Full mathematical audit]

---

*Last updated: 2026-03-17*
*This document represents the most thorough analysis to date of whether mu_0 = H_0/c can be derived from first principles. The answer is: not yet, but the convergent evidence from five independent approaches makes it the best-motivated assumption in the framework.*
