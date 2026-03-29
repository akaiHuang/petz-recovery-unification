# The Conservative g00-Only Strategy for Paper 2

## Date: 2026-03-11
## Status: Comprehensive research analysis
## Purpose: Rigorously investigate which predictions depend only on g00, which need g_rr, what g_rr options are allowed, and how to present this in Paper 2

---

## Executive Summary

Paper 2 derives g00 = -exp(-r_s/r) from quantum information theory (Petz recovery fidelity). The spatial component g_rr is NOT derived from QI -- it is obtained via the isometry assumption g00 * g_rr = -1 in isotropic coordinates. This document rigorously classifies every prediction as "robust" (g00-only) or "conditional" (requires g_rr), analyzes allowed g_rr options with explicit calculations, and argues that the conservative strategy makes Paper 2 STRONGER than the Yilmaz theory it superficially resembles.

**Critical correction to the task premise**: The claim that the "full exponential metric predicts only half the observed light bending and only 1/3 the perihelion precession" is **incorrect**. This confusion arises from conflating two distinct things:
1. The **Yilmaz field equations** (a different theory of gravity) -- which DO have problems with higher-order effects
2. The **exponential metric** (a specific geometry) -- which has **identical** PPN parameters (gamma = beta = 1) to Schwarzschild and gives **exactly the same** first-order light bending and perihelion precession

Mychelkin et al. (2024, arXiv:2403.11610) proved explicitly that the exponential metric gives identical 1PN light deflection and perihelion precession. Differences appear only at 2PN order (~10^{-15} for solar system tests). The exponential metric is NOT ruled out by classical tests.

---

## 1. Predictions That Depend ONLY on g00

### 1.1 No-Horizon Theorem

**Statement**: If Sigma_grav < infinity, then tau < 1, and event horizons are forbidden.

**Proof that g_rr does not enter**:

The temporal asymmetry parameter is defined as:
```
tau(r) = 1 - sqrt(-g_00(r))
```

For g_00 = -exp(-r_s/r):
```
tau(r) = 1 - exp(-r_s/(2r))
```

Since exp(-r_s/(2r)) > 0 for all r > 0, we have tau < 1 everywhere.

**Mathematical proof**: The no-horizon theorem requires only:
1. The Petz bound F >= exp(-Sigma/2) (Paper 1, model-independent)
2. Sigma_grav < infinity for any finite-mass object at finite distance
3. Therefore F > 0, hence tau < 1

None of these three steps reference g_rr, g_{theta theta}, or any spatial metric component. The spatial geometry is completely irrelevant.

**Status**: ROBUST -- holds for ANY g_rr.

### 1.2 Gravitational Redshift

**Statement**: The gravitational redshift between radii r_1 and r_2 is:
```
z = sqrt(g_00(r_1) / g_00(r_2)) - 1
```

**Proof that g_rr does not enter**:

The standard derivation proceeds via the geodesic equation for photons. For a static metric ds^2 = g_00 dt^2 + g_rr dr^2 + ..., the Killing vector xi^a = (d/dt)^a gives a conserved energy:
```
E = -g_00 (dt/dlambda)
```

The frequency measured by a static observer at radius r is:
```
omega(r) = E / sqrt(-g_00(r))
```

Therefore:
```
omega(r_1)/omega(r_2) = sqrt(-g_00(r_2)) / sqrt(-g_00(r_1))
```

This ratio depends exclusively on g_00 at the two points. The spatial metric g_rr does not appear because the derivation uses only the timelike Killing vector and the norm of the 4-velocity of static observers, both of which depend only on g_00.

**Status**: ROBUST -- independent of g_rr. This is a fundamental result in GR, not specific to any particular metric.

### 1.3 Gravitational Time Dilation

**Statement**: The proper time rate for a static observer at radius r relative to infinity is:
```
dtau/dt = sqrt(-g_00(r))
```

**Proof**: This follows directly from the metric for a static observer (dr = dtheta = dphi = 0):
```
ds^2 = g_00 c^2 dt^2
dtau = sqrt(-g_00) dt
```

The spatial metric components are absent because the observer is static.

**Status**: ROBUST -- independent of g_rr.

### 1.4 Coordinate Speed of Light (Refractive Index)

**Statement**: c_eff/c = exp(-Sigma_grav/2) = F_bound.

**Subtlety**: The coordinate speed of light DOES depend on g_rr in general. For a radial null ray (ds^2 = 0):
```
g_00 c^2 dt^2 + g_rr dr^2 = 0
c_eff = dr/dt = c * sqrt(-g_00/g_rr)
```

For the isotropic exponential metric with g_rr = exp(+r_s/r):
```
c_eff = c * sqrt(exp(-r_s/r) / exp(+r_s/r)) = c * exp(-r_s/r)
```

But in Schwarzschild-areal coordinates where g_rr = 1/(1-r_s/r):
```
c_eff = c * (1 - r_s/r)
```

**However**: The ISOTROPIC coordinate speed c_eff = c * sqrt(-g_00 / g_rr) involves g_rr. What Paper 2 actually identifies is:
```
sqrt(-g_00) = exp(-Sigma_grav/2) = F_bound
```

This is g_00-only. The refractive index interpretation n = 1/sqrt(-g_00) is the Dicke scalar refractive index, which depends only on g_00.

**Status**: The identification sqrt(-g_00) = F_bound is ROBUST. The full coordinate speed c_eff = c*sqrt(-g_00/g_rr) is CONDITIONAL on g_rr. Paper 2 should be careful to distinguish these.

### 1.5 Echo Existence (Yes/No)

**Statement**: The existence of echoes (not their time delay) follows from the no-horizon theorem.

**Proof**: If tau < 1 everywhere, then no event horizon forms, and the gravitational potential well is finite. Any compact object without an event horizon acts as a "partial reflector" for gravitational waves -- some radiation bounces back from the potential barrier, producing echoes.

This is a topological/qualitative argument: echoes exist if and only if there is no event horizon. Since tau < 1 follows from g_00 alone, echo existence is g_00-only.

**Caveat**: The echo AMPLITUDE depends on the shape of the effective potential, which involves g_rr. A very small potential barrier could produce undetectable echoes. The existence is robust; the detectability is conditional.

**Status**: Existence is ROBUST. Amplitude and time delay are CONDITIONAL.

### 1.6 Kodama Vector and Apparent Horizon Definition

**Statement**: tau_K = 1 defines the apparent horizon surface.

**Analysis**: The Kodama vector for a general spherically symmetric metric
```
ds^2 = -e^{2alpha(r,t)} dt^2 + e^{2beta(r,t)} dr^2 + R(r,t)^2 dOmega^2
```
is defined as K^a = epsilon^{ab} nabla_b R, where epsilon^{ab} is the Levi-Civita tensor on the (t,r) submanifold.

In the static case, |K| = sqrt(-g_00), independent of g_rr. This is because the Kodama vector reduces to the timelike Killing vector, whose norm is sqrt(-g_00).

**However**, in the dynamic case (e.g., Vaidya), the Kodama vector norm is:
```
|K|^2 = -g_00 - 2 g_{0r} (dR/dt)/(dR/dr) - g_{rr} (dR/dt)^2/(dR/dr)^2
```

For the standard Vaidya metric in Eddington-Finkelstein coordinates (ds^2 = -(1-2M(v)/r) dv^2 + 2 dv dr + r^2 dOmega^2), the Kodama vector norm simplifies to |K| = sqrt(1 - 2M(v)/r), which involves only the effective g_{vv} component.

**Key result**: The condition tau_K = 1 (apparent horizon) always corresponds to the vanishing of the outgoing null expansion theta_+ = 0. In the Schwarzschild case, this is r = r_s, where g_00 = 0 (for Schwarzschild) or g_00 -> 0 asymptotically (for exponential). The locus of the apparent horizon depends on g_00, not g_rr.

**Status**: In the static limit, ROBUST (g_00-only). In the dynamic case, the Kodama vector involves the full 2D metric on the (t,r) submanifold, but the apparent horizon condition theta_+ = 0 is still determined by the Misner-Sharp mass 2M/R = 1, which in the static limit reduces to a condition on g_00.

### 1.7 tau Locality Argument

**Statement**: tau(r) = 1 - sqrt(-g_00(r)) depends only on local data at point r.

**Proof**: By definition, tau involves only g_00(r) -- a single component of the metric at a single point. No derivatives, no integrals, no spatial metric components.

This is the core argument of Paper 2: information loss is local, determined by g_00, not by global causal structure. The event horizon requires knowledge of the entire future spacetime; tau requires only one number at one point.

**Status**: ROBUST -- this is a definitional statement.

### 1.8 Gravitational Decoherence Rate (Pikovski)

**Statement**: The decoherence rate is controlled by Delta_Phi/c^2 ~ Sigma_grav/2.

**Proof**: The Pikovski gravitational decoherence mechanism depends on the proper time difference between two heights:
```
Delta_tau/tau = sqrt(-g_00(r + Delta_r)) - sqrt(-g_00(r)) ~ (1/2) * g_00'(r)/g_00(r) * Delta_r
```

This involves only g_00 and its radial derivative. The spatial metric g_rr enters only through the physical separation between the paths, but in the weak field (Delta_Phi/c^2 << 1), the dominant effect is purely temporal.

**Status**: ROBUST to leading order. Sub-leading corrections may involve g_rr through the physical distance.

---

## 2. Predictions That NEED g_rr

### 2.1 Light Bending Angle

**Where g_rr enters**:

The null geodesic equation for light passing a mass M involves both metric components. For a general static spherically symmetric metric:
```
ds^2 = -A(r) c^2 dt^2 + B(r) dr^2 + C(r) r^2 dOmega^2
```

The orbit equation for light is:
```
(du/dphi)^2 = C(u)/B(u) * [C(u)/(b^2 A(u)) - 1]
```
where u = 1/r and b is the impact parameter.

In the PPN formalism, the metric is written as:
```
g_00 = -(1 - 2U/c^2 + 2 beta U^2/c^4 + ...)
g_ij = delta_ij (1 + 2 gamma U/c^2 + ...)
```
where U = GM/r is the Newtonian potential, and gamma, beta are PPN parameters.

**The light deflection angle at 1PN is**:
```
Delta_phi = (1 + gamma)/2 * 4GM/(bc^2)
```

For GR: gamma = 1, giving Delta_phi = 4GM/(bc^2) (the standard Einstein result).

**The crucial factor (1+gamma)/2** shows explicitly that light bending has two equal contributions:
- Factor of 1/2 from g_00 (time curvature, equivalent principle, Einstein 1911)
- Factor of gamma/2 from g_ij (spatial curvature)

For the exponential metric in isotropic coordinates:
```
g_00 = -exp(-r_s/r) = -1 + r_s/r - r_s^2/(2r^2) + ...
g_ij = delta_ij exp(+r_s/r) = delta_ij (1 + r_s/r + r_s^2/(2r^2) + ...)
```

Reading off PPN parameters: gamma = 1 (the coefficient of r_s/r in g_ij matches), beta = 1 (the coefficient of (r_s/r)^2 in g_00 matches).

**Therefore**: The exponential metric gives EXACTLY the same 1PN light bending as Schwarzschild. The "half bending" problem does NOT apply. Differences appear only at 2PN:

From Mychelkin et al. (2024):
- Schwarzschild impact parameter: J_S = 1 + 2M/r_0 + 7M^2/(4r_0^2) + O(M^3)
- Exponential impact parameter: J_P = 1 + 2M/r_0 + 2M^2/r_0^2 + O(M^3)

The 2PN coefficients differ (7/4 vs 2), but this difference is ~ (r_s/r)^2 ~ 10^{-12} in the solar system -- completely unmeasurable.

**Status**: CONDITIONAL on g_rr. But for any g_rr with gamma = 1, the 1PN light bending is identical to GR.

### 2.2 Perihelion Precession

**Where g_rr enters**:

The precession per orbit is:
```
Delta_phi = (2 + 2gamma - beta)/3 * 6pi GM/(a(1-e^2)c^2)
```

For GR (gamma = beta = 1): Delta_phi = 6pi GM/(a(1-e^2)c^2) (standard Einstein result).

Both the exponential metric and Schwarzschild have gamma = beta = 1, so they give IDENTICAL 1PN precession:
```
Delta_phi = 6pi M/L    [L = a(1-e^2)]
```

From Mychelkin et al. (2024), the 2PN difference:
```
Delta_phi_Schw = 6pi M/L + (25pi/2) M^2/L^2 - (9pi/2) M^2/N^2
Delta_phi_exp  = 6pi M/L + (41pi/3) M^2/L^2 - 4pi M^2/N^2
```

The difference: Delta_phi_exp - Delta_phi_Schw = (1/6) pi M^2 (7/L^2 + 3/N^2)

For Mercury: this is ~ 10^{-14} arcsec/century (unmeasurable; current precision ~ 0.002"/century).

**Status**: CONDITIONAL on g_rr. But identical to GR at 1PN for gamma = beta = 1.

### 2.3 Black Hole Shadow Size/Shape

**Where g_rr enters**:

The photon sphere radius and critical impact parameter depend on both A(r) = -g_00 and B(r) = g_rr through the effective potential:
```
V_eff(r) = A(r) * L^2 / (B(r) * r^2)    [for non-isotropic coords]
```
or equivalently:
```
V_eff(r) = A(r) / r^2    [in isotropic coords with C(r) = B(r)]
```

The photon sphere is at dV_eff/dr = 0. For the exponential metric (isotropic):
```
V_eff = exp(-r_s/r) / r^2
dV_eff/dr = exp(-r_s/r) * (r_s/r^2 - 2/r) / r^2 = 0
=> r_ps = r_s/2 = M    [isotropic coordinates]
```

The critical impact parameter:
```
b_crit = r_ps / sqrt(V_eff(r_ps)) = r_ps * exp(r_s/(2r_ps)) = M * e = 2.718M
```

In areal coordinates, this corresponds to a shadow angular diameter ~ 5.437 M, compared to Schwarzschild's 3sqrt(3) M ~ 5.196 M, a 4.6% difference.

**Important**: In isotropic coordinates, the effective potential for light involves BOTH g_00 and g_rr implicitly through the coordinate system. Even though V_eff = A(r)/r^2 appears to involve only g_00, the relationship between isotropic r and areal R depends on g_rr:
```
R = r * sqrt(g_rr(r) / g_00(r))^{1/2}    [schematic]
```

**Status**: CONDITIONAL on g_rr. Different g_rr choices give different shadow sizes.

### 2.4 Echo TIME DELAY

**Where g_rr enters**:

The echo time delay is determined by the tortoise coordinate:
```
r* = integral sqrt(g_rr / (-g_00)) dr
```

For the exponential metric:
```
r* = integral exp(r_s/r) dr
```

This integral is BOUNDED (converges as r -> 0), giving a finite delay:
```
Delta_t_echo ~ 4.17 r_s/c
```

For Schwarzschild:
```
r* = integral dr/(1 - r_s/r) = r + r_s ln|r/r_s - 1|
```

This DIVERGES as r -> r_s (horizon), giving infinite delay (no echoes -- absorbed by horizon).

**The tortoise coordinate explicitly involves g_rr**. The echo delay time is one of the most g_rr-sensitive predictions.

Different g_rr choices would give:
- g_rr = exp(+r_s/r): r* = integral exp(r_s/r) dr ~ 4.17 r_s/c delay (Paper 2 prediction)
- g_rr = 1/(1-r_s/r): r* diverges at r = r_s (no finite echo delay even with exp g_00)
- g_rr = 1 + r_s/r: r* = r + r_s ln(r) -- bounded, gives different delay

**Status**: CONDITIONAL on g_rr. The echo delay is highly sensitive to the spatial metric.

### 2.5 Gravitational Wave Propagation Speed

**Where g_rr enters**:

The speed of gravitational waves in a curved background is determined by the null condition ds^2 = 0 for the perturbation wavefront:
```
c_gw = dr/dt = c * sqrt(-g_00 / g_rr)
```

This is identical to the coordinate speed of light and explicitly involves both g_00 and g_rr.

**Status**: CONDITIONAL on g_rr.

### 2.6 ISCO Radius

**Where g_rr enters**:

The ISCO is determined by the effective potential for timelike geodesics. In isotropic coordinates:
```
V_eff = sqrt(-g_00) * sqrt(1 + L^2/(g_rr * r^2))
```

The ISCO conditions require V_eff' = 0 and V_eff'' = 0 simultaneously. These conditions involve both g_00 and g_rr.

For the exponential metric:
- R_ISCO ~ 3.17 r_s (areal radius), compared to Schwarzschild's 3 r_s (+5.6%)

From Ernazarov (2025): r_ISCO = 2M in isotropic coordinates (for pure exponential with l=0).

**Status**: CONDITIONAL on g_rr.

### 2.7 Photon Sphere Radius

**Where g_rr enters**: Same analysis as shadow (Section 2.3). The photon sphere is defined by the effective potential for null geodesics, which involves the full metric.

For the exponential metric: r_ps = r_s = 2M (isotropic) = M (Ernazarov's convention)
For Schwarzschild: r_ps = 3M (areal) ~ 1.87 M (isotropic)

**Status**: CONDITIONAL on g_rr.

### 2.8 Quasinormal Mode Frequencies

**Where g_rr enters**: QNMs are solutions to the wave equation on the spacetime background. The Regge-Wheeler/Zerilli equation involves the tortoise coordinate and effective potential, both of which depend on g_rr.

The eikonal QNM frequency is:
```
omega_QNM = Omega_ps * (l + 1/2) - i * lambda_Lyapunov * (n + 1/2)
```

where Omega_ps is the orbital frequency at the photon sphere. This depends on both g_00 and g_rr at r_ps.

**Status**: CONDITIONAL on g_rr.

---

## 3. Allowed g_rr Options and Their Predictions

### 3.1 Option A: g_rr = exp(+r_s/r) (Pure Exponential / Papapetrou)

**Combined metric** (isotropic coordinates):
```
ds^2 = -exp(-r_s/r) c^2 dt^2 + exp(+r_s/r) (dr^2 + r^2 dOmega^2)
```

**Properties**:
- Isometry condition: g_00 * g_rr = -exp(-r_s/r) * exp(+r_s/r) = -1 (satisfied exactly)
- Traversable wormhole (Boonserm et al. 2018)
- Phantom scalar field exact solution (Makukov & Mychelkin 2020)
- PPN parameters: gamma = 1, beta = 1

**1PN predictions** (identical to Schwarzschild):
- Light deflection: 4GM/(bc^2) = 1.75" for Sun
- Perihelion precession: 6piGM/(a(1-e^2)c^2) = 42.98"/century for Mercury
- Shapiro delay: identical to Schwarzschild

**Strong-field predictions** (differ from Schwarzschild):
- Shadow: b_crit = 2eM ~ 5.437M (vs 5.196M, +4.6%)
- ISCO: R_ISCO ~ 3.17 r_s (vs 3 r_s, +5.6%)
- QNM: f_QNM ~ 0.956 f_Schw (eikonal, -4.4%)
- Echo delay: ~4.17 r_s/c (~2.5 ms for 60 M_sun)
- No event horizon

**Compatibility with tests**: Passes all solar system tests. Compatible with EHT (within error bars). Marginal tension with GW250114 pSEOBNR analysis but consistent with pyRing analysis.

### 3.2 Option B: g_rr = (1 - r_s/r)^{-1} (Schwarzschild spatial part)

**Combined metric** (Schwarzschild-like coordinates):
```
ds^2 = -exp(-r_s/r) c^2 dt^2 + dr^2/(1 - r_s/r) + r^2 dOmega^2
```

**Properties**:
- NOT an isometry: g_00 * g_rr = exp(-r_s/r)/(1 - r_s/r) != -1
- g_rr diverges at r = r_s, creating a coordinate singularity
- But g_00 = -exp(-r_s/r) != 0 at r = r_s, so tau < 1 there
- This creates a peculiar geometry: a "frozen star" surface without an actual horizon

**Coordinate speed of light**:
```
c_eff = c * sqrt(exp(-r_s/r) * (1 - r_s/r))
```
This vanishes at r = r_s (coordinate singularity), but the redshift z = exp(r_s/(2r)) - 1 remains finite.

**PPN parameters**: gamma = 1, beta = 1 (same weak-field behavior)

**Light deflection and precession**: Identical to GR at 1PN. At 2PN, differs from both pure exponential and Schwarzschild.

**Shadow calculation**: Would need numerical computation. The photon sphere location changes because the effective potential has a different form:
```
V_eff = exp(-r_s/r) * (1 - r_s/r) * L^2/r^2
```
This has different critical points than either pure Schwarzschild or pure exponential.

**Assessment**: Physically odd -- mixes two different metric paradigms. Not well-motivated by any theory. The no-horizon theorem still holds (tau < 1), but the geometry has unusual features.

### 3.3 Option C: g_rr = 1/|g_00| = exp(+r_s/r) (Isotropic Condition)

This IS Option A. The isotropic condition g_00 * g_rr = -1 in isotropic coordinates uniquely gives g_rr = exp(+r_s/r) when g_00 = -exp(-r_s/r). This is the Dicke scalar refractive index condition n^2 = 1/(-g_00), which implies g_rr = n^2 = 1/(-g_00) = exp(+r_s/r).

### 3.4 Option D: g_rr = 1 + r_s/r + alpha_2 (r_s/r)^2 + ... (Generic Weak-Field Expansion)

**Combined metric**:
```
ds^2 = -(1 - r_s/r + r_s^2/(2r^2) - ...) c^2 dt^2 + (1 + r_s/r + alpha_2 (r_s/r)^2 + ...) dr^2 + ...
```

**PPN parameters**:
- The coefficient of r_s/r in g_rr is 1, matching gamma = 1 (required by Cassini: |gamma-1| < 2.3 x 10^{-5})
- The coefficient alpha_2 of (r_s/r)^2 in g_rr is a free parameter at 2PN

For the exponential metric: alpha_2 = 1/2
For Schwarzschild (in isotropic coords): alpha_2 = 3/4

**Light bending at 2PN**:
```
Delta_phi = 4GM/(bc^2) * [1 + (15pi/16 + f(alpha_2)) * (r_s/b)^2 + ...]
```
where f(alpha_2) depends on the specific value of alpha_2. The correction is ~ 10^{-12} for solar grazing light.

**Assessment**: Any alpha_2 gives identical 1PN predictions. Solar system tests cannot distinguish different alpha_2 values. The choice of alpha_2 is constrained only by:
1. Strong-field observations (EHT, GW)
2. Theoretical consistency (energy conditions, field equations)

### 3.5 Option E: Schwarzschild in Isotropic Coordinates

For completeness, Schwarzschild in isotropic coordinates:
```
g_00^Schw = -((1 - r_s/(4r))/(1 + r_s/(4r)))^2
g_rr^Schw = (1 + r_s/(4r))^4
```

Expansion:
```
g_00^Schw = -1 + r_s/r - r_s^2/(2r^2) + r_s^3/(4r^3) + ...
g_rr^Schw = 1 + r_s/r + 3r_s^2/(4r^2) + ...
```

Compare with exponential:
```
g_00^exp = -1 + r_s/r - r_s^2/(2r^2) + r_s^3/(6r^3) + ...
g_rr^exp = 1 + r_s/r + r_s^2/(2r^2) + ...
```

**First differences**:
- In g_00: at O((r_s/r)^3): coefficient -1/4 (Schw) vs -1/6 (exp)
- In g_rr: at O((r_s/r)^2): coefficient 3/4 (Schw) vs 1/2 (exp)

### 3.6 Summary Table: g_rr Options

| Option | g_rr | gamma | beta | Shadow b_crit | ISCO | Echo delay | Theory basis |
|--------|------|-------|------|---------------|------|------------|-------------|
| A (Papapetrou) | exp(+r_s/r) | 1 | 1 | 5.437M (+4.6%) | 3.17 r_s (+5.6%) | 4.17 r_s/c | Phantom scalar |
| B (hybrid) | 1/(1-r_s/r) | 1 | 1 | Needs calc | Needs calc | Divergent | None (ad hoc) |
| D (generic) | 1 + r_s/r + alpha_2 (r_s/r)^2 | 1 | 1 | Depends on alpha_2 | Depends on alpha_2 | Depends on alpha_2 | Free |
| E (Schwarzschild) | (1+r_s/(4r))^4 | 1 | 1 | 5.196M (GR) | 3 r_s (GR) | Infinite (horizon) | Einstein vacuum |

**All options have gamma = beta = 1 and pass all solar system tests.**

---

## 4. Why the g00-Only Strategy is STRONGER Than Yilmaz

### 4.1 The Yilmaz Theory Problem

The Yilmaz theory is a complete alternative theory of gravity with modified field equations:
```
G_ab + t_ab^{grav} = 8pi G T_ab
```
where t_ab^{grav} is Yilmaz's gravitational stress-energy tensor. This theory:
- NEEDS to specify the full metric (it predicts all components)
- Was criticized by Misner (1995, gr-qc/9504050) for issues with the stress-energy tensor
- Has difficulties with self-consistency in strong fields
- Is a COMPETITOR to GR with different field equations

### 4.2 Paper 2's Position

Paper 2 does NOT propose a new theory of gravity. It:
- Accepts GR (or any metric theory) as the framework
- Derives g_00 from quantum information (Petz recovery fidelity)
- Makes NO claim about g_rr
- States explicitly that g_rr requires additional physics

This is a fundamentally different epistemic position:
- **Yilmaz**: "My theory predicts the full metric, and it must pass all tests" -> FAILS in strong field
- **Paper 2**: "QI constrains g_00; g_rr is a separate question" -> CANNOT fail from g_rr tests

### 4.3 The Robustness Argument (Formalized)

**Theorem (Informal)**: Let M be a static spherically symmetric metric with g_00 = -exp(-r_s/r). The following predictions hold for ALL choices of g_rr:

1. tau(r) = 1 - exp(-r_s/(2r)) < 1 for all r > 0
2. Gravitational redshift: z = exp(r_s/(2r)) - 1
3. No event horizon (tau < 1 everywhere)
4. Gravitational decoherence rate ~ (Delta_Phi/c^2)^2 ~ (r_s/(2r))^2
5. The apparent horizon (if it exists) is determined by g_00 = 0, which never occurs

**Proof**: All five statements follow from tau(r) = 1 - sqrt(-g_00(r)) and the properties of g_00 = -exp(-r_s/r). No spatial metric component enters any of these derivations. QED.

**Corollary**: Paper 2's core results (no-horizon theorem, tau locality, refractive index identification) cannot be falsified by ANY measurement of spatial geometry (light bending, shadows, precession, etc.).

### 4.4 The Key Distinction

| Aspect | Yilmaz Theory | Paper 2 (tau framework) |
|--------|---------------|----------------------|
| Claims | Full metric from modified field equations | g_00 from QI; g_rr undetermined |
| Tests it must pass | ALL classical tests (bending, precession, shadow) | Only redshift, time dilation tests |
| Vulnerability | Fails if full metric disagrees with data | CANNOT fail from spatial tests |
| Strength | Ambitious but fragile | Conservative but robust |
| Nature | Alternative gravity theory | QI constraint on any metric theory |

---

## 5. Is There a Principled Reason to Leave g_rr Undetermined?

### 5.1 Why QI Constrains g_00 but Not g_rr

**Physical reason**: Quantum channels model TEMPORAL processes:
- State preparation at time t_1 -> state measurement at time t_2
- The channel N maps rho(t_1) to rho(t_2)
- Entropy production Sigma = D(rho||sigma) - D(N(rho)||N(sigma)) quantifies irreversibility of this temporal process
- The Petz recovery map attempts to "undo" the temporal evolution

All of these are inherently TEMPORAL operations. They probe the relationship between "now" and "later" -- which is encoded in g_00 (how clocks tick at different locations).

There is no analogous standard QI operation that directly probes spatial geometry:
- You cannot "send a state through space" without also evolving it in time
- Spatial entanglement (mutual information between regions) probes spatial structure, but the connection to g_rr is indirect

### 5.2 What QI Quantity Might Determine g_rr?

Several candidates exist, but none is as clean as the Petz fidelity -> g_00 connection:

**(a) Spatial mutual information**:
The mutual information I(A:B) between spatially separated regions A and B encodes the spatial metric through the Ryu-Takayanagi formula:
```
S_A = Area(gamma_A) / (4G)
```
The area of the minimal surface gamma_A depends on the full spatial metric. In principle, varying A probes different aspects of g_rr.

**Problem**: This works cleanly only in AdS/CFT (holographic systems). For asymptotically flat spacetimes, the connection is much less developed.

**(b) Entanglement structure (Cao-Carroll-Michalakis 2017)**:
The spatial metric can be reconstructed from the entanglement structure of the quantum state:
- Entanglement entropy across codimension-one surfaces defines areas
- A Radon transform converts area data into spatial geometry

**Problem**: Requires knowing the full quantum state, not just a channel property. This is a state-dependent reconstruction, not a channel-dependent constraint.

**(c) Bisognano-Wichmann spatial modular flow**:
The most promising candidate (identified in paper2_full_metric_from_tau.md). The Bisognano-Wichmann theorem generates boosts that mix t and x. The modular flow in the (t,r) plane constrains BOTH g_00 and g_rr simultaneously.

**Status**: OPEN -- this is the deepest approach but requires new formalism.

**(d) Jacobson/Dorau-Much program**:
These programs DO derive the full Einstein equation (all metric components) from information-theoretic premises. But their vacuum solution is Schwarzschild, not exponential. This suggests that if one takes the full Einstein equation seriously, g_rr is determined -- and it's NOT exp(+r_s/r).

**Tension**: The Jacobson/Dorau-Much result gives Schwarzschild g_rr from information theory. Paper 2's g_00 = exp(-r_s/r) is NOT the Schwarzschild g_00. This means either:
1. Paper 2's g_00 is an approximation (and the true metric is Schwarzschild), or
2. The exponential metric requires matter (phantom scalar) -- confirmed by the field equation analysis, or
3. The Jacobson program needs modification to account for quantum effects that produce the exponential correction

### 5.3 The Honest Assessment

The fact that QI constrains g_00 but not g_rr is not a bug -- it is a feature. It reflects a genuine physical distinction:

- **Temporal geometry** (g_00): directly probed by quantum channels, clocks, and information recovery
- **Spatial geometry** (g_rr): probed by light propagation, geodesics, and entanglement structure -- all of which involve both temporal AND spatial geometry simultaneously

The clean separation exists only for g_00. Any operational procedure that probes g_rr necessarily also involves g_00 (because probes must evolve in time). This asymmetry is fundamental, not a limitation of the framework.

---

## 6. Correcting the Task Premise: Yilmaz vs Exponential Metric

### 6.1 The "Half Bending" Myth

The task description states that the "FULL exponential metric predicts only half the observed light bending, only 1/3 the perihelion precession, and shadow size 48% smaller than Schwarzschild."

**This is WRONG.** The confusion comes from conflating:

1. **Yilmaz THEORY** (alternative field equations): Early versions of Yilmaz's theory had the linearized metric g_00 = -(1-2Phi/c^2), g_rr = 1 (NOT 1+2Phi/c^2). This gives gamma = 0, which gives half the bending: Delta_phi = (1+0)/2 * 4M/b = 2M/b. This was the version criticized.

2. **Exponential METRIC** (specific geometry): The metric ds^2 = -exp(-r_s/r)dt^2 + exp(+r_s/r)(dr^2 + r^2 dOmega^2) has BOTH components, with gamma = 1. This gives the FULL 4M/b bending.

The "half bending" applies to the linearized Yilmaz theory with gamma = 0, NOT to the exponential metric with gamma = 1.

### 6.2 What About the "48% Smaller Shadow"?

This also appears to be incorrect. The exponential metric gives a shadow that is 4.6% LARGER than Schwarzschild (b_crit = 5.437M vs 5.196M), not 48% smaller. The 48% figure does not correspond to any known calculation of the exponential metric shadow.

### 6.3 The Correct Comparison

| Observable | Schwarzschild | Exponential metric | Yilmaz (gamma=0) |
|-----------|--------------|-------------------|-------------------|
| Light bending (1PN) | 4GM/bc^2 | 4GM/bc^2 (SAME) | 2GM/bc^2 (HALF) |
| Perihelion precession (1PN) | 6piGM/(Lc^2) | 6piGM/(Lc^2) (SAME) | Different |
| Shadow | 5.196M | 5.437M (+4.6%) | Unknown |
| PPN gamma | 1 | 1 | 0 or context-dependent |

**The exponential metric passes all classical tests. The Yilmaz theory (with gamma=0) does not.**

Paper 2 uses the exponential METRIC, not the Yilmaz THEORY. This distinction is critical.

---

## 7. How Paper 2 Should Present This

### 7.1 Current Handling (Already Good)

Paper 2 already contains (lines 599-623):

> "The spatial metric is assumed, not derived. --- The information-theoretic framework constrains only the temporal component g_00 through sqrt(-g_00) = exp(-Sigma_grav/2). The spatial component g_rr = exp(+r_s/r) follows from the isometry condition g_00 * g_rr = -1 (isotropic coordinates), which is equivalent to Dicke's scalar refractive index. This is a physically motivated assumption --- not a derivation from quantum information theory."

And (lines 615-618):

> "Consequently, we classify our predictions as robust (depending only on g_00: no-horizon theorem, refractive index, gravitational decoherence rate) or conditional (depending on the full metric: shadow size, ISCO radius, QNM frequencies, echo delay)."

### 7.2 Recommended Additions

**(a) Explicit acknowledgment of prior exponential metric literature:**

> "The exponential metric was first proposed by Papapetrou (1954) and independently by Yilmaz (1958) and Rosen (1973). Makukov and Mychelkin (2020) showed that three independent scalar-field families within GR all reduce to this metric. Boonserm et al. (2018) classified the geometry as a traversable wormhole. Our contribution is not the metric itself, but its information-theoretic derivation: the specific prediction g_00 = -exp(-r_s/r) follows from the Petz recovery fidelity bound, not from modified field equations. The spatial component g_rr is not constrained by this derivation."

**(b) Clear distinction from Yilmaz theory:**

> "We emphasize that our framework is distinct from the Yilmaz theory of gravitation (Yilmaz 1958, 1971). Yilmaz proposed modified gravitational field equations that necessarily predict both metric components; criticisms of that theory (Misner 1995) do not apply here. Our framework makes no modification to Einstein's field equations and no prediction about g_rr. The exponential g_00 arises as an information-theoretic constraint within any metric theory of gravity."

**(c) Formalize the robust/conditional classification:**

> "To make the logical structure precise:
>
> Robust predictions (g_00-only, hold for ANY g_rr):
> - No-horizon theorem: tau < 1 for Sigma < infinity
> - Gravitational redshift: z = exp(r_s/(2r)) - 1
> - tau locality: information loss is local
> - Apparent horizon tracking via Kodama vector
> - Gravitational decoherence bound
>
> Conditional predictions (require g_rr assumption):
> - Shadow size (+4.6% vs Schwarzschild)
> - QNM frequencies (-4.4% to -7.8%)
> - Echo delay (~4.17 r_s/c)
> - ISCO radius (+5.6%)
>
> The conditional predictions assume the isometry condition g_rr = 1/(-g_00), equivalent to Dicke's scalar refractive index. This is the simplest choice consistent with the Newtonian limit and isotropy, but not the only one."

**(d) What would determine g_rr:**

> "The spatial metric component requires physics beyond the Petz recovery framework. Possible routes include: (i) the Jacobson/Dorau-Much program, which derives all components but yields Schwarzschild, not exponential; (ii) spatial modular flow (Bisognano-Wichmann), which constrains both temporal and spatial geometry simultaneously; (iii) holographic entanglement structure (Ryu-Takayanagi/Cao-Carroll-Michalakis), which reconstructs spatial geometry from entanglement data. Determining g_rr from quantum information remains an important open problem."

---

## 8. Summary: The Conservative Strategy Strengthens Paper 2

### 8.1 Three Levels of Confidence

| Level | Content | Confidence | Falsifiable by |
|-------|---------|-----------|---------------|
| **Layer 0** | Petz bound F >= exp(-Sigma/2) | Mathematical theorem | Nothing (it's proven) |
| **Layer 1** | tau(r) = 1 - sqrt(-g_00) is local, no horizon | Model-independent (g_00-only) | Finding tau = 1 with finite Sigma |
| **Layer 2** | g_00 = -exp(-r_s/r) (saturation) | Motivated by 3 routes, approximate | Strong-field g_00 measurements |
| **Layer 3** | g_rr = exp(+r_s/r) (isometry) | Physically motivated assumption | Shadow, QNM, echo measurements |

The conservative strategy explicitly separates Layer 1-2 (robust, g_00-only) from Layer 3 (conditional, assumes g_rr). This is intellectually honest and strategically sound:

- If Layer 3 is falsified (wrong g_rr), Layers 0-2 survive intact
- If Layer 2 is falsified (g_00 is not exactly exponential), Layer 0-1 still hold
- Layer 0-1 is essentially unfalsifiable given the mathematical rigor of the Petz bound

### 8.2 Comparison with Yilmaz

The Yilmaz theory was an all-or-nothing proposition: it predicted the full metric from modified field equations. When the full metric failed tests (in certain formulations), the entire theory was rejected.

Paper 2's conservative strategy avoids this trap. By deriving only g_00 and being explicit about what requires additional assumptions, the core results are immune to any future measurement of g_rr-dependent observables.

### 8.3 The Deepest Point

The reason g_00 can be derived from QI but g_rr cannot is not a failure of the framework -- it reflects a deep physical truth: time and space are not informationally equivalent.

- **Temporal information loss** (g_00): directly measurable by quantum channels
- **Spatial information structure** (g_rr): accessible only through operations that necessarily involve both time and space

This asymmetry between time and space at the information-theoretic level may be the most profound implication of the tau framework.

---

## References

### Exponential Metric / Yilmaz / Papapetrou
1. Papapetrou, A. (1954). Z. Phys. 139, 518.
2. Yilmaz, H. (1958). Phys. Rev. 111, 1417.
3. Rosen, N. (1973). Gen. Relativ. Gravit. 4, 435.
4. Makukov, M.A. & Mychelkin, E.G. (2020). Found. Phys. 50, 1346. arXiv:2009.08655.
5. Boonserm, P. et al. (2018). Phys. Rev. D 98, 084048. arXiv:1805.03781.
6. Misner, C.W. (1995). arXiv:gr-qc/9504050.

### Exponential Metric Observational Tests
7. Ernazarov, K.K. (2025). arXiv:2511.13471.
8. Mychelkin, E.G. et al. (2024). Gen. Relativ. Gravit. arXiv:2403.11610.
9. Nath, P.P. & Sarma, D. (2024). Eur. Phys. J. C 84, 1063. arXiv:2401.03738.

### PPN Formalism and Classical Tests
10. Will, C.M. (2014). Living Rev. Relativ. 17, 4.
11. Bertotti, B., Iess, L. & Tortora, P. (2003). Nature 425, 374.

### QI and Gravity
12. Jacobson, T. (1995). Phys. Rev. Lett. 75, 1260.
13. Dorau, J. & Much, A. (2025). Phys. Rev. Lett. arXiv:2510.24491.
14. Dicke, R.H. (1957). Rev. Mod. Phys. 29, 363.

### Holographic Spatial Reconstruction
15. Ryu, S. & Takayanagi, T. (2006). Phys. Rev. Lett. 96, 181602.
16. Cao, C., Carroll, S.M. & Michalakis, S. (2017). Phys. Rev. D 95, 024031.

### Kodama Vector
17. Kodama, H. (1980). Prog. Theor. Phys. 63, 1217.
18. Hayward, S.A. (1994). Phys. Rev. D 49, 6467.
19. Abreu, G. & Visser, M. (2010). Phys. Rev. D 82, 044027.
