# Paper 4 Critical Calculation: Can Promoting u^a to a Dynamical Khronon Field Solve the CMB Problem?

**Author**: Sheng-Kai Huang (with computational assistance)
**Date**: 2026-03-12
**Status**: Complete calculation with brutally honest verdict
**Purpose**: Determine whether identifying the tau framework's observer 4-velocity u^a with a dynamical Khronon/aether field can produce the c_s = 0 mode needed for CMB acoustic peaks.

---

## Executive Summary: The Verdict in Three Sentences

**The identification u^a = A^a (observer = aether) is structurally exact: the mathematical objects are identical.** Promoting u^a to dynamical does reproduce the Khronon/AeST field content and CAN achieve c_s = 0. **However, the tau framework does NOT constrain the free parameters -- it is a rewriting of AeST/Khronon, not a derivation from first principles.** The one genuinely new element: the tau framework provides a *physical interpretation* of why the aether field exists (it is the observer's time direction), and the DPI/Petz structure constrains the form of the action to be k-essence-like with a quadratic minimum -- but the numerical values of the coupling constants remain free.

**Classification: PARTIAL SUCCESS. Not a derivation, not a failure, but a physically motivated embedding with one constraining principle (DPI) that narrows the action space without fixing it completely.**

---

## Part I: The Identification u^a = A^a

### 1.1 The tau Framework's Observer

In Paper 2, the gravitational entropy production is:
```
Sigma_grav = -ln(-g_ab u^a u^b)
```
where u^a is the observer's 4-velocity. For a static observer at radius r in Schwarzschild:
```
u^a = (1/sqrt(-g_00), 0, 0, 0)
g_ab u^a u^b = -1   (by normalization)
Sigma_grav = -ln(1) = 0   ← WRONG if we use the normalized u^a
```

**Clarification**: The correct interpretation from Paper 2 is that Sigma_grav compares the observer's local metric to the asymptotic metric. More precisely:
```
Sigma_grav(r) = -ln(-g_00(r))
```
This is the entropy production experienced by an observer at r relative to infinity. The 4-velocity enters implicitly through the choice of time direction: g_00 = g_ab (dt)^a (dt)^b, where (dt)^a is the time basis vector aligned with the observer's worldline.

**The key object is not u^a itself but the 1-form aligned with it:**
```
n_a = -N partial_a T
```
where T is the time function defining the foliation and N = (-g^{ab} partial_a T partial_b T)^{-1/2} is the lapse. The unit normal to the T = const hypersurfaces is:
```
u_a = -N partial_a T,   u^a u_a = -1
```

### 1.2 The Einstein-Aether Field

In Einstein-aether theory (Jacobson & Mattingly 2001; Jacobson 2004, gr-qc/0410001), the aether is a unit timelike vector field A^a with:
```
g_ab A^a A^b = -1   (unit timelike constraint)
```

The action is:
```
S_ae = -(1/16piG) integral sqrt(-g) [R + K^{ab}_{cd} nabla_a A^c nabla_b A^d + lambda(g_ab A^a A^b + 1)] d^4x
```
where lambda is a Lagrange multiplier enforcing the unit constraint and:
```
K^{ab}_{cd} = c_1 g^{ab} g_{cd} + c_2 delta^a_c delta^b_d + c_3 delta^a_d delta^b_c + c_4 A^a A^b g_{cd}
```

The four dimensionless coupling constants c_1, c_2, c_3, c_4 parameterize all possible two-derivative terms for a unit vector field.

### 1.3 The Khronon Field

In the hypersurface-orthogonal limit (Jacobson 2010, arXiv:1001.4823; Blas, Pujolas & Sibiryakov 2010), the aether is:
```
A_a = -N partial_a T
```
where T is a scalar field (the "Khronon") and N = (-g^{ab} partial_a T partial_b T)^{-1/2}.

This is IDENTICAL to the tau framework's construction:
```
tau framework:    u_a = -N partial_a T    (observer's unit normal to T = const)
Khronon theory:   A_a = -N partial_a T    (aether aligned with T-foliation)
```

**THE IDENTIFICATION IS EXACT:**

| tau framework | Khronon/AeST | Mathematical object |
|---|---|---|
| Observer 4-velocity u^a | Aether field A^a | Unit timelike vector |
| Time function T | Khronon scalar T | Scalar field foliating spacetime |
| Lapse N | Lapse N | N = (-g^{ab} dT_a dT_b)^{-1/2} |
| Normalization u^a u_a = -1 | Constraint A^a A_a = -1 | Same equation |
| "Observer selects time direction" | "Aether defines preferred frame" | Same physics, different words |

### 1.4 What the Identification Means Physically

In the tau framework (Paper 5), the observer O is characterized by the degrees of freedom they access. The observer's worldline defines a preferred time direction u^a. Different observers see different tau_O.

Promoting u^a to dynamical means: **the time direction is not chosen by an external observer but is determined by the dynamics of the theory itself.** The "observer" becomes a field that pervades all of spacetime -- exactly what the aether is.

This is a deep conceptual shift:
- **Before**: u^a is external (observer chooses it); tau_O depends on the choice
- **After**: u^a is dynamical (determined by field equations); the "observer field" is part of spacetime

The physical content: **the existence of a time direction is not a free choice but a dynamical consequence of the theory.** Time is not imposed; it emerges from the action.

---

## Part II: What Action Does the tau Framework Give for u^a?

### 2.1 The Naive Candidate: S = integral Sigma_grav sqrt(-g) d^4x

The simplest guess from the tau framework:
```
S_Sigma = integral Sigma_grav sqrt(-g) d^4x = integral [-ln(-g_ab u^a u^b)] sqrt(-g) d^4x
```

But with the constraint u^a u_a = -1, this gives Sigma_grav = -ln(1) = 0 identically. This is TRIVIAL -- it tells us nothing.

**The problem**: Sigma_grav = -ln(-g_00) is defined relative to a REFERENCE frame at infinity, not in terms of the local u^a alone. In the Khronon language, Sigma_grav is:
```
Sigma_grav = -ln(-g_00) = -ln(g_ab (partial t)^a (partial t)^b / (-1))
```
where (partial t)^a is the coordinate time direction, NOT the normalized u^a. This distinction is crucial.

### 2.2 The Correct Candidate: Kinetic Terms from nabla_a u^b

The tau framework's Sigma_grav encodes the POTENTIAL (how much information is lost at a given point). But a dynamical field also needs KINETIC terms (how the field propagates). These come from the derivatives of u^a:
```
nabla_a u^b = sigma_{ab} + (1/3) theta h_{ab} + omega_{ab} + a_a u_b
```
where:
- theta = nabla_a u^a (expansion)
- sigma_{ab} = nabla_{(a} u_{b)} - (1/3) theta h_{ab} + a_{(a} u_{b)} (shear)
- omega_{ab} = nabla_{[a} u_{b]} + a_{[a} u_{b]} (twist/vorticity)
- a^a = u^b nabla_b u^a (acceleration)
- h_{ab} = g_{ab} + u_a u_b (spatial projector)

### 2.3 Connection to Einstein-Aether Parameters

The Einstein-aether kinetic tensor K^{ab}_{cd} nabla_a A^c nabla_b A^d decomposes as:
```
K^{ab}_{cd} nabla_a u^c nabla_b u^d = c_theta theta^2 + c_sigma sigma^2 + c_omega omega^2 + c_a a^2
```
where:
```
c_theta = (c_1 + 3c_2 + c_3) / 3
c_sigma = (c_1 + c_3) / 2
c_omega = (c_1 - c_3) / 2
c_a = c_1 + c_4
```

Or inversely:
```
c_1 = c_sigma + c_omega
c_2 = 3 c_theta - c_sigma - c_omega       [note: c_2 = c_theta - c_sigma for the standard decomposition]
c_3 = c_sigma - c_omega
c_4 = c_a - c_sigma - c_omega
```

The abbreviation conventions used in the literature:
```
c_{13} = c_1 + c_3 = 2 c_sigma
c_{14} = c_1 + c_4 = c_a + c_omega - c_sigma ...
```

(Note: different papers use slightly different conventions. I will use the Jacobson (2004) conventions below.)

### 2.4 Can the tau Framework Constrain These Parameters?

**The question**: Does ANY principle from the tau framework (DPI, Petz optimality, JRSWW bound, Crooks relation) constrain c_1, c_2, c_3, c_4?

**Attempt 1: Data Processing Inequality (DPI)**

The DPI states: for any CPTP map N, D(rho || sigma) >= D(N(rho) || N(sigma)).

Applied to the gravitational channel: Sigma >= 0 always. This constrains:
```
Sigma_grav = -ln(-g_00) >= 0   =>   -g_00 <= 1   =>   gravitational redshift
```

For the aether kinetic term, the DPI requires that the aether field equation does not GENERATE information -- it can only preserve or lose it. This means the aether stress-energy must satisfy the null energy condition (NEC):
```
T_ab^(ae) k^a k^b >= 0   for all null k^a
```

The NEC for Einstein-aether theory requires (Jacobson 2004):
```
c_1 + c_4 >= 0   (c_a >= 0)
c_1 - c_4 >= 0   (for subluminal propagation -- but DPI doesn't require this)
```

**Verdict**: DPI constrains the SIGN of certain combinations but not the VALUES. It requires c_a >= 0 (energy condition), but nothing more specific.

**Attempt 2: Petz Recovery Optimality**

The Petz map is the UNIQUE optimal retrodiction. In the gravitational context, this means the gravitational channel has a unique optimal recovery. The Petz map's structure constrains the channel to be "sufficiency-saturating" -- the closest to reversible.

For the aether field, this translates to: **the kinetic term should be the one that makes the gravitational channel maximally recoverable.** In k-essence language, this means the kinetic function K(X) should have a MINIMUM at the background value X_0 -- because at the minimum, perturbations have c_s = 0, meaning no information propagates in the scalar sector, which is the maximum-recovery condition.

**This is a genuine constraint from the tau framework!**

The argument:
1. Petz optimality => the channel should be maximally recoverable
2. Maximal recoverability => minimal entropy production for perturbations
3. Minimal entropy production => perturbations should not propagate (c_s = 0)
4. c_s = 0 => the kinetic function has a quadratic minimum at the background

This selects the k-essence form:
```
L_kinetic ~ mu^2 (X - X_0)^2
```
where X = -(1/2) g^{ab} nabla_a T nabla_b T.

This is EXACTLY the Blanchet-Skordis (2024) Khronon kinetic term K(Q) = mu^2 (Q - 1)^2!

**Attempt 3: JRSWW Bound**

The JRSWW bound F >= exp(-Sigma/2) constrains the recovery fidelity. For the aether perturbations:
```
F_perturbation >= exp(-Sigma_perturbation / 2)
```

If the perturbation is non-propagating (c_s = 0), then Sigma_perturbation = 0 and F = 1 (perfect recovery). This is consistent with the aether perturbation acting as a CONSERVED quantity (like CDM density contrast, which grows without dissipation).

**Attempt 4: Crooks Fluctuation Theorem**

The Crooks relation P(Sigma)/P(-Sigma) = exp(Sigma) constrains the forward/reverse ratio. For the aether field, this gives the detailed balance condition for fluctuations. It does NOT constrain the coupling constants directly.

### 2.5 Summary: What tau Constrains

| Principle | Constraint on aether parameters | Strength |
|---|---|---|
| DPI | c_a >= 0 (energy condition) | WEAK (sign only) |
| Petz optimality | K(X) has quadratic minimum => c_s = 0 at background | STRONG (form of action) |
| JRSWW bound | Consistent with c_s = 0 giving F = 1 | CONSISTENT (not constraining) |
| Crooks | Detailed balance for fluctuations | WEAK (no parameter constraint) |

**The tau framework constrains the FORM of the kinetic term (quadratic minimum => c_s = 0) but NOT the numerical values of the coupling constants.**

---

## Part III: Propagation Speeds and the c_s = 0 Condition

### 3.1 Einstein-Aether Propagation Speeds

In Einstein-aether theory around flat spacetime, the propagation speeds of the three modes are (Jacobson & Mattingly 2004, gr-qc/0402005; Jacobson 2007):

**Spin-2 (tensor/gravitational waves):**
```
c_T^2 = 1 / (1 - c_{13})     where c_{13} = c_1 + c_3
```

**Spin-1 (vector modes):**
```
c_V^2 = (2c_1 - c_1^2 + c_3^2) / (2 c_{14} (1 - c_{13}))     where c_{14} = c_1 + c_4
```

**Spin-0 (scalar mode):**
```
c_S^2 = (c_{123}) (2 - c_{14}) / (c_{14} (1 - c_{13}) (2 + c_{13} + 3c_2))
```
where c_{123} = c_1 + c_2 + c_3.

### 3.2 GW170817 Constraint

The multimessenger observation of GW170817 + GRB 170817A constrains (LIGO-Virgo 2017):
```
|c_T / c - 1| < 3 x 10^{-15}
```

From the spin-2 speed formula:
```
c_T^2 = 1/(1 - c_{13}) = 1   =>   c_{13} = c_1 + c_3 = 0
```

**This is essentially exact: c_1 = -c_3.**

### 3.3 The c_s = 0 Condition

Setting c_S = 0 requires:
```
c_{123} = c_1 + c_2 + c_3 = 0
```

Combined with c_1 + c_3 = 0 (from GW170817):
```
c_2 = 0
c_1 = -c_3   (free parameter)
```

The remaining free parameters are c_1 and c_4 (or equivalently, c_a = c_1 + c_4).

### 3.4 The Khronon Limit

In the hypersurface-orthogonal limit (Khronon), the twist omega_{ab} = 0 identically. This eliminates the c_omega dependence, and the effective kinetic term has only three parameters (c_theta, c_sigma, c_a) or equivalently:
```
alpha = c_{13} / 2 = c_sigma
beta = c_2
lambda_ae = c_{14} / 2 = c_a / 2
```
(using Blas-Pujolas-Sibiryakov notation)

With c_{13} = 0 (GW170817) and c_2 = 0 (c_s = 0):
```
alpha = 0
beta = 0
lambda_ae = c_4 / 2   (one free parameter)
```

This is a ONE-PARAMETER family of theories, all with:
- c_T = c (gravitational waves at light speed)
- c_S = 0 (scalar mode non-propagating)
- One free parameter controlling the coupling strength

### 3.5 The Blanchet-Skordis Khronon Theory

Blanchet & Skordis (2024, JCAP 11, 040; arXiv:2404.06584) use a different parameterization. Their Khronon field tau has kinetic function:
```
K(Q) = mu^2 (Q - 1)^2    where Q = c sqrt(-g^{ab} partial_a tau partial_b tau)
```

On Minkowski background with tau = t/c (comoving time), Q = 1 and K = 0.

**The second-order action gives a scalar mode with dispersion relation omega = 0 (non-propagating).**

The parameter mu sets the mass scale of the Khronon field:
```
mu^{-1} ~ 223 kpc    (DBI model with lambda_D ~ 30)
or
mu^{-1} ~ 22.3 Mpc   (DBI model with lambda_D ~ 1)
```

The DBI-inspired version:
```
K(Q) = mu^2 (sqrt(1 + (Q-1)^2/epsilon^2) - 1)
```
reduces to the quadratic form for small perturbations and gives DBI-like behavior for large perturbations (MOND regime).

### 3.6 Connection to the tau Framework

**The Petz optimality argument (Section 2.4, Attempt 2) selects exactly the quadratic minimum form:**
```
K(Q) ~ (Q - Q_0)^2
```
with Q_0 = 1 (comoving frame is the equilibrium).

This is a genuine prediction of the tau framework: the kinetic function must have a quadratic minimum at the background value, which is the unique form consistent with Petz-optimal retrodiction (minimal information loss for perturbations).

**What tau does NOT predict:**
- The mass scale mu (could be anything from Planck to Hubble)
- The DBI completion (the behavior away from the minimum)
- The coupling to the scalar field Phi_s (if one is needed)

---

## Part IV: The Scalar Field Question

### 4.1 AeST Requires Both u^a AND Phi_s

The Skordis-Zlosnik AeST theory (2021, PRL 127, 161302) requires TWO fields beyond the metric:
1. A unit timelike vector A^a (the aether)
2. A scalar field phi (controls the MOND interpolating function)

The Khronon theory (Blanchet-Skordis 2024) packages these differently:
1. A scalar field tau (the Khronon, whose gradient gives u^a)
2. A kinetic function K(Q) that plays the role of AeST's scalar field dynamics

In the newer Khronon-Tensor theory (Blanchet-Skordis 2025, arXiv:2507.00912), the vector field reappears alongside the Khronon.

### 4.2 What is Phi_s in the tau Framework?

**Natural candidate: Phi_s = Sigma_grav = -ln(-g_00)**

Let us check whether this identification is consistent.

In AeST, the scalar phi satisfies:
```
Y = u^a nabla_a phi    (time derivative along the aether)
```
where Y is the key dynamical variable.

If phi = Sigma_grav = -ln(-g_00), then in the static Schwarzschild case:
```
phi = -ln(1 - r_s/r) ≈ r_s/r + (1/2)(r_s/r)^2 + ...

Y = u^a nabla_a phi = (1/sqrt(-g_00)) partial_t phi = 0    (static case)
```

In the static case, Y = 0 and the dynamics reduce to the spatial gradient of phi:
```
nabla_i phi = (r_s/r^2) hat{r}_i + O((r_s/r)^3)
```

The MOND limit of AeST requires:
```
nabla . [mu(|nabla phi|/a_0) nabla phi] = 4 pi G rho_b
```

With phi = -ln(-g_00), in the weak field:
```
|nabla phi| ≈ |nabla (r_s/r)| = GM/(c^2 r^2) = a_N/c^2
```

For this to reproduce MOND, we need the interpolating function mu to depend on a_N/a_0, which requires:
```
|nabla phi| / a_0 ~ a_N / (a_0 c^2)
```

**This gives the right FUNCTIONAL form but the wrong normalization.** The acceleration scale a_0 must be put in by hand; the tau framework does not predict it.

### 4.3 The Static Limit Check

In the quasi-static, weak-field limit, AeST reduces to (Skordis & Zlosnik 2024, arXiv:2305.07742):
```
nabla^2 Phi_N = 4 pi G rho_b + nabla . [mu(|nabla phi_s|) nabla phi_s]
```

where the second term provides the MOND-like modification.

If phi_s = Sigma_grav, then in the tau framework:
```
nabla phi_s = nabla[-ln(-g_00)] ≈ nabla(2Phi_N/c^2) = 2 nabla Phi_N / c^2
```

So |nabla phi_s| = 2 a_N / c^2, and the MOND interpolating function becomes:
```
mu(2 a_N / (c^2 * a_0^{phi}))
```

where a_0^{phi} is the MOND scale in phi-units. If a_0^{phi} = 2 a_0/c^2, this reproduces standard MOND. This is self-consistent.

**Verdict on Phi_s = Sigma_grav**: The identification is dimensionally and functionally consistent in the quasi-static limit. It naturally connects the tau framework's entropy production to AeST's MOND scalar. However, it does NOT predict the value of a_0.

### 4.4 Cosmological Perturbations of Phi_s

On FRW background:
```
Phi_s^{(0)} = Sigma_grav^{(0)} = -ln(-g_00^{FRW}) = -ln(1) = 0
```

Perturbation:
```
delta Phi_s = -2Phi    (from delta g_00 = -2Phi in Newtonian gauge)
```

The time derivative along the aether:
```
Y = u^a nabla_a delta Phi_s = -2 dot{Phi}
```

In the matter-dominated era, Phi = const, so Y = 0 and delta Phi_s is FROZEN. This is exactly CDM-like behavior: the density perturbation grows in place without oscillating.

**This is a promising result.** The perturbation of Sigma_grav on FRW naturally has the right qualitative behavior: it is frozen when potentials are constant (matter domination) and evolves only when potentials decay (radiation domination).

### 4.5 But Does It Have the Right Quantitative Behavior?

The critical question: does delta Phi_s provide the correct AMOUNT of gravitational potential to match CMB observations?

In LCDM, the gravitational potential is sourced by:
```
k^2 Phi = -4 pi G a^2 (rho_b Delta_b + rho_c Delta_c + rho_gamma Delta_gamma + rho_nu Delta_nu)
```

The CDM term rho_c Delta_c dominates after z_eq. If we replace CDM with the Khronon/Sigma_grav field, we need:
```
rho_{Sigma} ~ Omega_DM rho_crit ~ 0.26 rho_crit
```

**Where does this energy density come from?**

In the Khronon theory, the background energy density of the Khronon field on FRW is:
```
rho_K = K(Q_0) + K'(Q_0) * (...) + ...
```

For K(Q) = mu^2(Q-1)^2 with Q_0 = 1 on FRW background: K(1) = 0 and rho_K = 0 at background level!

This means the Khronon field has ZERO background energy density. Its energy resides entirely in the perturbations. This is fine for producing CDM-like perturbative behavior, but it means:

**The Khronon field does NOT contribute to Omega_m h^2 at the background level.**

**This is a problem.** To get z_eq ~ 3400, we need Omega_m h^2 ~ 0.14. With baryons alone, Omega_b h^2 ~ 0.022. The Khronon must provide the missing Omega_DM h^2 ~ 0.12 at the background level.

In AeST, this is achieved through a non-trivial background configuration of the scalar field that provides effective dark matter energy density. In the Blanchet-Skordis Khronon theory, the effective dark matter energy comes from higher-order terms in K(Q) that give a non-zero contribution to the stress-energy tensor on the cosmological background.

**In the tau framework, Sigma_grav = 0 on FRW background, giving zero background energy density.**

This is a quantitative failure: the tau framework's Sigma_grav cannot provide the background dark matter density needed for z_eq ~ 3400.

---

## Part V: CMB Perturbation Theory with Dynamical u^a

### 5.1 Background: u^a = (1,0,0,0), Phi_s = 0

On the flat FRW background:
```
ds^2 = -dt^2 + a(t)^2 delta_{ij} dx^i dx^j
u^a = (1, 0, 0, 0)    (comoving observer)
Sigma_grav = 0
Q = 1,   K(Q) = 0
```

The Friedmann equations are UNMODIFIED (the Khronon contributes nothing to the background). This means:
```
H^2 = (8piG/3)(rho_b + rho_gamma + rho_nu)    [NO dark matter term]
```

Without CDM at the background level, z_eq ~ 530 instead of 3400. **The sound horizon, peak positions, and ALL background-dependent CMB observables are wrong by 20-30%.**

### 5.2 Perturbation: delta u^a, delta Phi_s

In Newtonian gauge:
```
ds^2 = -(1+2Phi)dt^2 + a^2(1-2Psi)delta_{ij}dx^i dx^j
```

The perturbed aether:
```
u^a = ((1-Phi), v^i/a)
```
where v^i is the velocity perturbation. The unit constraint u^a u_a = -1 fixes the time component.

The perturbed Khronon: T = t + delta_T(x,t), giving:
```
delta Q = c * delta[-g^{ab} partial_a T partial_b T]^{1/2}
```

At linear order:
```
delta Q = Phi + dot{delta_T}    (schematically)
```

### 5.3 The Scalar Mode Equation

For the Khronon field with K(Q) = mu^2(Q-1)^2, the linearized equation of motion for the scalar perturbation is:

On Minkowski:
```
omega^2 delta_T = 0   =>   omega = 0   (non-propagating!)
```

On FRW with perturbations:
```
ddot{delta_T} + 3H dot{delta_T} = source(Phi, Psi)
```

The effective sound speed:
```
c_s^2 = 0   (at leading order)
c_s^2 = O(k/mu)^2   (at subleading order, k-dependent)
```

For k << mu (which covers all CMB scales if mu^{-1} ~ kpc), c_s^2 ~ 0 to excellent approximation.

**This is the key CMB result: the Khronon perturbation is effectively pressureless on all CMB scales.**

### 5.4 The Effective Stress-Energy of Khronon Perturbations

The Khronon's stress-energy tensor perturbation acts as an effective dark matter:
```
delta T^0_0 (Khronon) = rho_{eff} Delta_{eff}
delta T^i_j (Khronon) = -p_{eff} delta^i_j = 0    (pressureless!)
```

The effective equation of state and sound speed:
```
w_eff ~ 0    (dust-like)
c_s^2 ~ 0   (non-propagating)
c_vis^2 ~ 0  (no viscosity)
```

**These are EXACTLY the GDM conditions needed for CMB fitting (see paper3_CMB_without_DM.md, Section 5.1).**

### 5.5 The Background Energy Problem, Revisited

Even though the perturbations behave like CDM, the BACKGROUND is wrong:

**Without CDM background energy:**
```
z_eq ~ 530 (instead of 3400)
s_* ~ 190 Mpc (instead of 144 Mpc)
l_1 ~ 170 (instead of 220)
```

**This means promoting u^a to dynamical gets the perturbative behavior right (c_s = 0) but the background wrong (missing Omega_DM h^2 ~ 0.12).**

**Can the Khronon provide background energy?**

Yes, IF K(Q) has a non-trivial minimum that contributes a cosmological-constant-like or dust-like term. For example, in the DBI version:
```
K(Q) = mu^2 [sqrt(1 + (Q-1)^2/epsilon^2) - 1]
```

On the FRW background with Q = 1, K = 0. But if the FRW background is not Q = 1 but Q = Q_0 != 1 (which happens if the Khronon is not perfectly aligned with cosmic time), then:
```
rho_K = K(Q_0) > 0
```

This is how AeST and Khronon theories actually work: the cosmological background has a non-trivial scalar field profile that provides the effective dark matter density. **But this requires tuning the initial conditions or the function K(Q) to give Omega_DM h^2 ~ 0.12 -- a free parameter, not a prediction.**

---

## Part VI: Does tau Constrain (c_1, c_2, c_3)?

### 6.1 What GW170817 + c_s = 0 Already Give

From purely observational/theoretical consistency:

```
c_T = c:     c_1 + c_3 = 0     (GW170817)
c_S = 0:     c_1 + c_2 + c_3 = 0     (CMB requirement)

Combined:    c_2 = 0,   c_1 = -c_3
```

This leaves TWO free parameters: c_1 and c_4.

### 6.2 What the tau Framework Adds

**Constraint from Petz optimality (Section 2.4):**

The kinetic function must have a quadratic minimum at the background value Q_0 = 1. This selects:
```
K(Q) = mu^2 (Q - 1)^2 + O((Q-1)^3)
```

In terms of the Einstein-aether parameters for the Khronon limit:
```
K = (alpha + beta/2) theta^2 + alpha sigma^2     (using c_{13} = 0, c_2 = 0)
```

Wait -- with c_1 + c_3 = 0 and c_2 = 0, the kinetic term becomes:
```
K^{ab}_{cd} nabla_a u^c nabla_b u^d = c_1(nabla_a u^b nabla^a u_b - nabla_a u^b nabla^b u_a) + c_4 a^2
                                      = c_1 * 2 sigma^2 + c_4 a^2
```
(using the decomposition and the fact that c_1 = -c_3 makes the theta^2 coefficient vanish when c_2 = 0).

Actually, let me be more careful. With c_1 = -c_3, c_2 = 0:
```
c_theta = (c_1 + 3*0 + c_3)/3 = (c_1 + c_3)/3 = 0
c_sigma = (c_1 + c_3)/2 = 0
c_omega = (c_1 - c_3)/2 = c_1    (but in hypersurface-orthogonal case, omega = 0)
c_a = c_1 + c_4
```

So the kinetic term reduces to:
```
L_kinetic = c_a * a^2 = (c_1 + c_4) * (u^b nabla_b u^a)(u^c nabla_c u_a)
```

This is purely acceleration-squared! The expansion theta and shear sigma have vanishing coefficients.

In the Khronon parameterization, this corresponds to:
```
K(Q) ~ (c_1 + c_4) * (Q - 1)^2
```
with mu^2 = (c_1 + c_4) * M_P^2 / (some factor).

**The tau framework's Petz constraint says mu^2 > 0, i.e., c_1 + c_4 > 0.** This is just the positive energy condition. It does NOT fix the VALUE of c_1 + c_4.

### 6.3 Can DPI Fix the Remaining Parameters?

The DPI gives monotonicity: coarse-graining always increases Sigma. For the Khronon field, this means the evolution of the Khronon must not decrease the QRE. Formally:
```
D(rho(t_2) || sigma(t_2)) <= D(rho(t_1) || sigma(t_1))    for t_2 > t_1
```

This is automatically satisfied by any CPTP evolution (it IS the DPI). For the Khronon, it means the field dynamics must be dissipative or at worst conservative. Since the Khronon has c_s = 0 (non-propagating), it is conservative by default -- perturbations grow but don't oscillate or radiate.

**DPI is automatically satisfied and gives NO additional parameter constraint.**

### 6.4 Summary of Parameter Constraints

| Source | Constraint | Parameters fixed |
|---|---|---|
| GW170817 | c_1 + c_3 = 0 | 1 relation |
| CMB (c_s = 0) | c_1 + c_2 + c_3 = 0 | 1 more relation (=> c_2 = 0) |
| Petz optimality | K(Q) has quadratic minimum | Form of action (not values) |
| DPI | c_a >= 0 | Sign constraint |
| Energy conditions | Various | Inequalities only |
| tau framework total | c_2 = 0, c_1 = -c_3, c_4 > -c_1 | 2 values + 1 inequality |
| FREE parameters | c_1 and c_4 (or mu and c_4) | 2 undetermined |

**The tau framework, combined with observations, leaves a 2-parameter family.** This is the same as the Khronon theory after GW170817 + CMB constraints. **No additional predictive power is gained.**

---

## Part VII: Honest Assessment -- Derivation, Rewriting, or Failure?

### 7.1 Scorecard

| Question | Answer | Evidence |
|---|---|---|
| Is u^a = A^a structurally exact? | **YES** | Same mathematical object (unit timelike vector) |
| Does promoting u^a to dynamical give Khronon/AeST content? | **YES** | Identical field content and constraint |
| Can this achieve c_s = 0? | **YES** | omega = 0 for Khronon perturbation on Minkowski |
| Does tau DERIVE the c_s = 0 condition? | **PARTIALLY** | Petz optimality argues for quadratic minimum (=> c_s = 0), but this is heuristic, not rigorous |
| Does tau fix the coupling constants? | **NO** | Two free parameters remain (c_1, c_4 or mu, c_4) |
| Does tau predict Omega_DM h^2 ~ 0.12? | **NO** | Background energy density is a free parameter |
| Does tau predict a_0? | **NO** | MOND scale is not determined |
| Does tau give a NEW physical interpretation? | **YES** | "DM = dynamical observer field = price of time" |
| Is this more than relabeling? | **MARGINALLY** | Petz optimality selects the form of K(Q) but not its parameters |
| Does the CMB work? | **YES, given parameters** | Same as Khronon/AeST by construction |

### 7.2 The Three Categories

**What is DERIVED (genuinely new):**
1. The identification u^a = A^a from the observer-dependence of tau
2. The Petz optimality argument for K(Q) having a quadratic minimum (c_s = 0)
3. The physical interpretation: dark matter = price of the observer's time direction
4. The connection Phi_s = Sigma_grav in the static MOND limit

**What is REWRITTEN (same physics, different words):**
1. The Khronon/AeST action (reproduced, not derived)
2. The c_s = 0 mechanism (from the Khronon structure, not from QRE)
3. The CMB fitting (identical to Blanchet-Skordis 2024)
4. The parameter constraints (from GW170817 + CMB, not from tau)

**What FAILS (doesn't work):**
1. Predicting Omega_DM h^2 ~ 0.12 from the framework
2. Predicting a_0 from the framework
3. Fixing c_1 and c_4 from information-theoretic principles
4. Explaining the background dark matter energy density from Sigma_grav (which is zero on FRW)

### 7.3 The Critical Distinction: Embedding vs. Derivation

The tau framework **embeds** Khronon/AeST physics by providing:
- A physical motivation for the aether field (the observer's time direction)
- A principle that selects the form of the action (Petz optimality => quadratic minimum)
- An information-theoretic interpretation (dark matter = temporal information cost)

But it does **not derive** Khronon/AeST because:
- The coupling constants are not determined
- The background energy density is not predicted
- The MOND scale a_0 is not derived
- The transition from cosmological (CDM-like) to galactic (MOND-like) behavior is put in through K(Q), not derived from QRE

**Analogy**: The tau framework's relation to AeST is like thermodynamics' relation to statistical mechanics. Thermodynamics provides the framework (laws, principles, constraints) but does not derive specific material properties (heat capacity, viscosity). The tau framework provides the framework (DPI, Petz, JRSWW) but does not derive specific gravitational parameters (c_i, mu, a_0).

This is VALUABLE but INCOMPLETE.

---

## Part VIII: The Background Energy Problem -- Can It Be Solved?

### 8.1 The Core Issue

The deepest problem is not the perturbative behavior (which works) but the background:

```
FRW background: Sigma_grav = 0, K(Q_0=1) = 0, rho_Khronon = 0
Required: rho_DM = Omega_DM rho_crit ~ 0.26 * 9.5e-27 kg/m^3
```

Zero background energy means z_eq is wrong, and ALL CMB observables shift by 20-30%.

### 8.2 Possible Solutions

**Solution A: Non-trivial cosmological Khronon background**

If the Khronon field tau is not perfectly aligned with cosmic time t (i.e., Q_0 != 1 on FRW), then K(Q_0) != 0 and there is non-zero background energy. AeST achieves this: the scalar field phi has a non-trivial time evolution on FRW that provides rho_DM.

**In the tau framework**: This would mean that the observer's time direction u^a does not exactly coincide with the cosmic rest frame. The misalignment provides an effective energy density. Physically: "the observer's clock runs slightly differently from cosmic time, and this discrepancy IS dark matter."

**Problem**: The misalignment Q_0 - 1 must be tuned to give Omega_DM h^2 = 0.12. This is a fine-tuning problem, not a prediction.

**Solution B: Additional contribution from Sigma_grav perturbation backreaction**

At second order in perturbation theory, the spatial average of Sigma_grav^2 terms could provide an effective background energy density:
```
<rho_Sigma^{(2)}> = <(delta Sigma_grav)^2> * (some coefficient)
```

With delta Sigma_grav = -2Phi and <Phi^2> ~ 10^{-10} (from primordial perturbations), this gives:
```
<rho_Sigma^{(2)}> / rho_crit ~ 10^{-10} << Omega_DM ~ 0.26
```

**Far too small by 9 orders of magnitude.** Backreaction cannot provide the dark matter density.

**Solution C: Topological/global contribution**

If the Khronon field has non-trivial global structure (e.g., winding number), the topological contribution could provide an effective background energy. This is speculative and has no calculation to support it.

**Solution D: Accept that the background requires an input**

The honest position: the tau framework provides the FORM of the action (Petz optimality => c_s = 0) and the INTERPRETATION (dark matter = observer dynamics), but the AMPLITUDE (Omega_DM h^2 = 0.12) is an input parameter, not a prediction.

This is analogous to how the Standard Model does not predict the Higgs VEV v = 246 GeV -- it is a free parameter that determines the scale of electroweak symmetry breaking.

### 8.3 The Cosmic Coincidence

There is one tantalizing hint. The dark matter density parameter:
```
Omega_DM h^2 ~ 0.12
```

And the MOND acceleration scale:
```
a_0 ~ c H_0 / (2 pi) ~ 1.2 x 10^{-10} m/s^2
```

These are connected through:
```
Omega_DM ~ a_0 / (c H_0) ~ O(1)
```

If the tau framework could derive the relation a_0 ~ c H_0 (which Milgrom has long argued is a cosmic coincidence or a deep fact), it would simultaneously explain both Omega_DM and a_0.

**The Petz bound F >= exp(-Sigma/2) with Sigma ~ H * L_Hubble / c ~ 1 gives tau ~ 0.4, which is O(1).** This might be connected to Omega_DM ~ O(1), but the connection is not rigorous.

---

## Part IX: Comparison with Blanchet-Skordis Khronon-Tensor Theory (2025)

### 9.1 The Khronon-Tensor Extension

Blanchet & Skordis (2025, arXiv:2507.00912) extend the 2024 Khronon theory by adding a vector (tensor) field B^a alongside the Khronon scalar tau. This is motivated by:
1. The need for a non-trivial cosmological background (Khronon alone gives zero background DM energy)
2. The desire to match both MOND and CMB precisely
3. Improved stability properties

The theory has:
- Khronon scalar tau (provides the time foliation and c_s = 0)
- Tensor/vector field B^a (provides the background DM energy density)

### 9.2 What This Means for the tau Framework

The Khronon-Tensor theory suggests that u^a ALONE is not enough. You also need a vector field B^a that provides the background energy. In the tau framework, this could be:

- u^a = observer's 4-velocity (the Khronon, providing c_s = 0 perturbations)
- B^a = spatial direction in the observer's reference frame (the additional vector, providing background energy)

Or alternatively:
- u^a = the Khronon-derived aether
- Sigma_grav = the scalar field that provides the MOND/CDM interpolation
- A separate mechanism for background energy

**The fact that state-of-the-art Khronon theories require ADDITIONAL fields beyond the Khronon itself suggests that the simple identification u^a = aether is necessary but NOT sufficient.**

---

## Part X: Final Verdict

### 10.1 Does tau -> Khronon Work?

**As a DERIVATION**: NO. The tau framework does not derive the Khronon action, coupling constants, or background energy density from information-theoretic principles. It provides a physical interpretation and a heuristic argument (Petz optimality => c_s = 0) but not a derivation.

**As a REWRITING**: MOSTLY. The identification u^a = A^a is exact, and the tau framework's language can describe the Khronon/AeST physics. But the translation adds no new parameters or predictions.

**As a PHYSICALLY MOTIVATED EMBEDDING**: YES, with genuine value. The tau framework:
1. Explains WHY there is an aether field (it is the observer's time direction)
2. Argues WHY c_s = 0 (Petz-optimal retrodiction at background)
3. Connects the CMB dark matter to the galactic MOND through Sigma_grav
4. Provides the "dark matter = price of time" interpretation

### 10.2 What Would Turn This into a Derivation

The embedding would become a derivation if ANY of these could be shown:

1. **D(g || G) perturbation theory on FRW naturally produces a Khronon-like field** with the correct kinetic term -- showing that the observer field EMERGES from the QRE structure rather than being put in.

2. **The JRSWW bound F >= exp(-Sigma/2) applied at the cosmological horizon determines Omega_DM** -- connecting the Petz fidelity to the dark matter density.

3. **The Crooks fluctuation theorem applied to cosmological expansion determines a_0** -- connecting the MOND scale to the arrow of time.

4. **The constrained optimization of D(rho_st || rho_m) produces the Khronon-Tensor field content** as Lagrange multipliers -- deriving the additional vector field B^a.

**None of these calculations exist.** They are the research program for Paper 4 (and possibly beyond).

### 10.3 The Honest One-Paragraph Summary

The tau framework's observer 4-velocity u^a is structurally identical to the Khronon/aether field. Promoting it to dynamical reproduces the field content of theories that successfully fit the CMB (AeST, Khronon theory). The Petz optimality principle provides a heuristic argument for why the scalar perturbation should have c_s = 0. However, the framework does not derive the coupling constants, the background dark matter density, or the MOND acceleration scale -- these remain free parameters, just as in Khronon/AeST. The deepest contribution is conceptual: it provides a physical interpretation of the aether field as the observer's time direction, making "dark matter" the dynamical cost of temporal ordering. Whether this interpretation can be promoted to a genuine derivation depends on calculations (constrained QRE perturbation theory, cosmological Petz bound) that do not yet exist.

### 10.4 Classification for the Paper 4 Architecture

```
Paper 4 Section on CMB:

PROVED:
- u^a ≡ A^a (structural identification, exact)
- Hypersurface orthogonal u^a ≡ Khronon field (exact)
- c_s = 0 follows from K(Q) = mu^2(Q-1)^2 form (known result, Blanchet-Skordis 2024)

NEW SYNTHESIS:
- Petz optimality => quadratic minimum of K(Q) => c_s = 0 (new argument)
- Phi_s = Sigma_grav identification in static MOND limit (new connection)
- "Dark matter = price of time" interpretation (new)

CONJECTURED:
- D(g||G) perturbation theory produces Khronon-like field
- Omega_DM h^2 determined by cosmological Petz bound
- a_0 determined by Crooks theorem at Hubble scale

FAILED/OPEN:
- Background dark matter density from Sigma_grav (zero on FRW)
- Fixing coupling constants from DPI/JRSWW
- Complete derivation of AeST from QRE
```

---

## Part XI: Implications for the Four-Paper Architecture

### 11.1 What This Calculation Changes

1. **Paper 3 (rotation curves)**: UNCHANGED. Running G at galactic scales is independent of the CMB problem.

2. **Paper 4 (grand unification)**: The CMB section should present the u^a = A^a identification as a structural result, the Petz optimality argument as motivation for c_s = 0, and honestly acknowledge the gap (background energy, coupling constants).

3. **The overall narrative**: The tau framework provides a UNIFYING INTERPRETATION (dark matter = observer's time direction = entropy production) but does NOT derive dark matter from first principles. This is valuable -- it connects galactic MOND, CMB peaks, and the arrow of time under one conceptual umbrella -- but it falls short of the "one equation, different solutions" dream.

### 11.2 Recommended Framing for Papers

**For Paper 3**: "The tau framework explains galactic-scale dark matter phenomenology through running G. The cosmological origin of the dark component is structurally compatible with the Khronon/AeST framework through the identification u^a = A^a, but requires additional work to derive from information-theoretic principles (Paper 4)."

**For Paper 4**: "The observer's 4-velocity u^a, when promoted to a dynamical field, is mathematically identical to the aether/Khronon field of theories that successfully fit the CMB. We present a Petz-optimality argument for why the kinetic function must have a quadratic minimum (c_s = 0). The background dark matter density and coupling constants remain free parameters in the current framework. We identify three concrete calculations that could close this gap."

---

## References

### Einstein-Aether Theory
- Jacobson & Mattingly 2001: arXiv:gr-qc/0007031 (gravity with dynamical preferred frame)
- Jacobson 2004: arXiv:gr-qc/0410001 (Einstein-aether theory review)
- Jacobson & Mattingly 2004: arXiv:gr-qc/0402005 (Einstein-aether waves -- propagation speeds)
- Jacobson 2010: arXiv:1001.4823 (extended Horava gravity and Einstein-aether theory)

### Khronon Theory
- Blas, Pujolas & Sibiryakov 2010: arXiv:0909.3525 (consistent extensions of Horava gravity)
- Blanchet & Skordis 2024: arXiv:2404.06584, JCAP 11, 040 (relativistic Khronon theory)
- Blanchet & Skordis 2025: arXiv:2507.00912 (Khronon-Tensor theory)

### AeST Theory
- Skordis & Zlosnik 2021: arXiv:2007.00082, PRL 127, 161302 (AeST -- CMB fit)
- Skordis & Zlosnik 2021: arXiv:2109.13287 (AeST linear stability)
- Skordis & Zlosnik 2024: arXiv:2305.07742 (AeST quasi-static limit)

### GW170817 Constraints
- LIGO-Virgo 2017: PRL 119, 161101 (GW170817 detection)
- Abbott et al. 2017: PRL 119, 251301 (strong constraints from GW170817)
- Oost, Mukohyama & Wang 2018: arXiv:1802.04303 (constraints on Einstein-aether after GW170817)

### tau Framework
- Paper 1: Huang (2026) -- Petz recovery unification
- Paper 2: Huang (2026) -- Information loss is local: gravitational refractive index
- Paper 5: Huang (2026) -- Observer-dependent temporal asymmetry

### Quantum Reference Frames
- De Vuyst et al. 2025: arXiv reference (quantum reference frames for general symmetry groups)
- Giacomini, Castro-Ruiz & Brukner 2019: Nat. Commun. 10, 494 (QRF and equivalence principle)

### k-Essence and Dust from Scalar Fields
- Scherrer 2004: PRL 93, 011301 (purely kinetic k-essence as unified DM)
- Kunz 2009: arXiv:0702615 (generalized dark matter)

### Other
- Dorau-Much 2025: arXiv:2510.24491 (QRE -> Einstein equations)
- Bianconi 2025: PRD 111, 066001 (gravity from entropy)

---

*Last updated: 2026-03-12*
*This calculation represents the CRITICAL assessment of whether the tau framework can address the CMB acoustic peak problem through the Khronon/aether route.*
*Classification: PARTIAL SUCCESS -- structural embedding with genuine physical insight, but not a first-principles derivation.*
