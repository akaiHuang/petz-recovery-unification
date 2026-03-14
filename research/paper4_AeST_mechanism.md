# Deep Technical Analysis: How AeST Achieves c_s = 0 Without CDM Particles

**Author**: Sheng-Kai Huang
**Date**: 2026-03-12
**Status**: Deep technical research for Paper 4
**Purpose**: Understand the PRECISE mechanism by which AeST (Aether Scalar Tensor theory) reproduces CDM-like behavior, and evaluate the connection to the tau framework

---

## Table of Contents

1. [Executive Summary](#1-executive-summary)
2. [The Complete AeST Action](#2-the-complete-aest-action)
3. [How c_s = 0 Is Achieved: The Precise Mechanism](#3-how-c_s--0-is-achieved-the-precise-mechanism)
4. [The Khronon Theory: A Simpler Realization](#4-the-khronon-theory-a-simpler-realization)
5. [CMB Fit Results](#5-cmb-fit-results)
6. [Connection to the tau Framework](#6-connection-to-the-tau-framework)
7. [Recent Developments (2023-2026)](#7-recent-developments-2023-2026)
8. [Honest Assessment: tau -> AeST](#8-honest-assessment-tau---aest)
9. [Key References](#9-key-references)

---

## 1. Executive Summary

### The Core Question

How does AeST (Skordis & Zlosnik 2021, PRL 127, 161302) produce a mode with zero sound speed (c_s = 0) that mimics cold dark matter at the perturbation level, enabling CMB fits comparable to LCDM -- without introducing dark matter particles?

### The Answer in One Paragraph

AeST introduces a unit timelike vector field A^mu (the "aether") and a shift-symmetric scalar field phi. The scalar field's kinetic structure is divided into two parts: Y = D^mu phi D_mu phi (spatial gradients, projected orthogonal to A^mu) and Q = A^mu nabla_mu phi (time derivative along the aether). The function F(Y,Q) in the action is chosen so that perturbations of the scalar field around the cosmological background have a dispersion relation omega = 0 -- the mode is NON-PROPAGATING. This happens because the kinetic function has a minimum at Q = Q_0 (the background value), making the quadratic perturbation of Q a STIFF POTENTIAL rather than a kinetic term. The scalar perturbation does not oscillate -- it grows in place, exactly like pressureless dust (CDM). The aether's unit-timelike constraint (A^mu A_mu = -1) breaks Lorentz invariance in the gravitational sector, which is the essential ingredient that allows c_s = 0 without introducing particles.

### The Bottom Line for Paper 4

| Question | Answer |
|----------|--------|
| Is tau -> AeST a derivation? | **No.** Not currently. |
| Is it a motivation? | **Yes.** The observer u^mu maps naturally to the aether A^mu. |
| Is it just notation? | **More than notation, less than derivation.** The structural parallel is deep but the mapping is not proven. |
| What new prediction would tau give? | The density Omega_DM h^2 = 0.12 from a QRE extremization condition. |
| What is the key gap? | Deriving F(Y,Q) from D(rho_spacetime || rho_matter). |

---

## 2. The Complete AeST Action

### 2.1 Full Action

From Skordis & Zlosnik (2021, arXiv:2007.00082; see also arXiv:2109.13287 and arXiv:2412.15395):

```
S = S_grav + S_matter

S_grav = integral d^4x sqrt(-g) / (16 pi G_tilde) *
         [ R
           - (K_B/2) F_{mu nu} F^{mu nu}
           + 2(2 - K_B) J^mu nabla_mu phi
           - (2 - K_B) Y
           - F(Y, Q)
           - lambda (A^mu A_mu + 1) ]

S_matter = S_m[psi, g_{mu nu}]
```

### 2.2 Field Content

| Field | Symbol | Nature | Role |
|-------|--------|--------|------|
| Metric | g_{mu nu} | Rank-2 tensor | Standard spacetime geometry |
| Aether | A^mu | Unit timelike vector | Selects preferred time direction |
| Scalar | phi | Shift-symmetric scalar | Sources MOND + CDM-like perturbations |
| Lagrange multiplier | lambda | Scalar | Enforces unit-norm constraint |

### 2.3 Key Definitions

**Aether constraint:**
```
g_{mu nu} A^mu A^nu = -1     (unit timelike)
```

**Vector field strength (Maxwell-like):**
```
F_{mu nu} = nabla_mu A_nu - nabla_nu A_mu
```

**Aether acceleration:**
```
J^mu = A^nu nabla_nu A^mu     (4-acceleration of aether worldlines)
```

**Spatial projection of scalar gradient:**
```
D_mu phi = (delta^nu_mu + A^nu A_mu) nabla_nu phi
         = nabla_mu phi + A_mu (A^nu nabla_nu phi)

Y = D^mu phi D_mu phi = g^{mu nu} D_mu phi D_nu phi
```

**THIS IS CRITICAL**: Y is the kinetic term of phi PROJECTED ORTHOGONAL to the aether. It measures the spatial gradient of phi relative to the preferred frame defined by A^mu.

**Time derivative of scalar along aether:**
```
Q = A^mu nabla_mu phi     (proper-time derivative of phi)
```

**THIS IS CRITICAL**: Q is the time derivative of phi in the aether frame. The function F(Y,Q) controls the dynamics through these two invariants.

### 2.4 Parameters

| Parameter | Role | Constraint |
|-----------|------|------------|
| G_tilde | Bare gravitational constant | Related to G_N by field redefinitions |
| K_B | Aether kinetic coupling | 0 < K_B < 2 (stability) |
| Q_0 | Background value of Q in cosmology | Set by cosmological evolution |
| lambda_s | Strong-field scalar coupling | lambda_s >= 0 (stability) |
| K_2 | Quadratic coupling around Q_0 | K_2 > 0 (ghost-free) |
| Free function F(Y,Q) | Controls MOND transition | Must satisfy multiple limits |
| a_0 | MOND acceleration scale | ~ 1.2 x 10^{-10} m/s^2 (empirical) |

**Total free parameters: ~8-10**, depending on how F(Y,Q) is parameterized.

### 2.5 Gravitational Wave Speed

AeST is specifically constructed so that tensor perturbations (gravitational waves) propagate at the speed of light:
```
c_T = c     (exact, by construction)
```

This satisfies the GW170817 constraint |c_T/c - 1| < 10^{-15}. Earlier theories (TeVeS) failed this test. AeST was designed from the start to satisfy it. The K_B parameter controls the vector kinetic term; the specific combination in the action ensures c_T = c regardless of K_B's value.

---

## 3. How c_s = 0 Is Achieved: The Precise Mechanism

### 3.1 The Cosmological Background

On an FLRW background with scale factor a(t):
- The aether aligns with the cosmic time direction: A^mu = (1/N, 0, 0, 0) where N is the lapse
- The scalar field depends only on time: phi = phi(t)
- The spatial gradient vanishes: Y = 0
- Only Q = A^mu nabla_mu phi = phi_dot/N survives

The background dynamics are governed by F(0, Q), which we write as F(Q) for brevity.

**Friedmann equation:**
```
H^2 + k/a^2 = (8 pi G_tilde / 3) rho_matter - (1/3)(F - Q F_Q) + Lambda/3
```

**Scalar field equation:**
```
d(F_Q)/dt + 3H F_Q = 0
```

where F_Q = dF/dQ.

**The scalar field equation has the immediate solution:**
```
F_Q = C / a^3     (C = constant)
```

This is the KEY: the "momentum" F_Q of the scalar field dilutes as a^{-3} -- exactly like the energy density of pressureless dust. The scalar field phi, through its kinetic coupling to the aether, behaves as a DUST-LIKE fluid.

### 3.2 The Effective Fluid Description

Define the effective energy density and pressure of the scalar-aether system:

```
rho_phi = -(1/(8 pi G_tilde)) (F - Q F_Q)
p_phi   =  (1/(8 pi G_tilde)) F
```

The equation of state:
```
w_phi = p_phi / rho_phi = F / (Q F_Q - F)
```

**For CDM-like behavior, we need w_phi ~ 0, which requires |F| << |Q F_Q - F|, i.e., |F| << |Q F_Q|.**

### 3.3 The c_s = 0 Mechanism: Step by Step

**Step 1: Choose F(Q) with a minimum at Q = Q_0.**

The function F(Q) is chosen to have the form (near Q_0):
```
F(Q) ~ -2 K_2 (Q - Q_0)^2 + higher order
```

At the minimum, F(Q_0) ~ 0 and F_Q(Q_0) = 0.

**Step 2: Linearize around the cosmological background.**

Write Q = Q_bg(t) + delta_Q(x, t). The perturbation equation for delta_Q follows from the second variation of the action.

**Step 3: The dispersion relation.**

On a Minkowski background (ignoring expansion for clarity), the linearized scalar perturbation has:

```
omega^2 = 0     (for the scalar mode)
```

This is NOT omega^2 = c_s^2 k^2 (which would give c_s > 0 for a standard scalar). The omega = 0 dispersion relation means the mode is NON-PROPAGATING.

**Step 4: Why omega = 0?**

The key is the structure of the kinetic term. For a standard scalar field with Lagrangian L = (1/2)(nabla phi)^2 - V(phi), the perturbation has omega^2 = k^2 + V''(phi). The spatial kinetic term k^2 comes from the spatial gradient (nabla phi)^2.

In AeST, the spatial and temporal kinetic terms are SEPARATELY controlled:
- Temporal: through Q = A^mu nabla_mu phi (the function F(Y=0, Q) controls temporal dynamics)
- Spatial: through Y = D^mu phi D_mu phi (the function F(Y, Q=Q_0) controls spatial dynamics)

When Y = 0 at the background level (no spatial gradients in FLRW), and F is expanded around this background, the spatial kinetic term for perturbations comes from:
```
delta^2 S / delta(nabla phi)^2 ~ partial F / partial Y |_{Y=0}
```

If the function F(Y, Q) is chosen such that its Y-dependence is MILD (specifically, (2-K_B)lambda_s Y with lambda_s controlling the coefficient), then the effective spatial kinetic energy for the scalar perturbation can be made ZERO or negligibly small. Without a spatial kinetic term, perturbations have no mechanism to propagate -- their sound speed vanishes.

**Step 5: Physical interpretation.**

The aether A^mu provides a preferred frame. In this frame, the scalar field phi has dynamics controlled by its time derivative Q but NOT by its spatial gradient Y (at the perturbative level). This is like a collection of dust particles that each evolve independently without communicating with neighbors -- they "fall" in the gravitational potential but don't generate pressure waves.

### 3.4 The Lorentz Invariance Breaking Is Essential

**Theorem (informal):** No Lorentz-invariant scalar field theory can have c_s = 0 perturbations around a homogeneous background.

**Proof sketch:** In a Lorentz-invariant theory, the kinetic term is (nabla phi)^2 = -(phi_dot)^2 + (nabla_spatial phi)^2. You cannot separate the temporal and spatial parts. The sound speed is:
```
c_s^2 = (dP/drho)|_entropy = (partial L / partial X) / (partial L / partial X + 2X partial^2 L / partial X^2)
```
where X = -(1/2)(nabla phi)^2. For c_s = 0, you need partial L / partial X = 0 at the background, which typically gives a degenerate kinetic structure.

In AeST, the aether A^mu BREAKS Lorentz invariance, allowing Y (spatial) and Q (temporal) to be independent variables. The function F(Y, Q) can have different dependence on each, enabling:
- F_Q ≠ 0 (non-trivial temporal dynamics)
- F_Y ~ 0 or F_Y small (suppressed spatial propagation)
- Result: c_s ~ 0

**This is the deep reason:** c_s = 0 requires Lorentz violation. AeST provides it through the aether. The tau framework would need to provide it through the observer's worldline.

### 3.5 The GDM Reduction

At cosmological scales, AeST perturbations reduce to the Generalized Dark Matter (GDM) parameterization (Kunz 2009, arXiv:0702615):

| GDM parameter | AeST value | CDM value | Meaning |
|---------------|------------|-----------|---------|
| w | ~ 0 (time-dependent) | 0 (exact) | Equation of state |
| c_s^2 | ~ 0 (k-dependent) | 0 (exact) | Sound speed squared |
| c_vis^2 | ~ 0 | 0 | Viscosity parameter |
| Omega h^2 | ~ 0.12 (input) | 0.12 (input) | Energy density |

**Important:** AeST does NOT predict Omega_DM h^2 = 0.12. This is an INPUT parameter, set by the initial conditions of the scalar field (the constant C in F_Q = C/a^3). The value 0.12 is fitted to observations, just as in LCDM.

### 3.6 Comparison with k-Essence

AeST's scalar sector is closely related to shift-symmetric k-essence:
```
L_k-essence = f(X)    where X = -(1/2) nabla_mu phi nabla^mu phi
```

A k-essence Lagrangian with a minimum at X = X_0 has c_s = 0 at the minimum:
```
c_s^2 = f'(X) / [f'(X) + 2X f''(X)] -> 0    when f'(X_0) = 0
```

The Khronon theory's K(Q) = mu^2 (Q-1)^2 is exactly this structure: it has a minimum at Q = 1, giving c_s = 0 at the minimum.

**AeST goes beyond k-essence** by having the SEPARATE spatial variable Y, controlled by the aether. This gives more freedom to ensure:
1. c_s = 0 for cosmological perturbations (via F(Y,Q))
2. MOND-like behavior for galactic dynamics (via J(Y) or the Y-dependence of F)
3. c_T = c for gravitational waves (via K_B)

---

## 4. The Khronon Theory: A Simpler Realization

### 4.1 The Khronon Action

Blanchet & Skordis (2024, JCAP 11, 040; arXiv:2404.06584):

```
S = (c^3 / (16 pi G)) integral d^4x sqrt(-g) [R - 2 J(Y) + 2 K(Q)] + S_matter
```

where:
- tau is the "Khronon" scalar field (foliates spacetime)
- n_mu = -(c/Q) nabla_mu tau is the unit normal to the foliation
- Q = c sqrt(-g^{mu nu} nabla_mu tau nabla_nu tau) (normalization)
- A_mu = c^2 n^nu nabla_nu n_mu = -c^2 q_mu^nu nabla_nu ln Q (acceleration)
- Y = A_mu A^mu / c^4 (acceleration-squared, dimensionless)
- q_mu^nu = delta_mu^nu + n_mu n^nu (spatial projector)

### 4.2 Key Simplification Over AeST

The Khronon theory achieves the SAME physics as AeST with fewer fields:
- **AeST**: metric + vector A^mu + scalar phi + Lagrange multiplier lambda (4 fields)
- **Khronon**: metric + scalar tau (2 fields)

The vector field is DERIVED from the scalar: n^mu = -(c/Q) g^{mu nu} nabla_nu tau. This eliminates the need for a separate aether field and its constraint.

### 4.3 How K(Q) Gives c_s = 0

**Baseline form:** K(Q) = mu^2 (Q - 1)^2

On the FLRW background, Q_bg ~ 1 (plus small corrections). Perturbations delta_Q around this background see the potential:
```
delta^2 K = 2 mu^2 (delta_Q)^2
```

This is a MASS TERM for the perturbation, not a kinetic term. The perturbation delta_Q is "frozen" by the stiff potential -- it does not propagate. Result:
```
omega = 0     (non-propagating mode)
```

**DBI variant:** K(Q) = mu^2 (sqrt(1 + (Q-1)^2/epsilon^2) - 1)

The DBI form provides better behavior across regimes and is preferred for CMB fitting.

### 4.4 How J(Y) Gives MOND

In the quasi-static, weak-field limit:
```
Y ~ |nabla Xi|^2 / c^4     where Xi = phi - sigma_dot + (1/2)|nabla sigma|^2
```

The function J(Y) is chosen so that at low accelerations (small Y):
```
J(Y) ~ Lambda - Y + c^2 a_0^{-1} Y^{3/2}
```

This gives the MOND interpolating function. The modified Poisson equation becomes:
```
nabla . [(1 + J_Y) nabla Xi] + mu^2 Xi = 4 pi G rho_m
```

where J_Y = dJ/dY. At low Y (low acceleration):
```
J_Y ~ -1 + (3/2) c^2 a_0^{-1} Y^{1/2} = -1 + (3/2)(|nabla Xi| / (c^2 a_0))
```

This reproduces the MOND potential equation:
```
nabla . [mu(|nabla Phi| / a_0) nabla Phi] = 4 pi G rho_b
```

where mu(x) = 1 + J_Y is the MOND interpolating function.

### 4.5 The Khronon Field as Time Itself

The scalar tau labels spacelike hypersurfaces -- it IS a time coordinate. The gradient nabla_mu tau defines the flow of time. In cosmology, tau aligns with cosmic time. In quasi-static situations, tau is the Newtonian universal time plus perturbative corrections.

**This is the deepest connection to the tau framework:** the Khronon field IS the physical realization of "time" in the tau framework's sense. The preferred foliation by tau is exactly the "observer's time" that defines Sigma = -ln(-g_00).

---

## 5. CMB Fit Results

### 5.1 AeST CMB Fit

From Skordis & Zlosnik (2021, PRL 127, 161302):

**Claim:** "A relativistic theory leading to Modified Newtonian Dynamics... [that] can lead to... a cosmological phenomenology that is consistent with observations of the cosmic microwave background and matter power spectrum."

**What this means concretely:**
- The TT, TE, EE CMB power spectra from AeST are visually indistinguishable from LCDM at the resolution of Planck
- The matter power spectrum P(k) reproduces the BAO wiggles and the correct shape
- No detailed chi^2/dof comparison has been published (this is a notable gap)

**How the fit is achieved:**
1. At the background level, the scalar field energy density rho_phi dilutes as a^{-3} (dust-like)
2. Set rho_phi to match Omega_DM h^2 ~ 0.12 -> this reproduces the correct z_eq ~ 3400
3. Perturbations of the scalar have c_s ~ 0 -> they grow like CDM perturbations
4. The combined effect on the Poisson equation, sound horizon, and damping tail is identical to CDM

### 5.2 Number of Free Parameters

| Parameter type | Count | Notes |
|---------------|-------|-------|
| Standard cosmological | 6 | Same as LCDM (H_0, Omega_b h^2, n_s, A_s, tau_reion, + implicit Omega_phi) |
| AeST-specific | 2-4 | K_B, mu (or Q_0), parameters of F(Y,Q) |
| MOND scale | 1 | a_0 (but this is fixed by galaxy data) |
| **Total** | **~8-10** | **vs LCDM's 6** |

The additional parameters are the PRICE of replacing CDM particles with fields. AeST has more parameters but covers more phenomenology (both CMB and galaxy rotation curves from the same theory).

### 5.3 What AeST Gets Right

1. **CMB TT spectrum**: All acoustic peaks in correct positions and heights
2. **CMB TE, EE spectra**: Correct polarization pattern
3. **Matter power spectrum**: Correct BAO scale and overall shape
4. **Galaxy rotation curves**: MOND limit reproduces flat rotation curves
5. **Gravitational wave speed**: c_T = c (satisfies GW170817)
6. **Solar system tests**: Reduces to GR in strong-field limit

### 5.4 What AeST Gets Wrong (or Hasn't Addressed)

1. **Galaxy clusters**: Negative effective mass density at certain radii (arXiv:2312.00889). The parameter m^2/f_G faces conflicting constraints from galactic and cluster scales (arXiv:2301.03499)
2. **Weak lensing at large radius**: Deviations from MOND at the radii probed by galaxy-galaxy lensing
3. **Detailed chi^2 comparison**: No published fit quality comparison with LCDM
4. **N-body structure formation**: No N-body simulation of AeST exists
5. **Bullet Cluster**: Not specifically addressed
6. **Non-linear regime**: Behavior during gravitational collapse, mergers unknown

---

## 6. Connection to the tau Framework

### 6.1 The Structural Parallel

| tau framework | AeST | Khronon | Match quality |
|---------------|------|---------|---------------|
| Observer's 4-velocity u^mu | Aether A^mu | Foliation normal n^mu | **EXACT structural match** |
| u^mu u_mu = -1 | A^mu A_mu = -1 | n^mu n_mu = -1 | **Identical constraint** |
| Sigma = -ln(-g_{00}) | F(Y,Q) | J(Y) + K(Q) | **Different but related** |
| tau = 1 - exp(-Sigma/2) | tau = 1 - F(rho, R_Petz) | tau = Khronon field | **Name collision!** |
| D(rho_spacetime \|\| rho_matter) | - | - | **No AeST analog** |
| Running G at galactic scales | MOND from J(Y) | MOND from J(Y) | **Same phenomenology, different mechanism** |

### 6.2 The Observer-Aether Identification

The single most important structural parallel:

**tau framework:** The observer with 4-velocity u^mu defines the reference frame in which Sigma = -ln(-g_ab u^a u^b) is computed. Different observers see different Sigma (Paper 5).

**AeST:** The aether A^mu defines a preferred reference frame. In this frame, the scalar field phi has its dynamics controlled by Q = A^mu nabla_mu phi.

**Identification:**
```
u^mu (observer) <-> A^mu (aether) <-> n^mu (Khronon normal)
```

In the tau framework, u^mu is treated as an EXTERNAL choice (the observer). In AeST, A^mu is a DYNAMICAL field (it has its own equation of motion from varying the action with respect to A^mu).

**The conceptual leap for Paper 4:** Promote the observer's 4-velocity from external to dynamical. The observer IS part of the physics, not outside it. This is exactly what AeST does with the aether.

### 6.3 Can Sigma Map to F(Y,Q)?

**Attempt 1: Direct identification Sigma_grav = -ln(-g_00) <-> Q = A^mu nabla_mu phi**

In the comoving frame (A^mu = (1,0,0,0)):
```
g_00 = -1 (FRW background)
Sigma = -ln(1) = 0

Q = phi_dot (background scalar time derivative)
```

These are different quantities. Sigma measures the gravitational redshift; Q measures the scalar field's rate of change. They do not directly map.

**Attempt 2: Sigma as the ACTION density, not a field variable**

In the tau framework, the proposed unified equation is:
```
Sigma = D(rho_spacetime || rho_matter)
```

In Bianconi's GfE, the gravitational action IS this QRE:
```
S_grav = integral sqrt(-g) D(rho_g || rho_m)
```

In AeST, the action is:
```
S_grav = integral sqrt(-g) [R - (K_B/2) F^2 + 2(2-K_B) J.nabla phi - (2-K_B)Y - F(Y,Q) - lambda(A^2+1)] / (16 pi G_tilde)
```

**The question:** Can the AeST action be rewritten as a QRE between spacetime and matter density matrices?

**Partial answer:** The Einstein-Hilbert term R can be derived from QRE via Dorau-Much (2025). The aether + scalar terms in AeST are ADDITIONAL to R. They represent degrees of freedom that GR does not have.

**Key insight:** If D(rho_spacetime || rho_matter) is computed for a spacetime WITH a preferred foliation (i.e., rho_spacetime includes the information about the foliation), then the additional terms might emerge from the foliation structure's contribution to the QRE.

### 6.4 The Natural Candidate: phi_s <-> Sigma_grav?

Consider the identification:
```
phi_s (AeST scalar) <-> Sigma_grav = -ln(-g_ab u^a u^b)
```

**In the static weak-field limit:**
```
g_00 = -(1 + 2 Phi_N)    (Newtonian gauge)
Sigma_grav = -ln(1 + 2 Phi_N) ~ -2 Phi_N

AeST scalar in quasi-static limit: Xi ~ Phi_N + corrections
```

So Sigma_grav ~ -2 Xi in weak field. This is a proportionality, not an identity, but it shows that the tau framework's Sigma_grav and AeST's scalar Xi carry SIMILAR information (both encode the gravitational potential).

**Does this mapping work quantitatively?**

At the perturbative level around FRW:
```
delta_Sigma = -2 Phi     (from tau framework)
delta_Xi = Phi + corrections     (from AeST quasi-static limit)
```

So delta_Sigma ~ -2 delta_Xi. The spatial gradient of Sigma maps to the spatial gradient of Xi (up to a factor of -2). The MOND equation in Khronon language:
```
nabla . [(1 + J_Y) nabla Xi] + mu^2 Xi = 4 pi G rho_m
```

Would translate to:
```
nabla . [(1 + J_Y) nabla(-Sigma/2)] + mu^2 (-Sigma/2) = 4 pi G rho_m
```

which is a modified Poisson equation for Sigma. This is consistent with the tau framework's structure but does NOT derive the MOND interpolating function J_Y from QRE principles.

### 6.5 What Would tau Add to AeST?

**Currently, AeST has these unexplained features:**

1. **The free function F(Y,Q) is ad hoc.** It is chosen to interpolate between MOND and GR. The tau framework could potentially derive F from the structure of D(rho_spacetime || rho_matter).

2. **Omega_DM h^2 = 0.12 is an input.** AeST does not predict the dark matter density; it is fitted to observations. The tau framework could potentially derive it from a QRE extremization condition (e.g., maximizing D or minimizing Sigma under constraints).

3. **The aether has no physical interpretation.** AeST introduces A^mu as a mathematical field with no deeper meaning. The tau framework identifies it as the observer's worldline -- time itself. This provides a physical interpretation: "dark matter" is the energetic cost of maintaining a time direction.

4. **The connection to quantum information is absent.** AeST is purely classical. The tau framework connects gravity to QI through Petz recovery, QRE, and DPI. This could provide quantum corrections to AeST and predictions for Planck-scale modifications.

5. **The acceleration scale a_0 is empirical.** The tau framework might derive a_0 from quantum information principles (e.g., a_0 = c H_0 from the cosmological QRE horizon).

### 6.6 The Dream Derivation for Paper 4

Starting from S_grav = integral sqrt(-g) D(rho_spacetime || rho_matter), show that:

1. At zeroth order in perturbation theory, the Einstein equations emerge (Dorau-Much)
2. The QRE structure naturally includes a preferred foliation (the modular flow generates time evolution)
3. The foliation normal n^mu = the aether A^mu = the observer u^mu
4. The modular Hamiltonian K = -ln(rho) generates a scalar field phi (the Khronon)
5. Perturbations of phi around the thermal background have omega = 0 (from the constrained nature of modular flow)
6. The effective fluid has w ~ 0, c_s ~ 0 (CDM-like)
7. The MOND interpolating function follows from the DPI structure at galactic scales

**Step 5 is the crucial one.** Why would the modular Hamiltonian perturbation have omega = 0? In standard QFT, modular flow generates a THERMAL state, and thermal perturbations propagate at the speed of sound of the thermal medium. For a relativistic thermal medium, c_s = c/sqrt(3). This is NOT zero.

**The escape route:** If the modular flow generates a CONSTRAINED perturbation (a Lagrange multiplier enforcing the first law of thermodynamics at each point), then the perturbation is non-propagating (omega = 0). This is the same mechanism as AeST: the constraint (A^mu A_mu = -1 or Q ~ Q_0) kills the propagating mode.

**Status:** This is a well-posed mathematical question but no calculation exists. It requires:
1. Computing the modular Hamiltonian for a QFT on FRW background
2. Perturbing the modular flow equations
3. Identifying the constraint structure
4. Determining the dispersion relation of constraint-derived modes

---

## 7. Recent Developments (2023-2026)

### 7.1 Blanchet & Skordis: Khronon Theory (2024)

**arXiv:2404.06584, JCAP 11 (2024) 040**

A simpler version of AeST using only a Khronon scalar field tau. Key improvements:
- Fewer degrees of freedom (eliminates independent vector field)
- DBI-type kinetic function naturally bridges MOND and cosmological regimes
- Claims CMB compatibility at linear perturbation order
- omega = 0 dispersion relation confirmed

### 7.2 Blanchet & Skordis: Khronon-Tensor Extension (2025)

**arXiv:2507.00912**

Extends the Khronon theory to include tensor modifications. Key developments:
- The action S = (c^3/16piG) integral sqrt(-g) [R - 2J(Y) + 2K(Q)] confirmed
- DBI-inspired K(Q) preferred over quadratic form
- Post-Newtonian limit gives phi = psi (equal gravitational potentials), ensuring consistent lensing
- Strong-field limit recovers GR + Lambda automatically

### 7.3 AeST Hamiltonian Formalism (Bataki, Skordis, Zlosnik 2023)

**arXiv:2307.15126**

Full Hamiltonian analysis of AeST:
- **6 physical degrees of freedom** at the fully nonlinear level
- 4 first-class constraints (gauge symmetries)
- 4 second-class constraints (reduce phase space)
- Confirms internal consistency of the theory

For comparison:
- GR has 2 physical DoF (two gravitational wave polarizations)
- AeST has 6 = 2 (tensor) + 2 (vector, massive) + 2 (scalar, one massive + one non-propagating)
- The non-propagating scalar mode is the c_s = 0 mode that mimics CDM

### 7.4 Stealth Black Holes in AeST (Skordis & Vokrouhlicky 2024)

**arXiv:2412.15395**

Major finding for strong-field regime:
- AeST admits "stealth" black hole solutions with EXACT GR geometry (Schwarzschild or Reissner-Nordstrom metric)
- The aether and scalar carry non-trivial "secondary hair" despite the metric being identical to GR
- One class of solutions can be continuously joined to the cosmological FLRW solution
- Stability conditions: 0 < K_B < 2, K_2 > 0, lambda_s >= 0

**Implication for tau framework:** AeST black holes LOOK like GR black holes at the metric level. The tau framework predicts exponential metric (no event horizon). If AeST is the correct embedding of the tau framework, this creates a TENSION: AeST gives Schwarzschild, tau gives exponential. Resolution options:
1. The tau framework's exponential metric is wrong (AeST's stealth BH is correct)
2. The stealth BH is one solution; others may give exponential metric
3. The tau framework should be embedded in a modification of AeST, not AeST itself

### 7.5 AeST Quasi-Static Analysis (Mistele, McGaugh, Hossenfelder 2023)

**arXiv:2305.07742**

Key findings:
- New mass scale m_x = Q_0 sqrt((2-K_B)/(2K_B)) affects non-spherical systems
- m_x has no effect in spherical symmetry (doesn't change rotation curves)
- Two limiting behaviors:
  - m_x -> infinity: reduces to AQUAL (single-field MOND)
  - m_x -> 0: reduces to two-field MOND (original Bekenstein-Milgrom)
- Wide binary tests may constrain m_x at percent level

### 7.6 AeST Galaxy Cluster Problem (2023)

**arXiv:2312.00889**

Serious problem discovered:
- "Negative mass density" artifacts appear in galaxy cluster fits
- Parameter m^2/f_G constrained to be > 1 Mpc^{-2} by clusters
- But galaxy-scale fits prefer m^2/f_G < 1 Mpc^{-2}
- **Conflicting requirements from different scales**
- This mirrors MOND's general failure at cluster scales

### 7.7 AeST Weak Lensing Confrontation (2023)

**arXiv:2301.03499**

Results:
- Inner galaxy regions: AeST matches MOND (and observations) well
- Large galactocentric distances: deviations from MOND emerge
- Ghost condensate mass becomes significant at ~Mpc scales
- **Tension between galaxy and cluster parameter requirements**

### 7.8 Dynamical Systems Analysis of AeST Cosmology (2023)

**arXiv:2309.06232**

Cosmological fixed points identified:
- Radiation-dominated era (t^{1/2}): exists and is a saddle point
- Matter-dominated era (t^{2/3}): exists and attracts trajectories
- de Sitter (Lambda): exists as a late-time attractor
- Scalar-dominated (t^{1/3}): exists for certain F(Q)

**Critical finding:** "The dust-like contribution [from the scalar field] generally occurs for a limited period of cosmic time." This means the scalar field does NOT behave as perfect CDM for all time -- it only mimics CDM during the matter era. This is different from real CDM particles, which behave as dust forever.

**For the tau framework:** This is actually a TESTABLE PREDICTION. If AeST is correct, there should be subtle deviations from LCDM at very early times (before the scalar field enters its dust-like phase) and at very late times (when it exits). These deviations could potentially be measured with future CMB experiments (CMB-S4, LiteBIRD).

---

## 8. Honest Assessment: tau -> AeST

### 8.1 Is tau -> AeST a Derivation?

**No.** There is no calculation showing that D(rho_spacetime || rho_matter), when properly defined and varied, produces the AeST action. The steps required:

1. Define rho_spacetime and rho_matter as density matrices on a Hilbert space [NOT DONE]
2. Compute D(rho_spacetime || rho_matter) for a spacetime with preferred foliation [NOT DONE]
3. Show that the variation delta D = 0 gives Einstein + aether + scalar field equations [NOT DONE]
4. Identify the free function F(Y,Q) from the QRE structure [NOT DONE]
5. Show that perturbations have omega = 0 [NOT DONE]

**Assessment: 0/5 steps completed. This is pure conjecture, not derivation.**

### 8.2 Is tau -> AeST a Motivation?

**Yes, and a strong one.** The structural parallels are deep:

| Structural parallel | Strength |
|--------------------|----------|
| u^mu <-> A^mu <-> n^mu | **Very strong** (identical mathematical object) |
| u^mu u_mu = -1 <-> A^mu A_mu = -1 | **Exact match** |
| Sigma = -ln(-g_00) is observer-dependent <-> phi depends on A^mu frame | **Strong parallel** |
| DPI/RG flow <-> MOND interpolating function | **Suggestive but not proven** |
| Petz recovery as retrodiction <-> Khronon as time foliation | **Conceptual match** |
| QRE as action principle <-> GfE (Bianconi) | **Direct analog** |

The motivation is genuine: the tau framework's mathematical structures (observer, QRE, DPI) naturally point toward the AeST field content (aether, scalar, preferred foliation).

### 8.3 Is It Just Notation?

**More than notation, less than derivation.** Specifically:

**What is more than notation:**
- The identification u^mu = A^mu has physical content: it says "dark matter is the cost of time"
- The QRE structure constrains the functional form of F(Y,Q) (if derivable)
- The DPI provides a reason for MOND's interpolating function (running of the gravitational channel)

**What is just notation (currently):**
- Writing Sigma = -ln(-g_00) is equivalent to writing phi ~ Phi_N in AeST
- The tau = 1 - exp(-Sigma/2) is a one-to-one map of the same information
- The "Petz recovery fidelity = coordinate speed of light" is a restatement of the gravitational redshift

### 8.4 What Specific New Prediction Would the tau Connection Give?

**Prediction 1: Omega_DM h^2 from QRE extremization**

If the dark matter density is set by extremizing D(rho_spacetime || rho_matter) under cosmological boundary conditions, then Omega_DM h^2 = 0.12 should follow from the QRE structure, not be fitted. This would reduce AeST's parameters by one.

**Status:** No calculation exists. This is the single most valuable prediction the tau framework could provide.

**Prediction 2: a_0 from the cosmological QRE**

If the MOND acceleration scale a_0 emerges from the QRE at the cosmological horizon (a_0 ~ c H_0, as observed empirically), this would provide a first-principles explanation for the MOND-cosmology coincidence.

**Status:** The coincidence a_0 ~ c H_0 is well-known but unexplained. If the tau framework derives it, this would be significant.

**Prediction 3: The free function F(Y,Q) from DPI constraints**

If the MOND interpolating function mu(x) = dJ/dY is constrained by the data processing inequality of the gravitational channel's RG flow, then F(Y,Q) would be PREDICTED rather than chosen. This would eliminate AeST's most ad hoc element.

**Status:** The DPI connection to the MOND function is suggestive (Casini-Huerta RG = DPI) but no explicit derivation exists.

**Prediction 4: Deviations from AeST at the Planck scale**

The tau framework includes quantum corrections (Paper 1's bound F >= exp(-Sigma/2)) that AeST does not have. These could produce observable effects:
- Modified dispersion relation for the c_s = 0 mode: omega ~ (l_Planck k)^2 instead of omega = 0
- Minimum sound speed c_s,min ~ l_Planck / t_Hubble ~ 10^{-61}
- Deviations from perfect CDM behavior at very small scales

**Status:** Completely speculative. But this is a UNIQUE prediction that AeST alone cannot make.

### 8.5 The Key Gap

**The single most important calculation for Paper 4:**

Show that D(rho_spacetime || rho_matter), computed on a spacetime with a dynamical foliation (n^mu from a scalar field tau), and varied to give field equations, produces:
1. The Einstein equations (known: Dorau-Much 2025)
2. A scalar field equation with K(Q) kinetic structure
3. The MOND interpolating function J(Y) from the DPI/RG structure
4. Perturbations with omega = 0 (from the constrained QRE optimization)

If this calculation succeeds, it would:
- Derive AeST/Khronon from quantum information
- Predict F(Y,Q) rather than assuming it
- Potentially predict Omega_DM h^2
- Provide a physical interpretation for all fields

If it fails, the tau framework must accept AeST as an independent theory and acknowledge that the structural parallel is suggestive but not derivable.

### 8.6 Comparison with Existing Unification Attempts

| Approach | Derives from? | Predicts F? | Predicts Omega_DM? | Status |
|----------|-------------|------------|-------------------|--------|
| AeST (Skordis-Zlosnik) | Constructed to fit observations | No (chosen) | No (fitted) | Works |
| Khronon (Blanchet-Skordis) | Simplified AeST | No (chosen) | No (fitted) | Works |
| Bianconi GfE | QRE of metrics | No | No | Gives Schwarzschild, not AeST |
| Dorau-Much | QRE of states | N/A (gives GR only) | N/A | Incomplete |
| **tau framework (proposed)** | **QRE of spacetime + observer** | **Possibly (from DPI)** | **Possibly (from QRE extremum)** | **No calculation** |

### 8.7 The Stealth Black Hole Tension

The stealth BH result (arXiv:2412.15395) creates a specific tension:

- **AeST predicts:** Schwarzschild metric for BHs (stealth solutions)
- **tau framework predicts:** Exponential metric g_00 = -exp(-r_s/r) (no event horizon)

These are DIFFERENT and observationally distinguishable (Paper 2's testable predictions).

**Resolution options:**
1. tau framework embeds in a MODIFIED AeST where stealth solutions do not exist
2. The exponential metric is the correct weak-field approximation; stealth BH is the strong-field completion
3. The identification tau <-> AeST breaks down in the strong-field regime
4. Both theories are wrong in the strong-field regime (need quantum gravity)

**Honest assessment:** This tension is REAL and must be acknowledged. Paper 4 cannot claim tau derives AeST if their strong-field predictions disagree.

---

## 9. Key References

### AeST Theory (Primary)
1. **Skordis & Zlosnik 2021a**: arXiv:2007.00082, PRL 127, 161302 -- Original AeST paper (CMB + MOND)
2. **Skordis & Zlosnik 2021b**: arXiv:2109.13287 -- Extended version with stability analysis, free function example
3. **Skordis & Vokrouhlicky 2024**: arXiv:2412.15395 -- Stealth black holes in AeST

### AeST Analysis
4. **Bataki, Skordis & Zlosnik 2023**: arXiv:2307.15126 -- Hamiltonian formalism (6 DoF)
5. **Mistele, McGaugh & Hossenfelder 2023**: arXiv:2305.07742 -- Quasi-static limit, m_x scale
6. **arXiv:2312.00889** (2023) -- Galaxy cluster problems (negative mass density)
7. **arXiv:2301.03499** (2023) -- Weak lensing confrontation
8. **arXiv:2309.06232** (2023) -- Dynamical systems analysis of AeST cosmology

### Khronon Theory
9. **Blanchet & Skordis 2024**: arXiv:2404.06584, JCAP 11 (2024) 040 -- Khronon theory
10. **Blanchet & Skordis 2025**: arXiv:2507.00912 -- Khronon-Tensor extension

### nuHDM and Alternatives
11. **Haslbauer, Banik & Kroupa 2020**: arXiv:2009.11292 -- nuHDM
12. **Samaras, Grandis & Kroupa 2025**: arXiv:2506.19196 -- nuHDM initial conditions problem

### Foundational
13. **Kunz 2009**: arXiv:0702615 -- Generalized Dark Matter framework
14. **Hu 2008**: arXiv:0802.3688 -- CMB Theory pedagogical reference
15. **Planck 2018**: arXiv:1807.06209 -- Cosmological parameters

### Entropic/Information-Theoretic Gravity
16. **Dorau-Much 2025**: arXiv:2510.24491 -- QRE -> Einstein equations
17. **Bianconi 2025**: PRD 111, 066001; arXiv:2408.14391 -- Gravity from Entropy
18. **Bianconi 2025**: arXiv:2510.22545 -- GfE Cosmological solutions

### k-Essence and CDM from Scalars
19. **Scherrer 2004**: PRL 93, 011301 -- k-essence as unified dark matter
20. **Berezhiani & Khoury 2015**: arXiv:1507.01530 -- Superfluid dark matter

### Einstein-Aether Theory
21. **Jacobson & Mattingly 2001**: PRD 64, 024028 -- Einstein-aether theory (original)
22. **Jacobson 2008**: arXiv:0801.1547 -- Einstein-aether gravity (review)

---

## Appendix A: Glossary of AeST Notation

| Symbol | Definition | tau framework analog |
|--------|-----------|---------------------|
| A^mu | Unit timelike aether vector | Observer 4-velocity u^mu |
| phi (or phi_s) | Shift-symmetric scalar field | ~Sigma_grav / (-2) |
| Y | Spatial scalar kinetic term D^mu phi D_mu phi | ~(nabla Sigma)^2 / 4 |
| Q | Temporal scalar derivative A^mu nabla_mu phi | ~d(Sigma)/d(proper time) |
| F(Y,Q) | Free function in action | Should be derived from D(rho_st \|\| rho_m) |
| K_B | Aether kinetic coupling | No direct analog |
| lambda | Lagrange multiplier for A^2 = -1 | Observer constraint |
| J(Y) | MOND function (Khronon) | Should follow from DPI/RG of gravitational channel |
| K(Q) | Khronon kinetic function | Should follow from QRE perturbation structure |
| mu | Inverse-length Khronon mass | Related to cosmological QRE? |
| a_0 | MOND acceleration scale | ~ c H_0 (cosmological QRE horizon?) |
| Q_0 | Background Q value | Set by cosmological Sigma evolution? |
| F_Q = C/a^3 | Scalar "momentum" dilution | Dark matter energy density dilution |

---

## Appendix B: The c_s = 0 Mechanism in k-Essence Language

For readers familiar with k-essence, here is the AeST mechanism translated:

**Standard k-essence:** L = f(X) where X = -(1/2)(nabla phi)^2

The background has X = X_0 = (1/2)(phi_dot)^2. Sound speed:
```
c_s^2 = f_X / (f_X + 2 X f_{XX})
```

For c_s = 0: need f_X(X_0) = 0. This gives a degenerate kinetic structure.

**Problem:** If f_X = 0, the kinetic energy vanishes and the field has no dynamics. Perturbations grow exponentially (Jeans instability) but the background solution is unstable.

**AeST solution:** Separate X into temporal (Q^2/2) and spatial (Y/2) parts using the aether:
```
X = Q^2/2 - Y/2     (schematically)
```

Now can have:
- F_Q ≠ 0 (temporal kinetics non-degenerate -> background dynamics healthy)
- F_Y ~ 0 or small (spatial kinetics suppressed -> c_s ~ 0)

The aether enables this split. Without the aether, X is a single variable and you cannot independently control temporal and spatial propagation.

**In Khronon language:** K(Q) = mu^2(Q-1)^2 gives:
- K_Q = 2 mu^2 (Q-1) ≠ 0 around Q ~ 1 + epsilon
- But the PERTURBATION dispersion is omega = 0 because the mode sits at the potential minimum

The non-propagating nature comes from the CONSTRAINT-like structure: Q is dynamically driven toward Q ~ 1 by the stiff potential, and perturbations around this minimum are "frozen."

---

## Appendix C: The Six Degrees of Freedom of AeST

From Bataki, Skordis & Zlosnik (2023, arXiv:2307.15126):

**Total DoF at the nonlinear level: 6**

| Mode | Count | Nature | Dispersion | Mass | Role |
|------|-------|--------|------------|------|------|
| Tensor (GW) | 2 | Propagating, massless | omega = c k | 0 | Standard gravitational waves |
| Vector | 2 | Propagating, massive | omega^2 = c^2 k^2 + m_V^2 | m_V | Constrained by GW observations |
| Scalar (massive) | 1 | Propagating, massive | omega^2 = c_s^2 k^2 + m_S^2 | m_S | Drives MOND phenomenology |
| Scalar (non-propagating) | 1 | Non-propagating | omega = 0 | - | **CDM-like mode** |

**The CDM-like mode is the 6th degree of freedom.** It is the one with omega = 0 that grows under gravitational instability without oscillating. This is the mode that mimics CDM.

The other 5 DoF (2 tensor + 2 vector + 1 massive scalar) are propagating modes with standard dispersion relations. They do not contribute to CDM-like behavior because they propagate and do not grow in place.

---

## Appendix D: Timeline of MOND Relativistic Theories

| Year | Theory | Authors | CMB? | c_T = c? | Clusters? |
|------|--------|---------|------|----------|-----------|
| 1984 | AQUAL | Bekenstein & Milgrom | No (non-relativistic) | N/A | No |
| 2004 | TeVeS | Bekenstein | Partial | **No** (ruled out by GW170817) | No |
| 2009 | Bimetric MOND | Milgrom | No | Unknown | No |
| 2010 | Generalized TeVeS | Skordis | Improved | No | No |
| 2015 | Superfluid DM | Berezhiani & Khoury | Yes (uses CDM) | Yes | Yes (uses CDM) |
| 2020 | nuHDM | Haslbauer et al. | Partial (misses 4th peak) | N/A | Partial |
| **2021** | **AeST** | **Skordis & Zlosnik** | **Yes** | **Yes** | **Problematic** |
| **2024** | **Khronon** | **Blanchet & Skordis** | **Claims yes** | **Yes** | **Unknown** |
| 2025 | Khronon-Tensor | Blanchet & Skordis | Claims yes | Yes | Unknown |

**AeST and Khronon are the ONLY theories that satisfy all three of: (1) CMB fit, (2) c_T = c, (3) MOND at galaxy scales.** Their shared failure point is galaxy clusters.

---

*Last updated: 2026-03-12*
*Research conducted for Paper 4 of the tau framework series: Sigma = D(rho_spacetime || rho_matter)*
*This document provides the technical foundation for understanding how AeST achieves CDM-like behavior and what the tau framework would need to reproduce it.*
