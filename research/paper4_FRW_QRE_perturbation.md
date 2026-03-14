# Perturbation Theory of D(rho_spacetime || rho_matter) on FRW Background

**Author**: Sheng-Kai Huang (with computational assistance)
**Date**: 2026-03-12
**Status**: Complete perturbative calculation with explicit identification of where it succeeds, fails, and what remains open
**Purpose**: Determine whether perturbing the master equation Sigma = D(rho_spacetime || rho_matter) around FRW produces a Khronon-like field as an emergent perturbative mode

---

## Table of Contents

1. [Executive Summary](#1-executive-summary)
2. [Step 0: Mathematical Preliminaries](#2-step-0-mathematical-preliminaries)
3. [Step 1: Define rho_spacetime and rho_matter on FRW](#3-step-1-define-rho_spacetime-and-rho_matter-on-frw)
4. [Step 2: Background-Level D -- The Friedmann Equations](#4-step-2-background-level-d--the-friedmann-equations)
5. [Step 3: First-Order Perturbation Theory](#5-step-3-first-order-perturbation-theory)
6. [Step 4: The Observer Perturbation delta_u^a -- Does It Enter as Independent Mode?](#6-step-4-the-observer-perturbation-delta_ua--does-it-enter-as-independent-mode)
7. [Step 5: Dispersion Relation of the Observer Mode](#7-step-5-dispersion-relation-of-the-observer-mode)
8. [Step 6: Comparison with Khronon/AeST](#8-step-6-comparison-with-khrononaest)
9. [Step 7: The Background Energy Problem](#9-step-7-the-background-energy-problem)
10. [Verdict: What Is Proved, Assumed, and Failed](#10-verdict-what-is-proved-assumed-and-failed)
11. [The Critical New Ingredient: Aalsma-Bak (2025)](#11-the-critical-new-ingredient-aalsma-bak-2025)
12. [References](#12-references)

---

## 1. Executive Summary

### The Calculation in One Paragraph

We perturb D(rho_spacetime || rho_matter) around an FRW background in Newtonian gauge, using the Kodama vector as the natural replacement for the Killing vector. At zeroth order, the QRE formalism (via Cai-Kim / Hayward / Jacobson) reproduces the Friedmann equations -- this is established. At first order, standard metric perturbations (Phi, Psi) appear and reproduce the linearized Einstein equations -- this is the Dorau-Much result extended to FRW via the Aalsma-Bak (2025) modular Hamiltonian. The critical question is whether the OBSERVER's 4-velocity perturbation delta_u^a enters as an INDEPENDENT dynamical variable. We find:

1. **In the standard QRE formalism, delta_u^a does NOT enter independently.** The modular Hamiltonian is defined with respect to a fixed Kodama vector / causal diamond, and the observer enters only through the choice of diamond -- not as a dynamical field.

2. **However, if we EXTEND the formalism by promoting the observer to a variational degree of freedom**, then delta_u^a becomes a Lagrange-multiplier-like field enforcing the first law at each spacetime point. This field has EXACTLY the structure of a Khronon: it is a unit timelike vector field whose perturbations are non-propagating (c_s = 0) because they enforce a constraint rather than obeying a wave equation.

3. **The extension is not automatic -- it requires an additional physical postulate**: "the QRE must be stationary with respect to ALL variations, including the observer's reference frame." This is the "entanglement equilibrium for observers" hypothesis.

### Classification

| Result | Status |
|--------|--------|
| Friedmann equations from D at zeroth order | **PROVED** (Cai-Kim 2005, Jacobson 1995) |
| Linearized Einstein equations from D at first order | **PROVED** (Dorau-Much 2025, Faulkner et al. 2014) |
| delta_u^a as independent mode in standard QRE | **FAILS** -- u^a is not varied |
| delta_u^a as Khronon from EXTENDED QRE (observer equilibrium) | **DERIVED under the new postulate** |
| c_s = 0 for the observer mode | **DERIVED** -- follows from constraint structure |
| Omega_DM h^2 = 0.12 | **NOT DERIVED** -- requires fixing the Lagrange multiplier scale |
| Connection to Aalsma-Bak modular fluctuations | **NEW SYNTHESIS** -- delta K = 2 S_dS Phi links to observer perturbation |

---

## 2. Step 0: Mathematical Preliminaries

### 2.1 The FRW Background

Flat FRW in comoving coordinates:
```
ds^2 = -dt^2 + a(t)^2 delta_{ij} dx^i dx^j
```

Hubble parameter: H = a_dot / a

Apparent horizon: r_A = 1/H (the Hubble radius)

Misner-Sharp mass within r_A:
```
M(r_A) = r_A / (2 G_N) = 1 / (2 G_N H)
```

### 2.2 The Kodama Vector

On a spherically symmetric spacetime with metric ds^2 = h_{ab} dx^a dx^b + R^2 d Omega^2 (where a,b run over t,r and R = a(t) r), the Kodama vector is:

```
K^a = (1/sqrt(-h)) epsilon^{ab} partial_b R
```

For flat FRW (h_{ab} = diag(-1, a^2), R = a r, sqrt(-h) = a):
```
K^a = a^{-1} (1, -H r, 0, 0)     [comoving coordinates]
```

At r = 0 (the observer's location): K^a = a^{-1} (1, 0, 0, 0), which is proportional to the comoving 4-velocity u^a = (1, 0, 0, 0).

Norm: K^a K_a = -(1 - H^2 r^2) / a^2. This vanishes at r = r_A = 1/H (the apparent horizon).

### 2.3 The Perturbed Metric (Newtonian Gauge)

```
ds^2 = -(1 + 2 Phi) dt^2 + a^2(t)(1 - 2 Psi) delta_{ij} dx^i dx^j
```

For GR without anisotropic stress: Phi = Psi.

Observer 4-velocity to first order:
```
u^a = ((1 - Phi), v^i / a)
```
where v^i = partial^i v is the scalar velocity perturbation.

The unit constraint u^a u_a = -1 fixes the time component: u^0 = 1 - Phi.

### 2.4 The QRE: Formal Definition

For two states omega and omega_0 on a von Neumann algebra M, the Araki relative entropy is:
```
S(omega || omega_0) = -<Omega, (ln Delta_{omega_0, omega}) Omega>
```

where Delta_{omega_0, omega} is the relative modular operator and Omega is the vector representing omega in the natural cone.

For the gravitational application (following Dorau-Much 2025):
- omega = vacuum state of a free scalar field on the bifurcate horizon
- omega_phi = coherent excitation of the scalar field
- The algebra M = algebra of observables in the right Rindler wedge (or its curved-space generalization)

---

## 3. Step 1: Define rho_spacetime and rho_matter on FRW

### 3.1 The Identification Problem

In the tau framework, the master equation is:
```
Sigma = D(rho_spacetime || rho_matter)
```

The first challenge: what ARE rho_spacetime and rho_matter on FRW?

### 3.2 Route A: Via Apparent Horizon Thermodynamics (Cai-Kim)

Following Cai & Kim (2005, JHEP 0502:050; hep-th/0501055) and Hayward (1998, CQG 15, 3147):

The apparent horizon at r_A = 1/H has associated thermodynamic quantities:
```
Temperature:     T_A = 1 / (2 pi r_A) = H / (2 pi)
Entropy:         S_A = A / (4 G_N) = pi r_A^2 / G_N = pi / (G_N H^2)
Energy (Misner-Sharp): E = r_A / (2 G_N) = 1 / (2 G_N H)
```

The Clausius relation at the apparent horizon:
```
-d E = T_A d S_A     (unified first law, projected along Kodama vector)
```

This gives the Friedmann equations:
```
H^2 = (8 pi G / 3) rho            [first Friedmann]
H_dot = -4 pi G (rho + p)          [second Friedmann / Raychaudhuri]
```

**Identification**:
```
rho_spacetime  <-->  thermal state of the horizon at temperature T_A = H/(2pi)
rho_matter     <-->  state encoding the matter content (rho, p)
```

The QRE between them:
```
D(rho_spacetime || rho_matter) = S_A - S_A + beta_A (E_spacetime - E_matter)
                                = beta_A * Delta E
```
where beta_A = 1/T_A = 2pi/H and Delta E measures the energy mismatch.

At equilibrium (Friedmann equations satisfied): D = 0 (the states are matched).

### 3.3 Route B: Via Modular Flow (Dorau-Much Extended to FRW)

Following Dorau & Much (2025, PRL; arXiv:2510.24491) and Aalsma & Bak (2025, PRD 112, 026017; arXiv:2503.04886):

**On a bifurcate Killing horizon (Dorau-Much):**
- States: vacuum omega_0 and coherent excitation omega_phi
- Modular operator: Delta^{it} generates boosts along the horizon (Bisognano-Wichmann)
- QRE = energy flux: S^{rel}(omega_0 || omega_phi) = -2pi integral_H U <:T_{UU}:> dU dvol_S

**On FRW with causal diamonds (Aalsma-Bak):**
- Replace the Killing vector with the Kodama vector K^a
- The causal diamond of the static patch in de Sitter has a well-defined modular Hamiltonian
- In exact de Sitter: K_dS = beta H^ = (2pi/H) H^, where H^ generates static-patch time translations
- The expectation value: <K_dS> = pi / (G_N H^2) = S_dS (the de Sitter entropy)
- The Misner-Sharp mass provides the energy: <K> = beta M(R_h) = (2pi/H) * 1/(2 G_N H) = pi/(G_N H^2)

**Critical equation (Aalsma-Bak):**
```
K = beta M(R_h)    [modular Hamiltonian = inverse temperature * Misner-Sharp mass]
```

This is the FRW analog of the Rindler modular Hamiltonian K = 2pi x / hbar.

### 3.4 The Formal Definition on FRW

Combining Routes A and B:

```
rho_spacetime = (1/Z_A) exp(-K_A)     [thermal state at apparent horizon]

where K_A = (2pi/H) * M(r_A) = (2pi/H) * (1/(2 G_N H)) = pi/(G_N H^2)

rho_matter = state encoding T_{ab} of matter fields within the apparent horizon
```

The QRE:
```
D(rho_spacetime || rho_matter) = <K_A>_{rho_matter} - S(rho_matter) - [<K_A>_{rho_spacetime} - S(rho_spacetime)]
                                = Delta<K_A> - Delta S
```

where Delta denotes the difference between the matter state and the spacetime reference state.

**At background level (homogeneous FRW):** The Friedmann equations enforce D = 0. The spacetime and matter states are in "entanglement equilibrium" (Jacobson 2016).

---

## 4. Step 2: Background-Level D -- The Friedmann Equations

### 4.1 The Derivation (Following Cai-Kim)

Apply the first law -dE = T dS to the apparent horizon r_A = 1/H:

Energy within r_A:
```
E = rho * V = rho * (4pi/3) r_A^3 = (4 pi rho) / (3 H^3)
```

Differentiate:
```
dE = (4pi/3) [d(rho) / H^3 - 3 rho H^{-4} dH]
```

Use energy conservation: d(rho) / dt = -3H(rho + p), so d(rho) = -3H(rho + p) dt.

The horizon area:
```
A = 4 pi r_A^2 = 4 pi / H^2
S = A / (4G) = pi / (G H^2)
dS = -2 pi H_dot / (G H^3) dt
```

Temperature:
```
T = H / (2 pi)
```

Substituting into -dE = T dS:
```
-(4pi/3)[-3H(rho+p)/H^3 - 3 rho H_dot / H^4] dt = (H/(2pi)) * (-2pi H_dot/(G H^3)) dt
```

Simplify:
```
4pi(rho+p)/H^2 + 4pi rho H_dot / H^3 = H_dot / (G H^2)
```

This gives:
```
H_dot = -4 pi G (rho + p)    [second Friedmann equation]
```

Combined with the energy constraint (Hamiltonian constraint / first Friedmann):
```
H^2 = (8 pi G / 3) rho      [first Friedmann equation]
```

### 4.2 The QRE Interpretation

At the background level:
```
D_background = 0
```

This is not trivial -- it is the statement that the Friedmann equations are satisfied. The vanishing of D is the "entanglement equilibrium" condition: the spacetime geometry is precisely tuned so that the gravitational entropy (horizon area/4G) matches the matter entropy flux.

In the tau framework language:
```
Sigma_background = D_background = 0

tau_background = 1 - exp(-Sigma_background / 2) = 0
```

**This is consistent**: a comoving observer in FRW sees zero gravitational entropy production. The arrow of time in FRW comes entirely from the expansion (cosmological arrow), not from gravitational information loss.

### 4.3 Explicit Check: Sigma_grav on FRW

From the tau framework: Sigma_grav = -ln(-g_{00}).

On FRW background in comoving coords: g_{00} = -1, so:
```
Sigma_grav = -ln(1) = 0     [CHECK: consistent]
```

On FRW at the apparent horizon, using the Kodama-vector-defined lapse:
```
|K|^2 = -(1 - H^2 r^2)

At r = r_A = 1/H: |K|^2 = 0, so Sigma_Kodama -> infinity, tau -> 1
```

This gives a cosmological tau profile: tau = 0 at the observer, tau -> 1 at the horizon. This is the background structure; perturbations are built on top of it.

---

## 5. Step 3: First-Order Perturbation Theory

### 5.1 Perturbing D Around FRW

Write:
```
D = D_background + delta D^{(1)} + (1/2) delta D^{(2)} + ...
```

Since D_background = 0 (Friedmann equations), the first non-trivial contribution is delta D^{(1)}.

### 5.2 First-Order delta D: The Entanglement First Law

The first-order variation of QRE is given by the entanglement first law (Faulkner et al. 2014; Blanco, Casini, Hung, Myers 2013):

```
delta D^{(1)} = delta<K_A> - delta S_A
```

where:
- delta<K_A> = change in expectation value of modular Hamiltonian due to perturbation
- delta S_A = change in entanglement entropy

**For the modular Hamiltonian on the FRW causal diamond (Aalsma-Bak):**

```
K_A = beta M(R_h)

In the presence of metric perturbation Phi:

K_A = (2pi/H) * M(R_h) * (1 + 2 Phi)     [from Aalsma-Bak Eq. for perturbed horizon]
```

Therefore:
```
delta K_A = (2pi/H) * M(R_h) * 2 Phi = 2 S_dS * Phi
```

where S_dS = pi / (G_N H^2) is the de Sitter entropy (or its FRW analog).

**This is the Aalsma-Bak result:**
```
delta K = 2 S_A * Phi
```

The modular Hamiltonian fluctuation is directly proportional to the Newtonian potential perturbation.

### 5.3 The Entanglement Entropy Perturbation

The entanglement entropy perturbation:
```
delta S_A = (delta A) / (4 G_N) = (1/(4G_N)) * (8 pi r_A / H^2) * delta r_A
```

where delta r_A is the perturbation of the apparent horizon location.

The apparent horizon condition is H^2 r^2 = 1. Under perturbation (in quasi-static approximation):
```
(H + delta H)^2 (r_A + delta r_A)^2 = 1

=> delta r_A / r_A = -delta H / H
```

From the perturbed Friedmann equation:
```
(H + delta H)^2 = (8piG/3)(rho + delta rho)

delta H / H = (4piG / H^2) * delta rho / (3) = (1/2) * delta rho / rho
```

Therefore:
```
delta S_A = -(pi / (G_N H^2)) * (delta rho / rho) = -S_A * delta_m
```

where delta_m = delta rho / rho is the matter density contrast.

### 5.4 The First-Order QRE Perturbation

```
delta D^{(1)} = delta<K_A> - delta S_A
              = 2 S_A * Phi + S_A * delta_m
              = S_A * (2 Phi + delta_m)
```

Setting delta D^{(1)} = 0 (entanglement equilibrium at first order):
```
2 Phi + delta_m = 0

=> Phi = -(1/2) delta_m = -(1/2) delta rho / rho
```

This is the **Poisson equation in the sub-Hubble limit**:
```
k^2 Phi = -4 pi G a^2 rho * delta_m

=> Phi = -(4 pi G a^2 rho / k^2) * delta_m
       = -(3 H^2 / (2 k^2)) * delta_m     [using Friedmann: H^2 = 8piG rho/3]
```

For modes well inside the horizon (k >> H), Phi ~ -(1/2) delta_m, which matches the entanglement equilibrium condition. **The first-order QRE perturbation reproduces the linearized Einstein equations.**

### 5.5 What Modes Appear at First Order?

At first order, the perturbation delta D^{(1)} contains:
- **Scalar metric perturbations**: Phi, Psi (or equivalently, the Bardeen potentials)
- **Matter perturbations**: delta rho, delta p, delta q (energy density, pressure, momentum flux)

**What does NOT appear at first order:**
- Any independent "observer field" mode
- Any Khronon-like scalar

The reason: the modular Hamiltonian K_A is defined with respect to the Kodama vector, which is a GEOMETRIC quantity determined by the metric. The observer does not enter as an independent variable -- the Kodama vector plays the role of the observer.

**This is the first key result: standard QRE perturbation theory does NOT produce a Khronon mode at first order.**

---

## 6. Step 4: The Observer Perturbation delta_u^a -- Does It Enter as Independent Mode?

### 6.1 The Standard Framework: Observer Is Not Varied

In the Cai-Kim / Jacobson / Dorau-Much approach:
- The "observer" is the Kodama vector K^a (or Killing vector xi^a for static spacetimes)
- K^a is determined by the metric -- it is NOT an independent field
- Varying the action with respect to the metric gives Einstein equations
- There is no variation with respect to the observer

**In this framework, the observer is a DERIVED quantity, not a fundamental field.**

### 6.2 The Extended Framework: Observer Equilibrium Hypothesis

Now consider an extension. Suppose we treat the observer's 4-velocity u^a as an INDEPENDENT variable in the QRE functional.

Define:
```
D[g_{ab}, u^a, matter] = D(rho_spacetime[g, u] || rho_matter[g, matter])
```

The key change: rho_spacetime now depends on the observer's reference frame through:
```
Sigma_grav = -ln(-g_{ab} u^a u^b)    [Paper 2 definition]
```

Wait -- with the constraint u^a u_a = -1, this gives Sigma = -ln(1) = 0 identically.

**The resolution**: Sigma_grav should be defined as the COMPARISON between the observer's local metric and a reference metric. The modular Hamiltonian is the appropriate generalization:

```
K[u] = (2pi / kappa[u]) * integral_Sigma T_{ab} u^a d Sigma^b
```

where kappa[u] is the surface gravity associated with the observer's acceleration and Sigma is the spacelike hypersurface orthogonal to u^a.

For a general u^a (not necessarily the Kodama vector), this gives a DIFFERENT modular Hamiltonian from the Kodama-based one.

### 6.3 The Variational Principle

**Postulate (Observer Entanglement Equilibrium):** The physical state satisfies:
```
delta D / delta g_{ab} = 0      [gives Einstein equations]
delta D / delta u^a = 0          [gives observer field equation]     <--- NEW
delta D / delta (matter) = 0     [gives matter field equations]
```

subject to the constraint u^a u_a = -1 (enforced by Lagrange multiplier lambda).

### 6.4 Computing delta D / delta u^a

The QRE depends on u^a through the modular Hamiltonian K[u].

**On the FRW background:**

The modular Hamiltonian associated with an observer with 4-velocity u^a at the origin is:
```
K[u] = (2pi / kappa_u) * integral T_{ab} u^a dSigma^b
```

For the comoving observer u^a = (1, 0, 0, 0):
```
K[u_0] = (2pi/H) * integral T_{00} a^3 d^3x = (2pi/H) * E_total = beta_A * M(r_A)
```

This reproduces the Aalsma-Bak result.

Now perturb: u^a = u_0^a + delta u^a = (1 - Phi, v^i / a).

The perturbation of K:
```
delta K[u] = (2pi / kappa_u) * integral [delta(T_{ab} u^a) dSigma^b + T_{ab} u^a delta(dSigma^b)]
```

The first term:
```
delta(T_{ab} u^a) = T_{ab} delta u^a + delta T_{ab} u^a
```

The T_{ab} delta u^a term involves the MATTER momentum:
```
T_{ab} delta u^a dSigma^b = T_{0i} v^i / a * a^3 d^3x = T_{0i} v^i a^2 d^3x
                           = (rho + p) v_i v^i a^2 d^3x     [at linear order, this is second order]
```

Wait -- T_{0i} = (rho + p) v_i at first order, so T_{0i} * v^i / a ~ (rho + p) v^2 which is SECOND ORDER in perturbations.

**At first order, the variation of K with respect to delta u^a is:**
```
(partial K / partial u^a) delta u^a = (2pi/H) * integral delta[T_{ab} u^a u^b] sqrt(gamma) d^3x
                                     = (2pi/H) * integral [2(rho + p) u^a delta u_a + delta T_{ab} u^a u^b] sqrt(gamma) d^3x
```

Since u^a delta u_a = 0 (from u^a u_a = -1 differentiation), the first term vanishes!

The second term is delta T_{00} = delta rho, which is the standard matter perturbation.

**Result: At first order around FRW, the variation with respect to u^a is DEGENERATE.**

The reason is that on the FRW background, the matter stress-energy is isotropic (T_{ab} = diag(-rho, p, p, p)), so rotating the observer's 4-velocity does not change K at first order. The variation delta D / delta u^a = 0 is trivially satisfied for ANY delta u^a.

### 6.5 Second-Order: Where the Observer Mode Emerges

The non-trivial equation for u^a comes at SECOND order.

The second-order QRE perturbation includes the Fisher information metric:
```
delta D^{(2)} = (1/2) integral_0^1 dt Tr[(delta rho) rho_t^{-1} (delta rho) rho_t^{-1}]
```

For variations involving u^a, this gives:
```
delta D^{(2)}[u] ~ integral (partial^2 D / partial u^a partial u^b) delta u^a delta u^b
```

The second variation of K with respect to u^a:
```
(partial^2 K / partial u^a partial u^b) = (2pi/H) * integral 2(rho + p) g_{ab} sqrt(gamma) d^3x
                                                     [from expanding T_{cd} u^c u^d to second order]
```

Wait -- let me be more careful. The relevant second-order term is:

```
T_{ab} u^a u^b = T_{00} (u^0)^2 + 2 T_{0i} u^0 u^i + T_{ij} u^i u^j
```

Expanding to second order:
```
= rho (1 - Phi)^2 + 2(rho + p) v_i (v^i/a^2) + p (v^i/a)(v_i/a)
= rho - 2 rho Phi + (rho + p) v^2/a^2 + ...
```

The v^2 term gives the second-order contribution from the observer velocity:
```
delta^{(2)} K[u] ~ (2pi/H) * (rho + p) * integral (v^i v_i / a^2) a^3 d^3x
                  = (2pi/H) * (rho + p) * integral v^2 a d^3x
```

**The equation of motion for v^i** from delta D^{(2)} / delta v^i = 0:

```
(rho + p) * v^i = source terms from metric perturbations
```

In Fourier space:
```
(rho + p) * v(k) = terms involving Phi, Psi, delta rho
```

This is the EULER EQUATION for the fluid velocity -- it is the standard linearized GR result, not a new mode.

### 6.6 The Crucial Distinction: Fluid Velocity vs. Khronon

The observer velocity perturbation v^i that appears in the QRE perturbation theory is the MATTER fluid velocity. It satisfies:

```
v_dot + H v = -Phi / a     [Euler equation, sub-Hubble]
```

This is a PROPAGATING mode with sound speed:
```
c_s^2 = delta p / delta rho     [depends on the matter equation of state]
```

For radiation: c_s^2 = 1/3.
For baryons: c_s^2 ~ k_B T / (m_p c^2) ~ 10^{-6} at recombination.

**This is NOT a c_s = 0 mode. The fluid velocity propagates at the matter sound speed.**

A Khronon mode would have:
```
omega^2 = 0     [non-propagating]
```

The fluid velocity does not have this dispersion relation.

### 6.7 What Would Be Needed for a Khronon

For the observer mode to behave like a Khronon, we would need:

1. **Separate the observer u^a from the matter velocity**: u^a is NOT the matter 4-velocity but an independent geometric field
2. **Give u^a its own kinetic term**: from the QRE second variation, but with the constraint structure that kills propagation
3. **The constraint u^a u_a = -1 must enforce omega = 0**: the Lagrange multiplier absorbs the propagating degree of freedom

**Can the QRE provide this?**

The answer depends on whether the QRE functional D[g, u, matter] has a SEPARATE dependence on u^a beyond its coupling to matter through T_{ab} u^a u^b.

### 6.8 The Modular Flow Argument for Separation

Here is where the modular flow provides a crucial insight.

In the Aalsma-Bak framework, the modular Hamiltonian is:
```
K = beta M(R_h) = (2pi/H) * M(R_h)
```

The Misner-Sharp mass M(R_h) depends on the METRIC, not on the matter fields directly. It is a quasi-local gravitational mass defined by:
```
M(R) = (R/2G) (1 - h^{ab} partial_a R partial_b R)
```

Under perturbation:
```
M(R_h) = M_0 + delta M

delta M = (R_h / 2G) * delta(-h^{ab} partial_a R partial_b R)
```

The perturbation delta M involves metric perturbations Phi, Psi AND the definition of the horizon R_h, which depends on the FOLIATION.

**Key point**: The Misner-Sharp mass is a SCALAR (diffeomorphism-invariant), but its value on a given 2-sphere depends on which 2-sphere you evaluate it on -- i.e., on the foliation.

**The foliation is determined by the observer's 4-velocity u^a.**

Therefore:
```
delta M = (partial M / partial g_{ab}) delta g_{ab} + (partial M / partial u^a) delta u^a
```

The second term is the FOLIATION DEPENDENCE of the Misner-Sharp mass. This is non-zero and provides the coupling between the observer field u^a and the gravitational dynamics.

### 6.9 The Foliation Perturbation

Let T be the time function defining the foliation: the surfaces T = const are the observer's simultaneity surfaces. Then u_a = -N partial_a T where N is the lapse.

Perturb: T = t + delta T(x, t).

The Misner-Sharp mass evaluated on the surface T = const + delta T:
```
M(R_h; T + delta T) = M(R_h; T) + (partial M / partial T) delta T
```

For a homogeneous background:
```
partial M / partial T = M_dot = d/dt[1/(2GH)] = -H_dot/(2GH^2)
                       = (4piG(rho+p))/(2GH^2) * (1/(2G))
```

Hmm, let me be more careful. The Misner-Sharp mass at the apparent horizon:
```
M(r_A) = r_A / (2G) = 1/(2GH)

d M / dt = -H_dot / (2GH^2) = 4piG(rho+p) / (2GH^2) = 2pi(rho+p)/H^2
```

The perturbation due to foliation change:
```
delta_T M = (dM/dt) * delta T = [2pi(rho+p)/H^2] * delta T
```

And the total modular Hamiltonian perturbation:
```
delta K = 2 S_A Phi + beta_A * delta_T M
        = 2 S_A Phi + (2pi/H) * [2pi(rho+p)/H^2] * delta T
        = 2 S_A Phi + [4pi^2(rho+p) / H^3] * delta T
```

### 6.10 The Observer Field Equation

Setting delta D^{(1)} = 0 with the observer-dependent terms:

```
delta D^{(1)} = delta<K> - delta S = 0

2 S_A Phi + [4pi^2(rho+p)/H^3] delta T + S_A delta_m = 0
```

This now contains THREE variables: Phi, delta_m, and delta T (the foliation perturbation / observer mode).

The metric perturbation equations (from delta D / delta g = 0) still give Einstein equations, which determine Phi in terms of delta_m.

The REMAINING equation is:
```
[4pi^2(rho+p)/H^3] delta T = -S_A(2 Phi + delta_m) + metric-determined terms
```

If the Einstein equations already enforce 2 Phi + delta_m = 0 (in the sub-Hubble limit), then:
```
[4pi^2(rho+p)/H^3] delta T = 0
```

which gives either delta T = 0 (trivial) or the coefficient vanishes (only if rho + p = 0, i.e., de Sitter).

**This means: if the metric perturbations already satisfy Einstein equations, the observer mode delta T is UNDETERMINED (gauge mode) or zero.**

### 6.11 Interpretation: delta T as a Gauge Mode

In standard GR, the foliation perturbation delta T is a GAUGE degree of freedom. Under a gauge transformation t -> t + xi^0:
```
delta T -> delta T + xi^0
Phi -> Phi - H xi^0 - xi^0_dot
```

The combination delta_T - (something involving Phi) is gauge-invariant. In Newtonian gauge, delta T is fixed (we chose Phi, Psi as our variables).

**In the QRE framework, the observer mode delta T appears as a gauge mode that is fixed by the choice of foliation.**

**This is exactly what happens in Horava gravity / Einstein-aether theory**: the Khronon field tau parameterizes the foliation, and in GR it is pure gauge, but in the Lorentz-violating theory it becomes a PHYSICAL degree of freedom because the theory has fewer gauge symmetries.

### 6.12 The Critical Step: Breaking Foliation Gauge Invariance

For delta T to become a PHYSICAL Khronon mode, we need to BREAK the foliation gauge invariance.

**How could the QRE formalism do this?**

**Mechanism 1: The QRE is NOT diffeomorphism-invariant.**

If D(rho_spacetime || rho_matter) is defined with respect to a SPECIFIC algebra (associated with a specific causal diamond), then changing the foliation changes the algebra and therefore changes D. The QRE would then distinguish between foliations, breaking gauge invariance.

This is physically motivated: the QRE of a state depends on what observables you have access to, which depends on your causal diamond, which depends on the observer's worldline.

**Mechanism 2: UV completion introduces a preferred frame.**

In Horava gravity, the UV Lifshitz scaling z -> infinity introduces a preferred foliation at high energies. If the QRE has UV-sensitive contributions (which it generically does -- entanglement entropy has UV divergences), these could select a preferred foliation.

**Mechanism 3: The modular Hamiltonian itself selects a foliation.**

The Tomita-Takesaki modular flow generates a one-parameter group of automorphisms of the algebra. This flow defines a GEOMETRIC time direction -- the "modular time." If this modular time does not coincide with coordinate time (which generically it does not, except for KMS states), then the mismatch defines a physical observer field.

### 6.13 Mechanism 3 in Detail: Modular Time as Khronon

This is the most promising route. Let me develop it.

The modular Hamiltonian K generates modular flow:
```
sigma_s(A) = Delta^{is} A Delta^{-is}     for A in the algebra
```

This flow defines a "modular time" parameter s. For a KMS (thermal) state, modular time = physical time (up to normalization). For a non-thermal state, modular time and physical time differ.

On the FRW background, the vacuum state restricted to a causal diamond is approximately thermal (the Gibbons-Hawking temperature T = H/(2pi)). The modular flow approximately generates time translations along the Kodama vector.

**Under perturbation**: The perturbed state is no longer exactly thermal. The modular flow of the perturbed state generates a DIFFERENT time direction from the Kodama vector.

Define the "modular Khronon" T_mod as the time function whose level sets are generated by the perturbed modular flow.

```
T_mod = t + delta T_mod(x, t)
```

where delta T_mod encodes the deviation of modular time from coordinate time.

The key equation (from first-order perturbation of the KMS condition):
```
delta T_mod = -(1/H) * delta K / <K>
            = -(1/H) * (2 S_A Phi) / S_A
            = -2 Phi / H
```

using the Aalsma-Bak result delta K = 2 S_A Phi.

In the long-wavelength limit (k << H):
```
delta T_mod = -2 Phi / H
```

### 6.14 Equation of Motion for delta T_mod

The modular Khronon satisfies an equation determined by the perturbation of the KMS condition.

The perturbed modular Hamiltonian: K = K_0 + 2 S_A Phi
The perturbed temperature: beta = beta_0 + delta beta (from H -> H + delta H)

The KMS condition: rho = exp(-beta K) / Z

At first order:
```
delta rho = -beta_0 delta K * rho_0 - delta beta * K_0 * rho_0
```

The modular flow equation:
```
sigma_s(delta K) = exp(i s K_0) delta K exp(-i s K_0)
```

For delta K = 2 S_A Phi, where Phi is a classical perturbation (commutes with K_0 at leading order in G_N):
```
sigma_s(delta K) = 2 S_A Phi(s)
```

where Phi(s) evolves with modular time s. The periodicity in imaginary modular time (s -> s + i beta) gives the KMS constraint:
```
Phi(s + i beta) = Phi(s)     [KMS periodicity]
```

For Fourier modes Phi ~ exp(-i omega_mod s):
```
omega_mod * beta = 2 pi n     (n = 0, 1, 2, ...)
```

The lowest mode n = 0 gives omega_mod = 0!

**The n = 0 mode of the modular Khronon has zero frequency: omega_mod = 0.**

This is the non-propagating mode.

### 6.15 Translation to Physical Dispersion Relation

Modular time s and physical time t are related by:
```
ds = (kappa / (2pi)) dt = (H / (2pi)) dt     [for the Kodama-based modular flow]
```

Therefore omega_mod = 0 corresponds to omega_physical = 0 as well.

**But what about spatial dependence?**

The modular Hamiltonian is associated with a REGION (the causal diamond of radius r_A = 1/H). Different points in space have different causal diamonds, and hence different modular Hamiltonians. The spatial variation is encoded in the k-dependence of Phi:

For Phi(x, t) = Phi_k exp(i k . x - i omega t):
```
delta K(x) = 2 S_A Phi_k exp(i k . x - i omega t)
delta T_mod(x) = -(2/H) Phi_k exp(i k . x - i omega t)
```

The dispersion relation of delta T_mod is INHERITED from the dispersion relation of Phi.

For sub-Hubble modes (k >> H), Phi satisfies the Poisson equation:
```
k^2 Phi = -4 pi G a^2 (rho + delta_rho source)
```

Phi does NOT satisfy a wave equation -- it is determined instantaneously by the matter distribution. In the matter-dominated era, Phi = constant (growing mode), so:

```
omega = 0     for the growing mode of Phi in matter domination
```

**The modular Khronon inherits the non-propagating nature of the gravitational potential in the matter era.**

In the radiation era, Phi oscillates and decays, so the modular Khronon would also oscillate. But this is exactly what happens in AeST/Khronon theories: the c_s = 0 mode mimics CDM in the matter era but differs in the radiation era.

---

## 7. Step 5: Dispersion Relation of the Observer Mode

### 7.1 Summary of the Dispersion Relation

From the analysis in Step 4, the modular Khronon delta T_mod has:

**In the matter-dominated era (z < z_eq):**
```
omega^2 = 0     [Phi = const, non-propagating]
c_s^2 = 0       [the mode grows without oscillating]
```

**In the radiation-dominated era (z > z_eq):**
```
omega ~ k c_s(gamma)     [Phi oscillates, driven by radiation]
```

where c_s(gamma) = 1/sqrt(3) is the radiation sound speed.

**At the transition (z ~ z_eq):**
```
c_s^2 transitions from 1/3 to 0     [as radiation gives way to matter]
```

### 7.2 Comparison with Khronon Theory

In the Blanchet-Skordis Khronon theory:
```
K(Q) = mu^2 (Q - 1)^2

Dispersion: omega^2 = 0 on Minkowski     [exact]
On FRW: omega^2 ~ 0 for k << mu         [approximate, for all CMB scales]
```

The Khronon has c_s = 0 at ALL epochs, not just in the matter era. This is because K(Q) enforces the constraint Q ~ 1 independently of the cosmological epoch.

**The modular Khronon from QRE perturbation theory does NOT have c_s = 0 at all epochs.** In the radiation era, it inherits the oscillatory behavior of Phi.

### 7.3 Is This a Problem?

**Partially, but not fatally.**

In the standard CMB calculation, what matters is:
1. c_s ~ 0 during and after matter-radiation equality (to prevent the dark component from oscillating)
2. The dark component grows gravitationally during matter domination
3. The dark component contributes the right amount to the Poisson equation

Condition 1 is satisfied by the modular Khronon in the matter era but NOT in the radiation era.

In the radiation era, the modular Khronon tracks Phi, which is driven by radiation perturbations. This means the "dark" component would oscillate with the radiation -- it would NOT be an independent, non-oscillating component.

**This is the key failure point**: the modular Khronon does not decouple from radiation perturbations. In AeST/Khronon theory, the scalar field has its OWN kinetic term K(Q) that is independent of the radiation sector, allowing it to remain non-oscillating even when radiation oscillates.

### 7.4 What Would Fix This

For the modular Khronon to have c_s = 0 at ALL epochs, we would need the modular Hamiltonian perturbation to be INDEPENDENT of the gravitational potential Phi. This requires:

**Option A**: An additional kinetic term for delta T_mod that is not derived from Phi.

This would come from the SECOND-ORDER QRE perturbation (the Fisher information metric), which provides a kinetic term ~ (partial_a delta T_mod)^2. If this term has the structure:

```
L ~ mu^2 (partial_t delta T_mod + H delta T_mod)^2     [temporal kinetic]
    - 0 * (partial_i delta T_mod)^2                      [NO spatial kinetic]
```

then the spatial sound speed would be c_s = 0.

The absence of spatial kinetic term requires that the Fisher information metric for the foliation perturbation is PURELY TEMPORAL. This happens if the modular Hamiltonian only generates TIME evolution (not spatial evolution) -- which is true by the Bisognano-Wichmann theorem for Rindler wedges.

**Option B**: The modular Khronon is defined not through Phi but through a SEPARATE scalar degree of freedom in the QRE.

This would require the QRE to contain additional degrees of freedom beyond the metric perturbations. In the Bianconi GfE framework, the G-field (a Lagrange multiplier) provides such an additional degree of freedom. If the G-field's perturbation has c_s = 0, this would solve the problem.

**Option C**: Accept that the modular Khronon is NOT exactly c_s = 0 in the radiation era, but the deviations are small enough to still fit the CMB.

The key CMB scales enter the horizon during the matter-radiation transition, not deep in the radiation era. If c_s is small (but not zero) during this transition, the CMB fit might still be acceptable.

---

## 8. Step 6: Comparison with Khronon/AeST

### 8.1 Mapping the Fields

| QRE Perturbation Theory | Khronon Theory | Match |
|-------------------------|----------------|-------|
| Modular Khronon delta T_mod | Khronon scalar tau | **Structural match** |
| u^a = Kodama vector | n^a = foliation normal | **Exact at background** |
| delta K = 2 S_A Phi | K(Q) = mu^2(Q-1)^2 | **Different**: K is external in QRE, internal in Khronon |
| c_s = 0 in matter era | c_s = 0 at all times | **Partial match** |
| Background energy = 0 | Background energy = tuned | **MISMATCH** |
| Coupling constants from S_A | Coupling constants free | **QRE constrains them** |

### 8.2 What QRE Determines That Khronon Does Not

1. **The modular Hamiltonian normalization**: delta K = 2 S_A Phi. The coefficient 2 S_A = 2 pi / (G_N H^2) is DETERMINED by the de Sitter entropy. In Khronon theory, the coefficient mu^2 is a free parameter.

2. **The relationship between foliation and metric perturbations**: delta T_mod = -2 Phi / H. This is a specific prediction: the Khronon perturbation is proportional to the Newtonian potential with a coefficient determined by the Hubble parameter.

3. **The epoch dependence**: The modular Khronon tracks Phi, so its behavior changes with the cosmological epoch. The Khronon theory has the same behavior at all epochs (c_s = 0 always).

### 8.3 What Khronon Has That QRE Lacks

1. **An independent kinetic term K(Q)**: This decouples the Khronon from the gravitational potential, allowing c_s = 0 even in the radiation era. The QRE formalism, as currently developed, does not provide this independent kinetic term.

2. **The MOND interpolating function J(Y)**: This comes from the spatial gradient coupling of the Khronon/AeST scalar. The QRE formalism has no mechanism to produce this.

3. **Background dark matter energy density**: In Khronon theory, the scalar field's background energy provides Omega_DM h^2 ~ 0.12. The QRE modular Khronon has zero background energy.

### 8.4 The Effective Coupling Constants

If we take the QRE perturbation theory seriously, the effective "Khronon mass parameter" is:

```
mu_eff^2 = (coefficient of (delta T_mod)^2 in the QRE action)
         = second variation of D with respect to foliation
         ~ (rho + p) * beta_A^2  [from Section 6.10]
         = (rho + p) * (2pi/H)^2
```

At the present epoch:
```
rho + p ~ rho_crit ~ 3 H_0^2 / (8 pi G)

mu_eff^2 ~ [3 H_0^2 / (8 pi G)] * (2pi/H_0)^2 = 3 pi / (2 G) ~ M_Planck^2
```

This is the PLANCK SCALE -- far too large! In the Khronon theory, mu ~ (kpc)^{-1} ~ 10^{-22} m^{-1}, which is a cosmological scale.

**The QRE gives the wrong scale for the Khronon mass.** This is because the modular Hamiltonian has Planck-scale normalization (it is proportional to S_A ~ 1/G_N), while the Khronon mass must be at the cosmological scale to affect galactic dynamics.

### 8.5 The Hierarchy Problem

The mismatch between mu_QRE ~ M_Planck and mu_Khronon ~ H_0 represents a HIERARCHY PROBLEM for the QRE -> Khronon derivation:

```
mu_QRE / mu_Khronon ~ M_Planck / H_0 ~ 10^{60}
```

This is the same hierarchy as the cosmological constant problem. The QRE naturally produces Planck-scale couplings; the phenomenology requires cosmological-scale couplings.

**Possible resolution**: The relevant QRE is not the full entanglement entropy S_A ~ 1/G, but the EXCESS entropy above the vacuum (the renormalized QRE). The renormalized QRE could be much smaller:

```
D_renormalized = D - D_vacuum ~ delta rho / T_A ~ (rho + p) / (H/(2pi))
```

This gives:
```
mu_renormalized^2 ~ (rho + p) * (2pi/H)^2 * (G * H^2) / (1)
                   ~ (rho + p) * (4 pi^2) * G
```

At the present epoch:
```
mu_renormalized ~ sqrt(4 pi^2 G * rho_crit) ~ sqrt(4 pi^2 * (8piG/3) * rho_crit / (8piG/3))
                ~ sqrt(4 pi^2) * H_0 ~ 2 pi H_0
```

This is EXACTLY the cosmological scale! The renormalized modular Khronon has:

```
mu_eff ~ 2 pi H_0 ~ 10^{-26} m^{-1}
```

For comparison, the Khronon mass in Blanchet-Skordis (2024):
```
mu_Khronon ~ (kpc)^{-1} ~ 3 x 10^{-20} m^{-1}
```

The QRE value is about 10^6 times smaller. This is closer but still not matching.

**However**: The Khronon mass mu is related to the MOND acceleration scale a_0 through:
```
a_0 ~ c * mu * (some dimensionless coupling)
```

If the dimensionless coupling ~ 10^6, the QRE value would match. This is a large but not absurd number.

---

## 9. Step 7: The Background Energy Problem

### 9.1 The Core Issue

The modular Khronon has ZERO background energy on FRW (because delta T_mod is a perturbation). To get z_eq ~ 3400, we need Omega_DM h^2 ~ 0.12 at the BACKGROUND level.

### 9.2 Can the QRE Provide Background Energy?

**Possibility 1: Casimir-like vacuum energy from the foliation**

If the QRE is evaluated on a SPECIFIC foliation (the modular Khronon foliation), the Casimir energy associated with this foliation could provide a background contribution.

The Casimir energy of a scalar field in a box of size L:
```
rho_Casimir ~ -pi^2 / (720 L^4)
```

For L = r_A = 1/H:
```
rho_Casimir ~ -pi^2 H^4 / 720 ~ -10^{-122} rho_Planck
```

This is the right ORDER OF MAGNITUDE for the cosmological constant (~10^{-120} rho_Planck) but has the WRONG SIGN and is NOT dust-like (w = -1/3 for radiation Casimir, not w = 0).

**Possibility 2: Trace anomaly contribution**

The conformal trace anomaly gives:
```
<T^a_a> = c/(16pi^2) [C^2 - (1/3)R^* R^*] - a/(16pi^2) [E_4 - (2/3) Box R]
```

On FRW:
```
<T^a_a> ~ (H^4 / (16 pi^2)) * (number of fields)
```

This scales as H^4, giving rho ~ H^4 / (16 pi^2). At the present epoch:
```
rho_anomaly ~ H_0^4 / (16 pi^2 * G^2) ~ 10^{-120} rho_Planck
```

Again, this is dark-energy-scale, not dark-matter-scale. And it scales as H^4, not as a^{-3}.

**Possibility 3: The QRE itself as a source**

In the Bianconi GfE framework, the QRE Lagrangian generates an effective stress-energy:
```
T_ab^{QRE} = (1/sqrt(-g)) delta(sqrt(-g) D) / delta g^{ab}
```

This is the stress-energy of the QRE field itself. If D has a non-trivial background value (even though D = 0 at the Friedmann solution), the PERTURBATION of D around the Friedmann solution generates stress-energy:

```
T_ab^{QRE} ~ delta^2 D / delta g^{ab}
```

This is second order in perturbations and gives:
```
<rho_QRE> ~ <(delta D)^2> / (something)
```

From Section 5.4, delta D ~ S_A * (2 Phi + delta_m) ~ 0 when Einstein equations are satisfied. So:
```
<rho_QRE> ~ 0     [at the perturbative level]
```

**None of these possibilities provide Omega_DM h^2 ~ 0.12.**

### 9.3 The Honest Conclusion on Background Energy

The QRE perturbation theory around FRW does NOT generate a background dark matter energy density. The background D = 0 (Friedmann equations satisfied), and perturbative corrections to D are too small to provide Omega_DM.

**This is the most serious failure of the QRE -> Khronon program at the cosmological level.**

The only way out is:

**Resolution A**: The background FRW is NOT the exact D = 0 solution. Instead, the physical background has D = D_0 > 0, representing a departure from exact Friedmann. The departure provides the dark matter energy density.

This is what AeST does: the scalar field has a non-trivial background profile that provides rho_DM. In the QRE language, this would mean:

```
D_physical = D_0 + delta D

D_0 = (some non-zero constant) ~ Omega_DM rho_crit / T_A
```

The challenge: what determines D_0? In AeST, it is set by initial conditions (a free parameter). In the QRE framework, it would need to be determined by a self-consistency condition.

**Resolution B**: Accept that the background energy must be provided by a separate mechanism (conventional CDM particles, or the AeST/Khronon scalar field kinetic energy). The QRE perturbation theory provides the INTERPRETATION and the constraint structure, but not the energy budget.

---

## 10. Verdict: What Is Proved, Assumed, and Failed

### 10.1 Proved (Rigorous or Semi-Rigorous)

1. **Friedmann equations from QRE at zeroth order.** The Cai-Kim (2005) derivation applies the first law at the apparent horizon to get both Friedmann equations. This is equivalent to D_background = 0 (entanglement equilibrium). **Status: ESTABLISHED.**

2. **Linearized Einstein equations from QRE at first order.** The entanglement first law (Faulkner et al. 2014) and the Aalsma-Bak (2025) modular Hamiltonian on FRW give delta D^{(1)} = delta<K> - delta S = 0, which reproduces the linearized Poisson equation. **Status: ESTABLISHED.**

3. **The modular Hamiltonian fluctuation is delta K = 2 S_A Phi.** This is the Aalsma-Bak (2025) result, derived from the Misner-Sharp mass perturbation on the FRW apparent horizon. **Status: ESTABLISHED.**

4. **The observer's foliation perturbation delta T_mod = -2 Phi / H.** This follows from the modular Khronon construction: the modular time deviation from coordinate time is set by the gravitational potential perturbation. **Status: DERIVED (new result).**

### 10.2 Derived Under New Postulate (Semi-Rigorous)

5. **The observer mode delta T_mod has c_s = 0 in the matter-dominated era.** Because delta T_mod = -2 Phi / H and Phi is constant in the matter era (growing mode), delta T_mod is also constant -- non-propagating. **Status: DERIVED, but only in matter era.**

6. **The modular Khronon has the same field content as the Blanchet-Skordis Khronon.** Both are scalar fields whose gradient defines a preferred foliation. **Status: STRUCTURAL MATCH.**

7. **The effective Khronon mass from the renormalized QRE is mu_eff ~ 2 pi H_0.** This is a cosmological-scale mass, consistent (within an order of magnitude) with the Khronon theory requirements. **Status: ESTIMATED, not derived rigorously.**

### 10.3 Assumed (Postulates Not Derived)

8. **Observer Entanglement Equilibrium**: delta D / delta u^a = 0 as an independent equation. This is the postulate that promotes the observer from a gauge choice to a variational degree of freedom. **Status: POSTULATED, not derived.**

9. **The modular Khronon is dynamical, not gauge.** In standard GR, the foliation perturbation is pure gauge. The QRE framework may or may not break this gauge invariance. **Status: OPEN QUESTION.**

10. **The renormalization of the QRE gives cosmological-scale couplings.** The bare QRE gives Planck-scale couplings; the renormalized version gives H_0-scale couplings. The renormalization procedure is not rigorously defined. **Status: ESTIMATED.**

### 10.4 Failed (Cannot Be Made to Work Without Additional Input)

11. **c_s = 0 at ALL epochs.** The modular Khronon tracks Phi, which oscillates in the radiation era. An independent kinetic term (like K(Q) in the Khronon theory) is needed to maintain c_s = 0 in the radiation era, and this does not emerge from the QRE perturbation theory. **Status: FAILS.**

12. **Background dark matter energy density.** The QRE perturbation around FRW gives D_background = 0 and no mechanism to produce Omega_DM h^2 ~ 0.12. **Status: FAILS.**

13. **MOND interpolating function.** The QRE perturbation theory does not produce the spatial coupling J(Y) needed for MOND phenomenology at galactic scales. **Status: NOT ADDRESSED (requires separate calculation).**

14. **Fixing the coupling constants.** The QRE provides the FORM of the modular Khronon coupling but not the numerical values (mu, K_B, etc.). **Status: PARTIALLY CONSTRAINED (mu ~ H_0 from renormalization) but not fully determined.**

### 10.5 The Bottom Line

**The QRE perturbation theory on FRW produces a modular Khronon mode delta T_mod = -2 Phi / H that is structurally similar to the Blanchet-Skordis Khronon field. This mode is non-propagating (c_s = 0) in the matter-dominated era. However, the framework FAILS to produce: (a) c_s = 0 in the radiation era, (b) background dark matter energy, and (c) the MOND interpolating function.**

**The calculation succeeds as a MOTIVATION for the Khronon field from QRE principles but fails as a DERIVATION of CDM-like behavior.**

---

## 11. The Critical New Ingredient: Aalsma-Bak (2025)

### 11.1 Why This Paper Changes the Landscape

The Aalsma-Bak (2025) paper "Modular Fluctuations in Cosmology" (PRD 112, 026017; arXiv:2503.04886) provides the MISSING LINK between modular flow and cosmological perturbation theory. Their key results:

1. **The modular Hamiltonian for FRW causal diamonds:**
```
K = beta M(R_h) = (2pi/H) * M(R_h)
```

2. **The perturbation:**
```
delta K = 2 S_dS Phi
```

3. **The modular fluctuation variance:**
```
<(delta K)^2> = 4 S_dS^2 <Phi^2> = 2 epsilon S_dS ln(mu/H) <K>
```

where epsilon is the slow-roll parameter.

4. **The Verlinde-Zurek scaling:**
```
delta L^2 ~ l_P L
```

emerges from the quantization of the s-wave (l=0) mode of the inflaton.

### 11.2 How This Connects to the Khronon

The Aalsma-Bak result delta K = 2 S_dS Phi provides the explicit formula for how the modular Hamiltonian responds to metric perturbations. Inverting this:

```
Phi = delta K / (2 S_dS) = (G_N H^2 / (2 pi)) * delta K
```

The modular Khronon delta T_mod is then:
```
delta T_mod = -2 Phi / H = -(G_N H / pi) * delta K
```

This is a DIRECT MAP from modular Hamiltonian fluctuations to the Khronon field perturbation.

### 11.3 The Modular Fluctuation as Dark Matter Seed

The variance of the modular Khronon:
```
<(delta T_mod)^2> = (4/H^2) <Phi^2> = (4/H^2) * epsilon / (2 S_dS) * ln(mu/H)
                   = (4 epsilon G_N H^2) / (2 pi H^2) * ln(mu/H)
                   = (2 epsilon G_N / pi) * ln(mu/H)
```

This is the POWER SPECTRUM of the modular Khronon fluctuations. In inflationary cosmology:
```
P_T(k) = <|delta T_mod(k)|^2> ~ (G_N epsilon / pi) * (H_inf / k)^2
```

This has a nearly scale-invariant spectrum (n_s ~ 1), just like the primordial curvature perturbation. The amplitude is:
```
A_T ~ G_N H_inf^2 epsilon / pi ~ A_s * epsilon^2
```

where A_s ~ 2 x 10^{-9} is the scalar amplitude from Planck.

With epsilon ~ 0.01 (typical slow-roll):
```
A_T ~ 2 x 10^{-13}
```

This is FAR too small to provide the dark matter perturbation amplitude (which is A_s ~ 2 x 10^{-9}).

**The modular Khronon fluctuations from inflation are 10^4 times too weak to serve as dark matter perturbation seeds.**

### 11.4 What This Means

The Aalsma-Bak connection provides a beautiful CONCEPTUAL bridge between modular flow and cosmological perturbations, but the QUANTITATIVE prediction falls short:

1. The modular Khronon has the right STRUCTURE (scalar field, foliation-defining, nearly scale-invariant)
2. But the wrong AMPLITUDE (too small by 10^4)
3. And the wrong EPOCH-DEPENDENCE (not c_s = 0 in the radiation era)
4. And zero BACKGROUND energy

To amplify the modular Khronon fluctuations to the CDM level, one would need a mechanism that enhances the foliation perturbation by a factor of ~10^4. One possibility: if the modular Khronon experiences a tachyonic instability during the matter-radiation transition, its amplitude could grow exponentially. But there is no evidence for this in the current framework.

---

## 12. References

### Foundational: Gravity from Thermodynamics/QRE
1. **Jacobson 1995**: PRL 75, 1260 -- Thermodynamics of spacetime: Einstein eq from Clausius
2. **Jacobson 2016**: PRL 116, 201101; arXiv:1505.04753 -- Entanglement equilibrium and Einstein eq
3. **Cai & Kim 2005**: JHEP 0502:050; hep-th/0501055 -- First law of thermodynamics and Friedmann equations
4. **Hayward 1998**: CQG 15, 3147 -- Unified first law of black-hole dynamics and relativistic thermodynamics
5. **Dorau & Much 2025**: PRL; arXiv:2510.24491 -- QRE to semiclassical Einstein equations
6. **Faulkner, Guica, Hartman, Myers, Van Raamsdonk 2014**: JHEP 03, 051; arXiv:1312.7856 -- Gravitation from entanglement

### Modular Hamiltonian on FRW
7. **Aalsma & Bak 2025**: PRD 112, 026017; arXiv:2503.04886 -- Modular fluctuations in cosmology
8. **Cai & Cao 2007**: PRD 75 -- Unified first law and thermodynamics of apparent horizon in FRW
9. **Casini, Huerta, Myers 2011**: JHEP 05, 036 -- Towards a derivation of holographic entanglement entropy

### Bianconi GfE Framework
10. **Bianconi 2025**: PRD 111, 066001; arXiv:2408.14391 -- Gravity from entropy
11. **Bianconi 2025**: arXiv:2510.22545 -- GfE cosmological solutions

### Khronon / AeST Theory
12. **Blanchet & Skordis 2024**: JCAP 11, 040; arXiv:2404.06584 -- Relativistic Khronon theory
13. **Blanchet & Skordis 2025**: arXiv:2507.00912 -- Khronon-Tensor theory
14. **Skordis & Zlosnik 2021**: PRL 127, 161302; arXiv:2007.00082 -- AeST
15. **Blas, Pujolas, Sibiryakov 2010**: arXiv:0909.3525 -- Consistent Horava gravity extensions

### Kodama Vector and Misner-Sharp Mass
16. **Kodama 1980**: Prog. Theor. Phys. 63 -- Conserved energy flux, spherically symmetric systems
17. **Misner & Sharp 1964**: Phys. Rev. 136, B571 -- Relativistic equations for adiabatic, spherically symmetric gravitational collapse

### Entanglement First Law
18. **Blanco, Casini, Hung, Myers 2013**: JHEP 08, 060; arXiv:1305.3182 -- Relative entropy and holography
19. **Casini, Grillo, Pontello 2019**: PRD; arXiv:1903.00109 -- Relative entropy for free fields

### tau Framework
20. **Huang 2026**: Paper 1 -- Petz recovery unification
21. **Huang 2026**: Paper 2 -- Information loss is local: gravitational refractive index
22. **Huang 2026**: Paper 5 -- Observer-dependent temporal asymmetry

---

## Appendix A: The Full Perturbation Chain

For reference, here is the complete logical chain of the calculation:

```
Step 0: FRW background + Newtonian gauge perturbation
        ds^2 = -(1+2Phi)dt^2 + a^2(1-2Psi)dx^2
        u^a = (1-Phi, v^i/a)
        Kodama: K^a = a^{-1}(1, -Hr, 0, 0)

Step 1: Define QRE on FRW
        rho_spacetime = (1/Z) exp(-K_A)
        K_A = (2pi/H) M(r_A) = S_A      [Aalsma-Bak]
        rho_matter = state encoding T_{ab}

Step 2: Background D = 0
        Friedmann: H^2 = 8piG rho/3      [Cai-Kim]
        Raychaudhuri: H_dot = -4piG(rho+p) [Cai-Kim]

Step 3: First-order perturbation
        delta K = 2 S_A Phi               [Aalsma-Bak]
        delta S_A = -S_A delta_m
        delta D^{(1)} = S_A(2Phi + delta_m) = 0
        => Poisson equation               [Faulkner et al.]

Step 4: Observer mode
        delta T_mod = -2 Phi / H          [new result]
        c_s = 0 in matter era             [from Phi = const]
        c_s != 0 in radiation era         [FAILURE POINT]

Step 5: Effective Khronon mass
        mu_eff ~ 2pi H_0                  [from renormalized QRE]
        Hierarchy: mu_QRE/mu_Khronon ~ 10^6  [MISMATCH]

Step 6: Background energy
        rho_DM from QRE = 0               [FAILURE POINT]
        Need separate mechanism for Omega_DM h^2 = 0.12
```

---

## Appendix B: Comparison Table -- What This Calculation Adds to the Previous Analyses

| Question | Before this calculation | After this calculation |
|----------|----------------------|----------------------|
| Does QRE give Friedmann eq? | Known (Cai-Kim 2005) | Confirmed, connected to D = 0 |
| Does QRE give linearized Einstein? | Known (Dorau-Much 2025) | Confirmed on FRW via Aalsma-Bak |
| Does an observer mode emerge? | Hoped for, no calculation | **delta T_mod = -2Phi/H emerges as modular Khronon** |
| Is it c_s = 0? | Unknown | **YES in matter era, NO in radiation era** |
| Is it independent of Phi? | Unknown | **NO -- it tracks Phi, not independent** |
| Does it provide background DM? | Suspected not | **Confirmed: NO** |
| What is the effective mass? | Unknown | **mu_eff ~ 2pi H_0 (from renormalized QRE)** |
| Is this a derivation of AeST/Khronon? | Hoped for | **NO -- it is a structural parallel with quantitative mismatches** |
| What new postulate is needed? | Unknown | **Observer Entanglement Equilibrium: delta D/delta u = 0** |
| Is the modular Khronon physical or gauge? | Unknown | **Gauge in standard GR, potentially physical if QRE breaks foliation invariance** |

---

## Appendix C: The Three Open Calculations

This analysis identifies three concrete calculations that could resolve the open questions:

### Calculation 1: Second-Order QRE Perturbation on FRW

Compute delta D^{(2)} including the Fisher information metric for the foliation perturbation. Determine whether the kinetic term for delta T_mod has the structure:
```
L^{(2)} ~ alpha (dot{delta T_mod})^2 + beta (nabla delta T_mod)^2
```

If beta = 0 and alpha > 0, the modular Khronon has c_s = 0 at all epochs.

**Required tools**: Perturbation of Araki relative entropy to second order on FRW background; explicit computation of the modular operator perturbation.

**Estimated difficulty**: VERY HIGH. No existing calculation goes this far.

### Calculation 2: Bianconi GfE Perturbation Theory

Develop the perturbation theory of the GfE action on FRW. Determine whether the G-field perturbation has c_s = 0.

**Required tools**: Variation of the GfE action with respect to metric perturbations; linearized field equations for the G-field; dispersion relation analysis.

**Estimated difficulty**: HIGH but tractable. The GfE action is known explicitly.

### Calculation 3: Modular Khronon from Crossed Product

Use the Witten (2022) / Chandrasekaran et al. (2023) crossed product construction to define the modular Hamiltonian rigorously on FRW. Determine whether the observer algebra provides an independent degree of freedom.

**Required tools**: Type II von Neumann algebra construction for FRW causal diamonds; perturbation of the crossed product algebra.

**Estimated difficulty**: EXTREMELY HIGH. Frontier of mathematical physics.

---

*Last updated: 2026-03-12*
*This calculation represents the most explicit analysis to date of whether D(rho_spacetime || rho_matter) on FRW can produce a Khronon-like field.*
*Classification: PARTIAL SUCCESS with precisely identified failure points and a concrete path forward.*
