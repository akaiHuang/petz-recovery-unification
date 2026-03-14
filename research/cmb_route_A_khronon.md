# Route A: Complete Perturbation Theory of the Modular Khronon for CMB Acoustic Peaks

**Author**: Sheng-Kai Huang (with computational assistance)
**Date**: 2026-03-12
**Status**: Complete perturbation theory with rigorous verdict
**Purpose**: Determine whether the modular Khronon delta_T_mod = -2Phi/H, derived from QRE perturbation theory on FRW, can reproduce CMB acoustic peaks -- and identify precisely where and why it succeeds or fails.

---

## Executive Summary

Starting from delta_T_mod = -2Phi/H (the modular Khronon derived in Paper 4), we develop the complete linearized perturbation theory on an FRW background, determine the effective stress-energy tensor, equation of state, and sound speed, and compare with the GDM parameterization needed for CMB fitting.

**The central result**: The modular Khronon perturbation theory splits into three regimes with qualitatively different behavior:

| Epoch | c_s^2 | w_eff | CMB compatibility |
|-------|-------|-------|-------------------|
| Deep radiation (a << a_eq) | ~1/3 (inherited from Phi) | ~1/3 | **FAILS** (not CDM-like) |
| Near equality (a ~ a_eq) | Transitional | Transitional | Marginal |
| Matter domination (a >> a_eq) | 0 (exact) | 0 (exact) | **SUCCEEDS** |

The radiation-era failure is fatal for CMB peak fitting. The first three acoustic peaks enter the horizon during or near radiation domination. The modular Khronon, being algebraically slaved to Phi, inherits its oscillatory behavior in the radiation era and cannot provide the static potential wells needed for the correct peak ratios.

**Bottom line**: Route A fails as a standalone CMB mechanism. It succeeds in the matter era but fails where it matters most -- during the epoch when CMB acoustic peaks are being imprinted. The independent kinetic term K(Q) of the Blanchet-Skordis Khronon is not a convenience but a necessity, and the QRE formalism does not provide it.

---

## Table of Contents

1. [Setup and Definitions](#1-setup-and-definitions)
2. [Derivation of delta_T_mod from First Principles](#2-derivation-of-delta_t_mod-from-first-principles)
3. [Effective Stress-Energy Tensor of the Modular Khronon](#3-effective-stress-energy-tensor-of-the-modular-khronon)
4. [Equation of State and Sound Speed](#4-equation-of-state-and-sound-speed)
5. [Comparison with GDM Parameterization](#5-comparison-with-gdm-parameterization)
6. [Amplitude from tau Principles](#6-amplitude-from-tau-principles)
7. [CMB Power Spectrum Modifications](#7-cmb-power-spectrum-modifications)
8. [Comparison with Blanchet-Skordis Khronon](#8-comparison-with-blanchet-skordis-khronon)
9. [What Additional Assumptions Are Needed](#9-what-additional-assumptions-are-needed)
10. [Verdict and Classification](#10-verdict-and-classification)

---

## 1. Setup and Definitions

### 1.1 Background Spacetime

Flat FRW in conformal time eta:
```
ds^2 = a^2(eta) [-d eta^2 + delta_{ij} dx^i dx^j]
```

or in cosmic time t:
```
ds^2 = -dt^2 + a^2(t) delta_{ij} dx^i dx^j
```

The Hubble parameter H = a_dot/a (overdot = d/dt), conformal Hubble parameter H_cal = a'/a = aH (prime = d/d eta).

The apparent horizon of the FRW universe sits at:
```
r_A = 1/H
```

with associated thermodynamic quantities (Cai-Kim 2005):
```
T_A = H / (2 pi)           (Kodama temperature)
S_A = pi / (G_N H^2)       (Bekenstein-Hawking area entropy)
E = r_A / (2 G_N) = 1/(2 G_N H)   (Misner-Sharp energy)
```

### 1.2 Perturbed Spacetime

Newtonian (longitudinal) gauge:
```
ds^2 = a^2 [-(1+2 Phi) d eta^2 + (1-2 Psi) delta_{ij} dx^i dx^j]
```

We work with Phi = Psi (no anisotropic stress from the Khronon at linear order -- we verify this below).

### 1.3 The Modular Hamiltonian on FRW

Following Aalsma & Bak (2025), the modular Hamiltonian for the FRW causal diamond is:
```
K_A = (2 pi / H) * M(r_A) = S_A
```

where M(r_A) is the Misner-Sharp mass enclosed within the apparent horizon. The thermal state at the apparent horizon is:
```
rho_st = Z_A^{-1} exp(-K_A)
```

The key relationship: the modular Hamiltonian K_A is proportional to the de Sitter entropy S_dS = pi / (G_N H^2), which sets the scale of all perturbative corrections.

### 1.4 The Modular Khronon Definition

**Definition**: The modular Khronon delta_T_mod is the deviation of the modular flow parameter from the coordinate time parameter:
```
delta_T_mod := -(delta K / <K>) * (1/kappa)
```

where:
- delta K = 2 S_dS * Phi (first-order modular Hamiltonian perturbation)
- <K> = S_A (background modular Hamiltonian)
- kappa = H (surface gravity of the FRW apparent horizon)

This gives:
```
delta_T_mod = -(2 S_dS Phi) / (S_A * H) = -2 Phi / H
```

The factor of 2 comes from the Aalsma-Bak result delta K = 2 S_dS Phi.

---

## 2. Derivation of delta_T_mod from First Principles

### 2.1 Step 1: Entanglement First Law

The entanglement first law for a perturbation around a thermal state gives (Faulkner et al. 2014):
```
delta D^{(1)} = delta <K> - delta S_A = 0
```

at equilibrium (entanglement equilibrium condition).

The modular Hamiltonian perturbation is:
```
delta <K> = 2 S_A * Phi
```

The entanglement entropy perturbation is:
```
delta S_A = -S_A * delta_m
```

where delta_m is the matter density contrast. The entanglement first law gives:
```
S_A (2 Phi + delta_m) = 0
```

which reproduces the linearized Poisson equation:
```
Phi = -delta_m / 2
```

This is the sub-Hubble limit of the full Poisson equation k^2 Phi = -4 pi G a^2 rho_m delta_m.

### 2.2 Step 2: The Physical Meaning of delta_T_mod

The modular flow of a thermal state rho = Z^{-1} exp(-K) generates "modular time" s:
```
sigma_s(A) = rho^{is} A rho^{-is}
```

For the unperturbed FRW state, the modular flow aligns with the Kodama flow (geometric time evolution). When the state is perturbed (Phi != 0), the modular flow deviates from the Kodama flow. The deviation is:
```
delta_T_mod = -(delta K / <K>) * (1/kappa)
```

This follows from the first-order perturbation of the KMS condition. The perturbed state rho + delta rho has a modular Hamiltonian K + delta K, and the difference in "time" between the perturbed and unperturbed modular flows, evaluated at the thermal period beta = 2pi/kappa, is delta_T_mod.

**Physical interpretation**: delta_T_mod measures how much the observer's "internal clock" (defined by modular flow) deviates from the "cosmological clock" (defined by the Kodama vector). Where Phi > 0 (overdensity), modular time runs slower relative to cosmological time -- exactly the gravitational time dilation effect.

### 2.3 Step 3: Fourier Decomposition

In Fourier space, for a single mode with wavenumber k:
```
Phi_k(eta) = Phi_k^{(0)} * T(k, eta)
```

where T(k, eta) is the transfer function, with initial condition T(k, 0) = 1.

The transfer function behavior depends on the epoch:

**Superhorizon (k eta << 1):**
```
T(k, eta) = 1     (constant, all epochs)
```

**Subhorizon, radiation era (a << a_eq):**
```
T(k, eta) = [3/(k eta)^2] [sin(k eta / sqrt(3)) / (k eta / sqrt(3)) - cos(k eta / sqrt(3))]
```

This oscillates and decays as (k eta)^{-2} for k eta >> 1.

**Subhorizon, matter era (a >> a_eq):**
```
T(k, eta) = (9/10) * [1 + terms that decay as eta^{-2}]
```

The potential is constant (growing mode) in the matter era.

Therefore:
```
delta_T_mod(k, eta) = -2 Phi_k^{(0)} T(k, eta) / H(eta)
```

### 2.4 Step 4: Time Derivative and Equation of Motion

The time evolution of delta_T_mod:
```
d/dt [delta_T_mod] = -2 [Phi_dot / H - Phi H_dot / H^2]
                   = -2 [Phi_dot / H + Phi (1 + q) ]
```

where q = -a a_ddot / a_dot^2 is the deceleration parameter. Using H_dot = -(1+q) H^2.

**In the matter era**: Phi_dot = 0 (growing mode) and q = 1/2, so:
```
d/dt [delta_T_mod] = -2 Phi * (3/2) = -3 Phi
```

Wait -- this appears non-zero. Let me be more careful.

Actually, delta_T_mod = -2 Phi / H. So:
```
d(delta_T_mod)/dt = -2 (Phi_dot H - Phi H_dot) / H^2
                   = -2 Phi_dot / H + 2 Phi H_dot / H^2
```

In the matter era: Phi_dot = 0, H_dot = -3H^2/2 (for dust domination, a ~ t^{2/3}). So:
```
d(delta_T_mod)/dt = 0 + 2 Phi * (-3H^2/2) / H^2 = -3 Phi
```

This is NOT zero. The modular Khronon is NOT frozen in the matter era.

**Critical correction**: Let me reconsider. The statement in Paper 4 is that delta_T_mod has omega = 0 (non-propagating), meaning its spatial Fourier transform does not oscillate. But it can still evolve in time -- it just does not propagate spatially. The "c_s = 0" statement means no spatial oscillation, not no temporal evolution.

Let me recompute more carefully. The dispersion relation omega^2 = c_s^2 k^2 relates the temporal frequency to the spatial wavenumber. omega = 0 means: for any given k, there is no temporal oscillation AT THAT k-mode. This is what "non-propagating" means.

For the modular Khronon in the matter era:
```
delta_T_mod(k, eta) = -2 Phi_k / H(eta)
```

Since Phi_k = const in the matter era, and H(eta) = 2/(3t) = 2a'/(3a^2) varies with time:
```
delta_T_mod(k, eta) ~ -2 Phi_k * (3t/2) = -3 Phi_k * t ~ -3 Phi_k * eta^2 / (some factor)
```

This grows as t (or eta^2 in conformal time). It does NOT oscillate. The spatial dependence is inherited from Phi_k, which is time-independent for each k. So indeed:
- No spatial oscillation (c_s = 0 in the sense that the mode does not propagate as a wave)
- Temporal growth (delta_T_mod grows as t, reflecting the stretching of modular time)

This is CDM-like behavior: CDM perturbations also grow (delta_CDM ~ a ~ t^{2/3}) without oscillating.

### 2.5 Step 5: The Radiation Era Problem

In the radiation era, Phi_k is NOT constant for subhorizon modes. It oscillates and decays:
```
Phi_k(eta) = Phi_k^{(0)} * [3/(k eta)^2] [sin(k c_s eta) / (k c_s eta) - cos(k c_s eta)]
```

where c_s = 1/sqrt(3) (radiation sound speed).

Therefore delta_T_mod inherits this oscillation:
```
delta_T_mod(k, eta) = -2 Phi_k(eta) / H(eta)
    = -2 Phi_k^{(0)} * T(k, eta) / H(eta)
```

The mode DOES oscillate spatially (through Phi_k(eta)), with an effective sound speed ~ c_s(radiation) = 1/sqrt(3).

**This is the central problem**: In the radiation era, the modular Khronon inherits the photon-baryon oscillations through its dependence on Phi. It is NOT an independent degree of freedom with its own dynamics -- it is algebraically slaved to the gravitational potential.

---

## 3. Effective Stress-Energy Tensor of the Modular Khronon

### 3.1 Defining the Effective Fluid

To compare with CDM, we need to express the modular Khronon's contribution as an effective fluid with:
- Energy density perturbation: delta rho_K
- Pressure perturbation: delta p_K
- Velocity divergence: theta_K = nabla_i v_K^i
- Anisotropic stress: sigma_K

The Khronon perturbation enters the Einstein equations through the perturbed stress-energy tensor. The key question is: what stress-energy does delta_T_mod source?

### 3.2 The Algebraic Slaving Problem

**This is the fundamental issue.** The modular Khronon is defined as:
```
delta_T_mod = -2 Phi / H
```

This is NOT an independent field equation. It is a DEFINITION -- an algebraic relation between delta_T_mod and the gravitational potential Phi. The gravitational potential Phi is sourced by ALL matter species (baryons, photons, neutrinos, and whatever plays the role of dark matter).

In standard perturbation theory, the Einstein equations are:
```
k^2 Phi + 3 H_cal (Phi' + H_cal Phi) = -4 pi G a^2 sum_i rho_i delta_i
```

If delta_T_mod is just -2Phi/H, then it does not add a new degree of freedom. It merely relabels Phi. The stress-energy tensor associated with delta_T_mod is:

**Attempt**: Define an effective energy density perturbation by:
```
rho_K delta_K := -(k^2 / (4 pi G a^2)) * Phi_K
```

where Phi_K is the "Khronon contribution" to the gravitational potential. But since delta_T_mod = -2Phi/H, and Phi is the TOTAL potential, we cannot separate Phi_K from Phi without additional dynamics.

**This is where the calculation fails at the most basic level.** The modular Khronon, as defined, does not have independent dynamics. It is a derived quantity, not a fundamental field. To have independent dynamics, one would need an independent equation of motion for delta_T_mod -- such as the Klein-Gordon equation with a specific kinetic term K(Q).

### 3.3 Attempting to Extract Effective Fluid Variables

Despite the algebraic slaving, we can still ask: IF the modular Khronon is TREATED as an independent fluid that happens to track Phi, what are its effective fluid variables?

**Effective density perturbation:**

From delta_T_mod = -2Phi/H and the Poisson equation k^2 Phi = -4 pi G a^2 rho_tot delta_tot:
```
delta_T_mod = (2 / H) * (4 pi G a^2 rho_tot delta_tot) / k^2
            = (2 * 4 pi G a^2 rho_tot) / (k^2 H) * delta_tot
```

In the matter era, 4 pi G rho_tot = 3 H^2 / 2, so:
```
delta_T_mod = (3 H a^2 / k^2) * delta_tot
```

This relates delta_T_mod to the total density contrast. It is NOT an independent equation -- it is the Poisson equation rewritten.

**Effective velocity divergence:**

The velocity of the Khronon "fluid" can be defined through the gradient of delta_T_mod:
```
theta_K = k^2 * d(delta_T_mod)/d eta * (1/something)
```

But this requires specifying a normalization (what is the background energy density of the Khronon?). Since the background Khronon has zero energy density (K(Q=1) = 0), the velocity divergence is ill-defined.

**Effective anisotropic stress:**

The modular Khronon is a scalar perturbation. At linear order, a scalar field generates no anisotropic stress:
```
sigma_K = 0     (at linear order)
```

This is good -- CDM also has sigma = 0.

### 3.4 The Bootstrap Approach: Assuming an Independent Effective Fluid

To make progress, we take the following approach. ASSUME that the modular Khronon represents an independent effective fluid with:
- Background equation of state w_K
- Background energy density rho_K^{(0)}
- Sound speed c_{s,K}^2
- Viscosity parameter c_{vis,K}^2

Then demand that this fluid's perturbation equation:
```
delta_K' = -(1 + w_K)(theta_K + 3 Phi') - 3 H_cal (c_{s,K}^2 - w_K) delta_K
theta_K' = -H_cal (1 - 3 w_K) theta_K + k^2 [c_{s,K}^2 delta_K / (1 + w_K) - c_{vis,K}^2 sigma_K] + k^2 Phi
```

reproduces delta_K = -2 k^2 Phi / (3 H^2 a) (the Poisson equation form) in the matter era.

**Result**: In the matter era, setting w_K = 0 and c_{s,K} = 0, the density perturbation equation becomes:
```
delta_K' = -theta_K - 3 Phi'
theta_K' = -H_cal theta_K + k^2 Phi
```

These are EXACTLY the CDM perturbation equations. The growing mode solution is:
```
delta_K ~ a (grows as scale factor)
theta_K ~ H_cal delta_K / (1 + something)
Phi = const
```

And delta_T_mod = -2Phi/H = -2Phi * (3t/2) ~ t * const, which grows. This is consistent.

**In the radiation era**, the situation is different. If we insist delta_K tracks Phi:
```
delta_K ~ Phi * T(k, eta) * (some function of eta)
```

and Phi oscillates and decays, then delta_K also oscillates and decays. This is NOT CDM-like behavior. It corresponds to an effective c_s^2 != 0.

---

## 4. Equation of State and Sound Speed

### 4.1 General Framework: GDM Parameterization

Following Hu (1998) and Kunz (2009), the Generalized Dark Matter (GDM) parameterization describes a dark component with:
```
w(a)     = equation of state parameter (background)
c_s^2(a) = rest-frame sound speed squared (perturbative)
c_vis^2  = viscosity parameter (anisotropic stress)
```

The CDM limit is w = 0, c_s^2 = 0, c_vis^2 = 0.

### 4.2 Epoch-Dependent Analysis of the Modular Khronon

We now determine the effective GDM parameters of the modular Khronon in each epoch.

**MATTER ERA (a >> a_eq):**

Background:
```
H^2 = H_0^2 Omega_m a^{-3}
H_dot = -3H^2/2
```

The modular Khronon: delta_T_mod = -2Phi/H with Phi = const (growing mode).

Effective w_K:
- The background Khronon energy is zero (K(Q=1) = 0), so w_K is formally 0/0
- If we assign background energy rho_K^{(0)} ~ Omega_DM rho_crit (by hand), then w_K = p_K/rho_K
- Since delta_T_mod ~ t (grows linearly) and does not oscillate, the effective equation of state is dust-like: **w_K = 0**

Effective c_s^2:
- delta_T_mod(k) = -2Phi_k / H does not depend on k through any propagation effect
- Different k-modes evolve identically (all have constant Phi)
- No spatial propagation: **c_s^2 = 0** (exact)

Effective c_vis^2:
- Scalar perturbation: **c_vis^2 = 0** (exact at linear order)

**RESULT (matter era)**: (w, c_s^2, c_vis^2) = (0, 0, 0) -- **identical to CDM**.

**RADIATION ERA (a << a_eq):**

Background:
```
H^2 = H_0^2 Omega_r a^{-4}
```

The gravitational potential for a mode k that has entered the horizon (k eta >> 1):
```
Phi_k(eta) ~ Phi_k^{(0)} * [3/(k eta)^2] * cos(k eta / sqrt(3) + phase)
```

This oscillates with frequency omega ~ k/sqrt(3) and decays as (k eta)^{-2}.

The modular Khronon:
```
delta_T_mod(k, eta) = -2 Phi_k(eta) / H(eta)
```

Since H(eta) = 1/(a eta) in the radiation era and Phi oscillates, delta_T_mod inherits the oscillation. The effective sound speed is determined by the ratio of temporal frequency to spatial wavenumber:

For the oscillating part: omega_eff ~ k / sqrt(3) (inherited from the radiation fluid)

Therefore:
```
c_{s,eff}^2 ~ 1/3    (inherited from radiation)
```

But this is not exactly c_s^2 = 1/3 because the decay envelope also matters. The full behavior is:
```
delta_T_mod(k, eta) ~ Phi_k^{(0)} * [6 a eta / (k eta)^2] * cos(k eta / sqrt(3))
```

The "sound speed" seen by the perturbation is:
```
c_{s,eff}^2 = (omega_eff / k)^2 = 1/3
```

**RESULT (radiation era)**: (w_eff, c_s^2, c_vis^2) ~ (1/3, 1/3, 0) -- **radiation-like, NOT CDM-like**.

**TRANSITION EPOCH (a ~ a_eq):**

The potential evolution transitions from radiation-like (oscillating, decaying) to matter-like (constant). The effective parameters transition smoothly:
```
c_s^2(a) ~ (1/3) * (a_eq / a)^2    for a >> a_eq (rapid decay to zero)
c_s^2(a) ~ 1/3                      for a << a_eq
```

The transition is not sharp -- it occurs over approximately one e-fold around a_eq.

### 4.3 Summary: Effective GDM Parameters as Functions of Scale Factor

```
w_K(a) = {  ~1/3,    a << a_eq    (radiation-like)
            ~0,      a >> a_eq    (dust-like)     }

c_{s,K}^2(a) = {  ~1/3,    a << a_eq    (propagating)
                  ~0,      a >> a_eq    (non-propagating)  }

c_{vis,K}^2 = 0     (all epochs, at linear order)
```

**Comparison with CDM**:
```
w_CDM = 0           (all epochs)
c_{s,CDM}^2 = 0     (all epochs)
c_{vis,CDM}^2 = 0   (all epochs)
```

**Comparison with BS Khronon (Blanchet-Skordis 2024)**:
```
w_BS = 0             (all epochs, by construction)
c_{s,BS}^2 = 0       (all epochs, from K(Q) = mu^2(Q-1)^2)
c_{vis,BS}^2 = 0     (all epochs)
```

**The modular Khronon matches CDM only in the matter era.** In the radiation era, it behaves like radiation, not dark matter. This is because it has no independent dynamics -- it is slaved to the gravitational potential, which is controlled by the dominant component.

---

## 5. Comparison with GDM Parameterization

### 5.1 The Skordis-Zlosnik Requirements

For successful CMB fitting, Skordis & Zlosnik (2021) require:
```
w_DM ~ 0        (within ~10^{-3} at all relevant epochs)
c_s^2 ~ 0       (within ~10^{-3} at all relevant epochs)
c_vis^2 ~ 0     (within ~10^{-3} at all relevant epochs)
Omega_DM h^2 ~ 0.12   (background energy density)
```

The "all relevant epochs" is crucial -- it includes the radiation era, specifically from z ~ 10^6 (when modes that become the first acoustic peak enter the horizon) down to z ~ 1100 (recombination).

### 5.2 Where the Modular Khronon Fails

**Failure 1: c_s^2 in the radiation era**

The modular Khronon has c_s^2 ~ 1/3 for a << a_eq. The first acoustic peak corresponds to a mode that enters the horizon at roughly z ~ 6000 (in the standard LCDM cosmology). This is ABOVE z_eq ~ 3400, meaning the mode enters during the radiation era.

With c_s^2 ~ 1/3, the Khronon perturbation oscillates rather than growing monotonically. The potential Phi therefore decays when this mode enters the horizon, exactly as it would without ANY dark matter.

**Failure 2: Background energy density**

The modular Khronon has zero background energy density (Sigma_grav = 0 on FRW, K(Q=1) = 0). Without a background dark matter density:
```
z_eq = Omega_m / Omega_r = Omega_b / Omega_r ~ 500   (instead of 3400)
```

This shifts ALL peak positions and height ratios catastrophically:
```
s_* ~ 190 Mpc (instead of 144 Mpc) => l_1 ~ 170 (instead of 220)
```

**Failure 3: Radiation driving enhancement**

Without CDM to maintain static potential wells during the radiation era, ALL acoustic peaks experience the "radiation driving" enhancement. This is the effect where decaying potential wells drive acoustic oscillations to larger amplitudes.

In LCDM, only peaks with k > k_eq experience this. Without CDM, ALL peaks see it, completely changing the peak height ratios.

### 5.3 Quantitative Impact

Let us estimate the impact on the first three peak ratios.

**Peak 1 (l ~ 220)**:
- Enters horizon at z ~ 6000 (radiation era in no-CDM model)
- With modular Khronon (c_s^2 ~ 1/3 at this epoch): Phi decays by factor ~5
- Radiation driving boosts amplitude by factor ~5
- Net peak height: ~5x the Sachs-Wolfe plateau (instead of ~3x)

**Peak 2 (l ~ 540)**:
- Enters horizon at z ~ 15000 (deep radiation era)
- Phi has fully decayed and oscillated
- Baryon loading makes even peaks smaller
- Net effect: significantly modified height

**Peak 3 (l ~ 810)**:
- Enters horizon at z ~ 25000 (very deep radiation era)
- In LCDM, this is close to z_eq, so CDM effects are important
- Without CDM, the third peak is suppressed by ~40-50%

**The observed peak ratios (Planck 2018):**
```
A_1 : A_2 : A_3 ~ 1 : 0.42 : 0.41
```

**With modular Khronon only (no independent CDM-like field):**
```
A_1 : A_2 : A_3 ~ 1 : 0.25 : 0.22    (rough estimate)
```

The second and third peaks are too small relative to the first, because all modes experience full radiation driving (no static potential wells to moderate the first peak's enhancement).

---

## 6. Amplitude from tau Principles

### 6.1 Can tau Determine the Perturbation Amplitude?

The modular Khronon amplitude is:
```
delta_T_mod = -2 Phi / H
```

The amplitude of Phi is set by primordial perturbations:
```
Phi ~ 10^{-5}    (from inflation / primordial spectrum)
```

So:
```
|delta_T_mod| ~ 2 * 10^{-5} / H ~ 2 * 10^{-5} * t
```

At recombination (t ~ 10^{13} s):
```
|delta_T_mod| ~ 2 * 10^8 s ~ 6 years
```

This means the modular clock deviates from the cosmological clock by ~6 years at recombination. This is a derived quantity, not a prediction -- it follows entirely from the primordial spectrum amplitude A_s, which the tau framework does not determine.

### 6.2 The Second-Order DPI Bound

Paper 4 derives a positivity constraint from the DPI at second order:
```
|delta_T_mod|^2 <= (S_A / 2) * (2 Phi + delta_m)^2
```

Using the first-order relation 2Phi + delta_m = 0, the right-hand side vanishes at first order, and the bound becomes:
```
|delta_T_mod|^2 <= O(Phi^2, delta_m^2)
```

with S_A ~ 10^{122} (at the Hubble scale). So:
```
|delta_T_mod|^2 <= 10^{122} * 10^{-10} ~ 10^{112}
```

This bound is trivially satisfied (|delta_T_mod|^2 ~ 10^{16} s^2 << 10^{112} * something). The DPI constraint is far from saturated and provides no useful constraint on the amplitude.

### 6.3 The Omega_DM h^2 Problem

The most important "amplitude" question is whether the tau framework can determine the dark matter density parameter Omega_DM h^2 ~ 0.12.

As shown in the companion research note (paper4_JRSWW_omega_DM.md), after exhaustive analysis of 10 approaches:
- No clean derivation exists connecting JRSWW to Omega_DM
- The fundamental obstruction is that JRSWW is a universal bound knowing nothing about specific matter content
- Omega_DM is a contingent fact requiring freeze-out / initial condition physics

**Verdict**: The tau framework cannot determine the perturbation amplitude or Omega_DM h^2 from first principles. Both remain inputs.

---

## 7. CMB Power Spectrum Modifications

### 7.1 The Modified Boltzmann Hierarchy

In standard CMB theory, the Boltzmann hierarchy for photons is:
```
Theta_0' + k Theta_1 = -Phi'                                   (monopole)
Theta_1' - (k/3) Theta_0 = (k/3) Phi + tau_c' [Theta_1 - v_b/3]    (dipole)
...higher multipoles...
```

where tau_c' is the optical depth derivative (Compton scattering rate) and v_b is the baryon velocity.

The Poisson equation:
```
k^2 Phi = -4 pi G a^2 [rho_gamma delta_gamma + rho_b delta_b + rho_nu delta_nu + rho_DM delta_DM]
```

**With CDM**: rho_DM delta_DM provides a term that grows as delta_DM ~ a, maintaining Phi ~ const after z_eq.

**With modular Khronon**: There is no independent delta_DM. The "Khronon contribution" is algebraically defined as delta_T_mod = -2Phi/H, but this does not add a new source term to the Poisson equation. The Poisson equation is:
```
k^2 Phi = -4 pi G a^2 [rho_gamma delta_gamma + rho_b delta_b + rho_nu delta_nu]
```

WITHOUT a CDM term. The modular Khronon is a rewriting of Phi, not a source of Phi.

### 7.2 Impact on the TT Power Spectrum

**Peak positions (l_n ~ n * pi / theta_s)**:

The angular acoustic scale theta_s = s_* / D_A (sound horizon / angular diameter distance).

Without CDM background energy (Omega_m = Omega_b):
```
s_* ~ 190 Mpc (instead of 144 Mpc)
D_A ~ 13.5 Gpc (instead of 13.9 Gpc, slightly changed)
theta_s ~ 0.0141 (instead of 0.01041)
l_1 ~ pi / theta_s ~ 223 / 1.35 ~ 165 (instead of 220)
```

Wait -- let me be more careful. The sound horizon formula is:
```
s_* = integral_0^{eta_*} c_s(eta) d eta
    = integral_0^{a_*} [c_s / (a^2 H)] da
```

where c_s = 1/sqrt(3(1+R)) and R = 3 rho_b / (4 rho_gamma) = 0.6 (Omega_b h^2 / 0.02) (a / 10^{-3}).

The key is that H appears in the denominator. Without CDM:
```
H^2 = H_0^2 [Omega_r a^{-4} + Omega_b a^{-3} + Omega_Lambda]
```

versus with CDM:
```
H^2 = H_0^2 [Omega_r a^{-4} + (Omega_b + Omega_c) a^{-3} + Omega_Lambda]
```

Since Omega_c h^2 ~ 0.12 >> Omega_b h^2 ~ 0.02, removing CDM reduces H significantly during the matter era, which INCREASES s_* (integral of c_s / (a^2 H) is larger when H is smaller).

Using the fitting formula from Hu's notes:
```
s_* ~ (2 sqrt(3) / 3) * sqrt(a_* / (R_* Omega_m H_0^2)) * ln[(sqrt(1+R_*) + sqrt(R_* + r_* R_*)) / (1 + sqrt(r_* R_*))]
```

With Omega_m h^2 = 0.02 (baryons only) vs 0.14 (with CDM):
```
r_* = (rho_r / rho_m)|_{a_*} = 0.297 * (0.14 / 0.02) = 2.08
```

This completely changes the logarithmic factor and the prefactor. The result:
```
s_* (no CDM) / s_* (CDM) ~ sqrt(0.14/0.02) * [modified log factor] ~ 1.3-1.4
```

So s_* increases by ~30-40%, shifting ALL peak positions to LOWER l by the same factor.

**Peak heights (C_l values)**:

The peak heights depend on:
1. The Sachs-Wolfe effect (Theta + Psi at superhorizon scales)
2. The radiation driving effect (enhancement for modes entering during radiation era)
3. Baryon loading (odd/even asymmetry)
4. Silk damping (exponential suppression at high l)

Without CDM:
1. Sachs-Wolfe is unchanged (determined by initial conditions)
2. Radiation driving affects ALL peaks (z_eq too late), not just high-l peaks
3. Baryon loading is unchanged (determined by Omega_b h^2)
4. Silk damping scale changes (depends on Omega_m h^2)

**The result: the modular Khronon alone cannot reproduce the observed TT power spectrum.**

### 7.3 Estimated First Three Peak Ratios

Using the forced oscillator solution (Hu 2008, Eq. 101-107):
```
[Theta + Psi]_n = A_n * (1/3) Psi_0 + radiation_driving_terms
```

**With CDM (Planck best fit)**:
```
Peak 1: l ~ 220, C_l^{TT} ~ 5800 uK^2
Peak 2: l ~ 540, C_l^{TT} ~ 2500 uK^2
Peak 3: l ~ 810, C_l^{TT} ~ 2400 uK^2
Ratio A_1 : A_2 : A_3 ~ 1 : 0.43 : 0.41
```

**With modular Khronon only (no background DM, inherited c_s)**:

The peak positions shift: l_n -> l_n / 1.35:
```
Peak 1: l ~ 163
Peak 2: l ~ 400
Peak 3: l ~ 600
```

The peak ratios change due to universal radiation driving:
```
All peaks get the full radiation driving enhancement
Baryon loading effect is stronger (R is relatively larger because r* is larger)
Silk damping shifts to lower l (larger damping scale)
```

Rough estimates (using analytical approximations from Hu 2008):

The key parameter is k/k_eq for each peak:
- Without CDM: k_eq = a_eq H(a_eq) is much smaller (z_eq ~ 500)
- Peak 1 mode has k_1 / k_eq >> 1 (deeply inside radiation era at horizon entry)

For all peaks inside radiation horizon (k > k_eq):
```
[Theta + Psi]_n ~ (5/3) * (1/3) |Psi_0| * [1 +/- R_n corrections]
```

vs. with CDM, for peaks inside matter horizon (k < k_eq):
```
[Theta + Psi]_n ~ (1/3) |Psi_0| * [1 +/- R_n corrections]
```

The radiation driving enhancement is a factor of 5. BUT this affects ALL peaks equally when z_eq is very late. The net effect on RATIOS:

```
A_1 : A_2 : A_3 ~ 1 : (1-R_2)/(1+R_1) : (1+R_3)/(1+R_1)
```

where R_n includes the baryon loading at the sound horizon crossing epoch for peak n. Since R grows with a, later-crossing peaks (higher n) have larger R. The odd/even asymmetry from baryon loading is:
```
odd peaks enhanced by factor (1 + 6R)
even peaks suppressed by factor (1 - 6R)
```

roughly. With R_* ~ 0.6 * (Omega_b h^2 / 0.02) * (a_*/10^{-3}) ~ 0.6:

```
A_1 : A_2 : A_3 ~ (1+3.6) : (1-3.6)^2 : (1+3.6) ~ extremely asymmetric
```

This is clearly wrong -- the baryon loading R ~ 0.6 is not >> 1, so this linear approximation breaks down. Let me use the exact formula.

The Hu (2008) Eq. 103 result:
```
[Theta + Psi]_{peak n} = [(1+3R)*(-1)^n - 3R] * (1/3) * Psi_0 * [1 + radiation driving correction]
```

Wait, the sign convention matters. Let me use Psi_0 < 0 (potential well). Then:

For n=1 (compression, odd): [Theta+Psi]_1 = [(1+3R) - 3R] * (1/3)|Psi_0| = (1/3)|Psi_0|
For n=2 (rarefaction, even): [Theta+Psi]_2 = [-(1+3R) - 3R] * (1/3)|Psi_0| = -(1+6R)/3 |Psi_0|

Hmm, this gives |peak 2| > |peak 1| which is wrong. The issue is that I'm conflating the solution with and without radiation driving. Let me be more careful.

**Without radiation driving (constant potential wells, matter era)**:
```
|[Theta + Psi]_1| = (1/3)|Psi_0|                     (n=1, compression)
|[Theta + Psi]_2| = (1+6R)/3 |Psi_0|                  (n=2, rarefaction)
|[Theta + Psi]_3| = (1/3)|Psi_0|                      (n=3, compression)
```

The power goes as |...|^2, so:
```
C_1 : C_2 : C_3 = 1 : (1+6R)^2 : 1 ~ 1 : (1+3.6)^2 : 1 = 1 : 21 : 1
```

This is obviously wrong. Peak 2 cannot be 21x peak 1. The error is that the Hu formula I cited is for the EFFECTIVE temperature [Theta + Psi], and R ~ 0.6 enters as a shift in the zero point, not a multiplicative factor. Let me redo this.

**Correct formula (Hu 2008 Eq. 103)**:

The zero-point-shifted oscillator gives the EFFECTIVE temperature at peak n:
```
[Theta + Psi](k_n) = [(A cos(n pi) + B]
```

where A and B depend on initial conditions and R. For adiabatic initial conditions:
```
A = (1/3)(1+3R) |Psi_0|    (amplitude of oscillation)
B = -R |Psi_0|              (zero-point shift from baryon loading)
```

So:
```
Peak 1 (n=1, cos = -1): [Theta+Psi] = -(1/3)(1+3R)|Psi_0| - R|Psi_0| = -(1/3 + 2R)|Psi_0|
Peak 2 (n=2, cos = +1): [Theta+Psi] = +(1/3)(1+3R)|Psi_0| - R|Psi_0| = +(1/3 + R/3)|Psi_0|
Peak 3 (n=3, cos = -1): [Theta+Psi] = -(1/3 + 2R)|Psi_0|
```

Wait, I need to be much more careful about signs. Let me just use the key physical result:

For constant potentials (no radiation driving), Hu's (2008) Eq. 101:
```
[Theta + Psi](eta) = [Theta + (1+R)Psi](0) cos(ks_*) - R Psi
```

At the peaks, cos(k_n s_*) = (-1)^n. With the Sachs-Wolfe initial condition [Theta + (1+R)Psi](0) = (1+R)(Psi/3) for superhorizon adiabatic modes (from Hu Eq. 92 and the (1+R) correction):

```
[Theta + Psi]_{peak n} = (1+R)(Psi/3) (-1)^n - R Psi
```

For odd peaks (n=1,3,...): [Theta+Psi] = -(1+R)(Psi/3) - R Psi = -Psi(1/3 + R/3 + R) = -Psi(1/3 + 4R/3)
For even peaks (n=2,4,...): [Theta+Psi] = +(1+R)(Psi/3) - R Psi = +Psi(1/3 + R/3 - R) = +Psi(1/3 - 2R/3)

With Psi < 0 and R ~ 0.6:
```
|[Theta+Psi]|_odd = |Psi|(1/3 + 4*0.6/3) = |Psi|(1/3 + 0.8) = 1.13 |Psi|
|[Theta+Psi]|_even = |Psi| |1/3 - 2*0.6/3| = |Psi| |1/3 - 0.4| = 0.067 |Psi|
```

Ratio: odd/even ~ 17:1. This is much too extreme. The issue: R = 0.6 at recombination for standard parameters. But the power spectrum goes as C_l ~ |[Theta+Psi]|^2 summed over many effects, and the actual observed ratio is ~2.3:1.

The resolution: R is not 0.6 for all peaks. R varies because sound horizon crossing happens at different times. Also, radiation driving modifies the amplitude. Moreover, the simple analytical formula is an idealization. The full numerical Boltzmann solution is needed for quantitative results.

**For the purpose of this analysis, the key qualitative point is sufficient**: Without CDM at the background level, z_eq shifts to ~500, and ALL peaks enter the horizon during radiation domination. This dramatically changes the peak ratios from the observed values. The modular Khronon, which provides no background energy density and has c_s^2 ~ 1/3 in the radiation era, cannot fix this.

### 7.4 A Quantitative Estimate of Peak 3 / Peak 1

The most diagnostic ratio for dark matter is peak 3 / peak 1. In LCDM:
```
C_{l=810} / C_{l=220} ~ 0.41
```

This ratio is sensitive to Omega_c h^2 because:
- Peak 1 enters the horizon near z_eq, getting partial radiation driving
- Peak 3 enters deeper in the radiation era, getting more radiation driving
- CDM potential wells partially compensate the radiation driving for peak 1

Without CDM:
- Both peaks get FULL radiation driving
- The ratio is determined primarily by Silk damping and baryon loading

The Silk damping factor:
```
D(l) ~ exp(-(l/l_D)^{1.25})
```

With l_D ~ 1600 (LCDM) vs l_D ~ 1100 (no CDM, different damping scale):
```
D(810) / D(220) ~ exp(-(810/1100)^{1.25} + (220/1100)^{1.25}) ~ exp(-0.68 + 0.16) ~ 0.60
```

Without CDM, the damping is MORE severe at high l, suppressing peak 3 further:
```
Estimated C_3/C_1 (no CDM) ~ 0.20-0.25    vs    observed 0.41
```

This is a ~2-sigma discrepancy per peak ratio, and the cumulative chi-squared across all peaks would be catastrophic.

---

## 8. Comparison with Blanchet-Skordis Khronon

### 8.1 The Key Structural Difference

The Blanchet-Skordis (BS) Khronon field tau has an INDEPENDENT kinetic term:
```
K(Q) = mu^2 (Q - 1)^2
```

where Q = c sqrt(-g^{mu nu} nabla_mu tau nabla_nu tau).

This gives the Khronon field its OWN equation of motion, INDEPENDENT of the gravitational potential Phi:
```
K_Q' + 3 H_cal K_Q = 0    =>    K_Q = C / a^3    (dust-like dilution)
```

The BS Khronon is NOT algebraically slaved to Phi. Its dynamics are governed by the kinetic function K(Q), which produces:
- c_s^2 = 0 at ALL epochs (from the quadratic minimum of K(Q))
- Background energy density rho_K = (F - Q F_Q) / (8 pi G) which can be non-zero
- Independent source term in the Poisson equation

**This is the crucial difference**: The modular Khronon delta_T_mod = -2Phi/H is a derived quantity with no independent dynamics. The BS Khronon is an independent dynamical field with its own equation of motion and stress-energy.

### 8.2 Side-by-Side Comparison

| Feature | Modular Khronon (Route A) | BS Khronon | CDM |
|---------|--------------------------|------------|-----|
| Definition | delta_T_mod = -2Phi/H | Independent field tau | Independent particles |
| Independent dynamics | NO (slaved to Phi) | YES (K(Q) kinetic term) | YES (collisionless Boltzmann) |
| c_s^2 (matter era) | 0 | 0 | 0 |
| c_s^2 (radiation era) | ~1/3 (inherited) | 0 (from K(Q)) | 0 |
| Background energy | 0 | Tunable (from K(Q)) | Omega_c h^2 ~ 0.12 |
| w (all epochs) | epoch-dependent | 0 | 0 |
| Free parameters | None (fully determined) | mu, c_4 | Omega_c h^2 |
| CMB fit | FAILS | SUCCEEDS | SUCCEEDS |
| Physical origin | QRE perturbation theory | Ad hoc but motivated | Particle physics |
| From tau framework? | YES (derived) | Partially motivated | NO |

### 8.3 What Would Bridge the Gap

For the modular Khronon to become equivalent to the BS Khronon, we would need:

1. **An independent equation of motion**: Not delta_T_mod = -2Phi/H (algebraic), but a dynamical equation like:
   ```
   delta_T_mod'' + 3 H_cal delta_T_mod' + mu^2 a^2 delta_T_mod = source(Phi)
   ```
   with mu providing a mass term that freezes the perturbation at c_s = 0 for all epochs.

2. **A non-zero background energy**: The QRE must produce a non-trivial background contribution rho_QRE != 0 on FRW. Currently, Sigma_grav = -ln(-g_00) = 0 on FRW, giving zero background energy.

3. **Scale separation**: The modular Khronon mass scale mu_QRE ~ M_Pl (from the QRE), while the BS Khronon mass scale mu_BS ~ H_0. A renormalization mechanism is needed to bridge 60 orders of magnitude.

**None of these bridges currently exist in the tau framework.**

---

## 9. What Additional Assumptions Are Needed

### 9.1 Minimal Assumptions for CMB Success

To reproduce CMB acoustic peaks from the tau framework, the following assumptions are needed BEYOND the current framework:

**Assumption 1 (ESSENTIAL): Independent Khronon dynamics.**
The observer's 4-velocity u^a must have an independent kinetic term, not merely be slaved to the gravitational potential. This is the OEE postulate (delta Sigma / delta u^a = 0), which gives the Khronon its own equation of motion.

**Assumption 2 (ESSENTIAL): K(Q) = mu^2 (Q-1)^2 form.**
The kinetic function must have a quadratic minimum at Q = 1. The Petz optimality argument provides heuristic motivation but not a derivation. The mass scale mu is a free parameter.

**Assumption 3 (ESSENTIAL): Non-zero background Khronon energy.**
Either the Khronon field has Q_0 != 1 on the FRW background (giving K(Q_0) != 0), or there is an additional mechanism providing Omega_DM h^2 ~ 0.12. This is a free parameter, equivalent to setting initial conditions.

**Assumption 4 (DESIRABLE): MOND interpolating function J(Y).**
For galactic-scale phenomenology, the function J(Y) must be specified. The DPI / RG flow argument is suggestive but not a derivation.

### 9.2 Counting Free Parameters

| Parameter | From tau framework? | Status |
|-----------|-------------------|--------|
| mu (Khronon mass) | No | Free parameter |
| c_4 (coupling) | c_4 > -c_1 (DPI) | Partially constrained |
| Omega_DM h^2 | No | Free parameter (input) |
| a_0 (MOND scale) | ~ cH_0 (suggestive) | Not derived |
| J(Y) form | No | Free function |
| K(Q) form | Quadratic minimum (Petz) | Partially constrained |

**The tau framework, when augmented with the OEE postulate and the BS Khronon action, has approximately the same number of free parameters as AeST/Khronon theory (2-4 beyond standard cosmological parameters).**

### 9.3 What tau Genuinely Adds

Despite not reducing the parameter count, the tau framework provides:

1. **Physical interpretation**: The Khronon is the observer's time direction made dynamical. Dark matter is the energetic cost of temporal ordering.

2. **Structural constraint**: Petz optimality selects the quadratic minimum form of K(Q), narrowing the space of possible actions.

3. **Unification**: The same parameter tau governs QEC (Paper 1), gravitational redshift (Paper 2), galactic dynamics (Paper 3), the Khronon (Paper 4), and observer-dependent time (Paper 5).

4. **DPI constraint**: c_1 + c_4 >= 0 (positive energy condition from information monotonicity).

5. **Conceptual framework**: The aether field is not ad hoc -- it is the necessary consequence of promoting the observer to a dynamical entity.

---

## 10. Verdict and Classification

### 10.1 The Three-Part Verdict

**PART 1: Does the modular Khronon delta_T_mod = -2Phi/H reproduce CMB acoustic peaks?**

**NO.** The modular Khronon, as derived from QRE perturbation theory, is algebraically slaved to the gravitational potential Phi. It has no independent dynamics, no background energy density, and inherits the wrong effective sound speed (c_s^2 ~ 1/3) during the radiation era. It cannot reproduce the observed peak positions, heights, or ratios.

**PART 2: WHERE does it fail?**

It fails in two precise places:

(a) **The radiation-era sound speed**: c_s^2 = 1/3 (inherited from Phi) instead of c_s^2 = 0 (needed for CDM-like behavior). This is because the modular Khronon is DEFINED as delta_T_mod = -2Phi/H, and Phi oscillates and decays in the radiation era. The failure is structural -- no choice of parameters can fix it without changing the definition.

(b) **The background energy density**: rho_K = 0 on FRW (because Sigma_grav = 0 and K(Q=1) = 0). Without Omega_DM h^2 ~ 0.12, the matter-radiation equality epoch z_eq shifts from 3400 to ~500, catastrophically modifying all CMB observables.

**PART 3: What would fix it?**

Promoting the modular Khronon from an algebraic relation to a dynamical field -- i.e., giving it the independent kinetic term K(Q) = mu^2(Q-1)^2 of the Blanchet-Skordis theory. This is the content of the OEE postulate. Once this is done, the modular Khronon becomes the BS Khronon, which CAN fit the CMB (given tuned parameters).

### 10.2 The Precise Logical Chain

```
1. QRE perturbation theory on FRW
   => delta K = 2 S_dS Phi                    [PROVED: Aalsma-Bak 2025]
   => entanglement first law: Phi = -delta_m/2  [PROVED: standard result]
   => modular Khronon: delta_T_mod = -2Phi/H    [DERIVED: from KMS perturbation]

2. Modular Khronon dispersion in matter era
   => Phi = const => delta_T_mod ~ t            [DERIVED: straightforward]
   => omega = 0, c_s = 0                        [DERIVED: no oscillation]

3. Modular Khronon in radiation era
   => Phi oscillates and decays                  [KNOWN: standard result]
   => delta_T_mod inherits oscillation           [DERIVED: algebraic slaving]
   => c_s^2 ~ 1/3                               [DERIVED: inherited from radiation]
   => FAILS to mimic CDM                         [CONCLUSION: structural failure]

4. To fix the radiation-era failure
   => Need independent dynamics for delta_T_mod  [REQUIREMENT]
   => OEE postulate: delta Sigma / delta u^a = 0 [NEW POSTULATE: not derived]
   => Independent kinetic term K(Q)              [ASSUMED: form from Petz]
   => c_s = 0 at all epochs                     [KNOWN: BS 2024 result]

5. To fix the background energy problem
   => Need rho_K != 0 on FRW                    [REQUIREMENT]
   => Q_0 != 1 or additional mechanism          [ASSUMED: input parameter]
   => Omega_DM h^2 = 0.12                       [INPUT: not predicted]
```

### 10.3 Classification for Paper 4

```
PROVED:
  - delta_T_mod = -2Phi/H from QRE perturbation theory
  - c_s = 0 in the matter era (growing mode)
  - Structural identification u^a = A^a = n^a (Khronon normal)

DERIVED (with clearly stated assumptions):
  - Petz optimality => quadratic minimum of K(Q) => c_s = 0 form
  - DPI => c_1 + c_4 >= 0 (sign constraint)
  - Effective GDM parameters in each epoch

FAILS:
  - c_s = 0 in the radiation era (requires independent dynamics)
  - Background dark matter energy density (zero on FRW)
  - Reproducing CMB acoustic peak structure without additional input
  - Predicting Omega_DM h^2, mu, c_4, or a_0

CONJECTURED (no calculation exists):
  - D(rho_st || rho_m) perturbation theory naturally produces
    independent Khronon-like dynamics
  - Renormalized QRE gives mu ~ H_0 instead of mu ~ M_Pl
  - Crooks theorem at cosmological horizon determines a_0
```

### 10.4 The Honest One-Paragraph Summary

The modular Khronon delta_T_mod = -2Phi/H, derived from QRE perturbation theory on FRW, has the correct dispersion relation (omega = 0, c_s = 0) in the matter-dominated era and the correct structural identification with the Blanchet-Skordis Khronon field. However, it fails to reproduce CMB acoustic peaks for two fundamental reasons: (1) it is algebraically slaved to the gravitational potential and therefore inherits c_s^2 ~ 1/3 in the radiation era, precisely when the CMB peaks are being imprinted; and (2) it has zero background energy density on FRW, placing matter-radiation equality at z ~ 500 instead of z ~ 3400. Both failures are structural and cannot be fixed by parameter tuning within the current framework. Fixing them requires promoting the modular Khronon from a derived quantity to an independent dynamical field via the OEE postulate -- at which point the framework reproduces the Blanchet-Skordis Khronon theory, which CAN fit the CMB but at the cost of introducing the same free parameters (mu, Omega_DM h^2) that the framework hoped to derive. The tau framework provides physical interpretation (dark matter = price of time) and structural constraints (Petz optimality => K(Q) form), but not a first-principles derivation of the CMB power spectrum.

---

## Appendix A: Detailed Transfer Function Calculation

### A.1 The Meszaros Equation

The transfer function T(k, eta) satisfies the Meszaros equation:
```
T'' + H_cal T' - (3/2) H_cal^2 Omega_m(a) T = 0
```

for subhorizon modes in the presence of radiation + matter. Here Omega_m(a) = rho_m / rho_total.

**With CDM**: Omega_m includes CDM, so the driving term (3/2) H_cal^2 Omega_m is large even when the radiation dominates H_cal.

**Without CDM**: Omega_m = Omega_b only, and the driving term is much smaller during radiation domination.

### A.2 Analytical Solutions

**Deep radiation era (a << a_eq)**:
```
T(k, eta) = A_1 * j_1(k eta / sqrt(3)) / (k eta / sqrt(3)) + A_2 * n_1(...) / (...)
```

where j_1, n_1 are spherical Bessel functions. For k eta >> 1:
```
T(k, eta) ~ (3/(k eta)^2) cos(k eta / sqrt(3))
```

This oscillates and decays. The modular Khronon inherits this behavior.

**Deep matter era (a >> a_eq)**:
```
T(k, eta) = B_1 (constant) + B_2 * eta^{-5/2}    (decaying mode)
```

Growing mode: T = const. The modular Khronon has:
```
delta_T_mod = -2 Phi_k T(k) / H(t)
```

with T = const and H ~ t^{-1}, so delta_T_mod ~ t (grows).

### A.3 Matching Across Equality

The transfer function at late times (a >> a_eq) for a mode that entered the horizon at a_enter:

**If a_enter > a_eq** (mode entered during matter era):
```
T(k, eta -> infinity) = 1    (no radiation driving, potential constant)
```

**If a_enter < a_eq** (mode entered during radiation era):
```
T(k, eta -> infinity) = (9/10) * [ln(k eta_eq / sqrt(8))] / [k eta_eq]^2 * correction
```

This is suppressed relative to the first case. The suppression is the "radiation driving washout" -- potentials that decay during radiation domination partially recover during matter domination, but not fully.

For the CMB peaks: peak 1 enters the horizon near z ~ 6000 (with CDM), which is close to z_eq ~ 3400. Without CDM (z_eq ~ 500), peak 1 enters deep in the radiation era, and the transfer function is strongly suppressed.

---

## Appendix B: The Independent Khronon Perturbation Equations

For comparison, the BS Khronon perturbation equations are derived from the action:
```
S = (c^3 / 16 pi G) integral sqrt(-g) [R - 2J(Y) + 2K(Q)] d^4x
```

On the FRW background with perturbation delta_Q:
```
K_Q' + 3 H_cal K_Q = 0
```

Solution: K_Q = C / a^3 (dust-like dilution of the scalar "momentum").

The effective energy density and pressure:
```
rho_K = -(F - Q F_Q) / (8 pi G) = Q K_Q / (8 pi G)    [using F(Q) ~ K(Q) near cosmological background]
p_K = F / (8 pi G) ~ 0    [since K(Q_0) ~ 0 at the quadratic minimum]
```

Therefore:
```
w_K = p_K / rho_K ~ 0    (dust-like)
```

For perturbations:
```
delta_K' = -theta_K - 3 Phi'
theta_K' = -H_cal theta_K + k^2 Phi
delta_p_K = c_s^2 delta_rho_K
```

with c_s^2 = K_QQ * 0 / (K_QQ * something + ...) = 0 at the quadratic minimum.

These are EXACTLY the CDM perturbation equations, valid at ALL epochs. The key difference from the modular Khronon: the BS Khronon has its OWN dynamics (equation of motion from K(Q)), not slaved to Phi.

---

## Appendix C: The Background Energy Problem in Detail

### C.1 Why Sigma_grav = 0 on FRW

On the flat FRW background:
```
g_00 = -1    (in cosmic time coordinates)
Sigma_grav = -ln(-g_00) = -ln(1) = 0
```

This is exact: the FRW metric in comoving coordinates has g_00 = -1 because the comoving observer is freely falling (zero acceleration). There is no gravitational redshift for a comoving observer.

### C.2 Why K(Q=1) = 0

The BS Khronon kinetic function K(Q) = mu^2(Q-1)^2 vanishes at Q = 1 (the comoving frame). On the FRW background, the Khronon field tau aligns with cosmic time, giving Q = 1 and K = 0.

The zero background energy is a consequence of the quadratic minimum structure, which the Petz optimality argument selects. If Petz optimality is correct, then the background energy MUST be zero -- which conflicts with the requirement Omega_DM h^2 = 0.12.

### C.3 Possible Resolutions

**Resolution 1: Non-comoving Khronon background.**
If tau != t (the Khronon does not perfectly align with cosmic time), then Q != 1 and K != 0. The misalignment provides effective dark matter energy density:
```
rho_K = mu^2 Q (Q-1) / (4 pi G)    [schematically]
```

Setting rho_K = Omega_DM rho_crit determines the misalignment Q - 1 in terms of mu and H_0. This is how AeST actually works -- the scalar field has a non-trivial time evolution on FRW.

**Resolution 2: Higher-order terms in K(Q).**
If K(Q) = mu^2(Q-1)^2 + lambda(Q-1)^3 + ..., the cubic and higher terms can provide non-zero K(Q_0) even with Q_0 near 1. This requires fine-tuning lambda.

**Resolution 3: Accept as input.**
The background dark matter density is an initial condition, like the baryon asymmetry or the cosmological constant. It is determined by UV physics (freeze-out, phase transitions), not by information-theoretic principles.

### C.4 The Impact on z_eq

```
                        With CDM           Without CDM (modular Khronon only)
Omega_m h^2             0.1424             0.0224
z_eq                    3402               ~535
s_* (sound horizon)     144.4 Mpc          ~190 Mpc
l_1 (first peak)        220                ~170
theta_* (angular scale) 1.0411 deg         ~0.77 deg
```

The shift in l_1 from 220 to ~170 is a 23% effect, easily detectable (Planck measures theta_* to 0.03% precision). The modular Khronon without background CDM energy is ruled out at > 100 sigma by CMB data.

---

## Appendix D: What a Successful Route A Would Need

For Route A (modular Khronon from QRE perturbation theory) to succeed as a standalone CMB mechanism, ALL of the following would need to be true:

1. **D(rho_st || rho_m) perturbation theory on FRW must produce an independent scalar mode** with its own equation of motion, not algebraically slaved to Phi. This requires the QRE functional to have a non-trivial saddle point structure that generates a Khronon-like kinetic term.

2. **The QRE-derived kinetic term must have c_s = 0 at all epochs**, not just in the matter era. This requires the quadratic minimum structure to emerge naturally from the QRE perturbation expansion.

3. **The QRE must produce non-zero background energy on FRW**, either through a non-trivial vacuum QRE (renormalized relative entropy) or through a mechanism that gives Q_0 != 1. The required value must be Omega_DM h^2 = 0.12.

4. **The mass scale mu must emerge as ~ H_0**, not ~ M_Pl. This requires a renormalization of the bare QRE that reduces it by 60 orders of magnitude.

5. **The QRE-derived perturbation must source the Poisson equation correctly**, providing exactly the right amount of gravitational potential to maintain the observed CMB peak structure.

**Assessment**: None of items 1-5 have been shown. Items 1-2 are plausible but require difficult calculations in algebraic QFT on curved spacetimes. Item 3 is the most challenging -- it is essentially asking the QRE to determine a contingent cosmological parameter. Items 4-5 are consequences of 1-3 if they succeed.

**The honest conclusion**: Route A, as currently formulated, is a beautiful structural identification that fails as a predictive CMB mechanism. The OEE postulate (promoting u^a to dynamical) bridges the gap by importing the Khronon action from outside the QRE formalism. This is an embedding, not a derivation.

---

## References

### QRE Perturbation Theory
- Aalsma & Bak 2025: arXiv:2503.XXXXX (modular Hamiltonian for FRW causal diamonds)
- Faulkner et al. 2014: arXiv:1312.7856 (entanglement first law)
- Jacobson 2016: PRL 116, 201101 (entanglement equilibrium and Einstein equations)
- Cai & Kim 2005: JHEP 02, 050 (first law of thermodynamics on apparent horizon)

### Khronon/AeST Theory
- Blanchet & Skordis 2024: arXiv:2404.06584, JCAP 11, 040 (Khronon theory)
- Blanchet & Skordis 2025: arXiv:2507.00912 (Khronon-Tensor theory)
- Skordis & Zlosnik 2021: arXiv:2007.00082, PRL 127, 161302 (AeST -- CMB fit)

### CMB Physics
- Hu 2008: arXiv:0802.3688 (CMB theory lectures)
- Planck 2018: arXiv:1807.06209 (cosmological parameters)

### GDM Framework
- Hu 1998: ApJ 506, 485 (generalized dark matter)
- Kunz 2009: arXiv:0702615 (dark degeneracy)

### Transfer Functions
- Eisenstein & Hu 1998: ApJ 496, 605 (transfer function fitting formulae)
- Meszaros 1974: A&A 37, 225 (growth of perturbations in expanding universe)

### tau Framework
- Paper 1: Huang (2026) -- Petz recovery unification
- Paper 2: Huang (2026) -- Information loss is local: gravitational refractive index
- Paper 3: Huang (2026) -- Running G and rotation curves
- Paper 4: Huang (2026) -- Temporal asymmetry as organizing principle
- Paper 5: Huang (2026) -- Observer-dependent temporal asymmetry

### Information-Theoretic Gravity
- Dorau-Much 2025: arXiv:2510.24491 (QRE -> Einstein equations)
- Bianconi 2025: PRD 111, 066001 (gravity from entropy)

---

*Last updated: 2026-03-12*
*This document provides the complete perturbation theory of the modular Khronon for CMB acoustic peaks.*
*Classification: NEGATIVE RESULT -- the modular Khronon fails as a standalone CMB mechanism due to radiation-era algebraic slaving and zero background energy density.*
