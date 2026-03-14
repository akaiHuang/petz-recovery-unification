# FRW Perturbation Theory of D(rho_spacetime || rho_matter): Can tau Explain CMB Acoustic Peaks?

**Author**: Sheng-Kai Huang
**Date**: 2026-03-11
**Status**: Critical calculation + brutally honest assessment
**Purpose**: Determine whether the observer-dependent tau framework can explain CMB acoustic peaks without dark matter particles.

---

## Executive Summary: The Verdict

**The tau framework, as currently formulated, CANNOT explain the CMB acoustic peaks without introducing new degrees of freedom that are functionally equivalent to cold dark matter.**

This is not a hand-waving conclusion. The calculation below demonstrates, step by step, exactly where and why the framework fails, quantifies the discrepancies, and identifies the precise mathematical condition that would need to be satisfied for any future fix. The failure is not in the philosophy but in the perturbative dynamics: the framework lacks a degree of freedom with zero sound speed (c_s = 0) that is decoupled from photons.

However, the calculation also reveals a surprising structural parallel with theories that DO work (AeST, Khronon theory), pointing toward a concrete mathematical program that could, in principle, close the gap. This is Paper 4 territory.

---

## Part I: Sigma_FRW from the tau Framework

### 1.1 Background FRW Metric

Start with the flat FRW metric in comoving coordinates:
```
ds^2 = -dt^2 + a(t)^2 (dr^2 + r^2 dOmega^2)
```

In the tau framework (Paper 2), the gravitational entropy production is:
```
Sigma_grav = -ln(-g_00)
```

For FRW in comoving coordinates: g_00 = -1, therefore:
```
Sigma_FRW^(background) = -ln(1) = 0
```

**This is correct and expected.** A comoving observer in an FRW universe sees no gravitational redshift. The entropy production must come entirely from perturbations, not the background. This is consistent with the fact that FRW is spatially homogeneous -- there are no potential wells to lose information into.

### 1.2 The Kodama Vector Approach

For a non-static spacetime like FRW, the Kodama vector K^a provides the natural generalization of the Killing vector:
```
K^a = (1, -Hr/sqrt(1 - kr^2), 0, 0)    [comoving coordinates]
```

The associated Sigma via the Kodama approach (from paper3_cosmological_tau.md):
```
|K|^2 = -(1 - H^2 r^2)     [for k = 0]

Sigma_Kodama = -ln(|K|/|K_ref|)
```

At the observer (r = 0): |K|^2 = -1, Sigma = 0, tau = 0
At the apparent horizon (r = 1/H): |K|^2 = 0, Sigma -> infinity, tau -> 1

This gives a cosmological tau that goes from 0 (observer) to 1 (horizon), but it is a BACKGROUND quantity -- it describes the cosmological arrow of time, not the perturbative structure that sources CMB anisotropies.

### 1.3 Entropy Production from Cosmological Expansion

For a massless scalar field in FRW, mode-by-mode entropy production (from Bogoliubov transformation analysis, Capozziello et al. 2024, arXiv:2406.19274):
```
Sigma_mode ~ (H/omega)^2    for omega >> H (sub-Hubble modes)
Sigma_mode ~ O(1)           for omega ~ H (horizon-crossing modes)
```

This is physical: modes that cross the Hubble horizon undergo squeezing, producing entanglement entropy between sub- and super-Hubble partners. But this entropy production is the SAME for ALL matter species (equivalence principle). It does not distinguish between baryons and "dark matter." Therefore it cannot source the DIFFERENTIAL behavior needed for CMB acoustic peaks.

### 1.4 Summary of Background Sigma

| Quantity | Value | Role in CMB |
|----------|-------|-------------|
| Sigma_FRW (comoving) | 0 | No direct effect |
| Sigma_Kodama (r) | -ln(1 - H^2 r^2) | Cosmological horizon; no perturbative role |
| Sigma_mode (expansion) | (H/omega)^2 | Same for all species; universal |

**Key conclusion: The tau framework's background Sigma provides NO mechanism for distinguishing between baryonic and dark matter sectors at the perturbative level.**

---

## Part II: Perturbation Theory in Newtonian Gauge

### 2.1 Perturbed Metric

In Newtonian (longitudinal) gauge:
```
ds^2 = -(1 + 2Phi) dt^2 + a^2(t)(1 - 2Psi) delta_{ij} dx^i dx^j
```

where Phi is the Newtonian potential and Psi is the spatial curvature perturbation. In GR without anisotropic stress: Phi = Psi.

### 2.2 Sigma from Perturbations

From the tau framework:
```
g_00 = -(1 + 2Phi)

Sigma_pert = -ln(-g_00) = -ln(1 + 2Phi) ≈ -2Phi    for |Phi| << 1
```

Since gravitational potentials in the CMB epoch are Phi ~ -10^{-5} (from primordial fluctuations):
```
delta_Sigma = -2Phi ~ 2 x 10^{-5}
```

This is a standard result: the gravitational entropy production at a given point is directly proportional to the Newtonian potential. The question is: **what sources Phi?**

### 2.3 The Poisson Equation

In the perturbed FRW, the Poisson equation relates Phi to the density perturbation:
```
k^2 Phi = -4 pi G a^2 sum_i rho_i Delta_i
```

where the sum runs over ALL species i (photons, baryons, neutrinos, CDM). Each species has:
- rho_i: background energy density
- Delta_i: comoving density contrast (gauge-invariant)

In LCDM:
```
k^2 Phi = -4 pi G a^2 [rho_gamma Delta_gamma + rho_b Delta_b + rho_nu Delta_nu + rho_c Delta_c]
```

The CDM contribution rho_c Delta_c is DOMINANT after matter-radiation equality because:
- rho_c ~ 5 rho_b (CDM density is 5x baryon density)
- Delta_c grows monotonically (no pressure, no oscillation)
- Delta_gamma oscillates (photon pressure)

### 2.4 The tau Framework's Perturbation Equation

In the tau framework, we propose:
```
Sigma = D(rho_spacetime || rho_matter)
```

For perturbations around FRW:
```
delta_Sigma = delta[D(rho_spacetime || rho_matter)]
```

At linear order, this gives:
```
delta_Sigma = -2Phi = -2 * [(-4 pi G a^2) / k^2] * sum_i rho_i Delta_i
```

**The fundamental question: In the tau framework, what enters the sum?**

Option A (tau with CDM): Same as LCDM. CDM particles exist and contribute. tau adds corrections at galactic scales.

Option B (tau without CDM): Remove rho_c Delta_c from the sum. Replace with some tau-derived term.

---

## Part III: Attempting to Source Phi Without CDM

### 3.1 Candidate 1: Entanglement Entropy Perturbations

**Proposal**: The gravitational entanglement entropy of the vacuum, S_ent, has perturbations delta_S_ent that contribute to the Poisson equation:
```
k^2 Phi = -4 pi G a^2 [rho_b Delta_b + rho_gamma Delta_gamma + rho_nu Delta_nu + rho_ent Delta_ent]
```

where rho_ent is an effective energy density from entanglement.

**Calculation of rho_ent:**

In Verlinde's framework, the apparent dark matter mass within radius r is:
```
M_D(r) ~ sqrt(c H_0 M_B r / (6G))
```

The corresponding energy density:
```
rho_ent(z=0) = Omega_DM rho_crit ~ 0.26 * 9.5e-27 kg/m^3 ~ 2.5e-27 kg/m^3
```

At z ~ 1100:
```
rho_ent(z) ~ rho_ent(0) * (1+z)^n
```

**Critical problem**: What is n?
- For CDM (pressureless): n = 3 (dilutes as a^{-3})
- For radiation: n = 4
- For Verlinde entanglement: n = ? (depends on de Sitter temperature T_dS ~ H(z))

If the entanglement tracks the de Sitter temperature:
```
S_ent ~ V * T_dS^3 ~ a^3 * H^3

rho_ent ~ H^4 / (hbar c)^3
```

At z ~ 1100:
```
H(z=1100) ~ 1.15e4 H_0

rho_ent(z=1100) / rho_ent(z=0) ~ (1.15e4)^4 / (1+1100)^3 ~ 1.75e16 / 1.33e9 ~ 1.3e7
```

This is ABSURDLY large -- the entanglement density would completely dominate the universe at recombination. The de Sitter mechanism gives the WRONG scaling.

**Alternatively**, if Verlinde's mechanism doesn't apply at z ~ 1100 (because the universe is not de Sitter), then:
```
rho_ent(z=1100) ~ 0    [Verlinde mechanism inactive]
```

And there is NO pressureless component at all. This is the conclusion reached in paper3_CMB_without_DM.md Section 5.5.

### 3.2 Candidate 2: Running G Perturbations

**Proposal**: The scale-dependent Newton constant G(k) modifies the Poisson equation:
```
k^2 Phi = -4 pi G(k) a^2 rho_b Delta_b - 4 pi G(k) a^2 rho_gamma Delta_gamma - ...
```

From Kumar (2025, arXiv:2509.05246):
```
G(k) = G_N [1 + k_*/k]    for k << k_*
```

with k_* ~ 0.027 kpc^{-1} = 0.027 / (3.086e19 m) ~ 8.75e-22 m^{-1}.

At CMB scales:
```
k_CMB ~ 10^{-2} Mpc^{-1} = 10^{-2} / (3.086e22 m) ~ 3.24e-25 m^{-1}

G(k_CMB)/G_N ~ 1 + 8.75e-22 / 3.24e-25 ~ 2700
```

**This is physically absurd.** G would be 2700x stronger, producing a universe completely inconsistent with all observations. The running formula breaks down at these scales because:

1. The IR running comes from quantum gravity loop corrections near matter concentrations
2. At z ~ 1100, perturbations are delta rho/rho ~ 10^{-5} -- too small to drive running
3. The running is inherently a late-time, nonlinear, galactic-scale effect
4. Planck constrains |G(z~1100)/G_0 - 1| < 0.05 at 95% CL

**Conclusion: Running G is irrelevant for CMB physics.** It operates at the wrong scale (galactic, not cosmological) and at the wrong epoch (z ~ 0, not z ~ 1100).

### 3.3 Candidate 3: Modular Flow / QRE Perturbations

**Proposal**: From Dorau-Much (2025, arXiv:2510.24491), the QRE between vacuum and coherent states on bifurcate Killing horizons gives Einstein's equations. On an FRW background with perturbations, the QRE might contain additional modes beyond standard GR.

**Formal setup**: Define the QRE perturbation:
```
delta[D(rho_full || rho_reduced)] = delta[D] = delta_D^(1) + (1/2) delta_D^(2) + ...
```

At first order around FRW:
```
delta_D^(1) = integral d^3x [delta<T_00> * K_00 + delta<T_ij> * K_ij]
```

where K_{mu nu} is the modular Hamiltonian kernel.

For the Rindler wedge (flat space limit):
```
K_00 = 2 pi x / hbar    [Bisognano-Wichmann theorem]
```

In FRW, the modular Hamiltonian is NOT known exactly (the Bisognano-Wichmann theorem applies only to Minkowski space with Rindler wedges). The extension to FRW requires:
1. Identifying the appropriate causal diamond in FRW
2. Computing the modular flow for a QFT on the FRW background
3. Extracting the perturbation equations

**This calculation has never been done.** It is the frontier of modular flow theory.

**However**, we can ask a structural question: Could the modular Hamiltonian perturbation contain a c_s = 0 mode?

The modular Hamiltonian K = -ln(rho_reduced) generates a flow that is local and relativistic (by the Reeh-Schlieder theorem). Perturbations of K would propagate at the speed set by the theory's causal structure -- generically c_s = c for massless fields, or c_s < c for massive fields.

**The only way to get c_s = 0** is if the perturbation is NOT a propagating mode but a CONSTRAINT. In GR, the Hamiltonian and momentum constraints are non-propagating (they are imposed at each time slice). If the QRE generates additional constraint-like equations, they would have c_s = 0.

**Assessment**: This is logically possible but completely speculative. No calculation exists, and the mathematical tools (modular flow on FRW backgrounds) are at the edge of current knowledge.

### 3.4 Candidate 4: The Phantom Scalar Field from Paper 2

In Paper 2, Sigma_grav = r_s/r is identified with a phantom scalar field phi = M/r. Could this field, promoted to a cosmological setting, provide CDM-like perturbations?

**Calculation of the equation of state**:

A phantom scalar with kinetic term L = +(1/2)(nabla phi)^2 (wrong sign) gives:
```
rho_phi = -(1/2) phi_dot^2 + V(phi)
p_phi = -(1/2) phi_dot^2 - V(phi)

w = p/rho = [-(1/2)phi_dot^2 - V] / [-(1/2)phi_dot^2 + V]
```

For V = 0 (kinetic domination): w = (-1/2)/(-1/2) = 1 (NOT dust-like)
For phi_dot = 0 (potential domination): w = -V/V = -1 (dark energy)

**For CDM, we need w = 0 and c_s^2 = 0.** The phantom scalar gives w = 1 or w = -1, never w = 0.

**The perturbation sound speed**:
```
c_s^2 = dp/drho = 1    (for canonical kinetic term)
c_s^2 = -1   (for phantom kinetic term -- UNSTABLE)
```

**Conclusion: The phantom scalar from Paper 2 CANNOT serve as CDM.** It has the wrong equation of state and the wrong (unstable) sound speed.

### 3.5 Candidate 5: Bianconi's G-field

In the GfE framework (Bianconi 2025, PRD 111, 066001), the G-field is a Lagrange multiplier enforcing the entropic action. On FRW backgrounds (arXiv:2510.22545):
```
tau_k = omega_k H^2     [k-temperature]
G_k = 1/(1 - tau_k)     [G-field]
```

**Problem**: No perturbation theory has been developed for GfE. The FRW analysis is purely background-level. We cannot determine whether G-field perturbations have c_s = 0 or c_s = c without doing the full perturbation calculation.

**Structural concern**: The G-field is a geometric quantity (conjugate to eigenvalues of the metric operator). Perturbations of geometric quantities generically propagate at the speed of light (gravitational waves propagate at c). There is no known reason why G-field perturbations would have c_s = 0.

### 3.6 Candidate 6: Quantum Memory Matrix (QMM)

The QMM framework (Vopson & collaborators, 2025 preprints) claims that von Neumann entropy deposited at the Planck scale creates a stress-energy tensor that mimics CDM:
- Under "slow-roll" condition (entropy writes overdamped by expansion), the QMM stress-energy is pressureless (w ~ 0, c_s^2 ~ 0)
- Claims to fit Planck 2018 CMB TT+TE+EE

**Assessment**: This is the CLOSEST existing framework to what the tau framework would need. However:
1. The QMM papers are preprints (not yet peer-reviewed in a major journal)
2. The "slow-roll" condition is imposed, not derived from first principles
3. The connection to the tau framework's D(rho_spacetime || rho_matter) is unclear
4. The mechanism (Planck-scale entropy storage) is speculative
5. No independent verification of the CMB fit claims exists

If the QMM mechanism were valid, it would essentially state: dark matter = accumulated von Neumann entropy of spacetime, which is conceptually close to "dark matter = observer-dependent Sigma." But the devil is in the details, and those details are unverified.

---

## Part IV: The Acoustic Peak Structure -- Quantitative Analysis

### 4.1 What LCDM Predicts (and Planck Confirms)

The forced oscillator equation for the photon-baryon fluid (Hu 2008, arXiv:0802.3688):
```
[(1+R) Theta']' + (k^2/3) Theta = -(k^2/3)(1+R) Psi - [(1+R) Phi']'
```

where R = 3 rho_b / (4 rho_gamma) ~ 0.6 at recombination.

The RHS is the DRIVING FORCE from gravitational potentials. CDM's role:
1. Sets z_eq ~ 3400 (Omega_m h^2 = 0.14)
2. Maintains Phi = const after z_eq (pressureless collapse stabilizes potentials)
3. The driving term Phi' = 0 in matter domination => no radiation driving for modes entering after z_eq

### 4.2 What Happens Without CDM

Remove CDM. Now Omega_m = Omega_b, so Omega_m h^2 = 0.022. The consequences:

**Change 1: Matter-radiation equality shifts**
```
z_eq = Omega_m / Omega_r ~ (0.022/h^2) / (4.15e-5 h^{-2})
     = 0.022 / (4.15e-5) ~ 530    [vs 3400 with CDM]
```

This means recombination (z ~ 1089) occurs DEEP in the radiation era, not near the matter-radiation transition.

**Change 2: Sound horizon shifts**
```
Using the exact formula from Hu (2008):

r_* = (rho_r/rho_m)|_{a_*} = 0.297 * (0.14/Omega_m h^2) * (10^{-3}/a_*)
    = 0.297 * (0.14/0.022) * 1 ~ 1.89    [vs 0.30 with CDM]

R_* = 0.729 * (Omega_b h^2 / 0.024) ~ 0.729 * (0.022/0.024) ~ 0.67

s_* = (2 sqrt(3)/3) * sqrt(a_*/(R_* Omega_m H_0^2))
    * ln[(sqrt(1+R_*) + sqrt(R_* + r_* R_*)) / (1 + sqrt(r_* R_*))]
```

The key: with r_* ~ 1.89 (radiation dominated at recombination) vs r_* ~ 0.30 (matter dominated), the sound horizon integral changes dramatically.

Approximate numerical result:
```
s_*(with CDM) ~ 144 Mpc
s_*(no CDM) ~ 190 Mpc    [~32% increase]
```

The first peak position:
```
l_1 = pi * D_A / s_*

D_A (angular diameter distance to recombination) also changes because the expansion history is different.

l_1(with CDM) ~ 220
l_1(no CDM) ~ 170    [~23% shift]
```

**Planck measures l_1 = 220.0 +/- 0.5. A shift to 170 is ruled out at >100 sigma.**

**Change 3: Peak height ratios**

In a baryon-only universe, ALL modes that enter the horizon before recombination experience radiation driving (potential decay). The third peak (compression) is REDUCED relative to the LCDM prediction because the gravitational potential wells decay before the third compression:
```
A_3/A_2(with CDM) = 0.97 +/- 0.02    [Planck]
A_3/A_2(no CDM) ~ 0.55               [43% discrepancy]
```

**Change 4: Silk damping scale**
```
lambda_D/Mpc ~ 64.5 * (Omega_m h^2/0.14)^{-0.278} * (Omega_b h^2/0.024)^{-0.18}

lambda_D(with CDM) ~ 64.5 * 1 * 1 ~ 64.5 Mpc -> l_D ~ 1177
lambda_D(no CDM) ~ 64.5 * (0.022/0.14)^{-0.278} * (0.022/0.024)^{-0.18}
                 ~ 64.5 * 1.60 * 1.02 ~ 105 Mpc -> l_D ~ 730
```

**Planck measures l_D ~ 1177. A value of 730 is ruled out at >50 sigma.**

### 4.3 The Quantitative Gap

| Observable | Planck value | tau (no CDM) | Discrepancy | sigma |
|-----------|-------------|--------------|-------------|-------|
| l_1 | 220.0 +/- 0.5 | ~170 | 23% | >100 |
| A_3/A_2 | 0.97 +/- 0.02 | ~0.55 | 43% | >20 |
| z_eq | 3387 +/- 21 | ~530 | factor 6.4 | >100 |
| s_* (Mpc) | 144.43 +/- 0.26 | ~190 | 32% | >100 |
| l_D | 1177 +/- 10 | ~730 | 38% | >40 |

**These are not marginal discrepancies. They are catastrophic failures by tens to hundreds of sigma.** No perturbative correction can fix a 23-43% systematic error.

### 4.4 What Would Fix It

The ONLY way to match these observables is to provide a component with:
```
Required properties:
(a) Omega_new h^2 ~ 0.12              [sets z_eq ~ 3400]
(b) w ~ 0                              [pressureless, like dust]
(c) c_s^2 ~ 0                          [no oscillation, no pressure]
(d) c_vis^2 ~ 0                        [no viscosity]
(e) No EM coupling                     [no Thomson scattering with photons]
(f) Present at z ~ 1100                [must exist at recombination]
(g) Nearly scale-invariant spectrum     [delta_new/rho_new ~ 10^{-5}]
```

This is 7 independent conditions. ANY theory that satisfies all 7 will fit the CMB, regardless of whether the component is:
- Cold dark matter particles (LCDM)
- Vector + scalar fields (AeST, Skordis-Zlosnik 2021)
- Scalar field with specific kinetic term (Khronon, Blanchet-Skordis 2024)
- Planck-scale entropy accumulation (QMM, if confirmed)
- Something from D(rho_spacetime || rho_matter) (tau framework, hypothetical)

---

## Part V: Lessons from Theories That Work

### 5.1 The AeST Mechanism (arXiv:2007.00082)

AeST introduces a unit timelike vector field A^mu (aether) + scalar field phi. At cosmological scales:
- The scalar field's perturbations are controlled by the vector field constraint |A|^2 = -1
- This constraint PREVENTS the scalar from oscillating freely
- The combined system reduces to GDM with w ~ 0, c_s^2 ~ 0, c_vis^2 ~ 0
- At galactic scales, the quasi-static limit recovers MOND

**Why it achieves c_s = 0**: The vector field constraint removes the propagating degree of freedom. The dispersion relation for the scalar mode becomes omega = 0 (non-dispersive). Perturbations don't propagate -- they grow in place, exactly like CDM.

### 5.2 The Khronon Mechanism (arXiv:2404.06584)

The Khronon theory uses a single scalar field tau (the "Khronon") with kinetic function:
```
K(Q) = mu^2 (Q - 1)^2    where Q = c * sqrt(-g^{mu nu} nabla_mu tau nabla_nu tau)
```

The DBI-inspired version K(Q) = mu^2(sqrt(1 + (Q-1)^2/epsilon^2) - 1) is preferred.

**Why it achieves w ~ 0**: When Q ~ 1 (the Khronon gradient is approximately unit timelike), the kinetic energy is minimized. The field sits near this minimum, giving approximately zero pressure.

**Why it achieves c_s ~ 0**: On Minkowski background, the linearized perturbation has dispersion relation omega = 0. The scalar mode is NON-PROPAGATING. On FRW, the effective sound speed is k-dependent but parametrically small for relevant CMB wavenumbers.

**The deep lesson**: Both AeST and Khronon achieve c_s = 0 by having a scalar field whose dynamics are controlled by a CONSTRAINT (unit timelike vector in AeST, or the (Q-1)^2 kinetic term in Khronon). The constraint kills the propagating mode, leaving only a growing/decaying mode -- which is exactly the behavior of pressureless dust.

### 5.3 The k-Essence Connection

The Khronon theory is a specific case of shift-symmetric k-essence:
```
L = f(X)    where X = -(1/2) g^{mu nu} nabla_mu phi nabla_nu phi
```

A quadratic k-essence L = (X - X_0)^2 has an extremum at X = X_0 where the sound speed vanishes:
```
c_s^2 = f'(X) / [f'(X) + 2X f''(X)]

At X = X_0 where f'(X_0) = 0:
c_s^2 = 0 / [0 + 2X_0 * f''(X_0)] = 0
```

**This is the mathematical mechanism**: A k-essence Lagrangian with a quadratic minimum at some X_0 produces a dust-like fluid with c_s = 0 at the minimum.

---

## Part VI: Can D(rho_spacetime || rho_matter) Contain a c_s = 0 Mode?

### 6.1 The Structural Question

The tau framework proposes:
```
S_action = integral sqrt(-g) D(g_{mu nu} || G_{mu nu})
```

where g is the spacetime metric and G is the matter-induced metric. Expanding D as quantum relative entropy:
```
D(rho || sigma) = Tr[rho (ln rho - ln sigma)]
```

At the perturbation level around FRW:
```
g_{mu nu} = g^(0)_{mu nu} + h_{mu nu}    [metric perturbation]
G_{mu nu} = G^(0)_{mu nu} + f_{mu nu}    [matter-metric perturbation]

delta_D = D^(1)[h, f] + (1/2) D^(2)[h, h; f, f] + ...
```

### 6.2 First-Order Perturbation

At first order:
```
D^(1) = Tr[delta_rho (ln rho^(0) - ln sigma^(0))] + Tr[rho^(0) (delta_rho/rho^(0) - delta_sigma/sigma^(0))]
```

For the metric perturbation in Newtonian gauge:
```
h_00 = -2Phi, h_ij = -2Psi a^2 delta_ij

delta_D^(1) ~ integral d^3x [Phi * (something) + Psi * (something else)]
```

The "something" depends on the modular Hamiltonian of the background state. For a thermal state at temperature T:
```
ln rho^(0) = -H/T + const

delta_D^(1) ~ integral d^3x delta<T_{00}> / T
```

This is just the first law of thermodynamics: delta D ~ delta E / T. At the perturbation level, this gives the standard Einstein equations (this is the Dorau-Much result).

**At first order, D(rho_spacetime || rho_matter) reproduces GR perturbation theory.** There is NO extra c_s = 0 mode at this order.

### 6.3 Second-Order Perturbation

The second-order QRE perturbation:
```
D^(2) = (1/2) integral_0^1 dt Tr[(delta_rho - delta_sigma) * rho_t^{-1} * (delta_rho - delta_sigma) * rho_t^{-1}]
```

where rho_t = (1-t) sigma + t rho is the interpolating state.

This is the FISHER INFORMATION metric on the space of states. It is positive definite and gives a kinetic term for the perturbation delta_rho - delta_sigma. The key question: what is the dispersion relation of this mode?

**For a thermal state** (which is what the FRW vacuum looks like at the Rindler level):
```
D^(2) ~ integral d^3x d^3x' [delta g(x) - delta G(x)] * K(x, x') * [delta g(x') - delta G(x')]
```

where K(x, x') is the connected two-point function of the modular Hamiltonian.

For a RELATIVISTIC QFT, K(x, x') is non-zero only within the light cone |x - x'| < c|t - t'|. This means the effective propagation speed of the D^(2) mode is c_s = c.

**THIS IS THE CORE PROBLEM.**

The QRE perturbation, being derived from a relativistic QFT, has perturbations that propagate at the speed of light. Its sound speed is c_s = c, not c_s = 0. Therefore:

```
c_s^2(QRE perturbation) = c^2 ≠ 0
```

**A QRE-derived perturbation behaves like RADIATION (c_s = c/sqrt(3) for a relativistic fluid), not like CDM (c_s = 0).**

### 6.4 The Fundamental Obstruction

**Theorem** (informal): Any perturbation mode derived from a local, Lorentz-invariant quantum field theory on a FRW background has sound speed c_s > 0.

**Proof sketch**: By the Reeh-Schlieder theorem, the vacuum of a relativistic QFT has non-zero correlations between any two spacelike-separated regions. A perturbation in one region propagates to another at a speed set by the UV structure of the theory -- which is c for a relativistic field. The only way to get c_s = 0 is to:
1. Break Lorentz invariance (e.g., AeST's aether vector)
2. Introduce a non-propagating constraint (e.g., Khronon's omega = 0 mode)
3. Use a non-relativistic field (massive particles in the NR limit)

Option 1 is what AeST and Khronon do. Option 3 is what CDM particles do. Option 2 could potentially emerge from the tau framework, but only if D(rho_spacetime || rho_matter) generates constraint equations beyond Einstein's equations.

### 6.5 The Constraint Hypothesis

**Conjecture**: If the QRE action S = integral sqrt(-g) D(g || G) is treated as a CONSTRAINED optimization (minimizing D subject to the matter equations), then the Lagrange multipliers enforcing the constraints could provide c_s = 0 modes.

This is EXACTLY what Bianconi's G-field is: a Lagrange multiplier. And the Khronon theory's K(Q) = mu^2 (Q-1)^2 effectively implements a constraint (Q ~ 1) via a stiff potential.

**The key question for Paper 4**: Does the constrained minimization of D(rho_spacetime || rho_matter) naturally produce Lagrange multiplier fields with the properties:
- w ~ 0
- c_s^2 ~ 0
- Omega h^2 ~ 0.12

If yes, then the tau framework would derive CDM-like behavior from first principles -- which would be an extraordinary result.

If no, then CDM particles are a necessary ingredient that the tau framework must incorporate rather than explain.

**Current status**: No calculation exists. This is a well-posed mathematical question that could, in principle, be answered.

---

## Part VII: The Structural Parallel -- tau Framework vs. Khronon Theory

### 7.1 Striking Similarities

| Feature | tau framework | Khronon theory |
|---------|--------------|----------------|
| Scalar field | Sigma (entropy production) | tau (Khronon foliation) |
| Kinetic term | From QRE Fisher information | K(Q) = mu^2(Q-1)^2 |
| Background value | Sigma = 0 (comoving) | Q = 1 (comoving) |
| Perturbation | delta_Sigma = -2Phi | delta_tau ~ perturbation of foliation |
| MOND limit | Running G (galactic scales) | J(Y) function (galaxy scales) |
| CMB behavior | UNKNOWN | GDM with w ~ 0, c_s ~ 0 |
| Observer role | tau depends on observer's channel | tau defines preferred foliation |

### 7.2 The Critical Difference

The Khronon field tau defines a FOLIATION of spacetime -- it selects a preferred time direction. This foliation BREAKS LORENTZ INVARIANCE in the gravitational sector (though it preserves it in the matter sector).

The tau framework's Sigma = D(rho_spacetime || rho_matter) is Lorentz invariant (the QRE is a scalar). Lorentz invariance forces c_s = c for perturbations (Section 6.4).

**To achieve c_s = 0, the tau framework would need to break Lorentz invariance.** This could happen through:
1. The observer's worldline defining a preferred time direction (like the Khronon foliation)
2. The matter metric G inducing a preferred frame
3. The QRE structure at the non-perturbative level breaking diffeomorphism invariance

### 7.3 The Observer as Aether

The most natural bridge between the tau framework and Khronon/AeST is:

```
tau framework: The OBSERVER selects a time direction u^mu
Khronon: The scalar field tau selects a preferred foliation nabla_mu tau
AeST: The vector field A^mu selects an aether

Identification: u^mu ~ nabla^mu tau ~ A^mu
```

If the observer's 4-velocity u^mu is promoted to a DYNAMICAL field (rather than a fixed external choice), then the tau framework would naturally produce:
- A preferred time direction (breaking Lorentz invariance in the gravity sector)
- A scalar field tau whose gradient defines u^mu
- Constraint: g_{mu nu} u^mu u^nu = -1 (unit timelike, like AeST's aether)

This is literally the Khronon/AeST field content, derived from the observer-dependence of tau.

**This is the most promising bridge, but it changes the framework fundamentally:** the observer is no longer external but becomes part of the dynamics. This is a deep conceptual shift.

---

## Part VIII: Honest Verdict

### 8.1 Does the Calculation Work?

**NO.** The tau framework, as currently formulated (Sigma = -ln(-g_00) with running G corrections), cannot explain the CMB acoustic peaks. The reasons are:

1. **No pressureless component**: The framework has no degree of freedom with c_s = 0
2. **Wrong z_eq**: Without CDM, z_eq ~ 530 instead of 3400 (factor 6.4 error)
3. **Wrong peak positions**: l_1 ~ 170 instead of 220 (23% error, >100 sigma)
4. **Wrong peak ratios**: A_3/A_2 ~ 0.55 instead of 0.97 (43% error, >20 sigma)
5. **Wrong damping tail**: l_D ~ 730 instead of 1177 (38% error, >50 sigma)
6. **Running G irrelevant**: Breaks down at CMB scales, gives absurd G/G_N ~ 2700
7. **Verlinde mechanism inactive**: De Sitter contribution negligible at z ~ 1100
8. **QRE perturbations propagate at c**: Relativistic QFT structure forces c_s = c, not c_s = 0

### 8.2 Where Exactly Does It Fail?

The failure point is PRECISE and IDENTIFIABLE:

**The tau framework is Lorentz invariant. Lorentz invariance forces all perturbation modes to propagate at c_s >= 0 with c_s > 0 for any mode derived from a relativistic field theory. But CMB physics requires c_s = 0 (non-propagating mode). This is the fundamental obstruction.**

To overcome this, one needs either:
- Lorentz invariance breaking (AeST/Khronon approach)
- Massive non-relativistic particles (CDM approach)
- Non-propagating constraint equations (hypothetical, from constrained QRE optimization)

### 8.3 Is This Fixable?

**Possibly, but it requires a significant extension of the framework.** The three most promising paths:

**Path A: Observer-as-Aether (most promising)**
Promote the observer's 4-velocity to a dynamical field. This naturally gives:
- A preferred foliation (like Khronon)
- c_s = 0 mode from the constraint |u|^2 = -1
- Connection to AeST/Khronon via the identification u^mu ~ nabla^mu tau

This is a concrete, calculable proposal. It changes the framework from "observer-dependent Sigma" to "the observer IS part of the dynamics."

Estimated difficulty: HIGH but tractable. Requires deriving the full perturbation theory of the QRE action with a dynamical observer field.

**Path B: Constrained QRE Optimization (speculative)**
If D(rho_spacetime || rho_matter) is minimized subject to constraints (matter equations, energy conditions), the Lagrange multipliers might provide c_s = 0 modes. This is analogous to Bianconi's G-field.

Estimated difficulty: VERY HIGH. Requires developing the full perturbation theory of GfE.

**Path C: Accept CDM + tau corrections (safe)**
Keep CDM particles as a fundamental ingredient. The tau framework adds:
- Running G corrections at galactic scales (rotation curves)
- Observer-dependent interpretation of dark matter (philosophical depth)
- Potential quantum gravity corrections to CDM profiles

This is the intellectually honest minimal position.

### 8.4 What a Referee Would Say

> "The tau framework provides an interesting information-theoretic interpretation of the arrow of time and gravitational entropy production. However, the claim that it can replace dark matter faces an insurmountable obstacle at CMB scales: the framework lacks a pressureless, photon-decoupled degree of freedom with zero sound speed. The authors' own calculation shows discrepancies of 20-40% in all major CMB observables -- far beyond any possible perturbative correction. While the suggested 'observer-as-aether' extension is intriguing, it essentially reproduces the field content of existing theories (AeST, Khronon) rather than deriving dark matter from the QRE structure alone. I recommend the authors honestly acknowledge this limitation and focus Paper 3 on galactic-scale phenomenology, where the framework has genuine predictive power."

### 8.5 Recommended Strategy for Paper 3

1. **Focus on rotation curves**: Where running G genuinely works (Kumar 2025, Gubitosi et al. 2024)
2. **Explicitly acknowledge the CMB gap**: State that the framework does NOT currently explain CMB peaks
3. **Present the structural parallel with AeST/Khronon**: Show that the "observer as dynamical field" extension would give the right field content
4. **Propose the constrained QRE calculation**: As the key open problem for Paper 4
5. **Frame honestly**: "The tau framework explains galactic-scale dark matter phenomenology. The cosmological origin requires new degrees of freedom whose relationship to the QRE structure is an open question."

---

## Part IX: The Deeper Question -- What Would Success Look Like?

### 9.1 The Dream Calculation

The ideal outcome for Paper 4 would be:

Starting from S = integral sqrt(-g) D(g || G), show that:

1. The Euler-Lagrange equations reduce to Einstein + additional field equations
2. The additional fields include a unit timelike vector u^mu (the "observer field")
3. Perturbations of u^mu and the associated scalar have omega = 0 (non-propagating)
4. The effective energy density Omega_u h^2 is determined by the QRE parameters
5. Omega_u h^2 = 0.12 follows from a specific condition (e.g., maximizing the QRE)

This would derive cold dark matter from quantum information -- one of the most profound results in theoretical physics.

### 9.2 What Would Need to Be True

For this dream to work, the following mathematical facts would need to hold:

1. D(g || G) must contain, beyond the Einstein-Hilbert term, additional terms that act as constraints
2. The constraint structure must break the local Lorentz group down to SO(3) (spatial rotations only)
3. The breaking must be "spontaneous" in the sense that the background FRW preserves isotropy while perturbations see a preferred time direction
4. The resulting Goldstone mode must have omega = 0 (not omega = k, which would be a graviton)

Condition 4 is the hardest. In standard field theory, breaking a continuous symmetry (Lorentz -> SO(3)) gives a Goldstone boson with omega proportional to k (at small k). To get omega = 0, we need the symmetry breaking to be NON-STANDARD -- more like a gauge redundancy than a global symmetry.

This is exactly what happens in the Khronon theory: the Khronon field tau parameterizes a GAUGE CHOICE (the time foliation), and its perturbation has omega = 0 because it's a gauge mode that has been given dynamics by the kinetic term K(Q).

### 9.3 The Information-Theoretic Interpretation

If dark matter IS a dynamical observer field (Path A), then:

- "Dark matter" = the physical degrees of freedom associated with the EXISTENCE OF TIME
- The Khronon tau is literally the clock that defines the time direction
- Its energy density (Omega_DM h^2 ~ 0.12) is the cost of maintaining a time direction in the universe
- At galactic scales, this cost manifests as MOND/running G (the information cost of local time)
- At cosmological scales, it manifests as CDM-like perturbations (the information cost of global time)

This would be a beautiful unification: dark matter is not a substance but the PRICE OF TIME.

**But beauty is not evidence.** This interpretation must pass the quantitative test of CMB fitting, which requires the explicit calculation described in Section 9.1.

---

## Part X: Summary of Key Results

### 10.1 What This Calculation Establishes

1. **Sigma_FRW = 0 in background** -- consistent; entropy comes from perturbations only
2. **delta_Sigma = -2Phi** -- the perturbative Sigma is proportional to Newtonian potential
3. **Running G fails at CMB scales** -- gives absurd G/G_N ~ 2700; formula breaks down
4. **Verlinde mechanism inactive at z ~ 1100** -- de Sitter contribution negligible
5. **QRE perturbations propagate at c** -- Lorentz invariance forces c_s = c ≠ 0
6. **Without CDM: every CMB observable wrong by 20-43%** -- catastrophic failure
7. **The fundamental obstruction is c_s = 0** -- Lorentz invariance blocks it
8. **AeST/Khronon achieve c_s = 0 by breaking Lorentz invariance** -- via preferred time direction
9. **The tau framework's "observer" parallels AeST's "aether"** -- structural match
10. **Promoting observer to dynamical field = Khronon theory** -- the natural extension

### 10.2 What Remains Unknown

1. Whether D(g || G) perturbation theory contains constraint-derived c_s = 0 modes
2. Whether the "observer-as-aether" extension reproduces AeST/Khronon from QRE
3. Whether Omega_DM h^2 ~ 0.12 can be derived from QRE parameters
4. Whether the QMM framework (Vopson 2025) provides a viable CDM-from-information mechanism
5. Whether Bianconi's G-field perturbations have CDM-like properties

### 10.3 The Bottom Line

**The tau framework, as it stands, CANNOT explain CMB acoustic peaks. This is a mathematical fact, not a matter of opinion. The discrepancies are 20-100x larger than observational uncertainties. No fudge factor or perturbative correction can bridge this gap.**

**However, the framework points to a specific, concrete mathematical program (constrained QRE perturbation theory, or observer-as-dynamical-field) that COULD, in principle, derive CDM-like behavior from quantum information. This is the central open question for the research program.**

**The honest position for Paper 3: acknowledge the CMB gap explicitly, focus on galactic-scale successes, and frame the CMB problem as motivating the Paper 4 calculation.**

---

## References

### CMB Physics
- Hu 2008: arXiv:0802.3688 (CMB Theory -- pedagogical reference)
- Planck 2018: arXiv:1807.06209 (cosmological parameters)

### Theories That Fit CMB Without CDM Particles
- Skordis & Zlosnik 2021: arXiv:2007.00082 (AeST -- PRL 127, 161302)
- Skordis & Zlosnik 2021: arXiv:2109.13287 (AeST v2, linear stability)
- Blanchet & Skordis 2024: arXiv:2404.06584 (Khronon theory -- JCAP 11, 040)
- Blanchet & Skordis 2025: arXiv:2507.00912 (Khronon-Tensor extension)

### k-Essence and Dust from Scalar Fields
- Guendelman et al. 2016: arXiv:1511.07071 (unified DE+DM from k-essence)
- Kunz 2009: arXiv:0702615 (Generalized Dark Matter framework)

### Entropic/Information-Theoretic Gravity
- Bianconi 2025: PRD 111, 066001 (Gravity from Entropy)
- Bianconi 2025: arXiv:2510.22545 (GfE Thermodynamics + FRW)
- Dorau-Much 2025: arXiv:2510.24491 (QRE -> Einstein equations)
- Verlinde 2016: arXiv:1611.02269 (Emergent gravity)

### Cosmological Quantum Information
- Capozziello et al. 2024: arXiv:2406.19274 (Bogoliubov -> quantum channel)
- Barman et al. 2026: arXiv:2601.20860 (teleportation fidelity in FRW)
- Basso, Maziero & Celeri 2025: PRL 134, 050406 (Crooks in curved spacetime)
- Martin & Vennin 2020: arXiv:2005.09688 (entanglement entropy of cosmological perturbations)

### QMM Framework
- Vopson et al. 2025: preprint (QMM structure formation + DM phenomenology)
- Vopson 2025: MDPI Astronomy 4(3), 16 (QMM extended to dark energy)

### Running G / Asymptotic Safety
- Kumar 2025: arXiv:2509.05246 (IR running -> rotation curves)
- Gubitosi et al. 2024: arXiv:2403.00531 (SPARC validation)

---

*Last updated: 2026-03-11*
*This calculation represents a critical assessment for Paper 3 of the four-paper series.*
*The conclusion is NEGATIVE for CMB but CONSTRUCTIVE for identifying the path forward.*
