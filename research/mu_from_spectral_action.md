# Attempt to Derive mu = (Sigma m_nu)^2 / M_Pl from the Spectral Action

**Author**: Sheng-Kai Huang (with systematic derivation analysis)
**Date**: 2026-03-19
**Status**: FIVE ROUTES EXPLORED -- one promising self-consistency argument, no clean first-principles derivation yet
**Purpose**: Determine whether the gravitational seesaw mu = (Sigma m_nu)^2 / M_Pl can be derived from the spectral action framework, building on the Paper 9 result g_QQ = (d-2)(d-3)/Q^2 = 2/Q^2 in d=4.

---

## Executive Summary

| Route | Method | Result | Verdict |
|-------|--------|--------|---------|
| 1 | Direct from D_F eigenvalues | mu ~ 10^{-31} eV not an eigenvalue of D_F | FAILS |
| 2 | One-loop dilaton mass (Coleman-Weinberg) | mu ~ m_nu^4/(8pi^2 M_Pl^2) ~ 10^{-32} eV | OFF by ~10x |
| 3 | Non-perturbative instanton | Lambda_NP = (mu M_Pl^3)^{1/4} ~ 27 TeV (unidentifiable) | NO MECHANISM |
| 4 | Ghost condensation self-consistency | M_GC = Sigma m_nu as IR cutoff -> only neutrinos contribute -> self-consistent | PROMISING but circular |
| 5 | Spectral action a_2 Fisher metric | g_QQ = 2/Q^2 constrains FORM of K(Q) but not mu | FORM MATCHES, mu undetermined |

**Bottom line**: The spectral action beautifully explains WHY K(Q) has the ghost condensation form (via the Fisher metric of the a_2 coefficient), and the self-consistency argument explains WHY neutrinos are singled out (they are the only fermions below M_GC). But the actual VALUE of mu = (Sigma m_nu)^2/M_Pl remains a numerical observation rather than a theorem. The closest approach (Route 4) is self-consistent but logically circular.

---

## 0. Known Inputs from Previous Analysis

From `paper9_d4_calculation.md`:
- g_QQ = (d-2)(d-3)/Q^2 is the Fisher metric of the a_2 (Einstein-Hilbert) sector
- In d=4: g_QQ = 2/Q^2, matching the channel result Sigma = 2 ln Q
- The a_4 sector (Yang-Mills, Higgs quartic) is conformally invariant, contributes g_QQ = 0
- This SELECTS d=4 as the unique dimension where g_QQ^{EH} = 2/Q^2

From `neutrino_mu_connection.md`:
- M_GC = sqrt(mu M_Pl) = 59.17 meV (with full Planck mass)
- Sigma m_nu (NH, min) = 58.21 meV (NuFIT 5.2), 58.68 meV (NuFIT 5.3)
- Discrepancy: 1.66% (NuFIT 5.2), 0.84% (NuFIT 5.3)
- Exact match at m_1 = 0.91 meV

From `sigma_khronon_verification.md`:
- CC sigma (2012) != Khronon (5 structural obstructions)
- CC dilaton (2006) ~ Khronon conformal mode (valid at Level 1)
- Correct picture: TWO fields (sigma from D_F + dilaton/Khronon from conformal mode)

From `mu_breakthrough_2026_03_19.md`:
- mu_0 = H_0/c EXCLUDED by CMB (> 55 sigma)
- BS2025 mu^{-1} = 22.3 Mpc passes CMB
- Three coincidences: mu c^2 = a_0(1+z_dec), mu = 2pi(1+Omega_b)/r_d, M_GC ~ Sigma m_nu

---

## 1. The Logical Structure of the Problem

### 1.1 What We Need to Derive

The gravitational seesaw:
```
mu = (Sigma m_nu)^2 / M_Pl
```

In the spectral action framework, this requires showing that the Khronon mass parameter mu is determined by the neutrino mass eigenvalues of D_F and the Planck mass (which comes from the a_2 coefficient).

### 1.2 Why This Is Hard

The CC dilaton (identified with the Khronon conformal mode) is **massless at tree level** due to its shift symmetry phi -> phi + const. Any mass must be generated dynamically, either through:
- Loop corrections (Coleman-Weinberg mechanism)
- Non-perturbative effects (instanton-like)
- Spontaneous symmetry breaking (ghost condensation)
- Self-consistent fixed point (bootstrap)

The gravitational seesaw formula mu = (Sigma m_nu)^2 / M_Pl has the structure of a **tree-level** relation (like the standard seesaw m_nu = m_D^2/M_R), but there is no obvious tree-level mechanism in the spectral action that produces it.

### 1.3 Dimensional Analysis Constraint

```
[mu] = eV (inverse length in natural units)
[Sigma m_nu] = eV
[M_Pl] = eV
```

Any derivation must produce the specific combination (Sigma m_nu)^2/M_Pl. This requires:
- Two powers of neutrino mass (quadratic coupling)
- One inverse power of M_Pl (gravitational suppression)
- No other dimensionful parameters

---

## 2. Route 1: Direct from D_F Eigenvalues

### 2.1 The Spectrum of D_F

In the standard NCG spectral triple (A_F = C + H + M_3(C)), the finite Dirac operator D_F has eigenvalues:
```
{0, y_e v, y_mu v, y_tau v, y_u v, y_c v, y_t v, y_d v, y_s v, y_b v, y_nu_i v, M_R}
```

where y_i are Yukawa couplings, v = 246 GeV is the Higgs VEV, and M_R ~ 10^{12-15} GeV is the Majorana mass.

### 2.2 Where Is mu?

```
mu = 2.87 x 10^{-31} eV
```

This is **52 orders of magnitude** below the smallest nonzero eigenvalue of D_F (y_e v ~ 0.5 MeV). It does not appear as any eigenvalue, ratio of eigenvalues, or simple algebraic combination of eigenvalues of D_F.

### 2.3 Verdict: FAILS

The Khronon mass mu is not directly readable from the spectral data of the finite geometry. It must be a DERIVED quantity involving both D_F and the manifold geometry.

---

## 3. Route 2: One-Loop Dilaton Mass (Coleman-Weinberg)

### 3.1 The Mechanism

The CC dilaton is massless at tree level. At one loop, fermion loops generate an effective mass through the conformal anomaly. The dilaton couples to fermion masses via the conformal factor: m_f(phi) = m_f * e^{-phi}.

### 3.2 Computation

The effective dilaton 'Yukawa' coupling to a fermion of mass m_f is:
```
g_f = m_f / f_phi
```
where f_phi ~ M_Pl is the dilaton decay constant (for a gravitational dilaton).

The one-loop Coleman-Weinberg contribution to the dilaton mass:
```
delta mu^2 ~ (1/(8pi^2)) * sum_i g_i^2 * m_i^2
           = (1/(8pi^2 M_Pl^2)) * sum_i m_i^4
```

### 3.3 Evaluation

For neutrinos only (m_1 ~ 0, m_2 = 8.68 meV, m_3 = 49.53 meV):
```
sum m_nu^4 = (8.68e-3)^4 + (49.53e-3)^4 = 5.67e-9 + 6.02e-6 = 6.03e-6 eV^4

delta mu = sqrt(6.03e-6 / (8pi^2 * (1.22e28)^2))
         = sqrt(6.03e-6 / (9.63e57))
         = sqrt(6.26e-64)
         = 2.50e-32 eV
```

**This is a factor of ~11 below mu_obs = 2.87e-31 eV.**

### 3.4 The Factor of 4pi Problem

The Coleman-Weinberg formula has scheme-dependent numerical prefactors. The factor 8pi^2 comes from the one-loop integral in dimensional regularization. Different regularization schemes (cutoff, zeta-function, spectral) give different O(1) prefactors.

In the spectral action framework, the one-loop effective potential uses the spectral zeta function regularization, which involves zeta values (zeta(3), zeta(5), ...) rather than 1/16pi^2. The CCSvS coefficients are:
```
c_0 ~ zeta(5)/zeta(4)  (for the a_0 term)
c_1 ~ zeta(3)/zeta(2)  (for the a_2 term)
```

If the relevant coefficient is c_1 ~ zeta(3)/zeta(2) = 1.202/(pi^2/6) = 7.21/pi^2 rather than 1/(8pi^2), the loop factor could be larger by a factor of ~7, partially closing the gap.

### 3.5 The "Almost Works" Formula

A simpler formula:
```
mu ~ m_nu_eff^2 / (4pi M_Pl)
```
where m_nu_eff is an effective neutrino mass scale. With m_nu_eff = 50 meV (atmospheric):
```
mu = (50e-3)^2 / (4pi * 1.22e28) = 1.63e-32 eV
```
Off by a factor of ~18.

With m_nu_eff = Sigma m_nu = 58 meV:
```
mu = (58e-3)^2 / (4pi * 1.22e28) = 2.19e-32 eV
```
Off by a factor of ~13.

### 3.6 Comparison with the Target

| Formula | Value (eV) | Ratio to mu_obs |
|---------|-----------|-----------------|
| m_nu^4 / (8pi^2 M_Pl^2) | 2.5e-32 | 0.087 |
| (Sigma m_nu)^2 / (4pi M_Pl) | 2.2e-32 | 0.077 |
| **(Sigma m_nu)^2 / M_Pl** | **2.78e-31** | **0.97** |
| mu_BS2025 | 2.87e-31 | 1.00 |

The target formula (Sigma m_nu)^2/M_Pl has **no** loop suppression factor 1/(4pi) or 1/(16pi^2). This is the hallmark of a **tree-level** mechanism, not a loop effect.

### 3.7 Verdict: CLOSE BUT OFF BY ~10x

The Coleman-Weinberg mechanism gives the right order of magnitude but with an unwanted factor of 1/(4pi) ~ 0.08. The gravitational seesaw formula has no such factor, suggesting it is NOT a one-loop effect.

---

## 4. Route 3: Non-Perturbative Effects

### 4.1 Instanton Mechanism

If the dilaton mass is generated by non-perturbative effects (analogous to the QCD axion mass from instantons):
```
mu ~ Lambda_NP^{n+4} / M_Pl^{n+3}
```
for some power n. The simplest case n=0 gives:
```
mu ~ Lambda_NP^4 / M_Pl^3
Lambda_NP = (mu M_Pl^3)^{1/4}
```

### 4.2 Evaluation

```
Lambda_NP = (2.87e-31 * (1.22e28)^3)^{1/4}
          = (2.87e-31 * 1.82e84)^{1/4}
          = (5.22e53)^{1/4}
          = 2.69e13 eV = 26.9 TeV
```

### 4.3 Is 27 TeV a Meaningful Scale?

This is above the electroweak scale (246 GeV) but below the GUT scale. It does not correspond to any known spectral action quantity. It is NOT related to Sigma m_nu in any obvious way.

### 4.4 Verdict: NO MECHANISM

The instanton route gives a scale that does not match any spectral data. There is no known non-perturbative effect in the NCG framework that operates at 27 TeV to generate the Khronon mass.

---

## 5. Route 4: Ghost Condensation Self-Consistency (MOST PROMISING)

### 5.1 The Key Idea

Ghost condensation creates an effective IR cutoff at the energy scale M_GC = sqrt(mu M_Pl). This cutoff determines which fermion species contribute to the Khronon mass. If only fermions with m < M_GC contribute, then:
- M_GC acts as an IR filter
- Only neutrinos (m_nu < 60 meV) pass through the filter
- All other fermions (m_e = 0.511 MeV and above) are decoupled
- The Khronon mass is set exclusively by the neutrino sector

### 5.2 Why Ghost Condensation Provides an IR Cutoff

At ghost condensation, the Khronon mode has quartic dispersion:
```
omega^2 = c_4^2 k^4 / M_GC^2 + ...
```

This means the Khronon propagator in the infrared is:
```
G(k) ~ M_GC^2 / (c_4^2 k^4)
```

For virtual fermion loops contributing to the Khronon mass, the loop integral involves the Khronon propagator. The quartic dispersion suppresses contributions from modes with energy E >> M_GC, because the Khronon cannot propagate efficiently at those scales. Conversely, modes with E << M_GC couple strongly to the Khronon condensate.

### 5.3 The Self-Consistency Equation

**Ansatz**: The Khronon mass is determined by the sum of masses of all fermions below the ghost condensation scale:
```
M_GC = sum_{i: m_i < M_GC} m_i    (*)
```

Combined with the definitional relation:
```
M_GC = sqrt(mu M_Pl) = sqrt(M_GC^2 / M_Pl * M_Pl) = M_GC    (**)
```

Equation (*) is the physical content; (**) is tautological. The self-consistency requires:

```
sqrt(mu M_Pl) = sum_{i: m_i < sqrt(mu M_Pl)} m_i
```

### 5.4 Verification with SM Fermion Spectrum

The SM fermion masses in ascending order:
```
nu_1:  ~0 meV          (below M_GC)
nu_2:  8.68 meV        (below M_GC)
nu_3:  49.53 meV       (below M_GC)
e:     0.511 MeV = 511,000 meV  (ABOVE M_GC)
u:     ~2.2 MeV        (ABOVE M_GC)
d:     ~4.7 MeV        (ABOVE M_GC)
...
```

For M_GC ~ 59 meV:
```
sum_{m_i < 59 meV} m_i = 0 + 8.68 + 49.53 = 58.21 meV
```

Check: sqrt(mu M_Pl) = 59.17 meV ~ 58.21 meV.  **Matches to 1.7%.**

The electron (m_e = 511 MeV >> 59 meV) is far above the cutoff. There is a **massive gap** between the heaviest neutrino (49.53 meV) and the electron (0.511 MeV), a factor of ~10,000. This gap ensures the self-consistency is **robust**: M_GC could vary by a factor of 100 and still select only the three neutrinos.

### 5.5 The Fixed-Point Structure

Define:
```
F(M) = sum_{i: m_i < M} m_i
```

This is a step function that jumps at each fermion mass. The self-consistency condition is:
```
M = F(M)
```

The function F(M) has the following structure:
```
M < 0 meV:       F(M) = 0
0 < M < 8.68 meV:  F(M) = 0 (only m_1 ~ 0)
8.68 < M < 49.53 meV: F(M) = 8.68 meV
49.53 < M < 511 MeV:  F(M) = 58.21 meV  <-- M_GC lives HERE
511 MeV < M < 2.2 MeV: F(M) = 58.21 + 0.511 MeV = ...
```

The fixed point M = F(M) requires M = 58.21 meV, which falls in the interval (49.53 meV, 511 MeV). In this interval, F(M) = 58.21 meV = constant. So:
```
M_GC = 58.21 meV (the unique fixed point)
```

**The fixed point is UNIQUE and STABLE.** Any small perturbation of M_GC within the interval (49.53 meV, 511 MeV) returns to the same value. This is a fixed-point bootstrap.

### 5.6 The Gravitational Seesaw Emerges

From the fixed point M_GC = Sigma m_nu:
```
mu = M_GC^2 / M_Pl = (Sigma m_nu)^2 / M_Pl
```

This IS the gravitational seesaw. It emerges as a **self-consistency condition** for the ghost condensation scale, not as a tree-level or loop-level formula.

### 5.7 Connection to the Spectral Action

In the spectral action framework:
1. The Dirac operator D_F determines the fermion mass spectrum {m_i}
2. The a_2 coefficient determines M_Pl^2 = f_2 Lambda^2 / pi
3. The conformal/dilaton mode (from CC 2006) becomes the Khronon
4. Ghost condensation creates the IR cutoff M_GC
5. The self-consistency M_GC = F(M_GC) determines mu from the spectral data

The spectral action provides BOTH ingredients:
- M_Pl from a_2 (gravitational sector)
- {m_nu_i} from D_F (fermion sector)

And the ghost condensation mechanism provides the self-consistency equation that combines them.

### 5.8 Honest Assessment: Circularity

**The argument is self-consistent but CIRCULAR.** The step from "ghost condensation creates IR cutoff at M_GC" to "only fermions with m < M_GC contribute to mu" is physically motivated but not derived from first principles.

The circularity is: we ASSUME M_GC acts as a cutoff, then DERIVE M_GC = Sigma m_nu, then CHECK self-consistency. But we have not PROVEN that ghost condensation creates this specific type of cutoff.

What would be needed for a first-principles derivation:
1. Start from the spectral action with the dilaton/Khronon
2. Compute the effective action at one loop, including all SM fermions
3. Show that the ghost condensation condition K'(Q_0) = 0 combined with the fermion loops gives mu = (Sigma m_nu)^2/M_Pl
4. This requires the fermion loop integral to have the self-consistency structure

### 5.9 Verdict: PROMISING BUT CIRCULAR

The self-consistency argument is beautiful and gives the correct numerical answer to 1.7%. It explains WHY neutrinos are singled out (they are the only fermions below M_GC) and WHY the sum (not individual masses) appears (the Khronon is a scalar that couples universally). But it is not a first-principles derivation.

---

## 6. Route 5: From g_QQ = 2/Q^2 (Paper 9 Result)

### 6.1 What Paper 9 Gives

From `paper9_d4_calculation.md`, the Fisher metric of the a_2 coefficient under conformal variation is:
```
g_QQ^{EH}(Q) = (d-2)(d-3)/Q^2
```

In d=4: g_QQ = 2/Q^2, matching the channel result Sigma = 2 ln Q.

### 6.2 Connection to K(Q)

The Khronon action is:
```
S_K = (M_Pl^2/2) integral sqrt(g) [R + 2K(Q)] d^4x
```

where K(Q) = mu^2(Q-1)^2. The second variation of the K(Q) term at Q=1:
```
K''(1) = 2 mu^2
```

The Fisher metric from the K(Q) term is:
```
g_QQ^{K}(1) = 2 mu^2 / (normalization)
```

### 6.3 Can g_QQ = 2/Q^2 Determine mu?

The Paper 9 result gives g_QQ = 2/Q^2 from the a_2 sector. This constrains the SHAPE of K(Q) near Q=1:
```
K(Q) = (g_QQ(1)/2) * (Q-1)^2 + O((Q-1)^3) = (1/Q^2)|_{Q=1} * (Q-1)^2 + ...
     = 1 * (Q-1)^2 + ...
```

So g_QQ = 2/Q^2 gives K(Q) ~ (Q-1)^2 near Q=1, which matches the ghost condensation form K = mu^2(Q-1)^2 with mu^2 absorbed into the normalization.

**But mu itself is NOT determined by the Fisher metric.** The Fisher metric fixes the SHAPE of K(Q) (quadratic at Q=1) but not the SCALE (the value of mu). This is because the Fisher metric is a dimensionless geometric quantity, while mu is dimensionful.

### 6.4 What Would Determine mu

To determine mu from the spectral action, we would need:
1. The NORMALIZATION of K(Q) relative to the Einstein-Hilbert term
2. This involves the ratio of the K(Q) coefficient to the a_2 coefficient
3. In the spectral action: S = (M_Pl^2/2) int [R + 2K(Q)] sqrt(g), so K(Q) is measured in units of R (curvature)
4. mu^2 has dimensions of (curvature) = (length)^{-2}

The spectral action tells us that K(Q) lives in the a_2 sector (same sector as R), and the Fisher metric matches. But to get mu, we need an additional input: the ratio of the K(Q) amplitude to the curvature R, which is set by the matter content (specifically, the neutrino sector through Route 4).

### 6.5 Verdict: FORM MATCHES, mu UNDETERMINED

The spectral action determines the form K(Q) ~ (Q-1)^2 (ghost condensation) through the Fisher metric constraint g_QQ = 2/Q^2. It does not determine mu.

---

## 7. The Spectral Action Two-Field Picture

### 7.1 Two Fields from the Spectral Triple

Following the analysis in `sigma_khronon_verification.md` Section 10:

**Field 1: CC sigma (from D_F)**
- Origin: Majorana sector of the finite Dirac operator
- Mass: m_sigma ~ M_R ~ 10^{12-15} GeV
- Potential: V(sigma) = lambda_sigma sigma^4 - mu_sigma^2 sigma^2 + lambda_{Hs} |H|^2 sigma^2
- Role: Provides seesaw mechanism (m_nu = m_D^2/M_R) and Higgs mass correction (m_H = 125 GeV)

**Field 2: CC dilaton / Khronon (from conformal mode)**
- Origin: Making the cutoff Lambda dynamical: Lambda -> Lambda e^{-phi}
- Mass: mu ~ 10^{-31} eV (generated dynamically)
- Action: K(Q) = mu^2(Q-1)^2 with Q = function of d phi
- Role: Ghost condensation -> dark matter / MOND

### 7.2 The Sigma-Dilaton Coupling

In the spectral action, the sigma and dilaton couple through the heat kernel expansion. The sigma enters through D_F (finite space), the dilaton through the conformal factor (manifold). The mixed terms in the spectral action create a portal:

```
L_portal ~ (sigma^2 / Lambda_GUT^2) * (d phi)^2
```

On the sigma condensation background (sigma = M_R):
```
L_eff ~ (M_R^2 / Lambda_GUT^2) * (d phi)^2
```

This modifies the dilaton kinetic term but does not directly generate the mass mu.

### 7.3 The Chain: D_F -> sigma -> m_nu -> mu

The chain that produces mu involves two levels:

**Level 1 (UV, at M_R)**: The sigma's VEV gives the seesaw
```
m_nu_i = y_nu_i^2 v^2 / M_R     (standard seesaw from spectral action)
```

**Level 2 (IR, at M_GC)**: Ghost condensation self-consistency (Route 4)
```
M_GC = Sigma m_nu_i = sum of neutrino masses below M_GC
mu = M_GC^2 / M_Pl = (Sigma m_nu)^2 / M_Pl    (gravitational seesaw)
```

The spectral action provides the input at Level 1 (through D_F and the sigma mechanism). Level 2 is a self-consistency condition that requires ghost condensation.

### 7.4 What the Spectral Action Provides (and What It Does Not)

**Provides**:
- The fermion mass spectrum {m_i} (from D_F)
- M_Pl (from a_2 coefficient)
- The form of K(Q) ~ (Q-1)^2 (from g_QQ = 2/Q^2)
- The dilaton as a natural degree of freedom (from conformal mode)
- The seesaw mechanism for neutrino masses (from sigma VEV)

**Does NOT provide (yet)**:
- The ghost condensation mechanism K'(Q_0) = 0 (not in the spectral action)
- The self-consistency equation M_GC = Sigma m_nu (requires ghost condensation)
- The value of mu directly

---

## 8. The Tr(m_nu^2)/M_Pl vs (Sigma m_nu)^2/M_Pl Distinction

### 8.1 The Two Natural Candidates

The spectral action naturally involves **Tr(D_F^2)**, which includes:
```
Tr(M_nu^dag M_nu) = sum_i m_nu_i^2
```

This gives:
```
mu_spectral = Tr(m_nu^2) / M_Pl = (8.68e-3)^2 + (49.53e-3)^2) / M_Pl
            = 2.53e-3 / 1.22e28 = 2.07e-31 eV
```

The gravitational seesaw gives:
```
mu_seesaw = (Sigma m_nu)^2 / M_Pl = (58.21e-3)^2 / M_Pl
          = 3.39e-3 / 1.22e28 = 2.78e-31 eV
```

### 8.2 Numerical Comparison

| Formula | mu (eV) | mu^{-1} (Mpc) | Discrepancy from BS2025 |
|---------|---------|---------------|------------------------|
| Tr(m_nu^2)/M_Pl | 2.07e-31 | 30.9 | +38.6% |
| **(Sigma m_nu)^2/M_Pl** | **2.78e-31** | **23.0** | **+3.3%** |
| mu_BS2025 | 2.87e-31 | 22.3 | reference |

### 8.3 The Algebraic Difference

```
(Sigma m_nu)^2 = sum_i m_i^2 + 2 * sum_{i<j} m_i m_j
               = Tr(m_nu^2) + 2 * sum_{i<j} m_i m_j
```

The cross terms 2 * sum_{i<j} m_i m_j = 2 * (0 * 8.68 + 0 * 49.53 + 8.68 * 49.53) meV^2 = 860 meV^2, contributing ~34% of the total.

### 8.4 Why (Sigma m_nu)^2 and Not Tr(m_nu^2)?

The spectral action's heat kernel naturally gives Tr(D_F^2) ~ sum m_i^2. To get (sum m_i)^2, we need a mechanism that takes the SQUARE OF THE TRACE rather than the TRACE OF THE SQUARE.

This distinction encodes whether the three neutrino species contribute:
- **Independently**: Tr(m^2) = sum m_i^2 (each contributes its own m_i^2)
- **Coherently**: (Tr m)^2 = (sum m_i)^2 (the SUM is squared, with cross terms)

The ghost condensation self-consistency argument (Route 4) naturally gives the coherent sum (sum m_i)^2, because M_GC is defined as the TOTAL mass below the cutoff, then squared.

The spectral action loop correction naturally gives the independent sum Tr(m^2), because each fermion species runs in the loop independently.

**The 34% difference between these two formulae is a discriminant**: the data slightly favor (sum m_nu)^2/M_Pl (3.3% off) over Tr(m_nu^2)/M_Pl (38.6% off), supporting the coherent (self-consistency) mechanism over the incoherent (loop) mechanism.

---

## 9. Attempt at a Clean Derivation

### 9.1 Setup

Consider the spectral action on M^4 x F with:
- Manifold M^4 with metric g_mu_nu, conformal mode Q (dilaton/Khronon)
- Finite space F with Dirac operator D_F having eigenvalues {m_nu_i, m_e, ...}
- Ghost condensation at Q_0 = 1

### 9.2 The Effective Action Below M_R

After integrating out the sigma field at the seesaw scale M_R:
- The neutrino masses m_nu_i = y_nu_i^2 v^2 / M_R are generated
- The dilaton remains massless at tree level
- The one-loop effective action for the dilaton includes contributions from all fermions

### 9.3 The Ghost Condensation Condition

K'(Q_0) = 0 at Q_0 = 1 requires:
```
sum_i (partial/partial Q) [delta S_i(Q)] |_{Q=1} = 0
```

where delta S_i(Q) is the Q-dependent part of the one-loop effective action from fermion species i.

### 9.4 The Loop Integral

For a single fermion of mass m in the Khronon background with Q = 1 + delta:
```
delta S_1-loop = (1/2) Tr ln(D^2(Q) / D^2(1))
```

Under Q-rescaling, D -> D/Q, so:
```
delta S_1-loop = (1/2) Tr ln(D^2/Q^2 / D^2) = -Tr ln Q = -N(Q*Lambda) * ln Q
```

where N(E) is the number of eigenvalues below E.

For the finite-temperature (spectral action) regularization:
```
delta S_1-loop(Q) = sum_n [s(lambda_n^2/(Q*Lambda)^2) - s(lambda_n^2/Lambda^2)]
```

### 9.5 The Fermion Mass Contribution

In the spectral action heat kernel expansion, the Q-dependent contribution from a fermion of mass m is (from the a_2 sector):
```
delta S_m(Q) = c_1 Lambda^2 * m^2 * (Q^2 - 1) / M_Pl^2 + ...
```

This has:
```
d/dQ [delta S_m(Q)]|_{Q=1} = 2 c_1 Lambda^2 m^2 / M_Pl^2
d^2/dQ^2 [delta S_m(Q)]|_{Q=1} = 2 c_1 Lambda^2 m^2 / M_Pl^2
```

### 9.6 The Ghost Condensation Condition from the Loop

If ghost condensation requires the TOTAL first derivative to vanish:
```
sum_i d/dQ [delta S_i(Q)]|_{Q=1} + d/dQ [S_grav(Q)]|_{Q=1} = 0
```

The gravitational contribution (from a_2):
```
d/dQ [S_grav(Q)]|_{Q=1} = 2 * (M_Pl^2/2) * int R sqrt(g) = M_Pl^2 * R_0
```

The fermion contribution:
```
sum_i d/dQ [delta S_i(Q)]|_{Q=1} = 2 c_1 Lambda^2 / M_Pl^2 * sum_i m_i^2
```

Setting these equal:
```
M_Pl^2 R_0 + 2 c_1 Lambda^2 sum m_i^2 / M_Pl^2 = 0
```

This gives R_0 in terms of the fermion masses -- but this is the COSMOLOGICAL CONSTANT problem, not the Khronon mass problem!

### 9.7 Where the Derivation Breaks Down

The issue is that the standard spectral action heat kernel expansion does not naturally produce K(Q) = mu^2(Q-1)^2 with a specific mu. It produces:
- From a_0: a cosmological constant term ~ Q^4
- From a_2: an Einstein-Hilbert term ~ Q^2 R
- From a_4: conformally invariant terms ~ Q^0

The K(Q) function must come from a DIFFERENT mechanism than the heat kernel expansion. The ghost condensation form K(Q) = mu^2(Q-1)^2 has no natural origin in the standard spectral action.

### 9.8 Verdict: DERIVATION INCOMPLETE

A clean first-principles derivation of mu = (Sigma m_nu)^2/M_Pl from the spectral action has not been achieved. The spectral action provides the ingredients (M_Pl, m_nu_i, the dilaton, g_QQ = 2/Q^2) but not the assembly instructions (ghost condensation, self-consistency equation).

---

## 10. Summary: What We Have and What We Need

### 10.1 What We Have (Solid)

1. **Numerical match**: M_GC = sqrt(mu M_Pl) ~ Sigma m_nu to 1.7% (0.8% with updated NuFIT)
2. **Fisher metric from spectral action**: g_QQ = 2/Q^2 in d=4 (Paper 9), matching the channel Sigma = 2 ln Q
3. **Self-consistency argument**: ghost condensation IR cutoff + SM fermion spectrum -> uniquely selects neutrinos -> M_GC = Sigma m_nu (Route 4)
4. **Two-field picture**: spectral action naturally contains sigma (from D_F) + dilaton/Khronon (from conformal mode)
5. **Seesaw chain**: D_F -> M_R -> m_nu (standard seesaw from spectral action, well-established)

### 10.2 What We Need (Open)

1. **First-principles derivation** of K(Q) = mu^2(Q-1)^2 from the spectral action
   - The form K ~ (Q-1)^2 follows from g_QQ = 2/Q^2, but the normalization mu is not determined
   - Ghost condensation K'(Q_0) = 0 has no analog in the heat kernel expansion

2. **First-principles derivation** of the self-consistency equation M_GC = Sigma m_nu
   - The physical argument (IR cutoff from ghost condensation) is compelling but not a theorem
   - Need to show: one-loop effective action with ghost condensation background reproduces this

3. **Why (Sigma m_nu)^2 and not Tr(m_nu^2)?**
   - The data prefer the coherent sum (3.3% off) over the incoherent sum (38.6% off)
   - This distinction encodes whether the mechanism is a self-consistency bootstrap or a loop correction
   - Need a principle that selects (sum m_i)^2 over sum m_i^2

4. **Which Planck mass?**
   - The match requires the FULL Planck mass M_Pl = 1.22e28 eV, not the reduced M_Pl_red = 2.44e27 eV
   - This is consistent with ghost condensation literature (Arkani-Hamed et al. 2004) but needs spectral action verification

### 10.3 The Honest Assessment

The gravitational seesaw mu = (Sigma m_nu)^2/M_Pl is:
- **Numerically verified** to 1.7% (possibly 0.8%)
- **Structurally motivated** by the seesaw analogy and the self-consistency argument
- **Compatible with** the spectral action framework (all ingredients are present)
- **NOT YET DERIVED** from first principles in any framework

The closest to a derivation is Route 4 (self-consistency), which is physically compelling and gives the exact right answer, but is logically circular in its current form.

### 10.4 Path Forward for Paper 9

For Paper 9 (spectral action from Sigma), the recommended strategy:

**Include (well-established)**:
- g_QQ = (d-2)(d-3)/Q^2 from the a_2 sector (PROVEN)
- d=4 selection from g_QQ = 2/Q^2 = g_QQ^{channel} (PROVEN)
- The spectral action naturally contains a dilaton/Khronon mode (ESTABLISHED)

**Discuss (suggestive)**:
- The gravitational seesaw mu = (Sigma m_nu)^2/M_Pl as a CONJECTURE
- The self-consistency argument as MOTIVATION
- The two-field picture (sigma + dilaton) as the CORRECT framework

**Do not claim (not yet proven)**:
- That mu is DERIVED from the spectral action
- That K(Q) = mu^2(Q-1)^2 FOLLOWS from the spectral action
- That ghost condensation is a CONSEQUENCE of the spectral action

---

## 11. Speculative: A Possible Path to the Derivation

### 11.1 The Missing Ingredient: Ghost Condensation from Spectral Geometry

The spectral action Tr f(D^2/Lambda^2) is defined with a FIXED cutoff Lambda. If Lambda becomes dynamical (CC 2006 dilaton), the spectral action gains a conformal degree of freedom.

The ghost condensation condition K'(Q_0) = 0 could arise if the ENTROPY FUNCTIONAL (CCSvS) has a SADDLE POINT in the conformal direction. Specifically:

```
d S_vN / d Q |_{Q_0} = 0    (entropy extremum)
d^2 S_vN / d Q^2 |_{Q_0} > 0    (minimum, not maximum)
```

From Section 5 of paper9_d4_calculation.md:
```
d S_vN / d Q |_{Q=1} = sum_k (4-2k) c_k Lambda^{4-2k} a_{2k}
```

Setting this to zero:
```
4 c_0 Lambda^4 a_0 + 2 c_1 Lambda^2 a_2 + 0 * a_4 + (-2) c_3 Lambda^{-2} a_6 + ... = 0
```

This is a condition on Lambda (or equivalently on the Hubble rate H, since on FRW background a_0 ~ Vol and a_2 ~ H^2 Vol). It determines the epoch at which Q_0 = 1 is an entropy extremum.

### 11.2 The Neutrino Condensation Epoch

When the most massive neutrino becomes non-relativistic:
```
z_nr(m_3 = 49.5 meV) ~ m_3 / (3 T_nu_0) ~ 98
```

At this epoch, the neutrino contribution to the a_2 coefficient changes (from radiation-like to matter-like). This shifts the entropy extremum condition and could trigger the ghost condensation.

If ghost condensation occurs at z ~ 100, then:
```
H(z=100) ~ H_0 * sqrt(Omega_m * (1+100)^3) ~ H_0 * 10^3
```

And the Khronon mass at that epoch:
```
mu(z=100) ~ H(z=100)/c ~ 10^3 * H_0/c ~ 10^3 * 10^{-33} eV ~ 10^{-30} eV
```

This is tantalizingly close to mu_BS = 2.87e-31 eV (off by ~3x).

### 11.3 Why This Might Work

If the ghost condensation is triggered by the neutrino non-relativistic transition:
1. The transition modifies the spectral action at the a_2 level (neutrino mass contribution)
2. The entropy extremum condition dS_vN/dQ = 0 is satisfied at the transition epoch
3. The Khronon mass is set by the neutrino masses at the transition: mu ~ (Sigma m_nu)^2 / M_Pl
4. This is a dynamical mechanism, not an arbitrary parameter choice

### 11.4 What Remains to Be Computed

To make this rigorous:
1. Compute the full a_2 coefficient including massive neutrinos (non-relativistic correction)
2. Show that dS_vN/dQ = 0 is satisfied at z ~ z_nr
3. Show that d^2S_vN/dQ^2 > 0 (ghost condensation, not anti-condensation)
4. Extract mu from the second derivative and show it equals (Sigma m_nu)^2/M_Pl

This is a well-defined computation that could be carried out. It would constitute the first-principles derivation.

---

## Key References

1. Chamseddine, A.H. & Connes, A. "Scale Invariance in the Spectral Action." J. Math. Phys. 47, 063504 (2006). [The dilaton paper]
2. Chamseddine, A.H. & Connes, A. "Resilience of the Spectral Standard Model." arXiv:1208.1030 (2012). JHEP 09, 104 (2012). [The sigma paper]
3. Chamseddine, A.H., Connes, A. & van Suijlekom, W.D. "Entropy and the Spectral Action." arXiv:1809.02944 (2018). Commun. Math. Phys. 373, 457-471 (2020). [CCSvS: entropy = spectral action]
4. Arkani-Hamed, N. et al. "Ghost Condensation and a Consistent Infrared Modification of Gravity." hep-th/0312099 (2004). JHEP 05, 074 (2004). [Ghost condensation]
5. Blanchet, L. & Skordis, C. "Khronon-Tensor theory." arXiv:2507.00912 (2025). [mu^{-1} = 22.3 Mpc]
6. Vassilevich, D.V. "Heat kernel expansion: user's manual." hep-th/0306138 (2003). [Heat kernel coefficients]

---

*Analysis completed 2026-03-19.*
*Status: The gravitational seesaw mu = (Sigma m_nu)^2/M_Pl is numerically verified (1.7%) and structurally motivated (self-consistency argument, spectral action ingredients). A first-principles derivation remains OPEN, with a concrete path identified (entropy extremum condition + neutrino non-relativistic transition). The closest approach is the ghost condensation self-consistency (Route 4), which gives the correct answer but is logically circular.*
*Classification: SUGGESTIVE -- upgraded from COINCIDENCE by the self-consistency argument and the spectral action connection.*
