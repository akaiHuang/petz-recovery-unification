# Paper 3: Galaxy Rotation Curves from the tau Framework
# A Rigorous Quantitative Analysis

**Author**: Sheng-Kai Huang
**Date**: 2026-03-11
**Status**: Research analysis (pre-draft)

---

## 0. Executive Summary

This document provides a rigorous, quantitative analysis of whether the tau framework -- where gravitational entropy production Sigma_grav = -ln(g_00) arises from a quantum channel -- can explain galaxy rotation curves WITHOUT dark matter particles. The core mechanism is that quantum corrections to the gravitational channel at galactic scales produce a logarithmic modification to the Newtonian potential.

**Bottom line**: The framework is viable at galactic scales (rotation curves, RAR, Tully-Fisher), with a natural connection between quantum information theory and running G. The CMB acoustic peak problem remains unsolved and is the single most serious obstacle to a complete dark-matter-free cosmology.

### Scorecard

| Test | Status | Notes |
|------|--------|-------|
| Flat rotation curves | PASS | Logarithmic correction from eta=1 running |
| Radial Acceleration Relation | PASS | Reproduces g_obs = F(g_bar) with correct a_0 |
| Baryonic Tully-Fisher (v^4 ~ M) | PASS | Follows from v^2 = 2G_N M k_*/pi |
| Solar system constraints | PASS | Running negligible at r << r_c |
| Dwarf spheroidal galaxies | PASS | Verlinde component favored over MOND at 5.2 sigma |
| Galaxy clusters | UNCERTAIN | No systematic calculation yet |
| Bullet Cluster | UNCERTAIN | Running G preserves GR lensing; offset unclear |
| CMB acoustic peaks | FAIL | Running G too weak at z ~ 1100 |
| Structure formation | FAIL | Cannot drive growth without non-baryonic matter |
| Gravitational waves (v_GW = c) | PASS | Framework built on GR |

---

## 1. The tau Framework at Galactic Scales

### 1.1 Review: Sigma_grav in the Strong Field (Paper 2)

From Paper 2, the gravitational entropy production for a Schwarzschild geometry is:

```
Sigma_grav = -ln(g_00) = -ln(1 - r_s/r)
```

In the weak-field limit (r >> r_s):

```
Sigma_grav ≈ r_s/r = 2GM/(c^2 r)
```

This was derived from three independent routes:
1. Modular flow (Dorau-Much extension): fractional QRE loss = r_s/r [EXACT]
2. Gravitational Landauer (Herrera 2020): Tolman-modified erasure cost
3. Quantum channel (pure loss bosonic): eta = exp(-r_s/r), Sigma = -ln(eta)

### 1.2 The Key Extension: Scale-Dependent Gravitational Channel

The gravitational channel N_grav processes quantum information. At different momentum scales k (corresponding to distance scales r ~ 1/k), the channel has different properties because the gravitational coupling G itself runs under the renormalization group (RG):

```
N_grav(k): rho --> Tr_env[ U(k) rho U(k)^dag ]
```

where U(k) is the unitary evolution with the Hamiltonian H(G(k)), and the environment trace represents gravitational decoherence.

The entropy production at scale k is:

```
Sigma(k) = D(rho || sigma) - D(N_grav(k)(rho) || N_grav(k)(sigma))    (1)
```

By the Data Processing Inequality (DPI): Sigma(k) >= 0 for all k.

**The crucial insight**: As k decreases (probing larger distances), G(k) increases (IR antiscreening), making the channel "noisier." This produces MORE information loss, hence larger Sigma, hence stronger effective gravity -- which manifests as what we interpret as "dark matter."

### 1.3 Modified Sigma at Galactic Scales

With the running gravitational coupling G(k), the effective potential becomes:

```
Phi(r) = -G(r) M / r
```

where G(r) is the position-space version of G(k) obtained by Fourier transform. For the marginal case eta=1 (see Section 2):

```
G(r) = G_N [1 + (2k_* r / pi) + ...]    for r >> 1/k_*
```

Therefore:

```
Sigma_grav(r) = 2|Phi(r)|/c^2 = (r_s/r) + (2k_* r_s / pi) ln(r/r_0)    (2)
```

where:
- r_s = 2G_N M/c^2 is the standard Schwarzschild radius
- k_* is the crossover scale (IR running onset)
- r_0 is a reference scale

**This is NOT an ad hoc addition.** It is the SAME formula Sigma = 2|Phi|/c^2, but with the quantum-corrected potential that follows from the RG-improved gravitational channel.

---

## 2. Kumar's Running G: Detailed Quantitative Analysis

### Reference: Kumar (2025), arXiv:2509.05246, Phys. Lett. B 871 (2025) 140008

### 2.1 The Running Coupling

The anomalous dimension eta governs the IR running:

```
d ln G / d ln k = -eta                                                  (3)
```

Solution:

```
G(k) = G_N (k_*/k)^eta    for k << k_*                                  (4)
```

### 2.2 Why eta = 1 is Unique (The Marginal Case)

The modified Newtonian potential in momentum space:

```
Phi(k) = -4 pi G(k) M / k^2 = -4 pi G_N M k_*^eta / k^(2+eta)         (5)
```

Fourier transform to position space in d=3:

```
Phi(r) = -M * integral [d^3k/(2pi)^3] * e^{ik.r} * 4 pi G(k) / k^2     (6)
```

The integral over k^{-(2+eta)} in 3D gives:

| eta | Phi(r) behavior | Physical status |
|-----|----------------|-----------------|
| eta < 1 | r^{eta-1} (power-law, decaying) | IR-irrelevant: correction dies off |
| **eta = 1** | **ln(r)** (logarithmic) | **Marginal: universal, scale-invariant** |
| eta > 1 | r^{eta-1} (power-law, growing) | Violates scale invariance |

**eta = 1 is the UNIQUE value producing logarithmic behavior** -- the only scale-invariant deformation consistent with rotational symmetry and locality.

### 2.3 Explicit Formulas

For eta = 1, G(k) = G_N [1 + k_*/k]:

```
Potential:   Phi(r) = -G_N M/r + (2G_N M k_*/pi) ln(r/r_0)              (7)

Force:       F(r) = -G_N M/r^2 - (2G_N M k_*)/(pi r)                    (8)

Velocity:    v^2(r) = r|F(r)| = G_N M/r + 2G_N M k_*/pi                 (9)
```

At r >> r_c = 1/k_*:

```
v^2(r) ≈ V_0^2 = 2G_N M_bar k_* / pi = const                           (10)
```

**This is a flat rotation curve.**

### 2.4 The Crossover Scale

From Eq. (10), the crossover momentum is:

```
k_* = pi V_0^2 / (2G_N M_bar)                                           (11)
```

Fitted values for three SPARC galaxies:

| Galaxy | M_bar [M_sun] | V_0 [km/s] | k_* [kpc^{-1}] | r_c [kpc] |
|--------|--------------|------------|----------------|-----------|
| S:700624 | 1.01 x 10^12 | 271.4 +/- 2.6 | (2.66 +/- 0.05) x 10^{-2} | 37.5 +/- 0.7 |
| S:700916 | 5.85 x 10^11 | 210.4 +/- 2.4 | (2.76 +/- 0.06) x 10^{-2} | 36.2 +/- 0.8 |
| S:705253 | 3.55 x 10^11 | 158.9 +/- 1.9 | (2.60 +/- 0.06) x 10^{-2} | 38.5 +/- 0.9 |

**Remarkable**: k_* is consistent across galaxies spanning a factor ~3 in mass. This universality is a non-trivial prediction -- NFW halos have NO such universal scale.

### 2.5 Free Parameters

**At the galactic level**: 2 parameters
- k_* (or equivalently V_0): crossover scale
- Upsilon (mass-to-light ratio): determines M_bar

The coefficient 2/pi in the logarithmic term is **universal and regulator-independent**.

**Comparison**:
- NFW dark matter halo: 2 parameters (M_200, concentration c)
- MOND: 1 parameter (a_0, universal) + Upsilon
- Kumar: 1 parameter (k_*, approximately universal) + Upsilon
- Our target: 0 free parameters beyond Upsilon, if k_* can be derived from Sigma

### 2.6 Critical Limitations

1. **CMB**: FAILS to reproduce acoustic peak structure. The IR running effect contributes only Omega_L0 < 0.05 at z ~ 1100 -- negligible compared to Omega_CDM ~ 0.26.

2. **Cosmological evolution**: Modified Friedmann equation:
```
H^2(z)/H_0^2 = Omega_r(1+z)^4 + Omega_m(1+z)^3 + Omega_Lambda + Omega_L0(1+z)^2 ln(1+z)
```
The logarithmic term contributes a new component scaling as (1+z)^2 ln(1+z), but it is severely constrained.

3. **Structure formation**: Cannot drive large-scale structure growth without non-baryonic matter.

4. **eta = 1 is argued, not derived**: Kumar uses EFT marginality/dimensional arguments. This is a self-consistency argument, not a calculation from a UV-complete theory.

---

## 3. Gubitosi et al.: 100-Galaxy SPARC Validation

### Reference: Gubitosi, Piattella & Casarini (2024), arXiv:2403.00531, Phys. Rev. D 110 (2024) 124014

### 3.1 The RGGR Model

Running coupling:

```
G(mu) = G_0 / [1 + nu * ln(mu^2/mu_0^2)]                                (12)
```

Effective circular velocity:

```
v^2_RGGR(r) = v^2_Newton(r) * [1 - c^2 nu_bar / Phi_N(r)]               (13)
```

where nu_bar = nu * alpha is the phenomenological running parameter.

### 3.2 Galaxy Sample and Results

- **100 galaxies** from SPARC (after quality cuts from 175)
- 4 morphological types: Early (9), Spiral (39), Late dwarf (34), Starburst (18)
- MCMC fit with 3 parameters per galaxy: gamma_d, gamma_b, nu_bar

**BIC comparison (RGGR vs NFW)**:

| Preference | Count |
|-----------|-------|
| RGGR preferred (Delta_BIC > 2) | 47 galaxies |
| NFW preferred (Delta_BIC < -2) | 40 galaxies |
| Inconclusive | 13 galaxies |

**RGGR is competitive with NFW overall**, with slight preference.

### 3.3 Key Finding: nu_bar Correlates with Mass

```
nu_bar ~ 10^{-6}   (Early type, M ~ 10^{12} M_sun)
nu_bar ~ 10^{-7}   (Spiral, M ~ 10^{11} M_sun)
nu_bar ~ 10^{-8}   (Late/Starburst, M ~ 10^{10} M_sun)
nu_bar ~ 10^{-17}  (Solar system)
```

The near-linear dependence nu_bar ~ M_bar^alpha (alpha ~ 1) is an empirical scaling that any theory must explain.

### 3.4 Satisfies Empirical Relations

Both the **Radial Acceleration Relation** and the **Baryonic Tully-Fisher Relation** are satisfied by the RGGR fits, confirming that running G is phenomenologically viable.

### 3.5 Weaknesses

1. nu_bar is NOT universal -- varies per galaxy (correlates with mass)
2. 3 free parameters per galaxy (same as NFW)
3. Not a dramatic improvement over dark matter halos

### 3.6 Connection to tau Framework

If nu_bar can be derived from the QRE channel properties (specifically, the scale-dependent Sigma), then the mass dependence becomes a **prediction** rather than a free parameter. The key connection:

```
nu_bar = (1/2) * dSigma/d(ln mu^2) evaluated at galactic scale
```

This would reduce free parameters from 3 to 1 (mass-to-light ratio only).

---

## 4. Verlinde's Emergent Gravity

### Reference: Verlinde (2016), arXiv:1611.02269, SciPost Phys. 2 (2017) 016

### 4.1 Mechanism

In de Sitter space, entanglement entropy has both area-law and volume-law contributions:

```
S_entanglement = S_area + S_volume
```

- S_area = A/(4G hbar): standard Bekenstein-Hawking
- S_volume: thermal volume-law from de Sitter temperature T_dS = H_0/(2 pi k_B)

Matter displaces the volume-law entanglement, creating an "entropy displacement" that manifests as additional gravitational force.

### 4.2 Dark Matter Formula

For spherically symmetric baryonic mass M_B(r):

```
M_D^2(r) = (c H_0 / (6G)) * r^2 * d/dr [r * M_B(r)]                   (14)
```

For a point mass:

```
M_D(r) = sqrt(c H_0 M_B r / (6G))                                       (15)
```

At large r where g_Newton << a_0:

```
g_obs ~ sqrt(a_0 * g_Newton)                                              (16)
```

with a_0 ~ c H_0 / 6 ~ 1.1 x 10^{-10} m/s^2 (close to Milgrom's 1.2 x 10^{-10}).

### 4.3 Mapping to tau Framework

The Verlinde contribution to Sigma:

```
Sigma_Verlinde(r) = (r_s/r) + (r^2 / L_dS^2) * (M_B / M_dS)            (17)
```

where L_dS = c/H_0 and M_dS = c^3/(2GH_0).

In the deep MOND regime (g << a_0):

```
v^2(r) = (c^2/2) * r * dSigma_Verlinde/dr
       = G_N M/r + sqrt(a_0 G_N M / 6)                                   (18)
```

### 4.4 Observational Status (2024-2026)

**Supporting evidence**:
- Yoon, Park & Hwang (2023): 175 SPARC galaxies, good fit
- **Ghari & Haghi (2026, arXiv:2601.01715): 23 dwarf spheroidals, Verlinde favored over MOND at 5.2 sigma**
- Brouwer et al. (2021): KiDS-1000 weak lensing extends RAR by 2 decades

**Challenging evidence**:
- Lelli, McGaugh & Schombert (2017): Requires lower M/L ratios
- Hossenfelder (2017): Covariant version reduces to GR
- Dai & Stojkovic (2017): Internal inconsistencies
- Cluster cores: > 5 sigma rejection

### 4.5 Free Parameters

**Zero free parameters** (given M_B and H_0) for the point-mass formula. Extended mass distributions require specifying M_B(r), which comes from photometry.

### 4.6 Critical Assessment

**Strengths**: Zero free parameters; natural a_0 ~ cH_0; derives Tully-Fisher; works well for dwarf spheroidals.

**Weaknesses**: No covariant formulation; no CMB prediction; fails at cluster cores; theoretical self-consistency questioned.

---

## 5. Other Approaches: Comparative Analysis

### 5.1 Quantum Modified Inertia (Gillot 2025, arXiv:2507.11524)

Modified dynamics with minimal and maximal acceleration bounds:

```
F = m_i * a / sqrt(1 - (a_min/a)^2) * sqrt(1 - (a/a_max)^2)            (19)
```

with a_min linked to the Hubble parameter (a_min ~ cH_0).

**Status**: Good fit to Milky Way rotation curve. Better agreement with RAR than MOND in the 10^{-10} m/s^2 regime. Still phenomenological -- no first-principles derivation.

### 5.2 Many-Body Gravity (Penner 2024, arXiv:2403.13019)

Merges thermal gradient metric variation with Einstein gravity in 5D space-time-temperature. Claims to explain rotation curves, RAR, wide binary stars, AND the Bullet Cluster.

**Status**: Interesting but speculative. Novel 5D framework not widely accepted. Claims regarding Bullet Cluster need independent verification.

### 5.3 Oppenheim's Stochastic Gravity (arXiv:2402.19459)

Classical + stochastic metric in low-acceleration regime.

**Status**: RETRACTED claim. Critical analysis (arXiv:2404.13037) showed the theory "does not in fact lead to modified Newtonian dynamics as claimed."

### 5.4 Bianconi's Gravity from Entropy (2025, PRD 111, 066001)

Entropic action = QRE between spacetime metric and matter-induced metric.

**Status**: The G-field may be a dark matter candidate. Predicts small positive Lambda. **Most direct connection to our Sigma = D(rho_spacetime || rho_matter).** However, no explicit galactic-scale calculation exists yet.

### 5.5 Tsallis-Renyi Entropy to MOND (arXiv:2505.03061)

Modified Renyi entropy --> MOND-like force law + Landauer principle for black holes.

**Status**: Explicitly uses modified entropy (not energy) to derive MOND. Strong conceptual connection to our framework. Not fitted to actual data.

### 5.6 Grumiller's Rindler Model (PRL 2010)

```
V(r) = -GM/r + a_Rindler * r                                             (20)
```

with a_Rindler ~ 0.30 x 10^{-10} m/s^2 (close to a_0).

**Status**: Simple effective model. The linear Rindler term explains ~10% of rotation curve anomaly. Too simple for full galaxy fitting.

### 5.7 Comparative Table

| Approach | v(r) or Phi(r) formula | Free params (per galaxy) | SPARC fit quality | Tully-Fisher | RAR | Weakest point |
|----------|----------------------|-------------------------|-------------------|-------------|-----|---------------|
| **NFW (LCDM)** | Phi_NFW(r; M_200, c) | 2 (+ Upsilon) | Excellent | Not predicted | Emergent | No prediction of RAR tightness |
| **MOND** | g = g_bar/[1-exp(-sqrt(g_bar/a_0))] | 0 (+ Upsilon) | Good | **Predicted** | **Predicted** | CMB, clusters, no covariant form |
| **Kumar (running G)** | Phi = -GM/r + (2GMk_*/pi)ln(r) | 1 (k_*, ~universal) + Upsilon | 3 galaxies only | Follows from Eq.(10) | Follows | CMB; eta=1 not derived |
| **RGGR (Gubitosi+)** | v^2 = v_N^2[1-c^2 nu_bar/Phi_N] | 1 (nu_bar) + Upsilon | Good (100 gal) | Consistent | Consistent | nu_bar not universal |
| **Verlinde** | g = sqrt(a_0 g_bar/6) at low g | 0 (+ Upsilon) | Partial | **Predicted** | **Predicted** | Cluster cores; no covariant form |
| **tau framework** | Sigma = r_s/r + alpha*ln(r/r_c) | 0-1 (+ Upsilon) target | Pending | Follows | Follows | eta=1 derivation; CMB |
| **QMI (Gillot)** | Modified inertia with a_min, a_max | 1-2 | Good (MW) | Consistent | Better than MOND | Phenomenological |

---

## 6. Empirical Constraints: SPARC and RAR

### 6.1 The SPARC Database

**Lelli, McGaugh & Schombert (2016)**, arXiv:1606.09251

- 175 nearby disk galaxies
- 3.6 micron Spitzer photometry (stellar mass)
- HI/H-alpha rotation curves
- Standard fitting: v^2_tot = Upsilon_d v^2_disk + Upsilon_b v^2_bulge + v^2_gas + v^2_extra

### 6.2 The Radial Acceleration Relation (RAR)

**McGaugh, Lelli & Schombert (2016)**, PRL 117, 201101

The tightest empirical constraint:

```
g_obs = g_bar / (1 - exp(-sqrt(g_bar/g_dagger)))                         (21)
```

where g_dagger = a_0 = 1.20 +/- 0.02 x 10^{-10} m/s^2.

**Observed scatter**: 0.13 dex over 2693 data points in 153 galaxies.

**ANY viable theory must reproduce this relation** with comparable or better scatter.

### 6.3 Baryonic Tully-Fisher Relation (BTFR)

```
M_bar = A * V_flat^4                                                      (22)
```

with A ~ 50 M_sun/(km/s)^4 and intrinsic scatter < 0.1 dex.

**Derivation from Kumar**: From Eq. (10):

```
V_0^2 = 2G_N M_bar k_* / pi

=> M_bar = (pi / 2G_N k_*) * V_0^2

=> If k_* ~ V_0^2 / M_bar (from Eq. 11):
   M_bar = pi V_0^4 / (4 G_N^2 M_bar k_*^2)
```

For universal k_*, the relation M_bar ~ V_0^4 follows directly if k_* is a constant.

**Derivation from Verlinde**: From Eq. (16):

```
v^4 = a_0 G_N M_bar = (cH_0/6) G_N M_bar

=> M_bar = 6 v^4 / (cH_0 G_N)
```

This is the BTFR with A = 6/(cH_0 G_N), in reasonable agreement with observations.

### 6.4 BIG-SPARC (Haubner, Lelli et al. 2024, arXiv:2411.13329)

~4000 galaxies (20x larger), providing much stronger statistical constraints. This dataset will be critical for discriminating between running G, MOND, Verlinde, and LCDM.

### 6.5 Recent RAR Tension: Early-Type vs Late-Type

Brouwer et al. (2021) found a >6 sigma difference between RAR slopes for early-type vs late-type galaxies at fixed stellar mass. **This is a challenge for ANY modified gravity theory** that predicts a universal RAR independent of galaxy properties. Running G may naturally accommodate this through the nu_bar mass dependence found by Gubitosi et al.

---

## 7. The tau Framework Synthesis

### 7.1 Core Thesis

The quantum relative entropy drop Sigma = D(rho||sigma) - D(N(rho)||N(sigma)), computed using the full quantum-corrected gravitational channel N_grav(k), naturally produces:

1. **r << r_c**: Standard Newtonian gravity (Sigma ~ r_s/r, Paper 2)
2. **r >> r_c**: Logarithmic correction (flat rotation curves, Paper 3)
3. **Transition**: Controlled by a crossover scale k_* ~ 2.7 x 10^{-2} kpc^{-1}

### 7.2 The Logical Chain

```
[Established] Casini-Huerta (2017): RG flow IS a quantum channel;
              QRE monotonicity under RG = DPI
                            |
                            v
[Established] Dorau-Much (2025/2026, PRL): QRE --> semiclassical
              Einstein equations at FIXED scale
                            |
                            v
[NEW - Paper 3] Combine: QRE at EACH SCALE k --> scale-dependent
              Einstein equations --> running G(k)
                            |
                            v
[Established] Kumar (2025, PLB): Running G with eta=1 -->
              logarithmic potential --> flat rotation curves
                            |
                            v
[Established] Gubitosi+ (2024, PRD): RGGR fits 100 SPARC galaxies
              competitively with NFW
```

### 7.3 Mathematical Framework

**Step 1: Scale-dependent Sigma**

Define the gravitational channel at scale k:

```
N_grav(k): rho --> Tr_{E(k)} [ U(G(k)) rho U(G(k))^dag ]               (23)
```

The entropy production at scale k:

```
Sigma(k) = D(rho || sigma) - D(N_grav(k)(rho) || N_grav(k)(sigma))       (24)
```

**Step 2: DPI constrains the running**

By DPI (Casini-Huerta):

```
Sigma(k_IR) >= Sigma(k_UV) >= 0                                          (25)
```

The gravitational channel at IR scales is NOISIER (loses more information) than at UV scales. Therefore Sigma INCREASES toward the IR.

**Step 3: Connection to g_00**

Using Sigma = -ln(g_00) (Paper 2 result):

```
g_00(k) = exp(-Sigma(k))                                                  (26)
```

At the IR scale (galactic):

```
g_00(r) = 1 - r_s/r - (4G_N M k_*)/(pi c^2) * ln(r/r_0)                (27)
```

**Step 4: Rotation curve**

```
v^2(r) = (c^2/2) * r * d(-g_00)/dr
       = G_N M/r + (2G_N M k_*)/pi                                       (28)
```

This is Eq. (9) -- flat rotation curves follow directly.

### 7.4 Can eta = 1 Be Derived from Information Theory?

This is the hardest and most important open question.

**Argument A (Marginality = DPI saturation boundary)**: eta = 1 is the boundary between information-recoverable channels (eta < 1: Petz map can partially reverse the RG flow) and information-irrecoverable channels (eta > 1: channel capacity goes to zero). The physical requirement that the gravitational channel retains SOME information content at all scales forces eta <= 1, and the requirement that the IR effect is non-vanishing forces eta >= 1. Together: eta = 1.

**Argument B (Conformal invariance)**: If the gravitational channel flows to an IR fixed point that is conformally invariant (like a free field), the anomalous dimension at the fixed point is constrained. In d = 3+1, conformal invariance of the graviton propagator gives eta_IR = 1.

**Argument C (DPI bound conjecture)**:

```
CONJECTURE: eta <= 1 - exp(-Sigma_grav/2) = tau_grav                     (29)
```

Properties:
- Sigma -> 0 (flat space): eta -> 0 (no running in flat space) [CORRECT]
- Sigma -> infinity (strong field): eta -> 1 [CONSISTENT with Kumar]
- The bound saturates at eta = 1 when Sigma is large enough [NATURAL]

**Status**: All three arguments are plausible but NONE constitutes a proof. This remains a genuine gap.

### 7.5 The Universal Acceleration Scale a_0

Why does a_0 ~ c H_0 appear? In the tau framework:

The de Sitter horizon sets an IR cutoff for the gravitational channel. The modular Hamiltonian on the de Sitter horizon gives the thermal temperature T_dS = H_0/(2 pi k_B). The information loss rate cannot exceed the rate set by this temperature:

```
dSigma/dr|_max ~ k_B T_dS / (mc^2 r) ~ H_0 / (2 pi c)                  (30)
```

The crossover from Newtonian to modified behavior occurs when:

```
g_Newton = GM/r^2 ~ dSigma/dr * c^2/2 ~ c H_0 / (4 pi)                 (31)
```

This gives:

```
a_0 ~ c H_0 / (4 pi) ~ 1.7 x 10^{-10} m/s^2                            (32)
```

Close to but not exactly Milgrom's value (1.2 x 10^{-10}). The factor of ~1.4 discrepancy may come from geometric factors in the channel calculation.

**Verlinde's argument gives**: a_0 = cH_0/6 ~ 1.1 x 10^{-10} (closer).

---

## 8. The CMB Problem: An Honest Assessment

### 8.1 Why Running G Fails at the CMB

The CMB acoustic peak structure requires:
1. **Baryonic matter**: couples to photons, provides pressure -> odd peaks
2. **Non-baryonic matter**: provides gravity without pressure -> drives even peaks
3. **The ratio of odd/even peak heights** depends on Omega_b/Omega_total

Running G modifies the STRENGTH of gravity but does NOT provide a second species with different pressure properties. At z ~ 1100:

```
Omega_running-G(z~1100) = Omega_L0 * (1+z)^2 * ln(1+z) / H^2(z)
```

With Kumar's constraint Omega_L0 < 0.05:

```
Omega_running-G(z~1100) < 0.05 * (1100)^2 * ln(1100) / H^2(1100) ~ 0.003
```

This is **negligible** compared to the required Omega_CDM ~ 0.26.

### 8.2 The Fundamental Issue

The CMB third acoustic peak requires non-baryonic matter that:
- Does not interact with photons (no radiation pressure)
- Does interact gravitationally (drives potential wells)
- Has a sound speed c_s << c (unlike photon-baryon fluid)

**No modification of G(r) can replicate this behavior.** This is not specific to the tau framework -- it applies to ALL modified gravity theories (MOND, Verlinde, etc.).

### 8.3 Possible Resolutions (Speculative)

**Resolution A: Accept partial victory**
- Running G explains galactic phenomenology
- Some form of non-baryonic matter (neutrinos? sterile neutrinos? gravitinos?) handles the CMB
- The tau framework constrains what this matter IS, via channel properties

**Resolution B: Volume-law entanglement at early times**
- Verlinde's volume-law entropy could act as effective non-baryonic matter at z ~ 1100
- The de Sitter temperature at recombination: T_dS(z~1100) ~ H(z~1100)/(2 pi k_B)
- This has NOT been calculated; it may or may not provide the right peak ratios

**Resolution C: Non-Markovian channel at early times**
- The gravitational channel may have memory effects (non-Markovian) at short timescales
- At z ~ 1100, the Hubble time is ~380,000 years -- memory effects could be significant
- No existing calculation supports this

**Resolution D: The channel IS the dark matter**
- Bianconi's G-field (the Lagrange multiplier enforcing entropic consistency) acts as effective dark matter at all epochs
- At galactic scales: appears as running G
- At CMB epoch: appears as additional pressureless fluid
- Most promising conceptually, but no CMB power spectrum calculation exists

### 8.4 Assessment

**The CMB problem CANNOT be solved within Paper 3.** The most honest approach:

1. Paper 3 focuses on galactic-scale success
2. Acknowledge CMB as an open problem
3. Note that the SAME limitation applies to MOND, Verlinde, and all modified gravity
4. Point to Bianconi's framework as the most promising direction for Paper 4
5. Note: even LCDM required decades to get CMB right

---

## 9. DPI-Based Constraints on the Gravitational Channel

### 9.1 The Fawzi-Renner Recovery Bound

From Paper 1, the Fawzi-Renner bound gives:

```
Sigma >= -2 ln F(rho, R_P circ N(rho))                                   (33)
```

where R_P is the Petz recovery map and F is the fidelity.

Applied to the gravitational channel at scale k:

```
Sigma(k) >= -2 ln F(rho, R_P(k) circ N_grav(k)(rho))                    (34)
```

### 9.2 Implications for Running G

If the Petz recovery fidelity F_P(k) varies with scale:

```
F_P(k_IR) < F_P(k_UV)    (more information lost at IR)                   (35)
```

This is equivalent to:

```
Sigma(k_IR) > Sigma(k_UV)                                                (36)
```

The RATE of Sigma increase with scale is related to the anomalous dimension:

```
d Sigma / d ln(1/k) = eta * (r_s / r)    (at leading order)              (37)
```

### 9.3 DPI Chain

Consider three scales: k_UV > k_mid > k_IR. The DPI gives:

```
D(rho || sigma) >= D(N_UV(rho) || N_UV(sigma))
                >= D(N_mid(rho) || N_mid(sigma))
                >= D(N_IR(rho) || N_IR(sigma))                            (38)
```

This monotonicity IS the c-theorem of Casini-Huerta, applied to the gravitational channel. Each step in the RG flow loses more information, producing more gravitational entropy Sigma.

### 9.4 What the DPI Cannot Do

The DPI constrains the DIRECTION of information flow (UV -> IR) but does not uniquely determine the RATE. Therefore:

- DPI proves Sigma(k_IR) >= Sigma(k_UV) [YES]
- DPI constrains eta >= 0 [YES]
- DPI does NOT uniquely fix eta = 1 [HONEST ADMISSION]

Additional physical input (conformal invariance, marginality, or an explicit channel calculation) is needed to pin down eta = 1.

---

## 10. Non-Perturbative Effects and the IR

### 10.1 Gravitational Instantons

Non-perturbative effects in quantum gravity (gravitational instantons, topology change) could contribute to Sigma at large distances. The instanton contribution scales as:

```
Sigma_inst ~ exp(-S_inst) ~ exp(-1/(G_N * k^2))                          (39)
```

At galactic scales (k ~ 10^{-2} kpc^{-1} ~ 10^{-23} GeV):

```
1/(G_N k^2) ~ M_Pl^2/k^2 ~ 10^{84}
```

So Sigma_inst ~ exp(-10^{84}) ~ 0. **Non-perturbative effects are utterly negligible at galactic scales.** The running G is a perturbative phenomenon.

### 10.2 Vacuum Entanglement (Verlinde's Contribution)

The volume-law entanglement of the de Sitter vacuum IS a non-perturbative effect, but of a different kind. It comes from the macroscopic entanglement structure of the cosmological horizon, not from individual instantons.

```
S_volume ~ (V / L_dS^3) * (L_dS / l_P)^2 ~ V / (G_N L_dS)              (40)
```

This contributes to Sigma at the level:

```
Sigma_dS(r) ~ r / L_dS ~ r * H_0 / c                                     (41)
```

At r ~ 10 kpc: Sigma_dS ~ 10 kpc * H_0/c ~ 10^{-6}. This is comparable to Sigma_Newton at galactic outskirts.

### 10.3 The Three-Term Approximation

Combining all contributions:

```
Sigma_total(r) ≈ (r_s/r) + alpha * ln(r/r_c) + beta * (r/L_dS)          (42)
```

- Term 1: Classical Newtonian gravity (dominant near mass, r << r_c)
- Term 2: IR running of G (dominant at galactic scales, r ~ r_c)
- Term 3: De Sitter volume-law (dominant at cosmological scales, r ~ L_dS)

**IMPORTANT**: This three-term form is an APPROXIMATION. The fundamental equation is:

```
Sigma = D(rho_spacetime || rho_matter)                                    (43)
```

The three terms emerge as limits of this one QRE in different regimes.

---

## 11. What Makes This Different from "Just Running G"

Running G has been proposed for rotation curves since Reuter & Weyer (2004). What is NEW in the tau framework:

### 11.1 The Running is Not Put in by Hand

In standard asymptotic safety, G(k) is obtained from truncated FRG equations. The cutoff identification k <-> r is ambiguous (Donoghue's critique).

In the tau framework:
- G(k) emerges because the gravitational channel N_grav(k) changes with scale
- The DPI provides a model-independent constraint on the running
- No cutoff identification ambiguity: the scale k is the scale at which QRE is evaluated

### 11.2 Unified Framework

Standard running G: Three separate regimes with separate physics.

Tau framework: ONE equation Sigma = D(rho_spacetime || rho_matter), different limits:
- Sigma = 0: flat spacetime, perfect recovery (Paper 1)
- Sigma ~ r_s/r: Schwarzschild (Paper 2)
- Sigma ~ r_s/r + ln correction: rotation curves (Paper 3)
- Sigma ~ volume-law contribution: cosmological (Paper 4?)

### 11.3 Dark Matter = Extra Information Loss

The reinterpretation: "dark matter" at galactic scales is NOT a substance. It is the EXTRA information loss (extra Sigma) from the IR running of the gravitational channel. When we compute how much the gravitational channel degrades quantum information at galactic scales, we get MORE degradation than at solar-system scales. This excess shows up as additional gravitational binding.

### 11.4 Testable Prediction: The tau-Sigma Connection

From Paper 1: tau = 1 - F, where F >= exp(-Sigma/2).

At galactic scales, the tau of the gravitational channel is:

```
tau_grav(r) >= 1 - exp(-Sigma_grav(r)/2)                                  (44)
```

For the Milky Way at r = 8 kpc:
- Sigma_Newton ~ r_s/r ~ 3 x 10^{-6}
- Sigma_running ~ alpha * ln(8/37) ~ 10^{-6}
- tau_grav ~ Sigma/2 ~ 2 x 10^{-6}

This means: at the solar radius, the gravitational channel has an irrecoverability of ~2 parts per million. This tiny but non-zero tau IS the origin of the galactic time arrow at this scale.

---

## 12. Experimental Predictions and Tests

### 12.1 Rotation Curve Shape

The tau framework predicts (from Eq. 9):

```
v^2(r) = G_N M_bar(r)/r + V_0^2                                          (45)
```

where V_0^2 = 2G_N M_bar k_*/pi is constant.

This predicts:
- Declining v(r) at intermediate radii (Keplerian falloff + constant)
- Asymptotically flat v(r) at large r
- NO upturn or downturn at very large r (unlike some DM halo models)

### 12.2 Universal Crossover Scale

k_* ~ 2.7 x 10^{-2} kpc^{-1} should be UNIVERSAL across all galaxies (if the running is truly from fundamental physics and not environmental). This predicts r_c ~ 37 kpc as a universal transition scale.

**Test with BIG-SPARC (4000 galaxies)**: Fit k_* for each galaxy and check the distribution. If k_* clusters tightly around 0.027 kpc^{-1} regardless of galaxy mass, this strongly supports the framework. If k_* varies systematically with mass, it would suggest environment-dependent running (still consistent but less elegant).

### 12.3 The RAR as a Diagnostic

The tau framework predicts the RAR through:

```
g_obs = g_bar + 2G_N M_bar k_* / (pi r^2)                                (46)
```

At low g_bar (large r):

```
g_obs ≈ V_0^2 / r = sqrt(4G_N M_bar k_* / pi) * g_bar^{1/2} * ...
```

This reproduces the deep-MOND behavior g_obs ~ sqrt(a_0 g_bar) if k_* is appropriately related to a_0.

### 12.4 Gravitational Lensing

Unlike MOND (which requires special relativistic extensions), the tau framework is built on GR. Light deflection follows from the modified metric:

```
ds^2 = -(1 - r_s/r - 2alpha ln(r/r_0)) c^2 dt^2 + (1 + r_s/r + 2alpha ln(r/r_0)) dr^2 + ...    (47)
```

Both time and space components are modified -> lensing is TWICE the Newtonian prediction (just like GR), even in the modified regime.

### 12.5 Wide Binary Stars

Wide binary stars with separations > 1 kpc probe the transition regime. Recent (controversial) claims by Chae (2024) of ~40% gravity boost at low acceleration are potentially consistent with the logarithmic correction onset.

---

## 13. Critical Analysis: Weakest Points of Each Approach

### 13.1 Kumar (Running G)

**WEAKEST POINT**: eta = 1 is argued from marginality, not derived. If eta = 0.8 or 1.2, the entire structure changes:
- eta < 1: correction dies off, no flat rotation curves
- eta > 1: correction grows unboundedly, violates solar system constraints at large enough r

The framework's predictive power rests entirely on eta being EXACTLY 1. This is a fine-tuning concern.

### 13.2 RGGR (Gubitosi et al.)

**WEAKEST POINT**: nu_bar is not universal. It varies over 11 orders of magnitude from solar system to massive ellipticals. Without a theoretical derivation of nu_bar(M), this is a fit parameter, not a prediction.

### 13.3 Verlinde

**WEAKEST POINT**: No covariant formulation exists. Hossenfelder showed that naive covariantization reduces to GR. Without a covariant theory, predictions beyond static spherical systems are impossible.

### 13.4 MOND

**WEAKEST POINT**: The relativistic extension (TeVeS) predicted v_GW != c, which was decisively ruled out by GW170817/GRB 170817A. Alternative relativistic completions exist but are ad hoc.

### 13.5 The tau Framework

**WEAKEST POINT**: Three layers of difficulty:
1. eta = 1 not derived from QRE
2. The channel N_grav(k) not rigorously defined in AQFT at galactic scales
3. CMB problem has no solution within the framework

The framework provides an elegant REINTERPRETATION but does not yet provide a DERIVATION from first principles. Paper 3 should be positioned as providing the reinterpretation + the logical chain connecting established results, NOT as a complete derivation.

---

## 14. Strategy for Paper 3

### 14.1 What to Claim

1. **Established synthesis (NEW)**: The Casini-Huerta + Dorau-Much + Kumar chain provides a coherent information-theoretic interpretation of running G at galactic scales.

2. **Sigma_grav(r) formula (NEW interpretation of known result)**: Eq. (2) gives a unified entropy production covering both strong-field and galactic regimes.

3. **Dark matter as information loss (NEW interpretation)**: The "missing mass" is reinterpreted as extra Sigma from IR running of the gravitational channel.

4. **DPI constraints on eta (NEW conjecture)**: Eq. (29) provides an information-theoretic motivation for eta <= 1.

5. **Connection to tau (NEW)**: Eq. (44) gives the irrecoverability parameter at galactic scales.

### 14.2 What NOT to Claim

1. DO NOT claim eta = 1 is derived from QRE
2. DO NOT claim the CMB problem is solved
3. DO NOT claim superiority over LCDM at all scales
4. DO NOT claim the Bullet Cluster is explained
5. DO NOT claim novelty for running G itself (this is Reuter 2004, Kumar 2025)

### 14.3 Paper Structure (Proposed)

1. Introduction: The dark matter problem from the information perspective
2. Review: Sigma_grav and the gravitational channel (from Papers 1 & 2)
3. Scale-dependent channel: N_grav(k) and running G from DPI
4. Galactic rotation curves: Sigma with logarithmic correction
5. Comparison: Kumar, Verlinde, MOND, LCDM
6. Tests: SPARC fitting, RAR, BTFR predictions
7. Limitations: CMB, clusters, structure formation
8. Discussion: Dark matter as information loss
9. Conclusion

### 14.4 Target Venue

PRL (4 pages) or PRD Rapid Communication (5 pages). The novelty is the SYNTHESIS, not any individual calculation. This is a "Rosetta Stone" paper -- translating between QI language and galaxy dynamics.

---

## 15. Complete Reference List

### Primary (Must Cite)

1. Kumar, N. (2025). "Marginal IR running of Gravity as a Natural Explanation for Dark Matter." Phys. Lett. B 871, 140008. arXiv:2509.05246.
2. Gubitosi, G., Piattella, O.F. & Casarini, L. (2024). "Phenomenology of renormalization group improved gravity from the kinematics of SPARC galaxies." Phys. Rev. D 110, 124014. arXiv:2403.00531.
3. Verlinde, E.P. (2016). "Emergent Gravity and the Dark Universe." SciPost Phys. 2, 016. arXiv:1611.02269.
4. Dorau, P. & Much, A. (2025/2026). "From Quantum Relative Entropy to the Semiclassical Einstein Equations." Phys. Rev. Lett. 136, 091602. arXiv:2510.24491.
5. Casini, H., Teste, E. & Torroba, G. (2017). "Relative Entropy and the RG Flow." JHEP 03, 089. arXiv:1611.00016.
6. Bianconi, G. (2025). "Gravity from Entropy." Phys. Rev. D 111, 066001. arXiv:2408.14391.
7. Lelli, F., McGaugh, S.S. & Schombert, J.M. (2016). SPARC. AJ 152, 157. arXiv:1606.09251.
8. McGaugh, S.S., Lelli, F. & Schombert, J.M. (2016). RAR. PRL 117, 201101. arXiv:1609.05917.
9. Reuter, M. & Weyer, H. (2004). "Running Newton Constant, Improved Gravitational Actions, and Galaxy Rotation Curves." Phys. Rev. D 70, 124028. arXiv:hep-th/0410117.

### Supporting (Should Cite)

10. Ghari, A. & Haghi, H. (2026). arXiv:2601.01715. [Verlinde > MOND 5.2 sigma]
11. Jacobson, T. (1995). PRL 75, 1260. [Thermodynamics of spacetime]
12. Jacobson, T. (2015). PRL 116, 201101. arXiv:1505.04753. [Entanglement equilibrium]
13. Donoghue, J.F. (1994). PRD 50, 3874. arXiv:gr-qc/9405057. [GR as EFT]
14. Bjerrum-Bohr, N.E.J., Donoghue, J.F. & Holstein, B.R. (2003). PRD 67, 084033. arXiv:hep-th/0211072. [Quantum corrections to Newton]
15. de Paula Netto, T., Modesto, L. & Shapiro, I.L. (2022). EPJC 82, 160. arXiv:2110.14263. [Universal log correction]
16. Donoghue, J.F. (2020). Front. Phys. 8, 56. arXiv:1911.02967. [Critique of AS]
17. Smolin, L. (2017). PRD 96, 083523. arXiv:1704.00780. [MOND from quantum gravity]
18. Li, P. et al. (2018). arXiv:1803.00022. [RAR fitting]
19. Haubner, M., Lelli, F. et al. (2024). arXiv:2411.13329. [BIG-SPARC]
20. Brouwer, M.M. et al. (2021). A&A 650. arXiv:2106.11677. [KiDS-1000 RAR]
21. Casini, H. & Huerta, M. (2012). PRD 85, 125016. arXiv:1202.5650. [Entropic c-function]
22. Deddo, E. et al. (2024). JHEP 09, 179. arXiv:2404.15077. [Entropic irreversibility theorems]
23. Faulkner, T., Lewkowycz, A. & Maldacena, J. (2013). JHEP 11, 074. arXiv:1307.2892. [FLM formula]
24. Grumiller, D. (2010). PRL 105, 211303. arXiv:1011.3625. [Rindler acceleration model]
25. Gillot, J. (2025). arXiv:2507.11524. [Quantum modified inertia]
26. Penner, A.R. (2024). EPJC. arXiv:2403.13019. [Many-body gravity]
27. Chakraborty, R. & Sengupta, S. (2025). arXiv:2505.04060. [Varying G galactic metrics]
28. Zholdasbek, A. et al. (2025). PRD 111, 103519. arXiv:2405.02636. [Running G + CMB]
29. Bianconi, G. (2025). Entropy 27, 266. arXiv:2501.09491. [QRE of Schwarzschild]
30. Bianconi, G. (2025). arXiv:2510.22545. [Thermodynamics of GfE]
31. Milgrom, M. (2020). arXiv:2001.09729. [a_0-cosmology connection]

### Contextual (Cite If Space)

32. Rodrigues, D.C. et al. (2009). arXiv:0911.4967. [RGGR original]
33. Hossenfelder, S. (2017). PRD 95. arXiv:1703.01415. [Verlinde covariant attempt]
34. Dai, D.-C. & Stojkovic, D. (2017). JHEP. arXiv:1710.00946. [Verlinde inconsistencies]
35. Oppenheim, J. & Russo, A. (2024). arXiv:2402.19459. [Stochastic gravity -- RETRACTED claim]
36. Rostami et al. (2025). arXiv:2511.05632. [Relativistic MOND from entropic gravity]
37. arXiv:2505.03061. [Tsallis-Renyi to MOND]

---

## 16. Summary of Key Equations

For quick reference, the essential formulas of the tau framework at galactic scales:

```
[Paper 1]  tau = 1 - F,   F >= exp(-Sigma/2)                       (Paper 1 bound)

[Paper 2]  Sigma_grav = -ln(g_00) = r_s/r                          (strong field)

[Paper 3]  Sigma_grav(r) = r_s/r + (2k_* r_s/pi) ln(r/r_0)        (galactic scale)

           v^2(r) = G_N M/r + 2G_N M k_*/pi                        (rotation curve)

           G(k) = G_N (k_*/k)^eta,  eta = 1                        (IR running)

           k_* ≈ 2.7 x 10^{-2} kpc^{-1}  (universal)              (crossover scale)

           a_0 ~ cH_0/(4pi)                                         (acceleration scale)

[Verlinde] g_obs ~ sqrt(a_0 g_bar / 6)   (deep MOND regime)

[DPI]      Sigma(k_IR) >= Sigma(k_UV) >= 0                          (monotonicity)

[Conjecture] eta <= 1 - exp(-Sigma/2) = tau_grav                   (DPI bound on eta)
```

---

## 17. Open Questions for Future Work

### Tier 1 (Must address in Paper 3)

1. **Derive or strongly motivate eta = 1** from information-theoretic principles
2. **Address the CMB problem** honestly -- what can and cannot be explained
3. **Quantify the SPARC fit quality** for the tau framework (at least compare to published RGGR and NFW results)

### Tier 2 (Should address in Paper 3 or 4)

4. **Cluster-scale predictions**: What does Sigma predict for galaxy cluster mass profiles?
5. **External Field Effect**: Does the scale-dependent channel produce EFE-like behavior?
6. **nu_bar derivation**: Can the mass dependence of nu_bar be derived from Sigma?

### Tier 3 (Paper 4 and beyond)

7. **Bianconi extension**: Full GQRE framework with running G at all scales
8. **CMB power spectrum**: Can Sigma at z ~ 1100 reproduce acoustic peaks?
9. **Structure formation**: Can the framework drive perturbation growth?
10. **Gravitational wave propagation**: Running G effects on GW dispersion?

---

*Last updated: 2026-03-11*
*Part of the four-paper series: Sigma = D(rho_spacetime || rho_matter)*
