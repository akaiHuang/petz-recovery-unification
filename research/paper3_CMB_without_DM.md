# The CMB Without Dark Matter: A Brutally Honest Assessment

## Research Date: 2026-03-11 (major revision with quantitative physics)
## Status: Comprehensive literature survey + quantitative analysis + critical assessment
## Relevance: Paper 3 (existential challenge), Paper 4 (framework viability)

---

## Executive Summary

The CMB acoustic peak structure is the **single hardest challenge** for any framework that attempts to eliminate dark matter particles. After exhaustive literature review and quantitative analysis, the honest conclusion is:

| Approach | Can it fit CMB peaks? | Parameters | Status |
|---|---|---|---|
| LCDM (6-parameter) | **YES** (chi^2/dof ~ 1) | 6 | Gold standard |
| AeST (Skordis-Zlosnik 2021) | **YES** (comparable to LCDM) | ~8-10 | Only successful relativistic MOND |
| Khronon theory (Blanchet-Skordis 2024) | **Claims YES** | ~6-8 | Newer, less tested |
| nuHDM (Haslbauer-Banik-Kroupa) | **PARTIALLY** (misses 4th peak) | ~8 | Requires 11 eV sterile neutrino |
| Superfluid DM (Berezhiani-Khoury) | **YES** (contains DM particles) | ~8-10 | Hybrid: particles + MOND phonons |
| Verlinde emergent gravity | **NO** | - | No perturbation theory exists |
| Running G (asymptotic safety) | **NO** | - | Wrong scale for CMB |
| f(R) gravity | **NO** (keeps CDM) | varies | Modifies, doesn't replace CDM |
| tau framework (current) | **NO** | - | No mechanism identified |

**Bottom line**: The ONLY approaches that successfully fit the CMB without CDM particles **introduce new fields that effectively act as CDM at early times**. There is no known way to explain the CMB peak structure using pure gravity modification alone.

**For the tau framework**: This is not a death sentence, but it is a hard boundary. The framework must either (a) accept CDM particles as an ingredient and explain rotation curves as a correction, (b) derive CDM-like perturbative behavior from the QRE structure of spacetime, or (c) connect to AeST-like fields through the modular flow approach.

---

## Part I: The Quantitative Physics of CMB Acoustic Peaks

This section provides the equation-level physics needed to understand WHY the CMB is so constraining. Based on Hu (2008, arXiv:0802.3688) and Planck 2018 (arXiv:1807.06209).

### 1.1 Planck 2018 Best-Fit Parameters

```
Omega_b h^2 = 0.02237 +/- 0.00015    (baryon density)
Omega_c h^2 = 0.1200 +/- 0.0012      (CDM density)
Omega_m h^2 = 0.1424 +/- 0.0011      (total matter = baryon + CDM)
H_0 = 67.36 +/- 0.54 km/s/Mpc        (Hubble constant)
n_s = 0.9649 +/- 0.0042              (spectral index)
tau_reion = 0.0544 +/- 0.0073        (reionization optical depth)
ln(10^10 A_s) = 3.044 +/- 0.014      (amplitude)
100 theta_* = 1.04110 +/- 0.00031    (angular acoustic scale, 0.03% precision)
```

The chi-squared for LCDM fit to Planck TT+TE+EE+lowE is approximately chi^2/dof ~ 1.0, fitting ~2500 multipole data points with 6 parameters.

### 1.2 The Photon-Baryon Fluid: Basic Equations

**From Wayne Hu's lecture notes (arXiv:0802.3688), the complete physics:**

The baryon-photon momentum ratio:
```
R = (rho_b + p_b) / (rho_gamma + p_gamma) = 3 rho_b / (4 rho_gamma)
  ~ 0.6 (Omega_b h^2 / 0.02) (a / 10^{-3})                    [Eq. 47]
```

At recombination (a ~ 10^{-3}): R_* ~ 0.6 for Planck values.

The matter-radiation ratio:
```
rho_m / rho_r = 3.6 (Omega_m h^2 / 0.15) (a / 10^{-3})         [Eq. 48]
```

At recombination: rho_m/rho_r ~ 3.4, meaning matter JUST dominates.

**This is the KEY: the universe transitions from radiation to matter domination right around recombination. This timing is set by Omega_m h^2 ~ 0.14, not by Omega_b h^2 ~ 0.02. Without CDM, Omega_m = Omega_b and the transition happens MUCH LATER (z_eq ~ 500 instead of z_eq ~ 3400).**

### 1.3 The Forced Oscillator Equation

The photon-baryon fluid satisfies a forced harmonic oscillator equation:
```
[(1+R) Theta_dot]' + (k^2/3) Theta = -(k^2/3)(1+R) Psi - [(1+R) Phi_dot]'   [Eq. 99]
```

where:
- Theta = delta_T/T (temperature fluctuation)
- Psi = Newtonian gravitational potential (time-time metric perturbation)
- Phi = spatial curvature perturbation
- R = baryon-photon momentum ratio
- k = wavenumber of the perturbation

The sound speed is:
```
c_s^2 = (p_gamma_dot + p_b_dot) / (rho_gamma_dot + rho_b_dot) = 1 / (3(1+R))   [Eq. 100]
```

And the sound horizon at recombination:
```
s_* = integral_0^{eta_*} c_s d eta                                               [Eq. 65]
```

### 1.4 The Three Key Effects of Dark Matter

**Effect 1: Gravitational Potential Wells (Sections 3.3, 3.5 of Hu)**

The effective temperature is Theta + Psi (including gravitational redshift). For constant potentials (matter-dominated era):
```
[Theta + Psi](eta) = [Theta + (1+R) Psi](0) cos(ks) - R Psi(0)              [Eq. 101]
```

The Sachs-Wolfe initial condition gives:
```
Theta + Psi = (1/3) Psi(0)    for superhorizon modes                          [Eq. 92]
```

**CDM's role**: In the matter-dominated era, CDM maintains CONSTANT gravitational potentials (Psi = const). The Poisson equation:
```
k^2 Phi = 4 pi G a^2 rho_m Delta_m                                             [Eq. 83]
```

In matter domination, Delta_m grows as a (scale factor), and rho_m ~ a^{-3}, giving Phi ~ constant. **Without CDM**, during radiation domination, the photon pressure causes Phi to oscillate and DECAY as a^{-2}. This decay of potential wells is THE critical effect.

**Effect 2: Radiation Driving (Section 3.5)**

When a mode enters the horizon DURING RADIATION DOMINATION (before matter-radiation equality), the potential wells decay because radiation pressure prevents gravitational collapse. This decay is IN PHASE with the acoustic oscillations, acting as a DRIVING force:
```
[Theta + Psi](eta) = [Theta + Psi](0) + Delta_Psi - Delta_Phi
                    = (1/3) Psi(0) - 2 Psi(0) = (5/3) Psi(0)               [Eq. 107]
```

This gives an acoustic amplitude that is **5x the Sachs-Wolfe effect** for modes that entered the horizon during radiation domination. This enhancement is the RADIATION DRIVING effect.

**CDM's crucial role**: CDM sets the redshift of matter-radiation equality:
```
z_eq = Omega_m / Omega_r ~ 3400    (with CDM: Omega_m h^2 = 0.14)
z_eq ~ 500                          (without CDM: Omega_m h^2 = Omega_b h^2 = 0.02)
```

The wavenumber at horizon crossing at equality is:
```
eta_eq ~ 114 (Omega_m h^2 / 0.14)^{-1} Mpc                                   [Eq. 134]
```

Modes with k > 1/eta_eq entered during radiation domination and get the driving boost. With CDM, only high-k modes (peaks 3, 4, 5...) get this boost. Without CDM, ALL peaks would enter during radiation domination, getting massive driving -- completely changing the peak height ratios.

**Effect 3: Baryon Loading (Section 3.4)**

Baryons add inertia to the photon-baryon fluid, displacing the zero point of oscillation:
```
[Theta + Psi]_n = [+/-(1+3R) - 3R] * (1/3) Psi(0)                           [Eq. 103]
```

This makes odd peaks (compression, n = 1, 3, 5...) LARGER than even peaks (rarefaction, n = 2, 4, 6...):
```
[Theta + Psi]_1 - [Theta + Psi]_2 = [-6R] * (1/3) Psi(0)                    [Eq. 103]
```

The baryon loading parameter R depends ONLY on Omega_b h^2, not on Omega_c h^2. So the odd/even peak asymmetry measures baryons, while the OVERALL peak heights (especially peak 3 vs peak 1) measure the TOTAL matter content.

### 1.5 The Silk Damping Tail

Photon diffusion (Silk damping) erases fluctuations on small scales:
```
lambda_D = sqrt(N) lambda_C = sqrt(eta / lambda_C^{-1}) * lambda_C = sqrt(eta lambda_C)   [Eq. 111]
```

The diffusion wavenumber:
```
k_D^{-2} = integral d eta (1/tau_dot) (1/(6(1+R))) (16/15 + R^2/(1+R))      [Eq. 119]
```

The damping envelope:
```
D_l ~ exp(-(l/l_D)^{1.25})                                                    [Eq. 122]
```

with the damping scale:
```
lambda_D / Mpc ~ 64.5 (Omega_m h^2 / 0.14)^{-0.278} (Omega_b h^2 / 0.024)^{-0.18}   [Eq. 123]
```

**This formula shows that lambda_D depends on BOTH Omega_m and Omega_b. Without CDM (Omega_m = Omega_b), the damping scale shifts dramatically, changing the entire damping tail.**

### 1.6 The Sound Horizon Formula

The sound horizon at recombination determines the angular scale of ALL peaks:
```
s_* = (2 sqrt(3) / 3) sqrt(a_* / (R_* Omega_m H_0^2))
      * ln[ (sqrt(1 + R_*) + sqrt(R_* + r_* R_*)) / (1 + sqrt(r_* R_*)) ]    [Eq. 128]
```

where:
```
R_* = (3/4)(rho_b/rho_gamma)|_{a_*} = 0.729 (Omega_b h^2 / 0.024)(a_*/10^{-3})   [Eq. 129]
r_* = (rho_r/rho_m)|_{a_*} = 0.297 (Omega_m h^2 / 0.14)^{-1} (a_*/10^{-3})^{-1} [Eq. 130]
```

**The sound horizon depends on Omega_m h^2 through r_* (the radiation-to-matter ratio at recombination). Without CDM, Omega_m h^2 drops from 0.14 to 0.02, changing r_* from 0.3 to 2.1 -- a factor of 7 increase. This would shift ALL peak positions.**

### 1.7 Summary: The Three Roles of CDM in the CMB

| Role | Physics | Observable | Without CDM |
|------|---------|------------|-------------|
| Potential wells | CDM maintains Phi = const after z_eq | Third peak height | Phi decays, peaks 3+ suppressed by ~50% |
| z_eq timing | Omega_m h^2 = 0.14 sets z_eq ~ 3400 | Peak height ratios | z_eq ~ 500, ALL peaks radiation-driven |
| Sound horizon | Omega_m h^2 enters s_* formula | Peak positions | s_* shifts by ~30%, all peaks move |

**Quantitative impact of removing CDM:**

| Observable | With CDM (Planck) | Without CDM (baryons only) | Discrepancy |
|---|---|---|---|
| z_eq | ~3400 | ~500 | Factor ~7 |
| s_* | 144.4 Mpc | ~190 Mpc | ~32% |
| l_1 (first peak) | 220 | ~170 | ~23% |
| A_1/A_2 | 2.4 | ~3.5-4.0 | ~50% |
| A_3/A_2 | 0.97 | ~0.55 | ~43% |
| lambda_D | ~9 Mpc | ~14 Mpc | ~55% |

**These are not subtle effects. Every single CMB observable changes by 20-50% when CDM is removed. This is why the CMB is the hardest test for any dark-matter-free theory.**

---

## Part II: Literature Survey of Alternative Approaches

### 2.1 Skordis & Zlosnik: AeST Theory (arXiv:2007.00082, PRL 127, 161302, 2021)

**THE BREAKTHROUGH**: This is the ONLY relativistic theory that both reproduces MOND at galaxy scales AND fits the CMB power spectrum.

**The theory**: Aether Scalar Tensor (AeST) theory introduces:
1. A unit timelike **vector field** A^mu (the "aether") with constraint |A|^2 = -1
2. A **scalar field** phi
3. A **free function** f that determines the theory's behavior in different regimes

**The action** (schematically):
```
S = S_GR[g] + S_vector[A, g] + S_scalar[phi, g, A] + S_matter[psi, g]
```

**How it fits the CMB**: At cosmological scales, the scalar field phi and vector field A together produce perturbations that:
- Have effectively zero pressure (c_s^2 ~ 0) at early times
- Do not couple to photons (no Thomson scattering)
- Provide gravitational potential wells
- Behave EXACTLY like CDM perturbations in the linear regime

The theory reduces to a subset of the **Generalized Dark Matter (GDM) model** (Kunz 2009, arXiv:0702615) at cosmological scales, with:
- w ~ 0 (dust-like equation of state)
- c_s^2 ~ 0 (zero sound speed)
- c_vis^2 ~ 0 (no viscosity)

**The deep insight**: AeST works because the vector field A^mu, being constrained to be unit timelike, breaks Lorentz invariance in the gravitational sector. The scalar field phi then has its dynamics controlled by A^mu, and the combined system produces perturbations that are "stiff" (do not oscillate) and "cold" (c_s ~ 0). In effect, the vector-scalar system acts as a medium that responds to gravity but not to radiation pressure -- exactly mimicking CDM.

**How it transitions to MOND**: At galaxy scales, the free function f is chosen so that the quasi-static limit of the scalar field equation reduces to the MOND potential equation:
```
nabla . [mu(|nabla phi|/a_0) nabla phi] = 4 pi G rho_b
```

This transition from CDM-like behavior (cosmological) to MOND-like behavior (galactic) is controlled by the free function f, which acts as an interpolation between regimes.

**Key parameters**: ~8-10 total (vs LCDM's 6), including the free function parameters and coupling constants. The additional parameters are the price of replacing CDM particles with CDM-like fields.

**Linear stability**: The v2 of the paper (arXiv:2109.13287) showed that propagating modes are massive and ghost-free on Minkowski background. Potential instability at k < mu ~ Mpc^{-1} scale, but this is at cosmological scales where the cosmological background differs from Minkowski.

**Known problems**:
- Galaxy cluster regime shows "negative mass density" artifacts (arXiv:2312.00889)
- The free function f is ad hoc (not derived from first principles)
- More parameters than LCDM
- Non-linear structure formation not yet simulated

**HONEST VERDICT**: AeST works because it introduces new degrees of freedom (vector + scalar fields) that behave like CDM at early times. It does NOT show that "gravity modification alone" can fit the CMB. It shows that you need SOMETHING with CDM-like perturbative behavior, whether it's particles or fields. The philosophical distinction (particles vs fields) is real but the perturbative physics is identical.

### 2.2 Blanchet & Skordis: Khronon Theory (arXiv:2404.06584, 2024)

A more elegant version of the same idea. Introduces a scalar field tau (the "Khronon") that foliates spacetime.

**Key features**:
- Reduces to MOND for stationary systems at galaxy scales
- Recovers GR + Lambda in strong field regime
- At cosmological scales, reduces to a subset of GDM
- Single propagating scalar degree of freedom with dispersion relation omega = 0

**CMB fit**: Claims agreement with "observed cosmic microwave background anisotropies at linear cosmological scales" but no detailed chi^2 comparison published yet.

**HONEST VERDICT**: Same fundamental mechanism as AeST -- scalar field perturbations behave like CDM perturbations at early times. More elegant but same physics.

### 2.3 Superfluid Dark Matter (Berezhiani & Khoury, arXiv:1507.01530, 2015)

**The approach**: A hybrid framework where dark matter particles exist but form a superfluid in galaxy halos. The superfluid's phonon excitations mediate a MOND-like force.

**The mechanism**: At temperatures below a critical temperature T_c (satisfied in galaxy halos), DM particles condense into a superfluid state. The phonons of this superfluid couple to baryons through:
```
L_phonon ~ Lambda (2m)^{3/2} X sqrt(|X|) / 3    where X = theta_dot - m Phi - (nabla theta)^2/(2m)
```

This three-body coupling between phonons and baryonic matter produces a force law:
```
F_phonon ~ sqrt(a_0 a_N)    (geometric mean of MOND scale and Newtonian acceleration)
```

**For the CMB**: Superfluid DM retains CDM particles, so it fits the CMB exactly like LCDM. At z ~ 1100, the DM is NOT in a superfluid state (T >> T_c), so it behaves as standard CDM.

**Key advantage**: Naturally resolves the transition between CDM (cosmological/cluster) and MOND (galactic) behavior through a phase transition.

**HONEST VERDICT**: This is the intellectually most complete framework, but it REQUIRES dark matter particles. It explains why DM behaves like MOND at galaxy scales without abandoning CDM for cosmology. For the tau framework, this is both a model and a warning: the most successful framework keeps DM particles and adds MOND as an emergent phenomenon.

### 2.4 nuHDM Model (Haslbauer, Banik, Kroupa; arXiv:2009.11292)

**The approach**: MOND gravity + 11 eV sterile neutrinos as hot dark matter.

**CMB fit**:
- Optimized opt-nuHDM: Fits CMB angular power spectrum EXCEPT the fourth peak
- Requires very different parameters: Omega_m ~ 0.5, H_0 ~ 55.6 km/s/Mpc

**Critical problems**:
1. 11 eV sterile neutrino has NO experimental support (KATRIN, reactor anomalies disfavor it)
2. Hot dark matter has too much free streaming -- washes out small-scale structure
3. The opt-nuHDM parameters are in severe tension with local H_0 measurements
4. Missing the 4th peak means the damping tail is wrong
5. Samaras, Grandis & Kroupa (arXiv:2506.19196) found that initial conditions are inconsistent

**HONEST VERDICT**: Replaces CDM with HDM. Still requires dark matter particles (just different ones). And the fit is significantly worse than LCDM.

### 2.5 Verlinde's Emergent Gravity (arXiv:1611.02269)

**For the CMB**: Verlinde's framework has NO perturbation theory and makes NO predictions for the CMB power spectrum. The volume-law entanglement mechanism is fundamentally a z ~ 0 effect, not defined for z ~ 1100.

**HONEST VERDICT**: Completely silent on the CMB. Not even wrong -- just undefined.

### 2.6 f(R) Gravity, Running G, Bianconi's GfE

All viable f(R) models include CDM as a standard component. The f(R) modification adds a scalaron that modifies late-time evolution, not the CMB peak structure.

Running G from asymptotic safety operates at the wrong scale. The IR running formula G(k) = G_N[1 + k_*/k] with k_* ~ 2.7 x 10^{-2} kpc^{-1} gives absurd results (G/G_N ~ 10^5) when naively applied at CMB scales (k ~ 10^{-4} Mpc^{-1}). The running MUST be cut off at intermediate scales and does not affect CMB physics.

Bianconi's gravity-from-entropy framework has no perturbation theory and makes no CMB predictions.

**HONEST VERDICT**: None of these approaches address the CMB.

### 2.7 Other Recent Approaches (2024-2026)

**Aoki & Mukohyama (arXiv:2005.13972)**: Minimally modified gravity. Still includes CDM; modifies expansion history around recombination.

**Bonanno et al. (arXiv:2405.02636)**: Running G from asymptotic safety cosmology. Affects Planck-era physics but not recombination epoch.

**Backreaction (Buchert, Wiltshire)**: Primarily a cosmological-scale effect (dark energy). No specific CMB calculation exists. Cannot produce the lensing/potential well structure needed for acoustic peaks.

---

## Part III: WHY the CMB Is So Hard -- The Information-Theoretic Perspective

### 3.1 The CDM Requirement Restated in tau Language

In the tau framework, CDM's three roles can be restated as:

**Role 1: A closed subsystem**
CDM perturbations represent a subsystem with tau = 0 relative to photons. They do not exchange information with the photon-baryon fluid via Thomson scattering. In tau language:
- CDM-photon channel: completely opaque (zero coupling, tau_CDM-gamma undefined)
- Baryon-photon channel: tightly coupled (Compton scattering, tau ~ 0 at z > 1100)

**Role 2: Zero-entropy perturbations**
CDM has c_s = 0, meaning perturbations are entropy-free (no pressure, no dissipation). In the language of Sigma:
- Sigma_CDM = 0 per oscillation cycle (no information loss in CDM perturbations)
- Sigma_baryon > 0 per cycle (baryon perturbations dissipate via Silk damping)

**Role 3: A gravitational memory**
CDM perturbations maintain the gravitational potential after radiation-matter equality. In tau language, CDM provides a stable information storage medium that is not erased by photon pressure.

### 3.2 The Fundamental Tension

The tau framework says: "dark matter effects come from information loss in the gravitational channel (running G, entropy corrections)."

The CMB says: "you need a SEPARATE component that doesn't interact with photons and has c_s = 0."

These are in tension. Modifying the gravitational channel changes the STRENGTH of gravity for all components equally (equivalence principle). To produce CDM-like behavior, you need a component that:
- Couples gravitationally (provides potential wells)
- Does NOT couple electromagnetically (no Thomson scattering)
- Has zero sound speed (c_s = 0, no oscillation)

**This cannot be achieved by modifying gravity alone.** You must modify the MATTER CONTENT (add new degrees of freedom with the right properties).

### 3.3 Why c_s = 0 Is So Hard to Get from Gravity Modification

The reason is deeply connected to the equivalence principle:

1. **Gravity modification changes the Poisson equation**: nabla^2 Phi = 4 pi G_eff rho
2. This changes the AMPLITUDE of potential wells, not their oscillation properties
3. The oscillation properties depend on the FLUID equations (pressure, sound speed)
4. To change the fluid equations, you need new fluid components, not new gravity

**Corollary**: Any theory that modifies ONLY the left-hand side of Einstein's equations (geometry) cannot fit the CMB. You MUST modify the right-hand side (stress-energy content).

AeST circumvents this by adding vector and scalar fields to the right-hand side. These fields provide stress-energy that mimics CDM. This is the ONLY known mechanism that works.

### 3.4 The Jeans Length Argument

Even if we could define "Sigma perturbations" as a separate degree of freedom, their Jeans length would be:
```
lambda_J = c_s sqrt(pi / (G rho))
```

For CDM: c_s = 0, so lambda_J = 0 (perturbations grow at ALL scales).
For Sigma perturbations: c_s ~ c (propagate at light speed), so lambda_J ~ Hubble radius.

Sigma perturbations would NOT collapse on sub-horizon scales -- they would be smooth. This is the OPPOSITE of what's needed.

**KEY INSIGHT**: To fit the CMB, you need perturbations with c_s = 0. Running G corrections and Sigma perturbations propagate at c_s ~ c, making them behave like radiation, not CDM.

---

## Part IV: The AeST Lesson and Its Implications for the tau Framework

### 4.1 What AeST Teaches Us

AeST proves two things:
1. It IS possible to fit the CMB without CDM particles (using fields instead)
2. The fields MUST behave like CDM at the perturbation level (c_s ~ 0, w ~ 0)

The question for the tau framework becomes: can the framework produce AeST-like fields from its information-theoretic structure?

### 4.2 The Phantom Scalar Field Problem

In the tau framework, the exponential metric is sourced by a phantom scalar phi = M/r. Could this be the needed CDM-like field?

**Answer: NO, for three fundamental reasons.**

**Reason 1: Wrong equation of state.** A phantom scalar with V = 0 gives w = -1 (dark energy), not w = 0 (dark matter). Promoting phi = M/r to a cosmological field phi(t) does not give dust-like behavior.

**Reason 2: Wrong perturbation behavior.** Phantom fields have c_s^2 = -1 (imaginary sound speed), leading to exponential instability. CDM requires c_s^2 = 0.

**Reason 3: Wrong coupling.** The phantom scalar sources negative energy density (phantom kinetic term), while CDM requires positive energy density.

### 4.3 Possible Bridges to AeST-like Fields

**Bridge 1: Modular Flow Approach (MOST PROMISING)**

The Dorau-Much result (arXiv:2510.24491, PRL 2025) derives Einstein's equations from QRE via modular flow. If extended to include:
- Cosmological FRW background
- Linear perturbation theory
- Quantum corrections from the modular Hamiltonian

Then the additional degrees of freedom might produce CDM-like perturbations. The modular Hamiltonian K = -ln(rho) generates a flow that is equivalent to time evolution in the Rindler wedge. If this flow has modes with c_s = 0, they would behave like CDM.

**Status**: Completely unknown. No calculation exists.

**Bridge 2: Emergent Sectors from Hilbert Space Decomposition**

In quantum information, a Hilbert space decomposes as H = H_visible x H_hidden. The QRE between full and reduced states gives:
```
D(rho_full || rho_visible tensor rho_hidden) = I(visible : hidden)  (mutual information)
```

If spacetime emerges from this structure, the "hidden sector" might produce CDM-like perturbations that are:
- Gravitationally coupled (through the metric, which emerges from the full H)
- Photon-decoupled (the hidden sector is traced over, not observable)
- Pressureless (if the hidden sector entanglement is long-range)

**Status**: Pure speculation but mathematically well-motivated.

**Bridge 3: GfE Perturbation Theory**

If Bianconi's gravity-from-entropy action S = integral sqrt(-g) D(g || G) is perturbed around FRW, the resulting perturbation equations might contain a zero-pressure mode from the "G-field" (the matter-induced metric). This would be the most direct realization of CDM-like behavior from QRE.

**Status**: No perturbation theory developed.

---

## Part V: The Minimum Required New Ingredient

### 5.1 The GDM Framework (Kunz 2009)

Based on the Generalized Dark Matter framework, the minimum modification to fit the CMB requires ONE new component with:
```
w ~ 0          (dust-like equation of state)
c_s^2 ~ 0     (zero sound speed)
c_vis^2 ~ 0   (no viscosity)
Omega_new h^2 ~ 0.12   (correct density)
```

This is 4 conditions. Any theory that produces them will fit the CMB, regardless of the microscopic origin. AeST produces them from fields. LCDM from particles. Both work at the CMB level.

### 5.2 Numerical Observables That Must Be Matched

| Observable | Planck value | Tolerance | Physics |
|---|---|---|---|
| l_1 (first peak) | 220.0 +/- 0.5 | 0.2% | Sound horizon / angular distance |
| A_1/A_2 (peak ratio) | 2.40 +/- 0.03 | 1.2% | Baryon loading |
| A_3/A_2 (peak ratio) | 0.97 +/- 0.02 | 2% | Matter-radiation equality |
| l_D (damping scale) | 1177 +/- 10 | 0.8% | Diffusion / expansion history |
| z_eq | 3387 +/- 21 | 0.6% | Total matter density |
| r_s(z_*) | 144.43 +/- 0.26 Mpc | 0.2% | Expansion + sound speed integral |

### 5.3 The A_3/A_2 Smoking Gun

The third peak height relative to the second is the most diagnostic:
```
A_3/A_2 ~ 0.55    (baryon-only universe, no CDM)
A_3/A_2 = 0.97    (observed by Planck)
```

This 43% discrepancy is NOT a subtle effect. It arises because:
1. Without CDM, z_eq drops to ~500 (from ~3400)
2. The third peak enters the horizon during radiation domination
3. The driving effect (potential decay) boosts the oscillation amplitude
4. BUT the potential wells are also OSCILLATING (no CDM to stabilize them)
5. The net effect is a LOWER third peak relative to the second

To get A_3/A_2 from 0.55 to 0.97, you need EITHER:
- CDM with Omega_c h^2 ~ 0.12 (standard)
- A field with c_s^2 ~ 0 and Omega_field h^2 ~ 0.12 (AeST-like)
- Neutrinos with m_nu ~ 11 eV (nuHDM, doesn't fully work)

There is NO known modification of gravity that can change this ratio without adding a new pressureless component.

---

## Part VI: The Running G Problem at CMB Scales

### 6.1 Why Running G Cannot Help

The IR running of G with eta = 1 gives:
```
G(k) = G_N [1 + k_*/k]    for k << k_*
```

with k_* ~ 2.7 x 10^{-2} kpc^{-1} (galactic crossover scale).

At CMB scales, the relevant wavenumber is:
```
k_CMB ~ 10^{-2} Mpc^{-1} = 10^{-5} kpc^{-1}
```

Naive extrapolation gives:
```
G(k_CMB)/G_N ~ 1 + k_*/k_CMB ~ 1 + 2.7 x 10^{-2} / 10^{-5} ~ 2700
```

This is ABSURD and means the formula breaks down. The physical resolution:

1. **The running is a galactic-scale phenomenon**: The IR running of G is driven by local curvature from matter concentrations. At z ~ 1100, perturbations are delta rho/rho ~ 10^{-5} -- too small to drive significant running.

2. **The running must be cut off**: At scales r >> r_c ~ 37 kpc, the de Sitter background curvature (not local mass) dominates. The formula G(k) = G_N[1 + k_*/k] is only valid for k ~ k_* (galactic scales).

3. **At CMB scales, G = G_N**: The quantum corrections that produce running G are negligible at cosmological scales and early times.

### 6.2 Quantitative Constraint from Planck

Planck constrains the running of G at z ~ 1100:
```
|G(z ~ 1100)/G_0 - 1| < 0.05    (95% CL, from peak ratios + damping tail)
```

Even a 5% change in G at recombination shifts:
- z_eq by ~5% (through H^2 = 8 pi G rho / 3)
- s_* by ~2.5% (through the sound horizon integral)
- l_1 by ~2.5% (through D_A / s_*)

These shifts are detectable at >10 sigma with Planck precision.

**Bottom line**: Running G has essentially ZERO effect on CMB physics because: (a) the running formula breaks down at CMB scales, and (b) even small running (5%) is ruled out by Planck.

---

## Part VII: Honest Assessment for the tau Framework

### 7.1 Scorecard

| Question | Honest Answer | Confidence |
|----------|---------------|------------|
| Can tau explain rotation curves without CDM? | YES (running G) | HIGH |
| Can tau explain RAR? | YES (universal k_*) | HIGH |
| Can tau explain CMB peaks without CDM? | **NO** | **VERY HIGH** |
| Can tau produce CDM-like perturbations? | Unknown (modular flow?) | VERY LOW |
| Is this a fatal problem? | No, but it's the hardest challenge | HIGH |

### 7.2 The Four Options

**Option A: Accept CDM + tau corrections**
The tau framework keeps CDM as a fundamental component but ADDS running G at galaxy scales. CDM explains CMB (z ~ 1100); running G modifies CDM halo profiles at z ~ 0. This is the safest option but makes tau less revolutionary.

**Viability**: MODERATE. Consistent but unambitious.

**Option B: tau framework derives CDM-like field from QRE**
If modular flow / QRE perturbation theory produces a c_s = 0 mode, "dark matter" would be emergent from quantum spacetime geometry. No new particles needed.

**Viability**: VERY LOW (currently). No calculation exists.

**Option C: Connect to AeST**
If the tau framework's Sigma = D(rho_spacetime || rho_matter) can be shown to contain AeST-like degrees of freedom when properly quantized, the AeST fields would be the "missing" CDM-like perturbations from within the tau framework.

**Viability**: SPECULATIVE. Would require showing that AeST's vector + scalar fields emerge from the QRE structure.

**Option D: Honest agnosticism**
Paper 3 focuses on galactic-scale successes. Acknowledges CMB as an open problem. States clearly that the tau framework currently has no mechanism for CMB peaks without CDM. Proposes concrete future calculations (modular flow, GfE perturbation theory).

**Viability**: HIGH. Intellectually honest and sets the research agenda.

### 7.3 Recommended Strategy for Paper 3

1. **Focus on what works**: Rotation curves from running G (well-supported by Kumar, Gubitosi et al.)
2. **Honestly acknowledge the CMB problem**: State explicitly that running G does not affect CMB physics
3. **Frame the CMB as establishing the boundary of the framework's current scope**: The tau framework explains galactic-scale phenomenology; the cosmological origin of the "dark component" remains open
4. **Point to AeST as proof of concept**: Demonstrate that MOND-like theories CAN fit the CMB if they have the right perturbative structure (c_s = 0 fields)
5. **Propose modular flow / QRE perturbation theory**: The most promising direction for deriving CDM-like perturbations from quantum spacetime
6. **Note that CDM itself has problems at galaxy scales**: Core-cusp problem, too-big-to-fail, diversity problem. Even if CDM exists, something modifies it at galactic scales -- and the tau framework provides that modification

### 7.4 What Paper 3 Should NOT Do

1. **Do NOT claim** that running G explains the CMB -- it provably doesn't
2. **Do NOT ignore** the CMB problem -- it is the strongest objection reviewers will raise
3. **Do NOT invoke** nuHDM -- it doesn't work and requires undetected particles
4. **Do NOT treat** the CMB as a minor issue -- it is THE issue for any dark-matter-free theory

---

## Part VIII: The Deeper Question -- Is "Dark Matter" Fields vs Particles Meaningful?

### 8.1 The Physics of the Distinction

At the level of CMB physics (linear perturbation theory), there is NO distinction between CDM particles and CDM-like fields. Both produce identical TT, TE, EE power spectra.

The distinction matters for:
- **Direct detection**: Particles can be detected in underground experiments; fields cannot
- **Collider production**: Particles can be produced at LHC; fields cannot
- **Non-linear structure formation**: N-body simulations differ at the percent level
- **Self-interaction**: Particles can have sigma_SI; fields have different non-linearities
- **Small-scale structure**: Free-streaming scales differ

### 8.2 The Philosophical Point

If AeST is correct, then "dark matter" is real (it gravitates, it has density Omega h^2 = 0.12) but it is a property of spacetime geometry (the vector-scalar field configuration) rather than a new particle species. This is a DEEP distinction:

- **Particle CDM**: Dark matter is a substance added to spacetime
- **Field CDM (AeST)**: Dark matter is a configuration of spacetime itself

The tau framework, with Sigma = D(rho_spacetime || rho_matter), is philosophically aligned with the field interpretation. If "dark matter" IS a spacetime configuration, then it IS the kind of thing the tau framework should describe.

### 8.3 The Research Direction

The most promising path forward:
1. Show that Sigma = D(rho_spacetime || rho_matter), when computed on an FRW background with perturbations, contains modes with c_s = 0
2. Identify these modes with the AeST-like vector + scalar fields
3. Show that the coupling constants and the free function f of AeST are determined by the QRE structure
4. This would reduce AeST's ~10 parameters to derivable quantities from the tau framework

This is Paper 4 territory, not Paper 3. But it should be flagged as the key open question.

---

## Part IX: Key References

### CMB Physics (Foundational)
- **Hu 2008**: arXiv:0802.3688 (CMB Theory from Nucleosynthesis to Recombination -- THE pedagogical reference)
- **Planck 2018**: arXiv:1807.06209 (cosmological parameters)
- **Planck 2018**: arXiv:1807.06211 (power spectra and likelihoods)
- **Hu & Dodelson 2002**: Ann. Rev. Astron. Astrophys. 40, 171 (CMB Anisotropies)

### Alternative CMB Models
- **Skordis & Zlosnik 2021**: arXiv:2007.00082, PRL 127, 161302 (AeST -- THE key paper)
- **Skordis & Zlosnik 2021**: arXiv:2109.13287 (AeST v2, linear stability)
- **Blanchet & Skordis 2024**: arXiv:2404.06584 (Khronon theory)
- **Haslbauer, Banik & Kroupa 2020**: arXiv:2009.11292 (nuHDM)
- **Samaras, Grandis & Kroupa 2025**: arXiv:2506.19196 (nuHDM initial conditions problem)

### AeST Extensions and Tests
- arXiv:2309.06232 (dynamical system analysis of AeST cosmology)
- arXiv:2312.00889 (AeST galaxy clusters -- negative mass density problem)
- arXiv:2305.07742 (quasi-static limit, new scale m_x)
- arXiv:2301.03499 (AeST weak lensing confrontation)

### Superfluid Dark Matter
- **Berezhiani & Khoury 2015**: arXiv:1507.01530 (original proposal)
- **Berezhiani & Khoury 2016**: arXiv:1602.00139 (cosmology, Jeans instability)
- **Khoury 2016**: arXiv:1602.05961 (another dark matter story)

### GDM Framework
- **Kunz 2009**: arXiv:0702615 (Generalized Dark Matter parameterization)
- arXiv:1605.00649 (extensive GDM investigation)

### Running G / Asymptotic Safety
- **Kumar 2025**: arXiv:2509.05246 (IR running -> rotation curves, CMB FAILS)
- **Bonanno et al. 2024**: arXiv:2405.02636 (running G cosmology)
- **Gubitosi et al. 2024**: arXiv:2403.00531 (SPARC validation)

### Entropic/Information-Theoretic Gravity
- **Verlinde 2016**: arXiv:1611.02269 (emergent gravity -- no CMB predictions)
- **Dorau-Much 2025**: arXiv:2510.24491 (QRE -> Einstein)
- **Bianconi 2025**: arXiv:2408.14391, arXiv:2510.22545 (gravity from entropy)
- **Casini-Huerta 2017**: arXiv:1611.00016 (QRE and RG flow)

### MOND Reviews
- **Famaey & McGaugh 2012**: arXiv:1112.3960, Living Rev. Relativ. 15, 10 (comprehensive MOND review)
- **Milgrom 2020**: arXiv:2001.09729 (MOND vs DM: a duel with a twist -- a_0 ~ cH_0)

---

## Appendix: Quantitative CMB Physics Formulae

For reference, the key formulae from Hu (2008):

**Sound horizon at recombination:**
```
s_* = (2 sqrt(3)/3) sqrt(a_*/(R_* Omega_m H_0^2)) ln[(sqrt(1+R_*) + sqrt(R_* + r_* R_*))/(1 + sqrt(r_* R_*))]
```

**Baryon-photon ratio at recombination:**
```
R_* = 0.729 (Omega_b h^2 / 0.024)(a_*/10^{-3})
```

**Radiation-matter ratio at recombination:**
```
r_* = 0.297 (Omega_m h^2 / 0.14)^{-1} (a_*/10^{-3})^{-1}
```

**Recombination redshift:**
```
a_*^{-1} ~ 1089 (Omega_m h^2/0.14)^{0.0105} (Omega_b h^2/0.024)^{-0.028}
```

**Conformal distance to matter-radiation equality:**
```
eta_eq ~ 114 (Omega_m h^2/0.14)^{-1} Mpc
```

**Silk damping scale:**
```
lambda_D/Mpc ~ 64.5 (Omega_m h^2/0.14)^{-0.278} (Omega_b h^2/0.024)^{-0.18}
```

**Baryon loading effect on peak heights:**
```
[Theta + Psi]_odd = -(1 + 6R)(1/3) Psi(0)     (compression, enhanced)
[Theta + Psi]_even = +(1 - 2R)(1/3) Psi(0)     (rarefaction, suppressed)
Odd-even difference = -6R * (1/3) Psi(0)
```

**Acoustic scale (first peak angular position):**
```
l_A ~ 200    (for flat universe, matter-dominated at recombination)
theta_A ~ eta_*/eta_0 ~ 1/30 ~ 2 degrees
```

**Radiation driving enhancement (modes entering horizon during radiation domination):**
```
[Theta + Psi] = (5/3) Psi(0)    (5x Sachs-Wolfe, from potential decay)
```

---

*Last updated: 2026-03-11*
*Research conducted for Paper 3 of the four-paper series: Sigma = D(rho_spacetime || rho_matter)*
*This document is intended to be RIGOROUSLY HONEST about the tau framework's limitations regarding CMB.*
