# Bridge Research: From Paper 2 (Strong Field) to Paper 3 (Weak Field)

## Research Date: 2026-03-10
## Goal: Connect Sigma_grav = r_s/r (strong field) to Sigma_grav = r_s/r + alpha*ln(r) (weak field)

---

## Executive Summary

This document surveys the literature connecting quantum relative entropy (QRE) to weak-field gravity, running Newton's constant G, and galactic dynamics. The central question is whether the same framework Sigma = D(rho_spacetime || rho_matter) that produces r_s/r at strong field can naturally produce r_s/r + alpha*ln(r) at weak field (galactic scales).

**Key findings:**
1. The logarithmic correction arises naturally from marginal IR running of G with anomalous dimension eta=1 (Kumar 2025)
2. The Casini-Huerta theorem proves QRE monotonically decreases under RG flow, which constrains running couplings -- but has NOT been directly applied to G
3. The DPI for QRE provides the mathematical structure to constrain eta, but this connection is NOVEL and must be constructed
4. The CMB problem remains serious: running G alone cannot reproduce acoustic peaks
5. The most promising route is: QRE --> Casini-Huerta monotonicity --> constrained running of G --> logarithmic correction at galactic scales

---

## 1. COMPLETE BIBLIOGRAPHY

### 1.1 Running Newton's Constant: Asymptotic Safety and Functional RG

**[AS1]** Reuter, M. & Saueressig, F.
"Quantum Einstein Gravity"
New J. Phys. 14 (2012) 055022
- Review of asymptotic safety program
- Running G(k) near UV fixed point: G(k) -> G_* k^{-2} as k -> infinity
- Key result: G_k runs with energy scale via Wetterich equation

**[AS2]** Reuter, M. & Saueressig, F.
"Quantum Gravity and the Functional Renormalization Group: The Road towards Asymptotic Safety"
Cambridge University Press, 2019
- Comprehensive monograph on asymptotic safety
- Full treatment of Wetterich equation applied to gravity
- UV fixed point: g_* = G_* k^2 = finite

**[AS3]** Reuter, M. & Weyer, H.
"Running Newton Constant, Improved Gravitational Actions, and Galaxy Rotation Curves"
Phys. Rev. D 70 (2004) 124028 [arXiv:hep-th/0410117]
- **KEY PAPER**: First application of running G to galaxy rotation curves
- RG improvement of Einstein-Hilbert action promotes G and Lambda to spacetime scalars
- Result: G grows at large (astrophysical) distances in the IR
- Modified Newtonian limit reproduces flat rotation curves without dark matter
- Cutoff identification k <-> r maps momentum running to position-space running

**[AS4]** Koch, B. & Ramirez, I.
"Exact Renormalization Group with Optimal Scale and its Application to Cosmology"
Class. Quantum Grav. 28 (2011) 055008 [arXiv:1010.2799]
- Derives consistency condition for local form of cutoff scale
- Applied to homogeneous cosmology with UV fixed point
- Running couplings with ultraviolet fixed point

**[AS5]** Donoghue, J.F.
"A Critique of the Asymptotic Safety Program"
Front. Phys. 8 (2020) 56 [arXiv:1911.02967]
- **IMPORTANT CRITIQUE**: No useful and universal definition of running G(E) below Planck scale
- Power-law running in AS conflicts with perturbative calculations
- The running of Lambda and G found in AS are not realized in the real world
- Key objection: cutoff identification k <-> r is not unique

**[AS6]** Zholdasbek, A., Chakrabarty, H., Malafarina, D. & Bonanno, A.
"An Emergent Cosmological Model from Running Newton Constant"
Phys. Rev. D 111 (2025) 103519 [arXiv:2405.02636]
- **RECENT**: Lagrangian approach preserving diffeomorphism invariance
- Universe emerges from quasi-de Sitter phase
- Constrains transition scale using CMB data
- Can probe antiscreening of G predicted by AS

### 1.2 Kumar (2025): QFT First-Principles Logarithmic Correction

**[K1]** Kumar, N.
"Marginal IR Running of Gravity as a Natural Explanation for Dark Matter"
Phys. Lett. B 871 (2025) 140008 [arXiv:2509.05246]
- **MOST IMPORTANT PAPER FOR PAPER 3**
- Anomalous dimension: eta(mu) = -mu d(ln G)/dmu
- Running G in IR: G(k) ~ G_N (k_*/k)^eta for k << k_*
- **KEY RESULT**: eta = 1 is the unique marginal case (logarithmic)
- Why eta=1: Fourier transform of k^{-(2+eta)} in d=3 gives ln(r) ONLY for eta=1
- eta < 1: IR-irrelevant (decays too fast)
- eta > 1: violates scale invariance
- **Modified potential:**
```
Phi(r) = -G_N M/r + (2G_N M k_*/pi) ln(r/r_0)
```
- **Force law:**
```
F(r) = -G_N M/r^2 - (2G_N M k_*/pi)(1/r)
```
- **Rotation curve:**
```
v^2(r) = G_N M/r + (2G_N M k_*/pi) = G_N M/r + V_0^2
```
- **Crossover scale:** k_* = pi V_0^2 / (2 G M_bar)
  - Fitted values: k_* ~ (2.6-2.8) x 10^{-2} kpc^{-1} (remarkably universal)
  - r_c = 1/k_* ~ 36-38 kpc
- **Cosmological Friedmann modification:**
```
H^2(z)/H_0^2 = Omega_r0(1+z)^4 + Omega_m0(1+z)^3 + Omega_Lambda0 + Omega_L0(1+z)^2 ln(1+z)
```
- **Logarithmic correction is universal and regulator-independent**
- **LIMITATIONS**: Fails to reproduce CMB acoustic peaks; requires Omega_L0 < 0.05

### 1.3 Casini-Huerta: Relative Entropy and RG Flow

**[CH1]** Casini, H., Teste, E. & Torroba, G.
"Relative Entropy and the RG Flow"
JHEP 03 (2017) 089 [arXiv:1611.00016]
- **KEY THEOREM**: Relative entropy between vacuum states of CFT and CFT+relevant perturbation is monotonically decreasing under RG flow
- By restricting to null Cauchy surface in causal domain of sphere:
  S_rel = Delta S_EE (difference of entanglement entropies)
- Inherits positivity and monotonicity of relative entropy
- **RESULTS**:
  - d=2: Alternative proof of Zamolodchikov c-theorem
  - d>2: Coefficient of area term in S_EE DECREASES along RG flow
- **CRUCIAL IMPLICATION**: If gravitational coupling G is related to area coefficient of entanglement entropy (via Adler-Zee), then RG flow constrains running of G

**[CH2]** Casini, H. & Huerta, M.
"On the RG Running of the Entanglement Entropy of a Circle"
Phys. Rev. D 85 (2012) 125016 [arXiv:1202.5650]
- Entropic c-function from entanglement entropy of intervals
- Strong subadditivity + Lorentz invariance proves c-theorem
- Area term coefficient is concave function of RG scale

**[CH3]** Deddo, E., Liu, J.T., Pando Zayas, L.A. et al.
"Explicit Entropic Proofs of Irreversibility Theorems for Holographic RG Flows"
JHEP 09 (2024) 179 [arXiv:2404.15077]
- **RECENT**: Reproves c-, F-, and a-theorems using only NEC + Ryu-Takayanagi
- Holographic c-function monotonic along flow
- Extends to flows across dimensions (AdS_{D+1} -> AdS_3, AdS_5 -> AdS_4)
- Key: monotonicity of entanglement entropy under coarse-graining = DPI in holographic context

### 1.4 Verlinde Emergent Gravity

**[V1]** Verlinde, E.P.
"Emergent Gravity and the Dark Universe"
SciPost Phys. 2 (2017) 016 [arXiv:1611.02269]
- In de Sitter space: S_entanglement = S_area + S_volume
- Volume-law entanglement overtakes area law at cosmological horizon
- Sub-Hubble: de Sitter states don't thermalize -> entropy displacement by matter
- **Dark matter formula:**
```
M_D^2(r) = (cH_0/(6G)) r^2 d/dr[r M_B(r)]
```
- For point mass: M_D(r) = sqrt(cH_0 M_B r / (6G))
- At low acceleration: g ~ sqrt(a_0 g_Newton) with a_0 ~ cH_0/6
- **LIMITATIONS**: Only valid for static, spherically symmetric; no CMB prediction; no cosmological evolution

### 1.5 SPARC Galaxy Data and Running G Fits

**[SP1]** Lelli, F., McGaugh, S.S., Schombert, J.M. & Pawlowski, M.S.
"One Law to Rule Them All: The Radial Acceleration Relation of Galaxies"
ApJ 836 (2017) 152
- SPARC: 175 disk galaxies with [3.6] photometry + HI/Ha rotation curves
- Radial Acceleration Relation (RAR): g_obs correlates tightly with g_bar over 4 decades
- Critical acceleration scale: a_0 ~ 10^{-10} m/s^2
- Universal across LTGs, ETGs, and classical dSphs

**[SP2]** Bhatia, E., Chakrabarti, S. & Chakraborty, S.
"Phenomenology of Renormalization Group Improved Gravity from the Kinematics of SPARC Galaxies"
[arXiv:2403.00531]
- Applied RGGR (logarithmic running G) to 100 SPARC galaxies
- **KEY RESULT**: RGGR fits observed kinematics consistently
- Parameter nu_bar has near-linear dependence on baryonic mass
- nu_bar decreases from early-type to starburst galaxies
- Satisfies both RAR and Baryonic Tully-Fisher relation
- 1-2 universal parameters vs NFW's 2-3 per galaxy

**[SP3]** Ghari, A. & Haghi, H. (actually: Yoon, Y., Han, S. & Hwang, H.S.)
"Comparison of MOND and Verlinde's Emergent Gravity in Dwarf Spheroidals"
[arXiv:2601.01715] (January 2026)
- 23 dwarf spheroidal galaxies
- **KEY RESULT**: Verlinde's EG favored over MOND at 5.2 sigma
- 21/23 samples: Verlinde follows observed values more closely
- Individual significance: -0.25 sigma to 3.41 sigma

### 1.6 DPI and Gravitational Information

**[DPI1]** Dorau, P. & Much, A.
"From Quantum Relative Entropy to the Semiclassical Einstein Equations"
Phys. Rev. Lett. 136 (2026) 091602 [arXiv:2510.24491]
- **MOST RECENT KEY PAPER**: QRE -> Einstein equations
- Uses modular theory (Araki-Uhlmann entropy) on bifurcate Killing horizons
- QRE between vacuum and coherent excitations = energy flux across horizon
- Under Bekenstein-Hawking assumption: energy flux proportional to area variation
- **Generalizes Jacobson 1995** from thermodynamic to quantum relative entropy
- Result: semiclassical Einstein equations follow from QRE + area-entropy

**[DPI2]** Jacobson, T.
"Thermodynamics of Spacetime: The Einstein Equation of State"
Phys. Rev. Lett. 75 (1995) 1260 [arXiv:gr-qc/9504004]
- Einstein equation from delta Q = T dS on local Rindler horizons
- Foundation for all "gravity from entropy" approaches

**[DPI3]** Jacobson, T.
"Entanglement Equilibrium and the Einstein Equation"
Phys. Rev. Lett. 116 (2016) 201101 [arXiv:1505.04753]
- Einstein equation from maximizing entanglement entropy in geodesic balls
- First-order variation of vacuum state: entanglement stationary iff Einstein holds

**[DPI4]** Faulkner, T., Lewkowycz, A. & Maldacena, J.
"Quantum Corrections to Holographic Entanglement Entropy"
JHEP 11 (2013) 074 [arXiv:1307.2892]
- FLM formula: S = A/(4G_N) + S_bulk
- Quantum corrections to Ryu-Takayanagi at order G^0

### 1.7 Bianconi: QRE as Gravitational Action

**[B1]** Bianconi, G.
"Gravity from Entropy"
Phys. Rev. D 111 (2025) 066001 [arXiv:2408.14391]
- **KEY FRAMEWORK**: Action = QRE between spacetime metric and matter-induced metric
- Metric of spacetime -> quantum operator (effective density matrix)
- Matter: Dirac-Kahler formalism (0-form + 1-form + 2-form)
- Low coupling: reduces to Einstein equations with Lambda=0
- **Predicts**: small positive cosmological constant
- Equations remain second-order

**[B2]** Bianconi, G.
"The Quantum Relative Entropy of the Schwarzschild Black Hole and the Area Law"
Entropy 27 (2025) 266 [arXiv:2501.09491]
- QRE of Schwarzschild metric obeys area law for large r_s
- Schwarzschild is approximate solution (valid in low coupling, small curvature)
- Geometric QRE (GQRE) framework

**[B3]** Bianconi, G.
"The Thermodynamics of the Gravity from Entropy Theory"
[arXiv:2510.22545] (October 2025)
- Thermodynamics of GfE: GQRE per unit volume is NOT increasing (consistent with relative entropy)
- Total entropy of GfE universes is non-decreasing
- Friedmann universes as solutions in low-energy limit
- **CRUCIAL**: provides natural thermodynamic interpretation for cosmological solutions

### 1.8 Weak-Field Quantum Gravity Corrections

**[QC1]** Donoghue, J.F.
"General Relativity as an Effective Field Theory: The Leading Quantum Corrections"
Phys. Rev. D 50 (1994) 3874 [arXiv:gr-qc/9405057]
- **FOUNDATIONAL**: Pioneered EFT approach to quantum gravity
- Leading corrections from massless particle propagation (including gravitons)
- Non-analytic/nonlocal contributions are parameter-free
- **Quantum-corrected potential:**
```
V(r) = -G m1 m2 / r [1 + 3G(m1+m2)/(c^2 r) + 41 G hbar / (10 pi c^3 r^2)]
```
- The r^{-2} term is the quantum correction (universal, UV-completion independent)

**[QC2]** Bjerrum-Bohr, N.E.J., Donoghue, J.F. & Holstein, B.R.
"Quantum Gravitational Corrections to the Nonrelativistic Scattering Potential of Two Masses"
Phys. Rev. D 67 (2003) 084033; erratum D 71 (2005) 069903 [arXiv:hep-th/0211072]
- Full one-loop nonanalytic scattering amplitude
- Post-Newtonian + quantum corrections to order G^2
- **Logarithmic quantum correction**: depends only on low-energy spectrum
- Related to Landau singularities of one-loop amplitude

**[QC3]** de Paula Netto, T., Modesto, L. & Shapiro, I.L.
"Universal Leading Quantum Correction to the Newton Potential"
Eur. Phys. J. C 82 (2022) 160 [arXiv:2110.14263]
- Recalculates corrections using functional methods
- Proves gauge- and parametrization-independence
- Universality: logarithmic correction depends only on low-energy spectrum
- **Key**: the log correction is universal across all UV completions of quantum gravity

### 1.9 Grumiller's Model

**[G1]** Grumiller, D.
"Model for Gravity at Large Distances"
Phys. Rev. Lett. 105 (2010) 211303 [arXiv:1011.3625]
- Effective model for gravity at large scales
- Leading order: cosmological constant + **Rindler acceleration** + Newtonian + subleading
- V(r) = -GM/r + a_Rindler * r (Rindler term)
- Rindler acceleration a ~ 0.30 x 10^{-10} m/s^2 (close to a_0!)
- Explains ~10% of Pioneer anomaly + ameliorates rotation curves
- Does not spoil solar system tests

### 1.10 Galactic Spacetime with Varying G

**[VG1]** Chakraborty, R. & Sengupta, S.
"Galactic Spacetime Solutions with a Varying Newton's Coupling"
[arXiv:2505.04060] (May 2025)
- **RECENT**: New family of galactic metrics with flat rotation curves
- Vacuum solutions in gravity theory with spatially varying G
- Effective "mass" receives negative non-baryonic contribution (purely geometric)
- Light deflection diminished (not enhanced) vs Einstein
- Observationally distinguishable from dark matter alternatives

### 1.11 CMB and Running G

**[CMB1]** Alvarez, P.D. et al.
"From Cavendish to Planck: Constraining Newton's Gravitational Constant with CMB"
JCAP (various years)
- Planck + DESI BAO: G constrained to 1.8% precision
- Effective gravitational constant on cosmological scales: G_eff = G_N (1.0 +/- 0.04)

**[CMB2]** Modified gravity and CMB acoustic peaks (general constraint):
- Three independent effects: (a) coupling to CDM changes peak amplitudes; (b) coupling to baryons shifts peaks; (c) both couplings -> Newton potential in Sachs-Wolfe effect
- **THE PROBLEM**: Running G alone cannot produce the correct ratio of odd/even CMB peaks (which requires non-baryonic matter with different sound speed)

### 1.12 Holographic RG and Entanglement

**[H1]** Hartman, T.
"Entanglement Entropy and the Renormalization Group" (lecture notes)
- Adler-Zee formula: 1/G_ren = 1/G_bare + sum of QFT contributions
- Holographic: renormalization of central charge = renormalization of G
- Area coefficient of S_EE encodes running of G

---

## 2. ANALYSIS: KEY QUESTIONS

### 2.1 Can Sigma = D(rho_spacetime || rho_matter) give both r_s/r and r_s/r + alpha*ln(r)?

**YES, through the gravitational channel's scale-dependent running.**

The argument proceeds in three steps:

**Step 1: Sigma = -ln(g_00) at all field strengths (Paper 2 result)**
- Dorau-Much [DPI1]: QRE -> Einstein equations via modular flow
- Our Paper 2: Sigma_grav = r_s/r for Schwarzschild (three independent derivations)

**Step 2: Include quantum loop corrections to the gravitational channel**
- At one loop (Donoghue [QC1], Bjerrum-Bohr [QC2]):
```
Phi(r) = -G_N M/r [1 + corrections at O(G/r) + O(G hbar/r^2)]
```
These are negligible at galactic scales (too small by ~60 orders of magnitude).

**Step 3: Include IR running of G (Kumar [K1])**
- The gravitational channel N_{grav} itself has scale-dependent coupling
- At momentum scale k: G(k) = G_N (k_*/k)^eta for k << k_*
- For eta = 1 (marginal, unique): G(k) ~ G_N (1 + k_*/k)
- This gives:
```
Sigma(r) = -ln(g_00) = 2|Phi|/c^2 = r_s/r + (4G_N M k_*)/(pi c^2) * ln(r/r_0)
```
- **The logarithmic correction emerges from the SAME QRE, not an additional term**
- The gravitational channel N changes at different scales -> Sigma changes

**ASSESSMENT**: This is the most natural route. The key insight is that Sigma = D(rho||sigma) - D(N(rho)||N(sigma)) depends on the channel N, and N itself runs with scale. In Paper 2, N is the strong-field gravitational channel (Schwarzschild). In Paper 3, N includes IR loop corrections that are negligible near black holes but dominate at galactic scales.

**Remaining gap**: The running G(k) is not derived from QRE first principles. Kumar uses EFT arguments. To make this fully self-consistent, we need:
- QRE monotonicity (Casini-Huerta) applied to gravitational coupling
- OR: Bianconi's framework extended to include running

### 2.2 Is there a QRE-based derivation of running G(k)?

**NOT YET, but the ingredients exist.**

**Existing pieces:**
1. **Casini-Huerta [CH1]**: QRE monotonically decreases under RG flow. The area coefficient of S_EE decreases along flow.
2. **Adler-Zee formula [H1]**: 1/G_ren = 1/G_bare + QFT contributions. The running of G is related to the running of the area coefficient.
3. **Dorau-Much [DPI1]**: QRE -> Einstein equations. This establishes QRE as the correct information-theoretic quantity for gravity.
4. **Bianconi [B1]**: Action = QRE between metrics. This gives modified Einstein equations with running.

**The missing link**: No one has combined Casini-Huerta + Dorau-Much to derive:
```
d G(k) / d ln(k) = -eta * G(k)   from QRE monotonicity
```

**Proposed route for Paper 3:**
- Start from Sigma = D(rho_spacetime || rho_matter) [Bianconi-type]
- Apply RG flow to the gravitational channel
- Use Casini-Huerta monotonicity to show that Sigma_grav decreases under flow
- The rate of decrease is eta, constrained by DPI

### 2.3 Can the DPI constrain the anomalous dimension eta?

**PLAUSIBLE, but not proven.**

The DPI states:
```
D(N(rho) || N(sigma)) <= D(rho || sigma)
```

If we identify:
- N = coarse-graining from UV to IR (the RG flow itself)
- D(rho || sigma) = Sigma_grav at UV scale
- D(N(rho) || N(sigma)) = Sigma_grav at IR scale

Then DPI says Sigma decreases under coarse-graining, which is Casini-Huerta's result.

**To constrain eta**: We need the RATE of decrease. The Fawzi-Renner bound gives:
```
D(rho || sigma) - D(N(rho) || N(sigma)) >= -2 ln F(rho, P circ N(rho))
```
where P is the Petz recovery map.

**Interpretation for gravity:**
- The "recovery deficit" = how much Sigma is lost under RG flow from UV to IR
- If the gravitational channel at scale k has Petz recovery fidelity F_k, then:
```
Delta Sigma = Sigma(k_UV) - Sigma(k_IR) >= -2 ln F_Petz
```

**KEY CONJECTURE FOR PAPER 3:**
```
eta <= 1 - exp(-Sigma/2)
```
This would make eta = 1 the maximal allowed value (saturating the bound when Sigma -> infinity), and would be a DPI-derived constraint on the anomalous dimension.

**Current status**: Speculative but physically motivated. No proof exists.

---

## 3. ASSESSMENT: MOST PROMISING ROUTE

### Route A: Pure QRE (Strongest but Hardest)
1. Bianconi's GQRE action -> modified Einstein equations
2. Include quantum corrections -> running G emerges
3. Casini-Huerta constrains the running via DPI
4. At galactic scales: logarithmic correction natural

**Pros**: Fully self-consistent, everything from one QRE
**Cons**: Bianconi's framework not yet developed enough; no explicit G(k)
**Likelihood**: Low for Paper 3 timeline, but ultimate goal for Paper 4

### Route B: EFT Running G Interpreted via QRE (Most Promising)
1. Kumar's result: eta=1 gives ln(r) correction (proven)
2. Reinterpret: running G = scale-dependent gravitational channel
3. Sigma(r) = -ln(g_00) with quantum-corrected g_00
4. DPI provides conceptual constraint on eta (new contribution)

**Pros**: Builds on proven results, adds QRE interpretation as new layer
**Cons**: Not fully derived from QRE; uses EFT as input
**Likelihood**: HIGH for Paper 3
**This is the recommended approach.**

### Route C: Verlinde's Volume-Law Entanglement (Complementary)
1. De Sitter background -> volume-law entanglement
2. Entropy displacement by matter -> apparent dark matter
3. Connection to Sigma: Sigma_dS includes volume-law contribution

**Pros**: Conceptually beautiful, works for dwarf spheroidals (5.2 sigma)
**Cons**: No CMB prediction, theoretical inconsistencies (Hossenfelder, Dai-Stojkovic)
**Likelihood**: Complementary to Route B, not standalone

### RECOMMENDED STRATEGY:
Paper 3 should use **Route B as primary + Route C as complementary**.
Paper 4 should aim for **Route A** (full QRE derivation).

---

## 4. ASSESSMENT: CAN THE CMB PROBLEM BE SOLVED?

### The Problem
CMB acoustic peaks require:
1. Baryonic matter (couples to photons, provides pressure)
2. Non-baryonic matter (provides gravity but no pressure -> different sound speed)
3. The ratio of odd/even peaks depends on the baryon-to-total-matter ratio

Running G alone modifies the gravitational strength but does NOT provide a second species with different pressure properties. This is the fundamental obstacle.

### Possible Solutions (from Literature)

**Solution 1: Scale-dependent G at CMB epoch (Partial)**
- Zholdasbek et al. [AS6]: Constrains transition scale using CMB data
- Running G affects early universe dynamics
- BUT: cannot replace the role of non-baryonic matter in determining peak ratios

**Solution 2: Kumar's Cosmological Modification (Insufficient)**
- Modified Friedmann: includes Omega_L0(1+z)^2 ln(1+z) term
- This acts like a new energy component
- BUT: Kumar finds Omega_L0 < 0.05 from H(z) data -> too small to affect CMB peaks significantly

**Solution 3: Volume-Law Entanglement (Untested)**
- Verlinde's thermal volume-law entropy could act as effective dark matter at CMB epoch
- Temperature-dependent: at z~1100, the de Sitter temperature is much lower
- NOT developed for CMB; Verlinde himself acknowledges this gap

**Solution 4: Separate the Roles (Most Honest)**
- Accept that running G explains galactic-scale anomalies
- CMB peaks require additional physics (possibly related to early-universe entanglement structure)
- Paper 3 should be HONEST about this limitation
- Connection to Paper 4: Sigma = D(rho_spacetime || rho_matter) in early universe may have different structure

**Solution 5: Non-Markovian Gravitational Channel (Speculative)**
- If the gravitational channel has memory effects (non-Markovian)
- At CMB epoch: shorter timescales -> stronger memory effects
- Could provide effective "mass" through information-theoretic mechanism
- No existing calculation supports this

### ASSESSMENT
**The CMB problem CANNOT be fully solved within Paper 3's framework.** The most honest approach:
1. Paper 3 focuses on galactic-scale success (rotation curves, dwarf spheroidals, RAR, BTFR)
2. Acknowledge CMB as a limitation
3. Point out that the SAME limitation applies to MOND, Verlinde, and all modified gravity theories
4. Suggest that Paper 4's full QRE framework (Sigma in early universe) may resolve this
5. Note: even LCDM required decades to get CMB right; modified gravity is younger

---

## 5. THE BRIDGE: How Paper 2 Connects to Paper 3

### Mathematical Bridge
```
Paper 2 (strong field):
  Sigma_grav = r_s/r = -ln(g_00) for Schwarzschild
  Valid for: r >> r_s (weak field limit of strong-field solution)

Paper 3 (weak field, galactic):
  Sigma_grav = r_s/r + alpha * ln(r/r_c)
  where alpha = (4 G_N M k_*)/(pi c^2)
  Valid for: r ~ 1-100 kpc (galactic scales)

Connection: SAME formula Sigma = -ln(g_00), but g_00 now includes
quantum loop corrections to the gravitational channel:
  g_00 = 1 - r_s/r  ->  g_00 = 1 - r_s/r - (4G_N M k_*)/(pi c^2) ln(r/r_0)
```

### Physical Bridge
- **Paper 2**: Sigma measures information loss to gravitational degrees of freedom (spacetime curvature)
- **Paper 3**: At galactic scales, the gravitational channel N has IR running -> more information loss at larger distances -> ADDITIONAL Sigma from running G -> this extra Sigma appears as "dark matter"
- **Unifying idea**: "Dark matter" = extra information loss from IR running of the gravitational channel

### Information-Theoretic Bridge
- **DPI**: D(N_IR(rho) || N_IR(sigma)) <= D(N_UV(rho) || N_UV(sigma))
- The gravitational channel at IR (galactic) scales loses MORE information than at UV (solar system) scales
- This extra loss manifests as additional gravitational effect (apparent dark matter)
- tau_IR > tau_UV: more time arrow at larger scales (more irreversibility)

---

## 6. KEY FORMULAS FOR PAPER 3

### 6.1 Running G in Sigma Language
```
Sigma(r) = 2 G(r) M / (c^2 r)
         = (r_s / r) * [1 + f(r)]
```
where f(r) encodes the running:
- f(r) = 0: Newtonian (no running)
- f(r) = (2 k_* r / pi) * ln(r/r_0): marginal IR running (eta=1)

### 6.2 DPI-Motivated Bound on eta
```
Proposed: eta <= 1 - exp(-Sigma_grav/2)
```
- Sigma -> 0 (flat space): eta -> 0 (no running in flat space) [CORRECT]
- Sigma -> infinity (strong field): eta -> 1 (maximal running) [CONSISTENT with Kumar]
- eta = 1 saturates when Sigma is large enough [NATURAL]

### 6.3 Verlinde's Formula in Sigma Language
```
Sigma_Verlinde(r) = r_s/r + Sigma_dS(r)
where Sigma_dS(r) = (r^2 / L_dS^2) * M_B / M_dS
```
- L_dS = Hubble radius = c/H_0
- M_dS = c^3/(2GH_0)

### 6.4 Rotation Curve from Sigma
```
v^2(r) = (c^2/2) * r * dSigma/dr
        = G_N M/r + (2G_N M k_*)/(pi)      [Kumar]
        = G_N M/r + sqrt(a_0 G_N M)         [Verlinde, at low g]
```

---

## 7. OPEN QUESTIONS AND NEXT STEPS

### For Paper 3:
1. **Derive eta=1 from QRE**: Can Casini-Huerta + Dorau-Much give eta=1 directly?
2. **Test the DPI bound on eta**: Does eta <= 1-exp(-Sigma/2) hold for known solutions?
3. **SPARC analysis**: Apply Sigma(r) = r_s/r + alpha*ln(r/r_c) to full SPARC catalog
4. **CMB section**: Honest discussion of limitations + future directions
5. **Bullet Cluster**: Can the framework handle cluster collisions?

### For Paper 4:
1. **Full Bianconi extension**: Include running G in GQRE framework
2. **Early universe Sigma**: What does Sigma = D(rho_spacetime || rho_matter) look like at z~1100?
3. **Unify Kumar + Verlinde**: Both as limits of the same QRE
4. **CMB from QRE**: Can the full framework (not just running G) reproduce acoustic peaks?

---

## 8. PAPER-BY-PAPER CONNECTION MAP

```
Paper 1: tau = 1-F, F >= exp(-Sigma/2)
  |
  | (Sigma_grav = r_s/r via 3 derivations)
  v
Paper 2: Sigma_grav = -ln(g_00) = r_s/r [Schwarzschild]
  |
  | (Include IR running of gravitational channel)
  v
Paper 3: Sigma_grav = r_s/r + alpha*ln(r/r_c) [Galactic]
  |        - Running G from EFT (Kumar): proven, tested
  |        - QRE interpretation: Casini-Huerta constrains running
  |        - DPI bound on eta: new contribution
  |        - Verlinde as complementary picture
  |
  | (Full QRE derivation, early universe)
  v
Paper 4: Sigma = D(rho_spacetime || rho_matter) [Universal]
           - One equation, different boundary conditions -> different physics
           - Strong field (Paper 2), weak field (Paper 3), cosmological (CMB)
```

---

## 9. SUMMARY TABLE OF ALL REFERENCES

| Ref ID | Authors | Year | arXiv | Key Result | Relevance |
|--------|---------|------|-------|------------|-----------|
| AS1 | Reuter & Saueressig | 2012 | - | AS review | Background |
| AS2 | Reuter & Saueressig | 2019 | - | AS monograph | Background |
| AS3 | Reuter & Weyer | 2004 | hep-th/0410117 | Running G -> rotation curves | Core for Paper 3 |
| AS4 | Koch & Ramirez | 2011 | 1010.2799 | Optimal scale RG | Background |
| AS5 | Donoghue | 2020 | 1911.02967 | Critique of AS | Must address |
| AS6 | Zholdasbek et al. | 2025 | 2405.02636 | Running G + CMB | CMB constraints |
| K1 | Kumar | 2025 | 2509.05246 | eta=1, ln correction | **PRIMARY** |
| CH1 | Casini, Teste, Torroba | 2017 | 1611.00016 | QRE monotone under RG | **KEY THEOREM** |
| CH2 | Casini & Huerta | 2012 | 1202.5650 | Entropic c-function | Supporting |
| CH3 | Deddo et al. | 2024 | 2404.15077 | Entropic irreversibility | Supporting |
| V1 | Verlinde | 2016 | 1611.02269 | Emergent gravity + DM | Complementary |
| SP1 | Lelli et al. | 2017 | - | SPARC + RAR | Data source |
| SP2 | Bhatia et al. | 2024 | 2403.00531 | RGGR fits SPARC | **VALIDATION** |
| SP3 | Yoon et al. | 2026 | 2601.01715 | Verlinde > MOND 5.2sigma | Supporting |
| DPI1 | Dorau & Much | 2026 | 2510.24491 | QRE -> Einstein | **KEY LINK** |
| DPI2 | Jacobson | 1995 | gr-qc/9504004 | Thermodynamic Einstein | Foundation |
| DPI3 | Jacobson | 2016 | 1505.04753 | Entanglement equilibrium | Foundation |
| DPI4 | Faulkner et al. | 2013 | 1307.2892 | FLM formula | Background |
| B1 | Bianconi | 2025 | 2408.14391 | QRE as gravity action | **KEY FRAMEWORK** |
| B2 | Bianconi | 2025 | 2501.09491 | QRE of Schwarzschild | Supporting |
| B3 | Bianconi | 2025 | 2510.22545 | Thermodynamics of GfE | Supporting |
| QC1 | Donoghue | 1994 | gr-qc/9405057 | GR as EFT | Foundation |
| QC2 | Bjerrum-Bohr et al. | 2003 | hep-th/0211072 | Quantum corrections | Foundation |
| QC3 | de Paula Netto et al. | 2022 | 2110.14263 | Universal log correction | Supporting |
| G1 | Grumiller | 2010 | 1011.3625 | Rindler acceleration | Alternative |
| VG1 | Chakraborty & Sengupta | 2025 | 2505.04060 | Varying G metrics | Recent |

---

## 10. CRITICAL ASSESSMENT: WHAT IS NEW IN OUR FRAMEWORK

### Already Known:
- Running G from asymptotic safety (Reuter 2004)
- Logarithmic correction from eta=1 (Kumar 2025)
- SPARC fits with running G (Bhatia et al. 2024)
- QRE monotonicity under RG flow (Casini-Huerta 2017)
- QRE -> Einstein equations (Jacobson 1995, Dorau-Much 2026)

### New in Our Paper 3:
1. **Sigma(r) = r_s/r + alpha*ln(r/r_c)** as a unified information-theoretic quantity bridging strong and weak field
2. **DPI-motivated bound on eta**: eta <= 1 - exp(-Sigma/2) [NEW CONJECTURE]
3. **Running G = scale-dependent gravitational channel**: reinterpretation of EFT running in quantum information language
4. **"Dark matter" = extra information loss**: the extra Sigma from IR running manifests as apparent dark matter
5. **Connection to Paper 1's tau**: at galactic scales, tau measures the degree to which gravitational information is lost, and the extra tau from running G is what we call "dark matter"

### What We Do NOT Claim:
- We do NOT derive running G from QRE (that's Paper 4's goal)
- We do NOT solve the CMB problem
- We do NOT explain the Bullet Cluster
- We do NOT claim this replaces dark matter at ALL scales

---

*Last updated: 2026-03-10*
*Research conducted for Paper 3 of the four-paper series: Sigma = D(rho_spacetime || rho_matter)*
