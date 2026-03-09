# Paper 3 Research: Flat Rotation Curves from Quantum Relative Entropy

**Research Date**: 2026-03-06
**Goal**: Show that the same Sigma = D(rho||sigma) - D(N(rho)||N(sigma)) (quantum relative entropy drop under gravitational channel), when evaluated at galactic scales, naturally produces flat rotation curves WITHOUT dark matter.

**Core philosophy**: NOT adding terms to Sigma (patchwork), but showing the SAME Sigma, computed with the full quantum-corrected gravitational channel, naturally gives different behavior at different scales. Like F=ma giving different solutions for different systems.

---

## 1. Kumar (2025, arXiv:2509.05246) -- DETAILED ANALYSIS

**Title**: "Marginal IR running of Gravity as a Natural Explanation for Dark Matter"
**Author**: Naman Kumar
**Published**: Physics Letters B 871 (2025) 140008

### 1.1 Core Derivation

The running of Newton's coupling in the IR:

```
G(k) ~ G_N (k_*/k)^eta   for k << k_*
```

The anomalous dimension eta governs the flow:

```
d ln G / d ln k = -eta
```

### 1.2 Why eta = 1 is Special (Marginal)

Fourier transform of k^(-(2+eta)) in 3 spatial dimensions:
- eta != 1: power-law behavior ~ r^(eta-1)
- **eta = 1: logarithmic behavior ~ ln(r)** (UNIQUE marginal case)

This is the ONLY scale-invariant deformation consistent with rotational symmetry and locality:
- eta < 1: decays too fast (IR-irrelevant)
- eta > 1: grows too fast (spoils scale invariance)
- **eta = 1: marginal, logarithmic, universal**

### 1.3 Derivation of Logarithmic Potential

Starting from Newtonian potential in momentum space:

```
Phi(r) = -M integral d^3k/(2pi)^3 e^{ik.r} 4piG(k)/k^2
```

With G(k) = G_N[1 + (k_*/k)] for eta=1:

```
Phi(r) = -G_N M/r + (2 G_N M k_*/pi) ln(r/r_0) + const
```

Force law:

```
F(r) = -dPhi/dr = -G_N M/r^2 - (2 G_N M k_*/pi) * (1/r)
```

At r >> r_c = 1/k_*, the **1/r** term dominates over **1/r^2**.

### 1.4 Crossover Scale k_* and r_c

```
k_* = pi V_0^2 / (2 G M_bar)
r_c = 1/k_*
```

For three sample galaxies:
- k_* ~ (2.6-2.8) x 10^{-2} kpc^{-1}
- **r_c ~ 36-38 kpc** (remarkably consistent!)

| Galaxy   | M_bar [M_sun]  | V_0 [km/s]    | k_* [kpc^{-1}]          | r_c [kpc]    |
|----------|----------------|----------------|--------------------------|--------------|
| S:700624 | 1.01 x 10^{12} | 271.4 +/- 2.6 | (2.66 +/- 0.05) x 10^{-2} | 37.5 +/- 0.7 |
| S:700916 | 5.85 x 10^{11} | 210.4 +/- 2.4 | (2.76 +/- 0.06) x 10^{-2} | 36.2 +/- 0.8 |
| S:705253 | 3.55 x 10^{11} | 158.9 +/- 1.9 | (2.60 +/- 0.06) x 10^{-2} | 38.5 +/- 0.9 |

### 1.5 Rotation Curve Formula

```
v^2(r) = G_N M_bar / r + (2 G_N M_bar k_*) / pi
```

At r >> r_c: v^2 ~ constant --> flat rotation curves!

### 1.6 Free Parameters

**Galactic scale**: 2 parameters
- k_* (crossover scale)
- G_N (Newton's constant)

The coefficient **2/pi** in the logarithmic term is **universal and regulator-independent**.

### 1.7 CRITICAL LIMITATIONS

1. **CMB**: FAILS to reproduce acoustic peak structure without dark matter. IR effect too weak at z ~ 1100.
2. **Cosmological**: Late-time H(z) data require Omega_L0 < 0.05.
3. **Cluster scales**: No systematic analysis of Bullet Cluster or cluster dynamics.
4. **Structure formation**: Cannot explain high-z structure formation without additional matter.

### 1.8 Connection to Asymptotic Safety

**Explicitly NOT discussed.** Kumar treats eta=1 as emerging from IR RG flow in EFT framework, NOT from UV fixed points. The paper references Reuter & Weyer but emphasizes EFT marginality arguments.

### 1.9 Connection to QRE/DPI

**None mentioned.** No information theory in this paper.

### 1.10 Assessment for Paper 3

**Strengths**:
- Clean derivation from first principles (RG flow)
- Universal coefficient (2/pi)
- Remarkably consistent k_* across galaxies
- Only 2 parameters at galactic scale

**Weaknesses**:
- No connection to quantum information theory (must be built)
- CMB problem is severe
- Not clear this is the SAME as "running from QRE channel"

**KEY INSIGHT for Paper 3**: If we can show that the gravitational channel N_grav has quantum corrections that produce G(k) with eta=1, then Kumar's logarithmic potential IS the natural consequence of computing Sigma with the full quantum channel. The 1/r force is not "added" -- it emerges from the RG flow of the channel itself.

---

## 2. Gubitosi, Piattella & Casarini (2024, arXiv:2403.00531) -- SPARC FITTING

**Title**: "Phenomenology of renormalization group improved gravity from the kinematics of SPARC galaxies"

### 2.1 RGGR Model Equations

Beta function for running G:

```
beta = mu dG^{-1}/dmu = 2nu G_0^{-1}
```

Running coupling:

```
G(mu) = G_0 / [1 + nu ln(mu^2/mu_0^2)]
```

Effective circular velocity:

```
v^2_RGGR(r) = v^2_N(r) [1 - c^2 nu_bar / phi_N(r)]
```

where nu_bar = nu * alpha (phenomenological parameter).

### 2.2 Galaxy Sample

- 100 galaxies from SPARC (out of 175), selected after quality cuts
- 4 morphological types: Early (9), Spiral (39), Late dwarf (34), Starburst (18)

### 2.3 Free Parameters Per Galaxy

Three parameters via MCMC:
- gamma_d (disk mass-to-light ratio): [0.3, 0.8]
- gamma_b (bulge mass-to-light ratio): [0.3, 0.8]
- nu_bar (RGGR parameter): 10^{-9} to 10^{-6}

### 2.4 Comparison: RGGR vs NFW

Using DELTA_BIC = BIC_NFW - BIC_RGGR:
- 47 galaxies prefer RGGR (DELTA_BIC > 2)
- 40 galaxies prefer NFW
- 13 inconclusive

**RGGR favored overall**, but not overwhelmingly.

No MOND comparison in this paper.

### 2.5 Key Finding: nu_bar Correlates with Mass

**nu_bar has near-linear dependence on galactic baryonic mass**:
- Early type: nu_bar ~ 10^{-6}
- Spiral: nu_bar ~ 10^{-7}
- Late/Starburst: nu_bar ~ 10^{-8}
- Solar system: nu_bar ~ 10^{-17}

### 2.6 Assessment for Paper 3

**Strengths**:
- 100 galaxies fitted systematically with SPARC data
- RGGR competitive with NFW (BIC analysis)
- Clear mass-dependent scaling of nu_bar
- Logarithmic running of G is the key ingredient

**Weaknesses**:
- nu_bar is NOT universal (varies per galaxy, correlates with mass)
- 3 free parameters per galaxy (same as NFW)
- Not a dramatic improvement over dark matter

**KEY INSIGHT**: The logarithmic running G(mu) = G_0 / [1 + nu ln(mu^2/mu_0^2)] is EXACTLY what appears when computing the gravitational channel's quantum corrections. If we can connect nu to the QRE drop Sigma, the mass-dependence of nu_bar becomes a prediction, not a free parameter.

---

## 3. Asymptotic Safety and Running G -- Status Report

### 3.1 Reuter's Program (Overview)

- **Founded**: 1998, Martin Reuter's seminal paper using FRG for gravity
- **Key claim**: Non-Gaussian UV fixed point (Reuter fixed point) exists, making gravity UV-complete
- **G(k) behavior**: G(k) ~ k^{-2} at high k (UV), grows at low k (IR)
- **Book**: Reuter & Saueressig, "Quantum Gravity and the Functional Renormalization Group" (Cambridge, 2019)

### 3.2 Application to Dark Matter

Reuter & Weyer (2004, hep-th/0410117):
- First paper connecting running G to rotation curves
- Power-law running G(r) with small exponent ~ 10^{-6}
- Can account for non-Keplerian rotation curves without dark matter

Rodrigues, Letelier, Fabris, Shapiro (2009-2024):
- Systematic development of RGGR model
- G varies ~ 10^{-7} of its value across a galaxy
- Series of papers culminating in Gubitosi et al. (2024)

### 3.3 Is eta = 1 Established?

**NOT firmly established from first principles.**

- Kumar (2025) derives eta=1 from EFT marginality/dimensional arguments -- this is a self-consistency argument, not a lattice/FRG derivation
- UV fixed point: eta_N = d-2 = 2 (well-established in AS)
- **IR running is much less understood** than UV running
- The flow FROM the UV fixed point TO the IR is the hard part

### 3.4 Criticisms of Asymptotic Safety

**Donoghue's critique** (Frontiers in Physics, 2020):
1. "Running of Lambda and G found in AS are not realized in physical processes"
2. Kinematic incompatibility: s and t channels carry opposite signs; g(s) and g(t) are wildly different
3. Cutoff identification G(k) is regularization-dependent, not physical
4. Power-divergent corrections (unlike logarithmic QCD running) are scheme-dependent

**Bonanno, Eichhorn et al.** (Frontiers in Physics, 2020):
- Acknowledge difficulty of computing observables
- Uncontrollable approximations in FRG
- Background field method concerns
- Difficult to compare with experiment/observation

**BRST symmetry** (Symmetry 2026):
- Fundamental limitations from BRST violation above gravitational cutoff
- Propagator depends on gauge choice and truncation

### 3.5 Current Status Assessment (Honest)

| Aspect | Status |
|--------|--------|
| UV fixed point exists | Strong evidence (FRG, lattice hints) |
| UV anomalous dimension | Well-established: eta_N = 2 in d=4 |
| IR running of G | Poorly understood, controversial |
| eta = 1 at galactic scales | **Argued by dimensionality, not derived** |
| Rotation curves from running G | Phenomenologically viable (Gubitosi+) |
| Observers from AS | Limited, criticized (Donoghue) |
| Connection to QI/QRE | **NON-EXISTENT in current literature** |

### 3.6 Assessment for Paper 3

**The gap**: Asymptotic safety people don't talk about quantum information. QI people don't talk about running couplings at galactic scales. **Paper 3 fills this gap.**

---

## 4. Can Running G Be Derived from QRE/DPI?

### 4.1 The Bridge: Casini-Huerta Theorem

**Key paper**: Casini & Huerta, "Relative entropy and the RG flow" (JHEP 2017, arXiv:1611.00016)

**Result**: S_rel(rho_UV || rho_IR) when restricted to null Cauchy surface in causal domain of a sphere equals the difference of entanglement entropies, which has positivity and monotonicity properties of relative entropy.

**Implications**:
1. Alternative proof of c-theorem in d=2 (Zamolodchikov)
2. For d>2: coefficient of area term in entanglement entropy DECREASES along RG flow
3. Monotonicity of relative entropy = irreversibility of RG flow = DPI!

**THIS IS THE MATHEMATICAL BRIDGE**: RG flow IS data processing. The monotonicity of QRE under quantum channels (DPI) IS the c-theorem for gravity.

### 4.2 Dorau-Much (2025, arXiv:2510.24491, PRL)

**Title**: "From Quantum Relative Entropy to the Semiclassical Einstein Equations"

**Derivation**:
1. Start with local Rindler horizon (equivalence principle)
2. Compute S_rel(omega_0 || omega_phi) between vacuum and coherent excitation using Araki-Uhlmann modular theory
3. Result: S_rel = -2pi integral_{H_B^R} U (partial_U phi)^2 dU dvol_S
4. Assume Bekenstein-Hawking: delta_A = 4 S_rel
5. Apply Raychaudhuri equation
6. **Einstein equations follow**: R_ab - (R/2)g_ab + Lambda g_ab = alpha <:T_ab:>

**Limitations**:
- Semiclassical only (no beyond-semiclassical corrections)
- Coherent states only
- NO mention of running G, loop corrections, or RG flow
- NO mention of DPI

### 4.3 The Logical Chain for Paper 3

**The argument we need to construct**:

1. **Dorau-Much**: QRE --> semiclassical Einstein equations (established, PRL 2025)
2. **Casini-Huerta**: QRE monotonicity under RG flow = c-theorem (established, JHEP 2017)
3. **NEW (Paper 3)**: Quantum corrections to the gravitational channel N_grav modify D(rho||sigma) - D(N(rho)||N(sigma)) at each scale k
4. **KEY STEP**: The RG-improved gravitational channel N_grav(k) has:
   - At small k (large r): enhanced information loss --> G grows --> logarithmic corrections
   - At large k (small r): standard Newton --> Sigma ~ r_s/r as in Paper 2
5. **RESULT**: The SAME Sigma formula, evaluated with N_grav(k), naturally gives flat rotation curves at galactic scales

### 4.4 Mathematical Sketch

Define the gravitational channel at scale k:

```
N_grav(k): rho --> e^{-iH(k)t} rho e^{iH(k)t}  (+ tracing over environment)
```

where H(k) uses G(k) instead of G_N.

Then:

```
Sigma(k) = D(rho || sigma) - D(N_grav(k)(rho) || N_grav(k)(sigma))
```

By DPI: Sigma(k) >= 0 always.

As k decreases (larger distances):
- G(k) increases (IR running with eta=1)
- N_grav(k) becomes a "noisier" channel
- Sigma(k) increases (MORE information loss)
- This enhanced Sigma produces extra gravitational potential

**Connection to Kumar**: If G(k) = G_N(k_*/k) for k << k_*, then:

```
Sigma_grav(r) = Sigma_Newton(r) + Sigma_quantum(r)
             = (r_s/r) + (2 k_* r_s / pi) ln(r/r_0)
```

The logarithmic correction IS the additional information loss due to quantum corrections.

### 4.5 Honest Assessment

**What is established**:
- QRE --> Einstein equations (Dorau-Much, PRL)
- QRE monotonicity = RG irreversibility (Casini-Huerta)
- Running G can fit rotation curves (Kumar, Gubitosi+)

**What needs to be constructed (our contribution)**:
- The explicit connection: QRE corrections AT EACH SCALE --> running G
- Showing eta=1 follows from QRE channel properties
- Showing the SAME Sigma formula gives Newton at small r and flat curves at large r

**Difficulty level**: HIGH. This is not a simple calculation. It requires combining:
- Algebraic QFT (modular theory)
- Functional RG (asymptotic safety)
- Quantum information (DPI/recovery maps)

---

## 5. Alternative Approaches: Dark Matter from Information Theory

### 5.1 Verlinde's Emergent Gravity (2016, arXiv:1611.02269)

**Mechanism**: Positive dark energy creates volume-law entropy competing with area-law entanglement entropy --> extra "dark gravity force" at sub-Hubble scales.

**Results**:
- Derives Tully-Fisher relation
- Reproduces MOND-like acceleration scale

**Status (2024-2025)**: Mixed reception.
- Inconsistent with some dwarf galaxy rotation curves (Pardo 2017)
- Requires decreased mass-to-light ratios, in tension with other estimates
- Active research continues but no major new developments

### 5.2 Oppenheim's Stochastic Gravity (2024, arXiv:2402.19459)

**Mechanism**: Spacetime metric is classical + stochastic. In low-acceleration regime, gravitational acceleration variance acts as entropic force.

**Key claim**: Can explain galactic rotation curves without dark matter.

**HOWEVER**: Critical analysis (arXiv:2404.13037) shows the theory "does not in fact lead to modified Newtonian dynamics as claimed." The fluctuation spectrum differs from observations.

### 5.3 Relativistic MOND from Entropic Gravity (2025, arXiv:2511.05632)

**Authors**: Rostami, Rezazadeh, Rostampour

**Mechanism**: Temperature-dependent corrections to equipartition law on holographic screen + Unruh relation --> modified Einstein equations with thermal corrections.

**Results**: Fits NGC 3198 rotation curve, RMOND performs as well as dark matter model at r > 20 kpc.

**Connection to our work**: Uses Unruh temperature, which connects to modular flow and Tomita-Takesaki theory -- same mathematical structure as Dorau-Much.

### 5.4 Tsallis-Renyi Entropy --> MOND (2025, arXiv:2505.03061)

**Result**: Modified Renyi entropy derived from Tsallis entropy --> MOND-like force law + Bekenstein bound + Landauer principle for black holes.

**KEY CONNECTION**: This explicitly uses modified ENTROPY (not energy) to derive MOND. If Tsallis entropy = non-extensive QRE, this supports our thesis that modified Sigma --> rotation curves.

### 5.5 Smolin: MOND as Regime of Quantum Gravity (2017, arXiv:1704.00780)

**Result**: Below deSitter temperature, equivalence principle weakens --> ratio of gravitational/inertial mass becomes environment-dependent --> MOND in classical limit.

**Interesting**: Connects MOND acceleration a_0 to cosmological constant: a_0 ~ cH_0 ~ c^2(Lambda/3)^{1/2}

### 5.6 Bianconi's Gravity from Entropy (2025, PRD 111, 066001)

**Mechanism**: Entropic action = QRE between spacetime metric (as density matrix) and matter-induced metric. Modified Einstein equations reduce to standard in low-coupling limit.

**Key claim**: The G-field might be a candidate for dark matter.

**Predictions**: Predicts small positive cosmological constant consistent with observations.

**CONNECTION TO OUR WORK**: Bianconi's action IS literally Sigma = D(g_spacetime || g_matter). If this has non-trivial galactic-scale solutions, they would be our rotation curves. **This is potentially the most direct connection.**

### 5.7 Information as Dark Matter/Energy (Vagnozzi direction, 2025)

**"Evidence for Dark Energy Driven by Star Formation: Information Dark Energy?"**

Key idea: Landauer principle --> each bit of information has energy kT ln(2). The total information content of the universe gives the right order of magnitude for dark energy. Mass of a bit replaces Planck mass --> reduces cosmological constant discrepancy by 120 orders of magnitude.

### 5.8 Assessment: Landscape of Approaches

| Approach | Mechanism | Fits Data? | Connection to Sigma? |
|----------|-----------|------------|---------------------|
| Kumar (running G) | RG flow, eta=1 | Yes (3 galaxies) | Indirect (RG = DPI) |
| RGGR (Gubitosi+) | Running G from AS | Yes (100 galaxies) | Indirect |
| Verlinde | Entanglement entropy volume/area | Partial | Moderate |
| Oppenheim | Stochastic metric | Disputed | Weak |
| Rostami+ | Debye corrections | 1 galaxy | Moderate (Unruh) |
| Tsallis-Renyi MOND | Modified entropy | Not fitted | **Strong** |
| Smolin | QG regime below T_dS | Qualitative | Moderate |
| Bianconi | QRE as action | Untested | **STRONGEST** |

---

## 6. SPARC Database Details

### 6.1 Original SPARC (Lelli, McGaugh & Schombert 2016)

**Reference**: arXiv:1606.09251, AJ 152, 157 (2016)
**URL**: https://astroweb.case.edu/SPARC/
**Content**: 175 nearby disk galaxies with:
- 3.6 micron Spitzer photometry (surface brightness profiles)
- High-quality H I/H-alpha rotation curves
- Newtonian mass models
- Bulge-disk decompositions

### 6.2 Data Files Available

1. **Table1.mrt**: Galaxy sample with rotation curve references
2. **Photometric Profiles**: Surface brightness in mag/arcsec^2 (zip)
3. **Bulge-Disk Decompositions**: (zip)
4. **Newtonian Mass Models**: Rotation curves + baryonic contributions + inclination-corrected density profiles (zip)
5. **WISE Stellar Masses**: Alternative mass estimates
6. **Scaling Relations**: BTFR, Central Surface Density, RAR
7. **Dark Matter Halo Fits**: NFW and other profiles (with MCMC chains)

### 6.3 Standard Fitting Procedure

From Li et al. (2018, arXiv:1803.00022):

1. Total velocity: v^2_tot = Upsilon_disk * v^2_disk + Upsilon_bulge * v^2_bulge + v^2_gas + v^2_DM (or modified gravity term)
2. Parameters fitted via MCMC: Upsilon_disk, Upsilon_bulge, (galaxy distance, inclination as nuisance)
3. Standard mass-to-light ratios: Upsilon_disk ~ 0.5 M_sun/L_sun, Upsilon_bulge ~ 0.7 M_sun/L_sun
4. Residual scatter: 0.057 dex (13%) around RAR

### 6.4 The Radial Acceleration Relation (RAR)

**McGaugh, Lelli & Schombert (2016, PRL 117, 201101, arXiv:1609.05917)**

The observed centripetal acceleration g_obs correlates tightly with the baryonic Newtonian acceleration g_bar:

```
g_obs = g_bar / (1 - e^{-sqrt(g_bar/g_dagger)})
```

where g_dagger = a_0 = 1.20 x 10^{-10} m/s^2.

**Observed scatter**: 0.13 dex over 2693 points in 153 galaxies.

**KEY FOR PAPER 3**: Any theory must reproduce this relation. The logarithmic correction from running G naturally does this because the correction becomes important precisely when g_bar << a_0.

### 6.5 BIG-SPARC (2024, arXiv:2411.13329)

**Authors**: Haubner, Lelli et al.
**Content**: ~4000 galaxies (20x larger than original SPARC)
- HI datacubes from public archives (APERTIF, ASKAP, ATCA, GMRT, MeerKAT, VLA, WSRT)
- NIR photometry from WISE
- Homogeneously derived rotation curves, surface brightness profiles, mass models

**Status**: Announced November 2024. Data may not be fully public yet.

---

## 7. SYNTHESIS: The Paper 3 Argument

### 7.1 Core Thesis

The quantum relative entropy drop Sigma = D(rho||sigma) - D(N(rho)||N(sigma)), when computed using the full quantum-corrected gravitational channel N_grav(k), naturally produces:

1. **Small scales (r << r_c)**: Standard Newton gravity (Sigma ~ r_s/r as in Paper 2)
2. **Large scales (r >> r_c)**: Logarithmic correction giving flat rotation curves
3. **The transition**: Controlled by a UNIVERSAL crossover scale k_* ~ 2.7 x 10^{-2} kpc^{-1}

The key insight: the gravitational channel N_grav has scale-dependent "noise" because of quantum corrections (running G). At large distances, the channel is "noisier" (more information loss), producing extra Sigma, which manifests as extra gravitational potential.

### 7.2 The Logical Chain

```
Paper 1: Sigma = 0 <==> perfect recovery (Petz map) <==> time reversal
Paper 2: Sigma_grav = r_s/r for Schwarzschild (semiclassical Einstein)
Paper 3: Sigma_grav(k) = r_s(k)/r with r_s(k) = 2G(k)M/c^2 and G(k) running
    --> At k << k_*: G(k) ~ G_N(k_*/k) gives extra ln(r) term
    --> Rotation curves are flat
Paper 4: Applications (QEC, cosmology, time)
```

### 7.3 What We Need to Prove/Argue

1. **The gravitational channel at scale k**: Define N_grav(k) precisely
   - Option A: RG-improved Hamiltonian H(k) with G(k)
   - Option B: Effective field theory one-loop corrected propagator
   - Option C: Modular channel with scale-dependent modular Hamiltonian

2. **Sigma(k) = Sigma_Newton(k) + delta_Sigma_quantum(k)**:
   - Show that delta_Sigma from quantum corrections gives logarithmic potential
   - Show that the coefficient is 2/pi (universal, matching Kumar)

3. **eta = 1 from information theory**:
   - This is the hardest step
   - Possible argument: eta = 1 is the ONLY value where the DPI bound is saturated by a recoverable channel in the IR (marginal = boundary between recoverable and non-recoverable)
   - Alternative: eta = 1 follows from conformal invariance of the IR fixed point

4. **Universal acceleration scale a_0**:
   - Show that a_0 ~ cH_0 emerges from the QRE framework
   - Connection: the de Sitter horizon entropy gives the IR cutoff for the gravitational channel

### 7.4 Strengths of This Approach

1. **NOT adding new physics**: Same Sigma, different channel
2. **Naturally scale-dependent**: Small scales = Newton, large scales = MOND-like
3. **Connected to established QI results**: DPI, Petz recovery, c-theorem
4. **Builds on Papers 1 & 2**: Direct extension of the framework
5. **Testable**: SPARC data with 175 galaxies, soon BIG-SPARC with 4000

### 7.5 Weaknesses/Gaps

1. **CMB problem**: Running G alone cannot explain CMB acoustic peaks (Kumar's limitation). Need to address this honestly.
2. **eta = 1 derivation**: Not yet derived from QRE first principles -- this is the main gap.
3. **Cluster scales**: Bullet Cluster and gravitational lensing not addressed.
4. **Cosmological perturbation theory**: Structure formation requires more than just modified static potential.
5. **The channel definition**: Defining N_grav(k) rigorously at galactic scales is non-trivial.

### 7.6 Honest Novelty Assessment

| Claim | Novelty | Support |
|-------|---------|---------|
| Sigma from QRE gives gravity (semiclassical) | Low (Dorau-Much 2025) | PRL |
| RG flow = DPI (monotonicity) | Low (Casini-Huerta 2017) | JHEP |
| Running G fits rotation curves | Low (Kumar, Gubitosi) | PLB, SPARC |
| **Sigma with quantum-corrected channel gives running G** | **HIGH** | New connection |
| **Same Sigma formula, different scales, different physics** | **HIGH** | New interpretation |
| **QRE framework unifying Newton, rotation curves, cosmology** | **VERY HIGH** | New synthesis |

---

## 8. Key References (Complete List)

### Primary References for Paper 3

1. Kumar (2025), arXiv:2509.05246, PLB 871, 140008 -- Marginal IR running, eta=1, rotation curves
2. Gubitosi, Piattella & Casarini (2024), arXiv:2403.00531 -- SPARC 100-galaxy RGGR fits
3. Dorau & Much (2025), arXiv:2510.24491, PRL -- QRE --> Einstein equations
4. Casini & Huerta (2017), arXiv:1611.00016, JHEP -- QRE, RG flow, c-theorem
5. Lelli, McGaugh & Schombert (2016), arXiv:1606.09251, AJ 152, 157 -- SPARC database
6. McGaugh, Lelli & Schombert (2016), arXiv:1609.05917, PRL 117, 201101 -- RAR
7. Bianconi (2025), arXiv:2408.14391, PRD 111, 066001 -- Gravity from entropy/QRE
8. Reuter & Weyer (2004), hep-th/0410117, PRD 70, 124028 -- Running G, rotation curves

### Supporting References

9. Rodrigues, Letelier, Fabris, Shapiro (2009), arXiv:0911.4967 -- RGGR original
10. Rodrigues et al. (2012), arXiv:1209.0504 -- Disk/elliptical galaxies in RGGR
11. Verlinde (2016), arXiv:1611.02269 -- Emergent gravity and dark universe
12. Smolin (2017), arXiv:1704.00780, PRD 96, 083523 -- MOND from quantum gravity
13. Jacobson (1995), gr-qc/9504004 -- Thermodynamics of spacetime
14. Jacobson (2015), arXiv:1505.04753, PRL 116, 201101 -- Entanglement equilibrium
15. Donoghue (2020), Frontiers in Physics, 10.3389/fphy.2020.00056 -- Critique of AS
16. Bonanno, Eichhorn et al. (2020), Frontiers in Physics, 10.3389/fphy.2020.00269 -- Critical reflections on AS
17. Li et al. (2018), arXiv:1803.00022 -- RAR fitting to individual SPARC galaxies
18. Haubner, Lelli et al. (2024), arXiv:2411.13329 -- BIG-SPARC
19. Oppenheim & Russo (2024), arXiv:2402.19459 -- Stochastic spacetime rotation curves
20. Rostami, Rezazadeh, Rostampour (2025), arXiv:2511.05632 -- Relativistic MOND from entropic gravity
21. arXiv:2505.03061 -- Tsallis-Renyi entropy to MOND-like force law

### Casini-Huerta & c-theorem

22. Casini & Huerta (2004), hep-th/0405111 -- Finite entanglement entropy and c-theorem
23. Casini & Huerta (2007), cond-mat/0610375 -- c-theorem for entanglement entropy
24. Casini & Huerta (2023), arXiv:2303.xxxxx, PRL 130, 111603 -- Entropic theorem in general dimensions

### 2024 Dirac Medal: Casini & Huerta (ICTP, 2024)
"For their insights on quantum entropy in quantum gravity and quantum field theories"

---

## 9. STRATEGIC NOTES

### 9.1 The Single Most Important Connection

**Casini-Huerta + Dorau-Much = Paper 3 foundation**

- Dorau-Much: QRE --> Einstein equations (at fixed scale)
- Casini-Huerta: QRE monotonicity = RG flow irreversibility
- Paper 3: QRE AT EACH SCALE --> scale-dependent Einstein equations --> running G --> rotation curves

This is NOT adding terms. This is saying: "The Einstein equations that follow from QRE at scale k have a G that depends on k, because QRE itself has RG flow."

### 9.2 The Philosophical Point

Standard approach: G is constant, add dark matter
AS approach: G runs, no dark matter (but ad hoc running)
**Our approach**: G runs BECAUSE the gravitational quantum channel processes information differently at different scales. The running IS the DPI applied scale by scale. Sigma naturally changes with scale because the channel changes with scale.

### 9.3 What Makes This Different from Just "Running G"

Running G has been proposed since 2004 (Reuter-Weyer). What's new:

1. **The running is not put in by hand** -- it follows from the channel's QRE at each scale
2. **Sigma = D(rho||sigma) - D(N_k(rho)||N_k(sigma))** gives a UNIFIED framework where:
   - Sigma = 0: flat spacetime, perfect recovery (Paper 1)
   - Sigma ~ r_s/r: Schwarzschild, semiclassical (Paper 2)
   - Sigma ~ r_s/r + quantum corrections: rotation curves (Paper 3)
3. **The universality of eta=1** may follow from QRE saturation properties
4. **The acceleration scale a_0** may emerge from the cosmological horizon (de Sitter temperature = IR cutoff for modular Hamiltonian)

### 9.4 Next Steps

1. **Detailed calculation**: Define N_grav(k) using Dorau-Much framework but at each RG scale
2. **Show Sigma(k) gives Kumar's potential**: This is the central calculation
3. **Fit SPARC data**: Use Kumar/Gubitosi formulas with our reinterpretation
4. **Address CMB**: Honestly discuss what Sigma-based running G can and cannot do
5. **Write Paper 3**: 5-page PRL letter + supplemental

---

## 10. APPENDIX: Radial Acceleration Relation Data

The RAR is the single tightest empirical constraint:

```
g_obs = F(g_bar) = g_bar / (1 - exp(-sqrt(g_bar / a_0)))
```

with a_0 = 1.20 x 10^{-10} m/s^2.

Scatter: 0.13 dex (2693 data points, 153 galaxies).

Any Paper 3 prediction must reproduce this with comparable or better scatter.

The RGGR model (Gubitosi+) achieves this with 3 parameters per galaxy.
Kumar achieves it with 2 parameters (k_*, G_N) but only tested on 3 galaxies.

**Our target**: Reproduce RAR with k_* derivable from Sigma, reducing to 1 effective free parameter (the mass-to-light ratio).
