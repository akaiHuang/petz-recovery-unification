# GfE Solutions Survey: Can Bianconi's Framework Produce the Exponential Metric?

**Research Date:** 2026-03-10
**Agent:** Physics research agent (comprehensive literature survey)
**Purpose:** Determine whether the Gravity from Entropy (GfE) field equations admit the exponential metric g₀₀ = −exp(−R_s/r) as a solution, and survey all related literature

---

## Table of Contents

1. [Executive Summary](#1-executive-summary)
2. [Complete Bianconi GfE Bibliography](#2-complete-bianconi-gfe-bibliography)
3. [Papers Citing/Extending GfE](#3-papers-citingextending-gfe)
4. [The GfE Field Equations: What We Know](#4-the-gfe-field-equations-what-we-know)
5. [CRITICAL FINDING: First Spherically Symmetric GfE Solution (2026)](#5-critical-finding-first-spherically-symmetric-gfe-solution-2026)
6. [Exponential Metric Literature: Complete Survey](#6-exponential-metric-literature-complete-survey)
7. [Phantom Scalar Field and the Exponential Metric](#7-phantom-scalar-field-and-the-exponential-metric)
8. [Modified Gravity Theories Producing the Exponential Metric](#8-modified-gravity-theories-producing-the-exponential-metric)
9. [Entropic/Information-Theoretic Metric Derivations](#9-entropicinformation-theoretic-metric-derivations)
10. [G-Field as Dark Matter: Galactic Scale Predictions](#10-g-field-as-dark-matter-galactic-scale-predictions)
11. [Numerical Methods for Solving GfE Equations](#11-numerical-methods-for-solving-gfe-equations)
12. [Assessment: GfE Solution — Schwarzschild vs Exponential vs New](#12-assessment-gfe-solution--schwarzschild-vs-exponential-vs-new)
13. [Implications for Paper 2 and Paper 3](#13-implications-for-paper-2-and-paper-3)
14. [Complete Bibliography](#14-complete-bibliography)

---

## 1. Executive Summary

### The Core Question
Does the GfE framework produce the exponential metric g₀₀ = −exp(−R_s/r) as its vacuum spherically symmetric solution?

### The Answer (as of March 2026)

**NO — but the story is more interesting than expected.**

The first paper to solve the GfE field equations for static spherically symmetric vacuum appeared on **February 14, 2026** (Thattarampilly, Zheng & Kakkat, arXiv:2602.13694). Their result:

> **The GfE solution is Schwarzschild + perturbative r⁻⁴ corrections:**
> - A(r) = 1 − r_S/r − β r_S²/(48 r⁴)
> - B(r) = (1 − r_S/r)⁻¹ + β r_S²/(12 r⁴)
> - Event horizon shifts to: r_h = r_S + β/(48 r_S) + O(β²)

This is **neither Schwarzschild nor the exponential metric** — it is a perturbative correction to Schwarzschild.

### Key Implications

1. **The GfE does NOT naturally produce the exponential metric.** At perturbative order, it gives Schwarzschild + r⁻⁴ corrections.
2. **However**, this is only a perturbative (small β) result. The full nonlinear GfE equations remain unsolved.
3. **The exponential metric exp(−R_s/r) = 1 − R_s/r + R_s²/(2r²) − ...** This Taylor expansion does NOT match the r⁻⁴ correction pattern found. The exponential metric has r⁻² corrections, not r⁻⁴.
4. **This is important for Paper 2**: If Paper 2 claims Σ_grav = r_s/r → exponential metric, the GfE framework does NOT support this at perturbative level. The GfE gives a different modification of Schwarzschild.

### Likelihood Assessment

| Solution | Likelihood from GfE | Evidence |
|----------|---------------------|----------|
| Schwarzschild (exact) | ❌ Ruled out | r⁻⁴ corrections exist |
| Exponential metric | ❌ Very unlikely | Corrections are r⁻⁴, not r⁻² |
| Schwarzschild + r⁻⁴ corrections | ✅ Confirmed | Thattarampilly et al. 2026 |
| Something new (full nonlinear) | ❓ Unknown | Full solution not yet computed |

---

## 2. Complete Bianconi GfE Bibliography

### 2.1 Core Papers by Ginestra Bianconi

**[B1] "Gravity from Entropy" (the founding paper)**
- Bianconi, G. (2025). Physical Review D 111, 066001
- arXiv: 2408.14391 (submitted Aug 26, 2024; revised Feb 8, 2025)
- DOI: 10.1103/PhysRevD.111.066001
- **Content**: Defines the GfE action S = ∫ √(-g) [Tr(g̃ · ln(G̃⁻¹)) − Λ] d⁴x, derives modified Einstein equations, introduces G-field as Lagrange multiplier, shows emergence of positive cosmological constant, proposes G-field as dark matter candidate.

**[B2] "The Quantum Relative Entropy of the Schwarzschild Black Hole and the Area Law"**
- Bianconi, G. (2025). Entropy 27(3), 266
- arXiv: 2501.09491 (submitted Jan 16, 2025)
- DOI: 10.3390/e27030266
- **Content**: Computes the GQRE for Schwarzschild metric, shows it obeys area law for large R_s. Notes Schwarzschild is an approximate (not exact) GfE solution valid at low coupling/small curvature.
- **Correction published**: Entropy 27(7), 724 (July 2025) — corrects minor errors.

**[B3] "The Thermodynamics of the Gravity from Entropy Theory"**
- Bianconi, G. (2026). arXiv: 2510.22545 (submitted Oct 29, 2025; revised Jan 7, 2026)
- **Content**: Develops GfE thermodynamics. Defines internal energy density ε = 2βΛ_G, k-temperatures θ_k = G_k − 1, k-pressures. Studies FRW cosmological solutions. Shows total GQRE per unit volume is non-increasing while total entropy is non-decreasing. Does NOT solve static vacuum equations.

**[B4] "Beyond Holography: The Entropic Quantum Gravity Foundations of Image Processing"**
- Bianconi, G. (2025). Physical Review E
- arXiv: 2503.14048 (submitted Mar 2025; published in PRE)
- DOI: 10.1103/gkkt-yysp
- **Content**: Shows that the Perona-Malik algorithm for image processing is the gradient flow maximizing the GfE action in its warm-up scenario. Bridge between GfE, machine learning, and brain research.

### 2.2 Status of Bianconi's GfE Program (Summary)

| Paper | Status | Vacuum solution? |
|-------|--------|-----------------|
| B1: Gravity from Entropy | Published PRD | No — states "beyond scope" |
| B2: Schwarzschild & Area Law | Published Entropy | Schwarzschild as approximate solution |
| B3: Thermodynamics | arXiv preprint | FRW cosmology only |
| B4: Beyond Holography | Published PRE | N/A (image processing) |

**Bianconi herself has NOT solved the vacuum spherically symmetric GfE equations.** She has only shown that Schwarzschild is an approximate solution at low coupling.

---

## 3. Papers Citing/Extending GfE

### 3.1 Thattarampilly & Zheng Group (Yangzhou University)

**[TZ1] "Inflation from Entropy"**
- Thattarampilly, U. & Zheng, Y. (2025). European Physical Journal C
- arXiv: 2509.23987 (submitted Sep 28, 2025; published in EPJC)
- **Content**: First external group to work with GfE equations. Solves modified Friedmann equations in GfE framework. Finds inflationary solutions without requiring additional scalar fields. High-entropy branch gives spectral index n_s = 0.962–0.964, consistent with CMB. Published in EPJC.

**[TZ2] "Spherically Symmetric Black Holes in Gravity from Entropy and Spontaneous Emission"** ⭐ CRITICAL
- Thattarampilly, U., Zheng, Y. & Kakkat, V. (2026)
- arXiv: 2602.13694 (submitted Feb 14, 2026)
- **Content**: **First solution of GfE field equations for static spherically symmetric vacuum.** See Section 5 for detailed analysis.

### 3.2 Obidi: Theory of Entropicity

**[O1] "On the Theory of Entropicity (ToE) and Ginestra Bianconi's Gravity from Entropy"**
- Obidi, J.O. (2025). SSRN preprint (Nov 12, 2025)
- Cambridge Open Engage, SSRN 5738123
- **Content**: Claims Bianconi's GfE emerges as a special limiting case within a broader "Theory of Entropicity" framework. The ToE's variational principle (Obidi Action) reduces to Bianconi's GQRE under quadratic expansion near equilibrium. Also claims ToE encompasses Jacobson, Padmanabhan, and Verlinde. **Note: This is a self-published preprint, not peer-reviewed. Treat with caution.**

### 3.3 Other citations

**[G1] Guo, X. & Yuan, F. (2025). "Quantum Effective Dynamics of Papapetrou Spacetime"**
- arXiv: 2506.08821; Published in EPJC
- **Content**: Studies quantum corrections to the Papapetrou (exponential) spacetime using loop quantum gravity methods. Finds quantum effects create a new wormhole throat while the classical throat disappears for very small masses. Not directly about GfE but relevant to quantum corrections of the exponential metric.

---

## 4. The GfE Field Equations: What We Know

### 4.1 The GfE Action

The Lagrangian density is:

```
L_GfE = Tr(g̃ · ln(G̃⁻¹)) − Λ
```

where:
- g̃ is the spacetime metric (treated as an operator)
- G̃ is the "induced metric" from matter fields and curvature
- Λ is a constant (NOT the cosmological constant — the cosmological constant emerges dynamically)
- The trace is over the fiber (zero-form, one-form, two-form) indices

### 4.2 The Modified Einstein Equations

The field equations (Eq. 45 in Bianconi 2025) have the form:

```
β R^{μν} [G₀]⁻¹ − ½ g^{μν} L + B^{μν} = T̃^{μν}
```

where:
- β is a positive coupling constant (small β → Einstein equations)
- B^{μν} contains derivative terms involving both G-field and Ricci tensor
- T̃^{μν} is the modified stress-energy tensor
- [G₀] is a background G-field component

### 4.3 The G-Field Equation

The G-field G̃ enters as a Lagrange multiplier but is promoted to a dynamical field. Its equation of motion relates it to both matter and curvature:

```
Θ̃ = G̃ · g̃⁻¹
```

The G-field "dresses" the metric and enables the separation of the action into gravitational + matter sectors.

### 4.4 Vacuum Behavior

In vacuum (no matter fields), the induced metric depends only on curvature terms. The equations become:

```
R_{(μν)} G⁻¹ − ½ g_{μν} (R_G − 2Λ_G) + D_{(μν)} = 0
```

where Λ_G is the emergent dynamical cosmological constant depending only on G.

**Key point**: In GfE, "vacuum" is NOT simply R_μν = 0 (as in GR). The G-field carries additional degrees of freedom even in vacuum.

---

## 5. CRITICAL FINDING: First Spherically Symmetric GfE Solution (2026)

### Paper Details

**"Spherically Symmetric Black Holes in Gravity from Entropy and Spontaneous Emission"**
- Authors: Udaykrishna Thattarampilly, Yunlong Zheng, Vishnu Kakkat
- arXiv: 2602.13694 (February 14, 2026)
- Affiliation: Yangzhou University; Luleå University of Technology

### Method

The authors:
1. Start from the full GfE field equations
2. Assume a static, spherically symmetric ansatz: ds² = −A(r)dt² + B(r)dr² + r²dΩ²
3. Solve perturbatively in the coupling parameter β around the Schwarzschild solution

### Results

**The metric functions receive r⁻⁴ corrections:**

```
A(r) = 1 − r_S/r − β r_S² / (48 r⁴)

B(r) = (1 − r_S/r)⁻¹ + β r_S² / (12 r⁴)
```

**Event horizon shifts:**

```
r_h = r_S + β/(48 r_S) + O(β²)
```

### Physical Predictions

1. **S2 star orbital precession**: Consistent with observations for −1 < β < 1
2. **EHT shadow diameter**: Consistent for −9.04 r_S² ≤ β ≤ 18.63 r_S²
3. **Spontaneous mass loss** (in Lemaître coordinates): Ṁ ≈ −β/24 − 0.17β/M²
   - This mimics Hawking radiation from a purely classical GfE framework!
   - Temperature scaling: T ∝ M⁻¹/² (different from Hawking's T_H ∝ M⁻¹)

### Comparison with Exponential Metric

The exponential metric has the Taylor expansion:

```
exp(−R_s/r) = 1 − R_s/r + R_s²/(2r²) − R_s³/(6r³) + R_s⁴/(24r⁴) − ...
```

The GfE correction is:

```
A(r) = 1 − R_s/r + 0 · r⁻² + 0 · r⁻³ − β R_s²/(48 r⁴)
```

**These are fundamentally different:**
- Exponential metric: corrections at ALL orders (r⁻², r⁻³, r⁻⁴, ...)
- GfE perturbative: correction ONLY at r⁻⁴ order

**Conclusion: The perturbative GfE solution does NOT converge to the exponential metric.** The r⁻² and r⁻³ terms are absent.

### Caveat

This is a **perturbative** result (small β). The full nonlinear GfE equations could in principle give a very different solution. However:
- If the full solution were the exponential metric, one would expect r⁻² corrections at leading order in perturbation theory
- The absence of r⁻² corrections is strong evidence against the exponential metric being the GfE solution
- The full nonlinear solution could still be "something new" that is neither Schwarzschild nor exponential

---

## 6. Exponential Metric Literature: Complete Survey

### 6.1 Classical References

**[P1] Papapetrou, A. (1954)**
- "Eine Theorie des Gravitationsfeldes mit einer Feldfunktion"
- Math. Nachr. 12, 129–141
- **First appearance** of the exponential metric. Proposed as a simpler alternative to GR.

**[Y1] Yilmaz, H. (1958, 1971)**
- "New Approach to General Relativity" Phys. Rev. 111, 1417 (1958)
- "New Theory of Gravitation" Phys. Rev. Lett. 27, 1399 (1971)
- Re-derived the exponential metric from a field-theoretic approach. Correctly predicts three classical tests of GR (perihelion precession, light deflection, gravitational redshift).

**[R1] Rosen, N. (1973)**
- Bimetric gravity theory
- The exponential metric appears as a solution in Rosen's bimetric framework.

### 6.2 Modern Exponential Metric Papers

**[MM1] Makukov, M.A. & Mychelkin, E.G. (2020). "Triple Path to the Exponential Metric"**
- Foundations of Physics 50, 1346–1355
- arXiv: 2009.08655
- **Key result**: Three independent families of static solutions (Fisher 1948, Janis-Newman-Winicour 1968, Xanthopoulos-Zannias 1989) ALL reduce to the exponential metric when scalar charge = central mass. This suggests a universal character of the scalar background.

**[MM2] Makukov, M.A. & Mychelkin, E.G. (2023). "Axially Symmetric Exponential Metric"**
- arXiv: 2309.02874
- Published in New Astronomy Reviews / Dark Universe
- **Content**: Constructs a two-parameter family of Axially Symmetric Exponential Metric (ASEM) solutions within GR, generalizing the spherically symmetric exponential metric. Interpolates between YEM and Curzon-Chazy geometry.

**[MM3] Makukov, M.A., Mychelkin, E.G. & Suliyeva, G. (2023). "Rotation in vacuum and scalar background: are there alternatives to Newman-Janis algorithm?"**
- arXiv: 2301.08118
- **Content**: Applies Newman-Janis algorithm to the Papapetrou antiscalar spacetime. Produces a rotational (Kerr-like) metric devoid of horizons and ergospheres.

**[MM4] Mychelkin, E.G., Makukov, M.A. & Suliyeva, G. (2024). "On the weak and strong field effects in antiscalar background"**
- General Relativity and Gravitation (2024)
- arXiv: 2403.11610
- **Content**: Compares vacuum (Schwarzschild) and antiscalar (exponential) metrics. Weak-field perihelion shift difference is observationally imperceptible even for S-cluster stars. Strong-field shadow effect is the most promising observational test.

**[MM5] Mychelkin, E.G., Suliyeva, G. & Makukov, M.A. (2025/2026). "The exponential metric: traversable wormhole and possible identification of scalar background"**
- arXiv: 2510.15391 (submitted Oct 2025; revised Jan 8, 2026)
- **Content**: Identifies the topological Gauss-Bonnet invariant sign swap at the wormhole throat (r = M). Associates the Ricci scalar at r = M/2 with extremal thermodynamic characteristics. Connects to the scalar background identification.

**[Bo1] Boonserm, P., Ngampitipan, T., Simpson, A. & Visser, M. (2018). "The exponential metric represents a traversable wormhole"**
- Physical Review D 98, 084048
- arXiv: 1805.03781
- **Key result**: The exponential metric ds² = −e^(−2m/r) dt² + e^(+2m/r){dr² + r²dΩ²} has no horizons and is actually a traversable wormhole.

**[N1] Nandkishore, A.K. et al. (2024). "A new class of traversable wormhole metrics"**
- European Physical Journal C 84, 1063
- **Content**: Generalizes the exponential wormhole metric family with exponential shape and redshift functions. Studies photon sphere, ISCO, QNMs, gravitational entropy.

**[NS1] Nath, P.P. & Sarma, D. (2025). "Exponential wormhole, its quasinormal modes, and study in logarithmic f(R) gravity"**
- Annals of Physics 479, 170067 (May 2025)
- **Content**: QNM analysis of exponential wormhole using 3rd-order WKB. Studies in logarithmic f(R) gravity with cubic gravity term. Finds stability through TOV equation.

**[E1] Ernazarov, K.K. (2025). "The Yilmaz-Rosen and JNW metric solutions in the scalar-Einstein-Gauss-Bonnet 4d gravitational model"**
- arXiv: 2510.21625 (Oct 24, 2025)
- **Key result**: In sEGB framework, the Yilmaz-Rosen metric requires a phantom-like scalar field with vanishing potential. All energy conditions are violated. The energy-momentum tensor implies exotic matter.

**[E2] Ernazarov, K.K. (2025). "The asymptotically Schwarzschild-like metric solutions"**
- arXiv: 2511.13471 (Nov 17, 2025)
- **Content**: The exponential metric is not a vacuum GR solution (Ricci tensor is non-zero). Computes photon sphere radius, ISCO, Regge-Wheeler potential, shadow size, and compares with Schwarzschild. Notes the exponential metric remains relevant for beyond-GR theories despite being ruled out by solar system tests as a GR vacuum solution.

**[M1] Manna, T., Rahman, F. & Chowdhury, T. (2023). "Strong lensing in exponential wormhole spacetimes"**
- New Astronomy Reviews 97, 101695
- **Content**: Computes deflection angle in exponential metric using Bozza's method. Studies relativistic image formation.

**[GY1] Guo, X. & Yuan, F. (2025). "Quantum Effective Dynamics of Papapetrou Spacetime"**
- EPJC (2025)
- arXiv: 2506.08821
- **Content**: Loop quantum gravity corrections to the Papapetrou (exponential) metric. Quantum effects create a new wormhole throat; classical throat disappears for extremely small masses.

---

## 7. Phantom Scalar Field and the Exponential Metric

### 7.1 The Key Connection: φ = M/r

The exponential metric g₀₀ = −exp(−2M/r) is an **exact solution** of Einstein's equations coupled to a massless **antiscalar (phantom) field** with φ = M/r.

An antiscalar field has a kinetic term with the "wrong" sign:

```
L_phantom = +½ g^{μν} ∂_μ φ ∂_ν φ   (note: + instead of −)
```

This violates the null energy condition (NEC), which is precisely why the exponential metric can be a traversable wormhole (wormholes require NEC violation by the topological censorship theorem).

### 7.2 The Three Paths (Makukov & Mychelkin 2020)

Three independent families of scalar field solutions in GR all converge to the exponential metric:

1. **Fisher (1948)**: Static spherically symmetric massless scalar field solutions
2. **Janis-Newman-Winicour (1968)**: Most general static spherically symmetric solution of Einstein-massless-scalar
3. **Xanthopoulos-Zannias (1989)**: Independent derivation

**All three give the exponential metric when scalar charge q = mass M.**

This is a remarkable convergence: three different mathematical approaches, spanning 41 years, all arriving at the same metric.

### 7.3 Ernazarov (2025): sEGB Reconstruction

Ernazarov's reconstruction of the Yilmaz-Rosen metric in the scalar-Einstein-Gauss-Bonnet framework confirms:
- The scalar field is phantom-like
- The potential vanishes (U = 0)
- All energy conditions are violated
- The energy-momentum tensor requires exotic matter

### 7.4 Relevance for Paper 2

In Paper 2, the claim is that Σ_grav = r_s/r is a phantom scalar field. This is **mathematically correct**: the exponential metric IS the solution of GR + phantom scalar field with φ = M/r. The three independent derivations (Makukov & Mychelkin) confirm this universality.

However, the question is whether this phantom scalar field emerges naturally from a more fundamental theory (like GfE), or must be put in by hand.

---

## 8. Modified Gravity Theories Producing the Exponential Metric

### 8.1 Theories Where the Exponential Metric Appears

| Theory | How it appears | Reference |
|--------|---------------|-----------|
| GR + phantom scalar (φ = M/r) | Exact solution | Fisher 1948; JNW 1968; Makukov-Mychelkin 2020 |
| Yilmaz gravity | Fundamental solution | Yilmaz 1958, 1971 |
| Rosen bimetric gravity | Solution | Rosen 1973 |
| sEGB with phantom scalar | Reconstruction | Ernazarov 2025 |
| Logarithmic f(R) gravity | Consistent wormhole | Nath & Sarma 2025 |
| Loop quantum gravity (corrected) | Quantum-corrected version | Guo & Yuan 2025 |

### 8.2 Theories Where it Does NOT Appear

| Theory | What appears instead | Reference |
|--------|---------------------|-----------|
| Standard GR (vacuum) | Schwarzschild | Birkhoff's theorem |
| GfE (perturbative) | Schwarzschild + r⁻⁴ corrections | Thattarampilly et al. 2026 |
| Brans-Dicke | Different corrections depending on ω | Various |
| Standard f(R) gravity | Modified Schwarzschild | Various |

### 8.3 Observational Distinguishability

From Mychelkin, Makukov & Suliyeva (2024):
- **Weak field**: Exponential metric vs Schwarzschild: perihelion shift difference is observationally imperceptible, even for S-cluster stars
- **Strong field**: Shadow diameter differs by detectable amount (most promising test)
- **Wormhole signature**: No event horizon → different ringdown signal

From Ernazarov (2025):
- Photon sphere radius: differs from Schwarzschild
- ISCO radius: differs from Schwarzschild
- Shadow size: larger than Schwarzschild shadow

---

## 9. Entropic/Information-Theoretic Metric Derivations

### 9.1 Comparison of Entropic Gravity Programs

| Program | Derives Einstein equations? | Derives specific metric? | Dark matter prediction? |
|---------|---------------------------|-------------------------|----------------------|
| Jacobson (1995) | Yes (from horizon thermodynamics) | No | No |
| Verlinde (2011/2016) | Yes (Newtonian limit + corrections) | No specific metric | Yes (apparent DM from volume-law entanglement) |
| Padmanabhan | Yes (from extremizing entropy functional) | No | No |
| Bianconi GfE (2025) | Yes (modified Einstein eq.) | Yes (Schwarzschild + r⁻⁴) | Yes (G-field) |
| Dorau & Much (2025) | Yes (semiclassical Einstein eq.) | No | No |

### 9.2 Dorau & Much (2025, PRL)

**"From Quantum Relative Entropy to the Semiclassical Einstein Equations"**
- arXiv: 2510.24491
- Published: Physical Review Letters
- DOI: 10.1103/lmq8-nsty

Key contribution: Using modular theory, they establish that the relative entropy between vacuum state and coherent excitations of a scalar QFT on a bifurcate Killing horizon equals the energy flux. Under the Bekenstein-Hawking area formula, this gives the semiclassical Einstein equations.

**This is a quantum field-theoretic generalization of Jacobson's 1995 thermodynamic derivation**, replacing classical thermodynamic entropy with Araki-Uhlmann quantum relative entropy.

### 9.3 Rostami, Rezazadeh & Rostampour (2025)

**"Relativistic MOND Theory from Modified Entropic Gravity"**
- arXiv: 2511.05632 (Nov 7, 2025)
- **Content**: Derives relativistic MOND from entropic gravity with Debye-like temperature corrections. Tests against NGC 3198 rotation curve. RMOND fits better than baryons-only prediction.

### 9.4 Key Insight: No Entropic Program Derives the Exponential Metric

None of the entropic gravity programs (Jacobson, Verlinde, Padmanabhan, Bianconi, Dorau-Much) produce the exponential metric. They all derive Einstein equations (or modified Einstein equations), and the specific metric solutions still depend on boundary conditions and matter content.

**The exponential metric requires exotic matter (phantom scalar field).** No entropic gravity program naturally generates this exotic matter content.

---

## 10. G-Field as Dark Matter: Galactic Scale Predictions

### 10.1 Bianconi's Claim

Bianconi suggests (in B1) that "the G-field might be a candidate for dark matter." The G-field is the Lagrange multiplier field that dresses the metric and depends on both matter and geometry.

### 10.2 Current Status: No Galactic-Scale Calculation Exists

**Nobody has computed the G-field profile for a galaxy.** The G-field dark matter claim is:
- A qualitative suggestion by Bianconi
- Not yet quantified with any galactic-scale calculation
- Not compared to rotation curve data
- Not compared to MOND or Verlinde predictions

### 10.3 What Would Be Needed

To test whether the G-field acts as dark matter at galactic scales:

1. Solve the GfE field equations for a non-vacuum case with a galactic mass distribution
2. Compute the G-field profile G̃(r) around a galaxy
3. Determine whether the G-field's gravitational effect matches the observed rotation curve
4. Compare with MOND and Verlinde predictions

### 10.4 Comparison with Verlinde's Approach

Verlinde (2016, arXiv:1611.02269):
- Derives apparent dark matter from volume-law entanglement entropy
- Makes specific predictions for rotation curves: M_D(r) = (c² H₀/6G) a_B r
- Successfully compared to SPARC data (Gubitosi et al. 2024)
- At 5.2σ better than MOND (Ghari & Haghi 2026)

**The GfE has not reached this level of quantitative prediction for galactic scales.**

### 10.5 Bridge to Paper 3

For Paper 3's framework (Σ → running G at galactic scales), the GfE's G-field could provide a concrete mechanism IF:
- The G-field naturally produces a radial dependence that mimics dark matter
- This can be related to Σ_grav

However, this connection is currently speculative. The perturbative GfE solution (Thattarampilly et al. 2026) gives r⁻⁴ corrections, which fall off too fast to explain flat rotation curves (which require r⁻¹ to r⁰ corrections at large r).

---

## 11. Numerical Methods for Solving GfE Equations

### 11.1 Why Numerical Methods Are Needed

The full nonlinear GfE equations have NOT been solved exactly for any spacetime (only perturbatively and in cosmological approximation). Key challenges:
- The equations couple the metric g_μν and the G-field G̃ simultaneously
- The G-field equation is nonlinear (involves logarithms of operators)
- Even in the spherically symmetric case, the system is a coupled nonlinear ODE system

### 11.2 Available Methods

**Shooting method:**
- Standard for two-point boundary value problems in GR (e.g., neutron star structure)
- Integrate from r → ∞ (asymptotic flatness) and from r → r_h (horizon regularity)
- Match at an intermediate point
- Suitable for the GfE equations IF proper asymptotic boundary conditions are known

**Relaxation method:**
- Discretize the ODE on a grid and solve the nonlinear algebraic system iteratively
- More robust than shooting for stiff problems
- Standard implementation: Newton-Raphson with band-diagonal Jacobian

**Spectral methods:**
- Errors decrease exponentially with number of grid points
- Standard in numerical relativity (e.g., KADATH, SpEC codes)
- Best for smooth solutions (may have issues near horizons)

### 11.3 Concrete Proposal for Solving Full Nonlinear GfE

1. **Start from Thattarampilly et al.'s perturbative solution** as initial guess
2. Use a **relaxation method** with the full nonlinear GfE equations
3. Gradually increase β from 0 (Schwarzschild) to finite values (parameter continuation)
4. Check whether the solution develops r⁻² corrections (would suggest convergence toward exponential metric) or maintains the r⁻⁴ pattern

**Required inputs:**
- Explicit form of the GfE vacuum equations in spherical symmetry (available from Thattarampilly et al.)
- Boundary conditions at infinity (asymptotic flatness) and at the horizon (or, if the horizon disappears, at the wormhole throat)
- The G-field profile as a function of r

---

## 12. Assessment: GfE Solution — Schwarzschild vs Exponential vs New

### 12.1 Evidence Summary

**For Schwarzschild (exact):**
- ❌ Ruled out: Thattarampilly et al. show r⁻⁴ corrections exist for any β ≠ 0

**For Exponential Metric:**
- ❌ Perturbative evidence strongly against it:
  - Exponential metric has r⁻² corrections at leading order
  - GfE gives r⁻⁴ corrections at leading order
  - These are inconsistent unless there is a miraculous cancellation at higher orders (extremely unlikely)
- ❌ The exponential metric requires exotic (phantom) matter in GR; the GfE vacuum equations have no such source
- ❌ No entropic gravity program has ever produced the exponential metric

**For Schwarzschild + r⁻⁴ corrections (perturbative GfE):**
- ✅ Explicitly derived by Thattarampilly et al. 2026
- ✅ Consistent with observational constraints (S2 star, EHT)
- ✅ Produces Hawking-like radiation from classical framework
- ❓ Only valid for small β; full nonlinear solution unknown

**For something entirely new (full nonlinear GfE):**
- ❓ Possible but unknown
- The full nonlinear solution could differ qualitatively from the perturbative expansion
- Possible scenarios: the horizon disappears (wormhole?), new singularity structure, etc.

### 12.2 Overall Assessment

**Likelihood table:**

| Solution | Probability | Reasoning |
|----------|------------|-----------|
| Schwarzschild (exact) | 0% | Explicitly ruled out |
| Exponential metric | ~2% | Leading-order corrections are wrong power |
| Schwarzschild + polynomial corrections | ~70% | Supported by perturbative calculation |
| Something qualitatively new | ~28% | Nonlinear effects could change picture |

### 12.3 What This Means for Paper 2

Paper 2 claims that Σ_grav = r_s/r leads to the exponential metric g₀₀ = exp(−r_s/r). This claim is:

1. **Mathematically valid as a standalone argument**: If you define Σ_grav = −ln(g₀₀), then g₀₀ = exp(−Σ_grav) follows by definition.

2. **Supported by the phantom scalar field connection**: The exponential metric is the universal attractor of scalar field solutions when q = M (Makukov & Mychelkin 2020). If Σ_grav behaves like a phantom scalar field, the exponential metric follows.

3. **NOT supported by the GfE framework**: The GfE equations do not produce the exponential metric at perturbative level. If Paper 2 invokes Bianconi's GfE as a foundational justification, this is problematic.

**Recommendation for Paper 2:**
- Justify the exponential metric via the phantom scalar field / Makukov-Mychelkin universality argument
- Do NOT claim GfE produces the exponential metric (it doesn't, based on current evidence)
- Cite GfE as an independently interesting entropic gravity framework, but note the GfE vacuum solution differs from the exponential metric

---

## 13. Implications for Paper 2 and Paper 3

### 13.1 Paper 2 Implications

**The good news:**
- The exponential metric is mathematically well-supported as a solution of GR + phantom scalar field
- Three independent derivations converge on it (Makukov-Mychelkin universality)
- Σ_grav = r_s/r → exp(−r_s/r) is a clean mathematical construction
- The no-horizon, traversable wormhole nature is established (Boonserm et al. 2018)

**The bad news:**
- GfE does NOT produce it (Thattarampilly et al. 2026)
- It requires exotic (phantom) matter — violates energy conditions
- Weak-field tests cannot distinguish it from Schwarzschild (Mychelkin et al. 2024)
- Only strong-field shadow measurements can potentially test it

**Strategy for Paper 2:**
- Present the exponential metric as emerging from Σ_grav = −ln(g₀₀) (information-theoretic definition)
- Use the Makukov-Mychelkin universality as supporting evidence
- Cite the three first-principles routes (modular flow, gravitational Landauer, quantum channel)
- Acknowledge the GfE result as a different (competing) prediction for the entropic gravity vacuum solution
- Propose the comparison between exponential metric and GfE r⁻⁴ solution as a testable prediction

### 13.2 Paper 3 Implications

**The G-field dark matter connection is currently too weak:**
- No galactic-scale G-field calculation exists
- The r⁻⁴ correction falls off too fast for rotation curves
- Verlinde's approach is more developed for galactic scales

**Better route for Paper 3:**
- Use the quantum channel / running G approach (Kumar 2025)
- Connect to Verlinde via volume-law entanglement
- Do NOT rely on GfE's G-field until galactic-scale calculations exist

### 13.3 The "Competition" Between Predictions

The GfE and the exponential metric now represent **two distinct predictions** for what entropic gravity does in the strong field:

| Observable | Exponential metric | GfE (perturbative) |
|-----------|-------------------|---------------------|
| Event horizon | NO (traversable wormhole) | YES (shifted by β) |
| g₀₀ at r = 2M | exp(−1) ≈ 0.368 | 1 − 1/2 − β/(48·16M²) |
| Photon sphere | Different from Schwarzschild | Close to Schwarzschild |
| QNMs | Wormhole echoes | Black hole QNMs + corrections |
| Hawking radiation | None (no horizon) | Modified (T ∝ M⁻¹/²) |

**This is actually good for Paper 2**: it makes the exponential metric a **testable alternative** to the GfE prediction.

---

## 14. Complete Bibliography

### Bianconi GfE Papers
1. Bianconi, G. (2025). "Gravity from entropy." PRD 111, 066001. arXiv:2408.14391
2. Bianconi, G. (2025). "The quantum relative entropy of the Schwarzschild black hole and the area law." Entropy 27(3), 266. arXiv:2501.09491
3. Bianconi, G. (2026). "The Thermodynamics of the Gravity from Entropy Theory." arXiv:2510.22545
4. Bianconi, G. (2025). "Beyond holography: the entropic quantum gravity foundations of image processing." PRE. arXiv:2503.14048
5. Bianconi, G. (2025). Correction to [2]. Entropy 27(7), 724

### GfE Extensions
6. Thattarampilly, U. & Zheng, Y. (2025). "Inflation from entropy." EPJC. arXiv:2509.23987
7. **Thattarampilly, U., Zheng, Y. & Kakkat, V. (2026). "Spherically symmetric black holes in Gravity from Entropy and spontaneous emission." arXiv:2602.13694** ⭐
8. Obidi, J.O. (2025). "On the Theory of Entropicity and Bianconi's GfE." SSRN 5738123

### Exponential Metric: Classical
9. Papapetrou, A. (1954). "Eine Theorie des Gravitationsfeldes mit einer Feldfunktion." Math. Nachr. 12, 129
10. Yilmaz, H. (1958). "New Approach to General Relativity." Phys. Rev. 111, 1417
11. Yilmaz, H. (1971). "New Theory of Gravitation." PRL 27, 1399
12. Rosen, N. (1973). Bimetric gravity theory
13. Fisher, I.Z. (1948). Zh. Eksp. Teor. Fiz. 18, 636

### Exponential Metric: Modern
14. Makukov, M.A. & Mychelkin, E.G. (2020). "Triple Path to the Exponential Metric." Found. Phys. 50, 1346. arXiv:2009.08655
15. Makukov, M.A. & Mychelkin, E.G. (2023). "Axially Symmetric Exponential Metric." arXiv:2309.02874
16. Makukov, M.A., Mychelkin, E.G. et al. (2023). "Rotation in vacuum and scalar background." arXiv:2301.08118
17. Mychelkin, E.G., Makukov, M.A. & Suliyeva, G. (2024). "On the weak and strong field effects in antiscalar background." GRG. arXiv:2403.11610
18. Mychelkin, E.G., Suliyeva, G. & Makukov, M.A. (2025/2026). "The exponential metric: traversable wormhole and possible identification of scalar background." arXiv:2510.15391
19. Boonserm, P. et al. (2018). "The exponential metric represents a traversable wormhole." PRD 98, 084048. arXiv:1805.03781
20. Nandkishore, A.K. et al. (2024). "A new class of traversable wormhole metrics." EPJC 84, 1063
21. Nath, P.P. & Sarma, D. (2025). "Exponential wormhole, QNMs, and logarithmic f(R) gravity." Ann. Phys. 479, 170067
22. Ernazarov, K.K. (2025). "Yilmaz-Rosen and JNW in sEGB." arXiv:2510.21625
23. Ernazarov, K.K. (2025). "Asymptotically Schwarzschild-like metric solutions." arXiv:2511.13471
24. Manna, T. et al. (2023). "Strong lensing in exponential wormhole spacetimes." New Astron. Rev. 97, 101695
25. Guo, X. & Yuan, F. (2025). "Quantum Effective Dynamics of Papapetrou Spacetime." EPJC. arXiv:2506.08821

### Entropic Gravity
26. Jacobson, T. (1995). "Thermodynamics of Spacetime: The Einstein Equation of State." PRL 75, 1260
27. Verlinde, E. (2011). "On the Origin of Gravity and the Laws of Newton." JHEP 1104, 029. arXiv:1001.0785
28. Verlinde, E. (2016). "Emergent Gravity and the Dark Universe." SciPost Phys. 2, 016. arXiv:1611.02269
29. Dorau, P. & Much, A. (2025). "From Quantum Relative Entropy to the Semiclassical Einstein Equations." PRL. arXiv:2510.24491
30. Rostami, A., Rezazadeh, K. & Rostampour, M. (2025). "Relativistic MOND Theory from Modified Entropic Gravity." arXiv:2511.05632

### Yilmaz Theory Criticism
31. Cooperstock, F.I. & Vollick, D.N. (1996). "The Yilmaz challenge to general relativity"
32. Misner, C.W. (1999). "Yilmaz Cancels Newton" (disputed by Alley & Yilmaz 1999)
33. Robertson, S.L. (various). Astrophysical applications of Yilmaz gravity

### JNW / Scalar Solutions
34. Janis, A.I., Newman, E.T. & Winicour, J. (1968). PRL 20, 878
35. Xanthopoulos, B.C. & Zannias, T. (1989). J. Math. Phys. 30, 2656

---

## Appendix A: GfE vs Exponential Metric — Side-by-Side Comparison

### A.1 Mathematical Structure

| Feature | GfE (Bianconi) | Exponential Metric |
|---------|----------------|-------------------|
| Action | S = ∫ Tr(g̃ ln G̃⁻¹) √(-g) d⁴x | S = ∫ (R − ½∂φ² + ...) √(-g) d⁴x |
| Field content | g_μν + G-field | g_μν + phantom scalar φ |
| Fundamental principle | Quantum relative entropy | Scalar charge = mass |
| Vacuum equations | G-field survives in vacuum | Phantom field φ = M/r |
| Energy conditions | Satisfied (standard matter) | ALL violated (phantom) |
| Cosmological constant | Emergent (positive) | Not addressed |
| Dark matter | G-field candidate | Not addressed |

### A.2 Predictions Comparison

| Observable | GfE perturbative (β small) | Exponential metric |
|-----------|---------------------------|-------------------|
| g₀₀(r) | 1 − r_S/r − βr_S²/(48r⁴) | exp(−r_S/r) |
| Horizon? | Yes (shifted) | No (wormhole) |
| Leading correction | r⁻⁴ | r⁻² |
| Perihelion precession | ~GR + tiny corrections | ~GR (indistinguishable) |
| Shadow | ~Schwarzschild + corrections | Larger than Schwarzschild |
| QNMs | BH + corrections | Wormhole echoes |
| Hawking radiation | T ∝ M⁻¹/² | None |

---

## Appendix B: Open Questions and Future Work

1. **Solve the full nonlinear GfE equations** for static spherically symmetric vacuum (not just perturbative)
2. **Compute the G-field profile** for a galaxy → test dark matter claim
3. **Compare GfE and exponential metric predictions** for EHT observations (shadow, photon ring)
4. **Determine whether GfE at large β** could produce a horizonless object (wormhole-like)
5. **Connect GfE to quantum channel** approach (is the G-field related to the channel's decoherence parameter?)
6. **Test whether GfE's Λ_G** gives the correct cosmological constant value
7. **Compare GfE inflation** (Thattarampilly & Zheng) with Starobinsky inflation
