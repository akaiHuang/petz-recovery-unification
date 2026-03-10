# τ Framework in Dynamic Spacetimes
## Date: 2026-03-10
## Author: Research synthesis for Sheng-Kai Huang
## Status: Deep research complete

---

## 1. Problem Statement

### 1.1 The Static Limitation

Paper 2 establishes the triple identification for static spacetimes:

```
c_eff/c = exp(-Σ_grav/2) = F_bound
```

with the covariant formulation (Paper 2, Eq. 18):

```
F_bound(A → B) = |ξ|_A / |ξ|_B = exp(-ΔΣ_grav/2)
```

where ξ^a is a **timelike Killing vector field**. This requires the spacetime to be stationary—a severe restriction.

### 1.2 Why Dynamic Extension Matters

The most physically important scenarios are dynamic:

| Scenario | Why it matters | Killing vector? |
|----------|---------------|-----------------|
| Gravitational collapse (Vaidya) | Black hole formation | NO |
| FRW cosmology | Expanding universe, dark energy | NO |
| Hawking radiation / evaporation | Information paradox resolution | NO (backreaction) |
| Binary merger (GW) | LIGO observations | NO |
| Primordial inflation | Cosmic perturbation generation | NO (approximate de Sitter) |

If τ only works in static spacetimes, it describes **equilibrium** situations but misses all the **non-equilibrium** physics where entropy production is most interesting.

### 1.3 What We Need

A covariant definition τ_cov that:
1. **Reduces** to τ = 1 − exp(−r_s/(2r)) in the static Schwarzschild/exponential limit
2. **Is well-defined** in arbitrary dynamic spacetimes (no Killing vector required)
3. **Has physical content**: gives new predictions beyond static case
4. **Is mathematically rigorous**: based on established frameworks (AQFT, quasi-local horizons, etc.)

---

## 2. Literature Survey (2024–2026)

### 2.1 Kodama Vector and Dynamic Horizons

| Paper | Key Result | Relevance |
|-------|-----------|-----------|
| Kodama (1980) | Divergence-free vector K^a in spherically symmetric spacetimes | Foundation |
| Hayward (1994, 1998) | Unified first law at trapping horizons using Kodama flow | Thermodynamics without Killing vectors |
| Abreu & Visser (2010, PRD 82, 044027; arXiv:1004.1456) | Geometrically preferred foliations via Kodama vector; surface gravity throughout spacetime | τ_K definition possible at all r |
| Dorau & Verch (2024, arXiv:2402.18993) | Kodama-like vector fields in axisymmetric (Kerr-Vaidya, Kerr-Vaidya-de Sitter) spacetimes | **Extension beyond spherical symmetry** |
| arXiv:2402.16484 (2024) | Geometrical origin of Kodama vector as obstruction to foliation by marginally trapped surfaces | Conceptual clarity |
| Ashtekar & Krishnan (2025, arXiv:2502.11825) | Comprehensive review: quasi-local horizons replace event horizons | **State of the art** on dynamical horizons |
| arXiv:2506.04469 (2025) | Generalized Misner-Sharp energy in f(R,G) gravity via Kodama vector | Extension to modified gravity |

### 2.2 Modular Flow / Algebraic QFT

| Paper | Key Result | Relevance |
|-------|-----------|-----------|
| Witten (2018, arXiv:1803.04993) | Entanglement in QFT: modular theory primer | Foundation |
| Witten (2022, arXiv:2112.12828) | Gravity and the crossed product: Type III → Type II | Entropy becomes well-defined |
| Chandrasekaran, Longo, Penington, Witten (2023, arXiv:2206.10780) | Algebra for de Sitter space via crossed product | **Dynamic spacetime algebra** |
| Dorau & Much (2025, PRL; arXiv:2510.24491) | QRE → semiclassical Einstein equations | Σ_grav derivation |
| Kudler-Flam et al. (2025, arXiv:2405.00114; JHEP 07, 2025, 146) | **Gravitational entropy is observer-dependent** | τ must be observer-dependent |
| arXiv:2412.15502 (2025; JHEP 07, 2025, 063) | Crossed products and QRFs: observer-dependence of gravitational entropy | Formalization of observer-dependence |
| Chandrasekaran (2025, JHEP 07, 2025, 192) | **Generalized second law beyond semiclassical regime** | GSL without Killing vectors |
| Cirafici (2024, arXiv:2408.04219; JHEP 11, 2024, 089) | **Fluctuation theorems + quantum channels + gravitational algebras in de Sitter** | Crooks theorem in dynamic de Sitter |
| arXiv:2503.14107 (2025) | Von Neumann algebraic approach to QFT on curved spacetime | Mathematical foundation |
| arXiv:2511.00622 (2025) | Algebra for covariant observers in de Sitter space | FRW-relevant |

### 2.3 Raychaudhuri + Entropy + Quantum Information

| Paper | Key Result | Relevance |
|-------|-----------|-----------|
| Wall (2012, PRD; arXiv:1105.3445) | GSL from QRE monotonicity | Foundation: DPI → area increase |
| arXiv:2412.07107 (2024) | Gravitational focusing + horizon entropy for higher-spin fields | Generalized Wall entropy |
| arXiv:2509.00628 (2025) | Generalised focusing theorem in diffeomorphism-invariant theories | **Wall entropy from generalised Raychaudhuri** |
| arXiv:2601.18860 (2026) | Towards proof of Improved QNEC | QNEC from modular Hamiltonian + Raychaudhuri |
| Moreira & Celeri (2024, arXiv:2407.21186) | Entropy production from spacetime fluctuations (graviton bath) | Microscopic Σ_grav |
| Basso, Maziero & Celeri (2025, PRL 134, 050406; arXiv:2405.03902) | **Quantum Crooks theorem in curved spacetime** | Σ from Riemann tensor, observer-dependent |

### 2.4 Quasi-Local Approaches

| Paper | Key Result | Relevance |
|-------|-----------|-----------|
| Brown & York (1993) | Quasi-local energy from boundary stress tensor | Foundation |
| Wang & Yau (2009) | Gauge-independent quasi-local mass | Improvement over Brown-York |
| arXiv:2406.10751 (2024) | Strong-field behavior of Wang-Yau quasi-local energy near apparent horizons | Behavior near τ → 1 |
| arXiv:2402.19310 (2024) | Remarks on Wang-Yau quasi-local mass: connection to Brown-York | Unification |
| Communications in Mathematical Physics (2024) | New gauge-independent quasi-local mass | Modern formulation |

### 2.5 FRW Cosmology + Quantum Information

| Paper | Key Result | Relevance |
|-------|-----------|-----------|
| Padmanabhan (2017, PLB 773; arXiv:1703.06144) | CosMIn: finite cosmic information requires Λ > 0 | τ_cosmological motivation |
| Capozziello et al. (2024, arXiv:2406.19274, EPJC) | Bogoliubov transformation = Gaussian bosonic thermal noise channel | Cosmological expansion as quantum channel |
| Barman et al. (2026, QIP; arXiv:2601.20860) | Quantum teleportation fidelity in FRW universe | F degradation from expansion |
| Bianconi (2025, arXiv:2510.22545) | GfE thermodynamics: FRW solutions with k-temperatures and dynamical Λ | GQRE cosmology |
| Bianconi (2025, PRD 111, 066001; arXiv:2408.14391) | Gravity from entropy: metric = density matrix | Unified equation framework |
| Cai, Wang & Zhao (2024, arXiv:2407.09912) | GR fluctuation theorems, covariant | Σ in general spacetimes |

---

## 3. Route Analysis

### 3a. Kodama Vector Route

#### 3a.1 Definition and Properties

The Kodama vector K^a is defined in any spherically symmetric spacetime (not necessarily static):

```
K^a = ε^(ab) ∇_b R
```

where R is the areal radius and ε^(ab) is the volume form on the 2D Lorentzian orbit space. In coordinates ds² = h_{AB} dx^A dx^B + R²(x^A) dΩ², the Kodama vector is:

```
K^A = (1/√(-det h)) ε^{AB} ∂_B R
```

**Key properties (Hayward 1994, Abreu-Visser 2010)**:
1. ∇_a (G^a_b K^b) = 0 — divergence-free when contracted with Einstein tensor
2. Conserved charge = Misner-Sharp mass E_MS = R(1 − h^{AB}∂_A R ∂_B R)/2
3. K^a is timelike outside the apparent horizon, null on it, spacelike inside
4. In static case: K^a = ξ^a (Killing vector), so it reduces correctly

**For FRW metric** ds² = −dt² + a(t)²(dr²/(1−kr²) + r²dΩ²):
```
K^a = (1, −Hr/√(1−kr²), 0, 0)   in comoving coordinates
```

where H = ȧ/a is the Hubble parameter.

#### 3a.2 Proposed τ_K Definition

Define the Kodama-based recovery fidelity between points A and B:

```
F_K(A → B) ≡ |K|_A / |K|_B
```

where |K| = √(−g_{ab} K^a K^b) for timelike K.

**Static limit check**: For Schwarzschild, K^a = ξ^a (the timelike Killing vector), so |K| = √(−g_{00}) and:
```
F_K = √(−g_{00}(A)) / √(−g_{00}(B)) = exp(−ΔΣ_grav/2)   ✓
```

This matches Paper 2, Eq. 18 exactly.

**FRW cosmology**: At the apparent horizon r_AH = 1/H (flat case):
- |K|² = −1 + H²R² where R = a(t)r is the areal radius
- r = 0 (observer): |K|² = −1, F_K = 1, τ = 0
- r = r_AH: |K|² = 0, F_K = 0, τ = 1 — consistent with horizon = complete information loss

**Kodama entropy production**:
```
Σ_K = −2 ln(F_K) = −2 ln(|K|_A / |K|_B)
```

For the exponential metric, this gives Σ_K = r_s/r exactly (static limit).

#### 3a.3 Dynamical Predictions

**Gravitational collapse (Vaidya spacetime)**: ds² = −(1 − 2M(v)/r) dv² + 2dvdr + r²dΩ²

The Kodama vector exists and is timelike outside the apparent horizon r_AH(v) = 2M(v). As M(v) increases (infalling radiation), the apparent horizon expands, and F_K at a fixed radius r changes:

```
F_K(r, v) = √(1 − 2M(v)/r)   [for Vaidya, reduces to Schwarzschild form]
```

**New prediction**: During collapse, τ_K at a fixed r **increases monotonically** with M(v). The surface τ_K = 1 tracks the apparent horizon, NOT the event horizon.

**Hawking evaporation**: As M(v) decreases, τ_K at fixed r **decreases** — information recovery improves as the black hole evaporates. This is physically correct: the Petz recovery fidelity should improve as the mass decreases.

#### 3a.4 Thermodynamic Support

Hayward's unified first law at the trapping horizon:

```
dE = A/(8πG) κ_K dR_AH + W dV
```

where κ_K is the **Kodama surface gravity**:

```
κ_K = (1/2√(-h)) ∂_A(√(-h) h^{AB} ∂_B R)
```

Projecting along the Kodama vector gives the second Friedmann equation (Hayward 1998):

```
Ḣ − k/a² = −4πG(ρ + p)
```

This means **Friedmann equations ARE the first law of thermodynamics at the apparent horizon** — exactly what the τ framework predicts: gravitational dynamics = information dynamics.

#### 3a.5 Limitation: Spherical Symmetry

The original Kodama vector requires spherical symmetry. However:

**Dorau & Verch (2024, arXiv:2402.18993)** have constructed Kodama-like vectors for **axisymmetric** spacetimes (Kerr-Vaidya, Kerr-Vaidya-de Sitter):
- Timelike vector field with covariantly conserved current
- Associated Brown-York mass is well-defined
- Reduces to Kodama in spherical limit and Killing in stationary limit

This partially addresses the limitation. However, for **generic** spacetimes without any symmetry (e.g., binary mergers), no Kodama-like vector exists.

#### 3a.6 Assessment

| Criterion | Score |
|-----------|-------|
| Static limit recovery | EXACT |
| Mathematical rigor | HIGH (well-established geometric object) |
| Physical predictions | YES (τ tracks apparent horizon, not event horizon) |
| Generality | MEDIUM (spherically symmetric + axisymmetric extensions) |
| Connection to QI | MODERATE (thermodynamic, not directly quantum channel) |

**Maturity: PROMISING** — ready for concrete calculations, but limited to symmetric spacetimes.

---

### 3b. Modular Flow (Algebraic QFT)

#### 3b.1 The Framework

Tomita-Takesaki theory provides time evolution from the state-algebra pair (ω, M) alone, without reference to any background symmetry:

1. Given a von Neumann algebra M and a cyclic separating state |Ω⟩:
2. Define S: a|Ω⟩ ↦ a†|Ω⟩ (the Tomita operator)
3. Polar decomposition: S = JΔ^{1/2} where Δ is the **modular operator**
4. The **modular flow** σ_t(a) = Δ^{it} a Δ^{-it} defines intrinsic "time" evolution

**No background geometry needed.** The modular flow is determined entirely by the algebra and the state.

#### 3b.2 Connection to Gravity

In specific cases, modular flow has geometric meaning:

| Setting | Modular flow = | Reference |
|---------|---------------|-----------|
| Rindler wedge (Minkowski vacuum) | Boost (Unruh effect) | Bisognano-Wichmann (1976) |
| Diamond (Minkowski vacuum) | Conformal Killing flow | Hislop-Longo (1982) |
| Schwarzschild (Hartle-Hawking state) | Killing flow (Hawking temperature) | Sewell (1982) |
| de Sitter (Bunch-Davies state) | Killing flow (Gibbons-Hawking temperature) | Chandrasekaran et al. (2023) |

**For dynamic spacetimes**: The modular flow **still exists** (Tomita-Takesaki is a theorem about operator algebras, not about spacetime), but it may **not** have a simple geometric interpretation.

#### 3b.3 Modular τ Definition

Define the modular entropy production:

```
Σ_mod = D(ω_1 || ω_2) − D(ω_1 ∘ α || ω_2 ∘ α)
```

where:
- ω_1, ω_2 are states on the algebra M
- α: M → N is the inclusion/restriction map (the quantum channel)
- D is the Araki relative entropy: D(ω_1 || ω_2) = −⟨Ω_1 | ln Δ_{ω_2, ω_1} | Ω_1⟩

Then:
```
τ_mod = 1 − F_mod,   F_mod ≥ exp(−Σ_mod/2)
```

This is the **JRSWW bound generalized to von Neumann algebras** (it holds for general von Neumann algebras, not just matrix algebras; see Junge et al. 2018).

#### 3b.4 Static Limit Recovery

In Schwarzschild with Hartle-Hawking state:
- Modular Hamiltonian K = (2π/κ) H at the bifurcate horizon
- Araki relative entropy S^rel(ω_0 || ω_φ) = −2π ∫_H U ⟨:T_{UU}:⟩ dU dvol_S
- Fractional QRE loss = r_s/r (Dorau-Much derivation)
- τ_mod = 1 − exp(−r_s/(2r)) in the saturation ansatz

**Static limit: EXACT** ✓

#### 3b.5 Dynamic Extension: Crossed Product Algebras

The breakthrough of Witten (2022) and Chandrasekaran et al. (2023):

1. QFT observables in a subregion form a Type III von Neumann algebra → **entropy is undefined** (all states have infinite entropy)
2. Including the observer (gravitational degrees of freedom via crossed product) promotes it to **Type II** → entropy is well-defined
3. The resulting entropy = generalized entropy S_gen = A/(4G) + S_matter in the semiclassical limit

**For dynamic spacetimes**: The crossed product construction does **not** require Killing vectors. It requires:
- A von Neumann algebra M (from QFT in the subregion)
- A modular automorphism group (from Tomita-Takesaki)
- An observer / clock degree of freedom

**Key recent result** (Kudler-Flam et al. 2025, arXiv:2405.00114):
> "Gravitational entropy is observer-dependent"
> Different quantum reference frames (observers) give different crossed product algebras → different entropies → **different τ**.

This is **consistent** with the τ framework: τ was always defined relative to a channel N (which depends on the observer).

#### 3b.6 Application to de Sitter

Cirafici (2024, arXiv:2408.04219) established **fluctuation theorems in gravitational algebras** for de Sitter space:

- Two-point measurement scheme → dynamical fluctuations of observables
- Quantum channels represented by subfactors (Jones theory)
- Fluctuation theorem: P(Σ)/P(−Σ) = exp(Σ) — the **Crooks relation holds** in de Sitter crossed product algebra

This means the Crooks → Petz → τ chain extends to de Sitter, which is a dynamical spacetime.

**Chandrasekaran et al. (2023, arXiv:2206.10780)**: The Type II algebra for static patch of de Sitter has well-defined entropy. The modular flow corresponds to the Killing flow of the static de Sitter metric.

**Beyond the static patch**: For a general FRW spacetime (not exactly de Sitter), the modular flow is not geometric. But it still exists as an algebraic object, and Σ_mod is still well-defined.

#### 3b.7 The Generalized Second Law

Chandrasekaran (2025, JHEP 07, 2025, 192): **Generalized second law beyond the semiclassical regime**.

The GSL ΔS_gen ≥ 0 can be proven using:
1. Crossed product algebras (no Killing vector needed)
2. Monotonicity of relative entropy (DPI)
3. Wall's argument: QRE decrease under restriction = entropy production

This is **exactly** the mathematical structure of the τ framework: DPI → Σ ≥ 0 → τ ≥ 0.

#### 3b.8 Assessment

| Criterion | Score |
|-----------|-------|
| Static limit recovery | EXACT (Bisognano-Wichmann → Killing flow) |
| Mathematical rigor | HIGHEST (operator algebras, theorem-level results) |
| Physical predictions | YES (observer-dependent τ, GSL without Killing vectors) |
| Generality | HIGHEST (works for any spacetime, any state, any algebra) |
| Connection to QI | DIRECT (Araki relative entropy = Umegaki generalization) |
| Computability | LOW (modular Hamiltonian unknown for most dynamic spacetimes) |

**Maturity: PROMISING (theoretical framework) / SPECULATIVE (concrete calculations)**

The modular approach is the **most mathematically rigorous** and **most general**, but the **hardest to compute with** in practice.

---

### 3c. Raychaudhuri Decomposition

#### 3c.1 The Raychaudhuri-τ Decomposition

Already identified in the MEMORY as a key result:

```
dΣ/ds = θ²/3 + σ² − ω² + R_{ab} u^a u^b
```

where:
- θ = expansion (volume change rate)
- σ = shear (shape change rate)
- ω = vorticity (rotation rate)
- R_{ab} u^a u^b = Ricci focusing (matter content via Einstein equations)
- s = affine parameter along the congruence

This is **already covariant** — no Killing vector required.

#### 3c.2 Physical Interpretation

Each term has an information-theoretic meaning:

| Term | Sign | Physical meaning | Information effect |
|------|------|-----------------|-------------------|
| θ²/3 | ≥ 0 | Expansion (stretching) | Dilution of information |
| σ² | ≥ 0 | Shear (distortion) | Scrambling |
| −ω² | ≤ 0 | Vorticity (rotation) | **Protection** (frame dragging) |
| R_{ab}u^au^b | Either | Matter/energy content | Focusing or defocusing |

For null energy condition (NEC): R_{ab}k^ak^b ≥ 0 → matter always focuses.
Vorticity **reduces** Σ — spinning objects protect information (consistent with Kerr having a smaller τ than Schwarzschild at the same mass).

#### 3c.3 Integration: Local vs Global

The Raychaudhuri equation gives a **local** (differential) Σ_Raych. To get the integrated Σ:

```
Σ_Raych(A → B) = ∫_A^B [θ²/3 + σ² − ω² + R_{ab}u^au^b] ds
```

**Static Schwarzschild vacuum**: R_{ab} = 0, and for a static congruence θ = σ = ω = 0, so Σ_Raych = 0.

This is **NOT** the same as Σ_grav = r_s/r! The two types of Σ measure different things:

| Quantity | What it measures | Static Schwarzschild |
|----------|-----------------|---------------------|
| Σ_grav = −ln(−g_{00}) | Redshift (integrated potential) | r_s/r |
| Σ_Raych | Geodesic focusing (local curvature) | 0 (vacuum, static) |

#### 3c.4 Resolution: Complementary Σ's

The two Σ's are **complementary**, not contradictory:

1. **Σ_grav** (redshift): measures the **cumulative cost** of maintaining a static observer at radius r. It's the entropy production of the **static channel** (stationary observer → asymptotic observer).

2. **Σ_Raych** (focusing): measures the **local rate** of information loss along a **freely falling congruence**. It's the entropy production of the **dynamical channel** (null/timelike geodesics).

In a dynamic spacetime (e.g., collapse), **both** are nonzero and contribute differently:
- Σ_grav increases as the mass grows (M(v) increases)
- Σ_Raych accumulates along the geodesic from the Ricci focusing

**Key insight**: In a fully dynamic scenario, Σ_total should include **both** contributions. The Raychaudhuri part captures the genuine non-equilibrium dynamics that Σ_grav alone misses.

#### 3c.5 Connection to Wall's Theorem and GSL

Wall (2012, arXiv:1105.3445) proved:

```
ΔS_gen = ΔA/(4G) + ΔS_matter ≥ 0
```

using monotonicity of relative entropy. The area increase comes from the **Raychaudhuri focusing**:

```
dA/dλ = ∫ θ dA_⊥
```

and θ satisfies the Raychaudhuri equation. So the GSL is essentially:

```
Σ_Wall = D(ρ‖σ) − D(N(ρ)‖N(σ)) = ΔA/(4G) + ΔS_matter ≥ 0
```

This is exactly the JRSWW entropy production! The connection:

```
Wall's GSL ←→ DPI for the restriction channel ←→ JRSWW bound ←→ τ framework
```

#### 3c.6 Recent Developments

**Gravitational focusing for higher-spin fields** (arXiv:2412.07107, 2024): Wall entropy extracted from generalised Raychaudhuri-type equations satisfies the first and second laws. Higher-spin fields introduce indefinite terms that obstruct the focusing theorem, but a "higher-spin focusing condition" restores it.

**Improved QNEC** (arXiv:2601.18860, 2026): Using modular Hamiltonian and relative entropy:

```
⟨T_{kk}⟩ ≥ (ℏ/2π) S''_{EE}
```

where S''_{EE} is the second null derivative of entanglement entropy. This is a **quantum version** of the Raychaudhuri focusing, directly relating energy-momentum (gravity) to entanglement entropy changes (quantum information).

#### 3c.7 Assessment

| Criterion | Score |
|-----------|-------|
| Static limit recovery | PARTIAL (gives Σ_Raych = 0 for vacuum, not Σ_grav = r_s/r) |
| Mathematical rigor | HIGH (covariant, well-established) |
| Physical predictions | YES (focusing, singularity theorems, GSL) |
| Generality | HIGH (any congruence in any spacetime) |
| Connection to QI | HIGH (via Wall → DPI → JRSWW) |
| Computability | HIGH (can integrate along any geodesic congruence) |

**Maturity: READY (as a component) / INCOMPLETE (as standalone τ definition)**

The Raychaudhuri approach is computationally tractable and covariant, but captures **different physics** (focusing) than the redshift-based Σ_grav. It should be a **component** of the full dynamical τ, not the whole story.

---

### 3d. Quasi-Local Approaches

#### 3d.1 Brown-York Quasi-Local Energy

The Brown-York quasi-local energy for a 2-surface S bounding a spacelike region:

```
E_BY = (1/8πG) ∮_S (k − k_0) dA
```

where k is the trace of extrinsic curvature and k_0 is the reference (flat space) subtraction.

**Proposed Σ_BY**:
```
Σ_BY = 2 |Φ_BY|/c² = 2 G E_BY / (R c²)
```

where R is the areal radius of S.

**Static limit**: For Schwarzschild, E_BY = M√(1 − r_s/r) [at areal radius r], so:
```
Σ_BY = 2GM√(1 − r_s/r) / (rc²) = (r_s/r)√(1 − r_s/r) ≠ r_s/r
```

This does NOT reproduce the target Σ_grav = r_s/r. The Brown-York energy includes a redshift factor that modifies the result.

#### 3d.2 Misner-Sharp Energy (Kodama-Associated)

The Misner-Sharp mass E_MS = R(1 − h^{AB}∂_A R ∂_B R)/2 is the conserved charge associated with the Kodama vector. For Schwarzschild: E_MS = M (exact, at all r).

**Proposed Σ_MS**:
```
Σ_MS = 2G E_MS / (Rc²) = r_s/r    ✓
```

This DOES reproduce the target! And it works in dynamic spherically symmetric spacetimes because E_MS is well-defined via the Kodama vector even without a Killing vector.

**Key insight**: The Misner-Sharp mass, not the Brown-York energy, is the correct quasi-local quantity for defining Σ in the τ framework.

#### 3d.3 Wang-Yau Quasi-Local Mass

Wang-Yau quasi-local mass (2009) is gauge-independent and embeds the 2-surface into **Minkowski space** (not R³ as Brown-York):

Recent results (arXiv:2406.10751, 2024): Near an apparent horizon:
- If the horizon cannot be isometrically embedded in R³: E_WY → ∞
- If it can: E_WY has a finite limit

This divergence behavior near horizons suggests that Wang-Yau mass could give τ → 1 at the apparent horizon, consistent with the framework.

However, no explicit connection between Wang-Yau mass and quantum channels/Petz recovery has been established.

#### 3d.4 Assessment

| Criterion | Score |
|-----------|-------|
| Static limit recovery | YES (via Misner-Sharp), NO (via Brown-York) |
| Mathematical rigor | HIGH (well-defined geometric quantities) |
| Physical predictions | LIMITED (energy, not directly entropy production) |
| Generality | MEDIUM (Misner-Sharp: spherically symmetric; Wang-Yau: general) |
| Connection to QI | WEAK (no direct quantum channel interpretation) |

**Maturity: SPECULATIVE** — The quasi-local energy approach provides useful geometric quantities but lacks a direct connection to quantum information theory. The Misner-Sharp mass is useful primarily because it is the Kodama charge — so this route essentially reduces to the Kodama vector route (3a).

---

### 3e. FRW Cosmology

#### 3e.1 Expansion as a Quantum Channel

**Capozziello et al. (2024, arXiv:2406.19274)**:
- Bogoliubov transformation of field modes in FRW → Gaussian bosonic thermal noise channel
- Channel parameters determined by particle production rate
- Modified gravity theories preserve more information (lower τ)

**Barman et al. (2026, arXiv:2601.20860)**:
- Teleportation fidelity F in FRW explicitly computed
- F decreases with expansion rate (H) and mode frequency
- Bogoliubov coefficients quantify fidelity degradation

**Key result**: Cosmological expansion IS a quantum channel N_cosmo, and fidelity degradation from Bogoliubov mixing provides a concrete, computable τ_cosmo.

#### 3e.2 Σ from Curvature (Basso-Celeri)

Basso, Maziero & Celeri (PRL 134, 050406, 2025):

Curved-spacetime Crooks relation:
```
P_fwd(W) / P_rev(−W) = exp[β(W − ΔF)]
```

Hamiltonian includes curvature coupling:
```
H_cm(τ) = m + p²/(2m) + m a_i(τ)x^i + (m/2) R_{τiτj}(τ) x^i x^j
```

The last term (tidal coupling to Riemann tensor) drives the system out of equilibrium.

**In de Sitter**: For QHO with frequency ω₀, transition probability ~ (H/ω₀)⁴. Per Hubble time, Σ ~ O(1) for modes at ω ~ H.

**Crucial property**: **Entropy production is observer-dependent**:
- Different worldlines → different R_{τiτj} → different Σ
- This matches the τ framework: τ depends on which channel N you observe through

#### 3e.3 Kodama Vector in FRW

As derived in Section 3a:
```
F_K(r) = √(1 − H²R²)  for flat FRW (R = areal radius)
```

- At observer (R = 0): F_K = 1, τ = 0
- At apparent horizon (R = 1/H): F_K = 0, τ = 1
- Smooth interpolation between

**Cosmological Σ_K**:
```
Σ_K = −2 ln(F_K) = −ln(1 − H²R²)
```

For R << 1/H (well inside the Hubble volume): Σ_K ≈ H²R², which is small.
At R = 1/H: Σ_K → ∞, τ → 1.

#### 3e.4 Three Independent Routes to Σ ~ O(1) at Hubble Scale

Already identified in paper3_cosmological_tau.md:

**Route A (Kodama)**: |K|² → 0 at apparent horizon → τ → 1

**Route B (Cumulative)**: Σ_cum ~ (4πGρ₀/3)R²/c² = (R/R_H)² Ω_m ~ 0.3 at R ~ R_H

**Route C (Bianconi GfE)**: GQRE per degree of freedom per Hubble time ~ O(1)

#### 3e.5 Bianconi's GfE Cosmology (arXiv:2510.22545)

FRW solutions in GfE framework:
- k-temperatures: θ_k ~ ω_k H² (5 distinct temperatures for scalar/vector/bivector modes)
- Entropy per unit volume: δs/δv = −ln(1 − ω_k H²), decreasing (as befits relative entropy)
- Total entropy S ∝ ω̄₁ t^{6/n−1}, increasing in time (volume growth dominates)
- **Dynamical cosmological constant** Λ^G = H/(2β) Σ_k z_k [G_k − 1 − ln G_k]

**Connection to τ**: In the GfE framework, the GQRE provides a natural Σ_cosmo that:
1. Equals zero in flat Minkowski (correct)
2. Increases with H (correct)
3. Has well-defined thermodynamics (first law, k-temperatures)

But recall the fundamental issue (Section 3 of paper4_bianconi_gfe_research.md): GQRE depends on R_s/r³ (curvature), while Σ_grav depends on R_s/r (potential). They are related by integration but NOT identical.

#### 3e.6 Crooks Theorem in Expanding Universe

**Status**: Yes, the Crooks theorem holds in expanding universes:
1. Basso-Celeri (PRL 2025): Fully general-relativistic quantum Crooks theorem, explicitly verified for FRW
2. Cirafici (JHEP 2024): Fluctuation theorems established in gravitational algebras of de Sitter space
3. Cai-Wang-Zhao (arXiv:2407.09912, 2024): Covariant GR fluctuation theorems

The Crooks → Petz → τ chain is valid in cosmological settings.

#### 3e.7 Assessment

| Criterion | Score |
|-----------|-------|
| Static limit recovery | N/A (cosmological regime, different physical setup) |
| Mathematical rigor | MEDIUM (Bogoliubov transformations well-defined; Σ computation ongoing) |
| Physical predictions | YES (F degradation vs H, observer-dependent τ) |
| Generality | HIGH (any FRW cosmology) |
| Connection to QI | STRONG (explicit quantum channels from Bogoliubov mixing) |

**Maturity: PROMISING** — Concrete calculations exist (teleportation fidelity, Crooks in FRW). The Kodama-based τ_cosmo has a clean definition. Main gap: connecting Bogoliubov-channel τ to Kodama-τ quantitatively.

---

## 4. Comparison Table

| Route | Maturity | Static Limit | Generality | QI Connection | Computability | Main Limitation |
|-------|----------|-------------|-----------|--------------|---------------|-----------------|
| **3a. Kodama** | PROMISING | EXACT | Spherically symmetric + axisymmetric | MODERATE | HIGH | No generic spacetimes |
| **3b. Modular flow** | PROMISING/SPECULATIVE | EXACT | ANY spacetime | DIRECT | LOW | Modular Hamiltonian unknown for most cases |
| **3c. Raychaudhuri** | READY (component) | PARTIAL (different Σ) | ANY spacetime | HIGH (via Wall/DPI) | HIGH | Measures focusing, not redshift |
| **3d. Quasi-local** | SPECULATIVE | YES (Misner-Sharp) | MEDIUM | WEAK | MEDIUM | Reduces to Kodama route |
| **3e. FRW cosmology** | PROMISING | N/A | FRW spacetimes | STRONG | MEDIUM | Specific to cosmology |

---

## 5. Static Limit Consistency Check

For each route, verify recovery of τ = 1 − exp(−r_s/(2r)) in static Schwarzschild/exponential:

### 5a. Kodama Vector
- K^a → ξ^a (Killing vector)
- |K| = √(−g_{00})
- F_K = √(−g_{00}(r))/√(−g_{00}(∞)) = exp(−Σ_grav/2)
- τ_K = 1 − exp(−r_s/(2r)) ✓ **EXACT**

### 5b. Modular Flow
- Modular Hamiltonian K = (2π/κ)H at bifurcate horizon
- Fractional QRE loss = r_s/r (Dorau-Much)
- Σ_mod = r_s/r in weak field
- τ_mod = 1 − exp(−r_s/(2r)) ✓ **EXACT** (with saturation ansatz)

### 5c. Raychaudhuri
- Static vacuum: R_{ab} = 0, θ = σ = ω = 0
- Σ_Raych = 0 ✗ **DOES NOT MATCH** (measures different physics)
- However, for null geodesics from r to ∞ through Schwarzschild:
  - The focusing integral gives the area theorem ΔA ≥ 0
  - This relates to generalized entropy, NOT to Σ_grav = r_s/r

### 5d. Quasi-Local (Misner-Sharp)
- E_MS = M for Schwarzschild
- Σ_MS = 2GM/(rc²) = r_s/r ✓ **EXACT**
- But this is just Σ = r_s/r by definition; the quasi-local framework adds no independent derivation

### 5e. FRW
- Static limit: H → 0, de Sitter becomes Minkowski
- Σ_K → 0, τ → 0 ✓ (trivial, but consistent)
- FRW regime does not directly recover the Schwarzschild limit

---

## 6. Recommendation

### 6.1 Primary Route: Kodama + Modular (Hybrid)

The optimal strategy combines the **computational tractability** of the Kodama vector with the **mathematical depth** of the modular flow approach:

**Layer 1 — Kodama (operational definition)**:
```
τ_K(A → B) = 1 − |K|_A / |K|_B
```
- Works for all spherically symmetric (and axisymmetric) spacetimes
- Immediately computable in Vaidya collapse, FRW cosmology, Kerr-Vaidya
- Gives concrete predictions (τ tracks apparent horizon, not event horizon)

**Layer 2 — Modular flow (foundational justification)**:
- Prove that the Kodama-based τ_K equals the modular Σ_mod in the static limit
- In the modular framework, the Kodama vector emerges as the geometric representative of the modular flow for spherically symmetric spacetimes
- The modular framework guarantees that the DPI/JRSWW structure extends to arbitrary spacetimes (even without Kodama)

**Layer 3 — Raychaudhuri (dynamical component)**:
- Include Σ_Raych as an additional contribution to total Σ in genuinely dynamic situations
- In static case, Σ_Raych = 0, so it does not interfere with the Kodama result
- In collapse/merger, Σ_Raych captures the focusing/shearing that Σ_K misses

### 6.2 Proposed Covariant τ Definition (for Paper 2 revision / Paper 4)

**For spherically symmetric spacetimes** (sufficient for Paper 2 + FRW):

```
Σ_cov = −2 ln(|K|_A / |K|_B)     [Kodama-based]
τ_cov = 1 − exp(−Σ_cov/2) = 1 − |K|_A / |K|_B
```

**For general spacetimes** (Paper 4 / future):

```
Σ_cov = D_Araki(ω_A || ω_0) − D_Araki(ω_B || ω_0)     [modular flow]
```

where:
- ω_A = state restricted to algebra at point/region A
- ω_0 = reference state (vacuum or thermal)
- D_Araki = Araki relative entropy on the von Neumann algebra

### 6.3 Concrete Next Steps (Priority Order)

1. **Immediate (for Paper 2 GRF Essay, deadline 3/31/2026)**:
   - Add one paragraph to Discussion/Open Questions: "The covariant formulation Eq. (18) requires a timelike Killing vector. For dynamical spacetimes, the Kodama vector provides a natural replacement in spherical symmetry (Kodama 1980, Dorau-Verch 2024). The extension to arbitrary spacetimes via modular flow (Witten 2022, Chandrasekaran et al. 2023) is an important open direction."
   - This is honest and forward-looking without overclaiming.

2. **Short-term (standalone calculation for Paper 2 revision)**:
   - Compute τ_K(r, v) for Vaidya collapse: explicit formula, plot τ vs r at different v
   - Show that τ_K = 1 surface = apparent horizon (NOT event horizon)
   - This is a new, concrete prediction that distinguishes the τ framework

3. **Medium-term (for Paper 4)**:
   - Compute τ_K for FRW with matter: Σ_K(R, t) = −ln(1 − H(t)²R²)
   - Compare with Bogoliubov-channel fidelity (Barman et al. 2026)
   - Establish quantitative connection between Kodama τ and channel τ

4. **Long-term (fundamental)**:
   - Prove that Kodama-based τ_K equals modular τ_mod in the static limit (should follow from Bisognano-Wichmann + Sewell)
   - Extend to axisymmetric case using Dorau-Verch Kodama-like vectors
   - Address binary merger: this requires a genuinely new approach (possibly numerical modular flow)

### 6.4 The Key Theoretical Prediction

**In the τ framework, the dynamical τ_K tracks the APPARENT HORIZON, not the EVENT HORIZON.**

This is physically significant because:
- The event horizon is teleological (requires knowledge of the future)
- The apparent horizon is quasi-local (determined by local geometry)
- Ashtekar & Krishnan (2025) argue that quasi-local horizons should replace event horizons
- τ_K = 1 at the apparent horizon is **consistent** with this modern view
- τ_K < 1 everywhere outside the apparent horizon (information always partially recoverable)

This prediction is **testable** in numerical relativity simulations of gravitational collapse.

---

## 7. Complete Reference List

### Foundational

1. **Kodama H** (1980). "Conserved quantities in the spherically symmetric gravitational field." Prog. Theor. Phys. 63, 1217.

2. **Hayward SA** (1994). "Quasilocal gravitational energy." Phys. Rev. D 49, 831.

3. **Hayward SA** (1998). "Unified first law of black-hole dynamics and relativistic thermodynamics." Class. Quant. Grav. 15, 3147. arXiv:gr-qc/9710089.

4. **Abreu G, Visser M** (2010). "Kodama time: Geometrically preferred foliations of spherically symmetric spacetimes." Phys. Rev. D 82, 044027. arXiv:1004.1456.

5. **Brown JD, York JW** (1993). "Quasilocal energy and conserved charges derived from the gravitational action." Phys. Rev. D 47, 1407. arXiv:gr-qc/9209012.

6. **Wang MT, Yau ST** (2009). "Quasilocal mass in general relativity." Phys. Rev. Lett. 102, 021101.

### Modular Flow / Algebraic QFT

7. **Bisognano JJ, Wichmann EH** (1976). "On the duality condition for quantum fields." J. Math. Phys. 17, 303.

8. **Witten E** (2018). "APS Medal for Exceptional Achievement in Research: Invited article on entanglement properties of quantum field theory." Rev. Mod. Phys. 90, 045003. arXiv:1803.04993.

9. **Witten E** (2022). "Gravity and the crossed product." JHEP 2022(10), 008. arXiv:2112.12828.

10. **Chandrasekaran V, Longo R, Penington G, Witten E** (2023). "An algebra of observables for de Sitter space." JHEP 2023(02), 082. arXiv:2206.10780.

11. **Dorau P, Much A** (2025). "From quantum relative entropy to the semiclassical Einstein equations." Phys. Rev. Lett. arXiv:2510.24491.

12. **Kudler-Flam J** et al. (2025). "Gravitational entropy is observer-dependent." JHEP 07, 146. arXiv:2405.00114.

13. **arXiv:2412.15502** (2025). "Crossed products and quantum reference frames: on the observer-dependence of gravitational entropy." JHEP 07, 063.

14. **Chandrasekaran V** (2025). "Generalised second law beyond the semiclassical regime." JHEP 07, 192.

15. **Cirafici M** (2024). "Fluctuation theorems, quantum channels and gravitational algebras." JHEP 11, 089. arXiv:2408.04219.

16. **arXiv:2503.14107** (2025). "A von Neumann algebraic approach to quantum theory on curved spacetime."

17. **arXiv:2511.00622** (2025). "An algebra for covariant observers in de Sitter space."

### Raychaudhuri + Focusing + Entropy

18. **Wall AC** (2012). "A proof of the generalized second law for rapidly-changing fields and arbitrary horizon slices." Phys. Rev. D 85, 104049. arXiv:1105.3445.

19. **arXiv:2412.07107** (2024). "Gravitational focusing and horizon entropy for higher-spin fields."

20. **arXiv:2509.00628** (2025). "Generalised focusing theorem and dynamical horizon entropy in diffeomorphism-invariant theories."

21. **arXiv:2601.18860** (2026). "Towards a proof of the Improved Quantum Null Energy Condition."

22. **Basso MLW, Maziero J, Celeri LC** (2025). "Quantum detailed fluctuation theorem in curved spacetimes." Phys. Rev. Lett. 134, 050406. arXiv:2405.03902.

23. **Moreira LC, Celeri LC** (2024). "Graviton-bath entropy production and decoherence." arXiv:2407.21186.

### Quasi-Local Horizons

24. **Ashtekar A, Krishnan B** (2025). "Quasi-local black hole horizons: recent advances." arXiv:2502.11825. (At press, Living Reviews in Relativity.)

25. **arXiv:2406.10751** (2024). "Strong field behavior of Wang-Yau quasi-local energy."

26. **arXiv:2402.19310** (2024). "Some remarks on Wang-Yau quasi-local mass."

27. **arXiv:2506.04469** (2025). "Generalized Misner-Sharp energy in f(R,G) gravity."

### Kodama Extensions

28. **Dorau P, Verch R** (2024). "Kodama-like vector fields in axisymmetric spacetimes." arXiv:2402.18993.

29. **arXiv:2402.16484** (2024). "Geometrical origin of the Kodama vector."

### FRW / Cosmological

30. **Padmanabhan T** (2017). "Cosmic information, the cosmological constant and the amplitude of primordial perturbations." Phys. Lett. B 773, 81. arXiv:1703.06144.

31. **Bianconi G** (2025a). "Gravity from entropy." Phys. Rev. D 111, 066001. arXiv:2408.14391.

32. **Bianconi G** (2025b). "The thermodynamics of the Gravity from Entropy Theory." arXiv:2510.22545.

33. **Capozziello S** et al. (2024). "Bogoliubov transformation as quantum channel in FRW." EPJC. arXiv:2406.19274.

34. **Barman S** et al. (2026). "Quantum teleportation in expanding FRW universe." QIP. arXiv:2601.20860.

35. **Cai RG, Wang B, Zhao L** (2024). "GR fluctuation theorems, covariant." arXiv:2407.09912.

### Information Theory / Petz Map

36. **Junge M, Renner R, Sutter D, Wilde MM, Winter A** (2018). "Universal recovery maps and approximate sufficiency of quantum relative entropy." Ann. Henri Poincaré 19, 2955.

37. **Hollands S, Longo R** (2021). "Relative entropy in curved spacetimes." arXiv:2107.06787.

38. **Casini H, Grillo AFG, Pontello D** (2019). "Relative entropy from Araki formula." Phys. Rev. D. arXiv:1903.00109.

---

## Appendix A: Key Equations Summary

### Static Framework (Paper 2)

```
Σ_grav = r_s/r = −ln(−g_{00})     [exponential metric]
F_bound = exp(−Σ_grav/2) = √(−g_{00})
τ = 1 − F_bound = 1 − exp(−r_s/(2r))

Covariant (static): F(A→B) = |ξ|_A / |ξ|_B     [ξ = Killing vector]
```

### Dynamic Extension (this document)

```
Kodama:     F_K(A→B) = |K|_A / |K|_B     [K = Kodama vector]
            Σ_K = −2 ln(|K|_A / |K|_B)
            τ_K = 1 − |K|_A / |K|_B

Modular:    Σ_mod = D_Araki(ω_A ‖ ω_0) − D_Araki(ω_B ‖ ω_0)
            F_mod ≥ exp(−Σ_mod/2)     [JRSWW on von Neumann algebras]

Raychaudhuri: dΣ_R/ds = θ²/3 + σ² − ω² + R_{ab}u^au^b

FRW (Kodama): Σ_K(R) = −ln(1 − H²R²)
              τ_K(R) = 1 − √(1 − H²R²)

Vaidya (Kodama): F_K(r,v) = √(1 − 2M(v)/r)
                 τ_K(r,v) = 1 − √(1 − 2M(v)/r)
```

### Hierarchy of Generality

```
Killing vector (static/stationary)
    ⊂ Kodama vector (spherically symmetric)
        ⊂ Kodama-like vector (axisymmetric; Dorau-Verch 2024)
            ⊂ Modular flow (arbitrary spacetime; Tomita-Takesaki)
```

Each level sacrifices computability for generality. The τ framework should be defined at the most general level (modular flow) but computed at the most specific level available for the problem at hand.

---

## Appendix B: Connection to Paper 2 Open Questions

Paper 2, Section VI.D lists as an open question:

> "Dynamic spacetimes. — The covariant formulation Eq. (18) requires a timelike Killing vector. For dynamical spacetimes (gravitational collapse, cosmological expansion), no such Killing vector exists, and the framework in its present form does not apply. Extension to dynamical settings — perhaps via approximate Killing vectors or the Kodama vector in spherical symmetry — is an important open direction."

This document provides the complete analysis of that open direction. The key findings:

1. **The Kodama vector IS the right replacement** for spherically symmetric dynamics (Vaidya collapse, FRW cosmology). It is not "approximate" — it is an exact geometric object.

2. **The modular flow approach provides the foundational justification** and extends to arbitrary spacetimes (Witten 2022, Chandrasekaran et al. 2023), though practical computation is difficult.

3. **The Raychaudhuri equation provides a complementary, computable quantity** that captures dynamical effects (focusing, shearing) not present in the redshift-based Σ.

4. **The key theoretical prediction** is that τ tracks the **apparent horizon** (quasi-local, no teleology) rather than the **event horizon** (global, teleological). This aligns with the modern quasi-local horizon program (Ashtekar-Krishnan 2025).

5. **The Crooks theorem holds in dynamic spacetimes** (Basso-Celeri 2025, Cirafici 2024), so the full Crooks → Petz → τ chain is valid beyond the static case.
