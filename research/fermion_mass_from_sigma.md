# Fermion Mass as Retrodiction Cost in Internal Space
## Paper 10 Research Note
## 2026-03-23
## Sheng-Kai Huang

---

## Executive Summary

We propose that **fermion mass is the cost of retrodiction failure in the internal (gauge) space**. In the Σ = 2 ln Q framework, the temporal asymmetry τ = 1 − F measures irreversibility. On spacetime, this produces gravity and dark matter (Papers 2–4). On the internal algebra A_F, the same structure should produce the fermion mass matrix.

The central formula:

```
m_f / M_ref = τ_f = 1 − F_f = 1 − exp(−Σ_f / 2)
```

where Σ_f = D(ρ_f ‖ ρ_vacuum) is the quantum relative entropy of fermion f in the internal space.

**Key finding**: Singh's parameter δ² = 3/8 in his J₃(O) mass ratio calculation numerically equals sin²θ_W at GUT scale. This is the **same algebraic origin** — the ratio of U(1) to total gauge degrees of freedom in A_F = C ⊕ H ⊕ M₃(C). If confirmed, this connects fermion masses to the electroweak mixing angle through Σ.

**Status labels**: Each step below is marked PROVEN, CONJECTURED, or SPECULATIVE.

---

## I. The Core Idea: Mass = Internal Retrodiction Cost

### I.1 Recap: Σ on spacetime

In the spacetime sector (Papers 1–4), we have:

- **PROVEN** (Paper 1): τ = 1 − F, where F ≥ exp(−Σ/2) is the Petz recovery fidelity
- **PROVEN** (Paper 2): On static backgrounds, Σ_grav = 2 ln Q = −ln(−g₀₀)
- **PROVEN** (Paper 9): The spectral Fisher metric g_QQ = (d−2)(d−3)/Q² matches Σ = 2 ln Q in d = 4

The pattern: **Σ_grav encodes how much information is lost by the gravitational channel. The metric is the consequence.**

### I.2 Extension to internal space

The Standard Model lives on a product geometry M × F (Connes):
- M = 4D spacetime with Dirac operator D_M (encodes gravity)
- F = internal finite space with Dirac operator D_F (encodes masses and mixings)

**CONJECTURE (Paper 10 core)**: The same Σ structure applies to F:

```
Σ_f = D(ρ_f ‖ ρ_0)_F
```

where:
- ρ_f = state of fermion f in the internal Hilbert space H_F
- ρ_0 = vacuum state (maximally mixed over accessible internal degrees of freedom)
- D(·‖·)_F = Umegaki QRE restricted to the internal algebra A_F

The fermion mass is then:

```
m_f = M_ref × (1 − exp(−Σ_f / 2))
```

where M_ref is a reference scale (see Section IV).

### I.3 Physical interpretation

**Heavy fermion (top quark, m_t ≈ 173 GeV)**:
- Σ_t ≫ 1
- τ_t ≈ 1
- The top quark state is "maximally irreversible" in the internal space
- Retrodicting the vacuum from the top state is nearly impossible
- F_t ≈ 0

**Light fermion (electron, m_e ≈ 0.511 MeV)**:
- Σ_e ≪ Σ_t
- τ_e ≪ 1
- The electron state is "nearly retrodictable" from the vacuum
- F_e close to 1

**Massless fermion (neutrino in the SM limit)**:
- Σ_ν = 0
- τ_ν = 0
- Perfect retrodiction = no mass = no asymmetry in the internal channel
- Neutrino mass ≠ 0 (DESI, oscillations) means Σ_ν > 0 but very small

This is satisfying: **mass hierarchy ↔ retrodiction hierarchy**.

---

## II. Σ_f from the Finite Dirac Operator

### II.1 The NCG setup [KNOWN]

The finite spectral triple of the Standard Model (Chamseddine–Connes 1996):

```
A_F = C ⊕ H ⊕ M₃(C)
H_F = C⁹⁶  (per generation × 3 generations = 96 × 3 = 288 total)
D_F = Yukawa coupling matrix (288 × 288)
```

For one generation, H_F decomposes as:
```
H_F^(1gen) = (ν_R, e_R, ν_L, e_L, u_R^{rgb}, d_R^{rgb}, u_L^{rgb}, d_L^{rgb}) + antiparticles
           = 16 Weyl spinors × 2 (particle/anti) × 3 (color for quarks) = 96
```

The Dirac operator D_F is block off-diagonal (it mixes L and R):
```
D_F = | 0      M† |
      | M      0  |
```

where M is the fermion mass matrix containing Yukawa couplings y_f:
```
M = diag(y_e, y_μ, y_τ, y_u, y_c, y_t, y_d, y_s, y_b) × v/√2
```
with v = 246 GeV the Higgs VEV.

### II.2 Defining Σ_f on A_F [CONJECTURED]

For a single fermion f, define the internal state as the Gibbs state of D_F restricted to the fermion's subspace:

```
ρ_f = exp(−β_f D_F²) / Tr exp(−β_f D_F²)|_{subspace f}
```

The vacuum reference state:
```
ρ_0 = I / dim(H_F)|_{subspace f}
```

Then:
```
Σ_f = D(ρ_f ‖ ρ_0) = Tr(ρ_f ln ρ_f) − Tr(ρ_f ln ρ_0)
     = −S(ρ_f) + ln(dim)
```

For a thermal state at inverse temperature β_f of a 2-level system (L-R mixing of one fermion species):
```
Σ_f = β_f m_f − ln(1 + e^{−β_f m_f}) + ln 2
```

In the low-temperature limit β_f m_f ≫ 1:
```
Σ_f ≈ β_f m_f − ln 2 + ln 2 = β_f m_f
```

In the high-temperature limit β_f m_f ≪ 1:
```
Σ_f ≈ (β_f m_f)²/4
```

**The question**: What determines β_f? There are two natural choices:

**Choice A**: Universal β = 1/(k_B T) at some reference scale (e.g., T = v = 246 GeV). Then:
```
Σ_f = m_f / v  (in the low-T limit)
```
This is just the Yukawa coupling: **Σ_f = y_f / √2**.

**Choice B**: Self-referential β_f = 1/m_f (each fermion at its own "internal temperature"). Then:
```
Σ_f = 1  for all fermions
```
This is universal and boring — it says all fermions are equally irreversible at their own scale.

**Choice C (preferred)**: β is set by the spectral geometry via the CCSvS relation. The von Neumann entropy of the spectral triple equals the spectral action:
```
S_vN(ρ_β) = Tr(h(β D_F))
```
The "natural temperature" is β = 1/Λ where Λ is the NCG cutoff. Then:
```
Σ_f = D(ρ_Λ,f ‖ I/dim) = function of y_f and Λ/v
```

**Status**: CONJECTURED. The precise form of Σ_f in terms of Yukawa couplings needs a careful calculation on the specific spectral triple. This is a concrete computation that can be done.

### II.3 The key constraint: Petz bound on masses

From Paper 1, F ≥ exp(−Σ/2). Applied to the internal channel:

```
F_f ≥ exp(−Σ_f / 2)
```

If the bound is saturated (as it is for the gravitational channel, Paper 2):

```
m_f = M_ref × (1 − exp(−Σ_f / 2))
```

For small Σ_f (light fermions):
```
m_f ≈ M_ref × Σ_f / 2
```

For large Σ_f (heavy fermions):
```
m_f ≈ M_ref × (1 − e^{−Σ_f/2}) ≈ M_ref
```

**This naturally caps the fermion mass at M_ref**, which should be identified with the electroweak scale v = 246 GeV (the Higgs VEV). The top quark (m_t ≈ 173 GeV ≈ 0.70 v) is close to saturation.

---

## III. Connection to Singh's J₃(O) Mass Ratios

### III.1 Singh's result [KNOWN, arXiv:2508.10131]

Singh (2025) showed that the complexified exceptional Jordan algebra J₃(O_C) determines fermion mass ratios with ONE free parameter δ², and the algebraic structure fixes δ² = 3/8.

The exceptional Jordan algebra J₃(O) consists of 3×3 Hermitian matrices over the octonions:
```
X = | α₁    a₃*   a₂  |
    | a₃    α₂    a₁* |
    | a₂*   a₁    α₃  |
```
where α_i ∈ R, a_i ∈ O.

Three eigenvalues of X → three generations. The "universal ladder" with parameter δ determines the eigenvalue spread:

```
λ₁ : λ₂ : λ₃ determined by δ² = 3/8
```

Key predictions (parameter-free):
- √(m_e) : √(m_u) : √(m_d) = 1 : 2 : 3
- √(m_τ/m_μ) = √(m_s/m_d) (from E₆ Dynkin diagram Z₂ automorphism / triality)
- CKM matrix elements from ladder state overlaps

### III.2 δ² = 3/8 = sin²θ_W [CONJECTURED — key new observation]

In the NCG framework, the Weinberg angle at unification scale is:

```
sin²θ_W = 3/8  [KNOWN: Chamseddine-Connes 1996]
```

This comes from the ratio of the U(1) coupling to the total:
```
sin²θ_W = Tr(Y²) / Tr(T₃² + Y²) = 3 / (5 + 3) = 3/8
```
where the trace is over one generation of SM fermions.

Singh's δ² = 3/8 comes from the octonionic algebra structure in J₃(O).

**CONJECTURE**: These are the same 3/8. Both arise from the algebraic structure of A_F = C ⊕ H ⊕ M₃(C), which is the maximal associative subalgebra of J₃(O) (Boyle 2020). The 3/8 reflects the ratio of abelian to total gauge content in the internal algebra.

If this is correct:
- Singh's mass ratios are not independent of the gauge structure — they ARE the gauge structure seen through J₃(O)
- The Weinberg angle determines the fermion mass hierarchy
- Both are manifestations of the Σ structure on J₃(O)

### III.3 Σ on J₃(O) [SPECULATIVE]

Define Σ on the Jordan algebra:

```
Σ_J(X) = D(ρ_X ‖ ρ_I)
```

where ρ_X is the state associated with the Jordan element X, and ρ_I is the trace state.

For a Jordan algebra, the natural "state" associated to X ∈ J₃(O) is:
```
ρ_X = exp(X) / Tr(exp(X))
```

The QRE is then:
```
Σ_J(X) = Tr(X · ρ_X) − ln Tr(exp(X))
```

For diagonal X = diag(λ₁, λ₂, λ₃):
```
Σ_J = Σᵢ pᵢ λᵢ − ln(Σᵢ exp(λᵢ))
```
where pᵢ = exp(λᵢ)/Σⱼ exp(λⱼ).

**The mass formula would be**:
```
m_i / M_ref = 1 − exp(−Σ_J^(i) / 2)
```

where Σ_J^(i) is the QRE contribution from the i-th eigenvalue.

### III.4 Can we derive δ² = 3/8 from Σ extremization? [SPECULATIVE]

**Proposal**: Extremize the total Σ on J₃(O) subject to:
1. Tr(X) = fixed (total "mass budget")
2. X ∈ J₃(O) (octonionic structure constraint)
3. Petz bound saturation on each eigenvalue

The extremal condition:
```
δΣ_J / δX = 0  subject to constraints
```

In the diagonal sector with eigenvalues (λ₁, λ₂, λ₃) and constraint λ₁ + λ₂ + λ₃ = S:

```
∂Σ_J/∂λᵢ = μ  (Lagrange multiplier for trace constraint)
```

This gives: all λᵢ equal (maximum entropy → all masses equal → unphysical).

**The symmetry breaking must come from the octonionic structure.** In J₃(O), the off-diagonal elements a_i ∈ O carry 8 real components each. The condition that X be a valid element of J₃(O) (with the octonionic multiplication) introduces correlations between eigenvalues that are absent in J₃(R) or J₃(C).

**Specifically**: The characteristic equation of X ∈ J₃(O) is:
```
λ³ − Tr(X)λ² + S₂(X)λ − det(X) = 0
```

where det(X) is the Freudenthal determinant (which uses octonionic multiplication and is NOT the ordinary determinant). The constraint that det(X) take a specific value — determined by the octonionic structure — restricts the eigenvalue ratios.

Singh shows that the "minimal universal ladder" in the flavor space of J₃(O) has:
```
δ² = dim(C ⊕ H) / dim(C ⊕ H ⊕ O) = (2 + 4) / (2 + 4 + 8 + 2) = 6/16 = 3/8
```

Wait — let me be more careful. Singh's δ² = 3/8 arises from the Sym³(3) representation of the flavor SU(3), specifically from the ratio of Casimir invariants. The numerical coincidence with sin²θ_W needs a precise algebraic proof.

**Status**: SPECULATIVE. The extremization program is well-defined but the computation is nontrivial. The coincidence δ² = sin²θ_W = 3/8 is suggestive but unproven to be more than numerological.

---

## IV. The Mass Formula

### IV.1 Attempt 1: Exponential suppression [CONJECTURED]

If Σ_f is determined by the Jordan eigenvalue λ_f and the Weinberg angle:

```
Σ_f = λ_f / sin²θ_W = (8/3) λ_f
```

Then:
```
m_f = M_ref × (1 − exp(−(4/3)λ_f))
```

For small λ_f:
```
m_f ≈ M_ref × (4/3)λ_f
```

The mass RATIOS are:
```
m_f / m_{f'} = (1 − exp(−(4/3)λ_f)) / (1 − exp(−(4/3)λ_{f'}))
             ≈ λ_f / λ_{f'}  (for small eigenvalues)
```

To get the observed hierarchy m_e : m_μ : m_τ = 1 : 207 : 3477, we need:
```
λ_e : λ_μ : λ_τ ≈ 1 : 207 : 3477
```

This is just restating the hierarchy in terms of Jordan eigenvalues — not explaining it.

### IV.2 Attempt 2: Double exponential [SPECULATIVE]

What if the Jordan eigenvalues themselves are exponentially spaced?

Singh's ladder structure gives eigenvalues of the form:
```
λ_n = λ_0 × (something)^n
```

If the "something" is exp(δ²) = exp(3/8):

```
λ₁ = λ₀
λ₂ = λ₀ × exp(3/8)
λ₃ = λ₀ × exp(3/4)
```

Then:
```
m₁ / m₂ = (1 − exp(−c λ₀)) / (1 − exp(−c λ₀ e^{3/8}))
```

For c λ₀ ≪ 1:
```
m₁ / m₂ ≈ e^{−3/8} ≈ 0.687
```

This gives mass ratios of order 1, not the observed 1:200. **Attempt 2 fails.**

### IV.3 Attempt 3: Power-law from Jordan trace formula [SPECULATIVE]

Singh's actual construction uses the Jordan trace formula differently. The key observation is that the squared masses (not masses) follow from the characteristic equation:

```
m²_i ∝ λ_i²
```

And the eigenvalue ratios come from the ladder operators in Sym³(3):

```
√(m_e) : √(m_u) : √(m_d) = 1 : 2 : 3
```

In the Σ framework, this would mean:
```
Σ_f ∝ m_f² / M_ref²  (quadratic, not linear)
```

or equivalently:
```
Σ_f = (y_f)² = (m_f / v)²
```

This is the **Yukawa squared**, which appears naturally in the spectral action! The CCSvS result gives:

```
S_vN ∝ Σ_k c_k (Λ/m_k)^{d-2k}
```

The fermion contribution to the spectral action at the a₀ level (cosmological constant) involves Tr(D_F⁰) = dim(H_F), while the a₂ level involves Tr(D_F²) ∝ Σ_f m_f².

**CONJECTURE**: Σ_f = (m_f / M_P)² × (M_P / v)² = y_f² × (M_P / v)².

This is the Yukawa coupling squared, multiplied by a hierarchy factor. For the top quark: Σ_t ≈ 1 × (10¹⁷)² ≈ 10³⁴. This is enormous and unphysical for the "retrodiction cost" interpretation.

### IV.4 Attempt 4: Logarithmic [CONJECTURED — most promising]

The Σ = 2 ln Q structure suggests a logarithmic relationship. In spacetime:

```
Σ_grav = −ln(−g₀₀) = 2 ln Q
```

By analogy, in the internal space:

```
Σ_f = −ln(1 − m_f² / M_ref²)
```

This has the right properties:
- Σ_f = 0 when m_f = 0 (massless = retrodictable)
- Σ_f → ∞ when m_f → M_ref (mass at reference scale = horizon)
- Σ_f ≈ m_f² / M_ref² for light fermions

**Inverting**:
```
m_f² / M_ref² = 1 − exp(−Σ_f) = τ_f²
```

So:
```
m_f = M_ref × τ_f = M_ref × √(1 − exp(−Σ_f))
```

For small Σ_f:
```
m_f ≈ M_ref × √(Σ_f)
```

**If M_ref = v = 246 GeV** (electroweak scale), then:
- m_e = 0.511 MeV → τ_e = m_e/v = 2.08 × 10⁻³ → Σ_e = τ_e² = 4.3 × 10⁻⁶
- m_μ = 106 MeV → τ_μ = 4.31 × 10⁻¹ → Σ_μ = 1.9 × 10⁻¹ (Wait, this doesn't work)

Actually m_μ/v = 4.31 × 10⁻⁴, giving Σ_μ = 1.86 × 10⁻⁷. Let me redo:

- m_e/v = 2.08 × 10⁻⁶ → Σ_e ≈ 4.3 × 10⁻¹²
- m_μ/v = 4.31 × 10⁻⁴ → Σ_μ ≈ 1.86 × 10⁻⁷
- m_τ/v = 7.22 × 10⁻³ → Σ_τ ≈ 5.22 × 10⁻⁵
- m_t/v = 0.703 → Σ_t ≈ 0.788

Ratios:
- Σ_μ/Σ_e ≈ 4.3 × 10⁴ ≈ (m_μ/m_e)² = (207)² = 4.3 × 10⁴ ✓ (trivially, since Σ ≈ m²/v²)
- Σ_τ/Σ_μ ≈ 281 ≈ (m_τ/m_μ)² = (16.8)² = 282 ✓

This is just restating m_f = v × y_f in logarithmic language. The question is whether the Σ_f values have a deeper origin in J₃(O).

### IV.5 The honest assessment of Attempts 1–4

**None of these attempts DERIVE the mass hierarchy.** They rewrite it in Σ language. The actual mass prediction requires knowing the eigenvalues of D_F on the specific spectral triple of the SM — which is the full content of the Yukawa sector.

**What the Σ framework CAN do**:
1. Provide a physical INTERPRETATION (mass = retrodiction cost)
2. Provide CONSTRAINTS (Petz bound → mass ≤ M_ref)
3. Provide a BRIDGE to J₃(O) (if Σ on J₃(O) = spectral action on A_F)

**What requires external input (Singh or similar)**:
- The specific eigenvalue spectrum of J₃(O) → mass ratios
- The value δ² = 3/8
- The CKM matrix elements

---

## V. The δ² = 3/8 = sin²θ_W Coincidence

### V.1 Why sin²θ_W = 3/8 in NCG [KNOWN]

In the Chamseddine-Connes framework, the gauge couplings at the GUT scale satisfy:

```
g₁² Tr(Y²) = g₂² Tr(T₃²) = g₃² Tr(T_a²)
```

where the traces are over one generation. For the SM fermion content:

```
Tr(Y²) = 2 × (1/3 × 3 × (4/9 + 1/9) + 1/3 × 1 × (1 + 0) + ...) = 10/3
Tr(T₃²) = 2 × (1/4 × 3 × 2 + 1/4 × 1 × 2) = 2
```

Wait, let me be more careful. The standard NCG normalization gives:

```
sin²θ_W = g₁²/(g₁² + g₂²) = Tr(T₃²)/(Tr(T₃²) + Tr(Y²))
```

With the specific SM representations:
```
sin²θ_W = 3/(3 + 5) = 3/8
```

This is a consequence of the algebra A_F = C ⊕ H ⊕ M₃(C) and nothing else.

### V.2 Why δ² = 3/8 in Singh [KNOWN from arXiv:2508.10131]

Singh's δ² parameterizes the "spread" of eigenvalues in the Jordan element. He determines it from the Sym³(3) representation of the flavor SU(3) acting on J₃(O_C).

His derivation (simplified):
1. The 27-dimensional representation of E₆ decomposes under SU(3)_flavor as 27 = 10 + 8 + 6 + 3
2. The "ladder" connecting generations lives in the Sym³(3) = 10 representation
3. The universal step δ is constrained by the octonionic multiplication to satisfy δ² = 3/8

### V.3 Are they the same 3/8? [CONJECTURED]

**Argument FOR**:
1. Both arise from A_F = C ⊕ H ⊕ M₃(C), which IS the maximal associative subalgebra of J₃(O) (Boyle 2020)
2. Both involve the ratio of abelian to total gauge content
3. The NCG framework UNIFIES gauge couplings and Yukawa couplings in D_F; sin²θ_W and mass ratios come from the SAME operator
4. The Σ framework says gauge forces AND masses come from the same QRE; if Σ on A_F = spectral action, then the 3/8 in gauge couplings and the 3/8 in mass ratios are the same Σ structure

**Argument AGAINST**:
1. Singh's derivation uses Sym³(3) of flavor SU(3), not the gauge SU(3)_color. The "3" in his 3/8 is the number of generations, not the number of weak isospin components
2. The NCG derivation uses the SM fermion content specifically; Singh's uses the abstract J₃(O) structure. They might give the same number by accident
3. sin²θ_W runs with energy. At low energy, sin²θ_W ≈ 0.231 ≠ 3/8. The mass ratios are measured at low energy. The 3/8 applies only at the GUT scale

**Assessment**: The coincidence is highly suggestive but the algebraic proof of equivalence does not exist. This is the key open question for Paper 10.

### V.4 What would a proof look like?

**Required**: Show that the characteristic equation of X ∈ J₃(O) subject to the constraint that Σ_J(X) be extremized produces eigenvalue ratios that match the SM fermion mass ratios, AND that the extremal Σ_J involves the factor 3/8 from the octonionic structure in a way that reduces to sin²θ_W when restricted to the gauge sector.

**Concretely**: Start with the spectral action on M × F:
```
S = Tr(f(D²/Λ²))
```

The Dirac operator is D = D_M ⊗ 1 + γ₅ ⊗ D_F, where D_F contains Yukawa couplings.

The CCSvS result: S_vN(ρ_β) = Tr(h(β D)).

If Σ on A_F ↔ first variation of S_vN on F, then:
```
δΣ_F = 0  ⟺  δTr(h(β D_F)) = 0
```

The extremal D_F is the physical Dirac operator with the observed Yukawa couplings.

**The 3/8 would appear as**: The ratio of Σ_gauge (from gauge bosons) to Σ_total (gauge + Yukawa) at the extremum equals 3/8 if the A_F structure forces this ratio.

---

## VI. Mass Predictions: What Can We Actually Get?

### VI.1 The Σ-Singh hybrid formula [SPECULATIVE]

Combining Σ interpretation with Singh's eigenvalues:

**Step 1**: Identify the J₃(O) eigenvalues with Σ_f:
```
Σ_f = f(λ_f)  for some monotonic f
```

**Step 2**: The simplest choice f = identity gives:
```
m_f = M_ref × (1 − exp(−λ_f / 2))
```

**Step 3**: Singh's eigenvalue ratios (for charged leptons):
```
λ_e : λ_μ : λ_τ = (1-δ)² : 1 : (1+δ)²
```
with δ² = 3/8, so δ = √(3/8) ≈ 0.6124.

```
λ_e ∝ (1 − 0.6124)² = 0.1503
λ_μ ∝ 1
λ_τ ∝ (1 + 0.6124)² = 2.5927
```

Ratios: λ_e : λ_μ : λ_τ = 1 : 6.653 : 17.25

**But observed**: m_e : m_μ : m_τ = 1 : 207 : 3477

These don't match at all. Singh's actual construction is more sophisticated — the ladder operates on SQUARED masses, not masses, and involves multiple steps.

### VI.2 Singh's actual mass formula [KNOWN]

Singh uses:
```
√(m_e) : √(m_u) : √(m_d) = 1 : 2 : 3
```

as the INPUT for the first generation, and then the generation ratios come from the Sym³(3) ladder:
```
m_{2nd gen} / m_{1st gen} = function of δ² = 3/8
```

His predictions (compared to data):

| Ratio | Singh prediction | Observed | Agreement |
|-------|-----------------|----------|-----------|
| m_μ/m_e | ~207 | 206.8 | ✓ excellent |
| m_τ/m_μ | ~16.8 | 16.8 | ✓ excellent |
| m_c/m_u | ~278 | ~277 | ✓ excellent |
| m_t/m_c | ~135 | ~136 | ✓ excellent |
| m_s/m_d | ~20 | ~20 | ✓ excellent |
| m_b/m_s | ~50 | ~49 | ✓ excellent |

These are remarkable. The key question: can Σ provide the PRINCIPLE that selects Singh's construction?

### VI.3 The Σ principle for Singh's mass formula [SPECULATIVE]

**Proposal**: Singh's "minimal universal ladder" IS the Petz recovery map on J₃(O).

The Petz map R recovers information. In J₃(O):
- The Petz map on the diagonal (eigenvalue) sector: R maps λ_i → λ_i (trivial)
- The Petz map on the off-diagonal (flavor-mixing) sector: R introduces the CKM-like overlaps

**The "minimal ladder"** = the channel that loses the least information while connecting different eigenvalue sectors. This is precisely what the Petz map does — it is the OPTIMAL retrodiction channel.

If this identification holds:
1. The ladder structure follows from Petz optimality
2. The step size δ follows from the algebraic constraints of J₃(O)
3. The 3/8 follows from the octonionic structure
4. The mass ratios follow from the Σ hierarchy of the Petz recovery on J₃(O)

**Status**: Highly speculative. This requires proving that the Petz map on J₃(O) reproduces Singh's ladder construction. No such proof exists.

---

## VII. Connection to Paper 9: Σ on A_F = Spectral Action

### VII.1 The bridge [CONJECTURED]

Paper 9 showed that the Fisher metric of the spectral action matches Σ = 2 ln Q in d = 4. Extending this to the internal space:

**Paper 9 result (gravity sector)**:
```
g_QQ^(EH) = (d-2)(d-3)/Q² = 2/Q²  in d=4
```

**Internal sector analog**: The Fisher metric of the spectral action on A_F with respect to the Yukawa couplings:
```
g_{yy}^(F) = ∂²S_F/∂y_f² = ?
```

The spectral action on the internal space is:
```
S_F = Tr(f(D_F²/Λ²)) ≈ f₀ Tr(D_F⁰) + f₂ Λ² Tr(D_F²) + f₄ Λ⁴ Tr(D_F⁴) + ...
```

The leading mass-dependent term:
```
Tr(D_F²) = Σ_f multiplicity(f) × m_f² = Σ_f N_f y_f² v²
```

This gives:
```
g_{yy}^(F) = ∂²(Σ_f N_f y_f² v²)/∂y_f² = 2 N_f v²
```

This is just a constant — the internal Fisher metric is flat in Yukawa space! This means Σ_f ∝ y_f² (quadratic, as in Attempt 4 of Section IV).

### VII.2 What this means for masses

If Σ_f = c × y_f² for some constant c (set by the NCG normalization), then:

```
τ_f = 1 − exp(−c y_f²/2)
```

The mass hierarchy is encoded in the Yukawa hierarchy, which IS the unsolved problem of the SM. The Σ framework recasts it as a retrodiction cost hierarchy but does not solve it without additional input (Singh's J₃(O) structure or equivalent).

### VII.3 Where Paper 10 adds value

Paper 10's contribution is NOT to derive fermion masses from nothing. It is to:

1. **Interpret** masses as retrodiction costs (new physical picture)
2. **Bridge** Singh's J₃(O) mass ratios to the Σ framework (connecting the algebraic and information-theoretic programs)
3. **Constrain** masses via the Petz bound (M_ref as upper bound)
4. **Unify** the 3/8 (if the coincidence is real)
5. **Predict** neutrino masses (the lightest = smallest Σ_internal)

---

## VIII. Neutrino Mass from Σ

### VIII.1 The seesaw in Σ language [CONJECTURED]

In the standard seesaw, light neutrino masses are:
```
m_ν ≈ y_ν² v² / M_R
```

In Σ language:
```
Σ_ν = y_ν² × (some constant) ≈ (m_ν M_R / v²) × c
```

If M_R ~ M_GUT ~ 10¹⁶ GeV:
```
Σ_ν ~ m_ν × 10¹⁶ / v² ~ m_ν × 10¹¹ / GeV
```

For m_ν ~ 0.05 eV: Σ_ν ~ 5 × 10⁻¹² × 10¹¹ ~ 5 × 10⁻¹. So Σ_ν ~ 0.5, which is O(1). Interesting — the neutrino is NOT the lightest in Σ space, because the seesaw amplifies the Σ back up.

### VIII.2 Alternative: Neutrino mass from cosmological Σ [SPECULATIVE]

In the mu framework (Papers 3-4):
```
μ·c² = a₀(1 + z_dec)
```

If μ ~ H₀/c is the fundamental mass scale, and neutrino mass is the "cosmological retrodiction cost":
```
m_ν ~ μ c² × (some factor from Σ_cosmological)
```

With μ c² ~ 10⁻³³ eV, this needs a large amplification. The seesaw provides it naturally.

### VIII.3 DESI constraint

DESI (2024-2025): Σ m_ν < 0.05-0.07 eV.

In the Σ framework: Σ_ν must be small enough that τ_ν < threshold. This constrains the minimum Σ in the internal space.

**Status**: SPECULATIVE for the cosmological neutrino mass; CONJECTURED for the seesaw in Σ language.

---

## IX. The CKM Matrix from Retrodiction Overlap

### IX.1 Singh's CKM from J₃(O) [KNOWN]

Singh obtains CKM matrix elements as overlaps between the "ladder states" of the up and down sectors in J₃(O):

```
V_{ij} = ⟨u_i | d_j⟩_J₃(O)
```

where |u_i⟩ and |d_j⟩ are the i-th and j-th generation states in the up and down sectors respectively.

### IX.2 CKM from Petz recovery [SPECULATIVE]

In the Σ framework, the CKM matrix should be the "fidelity matrix" between the Petz recovery channels of the up and down sectors:

```
|V_{ij}|² = F(R_up^(i), R_down^(j))
```

where R_up^(i) is the Petz recovery map for the i-th generation up-type quark, and F is the channel fidelity.

**Physical meaning**: |V_{ij}|² measures how well the retrodiction of the i-th up-type quark matches the retrodiction of the j-th down-type quark. If the retrodictions are "similar" (high fidelity), the CKM element is large (= the quarks mix). If they are "dissimilar" (low fidelity), the CKM element is small (= the quarks decouple).

This would explain:
- |V_ud| ≈ 1: up and down have nearly identical retrodiction channels (same generation, similar mass)
- |V_ub| ≈ 0.004: up and bottom have very different retrodiction channels (different generations, very different mass)
- CP violation: the complex phase of V = the phase acquired by the Petz map when acting across generations

**Status**: SPECULATIVE. No computation exists.

---

## X. Chirality from Σ Gradient

### X.1 The domain wall picture [CONJECTURED]

In the NCG framework, the Dirac operator D = D_M ⊗ 1 + γ₅ ⊗ D_F has a natural domain-wall structure:

- γ₅ is the chirality operator on M
- D_F is the "mass matrix" in the internal space
- The product γ₅ ⊗ D_F couples chirality to mass

In the Σ framework, the Σ gradient in the internal space defines a preferred direction. Fermions localized on the Σ = 0 surface (domain wall) are chiral:

```
Left-handed: ψ_L localized where ∇Σ_internal points "inward"
Right-handed: ψ_R localized where ∇Σ_internal points "outward"
```

The mass term is the tunneling amplitude across the domain wall:
```
m_f ∝ exp(−Σ_wall / 2)
```

where Σ_wall is the "width" of the domain wall in Σ space. Light fermions have large Σ_wall (wide domain wall, small tunneling). Heavy fermions have small Σ_wall (narrow domain wall, large tunneling).

This is the Kaplan mechanism (1992) rewritten in Σ language.

---

## XI. Summary: What's Proven, What's Conjectured, What's Speculative

### PROVEN
1. τ = 1 − F, Σ = D(ρ‖σ) ≥ 0, F ≥ exp(−Σ/2) — Paper 1
2. Σ_grav = 2 ln Q on static backgrounds — Paper 2
3. g_QQ = (d−2)(d−3)/Q² from spectral action — Paper 9
4. sin²θ_W = 3/8 at GUT scale in NCG — Chamseddine-Connes 1996
5. Spectral action = von Neumann entropy — CCSvS 2018
6. Singh's mass ratios from J₃(O) with δ² = 3/8 — Singh 2025

### CONJECTURED (testable with explicit computation)
7. Σ_f = D(ρ_f ‖ ρ_0) on A_F defines the fermion's retrodiction cost
8. m_f = M_ref × (1 − exp(−Σ_f/2)) with M_ref = v
9. Σ_f ∝ y_f² from the spectral action structure
10. δ² = sin²θ_W is the SAME algebraic quantity seen from J₃(O) vs A_F
11. The Petz bound constrains the maximum fermion mass to ≤ v

### SPECULATIVE (requires major new work)
12. Σ extremization on J₃(O) selects the physical mass ratios
13. Singh's ladder = Petz recovery map on J₃(O)
14. CKM matrix = Petz fidelity matrix between up/down channels
15. Chirality from Σ gradient (domain wall mechanism in internal space)
16. Neutrino mass from minimum Σ_internal

---

## XII. Concrete Next Steps

### Step 1: Compute Σ_f explicitly (Week 1-4)
On the SM spectral triple (A_F, H_F, D_F):
- Define ρ_f as thermal state of D_F restricted to fermion f
- Compute D(ρ_f ‖ I/dim) for each SM fermion
- Check if ratios match known mass ratios (they should, trivially, since D_F encodes masses)

### Step 2: Test the 3/8 connection (Week 4-8)
- Show algebraically whether Singh's δ² = 3/8 from Sym³(3) of J₃(O) is the same as sin²θ_W = 3/8 from Tr(Y²)/Tr(T₃²+Y²)
- Key question: does the Boyle identification A_F ⊂ J₃(O) map one 3/8 to the other?

### Step 3: Petz on J₃(O) (Month 3-6)
- Define the Petz recovery map on J₃(O)
- Check if it reproduces Singh's ladder
- This requires careful treatment of the Jordan algebra structure (non-associative)

### Step 4: CKM from fidelity (Month 6-12)
- Define channel fidelity between up/down Petz maps
- Compute |V_ij|²
- Compare with observed CKM elements

### Step 5: Neutrino prediction (Month 6-12)
- Use minimum Σ_internal principle to predict lightest neutrino mass
- Compare with DESI bounds

---

## XIII. Key References

1. **Singh** (2025), arXiv:2508.10131 — J₃(O) mass ratios with δ² = 3/8
2. **Furey & Hughes** (2024), arXiv:2409.17948 — Three generations from trialities
3. **Furey** (2025), arXiv:2505.07923 — Z₂⁵-graded superalgebra
4. **Chamseddine-Connes-van Suijlekom** (2018/2020), arXiv:1809.02944 — Spectral action = entropy
5. **Chamseddine-Connes** (1996), hep-th/9606001 — NCG spectral action
6. **Boyle** (2020), arXiv:2006.16265 — A_F ⊂ J₃(O)
7. **Farnsworth, Finster, Paganini, Singh** (2026), arXiv:2603.05018 — CFS + NCG + trace dynamics bridge
8. **Dorau-Much** (2025), arXiv:2510.24491 — QRE ⟹ Einstein equations (PRL)
9. **Bianconi** (2025), PRD 111, 066001 — Gravity from entropy
10. **Szangolies** (2025), arXiv:2512.17328 — n-qubit → SM gauge group
11. **Kaplan** (1992), Phys.Lett.B288:342 — Domain wall fermions

---

*Last updated: 2026-03-23*
*Status: Research note for Paper 10. Establishes the research program. Key open question: is δ² = sin²θ_W a coincidence or a theorem?*
