# Paper 9 to the Standard Model: What g_QQ = 2/Q^2 Implies for SU(3) and Fermions
## 2026-03-20
## Sheng-Kai Huang

---

## 0. Executive Summary

Paper 9 proved:

```
g_QQ^{spectral}(Q) = (d-2)(d-3)/Q^2   (for the a_2 sector)
g_QQ^{spectral}(Q) = 0                  (for the a_4 sector)
```

In d=4: the a_2 (Einstein-Hilbert + Higgs mass) sector gives g_QQ = 2/Q^2, matching Sigma = 2 ln Q. The a_4 (Yang-Mills + Higgs kinetic + Higgs quartic) sector gives g_QQ = 0 (conformally invariant).

**This document asks: what does this split tell us about SU(3), fermions, and the Standard Model?**

Key conclusions:

1. **g_QQ(a_4) = 0 is the conformal invariance of classical Yang-Mills.** This is not a constraint on the gauge group -- it holds for ANY gauge group in d=4. SU(3) is not selected by conformal invariance alone.

2. **g_QQ(a_2) = 2/Q^2 does NOT directly constrain A_F.** The Fisher metric of a_2 is 2/Q^2 regardless of which internal algebra A_F is chosen, because a_2 scales as Q^{d-2} universally.

3. **However, g_QQ = 2 SELECTS d=4.** And d=4 is the dimension where the NCG constraints on A_F are most restrictive. The combination (d=4) + (NCG axioms) + (Sigma framework) provides a three-pronged selection.

4. **Fermion generations: the condition (d-2)(d-3) = 2 selects d=4. No analogous single equation selects N_gen = 3.** But J_3(O) triality (3 eigenvalues) and SO(8) triality (3 inequivalent 8-dimensional representations) remain the best candidates, and these structures only work in d=4 + NCG.

5. **The real power of Paper 9 for the SM is INDIRECT**: it establishes Sigma as the master functional, then the NCG spectral triple (A_F, H_F, D_F) determines the SM content. Paper 9 does not derive SU(3) or fermions; it provides the framework within which those derivations (Papers 7-8, 10) become meaningful.

---

## 1. Question 1: Why is g_QQ(a_4) = 0, and What Does It Mean?

### 1.1 The Mathematical Fact

The a_4 Seeley-DeWitt coefficient in d=4 is conformally invariant:

```
a_4[Q^2 g] = Q^{d-2k} a_4[g] |_{k=2, d=4} = Q^0 a_4[g] = a_4[g]
```

Since a_4 is independent of Q, all Q-derivatives vanish, so g_QQ = 0.

### 1.2 Connection to Conformal Invariance of Classical Yang-Mills

The a_4 coefficient contains (Chamseddine-Connes 1996):

```
a_4 = (1/360)(4pi)^{-2} integral d^4 x sqrt{g} [
    curvature^2 terms (Gauss-Bonnet + Weyl^2)
    + alpha_i F_{mu nu}^{(i)} F^{mu nu (i)}   (Yang-Mills for each gauge factor)
    + |D_mu H|^2                                (Higgs kinetic)
    + lambda |H|^4                              (Higgs quartic)
]
```

The classical Yang-Mills action S_YM = integral sqrt{g} |F|^2 IS conformally invariant in d=4. Under g -> Q^2 g with constant Q:

```
sqrt{g-hat} |F-hat|^2 = Q^4 sqrt{g} * Q^{-4} |F|^2 = sqrt{g} |F|^2
```

The Q^4 from the volume exactly cancels the Q^{-4} from raising indices. This is a well-known result (Weyl 1918, Penrose 1963).

**This conformal invariance holds for ANY Yang-Mills theory in d=4: U(1), SU(2), SU(3), SU(5), SO(10), E_8, anything.** It does not constrain the gauge group.

### 1.3 What g_QQ(a_4) = 0 DOES Tell Us

The vanishing Fisher metric means: **the gauge coupling constants are not affected by uniform conformal rescaling.** This has two important consequences:

**Consequence A: Scale separation.** Gravity (a_2, g_QQ = 2/Q^2) responds to conformal deformation; gauge forces (a_4, g_QQ = 0) do not. This is the spectral-action version of the statement that gauge couplings are dimensionless while Newton's constant is dimensionful. The Sigma framework "sees" gravity but is "blind" to gauge forces under conformal variation.

**Consequence B: Gauge forces require a different probe.** To detect gauge structure via Sigma, one cannot use conformal variation. One must use INTERNAL fluctuations of the Dirac operator:

```
D -> D + A + J A J^{-1}
```

where A is a 1-form in the NCG sense (inner fluctuation) and J is the real structure operator. These inner fluctuations generate the gauge fields. The Fisher metric with respect to GAUGE parameters (not Q) is where SU(3) x SU(2) x U(1) enters.

This is already the content of Paper 6 for U(1): the EM sector enters through the PHASE of the complex Khronon, not through the conformal amplitude. For SU(2) and SU(3), the relevant parameters are the non-Abelian gauge potentials.

### 1.4 Summary for Question 1

| Statement | True? | Why |
|-----------|-------|-----|
| g_QQ(a_4) = 0 because Yang-Mills is conformally invariant in d=4 | **Yes** | Exact cancellation of volume and metric factors |
| This constrains the gauge group | **No** | Conformal invariance holds for ALL gauge groups in d=4 |
| Gauge structure is invisible to Sigma under conformal variation | **Yes** | Sigma "sees" a_2 (gravity) not a_4 (gauge) via Q |
| Gauge structure enters through inner fluctuations D -> D + A | **Yes** | NCG: gauge = inner automorphisms of A_F |

---

## 2. Question 2: Does g_QQ(a_2) = 2/Q^2 Constrain A_F?

### 2.1 The Direct Answer: No

The scaling a_2[Q^2 g] = Q^2 a_2[g] is purely geometric. It depends on the dimension d and the heat kernel coefficient index k, not on the internal algebra A_F.

For the full SM spectral triple (M^4 x F with A_F = C + H + M_3(C)):

```
a_2^{full} = (4pi)^{-2} integral d^4 x sqrt{g} [c_R R + c_H mu^2 |H|^2]
```

where c_R and c_H depend on the spectral data of A_F (traces of Yukawa couplings, etc.). But under Q^2 g:

```
a_2^{full}[Q^2 g] = Q^2 a_2^{full}[g]
```

The Q^2 scaling is UNIVERSAL -- the internal algebra only affects the overall normalization (the coefficients c_R, c_H), not the scaling exponent. Therefore:

```
g_QQ(a_2) = (d-2)(d-3)/Q^2 = 2/Q^2   for ALL choices of A_F in d=4
```

### 2.2 The Indirect Answer: Yes, Through d=4 Selection

While g_QQ = 2/Q^2 does not constrain A_F directly, it SELECTS d=4. And d=4 is the dimension where the NCG axioms most severely constrain A_F.

The NCG classification of irreducible finite spectral triples in d=4 (Chamseddine-Connes 2008; Krajewski 2000) shows that the 7 axioms (dimension, orientation, chirality, Poincare duality, reality, first order, regularity) allow only a very restricted set of algebras A_F. The Standard Model algebra A_F = C + H + M_3(C) is essentially unique (up to the number of generations, which is an input).

**The chain of inference:**

```
Sigma = 2 ln Q  (quantum channel)
       |
       v
g_QQ = 2/Q^2  (Fisher metric)
       |
       v
(d-2)(d-3) = 2  (from spectral action)
       |
       v
d = 4  (unique non-trivial solution)
       |
       v
NCG axioms in d=4  (severe constraints)
       |
       v
A_F = C + H + M_3(C)  (essentially unique)
       |
       v
G_SM = U(1) x SU(2) x SU(3) / Z_6
```

This is a SELECTION chain, not a derivation chain. Each step narrows the possibilities. The Sigma framework provides the first link (why d=4); the NCG axioms provide the rest.

### 2.3 What Would Make This a Derivation

For this to be a genuine derivation of A_F from Sigma, we would need:

**Missing step**: Show that the NCG axioms themselves follow from information-theoretic principles (QRE positivity, Petz recovery, DPI). This is an open problem.

Partial progress: Dorau-Much (2025, PRL) derived the Einstein field equations from QRE. If the NCG axioms can similarly be derived from QRE constraints on the spectral triple, then A_F = C + H + M_3(C) would follow from Sigma alone.

---

## 3. Question 3: Can Sigma Select N_gen = 3?

### 3.1 The Dimension Analogy

The condition (d-2)(d-3) = 2 is a quadratic equation with solutions d = 4 and d = 1 (trivial). This elegantly selects the physical spacetime dimension.

Is there an analogous condition that selects N_gen = 3?

### 3.2 Candidate Conditions

**Candidate A: From J_3(O) eigenvalues.**

The exceptional Jordan algebra J_3(O) consists of 3x3 Hermitian matrices over the octonions. It has:
- 3 real eigenvalues (corresponding to 3 generations)
- Automorphism group F_4 (52-dimensional)
- The identity component of the structure group is F_4

The number "3" in J_3(O) comes from the 3x3 matrix size. The question is: why 3x3 and not 2x2 or 4x4?

Answer (mathematical): J_2(O) is not exceptional (it is isomorphic to the spin factor V^{10}). J_4(O) does not exist as a Jordan algebra (the Jordan identity fails for 4x4 octonionic matrices due to non-associativity). So J_3(O) is the UNIQUE exceptional Jordan algebra.

**This gives a selection of "3" from algebraic uniqueness, not from a Fisher metric equation.**

**Candidate B: From SO(8) triality.**

SO(8) has a unique triality automorphism that permutes:
- 8_v (vector representation)
- 8_s (positive spinor)
- 8_c (negative spinor)

These three 8-dimensional representations are all inequivalent but related by an outer automorphism of order 3. In the Boyle (2020) framework, each representation corresponds to one generation of fermions.

The number "3" comes from the order of the triality automorphism, which is |Out(D_4)| = |S_3| = 6 (the symmetric group), but the cyclic part is Z_3.

**Candidate C: From Sigma extremal principle on J_3(O).**

Conjecture: The QRE Sigma = D(rho || sigma) on the state space of J_3(O) has exactly 3 extremal directions (corresponding to the 3 eigenvalues). These 3 directions correspond to the 3 fermion generations.

Status: UNVERIFIED. This requires computing the QRE landscape on J_3(O), which is a well-defined but computationally challenging problem.

**Candidate D: Fisher metric condition on the generation space.**

By analogy with (d-2)(d-3) = 2, one could ask: is there a Fisher metric condition f(N_gen) = (something fixed) that selects N_gen = 3?

If the generation structure enters through a heat kernel coefficient a_{2k} that scales as Q^{alpha(N_gen)}, then the Fisher metric would be alpha(alpha-1)/Q^2. Setting this equal to some target value could select N_gen.

However, in the standard NCG framework, N_gen enters as a MULTIPLICATIVE factor in the traces (tr_F -> N_gen * tr_{1-gen}), not as an exponent. So the Fisher metric is proportional to N_gen, not a quadratic function of it. This means:

```
g_QQ(N_gen) = N_gen * g_QQ(1-gen)
```

Setting g_QQ = 2 gives N_gen = 2/g_QQ(1-gen), which requires knowing g_QQ for a single generation -- and this is NOT an integer in general.

**Verdict: No clean analog of (d-2)(d-3) = 2 for N_gen has been found.**

### 3.3 The Honest Assessment

| Mechanism for N_gen = 3 | Status | Strength |
|--------------------------|--------|----------|
| J_3(O) uniqueness (algebraic) | Well-established math | Strong (algebraic necessity) |
| SO(8) triality | Well-established math | Strong (representation theory) |
| Sigma extremal on J_3(O) | Unverified conjecture | Unknown |
| Fisher metric condition | No viable candidate | Weak |
| Furey Z_2^5 grading | Recent (2025) | Promising but needs verification |
| Gourlay-Gresnigt Cl(10) + S_3 | Recent (2026) | Alternative route |

The most promising path to "3 generations from Sigma" is through J_3(O):

```
Sigma on J_3(O)  (Paper 10)
     |
     v
Spectral action on (M^4 x F) with A_F subset J_3(O)
     |
     v
D_F has 3 sets of eigenvalues (3 generations)
     |
     v
Mass ratios from D_F eigenvalue ratios (Singh 2025)
```

This route uses the UNIQUENESS of J_3(O) (not a Fisher metric condition) to get 3.

---

## 4. The a_2 vs a_4 Split: Structural Implications for the SM

### 4.1 The Two Sectors of the Spectral Action

Paper 9 reveals a clean division of the spectral action:

```
GRAVITY SECTOR (a_2):
  - Einstein-Hilbert: (1/16piG) integral sqrt{g} R
  - Higgs mass: mu^2 |H|^2
  - Fisher metric: g_QQ = 2/Q^2
  - Responds to conformal deformation
  - Sigma = 2 ln Q applies here

GAUGE SECTOR (a_4):
  - Yang-Mills: alpha_i |F_i|^2 for i = U(1), SU(2), SU(3)
  - Higgs kinetic: |D_mu H|^2
  - Higgs quartic: lambda |H|^4
  - Gauss-Bonnet: topological
  - Fisher metric: g_QQ = 0
  - Conformally invariant
  - Invisible to Sigma under conformal variation
```

### 4.2 The Physical Meaning

**The Higgs mass term lives in a_2 (gravity sector), not a_4 (gauge sector).**

This is the spectral-action perspective on the hierarchy problem: the Higgs mass is a GRAVITATIONAL quantity (it scales with the curvature), while the Higgs self-coupling is a GAUGE quantity (conformally invariant).

In the Sigma framework:
- The Higgs mass mu^2 is determined by the conformal mode Q (Fisher metric 2/Q^2)
- The gauge couplings g_1, g_2, g_3 are determined by the inner fluctuations A (Fisher metric from a_4 under A-variation, not Q-variation)
- The Higgs quartic lambda is determined by the a_4 sector

This suggests a natural decomposition of the SM Lagrangian:

```
L_SM = L_{Sigma-visible}(Q) + L_{Sigma-invisible}(A)
```

where "Sigma-visible" means contributing to g_QQ, and "Sigma-invisible" means conformally invariant.

### 4.3 Implications for Papers 7-8 (SU(2) and SU(3))

Since the gauge sector is invisible to conformal Sigma, Papers 7 and 8 CANNOT derive SU(2) and SU(3) from conformal variation of Sigma. They must use a different mechanism.

**The correct mechanism**: Sigma on INTERNAL (non-conformal) degrees of freedom.

For Paper 7 (SU(2)):
- The quaternionic structure of H subset A_F provides 3 non-conformal directions (i, j, k)
- Sigma decomposition under SU(2) twirling: D(rho || G_{SU(2)}(rho)) measures SU(2)-breaking
- Verified numerically (2026-03-19): D_break = ln(2j+1) for spin-j sector

For Paper 8 (SU(3)):
- The octonionic structure of M_3(C) subset A_F provides the color directions
- Sigma decomposition under SU(3) twirling: D(rho || G_{SU(3)}(rho)) measures color-breaking
- Confinement: D_color -> infinity for isolated quarks (Petz recovery failure)
- NOT YET COMPUTED

### 4.4 The Full Sigma Decomposition

Combining conformal and internal variations:

```
Sigma_total = Sigma_conformal(Q) + Sigma_internal(A)

Sigma_conformal = 2 ln Q                (from a_2, Papers 2-4)
Sigma_internal  = D(rho || G_G(rho))    (from inner fluctuations, Papers 6-8)
```

where G is the gauge group G_SM = U(1) x SU(2) x SU(3).

This decomposition is the FULL content of the spectral action from the Sigma perspective:
- The first term gives gravity + Higgs mass
- The second term gives gauge forces + Higgs kinetic + Higgs quartic

---

## 5. What g_QQ = 2/Q^2 Specifically Implies for Each SM Sector

### 5.1 For Gravity: CONFIRMED

g_QQ = 2/Q^2 from a_2 IS the content of Papers 2-4. The Einstein-Hilbert action, exponential metric, dark matter from Khronon -- all follow from the a_2 sector having Fisher metric 2/Q^2.

### 5.2 For the Higgs Mass: PREDICTION

The Higgs mass term mu^2 |H|^2 lives in a_2 and has g_QQ = 2/Q^2. This means:

```
mu^2(Q) = Q^2 * mu^2(1)
```

Under conformal rescaling, the Higgs mass scales linearly with Q^2. In the Sigma framework, this connects to the running of the Higgs mass with the cosmological scale factor:

```
mu^2(a) ~ a^2 * mu^2(a_0)   (if Q ~ a/a_0)
```

This is the spectral-action version of the statement that the Higgs mass is quadratically sensitive to the UV cutoff (hierarchy problem). The Sigma framework may provide a resolution through the ghost condensation mechanism (Paper 4, CW portal).

### 5.3 For Yang-Mills (SU(3), SU(2), U(1)): INVISIBLE TO Q

The gauge couplings are conformally invariant and do not appear in g_QQ. To constrain them, one must:

1. Compute the Fisher metric with respect to gauge parameters A (inner fluctuations)
2. This requires the full Dirac operator D + A + J A J^{-1} on M^4 x F
3. The result will depend on A_F and hence on the gauge group

**This is the key computation for Paper 9 -> Papers 7-8.**

### 5.4 For Fermions: REQUIRED BY NCG, NOT BY g_QQ

The Fisher metric g_QQ = 2/Q^2 does not directly involve fermions. Fermions enter through:

1. The Hilbert space H_F of the spectral triple (where fermions live)
2. The Dirac operator D_F on H_F (which determines fermion masses)
3. The spectral action Tr f(D^2/Lambda^2) which, through CCSvS, is the entropy that gives g_QQ

The fermion content is ENCODED in the spectral data but not DETERMINED by g_QQ alone.

---

## 6. Concrete Research Plan: Paper 9 -> SU(3) + Fermions

### Phase 1: Establish the Framework (Months 1-3) -- Paper 9

**Computation 1.1: Toy model S^1 x M_2(C)**
- Spectral triple: (C^inf(S^1) tensor M_2(C), L^2(S^1) tensor C^2, D = i d/dx tensor 1 + gamma tensor D_F)
- Compute S_vN (CCSvS formula) and Sigma = QRE between two Gibbs states
- Compare first variations: delta S_vN vs delta Sigma
- Expected: agreement at second order (Fisher level), disagreement at higher orders
- Go/No-Go: if the Fisher metrics match, PROCEED

**Computation 1.2: Fisher metric under inner fluctuations**
- On S^1 x M_2(C), add inner fluctuation A = a[D, b] with a, b in M_2(C)
- Compute g_{AA} (Fisher metric with respect to gauge parameter)
- Expected: g_{AA} ~ 1/(gauge coupling)^2, related to a_4 coefficient
- This tests whether Sigma "sees" gauge structure through internal (not conformal) variation

**Computation 1.3: Extension to A_F = C + M_2(C)**
- This gives U(1) x SU(2) gauge group (electroweak without strong)
- Compute a_2 and a_4 contributions to Fisher metric
- Verify: g_QQ(a_2) = 2/Q^2, g_QQ(a_4) = 0, g_{AA}(a_4) nonzero
- Timeline: Weeks 1-8

**Go/No-Go Criterion (Week 8)**:
- GREEN: Fisher metrics match for toy model, inner fluctuation metric is nonzero -> PROCEED to Phase 2
- YELLOW: Partial match, some discrepancies -> investigate, possibly modify the framework
- RED: Fundamental mismatch (e.g., QRE and S_vN give completely different Fisher metrics) -> Paper 9 fails, switch to alternative route

### Phase 2: The Gauge Sector Fisher Metric (Months 3-6) -- Toward Papers 7-8

**Computation 2.1: Fisher metric for SU(2) inner fluctuations**
- Full spectral triple: (C^inf(M^4) tensor (C + H), H, D)
- Inner fluctuations A = sum_alpha a_alpha [D, b_alpha] for a, b in C + H
- The SU(2) gauge field W_mu emerges from A
- Compute g_{W W} (Fisher metric with respect to W_mu)
- Expected: g_{W W} proportional to a_4 coefficient, giving the Yang-Mills action
- This IS Paper 7 (SU(2) from Sigma via inner fluctuations)

**Computation 2.2: Fisher metric for SU(3) inner fluctuations**
- Full spectral triple with A_F = C + H + M_3(C)
- Inner fluctuations generate the gluon field G_mu
- Compute g_{G G} (Fisher metric with respect to G_mu)
- Expected: g_{G G} proportional to a_4, giving SU(3) Yang-Mills
- This IS Paper 8 (SU(3) from Sigma via inner fluctuations)

**Computation 2.3: Verify A_F uniqueness**
- Compute Fisher metrics for alternative algebras (e.g., C + H + M_2(C), or C + M_4(C))
- Check whether the NCG axioms in d=4 uniquely select A_F = C + H + M_3(C)
- If YES: Sigma + d=4 + NCG -> SM gauge group (derivation chain complete)
- If NO: additional principle needed (possibly from QRE positivity or Petz recovery)
- Timeline: Weeks 8-24

**Go/No-Go Criterion (Week 16)**:
- GREEN: SU(2) and SU(3) Fisher metrics reproduce the correct Yang-Mills terms -> Papers 7-8 concrete
- RED: Inner fluctuation Fisher metric does not give gauge kinetic terms -> fundamental problem with the approach

### Phase 3: Fermion Content (Months 6-12) -- Paper 10

**Computation 3.1: Dirac operator eigenvalues on J_3(O)**
- Construct the Dirac operator D_F on the Hilbert space H_F = C^{96} (one generation) or C^{288} (three generations)
- D_F is a 96x96 (or 288x288) matrix whose entries are the Yukawa couplings
- Compute eigenvalues of D_F: these are the fermion masses
- Compare with Singh (2025) mass ratio predictions
- Timeline: Months 6-9

**Computation 3.2: Sigma extremal structure on J_3(O)**
- Define Sigma = D(rho || sigma) for states on the C*-algebra generated by J_3(O)
- Find the extremal points of Sigma under the constraint of fixed total trace
- Check whether the extremal structure has 3-fold structure (corresponding to 3 generations)
- If YES: Paper 10 has a Sigma-based explanation for 3 generations
- If NO: N_gen = 3 remains an algebraic input, not a Sigma prediction
- Timeline: Months 9-12

**Computation 3.3: Fermion mass hierarchy from Sigma**
- If Sigma on J_3(O) has extrema, compute the ratios of Sigma at different extrema
- These ratios may correspond to the fermion mass ratios (m_t/m_b, m_tau/m_mu, etc.)
- Compare with experimental values
- Timeline: Months 10-12

**Go/No-Go Criterion (Month 9)**:
- GREEN: J_3(O) eigenvalues match known fermion mass ratios to within 10% -> Paper 10 viable
- YELLOW: Qualitative agreement (correct ordering, wrong magnitudes) -> needs refinement
- RED: No relation between J_3(O) eigenvalues and fermion masses -> Paper 10 fails

### Phase 4: Grand Synthesis (Months 12-18) -- Paper 11

Combine all results:
- Sigma conformal (a_2) -> gravity + Higgs mass (Papers 2-4, 9)
- Sigma internal SU(2) (a_4) -> weak force (Paper 7)
- Sigma internal SU(3) (a_4) -> strong force (Paper 8)
- Sigma on J_3(O) -> 3 generations + mass ratios (Paper 10)
- Full spectral action = Sigma on M^4 x F = SM + GR

---

## 7. The Specific Computations Needed

### 7.1 Immediate (This Week)

| # | Computation | Input | Expected Output | Difficulty |
|---|-------------|-------|-----------------|------------|
| 1 | S^1 x M_2(C) Gibbs states | Dirac spectrum on S^1 x M_2(C) | S_vN(Q) and Sigma(Q) as functions of Q | Medium |
| 2 | Inner fluctuation on M_2(C) | A = a[D,b] for a,b in M_2(C) | g_{AA} = Fisher metric for gauge direction | Medium |

### 7.2 Short-Term (Months 1-3)

| # | Computation | Input | Expected Output | Difficulty |
|---|-------------|-------|-----------------|------------|
| 3 | A_F = C + M_2(C) spectral action | Full Dirac on S^1 x (C+M_2(C)) | a_2 and a_4 with matter | Hard |
| 4 | CCSvS formula verification | Standard heat kernel | S_vN = spectral action (coefficients) | Medium |
| 5 | QRE vs S_vN for inner fluctuations | Gibbs states with gauge field | Do they agree at Fisher level? | Hard |

### 7.3 Medium-Term (Months 3-6)

| # | Computation | Input | Expected Output | Difficulty |
|---|-------------|-------|-----------------|------------|
| 6 | g_{WW} for SU(2) | Full SM spectral triple (electroweak) | Yang-Mills action from Fisher metric | Very Hard |
| 7 | g_{GG} for SU(3) | Full SM spectral triple | QCD action from Fisher metric | Very Hard |
| 8 | NCG axiom verification for alternative A_F | Various algebras in d=4 | Uniqueness of A_F = C+H+M_3(C) | Hard |

### 7.4 Long-Term (Months 6-12)

| # | Computation | Input | Expected Output | Difficulty |
|---|-------------|-------|-----------------|------------|
| 9 | D_F eigenvalues on J_3(O) | J_3(O) structure constants | Fermion mass spectrum | Very Hard |
| 10 | Sigma extremal on J_3(O) | QRE on J_3(O) state space | 3-fold structure (?) | Research-level |
| 11 | Mass ratio predictions | D_F eigenvalue ratios | Comparison with experiment | Medium (if #9 succeeds) |

---

## 8. Go/No-Go Criteria Summary

### Paper 9 (Sigma = spectral action Fisher metric)

| Criterion | Test | Threshold | Current Status |
|-----------|------|-----------|----------------|
| g_QQ = 2/Q^2 from a_2 | Conformal scaling | Exact match | **PASSED** |
| a_4 conformally invariant | Gauge sector | g_QQ(a_4) = 0 | **PASSED** |
| d = 4 selection | (d-2)(d-3) = 2 | Unique solution | **PASSED** |
| Toy model Fisher match | S^1 x M_2(C) | Agree at O(epsilon^2) | PENDING (Week 8) |
| Inner fluctuation metric | A-direction Fisher | Nonzero | PENDING (Week 8) |
| Full Sigma != 2 ln Q | Q^2-1 vs 2 ln Q | Understand the discrepancy | **KNOWN (partial match)** |

### Papers 7-8 (SU(2), SU(3) from Sigma)

| Criterion | Test | Threshold | Current Status |
|-----------|------|-----------|----------------|
| SU(2) breaking from QRE | D(rho || G_{SU(2)}(rho)) | Correct quantum numbers | **PASSED** (Level 2) |
| SU(3) from octonionic QRE | D(rho || G_{SU(3)}(rho)) | Color quantum numbers | NOT STARTED |
| Yang-Mills from inner Fisher | g_{WW}, g_{GG} | Reproduce S_YM | NOT STARTED |

### Paper 10 (Fermion masses from Sigma)

| Criterion | Test | Threshold | Current Status |
|-----------|------|-----------|----------------|
| J_3(O) eigenvalues exist | Linear algebra | Well-defined | KNOWN (math) |
| Mass ratios within 10% | Compare with PDG | m_t/m_b, m_tau/m_mu, etc. | NOT STARTED |
| 3-fold structure from Sigma | QRE extremal on J_3(O) | 3 extremal directions | NOT STARTED |

---

## 9. Expected Timeline

```
2026 Q2 (Apr-Jun):
  [Week 1-4]  Computation #1-2: S^1 x M_2(C) toy model
  [Week 5-8]  Computation #3-5: A_F = C + M_2(C)
  [Week 8]    GO/NO-GO for Paper 9

2026 Q3 (Jul-Sep):
  [Week 9-16]  Computation #6-7: SU(2) and SU(3) Fisher metrics
  [Week 16]    GO/NO-GO for Papers 7-8
  [Week 17-24] Paper 9 writing + submission

2026 Q4 (Oct-Dec):
  [Week 25-36] Computation #9-10: J_3(O) eigenvalues and Sigma extremal
  [Week 36]    GO/NO-GO for Paper 10

2027 Q1 (Jan-Mar):
  [Week 37-48] Paper 10 writing + Paper 11 planning
  [Week 48]    Complete SM derivation chain (or identify remaining gaps)
```

**Total estimated time: 12 months from today (2026-03-20) to complete Papers 9-10.**

---

## 10. Risk Assessment

### High Risk (>30% failure probability)

1. **QRE vs S_vN mismatch** (~35%): The CCSvS bridge uses von Neumann entropy, while Sigma uses QRE. These are different quantities. The Fisher metrics may match (partial GO) but the full functions may not (preventing a complete identification).

2. **Inner fluctuation Fisher metric does not give Yang-Mills** (~25%): The Fisher metric with respect to gauge parameters A might not reproduce the standard Yang-Mills action. This would block Papers 7-8.

3. **J_3(O) fermion mass ratios are wrong** (~40%): The eigenvalue ratios of D_F on J_3(O) might not match experimental fermion mass ratios. This is the most speculative part.

### Medium Risk (10-30%)

4. **Lorentzian signature issues** (~20%): All spectral action calculations are in Euclidean signature. The Wick rotation to Lorentzian (needed for physical predictions) may introduce subtleties.

5. **Regularization scheme dependence** (~15%): The heat kernel expansion is an asymptotic series. Different regularization schemes may give different results at finite Lambda.

### Low Risk (<10%)

6. **g_QQ = 2/Q^2 is wrong** (~5%): This is a clean mathematical result from well-established Seeley-DeWitt theory. Very unlikely to be wrong.

7. **d = 4 selection fails** (~5%): The equation (d-2)(d-3) = 2 having d=4 as unique solution is a mathematical fact.

---

## 11. What Paper 9 CANNOT Do (Honest Limitations)

1. **Cannot derive A_F from Sigma alone.** The Fisher metric g_QQ = 2/Q^2 is universal (independent of A_F). A_F selection requires the NCG axioms, which are not (yet) derivable from Sigma.

2. **Cannot derive N_gen = 3 from a Fisher metric condition.** No equation analogous to (d-2)(d-3) = 2 selects 3 generations. The best route is through J_3(O) algebraic uniqueness.

3. **Cannot explain confinement.** While "Sigma_color -> infinity for isolated quarks" is a compelling interpretation, it is not a derivation of confinement from the spectral action.

4. **Cannot derive fermion masses from first principles.** The Dirac operator D_F on J_3(O) has eigenvalues, but their relation to physical fermion masses requires the full NCG machinery + renormalization group running.

5. **The full Sigma(Q) is Q^2-1, not 2 ln Q.** These agree at the Fisher (second-order) level but differ globally. The channel result Sigma = 2 ln Q and the spectral action result Sigma ~ Q^2-1 are different functions, suggesting they capture complementary aspects of the same physics.

---

## 12. The Philosophical Point

Paper 9's result reveals a clean DIVISION OF LABOR in the Sigma framework:

- **Gravity (a_2 sector)**: Sigma = 2 ln Q, Fisher metric 2/Q^2. SELECTS d=4. Governs spacetime geometry, dark matter, cosmology.

- **Gauge forces (a_4 sector)**: Conformally invariant, invisible to Sigma under conformal variation. Enters through INTERNAL fluctuations of the Dirac operator. The gauge group is determined by A_F (NCG input), not by Sigma directly.

- **Fermions (H_F)**: Required by the algebraic structure of A_F. The division algebra chain R -> C -> H -> O determines the representations. Sigma constrains their dynamics but does not create them.

This is consistent with Sheng-Kai's philosophy: Sigma does not CREATE new physics from nothing. It provides a unified principle that ORGANIZES known physics by revealing the information-theoretic structure underlying both gravity (conformal sector) and gauge forces (internal sector).

The deep question is whether the NCG axioms -- currently external input -- will eventually be shown to follow from information-theoretic principles. If so, the entire SM would be a consequence of "retrodiction in a quantum world" (Sigma = cost of retrodiction). This is the ultimate goal of Paper 11.

---

*Written: 2026-03-20*
*Status: Research plan for Paper 9 -> SM. Phase 1 toy model computation is the immediate next step.*
*Dependencies: CCSvS (2018), Chamseddine-Connes (1996), Vassilevich (2003), Furey (2025), Singh (2025)*
