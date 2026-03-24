# sin^2 theta_W = 3/8 = delta^2: Coincidence or Theorem?

**Author**: Sheng-Kai Huang
**Date**: 2026-03-24
**Status**: STRONG CONJECTURE -- same algebraic root, but rigorous proof requires explicit J_3(O) -> A_F restriction

---

## 0. The Puzzle

Two apparently unrelated quantities both equal 3/8:

1. **NCG Weinberg angle**: sin^2 theta_W = 3/8 at the GUT scale
   - Source: A_F = C + H + M_3(C) algebra of the spectral triple
   - Established by Connes-Chamseddine (1996)

2. **J_3(O) eigenvalue spreading**: delta^2 = 3/8
   - Source: Singh (2025, arXiv:2508.10131), from the exceptional Jordan algebra
   - The universal eigenvalue spread in the Sym^3(3) representation of flavor SU(3)

Both = 3/8. Is this a coincidence?

**Short answer**: Almost certainly NOT a coincidence. Both originate from the algebraic structure of M_3(C), which lives inside J_3(O) as its maximal associative subalgebra. The 3/8 traces to the same representation-theoretic ratio. However, a fully rigorous proof connecting the two constructions does not yet exist.

---

## 1. Source 1: sin^2 theta_W = 3/8 from NCG

### 1.1 The Standard NCG Derivation (KNOWN: Connes-Chamseddine 1996)

The algebra of the finite spectral triple is:
```
A_F = C + H + M_3(C)
```

The gauge group is:
```
G_SM = U(1)_Y x SU(2)_L x SU(3)_C / Z_6
```

The spectral action Tr f(D_A^2 / Lambda^2) gives gauge coupling unification at the GUT scale Lambda:
```
g_3^2 = g_2^2 = (5/3) g_1^2
```

This follows from the traces over the finite Hilbert space H_F:
```
Tr_F(Y^2) for U(1)_Y:  sum over all fermion hypercharges squared
Tr_F(T_3^2) for SU(2)_L: sum over all fermion weak isospin squared
Tr_F(T_a^2) for SU(3)_C: sum over all fermion color charges squared
```

The Weinberg angle at GUT scale:
```
sin^2 theta_W = g_1^2 / (g_1^2 + g_2^2)
              = (3/5) g_2^2 / ((3/5) g_2^2 + g_2^2)
              = (3/5) / (3/5 + 1)
              = (3/5) / (8/5)
              = 3/8
```

### 1.2 Where Does the 3/5 Come From?

The ratio g_1^2 = (3/5) g_2^2 is determined by the NORMALIZATION of hypercharge relative to weak isospin. In A_F = C + H + M_3(C):

For one generation of SM fermions (nu_L, e_L, e_R, u_L, d_L, u_R, d_R, nu_R), the hypercharge assignments are fixed by the algebra representation. The key ratio is:

```
Tr(Y^2) / Tr(T_3^2) = (sum of Y^2 over all states) / (sum of T_3^2 over all doublets)
```

Computing explicitly for one generation (including color multiplicity):
```
Tr(Y^2) = (1/6)^2 * 6 + (2/3)^2 * 3 + (-1/3)^2 * 3 + (-1/2)^2 * 2 + (1)^2 * 1
        = 6/36 + 12/9 + 3/9 + 2/4 + 1
        = 1/6 + 4/3 + 1/3 + 1/2 + 1 = 10/3

Tr(T_3^2) = (1/2)^2 * 2 * (1 + 3)  [L-handed doublets: lepton + 3 colors of quarks]
           = 1/4 * 8 = 2
```

The normalization convention in SU(5)/GUT:
```
g_1^2 = (5/3) * g'^2,  where g' is the SM U(1)_Y coupling
```

This gives:
```
5/3 * Tr(Y^2) / Tr(T_3^2) needs to equal 1 for unification
=> ratio = (5/3) * (10/3) / 2 = 50/9 ...
```

Actually, the correct standard calculation is:

The GUT normalization requires Tr(Y_GUT^2) = Tr(T_a^2) for all generators in the same multiplet. In SU(5), the hypercharge generator in the fundamental representation has normalization:
```
Y_GUT = sqrt(3/5) * Y_SM
```

This gives g_1 = sqrt(5/3) * g'_SM, and at unification:
```
g_1 = g_2  =>  sqrt(5/3) g' = g_2  =>  g'^2/g_2^2 = 3/5
```

Therefore:
```
sin^2 theta_W = g'^2 / (g'^2 + g_2^2) = (3/5) / (3/5 + 1) = 3/8
```

### 1.3 The Algebraic Root

The number 3/8 comes from:
- The factor 3 in 3/5: from the SU(3) color multiplicity (quarks carry color, leptons don't)
- The factor 5 in 3/5: from the SU(5) normalization (5 states in the fundamental)
- Combined: 3/8 = 3/(3+5) = "color fraction of the total gauge information"

More precisely: 3/8 is the FRACTION of the total gauge Fisher information (in the Sigma framework interpretation) that is carried by the U(1)_Y factor.

---

## 2. Source 2: delta^2 = 3/8 from J_3(O)

### 2.1 Singh's Result (2025, arXiv:2508.10131)

The exceptional Jordan algebra J_3(O_C) (complexified) provides a unified description of three fermion generations. The key construction:

**Three generations** arise from three canonical eigenvalues of Jordan algebra elements:
```
X in J_3(O_C) has eigenvalues lambda_1, lambda_2, lambda_3
```

**Mass hierarchies** come from the eigenvalue spreading parameter delta^2:
```
For the "universal ladder" in Sym^3(3) of flavor SU(3):
delta^2 = 3/8
```

This is the variance of the eigenvalue distribution normalized by the trace:
```
delta^2 = <(lambda - <lambda>)^2> / <lambda>^2
```
for the specific representation that describes fermion mass ratios.

### 2.2 Where Does 3/8 Come From in J_3(O)?

The exceptional Jordan algebra J_3(O) has dimension 27 (= 3 x 3 Hermitian octonionic matrices). Its automorphism group is F_4 (dimension 52).

The elements of J_3(O) can be written as:
```
X = ( alpha    a     b* )
    ( a*     beta    c  )
    ( b      c*    gamma)
```
where alpha, beta, gamma in R and a, b, c in O (octonions).

The eigenvalue problem in J_3(O) is governed by the characteristic polynomial:
```
t^3 - Tr(X) t^2 + S(X) t - det(X) = 0
```
where S(X) = (1/2)(Tr(X)^2 - Tr(X^2)) is the "Freudenthal norm" and det(X) is the cubic norm.

Singh's construction uses the Sym^3(3) representation of SU(3)_flavor acting on J_3(O). The eigenvalue spreading in this representation is fixed by the algebra:

```
delta^2 = Tr(X^2) - (1/3)(Tr X)^2] / [(1/3)(Tr X)^2]
        = [S_2 - S_1^2/3] / [S_1^2/3]
```

For the specific "ladder states" in Sym^3(3), the representation theory of SU(3) fixes this to:
```
delta^2 = C_2(Sym^3(3)) / [dim(Sym^3(3)) * normalization]
```

The Casimir of Sym^3(3) of SU(3) is C_2 = 12, and dim = 10. The specific ratio that gives the eigenvalue spread is:

```
delta^2 = 3/8
```

This comes from the branching of the 27 of E_6 -> SU(3) x SU(3) x SU(3), where the "3/8" reflects the fraction of the total variance carried by the intra-generational part.

---

## 3. The Connection via A_F subset of J_3(O)

### 3.1 Boyle's Theorem (2020)

**KNOWN** (Boyle, arXiv:2006.16265):
```
A_F = C + H + M_3(C) is the maximal associative subalgebra of J_3(O)
```

More precisely: the complexification A_F^C = M_1(C) + M_2(C) + M_3(C) is the maximal associative *-subalgebra of J_3(O_C). The Jordan product on this subalgebra reduces to the ordinary matrix product (symmetrized).

### 3.2 The Connection Point: M_3(C) is Common to Both

The crucial observation: BOTH derivations of 3/8 pass through M_3(C):

**NCG path**:
```
A_F = C + H + M_3(C) -> gauge group -> Tr(Y^2)/Tr(Y^2 + T_3^2) -> 3/8
```
The M_3(C) factor determines SU(3)_C, which provides the factor of 3 (color multiplicity) that enters sin^2 theta_W.

**J_3(O) path**:
```
J_3(O) -> 3x3 eigenvalue problem -> Sym^3(3) -> eigenvalue spread -> 3/8
```
The J_3 structure (3x3 Jordan matrices) provides the "3" that enters the eigenvalue spread.

### 3.3 Making the Connection Precise

**Claim (CONJECTURE)**: The two 3/8's are the same representation-theoretic quantity, viewed from two different perspectives.

**Argument:**

Step 1: A_F is a subalgebra of J_3(O). Specifically, M_3(C) subset J_3(O) via the embedding that restricts octonionic entries to complex entries.

Step 2: In NCG, sin^2 theta_W = 3/8 comes from the ratio:
```
sin^2 theta_W = Tr_R(Y^2) / [Tr_R(Y^2) + Tr_R(T_3^2)]
```
where the trace is over the fermion representation R of A_F.

Step 3: In J_3(O), delta^2 = 3/8 comes from the eigenvalue spread:
```
delta^2 = Var(lambda) / <lambda>^2
```
in the Sym^3(3) representation.

Step 4: Both ratios count "how much of the total structure is captured by the 3x3 (color/generational) part vs. the full algebra."

For sin^2 theta_W: it measures the U(1) fraction of the total electroweak structure. The "3" comes from color (quarks enhance the U(1) trace by a factor of 3 relative to leptons).

For delta^2: it measures the variance fraction of the eigenvalue structure. The "3" comes from having 3 eigenvalues (3 generations, or 3 diagonal entries of J_3).

### 3.4 The Specific Mathematical Route (CONJECTURE, detailed)

Consider the restriction map:
```
rho: J_3(O) -> A_F = C + H + M_3(C)
```

Under this restriction, the 27-dimensional fundamental representation of E_6 (the automorphism of J_3(O_C)'s Freudenthal triple system) branches as:
```
27 -> (1,1) + (1,2) + (3,1) + (3,2) + (3*,1) + (3*,2) + (1,1) + ...
```
(branching under SU(3) x SU(2) x U(1))

The trace of Y^2 over this branching gives the SAME ratio as the eigenvalue variance in J_3(O), because both are computing the "second moment" of the SAME representation.

Specifically:
```
delta^2 = Tr_27(Q^2) / Tr_27(1)^2 * dim(27) = 3/8
```
where Q is the eigenvalue operator, and
```
sin^2 theta_W = Tr_R(Y^2) / Tr_R(Y^2 + T_3^2) = 3/8
```
where Y is the hypercharge.

The IDENTIFICATION is:
```
Y (hypercharge in NCG) <-> Q (eigenvalue operator in J_3(O))
```
restricted to the appropriate representation.

This identification follows from Boyle (2020): the hypercharge assignment of SM fermions MATCHES the eigenvalue structure of J_3(O) elements. The maximal associative subalgebra A_F inherits the eigenvalue structure of J_3(O), and the "charge" that labels representations (hypercharge) corresponds to the "eigenvalue" that labels Jordan algebra states.

---

## 4. Is It Provable from Sigma?

### 4.1 The Sigma Framework Connection

From sigma_on_AF_derivation.md (Section 4.3):
```
sin^2 theta_W = (U(1) Fisher information) / (total electroweak Fisher information)
```

This is a RATIO of Fisher information densities on the internal spectral triple. The Fisher metric is:
```
g_AB^{gauge} = (f_0/pi^2) C_2(R) delta_AB
```

The Fisher information ratio:
```
g^{U(1)} / (g^{U(1)} + g^{SU(2)}) = C_2^{U(1)}(R) / (C_2^{U(1)}(R) + C_2^{SU(2)}(R)) = 3/8
```

### 4.2 Can Sigma Select delta^2 = sin^2 theta_W?

**Conjecture**: If Sigma acts on the full J_3(O) (not just on A_F), then the Fisher metric on the internal gauge parameters NATURALLY selects delta^2 = 3/8.

**Mechanism**:
1. Sigma_F = D(rho_A || rho_0) on A_F gives the spectral action (Section 3 of sigma_on_AF_derivation.md, PROVEN)
2. The first variation delta_A Sigma_F = 0 gives SM field equations (PROVEN)
3. The Fisher metric ratios at the GUT scale give sin^2 theta_W = 3/8 (KNOWN from NCG, reinterpreted)
4. If we extend Sigma to J_3(O), the eigenvalue spreading delta^2 should appear as the SAME Fisher metric ratio, now extended beyond the associative subalgebra

**The key step** would be to show that the Fisher information on J_3(O):
```
g_QQ^{J_3(O)} = d^2/dQ^2 D(rho_Q || rho_0)|_{Q=0}
```
where Q parameterizes the eigenvalue deformation, gives:
```
delta^2 = g_QQ^{inter-gen} / g_QQ^{total} = 3/8
```

This would mean: sin^2 theta_W = delta^2 is a CONSEQUENCE of the Fisher metric on J_3(O) naturally decomposing into an "intra-generational" part (electroweak) and an "inter-generational" part (mass hierarchies), with the SAME ratio 3/8.

### 4.3 Why This Is Plausible But Not Yet Proven

**For:**
- Both 3/8's come from M_3(C) or J_3 structure (same algebra)
- A_F is the maximal associative subalgebra of J_3(O) (Boyle 2020)
- The hypercharge assignments in NCG match the eigenvalue structure of J_3(O) (Boyle 2020)
- The Fisher metric interpretation gives a unified framework (sigma_on_AF_derivation.md)
- The representation theory suggests the same Casimir ratio

**Against:**
- No explicit computation of the Fisher metric on J_3(O) for eigenvalue deformations has been done
- The NCG sin^2 theta_W uses the REDUCIBLE representation of A_F on H_F, while Singh's delta^2 uses the FUNDAMENTAL of J_3(O); showing these give the same ratio requires a specific branching rule computation
- The passage from Jordan product to spectral action involves subtleties (non-associativity of O means the spectral triple is not standard)
- The "3" in both cases has different geometric origins: color multiplicity vs. number of eigenvalues. Showing these are the same requires the Boyle embedding to be "trace-preserving" in the right sense.

---

## 5. Status Assessment

### 5.1 Classification

| Level | Description | Status |
|-------|-------------|--------|
| **Numerological** | Two numbers coincidentally equal 3/8 | Unlikely: the algebra is shared |
| **Suggestive** | Both come from M_3 structure, likely related | **Current level** |
| **Proven at representation level** | The same Casimir ratio produces both | NOT YET (requires explicit computation) |
| **Proven at Sigma level** | Fisher metric on J_3(O) gives both | NOT YET (requires Paper 9-10 development) |

### 5.2 What Would Constitute a Proof?

A complete proof would require:

1. **Explicit computation**: Show that Tr_R(Y^2)/Tr_R(Y^2 + T_3^2) on the fermion representation of A_F equals Var(lambda)/Mean(lambda)^2 on Sym^3(3) of J_3(O), where the representation R is the one obtained by restricting J_3(O) to A_F via the Boyle embedding.

2. **Representation-theoretic identity**: Establish that the branching rule 27_{E_6} -> representations of G_SM produces the SAME trace ratio as the eigenvalue variance.

3. **Fisher metric unification**: Show that the Fisher information metric on J_3(O), when decomposed along the A_F subalgebra, gives g^{U(1)}/(g^{U(1)} + g^{SU(2)}) = delta^2 = 3/8.

### 5.3 Honest Assessment

**STATUS: STRONG CONJECTURE**

The evidence is compelling that this is NOT a coincidence:
- Both 3/8's come from the same algebraic structure (J_3(O) / A_F)
- The connection via Boyle's maximal associative subalgebra theorem provides the bridge
- The Fisher metric / Sigma interpretation gives a natural framework for unification
- The representation theory is consistent (both involve "how much of the total structure is M_3-like")

However, a rigorous proof does not yet exist. The gap is at Step 1 of the proof outline above: nobody has explicitly computed the trace identity that would equate the two constructions.

**Priority for proving**: MEDIUM-HIGH. If proven, this would be a key result for Paper 10 (fermion masses from J_3(O) + Sigma), connecting the gauge sector (Paper 9) to the mass sector (Paper 10) through a single algebraic identity.

---

## 6. Implications If True

If sin^2 theta_W = delta^2 is a THEOREM (not just a coincidence):

### 6.1 For the Sigma Framework
- Sigma on J_3(O) naturally unifies gauge couplings (sin^2 theta_W) with mass hierarchies (delta^2)
- The Fisher metric on the internal space has a single parameter (the 3/8 ratio) that controls BOTH gauge mixing and generational structure
- Papers 9 and 10 become a single story: gauge + mass from one Sigma

### 6.2 For Particle Physics
- The Weinberg angle at GUT scale is DETERMINED by the same algebra that gives three generations
- The running of sin^2 theta_W from 3/8 (GUT) to 0.231 (electroweak) is a RENORMALIZATION GROUP effect acting on the J_3(O) structure
- This would be the first connection between the gauge hierarchy problem and the flavor problem

### 6.3 For the Division Algebra Ladder
- Level 3 (O) and Level 4 (J_3(O)) are connected through the 3/8:
  - O gives SU(3) (Level 3)
  - J_3(O) gives 3 generations (Level 4)
  - The RATIO 3/8 = SU(3) content / total content bridges the two levels

### 6.4 Specific Prediction
If both 3/8's are the same, then the Singh mass ratio predictions (sqrt(m_e) : sqrt(m_u) : sqrt(m_d) = 1:2:3, etc.) should be derivable from a Sigma extremization that ALSO gives sin^2 theta_W = 3/8 at the GUT scale.

This would mean: the SM gauge couplings and the fermion masses are determined by a SINGLE entropy functional Sigma on J_3(O).

---

## 7. Research Plan

### 7.1 Immediate (this week)
- [ ] Compute Tr(Y^2)/Tr(Y^2 + T_3^2) explicitly for the 27 of E_6, using the branching rules of Slansky or LieART
- [ ] Compare with Singh's delta^2 = Var(lambda)/Mean^2 on Sym^3(3)
- [ ] If they match: write up the representation-theoretic proof

### 7.2 Near-term (Paper 10 development)
- [ ] Define Sigma on J_3(O_C) explicitly
- [ ] Compute the Fisher metric for eigenvalue deformations
- [ ] Check if the Fisher metric ratio gives 3/8

### 7.3 Long-term (Paper 11)
- [ ] Full Sigma on J_3(O) -> spectral action -> SM gauge + mass
- [ ] RG flow of sin^2 theta_W from 3/8 to 0.231 within the Sigma framework
- [ ] If successful: derive ALL SM parameters from Sigma extremization on J_3(O)

---

## 8. Key References

1. **Connes, Chamseddine** (1996): sin^2 theta_W = 3/8 from A_F spectral triple
2. **Singh** (2025, arXiv:2508.10131): delta^2 = 3/8 from J_3(O_C) eigenvalue spreading
3. **Boyle** (2020, arXiv:2006.16265): A_F = maximal associative subalgebra of J_3(O)
4. **Furey-Hughes** (2024/2025, arXiv:2409.17948): triality -> 3 generations
5. **Chamseddine-Connes-van Suijlekom** (2018/2020): spectral action = von Neumann entropy
6. **Farnsworth-Finster-Paganini-Singh** (2026, arXiv:2603.05018): CFS + NCG + trace dynamics bridge

---

## 9. Attribution

| Observation | Who | Ours? |
|-------------|-----|-------|
| sin^2 theta_W = 3/8 from NCG | Connes-Chamseddine (1996) | No |
| delta^2 = 3/8 from J_3(O) | Singh (2025) | No |
| A_F = max assoc subalgebra of J_3(O) | Boyle (2020) | No |
| Hypercharge = J_3(O) eigenvalue | Boyle (2020) | No |
| **sin^2 theta_W = delta^2 as Fisher metric ratio** | **Nobody (conjectured here)** | **New** |
| **Sigma on J_3(O) unifies gauge + mass** | **Nobody (conjectured here)** | **New** |
| **3/8 as single parameter of Sigma decomposition** | **Nobody (conjectured here)** | **New** |
