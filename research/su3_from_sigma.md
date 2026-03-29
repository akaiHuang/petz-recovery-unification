# SU(3) from Sigma: The Octonionic Challenge
## 2026-03-19
## Sheng-Kai Huang

---

## Executive Summary

This document investigates whether the Sigma framework (Sigma = D(rho_spacetime || rho_matter) = 2 ln Q) can derive SU(3) gauge symmetry (QCD / strong force). Five major routes are examined: (1) octonionic extension of the division algebra ladder, (2) Kaluza-Klein on S^7, (3) Connes NCG via A_F = C + H + M_3(C), (4) Farnsworth's nonassociative spectral geometry, and (5) confinement as Petz recovery failure.

**Key conclusions**:
1. The division algebra ladder at Level 3 (O -> SU(3)) FAILS at the algebraic level: QRE on 3-qubit systems under SU(2)_diagonal gives dim V_{3/2} = 4, not dim Im(O) = 7. The transition from SU(2) to SU(3) requires G_2 (the automorphism group of O), which does not arise from any n-qubit tensor product structure.
2. The Kaluza-Klein approach on S^7 gives Sigma_{11D} = Sigma_grav + Sigma_{SU(3)} + Sigma_{SU(2)} + ..., which is geometrically clean but does NOT derive SU(3) from information theory -- the gauge group is a geometric input from the choice of internal manifold.
3. The Connes NCG approach (A_F = C + H + M_3(C)) is the MOST VIABLE route. SU(3) emerges from the M_3(C) factor. If Paper 9 succeeds in connecting Sigma to the Fisher information of the spectral action, SU(3) comes for free from the algebra A_F.
4. Farnsworth's nonassociative spectral geometry (2025) is the MOST MATHEMATICALLY NOVEL approach but is in early development. If extended from G_2 x G_2 to the full SM, this could provide the first direct SU(3)-from-octonions derivation compatible with Sigma.
5. Confinement as Petz recovery failure (Sigma_color -> infinity for isolated quarks) is a compelling INTERPRETIVE framework but not a derivation. It connects naturally to area-law entanglement and the entanglement entropy of the QCD vacuum.

**Recommended strategy**: Pursue the Connes NCG route (Paper 9) as the primary path. This gives SU(3) as part of the full SM package. Treat the octonionic and confinement perspectives as complementary interpretive layers, not alternative derivations.

---

## 1. The Problem: Why SU(3) Is Fundamentally Harder Than SU(2)

### 1.1 Recap: What worked for U(1) and SU(2)

| Level | Algebra | Gauge Group | Method | Status |
|-------|---------|-------------|--------|--------|
| 1 | C | U(1) | 5 independent routes (Paper 6 Thm 1) | DONE |
| 2 | H | SU(2) | QRE chain rule, D_break = ln 3 = dim Im(H) | Level 2 DONE |

For U(1): the complex structure of Hilbert space (required by Sigma) has a natural U(1) phase symmetry. Gauging it gives electromagnetism.

For SU(2): the quaternionic structure of 2-qubit systems gives a natural SU(2) symmetry. The QRE chain rule under the SU(2) twirl decomposes D into 3 independent directions matching Im(H) = {i, j, k}.

### 1.2 Why the same approach fails for SU(3)

**Obstacle 1: Non-associativity.**
The octonions O are the natural division algebra at Level 3. But O is non-associative: (ab)c != a(bc). This means:
- No matrix representation of O exists (matrices are associative)
- O cannot form a C*-algebra (C*-algebras are associative by definition)
- QRE = D(rho || sigma) = Tr(rho ln rho - rho ln sigma) uses matrix operations, which require associativity
- Therefore, QRE CANNOT be directly defined on octonionic "density operators"

This is a structural impossibility, not a technical difficulty.

**Obstacle 2: No tensor product.**
For U(1): 1-qubit Hilbert space C^2 has natural U(1) structure.
For SU(2): 2-qubit Hilbert space C^2 tensor C^2 = C^4 has natural SU(2) structure.
For SU(3): 3-qubit Hilbert space C^2 tensor C^2 tensor C^2 = C^8 has SU(2)^3 x SU(2)_diagonal symmetry, NOT SU(3).

The 3-qubit system has 8 = 2^3 complex dimensions. The octonions also have 8 real dimensions. But:
- C^8 decomposes under SU(2)_diagonal as: j=3/2 (dim 4) + j=1/2 (dim 2) + j=1/2 (dim 2)
- The maximum irrep has dim 4, not 7 = dim Im(O)
- There is no natural SU(3) action on C^8 that comes from the tensor product structure

This was verified on 2026-03-19: Level 3 FAILS algebraically. The number 4 != 7.

**Obstacle 3: G_2, not SU(3), is the natural automorphism group.**
Aut(O) = G_2, the exceptional Lie group of dimension 14. SU(3) sits inside G_2 as a subgroup:

```
G_2 supset SU(3)    [index 2, dim 14 supset dim 8]
```

The decomposition: G_2 acts on Im(O) = R^7. Under the subgroup SU(3):
```
7 = 1 + 3 + 3_bar    (singlet + fundamental + anti-fundamental)
```

So SU(3) does not act on ALL 7 imaginary units -- it acts on a 6-dimensional subspace (3 + 3_bar), with one singlet direction.

The passage from G_2 to SU(3) requires CHOOSING a preferred imaginary unit (say e_7), which defines the embedding SU(3) = Stab_{G_2}(e_7). This choice is NOT natural in the Sigma framework -- it requires additional structure.

### 1.3 The octonionic Fano plane

The multiplication rules of octonions are encoded in the Fano plane (the projective plane over F_2):

```
        e_1
       / | \
      /  |  \
    e_2  |  e_4
    / \  |  / \
   /   \ | /   \
  e_3   e_6   e_5
   \   / | \   /
    \ /  |  \ /
    e_7--+--e_7
```

Each line (including the inscribed circle) represents a quaternionic triple: e_a * e_b = e_c (with appropriate signs). There are 7 lines, each giving a copy of Im(H) = {e_a, e_b, e_c}.

The Fano plane has PSL(2,7) = GL(3, F_2) as its automorphism group (168 elements). The continuous automorphism group is G_2 (14-dimensional). The SU(3) subgroup preserves one line (= one quaternionic subspace).

**Key insight for Sigma**: The QRE chain rule naturally decomposes along ONE quaternionic triple at a time (giving SU(2) chain rule). To get the FULL G_2 structure, one would need a decomposition that respects ALL 7 lines simultaneously. This is the "octonionic chain rule" that does not yet exist.

---

## 2. Route 1: Octonionic Extension of Division Algebra Ladder

### 2.1 The Szangolies hierarchy

Szangolies (2025, arXiv:2512.17328) proposed:

```
1-qubit entanglement -> Hopf S^3 -> S^2  -> fiber S^1 = U(1)
2-qubit entanglement -> Hopf S^7 -> S^4  -> fiber S^3 = SU(2)
3-qubit entanglement -> Hopf S^15 -> S^8 -> fiber S^7 => G_2 supset SU(3)
```

The fiber dimensions 1, 3, 7 match dim Im(C), dim Im(H), dim Im(O). This is the topological version of the division algebra ladder.

### 2.2 Why the algebraic version fails at Level 3

The algebraic QRE chain rule under the G twirl gives:
```
D(rho || I/d) = D(rho || G_twirl(rho)) + D(G_twirl(rho) || I/d)
```

For SU(2) twirl on 2 qubits: D_break = ln(2j+1) where j is the maximal irrep.
- j = 1 (triplet): D_break = ln 3 matches dim Im(H) = 3

For SU(2)_diagonal twirl on 3 qubits: the maximal symmetric irrep has j = 3/2, dim = 4.
- D_break = ln 4, NOT ln 7 or ln 8

For G_2 twirl on 3 qubits: G_2 does not act naturally on (C^2)^{tensor 3}. The fundamental representation of G_2 is 7-dimensional (real), while the 3-qubit space is 8-dimensional (complex). There is no natural embedding.

**The algebraic and topological pictures diverge at Level 3.** The Hopf fibration S^15 -> S^8 lives in the TOPOLOGICAL structure of 3-qubit entanglement, not in the algebraic structure of QRE decomposition.

### 2.3 What would be needed: octonionic QRE

A hypothetical "octonionic QRE" would require:
1. Density operators over O^n (n-tuples of octonions)
2. A trace operation Tr_O compatible with non-associativity
3. A logarithm operation compatible with non-associativity
4. D_O(rho || sigma) = Tr_O(rho *_O ln_O(rho) - rho *_O ln_O(sigma))

None of these exist in the standard mathematical literature. The non-associativity of O prevents defining matrix operations in the usual sense.

**Partial progress**: The exceptional Jordan algebra J_3(O) consists of 3x3 Hermitian matrices over O:

```
       / a_1   o_3*  o_2  \
X =   | o_3   a_2   o_1* |    a_i in R, o_i in O
       \ o_2*  o_1   a_3  /
```

with the Jordan product X o Y = (XY + YX)/2 (which IS well-defined even though O is non-associative, because the product only involves 2 elements at a time -- the Albert algebra is power-associative).

J_3(O) has dimension 27 = 3 * 1 + 3 * 8 = 3 real entries + 3 octonionic entries.

**Crucial fact**: The automorphism group of J_3(O) is F_4 (the exceptional group of dimension 52). F_4 contains G_2 (which preserves each octonionic entry) and also SU(3) (which permutes the three octonionic entries). In fact:

```
F_4 supset SU(3) x G_2
```

and the Standard Model gauge group G_SM = SU(3) x SU(2) x U(1) / Z_6 sits inside F_4.

### 2.4 The Boyle construction (arXiv:2006.16265)

Boyle (2020) showed that the FULL Standard Model structure can be extracted from the algebra C tensor H tensor O (or equivalently from J_3(O)):

1. One generation of SM fermions = a minimal left ideal of C tensor H tensor O
2. The gauge group G_SM = intersection of Aut(C tensor H tensor O) with the unbroken symmetries
3. Three generations arise from the 3 eigenvalues of a J_3(O) element (SO(8) triality)
4. Chirality arises from the Z_2 grading of the algebra

Boyle's construction gives the CORRECT particle content and gauge group, but does NOT give dynamics. It is a kinematic identification.

### 2.5 Assessment: Route 1

| Aspect | Status |
|--------|--------|
| Topological match (Hopf S^15 -> S^8, fiber S^7) | Correct |
| Algebraic QRE decomposition | FAILS (4 != 7) |
| Octonionic QRE | Does not exist (non-associativity) |
| J_3(O) Jordan algebra approach | Well-defined, F_4 supset G_SM |
| Dynamics from Sigma | NOT derived |

**Verdict**: The octonionic route identifies SU(3) through the algebraic structure of O and J_3(O), but CANNOT derive SU(3) from the QRE chain rule. The non-associativity of O is the fundamental obstruction. This route provides the correct GROUP but not the correct DYNAMICS.

---

## 3. Route 2: Kaluza-Klein on S^7 and 11D Supergravity

### 3.1 The geometric setup

The maximal KK construction for the Standard Model gauge group requires:

```
11D = 4D + 7D internal space
```

where the internal 7-manifold has isometry group containing G_SM = SU(3) x SU(2) x U(1).

The natural candidates for the internal space:
- S^7 (7-sphere): isometry group SO(8), which contains G_2, SU(3), SU(2), U(1)
- CP^2 x S^2 x S^1: isometry group SU(3) x SU(2) x U(1) (exactly G_SM!)
- Squashed S^7: can have isometry reduced to SU(3) x SU(2) x U(1)

The connection to octonions: S^7 = unit octonions, and:

```
S^7 = Spin(7) / G_2
G_2 / SU(3) = S^6
SU(3) / SU(2) = S^5
SU(2) / U(1) = S^3
```

So the Standard Model gauge groups form a natural chain of subgroups within the isometries of S^7.

### 3.2 Sigma in 11D decomposition

Extending the Paper 2 and Paper 6 KK decompositions, the 11D entropy production is:

```
Sigma_{11D} = -ln(-G_{00}^{(11)})
```

The 11D metric decomposes as:

```
ds_{11}^2 = g_{mu nu}(x) dx^mu dx^nu
           + r_1^2 [d Omega_5^2(CP^2) + 2 B^alpha_mu dx^mu (CP^2 Killing 1-forms)_alpha + ...]
           + r_2^2 [d Omega_2^2(S^2) + 2 A^a_mu dx^mu (S^2 Killing 1-forms)_a + ...]
           + r_3^2 [d theta^2 + 2 C_mu dx^mu d theta + ...]
```

where:
- CP^2 Killing vectors generate SU(3) (8 generators, but only gauge 8 bosons)
- S^2 Killing vectors generate SU(2) (3 generators)
- S^1 Killing vector generates U(1) (1 generator)

The temporal metric component picks up gauge field contributions:

```
G_{00}^{(11)} = g_{00} + r_1^2 |B_0|^2_{SU(3)} + r_2^2 |A_0|^2_{SU(2)} + r_3^2 |C_0|^2_{U(1)}
```

Therefore:

```
Sigma_{11D} = -ln(-g_{00} - r_1^2 |B_0|^2 - r_2^2 |A_0|^2 - r_3^2 |C_0|^2)
```

In the weak-field regime:

```
Sigma_{11D} = Sigma_grav + Sigma_{SU(3)} + Sigma_{SU(2)} + Sigma_{U(1)}
```

where:

```
Sigma_{SU(3)} = -ln(1 - r_1^2 |B_0|^2 / |g_{00}|)
              ~ r_1^2 |B_0|^2 / |g_{00}|   (weak field)
              ~ (1/2g_3^2) B^alpha_{0 mu} B^{alpha,0 mu}  (after KK reduction)
```

This is exactly the temporal component of the SU(3) Yang-Mills Lagrangian.

### 3.3 Connection to M-theory and 11D supergravity

The 11D KK construction is precisely the starting point of M-theory / 11D supergravity (Cremmer-Julia-Scherk 1978). In that context:

- The 11D theory has a 3-form potential C_{MNP} (the M-theory 3-form)
- Compactification on a 7-manifold Y_7 gives a 4D theory with gauge group = isometry group of Y_7
- For Y_7 = S^7 (round): gauge group is SO(8), enhanced to the full SU(8) supergravity
- For Y_7 = squashed S^7: gauge group reduces to SU(3) x SU(2) x U(1) (the Witten-Awada-Duff-Pope solution)

**The Sigma framework adds**: the entropy decomposition Sigma_{11D} = sum of gauge contributions, where each gauge sector contributes proportional to its temporal field strength. This is a NEW observation: the standard KK decomposition gives the Lagrangian, but the Sigma decomposition gives the information-theoretic cost of each gauge interaction.

### 3.4 Does this DERIVE SU(3)?

**No.** The KK approach assumes:
1. The dimension is 11 (not derived)
2. The internal manifold is S^7 or CP^2 x S^2 x S^1 (not derived)
3. The metric ansatz factorizes as 4D + 7D (an approximation)

SU(3) enters as the isometry group of the internal space, which is a GEOMETRIC INPUT. The Sigma framework provides the entropy decomposition once the geometry is given, but does not select the geometry.

**What Sigma COULD contribute**: If the variational principle delta Sigma_{11D} = 0 selects the internal manifold, then SU(3) would be derived. Specifically:

Conjecture: Among all 7-manifolds Y_7, the one that extremizes Sigma_{11D} (retrodictability) is the one whose isometry group contains G_SM.

This is a strong conjecture with no proof, but it would give a selection principle for the gauge group.

### 3.5 The gauge couplings from 11D

The three gauge couplings in the KK framework are:

```
g_3^2 = 16 pi G_{11} / (Vol(CP^2) * r_1^5)    -- SU(3) from CP^2
g_2^2 = 16 pi G_{11} / (Vol(S^2) * r_2^2)      -- SU(2) from S^2
g_1^2 = 16 pi G_{11} / (2 pi r_3)               -- U(1) from S^1
```

The Weinberg angle and coupling unification conditions relate the radii r_1, r_2, r_3. At the GUT scale, if all radii are equal (maximum symmetry):

```
sin^2(theta_W) = 3/8    (same as Connes NCG prediction)
alpha_3 = alpha_2        (SU(3) and SU(2) unify)
```

### 3.6 Assessment: Route 2

| Aspect | Status |
|--------|--------|
| 11D -> 4D + S^7 gives SU(3) | Standard KK (known) |
| Sigma_{11D} decomposition | NEW (explicit formula in Section 3.2) |
| SU(3) derived from Sigma | NO (geometric input) |
| Connection to M-theory | Standard |
| Gauge couplings | From geometry (r_1, r_2, r_3), not from Sigma |

**Verdict**: The KK route provides a clean geometric framework and a NEW entropy decomposition formula, but does not derive SU(3) from information theory. It is a "compatibility check," not a derivation.

---

## 4. Route 3: Connes NCG and the Spectral Action

### 4.1 The finite algebra A_F = C + H + M_3(C)

In Connes' noncommutative geometry (NCG), the Standard Model is described by an "almost-commutative" spectral triple:

```
(A, H, D) = (C^inf(M) tensor A_F, L^2(S) tensor H_F, D_M tensor 1 + gamma_5 tensor D_F)
```

where the finite algebra is:

```
A_F = C + H + M_3(C)
```

This algebra is associative (critical!) and the gauge group is:

```
G = U(A_F) / U(center(A_F)) = U(1) x SU(2) x SU(3)
```

The three factors come from:
- C: contributes U(1) (hypercharge)
- H: contributes SU(2) (weak isospin) -- unit quaternions are SU(2)
- M_3(C): contributes SU(3) (color) -- unitary 3x3 complex matrices

**SU(3) comes from M_3(C)**, the algebra of 3x3 complex matrices. This is a MATRIX algebra, associative and well-behaved.

### 4.2 Why M_3(C)?

Connes and Chamseddine showed that A_F = C + H + M_3(C) is essentially UNIQUE given the following axioms:
1. The spectral triple is "almost-commutative" (M x F with F finite)
2. The KO-dimension of F is 6 (mod 8) -- required for chirality
3. The first-order condition for the Dirac operator
4. Massless photon condition (unbroken U(1)_EM after symmetry breaking)
5. No more than one generation of fermions in the finite space

Under these conditions, the algebra is uniquely A_F = C + H + M_3(C) (or mild variants like M_2(H) + M_4(C) which reduce to the same physics).

**For the Sigma framework**: The question "why SU(3)?" reduces to "why M_3(C) in A_F?" which in turn reduces to "why these 5 axioms?" The Sigma framework could potentially motivate the axioms -- particularly the KO-dimension condition -- but this has not been done.

### 4.3 Connection to Sigma via spectral action

The CCSvS (2018) result states:

```
S_vN(rho_beta(D^2)) = Spectral Action = Tr f(D^2 / Lambda^2)
```

In the almost-commutative case M x F, the spectral action expands as:

```
S_spectral = integral d^4x sqrt(g) [
    (1/16 pi G) R                                    -- Einstein-Hilbert (from C^inf(M))
  + (1/4 g_1^2) B_{mu nu} B^{mu nu}                 -- U(1) hypercharge (from C)
  + (1/4 g_2^2) W^a_{mu nu} W^{a, mu nu}            -- SU(2) weak (from H)
  + (1/4 g_3^2) G^alpha_{mu nu} G^{alpha, mu nu}    -- SU(3) color (from M_3(C))
  + |D_mu H|^2 + mu_H^2 |H|^2 + lambda |H|^4       -- Higgs (from inner fluctuations)
  + Yukawa terms                                      -- fermion masses (from D_F)
]
```

**The SU(3) Yang-Mills term arises automatically from the M_3(C) factor in A_F.**

If Paper 9 succeeds in showing Sigma = Fisher information of spectral action (see spectral_action_sigma.md), then:

```
Sigma_SM = D(rho_{beta/Q}(D^2) || rho_beta(D^2))
```

naturally contains the SU(3) contribution through the M_3(C) factor. The SU(3) gauge dynamics would follow from delta Sigma_SM = 0.

### 4.4 The M_3(C) factor and the octonions

A deep mathematical connection: A_F = C + H + M_3(C) is the MAXIMAL ASSOCIATIVE SUBALGEBRA of J_3(O).

Specifically (Boyle-Farnsworth 2018):
```
J_3(O) contains A_F = C + H + M_3(C)
```

This means:
- The octonionic structure (J_3(O)) is the FULL algebraic structure
- The Standard Model algebra (A_F) is the part that is ASSOCIATIVE
- The non-associative part of J_3(O) encodes ADDITIONAL structure (3 generations? mass ratios?)

For the Sigma framework: QRE requires associativity (matrix trace, matrix log). Therefore:

```
Sigma can only "see" the associative part of J_3(O) = A_F
```

This gives a NATURAL EXPLANATION for why the gauge group is G_SM and not something bigger: G_SM is the gauge group of the maximal ASSOCIATIVE subalgebra, and Sigma (which requires associativity) can only probe this subalgebra.

### 4.5 The Krasnov characterization

Krasnov (2019, arXiv:1912.11282) proved:

```
G_SM = Stab_{Spin(9)}(J)
```

where J is a specific element in the Lie algebra of Spin(9) (the even part of Cl(9)). In English: the Standard Model gauge group is the stabilizer of a particular algebraic element within Spin(9).

The connection to division algebras: Spin(9) is the spin group in the same dimension as the octonionic projective plane OP^2, and Cl(9) contains the Moufang projective plane structure.

For the Sigma framework: If Sigma on the Clifford algebra Cl(9) has a natural "preferred direction" J, then G_SM = Stab(J) would follow as the symmetry group that preserves this direction. The preferred direction could be related to the Khronon field (which defines a preferred time direction).

### 4.6 Assessment: Route 3

| Aspect | Status |
|--------|--------|
| A_F = C + H + M_3(C) gives G_SM uniquely | PROVEN (Connes) |
| M_3(C) -> SU(3) | Standard (gauge group of unitaries) |
| Spectral action gives SU(3) Yang-Mills | PROVEN (Connes-Chamseddine) |
| Sigma = Fisher info of spectral action | CONJECTURED (Paper 9) |
| SU(3) from Sigma via spectral action | CONTINGENT on Paper 9 |
| A_F = maximal associative subalgebra of J_3(O) | PROVEN (Boyle-Farnsworth) |
| "QRE sees only associative part" -> G_SM | NEW OBSERVATION |

**Verdict**: This is the MOST VIABLE route. SU(3) comes packaged with the full SM from the spectral action. The Sigma framework's contribution is the Fisher information identification (Paper 9). If Paper 9 succeeds, SU(3) follows automatically. No separate "Paper 8" would be needed -- it would be absorbed into Paper 9.

---

## 5. Route 4: Farnsworth's Nonassociative Spectral Geometry

### 5.1 The program (arXiv:2506.21496)

Farnsworth (2025) develops spectral geometry for nonassociative algebras. The key idea: replace the associative algebra A in Connes' spectral triple (A, H, D) with a nonassociative algebra (specifically, the octonions O or their complexification C tensor O).

The challenges are formidable:
1. The representation theory of nonassociative algebras on Hilbert spaces is not standard
2. The Dirac operator must commute with the algebra action (first-order condition), which is complicated by non-associativity
3. The spectral action Tr f(D^2/Lambda^2) uses the trace, which requires some form of associativity

Farnsworth's partial results:
- Gauge group of the octonionic spectral triple: G_2 x G_2 (not the SM!)
- The reduction G_2 x G_2 -> G_SM requires additional structure (a "breaking pattern")
- The work is in early stages; the full SM has not been recovered

### 5.2 G_2 x G_2 -> SU(3) x SU(2) x U(1)

The breaking pattern:

```
G_2 x G_2  ->  SU(3) x SU(2) x U(1) / Z_6
```

requires choosing a preferred direction in each G_2 factor. Specifically:
- First G_2 -> SU(3): stabilize one imaginary unit in Im(O) (7 -> 6 = 3 + 3_bar)
- Second G_2 -> SU(2) x U(1): stabilize two imaginary units

The number of broken generators: 14 + 14 - 8 - 3 - 1 = 16. These become massive gauge bosons (analogous to W, Z masses from Higgs mechanism).

### 5.3 Connection to Sigma

If Farnsworth's program succeeds in recovering the SM, the connection to Sigma would be:

```
Sigma on nonassociative spectral triple
  = D(rho_1 || rho_2) with rho_i defined via the octonionic Dirac operator
```

The non-associativity of O means that standard QRE cannot be used directly. But the RESTRICTED QRE on the associative subalgebra (= A_F = C + H + M_3(C)) is well-defined. This circles back to Route 3.

The genuinely NEW element from Farnsworth would be: the selection of A_F from O by the requirement of well-defined QRE. In other words:

**Conjecture**: The maximal subalgebra of the octonionic spectral algebra on which QRE is well-defined is exactly A_F = C + H + M_3(C), and the gauge group of this subalgebra is G_SM.

This would give a deep reason for why the SM gauge group is what it is: it is the largest gauge group compatible with the quantum information structure (QRE).

### 5.4 Assessment: Route 4

| Aspect | Status |
|--------|--------|
| Nonassociative spectral geometry | IN DEVELOPMENT (Farnsworth 2025) |
| G_2 x G_2 recovered | Partial (not full SM) |
| SU(3) from octonionic Dirac operator | NOT YET (need breaking pattern) |
| Connection to Sigma | SPECULATIVE (QRE on associative subalgebra) |
| Timeline | 6-12 months (wait for Farnsworth's next paper) |

**Verdict**: Mathematically the most novel route, but too early to assess viability. Monitor Farnsworth's progress. The conjecture "QRE selects A_F from O" is worth pursuing as an independent research direction.

---

## 6. Route 5: Confinement as Petz Recovery Failure

### 6.1 The physical picture

In QCD, confinement means:
- Quarks cannot exist as free particles
- Color charge cannot be isolated
- The potential between quarks grows linearly with distance: V(r) ~ sigma * r (string tension)
- At low energies, only color-singlet states (mesons, baryons) are observed

In the Sigma framework, the natural interpretation is:

```
Sigma_color = D(rho_quark || rho_singlet)
```

where rho_quark is the state of a colored quark and rho_singlet is the nearest color-singlet state. Confinement would mean:

```
Sigma_color -> infinity for an isolated quark
```

This is the MAXIMUM possible Petz recovery failure: the color information is completely irrecoverable from the observable (color-singlet) sector.

### 6.2 Entanglement entropy and the area law

The entanglement entropy of a region in QCD obeys an area law:

```
S_EE = sigma_ent * A / (4 l_s^2) + ...
```

where sigma_ent is the entanglement string tension and l_s is the string length scale. This area law is a signature of confinement: the entanglement is concentrated on the boundary between the region and its complement, with a "string" connecting the quarks.

In the Sigma framework, the entanglement entropy is related to the QRE via:

```
D(rho_A || rho_A^{thermal}) >= S_EE(A) - S_{thermal}(A)
```

The area law for S_EE implies a corresponding behavior for the QRE: Sigma_color grows with the separation surface area, which for an isolated quark (infinite separation) gives Sigma_color -> infinity.

### 6.3 Lattice QCD evidence

Lattice QCD calculations (Buividovich-Polikarpov 2008, arXiv:0811.3824; Ratti 2018) show:

1. The color entanglement entropy across a plane grows as sigma * A (area law)
2. At the deconfinement transition (T = T_c ~ 150 MeV), the area law transitions to a volume law
3. The QCD string tension sigma ~ (440 MeV)^2 sets the confinement scale

In Sigma language:
```
T < T_c (confined):   Sigma_color ~ sigma * A    (area law, Petz recovery fails)
T > T_c (deconfined): Sigma_color ~ sigma' * V   (volume law, partial recovery)
```

The deconfinement transition is a PHASE TRANSITION in the Petz recoverability:
- Below T_c: color information is permanently lost (confinement = complete Petz failure)
- Above T_c: color information becomes partially accessible (deconfinement = partial Petz recovery)

### 6.4 Non-associativity and confinement: a speculative connection

The observation from imaginary_units_gauge_forces_2026_03_19.md:

```
Octonions are non-associative: (e_1 * e_2) * e_3 != e_1 * (e_2 * e_3)
Quarks cannot exist alone (confinement)
```

Is there a deeper connection? Consider:
- Non-associativity means that the ORDER OF OPERATIONS matters
- For three octonionic elements a, b, c: the "product" depends on which pair is computed first
- This is analogous to the fact that three quarks (baryon) cannot be decomposed into two quarks + one quark in a unique way
- The "associator" [a, b, c] = (ab)c - a(bc) measures the failure of associativity
- The "confinement potential" V(r) = sigma * r measures the failure to isolate a quark

**Speculative conjecture (Furey-inspired)**:

```
V_confine(r) ~ Tr |[a, b, c]|^2 * r
```

where a, b, c are octonionic fields evaluated at the quark positions. The string tension sigma is proportional to the squared associator of the octonionic algebra.

This is highly speculative and lacks rigorous mathematical support. But it suggests a direction: confinement could be the PHYSICAL MANIFESTATION of algebraic non-associativity, detected by the Sigma framework as the impossibility of defining a consistent QRE for the full octonionic algebra.

### 6.5 Color confinement from quantum information: literature

Several groups have studied confinement from quantum information perspectives:

**Klebanov-Kutasov-Murugan (2007, arXiv:0709.2140)**: Entanglement entropy in confining gauge theories shows area-law behavior. The "entanglement string tension" agrees with the Wilson loop string tension.

**Nishioka-Takayanagi (2007, arXiv:0611035)**: Holographic entanglement entropy in confining backgrounds (AdS soliton) shows a phase transition from connected to disconnected minimal surfaces, corresponding to confinement/deconfinement.

**Radicevic (2016, arXiv:1608.04732)**: Entanglement entropy of lattice gauge theories with explicit computation of color entanglement. Shows that the color entanglement structure is related to the center symmetry of the gauge group.

**Casini-Huerta-Rosabal (2014, arXiv:1312.1183)**: The QRE for gauge fields includes edge-mode contributions. For SU(3), the edge modes contribute:

```
S_edge = (N_surface / 2) * ln dim(SU(3)) = (N_surface / 2) * ln 8
```

This should be compared with the Level 2 result for SU(2): S_edge = (N/2) * ln 3. The SU(3) edge entropy is ln 8, which equals dim Im(O) + 1 = 8. This is NOT a coincidence -- it reflects the fact that dim(SU(3)) = 8 = dim(O).

### 6.6 The Petz map for color recovery

Consider the "color channel" N_color: rho_full -> Tr_color(rho_full), which traces out the color degrees of freedom. This is a CPTP map from the full QCD Hilbert space to the color-singlet sector.

The Petz recovery map would be:

```
R_Petz(sigma) = rho_0^{1/2} N_color^dagger(N_color(rho_0)^{-1/2} sigma N_color(rho_0)^{-1/2}) rho_0^{1/2}
```

where rho_0 is the QCD vacuum state.

For a color-singlet state |psi>: N_color(|psi><psi|) preserves full information. R_Petz recovers exactly. Sigma_color = 0.

For a colored state |q> (isolated quark): N_color(|q><q|) loses color information completely (the trace over color gives the maximally mixed state in the color sector). R_Petz fails maximally. Sigma_color -> infinity (formally).

**The connection to confinement**: Confinement means that the physical Hilbert space IS the color-singlet sector. The Petz recovery map is TRIVIALLY exact on physical states (Sigma_color = 0) and INFINITELY fails on unphysical states (Sigma_color = infinity). This is consistent with, but does not derive, confinement.

**What would constitute a derivation**: Show that the DYNAMICS of Sigma (delta Sigma = 0 or Sigma extremization) REQUIRE the physical Hilbert space to be the color-singlet sector. This would mean: confinement follows from the variational principle of the Sigma framework.

### 6.7 Assessment: Route 5

| Aspect | Status |
|--------|--------|
| Confinement = Sigma_color -> infinity | INTERPRETIVE (consistent, not derived) |
| Area law from Sigma | KNOWN (lattice QCD, holography) |
| Deconfinement transition = Petz phase transition | NEW INTERPRETATION |
| Non-associativity -> confinement | HIGHLY SPECULATIVE |
| Edge mode entropy ln 8 = dim O | NEW OBSERVATION (to verify) |
| Dynamics of confinement from delta Sigma = 0 | NOT DERIVED |

**Verdict**: Rich interpretive framework, but not a derivation. Confinement cannot be derived without having SU(3) gauge dynamics first (circular reasoning). This route is valuable as an INTERPRETATION of confinement within the Sigma framework, not as a derivation of SU(3).

---

## 7. Synthesis: The Three-Level Understanding

### 7.1 Level A: Algebraic (what the gauge group IS)

```
J_3(O) -- maximal associative subalgebra --> A_F = C + H + M_3(C)
                                                       |
                                         gauge group = U(A_F)/U(center)
                                                       |
                                              G_SM = U(1) x SU(2) x SU(3)
```

SU(3) = unitary group of M_3(C). This is the ALGEBRAIC origin of SU(3).

Sigma's contribution: QRE requires associativity, so Sigma can only "see" the associative part of J_3(O) = A_F. This selects G_SM as the maximal gauge group compatible with the QRE structure.

### 7.2 Level B: Geometric (how the dynamics work)

```
11D spacetime = 4D + 7D(CP^2 x S^2 x S^1)
                         |
              isometry = SU(3) x SU(2) x U(1)
                         |
          KK reduction -> Yang-Mills dynamics
                         |
     Sigma_{11D} = Sigma_grav + Sigma_{SU(3)} + Sigma_{SU(2)} + Sigma_{U(1)}
```

SU(3) gauge dynamics come from the SU(3) isometry of the internal space (CP^2 or its octonionic generalization S^7/G_2).

Sigma's contribution: the entropy decomposition Sigma_{11D} = sum over gauge sectors, with each sector contributing proportional to its field strength.

### 7.3 Level C: Information-theoretic (why confinement)

```
Color channel: N_color(rho) = Tr_color(rho)
                    |
  Sigma_color = D(rho || R_Petz(N_color(rho)))
                    |
  Confinement: Sigma_color = 0 for singlets, infinity for colored states
                    |
  Deconfinement = Petz recovery phase transition at T = T_c
```

Confinement is the statement that the color channel is COMPLETELY IRREVERSIBLE for isolated quarks.

Sigma's contribution: confinement is the extreme case of Petz recovery failure, connecting QCD to the deepest principle of the framework (tau = 1 = complete irreversibility).

### 7.4 How the three levels connect

```
Level A (algebra):  A_F = C + H + M_3(C)  --- selects the GROUP ---
         |
         | (spectral action on A_F)
         v
Level B (geometry): Spectral action = EH + YM + Higgs  --- gives DYNAMICS ---
         |
         | (Sigma = Fisher info of spectral action)
         v
Level C (info):     Confinement = Petz failure  --- gives INTERPRETATION ---
```

The three levels are not alternatives but LAYERS of the same structure:
1. A_F determines the gauge group (Level A)
2. The spectral action on A_F gives the gauge dynamics (Level B)
3. Sigma measures the information-theoretic content of the dynamics (Level C)

---

## 8. New Observations and Conjectures

### 8.1 "QRE selects A_F from J_3(O)" (NEW)

**Conjecture**: The requirement that QRE be well-defined (associativity, trace, logarithm) selects the maximal associative subalgebra of J_3(O), which is A_F = C + H + M_3(C). The gauge group G_SM is therefore the maximal gauge group compatible with quantum information theory applied to the octonionic exceptional Jordan algebra.

**Significance**: This would answer "why SU(3) and not something else?" -- because SU(3) is the largest simple Lie group that can arise from an associative algebra within the octonionic framework.

**Test**: Verify that no associative subalgebra of J_3(O) larger than A_F exists. (This is known to be true -- proven by Albert and Jacobson.)

### 8.2 Deconfinement as Petz phase transition (NEW)

**Conjecture**: The QCD deconfinement transition at T = T_c is a phase transition in the Petz recoverability of color information. Below T_c, the color channel is perfectly irreversible (Sigma_color = infinity for colored states). Above T_c, partial color recovery becomes possible (Sigma_color finite but large).

**Test**: Compute the Petz recovery fidelity for the color channel on the lattice at various temperatures. The fidelity should jump at T = T_c.

### 8.3 Edge mode entropy ln 8 = dim O (NEW OBSERVATION)

The Donnelly-Wall edge mode entropy for SU(3) is:

```
S_edge^{SU(3)} propto ln dim(SU(3)) = ln 8
```

And dim(O) = 8 (the octonions as a vector space). The coincidence ln 8 = ln(dim O) suggests that the SU(3) edge modes are counting octonionic degrees of freedom.

Compare with SU(2):
```
S_edge^{SU(2)} propto ln dim(SU(2)) = ln 3 = ln(dim Im(H))
```

For SU(2), the edge entropy counts the IMAGINARY quaternionic dimensions (3 = dim Im(H)).
For SU(3), the edge entropy counts the FULL octonionic dimensions (8 = dim O).

The discrepancy (imaginary vs full) arises because SU(3) has 8 generators while Im(O) has 7 dimensions. The 8th generator (the identity direction in the Gell-Mann basis, or equivalently lambda_8) corresponds to the REAL part of O.

**Refined observation**: dim(SU(N)) = N^2 - 1 for N >= 2.
- SU(2): 4 - 1 = 3 = dim Im(H)
- SU(3): 9 - 1 = 8 = dim O (total, not just imaginary!)

The pattern is:
```
SU(2) generators <-> imaginary quaternions (Im H)
SU(3) generators <-> full octonions (O, all 8 dimensions including real)
```

This is because SU(2) = unit quaternions = S^3 in Im(H), but SU(3) is NOT the unit octonions (which form S^7, not a Lie group). Instead, SU(3) is the subgroup of G_2 that preserves a quaternionic subspace, and its 8 generators span O as a vector space.

### 8.4 Octonionic Khronon and the 15D KK (SPECULATIVE)

Extending the complex Khronon (Paper 6) and quaternionic Khronon (Paper 7 plan):

```
Complex Khronon:     phi = f * e^{i theta}        in C  -- 2 real dof -- gravity + U(1)
Quaternionic Khronon: Phi = f * q                  in H  -- 4 real dof -- gravity + SU(2) + U(1)
Octonionic Khronon:   Omega = f * o               in O  -- 8 real dof -- gravity + G_2 => SU(3) + ...
```

But the octonionic Khronon has a FATAL PROBLEM: the non-associativity of O means that the kinetic term:

```
K(|d Omega|^2) = K(g^{mu nu} partial_mu Omega^* partial_nu Omega)
```

is ambiguous unless the product is carefully defined. For the complex and quaternionic cases, the products are associative (or at least alternative), so the kinetic term is unambiguous.

For octonions, one can use the ALTERNATIVITY property: a(ab) = (aa)b and (ab)b = a(bb). This means the norm |Omega|^2 = Omega^* Omega is well-defined, and therefore the kinetic term can be defined using the octonionic norm. But the gauge connection derived from the phase variation is more problematic.

**Alternative**: Use the SPLIT OCTONION approach. The split octonions O_s have the same dimension (8) but a different signature (4,4) and are ISOMORPHIC to 2x2 matrices over the quaternions: O_s = M_2(H). This is associative (being a matrix algebra), and could provide a "regularized" version of the octonionic Khronon.

This is highly speculative and requires significant mathematical development.

### 8.5 The asymptotic freedom connection (NEW INTERPRETATION)

QCD is asymptotically free: the coupling g_3(mu_RG) decreases at high energies. In the Sigma framework:

```
Sigma_{RG}^{QCD} = 2 ln(mu_RG / Lambda_QCD)
```

At high energies (mu_RG >> Lambda_QCD): Sigma_{RG}^{QCD} is large, meaning the high-energy QCD state is very different from the low-energy vacuum. The Petz recovery from high to low energies is poor.

At low energies (mu_RG ~ Lambda_QCD): Sigma_{RG}^{QCD} -> 0, meaning the state approaches the vacuum. But SIMULTANEOUSLY, confinement sets in, making Sigma_color -> infinity.

These are two DIFFERENT Sigma's:
- Sigma_{RG}: measures the "distance" from the vacuum in energy space (decreases at low energy)
- Sigma_color: measures the "distance" from color-singlet in color space (increases at low energy)

They are complementary: at low energy, the state is close to the vacuum (small Sigma_{RG}) but color is completely confined (large Sigma_color). At high energy, the state is far from the vacuum (large Sigma_{RG}) but color is liberated (small Sigma_color).

**Conjecture**: There is a DUALITY between Sigma_{RG} and Sigma_color:

```
Sigma_{RG}^{QCD} + Sigma_color >= C
```

where C is a constant related to the QCD scale. This would be an information-theoretic formulation of the complementarity between asymptotic freedom and confinement.

---

## 9. Comparison of Routes

| Criterion | Route 1 (Octonions) | Route 2 (KK 11D) | Route 3 (Connes NCG) | Route 4 (Farnsworth) | Route 5 (Confinement) |
|-----------|:-------------------:|:-----------------:|:--------------------:|:--------------------:|:---------------------:|
| Identifies SU(3) | via G_2 / J_3(O) | via isometry | via M_3(C) | via octonionic D | via color channel |
| Gives dynamics | No | Yes (KK) | Yes (spectral action) | In development | No (interpretation) |
| Connects to Sigma | Indirect (QRE selects A_F) | Direct (Sigma_{11D} decomposition) | Direct (Fisher info) | Speculative | Direct (Petz failure) |
| Non-associativity handled | Not resolved | Not needed (geometric) | Avoided (A_F associative) | Central challenge | Not needed |
| Mathematical maturity | High (division algebras well-studied) | High (KK textbook) | High (Connes 1996-2024) | Low (2025, in progress) | Medium (lattice QCD) |
| Novelty for Sigma framework | Medium | Low (standard KK + Sigma decomposition) | HIGH (Paper 9) | Very High | Medium |
| Feasibility | 6-12 months | Already done (Section 3.2) | 3-6 months (Paper 9) | 12+ months (wait for Farnsworth) | Interpretive only |
| Self-contained derivation | No | No (geometric input) | YES (if Paper 9 works) | Possibly | No |

**Winner: Route 3 (Connes NCG)**. It is the only route that provides BOTH the gauge group AND the dynamics, in a framework (spectral action) that already has a proven connection to entropy (CCSvS 2018) and therefore to Sigma.

---

## 10. Concrete Plan for Paper 8

### 10.1 Revised assessment

The original plan had Paper 8 as "SU(3) via O [6-12 months, wait for Farnsworth]". Based on this analysis, the plan should be revised:

**Option A (recommended)**: MERGE Paper 8 into Paper 9.

Paper 9 (Sigma on A_F = spectral action) gives the FULL SM including SU(3). There is no need for a separate Paper 8 if Paper 9 succeeds. The octonionic perspective and confinement interpretation can be included as sections within Paper 9.

**Option B (if Paper 9 fails)**: Write Paper 8 as an interpretive paper.

If Paper 9 fails (Sigma != Fisher info of spectral action), Paper 8 would be limited to:
- The "QRE selects A_F from J_3(O)" conjecture
- The confinement = Petz failure interpretation
- The edge mode entropy observation (ln 8 = dim O)
- The Sigma_{11D} KK decomposition

This would be a valuable conceptual paper but would NOT derive SU(3) from Sigma.

### 10.2 Timeline

**Months 1-3**: Paper 9 toy model (S^1 x M_2(C)). Go/No-Go.

**Month 3 decision point**:
- If Paper 9 GO: continue Paper 9 with full A_F. SU(3) comes automatically. No separate Paper 8.
- If Paper 9 NO-GO: pivot to Paper 8 as interpretive paper (Options below).

**Months 3-6 (if Paper 9 NO-GO)**:
- Formalize the "QRE selects A_F" conjecture with rigorous proof of maximal associativity
- Compute Sigma_color on lattice (collaboration needed)
- Write the Sigma_{11D} decomposition explicitly
- Explore Farnsworth's nonassociative spectral geometry

### 10.3 Mathematical prerequisites

For the Connes route (Paper 9):
1. Master the spectral action expansion (van Suijlekom textbook 2nd ed.)
2. Understand the CCSvS entropy-spectral action correspondence (arXiv:1809.02944)
3. Compute the QRE between Gibbs states on almost-commutative geometries
4. Identify where M_3(C) contributes to the Fisher information

For the octonionic route:
1. Read Boyle (arXiv:2006.16265) on SM from C tensor H tensor O
2. Read Furey (arXiv:2505.07923) on division algebra superalgebra
3. Read Farnsworth (arXiv:2506.21496) on nonassociative spectral geometry
4. Understand the J_3(O) -> A_F maximal associative subalgebra theorem
5. Read Todorov-Drenska (arXiv:1911.13124) on octonionic SM

For the confinement route:
1. Read Buividovich-Polikarpov (arXiv:0811.3824) on entanglement entropy in lattice gauge theory
2. Read Nishioka-Takayanagi (arXiv:0611035) on holographic entanglement in confining backgrounds
3. Read Donnelly-Wall (arXiv:1506.04267) on edge modes
4. Read Radicevic (arXiv:1608.04732) on lattice gauge entanglement

### 10.4 Key literature (complete list)

**SU(3) from octonions**:
- Gunaydin-Gursey 1973: original octonion -> SU(3) color
- Dixon 1994: "Division Algebras: Octonions, Quaternions, Complex Numbers and the Algebraic Design of Physics" (book)
- Furey 2012-2025: C tensor H tensor O minimal ideals -> one generation SM fermions
- Furey 2025 (arXiv:2505.07923): Z_2^5-graded superalgebra from division algebras
- Boyle 2020 (arXiv:2006.16265): SM from C tensor H tensor O and J_3(O)
- Singh 2025 (arXiv:2508.10131): fermion mass ratios from J_3(O)
- Krasnov 2019 (arXiv:1912.11282): G_SM = Stab_{Spin(9)}(J)
- Todorov-Drenska 2019 (arXiv:1911.13124): octonions, exceptional Jordan algebra, and SM
- Gourlay-Gresnigt 2026 (arXiv:2601.07857): Cl(10) + S_3

**Connes NCG**:
- Connes 1996 (hep-th/9606001): spectral action principle
- Chamseddine-Connes 2006: scale-invariant spectral action
- Chamseddine-Connes-van Suijlekom 2018 (arXiv:1809.02944): entropy = spectral action
- Van Suijlekom 2024: "Noncommutative Geometry and Particle Physics" 2nd ed.
- Boyle-Farnsworth 2018 (arXiv:1604.00847): differential graded algebra
- Connes 2025 (arXiv:2511.08159): spectral torsion of SM NCG

**Nonassociative spectral geometry**:
- Farnsworth 2025 (arXiv:2506.21496): nonassociative spectral geometry, G_2 x G_2
- Boyle-Farnsworth 2020 (arXiv:1910.11888): Jordan geometry

**Entanglement and gauge/confinement**:
- Donnelly-Wall 2015 (arXiv:1506.04267): edge modes and gauge entanglement
- Casini-Huerta-Rosabal 2014 (arXiv:1312.1183): QRE for gauge fields
- Buividovich-Polikarpov 2008 (arXiv:0811.3824): entanglement entropy in lattice gauge theory
- Nishioka-Takayanagi 2007 (arXiv:0611035): holographic entanglement in confining backgrounds
- Radicevic 2016 (arXiv:1608.04732): entanglement in lattice gauge theories
- Klebanov-Kutasov-Murugan 2007 (arXiv:0709.2140): entanglement in confining theories

**KK and M-theory**:
- Witten 1981: "Search for a Realistic Kaluza-Klein Theory"
- Cremmer-Julia-Scherk 1978: 11D supergravity
- Duff-Nilsson-Pope 1986: "Kaluza-Klein Supergravity" (review)

**Other**:
- Harlow 2018 (arXiv:1802.01040): gauge symmetry = QEC
- Szangolies 2025 (arXiv:2512.17328): n-qubit entanglement -> SM gauge group
- Bianconi 2025 (PRD 111, 066001): gravity from entropy

---

## 11. Summary: What Is Derived, What Is Assumed, What Is Open

| Result | Status | Route |
|--------|--------|-------|
| SU(3) structure from M_3(C) in A_F | KNOWN (Connes 1996) | 3 |
| A_F = maximal associative subalgebra of J_3(O) | KNOWN (Albert, Boyle) | 1/3 |
| G_SM = Stab_{Spin(9)}(J) | KNOWN (Krasnov 2019) | 1 |
| Spectral action gives SU(3) Yang-Mills | PROVEN (Connes) | 3 |
| Sigma_{11D} = Sigma_grav + Sigma_{SU(3)} + ... | NEW (Section 3.2) | 2 |
| Sigma = Fisher info of spectral action | CONJECTURED (Paper 9) | 3 |
| SU(3) from Sigma (Level 3 derivation) | NOT YET | -- |
| "QRE selects A_F from J_3(O)" | NEW CONJECTURE | 1/3 |
| Confinement = Sigma_color -> infinity (Petz failure) | NEW INTERPRETATION | 5 |
| Deconfinement = Petz phase transition | NEW INTERPRETATION | 5 |
| Edge mode entropy ln 8 = dim O | NEW OBSERVATION | 5 |
| Sigma_{RG} + Sigma_color duality | SPECULATIVE CONJECTURE | 5 |
| Non-associativity -> confinement | HIGHLY SPECULATIVE | 1/5 |
| Octonionic Khronon | BLOCKED (non-associativity) | 4 |

---

## 12. Attribution

| Discovery | Who | Ours? |
|-----------|-----|-------|
| Octonions -> SU(3) color | Gunaydin-Gursey (1973) | No |
| Division algebra ladder R,C,H,O -> SM | Dixon (1994), Furey (2012+) | No |
| A_F = C+H+M_3(C) -> SM spectral triple | Connes (1996) | No |
| A_F = maximal associative subalgebra of J_3(O) | Boyle-Farnsworth (2018) | No |
| G_SM = Stab_{Spin(9)}(J) | Krasnov (2019) | No |
| Entropy = spectral action | CCSvS (2018) | No |
| n-qubit -> Hopf -> SM | Szangolies (2025) | No |
| Nonassociative spectral geometry | Farnsworth (2025) | No |
| Edge modes for gauge fields | Donnelly-Wall (2015) | No |
| Entanglement in confining theories | Multiple groups | No |
| **Sigma_{11D} decomposition for SU(3)** | **New (this analysis)** | **Yes** |
| **"QRE selects A_F from J_3(O)" conjecture** | **New (this analysis)** | **Yes** |
| **Confinement = Sigma_color -> infinity (Petz failure)** | **New interpretation** | **Yes** |
| **Deconfinement = Petz recovery phase transition** | **New interpretation** | **Yes** |
| **Edge mode ln 8 = dim O observation** | **New observation** | **Yes** |
| **Sigma_{RG} + Sigma_color complementarity** | **New conjecture** | **Yes** |
| **Paper 8 merged into Paper 9 recommendation** | **New strategy** | **Yes** |

---

## Record
Created 2026-03-19 during systematic investigation of SU(3) from Sigma framework.
Builds on: su2_from_sigma.md, imaginary_units_gauge_forces_2026_03_19.md, matter_sector_strategy_2026_03_19.md, spectral_action_sigma.md, gap_analysis_2026_03_19.md, division_algebra_attribution_2026_03_19.md.
