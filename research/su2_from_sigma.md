# SU(2) from Sigma: Derivation Attempt, Weinberg Angle, and the Non-Abelian Gap
## 2026-03-19
## Sheng-Kai Huang

---

## Executive Summary

This document investigates whether the Sigma framework (Sigma = D(rho_spacetime || rho_matter)) can derive SU(2) gauge symmetry, and whether the Weinberg angle constrains the Khronon mass parameter mu. Five approaches are examined in detail. The conclusion is that **a complete Level 3 derivation of SU(2) from Sigma does not yet exist**, but three promising partial routes are identified, with concrete mathematical obstacles and potential resolutions mapped out.

**Key results**:
1. The fundamental obstacle to SU(2) from Sigma is identified precisely: QRE is defined on complex Hilbert space with tensor products, which structurally produces U(1) but not non-abelian structure without additional input.
2. The Weinberg angle does NOT directly constrain mu (the two mu's are unrelated), but a suggestive structural parallel exists through the Sigma_RG framework.
3. The most promising route is the **quaternionic Khronon doublet** combined with the **7D Kaluza-Klein** decomposition, where the fiber S^3 = SU(2) provides the gauge structure while Sigma provides the dynamics.
4. The Harlow QEC framework offers a second independent route: gauge symmetry as a quantum error-correcting code, with the Petz map as the recovery operator.

---

## 1. Why Level 3 Failed for SU(2): Precise Diagnosis

### 1.1 What Level 2 achieved

The QRE chain rule on 2-qubit systems gives:

```
D(rho || I/4) = D(rho || G_{SU(2)}(rho)) + D(G_{SU(2)}(rho) || I/4)
```

where G_{SU(2)} is the SU(2) twirl (Marvian-Spekkens 2014). For the triplet state |Phi+>:

- D(rho || G(rho)) = ln 3 (SU(2)-breaking part, 3 directions)
- 3 directions = dim Im(H) = {i, j, k}
- su(2) isomorphic to Im(H) as Lie algebras

This is an algebraic identity, verified numerically to 10^{-16} precision (2026-03-19). It says: "The QRE of a 2-qubit state naturally decomposes into an SU(2)-invariant part and an SU(2)-breaking part with exactly 3 independent directions."

### 1.2 What Level 3 requires but does not have

Level 3 would require showing that the **dynamics** of the Sigma framework (field equations from delta Sigma = 0) produce SU(2) gauge fields. Specifically:

**Missing link 1: From chain rule to gauge connection.**
The QRE chain rule is a kinematic identity on the state space. It tells us that SU(2) structure exists in the 2-qubit Hilbert space, but it does not produce a gauge connection A^a_mu (a = 1,2,3) with the correct Yang-Mills dynamics:

```
D_nu F^{a,mu nu} + g_2 epsilon^{abc} A^b_nu F^{c,mu nu} = J^{a,mu}
```

The non-abelian term (epsilon^{abc} A^b F^c) is the crux: it arises from the non-commutativity [T^a, T^b] = i epsilon^{abc} T^c, which is a property of the gauge connection, not of the QRE decomposition.

**Missing link 2: From Im(H) to gauge dynamics.**
The isomorphism su(2) = Im(H) is a Lie algebra identity. But a Lie algebra does not uniquely determine a gauge theory -- one also needs:
- A principal bundle P(M, SU(2))
- A connection 1-form on P
- Matter fields in specific representations
- The Yang-Mills action S = -(1/4g^2) int Tr(F wedge *F)

The Sigma framework provides none of these.

**Missing link 3: The tensor product obstruction.**
Paper 6 Theorem 1 excluded quaternionic Hilbert space precisely because H lacks tensor products. But SU(2) gauge theory is naturally quaternionic (the fundamental representation of SU(2) acts on C^2 = H^1, and the adjoint acts on Im(H)). This creates a tension:

- Sigma requires C (complex) to be defined
- SU(2) is naturally described by H (quaternions)
- The bridge between these two is: SU(2) = Sp(1) = unit quaternions, which acts on C^2 by left multiplication

This bridge exists but is not automatic from the Sigma framework.

### 1.3 The non-abelian gap: precise mathematical statement

**Theorem (Non-abelian gap).** Let Sigma = D(rho || sigma) be the quantum relative entropy on a complex Hilbert space H = H_1 tensor H_2 with dim H_1 = dim H_2 = 2. The QRE chain rule under the SU(2) twirl G decomposes D into SU(2)-covariant components, but the resulting decomposition is:
1. Kinematic (state-space identity), not dynamic (no field equation)
2. Linear in rho (QRE is convex), while the Yang-Mills action is quadratic in A
3. Abelian in structure: D(rho || G(rho)) decomposes as a sum over irreps, with no cross-terms between irreps corresponding to the non-abelian structure constants

The gap between the QRE chain rule and Yang-Mills dynamics is therefore **structural, not merely computational**. It cannot be bridged by a more careful calculation within the current framework. New mathematical structure is needed.

### 1.4 Why it worked for U(1) but fails for SU(2)

For U(1), the derivation in Paper 6 did not derive the gauge dynamics from QRE either. What it showed was:

1. Sigma requires complex Hilbert space (5 routes)
2. Complex Hilbert space has a U(1) phase symmetry
3. Gauging U(1) gives electromagnetism (standard argument)
4. KK decomposition with Sigma-dependent dilaton gives modified Maxwell

The key step (3) uses the **gauge principle** as an external input. The same strategy for SU(2) would require:

1. Sigma on 2-qubit systems has SU(2) structure (Level 2, verified)
2. SU(2) is the natural symmetry of the 2-qubit state space (Hopf fibration S^7 -> S^4)
3. Gauging SU(2) gives weak force (standard argument)
4. KK decomposition in 7D with Sigma-dependent moduli gives Yang-Mills

Steps 1-2 are done. Steps 3-4 require the gauge principle and KK geometry as external inputs. This is exactly the Level 2 situation: the Sigma framework identifies the correct gauge group but does not derive the gauge dynamics from first principles.

---

## 2. Literature: Gauge Symmetry from Information/Entropy

### 2.1 Harlow: Gauge symmetry and quantum error correction

Harlow (2018, arXiv:1802.01040) proved a remarkable theorem: **gauge symmetry is equivalent to a quantum error-correcting code structure**. Specifically:

- A gauge theory with gauge group G has a code subspace C = {|psi>: G|psi> = |psi>} (gauge-invariant states)
- The physical Hilbert space is the code subspace
- Gauss's law is the error-detection condition
- The Petz map recovers gauge-invariant information from local measurements

**Connection to Sigma framework**: The Petz recovery map is central to Papers 1-6. If gauge symmetry = QEC, and the Petz map = recovery operator, then:

```
tau_{gauge} = 1 - F(rho, R_{Petz}(N(rho)))
```

where N is the "gauging" channel (projecting onto the gauge-invariant subspace). For a non-abelian gauge theory, tau_{gauge} = 0 for gauge-invariant states (perfect recovery) and tau_{gauge} > 0 for gauge-variant states (information is "lost" to the gauge orbit).

This gives a new interpretation: **Sigma_{weak} = D(rho || G_{SU(2)}(rho))** is the "cost" of gauge-invariance, measured by QRE. The chain rule decomposition is the QEC structure of the weak force.

**Obstacle**: Harlow's framework assumes the gauge theory exists and then identifies the QEC structure. It does not derive the gauge theory from QEC. The arrow of reasoning needs to be reversed.

### 2.2 Swingle and Van Raamsdonk: Entanglement and gravity

Van Raamsdonk (2010) showed that entanglement builds spacetime geometry. Swingle (2012) connected this to tensor networks. But their work is about gravity (the geometric side of Sigma), not about gauge forces.

The key insight from their work for SU(2) is: **if entanglement builds geometry (gravity), perhaps multi-partite entanglement builds gauge connections (forces)**. This is precisely the Szangolies hierarchy:

```
1-qubit entanglement → U(1) geometry (Bloch sphere)
2-qubit entanglement → SU(2) geometry (quaternionic structure)
3-qubit entanglement → SU(3) geometry (octonionic structure)
```

But making this precise requires going from "entanglement structure" to "gauge dynamics," which is the non-abelian gap identified above.

### 2.3 Donnelly, Wall, Casini: Entanglement entropy and gauge fields

Donnelly and Wall (2015, arXiv:1506.04267) showed that the entanglement entropy for a gauge field has an additional "edge mode" contribution proportional to ln |G| (where G is the gauge group). For SU(2):

```
S_{edge} = (N_surface/2) ln dim(SU(2)) = (N_surface/2) ln 3
```

This is structurally identical to D(rho || G(rho)) = ln 3 from the Level 2 calculation. The Donnelly-Wall edge modes are the gauge degrees of freedom on the entangling surface, and they contribute exactly the same amount as the SU(2)-breaking QRE.

**This is a new observation**: the Level 2 result D_{break} = ln 3 for SU(2) is the Donnelly-Wall edge mode entropy. This connects the Sigma chain rule to the gauge theory entanglement entropy directly.

### 2.4 Casini, Huerta, Rosabal: Relative entropy and gauge

Casini, Huerta, and Rosabal (2014, arXiv:1312.1183) studied the relative entropy for gauge fields and showed that the QRE in a gauge theory decomposes into "electric center" and "magnetic center" contributions, with the gauge-invariant part being the physical relative entropy. Their framework is compatible with the Sigma decomposition.

### 2.5 Bianconi: Dirac-Kahler gauge from entropy

Bianconi (2025, PRD 111, 066001) constructs gravity from L = -Tr_F ln(G_tilde g_tilde^{-1}), which is the 0-form sector of a Dirac-Kahler construction. Opening up the 1-form sector naturally introduces gauge fields. The Sigma framework is the 0-form sector of Bianconi's construction (confirmed 2026-03-19).

If this identification holds, then: **Sigma on 1-forms = Yang-Mills gauge dynamics**. The 1-form in Bianconi's framework carries gauge indices, and the QRE between the 1-form field and its reference state should reproduce the Yang-Mills action in the weak-field limit.

This is the most concrete path to SU(2) from Sigma, but the detailed calculation has not been done.

---

## 3. The Weinberg Angle Connection

### 3.1 sin^2(theta_W) in the Standard Model

The Weinberg angle theta_W relates the SU(2)_L coupling g_2 to the U(1)_Y coupling g_1:

```
tan(theta_W) = g_1/g_2
sin^2(theta_W) = g_1^2 / (g_1^2 + g_2^2)
```

Measured value at the Z pole: sin^2(theta_W) = 0.23122 +/- 0.00003

At the GUT scale (from SU(5) or SO(10) unification): sin^2(theta_W) = 3/8 = 0.375

In Connes' NCG: sin^2(theta_W) = 3/8 at the unification scale is a **prediction**, not an input. It follows from the structure of A_F = C + H + M_3(C).

### 3.2 Two different mu's: no direct connection

**mu_{RG}** = renormalization scale (energy scale parameter in the running coupling):
```
d g_i / d(ln mu_{RG}) = beta_i(g_i)
```

**mu_{Khronon}** = Khronon mass parameter (mu^{-1} = 22.3 Mpc in Blanchet-Skordis):
```
K(Q) = mu^2 (Q-1)^2 + ...
```

These are completely different quantities:
- mu_{RG} has dimensions of energy, ranges from ~0 to ~10^{19} GeV
- mu_{Khronon} has dimensions of inverse length (or mass in natural units), with mu^{-1} ~ 22.3 Mpc ~ 10^{-33} eV

There is no a priori relationship. mu_{RG} is a running parameter in perturbative QFT; mu_{Khronon} is a fixed mass parameter in the Khronon Lagrangian.

### 3.3 But: a structural parallel through Sigma_RG

From the supplement (Section S4), the RG running of alpha has a Sigma structure:

```
alpha^{-1}(E) = alpha^{-1}(E_0) - (1/3pi) Sigma_{RG}
Sigma_{RG} = 2 ln(E/E_0)
```

This has the same form as Sigma_grav = 2 ln Q. The parallel suggests that RG flow is a form of "coarse-graining entropy production," which the Sigma framework should capture.

For the electroweak sector, the running of the Weinberg angle is:

```
sin^2(theta_W)(mu_{RG}) = sin^2(theta_W)(M_Z)
    + (alpha/4pi) * [sum over fermion charges] * ln(mu_{RG}/M_Z) + ...
```

In Sigma language:

```
sin^2(theta_W)(mu_{RG}) = sin^2(theta_W)(M_Z) + c * Sigma_{RG}(mu_{RG}, M_Z)
```

where c is a numerical coefficient determined by the matter content.

### 3.4 Could the Weinberg angle constrain mu_{Khronon}?

**Direct route: NO.** The two mu's live in different sectors (RG scale vs Khronon mass). There is no equation connecting them.

**Indirect route: MAYBE, through the KK/NCG connection.** If the Sigma framework at the KK level (Sections 4-5 below) relates the gauge couplings to the geometry of extra dimensions, and if the extra dimensions are related to the Khronon field, then:

```
sin^2(theta_W) = Vol(U(1)_Y fiber) / Vol(SU(2)_L x U(1)_Y fiber)
```

In standard KK, sin^2(theta_W) = 3/8 comes from the ratio of the U(1) and SU(2) fiber volumes in the Connes spectral triple. If the fiber size is set by mu_{Khronon} somehow, there could be a relation.

**Concrete check**: In the 7D KK decomposition (Section 5), the SU(2) gauge coupling is:

```
g_2^2 = (16 pi G) / (Vol(S^3) * R_{KK}^3)
```

where R_{KK} is the KK radius. If R_{KK} = mu_{Khronon}^{-1}, then:

```
g_2^2 = (16 pi G * mu^3) / (2 pi^2)
g_2^2 = 8 G mu^3 / pi
```

With mu^{-1} = 22.3 Mpc = 6.87 x 10^{23} m, mu ~ 1.46 x 10^{-24} m^{-1}:

```
g_2^2 ~ 8 * (6.67 x 10^{-11}) * (1.46 x 10^{-24})^3 / pi
      ~ 8 * 6.67e-11 * 3.11e-73 / 3.14
      ~ 5.3e-83
```

This is absurdly small (measured g_2 ~ 0.65, so g_2^2 ~ 0.42).

**Conclusion: mu_{Khronon} is far too small to set the SU(2) coupling constant through KK radius identification.** The KK radius for SU(2) must be at the Planck scale (~10^{-35} m), not at the Khronon scale (~10^{24} m). The Weinberg angle does not constrain mu_{Khronon}.

### 3.5 Connes NCG prediction: sin^2(theta_W) = 3/8

In the Connes framework, the 3/8 prediction at unification comes from:

```
A_F = C + H + M_3(C)
Tr_F(Y^2) / Tr_F(T_3^2) = 5/3  (hypercharge normalization)
sin^2(theta_W) = g_1^2/(g_1^2+g_2^2) = 3/(3+5) = 3/8
```

This is identical to the SU(5) GUT prediction. The factor 5/3 comes from the specific embedding of U(1)_Y in the finite algebra A_F.

**For our framework**: If Paper 9 succeeds in connecting Sigma to the spectral action on A_F, then sin^2(theta_W) = 3/8 at unification becomes a consequence of our framework (inherited from Connes). But it is not a constraint on mu. It is a constraint on the algebra A_F.

---

## 4. Quaternionic Structure of Sigma

### 4.1 Standard Sigma on complex Hilbert space

The standard QRE is defined on density operators over C^n:
```
D(rho || sigma) = Tr(rho ln rho - rho ln sigma)
```

This uses the complex trace, complex matrix logarithm, and complex eigendecomposition. The result is always real and non-negative (for states).

### 4.2 Quaternionic Hilbert space: Adler's framework

Adler (1995) developed quaternionic quantum mechanics (QQM) on a right H-module (quaternionic Hilbert space). Key features:

- States: |psi> in H^n (n-tuples of quaternions)
- Inner product: <psi|phi> = sum_i (psi_i)* phi_i in H (quaternion-valued)
- Operators: right-H-linear maps (A|psi>q = A(|psi>q) for q in H)
- Density operators: rho = sum_i p_i |psi_i><psi_i|, Hermitian (rho^dagger = rho) with Tr(rho) = 1

**The tensor product problem**: For two quaternionic systems H_1 and H_2:

```
(|psi>q) tensor |phi> != |psi> tensor (q|phi>) in general
```

because q does not commute with elements of H_2. The "tensor product" H_1 tensor_H H_2 is not a quaternionic Hilbert space -- it collapses to a real vector space.

Adler's resolution: use the **symplectic embedding** H^n -> C^{2n} via q = a + bj -> (a, b). This maps quaternionic QM to complex QM with an additional antilinear symmetry J (the quaternionic structure map, J^2 = -1).

### 4.3 QRE on quaternionic Hilbert space: indirect approach

Since direct tensor products don't exist in H, we cannot define D(rho_AB || sigma_A tensor sigma_B) directly. But we can define:

**Route A: Symplectic embedding.**
Map H^n -> C^{2n}. The quaternionic density operator rho_H maps to a complex rho_C with the additional constraint [rho_C, J] = 0, where J is the quaternionic structure map. Then:

```
D_H(rho || sigma) := D_C(rho_C || sigma_C)  subject to [rho_C, J] = 0, [sigma_C, J] = 0
```

The constraint [rho_C, J] = 0 is exactly the condition that rho_C is SU(2)-invariant (since J generates an SU(2) action on C^{2n}).

**Key insight**: Quaternionic QRE is complex QRE restricted to the SU(2)-invariant sector. The SU(2) structure emerges not from a derivation within Sigma, but from the restriction of Sigma to quaternionic states.

This reverses the logic: SU(2) is input (through the choice of quaternionic structure) rather than output.

**Route B: Watatani relative entropy for subfactors.**
For a von Neumann algebra M with subfactor N subset M, the relative entropy can be defined as:

```
D(omega|_N || phi|_N) = D(omega || phi) - D(omega|_{N'} || phi|_{N'})
```

where N' is the commutant. If N = M^{SU(2)} (the SU(2)-invariant subalgebra), this is precisely the Level 2 chain rule:

```
D(rho || I/d) = D(rho || G(rho)) + D(G(rho) || I/d)
D(rho || G(rho)) = Sigma_{SU(2)-breaking}
```

The subfactor approach naturally associates SU(2) to the Watatani index [M : N] = 4 (for 2-qubit systems with SU(2) invariance). The index [M : N] counts the number of "directions" that are lost when restricting to the subfactor.

### 4.4 Does quaternionic QRE naturally give SU(2)?

**Answer: Yes, but tautologically.** The quaternionic structure IS the SU(2) structure (unit quaternions = SU(2) = Sp(1)). Defining QRE on quaternionic Hilbert space is equivalent to defining QRE on the SU(2)-invariant sector of complex Hilbert space.

This means:
- **Level 2 is automatic**: any QRE calculation on a 2-qubit system restricted to J-invariant states will have SU(2) structure.
- **Level 3 still fails**: having the SU(2) structure in the state space does not produce SU(2) Yang-Mills dynamics.

**New result (from this analysis)**: The quaternionic Sigma can be defined as:

```
Sigma_H(rho, N) := D(rho_C || N(rho_C))  where rho_C = symplectic embedding, [rho_C, J] = 0
```

This satisfies:
1. Sigma_H >= 0 (DPI)
2. Sigma_H = 0 iff N|_{J-sector} is reversible
3. Sigma_H decomposes into SU(2) irreps exactly as in the Level 2 calculation

But this is a reformulation, not a derivation.

---

## 5. Kaluza-Klein Approach: 7D and SU(2)

### 5.1 Standard 7D Kaluza-Klein for SU(2) x U(1)

The standard KK construction for SU(2) uses a 7-dimensional manifold:

```
M^7 = M^4 x S^3
```

where S^3 = SU(2) is the internal space. The 7D metric decomposes as:

```
ds_7^2 = g_{mu nu}(x) dx^mu dx^nu
       + r^2 [sigma_1^2 + sigma_2^2 + sigma_3^2
              + 2 A^a_mu(x) sigma_a dx^mu + A^a_mu A^b_mu sigma_a sigma_b]
```

where sigma_a (a=1,2,3) are the left-invariant 1-forms on S^3:

```
sigma_1 = cos(psi) d(theta) + sin(psi) sin(theta) d(phi)
sigma_2 = -sin(psi) d(theta) + cos(psi) sin(theta) d(phi)
sigma_3 = d(psi) + cos(theta) d(phi)
```

and r is the radius of S^3. The vacuum Einstein equations R_{MN}^{(7)} = 0 in 7D yield:

1. 4D Einstein equations with SU(2) Yang-Mills stress-energy
2. SU(2) Yang-Mills equations: D_nu F^{a,mu nu} = 0
3. Scalar equation for r (the modulus, analogous to the 5D dilaton)

The SU(2) gauge coupling is:

```
g_2 = sqrt(16 pi G_4 / Vol(S^3)) = sqrt(16 pi G_4 / (2 pi^2 r^3))
```

### 5.2 Sigma in 7D: decomposition

In the Sigma framework, the 7D entropy production should decompose as:

```
Sigma_{7D} = -ln(-G_{00}^{(7)})
```

The 7D (0,0) metric component is:

```
G_{00}^{(7)} = g_{00}(x) + r^2 A^a_0 A^a_0
```

(where sum over a is implied). Therefore:

```
Sigma_{7D} = -ln(-g_{00} - r^2 A^a_0 A^a_0)
           = -ln(-g_{00}) - ln(1 - r^2 A^a_0 A^a_0 / |g_{00}|)
           = Sigma_grav + Sigma_{SU(2)}
```

where:

```
Sigma_{SU(2)} = -ln(1 - r^2 |A_0|^2 / |g_{00}|)
```

with |A_0|^2 = A^1_0 A^1_0 + A^2_0 A^2_0 + A^3_0 A^3_0 (the SU(2)-norm of the temporal gauge potential).

**This is the direct non-abelian generalization of the 5D decomposition from Paper 6 Supplement Section S4.**

### 5.3 Including both U(1) and SU(2): 8D = 4D + S^3 + S^1

For the electroweak sector SU(2) x U(1), we need 8D = 4D + S^3 + S^1:

```
Sigma_{8D} = Sigma_grav + Sigma_{SU(2)} + Sigma_{U(1)}
```

The Weinberg angle enters through the relative radii:

```
sin^2(theta_W) = g_1^2 / (g_1^2 + g_2^2)
               = R_{S^1}^{-2} / (R_{S^1}^{-2} + R_{S^3}^{-3} * [volume factor])
```

In standard KK, the gauge couplings are:

```
g_1^2 = 16 pi G / (2 pi R_1)    (from S^1)
g_2^2 = 16 pi G / (2 pi^2 R_3^3)  (from S^3)
```

Therefore:

```
sin^2(theta_W) = g_1^2/(g_1^2+g_2^2) = (2pi^2 R_3^3) / (2pi^2 R_3^3 + 2pi R_1)
```

This constrains the ratio R_3^3/R_1 but does not involve mu_{Khronon}.

### 5.4 Does 7D KK constrain g_2?

The 7D decomposition gives:

```
Sigma_{SU(2)} = -ln(1 - r^2 |A_0|^2 / |g_{00}|)
```

For weak fields (r |A_0| << sqrt(|g_{00}|)):

```
Sigma_{SU(2)} ~ r^2 |A_0|^2 / |g_{00}|
```

The non-abelian Yang-Mills equations emerge from the 7D Einstein equations. The gauge coupling g_2 is set by the geometry (r = radius of S^3), not by the Sigma framework. Sigma provides the entropy production decomposition; KK geometry provides the gauge coupling.

**Conclusion**: The 7D KK approach gives a clean decomposition Sigma_{7D} = Sigma_grav + Sigma_{SU(2)}, but does not derive g_2 from Sigma. The gauge coupling is a geometric input (radius of S^3), as in standard KK theory.

### 5.5 Sigma as the dynamical principle for KK moduli

There is one potential advantage of the Sigma framework over standard KK: the modulus r (radius of S^3) can be treated as a dynamical field whose value is determined by the Sigma variational principle delta Sigma = 0.

If we postulate that the total Sigma_{7D} is extremized:

```
delta/delta r [Sigma_{7D}(g, A, r)] = 0
```

this could fix r (and hence g_2) dynamically. The equation is:

```
d Sigma_{SU(2)} / dr = 0
=> 2r |A_0|^2 / (|g_{00}| - r^2 |A_0|^2) = 0
```

which gives r = 0 (trivial) or |A_0| = 0 (no gauge field). Neither is physical.

A more sophisticated variational principle would use the full 7D Sigma action:

```
S_{Sigma} = int sqrt(-G^{(7)}) Sigma_{7D} d^7x
```

The Euler-Lagrange equation for r from this action would be a nonlinear PDE coupling the modulus to the gauge field. This could in principle fix g_2, but the calculation has not been done.

### 5.6 Sigma complex Khronon doublet in 7D

The complex Khronon phi = f e^{i theta} from Paper 6 can be generalized to a **quaternionic Khronon**:

```
Phi = f exp(i theta_1 + j theta_2 + k theta_3) = f * q
```

where q in SU(2) = Sp(1) is a unit quaternion. The four real components of Phi are:

```
Phi = f (cos(alpha) + i sin(alpha) cos(beta) + j sin(alpha) sin(beta) cos(gamma) + k sin(alpha) sin(beta) sin(gamma))
```

This gives:
- 1 amplitude f -> gravity (Sigma_grav = 2 ln Q, Q determined by f)
- 3 angles (alpha, beta, gamma) -> SU(2) gauge field (3 components = 3 generators of SU(2))

The 7D KK metric with quaternionic Khronon is:

```
ds_7^2 = e^{2alpha phi} g_{mu nu} dx^mu dx^nu + f^2 (sigma_a + A^a_mu dx^mu)^2
```

where f = |Phi| is the Khronon amplitude and A^a_mu = partial_mu theta^a (in the leading approximation, before gauging).

**This is the quaternionic generalization of Paper 6's complex Khronon.** It provides the geometric structure for SU(2) from the phase degrees of freedom of the Khronon field.

**Critical question**: Does the ghost condensation K'(Q_0) = 0 still work for the quaternionic Khronon? The ghost condensation requires:

```
K(Q) with K'(1) = 0 (condensation at Q=1)
```

For a quaternionic Khronon, Q = |Phi|/N_Phi where N_Phi is the lapse computed from the quaternionic norm. The ghost condensation condition is still K'(Q_0) = 0, but now Q depends on the full quaternionic amplitude, not just the complex one.

The SU(2) gauge field lives in the "angular" part of the quaternionic Khronon, which is orthogonal to the "radial" part that determines gravity. Ghost condensation affects only the radial part. Therefore, **ghost condensation and SU(2) gauge dynamics are decoupled** at leading order.

---

## 6. Three Most Promising Routes to SU(2) from Sigma

### Route A: Bianconi 1-form sector (highest priority)

The Sigma framework = Bianconi's 0-form sector (confirmed 2026-03-19). Opening the 1-form sector introduces gauge fields. The construction:

```
L_{Bianconi} = -Tr_F ln(G_tilde g_tilde^{-1})
```

with Dirac-Kahler matter of degree p:
- p=0: scalar (Khronon) -> gravity -> our Sigma = 2 ln Q
- p=1: vector (gauge field) -> Yang-Mills -> Sigma_{gauge}
- p=2: 2-form (B field) -> higher gauge -> Sigma_{2-form}

The 1-form sector on a 4D manifold with SU(2) gauge group gives:

```
Sigma_{1-form} = -Tr_{SU(2)} ln(1 + F/b_W)
```

where F is the SU(2) field strength and b_W is the weak-force analogue of the Born-Infeld critical field.

In the weak-field limit:

```
Sigma_{1-form} ~ Tr(F^2) / (2 b_W^2) = (1/2) F^a_{mu nu} F^{a,mu nu} / b_W^2
```

which is the Yang-Mills action (up to normalization).

**Todo**: Work out the explicit 1-form Bianconi action with SU(2) structure and verify that the weak-field limit gives Yang-Mills.

### Route B: Quaternionic Khronon + 7D KK (Section 5.6 above)

The quaternionic Khronon Phi = f * q with q in SU(2) provides:
- Amplitude f -> gravity/DM (Papers 2-4)
- Phase q -> SU(2) gauge field (Paper 7)
- The 7D KK decomposition gives Sigma_{7D} = Sigma_grav + Sigma_{SU(2)}
- Modified Yang-Mills equations from the dynamical dilaton (f)

**Todo**: Write the full 7D Lagrangian with quaternionic Khronon, derive the 4D equations, and verify that SU(2) Yang-Mills emerges.

### Route C: Harlow QEC + Petz (longest term, most fundamental)

The Harlow framework identifies gauge symmetry = QEC. The Sigma framework provides the Petz recovery map. If we can show:

1. The optimal recovery channel for the SU(2) twirl is the Petz map
2. The Petz recovery fidelity gives F = e^{-Sigma_{SU(2)}/2}
3. The field equations from delta Sigma_{SU(2)} = 0 reproduce Yang-Mills

then SU(2) would follow from the Sigma framework at Level 3.

**The key difficulty**: Step 3. The variational equation delta Sigma = 0 gives equations for the state rho, not for the gauge connection A. One needs to show that the optimal state rho(A) depends on A in such a way that the state-space variational equation implies the configuration-space (gauge field) equation.

---

## 7. New Observations and Conjectures

### 7.1 D_break = ln 3 is the Donnelly-Wall edge mode entropy (NEW)

The Level 2 result D(rho || G_{SU(2)}(rho)) = ln 3 for the triplet coincides with the Donnelly-Wall edge mode contribution to entanglement entropy for SU(2) gauge fields. This suggests that the QRE decomposition under the twirl is computing the gauge-theoretic entanglement entropy, not just a generic Lie group decomposition.

If confirmed, this means: **Sigma_{weak} = edge mode entanglement entropy of the weak gauge field**.

### 7.2 Quaternionic Khronon conjecture (NEW)

**Conjecture**: The Khronon field in the full theory is quaternionic: Phi in H. The complex Khronon phi = f e^{i theta} of Paper 6 is a truncation to the C subset H sector. The full quaternionic Khronon gives:

```
amplitude |Phi| = f  -> gravity/DM  (Papers 2-4)
U(1) phase theta_1   -> EM          (Paper 6)
SU(2) phases theta_2, theta_3 -> weak force (Paper 7)
```

The gauge group SU(2) x U(1) emerges from the phase structure of the quaternionic Khronon, analogous to how U(1) emerged from the complex Khronon.

**Test**: Show that the ghost condensation K'(Q_0) = 0 for the quaternionic Khronon admits SU(2) fluctuations around the condensate.

### 7.3 The Weinberg angle from Khronon doublet geometry (SPECULATIVE)

If the quaternionic Khronon has the structure Phi = (phi_1, phi_2) in C^2 = H (doublet), the Weinberg angle could be related to the angle between the U(1) and SU(2) sectors in quaternionic space:

```
phi_1 = f_1 e^{i theta}     (U(1) sector)
phi_2 = f_2 e^{i alpha}     (SU(2) sector, simplified)
```

The mixing angle between these sectors at the condensation point could give:

```
sin^2(theta_W) = |phi_2|^2 / (|phi_1|^2 + |phi_2|^2) = f_2^2 / (f_1^2 + f_2^2)
```

At the GUT scale (maximum symmetry), f_1 = f_2 -> sin^2(theta_W) = 1/2 (wrong, need 3/8).

With the Connes normalization (5/3 hypercharge factor): sin^2(theta_W) = 3/8 (correct).

**This is highly speculative and likely wrong in detail.** The 3/8 prediction comes from the representation theory of A_F, not from a simple amplitude ratio. But it illustrates how the quaternionic Khronon could in principle encode the Weinberg angle.

---

## 8. Summary: What Is Derived, What Is Assumed, What Is Open

| Result | Status | Route |
|--------|--------|-------|
| SU(2) structure in 2-qubit QRE (Level 2) | VERIFIED (10^{-16}) | Chain rule |
| D_break = ln 3 = dim Im(H) for triplet | DERIVED | Marvian-Spekkens |
| D_break = ln 3 = Donnelly-Wall edge entropy | NEW OBSERVATION | This analysis |
| SU(2) Yang-Mills from Sigma (Level 3) | NOT YET DERIVED | -- |
| Quaternionic Khronon -> SU(2) gauge field | CONJECTURE | This analysis |
| Sigma_{7D} = Sigma_grav + Sigma_{SU(2)} | DERIVED (KK geometry) | Section 5 |
| sin^2(theta_W) from Sigma | NOT DERIVED | -- |
| mu_{Khronon} constrained by theta_W | NO (different mu's) | Section 3 |
| Bianconi 1-form -> SU(2) Yang-Mills | PLAUSIBLE (not calculated) | Route A |
| Harlow QEC + Petz -> SU(2) gauge | PLAUSIBLE (not calculated) | Route C |

---

## 9. Concrete Next Steps for Paper 7

### Week 1-2: Bianconi 1-form sector
- Read Bianconi (2025 PRD) in full
- Write out L = -Tr_F ln(G_tilde g_tilde^{-1}) for p=1 (1-form) with SU(2) indices
- Verify weak-field limit gives Yang-Mills action
- Compute Sigma_{SU(2)} from the 1-form action

### Week 3-4: Quaternionic Khronon Lagrangian
- Write the full Lagrangian for Phi in H with ghost condensation
- Verify that SU(2) fluctuations are allowed around the condensate
- Derive the 4D equations from 7D KK
- Check: does the dilaton coupling modify Yang-Mills?

### Week 5-6: Donnelly-Wall connection
- Compute the edge mode entanglement entropy for SU(2) gauge field on a ball
- Compare with D(rho || G(rho)) = ln 3
- If they match: write up the connection as a theorem

### Week 7-8: Go/No-Go
- Can we derive SU(2) Yang-Mills from any combination of Routes A-C?
- If yes: write Paper 7
- If no: identify the precise mathematical obstruction and document it

### Key literature to read
1. Harlow (2018): "TASI Lectures on the Emergence of Bulk Physics in AdS/CFT" (arXiv:1802.01040)
2. Harlow (2015): "The Ryu-Takayanagi Formula from Quantum Error Correction" (arXiv:1607.03901)
3. Donnelly-Wall (2015): "Geometric entropy and edge modes of the electromagnetic field" (arXiv:1506.04267)
4. Bianconi (2025): PRD 111, 066001 -- "Gravity from Entropy"
5. Adler (1995): "Quaternionic Quantum Mechanics and Quantum Fields" (Oxford UP)
6. Marvian-Spekkens (2014): "Extending Noether's theorem..." (arXiv:1404.3236)
7. Casini-Huerta-Rosabal (2014): "Remarks on entanglement entropy for gauge fields" (arXiv:1312.1183)

---

## 10. Attribution

| Discovery | Who | Ours? |
|-----------|-----|-------|
| SU(2) = unit quaternions | Hamilton (1843) | No |
| KK for SU(2): 7D with S^3 fiber | Witten (1981), Duff et al. | No |
| QRE chain rule under twirl | Marvian-Spekkens (2014) | No |
| Edge modes for gauge fields | Donnelly-Wall (2015) | No |
| Gauge = QEC | Harlow (2018) | No |
| Bianconi entropy → gravity | Bianconi (2025) | No |
| Szangolies qubit → gauge hierarchy | Szangolies (2025) | No |
| **D_break = ln 3 = Donnelly-Wall edge entropy** | **New (this analysis)** | **Yes** |
| **Quaternionic Khronon conjecture** | **New (this analysis)** | **Yes** |
| **Sigma_{7D} decomposition for SU(2)** | **New (explicit calculation)** | **Yes** |
| **mu_{Khronon} ≠ mu_{RG} (non-constraint)** | **New (explicit check)** | **Yes** |
| **Three-route strategy for Paper 7** | **New (this analysis)** | **Yes** |

---

## Record
Created 2026-03-19 during investigation of SU(2) from Sigma framework.
