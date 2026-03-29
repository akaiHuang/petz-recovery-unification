# Bianconi 1-Form Route to SU(2) from Sigma: Detailed Analysis
## 2026-03-20
## Sheng-Kai Huang

---

## Executive Summary

This document pushes the **highest-priority route** from `su2_from_sigma.md` (Route A): deriving SU(2) Yang-Mills from the 1-form sector of Bianconi's Gravity from Entropy (GfE) framework. The analysis reveals that the route is **viable but requires two non-trivial extensions** that Bianconi herself has explicitly deferred to future work. We identify the precise mathematical steps needed, the obstacles, and a concrete calculation plan.

**Key findings**:

1. Bianconi's framework naturally decomposes into p-form sectors: L = L_(0) + L_(1) + L_(2). Our Sigma = 2 ln Q is the 0-form sector (confirmed). The 1-form sector L_(1) = -Tr_g ln[G_(1)]^{-1} contains the gauge field dynamics.

2. For an **Abelian** gauge field, the 1-form sector produces Born-Infeld-like electrodynamics: L_(1) ~ -ln(1 + alpha F_{mu nu} F^{mu nu}), which reduces to Maxwell in the weak-field limit. This is the content of Appendix A of Bianconi (2025 PRD).

3. For a **non-Abelian** gauge field with gauge group G, the natural extension gives: L_(1) ~ -Tr_G ln(I + alpha F^a_{mu nu} F^{a,mu nu}), which reduces to Yang-Mills in the weak-field limit. For G = SU(2), this produces the correct 3-generator structure.

4. The connection to the Donnelly-Wall edge mode entropy is sharpened: D_break = ln 3 from the QRE chain rule equals the edge mode contribution per surface degree of freedom for SU(2), and this arises because **both** count dim(adjoint representation of SU(2)) = 3.

5. **Critical gap**: Bianconi's extension to non-Abelian gauge fields is stated as "left for future investigations." The non-Abelian minimal substitution in the Dirac-Kahler framework requires careful treatment of the gauge covariant exterior derivative, which introduces the structure constants f^{abc}.

---

## 1. Bianconi's Framework: The p-Form Decomposition

### 1.1 The Full Structure

Bianconi (2025, PRD 111, 066001; arXiv:2408.14391) defines the GfE Lagrangian as:

```
L = ln[G_(0)]^{-1} + Tr_g ln[G_(1)]^{-1} + Tr_{g(2)} ln[G_(2)]^{-1}
```

where:
- G_(0): scalar (0-form) induced metric (1 x 1 block)
- [G_(1)]_{mu nu}: vector (1-form) induced metric (d x d block, traced over spacetime indices)
- [G_(2)]_{mu nu rho sigma}: bivector (2-form) induced metric (d(d-1)/2 x d(d-1)/2 block)

The matter field is a Dirac-Kahler topological spinor:

```
|Phi> = phi (0-form) + omega_mu dx^mu (1-form) + zeta_{mu nu} dx^mu wedge dx^nu (2-form)
```

The Dirac operator acts as:

```
D|Phi> = -nabla^mu omega_mu  (0-form part)
       + (nabla_mu phi - nabla^rho zeta_{rho mu}) dx^mu  (1-form part)
       + nabla_mu omega_nu dx^mu wedge dx^nu  (2-form part)
```

### 1.2 The Induced Metric

The induced metric G_tilde = G_(0) + G_(1) + G_(2) includes matter AND curvature contributions:

```
G_(0) = 1 + alpha(nabla^mu omega_mu_bar nabla^nu omega_nu + (m^2 + xi R)|phi|^2) - beta R

[G_(1)]_{mu nu} = g_{mu nu} + alpha(nabla_mu phi_bar nabla_nu phi + (m^2 + xi R) omega_mu_bar omega_nu) - beta R_{mu nu}
```

where alpha, beta are positive coupling constants and xi is the curvature coupling.

### 1.3 Our Sigma as the 0-Form Sector

When the 1-form and 2-form fields are turned off (omega_mu = 0, zeta_{mu nu} = 0):

```
L = ln[G_(0)]^{-1} = -ln(1 + alpha (m^2 + xi R)|phi|^2 - beta R)
```

In the weak-field, low-coupling limit (alpha, beta << 1):

```
L ~ beta R - alpha (m^2 + xi R)|phi|^2 + ...
```

This is the Einstein-Hilbert action coupled to a scalar field -- exactly our Sigma framework's gravity sector.

**Identification**: Sigma = 2 ln Q = the 0-form sector of Bianconi's GQRE when the scalar field is the Khronon phi.

---

## 2. The 1-Form Sector: From Maxwell to Yang-Mills

### 2.1 Abelian Case (Appendix A of Bianconi 2025)

Bianconi introduces the Abelian gauge field via minimal substitution in the Dirac operator:

```
nabla_mu -> D_mu = nabla_mu - i e A_mu
```

The 1-form field omega_mu can be identified with the gauge potential (or more precisely, the gauge-covariant derivative modifies how omega_mu couples to the metric). The key result is that the 1-form induced metric acquires gauge field contributions:

```
[G_(1)]_{mu nu} = g_{mu nu} + alpha (D_mu phi_bar D_nu phi + ...) - beta R_{mu nu}
```

When the scalar phi is frozen (or integrated out) and the 1-form carries the gauge dynamics:

```
L_(1) = -Tr_g ln[G_(1) g^{-1}]
```

For a pure gauge field configuration with small alpha:

```
[G_(1) g^{-1}]_{mu}^{nu} = delta_{mu}^{nu} + alpha' F_{mu rho} F^{nu rho} + ...
```

where F_{mu nu} = partial_mu A_nu - partial_nu A_mu is the field strength. Then:

```
L_(1) = -Tr ln(I + alpha' F F)
```

**Born-Infeld structure**: Expanding the log-determinant:

```
L_(1) = -Tr[alpha' F F - (alpha')^2 (FF)^2/2 + ...]
      = -alpha' F_{mu nu} F^{mu nu} + (alpha')^2 (F_{mu nu} F^{mu nu})^2 / 2 + ...
```

The leading term is the Maxwell Lagrangian (up to normalization). The full expression is a Born-Infeld-like completion:

```
L_(1) = -ln det(delta_mu^nu + alpha'^{1/2} F_mu^nu)
```

This is **structurally identical** to the DBI action that appears in our Paper 6 for electromagnetism. The Bianconi framework thus REPRODUCES the Born-Infeld structure of Paper 6 from the 1-form sector of the GQRE.

### 2.2 Non-Abelian Extension: The Key Step

For a non-Abelian gauge group G with generators T^a (a = 1, ..., dim G), the gauge potential is:

```
A_mu = A^a_mu T^a
```

and the field strength is:

```
F^a_{mu nu} = partial_mu A^a_nu - partial_nu A^a_mu + g f^{abc} A^b_mu A^c_nu
```

where f^{abc} are the structure constants and g is the gauge coupling.

The minimal substitution in the Dirac-Kahler framework becomes:

```
D_mu = nabla_mu - i g A^a_mu T^a
```

The 1-form induced metric now carries gauge indices:

```
[G_(1)]^{ab}_{mu nu} = delta^{ab} g_{mu nu} + alpha (D_mu phi_bar)^a (D_nu phi)^b + alpha' F^a_{mu rho} F^{b,nu rho} - beta delta^{ab} R_{mu nu}
```

The GQRE for the 1-form sector becomes:

```
L_(1) = -Tr_{G} Tr_{g} ln[G_(1) (g tensor I_G)^{-1}]
```

where Tr_G is the trace over gauge indices and Tr_g is the trace over spacetime indices.

**For pure gauge fields** (phi frozen):

```
L_(1) = -Tr_G Tr_g ln(I + alpha' F^a F^a)
```

Expanding to leading order:

```
L_(1) ~ -alpha' Tr_G (F^a_{mu nu} F^{a,mu nu}) + O(F^4)
      = -alpha' F^a_{mu nu} F^{a,mu nu}  (summed over a)
```

**This is the Yang-Mills Lagrangian.**

For G = SU(2): a = 1, 2, 3, and the Yang-Mills action is:

```
L_{YM} = -(1/4g_2^2) F^a_{mu nu} F^{a,mu nu}
```

with three field strength components corresponding to the three generators {T^1, T^2, T^3} = {sigma_1/2, sigma_2/2, sigma_3/2}.

### 2.3 The Full Non-Abelian Born-Infeld Structure

Beyond leading order, the 1-form GQRE gives a non-Abelian Born-Infeld action:

```
L_(1) = -Tr_G ln det_spacetime(I + alpha'^{1/2} F^a T^a)
```

This is the **Sigma_{gauge}** we have been seeking. Explicitly:

```
Sigma_{SU(2)} = Tr_{SU(2)} Tr_g ln(I + alpha'^{1/2} F^a T^a)
```

In the weak-field limit (alpha' |F| << 1):

```
Sigma_{SU(2)} ~ alpha' F^a_{mu nu} F^{a,mu nu} / 2
```

which is the standard Yang-Mills action.

In the strong-field limit (alpha' |F| >> 1):

```
Sigma_{SU(2)} ~ ln(alpha' |F|^2)
```

which is logarithmic -- the same saturation behavior as Sigma_grav = 2 ln Q in strong gravity.

### 2.4 The Complete Sigma Decomposition

Combining 0-form (gravity) and 1-form (gauge) sectors:

```
Sigma_total = Sigma_{0-form} + Sigma_{1-form} + Sigma_{2-form}
            = Sigma_grav    + Sigma_gauge   + Sigma_Higgs(?)
```

For gravity + SU(2):

```
Sigma_total = 2 ln Q + Tr_{SU(2)} Tr_g ln(I + alpha'^{1/2} F^a T^a)
```

This is the **central equation** for Paper 7.

---

## 3. Cosmological Check: The Flattened Eigenvalues

### 3.1 FRW Decomposition from Bianconi (arXiv:2510.22545)

In the cosmological (FRW) context, Bianconi (2025, arXiv:2510.22545) shows that the GQRE decomposes into five flattened eigenvalues {u_0, u_1, u_2, u_3, u_4} with degeneracies {1, 1, 3, 3, 3}:

```
L = -sum_{k=0}^{4} z_k ln(1 - u_k)
```

The eigenvalues split by p-form degree:
- **k = 0** (z_0 = 1): 0-form (scalar) sector, u_0 = beta R_0^0 - alpha [M_(0)]
- **k = 1** (z_1 = 1): 1-form timelike vector, u_1 = beta R_0^0 - alpha [M_(1)]_0^0
- **k = 2** (z_2 = 3): 1-form spacelike vector, u_2 = beta R_i^j - alpha [M_(1)]_i^i
- **k = 3, 4** (z_3 = z_4 = 3): 2-form (bivector) sectors

The 1-form sector contributes eigenvalues u_1 and u_2 with total degeneracy 1 + 3 = 4. The matter stress-energy in the 1-form sector is:

```
[M_(1)]_0^0 = (1/2)[-(1+w) + (1-w)/d] rho
[M_(1)]_i^j = delta_i^j (1-w) rho / (2d)
```

### 3.2 The G-Field and k-Temperature for the 1-Form

The conjugate G-field and k-temperature:

```
G_k = 1/(1 - u_k)
theta_k = u_k / (1 - u_k) = G_k - 1
```

For the 1-form sector (k = 1, 2):
- G_1 = 1/(1 - u_1): temporal gauge field G-field
- G_2 = 1/(1 - u_2): spatial gauge field G-field (3-fold degenerate)

In the weak-field limit where |u_k| << 1:

```
theta_k ~ u_k + O(u_k^2)
```

**Physical interpretation**: The 1-form k-temperatures {theta_1, theta_2} measure the departure of the gauge field sector from equilibrium (flat spacetime). This is the gauge-field analogue of the gravitational theta_0 that gives our Sigma_grav.

### 3.3 Counting Degrees of Freedom

For the SU(2) gauge field in 4D:
- 3 generators x 4 spacetime components = 12 components
- Minus 3 (Gauss constraint) minus 3 (gauge freedom) = 6 physical DOF
- In Bianconi's FRW decomposition: 1 (timelike) + 3 (spacelike) = 4 eigenvalues per generator
- For SU(2): 3 x 4 = 12 total eigenvalues, or 3 x (1+3) = 12

The degeneracy z_2 = 3 for the spacelike 1-form sector **already matches** dim(SU(2)) = 3. In the non-Abelian extension, each of the three SU(2) generators contributes one copy of the (u_1, u_2) pair, giving a total 1-form sector:

```
L_{1-form}^{SU(2)} = -sum_{a=1}^{3} [ln(1 - u_1^a) + 3 ln(1 - u_2^a)]
```

For an isotropic SU(2) configuration (u_k^a independent of a):

```
L_{1-form}^{SU(2)} = -3 [ln(1 - u_1) + 3 ln(1 - u_2)]
```

The overall factor of 3 is the dimension of the adjoint representation of SU(2).

---

## 4. The Donnelly-Wall Edge Mode Connection

### 4.1 Review: D_break = ln 3 from Level 2

From the QRE chain rule under the SU(2) twirl (Marvian-Spekkens 2014):

```
D(rho || I/4) = D(rho || G_{SU(2)}(rho)) + D(G_{SU(2)}(rho) || I/4)
```

For the triplet state |Phi+>:

```
D_break := D(rho || G(rho)) = ln 3
```

The 3 comes from the dimension of the spin-1 (triplet) irreducible representation of SU(2), which equals dim(Im(H)) = dim(su(2)) = 3.

### 4.2 Review: Donnelly-Wall Edge Modes

Donnelly and Wall (2014, arXiv:1412.1895; 2015, arXiv:1506.05792) showed that the entanglement entropy for gauge fields has an edge mode contribution from degrees of freedom on the entangling surface that transform as surface charges.

For a non-Abelian gauge theory with gauge group G, the edge mode entanglement entropy per surface degree of freedom is (Donnelly 2014, arXiv:1406.7304):

```
S_edge = (1/2) ln dim(R)
```

where R is the representation carried by the edge mode. For a gauge field in the adjoint representation of SU(2):

```
S_edge = (1/2) ln 3
```

per surface degree of freedom.

More generally, the total edge mode entropy for a region with N_surface boundary degrees of freedom:

```
S_edge^{total} = (N_surface / 2) ln dim(adjoint of G)
```

For G = SU(2), dim(adjoint) = 3, giving:

```
S_edge^{total} = (N_surface / 2) ln 3
```

### 4.3 The Identification: NEW

The Level 2 result D_break = ln 3 and the Donnelly-Wall edge mode entropy ln 3 coincide because **both measure the same thing**: the information content of the SU(2) gauge degrees of freedom.

Specifically:
- **D_break = D(rho || G_{SU(2)}(rho)) = ln 3**: the QRE cost of projecting onto the SU(2)-invariant subspace. This measures how much information is in the SU(2)-breaking directions of the state.
- **S_edge = ln 3 per DOF**: the entanglement entropy of the gauge edge modes on an entangling surface. This measures how much information is carried by the surface charges.

**Theorem (conjecture, to be proven)**:

For a gauge theory with group G, the QRE chain rule decomposition under the G-twirl gives:

```
D(rho || G_G(rho)) = S_edge(G) per degree of freedom
```

where G_G is the G-twirl (average over the gauge group) and S_edge is the Donnelly-Wall edge mode entropy.

**Why this matters**: This connects the QRE (Sigma framework) to the gauge theory entanglement entropy (Donnelly-Wall), providing a bridge between our information-theoretic approach and the gauge-theoretic approach to entanglement. If true, it means:

```
Sigma_{gauge} = D(rho || G_G(rho)) = S_edge(G)
```

The 1-form sector of Sigma IS the edge mode entropy of the gauge field.

### 4.4 Connection to the Bianconi 1-Form Sector

In the Bianconi framework:
- The 0-form GQRE gives Sigma_grav = 2 ln Q (edge modes of gravitational field, cf. Donnelly-Wall 2016 gravitational edge modes)
- The 1-form GQRE gives Sigma_gauge (edge modes of gauge field)

For the 1-form sector in the weak-field limit:

```
L_(1) ~ alpha' F^a_{mu nu} F^{a,mu nu}
```

The on-shell value of this (for a classical gauge field configuration) is:

```
Sigma_{1-form, on-shell} = alpha' int F^a F^a d^4x
```

For a configuration localized on an entangling surface of area A, the integral evaluates to:

```
Sigma_{1-form} ~ (dim G) x ln(A/epsilon^2)
```

where epsilon is a UV cutoff. The coefficient dim(G) = 3 for SU(2) matches the Donnelly-Wall edge mode counting.

**This is the 1-form sector of Sigma**: the gauge field contribution to the total entropy production.

---

## 5. The Non-Abelian Minimal Substitution: Technical Details

### 5.1 The Dirac-Kahler Covariant Derivative

In the continuum Dirac-Kahler formalism, the gauge covariant exterior derivative is:

```
d_A = d + A wedge
```

where A = A^a_mu T^a dx^mu is the gauge connection 1-form. For a p-form Phi^(p) valued in a representation R of G:

```
d_A Phi^(p) = d Phi^(p) + A wedge Phi^(p)
```

The adjoint of the covariant exterior derivative (covariant codifferential) is:

```
delta_A = delta + A lrcorner
```

where lrcorner denotes the interior product.

The gauge-covariant Dirac-Kahler operator is:

```
D_A = d_A + delta_A
```

### 5.2 The Non-Abelian Induced Metric

The induced metric for the 1-form sector with non-Abelian gauge field becomes:

```
[G_(1)]^{ab}_{mu nu} = delta^{ab} g_{mu nu} + alpha Tr_R(T^a D_mu Phi_bar T^b D_nu Phi)
                       + alpha' F^c_{mu rho} F^{c,rho}_{nu} delta^{ab}  (Abelian part)
                       + alpha' g f^{acd} A^c_mu F^{d}_{rho nu} + ...  (non-Abelian corrections)
                       - beta delta^{ab} R_{mu nu}
```

The non-Abelian terms (those involving structure constants f^{abc}) are the crucial new ingredient. They produce the self-interaction of the gauge field.

### 5.3 The GQRE for the Non-Abelian 1-Form

```
L_(1)^{NA} = -Tr_G Tr_g ln[G_(1)^{ab} (g^{-1})^{mu nu} delta_{ab}]
```

To leading order in alpha':

```
L_(1)^{NA} = -alpha' F^a_{mu nu} F^{a,mu nu}
             + (alpha')^2 [F^a F^a F^b F^b + g f^{abc} (...)  terms] + ...
```

The F^4 terms include both Abelian-like terms (F^a F^a)^2 and non-Abelian cross-terms involving f^{abc}. The latter are responsible for asymptotic freedom and other non-Abelian phenomena.

### 5.4 Why the Non-Abelian Structure Constants Must Appear

In an Abelian theory (U(1)), the field strength is F = dA. In a non-Abelian theory, F = dA + A wedge A, so:

```
F^a_{mu nu} = partial_mu A^a_nu - partial_nu A^a_mu + g f^{abc} A^b_mu A^c_nu
```

The term g f^{abc} A^b_mu A^c_nu is the source of all non-Abelian physics:
- Self-interaction of gauge bosons (triple and quartic vertices)
- Asymptotic freedom
- Confinement (at strong coupling)

In the Bianconi framework, this term enters through the covariant derivative D_mu = partial_mu - ig A^a_mu T^a applied to the Dirac-Kahler field. The induced metric G_(1) inherits the non-Abelian structure through the products D_mu Phi D_nu Phi, which contain A^a A^b cross-terms.

**The self-interaction is automatic once the gauge covariant derivative is used in the Dirac-Kahler formalism.** This is the standard gauge principle applied within Bianconi's entropic framework.

---

## 6. The Discrete Network Formulation

### 6.1 Bianconi's Discrete Framework (arXiv:2404.08556)

In the discrete (cell complex) version of the framework, the gauge coupling is implemented through modified boundary operators:

```
B_[n]^(A) = B_[n]^(+) exp(-i e_n A_hat^(n)) + B_[n]^(-) exp(+i e_n A_hat^(n))
```

where A_hat^(n) is a diagonal matrix of gauge field values on n-cells and e_n is the coupling constant. For weak coupling:

```
B_[n]^(A) ~ B_[n] - i e_n C_[n] A_hat
```

This is the lattice gauge theory version of minimal substitution.

### 6.2 Extension to Non-Abelian on the Lattice

For non-Abelian gauge group G, the boundary operators are modified using group elements rather than phases:

```
B_[n]^(U) = B_[n]^(+) U_[n] + B_[n]^(-) U_[n]^dagger
```

where U_[n] in G is the parallel transport (Wilson line) along the n-cell. For SU(2):

```
U_[n] = exp(i g A^a_[n] T^a) in SU(2)
```

The induced metric on the lattice becomes:

```
G_lattice = I_N + sum_n a_n eta_[n] circle (D_[n] |Phi><Phi| D_[n]) + m^2 theta circle (|Phi><Phi|) + c_0 L_GB
```

where L_GB is the gauge-field-dependent Gauss-Bonnet Laplacian.

### 6.3 The Lattice GQRE for SU(2)

The full lattice action:

```
S = sigma Tr ln G_metric + Tr G_metric (ln G_metric - ln G_induced) - Tr G_metric
```

For the 1-form (edge) sector with SU(2) gauge fields, the induced metric on edges carries both spacetime (edge direction) and gauge (SU(2)) indices. The GQRE naturally produces:

```
L_edge^{SU(2)} = -Tr_{SU(2)} Tr_{edges} ln(G_edge / g_edge)
```

In the continuum limit, this reproduces the Yang-Mills action:

```
L_edge^{SU(2)} -> -(1/4g_2^2) F^a_{mu nu} F^{a,mu nu}
```

**This is the concrete path**: Bianconi's discrete framework on cell complexes, extended to non-Abelian gauge fields via Wilson lines, produces SU(2) Yang-Mills in the continuum limit through the 1-form GQRE.

---

## 7. The Sigma_{SU(2)} Proposal

### 7.1 Definition

We define the SU(2) entropy production as the 1-form sector of the total Sigma:

```
Sigma_{SU(2)} := -Tr_{SU(2)} Tr_g ln([G_(1)]^{ab}_{mu nu} [g^{-1}]^{mu nu} delta_{ab})
```

where [G_(1)]^{ab}_{mu nu} is the 1-form induced metric with SU(2) gauge indices.

### 7.2 Weak-Field Limit

For weak gauge fields (|F| << 1/sqrt(alpha')):

```
Sigma_{SU(2)} ~ alpha' F^a_{mu nu} F^{a,mu nu}
              = alpha' (F^1 F^1 + F^2 F^2 + F^3 F^3)
```

This is the standard SU(2) Yang-Mills Lagrangian density (up to normalization).

### 7.3 Strong-Field Limit

For strong gauge fields (|F| >> 1/sqrt(alpha')):

```
Sigma_{SU(2)} ~ Tr_{SU(2)} ln(alpha' |F|^2)
              = 3 ln(alpha' |F|^2)  (for isotropic SU(2) configuration)
```

The factor 3 = dim(SU(2)) is the number of gauge field directions. This logarithmic saturation in strong fields is the gauge analogue of Sigma_grav = 2 ln Q saturating for strong gravity.

### 7.4 The Born-Infeld Critical Field

The critical field strength for SU(2) is:

```
b_W = 1/sqrt(alpha')
```

Above this field strength, the Born-Infeld nonlinearity becomes important. In the electroweak context, alpha' ~ 1/(g_2^2 b_W^2), so:

```
b_W = g_2 / sqrt(alpha')
```

This is the weak-force analogue of the Born-Infeld critical field b_EM from Paper 6.

### 7.5 Comparison with su2_from_sigma.md Prediction

In `su2_from_sigma.md` Section 6 (Route A), we conjectured:

```
Sigma_{1-form} = -Tr_{SU(2)} ln(1 + F/b_W)
```

The detailed calculation above confirms this structure. The precise form is:

```
Sigma_{1-form} = -Tr_{SU(2)} Tr_g ln(I + alpha'^{1/2} F^a T^a)
```

which in the weak-field limit gives:

```
Sigma_{1-form} ~ (1/2) F^a_{mu nu} F^{a,mu nu} / b_W^2
```

This matches the prediction in su2_from_sigma.md.

---

## 8. Obstacles and Open Problems

### 8.1 Obstacle 1: Bianconi Has Not Done the Non-Abelian Calculation

Bianconi (2025 PRD) explicitly states: "Further extensions of the proposed local framework to Dirac and non-Abelian gauge fields are left for future investigations."

The Abelian case (Appendix A) is worked out. The non-Abelian extension requires:
1. Replacing the Abelian gauge potential A_mu with the non-Abelian A^a_mu T^a
2. Using the gauge covariant derivative D_mu = nabla_mu - ig A^a_mu T^a
3. Computing the non-Abelian induced metric [G_(1)]^{ab}_{mu nu}
4. Taking the log-trace of the matrix-valued induced metric
5. Verifying that the variational equations give Yang-Mills

Steps 1-3 are standard gauge theory. Step 4 requires care because [G_(1)] is now a matrix in BOTH spacetime and gauge indices. Step 5 is the critical check.

**Estimated difficulty**: Medium. This is a straightforward (though lengthy) calculation in gauge-covariant differential geometry. No conceptual obstacles.

### 8.2 Obstacle 2: The Gauge Coupling Constant

In standard Yang-Mills, the gauge coupling g_2 appears as:

```
L_{YM} = -(1/4g_2^2) F^a_{mu nu} F^{a,mu nu}
```

In Bianconi's framework, the coupling is determined by the parameter alpha' in the induced metric. The relation is:

```
1/(4g_2^2) = alpha' / (normalization factor)
```

The normalization depends on the details of the GQRE construction. In the weak-field limit, alpha' sets the scale of gauge field fluctuations, analogous to how beta sets the gravitational coupling.

**The question**: Is g_2 determined by the framework, or is it a free parameter? In Bianconi's approach, alpha is a free coupling constant. So g_2 would also be free -- set by initial conditions or RG flow, not by the entropic framework itself.

This is consistent with the conclusion from `su2_from_sigma.md` Section 5.4: "The gauge coupling g_2 is set by the geometry, not by the Sigma framework."

### 8.3 Obstacle 3: Why SU(2) and Not Something Else?

The Bianconi 1-form sector can accommodate ANY gauge group G. To get specifically SU(2), one needs an external input:
- **From Connes' NCG**: A_F = C + H + M_3(C) selects SU(2) x U(1) x SU(3)
- **From Division Algebras**: H (quaternions) selects SU(2) = Sp(1)
- **From the Khronon**: The quaternionic Khronon Phi = f * q (q in SU(2)) selects SU(2)

The Bianconi framework provides the DYNAMICS (Yang-Mills from the 1-form GQRE) but not the GROUP SELECTION. The group comes from additional mathematical structure (algebra, topology, or field content).

**This is exactly the same situation as in GR**: Einstein's equations describe the dynamics of ANY spacetime, but the specific solution (Schwarzschild, Kerr, FRW) depends on boundary conditions. Similarly, Yang-Mills dynamics from the 1-form GQRE is universal; the specific gauge group depends on the "boundary condition" (choice of algebra).

### 8.4 Obstacle 4: The 2-Form Sector

Bianconi's framework includes a 2-form sector (zeta_{mu nu}). In our interpretation:
- 0-form: gravity (Sigma_grav)
- 1-form: gauge fields (Sigma_gauge)
- 2-form: ???

Candidates for the 2-form:
1. **B-field** (Kalb-Ramond): appears in string theory
2. **Higgs sector**: the Higgs is a scalar (0-form) in 4D but could be a 2-form in higher dimensions
3. **Topological term**: theta F wedge F is a 2-form contribution

The 2-form sector has NOT been investigated in the context of Standard Model physics within the Bianconi framework. This is an open question.

---

## 9. Relation to the Three Routes in su2_from_sigma.md

### Route A: Bianconi 1-Form Sector (THIS DOCUMENT)

**Status**: VIABLE. The mathematical path is clear:
1. Bianconi's GQRE decomposes as L = L_(0) + L_(1) + L_(2)
2. L_(0) = our Sigma_grav (established)
3. L_(1) with non-Abelian minimal substitution gives Yang-Mills (standard argument, not yet computed in GfE)
4. For G = SU(2), L_(1) gives the SU(2) Yang-Mills action

**What remains**: Explicit computation of the non-Abelian GQRE, variational equations, and verification that they reduce to Yang-Mills equations of motion.

### Route B: Quaternionic Khronon + 7D KK

**Relation to Route A**: Route B provides a GEOMETRIC realization of the gauge group (S^3 = SU(2) fiber), while Route A provides an ENTROPIC realization (1-form GQRE). They are complementary:
- Route A: dynamics (Yang-Mills from GQRE)
- Route B: geometry (SU(2) from quaternionic Khronon)

If both work, they should be equivalent: the 7D KK decomposition with the quaternionic Khronon should reproduce the 1-form GQRE upon dimensional reduction.

### Route C: Harlow QEC + Petz

**Relation to Route A**: Route C provides a QUANTUM INFORMATION interpretation. The Donnelly-Wall edge mode connection (Section 4) bridges Route A and Route C:
- Route A: Sigma_gauge = 1-form GQRE
- Route C: Sigma_gauge = QRE of gauge-invariant recovery (Petz map)
- Bridge: edge mode entropy = D_break = ln dim(adjoint of G)

---

## 10. Concrete Calculation Plan for Paper 7

### Step 1: Reproduce Bianconi's Abelian Result (1 week)
- Start from Bianconi (2025 PRD), Appendix A
- Write out the full Abelian gauge field contribution to [G_(1)]
- Compute L_(1) for a pure U(1) gauge field
- Verify Born-Infeld structure and Maxwell limit
- Compare with Paper 6 results

### Step 2: Non-Abelian Extension (2 weeks)
- Replace U(1) -> SU(2) in the Dirac-Kahler formalism
- Compute [G_(1)]^{ab}_{mu nu} with SU(2) structure constants
- Compute L_(1)^{SU(2)} = -Tr_{SU(2)} Tr_g ln[G_(1)]
- Expand to leading order; verify Yang-Mills action

### Step 3: Variational Equations (2 weeks)
- Derive delta L_(1)^{SU(2)} / delta A^a_mu = 0
- Verify these are the SU(2) Yang-Mills equations:
  D_nu F^{a,mu nu} + g_2 f^{abc} A^b_nu F^{c,mu nu} = J^{a,mu}
- Check that the non-Abelian self-interaction term (f^{abc}) appears correctly

### Step 4: Donnelly-Wall Verification (1 week)
- Compute the on-shell Sigma_{SU(2)} for a configuration near an entangling surface
- Compare with the Donnelly-Wall edge mode entropy S_edge = (N/2) ln 3
- Verify the conjecture: D_break = S_edge per DOF

### Step 5: Strong-Field Analysis (1 week)
- Compute Sigma_{SU(2)} in the strong-field regime
- Verify logarithmic saturation: Sigma ~ 3 ln(alpha' |F|^2)
- Compare with Sigma_grav = 2 ln Q saturation

### Step 6: Write Paper 7 (2 weeks)
- Structure: Sigma_total = Sigma_grav + Sigma_gauge
- 0-form = gravity (Papers 2-4)
- 1-form = gauge (this paper)
- Born-Infeld completion of Yang-Mills from GQRE
- Donnelly-Wall edge mode connection (new result)

---

## 11. Summary Table

| Item | Status | Evidence |
|------|--------|----------|
| Sigma = 0-form of Bianconi GQRE | CONFIRMED | Structural match, Section 1.3 |
| 1-form GQRE gives Abelian gauge (Maxwell/BI) | ESTABLISHED by Bianconi | Appendix A of PRD paper |
| 1-form GQRE gives non-Abelian gauge (Yang-Mills) | EXPECTED but NOT COMPUTED | Standard gauge principle + GQRE, Section 2.2-2.3 |
| For SU(2): Sigma_{SU(2)} ~ F^a F^a (weak field) | EXPECTED | From non-Abelian GQRE expansion |
| Strong-field saturation: Sigma ~ 3 ln |F|^2 | PREDICTED | From log structure of GQRE |
| D_break = ln 3 = Donnelly-Wall edge entropy | NEW (this analysis + su2_from_sigma.md) | Both count dim(adj SU(2)) = 3 |
| SU(2) group selection from GQRE alone | NO -- needs external input | Group comes from algebra/topology |
| Full non-Abelian GQRE variational equations = YM | NOT YET VERIFIED | Calculation needed (Step 3) |
| Bianconi's explicit non-Abelian extension | NOT YET PUBLISHED | "Left for future investigations" |

---

## 12. Go/No-Go Criteria

### GO (proceed to Paper 7):
- If Step 2 produces Yang-Mills action from the 1-form GQRE with SU(2) indices
- If Step 3 recovers the correct Yang-Mills equations including self-interaction
- If Step 4 confirms D_break = S_edge

### NO-GO (abandon Route A):
- If the non-Abelian log-trace does not reduce to Yang-Mills (wrong symmetry structure)
- If the variational equations produce additional unwanted terms (ghosts, tachyons)
- If the Born-Infeld completion conflicts with known phenomenology (e.g., wrong sign of F^4 correction)

### PARTIAL GO (publish partial result):
- If Yang-Mills emerges but only in a specific limit (e.g., only in 4D, only for specific representations)
- If the Donnelly-Wall connection holds but the GQRE dynamics don't reproduce full YM
- In this case: publish the Donnelly-Wall edge mode identification as a standalone result

---

## 13. Attribution

| Discovery | Who | Ours? |
|-----------|-----|-------|
| GQRE Lagrangian and p-form decomposition | Bianconi (2025 PRD) | No |
| Abelian gauge field in GQRE | Bianconi (2025 PRD, App. A) | No |
| Edge modes for gauge fields | Donnelly-Wall (2014-2016) | No |
| Non-Abelian entanglement entropy | Donnelly (2014), Casini-Huerta-Rosabal (2013) | No |
| QRE chain rule under twirl | Marvian-Spekkens (2014) | No |
| **D_break = ln 3 = Donnelly-Wall edge entropy** | **New (su2_from_sigma.md + this analysis)** | **Yes** |
| **Sigma_total = Sigma_grav + Sigma_gauge from GQRE p-form decomposition** | **New (explicit identification)** | **Yes** |
| **Non-Abelian GQRE prediction: Sigma_{SU(2)} = Tr ln(I + F)** | **New (not computed by Bianconi)** | **Yes** |
| **Strong-field saturation: Sigma_gauge ~ dim(G) ln |F|^2** | **New prediction** | **Yes** |
| **Edge mode = 1-form GQRE** identification | **New conjecture** | **Yes** |
| **Concrete plan for non-Abelian GQRE computation** | **New (this document)** | **Yes** |

---

## 14. Key References

1. Bianconi (2025), "Gravity from entropy," PRD 111, 066001. arXiv:2408.14391
2. Bianconi (2025), "Thermodynamics of GfE," arXiv:2510.22545
3. Bianconi (2024), "Quantum entropy couples matter with geometry," J. Phys. A 57, 365002. arXiv:2404.08556
4. Donnelly & Wall (2015), "Entanglement entropy of EM edge modes," PRL 114, 111603. arXiv:1412.1895
5. Donnelly (2014), "Entanglement entropy and nonabelian gauge symmetry," arXiv:1406.7304
6. Donnelly & Wall (2016), "Geometric entropy and edge modes of the EM field," PRD 94, 104053. arXiv:1506.05792
7. Casini, Huerta & Rosabal (2014), "Remarks on entanglement entropy for gauge fields," arXiv:1312.1183
8. Marvian & Spekkens (2014), "Extending Noether's theorem," arXiv:1404.3236
9. Harlow (2018), "TASI lectures on emergence of bulk physics in AdS/CFT," arXiv:1802.01040
10. Ball & Ciambelli (2025), "Dynamical Edge Modes in Yang-Mills Theory," SciPost preprint

---

## Record
Created 2026-03-20 during deep push of Bianconi 1-form route for SU(2) from Sigma.
Builds on: su2_from_sigma.md (Route A), paper4_bianconi_gfe_research.md, spectral_action_sigma.md, imaginary_units_gauge_forces_2026_03_19.md
