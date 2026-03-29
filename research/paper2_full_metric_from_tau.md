# Can the Full Spacetime Metric g_ab Be Derived from the tau Framework?

## Date: 2026-03-10
## Status: Rigorous mathematical analysis (supersedes paper2_full_metric_derivation.md)
## Purpose: Systematically evaluate five approaches to deriving g_rr (and the full metric) from information-theoretic principles

---

## 0. Executive Summary

**The question**: Paper 2 derives g_00 = -exp(-r_s/r) from three independent information-theoretic routes, but the spatial component g_rr = exp(+r_s/r) is obtained via the isometry assumption g_00 * g_rr = -1 (Dicke's isotropic refractive index). Can we do better?

**The answer**: After rigorous analysis of five approaches, the honest assessment is:

| Approach | Can derive g_rr? | Gives exponential metric? | Mathematical status |
|----------|-----------------|--------------------------|-------------------|
| 1. Spatial channels | In principle | Possible but unproven | OPEN -- requires new formalism |
| 2. Modified Jacobson/Dorau-Much | YES (via field eqs) | NO (gives Schwarzschild) | PROVEN (for Einstein eqs) |
| 3. Entanglement structure (RT/MERA) | YES (in AdS) | Not tested | PROVEN in AdS; OPEN for asymptotically flat |
| 4. Modular geometry (Connes/Witten) | In principle | Unknown | OPEN -- deepest but hardest |
| 5. Consistency conditions (EFE) | YES | Only with phantom scalar | PROVEN but requires matter input |

**Key finding**: There IS a promising new pathway (Approach 1, via Bisognano-Wichmann spatial modular flow) that was not properly explored in the previous analysis. The Bisognano-Wichmann theorem generates boosts that mix t and x. The modular flow in the (t,r) plane constrains BOTH g_00 and g_rr simultaneously, not just g_00. This is the most promising route to a genuine derivation.

**However**: No approach currently provides a complete, rigorous derivation of g_rr = exp(+r_s/r) from purely information-theoretic principles. The gap is real.

---

## 1. The Problem: Precise Mathematical Statement

### 1.1 What We Have

Paper 2 establishes (via three independent routes):

```
Sigma_grav = -ln(-g_00)
```

Combined with the identification Sigma_grav = r_s/r (from modular flow, gravitational Landauer, and quantum channel arguments):

```
g_00 = -exp(-r_s/r)
```

The full exponential metric in isotropic coordinates is:

```
ds^2 = -exp(-r_s/r) c^2 dt^2 + exp(+r_s/r)(dr^2 + r^2 dOmega^2)
```

### 1.2 What We Need

A derivation of g_rr = exp(+r_s/r) from information-theoretic principles. Equivalently, any of:

(a) A "spatial entropy production" Sigma_spatial such that g_rr = exp(Sigma_spatial) with Sigma_spatial = r_s/r

(b) A derivation of the condition g_00 * g_rr = -1 from information theory

(c) A derivation of the full Einstein equations + the specific matter content (phantom scalar phi = M/r) that sources the exponential metric

(d) A reconstruction of the full spatial metric from entanglement entropy data

### 1.3 Why This Is Hard (The Fundamental Asymmetry)

Quantum mechanics treats time and space asymmetrically:

- **Time**: Evolution is described by CPTP maps (channels). The Petz recovery map recovers past from future. Entropy production Sigma = QRE drop under a channel is well-defined. Channels are fundamentally temporal.

- **Space**: There is no "spatial channel" in standard quantum mechanics. Spatial correlations are described by the state, not by a dynamical map. The distinction between "input" and "output" (essential for defining a channel) is temporal, not spatial.

This asymmetry is why Sigma_grav = -ln(-g_00) naturally constrains the TEMPORAL metric component but has no direct spatial analogue.

**However**: In relativistic QFT, the distinction between time and space is less sharp. The modular flow (Tomita-Takesaki theory) generates evolution that can mix temporal and spatial directions. This is the key insight that opens new possibilities.

---

## 2. Approach 1: Spatial Quantum Channels via Modular Flow

### 2.1 The Bisognano-Wichmann Theorem (The Key Tool)

**Theorem** (Bisognano-Wichmann 1975-76): For a free quantum field in Minkowski space, the modular automorphism group of the vacuum state restricted to the right Rindler wedge W_R = {x : x^1 > |x^0|} is:

```
sigma_s^{BW}(A) = U(Lambda_{2pi s}) A U(Lambda_{2pi s})^{-1}
```

where Lambda_{2pi s} is the Lorentz boost with rapidity 2pi s:

```
Lambda_{2pi s}: (t, x) -> (t cosh(2pi s) + x sinh(2pi s), x cosh(2pi s) + t sinh(2pi s))
```

**Critical observation**: The boost Lambda_{2pi s} MIXES t and x. It is NOT a pure time translation. The modular flow simultaneously acts on both temporal and spatial directions.

### 2.2 The Modular Hamiltonian and Both Metric Components

The modular Hamiltonian for the Rindler wedge is:

```
K_R = 2pi integral d^{d-1}x  x^1 T_{00}(0, x)
```

In curved spacetime (for a bifurcate Killing horizon with surface gravity kappa), the modular Hamiltonian generalizes (Sewell 1982, Jacobson-Wall):

```
K_xi = (2pi / kappa) integral_Sigma xi^a T_{ab} n^b sqrt(h) d^{d-1}x
```

where xi^a is the Killing vector that generates the horizon, n^b is the unit normal to the Cauchy surface Sigma, and h is the induced metric determinant.

**For the Schwarzschild/exponential metric**, the Killing vector is xi = partial_t, and on a t = const slice:

```
K = (2pi / kappa) integral sqrt(-g_00) T_{00} sqrt(g_rr) r^2 sin(theta) dr d(theta) d(phi)
```

The factor sqrt(-g_00) * sqrt(g_rr) appears naturally in K. This means the modular Hamiltonian depends on BOTH g_00 AND g_rr.

### 2.3 The Spatial Modular Flow Argument

**Proposition** (New): In Rindler space, the boost Killing vector xi = x partial_t + t partial_x generates flow in the (t,x) plane. The QRE associated with this boost can be decomposed:

```
S^rel = (2pi / kappa) integral_H (xi . T) dH
```

where the integral runs over the bifurcation surface. By the Bisognano-Wichmann theorem, this QRE encodes information about the full (t,x) sector of the metric, not just the temporal part.

**The critical calculation**: Consider two nested Rindler wedges, W_R(a) subset W_R, where W_R(a) = {x : x^1 > a + |x^0|}. The QRE between the vacuum restricted to W_R and W_R(a) involves the modular Hamiltonian of W_R:

```
S^rel(omega|_{W_R(a)} || omega|_{W_R}) = <omega| K_R - K_{R(a)} |omega>
```

where K_R - K_{R(a)} depends on the metric in the strip a < x < a + delta. By evaluating this in the curved spacetime case, one obtains constraints on both g_00 and g_rr.

### 2.4 Specific Equations for the Schwarzschild Geometry

For a static spherically symmetric metric ds^2 = -A(r) dt^2 + B(r) dr^2 + r^2 dOmega^2, consider a family of "Rindler-like" horizons at different radii. The local Rindler approximation at radius r_0 gives:

```
ds^2 approx -kappa^2 rho^2 dt^2 + d(rho)^2 + r_0^2 dOmega^2
```

where rho is the proper distance from the local horizon and kappa = A'(r_0) / (2 sqrt(A(r_0) B(r_0))) is the local surface gravity.

The modular Hamiltonian for this local Rindler approximation is:

```
K_local = 2pi integral rho T_{00} rho d(rho) dOmega
```

The modular flow at imaginary time s = i/(2pi) gives the KMS condition:

```
<A sigma_{i/(2pi)}(B)> = <BA>
```

at the Tolman temperature T(r) = kappa / (2pi sqrt(A(r))).

**The constraint on g_rr**: The proper distance rho is related to the coordinate distance by:

```
rho = integral_{r_H}^{r} sqrt(B(r')) dr'
```

The modular Hamiltonian involves rho * T_{00}, which depends on B(r) through the proper distance. Therefore, the modular flow constraints on the QRE involve both A(r) and B(r).

### 2.5 The Precise Constraint

**Claim**: The modular Hamiltonian of the local Rindler approximation at radius r constrains the combination:

```
kappa(r) = A'(r) / (2 sqrt(A(r) B(r)))
```

This is a SINGLE equation relating A(r) and B(r). Combined with the constraint A(r) = exp(-r_s/r) from the tau framework (Sigma_grav = -ln(A)), this gives:

```
A'(r) = (r_s/r^2) exp(-r_s/r)
kappa(r) = (r_s/r^2) exp(-r_s/r) / (2 sqrt(exp(-r_s/r) B(r)))
         = (r_s/(2r^2)) exp(-r_s/(2r)) / sqrt(B(r))
```

For the exponential metric: B(r) = exp(+r_s/r), so:

```
kappa(r) = (r_s/(2r^2)) exp(-r_s/(2r)) / exp(+r_s/(2r))
         = (r_s/(2r^2)) exp(-r_s/r)
```

For Schwarzschild (in Schwarzschild coordinates, A = 1 - r_s/r, B = 1/(1-r_s/r)):

```
kappa(r) = (r_s/r^2) (1/(2 * 1)) = r_s/(2r^2)  [at r >> r_s]
```

The surface gravity kappa involves A'/(2 sqrt(AB)). If modular flow can independently determine kappa(r) from information-theoretic principles, then knowing A(r) determines B(r).

### 2.6 Assessment of Approach 1

| Aspect | Status |
|--------|--------|
| Mathematical formulation | PARTIALLY DEFINED -- the modular flow argument is well-posed |
| Key calculation needed | Evaluate S^rel between nested Rindler wedges in curved spacetime |
| Physical motivation | STRONG -- Bisognano-Wichmann naturally mixes t and x |
| Can distinguish exp from Schwarzschild? | IN PRINCIPLE -- different B(r) gives different kappa(r) |
| Existing literature support | Sewell 1982, Jacobson 1995/2016, Dorau-Much 2025, Wall 2012 |
| Status | PROMISING but requires explicit calculation |

**The critical calculation that would resolve this**:

Compute the QRE between the vacuum state restricted to the exterior of radius r and the exterior of radius r + dr, for a general static spherically symmetric metric with A(r) and B(r). If this QRE depends on both A and B independently (not just through their ratio or product), then modular flow constrains both metric components.

**PROVEN**: The modular Hamiltonian depends on sqrt(A * B) through the proper distance factor.

**CONJECTURED**: This dependence, combined with Sigma_grav = -ln(A), uniquely determines B(r).

**OPEN**: Whether the QRE constraint gives B = 1/A (exponential metric) or some other relation.

---

## 3. Approach 2: Modified Jacobson/Dorau-Much

### 3.1 The Standard Jacobson Argument

Jacobson (1995) derives the full Einstein equation from the Clausius relation delta-Q = T delta-S applied to ALL local Rindler horizons. The key steps:

**Step 1**: For a pencil of null geodesic generators xi^a of a local Rindler horizon, the Raychaudhuri equation gives:

```
d(theta)/d(lambda) = -(1/2) theta^2 - sigma_{ab} sigma^{ab} + omega_{ab} omega^{ab} - R_{ab} xi^a xi^b
```

For a congruence that is initially non-expanding (theta_0 = 0), to first order:

```
delta-A = -integral R_{ab} xi^a xi^b d(lambda) dA
```

**Step 2**: The Bekenstein-Hawking entropy S = A/(4G) converts this to:

```
delta-S = -(1/(4G)) integral R_{ab} xi^a xi^b d(lambda) dA
```

**Step 3**: The Clausius relation delta-Q = T delta-S, with T = kappa/(2pi) (Unruh temperature) and delta-Q = integral T_{ab} xi^a xi^b d(lambda) dA, gives:

```
T_{ab} xi^a xi^b = (1/(8pi G)) R_{ab} xi^a xi^b
```

**Step 4 (The tensorialization)**: Since this holds for ALL null vectors xi^a (all boost directions at every point), the scalar equation upgrades to:

```
R_{ab} - (R/2) g_{ab} + Lambda g_{ab} = 8pi G T_{ab}
```

The cosmological constant Lambda appears as an integration constant.

### 3.2 Why This Gives Schwarzschild, Not Exponential

The Einstein equations in vacuum (T_{ab} = 0) with Lambda = 0 give:

```
R_{ab} = 0
```

For a static spherically symmetric metric ds^2 = -A(r) dt^2 + B(r) dr^2 + r^2 dOmega^2:

The (t,t) component: A''/(2B) - A'B'/(4B^2) + A'^2/(4AB) + A'/(rB) = 0
The (r,r) component: -A''/(2A) + A'B'/(4AB) + A'^2/(4A^2) + B'/(rB) = 0

Adding these: A'/(rB) + B'/(rB) = 0, which gives (AB)' = 0, hence AB = const. Boundary conditions fix const = 1: AB = 1.

Combined with the remaining equation, this gives A = 1 - r_s/r (Schwarzschild).

**The exponential metric does NOT satisfy R_{ab} = 0**. It requires matter:

```
R_{ab} - (R/2) g_{ab} = 8pi G T_{ab}^{phantom}
```

where T_{ab}^{phantom} is the stress-energy of a phantom (wrong-sign kinetic term) scalar field phi = M/r.

### 3.3 Can We Modify the Argument to Get the Exponential Metric?

**Idea**: Replace the Bekenstein-Hawking area law delta-S = delta-A/(4G) with a modified entropy-area relation motivated by the tau framework.

**Specifically**: If the Petz bound F >= exp(-Sigma/2) is the fundamental relation (rather than S = A/(4G)), we need to convert this into an entropy-area relation.

The Petz bound gives:

```
Sigma = -ln(-g_00) = r_s/r
```

which is related to the redshift factor, NOT the area. The area law S = A/(4G) is a statement about the ENTROPY of the horizon, while Sigma is a statement about the FIDELITY of information recovery.

**The mathematical question**: Is there a modified Clausius relation that uses Sigma instead of S?

**Attempt**: Replace delta-S with delta-Sigma:

```
delta-Q = T delta-Sigma
```

Then:

```
T_{ab} xi^a xi^b d(lambda) dA = T * delta(-ln(-g_00)) * dA
```

But -ln(-g_00) is a GLOBAL quantity (defined at a specific radius), not a LOCAL quantity like the area element. The Jacobson argument requires a LOCAL entropy that can be varied under deformations of the horizon. The redshift factor -ln(-g_00) does not have this property -- it is defined far from the horizon.

**Assessment**: The modification delta-S -> delta-Sigma is NOT mathematically well-posed because Sigma is not a local horizon quantity.

### 3.4 The Dorau-Much Variation

Dorau and Much (2025) replace the Clausius relation with the QRE directly:

```
S^rel(omega_0 || omega_phi) = -(2pi) integral_H U <:T_{UU}:> dU dA
```

This is more rigorous because it uses the well-defined Araki-Uhlmann relative entropy instead of the somewhat heuristic Clausius relation.

**Can we modify this to get the exponential metric?**

The Dorau-Much derivation assumes:
1. The Bekenstein-Hawking area law: delta-A/(4G) = entropy change
2. The Raychaudhuri equation: relates delta-A to R_{ab}

Both of these are standard GR results. To get the exponential metric instead of Schwarzschild, we would need to modify either the area law or the Raychaudhuri equation.

**Option A: Modify the area law**

Replace S = A/(4G) with S = f(A) for some function f. This changes the field equations to:

```
f'(A) R_{ab} xi^a xi^b = alpha T_{ab} xi^a xi^b
```

For f'(A) = const (linear), this gives the standard Einstein equations. For nonlinear f, this gives modified gravity (like f(R) gravity). The exponential metric could potentially arise from a specific f(A).

**However**: The choice of f(A) would be ad hoc and not motivated by information theory. We would be replacing one assumption (isometry condition) with another (modified area law).

**Option B: Modify the Raychaudhuri equation**

The Raychaudhuri equation is a KINEMATIC identity -- it does not depend on the field equations. It relates the expansion rate to the Ricci tensor:

```
d(theta)/d(lambda) = -(1/(d-2)) theta^2 - sigma^2 + omega^2 - R_{ab} xi^a xi^b
```

This is exact and cannot be modified without changing the geometric framework.

### 3.5 Assessment of Approach 2

| Aspect | Status |
|--------|--------|
| Can derive full Einstein equations? | YES -- PROVEN (Jacobson 1995, Dorau-Much 2025) |
| Gives exponential metric? | NO -- gives Schwarzschild (in vacuum) |
| Can be modified to give exponential? | NOT without ad hoc modifications |
| Information-theoretic basis? | YES (strongest available) |
| Mathematical rigor? | HIGH |
| Verdict | The approach works perfectly but gives the WRONG metric for Paper 2 |

**PROVEN**: Jacobson/Dorau-Much => Einstein equations => Schwarzschild.

**PROVEN**: The exponential metric does NOT satisfy the vacuum Einstein equations.

**CONCLUSION**: Approach 2 is a dead end for deriving the exponential metric. It provides strong evidence that the STANDARD information-theoretic route leads to Schwarzschild, and the exponential metric requires ADDITIONAL input beyond what Jacobson/Dorau-Much provide.

This creates a genuine tension: the most rigorous information-theoretic approach to gravity (Jacobson) gives Schwarzschild, while Paper 2's direct identification Sigma = -ln(-g_00) gives the exponential metric. They cannot both be right in strong field. This is a TESTABLE prediction.

---

## 4. Approach 3: Entanglement Structure (RT/MERA/Tensor Networks)

### 4.1 The Ryu-Takayanagi Formula

For a static spacetime in AdS/CFT, the entanglement entropy of a boundary region A is:

```
S_A = min_{gamma_A} Area(gamma_A) / (4G_N)
```

where gamma_A is the minimal surface in the bulk anchored to the boundary of A.

**How this constrains spatial geometry**: The area of gamma_A depends on the FULL bulk metric:

```
Area(gamma_A) = integral_{gamma_A} sqrt(det(h_{ij})) d^{d-2}sigma
```

where h_{ij} is the induced metric on gamma_A. For a bulk metric with both g_00 and g_rr components, the minimal surface probe different combinations depending on its orientation.

### 4.2 Metric Reconstruction from Entanglement

**Czech et al. (2015)** introduced "kinematic space" -- a Lorentzian geometry whose metric is defined by conditional mutual information:

```
ds^2_{kinematic} ~ partial^2 S_A / (partial u partial v) du dv
```

where u, v parametrize the endpoints of boundary intervals.

The bulk metric can be reconstructed from entanglement entropy data via an integral transform (analogous to the Radon transform):

```
Length(gamma) = integral_{kinematic space} (contribution from entanglement data)
```

**Key result**: For AdS_3/CFT_2, the spatial metric of a constant-time slice of AdS_3 can be FULLY reconstructed from the boundary entanglement entropy data. The reconstruction determines BOTH the radial metric component and the angular metric component.

### 4.3 Cao-Carroll-Michalakis (2017): Space from Hilbert Space

This is the most explicit demonstration that entanglement can determine spatial metric components.

**Method**:
1. Decompose Hilbert space into a tensor product H = tensor_{alpha} H_alpha
2. Define a "distance" between factors using mutual information: d(alpha, beta) ~ 1/I(alpha : beta)
3. Use multidimensional scaling to extract the best-fit spatial geometry
4. The induced spatial metric g_{ij}^{spatial} is determined by the entanglement pattern

**Result**: For "redundancy-constrained states" (states with area-law entanglement), the spatial geometry is recoverable, and perturbations of entanglement entropy give rise to curvature obeying the linearized Einstein equation.

### 4.4 Application to the tau Framework

**The question**: Can we use the Cao-Carroll-Michalakis framework to determine g_rr in the exponential metric?

**Requirements**:
1. Need a boundary theory (CFT) dual to the exponential metric -- DOES NOT EXIST (the exponential metric is asymptotically flat, not asymptotically AdS)
2. Need entanglement entropy data for the boundary theory -- NOT AVAILABLE
3. Need the RT formula to apply -- REQUIRES AdS/CFT

**The fundamental obstacle**: The RT formula and holographic reconstruction are proven only in the AdS/CFT context. The exponential metric is asymptotically flat. There is no known holographic dual.

**Possible resolution**: Flat-space holography (celestial holography, Carrollian CFT) is under active development. If a holographic duality for asymptotically flat spacetimes is established, the RT-like formula could be applied to reconstruct the spatial metric.

### 4.5 Can Entanglement Structure Distinguish Exponential from Schwarzschild?

Even without a full holographic duality, we can ask: does the entanglement structure of the vacuum differ between Schwarzschild and exponential geometries?

**For Schwarzschild**: The vacuum entanglement across the horizon is the Hawking entropy S = A/(4G). The entanglement at radius r scales as:

```
S(r) ~ A(r)/(4G) = 4pi r^2 / (4G)    [area law]
```

with corrections from the curvature.

**For the exponential metric**: There is no horizon, so the entanglement structure is fundamentally different. The entanglement entropy at radius r would be:

```
S(r) ~ 4pi r^2 exp(r_s/r) / (4G)    [modified by g_rr]
```

The factor exp(r_s/r) comes from the spatial metric affecting the proper area of the surface at coordinate radius r.

**The key point**: The entanglement entropy S(r) depends on g_rr through the proper area sqrt(det g_{angular}). In principle, measuring S(r) as a function of r determines g_rr(r).

### 4.6 Specific Equations

For a static spherically symmetric metric, the entanglement entropy across a sphere at radius r is (to leading order in the UV cutoff epsilon):

```
S(r) = c_0 * Area_proper(r) / epsilon^{d-2} + c_1 * integral_Sigma R_{induced} sqrt{h} d^{d-2}x + ...
```

where Area_proper(r) = 4pi r^2 sqrt{g_{theta theta}(r) g_{phi phi}(r)} depends on the angular metric components.

For the exponential metric in isotropic coordinates:

```
g_{theta theta} = exp(r_s/r) * r^2
```

so Area_proper(r) = 4pi r^2 exp(r_s/r).

For Schwarzschild in isotropic coordinates:

```
g_{theta theta} = (1 + r_s/(4r))^4 * r^2
```

so Area_proper(r) = 4pi r^2 (1 + r_s/(4r))^4.

In principle, measuring S(r) determines the angular metric, and hence g_rr (by spherical symmetry, g_rr is related to g_{theta theta} via the coordinate choice).

### 4.7 Assessment of Approach 3

| Aspect | Status |
|--------|--------|
| Mathematical framework | PROVEN in AdS/CFT; OPEN for flat space |
| Can determine g_rr in principle? | YES -- entanglement area law constrains spatial metric |
| Gives exponential metric? | NOT TESTED |
| Requires AdS/CFT? | Currently YES -- flat-space extension is open |
| Connection to tau framework? | INDIRECT -- entanglement entropy != Petz fidelity |
| Mathematical rigor | HIGH (in AdS); SPECULATIVE (for flat space) |

**PROVEN**: In AdS/CFT, spatial metric reconstruction from entanglement is possible (Czech 2015, Cao-Carroll 2017).

**CONJECTURED**: A similar reconstruction should work for asymptotically flat spacetimes, once flat-space holography is developed.

**OPEN**: Whether the tau framework's Sigma_grav is related to the entanglement entropy used in holographic reconstruction.

**Verdict**: This is a PROMISING but INCOMPLETE pathway. The key missing ingredient is a flat-space analogue of the RT formula. If one is found, it could determine g_rr from entanglement data.

---

## 5. Approach 4: Modular Geometry (Connes / Witten / Chandrasekaran)

### 5.1 Connes' Spectral Distance

In noncommutative geometry, Connes defines a distance on a (possibly noncommutative) space using the Dirac operator D:

```
d(p, q) = sup { |f(p) - f(q)| : ||[D, pi(f)]|| <= 1 }
```

For a Riemannian manifold with the standard Dirac operator, this recovers the geodesic distance. The entire metric is encoded in the spectral data (A, H, D) -- the algebra, Hilbert space, and Dirac operator.

**Application to our problem**: If we can identify the "gravitational Dirac operator" D_grav from information-theoretic principles, the full metric is determined by Connes' formula.

**The modular operator as Dirac operator**: In Tomita-Takesaki theory, the modular operator Delta satisfies:

```
Delta^{it} A Delta^{-it} = sigma_t(A)    [modular automorphism]
```

The "modular Dirac operator" would be ln(Delta). For the vacuum state of a QFT, the Bisognano-Wichmann theorem identifies:

```
ln(Delta) = -2pi K    [K = boost generator = modular Hamiltonian]
```

In curved spacetime, K depends on the full metric. Therefore, if we can independently determine K from information-theoretic principles, the metric is determined.

### 5.2 Witten's Crossed Product (2022)

Witten showed that in perturbative quantum gravity around a black hole background:

1. The algebra of observables on one side of the horizon is Type III_1 (no trace, no entropy)
2. Including 1/N corrections (gravitational dressing to an observer), the algebra becomes Type II_infinity (has a trace)
3. This Type II_infinity algebra is the crossed product of the Type III_1 algebra by its modular automorphism group

The entropy in the Type II_infinity algebra is:

```
S = -Tr(rho ln rho)    [well-defined for Type II]
```

and agrees with the generalized entropy S_gen = A/(4G) + S_outside.

### 5.3 Chandrasekaran-Longo-Penington-Witten (2023)

For de Sitter space, they construct a Type II_1 algebra (finite entropy) by dressing operators to the observer's worldline. The modular flow of this algebra generates:

```
sigma_t^{obs}(A) = e^{iHt} A e^{-iHt}
```

where H is the static Hamiltonian in the observer's causal diamond.

**The full metric is encoded in the modular structure**: The modular operator Delta depends on the state and the algebra. For the Hartle-Hawking state on the Schwarzschild background, the modular flow is the Killing time translation. The ENTIRE metric structure (g_00, g_rr, angular parts) determines the modular operator through:

```
Delta = exp(-beta H_static)
```

where beta = 1/T_Hawking and H_static = integral sqrt(-g_00) sqrt(g_rr) r^2 T_{00} dr dOmega.

### 5.4 Can the Modular Hamiltonian Determine Both g_00 AND g_rr?

**The key equation**: For a KMS state at temperature T = 1/beta on a static spacetime, the modular Hamiltonian is:

```
K = beta H_static = beta integral_Sigma xi^a T_{ab} n^b sqrt{|h|} d^{d-1}x
```

where:
- xi = partial_t is the Killing vector with norm |xi| = sqrt(-g_00)
- n^b = xi^b / |xi| is the unit normal to the t = const slice
- sqrt{|h|} = sqrt{g_rr * g_{theta theta} * g_{phi phi}} is the determinant of the induced spatial metric

Expanding:

```
K = beta integral sqrt(-g_00) * T_{00} * sqrt(g_rr * g_{theta theta} * g_{phi phi}) dr d(theta) d(phi)
```

For a spherically symmetric metric with g_{theta theta} = f(r) * r^2:

```
K = 4pi beta integral sqrt(-g_00(r)) * sqrt(g_rr(r)) * sqrt(f(r)) * r^2 * T_{00}(r) dr
```

**The modular Hamiltonian depends on the combination sqrt(g_00 * g_rr * f)**. This is ONE function of three unknowns (g_00, g_rr, f). Combined with the tau-framework constraint g_00 = -exp(-r_s/r), this gives:

```
K = 4pi beta integral exp(-r_s/(2r)) * sqrt(g_rr(r) * f(r)) * r^2 * T_{00}(r) dr
```

So K constrains the combination sqrt(g_rr * f), but NOT g_rr and f independently.

**For the exponential metric in isotropic coordinates**: g_rr = f(r) (both equal exp(r_s/r)), so g_rr * f = exp(2r_s/r), and sqrt(g_rr * f) = exp(r_s/r).

**For Schwarzschild in isotropic coordinates**: g_rr = (1 + r_s/(4r))^4 = f(r), so g_rr * f = (1 + r_s/(4r))^8.

These ARE different, so in principle the modular Hamiltonian distinguishes them.

### 5.5 The KMS Condition and Both Metric Components

The KMS condition states that for the modular flow sigma_t at inverse temperature beta:

```
<A sigma_{ibeta}(B)> = <BA>
```

In the Rindler-like approximation near radius r, the KMS condition is satisfied at the Tolman temperature:

```
T(r) = T_Hawking / sqrt(-g_00(r))
```

This depends only on g_00. HOWEVER, the modular flow sigma_t is generated by the modular Hamiltonian K, which depends on g_rr (as shown above). So:

- The TEMPERATURE (periodicity of modular flow) depends only on g_00
- The GENERATOR of the flow depends on both g_00 and g_rr
- The FLOW ITSELF (orbit structure in the (t,r) plane) depends on both

**The flow orbits**: For the Rindler-like approximation, the flow orbits are hyperbolas in the (t, rho) plane where rho = integral sqrt(g_rr) dr is the proper radial distance. The shape of these hyperbolas depends on g_rr through rho(r).

### 5.6 The Spectral Geometry Argument

**Claim**: If we know the FULL spectrum of the modular Hamiltonian K (not just its temperature), we can reconstruct both g_00 and g_rr.

**Reasoning**: The spectrum of K is determined by:

```
Spec(K) = { beta * E_n : n = 0, 1, 2, ... }
```

where E_n are the eigenvalues of the static Hamiltonian H_static. These eigenvalues depend on the potential well created by the curved geometry, which involves both A(r) and B(r).

For a scalar field of mass m in the background metric:

```
(-g)^{-1/2} partial_mu (sqrt{-g} g^{mu nu} partial_nu phi) = m^2 phi
```

The radial equation in the static case:

```
-(1/(r^2 sqrt(AB))) d/dr (r^2 sqrt(A/B) d(phi_l)/dr) + (l(l+1)/(r^2)) phi_l = (omega^2/A) phi_l
```

The eigenfrequencies omega_n depend on BOTH A(r) and B(r). Therefore, the spectrum of the modular Hamiltonian encodes both metric components.

**However**: Recovering A(r) and B(r) from the spectrum is an INVERSE problem (like hearing the shape of a drum). Such inverse problems are generically ill-posed (non-unique). Explicit conditions under which the metric is uniquely determined by the spectrum are not known.

### 5.7 Assessment of Approach 4

| Aspect | Status |
|--------|--------|
| Mathematical depth | DEEPEST of all approaches |
| Framework exists? | YES -- Tomita-Takesaki, crossed products, Connes' NCG |
| Can determine g_rr in principle? | YES -- modular Hamiltonian depends on g_rr |
| Gives specific g_rr? | UNKNOWN -- inverse spectral problem |
| Practical computability | VERY LOW -- requires full Type II algebra computation |
| Connection to tau framework | INDIRECT -- modular flow generates the tau-relevant dynamics |
| Existing results | Witten 2022, Chandrasekaran 2023 -- foundational but incomplete |

**PROVEN**: The modular Hamiltonian depends on both g_00 and g_rr.

**PROVEN**: The spectrum of K encodes information about both metric components.

**CONJECTURED**: The full modular structure (algebra + state + modular flow) uniquely determines the metric.

**OPEN**: Whether this determination is unique (inverse spectral problem) and whether it can be carried out explicitly.

**Verdict**: This is the DEEPEST approach and the one most likely to succeed in principle. But it is also the hardest to make concrete. The key obstacle is the inverse spectral problem. A breakthrough would be a theorem showing that the modular structure of a specific class of algebras (e.g., those arising from free fields on static spacetimes) uniquely determines the background metric.

---

## 6. Approach 5: Einstein Equations as Consistency Conditions

### 6.1 The Logic

If g_00 = -exp(-r_s/r) is fixed by the tau framework, the Einstein equations G_{mu nu} = 8pi G T_{mu nu} (plus specific matter content T_{mu nu}) determine g_rr.

**For vacuum** (T_{mu nu} = 0): G_{mu nu} = 0 forces AB = 1 (in Schwarzschild coordinates) AND A = 1 - r_s/r. Since g_00 = -exp(-r_s/r) != -(1 - r_s/r), the exponential metric is NOT a vacuum solution. The vacuum Einstein equations are INCONSISTENT with g_00 = -exp(-r_s/r).

**For matter** (specific T_{mu nu}): The exponential metric requires a specific matter source.

### 6.2 The Phantom Scalar Field

The exponential metric ds^2 = -exp(-r_s/r) dt^2 + exp(+r_s/r) dr^2 + r^2 exp(+r_s/r) dOmega^2 (in isotropic coordinates) satisfies the Einstein equations with a phantom scalar field:

```
phi(r) = sqrt(r_s / (8pi G)) * (1/r)     [in isotropic coordinates]
```

with the wrong-sign kinetic term:

```
L_phi = +(1/2) g^{ab} partial_a phi partial_b phi    [+ instead of -]
```

The stress-energy tensor:

```
T_{ab}^{phantom} = partial_a phi partial_b phi - (1/2) g_{ab} (partial phi)^2
```

with the sign flip in the kinetic term. This gives:

```
T_{00} = (1/2) exp(-r_s/r) (r_s/(8pi G r^2))^2 exp(-r_s/r) [from the gradient]
```

### 6.3 Can the Phantom Scalar Be Derived from the tau Framework?

**The question**: Is there an information-theoretic reason for a phantom scalar field phi = M/r?

**Observation (from Paper 2 research)**: Sigma_grav = r_s/r = 2GM/(c^2 r) has the same 1/r profile as the phantom scalar phi. This suggests a deep connection:

```
phi(r) propto Sigma_grav(r)
```

The entropy production Sigma IS (proportional to) the phantom scalar field. If we interpret Sigma as a physical field (not just an abstract information-theoretic quantity), then:

**Conjecture**: The gravitational entropy production Sigma_grav acts as a phantom scalar field sourcing the metric. The metric is self-consistently determined by:

```
g_00 = -exp(-Sigma_grav)    [from tau framework]
Sigma_grav solves the scalar field equation    [from EFE]
g_rr determined by EFE + phantom T_{ab}    [from consistency]
```

**Is this self-consistent?** Let us check:

1. Start with Sigma = r_s/r (from information theory)
2. This gives g_00 = -exp(-r_s/r) (tau framework)
3. The phantom scalar phi propto Sigma = r_s/r
4. Einstein equations with phantom source give g_rr = exp(+r_s/r)
5. The scalar field equation for phi = M/r in the background metric:

```
Box_g phi = (1/sqrt{-g}) partial_mu (sqrt{-g} g^{mu nu} partial_nu phi) = 0
```

For the exponential metric, this equation IS satisfied (Makukov-Mychelkin 2020). The exponential metric is an exact solution of the Einstein-phantom-scalar system.

**The circularity problem**: This argument is self-consistent but CIRCULAR:
- Step 1 assumes Sigma = r_s/r
- Step 3 identifies Sigma with the phantom scalar
- Step 5 verifies consistency

The circularity is: we assumed the metric (through Sigma = r_s/r) and then verified that it is a solution. We did NOT derive the metric from scratch.

### 6.4 Breaking the Circularity

To make Approach 5 non-circular, we would need to INDEPENDENTLY derive:

(a) Sigma_grav = -ln(-g_00) [from tau framework -- DONE]
(b) Sigma_grav = phi (phantom scalar identification) [from physics -- NEW CLAIM]
(c) The scalar field equation Box_g phi = 0 [from QFT -- standard]
(d) The Einstein equations G_{mu nu} = 8pi G T_{mu nu}^{phantom} [from Jacobson/Dorau-Much -- DONE]

If all four are independently derived, the system of equations determines g_00, g_rr, and phi simultaneously, with NO circularity.

**Status of each step**:

(a) ESTABLISHED -- three independent routes give Sigma = -ln(-g_00)

(b) CONJECTURED but NOT DERIVED -- the identification Sigma = phi is motivated by the 1/r profile but lacks a derivation from first principles. The physical meaning of "entropy production IS a scalar field" needs clarification.

(c) STANDARD -- the massless scalar field equation is well-known

(d) ESTABLISHED -- Jacobson (1995) and Dorau-Much (2025) derive Einstein equations from information theory. BUT: these give the STANDARD Einstein equations, which with a phantom scalar give the exponential metric; with NO matter, they give Schwarzschild.

**The key unresolved question for step (b)**: WHY should the entropy production Sigma_grav be identified with a phantom scalar field? What physical principle connects the information-theoretic quantity (QRE drop) to a dynamical scalar field?

### 6.5 A Possible Physical Argument for Sigma = phi

**Bianconi (PRD 2025)**: In the "Gravity from Entropy" framework, the metric g_ab IS treated as a density matrix, and the gravitational action is the QRE between the spacetime metric and the matter-induced metric. In this framework, entropy production naturally has the interpretation of a field.

**Information field theory (IFT)**: Ensslin et al. have developed a framework where physical fields are treated as information-theoretic quantities. In IFT, a scalar field phi(x) is a random variable, and the action is the negative log-likelihood. The entropy production Sigma could be viewed as a field in this sense.

**The Verlinde argument**: Verlinde (2010) argued that gravity is an entropic force, with the gravitational potential phi proportional to the entropy gradient. If phi propto Sigma, then the phantom scalar is literally the entropy field.

**However**: None of these arguments rigorously derive Sigma = phi from first principles. They provide physical motivation but not proof.

### 6.6 Assessment of Approach 5

| Aspect | Status |
|--------|--------|
| Can determine g_rr? | YES -- given specific matter content, EFE determine g_rr |
| Gives exponential metric? | YES -- with phantom scalar phi = M/r |
| Information-theoretic basis for matter? | PARTIAL -- Sigma propto phi is motivated but not proven |
| Self-consistency? | YES -- exponential metric solves Einstein-phantom system |
| Circularity-free? | NO -- currently circular; needs independent derivation of Sigma = phi |
| Mathematical rigor | HIGH (for the EFE part); LOW (for the Sigma = phi identification) |

**PROVEN**: The exponential metric is an exact solution of Einstein + phantom scalar (Makukov-Mychelkin 2020).

**PROVEN**: Sigma_grav = r_s/r has the same profile as the phantom scalar phi = M/r.

**CONJECTURED**: Sigma_grav IS the phantom scalar field (up to normalization).

**OPEN**: Physical derivation of Sigma = phi from information-theoretic principles.

**Verdict**: This approach CAN give g_rr, but only by importing the additional assumption that Sigma is a physical scalar field. This assumption is no better (or worse) than the isometry assumption g_00 * g_rr = -1. Both are physically motivated but underived.

---

## 7. Synthesis: The Landscape of Possibilities

### 7.1 Comparison of All Five Approaches

| Approach | g_rr derivable? | Exponential metric? | Key assumption | Most important gap |
|----------|----------------|--------------------|-----------------|--------------------|
| 1. Spatial modular flow | In principle | Possible | Modular flow constrains metric independently | Explicit curved-space calculation |
| 2. Modified Jacobson | Yes | No (Schwarzschild) | Standard area law | Cannot be modified without ad hoc choices |
| 3. Entanglement (RT/MERA) | Yes (in AdS) | Not tested | RT formula | No flat-space holographic dual |
| 4. Modular geometry | In principle | Unknown | Inverse spectral problem solvable | Extremely difficult computation |
| 5. EFE consistency | Yes | Yes (with phantom) | Sigma = phantom scalar | Circularity in current form |

### 7.2 The Most Promising Path Forward

**Approach 1 (Spatial modular flow)** is the most promising for a genuinely new result, for the following reasons:

1. It uses EXISTING mathematical tools (Bisognano-Wichmann, modular flow) that are well-understood
2. The modular flow in the (t,r) plane naturally constrains both g_00 and g_rr
3. The calculation is well-defined (QRE between states on nested Rindler wedges)
4. It does NOT require AdS/CFT (works in flat or curved spacetime directly)
5. It could potentially distinguish between Schwarzschild and exponential metrics

**The critical computation**: Evaluate the QRE drop for the modular flow in a general static spherically symmetric background with metric functions A(r) and B(r). Determine whether the combination of:

```
Sigma_temporal = -ln(A(r))    [from tau framework, established]
```

and

```
Sigma_total(A, B) = QRE drop from modular flow    [new computation needed]
```

uniquely determines B(r).

**If Sigma_total depends on A and B independently**: Then B is determined, and g_rr follows.

**If Sigma_total depends only on A**: Then g_rr remains undetermined, confirming the fundamental time-space asymmetry.

### 7.3 A New Conjecture

Based on the analysis above, I propose:

**Conjecture (Modular Flow Metric Determination)**: For a static spherically symmetric spacetime with metric ds^2 = -A(r) dt^2 + B(r) dr^2 + C(r) dOmega^2, the modular structure of the vacuum state restricted to the exterior of a sphere of radius r determines the combination:

```
kappa(r) = A'(r) / (2 sqrt(A(r) B(r)))    [local surface gravity]
```

Combined with the tau-framework constraint A(r) = exp(-Sigma_grav(r)), this determines B(r) up to a single function.

**If** the local surface gravity has an independent information-theoretic expression (e.g., through the Unruh temperature of an accelerated observer at r), then B(r) is FULLY determined:

```
B(r) = [A'(r)]^2 / (4 A(r) kappa(r)^2)
```

For the exponential metric: A = exp(-r_s/r), A' = (r_s/r^2) exp(-r_s/r), and:

```
kappa(r) = (r_s/(2r^2)) exp(-r_s/r) / sqrt(exp(-r_s/r) * exp(r_s/r)) = (r_s/(2r^2)) exp(-r_s/r)
```

So B = A'^2/(4A kappa^2) = (r_s/r^2)^2 exp(-2r_s/r) / (4 exp(-r_s/r) (r_s/(2r^2))^2 exp(-2r_s/r)) = exp(r_s/r). This is self-consistent.

For Schwarzschild: A = 1 - r_s/r, kappa = r_s/(2r^2) [r >> r_s], and B = 1/(1 - r_s/r). Also self-consistent.

**The discriminating question**: What determines kappa(r) from information theory? If kappa is the ASYMPTOTIC surface gravity (kappa = r_s/2), both metrics agree. If kappa is the LOCAL acceleration needed to maintain a static observer at r, then:

- Exponential: a(r) = A'/(2A sqrt(B)) = (r_s/(2r^2)) exp(-r_s/(2r)) / exp(r_s/(2r)) = (r_s/(2r^2)) exp(-r_s/r)
- Schwarzschild: a(r) = (r_s/(2r^2)) / sqrt(1 - r_s/r)

These differ in strong field. The information-theoretic determination of a(r) would discriminate between the two metrics.

---

## 8. The Deep Issue: Time-Space Asymmetry in Quantum Information

### 8.1 Why Time Is Special

The fundamental reason g_00 is accessible to information theory but g_rr is not is:

**Quantum channels are temporal objects.** A channel N: B(H) -> B(H) maps a state at time t_1 to a state at time t_2. The Petz recovery map R recovers the state at t_1 from the state at t_2. The entropy production Sigma = QRE drop quantifies information loss from t_1 to t_2.

There is NO analogous structure for spatial propagation. In quantum mechanics, spatial correlations are properties of the STATE (at a fixed time), not of a DYNAMICAL MAP.

### 8.2 How Relativistic QFT Partially Resolves This

In relativistic QFT, the sharp distinction between time and space is softened:

1. **Modular flow mixes t and x**: The Bisognano-Wichmann modular flow is a Lorentz boost, not a pure time translation. The modular "time" parameter s generates evolution in a direction that combines t and x.

2. **Causal structure replaces time order**: In QFT on curved spacetime, the relevant structure is the causal order (light cone), not the temporal order. Information propagation follows null geodesics, which have both temporal and spatial components.

3. **The Reeh-Schlieder theorem**: Any state can be approximated by acting with operators localized in an arbitrarily small region. This means spatial correlations ARE accessible through local operations -- the distinction between "temporal information loss" and "spatial information spreading" is blurred.

### 8.3 What This Means for the Full Metric

The partial resolution of the time-space asymmetry in QFT suggests:

1. **g_00 is directly accessible**: Through temporal channels (redshift, decoherence)
2. **g_rr is indirectly accessible**: Through the modular flow, which mixes temporal and spatial directions
3. **The full metric requires the full modular structure**: Not just the temperature (which gives g_00) but also the detailed form of the modular Hamiltonian (which encodes g_rr)

This is consistent with the findings from Approaches 1 and 4: the modular Hamiltonian contains information about g_rr, but extracting it requires more detailed knowledge than just the QRE drop (which gives Sigma = -ln(-g_00)).

---

## 9. Honest Conclusions

### 9.1 What Is PROVEN

1. The modular Hamiltonian K depends on both g_00 and g_rr (through the proper distance and induced metric)
2. The Jacobson/Dorau-Much program derives the FULL Einstein equation from entropy/QRE
3. The Einstein equations (given matter content) determine g_rr from g_00
4. The exponential metric is an exact solution of Einstein + phantom scalar
5. In AdS/CFT, the spatial metric IS reconstructible from entanglement data (Czech 2015)
6. The Bisognano-Wichmann modular flow mixes temporal and spatial directions

### 9.2 What Is CONJECTURED (Reasonable but Unproven)

1. The modular flow constrains both g_00 and g_rr simultaneously
2. The surface gravity kappa(r) has an independent information-theoretic expression
3. The phantom scalar field is related to the entropy production Sigma_grav
4. The isometry condition g_00 * g_rr = -1 follows from the isotropy of the modular flow

### 9.3 What Is OPEN (Could Go Either Way)

1. Whether the explicit modular flow calculation in curved spacetime gives an independent constraint on g_rr
2. Whether the inverse spectral problem (metric from modular spectrum) is well-posed
3. Whether flat-space holography can reconstruct the spatial metric
4. Whether the time-space asymmetry in QI is fundamental (g_rr is NOT determinable) or an artifact of current formulations

### 9.4 What Is PROVEN to be IMPOSSIBLE

1. The standard Jacobson argument CANNOT give the exponential metric (gives Schwarzschild)
2. The Petz recovery framework ALONE cannot constrain g_rr (it only defines temporal channels)
3. A scalar entropy production Sigma CANNOT be promoted to a tensor Sigma_ab in any natural way

### 9.5 Classification for Paper 2

**The honest statement for Paper 2**: The tau framework rigorously constrains g_00 = -exp(-Sigma_grav). The spatial metric g_rr requires additional input. The simplest consistent choice is the isometric refractive index (g_00 * g_rr = -1), which is physically motivated by the isotropy of the Petz recovery fidelity but not derived from the framework. The most promising avenue for a genuine derivation is the modular flow approach (Approach 1), which exploits the fact that the Bisognano-Wichmann modular flow mixes temporal and spatial directions, potentially constraining both metric components simultaneously.

---

## 10. Specific Calculations for Future Work

### 10.1 Priority 1: Modular Flow in Static Spherically Symmetric Spacetime

**Goal**: Compute the QRE between states on nested regions in a general static spherically symmetric background.

**Setup**: Consider a free scalar field on the background ds^2 = -A(r) dt^2 + B(r) dr^2 + r^2 dOmega^2. Define the algebra M(r_0) as the algebra of observables in the exterior region r > r_0.

**Computation**: For the vacuum state omega, the modular Hamiltonian of M(r_0) is:

```
K_{r_0} = (2pi/kappa(r_0)) integral_{r_0}^{infinity} sqrt(A(r)) * sqrt(B(r)) * r^2 * T_{00}(r) dr
```

where kappa(r_0) is the surface gravity at r_0 (if there is a horizon) or a regularization parameter (if there is none).

**The key quantity**: The QRE between the vacuum restricted to M(r_0) and M(r_0 + delta):

```
S^rel(omega|_{M(r_0)} || omega|_{M(r_0+delta)}) propto delta * sqrt(A(r_0) * B(r_0)) * r_0^2 * <T_{00}(r_0)>
```

This depends on both A(r_0) and B(r_0). If the tau framework independently fixes A(r_0), then this QRE constrains B(r_0).

### 10.2 Priority 2: Unruh-DeWitt Detector Response

**Goal**: Determine whether the response of an Unruh-DeWitt detector at radius r encodes both A(r) and B(r).

**Setup**: A static UDW detector at radius r with energy gap Omega. Its transition rate is:

```
R(Omega) = integral_{-infinity}^{infinity} e^{-i Omega tau} W(tau, 0) d(tau)
```

where W(x, x') = <0|phi(x)phi(x')|0> is the Wightman function and tau is proper time.

**The Wightman function depends on the full metric**: W(tau, 0) for a static detector at radius r is:

```
W(tau, 0) = <0|phi(t(tau), r) phi(0, r)|0>
```

where t(tau) = tau / sqrt(A(r)) relates proper time to coordinate time. The Wightman function in the Hartle-Hawking state involves a sum over modes:

```
W(tau, 0) = sum_l integral_0^{infinity} |u_{omega l}(r)|^2 e^{-i omega tau/sqrt(A)} (1/(e^{beta omega} - 1) + 1/2) d(omega)/(2pi)
```

The mode functions u_{omega l}(r) satisfy the radial equation that depends on BOTH A(r) and B(r). Therefore, the UDW detector response encodes both metric components.

**Status**: This calculation can be performed numerically for specific metrics. Comparing the exponential and Schwarzschild metrics would give a concrete prediction.

### 10.3 Priority 3: Information-Theoretic Surface Gravity

**Goal**: Derive kappa(r) = A'(r)/(2 sqrt(A(r) B(r))) from information-theoretic principles.

**Approach**: The surface gravity kappa determines the periodicity of the Euclidean section (beta = 2pi/kappa). For a detector at radius r, the Tolman temperature is T(r) = kappa/(2pi sqrt(A(r))). If T(r) can be determined from the QRE of the detector's state (without knowing the metric), then kappa is determined.

**The QRE of a thermalized detector**: A detector at radius r, after thermalization, has the Gibbs state:

```
rho_r = exp(-Omega/(T(r))) / Z = exp(-Omega sqrt(A(r)) / (kappa/(2pi))) / Z
```

The QRE between this state and the vacuum state:

```
D(rho_r || |0><0|) = S(|0><0|) - S(rho_r) + Tr(rho_r ln(|0><0|) - rho_r ln(rho_r))
```

This encodes kappa/sqrt(A). Combined with A from the tau framework, kappa is determined, and hence B.

---

## 11. References

### Primary (Used in Analysis)

1. Bisognano-Wichmann (1975-76): J. Math. Phys. 16, 985; 17, 303
2. Jacobson (1995): PRL 75, 1260 [gr-qc/9504004]
3. Jacobson (2016): PRL 116, 201101 [1505.04753]
4. Dorau-Much (2025): PRL [2510.24491]
5. Witten (2022): JHEP 2022(10), 008 [2112.12828]
6. Chandrasekaran et al. (2023): JHEP 2023(02), 082 [2206.10780]
7. Ryu-Takayanagi (2006): PRL 96, 181602 [hep-th/0603001]
8. Czech et al. (2015): JHEP 2015(10), 175 [1505.05515]
9. Cao-Carroll-Michalakis (2017): PRD 95, 024031 [1606.08444]
10. Van Raamsdonk (2010): Gen. Rel. Grav. 42, 2323 [1005.3035]
11. Swingle (2012): PRD 86, 065007 [1209.3304]
12. Connes (1994): Noncommutative Geometry (Academic Press)
13. Makukov-Mychelkin (2020): Found. Phys. 50, 1346 [2009.08655]
14. Dicke (1957): Rev. Mod. Phys. 29, 363
15. Sewell (1982): Ann. Phys. 141, 201
16. Bianconi (2025): PRD 111, 066001 [2408.14391]
17. Trejo-Calderon (2025): [2504.20457] -- Modular channels, thermal filtering
18. Verlinde (2010): JHEP 04, 029 [1001.0785]
19. Herrera (2020): Entropy 22, 340
20. Wall (2012): PRD 85, 104049 [1105.3445]

### Secondary (Contextual)

21. Petz (1988): Q. J. Math. 39, 97
22. JRSWW (2018): Ann. Henri Poincare 19, 2955
23. Padmanabhan (2010): Rep. Prog. Phys. 73, 046901 [0911.5004]
24. Boonserm et al. (2018): PRD 98, 084048 [1805.03781]
25. Xu-Wu (2023): [2305.01330] -- Metric reconstruction from entanglement and complexity
26. Kinoshita et al. (2024): [2410.07587] -- Spin systems as QFT simulators
27. Hollands-Longo (2021): [2107.06787] -- Relative entropy in curved spacetimes
28. Bao et al. (2020): JHEP 2020(12), 033 -- Tensor Radon transform
29. Cao-Carroll (2018): PRD 97, 086003 [1712.02803]

---

**Last updated**: 2026-03-10
**Supersedes**: paper2_full_metric_derivation.md (same date, earlier version)
