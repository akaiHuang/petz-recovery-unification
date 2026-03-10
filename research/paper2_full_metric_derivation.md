# From g00 to Full Metric: Can Information Theory Determine g_rr?

## Date: 2026-03-10
## Status: Comprehensive research report
## Purpose: Investigate whether information theory can determine spatial metric components beyond g00

---

## 1. Problem Statement

### 1.1 What Paper 2 Currently Achieves

Paper 2's triple identification constrains only the temporal metric component:

```
sqrt(-g_00) = exp(-Sigma_grav/2) = F_bound
```

With Sigma_grav = r_s/r, this gives:

```
g_00 = -exp(-r_s/r)
```

This is well-supported by three independent first-principles routes (modular flow, gravitational Landauer, bosonic channel model) and is the strongest result in Paper 2.

### 1.2 What Paper 2 Assumes

The full exponential metric in isotropic coordinates is:

```
ds^2 = -exp(-r_s/r) c^2 dt^2 + exp(+r_s/r) (dr^2 + r^2 dOmega^2)
```

The spatial component g_rr = exp(+r_s/r) is obtained in the paper via the statement (Section II.D, line ~411):

> "For the spatial components, we invoke the isotropy condition natural in the information-theoretic setting---the recovery fidelity is a scalar quantity and does not depend on the direction of signal propagation---together with asymptotic flatness, obtaining the full line element: g_rr = e^{+r_s/r} by the requirement g_00 * g_rr = -1 in isotropic coordinates."

### 1.3 The Core Question

Is g_rr = exp(+r_s/r) derivable from information theory, or is it an independent assumption? If the latter, how much does Paper 2 depend on it, and is the assumption justified?

### 1.4 Why This Matters

The spatial metric component g_rr determines:
- Photon sphere location
- Shadow size
- ISCO radius
- Geodesic structure
- Light deflection at strong field
- Wormhole throat geometry

If g_rr cannot be derived from information theory, then Paper 2's strong-field predictions (which depend on the full metric) rest on an assumption that goes beyond the information-theoretic framework.

---

## 2. Literature Survey

### 2.1 Programs That Derive Full Einstein Equations (Including All Metric Components)

#### Jacobson 1995 (PRL 75, 1260)

**Method**: delta-Q = T delta-S applied to all local Rindler horizons through each spacetime point.

**Key insight**: By demanding that the Clausius relation hold for ALL local Rindler causal horizons (not just one), and using the Raychaudhuri equation to relate the area change delta-A to the Ricci tensor R_ab xi^a xi^b, Jacobson obtains:

```
alpha <T_ab> xi^a xi^b = R_ab xi^a xi^b
```

Since this holds for ALL boost directions xi^a, the tensor equation follows:

```
R_ab - (R/2) g_ab + Lambda g_ab = 8pi G <T_ab>
```

**Crucially**: This gives the FULL tensorial Einstein equation, including all spatial components. The key step is that the argument works for boosts in ALL directions, not just time. The Raychaudhuri equation converts the scalar entropy-area relation into a tensor equation because it relates a scalar (expansion) to a contraction of R_ab with the null vector, and varying over all null vectors recovers the full tensor.

**Relevance to our problem**: The Jacobson program shows that thermodynamic/entropic arguments CAN determine all metric components -- but only if the argument is applied to ALL local horizons (all directions), not just radial ones.

#### Jacobson 2016 (PRL 116, 201101)

**Method**: Entanglement equilibrium -- the entanglement entropy in small geodesic balls (causal diamonds) is maximized at fixed volume in the vacuum.

**Result**: For conformal quantum fields, vacuum entanglement is stationary iff the Einstein equation holds. This is again the FULL tensorial equation.

**How spatial components emerge**: The causal diamond analysis probes entanglement across surfaces of all orientations. The conformal Killing vector of the diamond generates flow in all spacetime directions, not just time. By requiring entanglement equilibrium for causal diamonds of all shapes and orientations, the full R_ab is constrained.

**Limitation**: This is a consistency condition (entanglement equilibrium <=> Einstein equation), not a derivation that starts from a metric ansatz and determines it. It cannot select between Schwarzschild and exponential metrics -- it says whatever the correct metric is, it must satisfy the Einstein equation (or its modification).

#### Dorau-Much 2025 (PRL, arXiv:2510.24491)

**Method**: QRE between vacuum and coherent excitations on bifurcate Killing horizons, combined with Raychaudhuri equation and Bekenstein-Hawking area law.

**Result**: The full semiclassical Einstein equation:

```
R_ab - (R/2) g_ab + Lambda g_ab = alpha <:T_ab:>
```

**How all components are recovered**: The critical step (their eq. 17) starts as:

```
alpha <:T_ab:> xi^a xi^b = R_ab xi^a xi^b
```

This is a scalar equation for a specific Killing direction xi. By invoking the equivalence principle (the argument works for ANY local Rindler frame, i.e., for ANY xi^a), and applying local energy-momentum conservation (divergence-free condition), the FULL tensor equation is forced.

**Key point**: The tensorialization step -- going from R_ab xi^a xi^b = alpha T_ab xi^a xi^b (scalar, for one direction) to R_ab = alpha T_ab + N g_ab (tensor, for all directions) -- relies on:
1. The argument holding for ALL Killing vector directions
2. Local conservation: nabla^a T_ab = 0

This is why a SINGLE direction (like radial) is insufficient.

### 2.2 Entropic Gravity Programs

#### Verlinde 2010 (JHEP 04, 029)

**What it constrains**: Verlinde derives Newton's law F = ma and the gravitational force from holographic screens. The argument gives the FORCE, not the metric directly.

**Spatial metric**: Verlinde does NOT derive g_rr from entropy. His framework constrains the Newtonian potential (hence g_00 in weak field) through the entropic force argument, but the spatial metric is not addressed.

**Verlinde 2016**: Extends to the cosmological setting and derives apparent dark matter effects. Again, the focus is on the gravitational potential (related to g_00), not spatial metric components. The de Sitter entanglement entropy gives modifications to the potential, but g_rr is not independently constrained.

#### Padmanabhan 2010 (Rep. Prog. Phys. 73, 046901)

**Method**: Gravity as the free energy of spacetime. The gravitational action = surface entropy - bulk energy. Variation gives Einstein equations.

**Spatial metric**: Padmanabhan's approach recovers the FULL Einstein equations (all components), because the action principle naturally gives all components through variation of the full metric g_ab. The entropy functional involves the Lanczos-Lovelock Lagrangian, which depends on all metric components.

**Key distinction**: Padmanabhan does not derive a specific metric from information theory. He provides an entropic interpretation of the gravitational ACTION, which then determines the metric through field equations. The spatial components are determined by the field equations, not directly by entropy.

### 2.3 Holographic Approaches

#### Ryu-Takayanagi (2006) and Spatial Metric Reconstruction

The RT formula S_A = Area(gamma_A) / (4G) relates entanglement entropy to the area of minimal surfaces in the bulk. This constrains spatial geometry because:

1. The area of gamma_A depends on the FULL bulk metric (including spatial components)
2. By varying the boundary region A, one probes different minimal surfaces, which encode different aspects of the spatial metric

**Cao-Carroll-Michalakis 2017 (PRD 95, 024031)**: "Space from Hilbert Space"

They show that for "redundancy-constrained states" in Hilbert space:
1. Entanglement entropy across codimension-one surfaces defines areas
2. A Radon transform converts area data into a spatial metric
3. Under appropriate assumptions, the spatial geometry and Einstein's equation emerge in weak field

**Key result**: This is the most explicit demonstration that entanglement entropy CAN determine spatial metric components. However:
- It requires a specific class of quantum states
- It works in the bulk (AdS) setting, not asymptotically flat
- The Radon transform reconstruction is global, not local
- It has not been applied to derive g_rr for a specific spacetime like Schwarzschild or exponential

#### Van Raamsdonk (2010)

"Building up Spacetime with Quantum Entanglement": Argues that spatial connectivity = entanglement. Disconnecting two regions reduces entanglement to zero. This constrains the spatial metric qualitatively but not quantitatively.

### 2.4 Bianconi GfE (PRD 111, 066001, 2025)

**Framework**: The entropic action is the quantum relative entropy between the spacetime metric g_ab and the metric induced by matter fields G_ab:

```
S_GfE = integral sqrt(|g|) Tr(g_tilde ln(G_tilde^{-1})) d^4x
```

**Variation**: Variation with respect to g_ab gives modified Einstein equations that reduce to standard Einstein equations at low coupling. The variation is with respect to the FULL metric, so in principle all components (g_00, g_rr, etc.) are determined.

**Vacuum solution**: Thattarampilly, Zheng & Kakkat (arXiv:2602.13694, Feb 2026) solved the GfE equations perturbatively for static spherically symmetric vacuum:

```
A(r) = 1 - r_s/r - beta r_s^2 / (48 r^4)
B(r) = (1 - r_s/r)^{-1} + beta r_s^2 / (12 r^4)
```

This is Schwarzschild + r^{-4} corrections, NOT the exponential metric. The exponential metric has r^{-2} corrections to Schwarzschild. The GfE does NOT produce the exponential metric.

**Assessment**: The GfE framework CAN determine all metric components through its field equations, but the specific solution it gives is NOT the exponential metric.

---

## 3. Route Analysis (Each Route's Feasibility)

### Route A: The Isometry Condition g_00 * g_rr = -1

#### Statement of the Condition

In isotropic coordinates, if the spatial metric is conformally flat and determined by a SINGLE function:

```
ds^2 = -f(r) dt^2 + h(r) (dr^2 + r^2 dOmega^2)
```

then there exist multiple possible relationships between f and h. The specific condition g_00 * g_rr = -1 (i.e., f * h = 1) is a CHOICE, not a necessity.

#### When Does g_00 * g_rr = -1 Hold?

**In Schwarzschild coordinates**: g_tt * g_rr = -1 holds EXACTLY for the Schwarzschild vacuum solution. This follows from the vacuum Einstein equations R_tt = 0 and R_rr = 0, which together imply:

```
d/dr [ln(-g_tt) + ln(g_rr)] = 0  =>  g_tt * g_rr = const
```

Boundary conditions at infinity fix the constant to -1.

**CRITICAL POINT**: This g_tt * g_rr = -1 result is a consequence of the VACUUM Einstein equations. It holds for:
- Schwarzschild (exactly)
- Any vacuum solution of Einstein equations in Schwarzschild coordinates

It does NOT automatically hold for:
- Solutions with matter (the exponential metric requires phantom scalar matter)
- Solutions in isotropic coordinates (the transformation mixes the components)
- Modified gravity theories

#### For the Exponential Metric Specifically

In isotropic coordinates, the exponential metric has:

```
g_00 = -exp(-r_s/r)
g_rr = +exp(+r_s/r)
```

So g_00 * g_rr = -exp(-r_s/r) * exp(+r_s/r) = -1.

But wait -- this is in ISOTROPIC coordinates. The condition g_00 * g_rr = -1 in isotropic coordinates is NOT the same as in Schwarzschild coordinates.

**For Schwarzschild in isotropic coordinates**:

```
g_00^{Schw} = -((1 - r_s/(4r))/(1 + r_s/(4r)))^2
g_rr^{Schw} = (1 + r_s/(4r))^4
g_00 * g_rr = -((1 - r_s/(4r))/(1 + r_s/(4r)))^2 * (1 + r_s/(4r))^4
            = -(1 - r_s/(4r))^2 * (1 + r_s/(4r))^2
            = -(1 - r_s^2/(16r^2))^2
            != -1
```

So **Schwarzschild does NOT satisfy g_00 * g_rr = -1 in isotropic coordinates!**

The condition g_00 * g_rr = -1 in isotropic coordinates is SPECIFIC to the exponential metric (and other one-function metrics where the spatial part is the reciprocal of the temporal part).

#### Information-Theoretic Motivation

Paper 2's argument is: "the recovery fidelity is a scalar quantity and does not depend on the direction of signal propagation."

This is physically reasonable but mathematically incomplete:
1. The recovery fidelity F is indeed a scalar (direction-independent)
2. But the ENTROPY PRODUCTION Sigma could be direction-dependent in principle
3. The connection between Sigma and the metric component goes through a specific channel model
4. For time-direction: Sigma_t = -ln(-g_00)
5. For space-direction: what channel gives Sigma_r = ln(g_rr)?

**The gap**: Even if the fidelity is direction-independent, the entropy production for spatial propagation is a DIFFERENT quantity from the entropy production for temporal propagation. They need not be related by g_00 * g_rr = -1.

#### Viability Assessment

| Aspect | Assessment |
|--------|-----------|
| Mathematical well-definedness | YES -- the condition is clear |
| Derivation from information theory | NO -- it is an ADDITIONAL ASSUMPTION |
| Physical motivation | PARTIAL -- scalar fidelity argument is suggestive but not rigorous |
| Self-consistency | YES -- when combined with g_00, it gives a consistent metric |
| Uniqueness | NO -- other relations (e.g., g_rr = (-g_00)^{alpha} for alpha != -1) are equally consistent |

**Verdict**: The isometry condition is the SIMPLEST choice but cannot be derived from information theory alone. It is an independent physical assumption.

---

### Route B: Spatial Petz Recovery

#### Concept

Define a "spatial recovery" problem: instead of asking "can we recover a signal that climbed out of a potential well?" (temporal), ask "can we recover spatial information that was spread across a curved geometry?"

For a signal propagating in the radial direction at fixed time, the relevant metric component is g_rr. The spatial entropy production would be:

```
Sigma_spatial = -ln(g_rr)  [by analogy with Sigma_temporal = -ln(-g_00)]
```

This would give g_rr = exp(-Sigma_spatial). If Sigma_spatial = -r_s/r, then g_rr = exp(+r_s/r).

#### Analysis

**Problem 1: What is the spatial channel?**

For temporal propagation, the channel is clear: gravitational redshift degrades signals as they climb out. The Tolman relation gives the temperature gradient. The bosonic loss channel models the energy loss.

For spatial propagation at fixed time, the situation is DIFFERENT:
- There is no gravitational "energy loss" in the spatial direction
- A ruler at radius r has the same proper length regardless of where you view it from
- The spatial metric describes the GEOMETRY (proper distances), not an information degradation channel

**Conceptual problem**: Spatial metric components describe distances, not information loss channels. The temporal metric component g_00 has a natural interpretation as information degradation (redshift = loss channel). The spatial component g_rr describes how distances are stretched but does NOT naturally correspond to an information loss process.

**Problem 2: Directionality of recovery**

The Petz recovery map is defined for a quantum channel N. For gravitational redshift, N maps the state at radius r to the state at infinity. This is naturally a temporal/radial process.

A "spatial Petz recovery" would require a channel that maps states from one spatial location to another AT THE SAME TIME. In GR, simultaneity is coordinate-dependent, making this construction fundamentally observer-dependent.

**Problem 3: No analogue of Tolman relation for spatial direction**

The Tolman relation T(r) = T_inf / sqrt(-g_00) provides the temperature gradient that drives the temporal entropy production. There is NO spatial analogue: there is no "spatial temperature" that varies with angle or direction.

#### Viability Assessment

| Aspect | Assessment |
|--------|-----------|
| Mathematical formulation | UNCLEAR -- what is the spatial channel? |
| Physical motivation | WEAK -- spatial metric != information loss |
| Existing literature | NONE -- no one has attempted this |
| Connection to known results | NONE -- no analogue of Tolman relation |
| Self-consistency | UNKNOWN |

**Verdict**: Spatial Petz recovery is conceptually attractive but fundamentally problematic. The spatial metric does not naturally correspond to an information degradation channel. This route is NOT viable without major new insights.

---

### Route C: Jacobson Route (Entropy -> Field Equations -> Full Metric)

#### How It Would Work

1. Start with the Jacobson/Dorau-Much program: entropy/QRE on local horizons gives the FULL Einstein equation
2. The Einstein equation (plus boundary conditions) determines the full metric including g_rr
3. For vacuum: R_ab = 0 determines both g_00 and g_rr
4. For matter (phantom scalar): R_ab - (R/2)g_ab = 8pi G T_ab determines both

#### Analysis

**The good**: The Jacobson-Dorau-Much program DOES derive the full tensorial Einstein equation from information-theoretic principles. The key step is that the argument works for ALL boost directions (not just radial), converting the scalar QRE relation into a tensor equation. This then determines all metric components through the field equations.

**The problem for Paper 2**: The Jacobson program derives the EINSTEIN equations, whose vacuum solution is SCHWARZSCHILD, not exponential. To get the exponential metric, one needs:
- Either modified field equations (but the Jacobson program gives standard GR)
- Or matter content (phantom scalar, but this is not derived from information theory)

**The chain breaks**: Information theory -> Einstein equations -> Schwarzschild (NOT exponential)

If we want information theory -> exponential metric, we CANNOT go through the Jacobson route, because the Jacobson route gives the Einstein equations, which give Schwarzschild.

**A hybrid approach**: One could argue:
1. Jacobson/Dorau-Much: QRE -> Einstein equations (determines tensorial structure)
2. Paper 2's identification: Sigma_grav = r_s/r -> exponential g_00
3. Combine: the exponential g_00 is sourced by phantom matter, and the field equations then determine g_rr

But this requires knowing the matter content, which is derived from the exponential metric itself (circular).

#### Viability Assessment

| Aspect | Assessment |
|--------|-----------|
| Derives full metric? | YES (through field equations) |
| Gives exponential metric? | NO -- gives Schwarzschild |
| Information-theoretic basis? | YES -- strong (Jacobson, Dorau-Much) |
| Self-consistency? | YES (for Schwarzschild) |
| Compatible with Paper 2? | TENSION -- different metric |

**Verdict**: The Jacobson route CAN derive all metric components, but it gives the wrong metric (Schwarzschild, not exponential). This highlights a deep tension: the most rigorous information-theoretic approach to full metric derivation leads to Schwarzschild, not exponential.

---

### Route D: Dorau-Much Route (QRE -> Semiclassical Einstein -> All Components)

#### How It Works

Dorau & Much (2025 PRL) derive:

```
R_ab - (R/2) g_ab + Lambda g_ab = alpha <:T_ab:>_omega
```

This is the full tensorial equation. The derivation uses:
1. QRE on bifurcate Killing horizon = energy flux across horizon
2. Raychaudhuri equation converts energy flux to area change
3. Bekenstein-Hawking converts area to entropy
4. The equivalence principle makes the argument local and directional: valid for ALL Rindler frames
5. Energy-momentum conservation forces the scalar relation into a tensor equation

#### What It Constrains

The Dorau-Much result constrains the FULL metric (all components) indirectly: the Einstein equation, given matter content and boundary conditions, determines the metric uniquely (up to coordinate choice).

#### Relation to Paper 2's Problem

The Dorau-Much computation of fractional QRE loss = r_s/r (our Route A in modular flow notes) only probes the TEMPORAL direction (energy flux through the horizon, related to the boost Killing vector). This gives Sigma = -ln(-g_00).

To get information about spatial components, one would need to:
1. Apply the Dorau-Much argument to spatial boost directions
2. Compute the QRE loss for signals propagating spatially
3. Show this gives Sigma_spatial = ln(g_rr)

But the Dorau-Much computation uses the BOOST Killing vector, which is inherently temporal (it generates time translations for the accelerated observer). There is no natural "spatial boost" analogue.

**However**, the FULL Einstein equation that Dorau-Much derives DOES constrain all components. Once you have R_ab = alpha T_ab + N g_ab, the spatial components of R_ab are constrained. This works through the field equations, not through direct spatial entropy production.

#### Viability Assessment

| Aspect | Assessment |
|--------|-----------|
| Derives full metric? | YES (through Einstein equations) |
| Gives exponential metric? | NO (gives Schwarzschild in vacuum) |
| Information-theoretic basis? | YES -- the strongest available |
| Direct spatial constraint? | NO -- goes through field equations |

**Verdict**: Same as Route C. The Dorau-Much route gives full metric but via Einstein equations -> Schwarzschild.

---

### Route E: Bianconi GfE (Metric = Density Matrix)

#### How It Works

In Bianconi's framework, the metric g_ab IS (plays the role of) a density matrix. The entropic action is the quantum relative entropy between the spacetime metric and the metric induced by matter:

```
S_GfE = integral sqrt(|g|) D_QRE(g_tilde || G_tilde^{-1}) d^4x
```

Variation with respect to g_ab gives modified Einstein equations that determine ALL metric components.

#### What the Vacuum Solution Gives

The first explicit vacuum solution (Thattarampilly et al. 2026):
- Schwarzschild + r^{-4} perturbative corrections
- NOT the exponential metric
- The r^{-4} corrections are qualitatively different from the r^{-2} corrections that distinguish exponential from Schwarzschild

#### Can It Give the Exponential Metric?

**At perturbative level**: NO. The r^{-4} structure is incompatible with the exponential metric's r^{-2} structure.

**At full nonlinear level**: UNKNOWN. The full nonlinear GfE equations have not been solved. In principle, a non-perturbative solution could differ dramatically from the perturbative one. But there is no evidence for this.

**Fundamental issue**: The GfE action involves eigenvalues of the metric between bivectors. The exponential metric has specific eigenvalue structure: lambda_0 = exp(-r_s/r), lambda_1 = lambda_2 = lambda_3 = exp(+r_s/r). Whether this eigenvalue pattern is an extremum of the GfE action is an open question.

#### Viability Assessment

| Aspect | Assessment |
|--------|-----------|
| Determines all g_ab? | YES (through field equations) |
| Gives exponential metric? | NO (perturbative result gives Schwarzschild + r^{-4}) |
| Information-theoretic? | YES (QRE-based action) |
| Status | CONFIRMED NEGATIVE at perturbative level |

**Verdict**: The GfE framework determines all metric components but gives the WRONG answer for the exponential metric. This is a significant negative result for Paper 2.

---

### Route F: Entropy Density Tensor Sigma_ab

#### Concept

Define a tensor-valued entropy production:

```
Sigma_ab = some tensorial generalization of Sigma
```

such that:
- Sigma_00 = r_s/r (gives g_00 via exp(-Sigma_00))
- Sigma_rr = -r_s/r (gives g_rr via exp(-Sigma_rr) = exp(+r_s/r))
- Sigma_theta_theta, Sigma_phi_phi follow from spherical symmetry

#### Analysis

**The mathematical question**: Can the scalar entropy production Sigma be promoted to a tensor?

**Existing candidates**:

1. **Modular Hamiltonian**: The modular Hamiltonian K is associated with a specific region and state. It is an operator, not a tensor field on spacetime. It does not naturally give a tensor Sigma_ab.

2. **Relative modular operator**: Delta_{phi/omega} is an operator. Its logarithm gives the relative entropy S(omega || phi) = -<psi_omega| ln Delta |psi_omega>, which is a scalar.

3. **Fisher information metric**: The quantum Fisher information g_F^{ab} IS a tensor on state space. It could serve as a "metric on information geometry." But this is a metric on the SPACE OF STATES, not on spacetime. The connection to spacetime metric is indirect and not established.

4. **Stress-energy tensor of entropy production**: One could define:

```
S_ab = (1/alpha) R_ab  [using Einstein equations R_ab = alpha T_ab + ...]
```

This IS a tensor. For Schwarzschild: R_ab = 0 (vacuum), so S_ab = 0. This does not give Sigma_ab = diag(r_s/r, -r_s/r, ...).

5. **Weyl tensor components**: The Weyl tensor C_{abcd} contains the tidal information in vacuum. For Schwarzschild: C_{trtr} ~ r_s/r^3. But this is local (r^{-3}), not integrated (r^{-1}).

**The fundamental obstruction**: Entropy production Sigma is fundamentally a SCALAR quantity: it is the QRE drop D(rho||sigma) - D(N(rho)||N(sigma)), which is a real number. Promoting it to a tensor requires additional structure that does not exist in the standard JRSWW framework.

**A speculative direction**: In the Bianconi framework, the metric IS a density matrix. The "entropy" of the metric involves eigenvalues of g_tilde, which has a natural tensorial structure. One could define:

```
Sigma_ab = ln(-g_ab / eta_ab)  [ratio to flat metric]
```

For the exponential metric: Sigma_00 = -r_s/r, Sigma_rr = r_s/r, giving g_ab = eta_ab exp(Sigma_ab). But this is TAUTOLOGICAL -- it defines Sigma_ab from g_ab rather than deriving g_ab from Sigma_ab.

#### Viability Assessment

| Aspect | Assessment |
|--------|-----------|
| Mathematical well-definedness | PROBLEMATIC -- Sigma is naturally scalar |
| Physical motivation | WEAK -- no clear physical meaning for Sigma_rr |
| Existing literature | NONE -- no one has proposed this |
| Self-consistency | TAUTOLOGICAL if defined as ln(g_ab/eta_ab) |
| Predictive power | NONE unless independently motivated |

**Verdict**: The entropy tensor approach is currently ill-defined and tautological. Without a physical mechanism that produces direction-dependent entropy production, this route is NOT viable.

---

## 4. The Isometry Condition g_00 * g_rr = -1

### 4.1 Mathematical Status

In isotropic coordinates for a static spherically symmetric metric:

```
ds^2 = -A(r) dt^2 + B(r) (dr^2 + r^2 dOmega^2)
```

The condition A * B = 1 (equivalently g_00 * g_rr = -1) is NOT generic. It holds only for specific metrics.

### 4.2 Which Metrics Satisfy It?

| Metric | g_00 * g_rr (isotropic) | Satisfies? |
|--------|------------------------|------------|
| Minkowski | (-1)(+1) = -1 | YES |
| Exponential | -exp(-r_s/r) * exp(+r_s/r) = -1 | YES |
| Schwarzschild (isotropic) | -((1-m/r)/(1+m/r))^2 * (1+m/r)^4 = -(1-m^2/r^2)^2 | NO |
| de Sitter | Depends on coordinates | NO (generically) |
| Reissner-Nordstrom (isotropic) | | NO |

The condition is SPECIAL. It holds for metrics where the spatial factor is exactly the reciprocal of the temporal factor.

### 4.3 Physical Meaning

The condition A * B = 1 means:

```
ds^2 = -A(r) dt^2 + (1/A(r)) (dr^2 + r^2 dOmega^2)
```

This has a natural interpretation in terms of the coordinate speed of light. For radial light (ds^2 = 0, dOmega = 0):

```
(dr/dt)^2 = A(r) / B(r) = A(r)^2
=> dr/dt = A(r)
```

So the coordinate speed of light is sqrt(-g_00) in BOTH temporal and spatial directions. This is the "isotropic" property in the strongest sense: the effective refractive index n = c/c_eff = 1/sqrt(-g_00) is the SAME in all directions.

### 4.4 Connection to Dicke's Refractive Index

Dicke (1957) showed that GR can be reformulated as flat spacetime with a gravitational refractive index:

```
n(r) = c/c_eff(r) = 1/sqrt(-g_00(r))
```

For the exponential metric: n(r) = exp(r_s/(2r)).

The condition g_00 * g_rr = -1 means that this refractive index is ISOTROPIC -- the same for light propagating in any direction. This is the most natural condition in the Dicke picture.

### 4.5 Is Isotropic Refractive Index Derivable from Information Theory?

**Argument in favor**: The Petz recovery fidelity F is a scalar quantity. It does not depend on the direction in which information propagates. If we identify:

```
c_eff/c = F_bound = exp(-Sigma/2)
```

then c_eff should be direction-independent, giving isotropic refraction, hence g_00 * g_rr = -1.

**Argument against**: The entropy production Sigma_grav is defined for a specific quantum channel (radial propagation from r to infinity). A channel for LATERAL propagation (tangential to a sphere at radius r) would be fundamentally different -- there is no gravitational potential difference for tangential motion. So the entropy production for lateral propagation is ZERO, not r_s/r.

This means the effective speed of light is NOT isotropic if Sigma is channel-dependent:
- Radial: Sigma_radial = r_s/r => c_eff^radial/c = exp(-r_s/(2r))
- Tangential: Sigma_tangential = 0 => c_eff^tangential/c = 1

This would give ANISOTROPIC light propagation, inconsistent with the isometry condition.

**Resolution attempt**: The entropy production should be understood as the cost of maintaining a static observer at radius r (the "integrated informational cost"), not the cost of a specific propagation direction. In this interpretation, Sigma_grav = r_s/r is a PROPERTY OF THE LOCATION, not of the propagation direction. Both radial and tangential light at radius r experience the same informational cost.

This interpretation is consistent with the Tolman relation (temperature depends on position, not direction) and with the Unruh effect (the vacuum state appears thermal to the accelerated observer regardless of measurement direction).

**Assessment**: The isotropic refractive index is CONSISTENT with the information-theoretic framework but not DERIVABLE from it without the additional assumption that entropy production is location-dependent rather than direction-dependent.

---

## 5. Honest Assessment

### 5.1 What CAN Be Derived from Information Theory

| Quantity | Status | Method | Confidence |
|----------|--------|--------|------------|
| g_00 = -exp(-r_s/r) | DERIVED | Triple identification + 3 first-principles routes | MEDIUM-HIGH |
| g_rr = exp(+r_s/r) | ASSUMED | Isometry condition g_00 * g_rr = -1 | ASSUMPTION |
| Full Einstein equation | DERIVED | Jacobson/Dorau-Much programs | HIGH (but gives Schwarzschild) |
| No-horizon theorem | DERIVED | JRSWW bound + finite Sigma | HIGH (model-independent) |

### 5.2 The Fundamental Tension

There is a deep tension between two approaches:

**Approach 1 (Paper 2)**: Information theory -> Sigma_grav = r_s/r -> g_00 = exp(-r_s/r) -> [assume g_00*g_rr = -1] -> exponential metric

**Approach 2 (Jacobson/Dorau-Much)**: Information theory -> full Einstein equations -> Schwarzschild metric (including g_rr)

Approach 1 gives the exponential metric but requires the isometry assumption. Approach 2 gives all metric components from information theory but produces Schwarzschild, not exponential.

**The resolution**: These are DIFFERENT uses of information theory:
- Approach 1 identifies the metric with the Petz recovery fidelity bound (specific to the metric value, not the field equations)
- Approach 2 derives the field equations from entanglement/QRE, then solves them

The tension is not a contradiction but a different choice of what is "fundamental": the field equations (Jacobson) vs. the metric directly (Paper 2). They give different predictions in strong field, which is TESTABLE.

### 5.3 What g_rr Cannot Be

Even without deriving g_rr, we can constrain it:

1. **Asymptotic flatness**: g_rr -> 1 as r -> infinity
2. **Weak-field agreement**: g_rr = 1 + r_s/r + O((r_s/r)^2) (from PPN agreement)
3. **Spherical symmetry**: g_rr depends only on r
4. **Positive definiteness**: g_rr > 0 (signature requirement)
5. **Regularity**: g_rr should not have singularities (for the exponential metric philosophy)

These constraints are satisfied by BOTH:
- g_rr = exp(+r_s/r) (exponential, Paper 2)
- g_rr = (1 + r_s/(4r))^4 (Schwarzschild in isotropic coordinates)
- Many other functions

### 5.4 How Much Do Strong-Field Predictions Depend on g_rr?

| Observable | Depends on g_00? | Depends on g_rr? | Sensitivity to g_rr |
|-----------|-----------------|-----------------|---------------------|
| Gravitational redshift z | YES (directly) | NO | NONE |
| Photon sphere location | YES | YES | HIGH |
| Shadow size | YES | YES | HIGH |
| ISCO radius | YES | YES | HIGH |
| QNM frequencies | YES | YES | HIGH |
| Echo time delay | YES | YES | HIGH |
| Perihelion precession | YES | YES | MODERATE |
| Light deflection | YES | YES | MODERATE |

**The honest conclusion**: Paper 2's prediction of g_00 (no horizons, tau < 1) is robust and independent of g_rr. But ALL strong-field observational predictions (shadow, ISCO, QNM, echoes) depend on g_rr and therefore rest on the ASSUMED isometry condition.

### 5.5 What This Means for Paper 2

**Layer 1 (no-horizon argument)**: UNAFFECTED. This depends only on g_00 (more precisely, only on Sigma being finite). It does not need g_rr at all.

**Layer 2 (exponential g_00)**: UNAFFECTED. The derivation of g_00 = -exp(-r_s/r) via three independent routes does not require g_rr.

**Layer 3 (full exponential metric)**: WEAKENED. The g_rr component is ASSUMED via the isometry condition, not derived. All strong-field predictions that depend on the full metric (shadow, ISCO, QNM, echoes) are conditional on this assumption.

**Paper 2's predictions should be classified**:
- **Robust predictions** (independent of g_rr): no horizons, tau < 1, redshift at neutron star surface
- **Conditional predictions** (depend on g_rr assumption): photon sphere location, shadow size, ISCO, QNM, echo time delay

---

## 6. Recommendation for Paper 2

### 6.1 Honest Presentation

The current text (Section II.D, line ~411) should be expanded to honestly state:

```
The triple identification constrains g_00 but not the spatial
components of the metric. The full exponential metric
(Eq. 3) requires an additional assumption: in the Dicke
reformulation, gravity acts as a refractive index, and the
most natural information-theoretic assumption is that this
index is isotropic (direction-independent). Isotropy, together
with g_00 = -exp(-r_s/r) and asymptotic flatness, uniquely
determines g_rr = exp(+r_s/r).

This assumption is physically motivated---the Petz recovery
fidelity is a scalar and does not depend on the direction of
signal propagation---but it is not derived from information
theory. An alternative approach that does derive all metric
components (the Jacobson-Dorau-Much program) gives the
Einstein equations, whose vacuum solution is Schwarzschild,
not exponential. The two approaches represent different
choices of what information-theoretic quantity is fundamental:
the fidelity bound (our approach) versus the field equations
(Jacobson's approach). Strong-field observations will
discriminate between them.
```

### 6.2 Classification of Predictions

Add a table distinguishing robust vs. conditional predictions:

```
Table: Classification of predictions by theoretical robustness

Robust (g_00 only):
- No event horizons for finite-mass objects (Layer 1)
- Gravitational redshift at neutron-star compactness: 19% smaller than Schwarzschild
- tau < 1 everywhere

Conditional (full metric, assumes isotropic refraction):
- Photon sphere at r_ph = r_s (vs 1.87 m for Schwarzschild)
- Shadow 4.6% larger than Schwarzschild
- ISCO 5.6% larger
- QNM frequency shift ~4.4% (eikonal)
- Echo time delay ~4.2 r_s/c
```

### 6.3 The Limitations Section

The existing limitations section (Sec. V.C) already contains the correct statement:

> "One cannot invert [the relation] to obtain the full metric tensor g_ab from tau alone: tau constrains only g_00, not the spatial components. A complete reconstruction of g_ab would require additional information-theoretic input (e.g., directional recovery fidelities)."

This is honest and correct. The recommendation is to KEEP this statement and potentially strengthen it by explicitly noting that g_rr = exp(+r_s/r) is assumed, not derived.

### 6.4 Future Direction

The most promising path to deriving g_rr from information theory would be:

1. **Directional recovery fidelities**: Define separate Petz recovery problems for radial vs. tangential propagation. If these give different Sigma, g_rr != exp(+r_s/r). If they give the same Sigma, the isometry condition follows.

2. **Holographic reconstruction**: Use the Cao-Carroll-Michalakis framework to reconstruct the full spatial metric from entanglement entropy data across multiple surfaces.

3. **Modified Jacobson**: Modify the Jacobson entropy-area relation to use Sigma_grav = r_s/r instead of the standard area law. This would give modified field equations whose vacuum solution might be the exponential metric.

4. **Explicit gravitational channel**: Construct the full CPTP map N_grav using the crossed-product algebra, including its action on spatial degrees of freedom. This would determine g_rr from the channel's properties.

---

## 7. Complete Reference List

### Primary Papers

1. **Jacobson 1995**: T. Jacobson, "Thermodynamics of Spacetime: The Einstein Equation of State," PRL 75, 1260. [arXiv:gr-qc/9504004]

2. **Jacobson 2016**: T. Jacobson, "Entanglement Equilibrium and the Einstein Equation," PRL 116, 201101. [arXiv:1505.04753]

3. **Dorau-Much 2025**: P. Dorau and A. Much, "From Quantum Relative Entropy to the Semiclassical Einstein Equations," PRL (2026). [arXiv:2510.24491]

4. **Bianconi 2025**: G. Bianconi, "Gravity from entropy," PRD 111, 066001. [arXiv:2408.14391]

5. **Thattarampilly-Zheng-Kakkat 2026**: "First spherically symmetric GfE solution," [arXiv:2602.13694]

6. **Makukov-Mychelkin 2020**: M. A. Makukov and E. G. Mychelkin, "Triple Path to the Exponential Metric," Found. Phys. 50, 1346. [arXiv:2009.08655]

7. **Verlinde 2010**: E. Verlinde, "On the Origin of Gravity and the Laws of Newton," JHEP 04, 029. [arXiv:1001.0785]

8. **Verlinde 2016**: E. Verlinde, "Emergent Gravity and the Dark Universe." [arXiv:1611.02269]

9. **Padmanabhan 2010**: T. Padmanabhan, "Thermodynamical Aspects of Gravity: New insights," Rep. Prog. Phys. 73, 046901. [arXiv:0911.5004]

### Holographic Reconstruction

10. **Ryu-Takayanagi 2006**: S. Ryu and T. Takayanagi, "Holographic Derivation of Entanglement Entropy from AdS/CFT," PRL 96, 181602. [arXiv:hep-th/0603001]

11. **Cao-Carroll-Michalakis 2017**: C. Cao, S. M. Carroll, and S. Michalakis, "Space from Hilbert Space: Recovering Geometry from Bulk Entanglement," PRD 95, 024031. [arXiv:1606.08444]

12. **Cao-Carroll 2018**: C. Cao and S. M. Carroll, "Bulk Entanglement Gravity without a Boundary: Towards Finding Einstein's Equation in Hilbert Space," PRD 97, 086003. [arXiv:1712.02803]

13. **Van Raamsdonk 2010**: M. Van Raamsdonk, "Building up Spacetime with Quantum Entanglement," Gen. Rel. Grav. 42, 2323. [arXiv:1005.3035]

### Exponential Metric Properties

14. **Papapetrou 1954**: A. Papapetrou, "Eine Theorie des Gravitationsfeldes mit einer Feldfunktion," Z. Phys. 139, 518.

15. **Yilmaz 1958**: H. Yilmaz, "New Approach to General Relativity," Phys. Rev. 111, 1417.

16. **Boonserm et al. 2018**: P. Boonserm, T. Ngampitipan, A. Simpson, and M. Visser, "Exponential metric represents a traversable wormhole," PRD 98, 084048. [arXiv:1805.03781]

17. **Dicke 1957**: R. H. Dicke, "Gravitation without a Principle of Equivalence," Rev. Mod. Phys. 29, 363.

### Information Theory and Gravity

18. **JRSWW 2018**: M. Junge, R. Renner, D. Sutter, M. M. Wilde, and A. Winter, "Universal Recovery Maps and Approximate Sufficiency of Quantum Relative Entropy," Ann. Henri Poincare 19, 2955.

19. **Petz 1988**: D. Petz, "Sufficiency of channels over von Neumann algebras," Q. J. Math. 39, 97.

20. **Witten 2022**: E. Witten, "Gravity and the Crossed Product," JHEP 2022(10), 008. [arXiv:2112.12828]

21. **Chandrasekaran et al. 2023**: V. Chandrasekaran, R. Longo, G. Penington, and E. Witten, "An Algebra of Observables for de Sitter Space," JHEP 2023(02), 082. [arXiv:2206.10780]

22. **Herrera 2020**: L. Herrera, "Landauer Principle and General Relativity," Entropy 22, 340.

### GfE Extensions

23. **Bianconi 2025b**: G. Bianconi, "The Quantum Relative Entropy of the Schwarzschild Black Hole and the Area Law," Entropy 27(3), 266.

24. **Bianconi 2026**: G. Bianconi, "The Thermodynamics of the Gravity from Entropy Theory." [arXiv:2510.22545]

### Spatial Metric from Information

25. **Czech et al. 2015**: B. Czech, L. Lamprou, S. McCandlish, and J. Sully, "Integral Geometry and Holography," JHEP 2015(10), 175. [arXiv:1505.05515]

26. **Bao et al. 2020**: N. Bao, C. Cao, S. Fischetti, and C. Keeler, "Building bulk geometry from the tensor Radon transform," JHEP 2020(12), 033.

---

## Appendix A: Summary Table of All Routes

| Route | Can derive g_rr? | Gives exponential? | Information-theoretic? | Status |
|-------|-----------------|-------------------|----------------------|--------|
| A. Isometry g_00*g_rr=-1 | YES (assumed) | YES | Partially motivated | ASSUMPTION |
| B. Spatial Petz | No viable formulation | N/A | Would be, if feasible | DEAD END |
| C. Jacobson program | YES (via field eqs) | NO (Schwarzschild) | YES | WRONG METRIC |
| D. Dorau-Much | YES (via field eqs) | NO (Schwarzschild) | YES | WRONG METRIC |
| E. Bianconi GfE | YES (via field eqs) | NO (Schw + r^{-4}) | YES | WRONG METRIC |
| F. Entropy tensor | Ill-defined | N/A | No formulation | NOT VIABLE |
| G. Holographic (RT) | In principle | Not tested | YES (in AdS) | PROMISING BUT INCOMPLETE |
| H. Cao-Carroll | In principle | Not tested | YES | PROMISING BUT INCOMPLETE |
| I. Modified Jacobson | Possible | Unknown | Would be | UNEXPLORED |

## Appendix B: The Dicke Argument for Isotropic Refraction

Dicke (1957) showed that GR in isotropic coordinates can be written as:

```
ds^2 = n(r)^{-2} (-dt^2 + dr^2 + r^2 dOmega^2)   [conformal to Minkowski]
```

Wait -- this is NOT correct for general GR metrics. The Schwarzschild metric in isotropic coordinates is:

```
ds^2 = -((1-m/r)/(1+m/r))^2 dt^2 + (1+m/r)^4 (dr^2 + r^2 dOmega^2)
```

This is NOT conformally flat (g_00 * g_rr != -1). So Dicke's refractive index formulation does NOT automatically give isotropic refraction for Schwarzschild.

However, the exponential metric IS conformally flat in isotropic coordinates:

```
ds^2 = -exp(-r_s/r) dt^2 + exp(+r_s/r) (dr^2 + r^2 dOmega^2)
     = exp(+r_s/r) [-exp(-2r_s/r) dt^2 + dr^2 + r^2 dOmega^2]
```

Wait, this is NOT conformally flat either, because exp(-2r_s/r) != 1. Let me reconsider.

A metric is conformally flat if it can be written as:

```
g_ab = Omega^2 eta_ab
```

For the exponential metric:
```
g_00 = -exp(-r_s/r)
g_rr = +exp(+r_s/r)
```

If conformally flat: g_00 = -Omega^2 and g_rr = +Omega^2, so g_00/g_rr = -1 identically. But for the exponential metric: g_00/g_rr = -exp(-2r_s/r) != -1 (unless r_s = 0).

Therefore the exponential metric is NOT conformally flat. The condition g_00 * g_rr = -1 is DIFFERENT from conformal flatness.

**Correction**: g_00 * g_rr = -exp(-r_s/r) * exp(+r_s/r) = -1 for the exponential metric. But this does NOT mean the metric is conformally flat. Conformal flatness requires g_00 = -Omega^2 and g_rr = Omega^2 for the SAME Omega, which would give g_rr = -g_00, i.e., A(r) = B(r). The exponential metric has A = exp(-r_s/r) and B = exp(+r_s/r), so A != B. It satisfies A*B = 1 but NOT A = B.

**The correct geometric interpretation**: g_00 * g_rr = -1 means the metric has UNIT determinant in the (t, r) sector:

```
det(g_{tr}) = g_00 * g_rr - g_{0r}^2 = -1  [for diagonal metrics with g_00*g_rr = -1]
```

This means the (t, r) part of the metric is an element of SL(2, R) -- the area element in the (t, r) plane equals the Minkowski area element. This is a special geometric property that is weaker than conformal flatness.

**Physical interpretation**: In the (t, r) plane, proper time intervals and proper radial distances are "reciprocally related" -- if clocks slow down by a factor alpha, then radial rulers stretch by a factor 1/alpha. This is the content of the isometry condition.

---

## Appendix C: Why the Problem Is Hard

The fundamental difficulty in deriving g_rr from information theory is rooted in the asymmetry between time and space in quantum mechanics:

1. **Time**: Has a natural direction (past to future). Quantum channels model evolution in time. The Petz recovery map recovers past states from future states. All of this is inherently temporal.

2. **Space**: Has no natural direction in the same sense. Spatial propagation at fixed time is described by the Hamiltonian, not by a channel. The distinction between "input" and "output" (necessary for defining a channel) is temporal, not spatial.

3. **The Petz recovery framework**: Is designed for quantum channels (CPTP maps), which model temporal processes. There is no natural "spatial Petz recovery" because there is no natural "spatial channel."

4. **Entropy production**: Is defined as the drop in QRE under a channel. Channels are temporal objects. Spatial entropy production would require a "spatial channel" whose physical meaning is unclear.

This asymmetry between time and space is fundamental to quantum mechanics and is reflected in the metric structure: g_00 (temporal) is naturally connected to information channels and recovery, while g_rr (spatial) is connected to geometry and distances but not to information degradation.

**The deep insight**: The fact that information theory naturally constrains g_00 but not g_rr reflects the time-space asymmetry in quantum mechanics. This is not a bug but a feature: it shows that the information-theoretic approach to gravity has a natural domain of applicability (temporal metric components) and a natural boundary (spatial components require additional input).

---

**Last updated**: 2026-03-10
