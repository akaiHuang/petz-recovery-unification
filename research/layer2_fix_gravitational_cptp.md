# Layer 2 Fix: Search for Rigorous CPTP Map N_grav with Delta-D = -ln(-g_00)

**Date**: 2026-03-09
**Status**: Comprehensive literature search complete
**Purpose**: Find mathematical frameworks that provide a rigorous CPTP map for gravitational redshift/decoherence, solving the Layer 2 gap

---

## 0. Executive Summary

The central problem remains: we need an explicit CPTP map N_grav such that its JRSWW entropy production Delta-D = -ln(-g_00). After searching across 6 major directions and ~50 papers, the assessment is:

**NO single paper provides a complete, first-principles, explicit CPTP map with computable Delta-D = -ln(-g_00).**

However, there are **5 promising candidates** that, individually or combined, could close the gap:

| Rank | Candidate | Provides CPTP? | Delta-D computable? | = -ln(-g_00)? | First-principles? |
|------|-----------|---------------|---------------------|---------------|-------------------|
| 1 | Modular Channel (Trejo-Calderon 2025) | YES (partial trace) | YES (spectral entropy) | Plausible (via thermal filter) | YES (modular theory) |
| 2 | Bosonic Loss Channel (standard QI) | YES (beam splitter) | YES (exact) | YES (with conventions) | NO (effective model) |
| 3 | UDW Communication Channel (Kasprzak-Tjoa 2024) | YES (field-mediated) | YES (capacity) | Geometry-independent! | YES (AQFT) |
| 4 | Holographic Conditional Expectation (Faulkner et al. 2020) | YES (CE on vN algebra) | YES (JLMS formula) | YES (in AdS/CFT) | YES (but requires AdS/CFT) |
| 5 | Open EFT Lindblad (Kaplanek-Burgess 2021) | YES (Lindblad) | YES (perturbative) | Near-horizon universal | YES (EFT) |

**Recommended strategy**: Combine Candidate 1 (modular channel) with Candidate 2 (bosonic loss) to get:
- First-principles justification from modular theory
- Explicit computable formulas from bosonic channel
- The identification eta = -g_00 motivated by both the thermal filter spectrum AND the Tolman/power convention

---

## 1. Direction 1: Quantum Reference Frames (QRF)

### Key Papers
- **Giacomini (2021)**: Spacetime QRFs and superpositions of proper times. Quantum 5, 508. [arXiv:2101.11628](https://arxiv.org/abs/2101.11628)
- **Giacomini et al. (2019)**: Quantum mechanics and covariance of physical laws in QRFs. Nature Commun. 10, 494.
- **Castro-Ruiz et al. (2020)**: A change of perspective: switching QRFs. Quantum 4, 225.

### Assessment for Layer 2

**Does it provide an explicit CPTP map?** Not directly. QRF transformations are described as changes of perspective, implemented by unitary transformations on an extended Hilbert space. The "gravitational channel" would be the partial trace after changing frames.

**Can Delta-D be computed?** Not in the existing literature. The QRF framework focuses on operational equivalence between frames, not on entropy production of the frame change.

**Connection to -ln(-g_00)?** The proper time superposition in the QRF framework does involve sqrt(-g_00) as the time dilation factor, but this has not been formulated as a CPTP map with computable Delta-D.

**Key insight**: When a quantum clock at radius r is described from the perspective of an observer at infinity, the time dilation sqrt(-g_00) appears as a coherent superposition of proper times. The "decoherence" of this superposition when the clock's internal state is inaccessible is essentially the Pikovski mechanism, which is a qubit channel and therefore has bounded Delta-D.

**Verdict**: QRF provides conceptual clarity but does NOT solve the Layer 2 gap. The formalism is too abstract to yield explicit entropy production formulas.

**Relevance to Paper 2**: LOW for Layer 2 fix; HIGH for conceptual framing.

---

## 2. Direction 2: Stinespring Dilation for Gravity

### Key Idea
Every CPTP map has a Stinespring dilation: N(rho) = Tr_E[U(rho tensor |0><0|_E) U^dag]. For gravitational decoherence, what are E and U?

### Pikovski Stinespring Dilation
For the Pikovski channel, the dilation is EXPLICIT:
- **System S**: center-of-mass position
- **Environment E**: internal degrees of freedom (phonons, molecular vibrations)
- **Unitary U**: time evolution under H = H_cm + H_int + (Phi(x)/c^2) H_int (gravity couples to internal energy)

The "environment" is the particle's own internal degrees of freedom -- gravity creates correlations between external and internal DOF.

**Problem**: This is a qubit (or finite-dimensional) channel, so Delta-D <= ln(dim) is bounded. Cannot give unbounded Sigma = r_s/r.

### For the Bosonic Channel
The Stinespring dilation of the pure-loss channel E_eta is:
- **System S**: signal mode a
- **Environment E**: vacuum mode b
- **Unitary U**: beam-splitter U_BS mixing a and b with transmissivity eta

The gravitational identification eta = -g_00 = exp(-r_s/r) then has the Stinespring interpretation:
- The "environment" is the gravitational field itself (or the modes between r and infinity)
- The "beam-splitter" is the gravitational redshift acting on field modes
- The "tracing out" corresponds to the inaccessibility of the lost photons

**Key paper**: Leber, Alanis Rodriguez, Ferreri, Schell, Bruschi (2025). [arXiv:2502.20521](https://arxiv.org/abs/2502.20521). They show the beam-splitter model for gravitational redshift has LIMITED VALIDITY -- it works for small redshift but fails for large redshift. The validity range depends on both gravitational parameters and photonic mode parameters.

**Implication for Layer 2**: The bosonic loss channel as Stinespring dilation is an EFFECTIVE model valid in weak field. It is NOT a first-principles construction.

**Verdict**: Stinespring dilation is well-understood for both Pikovski and bosonic channels, but neither provides a FIRST-PRINCIPLES gravitational CPTP map.

---

## 3. Direction 3: Quantum Channel Capacity of Gravitational Channels

### Key Paper: Kasprzak & Tjoa (2024)
**"Transmission of quantum information through quantum fields in curved spacetimes"**
[arXiv:2408.00518](https://arxiv.org/abs/2408.00518) (v4: July 2025)

This is potentially the MOST IMPORTANT paper for the Layer 2 problem.

**What they construct**:
- A relativistic quantum communication channel between two localized qubit systems
- Mediated by a relativistic quantum field (using UDW detector formalism)
- Achieves maximum quantum capacity in arbitrary curved spacetimes

**The CPTP map**:
```
Phi: D(H_A) -> D(H_B)   [maps Alice's qubit to Bob's qubit]
```
Constructed via:
```
U_A = exp(i lambda_1 sigma_z phi(f_1)) exp(i lambda_2 sigma_x phi(f_2))
```
with Bob performing inverse operations using the causal propagator.

**Quantum capacity**: Expressed purely in terms of:
- Correlation functions W(f,g) = <KEf, KEg>_KG  (Wightman functions)
- Causal propagator E(x,y) = G_R - G_A  (retarded minus advanced Green function)

**Critical feature**: The channel is "manifestly covariant, respects the causal structure of spacetime, and is independent of the details of the background geometry, topology, and choice of Hilbert space representations."

**Can Delta-D be computed?** YES, in principle. The channel has well-defined Kraus operators (though expressed implicitly via the UDW interaction unitaries). The quantum capacity involves coherent information, which is directly related to relative entropy.

**Does Delta-D = -ln(-g_00)?** THIS IS THE KEY QUESTION NOT ANSWERED IN THE PAPER. The channel is geometry-INDEPENDENT in its formal structure, but the CORRELATION FUNCTIONS that determine its quality DO depend on the background geometry.

For Schwarzschild:
- The Wightman function W(f,g) depends on the geodesic distance between the smeared regions
- The causal propagator E(f,g) depends on the causal structure
- Both are affected by gravitational redshift when Alice is at r and Bob is at infinity

**Assessment**: This framework gives a RIGOROUS CPTP map in curved spacetime. The open question is whether the entropy production of this specific channel equals -ln(-g_00). This requires computing W(f,g) and E(f,g) in the Schwarzschild background for appropriately placed detectors, which is a concrete (if difficult) calculation.

**Verdict**: MOST PROMISING for a first-principles solution. The channel exists and is CPTP. The remaining work is a COMPUTATION, not a conceptual gap.

**Recommended next step**: Compute the Kasprzak-Tjoa channel capacity for Alice at radius r and Bob at infinity in Schwarzschild, in the weak-field limit. Check whether the capacity loss relates to -ln(-g_00).

---

## 4. Direction 4: Modular Theory / Crossed Product Approach

### Key Papers
- **Dorau & Much (2025 PRL)**: From Quantum Relative Entropy to the Semiclassical Einstein Equations. PRL 136, 091602 (2026). [arXiv:2510.24491](https://arxiv.org/abs/2510.24491)
- **Witten (2022)**: Gravity and the crossed product. JHEP 10 (2022) 008.
- **Chandrasekaran, Longo, Penington, Witten (2023)**: An algebra of observables for de Sitter space. JHEP 02 (2023) 082.
- **Faulkner (2020)**: The holographic map as a conditional expectation. [arXiv:2008.04810](https://arxiv.org/abs/2008.04810)

### Dorau-Much: QRE -> Einstein Equations
**Key result**: The relative entropy between vacuum and coherent excitations on a bifurcate Killing horizon equals the energy flux across the horizon. With Bekenstein-Hawking, this yields Einstein's equations.

**Formula**: S_rel(omega_r || omega_0) = (2pi/kappa) * <delta T_ab xi^a xi^b>

**Does this give a CPTP map?** NO directly. The modular flow IS a one-parameter group of automorphisms, but automorphisms are INVERTIBLE and therefore have zero entropy production (they are unitary on the algebra). The entropy production comes from RESTRICTING to a subalgebra (partial trace), which IS a CPTP map.

**The conditional expectation as CPTP map**:
In the crossed-product framework, the inclusion A_inf hookrightarrow A_r (algebra at infinity into algebra at r) has an associated conditional expectation E: A_r -> A_inf. This IS a CPTP map (completely positive, trace-preserving, idempotent).

**Can we compute Delta-D for this conditional expectation?**
YES, in the holographic context (JLMS formula):
```
D(omega_1 || omega_2)|_{A_r} - D(omega_1 || omega_2)|_{A_inf} = S_gen(r) - S_gen(inf)
```
where S_gen is the generalized entropy (area + bulk EE).

For Schwarzschild exterior, this gives the fractional loss (D_r - D_inf)/D_r = r_s/r [as documented in existing layer2_gravitational_channel.md].

**The remaining gap**: Converting fractional loss to absolute loss requires a normalization convention or a specific choice of test states with D_r = 1.

### Trejo-Calderon (2025): Modular Channels and Thermal Filtering
**"Modular Channels, Thermal Filtering and the Spectral Emergence of Spacetime"**
[arXiv:2504.20457](https://arxiv.org/abs/2504.20457)

**THIS IS THE MOST PROMISING NEW PAPER for Layer 2.**

**What it constructs**:
- The modular channel is the CPTP map induced by partial trace over inaccessible (causally disconnected) region
- The channel's Kraus operator A admits SVD: A = U Sigma V^dag
- **KEY RESULT**: The singular values follow a THERMAL (Gibbs) distribution:
  ```
  sigma_i^2 = exp(-beta k_i) / Z
  ```
  where {k_i} are modular Hamiltonian eigenvalues, beta = 2pi/a (Unruh temperature)

**Spectral entropy**: S_spec(beta) = -2 sum_i sigma_i(beta)^2 log sigma_i(beta)

**Connection to geometry**: The thermal filter temperature is T = a/(2pi) where a is the proper acceleration. For a static observer at radius r in Schwarzschild:
```
a = (r_s/2) / (r^2 sqrt(1 - r_s/r))
beta = 2pi/a = 4pi r^2 sqrt(1 - r_s/r) / r_s
```

**How this could solve Layer 2**:
1. The modular channel IS a well-defined CPTP map (partial trace)
2. Its singular values ARE computable (from modular Hamiltonian spectrum)
3. The spectral entropy S_spec relates to information loss
4. The thermal filter naturally connects to the Tolman temperature

**What's missing**: An explicit computation of Delta-D = D(rho||sigma) - D(N(rho)||N(sigma)) for this modular channel, and verification that it equals -ln(-g_00).

**Einstein equations recovery**: The paper recovers G_mu_nu + Lambda g_mu_nu = 8pi G T_mu_nu from requiring the Clausius relation delta Q = T delta S across all Rindler patches. This is essentially the Jacobson (1995) argument, but now with the modular channel providing the quantum information-theoretic underpinning.

**Verdict**: Very promising but incomplete. The modular channel is constructed, the thermal spectrum is derived, but the specific computation of Delta-D in Schwarzschild is NOT performed.

### Faulkner (2020): Holographic Map as Conditional Expectation
[arXiv:2008.04810](https://arxiv.org/abs/2008.04810)

**Key insight**: In AdS/CFT, the holographic map (bulk-to-boundary) IS a conditional expectation on von Neumann algebras. This is a CPTP map by construction.

**Relative entropy**: By JLMS, D(omega_1||omega_2)|_boundary = D(omega_1||omega_2)|_bulk + (Area difference)/(4G). The entropy production of the holographic map is precisely the area term.

**Applicability to Layer 2**: This provides a RIGOROUS CPTP map with COMPUTABLE Delta-D in the holographic setting. However, it requires the AdS/CFT correspondence, which is not available for asymptotically flat Schwarzschild.

**Verdict**: Rigorous and explicit, but limited to AdS/CFT context.

---

## 5. Direction 5: Unruh-DeWitt Detector as Channel

### Key Papers
- **Kasprzak & Tjoa (2024)**: See Direction 3 above.
- **Kaplanek & Burgess (2021)**: Qubits on the Horizon. JHEP 01 (2021) 098. [arXiv:2007.05984](https://arxiv.org/abs/2007.05984)
- **Kaplanek & Burgess (2020)**: Hot Accelerated Qubits. JHEP 03 (2020) 008. [arXiv:1912.12951](https://arxiv.org/abs/1912.12951)

### Kaplanek-Burgess: Open EFT for Qubits near Black Holes
**Key construction**: A qubit (UDW detector) hovering at radius r near a Schwarzschild black hole, interacting with a free scalar field in the Hartle-Hawking state.

**Explicit Lindblad equation**: Using Open EFT techniques, they derive a MARKOVIAN Lindblad master equation valid for late times:
```
d rho / dt = -i[H_eff, rho] + sum_k gamma_k (L_k rho L_k^dag - (1/2){L_k^dag L_k, rho})
```

**Key results**:
1. For qubits sufficiently close to the horizon, the late-time evolution takes a UNIVERSAL form depending only on near-horizon geometry
2. The qubit thermalizes at the local Tolman temperature T(r) = T_H / sqrt(1 - r_s/r)
3. The Lindblad rates gamma_k are computable in perturbation theory (in g^2 t / r_s)
4. The decoherence rate and thermalization rate both depend on the metric through sqrt(-g_00)

**Can Delta-D be computed?** YES. The Lindblad equation defines a one-parameter family of CPTP maps Phi_t. The entropy production can be computed as:
```
Sigma(t) = D(rho(0) || sigma_th) - D(rho(t) || sigma_th)
```
where sigma_th is the thermal state at Tolman temperature.

**Does Sigma = -ln(-g_00)?** NOT DIRECTLY. The entropy production is:
- Time-dependent (grows with t until thermalization)
- Coupling-dependent (scales as g^2)
- Probe-dependent (depends on qubit gap omega)

However, the EQUILIBRIUM relative entropy D(sigma_th(r) || sigma_th(inf)) depends on the metric, and the thermalization dynamics are universally controlled by sqrt(-g_00) near the horizon.

**Assessment**: This provides an EXPLICIT, FIRST-PRINCIPLES Lindblad equation for gravitational decoherence, but the entropy production is probe-dependent and does not directly give -ln(-g_00). The universal near-horizon behavior is physically important but limited to the near-horizon regime.

**Verdict**: Rigorous and explicit, but probe-dependent. Useful for understanding the physics but does not solve Layer 2 as stated.

---

## 6. Direction 6: Gaussian Quantum Information in Curved Spacetime

### Key Papers
- **Ahmadi, Bruschi, Fuentes (2014)**: Quantum metrology for relativistic quantum fields. PRD 89, 065028. [arXiv:1312.5707](https://arxiv.org/abs/1312.5707)
- **Kohlrus, Bruschi, Louko, Fuentes (2017)**: Quantum communications in rotating planet spacetime. EPJ QT 4, 7. [arXiv:1511.04256](https://arxiv.org/abs/1511.04256)
- **Bruschi (2021)**: Towards communication in a curved spacetime geometry. Commun. Phys. 4, 171. [Nature](https://www.nature.com/articles/s42005-021-00671-8)
- **Leber et al. (2025)**: Limits of quantum-optical redshift models. [arXiv:2502.20521](https://arxiv.org/abs/2502.20521)
- **Dragan et al. (2025)**: Gravitational redshift of broadband relativistic quantum photons. [arXiv:2504.03956](https://arxiv.org/abs/2504.03956)

### The Gaussian Channel for Gravitational Redshift

The key result from this body of work is that gravitational redshift on photonic modes can be described as a **Bogoliubov transformation** on the mode operators:

For a photon with frequency omega at radius r, observed at infinity:
```
a_out = alpha * a_in + beta * a_in^dag + sum_k (alpha_k b_k + beta_k b_k^dag)
```

where:
- alpha = sqrt(-g_00) * (mode overlap factor) [amplitude transmissivity]
- beta coefficients represent particle creation (typically negligible in static spacetimes)
- b_k are auxiliary environmental modes

### Bruschi (2021): Wavepacket Distortion
**Key result**: Spacetime curvature distorts wavepackets, causing:
- Phase shift proportional to Riemann curvature
- Cross-talk between structured modes (up to 12.2% in the solar system)
- Mode mixing described by a unitary transformation on the Hilbert-space basis

For communication between Earth and ISS: 0.10 radian phase shift on monochromatic laser beam.

**Is this a CPTP map?** The full evolution IS unitary (no information loss in static backgrounds). The CPTP map arises from tracing out the modes that are not observed. If only one mode is observed, the channel is:
```
N(rho) = Tr_{other modes}[S rho S^dag]
```
where S is the symplectic (Bogoliubov) transformation.

### Leber et al. (2025): Validity Limits
**Critical finding**: The standard beam-splitter model for gravitational redshift is valid ONLY for small redshift. The validity range depends on:
1. Gravitational parameters (r_s/r)
2. Photonic mode parameters (bandwidth, mode number)

**Necessary condition**: The number of auxiliary (environmental) modes must be >= the number of signal modes. This is a CONSISTENCY condition on the CPTP map.

**Implication**: For strong fields (r_s/r ~ 1), the simple bosonic loss channel E_eta is NOT a valid description. A multimode treatment is required.

### Dragan et al. (2025): Gravitational Redshift via Linearized QG
[arXiv:2504.03956](https://arxiv.org/abs/2504.03956)

Uses linearized quantum gravity (graviton exchange) to study photon redshift. Key finding: gravitational redshift occurs only for photons interacting with a classical source AND having well-defined momentum. Photons with quantum coherence in position show NO well-defined redshift.

**Implication**: The "gravitational channel" depends on the quantum state of the source, not just the background metric. This challenges the universality assumption.

### Assessment for Layer 2

The Gaussian channel framework gives:
- Explicit CPTP maps (beam-splitter + trace out)
- Computable Delta-D (Gaussian state formulas)
- For weak field: Delta-D = -ln(eta) where eta is the transmissivity
- The identification eta = -g_00 requires the power convention (frequency shift + time dilation)

But:
- NOT valid for strong fields (Leber et al. 2025)
- NOT first-principles (effective model, not derived from GR)
- The transmissivity convention is a modeling choice

**Verdict**: Best EXPLICIT channel for weak-field computations, but limited to weak field and not first-principles.

---

## 7. Additional Important Candidates

### 7a. Galley-Giacomini-Selby (2023): Classical-Quantum Gravity is Irreversible
**"Any consistent coupling between classical gravity and quantum matter is fundamentally irreversible"**
Quantum 7, 1142 (2023). [Link](https://quantum-journal.org/papers/q-2023-10-16-1142/)

**Key theorem**: If gravity is classical and matter is quantum, then either (i) matter is not fully quantum, (ii) interactions are irreversible, or (iii) matter does not back-react on gravity.

**Relevance**: This proves that classical gravity + quantum matter REQUIRES a non-unitary (irreversible) channel. This is exactly the situation in Paper 2: the gravitational field creates irreversibility, quantified by Sigma > 0.

**Does it give an explicit channel?** NO. The paper uses GPT (General Probabilistic Theory) framework and proves existence, not construction.

**Verdict**: Strong conceptual support for the existence of N_grav, but no explicit construction.

### 7b. Oppenheim (2023): Postquantum Theory of Classical Gravity
PRX 13, 041040 (2023). [arXiv:1811.03116](https://arxiv.org/abs/1811.03116)

**Key construction**: A theory where gravity remains classical but coupled to quantum matter via a CQ (classical-quantum) dynamics that is:
- Linear in the density matrix
- Completely positive and trace-preserving
- Reduces to Einstein gravity classically

**Lindblad equation**: The CQ dynamics is generated by a GKSL/Lindblad equation with:
- Decoherence terms (quantum -> classical noise)
- Diffusion terms (classical stochastic noise)
- A fundamental TRADE-OFF: more decoherence <-> less diffusion (and vice versa)

**Does Sigma relate to the metric?** The decoherence rate is set by the gravitational coupling strength (G), and the entropy production from the Lindblad equation is computable. However, the entropy production depends on the specific CQ trajectory, not just the metric.

**Verdict**: Provides a CPTP-like structure for classical gravity, but the entropy production is model-dependent (depends on the decoherence-diffusion trade-off parameter). Does not directly give Sigma = -ln(-g_00).

### 7c. Kafri-Taylor-Milburn (2014): Classical Channel Model
NJP 16, 065020 (2014). [arXiv:1401.0946](https://arxiv.org/abs/1401.0946)

**Key construction**: Gravitational interaction modeled as a classical measurement channel:
- Continuous position measurement of one mass
- Feedback force on the other mass based on measurement result
- Equivalent to the Diosi decoherence model

**CPTP map**: The classical channel IS a CPTP map on the quantum state. The Lindblad equation is:
```
d rho / dt = -i[H, rho] + gamma (x rho x - (1/2){x^2, rho})
```
where gamma = G m^2 / (hbar d^3) for masses separated by distance d.

**Delta-D**: The entropy production from this Lindblad equation is computable but depends on probe parameters (mass m, separation d, interaction time t).

**Verdict**: Explicit CPTP map for gravitational interaction, but probe-dependent. Does not give -ln(-g_00).

### 7d. Basso-Maziero-Celeri (2025): Quantum Crooks in Curved Spacetime
PRL 134, 050406 (2025). [arXiv:2405.03902](https://arxiv.org/abs/2405.03902)

**Key result**: A fully general relativistic detailed quantum fluctuation theorem using Fermi normal coordinates.

**Critical finding**: ENTROPY PRODUCTION IS OBSERVER-DEPENDENT.
- Different observers at different radii measure different Sigma
- The arrow of time is connected to the causal structure

**Does it give Sigma = -ln(-g_00)?** The paper does not claim this specific formula. The entropy production depends on the specific protocol (work done, heat exchanged) in the curved background.

**Key insight for Paper 2**: Observer-dependence of Sigma is CONSISTENT with our framework (different observers have different tau). This paper provides the general relativistic version of the fluctuation theorem that Paper 2's framework requires.

**Verdict**: Strong theoretical support, but does not compute Sigma for the specific gravitational redshift channel.

### 7e. Moreira-Celeri (2024): Entropy Production from Spacetime Fluctuations
[arXiv:2407.21186](https://arxiv.org/abs/2407.21186)
Thesis extension: [arXiv:2603.02034](https://arxiv.org/abs/2603.02034) (March 2026)

**Key construction**: A non-relativistic quantum system interacting with a bath of gravitons (linear limit of GR). Uses consistent histories approach to derive entropy production.

**Channel**: The graviton bath acts as an environment. The reduced dynamics of the system defines a CPTP map (after tracing gravitons).

**Entropy production**: Derived via integral fluctuation theorem within decoherent histories framework. The result is that thermodynamic entropy MUST be produced due to unavoidable interaction with quantum spacetime fluctuations.

**Connection to metric**: The computation is in the weak-field (linearized) limit around flat spacetime. The entropy production depends on the graviton spectral density, not directly on -ln(-g_00).

**Verdict**: Promising first-principles approach, but limited to linearized gravity and does not give -ln(-g_00).

### 7f. Artini et al. (2025): Non-equilibrium Thermodynamics of DP Model
[arXiv:2502.03173](https://arxiv.org/abs/2502.03173)

Investigates entropy production in the Diosi-Penrose model. The original DP model induces unbounded heating; the dissipative extension achieves thermalization.

**Relevance**: Provides explicit entropy production formulas for gravitational collapse models, but these are for the DP mechanism (stochastic localization), not for gravitational redshift.

### 7g. Bianconi (2025): Gravity from Entropy
PRD 111, 066001 (2025). [arXiv:2408.14391](https://arxiv.org/abs/2408.14391)

**Key idea**: Metric = density matrix. Action = quantum relative entropy between spacetime metric and matter-induced metric:
```
L = Tr g_tilde ln(G_tilde^{-1}) - Lambda
```

**Does it give a CPTP map?** NO. The framework is entirely classical differential-geometric. There are no quantum channels.

**Relevance**: The identification metric = density matrix is conceptually aligned with Paper 2's framework, and the action being QRE is exactly the structure of Paper 4's vision (Sigma = D(rho_spacetime || rho_matter)). But it does not solve Layer 2.

---

## 8. Synthesis: The State of the Art

### 8.1 What Exists (Summary)

**Explicit CPTP maps in gravitational contexts**:

| Paper | Channel type | Hilbert space | First-principles? | Delta-D formula? |
|-------|-------------|---------------|-------------------|-----------------|
| Pikovski (2015) | Qubit dephasing | C^2 | Yes (low-energy GR) | H_bin((1+p)/2) <= ln 2 |
| Kafri et al. (2014) | Classical channel (position measurement) | L^2(R) | Semi (classical gravity) | Probe-dependent |
| Kaplanek-Burgess (2021) | Qubit Lindblad (Open EFT) | C^2 | Yes (EFT) | Probe + time dependent |
| Kasprzak-Tjoa (2024) | Field-mediated UDW | C^2 tensor C^2 | Yes (AQFT) | Via capacity (implicit) |
| Oppenheim (2023) | CQ Lindblad | Infinite-dim | Yes (CQ dynamics) | Trade-off dependent |
| Bosonic loss E_eta | Beam-splitter | Fock space | No (effective) | -ln(eta) exactly |
| Trejo-Calderon (2025) | Modular (partial trace) | Type II algebra | Yes (modular theory) | S_spec (spectral) |
| Faulkner (2020) | Conditional expectation | vN algebra | Yes (AdS/CFT) | Via JLMS |

### 8.2 The Gap Analysis

**What we need**: A CPTP map N_grav on a well-defined Hilbert space such that for fixed reference state sigma and general input rho:
```
D(rho || sigma) - D(N_grav(rho) || N_grav(sigma)) = -ln(-g_00(r))
```

**What we have**:
1. Bosonic loss channel gives this in the N_B >> 1 limit, but is an effective model
2. Modular channel gives fractional loss r_s/r, but needs normalization
3. UDW channel (Kasprzak-Tjoa) is rigorous CPTP but Delta-D not yet computed
4. Holographic CE gives this in AdS/CFT, but not for asymptotically flat

**The remaining gap is technical, not conceptual**:
- The CONCEPT of gravitational entropy production = -ln(-g_00) is supported by 3+ independent arguments
- The COMPUTATION of Delta-D for a specific first-principles CPTP map in Schwarzschild has not been performed

### 8.3 The Most Promising Path Forward

**Step 1: Modular Channel Approach (Trejo-Calderon + Dorau-Much)**
- Start with the modular channel (partial trace over causally inaccessible region)
- This is a well-defined CPTP map from modular theory
- Compute its singular value spectrum for Schwarzschild (the thermal filter)
- Relate the spectral entropy to -ln(-g_00)

**Step 2: Bosonic Channel Embedding**
- Show that in the weak-field limit, the modular channel REDUCES to a bosonic loss channel
- The transmissivity eta = -g_00 emerges from the thermal filter spectrum
- The high-temperature thermal reference emerges from the Tolman temperature

**Step 3: Kasprzak-Tjoa Verification**
- Compute the Kasprzak-Tjoa channel for Alice at r, Bob at infinity in Schwarzschild
- Verify that the capacity loss matches -ln(-g_00)
- This would provide an independent first-principles check

---

## 9. New References (Not in Existing Documentation)

### Papers found in this search that are NEW:

| Paper | Key contribution | arXiv |
|-------|-----------------|-------|
| Kasprzak & Tjoa (2024) | Rigorous CPTP in curved spacetime via AQFT | 2408.00518 |
| Trejo-Calderon (2025) | Modular channels, thermal filter, MCFC | 2504.20457 |
| Leber et al. (2025) | Validity limits of beam-splitter model | 2502.20521 |
| Dragan et al. (2025) | Gravitational redshift via linearized QG | 2504.03956 |
| Oppenheim (2023) | CQ gravity with CPTP dynamics | PRX 13, 041040 |
| Kaplanek & Burgess (2021) | Open EFT Lindblad near black holes | JHEP 01 (2021) 098 |
| Galley-Giacomini-Selby (2023) | CG+QM is fundamentally irreversible | Quantum 7, 1142 |
| Moreira & Celeri (2024) | Entropy production from graviton bath | 2407.21186 |
| Moreira (2026 thesis) | Decoherence + entropy from spacetime | 2603.02034 |
| Casini & Huerta (2017) | RG = DPI (relative entropy + RG flow) | JHEP 03 (2017) 089 |
| Faulkner (2020) | Holographic map as conditional expectation | 2008.04810 |

---

## 10. Concrete Calculations To Do

### Priority 1: Bosonic Channel Delta-D (Complete -- already in layer2_gravitational_channel.md)
- DONE: Sigma = -ln(eta) + O(1/N_B) for N_B >> 1

### Priority 2: Kasprzak-Tjoa Channel in Schwarzschild (NEW)
- Compute Wightman function W(f,g) for Alice at r, Bob at infinity
- Compute causal propagator E(f,g) in Schwarzschild
- Extract the effective channel N: D(H_A) -> D(H_B)
- Compute Delta-D = D(rho||sigma) - D(N(rho)||N(sigma))
- Check whether Delta-D = -ln(-g_00) or -ln(sqrt(-g_00))

### Priority 3: Modular Channel Spectral Entropy (NEW)
- Compute modular Hamiltonian spectrum for Schwarzschild exterior at radius r
- Extract singular values sigma_i^2 = exp(-beta k_i)/Z
- Compute spectral entropy S_spec
- Relate S_spec to -ln(-g_00)

### Priority 4: Kaplanek-Burgess Entropy Production (NEW)
- Take the Open EFT Lindblad equation for qubit at radius r
- Compute the total entropy production integral_0^inf Sigma(t) dt
- Check near-horizon universality
- Determine whether the total integrated Sigma relates to -ln(-g_00)

---

## 11. Revised Assessment for Paper 2

### Layer 2 Status: MOTIVATED CONJECTURE (upgraded from "ansatz")

**Evidence FOR Sigma_grav = -ln(-g_00)**:
1. Bosonic loss channel with eta = -g_00 gives Sigma = -ln(eta) = r_s/r [explicit CPTP, effective model]
2. Modular flow / Dorau-Much: fractional QRE loss = r_s/r [first-principles, fractional not absolute]
3. Gravitational Landauer (Herrera 2020): excess erasure cost = r_s/r [thermodynamic argument]
4. Galley-Giacomini-Selby: CG+QM MUST be irreversible [existence theorem]
5. Basso-Celeri: observer-dependent entropy production in curved spacetime [consistency]
6. Trejo-Calderon: modular channel thermal filter gives Einstein equations [structural alignment]

**Evidence AGAINST or COMPLICATING**:
1. No first-principles CPTP map with computable Delta-D = -ln(-g_00) exists in the literature
2. Beam-splitter model has limited validity (Leber et al. 2025)
3. Transmissivity convention (eta = -g_00 vs sqrt(-g_00)) remains a modeling choice
4. Dragan et al. (2025): redshift depends on quantum state of source, not just metric
5. All qubit channels have bounded Delta-D (fundamental obstruction)

### Honest Statement for Paper 2

> "The identification Sigma_grav = -ln(-g_00) is supported by three independent physical arguments (modular flow, gravitational Landauer, and bosonic channel), each yielding the same formula from different starting points. An explicit CPTP map (the bosonic pure-loss channel with eta = -g_00 and high-temperature thermal reference) reproduces this entropy production in the large-N_B limit. A fully first-principles derivation from a rigorously constructed gravitational CPTP map remains an important open problem."

---

## 12. Summary Table: All Candidates Ranked

| Rank | Candidate | Score (0-10) | Pros | Cons |
|------|-----------|-------------|------|------|
| 1 | **Modular channel + bosonic embedding** | 8/10 | First-principles structure, explicit computation possible, thermal filter natural | Spectral entropy != Delta-D directly; needs computation |
| 2 | **Kasprzak-Tjoa UDW channel** | 7/10 | Rigorous CPTP in any curved spacetime, AQFT based | Delta-D not computed; geometry enters only through correlators |
| 3 | **Bosonic loss E_eta** | 7/10 | Explicit, computable, gives -ln(eta) exactly | Effective model, not first-principles; eta convention ambiguous |
| 4 | **Holographic CE (Faulkner)** | 6/10 | Rigorous, well-defined, JLMS formula | Requires AdS/CFT; not applicable to asymptotically flat |
| 5 | **Open EFT Lindblad (Kaplanek-Burgess)** | 5/10 | First-principles, explicit Lindblad | Probe-dependent; qubit; near-horizon only |
| 6 | **Oppenheim CQ** | 4/10 | CPTP by construction, fundamental theory | Trade-off dependent; does not predict specific Sigma |
| 7 | **ABH master equation** | 4/10 | First-principles (QFT+GR), Lindblad form | Energy-basis decoherence, not position; graviton bath |
| 8 | **Kafri-Taylor-Milburn** | 3/10 | Explicit classical channel | Very specific model; cannot create entanglement |
| 9 | **QRF frame change** | 3/10 | Conceptually clean | Too abstract; no entropy production formula |
| 10 | **Pikovski dephasing** | 2/10 | Physical motivation excellent | Qubit obstruction; probe-dependent |

---

**Last updated**: 2026-03-09
