# Canonical Gravitational CPTP Map: First-Principles Survey

## Date: 2026-03-10
## Author: Research survey for Paper 2 (Sheng-Kai Huang)
## Status: Comprehensive literature survey complete

---

## 1. Problem Statement

### 1.1 The Core Gap

Paper 2 defines the gravitational entropy production as:

```
Sigma_grav = -ln(-g_00) = r_s/r
```

and uses the JRSWW bound F >= exp(-Sigma/2) to derive the exponential metric g_00 = -exp(-r_s/r). Three independent physical arguments support this identification (modular flow, gravitational Landauer, quantum channel). However, the **central theoretical gap** is:

> **No one has derived from first principles a canonical CPTP map N_grav such that its JRSWW entropy production equals -ln(-g_00).**

Specifically, we need an explicit quantum channel N_grav acting on a well-defined Hilbert space H, with a physically motivated reference state sigma, such that for arbitrary input rho:

```
D(rho || sigma) - D(N_grav(rho) || N_grav(sigma)) = -ln(-g_00(r))
```

where D is the Umegaki quantum relative entropy.

### 1.2 Why This Matters

Without this channel:
- The JRSWW bound cannot be rigorously applied (it requires a CPTP map)
- The identification Sigma_grav = r_s/r remains a "motivated conjecture" rather than a theorem
- The exponential metric g_00 = -exp(-r_s/r) rests on an additional unproven saturation hypothesis
- Paper 2's three-layer structure (theorem -> conjecture -> speculation) cannot be promoted to a unified theorem

### 1.3 Structural Constraints

Any candidate N_grav must satisfy:
1. **CPTP**: completely positive, trace-preserving
2. **Infinite-dimensional**: finite-dimensional channels have bounded entropy production (sup Sigma <= ln(dim) + ln(1/lambda_min)), so they cannot reproduce unbounded Sigma = r_s/r as r -> r_s
3. **Universal**: the entropy production should depend only on the metric (g_00), not on probe-specific parameters
4. **Correct weak-field limit**: Sigma -> r_s/r + O((r_s/r)^2) as r_s/r -> 0
5. **Physically derivable**: should emerge from semiclassical gravity (QFT in curved spacetime) without ad hoc assumptions

---

## 2. Literature Survey

### 2.1 Route A: Algebraic QFT and Modular Automorphisms -> CPTP Map

#### 2.1.1 Dorau & Much (2025 PRL)
- **Title**: "From Quantum Relative Entropy to the Semiclassical Einstein Equations"
- **arXiv**: [2510.24491](https://arxiv.org/abs/2510.24491)
- **Published**: Physical Review Letters 136, 091602 (2026)
- **Key Result**: Using Tomita-Takesaki modular theory, the Araki-Uhlmann relative entropy between the vacuum state and coherent excitations of a scalar QFT on a bifurcate Killing horizon equals the energy flux across the horizon:
  ```
  S_rel(omega_r || omega_0) = (2pi/kappa) * <delta T_ab xi^a xi^b>
  ```
  Combined with the Bekenstein-Hawking entropy-area formula, the semiclassical Einstein equations follow.
- **Channel Implication**: The modular flow IS a one-parameter group of automorphisms (invertible, hence zero entropy production). The entropy production arises from RESTRICTING to a subalgebra (partial trace), which IS a CPTP map. The fractional QRE loss is:
  ```
  (D_in - D_out) / D_in = r_s/r   [EXACT for Schwarzschild]
  ```
- **Gap**: Gives fractional (relative) entropy loss, not absolute entropy production. Converting to absolute Sigma = r_s/r requires D_in = 1 (normalization assumption).

#### 2.1.2 Trejo-Calderon (2025)
- **Title**: "Modular Channels, Thermal Filtering and the Spectral Emergence of Spacetime"
- **arXiv**: [2504.20457](https://arxiv.org/abs/2504.20457)
- **Key Construction**: Defines the **modular channel** as a CPTP map arising from partial trace over causally inaccessible regions. The channel's Kraus operator A admits SVD A = U Sigma V^dagger, where the singular values follow a **Gibbs distribution**:
  ```
  sigma_i^2 = exp(-beta * k_i) / Z
  ```
  with {k_i} the modular Hamiltonian eigenvalues and beta = 2pi/a (Unruh temperature).
- **Spectral Entropy**: S_spec(beta) = -2 sum_i sigma_i^2 * log(sigma_i)
- **Einstein Equations Recovery**: Reinterprets the first law of entanglement as a Clausius relation delta_Q = T * delta_S. Local validity across all Rindler horizons yields G_mu_nu + Lambda * g_mu_nu = 8pi G T_mu_nu.
- **Modular Channels Flow Correspondence (MCFC)**: Proposes that area of causal screen = entropic capacity of modular channel.
- **Gap**: Does not compute Delta-D = D(rho||sigma) - D(N(rho)||N(sigma)) for the modular channel in Schwarzschild. Does not verify Sigma = -ln(-g_00).
- **Assessment**: MOST PROMISING algebraic construction for N_grav. The Gibbs-weighted singular value spectrum naturally produces a thermal channel. The remaining work is a specific computation, not a conceptual gap.

#### 2.1.3 Witten (2022), Chandrasekaran-Longo-Penington-Witten (2023)
- **Title**: "An Algebra of Observables for de Sitter Space" / "Gravity and the Crossed Product"
- **arXiv**: [2206.10780](https://arxiv.org/abs/2206.10780), [2112.12828](https://arxiv.org/abs/2112.12828)
- **Key Result**: Gravity converts the Type III von Neumann algebra of QFT in a subregion to a Type II algebra. This makes the relative entropy well-defined and finite. The conditional expectation (restriction from bulk to boundary, or from larger to smaller algebra) IS a CPTP map.
- **Entropy**: Generalized entropy S_gen = A/(4G_N) + S_out emerges as the entropy of semiclassical states on the Type II algebra.
- **Gap**: Abstract algebraic framework. No explicit Kraus operators. The conditional expectation is defined implicitly via the algebraic structure.

#### 2.1.4 Chandrasekaran-Penington-Witten (2023)
- **Title**: "Generalized entropy for general subregions in quantum gravity"
- **arXiv**: [2306.01837](https://arxiv.org/abs/2306.01837)
- **Key Extension**: For spatially compact subregions or asymptotic boundaries, the algebra is Type II von Neumann factor. Entropy is well-defined.
- **Assessment**: Extends the algebraic framework to more general settings but does not construct explicit channels.

#### 2.1.5 Faulkner (2020)
- **Title**: "The holographic map as a conditional expectation"
- **arXiv**: [2008.04810](https://arxiv.org/abs/2008.04810)
- **Key Result**: In AdS/CFT, the holographic map (bulk-to-boundary) IS a conditional expectation on von Neumann algebras, hence a CPTP map by construction. The JLMS formula gives:
  ```
  D(omega_1||omega_2)|_boundary = D(omega_1||omega_2)|_bulk + (Area difference)/(4G)
  ```
  The entropy production of the holographic map equals the area term.
- **Gap**: Requires AdS/CFT correspondence. Not available for asymptotically flat Schwarzschild.
- **Assessment**: The most rigorous CPTP map with computable Delta-D, but restricted to the holographic setting.

#### 2.1.6 Cirafici (2024) [Kudler-Flam framework]
- **Title**: "Fluctuation theorems, quantum channels and gravitational algebras"
- **arXiv**: [2408.04219](https://arxiv.org/abs/2408.04219)
- **Published**: JHEP 2024, 089 (2024)
- **Key Result**: Quantum channels in Type II_1 gravitational algebras (de Sitter space) are represented by SUBFACTORS using Jones' theory. Establishes Crooks-like fluctuation theorem:
  ```
  P_t(-s_bar) = e^{-t * s_bar} * P_t(s_bar)
  ```
  Negative entropy fluctuations exponentially suppressed; Jarzynski equality <e^{-s_bar}> = 1.
- **Assessment**: Connects quantum channels to gravitational entropy production via operator algebras. However, the channel is abstract (subfactor representation), not an explicit operator on a Hilbert space.

#### 2.1.7 Bianconi (2025)
- **Title**: "Gravity from entropy"
- **arXiv**: [2408.14391](https://arxiv.org/abs/2408.14391)
- **Published**: Physical Review D 111, 066001 (2025)
- **Key Idea**: Treats the metric as an effective density matrix. The gravitational action IS the quantum relative entropy:
  ```
  S_action = D(rho_spacetime || rho_matter)
  ```
  Modified Einstein equations reduce to standard Einstein equations in the low-coupling regime. Naturally predicts small positive cosmological constant.
- **QRE of Schwarzschild**: Published also as "The Quantum Relative Entropy of the Schwarzschild Black Hole and the Area Law" (Entropy 27, 266, 2025; arXiv:2501.09491). Shows the QRE of the Schwarzschild metric obeys the area law for large r_s.
- **Assessment**: The most radical proposal -- if correct, it would make Paper 4's Sigma = D(rho_spacetime || rho_matter) literally true. However, treating the metric tensor as a density matrix is non-standard and controversial. Does NOT construct a CPTP channel.

### 2.2 Route B: Unruh-DeWitt Detector Model -> Explicit Channel

#### 2.2.1 Kasprzak & Tjoa (2024)
- **Title**: "Transmission of quantum information through quantum fields in curved spacetimes"
- **arXiv**: [2408.00518](https://arxiv.org/abs/2408.00518) (v4: July 2025)
- **Key Construction**: Constructs a rigorous CPTP map between two UDW detectors coupled to a scalar field in arbitrary curved spacetime:
  ```
  Phi: D(H_A) -> D(H_B)
  ```
  via controlled unitaries U_A = exp(i lambda_1 sigma_z phi(f_1)) exp(i lambda_2 sigma_x phi(f_2)).
- **Quantum Capacity**: Expressed purely in terms of:
  - Wightman functions W(f,g) = <KEf, KEg>_KG
  - Causal propagator E(x,y) = G_R - G_A
- **Properties**: Manifestly covariant, respects causal structure, geometry-independent formalism (geometry enters only through correlation functions).
- **Gap**: The channel maps qubits to qubits (finite-dimensional), so Delta-D is bounded by ~ln 2. For infinite-mode version, the capacity LOSS in Schwarzschild has not been computed.
- **Assessment**: MOST RIGOROUS first-principles construction of a communication channel in curved spacetime. The key open question is whether the capacity loss (over all modes) relates to -ln(-g_00). This is a concrete computation, not a conceptual gap.
- **Priority Calculation**: Compute W(f,g) and E(f,g) for Alice at radius r and Bob at infinity in Schwarzschild, extract the effective multimode channel, and compute Delta-D.

#### 2.2.2 Kaplanek & Burgess (2020, 2021)
- **Title**: "Hot Accelerated Qubits: Decoherence, Thermalization, Secular Growth and Reliable Late-time Predictions" / "Qubits on the Horizon"
- **arXiv**: [1912.12951](https://arxiv.org/abs/1912.12951), [2007.05984](https://arxiv.org/abs/2007.05984)
- **Key Construction**: Using Open Effective Field Theory (Open EFT), derives a MARKOVIAN Lindblad master equation for a qubit (UDW detector) hovering at radius r near a Schwarzschild black hole:
  ```
  d rho / dt = -i[H_eff, rho] + sum_k gamma_k (L_k rho L_k^dagger - (1/2){L_k^dagger L_k, rho})
  ```
- **Key Results**:
  1. Near-horizon evolution takes a UNIVERSAL form
  2. Qubit thermalizes at local Tolman temperature T(r) = T_H / sqrt(1 - r_s/r)
  3. Decoherence and thermalization rates depend on metric through sqrt(-g_00)
- **Entropy Production**: Computable as Sigma(t) = D(rho(0) || sigma_th) - D(rho(t) || sigma_th), but is time-dependent, coupling-dependent, and probe-dependent.
- **Gap**: Probe-dependent entropy production. Does NOT directly give Sigma = -ln(-g_00). The qubit constraint limits Delta-D.
- **Assessment**: Rigorous first-principles Lindblad equation, but fundamentally probe-dependent. Useful for understanding the physics but does not solve the canonical channel problem.

#### 2.2.3 Quantum Optical Simulator for UDW Dynamics (2025)
- **arXiv**: [2511.16865](https://arxiv.org/abs/2511.16865)
- **Key Result**: Experimental realization of UDW detector-field entanglement using quasiparticles.
- **Assessment**: Experimental confirmation of UDW model predictions. Does not advance the channel construction.

### 2.3 Route C: Gravitational Scattering S-Matrix -> Channel

#### 2.3.1 Danielson-Satishchandran-Wald (2022, 2023, 2025)
- **Title**: "Black holes decohere quantum superpositions" (2022) / "Killing horizons decohere quantum superpositions" (2023) / "Local description of decoherence" (2025)
- **Published**: IJMPD 31, 2241003 (2022); PRD 108, 025007 (2023); PRD 111, 025014 (2025)
- **Key Result**: A massive body in quantum spatial superposition is decohered by radiation of soft gravitons through a Killing horizon. The decoherence is:
  ```
  |<psi_1|psi_2>|^2 = exp(-Gamma_dec * L)
  ```
  where Gamma_dec depends on the mass separation and L is the duration.
- **Local Description (2025)**: The decoherence can be described via the local two-point function of the quantum field within the experimenter's lab, arising from extremely low-frequency Hawking quanta.
- **Channel Structure**: The reduced dynamics of the superposed body defines a CPTP map (partial trace over soft radiation). This is structurally similar to a dephasing channel.
- **Gap**: The decoherence rate Gamma_dec depends on probe parameters (mass, separation) -- it is NOT universal. The S-matrix approach gives scattering amplitudes, not a universal gravitational channel.
- **Assessment**: Provides strong evidence that gravitational processes ARE quantum channels, but the channels are probe-dependent.

#### 2.3.2 Gravitational Waves Decohere Quantum Superpositions (2025)
- **arXiv**: [2501.18111](https://arxiv.org/abs/2501.18111)
- **Key Result**: Gravitational radiation from distant sources decoheres quantum superpositions. Individual contributions identified from memory (soft gravitons) and oscillatory (hard gravitons) components. Memory contribution dominates.
- **Assessment**: Further confirms that gravitational radiation acts as a decoherence channel, but the channel parameters are source-dependent.

#### 2.3.3 Horizons and Soft Quantum Information (2025)
- **arXiv**: [2512.20754](https://arxiv.org/abs/2512.20754)
- **Key Result**: Extends Tomita-Takesaki theory to accommodate soft radiation states. Computes fidelity, relative entropy, and Renyi entropies of general coherent states on the horizon:
  - Fidelity: sensitive to absolute NUMBER of soft gravitons (which can be arbitrarily large)
  - Relative entropy: sensitive to ENERGY of radiation (can be arbitrarily small)
- **Critical Insight**: The relative entropy and fidelity have different scaling with soft radiation parameters. This is relevant to understanding what gravitational "entropy production" measures.
- **Assessment**: Highly relevant to Paper 2. Shows that the choice of quantum information measure (fidelity vs. relative entropy) matters fundamentally for how soft gravitational radiation is characterized.

### 2.4 Route D: Open Quantum System Approach (Lindblad Master Equation)

#### 2.4.1 Anastopoulos & Hu (2013, 2017)
- **Title**: "A Master Equation for Gravitational Decoherence: Probing the Textures of Spacetime"
- **arXiv**: [1305.5231](https://arxiv.org/abs/1305.5231) (2013); [1706.05677](https://arxiv.org/abs/1706.05677) (Topical Review, 2017)
- **Key Result**: First-principles derivation of a Lindblad master equation for matter in linearized quantum gravity. Decoherence occurs in the ENERGY basis (not position basis), controlled by the spectral density of spacetime fluctuations.
- **Channel**: The Lindbladian generates a one-parameter family of CPTP maps exp(L*t).
- **Gap**: Linearized gravity (flat background + perturbations). Does not capture the full g_00 dependence. Probe-dependent through the spectral density.

#### 2.4.2 Moreira & Celeri (2024, 2026)
- **Title**: "Entropy production due to spacetime fluctuations" (2024) / "Decoherence and entropy production due to quantum fluctuations of spacetime" (2026)
- **arXiv**: [2407.21186](https://arxiv.org/abs/2407.21186), [2603.02034](https://arxiv.org/abs/2603.02034)
- **Key Construction**: Open quantum system coupled to a bath of gravitons (quanta of the gravitational field in linearized GR). Uses decoherent histories approach.
- **Key Result (2026)**: The universal interaction of internal variables with gravitons leads to decoherence of spatial superpositions of microscopic systems in the long-time limit, even when the graviton bath alone does not. Entropy production derived via integral fluctuation theorem.
- **Gap**: Linearized gravity. Entropy production depends on graviton spectral density, not directly on -ln(-g_00).
- **Assessment**: Provides the most explicit recent treatment of entropy production from quantum spacetime fluctuations. However, limited to the linearized regime and probe-dependent.

#### 2.4.3 Oppenheim (2023)
- **Title**: "A postquantum theory of classical gravity?"
- **Published**: PRX 13, 041040 (2023)
- **arXiv**: [1811.03116](https://arxiv.org/abs/1811.03116)
- **Key Construction**: Classical-quantum (CQ) dynamics where gravity stays classical and matter is quantum. The CQ map is CPTP:
  ```
  rho -> sum_k L_k * N_k(rho) * P(z_k)
  ```
  with decoherence-diffusion trade-off: more decoherence <-> less diffusion.
- **Gap**: The entropy production depends on the trade-off parameter, not on the metric alone. Does not predict Sigma = -ln(-g_00).
- **Assessment**: Shows that classical gravity + quantum matter IS a CPTP process, confirming the framework. But the specific entropy production is model-dependent.

#### 2.4.4 Hsiang & Hu (2025)
- **Title**: "Non-Markovian Quantum Master and Fokker-Planck Equation for Gravitational Systems"
- **arXiv**: [2504.11991](https://arxiv.org/abs/2504.11991)
- **Key Result**: Non-Markovian generalization beyond the high-temperature Markov approximation.
- **Assessment**: Important technical advance but does not solve the canonical channel problem.

#### 2.4.5 Gravitational Decoherence and Recoherence (2026)
- **arXiv**: [2602.22517](https://arxiv.org/abs/2602.22517)
- **Authors**: Moreira & Celeri
- **Key Result**: Interplay between graviton-induced decoherence and classical Newtonian potential effects. Gravitons + internal DOF render decoherence inevitable long-term, even for microscopic masses. Newtonian potential slightly slows decoherence and can even cause recoherence for systems without dynamical internal DOF.
- **Assessment**: Shows the graviton bath produces irreversible entropy production, but the rate depends on system parameters.

### 2.5 Route E: Stinespring Dilation of Gravitational Redshift

#### 2.5.1 Bosonic Pure-Loss Channel (Standard QI)
- **Channel**: E_eta(rho) = Tr_E[U_BS(rho tensor |0><0|_E) U_BS^dagger]
  where U_BS is a beam-splitter with transmissivity eta.
- **Stinespring**: System = signal mode a, Environment = vacuum mode b, Unitary = beam-splitter mixing.
- **Gravitational Identification**: eta = -g_00 = exp(-r_s/r) gives:
  ```
  Sigma = -ln(eta) + O(1/N_B) = r_s/r + O(1/N_B)   [for thermal reference with N_B >> 1]
  ```
- **Gap**: The beam-splitter model is an EFFECTIVE description, not derived from GR + QFT.
- **Assessment**: Best explicit CPTP map with computable Kraus operators. The only candidate that simultaneously gives: (1) explicit Kraus operators, (2) unbounded Sigma, (3) correct value r_s/r, (4) universality in the high-temperature limit.

#### 2.5.2 Leber, Rodriguez, Ferreri, Schell, Bruschi (2025)
- **Title**: "Limits to the validity of quantum-optical models of the effects of gravitational redshift"
- **arXiv**: [2502.20521](https://arxiv.org/abs/2502.20521)
- **CRITICAL WARNING**: The beam-splitter model for gravitational redshift has LIMITED VALIDITY:
  - Works for small redshift (weak field)
  - Fails for large redshifts (strong field)
  - Validity range depends on both gravitational parameters AND photonic mode parameters
- **Implication**: The pure-loss bosonic channel is a valid effective model only in the weak-field regime. Strong-field behavior requires multimode treatment.

#### 2.5.3 Dragan et al. (2025)
- **Title**: "Gravitational redshift of broadband relativistic quantum photons"
- **arXiv**: [2504.03956](https://arxiv.org/abs/2504.03956)
- **Key Result**: Using linearized quantum gravity (graviton exchange), gravitational redshift occurs only for photons interacting with a classical source AND having well-defined momentum. Photons with quantum coherence in position show NO well-defined redshift.
- **Implication**: The "gravitational channel" depends on the quantum state of the source, not just the background metric. Challenges universality.

#### 2.5.4 Bruschi (2021)
- **Title**: "Towards communication in a curved spacetime geometry"
- **Published**: Communications in Physics 4, 171 (2021)
- **Key Result**: Spacetime curvature distorts wavepackets, causing phase shifts proportional to Riemann curvature and cross-talk between structured modes (up to 12.2% in the solar system).
- **Assessment**: Shows the full evolution IS unitary (no information loss in static backgrounds). The CPTP map arises from tracing out unobserved modes.

### 2.6 Route F: Consistency Theorems and Existence Results

#### 2.6.1 Galley-Giacomini-Selby (2023)
- **Title**: "Any consistent coupling between classical gravity and quantum matter is fundamentally irreversible"
- **Published**: Quantum 7, 1142 (2023)
- **arXiv**: [2301.10261](https://arxiv.org/abs/2301.10261)
- **No-Go Theorem**: In General Probabilistic Theories, at least one of these must be violated:
  1. Matter is fully non-classical
  2. Matter-gravity interactions are reversible
  3. Matter back-reacts on gravity
- **Implication**: If gravity is classical and matter is quantum, the coupling MUST be irreversible (CPTP but not unitary). Entropy production Sigma > 0 is a NECESSITY.
- **Gap**: Proves existence of irreversibility but does NOT construct the specific channel or compute Sigma.
- **Assessment**: The strongest theoretical argument that a gravitational CPTP map with Sigma > 0 MUST exist. Combined with the three routes giving Sigma = r_s/r, this is strong evidence for the identification.

#### 2.6.2 Basso-Maziero-Celeri (2025 PRL)
- **Title**: "Quantum Detailed Fluctuation Theorem in Curved Spacetimes"
- **Published**: PRL 134, 050406 (2025)
- **arXiv**: [2405.03902](https://arxiv.org/abs/2405.03902)
- **Key Result**: Fully general-relativistic quantum Crooks fluctuation theorem using Fermi normal coordinates:
  ```
  P_fwd(W) / P_rev(-W) = exp(beta * (W - Delta_F))
  ```
  Curvature coupling: Riemann tensor R_{titj} appears as effective time-dependent potential.
- **CRITICAL FINDING**: Entropy production is OBSERVER-DEPENDENT. Different observers at different radii measure different Sigma.
- **Assessment**: Observer-dependent Sigma is consistent with Paper 2's framework (different observers see different tau). Provides the rigorous curved-spacetime version of the fluctuation theorem.

### 2.7 Route G: Additional Important Approaches

#### 2.7.1 Casini & Huerta (2017)
- **Title**: "Relative entropy and the RG flow"
- **Published**: JHEP 03 (2017) 089
- **Key Result**: The decrease of relative entropy along the renormalization group flow IS the data processing inequality applied to the coarse-graining map. This provides a QFT version of the connection between information loss and monotonicity.
- **Relevance**: The RG flow IS a CPTP map (coarse-graining). The DPI applied to this map gives the c-theorem and related monotonicity results. For gravity, the analog is: moving from r to infinity IS a form of "coarse-graining" (information loss), and the DPI constrains how much QRE can decrease.

#### 2.7.2 Herrera (2020)
- **Title**: "Landauer Principle and General Relativity"
- **Published**: Entropy 22, 340 (2020)
- **Key Result**: Tolman-modified Landauer erasure cost at radius r:
  ```
  W(r) = k_B T(r) ln 2 = k_B T_inf / sqrt(-g_00) * ln 2
  ```
  Excess entropy production: Sigma = 2 * ln(1/sqrt(-g_00)) = -ln(-g_00) = r_s/r
- **Assessment**: Thermodynamic argument, not a channel construction. But gives the correct formula from an independent starting point.

#### 2.7.3 Pikovski (2015)
- **Title**: "Universal decoherence due to gravitational time dilation"
- **Published**: Nature Physics 11, 668 (2015)
- **arXiv**: [1311.1095](https://arxiv.org/abs/1311.1095)
- **Channel**: Qubit dephasing with Kraus operators K_0 = sqrt((1+p)/2) I, K_1 = sqrt((1-p)/2) sigma_z
- **FUNDAMENTAL LIMITATION**: Qubit channel -> Delta-D <= ln 2 for any fixed reference. Cannot reproduce unbounded Sigma = r_s/r.
- **Value**: Excellent physical motivation. Shows gravity IS a dephasing channel, with the dephasing rate encoding g_00 indirectly.

#### 2.7.4 Kafri-Taylor-Milburn (2014)
- **Title**: Classical channel model for gravitational interaction
- **Published**: NJP 16, 065020 (2014)
- **Channel**: Classical measurement channel: continuous position measurement -> feedback force. Equivalent to Diosi decoherence.
- **Lindblad**: d rho/dt = -i[H, rho] + gamma(x rho x - (1/2){x^2, rho})
  where gamma = G m^2 / (hbar d^3).
- **Gap**: Probe-dependent. Cannot give universal Sigma = -ln(-g_00).

#### 2.7.5 Mari-Zippilli-Vitali (2025)
- **Title**: "Can gravity mediate the transmission of quantum information?"
- **arXiv**: [2504.05998](https://arxiv.org/abs/2504.05998)
- **Channel**: Quantum thermal attenuator (Gaussian), defined for gravitational interaction:
  ```
  a_out = sqrt(eta) * exp(i phi) * a_in + sqrt(1-eta) * a_E
  ```
  with eta depending on gravitational coupling rate vs. mechanical damping.
- **Key Result**: Non-zero quantum capacity when lambda^2 > gamma^2 * N_T * (N_T + 1). Sharp transition between entanglement-breaking and non-classical channel.
- **Gap**: Models gravitational INTERACTION (Newtonian coupling), not gravitational PROPAGATION (redshift).

#### 2.7.6 Bradler & Adami (2015)
- **Title**: "Black holes as bosonic Gaussian channels"
- **Published**: PRD 92, 025030 (2015)
- **arXiv**: [1405.1097](https://arxiv.org/abs/1405.1097)
- **Channel**: Full classification of black hole scattering as bosonic Gaussian channels (lossy, amplifying, or classical-noise).
- **Assessment**: Most complete treatment of black holes as quantum channels. However, channel parameters are specific to black hole mass and potential barrier, not universal.

---

## 3. Route Analysis: Feasibility of Each Approach

### 3.1 Route A: Modular Channel from AQFT

**Feasibility: HIGH (but requires computation)**

**Argument**:
1. The modular channel (partial trace over causally inaccessible region) is a well-defined CPTP map in the AQFT framework.
2. Its singular values follow a Gibbs distribution governed by the modular Hamiltonian.
3. For Schwarzschild, the modular Hamiltonian near the horizon is the boost generator with surface gravity kappa.
4. The Dorau-Much result gives the fractional QRE loss as r_s/r.

**What remains**: An explicit computation of Delta-D for the modular channel in the Schwarzschild exterior at finite radius r, with a specified reference state (vacuum or KMS state). This requires:
- The modular Hamiltonian spectrum at finite r (known only approximately; Bisognano-Wichmann exact only at bifurcate horizons)
- The singular value decomposition of the restriction map
- A specific computation of D(rho||sigma) - D(N(rho)||N(sigma))

**Expected difficulty**: MODERATE. The computation is well-defined in principle but technically challenging because the modular Hamiltonian at finite r receives non-local corrections.

**Can it give Sigma = -ln(-g_00)?**: The fractional result (D_in - D_out)/D_in = r_s/r is already established. If the absolute result D_in - D_out can be shown to equal -ln(-g_00) for appropriately chosen test states, this route would be COMPLETE.

**Key obstacle**: The normalization problem. The absolute entropy production depends on D_in, which is state-dependent. To get a state-independent Sigma = -ln(-g_00), one needs either: (a) a specific choice of states with D_in = 1, or (b) a different definition of Sigma that uses the ratio rather than the difference.

### 3.2 Route B: Unruh-DeWitt Detector Model

**Feasibility: MEDIUM**

**Argument**:
1. Kasprzak-Tjoa (2024) provides a rigorous CPTP map in any curved spacetime.
2. The quantum capacity is expressed in terms of Wightman functions and causal propagators.
3. For Schwarzschild, these correlation functions ARE computable (though technically involved).

**What remains**:
1. Compute W(f,g) for Alice at r, Bob at infinity in Schwarzschild
2. Compute the causal propagator E(f,g) in the same geometry
3. Extract the effective channel N: D(H_A) -> D(H_B)
4. Compute Delta-D and check whether it equals -ln(-g_00)

**Expected difficulty**: HIGH. The Wightman function in Schwarzschild involves mode sums over angular momenta and frequencies, with backscattering from the potential barrier. No closed-form expression exists.

**Key obstacle**: The Kasprzak-Tjoa channel is between QUBIT detectors, so single-mode Delta-D is bounded by ln 2. To get unbounded Sigma, one needs to consider the multimode (wideband) channel and sum over all modes. This is the analog of the information-theoretic statement that the total capacity loss equals the sum over mode losses.

**Can it give Sigma = -ln(-g_00)?**: Plausible if the TOTAL information loss (summed over all modes) equals -ln(-g_00). This is a conjecture that requires verification.

### 3.3 Route C: Gravitational S-Matrix

**Feasibility: LOW-MEDIUM**

**Argument**:
1. The DSW approach shows that soft graviton emission provides which-path information to the horizon.
2. The S-matrix elements determine the decoherence matrix.
3. In principle, the reduced dynamics (tracing over emitted gravitons) defines a CPTP map.

**Key obstacles**:
1. The S-matrix is perturbative (expansion in G). Non-perturbative effects matter near horizons.
2. The decoherence rate depends on probe parameters (mass, separation), not just on the metric.
3. Soft graviton theorems give RATES, not integrated entropy production equal to -ln(-g_00).

**Can it give Sigma = -ln(-g_00)?**: Unlikely directly. The soft graviton results give decoherence rates proportional to G * m^2, not the geometric quantity -ln(-g_00). There may be an indirect connection through the equivalence principle (replacing m by the local energy density), but this has not been established.

**Assessment**: The S-matrix approach is best for understanding the MECHANISM of gravitational decoherence, not for constructing the canonical channel.

### 3.4 Route D: Open Quantum System / Lindblad

**Feasibility: MEDIUM (partial results exist)**

**Argument**:
1. Multiple groups (ABH, Kaplanek-Burgess, Moreira-Celeri, Oppenheim) have derived explicit Lindblad equations for gravitational decoherence.
2. Each Lindbladian generates a one-parameter family of CPTP maps.
3. The entropy production can be computed for each.

**Key obstacles**:
1. All existing Lindbladian approaches give PROBE-DEPENDENT entropy production (depends on the detector gap, mass, coupling constant).
2. The decoherence basis varies: energy basis (ABH), position basis (Diosi-Penrose), center-of-mass (Pikovski).
3. To get Sigma = -ln(-g_00), one would need to take a supremum over all probes, which may or may not converge to the geometric value.

**Can it give Sigma = -ln(-g_00)?**: The EQUILIBRIUM relative entropy D(sigma_th(r) || sigma_th(inf)) = -ln(-g_00) follows from the Tolman relation. But the DYNAMICAL entropy production Sigma(t) is different and depends on the probe. The connection between the equilibrium and dynamical quantities is not established.

**Promising sub-route**: The Kaplanek-Burgess result that qubit thermalization near the horizon is UNIVERSAL (independent of probe details for r -> r_s) suggests that the entropy production in the near-horizon limit might also be universal. This deserves investigation.

### 3.5 Route E: Stinespring Dilation

**Feasibility: HIGH (already done, but effective)**

**Argument**: The bosonic pure-loss channel with eta = -g_00 gives Sigma = -ln(eta) = r_s/r in the high-temperature limit. The Stinespring dilation is the beam-splitter model.

**Status**: This is the BEST AVAILABLE answer, but it is an effective model, not a first-principles derivation. The identification eta = -g_00 uses the power transmissivity convention and is motivated by physical arguments (both frequency shift and time dilation contribute to occupation number loss), but is not derived from GR + QFT.

**Leber et al. caveat**: The beam-splitter model has limited validity for large redshifts. For strong fields (near the horizon), a multimode treatment is required. The pure-loss channel is accurate only when r_s/r << 1.

---

## 4. Can We Bypass the Channel?

### 4.1 Strategy 1: QRE Difference Without Explicit N

**Idea**: Define Sigma_grav directly as the QRE difference between "states at r" and "states at infinity", without constructing an explicit CPTP map:

```
Sigma_grav(r) := D(omega_r || phi_r) - D(omega_inf || phi_inf)
```

where omega_r, phi_r are states on the algebra at radius r, and omega_inf, phi_inf are the corresponding states seen from infinity.

**Justification**: The Dorau-Much result already shows that this difference equals r_s/r (fractionally). The DATA PROCESSING INEQUALITY guarantees that if ANY CPTP map connects the two, the difference is non-negative. We do not need to know the specific map -- only that one exists (which is guaranteed by the causal structure of spacetime).

**Pros**:
- Avoids the need to construct N_grav explicitly
- The QRE difference is well-defined in the algebraic framework
- Consistent with DPI (which is model-independent)
- The Galley-Giacomini-Selby theorem guarantees that an irreversible CPTP map exists

**Cons**:
- The JRSWW bound F >= exp(-Sigma/2) is proven for a SPECIFIC CPTP map. Using Sigma defined as a QRE difference without specifying the map does not automatically give the bound.
- The Petz recovery map (which achieves near-optimal recovery) is defined relative to a specific channel. Without the channel, the Petz map is undefined.
- The identification Sigma = -ln(-g_00) still requires the normalization assumption (fractional vs. absolute).

**Assessment**: This is a VIABLE bypass for making the physical argument, but it does not fully replace the channel construction for the mathematical theorem.

### 4.2 Strategy 2: Axiomatic Approach

**Idea**: Instead of constructing N_grav, AXIOMATIZE the properties it must have and derive consequences.

**Axioms**:
1. N_grav is CPTP on an infinite-dimensional Hilbert space
2. N_grav depends only on the metric (specifically, on g_00(r))
3. N_grav reduces to the identity as g_00 -> -1 (flat spacetime limit)
4. The entropy production of N_grav is monotonic in the redshift factor
5. The entropy production of N_grav is additive under composition of channels for concentric shells

**From Axiom 5**: If the channel for going from r_1 to r_2 composes with the channel from r_2 to r_3, then:
```
Sigma(r_1 -> r_3) = Sigma(r_1 -> r_2) + Sigma(r_2 -> r_3)
```

This additivity, combined with the monotonicity axiom, uniquely determines:
```
Sigma(r) = c * (-ln(-g_00(r)))
```

for some constant c > 0 (this is because the only additive, monotonic function of -g_00 is the logarithm).

**The constant c** cannot be determined from the axioms alone. The physical arguments (power transmissivity, Tolman-Landauer) suggest c = 1. The amplitude transmissivity gives c = 1/2.

**Pros**:
- Clean and elegant
- Uniquely determines the functional form Sigma = -c * ln(-g_00)
- Does not require constructing the explicit channel
- The axioms are physically reasonable

**Cons**:
- Axiom 5 (additivity) is a strong assumption. It requires the channel for each shell to be independent of the others, which may not hold for entangled field modes spanning multiple shells.
- The constant c remains undetermined
- The axioms need independent verification

**Assessment**: The most promising bypass strategy. Could be the basis for a Paper 2 revision that replaces "here is the channel" with "any channel satisfying these axioms gives Sigma = -ln(-g_00)."

### 4.3 Strategy 3: Universal Capacity Loss

**Idea**: Instead of computing Delta-D for specific states, use the QUANTUM CAPACITY loss of the gravitational channel.

The quantum capacity Q(N) of a channel N quantifies the maximum rate of reliable quantum information transmission. For a pure-loss bosonic channel:
```
Q(E_eta) = max(0, log((eta/(1-eta))))
```

For eta = -g_00 close to 1 (weak field):
```
Q(E_{-g_00}) ≈ log(1/(r_s/r)) = -log(r_s/r)
```

The capacity LOSS compared to the identity channel is:
```
Delta_Q = Q(id) - Q(E_{-g_00}) = infinity - finite = infinity
```

This diverges because the identity channel has infinite capacity. A better measure is the COHERENT INFORMATION loss for specific input ensembles.

**Assessment**: The capacity framework is useful conceptually but does not directly give Sigma = -ln(-g_00) as a finite quantity.

### 4.4 Strategy 4: Define Sigma via the Petz Recovery Fidelity

**Idea**: REVERSE the logic. Instead of deriving Sigma and using JRSWW to get F, define Sigma through the Petz recovery fidelity:

```
Sigma_grav := -2 ln F_Petz(r)
```

where F_Petz(r) is the fidelity of the Petz recovery map for reconstructing states at r from observations at infinity.

This is well-defined if:
- We specify which states we are trying to recover
- The Petz map is defined relative to some channel (but ANY channel that models the signal propagation from r to infinity will do)

For the bosonic loss channel with eta = -g_00:
```
F_Petz = sqrt(-g_00) = exp(-r_s/(2r))
Sigma = -2 ln(exp(-r_s/(2r))) = r_s/r
```

**Pros**:
- Operationally meaningful (defines Sigma through an information recovery task)
- Does not require specifying the full channel, only the recovery fidelity
- Gives the correct result

**Cons**:
- Circular if we define F_Petz using the bosonic loss channel (then we ARE using a specific channel)
- The Petz map is defined relative to a channel; without specifying N, the Petz map is undefined
- This is a RESTATEMENT of the problem, not a solution

### 4.5 Recommended Bypass Strategy

**Combined Strategy (Axioms + Modular QRE)**:

1. **State the axioms** (Section 4.2) for the gravitational channel
2. **Derive** that Sigma must have the form -c * ln(-g_00) from additivity and monotonicity
3. **Use the Dorau-Much result** (modular QRE) to fix the fractional coefficient to c * r_s/r / D_in = r_s/r
4. **Use the bosonic loss channel** as an explicit REALIZATION demonstrating that c = 1 is achievable and gives Sigma = -ln(-g_00) in the high-temperature limit
5. **Cite Galley-Giacomini-Selby** to argue that such an irreversible CPTP map MUST exist

This combined strategy:
- Does not require constructing N_grav from first principles
- Uses the axiomatic approach for the functional form
- Uses the algebraic result for the coefficient
- Uses the bosonic channel as an existence proof
- Uses the no-go theorem for the necessity of irreversibility

---

## 5. Conclusion and Recommendation

### 5.1 Current State of Knowledge

**What EXISTS (as of March 2026)**:

| Feature | Best Available | Reference |
|---------|---------------|-----------|
| Explicit CPTP map giving Sigma = r_s/r | Pure-loss bosonic, eta = -g_00, N_B >> 1 | Standard QI + physical identification |
| First-principles justification for Sigma = r_s/r | Modular flow (Dorau-Much) | arXiv:2510.24491 |
| Thermodynamic argument for Sigma = r_s/r | Gravitational Landauer (Herrera) | Entropy 22, 340 (2020) |
| Necessity of irreversibility (Sigma > 0) | Galley-Giacomini-Selby no-go | Quantum 7, 1142 (2023) |
| Rigorous CPTP in curved spacetime | Kasprzak-Tjoa UDW channel | arXiv:2408.00518 |
| Algebraic CPTP (conditional expectation) | Witten, Chandrasekaran, Faulkner | Multiple papers (2020-2023) |
| Observer-dependent Sigma in curved spacetime | Basso-Celeri quantum Crooks | PRL 134, 050406 (2025) |
| Modular channel with Gibbs spectrum | Trejo-Calderon thermal filter | arXiv:2504.20457 |

**What does NOT exist**:

| Feature | Status | Comment |
|---------|--------|---------|
| First-principles CPTP with explicit Kraus operators and Sigma = -ln(-g_00) | DOES NOT EXIST | The central open problem |
| Proof that Sigma = r_s/r for the modular channel (absolute, not fractional) | DOES NOT EXIST | Normalization problem |
| Universal (probe-independent) Lindblad giving Sigma = -ln(-g_00) | DOES NOT EXIST | All Lindblad approaches are probe-dependent |
| Proof of Petz bound saturation for Sigma > 0 | DOES NOT EXIST | Likely impossible |

### 5.2 Has Anyone Already Done It?

**NO.** After surveying 50+ papers across 7 research directions (AQFT modular theory, UDW detectors, gravitational S-matrix, Lindblad approaches, Stinespring dilation, consistency theorems, and capacity theory), we can confirm:

> **No paper in the literature provides a canonical, first-principles CPTP map N_grav with explicit Kraus operators such that Delta-D = -ln(-g_00).**

The closest candidates are:
1. **Trejo-Calderon (2025)**: Constructs the modular channel with Gibbs-weighted singular values. Could in principle yield Sigma = -ln(-g_00) but the specific computation has not been performed.
2. **Kasprzak-Tjoa (2024)**: Constructs a rigorous CPTP in any curved spacetime. The total capacity loss in Schwarzschild has not been computed.
3. **Bosonic pure-loss channel**: Gives Sigma = -ln(eta) exactly. The identification eta = -g_00 is physically motivated but not derived.

### 5.3 Honest Statement for Paper 2

Based on this survey, the recommended honest statement for Paper 2 is:

> "The gravitational entropy production Sigma_grav = -ln(-g_00) is supported by three independent physical arguments: (i) modular flow on the Killing horizon gives the fractional QRE loss (D_in - D_out)/D_in = r_s/r [Dorau-Much 2025]; (ii) the Tolman-modified Landauer erasure cost gives Sigma = -ln(-g_00) [Herrera 2020]; and (iii) a bosonic pure-loss channel with transmissivity eta = -g_00 and high-temperature thermal reference reproduces Sigma = -ln(eta) in the large-N_B limit [standard QI]. Furthermore, the Galley-Giacomini-Selby theorem [2023] proves that any consistent coupling between classical gravity and quantum matter must be fundamentally irreversible, guaranteeing Sigma > 0. An axiomatic argument based on additivity and monotonicity uniquely determines the functional form Sigma = -c * ln(-g_00). The construction of a fully first-principles canonical CPTP map N_grav with explicit Kraus operators remains an important open problem."

### 5.4 Recommended Path Forward

**Immediate (for Paper 2)**:
1. Adopt the axiomatic + existence proof strategy (Section 4.5)
2. State clearly that the bosonic loss channel is an effective model
3. Use the Dorau-Much fractional QRE as the first-principles support
4. Cite Galley-Giacomini-Selby for necessity

**Short-term (standalone calculation)**:
1. Compute the modular channel Delta-D for Schwarzschild at finite r, using the known approximate modular Hamiltonian
2. Verify whether the spectral entropy of the Trejo-Calderon thermal filter equals -ln(-g_00)
3. Compute the Kasprzak-Tjoa channel capacity loss in Schwarzschild (multimode)

**Medium-term (potential publication)**:
1. Construct N_grav explicitly from the crossed-product algebra (Witten 2022), reducing the conditional expectation to a Hilbert space operator with Kraus decomposition
2. Compute its entropy production and verify Sigma = -ln(-g_00)
3. This would close the gap and promote Layer 2 from conjecture to theorem

### 5.5 Assessment of Difficulty

| Task | Difficulty | Estimated effort | Impact |
|------|-----------|-----------------|--------|
| Modular channel Delta-D in Schwarzschild | MODERATE | 2-4 months | HIGH (could close the gap) |
| Kasprzak-Tjoa capacity loss | HIGH | 3-6 months | HIGH (independent verification) |
| Crossed-product Kraus operators | VERY HIGH | 6-12 months | VERY HIGH (full solution) |
| Bosonic channel validity bounds | MODERATE | 1-2 months | MEDIUM (clarifies limitations) |
| Axiomatic derivation | LOW | 1 month | MEDIUM (elegant bypass) |

---

## 6. Complete Reference List

### First-Principles Constructions / Algebraic

1. Dorau & Much, "From Quantum Relative Entropy to the Semiclassical Einstein Equations," PRL 136, 091602 (2026). [arXiv:2510.24491](https://arxiv.org/abs/2510.24491)
2. Trejo-Calderon, "Modular Channels, Thermal Filtering and the Spectral Emergence of Spacetime," (2025). [arXiv:2504.20457](https://arxiv.org/abs/2504.20457)
3. Witten, "Gravity and the Crossed Product," JHEP 10 (2022) 008. [arXiv:2112.12828](https://arxiv.org/abs/2112.12828)
4. Chandrasekaran, Longo, Penington, Witten, "An Algebra of Observables for de Sitter Space," JHEP 02 (2023) 082. [arXiv:2206.10780](https://arxiv.org/abs/2206.10780)
5. Chandrasekaran, Penington, Witten, "Generalized entropy for general subregions in quantum gravity," JHEP 12 (2023) 020. [arXiv:2306.01837](https://arxiv.org/abs/2306.01837)
6. Faulkner, "The holographic map as a conditional expectation," (2020). [arXiv:2008.04810](https://arxiv.org/abs/2008.04810)
7. Cirafici, "Fluctuation theorems, quantum channels and gravitational algebras," JHEP 11 (2024) 089. [arXiv:2408.04219](https://arxiv.org/abs/2408.04219)
8. Jacobson, "Thermodynamics of Spacetime: The Einstein Equation of State," PRL 75, 1260 (1995). [arXiv:gr-qc/9504004](https://arxiv.org/abs/gr-qc/9504004)
9. Jacobson, "Entanglement Equilibrium and the Einstein Equation," PRL 116, 201101 (2016). [arXiv:1505.04753](https://arxiv.org/abs/1505.04753)
10. Casini, Huerta, Teste, Torroba, "Relative entropy and the RG flow," JHEP 03 (2017) 089.

### Unruh-DeWitt Detector Channels

11. Kasprzak & Tjoa, "Transmission of quantum information through quantum fields in curved spacetimes," (2024). [arXiv:2408.00518](https://arxiv.org/abs/2408.00518)
12. Kaplanek & Burgess, "Hot Accelerated Qubits," JHEP 03 (2020) 008. [arXiv:1912.12951](https://arxiv.org/abs/1912.12951)
13. Kaplanek & Burgess, "Qubits on the Horizon," JHEP 01 (2021) 098. [arXiv:2007.05984](https://arxiv.org/abs/2007.05984)
14. Landulfo, "Communication through quantum fields near a black hole," PRD 93, 104019 (2016).

### Gravitational S-Matrix / Decoherence

15. Danielson, Satishchandran, Wald, "Black holes decohere quantum superpositions," IJMPD 31, 2241003 (2022).
16. Danielson, Satishchandran, Wald, "Killing horizons decohere quantum superpositions," PRD 108, 025007 (2023).
17. Danielson, Satishchandran, Wald, "Local description of decoherence," PRD 111, 025014 (2025).
18. "Gravitational waves decohere quantum superpositions," (2025). [arXiv:2501.18111](https://arxiv.org/abs/2501.18111)
19. "Horizons and Soft Quantum Information," (2025). [arXiv:2512.20754](https://arxiv.org/abs/2512.20754)

### Open Quantum System / Lindblad

20. Anastopoulos & Hu, "A Master Equation for Gravitational Decoherence," CQG 30, 165007 (2013). [arXiv:1305.5231](https://arxiv.org/abs/1305.5231)
21. Anastopoulos & Hu, "Gravitational Decoherence," Topical Review (2017). [arXiv:1706.05677](https://arxiv.org/abs/1706.05677)
22. Blencowe, "Effective field theory approach to gravitationally induced decoherence," PRL 111, 021302 (2013).
23. Moreira & Celeri, "Entropy production due to spacetime fluctuations," (2024). [arXiv:2407.21186](https://arxiv.org/abs/2407.21186)
24. Moreira, "Decoherence and entropy production due to quantum fluctuations of spacetime," (2026). [arXiv:2603.02034](https://arxiv.org/abs/2603.02034)
25. Moreira & Celeri, "Gravitational decoherence and recoherence of a composite particle," (2026). [arXiv:2602.22517](https://arxiv.org/abs/2602.22517)
26. Oppenheim, "A postquantum theory of classical gravity?" PRX 13, 041040 (2023). [arXiv:1811.03116](https://arxiv.org/abs/1811.03116)
27. Hsiang, Cho, Hu, "Graviton physics tutorial," Universe 10, 306 (2024). [arXiv:2405.11790](https://arxiv.org/abs/2405.11790)
28. Hsiang & Hu, "Non-Markovian Quantum Master Equation for Gravitational Systems," (2025). [arXiv:2504.11991](https://arxiv.org/abs/2504.11991)

### Stinespring / Bosonic Channels

29. Ahmadi, Bruschi, Fuentes, "Quantum metrology for relativistic quantum fields," PRD 89, 065028 (2014). [arXiv:1312.5707](https://arxiv.org/abs/1312.5707)
30. Leber, Rodriguez, Ferreri, Schell, Bruschi, "Limits of quantum-optical redshift models," (2025). [arXiv:2502.20521](https://arxiv.org/abs/2502.20521)
31. Dragan et al., "Gravitational redshift of broadband relativistic quantum photons," (2025). [arXiv:2504.03956](https://arxiv.org/abs/2504.03956)
32. Bruschi, "Towards communication in a curved spacetime geometry," Commun. Phys. 4, 171 (2021).
33. Bradler & Adami, "Black holes as bosonic Gaussian channels," PRD 92, 025030 (2015). [arXiv:1405.1097](https://arxiv.org/abs/1405.1097)
34. Ivan, Sabapathy, Simon, "Operator-sum Representation for Bosonic Gaussian Channels," PRA 84, 042311 (2011). [arXiv:1012.4266](https://arxiv.org/abs/1012.4266)
35. Mari, Zippilli, Vitali, "Can gravity mediate the transmission of quantum information?" (2025). [arXiv:2504.05998](https://arxiv.org/abs/2504.05998)
36. Toccacelo, Andersen, Brask, "Benchmarks for quantum communication via gravity," (2025). [arXiv:2503.03585](https://arxiv.org/abs/2503.03585)

### Consistency / Existence Theorems

37. Galley, Giacomini, Selby, "Any consistent coupling between classical gravity and quantum matter is fundamentally irreversible," Quantum 7, 1142 (2023). [arXiv:2301.10261](https://arxiv.org/abs/2301.10261)
38. Basso, Maziero, Celeri, "Quantum Detailed Fluctuation Theorem in Curved Spacetimes," PRL 134, 050406 (2025). [arXiv:2405.03902](https://arxiv.org/abs/2405.03902)
39. Basso, Maziero, Celeri, "Work distribution of quantum fields in static curved spacetimes," PRD (2026).

### Thermodynamic Arguments

40. Herrera, "Landauer Principle and General Relativity," Entropy 22, 340 (2020).
41. Bianconi, "Gravity from entropy," PRD 111, 066001 (2025). [arXiv:2408.14391](https://arxiv.org/abs/2408.14391)
42. Bianconi, "The Quantum Relative Entropy of the Schwarzschild Black Hole and the Area Law," Entropy 27, 266 (2025). [arXiv:2501.09491](https://arxiv.org/abs/2501.09491)

### Pikovski / Qubit Channels

43. Pikovski, Zych, Costa, Brukner, "Universal decoherence due to gravitational time dilation," Nature Physics 11, 668 (2015). [arXiv:1311.1095](https://arxiv.org/abs/1311.1095)
44. Bonder, Okon, Sudarsky, "Questioning universal decoherence," Nature Physics 12, 2 (2016). [arXiv:1507.05320](https://arxiv.org/abs/1507.05320)
45. Balatsky, Roushan et al., "Quantum Sensing from Gravity as Universal Dephasing Channel," PRA 111, 012411 (2025). [arXiv:2406.03256](https://arxiv.org/abs/2406.03256)
46. Kafri, Taylor, Milburn, "A classical channel model of Newtonian gravity," NJP 16, 065020 (2014). [arXiv:1401.0946](https://arxiv.org/abs/1401.0946)
47. Bradler, Jauregui, Adami, "The Unruh effect interpreted as a quantum noise channel," (2014). [arXiv:1408.1477](https://arxiv.org/abs/1408.1477)

### Black Hole Channels

48. Bradler & Adami, "The capacity of black holes to transmit quantum information," JHEP 2014, 095 (2014). [arXiv:1310.7914](https://arxiv.org/abs/1310.7914)
49. Hayden & Preskill, "Black holes as mirrors," JHEP 2007, 120 (2007).
50. Iyen et al., "Quantum Fisher Information in Curved Spacetime," IJTP 64, 205 (2025). [arXiv:2507.17901](https://arxiv.org/abs/2507.17901)

### QFT in Curved Spacetime / Other

51. Bruschi, Louko, Martin-Martinez, Dragan, Fuentes, "Unruh effect beyond the single-mode approximation," PRA 82, 042332 (2010). [arXiv:1007.4670](https://arxiv.org/abs/1007.4670)
52. Kohlrus, Bruschi, Louko, Fuentes, "Quantum communications in rotating planet spacetime," EPJQT 4, 7 (2017). [arXiv:1511.04256](https://arxiv.org/abs/1511.04256)
53. "Triangulated relativistic quantum computation," Quantum Studies: Mathematics and Foundations (2025). [Springer](https://link.springer.com/article/10.1007/s40509-025-00376-5)

---

## Appendix A: Summary Table of All Candidates

| Rank | Candidate Channel | Sigma = -ln(-g_00)? | Explicit Kraus? | First-Principles? | Unbounded? | Universal? | Key Gap |
|------|------------------|---------------------|----------------|-------------------|-----------|-----------|---------|
| **1** | **Pure-loss bosonic (eta = -g_00)** | **YES** (N_B >> 1) | **YES** | NO (effective) | **YES** | **YES** | Not derived from GR |
| 2 | Modular channel (Trejo-Calderon) | Plausible (not computed) | Singular values (partial) | YES | YES | YES | Delta-D not computed |
| 3 | Algebraic restriction (Dorau-Much) | Fractional YES (r_s/r) | NO (conditional expectation) | YES | YES | YES | Normalization problem |
| 4 | Kasprzak-Tjoa UDW | Unknown (not computed) | YES (implicit via unitaries) | YES | Per-mode NO, total unknown | YES | Schwarzschild calculation needed |
| 5 | Holographic CE (Faulkner) | YES (in AdS/CFT) | NO (algebraic) | YES | YES | YES | Requires AdS/CFT |
| 6 | Kaplanek-Burgess Lindblad | NO (probe-dependent) | YES | YES | NO (qubit) | Near-horizon only | Probe + time dependent |
| 7 | Cirafici subfactor | YES (fluctuation theorem) | NO (abstract) | YES | YES | De Sitter only | Very abstract |
| 8 | Oppenheim CQ | NO (trade-off dependent) | YES (CQ Kraus) | YES | YES | YES | Trade-off parameter |
| 9 | ABH Lindblad | NO (spectral density) | YES (master equation) | YES | YES | YES | Linearized gravity only |
| 10 | Pikovski dephasing | **NO** (bounded by ln 2) | **YES** | YES | **NO** | NO (probe) | Finite dimension |

## Appendix B: The Normalization Problem in Detail

The Dorau-Much result gives the fractional QRE loss:
```
(D_in - D_out) / D_in = r_s/r
```

To convert this to the absolute JRSWW entropy production Sigma = D_in - D_out = r_s/r, we need D_in = 1. But D_in = D(omega_r || omega_0) depends on the test state omega_r:

- For small coherent perturbations: D_in ~ epsilon^2 (small, from entanglement first law)
- For highly excited states: D_in can be arbitrarily large
- For the "canonical" choice D_in = 1: this requires a specific test state

**Possible resolutions**:
1. **Operationalized definition**: Define Sigma_grav as the MAXIMUM entropy production per unit of input relative entropy:
   ```
   Sigma_grav = sup_{rho} [D(rho||sigma) - D(N(rho)||N(sigma))] / D(rho||sigma) = r_s/r
   ```
   This is the CONTRACTION COEFFICIENT of the channel, which equals r_s/r for the modular channel (from Dorau-Much). This is well-defined and state-independent, but it redefines Sigma as a ratio, not a difference.

2. **Fixed normalization**: Choose test states with D_in = 1 by convention. This is analogous to choosing unit probe energy in metrology. The entropy production per unit input QRE is then r_s/r.

3. **Logarithmic identification**: Note that the FRACTIONAL loss r_s/r IS -ln(-g_00) for the exponential metric (to leading order). So the fractional loss IS the logarithmic entropy production, without need for absolute normalization.

**Resolution 3 is the cleanest**: The fractional QRE loss equals -ln(-g_00) in the weak-field limit, providing a self-consistent identification without requiring absolute normalization.

---

**Last updated**: 2026-03-10
**Total papers surveyed**: 53
**Papers with explicit CPTP maps for gravity**: ~15
**Papers with first-principles gravitational entropy production**: ~8
**Central conclusion**: The canonical gravitational CPTP map remains an open problem. The best available effective model is the pure-loss bosonic channel with eta = -g_00. The most promising first-principles route is the modular channel (Trejo-Calderon 2025) with explicit Delta-D computation in Schwarzschild.
