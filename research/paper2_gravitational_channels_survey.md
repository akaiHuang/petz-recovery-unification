# Survey: CPTP Quantum Channels for Gravitational Processes

**Date**: 2026-03-10
**Purpose**: Comprehensive literature survey of papers defining explicit CPTP quantum channels for gravitational processes, to identify the best candidate for the "canonical gravitational channel" N_grav with entropy production Sigma_grav = -ln(-g_00).

---

## Table of Contents

1. [Category A: Gravitational Dephasing Channels (Pikovski et al.)](#category-a)
2. [Category B: Bosonic Gaussian Channels in Curved Spacetime](#category-b)
3. [Category C: Unruh / Rindler Channels](#category-c)
4. [Category D: Hawking Radiation as Quantum Channel](#category-d)
5. [Category E: Relativistic Quantum Channels (Ahmadi, Fuentes, Mann et al.)](#category-e)
6. [Category F: Graviton-Mediated Decoherence Channels](#category-f)
7. [Category G: Channel Capacity in Gravitational Fields](#category-g)
8. [Category H: Modular Flow / Algebraic Channels](#category-h)
9. [Category I: Thermal Channels as Gravitational Analogs](#category-i)
10. [Category J: Quantum Crooks / Fluctuation Theorems in Curved Spacetime](#category-j)
11. [Category K: Gravity-Quantum Consistency and Irreversibility](#category-k)
12. [Category L: Quantum Metrology in Gravitational Fields](#category-l)
13. [Evaluation: Best Candidate Channels](#evaluation)
14. [Relevance to Paper 3 (Weak Field / Galactic Scale)](#paper3-relevance)
15. [Open Problems and Conclusions](#conclusions)

---

## Category A: Gravitational Dephasing Channels (Pikovski et al.) {#category-a}

### A1. Pikovski, Zych, Costa, Brukner (2015)
- **Title**: "Universal decoherence due to gravitational time dilation"
- **Journal**: Nature Physics 11, 668-672 (2015)
- **arXiv**: [1311.1095](https://arxiv.org/abs/1311.1095)
- **Channel Type**: Qubit dephasing channel
- **Kraus Operators**:
  ```
  K_0 = sqrt((1+p)/2) * I
  K_1 = sqrt((1-p)/2) * sigma_z
  ```
  where p = D(t) = exp(-Gamma * t^2) is the decoherence factor, with Gamma depending on the internal energy spread Delta_E and gravitational potential difference.
- **Action**: N_p(rho) = (1/2)(I + p*r_x*sigma_x + p*r_y*sigma_y + r_z*sigma_z)
  Off-diagonal elements multiplied by p; diagonal preserved.
- **Entropy Production**: Delta-D = H_bin((1+p)/2) <= ln 2 ~ 0.693 for sigma = I/2 (or any diagonal reference). **BOUNDED** -- cannot reproduce Sigma = r_s/r for strong fields.
- **Assessment for N_grav**: **FAILS** as canonical gravitational channel because:
  1. Qubit channel -- entropy production bounded by ln 2
  2. Probe-dependent (depends on Delta_E of the composite particle)
  3. Does not give universal Sigma = -ln(-g_00)
  4. Time-dependent (Gaussian decay), not a steady-state channel

### A2. Bonder, Okon, Sudarsky (2016) -- Critique
- **Title**: "Questioning universal decoherence due to gravitational time dilation"
- **Journal**: Nature Physics 12, 2 (2016)
- **arXiv**: [1507.05320](https://arxiv.org/abs/1507.05320)
- **Key Critique**: The Pikovski decoherence depends on the frame of reference. The "universality" claim is questioned -- the dephasing is relative to the chosen basis and the specific internal degrees of freedom.
- **Implication**: Even as an effective channel, the Pikovski CPTP map is frame-dependent.

### A3. Pang, Chen (2016)
- **Title**: "Universal Decoherence under Gravity: A Perspective through the Equivalence Principle"
- **Journal**: Physical Review Letters 117, 090401 (2016)
- **Key Result**: Clarifies the Pikovski mechanism via the equivalence principle. Does not resolve the fundamental limitation (qubit, bounded entropy production).

### A4. Balatsky, Roushan et al. (2025)
- **Title**: "Quantum Sensing from Gravity as Universal Dephasing Channel for Qubits"
- **Journal**: Physical Review A 111, 012411 (2025)
- **arXiv**: [2406.03256](https://arxiv.org/abs/2406.03256)
- **Channel Type**: Qubit dephasing, specialized to transmon qubits
- **Kraus Operators**:
  ```
  Sigma_k(t) = diag(1, e^{i*theta_k})
  ```
  where theta_k = (t/c^2) * delta_Phi_k * omega_k is the gravitational phase.
- **CPTP Map**: Lambda: rho -> Sigma(t) rho Sigma(t)^dagger (pure dephasing)
- **Connection to g_00**: The Hamiltonian includes H = hbar*omega*(1 + Phi(x)/c^2)|1><1|, so the dephasing rate encodes g_00 indirectly.
- **Assessment**: Same structural limitation as Pikovski (qubit, bounded Delta-D). However, provides an explicit, experimentally accessible gravitational channel with measured sensitivity delta_g/g ~ 10^-7.
- **Key value for Paper 2**: Confirms that gravity IS a dephasing channel in practice, even though qubit models cannot reproduce unbounded Sigma.

---

## Category B: Bosonic Gaussian Channels in Curved Spacetime {#category-b}

### B1. Pure-Loss Bosonic Channel (Standard QI Result)
- **References**: Holevo (2019), Serafini (2017), Wilde (2017)
- **Channel Definition**:
  ```
  E_eta(rho) = Tr_E[U_BS (rho tensor |0><0|_E) U_BS^dagger]
  ```
  where U_BS is a beam-splitter unitary with transmissivity eta:
  ```
  a_out = sqrt(eta) * a_in + sqrt(1-eta) * a_vac
  ```
- **Kraus Operators** (Fock basis):
  ```
  A_l = sum_{n=l}^{infinity} sqrt(C(n,l)) * eta^{(n-l)/2} * (1-eta)^{l/2} * |n-l><n|
  ```
  where C(n,l) = n! / (l! (n-l)!)
- **Entropy Production** (coherent input |alpha>, thermal reference sigma_{N_B}):
  ```
  Sigma = N_s * ln(1 + 1/N_B) - eta*N_s * ln(1 + 1/(eta*N_B))
        + ln(1+N_B) - ln(1+eta*N_B)
  ```
  For N_B >> 1: **Sigma = -ln(eta) + O(1/N_B)** independent of input state.
- **Gravitational Identification**: eta = -g_00 = exp(-r_s/r) gives Sigma = r_s/r.
- **Assessment**: **BEST CANDIDATE** for canonical gravitational channel. Explicit CPTP map, infinite-dimensional (can produce unbounded Sigma), and reproduces -ln(-g_00) in the high-temperature limit.
- **Caveats**:
  1. eta = -g_00 (power transmissivity) is a modeling choice, not a derivation
  2. High-temperature limit N_B >> 1 required (motivated by Tolman effect near horizons)
  3. Beam-splitter model is effective, not first-principles

### B2. Operator-Sum Representation for Bosonic Gaussian Channels
- **Authors**: Ivan, Sabapathy, Simon
- **Title**: "Operator-sum Representation for Bosonic Gaussian Channels"
- **Journal**: Physical Review A 84, 042311 (2011)
- **arXiv**: [1012.4266](https://arxiv.org/abs/1012.4266)
- **Key Result**: Provides the complete Choi-Jamiolkowski form and Kraus operators for ALL single-mode bosonic Gaussian channels, including thermal attenuator, amplifier, and classical noise channels.
- **Relevance**: Gives the explicit Kraus operators for the thermal attenuator channel (the gravitational channel candidate).

### B3. Leber, Rodriguez, Ferreri, Schell, Bruschi (2025)
- **Title**: "Limits to the validity of quantum-optical models of the effects of gravitational redshift on photonic quantum states"
- **arXiv**: [2502.20521](https://arxiv.org/abs/2502.20521)
- **Key Warning**: The beam-splitter/loss model for gravitational redshift has LIMITED validity. The model is consistent for small redshift but FAILS for larger redshifts. The range depends on both gravitational parameters AND photonic mode parameters.
- **Implication**: The pure-loss bosonic channel as a gravitational model is an effective description valid only for certain observational scenarios. This is a fundamental caveat for using it as the "canonical" gravitational channel.

---

## Category C: Unruh / Rindler Channels {#category-c}

### C1. Bradler, Adami (2014) -- Unruh Channel
- **Title**: "The capacity of black holes to transmit quantum information"
- **Journal**: Journal of High Energy Physics 2014, 095 (2014)
- **arXiv**: [1310.7914](https://arxiv.org/abs/1310.7914)
- **Channel Definition**: The Unruh channel emerges from tracing over the Rindler wedge inaccessible to an accelerated observer. For a Dirac field:
  ```
  N_Unruh(rho) = Tr_{Rindler II}[U rho U^dagger]
  ```
  where U implements the Bogoliubov transformation between Minkowski and Rindler modes.
- **Key Result**: Perfectly reflecting black hole channel is closely related to Unruh channel; capacity is non-vanishing.
- **Assessment**: Well-defined CPTP map but probe-dependent (depends on the field modes used). Temperature is the Unruh temperature T = a/(2*pi).

### C2. Bradler, Jauregui, Adami (2014) -- Unruh as Amplitude Damping
- **Title**: "The Unruh effect interpreted as a quantum noise channel"
- **arXiv**: [1408.1477](https://arxiv.org/abs/1408.1477)
- **Channel Type**: Amplitude-damping-like channel
- **Key Finding**: The Unruh effect produces an amplitude-damping-like channel associated with zero temperature, despite the Unruh temperature being non-zero. The Bloch sphere contracts by a finite factor rather than collapsing to a point.
- **Assessment**: Explicit CPTP map with Kraus representation. However, it acts on a qubit (modal projection), so it shares the bounded entropy production limitation.

### C3. Sinha, Ganapathy (2016)
- **Title**: "Characterization of Unruh Channel in the context of Open Quantum Systems"
- **arXiv**: [1603.05450](https://arxiv.org/abs/1603.05450)
- **Channel Characterization**: Uses the Choi matrix approach. Shows the Unruh channel effect on a qubit is similar to vacuum bath interaction.
- **Metrics Studied**: Nonlocality degradation, teleportation fidelity, entanglement loss.

---

## Category D: Hawking Radiation as Quantum Channel {#category-d}

### D1. Bradler, Adami (2015)
- **Title**: "Black holes as bosonic Gaussian channels"
- **Journal**: Physical Review D 92, 025030 (2015)
- **arXiv**: [1405.1097](https://arxiv.org/abs/1405.1097)
- **Channel Type**: Bosonic Gaussian channel (lossy, amplifying, or classical-noise)
- **Channel Definition**: The interaction of a Gaussian quantum state with a Schwarzschild black hole is identified as lying in the non-entanglement-breaking subset of:
  - Lossy channels C(loss)
  - Amplifying channels C(amp)
  - Classical-noise channels B2
- **Channel Parameters**: Depend on black hole mass M and properties of the potential barrier.
- **Key Result**: Classical and quantum capacities calculated. The classification enables capacity estimation even where no tractable expressions exist.
- **Assessment**: The most complete treatment of black holes as quantum channels. However, the channel parameters are black-hole-specific (mass, barrier), not universal gravitational parameters. Not directly giving Sigma = r_s/r.

### D2. Hayden, Preskill (2007)
- **Title**: "Black holes as mirrors: quantum information in random subsystems"
- **Journal**: Journal of High Energy Physics 2007, 120 (2007)
- **Key Result**: After the Page time, quantum information deposited in a black hole is revealed in Hawking radiation very rapidly. The black hole acts as a quantum information mirror.
- **Implication**: The black hole channel has fundamentally different character before and after the Page time -- phase transition in information transmission.

---

## Category E: Relativistic Quantum Channels (Ahmadi, Fuentes, Mann et al.) {#category-e}

### E1. Bruschi, Louko, Martin-Martinez, Dragan, Fuentes (2010)
- **Title**: "Unruh effect in quantum information beyond the single-mode approximation"
- **Journal**: Physical Review A 82, 042332 (2010)
- **arXiv**: [1007.4670](https://arxiv.org/abs/1007.4670)
- **Key Result**: The single-mode approximation commonly used in relativistic QI is not valid for arbitrary states. Provides corrections for both bosonic and fermionic cases.
- **Implication**: Any gravitational channel construction using single-mode approximation must be checked for validity.

### E2. Ahmadi, Bruschi, Fuentes (2014)
- **Title**: "Quantum metrology for relativistic quantum fields"
- **Journal**: Physical Review D 89, 065028 (2014)
- **arXiv**: [1312.5707](https://arxiv.org/abs/1312.5707)
- **Key Result**: Analytical formulas for optimal estimation precision in terms of Bogoliubov coefficients for single-mode and two-mode Gaussian channels in curved spacetime.
- **Channel**: Gaussian channels parametrized by Bogoliubov transformations.

### E3. Bruschi, Datta, Fuentes et al. (2014)
- **Title**: "Quantum estimation of the Schwarzschild spacetime parameters of the Earth"
- **Journal**: Physical Review D 90, 124001 (2014)
- **Key Result**: Uses quantum channel framework to estimate Schwarzschild parameters. The channel is defined through Bogoliubov transformations of field modes.
- **Relevance**: Explicit construction connecting Schwarzschild geometry to quantum channel parameters.

### E4. Kohlrus, Bruschi, Louko, Fuentes (2017)
- **Title**: "Quantum communications and quantum metrology in the spacetime of a rotating planet"
- **Journal**: EPJ Quantum Technology 4, 7 (2017)
- **arXiv**: [1511.04256](https://arxiv.org/abs/1511.04256)
- **Key Result**: Quantum channel for communication in the spacetime of a rotating planet. Shows effects of frame dragging on quantum channel parameters.

---

## Category F: Graviton-Mediated Decoherence Channels {#category-f}

### F1. Anastopoulos, Hu (2013)
- **Title**: "A Master Equation for Gravitational Decoherence: Probing the Textures of Spacetime"
- **Journal**: Classical and Quantum Gravity 30, 165007 (2013)
- **arXiv**: [1305.5231](https://arxiv.org/abs/1305.5231)
- **Channel Type**: Master equation (Lindblad form, can be converted to CPTP map via exponentiation)
- **Key Result**: First-principles derivation of master equation for matter field in linearly perturbed Minkowski spacetime. Predicts decoherence in the ENERGY basis, not position basis (unlike most gravitational decoherence proposals).
- **Decoherence Rate**: Depends on "textures of spacetime" -- temperature and spectral density of the gravitational environment.
- **Assessment**: Rigorous first-principles derivation but gives a DYNAMICAL channel (time-dependent), not a static gravitational transmissivity.

### F2. Anastopoulos, Hu (2017) -- Topical Review
- **Title**: "Gravitational Decoherence"
- **Journal**: Classical and Quantum Gravity 30, 165007 (2017)
- **arXiv**: [1706.05677](https://arxiv.org/abs/1706.05677)
- **Key Review**: Comprehensive review of all gravitational decoherence mechanisms: ABH (Anastopoulos-Blencowe-Hu) from graviton noise, Pikovski from time dilation, Penrose/Diosi from self-gravity.

### F3. Blencowe (2013)
- **Title**: "Effective field theory approach to gravitationally induced decoherence"
- **Journal**: Physical Review Letters 111, 021302 (2013)
- **Key Result**: Decoherence from linearized quantum gravity as effective field theory. Treats gravitons as a thermal bath.

### F4. Oniga, Wang (2017)
- **Title**: "Quantum coherence, radiance, and resistance of gravitational systems"
- **arXiv**: [1701.04122](https://arxiv.org/abs/1701.04122)
- **Key Result**: Framework for open dynamics of quantum particles subject to spacetime fluctuations. Includes spontaneous emission of gravitons and collectively amplified gravitational radiation.

### F5. Hsiang, Cho, Hu (2024)
- **Title**: "Graviton physics: Quantum field theory of gravitons, graviton noise and gravitational decoherence -- a concise tutorial"
- **Journal**: Universe 10, 306 (2024)
- **arXiv**: [2405.11790](https://arxiv.org/abs/2405.11790)
- **Key Value**: Pedagogical review covering graviton physics, graviton noise, and gravitational decoherence. Most comprehensive recent tutorial.
- **Channel Structure**: Master equation approach; the CPTP map is obtained by exponentiating the Lindbladian.

### F6. Lagouvardos, Anastopoulos (2021)
- **Title**: "Gravitational decoherence of photons"
- **arXiv**: [2011.08270](https://arxiv.org/abs/2011.08270)
- **Key Result**: Gravitational decoherence of photons specifically, extending the ABH framework.

### F7. Hsiang, Hu (2025)
- **Title**: "Non-Markovian Quantum Master and Fokker-Planck Equation for Gravitational Systems and Gravitational Decoherence"
- **arXiv**: [2504.11991](https://arxiv.org/html/2504.11991)
- **Key Result**: Non-Markovian generalization of gravitational decoherence master equation. Goes beyond high-temperature Markov approximation.

### F8. Recent (2026)
- **Title**: "Gravitational decoherence and recoherence of a composite particle: the interplay between gravitons and a classical Newtonian potential"
- **arXiv**: [2602.22517](https://arxiv.org/html/2602.22517)
- **Key Result**: Studies interplay between graviton-induced decoherence and classical Newtonian potential effects.

---

## Category G: Channel Capacity in Gravitational Fields {#category-g}

### G1. Landulfo (2016)
- **Title**: "Communication through quantum fields near a black hole"
- **Journal**: Physical Review D 93, 104019 (2016)
- **Key Result**: Studies the quantum channel between two systems communicating in Schwarzschild spacetime via a quantum field. The information carrying capacity includes contributions from direct null geodesics, black-hole-orbiting null geodesics, and timelike contributions.
- **Important Finding**: Non-direct-null and timelike contributions (unique to curved spacetime) can DOMINATE over direct null contributions.

### G2. Kasprzak, Tjoa (2024)
- **Title**: "Transmission of quantum information through quantum fields in curved spacetimes"
- **arXiv**: [2408.00518](https://arxiv.org/abs/2408.00518)
- **Channel Construction**: Uses Unruh-DeWitt detectors coupled to a scalar field in arbitrary curved spacetime.
- **CPTP Map**: Constructed through controlled unitaries coupling detectors to smeared field operators.
- **Quantum Capacity**: Expressed purely in terms of correlation functions and causal propagator.
- **Key Properties**:
  - Manifestly covariant
  - Respects causal structure
  - Independent of background geometry details
  - Can achieve theoretical maximum quantum capacity
- **Assessment**: Most general rigorous construction of a quantum channel in arbitrary curved spacetime. However, probe-dependent (Unruh-DeWitt detector parameters enter).

### G3. Mari, Zippilli, Vitali (2025)
- **Title**: "Can gravity mediate the transmission of quantum information?"
- **Journal**: Physical Review D (2025)
- **arXiv**: [2504.05998](https://arxiv.org/abs/2504.05998)
- **Channel Type**: Quantum thermal attenuator (Gaussian)
- **Channel Definition**:
  ```
  a_out2(omega) = sqrt(eta(omega)) * e^{i*phi(omega)} * a_in1(omega) + sqrt(1-eta(omega)) * a_E(omega)
  ```
  where a_E is a thermal environmental mode.
- **Transmissivity**:
  ```
  eta_opt = (2*lambda^2/gamma^2) / (1 + sqrt(1 + 4*lambda^2/gamma^2) + 2*lambda^2/gamma^2)
  ```
  where lambda = gravitational coupling rate, gamma = mechanical damping.
- **Thermal Noise**: N(omega) = N_T (environmental thermal occupation)
- **Quantum Capacity**: Non-zero when lambda^2 > gamma^2 * N_T * (N_T + 1). Sharp transition from entanglement-breaking to non-classical channel.
- **Assessment**: Explicitly defines gravity as a quantum thermal attenuator. Direct analog to the pure-loss channel. However, it models the gravitational INTERACTION (Newtonian coupling), not the gravitational PROPAGATION (redshift).

### G4. Toccacelo, Andersen, Brask (2025)
- **Title**: "Benchmarks for quantum communication via gravity"
- **Journal**: Physical Review A (2025)
- **arXiv**: [2503.03585](https://arxiv.org/abs/2503.03585)
- **Key Result**: Establishes bounds on quantum state transmission between gravitationally interacting oscillators under different gravity models. Does not require entanglement measurement.

---

## Category H: Modular Flow / Algebraic Channels {#category-h}

### H1. Witten (2022); Chandrasekaran, Longo, Penington, Witten (2023)
- **Title**: "An Algebra of Observables for de Sitter Space"
- **Journal**: Reviews of Modern Physics (2023)
- **arXiv**: [2206.10780](https://arxiv.org/abs/2206.10780)
- **Key Result**: Gravity converts the Type III algebra of QFT in a subregion to a Type II_1 algebra (de Sitter) or Type II_infinity (Schwarzschild). This makes von Neumann entropy and relative entropy well-defined.
- **Channel**: The restriction map (conditional expectation) from the larger algebra to a subalgebra IS a CPTP map. However, it is defined abstractly, not as an explicit Kraus representation.
- **Entropy**: The generalized entropy S_gen = A/(4G_N) + S_out emerges as the entropy of semiclassical states.

### H2. Chandrasekaran, Penington, Witten (2023)
- **Title**: "Generalized entropy for general subregions in quantum gravity"
- **Journal**: Journal of High Energy Physics 2023, 020 (2023)
- **arXiv**: [2306.01837](https://arxiv.org/abs/2306.01837)
- **Key Extension**: Extends H1 to general subregions (not just static patches).

### H3. Kudler-Flam (2024)
- **Title**: "Fluctuation theorems, quantum channels and gravitational algebras"
- **Journal**: Journal of High Energy Physics 2024, 089 (2024)
- **arXiv**: [2408.04219](https://arxiv.org/abs/2408.04219)
- **CRITICAL PAPER for Paper 2**
- **Channel Definition**: Quantum channels in Type II_1 gravitational algebras, represented by SUBFACTORS.
  ```
  rho -> rho' = sum_k E_k rho E_k^dagger,  sum_k E_k^dagger E_k = 1
  ```
- **Crooks-like Fluctuation Theorem in de Sitter**:
  ```
  P_t(-s_bar) = e^{-t*s_bar} * P_t(s_bar)
  ```
  Negative entropy fluctuations are exponentially suppressed.
- **Jarzynski Equality**: <e^{-s_bar}> = 1
- **Entropy Production**: From two-time measurement scheme. For semiclassical states:
  ```
  S_gen = A/(4G_N) + S_out
  ```
- **Key Innovation**: Quantum channels are represented by subfactors using Jones' theory. This connects quantum information to operator algebra in a gravitational context.
- **Assessment**: Most sophisticated algebraic framework connecting quantum channels to gravitational entropy. However, the CPTP map is abstract (conditional expectation on von Neumann algebras), not a Hilbert space operator with explicit Kraus operators.

### H4. Modular Channels, Thermal Filtering (2025)
- **Title**: "Modular Channels, Thermal Filtering and the Spectral Emergence of Spacetime"
- **arXiv**: [2504.20457](https://arxiv.org/abs/2504.20457)
- **CRITICAL PAPER for Paper 2**
- **Channel Definition**:
  ```
  E(rho) = sum_i A_i rho A_i^dagger,  sum_i A_i^dagger A_i = I
  ```
  Arises from partial tracing over inaccessible regions.
- **Singular Value Decomposition**: A = U Sigma V^dagger, where the singular values follow a Gibbs distribution:
  ```
  sigma_i^2 = e^{-beta * k_i} / Z
  ```
  with {k_i} eigenvalues of the modular Hamiltonian K = -ln(rho_in), beta = effective inverse temperature.
- **Unruh Temperature**: For uniformly accelerated observers, beta = 2*pi/a, giving T_Unruh = a/(2*pi).
- **Einstein Equations Derivation**: First law of entanglement delta_S_EE = delta_<K> reinterpreted as Clausius relation delta_Q = T*delta_S. Local validity across Rindler horizons gives:
  ```
  G_mu_nu + Lambda * g_mu_nu = 8*pi*G * T_mu_nu
  ```
- **Spectral Entropy**:
  ```
  S_spec(beta) = -2 * sum_i sigma_i(beta)^2 * log(sigma_i(beta))
  ```
- **Modular Channels Flow Correspondence (MCFC)**: Minimal holographic principle: area of causal screen = storage of filtered quantum information.
- **Assessment**: Provides the most unified framework connecting modular channels to spacetime emergence. The singular value spectrum being Gibbs-weighted gives a NATURAL thermal channel structure. However:
  1. The channel is defined through partial trace (abstract), not as an explicit operator on a specific Hilbert space
  2. The connection Sigma = -ln(-g_00) is not directly established
  3. Derivation of Einstein equations follows Jacobson's thermodynamic route, not a direct channel entropy calculation

### H5. Dorau, Much (2025)
- **Title**: "From Quantum Relative Entropy to the Semiclassical Einstein Equations"
- **Journal**: Physical Review Letters (2025)
- **arXiv**: [2510.24491](https://arxiv.org/abs/2510.24491)
- **Key Result**: Uses modular theory to establish that the Araki-Uhlmann relative entropy between vacuum and coherent excitations on a bifurcate Killing horizon equals the energy flux across the horizon. The semiclassical Einstein equations follow automatically.
- **Relation to N_grav**: The relative entropy DECREASE across a horizon is quantified. For Schwarzschild, the fractional QRE loss is:
  ```
  (D_in - D_out) / D_in = r_s/r  [EXACT]
  ```
- **Assessment**: Provides the most rigorous ALGEBRAIC justification that the "entropy production" of the gravitational channel is r_s/r. But does NOT construct an explicit CPTP map -- it demonstrates consistency with the data processing inequality.

### H6. Jacobson (1995, 2016)
- **Title**: "Thermodynamics of Spacetime: The Einstein Equation of State" (1995); "Entanglement Equilibrium and the Einstein Equation" (2016)
- **arXiv**: [gr-qc/9504004](https://arxiv.org/abs/gr-qc/9504004), [1505.04753](https://arxiv.org/abs/1505.04753)
- **Key Result**: Einstein equations follow from delta_Q = T*delta_S applied to local Rindler horizons. The 2016 version uses entanglement entropy: conformal field theory vacuum entanglement is stationary iff the Einstein equation holds.
- **Relation to Channels**: The Clausius relation can be reinterpreted as a constraint on the entropy production of the modular channel across the horizon.

---

## Category I: Thermal Channels as Gravitational Analogs {#category-i}

### I1. Bianconi (2025)
- **Title**: "Gravity from entropy"
- **Journal**: Physical Review D 111, 066001 (2025)
- **arXiv**: [2408.14391](https://arxiv.org/abs/2408.14391)
- **Key Idea**: The metric of Lorentzian spacetime is treated as a renormalizable effective density matrix. The gravitational action is the quantum relative entropy:
  ```
  S_action = D(rho_spacetime || rho_matter)
  ```
  where rho_spacetime = metric tensor treated as density matrix, and rho_matter = metric induced by matter fields.
- **Modified Einstein Equations**: Reduce to standard Einstein equations with zero cosmological constant in the low-coupling regime. Naturally predicts a small positive cosmological constant.
- **Assessment**: The most radical proposal -- gravity IS a quantum relative entropy. This would make Paper 4's equation Sigma = D(rho_spacetime || rho_matter) literally correct. However, the "density matrix" interpretation of the metric tensor is non-standard and controversial.

### I2. Herrera (2020)
- **Title**: "Landauer Principle and General Relativity"
- **Journal**: Entropy 22, 340 (2020)
- **Key Result**: The Landauer principle applied to gravitational systems shows:
  1. Information has mass; this has consequences for general relativity
  2. Gravitational radiation is irreversible because it conveys information
  3. The Tolman temperature gradient determines the local erasure cost
- **Channel Analog**: Gravity acts as a THERMAL channel with Tolman-temperature-dependent erasure cost.
- **Sigma Derivation**: Tolman-modified erasure cost gives Sigma = r_s/r (one of three routes in Paper 2).

---

## Category J: Quantum Crooks / Fluctuation Theorems in Curved Spacetime {#category-j}

### J1. Basso, Maziero, Celeri (2025)
- **Title**: "Quantum Detailed Fluctuation Theorem in Curved Spacetimes: The Observer Dependent Nature of Entropy Production"
- **Journal**: Physical Review Letters 134, 050406 (2025)
- **arXiv**: [2405.03902](https://arxiv.org/abs/2405.03902)
- **CRITICAL PAPER for Paper 2**
- **Channel**: Two-point measurement (TPM) scheme with time-ordered evolution operator.
- **Entropy Production**:
  ```
  Sigma = W - Delta_F
  ```
  with the detailed fluctuation theorem:
  ```
  P_fwd(W) / P_rev(-W) = e^{beta*(W - Delta_F)}
  ```
- **Curvature Coupling**: The Hamiltonian includes a curvature term:
  ```
  Z(tau) ~ 1 - p^2/(2m^2) + a_i(tau)*x^i + (1/2)*R_{titj}(tau)*x^i*x^j
  ```
  The Riemann tensor component R_{titj} acts as an effective time-dependent potential.
- **Observer Dependence**: "Entropy production is not an invariant quantity defined solely by the system. Rather, it depends on the observer who measures it." This deeply connects the arrow of time with causal structure.
- **Assessment**: EXACTLY what Paper 2 needs -- a rigorous Crooks theorem in curved spacetime with explicit entropy production. The observer dependence of Sigma is compatible with Paper 2's framework (different observers see different tau).

### J2. Basso, Maziero, Celeri (2026)
- **Title**: "Work distribution of quantum fields in static curved spacetimes"
- **Journal**: Physical Review D (2026)
- **Key Result**: Work distributions for quantum scalar fields in static curved spacetimes satisfy both Crooks fluctuation theorem and Jarzynski equality. Extends J1 to field theory.

### J3. Moreira, Celeri (2024)
- **Title**: "Entropy production due to spacetime fluctuations"
- **arXiv**: [2407.21186](https://arxiv.org/abs/2407.21186)
- **Key Result**: Models entropy production in a system interacting with a weak quantum gravitational field (graviton bath). Derives a fluctuation theorem where the first moment is entropy production.
- **Channel**: Open quantum system interacting with graviton bath. Langevin-like equation for geodesic deviation.

---

## Category K: Gravity-Quantum Consistency and Irreversibility {#category-k}

### K1. Galley, Giacomini, Selby (2023)
- **Title**: "Any consistent coupling between classical gravity and quantum matter is fundamentally irreversible"
- **Journal**: Quantum 7, 1142 (2023)
- **arXiv**: [2301.10261](https://arxiv.org/abs/2301.10261)
- **CRITICAL PAPER for Paper 2**
- **No-Go Theorem**: In general probabilistic theories (GPTs), at least one of these must be violated:
  1. Matter is fully non-classical
  2. Matter-gravity interactions are reversible
  3. Matter back-reacts on gravity
- **Implication**: If gravity is classical, the coupling to quantum matter MUST be irreversible (a non-unitary CPTP map). This provides independent justification that the gravitational channel must have Sigma > 0.
- **Assessment**: Proves that gravitational channels are NECESSARILY entropy-producing (irreversible) unless gravity itself is quantum. This is a fundamental consistency requirement that Paper 2's framework satisfies.

---

## Category L: Quantum Metrology in Gravitational Fields {#category-l}

### L1. Iyen et al. (2025)
- **Title**: "Quantum Fisher Information in Curved Spacetime: Dirac Particles in Noisy Channels around a Schwarzschild Black Hole"
- **Journal**: International Journal of Theoretical Physics 64, 205 (2025)
- **arXiv**: [2507.17901](https://arxiv.org/abs/2507.17901)
- **Channel**: Squeezed Generalized Amplitude Damping (SGAD) channel for Dirac particles near Schwarzschild black hole.
- **Kraus Operators**: Standard SGAD Kraus operators parametrized by Hawking temperature.
- **Key Finding**: Squeezing acts as error mitigation; QFI becomes resistant to Hawking temperature variations at strong squeezing (r=1).

### L2. Ahmadi, Bruschi, Sabini, Adesso, Fuentes (2014)
- **Title**: "Relativistic Quantum Metrology: Exploiting relativity to improve quantum measurement technologies"
- **Journal**: Scientific Reports 4, 4996 (2014)
- **Key Result**: Framework for quantum metrology in relativistic settings. Optimal estimation bounds given by quantum Fisher information in curved spacetime.

---

## Evaluation: Best Candidate Channels {#evaluation}

### Ranking of Candidates for Canonical Gravitational Channel N_grav

| Rank | Channel | Sigma | Explicit Kraus? | First-Principles? | Unbounded Sigma? | Universal? |
|------|---------|-------|----------------|-------------------|------------------|-----------|
| **1** | **Pure-loss bosonic (B1)** | -ln(eta) = r_s/r for eta = -g_00 | YES | NO (effective model) | YES | YES (N_B >> 1) |
| 2 | Modular channel (H4) | Gibbs-weighted spectral | Singular values, no explicit K_i | PARTIAL | YES | YES |
| 3 | Algebraic restriction (H1-H2) | Fractional QRE = r_s/r | NO (conditional expectation) | YES | YES | YES |
| 4 | Gravitational thermal attenuator (G3) | -ln(eta) | NO (input-output only) | PARTIAL | YES | PARTIAL |
| 5 | Pikovski dephasing (A1) | H_bin((1+p)/2) | YES | YES | NO (bounded by ln 2) | NO (probe-dependent) |
| 6 | Gravity as dephasing (A4) | Small phase shift | YES | YES | NO (qubit) | YES (universal rate) |
| 7 | Black hole Gaussian (D1) | Capacity-dependent | Classification only | YES | YES | NO (mass-dependent) |
| 8 | Unruh amplitude damping (C2) | Bounded (qubit) | YES | YES | NO (qubit) | NO |

### Assessment Summary

**Winner: Pure-Loss Bosonic Channel with eta = -g_00 and thermal reference (N_B >> 1)**

This is the ONLY candidate that satisfies ALL four requirements simultaneously:
1. Has explicit Kraus operators (beam-splitter model in Fock basis)
2. Gives entropy production Sigma = -ln(-g_00) = r_s/r
3. Reduces to identity for flat spacetime (eta = 1 => Sigma = 0)
4. Works for all static metrics via eta = -g_00(r) (not just Schwarzschild)

**Caveats**:
- The identification eta = -g_00 uses the POWER transmissivity convention (both frequency shift + time dilation), not the AMPLITUDE convention (sqrt(-g_00)). This is a modeling choice.
- The high-temperature limit N_B >> 1 is needed. Motivated by Tolman effect but introduces an approximation.
- The beam-splitter model has limited validity for large redshifts (Leber et al. 2025).

**Runner-Up: Modular Channel (arXiv:2504.20457)**

Most rigorous algebraic framework with:
- Singular values following Gibbs distribution
- Natural connection to Einstein equations
- Information-theoretic derivation of spacetime

But lacks explicit Kraus operators on a Hilbert space and does not directly yield Sigma = r_s/r.

**Key Insight**: There is NO "canonical" gravitational CPTP map that is:
1. Derived from first principles (general relativity + QFT)
2. Has explicit Kraus operators
3. Gives exactly Sigma = -ln(-g_00)
4. Works universally

This is the central **open problem** for Paper 2. The pure-loss bosonic channel is the best available EFFECTIVE model, and the algebraic framework provides the most rigorous STRUCTURAL justification. The gap between them is the gap between Layer 2 (conjecture) and a theorem.

---

## Relevance to Paper 3 (Weak Field / Galactic Scale) {#paper3-relevance}

### Which Results Carry Over?

1. **Pure-loss bosonic channel (B1)**: In the weak field limit r_s/r << 1, the transmissivity eta = 1 - r_s/r + O((r_s/r)^2), and Sigma ~ r_s/r. This gives the Newtonian potential as the leading-order entropy production. At galactic scales (r ~ kpc, M ~ 10^11 M_sun), r_s/r ~ 10^{-6}. The channel is nearly identity with very small entropy production.

2. **Basso-Celeri Crooks theorem (J1)**: The curvature coupling R_{titj} x^i x^j gives an effective time-dependent potential. At galactic scales, R_{titj} ~ M/r^3, providing the tidal field that sources entropy production. The observer dependence is relevant for rotating observers (galactic disk).

3. **Modular channel (H4)**: The spectral entropy S_spec depends on the modular Hamiltonian. For weak fields, the modular Hamiltonian is approximately the boost generator plus corrections. The correction terms could encode the running of G(r) predicted by entropic gravity models.

4. **Galley-Giacomini irreversibility (K1)**: The fundamental irreversibility of classical-gravity-quantum-matter coupling persists at all scales, including galactic. This provides a conceptual foundation for why Sigma > 0 even in weak fields.

5. **Gravitational Landauer (I2)**: The Tolman-modified erasure cost at galactic scales gives:
   ```
   Sigma(r) ~ r_s/r ~ GM/(c^2 * r) ~ 10^{-6}  [at galactic edge]
   ```
   This is too small to explain rotation curves directly. Paper 3 needs the RUNNING of G (logarithmic correction) to get enhanced effects at large r.

### Paper 3 Gap
None of the channels surveyed produce the logarithmic correction needed for rotation curves. The pure-loss channel gives Sigma ~ r_s/r = GM/(c^2 r), which decreases as 1/r. To explain flat rotation curves, one needs Sigma ~ ln(r/r_0) or similar. This requires a fundamentally different channel structure -- possibly the modular channel with scale-dependent spectral weights, or a renormalization group flow of channel parameters.

---

## Open Problems and Conclusions {#conclusions}

### Open Problem 1: First-Principles N_grav
No one has derived an explicit CPTP map from semiclassical gravity (GR + QFT) with Kraus operators giving Sigma = -ln(-g_00). The algebraic framework (Witten/Chandrasekaran) provides the correct STRUCTURE but not explicit operators. The bosonic channel provides explicit operators but not a derivation.

### Open Problem 2: Factor of 2
The amplitude transmissivity sqrt(-g_00) vs power transmissivity -g_00 gives Sigma = r_s/(2r) vs r_s/r. Which is physically correct depends on what observable is being tracked. No consensus exists.

### Open Problem 3: Saturation
No quantum channel with Sigma > 0 is known to saturate the JRSWW bound F = exp(-Sigma/2). This means the exponential metric (Paper 2) is at best a bound, not an exact prediction.

### Open Problem 4: Scale Dependence
All channels surveyed give Sigma ~ r_s/r at large r, which is insufficient for Paper 3 (rotation curves). A channel with logarithmic or sub-linear Sigma(r) behavior is needed.

### Open Problem 5: Reference State
The thermal reference state for the pure-loss channel must be chosen with N_B >> 1 (motivated by Tolman effect). But the Tolman temperature depends on the metric, creating a circularity: one needs the metric to define the channel, but the channel is supposed to determine the metric.

### Overall Conclusion
The state of the art as of 2026 is:
- **There exists an effective CPTP channel** (pure-loss bosonic) that reproduces Sigma_grav = -ln(-g_00) under reasonable physical assumptions.
- **There exists an algebraic framework** (modular channels / gravitational algebras) that rigorously connects QRE to Einstein equations and implies the correct fractional QRE decrease.
- **There does NOT exist** a single canonical gravitational CPTP map derived from first principles with explicit Kraus operators.
- The gap between the effective model and the algebraic framework is the central open problem for the tau framework.

---

## Complete Reference List

### Category A: Gravitational Dephasing
- [A1] Pikovski, Zych, Costa, Brukner, Nature Physics 11, 668 (2015), [arXiv:1311.1095](https://arxiv.org/abs/1311.1095)
- [A2] Bonder, Okon, Sudarsky, Nature Physics 12, 2 (2016), [arXiv:1507.05320](https://arxiv.org/abs/1507.05320)
- [A3] Pang, Chen, PRL 117, 090401 (2016)
- [A4] Balatsky et al., PRA 111, 012411 (2025), [arXiv:2406.03256](https://arxiv.org/abs/2406.03256)

### Category B: Bosonic Gaussian Channels
- [B1] Standard pure-loss channel; see Holevo, "Quantum Systems, Channels, Information" (2019); Wilde, "Quantum Information Theory" (2017)
- [B2] Ivan, Sabapathy, Simon, PRA 84, 042311 (2011), [arXiv:1012.4266](https://arxiv.org/abs/1012.4266)
- [B3] Leber, Rodriguez, Ferreri, Schell, Bruschi (2025), [arXiv:2502.20521](https://arxiv.org/abs/2502.20521)

### Category C: Unruh / Rindler Channels
- [C1] Bradler, Adami, JHEP 2014, 095 (2014), [arXiv:1310.7914](https://arxiv.org/abs/1310.7914)
- [C2] Bradler, Jauregui, Adami (2014), [arXiv:1408.1477](https://arxiv.org/abs/1408.1477)
- [C3] Sinha, Ganapathy (2016), [arXiv:1603.05450](https://arxiv.org/abs/1603.05450)

### Category D: Hawking Radiation as Channel
- [D1] Bradler, Adami, PRD 92, 025030 (2015), [arXiv:1405.1097](https://arxiv.org/abs/1405.1097)
- [D2] Hayden, Preskill, JHEP 2007, 120 (2007)

### Category E: Relativistic Quantum Channels
- [E1] Bruschi, Louko, Martin-Martinez, Dragan, Fuentes, PRA 82, 042332 (2010), [arXiv:1007.4670](https://arxiv.org/abs/1007.4670)
- [E2] Ahmadi, Bruschi, Fuentes, PRD 89, 065028 (2014), [arXiv:1312.5707](https://arxiv.org/abs/1312.5707)
- [E3] Bruschi, Datta, Fuentes et al., PRD 90, 124001 (2014)
- [E4] Kohlrus, Bruschi, Louko, Fuentes, EPJQT 4, 7 (2017), [arXiv:1511.04256](https://arxiv.org/abs/1511.04256)

### Category F: Graviton-Mediated Decoherence
- [F1] Anastopoulos, Hu, CQG 30, 165007 (2013), [arXiv:1305.5231](https://arxiv.org/abs/1305.5231)
- [F2] Anastopoulos, Hu, Topical Review (2017), [arXiv:1706.05677](https://arxiv.org/abs/1706.05677)
- [F3] Blencowe, PRL 111, 021302 (2013)
- [F4] Oniga, Wang (2017), [arXiv:1701.04122](https://arxiv.org/abs/1701.04122)
- [F5] Hsiang, Cho, Hu, Universe 10, 306 (2024), [arXiv:2405.11790](https://arxiv.org/abs/2405.11790)
- [F6] Lagouvardos, Anastopoulos (2021), [arXiv:2011.08270](https://arxiv.org/abs/2011.08270)
- [F7] Hsiang, Hu (2025), [arXiv:2504.11991](https://arxiv.org/html/2504.11991)
- [F8] Gravitational decoherence and recoherence (2026), [arXiv:2602.22517](https://arxiv.org/html/2602.22517)

### Category G: Channel Capacity in Gravitational Fields
- [G1] Landulfo, PRD 93, 104019 (2016)
- [G2] Kasprzak, Tjoa (2024), [arXiv:2408.00518](https://arxiv.org/abs/2408.00518)
- [G3] Mari, Zippilli, Vitali, PRD (2025), [arXiv:2504.05998](https://arxiv.org/abs/2504.05998)
- [G4] Toccacelo, Andersen, Brask, PRA (2025), [arXiv:2503.03585](https://arxiv.org/abs/2503.03585)

### Category H: Modular Flow / Algebraic Channels
- [H1] Chandrasekaran, Longo, Penington, Witten (2023), [arXiv:2206.10780](https://arxiv.org/abs/2206.10780)
- [H2] Chandrasekaran, Penington, Witten, JHEP (2023), [arXiv:2306.01837](https://arxiv.org/abs/2306.01837)
- [H3] Kudler-Flam, JHEP (2024), [arXiv:2408.04219](https://arxiv.org/abs/2408.04219)
- [H4] Modular Channels, Thermal Filtering (2025), [arXiv:2504.20457](https://arxiv.org/abs/2504.20457)
- [H5] Dorau, Much, PRL (2025), [arXiv:2510.24491](https://arxiv.org/abs/2510.24491)
- [H6] Jacobson, PRL (1995), [arXiv:gr-qc/9504004](https://arxiv.org/abs/gr-qc/9504004); PRL (2016), [arXiv:1505.04753](https://arxiv.org/abs/1505.04753)

### Category I: Thermal Channels as Gravitational Analogs
- [I1] Bianconi, PRD 111, 066001 (2025), [arXiv:2408.14391](https://arxiv.org/abs/2408.14391)
- [I2] Herrera, Entropy 22, 340 (2020)

### Category J: Quantum Crooks in Curved Spacetime
- [J1] Basso, Maziero, Celeri, PRL 134, 050406 (2025), [arXiv:2405.03902](https://arxiv.org/abs/2405.03902)
- [J2] Basso, Maziero, Celeri, PRD (2026)
- [J3] Moreira, Celeri (2024), [arXiv:2407.21186](https://arxiv.org/abs/2407.21186)

### Category K: Gravity-Quantum Consistency
- [K1] Galley, Giacomini, Selby, Quantum 7, 1142 (2023), [arXiv:2301.10261](https://arxiv.org/abs/2301.10261)

### Category L: Quantum Metrology
- [L1] Iyen et al., IJTP 64, 205 (2025), [arXiv:2507.17901](https://arxiv.org/abs/2507.17901)
- [L2] Ahmadi, Bruschi, Sabini, Adesso, Fuentes, Scientific Reports 4, 4996 (2014)

### Additional Key References
- [R1] Singh et al., "Realizing the Petz Recovery Map on an NMR Quantum Processor" (2025), [arXiv:2508.08998](https://arxiv.org/abs/2508.08998)
- [R2] Wilde, "Entropy of a quantum channel", PRR 3, 023096 (2021)
- [R3] "Quantum Relative Entropy of Schwarzschild Black Hole and the Area Law", Entropy 27, 266 (2025)

---

**Total papers surveyed**: 40+
**Papers with explicit CPTP maps**: ~15
**Papers with gravitational Kraus operators**: ~8
**Papers giving Sigma = r_s/r or equivalent**: 5 (via different routes)

**Last updated**: 2026-03-10
