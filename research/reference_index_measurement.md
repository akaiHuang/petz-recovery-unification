# Reference Index: Measurement Problem & tau Framework

> Compiled 2026-03-09. Papers connecting Petz recovery maps, quantum channels, recovery fidelity, retrodiction, and entropy production to the measurement problem, wave function collapse, and the quantum-to-classical transition.

---

## Category 1: Petz Recovery — Theory & Foundations

### [Parzygnat & Buscemi 2023] Axioms for retrodiction: achieving time-reversal symmetry with a prior
- arXiv: 2210.13531
- Journal: Quantum 7, 1013 (2023)
- Key result: Category-theoretic framework for retrodiction. Among all retrodiction families tested (including rotated and averaged Petz maps), the original Petz recovery map is the *only* one that satisfies the retrodiction functor axioms (compositionality). This establishes the Petz map as the unique quantum analogue of Bayes' theorem with full structural consistency.
- Relevance to tau: The Petz map as retrodiction functor is the mathematical backbone of the tau framework. tau = 1 - F quantifies how far the Petz map fails to perfectly retrodict, i.e., the irreversibility cost.

### [Bai, Buscemi & Scarani 2025] Quantum Bayes' rule and Petz transpose map from the minimum change principle
- arXiv: 2410.00319
- Journal: Physical Review Letters 135, 090203 (2025)
- Key result: Derives the Petz transpose map from a higher principle -- the quantum minimum change principle. When minimizing fidelity-based distance between forward and reverse processes, the solution is unique and recovers the Petz map. This is the first derivation of the Petz map from a variational principle rather than algebraic construction.
- Relevance to tau: Provides foundational justification for why the Petz map (and hence tau = 1 - F_Petz) is the "right" measure of irreversibility -- it arises from a principle of minimal informational disturbance.

### [Junge, Renner, Sutter, Wilde & Winter 2018] Universal recovery maps and approximate sufficiency of quantum relative entropy
- arXiv: 1509.07127
- Journal: Annales Henri Poincare 19(10), 2955-2978 (2018)
- Key result: Strengthened data processing inequality: D(rho || sigma) - D(N(rho) || N(sigma)) >= -log F(rho, R_sigma o N(rho)), where R_sigma is a universal recovery map depending only on sigma and N. The recovery map is explicit and channel-dependent but state-independent. This is the JRSWW bound that underpins the equivalence chain.
- Relevance to tau: The JRSWW bound is the core inequality in Paper 1. It gives F >= exp(-Sigma/2), which translates to tau <= 1 - exp(-Sigma/2). This is the quantitative link between entropy production and recoverability.

### [Surace & Scandi 2023] State retrieval beyond Bayes' retrodiction
- arXiv: 2201.09899
- Journal: Quantum 7, 990 (2023)
- Key result: Shows that Bayesian retrodiction is just one member of a whole class of reverse transformations. The class can be optimized to achieve better state retrieval than Bayes' rule alone. The Petz recovery map appears as a special case within this broader class, confirming its interpretation as the quantum Bayes inverse. Adding a single extra axiom isolates the Bayesian reverse uniquely.
- Relevance to tau: Establishes that while tau = 1 - F_Petz uses a specific (Bayesian) reverse, this choice is canonical and robust -- other choices in the class can only do better, making the Petz-based tau an upper bound on achievable irreversibility.

### [Sohail, Pandey, Singh & Das 2024] Fundamental limitations on the recoverability of quantum processes
- arXiv: 2403.12947
- Journal: Annales Henri Poincare (2025)
- Key result: Extends the data processing inequality from states to quantum channels under superchannels. Identifies a class of superchannels (analogue of subunital channels) under which channel entropy is non-decreasing. Strengthens the quantum DPI at the process level, establishing fundamental limits on how well physical transformations on channels can be reversed.
- Relevance to tau: Lifts the tau framework from state-level to channel-level. If tau quantifies state irreversibility via Petz, this work quantifies process irreversibility, potentially enabling a "second-order tau" for quantum error correction protocols themselves.

### [Song, Kwon & Scarani 2025] Exact and approximate conditions of tabletop reversibility: when is Petz recovery cost-free?
- arXiv: 2510.26895
- Journal: arXiv preprint (2025)
- Key result: Investigates when the Petz recovery map can be implemented using only the same hardware as the forward dynamics ("tabletop reversibility"). Exact TTR requires time-sensitive ancilla control; approximate TTR relaxes this. Derives Lindbladian TTR conditions under random-time collision models.
- Relevance to tau: Addresses the practical question: when can tau = 0 be achieved "for free"? The answer reveals that even when perfect recovery is mathematically possible (tau = 0), the physical implementation may carry costs -- connecting tau to thermodynamic work requirements.

### [Bhattacharya 2023] Approximate recoverability and the quantum data processing inequality
- arXiv: 2309.02074
- Journal: arXiv preprint (2023, revised 2025)
- Key result: Disproves a conjecture by Seshadreesan et al. on approximate recoverability. Establishes universal approximate recoverability inequalities using the Petz recovery map for sandwiched Renyi relative entropies at parameter t=2. Derives new convexity theorems for parametrized relative entropy and fidelity.
- Relevance to tau: Refines the mathematical foundations of the Petz-based recovery bounds, potentially tightening the tau <= 1 - exp(-Sigma/2) inequality for specific divergence families.

---

## Category 2: Petz Recovery — Experimental Implementations

### [Png & Scarani 2025] Petz recovery maps of single-qubit decoherence channels in an ion trap quantum processor
- arXiv: 2504.20399
- Journal: Physical Review A 112, 022613 (2025)
- Key result: First ion-trap implementation proposal for exact state-specific Petz maps. Circuit constructions require at most 1 ancilla qubit and 3 CNOT gates for Kraus rank 2 channels. Explicit circuits provided for depolarizing, dephasing, and amplitude damping channels. Precision requirements characterized for sub-0.01 recovery error.
- Relevance to tau: Enables direct experimental measurement of tau = 1 - F_recovered for standard decoherence channels. The recovery fidelity data can test F >= exp(-Sigma/2) on real hardware.

### [Singh, Sahani, Jagadish, Lautenbacher, Bernardes & Dorai 2025] Realizing the Petz Recovery Map on an NMR Quantum Processor
- arXiv: 2508.08998
- Journal: arXiv preprint (2025)
- Key result: First experimental realization of the Petz recovery map on an NMR quantum processor using the duality quantum computing (DQC) algorithm. Demonstrated for phase damping and amplitude damping channels; recovered states closely match theoretical predictions. Validates feasibility of Petz-based recovery in current quantum platforms.
- Relevance to tau: Provides experimental fidelity data that can be directly compared against the tau bound. Combined with Png's ion-trap proposal, these two platforms offer independent tests of F >= exp(-Sigma/2).

---

## Category 3: Retrodiction, Bayesian Inference & Time Symmetry

### [Parzygnat 2024] Reversing information flow: retrodiction in semicartesian categories
- arXiv: 2401.17447
- Journal: arXiv preprint (2024)
- Key result: Extends retrodiction to arbitrary semicartesian categories, going beyond classical and quantum probability. Incorporates Bayesian inference and Jeffrey's probability kinematics into a general categorical framework. Shows retrodiction generalizes time-reversal even when dynamics is not invertible.
- Relevance to tau: Provides the most abstract mathematical setting for retrodiction. The tau framework is a specific instance of this general theory -- this paper guarantees the conceptual coherence of using Petz-based retrodiction in any physical theory satisfying the semicartesian axioms.

### [Kuang, Torii & Buscemi 2025] Quantum measurement retrodiction and entropic uncertainty relations
- arXiv: 2511.20281
- Journal: arXiv preprint (2025)
- Key result: For quantum-to-classical measurement channels, ALL standard quantum divergences select the same retrodictive update -- a unique, divergence-independent quantum Bayesian inverse for any POVM and prior state. Introduces "mutual retrodictability" and derives two new entropic uncertainty relations that are consistently tighter than existing ones.
- Relevance to tau: The divergence-independence result is remarkable: it means tau for measurement channels is robust regardless of which divergence defines the entropy production Sigma. This strengthens the universality claim of the tau framework.

### [Liu, Bai & Scarani 2025] The state of a quantum system is not a complete description for retrodiction
- arXiv: 2502.10030
- Journal: Physical Review Letters 136, 060203 (2026)
- Key result: Mixed quantum states that give identical predictions for future measurements can give DIFFERENT retrodictions depending on their interpretation (proper mixture, improper mixture, or unspecified ignorance). This is a purely quantum feature -- classically, identical probability distributions always give identical retrodictions. Proves necessary and sufficient conditions for belief equivalence in retrodiction.
- Relevance to tau: Critical for the tau framework: tau = 1 - F depends not just on the state but on the *belief* about the state. For improper mixtures (from entanglement), retrodiction carries fundamentally different information than for proper mixtures. This connects to why decoherence (which converts improper to effectively proper mixtures) is irreversible at the retrodictive level.

### [Pinske & Molmer 2025] Retrodiction of measurement outcomes on a single quantum system reveals entanglement with its environment
- arXiv: 2504.16170
- Journal: Physical Review A 112, 032208 (2025)
- Key result: Retrodiction using both prior and posterior knowledge yields conditional probabilities that can witness whether a single quantum system is entangled with its environment. The degree of certainty in retroactively predicting multiple measurement outcomes quantifies both the existence and nature of environmental entanglement -- something the density matrix alone cannot distinguish.
- Relevance to tau: Directly connects retrodiction (the operational content of tau) to the decoherence/entanglement structure. A system with tau > 0 has leaked information to the environment; this paper shows that retrodiction can diagnose and quantify exactly this leakage from single-system measurements.

### [Liu, Bai & Scarani 2025] Unifying Quantum Smoothing Theories with Extended Retrodiction
- arXiv: 2510.08447
- Journal: arXiv preprint (2025)
- Key result: Develops a unified retrodictive framework showing all existing quantum smoothing approaches are special cases corresponding to different extended prior beliefs. Upper and lower bounds on smoothed state entropy are achieved by Petz-Fuchs and CLHS smoothed states respectively. Clarifies the role of priors in quantum retrodiction.
- Relevance to tau: Shows the Petz map (hence tau) naturally appears as the upper bound in the smoothing hierarchy. Quantum smoothing = optimal retrodiction given both past and future data, which is precisely the quantum eraser scenario that motivated the tau framework.

### [Bai 2024] Bayesian retrodiction of quantum supermaps
- arXiv: 2408.07885
- Journal: arXiv preprint (2024, revised 2025)
- Key result: Extends quantum Bayes' rule (Petz map) from channels to higher-order operations (supermaps). For specific families of initial beliefs, retrodiction supermaps can be constructed that implement the belief update via a "reverse" quantum supermap. Applications demonstrated for improved error correction in quantum cloud computing.
- Relevance to tau: Lifts tau to the supermap level. If tau quantifies channel irreversibility, this work enables quantifying the irreversibility of error correction protocols themselves -- a "meta-tau" for quantum computing architectures.

### [Aw & Liu 2025] Quantifying Irreversibility via Bayesian Subjectivity for Classical & Quantum Linear Maps
- arXiv: 2503.12112
- Journal: Physical Review E 112, 054123 (2025)
- Key result: Proposes quantifying irreversibility by the "Bayesian subjectivity" -- the sensitivity of retrodiction to one's choice of prior. More irreversible processes yield retrodictions that depend more heavily on the prior rather than on the data. Analytical results confirm physical intuitions at extremes of reversibility and irreversibility.
- Relevance to tau: Provides an alternative but complementary measure to tau. While tau = 1 - F_Petz measures recovery failure, Bayesian subjectivity measures prior-dependence. Both vanish for reversible processes and both increase with entropy production. The relationship between these two measures may yield new bounds.

---

## Category 4: Entropy Production, Irreversibility & Measurement

### [Bai, Buscemi & Scarani 2024] Fully quantum stochastic entropy production
- arXiv: 2412.12489
- Journal: arXiv preprint (2024)
- Key result: Defines entropy production for arbitrary quantum processes using Bayesian tools. The entropy production operator is Hermitian, has non-negative average (Second Law), satisfies both Jarzynski equality and Crooks fluctuation theorem, and recovers classical expressions. Identifies features where the quantum version fundamentally differs from straightforward translation of classical stochastic thermodynamics.
- Relevance to tau: This is the natural partner to the tau framework. Sigma in tau <= 1 - exp(-Sigma/2) needs a fully quantum definition of entropy production; this paper provides exactly that. The Crooks theorem connection means tau inherits the fluctuation theorem structure: P(tau < 0) / P(tau > 0) = exp(-|Sigma|).

### [Clarke & Ford 2023] Stochastic entropy production associated with quantum measurement in a framework of Markovian quantum state diffusion
- arXiv: 2301.08197
- Journal: Entropy 26, 1024 (2024)
- Key result: Models continuous quantum measurement as quantum state diffusion. Stochastic entropy production characterizes how a system asymptotically approaches an eigenstate during measurement. For simultaneous measurement of non-commuting observables, the system reaches a stationary distribution with zero further entropy production rather than collapsing to an eigenstate. Different measurement configurations produce distinct entropy signatures.
- Relevance to tau: Provides the dynamical picture of how tau evolves during continuous measurement. As the system collapses (tau -> 1 for the measured observable), entropy is produced; when non-commuting observables compete, tau reaches a nonzero steady state -- the system cannot fully collapse to either eigenstate.

### [Walls, Bloss & Ford 2025] Characterising quantum measurement through environmental stochastic entropy production in a two spin-1/2 system
- arXiv: 2504.10172
- Journal: Physical Review A (2025)
- Key result: Extends Clarke & Ford's framework to a two-spin system. The mean asymptotic rates of environmental stochastic entropy production during collapse depend on which eigenstate is selected AND on the initial state. Different measurement configurations (single-particle vs total spin) produce qualitatively different entropy signatures.
- Relevance to tau: The eigenstate-dependent entropy production rates mean tau is not just a property of the channel but couples to the measurement outcome. This is precisely what tau_signed captures: different outcomes produce different amounts of irreversibility.

### [Xue et al. 2024] Evidence of genuine quantum effects in nonequilibrium entropy production
- arXiv: 2402.06858
- Journal: Physical Review A 110, 042204 (2024)
- Key result: Experimentally demonstrates decomposition of entropy production into population-related (classical) and coherence-related (genuinely quantum) components using optical systems. The coherence contribution has no classical counterpart. Irreversibility at the quantum level can be reduced through proper management of both contributions.
- Relevance to tau: Shows that Sigma in the tau bound has two distinct sources: classical population imbalance and quantum coherence. This suggests tau itself might be decomposable: tau = tau_pop + tau_coh, with the coherence part being the genuinely quantum contribution to the arrow of time.

### [Ferri-Cortes et al. 2023] Conditional fluctuation theorems and entropy production for monitored quantum systems under imperfect detection
- arXiv: 2308.08491
- Journal: Physical Review Research (2024)
- Key result: Establishes universal fluctuation relations for open quantum systems monitored with imperfect detectors. Derives connection between thermodynamic entropy production and information-theoretical irreversibility at the trajectory level. Provides a practical estimator using incomplete measurement records to lower-bound entropy production.
- Relevance to tau: Addresses the realistic scenario where measurement of Sigma (and hence tau) is imperfect. The trajectory-level fluctuation relation means tau can be estimated even from incomplete data, making the framework experimentally accessible.

### [Reddy 2026] Information-Induced Quantum Measurement: Entropy Production and the Dynamical Origin of Wavefunction Collapse
- DOI: 10.21203/rs.3.rs-8810607/v1
- Journal: Research Square preprint (2026)
- Key result: Proposes a dynamical theory where wavefunction collapse emerges from irreversible information transfer and entropy production in system-detector-environment interactions. Derives stochastic nonlinear evolution from microscopic Hamiltonians. The collapse rate is uniquely determined by environmental entropy production. Claims to establish the Born rule dynamically and resolve foundational paradoxes.
- Relevance to tau: Most directly relevant to Paper 1b. If collapse rate ~ Sigma (as this paper claims), then tau = 1 - F provides the quantitative bound on how "collapsed" a system can be. The connection tau -> 1 iff Sigma -> infinity matches this paper's prediction that large entropy production ensures definite outcomes.

---

## Category 5: Gravitational Decoherence & Collapse Models

### [Galley, Giacomini & Selby 2023] Any consistent coupling between classical gravity and quantum matter is fundamentally irreversible
- arXiv: 2301.10261
- Journal: Quantum 7, 1142 (2023)
- Key result: No-go theorem using General Probabilistic Theories: if gravity remains classical and quantum matter back-reacts on it, the interaction must be irreversible. Specifically, if (i) matter is fully non-classical, (ii) matter back-reacts on gravity, and (iii) gravity is classical, then the coupling cannot be reversible. Either gravity is non-classical or the coupling is irreversible.
- Relevance to tau: This is the strongest theoretical argument for tau > 0 in gravitational contexts. If classical-gravity-quantum-matter coupling is fundamentally irreversible, then tau_grav > 0 is not approximate but *necessary*. This provides the foundational justification for Paper 2's Sigma_grav = r_s/r: the irreversibility is built into the structure of the coupling.

### [Pikovski, Zych, Costa & Brukner 2015] Universal decoherence due to gravitational time dilation
- arXiv: 1311.1095
- Journal: Nature Physics 11, 668-672 (2015)
- Key result: Gravitational time dilation couples internal degrees of freedom to center-of-mass position, causing decoherence of spatial superpositions without any external environment. The weak time dilation on Earth is already sufficient to decohere micron-scale objects. This provides a mechanism for gravity-induced classicality testable in matter-wave experiments.
- Relevance to tau: The Pikovski channel is the explicit CPTP map used in Paper 1b to connect tau to gravitational decoherence. tau_int = 2*tau_sp*(1 - tau_sp) connects internal decoherence to spatial decoherence, with Sigma = S_entanglement. This is the bridge between the abstract tau framework and concrete gravitational physics.

### [Danielson, Satishchandran & Wald 2022] Black Holes Decohere Quantum Superpositions
- arXiv: 2205.06279
- Journal: International Journal of Modern Physics D 31(14), 2241003 (2022)
- Key result: A nearby black hole destroys quantum superpositions by harvesting "which path" information through soft graviton radiation into the black hole. The decoherence rate depends on the superposition size and distance to the black hole. Won third prize in the 2022 Gravity Research Foundation Essay Competition.
- Relevance to tau: Black hole horizons act as perfect information sinks: tau -> 1 at the horizon (complete information loss, perfect decoherence). This is the strong-field limit of Paper 2's Sigma_grav = r_s/r, where r -> r_s gives Sigma -> 1 and tau -> 1 - exp(-1/2) ~ 0.39, approaching maximal decoherence.

### [Danielson, Satishchandran & Wald 2023] Killing Horizons Decohere Quantum Superpositions
- arXiv: 2301.00026
- Journal: Physical Review D 108, 025007 (2023)
- Key result: Generalizes the black hole result to ANY spacetime with a Killing horizon (Rindler, de Sitter, cosmological). Soft gravitons/photons propagate through the horizon carrying "which path" information, decohering superpositions at a steady rate. The decoherence is frame-dependent (Rindler vs inertial descriptions give equivalent physics) and occurs in finite time.
- Relevance to tau: Universal gravitational decoherence = universal tau_grav > 0 near horizons. This supports Paper 2's thesis that Sigma_grav is a geometric property of spacetime rather than a property of specific matter. The Rindler/inertial equivalence maps to the observer-dependence of tau in accelerated frames.

### [Artini, Lo Monaco, Donadi & Paternostro 2025] Non-equilibrium thermodynamics of gravitational objective-collapse models
- arXiv: 2502.03173
- Journal: Physical Review Research 7, 043017 (2025)
- Key result: Analyzes entropy production in the Diosi-Penrose collapse model. The original DP model produces unbounded heating, consistent with thermodynamics only under infinite-temperature noise assumption. The dissipative extension achieves physically consistent thermalization at low dissipation strength. Phase-space dynamics and entropy production rates fully characterized.
- Relevance to tau: Quantifies Sigma for objective collapse models. The unbounded heating of the original DP model means tau -> 1 (complete irreversibility), which is thermodynamically problematic. The dissipative extension that thermalizes corresponds to finite Sigma and hence finite tau -- more physically reasonable. This constrains which collapse models are compatible with the tau framework.

### [Dorau & Much 2025] From Quantum Relative Entropy to the Semiclassical Einstein Equations
- arXiv: 2510.24491
- Journal: Physical Review Letters 136, 091602 (2026)
- Key result: Using modular theory, proves that relative entropy between vacuum and coherent excitations on a bifurcate Killing horizon equals the energy flux across the horizon. Combined with Bekenstein-Hawking entropy-area formula, this directly yields the semiclassical Einstein equations. Generalizes Jacobson's thermodynamic derivation using quantum relative entropy instead of classical thermodynamic entropy.
- Relevance to tau: THE key paper for Paper 2. Sigma_grav = D(rho || sigma) directly implies Einstein equations. This closes the logical circle: entropy production (Sigma) -> curvature (G_ab) -> decoherence (tau > 0). The quantum relative entropy is the same QRE that appears in the JRSWW bound, making the connection between Papers 1 and 2 rigorous.

---

## Category 6: Time Symmetry, Retrocausality & Interpretations

### [Frank 2025] Time Symmetry, Retrocausality and the Emergent Arrow of Time: the Quantum Time-Symmetric Interpretation (QTSI)
- arXiv: 2508.19301
- Journal: arXiv preprint (2025)
- Key result: Proposes QTSI where isolated quantum systems have forward- and backward-propagating amplitudes. A "retrocausal coherence time" parameter governs their mixing. As systems couple to amplifying environments, backward components are dynamically suppressed, producing classical causality and effective collapse. Consistent with weak-measurement and quantum eraser experiments. Predicts testable signatures in temporal echoes and chaotic cavities.
- Relevance to tau: QTSI's retrocausal coherence parameter maps directly to tau. tau = 0 (isolated system) = both temporal components present = no arrow of time. tau > 0 (coupled to environment) = backward component suppressed = emergent time arrow. The "amplifying environment" that kills the backward component is precisely the mechanism that makes Sigma > 0.

### [Inage 2025] A Time-Symmetric Formulation of Quantum Measurement: Reinterpreting the Arrow of Time as Information Flow
- arXiv: 2511.22191
- Journal: arXiv preprint (2025)
- Key result: Models measurement as bidirectional informational update between forward-evolving state and backward-propagating effect. Preserves complete positivity, normalization, and no-signalling. In the classical limit, reduces to Kalman filter and RTS smoother. The temporal asymmetry of measurement arises from "the one-sided way in which measurement outcomes are incorporated into our description" rather than fundamental irreversibility.
- Relevance to tau: Provides the operational interpretation of tau as an information-flow asymmetry. tau = 0 means bidirectional information flow is balanced (time-symmetric). tau > 0 means the forward direction dominates -- not because physics is irreversible, but because our informational conditioning is one-sided. This aligns perfectly with the Bayesian interpretation of tau.

---

## Category 7: Quantum Smoothing & State Estimation

### [Liu, Bai & Scarani 2025] A Retrodictive Approach to Quantum State Smoothing
- arXiv: 2501.15986
- Journal: Physical Review A (2025)
- Key result: Applies the Bayesian retrodiction framework to quantum state smoothing, where both past preparation and future measurement data are used to estimate the present state. Connects quantum smoothing to the Petz map through the retrodictive framework.
- Relevance to tau: Quantum smoothing is the operational task behind the quantum eraser: using future measurement results to refine present state estimates. The smoothing fidelity improvement over filtering alone quantifies how much "future information" helps, which is precisely 1 - tau.

---

## Summary Statistics

- **Total papers indexed:** 27
- **Date range:** 2015-2026
- **Journals represented:** PRL (4), PRA (3), PRD (2), PRR (1), PRE (1), Nature Physics (1), IJMPD (1), Quantum (3), AHP (2), Entropy (1), arXiv preprints (8)
- **Key author clusters:**
  - Buscemi-Parzygnat-Bai-Scarani axis (retrodiction + Bayesian inference): 9 papers
  - Danielson-Satishchandran-Wald (gravitational decoherence): 2 papers
  - Clarke-Ford-Walls (stochastic entropy & measurement): 2 papers
  - Wilde-Junge-Renner-Sutter-Winter (universal recovery): 2 papers

## Cross-Reference with tau Framework Papers

| tau Paper | Most relevant references from this index |
|-----------|------------------------------------------|
| Paper 1 (equivalence chain) | JRSWW 2018, Parzygnat-Buscemi 2023, Bai-Buscemi-Scarani 2025 (PRL), Surace-Scandi 2023 |
| Paper 1b (collapse = Petz failure) | Galley et al. 2023, Pikovski et al. 2015, Artini et al. 2025, Reddy 2026, Clarke-Ford 2023 |
| Paper 2 (exponential metric) | Dorau-Much 2025, Danielson et al. 2022/2023, Galley et al. 2023 |
| Paper 4 (unification) | Dorau-Much 2025, Bai et al. 2024 (fully quantum Sigma), Inage 2025, Frank 2025 |

---

## Key Open Questions Identified

1. **Divergence-independence of tau for measurement channels** (Kuang et al. 2025): Does this extend beyond quantum-to-classical channels to general quantum channels?

2. **Decomposition of Sigma into classical + quantum parts** (Xue et al. 2024): Can tau be similarly decomposed as tau = tau_pop + tau_coh?

3. **Superchannel-level tau** (Sohail et al. 2024, Bai 2024): What does the second-order DPI imply for the tau of error correction protocols?

4. **Belief-dependence of retrodiction** (Liu et al. 2026 PRL): How does the improper/proper mixture distinction affect tau in gravitational decoherence scenarios?

5. **Cost of Petz recovery** (Song et al. 2025): When is tau = 0 achievable without thermodynamic cost? Connection to Landauer's principle?
