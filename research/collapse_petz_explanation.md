# Quantum Collapse as Petz Recovery Failure: The Full Argument

**Author**: Sheng-Kai Huang
**Date**: 2026-03-20
**Status**: Comprehensive research document
**Purpose**: Address head-on the question "Can we explain the REAL reason for collapse? Or do we only have kinematics?"

---

## The Honest Answer, Up Front

We have **more than kinematics but less than a complete ontology**.

Specifically:
- We have a **quantitative, parameter-free prediction** for when, how fast, and how completely collapse occurs (Theorems 1--3 of Paper 1b).
- We have a **mechanism**: the tensor-product structure of quantum mechanics + gravitational time dilation + the data-processing inequality, producing exponential suppression of recovery fidelity.
- We have a **conceptual clarification**: collapse is not a new process --- it is the failure of the unique Bayesian retrodiction functor (the Petz map) when the subsystem is open.
- We **do not** have: an explanation for the Born rule, a resolution of whether branches "cease to exist," or an answer to the hard problem of consciousness.

What follows is the full argument, with honest flags for where we have theorems and where we have interpretations.

---

## 1. What Does "Collapse" Mean Operationally?

### 1.1 The textbook story

Before measurement: system in superposition

    |psi> = alpha |0> + beta |1>

After measurement: system found in |0> or |1>, with probabilities |alpha|^2, |beta|^2.

The density matrix transitions from a pure state to an apparent mixture:

    rho_pure = |psi><psi| = |alpha|^2 |0><0| + alpha beta* |0><1| + alpha* beta |1><0| + |beta|^2 |1><1|

    --> rho_mixed = |alpha|^2 |0><0| + |beta|^2 |1><1|

The off-diagonal terms (coherences) vanish.

### 1.2 The three sub-problems (Schlosshauer 2007)

The measurement problem is actually three problems:

**(i) Preferred basis**: Why does the system collapse into position eigenstates rather than momentum eigenstates or some arbitrary basis?

**(ii) Outcomes**: Why does ONE outcome occur, rather than the superposition persisting?

**(iii) Timing**: When exactly does the transition happen?

Standard decoherence (Zurek 2003, Joos et al. 2003) answers (i) through einselection and gives a practical answer to (iii) via the decoherence timescale. It leaves (ii) open --- the diagonal density matrix is still a *proper* mixture "for all practical purposes" (FAPP), but other branches still exist in the full Hilbert space.

### 1.3 What people really want to know

The persistent dissatisfaction with decoherence is this: *explaining why off-diagonal terms become small is not the same as explaining why you experience one definite outcome*. Penrose, among others, insists that something *physical* must happen --- a real, objective, non-unitary process.

We address this dissatisfaction directly.

---

## 2. The Petz Recovery Answer

### 2.1 Measurement as a quantum channel

A measurement is a physical interaction. The system S couples to an apparatus/environment E:

    rho_S  -->  N(rho_S) = Tr_E [ U (rho_S tensor rho_E) U^dagger ]

This is a quantum channel N: the completely positive, trace-preserving map from the system's initial state to its reduced state after interaction. This is not an interpretation --- it is the Stinespring representation theorem, a mathematical fact.

### 2.2 The Petz map: the unique way to undo a channel

Given a channel N and a reference state sigma, the Petz recovery map is:

    R_{sigma,N}(rho) = sigma^{1/2} N^dagger( N(sigma)^{-1/2} rho N(sigma)^{-1/2} ) sigma^{1/2}

Parzygnat and Buscemi (2023, Quantum 7, 1013) proved that this is the **unique** retrodiction functor satisfying:
- Bayesian consistency (reduces to Bayes' theorem classically)
- Normalization (outputs a valid state)
- Correct classical limit

This is not a choice --- it is the ONLY mathematically consistent way to "undo" a quantum channel given a prior. When someone asks "can we reverse the measurement?", the Petz map IS the answer.

### 2.3 When the Petz map fails: collapse

Define the recovery fidelity:

    F = F(rho, R_{sigma,N} o N(rho))

and the temporal asymmetry parameter:

    tau = 1 - F

The equivalence chain (Paper 1, Theorem):

    tau = 0  <==>  Sigma = 0  <==>  I(A;E|B) = 0  <==>  Quantum Markov Chain

where Sigma = D(rho || sigma) - D(N(rho) || N(sigma)) is the entropy production and I(A;E|B) is the conditional mutual information of the Stinespring dilation.

The universal recovery bound (Junge-Renner-Sutter-Wilde-Winter 2018):

    F >= exp(-Sigma/2)

Equivalently:

    tau <= 1 - exp(-Sigma/2)

**The operational meaning of collapse**: When tau --> 1, the Petz recovery map fails completely. From the system's perspective, the pre-measurement state is irrecoverable. This IS what we call "collapse" --- not a new process, but the FAILURE of the unique retrodiction protocol.

When tau --> 0, the channel is reversible, quantum coherence is maintained, no collapse occurs. The quantum eraser is the experimental demonstration of this regime.

### 2.4 The critical insight: perspective-dependence

This is where we go beyond standard decoherence:

- **From the universe's perspective** (the total Hilbert space H_S tensor H_E): evolution is unitary. Sigma = 0 for the global state. No collapse. All branches exist.

- **From the subsystem's perspective** (after tracing out the environment): Sigma > 0. The Petz recovery fails. The state LOOKS collapsed.

The APPEARANCE of collapse is not an illusion --- it is a **mathematical consequence** of the partial trace. The information about coherences has not been destroyed; it has been distributed across system-environment correlations that are inaccessible to the subsystem observer.

tau measures the DEGREE of this inaccessibility. It is not binary (collapsed/not collapsed) but continuous.

---

## 3. Is This Just "Decoherence Rebranded"? What Is New?

### 3.1 What standard decoherence provides

- Off-diagonal elements of the reduced density matrix decay: rho_{01}(t) --> 0
- Preferred basis selected by einselection (pointer states)
- Decoherence timescale computed: t_d ~ (lambda_dB / Delta x)^2 / gamma

### 3.2 What standard decoherence does NOT provide

**(a) Quantification of irrecoverability.** Standard decoherence tells you the coherences are "practically zero." It does not tell you how hard it would be to recover them. We provide an exact bound: F >= exp(-Sigma/2). This is a theorem, not a heuristic.

**(b) Connection to thermodynamics.** Standard decoherence computes a timescale. We connect collapse to entropy production --- the same quantity that governs the thermodynamic arrow of time. tau = 0 <==> Sigma = 0 <==> time-reversal symmetry. Collapse and the arrow of time are the SAME phenomenon, measured by the same parameter.

**(c) Unification across scales.** The same tau parameter describes:
- Microscopic: quantum eraser (tau ~ 0, recovery succeeds)
- Mesoscopic: Pikovski gravitational decoherence (tau ~ 0.004 for proposed BMV experiment)
- Macroscopic: Schrodinger's cat (tau --> 1, recovery fails)
- Gravitational: black hole information (tau = 0.632 at Planck radius)
- Neural: consciousness states (tau bounded by brain entropy production)

Standard decoherence does not provide this unification.

**(d) The selection problem: why one outcome?** See Section 4 below.

### 3.3 The structural addition: Petz map is the unique retrodiction

The key mathematical fact that distinguishes this from "decoherence plus interpretation" is:

> The Petz recovery map is the UNIQUE Bayesian retrodiction functor (Parzygnat-Buscemi 2023).

This means there is no alternative recovery protocol that could do better. When tau = 1, the state is irrecoverable not just in practice but in the information-theoretic sense --- no protocol, no matter how clever, can reconstruct the pre-measurement state from the post-measurement state. This is a theorem about the structure of quantum theory itself, not a practical limitation.

---

## 4. The Selection Problem: Why One Outcome?

This is the hardest question. Three layers of answer, in increasing depth.

### 4.1 Layer 1: Tensor product amplification (THEOREM)

For an N-body system under product channels, Theorem 1 of Paper 1b gives:

    Sigma_total >= N * sigma_min

    F_total <= exp(-N * sigma_min / 2)

For Schrodinger's cat (N ~ 10^26 internal DOF), even sigma_min ~ 10^{-26} per mode gives Sigma ~ 1 and tau ~ 0.4 within ~100 nanoseconds. By milliseconds, tau --> 1.

**This is not the central limit theorem.** The tensor product structure gives exponential-in-N suppression of coherence, not sqrt(N) fluctuations. This is structural --- it follows from the axiom that composite systems are described by tensor products.

**What this proves:** No physical protocol can recover the superposition once N * sigma_min >> 1. The superposition is not just "hard to see" --- it is information-theoretically destroyed from the subsystem's perspective.

**What this does NOT prove:** That only one branch "really exists." The full Hilbert space still contains all branches. But no measurement performable by any subsystem observer can distinguish the situation from genuine collapse.

### 4.2 Layer 2: Quantum Darwinism and Spectrum Broadcast Structure (ESTABLISHED PHYSICS)

Zurek's quantum Darwinism (2009) and the spectrum broadcast structure (Korbicz et al. 2024) provide the selection mechanism:

When the system interacts with N environmental modes (e.g., gravitational field modes via Pikovski coupling), the which-branch information is **redundantly encoded** in the environment:

    rho_SE = sum_i p_i |i><i| tensor rho^1_i tensor rho^2_i tensor ... tensor rho^N_i

where {|i>} are the pointer states and rho^k_i are distinguishable states for each environmental fragment k.

Key facts:
- Any fraction f > f_critical ~ 1/N of the environment contains full which-branch information
- The redundancy R ~ N means the outcome is recorded independently ~N times
- Riedel and Zurek (2017) explicitly confirmed this for gravitational decoherence

**The answer to "why one outcome?":** From the perspective of any subsystem observer (including other environmental modes), only one branch is accessible because every fragment of the environment has independently recorded the same outcome. Observing a different outcome would require simultaneously contradicting ~10^23 independent records --- which is precisely what tau = 1 quantifies as impossible.

### 4.3 Layer 3: Beyond FAPP --- causal irrecoverability (RESEARCH FRONTIER)

Three mechanisms can make tau = 1 **genuinely** irreversible (not just FAPP):

**(a) Danielson-Satishchandran-Wald (DSW) mechanism:** A superposition of masses radiates soft gravitons. These soft gravitons cross the cosmological horizon. Information is *causally* lost --- it is beyond the observer's causal past, forever. This is not "hard to recover" but "impossible in principle given causal structure."

**(b) Infrared superselection sectors:** Soft graviton emission changes the superselection sector. Recovery would require acting across sectors, which is algebraically forbidden within the observable algebra. This is genuine irrecoverability.

**(c) Galley-Giacomini-Selby theorem (Quantum 7, 1142, 2023):** Any consistent coupling between classical gravity and quantum matter is *fundamentally* irreversible --- no Stinespring dilation exists. If gravity has any classical aspect, tau > 0 from gravitational coupling is genuine, not FAPP.

**Current status:** These three mechanisms are at the research frontier. The DSW mechanism is the most established (published in PRD). If any of them holds, then collapse in our framework is genuinely irreversible --- not just practically irreversible but physically, algebraically, or causally irreversible.

---

## 5. Penrose's Objection and Our Response

### 5.1 What Penrose argues

Penrose's objective reduction (OR) program:

(P1) Quantum mechanics permits superpositions of geometries.
(P2) Distinct mass distributions produce distinct metrics.
(P3) A massive superposition IS a superposition of spacetimes.
(P4) Such superpositions are "ill-defined" because one cannot identify points between different manifolds.
(P5) Therefore gravity causes objective, non-unitary collapse at rate Lambda = E_G / hbar.

### 5.2 Where Penrose is right

We agree with P1--P3 completely. Gravity DOES play a special role in collapse. Our framework confirms this:

- Paper 2: Sigma_grav = 2 ln Q = -ln(-g_00). More mass --> more Sigma --> faster tau --> 1.
- Paper 1b, Theorem 2: E_G = (c^4 / 32 pi G) integral |grad(Delta Sigma)|^2 d^3x.

The Penrose collapse energy IS the gradient energy of the Sigma field. We reproduce the Penrose timescale:

    t_DP = hbar / E_G = 32 pi G hbar / (c^4 integral |grad(Delta Sigma)|^2 d^3x)

This is long for microscopic superpositions (grad Sigma small) and short for macroscopic ones (grad Sigma large) --- exactly as Penrose predicted.

### 5.3 Where Penrose is wrong

**P4 commits a category error.** Penrose argues that superposing two metrics requires "identifying points between manifolds," which is meaningless in GR. But:

Kabel, Brukner, and Wieland (2024, arXiv:2402.10267) proved that "identification is pointless" --- diffeomorphism-invariant quantum gravity renders point identification meaningless. Superpositions of geometries live in the Hilbert space H_geom built over Wheeler's superspace, where each point represents an ENTIRE 3-geometry. The superposition alpha |g^(1)> + beta |g^(2)> is well-defined as a vector in this Hilbert space. No point identification is needed.

**The error:** Penrose confuses pointwise addition of tensor fields on a single manifold (which requires point identification) with superposition in the gravitational configuration space (which does not). P4 does not follow from P3. Therefore P5 does not follow from P4.

### 5.4 What we provide instead

We need no new collapse law. We need:
- Standard quantum mechanics (unitary evolution on the total system)
- Gravity coupling to mass-energy (via Pikovski-type Hamiltonian H = H_int (1 + Phi(x)/c^2))
- The data-processing inequality (a mathematical theorem)
- The Petz recovery bound F >= exp(-Sigma/2) (a mathematical theorem)

From these alone: tau --> 1 for macroscopic objects, with a timescale that reproduces E_G.

### 5.5 The equivalence principle test

Equation (14) of Paper 1b reveals something Penrose's formula does not: **a uniform gravitational field cannot cause collapse.**

If Delta Sigma = const everywhere, then grad(Delta Sigma) = 0 and E_G = 0. Only **tidal** (non-uniform) gravitational effects produce collapse. This is the equivalence principle at work --- a freely falling lab cannot distinguish gravity from no gravity, so it cannot collapse superpositions.

Penrose's formula E_G ~ integral delta_rho delta_rho / |x-x'| also satisfies this, but in his framework it is a coincidence. In ours, it is a structural consequence of Sigma being a relative entropy (which measures *difference* from a reference, not absolute properties).

---

## 6. The Experimental Distinction

### 6.1 The recoherence test

The sharpest difference between Penrose OR and our framework:

**Penrose OR predicts:** Collapse is a genuine, irreversible, non-unitary process. Once Lambda = E_G/hbar triggers collapse, the superposition is destroyed forever. Recoherence probability = 0.

**We predict:** Collapse is operationally irreversible (exponentially hard) but in-principle reversible. If ALL N subsystems and their gravitational environment could be coherently controlled, the Petz map would succeed:

    F_recoh = exp(-Sigma/2) > 0

This is exponentially small but **nonzero**.

Current experiments cannot reach the required control for mesoscopic objects, but the logical structure differs fundamentally.

### 6.2 Already-excluded Penrose predictions

- **Gran Sasso underground experiments** (Donadi et al., Nat. Phys. 2021): Excluded the parameter-free Diosi-Penrose model via absence of anomalous radiation from spontaneous collapse.
- **170 kDa nanoparticle interference** (Pedalino et al., Nature 649, 866, 2026): Demonstrated coherence at masses where simplest DP predictions require visible decoherence.

Both results are consistent with our framework, which predicts no anomalous decoherence beyond the standard Sigma mechanism.

### 6.3 Concrete predictions

| Experiment | Our prediction | Penrose OR prediction |
|-----------|---------------|----------------------|
| BMV experiment (10^14 atoms) | tau ~ 0.004, measurable entanglement | May or may not show entanglement (depends on E_G vs. experimental timescale) |
| Micro-oscillator (10^14 atoms) | tau --> 1 in ~microseconds | Collapse in ~milliseconds (parameter-dependent) |
| 170 kDa molecule interference | Coherence maintained (Sigma << 1) | Borderline (some versions predict decoherence) |
| EEG Petz bound | F >= exp(-Sigma_brain/2) satisfied ~75.7% of epochs | No prediction |

---

## 7. The Honest Comparison Table

| Question | Penrose (OR) | Us (Petz) | Standard decoherence | Many-worlds |
|----------|-------------|-----------|---------------------|-------------|
| Why collapse? | Gravity (new law) | Sigma > 0 (open system + Petz failure) | Einselection + decoherence | No collapse (all branches real) |
| When? | t_DP = hbar/E_G | When tau --> 1 (= when Sigma >> 1) | When off-diagonal --> 0 | Never |
| New physics needed? | YES (non-unitary) | NO | NO | NO |
| Explains Born rule? | Claims to (via gravitational measure) | No (assumed) | No (assumed) | Contested (decision-theoretic) |
| Quantitative? | Yes (t_DP) | Yes (tau, Sigma, F) | Semi (decoherence time) | N/A |
| Why one outcome? | Gravity selects | Tensor product amplification + quantum Darwinism | FAPP (diagonal density matrix) | All outcomes occur |
| Genuinely irreversible? | Yes (by postulate) | Operationally yes; in-principle depends on DSW/superselection | No (in principle reversible) | No (all branches persist) |
| Role of gravity | Fundamental (drives collapse) | Special (Sigma_grav = 2 ln Q, enhances decoherence) | Optional (environmental) | None |
| Testable? | Gravity-induced collapse experiments | tau predictions, EEG bound, recoherence | Already confirmed | Unclear |
| Consistent with Gran Sasso / 170 kDa? | Under tension | Yes | Yes | Yes |

---

## 8. What We CAN'T Explain

### 8.1 The Born rule

We use p_i = |alpha_i|^2 as an axiom. We do not derive it.

Parzygnat and Buscemi's retrodiction framework *suggests* that the Born rule may be derivable from the Bayesian structure of the Petz map --- if the Petz map is the unique retrodiction functor, and if retrodiction must be consistent with a prior, then the prior determines the probabilities. But this derivation has not been completed. It is an open direction flagged in Paper 1b.

### 8.2 Do other branches "cease to exist"?

Our framework says: from the subsystem's perspective, other branches are irrecoverable (tau = 1). Whether they "really exist" in some ontological sense is a question we cannot answer --- it lies outside the operational definition.

Three positions are compatible with our mathematics:
- **Everettian**: All branches exist; tau measures the observer's inability to access them.
- **Relational**: Branches are observer-relative; tau quantifies the relativity.
- **Objective collapse with DSW**: If the DSW mechanism holds, information genuinely crosses the cosmological horizon, and branches are causally destroyed.

We prefer not to commit to one. The mathematics is the same either way.

### 8.3 Why does consciousness experience one outcome?

The hard problem of consciousness (Chalmers 1995) remains. Paper 1b's consciousness paper (paper_consciousness_sigma.tex) shows that brain entropy production Sigma_brain tracks consciousness level, and the Petz bound constrains retrodiction fidelity. But this does not explain WHY subjective experience exists or why it is associated with one branch.

We note honestly: a thermostat has nonzero tau. It presumably is not conscious. Tau is necessary but not sufficient for consciousness.

### 8.4 The Markov assumption

The brain dynamics is modeled as Markov. Real brains have long-range temporal correlations, 1/f scaling, memory effects. The bound holds for the Markov approximation but may be loose for non-Markovian processes. Extension to process tensors is an open problem.

---

## 9. The Deep Answer to "Why Collapse?"

After all the formalism, here is the answer in plain language.

### 9.1 Collapse happens because the universe has parts

The total universe evolves unitarily (Sigma = 0). But any subsystem is an open system --- it interacts with the rest of the universe. Tracing out the rest of the universe is a quantum channel. The entropy production of this channel is Sigma > 0. The Petz recovery bound then guarantees F < 1 and tau > 0.

**Collapse is the price of being a part, not the whole.**

### 9.2 Gravity makes it fast

Why does collapse happen on laboratory timescales rather than cosmological ones? Because gravity couples universally to mass-energy via position-dependent time dilation (Pikovski 2015). Every atom in a macroscopic object experiences a slightly different tick rate depending on its position in the gravitational field. This difference, multiplied across N ~ 10^{23} atoms, produces Sigma >> 1 almost instantly.

Without gravity: decoherence still happens (via electromagnetic, thermal, etc. interactions), but gravity is special because:
1. It is UNIVERSAL (couples to everything)
2. It cannot be shielded
3. The coupling is via the metric component g_00, which IS Sigma_grav

### 9.3 The Petz map is why collapse is irreversible

Other decoherence mechanisms might, in principle, be reversed by clever experimental control (e.g., spin echo). The Petz map tells you exactly how hard this is: F = exp(-Sigma/2). For a macroscopic object, Sigma ~ 10^{26}, so F ~ exp(-10^{26}/2) --- a number so small it has no physical meaning. The recovery is "possible in principle" only in the mathematical sense that exp(-10^{26}/2) > 0.

### 9.4 The selection is not a process --- it is a perspective

From the total Hilbert space: no selection occurred. All branches coexist.

From the subsystem: one branch is experienced because tau = 1 means the observer CANNOT access information about other branches. The quantum Darwinism mechanism ensures that every environmental fragment independently confirms the same branch. Experiencing a "different" outcome would require simultaneously overwriting ~10^{23} independent records --- which tau = 1 certifies as impossible.

**The answer to "why one outcome?":** Because you are a subsystem, and tau = 1 means no subsystem can access the other branches.

---

## 10. Summary: What We Have and What We Don't

### THEOREMS (proven within standard QM + information theory):

1. **tau = 0 <==> Sigma = 0 <==> Quantum Markov Chain** (equivalence chain)
2. **F >= exp(-Sigma/2)** (universal recovery bound, JRSWW 2018)
3. **Sigma_total >= N * sigma_min** for N-body systems (extensive entropy production)
4. **E_G = (c^4/32piG) integral |grad(Delta Sigma)|^2 d^3x** (Penrose energy as Sigma gradient)
5. **K(Q) = mu^2 (Q-1)^2** from coherent-state QRE (ghost condensation from information geometry)
6. **The Petz map is the unique retrodiction functor** (Parzygnat-Buscemi 2023)

### ESTABLISHED PHYSICS (standard results we apply):

7. Gravitational time dilation causes decoherence (Pikovski et al. 2015)
8. Quantum Darwinism selects pointer states (Zurek 2009)
9. Spectrum broadcast structure for von Neumann interactions (Korbicz et al. 2024)
10. Wheeler-DeWitt equation is timeless; Page-Wootters mechanism recovers time (1983/2015)

### INTERPRETIVE PROPOSALS (well-motivated but not proven):

11. Collapse IS Petz recovery failure (identification, not derivation)
12. The five-level ontological structure (Level 0--4)
13. Consciousness requires intermediate Sigma (from paper_consciousness_sigma.tex)
14. Branches are genuinely destroyed via DSW mechanism / superselection (research frontier)

### OPEN PROBLEMS:

15. Born rule derivation from Petz structure
16. Ontological status of branches
17. Non-Markovian extension of the bound
18. Complete Level 0 derivation from UV-complete QG
19. Hard problem of consciousness

---

## 11. Connection to the Core Thesis

Recall Sheng-Kai Huang's core insight from quantum eraser experiments:

> In an open environment, the paired quantum's future observation affects the present --- but you must WAIT for the future to actually happen to verify this retroactively.

> In a **closed** (zero-entropy) environment, you don't need to wait, because the result is already determined.

Collapse is the open-environment case. The system interacts with the environment (Sigma > 0), and from the subsystem's perspective, the future state is determined --- you cannot undo it, you cannot access the other branches, you cannot retroduce the past. This is tau = 1.

The quantum eraser is the closed case. The system + idler form a closed system (Sigma = 0 for the total), and the "collapse" CAN be undone by measuring the idler appropriately. This is tau = 0 (or close to it).

**Collapse is what happens when you are in an open system (Sigma > 0) and cannot perform the quantum eraser because the environment has too many degrees of freedom to control (N >> 1).**

The Petz recovery map is the generalization of the quantum eraser to arbitrary channels. When it succeeds (tau ~ 0), you get "uncollapse" (quantum eraser). When it fails (tau ~ 1), you get "collapse" (measurement). There is no sharp boundary --- tau is continuous. The discreteness of outcomes emerges from the quantum Darwinism mechanism (redundant encoding in the environment).

---

## References

### Mathematical foundations
- Petz (1986, 1988): Recovery map
- Parzygnat & Buscemi (2023, Quantum 7, 1013): Unique retrodiction functor
- Junge, Renner, Sutter, Wilde, Winter (2018, Ann. Henri Poincare 19, 2955): Universal recovery bound
- Fawzi & Renner (2015, CMP 340, 575): Approximate quantum Markov chains

### Decoherence and measurement
- Zurek (2003, RMP 75, 715): Einselection and decoherence
- Schlosshauer (2004, RMP 76, 1267): Decoherence and measurement problem review
- Pikovski, Zych, Costa, Brukner (2015, Nat. Phys. 11, 668): Gravitational time dilation decoherence
- Korbicz et al. (2024): Spectrum broadcast structure

### Penrose and collapse
- Penrose (1996, GRG 28, 581): Gravity's role in quantum state reduction
- Diosi (1987, Phys. Lett. A 120, 377): Gravitational collapse model
- Kabel, Brukner, Wieland (2024, arXiv:2402.10267): "Identification is pointless"
- Donadi et al. (2021, Nat. Phys. 17, 74): Gran Sasso exclusion of DP model
- Pedalino et al. (2026, Nature 649, 866): 170 kDa interference

### Beyond FAPP
- Danielson, Satishchandran, Wald (2022, PRD 105, 086001): Gravitational decoherence via soft gravitons
- Galley, Giacomini, Selby (2023, Quantum 7, 1142): No Stinespring for classical-quantum coupling
- Oppenheim (2023, PRX 13, 041040): Post-quantum classical gravity

### This framework
- Huang, Paper 1: tau = 1 - F, equivalence chain, Petz recovery bound
- Huang, Paper 1b: Collapse = Petz recovery failure (Theorems 1--3)
- Huang, Paper 2: Sigma_grav = 2 ln Q = -ln(-g_00)
- Huang, paper_consciousness_sigma.tex: Brain dynamics as Petz recovery

---

*Last updated: 2026-03-20*
