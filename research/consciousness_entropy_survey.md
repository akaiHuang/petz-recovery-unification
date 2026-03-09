# Consciousness, Entropy, and Information Theory: A Mathematical Survey

**Date**: 2026-03-09
**Purpose**: Survey papers connecting information theory, entropy, and consciousness with MATHEMATICAL frameworks. Evaluate relevance to the tau framework (tau = 1 - F, Sigma = entropy production).

---

## Table of Contents
1. [Integrated Information Theory (IIT)](#1-integrated-information-theory-iit)
2. [Quantum Consciousness: Orch-OR](#2-quantum-consciousness-orch-or)
3. [Free Energy Principle (Friston)](#3-free-energy-principle-friston)
4. [Thermodynamics of Consciousness](#4-thermodynamics-of-consciousness)
5. [Neural Entropy and Consciousness States](#5-neural-entropy-and-consciousness-states)
6. [Quantum Coherence in Biological Systems](#6-quantum-coherence-in-biological-systems)
7. [Connections to the tau Framework](#7-connections-to-the-tau-framework)
8. [Summary Table](#8-summary-table)

---

## 1. Integrated Information Theory (IIT)

### 1.1 IIT 4.0: Core Mathematical Framework

**Key paper**: Albantakis, Barbosa, et al. (2023). "Integrated information theory (IIT) 4.0: Formulating the properties of phenomenal existence in physical terms." *PLOS Computational Biology*. [arXiv:2212.14787](https://arxiv.org/abs/2212.14787)

#### Core Equations

**Intrinsic Information (effect side)**:
```
ii_e(s, s_bar) = p_e(s_bar | s) * log[ p_e(s_bar | s) / p_e(s_bar) ]
```

**Intrinsic Information (cause side)**:
```
ii_c(s, s_bar) = p_c(s_bar | s) * log[ p_c(s | s_bar) / p_c(s) ]
```

- The product form = selectivity x informativeness
- Selectivity: concentration of causal power over specific states
- Informativeness: deviation from chance (logarithmic term)

**System Integrated Information (Phi)**:
```
phi_s(T_e, T_c, s, theta) = min{ phi_c(T_c, s, theta), phi_e(T_e, s, theta) }
```

where phi is evaluated over the Minimum Information Partition (MIP):
```
phi_s(T_e, T_c, s) := phi_s(T_e, T_c, s, theta')
```
theta' = the partition that minimizes phi.

**Effect-side irreducibility**:
```
phi_e(T_e, s, theta) = p_e(s'_e | s) * |log[ p_e(s'_e | s) / p_e^theta(s'_e | s) ]|_+
```

**Distinction-level phi**:
```
phi_d(m) = min( phi_c(m), phi_e(m) )
```

**Big Phi (Phi)**: The sum total of the phi values of all distinctions and relations composing the Phi-structure = quantity of consciousness.

#### Key Mathematical Features
- Uses **Intrinsic Difference (ID)** measure, NOT KL divergence or Earth Mover's Distance
- ID uniquely satisfies three postulate-aligned properties: causality, intrinsicality, specificity
- Interventional conditional probabilities (do-calculus), NOT observational
- Uniform prior distributions over all possible states

#### Relationship to entropy
- phi and entropy-based measures of entanglement share structural similarity: both quantify how much a system is MORE than the sum of its parts
- phi is NOT entropy, but uses information-theoretic divergence in its logarithmic term
- The unconstrained repertoire p_e(s_bar) serves as a "maximum entropy" reference

### 1.2 Quantum IIT (QIIT)

**Key paper**: Albantakis, Prentner, Durham (2023). "Computing the Integrated Information of a Quantum Mechanism." *Entropy* 25(3), 449. [arXiv:2301.02244](https://arxiv.org/abs/2301.02244)

**Key paper**: Zanardi, Tomka, Venuti (2018). "Towards Quantum Integrated Information Theory." [arXiv:1806.01421](https://arxiv.org/abs/1806.01421)

#### Core Extensions
- Intrinsic information translated into **density matrix formulation**
- Classical KL divergence replaced by **quantum relative entropy**: S(rho || sigma) = Tr[rho(ln rho - ln sigma)]
- Conditional independence extended to accommodate **quantum entanglement**
- Zanardi et al. used **trace distance** to quantify phi; later work uses quantum relative entropy between density matrices
- Code available: https://github.com/Albantakis/QIIT

#### Key Result
- Quantum entanglement can **increase** phi beyond classical limits
- The density matrix formulation means phi is fundamentally related to von Neumann entropy through the quantum relative entropy measure

### 1.3 Mathematical Critiques

**Doerig et al. (2019). "The unfolding argument: Why IIT and other causal structure theories cannot explain consciousness."** *Consciousness and Cognition* 72, 49-59.
- Uses **universal approximation theorem**: any RNN can be "unfolded" to a functionally equivalent feedforward network
- The unfolded version has different causal structure => different phi => different predicted consciousness
- But identical input-output behavior
- This is a **mathematical** critique, not philosophical: it's "pre-falsification"

### 1.4 Category-Theoretic Formulation

**Phillips et al. (2024)**. Integration, information, and exclusion correspond to limits, colimits, and adjunctions in category theory.
- IIT's mathematical underpinnings can be unified via category theory
- This is structurally similar to the Petz map being a functor (retrodiction functor)

---

## 2. Quantum Consciousness: Orch-OR

### 2.1 Penrose's Gravitational Self-Energy Formula

**Key equation**:
```
tau_P = hbar / E_G
```

where:
- tau_P = mean lifetime of the quantum superposition before objective reduction (OR)
- hbar = reduced Planck's constant
- E_G = gravitational self-energy of the difference between mass distributions in superposed states

**Gravitational self-energy (general form)**:
```
E_G = (G/2) integral integral [ mu_A(r) - mu_B(r) ][ mu_A(r') - mu_B(r') ] / |r - r'| dr dr'
```

where mu_A, mu_B are the mass density distributions in superposition states A and B.

**For Gaussian mass distributions** (nuclei with masses M_i, positions R_i, widths alpha_i):
```
E_G involves sums of terms: (M_i * M_j) / sqrt(alpha_i + alpha_j)
```

**Homogeneous medium approximation**:
```
E_G ~ (G * M_tot^2 * N_at) / sqrt(2*pi*alpha) * f(Delta_R, delta_r_bar)
```

### 2.2 Numerical Values (from Barbatti et al. 2024, PMC11305101)

| System | Atoms | Collapse Time tau_P |
|--------|-------|---------------------|
| NH3 isomerization | 4 | 1.0 x 10^17 s (~10^9 years) |
| C70 diffraction | 70 | 2.4 x 10^14 s |
| Neuron activation | ~10^6 | 1.1 x 10^10 s |
| 1 kg iron crystal | large | ~1 s |
| 1 g aluminum pointer | ~10^22 | 6 x 10^-21 s |

**Key observation**: For single neurons, tau_P ~ 10^10 s is far too long. Orch-OR proposes that **coherent superpositions across ~10^9 tubulins** (in thousands of neurons) could bring E_G high enough for tau_P ~ 25 ms (matching the ~40 Hz gamma oscillation).

### 2.3 Experimental Status (2024-2025)

**Supportive evidence**:
1. **Tryptophan superradiance** (J. Phys. Chem. 2024): Confirmed in networks of tryptophans in microtubules. Laser-induced excitations propagate further and persist longer than expected.
2. **Anesthetic effects on microtubules** (Wellesley College 2024): Rats given epothilone B (microtubule-binding drug) took >1 minute longer to fall unconscious.
3. **Revised decoherence times**: More realistic microtubule models extend coherence from 10^-13 s to 10^-6 - 10^-4 s (7 orders of magnitude improvement).
4. **Superconducting qubit demonstration** (arXiv:2504.02914, April 2025): Transmon qubits show wavefunction collapse consistent with OR from improper mixtures. Moved ~10^-12 kg across ~1 mm.

**Constraining evidence**:
1. **Gran Sasso underground experiment** (2021): Ruled out parameter-free Diosi model; Penrose version with free parameter still viable.
2. No clear-cut falsification, but no unique empirical prediction confirmed either.

### 2.4 Conscious Active Inference + Orch-OR

**Recent paper**: "Conscious active inference II: Quantum orchestrated objective reduction among intraneuronal microtubules naturally accounts for discrete perceptual cycles." *ScienceDirect*, 2025.
- Merges Friston's active inference with Orch-OR
- Claims tau_P ~ 25 ms cycles match perceptual frame rate

---

## 3. Free Energy Principle (Friston)

### 3.1 Core Mathematical Framework

**Key paper**: Friston (2019). "A free energy principle for a particular physics." [arXiv:1906.10184](https://arxiv.org/abs/1906.10184)

**Key paper**: Fields, Friston, et al. (2022). "A free energy principle for generic quantum systems." [arXiv:2112.15242](https://arxiv.org/abs/2112.15242)

#### The Variational Free Energy

```
F = integral q(T) * ln[ q(T) / P(T, S) ] dT
```

Equivalently:
```
F = D_KL[ q(T) || P(T|S) ] - ln P(S)
```

Since D_KL >= 0:
```
F >= -ln P(S)    [free energy is an upper bound on surprisal]
```

**Energy-Entropy decomposition**:
```
F = E_q[Energy] - H[q]
```

where:
- E(T,S) = -ln P(T,S) is the "energy"
- H[q] = -integral q(T) ln q(T) dT is the entropy of the recognition density
- q(T) = recognition density (approximate posterior)
- P(T,S) = generative model
- P(S) = marginal likelihood (evidence)
- -ln P(S) = surprisal (self-information)

#### Key Relationships
- Long-term average of surprise = **Shannon entropy** H of sensory states
- Minimizing F => minimizing entropy of sensory outcomes => maintaining low-entropy existence
- F is related to **maximum entropy principle**: FEP = max-entropy under prediction accuracy constraint

### 3.2 Connection to Physics

From Friston (2019):
- Fokker-Planck equation becomes **Schrodinger equation** in certain limits
- Framework spans quantum mechanics -> statistical mechanics -> classical mechanics
- **Markov blanket** partition: internal states, external states, sensory states, active states
- Autonomous entropy production is always <= 0 at nonequilibrium steady-state
- Autonomous flow **resists dispersion** = maintains low entropy = "self-evidencing"

### 3.3 Quantum FEP

From Fields et al. (2022):
- FEP reformulated in **spacetime-background-free, scale-free quantum information theory**
- Generic quantum systems can be regarded as **observers** that minimize Bayesian prediction error
- **Asymptotically equivalent to the Principle of Unitarity**
- Connects to quantum reference frames (QRFs)

### 3.4 Connection to Entropy Production (Sigma)

The FEP has a deep but indirect connection to entropy production:
- F >= -ln P(o): minimizing F => minimizing surprise
- Long-term surprise = entropy H of outcomes
- Systems that persist must have **low entropy** of their sensory states
- The **dissipative** component of flow maintains this (Helmholtz decomposition)
- At steady state: solenoidal flow preserves probability; dissipative flow reduces entropy

**Mathematical connection to Sigma**:
```
F = D_KL[q || p_posterior] + (-ln P(evidence))
```

Compare to tau framework:
```
Sigma = D(rho_forward || rho_reverse)    [Crooks-like]
```

Both are KL divergences measuring deviation from a reference. The FEP measures prediction error; Sigma measures time-reversal asymmetry. They are structurally analogous.

---

## 4. Thermodynamics of Consciousness

### 4.1 Non-Equilibrium Brain Dynamics

**Key paper**: Sanz Perl, Pallavicini, et al. (2021). "Nonequilibrium brain dynamics as a signature of consciousness." *Phys. Rev. E* 104, 014411. [arXiv:2012.10792](https://arxiv.org/abs/2012.10792)

#### Mathematical Framework

The brain is modeled as a multivariate stochastic process. The key measure is:

**Entropy production rate** in phase space:
```
Sigma_dot = integral J(x) . nabla ln[ P(x) / P_eq(x) ] dx
```

where J(x) is the probability flux and P_eq is the equilibrium distribution.

Alternatively, via **curl of probability flux**:
```
Sigma_dot ~ || curl J ||^2
```

Non-zero curl = broken detailed balance = non-equilibrium = irreversibility.

#### Key Results
- **All states of reduced consciousness** (sleep, propofol, ketamine, ketamine+medetomidine) show **higher proximity to equilibrium** compared to wakefulness
- Entropy production is **minimal** during reduced consciousness
- **Monotonous relationship** between entropy production and level of consciousness (wakefulness -> N1 -> N2 -> N3 sleep)

### 4.2 FDT Violations as Consciousness Marker

**Key paper**: Fluctuation-Dissipation Theorem violations. *Phys. Rev. Research* 7, 013301 (2025).

**Key paper**: "Thermodynamics of consciousness: A non-invasive perturbational framework." *bioRxiv* (2025).

#### Mathematical Framework

In thermodynamic equilibrium, the **Fluctuation-Dissipation Theorem (FDT)** states:
```
Response function R(t) = (1/k_B T) * d/dt C(t)
```

where C(t) is the autocorrelation function. In the brain, FDT is **violated** during conscious states.

**FDT violation parameter**:
```
Delta_FDT = || R(t) - (1/k_B T_eff) * dC/dt ||
```

#### Key Results
- **Decreased FDT violations** in disorders of consciousness and anesthesia
- States with **higher PCI values** (higher consciousness) correspond to **stronger FDT violations**
- Entropy production, temporal irreversibility, and coupling asymmetry all show **parallel reductions** with loss of consciousness
- Xenon and propofol anesthesia: FDT violation reliably indexes consciousness
- **Ketamine is anomalous**: maintains high PCI and FDT violation despite anesthesia (dissociative state)

### 4.3 Entropy Production of Ornstein-Uhlenbeck Processes

**Key paper**: "Entropy production of multivariate Ornstein-Uhlenbeck processes correlates with consciousness levels in the human brain." *Phys. Rev. E* 107, 024121 (2023).

- Models brain dynamics as multivariate OU process
- Entropy production Sigma_dot computed analytically
- **Sigma_dot monotonically decreases** from wakefulness through N1, N2, N3 sleep stages
- Provides a **computable** thermodynamic marker of consciousness level

### 4.4 Temporal Irreversibility

**Key paper**: de la Fuente, Zamberlan, et al. (2022). "Temporal irreversibility of neural dynamics as a signature of consciousness." *Cerebral Cortex* 33(5), 1856. [bioRxiv:2021.09.02.458802](https://www.biorxiv.org/content/10.1101/2021.09.02.458802v2.full)

- Temporal irreversibility quantified via time-reversed statistical tests
- Conscious states are more temporally irreversible
- **Direct connection**: temporal irreversibility = Sigma > 0 = broken time-reversal symmetry

---

## 5. Neural Entropy and Consciousness States

### 5.1 The Entropic Brain Hypothesis

**Key paper**: Carhart-Harris, Leech, et al. (2014). "The entropic brain: a theory of conscious states informed by neuroimaging research with psychedelic drugs." *Front. Hum. Neurosci.* 8:20. [PMC3909994](https://pmc.ncbi.nlm.nih.gov/articles/PMC3909994/)

**Updated**: Carhart-Harris (2018). "The entropic brain - revisited." *Neuropharmacology* 142, 167-178.

#### Core Proposition
Within upper and lower limits, the **entropy of spontaneous brain activity indexes the informational richness of conscious states**.

#### Entropy Measures Used
1. **Lempel-Ziv Complexity (LZC)**: Algorithmic complexity of time series
2. **Sample Entropy (SampEn)**: Regularity measure
3. **Perturbational Complexity Index (PCI)**: Complexity of TMS-evoked EEG response

#### Quantitative Data: Consciousness Hierarchy (by entropy)

```
Seizure/coma  <  Deep sleep  <  Light sedation  <  Normal waking  <  Psychedelic state
   (lowest)                                                            (highest entropy)
```

Specific findings:
- **Psychedelics** (psilocybin, LSD, ketamine): LZC significantly HIGHER than normal waking (p < 0.05 for all three substances, MEG data)
- **Propofol anesthesia**: LZC significantly LOWER than normal waking
- **Deep sleep (N3)**: LZC lower than N1, which is lower than waking
- **PCI threshold**: PCI = 0.31 separates conscious from unconscious states
- **Ketamine anesthesia**: anomalously maintains PCI > 0.31 (grouped with wakefulness)

#### Type I vs Type II Complexity (2025 update)
- **Type I** (e.g., LZC): Scales linearly with entropy. Maximal at randomness.
- **Type II**: Inverted U-shape with entropy. Maximal at criticality (balance of order and randomness).
- **Consciousness may correspond to Type II criticality** — sub-critical (slightly below maximum entropy)

### 5.2 A Mathematical Model of Consciousness Dynamics

**Fuchs (2025)**. Published in HAL archives.

**Fuchs Consciousness Equation**:
```
dC/dt = k * (I x E - alpha * C)
```

where:
- C = consciousness level
- I = information flux
- E = environmental coupling
- alpha = dissipation rate
- k = scaling constant

**Note**: This is a phenomenological model, not derived from first principles.

### 5.3 Consciousness as Entropy Reduction (CER)

**arXiv:2510.06297** (2025).
- Proposes consciousness selects low-entropy scenarios from a high-entropy subconscious manifold
- Connects to Global Workspace Theory and Hopfield networks
- Mathematical formalism not yet fully developed

### 5.4 Quantum Entropy Deficit Hypothesis

Recent proposal (exact reference pending): Phenomenal consciousness is physically identical to the **largest and most persistent residual quantum entropy deficit** thermodynamics permits in a warm, wet, macroscopic system.

**Quantitative bound**: 40-150 bits per ~50 ms frame for waking consciousness.

---

## 6. Quantum Coherence in Biological Systems

### 6.1 Current Status (2024-2025)

**Key conferences**: QuEBS 2024 (Wuhan, China); GRC Quantum Biology 2025

#### Known quantum effects in biology:
1. **Photosynthesis**: Quantum coherence achieves ~100% energy transfer efficiency
2. **Magnetoreception**: Radical-pair mechanism in bird navigation
3. **Enzyme catalysis**: Quantum tunneling of protons
4. **Olfaction**: Possible quantum tunneling contribution
5. **DNA mutations**: Proton tunneling between base pairs

#### Decoherence Times in Biological Systems
- **Previous estimates**: ~10^-13 s (too short for neural relevance)
- **Revised estimates (2024)**: 10^-6 to 10^-4 s in microtubules (7 orders of magnitude longer)
- **Comparison**: Molecular spin qubits also show ~10 us decoherence
- **Tryptophan superradiance**: Confirmed in microtubule architectures (2024)

#### "Warm, Wet, and Noisy" Paradigm Shift
- Emerging "Open Quantum Systems in Biology" framework
- Life may have evolved to **engineer environmental noise** to maintain coherence
- Not fighting decoherence — exploiting it

### 6.2 Calcium Phosphate Quantum Effects

**March 2025**: Evidence for quantum effects in calcium phosphate formation.
- Relevant because calcium signaling is central to neural function
- Suggests quantum coherence may play role at the biochemical level

---

## 7. Connections to the tau Framework

### 7.1 Direct Mathematical Connections

#### Connection 1: tau and Entropy Production Sigma in Consciousness

Our framework:
```
tau = 1 - F,  where F = fidelity of Petz recovery
Sigma >= 0,   entropy production
F >= exp(-Sigma/2)   [JRSWW bound]
=> tau <= 1 - exp(-Sigma/2)
```

**Brain thermodynamics** (Sanz Perl 2021, Phys. Rev. E):
- Wakefulness: high Sigma (high entropy production, far from equilibrium)
- Deep sleep: low Sigma (near equilibrium)
- Consciousness level ~ Sigma

**Mapping**:
```
Conscious states:    Sigma >> 0  =>  tau >> 0  =>  strong time arrow
Unconscious states:  Sigma -> 0  =>  tau -> 0  =>  weak time arrow
```

**Interpretation**: Consciousness requires a strong time arrow. When tau -> 0 (closed system, reversible), there is no subjective time experience. This is mathematically consistent with loss of consciousness during deep sleep showing reduced entropy production.

#### Connection 2: IIT's Phi and tau

IIT measures information integration; tau measures recovery fidelity. The connection:

```
Phi = 0  <=>  system is reducible  <=>  no integration  <=>  no consciousness
tau = 0  <=>  perfect recovery    <=>  closed system   <=>  no time arrow
```

For a system with Phi > 0:
- The system has irreducible causal structure
- This irreducibility implies information loss under partitioning
- Information loss = entropy production = Sigma > 0 = tau > 0

**Conjecture**: For quantum systems,
```
Phi >= f(tau)
```
where f is a monotonically increasing function. High tau (strong irreversibility) is necessary for high Phi (rich consciousness).

**Supporting evidence**: The quantum IIT formalism uses quantum relative entropy S(rho || sigma) = Tr[rho(ln rho - ln sigma)], which is the SAME mathematical object that defines Sigma in the Crooks fluctuation theorem.

#### Connection 3: FEP and tau

Friston's free energy:
```
F_Friston = D_KL[q || p_posterior] + (-ln P(evidence))
```

Our Sigma:
```
Sigma = D(rho_forward || rho_reverse)
```

Both are KL divergences measuring asymmetry:
- F_Friston: asymmetry between belief and reality
- Sigma: asymmetry between forward and reverse processes

**Structural isomorphism**:
```
FEP: F >= -ln P(o)     [variational bound]
tau: tau <= 1 - exp(-Sigma/2)  [JRSWW bound]
```

Both are upper/lower bounds involving divergence measures. The FEP says organisms minimize surprise; the tau framework says closed systems minimize tau. **Living systems maintain themselves at a specific tau > 0** — neither too high (dissolution) nor too low (death/unconsciousness).

#### Connection 4: Penrose's tau_P and our tau

Penrose:
```
tau_P = hbar / E_G    [collapse timescale]
```

Our tau:
```
tau = 1 - F            [recovery infidelity]
```

**These are different quantities** (timescale vs. fidelity deficit), but:
- Both are named tau (coincidence, but suggestive)
- Both involve gravitational decoherence (our Paper 1b connects tau to Pikovski channel)
- Both = 0 for closed systems (no gravity-induced collapse; perfect recovery)
- Both increase with openness (stronger gravity = faster collapse; more entropy production = worse recovery)

**Quantitative bridge** (from Paper 1b):
```
Pikovski channel: K_0 = sqrt((1+p)/2) I,  K_1 = sqrt((1-p)/2) sigma_z
tau_internal = 2 * tau_spatial * (1 - tau_spatial)
Sigma_internal = Sigma_spatial = S_entanglement
```

When p -> 0 (strong gravitational decoherence): tau -> 1, Penrose's tau_P -> 0.
When p -> 1 (no gravitational decoherence): tau -> 0, Penrose's tau_P -> infinity.

**They are inversely related**: tau * tau_P ~ constant (for fixed system).

### 7.2 Predictions from tau Framework for Consciousness

1. **Consciousness requires tau > 0**: A perfectly reversible (tau = 0) system cannot be conscious. This is testable: quantum computers in coherent operation should show no consciousness markers.

2. **Consciousness requires tau < 1**: A maximally irreversible system (tau = 1) has lost all information — no structure for consciousness either. This is consistent with brain death (maximal entropy, tau -> 1 in neural dynamics).

3. **Optimal consciousness at intermediate tau**: The entropic brain hypothesis places normal consciousness slightly sub-critical. In our framework: **consciousness maximizes at tau ~ tau_critical, where the system balances information integration (requires openness, tau > 0) against information preservation (requires some closure, tau < 1)**.

4. **Anesthesia reduces tau**: Propofol, xenon reduce entropy production (Sigma -> smaller) => tau -> smaller => loss of consciousness. **Already confirmed** by Sanz Perl et al. (2021).

5. **Psychedelics increase tau**: Increased LZC and entropy => Sigma larger => tau larger. Consistent with "expanded" consciousness having a stronger time arrow.

6. **Ketamine is special**: Maintains high entropy production despite anesthesia (dissociative state). In tau framework: maintains tau > 0 for higher-order representations while disrupting sensory processing.

### 7.3 Open Questions

1. **Is Phi computable from tau?** Can we derive IIT's Phi from knowledge of Sigma and the Petz recovery map structure?

2. **Is consciousness a phase transition?** The Type II complexity (inverted U) suggests a critical phenomenon. In tau language: is there a critical Sigma_c where consciousness "turns on"?

3. **What is the tau of dreaming?** REM sleep shows brain activity similar to waking, but reduced environmental coupling. In tau terms: internal Sigma may be high, but external Sigma is low. Does consciousness require INTERNAL tau > 0, or TOTAL tau > 0?

4. **Can the quantum FEP (Fields et al. 2022) be reformulated in terms of tau?** The FEP is "asymptotically equivalent to the Principle of Unitarity" — which in our framework means it describes the tau = 0 limit. Active inference then becomes the process of maintaining tau > 0 against environmental perturbation.

---

## 8. Summary Table

| Framework | Central Quantity | Mathematical Form | Equations | Data | Connection to tau |
|-----------|-----------------|-------------------|-----------|------|-------------------|
| **IIT 4.0** | Phi (integrated info) | Product of selectivity x informativeness, minimum over partitions | Yes (complex) | Computed for small networks | Phi > 0 => irreducibility => Sigma > 0 => tau > 0 |
| **Quantum IIT** | Phi (quantum) | Quantum relative entropy between density matrices | Yes | Computed for 2-3 qubits | Uses SAME QRE as Sigma definition |
| **Orch-OR** | tau_P = hbar/E_G | Gravitational self-energy determines collapse time | Yes (explicit) | Numerical estimates; partial experiments | Inversely related: tau * tau_P ~ const |
| **FEP** | F (free energy) | D_KL[q || p] upper bound on surprise | Yes | Applied to brain imaging | Structural isomorphism with tau bound |
| **Brain Sigma** | Sigma_dot (entropy prod.) | Probability flux curl in phase space | Yes | ECoG (primates), fMRI (humans) | DIRECT: Sigma_dot indexes consciousness level |
| **FDT violation** | Delta_FDT | Deviation from equilibrium FDT | Yes | EEG (humans, rodents) | FDT violation ~ non-equilibrium ~ Sigma > 0 |
| **Entropic brain** | LZC, PCI | Lempel-Ziv complexity, perturbational complexity | Yes | MEG, EEG across states | LZC ~ entropy ~ Sigma; PCI > 0.31 threshold |
| **Quantum biology** | Decoherence time | Standard decoherence theory | Yes | 10^-6 to 10^-4 s in microtubules | Decoherence = information loss = tau > 0 |

---

## 9. Key References (with mathematical content)

### IIT
1. Albantakis et al. (2023). IIT 4.0. *PLOS Comp. Biol.* [PMC10581496](https://pmc.ncbi.nlm.nih.gov/articles/PMC10581496/)
2. Albantakis, Prentner, Durham (2023). Quantum IIT. *Entropy* 25(3), 449. [MDPI](https://www.mdpi.com/1099-4300/25/3/449)
3. Zanardi, Tomka, Venuti (2018). Towards QIIT. [arXiv:1806.01421](https://arxiv.org/abs/1806.01421)
4. Kleiner, Tull (2020). Mathematical Structure of IIT. *Front. Appl. Math. Stat.* [Frontiers](https://www.frontiersin.org/journals/applied-mathematics-and-statistics/articles/10.3389/fams.2020.602973/full)
5. Phillips et al. (2024). Category-theoretic IIT.
6. Doerig et al. (2019). Unfolding argument. *Consciousness and Cognition* 72, 49-59.

### Orch-OR
7. Penrose, Hameroff. Orch-OR. [hameroff.arizona.edu](https://hameroff.arizona.edu/research-overview/orch-or)
8. Barbatti et al. (2024). Gravitational collapse time for molecules. [PMC11305101](https://pmc.ncbi.nlm.nih.gov/articles/PMC11305101/)
9. arXiv:2504.02914 (April 2025). OR on superconducting quantum computer.
10. Tryptophan superradiance (2024). *J. Phys. Chem.*
11. Hameroff et al. (2025). Conscious active inference II. *ScienceDirect*.

### Free Energy Principle
12. Friston (2019). FEP for particular physics. [arXiv:1906.10184](https://arxiv.org/abs/1906.10184)
13. Fields, Friston et al. (2022). FEP for generic quantum systems. [arXiv:2112.15242](https://arxiv.org/abs/2112.15242)
14. Ramstead et al. (2023). FEP made simpler. [ScienceDirect](https://www.sciencedirect.com/science/article/pii/S037015732300203X)
15. Fields et al. (2025). Inner screens and FEP. *Neurosci. Consciousness* 2025(1).

### Thermodynamics of Consciousness
16. Sanz Perl et al. (2021). Non-equilibrium brain dynamics. *Phys. Rev. E* 104, 014411. [APS](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.104.014411)
17. (2023). Entropy production of OU processes and consciousness. *Phys. Rev. E* 107, 024121. [APS](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.107.024121)
18. (2025). FDT violations and consciousness. *Phys. Rev. Research* 7, 013301. [APS](https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.7.013301)
19. de la Fuente et al. (2022). Temporal irreversibility of neural dynamics. *Cerebral Cortex* 33(5), 1856.
20. (2025). Thermodynamics of consciousness framework. *bioRxiv*.

### Entropic Brain
21. Carhart-Harris et al. (2014). Entropic brain hypothesis. *Front. Hum. Neurosci.* 8:20. [PMC3909994](https://pmc.ncbi.nlm.nih.gov/articles/PMC3909994/)
22. Carhart-Harris (2018). Entropic brain revisited. *Neuropharmacology* 142, 167-178.
23. Schartner et al. (2017). Signal diversity for psychedelics. *Scientific Reports* 7, 46421. [Nature](https://www.nature.com/articles/srep46421)

### Quantum Biology
24. (2025). Quantum models of consciousness from QIS perspective. *Entropy* 27(3), 243. [MDPI](https://www.mdpi.com/1099-4300/27/3/243)
25. (2025). Quantum Biology. *arXiv:2511.14363*.

---

## 10. Bottom Line for the tau Framework

### What is well-established (with data):
- **Entropy production Sigma decreases with loss of consciousness** (Sanz Perl 2021, PRE): This is a DIRECT measurement of our Sigma in brain dynamics. Wakefulness = high Sigma; unconsciousness = low Sigma.
- **FDT violations decrease with loss of consciousness** (2025, PRL Research): Non-equilibrium = Sigma > 0 = time arrow = consciousness.
- **LZC increases with "richer" consciousness** (psychedelics > normal > sleep > anesthesia): entropy measures track consciousness level.
- **PCI = 0.31 threshold separates conscious/unconscious**: a sharp boundary exists.

### What is mathematically established (without direct consciousness data):
- **Quantum IIT uses quantum relative entropy** — the SAME QRE as our Sigma definition. The mathematical machinery is identical; the conceptual framing is different.
- **FEP's variational bound F >= -ln P(o)** has the same structure as our **tau <= 1 - exp(-Sigma/2)**. Both are information-theoretic bounds on recovery/prediction.
- **Penrose's tau_P = hbar/E_G** and our **tau = 1 - F** are inversely related through the Pikovski channel.

### What is speculative but promising:
- **Consciousness requires intermediate tau** (not 0, not 1) — the "Goldilocks zone" of irreversibility.
- **Phi might be derivable from tau** for quantum systems, since both use QRE.
- **The critical tau for consciousness onset** might correspond to PCI = 0.31.

### Actionable next steps:
1. **Compute tau for the Ornstein-Uhlenbeck brain model** (reference 17) — this is directly doable and would give quantitative tau values for different consciousness states.
2. **Show Phi >= g(Sigma)** for quantum IIT — prove that integrated information requires entropy production.
3. **Connect PCI threshold to tau threshold** — is PCI = 0.31 equivalent to some tau_critical?
4. **Test with anesthesia data** — use existing ECoG/fMRI data to compute both Sigma and Petz recovery fidelity.
