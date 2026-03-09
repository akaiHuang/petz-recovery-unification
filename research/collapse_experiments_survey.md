# Experimental Survey: Wavefunction Collapse, Decoherence Dynamics, and Quantum-to-Classical Transition

**Date**: 2026-03-09
**Purpose**: Quantitative experimental data relevant to the tau framework (tau = 1 - F) and Petz recovery predictions
**Scope**: Papers from 2021-2026 with emphasis on 2024-2026

---

## Table of Contents

1. [Gravitational Decoherence Experiments](#1-gravitational-decoherence-experiments)
2. [Collapse Model Tests (Diosi-Penrose, CSL, GRW)](#2-collapse-model-tests)
3. [Decoherence Time Measurements](#3-decoherence-time-measurements)
4. [Quantum-to-Classical Transition Experiments](#4-quantum-to-classical-transition-experiments)
5. [Non-Markovian Dynamics / Information Backflow](#5-non-markovian-dynamics--information-backflow)
6. [Entropy Production in Quantum Systems](#6-entropy-production-in-quantum-systems)
7. [Petz Recovery Map Implementations](#7-petz-recovery-map-implementations)
8. [Connection to tau Framework](#8-connection-to-tau-framework)
9. [Summary: What Survives and What Is Ruled Out](#9-summary-what-survives-and-what-is-ruled-out)

---

## 1. Gravitational Decoherence Experiments

### 1.1 Pikovski Gravitational Time Dilation Decoherence

**Theory** (Pikovski et al., Nature Physics 11, 668, 2015):
- Gravitational time dilation induces universal decoherence of spatial superpositions
- Decoherence rate: Gamma_grav ~ (Delta_E / hbar) * (g * Delta_x / c^2)
- For micron-scale objects on Earth: already sufficient to decohere

**Experimental status (as of 2026)**:
- **No direct experimental test completed**
- Hu et al. (Quantum Sci. Technol. 2024): Analyzed feasibility of detecting quantum delocalization effects on relativistic time dilation in optical clocks
- Debski et al. (CQG 2024): Theoretical work on universality of quantum time dilation
- Pikovski group received 2024 grant for graviton detection research

**tau connection**: Pikovski channel K_0 = sqrt((1+p)/2) I, K_1 = sqrt((1-p)/2) sigma_z gives explicit CPTP map with tau_int = 2*tau_sp*(1 - tau_sp)

### 1.2 Aspelmeyer Group (Vienna) -- Optomechanical Tests

**Naber, Nimmrichter, Romero-Isart & Aspelmeyer (PNAS 121, 2024)**:
- **Platform**: Levitated 100-nm silica nanoparticle
- **Mass**: ~6 x 10^8 amu (50 nm radius, density 1850 kg/m^3)
- **Proposed superposition size**: Several nanometers, within ~2.1 ms
- **Fringe spacing**: 5 nm
- **Coherence length**: ~6.6 nm (713x zero-point motion)
- **Environmental conditions**: T = 300 K, P = 10^-10 mbar
- **Decoherence budget**:
  - Black-body radiation: dominant in free-evolution steps
  - Photon scattering: dominant during trapping steps
  - Gas collisions: P < 10^-10 mbar gives 90% completion probability
- **For coherent splitting (100 nm)**: Requires T < 50 K, P < 10^-13 mbar, coherence time > 1 s
- **Status**: Feasibility demonstrated; experiment not yet performed
- **Data availability**: Theoretical proposal with detailed noise modeling

**Aspelmeyer Group 2024 Additional Publications**:
- Ultrastrong linear optomechanical interaction (PRR 6, L042025, Oct 2024)
- Stabilizing nanoparticles in intensity minimum (Opt. Express 32, Nov 2024)
- Erasure beyond Landauer's bound using levitated optomechanics (PRR, 2024)

### 1.3 Nanodiamond Interferometry Progress

**arXiv:2508.15625 (2025)**:
- Fast neutralization of levitated charged nanodiamonds using focused UV laser
- Neutralization rate faster by > 1 order of magnitude vs state of art
- Critical step toward spatial interferometry with nanodiamonds
- Neutralization essential to avoid decoherence from stray electric fields

### 1.4 Diosi-Penrose Gravitationally Induced Entanglement

**Tilloy (PRD 111, L121101, June 2025)**:
- **Key finding**: The DP model of *classical* gravity can still generate entanglement between two separate particles
- **Condition**: Entanglement iff separation < d_max proportional to R_0 (free parameter)
- **Duration**: GIE can survive for > 1 day for reasonable experimental parameters
- **Implication**: GIE detection alone is NOT sufficient to validate quantum gravity
- **tau connection**: This constrains how GIE experiments should be interpreted within the tau framework

---

## 2. Collapse Model Tests

### 2.1 Diosi-Penrose (DP) Model

#### 2.1.1 Gran Sasso Underground Test (Donadi et al., Nature Physics, 2021)

**Platform**: High-Purity Germanium detector at LNGS underground lab
**Method**: Search for spontaneous X-ray/gamma-ray emission from protons in Ge nuclei

**Quantitative results**:
- **Ruled out**: Parameter-free version of Diosi model (original Diosi prediction with R_0 = 0)
- **Lower bound**: R_0 > 4 x 10^-10 m (from radiation experiments)
- **Heating-based lower bound**: 2 orders of magnitude weaker than radiation, but more robust to non-Markovian generalizations
- **Penrose version**: Still viable (has different parameter structure)

#### 2.1.2 Effectiveness of DP Collapse (Figurato et al., NJP 26, 113004, Nov 2024)

**Quantitative bounds on R_0**:
- **Lower bound** (from experiments): R_0 > 4 x 10^-10 m
- **Upper bound** (from requiring effective macroscopic collapse): R_0 < 10^-4 m
- **Test system**: Graphene plate, side length 25 um, ~2 x 10^10 carbon atoms
- **Key finding**: Not all macroscopic systems collapse effectively in DP model; effectiveness depends on geometry

**Allowed window**: 4 x 10^-10 m < R_0 < 10^-4 m (6 orders of magnitude)

#### 2.1.3 Gravitational Collapse Times for Molecules (Barbatti, PCCP, 2024)

**Quantitative predictions** (Diosi-Penrose model):

| System | Mass | Gravitational Self-Energy (J) | Collapse Time |
|--------|------|-------------------------------|---------------|
| NH_3 (isomerization) | 17 amu | 1.1 x 10^-51 | 10^17 s (~10^9 yr) |
| C_20 (fullerene) | 240 amu | -- | 1.3 x 10^15 s |
| C_70 | 840 amu | -- | 2.4 x 10^14 s |
| C_720 | 8640 amu | -- | 8.4 x 10^12 s |
| Neuron activation (Na+) | ~10^6 atoms | -- | 10^10 s |
| Pollen grain | 10^-11 kg | -- | 10^-7 s |
| 1 kg Al pointer (1 cm sep.) | 1 kg | -- | 6 x 10^-21 s |
| 1 kg homogeneous C | 1 kg | -- | 10^-27 s |

**Proposed experiment**: Osmium crystal, mass 5.6 x 10^-18 kg, radius ~0.04 um
- Predicted collapse time: ~1334 s
- Particle velocity: 10^-3 m/s
- Talbot length: ~1.3 m

### 2.2 CSL (Continuous Spontaneous Localization) Model

#### 2.2.1 LISA Pathfinder Rotational Bounds (arXiv:2501.08971, PRA 111, L020203, Feb 2025)

**Platform**: LISA Pathfinder space mission angular motion data
**Test masses**: Gold-platinum alloy cubes, L = 4.6 cm, m = 1.928 kg each

**Quantitative results**:
- **Minimum measured torque noise**: S_tau^exp = 5.7 x 10^-34 N^2 m^2 / Hz at f = 3 x 10^-3 Hz
- **Improvement**: Factor ~2 improvement over previous translational bounds
- **Range of improvement**: r_C from ~10^-5.5 m to ~10^-3.5 m
- **Key finding**: Rotational measurements stronger than translational for LISA geometry
- **CSL force/torque ratio**: alpha_CSL ~ L^2/6

**Theoretical CSL parameter targets**:
- GRW original: lambda = 10^-16 s^-1, r_C = 10^-7 m
- Adler enhanced: lambda ~ 10^-8 +/- 2 s^-1 (r_C = 10^-7 m) or lambda ~ 10^-6 +/- 2 s^-1 (r_C = 10^-6 m)

#### 2.2.2 Gran Sasso X-ray Bounds (2021-2022)

- Strongest bounds on CSL from spontaneous radiation emission
- For r_C < 10^-6 m: improved by > 1 order of magnitude
- Factor 13 gain over previous similar searches
- Nuclear contribution grows as Z^2 (vs Z for electrons) for E ~ 10-10^5 keV

#### 2.2.3 Astrophysical Bounds (arXiv:2406.04463, 2024)

- Revisiting astrophysical bounds on CSL models
- Complements laboratory constraints in different parameter regions

### 2.3 Current CSL Exclusion Summary

**What is ruled out**:
- Adler's enhanced CSL range (lambda ~ 10^-8 s^-1, r_C = 10^-7 m) is under severe pressure from X-ray and LISA Pathfinder data
- Original GRW values (lambda = 10^-16 s^-1) remain allowed but give negligible effects

**What survives**:
- CSL with lambda < 10^-10 s^-1 at r_C = 10^-7 m
- DP model with R_0 in range [4 x 10^-10, 10^-4] m
- Penrose version (not parameter-free Diosi) still viable after Gran Sasso

---

## 3. Decoherence Time Measurements

### 3.1 Record: 1400-Second Schrodinger Cat State (USTC, Nature Photonics, 2024)

**Platform**: Ytterbium-173 atoms (spin-5/2) in decoherence-free subspace
**Atom count**: ~10,000 atoms
**Temperature**: Few thousandths of a degree above absolute zero

**Quantitative results**:
- **Coherence time**: 1.4(1) x 10^3 s (= 23 minutes 20 seconds)
- **Measurement sensitivity**: Heisenberg-limit enhanced by 15(2) dB over coherent spin states
- **Mechanism**: Decoherence-free subspace isolates from magnetic noise
- **Significance**: Record coherence for multi-particle entangled state

**tau connection**: In DFS, tau -> 0 by construction (engineered closed system). The 1400 s lifetime measures residual tau from imperfect isolation.

### 3.2 Google Willow Superconducting Qubits (Nature, Dec 2024)

**Platform**: 105-qubit superconducting processor

**Quantitative results**:
- **T_1** (energy relaxation): Mean 68 us (range up to 100 us)
- **T_2,CPMG** (dephasing): Mean 89 us
- **Single-qubit gate fidelity**: 99.87%
- **Two-qubit gate fidelity**: 99.65%
- **Logical error suppression**: Lambda = 2.14 +/- 0.02 per code distance increase of 2
- **Logical memory**: Exceeded best physical qubit lifetime by factor 2.4 +/- 0.3
- **Error per cycle** (distance-7 code): 0.143% +/- 0.003%
- **5x improvement** in T_1 over Sycamore generation (20 -> 100 us)

**tau connection**: QEC = engineered local tau -> 0. The Lambda = 2.14 quantifies how effectively each additional code distance reduces tau. Physical tau from T_1 = 68 us sets the baseline; logical tau is exponentially suppressed.

### 3.3 Levitated Nanoparticle Ground-State Cooling

**Two-mode cooling (Nature Physics, 2024)**:
- Occupation numbers: n_x = 0.83, n_y = 0.81 (near ground state)
- Decoherence dominated by photon recoil

**Wavepacket expansion (2025)**:
- Coherence length expanded from ~10 pm to ~70 pm by weakening/restoring optical trap
- Enables future matter-wave interference with larger particles

**Mechanical quantum memory (Nature Physics, 2025)**:
- T_2 extended from 64 us to 1 ms using two-pulse dynamical decoupling

---

## 4. Quantum-to-Classical Transition Experiments

### 4.1 Nanoparticle Schrodinger Cats (Arndt Group, Nature, 2025)

**Platform**: Sodium nanoparticle clusters, far-field matter-wave interferometry

**Quantitative results**:
- **Mass**: 172 kDa mean (range 143-197 kDa), > 7000 atoms
- **Diameter**: ~8 nm
- **Velocity**: ~160 m/s (spread ~10 m/s)
- **Fringe visibility**: V = 0.10 +/- 0.01 (primary), V = 0.66 +/- 0.09 (400 kDa - 1 MDa subset)
- **Delocalization**: > 10x particle diameter
- **Macroscopicity measure**: mu = 15.5 (surpasses previous record by order of magnitude)
- **Coherence length**: 160-350 fm
- **Grating period**: d = 133 nm
- **Interferometer spacing**: L = 0.983 m
- **Vacuum**: ~9 x 10^-9 mbar
- **CSL exclusion**: "Most stringent exclusion limit for generic macrorealistic modifications"

**tau connection**: tau = 1 - V gives tau ~ 0.90 for 172 kDa particles. The question is whether this tau is entirely environmental or includes collapse contributions. The experiment constrains collapse-induced tau.

### 4.2 Quantum Darwinism in Superconducting Circuits (Science Advances, 2025)

**Platform**: Superconducting quantum processor, 12 qubits

**Quantitative results**:
- **T_1**: ~130 us
- **Single-qubit gate fidelity**: 0.9997
- **Two-qubit CZ gate fidelity**: 0.998
- **System**: 2-qubit system + 10 environment qubits
- **Mutual information**: I(S:F) reaches H_S ~ 1 bit for small environment fragments
- **At m=4, N=4, theta=pi/2**: I(S:F) ~ 1.83
- **Quantum discord**: D -> 0 in classical plateau (redundant encoding regime)
- **Discord peak**: D -> H_S when m = N = 4 at theta = pi/2

**tau connection**: Quantum Darwinism describes HOW tau > 0 arises from system-environment correlations. The discord -> 0 plateau corresponds to complete decoherence (tau -> 1 for off-diagonal elements). Redundant encoding = information about pointer states encoded in environment = tau approaching 1 for non-pointer superpositions.

### 4.3 Experimental Blueprint: Decoherence vs Objective Collapse (arXiv:2512.02838, Dec 2025)

**Proposed platform**: Levitated nanosphere in optical trap

**Quantitative parameters**:
- **Particle radius**: 50 nm
- **Mass**: 1.0 x 10^-17 kg (~6 x 10^9 amu)
- **Trap frequency**: Omega/2pi = 10^5 Hz
- **Zero-point width**: x_zpf ~ 7 x 10^-13 m
- **Temperature**: T < 5 K (ambient), T_int < 20 K (internal)
- **Pressure**: P < 10^-15 mbar
- **Laser power**: P_laser < 5 mW

**Decoherence budget**:
- Environmental: Gamma_env ~ 3 x 10^-3 s^-1
- CSL-induced: Gamma_CSL ~ 3 x 10^-4 s^-1
- Gas collisions: Gamma_gas < 10^-3 s^-1
- Blackbody radiation: < 10^-4 s^-1
- Trap photon recoil: Gamma_trap < 10^-3 s^-1

**Detection targets**:
- Superposition size: Delta_x = 100-200 nm
- Coherence time: tau_coh > 10 s
- Critical mass window: 10^7 - 10^8 amu
- Minimum detectable rate: Gamma_min = 0.05 s^-1
- Significance: > 5 sigma
- Measurement uncertainty: sigma_i / Gamma_i = 0.05 (5%)

**Key CSL predictions**:
- CSL test value: lambda_CSL = 10^-21 s^-1 at r_C ~ 100 nm
- CSL predicts: saturation of decoherence rate with superposition size, quadratic scaling with mass
- Environmental decoherence: linear scaling with both

**tau connection**: This experiment is specifically designed to distinguish tau_collapse (objective collapse contribution) from tau_env (environmental decoherence). In the tau framework: tau_total = tau_env + tau_collapse. The mass-scaling test directly probes whether tau has a collapse component beyond environmental Sigma.

---

## 5. Non-Markovian Dynamics / Information Backflow

### 5.1 High-Dimensional Non-Markovianity in Optical Fibers (Quantum, Aug 2024)

**Platform**: 4-core multicore optical fiber (MCF), d = 4 qudits

**Quantitative results**:
- **Photon number**: mu = 0.15 mean photons/pulse
- **Single photon purity**: 95.7% of non-vacuum pulses
- **State preparation fidelity**: > 99%
- **Qudit reconstruction fidelity**: 98.5%
- **Repetition rate**: 2 MHz
- **Phase stabilization interval**: 100 ms
- **Total pulses per evolution stage**: 200,000
- **Relative entropy of coherence** (maximally coherent): C_r,max = log(4)
- **Coherence minimum**: At t = 0.5 (when p_3 = p_0 = 1/2)
- **Non-Markovianity demonstrated**: Via Quantum Vault protocol

**Information backflow**: Non-monotonic decay of quantum capacities observed, with complete recovery possible in optimal scenarios.

**tau_signed connection**: This is a controlled demonstration of tau_signed < 0 (information flowing from environment back to system). The non-monotonic coherence = local Sigma < 0 episodes. The complete recovery in optimal scenarios = tau returning to 0.

### 5.2 Causal vs Non-Causal Information Revivals (arXiv:2405.05326, revised 2025)

- Provides operational condition to witness genuine backflow (not just non-causal revivals)
- Distinguishes true Sigma < 0 from apparent information recovery
- Important for correctly interpreting tau_signed measurements

### 5.3 Non-Markovianity and Teleportation (2025)

- Connection established between non-Markovianity and quantum teleportation protocols
- Information flows bidirectionally between system and environment during teleportation
- Extends tau_signed framework to communication protocols

---

## 6. Entropy Production in Quantum Systems

### 6.1 Quantum Dot Entropy Production (Shen et al., Nature Physics, 2026)

**Platform**: Quantum dot under stepwise laser driving
**Method**: Machine learning + physics-based model to extract trajectory-level entropy production

**Key achievements**:
- **First measurement** of full dynamics of trajectory-level entropy production in non-stationary, non-Markovian system
- Combined time-dependent driving with non-Markovian memory effects
- Measured entropy produced by quantum dot stochastically blinking under control protocol
- Quantifies reversibility of microscopic processes

**tau connection**: Direct measurement of Sigma at trajectory level in a non-Markovian system. The tau framework predicts F >= exp(-Sigma/2); this experiment provides the Sigma values needed to test the bound. Non-Markovian effects correspond to episodes of tau_signed < 0.

### 6.2 Landauer's Principle in Quantum Many-Body Regime (Aimet et al., Nature Physics, 2025)

**Platform**: Quantum field simulator using ultracold Bose gases (TU Wien + FU Berlin + UBC + U Crete + U Pavia)

**Key results**:
- Measured entropy production in Bose gas following global mass quench (massive -> massless Klein-Gordon)
- Resolved contributions from:
  - Quantum mutual information changes
  - Relative entropy differences of environment
  - Correlations vs dissipation
- Validated generalized Landauer principle: Delta S_sys >= -beta * Q + Delta I(S:E) + Delta D(rho_E || sigma_E)
- Tracked temporal evolution via dynamical tomographic reconstruction

**tau connection**: This is the most direct experimental test of the Retrodiction Landauer Principle framework. The decomposition into mutual information + relative entropy differences mirrors exactly the structure of the tau bound: tau <= 1 - exp(-Sigma/2) where Sigma includes both thermal and informational contributions.

### 6.3 NV Center Coherence and Entropy Production (npj Quantum Information, 2023)

**Platform**: Nitrogen-Vacancy center in diamond, driven qubit

**Key results**:
- Observed significant increase of irreversibility (entropy production) arising from quantum coherence in initial state
- Verified both detailed and integral fluctuation theorem
- Showed: measuring coherence contribution provides tight bound for average heat exchange
- Demonstrated that coherence in initial state is a genuine thermodynamic resource

**tau connection**: Coherence in initial state -> higher Sigma -> higher tau. This quantitatively confirms that quantum coherence (which Petz recovery must preserve) directly affects entropy production.

### 6.4 Iterative Quantum Information Transfer (Suleymanzade et al., PRX 15, 031054, Aug 2025)

**Platform**: Silicon-Vacancy (SiV) center in diamond nanocavity (Harvard + U Tokyo)

**Quantitative results**:
- **Temperature**: 200 mK
- **T_2***: ~4 us (dephasing)
- **T_1**: > 1 s (relaxation, negligible depolarization)
- **Feedback latency**: 50 ns (much shorter than coherence time)
- **Protocol**: 10 successive measurement-feedback cycles per trial
- **Markovian trials**: 7,000 repetitions
- **k-th order Markovian trials**: 10,000 repetitions (k = 1,2,3,4)
- **Key result**: Achieved Sigma < 0 (entropy reduction) through ground/excited state stabilization
- **Second law**: Sigma >= -i_QCT (quantum-classical transfer entropy)
- **Efficiency**: eta = -Sigma / i_QCT increases with non-Markovianity
- **Non-Markovian advantage**: k >= 2 feedback achieves Sigma < -i_FB^1 (beyond Markovian bound)
- **Fluctuation theorem**: <exp(-Sigma - i_QCT)> ~ 1 verified throughout all cycles

**tau connection**: This experiment explicitly achieves Sigma < 0 locally through measurement and feedback -- precisely the mechanism described by tau_signed < 0. The non-Markovian advantage (memory-dependent feedback outperforming Markovian) directly demonstrates that information backflow (tau_signed < 0 episodes) enables greater entropy reduction. The generalized fluctuation theorem verified here is the single-trajectory version of the Crooks relation underlying the tau bound.

### 6.5 Crooks Fluctuation Theorem in Single Nuclear Spin (PRA 109, L020401, 2024)

**Platform**: Single nuclear spin in diamond (NMR-style)

**Key results**:
- High-fidelity single-shot readout of nuclear spin
- Implemented two-point work measurement protocol
- Crooks theorem verified for different speeds of nonequilibrium processes
- Valid under various effective temperatures

**tau connection**: Direct verification of the Crooks relation P_forward/P_reverse = exp(Sigma) in a quantum system. This is the foundation of the tau bound derivation.

### 6.6 Generalized Quantum Fluctuation Theorems with Coherence (PMC, 2025)

**Platform**: Quantum photonic system

**Key results**:
- Validated generalized Crooks QFT for processes with nonclassical features
- Quasi-probabilistic entropy production distributions characterized
- Coherence induces imaginary components of quantum entropy production -> phase factors in QFT
- Ratio between quasi-probabilities of time-forward and time-reversed processes obeys generalized Crooks QFT

**tau connection**: Extends the Crooks foundation of the tau framework to include quantum coherence effects. The imaginary entropy production from coherence is a genuinely quantum feature not captured by classical tau.

---

## 7. Petz Recovery Map Implementations

### 7.1 NMR Implementation (Singh et al., arXiv:2508.08998, 2025)

**Platform**: 3-qubit NMR processor (1H, 19F, 13C nuclei), 600 MHz Bruker spectrometer, ~300 K

**Quantitative results**:
- **Pseudopure state fidelity**: 0.9799 +/- 0.0010
- **Channels tested**: Amplitude damping (AD) and phase damping (PD)
- **Reference states**: epsilon = 0.2, 0.5, 0.8
- **AD channel**: |0> maintained F ~ 1.0; |1> experienced significant decoherence with recovery improving at higher epsilon
- **PD channel**: |+> and |-> recovery highly dependent on reference state; epsilon = 0.5 (maximally mixed) gives NO recovery
- **Key finding**: Petz recovery performance strongly depends on reference state overlap with input state
- **Statistical precision**: Error bars smaller than marker size

### 7.2 Ion Trap Implementation (Pino et al., PRA 112, 022613, 2025)

**Platform**: Ion trap quantum processor

**Quantitative results**:
- **Ancilla qubits**: 1 (Kraus rank 2) or 2 (Kraus rank > 2)
- **CNOT gates**: 3 (Kraus rank 2) or 20 (Kraus rank > 2)
- **Recovery error threshold**: < 0.01
- **Channels**: Depolarizing, dephasing, amplitude damping
- **Noise modeling**: Includes residual spin-motion coupling
- **Status**: Circuit constructions validated; full experimental results in paper

### 7.3 Connection to tau Framework

The Petz recovery experiments directly test: F(rho, R_Petz(N(rho))) >= exp(-Sigma/2)

For the NMR experiment:
- When reference state matches input: high F -> low tau (successful retrodiction)
- When reference state mismatches: low F -> high tau (retrodiction fails)
- epsilon = 0.5 for PD channel: F ~ 0 -> tau ~ 1 (complete collapse, no retrodiction possible)

This is exactly the Paper 1b prediction: collapse = tau -> 1 = Petz recovery failure.

---

## 8. Connection to tau Framework

### 8.1 Direct Tests of F >= exp(-Sigma/2)

| Experiment | Platform | Sigma Measured? | F Measured? | tau Testable? |
|------------|----------|-----------------|-------------|---------------|
| Shen et al. 2026 | Quantum dot | Yes (trajectory) | No | Sigma available |
| Aimet et al. 2025 | Bose gas | Yes (decomposed) | No | Sigma decomposition matches tau structure |
| Suleymanzade et al. 2025 | SiV diamond | Yes (+ fluctuation thm) | No | Sigma < 0 achieved (tau_signed < 0) |
| Singh NMR 2025 | NMR Petz | Indirectly | Yes (recovery F) | **Directly testable** |
| Pino ion trap 2025 | Ion trap Petz | Indirectly | Yes (recovery F) | **Directly testable** |
| Arndt nanoparticles 2025 | Matter-wave | Environmental Sigma | V = 0.10 (tau ~ 0.90) | Collapse vs environmental Sigma |
| USTC cat state 2024 | Yb atoms | ~0 (DFS) | ~1 (1400 s) | tau ~ 0 by design |

### 8.2 tau_signed < 0 Tests

| Experiment | Platform | Evidence for Sigma < 0 | Mechanism |
|------------|----------|------------------------|-----------|
| Suleymanzade PRX 2025 | SiV + feedback | Yes (explicit) | Measurement + non-Markovian feedback |
| Optical fiber 2024 | MCF qudits | Yes (information backflow) | Non-Markovian channel dynamics |
| Micadei 2019 (NMR) | Nuclear spins | Yes (heat from cold to hot) | Quantum correlations |
| Crooks photonic 2025 | Photons | Quasi-probability Sigma < 0 | Quantum coherence effects |

### 8.3 Mass-Scaling Predictions

The tau framework predicts:
- Environmental decoherence: Gamma_env ~ m (linear)
- CSL collapse: Gamma_CSL ~ m^2 (quadratic)
- DP collapse: Gamma_DP ~ m^2 / R_0^3 (quadratic in mass, depends on R_0)

**Current experimental constraints on mass scaling**:
- Arndt 2025: mu = 15.5 at 172 kDa -- no excess decoherence beyond environmental
- Aspelmeyer proposal: 10^8 amu target mass would probe quadratic regime
- arXiv:2512.02838 blueprint: 10^7-10^8 amu critical mass window for CSL detection

### 8.4 QEC as Engineered tau = 0

Google Willow demonstrates:
- Physical tau: from T_1 = 68 us -> tau_physical ~ some finite value per cycle
- Logical tau: suppressed by Lambda^(d/2) where Lambda = 2.14
- Distance-7 code: tau_logical reduced by factor ~(2.14)^3.5 ~ 25x relative to distance-1
- This is exactly Paper 1's prediction: QEC = engineering local Sigma -> 0

---

## 9. Summary: What Survives and What Is Ruled Out

### Collapse Models Status (March 2026)

| Model | Status | Key Constraint |
|-------|--------|----------------|
| Diosi (parameter-free, R_0 = 0) | **RULED OUT** | Gran Sasso 2021 |
| Diosi-Penrose (R_0 > 0) | Allowed: R_0 in [4x10^-10, 10^-4] m | Gran Sasso + Figurato 2024 |
| Penrose (non-Diosi version) | **Still viable** | No definitive test yet |
| CSL (Adler enhanced) | **Under severe pressure** | LISA Pathfinder + X-ray bounds |
| CSL (GRW original) | Allowed but negligible effects | lambda = 10^-16 too small to detect |
| CSL (intermediate) | Partially constrained | Depends on (lambda, r_C) region |

### Key Experimental Milestones Needed

1. **Nanoparticle interferometry at 10^8 amu**: Would definitively test CSL quadratic mass scaling (Aspelmeyer proposal feasible)
2. **Coherent splitting > 100 nm at 10^8 amu**: Would test DP model in allowed R_0 window
3. **Direct Sigma measurement + simultaneous Petz recovery F**: Would test F >= exp(-Sigma/2) bound directly
4. **Pikovski gravitational decoherence test**: Still untested; requires molecular interferometry in Earth's gravity

### Most Promising Near-Term Experiments for tau Framework

1. **Re-analysis of Singh NMR 2025 data**: Petz recovery fidelities can be directly compared to channel entropy production to test F >= exp(-Sigma/2)
2. **Re-analysis of Pino ion trap 2025 data**: Same as above
3. **Suleymanzade PRX 2025**: Already measures Sigma at trajectory level with fluctuation theorem -- adding Petz recovery step would give complete tau test
4. **Arndt nanoparticle data**: Fringe visibility V vs mass scaling can be reinterpreted as tau vs Sigma plot
5. **Aimet Bose gas 2025**: Entropy production decomposition already matches Retrodiction Landauer structure

---

## References (Chronological)

1. Pikovski, I. et al. "Universal decoherence due to gravitational time dilation." Nature Physics 11, 668 (2015). [arXiv:1311.1095](https://arxiv.org/abs/1311.1095)

2. Micadei, K. et al. "Reversing the direction of heat flow using quantum correlations." Nature Communications 10, 2456 (2019). [DOI](https://www.nature.com/articles/s41467-019-10333-7)

3. Donadi, S. & Piscicchia, K. et al. "Underground test of gravity-related wave function collapse." Nature Physics 17, 74-78 (2021). [DOI](https://www.nature.com/articles/s41567-020-1008-4)

4. Experimental signature of initial quantum coherence on entropy production. npj Quantum Information (2023). [DOI](https://www.nature.com/articles/s41534-023-00738-0)

5. Naber, N., Nimmrichter, S., Romero-Isart, O. & Aspelmeyer, M. "Fast quantum interference of a nanoparticle via optical potential control." PNAS 121 (2024). [DOI](https://www.pnas.org/doi/abs/10.1073/pnas.2306953121)

6. Barbatti, M. "Gravitationally-induced wave function collapse time for molecules." PCCP (2024). [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC11305101/)

7. Figurato, A. et al. "On the effectiveness of the collapse in the Diosi-Penrose model." New J. Phys. 26, 113004 (2024). [DOI](https://iopscience.iop.org/article/10.1088/1367-2630/ad8c77)

8. Crooks fluctuation theorem in single nuclear spin. Phys. Rev. A 109, L020401 (2024). [arXiv:2401.17655](https://arxiv.org/abs/2401.17655)

9. Non-Markovianity in High-Dimensional Open Quantum Systems using Next-generation Multicore Optical Fibers. Quantum (Aug 2024). [DOI](https://quantum-journal.org/papers/q-2024-08-12-1436/)

10. Google Quantum AI. "Quantum error correction below the surface code threshold." Nature (Dec 2024). [DOI](https://www.nature.com/articles/s41586-024-08449-y)

11. Minutes-scale Schrodinger-cat state of spin-5/2 atoms. Nature Photonics (2024). [arXiv:2410.09331](https://arxiv.org/abs/2410.09331)

12. Improved bounds on collapse models from rotational noise of LISA Pathfinder. PRA 111, L020203 (Feb 2025). [arXiv:2501.08971](https://arxiv.org/abs/2501.08971)

13. Tilloy, A. "Diosi-Penrose model of classical gravity predicts gravitationally induced entanglement." PRD 111, L121101 (June 2025). [arXiv:2411.02287](https://arxiv.org/abs/2411.02287)

14. Aimet, S. et al. "Experimentally probing Landauer's principle in the quantum many-body regime." Nature Physics (2025). [DOI](https://www.nature.com/articles/s41567-025-02930-9)

15. Zhu, Q. et al. "Observation of quantum Darwinism and the origin of classicality with superconducting circuits." Science Advances (2025). [DOI](https://www.science.org/doi/10.1126/sciadv.adx6857)

16. Suleymanzade, A. et al. "Experimentally Probing Entropy Reduction via Iterative Quantum Information Transfer." PRX 15, 031054 (Aug 2025). [arXiv:2411.06709](https://arxiv.org/abs/2411.06709)

17. Generalized quantum fluctuation theorems with coherence. PMC (2025). [DOI](https://pmc.ncbi.nlm.nih.gov/articles/PMC12124356/)

18. Singh, R. et al. "Realizing the Petz Recovery Map on an NMR Quantum Processor." arXiv:2508.08998 (2025). [arXiv](https://arxiv.org/abs/2508.08998)

19. Pino, J. et al. "Petz recovery maps of single-qubit decoherence channels in an ion trap quantum processor." PRA 112, 022613 (2025). [arXiv:2504.20399](https://arxiv.org/abs/2504.20399)

20. Pedalino, S. et al. "Probing quantum mechanics with nanoparticle matter-wave interferometry." Nature (2025). [DOI](https://www.nature.com/articles/s41586-025-09917-9)

21. Experimental Blueprint for Distinguishing Decoherence from Objective Collapse. arXiv:2512.02838 (Dec 2025).

22. Shen, Y. et al. "Non-equilibrium entropy production and information dissipation in a non-Markovian quantum dot." Nature Physics (2026). [DOI](https://www.nature.com/articles/s41567-026-03177-8)

---

*Last updated: 2026-03-09*
