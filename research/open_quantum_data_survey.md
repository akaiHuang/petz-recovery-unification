# Open Quantum Data Survey: Testing F >= exp(-Sigma/2)

**Date**: 2026-03-09
**Purpose**: Identify open datasets and published experimental results containing measured decoherence/fidelity data that can verify the bound F >= exp(-Sigma/2).

---

## Executive Summary

This survey identifies **16 primary data sources** across 5 categories that are potentially usable for testing the recovery fidelity bound F >= exp(-Sigma/2). The most immediately actionable sources are:

1. **Singh et al. (2025) NMR Petz recovery** -- directly measures Petz recovery fidelity for amplitude damping and phase damping channels
2. **Png et al. (2025) ion trap Petz recovery** -- simulates Petz recovery under realistic ion trap noise
3. **Google Willow (2024) surface code** -- open data on Zenodo with logical error rates at multiple code distances
4. **Micadei et al. (2019) NMR entropy reversal** -- measures entropy production and mutual information in quantum heat flow reversal
5. **Photonic Crooks QFT (2025)** -- reconstructs quasi-probability distribution of quantum entropy production

---

## Category 1: Petz Recovery Map Experimental Implementations

### 1.1 Singh et al. -- NMR Petz Recovery [HIGHEST PRIORITY]

- **Paper**: "Realizing the Petz Recovery Map on an NMR Quantum Processor"
- **arXiv**: [2508.08998](https://arxiv.org/abs/2508.08998) (August 2025)
- **Platform**: NMR quantum processor (liquid-state)
- **What was measured**:
  - Petz recovery fidelity for **phase damping channel**
  - Petz recovery fidelity for **amplitude damping channel**
  - Comparison of recovered states with theoretical predictions
- **Method**: Duality quantum computing (DQC) algorithm to implement Petz map
- **Data availability**: Results shown in paper figures; check arXiv supplementary for raw data
- **Can we compute F >= exp(-Sigma/2)?**: **YES** -- The channel parameters (damping rate) directly give Sigma, and they measure recovery fidelity F. This is the most direct test.
- **Action items**:
  - Download paper and extract fidelity vs. damping parameter data from figures
  - Compute Sigma = -ln(eta) for each channel parameter
  - Check if F_measured >= exp(-Sigma/2) holds for all data points

### 1.2 Png et al. -- Ion Trap Petz Recovery

- **Paper**: "Petz recovery maps of single-qubit decoherence channels in an ion trap quantum processor"
- **arXiv**: [2504.20399](https://arxiv.org/abs/2504.20399) (April 2025)
- **Published**: Phys. Rev. A 112, 022613 (August 2025)
- **Platform**: Ion trap (simulation with realistic noise)
- **What was measured**:
  - Petz recovery circuits for depolarizing, dephasing, and amplitude damping channels
  - Recovery error vs. decoherence level
  - Effect of prior knowledge precision on recovery fidelity
  - Realistic noise from residual spin-motion coupling
- **Data availability**: Simulation results in paper; circuits specified explicitly
- **Can we compute F >= exp(-Sigma/2)?**: **YES** -- Channel parameters give Sigma, simulated recovery fidelity gives F
- **Action items**:
  - Extract recovery error (1-F) vs. decoherence level from paper figures
  - Compute Sigma for each channel instance
  - Verify bound

### 1.3 Chen, Song, Scarani -- Optical Petz Recovery

- **Paper**: "Recovery of optical losses with the Petz recovery map"
- **arXiv**: [2511.05941](https://arxiv.org/abs/2511.05941) (November 2025)
- **Platform**: Optical/photonic (theoretical with explicit implementations)
- **What was measured**:
  - Petz recovery fidelity for single-mode bosonic loss channel
  - Comparison with trivial recovery (replace with belief state)
  - Near-optimality of Petz map demonstrated
- **Data availability**: Analytical results with specific numerical examples
- **Can we compute F >= exp(-Sigma/2)?**: **YES** -- Loss parameter eta gives Sigma = -ln(eta), and they compute recovery fidelity
- **Note**: Gaussian loss channels are well-characterized; Sigma is directly computable

### 1.4 Song, Kwon, Scarani -- Tabletop Time-Reversibility

- **Paper**: "Exact and approximate conditions of tabletop reversibility: when is Petz recovery cost-free?"
- **arXiv**: [2510.26895](https://arxiv.org/abs/2510.26895) (October 2025)
- **Platform**: Theoretical framework with explicit conditions
- **What was measured**: Conditions under which Petz recovery requires no additional resources
- **Can we compute F >= exp(-Sigma/2)?**: Indirectly -- provides conditions for F = 1 (Sigma = 0 case)
- **Relevance**: Confirms the tau = 0 limit of the bound

---

## Category 2: Quantum Error Correction Data (F vs. Error Rate)

### 2.1 Google Willow -- Below-Threshold Surface Code [HIGH PRIORITY]

- **Paper**: "Quantum error correction below the surface code threshold"
- **arXiv**: [2408.13687](https://arxiv.org/abs/2408.13687) (August 2024)
- **Published**: Nature (December 2024)
- **Platform**: Superconducting (105-qubit Willow processor)
- **Open data**: **Zenodo: [10.5281/zenodo.13273331](https://doi.org/10.5281/zenodo.13273331)**
- **What was measured**:
  - Logical error rate per cycle: 0.143% +/- 0.003% (distance-7)
  - Error suppression factor: Lambda = 2.14 +/- 0.02 per distance step
  - Physical qubit T1 = 68 us, T2_CPMG = 89 us
  - Below-threshold performance for d=3, d=5, d=7
  - Real-time decoder with 63 us latency
- **Can we compute F >= exp(-Sigma/2)?**:
  - **PARTIALLY** -- Need to map logical error rate to recovery fidelity F and physical error rate to entropy production Sigma
  - F_logical = 1 - p_logical, where p_logical = 0.00143/cycle
  - Sigma requires careful definition: entropy produced per QEC cycle
  - The mapping is non-trivial but feasible: Sigma_per_cycle ~ (n-k) * physical_error_rate
- **Action items**:
  - Download HDF5 data from Zenodo
  - Extract per-cycle syndrome data and logical error rates at d=3, d=5, d=7
  - Define appropriate Sigma (e.g., syndrome entropy or physical entropy production)
  - Test whether F_logical >= exp(-Sigma/2) is saturated or has a gap

### 2.2 AlphaQubit -- ML Decoder on Sycamore [HIGH PRIORITY]

- **Paper**: "Learning high-accuracy error decoding for quantum processors"
- **Published**: Nature (November 2024), DOI: [10.1038/s41586-024-08148-8](https://www.nature.com/articles/s41586-024-08148-8)
- **Platform**: Superconducting (Sycamore processor, Google)
- **What was measured**:
  - Logical error rate: 2.8 x 10^-3 at distance-9
  - 6% fewer errors than tensor network decoders
  - 30% fewer errors than correlated matching
  - Performance on real Sycamore data for d=3 and d=5
  - Performance on simulated data for d up to 11
  - Error correction accuracy: 98.5% vs. 93% for best traditional decoders
- **Data availability**: Supplementary data likely with Nature paper; check for code/data repository
- **Can we compute F >= exp(-Sigma/2)?**:
  - **YES, with interpretation** -- The decoder hierarchy (MWPM < Correlated Matching < AlphaQubit) maps to a recovery hierarchy
  - Each decoder is a different recovery map; F_decoder can be compared with Sigma_physical
  - AlphaQubit as ML decoder = empirical approximation to optimal recovery ~ Petz map
- **Key insight for Paper 1**: AlphaQubit's decoder hierarchy validates the retrodiction hierarchy prediction

### 2.3 Harvard/QuEra -- Neutral Atom QEC

- **Paper**: "A fault-tolerant neutral-atom architecture for universal quantum computation"
- **Published**: Nature (November 2025)
- **Platform**: Neutral atoms (up to 448 Rb atoms)
- **What was measured**:
  - Below-threshold surface code performance
  - Transversal logical gates
  - Teleportation-based universality
  - Logical two-qubit gate fidelity > physical two-qubit fidelity
  - Logical error rate: a few percent (lower than physical 2Q error rate)
- **Data availability**: Supplementary in Nature paper
- **Can we compute F >= exp(-Sigma/2)?**: Similar mapping as Google Willow

### 2.4 Quantinuum -- Trapped Ion QEC

- **Platform**: Trapped ion (H-series, up to 56 qubits)
- **Key results (2024-2025)**:
  - Magic state infidelity: 5.1 x 10^-4 (28 physical ions)
  - Prepared state infidelity: 7 x 10^-5
  - Logical gate error: <= 2.3 x 10^-4
  - Physical 2Q gate fidelity: 99.921%
  - Quantum volume: > 2,000,000
- **Data availability**: Published results; no known open dataset
- **Can we compute F >= exp(-Sigma/2)?**: Yes with appropriate Sigma definition

---

## Category 3: Fluctuation Theorem Experiments (Direct Sigma Measurement)

### 3.1 Photonic Generalized Crooks QFT [HIGH PRIORITY]

- **Paper**: "Experimental demonstration of generalized quantum fluctuation theorems in the presence of coherence"
- **arXiv**: [2506.00524](https://arxiv.org/abs/2506.00524) (May 2025)
- **Published**: Science Advances (2025)
- **Platform**: Photonic (single-photon polarization qubit)
- **What was measured**:
  - Quasi-probability distribution of quantum entropy production
  - Verification of generalized Crooks QFT for both covariant and incovariant channels
  - Imaginary components of quantum entropy production from coherence
  - Forward and time-reversed process quasi-probabilities
- **Data availability**: Supplementary data with Science Advances paper
- **Can we compute F >= exp(-Sigma/2)?**:
  - **YES** -- They directly measure Sigma distributions; if they also characterize the channel, we can compute F_Petz
  - The Crooks ratio P_fwd/P_rev = exp(Sigma) is directly the exponential appearing in our bound
- **Action items**:
  - Extract Sigma distribution data
  - For each channel characterized, compute the Petz recovery fidelity bound exp(-Sigma/2)
  - Compare with any reported fidelity measures

### 3.2 Cheng et al. -- Crooks FT in Single Nuclear Spin

- **Paper**: "Experimental test of the Crooks fluctuation theorem in a single nuclear spin"
- **arXiv**: [2401.17655](https://arxiv.org/abs/2401.17655) (January 2024)
- **Published**: Phys. Rev. A 109, L020401 (2024)
- **Platform**: Diamond NV center (single nuclear spin)
- **What was measured**:
  - Two-point work measurement protocol
  - Crooks FT verified at different process speeds and temperatures
  - High-fidelity single-shot readout of nuclear spin
- **Data availability**: Results in paper figures
- **Can we compute F >= exp(-Sigma/2)?**: Yes -- work distribution gives Sigma distribution

### 3.3 NV Center Coherence and Entropy Production

- **Paper**: "Experimental signature of initial quantum coherence on entropy production"
- **Published**: npj Quantum Information (September 2023)
- **Platform**: Diamond NV center (electronic spin)
- **What was measured**:
  - Entropy production with and without initial quantum coherence
  - Generalized fluctuation theorem tracking coherence effects
  - Tight bound for average heat exchange
  - Quantum contribution to irreversibility quantified
- **Data availability**: Published figures with error bars from Monte Carlo sampling
- **Can we compute F >= exp(-Sigma/2)?**: Yes -- they measure Sigma directly; need to reconstruct channel to compute F_Petz

### 3.4 Micadei et al. -- Quantum Heat Flow Reversal

- **Paper**: "Reversing the direction of heat flow using quantum correlations"
- **arXiv**: [1711.03323](https://arxiv.org/abs/1711.03323)
- **Published**: Nature Communications 10, 2456 (2019)
- **Platform**: NMR (two coupled nuclear spins)
- **What was measured**:
  - Spontaneous energy flow from cold to hot system (enabled by entanglement)
  - Mutual information change Delta I(A:B)
  - Entropy production trajectory
  - Trade-off between correlations and entropy
- **Data availability**: Figures with numerical simulations and error bars
- **Can we compute F >= exp(-Sigma/2)?**:
  - **YES** -- This is a Sigma < 0 experiment. The negative entropy production means F should be close to 1.
  - Directly tests the tau_signed < 0 regime from Paper 1
  - The mutual information decrease compensates entropy production
- **Key significance**: Tests the **Sigma < 0 regime** where time arrow reverses locally

### 3.5 PRX Quantum -- Conditional Entropy Production on IBM Quantum

- **Paper**: "Conditional entropy production and quantum fluctuation theorem of dissipative information: Theory and experiments"
- **arXiv**: [2105.06419](https://arxiv.org/abs/2105.06419)
- **Published**: PRX Quantum 3, 030315 (2022)
- **Platform**: IBM quantum computers (cloud access)
- **What was measured**:
  - Conditional and unconditional quantum entropy production
  - Quantum dissipative information (correlation dissipation)
  - Quantum fluctuation theorem verified experimentally
  - Relationship: Sigma_conditional >= Sigma_unconditional + dissipative_info
- **Data availability**: Run on IBM Quantum; results in paper
- **Can we compute F >= exp(-Sigma/2)?**:
  - **YES** -- They measure Sigma on actual quantum hardware
  - The conditional vs. unconditional distinction maps to different recovery scenarios
  - Dissipative information = correlation-induced extra cost

---

## Category 4: IBM Quantum Calibration and Process Tomography Data

### 4.1 IBM Quantum Calibration Data [CONTINUOUSLY AVAILABLE]

- **Source**: IBM Quantum Platform (cloud.ibm.com)
- **Platform**: Superconducting (Eagle, Heron processors)
- **What is measured** (updated every ~24 hours):
  - T1 (energy relaxation time) per qubit
  - T2 (dephasing time, Hahn echo) per qubit
  - Single-qubit gate error (SX gate, from randomized benchmarking)
  - Two-qubit gate error (CX/ECR gates)
  - Readout assignment error
  - Median values across all qubits
- **Data availability**: **OPEN** via Qiskit API (`backend.properties()`)
- **Can we compute F >= exp(-Sigma/2)?**:
  - **YES** with modeling: T1 gives amplitude damping rate, T2 gives dephasing rate
  - Sigma = -ln(exp(-t/T1)) = t/T1 for amplitude damping
  - F_Petz can be computed analytically for these channels
  - Can do this for every IBM backend, every day
- **Action items**:
  - Write Qiskit script to pull calibration data for all available backends
  - For each qubit: compute Sigma from T1, T2, gate error rates
  - Compute theoretical F_Petz for amplitude damping / dephasing channels
  - Verify F >= exp(-Sigma/2) across all qubits
  - Generate scatter plot of F vs. exp(-Sigma/2)

### 4.2 Qiskit Process Tomography Framework

- **Source**: qiskit-experiments library
- **What it provides**:
  - Full quantum process tomography (QPT) for arbitrary circuits
  - State fidelity computation from reconstructed density matrices
  - ProcessTomography experiment class
- **Tools available**:
  - `qiskit_experiments.library.tomography.ProcessTomography`
  - `qiskit_experiments.library.tomography.StateTomography`
  - Gradient-descent QPT: [github.com/quantshah/gd-qpt](https://github.com/quantshah/gd-qpt)
  - Variational channel fidelity: [github.com/iitis/variational_channel_fidelity](https://github.com/iitis/variational_channel_fidelity)
- **Can we compute F >= exp(-Sigma/2)?**: Yes -- run QPT, extract Choi matrix, compute both F and Sigma
- **Note**: IEEE DataPort has a collection of adaptive QST measurement data (1500 entries)

---

## Category 5: Open Databases and Repositories

### 5.1 Awesome Quantum Computing Experiments

- **GitHub**: [francois-marie/awesome-quantum-computing-experiments](https://github.com/francois-marie/awesome-quantum-computing-experiments)
- **Associated paper**: arXiv [2507.03678](https://arxiv.org/abs/2507.03678)
- **Content**: Comprehensive database of QEC experiments across platforms
  - Repetition codes, color codes, surface codes
  - Physical qubit performance benchmarks
  - Error rates across different platforms and codes
  - Auto-generated figures (PNG, PDF, interactive JS)
- **Usefulness**: Cross-platform comparison of F vs. error rates

### 5.2 Zenodo Open Data Repositories

| Dataset | Zenodo DOI | Content |
|---------|-----------|---------|
| Google Willow surface code | [10.5281/zenodo.13273331](https://doi.org/10.5281/zenodo.13273331) | HDF5 files with syndrome data, measurement fidelity |
| Real-time QEC superconducting | [zenodo.org/records/15364358](https://zenodo.org/records/15364358) | Raw experimental data, confusion matrices |
| RL control of QEC | [10.5281/zenodo.17566522](https://doi.org/10.5281/zenodo.17566522) | Training data and error correction results |
| Qudits beyond break-even | [10.5281/zenodo.15009817](https://doi.org/10.5281/zenodo.15009817) | Qudit QEC data |
| Measurement-free logical QC | [zenodo.org/records/17375920](https://zenodo.org/records/17375920) | Circuits, simulation code, data |

### 5.3 Google Quantum AI Open Source

- **GitHub**: [github.com/quantumlib](https://github.com/quantumlib)
- **Cirq**: Framework for NISQ algorithms
- **Stim**: Fast Clifford circuit simulator for QEC
- **Note**: Sycamore/Willow calibration data sometimes included in paper supplements

---

## Priority Action Plan

### Tier 1: Immediate (can verify bound within days)

| # | Source | Why | Effort |
|---|--------|-----|--------|
| 1 | **Singh NMR Petz** (2508.08998) | Directly measures Petz F for known channels | Extract data from figures |
| 2 | **Png ion trap Petz** (2504.20399) | Recovery error vs. decoherence, explicit circuits | Extract data from figures |
| 3 | **IBM Quantum calibration** | T1/T2 -> Sigma, compute F_Petz analytically | Write Qiskit script |
| 4 | **Optical Petz** (2511.05941) | Analytical F for loss channel, Sigma = -ln(eta) | Read paper, compute |

### Tier 2: Moderate effort (1-2 weeks)

| # | Source | Why | Effort |
|---|--------|-----|--------|
| 5 | **Google Willow** Zenodo data | Open HDF5, need to define Sigma for QEC cycle | Download + analysis script |
| 6 | **Photonic Crooks QFT** (2506.00524) | Directly measures Sigma distribution | Extract from Science Advances |
| 7 | **AlphaQubit** Nature paper | Decoder hierarchy = recovery hierarchy | Need supplementary data |
| 8 | **PRX Quantum IBM** (2105.06419) | Sigma measured on IBM hardware | Extract from paper |

### Tier 3: Requires reanalysis (weeks)

| # | Source | Why | Effort |
|---|--------|-----|--------|
| 9 | **Micadei NMR** heat reversal | Tests Sigma < 0 regime (!) | Reanalyze published data |
| 10 | **NV center entropy** (npj QI 2023) | Direct Sigma measurement, need channel reconstruction | Significant reanalysis |
| 11 | **Harvard neutral atom** (Nature 2025) | Below-threshold, different platform | Need supplementary data |
| 12 | **Quantinuum** magic state | Lowest infidelity, but no Sigma measured directly | Need to model Sigma |

---

## Key Observations for Paper 1

### What makes a good test of F >= exp(-Sigma/2)?

A dataset is ideal if it provides **both**:
1. **Recovery fidelity F** (or logical fidelity, or process fidelity after recovery)
2. **Entropy production Sigma** (or enough information to compute it: channel parameters, noise rates, T1/T2)

### Gap in the literature

**No existing paper directly tests the bound F >= exp(-Sigma/2)**. The Petz recovery community measures F but doesn't compute Sigma. The fluctuation theorem community measures Sigma but doesn't compute F_recovery. **This is exactly the gap Paper 1 fills.**

### Most compelling demonstration strategy

1. **Analytical**: For amplitude damping channel with parameter gamma:
   - Sigma = -ln(1-gamma)
   - F_Petz can be computed analytically
   - Show F_Petz >= exp(-Sigma/2) = (1-gamma)^(1/2) across all gamma in [0,1]

2. **Experimental re-analysis**: Take Singh (NMR) and Png (ion trap) data:
   - Their measured F vs. their channel's Sigma
   - Plot F_measured, F_Petz_theoretical, and exp(-Sigma/2) on the same graph
   - Show the bound is respected and quantify the gap

3. **Large-scale statistical**: IBM Quantum calibration across hundreds of qubits:
   - Compute Sigma from T1, T2 for every qubit on every backend
   - Compute corresponding F_Petz bounds
   - Show the bound holds universally across 1000+ qubit instances

---

## Relevant Theoretical Papers (for context)

| Paper | arXiv | Key result |
|-------|-------|------------|
| Junge-Renner-Sutter-Wilde-Winter | [1509.07127](https://arxiv.org/abs/1509.07127) | Universal recovery maps, F >= exp(-Delta(S)) |
| Buscemi-Scarani-Bai (2025) | PRL (Aug 2025) | Quantum Bayes rule = Petz map from first principles |
| Optimality of Petz map | [2410.23622](https://arxiv.org/abs/2410.23622) | Necessary and sufficient conditions for Petz optimality |
| Triple trade-off | PRX Quantum 3, 020318 | Thermodynamic constraints on QEC and information gain |
| Thermodynamic limits on FTQC | [2411.12805](https://arxiv.org/abs/2411.12805) | Dynamical phase transition in QEC thermal management |
| Approximate recoverability & DPI | [2309.02074](https://arxiv.org/abs/2309.02074) | Disproves 2015 conjecture on sandwiched Renyi recovery |

---

## Summary Statistics

- **Total papers surveyed**: ~40
- **Papers with usable fidelity data**: 12
- **Papers with usable entropy production data**: 6
- **Papers with BOTH F and Sigma**: 3 (Singh NMR, Png ion trap, optical Petz)
- **Open datasets on Zenodo/GitHub**: 5+
- **Platforms covered**: NMR, ion trap, superconducting, photonic, NV center, neutral atom

---

*Survey conducted 2026-03-09. Sources verified via web search.*
