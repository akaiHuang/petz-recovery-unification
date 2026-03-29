# Negative Energy Production: State of the Art (2024--2026)

**Compiled: 2026-03-28**
**Purpose: Survey of technologies for negative energy density, NEC violation, and vacuum engineering -- prerequisites for warp drives and exotic spacetimes.**

---

## Table of Contents

1. [Executive Summary](#1-executive-summary)
2. [Squeezed Vacuum & Negative Energy Density](#2-squeezed-vacuum--negative-energy-density)
3. [Quantum Energy Teleportation (QET)](#3-quantum-energy-teleportation-qet)
4. [Casimir Effect Engineering](#4-casimir-effect-engineering)
5. [Dynamic Casimir Effect](#5-dynamic-casimir-effect)
6. [Topological Casimir Effect](#6-topological-casimir-effect)
7. [Casimir Repulsion](#7-casimir-repulsion)
8. [Casimir Torque](#8-casimir-torque)
9. [Unruh Effect Detection](#9-unruh-effect-detection)
10. [Schwinger Effect](#10-schwinger-effect)
11. [Analog Gravity Experiments](#11-analog-gravity-experiments)
12. [Null Energy Condition (NEC) Violation](#12-null-energy-condition-nec-violation)
13. [Vacuum Energy Extraction & Zero-Point Energy](#13-vacuum-energy-extraction--zero-point-energy)
14. [Quantum Vacuum Thrusters (EmDrive etc.)](#14-quantum-vacuum-thrusters-emdrive-etc)
15. [Warp Drive Theory Status](#15-warp-drive-theory-status)
16. [The Key Question: Most Negative Energy Ever Measured](#16-the-key-question-most-negative-energy-ever-measured)
17. [Ford-Roman Quantum Inequalities](#17-ford-roman-quantum-inequalities)
18. [Archimedes Experiment: Weighing the Vacuum](#18-archimedes-experiment-weighing-the-vacuum)
19. [Technology Readiness Summary](#19-technology-readiness-summary)
20. [Connection to Sigma Framework](#20-connection-to-sigma-framework)

---

## 1. Executive Summary

**Bottom line:** Negative energy density is REAL, routinely produced in quantum optics labs (squeezed states, Casimir cavities), and has been measured at up to -15 dB below vacuum. However:

- The amounts are **fantastically small** (sub-femtojoule scale, nanometer-scale geometries)
- **Quantum inequalities** (Ford-Roman) fundamentally constrain how much negative energy can exist and for how long
- **No macroscopic NEC violation** has been achieved or is on the horizon
- The gap between lab-scale negative energy and warp-drive requirements is ~60+ orders of magnitude

**Most promising near-term directions:**
1. Casimir nanostructure engineering (geometry-controlled forces, 10x suppression demonstrated)
2. Magnetic-field tunable Casimir force (attraction-to-repulsion switching, Nature Physics 2024)
3. Topological insulator multilayer Casimir cavities (enhanced repulsion)
4. Quantum Energy Teleportation on quantum computers (Ikeda 2023, Purdue 2024)
5. Unruh effect detection via Josephson junctions (Hiroshima 2025, potentially within reach)

---

## 2. Squeezed Vacuum & Negative Energy Density

### Status: REAL PHYSICS (routinely produced)

**What it is:** In squeezed states of light, quantum noise in one quadrature (amplitude or phase) is reduced BELOW the vacuum level. This suppression of vacuum fluctuations creates intervals where the energy density is formally negative relative to the vacuum.

**Record squeezing levels:**
- **-15 dB** (2016, Hannover/AEI): World record, generated at 1064 nm using a doubly resonant type-I OPA below threshold. Published in PRL 117, 110801 (2016).
- **-14 dB**: Mentioned as the largest degree of squeezing using standard nonlinear crystals in other contexts.
- **5.9 dB (lossless)**: Achieved in 2025 using polarization self-rotation in atomic vapor with AI-stabilized control (PMC 2025).

**Practical application -- LIGO (2024):**
- Frequency-dependent squeezed light injection into LIGO **doubled** the gravitational wave detection rate (Science, October 2024).
- New 300m optical cavities enable squeezing across the full LIGO frequency band.
- Expected 60% more merger detections vs. previous runs.

**Key 2025 developments:**
- **Bright Squeezed Vacuum (BSV)** driving strong-field photoemission at metal needle tips (Nature Physics, November 2025). BSV has zero mean field but large intensity fluctuations.
- **Ultrafast squeezed vacuum measurement** on lithium niobate nanophotonic circuits (CLEO 2025).
- **AI-controlled ultra-stable squeezed vacuum source** maintaining stability for hours (2025).

**Negative energy density interpretation:**
- In every squeezed vacuum state, the energy density falls below zero once every cycle for any nonzero squeezing parameter.
- The -15 dB record corresponds to a variance reduction factor of ~31.6x below vacuum.
- However, quantum inequalities ensure: the more negative the energy density, the shorter its duration.

**Classification: REAL PHYSICS** -- peer-reviewed, replicated worldwide, used in operational instruments (LIGO).

---

## 3. Quantum Energy Teleportation (QET)

### Status: REAL PHYSICS (demonstrated on quantum hardware)

**What it is:** Protocol (Hotta 2008) exploiting vacuum entanglement: Alice measures locally, sends classical information to Bob, who uses it to extract energy from his local vacuum -- energy that was "hidden" in vacuum correlations.

**Key experiment -- Ikeda 2023:**
- **First demonstration on real quantum hardware**: IBM superconducting quantum computers (ibmq_lima, ibm_cairo, ibmq_jakarta).
- Published: Phys. Rev. Applied 20, 024051 (2023).
- Bob received ~1-5% of Alice's input energy.
- Results consistent with exact solution after measurement error mitigation.

**Purdue 2024 breakthrough:**
- Enhanced QET protocol with additional qubit for energy **storage**.
- Previous protocols lost teleported energy to classical devices immediately.
- Tested on IBM quantum computers.
- Published as arXiv:2409.03973 (September 2024).
- Described as "extracting and storing energy from a quasi-vacuum."

**Key physics:**
- Bob extracts energy from vacuum fluctuations -- this is negative energy creation at Alice's site compensated by positive energy at Bob's.
- Energy conservation is preserved: Bob cannot extract more than Alice input.
- No FTL signaling: requires classical communication channel.
- Only significant at short distances.

**Classification: REAL PHYSICS** -- peer-reviewed, replicated, code open-sourced on GitHub.

---

## 4. Casimir Effect Engineering

### Status: REAL PHYSICS (precision measurements, active engineering)

**The Casimir energy density between parallel perfect conductors:**

$$\rho_{\text{Casimir}} = -\frac{\pi^2 \hbar c}{720\, a^4}$$

where $a$ is plate separation.

**Numerical example:**
- At $a = 100$ nm: $\rho \approx -4.4 \times 10^{-4}$ J/m$^3$ (negative!)
- At $a = 10$ nm: $\rho \approx -4.4 \times 10^{0}$ J/m$^3$
- Force: $F/A \approx 0.013$ dyn/cm$^2$ at $a = 1\,\mu$m

**Key 2024-2025 developments:**

1. **Magnetic-field tuning of Casimir force (Nature Physics, May 2024):**
   - USTC China: Gold sphere + silica plate in water-based ferrofluids.
   - Casimir force tuned from **attractive to repulsive** by varying magnetic field.
   - First reversible, in-situ switching of Casimir force sign.
   - Implications: switchable micromechanical devices.

2. **3D nanostructure control (Nano Letters, 2025):**
   - Experimentally demonstrated geometry engineering of Casimir interaction.
   - Cylindrical pillars, holes, periodic arrays.
   - **10x suppression** of Casimir force for single pillar geometry.
   - Dramatically modified force power law through geometry.

3. **Precision measurements (2024 review):**
   - AFM-based measurements between smooth non-magnetic and magnetic metals.
   - Lateral Casimir forces between corrugated surfaces.
   - Thermal Casimir effect with graphene at 250-700 nm separation.
   - Drude model conclusively excluded at 200-700 nm; plasma model agrees with data.

4. **Gapped metal Casimir switching (January 2025):**
   - Teflon + gapped metal surface across liquid media.
   - Switch between attractive and repulsive Casimir forces.
   - Zero-frequency Casimir effect demonstrated.

**Classification: REAL PHYSICS** -- decades of replication, sub-percent agreement with theory.

---

## 5. Dynamic Casimir Effect

### Status: REAL PHYSICS (demonstrated 2011, active development)

**What it is:** Photon creation from vacuum by rapidly changing boundary conditions (accelerating mirrors). Converts vacuum fluctuations into real photon pairs.

**Landmark experiment:**
- **2011 (Chalmers, Nature):** Superconducting circuit with SQUID acting as effective mirror moving at ~5% speed of light. First observation of dynamical Casimir photons.

**2024-2025 developments:**
- Comprehensive review of DCE publications 2020-2024 published March 2025.
- Optomechanical backreaction of quantum field processes studied (2024).
- Resonator frequency modulation approach (Optics Express, 2025).
- **Mechanical DCE with low-frequency oscillator** (Phys. Rev. A, 2025) -- bringing DCE to more accessible experimental regimes.
- DCE in superconducting cavities reviewed (April 2025).
- Gravitational-wave-induced DCE studied theoretically (2026).

**Key limitation:** Photon production rate is extremely small. Even in the 2011 experiment, the "mirror" was not physical but an effective boundary in a superconducting circuit.

**Classification: REAL PHYSICS** -- peer-reviewed, but only demonstrated in superconducting circuits, not with physical mirrors.

---

## 6. Topological Casimir Effect

### Status: REAL PHYSICS (theoretical + emerging experiments)

**What it is:** Topological materials (topological insulators) modify the Casimir force through their non-trivial surface states and broken symmetries.

**Key findings:**
- **Topological insulator multilayers** can enhance Casimir repulsion (established 2016, developed further 2025).
- **Photonic topological insulator multilayers with broken time-reversal symmetry (2025):**
  - Casimir repulsive force demonstrated in multilayer systems.
  - Casimir stable equilibrium and restoring force can be controlled.
  - Both attraction and repulsion enhanced through structural optimization (more layers, thicker layers).
- **Higher-temperature topological insulators (Wurzburg, 2025):** Working at -213 C (previously required much colder).

**Why it matters for negative energy:** Topological materials offer a route to engineer repulsive (effectively "negative") Casimir forces through material properties rather than geometry alone.

**Classification: REAL PHYSICS** -- theoretical framework solid, experimental verification emerging.

---

## 7. Casimir Repulsion

### Status: REAL PHYSICS (measured)

**What it is:** Under certain conditions (material-medium combinations), the Casimir force becomes repulsive -- the plates push apart rather than attract. This is effectively creating a region with "less negative" energy than expected, or equivalently, a repulsive quantum vacuum force.

**Established results:**
- **Munday & Capasso (Nature, 2009):** First measurement of repulsive Casimir-Lifshitz force between gold and silica surfaces separated by bromobenzene liquid.

**2024-2025 developments:**
1. **Magnetic-field tunable switching** (Nature Physics, 2024): Attraction-to-repulsion transition via ferrofluids.
2. **Gapped metal materials** (2025): Teflon/gapped-metal/liquid system enables Casimir force quantum switching.
3. **Magneto-electric materials theory** (2025): Phase diagram governing sign of Casimir force based on parity and time-reversal symmetry.

**Classification: REAL PHYSICS** -- measured, replicated, controllable.

---

## 8. Casimir Torque

### Status: REAL PHYSICS (measured)

**What it is:** Rotational torque arising from quantum vacuum fluctuations between anisotropic (birefringent) materials. The vacuum tries to align the optical axes.

**Key measurements:**
- **First measurement (Munday group, Nature 2018):** Between liquid crystal and solid birefringent crystal.
- **2024 floating gold flake experiment:** Flake achieved opposing orientation (~60 degrees), then rotated under Casimir torque until full alignment.
- **Thermal effects (PRL, July 2025):** Thermal corrections to Casimir torque vs. twist angle characterized. Dielectric function perpendicular to optic axis plays crucial role.

**Classification: REAL PHYSICS** -- peer-reviewed, replicated.

---

## 9. Unruh Effect Detection

### Status: CLAIMS (proposed experiments, approaching feasibility)

**What it is:** An accelerating observer sees the vacuum as a thermal bath at temperature $T_U = \hbar a / (2\pi c k_B)$. Requires acceleration ~$10^{20}$ m/s$^2$ for ~1 K.

**Breakthrough proposals and demonstrations (2025):**

1. **Hiroshima University -- Josephson Junction approach (PRL, July 2025):**
   - Uses circular motion of metastable fluxon-antifluxon pairs in coupled annular Josephson junctions.
   - Advances in superconducting microfabrication enable extremely small radii.
   - Achieves effective Unruh temperature of **a few kelvin** -- detectable with current technology!
   - **This may be the most promising near-term route to Unruh detection.**

2. **Stockholm + IISER Mohali -- Superradiance method (December 2025):**
   - Accelerating atoms collectively emit a flash of light.
   - Unruh effect shifts the timing of this flash -- a measurable signal.

3. **Trapped-ion demonstration (arXiv, October 2025):**
   - Proof-of-principle of **timelike Unruh effect** in trapped-ion system.
   - Two-level spin coupled to vibrational field encoding Minkowski vacuum.
   - Thermal response resembling Unruh effect demonstrated.

**Classification: CLAIMS** -- theoretical proposals are peer-reviewed but full experimental confirmation of genuine Unruh effect not yet achieved. The Hiroshima Josephson junction approach is the closest to feasibility.

---

## 10. Schwinger Effect

### Status: REAL PHYSICS (analog), CLAIMS (direct)

**What it is:** Spontaneous electron-positron pair creation from vacuum in a sufficiently strong electric field ($E_{\text{crit}} \sim 1.3 \times 10^{18}$ V/m). Never directly observed -- requires fields ~1000x beyond current lasers.

**Analog demonstrations:**
1. **Graphene (Nature Physics, 2023):** Quantitative measurement of Schwinger-pair production rate in doped graphene transistors at ENS Paris. Mesoscopic Klein-Schwinger effect.
2. **Graphene superlattice (2022):** Analog process between electrons and holes at Dirac point, National Graphene Institute (Geim et al.).
3. **Superfluid helium-4 (PNAS, September 2025):**
   - **UBC (Stamp, Desrochers, Marchand):** Thin film of superfluid He-4 cooled to frictionless vacuum state.
   - Flow acts as "electric field"; vortex/anti-vortex pairs appear spontaneously.
   - Key finding: **vortex mass is NOT fixed -- it changes with motion.**
   - This hints at corrections to Schwinger's original theory.

**Direct observation status:**
- Not yet achieved. Requires $\sim 10^{18}$ V/m.
- **Dynamically Assisted Schwinger Effect (DASE):** Combining strong low-frequency + weak high-frequency fields enhances pair production by orders of magnitude. May become feasible with next-generation laser facilities (ELI, XCELS).

**Classification:** Analog demonstrations = REAL PHYSICS. Direct vacuum pair production = NOT YET ACHIEVED.

---

## 11. Analog Gravity Experiments

### Status: REAL PHYSICS (multiple platforms)

**Platforms demonstrating analog Hawking radiation:**

1. **Superconducting circuits:** Josephson junction arrays creating effective moving mirrors (original DCE experiment, Chalmers 2011).
2. **Superconducting transmon qubits:** 10-qubit chain simulating black hole. Stimulated Hawking radiation verified via state tomography (2023).
3. **BEC (Bose-Einstein Condensate):** Steinhauer group (Technion) -- multiple observations of analog Hawking radiation, entanglement verification (2016, 2019).
4. **Polariton fluids:** Proposed as ideal platform (2025).

**2025 developments:**
- Tunneling method for Hawking quanta in analogue gravity (Comptes Rendus Physique, 2025).
- Path integral methods for correlation functions near acoustic horizons in superfluids.

**Classification: REAL PHYSICS** -- multiple independent groups, multiple platforms.

---

## 12. Null Energy Condition (NEC) Violation

### Status: REAL PHYSICS (quantum), NO LAB DEMONSTRATION (macroscopic)

**Theoretical status:**
- NEC states $T_{\mu\nu} k^\mu k^\nu \geq 0$ for all null vectors $k^\mu$.
- Quantum fields CAN violate NEC (squeezed states, Casimir effect).
- **Quantum Null Energy Condition (QNEC)** provides the quantum correction: $\langle T_{kk} \rangle \geq \frac{\hbar}{2\pi} S''$, where $S''$ is second derivative of entanglement entropy.

**Observational searches (2024):**
- LIGO/Virgo O1-O3 data searched for NEC-violating signatures in stochastic gravitational wave background. **No statistically significant evidence found** (arXiv:2404.07075).
- Pulsar timing arrays explored for NEC-violating inflation signatures (JHEP, 2024).

**Key result (2025):**
- "Null Impact of the Null Energy Condition in Current Cosmology" (arXiv:2511.07526): Current cosmological data cannot distinguish NEC-violating models.

**Laboratory status:** No direct macroscopic NEC violation has been demonstrated. All known laboratory negative energy (squeezed states, Casimir) is constrained by quantum inequalities to be too small and too brief for macroscopic effects.

**Classification: REAL PHYSICS** (quantum NEC violation exists in principle); **NO MACROSCOPIC DEMONSTRATION.**

---

## 13. Vacuum Energy Extraction & Zero-Point Energy

### Status: Mostly PSEUDOSCIENCE (extraction claims), REAL PHYSICS (concept)

**What is real:**
- Zero-point energy is a genuine quantum mechanical phenomenon (Lamb shift, van der Waals forces, Aharonov-Bohm effect, electronic noise).
- Casimir force is a measurable manifestation.
- QET demonstrates energy can be extracted from vacuum correlations (but Bob never gets more than Alice put in).

**What is NOT real:**
- "Zero-point energy generators" = perpetual motion machines = **PSEUDOSCIENCE**.
- ZPE represents the LOWEST energy state. You cannot extract energy from a ground state by definition.
- Any "over-unity" device claiming to harvest ZPE violates thermodynamics.

**Legitimate research directions:**
- **Casimir batteries/actuators:** Using Casimir force to do mechanical work (energy comes from changing plate separation, not from vacuum itself).
- **Vacuum fluctuation battery (thought experiment):** Casimir force doing work on charged plate stack.
- **Vacuum fluctuation engineering (Phys.org, June 2025):** "Opens the door to engineering novel quantum materials by reshaping vacuum" -- using optical cavities to confine fluctuations and modify material properties.

**DIA report:** The US Defense Intelligence Agency published "Concepts for Extracting Energy from the Quantum Vacuum" -- this is a review of theoretical ideas, NOT evidence that extraction works.

**Classification:** ZPE as physics = REAL. ZPE extraction devices = PSEUDOSCIENCE. Casimir-based mechanical work = REAL but not "free energy."

---

## 14. Quantum Vacuum Thrusters (EmDrive etc.)

### Status: PSEUDOSCIENCE (debunked)

**EmDrive:**
- No peer-reviewed evidence supports functionality as of 2025.
- **Definitively debunked:** Tajmar group (TU Dresden, 2018) measured null result on high-precision torsion balance. Follow-ups in 2021 and 2022 confirmed null results below photon thrust levels.
- All positive results traced to experimental artifacts (thermal effects, electromagnetic interference).

**IVO Quantum Drive:**
- Hardware reported "operational in space" by late 2025 (survived launch).
- However, **no thrust measurement** has been reported.
- Surviving launch is an engineering milestone, NOT evidence of new physics.

**Classification: PSEUDOSCIENCE** -- violates conservation of momentum, all careful experiments show null results.

---

## 15. Warp Drive Theory Status

### Status: ACTIVE THEORY (no experimental path yet)

**Key developments (2021-2025):**

1. **Lentz soliton (2021):** First superluminal warp drive solutions using ONLY positive energy densities. Published in CQG. However, energy requirements are astronomical.

2. **Bobrick-Martire subluminal warp (2021):** Physical warp drives at subluminal speeds without exotic matter.

3. **Warp Factory (Applied Physics, 2024):** Open-source computational tool for analyzing warp drive spacetimes. First numerical implementation of physical warp drives.

4. **Constant-Velocity Subluminal Warp Drive (2024):** New design complying with GR at constant subluminal speed WITHOUT negative energy.

5. **Electromagnetic field manipulation (March 2025):**
   - Energy requirement reduced to ~$4.9 \times 10^6$ J for 20m warp bubble.
   - Down from ~$10^{62}$ J in original Alcubierre metric.
   - Uses gauge transformations and EM field manipulation.
   - **Eliminates exotic matter requirement.**

6. **Black hole gravitational assist (2024):** BH gravitational field reduces negative energy requirement. Potentially instrumental for microscopic warp drives.

**Current theoretical status:**
- Superluminal warp drives still require NEC violation.
- Subluminal positive-energy warp drives exist mathematically.
- Energy requirements have been reduced by ~56 orders of magnitude theoretically.
- No experimental prototype or path to one.

**Classification: REAL PHYSICS (theoretical), but far from experimental realization.**

---

## 16. The Key Question: Most Negative Energy Ever Measured

### Answer: Squeezed vacuum at -15 dB below vacuum level

**Quantitative assessment:**

| Method | Negative Energy Achieved | Scale | Year |
|--------|------------------------|-------|------|
| Squeezed vacuum (OPA) | -15 dB below vacuum ($\sim 31.6\times$ suppression) | Single-mode optical field | 2016 |
| Casimir cavity (100 nm) | $\sim -4.4 \times 10^{-4}$ J/m$^3$ | Nano-gap between plates | Continuous |
| Casimir cavity (10 nm) | $\sim -4.4$ J/m$^3$ | Extreme nano-gap | Continuous |
| QET on IBM quantum computer | $\sim 1$-$5\%$ of Alice's input | Single qubits | 2023 |
| Dynamic Casimir (SQUID) | $\sim$ few photons | Microwave frequency | 2011 |

**In absolute terms:**
- For a Casimir cavity at $a = 10$ nm with area $1\,\text{cm}^2$: total negative energy $\sim -4.4 \times 10^{-12}$ J ($\sim -4$ pJ).
- For squeezed vacuum: the negative energy density exists in bursts lasting fractions of an optical cycle ($\sim 10^{-15}$ s), with magnitude bounded by quantum inequalities.

**The gap to warp drives:**
- Alcubierre metric (original): $\sim 10^{62}$ J of negative energy.
- Even "optimized" designs: $\sim 10^{6}$-$10^{30}$ J.
- Best laboratory negative energy: $\sim 10^{-12}$ J.
- **Gap: at least 18 orders of magnitude to the most optimistic warp designs, 74 orders to the original.**

---

## 17. Ford-Roman Quantum Inequalities

### The fundamental constraint

Ford and Roman (1995-1997) proved that negative energy density is subject to an **inverse relationship** between magnitude and duration:

$$|\langle \rho_{\text{neg}} \rangle| \cdot (\Delta t)^4 \lesssim \frac{\hbar}{c}$$

**Financial analogy (Ford & Roman):** Nature allows you to "borrow" negative energy, but:
1. The loan term is limited.
2. The repayment must always exceed the amount borrowed.
3. Larger loans have shorter terms.

**Implications:**
- You CAN have arbitrarily large negative energy density, but only for correspondingly short times.
- You CAN have long-duration negative energy, but only at correspondingly tiny magnitude.
- **Wormholes and warp drives require large negative energy for sustained periods -- quantum inequalities make this extraordinarily difficult.**
- For a wormhole throat of 1 m, the required negative energy density violates quantum inequality bounds by many orders of magnitude.

**Recent development (2024):**
- "Wormhole Restrictions from Quantum Energy Inequalities" (MDPI, 2024): Updated bounds on wormhole traversability from QEIs.

---

## 18. Archimedes Experiment: Weighing the Vacuum

### Status: REAL PHYSICS (prototype operational, results pending)

**Goal:** Determine if vacuum energy gravitates -- i.e., does the Casimir vacuum inside a cavity have weight?

**Method:**
- Beam balance measures weight change of a Casimir cavity as plates become superconducting.
- Superconducting plates are more reflective -> better vacuum expulsion -> lower energy -> weight change.
- Uses layered superconductor (YBCO, $T_c = 93.5$ K) as natural multi-layer Casimir cavity.

**Location:** SarGrav Laboratory, Sos Enattos former mine, Sardinia (seismically quiet).

**Current sensitivity (2025):**
- Torque sensitivity: $\sim 7 \times 10^{-13}$ Nm/$\sqrt{\text{Hz}}$ (improved from earlier $\sim 10^{-11}$).
- Operating band: 50-150 mHz.
- Final measurement target: detect weight variation from modulated Casimir energy.

**Why it matters:** If vacuum energy gravitates normally, it deepens the cosmological constant problem. If it doesn't, something fundamental is wrong with our understanding of vacuum energy and gravity.

**Classification: REAL PHYSICS** -- INFN-funded, operational prototype, published results on noise performance.

---

## 19. Technology Readiness Summary

| Technology | TRL | Status | Negative Energy? | Scalable? |
|-----------|-----|--------|-------------------|-----------|
| Squeezed vacuum states | 9 | Operational (LIGO) | Yes (quantum scale) | No (quantum inequalities) |
| Casimir force measurement | 7-8 | Routine lab measurements | Yes (relative to vacuum) | Limited |
| Casimir nanostructure engineering | 4-5 | Lab demonstrations | Tunable force | Possibly |
| Magnetic Casimir switching | 3-4 | Published (Nature Physics) | Reversible sign change | Unknown |
| Topological Casimir repulsion | 2-3 | Theory + early experiments | Enhanced repulsion | Unknown |
| QET on quantum computers | 3-4 | Demonstrated on IBM hardware | Yes (local vacuum extraction) | Limited |
| Dynamic Casimir effect | 3-4 | SC circuits only | Photon pair creation | Very limited |
| Unruh effect detection | 2-3 | Proposals approaching feasibility | N/A (detection only) | N/A |
| Schwinger effect (direct) | 1-2 | Analogs only, direct not achieved | Would create pairs | Requires $10^{18}$ V/m |
| Warp drive | 1 | Theory only | Requires macroscopic NEC violation | Not with current physics |

**TRL scale:** 1 = basic principles, 3 = proof of concept, 5 = lab validation, 7 = system prototype, 9 = operational.

---

## 20. Connection to Sigma Framework

The negative energy landscape connects to the Sigma framework in several ways:

1. **Casimir energy density = negative Sigma contribution**: The Casimir vacuum between plates has $\rho < 0$, corresponding to a region where $\Sigma_{\text{grav}} = D(\rho_{\text{spacetime}} \| \rho_{\text{matter}})$ acquires unusual structure due to the boundary-induced vacuum modification.

2. **Quantum inequalities and Petz recovery**: Ford-Roman bounds $|\rho_{\text{neg}}| \cdot (\Delta t)^4 \lesssim \hbar/c$ resemble the Petz recovery bound $F \geq e^{-\Sigma/2}$. Both express the fundamental cost of "departing from equilibrium" -- whether it is the vacuum state or a quantum channel's fixed point.

3. **QET and retrodiction**: QET explicitly uses vacuum entanglement to extract energy -- the "retrodiction" of Alice's measurement result by Bob is precisely the information-theoretic structure that Petz recovery formalizes.

4. **NEC and Sigma positivity**: The NEC $T_{\mu\nu} k^\mu k^\nu \geq 0$ is related to the positivity of relative entropy $D(\rho \| \sigma) \geq 0$. NEC violation corresponds to situations where the "gravitational channel" acquires features not captured by the standard Petz recovery picture.

5. **Warp drives and $\tau \to 1$**: A warp drive would require sustained $\tau \to 1$ (complete departure from Petz recoverability) in a macroscopic region -- precisely what quantum inequalities forbid.

---

## Key References

### Squeezed States & LIGO
- Vahlbruch et al., PRL 117, 110801 (2016) -- 15 dB squeezing record
- LIGO frequency-dependent squeezing, Science (October 2024)

### Quantum Energy Teleportation
- Hotta, Phys. Lett. A 372, 5671 (2008) -- original QET proposal
- Ikeda, Phys. Rev. Applied 20, 024051 (2023) -- first hardware demonstration
- Rodriguez-Briones et al., arXiv:2409.03973 (2024) -- energy storage QET (Purdue)

### Casimir Engineering
- Zhao et al., Nature Physics 20, 1282 (2024) -- magnetic-field tuning
- Nano Letters (2025) -- 3D nanostructure Casimir control

### Dynamic Casimir
- Wilson et al., Nature 479, 376 (2011) -- first DCE observation
- MDPI Physics 7(2), 10 (2025) -- 55-year DCE review

### Unruh Effect
- Hiroshima University, PRL (July 2025) -- Josephson junction proposal
- arXiv:2510.24163 (October 2025) -- trapped-ion timelike Unruh

### Schwinger Effect
- Schmitt et al., Nature Physics (2023) -- graphene Schwinger analog
- Desrochers, Marchand & Stamp, PNAS (September 2025) -- superfluid He-4 analog

### Warp Drives
- Lentz, CQG 38, 075015 (2021) -- positive-energy solitons
- Alcubierre-Bobrick-Martire, Warp Factory (2024)
- Applied Physics subluminal warp drive (2024)

### Quantum Inequalities
- Ford & Roman, Phys. Rev. D 51, 4277 (1995)
- Ford & Roman, Phys. Rev. D 55, 2082 (1997)

### Archimedes Experiment
- Calloni et al., INFN (2015-2025) -- vacuum weight measurement
- SarGrav Laboratory, Sardinia -- operational prototype

---

*This survey reflects the state of knowledge as of March 2026. The field is evolving rapidly, particularly in Casimir nanostructure engineering and Unruh effect detection proposals.*
