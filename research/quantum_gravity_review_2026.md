# Quantum Gravity Review (2024--2026)
## Relevant to the $\Sigma = 2\ln Q$ Framework

**Compiled**: 2026-03-23
**Purpose**: Survey of recent developments in quantum gravity with direct or indirect connections to the $\Sigma = D(\rho_{\text{spacetime}} \| \rho_{\text{matter}})$ programme.

---

## 1. Information-Theoretic Approaches to Gravity

### 1.1 QRE $\Rightarrow$ Einstein Equations (CRITICAL -- Direct Connection)

**Dorau & Much (2025 PRL)**
"From Quantum Relative Entropy to the Semiclassical Einstein Equations"
arXiv: [2510.24491](https://arxiv.org/abs/2510.24491) | PRL 136, 091602 (2026)

Using modular theory (Araki-Uhlmann relative entropy), they show: the relative entropy between the vacuum state and coherent excitations of a scalar QFT on a bifurcate Killing horizon equals the energy flux across the horizon. Assuming Bekenstein-Hawking $S = A/4G$, the semiclassical Einstein equations follow automatically. This is the rigorous QFT generalisation of Jacobson's 1995 thermodynamic derivation.

**Relevance to $\Sigma$**: This is the closest independent validation of our core claim that $\Sigma = D(\rho \| \sigma)$ generates gravitational dynamics. Dorau-Much use the *same* mathematical object (quantum relative entropy) and arrive at Einstein's equations. Our $\Sigma = 2\ln Q$ extends this by providing the explicit form of $\Sigma$ in terms of the Khronon lapse.

---

**Bianconi (2024--2025)**
"Gravity from Entropy"
arXiv: [2408.14391](https://arxiv.org/abs/2408.14391) | Phys. Rev. D 111, 066001 (2025)

Proposes an *entropic action* given by the quantum relative entropy between the spacetime metric (treated as a quantum operator / effective density matrix) and the metric induced by matter fields (Dirac-Kahler formalism). Modified Einstein equations emerge, reducing to standard GR at low coupling. Key bonus: an emergent small positive cosmological constant from a "G-field" (Lagrange multiplier), which she suggests as a dark matter candidate.

Follow-up: "The Quantum Relative Entropy of the Schwarzschild Black Hole and the Area Law" arXiv: [2501.09491](https://arxiv.org/abs/2501.09491) (2025) derives the area law for the QRE of the Schwarzschild solution.

**Relevance to $\Sigma$**: Bianconi's entropic action $S_{\text{entropic}} = D(\hat{g}_{\mu\nu} \| \hat{g}^{\text{matter}}_{\mu\nu})$ is structurally identical to our $\Sigma = D(\rho_{\text{spacetime}} \| \rho_{\text{matter}})$. The emergent cosmological constant and G-field dark matter are parallel to our $\mu_0$ and $\Omega_{\text{DM}}$ predictions. This is the closest existing framework to ours -- differences are in the specific parametrisation and the Khronon structure.

---

### 1.2 Holographic Entanglement and Emergent Spacetime

**Takayanagi (2025)**
"Emergent Holographic Spacetime from Quantum Information"
arXiv: [2506.06595](https://arxiv.org/abs/2506.06595) | PRL 134, 240001 (2025)

Vision essay proposing pseudo-entropy and timelike entanglement as tools for understanding how *time* emerges in holography. Extends holographic duality beyond AdS to general (cosmological) spacetimes.

**Relevance**: The timelike entanglement entropy concept maps onto our $\tau_{\text{signed}}$ analysis. In our framework, $\tau > 0$ (spacelike separation) vs $\tau < 0$ (timelike/retrocausal) determines the arrow of time. Takayanagi's complex-valued pseudo-entropy may encode the same physics.

---

**Timelike Entanglement Entropy Programme (2023--2025)**
- "Pseudoentropy in dS/CFT and Timelike Entanglement Entropy" arXiv: [2210.09457](https://arxiv.org/abs/2210.09457) | PRL 130, 031601 (2023)
- "Timelike entanglement entropy in non-conformal theories" JHEP 07 (2024) 243
- "Holographic timelike entanglement across dimensions" JHEP 11 (2025) 100
- "Black hole singularity and timelike entanglement" JHEP 10 (2024) 182
- "Timelike Entanglement Entropy in Higher Curvature Gravity" arXiv: [2509.04181](https://arxiv.org/abs/2509.04181) (2025)

The imaginary part of pseudo-entropy implies the *emergence of time* in dS/CFT. The Ryu-Takayanagi formula is being extended to Lorentzian (timelike) settings, with complex-valued entanglement entropy.

**Relevance**: Our framework predicts $\Sigma \to$ complex when time reversal is involved (retrocausal sector). This aligns with the timelike entanglement programme. Opportunity: show that $\Sigma = 2\ln Q$ with complex $Q$ reproduces timelike entanglement entropy.

---

**Celestial Quantum Error Correction (2025)**
"Celestial Quantum Error Correction. Part II. From qudits to celestial CFT"
JHEP 06 (2025) 121

Extends holographic QEC to asymptotically flat spacetimes via the celestial CFT framework. Vacuum degeneracies and IR divergences are encoded as correctable errors.

**Relevance**: Our Petz recovery map is the QEC decoder. This paper shows how QEC works beyond AdS, which is where our framework operates (asymptotically flat / cosmological backgrounds).

---

### 1.3 Modular Flow and Spacetime Emergence

**Trejo-Calderon (2025)**
"Modular Channels, Thermal Filtering and the Spectral Emergence of Spacetime"
arXiv: [2504.20457](https://arxiv.org/abs/2504.20457)

Proposes the Modular Channels Flow Correspondence (MCFC): the flow of entangled information across causal horizons is encoded in the spectral filtering structure of quantum channels. Derives Einstein's equations from requiring local validity of the first law of entanglement (reinterpreted as a Clausius relation for modular flow) across Rindler horizons. Gravity emerges as the thermodynamic response of modular information flow.

**Relevance**: This is a *channel-theoretic* derivation of Einstein's equations, directly paralleling our Paper 2 approach (gravity as a quantum channel with Kraus operators). The "spectral filtering" is related to our thermal Petz recovery structure. The Page curve emerges from an informational phase transition at Page time.

---

### 1.4 Petz Recovery Map: Recent Advances

**Bai, Buscemi & Scarani (2024--2025)**
"Quantum Bayes' rule and Petz transpose map from the minimum change principle"
arXiv: [2410.00319](https://arxiv.org/abs/2410.00319)

Shows that the Petz transpose map emerges uniquely from the "minimum change principle" when change maximises fidelity. This is the operational foundation for our use of the Petz map as retrodiction.

**Relevance**: Directly supports our Paper 1 claim that Petz recovery = retrodiction = Bayesian updating. The minimum change principle is new motivation for why the Petz map is the correct recovery channel.

---

**Optimality Condition for the Petz Map (2024)**
arXiv: [2410.23622](https://arxiv.org/abs/2410.23622) | PRL 134, 200602 (2025)

Establishes necessary and sufficient conditions for when the Petz map is optimal.

**Relevance**: Identifies exactly when $\exp(-\Sigma/2)$ bound is saturated, which corresponds to our "zero-entropy environment" condition where retrodiction is perfect.

---

**Petz map for long-range entangled states (2024)**
arXiv: [2408.00857](https://arxiv.org/abs/2408.00857)

Studies the Petz recovery fidelity as a diagnostic of quantum phases of matter, including long-range entangled (topologically ordered) states.

**Experimental realisation on NMR (2025)**
arXiv: [2508.08998](https://arxiv.org/abs/2508.08998)

First experimental implementation of the Petz recovery map on a nuclear magnetic resonance quantum processor.

**Relevance**: Experimental verification that the Petz map works in practice. A stepping stone toward testing our gravitational Petz recovery predictions in the lab.

---

### 1.5 Entropic Dynamics

**Caticha (2025)**
"Entropic Dynamics approach to Quantum Electrodynamics"
arXiv: [2511.19238](https://arxiv.org/abs/2511.19238)

Derives QED (including Maxwell equations and the Coulomb potential) from entropic inference on a statistical manifold, without any underlying action principle. The dynamics is driven by entropy subject to appropriate constraints.

**Relevance**: Parallel programme to ours. Both derive physical laws from information-theoretic principles. Key difference: Caticha uses Shannon/Fisher information geometry; we use quantum relative entropy and Petz recovery. Potential bridge: show both approaches converge in appropriate limits.

---

## 2. Experimental Progress

### 2.1 Nanoparticle Matter-Wave Interferometry (MAJOR 2026 Result)

**Pedalino, Ramirez-Galindo, Ferstl et al. (2026)**
"Probing quantum mechanics with nanoparticle matter-wave interferometry"
Nature 649, 866 (2026)

Demonstrated quantum interference of sodium nanoparticles with >7,000 atoms at masses >170 kDa. The nanoparticles propagate in a Schrodinger cat state with macroscopicity $\mu = 15.5$, surpassing previous experiments by an order of magnitude.

**Relevance**: This is the frontier of testing quantum superposition at increasingly large masses. Our framework predicts that gravitational decoherence (via $\Sigma$) should kick in at some mass scale. The current 170 kDa result sets a *lower bound* on where gravitational decoherence effects become important. Opportunity: compute the predicted $\tau$ (retrodiction parameter) for 170 kDa nanoparticles using $\Sigma = 2\ln Q$ and compare with the observed coherence.

---

### 2.2 Tabletop Quantum Gravity: BMV Experiment

**Bose-Marletto-Vedral (BMV) Programme (2024--2026)**

- "The Bose-Marletto-Vedral experiment with nanodiamond interferometers" arXiv: [2410.19601](https://arxiv.org/abs/2410.19601) (2024)
- "Absence of gravitationally induced entanglement in certain semi-classical theories" arXiv: [2510.20991](https://arxiv.org/abs/2510.20991) (2025)
- "Evolution of tripartite entanglement in QGEM with quantum decoherence" Scientific Reports (2026)

Status: The BMV experiment aims to detect gravity-mediated entanglement between two massive particles. Nanodiamond interferometers with NV centres are the leading implementation. Bose estimates >10 years before experimental realisation. Recent 2025 work shows that entanglement generation from *local classical* gravity models (not quantum) is also possible, complicating interpretation.

**Relevance**: Our framework predicts specific decoherence rates via $\Sigma$. If BMV observes entanglement, $\Sigma$ should quantify its magnitude. If it fails to observe entanglement, our $\exp(-\Sigma)$ suppression could explain why.

---

### 2.3 Single Graviton Detection

**Tobar, Manikandan, Beitel & Pikovski (2024)**
"Detecting single gravitons with quantum sensing"
arXiv: [2308.15440](https://arxiv.org/abs/2308.15440) | Nature Communications 15, 7229 (2024)

Proposes detecting single gravitons via quantum jumps in a cooled acoustic resonator (Weber bar) at the quantum ground state, correlated with LIGO gravitational wave detections. Requires cooling a 10-ton niobium bar to ~1 mK.

**Schutzhold (2025)**
"Stimulated Emission or Absorption of Gravitons by Light"
arXiv: [2502.10221](https://arxiv.org/abs/2502.10221) | PRL 135, 171501 (2025)

Proposes an "optical Weber bar" using extended laser interferometers (Mach-Zehnder/Sagnac) to observe graviton emission/absorption via phase shifts. Claims present-day technology may suffice with ~10^6 reflections creating ~10^6 km optical path in a 1 km setup.

**Relevance**: Single graviton detection = observing the quantum channel of gravity at the single-quantum level. In our framework, each graviton absorption/emission is a Kraus operator acting on the gravitational channel. The detection rate is controlled by $\Sigma$.

---

### 2.4 LIGO Quantum Noise and Beyond

**LIGO A+ (2024)**
"Squeezing the quantum noise of a gravitational-wave detector below the standard quantum limit"
arXiv: [2404.14569](https://arxiv.org/abs/2404.14569) | Science (2024)

LIGO achieved broadband sensitivity below the Standard Quantum Limit (SQL) by up to 3 dB using frequency-dependent squeezing. Detects ~60% more mergers than before.

**Relevance**: As LIGO pushes below the SQL, it enters the regime where quantum gravitational noise could in principle be distinguished from ordinary quantum noise. Our $\Sigma$ framework predicts specific spectral signatures of gravitational decoherence that could appear as excess noise at particular frequencies.

---

### 2.5 IceCube Neutrino Decoherence

**IceCube Collaboration (2024)**
"Search for neutrino decoherence from quantum gravity"
Nature Physics (2024)

Set the world's strongest constraints on neutrino-quantum gravity decoherence, over a million times stronger than previous bounds. No decoherence detected.

**Relevance**: This constrains gravitational decoherence models. Our framework must be consistent with these null results: for neutrinos, $\Sigma$ should be extremely small (consistent with $\Sigma \propto m^2/M_P^2$ scaling).

---

### 2.6 Quantum Spacetime and Galaxy Rotation

**Koch et al. (2025)**
"Geodesics in quantum gravity"
arXiv: [2510.00117](https://arxiv.org/abs/2510.00117) | Phys. Rev. D 112, 084056 (2025) -- TU Wien

Derives "q-desic" equations: quantum-corrected geodesics in quantised spacetime. Key result: deviations from classical geodesics are negligible at Solar System scales (~$10^{-35}$ m) but become significant at galactic scales (~$10^{21}$ m) due to the cosmological constant. Suggests potential explanation for galaxy rotation anomalies.

**Relevance**: This is an *independent* route to modified dynamics at galactic scales from quantum gravity, paralleling our Paper 3 prediction. Key difference: they quantise the metric directly; we derive the modification from $\Sigma = 2\ln Q$ through the Khronon field. Same phenomenology, different mechanism. Opportunity: show that our $\Sigma$-induced force at galactic scales reproduces the q-desic deviations.

---

## 3. Direct Connections to the $\Sigma$ Framework

### 3.1 Running $G$ from QFT (CRITICAL -- Validates Paper 3)

**Kumar et al. (2025)**
"Marginal IR running of Gravity as a Natural Explanation for Dark Matter"
arXiv: [2509.05246](https://arxiv.org/abs/2509.05246) | Physics Letters B (2026)

From the EFT framework, the marginal case $\eta = 1$ (anomalous dimension) yields a logarithmic gravitational potential $\Phi \sim \ln(r)$ and a $1/r$ force law at large distances, recovering Newtonian gravity at short scales. Fits galactic rotation curves from the Sofue S-sample with a single crossover scale.

**Relevance**: This is the QFT first-principles derivation that validates our Paper 3 structure. Our $\Sigma = 2\ln(1+\delta)$ in the FRW limit, when combined with $\mu(k) = k$ running, produces the same logarithmic potential. Kumar provides the independent QFT justification for $\eta = 1$.

---

### 3.2 Khronon-Tensor Theory (CRITICAL -- Same Field Content)

**Blanchet & Skordis (2024)**
"Relativistic Khronon Theory in agreement with Modified Newtonian Dynamics and Large-Scale Cosmology"
arXiv: [2404.06584](https://arxiv.org/abs/2404.06584) | JCAP 11 (2024) 040

**Blanchet & Skordis (2025)**
"Khronon-Tensor theory reproducing MOND and the cosmological model"
arXiv: [2507.00912](https://arxiv.org/abs/2507.00912)

A scalar-tensor theory based on the Khronon field $\tau$ with DBI kinetic structure $K(Q)$. Reproduces MOND at galactic scales, GR + $\Lambda$ in strong field, and standard cosmological model (GDM) at linear cosmological scales.

**Relevance**: This is *our* Khronon field. Blanchet-Skordis provide the relativistic completion that our Papers 2-4 use. The DBI kinetic function $K(Q)$ is directly related to our $\Sigma = 2\ln Q$. Their $\omega = 0$ scalar mode has $c_s^2 = 0$, which we proved in Paper 3 via $K'(Q_0) = 0$ at ghost condensation.

---

### 3.3 Kabel et al. -- Quantum Reference Frames and Spacetime Superposition

**Kabel, de la Hamette, Apadula et al. (2024)**
"Identification is Pointless: Quantum Coordinates, Localisation of Events, and the Quantum Hole Argument"
arXiv: [2402.10267](https://arxiv.org/abs/2402.10267) | Nat. Commun. Phys. (2025)

Shows that "sameness" and "difference" of configurations in superposition are quantum-reference-frame dependent. Applied to semi-classical spacetimes in superposition: event localisation is frame-dependent. Generalises Einstein's hole argument to quantum gravity.

**Relevance**: Validates our ontological position (Paper 1b) that measurement outcomes depend on the observer's reference frame. In our framework, $\tau$ (the Petz recovery parameter) is observer-dependent when spacetimes are in superposition, consistent with Kabel et al.

---

### 3.4 Jacobson's Programme: Updates

**Thermodynamic Gravity with Non-Extensive Horizon Entropy (2025)**
arXiv: [2602.20430](https://arxiv.org/abs/2602.20430)

Extends Jacobson's thermodynamic gravity by replacing Bekenstein-Hawking entropy with a stretched-exponential group-entropy functional through a deformed Legendre structure.

**Relevance**: Shows the Jacobson programme is alive and being generalised. Our $\Sigma$ provides a specific choice of entropy functional (quantum relative entropy) that unifies with Petz recovery.

---

## 4. Black Hole Information Problem

### 4.1 Island Formula and Page Curve

**Entanglement islands for Kerr black holes (2024)**
Phys. Rev. D 110, 066012 (2024)

First application of the island rule to *rotating* (Kerr) black holes, demonstrating the Page curve for non-extremal Kerr.

**Islands in Kerr-Newman Black Holes (2025)**
arXiv: [2510.24006](https://arxiv.org/abs/2510.24006)

Extends island prescription to charged rotating black holes.

**LQG resolution of information paradox (2025)**
Loop quantum gravity black holes yield quantum extremal surfaces that suppress late-time entropy growth and preserve unitarity, with island boundaries expanding in tandem with interior volume.

**Relevance**: The island formula $S_{\text{gen}} = \min_I \left[\frac{A(\partial I)}{4G} + S_{\text{matter}}(I \cup R)\right]$ is structurally related to our $\Sigma$. The "island" is the region where Petz recovery succeeds; outside, it fails ($\tau \to 1$). The Page time = the time when $\Sigma$ transitions from growing to shrinking, i.e., when retrodiction becomes possible.

---

### 4.2 Hawking Radiation as Quantum Channel

**Tunnelling approach to unitarity (2025)**
arXiv: [2502.09924](https://arxiv.org/abs/2502.09924)

Shows that black hole evaporation is governed by a time-dependent Schrodinger equation that maps pure states to pure states. The tunnelling amplitude encodes the unitary channel.

**Information recovery from gravitational waves (2024)**
CERN seminar, October 2024

Proposes that a large portion of the information of the initial collapsing state may be carried by gravitational waves rather than Hawking radiation.

**Relevance**: Both results frame black hole evaporation as a quantum channel problem. Our framework provides the explicit channel ($\Sigma$-dependent Kraus operators) and the Petz recovery map as the decoder. The GW information channel is interesting: our $\Sigma_{\text{grav}}$ would quantify the information content of the gravitational wave channel.

---

### 4.3 ER=EPR: Theoretical Progress and Tests

**ER=EPR as operational theorem (2024)**
arXiv: [2410.16496](https://arxiv.org/abs/2410.16496) | Physics Letters B (2025)

In the LOCC setting, ER=EPR follows as an operational theorem: non-traversability of ER bridges follows from the no-signalling theorem.

**Testing ER=EPR with Hydrogen (2025)**
arXiv: [2512.02156](https://arxiv.org/abs/2512.02156)

If gauge flux leaks through the wormhole connecting entangled particles, the hydrogen 21 cm hyperfine transition would be suppressed. Current precision measurements severely constrain this possibility.

**Relevance**: ER=EPR connects entanglement (our $\tau$) with geometry (our $\Sigma$). Our framework makes this precise: $\Sigma = 2\ln Q$ quantifies the "size" of the ER bridge, while $\tau = 1 - 1/Q$ measures the entanglement. As $Q \to 1$ (flat space), $\Sigma \to 0$ and the bridge vanishes.

---

## 5. Opportunities for the $\Sigma = 2\ln Q$ Framework

### 5.1 Unique Contributions Our Framework Can Make

| Gap in the literature | What $\Sigma = 2\ln Q$ offers |
|---|---|
| Dorau-Much derive Einstein eqs from QRE but have no explicit $\Sigma$ | We give $\Sigma = 2\ln Q$ with $Q = 1/N_\phi$ from Khronon geometry |
| Bianconi's entropic gravity has no specific dark matter prediction | We predict $\Omega_{\text{DM}} = 0.268$ from extremal principle |
| Takayanagi's pseudo-entropy lacks a retrodiction interpretation | Our $\tau = 1 - 1/Q$ provides the operational meaning |
| Island formula lacks a channel-theoretic derivation | Our Petz recovery framework derives islands as recovery regions |
| BMV experiment needs specific decoherence predictions | $\Sigma$ gives explicit rates via $\exp(-\Sigma)$ suppression |
| Koch's q-desics lack an information-theoretic foundation | $\Sigma$ provides the entropic origin of modified geodesics |
| Kumar's running $G$ is phenomenological ($\eta = 1$ assumed) | $\Sigma$ framework derives the same from DPI + extensivity |

### 5.2 Immediate Opportunities (Next 3 Months)

1. **Paper connecting $\Sigma$ to Dorau-Much**: Show that $\Sigma = 2\ln Q$ is the explicit solution to their QRE framework on a Khronon foliation. This would place our work in the PRL lineage.

2. **Comparison with Bianconi**: Her QRE action $D(\hat{g} \| \hat{g}^{\text{matter}})$ is structurally identical to our $\Sigma$. Write a brief note showing the precise mapping and where the frameworks diverge (she has a G-field; we have a Khronon).

3. **Pseudo-entropy connection**: Compute $\Sigma = 2\ln Q$ for complex $Q$ (Lorentzian continuation) and compare with Takayanagi's timelike entanglement entropy. This could yield a *prediction* for timelike EE in terms of Khronon data.

4. **Experimental prediction for nanoparticle interferometry**: Using the Pedalino et al. (2026) setup at 170 kDa, compute the gravitational decoherence rate from $\Sigma$ and determine the mass threshold where our prediction departs from standard quantum mechanics.

### 5.3 Medium-Term Opportunities (3--12 Months)

5. **Channel-theoretic Page curve**: Use $\Sigma$-dependent Kraus operators to derive the Page curve directly, without invoking the island formula. Show that the Page time corresponds to $\Sigma = \Sigma_{\text{critical}}$ where Petz recovery transitions from failure to success.

6. **Graviton detection rate from $\Sigma$**: Compute the single-graviton absorption cross-section in the $\Sigma$ framework and compare with Pikovski's (2024) and Schutzhold's (2025) predictions.

7. **CMB predictions via CLASS**: Complete the CLASS/Boltzmann integration and produce $C_\ell$ predictions that can be compared against Planck and future CMB-S4 data.

### 5.4 Key Risks and Challenges

- **Bianconi overlap**: Her framework is very close to ours. Need to clearly distinguish our contributions (Khronon structure, explicit $Q$, Petz recovery) from hers (Dirac-Kahler matter, G-field).
- **Koch overlap**: Their q-desic galaxy rotation prediction is qualitatively similar to ours. Need to show that $\Sigma$ provides the *quantitative* prediction that q-desics currently lack.
- **IceCube constraint**: Our decoherence predictions must be consistent with the null result for neutrinos.

---

## 6. Summary Table of Key Papers

| Paper | Year | arXiv | Connection to $\Sigma$ | Priority |
|---|---|---|---|---|
| Dorau & Much, QRE $\Rightarrow$ Einstein | 2025 | 2510.24491 | Direct: same QRE, same equations | ★★★ |
| Bianconi, Gravity from Entropy | 2025 | 2408.14391 | Direct: QRE action = our $\Sigma$ | ★★★ |
| Takayanagi, Emergent Holographic Spacetime | 2025 | 2506.06595 | Indirect: pseudo-entropy $\leftrightarrow$ $\tau$ | ★★ |
| Trejo-Calderon, Modular Channels | 2025 | 2504.20457 | Direct: channel derivation of Einstein eqs | ★★★ |
| Kumar, IR Running of $G$ | 2025 | 2509.05246 | Direct: validates Paper 3 logarithmic potential | ★★★ |
| Blanchet & Skordis, Khronon-Tensor | 2024/25 | 2404.06584, 2507.00912 | Direct: same Khronon field, DBI $K(Q)$ | ★★★ |
| Koch et al., Q-desics | 2025 | 2510.00117 | Indirect: same phenomenology, different mechanism | ★★ |
| Bai-Buscemi-Scarani, Petz = Bayes | 2024 | 2410.00319 | Direct: Petz map foundation | ★★ |
| Optimality of Petz Map | 2024 | 2410.23622 | Direct: when $\exp(-\Sigma/2)$ saturates | ★★ |
| Pedalino et al., 170 kDa interference | 2026 | Nature 649, 866 | Experimental: decoherence frontier | ★★★ |
| Tobar-Pikovski, Single Graviton | 2024 | 2308.15440 | Experimental: quantum channel of gravity | ★★ |
| Schutzhold, Optical Weber Bar | 2025 | 2502.10221 | Experimental: graviton as channel operation | ★★ |
| Kabel et al., Identification is Pointless | 2024 | 2402.10267 | Ontological: observer-dependent $\tau$ | ★★ |
| IceCube, Neutrino Decoherence | 2024 | Nature Physics | Constraint: $\Sigma$ must be small for neutrinos | ★★ |
| LIGO A+, Below SQL | 2024 | 2404.14569 | Future: quantum gravity noise signatures | ★ |
| Caticha, Entropic QED | 2025 | 2511.19238 | Parallel: entropy-driven dynamics | ★ |
| ER=EPR Operational Theorem | 2024 | 2410.16496 | Theoretical: $\Sigma \leftrightarrow$ bridge size | ★ |

---

## 7. Assessment

The field is converging on the idea that **quantum relative entropy is the fundamental object connecting quantum information to gravity**. Three independent groups (Dorau-Much, Bianconi, and our framework) have arrived at this conclusion from different directions. The time is ripe for a unifying paper that shows these are all limits of a single structure.

The experimental frontier is advancing rapidly: 170 kDa nanoparticle interference (2026), single graviton proposals (2024--2025), and IceCube constraints (2024) are all setting the stage for tests of quantum gravity in the information-theoretic regime. Our $\Sigma = 2\ln Q$ framework, with its explicit predictions for decoherence rates and dark matter, is well-positioned to make testable claims that distinguish it from competing approaches.

**Bottom line**: The $\Sigma = D(\rho_{\text{spacetime}} \| \rho_{\text{matter}})$ programme is *not* isolated. It sits at the convergence of at least three major independent research streams (Dorau-Much, Bianconi, holographic QEC) and connects to all major experimental proposals. The next 12 months will be critical for establishing priority and making falsifiable predictions.
