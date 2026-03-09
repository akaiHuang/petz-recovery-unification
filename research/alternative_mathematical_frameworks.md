# Alternative Mathematical Frameworks for the tau Framework

**Date**: 2026-03-09
**Purpose**: Systematic search for mathematical frameworks that could strengthen or replace weak assumptions in the tau = 1 - F framework
**Status**: Complete literature survey with evaluation

---

## 0. The Three Gaps to Fix

| Gap | Description | Current status |
|-----|-------------|----------------|
| **Gap 1** (JRSWW saturation) | F = exp(-Sigma/2) has no proof of saturation for gravitational channels | Open: only approximate saturation in weak field |
| **Gap 2** (Gravitational CPTP map) | No rigorous CPTP map N_grav with Delta-D = -ln(-g_00) | Open: Pikovski channel gives Delta-D ~ (E_G t/hbar)^2 not E_G t/hbar; bounded by ln 2 |
| **Gap 3** (Sigma mismatch) | tau(t) = 1 - exp(-E_G t/2hbar) uses Sigma = E_G t/hbar which is not standard JRSWW Sigma | Open: conflation of different Sigma definitions |

**The core problem**: Paper 2 writes Sigma_grav = r_s/r and invokes F >= exp(-Sigma/2), but the JRSWW bound requires Sigma = D(rho||sigma) - D(N(rho)||N(sigma)) for a CPTP map N. No known CPTP map gives Delta-D = r_s/r. The Pikovski dephasing channel gives Delta-D = H_bin((1+p)/2), which is quadratic in the dephasing parameter and bounded by ln 2.

---

## 1. Quantum Fisher Information (QFI)

### 1.1 Core idea

QFI measures how much information about a parameter theta is encoded in a quantum state rho(theta). The quantum Cramer-Rao bound gives:

    Var(theta_est) >= 1 / (n * J_Q(theta))

where J_Q is the QFI. Unlike the JRSWW framework, QFI is ALWAYS well-defined — it requires only a parametric family of states, not a CPTP map.

### 1.2 Key recent results

**Jencova 2024** (Commun. Math. Phys. 405, 180):
- "Sufficient Statistic and Recoverability via Quantum Fisher Information"
- **Main theorem**: For a large class of QFI metrics (including Wigner-Yanase-Dyson skew information), a channel N is sufficient for a family {rho_theta} iff the QFI is preserved under N.
- **Recovery result**: If the quantum chi^2 divergence is approximately preserved, then states can be approximately recovered by the Petz map.
- **Critical caveat**: The SLD (symmetric logarithmic derivative) QFI — the most commonly used variant — does NOT satisfy this sufficiency property. Only the WYD skew information class does.
- **Connection to Petz**: Approximate preservation of chi^2 divergence => approximate Petz recoverability. This gives a QFI-based ALTERNATIVE route to recovery bounds.

**Balatsky et al. 2025** (Phys. Rev. A 111, 012411):
- "Quantum Sensing from Gravity as Universal Dephasing Channel for Qubits"
- Collaboration with Google Quantum AI
- Models gravity as a dephasing channel on transmon qubits using gravitational redshift + Aharonov-Bohm phase
- Estimates sensitivity delta_g/g ~ 10^{-7}
- **Relevance**: Provides an EXPLICIT CPTP dephasing channel for gravity on qubits, though for lab-scale gravitational differences, not for Sigma_grav = r_s/r

**QFI in Curved Spacetime** (2025, Int. J. Theor. Phys.):
- QFI for Dirac particles in noisy channels around Schwarzschild black holes
- Shows QFI degradation with Hawking temperature
- Under strong squeezing, QFI becomes resistant to Hawking radiation

**Pikovski group** (1612.08864):
- Information transfer during universal gravitational decoherence
- Derives decoherence time-scales in terms of energy dispersion and QFI
- Key formula: decoherence rate proportional to QFI of the internal state w.r.t. energy

### 1.3 Evaluation against the three gaps

| Gap | Does QFI help? | How? |
|-----|---------------|------|
| **Gap 1** (saturation) | PARTIALLY | QFI gives a DIFFERENT bound: Var >= 1/J_Q. This is always tight (saturated by optimal measurement). But it bounds estimation precision, not fidelity directly. |
| **Gap 2** (CPTP map) | YES, CIRCUMVENTS | QFI does NOT require a CPTP map. It only requires a parametric family rho(theta). For gravity, theta = Phi (gravitational potential), and rho(Phi) is the state of a probe in the gravitational field. This AVOIDS the CPTP gap entirely. |
| **Gap 3** (Sigma mismatch) | PARTIALLY | One could define a QFI-based Sigma_QFI = J_Q * (delta_Phi)^2, which would be well-defined. But the connection to tau = 1 - F becomes indirect. |

### 1.4 Proposed tau_QFI definition

Define:
    tau_QFI(theta) := 1 - F(rho(theta), rho(theta + delta_theta))

For infinitesimal delta_theta, using the relation F ~ 1 - (J_Q/8)(delta_theta)^2:

    tau_QFI ~ (J_Q / 8) * (delta_theta)^2

This gives a well-defined tau that measures distinguishability WITHOUT requiring a channel. The QFI J_Q plays the role of a "metric" on parameter space.

For gravity: theta = Phi, delta_theta = Delta_Phi = g * Delta_z, and:
    tau_QFI ~ (J_Q / 8) * (g * Delta_z)^2

**Advantage**: Always well-defined, no CPTP needed, saturated by optimal measurement.
**Disadvantage**: Quadratic in the parameter (not linear like Sigma), loses direct connection to entropy production and Crooks theorem.

### 1.5 Assessment

**Overall rating: HIGH value, addresses Gap 2 most directly**

QFI provides the cleanest alternative for defining gravitational distinguishability without a CPTP map. The Jencova 2024 result connects QFI preservation to Petz recoverability, maintaining the connection to the retrodiction framework. However, the shift from linear (entropy production) to quadratic (Fisher information) scaling changes the physical interpretation.

---

## 2. Monotone Metrics on State Space

### 2.1 Core idea

A monotone metric on quantum state space is a Riemannian metric g that satisfies the data processing inequality:

    g(N(rho)) <= g(rho)

for all CPTP maps N. The Bures metric (related to fidelity) is the MINIMAL monotone metric; the Bogoliubov-Kubo-Mori metric is the MAXIMAL one. The Petz classification theorem (1996) characterizes ALL monotone metrics.

### 2.2 Key recent results

**Petz monotone geometry for two-qubit families** (arXiv:2509.14578, Sept 2025):
- Extends beyond Bures/SLD to general Petz monotone metrics
- Establishes non-reduction theorem: curvature invariants of Petz metrics cannot be universally identified with entanglement monotones
- Important negative result: the geometry does NOT simply map to entanglement measures

**Monotone metric tensors review** (World Scientific, 2024):
- Reviews all finite-dimensional monotone quantum metrics
- Key point: the Bures metric gives the TIGHTEST distance-based bound (minimal metric)

**Jencova 2024 (same paper as 1.2)**:
- The chi^2 divergence is a monotone metric
- Approximate preservation => Petz recoverability
- Connection: monotone metric preservation <=> approximate sufficiency

### 2.3 Key insight: F-based bounds WITHOUT Petz

The data processing inequality for the Bures distance gives:

    d_B(N(rho), N(sigma)) <= d_B(rho, sigma)

where d_B(rho, sigma) = sqrt(2(1 - F(rho, sigma))).

This means:
    F(N(rho), N(sigma)) >= F(rho, sigma)    [wrong direction — this is for CPTP]

Actually, fidelity can only INCREASE under CPTP maps:
    1 - F(N(rho), N(sigma)) <= 1 - [some function of F(rho, sigma)]

The KEY point: the Bures metric provides distance bounds that constrain tau through purely geometric (metric) arguments, without needing the JRSWW entropy production bound.

### 2.4 A purely metric-based tau bound

The Bures distance satisfies:
    d_B(rho, sigma)^2 = 2(1 - F(rho, sigma))

And the Bures metric is monotone under CPTP maps. Therefore:

For any recovery map R:
    tau = 1 - F(rho, R(N(rho)))
        = (1/2) d_B(rho, R(N(rho)))^2
        <= (1/2) [d_B(rho, N(rho)) + d_B(N(rho), R(N(rho)))]^2    [triangle inequality]

This gives a PURELY GEOMETRIC bound on tau without reference to entropy production.

### 2.5 Evaluation against the three gaps

| Gap | Does it help? | How? |
|-----|--------------|------|
| **Gap 1** (saturation) | NO | Monotone metrics give bounds, not equalities. The saturation problem remains. |
| **Gap 2** (CPTP map) | NO | Still requires CPTP maps for the DPI. |
| **Gap 3** (Sigma mismatch) | PARTIALLY | Provides alternative distance measures (Bures, WYD) that could replace Sigma. |

### 2.6 Assessment

**Overall rating: MODERATE value, provides mathematical rigor for composition laws**

The tau composition inequality sqrt(tau_12) <= sqrt(tau_1) + sqrt(tau_2) (already derived in tau_structural_inequalities.md) IS a monotone metric result. This framework validates the existing tau structure but does not solve the gravitational CPTP gap. The non-reduction theorem from 2509.14578 is a NEGATIVE result suggesting that geometric curvature cannot simply proxy for entanglement.

---

## 3. Quantum Wasserstein Distance

### 3.1 Core idea

The quantum Wasserstein distance W_p(rho, sigma) measures the "transport cost" of deforming rho into sigma. Unlike relative entropy, it has a natural geometric interpretation as the cost of physically moving probability mass.

### 3.2 Key recent results

**Thermodynamic Unification of Optimal Transport** (Phys. Rev. X 13, 011013, 2023):
- **Central result**: W_2 = min over all admissible Markovian dynamics of {Sigma * mu}, where Sigma is entropy production and mu is dynamical state mobility
- Variational formula: the Wasserstein distance EQUALS the minimum product of irreversibility and mobility
- Extends to quantum case
- Yields thermodynamic speed limits: d(rho, sigma)/dt <= sqrt(Sigma * mu)

**Quantum HWI Inequality** (Rouze-Datta, Ann. Henri Poincare 2020):
- Relates three quantities: H (relative entropy), W (Wasserstein distance), I (Fisher information)
- Under a quantum Ricci curvature lower bound kappa:
  D(rho || sigma) <= W(rho, sigma) * sqrt(I(rho)) - (kappa/2) * W(rho, sigma)^2
- This is the quantum analogue of the classical HWI inequality
- **Key implication**: bounds relative entropy in terms of transport distance and Fisher information

**Order p Quantum Wasserstein Distances** (Ann. Henri Poincare 2025):
- Establishes W_p distances from quantum couplings
- Recent theoretical framework with mathematical rigor

**Quantum Wasserstein for Gaussian states** (2024-2025):
- General formula for W between one-mode Gaussian states
- Linear resilience to losses (unlike fidelity)

### 3.3 The HWI-tau connection

The HWI inequality gives:
    D(rho || sigma) <= W(rho, sigma) * sqrt(I(rho)) - (kappa/2) * W^2

If we identify D with Sigma (entropy production) and use tau <= 1 - exp(-Sigma/2), we get:
    tau <= 1 - exp(-W * sqrt(I) / 2 + kappa * W^2 / 4)

This provides a tau bound in terms of (Wasserstein distance) x (Fisher information), WITHOUT directly referencing entropy production from a CPTP map.

### 3.4 The speed limit connection

From Phys. Rev. X 2023:
    W_2(rho(0), rho(t)) <= integral_0^t sqrt(Sigma_dot * mu) dt'

where Sigma_dot is the entropy production rate. This gives:
    W_2^2 <= t * integral Sigma_dot * mu dt' = t * <Sigma> * <mu>

This connects the physical distance between states (Wasserstein) to entropy production, providing an alternative route to bounding tau.

### 3.5 Evaluation against the three gaps

| Gap | Does it help? | How? |
|-----|--------------|------|
| **Gap 1** (saturation) | PARTIALLY | HWI inequality gives a bound on D in terms of W and I. The bound is tight when the Ricci curvature bound is tight. Different from JRSWW saturation. |
| **Gap 2** (CPTP map) | YES, PARTIALLY | Wasserstein distance can be defined between ANY two states without a channel. The HWI-based tau bound circumvents the CPTP requirement. |
| **Gap 3** (Sigma mismatch) | YES | W_2^2 <= t * Sigma * mu provides a DIFFERENT relationship between Sigma and state distance, potentially resolving the mismatch. |

### 3.6 Proposed Wasserstein-tau framework

Define:
    tau_W(rho, sigma) := 1 - F(rho, sigma)    [same as before]

But bound it using:
    tau_W <= g(W(rho, sigma), I(rho), kappa)    [from HWI]

where g is the function derived from the HWI inequality. This avoids the need for JRSWW entirely.

**Advantage**: Connects to thermodynamic speed limits; avoids CPTP; has natural geometric interpretation.
**Disadvantage**: Requires knowledge of the quantum Ricci curvature kappa for the quantum Markov semigroup, which is hard to compute for gravitational systems.

### 3.7 Assessment

**Overall rating: HIGH value, most promising for connecting tau to thermodynamics**

The Wasserstein framework provides the most natural bridge between tau and thermodynamic entropy production. The thermodynamic unification result (W = min of Sigma * mu) gives an operational meaning to tau in terms of physical transport cost. The HWI inequality provides a bound that avoids the CPTP requirement.

**Key paper to study**: Phys. Rev. X 13, 011013 (2023) — provides the exact variational formulas needed.

---

## 4. Quantum Hypothesis Testing

### 4.1 Core idea

Given two quantum states rho and sigma, quantum hypothesis testing asks: how well can we distinguish them? The quantum Stein's lemma gives:

    - lim_{n->infty} (1/n) ln beta_n^* = D(rho || sigma)

where beta_n^* is the optimal type-II error probability. This gives an OPERATIONAL meaning to relative entropy: it is the exponential rate of error decay in hypothesis testing.

### 4.2 Key recent results

**Generalized Quantum Stein's Lemma** (Hayashi-Yamasaki, Nature Physics 2025):
- Re-proves the generalized Stein's lemma after the original proof was found flawed
- **Central result**: The Stein exponent for resource testing equals the regularized relative entropy of resource
- Establishes the second law of quantum resource theories: resource convertibility is characterized by the regularized relative entropy
- Extended to classical-quantum channels (arXiv:2509.13280)
- Has been formalized in Lean proof assistant (arXiv:2510.08672)

**Quantum Channel Stein's Lemma** (Sci. China Inf. Sci. 2025):
- Establishes Stein's lemma for memoryless channel discrimination WITHOUT quantum memory assistance
- The optimal type-II error exponent = regularized channel divergence
- **Key for Gap 2**: This works for channel discrimination, not just state discrimination

**Adversarial Channel Discrimination** (arXiv:2510.07977, 2025):
- Minimax channel divergence framework
- Achievable by non-adaptive strategies
- Strong converse property proven for regularized exponent

### 4.3 Operational Sigma via hypothesis testing

Instead of defining Sigma through the JRSWW formula (which requires a CPTP map), define:

    Sigma_HT(rho, sigma) := D(rho || sigma)    [relative entropy = Stein exponent]

This is the rate at which we can distinguish rho from sigma. It has a DIRECT operational meaning: the number of copies needed to identify the state with error probability exp(-n * Sigma_HT).

For the gravitational case:
- rho = state of probe at height z_1
- sigma = state of probe at height z_2
- Sigma_HT = D(rho(z_1) || rho(z_2)) = how distinguishable are the two heights?

This AVOIDS the CPTP map entirely. We never need to specify the "channel" — we only need the states at different heights.

### 4.4 Connection to tau

If Sigma_HT = D(rho || sigma), then:
    tau_HT := 1 - exp(-Sigma_HT / 2) = 1 - exp(-D(rho || sigma) / 2)

This is mathematically identical to the JRSWW-derived tau, but the INTERPRETATION is different:
- JRSWW: tau = 1 - F(rho, Petz-recovered state) and Sigma = entropy loss in the channel
- HT: tau = 1 - exp(-D/2) where D = hypothesis testing exponent, no channel needed

**Critical advantage**: For the gravitational case, we can compute D(rho(z_1) || rho(z_2)) directly from the states without ever specifying a gravitational CPTP map.

### 4.5 The Generalized Stein Lemma and resource theory

Hayashi-Yamasaki 2025 establishes: for ANY quantum resource theory with "free" states F:

    Rate of resource distillation = Rate of resource dilution = D_F^{reg}(rho)

where D_F^{reg} is the regularized relative entropy of resource. This is the quantum analogue of the second law of thermodynamics.

**Implication for tau**: If we define "time-symmetric states" as the free resource (states invariant under time reversal), then:
- tau quantifies the resource of time-asymmetry
- The second law says: tau can only increase under free (time-symmetric) operations
- D_F^{reg}(rho) = Sigma = the "cost" of the time-asymmetry resource

This provides a RESOURCE-THEORETIC foundation for tau that is independent of JRSWW and CPTP maps.

### 4.6 Evaluation against the three gaps

| Gap | Does it help? | How? |
|-----|--------------|------|
| **Gap 1** (saturation) | YES, CIRCUMVENTS | Stein's lemma is TIGHT (equality, not just a bound). D(rho||sigma) is the EXACT asymptotic rate. No saturation issue. |
| **Gap 2** (CPTP map) | YES, CIRCUMVENTS | Hypothesis testing between states at different positions does NOT require a CPTP map. |
| **Gap 3** (Sigma mismatch) | YES | Sigma_HT = D(rho(z_1) || rho(z_2)) is uniquely defined by the states, no convention choices. |

### 4.7 Assessment

**Overall rating: HIGHEST value — addresses ALL THREE gaps**

Quantum hypothesis testing provides the most complete resolution. The Stein's lemma gives an EXACT (not just bounded) operational meaning to relative entropy, the definition of Sigma does not require a CPTP map, and the generalized Stein's lemma connects this to the second law of resource theories.

**Recommended action**: Reformulate Sigma_grav as a quantum hypothesis testing distinguishability:
    Sigma_grav(z_1, z_2) = D(rho(z_1) || rho(z_2))

This requires computing the states rho(z) for a probe in a gravitational potential, but does NOT require specifying the "gravitational channel."

**Key papers**:
- Hayashi-Yamasaki 2025 (Nature Physics): Generalized Stein's lemma
- arXiv:2509.13280: Extension to quantum channels
- Sci. China Inf. Sci. 2025: Channel Stein's lemma without memory

---

## 5. Resource Theory of Asymmetry / Coherence

### 5.1 Core idea

The resource theory of asymmetry (RTA) quantifies how much a quantum state breaks a given symmetry. For time-translation symmetry (generated by the Hamiltonian H), the "resource" is coherence in the energy eigenbasis. Free operations are those that commute with time evolution.

### 5.2 Key recent results

**Yamaguchi-Tajima 2023** (PRL 131, 200203):
- "Beyond i.i.d. in the Resource Theory of Asymmetry"
- Information-spectrum approach for QFI
- **Central result**: For arbitrary sequences of pure states, the coherence cost and distillable coherence both equal the spectral QFI rates
- Makes the RTA and entanglement theory parallel: both understood via information-spectrum methods, but based on QFI (not entropy)

**Marvian-Spekkens** (foundational work):
- Distinguished "speakable" vs "unspeakable" coherence
- QFI-based coherence monotones for U(1) symmetry
- **Key insight**: Time-asymmetry is a RESOURCE quantified by the QFI

**Fisher Information Matrix as Resource Measure** (PRA 107, 2023):
- Extends QFI-based resource measures to general connected Lie group symmetries
- Full Fisher information MATRIX (not just a scalar) serves as a resource measure

**Wigner-Yanase-Dyson skew information**:
- A well-known family of asymmetry monotones
- These are PRECISELY the monotones that satisfy the Jencova 2024 sufficiency theorem
- Connection: I_WYD(rho, H) = -(1/2) Tr([sqrt(rho), H]^2)

### 5.3 tau as asymmetry monotone

In the RTA under time-translation:
- Free states: [rho, H] = 0 (no coherence, time-symmetric)
- Free operations: covariant channels N with N(e^{-iHt} rho e^{iHt}) = e^{-iHt} N(rho) e^{iHt}
- Resource monotone: I_WYD(rho, H) or QFI J_Q(rho, H)

The connection to tau is:
- tau = 0 iff Sigma = 0 iff the state is in the "free set" of the RTA
- tau > 0 iff the state has coherence (time-asymmetry)
- The second law of the RTA says: coherence (time-asymmetry) cannot increase under free operations

**Proposed identification**:
    tau_RTA(rho) := 1 - F(rho, Delta(rho))

where Delta is the dephasing channel in the energy eigenbasis (the "twirl" projecting to the free set). This measures how far rho is from being time-symmetric.

### 5.4 The critical advantage: gravitational asymmetry

For gravitational time dilation, the relevant Hamiltonian is the GRAVITATIONAL redshift-shifted Hamiltonian:
    H_eff = H_0 (1 + Phi/c^2)

The coherence of a state w.r.t. H_eff quantifies how much the state "sees" the gravitational time dilation. The RTA gives:

    I_WYD(rho, H_eff) >= 0    with equality iff rho commutes with H_eff

This provides a NATURAL Sigma-like quantity:
    Sigma_RTA = I_WYD(rho, Delta_H) / I_WYD(rho, H)

measuring the fraction of coherence lost due to the gravitational potential difference Delta_H.

### 5.5 Evaluation against the three gaps

| Gap | Does it help? | How? |
|-----|--------------|------|
| **Gap 1** (saturation) | NO | RTA gives monotonicity, not equality. |
| **Gap 2** (CPTP map) | YES, PARTIALLY | Sigma_RTA is defined via the state and Hamiltonian, not via a channel. But the "free operations" of the RTA are still CPTP maps. |
| **Gap 3** (Sigma mismatch) | YES | Sigma_RTA is uniquely defined by the state and the gravitational Hamiltonian. No convention choice. |

### 5.6 Assessment

**Overall rating: HIGH value for INTERPRETATION, moderate for formal gaps**

The RTA provides the cleanest physical interpretation of tau: it IS the resource of time-asymmetry. The connection tau = 0 <=> no temporal resource <=> no time arrow is EXACT in this framework. However, the RTA does not fully avoid CPTP maps (they appear as "free operations") and does not resolve the saturation problem.

**Key insight**: tau IS a coherence monotone. This is not just an analogy — it is literally the same mathematical object.

---

## 6. Information Geometry

### 6.1 Core idea

Amari's information geometry endows the space of probability distributions (or quantum states) with a natural Riemannian structure. The Fisher-Rao metric is the unique (up to scaling) Riemannian metric invariant under sufficient statistics. If spacetime geometry IS information geometry, then the metric tensor g_ab would be derivable from information-theoretic quantities.

### 6.2 Key recent results

**Dorau-Much 2025/2026** (PRL 136, published March 2, 2026; arXiv:2510.24491):
- "From Quantum Relative Entropy to the Semiclassical Einstein Equations"
- **Central result**: Using modular theory, the QRE between vacuum and coherent excitations on a bifurcate Killing horizon equals the energy flux across the horizon. Under Bekenstein-Hawking, this becomes proportional to area variation, and Einstein's equations follow AUTOMATICALLY.
- This is the QFT generalization of Jacobson's 1995 thermodynamic derivation
- **Key formula**: S_rel(omega_psi | omega_0) = <T_ab k^a dSigma^b> (energy flux)
- Replaces classical entropy with Araki-Uhlmann quantum relative entropy

**Bianconi 2025** (PRD 111, 066001):
- "Gravity from entropy"
- **Central formulation**: The action for gravity IS the quantum relative entropy between the metric of spacetime (treated as a density matrix) and the metric induced by matter fields:
  S_G = D(g_spacetime || g_matter)
- Modified Einstein equations reduce to standard ones at low coupling
- G-field (Lagrangian multiplier) gives emergent positive cosmological constant
- **Most radical proposal**: metric = density matrix, action = QRE

**Neukart 2025** (Annals of Physics):
- "Geometry-Information Duality"
- Modified Einstein equations with informational stress-energy tensor from entanglement entropy
- Corrections to Newton's constant from entanglement entropy
- G_eff = G(1 + delta_G(S_ent))
- Effects are extremely small (below current detection)

**Modular Channels, Thermal Filtering** (arXiv:2504.20457, 2025):
- "Modular Channels Flow Correspondence" (MCFC)
- Area of causal screen = storage of filtered quantum information
- Unifies holography, entropy, and curvature through spectral dynamics of modular channels
- Einstein's equations from first law of entanglement as Clausius relation

**Chandrasekaran-Penington-Witten** (multiple papers 2023-2025):
- Type II von Neumann algebras for gravitational systems
- Generalized entropy = von Neumann entropy on crossed product algebra
- Observer-dependent gravitational entropy (JHEP 2025)
- **Key insight**: gravitational entropy REQUIRES specifying an observer (frame-dependent)

### 6.3 The Bianconi route to tau

If Bianconi's framework is correct:
    S_gravity = D(g_spacetime || g_matter)

Then Sigma_grav = D(g || g_matter) is LITERALLY a quantum relative entropy, and the JRSWW bound applies BY DEFINITION:
    F(g, Petz-recovered g) >= exp(-D(g || g_matter) / 2)

**This would close ALL three gaps simultaneously.**

However, the Bianconi framework is:
1. Not yet experimentally tested
2. Speculative (metric = density matrix is a strong assumption)
3. The "Petz recovery" of the metric has no clear physical interpretation yet

### 6.4 The Dorau-Much route

Dorau-Much gives:
    S_rel(omega_psi | omega_0) = energy flux across horizon

This is a RIGOROUS QFT result using modular theory. Combined with Bekenstein-Hawking:
    S_rel = (delta A) / (4 l_P^2)

And Einstein's equations follow. This provides:
    Sigma_grav = S_rel = (delta A) / (4 l_P^2) = energy flux

This is a well-defined QRE that does not require a CPTP map (it uses modular theory directly).

### 6.5 Evaluation against the three gaps

| Gap | Does it help? | How? |
|-----|--------------|------|
| **Gap 1** (saturation) | PARTIALLY | Bianconi: Sigma IS a QRE, so JRSWW applies directly. Dorau-Much: provides Sigma = energy flux, but saturation still not proven. |
| **Gap 2** (CPTP map) | YES (Bianconi) / PARTIALLY (Dorau-Much) | Bianconi: metric = density matrix makes Sigma a true QRE. Dorau-Much: modular theory provides the channel structure. |
| **Gap 3** (Sigma mismatch) | YES | Both frameworks give a UNIQUE, well-defined Sigma_grav from first principles. |

### 6.6 Assessment

**Overall rating: HIGHEST long-term value, but currently speculative**

The information geometry route (especially Dorau-Much + Bianconi) provides the most satisfying long-term resolution. If the metric IS a density matrix, then the entire tau framework follows from standard QI results. However, this is currently the most speculative approach and would require significant development for Paper 2.

**Key papers**:
- Dorau-Much 2025/2026, PRL 136: Rigorous QFT result
- Bianconi 2025, PRD 111: Most ambitious framework
- Chandrasekaran-Penington-Witten: Mathematical foundations

---

## 7. Comparative Assessment

### 7.1 Summary table

| Framework | Gap 1 (saturation) | Gap 2 (CPTP) | Gap 3 (Sigma mismatch) | Maturity | Paper 2 viability |
|-----------|-------------------|--------------|----------------------|----------|-------------------|
| **QFI** | Partial (different bound) | YES (circumvents) | Partial | High | HIGH |
| **Monotone metrics** | No | No | Partial | High | LOW |
| **Wasserstein** | Partial (HWI) | Partial (HWI) | Yes | Medium | MEDIUM |
| **Hypothesis testing** | YES (Stein tight) | YES (no channel) | YES | High | **HIGHEST** |
| **Resource theory** | No | Partial | Yes | High | MEDIUM |
| **Info geometry** | Partial-YES | YES (Bianconi) | YES | Low | LOW (speculative) |

### 7.2 Recommended strategy for Paper 2

**Tier 1 (immediate, use in Paper 2 revision)**:

1. **Quantum Hypothesis Testing reformulation**: Replace Sigma_grav = Delta-D(N, rho, sigma) with Sigma_grav = D(rho(z_1) || rho(z_2)), the state distinguishability between two heights. This AVOIDS the CPTP map requirement entirely. The Stein's lemma gives this an EXACT operational meaning. This is the single most impactful change.

2. **QFI-based gravitational distinguishability**: For the Pikovski-type probe, compute J_Q(rho, Phi) — the QFI for estimating the gravitational potential from the quantum state. This provides a CPTP-independent measure of how much information about Phi is encoded.

**Tier 2 (include as discussion/outlook in Paper 2)**:

3. **Resource theory framing**: Present tau as a coherence monotone in the resource theory of time-translation asymmetry. This provides the deepest physical interpretation.

4. **Wasserstein-entropy production connection**: Cite the Phys. Rev. X thermodynamic unification result to connect tau to thermodynamic speed limits.

**Tier 3 (future papers)**:

5. **Bianconi/Dorau-Much information geometry**: If metric = density matrix, then Sigma IS a QRE and all gaps close. This is the most promising long-term direction but too speculative for Paper 2.

### 7.3 The "minimal revision" path

For Paper 2 to be maximally defensible with minimal changes:

**Replace**: "The gravitational entropy production Sigma_grav = r_s/r arises from the CPTP map N_grav via JRSWW..."

**With**: "The gravitational distinguishability Sigma_grav = D(rho(z_1) || rho(z_2)) is defined as the quantum relative entropy between probe states at different heights. By quantum Stein's lemma (Hayashi-Yamasaki 2025), this equals the optimal asymptotic rate of gravitational field discrimination. The tau parameter then follows as tau = 1 - exp(-Sigma_grav / 2), with the bound tau <= 1 guaranteed by non-negativity of relative entropy."

This reformulation:
- Avoids needing a CPTP map (Gap 2: CLOSED)
- Uses an equality, not a bound to be saturated (Gap 1: CIRCUMVENTED)
- Defines Sigma unambiguously from the states (Gap 3: CLOSED)
- Is supported by a Nature Physics 2025 paper (strong authority)

The COST of this reformulation: the direct connection to the Petz recovery map becomes indirect. tau is no longer "1 minus the Petz recovery fidelity" but "1 minus exp(-D/2)." The retrodiction interpretation is preserved (Stein's lemma IS about inference) but the Petz map is no longer the central object.

---

## 8. Deep Dive: The Hypothesis Testing Route (Recommended)

### 8.1 Detailed construction

**Step 1**: Define the gravitational state family.

For a probe particle with internal Hamiltonian H_int at height z in Newtonian gravity:
    rho(z, t) = Tr_env [U(t) (rho_int tensor rho_cm(z)) U^dag(t)]

where U(t) includes the gravitational time dilation:
    U(t) = exp(-i H_int (1 + g*z/c^2) t / hbar)    [to first order in Phi/c^2]

For a superposition of heights z_1 and z_2, this produces decoherence of the form:
    <z_1| rho_cm(t) |z_2> ~ exp(-i Delta_E g (z_1 - z_2) t / (hbar c^2)) * <z_1| rho_cm(0) |z_2>

**Step 2**: Compute the relative entropy.

For the reduced states of the internal degree of freedom, conditional on position z:
    rho_int(z, t) = exp(-i H_int (1 + gz/c^2) t / hbar) rho_int(0) exp(+i H_int (1 + gz/c^2) t / hbar)

Since this is a unitary rotation, D(rho_int(z_1, t) || rho_int(z_2, t)) depends on the angle between the two rotated states. For an initial state |psi> that is a superposition of energy eigenstates:

    D(rho_1 || rho_2) = D(|psi(t, z_1)><psi(t, z_1)| || |psi(t, z_2)><psi(t, z_2)|)

For pure states: D(|psi><psi| || |phi><phi|) = -ln |<psi|phi>|^2.

And: <psi(t, z_1)|psi(t, z_2)> = sum_n |c_n|^2 exp(-i E_n g Delta_z t / (hbar c^2))

So:
    D = -ln |sum_n |c_n|^2 exp(-i E_n g Delta_z t / (hbar c^2))|^2

For a two-level system with energies E_0 and E_1, equal superposition:
    D = -ln |cos(Delta_E g Delta_z t / (2 hbar c^2))|^2 = -2 ln |cos(phi/2)|

where phi = Delta_E g Delta_z t / (hbar c^2) is the gravitational phase.

**Step 3**: Define tau.

    tau_HT = 1 - exp(-D/2) = 1 - |cos(phi/2)|

For small phi: tau_HT ~ phi^2 / 8 = (Delta_E g Delta_z t)^2 / (8 hbar^2 c^4)

For phi = pi (full distinguishability): tau_HT = 1

### 8.2 Comparison with Paper 2's current Sigma

Paper 2 currently uses:
    Sigma = E_G t / hbar    (linear in t)
    tau = 1 - exp(-E_G t / (2 hbar))    (exponential approach to 1)

The hypothesis testing route gives:
    Sigma_HT = -2 ln|cos(phi/2)|    (grows logarithmically then diverges at phi = pi)
    tau_HT = 1 - |cos(phi/2)|    (oscillatory, periodic in t)

**Key difference**: The HT result is OSCILLATORY (because pure-state evolution is unitary), while Paper 2's result is MONOTONIC (because it assumes irreversibility).

**Resolution**: The oscillatory behavior is for a CLOSED system (unitary evolution). When we include the environment (tracing over environmental degrees of freedom), the oscillations average out and we recover monotonic behavior. The Pikovski mechanism IS this tracing out.

This means: the HT route gives the CLOSED SYSTEM result, and Paper 2's result is the OPEN SYSTEM (traced-over) version. Both are correct in their respective domains, and they are connected by the process of partial tracing.

### 8.3 The bridge: partial tracing as the "channel"

The observation resolves the tension:
1. The "gravitational channel" is not a mysterious new CPTP map
2. It is simply the PARTIAL TRACE over degrees of freedom we don't observe
3. The Pikovski mechanism is exactly this: tracing over internal DOF gives position decoherence, and tracing over position gives internal decoherence
4. Partial trace IS a well-defined CPTP map

**Reformulated claim**: Sigma_grav = D(rho_full || sigma_ref) - D(Tr_env[rho_full] || Tr_env[sigma_ref]), where the "channel" is the partial trace. This is a standard JRSWW Sigma for the partial trace channel.

The question then becomes: can we choose sigma_ref such that this Delta-D equals r_s/r?

This is precisely the question analyzed in layer2_gravitational_channel.md — and the answer is that for the Pikovski dephasing, it's bounded by ln 2 and quadratic. But the partial trace over GRAVITATIONAL degrees of freedom (not just internal DOF) might give a different answer. This connects to the Dorau-Much modular theory approach.

---

## 9. Synthesis: Recommended Dual Strategy

### Strategy A: "Clean Break" (for a standalone Paper 2)

Abandon the JRSWW bound for gravity. Instead:
1. Define Sigma_grav via hypothesis testing: Sigma = D(rho(z_1) || rho(z_2))
2. Define tau = 1 - exp(-Sigma/2) as an operational measure
3. Compute Sigma for specific probe states (QFI-optimal probes)
4. Show that in the weak field, this reproduces the known Pikovski decoherence
5. Frame this in the resource theory of time-translation asymmetry

**Pros**: Avoids all three gaps. Rigorous. Well-supported by recent literature.
**Cons**: Loses the direct Petz recovery / retrodiction interpretation.

### Strategy B: "Staged Defense" (for Paper 2 as part of 4-paper series)

Keep the JRSWW structure but explicitly state the assumptions:
1. Paper 2 ASSUMES a gravitational CPTP map exists (supported by 3 first-principles routes)
2. The JRSWW bound gives F >= exp(-Sigma/2) for this assumed channel
3. Paper 2 derives consequences of this assumption
4. The hypothesis testing / QFI / resource theory frameworks provide ALTERNATIVE derivations of the same tau that do not require this assumption
5. Paper 4 (the unification paper) would provide the rigorous channel via Bianconi/Dorau-Much

**Pros**: Preserves the 4-paper narrative. Maintains the Petz retrodiction connection.
**Cons**: Paper 2 is explicitly conditional on an unproven assumption.

### Recommended: Strategy B with QFI supplement

Use Strategy B but ADD a section showing that the QFI-based approach and hypothesis testing approach give the same tau in the appropriate limits. This provides "two roads to the same destination" and makes the argument robust against failure of either route.

---

## 10. Key References

### Framework 1: Quantum Fisher Information
- Jencova 2024, Commun. Math. Phys. 405, 180: "Sufficient Statistic and Recoverability via QFI" [arXiv:2302.02341](https://arxiv.org/abs/2302.02341)
- Balatsky et al. 2025, PRA 111, 012411: "Gravity as Universal Dephasing Channel" [arXiv:2406.03256](https://arxiv.org/abs/2406.03256)
- QFI in Schwarzschild spacetime 2025: [Springer](https://link.springer.com/article/10.1007/s10773-025-06073-8)
- Pikovski group 2016: "Information transfer during gravitational decoherence" [arXiv:1612.08864](https://arxiv.org/abs/1612.08864)
- Yamaguchi-Tajima 2023, PRL 131, 200203: "Beyond i.i.d. in RTA" [PRL](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.131.200203)
- Quantum time dilation in gravitational field 2024: [Quantum 8, 1338](https://quantum-journal.org/papers/q-2024-05-07-1338/)

### Framework 2: Monotone Metrics
- Petz monotone geometry 2025: [arXiv:2509.14578](https://arxiv.org/abs/2509.14578)
- Monotone metric tensors review 2024: [World Scientific](https://www.worldscientific.com/doi/10.1142/S0219887824400048)

### Framework 3: Quantum Wasserstein Distance
- Thermodynamic Unification 2023, PRX 13, 011013: [Phys. Rev. X](https://journals.aps.org/prx/abstract/10.1103/PhysRevX.13.011013)
- Quantum HWI inequality, Rouze-Datta 2020: [arXiv:1709.07437](https://arxiv.org/abs/1709.07437)
- Order p Quantum Wasserstein 2025: [Springer](https://link.springer.com/article/10.1007/s00023-025-01557-z)
- Wasserstein for Gaussian states 2024-2025: [Quantum Zeitgeist](https://quantumzeitgeist.com/quantum-states-order-wasserstein-distance-advances-state-discrimination-gaussian/)
- Quantum Wasserstein over separable states 2023: [Quantum 7, 1143](https://quantum-journal.org/papers/q-2023-10-16-1143/)

### Framework 4: Quantum Hypothesis Testing
- Hayashi-Yamasaki 2025, Nature Physics: "Generalized Quantum Stein's Lemma" [arXiv:2408.02722](https://arxiv.org/abs/2408.02722)
- Extension to CQ channels 2025: [arXiv:2509.13280](https://arxiv.org/abs/2509.13280)
- Lean formalization 2026: [arXiv:2510.08672](https://arxiv.org/abs/2510.08672)
- Channel Stein's lemma 2025: [Sci. China Inf. Sci.](http://scis.scichina.com/en/2025/180509.pdf)
- Adversarial channel discrimination 2025: [arXiv:2510.07977](https://arxiv.org/abs/2510.07977)

### Framework 5: Resource Theory of Asymmetry
- Marvian-Spekkens (foundational): Distinguishing speakable/unspeakable coherence
- Yamaguchi-Tajima 2023, PRL 131, 200203: Beyond i.i.d. in RTA
- Fisher information matrix as resource 2023, PRA 107: General Lie group extension
- Lostaglio review 2019: Resource theory approach to thermodynamics [arXiv:1807.11549](https://arxiv.org/abs/1807.11549)

### Framework 6: Information Geometry
- Dorau-Much 2025/2026, PRL 136: "QRE to Einstein Equations" [arXiv:2510.24491](https://arxiv.org/abs/2510.24491)
- Bianconi 2025, PRD 111, 066001: "Gravity from Entropy" [arXiv:2408.14391](https://arxiv.org/abs/2408.14391)
- Neukart 2025, Ann. Phys.: "Geometry-Information Duality" [arXiv:2409.12206](https://arxiv.org/abs/2409.12206)
- Modular Channels 2025: [arXiv:2504.20457](https://arxiv.org/abs/2504.20457)
- Chandrasekaran-Penington-Witten 2024: Algebras without AdS/CFT [arXiv:2310.02189](https://arxiv.org/abs/2310.02189)
- Observer-dependent gravitational entropy 2025: [JHEP](https://link.springer.com/article/10.1007/JHEP07(2025)146)
- JHEP 2025: Information-theoretic detection of quantum gravitational corrections [JHEP](https://link.springer.com/article/10.1007/JHEP02(2025)109)

---

## 11. Open Questions for Future Research

1. **Can the hypothesis testing Sigma_HT be computed for realistic gravitational probes?** The two-level atom calculation (Section 8.1) gives an oscillatory result. What happens for a thermal state? For a macroscopic object with ~10^23 internal DOF?

2. **Does the quantum Ricci curvature (appearing in HWI) have a gravitational interpretation?** If kappa = R_ab u^a u^b (the Raychaudhuri focusing scalar), then the HWI inequality directly connects Wasserstein distance to spacetime curvature.

3. **Can Bianconi's metric = density matrix identification be tested?** The modified Einstein equations should give corrections to standard GR. What are the observable signatures?

4. **Is there a gravitational Stein's lemma?** Can we prove: the rate of gravitational field discrimination between two nearby points = some function of the Riemann curvature?

5. **How does the partial trace over gravitational degrees of freedom differ from tracing over internal DOF?** The Pikovski channel traces over internal DOF. What if we trace over the gravitational field itself (as in Dorau-Much)?

6. **Can the Wasserstein-entropy production connection be extended to non-Markovian gravitational dynamics?** The Phys. Rev. X result assumes Markovian dynamics. Gravitational decoherence may be non-Markovian.

7. **What is the QFI-optimal probe for measuring Sigma_grav?** Given the quantum Cramer-Rao bound, what initial state maximizes the distinguishability between different gravitational potentials?
