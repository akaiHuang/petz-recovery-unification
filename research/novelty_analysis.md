# Rigorous Novelty Analysis of the tau Framework

**Date**: 2026-03-11
**Purpose**: Brutally honest classification of what is genuinely new in Sheng-Kai Huang's tau framework versus what is known, repackaged, or synthesized from others' work.

## Classification Key

| Code | Category | Description |
|------|----------|-------------|
| **(A)** | Known result (just cited) | Published by others; Huang adds no new content |
| **(B)** | Known result, new notation | Same mathematics with different symbols/names |
| **(C)** | Genuinely new synthesis | Connecting previously unconnected results into a unified picture |
| **(D)** | Genuinely new interpretation | New physical meaning assigned to existing mathematics |
| **(E)** | Genuinely new mathematical result | New theorem, bound, or proof |
| **(F)** | Genuinely new physical prediction | Testable consequence not previously stated |

---

## Executive Summary

### What is genuinely new:
1. The **specific definition** tau = 1 - F(Petz) as a named "temporal asymmetry parameter" and the **interpretive framework** that this single number unifies the quantum eraser, QEC, thermodynamic irreversibility, and gravity **(C/D)**
2. The **explicit identification** of the quantum eraser as a Petz recovery problem with exact computation **(C)** -- though this is close to obvious once stated
3. The **triple identification** c_eff/c = exp(-Sigma_grav/2) = F_bound connecting coordinate speed of light to Petz recovery fidelity **(C/D)** -- genuinely novel as a stated equivalence
4. The **saturation theorem** (F^2 >= exp(-Delta D) saturated iff quantum sufficient statistic) **(E)** -- if proven correctly, this appears to be a new mathematical result
5. The **triangle inequality** for sqrt(tau) under channel composition **(E)** -- if proven correctly, this appears to be a new mathematical result
6. The **gravitational-wave echo prediction** at Delta t ~ 4 r_s/c from the exponential metric **(F)** -- though the exponential metric itself is known (Yilmaz/Papapetrou)
7. The claim that **tau = 1 at the apparent horizon** (not event horizon) as the operationally relevant boundary **(C/D)**

### What is NOT new (but may appear so due to presentation):
1. The Petz recovery map itself (Petz 1986, 1988)
2. The universal recovery bound F >= exp(-Sigma/2) (JRSWW 2018, Fawzi-Renner 2015)
3. Petz = retrodiction functor (Parzygnat-Buscemi 2023)
4. Petz = reverse channel in quantum Crooks (Kwon-Kim 2019)
5. Petz in entanglement wedge reconstruction (Chen-Penington-Salton 2020)
6. The exponential metric g_00 = -exp(-r_s/r) (Yilmaz 1958, Papapetrou, Rosen)
7. QRE -> Einstein equations (Jacobson 1995, Dorau-Much 2025)
8. HP decoding = Petz map for scrambling channels (Yoshida-Kitaev 2017, Nakayama et al. 2023)
9. Traversable wormholes as entanglement-assisted quantum channels (Bao-Chatwin-Davies 2018)
10. The gravitational Landauer principle (Herrera 2020, Daffertshoffer-Plastino)
11. Pikovski decoherence from gravitational time dilation (Pikovski et al. 2015)
12. Gravitational redshift as a quantum channel (Ahmadi-Fuentes 2014)
13. NEC violation and wormhole traversability (standard GR, Gao-Jafferis-Wall 2017)

---

## Detailed Classification Table

### Paper 1 Claims

| # | Claim | Classification | Priority Check | Assessment |
|---|-------|---------------|----------------|------------|
| 1 | **tau = 1 - F as measure of time's arrow** | **(B/D)** | Infidelity 1-F is a standard measure. Kwon-Kim (2019) already connected Petz recovery to entropy production and irreversibility. Bai-Buscemi-Scarani (2024, arXiv:2412.12489) constructed entropy production operators using Petz retrodiction. **However**, the specific naming "temporal asymmetry parameter" and the interpretive claim that tau = 0 iff no time arrow is Huang's packaging. | **Moderate novelty**. The definition is trivial (1 minus a known quantity), but the interpretive framework around it is new. The risk: a referee may say "this is just infidelity with a name." |
| 2 | **Equivalence chain: Petz <=> entropy production <=> Crooks <=> QEC** | **(C)** | Each pairwise connection existed: Petz-QEC (Barnum-Knill 2002), Petz-Crooks (Kwon-Kim 2019), Petz-retrodiction (Parzygnat-Buscemi 2023), entropy-production-irreversibility (standard thermodynamics). **But no single paper assembled all four into a unified chain.** | **This is the strongest contribution of Paper 1.** The synthesis is genuine -- connecting four communities (quantum foundations, QEC, quantum thermo, retrodiction) through one map. Closest prior: Bai-Buscemi-Scarani (2024) covers Crooks + Petz + entropy production but not QEC or eraser. |
| 3 | **tau <= 1 - exp(-Sigma/2) bound** | **(A/B)** | This is a direct restatement of the JRSWW bound (Junge-Renner-Sutter-Wilde-Winter 2018) and Fawzi-Renner (2015) in terms of tau = 1-F. No new mathematics whatsoever. The paper itself labels this "[restating the universal recovery bound of Junge et al.]". | **Not novel.** Honest labeling in the paper is good. |
| 4 | **Retrodiction Landauer Principle** | **(C/D)** | Classical Landauer principle: Landauer 1961. Quantum extensions: Reeb-Wolf 2014. Gravitational Landauer: Herrera 2020. Retrodiction + Petz: Parzygnat-Buscemi 2023. **The specific statement "the cost of retrodiction is bounded by entropy production via Landauer" appears to be a new synthesis.** | **Moderate novelty**. The name is new; the underlying mathematics combines known bounds. |
| 5 | **Quantum eraser = Petz recovery (explicit identification)** | **(C)** | No prior paper explicitly wrote Petz_{sigma,Tr_E}(rho_S) = |psi><psi|_{SE} for the quantum eraser setup with a worked computation. The eraser has been discussed in terms of entanglement and decoherence by Scully, Kim, etc., but not through the Petz map formalism. | **Genuine but incremental novelty.** Once you know the Petz map and the eraser, the identification is almost immediate (as the paper itself notes for pure bipartite states, DPI is always saturated). The contribution is making the obvious explicit. |
| 6 | **F^2 >= exp(-Delta D) saturated iff quantum sufficient statistic** | **(E)** | The DPI equality condition (saturation iff Petz exact recovery) is Petz's theorem (1988). The JRSWW bound strengthens this to approximate recovery. **The specific characterization of when the JRSWW bound is tight as "iff quantum sufficient statistic" may be new.** However, this is closely related to Petz's original sufficiency theorem. Need to check Sutter et al. (2016) and the "optimality condition for the Petz map" paper (arXiv:2410.23622, 2024). | **Potentially novel if the proof is correct**, but may overlap with known sufficiency conditions. Significance: moderate (characterizing tightness of known bounds). |
| 7 | **sqrt(tau) satisfies triangle inequality under channel composition** | **(E)** | The fact that sqrt(1-F) is related to metrics (Bures distance, purified distance) is well-known. However, the specific statement about channel composition (not just state comparison) and the claim that this property is **unique to the standard Petz map** (failing for rotated Petz) appears to be new. | **Potentially novel** if the proof is correct and the uniqueness claim holds. The Fuchs-van de Graaf inequalities relate sqrt(1-F) to trace distance, which is a metric. But the compositional statement for channels specifically is a different (and potentially new) claim. |

### Paper 2 Claims

| # | Claim | Classification | Priority Check | Assessment |
|---|-------|---------------|----------------|------------|
| 8 | **Sigma_grav = -ln(-g_00) = r_s/r (gravitational entropy production)** | **(B/D)** | The quantity -ln(-g_00) appears implicitly in multiple contexts: Tolman temperature T_loc = T_inf/sqrt(-g_00), gravitational redshift z+1 = 1/sqrt(-g_00), the Unruh-DeWitt detector response. The Dorau-Much (2025 PRL) framework derives Einstein equations from QRE, implicitly using the same mathematical structure. Herrera (2020) obtains the same via Landauer + Tolman. Ahmadi-Fuentes (2014) model redshift as a quantum channel with transmissivity g_00. **The specific naming "gravitational entropy production" and the definition Sigma_grav = r_s/r is Huang's, but the mathematical content is distributed across existing literature.** | **New notation/interpretation for known mathematics.** The value-add is giving this quantity a unified name and connecting it to the tau framework. A referee familiar with Jacobson (1995), Verlinde (2011), or Dorau-Much (2025) may find this incremental. |
| 9 | **c_eff/c = exp(-Sigma_grav/2) = F_bound (triple identification)** | **(C/D)** | c_eff/c = exp(Phi/c^2) is known (Einstein 1911, Dicke 1957). F >= exp(-Sigma/2) is JRSWW (2018). **The identification that these two expressions are the same quantity is new.** No prior paper states "the coordinate speed of light IS the Petz recovery fidelity bound." | **This is Paper 2's strongest claim.** The synthesis is genuine -- it connects a GR result (1911) to a QI result (2018) through the intermediary of Sigma_grav. The question is whether it's profound or tautological. The answer depends on whether Sigma_grav is well-motivated as "entropy production" or is just a definition chosen to make the identification work. The three independent routes (modular flow, Landauer, quantum channel) strengthen the case that Sigma_grav = r_s/r is natural, not ad hoc. |
| 10 | **Exponential metric g_00 = -exp(-r_s/r) from tau framework** | **(B/D)** | The exponential metric is due to Yilmaz (1958), Papapetrou, and Rosen. It is a known alternative to Schwarzschild that differs at strong-field. Recent work (arXiv:1805.03781, 2018; arXiv:2510.15391, 2025) establishes it as a traversable wormhole supported by a phantom scalar field. **Huang's contribution is interpreting this metric as "the one that saturates the Petz bound."** | **New interpretation of a known metric.** The exponential metric is well-studied. The claim that it "saturates the Petz bound" is Huang's interpretation. Whether this interpretation is physically meaningful depends on whether the gravitational CPTP map is well-defined (which remains an open problem, as the paper honestly acknowledges). |
| 11 | **Event horizon = tau = 1 (complete information loss)** | **(C/D)** | That event horizons represent maximal information loss is a qualitative statement going back to Hawking (1975). The specific quantification tau = 1 at the apparent (not event) horizon, with tau_K tracking the AH in Vaidya spacetime, is Huang's framework. The paper's demonstration that tau_K = 1 at AH but tau_K < 1 at EH during accretion is a concrete new calculation. | **Moderately novel synthesis.** The insight that tau should track the apparent horizon is physically well-motivated and consistent with the quantum extremal surface prescription (Engelhardt-Wall 2015). The Vaidya calculation with specific numerical values (Table II in the paper) is new. |
| 12 | **NEC violation <=> tau < 0 (via Raychaudhuri equation)** | **(C/D)** | The Raychaudhuri equation is standard GR. The connection between NEC violation and exotic phenomena (wormholes, etc.) is standard. The Raychaudhuri-tau decomposition dSigma/ds = theta^2/3 + sigma^2 - omega^2 + R_ab u^a u^b is Huang's rewriting using his notation. **The identification that NEC violation maps to tau < 0 in his framework is a new interpretive claim.** | **Moderate novelty.** This is primarily a notational translation of known GR into the tau language. The physical content (NEC violation enables information recovery) is known. |
| 13 | **Pikovski decoherence as Petz recovery problem** | **(D)** | Pikovski et al. (2015) showed gravitational time dilation causes decoherence. Reframing this as a Petz recovery problem (asking "how well can the Petz map undo Pikovski decoherence?") is Huang's interpretation. | **New interpretation, but no new calculation.** It would be significantly strengthened if the paper computed the actual Petz recovery fidelity for the Pikovski channel and showed it equals something meaningful. |

### Wormhole/QEC Claims (New in this session)

| # | Claim | Classification | Priority Check | Assessment |
|---|-------|---------------|----------------|------------|
| 14 | **Wormhole QEC: scrambling + entanglement = automatic error correction with tau_eff < 0** | **(C/D)** | That scrambling provides error correction is known: Hayden-Preskill (2007), Yoshida-Kitaev (2017). That wormholes are entanglement-assisted quantum channels is in Bao-Chatwin-Davies (2018). That HP decoding = Petz recovery is shown by Nakayama et al. (2023). **The specific combination "tau_eff < 0 means entanglement overcomes channel noise" is Huang's interpretive framework.** | **Moderate novelty as interpretation.** The underlying physics (entanglement-assisted QEC, wormhole traversability) is well-established. The tau_eff framing is new packaging. |
| 15 | **tau_eff = Sigma_channel - I(L;R) as figure of merit for entanglement-assisted QEC** | **(C)** | In entanglement-assisted QEC (Brun-Devetak-Hsieh 2006), the capacity enhancement from shared entanglement is well-characterized. The hashing bound Q_EA = (1/2)I(A;B) is known. **The specific formula tau_eff = Sigma - I(L;R) does not appear in the literature in this form.** However, the idea that mutual information offsets channel entropy is the entire point of entanglement-assisted capacity theorems. | **New notation for a known concept.** The formula repackages the EA capacity idea. The novelty is in naming this tau_eff and connecting it to the broader tau framework. A quantum information theorist would recognize this as a restatement of coherent information / EA capacity bounds. |
| 16 | **Lab wormholes (GJW protocol) as practical tau < 0 generators** | **(D)** | The GJW (Gao-Jafferis-Wall 2017) protocol and its lab realization (Jafferis et al. 2022 on Google Sycamore) are known. The interpretation as "tau < 0 generators" is Huang's framing. | **Pure interpretation.** No new physics or predictions. |
| 17 | **Connection between wormhole traversability and quantum chip optimization** | **(D)** | That holographic QEC insights can inform practical QEC is a general theme in the field (e.g., AlphaQubit decoder). **A specific, actionable connection between wormhole traversability and quantum chip design has not been established** by anyone, including Huang. | **Speculative interpretation.** No concrete results or predictions. |

---

## Priority Search: Has Anyone Published Similar Work?

### Closest competitors to the overall tau framework:

1. **Bai, Buscemi, and Scarani (arXiv:2412.12489, Dec 2024)**: "Fully quantum stochastic entropy production." They construct an entropy production OPERATOR using Petz retrodiction. This is the closest independent work -- it connects Crooks, Petz, Jarzynski, and the second law in one framework. **Key difference from Huang**: they do not define tau = 1-F, do not discuss QEC or the quantum eraser, and do not extend to gravity. They go deeper mathematically (operator-level) while Huang goes broader (more connections, less depth).

2. **Kwon and Kim (Phys. Rev. X, 2019)**: "Fluctuation Theorems for a Quantum Channel." They adopt the Petz map as the reverse channel and extend entropy production to quantum channels. This is a direct precursor to Huang's connection between Petz and thermodynamic irreversibility. **Key difference**: Kwon-Kim do not name tau = 1-F, do not discuss the quantum eraser, and do not connect to QEC.

3. **Parzygnat and Buscemi (Quantum, 2023)**: "Axioms for retrodiction." They prove the uniqueness of the Petz map as a retrodiction functor. This is the mathematical foundation Huang builds on. **Key difference**: Parzygnat-Buscemi prove the theorem; Huang interprets it physically.

4. **Nakayama et al. (PTEP, 2023)**: "The Petz (lite) recovery map for scrambling channel." They show HP decoding = Petz in specific models. **Key difference**: They work within black hole physics; Huang connects this to the broader equivalence chain.

5. **Dorau and Much (PRL, 2025)**: "From Quantum Relative Entropy to the Semiclassical Einstein Equations." They derive Einstein equations from QRE, building on Jacobson (1995). This is the closest prior work to Huang's Paper 2 claim about Sigma_grav. **Key difference**: Dorau-Much work at the level of Einstein equations; Huang works at the level of the Petz bound. Dorau-Much do not define Sigma_grav = r_s/r or connect to the speed of light.

6. **Bianconi (PRD 111, 2025)**: "Gravity from entropy." Proposes action = QRE between spacetime metric and matter-induced metric. **Key difference**: Bianconi derives modified Einstein equations; Huang works at the operational level of information recovery. The philosophy is similar (gravity = entropy) but the technical approach is entirely different.

### Verdict on priority:

**No single paper in the literature contains the tau framework as presented by Huang.** The closest is Bai-Buscemi-Scarani (2024), which covers Crooks+Petz+entropy production but not QEC, not the eraser, and not gravity. For Paper 2, no one has published the triple identification c_eff/c = exp(-Sigma_grav/2) = F_bound.

---

## Honest Assessment of Genuine Contribution

### What is truly original:

1. **The four-way equivalence chain (Paper 1)** -- This is the most defensible contribution. Each link was known; the chain was not assembled. Significance: **HIGH** for pedagogy and conceptual clarity, **MODERATE** for advancing the field technically.

2. **The triple identification c_eff/c = F_bound (Paper 2)** -- This is striking and, as far as I can find, genuinely new. Its significance depends critically on the open "channel problem" (no canonical gravitational CPTP map exists). Without solving this, the identification is suggestive rather than proven. Significance: **HIGH if the channel problem is solved**, **MODERATE as a suggestive framework**.

3. **The two mathematical results (saturation theorem, triangle inequality) in Paper 1** -- If correctly proven, these are genuine new mathematical contributions to quantum information theory. Significance: **MODERATE** (useful technical results, not breakthroughs).

4. **tau tracking the apparent (not event) horizon** -- This is a novel physical observation with concrete calculations (Vaidya spacetime). Significance: **MODERATE-HIGH** (connects nicely to modern quantum gravity via QES prescription).

### What is NOT original but is presented as if it might be:

1. **Sigma_grav = r_s/r** -- This is essentially -ln(-g_00) in different notation. The three "derivations" are really three previously-known results that all give the same answer. The paper honestly labels them as [KNOWN], which is good.

2. **The exponential metric** -- Known since Yilmaz (1958). The "interpretation" as Petz-bound-saturating is Huang's, but the metric and its properties (no event horizon, traversable wormhole) were known.

3. **tau_eff = Sigma - I(L;R) for wormholes** -- This repackages entanglement-assisted capacity theory in new notation.

### Fundamental vulnerability:

The entire Paper 2 framework rests on an **unproven assumption**: that there exists a canonical gravitational CPTP map N_grav such that the QRE decrease equals r_s/r. The paper honestly calls this "an open problem" and provides three circumstantial arguments. But until this channel is constructed from first principles, the triple identification remains a **conjecture**, not a theorem. The paper's intellectual honesty about this gap is commendable, but a skeptical referee will press hard on this point.

---

## Suggestions for Strengthening Originality

### For Paper 1:
1. **Emphasize the equivalence chain more**, the tau definition less. The chain is the real contribution; tau = 1-F is just notation.
2. **The mathematical results (Theorems 1 and 2) should be prominently featured** -- these are the hardest-to-dismiss contributions because they are provably new (if correct).
3. **Be explicit about what is known vs new** (the paper already does this well).
4. **Consider adding a concrete prediction or protocol** that the equivalence chain enables but that was not previously accessible.

### For Paper 2:
1. **The channel problem is the Achilles' heel.** If possible, construct N_grav explicitly (even for a toy model). The Trejo-Calderon modular channel approach mentioned in the paper seems promising.
2. **The echo prediction is the strongest selling point** -- quantify the expected signal strength and compare to existing LIGO constraints.
3. **The AH vs EH distinction and the Vaidya calculation are strong** -- expand these.
4. **Distance yourself from strong claims** about "the exponential metric IS the correct metric." Instead, present it as "the metric that saturates our bound" and note it is testable.

### For the wormhole/QEC claims:
1. **These are currently too speculative** for publication. They would need: (a) an explicit calculation of tau_eff for the GJW protocol, (b) comparison with the actual teleportation fidelity in the Jafferis et al. (2022) experiment, (c) a concrete protocol for quantum chip optimization.
2. **The connection to Bao-Chatwin-Davies (2018)** should be made explicit -- they already studied wormholes as quantum channels. What does the tau framing add that they did not have?

---

## Summary Scorecard

| Claim | Category | Novelty | Significance |
|-------|----------|---------|-------------|
| tau = 1-F definition | (B/D) | Low-Moderate | LOW (notation) |
| Four-way equivalence chain | (C) | **HIGH** | **HIGH** |
| tau <= 1-exp(-Sigma/2) | (A) | None | N/A (known bound) |
| Retrodiction Landauer Principle | (C/D) | Moderate | Moderate |
| Eraser = Petz | (C) | Moderate | Moderate |
| Saturation theorem | (E) | **HIGH** | Moderate |
| Triangle inequality for sqrt(tau) | (E) | **HIGH** | Moderate |
| Sigma_grav = r_s/r | (B/D) | Low | Low-Moderate |
| c_eff/c = F_bound | (C/D) | **HIGH** | **HIGH** (if channel problem solved) |
| Exponential metric from tau | (B/D) | Low | Moderate |
| tau = 1 at apparent horizon | (C/D) | **Moderate-High** | **Moderate-High** |
| NEC violation = tau < 0 | (C/D) | Moderate | Moderate |
| Pikovski as Petz problem | (D) | Low | Low (no calculation) |
| Wormhole QEC with tau_eff | (C/D) | Moderate | Low (speculative) |
| tau_eff = Sigma - I(L;R) | (B/C) | Low | Low (repackaged EA capacity) |
| Lab wormholes as tau<0 | (D) | Low | Low (pure interpretation) |
| GW echo prediction | (F) | Moderate | **HIGH** (testable) |

---

## Final Verdict

The tau framework's genuine contribution is **synthesis, not discovery**. It does not introduce new physics or new mathematics of the first rank, but it provides a **unified language** that connects quantum information recovery, thermodynamic irreversibility, quantum error correction, and (tentatively) gravitational physics through a single quantity.

This is a legitimate and potentially valuable contribution -- many important papers in physics are "Rosetta Stones" that translate between communities. The paper's own framing acknowledges this: "We have not discovered new things, but organized the nature of known things."

**The key risk**: a referee may dismiss this as "just notation." To defend against this, the paper should:
1. Lead with the mathematical theorems (which are provably new)
2. Show concrete utility of the equivalence chain (what can you calculate now that you couldn't before?)
3. For Paper 2, emphasize testable predictions over interpretive claims

**Bottom line**: The framework contains approximately 10-15% genuinely new content (the equivalence chain, two theorems, the triple identification, the echo prediction), 25% new synthesis, 25% new interpretation, and 35-40% known results in new notation. This is comparable to many successful conceptual physics papers, but the authors should be vigilant about not overclaiming.

---

## References

- Petz, D. (1986, 1988). Sufficient subalgebras and the relative entropy of states of a von Neumann algebra.
- Fawzi, O. and Renner, R. (2015). Quantum conditional mutual information and approximate Markov chains. Commun. Math. Phys. 340, 575-611.
- Junge, M., Renner, R., Sutter, D., Wilde, M.M., Winter, A. (2018). Universal recovery maps and approximate sufficiency of quantum relative entropy. Ann. Henri Poincare 19, 2955-2978.
- Parzygnat, A.J. and Buscemi, F. (2023). Axioms for retrodiction: achieving time-reversal symmetry with a prior. Quantum 7, 1013.
- Kwon, H. and Kim, M.S. (2019). Fluctuation Theorems for a Quantum Channel. Phys. Rev. X 9, 031029.
- Bai, G., Buscemi, F., and Scarani, V. (2024). Fully quantum stochastic entropy production. arXiv:2412.12489.
- Barnum, H. and Knill, E. (2002). Reversing quantum dynamics with near-optimal quantum and classical fidelity. J. Math. Phys. 43, 2097-2106.
- Chen, C.F., Penington, G., and Salton, G. (2020). Entanglement wedge reconstruction using the Petz map. JHEP 01, 168.
- Yoshida, B. and Kitaev, A. (2017). Efficient decoding for the Hayden-Preskill protocol. arXiv:1710.03363.
- Nakayama, Y. et al. (2023). The Petz (lite) recovery map for scrambling channel. PTEP 2023, 123B04.
- Jacobson, T. (1995). Thermodynamics of Spacetime: The Einstein Equation of State. Phys. Rev. Lett. 75, 1260.
- Dorau, J.B. and Much, A. (2025). From Quantum Relative Entropy to the Semiclassical Einstein Equations. Phys. Rev. Lett.
- Bianconi, G. (2025). Gravity from entropy. Phys. Rev. D 111, 066001.
- Herrera, L. (2020). Landauer Principle and General Relativity. Entropy 22, 407.
- Ahmadi, M. and Fuentes, I. (2014). Relativistic quantum metrology / gravitational channel models.
- Pikovski, I., Zych, M., Costa, F., and Brukner, C. (2015). Universal decoherence due to gravitational time dilation. Nature Physics 11, 668-672.
- Yilmaz, H. (1958). New approach to general relativity. Phys. Rev. 111, 1417.
- Einstein, A. (1911). On the Influence of Gravitation on the Propagation of Light.
- Dicke, R.H. (1957). Gravitation without a principle of equivalence. Rev. Mod. Phys. 29, 363.
- Gao, P., Jafferis, D.L., and Wall, A.C. (2017). Traversable Wormholes via a Double Trace Deformation. JHEP 12, 151.
- Jafferis, D. et al. (2022). Traversable wormhole dynamics on a quantum processor. Nature 612, 51-55.
- Bao, N. and Chatwin-Davies, A. (2018). Traversable Wormholes as Quantum Channels. JHEP 11, 071.
- Pastawski, F., Yoshida, B., Harlow, D., and Preskill, J. (2015). Holographic quantum error-correcting codes. JHEP 06, 149.
- Almheiri, A., Dong, X., and Harlow, D. (2015). Bulk Locality and Quantum Error Correction in AdS/CFT. JHEP 04, 163.
- Verlinde, E.P. (2016). Emergent Gravity and the Dark Universe. SciPost Phys. 2, 016.
- Engelhardt, N. and Wall, A.C. (2015). Quantum Extremal Surfaces. JHEP 01, 073.
- Brun, T.A. (2003). Computers with closed timelike curves can solve hard problems. Found. Phys. Lett. 16, 245-253.
- Kolchinsky, A. and Wolpert, D.H. (2017). Dependence of dissipation on the initial distribution over states. J. Stat. Mech. 083202.
- Reeb, D. and Wolf, M.M. (2014). An improved Landauer principle with finite-size corrections. New J. Phys. 16, 103011.
