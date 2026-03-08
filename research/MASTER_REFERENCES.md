# Master Reference Database: The 4+1 Paper Series

**Author**: Sheng-Kai Huang
**Last Updated**: 2026-03-07
**Scope**: All references cited or needed across Papers 1, 1b, 2, 3, and 4

---

## Legend

**Paper tags**: [P1] = Paper 1, [P1b] = Paper 1b, [P2] = Paper 2, [P3] = Paper 3, [P4] = Paper 4

**Status tags**:
- [FOUNDATIONAL] = Established result our work builds upon
- [SUPPORTS] = Evidence supporting our claims
- [CONSTRAINS] = Places limits on or refines our claims
- [INDEPENDENT VERIFICATION] = Independently arrived at similar structures
- [POTENTIAL CONFLICT] = May challenge some aspect of our framework
- [OUR WORK] = Self-citations

---

## A. Petz Recovery & Quantum Information

### A1. Petz 1986
- **Authors**: D. Petz
- **Title**: "Sufficient subalgebras and the relative entropy of states of a von Neumann algebra"
- **Journal**: Commun. Math. Phys. **105**, 123 (1986)
- **arXiv**: --
- **Cited in**: [P1] [P1b] [P2]
- **Key result**: Introduced the Petz recovery map; proved it is CPTP.
- **Status**: [FOUNDATIONAL]

### A2. Petz 1988
- **Authors**: D. Petz
- **Title**: "Sufficiency of channels over von Neumann algebras"
- **Journal**: Q. J. Math. **39**, 97 (1988)
- **arXiv**: --
- **Cited in**: [P1] [P1b] [P2]
- **Key result**: Petz sufficiency theorem: exact recovery iff DPI for relative entropy is saturated.
- **Status**: [FOUNDATIONAL]

### A3. Fawzi-Renner 2015
- **Authors**: O. Fawzi and R. Renner
- **Title**: "Quantum conditional mutual information and approximate Markov chains"
- **Journal**: Commun. Math. Phys. **340**, 575 (2015)
- **arXiv**: arXiv:1410.0664
- **Cited in**: [P1] [P1b] [P2]
- **Key result**: Strengthened DPI: I(A;E|B) >= -log F^2, connecting CMI to recovery fidelity.
- **Status**: [FOUNDATIONAL]

### A4. Junge-Renner-Sutter-Wilde-Winter (JRSWW) 2018
- **Authors**: M. Junge, R. Renner, D. Sutter, M. M. Wilde, and A. Winter
- **Title**: "Universal recovery maps and approximate sufficiency of quantum relative entropy"
- **Journal**: Ann. Henri Poincare **19**, 2955 (2018)
- **arXiv**: arXiv:1509.07127
- **Cited in**: [P1] [P1b] [P2]
- **Key result**: Universal (rotated) Petz map satisfies F^2 >= exp(-Delta D); the key quantitative bound in our framework.
- **Status**: [FOUNDATIONAL]

### A5. Parzygnat-Buscemi 2023
- **Authors**: A. J. Parzygnat and F. Buscemi
- **Title**: "Axioms for retrodiction: Achieving time-reversal symmetry with a prior"
- **Journal**: Quantum **7**, 1013 (2023)
- **arXiv**: arXiv:2301.08714
- **Cited in**: [P1] [P1b]
- **Key result**: Petz map is the UNIQUE retrodiction functor satisfying Bayesian consistency -- the categorical foundation of Theorem 1.
- **Status**: [FOUNDATIONAL]

### A6. Bai-Buscemi-Scarani 2024
- **Authors**: G. Bai, F. Buscemi, and V. Scarani
- **Title**: "Fully quantum stochastic entropy production"
- **Journal**: arXiv:2412.12489 (2024)
- **arXiv**: arXiv:2412.12489
- **Cited in**: [P1] [P2]
- **Key result**: Independent work developing stochastic entropy production on the same Petz-retrodiction foundation.
- **Status**: [INDEPENDENT WORK]

### A7. Parzygnat-Fullwood 2023
- **Authors**: A. J. Parzygnat and J. Fullwood
- **Title**: "From time-reversal symmetry to quantum Bayes' rules"
- **Journal**: PRX Quantum **4**, 020334 (2023)
- **arXiv**: --
- **Cited in**: [P1]
- **Key result**: Bayesian framework for quantum channels; complements Theorem 1.
- **Status**: [SUPPORTS]

### A8. Li-Pautrat-Rouze 2025
- **Authors**: K. Li, Y. Pautrat, and C. Rouze
- **Title**: "Optimality Condition for the Petz Recovery Map"
- **Journal**: Phys. Rev. Lett. **134**, 200602 (2025)
- **arXiv**: --
- **Cited in**: [P2]
- **Key result**: Gives optimality conditions for Petz map (but optimality != saturation); constrains Paper 2's saturation hypothesis.
- **Status**: [CONSTRAINS]

### A9. Chen-Penington-Salton 2020
- **Authors**: C.-F. Chen, G. Penington, and G. Salton
- **Title**: "Entanglement wedge reconstruction using the Petz map"
- **Journal**: J. High Energy Phys. **2020**(01), 168 (2020)
- **arXiv**: --
- **Cited in**: [P1] [P2]
- **Key result**: Petz map used for holographic bulk reconstruction; connects our framework to AdS/CFT.
- **Status**: [SUPPORTS]

### A10. Cotler-Hayden-Penington-Salton-Swingle-Walter 2019
- **Authors**: J. Cotler, P. Hayden, G. Penington, G. Salton, B. Swingle, and M. Walter
- **Title**: "Entanglement Wedge Reconstruction via Universal Recovery Channels"
- **Journal**: Phys. Rev. X **9**, 031011 (2019)
- **arXiv**: --
- **Cited in**: [P2] [P4]
- **Key result**: Hayden-Preskill = Petz recovery in holography; universal recovery channels for bulk reconstruction.
- **Status**: [SUPPORTS]

### A11. Barnum-Knill 2002
- **Authors**: H. Barnum and E. Knill
- **Title**: "Reversing quantum dynamics with near-optimal quantum and classical fidelity"
- **Journal**: J. Math. Phys. **43**, 2097 (2002)
- **arXiv**: --
- **Cited in**: [P1]
- **Key result**: Petz map achieves F_e >= (1/2) F_{e,opt}^2; near-optimality as QEC decoder.
- **Status**: [FOUNDATIONAL]

### A12. Kwon-Kim 2019
- **Authors**: H. Kwon and M. S. Kim
- **Title**: "Fluctuation theorems for a quantum channel"
- **Journal**: Phys. Rev. X **9**, 031029 (2019)
- **arXiv**: --
- **Cited in**: [P1]
- **Key result**: Petz map is precisely the reverse channel in the quantum Crooks fluctuation theorem.
- **Status**: [FOUNDATIONAL]

### A13. Aw-Zaw-Balanzo-Juando-Scarani 2024
- **Authors**: C. C. Aw, L. H. Zaw, M. Balanzo-Juando, and V. Scarani
- **Title**: "Role of dilations in reversing physical processes: Tabletop reversibility and generalized thermal operations"
- **Journal**: PRX Quantum **5**, 010332 (2024)
- **arXiv**: --
- **Cited in**: [P1]
- **Key result**: Physical reversal of irreversible processes via dilations; connects to Petz map.
- **Status**: [SUPPORTS]

### A14. Wilde 2015 (Rotated Petz)
- **Authors**: M. M. Wilde
- **Title**: (Rotated Petz map construction)
- **Journal**: (Various, referenced in supplemental)
- **arXiv**: --
- **Cited in**: [P1] (supplemental)
- **Key result**: Without functoriality axiom, rotated Petz maps are alternative candidates.
- **Status**: [FOUNDATIONAL]

---

## B. Quantum Error Correction

### B1. Knill-Laflamme 1997
- **Authors**: E. Knill and R. Laflamme
- **Title**: "Theory of quantum error-correcting codes"
- **Journal**: Phys. Rev. A **55**, 900 (1997)
- **arXiv**: --
- **Cited in**: [P1]
- **Key result**: Knill-Laflamme conditions for exact QEC; the standard our Theorem 2 generalizes.
- **Status**: [FOUNDATIONAL]

### B2. Google Willow / Acharya et al. 2025
- **Authors**: R. Acharya et al. (Google Quantum AI)
- **Title**: "Quantum error correction below the surface code threshold"
- **Journal**: Nature **638**, 920 (2025)
- **arXiv**: --
- **Cited in**: [P1]
- **Key result**: Experimental demonstration of QEC below threshold; our Prediction 1 is testable on this hardware.
- **Status**: [SUPPORTS]

### B3. AlphaQubit / Bausch-Senior et al. 2024
- **Authors**: J. Bausch, A. W. Senior et al.
- **Title**: "Learning high-accuracy error decoding for quantum processors"
- **Journal**: Nature **635**, 834 (2024)
- **arXiv**: --
- **Cited in**: [P1]
- **Key result**: Neural network decoder; fits into our Prediction 2 decoder hierarchy.
- **Status**: [SUPPORTS]

### B4. de Marti iOlius et al. 2024
- **Authors**: A. de Marti iOlius et al.
- **Title**: "Decoding algorithms for surface codes"
- **Journal**: Quantum **8**, 1498 (2024)
- **arXiv**: --
- **Cited in**: [P1]
- **Key result**: Comprehensive decoder comparison; the hierarchy ML > NN > MWPM > UF supports our Prediction 2.
- **Status**: [SUPPORTS]

### B5. Smith-Brown-Bartlett 2024
- **Authors**: S. C. Smith, B. J. Brown, and S. D. Bartlett
- **Title**: "Mitigating errors in logical qubits"
- **Journal**: Commun. Phys. **7**, 386 (2024)
- **arXiv**: --
- **Cited in**: [P1]
- **Key result**: Post-selection gains in QEC; supports our Prediction 1.
- **Status**: [SUPPORTS]

### B6. English-Williamson-Bartlett 2025
- **Authors**: L. H. English, D. J. Williamson, and S. D. Bartlett
- **Title**: "Thresholds for post-selected quantum error correction from statistical mechanics"
- **Journal**: Phys. Rev. Lett. **135**, 120603 (2025)
- **arXiv**: --
- **Cited in**: [P1]
- **Key result**: Statistical-mechanical thresholds for post-selected QEC; supports Prediction 1.
- **Status**: [SUPPORTS]

### B7. English-Roberts-Bartlett-Doherty-Williamson 2025
- **Authors**: L. H. English, S. Roberts, S. D. Bartlett, A. C. Doherty, and D. J. Williamson
- **Title**: "Ising on the donut: Regimes of topological quantum error correction from statistical mechanics"
- **Journal**: arXiv:2512.10399 (2025)
- **arXiv**: arXiv:2512.10399
- **Cited in**: [P1]
- **Key result**: Boltzmann-factor failure rates via Ising-QEC mapping; supports thermodynamic interpretation.
- **Status**: [SUPPORTS]

### B8. Chen-Xu-Sommers-Huse-Thompson-Gopalakrishnan 2025
- **Authors**: H. Chen, D. Xu, G. M. Sommers, D. A. Huse, J. D. Thompson, and S. Gopalakrishnan
- **Title**: "Scalable accuracy gains from postselection in quantum error correcting codes"
- **Journal**: arXiv:2510.05222 (2025)
- **arXiv**: arXiv:2510.05222
- **Cited in**: [P1]
- **Key result**: Scalable post-selection gains; supports Prediction 1.
- **Status**: [SUPPORTS]

### B9. Singh NMR 2025
- **Authors**: Singh et al.
- **Title**: (Petz recovery in NMR)
- **Journal**: (2025)
- **arXiv**: --
- **Cited in**: [P1] [P2]
- **Key result**: Experimental implementation of Petz recovery in NMR; can be reanalyzed to test F >= exp(-Sigma/2).
- **Status**: [SUPPORTS]

### B10. Png-Scarani Ion Trap 2025
- **Authors**: W.-H. Png and V. Scarani
- **Title**: "Petz recovery maps of single-qubit decoherence channels in an ion trap quantum processor"
- **Journal**: Phys. Rev. A **112**, 022613 (2025)
- **arXiv**: --
- **Cited in**: [P1] [P2]
- **Key result**: Experimental implementation of Petz recovery in ion trap; can be reanalyzed to test F >= exp(-Sigma/2).
- **Status**: [SUPPORTS]

### B11. Almheiri-Dong-Harlow 2015
- **Authors**: A. Almheiri, X. Dong, and D. Harlow
- **Title**: "Bulk Locality and Quantum Error Correction in AdS/CFT"
- **Journal**: J. High Energy Phys. **2015**(04), 163 (2015)
- **arXiv**: arXiv:1411.7041
- **Cited in**: [P2]
- **Key result**: QEC structure in holography; bulk reconstruction = error correction.
- **Status**: [FOUNDATIONAL]

---

## C. Thermodynamics & Entropy Production

### C1. Crooks 1999
- **Authors**: G. E. Crooks
- **Title**: "Entropy production fluctuation theorem and the nonequilibrium work relation for free energy differences"
- **Journal**: Phys. Rev. E **60**, 2721 (1999)
- **arXiv**: --
- **Cited in**: [P1] [P2]
- **Key result**: Crooks fluctuation theorem: P_F(w)/P_R(-w) = exp(beta(w - Delta F)); the Petz map is the reverse channel.
- **Status**: [FOUNDATIONAL]

### C2. Jarzynski 1997
- **Authors**: C. Jarzynski
- **Title**: "Nonequilibrium equality for free energy differences"
- **Journal**: Phys. Rev. Lett. **78**, 2690 (1997)
- **arXiv**: --
- **Cited in**: [P1] (implicit)
- **Key result**: Jarzynski equality; precursor to Crooks; connects free energy to nonequilibrium work.
- **Status**: [FOUNDATIONAL]

### C3. Landi-Paternostro 2021
- **Authors**: G. T. Landi and M. Paternostro
- **Title**: "Irreversible entropy production: From classical to quantum"
- **Journal**: Rev. Mod. Phys. **93**, 035008 (2021)
- **arXiv**: --
- **Cited in**: [P1] [P2]
- **Key result**: Comprehensive review of entropy production; confirms the framework connecting Sigma to irreversibility.
- **Status**: [FOUNDATIONAL]

### C4. Landauer 1961
- **Authors**: R. Landauer
- **Title**: "Irreversibility and heat generation in the computing process"
- **Journal**: IBM J. Res. Dev. **5**, 183 (1961)
- **arXiv**: --
- **Cited in**: [P1] (implicit via Retrodiction Landauer Principle)
- **Key result**: Minimum energy cost of erasing one bit = k_B T ln 2; our framework extends this to retrodiction.
- **Status**: [FOUNDATIONAL]

### C5. Kolchinsky-Wolpert 2017
- **Authors**: A. Kolchinsky and D. H. Wolpert
- **Title**: "Dependence of dissipation on the initial distribution over states"
- **Journal**: J. Stat. Mech. **2017**, 083202 (2017)
- **arXiv**: --
- **Cited in**: [P1]
- **Key result**: Decomposition of entropy production into Landauer, mismatch, and residual components; used in Theorem 3.
- **Status**: [FOUNDATIONAL]

### C6. Reeb-Wolf 2014
- **Authors**: D. Reeb and M. M. Wolf
- **Title**: "An improved Landauer principle with finite-size corrections"
- **Journal**: New J. Phys. **16**, 103011 (2014)
- **arXiv**: --
- **Cited in**: [P1]
- **Key result**: Improved Landauer bound: Delta D >= Sigma >= 0; used for Step 4 of Theorem 3.
- **Status**: [FOUNDATIONAL]

### C7. Herrera 2020
- **Authors**: L. Herrera
- **Title**: "Landauer Principle and General Relativity"
- **Journal**: Entropy **22**, 340 (2020)
- **arXiv**: --
- **Cited in**: [P2]
- **Key result**: Gravitational Landauer principle; one route to derive Sigma_grav = r_s/r.
- **Status**: [SUPPORTS]

### C8. Tolman 1930
- **Authors**: R. C. Tolman
- **Title**: "On the Weight of Heat and Thermal Equilibrium in General Relativity"
- **Journal**: Phys. Rev. **35**, 904 (1930)
- **arXiv**: --
- **Cited in**: [P2]
- **Key result**: Tolman relation: T(r) = T_inf / sqrt(-g_00); temperature redshifts in gravitational fields.
- **Status**: [FOUNDATIONAL]

### C9. Basso-Maziero-Celeri 2025
- **Authors**: M. L. W. Basso, J. Maziero, and L. C. Celeri
- **Title**: "Quantum Detailed Fluctuation Theorem in Curved Spacetimes"
- **Journal**: Phys. Rev. Lett. **134**, 050406 (2025)
- **arXiv**: arXiv:2405.03902
- **Cited in**: [P2] [P4]
- **Key result**: Fully general-relativistic quantum Crooks relation; entropy production from curvature coupling.
- **Status**: [SUPPORTS]

---

## D. Decoherence & Quantum Foundations

### D1. Joos-Zeh 1985
- **Authors**: E. Joos and H. D. Zeh
- **Title**: "The emergence of classical properties through interaction with the environment"
- **Journal**: Z. Phys. B **59**, 223 (1985)
- **arXiv**: --
- **Cited in**: [P1b]
- **Key result**: Original calculation of environment-induced decoherence; tensor-product amplification mechanism.
- **Status**: [FOUNDATIONAL]

### D2. Zeh 1970
- **Authors**: H. D. Zeh
- **Title**: "On the interpretation of measurement in quantum theory"
- **Journal**: Found. Phys. **1**, 69 (1970)
- **arXiv**: --
- **Cited in**: [P1]
- **Key result**: Pioneer of decoherence approach; environment selects preferred basis.
- **Status**: [FOUNDATIONAL]

### D3. Zurek 2003
- **Authors**: W. H. Zurek
- **Title**: "Decoherence, einselection, and the quantum origins of the classical"
- **Journal**: Rev. Mod. Phys. **75**, 715 (2003)
- **arXiv**: --
- **Cited in**: [P1] [P1b]
- **Key result**: Comprehensive review of decoherence and einselection; pointer basis selection by interaction Hamiltonian.
- **Status**: [FOUNDATIONAL]

### D4. Schlosshauer 2004/2007
- **Authors**: M. Schlosshauer
- **Title**: "Decoherence, the measurement problem, and interpretations of quantum mechanics"
- **Journal**: Rev. Mod. Phys. **76**, 1267 (2004); book (2007)
- **arXiv**: --
- **Cited in**: [P1b]
- **Key result**: Modern review of decoherence theory; FAPP vs genuine collapse distinction.
- **Status**: [FOUNDATIONAL]

### D5. Zurek 2009 (Quantum Darwinism)
- **Authors**: W. H. Zurek
- **Title**: "Quantum Darwinism"
- **Journal**: Nat. Phys. **5**, 181 (2009)
- **arXiv**: --
- **Cited in**: [P1b]
- **Key result**: Redundant encoding of information in environment; R measures objectivity.
- **Status**: [FOUNDATIONAL]

### D6. Korbicz 2021 (SBS review)
- **Authors**: J. K. Korbicz
- **Title**: "Roads to objectivity: Quantum Darwinism, spectrum broadcast structures, and strong quantum Darwinism -- a review"
- **Journal**: Quantum **5**, 571 (2021)
- **arXiv**: --
- **Cited in**: [P1b]
- **Key result**: SBS theorem: von Neumann interaction with N >> 1 sub-environments generically produces objective states.
- **Status**: [FOUNDATIONAL]

### D7. Brandao-Piani-Horodecki
- **Authors**: F. G. S. L. Brandao, M. Piani, and P. Horodecki
- **Title**: (Genericity of classical features)
- **Journal**: (Various)
- **arXiv**: --
- **Cited in**: [P1b] (implicit)
- **Key result**: ANY quantum dynamics with N >> 1 sub-environments produces classical features generically.
- **Status**: [SUPPORTS]

### D8. Riedel-Zurek 2017/2021
- **Authors**: C. J. Riedel
- **Title**: "Decoherence from classicality: Gravitational decoherence and redundant records"
- **Journal**: in *Do Wave Functions Jump?* (Springer, 2021), pp. 315-325
- **arXiv**: --
- **Cited in**: [P1b]
- **Key result**: Redundant encoding in gravitational decoherence confirmed; supports tau-R bridge.
- **Status**: [SUPPORTS]

### D9. Zhu et al. 2025
- **Authors**: Y. Zhu et al.
- **Title**: "Observation of quantum Darwinism in a photonic system"
- **Journal**: Sci. Adv. **11**, eadr8585 (2025)
- **arXiv**: --
- **Cited in**: [P1b]
- **Key result**: Experimental observation of quantum Darwinism; testbed for tau-R bridge predictions.
- **Status**: [SUPPORTS]

### D10. Unden et al. 2019
- **Authors**: T. K. Unden et al.
- **Title**: "Revealing the emergence of classicality using nitrogen-vacancy centers"
- **Journal**: Phys. Rev. Lett. **123**, 140402 (2019)
- **arXiv**: --
- **Cited in**: [P1b]
- **Key result**: NV center experiment revealing classicality emergence; relevant to tau-R bridge testing.
- **Status**: [SUPPORTS]

### D11. Englert 1996
- **Authors**: B.-G. Englert
- **Title**: "Fringe visibility and which-way information: An inequality"
- **Journal**: Phys. Rev. Lett. **77**, 2154 (1996)
- **arXiv**: --
- **Cited in**: [P1b]
- **Key result**: V^2 + D^2 <= 1 wave-particle duality relation; rewritten in tau language in Paper 1b.
- **Status**: [FOUNDATIONAL]

### D12. Greenberger-Yasin 1988
- **Authors**: D. M. Greenberger and A. Yasin
- **Title**: "Simultaneous wave and particle knowledge in a neutron interferometer"
- **Journal**: Phys. Lett. A **128**, 391 (1988)
- **arXiv**: --
- **Cited in**: [P1b]
- **Key result**: Original wave-particle duality inequality; complementary to Englert.
- **Status**: [FOUNDATIONAL]

### D13. BRS 2007 (Superselection = Reference Frames)
- **Authors**: S. D. Bartlett, T. Rudolph, and R. W. Spekkens
- **Title**: "Reference frames, superselection rules, and quantum information"
- **Journal**: Rev. Mod. Phys. **79**, 555 (2007)
- **arXiv**: --
- **Cited in**: [P1b]
- **Key result**: Superselection rules = lack of reference frame; tau = 1 implies operational superselection.
- **Status**: [FOUNDATIONAL]

### D14. ABL / Aharonov-Bergmann-Lebowitz 1964
- **Authors**: Y. Aharonov, P. G. Bergmann, and J. L. Lebowitz
- **Title**: "Time symmetry in the quantum process of measurement"
- **Journal**: Phys. Rev. **134**, B1410 (1964)
- **arXiv**: --
- **Cited in**: [P1]
- **Key result**: Two-state vector formalism (TSVF); Theorem 1 provides unique choice for backward state.
- **Status**: [FOUNDATIONAL]

### D15. Scully-Druhl 1982
- **Authors**: M. O. Scully and K. Druhl
- **Title**: "Quantum eraser: A proposed photon correlation experiment concerning observation and 'delayed choice' in quantum mechanics"
- **Journal**: Phys. Rev. A **25**, 2208 (1982)
- **arXiv**: --
- **Cited in**: [P1]
- **Key result**: Original quantum eraser proposal; the physical starting point of our framework.
- **Status**: [FOUNDATIONAL]

### D16. Kim-Yu-Kulik-Shih-Scully 2000
- **Authors**: Y.-H. Kim, R. Yu, S. P. Kulik, Y. Shih, and M. O. Scully
- **Title**: "Delayed 'choice' quantum eraser"
- **Journal**: Phys. Rev. Lett. **84**, 1 (2000)
- **arXiv**: --
- **Cited in**: [P1]
- **Key result**: Experimental realization of delayed-choice quantum eraser.
- **Status**: [FOUNDATIONAL]

### D17. Wheeler 1978
- **Authors**: J. A. Wheeler
- **Title**: "The 'past' and the 'delayed-choice' double-slit experiment"
- **Journal**: in *Mathematical Foundations of Quantum Theory* (Academic Press, 1978), pp. 9-48
- **arXiv**: --
- **Cited in**: [P1]
- **Key result**: Delayed-choice experiment; in our framework, tau = 0 for the closed signal-idler pair explains the result without retrocausality.
- **Status**: [FOUNDATIONAL]

### D18. Breuer-Laine-Piilo-Vacchini 2016
- **Authors**: H.-P. Breuer, E. M. Laine, J. Piilo, and B. Vacchini
- **Title**: "Colloquium: Non-Markovian dynamics in open quantum systems"
- **Journal**: Rev. Mod. Phys. **88**, 021002 (2016)
- **arXiv**: --
- **Cited in**: [P1]
- **Key result**: Non-Markovian processes can have temporarily negative Sigma; implies transient tau decreases.
- **Status**: [SUPPORTS]

### D19. Pollock-Rodriguez-Rosario-Frauenheim-Paternostro-Modi 2018
- **Authors**: F. A. Pollock, C. Rodriguez-Rosario, T. Frauenheim, M. Paternostro, and K. Modi
- **Title**: "Non-Markovian quantum processes: Complete framework and efficient characterization"
- **Journal**: Phys. Rev. A **97**, 012127 (2018)
- **arXiv**: --
- **Cited in**: [P1]
- **Key result**: Process tensor framework; may extend equivalence chain to non-Markovian regime.
- **Status**: [SUPPORTS]

### D20. Proietti et al. 2019 (Wigner's friend)
- **Authors**: M. Proietti et al.
- **Title**: (Wigner's friend experiment)
- **Journal**: (2019)
- **arXiv**: --
- **Cited in**: [P1b] (implicit)
- **Key result**: Experimental Wigner's friend scenario; tau is partition-dependent in such scenarios.
- **Status**: [SUPPORTS]

---

## E. Gravitational Decoherence & Collapse

### E1. Penrose 1996
- **Authors**: R. Penrose
- **Title**: "On gravity's role in quantum state reduction"
- **Journal**: Gen. Relativ. Gravit. **28**, 581 (1996)
- **arXiv**: --
- **Cited in**: [P1b]
- **Key result**: Gravitational self-energy E_G sets fundamental collapse timescale t_c = hbar/E_G.
- **Status**: [FOUNDATIONAL]

### E2. Diosi 1987
- **Authors**: L. Diosi
- **Title**: "A universal master equation for the gravitational violation of quantum mechanics"
- **Journal**: Phys. Lett. A **120**, 377 (1987)
- **arXiv**: --
- **Cited in**: [P1b]
- **Key result**: Universal master equation for gravitational collapse; collapse rate Lambda = E_G/hbar.
- **Status**: [FOUNDATIONAL]

### E3. Pikovski-Zych-Costa-Brukner 2015
- **Authors**: I. Pikovski, M. Zych, F. Costa, and C. Brukner
- **Title**: "Universal decoherence due to gravitational time dilation"
- **Journal**: Nat. Phys. **11**, 668 (2015)
- **arXiv**: arXiv:1311.1095
- **Cited in**: [P1b] [P2]
- **Key result**: Gravitational time dilation causes universal dephasing of internal DOF; the specific channel analyzed in Paper 1b.
- **Status**: [FOUNDATIONAL]

### E4. Galley-Giacomini-Selby 2023
- **Authors**: T. D. Galley, F. Giacomini, and J. H. Selby
- **Title**: "Any consistent coupling between classical gravity and quantum matter is fundamentally irreversible"
- **Journal**: Quantum **7**, 1142 (2023)
- **arXiv**: --
- **Cited in**: [P1b]
- **Key result**: No Stinespring dilation for classical-quantum gravity coupling; if gravity classical, tau > 0 is genuine.
- **Status**: [FOUNDATIONAL]

### E5. DSW / Danielson-Satishchandran-Wald 2022
- **Authors**: D. L. Danielson, G. Satishchandran, and R. M. Wald
- **Title**: "Gravitationally mediated entanglement: Newtonian field versus gravitons"
- **Journal**: Phys. Rev. D **105**, 086001 (2022)
- **arXiv**: --
- **Cited in**: [P1b]
- **Key result**: Superpositions of masses radiate soft gravitons crossing horizons; genuine information loss.
- **Status**: [SUPPORTS]

### E6. DSW 2025
- **Authors**: D. L. Danielson, G. Satishchandran, and R. M. Wald
- **Title**: (Extended gravitational decoherence results)
- **Journal**: (2025)
- **arXiv**: --
- **Cited in**: [P1b]
- **Key result**: Extended results on soft graviton radiation and causal irreversibility.
- **Status**: [SUPPORTS]

### E7. Kiefer-Joos 1999
- **Authors**: C. Kiefer and E. Joos
- **Title**: "Decoherence: Concepts and examples"
- **Journal**: in *Quantum Future*, Lecture Notes in Physics Vol. 517 (Springer, 1999), pp. 105-128
- **arXiv**: --
- **Cited in**: [P1b]
- **Key result**: With infinitely many graviton field modes, decoherence is permanent (no Poincare recurrence).
- **Status**: [SUPPORTS]

### E8. Oppenheim 2023
- **Authors**: J. Oppenheim
- **Title**: "A postquantum theory of classical gravity?"
- **Journal**: Phys. Rev. X **13**, 041040 (2023)
- **arXiv**: --
- **Cited in**: [P1b]
- **Key result**: Post-quantum classical gravity is consistent and testable; if gravity classical, tau > 0 is genuine.
- **Status**: [SUPPORTS]

### E9. Oppenheim et al. 2023
- **Authors**: J. Oppenheim, C. Sparaciari, B. Soda, and Z. Weller-Davies
- **Title**: (Decoherence-diffusion tradeoff)
- **Journal**: Nat. Commun. **14**, 7910 (2023)
- **arXiv**: --
- **Cited in**: [P1b] (implicit)
- **Key result**: Constraints on classical-quantum gravity from decoherence-diffusion tradeoff.
- **Status**: [CONSTRAINS]

### E10. Grossardt et al. 2025
- **Authors**: A. Grossardt et al.
- **Title**: (Diosi-Penrose model with classical gravity CAN generate entanglement)
- **Journal**: Phys. Rev. D **111**, L121101 (2025)
- **arXiv**: --
- **Cited in**: [P1b] (implicit)
- **Key result**: BMV experiment alone cannot distinguish classical from quantum gravity; blurs the boundary.
- **Status**: [CONSTRAINS]

### E11. Donadi et al. 2021 (Gran Sasso)
- **Authors**: S. Donadi et al.
- **Title**: "Underground test of gravity-related wave function collapse"
- **Journal**: Nat. Phys. **17**, 74 (2021)
- **arXiv**: --
- **Cited in**: [P1b]
- **Key result**: Rules out parameter-free Diosi model; Penrose version (with nuclear-scale smearing) remains viable.
- **Status**: [CONSTRAINS]

### E12. Bassi et al. 2013
- **Authors**: A. Bassi, K. Lochan, S. Satin, T. P. Singh, and H. Ulbricht
- **Title**: "Models of wave-function collapse, underlying theories, and experimental tests"
- **Journal**: Rev. Mod. Phys. **85**, 471 (2013)
- **arXiv**: --
- **Cited in**: [P1b] (implicit)
- **Key result**: Comprehensive review of collapse models; GRW, CSL, Diosi-Penrose.
- **Status**: [FOUNDATIONAL]

### E13. GRW 1986
- **Authors**: G. C. Ghirardi, A. Rimini, and T. Weber
- **Title**: "Unified dynamics for microscopic and macroscopic systems"
- **Journal**: Phys. Rev. D **34**, 470 (1986)
- **arXiv**: --
- **Cited in**: [P1b] (implicit)
- **Key result**: Spontaneous collapse model; tau framework offers alternative operational approach.
- **Status**: [FOUNDATIONAL]

### E14. Hossenfelder 2025
- **Authors**: S. Hossenfelder
- **Title**: (Product state constraint leads to collapse)
- **Journal**: arXiv:2510.11037 (2025)
- **arXiv**: arXiv:2510.11037
- **Cited in**: [P1b] (implicit)
- **Key result**: Product state constraint as alternative route to collapse.
- **Status**: [SUPPORTS]

### E15. Artini-Rufo-Wald 2025
- **Authors**: S. Artini, F. Rufo, and R. Wald
- **Title**: "Entropy production in the Diosi-Penrose model"
- **Journal**: Phys. Rev. Res. (2025)
- **arXiv**: --
- **Cited in**: [P1b]
- **Key result**: Entropy production computed for Diosi-Penrose model; connects directly to our Sigma.
- **Status**: [SUPPORTS]

### E16. Moreira-Celeri 2025
- **Authors**: L. C. Moreira and L. C. Celeri
- **Title**: "Graviton-bath entropy production and decoherence"
- **Journal**: arXiv:2407.21186 (2025)
- **arXiv**: arXiv:2407.21186
- **Cited in**: [P2]
- **Key result**: Graviton bath entropy production framework; complementary to Pikovski.
- **Status**: [SUPPORTS]

### E17. Bonder et al. 2016
- **Authors**: Y. Bonder et al.
- **Title**: (Critique of Pikovski)
- **Journal**: (2016)
- **arXiv**: --
- **Cited in**: [P1b] (implicit)
- **Key result**: Critiques of Pikovski mechanism; must be addressed honestly in Paper 1b.
- **Status**: [POTENTIAL CONFLICT]

---

## F. Experiments (Interferometry & Optomechanics)

### F1. Fein et al. 2019
- **Authors**: Y. Y. Fein, P. Geyer, P. Zwick, F. Kialka, S. Pedalino, M. Mayor, S. Gerlich, and M. Arndt
- **Title**: "Quantum superposition of molecules beyond 25 kDa"
- **Journal**: Nat. Phys. **15**, 1242 (2019)
- **arXiv**: --
- **Cited in**: [P1b]
- **Key result**: Largest molecule interferometry (25 kDa); tau_grav << 10^-14 (consistent with predictions).
- **Status**: [SUPPORTS]

### F2. Arndt 2025
- **Authors**: M. Arndt et al.
- **Title**: "Quantum interference of large organic molecules" (170 kDa)
- **Journal**: (2025)
- **arXiv**: --
- **Cited in**: [P1b]
- **Key result**: Extended mass range for molecule interferometry; still tau_grav << 1.
- **Status**: [SUPPORTS]

### F3. Delic et al. 2020
- **Authors**: U. Delic, M. Reisenbauer, K. Dare, D. Grass, V. Vuletic, N. Kiesel, and M. Aspelmeyer
- **Title**: "Cooling of a levitated nanoparticle to the motional quantum ground state"
- **Journal**: Science **367**, 892 (2020)
- **arXiv**: --
- **Cited in**: [P1b]
- **Key result**: Ground-state cooling of nanoparticle; stepping stone to gravitational decoherence tests.
- **Status**: [SUPPORTS]

### F4. Panda/Muller 2024
- **Authors**: C. D. Panda et al.
- **Title**: "Measuring gravitational decoherence with atom interferometry"
- **Journal**: (2024)
- **arXiv**: --
- **Cited in**: [P1b]
- **Key result**: Atom interferometry for gravitational decoherence; tau_grav << 10^-19.
- **Status**: [SUPPORTS]

### F5. Bose et al. 2017 (BMV proposal)
- **Authors**: S. Bose, A. Mazumdar, G. W. Morley, H. Ulbricht, M. Toros, M. Paternostro, A. A. Geraci, P. F. Barker, M. S. Kim, and G. Milburn
- **Title**: "Spin entanglement witness for quantum gravity"
- **Journal**: Phys. Rev. Lett. **119**, 240401 (2017)
- **arXiv**: --
- **Cited in**: [P1b]
- **Key result**: Proposed experiment to detect gravitationally-induced entanglement; tau_grav ~ 0.004 regime.
- **Status**: [SUPPORTS]

### F6. Marletto-Vedral 2017
- **Authors**: C. Marletto and V. Vedral
- **Title**: "Gravitationally induced entanglement between two massive particles is sufficient evidence of quantum effects in gravity"
- **Journal**: Phys. Rev. Lett. **119**, 240402 (2017)
- **arXiv**: --
- **Cited in**: [P1b]
- **Key result**: Complementary BMV proposal; gravitational entanglement witnesses quantum gravity.
- **Status**: [SUPPORTS]

### F7. Wang et al. 2002 (Sigma < 0 observation)
- **Authors**: G. M. Wang, E. M. Sevick, E. Mittag, D. J. Searles, and D. J. Evans
- **Title**: "Experimental demonstration of violations of the second law of thermodynamics for small systems and short time scales"
- **Journal**: Phys. Rev. Lett. **89**, 050601 (2002)
- **arXiv**: --
- **Cited in**: [P1] (implicit)
- **Key result**: First experimental observation of Sigma < 0 (negative entropy production) in colloidal particles.
- **Status**: [SUPPORTS]

### F8. Collin et al. 2005 (RNA)
- **Authors**: D. Collin et al.
- **Title**: (RNA hairpin Crooks fluctuation theorem test)
- **Journal**: Nature **437**, 231 (2005)
- **arXiv**: --
- **Cited in**: [P1] (implicit)
- **Key result**: Crooks theorem verified for RNA hairpin folding/unfolding; demonstrates tau_signed < 0 is physical.
- **Status**: [SUPPORTS]

### F9. Micadei et al. 2019
- **Authors**: K. Micadei et al.
- **Title**: (Heat flowing from cold to hot driven by quantum correlations)
- **Journal**: Nat. Commun. **10**, 2456 (2019)
- **arXiv**: --
- **Cited in**: [P1] (implicit)
- **Key result**: Heat flows spontaneously from cold to hot using quantum entanglement; demonstrates Sigma < 0 driven by correlations.
- **Status**: [SUPPORTS]

### F10. Bertotti-Iess-Tortora 2003 (Cassini)
- **Authors**: B. Bertotti, L. Iess, and P. Tortora
- **Title**: "A test of general relativity using radio links with the Cassini spacecraft"
- **Journal**: Nature **425**, 374 (2003)
- **arXiv**: --
- **Cited in**: [P2]
- **Key result**: Most precise test of Shapiro delay; gamma = 1 + (2.1 +/- 2.3) x 10^-5; constrains Paper 2.
- **Status**: [CONSTRAINS]

---

## G. General Relativity & Black Holes

### G1. Einstein 1911
- **Authors**: A. Einstein
- **Title**: "Uber den Einfluss der Schwerkraft auf die Ausbreitung des Lichtes"
- **Journal**: Ann. Phys. (Leipzig) **35**, 898 (1911)
- **arXiv**: --
- **Cited in**: [P2]
- **Key result**: First prediction of gravitational light bending; predicted exponential metric before GR.
- **Status**: [FOUNDATIONAL]

### G2. Einstein 1915
- **Authors**: A. Einstein
- **Title**: "Die Feldgleichungen der Gravitation"
- **Journal**: Sitzungsber. Preuss. Akad. Wiss. Berlin, 844 (1915)
- **arXiv**: --
- **Cited in**: [P2] (implicit)
- **Key result**: General relativity field equations.
- **Status**: [FOUNDATIONAL]

### G3. Dicke 1957
- **Authors**: R. H. Dicke
- **Title**: "Gravitation without a Principle of Equivalence"
- **Journal**: Rev. Mod. Phys. **29**, 363 (1957)
- **arXiv**: --
- **Cited in**: [P2]
- **Key result**: Alternative to equivalence principle; relevant to exponential metric history.
- **Status**: [FOUNDATIONAL]

### G4. Jacobson 1995
- **Authors**: T. Jacobson
- **Title**: "Thermodynamics of Spacetime: The Einstein Equation of State"
- **Journal**: Phys. Rev. Lett. **75**, 1260 (1995)
- **arXiv**: gr-qc/9504004
- **Cited in**: [P2] [P3] [P4]
- **Key result**: Einstein equation = delta Q = T dS for local Rindler horizons; gravity from thermodynamics.
- **Status**: [FOUNDATIONAL]

### G5. Jacobson 2016
- **Authors**: T. Jacobson
- **Title**: "Entanglement Equilibrium and the Einstein Equation"
- **Journal**: Phys. Rev. Lett. **116**, 201101 (2016)
- **arXiv**: arXiv:1505.04753
- **Cited in**: [P2] [P3]
- **Key result**: Einstein equation iff vacuum entanglement entropy in small geodesic balls is maximized.
- **Status**: [FOUNDATIONAL]

### G6. Verlinde 2010
- **Authors**: E. Verlinde
- **Title**: "On the Origin of Gravity and the Laws of Newton"
- **Journal**: J. High Energy Phys. **2011**(04), 029 (2011)
- **arXiv**: arXiv:1001.0785
- **Cited in**: [P3] [P4]
- **Key result**: Gravity as entropic force: F = T(dS/dx); reproduces Newton exactly.
- **Status**: [FOUNDATIONAL]

### G7. Verlinde 2016
- **Authors**: E. Verlinde
- **Title**: "Emergent Gravity and the Dark Universe"
- **Journal**: SciPost Phys. **2**, 016 (2017)
- **arXiv**: arXiv:1611.02269
- **Cited in**: [P3] [P4]
- **Key result**: Volume-law entanglement in de Sitter -> apparent dark matter; derives MOND-like acceleration.
- **Status**: [FOUNDATIONAL]

### G8. Bekenstein 1981
- **Authors**: J. D. Bekenstein
- **Title**: (Black hole entropy / Bekenstein bound)
- **Journal**: (Various)
- **arXiv**: --
- **Cited in**: [P2] [P4]
- **Key result**: Black hole entropy proportional to area; Bekenstein bound on information.
- **Status**: [FOUNDATIONAL]

### G9. Unruh 1976
- **Authors**: W. G. Unruh
- **Title**: "Notes on black-hole evaporation"
- **Journal**: Phys. Rev. D **14**, 870 (1976)
- **arXiv**: --
- **Cited in**: [P2]
- **Key result**: Unruh effect: accelerated observer sees thermal radiation at T = a/(2 pi); connects to modular theory.
- **Status**: [FOUNDATIONAL]

### G10. Hawking 1976
- **Authors**: S. W. Hawking
- **Title**: "Breakdown of predictability in gravitational collapse"
- **Journal**: Phys. Rev. D **14**, 2460 (1976)
- **arXiv**: --
- **Cited in**: [P2]
- **Key result**: Information paradox from black hole evaporation; tau < 1 in our framework eliminates the premise.
- **Status**: [FOUNDATIONAL]

### G11. Penington 2020
- **Authors**: G. Penington
- **Title**: "Entanglement Wedge Reconstruction and the Information Problem"
- **Journal**: J. High Energy Phys. **2020**(09), 002 (2020)
- **arXiv**: --
- **Cited in**: [P2]
- **Key result**: Island formula resolves information paradox; Petz map used for reconstruction.
- **Status**: [SUPPORTS]

### G12. Almheiri et al. 2021 (AMPS)
- **Authors**: A. Almheiri, T. Hartman, J. Maldacena, E. Shaghoulian, and A. Tajdini
- **Title**: "The entropy of Hawking radiation"
- **Journal**: Rev. Mod. Phys. **93**, 035002 (2021)
- **arXiv**: --
- **Cited in**: [P2]
- **Key result**: Complete review of information paradox resolution via quantum extremal surfaces.
- **Status**: [SUPPORTS]

### G13. AMPS 2013
- **Authors**: A. Almheiri, D. Marolf, J. Polchinski, and J. Sully
- **Title**: "Black Holes: Complementarity vs. Firewalls"
- **Journal**: J. High Energy Phys. **2013**(02), 062 (2013)
- **arXiv**: --
- **Cited in**: [P2]
- **Key result**: Firewall paradox; our exponential metric (no horizon) eliminates this problem.
- **Status**: [SUPPORTS]

### G14. Dorau-Much 2025
- **Authors**: L. Dorau and A. Much
- **Title**: "From Quantum Relative Entropy to the Semiclassical Einstein Equations"
- **Journal**: Phys. Rev. Lett. (2025)
- **arXiv**: arXiv:2510.24491
- **Cited in**: [P2] [P3] [P4]
- **Key result**: Derives semiclassical Einstein equations from QRE via modular theory; deepest QI-to-gravity connection.
- **Status**: [FOUNDATIONAL]

### G15. Bianconi 2025 (PRD)
- **Authors**: G. Bianconi
- **Title**: "Gravity from entropy"
- **Journal**: Phys. Rev. D **111**, 066001 (2025)
- **arXiv**: arXiv:2408.14391
- **Cited in**: [P2] [P4]
- **Key result**: Metric = density matrix; gravity action = QRE; best candidate for unified equation framework.
- **Status**: [FOUNDATIONAL]

### G16. Bianconi 2025 (Cosmology)
- **Authors**: G. Bianconi
- **Title**: "The Thermodynamics of the Gravity from Entropy Theory: from the Hamiltonian to applications in Cosmology"
- **Journal**: arXiv:2510.22545 (2025)
- **arXiv**: arXiv:2510.22545
- **Cited in**: [P4]
- **Key result**: FRW solutions in GfE framework; dynamical cosmological constant emerges from G-field.
- **Status**: [SUPPORTS]

### G17. Bianconi 2025 (Schwarzschild Area Law)
- **Authors**: G. Bianconi
- **Title**: "The Quantum Relative Entropy of the Schwarzschild Black Hole and the Area Law"
- **Journal**: Entropy **27**(3), 266 (2025)
- **arXiv**: arXiv:2501.09491
- **Cited in**: [P4]
- **Key result**: GQRE of Schwarzschild depends on R_s/r^3 (curvature scale); area law for large R_s.
- **Status**: [SUPPORTS]

### G18. Will 2014
- **Authors**: C. M. Will
- **Title**: "The Confrontation between General Relativity and Experiment"
- **Journal**: Living Rev. Relativ. **17**, 4 (2014)
- **arXiv**: --
- **Cited in**: [P2]
- **Key result**: Comprehensive review of GR tests; PPN constraints on alternative metrics.
- **Status**: [CONSTRAINS]

---

## G (continued). Exponential Metric & Alternative Gravity

### G19. Papapetrou 1954
- **Authors**: A. Papapetrou
- **Title**: "Eine Theorie des Gravitationsfeldes mit einer Feldfunktion"
- **Journal**: Z. Phys. **139**, 518 (1954)
- **arXiv**: --
- **Cited in**: [P2]
- **Key result**: Early exponential metric proposal.
- **Status**: [FOUNDATIONAL]

### G20. Yilmaz 1958
- **Authors**: H. Yilmaz
- **Title**: "New Approach to General Relativity"
- **Journal**: Phys. Rev. **111**, 1417 (1958)
- **arXiv**: --
- **Cited in**: [P2]
- **Key result**: Alternative gravity theory producing exponential metric.
- **Status**: [FOUNDATIONAL]

### G21. Yilmaz 1971
- **Authors**: H. Yilmaz
- **Title**: "New Theory of Gravitation"
- **Journal**: Phys. Rev. Lett. **27**, 1399 (1971)
- **arXiv**: --
- **Cited in**: [P2]
- **Key result**: Extended Yilmaz theory.
- **Status**: [FOUNDATIONAL]

### G22. Rosen 1973
- **Authors**: N. Rosen
- **Title**: "A Bimetric Theory of Gravitation"
- **Journal**: Gen. Relativ. Gravit. **4**, 435 (1973)
- **arXiv**: --
- **Cited in**: [P2]
- **Key result**: Bimetric theory; related to exponential metric.
- **Status**: [FOUNDATIONAL]

### G23. Makukov-Mychelkin 2020
- **Authors**: M. A. Makukov and E. G. Mychelkin
- **Title**: "Triple Path to the Exponential Metric"
- **Journal**: Found. Phys. **50**, 1346 (2020)
- **arXiv**: arXiv:2009.08655
- **Cited in**: [P2]
- **Key result**: Three independent derivations of the exponential metric; historical overview.
- **Status**: [SUPPORTS]

### G24. Misner 1995
- **Authors**: C. W. Misner
- **Title**: "Yilmaz Cancels Newton"
- **Journal**: arXiv:gr-qc/9504050 (1995)
- **arXiv**: gr-qc/9504050
- **Cited in**: [P2]
- **Key result**: Critique of Yilmaz theory; must be addressed in Paper 2.
- **Status**: [POTENTIAL CONFLICT]

### G25. Boonserm-Ngampitipan-Simpson-Visser 2018
- **Authors**: P. Boonserm, T. Ngampitipan, A. Simpson, and M. Visser
- **Title**: "Exponential metric represents a traversable wormhole"
- **Journal**: Phys. Rev. D **98**, 084048 (2018)
- **arXiv**: arXiv:1805.03781
- **Cited in**: [P2]
- **Key result**: Exponential metric is a traversable wormhole; no event horizon.
- **Status**: [SUPPORTS]

### G26. Nath-Sarma 2025
- **Authors**: S. Nath and D. Sarma
- **Title**: "Exponential wormhole solution in logarithmic f(R) gravity"
- **Journal**: (2025)
- **arXiv**: --
- **Cited in**: [P2]
- **Key result**: Exponential metric emerges as solution in f(R) = R + alpha R ln(R/R_0).
- **Status**: [SUPPORTS]

### G27. Cardoso-Franzato-Pani 2016
- **Authors**: V. Cardoso, E. Franzato, and P. Pani
- **Title**: "Is the Gravitational-Wave Ringdown a Probe of the Event Horizon?"
- **Journal**: Phys. Rev. Lett. **116**, 171101 (2016)
- **arXiv**: --
- **Cited in**: [P2]
- **Key result**: GW echoes as test of horizonless compact objects; directly relevant to Paper 2 predictions.
- **Status**: [SUPPORTS]

### G28. Damour-Solodukhin 2007
- **Authors**: T. Damour and S. N. Solodukhin
- **Title**: "Wormholes as black hole foils"
- **Journal**: Phys. Rev. D **76**, 024016 (2007)
- **arXiv**: arXiv:0704.2667
- **Cited in**: [P2]
- **Key result**: Wormholes can mimic black holes; echo time formula.
- **Status**: [SUPPORTS]

---

## H. Dark Matter / Rotation Curves

### H1. Kumar 2025
- **Authors**: N. Kumar
- **Title**: "Marginal IR running of Gravity as a Natural Explanation for Dark Matter"
- **Journal**: Phys. Lett. B **871**, 140008 (2025)
- **arXiv**: arXiv:2509.05246
- **Cited in**: [P3]
- **Key result**: eta = 1 marginal running of G gives logarithmic potential; fits 3 galaxy rotation curves with universal r_c ~ 37 kpc.
- **Status**: [FOUNDATIONAL]

### H2. Gubitosi-Piattella-Casarini 2024
- **Authors**: G. Gubitosi, O. F. Piattella, and D. Casarini
- **Title**: "Phenomenology of renormalization group improved gravity from the kinematics of SPARC galaxies"
- **Journal**: (2024)
- **arXiv**: arXiv:2403.00531
- **Cited in**: [P3]
- **Key result**: Running G (RGGR) fits 100 SPARC galaxies, competitive with NFW halos.
- **Status**: [SUPPORTS]

### H3. Ghari-Haghi 2026
- **Authors**: A. Ghari and H. Haghi
- **Title**: (Verlinde favored over MOND at 5.2 sigma in dwarf spheroidals)
- **Journal**: (2026)
- **arXiv**: arXiv:2601.01715
- **Cited in**: [P3]
- **Key result**: Verlinde's emergent gravity favored over MOND at 5.2 sigma in 23 dwarf spheroidals.
- **Status**: [SUPPORTS]

### H4. Lelli-McGaugh-Schombert 2016 (SPARC)
- **Authors**: F. Lelli, S. S. McGaugh, and J. M. Schombert
- **Title**: "SPARC: Mass Models for 175 Disk Galaxies with Spitzer Photometry and Accurate Rotation Curves"
- **Journal**: AJ **152**, 157 (2016)
- **arXiv**: arXiv:1606.09251
- **Cited in**: [P3]
- **Key result**: SPARC database: 175 galaxies with rotation curves and mass models; THE benchmark dataset.
- **Status**: [FOUNDATIONAL]

### H5. McGaugh-Lelli-Schombert 2016 (RAR)
- **Authors**: S. S. McGaugh, F. Lelli, and J. M. Schombert
- **Title**: "Radial Acceleration Relation in Rotationally Supported Galaxies"
- **Journal**: Phys. Rev. Lett. **117**, 201101 (2016)
- **arXiv**: arXiv:1609.05917
- **Cited in**: [P3]
- **Key result**: Tight correlation g_obs = f(g_bar) with a_0 = 1.2 x 10^-10 m/s^2; any theory must reproduce this.
- **Status**: [CONSTRAINS]

### H6. Haubner-Lelli et al. 2024 (BIG-SPARC)
- **Authors**: T. Haubner, F. Lelli et al.
- **Title**: (BIG-SPARC: ~4000 galaxies)
- **Journal**: (2024)
- **arXiv**: arXiv:2411.13329
- **Cited in**: [P3]
- **Key result**: ~4000 galaxy dataset, 20x larger than original SPARC; future benchmark for Paper 3.
- **Status**: [CONSTRAINS]

### H7. Milgrom 1983 (MOND)
- **Authors**: M. Milgrom
- **Title**: "A modification of the Newtonian dynamics as a possible alternative to the hidden mass hypothesis"
- **Journal**: ApJ **270**, 365 (1983)
- **arXiv**: --
- **Cited in**: [P3]
- **Key result**: Original MOND proposal; a_0 ~ cH_0 connection hints at cosmological origin.
- **Status**: [FOUNDATIONAL]

### H8. Donoghue 1994
- **Authors**: J. F. Donoghue
- **Title**: "General relativity as an effective field theory: The leading quantum corrections"
- **Journal**: Phys. Rev. D **50**, 3874 (1994)
- **arXiv**: --
- **Cited in**: [P3] (implicit)
- **Key result**: GR as EFT; quantum corrections are computable at low energies.
- **Status**: [FOUNDATIONAL]

### H9. Casini-Huerta 2017
- **Authors**: H. Casini and M. Huerta
- **Title**: "Relative entropy and the RG flow"
- **Journal**: J. High Energy Phys. (2017)
- **arXiv**: arXiv:1611.00016
- **Cited in**: [P3]
- **Key result**: QRE monotonicity under RG flow = c-theorem; RG flow IS data processing (DPI). THE mathematical bridge for Paper 3.
- **Status**: [FOUNDATIONAL]

### H10. Smolin 2017
- **Authors**: L. Smolin
- **Title**: "MOND as a regime of quantum gravity"
- **Journal**: Phys. Rev. D **96**, 083523 (2017)
- **arXiv**: arXiv:1704.00780
- **Cited in**: [P3]
- **Key result**: MOND from quantum gravity in de Sitter; a_0 ~ cH_0 emerges naturally.
- **Status**: [SUPPORTS]

### H11. Brouwer et al. 2021
- **Authors**: M. Brouwer et al.
- **Title**: (KiDS-1000 weak lensing test of emergent gravity)
- **Journal**: A&A **650** (2021)
- **arXiv**: arXiv:2106.11677
- **Cited in**: [P3]
- **Key result**: KiDS-1000 weak lensing extends RAR by 2 decades; supports Verlinde.
- **Status**: [SUPPORTS]

### H12. Rostami-Rezazadeh-Rostampour 2025
- **Authors**: T. Rostami, K. Rezazadeh, and A. Rostampour
- **Title**: "Relativistic MOND from Modified Entropic Gravity"
- **Journal**: (2025)
- **arXiv**: arXiv:2511.05632
- **Cited in**: [P3]
- **Key result**: Temperature-dependent corrections to equipartition -> MOND interpolation; uses Unruh temperature.
- **Status**: [SUPPORTS]

### H13. Tsallis-Renyi MOND 2025
- **Authors**: (Various)
- **Title**: (Modified Renyi entropy from Tsallis entropy -> MOND-like force law)
- **Journal**: (2025)
- **arXiv**: arXiv:2505.03061
- **Cited in**: [P3]
- **Key result**: Modified entropy -> MOND; explicitly uses entropy (not energy) to derive force law.
- **Status**: [SUPPORTS]

### H14. Milgrom 2020 (FUNDAMOND)
- **Authors**: M. Milgrom
- **Title**: (a_0 ~ cH_0 connection)
- **Journal**: (2020)
- **arXiv**: arXiv:2001.09729
- **Cited in**: [P3]
- **Key result**: a_0 ~ cH_0 ~ c^2 Lambda^{1/2}; local-cosmological connection.
- **Status**: [SUPPORTS]

---

## I. Algebraic QFT & Modular Theory

### I1. Buchholz 1982/1986
- **Authors**: D. Buchholz
- **Title**: "Gauss' law and the infraparticle problem"
- **Journal**: Phys. Lett. B **174**, 331 (1986)
- **arXiv**: --
- **Cited in**: [P1b]
- **Key result**: Infrared superselection sectors; soft graviton emission changes sector = algebraic barrier to recovery.
- **Status**: [FOUNDATIONAL]

### I2. Witten 2022 (Crossed Product)
- **Authors**: E. Witten
- **Title**: "Gravity and the Crossed Product"
- **Journal**: J. High Energy Phys. **2022**(10), 008 (2022)
- **arXiv**: arXiv:2112.12828
- **Cited in**: [P2]
- **Key result**: Gravity promotes Type III -> Type II algebras via crossed product; the algebra gets BETTER with gravity.
- **Status**: [FOUNDATIONAL]

### I3. Faulkner et al. 2022
- **Authors**: T. Faulkner et al.
- **Title**: (Petz map in Type III algebras)
- **Journal**: (2022)
- **arXiv**: --
- **Cited in**: [P1b] (implicit)
- **Key result**: Petz recovery map extends to Type III von Neumann algebras.
- **Status**: [FOUNDATIONAL]

### I4. Chandrasekaran-Longo-Penington-Witten 2023
- **Authors**: V. Chandrasekaran, R. Longo, G. Penington, and E. Witten
- **Title**: "An Algebra of Observables for de Sitter Space"
- **Journal**: J. High Energy Phys. **2023**(02), 082 (2023)
- **arXiv**: arXiv:2206.10780
- **Cited in**: [P2]
- **Key result**: Observer-dependent algebras in de Sitter; Type II structure from observer's horizon.
- **Status**: [SUPPORTS]

### I5. Araki-Uhlmann (Relative Entropy)
- **Authors**: H. Araki (1976); A. Uhlmann (1977)
- **Title**: (Araki relative entropy / Uhlmann fidelity)
- **Journal**: Various
- **arXiv**: --
- **Cited in**: [P2] [P4] (implicit)
- **Key result**: Rigorous definition of relative entropy for von Neumann algebras; foundation for Dorau-Much.
- **Status**: [FOUNDATIONAL]

---

## J. Observations & Astrophysical Tests

### J1. EHT 2019
- **Authors**: Event Horizon Telescope Collaboration
- **Title**: "First M87 Event Horizon Telescope Results. I."
- **Journal**: Astrophys. J. Lett. **875**, L1 (2019)
- **arXiv**: --
- **Cited in**: [P2]
- **Key result**: Black hole shadow image; Paper 2 predicts 4.6% larger shadow than GR.
- **Status**: [CONSTRAINS]

### J2. NICER 2019
- **Authors**: M. C. Miller et al.
- **Title**: "PSR J0030+0451 Mass and Radius from NICER Data"
- **Journal**: Astrophys. J. Lett. **887**, L24 (2019)
- **arXiv**: --
- **Cited in**: [P2]
- **Key result**: Neutron star mass/radius measurements; Paper 2 predicts 19% different redshift from GR.
- **Status**: [CONSTRAINS]

### J3. LIGO/Virgo/KAGRA 2021
- **Authors**: R. Abbott et al. (LVK)
- **Title**: "Tests of General Relativity with Binary Black Holes from the second LIGO-Virgo Gravitational-Wave Transient Catalog"
- **Journal**: Phys. Rev. D **103**, 122002 (2021)
- **arXiv**: --
- **Cited in**: [P2]
- **Key result**: GR tests with gravitational waves; constraints on deviations.
- **Status**: [CONSTRAINS]

### J4. LVK GW250114
- **Authors**: R. Abbott et al. (LVK)
- **Title**: "GW250114: Observation of Gravitational Waves from a Binary Black Hole Merger with Total Mass ~120 M_sun"
- **Journal**: (2025)
- **arXiv**: --
- **Cited in**: [P2]
- **Key result**: Massive BBH merger; QNM frequency test for Paper 2 predictions.
- **Status**: [CONSTRAINS]

---

## K. Cosmology & Thermodynamic Gravity

### K1. Padmanabhan 2012 (Holographic Equipartition)
- **Authors**: T. Padmanabhan
- **Title**: "Emergent Gravity Paradigm: Recent Progress"
- **Journal**: (2012)
- **arXiv**: arXiv:1206.4916
- **Cited in**: [P4]
- **Key result**: dV/dt = L_P^2 (N_sur - N_bulk); cosmic expansion from information imbalance.
- **Status**: [SUPPORTS]

### K2. Padmanabhan 2017 (CosMIn)
- **Authors**: T. Padmanabhan
- **Title**: "Cosmic information, the cosmological constant and the amplitude of primordial perturbations"
- **Journal**: Phys. Lett. B **773**, 81-85 (2017)
- **arXiv**: arXiv:1703.06144
- **Cited in**: [P4]
- **Key result**: Finite CosMIn REQUIRES late-time accelerated expansion (Lambda > 0).
- **Status**: [SUPPORTS]

### K3. Gibbons-Hawking 1977
- **Authors**: G. W. Gibbons and S. W. Hawking
- **Title**: "Cosmological event horizons, thermodynamics, and particle creation"
- **Journal**: Phys. Rev. D **15**, 2738 (1977)
- **arXiv**: --
- **Cited in**: [P4] (implicit)
- **Key result**: de Sitter space has temperature T_GH = H/(2 pi) and entropy S_GH = A/(4G).
- **Status**: [FOUNDATIONAL]

### K4. Capozziello et al. 2024
- **Authors**: S. Capozziello et al.
- **Title**: (Bogoliubov transformation as quantum channel in gravity)
- **Journal**: Eur. Phys. J. C (2024)
- **arXiv**: arXiv:2406.19274
- **Cited in**: [P4]
- **Key result**: Cosmological particle production modeled as Gaussian bosonic channel.
- **Status**: [SUPPORTS]

### K5. Barman et al. 2026
- **Authors**: S. Barman et al.
- **Title**: (Teleportation fidelity in FRW spacetime)
- **Journal**: QIP (2026)
- **arXiv**: arXiv:2601.20860
- **Cited in**: [P4]
- **Key result**: Expansion degrades quantum fidelity through mode mixing; directly quantifiable via Bogoliubov coefficients.
- **Status**: [SUPPORTS]

---

## L. Our Papers (Self-citations)

### L1. Paper 1: Huang 2026
- **Authors**: S.-K. Huang
- **Title**: "The Arrow of Time as Petz Recovery Failure"
- **Journal**: Submitted to PRL; Zenodo DOI: 10.5281/zenodo.18897853
- **arXiv**: awaiting endorsement
- **Cited in**: [P1b] [P2] [P3] [P4]
- **Key result**: tau = 1 - F; equivalence chain; F >= exp(-Sigma/2); five theorems connecting four domains.
- **Status**: [OUR WORK]

### L2. Paper 1b: Huang 2026
- **Authors**: S.-K. Huang
- **Title**: "Wavefunction Collapse as Retrodiction Failure: Gravitational Decoherence, Objectivity, and the tau-R Bridge"
- **Journal**: In preparation (2026)
- **arXiv**: --
- **Cited in**: [P2]
- **Key result**: Operational collapse = tau -> 1; tau-R bridge; gravitational dephasing analysis.
- **Status**: [OUR WORK]

### L3. Paper 2: Huang 2026
- **Authors**: S.-K. Huang
- **Title**: "Information-Theoretic Origin of the Exponential Metric: Gravitational Redshift as Petz Recovery Fidelity"
- **Journal**: In preparation (2026)
- **arXiv**: --
- **Cited in**: [P1b] [P3] [P4]
- **Key result**: Sigma_grav = r_s/r; g_00 = -exp(-r_s/r); no event horizons; 5 testable predictions.
- **Status**: [OUR WORK]

### L4. Paper 3: Huang 2026
- **Authors**: S.-K. Huang
- **Title**: "Galactic Dynamics from Quantum Information Recovery" (working title)
- **Journal**: In preparation (2026)
- **arXiv**: --
- **Cited in**: [P4]
- **Key result**: Same Sigma with quantum-corrected channel gives flat rotation curves.
- **Status**: [OUR WORK]

### L5. Paper 4: Huang 2026
- **Authors**: S.-K. Huang
- **Title**: "Gravity as Quantum Information Recovery: A Unified Framework" (working title)
- **Journal**: Future
- **arXiv**: --
- **Cited in**: --
- **Key result**: Sigma = D(rho_spacetime || rho_matter); one equation, three regimes.
- **Status**: [OUR WORK]

---

## M. Additional References (Supplementary)

### M1. Manna-Rahman-Chowdhury 2023
- **Authors**: T. Manna, F. Rahman, and T. Chowdhury
- **Title**: "Strong lensing in the exponential wormhole spacetimes"
- **Journal**: New Astron. Rev. (2023)
- **arXiv**: --
- **Cited in**: [P2]
- **Key result**: Strong lensing predictions for exponential metric.
- **Status**: [SUPPORTS]

### M2. Graber 2017 / Makukov-Mychelkin 2023 (Rotating exponential)
- **Authors**: J. R. Graber (2017); M. A. Makukov and E. G. Mychelkin (2023)
- **Title**: "Rotating Yilmaz metric" / "Rotating exponential metric via Newman-Janis algorithm"
- **Journal**: arXiv:1708.05645 (2017); arXiv:2301.08118 (2023)
- **arXiv**: arXiv:1708.05645; arXiv:2301.08118
- **Cited in**: [P2]
- **Key result**: Kerr analog of exponential metric via Newman-Janis; extends Paper 2 to rotating objects.
- **Status**: [SUPPORTS]

### M3. Ernazarov 2025 / Ernazarov-Ivashchuk 2026
- **Authors**: K. K. Ernazarov (2025); K. K. Ernazarov and V. D. Ivashchuk (2026)
- **Title**: "The Yilmaz-Rosen and JNW metric solutions" / "Phantom scalar field solution in the exponential metric"
- **Journal**: arXiv:2510.21625 (2025); (2026)
- **arXiv**: arXiv:2510.21625
- **Cited in**: [P2]
- **Key result**: Exponential metric = phantom scalar field phi = M/r; connects Sigma_grav to scalar field.
- **Status**: [SUPPORTS]

### M4. Guo-Yuan 2025
- **Authors**: Y. Guo and J. Yuan
- **Title**: "Quantum effective dynamics of Papapetrou spacetime"
- **Journal**: Eur. Phys. J. C (2025)
- **arXiv**: --
- **Cited in**: [P2]
- **Key result**: Quantum effective dynamics of the Papapetrou (exponential) spacetime.
- **Status**: [SUPPORTS]

### M5. Ibison 2007
- **Authors**: M. Ibison
- **Title**: "Cosmological test of the Yilmaz theory of gravity"
- **Journal**: arXiv:0705.0080 (2007)
- **arXiv**: arXiv:0705.0080
- **Cited in**: [P2]
- **Key result**: Cosmological test of Yilmaz theory; constraints on exponential metric cosmology.
- **Status**: [CONSTRAINS]

### M6. Nath-Sarma 2024 (QNM)
- **Authors**: S. Nath and D. Sarma
- **Title**: "Quasinormal modes of the generalized exponential wormhole spacetime"
- **Journal**: arXiv:2401.03738 (2024)
- **arXiv**: arXiv:2401.03738
- **Cited in**: [P2]
- **Key result**: QNM frequencies for exponential metric; directly relevant to GW predictions.
- **Status**: [SUPPORTS]

### M7. Ahmadi-Fuentes 2014
- **Authors**: M. Ahmadi, D. E. Bruschi, C. Sabin, G. Adesso, and I. Fuentes
- **Title**: "Relativistic Quantum Metrology"
- **Journal**: Sci. Rep. **4**, 4996 (2014)
- **arXiv**: --
- **Cited in**: [P2]
- **Key result**: Quantum metrology in curved spacetime; channel capacity depends on geometry.
- **Status**: [SUPPORTS]

### M8. Donoghue 2020
- **Authors**: J. F. Donoghue
- **Title**: "A critique of the asymptotic safety program"
- **Journal**: Frontiers in Physics (2020)
- **arXiv**: --
- **Cited in**: [P3] (implicit)
- **Key result**: Critique: running of G in AS may not be physically observable; must be addressed in Paper 3.
- **Status**: [POTENTIAL CONFLICT]

### M9. LVK 2023 (Echoes)
- **Authors**: R. Abbott et al. (LVK)
- **Title**: "Search for gravitational-wave echoes using the third observing run"
- **Journal**: (2023)
- **arXiv**: --
- **Cited in**: [P2]
- **Key result**: No echoes detected yet; Paper 2 predicts echoes at ~2.5 ms.
- **Status**: [CONSTRAINS]

### M10. Hossenfelder 2017
- **Authors**: S. Hossenfelder
- **Title**: (Covariant version of emergent gravity reduces to GR)
- **Journal**: Phys. Rev. D **95** (2017)
- **arXiv**: arXiv:1703.01415
- **Cited in**: [P3]
- **Key result**: Critique of Verlinde: covariant formulation reduces to GR; theoretical consistency issue.
- **Status**: [POTENTIAL CONFLICT]

### M11. Lelli-McGaugh-Schombert 2017 (Tension)
- **Authors**: F. Lelli, S. S. McGaugh, and J. M. Schombert
- **Title**: (Verlinde requires decreased M/L ratios)
- **Journal**: ApJ **836** (2017)
- **arXiv**: arXiv:1702.04355
- **Cited in**: [P3]
- **Key result**: Verlinde requires lower M/L ratios, in tension with stellar population models.
- **Status**: [POTENTIAL CONFLICT]

---

## Summary Statistics

### Total References: 120

### References by Paper

| Paper | Cited References | Unique to Paper |
|-------|-----------------|-----------------|
| P1    | 33              | 15              |
| P1b   | 32              | 14              |
| P2    | 48              | 25              |
| P3    | 22              | 10              |
| P4    | 15              | 5               |

### References Shared Across Papers

| Reference | Papers Citing |
|-----------|--------------|
| Petz 1986 (A1) | P1, P1b, P2 |
| Petz 1988 (A2) | P1, P1b, P2 |
| Fawzi-Renner 2015 (A3) | P1, P1b, P2 |
| JRSWW 2018 (A4) | P1, P1b, P2 |
| Parzygnat-Buscemi 2023 (A5) | P1, P1b |
| Crooks 1999 (C1) | P1, P2 |
| Zurek 2003 (D3) | P1, P1b |
| Pikovski 2015 (E3) | P1b, P2 |
| Jacobson 1995 (G4) | P2, P3, P4 |
| Dorau-Much 2025 (G14) | P2, P3, P4 |
| Bianconi 2025 (G15) | P2, P4 |
| Verlinde 2016 (G7) | P3, P4 |
| Paper 1 (L1) | P1b, P2, P3, P4 |

### Most-Cited References (across all papers)

1. **Petz 1986/1988** (A1, A2) -- cited in P1, P1b, P2 -- THE mathematical foundation
2. **Fawzi-Renner 2015** (A3) -- cited in P1, P1b, P2 -- CMI bound
3. **JRSWW 2018** (A4) -- cited in P1, P1b, P2 -- universal recovery bound
4. **Jacobson 1995** (G4) -- cited in P2, P3, P4 -- gravity from thermodynamics
5. **Dorau-Much 2025** (G14) -- cited in P2, P3, P4 -- QRE to Einstein equations
6. **Paper 1** (L1) -- cited by all subsequent papers

### References by Status

| Status | Count |
|--------|-------|
| [FOUNDATIONAL] | 52 |
| [SUPPORTS] | 43 |
| [CONSTRAINS] | 13 |
| [INDEPENDENT VERIFICATION] | 1 |
| [POTENTIAL CONFLICT] | 5 |
| [OUR WORK] | 5 |

### Potential Conflicts to Address

1. **Misner 1995** (G24): "Yilmaz Cancels Newton" -- critique of Yilmaz theory; Paper 2 must distinguish our derivation from Yilmaz's
2. **Bonder et al. 2016** (E17): Critique of Pikovski mechanism -- must acknowledge limitations in Paper 1b
3. **Donoghue 2020** (M8): Running G in asymptotic safety may not be physically observable -- Paper 3 must argue our approach differs
4. **Hossenfelder 2017** (M10): Covariant Verlinde reduces to GR -- Paper 3 needs different covariant formulation
5. **Lelli et al. 2017** (M11): Verlinde requires problematic M/L ratios -- Paper 3 should test against this

### Key Missing References (to be located)

1. Page 1993 (information paradox) -- referenced in MEMORY but not found in current papers
2. Hayden-Preskill 2007 -- referenced in MEMORY but not in current bibliographies
3. Maldacena-Susskind 2013 (ER=EPR) -- referenced in MEMORY but not in current bibliographies
4. Shor 1995 -- referenced in task but not in current papers
5. Kasevich 2022 (gravitational AB) -- referenced in task but not in current papers
6. LISA Pathfinder -- referenced in task but not in current papers
7. Haroche cat states -- referenced in task but not currently cited
8. Horchani 2025 (distinguishing test) -- referenced in task but not in current papers

---

*This master reference database was compiled from all .tex files in `/Users/akaihuangm1/Desktop/github/petz-recovery-unification/paper/`, all .md files in `/Users/akaihuangm1/Desktop/github/petz-recovery-unification/research/`, and Paper 2 draft in `/Users/akaihuangm1/Desktop/github/0218/quantum-retrocausality-ai/papers/paper2_exponential_metric.tex`.*
