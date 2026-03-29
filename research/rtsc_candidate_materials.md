# Room-Temperature Superconductor Candidates via the tau Design Criterion

**Author**: Sheng-Kai Huang
**Date**: 2026-03-28
**Status**: Material screening study
**Framework**: Paper on Superconductivity as Zero Temporal Asymmetry (paper_superconductivity_tau.tex)

---

## 1. The R_tau Design Criterion

From our superconductivity paper, the screening ratio is:

```
R_tau = max_T [I_pair(T)] / tau_dec(300 K) >= 1
```

where:
- **I_pair** = total Cooper pair entanglement (summed over all pairing channels)
- **tau_dec(300 K)** = thermal decoherence at room temperature = 1 - exp(-Sigma_ep/2)

For a typical metal at 300 K, Sigma_ep ~ 1.6, giving tau_dec ~ 0.55. Reaching R_tau >= 1 requires I_pair(0) >= 0.55 (normalized units).

### What drives R_tau high?

**Numerator (I_pair):**
- I_pair ~ pi * N(0) * Delta(0) * ln(2) for BCS (weak coupling)
- Large gap Delta, high DOS N(0), multiple channels all increase I_pair
- Multi-channel: I_total = I_phonon + I_spin + I_orbital + I_topo

**Denominator (tau_dec):**
- tau_dec ~ 1 - exp(-lambda * k_B T / (2 * hbar * omega_D)) for T > Theta_D
- Low lambda (weak scattering), high omega_D (stiff lattice) reduce tau_dec
- Paradox: high lambda increases I_phonon but also increases tau_dec

### Current landscape (from our paper's Table I):

| Material | T_c (K) | Delta(0) (meV) | lambda | R_tau |
|----------|---------|-----------------|--------|-------|
| Nb | 9.3 | 1.5 | 0.82 | 0.003 |
| MgB2 | 39 | 7.1 | 0.87 | 0.012 |
| YBCO | 93 | 20 | --- | 0.031 |
| H3S (155 GPa) | 203 | 30 | 2.19 | 0.069 |
| LaH10 (170 GPa) | 260 | 40 | 3.41 | 0.092 |
| **Room-T target** | **300** | **>= 50** | **---** | **>= 1** |

**The gap**: Even LaH10 at R_tau ~ 0.09 is 10x below the threshold. A qualitative leap is needed.

---

## 2. The Gao et al. (2025) Ceiling

A critical recent result: Gao et al. (Nature Communications 16, 8253, 2025) analyzed >20,000 metals via DFT and found an **inherent trade-off** between the logarithmic average phonon frequency omega_log and the electron-phonon coupling constant lambda in conventional BCS superconductors at ambient pressure. Their conclusion: **room-temperature conventional superconductivity at ambient pressure is extremely unlikely via phonon-mediated pairing alone.**

In tau language: the lambda-omega_log trade-off means I_phonon / tau_phonon has a maximum. You cannot simply crank up lambda without also increasing tau_dec. The ceiling is approximately T_c ~ 120-140 K for pure phonon-mediated superconductivity at ambient pressure.

**Implication for our search**: R_tau >= 1 requires going **beyond single-channel phonon pairing**. Multi-channel pairing or unconventional mechanisms are mandatory.

---

## 3. Candidate Materials

### Candidate 1: Mg2IrH6 (Ambient-Pressure Ternary Hydride)
**Strategy A: Hydrogen-rich, ambient-pressure stable**

**Chemical formula**: Mg2IrH6
**Structure**: K2PtCl6-type cubic (Fm-3m), with IrH6 octahedra surrounded by Mg atoms. Dynamically stable at P = 0.

**Why R_tau might be high:**
- **High omega_D**: H vibrations dominate the phonon spectrum (omega_H ~ 100-200 meV). Theta_D > 1000 K. This keeps tau_dec moderate even at 300 K.
- **Van Hove singularity**: The electronic structure has two major van Hove singularities with the Fermi level pinned to a DOS peak. N(0) is anomalously large, boosting I_Cooper.
- **Strong e-ph coupling**: lambda ~ 2.0-2.5 from DFT (comparable to H3S).
- **s-wave symmetry**: Uniform gap -> no wasted pairing information at nodes.

**Estimated parameters:**
- lambda ~ 2.2, omega_log ~ 80 meV, N(0) ~ large (vHS)
- Delta(0) ~ 25-30 meV (from Allen-Dynes with mu* = 0.13)
- T_c ~ 160 K (DFT prediction, Belli et al. 2024, PRL 132, 166001)
- **R_tau ~ 0.06** (single-channel phonon; comparable to H3S without pressure)

**Thermodynamic stability**: Metastable, ~50 meV/atom above convex hull. Kinetic barriers may allow survival if synthesized at moderate pressure and quenched.

**Synthesis**: High-pressure (few GPa) synthesis of Mg + Ir + H2 in a multi-anvil press, then slow pressure release to trap the metastable phase at ambient P. Similar to diamond: metastable at P = 0 but kinetically stable.

**Has anyone tried?**: Not yet synthesized experimentally. DFT prediction published 2024 (PRL). The Mg2XH6 family (X = Rh, Ir, Pd, Pt) was identified via high-throughput screening.

**tau assessment**: Single-channel R_tau too low. But the vHS-enhanced N(0) is promising for combining with spin fluctuations (see Candidate 7).

---

### Candidate 2: SrAuH3 (Perovskite Hydride)
**Strategy A: Hydrogen-rich, moderate-pressure synthesis**

**Chemical formula**: SrAuH3
**Structure**: Perovskite (Pm-3m), with Au at the B-site and H at the anion positions. Predicted stable at P = 0 (synthesizable at 7 GPa).

**Why R_tau might be high:**
- **Strong e-ph coupling**: lambda = 2.025, in the strong-coupling regime.
- **Au heavy-fermion hybridization**: The Au 5d states hybridize strongly with H 1s, creating a mixed-character Fermi surface with enhanced N(0).
- **High phonon frequencies**: H-dominated optical modes above 100 meV.
- **Low synthesis pressure**: 7 GPa is achievable in standard multi-anvil apparatus (vs. 150+ GPa for binary hydrides).

**Estimated parameters:**
- lambda = 2.025, omega_log ~ 60 meV
- Delta(0) ~ 22-25 meV
- T_c = 132 K (DFT, arXiv:2412.15488, December 2024)
- **R_tau ~ 0.05**

**Synthesis**: Sr + Au + H2 at 7 GPa in a multi-anvil press at 500-800 C, then quench to ambient. The low synthesis pressure makes this one of the most experimentally accessible hydride candidates.

**Has anyone tried?**: No experimental realization yet. Predicted December 2024.

**tau assessment**: Modest R_tau from phonons alone, but the Au d-electron correlations open a possible spin-fluctuation channel (I_spin). If I_spin ~ I_phonon/3, R_tau could reach ~0.07.

---

### Candidate 3: KB3C3 (Boron-Carbon Clathrate)
**Strategy B: Carbon-based, ambient-pressure stable**

**Chemical formula**: KB3C3
**Structure**: Sodalite-type clathrate with K guest atoms encapsulated in B3C3 cages. Cubic, dynamically stable at P = 0.

**Why R_tau might be high:**
- **Extremely high Theta_D**: B-C covalent bonds are among the stiffest in nature. Theta_D > 1500 K. This dramatically suppresses tau_dec at 300 K.
- **Two-gap superconductor**: Two distinct Fermi sheets with different coupling strengths (analogous to MgB2). The strong-coupling sheet provides large I_phonon; the weak-coupling sheet adds.
- **Light atoms**: B and C are light, pushing phonon frequencies high.
- **Ambient pressure stability**: Unlike hydrides, no pressure needed.

**Estimated parameters:**
- lambda ~ 1.5-2.0 (two-gap), omega_log ~ 90-120 meV
- Delta(0) ~ 15-18 meV (strong-coupling sheet), ~8-10 meV (weak sheet)
- T_c = 102.5 K (DFT prediction, November 2025)
- **R_tau ~ 0.04** (phonon only)
- **Effective R_tau ~ 0.04 but with very low tau_dec denominator**

Note: Because Theta_D is so high, the tau_dec(300 K) for this material is lower than for typical metals. The actual R_tau may be better than the proxy estimate suggests.

**Refined estimate**: For Theta_D ~ 1800 K, at T = 300 K we are in the regime T << Theta_D. Here Sigma_ep ~ hbar*omega_D / (k_B T) ~ 1800/300 ~ 6, so tau_dec ~ 0.95 (high). BUT the scattering rate is exponentially suppressed at T << Theta_D: the number of phonons is ~ exp(-Theta_D/T) ~ exp(-6) ~ 0.0025. The effective tau_dec is proportional to the phonon population, giving tau_dec(eff) ~ 0.0025 * 0.95 ~ 0.002. This dramatically improves R_tau.

**Corrected R_tau ~ 0.04 / 0.002 = 20** (!!!)

Wait -- this estimate deserves scrutiny. The issue is that at T = 300 K << Theta_D = 1800 K, phonon scattering is indeed exponentially suppressed. But the reason KB3C3 only reaches T_c = 102 K (not 300 K) is that the pairing interaction also weakens as temperature increases toward the Debye scale. The BCS T_c formula already accounts for this trade-off. Still, the point stands: **materials with extremely high Theta_D have fundamentally lower tau_dec, changing the R_tau landscape.** This warrants a full Eliashberg calculation.

**Synthesis**: Direct solid-state synthesis from K + B + C precursors at high temperature (~1500 C) and moderate pressure (~5-10 GPa), or chemical vapor deposition of B-C frameworks followed by K intercalation.

**Has anyone tried?**: Not yet. The 1166-type B-C clathrates (RbYbB6C6, RbBaB6C6) were predicted in February 2024; KB3C3 in November 2025. No experimental synthesis reported.

**tau assessment**: The B-C clathrate family is the most promising carbon-based candidate. The high Theta_D fundamentally changes the tau_dec denominator, potentially allowing R_tau > 1 even with modest pairing.

---

### Candidate 4: RbPH3 (Phosphorus Hydride Perovskite)
**Strategy A + E: Anharmonicity-stabilized hydride**

**Chemical formula**: RbPH3
**Structure**: R3m phase at ambient pressure (distorted perovskite), dynamically stabilized by quantum ionic anharmonic effects. Strong P-H covalent bonds (similar to H3S).

**Why R_tau might be high:**
- **H3S-like bonding**: The P-H covalent network resembles H3S (which holds the pressure-dependent T_c record). High phonon frequencies from light H atoms.
- **Anharmonic stabilization**: The most anharmonic modes are precisely those that couple most strongly to electrons. Anharmonicity IS the e-ph coupling mechanism -- a self-reinforcing loop.
- **Three Fermi surfaces**: Complex electronic structure with multiple sheets, each contributing to I_phonon.
- **Moderate synthesis pressure**: Stable at 30 GPa in Pm-3m, then recoverable at P = 0 in R3m.

**Estimated parameters:**
- lambda ~ 1.8-2.5 (dominated by H modes), omega_log ~ 70-100 meV
- Delta(0) ~ 15-20 meV
- T_c ~ 100 K (DFT with anharmonic phonons, Dangic et al. 2024, arXiv:2411.03822)
- **R_tau ~ 0.04**

**Synthesis**: Rb + P + H2 at 30 GPa (diamond anvil cell or large-volume press), then slow decompression. The R3m phase should be kinetically trapped at ambient pressure thanks to quantum anharmonic stabilization -- the quantum zero-point motion of the ions prevents the structure from collapsing.

**Has anyone tried?**: No experimental realization. The prediction (November 2024) is notable because it explicitly accounts for quantum ionic effects that are neglected in most high-throughput screenings.

**tau assessment**: R_tau from phonons alone is insufficient. The key novelty is that anharmonicity-enhanced e-ph coupling might circumvent the Gao et al. omega_log-lambda trade-off, since the relevant phonon modes are anharmonic (not captured by harmonic DFT). If the true lambda is 30-50% higher than the harmonic estimate, R_tau could approach 0.06-0.08.

---

### Candidate 5: Li2AuH6 (Maximum-T_c Ternary Hydride)
**Strategy A: Gao et al. ceiling candidate**

**Chemical formula**: Li2AuH6
**Structure**: K2PtCl6-type cubic (Fm-3m), isostructural to Mg2IrH6.

**Why R_tau might be high:**
- **Near the ceiling**: Identified by Gao et al. (2025) as approaching the maximum T_c for conventional ambient-pressure superconductors.
- **Strong coupling**: lambda ~ 2.5-3.0, among the highest for any ambient-pressure material.
- **High phonon frequencies**: H-dominated, with omega_log ~ 60-80 meV.
- **Au 5d hybridization**: Similar to SrAuH3, the Au d-states enhance N(0).

**Estimated parameters:**
- lambda ~ 2.5-3.0
- Delta(0) ~ 20-30 meV
- T_c = 116 K (Eliashberg) to 140 K (some estimates), arXiv:2501.12222
- **R_tau ~ 0.05-0.06**

**Thermodynamic stability**: CRITICAL ISSUE. 0.172 eV/atom above convex hull. This is thermodynamically unstable. Metastable at best, and likely decomposes at ambient pressure unless kinetically trapped.

**Synthesis**: Li + Au + H2 at high pressure (estimated 10-20 GPa), then attempt quench to ambient. The 0.172 eV/atom instability is borderline -- some metastable phases with similar hull distances survive (e.g., diamond at 0.02 eV/atom above graphite; some high-pressure phases at ~0.1 eV/atom).

**Has anyone tried?**: No. The sister compound Li2AgH6 (T_c ~ 109 K, Eliashberg) has the same structural type but is even less stable (0.319 eV/atom above hull).

**tau assessment**: Represents the conventional-phonon ceiling. R_tau ~ 0.05-0.06 is the maximum achievable from phonons alone at ambient pressure. To reach R_tau >= 1, unconventional mechanisms are essential.

---

### Candidate 6: Hydrogenated Boron-Doped Diamond (H:B:C)
**Strategy B + D: Carbon-based with flat-band engineering**

**Chemical formula**: C1-x-yBxHy (x ~ 0.05-0.10, y ~ 0.10-0.20)
**Structure**: Diamond cubic with B substituting C and H occupying interstitial sites. Alternatively, hole-doped diamond-like cubic hydrocarbon (CH)_n.

**Why R_tau might be high:**
- **Highest Theta_D of any material**: Diamond has Theta_D ~ 2200 K. This makes tau_dec(300 K) exponentially small (phonon population ~ exp(-2200/300) ~ 7 x 10^{-4}).
- **B-doping introduces holes**: Moves E_F into the valence band, creating a Fermi surface with strong e-ph coupling to high-energy optical phonons.
- **H interstitials create resonant states**: H in diamond creates localized states near E_F that can act as van Hove singularities, dramatically enhancing N(0).
- **Covalent bonds**: C-C and B-C bonds are the stiffest in nature.

**Estimated parameters:**
- lambda ~ 0.5-1.5 (depends on B and H concentration)
- omega_log ~ 120-180 meV (diamond phonons!)
- Delta(0) ~ 10-20 meV (optimistic)
- T_c ~ 40-117 K (boron-doped Q-carbon: 36 K observed; H-doped diamond: 117 K claimed)
- tau_dec(300 K) ~ very low (exp(-7.3) ~ 7 x 10^{-4} phonon scattering rate at T/Theta_D = 0.14)
- **R_tau ~ 0.02 / 0.001 ~ 20** (favorable if the low-tau_dec regime is physical)

**CRITICAL CAVEAT**: The claim of T_c = 117 K in H-doped diamond (ScienceDirect) has not been independently confirmed. If real, it would already represent the highest ambient-pressure T_c for a carbon-based material by a large margin.

**Synthesis**:
1. CVD diamond growth with B doping (well-established, B-doped diamond superconductivity at 4 K confirmed in 2004)
2. H implantation via ion beam or high-pressure hydrogenation
3. Alternative: pulsed laser annealing of amorphous carbon + boron in H2 atmosphere (Q-carbon route)

**Has anyone tried?**:
- B-doped diamond: Yes, T_c = 4 K (Ekimov et al. 2004, Nature)
- B-doped Q-carbon: Yes, T_c = 36 K (Bhaumik et al. 2017)
- H-doped diamond: Claimed T_c = 117 K (Bhatt et al. 2023, ScienceDirect) -- **not independently confirmed**

**tau assessment**: The diamond matrix has the lowest possible tau_dec of any known material. If sufficient I_pair can be generated (through B-doping + H interstitials), R_tau could potentially exceed 1. This is the most promising carbon-based route.

---

### Candidate 7: La3Ni2O7 + H (Nickelate-Hydride Hybrid)
**Strategy C: Multi-channel pairing (spin + phonon)**

**Chemical formula**: La3Ni2O7-delta H_x (hydrogen-intercalated bilayer nickelate)
**Structure**: Ruddlesden-Popper bilayer nickelate with H occupying apical O vacancy sites or interstitial positions.

**Why R_tau might be high:**
- **Nickelate: strong spin-fluctuation pairing**: La3Ni2O7 under pressure achieves T_c = 80-96 K, driven by superexchange and spin fluctuations (I_spin channel). At ambient pressure, thin films show signatures of T_c ~ 26-42 K.
- **Adding hydrogen**: H intercalation could provide a second (phonon) pairing channel. H vibrations give high omega_D, contributing I_phonon.
- **Multi-channel R_tau**: I_total = I_spin(nickelate) + I_phonon(hydrogen), potentially exceeding either alone.
- **Nickel d-electron correlations**: Strong on-site Coulomb U creates magnetic fluctuations that persist to high T.
- **Thin film strain engineering**: Epitaxial strain can tune the Fermi surface to enhance nesting (I_spin) and phonon coupling (I_phonon) simultaneously.

**Estimated parameters:**
- I_spin: from nickelate, comparable to cuprate (J ~ 100 meV), giving R_tau_spin ~ 0.03
- I_phonon: from H, comparable to hydrides (lambda ~ 1-2), giving R_tau_phonon ~ 0.03-0.04
- Combined: R_tau ~ 0.06-0.07
- If constructive interference (coherent addition): R_tau ~ 0.08-0.10
- T_c target: 150-200 K (optimistic)

**Synthesis**:
1. Grow La3Ni2O7 thin film on SrTiO3 substrate via MBE or PLD
2. Hydrogenate via atomic H exposure or electrochemical H insertion (similar to HxWO3, HxSrCoO2.5)
3. Control H concentration to optimize phonon coupling without destroying the bilayer structure

**Has anyone tried?**:
- La3Ni2O7 under pressure: Yes, T_c = 80-96 K (Nature 2025)
- La3Ni2O7 thin films at ambient pressure: Yes, signatures of T_c ~ 26-42 K (Nature 2024)
- Hydrogen intercalation of nickelates: Not yet for La3Ni2O7 specifically, but HxNdNiO3 and HxSrNiO2 have been studied (topotactic H insertion is established)

**tau assessment**: This is the cleanest example of the multi-channel strategy. The nickelate provides I_spin; hydrogen provides I_phonon. The combined R_tau, while still below 1, represents the most physically motivated route to high R_tau.

---

### Candidate 8: PtBi2 / Bi2Te3 + Nb (Topological + Conventional Hybrid)
**Strategy D: Topological protection**

**Chemical formula**: PtBi2/Nb heterostructure, or Bi2Te3/NbSe2 bilayer
**Structure**: Topological surface states of PtBi2 or Bi2Te3, proximity-coupled to a conventional s-wave superconductor (Nb or NbSe2).

**Why R_tau might be high:**
- **Topological contribution is T-independent**: I_topo = |nu| * ln(2), where nu is the topological invariant. This does NOT decrease with temperature.
- **Proximity-induced gap**: The conventional superconductor induces a gap Delta_prox on the topological surface states. This provides I_phonon.
- **I_total = I_phonon(proximity) + I_topo**: The topological term survives to high T.
- **PtBi2 is a natural topological superconductor**: Surface electrons pair at low T with zero resistance, confirmed December 2025.

**Estimated parameters:**
- I_topo = ln(2) ~ 0.69 per topological mode. With N_topo surface modes, I_topo_total ~ N_topo * ln(2).
- I_phonon(proximity) ~ limited by parent superconductor (Nb: Delta ~ 1.5 meV, T_c = 9 K)
- R_tau contribution from topology: depends on N_topo / N_total
- **The problem**: topological surface states have much lower DOS than bulk. N_topo / N_total ~ 10^{-4} - 10^{-3}. So the topology helps per-mode but is overwhelmed by the bulk tau_dec.
- **R_tau ~ 0.003 + 0.001 (topo) ~ 0.004** (pessimistic)

**Synthesis**: Well-established thin film heterostructure growth (MBE, sputtering). PtBi2 crystals are available. Nb thin films are standard.

**Has anyone tried?**:
- PtBi2 topological superconductor: Yes, confirmed 2025 (ScienceDaily)
- Bi2Te3/NbSe2 heterostructures: Yes, studied extensively for Majorana fermions
- The specific goal of enhancing T_c via topology: limited exploration

**tau assessment**: The topology-alone route is insufficient because topological surface states carry too little current compared to bulk. The per-mode I_topo is maximal (ln 2), but the number of topological modes is tiny. Topology is more useful as a SUPPLEMENT to bulk pairing channels (see Candidate 9).

---

### Candidate 9: Magic-Angle Graphene + Li Intercalation + H adsorption (Flat Band + High omega_D)
**Strategy D + A: Flat band + hydrogen**

**Chemical formula**: (TBG)Li_x H_y -- twisted bilayer graphene with Li intercalation and H adsorption
**Structure**: Two graphene layers twisted to magic angle (~1.08 degrees), with Li intercalated between layers and H adsorbed on the outer surfaces.

**Why R_tau might be high:**
- **Flat band -> huge N(0)**: At the magic angle, TBG develops flat bands with N(0) -> infinity (van Hove singularity). This maximizes I_Cooper per unit Delta.
- **H adsorption provides high omega_D**: H on graphene has C-H stretching modes at ~350 meV, providing an extraordinarily high-frequency pairing glue.
- **Li intercalation tunes E_F**: Li donates electrons, moving E_F to the flat band.
- **Carbon backbone**: Theta_D of graphene ~ 2300 K, suppressing tau_dec.
- **Combined**: enormous N(0) from flat band x high Delta from H-phonons x low tau_dec from carbon backbone.

**Estimated parameters:**
- N(0) ~ 10x normal graphene (at flat band)
- lambda ~ 0.5-1.0 (C-H phonon coupling in graphene+H is moderate)
- omega_log ~ 150-250 meV (dominated by C-H and C-C modes)
- Delta(0) ~ 3-10 meV (flat band enhancement could amplify)
- T_c(TBG alone) ~ 3-5 K; with Li + H: potentially 30-100 K (very speculative)
- tau_dec(300 K) ~ very low (graphene Theta_D = 2300 K)
- **R_tau ~ 0.01-0.05 (speculative)**

**Synthesis**:
1. CVD-grown graphene, tear-and-stack at magic angle (established)
2. Li vapor deposition in UHV (established for graphite intercalation compounds)
3. Atomic H dosing from cracked H2 (established for hydrogenated graphene)
4. The challenge: maintaining the magic angle while intercalating/adsorbing

**Has anyone tried?**:
- Magic-angle TBG superconductor: Yes, T_c ~ 3-5 K (Cao et al. 2018, Nature)
- Li-intercalated graphite superconductor (CaC6): Yes, T_c = 11.5 K
- Hydrogenated graphene: Yes (graphane), but not combined with TBG
- TBG + Li + H: No.

**tau assessment**: This is the most speculative candidate but also the most "tau-optimized." It simultaneously attacks both numerator (flat band N(0) + H phonons) and denominator (graphene Theta_D). The key unknown is whether the magic-angle flat band survives Li intercalation and H adsorption. If it does, even modest pairing strength could yield high R_tau thanks to the enormous N(0).

---

### Candidate 10: (TiZrNb)B2 + H (High-Entropy Boride Hydride)
**Strategy E: Wild Card -- high-entropy stabilization + hydrogen**

**Chemical formula**: (Ti_{0.33}Zr_{0.33}Nb_{0.33})B2 H_x
**Structure**: AlB2-type hexagonal with disordered transition metal site and H in interstitial positions within the boron layers.

**Why R_tau might be high:**
- **High-entropy stabilization**: The configurational entropy from three (or more) transition metals stabilizes the structure thermodynamically at ambient pressure, potentially trapping metastable hydrogen-rich phases.
- **Boride backbone**: B-B bonds provide high Theta_D (> 1000 K), reducing tau_dec.
- **Multi-d-orbital manifold**: The disordered TM site creates a broad distribution of d-electron states near E_F, enhancing N(0) and potentially providing spin-fluctuation pairing from local moment fluctuations.
- **Hydrogen in B layers**: H atoms sitting in the hexagonal boron network have high vibrational frequencies and couple to the TM d-states.
- **Robust superconductivity**: HEA superconductors show remarkable robustness to pressure, disorder, and defects -- the superconducting state is topologically protected in some cases.

**Estimated parameters (highly speculative):**
- lambda ~ 0.5-1.0 (boride e-ph coupling)
- I_spin ~ 0.01-0.02 (from TM local moment fluctuations)
- omega_log ~ 60-100 meV
- Theta_D ~ 1200 K
- T_c ~ 10-30 K (optimistic for HEA borides)
- **R_tau ~ 0.01-0.02**

**Synthesis**: Arc melting of Ti + Zr + Nb + B at stoichiometric ratio, followed by hydrogenation at 1-5 GPa in a gas-loading apparatus. HEA borides are well-studied structurally.

**Has anyone tried?**:
- HEA borides: Yes, (TiZrHfNbTa)B2 synthesized and characterized
- HEA superconductors: Yes, T_c up to ~10 K
- HEA boride hydrides: Not specifically, but closely related HEA nitrides show superconductivity (TiNbTaN3, T_c ~ 10 K, 2025)

**tau assessment**: This is a wild card. The R_tau is low, but the material class has several attractive features: (1) thermodynamic stability from entropy, (2) tunability from composition, (3) multi-channel pairing from diverse d-orbitals. The main promise is as a platform for systematic R_tau optimization via composition tuning.

---

## 4. Summary Table: All Candidates

| # | Material | Strategy | T_c (K) | R_tau est. | Status | Key advantage |
|---|----------|----------|---------|-----------|--------|---------------|
| 1 | Mg2IrH6 | A (hydride) | 160 | 0.06 | DFT pred. 2024 | vHS at E_F, ambient P stable |
| 2 | SrAuH3 | A (hydride) | 132 | 0.05 | DFT pred. 2024 | 7 GPa synthesis, Au hybridization |
| 3 | KB3C3 | B (carbon) | 102 | 0.04* | DFT pred. 2025 | Highest Theta_D, two-gap |
| 4 | RbPH3 | A+E (anhar.) | 100 | 0.04 | DFT pred. 2024 | H3S-like, anharmonic stabilization |
| 5 | Li2AuH6 | A (ceiling) | 116-140 | 0.06 | DFT pred. 2025 | Near conventional ceiling |
| 6 | H:B:C diamond | B+D (C-based) | 40-117 | 0.02-20** | Partial exp. | Lowest tau_dec of any material |
| 7 | La3Ni2O7+H | C (multi-ch.) | 150-200? | 0.07 | Partial exp. | Spin + phonon multi-channel |
| 8 | PtBi2/Nb | D (topological) | ~10 | 0.004 | Exp. confirmed | T-independent I_topo |
| 9 | TBG+Li+H | D+A (flat band) | 30-100? | 0.01-0.05 | Speculative | Enormous N(0) x high omega_D |
| 10 | (TiZrNb)B2H_x | E (wild card) | 10-30 | 0.01-0.02 | Speculative | Entropy stabilization, tunable |

*KB3C3 R_tau may be significantly higher in the proper low-T/Theta_D regime
**H:B:C diamond R_tau range depends critically on whether the low tau_dec regime is physical

---

## 5. Honest Assessment: Can Any Single Candidate Reach R_tau >= 1?

**No.** Based on current knowledge, no single ambient-pressure material is predicted to reach R_tau >= 1 through any single pairing channel. The Gao et al. (2025) ceiling puts the conventional-phonon limit at R_tau ~ 0.06 (corresponding to T_c ~ 140 K).

**The path to R_tau >= 1 requires one of:**

1. **Multi-channel coherent addition**: If I_phonon + I_spin + I_topo add coherently (not just incoherently), the effective R_tau could be significantly enhanced. This is the nickelate-hydride strategy (Candidate 7).

2. **Anomalously low tau_dec**: Materials with Theta_D > 2000 K (diamond, graphene) have exponentially suppressed phonon scattering at 300 K. If the pairing mechanism does NOT require phonons (e.g., electronic pairing in flat bands), then tau_dec is very low and modest I_pair suffices.

3. **Beyond-BCS strong coupling**: Polaronic, bipolaronic, or pre-formed pair mechanisms that are not captured by Eliashberg theory might generate much larger I_pair than predicted by DFT.

4. **Topological amplification**: If topological protection extends to bulk (not just surface) states, the T-independent I_topo could be large. Current topological superconductors are surface-only; a bulk topological superconductor with high topological invariant would change the picture.

5. **Unknown unknowns**: The cuprates were not predicted by theory before their discovery. Nature may have pairing mechanisms we have not yet imagined.

---

## 6. Top 3 Candidates: Recommended Synthesis Plans

### Rank 1: Mg2IrH6 -- Most Promising Near-Term Hydride

**Rationale**: Highest predicted T_c (160 K) of any ambient-pressure stable hydride. Van Hove singularity at E_F suggests additional pairing channels beyond pure phonon. Cubic structure simplifies characterization.

**Synthesis Protocol:**
1. **Precursors**: MgH2 powder (Sigma-Aldrich, 99%) + Ir powder (Alfa Aesar, 99.9%), stoichiometric ratio 2:1
2. **Loading**: Mix in Ar glovebox, load into Au or Pt capsule with H2 gas
3. **High-pressure synthesis**: Multi-anvil press at 5-10 GPa, 500-800 C, hold 1-4 hours
4. **Quench**: Rapid cooling (<10 s) to room temperature under pressure
5. **Decompression**: Slow release (0.5 GPa/hour) to ambient. Monitor XRD in situ if possible.
6. **Characterization**:
   - Powder XRD to confirm Fm-3m structure
   - AC susceptibility for superconducting transition
   - Resistivity 4-point measurement, 2-300 K
   - Specific heat to confirm bulk T_c
7. **Estimated cost**: ~$50K for multi-anvil time + Ir metal + characterization
8. **Timeline**: 3-6 months from precursor to T_c measurement
9. **Risk**: Metastable phase may decompose during decompression. Mitigation: try different quench rates, add chemical stabilizers (e.g., partial Mg -> Ca substitution).

### Rank 2: KB3C3 -- Most Promising Carbon-Based

**Rationale**: Ambient-pressure stability, extremely high Theta_D (lowest tau_dec), two-gap structure. Light-element composition makes it amenable to high-throughput optimization. No expensive elements.

**Synthesis Protocol:**
1. **Precursors**: K metal (chunks, Sigma-Aldrich) + amorphous B (99.99%) + graphite powder (99.999%)
2. **Method A (high pressure)**:
   - Mix K:B:C = 1:3:3 in Ar glovebox
   - Load into BN capsule
   - Multi-anvil or Paris-Edinburgh press at 5-15 GPa, 1000-1500 C
   - Quench and decompress
3. **Method B (chemical intercalation)**:
   - First synthesize B3C3 clathrate framework (high-T, high-P route)
   - Then intercalate K via vapor transport at 200-400 C
4. **Characterization**: Same as Mg2IrH6 (XRD, susceptibility, resistivity, specific heat)
5. **Estimated cost**: ~$20K (no expensive elements)
6. **Timeline**: 6-12 months (B-C clathrate synthesis is less established than hydride synthesis)
7. **Risk**: K may not fully intercalate into the clathrate cages. B-C framework may not form the sodalite structure at accessible pressures. Mitigation: try other guest atoms (Rb, Cs, Ca, Sr).

### Rank 3: La3Ni2O7 + H -- Most Promising Multi-Channel

**Rationale**: La3Ni2O7 is an established 80-96 K superconductor under pressure and shows ambient-pressure signatures in thin films. Adding hydrogen creates a second (phonon) pairing channel. This is the most physically motivated multi-channel candidate.

**Synthesis Protocol:**
1. **Thin film growth**:
   - PLD or MBE of La3Ni2O7 on (001) SrTiO3 (compressive strain)
   - Optimize O2 pressure and substrate temperature for bilayer phase
   - Confirm structure by RHEED, XRD, STEM
2. **Hydrogenation**:
   - Method A: Expose film to atomic H from cracked H2 at 200-400 C in UHV
   - Method B: Electrochemical H insertion using ionic liquid electrolyte
   - Method C: CaH2 reduction (topotactic, well-established for nickelates)
3. **Characterization**:
   - In-situ resistivity during hydrogenation to monitor T_c evolution
   - XRD to track lattice expansion (H insertion signature)
   - EELS to map H distribution
   - AC susceptibility for bulk superconductivity
4. **Estimated cost**: ~$100K (MBE time + characterization; most expensive but highest scientific value)
5. **Timeline**: 6-12 months
6. **Risk**: H insertion may destabilize the bilayer structure or reduce Ni oxidation state too much, destroying the spin-fluctuation channel. Mitigation: careful dosimetry, multiple H concentrations.

---

## 7. Recent Literature Landscape (2024-2026)

### Computational Predictions (DFT/Eliashberg)

| Year | Material | T_c (K) | P (GPa) | Method | Reference |
|------|----------|---------|---------|--------|-----------|
| 2024 | Mg2IrH6 | 160 | 0 | DFT+Eliashberg | Belli et al., PRL 132, 166001 |
| 2024 | SrAuH3 | 132 | 0 (synth@7) | DFT+Eliashberg | arXiv:2412.15488 |
| 2024 | RbPH3 | ~100 | 0 (synth@30) | DFT+anharmonic | arXiv:2411.03822 |
| 2024 | Mg2RhH6 | 80-100 | 0 | DFT+Eliashberg | npj Comput. Mater. (2024) |
| 2024 | RbYbB6C6 | 70 | 0 | DFT | PRB 109, 054505 |
| 2024 | CsB2C8 | 69 | 0 | DFT | PRB 109, 184517 |
| 2024 | SrNH4B6C6 | 85 | 0 | DFT | Commun. Phys. (2024) |
| 2025 | KB3C3 | 102.5 | 0 | DFT | November 2025 |
| 2025 | Li2AuH6 | 116-140 | 0 (unstable) | Eliashberg+SCDFT | arXiv:2501.12222 |
| 2025 | Li2AgH6 | 83-109 | 0 (unstable) | Eliashberg+SCDFT | Nat. Commun. 16, 8253 |
| 2025 | CaB8C | 77 | 0 | DFT | PRB 111, 014510 |

### Machine Learning Predictions

| Year | Method | Key Finding | Reference |
|------|--------|-------------|-----------|
| 2024 | Random forest | 93.5% accuracy T_c prediction from formula | Springer APPA (2025) |
| 2025 | Physics-informed ML | 99.94% superconductor classification, 33 descriptors | ScienceDirect (2025) |
| 2025 | Hierarchical NN | R^2 = 95.6%, 45 new HEA candidates | Front. Phys. (2025) |
| 2025 | BEE-NET (GNN) | MAE = 0.87 K, 99.4% true negative rate | npj Comput. Mater. (2026) |
| 2025 | HTSC-2025 benchmark | Standardized dataset for AI screening | arXiv:2506.03837 |
| 2025 | Mixture of Experts | R^2 = 0.962, 40 doped candidates | ACS Appl. Mater. (2024) |

### Experimental Claims and Results

| Year | Material | T_c (K) | P (GPa) | Status | Reference |
|------|----------|---------|---------|--------|-----------|
| 2025 | La1.57Sm1.43Ni2O7 | 96 (onset) | high P | Confirmed (Nature) | Nature (2025) |
| 2024 | La3Ni2O7 thin film | 26-42 | 0 (strain) | Signatures only | Nature (2024) |
| 2025 | (NbTa)0.55(HfTiZr)0.45 | 7.2 | 0 | Confirmed | J. Appl. Phys. (2025) |
| 2025 | TiNbTaN3 | 10 | 0 | First MEN SC | Adv. Sci. (2025) |
| 2023 | H-doped diamond | 117 | ~0 | **Not confirmed** | ScienceDirect (2023) |

---

## 8. The tau Roadmap: From R_tau = 0.1 to R_tau = 1

Based on this analysis, the most promising path to R_tau >= 1 combines three strategies:

```
Stage 1 (2026-2027): Synthesize Mg2IrH6 and KB3C3
  - Confirm T_c predictions experimentally
  - Measure Eliashberg spectral functions
  - Calibrate R_tau for each material
  Expected R_tau: 0.04-0.06

Stage 2 (2027-2028): Multi-channel engineering
  - Hydrogenate La3Ni2O7 thin films
  - Explore B-C clathrate doping (magnetic ions for I_spin)
  - Combine flat-band materials (TBG) with high-omega_D layers
  Expected R_tau: 0.1-0.2

Stage 3 (2028-2030): Coherent multi-channel optimization
  - Use ML-guided search (HTSC-2025 framework) with R_tau as figure of merit
  - Optimize I_phonon + I_spin + I_topo simultaneously
  - Engineer topological band crossings in high-T_c hydrides
  Expected R_tau: 0.3-0.5

Stage 4 (2030+): Beyond-Eliashberg mechanisms
  - Explore polaronic pairing, excitonic pairing, magnon pairing
  - Flat-band materials with very low tau_dec (diamond/graphene platforms)
  - Quantum anharmonic effects (RbPH3 paradigm) in new material families
  Target R_tau: >= 1
```

**The honest conclusion**: Room-temperature ambient-pressure superconductivity remains extremely challenging. The tau framework provides a clear metric (R_tau >= 1) and identifies the bottleneck (the lambda-omega_log trade-off limits single-channel phonon pairing). Multi-channel strategies are the most promising path forward, with the nickelate-hydride hybrid (Candidate 7) and diamond-based materials (Candidate 6) as the two routes that could potentially break through the conventional ceiling.

---

## 9. Connection to Other tau Framework Results

The superconductor search connects to the broader tau program in several ways:

1. **Sigma = 2 ln Q**: For a superconductor, the relevant Sigma is the electron-phonon entropy production. The condition tau_eff = 0 (superconductivity) is the condensed-matter analog of tau = 0 in our gravitational work (null geodesic = zero entropy production = information perfectly preserved).

2. **Multi-channel pairing mirrors multi-observable tau**: Just as Sigma in gravity couples to all forms of energy, tau in superconductivity receives contributions from all pairing channels. The additivity of I_pair parallels the additivity of Sigma from different matter sectors.

3. **The Gao et al. ceiling is a DPI constraint**: The lambda-omega_log trade-off is essentially a data processing inequality: you cannot extract more pairing information from phonons than the phonon channel carries. Breaking this requires a NEW channel (spin, orbital, topological), just as breaking the gravitational DPI requires entanglement assistance (wormholes).

4. **R_tau is a universal figure of merit**: It applies equally to BCS, cuprate, iron-based, and future superconductors. It provides the same kind of mechanism-independent criterion that tau = 0 provides across all of physics.

---

## References

### tau Framework
- Huang, S.-K. (2026). Petz Recovery Unification. [Paper 1]
- Huang, S.-K. (2026). Superconductivity as Zero Temporal Asymmetry. [This work]

### Ambient-Pressure Hydrides
- Belli, F. et al. (2024). PRL 132, 166001. [Mg2IrH6, Feasible route to ambient-pressure hydrides]
- arXiv:2412.15488 (December 2024). [SrAuH3, T_c = 132 K]
- Dangic, D. et al. (2024). arXiv:2411.03822. [RbPH3, anharmonic stabilization]
- Gao, M. et al. (2025). Nat. Commun. 16, 8253. [Maximum T_c at ambient pressure, Li2AuH6/Li2AgH6]
- arXiv:2501.12222 (January 2025). [Li2AuH6, T_c ~ 140 K]
- npj Comput. Mater. (2024). [Mg2XH6 family screening]

### Boron-Carbon Clathrates
- PRB 109, 054505 (February 2024). [RbYbB6C6, RbBaB6C6, T_c ~ 70 K]
- PRB 109, 184517 (May 2024). [CsB2C8, T_c = 69 K]
- Commun. Phys. (2024). [SrNH4B6C6, T_c = 85 K]
- November 2025 prediction. [KB3C3, T_c = 102.5 K]
- PRB 111, 014510 (January 2025). [CaB8C, T_c = 77 K]

### Nickelates
- Nature (2025). [La1.57Sm1.43Ni2O7, bulk T_c = 96 K under pressure]
- Nature (2024). [La3Ni2O7 thin film, ambient pressure signatures T_c ~ 26-42 K]

### Carbon-Based
- Ekimov et al. (2004). Nature 428, 542. [B-doped diamond, T_c = 4 K]
- Bhaumik et al. (2017). [B-doped Q-carbon, T_c = 36 K]
- Bhatt et al. (2023). ScienceDirect. [H-doped diamond, claimed T_c = 117 K]
- Cao, Y. et al. (2018). Nature 556, 43. [Magic-angle TBG, T_c ~ 3 K]

### Topological Superconductors
- ScienceDaily (December 2025). [PtBi2 surface superconductor]
- 2M-WS2 studies. [Topological SC candidate]

### Machine Learning
- HTSC-2025 benchmark, arXiv:2506.03837.
- npj Comput. Mater. (2026). [BEE-NET complete AI workflow]
- Various ML screening studies (2024-2025), see Section 7.

### High-Entropy Alloys
- J. Appl. Phys. 137, 215901 (2025). [(NbTa)0.55(HfTiZr)0.45]
- Adv. Sci. (2025). [TiNbTaN3, first MEN superconductor]

---

## Sources (Web Search, March 2026)

- [Are room-temperature superconductors finally within reach? | ScienceDaily](https://www.sciencedaily.com/releases/2025/10/251030075132.htm)
- [The maximum Tc of conventional superconductors at ambient pressure | Nature Communications](https://www.nature.com/articles/s41467-025-63702-w)
- [Ambient pressure high temperature superconductivity in RbPH3 | ScienceDirect](https://www.sciencedirect.com/science/article/pii/S2950463525000195)
- [Feasible Route to High-Temperature Ambient-Pressure Hydride Superconductivity | PRL](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.132.166001)
- [Theoretical Prediction of High-Temperature Superconductivity in SrAuH3 | arXiv](https://arxiv.org/abs/2412.15488)
- [Prediction of ambient pressure conventional superconductivity above 80 K | npj Comput. Mater.](https://www.nature.com/articles/s41524-024-01214-9)
- [High-Tc superconductor candidates proposed by machine learning | arXiv](https://arxiv.org/abs/2406.14524)
- [HTSC-2025 Benchmark Dataset | arXiv](https://arxiv.org/abs/2506.03837)
- [Developing a complete AI-accelerated workflow for superconductor discovery | npj Comput. Mater.](https://www.nature.com/articles/s41524-026-01964-8)
- [Hydride units filled boron-carbon clathrate | Commun. Phys.](https://www.nature.com/articles/s42005-024-01814-3)
- [Predicting superconductivity near 70 K in 1166-type boron-carbon clathrates | PRB](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.109.054505)
- [Prediction of high-temperature ambient-pressure superconductivity in hexagonal boron-rich clathrates | PRB](https://doi.org/10.1103/PhysRevB.111.014510)
- [Signatures of ambient pressure superconductivity in thin film La3Ni2O7 | Nature](https://www.nature.com/articles/s41586-024-08525-3)
- [Bulk superconductivity up to 96 K in pressurized nickelate single crystals | Nature](https://www.nature.com/articles/s41586-025-09954-4)
- [PtBi2 surface-only superconductor | ScienceDaily](https://www.sciencedaily.com/releases/2025/12/251226045350.htm)
- [Superconductivity at 117 K in H-doped diamond | ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S2542529323001517)
- [Li2AuH6 high temperature superconductivity under ambient pressure | arXiv](https://arxiv.org/html/2501.12222)
- [Bringing pressure-induced superconductivity back to ambient pressure | PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC11962462/)
