# Research: Verlinde / Entropic Gravity / MOND ↔ τ Framework

## Date: 2026-03-06
## Relevance: Paper 3 (weak field) + Paper 4 (unification)

---

## 1. Verlinde's Emergent Gravity

### 1.1 The Two Papers

**Verlinde 2010** (arXiv:1001.0785): Gravity as entropic force: F = T(dS/dx).
Reproduces Newton exactly but adds nothing new.

**Verlinde 2016** (arXiv:1611.02269, SciPost Phys. 2, 016 (2017)): The crucial extension.

Key insight: in de Sitter space, entanglement entropy has BOTH area law AND volume law:
```
S_entanglement = S_area + S_volume
```
- S_area = A/(4Gℏ) (standard Bekenstein-Hawking)
- S_volume = thermal volume-law from de Sitter entanglement at T_dS = H₀/(2π)

The competition between area and volume law means at sub-Hubble scales, de Sitter vacuum does NOT fully thermalize → **entropy displacement** when baryonic matter is present.

### 1.2 Verlinde's Dark Matter Formula

For spherically symmetric baryonic mass M_B(r):
```
M_D²(r) = (cH₀/(6G)) · r² · d/dr[r · M_B(r)]
```

For point mass:
```
M_D(r) = √(cH₀M_B r / (6G))
```

At large r where g_Newton ≪ a₀:
```
g ~ √(a₀ · g_Newton)    ← MOND interpolation!
```
with a₀ ~ cH₀/6 ~ 1.1×10⁻¹⁰ m/s² (close to Milgrom's 1.2×10⁻¹⁰)

### 1.3 Current Observational Status (2024-2026)

**Supporting:**
- Yoon, Park & Hwang (2023, CQGra; arXiv:2206.11685): 175 SPARC galaxies, good fit
- **Ghari & Haghi (2026, arXiv:2601.01715): 23 dwarf spheroidals, Verlinde favored over MOND at 5.2σ**
- Brouwer et al. (2017, MNRAS 466; arXiv:1612.03034): KiDS-450 weak lensing, no free parameters
- Brouwer et al. (2021, A&A 650; arXiv:2106.11677): KiDS-1000 extends RAR by 2 decades
- Chae (2024, JKAS): Wide binary stars show ~40% gravity boost at low acceleration

**Challenging:**
- Lelli, McGaugh & Schombert (2017, ApJ 836; arXiv:1702.04355): Requires lower M/L ratios, predicts unobserved residual correlations
- Hossenfelder (2017, PRD 95; arXiv:1703.01415): Covariant version reduces to GR
- Dai & Stojkovic (2017, JHEP; arXiv:1710.00946): Internal inconsistencies

**Bottom line**: Alive but controversial. Works well for dwarf spheroidals and weak lensing, theoretical consistency issues.

---

## 2. Mathematical Bridge: Verlinde ↔ τ Framework

### 2.1 Structural Parallel

| τ Framework | Verlinde Framework |
|---|---|
| Σ (entropy production) | S_displacement (entropy displaced by matter) |
| F = exp(−Σ/2) | c_eff/c = exp(−Σ_grav/2) |
| τ > 0 = information loss | Entropy displacement = "apparent dark matter" |
| Petz recovery map | "Elastic" restoration of entanglement |

### 2.2 The Unified Σ Approach (NOT patchwork)

**Wrong approach**: Σ_total = r_s/r + Σ_dS (patchwork)

**Right approach**: Σ = D(ρ_spacetime ‖ ρ_matter)
- In asymptotically flat spacetime → only r_s/r behavior
- In de Sitter background → naturally includes volume-law contribution
- The "extra term" emerges from the same QRE, not added by hand

### 2.3 Verlinde's Formula in τ Language

The additional gravitational potential from volume-law entanglement:
```
Φ_dark(r) = −√(cH₀GM_B/6) · ln(r/r₀) + const
```

In τ language: this is additional information loss from the de Sitter background being disturbed by matter.

**Dark matter = additional temporal asymmetry (τ > 0) from de Sitter vacuum's volume-law entanglement being displaced by baryonic matter.**

### 2.4 Scale Comparison

For Milky Way (M_B ~ 6×10¹⁰ M_sun):
- At r ~ 8 kpc: Σ_grav ~ 3×10⁻⁶, Σ_dS ~ 10⁻⁶ (comparable!)
- At r ~ 30 kpc: Σ_grav ~ 10⁻⁶, Σ_dS ~ 2×10⁻⁶ (dS dominates → flat curve)
- At solar system: Σ_dS ~ 10⁻¹⁴ (negligible → passes all tests)

---

## 3. Deeper Connections via Quantum Relative Entropy

Three independent lines converge on QRE:
1. **Jacobson (1995, PRL 75, 1260)**: Einstein eq = δQ = TdS for local Rindler horizons
2. **Jacobson (2016, PRL 116, 201101; arXiv:1505.04753)**: Einstein eq = maximal vacuum entanglement
3. **Dorau-Much (2025, PRL; arXiv:2510.24491)**: Einstein eq from QRE via modular theory
4. **Bianconi (2025, PRD 111, 066001; arXiv:2408.14391)**: Gravity action = QRE, metric = density matrix

The τ framework's Σ IS the QRE drop. The DPI is the universal structure.

---

## 4. MOND from Information Theory

- **Smolin (2017, PRD 96; arXiv:1704.00780)**: MOND from quantum gravity regime in de Sitter
- **Milgrom (2020, arXiv:2001.09729)**: a₀ ~ cH₀ ~ c²Λ^{1/2} ("FUNDAMOND")
- **Relativistic MOND from modified entropic gravity (2025, arXiv:2511.05632)**: Temperature-dependent corrections → MOND interpolation

---

## 5. Padmanabhan's Thermodynamic Gravity

**Holographic Equipartition** (2012, arXiv:1206.4916):
```
dV/dt = L_P² (N_sur − N_bulk)
```
N_sur = N_bulk → de Sitter equilibrium. N_sur > N_bulk → expansion.

**CosMIn** (2017, arXiv:1703.06144): Finite cosmic information REQUIRES Λ > 0.

**Connection to τ**:
- N_sur − N_bulk > 0 ↔ τ > 0
- dV/dt > 0 (expansion) ↔ Σ_total > 0
- CosMIn = total τ integrated over cosmic history

---

## 6. Key Open Questions

1. Can D(ρ_spacetime ‖ ρ_matter) in de Sitter background reproduce Verlinde's M_D² formula?
2. Does Bianconi's G-field correspond to Verlinde's apparent dark matter?
3. Can Dorau-Much's QRE → Einstein framework be extended to include volume-law corrections?
4. What is the covariant formulation? (Hossenfelder's attempt failed)

---

## References (Key Papers)

- Verlinde 2016: arXiv:1611.02269
- Ghari & Haghi 2026: arXiv:2601.01715 (5.2σ over MOND)
- Brouwer et al. 2021: arXiv:2106.11677 (KiDS-1000)
- Lelli et al. 2017: arXiv:1702.04355 (tension)
- Smolin 2017: arXiv:1704.00780
- Milgrom 2020: arXiv:2001.09729
- Hossenfelder 2017: arXiv:1703.01415
- Jacobson 1995: PRL 75, 1260
- Jacobson 2016: arXiv:1505.04753
- Dorau-Much 2025: arXiv:2510.24491
- Bianconi 2025: arXiv:2408.14391
- Padmanabhan 2012: arXiv:1206.4916
- McGaugh et al. 2016: arXiv:1609.05917 (RAR)
- Relativistic MOND 2025: arXiv:2511.05632
