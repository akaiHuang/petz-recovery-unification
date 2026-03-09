# Research: Cosmological Extension of τ Framework

## Date: 2026-03-06
## Relevance: Paper 3 (background) + Paper 4 (core)

---

## 1. τ in de Sitter Space

### 1.1 De Sitter Thermodynamics
- Gibbons-Hawking temperature: T_GH = H/(2π)
- Local temperature (Akhmedov & Singleton, arXiv:2312.02292): T_local = H/π (2× T_GH)
- Gibbons-Hawking entropy: S_GH = A/(4G) = π/(GH²)

### 1.2 Natural Σ_cosmological
Parallel to black hole:
- BH: Σ_grav = r_s/r → τ = 1 at event horizon
- dS: Σ_dS = r/r_H → τ = 1 at cosmological horizon

For observer at center:
- r → 0: Σ → 0, τ → 0 (good recovery)
- r → r_H: Σ → 1, τ → 1 − exp(−1/2) ≈ 0.39
- Beyond horizon: τ → 1 (complete information loss)

### 1.3 Bianconi's GfE for de Sitter (arXiv:2510.22545, 2025)
- FRW solutions exist in Gravity-from-Entropy (GfE) framework
- k-temperature: θ_k ~ ω_k H²
- Total entropy ∝ H⁻², matching Gibbons-Hawking
- GQRE per unit volume decreases, but total entropy increases (volume grows as exp(3Ht))

---

## 2. Cosmological Expansion as a Quantum Channel

### 2.1 Bogoliubov Transformation = Quantum Channel

**Capozziello et al. (2024, arXiv:2406.19274, EPJC)**:
- Bogoliubov transformation of field modes → Gaussian bosonic thermal noise channel
- Information preservation depends on particle production rate
- Key: particle production depends on gravity theory → modified gravity preserves more information

**Barman et al. (2026, arXiv:2601.20860, QIP)**:
- Teleportation fidelity in FRW spacetime computed
- Expansion degrades fidelity through mode mixing
- Directly quantifiable via Bogoliubov coefficients

### 2.2 Entropy Production from Expansion
For massless scalar in de Sitter:
```
Σ_mode ~ (H/ω)² for ω ≫ H
```
Modes at ω ~ H: **Σ ~ O(1) per Hubble time**

---

## 3. Quantum Crooks Theorem in Curved Spacetime

### 3.1 Basso, Maziero & Celeri (PRL 134, 050406, 2025; arXiv:2405.03902)

**Main result**: Fully general-relativistic quantum detailed fluctuation theorem.

Curved-spacetime Crooks relation:
```
P_fwd(W) / P_rev(−W) = exp[β(W − ΔF)]
```

Hamiltonian includes curvature:
```
H_cm(τ) = m + p²/(2m) + m a_i(τ)x^i + (m/2) R_{τiτj}(τ) x^i x^j
```
Last term (Riemann tensor) drives system out of equilibrium → entropy production.

### 3.2 In de Sitter
QHO in de Sitter: transition probability scales as (H/ω₀)⁴.
For H ~ 10⁻⁶¹ t_P⁻¹, atomic oscillators: H/ω₀ ~ 10⁻³¹ → vanishingly small per oscillation.
But nonzero whenever H ≠ 0.

### 3.3 Observer Dependence
**Entropy production is observer-dependent**:
- Different worldlines → different R_{τiτj}
- Arrow of time rooted in causal structure
- This matches τ framework: τ depends on which channel N you're looking at

### 3.4 Related Papers
- Cai, Wang & Zhao (2024, arXiv:2407.09912): GR fluctuation theorems, covariant
- Moreira & Celeri (2024, arXiv:2407.21186): Entropy from spacetime fluctuations (graviton bath)
- PRL 134, 232302 (2025): Effective action for relativistic hydro from Crooks

---

## 4. Dark Energy from Quantum Information

### 4.1 Padmanabhan's CosMIn (arXiv:1703.06144, PLB 773, 81-85)
- CosMIn = total information transferred quantum→classical, accessible to eternal observer
- Conserved quantity
- **Finite CosMIn REQUIRES late-time accelerated expansion (Λ > 0)**
- Both Λ and primordial perturbation amplitude determined by one parameter
- Connection to τ: CosMIn = ∫τ dt; finite CosMIn requires τ → finite limit → Λ > 0

### 4.2 Bianconi's GfE Cosmology (arXiv:2510.22545, 2025)
- FRW solutions exist and are thermal
- First law: δε_k = θ_k δs_k − π_k δv
- Low-energy limit: approximates standard Friedmann
- **Emergent Λ**: dynamical, from Hamiltonian H = 2βΛ^G
- Hubble dependence: θ_k ~ ω_k H², π_k ~ (1/2)ω_k² H⁴

### 4.3 Dorau-Much (arXiv:2510.24491, PRL 2025)
- QRE between vacuum and coherent excitations on bifurcate Killing horizon = energy flux
- Under Bekenstein-Hawking → area variation → Einstein equations
- **Deepest QI→gravity connection: uses exactly the QRE underlying Petz map and τ**

### 4.4 Other Approaches
- Neukart (2024, arXiv:2409.12206): Spacetime curvature ∝ entanglement entropy
- Gomez (2020, arXiv:2002.04294): Inflation via time-dependent density matrices

---

## 5. Kodama Vector for Cosmological τ

### 5.1 Definition
In FRW (no Killing vector), the Kodama vector K^a:
- Divergence-free in any spherically symmetric spacetime
- Conserved charge = Misner-Sharp mass E_MS
- Surface gravity at apparent horizon (not event horizon)

For FRW: K^a = (1, −Hr/√(1−kr²), 0, 0) in comoving coordinates

### 5.2 F_bound via Kodama
```
F_K = |K| / |K_ref|
```
- r = 0 (observer): |K|² = −1, F_K = 1, τ = 0
- r = r_AH (apparent horizon): |K|² = 0, F_K = 0, τ = 1
- Smooth interpolation → exactly the desired cosmological τ

### 5.3 Thermodynamic Foundation
First law at FRW apparent horizon: dE = TdS + WdV
Projecting unified first law along Kodama vector → second Friedmann equation
**Friedmann equations ARE the first law of thermodynamics**

### 5.4 Extensions
- Dorau & Verch (2024, arXiv:2402.18993): Kodama-like vectors in axisymmetric spacetimes (Kerr-Vaidya-de Sitter)
- arXiv:2402.16484 (2024): Geometrical origin of Kodama vector

---

## 6. VSL Cosmology Status (2024-2025)

### meVSL Model (Lee, arXiv:2011.09274, JCAP 2021)
- c varies with scale factor a(t), related constants co-vary to preserve Lorentz invariance
- H(z) = (1+z)^{−b/4} H^{GR}(z), single parameter b
- May address Hubble tension

### Observational Tests
| Paper | Data | Result |
|---|---|---|
| Hong et al. (arXiv:2409.03248, MNRAS 534) | BAO + CC | Barrow VSL ruled out; CPL-VSL acceptable |
| arXiv:2505.15768 (2025) | DESI BAO + Pantheon+ | ~4σ tension with 2D BAO only |
| arXiv:2509.08840 (2025) | Hubble tension | meVSL can alleviate H₀ tension |

---

## 7. Key Result: Σ ~ O(1) at Hubble Scale

Three independent routes:

**Route A: Cosmological Horizon (MOST PROMISING)**
Kodama: |K|² → 0 at apparent horizon → τ → 1.
Modes at ω ~ H: Σ ~ O(1) per Hubble time.

**Route B: Cumulative Weak Fields**
```
Σ_cumulative ~ (4πGρ₀/3)r²/c² = (r/r_H)² · Ω_m
```
At r ~ r_H with Ω_m ~ 0.3: **Σ ~ 0.3**, which is O(1).

**Route C: Bianconi's GfE**
Per degree of freedom per Hubble time: Σ ~ O(1).

---

## 8. Recommendations

### Immediate
1. Define τ_cosmological via Kodama vector
2. Connect to Basso et al. Crooks theorem (Σ from curvature coupling)
3. Use Bianconi's GQRE as first-principles Σ_cosmological

### Medium-term
4. Cosmological Crooks for FRW: derive τ(z) as function of redshift
5. Model expansion as Gaussian channel, compute Petz recovery map
6. Bridge CosMIn ↔ τ: finite CosMIn requiring Λ > 0

### Speculative
7. Dark energy as self-consistent τ
8. Observer-dependent cosmological arrows of time
