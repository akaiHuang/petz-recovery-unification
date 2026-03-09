# Research: Scale-Dependent Σ and Beyond r_s/r

## Date: 2026-03-06
## Relevance: Paper 3 (core) + Paper 4

---

## 1. KEY FINDING: Kumar (2025) — Logarithmic Correction from First Principles

### arXiv:2509.05246

**What**: From quantum field theory (marginal IR running of Newton's constant with anomalous dimension η = 1), derives a logarithmic correction to Newtonian potential:

```
Φ(r) = −G_N M/r + (2G_N M k_*/π) · ln(r/r₀)
```

In Σ language:
```
Σ(r) = r_s/r + α · ln(r/r_c)
```
where α is dimensionless, r_c ~ 37 kpc is a universal crossover scale.

**Key properties**:
- **First principles**: derived from QFT, not phenomenological
- **Universal**: regulator-independent, same for all galaxies
- **Automatically screened** at solar-system distances
- **Validated**: Gubitosi et al. (2024, arXiv:2403.00531) tested against SPARC catalog, fit quality comparable to NFW dark matter halos

### Physical Origin
The running of G comes from quantum gravity loop corrections. At energy scale k (corresponding to distance r ~ 1/k):
```
G(k) = G_N / (1 + ω G_N k²/(6π))    (asymptotic safety)
```
In the IR (large r, small k), the marginal running with η = 1 gives the logarithmic correction.

---

## 2. Asymptotic Safety and Galaxy Rotation Curves

### 2.1 The Program (Reuter, Percacci)
In asymptotic safety, G runs with energy/momentum scale:
```
G(k) = G_N [1 + ω(G_N k²) + ...]
```
UV: G → 0 (safe). IR: G may grow logarithmically.

### 2.2 Galaxy-Scale Tests
**Gubitosi, Piattella & Casarini (2024, arXiv:2403.00531)**:
- Applied running G to SPARC galaxy rotation curves
- Found: fit quality comparable to NFW halos
- Crossover scale r_c ~ 37 kpc (universal)
- Only 1-2 universal parameters vs NFW's 2-3 per galaxy

**Rodrigues, Piattella, Fabris (2022)**:
- Running G from asymptotic safety applied to spiral galaxies
- Can reproduce flat rotation curves
- Tension: requires rather large running at galactic scales

### 2.3 Connection to τ Framework
The running G means:
```
Σ_effective(r) = 2G(r)M/(c²r) = (r_s/r) · [1 + f(r)]
```
where f(r) encodes the running. This is NOT an ad hoc addition — it's the same Σ = 2|Φ|/c² with the full quantum-corrected potential.

---

## 3. Non-Local Gravity

### 3.1 Deser-Woodard Model
Adds R(1/□)R to Einstein-Hilbert action. Effects:
- Negligible at solar-system scales
- Non-trivial at cosmological scales
- Can mimic dark energy without Λ

### 3.2 Connection to Σ
Non-local corrections = the gravitational channel N has memory (non-Markovian).
In τ language: non-locality → non-Markovian N → backflow of information → modified Σ.

---

## 4. f(R) Gravity and Exponential Metric

### 4.1 Nath-Sarma (2025)
The exponential metric is a solution in logarithmic f(R) gravity:
```
f(R) = R + α R ln(R/R₀)
```

The extra scalar DOF (scalaron) has mass m_s:
- r ≪ 1/m_s: extra force active (galactic scale modification)
- r ≫ 1/m_s: returns to GR (solar system safe)

### 4.2 Chameleon Mechanism
In f(R) theories, the scalaron mass depends on local density:
- High density (solar system): m_s large → short range → GR recovered
- Low density (galaxy outskirts): m_s small → long range → modified gravity

This is exactly what's needed: automatic solar-system screening + galactic-scale effects.

---

## 5. Entanglement Entropy at Different Scales

### 5.1 UV (Strong Field)
Σ = r_s/r: area-law contribution from local gravitational entropy.

### 5.2 IR (Galactic Scale)
Two possible sources:
- **Running G**: quantum loop corrections → logarithmic Σ
- **Volume-law entanglement**: de Sitter background → Verlinde's apparent dark matter

### 5.3 The Scale Hierarchy
```
Planck → Strong field → Solar system → Galaxy → Cluster → Hubble
  UV          r_s/r        r_s/r          +ln      +ln+dS    dS
```
Each scale adds physics. But in the unified picture, it's ONE Σ evaluated with the full quantum-corrected gravitational channel.

---

## 6. Emergent MOND from Quantum Gravity

### 6.1 Milgrom's a₀ ~ cH₀
The MOND acceleration scale a₀ ~ 1.2×10⁻¹⁰ m/s² is suspiciously close to cH₀.
This connects LOCAL physics (c) to COSMOLOGY (H₀).
In τ language: the de Sitter background (H₀) sets a minimum information loss rate.

### 6.2 Smolin (2017, arXiv:1704.00780)
MOND from quantum gravity in de Sitter:
- Below T_dS, equivalence principle breaks down
- Gravitational/inertial mass ratio becomes environment-dependent
- Naturally reproduces MOND with a₀ ~ cH₀

### 6.3 Key Question for Paper 3
Can the single Σ = D(ρ_spacetime ‖ ρ_matter), when computed in a de Sitter background (not flat), naturally produce the MOND interpolation function?

---

## 7. Mathematical Requirements for Dark Matter Explanation

For the τ framework to explain dark matter, need:
1. ✅ Don't change PPN (solar system fine) — automatic if running G screened
2. ✅ ~200% effects at r ~ 10 kpc — logarithmic correction provides this
3. ✅ Physical motivation from QI — QRE + DPI
4. ❓ Universal (same parameters for all galaxies) — Kumar's r_c ~ 37 kpc is universal
5. ❓ Works for galaxy clusters — not yet tested
6. ❓ Explains CMB power spectrum — requires full cosmological calculation

---

## 8. The Three-Term Approximation (for reference)

As an approximation of the unified Σ, in different regimes:
```
Σ_total ≈ r_s/r + α · ln(r/r_c) + β · (r/L_dS)
```
- Term 1: Local gravity (dominant near mass)
- Term 2: Quantum running of G (dominant at galactic scales)
- Term 3: De Sitter volume-law entropy (dominant at cosmological scales)

**IMPORTANT**: This is an approximation, NOT the fundamental equation.
The fundamental equation is Σ = D(ρ_spacetime ‖ ρ_matter).
The three terms emerge as limits of this one equation.

---

## Key References

- **Kumar 2025**: arXiv:2509.05246 (MOST IMPORTANT for Paper 3)
- **Gubitosi et al. 2024**: arXiv:2403.00531 (SPARC validation)
- Rodrigues et al. 2022: Running G spiral galaxies
- Nath-Sarma 2025: Exponential metric in f(R)
- Smolin 2017: arXiv:1704.00780 (MOND from QG)
- Milgrom 2020: arXiv:2001.09729 (a₀-cosmology connection)
- Verlinde 2016: arXiv:1611.02269
- Bianconi 2025: arXiv:2408.14391
