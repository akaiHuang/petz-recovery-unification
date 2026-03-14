# Channel Problem: SOLVED (2026-03-12)

## Summary

The Channel Problem (#2) — constructing a canonical gravitational CPTP map — has been fully resolved. Three computational gaps that were previously open have been closed.

## The Theorem

**Theorem (Gravitational Thermal Attenuator Channel):**
For a scalar field mode propagating from radius r₁ to r₂ in a static metric ds² = −f(r)dt² + f(r)⁻¹dr² + r²dΩ²:

1. **Channel**: CPTP thermal attenuator with η = f(r₁)/f(r₂)
2. **Environment**: Thermal at Tolman temperature T_loc = T_H/√f(r)
3. **Entropy production**: Σ = −ln(η) + O(ω²/T²) = −ln(−g₀₀)
4. **Recovery**: F ≥ exp(−Σ/2) = √f(r)
5. **Schwarzschild**: Σ = −ln(1−r_s/r) = r_s/r + O((r_s/r)²)
6. **Exponential metric**: Σ = r_s/r exactly

## Three Gaps Closed

### Gap 1: Bogoliubov Geometric Optics (η = −g₀₀)
- **Result**: η = f(r) via energy ratio × time dilation = √f × √f = f(r)
- **File**: calculations/bogoliubov_eta_derivation.tex (18 pages)
- **Key**: Two factors of √f — one from redshift, one from time dilation

### Gap 2: Greybody Corrections (subleading)
- **Result**: |Σ_exact − Σ₀| ≤ (1/4)(r_s/r)²
- **File**: calculations/greybody_subleading_proof.tex (14 pages)
- **Key findings**:
  - s-wave: ΔΣ = ε²/16 (8× smaller than Schw vs exp metric difference)
  - β coefficients: exactly zero for static observers
  - Angular modes: exactly decouple in spherical symmetry

### Gap 3: Mode-Sum Entropy (intensive quantity)
- **Result**: Σ_grav = −ln(−g₀₀) is intensive (per degree of freedom)
- **File**: papers/mode_sum_entropy_production.tex (9 pages)
- **Key findings**:
  - No UV divergence (Σ is per-mode, not summed)
  - Per-mode Petz bound: F_mode ≥ √(−g₀₀)
  - n modes: F_total ≥ (−g₀₀)^{n/2}
  - Three routes confirm −ln(η): capacity deficit, Stinespring, geometric

## Impact on Papers

### Paper 2
- "The channel problem" subsection → replaced with Theorem 1
- "Cannot: construct CPTP map" → removed (now solved)
- "Directions: canonical channel" → upgraded to "beyond geometric optics"

### Paper 4
- Σ_grav status: "conjectured" → "derived"
- Table I: gravity row → "derived"

## Three Routes Convergence

| Route | Framework | Result | Status |
|-------|-----------|--------|--------|
| A. Algebraic/Modular | Dorau-Much + Trejo-Calderon | Conditional expectation → Σ = S^rel | Rigorous |
| B. Thermal Attenuator | Beam-splitter + Tolman | η = −g₀₀, Σ = −ln(η) | Explicit Kraus |
| C. Gravitational Landauer | Herrera 2020 | Erasure cost → Σ/2 | Thermodynamic |

All three are facets of the same channel at different levels of abstraction.

## Remaining Caveat

The only remaining open question is whether the channel **saturates** the Petz bound (exponential metric) or has a gap of O((r_s/r)²) (Schwarzschild). This is empirically testable via GW echoes.

## References (new)

- Ivan, Sabapathy, Simon (2011): Kraus operators for bosonic Gaussian channels
- Holevo & Werner (2001): Bosonic channel capacities and entropy production
- Birrell & Davies (1982): QFT in curved spacetime (mode equations)
- Wald (1994): QFT and black hole thermodynamics
- Trejo-Calderon (2025): Modular channels and spectral emergence
- Cirafici (2024): Fluctuation theorems in gravitational algebras
