# Four-Paper Architecture: From τ to Unified Gravity

## Author: Sheng-Kai Huang
## Date: 2026-03-06
## Status: Active research program

---

## Core Philosophy

> 比起找出新的數學或方程式，提供一個新的論點串起已知的東西，讓事情變「對」，更有意義。
> 我們不是在不同尺度加不同項（補丁）。我們要找出一個統一的 Σ，在不同尺度自然展現不同行為。

## The Unified Equation (Goal)

```
Σ = D(ρ_spacetime ‖ ρ_matter)
```

One equation. Different environments → different solutions → different physics.

- Near a black hole → exponential metric (Paper 2)
- At galactic scales → flat rotation curves without dark matter (Paper 3)
- At cosmological scales → accelerating expansion (Paper 4)

NOT: Σ = term1 + term2 + term3 (patchwork)
YES: Σ = one QRE, solved in different boundary conditions

---

## Paper 1: The Tool (COMPLETED)

**Title**: "Universal Recovery: The Petz Map as Retrodiction Functor, Error Corrector, and Entropy Bound"

**Core**: τ = 1 − F, equivalence chain, F ≥ exp(−Σ/2)
**Status**: Written, on GitHub, awaiting arXiv endorsement (Wilde contacted)
**Risk**: Low (mathematically rigorous)
**Validation**: 1620/1620 numerical points, 3 falsifiable predictions all supported

## Paper 2: Strong Field (IN PROGRESS)

**Title**: "Information-Theoretic Origin of the Exponential Metric: Gravitational Redshift as Petz Recovery Fidelity"

**Core**: Σ_grav = r_s/r → g₀₀ = −exp(−r_s/r) → no event horizons
**Status**: 9-page draft exists (`papers/paper2_exponential_metric.tex`), needs revision to ~6 pages
**Risk**: Medium (saturation hypothesis)
**Key predictions**:
1. GW echoes at ~2.5ms (89× shorter than Planck-wall models)
2. Neutron star redshift 19% different from GR
3. QNM frequency shifted ~4.4%
4. Shadow 4.6% larger
5. ISCO 5.6% larger

**What Paper 2 does NOT do**: dark matter, dark energy, cosmology

## Paper 3: Weak Field (PLANNED)

**Title** (working): "Galactic Dynamics from Quantum Information Recovery"

**Core**: The same Σ = D(ρ‖σ) − D(N(ρ)‖N(σ)), solved at galactic scales, naturally produces flat rotation curves.

**Approach**: NOT "add a logarithmic term." Instead: compute Σ for the complete gravitational channel (including quantum corrections and de Sitter background) at galactic scales.

**Key resources identified**:
- Kumar (2025, arXiv:2509.05246): QFT first-principles logarithmic correction to Newton's potential
- Gubitosi et al. (2024, arXiv:2403.00531): SPARC galaxy fitting with running G
- Verlinde (2016, arXiv:1611.02269): Volume-law entanglement → apparent dark matter
- Ghari & Haghi (2026, arXiv:2601.01715): Verlinde favored over MOND at 5.2σ in dwarf spheroidals
- Brouwer et al. (2021, arXiv:2106.11677): KiDS-1000 weak lensing support

**The key question**: Can D(ρ_spacetime ‖ ρ_matter) at galactic scales, computed from first principles (not as an ansatz), reproduce the observed rotation curves?

**Risk**: Medium-high (requires significant new calculation)
**Validation**: SPARC database (175 galaxies, public data)

## Paper 4: Grand Unification (FUTURE)

**Title** (working): "Gravity as Quantum Information Recovery: A Unified Framework"

**Core**:
```
Σ = D(ρ_spacetime ‖ ρ_matter)
```
One equation that:
- Reproduces Paper 2 (strong field) as a limit
- Reproduces Paper 3 (weak field) as a limit
- Extends to cosmology (FRW via Kodama vector)
- Connects to Bianconi's "metric = density matrix" framework

**Key resources identified**:
- Bianconi (2025, PRD 111, 066001): gravity action = QRE
- Bianconi (2025, arXiv:2510.22545): FRW cosmological solutions in GfE
- Dorau-Much (2025 PRL, arXiv:2510.24491): QRE → Einstein equations
- Basso et al. (2025, PRL 134, 050406): Quantum Crooks in curved spacetime
- Padmanabhan CosMIn: finite cosmic information requires Λ > 0
- Kodama vector for dynamic spacetimes

**Risk**: High (most ambitious)
**Timeline**: After Papers 1-3 are accepted

---

## The Narrative Arc

> Paper 1: "Time's asymmetry can be quantified by a single number τ."
> Paper 2: "τ determines the geometry of spacetime near black holes."
> Paper 3: "τ also explains why galaxies rotate too fast — no dark matter needed."
> Paper 4: "From black holes to cosmic expansion, it's all one Σ."

Each step is bolder than the last, but each stands on the previous one's foundation.

---

## Risk Management

- Paper 1 is independent: even if 2-4 fail, Paper 1 stands
- Paper 2 Layer 1 (τ < 1, no event horizons) is independent of Layer 2 (specific metric)
- Paper 3 is independently testable against SPARC data
- Paper 4 only attempted after 1-3 are validated

---

## Key Distinction: NOT Patchwork

The three-term approximation Σ ≈ r_s/r + α·ln(r/r_c) + β·(r/L_dS) is an APPROXIMATION of the unified Σ = D(ρ_spacetime ‖ ρ_matter) in different regimes. Like how:

- F = ma is ONE equation
- Applied to planets → Kepler orbits
- Applied to springs → SHM
- Applied to fluids → Navier-Stokes

Not three forces. One force, three solutions.

Similarly: not three Σ terms. One Σ, three regimes.
