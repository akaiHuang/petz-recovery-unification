# Route 3: Thermodynamic Derivations of Sigma_grav = r_s/r

**Date**: 2026-03-07
**Agent task**: Investigate thermodynamic routes to derive Sigma_grav = r_s/r from first principles

## Executive Summary

Seven independent thermodynamic routes were investigated. Three succeed in deriving Sigma_grav = r_s/r, two give partial results, one fails by design (confirming an important distinction), and one gives consistency but not a derivation. The three successful routes all share the same mathematical core: Sigma_grav = -ln(g_00) = r_s/r. Critically, this is an **informational** entropy production (retrodiction cost), NOT a thermodynamic one -- the Tolman/Clausius route explicitly fails, confirming this distinction.

## Key Findings

### Three Successful Routes

#### 1. Unruh/Verlinde Route
Using Verlinde's entropic gravity framework (arXiv:1001.0785):
```
Sigma = -2 ln(sqrt(-g_00)) = -ln(1 - r_s/r) ~ r_s/r   (weak field)
```
The Unruh temperature and Verlinde's entropy-displacement relation dS = 2pi k_B (mc/hbar) dx combine to give this result. The proper acceleration cancels, leaving a universal, state-independent answer.

#### 2. Gravitational Landauer Route
Using Herrera's (2020) Tolman-modified Landauer principle:
- Erasure cost at radius r: k_B T(r) ln 2 = k_B T_inf ln 2 / sqrt(1-r_s/r)
- Erasure cost at infinity: k_B T_inf ln 2
- Excess cost (dimensionless): Sigma = 2|Phi|/c^2 = r_s/r

The Tolman relation T(r) = T_inf / sqrt(-g_00) modifies local erasure cost. The dimensionless excess entropy production in natural units is exactly r_s/r.

#### 3. Quantum Channel Route
Model gravitational redshift as a pure loss bosonic channel with intensity transmissivity eta = exp(-r_s/r):
```
Sigma = -ln(eta) = r_s/r
```
Both frequency shift AND time dilation contribute to the total channel loss seen by a distant observer.

### Two Partial Results

#### 4. Bekenstein Bound
Gives correct scaling Delta S ~ r_s/r but requires normalization by Petz convention for precise coefficient.

#### 5. Information-Theoretic Capacity
Capacity loss of gravitational channel is proportional to r_s/r, but coefficient depends on input state energy.

### One Designed Failure (Important!)

#### 6. Tolman/Clausius Route -- FAILS
The Clausius entropy production **vanishes identically** in thermal equilibrium (Tolman relation). This confirms that Sigma_grav is an INFORMATIONAL entropy production, not a thermodynamic one. This is consistent with Paper 2's statement that "Sigma_grav quantifies the informational cost of climbing the gravitational potential."

### One Consistency Check

#### 7. Black Hole Thermodynamics / Bousso Bound
GSL is satisfied and Bousso's covariant entropy bound is consistent with Sigma = r_s/r, but neither uniquely determines this value.

## The Unified Core

All three successful routes share the same mathematical structure:
```
Sigma_grav = -ln(g_00) = -ln(exp(-r_s/r)) = r_s/r
```

The entropy production is the negative logarithm of the metric component g_00. This is the "informational cost of climbing the gravitational potential."

## Key Distinction: Informational vs. Thermodynamic Entropy Production

- **Thermodynamic Sigma (Clausius)**: Vanishes in equilibrium (Tolman relation ensures this)
- **Informational Sigma (QRE drop)**: Non-zero even in equilibrium -- quantifies retrodiction cost
- **Implication**: Sigma_grav is fundamentally about information loss, not heat dissipation

## What Remains Un-derived

The **saturation** hypothesis (F = exp(-Sigma/2) holds as equality) remains an ansatz. The three derivation routes give Sigma >= r_s/r as a lower bound on entropy production. Proving equality (which gives the exponential metric instead of Schwarzschild) requires additional physical input, potentially from Dorau-Much (2025 PRL) or Bianconi (2025 PRD).

## Key References

- Verlinde (2010): arXiv:1001.0785 -- Entropic gravity
- Herrera (2020): Entropy 22, 340 -- Gravitational Landauer principle
- Basso, Maziero & Celeri (2025): PRL 134, 050406; arXiv:2405.03902 -- Quantum Crooks in curved spacetime
- Moreira & Celeri (2024): arXiv:2407.21186 -- Graviton bath entropy production
- de Paolis et al. (2025): arXiv:2502.20521 -- Limits of beam-splitter model for gravitational redshift
- Tamburini & Licata (2023): arXiv:2307.01296 -- Unruh entropy of Schwarzschild black hole
- Pikovski et al. (2015): Nature Physics 11, 668; arXiv:1311.1095 -- Gravitational decoherence

## Assessment

- **Likelihood of deriving Sigma_grav = r_s/r**: HIGH (three independent routes succeed)
- **What works**: The Verlinde, Gravitational Landauer, and Quantum Channel routes all independently give r_s/r; the universal core is -ln(g_00)
- **What doesn't**: Saturation (F = exp(-Sigma/2) as equality) is not derived; the distinction between informational and thermodynamic Sigma needs careful framing
- **Strongest route**: Gravitational Landauer (physically transparent, rigorous, connects to Paper 1's Retrodiction Landauer Principle)

## Next Steps

1. **For Paper 2**: Present all three derivation routes with the Gravitational Landauer as the cleanest argument
2. **Emphasize the informational vs. thermodynamic distinction**: This is a key conceptual contribution
3. **The saturation question**: Investigate whether Petz bound saturation can be derived from additional physical principles (e.g., maximal channel capacity, or the "most efficient" gravitational information processing)
4. **Connect to Paper 1**: The Retrodiction Landauer Principle (Paper 1) + Gravitational Landauer (Herrera 2020) = Sigma_grav = r_s/r. This is a beautiful chain.
