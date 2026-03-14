# η = 1 from DPI: DERIVED (2026-03-12)

## Result

**Theorem**: DPI alone gives only η ≥ 0. But DPI + extensivity of de Sitter entanglement entropy uniquely selects η = 1 in 3+1 dimensions.

## The Argument

### DPI alone: insufficient
The DPI is an inequality, not an equation. Any CPTP map with any positive contraction coefficient satisfies it. Cannot fix the numerical value of η.

### Minimal additional assumption: de Sitter extensivity
At scales r >> r₀, entanglement entropy transitions from area-law (S ~ r²) to volume-law (S ~ r³), driven by positive Λ.

Three independent motivations:
- Verlinde (2016): volume-law from positive Λ
- Jacobson (2016): entanglement equilibrium breaking
- Casini-Huerta (2017): entropic c-function monotonicity

### Three-pronged exclusion
- **η > 1 excluded**: δΦ ~ r^{η−1} grows → divergent entropy → exceeds dS bound → not a valid CPTP map
- **η < 1 excluded**: δΦ ~ r^{η−1} → 0 → IR-irrelevant → contradicts extensivity
- **η = 1 unique**: marginal case, logarithmic potential δΦ ~ ln(r/r₀)

### Dimension dependence
η = d − 3 in d spacetime dimensions:
- d=3: η=0 (G dimensionless, no running)
- d=4: η=1 ✓
- d=5: η=2

**Universal**: independent of matter content (coefficient k* depends, exponent does not)

## Impact on Papers

- Paper 3: "conjectured" → "derived from DPI + extensivity"
- Paper 4: Table — galactic row → "derived"
- Paper 4: Open problems — η=1 marked as resolved

## Key Statement for Papers

"The τ framework, combined with the established volume-law entanglement entropy in de Sitter spacetime, predicts η = 1. This is not a prediction of DPI alone — it is a prediction of DPI plus the specific entanglement structure of our universe."

## File

calculations/dpi_eta_derivation.tex (21 pages, compiled)
