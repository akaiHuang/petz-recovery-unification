# Complementary Uncertainty: PROVED (2026-03-12)

## Result

**Theorem (Complementary Observer Uncertainty):**
For d+1 MUB observers in prime-power dimension d and any pure state |ψ⟩:

Σᵢ τ_{Oᵢ} ≥ (d+1) − √(2(d+1)) =: g(d) > 0

The bound is **tight** (achievable).

## Proof (3 elementary ingredients)

1. **Fidelity formula**: For dephasing in basis i, F_i = √P_i where P_i = Σ_k |⟨e_k|ψ⟩|⁴
   Therefore τ_i = 1 − √P_i

2. **2-design sum rule**: Complete MUBs form a 2-design (Klappenecker-Rotteler 2005)
   → Σᵢ P_i = 2 for ANY pure state (state-independent!)

3. **Cauchy-Schwarz**: Σᵢ √P_i ≤ √((d+1)·2) = √(2(d+1))
   → Σᵢ τ_i ≥ (d+1) − √(2(d+1))

## Values

| d | K=d+1 | g(d) | Numerical |
|---|-------|------|-----------|
| 2 | 3 | 3−√6 | 0.5505 |
| 3 | 4 | 4−2√2 | 1.1716 |
| 4 | 5 | 5−√10 | 1.8377 |
| 7 | 8 | 4 | 4.0000 |
| ∞ | — | ~d | g(d)/d → 1 |

**Note**: g(2) = 0.5505, NOT 1/2 as originally conjectured. The exact bound is 10% stronger.

## Physical Interpretation

Temporal analogue of Maassen-Uffink: total irreversibility across complementary observers has a positive floor. The arrow of time cannot be eliminated from all perspectives simultaneously.

"A state timeless for one observer is necessarily time-full for all complementary observers."

## Impact on Papers

- Paper 5: Conjecture 1 → Theorem (thm:comp)
- Paper 5 supplemental: S11 "Evidence for Conjecture 1" → "Proof of Theorem"
- Paper 4: conjecture → theorem in observer-dependent section

## File

calculations/complementary_observer_uncertainty_proof.tex (4 pages, compiled)
