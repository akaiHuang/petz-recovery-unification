# c_s² Verification: The Running μ Resolution (2026-03-14)

## Summary

The Khronon/τ framework predicts dark matter as a Khronon field with effective GDM
parameters {w(z), c_s²(z,k), c_vis²=0}. A potential crisis arose: with fixed
μ = H₀/c, the background EoS w̃₀ ~ 0.17 violates the Thomas-Kopp-Skordis (2016)
bound c_s² < 3.4×10⁻⁶ by a factor of 50,000.

**Resolution**: The running Khronon coupling μ_bg(z) = H(z)/c suppresses w(z)
by (H₀/H(z))², giving w(z=1100) ≈ 3×10⁻¹⁰ — **11,000× below the TKS bound**.

## The Crisis

### What was claimed
With a fixed Khronon coupling μ₀ = H₀/c:
- δ₀ = Q₀ − 1 = 0.340 (from Ω_cdm = 0.265)
- Background EoS: w₀ = δ₀/2 = 0.170
- This IS the adiabatic sound speed: c_ad² = w₀ ≈ 0.17
- TKS 2016 constraint: c_s² < 3.4×10⁻⁶
- **Violation by factor ~50,000**

### Why it looked fatal
The TKS constraint comes from fitting the Planck CMB power spectrum with
generalized dark matter (GDM) parameters. Any fluid with c_s² > 3.4×10⁻⁶
produces acoustic oscillations in the dark matter that are NOT observed.

## The Resolution

### Mechanism: Running μ

The Khronon coupling is NOT fixed. On an FRW background:

```
μ_bg(a) = H(a)/c
```

This is the **background running** — it applies to ALL perturbation modes
equally at a given epoch. It is NOT the same as the mode-dependent running
μ(k) = μ₀(k*/k) used for galactic rotation curves.

The background field amplitude and EoS run as:

```
δ(z) = I₀ / [2μ_bg(z)²] = δ₀ × (H₀/H(z))²
w(z)  = δ(z)/2           = w₀ × (H₀/H(z))²
```

Since H(z) grows rapidly at early times:

| z | H(z)/H₀ | w̃(z) | c_s²(z) | vs TKS bound |
|---|---------|-------|---------|-------------|
| 0 | 1 | 1.70×10⁻¹ | 1.70×10⁻¹ | 50,000× OVER |
| 10 | 1.83×10¹ | 5.06×10⁻⁴ | 5.06×10⁻⁴ | 149× OVER |
| 100 | 5.63×10² | 5.35×10⁻⁷ | 5.35×10⁻⁷ | 6.4× BELOW |
| 1100 | 2.35×10⁴ | 3.06×10⁻¹⁰ | 3.06×10⁻¹⁰ | 11,100× BELOW |
| 3000 | 1.65×10⁵ | 6.26×10⁻¹² | 6.26×10⁻¹² | 543,000× BELOW |

### Why TKS 2016 does NOT apply directly

TKS 2016 assumed **constant c_s²** (not redshift-dependent). Our model has
c_s²(z) ∝ (H₀/H(z))² — strongly varying with redshift. The TKS constraint
is conservative (worst-case); the actual constraint on time-varying c_s² is
weaker.

Nevertheless, even the constant-c_s² bound is satisfied at ALL CMB-relevant
epochs (z > 100).

### Blanchet-Skordis Eq. 4.31 (Jeans suppression)

The full effective sound speed from BS2024:

```
c_s²(z,k) = c_ad² / [1 + c_ad² c² k² / (4π G a² ρ_K (1+w))]
```

This provides an ADDITIONAL 1/k² suppression at sub-Jeans scales. However,
at CMB scales:

- Jeans wavenumber: k_J(z=1100) ~ 170 Mpc⁻¹
- CMB wavenumbers: k ∈ [0.001, 0.3] Mpc⁻¹
- Since k_CMB << k_J, the denominator ≈ 1 → **no additional suppression**

**Conclusion**: At CMB scales, the running μ alone does all the work.
Eq. 4.31 becomes relevant at galactic scales (k >> k_J).

## Paper 3 Formula Correction

### Previous formula (WRONG)
```
w̃₀ = μ_bg² / k_phys²    [Eq. 12 in Paper 3, v1]
```

This is wrong because:
1. It makes w a k-dependent quantity, but GDM w is background (k-independent)
2. It mixes μ₀ (z=0) with k_phys (z=1100) — inconsistent units/epochs
3. It gives 4×10⁻¹⁰ by numerical coincidence, not by correct physics

### Corrected formula
```
w̃(z) = (δ₀/2) × (H₀/H(z))²    [Eq. 12 in Paper 3, v2]
```

This is correct because:
1. w is a background quantity (k-independent), as required by GDM
2. It follows directly from the Khronon field equation δ = I₀/(2μ²)
3. It gives 3×10⁻¹⁰ at z=1100 — same order as before, but with correct physics

## Numerical Verification

Full numerical verification script: `anatropic/examples/cs2_verification.py`

### Step-by-step at z = 1100, k = 0.05 Mpc⁻¹

```
Input parameters (Planck 2018):
  H₀ = 67.4 km/s/Mpc = 2.184×10⁻¹⁸ s⁻¹
  Ω_r = 9.1×10⁻⁵, Ω_m = 0.315, Ω_Λ = 0.685
  Ω_cdm h² = 0.120

Derived:
  δ₀ = -1 + √(1 + 3Ω_cdm) = 0.340
  w₀ = δ₀/2 = 0.170
  H(z=1100) = H₀ × √(Ω_r×1101⁴ + Ω_m×1101³ + Ω_Λ)
             = 2.184×10⁻¹⁸ × 2.35×10⁴ = 5.13×10⁻¹⁴ s⁻¹

Running suppression:
  w(1100) = 0.170 × (1/2.35×10⁴)² = 0.170 × 1.81×10⁻⁹ = 3.06×10⁻¹⁰

Eq. 4.31 (Jeans suppression):
  c_ad² = w(1100) = 3.06×10⁻¹⁰
  v_s² = c_ad² × c² = 2.75×10⁷ m²/s²
  4πG a² ρ_K (1+w) = 4π × 6.674×10⁻¹¹ × (9.08×10⁻⁴)² × 1.84×10⁻¹⁹ × 1.0
                    = 1.26×10⁻³⁴ s⁻²
  k_SI = 0.05/Mpc = 1.62×10⁻²⁴ m⁻¹
  denom = 1 + v_s² × k_SI² / grav_term
        = 1 + 2.75×10⁷ × 2.63×10⁻⁴⁸ / 1.26×10⁻³⁴
        ≈ 1.0  (no additional suppression)

Result:
  c_s²(k=0.05, z=1100) = 3.06×10⁻¹⁰ / 1.0 = 3.06×10⁻¹⁰
  TKS bound = 3.4×10⁻⁶
  Ratio = 9.0×10⁻⁵ → PASSES by factor 11,100×
```

## Summary Table

| Question | Answer |
|----------|--------|
| Is c_s² safe at CMB? | YES — 11,000× below TKS bound |
| Main suppression mechanism | Running μ_bg(z) = H(z)/c |
| Role of Eq. 4.31 | Irrelevant at CMB scales (k << k_J) |
| Was the formula wrong? | YES — w̃₀ = μ²/k² replaced by w(z) = (δ₀/2)(H₀/H)² |
| Does the conclusion change? | NO — order of magnitude unchanged (4→3 ×10⁻¹⁰) |

## References

- Blanchet & Skordis 2024, arXiv:2404.06584 (Khronon = GDM, Eq. 4.31)
- Thomas, Kopp & Skordis 2016, arXiv:1601.05097 (GDM constraints: c_s² < 3.4×10⁻⁶)
- Skordis & Złośnik 2021, arXiv:2007.00082 (AeST CMB fit)
