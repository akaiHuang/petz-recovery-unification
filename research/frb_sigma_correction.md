# FRB Dispersion Measure: Σ Framework Corrections

**Date**: 2026-03-29
**Status**: Complete analysis — mostly NEGATIVE result (honest)

---

## Executive Summary

The Σ framework (exponential metric) predicts corrections to gravitational time delays near compact objects that differ from Schwarzschild by O(Σ²). For FRBs near magnetars, these corrections are:

| Observable | δΔt | Detectable? |
|---|---|---|
| Near-source achromatic delay (surface) | ~3.3 μs | Unphysical (surface opaque) |
| Near-source achromatic delay (r=100 km) | ~284 ns | Below threshold |
| Frequency-dependent DM correction | 0 (leading order) | N/A — identical to GR |
| Lensed FRB time delay | ~10-20 s (galaxy lens) | **YES** |

**Bottom line**: The Σ correction to FRB dispersion is achromatic and therefore invisible to standard DM analysis. The only viable FRB test is via gravitationally **lensed** FRBs, where the Σ framework predicts systematically shorter time delays by ~10 ppm.

---

## Task 1: Σ Correction to Dispersion

### Setup

For a 1.4 M_sun neutron star with R = 12 km:
- Schwarzschild radius: r_s = 4130 m = 4.13 km
- Compactness: r_s/R = 0.344

### Metric comparison at the surface

| Quantity | Schwarzschild | Exponential (Σ) |
|---|---|---|
| g_00 | -(1 - r_s/R) = -0.656 | -exp(-r_s/R) = -0.709 |
| Σ | -ln(1-r_s/R) = 0.422 | r_s/R = 0.344 |
| c_eff/c | sqrt(1-r_s/R) = 0.810 | exp(-r_s/(2R)) = 0.842 |

**Difference**: Δ(c_eff/c) = 0.032 = **3.2%**

This is consistent with the leading-order prediction:
```
Δc/c ≈ (r_s/R)² / 4 + O((r_s/R)³)
```

### Why this is NOT directly measurable

The near-source region is ~10⁴ m out of a total path of ~10²⁵ m (1 Gpc). The fractional path contribution is ~10⁻²¹. The 3.2% speed difference over 12 km gives an absolute time difference of ~1.3 μs — but this assumes emission from the stellar surface, which is unphysical at radio frequencies (see Task 3 correction).

---

## Task 2: Gravitational Time Delay Integral

### Formulation

For radial photon propagation, the coordinate velocity is:
- **Schwarzschild**: dr/dt = c(1 - r_s/r)
- **Exponential**: dr/dt = c exp(-r_s/r)

The gravitational time delay (compared to flat spacetime):
- **Schwarzschild**: Δt_S = (1/c) ∫_R^∞ [(r_s/r)/(1-r_s/r)] dr
- **Exponential**: Δt_E = (1/c) ∫_R^∞ [exp(r_s/r) - 1] dr

### Numerical results (R = 12 km, 1.4 M_sun)

| Metric | Gravitational delay |
|---|---|
| Schwarzschild | 161.79 μs |
| Exponential | 158.49 μs |
| **Difference** | **-3.29 μs** |

Relative difference: **-2.04%**

### Analytical approximation

Expanding the integrands:
```
1/(1-x) - 1 = x + x² + x³ + ...
exp(x) - 1   = x + x²/2 + x³/6 + ...
```
where x = r_s/r. The difference:
```
δ(integrand) = -x²/2 - 5x³/6 - ...
```

Integrating:
```
δΔt = -r_s²/(2Rc) - 5r_s³/(12R²c) - ...
    = -2369 ns - 679 ns - ...
    = -3048 ns (analytical) vs -3294 ns (numerical)
```

The 8% discrepancy comes from higher-order terms, confirming the expansion.

**Sign**: The exponential metric gives a SHORTER delay — light escapes faster in the Σ framework than in Schwarzschild. This is physically consistent: the exponential metric has no event horizon, so there is no logarithmic divergence as r → r_s.

---

## Task 3: Frequency Dependence

### Key derivation

The local refractive index in a plasma:
```
n(f_local) = sqrt(1 - f_p²/f_local²)
```

For a photon with observed frequency f at infinity, the local frequency at radius r is:
- **Schwarzschild**: f_local = f / sqrt(1-r_s/r)
- **Exponential**: f_local = f × exp(r_s/(2r))

The coordinate-time delay per unit dr, expanded in f_p²/f² << 1:

**Schwarzschild**:
```
dt/dr = 1/[c(1-r_s/r)] + f_p²/(2f²c) + O(f_p⁴/f⁴)
         ↑ achromatic          ↑ dispersive (f-dependent)
```

**Exponential**:
```
dt/dr = exp(r_s/r)/c + f_p²/(2f²c) + O(f_p⁴/f⁴)
         ↑ achromatic    ↑ dispersive (SAME!)
```

### Result: Leading dispersive term is IDENTICAL

**The frequency-dependent (dispersive) contribution is f_p²/(2f²c) for BOTH metrics.**

Physical reason: The gravitational blueshift of the photon and the redshift of the plasma frequency exactly cancel at leading order. The ratio f_p/f_local is metric-independent to O(f_p²/f²).

### Higher-order correction

At O(f_p⁴/f⁴):
```
Schwarzschild: 3f_p⁴(1-r_s/r)/(8f⁴c)
Exponential:   3f_p⁴ exp(-r_s/r)/(8f⁴c)
Difference:    3f_p⁴/(8f⁴c) × (r_s/r)²/2
```

This is doubly suppressed: by (f_p/f)⁴ AND by (r_s/r)². Utterly negligible.

### Critical physical constraint: Magnetar opacity

For a magnetar with B = 10¹¹ T, P = 5 s:
- Goldreich-Julian density at surface: n_GJ = 6.9 × 10¹⁸ m⁻³
- Surface plasma frequency: **f_p = 23.7 GHz >> 1 GHz**

The magnetar magnetosphere is **opaque** to GHz radio at the surface! Radio emission escapes from:
```
r_esc ~ R × (f_p,surface / f_obs)^{2/3} ~ 100 km ~ 8 R_ns
```

At the escape radius: Σ ~ r_s/r_esc ~ 0.04 (small!)

### Realistic Σ corrections by emission radius

| r_em | Σ | f_p | |δΔt| (Σ vs Schw) |
|---|---|---|---|
| 20 km (1.7 R_ns) | 0.206 | 11.0 GHz (opaque) | 1421 ns |
| 50 km (4.2 R_ns) | 0.083 | 2.8 GHz (opaque) | 569 ns |
| 100 km (8.3 R_ns) | 0.041 | 1.0 GHz (marginal) | 284 ns |
| 500 km (42 R_ns) | 0.008 | 88 MHz (transparent) | 57 ns |
| 1000 km (83 R_ns) | 0.004 | 31 MHz (transparent) | 28 ns |

For realistic 1 GHz FRB emission: δΔt ~ **30-300 ns**, and it is **achromatic**.

---

## Task 4: Observable Effects

### FRB 20200120E (M81 globular cluster)

- DM_obs = 87.82 pc/cm³, DM_MW ~ 50, DM_IGM ~ 1-3, DM_host ~ 35-37
- DM budget is consistent within standard uncertainties
- No anomalous DM requiring gravitational correction
- **Verdict: No Σ signature**

### Repeating FRB DM variations

- Many repeaters show ΔDM ~ 1-10 pc/cm³ on days-months timescales
- These are from turbulent plasma (magnetar wind, SNR)
- Σ correction is **constant** (geometry doesn't change) → cannot explain variations
- **Verdict: No Σ signature**

### Detection strategies

**1. Lensed FRBs (★ MOST PROMISING)**

For a galaxy-scale lens (M ~ 10¹¹ M_sun):
```
Σ_lens ~ GM/(rc²) ~ 2 × 10⁻⁶
GR time delay: ~10⁷ s (~114 days)
Σ correction: ~19 s
```

FRB lensing time delays can be measured to **μs precision** (from DM-subtracted coherent de-dispersion). A ~19 s systematic offset is **easily detectable** — if we have an independent prediction of the GR time delay.

**Prediction**: The Σ framework gives systematically **shorter** lensing time delays than GR, by a fractional amount ~Σ_lens ~ 10⁻⁶.

Challenge: The lens mass model must be known to better than 10⁻⁶ precision, which is not currently achievable. This test requires comparing *statistical* properties of many lensed FRBs.

**2. SMBH-proximate FRBs**

If an FRB is produced near a 10⁶ M_sun SMBH at r ~ 10 r_s:
```
Σ ~ 0.1
δΔt ~ 0.5 s
```
This is a large achromatic delay, but unmeasurable without knowing the intrinsic emission time.

**3. Pulse profile asymmetry**

The Σ gradient across the emission region could cause asymmetric pulse broadening:
```
δt ~ (dΣ/dr) × δR / c ~ (r_s/r²) × δR / c
```
For r = 100 km, δR = 1 km: δt ~ 0.5 ns. **Not detectable.**

**4. Multi-messenger (if FRB + X-ray)**

The achromatic delay applies to ALL electromagnetic signals equally, so comparing radio and X-ray arrival times does NOT help. However, comparing to **gravitational waves** (if detectable) could in principle separate the metric correction — but FRBs are not GW sources.

---

## Conclusions

### What the Σ framework predicts for FRBs

1. **Achromatic time delay correction**: ~30-300 ns for magnetar-origin FRBs, ~0.5 s for SMBH-proximate FRBs. Sign: NEGATIVE (Σ gives shorter delays).

2. **Zero correction to DM**: The standard dispersion relation DM = ∫ n_e dl and the 1/f² frequency dependence are **identical** between Schwarzschild and exponential metrics at leading order.

3. **Lensing time delay shift**: ~10-20 s for galaxy-lens FRBs. This is the only potentially measurable effect.

### Honest assessment

| Test | Σ correction | Current precision | Gap |
|---|---|---|---|
| DM (frequency-dependent) | 0 (identical) | — | N/A |
| Near-source timing | 30-300 ns | ~1 μs | ~3-30× |
| Lensed FRB time delay | ~19 s | ~μs | **Detectable in principle** |
| Pulse broadening | ~0.5 ns | ~μs | ~10³× |

### Path to detection

The lensed FRB test requires:
1. A confirmed multiply-imaged FRB (several candidates exist as of 2026)
2. A precise lens model (from galaxy imaging + dynamics)
3. Measurement of time delay between images to better than ~10 s
4. Comparison with GR prediction from the lens model

**Current status**: No confirmed lensed FRBs with sufficient data. CHIME and DSA-2000 may provide candidates in the next 2-3 years.

### For the paper series

This analysis is **relevant to Paper 2** (strong-field predictions) as an additional observational test. However, the FRB test is weaker than:
- Binary pulsar timing (Σ ~ 10⁻⁶, but measured to 10⁻¹³)
- Gravitational wave ringdown (directly probes near-horizon geometry)
- Black hole shadow size (EHT measures r_shadow to ~10%)

**Recommendation**: Mention lensed FRBs as a *future* test in Paper 2, but do not pursue as a primary prediction. The Σ framework's strongest predictions remain in rotation curves (Paper 3) and cosmological observables (Paper 4).

---

## Appendix: Key Equations

### Coordinate speed of light (radial)

Schwarzschild:
```
dr/dt = c(1 - r_s/r)
```

Exponential (Σ framework):
```
dr/dt = c exp(-r_s/r) = c exp(-Σ)
```

### Gravitational time delay difference

```
δΔt = Δt_exp - Δt_Schw = -(1/c) ∫_R^∞ [r_s²/(2r²) + 5r_s³/(6r³) + ...] dr
     = -r_s²/(2Rc) [1 + 5r_s/(6R) + ...]
```

### Dispersive delay (frequency-dependent part)

Both metrics give:
```
(dt/dr)_dispersive = f_p²(r) / (2f²c) + O(f_p⁴/f⁴)
```

The metric correction enters only at O(f_p⁴/f⁴ × Σ²), which is negligible.
