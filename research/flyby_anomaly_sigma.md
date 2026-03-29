# Flyby Anomaly vs Sigma Framework: Honest Calculation

**Date**: 2026-03-29
**Result**: NEGATIVE — Sigma framework cannot explain the flyby anomaly
**Status of anomaly itself**: Largely resolved by conventional effects (post-2013)

---

## 1. Background

The flyby anomaly (Anderson et al. 2008, arXiv:0809.4351): spacecraft performing Earth gravity assists showed asymptotic velocity changes differing from GR predictions by ~1-13 mm/s.

**Empirical formula** (Anderson et al.):

$$\frac{\Delta v}{v} = \frac{2\omega_E R_E}{c}(\cos\delta_{\rm in} - \cos\delta_{\rm out})$$

Known flyby data:

| Flyby | Perigee alt (km) | v_inf (km/s) | Δv_obs (mm/s) |
|-------|-----------------|---------------|----------------|
| Galileo I (1990) | 960 | 8.949 | +3.92 |
| NEAR (1998) | 539 | 6.851 | +13.46 |
| Cassini (1999) | 1175 | 16.01 | -2.0 |
| Rosetta I (2005) | 1956 | 3.863 | +1.80 |
| Messenger (2005) | 2347 | 4.056 | +0.02 |
| Rosetta II (2007) | 5322 | 9.36 | 0.0 |
| **Juno (2013)** | 559 | — | **0 ± 0.1** |

---

## 2. Earth Parameters

| Quantity | Value |
|----------|-------|
| M_E | 5.9722 × 10²⁴ kg |
| R_E | 6.371 × 10⁶ m |
| ω_E | 7.2921 × 10⁻⁵ rad/s |
| I_E | 8.017 × 10³⁷ kg·m² (factor 0.3307) |
| J_E | 5.846 × 10³³ kg·m²/s |
| **r_s = 2GM_E/c²** | **8.870 × 10⁻³ m = 8.87 mm** |

---

## 3. Task 1: Static Sigma Correction

### Metric difference

Exponential metric: g₀₀ = -exp(-r_s/r)

Schwarzschild: g₀₀ = -(1 - r_s/r)

Expansion:
```
exp(-r_s/r) = 1 - (r_s/r) + (1/2)(r_s/r)² - (1/6)(r_s/r)³ + ...
Schwarzschild = 1 - (r_s/r)
```

**Leading difference**: Δg₀₀ = (1/2)(r_s/r)² [at O((r_s/r)²), NOT O((r_s/r)³)]

> **Correction to the problem statement**: The leading difference between the exponential and Schwarzschild metrics is at second order in r_s/r, not third order. This makes the effect slightly larger than initially estimated, but still hopelessly small.

### Numerical results

At Galileo I perigee (r = 7.331 × 10⁶ m):
```
r_s/r = 8.87e-3 / 7.33e6 = 1.21 × 10⁻⁹
(r_s/r)² = 1.46 × 10⁻¹⁸
Δg₀₀ = (1/2)(r_s/r)² = 7.32 × 10⁻¹⁹
Δv = v × Δg₀₀/2 = 8949 × 3.66e-19 = 3.28 × 10⁻¹⁵ m/s
```

**Result: Δv_static ~ 10⁻¹⁵ m/s vs observed 3.92 × 10⁻³ m/s**

**DEFICIT: 12 orders of magnitude** (corrected from 18 in preliminary estimate; still catastrophically small)

Full table:

| Flyby | r_s/r | (r_s/r)² | Δv_static (m/s) | Δv_obs (mm/s) |
|-------|-------|----------|-----------------|----------------|
| Galileo I | 1.21e-9 | 1.46e-18 | 3.28e-15 | 3.92 |
| NEAR | 1.28e-9 | 1.65e-18 | 2.82e-15 | 13.46 |
| Cassini | 1.18e-9 | 1.38e-18 | 5.53e-15 | -2.0 |
| Rosetta I | 1.07e-9 | 1.13e-18 | 1.10e-15 | 1.80 |
| Messenger | 1.02e-9 | 1.04e-18 | 1.05e-15 | 0.02 |
| Rosetta II | 7.59e-10 | 5.75e-19 | 1.35e-15 | 0.0 |

---

## 4. Task 2: Frame-Dragging (Lense-Thirring)

This is NOT a Sigma effect, but sets the scale for rotational gravitomagnetic corrections.

```
v_LT = 2GJ/(c²r²)
```

At Galileo perigee:
```
v_LT = 2 × 6.674e-11 × 5.846e33 / (9e16 × 5.374e13)
     = 1.62 × 10⁻⁷ m/s = 1.62 × 10⁻⁴ mm/s
```

**4 orders of magnitude too small** for the observed anomaly.

Note: the Anderson empirical prefactor 2ω_E R_E/c = 3.10 × 10⁻⁶ is much larger than GJ/(c²R²) ~ 10⁻⁷/v. The empirical formula is kinematic (involving Earth's surface velocity), not gravitational.

---

## 5. Task 3: Exponential Rotating Metric vs Kerr

If the exponential metric generalizes to a rotating solution, the off-diagonal term would be:
```
g₀φ^exp = g₀φ^Kerr × [1 + O(r_s/r)]
```

The correction to frame-dragging:
```
Δv_rot = (r_s/r) × v_LT
       = 1.21e-9 × 1.62e-7
       = 1.95 × 10⁻¹⁶ m/s
```

**13 orders of magnitude too small.**

Even optimistic scenarios (second-order corrections, coupling terms) cannot bridge this gap:

| Effect | Magnitude (m/s) | Orders below observed |
|--------|-----------------|----------------------|
| Static Σ correction | 10⁻¹⁵ | 12 |
| Lense-Thirring (GR) | 10⁻⁷ | 4 |
| Σ × Lense-Thirring | 10⁻¹⁶ | 13 |
| (r_s/r)² × LT | 10⁻²⁵ | 22 |

---

## 6. Task 4: Is the Flyby Anomaly Still Anomalous?

**Short answer: Largely NO.**

Key developments post-Anderson (2008):

1. **Juno Earth flyby (2013)**: Δv = 0 ± 0.1 mm/s — NO anomaly detected (Thompson et al. 2014). This was a high-precision test with modern tracking.

2. **OSIRIS-REx Earth flyby (2017)**: Consistent with zero anomaly.

3. **Anderson formula fails**: Predicts Δv ~ 7 mm/s for Juno, but observed Δv = 0. The empirical formula does not hold for modern flybys.

4. **Conventional explanations identified**:
   - Thermal radiation recoil from spacecraft asymmetries
   - Earth albedo radiation pressure
   - Atmospheric drag at low perigee (especially NEAR at 539 km)
   - Tracking station geometry systematics (incomplete Doppler coverage)
   - Improved orbit determination software eliminates some residuals

5. **Community consensus (circa 2020-2025)**: The flyby anomaly is either resolved or too marginal to constitute evidence for new physics. No new flybys have reproduced the effect with modern tracking capabilities.

---

## 7. Task 5: Fundamental Impossibility

For the flyby anomaly to be explained by a metric correction:
```
Required: Δv/v ~ 10⁻⁶
Available from (r_s/r)^n at r ~ R_E:
  r_s/r ~ 10⁻⁹

  n=1: 10⁻⁹   (3 orders too small — this IS Newtonian gravity itself)
  n=2: 10⁻¹⁸  (12 orders too small)
  n=3: 10⁻²⁷  (21 orders too small)
```

To achieve 10⁻⁶ from (r_s/r)^n, we need n = log(10⁻⁶)/log(10⁻⁹) = 0.67.

**A fractional power of r_s/r is not physical.** No Taylor expansion of any metric produces fractional powers.

The observed scale 2ω_E R_E/c ~ 10⁻⁶ corresponds to the ratio of Earth's equatorial velocity to c. This is:
- NOT proportional to any power of r_s/r
- A purely kinematic quantity
- Unrelated to the PPN expansion or metric corrections

**No metric theory** — not just Sigma, but ANY theory based on post-Newtonian corrections — can produce an effect at this scale for Earth. The effect would need to be proportional to v_rot/c, which is not how gravitational metric corrections work.

---

## 8. Conclusions

### Honest Assessment

| Question | Answer |
|----------|--------|
| Can static Σ explain the flyby anomaly? | **NO** (12 orders too small) |
| Can rotating Σ explain it? | **NO** (13 orders too small for the correction) |
| Can ANY metric correction explain it? | **NO** (would need fractional power n=0.67) |
| Is the anomaly even real? | **Probably NOT** (Juno 2013 null result; conventional explanations) |
| Does this hurt the Σ framework? | **NO** — this is actually good news |

### Why This is Good News

The Sigma framework's exponential metric corrections are too small at Earth scales to produce any detectable flyby effect. This means:

1. **No conflict**: The framework does not predict an anomaly that doesn't exist.
2. **Consistency**: The exponential metric agrees with Schwarzschild to extraordinary precision at Earth-scale fields (r_s/r ~ 10⁻⁹), as it should.
3. **Testability preserved**: The framework's distinctive predictions (Papers 2-3) live in the strong-field regime (r_s/r not too small) and the cosmological regime (Sigma integrated over Hubble volume), NOT in Earth-scale weak fields.

### Lesson

The flyby anomaly is a cautionary tale: not every unexplained observation requires new physics. When the anomaly was first reported, dozens of papers proposed exotic explanations (dark matter halos, modified gravity, etc.). The mundane explanation — incomplete thermal and tracking modeling — turned out to be correct. The Sigma framework correctly predicts zero anomaly at this scale, consistent with modern observations.

---

## Appendix: Calculation Cross-Checks

**Check 1**: r_s for Earth
```
r_s = 2 × 6.674e-11 × 5.972e24 / (2.998e8)² = 8.870e-3 m ✓
```

**Check 2**: Anderson prefactor
```
2ω_E R_E/c = 2 × 7.292e-5 × 6.371e6 / 2.998e8 = 3.10e-6 ✓
```

**Check 3**: Lense-Thirring at R_E
```
2GJ/(c²R_E²) = 2 × 6.674e-11 × 5.846e33 / (8.988e16 × 4.059e13)
             = 7.801e23 / 3.649e30 = 2.14e-7 m/s ✓
```

**Check 4**: (r_s/R_E)² sanity
```
(8.87e-3 / 6.37e6)² = (1.39e-9)² = 1.94e-18 ✓
```

All calculations verified independently.
