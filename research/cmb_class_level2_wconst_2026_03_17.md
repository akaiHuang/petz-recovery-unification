# CLASS Level 2: Constant-w Exclusion Table (cs2 = 1e-10 Regularization)

**Author**: Sheng-Kai Huang
**Date**: 2026-03-17
**Status**: Complete -- all 7 runs succeeded (no timeouts)
**Code**: gdm_class_public (Ilic, Kopp, Thomas & Skordis, arXiv:2004.09572)

---

## 1. Executive Summary

This analysis maps the CMB exclusion boundary for constant equation-of-state dark matter using the GDM framework. By scanning constant w from 0 to 0.14 with regularized sound speed cs2 = 1e-10, we establish:

1. **w = 0 (CDM baseline)**: SAFE. GDM reproduces LCDM to 0.002% -- confirms cs2 = 1e-10 regularization is invisible.
2. **w = 0.001**: MARGINAL (1.3 sigma per-ell). This is the edge of detectability.
3. **w >= 0.01**: EXCLUDED at high significance (>= 16 sigma).
4. **w = 0.14 (Khronon z=0 value)**: Catastrophically excluded at 281 sigma.

**The constant-w exclusion threshold is w < ~0.0015 (2-sigma), or w < ~0.001 (0.5-sigma, "safe").**

---

## 2. Technical Fix: cs2 = 1e-10 Regularization

Previous Level 2 runs with w > 0 and cs2 = 0 caused CLASS to hang indefinitely. The fix:

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| cs2_gdm | 1e-10 | Numerical regularization |
| TKS bound | 3.4e-6 | Observational limit from TKS (2016) |
| Safety margin | 4 orders of magnitude below TKS | Completely undetectable |
| Confirmed impact | < 0.002% on C_l (Level 1) | Negligible |

The w=0 run confirms this: GDM with w=0, cs2=1e-10 reproduces standard LCDM to 0.002% (max |dTT/TT|). This is well below cosmic variance at all multipoles.

---

## 3. Cosmological Parameters

All runs use Planck 2018 best-fit:

| Parameter | Value |
|-----------|-------|
| omega_b | 0.02237 |
| omega_cdm (LCDM) / omega_gdm (GDM) | 0.1200 |
| h | 0.6736 |
| tau_reio | 0.0544 |
| A_s | 2.1e-9 |
| n_s | 0.9649 |
| N_eff | 3.046 (N_ur=2.0328, N_ncdm=1, m_ncdm=0.06 eV) |

GDM runs: omega_cdm = 0, omega_gdm = 0.12, constant w and cs2 across 6 time bins.

---

## 4. Results: Constant-w Exclusion Table

| w | max\|dTT/TT\| | rms\|dTT/TT\| | max\|dEE/EE\| | N_sigma | Status | Runtime |
|------|----------------|----------------|----------------|---------|----------|---------|
| 0.000 | 0.002% | 0.001% | 0.017% | 0.0 | **SAFE** | 5.2s |
| 0.001 | 2.81% | 1.23% | 5.67% | 1.3 | **MARGINAL** | 5.2s |
| 0.010 | 33.9% | 13.0% | 75.9% | 15.9 | **EXCLUDED** | 5.2s |
| 0.050 | 134% | 45.1% | 186% | 67.2 | **EXCLUDED** | 4.8s |
| 0.100 | 375% | 87.1% | 440% | 187.5 | **EXCLUDED** | 4.4s |
| 0.140 | 591% | 214% | 744% | 281.2 | **EXCLUDED** | 4.1s |

### Scaling Law

The maximum TT deviation scales as:

```
max|dTT/TT| ~ w^{1.06}
```

Nearly linear in w. This makes physical sense: the constant w shifts the entire matter density history, and the CMB response is roughly proportional to the fractional change in dark matter density at recombination, which is ~3w * ln(1/a_rec) ~ 20w.

### Exclusion Thresholds (Interpolated)

| Criterion | N_sigma threshold | w_max |
|-----------|-------------------|-------|
| SAFE (below cosmic variance) | < 0.5 | w < 0.0010 |
| MARGINAL (detectable but not excluded) | < 2.0 | w < 0.0015 |
| EXCLUDED (definitively ruled out) | > 5.0 | w > 0.0034 |

---

## 5. Detailed Per-Run Analysis

### 5.1 w = 0.000 (CDM Baseline)

- max |dTT/TT| = 1.66e-05 (0.002%)
- Dominated by numerical noise from cs2 = 1e-10 regularization
- Confirms the regularization is invisible to any observation
- **Verdict: SAFE -- cs2 = 1e-10 has no observable effect**

### 5.2 w = 0.001

- max |dTT/TT| = 2.81% at l = 2175 (damping tail)
- Low-l ISW: 0.38%
- Acoustic peaks: 2.19%
- Damping tail: 2.81% (dominant)
- **Verdict: MARGINAL -- at the edge of Planck sensitivity (1.3 sigma per-ell)**
- Note: Per-ell cosmic variance is generous; a full chi2 analysis with Planck data would likely detect this.

### 5.3 w = 0.010

- max |dTT/TT| = 33.9% at l = 2202
- Low-l ISW: 3.2% (already substantial)
- Acoustic peaks: 21.1% (peak heights shifted)
- Damping tail: 33.9% (gross mismatch)
- **Verdict: EXCLUDED at 16 sigma**

### 5.4 w = 0.050

- max |dTT/TT| = 134% at l = 2500
- The TT spectrum is more than doubled at the damping tail
- All multipole regions grossly wrong
- **Verdict: EXCLUDED at 67 sigma**

### 5.5 w = 0.100

- max |dTT/TT| = 375% at l = 2500
- TT spectrum off by a factor of ~5 in the tail
- **Verdict: EXCLUDED at 188 sigma**

### 5.6 w = 0.140 (Khronon z=0 Value)

- max |dTT/TT| = 591% at l = 2262
- TT spectrum off by a factor of ~7 in the tail
- Even at low-l (ISW): 53% deviation
- **Verdict: EXCLUDED at 281 sigma -- catastrophic**

---

## 6. Multipole Structure of Deviations

The constant-w effect has a clear multipole structure:

1. **Low-l (l < 30, ISW effect)**: Deviation grows roughly as ~50*w. At w=0.14, this is ~53%. The late-time integrated Sachs-Wolfe effect is sensitive to the expansion history, which w modifies.

2. **Acoustic peaks (100 < l < 1000)**: Deviation grows as ~20*w. The peak heights and positions shift because the matter-radiation equality redshift changes, and the acoustic driving changes.

3. **Damping tail (l > 1000)**: Largest deviations, growing faster than linearly (~40*w for w=0.01 but much steeper for larger w). The Silk damping scale is modified because the DM density at recombination differs.

---

## 7. Implications for the Khronon Theory

### 7.1 The Key Question

The Khronon f(delta) correction gives an effective w_eff(a) that varies from ~0 at z >> 1 to ~0.14 at z = 0. How does this compare to the constant-w exclusion?

### 7.2 Constant-w is a WORST CASE for Late-Time Effects

A constant w=0.14 means w=0.14 at ALL epochs, including recombination. This changes the DM density at z=1100 by a factor of ~(1100)^{3*0.14} ~ 7, which is catastrophic.

The Khronon w_eff(a) is ~0 at z=1100, so the effect on recombination physics is MUCH smaller. The constant-w table is therefore an **upper bound** for the damping tail and acoustic peak effects.

### 7.3 But the Khronon Also Has a Density Normalization Problem

From the previous Level 2 analysis (11-bin w(a), Run B):
- At z=1100, the Khronon has rho_K = 1.287 * rho_CDM (with delta_0 = 0.34)
- This 28.7% density excess at recombination is the DOMINANT effect
- It produces 33% TT deviation even though w_eff(z=1100) ~ 0

The density normalization problem is NOT captured by the constant-w table. It comes from normalizing omega_gdm = 0.12 at z=0 while f(delta) > 1 at z > 0.

### 7.4 Combined Interpretation

| Effect | Constant-w captures? | Magnitude | Notes |
|--------|---------------------|-----------|-------|
| Late-time ISW | Partially (overestimates) | ~0.5-5% | Depends on w_eff(z<2) |
| Acoustic peak shifts | No (overestimates for peaks) | ~20% from density | Dominated by density normalization |
| Damping tail | No (overestimates) | ~30% from density | Dominated by density normalization |
| Background density shift | Not directly | ~29% at z=1100 | This is the real problem |

### 7.5 What the Constant-w Table DOES Tell Us

The table establishes a clean, model-independent exclusion:

**Any dark matter component with constant w > 0.001 and the same omega as CDM is at the edge of CMB exclusion. Constant w > 0.01 is definitively excluded.**

For the Khronon theory, the relevant comparison is:
- The **average** w over cosmic history (weighted by the CMB kernel) is what matters
- If the time-averaged effective w is < 0.001, the CMB is safe
- If it is > 0.01, it is excluded

From the 11-bin analysis:
- The average w weighted by the recombination kernel is essentially 0 (w_eff ~ 0 for z > 100)
- But the density normalization effect is equivalent to w_avg ~ 0.03 in terms of CMB impact
- This is why the 11-bin runs showed 33% deviation: the density normalization is equivalent to a constant w ~ 0.03

### 7.6 Conclusion for the Theory

The constant-w exclusion table, combined with the previous 11-bin w(a) analysis, establishes:

1. **The Khronon c_s^2 = 0 perturbation prediction is SAFE** (Level 1 result, confirmed).
2. **The Khronon background f(delta) correction with delta_0 = 0.34 is EXCLUDED** (Level 2 result, confirmed).
3. **Any rescue requires delta_0 < 0.003** (from the delta_0 scan) OR a different cosmological f(delta).
4. **The constant-w exclusion threshold is w < 0.001 (safe) or w < 0.0015 (marginal).**
5. **The regularization cs2 = 1e-10 works perfectly** -- no timeouts, invisible to observations.

---

## 8. Technical Details

### 8.1 GDM Bin Configuration

All GDM runs use 6 time bins with identical w and cs2 in each bin:
```
time_values_gdm = 1e-05, 0.0001, 0.001, 0.01, 0.1
w_values_gdm = w, w, w, w, w, w     (6 values for 6 bins)
cs2_values_gdm = 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10
cv2_values_gdm = 0, 0, 0, 0, 0, 0
```

### 8.2 Files Generated

- INI files: `/Users/akaihuangm1/Desktop/github/gdm_class_public/level2_wconst_*.ini`
- C_l output: `/Users/akaihuangm1/Desktop/github/gdm_class_public/output/level2_wconst_*_cl.dat`
- Analysis script: `/Users/akaihuangm1/Desktop/github/gdm_class_public/level2_wconst_scan.py`
- Summary: `/Users/akaihuangm1/Desktop/github/gdm_class_public/output/level2_wconst_summary.txt`

### 8.3 Runtime

All 7 CLASS runs completed in ~34 seconds total. No timeouts. The cs2 = 1e-10 regularization completely eliminates the numerical instability that caused previous runs to hang.

---

## References

- Ilic, Kopp, Thomas & Skordis 2020, arXiv:2004.09572 (gdm_class_public)
- Thomas, Kopp & Skordis 2016, arXiv:1601.05097 (TKS bound: c_s^2 < 3.4e-6)
- Blanchet & Skordis 2024, arXiv:2404.06584 (Khronon dark matter)
- Planck 2018, arXiv:1807.06209 (cosmological parameters)

---

*Generated: 2026-03-17. All CLASS runs completed successfully with cs2 = 1e-10 regularization.*
