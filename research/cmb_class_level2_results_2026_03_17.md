# CLASS Boltzmann Solver Results: Level 2 -- f(delta) Background Correction

**Author**: Sheng-Kai Huang
**Date**: 2026-03-17
**Status**: Complete
**Code**: gdm_class_public (Ilic, Kopp, Thomas & Skordis, arXiv:2004.09572)

---

## 1. Executive Summary

Level 2 tests the Khronon f(delta) background energy density correction using the GDM binned w(a) approach. The key physics:

- rho_K(a) = rho_CDM(a) x f(delta(a)) / f(delta_0)
- f(delta) = (4 + 3*delta) / 4
- delta(a) = delta_0 x (H_0/H(a))^2 / a^3 (from running mu = H(a)/c)

### Critical Findings

1. **At high redshift, the Khronon has MORE energy density than CDM** (when normalized to match at z=0)
   - At z=1100: rho_K = 1.2871 x rho_CDM (+28.7%)
   - Peak ratio: 1.4413 at z ~ 11

2. **The effective w(a) is positive at late times** (~0.14 at z=0), meaning the Khronon DM redshifts FASTER than a^{-3} at late times.

3. **The w(a) is negligible at early times** (<10^{-3} at z > 100), so the early universe is CDM-like.

---

## 2. Physics of the f(delta) Correction

### 2.1 delta(a) Evolution

| Epoch | z | delta | f(delta) | rho_K/rho_CDM | w_eff |
|-------|---|-------|----------|---------------|-------|
| Radiation | 10000 | 0.2769 | 1.2077 | 0.962297 | -4.267e-02 |
| Mat-rad eq | 3400 | 0.5443 | 1.4082 | 1.122102 | -4.809e-02 |
| Recombination | 1100 | 0.8205 | 1.6153 | 1.287123 | -3.084e-02 |
| z=500 | 500 | 0.9456 | 1.7092 | 1.361906 | -1.761e-02 |
| z=100 | 100 | 1.0526 | 1.7895 | 1.425867 | -4.202e-03 |
| z=10 | 10 | 1.0784 | 1.8088 | 1.441252 | +2.558e-04 |
| z=3 | 3 | 1.0466 | 1.7850 | 1.422278 | +1.435e-02 |
| z=1 | 1 | 0.8506 | 1.6379 | 1.305139 | +8.351e-02 |
| z=0.5 | 0.5 | 0.6574 | 1.4930 | 1.189660 | +1.298e-01 |
| Today | 0 | 0.3400 | 1.2550 | 1.000000 | +1.394e-01 |

### 2.2 Key Insight: Non-Monotonic delta(a)

delta(a) is NOT monotonically decreasing with redshift:
- Radiation era (z >> 3000): delta ~ delta_0 * a / Omega_r -> very small
- Matter era (z ~ 10-1000): delta ~ delta_0 / Omega_m ~ 1.08 (LARGER than delta_0!)
- Lambda era (z < 1): delta decreases back to delta_0 at z=0

Peak: delta_max = 1.0785 at z ~ 11

### 2.3 The w(a) Approach

The GDM code parameterizes dark matter via w(a), c_s^2(a), c_vis^2(a). We set c_s^2 = c_vis^2 = 0 (omega=0 dispersion, Level 1 result) and model the f(delta) correction entirely through w(a).

The effective EOS is:
w_eff(a) = -(1/3) * d ln f(delta(a)) / d ln a

This is ~0 at high z and grows to ~0.14 at z=0.

---

## 3. CLASS Run Configurations

| Run | Description | omega_gdm | omega_cdm | Bins | Key Feature |
|-----|-------------|-----------|-----------|------|-------------|
| A | 6-bin w(a), delta_0=0.34 | 0.1200 | 0.0000 | 6 | ... |
| B | 11-bin w(a), delta_0=0.34 | 0.1200 | 0.0000 | 11 | ... |
| C | 6-bin w(a), delta_0=0.10 | 0.1200 | 0.0000 | 6 | ... |
| D | constant w=0.001 (sanity) | 0.1200 | 0.0000 | 6 | ... |
| E | 11-bin w(a), omega_gdm=0.0932 (renormalized) | 0.0932 | 0.0000 | 11 | ... |
| F | Two-comp: omega_cdm=0.0932 + omega_gdm=0.0268 | 0.0268 | 0.0932 | 11 | ... |
| G | constant w=0.01 (sanity) | 0.1200 | 0.0000 | 6 | ... |

### Run descriptions:
- **A**: 6-bin w(a) from f(delta) with delta_0=0.34. Late-time only (w=0 for a<0.01).
- **B**: 11-bin w(a) from f(delta) with delta_0=0.34. Better late-time resolution.
- **C**: 6-bin w(a) from f(delta) with delta_0=0.10. Conservative correction.
- **D**: Constant w=0.001 everywhere. Sanity check.
- **E**: 11-bin w(a) with REDUCED omega_gdm=0.0932 to match omega_cdm=0.12 at recombination.
- **F**: Two-component: omega_cdm=0.0932 (standard CDM) + omega_gdm=0.0268 (f(delta) excess).
- **G**: Constant w=0.01 everywhere. Larger sanity check.

---

## 4. Results

| Run | max|dTT/TT| | rms|dTT/TT| | max N_sigma | Status |
|-----|-------------|-------------|-------------|--------|
| A | 1.8631e-02 (1.86%) | 1.0920e-02 (1.09%) | 0.6 | **MARGINAL** |
| B | 3.3047e-01 (33.05%) | 1.6969e-01 (16.97%) | 12.8 | **EXCLUDED** |
| C | 7.8541e-03 (0.79%) | 4.6075e-03 (0.46%) | 0.3 | **MARGINAL** |
| D | 2.8119e-02 (2.81%) | 1.2263e-02 (1.23%) | 1.1 | **MARGINAL** |
| E | 5.8426e-02 (5.84%) | 3.3795e-02 (3.38%) | 2.3 | **EXCLUDED** |
| F | 1.0358e-01 (10.36%) | 4.9159e-02 (4.92%) | 4.0 | **EXCLUDED** |
| G | 3.3913e-01 (33.91%) | 1.3038e-01 (13.04%) | 13.3 | **EXCLUDED** |

### Run A: 6-bin w(a), delta_0=0.34
- Max |Delta C_l^TT / C_l^TT| = 1.8631e-02 (1.863%)
- RMS |Delta C_l^TT / C_l^TT| = 1.0920e-02 (1.092%)
- Max |Delta C_l^EE / C_l^EE| = 3.0077e-02 (3.008%)
- Max significance = 0.6 sigma (per-ell cosmic variance)

### Run B: 11-bin w(a), delta_0=0.34
- Max |Delta C_l^TT / C_l^TT| = 3.3047e-01 (33.047%)
- RMS |Delta C_l^TT / C_l^TT| = 1.6969e-01 (16.969%)
- Max |Delta C_l^EE / C_l^EE| = 5.3015e-01 (53.015%)
- Max significance = 12.8 sigma (per-ell cosmic variance)

### Run C: 6-bin w(a), delta_0=0.10
- Max |Delta C_l^TT / C_l^TT| = 7.8541e-03 (0.785%)
- RMS |Delta C_l^TT / C_l^TT| = 4.6075e-03 (0.461%)
- Max |Delta C_l^EE / C_l^EE| = 1.2749e-02 (1.275%)
- Max significance = 0.3 sigma (per-ell cosmic variance)

### Run D: constant w=0.001 (sanity)
- Max |Delta C_l^TT / C_l^TT| = 2.8119e-02 (2.812%)
- RMS |Delta C_l^TT / C_l^TT| = 1.2263e-02 (1.226%)
- Max |Delta C_l^EE / C_l^EE| = 5.6741e-02 (5.674%)
- Max significance = 1.1 sigma (per-ell cosmic variance)

### Run E: 11-bin w(a), omega_gdm=0.0932 (renormalized)
- Max |Delta C_l^TT / C_l^TT| = 5.8426e-02 (5.843%)
- RMS |Delta C_l^TT / C_l^TT| = 3.3795e-02 (3.379%)
- Max |Delta C_l^EE / C_l^EE| = 1.0674e-01 (10.674%)
- Max significance = 2.3 sigma (per-ell cosmic variance)

### Run F: Two-comp: omega_cdm=0.0932 + omega_gdm=0.0268
- Max |Delta C_l^TT / C_l^TT| = 1.0358e-01 (10.358%)
- RMS |Delta C_l^TT / C_l^TT| = 4.9159e-02 (4.916%)
- Max |Delta C_l^EE / C_l^EE| = 1.8541e-01 (18.541%)
- Max significance = 4.0 sigma (per-ell cosmic variance)

### Run G: constant w=0.01 (sanity)
- Max |Delta C_l^TT / C_l^TT| = 3.3913e-01 (33.913%)
- RMS |Delta C_l^TT / C_l^TT| = 1.3038e-01 (13.038%)
- Max |Delta C_l^EE / C_l^EE| = 7.5877e-01 (75.877%)
- Max significance = 13.3 sigma (per-ell cosmic variance)


---

## 5. Interpretation

### 5.1 Critical Finding: 6-bin vs 11-bin Discrepancy

Run A (6-bin) shows only 1.86% TT deviation, while Run B (11-bin, same physics) shows 33%. **The 6-bin version drastically under-resolves the late-time w(a) evolution.** The 6-bin structure averages w over the entire range a=0.1 to a=1.0, getting w_avg ~ 0.005, while the 11-bin version correctly captures w growing from ~0.002 at a=0.1 to ~0.146 at a=0.9. **Run B is the physically correct result.**

### 5.2 The f(delta) Effect on CMB

With delta_0 = 0.34, the Khronon f(delta) correction produces:
1. **+28.7% more DM at recombination** (z=1100): shifts acoustic peak positions and heights
2. **Rapidly increasing w(a) at z < 3**: up to w ~ 0.14 at z=0, meaning the DM density dilutes much faster than a^{-3} at late times
3. **Combined effect: 33% max TT deviation, 12.8 sigma per-ell** -- **DEFINITIVELY EXCLUDED by Planck**

### 5.3 Comparison with Planck Sensitivity

Planck constrains omega_cdm to ~1% precision (sigma ~ 0.0012).
The f(delta) correction at z=1100 changes effective omega_dm by 28.7%.
This is ~24 sigma beyond Planck sensitivity -- catastrophically excluded.

### 5.4 Re-normalization Strategy (Run E)

If omega_gdm is reduced from 0.12 to 0.0932 to match omega_cdm=0.12 at recombination:
- At z=1100: omega_dm_eff = 0.0932 x 1.2871 = 0.1200 (matches Planck)
- At z=0: omega_dm = 0.0932 (-22.3% less than CDM)
- **Result: still 5.84% max TT deviation (EXCLUDED)**
- The late-time w(a) effect changes the expansion history and ISW effect enough to be detected

### 5.5 What This Means for the Theory

The f(delta) correction with delta_0 = 0.34 is excluded by the CMB regardless of normalization strategy. This implies one or more of:

1. **delta_0 must be much smaller**: Run C (delta_0 = 0.10) gives 0.79% deviation (MARGINAL). Even smaller delta_0 ~ 0.03-0.05 might be safe.
2. **The f(delta) = (4+3*delta)/4 formula may not apply to the full Khronon action**: This was derived for the weak-field limit (rotation curves). The cosmological context may have a different f(delta) that is closer to 1 at all epochs.
3. **The omega_gdm normalization needs reinterpretation**: In the Khronon theory, omega_gdm is not a free parameter -- it is determined by the theory. The correct interpretation may give a much smaller effective f(delta) correction at the background level.
4. **The perturbation-level prediction (c_s^2 = 0) remains SAFE** (Level 1 result). The issue is purely at the background level.

### 5.6 Quantitative Exclusion Thresholds

| delta_0 | Max TT dev (11-bin) | Status | Notes |
|---------|---------------------|--------|-------|
| 0.34    | ~33%                | EXCLUDED | Current prediction |
| 0.10    | ~0.8% (6-bin)       | MARGINAL | Needs 11-bin check |
| ~0.05   | ~0.2% (estimated)   | DETECTABLE | Might survive |
| ~0.01   | ~0.01% (estimated)  | SAFE | Below cosmic variance |

---

## 6. Technical Discovery: c_s^2 = 0 Numerical Instability

**CLASS (gdm_class_public) hangs indefinitely when w > 0 and c_s^2 = 0.** The perturbation solver becomes singular. The fix is to use c_s^2 = 1e-10 as numerical regularization. This value is:
- 4 orders of magnitude below the TKS exclusion bound (3.4e-6)
- Completely undetectable by any current or planned CMB experiment
- Confirmed to have < 0.002% effect on C_l (from Level 1 comparison of c_s^2 = 0 vs 1e-8)

This is a numerical issue, not a physical one.

---

## 7. Plots Generated

All saved to: `/Users/akaihuangm1/Desktop/github/petz-recovery-unification/research/cmb_plots/`

1. **level2_physics_summary.png** -- 4-panel: delta(a), density ratio, w_eff, cumulative effect
2. **level2_weff_bins.png** -- w_eff(a) with 6-bin and 11-bin overlays
3. **level2_tt_comparison.png** -- C_l^TT spectra with residuals
4. **level2_nsigma.png** -- Significance in cosmic variance units
5. **level2_residuals_zoomed.png** -- Residuals by multipole region (low-l, peaks, tail)
6. **level2_all_spectra.png** -- TT, EE, TE, phiphi
7. **level2_deviation_summary.png** -- Bar chart of max deviations

---

## 8. GDM Bin Values Used

### Run A (6-bin, delta_0=0.34):
```
time_values_gdm = 1e-05, 0.0001, 0.001, 0.01, 0.1
w_values_gdm = 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 4.86853552e-03
```

### Run B (11-bin, delta_0=0.34):
```
time_values_gdm = 1e-05, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8
w_values_gdm = 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 1.95003459e-03, 1.68895627e-02, 7.46944639e-02, 1.30755102e-01, 1.45824926e-01
```

---

## 9. Technical Notes

### 8.1 The w(a) Mapping

The density of the GDM fluid evolves as:
rho_gdm(a) = rho_gdm_0 * a^{-3} * exp(-3 * int_1^a w(a')/a' da')

We want: rho_K(a) = rho_K_0 * a^{-3} * f(delta(a))/f(delta_0)

So: exp(-3 int w/a' da') = f/f_0
=> w(a) = -(1/3) * d ln(f/f_0) / d ln a

For the binned GDM, we average w(a) within each bin, weighted by a^{-3} (density weighting).

### 8.2 Limitations

1. The w(a) binning is an approximation -- the true evolution is continuous
2. For Run A/B/C: omega_gdm = 0.12 at z=0, so at z=1100 we have 1.2871x more DM
3. Runs E/F attempt to compensate by adjusting the normalization
4. A proper implementation would modify CLASS source code (background.c)

### 8.3 Reproducibility

All .ini files are in: `/Users/akaihuangm1/Desktop/github/gdm_class_public/level2_*.ini`
Output data in: `/Users/akaihuangm1/Desktop/github/gdm_class_public/output/level2_*_cl.dat`
Analysis script: `/Users/akaihuangm1/Desktop/github/gdm_class_public/level2_run.py`

---

## 10. Next Steps

1. **Immediate: delta_0 scan with 11-bin resolution** -- Run B shows that 11-bin resolution is essential. Scan delta_0 = {0.01, 0.03, 0.05, 0.10, 0.15, 0.20} with 11-bin w(a) to find the exclusion threshold.
2. **Re-examine the f(delta) derivation**: The weak-field formula f(delta) = (4+3*delta)/4 may not be the correct cosmological f(delta). Check whether the full Khronon action gives a different background correction.
3. **Level 3**: Modify CLASS `background.c` to implement rho_K(a) directly (avoids binning approximation).
4. **Parameter scan**: Grid over (omega_gdm, delta_0) with omega_gdm as a free parameter to find the best fit to Planck.
5. **Planck data comparison**: Download binned Planck 2018 data for formal chi^2.

---

## 11. Addendum: delta_0 Scan (11-bin resolution)

Performed a systematic scan over delta_0 = {0.01, 0.03, 0.05, 0.10, 0.15, 0.20, 0.25, 0.34} with 11-bin w(a) resolution.

### 11.1 Results Table

| delta_0 | rho_K/CDM(z=1100) | max w(a) | max\|dTT/TT\| | rms\|dTT/TT\| | N_sigma | Status |
|---------|-------------------|----------|---------------|---------------|---------|--------|
| 0.01 | 1.0105 | 0.006 | 1.88% | 0.87% | 0.7 | MARGINAL |
| 0.03 | 1.0311 | 0.017 | 5.37% | 2.51% | 2.1 | **EXCLUDED** |
| 0.05 | 1.0511 | 0.027 | 8.52% | 4.02% | 3.3 | **EXCLUDED** |
| 0.10 | 1.0986 | 0.052 | 15.16% | 7.30% | 5.9 | **EXCLUDED** |
| 0.15 | 1.1429 | 0.075 | 20.44% | 10.00% | 7.9 | **EXCLUDED** |
| 0.20 | 1.1843 | 0.095 | 24.69% | 12.26% | 9.6 | **EXCLUDED** |
| 0.25 | 1.2231 | 0.115 | 28.18% | 14.18% | 10.9 | **EXCLUDED** |
| 0.34 | 1.2871 | 0.146 | 33.05% | 16.97% | 12.8 | **EXCLUDED** |

### 11.2 Scaling Law

The maximum TT deviation scales as a power law:
```
max|dTT/TT| ~ delta_0^{0.82}
```

This is slightly sub-linear, consistent with the logarithmic nature of f(delta) = (4+3*delta)/4.

### 11.3 Exclusion Thresholds (Interpolated)

| Threshold | max\|dTT/TT\| | delta_0 limit |
|-----------|---------------|---------------|
| SAFE (< cosmic variance) | < 0.05% | delta_0 < 0.0002 |
| DETECTABLE (Planck edge) | < 0.5% | delta_0 < 0.0025 |
| EXCLUDED (> 5%) | < 5% | delta_0 < 0.028 |

### 11.4 Critical Conclusion

**The f(delta) = (4+3*delta)/4 background correction with delta_0 = 0.34 is ruled out at ~13 sigma by the CMB.** Even delta_0 = 0.03 is excluded at >5%. To be completely safe, delta_0 < 0.003 is required.

This means one of:
1. **The f(delta) formula does not apply at the cosmological background level** -- it may only be valid for the weak-field perturbative regime (rotation curves)
2. **The background evolution has a different normalization** -- the Khronon action may give f_cosmo(delta) much closer to 1 than f_gal(delta)
3. **delta_0 is much smaller than 0.34** -- which would reduce the rotation curve effect proportionally, potentially undermining the galaxy-scale predictions

**This is the most important constraint from Level 2: the background f(delta) correction, if taken at face value, is catastrophically excluded by the CMB. The theory must address this.**

### 11.5 Additional Plots

8. **level2_delta0_exclusion.png** -- Exclusion curve: max deviation and N_sigma vs delta_0
9. **level2_delta0_scan_residuals.png** -- TT residuals for all delta_0 values
10. **level2_delta0_scaling.png** -- Power law scaling with exclusion thresholds

---

## References

- Ilic, Kopp, Thomas & Skordis 2020, arXiv:2004.09572 (gdm_class_public)
- Thomas, Kopp & Skordis 2016, arXiv:1601.05097 (TKS bound)
- Blanchet & Skordis 2024, arXiv:2404.06584 (Khronon dark matter)
- Planck 2018, arXiv:1807.06209 (cosmological parameters)

---

*Last updated: 2026-03-17*
*This extends Level 0+1 by including the Khronon f(delta) background correction and a systematic delta_0 scan.*
