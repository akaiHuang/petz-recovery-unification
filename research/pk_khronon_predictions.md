# P(k) Predictions for Khronon DBI Dark Matter

**Date**: 2026-03-24
**Code**: CLASS + Khronon DBI (alpha_dbi parameter)
**Cosmology**: Planck 2018 best-fit (h=0.6736, omega_b=0.02237, omega_cdm=0.12)

## Executive Summary

We computed the matter power spectrum P(k) for the Khronon DBI dark matter model
across alpha_dbi = 0 (LCDM) to 0.5, spanning 10 orders of magnitude.

**KEY FINDING**: The DBI pressure mechanism is MUCH more powerful than naively expected.
Even alpha_dbi = 1e-7 produces 4.5% suppression at k = 0.2 h/Mpc and reduces sigma_8
by 2.6%. The original scan range (alpha = 0.001 to 0.5) is catastrophically excluded --
those values destroy ALL small-scale structure.

**Upper bound**: alpha_dbi < ~1e-7 from existing P(k) data (5% at k=0.2 h/Mpc).

---

## 1. P(k) Ratio vs LCDM

### DeltaP/P = P(k, Khronon) / P(k, LCDM) - 1

| k [h/Mpc] | alpha=1e-10 | alpha=1e-9 | alpha=1e-8 | alpha=5e-8 | alpha=1e-7 |
|------------|-------------|------------|------------|------------|------------|
| 0.001      | -9.3e-10    | -9.3e-9    | -9.3e-8    | -4.7e-7    | -9.3e-7    |
| 0.01       | -1.1e-7     | -1.1e-6    | -1.1e-5    | -5.5e-5    | -1.1e-4    |
| 0.05       | -2.8e-6     | -2.8e-5    | -2.8e-4    | -1.4e-3    | -2.8e-3    |
| 0.1        | -9.0e-6     | -1.1e-4    | -1.1e-3    | -5.7e-3    | -1.1%      |
| 0.2        | -4.6e-5     | -4.6e-4    | -4.5e-3    | -2.3%      | -4.5%      |
| 0.5        | -2.9e-4     | -2.8e-3    | -2.8%      | -13.3%     | -24.9%     |
| 1.0        | -1.0e-3     | -1.0%      | -9.6%      | -40.0%     | -64.5%     |
| 5.0        | -3.4%       | -29.1%     | -97.5%     | -99.8%     | -99.9%     |
| 10.0       | -12.7%      | -75.6%     | -99.7%     | -100%      | -100%      |

### Jeans Scale (where suppression reaches 1%, 10%, 50%)

| alpha_dbi | k(1%) [h/Mpc] | k(10%) [h/Mpc] | k(50%) [h/Mpc] |
|-----------|---------------|-----------------|-----------------|
| 1e-10     | 3.46          | 10.96           | >11             |
| 1e-9      | 0.94          | 3.46            | 8.70            |
| 1e-8      | 0.30          | 1.11            | 2.75            |
| 5e-8      | 0.14          | 0.43            | 1.11            |
| 1e-7      | 0.094         | 0.31            | 0.84            |
| 1e-6      | 0.030         | 0.098           | 0.25            |
| 1e-5      | 0.010         | 0.031           | 0.079           |

### Why suppression is so strong

The Khronon DBI sound speed from the code is:

```
c_s^2(k,a) = alpha_dbi * k^2 / (k_J(a)^2 + k^2)
k_J^2(a) = 1.5 * Omega_cdm * H_0^2 / a
```

At z=0: **k_J = 0.00021 h/Mpc** (lambda_J = 29,906 Mpc/h ~ Hubble horizon)

This means c_s^2 = alpha_dbi for essentially ALL observable scales (k > 0.001 h/Mpc).
The pressure acts on every sub-horizon mode, and the effect accumulates over ~10 e-folds
from matter-radiation equality to today. The cumulative suppression is exponential
in alpha_dbi, not linear.

---

## 2. sigma_8 as Function of alpha_dbi

| alpha_dbi | sigma_8  | Delta%  | S_8      |
|-----------|----------|---------|----------|
| 0 (LCDM) | 0.8317   | --      | 0.8525   |
| 1e-10     | 0.8317   | -0.003% | 0.8525   |
| 1e-9      | 0.8315   | -0.028% | 0.8523   |
| 1e-8      | 0.8294   | -0.28%  | 0.8502   |
| 5e-8      | 0.8207   | -1.33%  | 0.8412   |
| **1e-7**  | **0.8104** | **-2.57%** | **0.8306** |
| 1e-6      | 0.6781   | -18.5%  | 0.6950   |
| 1e-5      | 0.3421   | -58.9%  | 0.3506   |

Reference values:
- Planck 2018: sigma_8 = 0.811 +/- 0.006, S_8 = 0.832 +/- 0.013
- KiDS-1000:  S_8 = 0.759 +/- 0.024
- DES Y3:     S_8 = 0.776 +/- 0.017

### S_8 tension resolution

To bring sigma_8 from 0.832 down to the weak lensing values:
- sigma_8 = 0.80 requires alpha_dbi ~ 1.2e-7
- sigma_8 = 0.78 requires alpha_dbi ~ 1.7e-7
- sigma_8 = 0.77 requires alpha_dbi ~ 2.0e-7
- sigma_8 = 0.76 requires alpha_dbi ~ 2.4e-7

**HOWEVER**: alpha_dbi ~ 2e-7 gives |DeltaP/P| ~ 8% at k=0.2 h/Mpc, which is
EXCLUDED by current BOSS/eBOSS data. The constant-alpha_dbi model cannot resolve
the S_8 tension without violating P(k) constraints at BAO scales.

The fundamental issue: the DBI pressure suppresses ALL sub-horizon modes equally.
To resolve S_8, one needs scale-DEPENDENT suppression (small at large scales, large
at sigma_8 scales). This requires the time-binned GDM approach or k-dependent alpha.

---

## 3. Predictions for Surveys

### Observational Upper Bound (current data)

**alpha_dbi < 1.1e-7** (from requiring |DeltaP/P| < 5% at k = 0.2 h/Mpc)

### Rubin LSST (Y10, 3x2pt)
- k = 0.01 -- 1.0 h/Mpc, z = 0.2-1.2, P(k) precision ~1%
- **Detectable for alpha_dbi > 1e-9** (marginal) to **> 1e-8** (definitive)
- At alpha = 1e-8: max suppression 9.6% at k=0.94 h/Mpc
- At alpha = 5e-8: max suppression 40% at k=0.94 h/Mpc

### Euclid (spectroscopic)
- k = 0.01 -- 0.3 h/Mpc, z = 0.9-1.8, P(k) precision ~0.5%
- **Detectable for alpha_dbi > 1e-8** (marginal) to **> 5e-8** (definitive)
- Best sensitivity at k ~ 0.2-0.3 h/Mpc

### Euclid (photometric)
- k = 0.001 -- 5.0 h/Mpc, z = 0.2-2.5, P(k) precision ~1%
- **Detectable for alpha_dbi > 1e-10** due to wide k-range
- At k > 1 h/Mpc, even alpha = 1e-10 gives >1% suppression

### DESI (5yr, complete)
- k = 0.01 -- 0.25 h/Mpc, z = 0.1-4.2, P(k) precision ~0.3%
- **Detectable for alpha_dbi > 1e-8**
- BAO measurements at 0.3% give the tightest constraint per mode

### SPHEREx
- k = 0.003 -- 1.0 h/Mpc, z = 0.0-5.0, P(k) precision ~1%
- **Detectable for alpha_dbi > 1e-9** (marginal) to **> 1e-8** (definitive)

### Summary: Detectability Thresholds

| Survey              | Marginal detection | Definitive detection |
|--------------------|--------------------|---------------------|
| Euclid (photo)     | alpha > 1e-10      | alpha > 1e-9        |
| Rubin LSST         | alpha > 1e-9       | alpha > 1e-8        |
| SPHEREx            | alpha > 1e-9       | alpha > 1e-8        |
| Euclid (spectro)   | alpha > 1e-8       | alpha > 5e-8        |
| DESI complete      | alpha > 1e-8       | alpha > 5e-8        |

---

## 4. The Observationally Interesting Window

The window alpha_dbi ~ 1e-8 to 1e-7 is particularly interesting:

| Property | alpha = 1e-8 | alpha = 5e-8 | alpha = 1e-7 |
|----------|-------------|-------------|-------------|
| sigma_8  | 0.8294      | 0.8207      | 0.8104      |
| DeltaP/P at k=0.1 | -0.11% | -0.57% | -1.14% |
| DeltaP/P at k=0.2 | -0.45% | -2.25% | -4.46% |
| DeltaP/P at k=0.5 | -2.8%  | -13.3% | -24.9% |
| DeltaP/P at k=1.0 | -9.6%  | -40.0% | -64.5% |
| Current constraint | ALLOWED | ALLOWED | MARGINAL |
| Euclid detection   | marginal | YES | YES |

This window:
- Is ALLOWED by current data (barely, for alpha=1e-7)
- DETECTABLE by upcoming surveys
- Gives sigma_8 shift of 0.3-2.6% (small compared to S_8 tension but measurable)
- Provides a clear prediction: P(k) suppression increasing monotonically with k

---

## 5. Implications for the Theory

### What this analysis reveals

1. **The constant-alpha_dbi parametrization is too simple**: Because k_J is at the
   horizon scale, alpha_dbi acts as a constant c_s^2 at all observable scales. This
   gives warm-dark-matter-like suppression, which is already tightly constrained.

2. **The BS2025 GDM approach is essential**: The time-binned w(tau), cs2(tau) approach
   allows the effective sound speed to vary with time, which translates to scale-dependent
   effects via the growth history. This is physically distinct from constant alpha_dbi.

3. **The observational bound alpha_dbi < 1e-7 translates to c_s^2 < 1e-7**: At z=0,
   the DBI sound speed is c_s ~ 3e-4 * c. This is already much tighter than the
   GDM constraints from Planck (cs2 < 1e-6 at 2sigma for constant cs2).

4. **S_8 tension cannot be addressed**: To resolve S_8, one needs to reduce sigma_8
   without affecting P(k) at k < 0.1 h/Mpc. The constant-alpha_dbi model suppresses
   all scales equally, making this impossible.

### Prediction for the DBI Khronon physical model

In the full DBI Khronon theory (Blanchet & Skordis 2025), the effective alpha_dbi
is NOT constant but depends on the ghost condensation physics. The time-dependent
cs2(tau) from the DBI structure naturally gives small cs2 at early times (CDM-like)
and potentially larger cs2 at late times, creating the scale-dependent suppression
needed. The GDM-binned approach captures this correctly.

The constant alpha_dbi parameter should be viewed as a diagnostic tool for understanding
the CUMULATIVE effect of DBI pressure, not as the physical model prediction.

---

## 6. Files Generated

- `output/pk_scan_a{X}_pk.dat`: P(k) output files for each alpha value
- `analyze_pk_scan.py`: Full analysis script
- `pk_scan_a{X}.ini`: CLASS input files

Alpha values scanned: 0, 1e-10, 1e-9, 1e-8, 5e-8, 1e-7, 1e-6, 1e-5, 1e-4,
0.001, 0.005, 0.01, 0.05, 0.5
