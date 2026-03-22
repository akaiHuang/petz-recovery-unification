# Precision Audit: mlx_class vs CLASS

**Date**: 2026-03-23
**Method**: Systematic comparison of every intermediate quantity against CLASS v3.x output.
**Reference files**:
- CLASS transfer functions: `gdm_class_public/output/lcdm_tk_tk.dat` (z=1100, 114 k-modes)
- CLASS C_l: `gdm_class_public/output/mlx_ref_lcdm_cl.dat` (l=2..2500, TT only)
- CLASS parameters: `gdm_class_public/mlx_ref_lcdm.ini` (Planck 2018, N_ur=3.046, no massive nu)

---

## 1. Background

| Quantity | Ours | CLASS | Error |
|----------|------|-------|-------|
| tau_rec [Mpc] | 280.43 | ~280.3 | 0.05% |
| r_s [Mpc] | 144.39 | ~144.39 | <0.01% |
| tau_0 [Mpc] | 14174.5 | ~14168 | 0.04% |
| D_A [Mpc] | 13894.0 | ~13888 | 0.04% |
| R_rec | 0.622 | 0.622 | <0.01% |
| l_1 = pi*D_A/r_s | 302 | 302 | <0.1% |
| Visibility FWHM [Mpc] | 22.5 | ~22 | ~2% |
| Visibility peak z | 1090 | 1089 | 0.1% |

**Verdict**: Background matches CLASS to 0.05% or better. Not a source of significant error.

---

## 2. Transfer Functions at z=1100 (Synchronous Gauge)

Normalization verified: For deeply superhorizon modes (k < 10^-3 Mpc^-1), our delta_c and delta_g match CLASS to 0.05%.

### 2a. Density perturbations (sync gauge)

| k [Mpc^-1] | delta_g ratio | delta_c ratio | delta_b ratio |
|-------------|:-------------:|:-------------:|:-------------:|
| 0.005 | 0.990 | 0.987 | 0.990 |
| 0.010 | 0.998 | 0.987 | 0.998 |
| 0.015 | 1.008 | 0.988 | 1.008 |
| 0.020 | 1.022 | 0.986 | 1.022 |
| 0.025 | 1.014 | 0.984 | 1.015 |
| 0.030 | 0.785* | 0.980 | 0.816* |
| 0.040 | 1.033 | 0.971 | 1.035 |
| 0.060 | 1.026 | 0.959 | 1.025 |
| 0.100 | 1.103 | 0.946 | 1.102 |

*Near zero crossing of delta_g -- ratio diverges at nodes, not a real error.

### 2b. Newtonian gauge potentials

| k [Mpc^-1] | Phi_N ratio | Psi_N ratio (current) | Psi_N ratio (factor/2 fix) |
|-------------|:-----------:|:---------------------:|:--------------------------:|
| 2.23e-5 | 1.0082 | 0.904 | 0.993 |
| 7.05e-5 | 1.0082 | 0.904 | 0.993 |
| 2.23e-4 | 1.0082 | 0.904 | 0.993 |
| 7.05e-4 | 1.0081 | 0.904 | 0.993 |
| 2.23e-3 | 1.0074 | 0.905 | 0.992 |
| 7.05e-3 | 1.0032 | 0.924 | 0.990 |
| 4.07e-2 | 0.966 | 0.958 | 0.962 |
| 8.73e-2 | 0.947 | 0.944 | 0.946 |

---

## 3. Root Causes (Ranked by Impact)

### CRITICAL: Anisotropic Stress Factor of 2

**File**: `perturbations_sync.py`, line 362-363
**Current code**:
```python
aniso = 12.0 * H02 / (a * a * k2) * (Og * sigma_g + On * sigma_n)
```
**Should be**:
```python
aniso = 6.0 * H02 / (a * a * k2) * (Og * sigma_g + On * sigma_n)
```

**Derivation**: The trace-free Einstein equation gives:
```
Psi = Phi - (9/2) * (a^2/k^2) * sum_i (1+w_i)*rho_i*sigma_i
```
where `rho_i = H0^2 * Omega_i / a^{3(1+w)}` in CLASS units.
Expanding: `(9/2) * (4/3) * H0^2 * Omega / (a^2 * k^2) * sigma = 6 * H0^2 * Omega / (a^2 * k^2) * sigma`.
Our code uses 12 instead of 6.

**Verification**: Confirmed against CLASS source code (`perturbations.c` line 6491):
```c
psi = phi - 4.5 * (a2/k2) * rho_plus_p_shear;
```
where `rho_plus_p_shear = (4/3)*rho*sigma` for radiation, giving factor = 4.5 * 4/3 = 6.

**Impact**:
- Psi_N is ~10% wrong for all superhorizon modes (Psi/Phi = 0.836 vs CLASS 0.932)
- Sachs-Wolfe source (Theta0+Psi) is **39% wrong** at k=0.005 and **13% wrong** at k=0.01
- Affects SW plateau (l < 30) by ~5-10%
- Affects first peak (l~220) by ~2-5% (these modes are borderline superhorizon)
- Negligible for k > 0.03 (1/k^2 suppression makes aniso stress irrelevant)

### CRITICAL: Neutrino l=2 Metric Source Wrong (eta' coefficient)

**File**: `perturbations_sync.py`, line 213-214
**Current code**:
```python
dFn[2] = ((k / 5.0) * (2.0 * Fn[1] - 3.0 * F3n)
          + (4.0 / 15.0) * h_prime + (8.0 / 15.0) * eta_prime)
```
**Should be**:
```python
dFn[2] = ((k / 5.0) * (2.0 * Fn[1] - 3.0 * F3n)
          + (4.0 / 15.0) * h_prime + (8.0 / 5.0) * eta_prime)
```

**Root cause**: The Ma & Bertschinger Eq. (28) neutrino l=2 source `(4/15)*h' + (8/15)*eta'` is **missing the metric shear contribution**. In synchronous gauge, the neutrino hierarchy l=2 equation requires the full metric source `(8/15)*(theta_nu + metric_shear)` where `metric_shear = (h' + 6*eta')/2`. This gives a source of `(4/15)*h' + (8/5)*eta'`, i.e., the eta' coefficient should be `8/5 = 1.6` not `8/15 = 0.533` -- a factor of 3 difference.

**Verification**: Confirmed against CLASS source code (`perturbations.c` line 9347-9352):
```c
dy[shear_ur] = 0.5*(8./15.*(theta_ur + metric_shear) - 3./5.*k*l3_ur);
```

**Numerical verification**: With the CLASS-like source, the CDM deficit is COMPLETELY FIXED:

| k [Mpc^-1] | delta_c ratio (M&B) | delta_c ratio (CLASS-like) |
|-------------|:-------------------:|:--------------------------:|
| 0.01 | 0.987 | 0.991 |
| 0.03 | 0.980 | **1.001** |
| 0.06 | 0.959 | **1.000** |
| 0.10 | 0.946 | **1.000** |
| 0.15 | 0.936 | **1.000** |

**Impact**: Fixes the 3-7% CDM transfer function deficit at all k. Also fixes the 0.8% Phi error (since Phi depends on delta_c through the constraint equations). **This is the single most impactful fix.**

Note: The same fix should also be applied to the PHOTON l=2 equation (line 185-187), changing `(8.0 / 15.0) * eta_prime` to `(8.0 / 5.0) * eta_prime`.

### MODERATE: Phi_N 0.8% Too High (CAUSED BY Bug #2)

**Symptom**: Phi_N/Phi_CLASS = 1.008 for all superhorizon modes.

**Cause**: The neutrino quadrupole error (Bug #2) propagates: wrong F_n,2 -> wrong neutrino contribution to eta' and h' -> wrong alpha -> wrong Phi_N. **This will be fixed automatically by fixing Bug #2.**

**Impact**: ~0.8% systematic error on all Newtonian gauge potentials.

### MINOR: Photon Phase Shift (~1%)

**Zero crossing positions**:

| Node | CLASS k [Mpc^-1] | Ours k [Mpc^-1] | Shift |
|------|:-----------------:|:----------------:|:-----:|
| 1st | 0.0314 | 0.0311 | -0.9% |
| 2nd | 0.0506 | 0.0503 | -0.5% |
| 3rd | 0.0748 | 0.0744 | -0.5% |
| 4th | 0.0954 | 0.0952 | -0.2% |
| 5th | 0.1186 | 0.1183 | -0.3% |

**Cause**: Effective sound speed ~1% too fast, possibly from the CDM deficit changing the photon-baryon coupling through the gravitational potential.

**Impact**: Creates large ratio errors near acoustic nodes (irrelevant for C_l since the contribution is zero at nodes) but slightly shifts peak positions.

### MINOR: Missing Baryon Sound Speed

**File**: `perturbations_sync.py`, line 170
**Current**: `c_s^2 * k^2 * delta_b` term omitted (cold baryon approximation)
**Effect**: At z~1100, c_b^2 ~ 2.6e-7, so `c_b^2 * k^2 * delta_b << |kd|*(theta_g-theta_b)/R`. Negligible for k < 0.3.

### MINOR: Missing Photon Polarization Hierarchy

Our photon l=2 equation uses `-(9/10)*|kd|*F_2` for the Thomson damping.
CLASS includes polarization: `-(9/10)*|kd|*(F_2 - sqrt(6)*Pi)` where `Pi = (F_2 + G_0 + G_2)/8`.
This changes the effective Silk damping scale, affecting the damping tail (l > 1000).

---

## 4. Source Functions at Visibility Peak

| k [Mpc^-1] | (Th0+Psi) current | (Th0+Psi) with factor/2 fix | Error (current) |
|-------------|:------------------:|:---------------------------:|:---------------:|
| 0.005 | 0.070 | 0.115 | 39% |
| 0.010 | -0.213 | -0.188 | 13.4% |
| 0.015 | -0.489 | -0.484 | 1.1% |
| 0.020 | -0.580 | -0.586 | 1.0% |
| 0.030 | -0.051 | -0.055 | 6.3% |
| 0.040 | 0.499 | 0.500 | 0.1% |
| 0.060 | -0.629 | -0.628 | 0.1% |
| 0.100 | -0.327 | -0.326 | 0.1% |

---

## 5. C_l Comparison

| l | D_l CLASS [uK^2] | Feature |
|---|:-----------------:|---------|
| 2-30 | ~1000-1200 | SW plateau (affected by Bug #1: ~5-10%) |
| 221 | 5740 | 1st peak |
| 538 | 2606 | 2nd peak |
| 2nd/1st ratio | 0.454 | Baryon loading diagnostic |

Current mlx_class: first peak 1.03x, second peak 0.59x CLASS (from prior run).

---

## 6. Recommended Fix Priority

1. **Fix neutrino (AND photon) l=2 metric source** (`perturbations_sync.py` lines 186, 214): Change eta' coefficient from `8/15` to `8/5` in both photon and neutrino l=2 equations. This is the highest-impact fix: **completely eliminates the 3-7% CDM deficit, fixes the 0.8% Phi error, and likely fixes the 1% photon phase shift.** Estimated improvement: **5-7% in peak-normalized RMS**.

2. **Fix anisotropic stress factor** (12 -> 6 in `perturbations_sync.py` line 362): Fixes Psi_N for superhorizon modes, SW plateau, low-l C_l. **Estimated improvement: 2-3% in peak-normalized RMS** (mostly at low l).

3. **Add photon polarization hierarchy** (G_0, G_2): Improves Silk damping accuracy for l > 800. **Estimated improvement: 1-2% in damping tail.**

4. **Add baryon sound speed term**: `c_b^2 * k^2 * delta_b` in theta_b equation. Minor but easy fix.

**Combined estimate**: Fixes #1 and #2 together should reduce peak-normalized RMS from 7.4% to ~1-2%.

---

## 7. Raw Data: CLASS Transfer Functions at z=1100

CLASS columns: k[h/Mpc], d_g, d_b, d_cdm, d_ur, d_tot, phi, psi
All normalized to initial curvature perturbation = 1.
Our normalization: C = 1 in adiabatic ICs (eta -> 1 for superhorizon).
Normalization match verified to 0.05% for superhorizon modes.
