# Self-Consistent delta(z) for Khronon Dark Matter

**Date**: 2026-03-17
**Status**: COMPLETED — reveals a fundamental 17% excess at CMB from the nonlinear delta^2 term

## Executive Summary

We solved the coupled Friedmann + Khronon equations self-consistently. The key findings:

1. **Self-consistency reduces the naive 20.5% CMB excess to 17.1%** — improvement is real but insufficient
2. **Root cause**: the `delta^2` term in `rho_K = H^2 delta(2+delta)/(8piG)` makes the Khronon energy density grow *faster* than `(1+z)^3`
3. **The excess is intrinsic to the nonlinear Khronon kinetic structure** — not an artifact of approximation
4. **z_eq shifts from 3434 (LCDM) to 3645 (self-consistent)**, a 6.1% change detectable by Planck at >>100 sigma
5. **w_eff(z) is NOT exactly zero**: ranges from +0.098 (z=0) to -0.030 (z~5000)

## The Coupled System

### Equations

Friedmann (without CDM, with Khronon):
```
H^2 = H^2_std / [1 - delta(2+delta)/3]
```
where `H^2_std = H_0^2 [Omega_b(1+z)^3 + Omega_r(1+z)^4 + Omega_Lambda]`

Conservation law (from Khronon EOM):
```
H^2 delta = I_0 (1+z)^3 / 2    (exact)
```

Combined: at each z, solve for delta from
```
delta [1 + x(2+delta)/3] = x
```
where `x(z) = I_0 c^2 / (2 H^2_std(z) a^3(z))` is the "naive" delta ignoring self-consistency.

### Solution

This reduces to a **quadratic** (not cubic):
```
(x/3) delta^2 + (1 + 2x/3) delta - x = 0
```
Positive root:
```
delta = [-(1+2x/3) + sqrt((1+2x/3)^2 + 4x^2/3)] / (2x/3)
```

### Calibration

At z=0: `Omega_DM = delta_0(2+delta_0)/3 = 0.265` gives:
- **delta_0 = 0.339776**
- **x_0 = delta_0/(1-Omega_DM) = 0.462280**

## Results

### Table 1: delta(z) and rho_K/rho_CDM at key redshifts

| z | delta_SC | delta_naive | Omega_K | rho_K/rho_CDM | H_SC/H_LCDM | Note |
|---:|---:|---:|---:|---:|---:|---|
| 0 | 0.3398 | 0.3398 | 0.265 | 1.000 | 1.000 | by construction |
| 1 | 0.759 | 0.850 | 0.698 | 1.179 | 1.058 | transition |
| 10 | 0.897 | 1.077 | 0.867 | 1.238 | 1.096 | matter era |
| 100 | 0.883 | 1.051 | 0.849 | 1.232 | 1.091 | matter era |
| 500 | 0.820 | 0.944 | 0.771 | 1.205 | 1.073 | matter era |
| **1100** | **0.739** | **0.819** | **0.675** | **1.171** | **1.053** | **CMB** |
| 3400 | 0.526 | 0.544 | 0.443 | 1.080 | 1.017 | z_eq |
| 10000 | 0.278 | 0.277 | 0.211 | 0.974 | 0.997 | radiation era |

### Table 2: Energy budget (self-consistent)

| z | Omega_K | Omega_b | Omega_r | Omega_Lambda | Total |
|---:|---:|---:|---:|---:|---:|
| 0 | 0.265 | 0.049 | 0.0001 | 0.686 | 1.000 |
| 100 | 0.849 | 0.127 | 0.024 | 0.000 | 1.000 |
| 1100 | 0.675 | 0.107 | 0.219 | 0.000 | 1.000 |
| 3400 | 0.443 | 0.076 | 0.481 | 0.000 | 1.000 |
| 10000 | 0.211 | 0.040 | 0.749 | 0.000 | 1.000 |

### Table 3: Effective equation of state

| z | w_eff |
|---:|---:|
| 0 | +0.098 |
| 1 | +0.044 |
| 10 | +0.000 |
| 100 | -0.002 |
| 1100 | -0.017 |
| 3400 | -0.031 |
| 10000 | -0.030 |

## Root Cause Analysis

### Why rho_K/rho_CDM != 1

The conservation law `H^2 delta = C (1+z)^3` guarantees that the **linear** part of rho_K scales correctly:
```
rho_K = H^2 delta(2+delta)/(8piG) = C(1+z)^3 (2+delta) / (8piG)
```

For rho_K to scale as (1+z)^3 (i.e., w=0), we would need `2+delta(z) = const`. But delta **varies**:
- z=0: delta = 0.34, so 2+delta = 2.34
- z=100: delta = 0.88, so 2+delta = 2.88
- z=10000: delta = 0.28, so 2+delta = 2.28

The (2+delta) factor varies by ~25% across the relevant range, producing:
```
w_eff = -(1/3) d ln(2+delta) / d ln(1+z)
```

### Quantitative: the delta^2 contribution

In the deep matter era:
- x -> delta_0/Omega_b = 6.93 (large!)
- delta -> 0.900
- delta^2/(2 delta) = 0.45 — the quadratic term is **45% of the linear term**
- Omega_K = 0.870, but need Omega_K = 0.844 for H_SC = H_LCDM
- **Mismatch: 3.1%** in Omega_K, amplified to **24% excess in rho_K/rho_CDM**

### Analytic matter-era formula

In pure matter era (z >> 1, z << z_eq):
```
rho_K / rho_CDM = Omega_b * Omega_K / [(1 - Omega_K) * Omega_DM]
```
where Omega_K = delta_m(2+delta_m)/3 and delta_m solves the quadratic with x = delta_0/Omega_b.

**Result: rho_K/rho_CDM = 1.2395 in the deep matter era** (z ~ 10-50).

This is a **constant** ratio — the Khronon tracks CDM with a fixed 24% offset during matter domination.

## CMB Impact Assessment

### Direct impact at z=1100

| Quantity | Self-consistent | LCDM | Difference |
|---|---:|---:|---:|
| rho_K/rho_CDM | 1.171 | 1.000 | +17.1% |
| H(z=1100) | 1.053 H_LCDM | H_LCDM | +5.3% |
| z_eq | 3645 | 3434 | +6.1% |
| theta_s shift | ~5.3% | — | ~177 sigma |

### What this means for CMB power spectrum

1. **Sound horizon**: rs is ~5% smaller (H is larger), shifting all peaks by ~5% in multipole space
2. **z_eq shift**: the 1st-to-3rd peak ratio depends sensitively on z_eq; a 6% shift is catastrophic
3. **Silk damping**: higher H at decoupling means more damping at high l
4. **ISW effect**: modified late-time evolution changes the low-l spectrum

All effects are many orders of magnitude above Planck sensitivity.

## Parameter scan: can varying Omega_DM help?

| Omega_DM | delta_0 | rho_K/rho_CDM(z=1100) | Excess |
|---:|---:|---:|---:|
| 0.10 | 0.140 | 1.179 | +17.9% |
| 0.15 | 0.204 | 1.190 | +19.0% |
| 0.20 | 0.265 | 1.185 | +18.5% |
| 0.265 | 0.340 | 1.171 | +17.1% |
| 0.30 | 0.378 | 1.161 | +16.1% |
| 0.35 | 0.432 | 1.145 | +14.5% |
| 0.40 | 0.483 | 1.129 | +12.9% |

**The excess is 13-19% regardless of Omega_DM**. The problem is structural.

## Comparison: self-consistent vs naive

| Method | rho_K/rho_CDM at z=1100 | w_eff(z=1100) |
|---|---:|---:|
| Naive (LCDM H(z)) | 1.205 (+20.5%) | ~-0.02 |
| Self-consistent | 1.171 (+17.1%) | -0.017 |

Self-consistency improves things by ~3.4%, but the dominant effect — the nonlinear delta^2 term — remains.

## Critical Insight: Three Possible Resolutions

### (a) Modified kinetic function K(Y)

The ghost condensation kinetic function K(Y) determines the exact form of rho_K. If K(Y) has **higher-order terms beyond Y^2**, they could modify the `delta(2+delta)` structure:
```
rho_K = f(delta) * H^2 / (8piG)
```
where f(delta) could differ from delta(2+delta) at large delta. This is the most natural resolution within the Khronon framework.

### (b) Perturbation-level treatment

If the Khronon perturbation is treated at the linearized level (delta << 1), then rho_K ~ 2H^2 delta/(8piG) exactly, and the delta^2 problem disappears. However:
- At z=0, we need Omega_DM = 2 delta_0/3 = 0.265, requiring delta_0 = 0.40 (not small)
- This is inconsistent with a perturbative expansion

### (c) Running mu(k) = k absorbs the excess

If mu has scale-dependent running, the background mu_0 = H_0/c could set the z=0 value, while the perturbation spectrum could differ. The excess in background rho_K would then need to be compensated by perturbation-level effects in the Boltzmann hierarchy.

## Plots

All plots saved in `cmb_plots/`:
- `self_consistent_delta_v2.png` — 6-panel overview (delta, ratio, w_eff, H ratio, energy fractions, Omega comparison)
- `self_consistent_ratio_cmb.png` — Focused CMB-era rho_K/rho_CDM ratio

## Scripts

- `self_consistent_delta.py` — v1 basic computation
- `self_consistent_delta_v2.py` — v2 with root cause analysis, parameter scan, plots

## Summary

**Self-consistency reduces but does NOT resolve the CMB tension.**

The fundamental issue is that `rho_K = H^2 delta(2+delta)/(8piG)` contains a nonlinear `delta^2` term that makes the Khronon energy density scale differently from (1+z)^3. At z=1100, this produces a 17% excess in rho_K over what CDM would give, resulting in a ~5% shift in H(z) and ~6% shift in z_eq — both catastrophically detectable by Planck.

The resolution likely requires either:
1. A modified K(Y) that linearizes the delta dependence, or
2. A perturbation-level Boltzmann treatment where the background mismatch is absorbed into the perturbation spectrum, or
3. Accepting that the pure background-level Khronon is only an approximation, and the full theory operates at the perturbation level (delta k-modes, not a homogeneous delta)

**This is a CRITICAL gap that must be addressed before claiming CMB compatibility.**
