# BAO Compatibility Analysis for Khronon/tau Framework

**Author**: Sheng-Kai Huang (with computational analysis)
**Date**: 2026-03-16
**Status**: Complete quantitative analysis
**Assessment**: SAFE (with caveats on late-time DM density)

---

## Executive Summary

The Khronon/tau framework is **fully compatible** with Baryon Acoustic Oscillation measurements, including the latest DESI DR1 results. The sound horizon r_s is identical to LCDM to 10^{-11} precision. Surprisingly, the Khronon model with its late-time w_tilde correction provides a **better fit** to DESI BAO distance measures than flat LCDM (Delta chi^2 = -4.72, Khronon preferred).

The key mechanism: The running mu_bg(z) = H(z)/c ensures dust-like behavior (w_tilde ~ 10^{-10}) at z=1100, preserving the sound horizon and CMB. At z=0, w_tilde ~ 0.20, making the Khronon density 20% higher than CDM. This requires 5.2% less Lambda to maintain flatness, producing an effective quintessence-like expansion (w_eff ~ -0.95) that happens to better fit the DESI data at z ~ 0.5-1.0.

---

## 1. Sound Horizon

### Formula
```
r_s = integral_0^{z_drag} c_s(z) / H(z) dz

c_s(z) = c / sqrt(3(1 + R_b(z)))
R_b(z) = 3 rho_b / (4 rho_gamma) = (3 Omega_b)/(4 Omega_gamma) * 1/(1+z)
```

### Result
```
r_s (LCDM):    147.197 Mpc
r_s (Khronon): 147.197 Mpc
Delta r_s / r_s = -1.2 x 10^{-11}
Planck uncertainty: +/- 1.8 x 10^{-3}
```

**SAFE by a factor of ~10^8 relative to Planck uncertainty.**

### Why r_s is unchanged
1. The sound horizon integral is dominated by z > 1000
2. At z=1100: w_tilde = 3.6 x 10^{-10} (negligible)
3. The Khronon background density matches CDM density to 10^{-10} at high z
4. H(z) at z > 100 is identical to LCDM to better than 10^{-7}

---

## 2. Distance Measures vs DESI DR1

### DESI BAO Measurements (arXiv:2404.03002)

| Tracer | z_eff | Quantity | DESI | LCDM | Khr | LCDM pull | Khr pull |
|--------|-------|----------|------|------|-----|-----------|----------|
| BGS | 0.295 | D_V/r_d | 7.93 +/- 0.15 | 8.06 | 8.00 | +0.87sigma | +0.45sigma |
| LRG1 | 0.510 | D_M/r_d | 13.62 +/- 0.25 | 13.51 | 13.39 | -0.45sigma | -0.92sigma |
| LRG1 | 0.510 | D_H/r_d | 20.98 +/- 0.61 | 22.76 | 22.45 | +2.92sigma | +2.41sigma |
| LRG2 | 0.706 | D_M/r_d | 16.85 +/- 0.32 | 17.71 | 17.54 | +2.70sigma | +2.15sigma |
| LRG2 | 0.706 | D_H/r_d | 20.08 +/- 0.60 | 20.20 | 19.92 | +0.19sigma | -0.27sigma |
| LRG3+ELG1 | 0.930 | D_M/r_d | 21.71 +/- 0.28 | 21.94 | 21.71 | +0.83sigma | +0.00sigma |
| LRG3+ELG1 | 0.930 | D_H/r_d | 17.88 +/- 0.35 | 17.64 | 17.42 | -0.69sigma | -1.33sigma |
| ELG2 | 1.317 | D_M/r_d | 27.79 +/- 0.69 | 28.05 | 27.75 | +0.38sigma | -0.05sigma |
| ELG2 | 1.317 | D_H/r_d | 13.82 +/- 0.42 | 14.12 | 13.99 | +0.72sigma | +0.39sigma |
| QSO | 1.491 | D_V/r_d | 26.07 +/- 0.67 | 26.07 | 25.81 | -0.00sigma | -0.39sigma |
| Lya | 2.330 | D_M/r_d | 39.71 +/- 0.94 | 39.23 | 38.86 | -0.51sigma | -0.91sigma |
| Lya | 2.330 | D_H/r_d | 8.52 +/- 0.17 | 8.63 | 8.60 | +0.67sigma | +0.45sigma |

### Chi-squared comparison
```
chi^2 (LCDM):    19.36  (12 data points)
chi^2 (Khronon): 14.63  (12 data points)
Delta chi^2 = -4.72 (Khronon PREFERRED)
```

**Note**: This is a naive chi^2 without the full D_M-D_H covariance matrix. The DESI measurements have anti-correlations (rho ~ -0.4) between D_M and D_H in the same bin that would modify this comparison.

### Summary: Khronon better in 7/12 measurements, worse in 4/12, similar in 1/12.

---

## 3. Khronon Background Model

### The w_tilde correction
With running mu_bg(z) = H(z)/c:
```
w_tilde(z) = (3 Omega_c / 4) * (H_0 / H(z))^2

w_tilde(z=0)    = 0.199
w_tilde(z=0.5)  = 0.112
w_tilde(z=1.0)  = 0.052
w_tilde(z=10)   = 4.7 x 10^{-4}
w_tilde(z=1100) = 3.6 x 10^{-10}
```

### Impact on Friedmann equation
```
H^2 = H_0^2 [Omega_r(1+z)^4 + Omega_b(1+z)^3 + Omega_c(1+z)^3(1+w_tilde(z)) + Omega_Lambda,Khr]

Omega_Lambda,Khr = 1 - Omega_b - Omega_c(1+w_tilde(0)) - Omega_r = 0.634
Omega_Lambda,LCDM = 0.686
Delta Omega_Lambda = -0.052
```

### H(z) comparison
```
z     | Delta H/H (Khronon vs LCDM)
0.3   | +1.1%
0.5   | +1.4%
0.7   | +1.4%
1.0   | +1.3%
1.5   | +0.8%
2.3   | +0.4%
10    | +0.01%
1100  | +6.7 x 10^{-11}
```

The 1-1.4% excess in H(z) at z ~ 0.5-1.0 is what drives the improved DESI fit.

---

## 4. Effective Dark Energy Equation of State

An observer fitting w_0-w_a to the Khronon expansion history would infer:
```
z     | w_eff
0.1   | -0.944
0.3   | -0.942
0.5   | -0.946
0.7   | -0.953
1.0   | -0.963
1.5   | -0.977
2.0   | -0.985
```

This is a nearly constant quintessence-like w ~ -0.95.

### Comparison with DESI dark energy signal
- DESI + CMB + PantheonPlus: w_0 = -0.727, w_a = -1.05
- The DESI signal is MUCH stronger (w crossing -1, phantom at high z)
- Khronon predicts a milder, monotonic w_eff(z)
- **Khronon cannot explain the DESI dark energy signal**
- But Khronon is orthogonal to it: Khronon replaces CDM, not Lambda
- If DESI's signal is real, Khronon + evolving DE works fine

---

## 5. BAO Damping and Peak Broadening

### Silk damping
```
Silk damping scale: k_S ~ 0.15 Mpc^{-1}
Khronon Jeans scale: k_J ~ 3.3 x 10^8 Mpc^{-1}
Ratio: k_J / k_S ~ 2 x 10^9 --> NO effect
```

### BAO peak broadening
```
BAO scale: k_BAO ~ 2pi/r_s ~ 0.04 Mpc^{-1}
Khronon Jeans scale: k_J ~ 3.3 x 10^8 Mpc^{-1}
Ratio: k_J / k_BAO ~ 10^{10} --> NO effect
```

### Non-linear broadening (z < 1)
The BAO peak broadening from non-linear structure growth is controlled by the displacement field Sigma_nl, which depends on the linear growth rate. Khronon modifies this by O(1%) through its slightly different growth factor. This is negligible compared to the ~5 Mpc/h non-linear broadening scale.

**All BAO damping mechanisms: SAFE.**

---

## 6. Alcock-Paczynski Effect

```
F_AP = D_M(z) * H(z) / c

z     | Delta F_AP / F_AP
0.3   | +0.5%
0.5   | +0.5%
0.7   | +0.4%
0.9   | +0.2%
1.3   | -0.1%
1.5   | -0.2%
2.3   | -0.5%
```

Current DESI precision on F_AP: ~1-3%. Khronon effect is sub-percent.

**SAFE.**

---

## 7. New Interpretation: BAO as a tau Phenomenon

The BAO imprint represents the maximum distance sound waves traveled in the baryon-photon plasma before recombination. In the tau framework:

1. **Baryon-photon plasma**: Open system with Sigma > 0 (photon diffusion, Silk damping = irreversible entropy production). The plasma has tau > 0 -- time flows, acoustic waves propagate, information is lost to the photon diffusion scale.

2. **Khronon (dark sector)**: At z > 1000, w_tilde ~ 10^{-10}, c_s^2 ~ 10^{-10}. This means tau_dark ~ 0 -- the dark sector is essentially a CLOSED system (no entropy production, no information loss, no sound propagation). Dark matter "sits there" and gravitates without participating in the acoustic dynamics.

3. **The BAO ruler is a tau_baryon ruler**: The sound horizon r_s measures the causal reach of the baryon-photon tau (the time arrow of the luminous sector). The dark sector's tau is zero -- it has no "arrow of time" at the perturbation level -- so it does not contribute to the acoustic scale.

4. **Why BAO is preserved**: Any dark matter replacement that has tau_dark ~ 0 (equivalently, c_s^2 ~ 0) at the CMB epoch will preserve the BAO ruler. The Khronon achieves this through its quadratic Lagrangian K(Q) = mu^2(Q-1)^2, which gives omega = 0 (non-propagating scalar mode) -- the mathematical expression of tau = 0 for perturbations.

---

## 8. Key Papers

### BAO measurements
- **DESI 2024 VI** (arXiv:2404.03002): BAO from DESI DR1, 7 tracers, z = 0.1-4.2
- **DESI Lya BAO** (arXiv:2404.03001): H(2.33) and D_M(2.33) from Lyman-alpha forest
- **BOSS DR12** (arXiv:1607.03155): Previous BAO standard (superseded by DESI)

### GDM constraints relevant to BAO
- **Thomas, Kopp, Skordis 2016** (arXiv:1601.05097): |w_DM| < 2.4 x 10^{-3}, c_s^2 < 3.2 x 10^{-6} from Planck CMB (99.7% CL). These bounds apply at the CMB epoch; the Khronon satisfies them by factors of 10^4 - 10^7.

### Khronon/AeST theory
- **Skordis & Zlosnik 2021** (PRL 127, 161302; arXiv:2007.00082): AeST theory reproducing CDM-like CMB
- **Jacobson 2010** (arXiv:1001.4823): Khronon field as IR limit of Horava gravity
- **Blanchet & Skordis** (in prep): Khronon as GDM with c_s = 0

### DESI dark energy constraints
- DESI + CMB + PantheonPlus: w_0 = -0.727 +/- 0.067, w_a = -1.05 (+0.31/-0.27)
- DESI + CMB: 2.6sigma preference for evolving DE
- DESI + CMB + SN: 3.5-3.9sigma preference for evolving DE

---

## 9. Assessment Summary

| Test | Status | Detail |
|------|--------|--------|
| Sound horizon r_s | SAFE | Delta = 10^{-11}, undetectable |
| BAO peak position | SAFE | Identical to LCDM |
| BAO damping | SAFE | k_J >> k_BAO by 10^{10} |
| DESI distance measures | SAFE+ | chi^2 BETTER than LCDM by 4.72 |
| Alcock-Paczynski | SAFE | Delta F/F < 0.5% |
| CMB compatibility | SAFE | w_tilde << 1 at z=1100 |
| Late-time DM density | CAUTION | 20% higher than CDM at z=0 |

### OVERALL: SAFE

### Caveats requiring further investigation
1. **w_tilde(z=0) = 0.20**: The Khronon energy density today is 20% higher than CDM. While this improves BAO fits, it changes the late-time expansion and requires 5.2% less Lambda. This could affect Type Ia supernova constraints, weak lensing, and the growth rate f*sigma8.

2. **TKS 2016 bounds are CMB-era only**: The bound |w_DM| < 2.4 x 10^{-3} applies at z ~ 1100. At z < 1, the effective w_K ~ 0.05-0.20 is large. Low-z constraints from galaxy clustering and RSD need separate analysis.

3. **Covariance matrix needed**: The Delta chi^2 = -4.72 is computed without D_M-D_H correlations. A proper MCMC analysis with the full DESI covariance matrix is needed.

---

## 10. Surprising Finding: Khronon as Phantom Dark Energy Mimic

The most unexpected result: the Khronon model naturally produces a quintessence-like effective dark energy (w_eff ~ -0.95) without any dark energy modification. This arises purely from the dark matter sector having w_tilde(z=0) = 0.20, which requires less Lambda and produces a slightly different expansion history.

This is qualitatively in the same direction as the DESI dark energy signal (w > -1 at low z), though the Khronon effect is milder (w ~ -0.95 vs DESI's w_0 ~ -0.73). The Khronon prediction is a ZERO-parameter prediction (given Omega_c h^2 from CMB and the running mu_bg = H/c), making this a genuine test of the framework.

**If future DESI data confirm w_eff ~ -0.95 (rather than the current w_0 ~ -0.73), this would be strong evidence for the Khronon model.**
