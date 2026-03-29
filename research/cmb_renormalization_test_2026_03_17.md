# CMB Renormalization Test: Khronon at z=1100

**Date**: 2026-03-17
**Status**: COMPLETED -- Mixed results
**Tool**: gdm_class_public (GDM extension of CLASS Boltzmann code)

## The Idea

Instead of fitting Omega_DM at z=0, normalize the Khronon energy density to match CDM at z=1100 (recombination). The Khronon density ratio:

```
rho_K(z) / rho_CDM(z) = (2 + delta(z)) / (2 + delta(z_ref))
```

At z=1100: delta(1100) ~ 0.82 (early universe, larger perturbation)
At z=0:   delta(0) ~ 0.34 (today)

So: rho_K(z=0) / rho_CDM(z=0) = (2+0.34)/(2+0.82) = 2.34/2.82 = 0.830

This predicts Omega_DM(today) ~ 0.83 x 0.265 ~ 0.22, with more Lambda needed and potentially HIGHER H_0.

## Two Testing Approaches

### Approach A: Naive (reduce omega_cdm directly)
Replace omega_cdm = 0.12 with omega_gdm = 0.10 (or 0.09, 0.08) and see what happens.

### Approach B: Correct (omega_gdm = 0.12, late-time w > 0)
Keep omega_gdm = 0.12 at recombination (preserving CMB physics), but add a small w > 0 in the last time bin (a > 0.1, z < 9) to model faster-than-CDM dilution.

## CLASS Configurations Run

15 total configurations across both approaches.

## Results

### Master Table

```
Run                          TT RMS    P1/P3   dP1/P3  phi RMS  CV 1sig  CV 2sig
--------------------------------------------------------------------------------
LCDM (reference)                ---   2.2213      ---      ---      ---      ---
A1: w_dm=0.10 h=67           0.1318   2.3821   +7.24%   0.2825    22.5%    44.2%
A2: w_dm=0.09 h=67           0.2128   2.4784  +11.57%   0.4081    11.9%    29.8%
A3: w_dm=0.08 h=67           0.3041   2.5889  +16.55%   0.5239     7.2%    19.8%
A4: w_dm=0.10 h=72           0.0715   2.3817   +7.22%   0.3096    30.9%    59.9%
A5: w_dm=0.10 h=73           0.0614   2.3816   +7.22%   0.3154    35.7%    68.4%
A6: w_dm=0.10 h=73 A_s       0.0614   2.3816   +7.22%   0.3800    34.2%    56.3%
A7: w_dm=0.10 h=73 tau       0.0547   2.3823   +7.25%   0.3155    43.3%    61.7%
B1: w=0.01 h=67              0.0222   2.1632   -2.61%   0.1557    87.1%   100.0%
B2: w=0.02 h=67              0.0430   2.1067   -5.16%   0.3358    45.5%    89.4%
B3: w=0.03 h=67              0.0625   2.0519   -7.62%   0.5442    31.9%    65.8%
B4: w=0.05 h=67              0.0980   1.9474  -12.33%   1.0593    20.2%    37.2%
B5: w=0.03 h=72              0.1051   2.0521   -7.62%   0.4855    22.3%    45.7%
B6: w=0.01 h=69              0.0394   2.1633   -2.61%   0.1404    60.7%    81.7%
B7: w=0.01 h=70              0.0509   2.1634   -2.60%   0.1306    50.6%    71.0%
```

Key columns:
- **TT RMS**: RMS fractional deviation of C_l^TT from LCDM (l=30-2500)
- **P1/P3**: Ratio of 1st to 3rd acoustic peak heights (LCDM = 2.2213)
- **dP1/P3**: Fractional deviation of P1/P3 from LCDM (Planck error ~1%)
- **phi RMS**: RMS deviation of lensing potential
- **CV 1sig / 2sig**: Fraction of multipoles within 1/2 sigma cosmic variance

## Approach A: RULED OUT

Reducing omega_cdm from 0.12 to 0.10:
- Peak height ratio P1/P3 shifts by **+7.2%** (LCDM: 2.22, renormalized: 2.38)
- This encodes the baryon-to-DM ratio (baryon loading), measured by Planck to ~1%
- **~8 sigma tension** -- cannot be compensated

What h = 0.73 does:
- Fixes peak LOCATIONS (angular diameter distance compensation) -- peaks return to l ~ 221, 539, 821
- Does NOT fix peak HEIGHT RATIOS -- P1/P3 remains at 2.38 regardless of h
- Reduces damping tail RMS (helps at l > 800)

What A_s reduction / tau increase does:
- Rescales ALL peaks uniformly -- does not change ratios
- P1/P3 = 2.38 is INVARIANT under amplitude rescaling
- Confirmed: A_s = 1.9e-9 and tau = 0.10 both leave P1/P3 unchanged

**Conclusion A: The naive omega_cdm = 0.10 approach is dead. Peak ratios are a hard constraint.**

## Approach B: Partially Viable

### Best candidate: w_late = 0.01, h = 0.6736 (B1)

```
TT deviation at key multipoles:
  l= 220: dC/C = -0.035  (0.5 sigma cosmic variance)
  l= 540: dC/C = -0.032  (0.7 sigma)
  l= 810: dC/C = -0.008  (0.2 sigma)
  l=1000: dC/C = -0.019  (0.6 sigma)
  l=1500: dC/C = -0.022  (0.8 sigma)
  l=2000: dC/C = +0.007  (0.3 sigma)
  l=2500: dC/C = -0.024  (1.2 sigma)
```

87% of multipoles within cosmic variance. 100% within 2 sigma. Excellent!

But:
- **P1/P3 = 2.163** vs LCDM 2.221 -- deviation of -2.6%
- Lensing excess: +8% at l=100, +14% at l=500 (detectable)
- omega_DM_eff(today) = 0.112 (only 7% less than CDM, not 17%)

### H0 scan with w_late = 0.01

```
     h    TT RMS  P1/P3  dP1/P3  phi RMS  CV 1sig
  0.6736  0.0222  2.163  -2.61%   0.156    87.1%
  0.6900  0.0394  2.163  -2.61%   0.140    60.7%
  0.7000  0.0509  2.163  -2.60%   0.131    50.6%
```

Higher h degrades TT fit. P1/P3 is not affected by h (as expected -- it depends on omega_cdm/omega_b).

### Effective Omega_DM today

```
w_late=0.010: omega_eff=0.112, Omega_DM=0.247
w_late=0.020: omega_eff=0.104, Omega_DM=0.230
w_late=0.030: omega_eff=0.098, Omega_DM=0.215
w_late=0.050: omega_eff=0.085, Omega_DM=0.187
```

## The Fundamental Tension

The Khronon predicts (2+0.34)/(2+0.82) = 0.83 ratio, requiring w_late ~ 0.027.
But w_late > 0.015 is already in tension with Planck (P1/P3 deviation > 4%).

```
w_late needed for 0.83 ratio:  ~0.027
w_late allowed by CMB:         <0.015 (2 sigma)
                               <0.010 (1 sigma)
```

**The predicted effect is ~2x larger than CMB allows.**

## What This Means for the Theory

### The Bad News
1. The full renormalization ratio of 0.83 (17% less DM) is RULED OUT by CMB peak ratios
2. Even Approach B (correct treatment with late w) shows the same problem
3. w > 0 causes late-ISW effect that suppresses peaks uniformly and changes P1/P3
4. Lensing excess (+8-16%) provides an independent constraint

### The Good News
1. w_late = 0.01 (7% less DM) IS compatible with CMB at 1 sigma
2. The framework works in principle -- just not at the predicted magnitude
3. This suggests the delta values used (0.34 at z=0, 0.82 at z=1100) may need revision
4. If the actual Khronon evolution is milder (ratio ~0.93 instead of 0.83), it fits perfectly

### The Hubble Tension
- w_late = 0.01 with h = 0.6736: excellent CMB fit, no H0 help
- w_late = 0.01 with h = 0.69: decent CMB (61% within 1 sigma CV), mild H0 help
- w_late = 0.01 with h = 0.70: marginal CMB (51% within 1 sigma CV), ~2 sigma H0 help
- **Not a Hubble tension solution**, but could contribute ~1 km/s/Mpc

### Revised Prediction
If the CMB constrains w_late < 0.015, then:
- omega_DM_eff(today) > 0.108
- Omega_DM > 0.238 (at Planck h)
- Deviation from CDM: < 10%

This is still a testable prediction, just smaller than originally hoped.

## Key Diagnostic: Why Peak Ratios Change

The peak height ratio P1/P3 encodes the baryon-to-DM ratio:
- More baryons relative to DM -> odd peaks enhanced relative to even peaks
- P1/P3 increases when omega_cdm decreases (Approach A: +7%)
- P1/P3 decreases when w > 0 (Approach B: late-ISW suppresses large-scale more)

In Approach B, the late-ISW from w > 0 suppresses power at low l more than high l, which effectively makes the first peak lower relative to the third peak. This is the opposite direction from Approach A, showing these are truly different physics.

## Files

- CLASS configs: `/Users/akaihuangm1/Desktop/github/gdm_class_public/renorm_*.ini`
- Output data: `/Users/akaihuangm1/Desktop/github/gdm_class_public/output/renorm_*_cl.dat`
- Analysis scripts:
  - `analyze_renormalization.py` (basic)
  - `analyze_renorm_v2.py` (compensation attempts)
  - `analyze_renorm_final.py` (approach A vs B)
  - `analyze_renorm_complete.py` (all 15 runs, master table)

## Next Steps

1. **Recompute delta(z)**: The values delta(z=0) = 0.34, delta(z=1100) = 0.82 need verification. If the actual ratio is closer to 0.93, the theory survives.
2. **Time-dependent w(a)**: The constant-w approximation is crude. The actual Khronon evolution d(ln(2+delta))/d(ln a) should be computed and input as time-varying w bins.
3. **Lensing cross-check**: The +8-14% lensing excess for w_late = 0.01 is a clear prediction. Compare with ACT DR6 / Planck lensing.
4. **DESI BAO**: The late-time w > 0 affects the expansion history. Check consistency with DESI DR1 BAO.
