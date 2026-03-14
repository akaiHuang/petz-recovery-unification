# CMB Boltzmann Code Verification Plan for mu_0 = H_0/c

**Date**: 2026-03-14
**Status**: Research complete, ready for implementation
**Purpose**: Determine how to verify the tau framework's CMB prediction using modified CLASS (GDM_CLASS)

---

## 1. Executive Summary

The tau framework claims that the Khronon field with mu_0 = H_0/c produces dust-like perturbations (w_tilde << 1) at CMB scales. Verification requires running a modified Boltzmann code (CLASS or CAMB) with the predicted GDM parameters and comparing the resulting C_l spectrum with Planck 2018 data.

**Key findings from this research:**

1. **gdm_class_public** (github.com/s-ilic/gdm_class_public) already supports k-dependent sound speed, but with the WRONG k-dependence (k^2 instead of 1/k^2). A minor code modification (~5 lines) is needed.

2. **~~The paper's numerical claim has an error~~** (**FIXED 2026-03-14**): Eq.(11) of paper3_dark_matter.tex previously claimed w_tilde(z~1100) ~ 6e-11; corrected to ~4e-10. The qualitative conclusion (w << 1) still holds.

3. **The naive 1/k^2 sound speed may already be excluded** by existing TKS2016 constraints (c_s^2 < 3.4e-6) at low k, because c_s^2(k=0.01) ~ 5e-4 >> 3.4e-6. This needs careful checking because TKS2016 used constant c_s^2, not k-dependent.

4. **The Khronon is NOT a simple GDM fluid** -- it has omega = 0 dispersion (non-propagating mode), which means the mapping to GDM parameters may be more subtle than c_s^2 = (mu_0/k)^2.

---

## 2. The GDM_CLASS Code

### 2.1 Repository
- **URL**: https://github.com/s-ilic/gdm_class_public
- **Authors**: Ilic, Kopp, Thomas, Skordis (arXiv:2004.09572)
- **Based on**: CLASS (Lesgourgues & Tram)
- **License**: MIT

### 2.2 GDM Parameters (from GDM.ini)

```ini
omega_gdm = 0.12                           # density (replaces omega_cdm)
type_gdm = time_only_bins                  # binned in scale factor
smooth_bins_gdm = yes                      # smooth transitions
time_transition_width_gdm = 8.             # transition width
time_values_gdm = 0.00001, 0.0001, 0.001, 0.01, 0.1  # bin edges (scale factor)
w_values_gdm = 0., 0., 0., 0., 0., 0.     # equation of state per bin
cs2_values_gdm = 0., 0., 0., 0., 0., 0.   # sound speed per bin
cv2_values_gdm = 0., 0., 0., 0., 0., 0.   # viscosity per bin
```

### 2.3 k-Dependent Sound Speed Support

GDM_CLASS has a built-in flag `k2_cs2_gdm`:

```ini
k2_cs2_gdm = yes   # enables k^2-dependent sound speed
```

**Implementation** (perturbations.c line 10442-10474):
```c
double cs2_gdm_of_a_and_k(...) {
    cs2 = twoD_pixel(pba, a, k, pba->cs2_values_gdm);  // time-dependent bin value

    if (pba->k2_cs2_gdm == _TRUE_) {
        double k_pivot = 0.01;  // hardcoded pivot in 1/Mpc
        double ca2_gdm = ppw->pvecback[pba->index_bg_ca2_gdm];
        cs2 = ca2_gdm + k*k/k_pivot/k_pivot * cs2;  // ← k^2 dependence
        if(cs2 > 1.) cs2 = 1.;
        if(cs2 < 0.) cs2 = 0.;
    }
    return cs2;
}
```

**Problem**: This gives cs2 ~ k^2 (INCREASING with k).
The tau framework predicts cs2 ~ 1/k^2 (DECREASING with k).

### 2.4 GDM Perturbation Equations (perturbations.c line 9173-9218)

```
delta' = -(1+w)(theta + h'/2) + 3aH[(w-ca2)*delta - Pi_nad]
theta' = -(1-3ca2)*aH*theta + k^2/(1+w)*[ca2*delta + Pi_nad] + metric_euler - k^2*shear
shear' = -3aH*shear + (8/3)*cv2/(1+w)*(theta + metric_shear)
```

where Pi_nad = (cs2 - ca2)*(delta + 3aH(1+w)*theta/k^2).

These are the STANDARD Hu (1998) GDM equations. The Khronon perturbation equation is different (has omega = 0 dispersion), but MAY map to these with appropriate parameter choices.

---

## 3. The tau Framework's Prediction

### 3.1 The Khronon GDM Parameters

From Blanchet-Skordis 2024 (arXiv:2404.06584):
- **w = 0** (dust-like in the relevant regime)
- **c_s^2 = 0** on Minkowski (omega = 0 dispersion)
- **c_s^2 ~ (mu_0/k)^2** on FRW (small correction from Hubble expansion)
- **c_vis^2 = 0** (no viscosity)

### 3.2 Numerical Values of c_s^2(k)

Using mu_0 = H_0/c = 2.248e-4 Mpc^{-1}:

| k (Mpc^-1) | c_s^2 = (mu_0/k)^2 | CMB relevance |
|-------------|---------------------|---------------|
| 0.0001 | 5.05 | Super-horizon |
| 0.001 | 0.0505 | Large scale ISW |
| 0.005 | 0.00202 | Low ell |
| 0.01 | 5.05e-4 | First peak |
| 0.05 | 2.02e-5 | Higher peaks |
| 0.1 | 5.05e-6 | High ell |
| 0.2 | 1.26e-6 | Damping tail |
| 0.5 | 2.02e-7 | Small scale |

### 3.3 Comparison with Existing Constraints

TKS2016 (arXiv:1605.00649) obtained from Planck 2015:
- c_s^2 < 3.4e-6 (95% CL) for **constant** c_s^2
- c_vis^2 < 3.3e-6 (95% CL) for **constant** c_vis^2

**Important caveat**: These constraints assume CONSTANT c_s^2. The tau prediction has c_s^2 ~ 1/k^2, which is large at low k but small at high k. The CMB is sensitive to a weighted average over k, so the effective constraint may be different.

At k = 0.01 (first peak): c_s^2 ~ 5e-4 >> 3.4e-6 -- appears excluded
At k = 0.1 (higher peaks): c_s^2 ~ 5e-6 -- marginally excluded
At k > 0.2: c_s^2 < 3e-6 -- within constraint

**This is a potential FALSIFICATION of the naive c_s^2 = (mu_0/k)^2 mapping.**
But the mapping itself is uncertain -- the Khronon has omega = 0 dispersion (non-propagating), which is fundamentally different from a GDM fluid with c_s^2 > 0. The correct translation may involve the Khronon being treated as having c_s^2 = 0 with a modified growth rate, not as c_s^2 > 0 with pressure support.

### 3.4 Error in Paper3 Eq.(11)

The paper claims w_tilde_0(z~1100) ~ 6e-11, using:
- mu_bg = H(z_rec)/c ~ 2.3e-14 s^{-1} (but this has the WRONG UNITS for the formula)
- k_phys ~ 0.01 Mpc^{-1}

**Correct computation**: With mu_0 = H_0/c = 2.248e-4 Mpc^{-1} and k_phys = k*(1+z) = 11 Mpc^{-1}:
- w_tilde = (mu_0/k_phys)^2 = (2.248e-4/11.01)^2 = 4.16e-10

The paper's 6e-11 is off by a factor of ~7 (likely a unit conversion error mixing s^{-1} with Mpc^{-1}). The qualitative conclusion (w_tilde << 1) still holds.

---

## 4. CAMB Alternative

### 4.1 CAMB Dark Energy Fluid
CAMB supports `DarkEnergyFluid` with constant c_s^2. It does NOT support:
- k-dependent c_s^2
- GDM (generalized dark matter)
- Custom dark matter sound speed

**Verdict**: CAMB is NOT suitable without significant Fortran modifications. GDM_CLASS is the better starting point.

---

## 5. Step-by-Step Verification Plan

### Level 0: Sanity Check (1-2 hours)
**Goal**: Confirm gdm_class reproduces LCDM when cs2 = 0.

1. Clone gdm_class_public
2. Compile: `make class` (requires gcc, no special dependencies)
3. Run with GDM.ini parameters:
   ```
   omega_gdm = 0.12, omega_cdm = 0 (or 1e-10)
   w_values_gdm = 0,...,0
   cs2_values_gdm = 0,...,0
   cv2_values_gdm = 0,...,0
   ```
4. Compare output C_l with standard CLASS (omega_cdm = 0.12)
5. **Expected**: Identical spectra (GDM with w=0, cs2=0, cvis2=0 IS CDM)

### Level 1: Constant c_s^2 Bracketing (2-4 hours)
**Goal**: Run constant c_s^2 values that bracket the prediction at different scales.

Run gdm_class with:
- cs2 = 0 (CDM baseline)
- cs2 = 1e-8 (well below constraint)
- cs2 = 1e-6 (near constraint boundary)
- cs2 = 1e-5 (above constraint)
- cs2 = 1e-4 (similar to prediction at k ~ 0.015)
- cs2 = 5e-4 (similar to prediction at first peak)

**Expected**: Increasing cs2 suppresses power at small scales (high ell), shifts peak positions.

Compare each with Planck 2018 TT power spectrum to see at what cs2 the spectrum becomes incompatible.

### Level 2: k-Dependent Sound Speed (1-2 days)
**Goal**: Implement c_s^2 ~ 1/k^2 in gdm_class.

**Code modification**: In `/source/perturbations.c`, add a new flag `kinv2_cs2_gdm`:

```c
// After line 10470 in cs2_gdm_of_a_and_k():
if (pba->kinv2_cs2_gdm == _TRUE_) {
    double k_pivot = 0.01;  // 1/Mpc
    double ca2_gdm = ppw->pvecback[pba->index_bg_ca2_gdm];
    cs2 = ca2_gdm + cs2 * k_pivot*k_pivot / (k*k);
    if(cs2 > 1.) cs2 = 1.;
    if(cs2 < 0.) cs2 = 0.;
}
```

Changes needed:
1. `/include/background.h`: Add `short kinv2_cs2_gdm;`
2. `/source/input.c`: Add parser for `kinv2_cs2_gdm`
3. `/source/input.c`: Add default `pba->kinv2_cs2_gdm = _FALSE_;`
4. `/source/perturbations.c`: Add the code block above

Then run with:
```ini
omega_gdm = 0.12
kinv2_cs2_gdm = yes
cs2_values_gdm = 5.05e-4,...  # = mu_0^2/k_pivot^2 = (2.248e-4/0.01)^2
```

This gives cs2(k) = 5.05e-4 * (0.01/k)^2 = mu_0^2/k^2, exactly the tau prediction.

### Level 3: Time-Dependent w(a) + k-Dependent c_s^2 (3-5 days)
**Goal**: Implement the full Khronon GDM prediction.

The Khronon background has a stiff-to-dust transition at:
- w(a) transitions from ~1 (stiff, a << a_transition) to ~0 (dust, a >> a_transition)
- The transition epoch depends on w_tilde_0 = delta_0/2

For mu >> H_0/c (the running case): the transition happens very early, so w = 0 throughout the CMB epoch.

Use time bins:
```ini
time_values_gdm = 1e-6, 1e-5, 1e-4, 1e-3, 0.01, 0.1
w_values_gdm = 0, 0, 0, 0, 0, 0, 0     # dust throughout
cs2_values_gdm = [computed from mu_bg(a)]  # time-dependent
```

The time-dependence of cs2 adds:
- cs2(a, k) = (mu_bg(a))^2 / k^2 = (H(a)/c)^2 / k^2

This requires implementing a custom function, not just bin values.

### Level 4: Full Khronon Perturbation Equations (1-2 weeks)
**Goal**: Implement the actual Khronon scalar field equation in CLASS.

The Khronon perturbation is NOT a simple GDM fluid. It has:
- omega = 0 dispersion relation on Minkowski
- A constraint (unit timelike vector) that kills the propagating mode
- Coupling to metric perturbations through the modified Einstein equations

This would require:
1. Adding the Khronon field as a new species in CLASS (not GDM)
2. Implementing the linearized scalar field equation from Blanchet-Skordis 2024
3. Modifying the Einstein constraint equations
4. Implementing the modified Friedmann equation for the background

This is essentially what Skordis-Zlosnik did for AeST, and it took years of development.

---

## 6. Effort Estimates

| Level | Description | Effort | Prerequisite |
|-------|-------------|--------|-------------|
| 0 | Sanity check (GDM with cs2=0 = CDM) | 1-2 hours | Compile gdm_class |
| 1 | Constant c_s^2 bracketing | 2-4 hours | Level 0 |
| 2 | k-dependent c_s^2 ~ 1/k^2 | 1-2 days | Level 1 + C coding |
| 3 | Time-dependent + k-dependent | 3-5 days | Level 2 |
| 4 | Full Khronon equations | 1-2 weeks | Expert-level CLASS knowledge |

### Minimum Viable Verification: Levels 0 + 1 (half a day)

Running GDM_CLASS with constant c_s^2 at several values immediately tells us:
1. Whether c_s^2 = 5e-4 (the prediction at first peak) is compatible with Planck
2. At what c_s^2 the spectrum becomes visibly different from CDM
3. Whether the TKS2016 bound c_s^2 < 3.4e-6 indeed excludes the naive prediction

**If cs2 = 5e-4 is excluded by the CMB (very likely based on TKS2016), this immediately constrains the tau framework: the effective GDM c_s^2 must be < 3.4e-6 at all CMB-relevant scales, meaning the Khronon cannot be a simple fluid with c_s^2 = (mu_0/k)^2.**

---

## 7. Critical Questions Revealed by This Analysis

### 7.1 Is c_s^2 = (mu_0/k)^2 the correct GDM mapping?

The Khronon has omega = 0 dispersion (non-propagating mode). This is qualitatively different from a fluid with c_s^2 > 0 and Jeans-like oscillation. The correct mapping may be:
- c_s^2 = 0 (the mode doesn't propagate)
- Modified growth rate (the mode grows differently from CDM)

In AeST (Skordis-Zlosnik 2021), the scalar field IS treated as having c_s = 0 on the FRW background. The dust-like behavior comes from the constraint structure, not from a small but nonzero sound speed. If the Khronon works the same way, then c_s^2 = 0 is the correct GDM parameter, and the verification at Level 0 (which recovers CDM exactly) would be the appropriate test.

### 7.2 What does w_tilde actually control?

The Khronon has two contributions to its energy density:
- Dust-like (rho ~ a^{-3}): dominates when w_tilde << 1
- Stiff-like (rho ~ a^{-6}): dominates when w_tilde >> 1/a^3

The w_tilde parameter controls the BACKGROUND equation of state (the ratio of stiff to dust energy), NOT the perturbation sound speed. The equation w_tilde = mu_bg^2/k_phys^2 in the paper may refer to the effective equation of state for a specific PERTURBATION MODE, not the GDM c_s^2.

### 7.3 The Catch-22 and its resolution

Without running (mu = mu_0 = H_0/c):
- Omega_DM ~ 0.26 requires delta_0 ~ 0.34
- Dust-like at z=1100 requires delta_0 << 1e-9
- INCOMPATIBLE (by 8 orders of magnitude)

With running mu_bg(a) = H(a)/c:
- The running increases the effective mu at early times
- This suppresses w_tilde at high redshift
- BUT: the computation shows w_tilde ~ O(0.5) even with running (see Section 3.4 of omega_dm_prediction.md)

**The running mu_bg(a) = H(a)/c alone does NOT solve the Catch-22.** The paper's claim relies on the scale-dependent running mu(k) = k, but this enters at the perturbation level, not the background level. The background Catch-22 remains.

### 7.4 What would a positive verification look like?

If the Khronon truly has c_s^2 = 0 (as expected from omega = 0 dispersion):
- GDM with w = 0, cs2 = 0 IS CDM
- The C_l spectrum reproduces Planck exactly
- No new prediction; the theory is compatible but not falsifiable at this level

If the Khronon has c_s^2 = (mu_0/k)^2:
- GDM with k-dependent cs2 gives a DIFFERENT C_l from CDM
- This is potentially falsifiable (and likely falsified at low k)
- Interesting only if it reproduces Planck better than CDM (unlikely)

The most valuable outcome: running Level 1 to determine the exact c_s^2 threshold where the CMB becomes incompatible, then comparing with the tau prediction.

---

## 8. Immediate Next Steps

1. ~~**Fix paper3 eq.(11)**~~: **DONE (2026-03-14)** — corrected to ~4e-10 with proper units (Mpc^{-1} throughout).

2. **Compile and run gdm_class** (Level 0):
   ```bash
   cd /tmp/gdm_class_public
   make class
   ./class GDM.ini
   ```
   Compare output with standard CLASS.

3. **Run Level 1 bracketing**: Modify GDM.ini with several constant cs2 values, generate C_l spectra, compare with Planck.

4. **Clarify the GDM mapping**: Read the full text of Blanchet-Skordis 2024 (arXiv:2404.06584) Section 4.2 to determine whether c_s^2 = 0 or c_s^2 = (mu/k)^2 is the correct perturbation sound speed for the Khronon on FRW.

5. **If c_s^2 = 0**: The verification is trivial (GDM with cs2=0 IS CDM), but the interesting question becomes whether the BACKGROUND w_tilde is consistent (the Catch-22 problem).

6. **If c_s^2 > 0**: Implement Level 2 (1/k^2 sound speed) and run the comparison with Planck.

---

## Appendix: GDM_CLASS Quick Start

```bash
# Clone
git clone https://github.com/s-ilic/gdm_class_public
cd gdm_class_public

# Compile (requires gcc)
make class

# Run with default GDM parameters
./class GDM.ini

# Output in output/ directory:
# output/kindPinad_cl.dat          # C_l spectrum
# output/kindPinad_cl_lensed.dat   # Lensed C_l
# output/kindPinad_background.dat  # Background evolution

# Python wrapper:
make   # compiles both C code and Python wrapper
python -c "from classy import Class; M = Class(); M.set({'omega_gdm':0.12}); M.compute()"
```

### Key files to modify for custom c_s^2(k):
- `/source/perturbations.c` line 10442-10474: `cs2_gdm_of_a_and_k()` function
- `/include/background.h` line 79: k2_cs2_gdm flag declaration
- `/source/input.c` line 982-989: parser for k2_cs2_gdm
- `/GDM.ini`: configuration file
