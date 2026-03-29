# Feasibility Report: Full Boltzmann Solver CMB Computation for Khronon Dark Matter

**Author**: Sheng-Kai Huang (with computational analysis)
**Date**: 2026-03-16
**Status**: Complete feasibility assessment
**Purpose**: Determine practical requirements for running a full C_l computation using a modified Boltzmann solver, compare with Planck 2018, and assess what the results would prove or disprove.

---

## 1. Recommended Approach: CLASS via gdm_class_public

### 1.1 Why CLASS over CAMB

| Criterion | CLASS | CAMB | Winner |
|-----------|-------|------|--------|
| Language | C (core) + Python wrapper | Fortran + Python wrapper | CLASS (C easier to modify) |
| GDM support | gdm_class_public exists | No GDM module | CLASS |
| k-dependent c_s^2 | Built-in flag (k^2); needs ~5 lines for 1/k^2 | Not supported | CLASS |
| Time-dependent w(a) | Binned support built-in | DarkEnergyFluid (constant c_s^2 only) | CLASS |
| Community forks | Many (axionCAMB, MGCAMB, class.VFDM, etc.) | Some | Tie |
| Speed | ~2.5x faster than CAMB for LCDM | Slower | CLASS |
| Documentation | Excellent (Lesgourgues lectures, source comments) | Good | CLASS |
| Python wrapper | `classy` module, pip-installable after compile | `camb` pip package | Tie |
| Differentiability | No (but SymBoltz.jl offers this in Julia) | No | Tie |

**Verdict**: CLASS via `gdm_class_public` is the clear winner. It already supports GDM with binned w(a), c_s^2(a), c_vis^2(a) and has a k-dependent sound speed flag that needs only a minor modification (~5 lines of C code) to implement c_s^2 ~ 1/k^2 instead of the built-in k^2 dependence.

### 1.2 Alternative: SymBoltz.jl

A newer option is [SymBoltz.jl](https://github.com/hersle/SymBoltz.jl) (Sletmoen 2025, arXiv:2509.24740), a symbolic-numeric Julia package that:
- Allows models to be defined from symbolic equations
- Compiles to efficient numerical code
- Supports automatic differentiation for parameter derivatives
- Is approximation-free (implicit stiff solvers)

**Advantage**: Custom species can be defined symbolically without modifying C source code. The Khronon perturbation equations could be entered directly.

**Disadvantage**: Newer, less tested, smaller community. No GDM module exists. Would require writing the full Khronon perturbation system from scratch.

**Recommendation**: Start with gdm_class_public (proven, well-tested). Consider SymBoltz.jl for Level 4 (full Khronon equations) if the symbolic interface proves more productive than modifying CLASS C code.

### 1.3 Repository and Dependencies

**gdm_class_public**:
- URL: https://github.com/s-ilic/gdm_class_public
- Authors: Ilic, Kopp, Thomas, Skordis (arXiv:2004.09572)
- License: MIT
- Dependencies: gcc (C compiler), Python 3 (for wrapper), numpy, scipy, cython
- Based on: CLASS v2.x (Lesgourgues & Tram)

**System requirements**:
- macOS: Xcode command line tools (for gcc/clang)
- No GPU needed
- No MPI needed for single runs
- ~200 MB disk space for code + output

---

## 2. What Modifications Are Needed

The Khronon framework has a layered structure. Different levels of modification test different aspects of the theory. The key insight from existing research (cmb_verification_result.md, cs2_verification_supplement.md) is:

**The Khronon at CMB scales is effectively CDM** (w ~ 0, c_s^2 ~ 0), with deviations controlled by:
- Background: w(z) = (delta_0/2) x (H_0/H(z))^2 ~ 3 x 10^{-10} at z=1100
- Perturbation: c_s^2 = 0 (omega = 0 dispersion on Minkowski)
- The f(delta) = (4+3delta)/4 correction to rho_K modifies H(z) at late times

### 2.1 Level 0: Sanity Check (GDM = CDM baseline)

**Modification**: None (use gdm_class_public as-is).

**Configuration** (GDM.ini):
```ini
# Replace CDM with GDM having identical parameters
omega_cdm = 0        # or 1e-15
omega_gdm = 0.1200   # Planck 2018 best-fit

# GDM parameters: dust (= CDM)
type_gdm = time_only_bins
smooth_bins_gdm = yes
time_transition_width_gdm = 8.
time_values_gdm = 0.00001, 0.0001, 0.001, 0.01, 0.1
w_values_gdm = 0., 0., 0., 0., 0., 0.
cs2_values_gdm = 0., 0., 0., 0., 0., 0.
cv2_values_gdm = 0., 0., 0., 0., 0., 0.

# Standard cosmological parameters (Planck 2018)
h = 0.6736
omega_b = 0.02237
T_cmb = 2.7255
A_s = 2.1e-9
n_s = 0.9649
tau_reio = 0.0544

# Output
output = tCl, pCl, lCl, mPk
l_max_scalars = 2500
```

**Expected result**: C_l spectrum identical to standard CLASS with omega_cdm = 0.12.

**Purpose**: Confirms the code compiles and runs correctly, and that GDM(w=0, c_s^2=0) = CDM.

### 2.2 Level 1: Constant c_s^2 Bracketing

**Modification**: None (just change parameter values in GDM.ini).

**Runs**:
```
cs2 = 0         (CDM baseline)
cs2 = 1e-8      (well below TKS bound)
cs2 = 1e-7      (below TKS bound)
cs2 = 1e-6      (near TKS bound: 3.4e-6)
cs2 = 5e-6      (above TKS bound)
cs2 = 1e-5      (above TKS bound)
cs2 = 1e-4      (Khronon prediction at k ~ 0.015 IF c_s^2 = (mu_0/k)^2)
cs2 = 5e-4      (Khronon prediction at first peak k ~ 0.01)
```

**Purpose**: Determines the exact c_s^2 threshold where the CMB spectrum becomes incompatible with Planck. This directly tests whether the naive mapping c_s^2 = (mu_0/k)^2 is viable.

**Key test**: Compare each C_l^TT with Planck 2018 binned data. Compute chi^2 residuals.

### 2.3 Level 2: The f(delta) Background Correction

**Modification**: Custom background H(z) incorporating the Khronon energy density correction.

The Khronon energy density is:
```
rho_K(a) = (c^2 I_0 / (16 pi G)) * f(delta(a)) / a^3
```
where:
```
delta(a) = I_0 / (2 mu_bg(a)^2)
mu_bg(a) = H(a)/c           (running coupling)
f(delta) = (4 + 3 delta) / 4
I_0 = delta_0 * 2 * mu_0^2  (initial condition, delta_0 = 0.340)
```

At z=0: delta_0 = 0.340, f = 1.255, so rho_K(today) = 1.255 * rho_CDM
This means 25.5% MORE dark matter energy density than in standard CDM!

To maintain Omega_total = 1, this requires:
- Less Lambda: Omega_Lambda decreases by ~5.2%
- Different H(z): late-time expansion history changes

**Implementation in gdm_class_public**:

The simplest approach is to use the time-binned w(a) feature:

```ini
# The effective w(a) from f(delta) correction
# w_eff(a) = p_K / rho_K ≈ delta(a) / 2 for small delta

time_values_gdm = 1e-6, 1e-5, 1e-4, 1e-3, 0.01, 0.03, 0.1, 0.3, 0.5, 0.8
w_values_gdm = 0, 0, 0, 0, 0, 2e-5, 3e-4, 0.01, 0.04, 0.12, 0.17
```

Actually, the more rigorous approach is to modify the background module. The key file is `/source/background.c`, specifically the function `background_functions()` where the dark matter density is computed. We would change:

```c
// Standard CDM: rho_cdm = rho_cdm_0 / a^3
// Khronon: rho_K = rho_K_0 * f(delta(a)) / (a^3 * f(delta_0))
// where delta(a) = delta_0 * (H_0/H(a))^2

double delta_a = delta_0 * pow(H_0_over_Ha, 2);  // H_0/H(a) squared
double f_a = (4.0 + 3.0 * delta_a) / 4.0;
double f_0 = (4.0 + 3.0 * delta_0) / 4.0;
pvecback[pba->index_bg_rho_gdm] = rho_gdm_0 * f_a / (a3 * f_0);
```

**Challenge**: H(a) itself depends on rho_K(a), creating a self-consistency loop. This must be solved iteratively or via a shooting method. CLASS already handles this for dark energy models, so the infrastructure exists.

**Simpler alternative**: Since delta(a) = delta_0 * (H_0/H(a))^2 and H(a) is dominated by radiation at early times, we can pre-compute delta(a) using the LCDM H(a) as a first approximation, then iterate.

### 2.4 Level 3: Time-Dependent w(a) + k-Dependent c_s^2(a,k)

**Modification**: Implement the full GDM mapping of the Khronon.

From Blanchet-Skordis 2024 (arXiv:2404.06584), the Khronon maps to GDM with:
```
w(a) = delta(a) / 2 = (delta_0 / 2) * (H_0 / H(a))^2
c_s^2(a, k) = 0    (leading order, omega = 0 dispersion)
c_vis^2 = 0
```

The c_s^2 = 0 at leading order is the crucial result. The subleading correction from FRW expansion is:
```
c_s^2(a, k) = (mu_bg(a) / k)^2 = (H(a) / (ck))^2
```

For CMB-relevant k and z=1100:
- k = 0.05 Mpc^{-1}: c_s^2 ~ (2.35e4 * 2.25e-4 / 0.05)^2 ~ (0.106)^2 ~ 0.011
  WAIT: This uses mu_bg = H(z)/c at z=1100 with H(z=1100) ~ 2.35e4 * H_0

Actually, from the cs2_verification_supplement.md, the correct computation is:
```
c_s^2(z, k) = w(z) = (delta_0/2) * (H_0/H(z))^2  [BACKGROUND sound speed]
```
NOT c_s^2 = (mu/k)^2.

The distinction is critical:
- Background EoS: w(z=1100) ~ 3e-10 (safe, from running mu)
- Perturbation sound speed: c_s^2 = 0 (omega = 0 dispersion on Minkowski)
- The effective sound speed from Jeans analysis (BS2024 Eq. 4.31): suppressed at CMB scales

**Code changes needed** (in perturbations.c):
```c
// In cs2_gdm_of_a_and_k():
// For Khronon: c_s^2 = 0 at all scales (leading order)
// Subleading: c_s^2 = c_ad^2 / (1 + c_ad^2 * c^2 * k^2 / (4*pi*G*a^2*rho_K*(1+w)))
// At CMB scales k << k_J, this reduces to c_s^2 ≈ c_ad^2 = w(z) ~ 3e-10
cs2 = w_gdm;  // adiabatic sound speed = w for dust-like fluid
```

### 2.5 Level 4: Full Khronon Perturbation Equations (Expert Level)

**Modification**: Implement the actual Khronon scalar field as a new species in CLASS.

From Blanchet-Skordis 2024, the linearized Khronon equation on FRW is:
```
delta_K'' + 2 H delta_K' = 4 pi G a^2 sum_i rho_i delta_i   (schematically)
```

where delta_K is the Khronon density contrast, and the source includes ALL species (including the Khronon itself). The key features:
- No pressure term (omega = 0 dispersion => no k^2 delta_K term)
- Couples to metric perturbations through modified Poisson equation
- The Khronon modifies the Poisson equation: k^2 Phi = -4 pi G a^2 (rho_b delta_b + rho_gamma delta_gamma + rho_nu delta_nu + rho_K delta_K)

**Files to modify**:
1. `/include/perturbations.h`: Add Khronon perturbation variables
2. `/source/perturbations.c`: Add Khronon evolution equations in `perturbations_derivs()`
3. `/source/perturbations.c`: Add Khronon initial conditions in `perturbations_initial_conditions()`
4. `/source/background.c`: Add Khronon background evolution
5. `/source/thermodynamics.c`: Verify Khronon does not couple to photons/baryons
6. `/source/input.c`: Add Khronon parameters

**This is essentially what Skordis-Zlosnik did for AeST.** Their code is NOT publicly available. The closest public implementation is `class.VFDM` (ultralight vector field dark matter, arXiv:2408.12052), which added a custom vector field species to CLASS. This could serve as a template.

**Estimated effort**: 1-3 weeks for a physicist comfortable with CLASS internals.

---

## 3. Expected Runtime and Resources

### 3.1 Hardware Requirements

| Level | CPU | RAM | Disk | GPU |
|-------|-----|-----|------|-----|
| 0-1 | Any modern laptop | 1 GB | 500 MB | Not needed |
| 2-3 | Any modern laptop | 2 GB | 1 GB | Not needed |
| 4 | Any modern laptop | 4 GB | 2 GB | Not needed |
| MCMC | Multi-core desktop or cluster | 8+ GB | 10+ GB | Not needed |

**A MacBook (M1/M2/M3) is more than sufficient for all levels.**

### 3.2 Computation Time Per Run

CLASS is highly optimized. Typical run times on a modern laptop:

| Computation | l_max | Time | Notes |
|-------------|-------|------|-------|
| C_l^TT only | 2500 | ~3-5 seconds | Standard LCDM |
| C_l^TT + TE + EE | 2500 | ~5-10 seconds | With polarization |
| C_l + lensing | 2500 | ~10-15 seconds | Full lensed spectra |
| C_l + P(k) | 2500 | ~15-20 seconds | Plus matter power spectrum |
| GDM (constant params) | 2500 | ~5-15 seconds | Similar to LCDM |
| GDM (time-binned) | 2500 | ~10-30 seconds | Interpolation overhead |
| Modified background | 2500 | ~10-30 seconds | Iterative H(z) |

### 3.3 Total Time Estimates by Level

| Level | Description | Runs needed | Total compute | Human effort |
|-------|-------------|-------------|---------------|--------------|
| 0 | Sanity check | 2 | ~30 seconds | 1-2 hours (compile + run) |
| 1 | Constant c_s^2 bracketing | 8 | ~2 minutes | 2-4 hours |
| 2 | f(delta) background | 10 | ~5 minutes | 1-2 days (C coding) |
| 3 | Time + k dependent | 20 | ~10 minutes | 3-5 days (C coding) |
| 4 | Full Khronon equations | 50+ | ~30 minutes | 1-3 weeks |
| MCMC | Parameter estimation | 10^4-10^5 | 1-7 days | 1-2 weeks |

**The computational bottleneck is NEVER the run time. It is the human effort to implement and debug the code modifications.**

### 3.4 Comparison Data

Planck 2018 data products are freely available:
- Planck Legacy Archive: https://pla.esac.esa.int/
- Binned C_l^TT, C_l^TE, C_l^EE with error bars
- Full likelihood code (Plik, Commander, etc.)
- Best-fit LCDM parameters

For a first comparison, the binned power spectrum (30 points) is sufficient. For a rigorous chi^2, the full likelihood code is needed (available as a CLASS/CAMB plugin via `clik`/`plc`).

---

## 4. What the Output Would Tell Us (Pass/Fail Criteria)

### 4.1 The Three Scenarios

**Scenario A: c_s^2 = 0 (Khronon IS CDM at perturbation level)**

If the Khronon truly has c_s^2 = 0 (omega = 0 dispersion), then:
- GDM(w=0, c_s^2=0, c_vis^2=0) IS CDM by construction
- The C_l spectrum will be identical to LCDM
- **The theory is compatible but not yet distinguishable from CDM**
- The distinguishing signal comes from the f(delta) background correction (Level 2)

**Pass criteria**: C_l residuals < 0.1% relative to LCDM (within cosmic variance)
**Fail criteria**: Cannot fail in this scenario (it is CDM by definition)

**What we learn**: The Khronon framework passes the CMB test trivially. The interesting physics is in the BACKGROUND correction (f(delta) effect) and the late-time ISW.

**Scenario B: c_s^2 = w(z) ~ 3 x 10^{-10} (adiabatic sound speed from running mu)**

If the effective perturbation sound speed equals the adiabatic sound speed:
- At z=1100: c_s^2 ~ 3e-10 (11,000x below TKS bound of 3.4e-6)
- This is indistinguishable from c_s^2 = 0 at Planck precision
- Even CMB-S4 (sigma(c_s^2) ~ 10^{-7}) would not distinguish this from zero

**Pass criteria**: Same as Scenario A (effectively identical to CDM)
**What we learn**: The running mu successfully suppresses the sound speed to undetectable levels.

**Scenario C: c_s^2 = (mu_0/k)^2 (naive mapping, NO running)**

If someone incorrectly uses the FIXED mu_0 = H_0/c (without running):
- At k=0.01 (first peak): c_s^2 ~ 5e-4 >> 3.4e-6 (TKS bound)
- This would be EXCLUDED by Planck at high significance

**Pass criteria**: This scenario SHOULD fail (it represents the incorrect fixed-mu case)
**What we learn**: Confirms that the running mu_bg(z) = H(z)/c is ESSENTIAL for CMB compatibility.

### 4.2 The Key Diagnostic: f(delta) Background Effect

The most scientifically interesting output is from Level 2 (the background correction):

```
rho_K(z=0) = 1.255 * rho_CDM  (25.5% more DM energy today)
Omega_Lambda must decrease by ~5.2% to maintain flatness
Effective late-time expansion: w_eff ~ -0.95 (quintessence-like)
```

This changes:
1. **Late-time ISW effect**: Enhanced ISW at low ell (l < 30) due to faster potential decay
2. **Distance to last scattering**: Slightly different D_A (shifts peak positions by ~0.1%)
3. **BAO distances**: Better fit to DESI DR1 data (Delta chi^2 = -4.72 from bao_compatibility_2026_03_16.md)

**Pass criteria for f(delta) effect**:
- Peak positions shift < 0.5 sigma from Planck best-fit
- Chi^2 relative to Planck within delta chi^2 < 5 of LCDM
- BAO compatibility maintained (already shown to be better than LCDM)

**Fail criteria**:
- Peak positions shift > 2 sigma
- Overall chi^2 degradation > 10 relative to LCDM
- Inconsistency between CMB and BAO distance measures

### 4.3 The Secondary Diagnostic: Matter Power Spectrum

The Khronon also predicts a modified matter power spectrum P(k):
- At large k (small scales): identical to CDM (c_s^2 = 0)
- At small k (large scales): possible modification from f(delta) and late-time growth

The matter power spectrum comparison with SDSS/BOSS data would be a secondary test. The Fagin SLACS measurement already shows a 2.4-sigma preference for Khronon over CDM at lensing scales (fagin_comparison_result.md).

### 4.4 Summary: What Each Level Proves

| Level | If PASS | If FAIL |
|-------|---------|---------|
| 0 | Code works; GDM = CDM confirmed | Code bug (fix and retry) |
| 1 | TKS bound reproduced; c_s^2 threshold mapped | Unexpected; would indicate code issue |
| 2 | f(delta) correction compatible with CMB; possible better BAO fit | **Theory needs modification**: f(delta) too large |
| 3 | Full Khronon GDM mapping works; CMB compatible | **Serious problem**: time-dependent w(z) or k-dependent c_s^2 excluded |
| 4 | Full Khronon perturbation theory matches Planck | **Would be a major result**: first full CMB computation for Khronon framework |

---

## 5. Step-by-Step Implementation Plan

### Phase 1: Setup and Sanity Check (Day 1)

```bash
# Step 1: Clone the repository
git clone https://github.com/s-ilic/gdm_class_public
cd gdm_class_public

# Step 2: Compile
# On macOS, may need: export CC=gcc (or clang)
make clean
make class

# Step 3: Run standard LCDM
cp explanatory.ini test_lcdm.ini
# Edit test_lcdm.ini: set output = tCl, pCl, lCl
./class test_lcdm.ini

# Step 4: Run GDM = CDM
cp GDM.ini test_gdm_cdm.ini
# Edit: omega_gdm = 0.12, all w/cs2/cv2 = 0
./class test_gdm_cdm.ini

# Step 5: Compare outputs
# Use Python to load and compare C_l files
python3 -c "
import numpy as np
lcdm = np.loadtxt('output/lcdm_cl.dat')
gdm = np.loadtxt('output/gdm_cl.dat')
diff = np.max(np.abs((gdm[:,1] - lcdm[:,1]) / lcdm[:,1]))
print(f'Max relative difference: {diff:.2e}')
assert diff < 1e-6, 'GDM != CDM: bug detected'
print('PASS: GDM(w=0, cs2=0) = CDM')
"
```

### Phase 2: Constant c_s^2 Bracketing (Day 1-2)

```bash
# Run 8 configurations with different constant cs2
for cs2 in 0 1e-8 1e-7 1e-6 5e-6 1e-5 1e-4 5e-4; do
    # Create config file with this cs2
    # Run CLASS
    # Save output
done

# Compare all with Planck 2018 binned TT spectrum
# Plot residuals
# Determine exclusion threshold
```

**Python analysis script** (skeleton):
```python
import numpy as np
import matplotlib.pyplot as plt

# Load Planck 2018 binned power spectrum
# Available from: https://pla.esac.esa.int/
# Or use the Planck bestfit file included with CLASS

planck = np.loadtxt('planck_2018_binned_TT.dat')
ell_planck = planck[:, 0]
cl_planck = planck[:, 1]
sigma_planck = planck[:, 2]

cs2_values = [0, 1e-8, 1e-7, 1e-6, 5e-6, 1e-5, 1e-4, 5e-4]
for cs2 in cs2_values:
    cl_model = np.loadtxt(f'output/gdm_cs2_{cs2}_cl.dat')
    # Interpolate to Planck ell values
    # Compute chi^2
    chi2 = np.sum(((cl_model_interp - cl_planck) / sigma_planck)**2)
    print(f'cs2 = {cs2:.0e}: chi^2 = {chi2:.1f}')
```

### Phase 3: f(delta) Background Correction (Day 3-7)

**Option A: Binned w(a) approximation** (simpler, 1-2 days)

Pre-compute w(a) = delta(a)/2 at discrete scale factor values and use the binned GDM feature:

```python
import numpy as np

# Planck 2018 parameters
H0 = 67.36  # km/s/Mpc
Omega_m = 0.3153
Omega_r = 9.14e-5
Omega_L = 0.6847
delta_0 = 0.340

# Compute w(a) at bin edges
a_bins = np.array([1e-6, 1e-5, 1e-4, 1e-3, 0.003, 0.01, 0.03, 0.1, 0.3, 0.5, 0.8, 1.0])
for a in a_bins:
    z = 1/a - 1
    Ha_over_H0 = np.sqrt(Omega_r/a**4 + Omega_m/a**3 + Omega_L)
    delta_a = delta_0 * (1/Ha_over_H0)**2
    w_a = delta_a / 2
    print(f'a = {a:.0e}, z = {z:.0f}, w = {w_a:.3e}')
```

Expected output:
```
a = 1e-06, z = 999999, w ~ 5e-13
a = 1e-04, z = 9999, w ~ 5e-10
a = 9.1e-04, z = 1100, w ~ 3e-10
a = 0.01, z = 99, w ~ 3e-7
a = 0.1, z = 9, w ~ 5e-4
a = 0.5, z = 1, w ~ 0.04
a = 1.0, z = 0, w ~ 0.17
```

The w(a) is negligibly small at all CMB-relevant epochs (z > 100) and only becomes significant at z < 10.

**Option B: Modified background.c** (more rigorous, 3-5 days)

Directly modify the background integration to use the Khronon energy density formula. This requires understanding the CLASS background module structure, which iterates:
1. Guess H(a) from previous step
2. Compute all energy densities (including Khronon with f(delta))
3. Update H(a) from Friedmann equation
4. Check convergence

### Phase 4: Full k-Dependent c_s^2 (Day 7-12, if needed)

**Only needed if Level 1 shows sensitivity to c_s^2 at the predicted levels.**

From the cs2_verification_supplement.md analysis, c_s^2 ~ 3e-10 at z=1100 is 11,000x below the TKS bound. This is almost certainly undetectable.

If nevertheless needed, implement the 1/k^2 sound speed in perturbations.c:

```c
// In cs2_gdm_of_a_and_k(), add after existing k2_cs2_gdm block:
if (pba->kinv2_cs2_gdm == _TRUE_) {
    double k_pivot = 0.01;  // 1/Mpc
    double ca2_gdm = ppw->pvecback[pba->index_bg_ca2_gdm];
    cs2 = ca2_gdm + cs2 * k_pivot * k_pivot / (k * k);
    if (cs2 > 1.) cs2 = 1.;
    if (cs2 < 0.) cs2 = 0.;
}
```

Files to modify (4 changes total):
1. `/include/background.h`: Add `short kinv2_cs2_gdm;`
2. `/source/input.c` (parser): Add `class_read_flag("kinv2_cs2_gdm", pba->kinv2_cs2_gdm);`
3. `/source/input.c` (default): Add `pba->kinv2_cs2_gdm = _FALSE_;`
4. `/source/perturbations.c`: Add the code block above

### Phase 5: Comparison with Planck (Ongoing throughout)

For each level, produce:

1. **C_l^TT plot**: Model vs Planck 2018 best-fit, with Planck error bars
2. **Residual plot**: (C_l^model - C_l^LCDM) / sigma_Planck vs multipole l
3. **Chi^2 table**: Comparing model chi^2 with LCDM chi^2
4. **Parameter shifts**: How much do derived parameters (H_0, sigma_8, etc.) shift?

### Phase 6: Full Khronon Equations (Week 2-3, optional)

This is the expert-level implementation. Use `class.VFDM` (arXiv:2408.12052) as a template for adding a custom field species to CLASS. The key equations to implement:

**Background**:
```
rho_K = (mu^2 c^2 / (8 pi G)) * delta^2 * (1 + delta/2) / (1 + delta)^2
p_K = (mu^2 c^2 / (8 pi G)) * delta^3 / (2 * (1 + delta)^2)
```

**Perturbation** (synchronous gauge):
```
delta_K' = -(1 + w_K)(theta_K + h'/2) + 3 aH (w_K - c_ad^2) delta_K
theta_K' = -(1 - 3 c_ad^2) aH theta_K + c_s^2 k^2 delta_K / (1 + w_K)
```
with c_s^2 = 0 (at leading order) and c_ad^2 = w_K.

---

## 6. Risk Assessment: What Could Go Wrong?

### 6.1 Technical Risks

| Risk | Probability | Impact | Mitigation |
|------|------------|--------|------------|
| Compilation failure on macOS | MEDIUM | LOW | Use Docker with Linux image; or conda environment |
| GDM(w=0) != CDM (code bug) | LOW | HIGH | Start with Level 0 sanity check |
| Numerical instability at c_s^2 = 0 | LOW | MEDIUM | CLASS handles this; set cs2 = 1e-15 instead of 0 |
| f(delta) self-consistency loop fails | MEDIUM | MEDIUM | Use iterative approach; start from LCDM solution |
| Memory issues for high l_max | LOW | LOW | l_max = 2500 is standard; no issue on modern hardware |

### 6.2 Physics Risks

| Risk | Probability | Impact | Mitigation |
|------|------------|--------|------------|
| f(delta) correction too large at low z | MEDIUM | HIGH | If chi^2 degrades by > 10, theory needs revision |
| Running mu not self-consistent with H(z) | LOW | HIGH | Verify: delta(a) computed from LCDM H(a) vs self-consistent H(a) |
| ISW enhancement at low ell | MEDIUM | MEDIUM | Compare with cosmic variance (large errorbars at low ell) |
| Incorrect GDM mapping of Khronon | MEDIUM | HIGH | The GDM mapping is from Blanchet-Skordis 2024; should be reliable |
| Need c_vis^2 != 0 | LOW | MEDIUM | BS2024 predicts c_vis^2 = 0; if needed, add as free parameter |

### 6.3 The Three Most Likely Failure Modes

**Failure 1: f(delta) overcorrects H(z) at z < 2**

The 25.5% excess Khronon energy at z=0 requires 5.2% less Lambda. If the transition is too sharp, it could create a "dip" in H(z) around z ~ 0.5-1 that conflicts with supernovae or BAO. The bao_compatibility analysis (already done) suggests this is actually a FEATURE (better fit to DESI), but the full CMB computation would confirm.

**Failure 2: Omega_gdm h^2 needs re-tuning**

If the f(delta) correction changes the effective Omega_m h^2 at recombination, the peak positions would shift. The running mu ensures w(z=1100) ~ 3e-10, which should be safe, but the integrated effect on the sound horizon needs numerical verification.

**Failure 3: Initial conditions for Khronon perturbations**

At Level 4, the Khronon perturbation initial conditions (how they are set in the early radiation era) matter. CDM starts with adiabatic initial conditions (delta_c = (3/4) delta_gamma). The Khronon might have different initial conditions if its early-time dynamics differ from CDM. Getting this wrong would shift all peak heights.

### 6.4 What "Failure" Actually Means

**A failure at Level 2 (f(delta) excluded) would NOT kill the theory.** It would mean:
- The initial condition delta_0 needs to be re-fit (perhaps delta_0 < 0.340)
- Or the running prescription mu_bg(a) = H(a)/c needs modification
- Or the Khronon background energy formula needs correction

**A failure at Level 4 (full Khronon perturbations excluded) WOULD be serious.** It would mean:
- The omega = 0 dispersion is not sufficient for CMB compatibility
- The Khronon perturbation equations differ from CDM in a way that Planck detects
- The theory would need fundamental revision at the perturbation level

---

## 7. Existing Reference Implementations

### 7.1 Closest Analogues in the Literature

| Theory | Code | Public? | Reference | Relevance |
|--------|------|---------|-----------|-----------|
| GDM | gdm_class_public | YES | Ilic+ 2020, arXiv:2004.09572 | Direct starting point |
| Ultralight vector DM | class.VFDM | YES | Gomez+ 2024, arXiv:2408.12052 | Template for custom species |
| Dark energy radiation | class_der | YES | Berghaus+ 2023, arXiv:2311.08638 | Template for modified background |
| Scalar field DM | class_sfdm | YES | Urena-Lopez 2024 | Similar conceptual structure |
| AeST (Skordis-Zlosnik) | Custom (private) | NO | Skordis+ 2021, arXiv:2007.00082 | Exact theory; code unavailable |
| Khronon (Blanchet-Skordis) | Not implemented | NO | BS2024, arXiv:2404.06584 | The target theory |
| Khronon-Tensor | Not implemented | NO | BS2025, arXiv:2507.00912 | Extended version |

### 7.2 AeST Code Availability

Skordis and Zlosnik's AeST CMB code (which produced the PRL 2021 CMB fit) is NOT publicly available. However, Skordis is one of the authors of gdm_class_public, and the GDM framework was specifically designed to test AeST-like theories at the effective fluid level. This strongly suggests that gdm_class_public captures the essential physics.

### 7.3 What Blanchet-Skordis 2024 Actually Computed

From arXiv:2404.06584, Blanchet and Skordis:
1. Derived the Khronon perturbation equations on FRW
2. Showed the Khronon maps to GDM with specific {w(z), c_s^2(z,k)}
3. Argued (but did not numerically demonstrate) CMB compatibility
4. Their Fig. 3 shows qualitative C_l curves but NOT a quantitative Planck comparison

**Our computation would be the FIRST quantitative C_l comparison of the Khronon theory with Planck data.** This alone would be a publishable result.

---

## 8. Recommended Execution Strategy

### 8.1 Minimum Viable Product (MVP): Half a Day

**Do Level 0 + Level 1 only.**

This requires NO code modification -- just running gdm_class_public with different parameter values in the .ini file.

The output immediately answers:
- Is GDM(w=0, c_s^2=0) = CDM? (Should be yes)
- At what constant c_s^2 does the CMB become incompatible?
- Is c_s^2 = 5e-4 (naive prediction at first peak) excluded?

**If c_s^2 = 5e-4 is excluded** (very likely): Confirms that the omega = 0 dispersion (c_s^2 = 0) is the correct interpretation, not c_s^2 = (mu_0/k)^2.

**If c_s^2 = 5e-4 is NOT excluded** (very unlikely): Would suggest the k-dependent c_s^2 might be viable, motivating Level 2-3.

### 8.2 Full Program: 2-3 Weeks

```
Week 1: Levels 0-2
  Day 1: Setup, compile, Level 0 sanity check, Level 1 bracketing
  Day 2-3: Level 2 (binned w(a) approximation)
  Day 4-5: Level 2 (modified background.c, if binned approach insufficient)

Week 2: Level 3 + Analysis
  Day 6-7: Level 3 (k-dependent c_s^2, if needed)
  Day 8-9: Full Planck comparison (chi^2, residuals, parameter shifts)
  Day 10: Write up results

Week 3: Level 4 (Optional, for publication)
  Day 11-15: Full Khronon perturbation implementation
  Day 16-17: Validation and Planck comparison
  Day 18: Draft paper section
```

### 8.3 What to Publish

A paper (or paper section) containing:
1. First quantitative C_l^TT, C_l^TE, C_l^EE comparison of Khronon theory with Planck 2018
2. The f(delta) background effect on late-time observables (ISW, BAO)
3. Constant c_s^2 exclusion threshold (reproducing/refining TKS 2016)
4. Demonstration that running mu_bg(z) = H(z)/c ensures CMB compatibility
5. Prediction for CMB-S4 distinguishability

This would be the first numerically verified CMB computation for the Khronon dark matter framework.

---

## 9. Summary

| Question | Answer |
|----------|--------|
| Which solver? | CLASS via gdm_class_public |
| Can run on laptop? | YES (M1 MacBook is more than sufficient) |
| Time per C_l computation? | 5-30 seconds |
| Minimum viable test? | Half a day (Level 0+1, no code modification) |
| Full program? | 2-3 weeks |
| Key modification needed? | f(delta) background correction (Level 2) |
| Expected outcome? | CMB compatible at perturbation level (c_s^2 = 0 = CDM) |
| Distinguishing signal? | f(delta) background effect at z < 2 (ISW + BAO) |
| Publication potential? | YES (first quantitative Khronon CMB computation) |
| Biggest risk? | f(delta) overcorrecting H(z) at low redshift |

### Critical Path

```
Compile gdm_class → Level 0 sanity check → Level 1 c_s^2 bracketing
                                                     ↓
                              Confirm c_s^2 = 0 is the correct regime
                                                     ↓
                              Level 2: f(delta) background correction
                                                     ↓
                              Full Planck comparison + BAO cross-check
                                                     ↓
                              Level 4 (optional): Full Khronon equations
```

The minimum viable product (Level 0+1) can be completed in half a day with no code modification. The most scientifically interesting result (the f(delta) background effect, Level 2) requires 3-5 days of C coding. The full Khronon implementation (Level 4) is a 1-3 week project but would yield a publishable result.

---

## References

### Boltzmann Solvers
- CLASS: Lesgourgues 2011, arXiv:1104.2932; Blas, Lesgourgues & Tram 2011, arXiv:1104.2933
- gdm_class_public: Ilic, Kopp, Thomas & Skordis 2020, arXiv:2004.09572
- CAMB: Lewis, Challinor & Lasenby 2000, arXiv:astro-ph/9911177
- SymBoltz.jl: Sletmoen 2025, arXiv:2509.24740
- class.VFDM: Gomez et al. 2024, arXiv:2408.12052
- class_der: Berghaus et al. 2023, arXiv:2311.08638

### GDM Constraints
- Thomas, Kopp & Skordis 2016, arXiv:1601.05097 (c_s^2 < 3.4e-6)
- Kopp, Skordis, Thomas & Ilic 2016, arXiv:1605.00649 (extensive GDM investigation)
- Ilic et al. 2021, arXiv:1802.09541 (DM EoS through cosmic history)

### Khronon Theory
- Blanchet & Skordis 2024, arXiv:2404.06584, JCAP 11, 040 (Khronon theory)
- Blanchet & Skordis 2025, arXiv:2507.00912 (Khronon-Tensor extension)
- Skordis & Zlosnik 2021, arXiv:2007.00082, PRL 127, 161302 (AeST CMB fit)

### Planck Data
- Planck 2018: arXiv:1807.06209 (parameters)
- Planck Legacy Archive: https://pla.esac.esa.int/

### Internal References
- cmb_boltzmann_verification_plan.md (2026-03-14)
- cmb_verification_result.md (2026-03-13)
- cs2_verification_supplement.md (2026-03-14)
- bao_compatibility_2026_03_16.md (2026-03-16)
- paper4_khronon_CMB_calculation.md (2026-03-12)

---

*Last updated: 2026-03-16*
*This report assesses the practical feasibility of a full CMB Boltzmann solver computation for the Khronon dark matter theory.*
*Verdict: FEASIBLE on a laptop within 2-3 weeks, with minimum viable result in half a day.*
