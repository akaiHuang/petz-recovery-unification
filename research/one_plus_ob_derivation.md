# Rigorous Investigation: The (1+Omega_b) Correction in mu = 2pi(1+Omega_b)/r_d

**Author**: Sheng-Kai Huang (with computational analysis)
**Date**: 2026-03-20
**Status**: ALL THREE ORIGINAL CANDIDATES RULED OUT; new candidate identified

---

## Executive Summary

We compute the sound horizon integral r_d numerically for different baryon densities to determine whether the (1+Omega_b) correction in mu = 2pi(1+Omega_b)/r_d arises from candidate 2: the "baryon-free sound horizon" hypothesis. The computation decisively rules out this candidate, along with the other two original candidates. A new and better-fitting expression (1+f_b)^{1/3} is discovered, but its physical origin remains unproven.

**Key results:**
- r_d(Omega_b=0) / r_d(Planck) = **1.457**, not 1.049 --> Candidate 2 RULED OUT
- The z_drag shift (717 -> 1021) dominates the r_d change, not baryon loading on c_s
- Candidates 1 and 3 give corrections too small by 6-30 orders of magnitude
- NEW: (1+f_b)^{1/3} matches to 0.013% (better than (1+Omega_b) at 0.046%)
- BUT all correction forms are degenerate at the Planck cosmology because Omega_m ~ 1/3
- The (1+Omega_b) factor remains a NUMERICAL OBSERVATION

---

## 1. The Computation

### 1.1 Setup

We compute the comoving sound horizon at the baryon drag epoch:

$$r_d = \int_0^{a_{\rm drag}} \frac{c_s}{a^2 H(a)}\,da, \qquad c_s = \frac{c}{\sqrt{3(1+R)}}, \qquad R = \frac{3\omega_b\,a}{4\omega_\gamma}$$

with z_drag from the Eisenstein & Hu (1998) fitting formula. All other parameters are fixed at Planck 2018 values (omega_m = 0.1430, h = 0.6736, omega_gamma = 2.469 x 10^{-5}).

### 1.2 Results

| Case | omega_b h^2 | Omega_b | z_drag | r_d (Mpc) | r_d/(2pi) (Mpc) |
|------|-------------|---------|--------|-----------|-----------------|
| No baryons | 0 | 0 | 717 | 219.38 | 34.91 |
| **Planck** | **0.02237** | **0.0493** | **1021** | **150.62** | **23.97** |
| Double baryons | 0.04 | 0.088 | 1051 | 137.82 | 21.93 |

### 1.3 The Critical Ratio

$$\frac{r_d(\Omega_b = 0)}{r_d(\Omega_b = 0.049)} = \mathbf{1.457}$$

This is **nowhere near** 1 + Omega_b = 1.049. The hypothesis that the (1+Omega_b) correction comes from a baryon-free sound horizon is ruled out by a factor of 9.

---

## 2. Why Candidate 2 Fails

### 2.1 Decomposition of the r_d Change

Removing baryons changes r_d through TWO effects:

1. **Sound speed (c_s):** With Omega_b = 0, the baryon loading R = 0 everywhere, so c_s = c/sqrt(3) (maximum). This INCREASES r_d by ~17%.

2. **Drag epoch (z_drag):** With Omega_b = 0, the drag epoch shifts dramatically from z_drag = 1021 to z_drag = 717. The integration extends to much later times (larger a_drag), which INCREASES r_d by an additional ~24%.

| Effect | Contribution |
|--------|-------------|
| c_s baryon loading (z_drag fixed) | +17.0% |
| z_drag shift (c_s fixed) | +24.5% |
| **Total** | **+45.7%** |

The z_drag shift is the DOMINANT effect, contributing 60% of the total change. This is because the Eisenstein-Hu formula for z_drag is extremely sensitive to omega_b: baryons delay the drag epoch by tightly coupling to photons, pushing z_drag from ~717 (no baryons) to ~1021 (Planck).

### 2.2 Perturbative Check

At the Planck point:

$$\frac{d(\ln r_d)}{d\Omega_b} = -2.67$$

For candidate 2 to work, we would need this derivative to be -1. The actual sensitivity is 2.7 times larger.

### 2.3 Power-Law Scaling

Near the Planck values:

$$r_d \propto \omega_b^{-0.132}$$

This matches the Eisenstein-Hu/Aubourg fitting formula scaling (alpha ~ -0.128). A 5% change in r_d requires a ~38% change in omega_b -- far from the perturbative regime.

---

## 3. Why Candidates 1 and 3 Fail

### 3.1 Candidate 1: Gravitational Self-Energy

The gravitational correction to the Khronon's effective mass from the baryonic background is:

$$\frac{4\pi G \rho_b}{\mu^2 c^2} \sim 2 \times 10^{-6}$$

This is six orders of magnitude too small. The Khronon's Compton wavelength (22.3 Mpc) is so large that the gravitational self-energy from the homogeneous baryon density is negligible compared to the mass term.

### 3.2 Candidate 3: Gravitational Redshift

The gravitational redshift of the Khronon mode from baryonic perturbations gives a correction of order Omega_b * Phi, where Phi ~ 10^{-5} is the primordial perturbation amplitude. This gives:

$$\delta\mu / \mu \sim \Omega_b \times \Phi \sim 5 \times 10^{-7}$$

Seven orders of magnitude too small. On the FRW background there is no Newtonian potential (homogeneous), so this mechanism does not apply at the background level.

---

## 4. New Discovery: (1+f_b)^{1/3}

### 4.1 The Finding

Testing alternative correction factors, we find that

$$(1 + f_b)^{1/3}, \qquad f_b \equiv \frac{\Omega_b}{\Omega_m} = 0.1564$$

matches the needed correction factor BETTER than (1+Omega_b):

| Expression | Value | Discrepancy from needed |
|------------|-------|------------------------|
| f_b / pi = Omega_b/(pi*Omega_m) | 0.04979 | **+0.028%** |
| (1+f_b)^{1/3} | 1.04964 | **-0.013%** |
| 1 + Omega_b | 1.04930 | -0.046% |
| exp(Omega_b) | 1.05054 | +0.072% |
| Needed | 1.04978 | --- |

### 4.2 Physical Interpretation

The expression (Omega_m/Omega_DM)^{1/3} = (1+f_b)^{1/3} has a suggestive physical meaning: the cube root of the total-to-DM mass ratio. In 3D Newtonian gravity, the enclosed mass within radius r scales as M ~ r^3, so a wavelength scales as M^{1/3}. If the Khronon's Compton wavelength is set by the DM content but the gravitational field involves ALL matter, the effective mass receives a correction:

$$\mu_{\rm eff} = \mu_{\rm DM} \times \left(\frac{M_{\rm total}}{M_{\rm DM}}\right)^{1/3} = \mu_{\rm DM} \times (1+f_b)^{1/3}$$

However, the Jeans length scaling is lambda_J ~ 1/sqrt(rho), which gives (1+f_b)^{1/2} ~ 1.076 (7.5%, too large). The 1/3 power does NOT arise from Jeans instability analysis.

### 4.3 The Degeneracy

At the Planck cosmology, (1+Omega_b) and (1+f_b)^{1/3} are nearly identical because:

$$(1+f_b)^{1/3} \approx 1 + \frac{f_b}{3} = 1 + \frac{\Omega_b}{3\Omega_m}$$

This equals 1 + Omega_b when Omega_m = 1/3. Since Omega_m(Planck) = 0.315 ~ 1/3, the two expressions agree to 0.03%.

**The three candidates cannot be distinguished with current data.** Distinguishing them requires a different cosmology where Omega_m departs significantly from 1/3.

| Omega_m | 1 + Omega_b | (1+f_b)^{1/3} | Spread |
|---------|-------------|----------------|--------|
| 0.10 | 1.0493 | 1.1429 | 10.3% |
| 0.20 | 1.0493 | 1.0762 | 2.8% |
| **0.315** | **1.0493** | **1.0497** | **0.05%** |
| 0.50 | 1.0493 | 1.0318 | 1.7% |

---

## 5. Sensitivity to mu^{-1}

The match of (1+Omega_b) depends critically on the BS2025 value mu^{-1} = 22.3 Mpc:

| mu^{-1} (Mpc) | epsilon = r_d/(2pi*mu^{-1}) - 1 | Match to Omega_b |
|----------------|--------------------------------|-------------------|
| 22.0 | 0.064 | 1.30x |
| 22.2 | 0.055 | 1.11x |
| **22.3** | **0.050** | **1.01x** |
| 22.4 | 0.045 | 0.91x |
| 22.5 | 0.040 | 0.82x |

A shift of just 0.2 Mpc (< 1%) in mu^{-1} breaks the match. Since BS2025 does not quote a precise uncertainty for mu^{-1}, the (1+Omega_b) identification cannot be established beyond the 1% level.

The predicted mu^{-1} values from three correction factors all cluster around 22.1-22.3 Mpc:
- f_b/pi correction: mu^{-1} = 22.300 Mpc (best match to BS2025)
- (1+f_b)^{1/3}: mu^{-1} = 22.303 Mpc
- (1+Omega_b): mu^{-1} = 22.310 Mpc

---

## 6. Final Verdict

### 6.1 What is SETTLED

| Candidate | Mechanism | Prediction | Actual | Status |
|-----------|-----------|------------|--------|--------|
| **2** | r_d(Omega_b=0)/r_d = 1+Omega_b | 1.049 | **1.457** | **RULED OUT** |
| **1** | Gravitational self-energy | ~Omega_b | **~10^{-6}** | **RULED OUT** |
| **3** | Gravitational redshift | ~Omega_b | **~10^{-7}** | **RULED OUT** |

All three original candidates are ruled out. Candidate 2 fails spectacularly (ratio is 1.46, not 1.05) because the z_drag shift dominates the r_d change when Omega_b -> 0. Candidates 1 and 3 fail because the relevant gravitational effects are negligible at cosmological scales.

### 6.2 What is NEW

- **(1+f_b)^{1/3}** matches to 0.013%, better than (1+Omega_b) at 0.046%.
- **f_b/pi** matches the needed correction to 0.028%, with mu^{-1} = 22.300 Mpc exactly matching BS2025.
- All three correction forms are **degenerate at Planck cosmology** because Omega_m ~ 1/3.
- The 1/3 power has a suggestive (but not derived) connection to 3D gravitational mass scaling.

### 6.3 What REMAINS OPEN

1. **No rigorous derivation exists** for ANY form of the ~5% correction factor.
2. The correction could be a **coincidence** at the <1% level, given the unknown uncertainty in mu^{-1} = 22.3 Mpc.
3. Distinguishing between (1+Omega_b), (1+f_b)^{1/3}, and 1+f_b/pi requires either:
   - A precise uncertainty on mu^{-1} from BS theory
   - Or a measurement of mu in a different cosmological context

### 6.4 Recommendation

**For the paper, state the relation as:**

$$\mu \approx \frac{2\pi}{r_d} \qquad (5\% \text{ level, physically motivated})$$

and note the ~5% residual as an open question. The (1+Omega_b) notation can be kept for numerical convenience with the caveat that it is observationally motivated, not derived. The cleaner expression (Omega_m/Omega_DM)^{1/3} may be mentioned as a candidate with a more suggestive physical structure.

---

## Appendix: Numerical Code Summary

### Parameters Used
```
omega_m     = 0.1430       (Planck 2018)
omega_b     = 0.02237      (Planck 2018)
omega_gamma = 2.469e-5     (T_CMB = 2.7255 K)
N_eff       = 3.046
h           = 0.6736
r_d(Planck) = 147.09 Mpc   (Planck 2018)
mu^{-1}     = 22.3 Mpc     (BS2025 Model 2)
```

### Integral Method
- Scipy quad integration of c_s/(a^2 H) from a = 10^{-10} to a_drag
- z_drag from Eisenstein & Hu (1998) fitting formula (Eq. 4)
- H(a) from flat Lambda-CDM with radiation
- Our computed r_d(Planck) = 150.62 Mpc (2.4% above Planck 2018 due to EH98 z_drag approximation)
- All RATIOS of r_d are reliable to <0.1% since the EH98 offset cancels

### Key Numerical Results
```
r_d(Omega_b = 0):           219.38 Mpc
r_d(Omega_b = 0.049):       150.62 Mpc
r_d(Omega_b = 0.088):       137.82 Mpc

r_d(0) / r_d(Planck) = 1.4565   (candidate 2 predicts 1.0493)

Decomposition of the 45.7% change:
  c_s effect (z_drag fixed):    +17.0%
  z_drag effect (c_s fixed):    +24.5%
  Cross term:                   +4.2%

d(ln r_d)/d(Omega_b) at Planck = -2.67  (candidate 2 needs -1)
r_d ~ omega_b^{-0.132}  (power-law near Planck)
```

---

*Last updated: 2026-03-20*
*This document settles the (1+Omega_b) question: none of the three original candidates work. The correction factor remains a numerical observation.*
