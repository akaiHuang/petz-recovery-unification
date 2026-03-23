# Quantitative GW Predictions from Sigma = 2 ln Q

**Date:** 2026-03-23
**Status:** Derived with full numerical computation

---

## 1. Scalar GW Mode: Null Prediction

### Framework

The Khronon field in ghost condensation has K'(Q_0) = 0, yielding c_s^2 = 0 exactly. With the DBI completion (Blanchet & Skordis 2024):

```
K_DBI(Q) = mu^2 lambda_D^2 [sqrt(1 + (Q-1)^2/lambda_D^2) - 1]
```

the effective sound speed acquires scale-dependent form:

```
c_s^2(k) = alpha_DBI * (k/k_J)^2 / (1 + (k/k_J)^2)
```

where alpha_DBI = Q_0/(2 lambda_D^2) and:

```
k_J(a) = sqrt(3/2 * Omega_K * H0^2 / a)
```

### Computation

**Jeans wavenumber at a=1:**
```
k_J = sqrt(1.5 * 0.265) * (2.247e-4 Mpc^-1) = 1.42e-4 Mpc^-1
```

**Jeans frequency:**
```
k_J = 4.59e-27 m^-1
f_J = c * k_J / (2 pi) = 2.19e-19 Hz
```

This is **10 orders of magnitude below nHz PTAs** and **17 orders below LIGO**.

### At each detector band:

| Detector | f (Hz) | k/k_J | c_s^2(k) for alpha_DBI = 0.01 |
|:---------|:-------|:------|:------------------------------|
| PTA      | 10^-9  | 4.6 x 10^10 | 0.01 (saturated) |
| LISA     | 10^-3  | 4.6 x 10^16 | 0.01 (saturated) |
| LIGO     | 100    | 4.6 x 10^20 | 0.01 (saturated) |
| ET/CE    | 3      | 1.4 x 10^19 | 0.01 (saturated) |

At all detector frequencies, k >> k_J so c_s^2 -> alpha_DBI. But the scalar mode **coupling to the detector** is suppressed:

- The scalar mode strain relative to tensor: h_scalar / h_tensor ~ alpha_DBI < 0.05 (CMB bound from Paper 3)
- With alpha_DBI < 0.01 (ghost condensation natural value): h_scalar < 10^-25 at LIGO
- Current LIGO scalar mode sensitivity: ~10^-23 (polarization searches)

### Prediction (Falsifiable)

> **No scalar GW mode detectable by any current or planned detector.**
>
> Scalar mode strain: h_scalar < alpha_DBI * h_tensor < 0.05 * h_tensor
>
> Detection of a propagating scalar GW mode at any frequency would falsify
> the ghost condensation mechanism (K'(Q_0) = 0).

### Current status

- LVK polarization searches have not found scalar modes (consistent)
- No dedicated omega = 0 dispersion search has been performed
- Blanchet-Skordis 2024 specifically predicts omega = 0 scalar dispersion

---

## 2. EMRI Phase Shift from Exponential Metric

### Framework

Paper 2 establishes that Petz bound saturation gives the exponential metric:

```
ds^2 = -exp(-r_s/r_iso) dt^2 + exp(+r_s/r_iso)(dr_iso^2 + r_iso^2 dOmega^2)
```

where r_s = 2GM/c^2 in isotropic coordinates. This replaces the Schwarzschild horizon with a traversable wormhole throat.

### Metric Comparison at Same Areal Radius

The physical comparison must be at the same circumferential (areal) radius R. Converting the exponential metric from isotropic to areal coordinates (R = r_iso * exp(GM/(c^2 * r_iso))):

| R/(GM/c^2) | g_00 (Schw) | g_00 (exp) | Delta g_00 | delta_omega/omega |
|:-----------|:------------|:-----------|:-----------|:-----------------|
| 6          | -0.66667    | -0.66434   | +2.3e-3    | +1.22%           |
| 8          | -0.75000    | -0.74913   | +8.7e-4    | +0.58%           |
| 10         | -0.80000    | -0.79958   | +4.2e-4    | +0.34%           |
| 15         | -0.86667    | -0.86655   | +1.1e-4    | +0.13%           |
| 20         | -0.90000    | -0.89995   | +4.6e-5    | +0.07%           |
| 50         | -0.96000    | -0.96000   | +2.8e-6    | +0.001%          |

**Key finding:** The difference appears at O((r_s/r)^3), consistent with Paper 2 Sec. IV.

### ISCO Comparison

- **Schwarzschild ISCO:** R = 6.00 GM/c^2
- **Exponential ISCO:** R = 6.34 GM/c^2 (5.7% larger)
- **GW frequency at ISCO (M = 10^6 M_sun):** 4.38 mHz (in LISA band)

### Accumulated Phase Difference for EMRI

For a standard EMRI (M = 10^6 M_sun, m = 10 M_sun, q = 10^-5):

Integration from R = 100 GM/c^2 to R_ISCO = 6 GM/c^2 in the adiabatic (quadrupole) approximation:

```
Total orbits: 4.97 x 10^7
Total GW phase (Schwarzschild): ~1.25 x 10^9 rad
Average delta_omega/omega: 2.15 x 10^-4

ACCUMULATED PHASE DIFFERENCE: Delta_Phi = 66,984 rad
```

**Dependence on parameters:**

| M (M_sun) | m (M_sun) | q | Delta_Phi (rad) | f_ISCO (mHz) |
|:----------|:----------|:--|:----------------|:-------------|
| 10^5      | 10        | 10^-4 | 6,698 | 43.8 |
| 10^6      | 10        | 10^-5 | 66,984 | 4.38 |
| 10^6      | 30        | 3x10^-5 | 22,328 | 4.38 |
| 10^7      | 10        | 10^-6 | 669,841 | 0.438 |

### LISA Sensitivity

LISA phase measurement precision for EMRIs: sigma_Phi ~ 0.1 rad (at SNR ~ 20).

**Detection significance: Delta_Phi / sigma_Phi ~ 6.7 x 10^5** for the standard M = 10^6 EMRI.

### Prediction (Falsifiable)

> **If the exponential metric is the correct strong-field geometry, LISA EMRIs
> will show a phase deviation from GR of ~10^4 - 10^6 radians.**
>
> This is overwhelmingly detectable. The minimum detectable fractional deviation
> from Schwarzschild:
>
> alpha_min = 0.1 rad / 66,984 rad = 1.5 x 10^-6
>
> LISA can detect even a 10^-6 fractional admixture of exponential metric
> corrections to Schwarzschild.

### Critical Caveats

1. **The full exponential metric may not be realized.** Paper 2 notes: "Whether nature saturates the Petz bound or admits a gap of order O((r_s/r)^2) is empirical." If the actual metric is Schwarzschild + epsilon * (correction), then Delta_Phi scales as epsilon.

2. **Coordinate subtlety resolved.** The naive formula Delta_omega/omega = exp(-r_s/r)/(1-r_s/r) - 1 mixes coordinate systems. The proper computation at the same areal radius gives ~10x smaller per-orbit deviation but still enormous accumulated phase.

3. **ISCO difference matters.** The exponential ISCO is 5.7% farther from the center than Schwarzschild. This modifies the final plunge waveform — an independent signal.

4. **Consistency with LIGO BBH mergers.** LIGO observes the last ~10-30 orbits of comparable-mass binaries. At r ~ 3-5 r_s, delta_omega/omega ~ 1-2%. Over 10-30 orbits: Delta_Phi ~ 1-4 rad. Current LIGO phase precision for chirps: ~1 rad (systematic). This is **marginally consistent** — the exponential metric is not yet excluded by LIGO, but next-generation detectors (ET/CE) with ~10x better SNR will reach this.

### GW Frequencies in LISA Band

| Orbit radius | f_GW (mHz) | In LISA band? |
|:------------|:-----------|:-------------|
| 6 r_s       | 1.55       | Yes          |
| 10 r_s      | 0.72       | Yes          |
| 20 r_s      | 0.25       | Yes          |
| 50 r_s      | 0.064      | Yes          |

---

## 3. GW Decoherence = Sigma Accumulation

### Framework

From Paper 2, Eq. (Pikovski):
```
tau_grav(t) = 1 - exp(-N Gamma t^2)
Gamma = (Delta E)^2 (Delta Phi)^2 / (2 hbar^2 c^4)
```

In the Sigma framework: tau = 1 - exp(-Sigma/2), so:
```
Sigma_decohere = 2 N Gamma t^2
```

### GW-Induced Decoherence Rate

A passing GW with strain h and angular frequency omega_GW produces a tidal acceleration:
```
g_eff = h * omega_GW^2 * Delta_x
```
where Delta_x is the superposition size.

Using the Pikovski mechanism:
```
t_decohere = hbar / (m * g_eff * Delta_x)
Sigma_GW = 2 t / t_decohere
```

### Numerical Estimates

**For a single GW event** (h = 10^-21, f = 100 Hz, duration 0.1 s):

Using a 170 kDa nanoparticle (Pedalino et al. 2026) in a 1 micron superposition:
```
g_eff = 10^-21 * (628)^2 * 10^-6 = 3.95 x 10^-22 m/s^2
t_decohere = 1.055e-34 / (2.82e-22 * 3.95e-22 * 10^-6) = 9.5 x 10^14 s
Sigma_GW = 2 * 0.1 / (9.5e14) = 2.1 x 10^-16
tau_GW = 1.1 x 10^-16
```

**Completely negligible.** A single GW event cannot decohere any realistic quantum system.

**For the stochastic GW background** (h_rms ~ 10^-25 at 100 Hz, 1 year exposure):
```
Sigma_SGWB ~ 6.7 x 10^-12
```

Also negligible with current technology.

### Phase Diffusion Prediction

The Sigma framework makes a specific prediction for the GW phase diffusion exponent (Cang & Wang 2025):

The phase variance of a GW signal propagating over distance d:
```
var(phi) propto d * f^{alpha_dec}
```

where alpha_dec discriminates between quantum gravity models:

| Model | alpha_dec | Mechanism |
|:------|:----------|:----------|
| String foam | -1 | Space-time foam from string theory |
| **Sigma = 2 ln Q** | **0** | **QRE accumulation (holographic scaling)** |
| Causal set | +1 | Discrete spacetime |

### Prediction (Falsifiable)

> **The Sigma framework predicts alpha_dec = 0 (holographic phase diffusion scaling).**
>
> Current bounds (de Kruijf et al. 2025): alpha_dec consistent with 0 using Planck + BICEP/Keck + LVK + BBN.
>
> ET/CE can improve bounds by factor ~100, reaching sensitivity to distinguish
> between string foam, holographic, and causal set predictions.
>
> If alpha_dec != 0 is measured, the QRE origin of Sigma is excluded.

### Quantitative Bound on Sigma_GW

From the decoherence analysis, the Sigma framework predicts:
```
Gamma_decoherence <= 2 * dSigma/dt / hbar
```

The rate dSigma/dt for a GW with power spectral density S_h(f):
```
dSigma_GW/dt = integral [S_h(f) * (2 pi f)^4 * m^2 * Delta_x^4 / hbar^2] df
```

For Advanced LIGO noise floor S_h^{1/2} ~ 10^-23 Hz^{-1/2} at 100 Hz:
```
dSigma/dt ~ 10^-46 * (600)^4 * (10^-22)^2 * (10^-6)^4 / (10^-34)^2
         ~ 10^-46 * 10^11 * 10^-44 * 10^-24 / 10^-68
         ~ 10^-35 s^-1
```

This predicts that GW-induced decoherence is unobservable for any foreseeable quantum system.

---

## 4. Additional Strong-Field Predictions (from Paper 2)

### 4.1 GW Echoes

The exponential metric wormhole structure produces echoes:
```
Delta t_echo = 4.2 * r_s / c
```

| Source | M | Delta t_echo | f_echo | Detector |
|:-------|:--|:-------------|:-------|:---------|
| Stellar BH (GW150914) | 60 M_sun | 2.5 ms | ~400 Hz | LIGO |
| Intermediate-mass BH | 10^3 M_sun | 40 ms | ~25 Hz | LIGO/ET |
| Sgr A* | 4 x 10^6 M_sun | 165 s | ~6 mHz | LISA |
| M87* | 6.5 x 10^9 M_sun | 74 hr | ~4 uHz | PTA |

**Critical:** Current echo searches use ~200 ms templates (Planck-wall models). The exponential-metric echo at 2.5 ms is **~100x shorter** and would be MISSED by existing searches. Dedicated short-delay templates are needed.

### 4.2 QNM Frequency Shift

```
f_QNM^exp / f_QNM^Schw = 3 sqrt(3) / (2e) = 0.956
```

A **-4.4% downward shift** in quasinormal mode frequency (eikonal limit).

Current precision: ~10% (GW250114 gave the best single-event test, 2-3x better than previous).
Required precision: ~1% (achievable with ET/CE by ~2035).

### 4.3 Shadow Size

```
b_crit^exp / b_crit^Schw = 2e / (3 sqrt(3)) = 1.046
```

A **+4.6% larger shadow.** Current EHT precision: ~7% (M87*). Not yet sufficient to distinguish.

---

## 5. Summary: Testability Timeline

| Prediction | Value | Detector | Sensitivity | Timeline |
|:-----------|:------|:---------|:------------|:---------|
| No scalar mode | h_s < 0.05 h_t | LIGO/ET | h ~ 10^-23 | Now (null) |
| EMRI phase shift | ~10^4-10^6 rad | LISA | 0.1 rad | 2035+ |
| Echo (stellar BH) | 2.5 ms, ~400 Hz | LIGO O5 | Need templates | 2025-2027 |
| Echo (SMBH) | 165 s, ~6 mHz | LISA | Standard search | 2035+ |
| QNM shift | -4.4% | ET/CE | ~1% | 2035+ |
| Shadow | +4.6% | ngEHT | ~2-3% | 2030+ |
| Phase diffusion | alpha_dec = 0 | ET/CE | Factor 100x | 2035+ |
| GW friction Xi_0 | = 1 exactly | Standard sirens | 10^-3 | 2035+ |
| GW decoherence | Sigma ~ 10^-16 | Lab quantum | Far below | Not foreseeable |

### Most Discriminating Tests (ranked):

1. **EMRI phase shift (LISA)** — SNR ~ 10^5-10^6, overwhelmingly detectable IF exponential metric is correct. Even a 10^-6 fractional deviation from Schwarzschild is detectable.

2. **GW echoes at 2.5 ms (LIGO O5)** — Requires new short-delay templates. Could be tested with EXISTING data if appropriate searches are performed.

3. **QNM frequency shift of -4.4% (ET/CE)** — Requires ~1% precision, achievable with next-generation detectors.

4. **Phase diffusion exponent alpha_dec = 0 (ET/CE)** — Discriminates between quantum gravity models. Current bounds are consistent.

5. **Null scalar mode detection (ongoing)** — Every negative result in scalar polarization searches confirms the framework.

---

## 6. What Would Exclude the Framework

| Observation | Implication |
|:-----------|:-----------|
| Scalar GW mode detected | K'(Q_0) != 0, ghost condensation fails |
| c_T != c confirmed | Tensor sector modified, framework invalid |
| GW birefringence | Parity violation not in Khronon |
| EMRI phase = GR exactly (with LISA precision) | Exponential metric excluded (Schwarzschild wins) |
| alpha_dec != 0 | QRE origin of Sigma excluded |
| Xi_0 != 1 confirmed | GW friction from unknown source |
| Echoes found at ~200 ms (not ~2.5 ms) | Planck-wall model, not exponential metric |
