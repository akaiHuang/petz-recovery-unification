# Exponential vs Schwarzschild: Experimental Distinguishability
## Date: 2026-03-10

## Executive Summary

The exponential metric g_00 = -exp(-r_s/r) and Schwarzschild g_00 = -(1-r_s/r) are indistinguishable through second post-Newtonian order (identical PPN parameters gamma = beta = 1). Their first difference appears at O((r_s/r)^3) ~ 10^{-18} in the solar system. In strong fields (r_s/r > 0.3), five predictions diverge at the 4-20% level. **As of March 2026, the exponential metric is NOT definitively ruled out**, but GW250114 QNM data is beginning to probe the relevant regime. The closest constraint comes from the pSEOBNR ringdown analysis of GW250114, which constrains delta_f_220 to [0.00, 0.04] at 90% credibility -- this would exclude the eikonal -4.4% shift IF the eikonal estimate is exact, but the ringdown-only (pyRing) analysis gives [-0.13, +0.43], within which even -7.8% remains viable.

---

## 1. Current Status Summary

### What Paper 2 predicts (five testable signatures):

| # | Observable | Exponential prediction | Schwarzschild | Difference | Current precision |
|---|-----------|----------------------|---------------|------------|------------------|
| 1 | Shadow diameter (M87*) | 41.5 uas | 39.7 uas | +4.6% | +/-7% (EHT 2019) |
| 2 | QNM fundamental frequency | 0.956 f_Schw (eikonal) | f_Schw | -4.4% | ~tens of % (GW250114) |
| 3 | GW echo delay (60 M_sun) | ~2.5 ms | No echo (horizon) | qualitative | searches at ~200 ms |
| 4 | ISCO radius | 3.17 r_s | 3 r_s | +5.6% | ~10% (NICER) |
| 5 | Gravitational redshift (NS) | z = 0.181 | z = 0.225 | -19% | ~5% (NICER) |

### Bottom line:
- **Solar system**: Completely degenerate with Schwarzschild. No constraint possible.
- **EHT shadows**: Current ~7-17% uncertainty; 4.6% difference NOT yet resolvable.
- **LIGO QNMs**: GW250114 is approaching but NOT yet decisively constraining.
- **GW echoes**: Current searches target wrong delay range (100-200 ms vs predicted 1-10 ms).
- **NICER**: 19% redshift difference is in principle the largest discriminator, but current ~5% compactness precision does not yet resolve the metric difference due to EOS degeneracies.

---

## 2. EHT Shadow Constraints

### 2.1 M87* Measurements

**EHT 2019 (original)**:
- Shadow diameter: 42 +/- 3 uas
- Schwarzschild deviation parameter: delta = -0.01 +/- 0.17
- This means ~17% fractional uncertainty on the deviation from Schwarzschild
- The exponential metric predicts delta ~ +0.046 (4.6% larger shadow)
- **Status: Well within error bars. Cannot distinguish.**

**EHT 2024 (persistent shadow, 2018 data)**:
- Confirmed persistent shadow with angular diameter consistent with 2017 observations
- Achieved 345 GHz observations in August 2024 (highest diffraction-limited angular resolution on Earth)
- Shadow size measurement is consistent across epochs, reducing systematic concerns
- **Precision has NOT significantly improved beyond the original ~17% deviation uncertainty**

**Exponential vs Schwarzschild for M87***:
- Schwarzschild prediction: theta_sh = 39.7 uas
- Exponential prediction: theta_sh = 41.5 uas
- EHT measurement: 42 +/- 3 uas
- Both predictions lie within 1-sigma of the measurement
- The exponential value (41.5) is actually closer to the central measured value (42)

### 2.2 Sgr A* Measurements

**EHT 2022**:
- Shadow diameter: 48.7 +/- 7 uas
- Schwarzschild deviation:
  - VLTI: delta = -0.08 (+0.09, -0.09)
  - Keck: delta = -0.04 (+0.09, -0.10)
- Fractional uncertainty: ~9-10% in the deviation parameter
- **Status: Even larger error bars than M87*. Cannot distinguish.**

### 2.3 Next-Generation EHT (ngEHT)

- **Phase 1 (2023-2026)**: Add 5 new VLBI stations, enhance current array
- **Phase 2 (2026-2030)**: Add ~10 total new sites worldwide
- Expected angular resolution improvement: ~50% over current EHT
- **At 345 GHz**: resolution ~10 uas
- The ngEHT may achieve shadow size precision sufficient to probe 4.6% deviations, but this requires reducing both statistical AND systematic uncertainties (astrophysical modeling of accretion flow, jet contamination, etc.)
- **Lunar-based VLBI (proposed)**: Could achieve ~0.85 uas angular resolution at 230 GHz

### 2.4 Assessment for Paper 2

The EHT shadow test is **not currently constraining** and will likely not become so until ngEHT Phase 2 (~2030). The 4.6% shadow size difference requires <2% fractional precision, which is challenging given astrophysical systematics. The shadow size test alone is unlikely to be the first to discriminate the metrics.

**Key reference**: EHT Collaboration 2024, A&A 686, A154 (M87* persistent shadow)

---

## 3. LIGO/Virgo GW Echo Searches

### 3.1 Current Echo Search Results (O1-O3)

**LVK Echo Search (2023, GWTC-3 data)**:
- Analyzed 33 BBH merger events with network SNR >= 10
- Search window: <1 second post-merger
- Echo delay times searched: t_echo up to 200 ms (blind time 50 ms, post-merger window 300 ms)
- **Result: No evidence of echoes found**
- Upper limits on echo amplitude ratio A: 0.11 to 0.75 (event-dependent)
- Specific examples:
  - GW150914: h_rss_UL = 3.4 x 10^{-23}/sqrt(Hz), amplitude ratio 0.21
  - GW190814: h_rss_UL = 1.5 x 10^{-23}/sqrt(Hz), amplitude ratio 0.75
  - GW190828_063405: h_rss_UL = 1.6 x 10^{-23}/sqrt(Hz), amplitude ratio 0.11

### 3.2 Critical Mismatch with Exponential Metric Prediction

**The current LVK echo searches are optimized for the WRONG delay range for the exponential metric.**

| Model | Echo delay (60 M_sun) | Frequency | Currently searched? |
|-------|----------------------|-----------|-------------------|
| Planck-wall (Damour-Solodukhin) | ~220 ms | ~4.5 Hz | YES |
| Quantum firewall | ~100-200 ms | ~5-10 Hz | YES |
| **Exponential metric** | **~2.5 ms** | **~400 Hz** | **NO** |

Paper 2's echo prediction (Delta_t ~ 4.17 r_s/c):
- 20 M_sun remnant: 0.82 ms, 1218 Hz
- 60 M_sun remnant: 2.46 ms, 406 Hz
- 100 M_sun remnant: 4.10 ms, 244 Hz
- 150 M_sun remnant: 6.16 ms, 162 Hz

**All exponential-metric echoes fall in the LIGO sensitivity band (10-2000 Hz)** but at delays ~100x shorter than current search templates.

### 3.3 O4 Results (2023-2025)

- O4 ran from May 24, 2023 to November 18, 2025
- ~250 candidate signals detected in real time
- GWTC-4.0 released August 2025 with 128 significant events from O4a
- O4b included GW250114 (SNR ~77-80)
- **No dedicated short-delay echo search has been published from O4 data as of March 2026**
- A model-independent search targeting long-lived QNMs from strong interior reflection was applied to GW150914, GW231226, and GW250114 -- no statistically significant evidence for postmerger echoes found

### 3.4 Assessment

**The exponential metric echo prediction has NOT been tested.** This is because:
1. Current templates are optimized for ~100-200 ms delays
2. The ~2.5 ms delay falls in a regime that overlaps with the tail of the ringdown itself
3. Separating a 2.5 ms echo from the ringdown requires careful template modeling

**This is the most actionable gap**: A dedicated search for short-delay echoes (1-10 ms) in O3/O4 data could either detect or constrain the exponential metric. Paper 2 correctly calls for this reanalysis.

---

## 4. QNM Measurements

### 4.1 GW250114: The Gold Standard

GW250114 (January 14, 2025) is the loudest GW signal ever detected:
- Binary black hole merger, total mass ~120 M_sun (remnant ~62.7 M_sun)
- Signal-to-noise ratio: ~77-80
- First detection of the l=m=2 overtone at 4.1 sigma significance

**QNM frequency constraints from GW250114**:

| Analysis method | Constraint on delta_f_220 | Exponential prediction | Consistent? |
|----------------|--------------------------|----------------------|-------------|
| pSEOBNR (full IMR) | [0.00, 0.04] at 90% | -0.044 (eikonal) | **Marginal exclusion** |
| pyRing (ringdown only) | [-0.13, +0.43] | -0.044 to -0.078 | **YES, consistent** |
| TIGER (parametrized) | "tens of percent" | -0.044 | YES |

**Key subtlety**: The pSEOBNR analysis uses the full inspiral-merger-ringdown waveform and constrains delta_f_220 to [0.00, 0.04], which would marginally exclude the eikonal prediction of -4.4%. However:

1. The exponential metric's QNM frequency has NOT been computed exactly (only the eikonal estimate exists)
2. Nath-Sarma (2024) found -7.8% for a generalized exponential wormhole using WKB -- larger than eikonal
3. The pyRing (ringdown-only) analysis, which is more model-independent, gives [-0.13, +0.43], well encompassing both estimates
4. The pSEOBNR constraint assumes a Kerr template for the inspiral, which biases the result if the inspiral itself deviates

### 4.2 Previous Constraints (O1-O3)

- Combined O1-O3 constraints from GWTC-3: QNM frequencies consistent with Kerr to ~10-30% for individual events
- GW150914 alone: ~20% precision on dominant QNM frequency
- Stacking multiple events: ~few percent precision, but model-dependent

### 4.3 Future: Einstein Telescope and Cosmic Explorer

- **Einstein Telescope (ET)**: Expected to measure two independent QNM modes to O(1)% relative uncertainty
- **Cosmic Explorer (CE)**: Similar or better precision for loud events
- **Combined ET+CE**: ~1000 events/year with <10% QNM precision individually
- **Timeline**: ET first science ~2035; CE ~2035-2040

**At O(1)% QNM precision, a 4.4% shift would be clearly detectable or excluded.**

### 4.4 Assessment

The QNM test is currently the **closest to constraining** the exponential metric. The pSEOBNR analysis of GW250114 is marginally in tension with the eikonal estimate. However:
- The exact QNM spectrum of the pure exponential metric remains uncomputed
- Ringdown-only analyses remain consistent
- Definitive constraint requires either exact QNM computation OR ET/CE precision

**Priority for Paper 2**: Computing the exact l=2, n=0 QNM frequency of the pure exponential metric (not the generalized version of Nath-Sarma) is the single most impactful calculation to perform.

---

## 5. NICER Neutron Star Constraints

### 5.1 Current NICER Measurements

**PSR J0030+0451** (Riley et al. 2019, Miller et al. 2019, updated 2024):
- Mass: M = 1.44 +/- 0.15 M_sun
- Radius: R = 13 (+1.2, -1.0) km
- Compactness: r_s/R ~ 0.31

**PSR J0740+6620** (Miller et al. 2021, updated 2024):
- Mass: M = 2.08 +/- 0.07 M_sun
- Radius: R = 13.7 (+2.6, -1.5) km (68% CI)
- Compactness: r_s/R ~ 0.45

**PSR J0437-4715** (Choudhury et al. 2024):
- New radius measurement for this 1.42 M_sun pulsar
- Adds additional constraint on EOS

**Overall precision**: Radius known to ~5% for a 1.4 M_sun neutron star.

### 5.2 Gravitational Redshift Comparison

At neutron star compactness r_s/r ~ 0.33:
- Schwarzschild: z = 0.225
- Exponential: z = 0.181
- **Difference: 19% (the largest fractional difference among all observables)**

However, NICER does not directly measure z. It measures:
- X-ray pulse profiles -> mass and radius
- The redshift is inferred from M/R
- EOS degeneracies dominate the error budget

### 5.3 Can NICER Distinguish the Metrics?

**Not directly** with current data. The issue:
1. NICER constrains M and R through pulse profile modeling
2. The metric enters through light bending and time delay in the pulse profile
3. At r_s/R ~ 0.3, the exponential and Schwarzschild metrics differ by:
   - Light bending: ~O((r_s/r)^3) correction
   - Time delay: ~O((r_s/r)^3) correction
   - These are ~3% effects, comparable to current systematic uncertainties

4. The 19% redshift difference is a **derived quantity** that requires knowing the full metric, not just M and R independently

### 5.4 Future Missions

- **eXTP (enhanced X-ray Timing and Polarimetry)**: Chinese-European mission
  - Expected M-R precision: ~2% ("best case" scenario)
  - This would start to probe the metric difference at neutron star surfaces
  - Launch timeline: ~2029

- **STROBE-X**: Proposed NASA mission
  - Similar precision goals to eXTP
  - Would provide independent constraints

### 5.5 NICER Science Operations Suspended

- NICER science operations were suspended in June 2025
- All published data remain available
- Future NICER data releases may still occur from already-collected observations

### 5.6 Assessment

The neutron star gravitational redshift is in principle the **strongest discriminator** (19% difference). However:
- Current observations cannot yet directly test the metric at neutron star surfaces
- EOS degeneracies are the dominant systematic
- eXTP (~2029) is the most promising near-future mission for this test
- **Key insight**: A simultaneous measurement of M, R, AND z for the same neutron star would be the cleanest test

---

## 6. Other Observational Tests

### 6.1 Gravitational Wave Memory Effect

**Current status**:
- GW memory has never been detected
- LIGO/Virgo are mostly insensitive to the memory effect (low-frequency regime)
- **LISA** (launch ~2035): Expected to detect GW memory from massive BH binary mergers (10^5-10^7 M_sun)
  - Recent analysis (Inchauspe et al. 2024, PRD 111, 044044) shows LISA can measure GW memory for several events
  - Bayesian analysis (2026, arXiv:2601.23230) confirms detectability

**For exponential metric**:
- The memory effect depends on the asymptotic falloff of the metric
- Since exponential and Schwarzschild metrics share the same weak-field/asymptotic behavior, the **memory effect is identical** at leading order
- Differences would appear only at O((r_s/r)^3) corrections to the radiation
- **Conclusion: GW memory is NOT a useful discriminator for these two metrics**

### 6.2 Pulsar Timing: Double Pulsar PSR J0737-3039

**Current status (Kramer et al. 2021, Phys. Rev. X 11, 041050)**:
- 16 years of timing data
- 7 post-Keplerian parameters measured
- GR verified to 0.05% precision (via parameter s)
- First observation of relativistic light propagation effects in strong fields
- 2PN periastron advance contribution is being probed (~same order as Lense-Thirring precession)

**For exponential vs Schwarzschild**:
- Both metrics have identical PPN parameters (gamma = beta = 1)
- All timing effects through 1PN are identical
- **2PN periastron advance differs**: The coefficient of (r_s/r)^3 term differs (-1/6 vs -1/4)
- For the double pulsar: r_s/r ~ 10^{-5} per orbit
- 2PN correction: ~ (r_s/r)^3 ~ 10^{-15} timing effect
- **Current timing precision: ~microseconds = 10^{-6} of orbital period**
- **Need: ~10^{-15} timing precision to detect the metric difference**
- **Conclusion: Pulsar timing CANNOT distinguish the metrics with current or foreseeable technology**

### 6.3 Solar System Tests

**Cassini Shapiro delay** (Bertotti et al. 2003):
- Constrains |gamma - 1| < 2.3 x 10^{-5}
- Both metrics have gamma = 1 exactly
- **No constraint**

**Mercury perihelion precession**:
- Both metrics predict identical precession through O(r_s/r)
- The O((r_s/r)^3) correction to Mercury's precession is ~10^{-18} of the classical value
- **No constraint**

**Lunar Laser Ranging**:
- Tests Nordtvedt effect, PPN parameters
- Both metrics pass identically
- **No constraint**

### 6.4 Strong Gravitational Lensing

**Manna, Rahman, Chowdhury (2023)**:
- Computed strong lensing in exponential wormhole spacetime
- Found differences in deflection angle compared to Schwarzschild at high compactness
- Relevant for future micro-lensing surveys of stellar-mass compact objects
- **Current data: Cannot yet distinguish** (requires resolving individual Einstein rings of stellar-mass objects)

### 6.5 X-ray Reflection Spectroscopy

**Iron K-alpha line profile**:
- Depends on ISCO location (5.6% difference) and metric geometry near ISCO
- Current XMM-Newton/Chandra/NuSTAR data have ~10-20% systematic uncertainties in ISCO determination
- **eXTP/Athena** may reach ~5% precision on iron line profiles
- **Conclusion: Not yet constraining, but within reach of next-generation X-ray observatories**

---

## 7. Calculated Predictions (ISCO, precession, Shapiro delay)

### 7.1 ISCO Radius and Orbital Frequency

**Schwarzschild** (in Schwarzschild coordinates):
- r_ISCO = 6M = 3r_s
- Orbital frequency: f_ISCO = c^3 / (6^{3/2} pi r_s c) = 1/(6sqrt(6) pi) * c/r_s

**Exponential metric** (Paper 2 computation):
- R_ISCO ~ 3.17 r_s (areal radius)
- **Difference: +5.6% outward shift**

**Ernazarov 2025 (arXiv:2511.13471) independent calculation**:
- For the pure exponential metric with l=0: r_ISCO = 2M (in isotropic coordinates)
- This corresponds to the same result in areal coordinates

**Orbital frequency at ISCO**:
- Exponential: f_ISCO^exp ~ 0.97 f_ISCO^Schw (due to larger ISCO radius)
- Difference: ~3% in orbital frequency at ISCO
- For a 10 M_sun object: f_ISCO^Schw ~ 220 Hz, f_ISCO^exp ~ 213 Hz

**Radiative efficiency**:
- Schwarzschild: eta = 1 - sqrt(1 - 2/(3*2)) = 1 - sqrt(8/9) ~ 5.72%
- Exponential: eta ~ 5.48%
- Difference: 0.24 percentage points (4.2% relative)

### 7.2 Photon Sphere and Shadow

**Schwarzschild**:
- Photon sphere: r_ph = 3M (Schwarzschild coords) = 1.87m (isotropic coords)
- Critical impact parameter: b_crit = 3sqrt(3)M ~ 5.196M

**Exponential metric**:
- Photon sphere: r_ph = r_s = 2M (isotropic coords)
- The photon sphere coincides with the wormhole throat
- Critical impact parameter: b_crit = 2eM ~ 5.437M

**Shadow diameter difference**: (5.437 - 5.196)/5.196 = +4.64%

**Nath-Sarma (2024) results for generalized exponential wormhole**:
- Photon sphere radius is less than Schwarzschild (consistent)
- Shadow size is greater than Schwarzschild (consistent with 4.6% prediction)
- QNMs computed via 3rd-order WKB: fundamental l=2 mode shifted by -7.8% from Schwarzschild

### 7.3 Perihelion Precession

**General formula for precession per orbit**:

For a metric with g_00 = -f(r), the precession per orbit to leading order is:
Delta_phi = 6pi M / (a(1-e^2))

Both metrics give identical precession at this order. The difference appears at the NEXT order.

**Schwarzschild** (exact to 2PN):
Delta_phi = (6pi M) / (a(1-e^2)) * [1 + (3/4 - e^2/4) * (r_s/(a(1-e^2)))^2 + ...]

The coefficient at 2PN: 3/4 - e^2/4 (depends on eccentricity)

**Exponential metric** (to next order):
The difference enters at O((r_s/a)^3):
- Schwarzschild coefficient of (r_s/r)^3 in g_00: -1/4 (isotropic coords)
- Exponential coefficient of (r_s/r)^3 in g_00: -1/6

**For Mercury** (a = 5.79 x 10^10 m, e = 0.205, r_s = 2.95 km):
- r_s/a ~ 5.1 x 10^{-8}
- (r_s/a)^3 ~ 1.3 x 10^{-22}
- Precession difference: ~ (1/4 - 1/6) * (r_s/a)^2 * 42.98"/century
- = (1/12) * (2.6 x 10^{-15}) * 42.98"/century
- ~ 10^{-14} arcseconds/century
- **Completely unmeasurable** (current precision: ~0.002"/century)

### 7.4 Shapiro Time Delay

**Schwarzschild** (for light passing near the Sun):
Delta_t = (2r_s/c) * [1 + ln(4r_1 r_2 / d^2)]

where d is the closest approach distance, r_1, r_2 are distances to emitter/receiver.

**Exponential metric**:
At leading order (r_s/r), the Shapiro delay is identical. The correction appears at O((r_s/r)^2):

Schwarzschild: Delta_t includes a term ~ (r_s/d)^2 with coefficient 15pi/16
Exponential: coefficient differs by O(1) factor

For the Cassini experiment (d ~ R_sun = 7 x 10^8 m, r_s = 2.95 km):
- (r_s/d)^2 ~ 1.8 x 10^{-11}
- Correction to Shapiro delay: ~ (r_s/c) * (r_s/d)^2 ~ 10^{-22} seconds
- **Completely unmeasurable** (Cassini timing precision: ~10^{-5} seconds)

### 7.5 Light Deflection

**Both metrics predict identical deflection at 1PN**:
Delta_theta = 4M/b (1 + cos(phi))

The difference at 2PN involves the coefficient of (M/b)^2:
- Schwarzschild: 15pi/4 * (M/b)^2
- Exponential: same at this order (both agree through 2PN)

Difference appears at 3PN ~ (M/b)^3 ~ 10^{-18} for solar-system light bending.

---

## 8. Is Exponential Metric Ruled Out?

### 8.1 Summary of Constraints

| Test | Exponential prediction | Current data | Ruled out? |
|------|----------------------|-------------|------------|
| Solar system (PPN) | gamma = beta = 1 | gamma = 1 +/- 2.3e-5 | NO |
| Mercury precession | 42.98"/century + 10^{-14} | 42.98 +/- 0.002 | NO |
| Shapiro delay | Identical to 1PN | Cassini bound | NO |
| M87* shadow | 41.5 uas | 42 +/- 3 uas | NO |
| Sgr A* shadow | larger by 4.6% | 48.7 +/- 7 uas | NO |
| QNM (GW250114 pSEOBNR) | delta_f ~ -4.4% to -7.8% | [0.00, 0.04] | **MARGINAL** |
| QNM (GW250114 pyRing) | delta_f ~ -4.4% to -7.8% | [-0.13, +0.43] | NO |
| GW echoes | ~2.5 ms for 60 M_sun | Not searched at this delay | NOT TESTED |
| NS redshift | 19% smaller at r_s/r=0.33 | ~5% radius precision | NO (indirect) |
| Pulsar timing | O((r_s/r)^3) difference | 0.05% GR agreement | NO |

### 8.2 The Closest Constraint

**GW250114 pSEOBNR analysis** provides the tightest constraint:
- Constrains delta_f_220 in [0.00, 0.04] at 90% credibility
- The eikonal estimate of -4.4% lies just outside this range
- **HOWEVER**: This is NOT a definitive exclusion because:
  1. The exact QNM of the exponential metric is unknown (eikonal is approximate)
  2. The pSEOBNR model assumes Kerr inspiral (circular reasoning if testing the metric)
  3. The ringdown-only analysis (pyRing) gives [-0.13, +0.43], consistent with exponential
  4. The eikonal approximation may not be valid when photon sphere = throat

### 8.3 What Would Definitively Rule It Out

1. **Short answer**: Compute the exact QNM of the exponential metric, then compare with GW250114 pyRing constraints
2. **If exact QNM < -13%**: Already excluded by pyRing
3. **If exact QNM in [-13%, +4%]**: Consistent with current data; wait for ET/CE
4. **If neutron star z measured to <5%**: The 19% difference would be resolvable
5. **If short-delay echoes searched and not found**: Would constrain the reflectivity of the wormhole throat (but does NOT exclude the metric -- the metric can exist without producing detectable echoes if throat reflectivity is low)

### 8.4 Answer to the Key Question

**The exponential metric is NOT ruled out as of March 2026.** It remains viable because:
1. All solar system tests are passed identically
2. EHT shadow measurements have insufficient precision
3. QNM constraints are model-dependent and the exact prediction is unknown
4. The most distinctive prediction (short-delay echoes) has not been searched for
5. Neutron star constraints are indirect and EOS-degenerate

---

## 9. Most Promising Near-Term Test

### Rank-ordered by feasibility and discriminating power:

**1. Short-delay GW echo search in O3/O4 data (MOST ACTIONABLE)**
- Can be done NOW with existing data
- Requires dedicated templates at 1-10 ms delays
- Would test the most distinctive prediction of the exponential metric
- No new observations needed -- just new analysis
- Recommended: Reanalyze GW150914, GW190521, GW250114 with short-delay templates

**2. Exact QNM computation of pure exponential metric (MOST IMPACTFUL THEORY)**
- Single most important calculation for Paper 2
- Would convert the "marginal" GW250114 constraint to a definitive test
- Methods: Direct integration, continued fraction, or time-domain evolution
- The Nath-Sarma (2024) WKB result for the generalized version gives -7.8%, but the pure Papapetrou metric may differ

**3. NICER/eXTP neutron star gravitational redshift (~2029)**
- 19% is the largest fractional difference
- eXTP aims for ~2% M-R precision
- Combined with improved EOS constraints from GW170817-like events
- Key: Need simultaneous M, R, AND spectroscopic z measurement

**4. ngEHT shadow precision (~2030)**
- Phase 2 may achieve <2% shadow size precision
- Would clearly resolve 4.6% difference
- Long timeline

**5. Einstein Telescope QNM spectroscopy (~2035)**
- O(1)% QNM frequency precision
- Would definitively test the 4.4% prediction
- Very long timeline

### Recommended Action Items for Paper 2:

1. **Immediate**: Call for short-delay echo searches in existing LVK data
2. **Short-term**: Compute exact QNMs of the pure exponential metric
3. **In paper**: Emphasize the echo prediction as the most distinctive and currently untested signature
4. **Acknowledge**: The pSEOBNR constraint from GW250114 is approaching the relevant regime

---

## 10. Complete Reference List

### Exponential/Yilmaz Metric Theory
1. Papapetrou, A. (1954). "Eine Theorie des Gravitationsfeldes mit einer Feldfunktion." Z. Phys. 139, 518.
2. Yilmaz, H. (1958). "New Approach to General Relativity." Phys. Rev. 111, 1417.
3. Yilmaz, H. (1971). "New Theory of Gravitation." Phys. Rev. Lett. 27, 1399.
4. Rosen, N. (1973). "A Bimetric Theory of Gravitation." Gen. Relativ. Gravit. 4, 435.
5. Makukov, M.A. & Mychelkin, E.G. (2020). "Triple Path to the Exponential Metric." Found. Phys. 50, 1346. arXiv:2009.08655.
6. Boonserm, P., Ngampitipan, T., Simpson, A. & Visser, M. (2018). "Exponential metric represents a traversable wormhole." Phys. Rev. D 98, 084048. arXiv:1805.03781.
7. Misner, C.W. (1995). "Yilmaz Cancels Newton." arXiv:gr-qc/9504050.
8. Ibison, M. (2007). "Cosmological test of the Yilmaz theory of gravity." arXiv:0705.0080.
9. Ernazarov, K.K. (2025). "The asymptotically Schwarzschild-like metric solutions." arXiv:2511.13471.

### QNMs and Wormhole Properties
10. Nath, P.P. & Sarma, D. (2024). "A new class of traversable wormhole metrics." Eur. Phys. J. C 84, 1063. arXiv:2401.03738.
11. Nath, P.P. & Sarma, D. (2025). "Exponential wormhole, its quasinormal modes, and study in logarithmic f(R) gravity." Ann. Phys. (2025). (ScienceDirect)
12. Cardoso, V., Franzato, E. & Pani, P. (2016). "Is the Gravitational-Wave Ringdown a Probe of the Event Horizon?" Phys. Rev. Lett. 116, 171101.
13. Damour, T. & Solodukhin, S.N. (2007). "Wormholes as black hole foils." Phys. Rev. D 76, 024016.

### EHT Observations
14. EHT Collaboration (2019). "First M87 Event Horizon Telescope Results. I." ApJ Lett. 875, L1.
15. EHT Collaboration (2022). "First Sagittarius A* Event Horizon Telescope Results." ApJ Lett. (multiple papers).
16. EHT Collaboration (2024). "The persistent shadow of the supermassive black hole of M87. I." A&A 686, A154.
17. EHT at 345 GHz (2024). Highest diffraction-limited angular resolution achieved.

### LIGO/Virgo/KAGRA O4 and Echoes
18. LVK Collaboration (2023). "Constraints on the amplitude of gravitational wave echoes from black hole ring-down using minimal assumptions." Phys. Rev. D 108, 064018. arXiv:2302.12158.
19. LVK Collaboration (2025). "GWTC-4.0: Updated Gravitational-Wave Catalog." (August 2025).
20. LVK Collaboration (2025). "GW250114: Observation of Gravitational Waves from a Binary Black Hole Merger with Total Mass ~120 M_sun."
21. LVK Collaboration (2025). "Black Hole Spectroscopy and Tests of General Relativity with GW250114." arXiv:2509.08099.

### NICER and Neutron Stars
22. Miller, M.C. et al. (2019). "PSR J0030+0451 Mass and Radius from NICER Data." ApJ Lett. 887, L24.
23. Miller, M.C. et al. (2021). "The Radius of PSR J0740+6620 from NICER and XMM-Newton Data." ApJ Lett. 918, L28.
24. Choudhury, D. et al. (2024). "The Radius of PSR J0437-4715 from NICER Data."
25. Salmi, T. et al. (2024). "Constraining the dense matter equation of state with new NICER mass-radius measurements." arXiv:2407.06790.

### Pulsar Timing
26. Kramer, M. et al. (2021). "Strong-Field Gravity Tests with the Double Pulsar." Phys. Rev. X 11, 041050.

### GW Memory
27. Inchauspe, H. et al. (2024). "Measuring gravitational wave memory with LISA." Phys. Rev. D 111, 044044. arXiv:2406.09228.
28. (2026). "Detectability of Gravitational-Wave Memory with LISA: A Bayesian Approach." arXiv:2601.23230.

### Future Detectors
29. Punturo, M. et al. (2023). "Landscape of stellar-mass black-hole spectroscopy with third-generation gravitational-wave detectors." Phys. Rev. D 108, 043019. arXiv:2304.02283.
30. (2025). "Einstein Telescope and Cosmic Explorer." arXiv:2505.11033.
31. (2023). "Fundamental Physics Opportunities with the Next-Generation Event Horizon Telescope." arXiv:2312.02130.

### Strong Lensing
32. Manna, T., Rahman, F. & Chowdhury, T. (2023). "Strong lensing in the exponential wormhole spacetimes." New Astron. Rev.

### Solar System Tests
33. Bertotti, B., Iess, L. & Tortora, P. (2003). "A test of general relativity using radio links with the Cassini spacecraft." Nature 425, 374.
34. Will, C.M. (2014). "The Confrontation between General Relativity and Experiment." Living Rev. Relativ. 17, 4.
35. PDG (2024). "Experimental Tests of Gravitational Theory." Review of Particle Physics.

### Additional Relevant Papers
36. Ernazarov, K.K. (2025). "The Yilmaz-Rosen and JNW metric solutions in scalar-EGB 4d gravitational model." arXiv:2510.21625.
37. Graber, J.R. (2017). "Rotating Yilmaz metric." arXiv:1708.05645.
38. Makukov, M.A. & Mychelkin, E.G. (2023). "Rotating exponential metric via Newman-Janis algorithm." arXiv:2301.08118.
39. Robertson, S.L. (2018). "Astrophysical Applications of Yilmaz Gravity Theory." (vixra preprint)

---

## Appendix A: Key Formulas for the Exponential Metric

### Metric (isotropic coordinates)
ds^2 = -exp(-r_s/r) c^2 dt^2 + exp(+r_s/r) [dr^2 + r^2 dOmega^2]

### Photon sphere
r_ph = r_s (isotropic) -- coincides with wormhole throat

### Critical impact parameter (shadow boundary)
b_crit = 2e M ~ 5.437 M (vs 3sqrt(3) M ~ 5.196 M for Schwarzschild)

### ISCO
R_ISCO ~ 3.17 r_s (areal radius) vs 3 r_s for Schwarzschild

### Gravitational redshift
z_exp = exp(r_s/2r) - 1
z_Schw = (1+r_s/4r)^2 / (1-r_s/4r)^2 - 1

### Echo delay time
Delta_t_echo = (2r_s/c) * [(e - Ei(1)) - (e^2/2 - Ei(2))] ~ 4.17 r_s/c

### QNM frequency (eikonal estimate)
f_QNM^exp / f_QNM^Schw = 3sqrt(3)/(2e) ~ 0.956

### Wormhole throat
R_throat = (e/2) r_s ~ 1.36 r_s

### PPN parameters
gamma = beta = 1 (identical to Schwarzschild)

### First metric difference
g_00 expansion coefficient at O((r_s/r)^3): -1/6 (exp) vs -1/4 (Schw)

---

## Appendix B: Ernazarov 2025 Results for General Exponential Metrics

Ernazarov (arXiv:2511.13471, November 2025) systematically studied a family of asymptotically Schwarzschild-like metrics including the exponential metric. Key numerical results for the pure exponential case (l=0 limit):

- Photon sphere radius: r* = M (isotropic coordinates)
- ISCO radius: r_ISCO = 2M (isotropic coordinates)
- Shadow radius (critical impact parameter): b* = M * e ~ 2.718 M (isotropic)

For comparison, Schwarzschild in isotropic coordinates:
- Photon sphere: r* ~ 1.87 M
- ISCO: r_ISCO ~ 2.91 M
- Shadow: b* = 3sqrt(3) M / 2 ~ 2.598 M (NOTE: isotropic b differs from areal b)

Note: The shadow is LARGER for exponential despite the photon sphere being CLOSER -- this is because the stronger gravitational lensing (exp(r_s/r) factor) bends more light around the photon sphere.

---

## Appendix C: Timeline of Observational Milestones

| Year | Observatory/Event | What it constrains | Precision for exp vs Schw |
|------|------------------|-------------------|--------------------------|
| 2019 | EHT M87* shadow | Shadow size | ~17% -- insufficient |
| 2022 | EHT Sgr A* shadow | Shadow size | ~9-10% -- insufficient |
| 2024 | EHT 345 GHz | Shadow resolution | improved but TBD |
| 2025 | GW250114 (O4b) | QNM frequency | ~few % (pSEOBNR), ~tens % (pyRing) |
| 2025 | O4 completes | Echoes, QNMs | Need dedicated short-delay search |
| ~2028 | ngEHT Phase 1 | Shadow size | ~5% -- marginal |
| ~2029 | eXTP launch | NS M-R-z | ~2% -- approaching |
| ~2030 | ngEHT Phase 2 | Shadow size | <2% -- sufficient |
| ~2035 | Einstein Telescope | QNM, echoes | ~1% -- definitive |
| ~2035 | LISA | GW memory, massive BH | Not a useful discriminator |
| ~2040 | Cosmic Explorer | QNM, echoes | <1% -- definitive |

**The decade 2025-2035 is the critical window for testing the exponential metric.**
