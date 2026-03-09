# Astronomical Observations and Open Datasets for Testing the Exponential Metric

**Date**: 2026-03-09
**Purpose**: Survey of observational data that could distinguish the exponential metric g_00 = exp(-r_s/r) from the Schwarzschild metric g_00 = 1 - r_s/r in the strong-field regime.

---

## Executive Summary

The two metrics agree to first order in r_s/r (weak field) but differ at O((r_s/r)^2) in the strong field. The key question is whether current or near-future astronomical observations can distinguish them. This survey covers six observational channels and assesses each one's discriminating power.

**Bottom line**: Current data cannot yet definitively distinguish the two metrics at the O((r_s/r)^2) level in most channels. However, several near-future observations (ngEHT, BHEX, next-generation gravitational wave detectors, NICER continued operation) will reach the necessary precision within the next 5-10 years. The most promising near-term test is neutron star compactness via NICER, where r_s/r ~ 0.3 gives ~5% metric differences.

---

## Important Caveat: The Yilmaz/Papapetrou Exponential Metric Literature

The exponential metric g_00 = exp(-r_s/r) has been studied before under different names:

- **Yilmaz metric** (Huseyin Yilmaz, 1958): Proposed as a full gravitational theory with an exponential metric. Controversial — Misner published "Yilmaz Cancels Newton" (Nuovo Cimento B 114, 1079, 1999), which was rebutted by Alley & Yilmaz (arXiv:gr-qc/9506082).
- **Papapetrou metric** (1954): Proposed as a solution to Einstein-scalar + Klein-Gordon equations in anti-scalar regime.
- **Rosen bimetric gravity**: Related exponential metric structure.

**Critical distinction for Paper 2**: Paper 2's exponential metric arises from a *different* physical motivation (Sigma_grav = r_s/r from quantum relative entropy loss). It is NOT the full Yilmaz theory (which has its own stress-energy tensor for gravity). Paper 2 proposes the exponential metric as a phenomenological consequence of information loss in the gravitational channel, not as a new field theory of gravity.

**Known objections to the Yilmaz theory** (which may or may not apply to Paper 2):
- Predicts only half the observed light bending at first post-Newtonian order (Ernazarov 2025)
- Predicts one-third the perihelion precession (Ernazarov 2025)
- Different gravitational wave energy loss in binary systems

**Paper 2 must address**: Whether the quantum-channel-derived exponential metric avoids these problems by modifying only g_00 while keeping the spatial metric compatible with GR predictions to 1PN order.

### Key Reference: Ernazarov (2025)
- **arXiv**: 2511.13471, "The asymptotically Schwarzschild-like metric solutions"
- Comprehensive study of exponential metrics, computing photon sphere, ISCO, shadow, and QNMs
- **Photon sphere** (l=0 case): r_ph = M (vs 3M Schwarzschild) — dramatically different
- **ISCO** (l=0 case): r_ISCO = 2M (vs 6M Schwarzschild)
- **Shadow** (l=0 case): critical impact parameter b = Me ~ 2.718M (vs 3sqrt(3)M ~ 5.196M Schwarzschild)
- At l ~ 0.236M, shadow coincides with Schwarzschild
- **Conclusion**: "ruled out by solar system tests" — but this applies to the pure exponential metric with ALL metric components exponential. Paper 2's version may differ.

---

## 1. Event Horizon Telescope (EHT) Data

### 1.1 Current Measurements

**M87* (2017, 2018 observations)**:
- Angular shadow diameter: 42 +/- 3 microarcseconds (uas)
- Mass: ~6.5 x 10^9 M_sun
- Persistent shadow confirmed across years (2024 paper: A&A)
- Brightness peak shifted ~30 degrees between 2017 and 2018 (turbulent accretion)

**Sgr A* (2017 observations)**:
- Angular shadow diameter: 48.7 +/- 7 uas
- Shadow consistent with GR predictions

### 1.2 Can EHT Distinguish the Two Metrics?

**Shadow radius comparison**:
- Schwarzschild: r_shadow = 3*sqrt(3)*M ~ 5.196*M = 5.196 * r_s/2
- For the *pure* exponential metric (Ernazarov): r_shadow ~ 2.718*M (factor of ~2 smaller)
- This is a ~48% difference — easily detectable if the pure exponential metric were the correct one

**However**: Paper 2's exponential metric is g_00 = exp(-r_s/r) but with spatial components that may differ from the pure Yilmaz form. The shadow radius depends on BOTH g_00 and g_rr. If g_rr = 1/g_00 (as in Schwarzschild), the effective potential changes dramatically. If g_rr = exp(+r_s/r), we get the Yilmaz/Papapetrou values. Paper 2 needs to specify the full line element.

**Current EHT precision**: ~7% on shadow diameter (M87*), ~14% (Sgr A*)
- Sufficient to rule out the *pure* Yilmaz exponential metric
- NOT sufficient to detect O((r_s/r)^2) corrections if the leading shadow term matches GR

### 1.3 Open Data Access

**All EHT data is publicly available**:
- **CyVerse Data Commons**: https://datacommons.cyverse.org/
  - M87* 2017 calibrated data: `EHTC_FirstM87Results_Apr2019`
  - M87* 2017 polarized data: `EHTC_M87pol2017_Nov2023`
  - M87* 2018 calibrated data: `EHTC_M87-2018_Mar2024`
  - Sgr A* 2017 calibrated data: `EHTC_FirstSgrAResults_May2022`
- **GitHub repositories**: https://github.com/eventhorizontelescope/
  - `2019-D01-01`: First M87 calibrated data
  - `2019-D01-02`: Imaging pipelines
  - `2022-D02-01`: First Sgr A* calibrated data
- **Data formats**: uvfits, csv, text files, lightcurve-normalized variants
- **EHT data products page**: https://eventhorizontelescope.org/for-astronomers/data
- **ALMA Science Portal**: https://almascience.nrao.edu/alma-data/eht-2021

### 1.4 Future EHT-Related Missions

**Next-Generation EHT (ngEHT)**:
- Enhanced ground-based array with more stations
- Angular resolution: a few uas (vs ~20 uas current)
- Can distinguish alternative metrics when image mismatch exceeds **~4.8%**
- Timeline: Under development, expected operational in late 2020s

**Black Hole Explorer (BHEX)**:
- Space-based VLBI mission (3.4m antenna in medium Earth orbit)
- Angular resolution: <=6 uas at 86-320 GHz
- Can detect **photon rings** (not just the shadow)
- Can distinguish alternative metrics with image mismatch >= **~8%**
- **Launch planned: 2031**
- Antenna technology at TRL 5 (target TRL 6 by 2026)
- **Key capability**: Photon ring detection would test the metric at the photon sphere, exactly where exponential vs Schwarzschild differ most

**Assessment**: ngEHT and BHEX together will achieve percent-level image fidelity, sufficient to test O((r_s/r)^2) deviations for supermassive black holes where r_s/r ~ 1 near the photon sphere.

### References
- EHT Collaboration (2024), "The persistent shadow of M87*", A&A
- Younsi et al. (2025), "The future ability to test theories of gravity with black-hole shadows", Nature Astronomy (arXiv:2511.03789)
- BHEX Mission: https://www.cfa.harvard.edu/research/black-hole-explorer-bhex

---

## 2. Gravitational Lensing

### 2.1 Strong Lensing by Compact Objects

**Key observable**: Deflection angle in the strong-field regime

For the exponential metric, strong lensing has been studied in the "exponential wormhole spacetime" framework:
- ISCO, MBO, and photon sphere radii are ALL less than Schwarzschild values
- But the shadow is LARGER than Schwarzschild (paradoxically, for the Papapetrou form)
- Lensing coefficients and deflection angles differ significantly from Schwarzschild

**Manna, Rahman & Chowdhury (2023)**: "Strong Lensing in the Exponential Wormhole Spacetimes" (SSRN/ScienceDirect)
- Detailed calculation of lensing observables in exponential metric
- Comparison with Schwarzschild lensing coefficients

### 2.2 Weak Lensing / Solar System Light Deflection

Both metrics predict the same first-order deflection angle: Delta_phi = 4GM/(c^2 b)
- Differences appear at second order: O((GM/c^2 b)^2)
- For solar limb grazing: GM/c^2 b ~ 4 x 10^-6 — second-order correction is ~10^-11 radians
- Current astrometric precision (~10 uas ~ 5 x 10^-11 rad) is close but not quite sufficient
- **Future**: 1 nano-arcsecond astrometry would reach this level

### 2.3 Galaxy Cluster Lensing

Galaxy cluster lensing tests are primarily used for dark matter distribution mapping, not metric tests. However:
- Statistical strong lensing has been used to test conformal gravity (arXiv:2506.21019)
- Abell 370 and Abell 2390 used to constrain alternative gravity deflection angles to second order
- **Potential**: If Paper 3 predicts modified lensing at galaxy scales, cluster lensing data could test this

### 2.4 Open Data
- Hubble Space Telescope gravitational lens images: publicly available via MAST
- Euclid, LSST (Vera Rubin Observatory): future surveys will dramatically increase lens sample sizes
- VLBI lens models: ~3 confirmed strong lensing systems currently

### Assessment
**Discriminating power: LOW (weak field) to MODERATE (strong field near compact objects)**
- Weak-field lensing cannot distinguish the two metrics with current technology
- Strong-field lensing near the photon sphere could distinguish them but requires EHT-class resolution
- Galaxy-cluster lensing is more relevant for Paper 3 (dark matter) than Paper 2

### References
- Manna et al. (2023), "Strong lensing in exponential wormhole spacetimes"
- Makarov et al. (2024), "Observational tests in scale invariance II: gravitational lensing" (arXiv:2403.08379)

---

## 3. Gravitational Wave Data (LIGO/Virgo/KAGRA)

### 3.1 Current State of the Art

**GW250114** (January 14, 2025): The loudest gravitational-wave signal ever detected
- **Black hole spectroscopy**: At least 3 QNM modes identified (220 fundamental, 221 first overtone, 440 hexadecapolar)
- **Precision**: Quadrupolar frequency constrained to within **a few percent** of Kerr prediction
- **Alternative gravity constraints**: Coupling length scales constrained below 34-49 km for various modified gravity theories
- **Key result**: "The most stringent single-event verification of GR and the Kerr nature of black holes"

**GWTC-4.0 Catalog** (released August 2025):
- 128 new significant gravitational wave candidates from O4a
- Includes GW231123 — likely highest-mass binary black hole merger observed
- Combined constraints from dozens of events

### 3.2 Can QNM Frequencies Distinguish the Two Metrics?

**QNM frequencies depend on the metric near the photon sphere**. For Schwarzschild:
- omega_n ~ (l + 1/2) * Omega_photon_sphere - i*(n + 1/2)*lambda_Lyapunov

For the exponential metric, the photon sphere is at a different radius, so:
- QNM frequencies would differ significantly from Schwarzschild
- In the eikonal limit, QNM frequency ~ 1/r_ph, and r_ph differs by ~3x between metrics

**Current precision**: QNM frequencies measured to ~10-30% for the loudest events
- This is NOT sufficient to distinguish the two metrics at the O((r_s/r)^2) level
- But it IS sufficient to rule out the pure exponential metric (which predicts ~3x different QNMs)

**METRICS framework** (2025): A new tool for computing QNM frequencies in modified gravity:
- Uses spectral methods with exponential convergence
- Has been applied to Chern-Simons, scalar-Gauss-Bonnet, and axi-dilaton gravity
- Could potentially be applied to the exponential metric

### 3.3 Open Data

**Gravitational-Wave Open Science Center (GWOSC)**: https://gwosc.org/
- **All LIGO/Virgo/KAGRA strain data** publicly available
- O4a data (May 2023 - Jan 2024): released
- GWTC-4.0 parameter estimation: https://zenodo.org/records/16053484
- GWTC-4.0 candidate data: https://zenodo.org/records/16053641
- Tutorials, tools, and software provided
- **Strain time series**: Calibrated data for each detector
- **Parameter estimation posteriors**: Mass, spin, distance, sky location

### 3.4 Future Prospects

**LIGO O5** (expected late 2020s): Improved sensitivity
**LISA** (2030s): Space-based, sensitive to supermassive black hole mergers
**Einstein Telescope / Cosmic Explorer** (2030s-2040s): 10x better sensitivity

With improved SNR, QNM frequency precision will reach ~1%, potentially sufficient to test O((r_s/r)^2) corrections.

### Assessment
**Discriminating power: MODERATE (current) to HIGH (future)**
- Current: Can rule out pure exponential metric but not O((r_s/r)^2) corrections
- Future: Einstein Telescope / Cosmic Explorer could reach the needed precision
- **Key advantage**: QNMs probe the geometry very close to the would-be horizon

### References
- LIGO/Virgo/KAGRA Collaboration (2025), "Black Hole Spectroscopy and Tests of GR with GW250114" (arXiv:2509.08099)
- Wagle et al. (2025), "Probing quadratic gravity with black-hole ringdown" (arXiv:2506.14695)
- GWTC-4.0 documentation: https://gwosc.org/GWTC-4.0/

---

## 4. Shapiro Delay and Solar System Precision Tests

### 4.1 Current Best Constraint

**Cassini spacecraft (Bertotti et al. 2003)**:
- PPN parameter gamma = 1.000021 +/- 0.000023
- Precision: ~2 x 10^-5
- Measured via Shapiro time delay of radio signals passing near the Sun

**Both metrics predict gamma = 1 to first order**. The difference appears at 2PN:
- Schwarzschild: gamma = 1 (exact in isotropic coordinates)
- Exponential metric: gamma = 1 (to 1PN), but differs at 2PN
- The 2PN correction is of order (GM/c^2 R_sun)^2 ~ (2 x 10^-6)^2 ~ 4 x 10^-12
- Current Cassini precision (~10^-5) is **7 orders of magnitude** too coarse

### 4.2 Perihelion Precession

Mercury precession measured to ~0.1% precision:
- Both metrics agree at 1PN (43 arcsec/century)
- 2PN correction: ~10^-4 arcsec/century — far below measurement precision
- **Not useful** for distinguishing the two metrics in the solar system

### 4.3 Future Missions

**BepiColombo** (ESA/JAXA, at Mercury):
- Goal: Measure PPN gamma and beta with improved precision via Mercury Orbiter Radio science Experiment (MORE)
- Expected improvement over Cassini: possibly ~10^-6 on gamma
- Still insufficient for 2PN distinction (~10^-12 needed)

**Future astrometry**: ~1 nano-arcsecond accuracy planned
- 2PN effects in solar system expected to be one order of magnitude below this capability
- **Not sufficient** even with futuristic solar system tests

### 4.4 Galactic Center S-Stars

**S2 and S62 orbits around Sgr A***:
- Future observations can measure 1PN parameters at the Galactic Center
- r_s/r ~ 10^-4 for S2 at pericenter — still weak field
- 2PN effects: ~10^-8, potentially measurable with 30m-class telescopes

### Assessment
**Discriminating power: VERY LOW**
- Solar system tests cannot reach the 2PN level needed to distinguish the metrics
- The strong-field regime (r_s/r ~ 1) is where the metrics diverge, and solar system tests are at r_s/r ~ 10^-6
- Even future solar system missions will not reach the necessary precision
- S-star observations are more promising but still ~4 orders of magnitude short

### References
- Bertotti et al. (2003), Nature 425, 374
- BepiColombo MORE experiment documentation
- Grould et al. (2025), "Future prospects for measuring 1PPN parameters using S2 and S62", A&A

---

## 5. Neutron Star Observations (NICER)

### 5.1 Why Neutron Stars Are the Best Test

Neutron stars have **r_s/r ~ 0.2-0.5**, putting them in the regime where the two metrics differ by ~2-12%:
- At r_s/r = 0.3: g_00(exp) = exp(-0.3) = 0.741 vs g_00(Schw) = 1 - 0.3 = 0.700 → **6% difference**
- At r_s/r = 0.5: g_00(exp) = exp(-0.5) = 0.607 vs g_00(Schw) = 1 - 0.5 = 0.500 → **21% difference**

This means the *surface redshift* of neutron stars differs between the two metrics at the percent level.

### 5.2 Current NICER Measurements

**PSR J0437-4715** (nearest millisecond pulsar):
- Mass: M = 1.418 +/- 0.037 M_sun
- Radius: R = 11.36 (+0.95/-0.63) km
- Compactness: GM/Rc^2 ~ 0.185 → r_s/r ~ 0.37
- **Most stringent EOS constraint** from NICER
- Reference: Choudhury et al. (2024), ApJ 971

**PSR J0030+0451**:
- Compactness: GM/Rc^2 <= 0.171
- Recent re-analysis primarily constrains compactness (not individual M, R)
- Reference: arXiv:2403.14105

**PSR J0740+6620** (most massive known neutron star):
- Mass: ~2.08 M_sun
- Radius: ~12-13 km (similar to less massive pulsars!)
- Compactness: GM/Rc^2 ~ 0.24 → r_s/r ~ 0.48
- At this compactness, the two metrics differ by ~12%

**XTE J1814-338** (extreme compactness candidate):
- Mass: M = 1.21 +/- 0.05 M_sun, Radius: R = 7.0 +/- 0.4 km
- Compactness: GM/Rc^2 ~ 0.26 → r_s/r ~ 0.51
- **If confirmed**: 13% metric difference — potentially detectable
- Currently debated: May require exotic matter (dark matter admixture, strange quark matter)

### 5.3 How to Test with Neutron Stars

The metric choice affects:
1. **Surface gravitational redshift**: z = (g_00)^{-1/2} - 1
   - Schwarzschild: z = (1 - r_s/r)^{-1/2} - 1
   - Exponential: z = exp(r_s/2r) - 1
   - At r_s/r = 0.5: z_Schw = 0.414, z_exp = 0.284 → **30% difference in redshift**

2. **Pulse profile modeling**: NICER fits X-ray pulse profiles accounting for light bending
   - The bending depends on the metric
   - Current analysis assumes Schwarzschild — could be re-analyzed with exponential metric

3. **Maximum mass**: The TOV equation changes with the metric
   - Exponential metric (no horizon) → different maximum mass prediction
   - Test: Does the observed maximum NS mass agree with exponential metric TOV?

4. **I-Love-Q relations**: Moment of inertia, Love number, and quadrupole moment
   - Universal relations depend on the background metric
   - Future GW observations of binary NS can constrain these

### 5.4 Open Data

**NICER data**: Available through NASA's HEASARC archive
- https://heasarc.gsfc.nasa.gov/docs/nicer/
- X-ray timing data for all observed pulsars
- Reduced data products available for key pulsars

**GW170817 tidal deformability**: Available through GWOSC
- Constrains neutron star EOS and compactness

### 5.5 Assessment
**Discriminating power: HIGH — the most promising near-term test**
- Neutron stars are in the strong-field regime (r_s/r ~ 0.3-0.5)
- Metric differences are 5-20% — within reach of current precision
- **Concrete proposal**: Re-analyze NICER pulse profile data assuming exponential metric
- **Challenge**: Disentangling metric effects from equation-of-state uncertainties
- **Key advantage**: Multiple neutron stars with different compactnesses provide a "ladder" of tests

### References
- Choudhury et al. (2024), "A NICER View of PSR J0437-4715", ApJ 971
- Miller et al. (2025), "The Radius of PSR J0437-4715 from NICER Data" (arXiv:2512.08790)
- Riley et al. (2024), "Bulk properties of PSR J0030+0451" (arXiv:2403.14105)
- Kini et al. (2024), XTE J1814-338 measurement

---

## 6. Summary: Observational Channels Ranked by Discriminating Power

| Channel | r_s/r regime | Metric difference | Current precision | Can distinguish? | Timeline |
|---------|-------------|-------------------|-------------------|-----------------|----------|
| **Neutron stars (NICER)** | 0.3-0.5 | 5-20% | ~5-10% on redshift | **Possibly NOW** | Now-2030 |
| **EHT shadow** | ~1 (photon sphere) | ~48% (pure exp) | ~7-14% | Yes (pure exp) | Now |
| **GW ringdown (QNM)** | ~1 (photon sphere) | ~3x (pure exp) | ~10-30% | Yes (pure exp) | Now-2030 |
| **ngEHT / BHEX** | ~1 | few % (if only g_00 differs) | ~2-5% (planned) | **Yes** | 2028-2031 |
| **Future GW detectors** | ~1 | ~% level | ~1% (planned) | **Yes** | 2030s |
| **Shapiro delay** | ~10^-6 | ~10^-12 | ~10^-5 | **No** | Never (solar system) |
| **Perihelion precession** | ~10^-7 | ~10^-14 | ~10^-3 | **No** | Never (solar system) |
| **S-star orbits** | ~10^-4 | ~10^-8 | ~10^-4 | Marginal | 2030s |

---

## 7. Recommended Strategy for Paper 2

### Immediate Actions
1. **Specify the full line element**: Paper 2 must define both g_00 AND g_rr (and g_theta_theta, g_phi_phi). The observational predictions depend critically on whether g_rr = exp(+r_s/r), g_rr = 1/g_00, or something else.

2. **Compute the photon sphere and shadow radius** for Paper 2's specific metric. Compare with EHT M87* constraint (42 +/- 3 uas).

3. **Compute the QNM frequencies** using WKB or the METRICS framework. Compare with GW250114.

4. **Address the Yilmaz objections**: Why does Paper 2's metric not suffer from the "half the light bending" problem? (Answer likely: it's only g_00 that's exponential, not the full metric.)

### Re-Analysis Proposals (for Paper 2 or follow-up)
5. **NICER pulse profile re-analysis**: Use existing NICER data but fit with exponential metric instead of Schwarzschild. This is the most feasible independent test.

6. **EHT shadow fitting**: Use public EHT visibility data to fit ring diameter assuming exponential metric. Compare chi^2 with Schwarzschild fit.

### Future Observations to Propose
7. **BHEX photon ring measurement** (2031+): Photon ring spacing depends sensitively on the metric at the photon sphere.

8. **Binary neutron star GW observations**: Tidal deformability Lambda depends on the metric through the TOV equation.

---

## 8. Key Open Data Portals

| Dataset | URL | Content |
|---------|-----|---------|
| EHT M87* 2017 | https://datacommons.cyverse.org/ | Calibrated visibility data |
| EHT Sgr A* 2017 | https://datacommons.cyverse.org/ | Calibrated visibility data |
| EHT GitHub | https://github.com/eventhorizontelescope/ | Data + imaging pipelines |
| GWOSC | https://gwosc.org/ | All LIGO/Virgo/KAGRA strain data |
| GWTC-4.0 PE | https://zenodo.org/records/16053484 | Parameter estimation posteriors |
| GWTC-4.0 candidates | https://zenodo.org/records/16053641 | Search results and candidates |
| NICER data | https://heasarc.gsfc.nasa.gov/docs/nicer/ | X-ray timing data |
| GW170817 | https://gwosc.org/ | Binary NS merger strain data |

---

## 9. Key Literature for Exponential Metric Observational Tests

### Exponential Metric Theory
- Yilmaz (1958): Original exponential metric proposal
- Papapetrou (1954): Einstein-scalar solution with exponential metric
- Alley & Yilmaz (1995): arXiv:gr-qc/9506082, Rebuttal of Misner's criticism
- Ernazarov (2025): arXiv:2511.13471, Comprehensive analysis of exponential metric observables
- Robertson: "Astrophysical Applications of Yilmaz Gravity Theory"

### EHT and Shadow Tests
- EHT Collaboration (2019): First M87 image papers (I-VI)
- EHT Collaboration (2022): First Sgr A* image papers
- EHT Collaboration (2024): Persistent M87* shadow, A&A
- Younsi et al. (2025): arXiv:2511.03789, Future ability to test gravity with shadows

### Gravitational Wave Tests
- LVK Collaboration (2025): arXiv:2509.08099, GW250114 spectroscopy and GR tests
- Wagle et al. (2025): arXiv:2506.14695, Probing quadratic gravity with ringdown
- GWTC-4.0 documentation: https://gwosc.org/GWTC-4.0/

### Neutron Star Tests
- Choudhury et al. (2024): ApJ 971, PSR J0437-4715 from NICER
- Miller et al. (2025): arXiv:2512.08790, Radius of PSR J0437-4715
- Riley et al. (2024): arXiv:2403.14105, Compactness of PSR J0030+0451
- Kini et al. (2024): XTE J1814-338 extreme compactness

### Solar System Tests
- Bertotti et al. (2003): Nature 425, 374, Cassini Shapiro delay
- Grould et al. (2025): A&A, Future 1PPN parameters from S2 and S62

### Strong Lensing
- Manna et al. (2023): "Strong Lensing in the Exponential Wormhole Spacetimes"
- Makarov et al. (2024): arXiv:2403.08379, Observational tests in scale invariance

---

## 10. Conclusion

The exponential metric g_00 = exp(-r_s/r) makes testable predictions that differ from Schwarzschild at the O((r_s/r)^2) level. The most promising observational tests are:

1. **Neutron star surface redshift and pulse profiles** (NICER data, available now)
2. **Black hole shadow morphology** (ngEHT and BHEX, 2028-2031)
3. **Quasinormal mode frequencies** (next-generation GW detectors, 2030s)

The pure Yilmaz/Papapetrou exponential metric (where ALL metric components are exponential) is already ruled out by solar system tests. However, Paper 2's version — where only g_00 = exp(-r_s/r) arises from quantum information loss — could potentially survive if the spatial metric components remain compatible with GR at 1PN. This distinction must be clearly articulated in Paper 2.
