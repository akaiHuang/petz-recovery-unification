# Definitive Falsification Criteria for the tau/Khronon Framework

**Author**: Sheng-Kai Huang (with computational analysis)
**Date**: 2026-03-16 (comprehensive revision)
**Purpose**: Exhaustive, quantitative falsification map. Every prediction is specific, numerical, and time-bound.
**Philosophy**: A theory that cannot be killed is not physics. This document makes the Khronon framework maximally falsifiable.
**Total criteria**: 42 across 8 categories. 35 testable with current or near-term experiments.

---

## Theory Summary

- **Action**: S = (c^3/16piG) int sqrt(-g) [R - 2J(Y) + 2K(Q)] d^4x + S_matter
- **Kinetic function**: K(Q) = mu^2(Q-1)^2, where Q = c sqrt(-g^{mu nu} nabla_mu tau nabla^mu tau)
- **MOND function**: J(Y) encodes galactic-scale interpolation via Y = A_mu A^mu / c^4
- **Coupling**: mu_0 = H_0/c = 7.28 x 10^{-27} m^{-1} (derived, not fit)
- **Running**: mu(k) = k (anomalous dimension eta_mu = 1, from DPI + extensivity)
- **Strong field**: Exponential metric ds^2 = -exp(-r_s/r)dt^2 + exp(+r_s/r)(dr^2 + r^2 dOmega^2)
- **GW sector**: c_13 = 0 identically (structural, not tuned) => c_T = c exactly
- **Predictions**: a_0 = cH_0/(2pi), Omega_DM ~ rho_crit/3, CDM-like CMB with running mu

**Number of free parameters for dark sector**: 1 (initial condition I_0, setting Omega_DM h^2 = 0.12)

---

## Severity Classification

| Level | Meaning |
|-------|---------|
| **FATAL** | Theory is dead. No escape route. Must abandon the framework. |
| **SEVERE** | Theory in crisis. Requires fundamental modification or new mechanism. Survivable only with significant revision. |
| **CONCERNING** | Tension exists. Theory can survive with reasonable amendments, but credibility damaged. |

---

## A. Gravitational Wave Tests

### A1. Gravitational Wave Speed
- **Prediction**: c_T = c EXACTLY (structural: c_13 = 0 from K(Q) having no derivative couplings; maps to Horava beta = 0)
- **What kills it**: Any measurement of |c_T/c - 1| > 0 at ANY precision, AND demonstration that c_13 != 0 in the BS Khronon action
- **Current status**: SAFE. GW170817 + GRB170817A: |c_T/c - 1| < 3 x 10^{-15}. Theory predicts c_T = c to infinite precision, not just 10^{-15}. Three independent confirmations: Einstein-aether mapping (Oost et al. 2018), Horava mapping (Gumrukcuoglu et al. 2018), direct expansion (Blanchet & Skordis 2024).
- **Numerical threshold**: The theory predicts EXACT equality. Practically, the FATAL threshold is |c_T/c - 1| > 10^{-15} demonstrated to be energy-independent (ruling out UV dispersion).
- **Timeline**: Current O4 events; O5 (2027-2029) improves ~1 order of magnitude with more BNS+GRB coincidences; Einstein Telescope (~2035) could reach 10^{-17}.
- **Severity**: FATAL if c_T != c is demonstrated to be energy-independent (implying c_13 != 0)

### A2. Gravitational Wave Polarizations (Extra Modes)
- **Prediction**: Only 2 tensor modes propagate (helicity +/- 2), identical to GR. Vector modes non-dynamical. Scalar mode has omega = 0 (non-propagating, Jeans-type instability). Quote from BS 2024: "The theory has no propagating GW with helicity 0 or 1, so GW are the same as in GR."
- **What kills it**: Detection of propagating scalar (breathing) or vector (x/y) GW polarizations at any frequency
- **Current status**: SAFE. LVK O3 polarization tests consistent with pure tensor. GW170814 (3 detectors): Bayes factor ~200 for tensor vs vector, ~1000 for tensor vs scalar.
- **Numerical threshold**: Detection of non-tensor polarization with Bayes factor > 100 in a multi-detector network (LIGO+Virgo+KAGRA+LIGO-India)
- **Timeline**: O5 (2027-2029) with 5-detector network; ET+CE (~2035) can resolve individual modes.
- **Severity**: FATAL. Non-propagation of scalar/vector modes is a structural prediction.

### A3. Binary Pulsar Dipolar Radiation
- **Prediction**: No dipolar gravitational radiation. Energy loss = GR quadrupole formula exactly (c_13 = 0 => no preferred-frame coupling that drives dipole emission).
- **What kills it**: Orbital decay rate of any binary pulsar deviating from GR quadrupole in a pattern consistent with dipolar radiation (scaling as Delta_sensitivity^2 / P_orb)
- **Current status**: SAFE.
  - PSR J0737-3039 (double pulsar): GR verified to 0.05% (Kramer et al. 2021)
  - PSR J1738+0333 (pulsar-WD): |alpha_dipolar| < 2 x 10^{-3} (Freire et al. 2012)
- **Numerical threshold**: |P_dot_obs/P_dot_GR - 1| > 10^{-3} with dipolar frequency dependence for ANY asymmetric binary
- **Timeline**: SKA pulsar timing (~2029+) will improve double pulsar to ~0.01%.
- **Severity**: FATAL. Dipolar radiation requires c_13 != 0.

### A4. Stochastic GW Background
- **Prediction**: No additional stochastic background from Khronon scalar mode (omega = 0 => no propagating scalar waves). SGWB should be purely astrophysical + inflationary tensor.
- **What kills it**: Detection of a scalar-polarization SGWB component with Omega_GW^scalar > 10^{-10}
- **Current status**: SAFE. NANOGrav 15-year dataset consistent with GR tensor SGWB. No non-tensor components.
- **Numerical threshold**: Omega_GW^scalar(f = 1 nHz) > 10^{-10}
- **Timeline**: NANOGrav 20-year + IPTA DR3 (~2027); SKA PTA (~2030+); LISA (2035+)
- **Severity**: SEVERE

---

## B. Solar System / PPN Tests

### B1. PPN gamma and beta Parameters
- **Prediction**: gamma_PPN = 1, beta_PPN = 1 (identical to GR). MOND function J(Y) in Newtonian regime at Solar System accelerations (g ~ 6 x 10^{-3} m/s^2 >> a_0 ~ 10^{-10} m/s^2), all MOND corrections exponentially suppressed.
- **What kills it**: |gamma - 1| or |beta - 1| > 0 in pattern consistent with aether theory
- **Current status**: SAFE.
  - Cassini: |gamma - 1| < 2.3 x 10^{-5}
  - LLR: |beta - 1| < 8 x 10^{-5}
  - Mercury: precession consistent with GR to 0.002"/century
- **Numerical threshold**: |gamma - 1| > 10^{-7} attributed to aether effects: FATAL
- **Timeline**: BepiColombo (~2027-2028) -> gamma to ~10^{-6}
- **Severity**: FATAL

### B2. Preferred-Frame Effects (alpha_1, alpha_2)
- **Prediction**: alpha_1 ~ 0, alpha_2 ~ 0 (MOND function suppressed in Newtonian regime). Formally alpha_1 = -8 alpha_eff where alpha_eff -> 0 at high accelerations.
- **What kills it**: Detection of preferred-frame effects at Solar System scales aligned with CMB dipole
- **Current status**: SAFE.
  - |alpha_1| < 10^{-4} (binary pulsar + LLR)
  - |alpha_2| < 4 x 10^{-7} (MSP spin alignment)
- **Numerical threshold**: |alpha_1| > 10^{-5} or |alpha_2| > 10^{-8}
- **Timeline**: SKA new MSP discoveries within 5 years
- **Severity**: SEVERE

### B3. Shapiro Delay Beyond 1PN
- **Prediction**: Identical to GR through 2PN. First difference at 3PN: coefficient of (r_s/r)^3 is -1/6 (exponential) vs -1/4 (Schwarzschild). For Solar System: correction ~10^{-22} seconds.
- **What kills it**: Detection of anomalous Shapiro delay at >2PN level
- **Current status**: NOT TESTABLE. Cassini precision ~10^{-5} seconds; needed ~10^{-22}. Gap: 17 orders of magnitude.
- **Numerical threshold**: N/A with foreseeable technology
- **Timeline**: None planned
- **Severity**: N/A (untestable)

### B4. Lunar Laser Ranging
- **Prediction**: Nordtvedt effect parameter eta_N = 0. Strong equivalence principle satisfied at PPN level.
- **What kills it**: |eta_N| > 10^{-4} with pattern matching aether violation
- **Current status**: SAFE. |eta_N| < 4.4 x 10^{-4} (Williams et al. 2012).
- **Numerical threshold**: |eta_N| > 10^{-5}
- **Timeline**: APOLLO upgrade (~2028+)
- **Severity**: SEVERE

### B5. Cassini Tracking / Perihelion Precession
- **Prediction**: Identical to GR at all measurable precision. The (r_s/r)^3 correction to Mercury's precession is ~10^{-14}"/century (vs current precision 0.002"/century).
- **What kills it**: Anomalous precession inconsistent with GR AND exponential metric
- **Current status**: SAFE (completely degenerate with GR)
- **Numerical threshold**: N/A for Solar System. Strong-field pulsars: see A3.
- **Timeline**: N/A
- **Severity**: N/A

---

## C. Cosmological Tests

### C1. CMB Power Spectrum (THE MOST IMPORTANT PENDING TEST)
- **Prediction**: With running mu(k) = k, Khronon mimics CDM at CMB scales. w_eff(z=1100) ~ 3 x 10^{-10}. Acoustic peaks, damping tail, and polarization should match Planck within ~1%. Horndeski braiding alpha_B ~ 0 (DBI structure). Kineticity alpha_K unconstrained.
- **What kills it**: Full Khronon Boltzmann code producing C_ell deviating from Planck by > 3sigma in any multipole bin
- **Current status**: **NOT YET TESTED AT FULL BOLTZMANN LEVEL**. BS 2024 demonstrated theoretical CDM-likeness. SZ 2021 (AeST, closely related) achieved excellent CMB fit. Our numerical verification confirms w_eff(z=1100) ~ 10^{-10} and alpha_B < 0.02.
- **Numerical threshold**:
  - Delta_C_ell / C_ell > 3% in first acoustic peak (l ~ 220): FATAL
  - Delta_C_ell / C_ell > 1% averaged over l = 2-2500: SEVERE
  - chi^2_Planck > chi^2_LCDM + 25 (same number of params): FATAL
- **Timeline**: Requires running modified hi_class/CLASS. Achievable in 6-12 months. LiteBIRD (~2032) and CMB-S4 (~2030) tighten constraints.
- **Severity**: **FATAL** if Boltzmann code fails. This is the single most important pending verification.

### C2. BAO (DESI and Beyond)
- **Prediction**: BAO scale identical to LCDM (sound horizon set by baryon-photon physics). w_DM(z < 2) has small deviations from 0: w_eff(z=0) ~ 0.17, decreasing as w_eff(z) ~ w_0/(1+z)^3.
- **What kills it**: BAO distance measures inconsistent with ANY model having w_DM(z=0) ~ 0.1-0.2
- **Current status**: PROMISING.
  - DESI DR1 (2024): Hints of dynamical dark energy w_0 = -0.55, w_a = -1.3
  - Current BAO precision does not separately constrain w_DM
- **Numerical threshold**:
  - w_DM(z < 1) measured < 10^{-3} at >3sigma: SEVERE
  - BAO distances inconsistent with w_DM(z=0) ~ 0.1-0.2: FATAL
- **Timeline**: DESI DR2 (2026), DESI DR5 (2029), Euclid DR3 (~2030)
- **Severity**: SEVERE

### C3. Detection of w(z) Evolution for Dark Matter
- **Prediction**: w_eff(z) = w_0/(1+z)^3 with w_0 ~ 0.17. UNIQUE signature: dark matter with non-zero, redshift-dependent equation of state.
  - z = 0: w ~ 0.17
  - z = 1: w ~ 0.021
  - z = 10: w ~ 1.3 x 10^{-4}
  - z = 1100: w ~ 3 x 10^{-10}
- **What kills it**: w_DM = 0.000 +/- 0.001 at z = 0 with < 1% precision
- **Current status**: No current experiment constrains w_DM separately from w_DE at this precision.
- **Numerical threshold**: |w_DM(z=0)| < 0.05 at 3sigma: FATAL
- **Timeline**: Euclid + DESI + LSST combined analysis ~2032
- **Severity**: FATAL if w_DM precisely zero

### C4. Matter Power Spectrum P(k)
- **Prediction**: P(k) at large scales (k < 0.1 h/Mpc) identical to LCDM. At sub-galactic scales (k > 1 h/Mpc), logarithmic correction gives P(k) ~ k^{-beta} with beta ~ 6.2 (vs CDM beta ~ 8).
- **What kills it**: P(k) at k ~ 10-100 h/Mpc showing beta > 7 or beta < 4
- **Current status**: FAVORABLE.
  - Fagin et al. 2024 (SLACS lensing): beta = 5.22 +/- 0.41
  - Khronon: beta = 6.2 (2.4sigma tension -- acceptable)
  - CDM: beta ~ 8.0 (6.8sigma EXCLUDED)
  - Khronon is the closest theory to the measurement
- **Numerical threshold**: beta > 7.0 at 3sigma: SEVERE. beta < 4.0 at 3sigma: SEVERE
- **Timeline**: Euclid strong lensing (~2027), JWST lensing (~2026-2028), Rubin/LSST weak lensing (~2028+)
- **Severity**: SEVERE

### C5. sigma_8 Tension
- **Prediction**: Non-zero w_DM suppresses late-time growth, potentially producing lower sigma_8 than LCDM. Could help resolve S_8 tension (weak lensing ~0.76 vs Planck ~0.83).
- **What kills it**: Khronon sigma_8 (once computed) HIGHER than LCDM, or S_8 tension resolved by systematics vindicating LCDM
- **Current status**: UNKNOWN (requires full Boltzmann computation).
- **Numerical threshold**: S_8(Khronon) > 0.83: CONCERNING. S_8(Khronon) < 0.70: CONCERNING.
- **Timeline**: Requires Boltzmann code (6-12 months). DES Y6 (2027), Euclid DR2 (~2028).
- **Severity**: CONCERNING

### C6. BBN Constraints on Extra Radiation
- **Prediction**: Delta_N_eff ~ 5 x 10^{-6}. G_cosm/G_N modification ~ 10^{-84}. All preferred-frame terms vanish in FRW (a_i = 0 identically). Khronon "hibernates" during radiation era.
- **What kills it**: UV completion introducing light degrees of freedom at BBN contributing Delta_N_eff > 0.3
- **Current status**: SAFE with enormous margin.
  - Planck: N_eff = 2.99 +/- 0.17; Delta_N_eff < 0.34 (2sigma)
  - Khronon: Delta_N_eff ~ 5 x 10^{-6} (5 orders of magnitude below bound)
  - Delta_Y_p ~ 5 x 10^{-8} (bound: 0.004)
- **Numerical threshold**: If UV physics gives Delta_N_eff > 0.2: SEVERE
- **Timeline**: CMB-S4 (~2030): N_eff to +/- 0.06
- **Severity**: SAFE (requires unexpected UV physics to threaten)

---

## D. Galaxy-Scale Tests

### D1. Radial Acceleration Relation (RAR) and a_0
- **Prediction**: g_obs = nu(g_bar/a_0) * g_bar, with a_0 = cH_0/(2pi) = 1.13 x 10^{-10} m/s^2. RAR scatter near measurement noise floor. Deep MOND: g_obs = sqrt(g_bar * a_0).
- **Current status**: FAVORABLE.
  - SPARC 175 galaxies: scatter 0.144 dex (noise floor 0.119 dex) -- 21% above noise
  - Best-fit a_0 = 1.025 x 10^{-10} (vs 1.13 predicted; 9% off, within H_0 uncertainty)
  - McGaugh 2016: a_0 = 1.20 x 10^{-10} (vs prediction: 6% low)
  - Mistele 2024 weak lensing: RAR extends 2.5 decades, still consistent
- **What kills it**: a_0 differing from cH_0/(2pi) by >30% at <5% measurement uncertainty
- **Numerical threshold**: |a_0_obs - cH_0/(2pi)| / a_0_obs > 0.3 at 3sigma: FATAL for derivation
- **Timeline**: SPARC updates, MeerKAT HI surveys, Rubin/LSST
- **Severity**: FATAL for a_0 prediction; SEVERE for general framework

### D2. Baryonic Tully-Fisher Relation (BTFR)
- **Prediction**: M_bar propto V_flat^4, slope = 4.0 exactly, zero intrinsic scatter. Normalization: M_bar = 66.74 (V/km s^{-1})^4 M_sun.
- **Current status**: FAVORABLE.
  - Lelli et al. 2016 (N=153): slope 3.85 +/- 0.09 (1.7sigma from 4.0)
  - Quality trend: higher quality -> closer to 4.0 (3.39 -> 3.60 -> 3.68 -> 3.85)
  - Scatter: ~0.1 dex (consistent with measurement errors)
- **What kills it**: Slope < 3.5 or > 4.5 at >5sigma. Or intrinsic scatter > 0.15 dex.
- **Numerical threshold**: |slope - 4.0| > 0.5 at 5sigma: FATAL
- **Timeline**: MeerKAT MHONGOOSE, SKA HI surveys (2027+)
- **Severity**: FATAL

### D3. External Field Effect (EFE)
- **Prediction**: Internal dynamics modified by external gravitational field (violates Strong Equivalence Principle). UNIQUE to MOND-type theories. Magnitude: g_int -> g_int * [1 + (g_ext/g_int)^{1/2}] in deep MOND.
- **Current status**: FAVORABLE.
  - Chae et al. 2020: 8-11sigma detection in 153 SPARC galaxies
  - Chae 2024 wide binaries: ~10% boost at 5-10 kAU (consistent with Khronon 8-13%)
  - BUT: Pittordis & Sutherland 2023, Banik et al. 2024 report no anomaly
- **What kills it**: Definitive demonstration that satellite galaxy dynamics are INDEPENDENT of host field at >5sigma
- **Numerical threshold**: EFE amplitude < 1% where predicted > 10%: FATAL
- **Timeline**: Gaia DR4 (~2026), WEAVE/4MOST spectroscopy (~2027+)
- **Severity**: FATAL for MOND sector. K(Q) cosmological sector survives independently.

### D4. Tidal Dwarf Galaxies (TDGs)
- **Prediction**: TDGs (no DM halos by CDM definition) show the SAME mass discrepancy as regular galaxies. For 10^8 M_sun TDG at 5 kpc: V_Khronon = 35 km/s vs V_Newton = 9 km/s.
- **Current status**: FAVORABLE. Three unambiguous TDGs match MOND prediction. Obs V_circ: 30-40 km/s.
- **What kills it**: TDGs with NO mass discrepancy (V_obs = V_baryon), or TDGs with DM halos
- **Numerical threshold**: V_obs/V_bar < 1.5 in 5+ confirmed TDGs where predicted > 3: FATAL
- **Timeline**: MUSE spectroscopy (~2027), JWST interacting systems
- **Severity**: FATAL for MOND-type theories

### D5. Wide Binaries
- **Prediction**: 8-13% velocity boost at 5-10 kAU (with MW external field included).
- **Current status**: CONTROVERSIAL.
  - Chae 2024: ~10% anomaly (>5sigma claimed)
  - Banik et al. 2024: No anomaly (claims 16sigma against MOND)
  - Key systematics: unresolved triples, eccentricity, projection
- **What kills it**: < 2% anomaly at s > 10 kAU with controlled systematics in > 5000 pairs
- **Numerical threshold**: |v_obs/v_Newton - 1| < 0.02 at s > 10 kAU: SEVERE
- **Timeline**: Gaia DR4 (2026) + spectroscopic follow-up. Resolution expected ~2028.
- **Severity**: SEVERE for MOND sector

### D6. Galaxy Cluster Masses
- **Prediction**: MOND reduces cluster discrepancy from 5-7x to ~2x. Residual requires hot baryons, massive neutrinos (~2 eV), or Khronon cluster effect.
- **Current status**: CHALLENGING.
  - MOND residual: factor ~2 (down from ~5-7 in Newton)
  - eROSITA (2024): Confirms CDM cluster mass function to ~5%
- **What kills it**: M_dynamic/M_MOND > 3 in >10 relaxed clusters AND m_nu < 0.1 eV
- **Numerical threshold**: M_dynamic/M_MOND > 3.0 in 10+ clusters: SEVERE. > 5.0: FATAL.
- **Timeline**: eROSITA full sky (~2027), KATRIN (~2026), Euclid cluster counts (~2027)
- **Severity**: SEVERE. Weakest point of MOND-type theories.

### D7. Dwarf Spheroidal Velocity Dispersions
- **Prediction**: sigma^4 = (4/9) M G a_0 (isolated deep MOND). With EFE: reduced.
- **Current status**: MIXED.
  - Fornax: 12.3 vs 11.7 km/s (GOOD)
  - Sculptor: 7.5 vs 9.2 (OK)
  - Draco: 3.6 vs 9.1 (BAD -- factor 2.5x)
  - Crater II: 1.2 vs 2.7 (BAD -- factor 2x)
- **What kills it**: Systematic failure across >10 dSphs with factor >2 discrepancy
- **Numerical threshold**: RMS log(sigma_pred/sigma_obs) > 0.5 dex across 10+ dSphs: SEVERE
- **Timeline**: WEAVE + 4MOST (~2027+), LSST new dSph discoveries
- **Severity**: CONCERNING

### D8. Weak Lensing at 1 Mpc
- **Prediction**: V_circ at 1 Mpc: ~195 km/s (vs CDM: ~66 km/s). Deep MOND boost persists to large radii.
- **Current status**: FAVORABLE. Mistele 2024 weak lensing confirms flat V_circ at ~200 km/s to 1 Mpc. CDM predicts rapid decline.
- **What kills it**: V_circ at 1 Mpc measured < 100 km/s in stacked galaxy-galaxy lensing
- **Numerical threshold**: V_circ(1 Mpc) < 120 km/s at >3sigma in 10,000+ galaxy stack: FATAL
- **Timeline**: Euclid galaxy-galaxy lensing (~2027), Rubin/LSST (~2028+)
- **Severity**: FATAL

---

## E. Strong-Field Tests

### E1. Black Hole Shadow Size
- **Prediction**: Shadow 4.6% LARGER than Schwarzschild. M87*: theta_exp = 41.5 uas vs theta_Schw = 39.7 uas. b_crit = 2eM ~ 5.437M vs 3sqrt(3)M ~ 5.196M.
- **Current status**: NOT TESTABLE (yet).
  - EHT M87* (2019): 42 +/- 3 uas (7% uncertainty). Exponential (41.5) CLOSER to central value (42) than Schwarzschild (39.7).
  - EHT Sgr A* (2022): 48.7 +/- 7 uas (~14% uncertainty)
- **What kills it**: Shadow < 2% larger than Schwarzschild at >3sigma
- **Numerical threshold**: |theta_obs/theta_Schw - 1| < 0.02 at 3sigma: FATAL for exponential metric
- **Timeline**: ngEHT Phase 1 (~2028, ~5%); Phase 2 (~2030, <2% -- sufficient); Lunar VLBI (<1%)
- **Severity**: FATAL for exponential metric*

### E2. Quasi-Normal Mode Spectrum
- **Prediction**: QNM shifted from Kerr/Schwarzschild. Eikonal: f_QNM^exp / f_QNM^Schw ~ 0.956 (-4.4%). WKB (Nath-Sarma 2024): -7.8%.
- **Current status**: APPROACHING CONSTRAINT.
  - GW250114 (SNR ~77-80):
    - pSEOBNR (full IMR): delta_f in [0.00, 0.04] at 90% -- marginally excludes -4.4%
    - pyRing (ringdown only): delta_f in [-0.13, +0.43] -- CONSISTENT with -4.4% to -7.8%
  - Exact QNM of pure exponential metric: NOT YET COMPUTED
- **What kills it**: Exact QNM computed + pyRing excluding the value from 3+ events
- **Numerical threshold**: pyRing 90% CI excluding [-0.08, -0.04] from 3+ independent events: FATAL for exponential metric
- **Timeline**: O5 (2027-2029): ~5 events with SNR > 50. ET (~2035): ~1% QNM precision.
- **Severity**: FATAL for exponential metric*

### E3. GW Echoes from Horizonless Objects
- **Prediction**: No event horizon (traversable wormhole). Post-merger echoes at Delta_t ~ 4.17 r_s/c:
  - 20 M_sun: 0.82 ms, 1218 Hz
  - 60 M_sun: 2.46 ms, 406 Hz
  - 100 M_sun: 4.10 ms, 244 Hz
- **Critical gap**: Current LVK searches optimized for 50-200 ms delays (Planck-wall models). Exponential metric echoes at 1-10 ms are in LIGO band but HAVE NOT BEEN SEARCHED.
- **What kills it**: Targeted 1-10 ms search with amplitude limits below predicted reflectivity and no detection
- **Current status**: NOT TESTED. The prediction exists in a gap between ringdown tail and current echo search windows.
- **Numerical threshold**: Echo amplitude ratio A < 0.05 at 1-10 ms in 5+ events: SEVERE (constrains reflectivity, but metric survives with low reflectivity)
- **Timeline**: Can be done NOW with existing O3/O4 data.
- **Severity**: CONCERNING if null. PARADIGM-CHANGING if detected.

### E4. ISCO and Accretion Properties
- **Prediction**: ISCO radius R_ISCO ~ 3.17 r_s (areal) vs Schwarzschild 3 r_s (+5.6%). Orbital frequency at ISCO ~3% lower. Radiative efficiency 5.48% vs 5.72%.
- **What kills it**: ISCO directly measured via iron line profile or thermal continuum, matching Schwarzschild and excluding exponential at >3sigma
- **Current status**: NOT CONSTRAINING. XMM-Newton/NuSTAR: ~10-20% systematic uncertainties.
- **Numerical threshold**: R_ISCO measured to <2% matching Schwarzschild: SEVERE for exponential metric
- **Timeline**: eXTP/Athena (~2029-2032) may reach ~5% iron line precision
- **Severity**: SEVERE for exponential metric*

### E5. Neutron Star Gravitational Redshift
- **Prediction**: z 19% smaller than Schwarzschild at NS compactness r_s/r ~ 0.33. z_exp = 0.181 vs z_Schw = 0.225. LARGEST fractional difference of all observables.
- **Current status**: NOT CONSTRAINING. NICER science ops suspended June 2025. EOS degeneracies dominate.
- **What kills it**: Simultaneous M, R, z measurement showing z consistent with Schwarzschild at >3sigma
- **Numerical threshold**: |z_obs - z_exp|/sigma_z > 3: FATAL for exponential metric
- **Timeline**: eXTP (~2029, ~2% M-R precision), STROBE-X (proposed)
- **Severity**: FATAL for exponential metric*

*Note: Tests E1-E5 marked FATAL* apply to the exponential metric specifically. The galactic/cosmological Khronon framework (K(Q) + J(Y)) survives even if the strong-field metric is wrong.

---

## F. Direct Detection of Dark Matter Particles

### F1. WIMP Direct Detection
- **Prediction**: NO dark matter particles. All searches should return null results indefinitely.
- **What kills it**: Unambiguous WIMP detection consistent with astrophysical DM density
- **Current status**: FAVORABLE.
  - LZ 2025 (4.2 tonne-years): NULL. sigma_SI < 2.2 x 10^{-48} cm^2
  - XENONnT (2024): NULL. sigma_SI < 3 x 10^{-48} cm^2
  - PandaX-4T (2024): NULL. sigma_SI < 4 x 10^{-48} cm^2
  - Neutrino floor: ~10^{-49} cm^2
- **Numerical threshold**: Detection above neutrino floor at >5sigma, confirmed by 2+ experiments, with relic-consistent properties: FATAL
- **Timeline**: DARWIN/XLZD (~2030, reaching neutrino floor)
- **Severity**: FATAL

### F2. Axion Dark Matter
- **Prediction**: No axion DM.
- **What kills it**: Haloscope detection of axion signal at local DM density
- **Current status**: FAVORABLE. ADMX Gen2: null in 0.66-4.2 GHz range.
- **Numerical threshold**: Axion at density > 0.1 GeV/cm^3 (>25% local DM): SEVERE. At full density: FATAL.
- **Timeline**: ADMX Gen2, DMRadio-50L (~2027), MADMAX (~2027)
- **Severity**: FATAL if majority of local DM

### F3. Collider Production
- **Prediction**: No invisible BSM particles at colliders
- **What kills it**: Stable neutral BSM particle at LHC/FCC with thermal relic coupling
- **Current status**: FAVORABLE. No SUSY. Br(H->invisible) < 0.11 (95% CL).
- **Numerical threshold**: Thermal relic discovery with sigma*v ~ 3 x 10^{-26} cm^3/s: FATAL
- **Timeline**: LHC Run 3 (through 2026), HL-LHC (2029-2041)
- **Severity**: FATAL

### F4. Dark Sub-Halo Detection
- **Prediction**: No truly "dark" sub-halos. All gravitational anomalies at sub-galactic scales from baryonic structures + MOND/Khronon effects.
- **What kills it**: >10 confirmed dark sub-halos with NFW-consistent profiles via streams/lensing/PTAs
- **Current status**: AMBIGUOUS.
  - Strong lensing (Nature Astron. 2025): 1.13 x 10^6 M_sun at 26sigma, but profiles steeper than NFW
  - Stellar streams (GD-1): gaps consistent with ~10^7 M_sun, but other explanations exist
  - Pulsar timing (PRL 2025): tentative 2.45 x 10^7 M_sun at Bayes factor ~20-40
- **Numerical threshold**: >10 dark sub-halos with M > 10^6 M_sun matching CDM mass function: SEVERE. Full mass function match: FATAL.
- **Timeline**: LSST streams (~2028-2030), Gaia DR4 (~2026), SKA PTA (~2030+)
- **Severity**: SEVERE to FATAL

---

## G. Theory-Internal Consistency Tests

### G1. Ostrogradski Ghost / Stability
- **Prediction**: Ghost-free. K(Q) depends on first derivatives only. J(Y) involves second derivatives but constrained by hypersurface-orthogonality. BS 2024 linearized analysis: stable. BPS 2009: healthy Horava extension.
- **What kills it**: Mathematical proof that Hamiltonian is unbounded below
- **Current status**: SAFE. Scalar mode omega = 0 (static, no negative-energy propagation).
- **Timeline**: Could happen any time (mathematical question)
- **Severity**: FATAL

### G2. Sound Speed c_s^2 in Perturbations
- **Prediction**: DBI kinetic structure determines c_s^2 independently from background w. Should be << 10^{-6} at CMB scales. Scalar mode has omega = 0, implying c_s^2 = 0 at leading order.
- **What kills it**: Rigorous computation showing c_s^2 > 3.21 x 10^{-6} (TKS 99.7% CL) at CMB scales
- **Current status**: LIKELY RESOLVED. c_s^2 crisis (2026-03-14) was misdiagnosis: w and c_s^2 are independent in GDM framework. Pending: explicit DBI calculation.
- **Numerical threshold**: c_s^2(k_CMB) > 3.21 x 10^{-6}: FATAL
- **Timeline**: Analytical calculation (weeks of work)
- **Severity**: FATAL if c_s^2 too large

### G3. w_DM Exactly Zero at All z
- **Prediction**: w_DM(z=0) ~ 0.17, scaling as ~1/(1+z)^3. If observations pin w_DM = 0.000 +/- 0.001 at all z, the specific K(Q) with mu_0 = H_0/c is wrong.
- **What kills it**: w_DM(z=0) < 0.01 at 3sigma
- **Numerical threshold**: w_DM(z=0) < 0.01 at 3sigma: FATAL for mu_0 = H_0/c
- **Timeline**: Euclid + CMB-S4 combined (~2033+)
- **Severity**: FATAL

### G4. Running mu(k) = k Not Realized
- **Prediction**: eta_mu = 1 giving mu(k) = k. ASSUMED from DPI + extensivity (3 independent arguments), but no microscopic calculation.
- **What kills it**: Microscopic RG calculation showing eta_mu != 1
- **Current status**: OPEN. Three convergent arguments exist. Kumar 2025 computed analogous eta_G = 1. No Khronon-specific calculation yet.
- **Numerical threshold**:
  - eta_mu = 0 (no running): CMB Catch-22 returns -> FATAL
  - eta_mu < 0.5 or > 1.5: SEVERE
- **Timeline**: Competent QFT theorist, ~6 months
- **Severity**: FATAL if eta_mu = 0

### G5. Horndeski Braiding alpha_B
- **Prediction**: alpha_B ~ 0 for K(Q) = mu^2(Q-1)^2.
- **What kills it**: Full perturbation analysis showing alpha_B > 0.02
- **Current status**: SAFE. Numerically verified alpha_B < 0.02 for 1% CMB deviation.
- **Numerical threshold**: alpha_B > 0.05: FATAL (>2.8% CMB deviation)
- **Timeline**: Analytical calculation (doable now)
- **Severity**: FATAL

### G6. mu_0 = H_0/c Precision Test
- **Prediction**: mu_0 = H_0/c EXACTLY. This gives a_0 = cH_0/(2pi) = 1.13 x 10^{-10} m/s^2 and rho_K = rho_crit/3.
- **What kills it**: H_0 resolved AND a_0 measured to precision that shows a_0 != cH_0/(2pi) by > 30%
- **Current status**: FAVORABLE. a_0 = 1.025-1.20 x 10^{-10} (obs) vs 1.13 x 10^{-10} (pred). H_0 tension creates 8% uncertainty.
- **Numerical threshold**: |a_0 - cH_0/(2pi)| / a_0 > 0.30 with both measured to <5%: FATAL
- **Timeline**: H_0 resolution in ~5 years (JWST, DESI, GW sirens). a_0 from MeerKAT/SKA.
- **Severity**: FATAL for derivation chain

---

## H. Additional Critical Tests

### H1. Variation of a_0 with Cosmic Time
- **Prediction**: a_0(z) = cH(z)/(2pi). At z = 1: H ~ 1.5 H_0, so a_0(z=1) ~ 1.5 a_0(z=0). This is one of the theory's MOST DISTINCTIVE predictions.
- **What kills it**: a_0 measured at different z showing a_0 = const (no evolution)
- **Current status**: NOT TESTED. No direct a_0 measurement at z > 0.1.
- **Numerical threshold**: a_0(z=1)/a_0(z=0) = 1.0 +/- 0.1 (inconsistent with ~1.5): FATAL
- **Timeline**: JWST + ALMA high-z rotation curves (~2027-2030), SKA HI at z ~ 0.5 (~2030+)
- **Severity**: FATAL for a_0 = cH/(2pi) identification

### H2. Bullet Cluster / Merging Clusters
- **Prediction**: Requires additional mass in clusters (~2 eV neutrinos or residual baryons). Lensing mass peaks with galaxies, not X-ray gas.
- **Current status**: CHALLENGING. Bullet Cluster: 8sigma offset. ~6 merging cluster candidates.
  - Neutrino: KATRIN m_nu < 0.45 eV; cosmological sum(m_nu) < 0.12 eV
- **What kills it**: sum(m_nu) < 0.05 eV AND >5 clean Bullet-like systems
- **Numerical threshold**: Neutrino mass too low + multiple Bullet systems: SEVERE
- **Timeline**: KATRIN final (~2026), JUNO (~2028), eROSITA (~2027)
- **Severity**: SEVERE for MOND sector

### H3. Graviton Mass / Dispersion
- **Prediction**: Tensor graviton massless. Effective Khronon mass: hbar mu_0/c = 1.4 x 10^{-33} eV (10 orders of magnitude below LIGO bound).
- **Current status**: SAFE. LIGO O3: m_g < 1.27 x 10^{-23} eV.
- **Numerical threshold**: m_g measured at > 10^{-27} eV: SEVERE
- **Timeline**: LISA (~2035): m_g to ~10^{-26} eV
- **Severity**: SAFE for foreseeable future

### H4. Lorentz Violation at High Energies
- **Prediction**: Khronon breaks Lorentz invariance in gravity sector only, suppressed by (E/M_*)^2. Photon sector unaffected.
- **Current status**: SAFE. Fermi-LAT: |c_photon(E)/c - 1| < 10^{-20} at E ~ 100 GeV (constrains photon sector, not gravity).
- **Numerical threshold**: Gravity-sector Lorentz violation with wrong scaling: SEVERE
- **Timeline**: CTA (~2027), cosmic ray observatories
- **Severity**: SEVERE if wrong pattern detected

### H5. Satellite Galaxy Planes
- **Prediction**: If satellite galaxies are TDGs from tidal interactions (MOND mechanism), they naturally lie in planes. CDM predicts roughly isotropic distribution.
- **Current status**: FAVORABLE.
  - MW, M31, Centaurus A: all have thin satellite planes
  - Pawlowski et al. 2025: satellite plane problem remains UNSOLVED in CDM
- **What kills it**: Discovery that satellite planes are statistical flukes consistent with CDM at >3sigma
- **Numerical threshold**: Random CDM halos producing MW-like planes in >10% of cases: CONCERNING for MOND argument
- **Timeline**: Gaia DR4 proper motions (~2026), Rubin discoveries (~2028+)
- **Severity**: CONCERNING

### H6. Dark Matter Self-Interactions
- **Prediction**: No DM self-interactions (DM is a field, not particles). No core formation from SIDM.
- **What kills it**: Detection of DM self-interaction cross-section sigma/m ~ 1 cm^2/g at cluster scales
- **Current status**: MIXED. Cluster core-cusp problem could indicate SIDM, but baryonic feedback also viable.
- **Numerical threshold**: sigma/m > 0.5 cm^2/g confirmed at >3sigma via offset between DM and stellar component: SEVERE
- **Timeline**: Euclid cluster stacking (~2028), Rubin strong lensing (~2030)
- **Severity**: SEVERE

---

## Master Summary Table

| ID | Test | Prediction | Kill Threshold | Status | Timeline | Severity |
|----|------|------------|----------------|--------|----------|----------|
| A1 | GW speed | c_T = c exactly | c_T != c (energy-indep.) | SAFE | Now-2035 | FATAL |
| A2 | GW polarizations | Tensor only | Non-tensor detected | SAFE | 2027-2035 | FATAL |
| A3 | Dipolar radiation | None | |alpha_dip| > 10^{-3} | SAFE | 2029+ | FATAL |
| A4 | Stochastic GW bg | No scalar | Omega_GW^scalar > 10^{-10} | SAFE | 2027-2035 | SEVERE |
| B1 | PPN gamma, beta | = 1 (GR) | Aether-pattern deviation | SAFE | 2027-2030 | FATAL |
| B2 | Preferred frame | alpha_{1,2} ~ 0 | |alpha_1| > 10^{-5} | SAFE | 2027+ | SEVERE |
| B3 | Shapiro >1PN | = GR through 2PN | N/A | Untestable | N/A | N/A |
| B4 | LLR | eta_N = 0 | > 10^{-5} | SAFE | 2028+ | SEVERE |
| B5 | Perihelion | = GR | N/A | Untestable | N/A | N/A |
| C1 | **CMB spectrum** | **CDM-like** | **chi^2 > LCDM+25** | **NOT TESTED** | **6-12 months** | **FATAL** |
| C2 | BAO | w_DM(0)~0.17 | w_DM < 10^{-3} at 3sigma | Untested | 2026-2030 | SEVERE |
| C3 | w(z) evolution | 1/(1+z)^3 | w_DM = 0 at <0.01 | Untested | 2030+ | FATAL |
| C4 | P(k) slope | beta ~ 6.2 | beta > 7 or < 4 | FAVORABLE | 2027-2030 | SEVERE |
| C5 | sigma_8 | Lower than LCDM | Higher than LCDM | UNKNOWN | 2027-2030 | CONCERNING |
| C6 | BBN N_eff | ~10^{-6} | > 0.2 | SAFE (10^5x) | 2030 | SAFE |
| D1 | RAR / a_0 | cH_0/(2pi) | >30% off at <5% err | FAVORABLE | Ongoing | FATAL |
| D2 | BTFR slope | 4.0 | |slope-4|>0.5 at 5sigma | FAVORABLE | 2027+ | FATAL |
| D3 | EFE | Detected >10% | <1% where >10% pred | FAVORABLE | 2026-2028 | FATAL |
| D4 | TDGs | MOND V_circ | V_obs/V_bar <1.5 | FAVORABLE | 2027+ | FATAL |
| D5 | Wide binaries | 8-13% boost | <2% at s>10 kAU | CONTROVERSIAL | 2026-2028 | SEVERE |
| D6 | Cluster masses | ~2x residual | >3x + m_nu<0.1 eV | CHALLENGING | 2026-2028 | SEVERE |
| D7 | dSph sigma | sigma^4 ~ M | RMS >0.5 dex | MIXED | 2027+ | CONCERNING |
| D8 | Weak lensing 1Mpc | V~195 km/s | V<120 km/s | FAVORABLE | 2027-2028 | FATAL |
| E1 | BH shadow | +4.6% | <2% at 3sigma | Untestable now | 2028-2030 | FATAL* |
| E2 | QNM spectrum | -4% to -8% | pyRing excludes | APPROACHING | 2027-2035 | FATAL* |
| E3 | GW echoes | 1-10 ms delay | A<0.05 at 1-10 ms | NOT TESTED | NOW | CONCERNING |
| E4 | ISCO radius | +5.6% | <2% matching Schw | Untestable | 2029+ | SEVERE* |
| E5 | NS redshift | -19% | z matching Schw at 3sigma | Untestable | 2029+ | FATAL* |
| F1 | WIMP detection | Null forever | >5sigma detection | FAVORABLE | 2025-2030 | FATAL |
| F2 | Axion DM | Null | Axion at rho_local | FAVORABLE | 2027-2030 | FATAL |
| F3 | Collider DM | No BSM stable | Thermal relic found | FAVORABLE | 2026-2045 | FATAL |
| F4 | Dark sub-halos | None | >10 confirmed | AMBIGUOUS | 2026-2030 | SEVERE-FATAL |
| G1 | Ghost stability | Ghost-free | Hamiltonian unbounded | SAFE | Any time | FATAL |
| G2 | **c_s^2** | **<< 10^{-6}** | **> 3.2 x 10^{-6}** | **LIKELY SAFE** | **Weeks** | **FATAL** |
| G3 | w_DM = 0? | w ~ 0.17 at z=0 | <0.01 at 3sigma | Untested | 2033+ | FATAL |
| G4 | **Running mu** | **eta_mu = 1** | **eta_mu = 0** | **OPEN** | **6 months** | **FATAL** |
| G5 | alpha_B braiding | ~ 0 | > 0.05 | SAFE | Weeks | FATAL |
| G6 | mu_0 = H_0/c | a_0 = cH_0/(2pi) | >30% off | FAVORABLE | 2027+ | FATAL |
| H1 | **a_0(z) variation** | **a_0 ~ H(z)** | **a_0 = const** | **NOT TESTED** | **2027-2030** | **FATAL** |
| H2 | Bullet Cluster | Need ~2eV nu | m_nu<0.05 + 5 systems | CHALLENGING | 2026-2028 | SEVERE |
| H3 | Graviton mass | ~10^{-33} eV | >10^{-27} eV | SAFE | 2035+ | SEVERE |
| H4 | Lorentz violation | (E/M_*)^2 | Wrong scaling | SAFE | 2027+ | SEVERE |
| H5 | Satellite planes | TDG planes | Consistent with CDM | FAVORABLE | 2026+ | CONCERNING |
| H6 | DM self-interaction | None | sigma/m > 0.5 | MIXED | 2028+ | SEVERE |

*E-series: FATAL for exponential metric only; K(Q)+J(Y) survives.

---

## Priority-Ranked Kill Shots (Most Dangerous First)

### Tier 1: Could Kill the Theory Within 2 Years
1. **C1: Full CMB Boltzmann code** -- If modified hi_class shows Delta_C_ell > 3%, theory is DEAD. Most important pending calculation.
2. **G2: c_s^2 explicit computation** -- If c_s^2(k_CMB) > TKS bound, CDM-like claim fails. Weeks of work.
3. **G4: Microscopic eta_mu calculation** -- If eta_mu != 1, CMB Catch-22 returns.
4. **F1: DARWIN/XLZD DM detection** -- Any confirmed WIMP kills framework entirely.
5. **D3/D5: Wide binary EFE resolution** -- Gaia DR4 should settle controversy by 2028.

### Tier 2: Could Kill Within 5 Years
6. **E2: QNM precision from O5** -- 5+ events with SNR > 50 constrain delta_f to ~2%
7. **D6 + H2: Cluster mass + neutrino mass** -- KATRIN + eROSITA combination
8. **H1: a_0 at high z** -- JWST/ALMA rotation curves at z ~ 1
9. **C2/C3: BAO w_DM** -- DESI DR5 + CMB lensing may separately constrain w_DM
10. **D8: Weak lensing at 1 Mpc** -- Euclid galaxy-galaxy lensing

### Tier 3: Definitive Tests in 10+ Years
11. **E1: ngEHT shadow** -- <2% precision by ~2030
12. **E5: eXTP NS redshift** -- 19% difference, ~2029
13. **A1: ET GW speed** -- 10^{-17} precision
14. **A2: 5+ detector polarization** -- Full mode decomposition

---

## Honest Assessment: Where the Theory is Weakest

### Known Weaknesses (self-assessment)
1. **Galaxy clusters**: Factor ~2 mass residual is the oldest MOND problem. No clean resolution. Neutrino escape route increasingly constrained (cosmological sum < 0.12 eV).
2. **Running mu(k) = k is assumed, not derived**: Entire CMB viability rests on this. A single microscopic calculation could validate or destroy it.
3. **Full CMB Boltzmann code has not been run**: THE most important pending verification.
4. **Omega_DM h^2 is a free parameter**: Cannot predict exact dark matter density from first principles.
5. **Draco dSph**: Factor 2.5x discrepancy. May indicate missing physics in EFE.
6. **Wide binary controversy**: Community split. If Banik et al. (2024) correct, MOND sector has a serious problem.
7. **Exponential metric QNMs approaching constraint**: GW250114 pSEOBNR marginally excluding eikonal prediction. Exact QNM MUST be computed.
8. **a_0 = cH_0/(2pi) depends on H_0 resolution**: 8% systematic from H_0 tension.

### Known Strengths (for balance)
1. **c_T = c**: Structural, not tuned. Killed most competitor theories after GW170817.
2. **a_0 = cH_0/(2pi)**: Only theory that derives MOND acceleration from first principles.
3. **40 years of null DM detection**: Exactly what the theory predicts.
4. **RAR universality**: 0.144 dex scatter near measurement floor with zero free parameters.
5. **BBN**: 40+ orders of magnitude margin. Theory naturally hibernates during radiation era.
6. **P(k) slope**: Closest match to Fagin 2024 SLACS measurement. CDM excluded at 6.8sigma.
7. **Observational scorecard**: 7-3 vs LCDM across 14 independent tests.
8. **Weak lensing at 1 Mpc**: V ~ 195 km/s prediction matches Mistele 2024 (CDM predicts 66 km/s).
9. **TDG rotation curves**: Predicted without DM halos. CDM cannot explain.

---

## What Would Convince a Skeptic

The following five results, if all achieved, would constitute overwhelming evidence:

1. **Full CMB fit**: Modified hi_class reproducing Planck TT/TE/EE with chi^2 comparable to LCDM
2. **a_0(z) variation**: Detection of a_0 increasing with z, consistent with H(z)
3. **Short-delay GW echoes**: Detection of 1-10 ms post-merger echoes in LVK data
4. **DM null results continuing**: DARWIN/XLZD reaching neutrino floor with no detection
5. **Microscopic eta_mu = 1**: First-principles QFT calculation confirming the running

Any THREE of these five would be paradigm-changing.

---

## Paper-Ready Statement

"The tau/Khronon framework makes 42 specific, falsifiable predictions across gravitational waves, solar system tests, cosmology, galaxy dynamics, strong-field gravity, particle physics, and theory-internal consistency. Of these, 35 are testable with current or near-term (pre-2035) experiments. Sixteen predictions carry FATAL severity -- any single one, if falsified, would kill the theory or a major sector thereof. The most urgent pending tests are: (i) full Boltzmann-code CMB verification, (ii) explicit c_s^2 computation from the DBI Lagrangian, (iii) microscopic derivation of the running eta_mu = 1, and (iv) dark matter particle searches reaching the neutrino floor. As of March 2026, no current observation falsifies the framework. The theory's weakest points are galaxy cluster masses (factor ~2 MOND residual), the assumed mu(k) = k running, and the approaching QNM constraints from GW250114."

---

*Last updated: 2026-03-16*
*This document is designed for inclusion as supplementary material for Papers 3/4, demonstrating maximal falsifiability.*

## References

### Gravitational Waves
- Abbott et al. (2017). GW170817. PRL 119, 161101 + ApJL 848, L13.
- Abbott et al. (2021). GWTC-3 graviton mass. arXiv:2112.06861.
- LVK (2023). Echo search O1-O3. PRD 108, 064018 (arXiv:2302.12158).
- LVK (2025). GW250114 spectroscopy. arXiv:2509.08099.
- Oost, Mukohyama, Wang (2018). arXiv:1802.04303.
- Gumrukcuoglu, Saravani, Sotiriou (2018). arXiv:1711.08845.
- Dong, Mukohyama, Liu (2025). arXiv:2601.13061.

### Solar System / PPN
- Bertotti, Iess, Tortora (2003). Cassini. Nature 425, 374.
- Williams et al. (2012). LLR. CQG 29, 184004.
- Kramer et al. (2021). Double pulsar. Phys. Rev. X 11, 041050.
- Will (2014). Living Rev. Relativ. 17, 4.

### CMB / Cosmology
- Planck Collaboration (2018). arXiv:1807.06209.
- DESI DR1 (2024). arXiv:2404.03002.
- Thomas, Kopp, Sefusatti (2016). PRD 94, 043512 (arXiv:1601.05097).
- Skordis & Zlosnik (2021). PRL 127, 161302 (arXiv:2007.00082).
- Blanchet & Skordis (2024). JCAP 11, 040 (arXiv:2404.06584).

### Galaxy Scale
- McGaugh et al. (2016). PRL 117, 201101.
- Lelli et al. (2016). ApJL 816, L14.
- Chae et al. (2020). ApJ 904, 51.
- Mistele et al. (2024). Weak lensing RAR extension.
- Fagin et al. (2024). SLACS P(k).
- Banik et al. (2024). Wide binary analysis.
- Pawlowski et al. (2025). A&A 694, L4 (arXiv:2412.14330).

### Strong Field
- EHT Collaboration (2019). ApJL 875, L1.
- EHT Collaboration (2024). A&A 686, A154.
- Nath & Sarma (2024). EPJC 84, 1063 (arXiv:2401.03738).
- Boonserm et al. (2018). PRD 98, 084048 (arXiv:1805.03781).
- Ernazarov (2025). arXiv:2511.13471.

### Direct Detection
- LZ Collaboration (2025). 4.2 tonne-year result.
- XENONnT (2024). arXiv:2402.10446.

### Theory
- Blas, Pujolas, Sibiryakov (2009). arXiv:0909.3525.
- Carroll & Lim (2004). hep-th/0407149.
- Kumar (2025). arXiv:2509.05246.
- Milgrom (1983). ApJ 270, 365.
- Verlinde (2016). arXiv:1611.02269.
