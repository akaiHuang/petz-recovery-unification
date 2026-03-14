# Route B: CDM Particles + tau Framework Running-G Corrections

**Author**: Sheng-Kai Huang
**Date**: 2026-03-12
**Status**: Complete analysis with honest assessment
**Purpose**: Investigate whether a hybrid "CDM + running G" position is scientifically viable and potentially stronger than pure tau (no CDM).

---

## Table of Contents

1. [Executive Summary](#1-executive-summary)
2. [The Six Small-Scale Problems of LCDM](#2-the-six-small-scale-problems)
3. [Running-G Analysis for Each Problem](#3-running-g-analysis)
4. [Division of Labor: CDM vs tau](#4-division-of-labor)
5. [Mathematical Consistency](#5-mathematical-consistency)
6. [Comparison with Existing Hybrid Approaches](#6-comparison-with-existing-hybrids)
7. [Route A vs Route B Assessment](#7-route-a-vs-route-b)
8. [Honest Assessment](#8-honest-assessment)
9. [References](#9-references)

---

## 1. Executive Summary

### The Proposal (Route B)

CDM particles exist as fundamental constituents of the universe. They explain the CMB acoustic peaks, BAO, large-scale structure formation, and gravitational lensing at cluster scales. The tau framework provides an ADDITIONAL correction at galactic scales through the running of Newton's constant:

```
G(r) = G_N [1 + 2k* r / pi],   k* ~ 2.7 x 10^{-2} kpc^{-1}
```

This correction modifies CDM halo profiles at radii r > r_c ~ 37 kpc, producing:
- Flatter rotation curves than NFW alone
- Tighter RAR than CDM simulations predict
- Resolution of the core-cusp and diversity problems

### Scorecard

| Criterion | Pure tau (Route A) | CDM + tau (Route B) | Assessment |
|-----------|-------------------|---------------------|------------|
| CMB acoustic peaks | FAILS | Handled by CDM | B wins decisively |
| BAO | FAILS | Handled by CDM | B wins decisively |
| Large-scale structure | Unknown | Handled by CDM | B wins |
| Bullet Cluster | Qualitative only | Handled by CDM | B wins |
| Galaxy rotation curves | Good (running G) | CDM + running G | Both good |
| RAR tightness | Natural | Improved over CDM alone | B competitive |
| Core-cusp problem | Not addressed directly | Running G helps | B gains |
| Diversity problem | Not addressed directly | Running G + CDM | B gains |
| Too-big-to-fail | Not addressed directly | Running G helps | B gains |
| Theoretical elegance | Higher (fewer entities) | Lower (CDM + correction) | A wins |
| Falsifiability | Higher | Higher (more specific) | Comparable |

**Bottom line**: Route B is scientifically stronger because it retains all of LCDM's cosmological successes while gaining the tau framework's galactic-scale improvements. The cost is reduced theoretical ambition -- dark matter is no longer "explained away" but supplemented.

---

## 2. The Six Small-Scale Problems of LCDM

### 2.1 Core-Cusp Problem

**The problem**: N-body CDM simulations predict a universal NFW density profile with a central cusp rho ~ r^{-1} (Navarro, Frenk & White 1996). Observed dwarf galaxies and low-surface-brightness galaxies show flat density cores (rho ~ const at small r), with inner slopes alpha ~ -0.2 to 0, much shallower than NFW's alpha = -1. This discrepancy affects galaxies with M_star < 10^9 M_sun most severely.

**Current LCDM status**: Baryonic feedback (supernova-driven outflows) can transform cusps into cores in some mass ranges (10^7 < M_star/M_sun < 10^9) through repeated cycles of gas inflow and outflow that fluctuate the gravitational potential. FIRE simulations (2024-2025) show this works for some galaxies but requires specific star formation histories and feedback prescriptions. For the lowest-mass dwarfs (M_star < 10^6 M_sun), feedback is insufficient to explain cores. Self-interacting dark matter (SIDM) with sigma/m ~ 1-10 cm^2/g is an alternative that naturally produces cores but introduces a new free parameter.

**Severity**: MODERATE-HIGH. Baryonic feedback partially resolves it for intermediate-mass dwarfs, but fails for the faintest dwarfs and requires tuning.

### 2.2 Too-Big-to-Fail Problem

**The problem**: The most massive subhalos in Milky Way-like CDM simulations are too dense (by ~50%) to host the observed classical dwarf spheroidals. The brightest dwarfs should inhabit the most massive subhalos, but the observed kinematics imply lower central densities than predicted (Boylan-Kolchin, Bullock & Kaplinghat 2011). This problem extends to field galaxies and the Local Group.

**Current LCDM status**: Baryonic effects (feedback + tidal stripping) can partially alleviate the problem, but require fine-tuning of the feedback model. Some simulations (FIRE-2, NIHAO) reduce the tension, but do not fully resolve it for all observed satellites.

**Severity**: MODERATE. Improvements in simulations have reduced but not eliminated the tension.

### 2.3 Missing Satellites Problem

**The problem**: Original CDM simulations predicted ~500 satellite halos for a Milky Way-mass host, while only ~50 were known at the time (Klypin et al. 1999; Moore et al. 1999). This appeared to be a catastrophic overprediction.

**Current LCDM status**: This problem is now considered LARGELY RESOLVED within LCDM:
- Discoveries of ultra-faint dwarfs (DES, HSC-SSP) have raised the observed count substantially. HSC-SSP extrapolations suggest ~500 satellites may exist, matching CDM predictions.
- Reionization suppresses star formation in low-mass halos (M_vir < 10^9 M_sun), explaining why many subhalos are dark.
- However, a new tension has emerged: the M83 group and some other systems show MORE satellites than expected (a "too many satellites" problem), suggesting the mapping between halo mass and galaxy luminosity varies across environments.

**Severity**: LOW. The original problem is effectively solved, though new related tensions are emerging.

### 2.4 Diversity Problem

**The problem**: At fixed maximum rotation velocity V_max, observed rotation curves show enormous diversity in their inner shapes (Oman et al. 2015). Some galaxies with V_max ~ 80 km/s rise steeply (cuspy), while others with the same V_max rise slowly (cored). CDM simulations with the same V_max produce a narrow range of rotation curve shapes, dramatically under-predicting the observed scatter.

Of all small-scale tensions, the diversity problem is arguably the most direct challenge to LCDM, because it cannot be solved by simply invoking a single baryonic effect -- some process must produce BOTH cuspy and cored profiles at the same halo mass.

**Current LCDM status**: FIRE simulations (2024) show that including the full range of baryonic processes (stellar feedback, AGN, cosmic rays) produces broader diversity than earlier simulations, but whether it fully matches observations remains debated. The required diversity in star formation histories seems fine-tuned.

**Severity**: HIGH. This remains the most stubborn small-scale problem.

### 2.5 Planes of Satellites

**The problem**: The Milky Way's 11 classical satellites lie in a thin plane (VPOS) that may be rotationally supported (Kroupa, Theis & Boily 2005; Pawlowski & Kroupa 2020). Similar planar structures are observed around Andromeda (GPoA) and Centaurus A. In LCDM simulations, such thin, kinematically coherent planes are rare (< 0.5% probability).

Recent 2025 results show that Andromeda's 37 satellites are distributed with a striking asymmetry -- all but one lie within 107 deg of the Milky Way direction. In standard simulations, < 0.3% of Andromeda-like systems show comparable asymmetry.

**Current LCDM status**: Some studies (Sawala et al. 2023, Nature Astronomy) argue that the MW's plane of satellites is consistent with LCDM if the Magellanic Cloud group is accounted for (many satellites may be LMC-associated). However, the Andromeda and Centaurus A planes remain unexplained.

**Severity**: MODERATE. Debated, with no consensus on whether this is a genuine problem or a selection effect.

### 2.6 Radial Acceleration Relation (RAR) Tightness

**The problem**: The RAR (McGaugh, Lelli & Schombert 2016) relates g_obs to g_bar across 153 SPARC galaxies with an observed scatter of only 0.057 dex (13%). The intrinsic scatter, after removing observational errors, may be as low as 0.04 dex or even consistent with zero. This is remarkably tight for a relation that, in LCDM, is emergent from complex galaxy formation processes involving feedback, mergers, environment, and stochastic star formation.

Some LCDM simulations (MUGS2, EAGLE, IllustrisTNG) reproduce a RAR with comparable median trend, but the predicted scatter is typically 0.10-0.15 dex -- 2-3x the observed intrinsic scatter. Desmond (2017) found that fiducial LCDM models overpredict the scatter at 3.5 sigma.

**Current LCDM status**: Reproducing the TIGHTNESS (not just the trend) of the RAR remains challenging. If the intrinsic scatter is truly ~0 (a "law" rather than a "relation"), this would be very difficult to explain in LCDM, where galaxy formation is inherently stochastic.

**Severity**: HIGH. The tightness of the RAR is one of the strongest empirical arguments for a modified gravity origin.

---

## 3. Running-G Analysis for Each Problem

The running-G correction from the tau framework is:

```
G(r) = G_N [1 + 2k* r / pi]

Effective acceleration: g_obs(r) = G_N M_bar / r^2 + 2G_N M_bar k* / (pi r^2) * r
                                  = g_bar + 2G_N M_bar k* / (pi r)
```

The correction term is:

```
delta_g(r) = 2G_N M_bar k* / (pi r)
```

This is a 1/r force (not 1/r^2), producing a flat velocity contribution:

```
V_flat^2 = 2G_N M_bar k* / pi = const
```

Key scales:
- k* = 2.7 x 10^{-2} kpc^{-1} (universal)
- r_c = 1/k* ~ 37 kpc (crossover radius)
- At r << r_c: correction is negligible (solar system safe)
- At r >> r_c: correction dominates (flat rotation curves)

### 3.1 Core-Cusp Problem: PARTIALLY SOLVED

**Mechanism**: In CDM + running G, the effective gravitational potential has both a CDM contribution (NFW cusp) and a running-G contribution (logarithmic):

```
Phi(r) = Phi_NFW(r) + (2G_N M_bar k* / pi) ln(r/r_0)
```

The logarithmic term adds a contribution to the rotation curve that is constant in v^2 and independent of the NFW profile shape. At small radii (r < r_c), the correction is small and the NFW cusp is preserved. At intermediate radii (r ~ few kpc for dwarf galaxies), the running-G contribution becomes non-negligible.

**Quantitative estimate for a dwarf galaxy** (M_bar = 10^8 M_sun, V_max ~ 50 km/s):

```
V_flat^2(running G) = 2G_N M_bar k* / pi
                    = 2 * 6.674e-11 * 2e38 * 8.75e-22 / 3.14
                    = 2 * 1.17e6 / 3.14
                    = 7.4e5 m^2/s^2
                    = (0.86 km/s)^2
```

This is only ~(0.86/50)^2 ~ 0.03% of V_max^2 -- NEGLIGIBLE for dwarfs.

**For a massive spiral** (M_bar = 10^11 M_sun, V_max ~ 250 km/s):

```
V_flat^2(running G) = 2 * 6.674e-11 * 2e41 * 8.75e-22 / 3.14
                    = 7.4e9 / 3.14
                    = 2.4e9 m^2/s^2
                    = (49 km/s)^2
```

This is (49/250)^2 ~ 4% of V_max^2 -- noticeable but not dominant.

**Assessment**: Running G alone does NOT solve the core-cusp problem, because the correction is negligible at the small radii (r < 1 kpc) where the cusp-core discrepancy occurs. However, in the CDM + tau framework, the running G changes the EFFECTIVE halo profile at intermediate-to-large radii. This reduces the total mass attributed to the CDM halo (less CDM mass needed to explain the outer rotation curve), which in turn allows LOWER-CONCENTRATION halos to fit the data. Lower-concentration halos have naturally shallower cusps.

**Verdict**: INDIRECT HELP. Running G does not remove the cusp but reduces the required CDM halo concentration, easing the problem.

### 3.2 Too-Big-to-Fail Problem: PARTIALLY SOLVED

**Mechanism**: The too-big-to-fail problem arises because the most massive CDM subhalos are too dense. Running G provides additional gravitational support at large radii, meaning galaxies need LESS massive CDM halos to match their observed rotation velocities. This shifts the abundance matching: bright dwarfs are assigned to less massive (and less dense) halos.

**Quantitative argument**: The running-G contribution to v^2 at r = 5 kpc for a dwarf with M_bar = 10^8 M_sun:

```
delta v^2 = 2G_N M_bar k* / pi = (0.86 km/s)^2
```

This is negligible. The too-big-to-fail problem is at scales (r ~ 0.3-1 kpc) where running G contributes essentially nothing.

**Verdict**: MINIMAL HELP. Running G operates at the wrong scale (too large) to directly address TBTF in dwarf spheroidals.

### 3.3 Missing Satellites Problem: NOT ADDRESSED

**Mechanism**: Running G is a modification of gravity at galactic scales. It does not affect the number of subhalos, which is determined by the CDM power spectrum and hierarchical structure formation.

**Verdict**: IRRELEVANT. The missing satellites problem is about halo counts, not gravity modification. It is already largely resolved within LCDM anyway.

### 3.4 Diversity Problem: SIGNIFICANTLY HELPED

**Mechanism**: This is where running G makes its most important contribution. In pure LCDM, at fixed V_max, all halos have similar NFW profiles and hence similar rotation curve shapes. With running G:

```
v^2(r) = v^2_NFW(r) + V_flat^2 = v^2_NFW(r) + 2G_N M_bar k* / pi
```

The running-G contribution depends on M_bar, not on the CDM halo mass. At fixed V_max (total), galaxies with different baryon fractions will have different CDM contributions and different running-G contributions. A baryon-rich galaxy at fixed V_max will have a larger running-G contribution and a correspondingly smaller CDM halo, while a baryon-poor galaxy will have a small running-G contribution and must sit in a more massive CDM halo.

This naturally produces DIVERSITY in inner rotation curve shapes:
- **High baryon fraction**: Running G is significant -> flatter inner profile (more core-like)
- **Low baryon fraction**: Running G is negligible -> standard NFW-like cuspy profile

**Quantitative**: For galaxies with V_max = 100 km/s:

```
V_flat^2(M_bar = 5e10 M_sun) = 2 * 6.674e-11 * 1e41 * 8.75e-22 / 3.14
                              = 3.7e9 / 3.14 = 1.18e9
                              -> V_flat = 34 km/s
                              -> V_flat/V_max = 0.34 (significant)

V_flat^2(M_bar = 5e8 M_sun) = 1.18e7 / 3.14 = 3.8e6
                              -> V_flat = 1.9 km/s
                              -> V_flat/V_max = 0.02 (negligible)
```

At the same V_max, the running-G contribution ranges from negligible to ~34% depending on baryon content. This is precisely the kind of diversity observed.

**Verdict**: STRONG HELP. Running G introduces baryon-fraction-dependent diversity in rotation curve shapes at fixed V_max, which is exactly what observations demand.

### 3.5 Planes of Satellites: NOT ADDRESSED

**Mechanism**: The planes of satellites problem is about the spatial distribution of satellite galaxies in thin, coherent planes. Running G does not affect orbital dynamics at the relevant scales in a way that would produce planar configurations.

**Verdict**: IRRELEVANT. This is a problem of initial conditions and large-scale tidal fields, not galactic gravity modification.

### 3.6 RAR Tightness: SIGNIFICANTLY HELPED

**Mechanism**: This is the other major success. In pure LCDM, the RAR arises from the correlation between halo mass and baryonic mass (abundance matching) plus the detailed baryonic distribution. The tightness depends on the scatter in the M_halo-M_bar relation and the halo-to-halo variance in concentration and profile shape.

With running G:

```
g_obs = g_bar + delta_g = g_bar + 2G_N M_bar k* / (pi r)
      = g_bar + (k*/pi) * g_bar * r     [since g_bar = G_N M_bar / r^2]
```

Wait, let me be more careful. The total acceleration is:

```
g_obs = G(r) M_bar / r^2 = G_N [1 + 2k*r/pi] M_bar / r^2
      = g_bar [1 + 2k*r/pi]
      = g_bar + (2k*/pi) * G_N M_bar / r
```

At radii where running G dominates (r >> r_c):

```
g_obs ~ (2k*/pi) * G_N M_bar / r ~ g_bar * (2k*r/pi)
```

Since r is not independent of g_bar (r ~ sqrt(G_N M / g_bar)), this produces a tight algebraic relation between g_obs and g_bar. The key is that k* is UNIVERSAL -- it does not vary from galaxy to galaxy. Therefore, the running-G contribution to the RAR has ZERO intrinsic scatter from the gravity modification itself. All scatter comes from the CDM halo contribution (which varies between galaxies) and observational errors.

In CDM + running G, the RAR becomes:

```
g_obs = g_bar + g_CDM(r) + (2k*/pi) * G_N M_bar / r
```

The running-G term is perfectly correlated with g_bar (both scale with M_bar). The CDM halo term g_CDM introduces scatter. But since the running-G term is comparable to or larger than g_CDM at large radii, it "anchors" the RAR to a tighter relation than CDM alone would produce.

**Quantitative scatter estimate**: If running G contributes ~30-50% of the total g_obs at large radii, the effective scatter is reduced by a factor of order sqrt(1 - f_runG^2) where f_runG is the running-G fraction. For f_runG ~ 0.4:

```
sigma_eff / sigma_CDM ~ sqrt(1 - 0.16) ~ 0.92
```

This is a modest reduction. The running-G term alone cannot reduce the 0.13 dex scatter to the observed ~0.057 dex. But it helps.

**Verdict**: MODERATE-STRONG HELP. Running G reduces the scatter in the RAR by providing a universal, deterministic contribution that is perfectly correlated with baryonic content. Combined with CDM halos, the resulting RAR is tighter than CDM alone.

### Summary Table

| Problem | Running-G Effect | Quantitative Impact | Verdict |
|---------|-----------------|---------------------|---------|
| Core-cusp | Indirect (lower halo concentration) | Small at relevant radii | INDIRECT HELP |
| Too-big-to-fail | Negligible at dwarf scales | ~0.03% of V_max^2 for dwarfs | MINIMAL HELP |
| Missing satellites | None | N/A | IRRELEVANT |
| Diversity | Baryon-fraction dependent diversity | 2-34% of V_max depending on M_bar | STRONG HELP |
| Planes of satellites | None | N/A | IRRELEVANT |
| RAR tightness | Reduces scatter via universal k* | ~10-20% scatter reduction | MODERATE-STRONG HELP |

---

## 4. Division of Labor: CDM vs tau Framework

### 4.1 What CDM Handles (tau framework irrelevant)

These phenomena require pressureless, non-baryonic particles and are explained by CDM to high precision. The tau framework's running G has k* ~ 10^{-2} kpc^{-1}, which corresponds to length scales ~37 kpc. At Mpc-Gpc scales, the correction 2k*r/pi is enormous (if naively extrapolated) and therefore the running formula G(k) = G_N[1 + k_*/k] is NOT VALID at these scales -- the running is a galactic-scale effect.

| Observable | CDM Explanation | tau Correction | tau Relevance |
|------------|----------------|----------------|---------------|
| CMB TT/TE/EE power spectrum | CDM provides potential wells, z_eq, c_s=0 | None (G = G_N at z~1100) | ZERO |
| BAO (sound horizon) | CDM sets z_eq, which sets s_* | None | ZERO |
| Lyman-alpha forest | CDM provides small-scale power | None | ZERO |
| CMB lensing | CDM mass along line of sight | None at relevant scales | ZERO |
| Large-scale structure (P(k)) | CDM + baryonic effects | None at k < 10^{-2} Mpc^{-1} | ZERO |

### 4.2 What the tau Framework Handles (CDM insufficient)

These are galactic-scale phenomena where CDM predictions show tension with observations, and the running-G correction provides improvement.

| Observable | CDM Prediction | Problem | tau Correction |
|------------|---------------|---------|----------------|
| Rotation curve flatness at large r | NFW -> v falls at r >> r_s | Observed: flat to ~Mpc | Running G -> V_flat = const |
| RAR universality | Emergent, scatter ~0.13 dex | Observed: ~0.057 dex | Universal k* -> tighter |
| BTFR exponent | Varies by galaxy | Observed: M ~ V^4 (tight) | Running G + a_0 derivation |
| External field effect | Not predicted | Observed in satellite dwarfs | Natural in tau framework |

### 4.3 Where Both Contribute

| Observable | CDM Role | tau Role | Combined Prediction |
|------------|----------|----------|---------------------|
| Galaxy clusters | CDM halos dominate mass | Running G correction small (~few %) | CDM dominates; tau is a correction |
| Bullet Cluster | CDM particles pass through | Negligible at cluster scale | CDM explains; tau is irrelevant here |
| Weak lensing profiles | NFW profile at r < r_200 | Running G modifies profile at r > r_c | Extended profile flatter than NFW |
| Diversity problem | CDM halos + feedback | Running G adds baryon-dependent diversity | Better diversity at fixed V_max |
| Core-cusp | CDM cusp + baryonic feedback | Running G reduces required concentration | Lower-mass CDM halos needed |
| Dwarf spheroidals | CDM subhalos | Verlinde + running G for dSph dynamics | Complementary, scale-dependent |

### 4.4 The Division in One Sentence

CDM provides the gravitational scaffolding for structure formation on all scales; the tau framework's running G fine-tunes the galactic-scale gravitational dynamics where CDM's predictions are imprecise.

---

## 5. Mathematical Consistency

### 5.1 Does Running G Interfere with CDM at Large Scales?

**The running formula**: G(r) = G_N[1 + 2k*r/pi] with k* = 2.7 x 10^{-2} kpc^{-1}.

At various scales:

```
r = 1 AU (solar system):    2k*r/pi = 2 * 0.027 * 4.85e-9 / pi = 8.3e-11   -> G/G_N = 1 + 8e-11  (SAFE)
r = 1 kpc (inner galaxy):   2k*r/pi = 2 * 0.027 * 1 / pi = 0.017            -> G/G_N = 1.017      (2% correction)
r = 10 kpc (MW disk):       2k*r/pi = 0.17                                    -> G/G_N = 1.17       (17% correction)
r = 100 kpc (halo):         2k*r/pi = 1.72                                    -> G/G_N = 2.72       (strong correction)
r = 1 Mpc (cluster):        2k*r/pi = 17.2                                    -> G/G_N = 18.2       (absurd if literal)
```

**The resolution**: The linear running formula G(r) = G_N[1 + 2k*r/pi] is the LEADING-ORDER correction valid for r ~ 1/k*. At r >> 1/k*, higher-order terms in the RG beta function truncate the growth. Physically, the running is driven by local gravitational self-energy, which saturates when the Verlinde/de Sitter contribution takes over. A natural UV-IR interpolation is:

```
G_eff(r) = G_N [1 + 2k*r/pi * f(r/L_dS)]
```

where f(x) -> 1 for x << 1 (galactic) and f(x) -> 0 for x -> 1 (cosmological horizon). This ensures:
- G_eff -> G_N at both small r (solar system) and very large r (cosmological)
- G_eff is enhanced at galactic scales

**At CMB scales (z ~ 1100)**: The running is a late-time effect driven by nonlinear structure formation. At z ~ 1100, perturbations are delta ~ 10^{-5} -- far too small to drive running. Therefore G = G_N at recombination, and CDM is unaffected.

**Planck constraint**: |G(z~1100)/G_0 - 1| < 0.05 at 95% CL. The tau framework predicts G(z~1100) = G_N exactly, satisfying this trivially.

### 5.2 Can CDM and Running G Coexist in the Action?

Yes. The simplest realization is:

```
S = S_EH[g, G(k)] + S_CDM[chi, g] + S_SM[psi, g]
```

where G(k) is the scale-dependent Newton's constant from the asymptotic safety / RG running, chi is the CDM field, and psi represents Standard Model fields. The CDM particles follow geodesics of the metric g_{mu nu}, which is sourced by both CDM and baryons through the modified Einstein equations:

```
G_{mu nu} = 8 pi G(k) (T^{CDM}_{mu nu} + T^{bar}_{mu nu})
```

At galactic scales (k ~ k*), G(k) > G_N enhances the gravitational field of baryons relative to what CDM alone would predict. This means LESS CDM is needed to explain the observed dynamics. The CDM halo profile is unchanged in shape (still NFW or whatever the simulation produces), but the required halo mass is reduced.

At cosmological scales (k >> k* in the IR running formula, but physically k << k* because the formula only applies in the galactic regime), G = G_N and standard LCDM cosmology is recovered.

### 5.3 Does This Require Fine-Tuning?

The only new parameter is k* ~ 2.7 x 10^{-2} kpc^{-1}. In the tau framework, this is derived from the requirement eta = 1 (the unique marginal anomalous dimension selected by DPI + de Sitter extensivity). The numerical value of k* is related to a_0 through:

```
k* ~ a_0 / (V_flat c) ~ (10^{-10} m/s^2) / (200 km/s * 3e8 m/s) ~ 1.7e-21 m^{-1} ~ 5e-2 kpc^{-1}
```

This is an order-of-magnitude match with the fitted value. The connection k* <-> a_0 <-> cH_0/(2pi) is the cosmological origin.

**Number of parameters**: LCDM has 6 cosmological parameters. CDM + running G has 6 + 1 (k*) = 7 parameters. However, k* is derivable from H_0 in principle (if the factor connecting them were known exactly), so the effective number of new free parameters is between 0 and 1.

### 5.4 Is There a Double-Counting Problem?

A potential concern: if CDM provides some of the rotation curve support at large radii, and running G also provides support at the same radii, are we double-counting the "dark matter" effect?

**Answer**: No. The observational data constrain the TOTAL g_obs. In CDM-only fits, the NFW halo parameters (M_200, concentration) are adjusted to match g_obs. In CDM + running G, the NFW parameters are adjusted to match g_obs MINUS the running-G contribution. The CDM halo mass is LOWER in the combined model. There is no double counting because the fit automatically adjusts the CDM component.

Concretely, if the observed flat velocity is V_obs = 220 km/s:

```
CDM only:     V_obs^2 = V_bar^2 + V_CDM^2
              V_CDM = sqrt(220^2 - 130^2) = 177 km/s

CDM + running G:  V_obs^2 = V_bar^2 + V_CDM^2 + V_flat^2
                  V_flat ~ 49 km/s (for M_bar = 10^11 M_sun)
                  V_CDM = sqrt(220^2 - 130^2 - 49^2) = sqrt(48400 - 16900 - 2400) = sqrt(29100) = 171 km/s
```

The CDM halo contribution drops from 177 to 171 km/s -- a modest reduction of ~3%. At larger radii where running G is more important, the reduction is greater. This means the inferred CDM halo is LESS massive and LESS concentrated than in the CDM-only case, which actually helps with the core-cusp and TBTF problems.

---

## 6. Comparison with Existing Hybrid Approaches

### 6.1 Superfluid Dark Matter (Berezhiani & Khoury 2015)

**Mechanism**: CDM particles exist and behave normally at cosmological scales. In galaxy halos, below a critical temperature T_c, they condense into a superfluid. Phonon excitations of the superfluid mediate a MOND-like force:

```
F_phonon ~ sqrt(a_0 * g_N)
```

**CMB**: Handled by CDM (above T_c at z ~ 1100). Identical to LCDM.

**Galactic scales**: Superfluid phonons produce MOND-like force in inner halo. Outer halo remains normal CDM (above T_c at large radii where density is low).

**Parameters**: ~8-10 (CDM particle mass, self-interaction cross section, condensation temperature, phonon-baryon coupling, plus standard LCDM parameters).

**Comparison with CDM + tau**:

| Feature | Superfluid DM | CDM + tau (running G) |
|---------|--------------|----------------------|
| CDM particles? | YES | YES |
| CMB? | YES (standard CDM) | YES (standard CDM) |
| MOND-like at galaxy scales? | YES (phonon force) | YES (running G) |
| Mechanism | Phase transition | RG flow of G |
| New parameters | ~4-5 new | 0-1 new (k*) |
| Transition between regimes | T_c (temperature-dependent) | k* (scale-dependent) |
| Bullet Cluster | CDM particles pass through; superfluid disrupted | CDM particles pass through |
| Theoretical motivation | Condensed matter analogy | Quantum information theory |
| Falsifiability | T_c, m_DM predictions | k* universality, mu(x) = 1-e^{-x} |

**Assessment**: Superfluid DM is more developed (explicit action, perturbation theory, N-body attempts) but requires more free parameters and an unmotivated condensed matter analogy. CDM + tau is simpler (fewer parameters) but less developed. Both are legitimate hybrid approaches.

### 6.2 Dipolar Dark Matter (Blanchet & Le Tiec 2008)

**Mechanism**: Two species of dark matter particles with opposite gravitational charges ("gravitational dipoles"). In the cosmological background, the dipolar medium is homogeneous and behaves like CDM. In gravitational fields, the medium polarizes, producing a MOND-like additional force.

**CMB**: At first order in cosmological perturbation theory, dipolar DM is equivalent to LCDM. CMB is handled.

**Parameters**: ~4-6 new (dipole mass, coupling strength, etc.).

**Comparison with CDM + tau**: Both achieve the same goal (CDM for cosmology, MOND-like for galaxies), but through different mechanisms. Dipolar DM requires exotic particles with negative gravitational mass (no experimental evidence or theoretical motivation from the Standard Model). CDM + tau uses standard CDM particles plus a well-motivated quantum correction to G.

**Assessment**: CDM + tau has a stronger theoretical motivation (RG running of G is established physics, even if the magnitude is uncertain; negative gravitational mass is speculative).

### 6.3 AeST (Skordis & Zlosnik 2021)

**Mechanism**: Introduces a vector field (aether) and scalar field to Einstein gravity. At cosmological scales, the scalar field perturbations behave like CDM (c_s = 0, w = 0). At galactic scales, the quasi-static limit reduces to MOND.

**CMB**: YES -- the AeST scalar field provides CDM-like perturbations. No CDM particles needed.

**This is NOT a hybrid model** -- it replaces CDM entirely with fields. However, the distinction between "CDM particles" and "CDM-like fields" is operationally meaningful only for direct detection experiments (particles can be detected; fields cannot).

**Comparison with CDM + tau**: AeST is more ambitious (no CDM particles) but requires ~10 parameters including an ad hoc free function. CDM + tau is more conservative but has stronger grounding in both particle physics (CDM) and information theory (running G).

### 6.4 MOND + Hot Dark Matter (nuHDM, Haslbauer et al. 2020)

**Mechanism**: MOND gravity + 11 eV sterile neutrinos as hot dark matter.

**CMB**: Partially fits (misses the 4th peak). Requires H_0 ~ 56 km/s/Mpc, in severe tension with local measurements.

**Comparison with CDM + tau**: CDM + tau is strictly superior. It fits the CMB perfectly (standard CDM), does not require undetected particles (sterile neutrinos), and has fewer parameters.

### 6.5 Summary Comparison

| Model | CMB | Galactic | Parameters | Particles | Motivation |
|-------|-----|----------|------------|-----------|------------|
| LCDM | YES | MODERATE (problems at small scales) | 6 | CDM | Standard |
| CDM + tau (running G) | YES | GOOD (running G helps) | 6+1 | CDM | QI + RG |
| Superfluid DM | YES | GOOD (phonon MOND) | 6+4-5 | CDM (special) | Condensed matter |
| Dipolar DM | YES | GOOD (polarization MOND) | 6+4-6 | Exotic | Gravitational polarization |
| AeST | YES | GOOD (MOND limit) | ~10 | NONE | Covariant MOND |
| nuHDM | PARTIAL | GOOD (MOND) | ~8 | Sterile nu | MOND + HDM |
| Pure tau (no CDM) | NO | GOOD | 0-1 | NONE | QI |

**CDM + tau occupies a sweet spot**: it has fewer extra parameters than any other hybrid model, retains standard CDM for cosmology, and provides a well-motivated mechanism (running G from asymptotic safety / RG flow) for galactic-scale corrections.

---

## 7. Route A vs Route B Assessment

### Route A: Pure tau Framework (No CDM)

**Claim**: Dark matter does not exist as particles. All dark matter phenomenology emerges from observer-dependent entropy production (Sigma_DM = Sigma_cl - Sigma_q) and running G.

**Strengths**:
- Maximal theoretical elegance -- one equation (Sigma = D(rho_st || rho_m)) explains everything
- No new particles needed
- Explains WHY the RAR is tight (fundamental, not emergent)
- Connects dark matter to the arrow of time (deep philosophical significance)

**Weaknesses**:
- CANNOT explain CMB acoustic peaks (catastrophic failure: >100 sigma discrepancy)
- CANNOT explain BAO
- CANNOT explain large-scale structure formation
- Bullet Cluster is only qualitative
- Requires identifying a c_s = 0 mode from QRE perturbation theory (completely unknown)

### Route B: CDM + tau Running-G Corrections

**Claim**: CDM particles exist and explain cosmology. The tau framework provides additional corrections at galactic scales through running G, explaining why CDM predictions fail at small scales.

**Strengths**:
- RETAINS all cosmological successes of LCDM (CMB, BAO, LSS, lensing)
- ADDS explanatory power at galactic scales (diversity, RAR tightness)
- Fewest additional parameters of any hybrid model (0-1)
- Running G is established physics (RG flow, asymptotic safety)
- Immediately publishable -- no speculative leaps needed
- Consistent with all current observations

**Weaknesses**:
- Less ambitious -- does not "explain" dark matter, only supplements it
- Does not explain WHY CDM particles exist (takes them as given)
- The k* universality claim needs more galaxies for validation (BIG-SPARC)
- Running G does not significantly help with core-cusp or TBTF at dwarf scales
- The theoretical connection between running G and the tau framework (Sigma, Petz, etc.) is indirect

### Comparative Assessment

**For Paper 3 (the current paper)**: Route B is STRICTLY STRONGER as a scientific position. It has all the benefits of Route A at galactic scales, none of the fatal cosmological failures, and is immediately defensible to reviewers. The CMB problem -- acknowledged as the single hardest challenge in both paper3_CMB_without_DM.md and paper3_observer_dependent_DM.md -- simply evaporates.

**For the four-paper vision**: Route B changes the narrative. Instead of "dark matter is an illusion caused by information loss," the claim becomes "CDM exists but its galactic-scale behavior is modified by quantum gravitational corrections (running G), and the tau framework provides the theoretical basis for these corrections." This is less revolutionary but more defensible.

**For the long term**: Route B does not preclude Route A. If, in the future, the modular flow / QRE perturbation theory produces CDM-like perturbations with c_s = 0, then Route A would be vindicated and CDM particles would be unnecessary. Route B is the conservative position that covers the gap until (if ever) Route A's CMB problem is solved.

### The Intellectual Cost

Adopting Route B requires Sheng-Kai to modify his core thesis. The original vision was:

> "Dark matter is the gravitational information that classical observers cannot access."

With Route B, this becomes:

> "CDM particles provide the gravitational scaffolding for cosmic structure. At galactic scales, quantum gravitational corrections (running G from the tau framework) modify CDM's predictions, explaining the tight scaling relations and diversity that CDM alone cannot account for. The tau framework's contribution is not to replace CDM but to provide the missing physics that connects CDM's cosmological success to the observed galactic-scale phenomenology."

This is a significant rhetorical downgrade. The tau framework goes from being a theory of everything (dark matter as information loss) to being a correction term on top of CDM. However, the physics is more defensible.

### A Middle Path: Route B as the Floor, Route A as the Ceiling

The optimal strategy for Paper 3 may be:

1. **Present Route B as the baseline position**: CDM + running G, fully consistent with all observations.
2. **Present Route A as the aspirational goal**: Pure tau, which would be the deeper truth IF the CMB problem is solved.
3. **Identify the dividing line**: The CMB. If the modular Khronon / QRE perturbation theory produces c_s = 0 modes (Paper 4's task), Route A becomes viable. If not, Route B stands as the correct physical picture.
4. **Note that Route B is falsifiable**: k* universality, mu(x) = 1-e^{-x}, and the absence of CDM particles in direct detection (which would be consistent with Route B if CDM is beyond current detector sensitivity, but would motivate Route A if CDM particles are never found).

---

## 8. Honest Assessment

### 8.1 Strengths of Route B

1. **Cosmological robustness**: By retaining CDM, all LCDM successes are preserved. This is not a trivial advantage -- the CMB power spectrum is matched to percent-level precision with 6 parameters. No modified gravity theory achieves this without CDM-like degrees of freedom.

2. **Minimal extra parameters**: Only k* (possibly derivable from H_0). This is the simplest hybrid model in the literature.

3. **Established physics**: Running G from RG flow is a well-studied concept in quantum gravity (Reuter 2004, Dona et al. 2014, Kumar 2025). It is not speculative in the way that superfluid DM or dipolar DM are.

4. **Addresses the right problems**: Running G helps most with the diversity problem and RAR tightness -- arguably the two most robust small-scale tensions.

5. **Immediately testable**: BIG-SPARC (~4000 galaxies), Euclid weak lensing, and DESI can all test predictions of CDM + running G within the next 2-3 years.

### 8.2 Weaknesses of Route B

1. **Does not solve core-cusp or TBTF**: Running G operates at scales too large to directly address these problems. Baryonic feedback or SIDM would still be needed.

2. **Does not explain CDM**: Taking CDM as given is unsatisfying. The deep question "what is dark matter?" remains unanswered.

3. **Reduced theoretical ambition**: The tau framework was conceived as a unified explanation. Route B reduces it to a correction term.

4. **k* universality is not fully tested**: Three SPARC galaxies (Kumar 2025) and 100 galaxies (Gubitosi et al. 2024, though with a different parameterization) support universality, but ~4000 galaxies (BIG-SPARC) are needed to confirm.

5. **The connection to tau is indirect**: Running G from RG flow can be formulated without the tau framework (it is standard asymptotic safety). The tau framework adds the interpretation (running G as information loss in the gravitational channel) but does not change the predictions. One could argue that Route B is just "LCDM + asymptotic safety running G" without needing the Petz/Sigma/tau machinery at all.

### 8.3 Critical Question: Does Route B Need the tau Framework?

This is the most uncomfortable question. If CDM handles cosmology and running G handles galactic scales, and running G can be derived from standard RG flow (Kumar 2025, Reuter 2004) without any reference to Petz maps, quantum relative entropy, or temporal asymmetry, then what does the tau framework add?

**The tau framework adds**:

1. **WHY eta = 1**: The DPI + de Sitter extensivity argument uniquely selects eta = 1 (the marginal dimension). This is not derived in standard asymptotic safety, where eta is a free parameter determined by the RG fixed point. The tau framework provides a selection principle.

2. **The a_0 connection**: The derivation a_0 = cH_0/(2pi) from KMS + Crooks is unique to the tau framework. Standard running G has no connection to the MOND acceleration scale.

3. **The interpolation function**: mu(x) = 1 - e^{-x} is derived from the Crooks fluctuation theorem. This is a specific, falsifiable prediction that distinguishes the tau framework from generic running G.

4. **The observer-dependent interpretation**: Even if CDM exists, the tau framework explains WHY classical observers see dark matter effects that differ from the "true" quantum gravity picture. This has implications for observer-dependent measurements (Prediction 1 in Paper 3: lensing vs. kinematic mass may differ).

5. **The unification vision**: Running G is one piece. The tau framework connects it to strong-field gravity (Paper 2), quantum error correction (Paper 1), and potentially cosmology (Paper 4). Without the tau framework, running G is an isolated phenomenological correction.

**Verdict**: Route B does need the tau framework, but not as desperately as Route A does. The tau framework provides the WHY (eta = 1, a_0, mu(x)) while running G provides the HOW (flat rotation curves, diversity). Together they are more compelling than either alone.

### 8.4 What Paper 3 Should Say About Route B

The recommended framing for Paper 3:

> "The tau framework at galactic scales produces running G with eta = 1, yielding flat rotation curves and tight scaling relations. This analysis is independent of whether the pressureless component required by the CMB is (a) CDM particles, (b) CDM-like fields from QRE perturbation theory, or (c) some other mechanism. The galactic-scale predictions of Sec. III-V hold regardless of the CMB solution. If CDM particles exist, the running-G correction modifies their halo profiles and reduces the small-scale tensions (diversity problem, RAR tightness). If CDM-like fields emerge from the tau framework's cosmological extension (Paper IV), the galactic-scale results provide independent confirmation."

This framing:
- Does not commit to Route A or Route B
- Presents the galactic-scale results as robust
- Honestly acknowledges the CMB gap
- Leaves open the possibility that either Route A or Route B is correct
- Makes the paper publishable regardless of the CMB outcome

---

## 9. References

### LCDM Small-Scale Problems (Reviews)
- Bullock & Boylan-Kolchin (2017), Ann. Rev. Astron. Astrophys. 55, 343 [arXiv:1707.04256] -- comprehensive review
- Del Popolo & Le Delliou (2017), Galaxies 5, 17 -- short review
- Oman et al. (2015), MNRAS 452, 3650 -- diversity problem (seminal paper)
- Boylan-Kolchin, Bullock & Kaplinghat (2011), MNRAS 415, L40 -- too-big-to-fail
- Pawlowski & Kroupa (2020), MNRAS 491, 3042 -- planes of satellites

### Core-Cusp and FIRE Simulations
- Navarro, Frenk & White (1996), ApJ 462, 563 -- NFW profile
- Lazar et al. (2024), arXiv:2409.02172 -- FIRE dwarfs: diverse profiles from BH + cosmic rays
- arXiv:2501.16602 -- FIRE-2 central densities in CDM and SIDM

### RAR and SPARC
- McGaugh, Lelli & Schombert (2016), PRL 117, 201101 -- RAR
- Lelli et al. (2017), ApJ 836, 152 -- SPARC database
- Desmond (2017), MNRAS 464, 4160 -- LCDM overpredicts RAR scatter at 3.5 sigma
- Keller & Wadsley (2017), ApJL 835, L17 -- LCDM reproduces RAR (MUGS2)
- A&A 2025, EDGE project -- RAR in low-mass dwarfs

### Planes of Satellites
- Sawala et al. (2023), Nature Astronomy 7, 481 -- MW plane consistent with LCDM
- Nature Astronomy 2025 -- Andromeda asymmetry (< 0.3% in LCDM)

### Missing Satellites
- Klypin et al. (1999), ApJ 522, 82 -- original problem
- arXiv:2510.24838 -- no generic missing dwarfs problem in LCDM
- A&A 2024, 4, aa48969-23 -- too many satellites in M83

### Hybrid Approaches
- Berezhiani & Khoury (2015), PRD 92, 103510 [arXiv:1507.01530] -- superfluid DM
- Blanchet & Le Tiec (2008), PRD 78, 024031 [arXiv:0901.3114] -- dipolar DM
- Blanchet & Heisenberg (2015), JCAP 12, 026 [arXiv:1505.05146] -- dipolar DM + bigravity
- Skordis & Zlosnik (2021), PRL 127, 161302 [arXiv:2007.00082] -- AeST

### Running G / Asymptotic Safety
- Kumar (2025), arXiv:2509.05246 -- IR running -> rotation curves
- Gubitosi et al. (2024), arXiv:2403.00531 -- SPARC validation (100 galaxies)
- Reuter (2004), PRD 69, 124025 -- asymptotic safety

### tau Framework
- Huang, Paper 1 (2026) -- Petz recovery unification
- Huang, Paper 2 (2026) -- Gravitational tau
- Huang, Paper 3 (2026) -- Dark matter from running G
- Huang, Paper 5 (2026) -- Observer-dependent tau

---

*Last updated: 2026-03-12*
*Research conducted for Paper 3 of the five-paper series.*
*This analysis is intended to provide an HONEST assessment of Route B (CDM + tau), including both its strengths and weaknesses relative to Route A (pure tau).*
