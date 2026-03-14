# Can the JRSWW Bound Predict Omega_DM h^2?

**Author**: Sheng-Kai Huang (with computational analysis)
**Date**: 2026-03-12
**Status**: NEGATIVE RESULT -- honest assessment
**Relevance**: Paper 4 (grand unification), Paper 3 (dark matter)

---

## Executive Summary

**VERDICT: This path is a dead end for predicting Omega_DM h^2 from first principles.**

We attempted to use the JRSWW bound (F >= exp(-Sigma/2)) at the cosmological horizon to predict the dark matter density parameter Omega_DM h^2 ~ 0.12. After exhaustive analysis of 10 distinct approaches, the conclusion is unambiguous:

- **No clean derivation exists** connecting JRSWW to Omega_DM.
- **The closest numerical coincidence** (volume-averaged tau with linear Sigma giving 0.309 vs Omega_m = 0.314, 1.5% discrepancy) is unmotivated and constitutes p-hacking.
- **The fundamental obstacle** is that JRSWW is a universal information-theoretic bound that knows nothing about the specific matter content of the universe. Omega_DM is a contingent fact requiring specific physics (freeze-out, baryogenesis, etc.) that no universal bound can reproduce.

The tau framework's value for dark matter lies in **reinterpretation** (dark matter = inaccessible information, as in Paper 3), not in **prediction** of numerical values.

---

## Table of Contents

1. [The Setup](#1-the-setup)
2. [Approach 1: Sigma = -ln(g_00) at Fixed Radius](#2-approach-1)
3. [Approach 2: tau = Omega_DM/Omega_total](#3-approach-2)
4. [Approach 3: Numerology with Universal Constants](#4-approach-3)
5. [Approach 4: Padmanabhan's N_sur vs N_bulk](#5-approach-4)
6. [Approach 5: Cumulative Sigma at Hubble Scale](#6-approach-5)
7. [Approach 6: Bekenstein Bound](#7-approach-6)
8. [Approach 7: Sigma = 1 at Horizon (from paper3_cosmological_tau)](#8-approach-7)
9. [Approach 8: Trivedi-Scherrer Holographic DM](#9-approach-8)
10. [Approach 9: Padmanabhan N_bulk/N_sur Ratio](#10-approach-9)
11. [Approach 10: Baryogenesis Entropy Production](#11-approach-10)
12. [The Best Numerological Coincidence](#12-best-coincidence)
13. [Fundamental Obstruction: Why This Cannot Work](#13-fundamental-obstruction)
14. [What the tau Framework CAN Do About Dark Matter](#14-what-tau-can-do)
15. [Comparison with Existing Information-Theoretic Approaches](#15-comparison)
16. [Recommendations](#16-recommendations)
17. [References](#17-references)

---

## 1. The Setup

### 1.1 The JRSWW Bound

The Junge-Renner-Sutter-Wilde-Winter bound (Annales Henri Poincare 19, 2955, 2018):

```
F(rho, R_Petz(N(rho))) >= exp(-Sigma/2)
```

where:
- F = fidelity between original state and recovered state
- R_Petz = the Petz recovery map
- N = quantum channel
- Sigma = D(rho || sigma) - D(N(rho) || N(sigma)) = decrease in quantum relative entropy

In the tau framework: tau = 1 - F, so:

```
tau <= 1 - exp(-Sigma/2)
```

### 1.2 The Question

Can we identify:
- N = cosmological horizon channel (traces out trans-horizon DOF)
- Sigma = entropy production at/of the horizon
- tau_horizon = irrecoverable information fraction
- Omega_DM = tau_horizon (or some function thereof)

and thereby PREDICT Omega_DM h^2 ~ 0.12?

### 1.3 Known Cosmological Parameters

| Parameter | Value | Source |
|-----------|-------|--------|
| H_0 | 67.4 km/s/Mpc = 2.184 x 10^-18 s^-1 | Planck 2018 |
| Omega_Lambda | 0.685 | Planck 2018 |
| Omega_DM | 0.265 | Planck 2018 |
| Omega_b | 0.049 | Planck 2018 |
| Omega_r | 9 x 10^-5 | CMB + BBN |
| r_H = c/H_0 | 1.37 x 10^26 m | derived |
| T_H = hbar H/(2 pi k_B) | 2.65 x 10^-30 K | Gibbons-Hawking |
| S_H = A/(4 l_P^2) | ~2.27 x 10^122 | Bekenstein-Hawking |

---

## 2. Approach 1: Sigma = -ln(-g_00) at Fixed Radius

### Setup
For static de Sitter: ds^2 = -(1 - H^2 r^2/c^2) dt^2 + ...

Using Paper 2's ansatz Sigma_grav = -ln(-g_00):

```
Sigma_dS(r) = -ln(1 - r^2/r_H^2)
```

### Results

| r/r_H | Sigma | F >= exp(-Sigma/2) | tau <= |
|--------|-------|-------------------|--------|
| 0.01 | 0.0001 | 0.9999 | 0.0001 |
| 0.10 | 0.0101 | 0.9950 | 0.0050 |
| 0.30 | 0.0943 | 0.9539 | 0.0461 |
| 0.50 | 0.2877 | 0.8660 | 0.1340 |
| 0.70 | 0.6733 | 0.7141 | 0.2859 |
| 0.90 | 1.6607 | 0.4359 | 0.5641 |
| 0.99 | 3.9170 | 0.1411 | 0.8589 |

To get tau = Omega_DM = 0.265, we need r/r_H = 0.678.

### Assessment
**DEAD END.** No physical principle selects r/r_H = 0.678 as special. The identification is completely arbitrary.

---

## 3. Approach 2: tau = Omega_DM/Omega_total

### Setup
Identify tau_horizon = Omega_DM/Omega_total = 0.265.
This requires Sigma = -2 ln(1 - 0.265) = 0.616.

### Assessment
**DEAD END.** This is the answer, not a derivation. We need to DERIVE Sigma = 0.616 from horizon physics, which we cannot do.

---

## 4. Approach 3: Numerology with Universal Constants

### Checks

| Expression | Value | Compare with |
|------------|-------|-------------|
| 1 - 1/e | 0.632 | Omega_DM = 0.265 |
| 1/e | 0.368 | Omega_DM = 0.265 |
| exp(-1/2) | 0.607 | Omega_DM = 0.265 |
| 1 - exp(-1/2) | 0.394 | Omega_m = 0.314 (25% off) |
| exp(-Omega_Lambda/2) | 0.710 | -- |
| 1 - exp(-Omega_Lambda/2) | 0.290 | Omega_DM = 0.265 (9.5% off) |

The last one (1 - exp(-Omega_Lambda/2) = 0.290 vs 0.265) is numerically interesting but CIRCULAR: it uses Omega_Lambda as input to "predict" Omega_DM.

### Assessment
**DEAD END.** No simple universal-constant expression matches Omega_DM.

---

## 5. Approach 4: Padmanabhan's Holographic Equipartition

### Framework (arXiv:1206.4916, arXiv:1703.06144)
Padmanabhan's expansion law:

```
dV/dt = L_P^2 (N_sur - N_bulk)
```

where N_sur = surface DOF on Hubble sphere, N_bulk = bulk DOF in thermal equilibrium with Hubble temperature.

At equilibrium (de Sitter endpoint): N_sur = N_bulk, which determines Lambda.

### Assessment
**NOT APPLICABLE.** Padmanabhan's framework predicts Lambda (the cosmological constant), not Omega_DM. Matter is the non-equilibrium excess that drives expansion. The framework does not partition matter into dark and baryonic components.

---

## 6. Approach 5: Cumulative Sigma at Hubble Scale

### From paper3_cosmological_tau.md, Route B

```
Sigma_cumulative ~ (4 pi G rho_0 / 3) r^2/c^2 = (r/r_H)^2 x Omega_m
```

At r = r_H: Sigma ~ Omega_m = 0.314
tau = 1 - exp(-0.314/2) = 0.145

### Assessment
**DEAD END (circular).** This uses Omega_m as input. The entire point of the calculation is that Sigma ~ Omega_m at Hubble scale -- the matter content determines the entropy production, not the other way around.

---

## 7. Approach 6: Bekenstein Bound

### Setup
- Total energy in observable universe: E_total ~ rho_crit c^2 V_H ~ 8.3 x 10^69 J
- Bekenstein bound: S_Bek = 2 pi E R / (hbar c) ~ 2.27 x 10^122
- Holographic bound: S_holo = A/(4 l_P^2) ~ 2.27 x 10^122
- Actual entropy (dominated by SMBHs): S_actual ~ 3.1 x 10^104

### Results
- S_actual/S_holo ~ 10^-18

### Assessment
**DEAD END.** The ratio 10^-18 is far too small to have anything to do with Omega_DM ~ 0.265. The universe is nowhere near saturating its Bekenstein/holographic bound.

---

## 8. Approach 7: Sigma = 1 at Horizon

### From paper3_cosmological_tau.md
The ansatz Sigma_dS = r/r_H (by analogy with Sigma_BH = r_s/r for black holes) gives:

```
Sigma = 1 at r = r_H
tau(r_H) = 1 - exp(-1/2) = 0.394
```

### Comparison
- Omega_DM = 0.265 (48% off from tau)
- Omega_m = 0.314 (25% off from tau)

### Assessment
**SUGGESTIVE NUMEROLOGY at best.** The number 0.394 is in the right ballpark of Omega_m = 0.314, but:
1. 25% discrepancy is not a "prediction"
2. The ansatz Sigma = r/r_H is from analogy, not derivation
3. The identification tau(r_H) = Omega_m has no physical justification

---

## 9. Approach 8: Trivedi-Scherrer Holographic DM (arXiv:2511.10617)

### Their Mechanism
Trivedi and Scherrer (Nov 2025, Vanderbilt) propose holographic dark matter where the IR cutoff is the Ricci scalar:

```
rho_holo = 3c^2 / (8 pi G L^2)
L^(-2) = (2/c^2)(H_dot + 2H^2)    [Ricci cutoff]
```

In a universe with only baryons and radiation, this automatically generates a dark-matter-like component. The model also reverses the sign of a pre-existing negative vacuum energy (relevant for string landscape).

### Connection to tau Framework
The holographic DM density arises from the information capacity of the cosmic horizon, which resonates with the tau reinterpretation. However:
- It is NOT derived from JRSWW
- The Ricci cutoff is a phenomenological choice
- The connection to tau would require showing that the Ricci cutoff = Sigma_horizon, which has not been done

### Assessment
**INDIRECT CONNECTION only.** The most promising existing approach to holographic DM, but fundamentally different from a JRSWW derivation.

---

## 10. Approach 9: N_bulk/N_sur Ratio

### Setup
In Padmanabhan's framework: N_bulk,matter / N_sur = Omega_m

If we define tau = N_bulk,m / (N_sur + N_bulk,m):

```
tau = Omega_m / (1 + Omega_m) = 0.239
```

### Assessment
**CIRCULAR.** This is just the Friedmann equation in disguise. The "tau" here is not independently determined.

---

## 11. Approach 10: Baryogenesis Entropy Production

### Check
- Omega_DM/Omega_b = 5.41
- eta_b (baryon-to-photon ratio) = 6.1 x 10^-10
- ln(1/eta_b) = 21.22
- ln(Omega_DM/Omega_b) = 1.688

### Assessment
**DEAD END.** No clean relationship between the DM/baryon ratio and any information-theoretic quantity. The ratio ~5.4 has no obvious origin from JRSWW or any entropy bound.

---

## 12. The Best Numerological Coincidence

Among all attempts, the closest match found is:

### Volume-Averaged tau with Linear Sigma

Using the de Sitter ansatz Sigma_dS = r/r_H:

```
<tau>_linear = 3 integral_0^1 [1 - exp(-u/2)] u^2 du = 0.3094
```

Compare: **Omega_m = Omega_DM + Omega_b = 0.314**

**Discrepancy: 1.5%**

### Why This Cannot Be Taken Seriously

1. **Unmotivated identification**: Why should the volume-averaged tau equal Omega_m? No physical argument connects "average irrecoverable information fraction inside the Hubble sphere" to "total matter density fraction." These are different concepts.

2. **Ad hoc ansatz**: Sigma = r/r_H is chosen by analogy with Sigma_BH = r_s/r for black holes. The black hole case has a first-principles derivation (modular flow, gravitational Landauer, etc.). The de Sitter case does not -- it is pure analogy.

3. **Predicts Omega_m, not Omega_DM**: Even if taken seriously, this gives the total matter fraction, not the dark matter fraction separately. There is no mechanism to separate Omega_DM from Omega_b.

4. **p-hacking**: We tried ~15 different formulas. Finding one with 1.5% agreement is statistically expected by chance. With N ~ 15 independent trials and each spanning a range of ~2 orders of magnitude, a ~2% match somewhere is likely.

5. **Not a prediction**: The value 0.309 depends on the specific form of Sigma(r). Changing the ansatz (e.g., Sigma = (r/r_H)^2 or Sigma = -ln(1-r^2/r_H^2)) gives completely different answers (0.145 or 0.411 respectively).

**GRADE: Suggestive numerology. Not a derivation. Not publishable as a prediction.**

---

## 13. Fundamental Obstruction: Why This Cannot Work

### The Core Problem: Universality vs Contingency

The JRSWW bound F >= exp(-Sigma/2) is a **UNIVERSAL** statement from quantum information theory. It holds for:
- Every quantum channel N
- Every pair of states (rho, sigma)
- Every Hilbert space dimension

It knows nothing about:
- The specific matter content of the universe
- The baryon-to-photon ratio
- Big Bang nucleosynthesis
- The freeze-out temperature of any particle
- Whether dark matter is WIMPs, axions, or something else

Meanwhile, Omega_DM h^2 ~ 0.12 is a **CONTINGENT** fact about our universe, determined by:
- The DM particle mass and cross-section (in standard LCDM)
- OR the thermal history of the universe
- OR the specific dynamics of whatever mechanism generates DM

**A universal bound cannot determine a contingent parameter.**

### Analogy
The second law of thermodynamics says Delta S >= 0. This CANNOT predict the specific entropy of a glass of water -- you need the heat capacity, temperature, amount of substance, etc.

Similarly, F >= exp(-Sigma/2) CANNOT predict Omega_DM. You need the specific entropy production Sigma of the cosmological channel, which requires knowing the matter content -- the very thing we are trying to predict.

**THE ARGUMENT IS CIRCULAR.**

### The Bound vs The Exact Value

Even if we could determine Sigma_horizon from first principles, the JRSWW bound gives:

```
tau <= 1 - exp(-Sigma/2)    [UPPER BOUND]
```

To predict Omega_DM, we need the EXACT tau, not just a bound. The actual tau could be anywhere from 0 to the upper bound. The bound alone is insufficient unless we can prove saturation.

Is there a reason to expect saturation at the cosmological horizon? **No.** In fact, Petz saturation (tau = 0 exactly) requires the channel to be a sufficient statistic -- a very special condition. Saturation of the JRSWW bound (tau = 1 - exp(-Sigma/2) exactly) has no known mechanism for gravitational channels.

---

## 14. What the tau Framework CAN Do About Dark Matter

### Reinterpretation (Paper 3 approach)

The tau framework provides a rigorous reinterpretation of dark matter:

```
LCDM:         Omega_DM = mass density of unknown particles
tau framework: Omega_DM = information inaccessible to EM observer
              = Sigma_classical - Sigma_quantum
              = gravitational entanglement entropy that classical astronomy cannot access
```

This is physically meaningful and connects to Verlinde's emergent gravity program:
- Volume-law entanglement of the de Sitter vacuum provides "apparent dark matter"
- The tau framework gives this a precise information-theoretic formulation
- Observer dependence (Theorem 1 from observer_dependent_tau.md) explains why different experiments see different amounts of "dark matter"

### What It Does NOT Do

The reinterpretation takes Omega_DM = 0.265 as **observational input** and reinterprets its **meaning**. It does not predict the **value**.

This is analogous to how thermodynamics reinterprets heat as molecular kinetic energy without predicting the specific heat capacity of iron.

---

## 15. Comparison with Existing Information-Theoretic Approaches

### Trivedi-Scherrer (arXiv:2511.10617, Nov 2025)
- **What they do**: Holographic DM from Ricci scalar cutoff
- **Can they predict Omega_DM?**: They generate a DM-like component, but its amplitude depends on the cutoff choice. It's phenomenological.
- **Relation to JRSWW**: None. Different framework entirely.

### Verlinde (arXiv:1611.02269, 2016)
- **What he does**: Apparent DM from volume-law entanglement entropy exceeding area law
- **Can he predict Omega_DM?**: He predicts galaxy rotation curves and the acceleration scale a_0. The cosmological Omega_DM is not predicted from first principles.
- **Relation to JRSWW**: Indirect. The entropy displacement is related to Sigma, but the connection has not been made precise. Ghari & Haghi (arXiv:2601.01715) show Verlinde beats MOND at 5.2 sigma.

### Padmanabhan (arXiv:1703.06144, 2017)
- **What he does**: Predicts Lambda from cosmic information (CosMIn)
- **Can he predict Omega_DM?**: No. His framework determines Omega_Lambda, not Omega_DM.
- **Relation to JRSWW**: CosMIn is related to total integrated tau (see paper3_cosmological_tau.md, Section 4.1).

### Asymmetric Dark Matter / Dynamical Adjustment
- Several particle physics models explain Omega_DM/Omega_b ~ 5 through shared production mechanisms (arXiv:2310.07777, PRL 132, 201001, 2024; arXiv:2410.22412)
- These are NOT information-theoretic approaches
- They address the "coincidence problem" (why DM ~ 5x baryon) through particle physics

### An Interesting Numerical Coincidence (NOT a prediction)

```
cH_0 / a_0(Milgrom) = 5.46
Omega_DM / Omega_b = 5.41
```

These are suspiciously close but involve **different physics** (ratio of accelerations vs ratio of densities). The near-equality is a known coincidence related to a_0 ~ cH_0/6.

---

## 16. Recommendations

### For the 4-Paper Program

1. **Do NOT pursue JRSWW -> Omega_DM as a Paper 4 result.** It is a dead end.

2. **Paper 3 should focus on the reinterpretation approach** (dark matter = inaccessible information), using:
   - Verlinde's volume-law entanglement mechanism
   - The observer-dependent tau formalism (paper3_observer_dependent_DM.md)
   - Running G / EFT approach (Kumar 2025, Gubitosi et al. 2024)

3. **Paper 4's unification equation** Sigma = D(rho_spacetime || rho_matter) does NOT need to predict Omega_DM from first principles. What it needs to do:
   - Show that the SAME Sigma governs QI recovery (Paper 1), strong gravity (Paper 2), galactic dynamics (Paper 3), and cosmological evolution
   - The specific VALUE of Omega_DM is an initial/boundary condition, not determined by the equation itself (just as F = ma doesn't tell you the mass of the Earth)

4. **The volume-averaged <tau> ~ Omega_m coincidence (1.5%)** could be mentioned as a "curious observation" in a discussion section, clearly labeled as numerology. It is not a result.

### For Future Research

5. **Trivedi-Scherrer (arXiv:2511.10617)** is the most promising existing holographic DM approach. Investigate whether their Ricci cutoff can be derived from tau-framework Sigma:
   - Their L^(-2) ~ (H_dot + 2H^2) ~ Ricci scalar
   - In tau framework: is there a natural channel N whose Sigma reduces to the Ricci scalar?
   - This is a well-defined mathematical question.

6. **Padmanabhan's CosMIn + tau** (paper3_cosmological_tau.md, Section 4.1): the identification CosMIn ~ integral(tau dt) constraining Lambda is more promising than trying to get Omega_DM.

7. **The CMB problem remains the hardest obstacle** for any non-particle DM explanation. Until the tau framework can reproduce the CMB acoustic peak heights without particle DM, Omega_DM h^2 = 0.12 cannot be "explained away" by information theory.

---

## 17. References

### JRSWW and Recovery Bounds
- Junge, Renner, Sutter, Wilde, Winter (2018). Universal Recovery Maps and Approximate Sufficiency of Quantum Relative Entropy. Ann. Henri Poincare 19, 2955-2978. arXiv:1509.07127
- Fawzi, Renner (2015). Quantum Conditional Mutual Information and Approximate Markov Chains. Commun. Math. Phys. 340, 575-611. arXiv:1410.0664

### Cosmological Information
- Padmanabhan, Padmanabhan (2017). Cosmic Information, the Cosmological Constant and the Amplitude of Primordial Perturbations. PLB 773, 81-85. arXiv:1703.06144
- Padmanabhan (2012). Emergence and Expansion of Cosmic Space as due to the Quest for Holographic Equipartition. arXiv:1206.4916

### Holographic Dark Matter
- Trivedi, Scherrer (2025). Dark Matter from Holography. arXiv:2511.10617

### Emergent Gravity and Dark Matter
- Verlinde (2016). Emergent Gravity and the Dark Universe. SciPost Phys. 2, 016 (2017). arXiv:1611.02269
- Ghari, Haghi (2026). arXiv:2601.01715 (Verlinde > MOND at 5.2 sigma)

### Dark Matter Coincidence Problem
- PRL 132, 201001 (2024). Dynamical Explanation of the Dark Matter and Baryon Energy Density Coincidence. arXiv:2310.07777
- arXiv:2410.22412 (2024). Predicting the Dark Matter -- Baryon Abundance Ratio
- arXiv:2506.03269 (2025). Parametric Coincidence in the Baryon to DM Ratio

### de Sitter Horizon and Information
- Gibbons, Hawking (1977). Cosmological Event Horizons, Thermodynamics, and Particle Creation. PRD 15, 2738
- Basso, Maziero, Celeri (2025). Curved-spacetime Crooks relation. PRL 134, 050406. arXiv:2405.03902
- Bianconi (2025). Gravity from Entropy cosmological solutions. arXiv:2510.22545

### Related tau Framework Documents
- paper3_cosmological_tau.md (de Sitter tau, Kodama vector, cosmological expansion as channel)
- paper2_dark_matter_energy.md (comprehensive DM/DE assessment)
- paper3_observer_dependent_DM.md (observer-dependent dark matter reinterpretation)
- paper4_grand_unification.md (master equation Sigma = D(rho_spacetime || rho_matter))
