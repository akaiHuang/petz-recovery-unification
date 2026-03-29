# Sigma = 2 ln Q at the Planck Scale: Implications for Quantum Gravity

**Author**: Sheng-Kai Huang (with systematic analysis)
**Date**: 2026-03-19
**Status**: Research investigation -- five tasks completed
**Purpose**: Investigate the behavior of Sigma = D(rho_spacetime || rho_matter) at the Planck scale and its implications for quantum gravity, the Page curve, the firewall problem, and Planck-scale completion of the framework.

---

## Executive Summary

The Sigma = 2 ln Q framework, established in Papers 1-2 for the gravity sector, has a natural quantum gravity regime where Sigma ~ O(1). This investigation reveals:

1. **Planck-mass black hole**: Sigma = 2 at the Planck radius, giving tau = 1 - e^{-1} = 0.632. The "63% information loss" has a precise operational meaning: a Petz recovery protocol applied to quantum information at the Planck radius recovers at most 37% fidelity. This is the **boundary between the quantum and classical gravity regimes** -- neither fully recoverable nor fully lost.

2. **Known quantum gravity approaches**: Sigma connects to five established programs: Jacobson (1995), Bousso bound, holographic entanglement (JLMS/RT), asymptotic safety, and loop quantum gravity. The deepest connection is to the Dorau-Much (2025) QRE derivation of Einstein equations, where Sigma = 0 gives GR and Sigma > 0 gives GR + Khronon.

3. **Page curve**: The tau(t) trajectory for an evaporating black hole can be described within the framework, but the identification S_Page ~ Sigma is heuristic. A rigorous derivation would require the explicit Hawking channel and its entropy production.

4. **Firewall problem**: The Sigma framework provides a **third option** beyond AMPS's "information loss OR firewall": the exponential metric has no horizon (Sigma remains finite everywhere), so neither premise of the AMPS argument applies. The Petz map mediates between interior and exterior with finite (not zero) fidelity at the throat.

5. **Planck-scale completion**: The exponential metric ds^2 = -e^{-Sigma}dt^2 + ... with Sigma = r_s/r enters a **transition regime** at Sigma ~ O(1) where the semiclassical approximation breaks down. Three candidate completions are identified: noncommutative geometry (Connes spectral triple), loop quantum gravity (polymer quantization of the throat), and asymptotic safety (non-Gaussian fixed point). The strongest connection is to Paper 9 (spectral action), where Sigma = Fisher information of the spectral action is UV-finite in d = 4.

---

## Task 1: Sigma at the Planck Scale

### 1.1 Planck-Mass Black Hole

For a Planck-mass object: M = M_Pl = sqrt(hbar c / G).

The Schwarzschild radius is:
```
r_s = 2GM_Pl/c^2 = 2G sqrt(hbar c / G) / c^2 = 2 sqrt(G hbar / c^3) = 2 l_Pl
```

where l_Pl = sqrt(G hbar / c^3) ~ 1.616 x 10^{-35} m is the Planck length.

The gravitational entropy production at the Planck radius r = l_Pl:
```
Sigma_grav = r_s / r = 2 l_Pl / l_Pl = 2
```

The Petz recovery fidelity:
```
F = exp(-Sigma/2) = exp(-1) = e^{-1} ~ 0.368
```

The temporal asymmetry parameter:
```
tau = 1 - F = 1 - e^{-1} ~ 0.632
```

### 1.2 The Meaning of "63% Information Loss"

This tau = 0.632 has a precise operational interpretation within the framework:

**Operational meaning**: Consider a quantum state rho prepared at the Planck radius r = l_Pl from a Planck-mass object. The gravitational channel N_grav (thermal attenuator with eta = -g_00 = e^{-2}) degrades this state. The optimal Petz recovery protocol R_sigma applied to the degraded state achieves fidelity:

```
F(rho, R_sigma o N_grav(rho)) >= e^{-1} ~ 0.368
```

This means: of all quantum information content at the Planck radius, at most 37% can be faithfully recovered by any protocol. The remaining 63% is irreversibly lost to the gravitational environment.

**Physical regimes defined by Sigma**:

| Sigma | tau | F | Physical regime |
|-------|-----|---|-----------------|
| << 1 | << 1 | ~ 1 | Weak field: nearly complete recoverability, classical gravity |
| ~ 1 | ~ 0.39 | ~ 0.61 | Transitional: quantum gravity onset, partial recoverability |
| = 2 | = 0.63 | = 0.37 | Planck scale: critical regime, significant but not total loss |
| >> 1 | -> 1 | -> 0 | Strong field: near-total loss, classical horizon limit |

### 1.3 Why Sigma = 2 is Special

The value Sigma = 2 at the Planck radius coincides with the wormhole throat value in the exponential metric (Section 2.2 of paper2_wormhole_sigma.md). At the throat:

```
r_throat = r_s/2 = m    (in isotropic coordinates)
Sigma_throat = r_s / r_throat = 2
```

This is NOT a coincidence. For a Planck-mass object, the Schwarzschild radius IS 2 l_Pl, so the throat IS at l_Pl. The Planck scale is where the exponential wormhole throat sits for a Planck-mass black hole.

**Implication**: The Planck scale is the **minimum achievable throat** -- the narrowest wormhole passage permitted by the framework. Objects smaller than the Planck mass have r_s < 2 l_Pl, and for them the semiclassical approximation behind Sigma = r_s/r breaks down before reaching Sigma = 2.

### 1.4 The Sigma = 2 Bound and Minimum Length

A natural conjecture emerges:

**Conjecture (Planck Sigma Bound)**: For any semiclassical geometry, Sigma_max <= 2 at the smallest resolvable scale. Equivalently, the minimum resolvable length is:

```
r_min = r_s / Sigma_max = r_s / 2
```

For a Planck-mass object, r_min = l_Pl. For any mass M, the minimum resolvable radius is r_min = GM/c^2, which is exactly the wormhole throat.

This is structurally similar to -- but distinct from -- the generalized uncertainty principle (GUP):

```
Delta x >= hbar / (2 Delta p) + alpha l_Pl^2 Delta p / hbar
```

In the Sigma framework, the minimum length does not arise from modifying the uncertainty principle, but from the **maximum tolerable entropy production** for a semiclassical description. When Sigma exceeds O(1), the quantum channel N_grav becomes so noisy that the classical spacetime description itself becomes unreliable.

### 1.5 Comparison with Bekenstein-Hawking Entropy

For a Schwarzschild black hole of mass M:

```
S_BH = A / (4 l_Pl^2) = 4 pi r_s^2 / (4 l_Pl^2) = pi r_s^2 / l_Pl^2
```

For a Planck-mass BH: S_BH = pi (2 l_Pl)^2 / l_Pl^2 = 4 pi ~ 12.6.

Meanwhile, Sigma at the horizon (Schwarzschild) diverges, but Sigma at the Planck radius is 2. The relationship:

```
S_BH ~ (Sigma_throat)^2 x (geometric factor) x (M/M_Pl)^2
```

suggests that the extensive Bekenstein-Hawking entropy counts the number of independent "Sigma units" (each contributing O(1) to total entropy) over the horizon area. This is consistent with the brick-wall/entanglement entropy interpretation (t'Hooft 1985, Bombelli-Koul-Lee-Sorkin 1986).

---

## Task 2: Connections to Known Quantum Gravity Approaches

### 2.1 Jacobson (1995): Einstein Equations from Thermodynamics

**Connection strength: STRONG**

Jacobson derived Einstein's equations from the Clausius relation delta Q = T dS applied to local Rindler horizons. The key equation:

```
delta Q = T_Unruh delta S_entanglement
=> T_ab xi^a xi^b = (1/8piG) R_ab xi^a xi^b    (for all xi^a)
=> G_ab + Lambda g_ab = 8piG T_ab
```

In the Sigma framework, Dorau-Much (2025 PRL, arXiv:2510.24491) extended this:

```
Sigma = D(rho_coherent || rho_vacuum) = QRE across bifurcate horizon
```

The Jacobson equation delta Q = T dS becomes, in QRE language:

```
Sigma = 0  =>  Einstein's equations (entanglement equilibrium)
Sigma > 0  =>  Einstein + extra entropy production (our framework)
```

**Key insight**: Jacobson's derivation is the Sigma = 0 limit of our framework. The Khronon field is the source of Sigma > 0 -- it represents **departure from entanglement equilibrium**. At the Planck scale, Sigma ~ O(1), meaning the departure is maximal and the Jacobson equilibrium approximation breaks down.

### 2.2 Bousso Bound and Sigma

**Connection strength: MODERATE**

The covariant Bousso bound states that the entropy flux S through any light sheet L(B) of a codimension-2 surface B of area A satisfies:

```
S[L(B)] <= A / (4G hbar)
```

In Sigma language, consider a light sheet from radius r_1 to r_2 in the exponential metric. The entropy production along the sheet is:

```
Delta Sigma = Sigma(r_2) - Sigma(r_1) = r_s/r_2 - r_s/r_1
```

The Bousso bound constrains the total entropy on the light sheet by the area. For the exponential metric:

```
A(r) = 4 pi r^2 exp(r_s/r) = 4 pi r^2 exp(Sigma)
```

At the throat (Sigma = 2):
```
A_throat = 4 pi (r_s/2)^2 e^2 ~ 4 pi e^2 l_Pl^2 ~ 93 l_Pl^2    (for M = M_Pl)
S_Bousso <= A_throat / (4 l_Pl^2) ~ 23
```

This is consistent with having O(10) Sigma "units" (each O(1)) fitting within the throat area -- a self-consistency check.

**Open question**: Does the Petz bound F >= exp(-Sigma/2) IMPLY the Bousso bound, or are they independent? If the Petz bound is more fundamental, the Bousso bound would follow from information-theoretic principles rather than from the generalized second law.

### 2.3 Holographic Entanglement: JLMS and Sigma

**Connection strength: STRONG (in AdS)**

The JLMS (Jafferis-Lewkowycz-Maldacena-Suh 2016) formula equates boundary QRE with bulk QRE:

```
D(rho_R || sigma_R)_boundary = D(rho_r || sigma_r)_bulk
```

to leading order in G_N. Here rho_R and sigma_R are boundary states on region R, and rho_r, sigma_r are bulk states on the entanglement wedge r of R.

Paper 2 already notes (Section IV.D) that entanglement wedge reconstruction satisfies:

```
F^2 >= exp(-Delta S_gen)
```

where Delta S_gen is the change in generalized entropy (Penington 2020, AEMM 2021). This IS the Petz bound with Sigma = Delta S_gen.

**At the Planck scale**: The JLMS formula receives quantum corrections proportional to G_N. The Sigma framework predicts:

```
Sigma_quantum = Sigma_classical + O(G_N / r^2) = Sigma_classical + O(l_Pl^2 / r^2)
```

When r ~ l_Pl, the quantum corrections become O(1) -- the same order as the classical term. This is the regime where the semiclassical JLMS formula breaks down and a full quantum gravity treatment is needed.

**The quantum extremal surface (QES) prescription** (Engelhardt-Wall 2015) replaces the RT surface with the surface that extremizes generalized entropy S_gen = A/(4G) + S_bulk. In Sigma language:

```
Sigma_QES = Delta A / (4G) + Sigma_bulk
```

At the Planck scale, Delta A ~ l_Pl^2 and Sigma_bulk ~ O(1), so both terms contribute equally. The QES prescription is the natural Planck-scale extension of the Sigma framework.

### 2.4 Asymptotic Safety and Sigma

**Connection strength: MODERATE (speculative)**

If gravity is asymptotically safe (Weinberg 1979, Reuter 1998), there exists a non-Gaussian UV fixed point where the dimensionless gravitational coupling g* = G(k) k^2 / (hbar c) approaches a finite value as k -> infinity.

In the Sigma framework, the running gravitational coupling implies:

```
Sigma(k) = 2 ln Q(k) = r_s(k) / r
```

where r_s(k) = 2 G(k) M / c^2 is the running Schwarzschild radius. At the fixed point:

```
G(k) ~ g* / k^2    for k -> infinity
r_s(k) ~ 2 g* M / (k^2 c^2)  ->  0    as k -> infinity
```

Therefore Sigma(k) -> 0 at the UV fixed point: **gravity becomes information-lossless at trans-Planckian energies**. This is consistent with the asymptotic safety prediction that gravity becomes conformally invariant (and hence entropy-free) at short distances.

The connection to running mu(k) = k (gap B7 in the gap analysis) is suggestive: if both G(k) and mu(k) run, the Sigma at scale k is:

```
Sigma(k) ~ G(k) mu(k)^2 M / c^2 ~ g* M / c^2 = const    (at UV fixed point)
```

The constancy of Sigma at the fixed point would mean the theory has a **finite, nonzero entropy production at all scales** -- a prediction that differs from pure asymptotic safety (where Sigma -> 0).

### 2.5 Loop Quantum Gravity and the Planck-Scale Throat

**Connection strength: STRONG (structural)**

Loop quantum gravity (LQG) predicts that the Schwarzschild singularity is replaced by a bounce at a minimum radius r_min ~ l_Pl, producing a Planck-mass star or a black-to-white-hole transition (Rovelli-Vidotto 2014, Haggard-Rovelli 2015).

This is strikingly similar to the exponential metric's wormhole structure:

| Feature | LQG Planck star | Exponential metric wormhole |
|---------|----------------|---------------------------|
| Singularity? | Resolved (bounce) | Resolved (throat) |
| Minimum radius | r_min ~ l_Pl (polymer effects) | r_throat = r_s/2 (information theory) |
| Sigma at minimum | Not defined in LQG | Sigma = 2 |
| Horizon? | No (effective repulsion) | No (exponential metric) |
| Second asymptotic region? | White hole region | "Other universe" |

Guo and Yuan (2025, arXiv:2506.08821) explicitly computed LQG corrections to the exponential metric ("Papapetrou spacetime"). Their main result: the LQG-corrected exponential metric has a modified throat at:

```
r_throat^{LQG} = r_s/2 (1 + O(l_Pl^2 / r_s^2))
```

For M >> M_Pl, the correction is negligible. For M ~ M_Pl, the correction is O(1), significantly modifying the throat geometry.

**Key result**: LQG and the Sigma framework agree on the qualitative picture (no singularity, finite minimum radius, wormhole-like structure) and differ quantitatively only at the Planck scale. The Sigma framework may be the semiclassical limit of LQG's full quantum description.

---

## Task 3: The Page Curve and Sigma

### 3.1 The Page Curve in tau(t) Language

Paper 2 (Section IV.B) established the heuristic description:

```
tau(t) <= 1 - exp(-S_P(t)/2)
```

where S_P(t) = min(S_rad(t), S_BH(t)) is the Page entropy.

Let us make this more precise. For an evaporating black hole of initial mass M_0:

**Phase I (before Page time, t < t_Page)**:
- S_rad grows as more Hawking quanta are emitted
- S_BH = A(t)/(4G) decreases as the black hole shrinks
- S_P = S_rad (< S_BH): Hawking radiation is nearly thermal
- Sigma ~ S_BH >> 1: extremely high entropy production
- tau ~ 1: initial state is essentially unrecoverable from radiation alone

**At the Page time (t = t_Page)**:
- S_rad = S_BH: the "crossing point"
- S_P = S_rad = S_BH: maximum of the Page entropy
- tau is at its maximum value
- For a Planck-mass BH remnant: S_BH ~ O(1), so tau ~ O(1) at the Page time

**Phase II (after Page time, t > t_Page)**:
- S_P = S_BH (now < S_rad)
- As BH continues evaporating, S_BH -> 0
- Sigma -> 0: entropy production vanishes
- tau -> 0: initial state becomes recoverable from the full radiation

**Complete evaporation (t = t_evap)**:
- S_BH = 0, S_P = 0
- tau = 0: unitarity fully restored

### 3.2 Sigma(t) for an Evaporating Black Hole

The mass decreases via Hawking evaporation:

```
dM/dt = -alpha / M^2    (Stefan-Boltzmann law for BH)
```

where alpha = hbar c^4 / (15360 pi G^2). Solution:

```
M(t) = M_0 (1 - t/t_evap)^{1/3}
t_evap = 5120 pi G^2 M_0^3 / (hbar c^4)
```

The Sigma at the would-be horizon (r = r_s(t) in Schwarzschild, or r = r_s(t)/2 at the throat in exponential):

```
Sigma_throat(t) = 2    (constant in exponential metric -- but r_s(t) shrinks!)
```

The physical content is in the areal radius of the throat:

```
R_throat(t) = (e/2) r_s(t) = e G M(t) / c^2
```

As M(t) -> 0, the throat shrinks to zero size. When R_throat ~ l_Pl:

```
M ~ M_Pl: the Planck scale is reached
Sigma_throat = 2 still, but the semiclassical description breaks down
```

### 3.3 Can We Derive the Page Curve from Sigma?

**Honest assessment: NOT YET.**

The Page curve describes S_rad(t) as a function of time. Deriving it requires:

1. Constructing the Hawking radiation as an explicit quantum channel N_Hawk(t)
2. Computing the entropy production Sigma_Hawk(t) of this channel
3. Showing that tau(t) = 1 - F(rho_0, R o N_Hawk(t)(rho_0)) traces the Page curve

Step 1 is the hardest: the Hawking channel is not a simple thermal attenuator (it involves pair creation at the horizon, entanglement between partners, and backreaction). The recent island formula (Penington 2020, AEMM 2021) provides the entropy via:

```
S_rad = min_{extremal} [A(X) / (4G) + S_bulk(Sigma_rad union Island)]
```

Identifying Sigma with the QRE component of S_gen:

```
Sigma_island = D(rho_rad+island || rho_vacuum) = S_bulk(Sigma_rad union Island)
```

The Page curve would follow if we can show:

```
tau(t) = 1 - exp(-Sigma_island(t)/2)
```

with Sigma_island transitioning from ~ S_BH (before Page time) to ~ 0 (after Page time) as the island "turns on."

**This is a well-defined calculation that could be performed in a solvable model** (e.g., JT gravity + matter, following Almheiri-Engelhardt-Marolf-Maxfield 2020). It is flagged as a concrete project.

### 3.4 The Scrambling Time and Sigma

The scrambling time t_scr ~ (beta/2pi) ln S_BH (Sekino-Susskind 2008) is the time for a perturbation to be "mixed" across all BH degrees of freedom.

In the Sigma framework, the MSS (Maldacena-Shenker-Stanford 2016) bound on scrambling:

```
|d tau/dt| <= (2 pi / beta) tau
```

gives:
```
tau(t) <= tau_0 exp(2 pi t / beta)
```

For initial perturbation tau_0 ~ 1/S_BH (one qubit worth of information):

```
tau ~ 1  at  t ~ (beta/2pi) ln S_BH = t_scr
```

This reproduces the scrambling time from the tau dynamics. **The MSS bound on chaos IS a bound on the rate of tau growth** -- a new interpretation that connects quantum chaos to information recovery.

---

## Task 4: The Firewall Problem

### 4.1 The AMPS Argument

Almheiri-Marolf-Polchinski-Sully (2013) argued that for an old black hole (post-Page time), three assumptions cannot simultaneously hold:

**(A1)** Unitarity: Information is recovered in Hawking radiation
**(A2)** No drama: The horizon is smooth for an infalling observer
**(A3)** EFT validity: Effective field theory holds outside the stretched horizon

The conflict: (A1) requires the late Hawking quantum to be maximally entangled with the early radiation (for unitarity). (A2) requires the same quantum to be maximally entangled with its interior partner (for a smooth horizon). Monogamy of entanglement forbids both.

Conclusion: Either (A1) fails (information loss), or (A2) fails (firewall at horizon), or (A3) fails (new physics).

### 4.2 The Sigma Framework's Resolution: No Horizon, No Paradox

The exponential metric provides a **structural dissolution** of the AMPS argument:

**There is no event horizon.**

In the exponential metric:
- g_00 = -exp(-r_s/r) never vanishes for r > 0
- Sigma remains finite everywhere (Sigma = 2 at the throat)
- tau < 1 everywhere (tau = 0.632 at the throat)

The AMPS argument's premise (A2) assumes a horizon where the vacuum looks smooth. If there is no horizon, (A2) is vacuously satisfied -- there is no surface where one needs to choose between "smooth" and "firewall."

**What replaces the horizon?** The exponential wormhole throat at r = r_s/2 (Sigma = 2). At the throat:

```
F_throat = exp(-1) ~ 0.368
```

Information is degraded (63% loss) but NOT completely destroyed. The Petz map at the throat achieves non-zero fidelity, meaning there IS a recovery protocol -- it just works with 37% efficiency.

### 4.3 Does Q -> infinity Predict a Firewall?

In the Schwarzschild metric: Q = 1/sqrt(-g_00) = 1/sqrt(1 - r_s/r) diverges at r = r_s. This gives Sigma = 2 ln Q -> infinity and F -> 0, which would imply **complete** information loss -- consistent with a firewall.

In the exponential metric: Q = 1/sqrt(-g_00) = exp(r_s/(2r)) diverges only as r -> 0 (the "other side's infinity"), which is infinitely far away in proper distance. At the throat:

```
Q_throat = exp(1) ~ 2.718
Sigma_throat = 2 ln Q_throat = 2
```

Q is large but finite. No firewall.

**The key distinction**: The Sigma framework predicts that **Petz bound saturation** (exponential metric) prevents Q from reaching infinity at any finite proper distance. The firewall arises in Schwarzschild because the bound is NOT saturated -- there is a "gap" between F and exp(-Sigma/2) that allows Sigma -> infinity at finite r.

### 4.4 The Petz Map as Firewall Mediator

In the exponential wormhole:

```
Interior (other side) <--[Petz map at throat]--> Exterior (our side)
```

The Petz recovery map R_sigma at the throat connects the interior and exterior descriptions with fidelity F = e^{-1}. This is neither:
- A perfectly smooth transition (which would require F = 1), nor
- A firewall (which would require F = 0)

It is a **partial recovery** -- a "fuzzy" transition zone where information is degraded but not destroyed. The fidelity e^{-1} is a universal value at the throat (independent of mass, charge, or spin of the object).

### 4.5 Complexity and the Firewall

Paper 2 (Section IV.E) raised the important distinction between fidelity and complexity:

- F > 0 at the throat means recovery IS possible in principle
- The COMPLEXITY of the Petz map may be exponential: C(R_sigma) ~ exp(S_BH)
- A "complexity-weighted tau" would distinguish between:
  - Easy recovery: F ~ 1, C ~ poly(S)
  - Hard recovery: F ~ 1, C ~ exp(S) (the ER=EPR regime)
  - Degraded recovery: F < 1, any C (the throat regime)

**For the firewall problem**: The resolution is not just that F > 0, but that the exponential metric structurally prevents F = 0 at any finite proper distance. The AMPS "firewall" is an artifact of the Schwarzschild geometry's singular behavior at the horizon.

### 4.6 Comparison with Other Resolutions

| Resolution | Mechanism | Sigma language | Status |
|-----------|-----------|----------------|--------|
| Information loss (Hawking) | Unitarity violated | tau = 1 permanent | Excluded (unitarity verified) |
| Complementarity (Susskind) | Observer-dependent descriptions | tau is observer-dependent | Compatible with Sigma |
| ER=EPR (Maldacena-Susskind) | Entanglement = geometry | tau_ER = 0 (wormhole channel) | Compatible; wormhole IS the channel |
| Firewall (AMPS) | Drama at horizon | F = 0 at r = r_s | Only in Schwarzschild, not exponential |
| Fuzzball (Mathur) | No interior; stringy surface | tau ~ 1 at fuzzball surface | Compatible at qualitative level |
| **Exponential metric (this work)** | **No horizon; throat at Sigma=2** | **F = e^{-1} at throat** | **Structural resolution** |
| Island formula (PSSY) | Islands in radiation entropy | Sigma_island -> 0 after Page time | Compatible; deeper derivation needed |

---

## Task 5: Planck-Scale Completion

### 5.1 Where the Semiclassical Description Breaks Down

The exponential metric ds^2 = -e^{-Sigma}dt^2 + e^{+Sigma}(dr^2 + r^2 d Omega^2) is a classical (semiclassical at best) geometry. It breaks down when:

1. **Sigma ~ O(1)**: Quantum fluctuations of the metric are comparable to the classical value. For Sigma = 2 at the throat, delta g / g ~ l_Pl / r_throat. When r_throat ~ l_Pl (i.e., M ~ M_Pl), delta g / g ~ 1 and the metric is no longer well-defined.

2. **Curvature ~ l_Pl^{-2}**: The Kretschner scalar of the exponential metric at the throat is K ~ M^4/r_throat^8 * exp(-4M/r_throat) = M^4 / (M/2)^8 * e^{-8} ~ 256 e^{-8} / M^4. Setting K ~ l_Pl^{-4}: M^4 ~ 256 e^{-8} l_Pl^4, giving M ~ (256 e^{-8})^{1/4} M_Pl ~ 0.7 M_Pl. The semiclassical description breaks down for M < M_Pl, as expected.

3. **The Khronon EFT breaks down**: The strong coupling scale of the ghost condensation is Lambda_3 = (mu^2 M_Pl)^{1/3} ~ 10^{-13} eV (using mu ~ H_0/c). Above this scale, higher-order operators dominate. In energy terms, Lambda_3 corresponds to a length scale much larger than l_Pl, so the Khronon EFT actually breaks down BEFORE the Planck scale. Between Lambda_3^{-1} and l_Pl, the Sigma framework needs UV completion.

### 5.2 Candidate Completion 1: Noncommutative Geometry (Connes)

**Connection to Paper 9 (spectral action): STRONGEST**

In Connes' noncommutative geometry, spacetime at the Planck scale is described by a spectral triple (A, H, D) where:
- A = C^infinity(M) tensor A_F is the algebra (smooth functions tensor finite algebra)
- H is the Hilbert space of spinors
- D is the Dirac operator encoding the metric

The key insight from the spectral action principle:

```
S = Tr(f(D^2 / Lambda^2))
```

where f is a cutoff function at scale Lambda. The spectral action is AUTOMATICALLY UV-finite in d = 4 due to the Seeley-DeWitt expansion:

```
S ~ f_4 Lambda^4 a_0 + f_2 Lambda^2 a_2 + f_0 a_4 + ...
```

If Sigma = Fisher information of the spectral action (Paper 9 conjecture), then Sigma is also UV-finite:

```
Sigma = delta^2 S_spectral / delta g^2 = UV-finite quantity
```

**At the Planck scale**: The spectral triple is modified by the finite algebra A_F. The Dirac operator on A_F has discrete eigenvalues, and the "minimum distance" on A_F is:

```
d_NCG(p, q) = sup { |a(p) - a(q)| : ||[D, a]|| <= 1 }
```

For the Standard Model spectral triple, this gives a minimum distance ~ 1/M_H (the Higgs mass), not l_Pl. The Planck length enters through the gravitational part of D, not the finite part.

**Prediction**: In the NCG framework, Sigma at the Planck scale is replaced by:

```
Sigma_NCG = D(rho_1 || rho_2) on (A tensor A_F)
```

where rho_1 and rho_2 are states on the full spectral triple. The finite algebra A_F provides automatic regularization, and Sigma_NCG remains finite even at r = 0 (where the classical Sigma diverges).

### 5.3 Candidate Completion 2: Loop Quantum Gravity

**Connection: STRUCTURAL**

In LQG, the area operator has discrete spectrum:

```
A = 8 pi gamma l_Pl^2 sum_i sqrt(j_i(j_i + 1))
```

where gamma is the Barbero-Immirzi parameter and j_i are spin labels.

The minimum nonzero area is:

```
A_min = 4 pi sqrt(3) gamma l_Pl^2 ~ 20 l_Pl^2    (for gamma ~ 0.274)
```

**In the Sigma framework**: The Bekenstein-Hawking entropy S = A/(4 l_Pl^2) has a minimum value:

```
S_min = A_min / (4 l_Pl^2) ~ 5
```

This means the minimum Sigma for a black hole throat:

```
Sigma_min ~ O(1)    (consistent with Sigma_throat = 2)
```

The LQG area gap provides a natural explanation for why Sigma_throat = 2 is the minimum: the throat cannot shrink below the area gap, and Sigma ~ O(1) at the gap.

Guo-Yuan (2025) computed the LQG-corrected exponential metric explicitly, finding:

```
g_00^{LQG} = -exp(-r_s/r) * [1 + O(l_Pl^2 / r^2)]
Sigma^{LQG} = Sigma_classical + O(l_Pl^2 / r^2)
```

At the throat: Sigma^{LQG} = 2 + O(1) for M ~ M_Pl, significantly different from the classical value.

### 5.4 Candidate Completion 3: Asymptotic Safety

**Connection: MODERATE**

In the asymptotic safety scenario, the running gravitational coupling G(k) approaches a non-Gaussian fixed point g* at high energies. The running Sigma:

```
Sigma(k) = 2 G(k) M / (c^2 r)
```

At the fixed point (k -> infinity):

```
G(k) -> g* hbar c / k^2
Sigma(k) -> 2 g* hbar M / (c k^2 r)  ->  0
```

Sigma vanishes at trans-Planckian energies, meaning gravity becomes information-lossless. This is the "conformal phase" of gravity, where the theory is scale-invariant.

**Implication for the Planck scale**: At k ~ 1/l_Pl:

```
Sigma(k_Pl) = 2 g* M / (M_Pl r/l_Pl) = 2 g* (M/M_Pl) (l_Pl/r)
```

For M = M_Pl, r = l_Pl: Sigma ~ 2g*. If g* ~ O(1) (as typical in asymptotic safety), this gives Sigma ~ O(1) at the Planck scale -- consistent with our earlier estimate.

### 5.5 Synthesis: The Three Completions Compared

| Feature | NCG (Connes) | LQG | Asymptotic Safety |
|---------|-------------|-----|-------------------|
| UV behavior of Sigma | Finite (spectral regularization) | Discrete (area gap) | -> 0 (conformal phase) |
| Minimum Sigma | O(1) from A_F | O(1) from area gap | 0 (approaches from above) |
| Minimum length | Dirac distance on A_F | sqrt(A_min) ~ l_Pl | None (but effective) |
| Connection to SM | DIRECT (A_F = C+H+M_3(C)) | Indirect (through coupling to matter) | Indirect (matter as perturbations) |
| Connection to Paper 9 | STRONGEST (spectral action = entropy) | Moderate (BH entropy) | Moderate (fixed point structure) |
| Status in framework | Highest priority (Paper 9) | Structural agreement | Running G connection |

**The most promising completion is NCG**, because:
1. It directly connects to the Standard Model (A_F determines the gauge group)
2. The spectral action provides a UV-finite Sigma
3. The Fisher information interpretation bridges QRE and field equations
4. It is the subject of Paper 9 (highest priority in the publication plan)

### 5.6 Minimum Length from Sigma

Combining all three approaches, a consistent picture emerges:

```
r_min ~ l_Pl (from Sigma ~ O(1) at the transition)
```

but the mechanism differs:
- NCG: spectral distance on the noncommutative algebra sets a minimum distance
- LQG: area gap prevents the throat from shrinking below ~ l_Pl
- AS: running G -> 0 at high k prevents Sigma from growing at short distances

In the Sigma framework, the minimum length can be expressed as:

```
r_min = r_s / Sigma_max
```

If Sigma_max = 2 (the throat value), then r_min = r_s/2. For M = M_Pl, this is l_Pl. For M > M_Pl, this is the classical throat. For M < M_Pl, the semiclassical framework breaks down and one of the three completions takes over.

**This is a falsifiable prediction**: the minimum resolvable scale is set by information-theoretic saturation (Sigma = 2), not by a fundamental length scale postulated by hand.

---

## Open Questions and Future Directions

### Immediate (next 2-4 weeks)

1. **Page curve from JT gravity + Sigma**: Compute tau(t) explicitly in the JT gravity model (solvable), where the island formula gives the Page curve analytically. Verify that tau(t) = 1 - exp(-S_P(t)/2) is satisfied.

2. **MSS bound as tau bound**: Formalize the scrambling interpretation: |d tau/dt| <= (2pi/beta) tau. Check whether this follows from the Petz bound alone or requires additional input (e.g., chaos assumption).

### Short term (1-3 months)

3. **QES in the exponential metric**: Compute the quantum extremal surface for the exponential metric in a solvable model. Does the QES coincide with the throat (r = r_s/2)?

4. **LQG corrections to Sigma**: Use the Guo-Yuan (2025) results to compute Sigma^{LQG}(r) near the throat. How does Sigma_max change?

### Medium term (3-6 months, connects to Paper 9)

5. **Sigma as Fisher information of spectral action**: The toy model S^1 x M_2(C). Compute delta^2 S_spectral / delta g^2 and compare with Sigma = 2 ln Q. This is the go/no-go test for the NCG completion.

6. **Area-entropy connection**: Derive S_BH = A/(4G) from Sigma = 2 at the throat, using the extensivity argument (Sigma_total >= N * sigma_min with N ~ A/l_Pl^2 independent modes).

### Long term (6-18 months)

7. **Full quantum gravity regime**: Extend the Sigma framework to the Planck scale using the NCG (or LQG, or AS) completion. The goal: a finite, well-defined Sigma at all scales, reproducing both the semiclassical limit (Papers 1-4) and the full quantum gravity regime.

---

## Summary of Key Results

| Task | Finding | Status | Confidence |
|------|---------|--------|------------|
| 1. Sigma at Planck scale | Sigma = 2, tau = 0.632, F = 0.368 at Planck radius | Calculation | HIGH |
| 1b. Minimum length | r_min = r_s/Sigma_max = r_s/2 from information saturation | Conjecture | MODERATE |
| 2. Known QG approaches | Strong connections to Jacobson, JLMS, LQG; moderate to AS, Bousso | Analysis | HIGH |
| 3. Page curve | tau(t) traces Page curve qualitatively; rigorous derivation needs explicit Hawking channel | Heuristic | MODERATE |
| 3b. Scrambling time | MSS bound = tau growth bound; t_scr from tau dynamics | New interpretation | MODERATE |
| 4. Firewall problem | Exponential metric dissolves AMPS: no horizon => no paradox; F = e^{-1} at throat | Structural resolution | HIGH |
| 5. Planck-scale completion | Three candidates: NCG (strongest, Paper 9), LQG (structural match), AS (running G) | Analysis | MODERATE |
| 5b. NCG connection | Sigma = Fisher info of spectral action; UV-finite in d=4; go/no-go test via toy model | Conjecture (Paper 9) | LOW-MODERATE |

---

## References

### Primary (used directly in this analysis)

1. Jacobson (1995) PRL 75, 1260 -- Thermodynamics of spacetime
2. Dorau-Much (2025) PRL -- QRE to Einstein equations
3. Chen-Penington-Salton (2020) JHEP -- Entanglement wedge via Petz map
4. Almheiri-Marolf-Polchinski-Sully (2013) JHEP -- Firewalls
5. Penington (2020) JHEP -- Entanglement wedge and information problem
6. Almheiri-Engelhardt-Marolf-Maxfield (2021) RMP -- Entropy of Hawking radiation
7. Maldacena-Shenker-Stanford (2016) JHEP -- Bound on chaos
8. CCSvS (2018) arXiv:1809.02944 -- Entropy = spectral action
9. Lashkari-Van Raamsdonk (2015) arXiv:1508.00897 -- Fisher info = canonical energy
10. Guo-Yuan (2025) arXiv:2506.08821 -- LQG corrections to Papapetrou metric
11. Boonserm-Ngampitipan-Simpson-Visser (2018) PRD 98, 084048 -- Exponential metric = traversable wormhole

### Supporting (context and connections)

12. Engelhardt-Wall (2015) -- Quantum extremal surfaces
13. Jafferis-Lewkowycz-Maldacena-Suh (2016) -- JLMS formula
14. Maldacena-Susskind (2013) -- ER=EPR
15. Cotler et al. (2019) PRX -- HP protocol = Petz map
16. Rovelli-Vidotto (2014) -- Planck stars
17. Haggard-Rovelli (2015) -- Black hole to white hole transition
18. Reuter (1998) -- Asymptotic safety
19. Sekino-Susskind (2008) -- Fast scramblers
20. Dong-Khalkhali-van Suijlekom (2019) arXiv:1903.09624 -- Chemical potential and spectral action

---

*Last updated: 2026-03-19*
*This analysis addresses gap E1 (Sigma at Planck scale) and E2 (Black hole information problem) in the comprehensive gap analysis.*
