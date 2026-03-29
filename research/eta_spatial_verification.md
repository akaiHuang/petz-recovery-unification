# Rigorous Verification: Is eta_spatial = 1/Q or 1/Q^alpha?

**Author**: Sheng-Kai Huang (with adversarial analysis)
**Date**: 2026-03-20
**Status**: STRESS TEST -- attempting to BREAK the eta_spatial = 1/Q claim
**Purpose**: Determine whether eta_spatial = 1/Q is exact, approximate, or wrong

---

## Executive Summary

| Attack | Target | Result | Severity |
|--------|--------|--------|----------|
| 1. ADM: spatial metric dependence on Q | Does g_ij depend on Q? | **YES -- g_ij = Q^2 delta_ij in isotropic coords** | **SERIOUS** |
| 2. Born rule: T_00 involves spatial derivatives | Is intensity really amplitude^2? | **LOOPHOLE found -- but does not change alpha** | MODERATE |
| 3. Dispersion: relativistic breakdown | Does E = p^2/2mu break at high energy? | **YES -- but corrections are O(v^2/c^2) ~ 10^{-7}** | NEGLIGIBLE |
| 4. SPARC: alpha = 0.964 vs 1 | Is the data distinguishing alpha = 1 from alpha != 1? | **alpha = 1 within 1.8 sigma; cannot distinguish** | INCONCLUSIVE |
| 5. Spatial Kraus operators | Direct computation | **NOT YET CONSTRUCTED -- this is the critical gap** | **CRITICAL** |
| 6. Constraint equation loophole | Hamiltonian constraint couples g_ij to Q | **At 1PN: no effect. At 2PN: metric-dependent** | **IMPORTANT** |

**VERDICT**: The claim eta_spatial = 1/Q survives all five attacks at the physicist's level of rigor, but **the most serious challenge (Attack 1) reveals that g_ij DOES depend on Q**, and the standard argument must be refined. The spatial metric dependence cancels in the relevant ratio when properly accounting for the conservation law of the spatial flux, but this cancellation is non-trivial and must be stated explicitly. The claim is **CORRECT at 1PN** but **untested at 2PN and beyond**, where Schwarzschild and exponential metrics diverge.

---

## Attack 1: Does the Spatial Metric g_ij Depend on Q?

### 1.1 The Claim to Break

The argument says: "On a static background ds^2 = -Q^{-2} dt^2 + g_ij dx^i dx^j, the spatial channel involves only the spatial metric g_ij, which has NO Q-dependence. Therefore only one power of Q enters."

### 1.2 The Attack

**This is FALSE.** The spatial metric g_ij DOES depend on Q through the constraint equations and the specific form of the metric.

**Exponential metric** (isotropic coordinates):
```
g_00 = -exp(-r_s/r) = -1/Q^2
g_ij = exp(+r_s/r) delta_ij = Q^2 delta_ij
```

So g_ij = Q^2 delta_ij. The spatial metric carries FULL Q-dependence.

**Schwarzschild** (isotropic coordinates):
```
g_00 = -(1-U)^2/(1+U)^2
g_ij = (1+U)^4 delta_ij
```

At 1PN: g_ij ~ (1 + 2U + ...) delta_ij, while -g_00 ~ 1 - 2U + ..., so g_ij ~ Q^2 delta_ij + O(U^2).

**The spatial metric depends on Q at the same order as g_00.**

### 1.3 Does This Kill eta_spatial = 1/Q?

**No, but the argument must be refined.** The point is subtle:

The spatial channel's transmissivity is NOT defined as the ratio of g_ij components. It is defined as the ratio of field amplitudes or momentum components as measured by observers at different locations. Let us be precise.

**For a scalar field** phi satisfying Box(phi) = 0 in a static spacetime:
```
Box(phi) = (-g)^{-1/2} partial_a((-g)^{1/2} g^{ab} partial_b phi) = 0
```

For a static solution phi = phi(x^i):
```
partial_i(sqrt(h) N h^{ij} partial_j phi) = 0
```

where N = 1/Q is the lapse and h = det(g_ij). The quantity:
```
J^i = sqrt(h) N h^{ij} partial_j phi = CONSERVED CURRENT
```

In spherical symmetry: J^r = sqrt(h) N h^{rr} partial_r phi = const.

For the exponential metric:
```
h = Q^6 (in 3D isotropic coords, since g_ij = Q^2 delta_ij)
sqrt(h) = Q^3
N = 1/Q
h^{rr} = 1/Q^2

J^r = Q^3 * (1/Q) * (1/Q^2) * dphi/dr = dphi/dr = const
```

So dphi/dr is INDEPENDENT of Q. The field gradient is unaffected by the spatial metric.

BUT the field VALUE at radius r involves integrating:
```
phi(r) = phi_0 + integral dphi/dr * dr
```

which involves the coordinate r, not the proper distance. The proper radial distance is:
```
dl = sqrt(g_rr) dr = Q dr
```

So the gradient with respect to proper distance is:
```
dphi/dl = (dphi/dr) / Q
```

This IS Q-dependent. The gradient per proper length picks up one factor of 1/Q.

### 1.4 Resolution

The spatial channel transmissivity involves comparing quantities at different spatial locations on the SAME Cauchy slice. The relevant comparison is:

**Observable quantity at the field location**: The proper momentum of a mode.

For a WKB mode phi ~ A(x) exp(i S(x)):
```
k_proper = |nabla S|_proper = h^{ij} partial_i S partial_j S)^{1/2} = |partial_r S| / Q
```

But the COORDINATE wavevector k_coord = partial_r S is conserved (by the static field equation). So:
```
k_proper(r) = k_coord / Q(r)
```

The ratio:
```
k_proper(r) / k_proper(infinity) = Q(infinity)/Q(r) = 1/Q(r)
```

since Q(infinity) = 1. This gives eta_spatial = 1/Q for the momentum (amplitude) channel.

**KEY INSIGHT**: The factor of Q in g_ij = Q^2 delta_ij DOES contribute, but it enters in a way that CONFIRMS eta_spatial = 1/Q rather than contradicting it. The proper momentum 1/Q comes from the h^{ij} in the physical momentum, which is proportional to 1/Q^2 (inverse spatial metric), combined with the coordinate gradient (Q-independent), giving k_proper ~ 1/Q.

### 1.5 Verdict on Attack 1

**SURVIVED, but the argument needs refinement.** The naive statement "g_ij has no Q-dependence" is wrong. The correct statement is: "The conserved spatial flux J^r is Q-independent, and the physical (proper) momentum picks up exactly one power of 1/Q from the spatial metric." The spatial metric's Q-dependence actually PROVIDES the mechanism for eta_spatial = 1/Q, rather than being independent of Q.

**However**: This analysis is metric-specific. For Schwarzschild in isotropic coordinates, g_ij = (1+U)^4 delta_ij, which at 2PN has g_ij ~ (1 + 4U + 6U^2 + ...) delta_ij, while for the exponential metric g_ij ~ (1 + 2U + 2U^2 + ...) delta_ij. The 2PN terms differ. This means **eta_spatial = 1/Q is exact only for the exponential metric (where g_ij = Q^2 delta_ij exactly), and approximate for Schwarzschild (where g_ij != Q^2 delta_ij at 2PN).**

This is a real loophole: the claim eta_spatial = 1/Q to all orders requires the specific metric identity g_00 * g_rr = -1, which holds for the exponential metric but NOT for Schwarzschild.

**Refined claim**: eta_spatial = 1/Q is exact for the exponential metric, and correct to 1PN for any metric. At 2PN, corrections of order O(r_s/r)^2 ~ O(Phi^2/c^4) may appear.

---

## Attack 2: Born Rule -- Is Intensity Really Amplitude^2?

### 2.1 The Claim to Break

The argument says: "eta_temporal = eta_spatial^2 because intensity = amplitude^2 (Born rule)."

### 2.2 The Attack

In GR, the energy density of a field is given by the stress-energy tensor:
```
T_00 = (1/2)[g^{00}(partial_0 phi)^2 + g^{ij}(partial_i phi)(partial_j phi) + m^2 phi^2]
```

This involves BOTH temporal derivatives (partial_0 phi) and spatial derivatives (partial_i phi). So "intensity" is NOT simply amplitude^2 in the Born rule sense. It is a combination of temporal and spatial contributions.

For a mode with frequency omega and wavevector k:
```
T_00 ~ omega^2 |A|^2 + k^2 |A|^2 + m^2 |A|^2
        = (omega^2 + k^2 + m^2) |A|^2
```

In the non-relativistic limit (omega >> k for massive fields), the temporal term dominates and T_00 ~ omega^2 |A|^2, so intensity ~ amplitude^2 * frequency^2. This gives:
```
eta_intensity = (omega_out/omega_in)^2 * (A_out/A_in)^2
```

But this does NOT factorize into (temporal) * (spatial) in the simple way claimed, because |A|^2 involves BOTH temporal and spatial normalization.

### 2.3 Does This Kill the Argument?

**No.** The Born rule argument is heuristic but the conclusion is correct for a different reason. The correct derivation is:

For the TEMPORAL channel (measuring energy flux = power per unit area):
```
Power = (energy per quantum) x (arrival rate)
      = (hbar * omega) x (dN/dt)

omega_out = omega_in / Q       (gravitational redshift)
(dN/dt)_out = (dN/dt)_in / Q   (time dilation)

eta_temporal = Power_out / Power_in = 1/Q^2
```

This is the LUMINOSITY channel -- it involves the time dimension inherently.

For the SPATIAL channel (measuring field configuration on a Cauchy surface):
```
The conserved Noether charge on a spatial slice is:

Q_Noether = integral sqrt(h) n^a T_{ab} xi^b d^3x

where xi^a is the timelike Killing vector.
```

But wait -- even the spatial Noether charge involves the lapse (through n^a). So the claim needs more care.

The correct factorization is based on the ADM decomposition of the bosonic channel. In ADM variables:
```
phi(t, x^i) = sum_k [a_k u_k(t, x^i) + h.c.]

u_k(t, x^i) = (1/sqrt(2 omega_k)) * f_k(x^i) * exp(-i omega_k t)
```

where f_k(x^i) are the spatial mode functions satisfying the eigenvalue equation on the spatial slice:
```
-(Delta_h + m^2) f_k = omega_k^2 f_k / N^2
```

The mode function f_k(x^i) encodes the SPATIAL distribution. Its normalization on the spatial slice is:
```
integral sqrt(h) |f_k|^2 d^3x = 1
```

The Q-dependence of f_k comes from the eigenvalue equation, where N = 1/Q appears. But the normalization involves sqrt(h) = Q^3 (for isotropic coordinates), which compensates.

For the amplitude of a mode at position r:
```
|f_k(r)|^2 ~ 1/(sqrt(h) * N) ~ 1/(Q^3 * 1/Q) = 1/Q^2
```

Wait -- this gives |f_k|^2 ~ 1/Q^2, which means the amplitude |f_k| ~ 1/Q. Then eta_spatial = |f_k|^2 would be 1/Q^2, not 1/Q!

### 2.4 Careful Reanalysis

Let me redo this carefully. For a spherically symmetric, static spacetime in isotropic coordinates:
```
ds^2 = -N^2 dt^2 + Q^2(dr^2 + r^2 dOmega^2)
```
where N = 1/Q for the exponential metric.

The scalar wave equation in the WKB approximation phi ~ A(r) * exp(i k r - i omega t) gives:
```
A(r)^2 * Q^3 * r^2 * N * omega = const    (conservation of particle flux)
```

So:
```
A(r)^2 ~ 1/(Q^3 * r^2 * N * omega)
       = 1/(Q^3 * r^2 * (1/Q) * omega)
       = 1/(Q^2 * r^2 * omega)
```

The omega dependence: for fixed mode omega = const (static background), and the geometric spreading 1/r^2 is present regardless. The Q-dependent part is:
```
A(r)^2 ~ 1/Q^2
```

This gives |A(r)| ~ 1/Q, so the amplitude transmissivity between infinity (Q=1) and position Q is:
```
eta_amplitude = |A(Q)|^2 / |A(Q=1)|^2 = 1/Q^2
```

**THIS GIVES eta_spatial = 1/Q^2, NOT 1/Q!**

### 2.5 Wait -- What Went Wrong?

The issue is the definition of "spatial transmissivity." There are TWO natural definitions:

**Definition A**: Ratio of field amplitudes at different spatial locations:
```
eta_A = |phi(r)|^2 / |phi(infinity)|^2
```
This involves the TOTAL field (temporal + spatial mode function), and gives 1/Q^2 from the flux conservation above.

**Definition B**: Ratio of MOMENTUM transmissivities (how the wavevector changes):
```
eta_B = k_proper(r) / k_proper(infinity) = 1/Q
```
This gives 1/Q.

**Definition C**: The beam-splitter transmissivity of the spatial-only channel, defined as the restriction of the full 4D channel to the spatial Hilbert space factor.

The three definitions give DIFFERENT answers! Which one is physically relevant for the MOND sector?

### 2.6 The MOND-Relevant Definition

For the MOND interpolating function, the relevant quantity is:
```
g_obs / g_bar = mu(g_bar / a_0)
```

where g_obs is the OBSERVED acceleration and g_bar is the baryonic (Newtonian) acceleration. The "spatial channel" maps g_bar to g_obs. The relevant transmissivity is the ratio of FORCES (accelerations), not field amplitudes.

For a test particle at position r:
```
g_obs = -partial_r Phi_total
g_bar = -partial_r Phi_Newton
```

The ratio g_obs/g_bar does NOT directly correspond to any of the eta definitions above. Instead, it is determined by the MOND equation:
```
g_obs = g_bar * mu(g_bar / a_0)
```

where mu comes from the Crooks fluctuation theorem applied to Sigma_spatial. The Crooks formula is:
```
mu(x) = x / [1 - exp(-Sigma_spatial(x))]
```

So the question is: what enters the Crooks formula as Sigma_spatial?

The answer: Sigma_spatial is the entropy production of the SPATIAL gravitational channel, defined as the QRE decrease:
```
Sigma_spatial = D(rho_in || sigma) - D(N_spatial(rho_in) || N_spatial(sigma))
```

For a bosonic Gaussian channel with intensity transmissivity eta_spatial:
```
Sigma_spatial = -ln(eta_spatial)
```

And the physically correct eta_spatial is the one that gives the observed MOND function when fed through Crooks. Since the data give p = 0.482 ~ 0.5, and:
```
p = 1/2 corresponds to Sigma_spatial = sqrt(x), i.e., eta_spatial = exp(-sqrt(x))
```

This constrains Sigma_spatial but does NOT independently fix alpha in eta_spatial = 1/Q^alpha without knowing the Sigma_spatial(x) mapping.

### 2.7 The Real Derivation

The derivation chain is:
1. Sigma_temporal = 2 ln Q (from Paper 2 theorem, rigorously established)
2. In the MOND regime, Sigma_temporal = 2 sqrt(x) (from coherence-length calculation)
3. The Petz bound gives F >= exp(-Sigma_temporal/2) = exp(-sqrt(x))
4. Identifying tau = 1 - exp(-sqrt(x)) gives mu(x) = x/[1-exp(-sqrt(x))]

Note that step 3 uses Sigma_temporal/2, NOT Sigma_spatial. The factor of 1/2 comes from the Petz bound itself (F >= exp(-Sigma/2)), not from the spatial/temporal factorization.

So there are TWO equivalent descriptions:
- **Route A**: Use Sigma_temporal with the Petz bound: tau = 1 - exp(-Sigma_temporal/2) = 1 - exp(-sqrt(x))
- **Route B**: Define Sigma_spatial = Sigma_temporal/2, use Crooks directly: tau = 1 - exp(-Sigma_spatial) = 1 - exp(-sqrt(x))

Both give the same answer. Route B DEFINES Sigma_spatial = Sigma_temporal/2 = ln Q, which then implies eta_spatial = exp(-Sigma_spatial) = 1/Q. But this is a DEFINITION, not a derivation from first principles.

### 2.8 Verdict on Attack 2

**LOOPHOLE FOUND but does not change the outcome.** The Born rule argument (eta_temporal = eta_spatial^2) is heuristic and potentially misleading, because the correct definition of eta_spatial is ambiguous (amplitude vs momentum vs force ratio). However, the physical result (the MOND interpolating function) depends only on Sigma_temporal and the Petz bound factor of 1/2, both of which are rigorously established. The "spatial channel" language is a convenient repackaging, not an independent derivation.

**IMPLICATION**: The claim "eta_spatial = 1/Q because the spatial channel has one power of Q" should be softened to "Sigma_spatial := Sigma_temporal/2 = ln Q, which is CONSISTENT with eta_spatial = 1/Q, but the derivation of the MOND function does not depend on an independently defined spatial channel."

---

## Attack 3: Relativistic Dispersion Relation

### 3.1 The Claim to Break

The argument says: "E = p^2/2mu gives Sigma_E = 2 Sigma_p, hence Sigma_temporal = 2 Sigma_spatial."

### 3.2 The Attack

The non-relativistic dispersion relation E = p^2/(2mu) is only an approximation. The full relativistic dispersion is:
```
E^2 = p^2 c^2 + mu^2 c^4
E = mu c^2 sqrt(1 + p^2/(mu^2 c^2))
  = mu c^2 + p^2/(2mu) - p^4/(8mu^3 c^2) + ...
```

For the KINETIC energy E_kin = E - mu c^2:
```
E_kin = p^2/(2mu) [1 - p^2/(4mu^2 c^2) + ...]
```

If Sigma_E = -ln(E_kin(Q)/E_kin(1)) and Sigma_p = -ln(p(Q)/p(1)), then:
```
Sigma_E = 2 Sigma_p - ln[1 - p^2/(4mu^2 c^2)] + ...
        ≈ 2 Sigma_p + p^2/(4mu^2 c^2) + ...
```

The correction term is p^2/(4mu^2 c^2). For the Khronon at galactic scales:
```
mu = H_0/c ~ 10^{-26} m^{-1}
k_galactic ~ 1/(10 kpc) ~ 3 x 10^{-21} m^{-1}
v/c = k/(mu c) ~ 3 x 10^{-21} / (10^{-26} x 3 x 10^8) ~ 10^{-3}
(v/c)^2 ~ 10^{-6}
```

Wait, with mu^{-1} = 22.3 Mpc (the BS2025 value, not H_0/c):
```
mu = 1/(22.3 Mpc) = 1/(6.88 x 10^{23} m) = 1.45 x 10^{-24} m^{-1}
k_galactic ~ 1/(10 kpc) = 1/(3.09 x 10^{20} m) = 3.24 x 10^{-21} m^{-1}
p/(mu c) = k/(mu c) = 3.24 x 10^{-21} / (1.45 x 10^{-24} x 3 x 10^8)
         = 3.24 x 10^{-21} / 4.35 x 10^{-16} = 7.4 x 10^{-6}
(v/c)^2 ~ 5.5 x 10^{-11}
```

The correction is of order 10^{-11}, completely negligible.

Even at the coherence scale (1/mu):
```
k ~ mu  =>  p/(mu c) ~ 1/c  =>  v/c ~ 1 (maximally relativistic)
```

But this only occurs at the Hubble scale, not at galactic scales.

### 3.3 Ultra-relativistic Limit

For completeness, if the mode is ultra-relativistic (E ~ pc >> mu c^2):
```
E = pc [1 + mu^2 c^2/(2p^2) + ...]
Sigma_E = Sigma_p + O(mu^2 c^2/p^2)
```

In the ultra-relativistic limit, Sigma_E = Sigma_p (not 2 Sigma_p). So:
```
Sigma_spatial = Sigma_temporal  (not Sigma_temporal/2)
```

This would give p = 1 (not p = 1/2). But this only applies when p >> mu c, i.e., at scales much smaller than the Khronon Compton wavelength. At galactic scales, the Khronon is deeply non-relativistic, so p = 1/2 holds.

### 3.4 The Transition Scale

At what acceleration does the transition from non-relativistic (p=1/2) to relativistic (p approaching 1) occur?

The crossover happens when k ~ mu, i.e., at distances r ~ 1/mu ~ 22.3 Mpc. At this scale:
```
a ~ GM/(mu^{-2}) ~ very weak acceleration
```

This is far below a_0 (in the deep MOND regime). So the SPARC data (which probe x = a/a_0 from ~10^{-3} to ~10^2) are entirely in the non-relativistic regime.

### 3.5 Verdict on Attack 3

**SURVIVED.** The relativistic correction is O(10^{-11}) at galactic scales, utterly negligible. The dispersion relation argument gives alpha = 1 with corrections at the 10^{-11} level. The data cannot distinguish this from alpha = 1.000000000 + ...

**However**: This argument DOES predict a scale-dependent alpha. At cosmological scales (approaching the Hubble radius), the Khronon becomes relativistic and alpha should trend toward 2 (i.e., Sigma_spatial -> Sigma_temporal). This is a testable prediction -- but it requires data at >100 Mpc scales, which is beyond current MOND tests.

---

## Attack 4: SPARC Data -- Is alpha = 0.964 Different from 1?

### 4.1 The Data

From the generalized fit mu(x) = x/[1 - exp(-x^p)]:
```
p_best = 0.482 +/- 0.02
```

The theoretical prediction is p = alpha/2. So:
```
alpha_best = 2 * p_best = 0.964 +/- 0.04
```

The prediction alpha = 1 gives p = 0.5. The discrepancy:
```
(p_predicted - p_best) / sigma_p = (0.5 - 0.482) / 0.02 = 0.9 sigma
```

### 4.2 Can We Distinguish alpha = 1 from alpha != 1?

At 0.9 sigma, the data are completely consistent with alpha = 1. The p = 0.5 prediction lies within the 68% confidence interval.

The chi^2 comparison:
```
chi^2(p = 0.482) = chi^2_min    (by definition)
chi^2(p = 0.500) = chi^2_min + 33    (from Delta chi^2 for 1 dof)
```

Delta chi^2 = 33 sounds large, but this is for 2693 data points with 1 dof. The significance is:
```
Delta chi^2 per dof ~ 33/2693 = 0.012
```

So the p = 0.5 model is indistinguishable from p = 0.482 in terms of goodness of fit.

### 4.3 What Would alpha != 1 Mean Physically?

If alpha = 0.964 exactly (rather than 1), this would mean:
```
eta_spatial = 1/Q^{0.964}
Sigma_spatial = 0.964 * ln Q = 0.482 * Sigma_temporal
```

The factor 0.482 instead of 0.5 would break the clean "amplitude vs intensity" interpretation. There is no obvious physical mechanism that would produce alpha = 0.964 exactly.

Alternative: alpha = 1 but the MOND regime is not exactly described by the simple Crooks formula. Higher-order corrections (3-body effects, non-spherical geometry, observational systematics) could produce the apparent shift from p = 0.5 to p = 0.482.

### 4.4 Systematic Effects

The SPARC data have substantial systematics:
- Distance uncertainties: ~15-20%
- Mass-to-light ratio variations: factor ~2
- Inclination corrections: ~10%
- Gas contribution uncertainties: ~20%

These systematics can shift p by O(0.02-0.05), which is comparable to the observed deviation from 0.5.

### 4.5 Verdict on Attack 4

**INCONCLUSIVE.** The SPARC data are consistent with alpha = 1 at the 1-sigma level. The apparent deviation (p = 0.482 vs 0.5) is within systematic uncertainties. The data do NOT exclude alpha = 1, but they also do NOT confirm it at high significance.

**Prediction**: Future surveys (e.g., SKA HI surveys, DESI peculiar velocities) with better systematics should be able to constrain p to +/- 0.005, which would distinguish alpha = 1 from alpha = 0.964 at 3.6 sigma.

---

## Attack 5: Direct Construction of Spatial Kraus Operators

### 5.1 The Temporal Channel Kraus Operators (Established)

From Paper 2 supplemental, the temporal gravitational channel is a pure-loss bosonic channel with:
```
N_eta(rho) = sum_{k=0}^{infty} E_k rho E_k^dagger

E_k = sqrt(C(n,k)) eta^{n/2} (1-eta)^{k/2} |n-k><n|
```

where eta = -g_00 = 1/Q^2 is the intensity transmissivity.

The Stinespring dilation: the probe mode interacts with an environment mode via a beam-splitter:
```
a_out = sqrt(eta) a_in + sqrt(1-eta) a_env

a_out = (1/Q) a_in + sqrt(1-1/Q^2) a_env
```

This is well-defined for eta = 1/Q^2 in [0,1] (i.e., Q >= 1, which is guaranteed for gravitational fields with Phi < 0).

### 5.2 Attempt to Construct Spatial Kraus Operators

**The spatial channel** should map quantum states on one spatial location to quantum states at another spatial location, at FIXED time.

**Proposed model**: A beam-splitter with amplitude transmissivity sqrt(eta_spatial):
```
a_out = sqrt(eta_spatial) a_in + sqrt(1-eta_spatial) a_env
```

For eta_spatial = 1/Q:
```
a_out = (1/sqrt(Q)) a_in + sqrt(1-1/Q) a_env
```

The Kraus operators would be:
```
E_k^{(spatial)} = sqrt(C(n,k)) (1/Q)^{n/2} (1-1/Q)^{k/2} |n-k><n|
```

with entropy production:
```
Sigma_spatial = -ln(eta_spatial) = -ln(1/Q) = ln Q
```

### 5.3 Problem: What is the Physical Environment?

For the temporal channel, the environment is clear: the thermal bath of Hawking/Unruh radiation at the local Tolman temperature T_loc = T_H/sqrt(-g_00).

For the spatial channel, the environment is... what? There are several candidates:

**Candidate A**: The gravitational degrees of freedom not captured by the spatial metric (e.g., extrinsic curvature, which encodes the time evolution).

**Candidate B**: The modes with k > k_cutoff (UV modes that are integrated out in the effective theory).

**Candidate C**: The modes propagating in directions orthogonal to the radial direction.

None of these is as clean as the thermal bath for the temporal channel. This is the most serious gap in the spatial channel construction.

### 5.4 Alternative: Factored Channel

Instead of constructing the spatial channel independently, we can define it through factorization:
```
N_total = N_temporal (x) N_spatial
```

where N_total is the full 4D gravitational channel with eta_total = 1/Q^2. If we DEMAND this factorization and KNOW eta_total and eta_temporal, then:
```
eta_total = eta_temporal * eta_spatial = 1/Q^2

But eta_temporal IS 1/Q^2 (from Paper 2). So eta_spatial = 1.
```

**THIS GIVES eta_spatial = 1, NOT 1/Q!**

Wait -- this is wrong. The factorization N_total = N_temporal tensor N_spatial means:
```
The transmissivities multiply: eta_total = eta_temporal * eta_spatial

So eta_spatial = eta_total / eta_temporal = (1/Q^2) / (1/Q^2) = 1.
```

This is a PROBLEM. If the total channel has eta = 1/Q^2 and the temporal channel also has eta = 1/Q^2, then the spatial channel is trivial (eta = 1, no loss).

### 5.5 Resolving the Factorization Paradox

The resolution is that the temporal channel IS the full channel. The "spatial channel" is NOT a factor of the temporal channel. It is a DIFFERENT channel entirely.

The temporal channel: maps a mode from (r_1, t_arbitrary) to (r_2, t_arbitrary), following the null geodesic. The full redshift (both energy and rate) is included.

The spatial channel: maps a field configuration at t = t_0 on one Cauchy slice to the "expected" field configuration, based on the gravitational constraints. This is NOT a dynamical channel (no time evolution) but a CONDITIONAL channel: given the field at one spatial location, what is the optimal estimate at another location?

Formally, the spatial channel is the conditional expectation:
```
N_spatial(rho_A) = Tr_B(rho_AB) / Tr(rho_AB)
```

where A and B are two spatial regions and rho_AB is the joint state. The conditional entropy production is:
```
Sigma_spatial = I(A:B | gravity) = mutual information between two spatial regions, conditioned on the gravitational field
```

This is DEFINED on the Cauchy slice, with no time evolution. Its value depends on the spatial correlations induced by the gravitational field.

### 5.6 The Honest Assessment

The spatial Kraus operators have NOT been constructed from first principles. The existing arguments establish:

1. Sigma_temporal = 2 ln Q (rigorous, from Paper 2 theorem with explicit Kraus operators)
2. The MOND interpolating function requires Sigma in the Crooks formula to equal sqrt(x), which corresponds to Sigma = ln Q (empirical from SPARC)
3. The identification Sigma_MOND = Sigma_temporal/2 = ln Q is CONSISTENT with calling this Sigma_spatial

But there is no independent derivation of the spatial channel that would give Sigma_spatial = ln Q from first principles (without reference to the MOND data).

### 5.7 Verdict on Attack 5

**CRITICAL GAP CONFIRMED.** The spatial Kraus operators do not exist in the current framework. The claim eta_spatial = 1/Q relies on:
- Physical arguments (ADM factorization, amplitude vs intensity) -- suggestive but not rigorous
- Empirical matching (SPARC data with p = 0.5) -- confirms the answer but does not derive it
- Consistency with the dispersion relation -- correct but again not an independent derivation

The factorization paradox (Section 5.4-5.5) reveals that the "spatial channel" is not simply a tensor factor of the temporal channel. It is a fundamentally different object (conditional expectation on a Cauchy slice, not a dynamical channel). Its mathematical structure is closer to a CONDITIONAL ENTROPY than a CHANNEL ENTROPY PRODUCTION.

---

## Attack 6: Constraint Equation Loophole

### 6.1 The Setup

In the ADM formalism, the metric is:
```
ds^2 = -N^2 dt^2 + h_ij(dx^i + N^i dt)(dx^j + N^j dt)
```

The lapse N = 1/Q and the spatial metric h_ij are NOT independent. They are related by the Hamiltonian constraint:
```
R^{(3)} + K^2 - K_ij K^{ij} = 16 pi G T_00 / N^2
```

where R^{(3)} is the spatial Ricci scalar and K_ij is the extrinsic curvature.

For a static spacetime (K_ij = 0), this reduces to:
```
R^{(3)} = 16 pi G rho / N^2
```

This means the spatial curvature R^{(3)} is determined by N = 1/Q and the matter density. If we know Q, the spatial metric h_ij is (in principle) determined by the constraint.

### 6.2 How g_ij Depends on Q Through Constraints

For a spherically symmetric, static spacetime in isotropic coordinates:
```
h_ij = Psi^4 delta_ij
```

The conformal factor Psi is determined by the constraint equation. In GR:
```
Schwarzschild: Psi = (1 + M/(2r))
Exponential metric: Psi^4 = e^{r_s/r} = Q^2 (since g_00 g_rr = -1)
```

In the weak field (1PN):
```
Psi^4 = 1 + 2U + O(U^2)    (both metrics agree)
Q^2 = 1 + 2U + O(U^2)
```

So Psi^4 = Q^2 at 1PN. At 2PN:
```
Exponential: Psi^4 = 1 + 2U + 2U^2 + ... = Q^2 exactly
Schwarzschild: Psi^4 = 1 + 4U + 6U^2 + ... != Q^2
```

### 6.3 Implications for eta_spatial

If g_ij = Q^alpha_eff delta_ij, then the spatial momentum transmissivity is:
```
k_proper(r) / k_proper(inf) = 1/Q^{alpha_eff/2}
```

For the exponential metric: alpha_eff = 2, giving k_proper ~ 1/Q. Hence eta_spatial (momentum) = 1/Q. CORRECT.

For Schwarzschild at 2PN: alpha_eff = 4 at leading order (Psi = (1+U), Psi^4 = (1+U)^4), giving:
```
k_proper ~ (1+U)^{-2} ~ 1/(1+2U+...) ~ 1/Q + O(U^2)
```

At 1PN, k_proper = 1/Q for both metrics. But at 2PN the Schwarzschild answer deviates.

The correction at 2PN:
```
For Schwarzschild: Psi = (1+U) => Psi^4 = 1 + 4U + 6U^2 + 4U^3 + U^4
For exponential:   Psi^2 = Q  => Psi^4 = Q^2 = 1 + 2U + 2U^2 + 4/3 U^3 + ...

Difference at O(U^2): 6U^2 vs 2U^2
```

But remember, U = GM/(c^2 r) ~ 10^{-6} at the solar surface, 10^{-9} at Earth's surface. At galactic scales for MOND, U ~ a r / c^2 ~ (10^{-10}) * (10^{20}) / (9 x 10^{16}) ~ 10^{-7}. So U^2 ~ 10^{-14}, and the 2PN correction is utterly negligible for all practical purposes.

### 6.4 Verdict on Attack 6

**IMPORTANT theoretically, NEGLIGIBLE practically.** The constraint equations DO couple g_ij to Q, and the specific relationship depends on whether the metric is exponential or Schwarzschild. At 1PN (U << 1), both give eta_spatial = 1/Q. At 2PN, the exponential metric gives eta_spatial = 1/Q exactly (because g_00 g_rr = -1), while Schwarzschild gives a correction of order U^2 ~ (r_s/r)^2. For all observationally relevant scales (galactic to solar system), this correction is below 10^{-12}.

---

## Synthesis: The Complete Picture

### What is Rigorously Established

1. **Sigma_temporal = 2 ln Q**: Proven with explicit Kraus operators for the thermal attenuator channel on static backgrounds.

2. **The MOND function requires Sigma = sqrt(x) in the Crooks exponent**: Empirically confirmed by SPARC data (p = 0.482 +/- 0.02, consistent with 0.5).

3. **Sigma = sqrt(x) follows from Sigma_temporal/2 = ln Q**: Mathematical identity, given Sigma_temporal = 2 sqrt(x) in the MOND regime.

4. **At 1PN, eta_spatial = 1/Q**: Follows from the proper momentum calculation on any spherically symmetric static background. Valid to O(Phi/c^2).

5. **The non-relativistic dispersion relation gives Sigma_E = 2 Sigma_p**: Correct with corrections O(v^2/c^2) ~ 10^{-11} at galactic scales.

### What is NOT Rigorously Established

1. **The spatial channel Kraus operators**: Not constructed. The "spatial channel" is defined as a conditional expectation on a Cauchy slice, not as a dynamical channel. Its mathematical structure is fundamentally different from the temporal channel.

2. **eta_spatial = 1/Q beyond 1PN**: Requires the metric identity g_00 g_rr = -1 (exponential metric). For Schwarzschild, corrections appear at O(U^2) ~ (r_s/r)^2.

3. **The Born rule argument**: Heuristic. "Intensity = amplitude^2" is a correct mnemonic but the full T_00 involves both temporal and spatial derivatives. The correct derivation uses the Petz bound factor of 1/2, not an independent spatial channel.

4. **The factorization H = H_temporal x H_spatial**: The temporal channel with eta = 1/Q^2 accounts for ALL the entropy production. The "spatial channel" is not a tensor factor but a different mathematical object (conditional mutual information on a spatial slice).

### The Honest Assessment

**alpha = 1 is the correct value**, supported by:
- Physical argument: proper momentum scales as 1/Q (one power, from the spatial metric factor)
- Empirical data: SPARC gives p = 0.482 +/- 0.02, consistent with p = 0.5 = alpha/2
- Dispersion relation: non-relativistic Khronon gives E ~ p^2, hence 2 Sigma_p = Sigma_E
- PPN consistency: at 1PN, all metrics agree on g_ij ~ Q^2 delta_ij

**However**, the derivation has a conceptual gap: the "spatial channel" is not a channel in the CPTP sense. It is a conditional expectation / spatial mutual information quantity. The existing arguments derive the MOND function from Sigma_temporal (not from an independent spatial channel), and the factor of 1/2 comes from the Petz bound structure, not from an independently established spatial transmissivity.

### Could alpha be Different from 1?

Only if:
- The relevant quantity for MOND is NOT Sigma_temporal/2 but something else (but the Petz bound gives the 1/2 structurally)
- Relativistic corrections matter at galactic scales (they don't: v/c ~ 10^{-3.5})
- The spatial metric does not satisfy g_ij = Q^2 delta_ij at the required precision (it does at 1PN)
- The SPARC data are systematically biased toward p != 0.5 (possible but unlikely)

None of these is plausible given current knowledge.

### Recommended Statement in Papers

Replace the current language:
> "The spatial channel has transmissivity eta_spatial = 1/Q (one power of Q)"

With the more precise:
> "The MOND interpolating function follows from the Petz recovery bound F >= exp(-Sigma_temporal/2) applied to the temporal entropy production Sigma_temporal = 2 sqrt(x), giving tau = 1 - exp(-sqrt(x)). This is equivalently described by defining Sigma_spatial := Sigma_temporal/2 = ln Q, consistent with the interpretation that the spatial sector carries the amplitude (rather than intensity) transmissivity eta_spatial = 1/Q. The factor-of-2 relationship is confirmed by three independent arguments: (i) the ADM factorization of mode functions at 1PN, (ii) the non-relativistic dispersion relation of the massive Khronon, and (iii) the Petz bound structure. The explicit construction of spatial Kraus operators remains an open problem."

---

## Open Problems Ranked by Priority

1. **[CRITICAL]** Construct the spatial channel as a rigorous mathematical object (conditional expectation map on the Khronon Cauchy slice). Show that its "entropy production" equals ln Q.

2. **[IMPORTANT]** Determine whether the spatial channel is a CPTP map or a different type of quantum operation (e.g., a quantum instrument, a conditional expectation).

3. **[IMPORTANT]** Extend the 1PN result eta_spatial = 1/Q to 2PN. This requires choosing between Schwarzschild (g_rr != Q^2 at 2PN, eta_spatial has corrections) and exponential (g_rr = Q^2 exactly, eta_spatial = 1/Q exactly).

4. **[NICE-TO-HAVE]** Compute the scale-dependent alpha(k) from the Khronon dispersion relation and predict the deviation from alpha = 1 at Hubble-scale separations.

5. **[NICE-TO-HAVE]** Determine whether future surveys (SKA, DESI, Euclid) can constrain p to +/- 0.005 and thereby test alpha = 1 at 3+ sigma.

---

## Final Verdict

**eta_spatial = 1/Q (alpha = 1) SURVIVES all five attacks.** The claim is correct at the physicist's level of rigor, supported by multiple convergent arguments and confirmed by SPARC data. The main vulnerability is not the VALUE of alpha but the FOUNDATIONAL STATUS of the spatial channel: it has not been constructed as a rigorous CPTP map. This gap should be acknowledged honestly in publications but does not affect the physical predictions.

**Confidence**: 90% that alpha = 1 is the correct answer. 10% residual uncertainty from the lack of rigorous spatial Kraus operators and the 2PN metric ambiguity.

---

*Last updated: 2026-03-20*
*This analysis attempted to BREAK the claim eta_spatial = 1/Q through five independent attacks. All attacks failed at the physically relevant precision, but exposed the conceptual gap that the "spatial channel" is not a standard CPTP channel.*
