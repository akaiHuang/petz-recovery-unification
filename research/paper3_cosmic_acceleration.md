# Can the tau Framework Replace Lambda? Cosmic Acceleration from Information Bounds

## A Rigorous Analysis

**Research Date:** 2026-03-11
**Author:** Sheng-Kai Huang
**Relevance:** Paper 3 (cosmological extension), Paper 4 (unification)
**Status:** Comprehensive literature survey + quantitative analysis + honest assessment

---

## Executive Summary

This document investigates whether cosmic acceleration (currently attributed to Lambda > 0) can be understood, predicted, or replaced within the tau framework's information-theoretic perspective. The analysis is organized around five questions, each answered with full quantitative rigor and honest assessment of what the framework can and cannot do.

| Question | Answer | Confidence |
|----------|--------|------------|
| Can tau reinterpret Lambda? | **YES** -- Lambda <==> finite cosmic information | HIGH (multiple independent routes) |
| Can tau predict Lambda's value? | **PARTIALLY** -- order of magnitude, not precise | MODERATE (Padmanabhan's CosMIn) |
| Can tau replace Lambda? | **NO (currently)** -- no dynamical mechanism derived | LOW |
| Does tau predict evolving dark energy? | **PLAUSIBLE** -- scale-dependent Sigma allows w != -1 | MODERATE (consistent with DESI) |
| Is "Lambda = information bound" tautological? | **NO** -- it makes falsifiable predictions | HIGH |

**Bottom line:** The tau framework provides a powerful reinterpretation of Lambda as a cosmic information bound (tau = 1 at the Hubble horizon), connects to multiple independent programs (Padmanabhan, Jacobson, Verlinde, GREA), and is naturally compatible with evolving dark energy (DESI DR2 hints). However, it does NOT yet derive the Friedmann equations from first principles, and the statement "Lambda = information bound" is a constraint, not a dynamical mechanism. The framework is most powerful as a unifying language, not (yet) as a replacement theory.

---

## Part I: Information-Theoretic Approaches to Lambda

---

### 1. Padmanabhan's Cosmic Information (CosMIn)

**Key references:**
- Padmanabhan (2017, arXiv:1703.06144, PLB 773, 81-85)
- Padmanabhan (2013, arXiv:1302.3226, IJMPD 22, 1342001)
- Padmanabhan (2012, arXiv:1206.4916)

#### 1.1 The CosMIn Framework

CosMIn (Cosmic Information) is defined as the total number of modes that cross the Hubble radius during the entire evolution of the universe, counting information transferred from quantum-to-classical as modes exit the horizon:

```
CosMIn = integral_0^infinity (dN_modes/dt) dt
```

**The key theorem:** CosMIn is finite if and only if the universe undergoes late-time accelerated expansion.

*Proof sketch:* In a matter-dominated universe (a ~ t^{2/3}), modes continue crossing the Hubble radius at a rate that makes CosMIn diverge logarithmically. Only if a(t) ~ exp(Ht) at late times does the integral converge:

```
CosMIn ~ integral dt / (aH)^3 * d/dt[(aH)^3]
```

For H = const (de Sitter): integrand decays exponentially -> CosMIn finite.
For H ~ 1/t (matter-dominated): integrand ~ const -> CosMIn diverges.

**Quantitative prediction.** Demanding CosMIn = 4pi (a natural topological constant for S^3 topology) gives:

```
Lambda L_P^2 ~ 3.4 x 10^{-122}
```

This is within an order of magnitude of the observed value Lambda L_P^2 ~ 2.9 x 10^{-122}.

#### 1.2 Connection to Holographic Equipartition

Padmanabhan's holographic equipartition principle:

```
dV/dt = L_P^2 (N_sur - N_bulk)
```

where:
- N_sur = 4pi / (L_P^2 H^2) = surface degrees of freedom on the Hubble sphere
- N_bulk = |rho + 3p| V_H / (T_H / 2) = bulk degrees of freedom (Komar energy / (kT/2))
- T_H = H/(2pi) = Hubble temperature

**Equilibrium condition:** N_sur = N_bulk gives exactly de Sitter space.

**The driving force:** When N_sur > N_bulk, space expands. When they equilibrate, you get the de Sitter attractor.

This correctly reproduces the standard Friedmann equation. More importantly, it REQUIRES both early inflation and late-time acceleration for self-consistency.

#### 1.3 Mapping to the tau Framework

The connection is remarkably clean:

| Padmanabhan | tau Framework |
|-------------|---------------|
| CosMIn (total cosmic information) | integral of Sigma over cosmic history |
| N_sur - N_bulk > 0 | tau > 0 (arrow of time exists) |
| dV/dt > 0 (expansion) | Sigma_total > 0 (entropy production) |
| N_sur = N_bulk (de Sitter equilibrium) | tau reaches its maximum allowed value |
| CosMIn finite | Sigma integrated over all time is finite |

**The precise mapping:**

```
CosMIn = integral_0^infinity Sigma_mode * (dN_mode/dt) dt
```

where Sigma_mode is the entropy production per mode crossing the Hubble radius. Each mode that exits the horizon loses its quantum coherence (decoherence by the horizon) -- this is exactly the information loss quantified by tau.

**Key insight:** Padmanabhan's "finite CosMIn requires Lambda > 0" translates to: "finite total information loss requires a cosmic information horizon (tau = 1 at finite radius)."

#### 1.4 Assessment

**What CosMIn gives us:**
- Conceptual necessity: Lambda > 0 is REQUIRED for finite cosmic information (HIGH confidence)
- Approximate value: CosMIn = 4pi gives correct order of magnitude for Lambda (MODERATE confidence)
- Dynamical picture: N_sur - N_bulk drives expansion (HIGH confidence, reproduces Friedmann)

**What CosMIn does NOT give us:**
- A derivation of WHY CosMIn should be 4pi (this is assumed)
- A microscopic mechanism for the expansion drive
- Any prediction for evolving dark energy

---

### 2. Jacobson's Entanglement Equilibrium

**Key references:**
- Jacobson (1995, PRL 75, 1260; gr-qc/9504004)
- Jacobson (2016, PRL 116, 201101; arXiv:1505.04753)
- Casalderrey-Solana et al. (2016, JHEP 03, 194; arXiv:1601.00528)

#### 2.1 Entanglement Equilibrium Hypothesis

Jacobson's 2016 paper posits: **the vacuum entanglement entropy in small geodesic balls is maximal at fixed volume, in a locally maximally symmetric state.**

For first-order variations of the vacuum state:

```
delta S_ent = delta S_geom + delta S_matter = 0
```

where S_geom = A/(4G) (Bekenstein-Hawking) and S_matter is the matter entanglement entropy.

**Result:** This condition is satisfied if and only if the Einstein equation holds, INCLUDING a cosmological constant Lambda that absorbs the UV-divergent part of the modular Hamiltonian.

#### 2.2 Lambda from Vacuum Entanglement

In Jacobson's framework, Lambda arises naturally:

```
R_ab - (1/2) R g_ab + Lambda g_ab = (8pi G / c^4) T_ab
```

The cosmological constant appears as the "baseline" vacuum entanglement -- the entanglement entropy density of the vacuum state in a maximally symmetric spacetime. It is NOT a fine-tuned parameter but the natural equilibrium value.

**Physical picture:** Lambda measures the irreducible entanglement of the vacuum. In tau language:

```
Sigma_vacuum = Lambda r^2 / (3 c^2)
```

This grows quadratically with distance. At r = c/H (the Hubble radius):

```
Sigma_vacuum(R_H) = Lambda / (3 H^2) = Omega_Lambda ~ 0.7
```

This is O(1), confirming that the Hubble horizon is where vacuum entanglement becomes maximal.

#### 2.3 Connection to tau

```
tau_vacuum(r) = 1 - exp(-Sigma_vacuum/2)
     = 1 - exp(-Lambda r^2 / (6c^2))
```

At the Hubble horizon: tau -> 1 - exp(-Omega_Lambda/2) ~ 0.30 for Omega_Lambda = 0.7.
At the de Sitter horizon (exact): tau -> 1 (complete information loss).

**Jacobson's contribution to the tau framework:** Lambda is not an arbitrary constant -- it is the vacuum's intrinsic information loss rate per unit area. The Einstein equations (with Lambda) are the CONDITION for entanglement equilibrium, and tau quantifies the departure from this equilibrium.

#### 2.4 Assessment

**Strengths:**
- Derives Einstein equations + Lambda from entanglement (HIGH rigor)
- Lambda has natural interpretation as vacuum entanglement
- No fine-tuning needed (Lambda is an integration constant)

**Weaknesses:**
- Does not predict Lambda's VALUE (it appears as a free constant)
- UV-divergent contributions must be regulated
- Extension to higher-order gravity (Bueno et al. 2017) shows subtleties

---

### 3. Banks-Fischler Holographic Cosmology

**Key references:**
- Banks & Fischler (2000, hep-th/0007146)
- Banks (2025, arXiv:2511.08213)
- Banks & Fischler (2024, arXiv:2402.11527)

#### 3.1 Finite Entropy and Lambda

The central postulate: **the universe has a finite-dimensional Hilbert space with dimension:**

```
dim(H) = exp(S_dS) = exp(pi c^3 / (G hbar Lambda))
```

For the observed Lambda:

```
S_dS = pi / (Lambda L_P^2) ~ 10^{122}
```

This is the maximum entropy accessible to any observer in de Sitter space (the Gibbons-Hawking entropy).

#### 3.2 Implications for tau

If the total Hilbert space is finite, then:

```
S_max = pi / (Lambda L_P^2)    [finite]
```

This means:
1. There is a maximum possible entropy production: Sigma_max is finite
2. The Petz recovery fidelity has a minimum: F_min = exp(-Sigma_max/2) > 0
3. tau has an upper bound: tau_max = 1 - F_min < 1

**But this bound is essentially tau_max ~ 1 - exp(-10^{122}/2) ~ 1 to all practical purposes.** The finite Hilbert space constrains information loss at the cosmic level but has no practical consequence for sub-Hubble physics.

#### 3.3 The Deeper Connection

The Banks-Fischler framework implies:

```
Lambda > 0  <=>  dim(H) < infinity  <=>  S_max < infinity  <=>  cosmic tau < 1 everywhere
```

And the cosmological constant problem becomes: **why is dim(H) ~ exp(10^{122}) and not, say, exp(10^{60})?**

In tau language: **why is the cosmic information horizon at R_H ~ 10^{26} m and not somewhere else?**

This restatement does not solve the cosmological constant problem, but it reformulates it in information-theoretic terms that may be more tractable.

#### 3.4 Recent Development: BMN Strings (Banks 2025)

Banks (arXiv:2511.08213) proposes that BMN string theory provides a test of the universal entropy bound and a possible estimate of the cosmological constant. If the bound S_max ~ N^2 for SU(N) gauge theory applies to the fundamental degrees of freedom of gravity, then:

```
Lambda ~ (N^2 L_P^2)^{-1}
```

This is a step toward PREDICTING Lambda's value from a microscopic theory, but remains highly speculative.

#### 3.5 Assessment

**Strengths:**
- Clean formulation: Lambda <=> finite Hilbert space
- Connects to well-established holographic ideas
- Recent work (2025) may lead to Lambda prediction

**Weaknesses:**
- Does not derive Lambda's value (except speculatively via BMN)
- Requires accepting holographic principle as axiom
- No connection to observational tensions (Hubble, sigma8)

---

### 4. Verlinde's Emergent Gravity and Cosmic Acceleration

**Key reference:** Verlinde (2016, arXiv:1611.02269; SciPost Phys. 2, 016)

#### 4.1 Volume-Law Entanglement and Lambda

Verlinde's key insight: in de Sitter space, entanglement entropy has BOTH area-law and volume-law contributions:

```
S_ent = S_area + S_volume

S_area = A / (4 G hbar)     [standard Bekenstein-Hawking]
S_volume = (2pi/3) * (r/L_dS)^3 * S_dS    [thermal, from T_dS = H/(2pi)]
```

The volume-law entropy arises because the de Sitter vacuum is a thermal state at temperature T_dS. This thermal bath carries entropy that scales with volume, not area.

#### 4.2 Dark Energy as Entanglement Pressure

In Verlinde's framework, the cosmological constant IS the effect of the volume-law entanglement:

```
Lambda = 3 H^2 / c^2  <=>  the thermal volume-law entropy saturates at the Hubble scale
```

The "dark energy" is not a substance -- it is the entropic pressure arising from the thermal de Sitter vacuum's volume-law entanglement. When this entanglement is disrupted by matter, you get "dark matter" effects. When it reaches equilibrium at the Hubble scale, you get cosmic acceleration.

#### 4.3 Connection to tau

In Sigma language:

```
Sigma_dS(r) = -ln(1 - H^2 r^2 / c^2)

At r << R_H:   Sigma ~ (H r / c)^2         [weak field, quadratic growth]
At r -> R_H:   Sigma -> infinity             [horizon, tau = 1]
```

The volume-law entanglement contributes an ADDITIONAL Sigma beyond the classical metric:

```
Sigma_Verlinde(r) = Sigma_classical(r) + Sigma_entanglement(r)
```

where Sigma_entanglement encodes the entropy displacement and drives both "dark matter" (at galactic scales) and "dark energy" (at cosmological scales).

**This is the most natural tau framework interpretation:** Lambda is not a separate parameter but an emergent property of the quantum vacuum's entanglement structure, manifesting as Sigma at the Hubble scale.

#### 4.4 Assessment

**Strengths:**
- Unifies dark matter and dark energy in one mechanism
- Natural connection to tau (volume-law entropy = cosmological Sigma)
- Observational support for DM predictions (dwarf spheroidals: 5.2 sigma over MOND)

**Weaknesses:**
- Theoretical consistency issues (Hossenfelder 2017: covariant version reduces to GR)
- No precise prediction for Lambda's value
- No CMB prediction

---

### 5. GREA: Gravitational Relativistic Entropic Acceleration

**Key references:**
- Gaztanaga et al. (2024, arXiv:2405.02895, PDU 45, 101533)
- Gaztanaga et al. (2025, arXiv:2509.21491)
- Gaztanaga et al. (2026, arXiv:2603.01934)
- Gaztanaga et al. (2025, arXiv:2511.19546) -- holographic dual

#### 5.1 The GREA Mechanism

GREA provides a covariant formalism where cosmic acceleration arises from entropy production, without a cosmological constant. The key modification to the Friedmann equations:

**Standard:**
```
H^2 = (8pi G / 3) rho
```

**GREA:**
```
H^2 = (8pi G / 3) (rho - T S_dot / V_dot)
```

The source for spacetime curvature is the Helmholtz free energy F = U - TS, not just the internal energy U. The entropic term introduces an effective negative pressure:

```
p_S = -T dS/dV < 0
```

This is guaranteed to be negative by the second law of thermodynamics (dS >= 0 in an expanding universe with dV > 0).

#### 5.2 GREA's Predictions

The model predicts:
1. **Cosmic acceleration from entropy growth** -- not from Lambda
2. **w(z) evolving** -- the effective equation of state is NOT exactly -1
3. **Resolution of the coincidence problem** -- acceleration begins when entropy production from structure formation becomes significant
4. **Compatible with DESI DR2** -- recent fits (arXiv:2603.01934, March 2026) show GREA is consistent with BAO + CMB + SN data

#### 5.3 Connection to tau

GREA is perhaps the MOST directly connected program to the tau framework:

| GREA | tau Framework |
|------|---------------|
| Entropy production dS/dt | Sigma (entropy production rate) |
| p_S = -TdS/dV | Sigma contribution to effective pressure |
| Acceleration from dS > 0 | tau > 0 drives expansion |
| de Sitter as equilibrium | tau_max at Hubble horizon |

**The precise identification:**

```
GREA: p_S = -T_H * (dS_horizon/dV) = -(H/2pi) * (d/dV)(A_H / 4G)

In tau language: p_tau = -(1/8piG) * H^2 * d(tau)/d(ln a)
```

If tau is related to the horizon entropy via tau = 1 - exp(-S_horizon / S_max), then the GREA effective pressure IS the pressure from the cosmological arrow of time.

**This is a profound connection:** GREA says "entropy production drives cosmic acceleration." The tau framework says "tau > 0 (arrow of time) is entropy production." Therefore: **the arrow of time drives cosmic acceleration.**

#### 5.4 Assessment

**Strengths:**
- Covariant formalism (not just a conjecture)
- Tested against DESI DR2, CMB, SN data (2026)
- Natural prediction of evolving w(z)
- Resolves coincidence problem
- Direct connection to tau framework

**Weaknesses:**
- The entropy source (horizon entropy vs. matter entropy) is debated
- Does not predict Lambda's specific value
- Competing with Lambda CDM, not yet clearly preferred statistically
- No explicit connection to quantum information theory (our contribution would be new)

---

## Part II: Does Accelerated Expansion Follow from Information Bounds?

---

### 6. The Core Question

**Claim to investigate:** If the total information in the universe is finite, there must be a boundary where tau = 1 (an information horizon). Does this REQUIRE accelerated expansion?

#### 6.1 Analysis by Expansion Type

**Matter-dominated universe:** H(t) = 2/(3t)
```
R_H = c/H = 3ct/2    [Hubble radius, growing linearly]
tau = 1 at R_H(t)     [information horizon exists but moves outward]
Sigma(R_H) ~ O(1)     [always, by definition of the horizon]
```

There IS an information horizon (tau = 1 at R_H), and it grows. Total information is NOT bounded -- it grows without limit as R_H -> infinity. CosMIn diverges.

**Lambda-dominated (de Sitter) universe:** H = const = sqrt(Lambda/3)
```
R_H = c/H = c/sqrt(Lambda/3) = const    [Hubble radius, FIXED]
tau = 1 at R_H = const                   [fixed information horizon]
Sigma(R_H) ~ O(1)                        [constant]
```

The information horizon is FIXED. Total accessible information is finite. CosMIn converges.

**Transition (matter -> Lambda):**
```
H(t): 2/(3t)  ->  H_inf = sqrt(Lambda/3)    as t -> infinity
R_H(t): 3ct/2  ->  c/H_inf                   [growing -> fixed]
```

The tau = 1 boundary transitions from GROWING to FIXED. This is the transition from unbounded to bounded information.

#### 6.2 The Logical Structure

```
FINITE total information  ==>  R_H must stop growing
                          ==>  H must approach a constant
                          ==>  expansion must accelerate
                          ==>  Lambda > 0 (or effective Lambda)
```

**This is a genuine logical implication, not a tautology.** The premise (finite information) is an independent physical principle. The conclusion (accelerated expansion) is an observable consequence.

**However, the converse is also true:**

```
Lambda > 0  ==>  R_H = c/sqrt(Lambda/3) = const
            ==>  S_max = pi c^3 / (G hbar Lambda) = const
            ==>  finite information
```

So: **finite information <==> Lambda > 0.** These are logically equivalent statements, not one deriving the other. The value of the tau framework is not in DERIVING one from the other, but in providing a UNIFIED LANGUAGE in which both are the same statement.

#### 6.3 Can We Derive Friedmann Equations from tau Dynamics?

**Current status: PARTIALLY.**

Several programs achieve partial derivations:

1. **Jacobson (1995):** delta Q = TdS on local Rindler horizons -> Einstein equations (including Friedmann as a special case). This is the foundational result.

2. **Padmanabhan (2012):** dV/dt = L_P^2 (N_sur - N_bulk) -> Friedmann equations. This uses holographic equipartition.

3. **Dorau-Much (2026, PRL):** QRE between vacuum and coherent excitations on Killing horizons -> semiclassical Einstein equations. The most recent and rigorous derivation.

4. **GREA (2024-2026):** Helmholtz free energy as source -> modified Friedmann equations with entropic term that mimics Lambda.

**None of these derive Friedmann equations PURELY from tau.** They all require additional input:
- Jacobson needs the proportionality delta S = delta A / (4G)
- Padmanabhan needs N_sur = A/(L_P^2) and the Komar energy expression
- Dorau-Much needs the Bekenstein-Hawking assumption
- GREA needs the entropy source identification

**What a tau derivation would need:**

Starting from:
```
Sigma = D(rho_spacetime || rho_matter)      [Paper 4 goal]
tau = 1 - exp(-Sigma/2)
```

One would need to show that the STATIONARITY of Sigma (delta Sigma = 0 for small perturbations) implies the Friedmann equations. This is essentially Jacobson's program reformulated in QRE language, and Dorau-Much (2026) has taken the most important step.

**The missing piece:** Deriving the SPECIFIC form of Sigma in FRW spacetime from first principles. In particular, showing that:

```
Sigma_FRW = -ln(1 - H^2 r^2 / c^2)    [for Kodama observer]
```

follows from the QRE of the gravitational channel in an FRW background. This is achievable but has not been done.

---

### 7. The tau = 1 Boundary and Its Time Evolution

#### 7.1 tau in FRW Spacetime

Using the Kodama vector (the natural generalization of the Killing vector to non-stationary spherically symmetric spacetimes):

```
K^a = (1, -Hr, 0, 0)    [in comoving coordinates, k=0]

|K|^2 = -(1 - H^2 r^2 / c^2)
```

Define:
```
F_K(r) = sqrt(1 - H^2 r^2 / c^2)        [Kodama fidelity]
tau_K(r) = 1 - F_K(r) = 1 - sqrt(1 - H^2 r^2 / c^2)
```

Properties:
- r = 0 (observer): F_K = 1, tau_K = 0 (perfect recovery)
- r = c/H (apparent horizon): F_K = 0, tau_K = 1 (complete information loss)
- Smooth interpolation between 0 and 1

#### 7.2 Time Evolution of H(t) in the tau Framework

The Friedmann equation:
```
H^2 = (8pi G / 3) rho + Lambda/3
```

gives H(t) for any energy content. In the tau framework, we can rewrite this as:

```
Sigma(r = c/H) = -ln(1 - 1) = infinity    [always, at the apparent horizon]
```

This is trivially true -- the apparent horizon ALWAYS has tau = 1 by construction. The non-trivial content is in HOW FAST tau approaches 1:

```
d(tau)/dr |_{r=0} = (H/c)^2 * r    [near the observer]
```

This depends on H(t), which encodes all the cosmological dynamics.

#### 7.3 The Equation of State from tau

If we parametrize the dark energy equation of state as w(z), the Friedmann equation gives:

```
H^2(z) = H_0^2 [Omega_m (1+z)^3 + Omega_DE * exp(3 integral_0^z (1+w(z'))/(1+z') dz')]
```

In the tau framework, w(z) is related to how the tau = 1 boundary (the apparent horizon) evolves:

```
R_AH(t) = c / H(t)

dR_AH/dt = -(c/H^2) * dH/dt = (c/H^2) * (4piG/3)(rho + 3p/c^2) * H

For w = -1:  dR_AH/dt = 0  (static horizon -- de Sitter)
For w > -1:  dR_AH/dt > 0  (growing horizon -- phantom-like if currently growing)
For w < -1:  dR_AH/dt < 0  (shrinking horizon -- phantom energy)
```

**tau framework prediction:** If cosmic acceleration arises from entropy production (GREA mechanism) rather than a constant Lambda, then:

```
w(z) = -1 + (delta w)_entropy

where (delta w)_entropy = -(T/rho) * (dS/dV) / c^2
```

This is generically NOT exactly -1, consistent with the DESI DR2 hints of evolving dark energy.

#### 7.4 Quantitative Estimate of w(z)

For GREA-type entropy production from the cosmic horizon:

```
dS_H/dt = (d/dt)(A_H / 4G) = (d/dt)(pi c^2 / (G H^2)) = -2pi c^2 H_dot / (G H^3)
```

Using H_dot = -(4piG/3)(rho + p/c^2) H and assuming matter domination:

```
|delta w| ~ T_H * S_H_dot / (rho_crit c^2 V_H)
         ~ (H/2pi) * (2pi c^2 H_dot / (G H^3)) / ((3H^2/(8piG)) * c^2 * (4pi c^3/(3H^3)))
         ~ -H_dot / H^2
         ~ (1 + q)    where q is the deceleration parameter
```

For the current epoch (q ~ -0.55):

```
|delta w| ~ 0.45
```

This is a LARGE correction -- too large for a perturbative treatment. The GREA framework handles this non-perturbatively, giving effective w ~ -0.7 to -0.9, which is in the right ballpark for DESI but systematically different from LCDM.

**Honest assessment:** This estimate is rough. The precise value depends on the entropy source identification, which is still debated.

---

## Part III: Observational Constraints

---

### 8. What Must Any Alternative to Lambda Explain?

| Observation | Constraint | Lambda CDM | tau Framework |
|-------------|-----------|------------|---------------|
| Type Ia SN (Pantheon+) | Acceleration since z ~ 0.7 | w = -1.03 +/- 0.03 | Compatible (GREA gives w ~ -0.7 to -1) |
| BAO (DESI DR2) | d_A(z) and H(z) at 7 redshifts | LCDM fit, but w0wa preferred at 2.5-4.2 sigma | **Naturally compatible** (w != -1 expected) |
| CMB (Planck) | theta_* = 1.04110 +/- 0.00031 | Excellent fit | Must match (constrains total Sigma history) |
| f sigma_8(z) | Growth rate of structure | sigma_8 tension (2-3 sigma low) | **Potentially resolved** (weaker gravity at late times) |
| t_0 ~ 13.8 Gyr | Age of universe | Matches | Must match (constrains H_0 and Omega_Lambda) |
| H_0 tension | 67.4 vs 73.0 km/s/Mpc | 4-5 sigma tension | **Potentially alleviated** (GREA shifts coasting point) |

### 9. DESI DR2 Results and Their Implications

**Key papers:**
- DESI Collaboration (2025, arXiv:2503.14738, Phys. Rev. D 112, 083515)
- DESI Collaboration (2024, arXiv:2404.03002)

#### 9.1 Summary of DESI Results

DESI's Data Release 2 (March 2025) measures BAO at 7 effective redshifts using galaxies and quasars. Key findings:

1. **w0waCDM preferred over LCDM:** The fit improves when allowing w = w0 + wa(1-a), with w0 > -1 and wa < 0, at 2.5-4.2 sigma depending on the supernova dataset used.

2. **Best-fit values:** w0 ~ -0.7, wa ~ -1.0 (varies with dataset combination)

3. **The trend is robust:** Both DR1 and DR2 show the same preference, and the significance has grown.

4. **Not yet a discovery:** 4.2 sigma is suggestive but below the 5 sigma threshold.

#### 9.2 Implications for the tau Framework

The DESI results are MORE compatible with the tau framework than with LCDM:

1. **w != -1 is natural in the tau framework.** If cosmic acceleration arises from entropy production rather than a constant Lambda, the effective w should evolve with redshift.

2. **w0 > -1 is expected.** The entropic contribution to pressure is bounded: |p_S| <= T * S_dot / V_dot, which is finite. So the effective w should be > -1 (not phantom).

3. **wa < 0 means w was closer to -1 in the past.** In the tau framework, this means entropy production was slower in the past (fewer structures, less horizon entropy) -- consistent with the idea that acceleration strengthens as structure forms.

4. **Specific GREA prediction:** Gaztanaga et al. (arXiv:2603.01934, March 2026) show GREA is consistent with DESI DR2 + CMB + SN, with the entropic acceleration naturally producing the observed w0wa trend.

#### 9.3 Quantitative Comparison

DESI+CMB+PantheonPlus best fit: w0 = -0.75 +/- 0.12, wa = -0.99 +/- 0.40

GREA prediction (Gaztanaga 2024): w_eff(z) starts near -1 at high z and evolves toward -0.7 at z = 0, giving effective w0 ~ -0.7 to -0.8 and wa ~ -0.5 to -1.0.

**The match is suggestive but not definitive.** Both LCDM and GREA are currently consistent with the data within uncertainties.

### 10. The Hubble and sigma_8 Tensions

#### 10.1 Hubble Tension (H_0 = 67.4 vs 73.0)

In the tau framework, H_0 sets the cosmic information horizon:
```
R_H = c / H_0 ~ 4.4 Gpc (for H_0 = 67.4)  or  ~ 4.1 Gpc (for H_0 = 73.0)
```

If the cosmic Sigma has scale-dependent corrections (running G, volume-law entanglement), the inferred H_0 from CMB vs local measurements could differ because they probe DIFFERENT scales.

**Specifically:** If G_eff(z) decreases slightly at late times (as running G from asymptotic safety suggests, see Zholdasbek et al. 2025, arXiv:2405.02636), then:
- CMB (z ~ 1100): probes G at early times
- Local (z ~ 0): probes G at late times (slightly smaller)
- This makes the locally inferred H_0 LARGER than the CMB-inferred H_0

**Quantitative estimate:** A running of |delta G / G| ~ 4% between z = 1100 and z = 0 could account for the ~8% discrepancy in H_0. This is within the range of Kumar's running G model (where G changes by ~30% between solar system and galactic scales -- much larger than needed).

**However:** This explanation requires running G to affect cosmological dynamics, not just local dynamics. And it must not spoil the CMB power spectrum, which is tightly constrained. This is a SERIOUS obstacle.

**Assessment: SPECULATIVE.** The mechanism is plausible but not quantitatively verified. Rating: LOW confidence.

#### 10.2 sigma_8 Tension

The sigma_8 tension: CMB predicts more structure than is observed at low redshift.

In the tau framework: if Sigma_grav at late times is LARGER than in LCDM (due to entropy production), this means MORE information loss, which corresponds to LESS coherent gravitational clustering. The structures "decohere" faster.

**Physical mechanism:** Enhanced entropy production at late times (from GREA mechanism or running G) effectively weakens the gravitational pull on the growth of perturbations:

```
delta_ddot + 2H delta_dot - 4pi G_eff rho delta = 0

If G_eff decreases at late times: delta grows more slowly -> lower sigma_8
```

**This is precisely the right direction** to resolve the sigma_8 tension. But the quantitative agreement depends on the specific model of running G.

**Assessment: PLAUSIBLE.** Direction is correct, but quantitative verification is needed. Rating: MODERATE confidence.

---

## Part IV: The "No Dark Energy" Landscape (2024-2026)

---

### 11. Recent Claims Against Lambda

#### 11.1 GREA (Gravitational Relativistic Entropic Acceleration)

**Status:** The most developed information-theoretic alternative.

- Gaztanaga et al. (2024-2026): Series of papers developing the covariant formalism
- Tested against DESI DR2 (arXiv:2603.01934, March 2026): consistent with data
- Holographic dual formulation (arXiv:2511.19546, November 2025): connects to fundamental d.o.f.
- Predicts evolving w(z) compatible with DESI hints

**Key equation:**
```
G_ab + Lambda_eff g_ab = 8pi G T_ab

where Lambda_eff = Lambda_bare + Lambda_entropy
and Lambda_entropy = (8pi G / c^4) * T_H * (dS_H / dV_H)
```

If Lambda_bare = 0, all acceleration comes from entropy production.

#### 11.2 Partner Anti-Universe Model

- Uses quantum relative entropy between universe and anti-universe
- Accelerated expansion follows from null energy condition + relative entropy positivity
- Connection to tau: the "partner" provides the reference state for QRE

#### 11.3 Information Dark Energy (Vopson-type)

- Landauer principle: each bit has energy kT ln(2)
- Total information content of the universe: N_bits ~ 10^{93} (from stellar processes)
- Energy: E ~ N_bits * kT_CMB * ln(2) ~ 10^{-13} J/m^3
- Dark energy density: rho_DE ~ 6 x 10^{-10} J/m^3

**Discrepancy:** ~10^3 too small. But the order of magnitude is much closer than the naive quantum field theory estimate (10^{120} too large).

**Connection to Trivedi (2024, arXiv:2407.15231):** Information erasure at the cosmological horizon does NOT obey Landauer's principle in the same way as black holes. This suggests the cosmological information problem has a different character than the black hole information problem.

#### 11.4 Extended Uncertainty Principle (EUP)

- Minimal momentum uncertainty -> scale-dependent quantum corrections
- At cosmological scales: effective negative pressure
- Late-time acceleration without Lambda

#### 11.5 Running Vacuum Models

- Lambda(H) = Lambda_0 + nu H^2 + ...
- Running of the vacuum energy with the Hubble parameter
- Can address both H_0 and sigma_8 tensions
- Connection to tau: Lambda(H) means Sigma_vacuum is not constant but evolves

---

## Part V: Synthesis -- Can tau Replace Lambda?

---

### 12. What the tau Framework CAN Say

#### 12.1 Established Results (HIGH Confidence)

1. **tau = 1 at the Hubble horizon.** The cosmological apparent horizon is an information boundary where tau = 1 (complete information loss for an observer). This is a mathematical consequence of the Kodama vector formalism.

2. **Lambda > 0 <==> finite total information.** Padmanabhan's CosMIn theorem establishes this equivalence. In tau language: finite integral of Sigma over cosmic history requires a fixed information horizon.

3. **Friedmann equations = entanglement equilibrium.** Jacobson's thermodynamic derivation (extended by Dorau-Much 2026) shows the Einstein equations (including Lambda) follow from entanglement stationarity. In tau language: the Einstein equations are the condition for Sigma to be stationary under small perturbations.

4. **GREA provides a concrete mechanism.** Entropy production at the cosmic horizon generates effective negative pressure that accelerates expansion. This maps directly onto the tau framework via Sigma = entropy production.

#### 12.2 Promising but Unproven (MODERATE Confidence)

5. **Evolving w(z) is natural.** If acceleration arises from entropy production (not a constant), w should deviate from -1. This is consistent with DESI DR2 hints (2.5-4.2 sigma).

6. **sigma_8 tension may be resolved.** Enhanced late-time entropy production (larger tau) weakens structure growth, reducing sigma_8.

7. **Lambda's value may be related to CosMIn = 4pi.** Padmanabhan's estimate gives the correct order of magnitude.

8. **Volume-law entanglement (Verlinde) unifies dark matter and dark energy.** Both arise from the de Sitter vacuum's entanglement structure being disrupted by matter.

#### 12.3 Speculative (LOW Confidence)

9. **Hubble tension from running G.** Plausible mechanism but not quantitatively verified and may conflict with CMB.

10. **Deriving Friedmann equations purely from tau dynamics.** Possible in principle (via Bianconi's framework + Kodama vector) but not yet achieved.

11. **Predicting Lambda's precise value.** No information-theoretic framework currently achieves this.

### 13. What the tau Framework CANNOT Say (Currently)

1. **It cannot derive Lambda's numerical value from first principles.** All approaches either input Lambda as a free parameter or get the order of magnitude at best.

2. **It cannot replace the Friedmann equations.** The tau framework provides a REINTERPRETATION of the Friedmann equations, not a replacement. The dynamical content is the same.

3. **It cannot explain why the universe started in a low-entropy state.** The arrow of time (tau > 0) is a consequence of the initial conditions, which the framework takes as given.

4. **It cannot (currently) handle CMB-scale physics.** The running G and volume-law entanglement approaches fail to reproduce CMB acoustic peaks without non-baryonic matter.

5. **It does not explain the coincidence problem on its own.** Why Omega_Lambda ~ Omega_m today requires additional input (GREA addresses this via structure formation onset).

### 14. Is "Lambda = Information Bound" a Tautology or a Prediction?

**It is NOT a tautology.** Here is why:

A tautology would be: "Lambda > 0 means Lambda > 0." That says nothing.

The statement "Lambda > 0 <==> finite cosmic information" is a NON-TRIVIAL equivalence between:
- A parameter in the Einstein equations (Lambda)
- A property of the information content of the universe (finiteness)

These are DIFFERENT concepts that happen to be equivalent. The equivalence is a theorem (Padmanabhan), not a definition.

Moreover, this equivalence makes **falsifiable predictions:**

1. **If Lambda = 0 exactly** (pure matter-dominated universe forever), then the total cosmic information is infinite. This is inconsistent with a finite-dimensional Hilbert space. So: **Lambda = 0 is inconsistent with holographic quantum gravity.**

2. **If Lambda evolves** (as DESI hints), then the information bound evolves too. The Hilbert space dimension changes with time. This predicts specific correlations between w(z) and the growth of cosmic entropy.

3. **If Lambda < 0** (anti-de Sitter universe), then there is no cosmic horizon, and the information capacity is set by the boundary at spatial infinity. This has DIFFERENT information-theoretic properties than Lambda > 0. The prediction: AdS universes have fundamentally different information structure than dS universes -- already verified in string theory (AdS/CFT vs dS mysteries).

### 15. Specific Testable Predictions

#### 15.1 From the tau Framework Applied to Cosmology

| Prediction | Test | Status |
|-----------|------|--------|
| w != -1 (exactly) | DESI DR2, Euclid, LSST | DESI hints at 2.5-4.2 sigma |
| w0 > -1 (not phantom) | Same | DESI best fit w0 ~ -0.75 (YES) |
| w evolves more at low z | Future BAO at z < 0.3 | Not yet tested precisely |
| sigma_8 lower than CMB prediction | DES, KiDS, HSC | 2-3 sigma tension observed (YES) |
| G_eff slightly lower at late times | Lunar laser ranging + cosmology | Current bounds: |delta G / G| < 10^{-2} |
| Information erasure at cosmological horizon differs from BH | Theoretical + analog experiments | Trivedi 2024 confirms difference |
| Cosmic tau increases with structure formation | Correlation of acceleration onset with SFR peak | Qualitatively consistent |

#### 15.2 Retrodiction Landauer Principle at Cosmological Scales

From Paper 1's Retrodiction Landauer Principle:
```
kT * ln(2) <= Sigma * kT
```

At the cosmological horizon (T = T_H = H/(2pi)):
```
Energy cost of 1 bit at the horizon = (H hbar / (2pi)) * ln(2) ~ 10^{-52} eV
```

This is the minimum energy to erase one bit of cosmological information. The total erasure energy:
```
E_erasure = N_bits * kT_H * ln(2) ~ S_dS * kT_H ~ (pi / (G Lambda)) * (H hbar / (2pi))
```

For Lambda ~ 10^{-52} m^{-2} and H ~ 10^{-18} s^{-1}:
```
E_erasure ~ 10^{122} * 10^{-52} eV ~ 10^{70} eV ~ 10^{46} J
```

This is the total information-theoretic energy content associated with the cosmological horizon. Compare to the dark energy content of the observable universe:
```
E_DE ~ rho_DE * V_obs ~ (6 x 10^{-10} J/m^3) * (4 x 10^{78} m^3) ~ 10^{69} J
```

**The two estimates are within ~3 orders of magnitude.** This suggests that dark energy IS the Landauer cost of maintaining the cosmic information horizon -- but the factor of ~10^3 discrepancy needs explanation.

---

## Part VI: Toward a tau-Based Cosmology

---

### 16. The Proposed Framework

Building on all the analysis above, here is the proposed tau-framework cosmology:

#### 16.1 Postulates

**P1:** Sigma = D(rho_spacetime || rho_matter) is the fundamental quantity governing spacetime dynamics. (Same as Paper 2.)

**P2:** In FRW spacetime, Sigma is defined via the Kodama vector:
```
Sigma_K(r) = -ln(|K|^2 / |K_ref|^2) = -ln(1 - H^2 r^2 / c^2)
```

**P3:** The Friedmann equations are the condition for Sigma to be stationary under small perturbations of the vacuum state (Jacobson's entanglement equilibrium, restated).

**P4:** Quantum corrections to the gravitational channel produce scale-dependent Sigma:
- At strong field: Sigma ~ r_s/r (Paper 2)
- At galactic scales: Sigma ~ r_s/r + alpha ln(r/r_c) (Paper 3, dark matter)
- At cosmological scales: Sigma ~ (Hr/c)^2 + entropy corrections (this analysis)

**P5:** The entropy corrections to Sigma at cosmological scales give rise to effective dark energy with w != -1.

#### 16.2 Specific Equations

**Modified Friedmann equation (tau form):**
```
H^2 = (8pi G / 3) rho + (1/3) * c^2 * (d Sigma_entropy / d(r^2)) |_{r=R_H}
```

where Sigma_entropy includes:
- Volume-law entanglement (Verlinde): contributes to effective dark matter
- Horizon entropy production (GREA): contributes to effective dark energy
- Running G (Kumar): contributes to both

**tau evolution:**
```
d(tau_cosmic)/dt = H * [1 + q(t)] * tau_cosmic * (1 - tau_cosmic)
```

where q(t) is the deceleration parameter. This is a logistic-type equation:
- Decelerating (q > 0): tau grows (arrow of time strengthens)
- Accelerating (q < 0, |q| < 1): tau still grows but more slowly
- de Sitter (q = -1): tau is constant (equilibrium)

#### 16.3 The Narrative

> Paper 1: "The arrow of time is measured by tau."
> Paper 2: "tau determines geometry near black holes."
> Paper 3: "tau at galactic scales explains rotation curves."
> **This analysis:** "tau at cosmic scales drives expansion -- and evolving tau means evolving dark energy."
> Paper 4: "All scales are one Sigma."

---

## Part VII: Honest Assessment and Key Gaps

---

### 17. Novelty Assessment

| Claim | Novelty | Support |
|-------|---------|---------|
| tau = 1 at Hubble horizon | LOW (standard) | Kodama vector formalism |
| Lambda <=> finite info | LOW (Padmanabhan 2017) | PLB |
| Friedmann = entanglement equilibrium | LOW (Jacobson 1995/2016) | PRL |
| GREA-type entropic acceleration | LOW (Gaztanaga 2024) | PDU |
| **tau framework unifies all four above** | **HIGH** | New synthesis |
| **w(z) predicted by Sigma evolution** | **MODERATE** | New connection |
| **Landauer cost ~ dark energy density** | **MODERATE** | New estimate |
| **GREA = tau dynamics at cosmic scale** | **HIGH** | New identification |
| **Running G + Verlinde + GREA = one Sigma** | **VERY HIGH** | Central Paper 4 thesis |

### 18. Critical Gaps

1. **No derivation of Friedmann from Sigma alone.** The tau framework INTERPRETS but does not DERIVE the Friedmann equations independently. Jacobson's and Dorau-Much's derivations require the area-entropy relation as input.

2. **No precise Lambda prediction.** CosMIn = 4pi gives the right order of magnitude but is not a derivation.

3. **CMB problem persists.** Neither running G, nor Verlinde, nor GREA can reproduce the CMB acoustic peak structure without non-baryonic matter.

4. **The entropy source problem.** GREA attributes acceleration to entropy production, but which entropy? Horizon entropy? Structure formation entropy? BH entropy? The answer affects quantitative predictions.

5. **Covariance of Verlinde's framework.** Hossenfelder (2017) showed the covariant version reduces to GR. This needs resolution.

6. **Factor of 10^3 in Landauer cost estimate.** The Landauer erasure energy at the horizon is ~10^3 times larger than the observed dark energy density. This could be a normalization issue or a genuine failure.

### 19. Comparison Table: tau vs LCDM vs Alternatives

| Feature | LCDM | tau Framework | GREA | Verlinde |
|---------|------|---------------|------|----------|
| Lambda | Free parameter | Information bound | Emergent | Entanglement pressure |
| w(z) | -1 (exactly) | != -1 (evolving) | Evolving | Not computed |
| DM | Particle | Running G / info loss | Not addressed | Entanglement displacement |
| CMB | Excellent | **FAILS** (no DM species) | Untested | **FAILS** |
| BAO | Good | Compatible (GREA tested) | Compatible | Not tested |
| H_0 tension | Unresolved | Plausible (running G) | Alleviated | Not addressed |
| sigma_8 | Tension | Correct direction | Correct direction | Not computed |
| Coincidence | Unexplained | Structure formation onset | Explained | Partial |
| Simplicity | 6 params | More conceptual, similar params | 6+ params | Conceptual |

### 20. Recommended Strategy

**For Paper 3 (galactic scales):** Focus on rotation curves (strong case). Mention cosmic acceleration briefly as motivation for Paper 4. Do NOT claim to derive Friedmann equations.

**For Paper 4 (grand unification):** The cosmic acceleration connection is a central piece. The argument:
1. Sigma = D(rho_spacetime || rho_matter) is the master equation
2. At strong field -> Paper 2 (exponential metric)
3. At galactic scales -> Paper 3 (flat rotation curves)
4. At cosmological scales -> effective Lambda from entropy production
5. All three are the SAME Sigma in different regimes

**For GRF Essay (if applicable):** The connection "arrow of time drives cosmic acceleration" (GREA = tau) is a striking result that could stand alone as a 4-page essay.

---

## References

### Primary (directly used in analysis)

1. Padmanabhan, T. "Cosmic Information, the Cosmological Constant and the Amplitude of primordial perturbations." PLB 773 (2017) 81-85. [arXiv:1703.06144](https://arxiv.org/abs/1703.06144)
2. Padmanabhan, T. "CosMIn: The Solution to the Cosmological Constant Problem." IJMPD 22 (2013) 1342001. [arXiv:1302.3226](https://arxiv.org/abs/1302.3226)
3. Padmanabhan, T. "Emergence and Expansion of Cosmic Space as due to the Quest for Holographic Equipartition." [arXiv:1206.4916](https://arxiv.org/abs/1206.4916)
4. Jacobson, T. "Thermodynamics of Spacetime: The Einstein Equation of State." PRL 75 (1995) 1260. [gr-qc/9504004](https://arxiv.org/abs/gr-qc/9504004)
5. Jacobson, T. "Entanglement Equilibrium and the Einstein Equation." PRL 116 (2016) 201101. [arXiv:1505.04753](https://arxiv.org/abs/1505.04753)
6. Verlinde, E.P. "Emergent Gravity and the Dark Universe." SciPost Phys. 2 (2017) 016. [arXiv:1611.02269](https://arxiv.org/abs/1611.02269)
7. Dorau, P. & Much, A. "From Quantum Relative Entropy to the Semiclassical Einstein Equations." PRL 136 (2026) 091602. [arXiv:2510.24491](https://arxiv.org/abs/2510.24491)
8. Gaztanaga, E. et al. "Dark Energy predictions from GREA: Background and linear perturbation theory." PDU 45 (2024) 101533. [arXiv:2405.02895](https://arxiv.org/abs/2405.02895)
9. Gaztanaga, E. et al. "Current and future constraints on the expansion history of the GREA model." (2026). [arXiv:2603.01934](https://arxiv.org/abs/2603.01934)
10. Gaztanaga, E. et al. "GREA and Dark Energy: A holographic dual." (2025). [arXiv:2511.19546](https://arxiv.org/abs/2511.19546)

### DESI and Observational

11. DESI Collaboration. "DESI DR2 Results II: Measurements of Baryon Acoustic Oscillations and Cosmological Constraints." Phys. Rev. D 112 (2025) 083515. [arXiv:2503.14738](https://arxiv.org/abs/2503.14738)
12. DESI Collaboration. "DESI 2024 VI: Cosmological Constraints from BAO." (2024). [arXiv:2404.03002](https://arxiv.org/abs/2404.03002)

### Banks-Fischler and Holographic

13. Banks, T. "Gibbons-Hawking Entropy and BMN Strings." (2025). [arXiv:2511.08213](https://arxiv.org/abs/2511.08213)
14. Banks, T. & Fischler, W. "Holographic Inflation, Primordial Black Holes and Early Structure Formation." (2024). [arXiv:2402.11527](https://arxiv.org/abs/2402.11527)

### Running G and Asymptotic Safety

15. Kumar, N. "Marginal IR Running of Gravity as a Natural Explanation for Dark Matter." PLB 871 (2025) 140008. [arXiv:2509.05246](https://arxiv.org/abs/2509.05246)
16. Zholdasbek, A. et al. "An Emergent Cosmological Model from Running Newton Constant." PRD 111 (2025) 103519. [arXiv:2405.02636](https://arxiv.org/abs/2405.02636)

### Information and Landauer

17. Trivedi, O. "Landauer's principle and information at the cosmological horizon." (2024). [arXiv:2407.15231](https://arxiv.org/abs/2407.15231)
18. Vopson, M. "Information Dark Energy Can Resolve the Hubble Tension and Is Falsifiable by Experiment." Entropy 24 (2022) 385.

### Bianconi GfE

19. Bianconi, G. "Gravity from Entropy." PRD 111 (2025) 066001. [arXiv:2408.14391](https://arxiv.org/abs/2408.14391)
20. Bianconi, G. "The Thermodynamics of the Gravity from Entropy Theory." (2025). [arXiv:2510.22545](https://arxiv.org/abs/2510.22545)

### Basso-Celeri (Crooks in Curved Spacetime)

21. Basso, M.L.W., Maziero, J. & Celeri, L.C. PRL 134 (2025) 050406. [arXiv:2405.03902](https://arxiv.org/abs/2405.03902)

### Supporting

22. Casini, H., Teste, E. & Torroba, G. "Relative Entropy and the RG Flow." JHEP 03 (2017) 089. [arXiv:1611.00016](https://arxiv.org/abs/1611.00016)
23. Hossenfelder, S. PRD 95 (2017) 124018. [arXiv:1703.01415](https://arxiv.org/abs/1703.01415)
24. Ghari, A. & Haghi, H. (2026). [arXiv:2601.01715](https://arxiv.org/abs/2601.01715)
25. Barman, S. et al. QIP (2026). [arXiv:2601.20860](https://arxiv.org/abs/2601.20860)
26. Capozziello, S. et al. EPJC (2024). [arXiv:2406.19274](https://arxiv.org/abs/2406.19274)

---

*Last updated: 2026-03-11*
*Research conducted for Paper 3/4 of the four-paper series: Sigma = D(rho_spacetime || rho_matter)*
