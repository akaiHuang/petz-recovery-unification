# How AeST Passes CMB Tests: Lessons for the Khronon Theory

**Author**: Sheng-Kai Huang (with research assistance)
**Date**: 2026-03-17
**Status**: Research analysis -- Key strategic findings
**Classification**: Critical comparison: AeST vs simplified Khronon

---

## Executive Summary

The simplified Khronon theory with K(Q) = mu^2(Q-1)^2 and running mu = H/c is excluded by Planck CMB at >10 sigma (17-20% excess dark matter at z=1100). Meanwhile, AeST (Skordis & Zlosnik 2021, PRL 127, 161302) explicitly passes the full CMB test. This document identifies the SPECIFIC mechanisms that allow AeST to succeed where our simplified model fails.

### The Answer in One Sentence

**AeST succeeds because the vector field A^mu provides SEPARATE degrees of freedom that decouple the cosmological background energy density from the MOND sector**, while our model tries to get BOTH from a single K(Q) function.

### Key Findings

| Question | Answer |
|----------|--------|
| Why does AeST pass CMB? | The vector field acts as "ersatz dark matter" at the background level, independently tunable |
| Why does our K(Q)-only model fail? | The single K(Q) cannot simultaneously give w=0 (CMB) and MOND (galaxies) |
| Does AeST predict Omega_DM? | **NO** -- it is an input parameter (fitted to 0.12) |
| Does AeST have c_s^2 = 0? | **YES** -- non-propagating mode with omega = 0 dispersion |
| How many extra parameters? | AeST: ~8-10 total (vs our 1-2) |
| Is tau compatible with AeST? | **YES** -- u^mu <-> A^mu is an exact structural match |
| Does Omega_DM = 0.268 survive? | Unclear -- the extremal principle may apply at a different level |
| What must we adopt? | The TWO-SECTOR architecture: separate background (CDM-like) from galactic (MOND) |

---

## 1. The Fundamental Difference: One Field vs Two Fields

### 1.1 Our Simplified Khronon Model

We use the Blanchet-Skordis (BS2024) action with only one dynamical scalar (the Khronon tau):

```
S = (c^3 / 16piG) integral sqrt(-g) [R - 2J(Y) + 2K(Q)] d^4x + S_matter
```

where:
- K(Q) = mu^2(Q-1)^2 controls the COSMOLOGICAL sector (dark matter energy density)
- J(Y) controls the GALACTIC sector (MOND rotation curves)
- Q = c * sqrt(-g^{mu nu} nabla_mu tau nabla_nu tau) (Khronon normalization)
- Y = A_mu A^mu / c^4 (acceleration-squared)

The Khronon field tau alone generates BOTH sectors. There is no independent vector field -- the normal n^mu = -(c/Q) g^{mu nu} nabla_nu tau is DERIVED from the scalar.

**The problem**: K(Q) must simultaneously:
1. Give rho_K ~ a^{-3} (CDM-like background for CMB)
2. Give c_s^2 = 0 (dust-like perturbations)
3. Connect to J(Y) at galactic scales (MOND)
4. Have K'(Q_0) = 0 at Q = 1 (ghost condensation / c_s^2 = 0)

Requirements 1 and 4 are incompatible: K'(1) = 0 means the minimum is at Q = 1, but the cosmological background has Q = 1 + delta with delta ~ 0.34, making the equation of state w = delta/(2+delta) = 0.145, which is far from the CDM value w = 0.

### 1.2 The Full AeST Theory

AeST (Skordis & Zlosnik 2021, arXiv:2007.00082) has THREE dynamical fields:

```
S = integral d^4x sqrt(-g) / (16pi G_tilde) *
    [ R
      - (K_B/2) F_{mu nu} F^{mu nu}
      + 2(2 - K_B) J^mu nabla_mu phi
      - (2 - K_B) Y
      - F(Y, Q)
      - lambda (A^mu A_mu + 1) ]
    + S_matter
```

where:
- **g_{mu nu}**: metric (as in GR)
- **A^mu**: unit timelike vector field (the "aether") -- INDEPENDENT dynamical field
- **phi**: shift-symmetric scalar field -- SEPARATE from the aether
- **lambda**: Lagrange multiplier enforcing A^mu A_mu = -1

The critical difference: A^mu and phi are INDEPENDENT fields. The scalar phi's time derivative Q = A^mu nabla_mu phi couples them, but they have separate dynamics and separate contributions to the stress-energy tensor.

### 1.3 The Key Structural Difference

| Feature | Our K(Q)-only model | AeST |
|---------|---------------------|------|
| Dynamical fields | metric + scalar tau | metric + vector A^mu + scalar phi |
| Independent DoF | 2 (tensor) + 1 (scalar) = 3 | 2 (tensor) + 2 (vector) + 2 (scalar) = 6 |
| Background DM source | K(Q) energy density from tau | Scalar phi "momentum" F_Q = C/a^3 |
| MOND source | J(Y) from tau's acceleration | J(Y) or F(Y,Q) with Y = spatial kinetic |
| Tunable background | NO -- locked by K(Q) form and delta_0 | YES -- F(Q) can be chosen to give exact w=0 |
| c_s^2 = 0 mechanism | K'(Q_0=1) = 0 (ghost condensation) | Lorentz violation: Y and Q independently controlled |
| Number of free parameters | 1-2 (mu, possibly Lambda_DBI) | 8-10 (K_B, Q_0, F(Y,Q) params, mu, lambda_s, ...) |

---

## 2. How AeST Passes the CMB: The Specific Mechanism

### 2.1 Background Level: Scalar Field as Perfect Dust

On the FLRW background:
- A^mu aligns with cosmic time: A^mu = (1/N, 0, 0, 0)
- phi depends only on time: phi = phi(t)
- Spatial gradients vanish: Y = 0
- Only Q = A^mu nabla_mu phi = phi_dot / N survives

The scalar field equation gives the CRITICAL conservation law:

```
F_Q = d F / d Q = C / a^3     (C = constant)
```

This is EXACTLY the dilution law for pressureless dust: the "momentum" F_Q scales as a^{-3}. The scalar field phi, through its kinetic coupling to the aether, behaves as a DUST-LIKE fluid.

The effective energy density and pressure of the scalar-aether system are:

```
rho_phi = -(1/(8pi G_tilde)) * (F - Q F_Q)
P_phi   =  (1/(8pi G_tilde)) * F
w_phi   = F / (Q F_Q - F)
```

**For CDM-like behavior**: need |F| << |Q F_Q|, so that rho_phi ~ Q F_Q / (8pi G_tilde) and P_phi ~ 0.

AeST achieves this by choosing F(Q) to have a MINIMUM near Q = Q_0 (the background value). Near the minimum, F ~ 0 while F_Q is determined by the conservation law. This makes w ~ 0 to high precision.

### 2.2 The Crucial Choice of F(Q)

In the original AeST paper (Skordis & Zlosnik 2021), the function F(Q) is chosen so that:

1. **At Q = Q_0 (cosmological background)**: F(Q_0) ~ 0, F_Q(Q_0) != 0. This gives rho_phi ~ Q_0 * C / (8pi G_tilde * a^3) -- exact a^{-3} scaling (dust).

2. **The equation of state**: w = F / (Q F_Q - F) ~ 0 because F ~ 0 at the background.

3. **The dust-like behavior persists for a limited cosmic time** (Dynamical systems analysis, arXiv:2309.06232). It is NOT exact dust forever, but it is close enough during the matter era and recombination to pass CMB tests.

Different choices of F(Q) give different cosmological histories:
- n = 2 models: approach dust-like behavior at LATE times
- n = 1/2 models: approximate dust at EARLY times

The flexibility to choose F(Q) is what gives AeST its CMB-fitting power.

### 2.3 Comparison with Our Model

Our model uses K(Q) = mu^2(Q-1)^2. On the FRW background:

```
rho_K = mu^2 * delta * (2+delta) / (8piG)     [BS2024 Eq. 4.14]
P_K   = mu^2 * delta^2 / (8piG)
w     = delta / (2+delta)                       [BS2024 Eq. 4.17]
```

**The problem**: w = delta/(2+delta) is NOT zero. For delta_0 = 0.34 (needed for Omega_DM = 0.265):
- w(z=0) = 0.145
- w(z=11) ~ 0.35 (peak, during matter era)
- w(z=1100) ~ 0.29

This is catastrophically far from CDM (w = 0). The energy density has a factor (2+delta) that varies by ~30% across the matter-radiation transition.

**In AeST**: The function F(Q) is CHOSEN to make w ~ 0. In our model, K(Q) = mu^2(Q-1)^2 is FIXED by the ghost condensation / Fisher metric argument, and it gives w >> 0.

### 2.4 The Vector Field's Role

The vector field A^mu in AeST plays THREE critical roles:

**(a) Providing the preferred frame**: A^mu defines the "cosmic time" direction. This is the same role as our Khronon normal n^mu, but A^mu is an INDEPENDENT dynamical field, not derived from a scalar.

**(b) Breaking Lorentz invariance**: The unit-norm constraint A^mu A_mu = -1 breaks Lorentz symmetry, allowing the temporal (Q) and spatial (Y) kinetic terms to be SEPARATELY controlled. This is ESSENTIAL for c_s = 0 (Section 3 below).

**(c) Acting as a "driving term" in perturbations**: At the perturbation level, the vector field provides an additional force that compensates baryon drag in the acoustic oscillations (Triton Station blog analysis). The vector field "acts like ersatz Dark Matter to closely reproduce the peaks in the CMB power spectrum." Where baryons dampen oscillations, the vector field maintains amplitude, enabling the third CMB peak to have the correct height.

**In our model**: The Khronon normal n^mu provides role (a) but is not independently dynamical. There is no separate vector degree of freedom to provide roles (b) and (c).

---

## 3. How AeST Achieves c_s^2 = 0: The Perturbation Mechanism

### 3.1 The Non-Propagating Mode

AeST has 6 physical degrees of freedom at the nonlinear level (Bataki, Skordis & Zlosnik 2023, arXiv:2307.15126):

| Mode | Count | Dispersion | Role |
|------|-------|------------|------|
| Tensor (GW) | 2 | omega = c*k (massless) | Standard gravitational waves |
| Vector | 2 | omega^2 = c^2*k^2 + m_V^2 (massive) | Constrained by GW observations |
| Scalar (massive) | 1 | omega^2 = c_s^2*k^2 + m_S^2 | Drives MOND phenomenology |
| Scalar (non-propagating) | 1 | **omega = 0** | **CDM-like mode** |

The CDM-like mode has omega = 0 -- it is NON-PROPAGATING. It grows under gravitational instability without oscillating, exactly like pressureless dust.

### 3.2 Why omega = 0?

The mechanism relies on the Lorentz-violating split between temporal (Q) and spatial (Y) kinetic terms:

1. The function F(Y, Q) is expanded around the cosmological background (Y=0, Q=Q_bg)
2. The TEMPORAL kinetic term comes from F_Q (non-zero -- gives healthy dynamics)
3. The SPATIAL kinetic term comes from F_Y (which can be made zero or very small)
4. Without a spatial kinetic term, perturbations cannot propagate -- their sound speed vanishes

**Theorem (informal)**: No Lorentz-invariant scalar field theory can have c_s = 0 perturbations around a homogeneous background. The aether A^mu provides the Lorentz violation needed to decouple Y from Q.

### 3.3 Our Model's c_s^2 = 0 Mechanism

Our model achieves c_s^2 = 0 through a DIFFERENT mechanism: ghost condensation. K(Q) = mu^2(Q-1)^2 has K'(Q=1) = 0, making the perturbation around Q = 1 a "stiff potential" rather than a kinetic term.

**This mechanism is analogous to k-essence**: In k-essence with L = f(X), c_s^2 = f'(X) / [f'(X) + 2X f''(X)]. At a minimum f'(X_0) = 0, we get c_s^2 = 0.

**The critical difference**: In our model, the cosmological background has Q = 1 + delta_0 = 1.34, NOT Q = 1. The perturbation is around Q_bg ~ 1.34, where K'(Q_bg) != 0. The c_s^2 = 0 property is exact only AT the ghost condensation point Q = 1.

BS2024 addresses this by noting that the perturbation sound speed is k-dependent:

```
c_s^2 = c_s^2(k)    [k-dependent, not just time-dependent]
```

For the GDM reduction, the Khronon produces:
- Zero bulk viscosity: c_vis^2 = 0
- Time-dependent equation of state: w(a) = delta(a) / (2 + delta(a))
- k-dependent sound speed: c_s^2(a, k)

**BS2024's key finding**: The quadratic K(Q) = mu^2(Q-1)^2 produces TENSION between CMB and MOND. The DBI form resolves this tension.

### 3.4 The DBI Resolution (BS2024)

BS2024 explicitly states that the quadratic K alone has "tension between the cosmology and MOND." The DBI kinetic function:

```
K_DBI(Q) = mu^2 * Lambda_D^2 * [sqrt(1 + (Q-1)^2/Lambda_D^2) - 1]
```

provides better behavior. BS2024 identifies two specific DBI parameter sets compatible with both MOND and CMB:

1. **mu^{-1} ~ 223 kpc, Lambda_D ~ 30** (transition at galactic scales)
2. **mu^{-1} ~ 22.3 Mpc, Lambda_D ~ 1** (transition at cosmological scales)

**However**: Our analysis (dbi_completion_2026_03_17.md) shows that the DBI completion trades one problem for another:
- In the deeply DBI regime, rho_DBI tracks rho_crit (dark energy-like) rather than a^{-3} (CDM-like)
- The DBI completion BREAKS c_s^2 = 0 away from the ghost condensation point
- The Gaussian Petz saturation (F = 1/Q) that gives Omega_DM = 0.268 is destroyed

**The DBI is NOT the answer for our model's CMB problem.** BS2024's CMB compatibility claim for DBI relies on the PERTURBATION-level behavior (c_s^2 near zero), not on the background equation of state being exactly w = 0.

---

## 4. What AeST Does That We Cannot Do with K(Q) Alone

### 4.1 The Two-Sector Architecture

AeST achieves CMB compatibility through a TWO-SECTOR architecture:

**Sector 1 -- Cosmological (background + perturbations):**
- The scalar field phi provides dust-like energy density rho ~ a^{-3}
- The function F(Q) is tuned so that w ~ 0 at all relevant epochs
- The non-propagating mode (omega = 0) mimics CDM perturbation growth
- The energy density Omega_DM h^2 ~ 0.12 is an INPUT (fitted to Planck)

**Sector 2 -- Galactic (quasi-static, weak-field):**
- The Y-dependent terms (spatial gradients of phi, acceleration of A^mu) become important
- The function F(Y, Q) or J(Y) interpolates to MOND at low accelerations
- This sector is completely independent of the cosmological background
- The MOND acceleration scale a_0 ~ 1.2 x 10^{-10} m/s^2 is empirical

**The key insight**: These two sectors are DECOUPLED in AeST. The cosmological parameters (Omega_DM, w, c_s^2) are controlled by F(Q) and can be freely adjusted without affecting the galactic MOND behavior (controlled by F(Y,Q)). Our model has only ONE function K(Q) controlling BOTH sectors, creating the fatal tension.

### 4.2 The "Dust-Like Behavior for a Limited Period" Issue

The dynamical systems analysis (arXiv:2309.06232) reveals an important subtlety: in AeST, "the dust-like contribution generally occurs for a limited period of cosmic time." The scalar field does NOT behave as perfect CDM forever -- it only approximates CDM during the matter era.

This means:
- AeST is NOT exactly CDM at all epochs
- There are small deviations at very early times (before dust-like phase) and late times (after)
- These deviations are small enough to be compatible with current CMB data
- They could potentially be detected by future experiments (CMB-S4, LiteBIRD)

**Implication for our model**: Even AeST only APPROXIMATELY mimics CDM. The question is whether the approximation is good enough (< 1% deviation during recombination). AeST achieves this; our model's 17-20% deviation is catastrophically too large.

### 4.3 The Price of AeST's Success

AeST pays for CMB compatibility with:

1. **More free parameters**: 8-10 total vs our 1-2. The function F(Y,Q) contains multiple tunable coefficients.

2. **The free function F(Y,Q) is ad hoc**: It is CHOSEN to interpolate between CDM cosmology and MOND galaxies. There is no first-principles derivation.

3. **Omega_DM is an input, not a prediction**: AeST does not explain WHY Omega_DM h^2 = 0.12. This value is set by the integration constant C in F_Q = C/a^3, which is fitted to observations.

4. **Galaxy cluster problems persist**: AeST has "negative mass density" artifacts at cluster scales (arXiv:2312.00889), and conflicting parameter requirements from galactic vs cluster scales (arXiv:2301.03499).

5. **No physical interpretation for the fields**: The aether A^mu and scalar phi are mathematical constructs without deeper physical meaning in AeST alone.

---

## 5. Compatibility of the tau Framework with AeST's Solution

### 5.1 The Structural Parallel Is Deep

The tau framework maps naturally onto AeST's field content:

| tau framework | AeST | Match quality |
|---------------|------|---------------|
| Observer u^mu | Aether A^mu | **Exact structural match** |
| u^mu u_mu = -1 | A^mu A_mu = -1 | **Identical constraint** |
| Sigma = -ln(-g_{00}) | Related to F(Y,Q) | **Strong parallel** |
| Petz recovery failure = decoherence | Non-propagating mode = constraint | **Conceptual match** |
| "Dark matter = cost of time" | "Dark matter = scalar field momentum" | **Complementary interpretations** |

The identification u^mu <-> A^mu <-> n^mu (Khronon normal) has PHYSICAL content: it says the preferred foliation that defines time IS the field that generates dark-matter-like energy density.

### 5.2 What tau Could Add to AeST

The tau framework could potentially provide what AeST currently lacks:

**(a) Prediction of Omega_DM h^2**:
AeST fits Omega_DM to observations. The tau framework's extremal principle (maximize R = integral F * Sigma dt) could predict the scalar field integration constant C, thereby predicting Omega_DM from first principles. This would reduce AeST's parameter count by one.

**Status**: The current prediction Omega_DM = 0.268 (1.1% from observed) uses the simplified K(Q)-only model that fails CMB. Whether a similar extremal principle works in AeST is an open question.

**(b) Derivation of a_0 = cH_0/(2pi)**:
AeST treats a_0 as an empirical input. The tau framework derives it from the de Sitter temperature T_dS = hbar*H_0 / (2pi*k_B), connecting the MOND scale to cosmological QRE.

**(c) Physical interpretation of the aether**:
AeST introduces A^mu as a mathematical field. The tau framework identifies it as the observer's worldline -- "time itself." This provides a physical interpretation: the aether IS the direction of time flow, and dark matter IS the energetic cost of maintaining temporal order.

**(d) Derivation of F(Y,Q) from QRE**:
The dream derivation: show that S_grav = integral sqrt(-g) D(rho_spacetime || rho_matter) produces the AeST action when rho_spacetime includes foliation information. This would make AeST a CONSEQUENCE of quantum information geometry rather than an ad hoc construction.

**Status**: 0/4 steps of this derivation completed. This is conjecture, not result.

### 5.3 Tension Points

**(a) Strong-field regime**:
AeST admits "stealth" black hole solutions with EXACT Schwarzschild metric (arXiv:2412.15395). The tau framework predicts exponential metric g_00 = -exp(-r_s/r) with NO event horizons (Paper 2). These are incompatible. Resolution requires understanding which theory applies where.

**(b) Number of fields**:
The tau framework's elegance comes from using a single Sigma = 2 ln Q to unify everything. AeST's success comes from having MULTIPLE fields with SEPARATE roles. Adopting AeST means accepting that "one equation for everything" may need to become "one framework organizing multiple equations."

**(c) Parameter count**:
The tau framework aims for minimal parameters. AeST has 8-10. The tau framework can potentially REDUCE AeST's parameters (by predicting Omega_DM, a_0, and constraining F(Y,Q)), but it cannot eliminate the need for the two-sector architecture.

---

## 6. What We Need to Adopt from AeST

### 6.1 Non-Negotiable Changes

Based on this analysis, the following modifications to our framework are NECESSARY:

**(1) Separate the background DM from the galactic MOND**:
The single K(Q) CANNOT provide both CDM-like background and MOND phenomenology. We need either:
- The full AeST action with independent vector and scalar fields, OR
- A modified Khronon theory where the background and perturbation contributions are architecturally separated

**(2) Accept that Omega_DM is (currently) a fitted parameter**:
Until the extremal principle is reformulated for the AeST context, Omega_DM h^2 = 0.12 should be treated as an input, not a prediction. The 0.268 prediction from the simplified model is associated with a cosmological sector that is excluded by CMB.

**(3) The equation of state must be w ~ 0 to better than 1% at recombination**:
Any viable theory must satisfy the Planck constraint on the dark matter equation of state. Our K(Q) = mu^2(Q-1)^2 gives w = 0.29 at z = 1100, which is excluded at >10 sigma. The replacement must give |w| < 0.01 at recombination.

### 6.2 Negotiable Changes (Strategic Choices)

**(4) Whether to use the full AeST or a simpler subset**:
The BS Khronon theory (BS2024) is simpler than AeST (one scalar vs scalar + vector), and BS2024 claims CMB compatibility with DBI kinetic terms. However, our analysis shows the DBI route has serious issues (loss of c_s^2 = 0, tracking dark energy). The full AeST with its 6 DoF may be necessary.

**(5) Whether the tau framework provides the organizing principle or the mechanism**:
Option A: tau derives AeST from QRE principles (dream scenario, no calculation exists)
Option B: tau provides conceptual interpretation of AeST (u^mu = aether = time, DM = temporal cost)
Option C: tau and AeST are complementary but independent theories

Current evidence supports Option B as the realistic near-term position.

**(6) Whether to keep the Omega_DM prediction as a "suggestive result"**:
The 0.268 prediction is beautiful and may point to deep physics, even if the specific derivation path (quadratic K + running mu) is wrong. It could be presented as a suggestive result motivating further work, rather than a firm prediction.

---

## 7. Does the Omega_DM = 0.268 Prediction Survive?

### 7.1 The Prediction's Logical Structure

The prediction follows from:
```
Maximize R = integral [F(Q) * Sigma_SE(Q)] dt
           = integral [(Q^2-1)/Q] dt

=> delta_0 = 0.343, Omega_DM = delta_0*(2+delta_0)/3 = 0.268
```

This requires:
1. K(Q) = mu^2(Q-1)^2 (quadratic, for Sigma_SE = Q^2 - 1)
2. F = 1/Q (Gaussian Petz saturation)
3. mu = H/c (running)
4. The extremal principle itself

### 7.2 Impact of AeST Adoption

If we adopt AeST's two-sector architecture:

- Items 1-2 apply ONLY if the Khronon scalar still has K(Q) = mu^2(Q-1)^2. In AeST, the function is F(Q), which has a different form. The extremal principle would need to be reformulated for AeST's F(Q).

- Item 3 (running mu) is conceptually similar to AeST's scalar field evolution. The conservation law F_Q = C/a^3 in AeST is analogous to our K' = I_0/a^3.

- Item 4 (the extremal principle) could potentially apply to AeST's scalar sector. The key question: does maximizing some QRE-related functional of F(Q) give the correct integration constant C that matches Omega_DM = 0.265?

### 7.3 Possible Survival Routes

**(a) Direct reformulation**: Define a new extremal principle for AeST's F(Q) that reduces to the Omega_DM = 0.268 result. This would require:
- Identifying the QRE analog of F * Sigma_SE in AeST
- Showing that the extremum gives C = specific value
- Checking that this value corresponds to Omega_DM ~ 0.265

**(b) Universality argument**: The Omega_DM = 0.268 result might be UNIVERSAL -- independent of the specific K(Q) or F(Q) form, depending only on the extremal principle structure. If so, it would survive the transition to AeST.

**(c) Coincidence**: The 1.1% accuracy could be a numerical coincidence. The prediction lies in the range 0.2-0.4 where many O(1) expressions naturally land.

### 7.4 Honest Assessment

**Probability of survival in current form: 15%**
**Probability of survival in modified form: 30%**
**Probability of coincidence: 55%**

The prediction is suggestive but cannot be considered robust until reformulated in an AeST-compatible context.

---

## 8. Comparison Table: Our Model vs AeST vs What We Need

| Feature | Our K(Q) model | AeST (S&Z 2021) | What we need |
|---------|---------------|-----------------|--------------|
| Background w at z=1100 | 0.29 (EXCLUDED) | ~0 (SAFE) | |w| < 0.01 |
| Background rho scaling | a^{-3} * (2+delta) (varies ~30%) | a^{-3} (exact to ~1%) | a^{-3} to <1% |
| c_s^2 at perturbation | 0 at Q=1, k-dependent elsewhere | 0 (non-propagating mode) | 0 or < 10^{-6} |
| c_T (GW speed) | c (exact) | c (exact) | c (exact) -- BOTH pass |
| MOND at galaxy scales | J(Y) sector (independent, works) | F(Y,Q) (works) | YES |
| Galaxy clusters | Unknown | Problematic (negative rho) | Major open problem |
| Omega_DM | Predicted: 0.268 | Input: 0.12 (fitted) | Predicted OR well-motivated input |
| a_0 | Derived: cH_0/(2pi) | Input: empirical | Derived OR well-motivated |
| Free parameters | 1-2 | 8-10 | Minimal (tau reduces AeST?) |
| Physical interpretation | "DM = cost of time" | None (fields are ad hoc) | tau provides interpretation |
| CMB fit | EXCLUDED (>10 sigma) | PASSES (Planck compatible) | Must pass |
| BBN | SAFE (40 OoM margin) | SAFE | SAFE |
| GW170817 | SAFE (structural) | SAFE (structural) | SAFE |

---

## 9. Strategic Conclusions

### 9.1 The Core Lesson

**A single K(Q) function CANNOT simultaneously produce:**
1. CDM-like background evolution (w = 0, rho ~ a^{-3})
2. Ghost condensation (K'(Q=1) = 0) giving c_s^2 = 0
3. Connection to MOND at galactic scales
4. The correct Omega_DM

AeST solves this by having MULTIPLE fields with SEPARATE roles. The tau framework must accept this architectural lesson while preserving its conceptual contributions (physical interpretation, potential parameter predictions, QRE connection).

### 9.2 The Path Forward

**Immediate** (for Paper 3):
- Split Paper 3 into galactic sector (solid, publishable) and cosmological sector (needs work)
- The galactic sector (J(Y) + MOND + RAR) is INDEPENDENT of the CMB problem
- State explicitly: "The cosmological completion requires additional structure, as demonstrated by the AeST framework"

**Short-term** (1-3 months):
- Study whether the tau extremal principle can be reformulated for AeST's F(Q)
- Investigate whether the Khronon-Tensor extension (BS 2025, arXiv:2507.00912) resolves any issues
- Explore whether the tau framework can PREDICT AeST's free function F(Y,Q) from DPI constraints

**Long-term** (3-12 months):
- Full integration of tau framework with AeST
- Potential reduction of AeST's parameter count through QRE-derived constraints
- CLASS/Boltzmann implementation of the combined theory
- Confrontation with Planck, DESI, CMB-S4 data

### 9.3 What Survives Regardless

The following components of the tau framework are INDEPENDENT of the CMB problem and survive regardless of the cosmological sector's resolution:

1. **tau = 1 - F** (Paper 1): Pure QIT, no gravity assumptions
2. **F >= exp(-Sigma/2)** (Paper 1): Petz bound, proven
3. **Sigma = 2 ln Q on adiabatic backgrounds** (Paper 2): Valid for any K(Q), independent of CMB
4. **Exponential metric in strong field** (Paper 2): Independent of cosmological K(Q)
5. **c_T = c** (all papers): Structural from the BS action
6. **c_s^2 = 0 at ghost condensation** (Paper 3): Local property at Q = 1
7. **J(Y) sector giving MOND** (Paper 3): Independent of K(Q) cosmology
8. **u^mu <-> A^mu identification**: Provides physical meaning to the aether

---

## 10. Technical Appendix: Why the f(delta) Problem Is Structural

### 10.1 The Mathematical Proof

For K(Q) = mu^2(Q-1)^2 with running mu = H(a)/c:

```
rho_K = mu^2 * delta * (2+delta) / (8piG)
      = [H^2(a)/c^2] * delta(a) * [2+delta(a)] / (8piG)
```

The conservation law a^3 * Q * K' = I_0 gives:
```
a^3 * (1+delta) * 2mu^2 * delta = I_0
=> mu^2 * delta = I_0 / [2 a^3 * (1+delta)]
```

Therefore:
```
rho_K = I_0 * (2+delta) / [16piG * a^3 * (1+delta)]
      = (I_0 / 16piG) * [(2+delta)/(1+delta)] * a^{-3}
```

For exact CDM: rho = rho_0 * a^{-3}. The ratio:
```
r(z) = rho_K(z) / [rho_K(0) * (1+z)^3] = [(2+delta(z))/(1+delta(z))] / [(2+delta_0)/(1+delta_0)]
```

Since delta(z) varies from 0.34 (today) to ~1.08 (matter era) to ~0.82 (z=1100):
- r(z=0) = 1 (by definition)
- r(z=11) = 1.32 (matter era peak)
- r(z=1100) = 1.205 (recombination)

**This 20.5% excess at recombination is the fundamental problem.** It requires (2+delta)/(1+delta) to be constant across all redshifts, which is impossible unless delta = const (requires constant mu, but that gives a^{-9} stiff matter).

### 10.2 Why AeST Avoids This

In AeST, the effective energy density is:
```
rho_phi = (Q * F_Q - F) / (8pi G_tilde)
```

The conservation law F_Q = C/a^3 gives:
```
rho_phi = [Q * C/a^3 - F(Q)] / (8pi G_tilde)
```

If F(Q_bg) ~ 0 (minimum), then:
```
rho_phi ~ Q_bg * C / (8pi G_tilde * a^3)
```

This is EXACT a^{-3} scaling (assuming Q_bg varies slowly). There is no (2+delta) factor because F(Q) is not quadratic but is chosen to have F ~ 0 at the background.

**The structural difference**: In our model, K(Q=1+delta) = mu^2*delta^2 is always NONZERO at the background (delta=0.34 means K is O(mu^2)). In AeST, F(Q_bg) is TUNED to be ~ 0 at the background. This tuning is what buys CDM-like behavior -- at the cost of having a free function rather than a derivation.

---

## References

### Primary
1. Skordis, C. & Zlosnik, T. (2021). "New relativistic theory for Modified Newtonian Dynamics." PRL 127, 161302. arXiv:2007.00082.
2. Blanchet, L. & Skordis, C. (2024). "Relativistic Khronon Theory..." JCAP 11, 040. arXiv:2404.06584.
3. Skordis, C. & Zlosnik, T. (2021). "AeST: Linear stability on Minkowski space." PRD 106, 104041. arXiv:2109.13287.

### AeST Analysis
4. Bataki, Skordis & Zlosnik (2023). "AeST: Hamiltonian Formalism." arXiv:2307.15126.
5. Mistele, McGaugh & Hossenfelder (2023). "AeST: Quasistatic spherical solutions." MNRAS 531, 272. arXiv:2305.07742.
6. Dynamical systems analysis (2023). arXiv:2309.06232.
7. Galaxy cluster problems (2023). arXiv:2312.00889.
8. Weak lensing confrontation (2023). arXiv:2301.03499.
9. Skordis & Vokrouhlicky (2024). "Stealth black holes in AeST." arXiv:2412.15395.

### Khronon Extensions
10. Blanchet & Skordis (2025). "Khronon-Tensor theory." arXiv:2507.00912.

### Crisis Analysis (This Repository)
11. fdelta_crisis_analysis_2026_03_17.md -- BS formula correction and CMB crisis proof
12. dbi_completion_2026_03_17.md -- DBI completion analysis (FAILS)
13. honest_assessment_2026_03_17.md -- Strategic assessment
14. paper4_AeST_mechanism.md -- Detailed AeST technical analysis

---

*Last updated: 2026-03-17*
*This document identifies the specific mechanisms by which AeST passes CMB tests and our simplified Khronon model fails. The core lesson: a single K(Q) cannot do both CDM cosmology and MOND galaxies. AeST's two-sector architecture (independent vector + scalar) is the proven solution. The tau framework's role shifts from "mechanism" to "organizing principle" -- providing physical interpretation, potential parameter predictions, and QRE-derived constraints for the richer AeST structure.*
