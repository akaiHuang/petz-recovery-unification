# J(Y) Functional Form Analysis: Rigorous Assessment

**Author**: Sheng-Kai Huang (with computational analysis)
**Date**: 2026-03-16
**Status**: Complete analysis with honest parameter counting
**Sources**: Blanchet-Skordis 2024 (arXiv:2404.06584), Blanchet-Skordis 2025 (arXiv:2507.00912), Skordis-Zlosnik 2021 (arXiv:2007.00082, arXiv:2109.13287), Skordis-Vokrouhlicky 2024 (arXiv:2412.15395)

---

## Executive Summary

**J(Y) is a FREE FUNCTION in Blanchet-Skordis 2024.** It is never specified as a unique closed-form expression. Only its asymptotic behavior at low Y (MOND regime) and high Y (GR regime) is constrained. This is the single most important fact for honest parameter counting: the BS Khronon theory contains an entire free function, not just a few parameters.

The full theory has: **1 free function (J(Y)) + 3 genuine free parameters (a_0, mu, Omega_DM h^2) + 1 semi-free function (K(Q), constrained to have a quadratic minimum)**, plus the standard 6 cosmological parameters. The J and K sectors are structurally independent -- there is no known theoretical mechanism connecting them.

---

## 1. Is J(Y) Specified in BS 2024?

### Answer: NO -- J(Y) is a free function

From Blanchet-Skordis 2024 (arXiv:2404.06584, JCAP 11, 040):

The action is (their Eq. 2.6):
```
S = (c^3 / 16piG) integral d^4x sqrt(-g) [R - 2 J(Y) + 2 K(Q)] + S_m[Psi, g]
```

The paper explicitly states that J(Y) "has already been postulated in [18] and [19] and is part of the original BM theory." References [18] and [19] are Blanchet & Marsat (2011) and Blanchet (2021). The section heading reads "Choice of function J(Y)" -- the word "choice" indicating it is not uniquely determined.

**BS 2024 provides ONLY the deep-MOND asymptotic expansion (their Eq. 3.23):**
```
J(Y) ~ Lambda - Y + (c^2 / a_0) Y^{3/2} + O(Y^2)     when Y << a_0^2/c^4
```

And the high-acceleration limit:
```
J(Y) -> 0      when Y >> a_0^2/c^4     (recovers GR)
```

**No closed-form J(Y) valid for all Y is given.** The paper focuses on the K(Q) sector (which is new to BS 2024) and takes J(Y) as inherited from the BM/MOND tradition.

### What BS 2025 (Khronon-Tensor) adds

Blanchet-Skordis 2025 (arXiv:2507.00912) confirms the same structure:
- High acceleration: J = 0 when Y >> a_0^2/c^4
- Low acceleration: J = Lambda - Y + (2c^2/3a_0) Y^{3/2} + O(Y^2) when Y << a_0^2/c^4

Note the slightly different numerical coefficient (2/3 vs 1) between 2024 and 2025 versions -- this reflects different conventions/normalizations, not physics.

### Historical context

The function J(Y) traces back to the original Bekenstein-Milgrom (1984) AQUAL theory, where the MOND Lagrangian density contains a free function f(|nabla Phi|^2/a_0^2). In the relativistic setting of BS Khronon, this becomes J(Y). **The free function has been free since 1984.** Forty-two years of MOND research have not uniquely determined it.

---

## 2. Constraints on J(Y)

### 2.1 Deep-MOND limit (Y -> 0, low acceleration)

For Y << a_0^2/c^4:
```
J(Y) = Lambda - Y + (c^2/a_0) Y^{3/2} + O(Y^2)
```

This gives the MOND interpolating function behavior:
- J_Y = dJ/dY = -1 + (3c^2/(2a_0)) Y^{1/2} + ...
- In the modified Poisson equation: Delta phi + nabla . (J_Y nabla Xi) + mu^2 Xi = 4piG rho_m
- At low acceleration: J_Y -> -1, so the equation becomes: nabla . [(1-1) nabla Xi] + ... which reduces to the MOND equation with sqrt(a_0 a_N) behavior

**This constraint fixes the O(Y^{3/2}) coefficient in terms of a_0.**

### 2.2 Newtonian/GR limit (Y -> infinity, high acceleration)

For Y >> a_0^2/c^4:
```
J(Y) -> 0    (or J -> Lambda for the cosmological constant piece)
J_Y -> 0
```

When J_Y = 0, the modified Poisson equation reduces to:
```
Delta phi + mu^2 Xi = 4piG rho_m
```

At scales r << 1/mu (solar system, stellar systems), the mu^2 Xi term is negligible compared to Delta phi, and standard Newtonian gravity is recovered.

**This constraint requires J to asymptotically vanish faster than Y at large Y.**

### 2.3 PPN constraints (Solar system)

BS 2024 shows that in the quasi-static weak-field limit, the theory gives:
```
phi = psi     (equal gravitational potentials in Newtonian gauge)
```

This means gamma_{PPN} = 1 exactly, regardless of the choice of J(Y), as long as the high-Y limit is satisfied. The PPN constraint does NOT further constrain the functional form of J(Y).

### 2.4 Transition region

The region Y ~ a_0^2/c^4 (transition between MOND and Newtonian) is where J(Y) is least constrained. This corresponds to the MOND "interpolating function" mu(x), where x = a/a_0:
```
mu(x) ~ 1 + J_Y(Y(x))
```

Different choices of J(Y) in the transition region correspond to different interpolating functions. Common choices in the MOND literature:
- "Simple": mu(x) = x/(1+x)
- "Standard": mu(x) = x/sqrt(1+x^2)
- "RAR": mu(x) = 1/(1 - exp(-sqrt(x)))

**These are phenomenologically distinguishable** with high-quality rotation curve data (e.g., SPARC), but current data does not uniquely select one.

### 2.5 Are the constraints enough to determine J(Y)?

**NO.** The asymptotic constraints at Y -> 0 and Y -> infinity are necessary conditions but leave the transition region completely free. Any smooth function interpolating between the two limits is viable. This is equivalent to the well-known freedom in choosing the MOND interpolating function mu(x).

---

## 3. Complete Free Parameter/Function Count

### 3.1 The BS Khronon action

```
S = (c^3/16piG) integral d^4x sqrt(-g) [R - 2J(Y) + 2K(Q)] + S_m
```

### 3.2 Itemized freedom

| Item | Type | Status | Source |
|------|------|--------|--------|
| **J(Y)** | Free function | Constrained at Y -> 0 and Y -> infinity only | MOND tradition since 1984 |
| **a_0** | Parameter of J | Input, a_0 ~ 1.2 x 10^{-10} m/s^2 | Empirical (rotation curves) |
| **Lambda** | Parameter of J | Fixed by cosmological constant measurement | Observation |
| **K(Q)** | Function | Semi-fixed: K(Q) = mu^2(Q-1)^2 + higher-order | Ghost condensation argument |
| **mu** | Parameter of K | Free; mu^{-1} sets the Khronon Compton wavelength | Dimensional analysis suggests mu_0 = H_0/c |
| **Q_0** | Background value | Set by cosmological initial conditions | Not predicted |
| **Omega_DM h^2** | Cosmological | Input (initial condition of Khronon field) | 0.12 from Planck; not predicted |
| **G** | Coupling | Newton's constant, measured | Cavendish experiment |
| **6 standard** | Cosmological | H_0, Omega_b h^2, n_s, A_s, tau_reion, etc. | Standard cosmology |

### 3.3 Honest count

**Beyond standard cosmological parameters (6), the BS Khronon theory has:**

1. **1 free function**: J(Y) -- the MOND sector. Contains effectively infinite freedom, though in practice parameterized by a_0 plus the shape of the interpolating function (1 parameter + 1 shape function)
2. **1 free parameter**: mu (or equivalently, the Khronon mass scale)
3. **1 free initial condition**: Omega_DM h^2 (the Khronon field amplitude on the cosmological background)
4. **1 semi-free function**: K(Q), constrained to be quadratic near Q = 1 but with freedom in higher-order terms (DBI extension adds lambda_D)

**Minimal count**: 6 (standard) + a_0 + mu + Omega_DM h^2 + shape of J = **9 parameters + 1 free function**

**Comparison with LCDM**: 6 standard + Omega_DM h^2 = **7 parameters, no free functions**

**Comparison with MOND**: a_0 + shape of mu(x) = **1 parameter + 1 free function** (but MOND has no cosmology)

**Comparison with AeST (SZ 2021)**: 6 standard + K_B + K_2 + lambda_s + Q_0 + a_0 + shape of F(Y,Q) = **~10 parameters + 1 free function**

### 3.4 The DBI extension of K(Q)

BS 2024 Section 4.3 discusses generalizations of K(Q) beyond the quadratic form:
- Cubic: K(Q) = mu^2(Q-1)^2 + alpha_3(Q-1)^3
- DBI-inspired: K(Q) = mu^2 * f_DBI((Q-1)/epsilon) where f_DBI is a DBI-type function
- The DBI form introduces at least one additional parameter (lambda_D or epsilon)

BS 2025 explicitly mentions the DBI form as preferred for passing "both constraints from the CMB and MOND for reasonable values of parameters." This suggests the quadratic K(Q) alone may not be sufficient for full observational compatibility, and the DBI extension (with its additional parameter) may be necessary.

---

## 4. Connection Between J(Y) and a_0

### 4.1 How a_0 enters

a_0 appears explicitly in the deep-MOND expansion of J(Y):
```
J(Y) ~ Lambda - Y + (c^2/a_0) Y^{3/2} + ...
```

It sets the transition scale between MOND and Newtonian regimes:
- Y ~ a_0^2/c^4: transition region
- Y >> a_0^2/c^4: Newtonian
- Y << a_0^2/c^4: deep MOND

### 4.2 Is a_0 emergent or input?

**a_0 is an INPUT parameter of J(Y).** It is not derived from the theory. It is set by fitting to galaxy rotation curve data (empirically a_0 ~ 1.2 x 10^{-10} m/s^2).

### 4.3 The tantalizing coincidence a_0 ~ cH_0/(2pi)

Numerically:
```
cH_0/(2pi) = (3 x 10^8)(2.27 x 10^{-18})/(2pi) = 1.08 x 10^{-10} m/s^2
```

This is within 10% of the empirical a_0 ~ 1.2 x 10^{-10} m/s^2. This "Milgrom coincidence" has been known since the 1980s but has never been rigorously derived.

**In the tau framework**, the argument is:
1. mu_0 = H_0/c (dimensional analysis)
2. a_0 = c^2 mu_0 / (2pi) = cH_0/(2pi) (if a_0 and mu are connected)

But as established in Section 5 below, a_0 and mu live in INDEPENDENT sectors of the action. The connection a_0 ~ mu c^2 is a suggestive numerological relation, not a theoretical derivation.

### 4.4 Could a_0 emerge from K(Q)?

No. a_0 controls the Y-dependence of J, while mu controls the Q-dependence of K. The variables Y and Q are:
- Y = A_mu A^mu / c^4 (spatial acceleration squared, relevant in quasi-static systems)
- Q = c sqrt(-g^{mu nu} partial_mu tau partial_nu tau) (temporal gradient of Khronon, relevant in cosmology)

In the quasi-static limit, Q -> 1 and K -> 0, so K(Q) contributes nothing to galactic dynamics. In the cosmological limit, Y -> 0 and J -> Lambda, so J contributes nothing to cosmological dark matter density. **The two sectors decouple.**

---

## 5. Can J(Y) Be Connected to K(Q)?

### 5.1 Structural analysis

The action decomposes cleanly:
```
S = S_GR + S_MOND[J(Y)] + S_cosmo[K(Q)] + S_matter
```

where:
- S_MOND = -(c^3/8piG) integral sqrt(-g) J(Y) d^4x controls galactic phenomenology
- S_cosmo = (c^3/8piG) integral sqrt(-g) K(Q) d^4x controls cosmological DM

The variables Y and Q have different physical origins:
- Y measures the acceleration of the Khronon congruence (spatial, non-perturbative)
- Q measures the temporal gradient of the Khronon (temporal, perturbative around Q = 1)

### 5.2 Where could a connection arise?

**Possibility A: Through the Khronon field itself**

Both J and K depend on the same scalar field tau:
- Q = c sqrt(-g^{mu nu} partial_mu tau partial_nu tau)
- A_mu = -c^2 q_mu^nu nabla_nu ln Q -> Y = A_mu A^mu / c^4

So Y is actually constructed from Q through spatial gradients of ln Q. In principle, a nonlinear field equation could couple the two sectors. However, in the regimes where each is important (galactic for J, cosmological for K), the other sector's contribution is negligible. The coupling is formally present but practically irrelevant.

**Possibility B: Through a unified function**

One could imagine a single function G(Y, Q) that reduces to J(Y) when Q -> 1 and to K(Q) when Y -> 0:
```
G(Y, Q) -> J(Y)     when Q = 1 (quasi-static)
G(Y, Q) -> K(Q)     when Y = 0 (cosmological)
```

This is essentially what AeST's F(Y, Q) does. However, even in AeST, the Y-dependence and Q-dependence are effectively independent -- no constraint from the theory forces them to be related.

**Possibility C: Through the tau framework's DPI**

The Data Processing Inequality constrains both sectors via monotonicity of quantum relative entropy. But DPI only gives inequality constraints (eta >= 0, c_s^2 >= 0), not equalities connecting J to K.

### 5.3 Verdict

**J(Y) and K(Q) are genuinely independent sectors.** There is no known theoretical principle -- within the BS Khronon theory, AeST, or the tau framework -- that connects them. The numerical coincidence a_0 ~ c^2 mu_0 is suggestive but not derived.

This is arguably the **most significant theoretical gap** in the BS Khronon framework: two sectors, each with its own freedom, unified only by sharing the same Khronon field but not by any deeper principle connecting their functional forms.

---

## 6. Comparison with RMOND/AeST (Skordis-Zlosnik 2021)

### 6.1 AeST action

From Skordis-Zlosnik 2021 (arXiv:2007.00082, PRL 127, 161302) and Skordis-Vokrouhlicky 2024 (arXiv:2412.15395):

```
S_g = integral d^4x sqrt(-g)/(16pi G_tilde) [R - (K_B/2) F_{mu nu} F^{mu nu}
      + 2(2-K_B) J^mu nabla_mu phi - (2-K_B) Y - F(Y, Q) - lambda(A^mu A_mu + 1)]
```

where:
- A^mu: unit timelike vector field (aether)
- phi: shift-symmetric scalar field
- Y = D^mu phi D_mu phi (spatial kinetic term)
- Q = A^mu nabla_mu phi (temporal derivative along aether)
- F(Y,Q): free function controlling both MOND and cosmological behavior

### 6.2 The free function F(Y,Q) in AeST

**F(Y,Q) is ALSO a free function in AeST.** Skordis & Zlosnik give a specific example for their CMB fit (in arXiv:2109.13287, the extended version), but emphasize that the function is not uniquely determined.

For the strong-field/black-hole analysis, Skordis & Vokrouhlicky (2024) use only the leading terms:
```
F = (2-K_B) lambda_s Y - 2 K_2 (Q - Q_0)^2 + ...
```
explicitly excluding MOND terms (Y^{3/2}) since those are irrelevant at compact-object scales.

### 6.3 AeST free parameters

From the Skordis-Vokrouhlicky analysis:

| Parameter | Role | Constraint |
|-----------|------|------------|
| K_B | Aether kinetic coupling | 0 < K_B < 2 (stability) |
| lambda_s | Spatial scalar coupling | lambda_s >= 0 (stability) |
| K_2 | Temporal scalar coupling | K_2 > 0 (ghost-free) |
| Q_0 | Background Q value | Set by cosmological solution |
| F(Y,Q) | Full free function | MOND + CMB matching |
| a_0 | MOND scale | Parameter of F |

**Total AeST-specific parameters: K_B + lambda_s + K_2 + Q_0 + a_0 = 5 numerical parameters + 1 free function F(Y,Q)**

### 6.4 BS Khronon vs AeST: comparison

| Feature | BS Khronon (2024) | AeST (SZ 2021) |
|---------|-------------------|-----------------|
| Fields | metric + scalar tau | metric + vector A^mu + scalar phi |
| Free function | J(Y) (MOND only) | F(Y,Q) (MOND + cosmo combined) |
| Cosmo kinetic | K(Q) = mu^2(Q-1)^2 (specific) | Part of F(Y,Q) (within free function) |
| c_T = c | Yes (by construction) | Yes (by construction) |
| c_s = 0 mode | Yes (from K(Q) minimum) | Yes (from F minimum in Q) |
| Free parameters | a_0, mu (+higher-order) | K_B, lambda_s, K_2, Q_0, a_0 |
| Degrees of freedom | 3 (2 tensor + 1 scalar) | 6 (2 tensor + 2 vector + 2 scalar) |

**Key difference**: BS Khronon is simpler (fewer fields, fewer parameters) but separates the MOND sector (J) from the cosmological sector (K) more explicitly. AeST lumps both into a single free function F(Y,Q) with more internal parameters.

### 6.5 Did SZ 2021 specify J(Y) for CMB?

SZ 2021 demonstrated CMB + matter P(k) fits with a specific choice of F(Y,Q), showing the framework CAN work. However:
1. The specific F used for the CMB fit was given in the extended paper (arXiv:2109.13287), not in the PRL letter
2. The choice is an EXAMPLE, not claimed to be unique
3. The CMB fit is primarily sensitive to the Q-sector (cosmological perturbations with c_s = 0), not the Y-sector (MOND)
4. The MOND behavior enters only through galactic-scale predictions, which are separate from the CMB fit

**Therefore: SZ 2021's CMB fit does NOT uniquely determine J(Y) / the MOND function.** It only constrains the cosmological sector (analogous to K(Q) in BS language).

---

## 7. What the tau Framework Contributes (and What It Does Not)

### 7.1 What the tau framework adds

1. **Physical interpretation of the Khronon**: tau = observer's time direction, not just a mathematical field. "Dark matter" = price of temporal ordering.

2. **Petz optimality argument for K(Q)**: The kinetic function should have a quadratic minimum (c_s = 0) because this is the maximally recoverable channel configuration. This provides a principle selecting K(Q) = mu^2(Q-1)^2 rather than arbitrary functions.

3. **Dimensional analysis for mu_0 = H_0/c**: The unique choice without Planck-scale physics.

4. **Crooks fluctuation theorem for mu(x)**: A proposed MOND interpolating function mu(x) = 1 - exp(-x) from the Crooks relation P(Sigma)/P(-Sigma) = exp(Sigma). This would constrain J(Y) if validated.

### 7.2 What the tau framework does NOT add

1. **Does not derive J(Y)**: The free function remains free. The Crooks interpolation is a proposal, not a derivation.

2. **Does not connect J to K**: No principle from DPI, JRSWW, or Petz connects the MOND sector to the cosmological sector.

3. **Does not predict a_0**: The relation a_0 = cH_0/(2pi) is numerology/naturality, not a theorem.

4. **Does not predict Omega_DM h^2**: This remains an initial condition.

5. **Does not reduce the free function count**: BS Khronon with tau framework still has 1 free function (J) + ~3 free parameters.

---

## 8. Summary Table: What Is Determined vs What Is Free

| Quantity | Status | Determined by | Confidence |
|----------|--------|---------------|------------|
| K(Q) form: quadratic minimum | DETERMINED | Ghost condensation + Petz optimality | Strong |
| K(Q) = mu^2(Q-1)^2 specifically | ASSUMED | Leading-order truncation | Moderate |
| mu_0 = H_0/c | DETERMINED | Dimensional analysis (unique) | Strong |
| mu(k) = k (running) | ASSUMED | DPI + extensivity + modular flow | Moderate (5 arguments converge) |
| J(Y) at Y -> 0 | CONSTRAINED | MOND phenomenology | Strong (observations) |
| J(Y) at Y -> infinity | CONSTRAINED | GR recovery | Strong (solar system tests) |
| J(Y) transition region | FREE | Not determined by any principle | N/A |
| a_0 value | INPUT | Rotation curve fits | Empirical |
| a_0 ~ cH_0/(2pi) relation | SUGGESTIVE | Numerology / naturality | Weak |
| J(Y) connected to K(Q) | NO | No known principle | N/A |
| Omega_DM h^2 = 0.12 | INPUT | CMB + BAO fit | Not predicted |
| Lambda (cosmological constant) | INPUT | Supernova + CMB | Not predicted |

---

## 9. Honest Assessment for Paper Writing

### What can be claimed

1. "The BS Khronon theory unifies MOND and cosmological dark matter in a single action with one scalar field."
2. "The cosmological sector K(Q) = mu^2(Q-1)^2 is fixed by ghost condensation / Petz optimality, with one parameter mu_0 = H_0/c from dimensional analysis."
3. "The theory has fewer free parameters than AeST (no K_B, lambda_s, K_2) while maintaining the same observational reach."
4. "The tau framework provides a physical interpretation: the Khronon IS the observer's time direction."

### What cannot be claimed

1. ~~"One equation, one parameter"~~ -- The theory has 1 free function + ~3 free parameters beyond standard cosmology.
2. ~~"J(Y) is derived from information theory"~~ -- It is free, inherited from MOND tradition.
3. ~~"a_0 and mu are connected"~~ -- They live in independent sectors.
4. ~~"The theory is more predictive than LCDM"~~ -- It has MORE freedom (a free function) than LCDM.
5. ~~"The CMB fit uniquely determines the MOND function"~~ -- CMB constrains K, not J.

### Recommended framing

"The BS Khronon theory contains two independent functional sectors: J(Y) controlling galactic MOND dynamics, and K(Q) controlling cosmological dark matter density. The K(Q) sector is well-constrained (quadratic form from ghost condensation, mu_0 = H_0/c from dimensional analysis, running mu(k) = k from DPI + extensivity). The J(Y) sector inherits the traditional MOND free function, constrained at both asymptotic limits but free in the transition region. The tau framework contributes a physical interpretation (Khronon = observer's time) and constrains K(Q) but does not determine J(Y). The numerical coincidence a_0 ~ c^2 mu_0/(2pi) suggests a deeper connection between the two sectors, but no derivation exists."

---

## 10. Comparison: Parameter Counts Across Theories

| Theory | Standard cosmo params | Theory-specific params | Free functions | Total |
|--------|----------------------|----------------------|----------------|-------|
| LCDM | 6 | 1 (Omega_DM h^2) | 0 | 7 |
| MOND (non-relativistic) | N/A | 1 (a_0) | 1 (mu(x)) | N/A |
| TeVeS (Bekenstein 2004) | 6 | ~3-4 | 1 (f(Y)) | ~10 + 1 func |
| AeST (SZ 2021) | 6 | ~5 (K_B, K_2, lambda_s, Q_0, a_0) | 1 (F(Y,Q)) | ~11 + 1 func |
| BS Khronon (2024) | 6 | ~3 (a_0, mu, Omega_DM h^2) | 1 (J(Y)) | ~9 + 1 func |
| BS Khronon + tau | 6 | ~2 (a_0, Omega_DM h^2; mu fixed) | 1 (J(Y)) | ~8 + 1 func |
| BS Khronon + tau + Crooks | 6 | ~2 (a_0, Omega_DM h^2) | 0 (J constrained!) | ~8 |

**Key observation**: The Crooks interpolation function mu(x) = 1 - exp(-x), if validated against SPARC data, would eliminate the free function entirely. This is the single most impactful prediction the tau framework could make: reducing "1 free function + 2 parameters" to "0 free functions + 2 parameters."

---

## 11. Open Questions for Future Work

1. **Can J(Y) be derived from QRE?** If D(rho_spacetime || rho_matter) generates both J and K from a single variational principle, this would be a genuine derivation.

2. **Does the Crooks interpolation fit SPARC?** A systematic fit of mu(x) = 1 - exp(-x) to the SPARC 175-galaxy sample would test the tau framework's most specific prediction for J(Y).

3. **Can a_0 be derived from mu_0?** The coincidence a_0 ~ c^2 mu_0/(2pi) needs either a proof or an explanation for why it is accidental.

4. **What is the DBI form of K(Q)?** BS 2024 discusses DBI extensions but does not give the specific DBI K(Q). This introduces at least one more parameter (lambda_D or epsilon).

5. **Can the K-J coupling be measured?** In the full nonlinear theory, J and K do couple through the Khronon field equations. Are there regimes where this coupling is observable?

---

## References

- Blanchet & Skordis 2024: arXiv:2404.06584, JCAP 11 (2024) 040 -- "Khronon theory: the original Blanchet-Marsat nonrelativistic modified gravity as a relativistic theory"
- Blanchet & Skordis 2025: arXiv:2507.00912 -- "Khronon-Tensor theory"
- Skordis & Zlosnik 2021a: arXiv:2007.00082, PRL 127, 161302 -- AeST (CMB + MOND)
- Skordis & Zlosnik 2021b: arXiv:2109.13287 -- AeST extended version
- Skordis & Vokrouhlicky 2024: arXiv:2412.15395 -- Stealth black holes in AeST
- Bekenstein & Milgrom 1984: ApJ 286, 7 -- AQUAL (original MOND free function)
- Scherrer 2004: PRL 93, 011301 -- K-essence as unified dark matter
- Arkani-Hamed et al. 2004: JHEP 0405:074 -- Ghost condensation
- Thomas, Kopp & Skordis 2016: arXiv:1601.05097 -- GDM constraints from Planck
- Bataki, Skordis & Zlosnik 2023: arXiv:2307.15126 -- AeST Hamiltonian (6 DoF)
