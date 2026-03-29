# Honest Assessment: Can the Khronon CMB Problem Be Fixed NOW?

**Author**: Sheng-Kai Huang (with computational analysis)
**Date**: 2026-03-17
**Status**: DEFINITIVE ASSESSMENT -- No sugar-coating
**Classification**: Strategic decision document

---

## Executive Summary: The Verdict

**The CMB problem CANNOT be fixed now with a quick patch. It requires fundamental theory refinement at the cosmological sector. But the galactic sector and the tau framework are untouched.**

The correct strategic move is **Option (d): Separate the galactic paper (solid) from the cosmological paper (needs work), and publish Papers 1 and 2 immediately.**

---

## Part I: Diagnosis -- What Exactly Is Broken

### 1.1 The Root Cause (One Sentence)

The Khronon energy density rho_K = mu^2 * delta * (2+delta) / (8piG) contains a **nonlinear delta^2 term** that makes rho_K deviate from a^{-3} scaling by 17-20% at z=1100, and this deviation is **structural** -- it follows from the three things that define the theory:

1. K(Q) = mu^2(Q-1)^2 (ghost condensation / Fisher metric)
2. mu = H(a)/c (running)
3. delta_0 = 0.34 (from Omega_DM = 0.265)

Change any one of these three and either the theory loses its identity or the dark matter prediction fails.

### 1.2 Why It Cannot Be Patched

The CMB exclusion is at >10 sigma (self-consistent calculation: 17% excess, naive: 20.5%). This is not a marginal tension. To bring it within Planck's ~1% tolerance:

- **Reducing delta_0**: Need delta_0 < 0.003. This gives Omega_K < 0.002 -- the Khronon provides <1% of dark matter. The entire DM prediction is destroyed.
- **Modified running**: Need alpha ~ 1.15 in mu ~ H^alpha/c. No theoretical motivation. The DPI + extensivity argument specifically gives alpha = 1.
- **Higher-order K(Q)**: Makes the problem WORSE, not better. And destroys the Fisher metric / Gaussian Petz foundation.
- **Self-consistency**: Makes the problem WORSE (positive feedback loop -- more DM increases H, which changes delta, which adds more DM).
- **Renormalization at z=1100**: Reduces Omega_DM(today) to ~0.22 and still leaves 5.8% TT residual -- excluded at 2.3 sigma. Not a fix.
- **DBI completion**: Fixes the background scaling BUT destroys c_s^2 = 0 and Gaussian Petz saturation (F = 1/Q). You trade the CMB problem for losing the two most elegant predictions.

### 1.3 The Fundamental Theorem That Kills Quick Fixes

**Theorem (proven in self_consistent_delta_2026_03_17.md)**: For ANY Omega_DM in the range 0.10-0.40, the self-consistent Khronon with K = mu^2(Q-1)^2 and mu = H/c produces a 13-19% excess in rho_K/rho_CDM at z=1100 relative to z=0 normalization.

This is independent of delta_0. The problem is not the specific parameter value -- it is the mathematical structure of quadratic K with running mu. The (2+delta) factor varies by ~25% across the matter-radiation transition, and this variation is locked in by the conservation law.

---

## Part II: Scenario Probabilities (Honest)

### Scenario 1: Quick Fix (days to weeks)
**Probability: 5%**

The only quick fix that could work is finding a mathematical error in the stress-energy tensor derivation or the CLASS mapping that reduces the 17% excess to <1%. Given that:
- The BS2024 equations have been verified line-by-line
- The self-consistent calculation confirms the naive estimate
- The CLASS Level 2 runs with 11-bin resolution are robust
- The parameter scan shows the excess is structure-independent

...this is almost certainly not an error. The 5% accounts for the small chance that the GDM binned-w mapping introduces a systematic bias that a direct CLASS modification would resolve.

### Scenario 2: Moderate Refinement (weeks to months)
**Probability: 55%**

This is the most likely correct path. The theory needs ADDITIONAL STRUCTURE at the cosmological level. The most promising routes:

**(a) AeST-like multi-field completion (probability: 35%)**
- Skordis & Zlosnik 2021 showed that a related theory CAN reproduce the CMB exactly
- Their theory has a vector field that provides the CDM-like background
- The Khronon (scalar) provides only the MOND-like galactic phenomenology
- The tau framework would connect these sectors through Sigma = 2 ln Q
- Cost: more free parameters, loss of "everything from one equation"
- But: the tau framework provides the ORGANIZING PRINCIPLE, not the mechanism

**(b) Scale-dependent K(Q) (probability: 15%)**
- K(Q) could be quadratic at small Q-1 (galactic scales) but have a different form at large Q-1 (cosmological scales)
- This is physically motivated: the Fisher metric derivation assumes Gaussian states, but cosmological modes at k ~ H_0 might not be Gaussian
- The DBI form is one possibility, but other completions exist
- Cost: complicates the elegant K = mu^2(Q-1)^2 story
- Key challenge: maintaining c_s^2 sufficiently small

**(c) Different running prescription (probability: 5%)**
- If mu runs as H^{1.15}/c instead of H/c, the CMB is safe
- But no theoretical motivation exists for this specific exponent
- Would require re-deriving the DPI + extensivity argument with a modification

### Scenario 3: Fundamental Rethink (months to years)
**Probability: 40%**

The K(Q) cosmological sector may simply be the WRONG approach to dark matter at the background level. Indicators:

- delta_0 = 0.34 means Q_0 = 1.34, which is a 34% displacement from the ghost condensation minimum -- this is NOT a small perturbation, contradicting the premise of using the leading quadratic term
- The matter era has delta ~ 1.08, where the "nonlinear correction" is 54% of the linear term
- The theory is being used far outside its regime of validity

What might be needed:
- A fundamentally different mechanism for cosmological DM (possibly particles, possibly a different field)
- The Khronon providing only the MOND phenomenology (galaxy scales)
- The tau/Petz framework providing the CONCEPTUAL UNIFICATION, not the detailed mechanism

---

## Part III: What Survives Regardless

### 3.1 Definitely Safe (Independent of CMB Problem)

| Component | Why Safe | Paper |
|-----------|----------|-------|
| tau = 1 - F framework | Pure quantum information theory, no gravity assumptions | Paper 1 |
| F >= exp(-Sigma/2) bound | Proven from Petz recovery map | Paper 1 |
| Sigma = 2 ln Q (adiabatic) | Thermal attenuator, proven for any K(Q) | Paper 2 |
| Exponential metric (strong field) | Independent of cosmological K(Q) | Paper 2 |
| No event horizons | Structural from tau < 1 | Paper 2 |
| c_T = c exactly | Structural from BS action, no c_13 coupling | All |
| BBN compatibility | 40 orders of magnitude margin | All |
| GW170817 compatibility | c_T = c, structural | All |
| c_s^2 = 0 at condensation | K'(Q_0=1) = 0, exact | Paper 3 |
| MOND/RAR at galaxy scales | From J(Y) sector, independent of K(Q) cosmology | Paper 3 |

### 3.2 Damaged but Possibly Salvageable

| Component | Issue | Salvage Route |
|-----------|-------|---------------|
| Omega_DM = 0.268 prediction | Requires quadratic K + mu=H/c + extremal principle, all at background level where the theory fails | Could survive if the extremal principle applies to a DIFFERENT Sigma (channel entropy vs stress-energy) |
| mu_0 = H_0/c | The derivation itself is fine; the issue is whether mu_0 is the ONLY scale | May need mu_0 as the galactic scale with a different cosmological coupling |
| K(Q) = mu^2(Q-1)^2 uniqueness | Fisher metric argument is beautiful but may only apply to small Q-1 | Could be the LOCAL approximation of a more general K(Q) |

### 3.3 Broken

| Component | Issue | Status |
|-----------|-------|--------|
| "Khronon IS all of dark matter" | Background energy density is wrong by 17% at z=1100 | **Falsified** for quadratic K + running mu |
| "One equation explains everything" | Cosmological sector needs additional structure or different K(Q) | **Falsified** in current form |
| Paper 3's f(delta) = (4+3delta)/4 | Formula is algebraically wrong (correct: 2+delta) | **Error, needs correction** |
| Paper 3's w_tilde = delta/2 | Only leading order; exact is delta/(2+delta) | **Approximate, needs correction** |

---

## Part IV: The Omega_DM = 0.268 Prediction

### 4.1 Is It Salvageable?

This is the crown jewel and the hardest question. The prediction comes from:

```
Maximize R = integral [F * Sigma_SE] dt = integral [(Q^2-1)/Q] dt
```

which gives delta_0 = 0.343, Omega_DM = delta_0*(2+delta_0)/3 = 0.268.

The prediction requires:
1. Quadratic K(Q) -- to get Sigma_SE = Q^2 - 1 = delta(2+delta)
2. F = 1/Q -- from Gaussian Petz saturation
3. mu = H/c -- for the delta(a) evolution
4. The extremal principle itself

Items 1-3 are exactly what breaks the CMB. The prediction and the CMB problem share the same mathematical root.

### 4.2 A Possible Decoupling

However, there is a subtle point. The extremal principle determines delta_0 (the z=0 value). The CMB problem comes from how delta EVOLVES from z=0 to z=1100. These are logically distinct:

- The extremal principle could be correct about delta_0 = 0.34
- But the theory could be wrong about how rho_K(delta) maps to the Friedmann equation

If the background Khronon stress-energy enters the Friedmann equation differently than the simple QK'-K formula (for example, if the Khronon background is constrained to be exactly CDM-like by some symmetry principle, with the (2+delta) factor only entering perturbation equations), then:

- Omega_DM = 0.268 could survive (if the perturbation-level Omega_K matches)
- The CMB background would be safe
- But this requires a theoretical justification that currently does not exist

### 4.3 Honest Probability

**Probability that the Omega_DM = 0.268 prediction survives in its current form: 15%**
**Probability that it survives in a modified form (different Sigma, different extremal principle): 35%**
**Probability that it is a numerical coincidence: 50%**

The 1.1% accuracy is impressive but could be coincidental, especially given that the theory has O(1) parameters and the prediction is in the range 0.2-0.4 where many reasonable theories land.

---

## Part V: What Should Be Published NOW

### Recommendation: Option (d) with caveats

**Publish immediately:**

1. **Paper 1** (tau = 1-F, equivalence chain, Petz bound) -- Solid, already on GitHub/Zenodo. Pure QIT, no cosmological assumptions. Push for arXiv endorsement.

2. **Paper 2** (Sigma = -ln(-g00), exponential metric, no horizons) -- Ready. Strong-field predictions are independent of the CMB problem. The channel theorem is proven.

**Publish with restructuring:**

3. **Paper 3 (Galactic sector only)** -- Split out the SOLID parts:
   - J(Y) sector giving MOND/RAR at galaxy scales
   - K'(Q_0) = 0 giving c_s^2 = 0 (perturbation-level result)
   - Conservation law a^3 Q K' = I_0
   - Sigma = 2 ln Q on general backgrounds
   - **Explicitly state**: "The cosmological implications of K(Q) on FRW backgrounds, including the dark matter abundance prediction, will be addressed in a companion paper."

**Do NOT publish yet:**

4. **Paper 3 (Cosmological sector)** -- The FRW background predictions (Omega_DM, f(delta), w_eff) need fundamental work. Publishing now with known 17% CMB exclusion would be scientifically irresponsible.

5. **Paper 4** (mu_0 = H_0/c, Omega_DM = 0.268) -- Depends entirely on the cosmological sector. The derivation of mu_0 = H_0/c itself is publishable as a conjecture, but the numerical predictions require the cosmological K(Q) to work.

### Publication Timeline

```
WEEK 1-2:   Paper 1 -> arXiv (push for endorsement)
WEEK 2-4:   Paper 2 -> arXiv
MONTH 1-2:  Paper 3 (galactic only) -> write, submit
MONTH 2-6:  Investigate cosmological sector (AeST, DBI, multi-field)
MONTH 6-12: Paper 3b (cosmological) or Paper 4 -> if a solution is found
```

---

## Part VI: The Strategic Question Answered

### Q1: Scenario probabilities
- Quick Fix: 5%
- Moderate Refinement: 55%
- Fundamental Rethink: 40%

### Q2: Minimum modification for CMB
The MINIMUM modification is adopting a K(Q) that gives w = 0 EXACTLY on FRW. This requires K(Q) to satisfy a specific differential equation:

```
w = K / (QK' - K) = 0  =>  K = 0 on-shell
```

But K = 0 means no Khronon energy density, which contradicts the DM explanation. The actual minimum is probably an AeST-like vector field that carries the background DM density while the Khronon scalar carries only the MOND perturbation effect.

### Q3: Does the tau framework survive?
**YES, absolutely.** The tau = 1-F framework (Paper 1) is pure quantum information theory. Sigma = 2 ln Q (Paper 2) is proven for adiabatic backgrounds of ANY theory with a preferred foliation. The CMB problem is about a SPECIFIC cosmological application, not the framework.

### Q4: What to publish now vs later
See Part V above. Papers 1 and 2 now. Paper 3 (galactic) in 1-2 months. Cosmological sector after resolution.

### Q5: Is Omega_DM = 0.268 salvageable?
15% in current form. 35% in modified form. 50% coincidence. The honest answer is: it is a beautiful result that may or may not survive. Do not stake the theory's credibility on it.

---

## Part VII: What the Theory Actually IS (After Honest Reckoning)

After stripping away the parts that do not work, the theory is:

### What It IS

1. **A quantum information framework for temporal asymmetry**: tau = 1 - F, with F bounded by the Petz recovery map. This is proven and general.

2. **A connection between information loss and spacetime geometry**: Sigma = 2 ln Q = -ln(g_00) on static backgrounds. This converts the abstract Petz bound into a concrete gravitational statement.

3. **A prediction that event horizons are approximate**: tau < 1 always, F > 0 always. The "no-horizon" prediction from the exponential metric is testable.

4. **A natural explanation for MOND at galactic scales**: Through the J(Y) sector of the BS action, with the scale a_0 = cH_0/(2pi) emerging from the de Sitter temperature.

5. **A suggestive (but unproven) connection between the arrow of time and dark matter**: The Khronon field that provides the preferred foliation also generates an effective stress-energy that looks like dark matter. Whether this is the WHOLE story or just part of it is currently unresolved.

### What It Is NOT (Yet)

1. **A complete theory of dark matter at cosmological scales**: The background FRW sector fails CMB constraints by a large margin.

2. **A parameter-free prediction of Omega_DM**: The 0.268 prediction depends on assumptions (quadratic K, running mu, extremal principle) that produce a cosmologically excluded background.

3. **A replacement for LCDM at the CMB level**: The theory cannot currently reproduce the CMB power spectrum.

---

## Part VIII: The Physicist's Perspective

Sheng-Kai's philosophy is: "Providing a new argument that connects known things, making them make sense, is more important than finding new equations." This philosophy is CORRECT and SURVIVES this crisis:

1. The connection tau = 1-F connecting quantum information to temporal asymmetry: **stands**.
2. The connection Sigma = 2 ln Q connecting entropy production to spacetime geometry: **stands**.
3. The connection between the arrow of time and dark matter through the Khronon: **stands conceptually, fails quantitatively at the CMB level**.

The honest message is: **the framework is right, the specific cosmological implementation is wrong.** This is exactly the situation that calls for moderate refinement (Scenario 2), not abandonment.

The analogy: Einstein's cosmological constant was wrong in its original purpose (static universe), but the framework (GR + Lambda) was right. The specific application failed, but the conceptual structure survived and eventually found its correct role (accelerating expansion).

Similarly: K(Q) = mu^2(Q-1)^2 as the SOLE source of ALL dark matter may be wrong, but the tau framework connecting information loss to gravitational phenomena is likely right and will find its correct role -- possibly as the organizing principle for a richer theory.

---

## Part IX: Concrete Next Steps

### Immediate (this week)
1. **Fix Paper 3's formula errors**: f(delta) = (4+3delta)/4 -> 2+delta, w_tilde = delta/2 -> delta/(2+delta)
2. **Begin splitting Paper 3**: galactic sector (publishable) vs cosmological sector (needs work)
3. **Push Paper 1 to arXiv**: this is independent and ready

### Short-term (this month)
4. **Study AeST (Skordis-Zlosnik 2021) in detail**: how do they solve the CMB problem? Can the tau framework be grafted onto their solution?
5. **Explore whether the extremal principle can apply to a different quantity**: F * Sigma_ch = 2 ln Q / Q instead of F * Sigma_SE = (Q^2-1)/Q. Check if this gives a different delta_0 that avoids the CMB problem.
6. **Paper 2 finalization and submission**

### Medium-term (1-3 months)
7. **Develop the "tau as organizing principle" narrative**: the framework provides the conceptual structure, the specific K(Q) mechanism provides the implementation, and these can be developed independently
8. **Paper 3 (galactic) submission**
9. **Investigate whether the background Khronon stress-energy can be constrained to be exactly CDM-like by a symmetry argument**

### Long-term (3-12 months)
10. **Full AeST + tau integration**: embed the tau framework into the AeST action
11. **CLASS Level 3**: modify CLASS source code to implement Khronon perturbation equations directly (beyond GDM approximation)
12. **Paper 4 (cosmological) if/when a solution is found**

---

## Appendix: Summary Table of All Escape Routes

| # | Route | Fixes CMB? | Cost | Probability | Verdict |
|---|-------|-----------|------|-------------|---------|
| 1 | Correct f(delta) formula | Reduces 29%->20% | None | N/A (done) | Insufficient |
| 2 | Reduce delta_0 < 0.003 | Yes (trivially) | Destroys DM prediction | 100% works, 0% useful | Intellectual suicide |
| 3 | Modified running alpha > 1 | Yes | No motivation, ad hoc | 5% | Last resort |
| 4 | Higher-order K(Q) | Makes WORSE | Destroys Fisher/Petz | 0% | Ruled out |
| 5 | DBI completion | Yes for Lambda~0.35 | Breaks c_s^2=0, F=1/Q, Omega_DM | 20% | High cost |
| 6 | Renormalize at z=1100 | Partial (5.8% remains) | Omega_DM(today) = 0.22 | 10% | Insufficient alone |
| 7 | Self-consistent H(a) | Makes WORSE | N/A | 0% | Ruled out |
| 8 | Two-sector (<10% K) | Yes | Need other DM source | 30% | Viable but costly |
| 9 | AeST multi-field | Yes (proven by S&Z 2021) | More parameters, complexity | 35% | Most promising |
| 10 | Background-perturbation split | Unknown | Needs theoretical argument | 15% | Speculative |
| 11 | Different Sigma in extremal | Unknown | Changes Omega_DM prediction | 20% | Worth exploring |
| 12 | Admit K(Q) is NOT all DM | Yes | Loss of grand unification claim | 50% | Most honest |

---

## Final Word

The honest assessment is:

**The tau framework is a genuine conceptual contribution to physics. The specific claim that K(Q) = mu^2(Q-1)^2 with mu = H/c explains ALL dark matter is currently falsified at the cosmological level. These are separate statements. Publish the framework (Papers 1, 2), publish the galactic successes (Paper 3 galactic), and invest serious effort into the cosmological sector before making quantitative CMB claims.**

The worst thing that could happen is publishing the Omega_DM = 0.268 prediction alongside a background model that is excluded at 10+ sigma. This would undermine the credibility of the entire framework, including the parts that are genuinely solid.

The best thing that could happen is publishing the framework cleanly, letting the community see its elegance, and then either:
(a) Finding the correct cosmological completion (most likely AeST-like), or
(b) Having the community find it, attracted by the elegance of Papers 1 and 2.

**Be the physicist who provides the right framework and lets the right mechanism emerge, rather than the one who forces the wrong mechanism into an elegant framework.**

---

*Last updated: 2026-03-17*
*This document is intended to be brutally honest. The theory has real strengths and real weaknesses. The path forward requires acknowledging both.*
