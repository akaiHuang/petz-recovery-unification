# f(delta) Crisis: Systematic Escape Route Analysis

**Author**: Sheng-Kai Huang (with computational analysis)
**Date**: 2026-03-17
**Status**: Complete analysis of 12 escape routes
**Severity**: CRITICAL -- the f(delta) background correction is excluded at ~13 sigma by Planck

---

## 0. Executive Summary

The Khronon dark matter theory predicts that the dark matter energy density deviates from pure a^{-3} dilution due to the f(delta) background correction. With delta_0 = 0.34 (needed for Omega_DM = 0.265), the correction produces ~20-29% more DM at z=1100 than CDM, which is catastrophically excluded by Planck.

**Key Discovery**: The Paper 3 formula f(delta) = (4+3delta)/4 is INCONSISTENT with the standard k-essence stress-energy tensor, which gives g(delta) = (2+delta). This reduces the excess from 29% to 20%, but 20% is still fatally excluded.

### Viability Ranking

| Rank | Route | Description | Verdict | Action Required |
|------|-------|-------------|---------|-----------------|
| 1 | **Formula correction** | f(delta)=(4+3d)/4 -> g(delta)=2+d | Reduces crisis, doesn't solve | Verify against BS2024 Eq. 3.8 |
| 2 | **Re-normalize + H_0** | Fit omega_K at z=1100, not z=0 | Partial fix + H_0 tension help | Re-run CLASS with corrected formula |
| 3 | **Perturbation-only** | Background is EXACTLY CDM; f only affects perturbations | Best if justifiable | Needs theoretical argument |
| 4 | **Two-sector** | K(Q) provides <10% of DM | Safe if K fraction < 10% | Needs J(Y) to provide cosmological DM |
| 5 | **Small delta_0** | delta_0 << 0.34 | Safe but no DM prediction | Destroys the Omega_DM prediction |
| 6 | **DBI completion** | K_DBI linearizes at large delta | Promising but changes conservation law | Full DBI cosmology needed |
| 7 | **Modified running** | mu ~ H^alpha/c^alpha, alpha > 1 | Works for alpha ~1.25 | No theoretical motivation |
| 8 | **K(Q) = mu^2(Q-1)^n** | Higher n reduces excess | Loses Omega_DM and Petz saturation | Unacceptable cost |
| 9 | **Self-consistent H(a)** | Backreaction of rho_K on H | Makes crisis WORSE | Ruled out |

---

## 1. The Problem in Detail

### 1.1 The f(delta) Correction

With running mu = H(a)/c, the Khronon condensate parameter evolves as:

```
delta(a) = delta_0 * (H_0/H(a))^2 / a^3
```

This gives:
- Radiation era: delta ~ 0 (safe)
- Matter era: delta ~ delta_0/Omega_m ~ 1.08 (problematic)
- Today: delta_0 = 0.34

The Khronon energy density is NOT exactly proportional to a^{-3}. It includes a correction factor:

```
rho_K(a) propto g(delta(a)) / a^3
```

where g(delta) encodes the non-trivial stress-energy of the K(Q) sector.

### 1.2 Two Competing Formulas

**Paper 3 formula** (citing BS2024):
```
rho_K = c^2 mu^2 delta f(delta) / (8piG)
f(delta) = (4 + 3*delta) / 4
```

**K-essence stress-energy tensor** (standard derivation):
```
rho = 2X L_X - L  where  L = 2K(Q),  X = Q^2/2
rho_K = [c^4/(8piG)] * [Q K'(Q) - K(Q)]
     = [c^2 mu^2/(8piG)] * delta * (2 + delta)
```

These give:
- Paper 3: Omega_K = delta*(4+3delta)/12 = 0.142 at delta_0=0.34
- K-essence: Omega_K = delta*(2+delta)/3 = 0.265 at delta_0=0.34

**The k-essence formula matches Omega_DM = 0.265, the Paper 3 formula does not.**

The discrepancy is roughly a factor of 2 at small delta. This means either:
1. Paper 3 has a formula error (most likely), OR
2. The BS stress-energy decomposition uses a different convention that absorbs a factor into the definition of rho_K

**URGENT**: Verify against BS2024 (arXiv:2404.06584) Equations 3.7-3.9 directly.

### 1.3 The Crisis Quantified

Using the CORRECT g(delta) = 2+delta:

| z | delta | g(delta) | g/g_0 | Excess over CDM |
|---|-------|----------|-------|-----------------|
| 0 | 0.340 | 2.340 | 1.000 | 0% (normalized) |
| 0.5 | 0.657 | 2.657 | 1.135 | +13.5% |
| 1 | 0.850 | 2.850 | 1.218 | +21.8% |
| 10 | 1.077 | 3.077 | 1.315 | +31.5% |
| 100 | 1.051 | 3.051 | 1.304 | +30.4% |
| 1100 | 0.819 | 2.819 | 1.205 | **+20.5%** |
| 3400 | 0.543 | 2.543 | 1.087 | +8.7% |
| 10000 | 0.277 | 2.277 | 0.973 | -2.7% |

At z=1100: +20.5% more DM than CDM when normalized to match at z=0.

Planck constrains omega_cdm to ~1% precision. A 20% excess at recombination is excluded at ~20 sigma. Even after renormalization at z=1100, the w_eff(a) shape produces ~6% TT residuals (>2 sigma, from CLASS Level 2 Run E).

---

## 2. Route-by-Route Analysis

### Route 1: Re-normalize at High Redshift

**Idea**: Fit rho_K to match CDM at z=1100 instead of z=0.

**Calculation**:
- ratio at z=1100: g(0.819)/g(0.340) = 2.819/2.340 = 1.205
- omega_K,0 = 0.120/1.205 = 0.0996 (as h^2) -> Omega_DM ~ 0.220
- At z=0: 17% less DM than CDM

**Does it solve the crisis?** Partially. The normalization at z=1100 matches, but the w_eff(a) shape still produces detectable ISW and peak-height shifts. CLASS Level 2 Run E (similar strategy with old formula) gave 5.8% max TT deviation -- still excluded at 2.3 sigma.

**Does it break anything?** It gives Omega_DM ~ 0.220 at z=0, which is 17% less than the CDM value. This:
- HELPS the S_8 tension (less DM -> less clustering)
- HELPS the H_0 tension: lower omega_m implies higher H_0 ~ 72.8 km/s/Mpc (from omega_m h^2 ~ const)
- May CONFLICT with galaxy cluster counts (which prefer omega_m ~ 0.3)
- Needs re-checking against BAO (but the previous DESI result preferred lower omega_m)

**Elegance cost**: Medium. Treating omega_K as a free parameter is already done in LCDM (omega_cdm is fitted). The cost is that the delta_0 value is now a derived quantity rather than giving a clean prediction.

**Testable?** Yes -- the prediction is Omega_DM(z=0) = 0.22 and H_0 ~ 73 km/s/Mpc.

**VERDICT**: Partial fix. Helps with cosmological tensions but doesn't fully resolve the w_eff shape problem.

---

### Route 2: Correct the f(delta) Formula

**Idea**: The Paper 3 formula f(delta) = (4+3delta)/4 may be incorrect. The standard k-essence derivation gives g(delta) = 2+delta.

**Evidence**:
1. rho_K = [Q K'(Q) - K(Q)] * normalization = delta*(2+delta) * normalization
2. This gives Omega_K = delta*(2+delta)/3 = 0.265 at delta_0 = 0.34 -- matches observations
3. The Paper 3 formula gives Omega_K = 0.142 -- doesn't match
4. The deep_intersection document independently uses rho_K = (rho_crit/3)(Q^2-1) = (rho_crit/3)*delta*(2+delta) -- consistent with the k-essence derivation

**Does it solve the crisis?** Reduces the excess from 29% to 20%. A 30% relative improvement, but 20% is still catastrophically excluded.

**Does it break anything?** No -- it fixes a potential formula error. All other results (c_s^2 = 0, conservation law, Sigma = 2 ln Q) are unaffected.

**Elegance cost**: Positive -- this IMPROVES the theory by fixing an error.

**Testable?** The corrected formula changes the w_eff bins for CLASS. This needs to be recomputed.

**CRITICAL ACTION**: Check BS2024 (arXiv:2404.06584) Eq. 3.8 to determine the correct formula. The f(delta) = (4+3delta)/4 factor may come from a different stress-energy decomposition convention in the BS paper that differs from the standard k-essence T_mu_nu.

**VERDICT**: Important correction, but insufficient alone. The crisis shrinks from ~13 sigma to ~8 sigma.

---

### Route 3: Modified Running Prescription

**Idea**: If mu runs faster than H/c, delta stays smaller at high z.

**Calculation** (mu propto H^alpha/c^alpha):

| alpha | delta(z=1100) | g-ratio | Excess |
|-------|---------------|---------|--------|
| 0.50 | 19283 | 8241 | Catastrophic |
| 0.75 | 126 | 55 | Catastrophic |
| 1.00 | 0.819 | 1.205 | +20.5% (standard) |
| 1.25 | 0.005 | 0.857 | -14.3% |
| 1.50 | 0.000 | 0.855 | -14.5% |

**Does it solve the crisis?** Yes, for alpha > 1.0. At alpha ~ 1.15, the excess at z=1100 would be ~0, and the theory would be safe.

**Does it break anything?**
- alpha > 1 means mu runs faster than H/c, so at early times mu_bg is much larger
- This changes delta_0 from 0.34 to a different value (need to refit)
- The formula delta = delta_0/(H^{2alpha}/H_0^{2alpha} * a^3) changes the cosmological predictions
- alpha = 1 is motivated by DPI + extensivity; alpha != 1 needs new theoretical motivation

**Elegance cost**: High. The running mu = H/c is motivated by three convergent arguments. Changing to H^{1.15}/c^{1.15} has no known motivation and introduces an ad hoc exponent.

**Testable?** In principle, different alpha values give different late-time w_eff(z) profiles, testable by DESI DR2.

**VERDICT**: Works numerically but theoretically unmotivated. Last resort.

---

### Route 4: Different K(Q)

**Idea**: K(Q) = mu^2(Q-1)^n with n > 2.

**Calculation**: For n > 2, the conservation law gives delta^{n-1} propto 1/(mu^2 a^3), so delta grows more slowly going backwards. But to maintain Omega_DM = 0.265, delta_0 must increase.

| n | delta_0 for Omega=0.265 | delta(z=1100) | Excess |
|---|------------------------|---------------|--------|
| 2 | 0.340 | 0.819 | +20.5% |
| 3 | 0.451 | 0.701 | +172% (WORSE!) |
| 4 | 0.523 | 0.701 | +164% (WORSE!) |
| 6 | 0.615 | 0.733 | +157% (WORSE!) |

**Does it solve the crisis?** NO. Higher n makes the excess WORSE when delta_0 is adjusted to match Omega_DM.

**Does it break anything?** Yes:
- n > 2 breaks the Gaussian Petz saturation (F = 1/(1+delta) only for n=2)
- n > 2 breaks the Fisher metric interpretation (K = QRE only for n=2)
- n > 2 may affect c_s^2 = 0 (needs checking)

**Elegance cost**: Extreme. Destroys the information-geometric foundation.

**Testable?** Different n gives different perturbation growth rates.

**VERDICT**: Unacceptable. Makes crisis worse and destroys theoretical elegance.

---

### Route 5: f(delta) as Perturbation Only

**Idea**: The g(delta) = 2+delta correction affects only perturbation equations, not the background Friedmann equation. The background rho_K is exactly proportional to a^{-3}.

**Analysis**: For this to work, the Khronon energy density in the Friedmann equation must be:
```
rho_K,background = c^2 I_0 / (16piG a^3)  [no g factor]
```
This gives Omega_K = delta_0/3 = 0.113. TOO SMALL for Omega_DM = 0.265.

To get Omega_DM = 0.265, need delta_0 = 0.795, giving delta(z=1100) ~ 1.9 -- the perturbative treatment breaks down.

**Does it solve the crisis?** No, because Omega_DM cannot match observations without the g(delta) factor.

**Alternative interpretation**: Perhaps the PHYSICAL rho_K includes the g factor (which correctly gives Omega_DM), but this rho_K should be DEFINED as CDM-like in the Friedmann equation (i.e., absorb the g correction into the definition of the "bare" density). This is a semantic trick: it amounts to redefining what you call rho_CDM vs rho_correction.

**Elegance cost**: Medium. It's not wrong per se, but it's an unusual decomposition.

**Testable?** If the background is exactly CDM, all CMB effects come from perturbation-level modifications only.

**VERDICT**: Cannot match Omega_DM. Dead end.

---

### Route 6: Small delta_0 + Other DM Source

**Idea**: If delta_0 << 0.34, the Khronon provides negligible cosmological DM. The bulk of Omega_DM comes from elsewhere.

**Calculation**: For delta_0 = 0.001: Omega_K = 0.0007 (0.3% of DM). CMB deviation < 0.01%. Completely safe.

**Does it solve the crisis?** Yes -- trivially.

**Does it break anything?** Destroys the Omega_DM prediction (the claim that Khronon naturally explains DM density). The bulk of DM must be:
- Actual CDM particles (returns to the standard problem)
- J(Y) sector contribution (but J -> Lambda at Y -> 0 on FRW)
- A different mechanism entirely

**Elegance cost**: Very high. The entire point of the theory is that K(Q) with mu_0 = H_0/c gives Omega_DM ~ rho_crit/3. Reducing delta_0 by factor 100 abandons this.

**Testable?** Galaxy rotation curves still work (from J(Y) sector), but cosmological DM is no longer explained.

**VERDICT**: Technically safe but intellectually devastating. Only acceptable if ALL other routes fail.

---

### Route 7: Two-Sector Model

**Idea**: K(Q) provides a fraction of cosmological DM; the rest comes from a different sector.

**Calculation**:

| K fraction | delta_0 | Excess at z=1100 | Total DM excess | Status |
|-----------|---------|-------------------|-----------------|--------|
| 10% | 0.039 | 2.7% | 0.3% | MARGINAL |
| 20% | 0.077 | 5.2% | 1.0% | MARGINAL |
| 30% | 0.113 | 7.5% | 2.3% | EXCLUDED |
| 50% | 0.182 | 11.8% | 5.9% | EXCLUDED |
| 100% | 0.340 | 20.5% | 20.5% | EXCLUDED |

**Does it solve the crisis?** Yes, if K(Q) provides < 10% of DM. But needs another DM source for the remaining 90%.

**Does it break anything?** Requires J(Y) or another sector to contribute ~90% of cosmological DM. In the BS framework, J(Y) -> Lambda on FRW backgrounds, so it cannot provide CDM-like matter. This would require:
- A new sector not in the BS action
- Or a non-trivial FRW limit of J(Y) that acts like dust
- Or actual CDM particles coexisting with the Khronon

**Elegance cost**: High. Defeats the purpose of explaining DM without particles.

**Testable?** The two sectors would have different perturbation properties, potentially distinguishable by CMB-S4.

**VERDICT**: Viable but costly. The K(Q) sector would be a small correction, not the DM.

---

### Route 8: Self-Consistent H(a)

**Idea**: The LCDM H(a) used to compute delta(a) is modified by the f(delta) correction itself. Maybe the self-consistent solution differs.

**Calculation**: Solving the coupled system iteratively:

| z | delta(LCDM) | delta(self-consistent) | g-ratio(LCDM) | g-ratio(self) |
|---|-------------|----------------------|----------------|---------------|
| 0 | 0.340 | 0.340 | 1.000 | 1.000 |
| 1 | 0.850 | 1.172 | 1.218 | 1.356 |
| 10 | 1.077 | 5.546 | 1.315 | 3.225 |
| 100 | 1.051 | 5.807 | 1.304 | 3.336 |
| 1100 | 0.819 | 2.268 | 1.205 | 1.824 |

**Does it solve the crisis?** NO. It makes the crisis DRAMATICALLY WORSE. The self-consistent solution has MORE DM at high z because the extra DM increases H(a), which reduces delta(a) less than expected, leading to even more DM in a positive feedback loop.

**Does it break anything?** N/A -- it doesn't work.

**VERDICT**: Definitively ruled out. Self-consistency exacerbates the problem.

---

### Route 9: w_tilde / w_EOS / w_eff Disambiguation

**Idea**: There are three different "w" quantities that have been conflated.

**The three w's**:

1. **w_tilde = delta/2** (BS effective EOS, from their stress-energy decomposition)
2. **w_EOS = p/rho = delta/(2+delta)** (thermodynamic EOS, pressure/density ratio)
3. **w_eff = -(1/3) d ln[rho_K a^3] / d ln a** (GDM effective w from density evolution)

These are ALL DIFFERENT:

| z | w_tilde | w_EOS | w_eff |
|---|---------|-------|-------|
| 0 | 0.170 | 0.145 | 0.785 |
| 1 | 0.425 | 0.298 | 0.278 |
| 10 | 0.538 | 0.350 | 0.001 |
| 100 | 0.525 | 0.344 | -0.013 |
| 1100 | 0.410 | 0.291 | -0.104 |

**Key insight**: w_eff (what CLASS needs) is NEGATIVE at z > 10 and POSITIVE at z < 10. This means:
- At early times: rho_K grows FASTER than a^{-3} going backward (w_eff < 0)
- At late times: rho_K dilutes FASTER than a^{-3} going forward (w_eff > 0)
- The net effect is MORE DM at z=1100 than CDM

**Does it solve the crisis?** No, but it CLARIFIES the problem. The CLASS Level 2 computation used w_eff derived from the OLD f formula. With the CORRECT g formula, the w_eff bins change:

| z | w_eff(old f) | w_eff(correct g) |
|---|-------------|-----------------|
| 0 | 0.825 | 0.785 |
| 1 | 0.297 | 0.278 |
| 100 | -0.014 | -0.013 |
| 1100 | -0.112 | -0.104 |

The values are slightly different (~7% at z=0). This means the CLASS results need to be rerun with the corrected w_eff.

**Does it break anything?** No. It's a clarification, not a modification.

**Elegance cost**: Positive -- resolves a conceptual confusion.

**VERDICT**: Essential clarification. Does not solve the crisis but properly frames it. All future CLASS runs must use w_eff, not w_tilde or w_EOS.

---

### Route 10: w_eff Sign Change -- Late-Time Phenomenology

**Observation**: w_eff transitions from negative (z > 10) to strongly positive (z < 3). At z=0, w_eff ~ 0.79, which is larger than any previous estimate.

This means:
- DM dilutes as a^{-3(1+0.79)} = a^{-5.4} at z=0 -- MUCH faster than CDM
- At z > 10, DM dilutes as a^{-3(1-0.01)} = a^{-2.97} -- slightly slower than CDM
- The cumulative effect: 20% more DM at z=1100

**Late-time implications**: w_eff ~ 0.8 at z=0 means the DM density TODAY is dropping rapidly. This:
- Reduces S_8 tension (less clustering at low z)
- Affects BAO distances (faster expansion at low z)
- Creates a detectable ISW effect
- Makes the theory look like "dark matter + dark energy" combined

**But**: w_eff = 0.79 at z=0 is extreme. Even w = 0.01 constant is excluded at 13 sigma by CLASS (Run G). The only reason this might survive is that w_eff is very small at z > 10 (where the CMB acoustic physics happens).

---

### Route 11: DBI Completion

**Idea**: Replace K(Q) = mu^2(Q-1)^2 with the DBI form:
```
K_DBI = mu^2 lambda^2 [sqrt(1 + (Q-1)^2/lambda^2) - 1]
```

For delta << lambda: K_DBI -> mu^2 delta^2/2 (recovers quadratic)
For delta >> lambda: K_DBI -> mu^2 lambda |delta| (linearizes)

In the linear regime, K' ~ const, so the conservation law gives a^3 * Q * K' = const, which means Q ~ const/a^3 and rho ~ const. This is EXACTLY CDM-like behavior when the DBI kicks in.

**Key insight**: If lambda_D is chosen so that delta > lambda_D during the matter era, the DBI completion AUTOMATICALLY regularizes the f(delta) correction, making rho_K truly proportional to a^{-3}.

**Challenge**: The DBI form changes the conservation law from mu^2 delta = I_0/(2a^3) to a different relation. The entire cosmological calculation needs to be redone. Also, K'(Q_0) != 0 for the DBI form, which means c_s^2 != 0 -- this may conflict with the TKS bound.

**Does it solve the crisis?** Potentially yes, if lambda_D ~ 0.1-0.5.

**Does it break anything?**
- c_s^2 may become non-zero (need to compute)
- The Petz saturation F = 1/(1+delta) may break for non-Gaussian K
- The Omega_DM prediction needs recalculation

**Elegance cost**: Medium. BS2024 already discusses the DBI completion as a natural UV regularization. The parameter lambda_D is an additional free parameter.

**VERDICT**: Promising but needs full calculation. This is the most promising STRUCTURAL escape route.

---

### Route 12: Reconsidering the BS Stress-Energy Formula

**Idea**: The discrepancy between f(delta)=(4+3delta)/4 and g(delta)=2+delta might be due to the BS stress-energy decomposition including extra terms (e.g., spatial gradient contributions) that are absent in the standard k-essence derivation.

**Analysis**: The BS action includes both K(Q) and J(Y) sectors. On FRW:
- Q != 1 (Khronon condensate), Y = 0 (no acceleration)
- The K sector stress-energy could include cross-terms with the metric that don't appear in the simple P(X) formalism

Specifically, the BS stress-energy tensor is:

```
T^K_{ab} = (c^4/(8piG)) * [2K'(Q) n_a n_b (1 - something) + 2K g_{ab} + ...]
```

The exact form involves the unit normal n_a and may have additional terms proportional to the extrinsic curvature. On FRW, these additional terms could contribute to T_00, modifying the energy density from the simple QK'-K formula.

**The f=(4+3delta)/4 factor MIGHT be the correct BS result** if:
```
rho_K^{BS} = (c^2 mu^2/(8piG)) * delta * (4+3delta)/4
```
with a different overall normalization. Then Omega_K = delta*(4+3delta)/12 = 0.142 at delta_0=0.34, which would require delta_0 = 0.576 to get Omega_K = 0.265.

At delta_0 = 0.576:
- delta(z=1100) ~ 1.39
- f_ratio = (4+3*1.39)/4 / (4+3*0.576)/4 = 8.17/5.73 = 1.43
- Excess at z=1100: 43% -- EVEN WORSE

**VERDICT**: If f=(4+3delta)/4 is correct, the crisis is worse (43% with delta_0=0.576), not better. The g(delta)=2+delta formula at least reduces the crisis. Either way, the fundamental problem remains.

---

## 3. The Three Most Promising Paths Forward

### Path A: Formula Correction + Renormalization + H_0 (Routes 1+2+9)

1. Verify the correct stress-energy formula by checking BS2024 directly
2. Use g(delta)=2+delta (reducing excess from 29% to 20%)
3. Renormalize omega_K at z=1100 (giving Omega_DM,0 ~ 0.22)
4. Accept lower Omega_DM as a FEATURE (helps H_0 and S_8 tensions)
5. Re-run CLASS with corrected w_eff bins
6. Check if the w_eff shape is compatible after full parameter refit

**Probability of success**: ~30%. The w_eff shape may still be excluded even after renormalization.

**H_0 prediction**: H_0 ~ 72.8 km/s/Mpc (in the right direction!)

### Path B: DBI Completion (Route 11)

1. Work out the full DBI cosmology: K_DBI with running mu = H/c
2. Compute the modified conservation law and delta(a) evolution
3. Find the lambda_D value that makes the CMB safe
4. Check whether c_s^2 remains sufficiently small
5. Compute the DBI Omega_DM prediction

**Probability of success**: ~40%. The DBI completion is physically motivated and may naturally regularize the problem.

**Cost**: One additional free parameter (lambda_D). Possible loss of c_s^2 = 0.

### Path C: Background-Perturbation Separation (Route 5, revised)

**New idea**: Argue that the Khronon background should be treated DIFFERENTLY from the perturbation level. Specifically:

In the Friedmann equation, the Khronon contributes through:
```
rho_K^{bg} = (c^2 I_0)/(16piG a^3)  [pure a^{-3}, no g correction]
```

The g(delta) = 2+delta factor enters only at the perturbation level, affecting how the Khronon density contrast grows.

**Justification**: The conservation law a^3 Q K' = I_0 gives mu^2 delta = I_0/(2a^3). The "bare" energy is proportional to mu^2 delta, which scales as a^{-3} exactly. The (2+delta) factor comes from the FULL T_00, which includes the "kinetic + potential" decomposition. Perhaps the "potential" part (the K(Q) = mu^2 delta^2 piece) should be treated as a correction.

This would give:
```
rho_K^{bg} = c^2 mu^2 delta / (4piG)    [= 2 * (rho_crit/3) * delta]
```
with Omega_K = 2*delta/3 = 0.227 at delta_0=0.34.

The excess at z=1100 would be:
```
ratio = delta(1100)/delta(0) * (at constant mu^2) = no excess if mu^2 delta = const
```

Wait -- if we use only the "kinetic" piece (proportional to mu^2 delta), and mu^2 delta = I_0/(2a^3) is exactly conserved, then rho_K propto 1/a^3 EXACTLY with NO correction.

**The excess comes entirely from the "potential" piece K(Q) = mu^2 delta^2.**

If we argue that the potential piece is a "perturbative correction" (valid for delta << 1, which holds at z > 3400), then the background is exactly CDM and the correction is small.

**Problem**: At z=0-1100, delta ~ 0.3-1.1, so the "perturbative" argument is marginal.

**Probability of success**: ~25%. Depends on a defensible theoretical argument for separating the two pieces.

---

## 4. What Needs to Be Done Immediately

### Priority 1: Verify the f(delta) formula
- Read BS2024 (arXiv:2404.06584) Eq. 3.7-3.9 directly
- Determine whether f(delta)=(4+3delta)/4 or g(delta)=(2+delta) is the correct BS formula
- This changes everything downstream

### Priority 2: Re-run CLASS with corrected w_eff
- Compute w_eff bins from the CORRECT formula
- Run CLASS Level 2 with 11-bin resolution
- Include renormalization (fit omega_K at z=1100)
- Full parameter refit (H_0, n_s, A_s, tau_reio as free)

### Priority 3: Explore the DBI completion
- Derive the DBI conservation law on FRW
- Compute delta(a) with K_DBI
- Find lambda_D that satisfies CMB constraints
- Check c_s^2 for DBI

### Priority 4: Assess the background-perturbation separation
- Is there a principled way to argue that only the "kinetic" piece enters the Friedmann equation?
- In Einstein-Aether theory and Horava gravity, how is the stress-energy decomposed?
- Check whether the standard Noether current gives the a^{-3} scaling without the g correction

---

## 5. Summary Table

| Route | Solves? | Breaks? | Elegance | Testable? | Priority |
|-------|---------|---------|----------|-----------|----------|
| 1. Re-normalize | Partial | H_0 help | Medium | Yes (H_0, S_8) | HIGH |
| 2. Formula fix | Reduces 29->20% | Nothing | High | Yes (CLASS) | HIGHEST |
| 3. Modified running | Yes (alpha>1) | Motivation | Low | Yes (w_eff) | LOW |
| 4. Higher n | No (worse) | Petz, Fisher | Very low | N/A | NONE |
| 5. Perturbation only | Omega_DM fails | Structure | Medium | Indirect | MEDIUM |
| 6. Small delta_0 | Trivially | DM prediction | Very low | N/A | LAST RESORT |
| 7. Two-sector | If frac<10% | J(Y) needed | Medium | CMB-S4 | MEDIUM |
| 8. Self-consistent | NO (worse!) | N/A | N/A | N/A | RULED OUT |
| 9. w disambiguation | Clarifies | Nothing | High | CLASS rerun | HIGH |
| 10. w_eff late-time | Clarifies | Nothing | High | DESI DR2 | HIGH |
| 11. DBI completion | Promising | c_s^2? | Medium | Full calc | HIGH |
| 12. BS decomposition | Unclear | Unclear | High | BS2024 check | HIGHEST |

---

## 6. The Deepest Question

The f(delta) crisis reveals a fundamental tension:

> **The Khronon theory needs delta_0 ~ O(1) to explain Omega_DM, but delta ~ O(1) means the Khronon is NOT a small perturbation of Minkowski. The (2+delta) correction IS the leading-order contribution to the energy density, not a perturbative correction.**

This means the "ghost condensation" regime (delta << 1) is NOT where the theory actually lives. At delta ~ 1, we are in the strongly nonlinear regime of K(Q), where:
- The quadratic approximation K ~ mu^2 delta^2 gets comparable corrections from higher-order terms
- The DBI completion becomes important
- The perturbative treatment of f(delta) breaks down

**The real question is not "how to escape the f(delta) crisis" but "what does the Khronon theory predict in the strongly nonlinear regime (delta ~ 1)".**

The quadratic K(Q) was derived from ghost condensation (small fluctuations around Q=1). But the cosmological solution has Q_0 = 1.34, which is a 34% displacement from the condensation point. At this displacement, higher-order terms in K(Q) are no longer negligible.

This suggests that the resolution may lie in the UV completion of K(Q), which is precisely the DBI form discussed in Route 11.

---

*Last updated: 2026-03-17*
*This analysis systematically evaluates 12 escape routes from the f(delta) CMB crisis.*
*The most promising paths are: formula verification (Route 12), DBI completion (Route 11), and renormalization + H_0 (Route 1).*
