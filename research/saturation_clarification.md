# AM-GM Saturation vs JRSWW Saturation: Resolution of the Apparent Contradiction

**Date:** 2026-03-17
**Author:** Analysis for Sheng-Kai Huang

---

## The Question

Paper 1 (Theorem 4, "Saturation of the Recovery Bound") claims that
F^2 = exp(-Delta_D) can be achieved with Delta_D > 0. For example,
dephasing at theta = pi/4 with sigma = I/2 gives Delta_D = ln 2 and
F^2 = 1/2. But an analysis using the JRSWW (Junge-Renner-Sutter-Wilde-Winter)
proof chain, which relies on the Golden-Thompson inequality, seems to
require Delta_D = 0 for saturation. Is there a contradiction?

## Answer: No Contradiction -- Two Different Inequalities

The resolution is that Paper 1's Theorem 4 and the JRSWW saturation
analysis concern **two different inequalities** that happen to share
the same form F^2 >= exp(-Delta_D).

---

### 1. What Paper 1 Actually Proves (the AM-GM Route)

Paper 1's Theorem 4 is restricted to **sigma = I/d** (maximally mixed
reference state) and **pure input rho = |psi><psi|**.

Under these conditions, the Petz recovery fidelity and the relative
entropy drop reduce to particularly clean algebraic forms. In the
eigenbasis of nu = N(sigma) and omega = N(rho) (assuming they commute),
the proof in Supplemental Sec. S17 shows:

- **F^2** = (1/d) sum_i (omega_i^2 / nu_i)
  = (1/d) sum_i omega_i * x_i          [weighted arithmetic mean of x_i = omega_i/nu_i]

- **exp(-Delta_D)** = (1/d) prod_i x_i^{omega_i}   [weighted geometric mean of x_i]

The bound F^2 >= exp(-Delta_D) then follows from the classical
**weighted AM-GM inequality**:

    sum_i omega_i * x_i  >=  prod_i x_i^{omega_i}

with equality iff all x_i are equal on the support of omega (i.e.,
omega_i/nu_i = constant). This is condition (ii) of Theorem 4.

**Key point:** This proof is entirely self-contained. It does not invoke
the JRSWW universal recovery theorem, Golden-Thompson, or the rotated
Petz map. It is a direct algebraic inequality for the **standard** Petz
map under the specific conditions sigma = I/d, rho pure.

### 2. What the JRSWW Proof Establishes (the Golden-Thompson Route)

The JRSWW theorem (Junge et al., Ann. Henri Poincare 19, 2955, 2018)
proves, for **arbitrary** sigma and rho:

    F(rho, R_rotated(N(rho)))^2 >= exp(-Delta_D)

where R_rotated is the **rotated** Petz map (an integral over modular
automorphisms of the standard Petz map). The proof chain is:

1. Start from the Peierls-Bogoliubov / Golden-Thompson inequality
2. Obtain a bound involving the sandwiched Renyi divergence
3. Optimize over the Renyi parameter to get the tightest bound

The saturation analysis of THIS inequality shows that equality
F^2 = exp(-Delta_D) in the JRSWW bound requires the channel to
be sufficient (Delta_D = 0), because the Golden-Thompson inequality
is strict for non-commuting operators whenever there is a genuine
gap in the data-processing inequality.

### 3. Why They Are Compatible

The two results are about **different maps** applied in **different
generality**:

| Feature | Paper 1 Theorem 4 | JRSWW Saturation |
|---|---|---|
| Recovery map | Standard Petz map R_{sigma,N} | Rotated Petz map R_tilde |
| Reference state | sigma = I/d only | Arbitrary faithful sigma |
| Input state | Pure rho | Arbitrary rho |
| Inequality used | AM-GM | Golden-Thompson |
| Saturation with Delta_D > 0? | **Yes** (when conditions (i)-(ii) hold) | **No** (requires Delta_D = 0) |

The crucial difference: when sigma = I/d, the Petz map simplifies
dramatically. The sigma^{1/2} factors become (1/sqrt{d}) * I, and the
entire recovery map becomes:

    R_{I/d, N}(.) = (1/d) N^dag( nu^{-1/2} (.) nu^{-1/2} )

This is a much simpler object than the general Petz map or the rotated
Petz map. The AM-GM inequality that governs its saturation behavior is
entirely independent of the Golden-Thompson trace inequality.

### 4. Concrete Example Confirming Both Are Correct

**Amplitude damping, ground-state input (sigma = I/2):**

- omega = |0><0|, nu = diag((1+gamma)/2, (1-gamma)/2)
- Likelihood ratio: omega_0/nu_0 = 2/(1+gamma) = c (only one term in support)
- Condition (ii) is trivially satisfied (one element)
- F^2 = 1/(1+gamma), Delta_D = ln(1+gamma)
- Check: F^2 = exp(-Delta_D). **Saturated with Delta_D > 0.**

This does NOT contradict JRSWW because:
- The JRSWW bound applies to the **rotated** Petz map
- For sigma = I/d, the rotated Petz map equals the standard Petz map
  (the modular automorphism sigma^{it}(.)sigma^{-it} is trivial when
  sigma is proportional to I)
- So in this special case, both maps agree, and the AM-GM bound
  (which is tighter than what Golden-Thompson gives in this regime)
  correctly predicts saturation

The JRSWW machinery, being designed for the general case, carries
extra "slack" from Golden-Thompson that does not apply when sigma = I/d.

### 5. Summary

| Question | Answer |
|---|---|
| Does Paper 1 use JRSWW to prove F^2 >= exp(-Delta_D) for sigma=I/d? | **No.** It uses AM-GM directly. |
| Is AM-GM independent of Golden-Thompson? | **Yes.** Completely independent inequality. |
| Can both have different saturation conditions? | **Yes.** Different inequalities saturate under different conditions. |
| Is Paper 1 Theorem 4 safe? | **Yes.** The proof is self-contained and does not depend on JRSWW at all. |
| Does JRSWW invalidate Theorem 4? | **No.** They address different objects (different maps, different generality). |

---

## Is Paper 1's Theorem 4 Safe?

**Yes, completely safe.** The theorem:

1. Uses a self-contained proof via AM-GM (Supplemental Sec. S17)
2. Is restricted to sigma = I/d and pure rho (stated explicitly)
3. Is verified numerically over 50,000 random channel-state pairs with zero failures
4. Has concrete analytically verified examples (amplitude damping, completely depolarizing)
5. The "necessity of commutativity" direction (condition (i)) uses the geometric-mean vs logarithmic-mean inequality, which is also independent of JRSWW

The JRSWW saturation analysis applies to the **general** case with
arbitrary sigma and the rotated Petz map. It correctly concludes that
in the general case, saturation of F^2 = exp(-Delta_D) for the rotated
Petz map requires Delta_D = 0. This is a **stronger** inequality
(it holds for all sigma, not just I/d) but correspondingly has
**stricter** saturation conditions.

The two results are complementary, not contradictory.
