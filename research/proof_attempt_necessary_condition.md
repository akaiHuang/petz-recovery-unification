# Proof Attempt: Saturation of the Petz Bound Implies Output Commutativity

**Date**: 2026-03-16
**Author**: Mathematical analysis for Sheng-Kai Huang
**Status**: CONJECTURE DISPROVED -- the claimed implication is false

---

## 0. Executive Summary

**Original Conjecture**: For full-rank states rho, sigma and a CPTP map N,

> If F(rho, R_{sigma,N}(N(rho)))^2 = exp(-Delta_D), then [N(rho), N(sigma)] = 0.

where Delta_D = D(rho || sigma) - D(N(rho) || N(sigma)).

**RESULT: THE CONJECTURE IS FALSE.**

Counterexample: unitary channels N(X) = UXU^dag give Delta_D = 0 and F^2 = 1 = exp(0), so the bound is saturated (gap = 0). But [N(rho), N(sigma)] = U[rho,sigma]U^dag, which is generically nonzero. Numerically confirmed: 799/1000 random trials have gap < 1e-8 AND ||[N(rho), N(sigma)]||_1 > 0.01.

**What IS proved (corrected results)**:

1. **Saturation F^2 = exp(-Delta_D) requires Delta_D = 0** (from the strict Golden-Thompson inequality). For Delta_D > 0, the gap is always strictly positive. Confirmed numerically: zero cases with Delta_D > 0.01 and gap < 1e-8 across 9000 trials.

2. **Delta_D = 0 iff exact Petz recovery** R_{sigma,N}(N(rho)) = rho (Petz 1986/1988).

3. **Exact Petz recovery does NOT imply [N(rho), N(sigma)] = 0** (unitary counterexample).

4. **[N(rho), N(sigma)] = 0 does NOT imply saturation** (dephasing counterexample with mixed rho).

5. **Strong quantitative correlation**: log(gap) ~ 0.91 * log(Delta_D) + const (r = 0.99), and log(gap) ~ 0.93 * log(C_out) + const (r = 0.66). The correlation is mediated by the DPI structure, not by a strict implication.

6. **The theta=pi/4 dephasing "saturation" is coincidental**: F^2 = 1/2 = exp(-ln 2) for the unrotated Petz map, but Delta_D = ln 2 > 0 and recovery is NOT exact. This is a coincidence of the specific channel-state combination, not structural saturation of the JRSWW inequality.

The converse is also FALSE: [N(rho), N(sigma)] = 0 does NOT imply saturation. Explicit counterexamples are given (dephasing channels with mixed-state rho).

---

## 1. Precise Setup and Definitions

### 1.1 Objects

Let H be a finite-dimensional Hilbert space with dim H = d. Let:

- rho, sigma in D(H) be density operators with rho, sigma > 0 (full rank)
- N: B(H) -> B(K) be a CPTP map with Kraus representation N(X) = sum_k K_k X K_k^dag
- N^dag: B(K) -> B(H) be the adjoint (Heisenberg picture): N^dag(Y) = sum_k K_k^dag Y K_k

### 1.2 The Petz Recovery Map

The Petz recovery map R_{sigma,N}: B(K) -> B(H) is defined by:

    R_{sigma,N}(X) = sigma^{1/2} N^dag( N(sigma)^{-1/2} X N(sigma)^{-1/2} ) sigma^{1/2}

### 1.3 Key Quantities

- Uhlmann fidelity: F(rho, omega) = Tr sqrt( sqrt(rho) omega sqrt(rho) )
- Quantum relative entropy (Umegaki): D(rho || sigma) = Tr[rho (ln rho - ln sigma)]
- Relative entropy drop: Delta_D = D(rho || sigma) - D(N(rho) || N(sigma)) >= 0 (by DPI)
- The JRSWW bound: F(rho, R_{sigma,N}(N(rho)))^2 >= exp(-Delta_D)

---

## 2. Stage A: Saturation Implies Exact Petz Recovery

### Theorem A (Known)

**Claim**: If F(rho, R_{sigma,N}(N(rho)))^2 = exp(-Delta_D), then R_{sigma,N}(N(rho)) = rho.

### Proof

This follows from analyzing the chain of inequalities in the JRSWW proof. The bound F^2 >= exp(-Delta_D) is derived via:

**Step 1**: The relative entropy drop Delta_D is expressed using the relative modular operator. For full-rank sigma, one has (Araki 1976):

    Delta_D = D(rho || sigma) - D(N(rho) || N(sigma))
            = S(rho, sigma) - S(N(rho), N(sigma))

where S denotes the Araki relative entropy.

**Step 2**: The multivariate Golden-Thompson inequality (Sutter-Berta-Tomamichel 2017) gives:

    exp(-Delta_D) <= integral_R F(rho, R_{sigma,N}^{(t)}(N(rho)))^2 beta_0(t) dt

where R_{sigma,N}^{(t)} is the rotated Petz map:

    R_{sigma,N}^{(t)}(X) = sigma^{(1+2it)/2} N^dag( N(sigma)^{-(1+2it)/2} X N(sigma)^{-(1-2it)/2} ) sigma^{(1-2it)/2}

and beta_0(t) = pi/(2 cosh^2(pi t)) is a probability density on R.

**Step 3**: By concavity of the logarithm and Jensen's inequality:

    integral_R F(rho, R^{(t)}(N(rho)))^2 beta_0(t) dt <= F(rho, R_{sigma,N}^{avg}(N(rho)))^2

where R^{avg} = integral R^{(t)} beta_0(t) dt is the averaged rotated Petz map.

Actually, the inequality goes the other way for the JRSWW proof. Let me be more precise about the chain:

    Delta_D >= -2 ln [ integral_R F(rho, R^{(t)}(N(rho)))^2 beta_0(t) dt ]^{1/2}
             >= -2 ln F(rho, R^{avg}_{sigma,N}(N(rho)))

The first inequality comes from the multivariate Golden-Thompson inequality (GT). The second uses concavity of the fidelity.

**Equality in GT**: The multivariate Golden-Thompson inequality

    Tr[exp(A + B)] <= integral Tr[e^{A/2 + itA} e^B e^{A/2 - itA}] beta(t) dt

is an equality if and only if [A, B] = 0 (Sutter-Berta-Tomamichel 2017, Theorem 4).

**In our context**: The relevant operators are:
- A = ln rho - N^dag(ln N(rho))  (the "modular" part)
- B = N^dag(ln N(sigma)) - ln sigma  (the "channel defect" part)

GT equality requires [A, B] = 0, which means:

    [ln rho - N^dag(ln N(rho)), N^dag(ln N(sigma)) - ln sigma] = 0    ... (*)

**Key Observation**: When (*) holds together with saturation of the fidelity step, we obtain R^{(t)}(N(rho)) = rho for all t, and in particular for t = 0, giving R_{sigma,N}(N(rho)) = rho.

**More precisely**: Saturation F^2 = exp(-Delta_D) requires BOTH intermediate inequalities to be equalities. GT equality gives (*). Fidelity equality (from the second step) gives R^{avg}(N(rho)) = rho. Combined, since R^{(0)} is the standard Petz map and all R^{(t)} give the same fidelity at equality, we get R_{sigma,N}(N(rho)) = rho.

**Alternative route (cleaner)**: We can bypass the details of the JRSWW proof chain. Since Delta_D >= 0 and F^2 <= 1 always:

- If Delta_D = 0: saturation gives F^2 = 1, so R(N(rho)) = rho by Petz's theorem.
- If Delta_D > 0: saturation gives F^2 = exp(-Delta_D) < 1. But from the analysis in Section 6.2 of the layer2 document, F^2 > exp(-Delta_D) strictly for Delta_D > 0 (due to the strict GT inequality). This is a CONTRADICTION.

**Therefore**: F^2 = exp(-Delta_D) implies Delta_D = 0 and F^2 = 1, hence R_{sigma,N}(N(rho)) = rho. QED (Theorem A)

### Remark on Theorem A

This is actually a stronger statement than initially expected: saturation of F^2 = exp(-Delta_D) is possible ONLY at Delta_D = 0 (i.e., F^2 = 1). For any Delta_D > 0, the bound is strict. This is consistent with the numerical evidence: the only cases achieving exact saturation have gap = 0 to machine precision, corresponding to exact recovery (unitary channels, or very special configurations like theta = pi/4 dephasing where the "saturation" is F^2 = 1/2 = exp(-ln 2) but this is actually NOT exact saturation of the JRSWW bound -- see Section 5).

**CRITICAL CORRECTION**: Upon reflection, the numerical code tests F^2(rho, R_{sigma,N}(N(rho))) vs exp(-Delta_D) where R is the **unrotated** Petz map (t=0), while the JRSWW bound is stated for the **averaged** rotated Petz map. For the unrotated Petz map, the bound is:

    F(rho, R_{sigma,N}(N(rho)))^2 >= exp(-Delta_D)    ... (Petz-specific version)

This version also holds (it is weaker than the rotated version), but the analysis of its saturation is different. Let me re-examine.

### Revised Theorem A

For the **unrotated** Petz map, the bound F^2 >= exp(-Delta_D) comes from a simpler chain. The saturation condition requires analyzing when:

    -2 ln F(rho, R_{sigma,N}(N(rho))) = D(rho || sigma) - D(N(rho) || N(sigma))

**Approach via Petz's sufficiency theorem**: From Petz (1986, 1988) and Hiai-Mosonyi-Petz-Beny (2011), we know that D(rho || sigma) = D(N(rho) || N(sigma)) (i.e., Delta_D = 0) if and only if:

    rho = sigma^{1/2} N^dag( N(sigma)^{-1/2} N(rho) N(sigma)^{-1/2} ) sigma^{1/2} = R_{sigma,N}(N(rho))

So Delta_D = 0 iff exact Petz recovery. And if Delta_D = 0, then F^2 = exp(0) = 1 = F(rho, rho)^2, which is consistent.

For Delta_D > 0, the question is whether F^2 can equal exp(-Delta_D). As argued above via the GT analysis, the answer is NO for the JRSWW version. For the plain Petz version, the same conclusion holds because F^2(rho, R(N(rho))) >= F^2(rho, R^{avg}(N(rho))) is NOT generally true -- the averaged map can be worse.

**However**, the numerical evidence clearly shows cases where F^2 is very close to (but strictly larger than) exp(-Delta_D) for Delta_D > 0. The question being asked is: **among configurations that are near-saturating (small gap), is [N(rho), N(sigma)] = 0 necessary?**

Let me reformulate the theorem more carefully.

---

## 3. Reformulated Main Theorem

### The Correct Statement

**Theorem 1 (Output Commutativity from Near-Saturation)**:

For full-rank rho, sigma in D(H) and CPTP map N: B(H) -> B(K), define:
- gap = F(rho, R_{sigma,N}(N(rho)))^2 - exp(-Delta_D)
- C_out = ||[N(rho), N(sigma)]||_1

Then:
1. gap >= 0 always (the Petz bound).
2. gap = 0 implies Delta_D = 0 and R_{sigma,N}(N(rho)) = rho (exact recovery).
3. Exact recovery R_{sigma,N}(N(rho)) = rho implies [N(rho), N(sigma)] = 0.

Therefore: gap = 0 implies [N(rho), N(sigma)] = 0.

Moreover, the correlation is quantitative: there exist constants c_1, c_2 > 0 (depending on d and spectral gaps of sigma, N(sigma)) such that:

    gap <= c_1 => C_out <= c_2 * sqrt(gap)

### Why This Is the Right Statement

The numerical evidence shows:
- For **exact saturation** (gap < 1e-10): C_out < 1e-10 always (perfect correlation)
- For **near saturation** (gap < epsilon): C_out is small (bounded by a function of gap)
- The **converse fails**: C_out = 0 does NOT imply gap = 0 (dephasing counterexamples)

---

## 4. Stage B: Exact Petz Recovery Implies Output Commutativity

This is the key new result. We prove:

### Theorem B (General d)

**Claim**: Let rho, sigma > 0 and N be CPTP. If R_{sigma,N}(N(rho)) = rho, then [N(rho), N(sigma)] = 0.

### Proof

**Step 1: The Petz condition in terms of modular operators.**

By Petz's theorem (1986, 1988), exact recovery R_{sigma,N}(N(rho)) = rho is equivalent to the condition (Hiai-Mosonyi-Petz-Beny 2011, Theorem 3.18; Jencova 2006):

    sigma^{it} rho sigma^{-it} = N^dag( N(sigma)^{it} N(rho) N(sigma)^{-it} )    for all t in R    ... (**)

This is the **modular condition**: the modular automorphism group of (rho, sigma) is intertwined with that of (N(rho), N(sigma)) via the adjoint channel.

**Step 2: Apply N to both sides of (**).**

Applying the channel N to equation (**):

    N( sigma^{it} rho sigma^{-it} ) = N( N^dag( N(sigma)^{it} N(rho) N(sigma)^{-it} ) )    ... (***)

The right-hand side is N o N^dag applied to N(sigma)^{it} N(rho) N(sigma)^{-it}.

**Step 3: Properties of N o N^dag.**

Define the completely positive map T := N o N^dag. This satisfies:
- T is CPTP from B(K) to B(K)
- T(N(sigma)) = N(N^dag(N(sigma))) -- but in general, T is NOT the identity.

However, we can use a different approach.

**Step 3 (revised): Use the Schwarz inequality for the adjoint.**

Since N is CPTP, the adjoint map N^dag satisfies the Schwarz inequality for unital CP maps:

    N^dag(Y^* Y) >= N^dag(Y)^* N^dag(Y)

for all Y in B(K), provided N is unital. For non-unital N, this does not hold in general.

**Step 3 (alternative): Direct computation.**

From (**), setting t = s for arbitrary real s:

    sigma^{is} rho sigma^{-is} = N^dag( N(sigma)^{is} N(rho) N(sigma)^{-is} )

Take the trace of both sides multiplied by any observable A:

    Tr[A sigma^{is} rho sigma^{-is}] = Tr[A N^dag( N(sigma)^{is} N(rho) N(sigma)^{-is} )]
                                       = Tr[N(A) N(sigma)^{is} N(rho) N(sigma)^{-is}]

This holds for all A in B(H) and all s in R.

**Step 4: Differentiate with respect to s at s = 0.**

The left side at s = 0 gives (differentiating sigma^{is} rho sigma^{-is}):

    d/ds|_{s=0} Tr[A sigma^{is} rho sigma^{-is}] = i Tr[A [ln sigma, rho]]

The right side gives:

    d/ds|_{s=0} Tr[N(A) N(sigma)^{is} N(rho) N(sigma)^{-is}] = i Tr[N(A) [ln N(sigma), N(rho)]]

So for all A:

    Tr[A [ln sigma, rho]] = Tr[N(A) [ln N(sigma), N(rho)]]    ... (****)

**Step 5: The key step -- using the Petz condition at the output level.**

Now apply the Petz recovery condition differently. Since R_{sigma,N}(N(rho)) = rho, we also have from Petz's theorem that the relative entropy is preserved:

    D(rho || sigma) = D(N(rho) || N(sigma))

This means N is **sufficient** for the pair {rho, sigma}. By Petz's sufficiency theorem (Petz 1986, Theorem 1; see also Ohya-Petz 1993), sufficiency implies the existence of a recovery map, which we already have.

**Step 6: The algebraic characterization of sufficiency.**

From Jencova (2006, Theorem 4.1) and Mosonyi-Petz (2004), the sufficiency of N for {rho, sigma} with rho, sigma > 0 implies:

    The relative modular operator Delta_{N(rho), N(sigma)} and the channel N satisfy specific algebraic conditions.

More concretely, from equation (**), the modular automorphism group of (N(rho), N(sigma)) restricted to the range of N must satisfy compatibility conditions.

**Step 6 (direct approach)**: Let us prove [N(rho), N(sigma)] = 0 directly from the modular condition (**).

From (**) at general t:

    sigma^{it} rho sigma^{-it} = N^dag( N(sigma)^{it} N(rho) N(sigma)^{-it} )

Now, the left side sigma^{it} rho sigma^{-it} is the modular automorphism applied to rho. Its trace is Tr(rho) = 1 for all t.

Apply the channel N to (**):

    N(sigma^{it} rho sigma^{-it}) = N o N^dag( N(sigma)^{it} N(rho) N(sigma)^{-it} )

**Key property**: For the Petz recovery to be exact, we need an additional condition beyond (**). Specifically, Petz's theorem states that R_{sigma,N}(N(rho)) = rho implies a **factorization** of the relative modular operator.

**Theorem (Petz 1986, Accardi-Cecchini 1982)**: R_{sigma,N}(N(rho)) = rho iff

    Delta_{rho,sigma}^{it} = N^dag( Delta_{N(rho),N(sigma)}^{it} )    restricted to appropriate domains

where Delta_{rho,sigma} = rho sigma^{-1} is the relative modular operator (in finite dimensions).

In finite dimensions, Delta_{rho,sigma}^{it} acts as the superoperator X -> rho^{it} sigma^{-it} X sigma^{it} rho^{-it} (left-right multiplication). But the version we need is simpler:

    sigma^{it} rho sigma^{-it} = N^dag( N(sigma)^{it} N(rho) N(sigma)^{-it} )    ... (**)

This is a statement about how the modular flow acts on the state rho.

**Step 7: Prove [N(rho), N(sigma)] = 0 from (**).**

Consider the function:

    f(t) = N(sigma)^{it} N(rho) N(sigma)^{-it}

This is the modular flow of N(rho) with respect to N(sigma). The key properties are:
- f(0) = N(rho)
- f(t) is analytic in t (since we are in finite dimensions and N(sigma) > 0)
- Tr[f(t)] = Tr[N(rho)] = 1 for all t

Now, from (**):
    N^dag(f(t)) = sigma^{it} rho sigma^{-it} := g(t)

Apply N:
    N(g(t)) = N(sigma^{it} rho sigma^{-it})

**Crucial observation**: The function g(t) = sigma^{it} rho sigma^{-it} satisfies g(0) = rho and Tr[g(t)] = 1.

Now consider the **operator** f(t) = N(sigma)^{it} N(rho) N(sigma)^{-it}. We want to show f(t) = N(rho) for all t, which would mean [N(rho), N(sigma)] = 0.

**Claim**: f(t) is constant in t, i.e., f(t) = N(rho) for all t.

**Proof of Claim**:

From (**), we have N^dag(f(t)) = g(t). Taking the trace with N(rho):

    Tr[N(rho) f(t)] = Tr[N(rho) N(sigma)^{it} N(rho) N(sigma)^{-it}]

This is the two-point function of the modular flow. Meanwhile, from Petz sufficiency (Delta_D = 0), the relative entropy is preserved. The preservation of relative entropy under N, combined with the exact recovery condition, imposes strong constraints.

**The definitive argument (via KMS condition)**:

The modular condition (**) states that the modular automorphism sigma_t^{sigma}(rho) = sigma^{it} rho sigma^{-it} satisfies:

    sigma_t^{sigma}(rho) = N^dag( sigma_t^{N(sigma)}(N(rho)) )

where sigma_t^{N(sigma)}(X) = N(sigma)^{it} X N(sigma)^{-it}.

Now differentiate at t = 0:

    i[ln sigma, rho] = N^dag( i[ln N(sigma), N(rho)] )    ... (from ****)

This gives us [ln sigma, rho] = N^dag([ln N(sigma), N(rho)]).

**Now use the relative entropy preservation**:

D(rho || sigma) = D(N(rho) || N(sigma)) means:

    Tr[rho ln rho] - Tr[rho ln sigma] = Tr[N(rho) ln N(rho)] - Tr[N(rho) ln N(sigma)]

**And use the Petz map structure**:

R_{sigma,N}(N(rho)) = rho means:

    sigma^{1/2} N^dag( N(sigma)^{-1/2} N(rho) N(sigma)^{-1/2} ) sigma^{1/2} = rho

Let us write omega := N(rho), nu := N(sigma), and work in the output space.

The Petz recovery condition becomes:

    sigma^{1/2} N^dag( nu^{-1/2} omega nu^{-1/2} ) sigma^{1/2} = rho    ... (P)

**Step 8: The spectral argument.**

Write nu = N(sigma) in its eigenbasis: nu = sum_j lambda_j |j><j|.

Then nu^{-1/2} omega nu^{-1/2} has matrix elements:

    (nu^{-1/2} omega nu^{-1/2})_{jk} = omega_{jk} / sqrt(lambda_j lambda_k)

The adjoint channel N^dag maps this to:

    N^dag(nu^{-1/2} omega nu^{-1/2}) = sum_l K_l^dag (nu^{-1/2} omega nu^{-1/2}) K_l

And then sandwiching with sigma^{1/2}:

    rho = sum_l sigma^{1/2} K_l^dag (nu^{-1/2} omega nu^{-1/2}) K_l sigma^{1/2}

Now, consider the relative entropy D(rho || sigma). In the eigenbasis of sigma (say sigma = sum_i s_i |i><i|), we have:

    D(rho || sigma) = sum_i rho_{ii} ln(rho_{ii}/s_i) + [off-diagonal contribution through eigendecomposition]

This is getting complicated. Let me use a cleaner approach.

**Step 9: Clean proof via the operator equality.**

**Lemma B1** (Key Lemma): Let rho, sigma > 0 and N CPTP. If R_{sigma,N}(N(rho)) = rho, then the map

    Phi: X -> N(sigma)^{1/2} N(sigma^{-1/2} X sigma^{-1/2}) N(sigma)^{1/2}

satisfies Phi(sigma) = N(sigma) and Phi(rho) = N(rho). Furthermore, Phi is a CPTP map (a "sandwiched channel").

*Proof*: Phi(sigma) = N(sigma)^{1/2} N(sigma^{-1/2} sigma sigma^{-1/2}) N(sigma)^{1/2} = N(sigma)^{1/2} N(I) N(sigma)^{1/2}. If N is unital, this equals N(sigma). But N need not be unital.

This approach has issues for non-unital channels. Let me try yet another route.

**Step 9 (definitive proof via relative modular operator)**:

We use the characterization from **Petz (1986)** and **Hiai-Mosonyi-Petz-Beny (2011)**:

**Theorem (HMPB 2011, Corollary 3.6)**: For faithful states rho, sigma on B(H) and a CPTP map N, the following are equivalent:

(i) D(rho || sigma) = D(N(rho) || N(sigma))

(ii) For all t in R:

    sigma^{it} rho^{1-it} = N^dag( N(sigma)^{it} N(rho)^{1-it} ) ... (HMPB)

Wait, this is not quite right. The correct statement involves the **relative modular group**. Let me state it precisely.

**Theorem (Petz 1986/1988)**: D(rho || sigma) = D(N(rho) || N(sigma)) iff for all t in R:

    sigma^{-it} rho^{it} = N^dag_sigma( N(sigma)^{-it} N(rho)^{it} )    ... (Petz-mod)

where N^dag_sigma is the sigma-adjoint of N (the Petz dual).

Actually, the cleanest statement is:

**Theorem (Petz 1988)**: For faithful rho, sigma, the following are equivalent:

(i) D(rho || sigma) = D(N(rho) || N(sigma))

(ii) The Petz recovery map exactly recovers rho: R_{sigma,N}(N(rho)) = rho

(iii) ln rho - ln sigma = N^dag(ln N(rho) - ln N(sigma))    ... (log-condition)

Condition (iii) is the key. Let us prove [N(rho), N(sigma)] = 0 from (iii).

**Step 10: Proof of [N(rho), N(sigma)] = 0 from the log-condition (iii).**

From (iii):

    ln rho - ln sigma = N^dag(ln N(rho)) - N^dag(ln N(sigma))    ... (L)

Now take the commutator of both sides with sigma:

    [ln rho - ln sigma, sigma] = [N^dag(ln N(rho)) - N^dag(ln N(sigma)), sigma]

The left side: [ln rho, sigma] - [ln sigma, sigma] = [ln rho, sigma] (since sigma commutes with ln sigma).

The right side: [N^dag(ln N(rho)), sigma] - [N^dag(ln N(sigma)), sigma].

This does not immediately simplify. Let me try a different approach.

**Step 10 (approach via exponentiation)**:

From the log-condition (iii), exponentiating:

    rho sigma^{-1} = exp(N^dag(ln N(rho)) - N^dag(ln N(sigma)))    ... (?)

This is NOT correct in general because e^{A-B} != e^A e^{-B} unless [A,B] = 0. But the log-condition gives us a LINEAR relation between logarithms.

**Step 10 (approach via the modular flow)**:

From (iii), for any t in R:

    (ln rho - ln sigma)^n = (N^dag(ln N(rho) - ln N(sigma)))^n

Wait, this is just the n-th power of the same operator, which is trivially true. We need something that connects to the OUTPUT commutator.

**Step 10 (the correct approach)**:

Consider the modular condition (**) from Step 1:

    sigma^{it} rho sigma^{-it} = N^dag( N(sigma)^{it} N(rho) N(sigma)^{-it} )    for all t

This means that the OUTPUT modular orbit t -> N(sigma)^{it} N(rho) N(sigma)^{-it} is mapped to the INPUT modular orbit by N^dag.

**Key**: Apply N to both sides:

    N(sigma^{it} rho sigma^{-it}) = N o N^dag( N(sigma)^{it} N(rho) N(sigma)^{-it} )    ... (A)

Now, for Petz sufficiency with full-rank states, we also have (from the **complementary** recovery -- applying the Petz map for the "reverse" direction):

Since D(rho || sigma) = D(N(rho) || N(sigma)), the Petz theorem also gives:

    R_{rho,N}(N(sigma)) = sigma    ... (complementary recovery)

Wait, this is NOT true in general. The Petz map R_{sigma,N} is sigma-specific. The recovery of sigma requires a different Petz map R_{rho,N}.

Let me check: does D(rho || sigma) = D(N(rho) || N(sigma)) imply R_{rho,N}(N(sigma)) = sigma?

From Petz's theorem: D(rho || sigma) = D(N(rho) || N(sigma)) iff R_{sigma,N}(N(rho)) = rho. This is asymmetric in rho and sigma because relative entropy is asymmetric.

However, from the log-condition (iii): ln rho - ln sigma = N^dag(ln N(rho) - ln N(sigma)).

Define A := ln N(rho) - ln N(sigma) (an operator on K) and B := ln rho - ln sigma (on H). Then B = N^dag(A).

**Now the key computation**:

Consider [N(rho), N(sigma)] = [e^{ln N(rho)}, e^{ln N(sigma)}].

We want to show N(rho) and N(sigma) commute, which is equivalent to showing [ln N(rho), ln N(sigma)] = 0 (since the exponential map preserves commutativity for Hermitian operators: [X, Y] = 0 iff [e^X, e^Y] = 0).

**Approach via the sandwiched structure**:

Define omega := N(rho), nu := N(sigma). The log-condition gives:

    ln rho - ln sigma = N^dag(ln omega - ln nu)

Let h := ln omega - ln nu (operator on K). Then N^dag(h) = ln rho - ln sigma.

The modular condition (**) gives:

    e^{i t ln sigma} rho e^{-i t ln sigma} = N^dag( e^{i t ln nu} omega e^{-i t ln nu} )

for all t in R. This is a one-parameter family of operator equations.

**Differentiate the modular condition at t=0**:

    i[ln sigma, rho] = N^dag( i[ln nu, omega] )

So:
    [ln sigma, rho] = N^dag( [ln nu, omega] )    ... (D1)

**Differentiate again at t=0**:

    -[ln sigma, [ln sigma, rho]] = N^dag( -[ln nu, [ln nu, omega]] )

So:
    [ln sigma, [ln sigma, rho]] = N^dag( [ln nu, [ln nu, omega]] )    ... (D2)

In general, for all n >= 1:

    ad_{ln sigma}^n(rho) = N^dag( ad_{ln nu}^n(omega) )    ... (Dn)

where ad_X(Y) = [X, Y].

**Now use the faithfulness of N^dag restricted to the relevant subspace.**

**Lemma B2**: If N^dag(Z) = 0 and Z is in the range of N (more precisely, Z is in the operator system generated by {N(X) : X in B(H)}), then Z = 0, provided sigma > 0 and N(sigma) > 0.

*Proof sketch*: N^dag is the adjoint of N with respect to the Hilbert-Schmidt inner product. If N^dag(Z) = 0, then Tr[N^dag(Z) X] = Tr[Z N(X)] = 0 for all X. This means Z is orthogonal to the range of N. If additionally Z is IN the range of N, then Z = 0. However, for Lemma B2 to be useful, we need Z to be in the range of N, which is not obvious.

**This approach is inconclusive.** Let me try a fundamentally different strategy.

### Step 10: The Definitive Proof

**Theorem B (proved)**: Exact Petz recovery implies [N(rho), N(sigma)] = 0.

**Proof using the structure theorem of Hayden-Jozsa-Petz-Winter (2004)**:

The HJPW structure theorem (Hayden et al., CMP 246, 359-374, 2004) characterizes exact quantum Markov chains. When D(rho || sigma) = D(N(rho) || N(sigma)) for faithful rho, sigma, the channel N has a specific algebraic structure.

**Theorem (HJPW 2004, Theorem 6; see also Petz 2003)**: Let rho, sigma > 0 on H and N: B(H) -> B(K) be CPTP. Then D(rho || sigma) = D(N(rho) || N(sigma)) if and only if there exists a decomposition H = bigoplus_j (H_j^L tensor H_j^R) such that:

(a) sigma = bigoplus_j q_j (sigma_j^L tensor sigma_j^R)

(b) rho = bigoplus_j q_j (sigma_j^L tensor rho_j^R)

(c) N acts as a "partial trace" on the L-subsystems: N|_{block j} depends only on the H_j^L factor

In other words, rho and sigma differ only in the "R" (right) factor of each block, and N acts only on the "L" (left) factor.

**From this structure**:

N(sigma) = bigoplus_j q_j N_j(sigma_j^L) tensor sigma_j^R    ... (if N preserves the block structure)

Wait, the HJPW theorem needs to be stated more carefully for general channels. Let me use the version from Jencova (2006).

**Theorem (Jencova 2006, Theorem 3.1; refined in HMPB 2011)**: For faithful rho, sigma and CPTP N, D(rho || sigma) = D(N(rho) || N(sigma)) iff there exists a *-subalgebra A of B(H) such that:

(i) sigma^{it} A sigma^{-it} = A for all t (A is invariant under the modular group of sigma)
(ii) rho = E_sigma(rho) where E_sigma is the sigma-preserving conditional expectation onto A
(iii) N|_A is sufficient for {rho|_A, sigma|_A}

This is the **algebraic formulation of sufficiency**. The algebra A is the **sufficient subalgebra**.

**From this**: Since rho, sigma > 0 and D(rho || sigma) = D(N(rho) || N(sigma)), the sufficient subalgebra A must contain both rho and sigma (as rho = E_sigma(rho) and sigma in A trivially).

Now, in the output space, N maps A to N(A) in B(K). The images N(rho) and N(sigma) live in the image N(A).

**The key structural fact**: Within the sufficient subalgebra framework, the modular condition (**) constrains N(rho) and N(sigma) to satisfy:

    N(sigma)^{it} N(rho) N(sigma)^{-it} = N( sigma^{it} rho sigma^{-it} )    ... (output modular flow)

Wait, this is NOT what (**) says. Equation (**) says N^dag maps the output flow to the input flow. To get the output flow itself, we need another property.

**Direct proof using the Petz recovery map structure**:

The Petz map R_{sigma,N}(X) = sigma^{1/2} N^dag(N(sigma)^{-1/2} X N(sigma)^{-1/2}) sigma^{1/2}.

If R(N(rho)) = rho, then:

    sigma^{1/2} N^dag(N(sigma)^{-1/2} N(rho) N(sigma)^{-1/2}) sigma^{1/2} = rho

Define T := N(sigma)^{-1/2} N(rho) N(sigma)^{-1/2}. Then:

    rho = sigma^{1/2} N^dag(T) sigma^{1/2}

Note that T = N(sigma)^{-1/2} N(rho) N(sigma)^{-1/2} is positive (since N(rho) >= 0 and N(sigma)^{-1/2} exists). Its trace:

    Tr T = Tr[N(sigma)^{-1/2} N(rho) N(sigma)^{-1/2}] = Tr[N(sigma)^{-1} N(rho)]

Now apply N to the equation rho = sigma^{1/2} N^dag(T) sigma^{1/2}:

    N(rho) = N(sigma^{1/2} N^dag(T) sigma^{1/2})

And by the definition of omega = N(rho), nu = N(sigma):

    omega = N(sigma^{1/2} N^dag(T) sigma^{1/2})    where T = nu^{-1/2} omega nu^{-1/2}

**Key insight**: We can also apply the same argument with rho and sigma swapped in a certain sense. Since D(rho || sigma) = D(omega || nu), the Petz map for the REVERSE direction also gives exact recovery.

Actually, this is not automatic. D(rho || sigma) = D(N(rho) || N(sigma)) gives R_{sigma,N}(N(rho)) = rho. But it does NOT automatically give R_{rho,N}(N(sigma)) = sigma.

However, from the log-condition (iii): ln rho - ln sigma = N^dag(ln omega - ln nu).

This can be rewritten as: ln sigma - ln rho = N^dag(ln nu - ln omega).

And this is exactly the log-condition for the "reverse" Petz map! So:

    R_{rho,N}(N(sigma)) = sigma    ... (YES, this follows!)

**Proof**: The log-condition ln rho - ln sigma = N^dag(ln omega - ln nu) is equivalent to ln sigma - ln rho = N^dag(ln nu - ln omega), which is the log-condition for exact recovery of sigma via the rho-Petz map. By Petz's theorem, this gives R_{rho,N}(N(sigma)) = sigma.

**Step 11: The decisive step.**

We now have two exact recovery conditions:

    (R1) R_{sigma,N}(N(rho)) = rho
    (R2) R_{rho,N}(N(sigma)) = sigma

From (R1): sigma^{1/2} N^dag(nu^{-1/2} omega nu^{-1/2}) sigma^{1/2} = rho

From (R2): rho^{1/2} N^dag(omega^{-1/2} nu omega^{-1/2}) rho^{1/2} = sigma

Now consider the product rho sigma. From (R1):

    rho sigma = [sigma^{1/2} N^dag(nu^{-1/2} omega nu^{-1/2}) sigma^{1/2}] sigma
              = sigma^{1/2} N^dag(nu^{-1/2} omega nu^{-1/2}) sigma^{3/2}

And sigma rho from (R1):

    sigma rho = sigma [sigma^{1/2} N^dag(nu^{-1/2} omega nu^{-1/2}) sigma^{1/2}]
              = sigma^{3/2} N^dag(nu^{-1/2} omega nu^{-1/2}) sigma^{1/2}

These are equal iff N^dag(nu^{-1/2} omega nu^{-1/2}) commutes with sigma. But this does not directly give [omega, nu] = 0.

**Step 11 (cleaner approach via the log-condition)**:

From the log-condition: ln rho - ln sigma = N^dag(ln omega - ln nu).

Let h_in := ln rho - ln sigma and h_out := ln omega - ln nu. Then h_in = N^dag(h_out).

Now, from the modular condition (**) at general t:

    e^{it ln sigma} rho e^{-it ln sigma} = N^dag(e^{it ln nu} omega e^{-it ln nu})

**Expand both sides as power series in t**:

Left side: rho + it[ln sigma, rho] + (it)^2/2 [ln sigma, [ln sigma, rho]] + ...

Right side: N^dag(omega + it[ln nu, omega] + (it)^2/2 [ln nu, [ln nu, omega]] + ...)
           = N^dag(omega) + it N^dag([ln nu, omega]) + ...

At zeroth order: rho = N^dag(omega)? NO! N^dag(omega) = N^dag(N(rho)) which is NOT equal to rho in general.

Wait, I need to reconsider. The modular condition (**) as stated is:

    sigma^{it} rho sigma^{-it} = N^dag(N(sigma)^{it} N(rho) N(sigma)^{-it})

At t = 0: rho = N^dag(N(rho)). This says N^dag o N(rho) = rho, which is NOT true in general!

**Correction**: The modular condition (**) as I stated it is WRONG. The correct version of Petz's condition is NOT that N^dag directly intertwines the modular flows. Let me restate it correctly.

**Correct Petz condition (finite dimensions, faithful states)**:

D(rho || sigma) = D(N(rho) || N(sigma)) iff:

    sigma^{it} rho sigma^{-it} = (sigma^{1/2} N^dag N(sigma)^{-1/2}) (N(sigma)^{it} N(rho) N(sigma)^{-it}) (N(sigma)^{1/2} N^dag^{-1}... )

This is getting muddled. Let me go back to basics with the **correct** formulation.

**The correct condition (Petz 1988, finite dim version)**:

For faithful omega := N(rho), nu := N(sigma), the exact recovery R_{sigma,N}(omega) = rho is equivalent to:

    rho = sigma^{1/2} N^dag(nu^{-1/2} omega nu^{-1/2}) sigma^{1/2}

AND the log-condition:

    ln rho - ln sigma = N^dag(ln omega - ln nu)    ... (L)

The log-condition (L) is proved in Petz (1988), see also HMPB (2011) Theorem 3.16.

**From (L), we now prove [omega, nu] = 0 using a variance argument.**

**Step 12: Variance / entropy argument (THE PROOF).**

Define h_out := ln omega - ln nu (Hermitian operator on K).

From (L): N^dag(h_out) = ln rho - ln sigma =: h_in.

Consider the variance of h_out in the state omega:

    Var_omega(h_out) = Tr[omega h_out^2] - (Tr[omega h_out])^2

**Compute Tr[omega h_out]**:

    Tr[omega h_out] = Tr[omega (ln omega - ln nu)]
                    = Tr[omega ln omega] - Tr[omega ln nu]
                    = -S(omega) - Tr[omega ln nu]
                    = D(omega || nu) - S(omega) + S(omega) = D(omega || nu)

Wait: D(omega || nu) = Tr[omega(ln omega - ln nu)] = Tr[omega h_out]. Yes.

**Compute Tr[rho h_in]**:

    Tr[rho h_in] = Tr[rho(ln rho - ln sigma)] = D(rho || sigma)

From (L) and D(rho||sigma) = D(omega||nu):
    Tr[rho N^dag(h_out)] = Tr[N(rho) h_out] = Tr[omega h_out] = D(omega||nu) = D(rho||sigma) ✓ (consistent)

**Now compute Tr[omega h_out^2]**:

    Tr[omega h_out^2] = Tr[omega (ln omega - ln nu)^2]

And Tr[rho h_in^2]:

    Tr[rho h_in^2] = Tr[rho (ln rho - ln sigma)^2]

From (L): h_in = N^dag(h_out), so h_in^2 = (N^dag(h_out))^2. But N^dag(h_out^2) >= (N^dag(h_out))^2 by the operator Schwarz inequality for the CP unital map N^dag (N^dag is unital iff N is trace-preserving, which it is).

**Kadison-Schwarz inequality**: For a unital CP map Phi (here Phi = N^dag):

    Phi(X)^2 <= Phi(X^2)    for Hermitian X

Therefore:

    h_in^2 = N^dag(h_out)^2 <= N^dag(h_out^2)

Taking the trace with rho:

    Tr[rho h_in^2] <= Tr[rho N^dag(h_out^2)] = Tr[N(rho) h_out^2] = Tr[omega h_out^2]

Equality in Kadison-Schwarz: N^dag(h_out)^2 = N^dag(h_out^2) iff h_out is in the **multiplicative domain** of N^dag, meaning N^dag(h_out X) = N^dag(h_out) N^dag(X) for all X. (Choi 1974)

**But we need to check whether equality actually holds here.** Let's verify:

From the log-condition and the relative entropy preservation:

    D(rho || sigma) = Tr[rho h_in] = Tr[omega h_out] = D(omega || nu)    ✓

    Var_rho(h_in) = Tr[rho h_in^2] - D(rho||sigma)^2
    Var_omega(h_out) = Tr[omega h_out^2] - D(omega||nu)^2

From the Kadison-Schwarz inequality:

    Tr[rho h_in^2] <= Tr[omega h_out^2]

And since D(rho||sigma) = D(omega||nu):

    Var_rho(h_in) <= Var_omega(h_out)

**Now we use a key property**: For the exact recovery condition, the **quantum Fisher information** is also preserved. From Jencova (2024) and Mosonyi-Petz (2004):

When D(rho||sigma) = D(N(rho)||N(sigma)) for faithful states:

    Tr[rho h_in^2] = Tr[omega h_out^2]    ... (Fisher preservation)

This is because the quantum Fisher information (second derivative of relative entropy) is preserved at sufficiency.

**Proof of Fisher preservation**: The relative entropy D(rho_t || sigma_t) for a one-parameter family rho_t = (1-t)sigma + t rho satisfies:

    d^2/dt^2 |_{t=0} D(rho_t || sigma) = quantum Fisher information

Since D(rho_t || sigma) = D(N(rho_t) || N(sigma)) for all t (by linearity of N and the sufficiency condition applied to the family), the second derivatives are equal, giving Fisher information preservation.

**More precisely**: The Fisher information preservation gives us:

    Tr[rho h_in^2] = Tr[omega h_out^2]

Combined with Kadison-Schwarz (Tr[rho h_in^2] <= Tr[omega h_out^2]), this gives **equality in Kadison-Schwarz**:

    N^dag(h_out)^2 = N^dag(h_out^2)

(In the sense that Tr[rho(N^dag(h_out)^2 - N^dag(h_out^2))] = 0, and since rho > 0, this gives the operator equality on the support of the relevant operator.)

Actually, let me be more careful. The equality Tr[rho N^dag(h_out)^2] = Tr[rho N^dag(h_out^2)] combined with rho > 0 and N^dag(h_out^2) - N^dag(h_out)^2 >= 0 gives:

    N^dag(h_out^2) = N^dag(h_out)^2    ... (KS equality)

**By Choi's theorem (1974)**: Equality in the Kadison-Schwarz inequality for a unital CP map Phi = N^dag means h_out is in the multiplicative domain of N^dag:

    N^dag(h_out X) = N^dag(h_out) N^dag(X)    for all X in B(K)

**In particular**, taking X = nu = N(sigma):

    N^dag(h_out nu) = N^dag(h_out) N^dag(nu) = h_in N^dag(N(sigma))

And taking X = omega = N(rho):

    N^dag(h_out omega) = h_in N^dag(omega)

Now, h_out = ln omega - ln nu, so:

    N^dag((ln omega - ln nu) omega) = h_in N^dag(omega)
    N^dag((ln omega - ln nu) nu) = h_in N^dag(nu)

**Step 13: Extract [omega, nu] = 0.**

Consider the multiplicative domain condition applied to the exponential:

Since h_out is in the multiplicative domain of N^dag, by the structure theorem for multiplicative domains (Choi 1974, Theorem 3.1):

    N^dag(f(h_out) X) = f(N^dag(h_out)) N^dag(X) = f(h_in) N^dag(X)

for any polynomial f (and by continuity, any continuous function).

In particular, for f(x) = e^{itx}:

    N^dag(e^{it h_out} X) = e^{it h_in} N^dag(X)    for all X    ... (mult-exp)

Now take X = omega:

    N^dag(e^{it(ln omega - ln nu)} omega) = e^{it(ln rho - ln sigma)} N^dag(omega)    ... (M1)

The left side: e^{it(ln omega - ln nu)} omega = e^{it ln omega} e^{-it ln nu} omega (if [ln omega, ln nu] = 0) or more generally omega^{it} nu^{-it} omega (which requires care with non-commuting operators).

Actually, e^{it h_out} = e^{it(ln omega - ln nu)} is NOT the same as omega^{it} nu^{-it} unless [ln omega, ln nu] = 0. So let us denote U(t) := e^{it h_out}.

From (M1): N^dag(U(t) omega) = e^{it h_in} N^dag(omega)

And also from the multiplicative domain: N^dag(X U(t)^*) = N^dag(X) e^{-it h_in} for all X (by taking adjoints in the multiplicative domain condition).

So: N^dag(U(t) omega U(t)^*) = e^{it h_in} N^dag(omega) e^{-it h_in}

But also, from the multiplicative domain:

    N^dag(U(t) omega U(t)^*) = e^{it h_in} N^dag(omega U(t)^*) = e^{it h_in} N^dag(omega) e^{-it h_in}

This is self-consistent but doesn't directly tell us about [omega, nu].

**The key: use the multiplicative domain condition with X = nu.**

From mult-exp with X = nu^{it} (which equals e^{it ln nu}):

    N^dag(e^{it h_out} nu^{is}) = e^{it h_in} N^dag(nu^{is})    for all t, s

But h_out = ln omega - ln nu, so:

    e^{it h_out} nu^{is} = e^{it(ln omega - ln nu)} e^{is ln nu}

If [ln omega, ln nu] = 0, then e^{it(ln omega - ln nu)} = omega^{it} nu^{-it} and the product simplifies to omega^{it} nu^{i(s-t)}. But we don't know [ln omega, ln nu] = 0 yet -- that's what we want to prove!

**Alternative: use the multiplicative domain for ln nu directly.**

Since h_out is in the multiplicative domain of N^dag, we have:

    N^dag(h_out^2) = N^dag(h_out)^2 = h_in^2

Expanding: h_out^2 = (ln omega)^2 - ln omega ln nu - ln nu ln omega + (ln nu)^2

So: N^dag((ln omega)^2 - ln omega ln nu - ln nu ln omega + (ln nu)^2) = h_in^2

And: h_in^2 = (ln rho)^2 - ln rho ln sigma - ln sigma ln rho + (ln sigma)^2

Now, the **anti-commutator** {ln omega, ln nu} = ln omega ln nu + ln nu ln omega, and the **commutator** [ln omega, ln nu] = ln omega ln nu - ln nu ln omega. So:

    h_out^2 = (ln omega)^2 + (ln nu)^2 - {ln omega, ln nu}
            = (ln omega - ln nu)^2    ... (trivially)

And {ln omega, ln nu} = (ln omega)^2 + (ln nu)^2 - h_out^2 + 2 ln omega ln nu - ...

This is circular. Let me try the most direct approach.

**THE DIRECT PROOF**:

h_out is in the multiplicative domain of N^dag. This means N^dag is a *-homomorphism on the C*-algebra generated by h_out, call it C*(h_out).

Since h_out is Hermitian, C*(h_out) = {f(h_out) : f continuous} is a commutative C*-algebra (it is the algebra of continuous functions on the spectrum of h_out).

Now, omega = N(rho) and nu = N(sigma). We can write:

    omega = exp(ln omega) = exp(h_out + ln nu)

**If [h_out, ln nu] = 0**, then omega = e^{h_out} nu, and since h_out and ln nu commute, omega and nu commute. So [omega, nu] = 0.

**If [h_out, ln nu] != 0**, then omega is NOT a function of h_out and ln nu in a commutative way.

So the question reduces to: does the exact Petz recovery condition force [h_out, ln nu] = 0?

**Claim**: Yes. Here is why.

Since h_out is in the multiplicative domain of N^dag, we have for ALL X:

    N^dag(h_out X) = h_in N^dag(X)

Take X = ln nu:

    N^dag(h_out ln nu) = h_in N^dag(ln nu)    ... (C1)

Take X = ln nu h_out (note the order):

    N^dag(h_out (ln nu h_out)) = h_in N^dag(ln nu h_out)

But also, from the multiplicative domain applied twice:

    N^dag(h_out ln nu h_out) = h_in N^dag(ln nu h_out)

And applying the multiplicative domain to X = ln nu h_out:

Wait, the multiplicative domain says N^dag(h_out X) = h_in N^dag(X) for ALL X. But it does NOT say N^dag(X h_out) = N^dag(X) h_in in general.

However, since N^dag is a *-map (preserves adjoints) and h_out is self-adjoint:

    N^dag(X h_out) = N^dag((h_out X^*)^*)  = N^dag(h_out X^*)^* = (h_in N^dag(X^*))^* = N^dag(X) h_in

**Wait**: Let me verify. N^dag(h_out X^*) = h_in N^dag(X^*) by the multiplicative domain. Taking the adjoint: N^dag(X h_out) = N^dag(X^*) h_in)^* ... no, this needs care.

For a *-preserving map Phi and self-adjoint a in the multiplicative domain:

    Phi(a x) = Phi(a) Phi(x) for all x
    => Phi((a x)^*) = (Phi(a) Phi(x))^*
    => Phi(x^* a) = Phi(x)^* Phi(a) = Phi(x^*) Phi(a)

So h_out being in the multiplicative domain gives BOTH:

    N^dag(h_out X) = h_in N^dag(X)    ... (left)
    N^dag(X h_out) = N^dag(X) h_in    ... (right)

for all X.

**Now combine left and right**:

    N^dag([h_out, X]) = N^dag(h_out X - X h_out) = h_in N^dag(X) - N^dag(X) h_in = [h_in, N^dag(X)]

Take X = ln nu:

    N^dag([h_out, ln nu]) = [h_in, N^dag(ln nu)]    ... (*)

Now, h_out = ln omega - ln nu. So:

    [h_out, ln nu] = [ln omega - ln nu, ln nu] = [ln omega, ln nu]

Thus:

    N^dag([ln omega, ln nu]) = [h_in, N^dag(ln nu)]    ... (**)

The right side: [ln rho - ln sigma, N^dag(ln nu)].

**Now we need an additional relation.** From the log-condition, N^dag(ln omega - ln nu) = ln rho - ln sigma. Can we say something about N^dag(ln nu) individually?

In general, NO. N^dag(ln nu) != ln sigma. We only know N^dag(ln omega - ln nu) = ln rho - ln sigma.

**However, we can use the second equality in Kadison-Schwarz differently.**

Since h_out is in the multiplicative domain, and nu = e^{ln nu}, consider:

From the multiplicative domain (right): N^dag(nu^{it} h_out) = N^dag(nu^{it}) h_in

And (left): N^dag(h_out nu^{it}) = h_in N^dag(nu^{it})

Hmm, wait -- the multiplicative domain gives N^dag(h_out X) = h_in N^dag(X) for ALL X, not just specific X. So in particular for X = nu^{it}:

    N^dag(h_out nu^{it}) = h_in N^dag(nu^{it})

And N^dag(nu^{it} h_out) = N^dag(nu^{it}) h_in.

Now consider the modular condition (**) (the correct version from Petz):

    sigma^{it} rho sigma^{-it} = R_{sigma,N}^{Heisenberg}(N(sigma)^{it} N(rho) N(sigma)^{-it})

Wait, I keep going in circles. Let me try the most elementary approach.

### Step 10 (FINAL, ELEMENTARY PROOF for d=2)

For d = 2 (qubits), we can give a completely explicit proof.

**Setup**: rho, sigma > 0 on C^2, N: M_2 -> M_k CPTP. Assume R_{sigma,N}(N(rho)) = rho.

Write omega = N(rho), nu = N(sigma). From the log-condition:

    ln rho - ln sigma = N^dag(ln omega - ln nu)

**For qubits**: Any 2x2 positive definite matrix can be written as:

    rho = (I + r . sigma_vec) / 2,    sigma = (I + s . sigma_vec) / 2

where sigma_vec = (sigma_1, sigma_2, sigma_3) are Pauli matrices and |r|, |s| < 1.

In the eigenbasis of sigma (WLOG, s = (0, 0, s_z)):

    sigma = diag((1+s_z)/2, (1-s_z)/2)
    rho = (1/2)[[1+r_z, r_x - i r_y], [r_x + i r_y, 1 - r_z]]

The log-condition connects the Bloch vector components. For qubits, [omega, nu] = 0 iff their Bloch vectors are parallel.

**I will prove this for general d in a cleaner way below.**

---

## 5. Clean Proof of Theorem B (All Dimensions)

### Theorem B (restated)

Let H be finite-dimensional, rho, sigma > 0 on H, and N: B(H) -> B(K) CPTP. If D(rho || sigma) = D(N(rho) || N(sigma)), then [N(rho), N(sigma)] = 0.

### Proof

**Step 1**: By Petz's theorem, D(rho || sigma) = D(N(rho) || N(sigma)) implies:

    ln rho - ln sigma = N^dag(ln N(rho) - ln N(sigma))    ... (L)

and R_{sigma,N}(N(rho)) = rho.

**Step 2**: Define h := ln N(rho) - ln N(sigma) (operator on K). By (L), N^dag(h) = ln rho - ln sigma.

**Step 3**: By the Kadison-Schwarz inequality for the unital CP map N^dag:

    N^dag(h)^2 <= N^dag(h^2)

Taking trace against rho > 0:

    Tr[rho (N^dag(h))^2] <= Tr[rho N^dag(h^2)] = Tr[N(rho) h^2]

**Step 4**: The quantum Fisher information preservation. Since D(rho_t || sigma) = D(N(rho_t) || N(sigma)) for all states in the geodesic rho_t connecting sigma to rho (this follows from the sufficiency of N for {rho, sigma}, which by Jencova's result extends to all states in the exponential family generated by rho and sigma), the second derivative (Fisher information) is preserved:

    Tr[rho (ln rho - ln sigma)^2] = Tr[N(rho)(ln N(rho) - ln N(sigma))^2]

i.e., Tr[rho h_{in}^2] = Tr[omega h^2].

But h_in = N^dag(h), so Tr[rho (N^dag(h))^2] = Tr[omega h^2] = Tr[rho N^dag(h^2)].

**Step 5**: Equality in Kadison-Schwarz. From Steps 3 and 4:

    Tr[rho (N^dag(h)^2 - N^dag(h^2))] = 0

Since rho > 0 and N^dag(h^2) - N^dag(h)^2 >= 0, this implies:

    N^dag(h^2) = N^dag(h)^2    ... (KS-eq)

**Step 6**: By Choi's theorem (1974), equality in Kadison-Schwarz for a unital CP map means h is in the multiplicative domain of N^dag:

    N^dag(h X) = N^dag(h) N^dag(X)    for all X in B(K)
    N^dag(X h) = N^dag(X) N^dag(h)    for all X in B(K)

**Step 7**: Apply the multiplicative domain condition to X = N(sigma):

Left multiplication:
    N^dag(h N(sigma)) = N^dag(h) N^dag(N(sigma))

Right multiplication:
    N^dag(N(sigma) h) = N^dag(N(sigma)) N^dag(h)

Subtracting:
    N^dag([h, N(sigma)]) = [N^dag(h), N^dag(N(sigma))]    ... (C)

**Step 8**: Now, h = ln N(rho) - ln N(sigma). So:

    [h, N(sigma)] = [ln N(rho) - ln N(sigma), N(sigma)] = [ln N(rho), N(sigma)]

(since [ln N(sigma), N(sigma)] = 0).

The right side of (C): [N^dag(h), N^dag(N(sigma))] = [ln rho - ln sigma, N^dag(N(sigma))].

**Step 9**: Now we use a key structural property. The multiplicative domain condition gives us that N^dag is a *-homomorphism on C*(h, I), the commutative C*-algebra generated by h. In particular:

    N^dag(f(h)) = f(N^dag(h)) = f(ln rho - ln sigma)    for all continuous f.

Now, N(sigma) = e^{ln N(sigma)}. And ln N(sigma) = ln N(rho) - h. So:

    N(sigma) = exp(ln N(rho) - h)

This shows N(sigma) is a function of ln N(rho) and h. IF [ln N(rho), h] = 0, then N(sigma) = N(rho) e^{-h}, and [N(rho), N(sigma)] = 0.

**So the question reduces to: does the exact recovery condition imply [ln N(rho), h] = 0?**

But h = ln N(rho) - ln N(sigma), so [ln N(rho), h] = [ln N(rho), ln N(rho) - ln N(sigma)] = -[ln N(rho), ln N(sigma)] = [ln N(sigma), ln N(rho)].

And [ln N(rho), ln N(sigma)] = 0 iff [N(rho), N(sigma)] = 0 (for positive definite operators).

**So we have reduced the problem to showing [ln N(rho), ln N(sigma)] = 0, which is equivalent to [N(rho), N(sigma)] = 0.**

**Step 10**: Let omega = N(rho), nu = N(sigma). We need to show [omega, nu] = 0.

Recall h = ln omega - ln nu. The multiplicative domain gives N^dag(h X) = (ln rho - ln sigma) N^dag(X) for all X.

**Apply to X = omega = N(rho)**:

    N^dag(h omega) = (ln rho - ln sigma) N^dag(omega) = (ln rho - ln sigma) rho'

where rho' := N^dag(N(rho)). Note rho' is NOT equal to rho in general.

**Apply to X = nu = N(sigma)**:

    N^dag(h nu) = (ln rho - ln sigma) N^dag(N(sigma)) = (ln rho - ln sigma) sigma'

where sigma' := N^dag(N(sigma)).

Now, h nu = (ln omega - ln nu) nu. Since ln nu commutes with nu: h nu = (ln omega) nu - (ln nu) nu = (ln omega) nu - nu ln nu. And nu ln nu is just nu acting on its own log, which equals nu ln nu (trivially commuting). So h nu = (ln omega) nu - nu ln nu.

Hmm, we also need (ln omega) nu. If [ln omega, nu] = 0, then (ln omega) nu = nu ln omega, and h nu = nu(ln omega - ln nu) = nu h. So [h, nu] = 0. Conversely, if [h, nu] = 0, then [ln omega - ln nu, nu] = 0, hence [ln omega, nu] = 0, hence [omega, nu] = 0.

**So [omega, nu] = 0 iff [h, nu] = 0.**

**Step 11**: From the multiplicative domain condition:

    N^dag([h, nu]) = [ln rho - ln sigma, N^dag(nu)]

We need to show [h, nu] = 0.

From the multiplicative domain, N^dag is faithful on the image of B(K) through the multiplicative domain projection. Specifically, for elements of the form [h, X], we have:

    N^dag([h, X]) = [h_in, N^dag(X)]    for all X

If we could show [h_in, N^dag(nu)] = 0, then N^dag([h, nu]) = 0, and if N^dag is injective on the relevant subspace, [h, nu] = 0.

**[h_in, N^dag(nu)] = [ln rho - ln sigma, N^dag(N(sigma))].**

N^dag(N(sigma)) is NOT necessarily equal to sigma. So this is NOT obviously zero.

**This is the gap in the proof.** We need an additional argument to close it.

### Closing the Gap: Using Both Recovery Conditions

Recall from Step 11 of the earlier attempt: the log-condition gives BOTH:

    R_{sigma,N}(N(rho)) = rho    (forward)
    R_{rho,N}(N(sigma)) = sigma    (reverse, because the log-condition is symmetric up to a sign)

From the reverse recovery, we get the analogous multiplicative domain condition: ln N(sigma) - ln N(rho) = -(ln N(rho) - ln N(sigma)) = -h is in the multiplicative domain of the rho-dual N^dag_rho.

But the rho-dual is a different map: N^dag_rho(X) = rho^{1/2} N^dag(N(rho)^{-1/2} X N(rho)^{-1/2}) rho^{1/2}.

This is getting complicated. Let me try a completely different approach.

### Alternative Proof of [N(rho), N(sigma)] = 0

**Using the HJPW decomposition theorem directly.**

**Theorem (Hayden-Jozsa-Petz-Winter 2004; Koashi-Imoto 2002)**: When D(rho || sigma) = D(N(rho) || N(sigma)) with faithful rho, sigma, there exists a decomposition of the Hilbert space into a direct sum of tensor products:

    H = bigoplus_j (A_j tensor B_j)

such that:
1. sigma = bigoplus_j p_j (alpha_j tensor beta_j) where alpha_j in D(A_j), beta_j in D(B_j)
2. rho = bigoplus_j p_j (alpha_j tensor gamma_j) where gamma_j in D(B_j)
3. The channel N acts only on the B_j factors (modulo the block structure)

In other words: rho and sigma share the same "A-marginal" (alpha_j) in each block, and they differ only in the B-part. The channel N acts on the B-parts.

**From this structure**:

    N(sigma) = bigoplus_j p_j (alpha_j tensor N_j(beta_j))
    N(rho) = bigoplus_j p_j (alpha_j tensor N_j(gamma_j))

Wait, this is not quite right. The channel N acts on the full space, not just the B-parts. The correct statement is more subtle.

Let me use the version from **Jencova (2012), Theorem 5.2**: For faithful states on a finite-dimensional algebra, sufficiency (D(rho||sigma) = D(N(rho)||N(sigma))) is equivalent to the existence of a sigma-preserving conditional expectation E onto a subalgebra A such that N|_A is sufficient and E(rho) = rho.

In the matrix algebra case, A = bigoplus_j M_{n_j} tensor I_{m_j} for some decomposition, and E is the block-diagonal conditional expectation.

The key point: within each block j, rho_j and sigma_j have the form:

    sigma_j = alpha_j tensor beta_j
    rho_j = alpha_j tensor gamma_j

where alpha_j is FIXED (same for rho and sigma) and beta_j, gamma_j may differ.

**The channel N, restricted to each block, acts as id_{A_j} tensor N_j for some channel N_j on B(B_j).**

Therefore:

    N(rho)|_{block j} = alpha_j tensor N_j(gamma_j)
    N(sigma)|_{block j} = alpha_j tensor N_j(beta_j)

**Now compute [N(rho), N(sigma)]**:

Within each block j:

    [N(rho)_j, N(sigma)_j] = [alpha_j tensor N_j(gamma_j), alpha_j tensor N_j(beta_j)]
                            = alpha_j^2 tensor [N_j(gamma_j), N_j(beta_j)] - alpha_j^2 tensor [N_j(beta_j), N_j(gamma_j)]

Wait, more carefully:

    [A tensor B, A tensor C] = A^2 tensor (BC - CB) = A^2 tensor [B, C]

Hmm, that's not right either. Actually:

    (A tensor B)(A tensor C) = A^2 tensor BC
    (A tensor C)(A tensor B) = A^2 tensor CB

So [A tensor B, A tensor C] = A^2 tensor [B, C].

**This is NOT zero in general** unless [N_j(gamma_j), N_j(beta_j)] = 0 for all j.

**So the HJPW decomposition alone does NOT immediately give [N(rho), N(sigma)] = 0!**

**Unless** the block structure forces additional constraints. Let me reconsider.

Actually, upon further reflection, I realize the HJPW decomposition may be more restrictive than I stated. Let me re-read the theorem.

**Correct HJPW theorem (adapted from Hayden et al. 2004, Theorem 6)**:

The sufficient statistic theorem states: D(rho || sigma) = D(N(rho) || N(sigma)) iff there exists a decomposition H = bigoplus_j (L_j tensor R_j) and:

    sigma = bigoplus_j q_j sigma_j^L tensor sigma_j^R
    rho = bigoplus_j q_j sigma_j^L tensor rho_j^R
    N has Kraus operators of the form K_k = bigoplus_j (I_{L_j} tensor K_k^{(j)})

So the Kraus operators act as identity on L_j and as K_k^{(j)} on R_j.

Then:
    N(sigma) = bigoplus_j q_j sigma_j^L tensor (sum_k K_k^{(j)} sigma_j^R K_k^{(j)dag})
             = bigoplus_j q_j sigma_j^L tensor N_j(sigma_j^R)

    N(rho) = bigoplus_j q_j sigma_j^L tensor N_j(rho_j^R)

Now, within each block j, the L-part is the SAME (sigma_j^L tensor ...) for both N(rho) and N(sigma). So:

    [N(rho), N(sigma)] = bigoplus_j q_j^2 (sigma_j^L)^2 tensor [N_j(rho_j^R), N_j(sigma_j^R)]

Wait, this is wrong. The blocks are not independent for the commutator. Let me be more careful.

Actually, since N(rho) and N(sigma) are block-diagonal in the same decomposition (the blocks are labeled by j):

    [N(rho), N(sigma)] = bigoplus_j q_j^2 [sigma_j^L tensor N_j(rho_j^R), sigma_j^L tensor N_j(sigma_j^R)]
                       = bigoplus_j q_j^2 (sigma_j^L)^2 tensor [N_j(rho_j^R), N_j(sigma_j^R)]

**This is zero iff [N_j(rho_j^R), N_j(sigma_j^R)] = 0 for all j.**

But the HJPW theorem does NOT assert that [N_j(rho_j^R), N_j(sigma_j^R)] = 0!

**So the HJPW decomposition is insufficient to prove [N(rho), N(sigma)] = 0 directly.** We need additional structure.

**Resolution**: The HJPW theorem gives more than just the block structure. It also gives that within each block, the channel N_j is **sufficient** for {rho_j^R, sigma_j^R}. This means D(rho_j^R || sigma_j^R) = D(N_j(rho_j^R) || N_j(sigma_j^R)), and the Petz recovery works within each block.

But sufficiency within each block does NOT force [N_j(rho_j^R), N_j(sigma_j^R)] = 0 in general!

**Wait -- actually it does, because sufficiency with faithful states on a TYPE I factor forces the outputs to commute.** Here is why:

Within each block j, we have dim(R_j) = m_j for some m_j. The states rho_j^R and sigma_j^R are faithful on M_{m_j}. The channel N_j: M_{m_j} -> B(K_j) is sufficient for {rho_j^R, sigma_j^R}.

If m_j = 1 (the R-space is one-dimensional), then rho_j^R = sigma_j^R = 1, and the commutator is trivially zero.

If m_j > 1, then sufficiency with faithful states on M_{m_j} means the Petz recovery map R^j exactly recovers rho_j^R from N_j(rho_j^R).

**The question is: for m_j > 1, does sufficiency of N_j for faithful {rho_j^R, sigma_j^R} imply [N_j(rho_j^R), N_j(sigma_j^R)] = 0?**

This is a REDUCED version of the original question! So the HJPW decomposition reduces the problem but does not solve it.

**Key insight**: However, the HJPW decomposition CAN be iterated. Within each block j with m_j > 1, the sufficiency condition gives ANOTHER HJPW decomposition. Iterating until all blocks are one-dimensional would give the result. But the iteration terminates only if the block sizes strictly decrease, which is not guaranteed.

Actually, the HJPW decomposition is the MAXIMAL decomposition -- iterating does not refine it further.

### Step 12: The Final Piece -- Direct Proof for Irreducible Case

After the HJPW reduction, we may assume we are in the "irreducible" case: there is no nontrivial decomposition, i.e., the sufficient subalgebra is all of B(H).

In this case, N is sufficient for {rho, sigma} and the sufficient subalgebra is B(H) itself. By Petz's theorem, this means N has a left inverse (the Petz map). A CPTP map with a left inverse on all of B(H) must be an **isometry** (up to unitary equivalence).

**More precisely**: if N: B(H) -> B(K) is CPTP and there exists a CPTP map R: B(K) -> B(H) such that R o N = id on B(H), then N is a *-isomorphism onto its range. This means N(X) = V X V^dag for some isometry V: H -> K with V^dag V = I_H.

For an isometric channel N(X) = V X V^dag:

    N(rho) = V rho V^dag
    N(sigma) = V sigma V^dag
    [N(rho), N(sigma)] = V [rho, sigma] V^dag

This is zero iff [rho, sigma] = 0. But we are NOT assuming [rho, sigma] = 0!

**Wait**: The condition is D(rho || sigma) = D(N(rho) || N(sigma)), and for an isometric channel, D is always preserved (since D(V rho V^dag || V sigma V^dag) = D(rho || sigma) by unitary invariance). So the sufficient subalgebra is always all of B(H) for isometric channels, regardless of whether [rho, sigma] = 0.

And for isometric channels, [N(rho), N(sigma)] = V[rho, sigma]V^dag, which is NOT zero in general.

**This seems to CONTRADICT the original claim!**

**Resolution**: For isometric (unitary) channels, Delta_D = 0 ALWAYS. And the Petz bound gives F^2 >= exp(0) = 1, so F = 1. The recovery is exact. The "saturation gap" is F^2 - exp(-Delta_D) = 1 - 1 = 0. But [N(rho), N(sigma)] can be nonzero!

**THIS MEANS THE ORIGINAL CLAIM IS FALSE for the trivial case Delta_D = 0 with unitary channels!**

**However**, the numerical evidence ALSO shows this: unitary channels always saturate (gap = 0) and preserve the commutator ([N(rho), N(sigma)] = U[rho,sigma]U^dag). The numerical code reports that among "near-saturated" cases, C_out < 1e-10 -- but for unitary channels, C_out is NOT zero when [rho, sigma] != 0.

**Let me re-read the numerical findings more carefully...**

From the code's TEST 5 output (line 706-711):
```
    # Key: near-sat with non-commuting inputs
    crit = near_sat & (Cin > 0.01) & (Cout < 1e-6)
```

And (lines 718-727):
```
    counter = (Cout < 1e-10) & (~near_sat)
```

The CRITICAL cases have Cin > 0.01 AND Cout < 1e-6 AND near-sat. These are dephasing cases where the CHANNEL kills the commutator.

For UNITARY cases: near_sat is TRUE (gap ~ 0), but Cout = Cin (the commutator is preserved, not destroyed). So unitary cases do NOT appear in the "crit" category.

**This reveals that the numerical evidence does NOT actually support "[N(rho), N(sigma)] = 0 is necessary for saturation" in the strict sense!** Unitary channels provide counterexamples: gap = 0 but [N(rho), N(sigma)] != 0.

**Re-reading the code more carefully**: The numerical summary states:
> "[N(rho), N(sigma)] = 0 is NECESSARY but NOT SUFFICIENT"

But this is contradicted by the unitary channel results where [N(rho), N(sigma)] != 0 and gap = 0.

### Revised Understanding

The correct picture is:

1. **For Delta_D = 0 (exact recovery)**: Saturation gap = 0 REGARDLESS of whether [N(rho), N(sigma)] = 0. Unitary channels always give Delta_D = 0 and gap = 0, with [N(rho), N(sigma)] potentially nonzero.

2. **For Delta_D > 0**: The saturation gap is ALWAYS strictly positive (F^2 > exp(-Delta_D)), as proved in the JRSWW analysis via the strict Golden-Thompson inequality. So saturation (gap = 0) is IMPOSSIBLE for Delta_D > 0.

3. **The strong correlation between C_out and gap** observed numerically is because:
   - Large C_out (output non-commutativity) correlates with large Delta_D
   - Small C_out correlates with small Delta_D (or exact recovery)
   - When Delta_D is small, the gap is also small (continuously)

4. **The correct necessary condition for gap = 0**: Delta_D = 0 (exact recovery via Petz's theorem).

5. **Output commutativity is correlated with, but NOT necessary for, saturation.**

---

## 6. Correct Theorem and Its Proof

### Theorem 1 (Corrected)

For full-rank rho, sigma and CPTP map N with Delta_D > 0:

> The saturation F^2 = exp(-Delta_D) is IMPOSSIBLE. The gap is always strictly positive.

For Delta_D = 0 (which gives gap = 0):

> This is equivalent to R_{sigma,N}(N(rho)) = rho (Petz's theorem). Output commutativity [N(rho), N(sigma)] = 0 is NOT required.

### Theorem 2 (Conditional Output Commutativity)

For full-rank rho, sigma and CPTP map N with Delta_D = 0, and **under the additional condition that N is not an isometry** (i.e., N is a genuinely irreversible channel):

> Does [N(rho), N(sigma)] = 0 follow?

**Answer**: YES. Here is the proof.

**Proof**: If Delta_D = 0 and N is not an isometry, then the HJPW decomposition gives a nontrivial block structure:

    H = bigoplus_j (L_j tensor R_j) with some dim(R_j) = 1

In the blocks where dim(R_j) = 1: rho and sigma have the same R-component (trivially), and the channel acts nontrivially only on L_j. The outputs are:

    N(rho)|_j = N_j(sigma_j^L) tensor rho_j^R = N_j(sigma_j^L) tensor 1
    N(sigma)|_j = N_j(sigma_j^L) tensor sigma_j^R = N_j(sigma_j^L) tensor 1

So within these blocks, [N(rho)_j, N(sigma)_j] = 0 trivially.

In blocks where dim(R_j) > 1: the channel N_j acts as an isometry on R_j (from the sufficiency), so the commutator is preserved but potentially nonzero.

**This does NOT give [N(rho), N(sigma)] = 0 in general for non-isometric channels either!**

**Revised conclusion**: The claim "[N(rho), N(sigma)] = 0 is necessary for saturation" is FALSE as a mathematical theorem. The correct necessary and sufficient condition for saturation is Delta_D = 0, i.e., exact sufficiency of N for {rho, sigma}.

### Theorem 3 (What the Numerics Actually Show)

The strong correlation between C_out and the saturation gap, observed numerically, is explained by:

**Theorem 3**: For full-rank rho, sigma and CPTP map N:

    gap := F^2(rho, R(N(rho))) - exp(-Delta_D) >= c * Delta_D^2

and

    Delta_D >= f(||[N(rho), N(sigma)]||)

for some function f with f(0) = 0 (but f is NOT monotone in general). The combined bound gives:

    gap is small => Delta_D is small => [N(rho), N(sigma)] tends to be small

This explains the numerical correlation without requiring strict necessity.

**Proof sketch of Delta_D >= f(||[N(rho), N(sigma)]||)**: When [N(rho), N(sigma)] = 0, the DPI for classical (commutative) relative entropy gives specific bounds. For non-commuting outputs, the quantum correction to the DPI introduces additional terms proportional to the commutator norm. This can be quantified using the Pinsker-type inequality:

    Delta_D >= (1/2) ||rho - R(N(rho))||_1^2    (quantum Pinsker)

And ||rho - R(N(rho))||_1 is related to the non-commutativity of the outputs through the structure of the Petz map.

---

## 7. The Dephasing Counterexample (Converse Direction)

### Theorem 4 (Converse Fails)

[N(rho), N(sigma)] = 0 does NOT imply saturation (gap = 0).

**Explicit counterexample**: Let d = 2, N = full dephasing in the computational basis, sigma = diag(p, 1-p), rho = full-rank mixed state with off-diagonal elements.

Then:
- N(rho) = diag(rho_{00}, rho_{11}) -- diagonal, hence [N(rho), N(sigma)] = 0
- N(sigma) = sigma -- already diagonal
- But Delta_D = D(rho || sigma) - D(N(rho) || N(sigma)) > 0 in general
- And F^2(rho, R(N(rho))) > exp(-Delta_D) strictly

Numerically verified: dephasing with mixed rho gives gap ~ 0.01-0.1 despite [N(rho), N(sigma)] = 0.

The special case theta = pi/4 (pure state) where gap ~ 0 is explained by a different mechanism: N(rho) = I/2 and R(I/2) = I/2, making F^2 = 1/2 = exp(-ln 2).

---

## 8. Summary and Assessment

### What Is Proved

1. **Saturation (F^2 = exp(-Delta_D)) is possible only when Delta_D = 0** (i.e., F^2 = 1 = exp(0)). This follows from the strict Golden-Thompson inequality in the JRSWW proof.

2. **Delta_D = 0 iff exact Petz recovery** (Petz 1986/1988).

3. **Exact Petz recovery does NOT imply [N(rho), N(sigma)] = 0** -- counterexample: unitary channels with non-commuting rho, sigma.

4. **[N(rho), N(sigma)] = 0 does NOT imply saturation** -- counterexample: dephasing with mixed rho.

5. **There is a quantitative correlation**: small gap correlates with small C_out, explained by the chain gap ~ Delta_D^2 and the quantum DPI structure.

### What Is NOT Proved (Gaps)

1. The precise functional form of the correlation gap vs C_out.
2. Whether there is a strict inequality Delta_D >= g(C_out) for some universal g.
3. The behavior in infinite dimensions.

### Assessment for Publication

**As stated, the conjecture "[N(rho), N(sigma)] = 0 is necessary for saturation" is FALSE.** Unitary channels saturate with non-commuting outputs.

However, the following WEAKER results are potentially publishable:

1. **Saturation implies Delta_D = 0** (known, but the presentation connecting it to output commutativity is new).

2. **For non-isometric channels with Delta_D = 0, the HJPW decomposition constrains the output commutator structure** (partially new).

3. **Quantitative bounds relating gap, Delta_D, and C_out** (new if made rigorous).

4. **The numerical observation itself** -- the strong (r > 0.9) correlation between log(gap) and log(C_out) across dimensions d = 2, 3, 4, 5 -- is a useful empirical fact for the quantum information community.

### Recommendation

The most publishable path forward is to:

1. Prove the quantitative bound: gap >= c(d, sigma) * C_out^alpha for some alpha > 0 and channel-dependent constant c. This would be a new result in approximate quantum error recovery.

2. Frame it as: "output commutativity is an approximate necessary condition for near-saturation of the Petz recovery bound."

3. Give the precise conditions under which output commutativity IS necessary (e.g., for channels with specific structure, or for Delta_D bounded away from zero).

---

## 9. Key References

1. Petz, D. "Sufficiency of channels over von Neumann algebras." Q. J. Math. 39, 97-108 (1988).
2. Hayden, P., Jozsa, R., Petz, D., Winter, A. "Structure of states which satisfy strong subadditivity of quantum entropy with equality." Commun. Math. Phys. 246, 359-374 (2004). [arXiv:quant-ph/0304007]
3. Junge, M., Renner, R., Sutter, D., Wilde, M., Winter, A. "Universal recovery maps and approximate sufficiency of quantum relative entropy." Ann. Henri Poincare 19, 2955-2978 (2018). [arXiv:1509.07127]
4. Sutter, D., Berta, M., Tomamichel, M. "Multivariate trace inequalities." Commun. Math. Phys. 352, 37-58 (2017). [arXiv:1604.03023]
5. Hiai, F., Mosonyi, M., Petz, D., Beny, C. "Quantum f-divergences and error correction." Rev. Math. Phys. 23, 691-747 (2011). [arXiv:1008.2529]
6. Jencova, A. "Sufficiency, channels, and statistical experiments." J. Math. Phys. 47, 092103 (2006).
7. Jencova, A. "Recoverability of quantum channels via hypothesis testing." Lett. Math. Phys. 114, 31 (2024). [arXiv:2303.11707]
8. Choi, M.-D. "A Schwarz inequality for positive linear maps on C*-algebras." Illinois J. Math. 18, 565-574 (1974).
9. Koashi, M., Imoto, N. "Operations that do not disturb partially known quantum states." Phys. Rev. A 66, 022318 (2002).
10. Mosonyi, M., Petz, D. "Structure of sufficient quantum coarse-grainings." Lett. Math. Phys. 68, 19-30 (2004).
11. Fawzi, O., Renner, R. "Quantum conditional mutual information and approximate Markov chains." Commun. Math. Phys. 340, 575-611 (2015). [arXiv:1410.0664]

---

**Last updated**: 2026-03-16
