# Structural Inequalities from the tau = 1 - F Framework

## Mathematical Setting

**Definition.** For a quantum channel N: A -> B, input state rho on A, and reference state sigma on A:

    tau(rho, N, sigma) := 1 - F(rho, R_sigma o N(rho))

where R_sigma is the Petz recovery map for (N, sigma), and F is the Uhlmann fidelity.

**Known bounds:**
- tau <= 1 - exp(-Sigma/2) where Sigma is entropy production (from Fawzi-Renner + Crooks)
- tau = 0 iff Sigma = 0 iff quantum Markov chain iff perfect recovery
- The Petz map is the unique retrodiction functor (Parzygnat-Buscemi 2023)

---

## Direction 1: Composition Law for Sequential Channels

### Setup

Let rho be a state on A, with channels N1: A -> B and N2: B -> C.
Define:
- tau_1 := tau(rho, N1, sigma) = 1 - F_1
- tau_2 := tau(N1(rho), N2, N1(sigma)) = 1 - F_2
- tau_{12} := tau(rho, N2 o N1, sigma) = 1 - F_{12}

where F_i are the respective recovery fidelities.

### Key Ingredient: Petz Map Functoriality

From Parzygnat-Buscemi (2023), the Petz map is the unique retrodiction functor:

    R_sigma(N2 o N1) = R_sigma,N1 o R_{N1(sigma),N2}

That is, the Petz recovery map for the composition equals the composition of individual Petz maps in reverse order. This is exactly the functorial (contravariant) property.

### Derivation

**Step 1.** The recovery fidelity for the composed channel is:

    F_{12} = F(rho, R_sigma(N2 o N1) o (N2 o N1)(rho))
           = F(rho, R_{sigma,N1} o R_{N1(sigma),N2} o N2 o N1(rho))

**Step 2.** Now consider applying R_{N1(sigma),N2} o N2 to the state N1(rho). By definition:

    F(N1(rho), R_{N1(sigma),N2} o N2 o N1(rho)) = F_2

So R_{N1(sigma),N2} o N2(N1(rho)) has fidelity F_2 with N1(rho).

**Step 3.** Now apply R_{sigma,N1} to both sides. By the data processing inequality for fidelity (F is non-decreasing under CPTP maps):

    F(R_{sigma,N1}(N1(rho)), R_{sigma,N1} o R_{N1(sigma),N2} o N2 o N1(rho)) >= F_2

But R_{sigma,N1}(N1(rho)) is the Petz-recovered state from step 1, with:

    F(rho, R_{sigma,N1}(N1(rho))) = F_1

**Step 4.** By the triangle inequality for the Bures metric (or Bures angle):

    d_B(rho, R_{12}(N2 o N1)(rho)) <= d_B(rho, R_1(N1(rho))) + d_B(R_1(N1(rho)), R_{12}(N2 o N1)(rho))

where d_B is the Bures distance, d_B(rho, sigma) = sqrt(2(1 - F(rho, sigma))).

**Step 5.** Using d_B = sqrt(2(1-F)) = sqrt(2*tau):

    sqrt(2*tau_{12}) <= sqrt(2*tau_1) + sqrt(2*tau_2)

Therefore:

### Result 1: Composition Inequality (Bures Triangle)

    sqrt(tau_{12}) <= sqrt(tau_1) + sqrt(tau_2)

Or equivalently:

    tau_{12} <= tau_1 + tau_2 + 2*sqrt(tau_1 * tau_2)
             = (sqrt(tau_1) + sqrt(tau_2))^2

**Physical interpretation:** The temporal asymmetry of a composed process is bounded by the sum of individual asymmetries *measured in the Bures metric*. The sqrt-metric is the natural one.

**But wait -- is Step 3 valid?** Actually, we need to be more careful. Let me redo this properly.

### Rigorous Derivation via Fidelity Submultiplicativity

**Step 1 (clean).** Define:
- rho_0 = rho (initial state)
- rho_1 = N1(rho) (after first channel)
- rho_2 = N2(N1(rho)) (after both channels)
- rho_1' = R_{N1(sigma),N2}(rho_2) (partial recovery: undo N2)
- rho_0' = R_{sigma,N1}(rho_1') (full recovery: undo N1 after undoing N2)

By functoriality, rho_0' = R_{sigma}(N2 o N1)(rho_2).

**Step 2.** F_{12} = F(rho_0, rho_0') by definition.

**Step 3.** By data processing inequality applied to R_{sigma,N1}:

    F(rho_0, rho_0') = F(rho_0, R_{sigma,N1}(rho_1'))

Now, we want to relate this to F_1 and F_2.

**Step 4.** We know F_2 = F(rho_1, rho_1') and F_1 = F(rho_0, R_{sigma,N1}(rho_1)).

**Step 5.** By the Bures metric triangle inequality on the triple (rho_0, R_{sigma,N1}(rho_1), R_{sigma,N1}(rho_1')):

    d_B(rho_0, R_{sigma,N1}(rho_1')) <= d_B(rho_0, R_{sigma,N1}(rho_1)) + d_B(R_{sigma,N1}(rho_1), R_{sigma,N1}(rho_1'))

**Step 6.** The first term gives d_B(rho_0, R_{sigma,N1}(rho_1)) = sqrt(2*tau_1).

**Step 7.** For the second term, by data processing under R_{sigma,N1}:

    d_B(R_{sigma,N1}(rho_1), R_{sigma,N1}(rho_1')) <= d_B(rho_1, rho_1') = sqrt(2*tau_2)

Wait, this inequality goes the WRONG way. Data processing says fidelity increases (distance decreases) under CPTP maps. So:

    F(R_{sigma,N1}(rho_1), R_{sigma,N1}(rho_1')) >= F(rho_1, rho_1') = F_2

which means:

    d_B(R_{sigma,N1}(rho_1), R_{sigma,N1}(rho_1')) <= d_B(rho_1, rho_1') = sqrt(2*tau_2)

Yes, this IS the correct direction! CPTP maps are contractions in Bures distance.

**Step 8.** Combining:

    sqrt(2*tau_{12}) = d_B(rho_0, rho_0') <= sqrt(2*tau_1) + sqrt(2*tau_2)

### Theorem 1 (Composition Inequality for tau)

For sequential channels N1: A -> B, N2: B -> C:

    sqrt(tau(rho, N2 o N1, sigma)) <= sqrt(tau(rho, N1, sigma)) + sqrt(tau(N1(rho), N2, N1(sigma)))

Equivalently:

    tau_{12} <= (sqrt(tau_1) + sqrt(tau_2))^2 = tau_1 + tau_2 + 2*sqrt(tau_1*tau_2)

**Corollary 1.1** (Weak sub-additivity): tau_{12} <= 4*max(tau_1, tau_2).

**Corollary 1.2** (Near-perfect regime): If tau_1, tau_2 << 1, then tau_{12} ~ tau_1 + tau_2 (to leading order).

### Is This Known?

**Assessment:** The Bures triangle inequality is well known. The data processing contractivity is well known. But the specific application to the Petz recovery map using functoriality to get an inequality specifically for tau = 1 - F_Petz does NOT appear in the literature. The key point is:

- Generic recovery maps do NOT compose functorially (Parzygnat-Buscemi 2023 proved this)
- The JRSWW universal recovery map fails compositionality
- Only the Petz map allows this clean decomposition

So while each ingredient is known, the *specific combination* -- that tau_Petz satisfies a Bures-metric triangle inequality for sequential channels, and that this FAILS for other recovery maps -- appears to be new.

**Wilde assessment:** He would likely say "this follows from Bures triangle + DPI + Petz functoriality" (which is true), but might find it interesting that it specifically fails for the JRSWW map.

### Connection to Entropy Production

From the bound tau <= 1 - exp(-Sigma/2):

    sqrt(tau_{12}) <= sqrt(1 - exp(-Sigma_1/2)) + sqrt(1 - exp(-Sigma_2/2))

For small Sigma: sqrt(tau) ~ sqrt(Sigma/2), so:

    sqrt(Sigma_{12}/2) <= sqrt(Sigma_1/2) + sqrt(Sigma_2/2)

which gives Sigma_{12} <= 2*(sqrt(Sigma_1) + sqrt(Sigma_2))^2.

But we also know Sigma_{12} = Sigma_1 + Sigma_2 (entropy production is additive for sequential processes). So this is weaker than the exact result. The tau inequality is weaker but works directly with the operationally measurable quantity (fidelity).

---

## Direction 2: Tensor Product Scaling

### Setup

Consider N^{otimes n} acting on rho^{otimes n} with reference sigma^{otimes n}.

### Key Facts

**Fact (Berta-Tomamichel 2015):** The fidelity of recovery is multiplicative:

    F(A;B|C)_{rho^{otimes n}} = [F(A;B|C)_rho]^n

This was proven in "The Fidelity of Recovery is Multiplicative" using SDP duality.

### Derivation

For n copies of channel N:

    F_n = F(rho^{otimes n}, R_{sigma^{otimes n}}(N^{otimes n}) o N^{otimes n}(rho^{otimes n}))

**Claim:** F_n = F_1^n (by multiplicativity of fidelity of recovery).

This requires a slight translation: the Berta-Tomamichel result is for tripartite states ρ_ABC and recovery of A from C. Our setting has ρ going through channel N. The connection: by the Stinespring dilation, N(rho)_B ~ Tr_E[V rho V^dagger]_{BE} with a tripartite purification rho_{ABE}. The Petz recovery from B is related to the fidelity of recovery F(A;R|B) where R is the reference.

Under this identification, for product states and product channels:

    F_n = F_1^n

Therefore:

    tau_n = 1 - F_1^n = 1 - (1 - tau_1)^n

### Theorem 2 (Tensor Product Scaling)

    tau(rho^{otimes n}, N^{otimes n}, sigma^{otimes n}) = 1 - (1 - tau(rho, N, sigma))^n

**Corollary 2.1** (Linear regime): For tau_1 << 1:

    tau_n ~ n * tau_1    (first order)

**Corollary 2.2** (Saturation): For large n:

    tau_n -> 1    exponentially fast when tau_1 > 0

**Corollary 2.3** (Concentration):

    tau_n = 1 - exp(n * ln(1 - tau_1)) = 1 - exp(-n * tau_eff)

where tau_eff = -ln(1 - tau_1) ~ tau_1 + tau_1^2/2 + ...

### Physical Interpretation: The Arrow of Time is Exponentially Amplified

This is a key insight:

**For a single channel with tiny irreversibility tau_1 = epsilon, after n applications:**

    tau_n = 1 - (1 - epsilon)^n ~ 1 - exp(-n*epsilon)

This means: even a TINY per-step irreversibility accumulates exponentially toward perfect irreversibility. After n ~ 1/epsilon steps, tau_n ~ 1 - 1/e ~ 0.63.

**This is a quantitative version of "the arrow of time emerges from microscopic irreversibility."**

In the Huang framework: microscopic tau ~ 0 (nearly reversible) but macroscopic tau ~ 1 (completely irreversible) after O(1/tau) steps. The transition from "time doesn't exist" to "time is absolute" is controlled by n*tau_1.

### Is This Known?

**Assessment:** The multiplicativity of the fidelity of recovery IS known (Berta-Tomamichel 2015). But the translation into the tau language and the physical interpretation as "exponential amplification of the time arrow" is new and non-trivial.

The formula tau_n = 1 - (1-tau_1)^n is elementary given multiplicativity, but:
1. The connection to the emergence of the macroscopic arrow of time is new
2. The characteristic scale n* = 1/tau_1 as a "decoherence time" in the tau framework is a new interpretation
3. The exponential saturation tau_n ~ 1 - exp(-n*tau_eff) mirrors the Crooks bound tau <= 1 - exp(-Sigma/2), suggesting tau_eff plays the role of a single-step entropy production

**Wilde assessment:** He knows the multiplicativity result (it's Berta's paper, Wilde was involved in related work). The formula itself would be immediate for him. The physical interpretation might be interesting but is more physics than math. Verdict: **known math, potentially new physics interpretation**.

---

## Direction 3: Differential Structure -- Speed Limit for tau

### Setup

Consider a Lindbladian evolution: d(rho)/dt = L(rho), generating the family of channels N_t = exp(tL).

Define tau(t) := tau(rho, N_t, sigma) = 1 - F(rho, R_{sigma,N_t} o N_t(rho)).

### Question: What is d(tau)/dt and is there a speed limit?

### Derivation

**Step 1.** tau(t) = 1 - F(t), so d(tau)/dt = -dF/dt.

**Step 2.** The time derivative of fidelity F(rho, rho_recovered(t)). Using the Bures metric:

    d_B(t)^2 = 2(1 - F(t)) = 2*tau(t)

So d_B(t) = sqrt(2*tau(t)), and:

    d(d_B)/dt = (1/sqrt(2*tau)) * d(tau)/dt

**Step 3.** From the quantum speed limit literature (Deffner-Lutz, Funo-Shiraishi-Watanabe):

For Lindbladian dynamics, the Bures distance satisfies:

    |d(d_B)/dt| <= v_QSL

where v_QSL depends on the Lindbladian generator. For the Mandelstam-Tamm type:

    v_QSL = Delta_H / hbar  (for unitary part)

plus dissipative corrections from the Lindblad operators.

**Step 4.** More precisely, for the full Lindbladian L = -i[H,*] + sum_k (L_k * L_k^dagger - {L_k^dagger L_k, *}/2):

    |d(tau)/dt| = sqrt(2*tau) * |d(d_B)/dt| <= sqrt(2*tau) * v_QSL(L)

### Theorem 3 (Speed Limit for tau)

    |d(tau)/dt| <= 2*sqrt(tau) * v_QSL(L)

where v_QSL(L) is the quantum speed limit velocity for the Lindbladian L.

For small tau (near-perfect recovery regime):

    |d(tau)/dt| <= 2*sqrt(tau) * v_QSL

This means: **the rate at which time's arrow can grow is bounded by both the current arrow strength AND the speed limit.**

**Corollary 3.1** (Initial growth): At t=0 where tau(0) = 0 (unitary at t=0):

    d(tau)/dt |_{t=0} = 0

This is because sqrt(tau(0)) = 0. The arrow of time starts growing at zero rate!

Wait -- this needs more care. At t=0, rho is unchanged and the channel is the identity, so tau=0 exactly. The derivative:

    d(tau)/dt = -dF/dt = -(d/dt)F(rho, R_{sigma,N_t} o N_t(rho))|_{t=0}

For N_0 = id, R_{sigma,id} = id, so:

    F(rho, N_0(rho)) = F(rho, rho) = 1

and:

    dF/dt|_{t=0} = d/dt F(rho, R_{sigma,N_t}(N_t(rho)))|_{t=0}

This requires computing the derivative of the Petz map itself. Let me approach this differently.

### Alternative: Using the entropy production rate

Since tau <= 1 - exp(-Sigma/2), and for small tau this gives tau <= Sigma/2:

    d(tau)/dt <= (1/2) * d(Sigma)/dt * exp(-Sigma/2)

The entropy production rate for Lindbladian dynamics is:

    sigma_dot := d(Sigma)/dt = -d/dt D(N_t(rho) || N_t(sigma))

where D is the relative entropy. For detailed-balance Lindbladians:

    sigma_dot = EP(rho(t))  (the entropy production rate)

### Theorem 3' (Entropy Production Speed Limit for tau)

    d(tau)/dt <= (1/2) * exp(-Sigma(t)/2) * sigma_dot(t)

where sigma_dot is the instantaneous entropy production rate.

**Corollary 3'.1:** At early times (Sigma ~ 0):

    d(tau)/dt |_{t=0+} <= sigma_dot(0)/2

So the initial growth rate of tau is bounded by half the initial entropy production rate.

### A New Speed-Limit Inequality

Combining with the Lindbladian speed limit on entropy production (from the literature):

    sigma_dot <= 2 * Delta_H^2 / hbar^2   (plus dissipative terms)

we get:

    d(tau)/dt <= Delta_H^2 / hbar^2 * exp(-Sigma/2)

### Is This Known?

**Assessment:** The individual pieces (quantum speed limit for fidelity, entropy production rate bounds) are known. But:

1. The specific inequality d(tau)/dt <= (1/2) exp(-Sigma/2) * sigma_dot is new as stated
2. The observation that d(tau)/dt = 0 at t=0 (the arrow of time starts growing at zero velocity) is a new physical insight
3. The connection between speed limits and the rate of emergence of time's arrow is new

However, Theorem 3' is essentially just differentiating the known bound tau <= 1 - exp(-Sigma/2). So it's really a corollary of Fawzi-Renner + Crooks, not a new result per se.

**Wilde assessment:** Would recognize this as a straightforward differentiation of known bounds. Not impressed mathematically, but might find the physical framing interesting. Verdict: **not genuinely new mathematics**.

### What WOULD Be New: A Reverse Speed Limit

The interesting question is: is there a LOWER bound on d(tau)/dt? That is:

**If the system is dissipative (L has Lindblad operators), must tau grow at some minimum rate?**

This would be related to the modified log-Sobolev inequality (MLSI):

    D(N_t(rho) || pi) <= exp(-alpha_MLSI * t) * D(rho || pi)

where pi is the fixed point. The MLSI constant alpha_MLSI gives a minimum rate of entropy production, which would translate to a minimum rate of tau growth.

### Theorem 3'' (Lower Speed Limit for tau, conditional)

If the Lindbladian L satisfies an MLSI with constant alpha, then for t >= 0:

    tau(t) >= 1 - exp(-(1/2) * D(rho||pi) * (1 - exp(-alpha*t)))

This gives a MINIMUM rate at which the arrow of time must grow, governed by the MLSI constant.

**Is this new?** The MLSI bound on relative entropy decay is known. Translating it to a lower bound on tau is immediate. But the PHYSICAL STATEMENT is new: "the arrow of time must emerge at a rate controlled by the log-Sobolev constant."

---

## Direction 4: tau and Channel Capacity

### Setup

For a quantum channel N, the quantum capacity is:

    Q(N) = lim_{n->inf} (1/n) max_rho I_c(rho, N^{otimes n})

where I_c is the coherent information.

### Relation to tau

**Key observation:** For a degradable channel N with complement N^c:

    I_c(rho, N) = S(N(rho)) - S(N^c(rho))

The Petz recovery fidelity tau is related to how well we can recover the input from the output. The coherent information measures how much quantum information survives.

**Claim:** For degradable channels, there should be a direct relation between tau and Q(N).

### Derivation Attempt

For degradable channels, N^c = D o N for some degrading channel D. Then:

    I_c(rho, N) = S(N(rho)) - S(D o N(rho))

The recovery fidelity tau measures something different: how well R_Petz o N recovers rho, not how much information survives in N(rho).

**Connection via the relative entropy:**

    D(rho || sigma) - D(N(rho) || N(sigma)) = Sigma (entropy production)

and I_c(rho, N) = S(B)_psi - S(BE)_psi = -S(A|B)_psi for Stinespring isometry |psi>_{ABE}.

### Theorem 4 (tau Bound on Coherent Information)

From Fawzi-Renner:

    I(A;B|C) >= -log F(A;B|C) = -log(1 - tau)

For the channel setting, identifying C = B (output), B = R (reference), A = A (input):

    I_c(rho, N) is NOT directly equal to I(A;R|B)

but rather I_c(rho, N) = I(A;B)_psi - I(A;E)_psi where psi is the Stinespring dilation.

Actually, the correct identification is:

    I(R;A|B)_psi >= -log F(R;A|B)

where R is the purifying system, A is the input, B is the output. And I(R;A|B) = I(R;A) - I(R;A;B) ... this is getting complicated.

### Cleaner Approach

The Fawzi-Renner bound states:

    D(N(rho) || N(sigma)) >= D(rho || sigma) - D(rho || R o N(rho))   ... not quite

Actually, the relevant form is:

    D(rho || sigma) - D(N(rho) || N(sigma)) >= -2 log F(rho, R_{sigma,N} o N(rho))
                                                = -2 log(1 - tau)

Rearranging:

    tau >= 1 - exp(-(Sigma)/2)

This we already know. The question is whether there's a direct tau-capacity relation.

### No Direct Clean Relation

After careful analysis: tau and channel capacity measure different things.
- tau measures how well a SPECIFIC input state can be recovered
- Q(N) measures the maximum rate of reliable quantum information transmission

There is no clean formula "Q(N) = f(tau)" because:
1. tau depends on the choice of rho and sigma, while Q involves optimization
2. Degradable channels can have the same tau but different Q
3. The relation goes through the coherent information, which involves entropy balances, not fidelity directly

**What IS true:**
- If Q(N) > 0, then there exists a code such that tau_code -> 0 (by definition of capacity)
- If tau(rho, N, sigma) = 0 for all rho, then N is reversible and Q(N) = log(dim A)
- For the maximally mixed input, tau is related to the entanglement fidelity, which IS related to capacity through the quantum Gilbert-Varshamov bound

### Is This Direction Productive?

**Verdict: No.** The relationship between tau and channel capacity is indirect and goes through well-known intermediaries (coherent information, entanglement fidelity). No clean new inequality emerges.

**Wilde assessment:** Would immediately point out that channel capacity involves regularization and optimization over inputs, while tau is a fixed-input quantity. Not the same beast. Verdict: **dead end**.

---

## Direction 5: Conditional tau and CMI

### Setup

For a tripartite state rho_ABC, define:

    tau(A|B) := 1 - max_{R: B->AB} F(rho_AB, R(rho_B))

This is precisely 1 - F(A;_|B) in the Fawzi-Renner notation (with the roles adjusted).

### Known Results

**Fawzi-Renner (2015):**

    I(A:C|B) >= -log F(A:C|B) = -log(1 - tau_{A:C|B})

**Improved bound (Wilde, 2015):**

    D(rho_ABC || sigma_ABC) - D(N(rho) || N(sigma)) >= -2 log F(rho, R o N(rho))

and various rotated versions.

### What Would Be New?

**Question:** Is there an UPPER bound on I(A:C|B) in terms of tau?

    I(A:C|B) <= f(tau_{A:C|B}) ?

Since I(A:C|B) can be arbitrarily large (for highly entangled states) while tau in [0,1], any such bound would need dimension dependence:

    I(A:C|B) <= g(tau, d_A, d_B, d_C)

**Derivation attempt:**

    I(A:C|B) = S(AB) + S(BC) - S(B) - S(ABC)

Using the Fannes-Audenaert inequality:

    |S(rho) - S(sigma)| <= T * log(d-1) + H(T)

where T = (1/2)||rho - sigma||_1 and H is binary entropy.

From the relation between fidelity and trace distance:

    1 - F <= T <= sqrt(1 - F^2)

So: T <= sqrt(2*tau - tau^2) <= sqrt(2*tau)

But this gives a bound on S(rho_AB) - S(rho_AB^recovered), not on I(A:C|B) directly.

### Theorem 5 (Conditional tau Quantifies CMI)

Combining Fawzi-Renner lower bound with the Fannes upper bound:

    -log(1 - tau) <= I(A:C|B) <= 4*sqrt(2*tau) * log(d_A) + 4*H(sqrt(2*tau))

where d_A = dim(A) and H is binary entropy.

**Assessment:** The lower bound is Fawzi-Renner. The upper bound follows from Fannes applied four times to the four entropy terms in CMI. This is a valid new inequality but it's essentially combining known bounds.

### Is This Known?

The individual bounds (Fawzi-Renner and Fannes-Audenaert) are well known. The combination giving a two-sided bound on CMI in terms of tau is, to my knowledge, not explicitly stated in the literature, but would be considered straightforward by experts.

**Wilde assessment:** Would say "this is just Fannes + Fawzi-Renner, nothing deep." Verdict: **not genuinely new**.

### More Interesting: Chain Rule for Conditional tau

For four systems A, B, C, D, the chain rule for CMI gives:

    I(A:CD|B) = I(A:C|B) + I(A:D|BC)

Does this translate to a chain rule for tau?

    tau_{A:CD|B} vs tau_{A:C|B} and tau_{A:D|BC}

From the lower bounds:

    -log(1 - tau_{A:CD|B}) <= -log(1-tau_{A:C|B}) + (-log(1-tau_{A:D|BC}))...

NO, this goes the wrong way. The chain rule is:

    I(A:CD|B) = I(A:C|B) + I(A:D|BC)

and:

    -log(1-tau_{A:CD|B}) <= I(A:CD|B) = I(A:C|B) + I(A:D|BC)

We want to relate this to tau_{A:C|B} and tau_{A:D|BC}. We have:

    I(A:C|B) >= -log(1 - tau_{A:C|B})
    I(A:D|BC) >= -log(1 - tau_{A:D|BC})

So:
    -log(1 - tau_{A:CD|B}) <= -log(1 - tau_{A:C|B}) + (-log(1 - tau_{A:D|BC}))

Wait, this says nothing useful because the inequality goes the wrong way on both sides.

However, from the MULTIPLICATIVITY of fidelity of recovery:

    F(A:CD|B) <= F(A:C|B) * F(A:D|BC)   ... is this true?

If the fidelity is SUB-multiplicative in this sense, then:

    1 - tau_{A:CD|B} <= (1 - tau_{A:C|B})(1 - tau_{A:D|BC})

which gives:

    tau_{A:CD|B} >= tau_{A:C|B} + tau_{A:D|BC} - tau_{A:C|B}*tau_{A:D|BC}

This would be a SUPER-ADDITIVITY result for tau under the chain rule.

**But is F sub-multiplicative in this sense?** This is NOT the same as the Berta-Tomamichel multiplicativity (which is about tensor products). This is about conditioning on more systems.

Actually, by data processing: tracing out D from a recovery of ACD gives a recovery of AC, but the converse is not true. So we get:

    F(A:CD|B) <= F(A:C|B)   (recovering more is harder)

but NOT necessarily F(A:CD|B) <= F(A:C|B) * F(A:D|BC).

**Verdict:** The chain rule for conditional tau does not yield clean new results.

---

## Direction 6: tau Under Channel Mixing

### Setup

Let N = sum_i p_i N_i be a convex mixture of channels. How does tau(N) relate to tau(N_i)?

### Derivation

**Step 1.** The Petz recovery map for a mixture: R_{sigma, sum p_i N_i} is NOT equal to sum p_i R_{sigma, N_i} in general. This is because R_sigma depends nonlinearly on N.

**Step 2.** The fidelity F(rho, R_N o N(rho)):

    N(rho) = sum_i p_i N_i(rho)

This is a mixture of states. The Petz map for the mixture is:

    R_{sigma,N}(X) = sigma^{1/2} N^dagger(N(sigma)^{-1/2} X N(sigma)^{-1/2}) sigma^{1/2}

where N^dagger is the adjoint of N = sum p_i N_i.

**Step 3.** Since N^dagger = sum p_i N_i^dagger and N(sigma) = sum p_i N_i(sigma), we see that R_{sigma,N} depends on the mixture in a complicated way.

**Step 4.** Can we bound F(rho, R_N o N(rho)) in terms of F(rho, R_{N_i} o N_i(rho))?

**Concavity argument:** F(rho, *) is concave in the second argument (known property of fidelity). However, R_N o N is not a convex combination of R_{N_i} o N_i. So joint concavity of fidelity does not directly help.

**Step 5.** Data processing approach: Each N_i can be obtained from N by a "filter" (post-selection on which N_i was applied). But post-selection is not a CPTP map, so data processing doesn't apply.

### Partial Result

**Lower bound via joint concavity of fidelity:**

For any fixed recovery map R (not necessarily the Petz map):

    F(rho, R o N(rho)) = F(rho, R(sum p_i N_i(rho))) = F(rho, sum p_i R(N_i(rho)))
                       >= (sum_i p_i sqrt(F(rho, R(N_i(rho)))))^2  ... NO, this is wrong

Actually, F(rho, sum p_i sigma_i) is related to sum p_i F(rho, sigma_i) by:

    F(rho, sum p_i sigma_i) >= sum_i p_i F(rho, sigma_i)   (by concavity of fidelity in second arg)

... wait, fidelity is CONCAVE in each argument separately? Let me check.

F(rho, sigma) = (Tr sqrt(sqrt(rho) sigma sqrt(rho)))^2 is actually CONCAVE in sigma when rho is fixed (Uhlmann 1976).

But here, both the recovery map R_N and the channel output N(rho) depend on the mixture. So:

    R_N o N(rho) is NOT a convex combination of R_{N_i} o N_i(rho)

### Theorem 6 (No Clean Convexity for tau)

tau(N) = tau(sum p_i N_i) is NEITHER convex nor concave in the mixing weights {p_i} in general.

**Proof sketch (by counterexample reasoning):**
- If all N_i are unitary (hence reversible), tau(N_i) = 0 for all i
- But N = sum p_i U_i(*) U_i^dagger can be a highly irreversible channel (e.g., depolarizing channel from random unitaries)
- So tau(N) >> sum p_i tau(N_i) = 0
- This shows tau is NOT convex (the mixture has higher tau than the average)

Wait -- this shows tau IS convex-LIKE (mixture has larger tau). Let me reconsider:
- tau(N_i) = 0 for unitaries
- tau(sum p_i N_i) > 0 for non-trivial mixture
- So tau(mixture) > sum p_i tau(N_i)

This means tau is NOT CONCAVE (since concavity would require tau(mixture) >= sum p_i tau(N_i), but we want the OPPOSITE direction).

Actually wait: tau = 0 for all components, and tau(mixture) > 0. This means:
- Concavity: tau(mixture) >= sum p_i tau(N_i) = 0. This IS satisfied.
- Convexity: tau(mixture) <= sum p_i tau(N_i) = 0. This is VIOLATED.

So mixing UNITARIES shows that tau is **NOT convex** (mixing can increase tau).

Is tau concave? Consider mixing a very noisy channel (tau ~ 1) with the identity (tau = 0):
- N = p * id + (1-p) * N_noisy
- For p close to 1, tau(N) ~ 0 (close to identity)
- sum p_i tau(N_i) = (1-p) * 1 ~ 0
- So tau(N) ~ (1-p)*tau(N_noisy) ~ sum p_i tau(N_i)
- For intermediate p, tau(N) could be less than (1-p)*tau(N_noisy)

Actually, adding identity to a channel improves recovery, so tau could decrease faster than linearly. This would violate concavity.

### Is This Known?

**Assessment:** The non-convexity of infidelity under channel mixing is probably known but I cannot find it stated explicitly for the Petz recovery map. The counterexample with random unitaries is folklore.

**Wilde assessment:** Would consider this obvious. Verdict: **known/obvious**.

---

## ============================================================
## SYNTHESIS: What Is Genuinely New?
## ============================================================

After systematic exploration of all six directions, here is my assessment:

### Tier 1: Genuinely New Results (would interest Wilde)

**NONE of the individual results are genuinely new mathematics.** Each direction either:
(a) Reduces to combining known results (Bures triangle + DPI, multiplicativity, speed limits)
(b) Has no clean answer (capacity, mixing)
(c) Is a reformulation of known bounds (conditional tau, differential structure)

### Tier 2: New Combinations/Perspectives That Have Value

1. **Theorem 1 (Composition)**: sqrt(tau_{12}) <= sqrt(tau_1) + sqrt(tau_2)
   - New as explicitly stated for the Petz map
   - Relies on functoriality (UNIQUE to Petz, fails for JRSWW)
   - The fact that this FAILS for universal recovery maps is potentially interesting

2. **Theorem 2 corollary**: tau_n = 1 - (1-tau_1)^n as "exponential emergence of arrow"
   - The formula is immediate from multiplicativity
   - The interpretation as arrow-of-time emergence is new
   - The characteristic time n* = 1/tau is the "decoherence time" in the tau framework

3. **Theorem 3'' (Lower speed limit)**: MLSI constant controls minimum rate of tau growth
   - Translates a known mathematical bound into a new physical statement
   - "The arrow of time MUST emerge at a rate controlled by the log-Sobolev constant"

### Tier 3: What WOULD Impress Wilde

To genuinely impress Wilde, one would need:

**Option A: A tight composition inequality**
Show that the Bures triangle inequality for tau is TIGHT for some class of channels. That is, find channels where sqrt(tau_{12}) = sqrt(tau_1) + sqrt(tau_2). This would require a new mathematical result about when equality holds.

**Option B: A non-obvious duality**
Show that tau(rho, N, sigma) = tau(sigma, R_N, rho) in some sense -- a duality between forward and backward processes. This is almost true by the functorial structure but needs careful formulation.

**Option C: An incompatibility result**
Show that no recovery map can simultaneously satisfy:
(i) The Bures composition inequality
(ii) A tighter-than-Fawzi-Renner bound
(iii) Some natural operational property
This would give an information-theoretic "no-go theorem."

**Option D: A new monotone**
Define M(N) = -log(1 - max_rho tau(rho, N, sigma)) and show it has properties that distinguish it from known channel monotones (like coherent information, private information, etc.).

---

## ============================================================
## MOST PROMISING DIRECTION: The Functorial Composition Failure
## ============================================================

The most promising direction, which I believe IS genuinely new, comes from combining Direction 1 with the Parzygnat-Buscemi uniqueness result:

### Theorem (Functorial Composition Gap)

**Statement:** The composition inequality sqrt(tau_{12}) <= sqrt(tau_1) + sqrt(tau_2) holds for the Petz recovery map but FAILS in general for the JRSWW universal recovery map.

**Significance:** The JRSWW map achieves a tighter bound on tau (it gives the strongest known recovery), but it sacrifices compositionality. The Petz map satisfies a clean composition law but gives weaker recovery.

**This is a genuine trade-off that has not been formulated in the literature:**

There is a Pareto frontier between:
- RECOVERY QUALITY (how small is tau for a single channel)
- COMPOSITIONAL CONSISTENCY (how well tau composes for sequential channels)

The Petz map sits at the extreme of perfect compositionality, while the JRSWW map optimizes single-channel recovery.

### Why This Matters for Physics

In the tau framework, the Petz map is the "natural" retrodiction because it's the unique functor. The composition inequality then says:

**The natural measure of temporal asymmetry satisfies a triangle inequality in the Bures metric.**

This means tau has a natural GEOMETRIC structure: the space of quantum processes, equipped with the Petz recovery infidelity, has a (pseudo-)metric structure where:
- The "distance from reversibility" of a composed process is bounded by the sum of individual distances
- This distance is measured in the SQUARE ROOT of infidelity (Bures distance), not infidelity itself

The fact that sqrt(tau), not tau, is the natural metric quantity is itself an insight: **the natural measure of time's arrow is sqrt(tau), the Bures distance from reversibility.**

---

## ADDENDUM: One More Attempt -- tau and the Petz-Renyi Divergences

### Setup

The Renyi generalization: for alpha in (0,1) u (1,2]:

    D_alpha(rho || sigma) := (1/(alpha-1)) log Tr[rho^alpha sigma^{1-alpha}]

The Petz-Renyi relative entropy satisfies a DPI with remainder:

    D_alpha(rho||sigma) - D_alpha(N(rho)||N(sigma)) >= ... (various bounds)

### Generalized tau

Define:

    tau_alpha(rho, N, sigma) := 1 - F_alpha(rho, R_Petz o N(rho))

where F_alpha is the alpha-fidelity (or alpha-affinity):

    F_alpha(rho, sigma) = Tr[rho^alpha sigma^{1-alpha}]

For alpha = 1/2: F_{1/2} = F (Uhlmann fidelity squared... no, F_{1/2} = Tr[rho^{1/2} sigma^{1/2}] is the Bhattacharyya coefficient).

### One-Parameter Family

tau_alpha gives a one-parameter family of "arrow of time" measures, each with its own composition law and entropy production bound.

**Key question:** Do different alpha values give different ordering of channels by "temporal asymmetry"?

If tau_{alpha_1}(N1) < tau_{alpha_1}(N2) but tau_{alpha_2}(N1) > tau_{alpha_2}(N2), this would show that the arrow of time is NOT uniquely defined -- it depends on which notion of distinguishability you use.

### Is This Direction Known?

The Petz-Renyi divergences and their data processing properties are well studied (Muller-Lennert et al., Wilde-Winter-Yang 2014). The idea of alpha-dependent recovery has been explored. But the specific question of whether different alpha give different orderings of the "arrow of time" is, to my knowledge, new.

**Wilde assessment:** He is an expert on Renyi divergences. This direction would need a specific new inequality to be interesting, not just a general observation. Verdict: **potentially interesting if a specific result emerges**.

---

## ============================================================
## DIRECTION 7: THE GENUINELY NEW RESULT
## The Composition-Recovery Trade-off Theorem
## ============================================================

After the systematic exploration above, the most promising new result emerges from the INTERSECTION of Directions 1 and 2, combined with a deeper analysis of the Parzygnat-Buscemi uniqueness theorem. Here is the full derivation.

### Setup and Motivation

For any recovery map family R, define:
- tau_R(rho, N, sigma) := 1 - F(rho, R_{sigma,N} o N(rho))
- tau_R^{JRSWW}: uses the Junge-Renner-Sutter-Wilde-Winter universal recovery map
- tau_R^{Petz}: uses the Petz recovery map

**Known:** tau^{JRSWW} <= tau^{Petz} (JRSWW gives tighter/better recovery for individual channels).

### The Key Observation

From Parzygnat-Buscemi (2023): the Petz map is the UNIQUE retrodiction functor. This means it is the ONLY recovery map satisfying:

    R_{sigma}(N2 o N1) = R_{sigma,N1} o R_{N1(sigma),N2}   ... (Functoriality)

The JRSWW map does NOT satisfy this. Nor does any other known recovery map.

### Theorem 7 (Composition-Recovery Trade-off)

**Part (a) -- Composition bound for Petz tau:**

For sequential channels N1, N2:

    sqrt(tau^Petz_{12}) <= sqrt(tau^Petz_1) + sqrt(tau^Petz_2)

This is Theorem 1, proven above via functoriality + Bures triangle + DPI.

**Part (b) -- No such bound for JRSWW tau:**

There exist channels N1, N2 and states rho, sigma such that:

    sqrt(tau^JRSWW_{12}) > sqrt(tau^JRSWW_1) + sqrt(tau^JRSWW_2)

**Proof sketch of Part (b):**

The proof of Part (a) uses three ingredients:
1. Bures triangle inequality (holds for ANY states)
2. DPI contractivity (holds for ANY CPTP map)
3. Functoriality: R(N2 o N1) = R(N1) o R(N2)

Ingredient 3 is what fails for JRSWW. Concretely: the JRSWW recovery for the composed channel N2 o N1 is:

    R^{JRSWW}_{sigma}(N2 o N1) = integral over t of [rotated Petz map for (N2 o N1)]

This is NOT equal to the composition of individual JRSWW maps. Therefore, when we insert the intermediate state R^{JRSWW}_{sigma,N1}(rho_1) in the Bures triangle inequality, the third vertex of the triangle is NOT equal to the full JRSWW recovery of the composed channel.

More precisely, the recovery for the composed channel involves AVERAGING over rotations t of the composed channel's modular operator, while composing individual recovery maps averages over SEPARATE modular operators. These are generically different.

**Construction of counterexample:**

Consider N1 = partial trace over E (environment coupling) and N2 = amplitude damping with high damping parameter. Choose rho to be entangled with the environment in a specific way such that:
- The JRSWW map for N2 o N1 "sees" the composed structure and optimizes globally
- But composing individual JRSWW maps uses intermediate modular operators that are suboptimal for the global recovery

In such cases, the JRSWW composed recovery can outperform the composed individual recoveries by enough to violate the Bures triangle bound.

[Note: A rigorous counterexample with explicit calculation is needed here. The argument above shows WHY it should fail but not the specific numerics.]

### Theorem 7 (cont.) -- The Trade-off

**Part (c) -- Impossibility:**

There is NO recovery map family R that simultaneously satisfies:
(i) Composition: sqrt(tau^R_{12}) <= sqrt(tau^R_1) + sqrt(tau^R_2) for all N1, N2
(ii) Optimality: tau^R(rho, N, sigma) <= tau^{Petz}(rho, N, sigma) for all rho, N, sigma
(iii) Strict improvement: tau^R < tau^{Petz} for at least one (rho, N, sigma) triple

In other words: ANY recovery map that improves on the Petz map for even one channel MUST violate the composition inequality for some pair of channels.

**Proof sketch:**

By Parzygnat-Buscemi, the Petz map is the UNIQUE functor. Any map satisfying (i) must be compositional. If it also satisfies (ii), then it achieves recovery at least as good as Petz for every single channel AND composes correctly. But the Petz map is the unique functor that is consistent with classical Bayesian inversion. So any map satisfying (i) must agree with the Petz map on all classical channels. By the uniqueness theorem's logic (which bootstraps from classical to quantum), it must agree with Petz everywhere. This contradicts (iii).

[Note: This proof sketch uses the Parzygnat-Buscemi uniqueness in a strong form. The precise conditions under which their uniqueness theorem applies need to be checked carefully.]

### Why This Is Interesting

**For mathematics (Wilde's perspective):**

This is a structural impossibility theorem. It says that the "metric structure" of temporal asymmetry (the triangle inequality for sqrt(tau)) is INCOMPATIBLE with having a tighter-than-Petz recovery bound. You must choose one or the other.

This is analogous to known impossibility results in quantum information, such as:
- No-cloning theorem: can't copy and maintain coherence
- No-broadcasting theorem: can't distribute quantum states
- Here: can't compose recoveries AND beat the Petz bound

**For physics (Huang's perspective):**

The Petz map IS retrodiction (uniquely). So tau^Petz is the "true" arrow of time. The fact that it satisfies a triangle inequality means the arrow of time has a natural metric structure. Optimized recovery maps (JRSWW) can achieve better single-step recovery but lose this metric structure. In physical terms:

"The natural arrow of time is geometrically consistent (triangle inequality) but not optimally sharp. Any sharper measure of irreversibility loses geometric consistency."

This is a fundamental statement about the nature of temporal asymmetry.

---

## ============================================================
## DIRECTION 8: Forward-Backward Duality of tau
## (Bonus: Connection to Kwon-Kim Complex Entropy Production)
## ============================================================

### Setup

The Petz map satisfies R_{sigma,N}(N(sigma)) = sigma (exact recovery of the reference). What about tau(sigma, R_N, N(sigma))? That is: if we go BACKWARD through the Petz map and then try to recover by going forward through N, what happens?

### Derivation

Define:
- tau_forward := tau(rho, N, sigma) = 1 - F(rho, R_Petz o N(rho))
- tau_backward := tau(N(rho), R_Petz, N(sigma)) -- the infidelity of recovering N(rho) after applying R_Petz, using N(sigma) as reference

For the backward process: we start with N(rho) in the output space, apply R_Petz to get (approximately) rho, and then need to recover N(rho). The "recovery of the recovery" is:

    R_{N(sigma), R_Petz} -- the Petz map for the Petz map itself

By functoriality applied to the identity channel id = N o R_Petz (approximately):

If N o R_Petz were exactly the identity, then tau_backward = 0. But it's not exact unless tau_forward = 0.

### Key Result

By the Petz map construction:

    R_{sigma,N}(X) = sigma^{1/2} N^*(N(sigma)^{-1/2} X N(sigma)^{-1/2}) sigma^{1/2}

where N* is the adjoint (Heisenberg picture) of N.

The fidelity F(N(rho), N o R_Petz o R_Petz o N(rho)):
- If we define the "round trip" map T := N o R_Petz (forward after backward)
- Then T(N(rho)) should be close to N(rho), and the infidelity measures how much is lost

### Theorem 8 (Forward-Backward Symmetry)

For the Petz recovery map with reference sigma:

    tau_forward(rho, N, sigma) = tau_backward(N(rho), R_Petz, N(sigma))

when sigma = rho (the Petz map with the true state as reference).

**Proof:** When sigma = rho, the Petz map R_{rho,N} satisfies:

    F(rho, R_{rho,N}(N(rho))) = F_forward

Now consider the backward direction:

    F(N(rho), N(R_{rho,N}(N(rho)))) ...

Hmm, this is NOT automatically equal. By data processing under N:

    F(N(rho), N(R_{rho,N}(N(rho)))) >= F(rho, R_{rho,N}(N(rho))) = F_forward

So tau_backward <= tau_forward: the backward process has LESS irreversibility than the forward process (or equal). This makes physical sense: applying a channel and then undoing it is a contraction, and the reverse process (undoing and then reapplying) contracts further.

### Refined Statement

    tau_backward(N(rho), R_Petz, N(sigma)) <= tau_forward(rho, N, sigma)

with equality iff N is a unitary (tau = 0).

This says: **the backward arrow of time is always WEAKER than the forward arrow.** Applying retrodiction and then prediction again loses more information than the initial prediction alone recovered.

**Is this known?** The DPI-based inequality F(N(A), N(B)) >= F(A,B) is trivially known. The interpretation in terms of forward/backward arrows of time is new but the math is trivial. Verdict: **trivial math, nice packaging**.

### Connection to Kwon-Kim Complex Entropy Production (PRX 2019)

Kwon and Kim (2019) showed that the entropy production has an IMAGINARY part when the channel involves coherence transfers:

    Sigma_complex = Sigma_real + i * Sigma_imag

where:
- Sigma_real = D(rho || R o N(rho)) >= 0 (the standard entropy production)
- Sigma_imag witnesses the "broken time-reversal symmetry" of the channel

In the tau framework:
- tau is sensitive to |Sigma_real| (through the bound tau <= 1 - exp(-Sigma_real/2))
- The imaginary part Sigma_imag has NO effect on tau (since F is a real quantity)

This means: **tau is blind to the imaginary part of entropy production.** The "arrow of time" measured by tau captures thermodynamic irreversibility but NOT coherence-based asymmetry.

This is actually an important limitation of the tau framework that deserves discussion in Paper 1.

---

## ============================================================
## FINAL VERDICT: RANKED LIST OF RESULTS
## ============================================================

### Rank 1: Composition-Recovery Trade-off (Theorem 7)
**Novelty: HIGH | Math difficulty: MEDIUM | Requires: Explicit counterexample**
- The impossibility theorem (Part c) would be the strongest result IF it can be made fully rigorous
- Needs the Parzygnat-Buscemi uniqueness theorem in a specific strengthened form
- Would be publishable as a standalone result

### Rank 2: Exponential Emergence of Arrow of Time (Theorem 2)
**Novelty: MEDIUM | Math difficulty: LOW | Interpretation: HIGH**
- tau_n = 1 - (1-tau_1)^n with characteristic scale n* = 1/tau_1
- Known math (Berta-Tomamichel multiplicativity), new physics interpretation
- Perfect for Paper 1 as a discussion point

### Rank 3: MLSI Lower Bound on tau Growth (Theorem 3'')
**Novelty: MEDIUM | Math difficulty: LOW | Physical impact: MEDIUM**
- "The arrow of time must emerge at a minimum rate controlled by the log-Sobolev constant"
- Translation of known MLSI results into tau language
- Good for a paragraph in Paper 1

### Rank 4: Bures Triangle Inequality for tau (Theorem 1)
**Novelty: LOW-MEDIUM | Math difficulty: LOW | Physical impact: MEDIUM**
- sqrt(tau_{12}) <= sqrt(tau_1) + sqrt(tau_2)
- Each ingredient is known; combination is natural but unstated
- Key insight: sqrt(tau) is the natural metric, not tau

### Rank 5: tau Blind to Imaginary Entropy Production
**Novelty: LOW | Math difficulty: TRIVIAL | Physical impact: MEDIUM-HIGH**
- Important limitation of the framework
- Connects to Kwon-Kim PRX 2019
- Should be acknowledged in Paper 1

### What Would Make Wilde Pay Attention

The ONLY result that would genuinely make Wilde think "this is new and non-trivial" is Theorem 7 Part (c), IF it can be made into a rigorous impossibility theorem with an explicit counterexample. Everything else, he would classify as "nice observation" at best and "obvious" at worst.

The honest assessment is: the tau = 1 - F framework is a useful REPACKAGING of known results, not a generator of new mathematics. Its value is in physical interpretation, not mathematical novelty.
