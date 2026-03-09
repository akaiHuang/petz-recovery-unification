# Layer 2 Analysis: Does a CPTP Map N_grav Exist with Delta-D = -ln(-g_00)?

**Date**: 2026-03-09
**Status**: Complete rigorous analysis
**Purpose**: Determine whether the JRSWW entropy production Sigma can be rigorously identified with -ln(-g_00) via an explicit CPTP map

---

## 0. The Central Problem

Paper 2 claims Sigma_grav = r_s/r and uses the JRSWW bound F >= exp(-Sigma/2). But the JRSWW bound is proven for CPTP maps. To make this connection rigorous, we need:

1. An explicit CPTP map N_grav
2. A reference state sigma
3. A proof that Delta-D(N_grav, rho, sigma) = D(rho || sigma) - D(N_grav(rho) || N_grav(sigma)) = -ln(-g_00)

Without this, Paper 2's "Sigma_grav" and the JRSWW "Sigma" are formally different quantities that happen to share a symbol.

**Notation**: Throughout, D(rho || sigma) = Tr[rho(ln rho - ln sigma)] is the Umegaki relative entropy. H_bin(x) = -x ln x - (1-x) ln(1-x) is the binary Shannon entropy (in nats). S(rho) = -Tr[rho ln rho] is the von Neumann entropy.

---

## 1. Approach 1: Pikovski Dephasing Channel

### 1.1 Channel Definition

The Pikovski mechanism (Nature Physics 11, 668, 2015) gives an explicit qubit dephasing channel. For a qubit encoding (keeping only the two energy eigenstates relevant to the superposition):

**Kraus operators**:
```
K_0 = sqrt((1+p)/2) * I
K_1 = sqrt((1-p)/2) * sigma_z
```

where p = D(t) is the decoherence factor, 0 <= p <= 1. For Pikovski:
```
p = exp(-Gamma * t^2)   [Gaussian decay, NOT exponential]
```
where Gamma depends on Delta_E and the gravitational potential difference.

**Action on a general qubit state** rho = (1/2)(I + r . sigma):
```
N_p(rho) = (1/2)(I + p*r_x*sigma_x + p*r_y*sigma_y + r_z*sigma_z)
```

The off-diagonal elements are multiplied by p; diagonal elements are preserved.

### 1.2 Computation of Delta-D with sigma = I/2

**Claim**: For sigma = I/2 (maximally mixed state), Delta-D = S(N(rho)) - S(rho).

**Proof**:
```
D(rho || I/2) = Tr[rho(ln rho - ln(I/2))]
              = Tr[rho ln rho] + Tr[rho] * ln 2
              = -S(rho) + ln 2
```

Similarly D(N(rho) || I/2) = -S(N(rho)) + ln 2.

Therefore:
```
Delta-D = D(rho || I/2) - D(N(rho) || I/2)
        = [-S(rho) + ln 2] - [-S(N(rho)) + ln 2]
        = S(N(rho)) - S(rho)
```

This is just the entropy increase under the channel. QED.

**Evaluation for pure input |+> = (|0> + |1>)/sqrt(2)**:

S(rho) = 0 (pure state).

After dephasing:
```
N_p(|+><+|) = (1/2)[(1+p)|+><+| + (1-p)|->< -|]
```

with eigenvalues (1+p)/2 and (1-p)/2. Therefore:

```
Delta-D = S(N(|+><+|)) = H_bin((1+p)/2)
        = -((1+p)/2) ln((1+p)/2) - ((1-p)/2) ln((1-p)/2)
```

**Key problem**: H_bin((1+p)/2) is bounded above by ln 2 approx 0.693 (achieved at p = 0). But the target -ln(-g_00) = r_s/r grows without bound as r -> 0.

**For small dephasing** (p close to 1, weak field): Let p = 1 - epsilon. Then:
```
H_bin((2-epsilon)/2) approx (epsilon^2)/2 + O(epsilon^3)    [expansion around p=1]
```

Meanwhile, we want Sigma = r_s/r, so for the Pikovski channel we would need p = exp(-E_G t / hbar) and thus:
```
H_bin((1+p)/2) approx (1-p)^2 / 2 approx (E_G t / hbar)^2 / 2
```

But this is (E_G t / hbar)^2 / 2, NOT E_G t / hbar.

**Verdict**: Delta-D for dephasing with sigma = I/2 is QUADRATIC in the dephasing parameter and BOUNDED by ln 2. It cannot equal -ln(-g_00) = r_s/r.

### 1.3 Can a Different Reference State sigma Fix This?

We need to find sigma such that:
```
D(rho || sigma) - D(N_p(rho) || N_p(sigma)) = -ln(p)
```

for the case p = exp(-r_s/r), so that Delta-D = r_s/r.

**General formula for qubit dephasing**: Let rho = |+><+| and sigma = diag(q, 1-q) in the {|0>, |1>} basis.

```
D(rho || sigma) = Tr[rho(ln rho - ln sigma)]
```

For rho = |+><+| = (1/2)(I + sigma_x):
```
rho = [[1/2, 1/2], [1/2, 1/2]]
ln sigma = diag(ln q, ln(1-q))
```

Therefore:
```
D(rho || sigma) = Tr[rho ln rho] - Tr[rho ln sigma]
                = 0 - (1/2)(ln q + ln(1-q))
                = -(1/2) ln(q(1-q))
                = ln 2 - (1/2) ln(4q(1-q))
```

Wait, let me be more careful. For a pure state rho = |psi><psi|, we have Tr[rho ln rho] = 0, so:
```
D(rho || sigma) = -Tr[rho ln sigma] = -<psi| ln sigma |psi>
```

For |+> and sigma = diag(q, 1-q):
```
D(|+><+| || sigma) = -(1/2)[ln q + ln(1-q)] = -(1/2) ln[q(1-q)]
```

After dephasing: N_p(|+><+|) = diag(1/2, 1/2) + p * [[0, 1/2],[1/2, 0]] has eigenvalues (1+p)/2 and (1-p)/2.

And N_p(sigma) = sigma (since sigma is diagonal, dephasing does not change it).

So:
```
D(N_p(rho) || N_p(sigma)) = D(diag((1+p)/2, (1-p)/2) + off-diag ... || diag(q, 1-q))
```

Actually, N_p(rho) in the {|0>,|1>} basis is:
```
N_p(|+><+|) = [[1/2, p/2], [p/2, 1/2]]
```

The eigenvalues of this matrix are (1+p)/2 and (1-p)/2 with eigenvectors |+> and |->. But D(A || B) when B is diagonal in a different basis is more complex.

Let me compute this directly. In the {|0>, |1>} basis:

```
N_p(rho) has matrix [[1/2, p/2], [p/2, 1/2]]
N_p(sigma) = sigma = diag(q, 1-q)
```

The relative entropy D(N_p(rho) || N_p(sigma)) = Tr[N_p(rho) * (ln N_p(rho) - ln N_p(sigma))].

Since N_p(sigma) is diagonal: ln N_p(sigma) = diag(ln q, ln(1-q)).

So:
```
Tr[N_p(rho) ln N_p(sigma)] = (1/2) ln q + (1/2) ln(1-q) = (1/2) ln[q(1-q)]
```

(The off-diagonal terms of N_p(rho) vanish when multiplied with the diagonal ln N_p(sigma).)

And:
```
Tr[N_p(rho) ln N_p(rho)] = -S(N_p(rho)) = -H_bin((1+p)/2)
```

Therefore:
```
D(N_p(rho) || N_p(sigma)) = -H_bin((1+p)/2) - (1/2) ln[q(1-q)]
```

Now compute Delta-D:
```
Delta-D = D(rho || sigma) - D(N_p(rho) || N_p(sigma))
        = [-(1/2) ln(q(1-q))] - [-H_bin((1+p)/2) - (1/2) ln(q(1-q))]
        = H_bin((1+p)/2)
```

**The reference state sigma (diagonal) CANCELS OUT completely!**

**This is a general result for qubit dephasing**: When sigma is diagonal in the dephasing basis, Delta-D is independent of sigma and equals the entropy increase S(N(rho)) - S(rho).

### 1.4 Non-Diagonal Reference States

What if sigma is NOT diagonal in the dephasing basis? Let sigma = |phi><phi| be a pure state with |phi> = cos(theta)|0> + sin(theta)|1>.

Then D(rho || sigma) could be infinite if rho is not in the support of sigma (which it won't be generically for pure sigma). So pure reference states are problematic.

Let sigma = (1/2)(I + s . vector-sigma) be a general qubit state with Bloch vector s. Then:

```
N_p(sigma) = (1/2)(I + p*s_x*sigma_x + p*s_y*sigma_y + s_z*sigma_z)
```

For D(rho || sigma) to be finite, we need supp(rho) subset supp(sigma), i.e., sigma must be full-rank, i.e., |s| < 1.

The computation becomes very involved for general Bloch vectors. However, we can note a key structural point:

**For any phase-covariant reference state** sigma (i.e., sigma commuting with the dephasing, meaning s_x = s_y = 0), we get the same result: Delta-D = H_bin((1+p)/2), independent of sigma.

**For general sigma**: The formula changes, but the fundamental limitation remains. The dephasing channel is a finite-dimensional (qubit) channel. The MAXIMUM entropy production over ALL input states rho and ALL reference states sigma is:

```
sup_{rho, sigma} Delta-D = sup_{rho, sigma} [D(rho || sigma) - D(N_p(rho) || N_p(sigma))]
```

By the data processing inequality, Delta-D >= 0 always. An upper bound comes from the contractivity coefficient of the channel.

**For the dephasing channel N_p, the contraction coefficient is**:
```
eta(N_p) = sup_{rho != sigma} D(N_p(rho) || N_p(sigma)) / D(rho || sigma)
```

For qubit dephasing channels, this equals p^2 (proven by Hiai, Ruskai, Petz). Therefore:
```
Delta-D = D(rho || sigma) * (1 - D(N(rho)||N(sigma))/D(rho||sigma))
        <= D(rho || sigma) * (1 - p^2)
```

This can be made arbitrarily large by choosing rho and sigma with large D(rho || sigma) (e.g., nearly orthogonal states). But the entropy production Sigma in the JRSWW sense requires a FIXED reference state sigma, and then considers the supremum over rho. The divergence D(rho || sigma) for a fixed full-rank sigma on a qubit is bounded:

```
D(rho || sigma) <= -ln(lambda_min(sigma))
```

where lambda_min is the smallest eigenvalue of sigma. This can be large but finite for any fixed sigma.

**Critical finding**: By choosing sigma with very small lambda_min, one can make Delta-D large. Specifically, if sigma = diag(epsilon, 1-epsilon) in an appropriate basis, then D(rho || sigma) ~ -ln(epsilon), and:

```
Delta-D ~ -ln(epsilon) * (1 - p^2)    [crude upper bound]
```

But this is (1-p^2) * ln(1/epsilon), not -ln(p). There is no choice of sigma that gives Delta-D = -ln(p) for ALL input states rho simultaneously.

### 1.5 Amplitude Damping as Alternative Qubit Channel

Consider an amplitude damping channel instead of dephasing:
```
K_0 = [[1, 0], [0, sqrt(p)]]
K_1 = [[0, sqrt(1-p)], [0, 0]]
```

This maps |1> -> sqrt(p)|1> + sqrt(1-p)|0>, modeling energy loss.

For sigma = |0><0| (ground state) and rho = |1><1| (excited state):

```
N(rho) = p|1><1| + (1-p)|0><0|
N(sigma) = |0><0|
```

```
D(rho || sigma) = infinity  [since <1|0><0|1> = 0]
```

So this particular choice diverges. Instead, take sigma = diag(1-epsilon, epsilon):

```
D(|1><1| || sigma) = -ln(epsilon)
D(N(|1><1|) || N(sigma)) = D(diag(1-p, p) || diag(1-epsilon, epsilon))
                          = (1-p) ln((1-p)/(1-epsilon)) + p ln(p/epsilon)
```

For small epsilon (nearly pure ground state reference):
```
D(N(rho) || N(sigma)) approx p * ln(p/epsilon) + (1-p) * ln(1)
                      = p * ln(p/epsilon)
                      = p * [ln p - ln epsilon]
                      = p * ln p + p * (-ln epsilon)
```

Therefore:
```
Delta-D = (-ln epsilon) - p * ln p - p * (-ln epsilon) - (1-p) * ln((1-p)/(1-epsilon))
        = (-ln epsilon)(1 - p) - p * ln p - (1-p) * ln((1-p)/(1-epsilon))
```

For epsilon -> 0:
```
Delta-D approx (1-p)(-ln epsilon) - p ln p - (1-p) ln(1-p)
             = (1-p)(-ln epsilon) + H_bin(p)     [where H_bin includes signs]
```

This diverges as ln(1/epsilon), multiplied by (1-p). NOT proportional to -ln(p).

**Special choice**: What if we set epsilon = p? Then sigma = diag(1-p, p).

```
D(|1><1| || diag(1-p, p)) = -ln(p) = r_s/r   [if p = exp(-r_s/r)]
```

And:
```
D(N(|1><1|) || N(diag(1-p, p)))
= D(diag(1-p, p) || diag(1-p, p))
= 0
```

So Delta-D = -ln(p) - 0 = -ln(p).

**THIS GIVES Delta-D = -ln(p) = r_s/r !**

But there is a critical problem: **the reference state sigma = diag(1-p, p) depends on p, which depends on r**. This means sigma is not a fixed reference state -- it changes with the gravitational field strength.

### 1.6 Assessment of Approach 1

**Summary of Results**:

| Channel | Reference sigma | Input rho | Delta-D | Equals -ln(p)? |
|---------|----------------|-----------|---------|-----------------|
| Dephasing | I/2 | \|+> | H_bin((1+p)/2) <= ln 2 | NO |
| Dephasing | diag(q, 1-q) | \|+> | H_bin((1+p)/2) <= ln 2 | NO |
| Dephasing | any full-rank | any | <= -ln(lam_min) * (1-p^2) | NO (wrong functional form) |
| Amp. damping | diag(1-p, p) | \|1> | -ln(p) | YES, but sigma depends on p |
| Amp. damping | fixed diag(q, 1-q) | \|1> | (1-p)ln((1-q)/q) + ... | NO (generically) |

**Fundamental obstruction**: For any FIXED reference state sigma on a finite-dimensional system, Delta-D for a qubit channel is bounded by 2 ln(dim) = 2 ln 2 approx 1.386. But -ln(-g_00) = r_s/r can be arbitrarily large.

*Proof*: Delta-D = D(rho || sigma) - D(N(rho) || N(sigma)) <= D(rho || sigma) <= S(rho || sigma) <= ln(dim/lam_min(sigma)). For fixed sigma, this is a finite constant.

**Actually, the bound is not so tight.** More precisely:
```
D(rho || sigma) <= ln(1/lambda_min(sigma))
```
for dim-dimensional systems. On a qubit this is indeed bounded for any fixed sigma. But the gravitational problem involves a CONTINUOUS variable (the field mode), not a qubit.

**Conclusion for Approach 1**: A qubit dephasing channel CANNOT reproduce Delta-D = -ln(-g_00) with a fixed reference state, because the qubit Hilbert space is too small. The amplitude damping channel can do it with a p-dependent reference state, but this is physically questionable (the reference state should not depend on the quantity being measured). The Pikovski channel, being a qubit channel, fundamentally cannot produce unbounded entropy production.

---

## 2. Approach 2: Bosonic Pure-Loss Channel

### 2.1 Channel Definition

The pure-loss bosonic channel E_eta with transmissivity eta in [0,1] acts on a single bosonic mode:

```
E_eta(rho) = Tr_E[U_BS (rho tensor |0><0|_E) U_BS^dagger]
```

where U_BS is the beam-splitter unitary mixing signal and environment modes with transmissivity eta:
```
a_out = sqrt(eta) * a_in + sqrt(1-eta) * a_vac
```

### 2.2 Entropy Production for Coherent State Inputs

**Input**: coherent state |alpha> with mean photon number N_s = |alpha|^2.
**Reference**: thermal state sigma_N_B with mean photon number N_B.

```
sigma_N_B = sum_{n=0}^{infinity} [N_B^n / (1+N_B)^{n+1}] |n><n|
```

**Output**: E_eta(|alpha><alpha|) = |sqrt(eta) alpha><sqrt(eta) alpha| (coherent state, pure).
**Output reference**: E_eta(sigma_N_B) = sigma_{eta * N_B} (thermal state with reduced photon number).

**Relative entropy computation**:

D(|alpha><alpha| || sigma_N_B) is the relative entropy of a coherent state with respect to a thermal state. Using the formula:

```
D(|alpha><alpha| || sigma_N_B) = N_s * ln(1 + 1/N_B) + ln(1 + N_B) - H(|alpha>)
```

Since |alpha> is pure, H(|alpha>) = S(|alpha><alpha|) = 0:

```
D(|alpha><alpha| || sigma_N_B) = |alpha|^2 / N_B - ln(1 + 1/N_B) * N_B + ln(1 + N_B)
```

Wait, let me derive this more carefully. For a coherent state displaced from the origin by alpha, with reference thermal state sigma with mean N_B:

The moment-generating function approach gives:
```
D(|alpha><alpha| || sigma_N_B) = |alpha|^2 / (N_B + 1) + g(N_B)
```
where g(N_B) = (N_B + 1) ln(N_B + 1) - N_B ln N_B.

Actually, the standard result (see Holevo, "Quantum Systems, Channels, Information", or Serafini "Quantum Continuous Variables"):

For a Gaussian state rho with covariance matrix V_rho and first moments d_rho, relative to a thermal state sigma with covariance matrix V_sigma = (2N_B + 1)I:

```
D(rho || sigma) = [Tr(V_sigma^{-1} V_rho) + (d_rho - d_sigma)^T V_sigma^{-1} (d_rho - d_sigma) - 2] / 2
                  + ln(det(V_sigma) / det(V_rho)) / 2
```

Hmm, this formula is for the general Gaussian case but the sign conventions vary. Let me use a direct computation instead.

**Direct computation for single-mode case**:

A coherent state |alpha> in the Fock basis is:
```
|alpha> = e^{-|alpha|^2/2} sum_n (alpha^n / sqrt(n!)) |n>
```

The thermal state is:
```
sigma_N = sum_n [N^n / (1+N)^{n+1}] |n><n|
```

Therefore:
```
<n| sigma_N |n> = N^n / (1+N)^{n+1}
ln <n|sigma_N|n> = n ln N - (n+1) ln(1+N)
```

```
D(|alpha><alpha| || sigma_N) = -S(|alpha><alpha|) - Tr[|alpha><alpha| ln sigma_N]
= 0 - sum_n |<n|alpha>|^2 * [n ln N - (n+1) ln(1+N)]
= -<n> * ln N + (<n> + 1) * ln(1+N)
= -|alpha|^2 ln N + (|alpha|^2 + 1) ln(1+N)
```

where I used <n> = |alpha|^2 and <n+1> = |alpha|^2 + 1 (expectation of n+1 in coherent state).

Wait, more carefully:
```
sum_n |<n|alpha>|^2 * n = |alpha|^2
sum_n |<n|alpha>|^2 * (n+1) = |alpha|^2 + 1
```

So:
```
D(|alpha><alpha| || sigma_N) = -|alpha|^2 ln(N) + (|alpha|^2 + 1) ln(1+N)
                              = |alpha|^2 ln((1+N)/N) + ln(1+N)
                              = |alpha|^2 ln(1 + 1/N) + ln(1+N)
```

Let me denote N_s = |alpha|^2. Then:
```
D_in := D(|alpha><alpha| || sigma_{N_B}) = N_s * ln(1 + 1/N_B) + ln(1 + N_B)   ... (*)
```

**After the loss channel**: The output is |sqrt(eta) alpha><sqrt(eta) alpha|, and the output reference is sigma_{eta N_B} (a thermal state with mean photon number eta * N_B, since the loss channel maps thermal states to thermal states with attenuated mean).

```
D_out := D(|sqrt(eta) alpha><sqrt(eta) alpha| || sigma_{eta N_B})
       = eta * N_s * ln(1 + 1/(eta N_B)) + ln(1 + eta N_B)
```

**Entropy production**:
```
Sigma = D_in - D_out
      = N_s * ln(1 + 1/N_B) - eta N_s * ln(1 + 1/(eta N_B))
        + ln(1 + N_B) - ln(1 + eta N_B)
```

### 2.3 Special Case: Vacuum Reference (N_B -> 0)

In the limit N_B -> 0, sigma_0 = |0><0|, but D(rho || |0><0|) is infinite unless rho is supported on |0>. So this limit is singular.

Take N_B small but nonzero:
```
ln(1 + 1/N_B) approx ln(1/N_B) = -ln(N_B)    [for N_B << 1]
ln(1 + 1/(eta N_B)) approx -ln(eta N_B) = -ln(eta) - ln(N_B)
```

Then:
```
Sigma approx N_s * (-ln N_B) - eta N_s * (-ln eta - ln N_B) + ...
     = N_s(-ln N_B) + eta N_s ln eta + eta N_s ln N_B + ...
     = N_s ln N_B (eta - 1) + eta N_s ln eta + ...
     = -(1-eta) N_s ln N_B + eta N_s ln eta + ...
```

As N_B -> 0, this diverges. So the vacuum reference does not give a finite answer.

### 2.4 Special Case: Large Thermal Reference (N_B >> 1, N_B >> N_s)

```
ln(1 + 1/N_B) approx 1/N_B
ln(1 + 1/(eta N_B)) approx 1/(eta N_B)
ln(1 + N_B) approx ln N_B
ln(1 + eta N_B) approx ln(eta N_B) = ln eta + ln N_B
```

Then:
```
Sigma approx N_s/N_B - eta N_s/(eta N_B) + ln N_B - ln eta - ln N_B
     = N_s/N_B - N_s/N_B - ln eta
     = -ln eta
```

**RESULT**: For a large thermal reference (N_B >> 1), the entropy production of the pure-loss channel with coherent state input is:

```
Sigma = -ln(eta) + O(1/N_B)
```

**This is the key result.** If eta = -g_00 = exp(-r_s/r), then Sigma = r_s/r.

### 2.5 Exact Computation for General N_B

Let us compute more carefully. Define x = 1/N_B. Then:

```
Sigma = N_s * ln(1 + x) - eta N_s * ln(1 + x/eta) + ln((1 + 1/x)/(1 + eta/x))
```

Hmm, let me just keep the original form and expand differently. Define f(eta, N_s, N_B):

```
f = N_s [ln(1 + N_B) - ln N_B] - eta N_s [ln(1 + eta N_B) - ln(eta N_B)]
  + ln(1 + N_B) - ln(1 + eta N_B)
```

Simplify:
```
f = N_s ln((1+N_B)/N_B) - eta N_s ln((1+eta N_B)/(eta N_B))
  + ln((1+N_B)/(1+eta N_B))
```

For N_B >> 1:
```
ln((1+N_B)/N_B) = ln(1 + 1/N_B) approx 1/N_B - 1/(2N_B^2) + ...
ln((1+eta N_B)/(eta N_B)) = ln(1 + 1/(eta N_B)) approx 1/(eta N_B) - 1/(2 eta^2 N_B^2) + ...
ln((1+N_B)/(1+eta N_B)) = ln(1 + (1-eta)N_B/(1+eta N_B)) approx ln(1/eta) - (1-eta)/(eta N_B) + ...
```

Wait, let me redo the last one:
```
(1+N_B)/(1+eta N_B) = (1 + N_B)/(1 + eta N_B)
```
For N_B >> 1: approx N_B/(eta N_B) = 1/eta.

So ln((1+N_B)/(1+eta N_B)) approx -ln(eta) + (1-eta)/(eta N_B) + O(1/N_B^2).

Hmm, more precisely: (1+N_B)/(1+eta N_B) = (1/eta)(1 + (1-eta)/(1+eta N_B)) * (eta N_B/(N_B)) * ...

Let me just do it directly:
```
(1+N_B)/(1+eta N_B) = [1/eta] * [eta + eta N_B] / [1 + eta N_B]
                     = [1/eta] * [1 + (eta - 1)/(1 + eta N_B)]
```

For large N_B:
```
approx (1/eta) * [1 + (eta-1)/(eta N_B)]
     = (1/eta) * [1 - (1-eta)/(eta N_B)]
```

So:
```
ln((1+N_B)/(1+eta N_B)) approx -ln(eta) - (1-eta)/(eta N_B) + O(1/N_B^2)
```

Now substituting:
```
Sigma approx N_s/N_B - eta N_s/(eta N_B) + [-ln eta - (1-eta)/(eta N_B)] + O(1/N_B^2)
     = N_s/N_B - N_s/N_B - ln eta - (1-eta)/(eta N_B) + ...
     = -ln(eta) - (1-eta)/(eta N_B) + ...
```

So to leading order in 1/N_B:
```
Sigma = -ln(eta) - (1-eta)/(eta N_B) + O(1/N_B^2)
```

The correction is O(1/N_B) and **negative**, meaning the actual Sigma is slightly less than -ln(eta) for finite N_B.

### 2.6 The Physical Identification eta = -g_00

**Question**: Why should the transmissivity equal -g_00?

There are two natural conventions:

**(A) Amplitude transmissivity**: A photon's frequency redshifts by the factor sqrt(-g_00). Since energy E = hbar omega, the energy transmissivity is sqrt(-g_00). As a beam-splitter model, this gives:
```
eta_A = sqrt(-g_00) = exp(-r_s/(2r))
```
Then Sigma = -ln(eta_A) = r_s/(2r). **Factor of 2 too small.**

**(B) Power/intensity transmissivity**: The power received at infinity from a source at r is reduced by a factor of (-g_00):
- Factor of sqrt(-g_00) from frequency redshift (energy per photon)
- Factor of sqrt(-g_00) from time dilation (photon arrival rate)
- Combined: power ratio = -g_00

```
eta_B = -g_00 = exp(-r_s/r)
```
Then Sigma = -ln(eta_B) = r_s/r. **Correct.**

**(C) Quantum state fidelity convention**: For a bosonic mode, the fidelity of a coherent state after the loss channel is:
```
F(|alpha>, |sqrt(eta) alpha>) = |<alpha|sqrt(eta) alpha>|^2 = exp(-|alpha|^2(1-sqrt(eta))^2)
```
This depends on |alpha|^2 and is NOT simply eta.

**Assessment of conventions**: Convention (B) gives the right answer, but the justification requires careful thought:

1. The gravitational redshift acts on EACH photon's frequency: omega_out = sqrt(-g_00) * omega_in
2. The number of photons is also reduced per unit coordinate time: dN/dt_out = sqrt(-g_00) * dN/dt_in (time dilation)
3. For a single-mode bosonic channel, the relevant quantity is the OCCUPATION NUMBER (mean photon number) of the mode: n_out = (-g_00) * n_in (both effects combined)
4. A pure-loss channel with eta = -g_00 correctly captures this

**HOWEVER**, there is a subtlety: de Paolis et al. (arXiv:2502.20521, 2025) have shown that the beam-splitter model for gravitational redshift has limitations. The gravitational redshift is NOT literally a beam-splitter interaction -- it is a unitary transformation of the mode structure. The beam-splitter model is an effective description valid for specific observational scenarios (e.g., photodetection at infinity).

### 2.7 Exact Formula for the Bosonic Channel

Collecting the results for the pure-loss bosonic channel:

**Channel**: E_eta (pure loss, transmissivity eta)
**Input**: coherent state |alpha> with N_s = |alpha|^2
**Reference**: thermal state sigma_{N_B} with mean photon number N_B

```
Sigma(eta, N_s, N_B) = N_s ln((1+N_B)/N_B) - eta N_s ln((1+eta N_B)/(eta N_B))
                     + ln(1+N_B) - ln(1+eta N_B)
```

**Limits**:
- N_B >> 1: Sigma -> -ln(eta)        [INDEPENDENT of input state N_s!]
- N_B << 1: Sigma -> divergent
- N_s = 0 (vacuum input): Sigma = ln(1+N_B) - ln(1+eta N_B) = ln((1+N_B)/(1+eta N_B))
  For N_B >> 1: Sigma -> -ln(eta). For N_B << 1: Sigma -> (1-eta) N_B.

**Key observation**: In the high-temperature limit (N_B >> 1), the entropy production becomes -ln(eta) regardless of the input state. This universality is physically appealing: it means the gravitational entropy production is a property of the CHANNEL, not the probe.

### 2.8 The Thermal Reference as Tolman-Gibbs State

The choice N_B >> 1 corresponds physically to a high-temperature reference. In the gravitational context, this can be motivated by the Unruh/Tolman effect: a static observer at radius r perceives the vacuum as a thermal state with temperature T(r). Taking the reference state to be the local thermal state is the standard choice in the JRSWW framework when applied to thermodynamic entropy production.

Specifically, if we identify the reference state with the Hartle-Hawking vacuum restricted to the region r < R (which appears thermal at the local Tolman temperature), then N_B grows as we approach the horizon, and the high-N_B limit becomes exact near strong-field regions.

**Assessment**: The identification Sigma = -ln(eta) = -ln(-g_00) = r_s/r works for the bosonic loss channel with thermal reference in the high-temperature limit. The remaining questions are:

1. Why should the reference be high-temperature? (Answer: Tolman effect gives local temperature T(r) = T_H/sqrt(1-r_s/r), which diverges near the horizon)
2. Why should eta = -g_00 rather than sqrt(-g_00)? (Answer: power transmissivity, accounting for both frequency shift and time dilation)
3. Is the bosonic loss channel the correct model for gravitational redshift? (Answer: it is an effective model, not a first-principles derivation)

---

## 3. Approach 3: Algebraic (Witten/Chandrasekaran)

### 3.1 Setup

In the crossed-product algebra framework (Witten, arXiv:2112.12828; Chandrasekaran et al., arXiv:2206.10780), gravitational observables in a subregion form a Type II_infinity von Neumann algebra.

For a static Schwarzschild black hole, the algebra of observables outside the horizon is:
```
A = A_QFT x_sigma R
```
where x_sigma denotes the crossed product with the modular automorphism group, and R accounts for the observer's clock.

The key feature: in this algebra, the relative entropy D(omega_1 || omega_2) is well-defined and finite for normal states omega_1, omega_2.

### 3.2 Modular Operator and Relative Entropy

For two states omega and phi on the algebra, the Araki relative entropy is:
```
S(omega || phi) = -<Psi_omega| ln Delta_{phi/omega} |Psi_omega>
```

where Delta_{phi/omega} is the relative modular operator.

For the Rindler wedge (the near-horizon approximation of Schwarzschild), the modular flow of the vacuum is the boost. The modular Hamiltonian is K = (2pi/kappa) * H_boost.

### 3.3 QRE for States at Different Radii

Consider two states:
- omega_r: a state describing excitations at radius r
- omega_inf: the corresponding state as seen from infinity

The gravitational channel N_grav maps omega_r to omega_inf. In the algebraic setting, this is the restriction from the algebra at r to the algebra at infinity (more precisely, the inclusion map of the subalgebra accessible from infinity into the full algebra at r).

**The Dorau-Much computation** (Section 1 of Route 2 notes):

Starting from the entanglement first law on a bifurcate Killing horizon:
```
D(omega_r || omega_0) = (2pi/kappa) * <delta T_{ab} xi^a xi^b>
```

With the Tolman relation and energy redshift:
```
D_r / D_inf = 1 / (1 - r_s/r)
```

leading to:
```
Sigma = D_r - D_inf = D_inf * r_s/r / (1 - r_s/r)
```

and the fractional loss:
```
(D_r - D_inf) / D_r = r_s/r     [EXACT for Schwarzschild]
```

### 3.4 Does This Give a CPTP Map?

**Critical issue**: The Dorau-Much framework does NOT explicitly construct a CPTP map. What it gives is:

1. A well-defined relative entropy D(omega_r || omega_0) at each radius
2. A monotonic decrease as we go from r to infinity
3. The fractional decrease equals r_s/r

The decrease of relative entropy is COMPATIBLE with the existence of a CPTP map (by the data processing inequality), but it does not CONSTRUCT one.

**What would be needed**: A CPTP map N: A(r) -> A(inf) such that:
```
D(omega_1 || omega_2) - D(N(omega_1) || N(omega_2)) = -ln(-g_00(r))
```

for appropriate states omega_1, omega_2.

**The restriction map** from A(r) to A(inf) IS a completely positive trace-preserving map (it's the partial trace over degrees of freedom between r and infinity). In the algebraic setting, it corresponds to the conditional expectation.

But the entropy production of this restriction map is:
```
Sigma = D_r - D_inf = D_inf * r_s/r / (1 - r_s/r)     [for Schwarzschild]
```

This is NOT exactly -ln(1 - r_s/r). The absolute Sigma depends on D_inf, which depends on the test state.

**The fractional loss** (D_r - D_inf)/D_r = r_s/r is state-independent, but it is a RATIO, not an absolute entropy production.

### 3.5 The Normalization Problem

The JRSWW bound involves an ABSOLUTE entropy production:
```
F >= exp(-Sigma/2)
```

But the algebraic computation gives a RELATIVE (fractional) entropy production. To convert:

If (D_r - D_inf)/D_r = r_s/r, then D_r - D_inf = D_r * r_s/r.

For D_r = 1 (unit QRE at the source), Sigma = r_s/r.

**But D_r is not generically equal to 1.** It depends on the test state. The identification Sigma_grav = r_s/r thus requires either:
(a) A normalization convention (divide by D_r), or
(b) Choosing test states with D_r = 1

Neither is automatic. This is the deepest unresolved issue.

### 3.6 Connection to Finite Radius

For an observer at finite radius r (not at the horizon), the modular flow is NOT a simple boost. The Bisognano-Wichmann theorem applies only at bifurcate Killing horizons. At finite r, the modular Hamiltonian receives corrections:

```
K_r = K_horizon + corrections
```

These corrections are known in the literature (see Casini, Huerta, Myers, "Towards a derivation of holographic entanglement entropy") but make exact computation difficult.

**Assessment**: The algebraic approach provides the most rigorous framework but falls short of constructing an explicit CPTP map with Sigma = -ln(-g_00). It gives a FRACTIONAL decrease of r_s/r for Schwarzschild, which becomes the absolute decrease only under a normalization assumption.

---

## 4. The Factor-of-2 Problem: Systematic Analysis

### 4.1 Catalog of Results

| Approach | Formula | Value for Schwarzschild | Factor |
|----------|---------|------------------------|--------|
| Frequency redshift | omega_out/omega_in = sqrt(-g_00) | exp(-r_s/(2r)) | 1/2 |
| Energy per photon | E_out/E_in = sqrt(-g_00) | exp(-r_s/(2r)) | 1/2 |
| Amplitude transmissivity | eta_A = sqrt(-g_00) | exp(-r_s/(2r)) | 1/2 |
| Power transmissivity | eta_P = -g_00 | exp(-r_s/r) | 1 |
| Tolman temperature ratio | T_inf/T_r = sqrt(-g_00) | exp(-r_s/(2r)) | 1/2 |
| Gravitational Landauer | Sigma = -ln(-g_00) | r_s/r | 1 |
| Bosonic loss (eta = -g_00) | Sigma = -ln(eta) | r_s/r | 1 |
| Bosonic loss (eta = sqrt(-g_00)) | Sigma = -(1/2)ln(-g_00) | r_s/(2r) | 1/2 |
| Modular flow (fractional) | (D_r - D_inf)/D_r | r_s/r | 1 |
| c_eff/c = sqrt(-g_00) | -ln(c_eff/c) | r_s/(2r) | 1/2 |
| g_00 = -exp(-r_s/r) [Paper 2] | -ln(-g_00) | r_s/r | 1 |

### 4.2 Where Each Factor of 2 Comes From

**The fundamental issue**: There are two natural "gravitational transmissivities":

1. **Single-particle energy ratio**: sqrt(-g_00) = exp(-r_s/(2r)). This is what a single photon experiences. It gives Sigma = r_s/(2r).

2. **Phase space / power ratio**: -g_00 = exp(-r_s/r). This accounts for BOTH the energy per photon AND the time dilation affecting the rate. It gives Sigma = r_s/r.

**Physical argument for using -g_00**: Consider a quantum channel that maps a field mode at radius r to a field mode at infinity. The channel must account for:
- The change in mode frequency: omega_inf = sqrt(-g_00) * omega_r
- The change in mode normalization: modes are normalized per unit proper time, so the normalization factor is also sqrt(-g_00)
- Combined effect on the occupation number: n_inf = -g_00 * n_r

This combined effect gives eta = -g_00, hence Sigma = -ln(-g_00) = r_s/r.

**Alternatively**: The Tolman-Oppenheimer-Volkoff equilibrium requires the local temperature T(r) = T_inf / sqrt(-g_00). The Landauer erasure cost at r is k_B T(r) ln 2. The RATIO of erasure costs is:
```
W(r)/W(inf) = T(r)/T(inf) = 1/sqrt(-g_00)
```

The dimensionless entropy production per bit is:
```
Sigma_bit = ln(W(r)/W(inf)) = -(1/2) ln(-g_00) = r_s/(2r)
```

This gives the HALF value! The full value r_s/r requires counting BOTH the erasure cost AND the transport cost (moving the information from r to infinity), which doubles the result.

### 4.3 Resolution

The factor-of-2 problem reflects a genuine ambiguity in what "gravitational entropy production" means:

**(i) Per-particle Sigma = r_s/(2r)**:
- Natural if measuring the degradation of a single quantum state
- Corresponds to sqrt(-g_00) as the fidelity
- Gives the CORRECT weak-field metric g_00 = -(1 - r_s/r + ...) when sqrt(-g_00) = 1 - r_s/(2r)

**(ii) Per-mode Sigma = r_s/r**:
- Natural if measuring the informational capacity of the gravitational channel
- Corresponds to -g_00 as the "power transmissivity"
- Gives g_00 = -exp(-r_s/r), the exponential metric

**IMPORTANT**: In the JRSWW framework, Sigma is defined as D(rho||sigma) - D(N(rho)||N(sigma)). This is a property of the PAIR (rho, sigma) and the channel N. Different physical quantities (energy per photon, power, occupation number) correspond to different choices of what the channel N represents.

**Paper 2's choice**: Paper 2 implicitly uses convention (ii) by defining Sigma_grav = r_s/r = 2|Phi|/c^2 and then setting sqrt(-g_00) = exp(-Sigma/2). This is self-consistent, but the factor of 2 in the exponent is a CONVENTION that corresponds to using the power transmissivity.

**The honest statement**: There exist CPTP channels for which Sigma = r_s/(2r) AND channels for which Sigma = r_s/r. The question of which is "the gravitational channel" depends on what physical process is being modeled. Paper 2's choice of Sigma = r_s/r corresponds to modeling the FULL phase-space degradation (both energy and timing), not just the energy shift.

---

## 5. Synthesis: The Key Mathematical Question

### 5.1 Precise Statement of What We Need

For which quantum channel N and reference state sigma does:
```
D(rho || sigma) - D(N(rho) || N(sigma)) = -ln(-g_00(r))
```

### 5.2 Explicit Constructions That Work

**Construction A: Bosonic loss channel with high-temperature thermal reference**

- Channel: E_eta with eta = -g_00 = exp(-r_s/r)
- Reference: sigma = sigma_{N_B} with N_B >> 1
- Input: any coherent state |alpha>
- Result: Sigma = -ln(eta) + O(1/N_B) = r_s/r + O(1/N_B)
- Status: WORKS in the large N_B limit

**Assumptions required**:
1. Gravitational redshift IS modeled by a pure-loss bosonic channel (effective model, not derived)
2. The transmissivity is eta = -g_00 (power convention, not amplitude convention)
3. The reference state is a high-temperature thermal state (motivated by Tolman/Unruh)
4. The high-temperature limit N_B >> 1 is a good approximation (holds near strong-field regions)

**Construction B: Amplitude damping with state-dependent reference** (Section 1.5)

- Channel: amplitude damping with parameter p = -g_00
- Reference: sigma = diag(1-p, p) [depends on p!]
- Input: |1><1|
- Result: Sigma = -ln(p) = r_s/r
- Status: WORKS but reference state depends on the field strength (physically questionable)

**Construction C: Algebraic (modular channel)**

- Channel: restriction from A(r) to A(infinity) (conditional expectation)
- Reference: vacuum/KMS state
- Input: coherent perturbation with unit QRE
- Result: Sigma = r_s/r (for appropriately normalized states)
- Status: WORKS formally but requires normalization assumption and is not an explicit CPTP map on a Hilbert space

### 5.3 What Does NOT Work

**Qubit dephasing (Pikovski)**: Cannot produce unbounded Sigma because the Hilbert space is too small. Delta-D <= ln 2 for any fixed reference state on a qubit.

**Amplitude damping with fixed reference**: Does not give -ln(p) for generic input states; the result depends on both input and reference.

**Bosonic loss with vacuum reference**: Sigma diverges (the vacuum is not a full-rank state in the Fock space sense needed for JRSWW).

**Any finite-dimensional channel with fixed reference**: Sigma is bounded by ln(dim), which is finite. The gravitational Sigma = r_s/r is unbounded.

---

## 6. The Factor-of-2 Problem: Can We Have a Single Consistent Definition?

### 6.1 The Two Natural Definitions

**Definition 1** (Paper 2 convention):
```
Sigma_grav^(1) = 2|Phi|/c^2 = r_s/r = -ln(-g_00^{exp})
```
where g_00^{exp} = -exp(-r_s/r) is the exponential metric.

**Definition 2** (single-particle redshift):
```
Sigma_grav^(2) = |Phi|/c^2 = r_s/(2r) = -ln(sqrt(-g_00^{exp}))
```

### 6.2 Relationship to the JRSWW Bound

The JRSWW bound is:
```
F(rho, R o N(rho)) >= exp(-Sigma/2)
```

Under Definition 1 with Sigma = r_s/r:
```
F >= exp(-r_s/(2r)) = sqrt(-g_00)
```
This identifies sqrt(-g_00) as the FIDELITY bound.

Under Definition 2 with Sigma = r_s/(2r):
```
F >= exp(-r_s/(4r)) = (-g_00)^{1/4}
```
This is a weaker bound.

Paper 2's identification sqrt(-g_00) = F_bound is precisely the statement that Sigma = r_s/r and the bound is saturated.

### 6.3 Is There a Unique Correct Choice?

**No.** The "correct" Sigma depends on the physical observable being tracked:

- If tracking **quantum state fidelity** of a single photon: the natural channel has Sigma = r_s/(2r), and F >= (-g_00)^{1/4}.
- If tracking **total information throughput** of a communication channel: the natural channel has Sigma = r_s/r, and F >= sqrt(-g_00).

Paper 2 chooses the latter, which gives the cleaner identification sqrt(-g_00) = F_bound. This is a CHOICE, not a theorem.

---

## 7. Conclusions

### 7.1 Does a CPTP Map Exist with Delta-D = -ln(-g_00)?

**Answer: CONDITIONAL YES.**

A bosonic pure-loss channel with transmissivity eta = -g_00, using a high-temperature thermal reference state, gives Delta-D = -ln(-g_00) in the limit N_B >> 1. This is an explicit, well-defined CPTP map.

**The conditions are**:
1. Transmissivity convention: eta must be the power transmissivity (-g_00), not the amplitude transmissivity (sqrt(-g_00)). This is physically motivated by accounting for both frequency shift and time dilation, but is a modeling choice.
2. Reference state: must be a high-temperature thermal state. This is motivated by the Tolman/Unruh effect but requires N_B >> 1.
3. The bosonic loss channel is an EFFECTIVE model of gravitational redshift, not a first-principles derivation.

### 7.2 What Channel and What Reference State?

**Best candidate**: Pure-loss bosonic channel E_eta with:
- eta = -g_00(r) = exp(-r_s/r) [for exponential metric] or (1-r_s/r) [for Schwarzschild]
- sigma = thermal state with mean photon number N_B >> 1 (Tolman/Unruh reference)
- Result: Sigma = -ln(eta) = r_s/r [exponential] or -ln(1-r_s/r) [Schwarzschild]

**Runner-up**: The algebraic (modular channel) construction via the crossed-product algebra. This is more rigorous but less explicit.

### 7.3 What Is the Closest We Can Get Without the Full Identification?

Even without specifying the channel:
- The DATA PROCESSING INEQUALITY guarantees that if ANY channel describes the propagation from r to infinity, then D(rho||sigma) decreases monotonically.
- The Dorau-Much computation shows the fractional decrease is r_s/r for Schwarzschild, independent of the test state.
- This is consistent with Sigma = -ln(-g_00) but does not uniquely determine it (because the fractional vs. absolute distinction remains).

### 7.4 What Additional Assumptions Are Needed?

To go from the existing results to a complete proof that Sigma_grav = -ln(-g_00), one needs:

1. **Channel construction**: An explicit CPTP map on a well-defined Hilbert space that represents gravitational propagation from r to infinity. The bosonic loss channel is a candidate but is not derived from first principles. The crossed-product construction (Witten 2022) provides the algebraic framework but has not been reduced to an explicit CPTP map.

2. **Reference state selection**: The choice of sigma must be physically motivated and r-independent. The Tolman thermal state is a natural candidate, but its temperature depends on r (through the Tolman relation), raising circularity concerns.

3. **Transmissivity convention**: The identification eta = -g_00 (power) rather than eta = sqrt(-g_00) (amplitude) must be justified from the physical process being modeled.

4. **High-temperature limit**: The result Sigma = -ln(eta) is exact only for N_B -> infinity. At finite N_B, there are O(1/N_B) corrections.

### 7.5 Classification: Theorem, Conjecture, or Wishful Thinking?

**Layer 1 (No-horizon argument)**: This is a THEOREM, conditional on one assumption: that gravitational redshift can be modeled as a quantum channel with finite entropy production. The conclusion (F > 0, hence tau < 1) follows rigorously from the JRSWW inequality.

**Layer 2 (Exponential metric from Sigma = r_s/r)**: This is a MOTIVATED CONJECTURE. The motivation is strong:
- Three independent routes converge on Sigma = -ln(-g_00)
- The bosonic loss channel provides an explicit CPTP map that reproduces this
- The fractional QRE decrease in the algebraic framework is exactly r_s/r

But the conjecture is not a theorem because:
- No first-principles derivation of N_grav as a CPTP map exists
- The transmissivity convention (eta = -g_00 vs sqrt(-g_00)) is a modeling choice
- The high-temperature limit for the reference state requires justification
- The saturation hypothesis (F = exp(-Sigma/2) as equality) is an additional assumption with no proof

**Layer 3 (Saturation)**: The claim that the Petz bound is SATURATED (F = exp(-Sigma/2) exactly) is the WEAKEST part. No known quantum channel saturates the JRSWW bound for Sigma > 0 (Petz proved saturation only for Sigma = 0). This is better described as SPECULATIVE: physically appealing but mathematically unsupported.

### 7.6 Honest Summary for Paper 2

The connection between gravity and JRSWW rests on three layers of decreasing rigor:

```
Layer 1 (THEOREM):     Finite Sigma + JRSWW  =>  F > 0  =>  tau < 1  =>  no horizon
Layer 2 (CONJECTURE):  Sigma = -ln(-g_00)    =>  F >= sqrt(-g_00)    =>  specific bound
Layer 3 (SPECULATION):  Saturation F = F_bound  =>  g_00 = -exp(-r_s/r)  =>  exponential metric
```

Each layer adds assumptions:
- Layer 1: gravitational redshift IS a quantum channel (widely accepted)
- Layer 2: the channel's entropy production equals -ln(-g_00) (supported by 3+ independent arguments, one explicit CPTP construction, but not derived from first principles)
- Layer 3: the Petz bound is saturated (no mathematical proof; physically motivated but not demonstrated for any Sigma > 0 channel)

---

## 8. Technical Appendix: Detailed Calculations

### A.1 Dephasing Channel Delta-D (Complete Derivation)

**Setup**: N_p: B(C^2) -> B(C^2) is the qubit dephasing channel with Kraus operators K_0 = sqrt((1+p)/2) I, K_1 = sqrt((1-p)/2) sigma_z.

In the Bloch ball representation:
```
N_p: (r_x, r_y, r_z) -> (p r_x, p r_y, r_z)
```

**Let rho = (1/2)(I + a sigma_x + b sigma_y + c sigma_z) with a^2 + b^2 + c^2 <= 1.**
**Let sigma = (1/2)(I + u sigma_x + v sigma_y + w sigma_z) be full-rank (u^2+v^2+w^2 < 1).**

**Eigenvalues**:
```
rho: lambda_pm = (1 pm r)/2 where r = sqrt(a^2 + b^2 + c^2)
sigma: mu_pm = (1 pm s)/2 where s = sqrt(u^2 + v^2 + w^2)
```

**N_p(rho)**: eigenvalues (1 pm r')/2 where r' = sqrt(p^2 a^2 + p^2 b^2 + c^2) = sqrt(p^2(a^2+b^2) + c^2)
**N_p(sigma)**: eigenvalues (1 pm s')/2 where s' = sqrt(p^2(u^2+v^2) + w^2)

**The relative entropy** D(rho || sigma) for qubit states is:
```
D(rho || sigma) = sum_{i,j} <e_i| rho |e_i> * |<e_i|f_j>|^2 * ln(<e_i|rho|e_i> / <f_j|sigma|f_j>)
```

This is complex for general non-commuting states. However, for the special case where rho and sigma share an eigenbasis (commuting case):

```
D(rho || sigma) = ((1+r)/2) ln((1+r)/(1+s)) + ((1-r)/2) ln((1-r)/(1-s))
```

For the pure state rho = |+><+| (r = 1, along x-axis) and sigma diagonal in z-basis (s = |w|, along z-axis):
These do NOT commute, and the computation requires the full formula.

However, as computed in Section 1.3, for sigma diagonal in the dephasing basis:
```
Delta-D = H_bin((1+p)/2)     [independent of sigma, proven above]
```

### A.2 Bosonic Loss Channel (Complete Derivation)

**The quantum relative entropy** for two Gaussian states on a single bosonic mode.

For coherent state |alpha> vs thermal state sigma_N:
```
D(|alpha><alpha| || sigma_N) = |alpha|^2 ln(1 + 1/N) + ln(1 + N)
```

This was derived in Section 2.2 using the Fock basis expansion.

**Verification**: In the limit N -> infinity, D -> |alpha|^2/N + ln N + ... which goes to infinity. In the limit |alpha| -> 0 (vacuum vs thermal):
```
D(|0><0| || sigma_N) = ln(1 + N)
```
which is the relative entropy between vacuum and thermal state.

**After the loss channel**:
```
E_eta(|alpha><alpha|) = |sqrt(eta) alpha><sqrt(eta) alpha|
E_eta(sigma_N) = sigma_{eta N}
```

The second equality follows from: the beam splitter maps a thermal state with mean N to a thermal state with mean eta N (when the other input port is vacuum).

Therefore:
```
D(E_eta(|alpha><alpha|) || E_eta(sigma_N)) = eta |alpha|^2 ln(1 + 1/(eta N)) + ln(1 + eta N)
```

**Entropy production**:
```
Sigma = |alpha|^2 [ln(1 + 1/N) - eta ln(1 + 1/(eta N))] + [ln(1+N) - ln(1+eta N)]
```

**Large N expansion**: Let x = 1/N -> 0.
```
ln(1 + x) - eta ln(1 + x/eta) = x - x^2/2 + ... - eta(x/eta - x^2/(2eta^2) + ...)
= x - x^2/2 - x + x^2/(2eta) + ...
= x^2(1/(2eta) - 1/2) + ...
= x^2(1-eta)/(2eta) + ...
```

So the first term is |alpha|^2 (1-eta)/(2eta N^2) + ... which vanishes as N -> infinity.

```
ln(1+N) - ln(1+eta N) = ln(N) + ln(1+1/N) - ln(eta N) - ln(1+1/(eta N))
= -ln(eta) + O(1/N)
```

Therefore Sigma = -ln(eta) + O(1/N). QED.

### A.3 Why No Finite-Dimensional Channel Suffices for Unbounded Sigma

**Theorem**: For any CPTP map N: M_d -> M_d' (finite dimensional) and any fixed full-rank reference state sigma, the entropy production Sigma(rho) = D(rho||sigma) - D(N(rho)||N(sigma)) satisfies:

```
sup_rho Sigma(rho) <= ln(d) + ln(1/lambda_min(sigma))
```

where d = dim of the input Hilbert space and lambda_min(sigma) is the smallest eigenvalue of sigma.

**Proof sketch**: D(rho || sigma) <= ln(d/lambda_min(sigma)) for all rho (standard bound on relative entropy). Since D(N(rho)||N(sigma)) >= 0, we get Sigma <= D(rho||sigma) <= ln(d/lambda_min(sigma)). QED.

**Consequence**: Since -ln(-g_00) = r_s/r can be arbitrarily large, NO finite-dimensional channel with a fixed reference state can reproduce Sigma_grav = -ln(-g_00) for all r. The gravitational channel, if it exists as a CPTP map, must act on an infinite-dimensional (bosonic) Hilbert space.

**This rules out**: all qubit models (Pikovski dephasing, qubit amplitude damping, etc.) as candidates for the exact gravitational channel. They can only serve as effective models in regimes where r_s/r < ln 2.

---

## 9. Open Problems and Directions

### 9.1 First-Principles Construction of N_grav

The most important open problem is to construct N_grav as an explicit CPTP map from first principles, without assuming a bosonic loss model. Possible approaches:

1. **Crossed-product algebra** (Witten 2022, Chandrasekaran 2023): The most promising rigorous framework. The conditional expectation from A(r) to A(inf) is a well-defined completely positive map. The challenge is computing its entropy production explicitly.

2. **Unruh-DeWitt detector model**: An explicit model of a quantum system interacting with a field in curved spacetime. The channel can be computed perturbatively, but it involves probe-dependent parameters.

3. **Functional integral approach**: Integrating out field modes between r and infinity defines a channel. This is related to the Wilsonian renormalization group, where DPI = RG monotonicity (Casini-Huerta 2017).

### 9.2 Saturation Problem

Prove or disprove: does ANY quantum channel N with Sigma > 0 saturate the JRSWW bound F = exp(-Sigma/2)?

If no such channel exists (likely), then the exponential metric is at best an APPROXIMATION to the true gravitational geometry, with the actual metric sitting between Schwarzschild and exponential.

### 9.3 The Reference State Problem

In thermodynamic applications of JRSWW, the reference state is naturally the thermal/Gibbs state. In the gravitational context, the natural reference is the Hartle-Hawking/KMS state. But this state depends on the background geometry (through the Tolman temperature), creating a potential circularity: to define Sigma_grav, one needs a reference state, which depends on the metric, which is supposed to be determined by Sigma_grav.

Breaking this circularity requires either:
(a) A reference state that does not depend on the metric (e.g., Minkowski vacuum), or
(b) A self-consistency argument (the metric is the fixed point of the map Sigma -> g_00 -> sigma -> Sigma)

### 9.4 Experimental Discrimination

The factor-of-2 ambiguity (Sigma = r_s/r or r_s/(2r)) corresponds to different physical observables:
- Sigma = r_s/r: predicts sqrt(-g_00) = F_bound (exponential metric if saturated)
- Sigma = r_s/(2r): predicts (-g_00)^{1/4} = F_bound (different metric)

These give different strong-field predictions that are in principle observable.

---

## 10. Summary Table

| Question | Answer | Confidence |
|----------|--------|------------|
| Does a CPTP map with Delta-D = -ln(-g_00) exist? | Yes, the bosonic loss channel | HIGH (conditional on modeling choices) |
| What channel? | Pure-loss bosonic, eta = -g_00 | MEDIUM (effective model, not first-principles) |
| What reference? | High-temperature thermal | MEDIUM (motivated by Tolman, requires N_B >> 1) |
| Is Sigma_grav = r_s/r? | Yes, via 3+ independent routes | MEDIUM-HIGH |
| Is the Petz bound saturated? | Unknown, likely not exactly | LOW (no proof for any Sigma > 0) |
| Is the connection a theorem? | Layer 1 yes, Layer 2 conjecture, Layer 3 speculation | -- |
| Can a qubit channel suffice? | NO (fundamental dimensional obstruction) | HIGH (proven) |
| Factor of 2 resolved? | Convention-dependent; eta = -g_00 gives r_s/r | MEDIUM |

---

**Last updated**: 2026-03-09
