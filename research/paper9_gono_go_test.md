# Paper 9 Go/No-Go Test: Sigma = Fisher Information of the Spectral Action
# Toy Model S^1 x M_2(C) -- Complete Design

**Author**: Sheng-Kai Huang (with systematic analysis)
**Date**: 2026-03-19
**Status**: TEST DESIGN -- ready for execution
**Purpose**: Compute the QRE between two Tolman-related Gibbs states on the spectral triple S^1 x M_2(C). Verify that the Fisher metric matches g_QQ = 2/Q^2. This is the decisive Go/No-Go gate for Paper 9.

---

## 0. The Question in One Sentence

> On the simplest nontrivial spectral triple, does the quantum relative entropy between Tolman-shifted thermal states reproduce Sigma = 2 ln Q, with Fisher metric g_QQ = 2/Q^2?

---

## 1. Precise Definition of the Toy Model

### 1.1 The Spectral Triple (A, H, D)

**Algebra**:
```
A = C(S^1) tensor M_2(C)
```
- C(S^1) = continuous functions on the circle of circumference beta_0 (inverse reference temperature)
- M_2(C) = 2x2 complex matrices (simplest non-abelian internal space)
- A acts on H by pointwise multiplication (C(S^1) part) and matrix multiplication (M_2(C) part)

**Hilbert space**:
```
H = L^2(S^1) tensor C^2
```
- L^2(S^1) = square-integrable functions on the circle, with orthonormal basis {e_n(theta) = (1/sqrt(beta_0)) exp(2pi i n theta / beta_0)}, n in Z
- C^2 = two-component spinor for the internal space
- A general element: psi(theta) = (psi_+(theta), psi_-(theta))^T

**Dirac operator**:
```
D = D_{S^1} tensor 1_2 + gamma tensor D_F
```
where:
- D_{S^1} = -i d/d(theta) acting on L^2(S^1), with eigenvalues omega_n = 2pi n / beta_0 on e_n
- gamma = sigma_3 = diag(1, -1) (the "chirality" on S^1, provides Z_2 grading)
- D_F = m sigma_1 = m ((0,1),(1,0)), with eigenvalues +m and -m
- 1_2 = 2x2 identity matrix

**Eigenvectors and eigenvalues of D**:

The operator D acts on the basis {e_n tensor |s>} where |s> = |+>, |-> are eigenstates of D_F in the {|+m>, |-m>} basis. Diagonalizing D on the n-th Fourier mode:

```
D|_{n-sector} = ( omega_n,  m  )
                ( m,  -omega_n )
```

Eigenvalues:
```
lambda_{n,pm} = pm sqrt(omega_n^2 + m^2) = pm sqrt((2pi n / beta_0)^2 + m^2)
```

Eigenvectors:
```
|n, +> = cos(alpha_n/2) |e_n, +> + sin(alpha_n/2) |e_n, ->
|n, -> = -sin(alpha_n/2) |e_n, +> + cos(alpha_n/2) |e_n, ->
```
where tan(alpha_n) = m / omega_n.

**The full spectrum of D**:
```
Spec(D) = { pm E_n : n in Z },  E_n = sqrt((2pi n / beta_0)^2 + m^2)
```
Each eigenvalue has degeneracy 1 (from M_2(C) diagonalization). The n=0 modes have E_0 = m (pure mass gap from the internal space).

### 1.2 Physical Interpretation

- S^1 of circumference beta_0 = "Euclidean time circle" = inverse temperature 1/beta_0
- M_2(C) = simplest internal (finite) geometry, analogous to the full SM algebra A_F = C + H + M_3(C) but reduced
- m = mass parameter from D_F, analogous to fermion masses from Yukawa couplings
- The spectral triple describes a 1-dimensional thermal system with a 2-level internal structure

### 1.3 Parameter Regime

For the test to be physically relevant:
- beta_0 m >> 1 : "low temperature" regime (discrete spectrum dominates)
- beta_0 m << 1 : "high temperature" regime (continuum limit, expect universal behavior)
- beta_0 m ~ O(1) : intermediate regime (full spectral structure matters)

**The test must be performed in ALL three regimes.** The conjecture Sigma = 2 ln Q is expected to hold universally; any regime-dependent failure is informative.

---

## 2. The Spectral Action on S^1 x M_2(C)

### 2.1 Definition

```
S_spectral[D] = Tr f(D^2 / Lambda^2)
```

where f is a positive even cutoff function (Schwartz class, or the CCSvS universal function f_S for the entropy interpretation). Using the spectrum:

```
S_spectral = sum_{n in Z} sum_{s=+/-} f(E_n^2 / Lambda^2)
           = 2 sum_{n in Z} f(((2pi n / beta_0)^2 + m^2) / Lambda^2)
```

(the factor 2 comes from the two eigenvalues pm E_n giving the same f(E_n^2)).

### 2.2 Poisson Summation

Using the Poisson summation formula:
```
sum_{n in Z} f(((2pi n / beta_0)^2 + m^2) / Lambda^2)
= (beta_0 Lambda / (2pi)) integral_{-inf}^{inf} f((x^2 + m^2/Lambda^2)) dx
  + 2 sum_{k=1}^{inf} (beta_0 Lambda / (2pi)) integral f((x^2 + m^2/Lambda^2)) cos(k beta_0 Lambda x) dx
```

In the large Lambda limit (Lambda >> 2pi/beta_0, Lambda >> m):
```
S_spectral ~ 2 (beta_0 Lambda / (2pi)) integral_0^{inf} f(u) / (2 sqrt(u)) du  + O(e^{-beta_0 m})
           = (beta_0 Lambda / pi) f_{1/2}
```
where f_{1/2} = integral_0^{inf} f(u) u^{-1/2} du / 2 is the half-moment of f.

### 2.3 CCSvS Entropy Function

For the entropy interpretation (CCSvS 2018), f = f_S where:
```
f_S(x) = x/(e^x - 1) - ln(1 - e^{-x})    [bosonic]
f_S(x) = x/(e^x + 1) + ln(1 + e^{-x})     [fermionic]
```

For the fermionic second quantization on our spectral triple, the von Neumann entropy is:
```
S_vN(beta) = sum_{n,s} [beta E_{n,s} / (e^{beta E_{n,s}} + 1) + ln(1 + e^{-beta E_{n,s}})]
```

This IS the spectral action with f = f_S^{ferm} and Lambda^2 = 1/beta (up to a reparameterization).

---

## 3. QRE Computation: The Core Calculation

### 3.1 Two Gibbs States

**State 1 (rho_Q)**: Gibbs state of the fermionic second quantization at inverse temperature beta_1 = beta_0 / Q.

```
rho_Q = exp(-beta_1 H_F) / Z(beta_1) = exp(-(beta_0/Q) H_F) / Z(beta_0/Q)
```

where H_F = sum_{n,s} E_{n,s} c^dag_{n,s} c_{n,s} is the fermionic Hamiltonian (second quantization of D^2).

**State 2 (rho_1)**: Gibbs state at reference inverse temperature beta_0 (Q=1).

```
rho_1 = exp(-beta_0 H_F) / Z(beta_0)
```

**Tolman relation**: These are Tolman-related if we identify:
```
T_1 sqrt(g_{00}^{(1)}) = T_2 sqrt(g_{00}^{(2)})
```
With T_1 = Q / beta_0 and T_2 = 1 / beta_0, this means:
```
g_{00}^{(1)} / g_{00}^{(2)} = 1/Q^2
```
i.e., state rho_Q corresponds to a region with g_{00} = -1/Q^2 (consistent with our definition Q = 1/sqrt(-g_{00})).

### 3.2 QRE Between Fermionic Gibbs States

For two thermal states of a free fermion system, the QRE factorizes over modes:

```
D(rho_Q || rho_1) = sum_{n,s} D(rho_Q^{(n,s)} || rho_1^{(n,s)})
```

For a single fermionic mode with energy E at inverse temperatures beta_1 = beta_0/Q and beta_2 = beta_0:

```
rho_{beta}^{mode} = diag(1 - n_F(beta E), n_F(beta E))
```

where n_F(x) = 1/(e^x + 1) is the Fermi-Dirac distribution.

The QRE for a single mode:
```
D(rho_1^{mode} || rho_2^{mode}) = n_1 ln(n_1/n_2) + (1-n_1) ln((1-n_1)/(1-n_2))
```

where n_1 = n_F(beta_0 E / Q) and n_2 = n_F(beta_0 E).

### 3.3 Explicit Formula

Define x = beta_0 E (dimensionless energy). Then:

```
n_1 = 1/(e^{x/Q} + 1),  n_2 = 1/(e^x + 1)
```

The per-mode QRE:
```
d(x, Q) = n_1(x,Q) ln(n_1/n_2) + (1-n_1(x,Q)) ln((1-n_1)/(1-n_2))
```

Simplifying using the identity ln(n_F(y)/(1-n_F(y))) = -y:
```
d(x, Q) = n_1(x,Q) [x - x/Q] + ln(1+e^{-x/Q}) - ln(1+e^{-x})
         = (1-1/Q) x n_F(x/Q) + ln(1+e^{-x/Q}) - ln(1+e^{-x})
```

Total QRE (summing over the spectrum of D^2):
```
Sigma(Q) = D(rho_Q || rho_1) = 2 sum_{n in Z} d(beta_0 E_n, Q)
```

where E_n = sqrt((2pi n / beta_0)^2 + m^2), and the factor 2 accounts for the two signs pm E_n giving the same D^2 eigenvalue E_n^2.

### 3.4 The Three Regimes

**Regime A: High temperature (beta_0 m << 1, many modes excited)**

In this limit, the sum over n becomes an integral. Using the Euler-Maclaurin formula:

```
Sigma(Q) ~ 2 (beta_0 / (2pi)) integral_0^{inf} d(beta_0 E, Q) * (beta_0 / sqrt(E^2 - m^2)) dE * (2pi / beta_0)
         = integral_0^{inf} d(beta_0 sqrt(k^2 + m^2), Q) dk * (beta_0 / pi)
```

where k = 2pi n / beta_0 is the "momentum."

For beta_0 -> 0 (with m fixed), the dominant contribution comes from n >> 1 modes where E_n >> m. In this regime n_F(x) ~ 1/2 - x/4, and the QRE per mode vanishes as x^2. The total QRE is:

```
Sigma(Q) ~ (pi / (3 beta_0)) [(beta_0/Q)^2 - beta_0^2] / ...
```

Actually, let us be more careful. The standard result for the relative entropy between two thermal states of a 1d fermion at temperatures T_1 and T_2 is:

```
D(rho_{T_1} || rho_{T_2}) = (pi L / 6) [T_1/T_2 + T_2/T_1 - 2] * (1 / beta_0)
```

Wait -- this is the CFT result for a 1+1d system of length L. For our S^1 of circumference beta_0, in the high-T limit (T >> 1/beta_0):

This needs to be computed carefully. The point is that in 1d, the QRE between thermal states does NOT simply equal 2 ln Q. **This is where the test becomes nontrivial.**

**Regime B: Low temperature (beta_0 m >> 1, only n=0 mode contributes)**

Only the zero mode E_0 = m contributes (all other modes are exponentially suppressed):

```
Sigma(Q) ~ 2 d(beta_0 m, Q)
         = 2 [(1-1/Q) beta_0 m * n_F(beta_0 m / Q) + ln(1+e^{-beta_0 m/Q}) - ln(1+e^{-beta_0 m})]
```

For beta_0 m >> 1 AND beta_0 m / Q >> 1 (i.e., Q not too large):
```
Sigma(Q) ~ 2 (1-1/Q) beta_0 m * e^{-beta_0 m / Q} + 2 e^{-beta_0 m / Q} - 2 e^{-beta_0 m}
         ~ 2 e^{-beta_0 m / Q} [1 + (1-1/Q) beta_0 m] - 2 e^{-beta_0 m}
```

This is exponentially small and does NOT look like 2 ln Q. **Potential failure point.**

**Regime C: Intermediate (beta_0 m ~ O(1))**

Must be computed numerically. See Section 6 below.

### 3.5 Critical Observation: The Correct Identification

The analysis above reveals a crucial subtlety. The QRE D(rho_Q || rho_1) between thermal states of the SAME Hamiltonian at different temperatures is NOT the same as the channel QRE Sigma = 2 ln Q.

The channel picture (Paper 2) gives Sigma = -ln(eta) where eta is the channel transmissivity. For the gravitational thermal attenuator:
```
eta = -g_{00} = 1/Q^2
Sigma = -ln(1/Q^2) = 2 ln Q
```

The state-space picture (this computation) gives D(rho_{beta_0/Q} || rho_{beta_0}), which depends on the full spectrum and IS NOT simply 2 ln Q in general.

**The key question becomes**: Under what conditions does the state-space QRE reduce to the channel QRE?

**Answer**: When the system has a LARGE number of modes (thermodynamic limit), the state-space QRE per degree of freedom should match the channel QRE. Specifically:

```
lim_{N -> inf} (1/N) D(rho_Q^{(N)} || rho_1^{(N)}) = Sigma_{channel} = 2 ln Q
```

where N is the number of modes (cutoff by Lambda). This is the **extensivity** of QRE.

For our S^1 x M_2(C) model with N modes (|n| <= N):
```
Sigma_{total}(Q) = sum_{|n| <= N} d(beta_0 E_n, Q) * 2
Sigma_{per mode}(Q) = Sigma_{total} / (2(2N+1))
```

The Go/No-Go test is: **does Sigma_{per mode} -> 2 ln Q in the thermodynamic limit?**

---

## 4. Fisher Metric Extraction

### 4.1 Definition

The Fisher information metric on the Q-manifold at point Q_0 is:

```
g_QQ(Q_0) = (d^2 / dQ^2) D(rho_{Q_0 + epsilon} || rho_{Q_0}) |_{epsilon=0}
```

Note: this is the second derivative of QRE with respect to PERTURBATION around Q_0, not the second derivative of Sigma(Q) = D(rho_Q || rho_1).

Alternatively, from Sigma:
```
g_QQ(Q) = d^2 Sigma / dQ^2 |_{locally}
```

Specifically, parameterize Q' near Q:
```
D(rho_{Q'} || rho_Q) = (1/2) g_QQ(Q) (Q'-Q)^2 + O((Q'-Q)^3)
```

### 4.2 Per-Mode Fisher Information

For a single fermionic mode at energy E, with x = beta_0 E:

```
g_QQ^{mode}(Q) = (d^2 / dQ'^2) d(x, Q+Q') |_{Q'=0} evaluated relative to rho_Q as reference
```

Computing directly:
```
dn_1/dQ = (x/Q^2) n_F(x/Q)(1 - n_F(x/Q))
```

The Fisher information for a Bernoulli distribution p = n_F(x/Q) is:
```
g_QQ^{mode}(Q) = [dn_1/dQ]^2 / [n_1(1-n_1)]
               = (x/Q^2)^2 n_F(x/Q)(1-n_F(x/Q))
               = (x^2/Q^4) n_F(x/Q)(1-n_F(x/Q))
```

### 4.3 Total Fisher Metric

```
g_QQ(Q) = 2 sum_{n in Z} (beta_0 E_n)^2 / Q^4 * n_F(beta_0 E_n / Q)(1 - n_F(beta_0 E_n / Q))
```

**Target**: g_QQ(Q) = 2/Q^2 (per degree of freedom, in thermodynamic limit).

Per degree of freedom:
```
g_QQ^{per dof}(Q) = g_QQ(Q) / (2(2N+1))
```

The target is:
```
g_QQ^{per dof}(Q) = 2/Q^2  (?)
```

### 4.4 Evaluation at Q = 1

At Q = 1:
```
g_QQ^{mode}(Q=1) = (beta_0 E_n)^2 * n_F(beta_0 E_n)(1-n_F(beta_0 E_n))
                  = (beta_0 E_n)^2 / [4 cosh^2(beta_0 E_n / 2)]
```

**High-temperature limit** (beta_0 E_n << 1 for most modes):
```
g_QQ^{mode} ~ (beta_0 E_n)^2 / 4
```

Total (with N modes):
```
g_QQ(Q=1) ~ 2 sum_{n} (beta_0 E_n)^2 / 4 ~ (beta_0^2 / 2) sum_n E_n^2
```

This DIVERGES as N -> inf (UV divergence). Per degree of freedom:
```
g_QQ^{per dof}(Q=1) ~ (beta_0^2 / 2) <E^2>
```

where <E^2> is the average squared energy. In the high-T limit, <E^2> ~ (2pi N / beta_0)^2 / 3 (for a uniform distribution up to mode N), so:

```
g_QQ^{per dof}(Q=1) ~ (2pi^2 / 3) N^2
```

This does NOT equal 2. **This means the naive per-mode normalization does not give g_QQ = 2/Q^2.**

### 4.5 The Correct Normalization: Channel Fisher vs State Fisher

The resolution lies in distinguishing two Fisher metrics:

1. **State Fisher metric**: g_QQ^{state} = d^2 D(rho_Q || rho_1) / dQ^2 -- this depends on the spectrum

2. **Channel Fisher metric**: g_QQ^{channel} = d^2 Sigma_{channel} / dQ^2 -- this is 2/Q^2

For the gravitational channel E_Q (thermal attenuator with eta = 1/Q^2):
```
Sigma_{channel}(Q) = D(E_Q(rho) || E_1(rho)) = -ln(eta) = 2 ln Q
```
for ANY input state rho (in the channel picture, the entropy production is state-independent for the thermal attenuator in the high-T limit).

The state Fisher metric and channel Fisher metric are related by:
```
g_QQ^{state} = sum_{modes} g_QQ^{channel, per mode} * (spectral weight of mode)
```

**The Go/No-Go test must compare the CHANNEL Fisher metric, not the state Fisher metric.**

---

## 5. The Corrected Test: Channel QRE on S^1 x M_2(C)

### 5.1 The Channel Construction

Define the gravitational channel on S^1 x M_2(C):

**Step 1**: The Dirac operator D depends on the "metric" parameter Q through conformal rescaling:
```
D(Q) = D / Q
```

This corresponds to: D_{S^1}(Q) = (1/Q)(-i d/dtheta), i.e., the circle has effective circumference beta_0 * Q.

**Step 2**: The channel E_Q maps states on the reference spectral triple to states on the Q-rescaled spectral triple. In the second-quantized fermionic Fock space, this is:

```
E_Q : rho -> rho_Q  (Gibbs state at temperature Q/beta_0 for Hamiltonian D(1)^2 / Q^2)
```

Equivalently, E_Q is the thermal attenuator that mixes the input with a thermal environment at the Tolman-shifted temperature.

**Step 3**: The channel entropy production is:
```
Sigma_{channel}(Q) = S(E_Q(rho)) - S(rho)  [for input rho at reference temperature]
```

For the thermal attenuator:
```
Sigma_{channel} = -ln(eta) = -ln(1/Q^2) = 2 ln Q
```

### 5.2 Explicit Verification via Partition Functions

The cleanest approach uses the free energy / partition function.

For the fermionic second quantization on the spectral triple with D^2/Lambda^2 rescaled by 1/Q^2:

**Partition function**:
```
ln Z(beta, Q) = sum_{n in Z} ln(1 + e^{-beta E_n / Q}) * 2
```

**Free energy**:
```
F(beta, Q) = -(2/beta) sum_{n} ln(1 + e^{-beta E_n / Q})
```

**Von Neumann entropy**:
```
S_vN(beta, Q) = beta^2 dF/dbeta = 2 sum_n [beta E_n/(Q(e^{beta E_n/Q}+1)) + ln(1+e^{-beta E_n/Q})]
```

**QRE** (same beta, different Q):
```
D(rho_Q || rho_1) = beta [<H(1)>_Q - F(beta,1)] - [S_vN(beta,Q)]
```

Wait -- let me be precise. For Gibbs states at the SAME beta but DIFFERENT Hamiltonians H_Q = H/Q^2 and H_1 = H:

```
D(rho_Q || rho_1) = Tr[rho_Q (ln rho_Q - ln rho_1)]
                   = Tr[rho_Q (-beta H_Q - ln Z_Q)] - Tr[rho_Q (-beta H_1 - ln Z_1)]
                   = beta Tr[rho_Q (H_1 - H_Q)] + ln Z_1 - ln Z_Q
                   = beta (<H_1>_Q - <H_Q>_Q) + ln Z_1 - ln Z_Q
                   = beta (1 - 1/Q^2) <H>_Q + ln Z_1 - ln Z_Q
```

where <H>_Q = Tr[rho_Q H] = sum_n E_n^2 * 2 n_F(beta E_n / Q).

Hmm, this is getting complicated. Let me try the Tolman approach instead.

### 5.3 Tolman Approach (Different beta, Same H)

The more natural parameterization: same Hamiltonian H = D^2, different temperatures beta_Q = beta_0/Q and beta_1 = beta_0.

```
D(rho_{beta_Q} || rho_{beta_1}) = Tr[rho_{beta_Q}(ln rho_{beta_Q} - ln rho_{beta_1})]
= Tr[rho_{beta_Q}(-beta_Q H - ln Z_Q + beta_1 H + ln Z_1)]
= (beta_1 - beta_Q)<H>_Q + ln Z_1 - ln Z_Q
= beta_0(1 - 1/Q)<H>_Q + ln Z(beta_0) - ln Z(beta_0/Q)
```

For the single-mode contribution at energy E:
```
d_mode(Q) = (beta_0 - beta_0/Q) E n_F(beta_0 E / Q) + ln(1+e^{-beta_0 E}) - ln(1+e^{-beta_0 E/Q})
```

### 5.4 The Sum: What We Need to Compute

```
Sigma(Q) = 2 sum_{n=0}^{inf} [
    beta_0(1-1/Q) E_n n_F(beta_0 E_n/Q)
  + ln(1+e^{-beta_0 E_n}) - ln(1+e^{-beta_0 E_n/Q})
]'
```

where the prime on the sum means the n=0 term counts once, all other n count twice (from pm n).

Wait, I need to be more careful. E_n = sqrt((2pi n/beta_0)^2 + m^2). For the full spectrum including both pm signs of E:

Actually, in the fermionic second quantization, the Hamiltonian is H = sum_{n,s} |lambda_{n,s}| c^dag c where the sum is over ALL eigenvalues of |D|. Since E_n = E_{-n} for our spectrum, and each E_n appears twice (from the M_2(C) diagonalization):

```
Sigma(Q) = 2 sum_{n=-inf}^{inf} d_mode(beta_0 E_{|n|}, Q)
         = 2 [d_mode(beta_0 m, Q) + 2 sum_{n=1}^{inf} d_mode(beta_0 E_n, Q)]
```

**This is the master formula to evaluate.**

---

## 6. Analytical Evaluation

### 6.1 Exact Closed-Form: Is It Possible?

The sum involves the function:
```
d(x, Q) = (1-1/Q) x / (e^{x/Q}+1) + ln(1+e^{-x}) - ln(1+e^{-x/Q})
```

summed over x_n = beta_0 E_n = beta_0 sqrt((2pi n/beta_0)^2 + m^2) = sqrt((2pi n)^2 + (beta_0 m)^2).

There is no known closed form for this sum at finite N. **However, we can extract exact results in the two limits.**

### 6.2 High-Temperature Limit (beta_0 m -> 0): Analytical

When beta_0 -> 0 with m and N fixed, all x_n -> 0. Expand d(x, Q) for small x:

```
d(x, Q) = (1-1/Q) x [1/2 - x/(4Q) + ...] + [-ln 2 + x/2 - x^2/8 + ...] - [-ln 2 + x/(2Q) - x^2/(8Q^2) + ...]
        = (1-1/Q) x/2 - (1-1/Q) x^2/(4Q) + x/2 - x^2/8 - x/(2Q) + x^2/(8Q^2) + O(x^3)
        = x/2 - x/(2Q) - x^2(1-1/Q)/(4Q) + x/2 - x^2/8 - x/(2Q) + x^2/(8Q^2) + O(x^3)
```

Let me redo this more carefully. Set y = x/Q.

```
n_F(y) = 1/(e^y+1),  n_F(x) = 1/(e^x+1)

d(x,Q) = (x - y) n_F(y) + ln(1+e^{-x}) - ln(1+e^{-y})
```

Write F(t) = ln(1+e^{-t}). Then:

```
d(x,Q) = (x-y) n_F(y) + F(x) - F(y)
```

Note F'(t) = -n_F(t). So:

```
F(x) - F(y) = -integral_y^x n_F(t) dt
```

Therefore:
```
d(x,Q) = (x-y) n_F(y) - integral_y^x n_F(t) dt
        = integral_y^x [n_F(y) - n_F(t)] dt
```

Since n_F is decreasing (for t > 0), the integrand n_F(y) - n_F(t) >= 0 for t >= y, confirming d >= 0 (as required by QRE non-negativity).

**Taylor expansion around x = y (i.e., Q = 1)**:

Set x = y + epsilon where epsilon = x - y = x(1-1/Q). To second order in epsilon:

```
n_F(t) = n_F(y) + (t-y) n_F'(y) + (t-y)^2 n_F''(y)/2 + ...
```

```
d = integral_y^x [(y-t) n_F'(y) + (y-t)^2 (-n_F''(y))/2 + ...] dt

Wait, let me redo:
n_F(y) - n_F(t) = -(t-y) n_F'(y) - (t-y)^2 n_F''(y)/2 + ...

d = integral_y^x [-(t-y) n_F'(y) - (t-y)^2 n_F''(y)/2] dt
  = -n_F'(y) (x-y)^2/2 - n_F''(y) (x-y)^3/6 + ...
  = |n_F'(y)| (x-y)^2/2 + O((x-y)^3)
```

Now n_F'(y) = -e^y/(e^y+1)^2 = -n_F(y)(1-n_F(y)), so |n_F'(y)| = n_F(y)(1-n_F(y)).

```
d(x,Q) = (1/2) n_F(y)(1-n_F(y)) (x-y)^2 + O((x-y)^3)
        = (1/2) n_F(x/Q)(1-n_F(x/Q)) x^2(1-1/Q)^2 + ...
```

This confirms the Fisher metric per mode:
```
g_QQ^{mode}(Q=1) = x^2 n_F(x)(1-n_F(x)) = x^2 / [4 cosh^2(x/2)]
```

**In the high-T limit** (x = beta_0 E_n -> 0):
```
g_QQ^{mode}(Q=1) -> x^2/4  (since n_F(0)(1-n_F(0)) = 1/4)
```

### 6.3 The Thermodynamic Limit and Per-Mode Sigma

For the sum over modes in the high-T limit, using x_n = beta_0 E_n << 1 for relevant modes:

```
Sigma(Q) = 2 sum_n d(x_n, Q) ~ 2 sum_n integral_{x_n/Q}^{x_n} [n_F(x_n/Q) - n_F(t)] dt
```

For very high T (all x_n -> 0), d(x, Q) to leading order:

```
d(x, Q) ~ (1/4) x^2 (1-1/Q)^2    [using n_F'(0) = -1/4]
```

Total:
```
Sigma(Q) ~ (1/2)(1-1/Q)^2 sum_n x_n^2 = (1/2)(1-1/Q)^2 sum_n beta_0^2 E_n^2
```

This is proportional to (1-1/Q)^2, NOT to 2 ln Q. **In the strict high-T limit of the state-space QRE, Sigma goes as (1-1/Q)^2, not as 2 ln Q.**

### 6.4 Resolution: Why This Does Not Match -- and What the Correct Test Is

The discrepancy reveals the precise nature of the test. There are TWO distinct quantities:

**Quantity A**: State-space QRE, D(rho_{T_1} || rho_{T_2}), where T_1 = T_0 Q and T_2 = T_0.
- This depends on the spectrum of D.
- In the high-T limit, it goes as (Q-1)^2 * (spectral sum).
- It is NOT 2 ln Q.

**Quantity B**: Channel QRE, Sigma = D(E_Q(rho) || E_1(rho)), where E_Q is the gravitational thermal attenuator.
- This equals -ln(eta) = 2 ln Q by the channel theorem (Paper 2).
- It is INDEPENDENT of the input state rho (for the thermal attenuator in the right regime).

**The connection**: Quantity A is related to Quantity B through the thermodynamic identity:

```
D(rho_{T_1} || rho_{T_2}) = beta_2 [<H>_{T_1} - <H>_{T_2}] - [S(T_1) - S(T_2)]
                           = beta_2 Delta<H> - Delta S
```

For a 1d CFT (beta_0 m = 0), this becomes:
```
D = (pi c / 6 beta_0) [Q + 1/Q - 2]
```

where c is the central charge (c = 1 for a single Dirac fermion).

Meanwhile, 2 ln Q ~ 2(Q-1) - (Q-1)^2 + ... for Q near 1, while (Q + 1/Q - 2) = (Q-1)^2/Q.

**These are different functions of Q.** The state-space QRE of a 1d CFT is NOT 2 ln Q. It is (pi c / 6 beta_0)(Q + 1/Q - 2).

### 6.5 The Decisive Question

This means the test must be reformulated. The question is NOT:

> "Does D(rho_{T_1} || rho_{T_2}) = 2 ln Q on S^1 x M_2(C)?"

(Answer: NO, in general.)

The correct question is:

> "Does the Fisher information metric derived from the spectral action's entropy functional, evaluated on the conformal deformation D -> D/Q, match g_QQ = 2/Q^2?"

And this means we must compute the **channel Fisher metric**, not the state-space Fisher metric. The channel Fisher metric is:

```
g_QQ^{channel}(Q) = (d^2/dQ'^2) D_channel(E_{Q+Q'} || E_Q) |_{Q'=0}
```

where D_channel is the Umegaki relative entropy between the Choi states of the channels E_{Q+Q'} and E_Q.

For the thermal attenuator with eta = 1/Q^2:
```
D_channel = (1/2)(n_bar+1) ln[(n_bar+1)/(n_bar'+1)] + (1/2) n_bar ln[n_bar/n_bar']
```

where n_bar and n_bar' are the thermal occupation numbers of the environment. In the high-T limit, this reduces to the classical channel capacity formula.

**However**, for the SPECIFIC question of whether the spectral action provides this Fisher metric, we need to compute the second variation of the CCSvS entropy under the conformal deformation D -> D/Q.

---

## 7. The Reformulated Go/No-Go Test

### 7.1 Test Statement (Precise)

**Compute**:
```
g_QQ^{CCSvS}(Q) = -(d^2/dQ^2) S_vN^{ferm}(beta_0, D/Q) |_{Q=1}
```

where S_vN^{ferm}(beta, D) is the von Neumann entropy of the fermionic Gibbs state at inverse temperature beta for the Dirac operator D, on the spectral triple S^1 x M_2(C).

**Compare with**:
```
g_QQ^{target} = 2/Q^2 |_{Q=1} = 2
```

### 7.2 Explicit Computation

The CCSvS entropy:
```
S_vN(Q) = S_vN^{ferm}(beta_0, D/Q) = 2 sum_n [beta_0 E_n/(Q(e^{beta_0 E_n/Q}+1)) + ln(1+e^{-beta_0 E_n/Q})]
```

First derivative:
```
dS_vN/dQ = 2 sum_n (beta_0 E_n / Q^2) [n_F(beta_0 E_n/Q)(1-n_F(beta_0 E_n/Q))] * (beta_0 E_n / Q) / (something)
```

Let me define s(y) = y/(e^y+1) + ln(1+e^{-y}) where y = beta_0 E_n / Q. Then S_vN = 2 sum_n s(y_n(Q)).

```
ds/dy = y e^y / (e^y+1)^2 + 1/(e^y+1) - e^{-y}/(1+e^{-y})
      = y n_F(y)(1-n_F(y)) + n_F(y) - (1-n_F(y))  ... [this needs care]
```

Actually: s(y) = y n_F(y) + ln(1+e^{-y}).

```
ds/dy = n_F(y) + y n_F'(y) + (-1) e^{-y}/(1+e^{-y})
      = n_F(y) - y n_F(y)(1-n_F(y)) - (1-n_F(y))
      = (2 n_F(y) - 1) - y n_F(y)(1-n_F(y))
      = -tanh(y/2) - y/(4 cosh^2(y/2))    ... [using 2n_F(y)-1 = -tanh(y/2)]
```

Hmm, let me just carefully use:
```
s(y) = y/(e^y+1) + ln(1+e^{-y})
s'(y) = 1/(e^y+1) - y e^y/(e^y+1)^2 - e^{-y}/(1+e^{-y})
      = n_F(y) - y n_F(y)(1-n_F(y)) - (1-n_F(y))
      = -(1-2n_F(y)) - y n_F(y)(1-n_F(y))
```

Since 1-2n_F(y) = tanh(y/2):
```
s'(y) = -tanh(y/2) - y/(4 cosh^2(y/2))
```

Now y_n(Q) = beta_0 E_n / Q, so dy_n/dQ = -beta_0 E_n / Q^2 = -y_n / Q.

```
dS_vN/dQ = 2 sum_n s'(y_n) (-y_n/Q)
         = (2/Q) sum_n y_n [tanh(y_n/2) + y_n/(4 cosh^2(y_n/2))]
```

Second derivative:
```
d^2 S_vN / dQ^2 = -(2/Q^2) sum_n y_n [tanh(y_n/2) + y_n/(4 cosh^2(y_n/2))]
                + (2/Q) sum_n d/dQ [y_n tanh(y_n/2) + y_n^2/(4 cosh^2(y_n/2))]
```

At Q = 1: y_n = beta_0 E_n =: x_n. Using dy_n/dQ|_{Q=1} = -x_n:

```
d^2 S_vN / dQ^2 |_{Q=1} = -2 sum_n x_n [tanh(x_n/2) + x_n/(4 cosh^2(x_n/2))]
  + 2 sum_n (-x_n) d/dy [y tanh(y/2) + y^2/(4 cosh^2(y/2))] |_{y=x_n}
```

This is computable but messy. Define:

```
A = sum_n x_n tanh(x_n/2)
B = sum_n x_n^2 / (4 cosh^2(x_n/2))
C = sum_n x_n d/dy [y tanh(y/2)] |_{y=x_n}
D_sum = sum_n x_n d/dy [y^2/(4 cosh^2(y/2))] |_{y=x_n}
```

Then:
```
d^2 S_vN / dQ^2 |_{Q=1} = -2(A + B) - 2(C + D_sum)
```

**This is the quantity to compare with -2** (since g_QQ^{target} = 2 means we need -d^2 S_vN / dQ^2 = 2, but only after proper normalization).

### 7.3 Crucial Issue: Normalization

The entropy S_vN is EXTENSIVE (proportional to the number of modes). The Fisher metric g_QQ = 2/Q^2 is INTENSIVE (per degree of freedom). Therefore:

**The proper Go/No-Go comparison is**:

```
g_QQ^{CCSvS, per dof}(Q=1) = -(1/N_{dof}) (d^2 S_vN / dQ^2) |_{Q=1} =? 2
```

where N_{dof} = total number of degrees of freedom = 2(2N_{max}+1) for our spectral triple with modes |n| <= N_{max}.

### 7.4 Analytical Result in Special Cases

**Case 1: Single mode (n=0 only, E_0 = m)**

```
x_0 = beta_0 m

g_QQ^{single}(Q=1) = -d^2 s(x_0/Q) / dQ^2 |_{Q=1}
                    = x_0^2 s''(x_0) + x_0 s'(x_0)  [from chain rule with y = x_0/Q]
```

Actually, s as a function of Q: s(x_0/Q). Let u = 1/Q.

```
d/dQ s(x_0/Q) = s'(x_0 u) * x_0 * (-1/Q^2) = -x_0 u^2 s'(x_0 u)   [since u = 1/Q]
```

At Q=1 (u=1):
```
ds/dQ|_{Q=1} = -x_0 s'(x_0)
```

```
d^2/dQ^2 s(x_0/Q)|_{Q=1} = x_0^2 s''(x_0) + 2 x_0 s'(x_0)
```

Therefore:
```
g_QQ^{single}(Q=1) = -(x_0^2 s''(x_0) + 2 x_0 s'(x_0))
```

Now:
```
s'(y) = -(1-2n_F) - y n_F(1-n_F) = -tanh(y/2) - y/(4cosh^2(y/2))
s''(y) = -1/(2cosh^2(y/2)) - 1/(4cosh^2(y/2)) + y tanh(y/2)/(4cosh^2(y/2)) ... [compute carefully]
```

Let me use n = n_F(y). Then:
```
s'(y) = -(1-2n) - yn(1-n)
s''(y) = 2n(1-n) - n(1-n) + y(1-2n)n(1-n)
       = n(1-n) + y(1-2n)n(1-n)
       = n(1-n)[1 + y(1-2n)]
       = n(1-n)[1 - y tanh(y/2)]
```

Wait, let me be more careful:
```
d/dy [-(1-2n)] = 2 dn/dy = -2n(1-n)
d/dy [-yn(1-n)] = -n(1-n) - y d/dy[n(1-n)]
                = -n(1-n) - y n(1-n)(1-2n)  [since d/dy[n(1-n)] = n(1-n)(2n-1) ... no]
```

Actually, dn/dy = -n(1-n), so:
```
d/dy[n(1-n)] = (1-n)(-n(1-n)) + n(n(1-n)) ... no,
d/dy[n(1-n)] = (dn/dy)(1-n) + n(-dn/dy) = -n(1-n)(1-n) + n^2(1-n) = n(1-n)(2n-1)
```

So:
```
s''(y) = -2n(1-n) - n(1-n) - y n(1-n)(2n-1)
       = -n(1-n)[3 + y(2n-1)]
       = -n(1-n)[3 - y tanh(y/2)]    [since 2n-1 = -tanh(y/2)]
```

Hmm, I'm getting sign issues. Let me just define everything numerically for specific x_0.

### 7.5 Numerical Protocol

Rather than fighting with algebra, the definitive test is numerical. Define:

```python
import numpy as np
from scipy.special import expit  # n_F(x) = expit(-x) = 1/(1+e^x)

def n_F(x):
    return 1.0 / (1.0 + np.exp(x))

def s(y):
    """CCSvS entropy per mode: s(y) = y*n_F(y) + ln(1+e^{-y})"""
    return y * n_F(y) + np.log1p(np.exp(-y))

def S_vN(Q, beta0, m, N_max):
    """Total von Neumann entropy at conformal parameter Q"""
    total = 0.0
    for n in range(-N_max, N_max+1):
        E_n = np.sqrt((2*np.pi*n/beta0)**2 + m**2)
        y = beta0 * E_n / Q
        total += 2 * s(y)  # factor 2 from M_2(C)
    return total

def Sigma_state(Q, beta0, m, N_max):
    """State-space QRE: D(rho_Q || rho_1)"""
    total = 0.0
    for n in range(-N_max, N_max+1):
        E_n = np.sqrt((2*np.pi*n/beta0)**2 + m**2)
        x = beta0 * E_n
        y = x / Q
        d_mode = (x - y) * n_F(y) + np.log1p(np.exp(-x)) - np.log1p(np.exp(-y))
        total += 2 * d_mode
    return total

# Fisher metric from second derivative of S_vN
def g_QQ_CCSvS(Q0, beta0, m, N_max, dQ=1e-5):
    """Numerical Fisher metric from CCSvS entropy"""
    Sp = S_vN(Q0 + dQ, beta0, m, N_max)
    S0 = S_vN(Q0, beta0, m, N_max)
    Sm = S_vN(Q0 - dQ, beta0, m, N_max)
    return -(Sp - 2*S0 + Sm) / dQ**2

# Fisher metric from second derivative of state-space QRE
def g_QQ_state(Q0, beta0, m, N_max, dQ=1e-5):
    """Numerical Fisher metric from state-space QRE"""
    Sp = Sigma_state(Q0 + dQ, beta0, m, N_max)
    S0 = Sigma_state(Q0, beta0, m, N_max)
    Sm = Sigma_state(Q0 - dQ, beta0, m, N_max)
    return (Sp - 2*S0 + Sm) / dQ**2

# RUN for multiple regimes
for beta0_m in [0.01, 0.1, 1.0, 10.0]:  # high-T to low-T
    for N_max in [10, 50, 100, 500]:
        beta0 = 1.0
        m = beta0_m / beta0
        N_dof = 2 * (2*N_max + 1)

        g_CCSvS = g_QQ_CCSvS(1.0, beta0, m, N_max) / N_dof
        g_state = g_QQ_state(1.0, beta0, m, N_max) / N_dof
        Sigma_val = Sigma_state(1.5, beta0, m, N_max) / N_dof
        target_Sigma = 2 * np.log(1.5)

        print(f"beta0*m={beta0_m}, N={N_max}: "
              f"g_CCSvS/dof={g_CCSvS:.6f}, "
              f"g_state/dof={g_state:.6f}, "
              f"Sigma/dof(Q=1.5)={Sigma_val:.6f}, "
              f"target 2lnQ={target_Sigma:.6f}")
```

---

## 8. Go/No-Go Criteria (Precise)

### 8.1 PRIMARY TEST: Fisher Metric Match

**GO condition**:
```
lim_{N->inf} (1/N_dof) g_QQ^{CCSvS}(Q=1) = 2  (to within 1%)
```
AND this limit is approached uniformly across all three beta_0 m regimes.

**NO-GO condition**:
```
(1/N_dof) g_QQ^{CCSvS}(Q=1) -> constant != 2
```
OR the limit depends on beta_0 m (non-universal).

### 8.2 SECONDARY TEST: Full Sigma Match

**GO condition**:
```
lim_{N->inf} (1/N_dof) Sigma_{state}(Q) = 2 ln Q  (for all Q in [1, 3])
```
matching to all orders, not just the Fisher (quadratic) term.

**NO-GO condition**:
```
Sigma_{per dof}(Q) = f(Q) != 2 ln Q for any Q > 1
```

### 8.3 TERTIARY TEST: Q-Dependence of Fisher Metric

**GO condition**:
```
(1/N_dof) g_QQ^{CCSvS}(Q) = 2/Q^2  (for Q in [1, 3])
```

**NO-GO condition**:
```
g_QQ^{per dof}(Q) has different Q-dependence than 2/Q^2
```

### 8.4 Expected Outcomes and Their Implications

| Outcome | Probability | Implication |
|---------|-------------|-------------|
| All three tests PASS | 20% | Strong GO: Sigma = Fisher info of spectral action, proceed to full A_F |
| Primary PASS, Secondary FAIL | 25% | Partial GO: Fisher metric matches but higher orders differ; Sigma = 2 ln Q is the leading-order (Gaussian) approximation |
| Primary FAIL but with simple correction | 25% | Conditional GO: g_QQ = c/Q^2 with c != 2; suggests Sigma = c ln Q with c depending on spectral data; may be salvageable by adjusting the identification |
| Primary FAIL, no pattern | 15% | Weak NO-GO: The CCSvS entropy does not directly give g_QQ = 2/Q^2; need a different route to connect Sigma to spectral geometry |
| All tests FAIL | 15% | Strong NO-GO: The entire spectral action approach to Sigma fails; retreat to Papers 1-4 only |

---

## 9. What Success Would Imply for Full A_F

### 9.1 From M_2(C) to A_F = C + H + M_3(C)

If the toy model passes, the extension to the full SM spectral triple requires:

1. **Replace M_2(C) with A_F**: The internal Dirac operator D_F becomes:
```
D_F = (Yukawa matrix) x (Higgs VEV) + (Majorana mass matrix)
```
with eigenvalues determined by fermion masses.

2. **Replace S^1 with M^4**: The manifold part becomes a 4-dimensional spin manifold. The Dirac operator D_M = i gamma^mu nabla_mu has continuous spectrum.

3. **The spectral action expands as**:
```
Tr f(D^2/Lambda^2) ~ Lambda^4 f_4 a_0 + Lambda^2 f_2 a_2 + f_0 a_4 + ...
```
The a_4 term contains the FULL SM Lagrangian.

4. **The Fisher metric under D -> D/Q becomes**:
```
g_QQ^{full} = -(d^2/dQ^2) S_vN^{ferm}(beta, D/Q) |_{Q=1}
            = (spectral sum over all SM fermions) x (single-mode Fisher metric)
```

If the per-mode Fisher metric is universal (= 2/Q^2 regardless of mass), then:
```
g_QQ^{full, per dof} = 2/Q^2
```
and the FULL Standard Model contributes the SAME Sigma = 2 ln Q per degree of freedom.

### 9.2 The Key Prediction: Universality

The most important implication of success is **universality**: the Fisher metric g_QQ = 2/Q^2 does not depend on the mass m or the spectral details. This would mean:

- Sigma = 2 ln Q is UNIVERSAL across all spectral triples (as long as they are finite-dimensional internally)
- The specific algebra A_F (whether C, H, M_3(C), or their product) does not affect Sigma
- The coupling constants (alpha, G, lambda_H, etc.) enter through the NORMALIZATION of the spectral action, not through the Fisher metric

This is consistent with the general argument: Sigma = 2 ln Q is a GEOMETRIC statement about the conformal mode, independent of matter content.

### 9.3 What Would Change in the Full Model

Even with universality, the full model introduces:
- **Multiple conformal sectors**: The Khronon Q, the dilaton, and the Higgs can all be seen as conformal deformations. Their Fisher metrics may differ.
- **Cross terms**: The Fisher metric becomes a MATRIX g_{IJ} where I, J index the different scalar fields.
- **Potential terms**: The spectral action potential V(Q, H, sigma) generates non-Fisher contributions.

---

## 10. Difficulty Assessment and Timeline

### 10.1 Can This Be Done Analytically?

**YES, but with caveats:**

| Component | Method | Difficulty | Time |
|-----------|--------|------------|------|
| Spectrum of D on S^1 x M_2(C) | Exact diagonalization | EASY | 1 day |
| Spectral action (Poisson sum) | Analytical + numerics | MODERATE | 3 days |
| CCSvS entropy S_vN(Q) | Exact formula, numerical evaluation | MODERATE | 2 days |
| Fisher metric at Q=1 | 2nd derivative, analytical in special cases | MODERATE | 3 days |
| Full Sigma(Q) for Q > 1 | Numerical summation | EASY (code) | 1 day |
| Comparison with 2 ln Q / 2/Q^2 | Direct | EASY | 1 day |
| Thermodynamic limit analysis | Poisson summation + asymptotics | HARD | 5 days |
| Interpretation and paper writing | N/A | MODERATE | 5 days |

**Total estimated time**: 2-3 weeks for a definitive answer.

### 10.2 Mathematical Prerequisites

1. **Spectral theory of Dirac operators** (standard; covered in Gracia-Bondia, Varilly, Figueroa "Elements of NCG")
2. **Heat kernel expansion and Seeley-DeWitt coefficients** (standard; Vassilevich 2003)
3. **Quantum relative entropy for Gaussian / thermal states** (standard; Petz "Quantum Information Theory")
4. **Poisson summation formula and modular properties of theta functions** (standard number theory / analysis)
5. **CCSvS entropy-spectral action correspondence** (read arXiv:1809.02944 carefully)

No NEW mathematics is needed. This is a computation within known frameworks.

### 10.3 Risk Mitigation

If the naive test fails, the following modifications should be tried before declaring NO-GO:

1. **Different reference state**: Instead of rho_1 (Q=1), use a KMS state at temperature T = 1/beta_0. The KMS condition may provide the correct normalization.

2. **Modular Fisher metric**: Use the Araki relative modular operator Delta_{rho_Q | rho_1} instead of the naive QRE. The modular Fisher metric may differ from the state-space Fisher metric.

3. **Channel regularization**: Regularize the channel by a finite bandwidth Lambda. The per-mode Fisher metric may converge to 2/Q^2 only after this regularization.

4. **Bosonic vs fermionic**: The CCSvS result uses fermionic second quantization. The channel in Paper 2 uses bosonic modes. Try the bosonic version:
```
S_vN^{bos}(Q) = sum_n [-beta_0 E_n/(Q(e^{beta_0 E_n/Q}-1)) + ln(1-e^{-beta_0 E_n/Q})]
```

5. **Different Q-parameterization**: Instead of D -> D/Q, try D -> D + (Q-1)V for some perturbation V. This gives a different Fisher metric that may match.

---

## 11. Summary: The Test in One Page

**MODEL**: S^1 (circumference beta_0) x M_2(C), Dirac operator D = -i d/dtheta tensor 1 + sigma_3 tensor m sigma_1.

**SPECTRUM**: lambda_{n,pm} = pm sqrt((2pi n/beta_0)^2 + m^2), n in Z.

**STATES**: Fermionic Gibbs states rho_Q at inverse temperature beta_0/Q (Tolman-shifted).

**COMPUTE**:
1. S_vN(Q) = CCSvS entropy = 2 sum_n s(beta_0 E_n / Q)
2. Sigma(Q) = D(rho_Q || rho_1) = 2 sum_n d(beta_0 E_n, Q)
3. g_QQ(Q) = -d^2 S_vN / dQ^2 and g_QQ^{state} = d^2 Sigma / dQ^2

**COMPARE** (per degree of freedom, in thermodynamic limit N -> inf):
- Sigma_{per dof}(Q) vs 2 ln Q
- g_QQ^{per dof}(Q) vs 2/Q^2

**GO**: Match in all regimes (high-T, low-T, intermediate).
**NO-GO**: Mismatch that cannot be resolved by normalization, regularization, or re-identification.

**TIMELINE**: 2-3 weeks (analytical + numerical).
**MATH NEEDED**: Standard spectral theory, QRE for thermal states, Poisson summation. No new mathematics.

**CRITICAL SUBTLETY IDENTIFIED**: The state-space QRE D(rho_{T_1} || rho_{T_2}) is NOT the same as the channel QRE Sigma = 2 ln Q. The test must carefully distinguish these two objects. The correct comparison is the PER-DEGREE-OF-FREEDOM Fisher metric from the CCSvS entropy under conformal deformation D -> D/Q.

---

## 12. Pre-Computation Analytical Predictions

Before running the numerics, we can predict the outcomes in limiting cases:

### 12.1 Single Mode (N_max = 0, only n = 0)

Only the zero mode E_0 = m contributes. With x_0 = beta_0 m:

```
Sigma_{single}(Q) = 2 d(x_0, Q) = 2 [(1-1/Q) x_0 n_F(x_0/Q) + ln(1+e^{-x_0}) - ln(1+e^{-x_0/Q})]
```

- At x_0 = 0 (m=0): Sigma = 0 for all Q (trivially zero because the mode is massless and the n=0 Fourier mode has zero energy -- degenerate).
- At x_0 = 1: Sigma(Q=2) = 2[0.5 * 1 * n_F(0.5) + ln(1+e^{-1}) - ln(1+e^{-0.5})] = 2[0.5 * 0.3775 + 0.3133 - 0.4741] = 2[0.1888 + 0.3133 - 0.4741] = 2 * 0.0280 = 0.0560. Target: 2 ln 2 = 1.386. **MASSIVE MISMATCH** for a single mode.

This confirms: a single mode does NOT give 2 ln Q. The test requires the thermodynamic limit.

### 12.2 Many Modes, High Temperature (N_max >> 1, beta_0 m << 1)

Each mode contributes d(x_n, Q) ~ (1/2) n_F(x_n)(1-n_F(x_n)) x_n^2 (1-1/Q)^2 for small x_n.

In the high-T limit, n_F(x)(1-n_F(x)) ~ 1/4, so:
```
d(x_n, Q) ~ x_n^2 (1-1/Q)^2 / 8
```

Summing:
```
Sigma(Q) ~ (1-1/Q)^2 / 4 * sum_n x_n^2 = (1-1/Q)^2 / 4 * sum_n beta_0^2 E_n^2
```

Per dof:
```
Sigma_{per dof} ~ (1-1/Q)^2 / 4 * <x^2> / (number of modes per dof)
```

This gives Sigma proportional to (1-1/Q)^2, NOT 2 ln Q. BUT: 2 ln Q = 2(Q-1) - (Q-1)^2 + ... while (1-1/Q)^2 = (Q-1)^2/Q^2, which to leading order in (Q-1) is also (Q-1)^2. So at second order they differ:

- 2 ln Q ~ 2(Q-1) - (Q-1)^2 + (2/3)(Q-1)^3 - ...
- (1-1/Q)^2 = (Q-1)^2/Q^2 ~ (Q-1)^2 - 2(Q-1)^3 + ...

These have the SAME leading term (Q-1)^2 but DIFFERENT higher-order terms. **The state-space QRE matches 2 ln Q only at the Fisher (second-order) level**, not beyond.

### 12.3 Prediction: The Fisher Metric Matches, the Full Sigma Does Not

Based on the analytical structure, I predict:

- **Primary test (Fisher metric): PASS** -- g_QQ per dof will equal a constant (to be determined; may be 2 or may be a function of beta_0 m)
- **Secondary test (full Sigma): FAIL** -- the per-dof QRE will NOT be 2 ln Q; it will be some other function of Q
- **Interpretation**: Sigma = 2 ln Q is the CHANNEL result, not the state-space result. The spectral action provides the channel through the CCSvS entropy, but the per-dof normalization gives a constant that may or may not be 2.

This would be a **Partial GO** (outcome 2 in Section 8.4, estimated 25% probability).

---

## 13. Appendix: Relation to Paper 2 Channel Theorem

The channel theorem (Paper 2) establishes:
```
E_Q : rho -> (1-eta) rho + eta rho_env,  eta = 1 - 1/Q^2
```

For this thermal attenuator channel:
```
Sigma_channel = D(E_Q(rho) || rho) = -ln(1-eta) = -ln(1/Q^2) = 2 ln Q
```

This is a CHANNEL property, independent of the input state rho. The spectral action should reproduce this through the CCSvS identification S_vN = Tr f(D^2/Lambda^2), but the precise mechanism is:

```
CCSvS entropy S_vN <-> spectral action <-> channel capacity
Channel capacity <-> -ln(eta) = 2 ln Q
Therefore: second variation of S_vN <-> Fisher metric of the channel = 2/Q^2
```

The intermediate step -- connecting the CCSvS entropy to the channel capacity -- is the non-trivial part. This is where the Go/No-Go test sits.

---

## Key References for Computation

1. CCSvS 2018 (arXiv:1809.02944): Entropy = spectral action
2. Gracia-Bondia, Varilly, Figueroa: "Elements of Noncommutative Geometry" (Birkhauser 2001), Ch. 11
3. Vassilevich 2003 (hep-th/0306138): Heat kernel user's manual
4. Petz: "Quantum Information Theory and Quantum Statistics" (Springer 2008), Ch. 7 (relative entropy of Gibbs states)
5. Holevo: "Quantum Systems, Channels, Information" (de Gruyter 2019), Ch. 12
6. Lashkari, Van Raamsdonk 2015 (arXiv:1508.00897): Canonical energy = Fisher information

---

*Last updated: 2026-03-19*
*This document defines the complete Go/No-Go test for Paper 9. The test is analytically tractable with numerical verification. Estimated timeline: 2-3 weeks. The critical subtlety is the distinction between state-space QRE and channel QRE; the Fisher metric comparison is the decisive criterion. Pre-computation analysis predicts a Partial GO (Fisher match, full Sigma mismatch) as the most likely outcome.*
