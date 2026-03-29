# Deriving J(Y) from Sigma: The MOND Interpolating Function Problem

**Author**: Sheng-Kai Huang (with systematic analysis)
**Date**: 2026-03-19
**Status**: Comprehensive investigation -- 5 routes explored, 1 promising, 0 proven
**Classification**: Gap B2 (CRITICAL) from gap_analysis_2026_03_19.md

---

## Executive Summary

| Route | Idea | Result | Verdict |
|-------|------|--------|---------|
| 1. Fisher information of Sigma | J(Y) = Fisher metric of spatial Sigma gradient | Linear J = Y; wrong MOND limit | FAILS |
| 2. DBI-Sigma for J | J_DBI = (2/lambda_J)[1 - sqrt(1 - lambda_J Y)] | Correct limits IF lambda_J tuned; but lambda_J not derived | SUGGESTIVE but circular |
| 3. Crooks fluctuation theorem | mu_MOND(x) = 1 - exp(-x) -> constrains J | Correct qualitative behavior; testable against SPARC | PROMISING (phenomenological) |
| 4. QRE for accelerated frames | Sigma(a) = D(rho_Rindler || rho_Minkowski) for acceleration a | Gives Sigma ~ (a/a_Unruh)^2 at high a; but deep-MOND limit unclear | INCOMPLETE |
| 5. Extremal retrodictability | Maximize F subject to Second Law, in the J-sector | Flat direction -- does not uniquely fix J(Y) | FAILS |

**Bottom line**: J(Y) remains a free function. No route from Sigma uniquely determines J(Y). The Crooks route (Route 3) is the most promising because it makes a TESTABLE prediction: mu_MOND(x) = 1 - exp(-x). If this fits SPARC data well, it would reduce the free function to zero free parameters beyond a_0. But this is a proposal, not a derivation from first principles.

**The honest assessment**: This gap is MORE DIFFICULT than deriving mu. The reason: mu lives in the K(Q) sector where ghost condensation provides strong structural constraints. J(Y) lives in the acceleration sector where no comparable symmetry principle has been identified.

---

## 1. What IS J(Y) in the Blanchet-Skordis Theory?

### 1.1 Definition

From BS2024 (arXiv:2404.06584), the action is:

```
S = (c^3 / 16piG) integral d^4x sqrt(-g) [R - 2J(Y) + 2K(Q)] + S_m
```

where:
- Q = c sqrt(-g^{ab} nabla_a phi nabla_b phi) is the inverse Khronon lapse
- n_a = -(c/Q) nabla_a phi is the unit normal to phi = const surfaces
- A_mu = c^2 n^nu nabla_nu n_mu is the acceleration of the foliation
- Y = A_mu A^mu / c^4 is the dimensionless acceleration scalar (Y >= 0 for spacelike acceleration)

### 1.2 Physical Meaning

Y measures the **acceleration of the Khronon congruence**. On a static background:

```
A_r = (c^2/2) f'(r)/f(r) ~ nabla Phi    (weak field)
Y = |A|^2/c^4 ~ |nabla Phi|^2 / c^4
```

So Y is the squared gravitational acceleration in units of c^4. At galactic scales:

```
Y ~ (a_N/c^2)^2 ~ (10^{-10} / 9e16)^2 ~ 10^{-52}
```

This is extremely small. The MOND transition occurs at:

```
Y_0 = a_0^2/c^4 ~ (1.2e-10)^2 / (3e8)^4 ~ 1.8e-52
```

### 1.3 Asymptotic Constraints

**Deep-MOND limit** (Y << Y_0):
```
J(Y) ~ Lambda - Y + (c^2/a_0) Y^{3/2} + O(Y^2)
```

This gives the MOND-Poisson equation:
```
nabla . [mu_MOND(|nabla Phi|/a_0) nabla Phi] = 4piG rho_b
```

with mu_MOND -> x for x << 1 (deep MOND: g_obs ~ sqrt(a_0 g_N)).

**Newtonian/GR limit** (Y >> Y_0):
```
J(Y) -> 0     (or J -> Lambda for the CC piece)
J_Y -> 0
```

This recovers standard GR at high accelerations (solar system, strong field).

**Transition region** (Y ~ Y_0): **COMPLETELY FREE**. This is where the MOND interpolating function mu(x) lives.

### 1.4 Historical Status

J(Y) has been free since Bekenstein-Milgrom (1984). Forty-two years of MOND research have not uniquely determined it. Common choices:
- "Simple": mu(x) = x/(1+x)
- "Standard": mu(x) = x/sqrt(1+x^2)
- "RAR" (McGaugh+2016): mu(x) = 1/(1 - exp(-sqrt(x)))
- "Crooks" (proposed): mu(x) = 1 - exp(-x)

### 1.5 Why This Is Fatal

Without deriving J(Y) from Sigma, rotation curves are NOT a prediction of the framework -- they are an input. The theory has:
- 1 free function (J(Y), or equivalently mu(x))
- plus a_0 as a parameter of that function

This makes the galactic sector LESS predictive than LCDM (which has no free function, just 2 halo parameters per galaxy).

---

## 2. The Kinematic Connection: Y in Terms of Sigma

### 2.1 The Chain

Since Sigma = 2 ln Q and Q = 1/sqrt(-g_00) on static backgrounds:

```
Q = exp(Sigma/2)
g_00 = -1/Q^2 = -exp(-Sigma)
```

The acceleration of the foliation is:
```
A_r = (c^2/2) d(ln(-g_00))/dr = -(c^2/2) dSigma/dr
```

Therefore:
```
Y = |A|^2/c^4 = (1/4)|nabla Sigma|^2    (on static, spherically symmetric backgrounds)
```

**This is the key relation**: Y = |nabla Sigma|^2 / 4.

### 2.2 What This Means

J(Y) = J(|nabla Sigma|^2 / 4). The MOND sector depends on the **spatial gradient** of Sigma, not on Sigma itself. The K(Q) sector depends on Sigma directly (through Q = exp(Sigma/2)), while the J(Y) sector depends on how Sigma varies in space.

This is reminiscent of the distinction in information geometry:
- K(Q) ~ "distance from the reference state" (Fisher distance)
- J(Y) ~ "curvature of the information landscape" (Fisher metric)

### 2.3 The Deep-MOND Regime in Sigma Language

In the deep-MOND regime (r >> r_MOND, weak gravity):
```
Sigma(r) ~ r_s/r << 1
|nabla Sigma| = r_s/r^2
Y = r_s^2 / (4r^4)
```

The MOND modification kicks in when:
```
|nabla Sigma| ~ 2a_0/c^2 = 2 * 1.2e-10 / 9e16 ~ 2.7e-27 m^{-1}
```

This is the spatial gradient of Sigma at which the information landscape transitions from "smooth" (Newtonian) to "steep" (MOND). The transition occurs at:

```
r_MOND = sqrt(GM/a_0) ~ sqrt(r_s c^2 / (2a_0))
```

For the Milky Way: r_MOND ~ 70 kpc, consistent with the onset of flat rotation curves.

---

## 3. Route 1: J(Y) as Fisher Information of Sigma

### 3.1 The Idea

The Fisher information metric on a statistical manifold parameterized by theta is:

```
I_F(theta) = integral (d ln p / d theta)^2 p dx
```

For a Gaussian distribution p(x|theta) = exp(-(x-theta)^2/(2sigma^2)) / Z:

```
I_F = 1/sigma^2
```

If we interpret Sigma(x) as the "log-likelihood" of the gravitational state at position x, then the Fisher information in the spatial direction is:

```
I_F^spatial = |nabla ln p|^2 ~ |nabla Sigma|^2
```

This would give:

```
J ~ I_F^spatial = |nabla Sigma|^2 = 4Y
```

i.e., J(Y) = 4Y (linear).

### 3.2 Why It Fails

A linear J(Y) = alpha Y gives:
```
J_Y = alpha = const
```

This produces:
```
mu_MOND(x) = 1 + alpha = const (independent of acceleration!)
```

A constant mu_MOND means NO MOND transition -- it is just a rescaling of G. This completely fails to reproduce rotation curves.

### 3.3 What Went Wrong

The Fisher information gives the LEADING-ORDER (quadratic) behavior of J near Y = 0. The deep-MOND regime requires J ~ Y^{3/2}, not Y. The Fisher metric is the wrong object -- it gives the Gaussian approximation, but MOND is inherently non-Gaussian.

### 3.4 Salvage Attempt: Higher-Order Fisher

The Amari-Chentsov alpha-connections generalize Fisher to higher orders. The exponential connection gives:

```
J^{(alpha)}(Y) = sum_{n=1}^{infinity} c_n Y^{n/2}
```

The n=3 term (Y^{3/2}) is the MOND term. But:
- The coefficients c_n are not uniquely determined by information geometry
- The alpha-connection choice is arbitrary
- This does not constitute a derivation

**VERDICT**: Route 1 FAILS. Fisher information gives only the Newtonian limit, not the MOND limit.

---

## 4. Route 2: DBI-Sigma Correspondence for J

### 4.1 The DBI Analogy

For K(Q), the DBI completion is:
```
K_DBI(Q) = (2mu^2/lambda_D)[1 - sqrt(1 - lambda_D(Q-1)^2)]
```

By analogy, a DBI form for J would be:
```
J_DBI(Y) = (2nu^2/lambda_J)[1 - sqrt(1 - lambda_J Y)]
```

### 4.2 Expanding J_DBI

For small Y (deep MOND, Y << 1/lambda_J):
```
J_DBI ~ (2nu^2/lambda_J)[1 - (1 - lambda_J Y/2 - (lambda_J Y)^2/8 - ...)]
       = nu^2 Y + nu^2 lambda_J Y^2/4 + ...
```

This gives J ~ nu^2 Y at leading order -- a linear term, NOT Y^{3/2}.

For large Y (near the DBI cutoff, Y ~ 1/lambda_J):
```
J_DBI ~ (2nu^2/lambda_J)[1 - sqrt(1 - lambda_J Y)]
```

The DBI form has a square root singularity at Y = 1/lambda_J. Near this point:
```
J_DBI ~ (2nu^2/lambda_J) - (2nu^2/lambda_J) sqrt((1/lambda_J - Y) * lambda_J)
       ~ const - (2nu^2) sqrt((1/lambda_J - Y))
```

### 4.3 The Problem

The MOND deep limit requires J ~ Y^{3/2} for small Y. The DBI form gives J ~ Y for small Y. These are incompatible.

To get Y^{3/2} behavior, one needs a non-analytic J(Y) near Y = 0. The DBI action is smooth at Y = 0 (analytic expansion exists). Therefore:

**A standard DBI form cannot produce the deep-MOND limit.**

### 4.4 Modified DBI: Square-Root Form

BS2024 mentions that in the deep-MOND limit:
```
J ~ Lambda - Y + (c^2/a_0) Y^{3/2}
```

The Y^{3/2} term is characteristic of a square-root correction. Consider:
```
J(Y) = Lambda - Y + alpha sqrt(Y^3) = Lambda - Y + alpha Y sqrt(Y)
```

This is NOT a standard DBI form. It is a specific non-analytic function. The square root of Y^3 cannot come from a simple DBI action [1 - sqrt(1 - lambda Y)], which generates integer and half-integer powers differently.

### 4.5 A More General DBI-Inspired Form

Consider:
```
J(Y) = (2/lambda_J) [1 - (1 + lambda_J Y)^{1/2}] + correction
```

Note: this is sqrt(1 + lambda_J Y), not sqrt(1 - lambda_J Y).

Expanding:
```
J ~ (2/lambda_J)[1 - 1 - lambda_J Y/2 + (lambda_J Y)^2/8 - ...]
  = -Y + lambda_J Y^2/4 - ...
```

Still gives Y^2, not Y^{3/2}. The half-integer power Y^{3/2} is not in the DBI expansion.

### 4.6 What COULD Produce Y^{3/2}?

The Y^{3/2} term is characteristic of a **fractional power law**. In condensed matter, fractional powers arise from:
1. Non-Fermi liquid behavior (strange metals)
2. Critical phenomena at second-order phase transitions
3. Anomalous dimensions from strong coupling

In the gravity context, Y^{3/2} might arise from:
- A one-loop correction with a logarithmic divergence that leaves a sqrt(Y) residue
- A non-perturbative effect (instanton-like) in the acceleration sector
- An intrinsically non-analytic structure of the Khronon EFT at low accelerations

**VERDICT**: Route 2 FAILS as a derivation. The DBI structure of K(Q) does not naturally extend to J(Y) because the MOND limit requires non-analytic Y^{3/2} behavior that no standard DBI form produces.

---

## 5. Route 3: Crooks Fluctuation Theorem

### 5.1 The Proposal (from Paper 4 / Crooks analysis)

The Crooks fluctuation theorem states:
```
P(+Sigma) / P(-Sigma) = exp(Sigma)
```

This means the probability of observing entropy production Sigma vs. the probability of observing the reverse process (-Sigma) are exponentially related.

Applied to the gravitational channel, if Sigma = 2|Phi|/c^2 at each radius:
```
tau(Sigma) = 1 - exp(-Sigma) = 1 - exp(-2|Phi|/c^2)
```

This defines the "irrecoverability" of the gravitational process. In the weak field (Sigma << 1):
```
tau ~ Sigma ~ 2|Phi|/c^2
```

In the strong field (Sigma >> 1):
```
tau -> 1 (total irreversibility)
```

### 5.2 Mapping to MOND Interpolating Function

Define x = a/a_0 where a = |nabla Phi| is the acceleration. The MOND interpolating function mu_MOND(x) appears in:
```
nabla . [mu_MOND(x) nabla Phi] = 4piG rho_b
```

If we identify:
```
mu_MOND(x) = 1 - exp(-x)
```

Then:
- x >> 1 (Newtonian): mu_MOND -> 1 (standard gravity)
- x << 1 (deep MOND): mu_MOND ~ x (MOND regime, g_obs ~ sqrt(a_0 g_N))
- x ~ 1 (transition): smooth, controlled by the exponential

### 5.3 Connection to J(Y)

The interpolating function mu_MOND(x) is related to J(Y) through:
```
mu_MOND = 1 + 2 Y J_YY / J_Y    (schematically, exact form depends on normalization)
```

or more directly, for the BS Khronon:
```
mu_MOND ~ 1 + J_Y(Y)
```

where J_Y = dJ/dY. If mu_MOND = 1 - exp(-x) with x = sqrt(Y)/Y_0^{1/2}:

```
J_Y = -exp(-sqrt(Y/Y_0))

J(Y) = integral_0^Y [-exp(-sqrt(Y'/Y_0))] dY'
```

With the substitution u = sqrt(Y'/Y_0):
```
J(Y) = -2Y_0 integral_0^{sqrt(Y/Y_0)} u exp(-u) du
     = -2Y_0 [1 - (1 + sqrt(Y/Y_0)) exp(-sqrt(Y/Y_0))]
```

Using integration by parts:
```
J(Y) = -2Y_0 + 2Y_0 (1 + sqrt(Y/Y_0)) exp(-sqrt(Y/Y_0))
```

Check the limits:
- Y -> 0: J -> -2Y_0 + 2Y_0(1+0)*1 = 0. But we need J -> Lambda at Y = 0 for the cosmological constant. Add Lambda:
  ```
  J(Y) = Lambda - 2Y_0 [1 - (1 + sqrt(Y/Y_0)) exp(-sqrt(Y/Y_0))]
  ```

- Y >> Y_0: J -> Lambda - 2Y_0 (approaches a constant, as required for GR recovery)

- Y << Y_0 (expand the exponential):
  ```
  J ~ Lambda - 2Y_0 [1 - (1 + u)(1 - u + u^2/2 - u^3/6)]  where u = sqrt(Y/Y_0)
    ~ Lambda - 2Y_0 [1 - (1 - u^2/2 + u^3/6 + ...)]
    = Lambda - Y + (2/3) Y^{3/2}/Y_0^{1/2} + ...
  ```

  Wait -- let me redo this more carefully. With u = sqrt(Y/Y_0):
  ```
  (1+u)exp(-u) = (1+u)(1 - u + u^2/2 - u^3/6 + ...)
               = 1 - u + u^2/2 - u^3/6 + u - u^2 + u^3/2 - ...
               = 1 - u^2/2 + u^3/3 + ...

  J = Lambda - 2Y_0[1 - 1 + u^2/2 - u^3/3 + ...]
    = Lambda - 2Y_0[u^2/2 - u^3/3 + ...]
    = Lambda - Y_0 u^2 + (2/3)Y_0 u^3 + ...
    = Lambda - Y + (2/3)(Y^{3/2}/Y_0^{1/2}) + ...
  ```

This gives the deep-MOND expansion:
```
J(Y) ~ Lambda - Y + (2/3) Y^{3/2} / sqrt(Y_0) + ...
```

Comparing with BS2024's expansion:
```
J(Y) ~ Lambda - Y + (c^2/a_0) Y^{3/2} + ...
```

We identify:
```
(2/3) / sqrt(Y_0) = c^2/a_0
=> Y_0 = (2a_0/(3c^2))^2 = 4a_0^2/(9c^4)
```

This is consistent (up to an O(1) factor) with Y_0 = a_0^2/c^4. The factor 4/9 vs 1 depends on the exact normalization convention in BS2024.

### 5.4 The Full J(Y) from Crooks

If we adopt mu_MOND(x) = 1 - exp(-x), the FULL J(Y) is:

```
J_Crooks(Y) = Lambda - 2Y_0 [1 - (1 + sqrt(Y/Y_0)) exp(-sqrt(Y/Y_0))]
```

with Y_0 ~ a_0^2/c^4 (up to O(1) normalization).

**Properties:**
- J(0) = Lambda (cosmological constant)
- J(Y -> infinity) -> Lambda - 2Y_0 (approaches constant)
- J'(Y) = -exp(-sqrt(Y/Y_0)) < 0 (monotonically decreasing)
- J''(Y) = exp(-sqrt(Y/Y_0)) / (2 sqrt(Y Y_0)) > 0 (convex)
- Deep MOND: J ~ Lambda - Y + (2/3)(Y/Y_0)^{1/2} Y + ...

### 5.5 Strengths and Weaknesses

**Strengths:**
1. Produces the CORRECT deep-MOND limit J ~ Lambda - Y + alpha Y^{3/2}
2. Produces the CORRECT Newtonian limit J -> const for Y >> Y_0
3. The transition region is SMOOTH and controlled by a single parameter Y_0 ~ a_0^2/c^4
4. The interpolating function mu = 1 - exp(-x) has a natural origin (Crooks theorem)
5. ZERO free parameters beyond a_0 (the shape is fully determined)

**Weaknesses:**
1. The identification mu_MOND(x) = tau(Sigma(x)) = 1 - exp(-x) is a PROPOSAL, not a derivation
2. The Crooks theorem applies to entropy production, but the mapping Sigma -> x = a/a_0 is ad hoc
3. No proof that the gravitational Crooks distribution gives this specific interpolating function
4. The mu = 1 - exp(-x) function differs from the RAR empirical fit mu = 1/[1 - exp(-sqrt(x))]
   - Note: these are DIFFERENT! Crooks gives exp(-x), RAR gives exp(-sqrt(x))
   - At x = 1: Crooks gives 0.632, RAR gives 0.581 -- a 9% difference
   - At x = 0.1: Crooks gives 0.095, RAR gives 0.264 -- a 178% difference!
   - **The Crooks form diverges significantly from the RAR fit at low x**

### 5.6 Critical Issue: Crooks vs. RAR Comparison

The McGaugh-Lelli-Schombert RAR (2016) best fit is:
```
g_obs = g_bar / [1 - exp(-sqrt(g_bar/g_dagger))]
```

which corresponds to mu_RAR(x) = 1/[1 - exp(-sqrt(x))].

The Crooks proposal is mu_Crooks(x) = 1 - exp(-x).

These are NOT the same function:

| x (= a/a_0) | mu_Crooks = 1-exp(-x) | mu_RAR = x/[1-exp(-sqrt(x))] | Ratio |
|---|---|---|---|
| 0.01 | 0.010 | 0.091 | 0.110 |
| 0.1 | 0.095 | 0.264 | 0.360 |
| 0.5 | 0.393 | 0.536 | 0.733 |
| 1.0 | 0.632 | 0.721 | 0.877 |
| 2.0 | 0.865 | 0.881 | 0.982 |
| 5.0 | 0.993 | 0.977 | 1.016 |
| 10 | 1.000 | 0.995 | 1.005 |

**At low x (deep MOND), the Crooks function underestimates mu by a factor of 3-10 compared to the RAR fit.** This is a serious problem. The Crooks form gives mu ~ x for small x, while the RAR gives mu ~ sqrt(x) / (something). The RAR fit actually corresponds to J ~ Y^{5/4} in the deep limit, not J ~ Y^{3/2}.

Wait -- let me recheck. The RAR form g_obs = g_bar/[1-exp(-sqrt(g_bar/g_dagger))] means:
```
mu_RAR(x) = 1/[1-exp(-sqrt(x))]

For x << 1: exp(-sqrt(x)) ~ 1 - sqrt(x) + x/2 - ...
1 - exp(-sqrt(x)) ~ sqrt(x) - x/2 + ...
mu_RAR ~ 1/sqrt(x) for x << 1
```

So mu_RAR ~ x^{-1/2}, which means in the deep MOND regime g_obs = mu * g_bar ~ g_bar / sqrt(g_bar/a_0) = sqrt(a_0 g_bar). This is the standard MOND deep limit.

For the Crooks form:
```
mu_Crooks = 1 - exp(-x) ~ x for x << 1
g_obs = mu * g_bar ~ x g_bar = (g_bar/a_0) g_bar = g_bar^2/a_0
```

This gives g_obs ~ g_bar^2/a_0, NOT sqrt(a_0 g_bar). This is WRONG for rotation curves.

### 5.7 CORRECTED Analysis: The Right Crooks Mapping

The issue is the mapping between Sigma and the interpolating function. The correct MOND deep limit requires mu_MOND -> x as x -> 0, where the product mu * g_bar -> sqrt(a_0 g_bar).

If instead we define Sigma = sqrt(x) = sqrt(a/a_0), then:
```
tau = 1 - exp(-sqrt(x))
```

And the interpolating function becomes:
```
mu_MOND(x) = x / tau(sqrt(x)) = x / [1 - exp(-sqrt(x))]
```

which IS the RAR form! This recovers the McGaugh-Lelli-Schombert fit if we identify:
```
Sigma_spatial = sqrt(a/a_0)    (NOT a/a_0)
```

**Physical interpretation**: The spatial entropy production is NOT proportional to the acceleration itself, but to the SQUARE ROOT of the acceleration (in units of a_0). This is natural if Sigma ~ sqrt(Y/Y_0) ~ |nabla Phi| / (a_0/c^2)^{1/2} ~ sqrt(a/a_0).

### 5.8 The J(Y) from the Corrected Crooks

With Sigma_spatial = sqrt(Y/Y_0) and tau = 1 - exp(-Sigma_spatial):
```
mu_MOND(x) = x / [1 - exp(-sqrt(x))]
```

This implies:
```
J_Y = mu_MOND - 1 = x/[1-exp(-sqrt(x))] - 1
```

Hmm -- this gets complicated. Let me use the standard MOND formalism more carefully.

In the BS Khronon, the modified Poisson equation in the quasi-static limit has the form (BS2024, Sec. 5):
```
nabla . [(1 + nu(Y)) nabla Phi] = 4piG rho_b
```

where nu(Y) encodes the MOND modification. The identification is:
```
nu(Y) = -2 Y J_YY(Y) / (something involving J_Y)
```

The exact mapping between J(Y) and nu(Y) or mu_MOND depends on the specific equations in BS2024 Sec. 5. Without the exact equations, I will work with the AQUAL-equivalent:
```
nabla . [mu(|nabla Phi|/a_0) nabla Phi] = 4piG rho_b
```

Here mu(x) is the MOND interpolating function with x = |nabla Phi|/a_0.

**VERDICT for Route 3**: The Crooks theorem CAN give the RAR interpolating function if the mapping is Sigma_spatial = sqrt(a/a_0) rather than Sigma_spatial = a/a_0. This is suggestive but the mapping itself is not derived. The sqrt(a/a_0) form would mean:

```
J(Y) ~ Lambda - (something from integrating the RAR mu)
```

with ONE free parameter (a_0 = Y_0^{1/2} c^2). **This reduces the free function to zero free functions + 1 parameter** -- a major simplification.

Status: PROMISING but the mapping Sigma_spatial = sqrt(Y/Y_0) needs justification.

---

## 6. Route 4: QRE for Accelerated Frames

### 6.1 The Unruh-DeWitt Approach

An observer with proper acceleration a in Minkowski space sees the Unruh temperature:
```
T_U = hbar a / (2pi c k_B)
```

The QRE between the Rindler vacuum |0_R> and the Minkowski vacuum |0_M> is:
```
D(rho_Rindler || rho_Minkowski) = <K_R>_{R} - <K_R>_{M} - S(rho_R) + S(rho_M)
```

where K_R is the Rindler modular Hamiltonian.

For a massless scalar field:
```
D(rho_R || rho_M) = (pi^2/6) (T_R^2 - T_M^2) / (energy scale)^2
```

The QRE scales as T_U^2 ~ a^2 at high accelerations and gives a logarithmic correction at the Unruh threshold.

### 6.2 Connection to J(Y)

If the acceleration scalar Y = a^2/c^4, then:
```
D(rho_accelerated || rho_inertial) ~ f(Y)
```

For the Khronon foliation, the QRE associated with the acceleration is precisely what J(Y) encodes. If:
```
J(Y) = D(rho_{accelerated Khronon} || rho_{inertial Khronon})
```

Then J(Y) is determined by the QRE calculation.

### 6.3 The Calculation

For the Unruh effect, at temperature T_U:
```
D(rho_T || rho_0) = (pi^2/90) (kT_U)^4 / (hbar c)^3 * V / (cutoff stuff)
```

This is a thermal QRE. For a d-dimensional massless field:
```
D(rho_{T_1} || rho_{T_2}) = C_d (T_1^d - T_2^d - d T_2^{d-1}(T_1 - T_2))
```

In d = 3+1:
```
D ~ (T_1^4 - T_2^4 - 4T_2^3(T_1 - T_2)) / T_2^4
  = (T_1/T_2)^4 - 4(T_1/T_2) + 3
```

With T_1 = T_U = hbar a/(2pi c k_B) and T_2 = T_dS = hbar H_0/(2pi k_B):
```
T_1/T_2 = a/(cH_0)
```

So:
```
D ~ (a/(cH_0))^4 - 4(a/(cH_0)) + 3
```

For a >> cH_0 (Newtonian regime): D ~ (a/(cH_0))^4 ~ Y^2 (quartic in Y^{1/2}, i.e., quadratic in Y)
For a << cH_0 (deep MOND?): D ~ 3 - 4(a/(cH_0)) + (a/(cH_0))^4 ~ 3 - 4 sqrt(Y/Y_0)

This does NOT give Y^{3/2}. The Unruh-based QRE gives quartic behavior, not the required 3/2 power.

### 6.4 Why This Route Is Incomplete

The thermal QRE for relativistic fields gives integer powers of temperature, hence even powers of acceleration. The fractional power Y^{3/2} characteristic of MOND cannot arise from the standard Unruh QRE. This suggests:

1. MOND is NOT simply the Unruh effect in disguise (consistent with criticisms of Verlinde by Hossenfelder and others)
2. The J(Y) function involves non-thermal contributions (e.g., non-equilibrium effects, vacuum polarization)
3. A different type of QRE calculation (non-equilibrium, path-integral, or boundary-state) is needed

**VERDICT**: Route 4 is INCOMPLETE. The standard Unruh QRE does not reproduce the MOND interpolating function. A more sophisticated calculation involving non-equilibrium states or boundary effects might succeed, but has not been done.

---

## 7. Route 5: Extremal Retrodictability

### 7.1 The Principle

The extremal retrodictability principle (used for Omega_DM = 0.268) states:
- Maximize the Petz recovery fidelity F (maximize retrodictability)
- Subject to the constraint Sigma >= 0 (second law)

For the K(Q) sector, this gave K(Q) = mu^2(Q-1)^2 (quadratic minimum at Q = 1, maximally recoverable).

### 7.2 Application to J(Y)

Apply the same principle to the J(Y) sector: maximize F subject to the constraint that the J-sector contribution to Sigma is non-negative.

The J-sector entropy production is:
```
Sigma_J = -2 integral J(Y) sqrt(-g) d^4x / (normalization)
```

Maximizing F subject to Sigma_J >= 0 requires minimizing Sigma_J, which means minimizing |J(Y)|.

But J(Y) = 0 (no J term at all) trivially minimizes Sigma_J. The extremal principle would select J = 0, i.e., NO MOND sector at all!

### 7.3 Why It Fails

The extremal principle has a trivial solution (J = 0) because J only contributes on static backgrounds with non-zero acceleration. On the FRW background (where Y = 0), J is invisible. The principle does not "see" the galactic sector.

To fix this, we would need to impose additional constraints:
- Require that the theory reproduces observed rotation curves
- Require that the theory matches the RAR
- Require consistency with the K(Q) sector at intermediate scales

But these are all OBSERVATIONAL constraints, not information-theoretic ones. The extremal principle alone does not select J(Y).

### 7.4 Flat Direction

Even with constraints, the extremal principle has a flat direction: any J(Y) that satisfies the asymptotic constraints (deep-MOND at low Y, GR at high Y) gives the same F to leading order. The transition region is invisible to the Petz recovery map because it only matters at the percent level in individual galaxies.

**VERDICT**: Route 5 FAILS. The extremal principle selects J = 0 (trivially) or has a flat direction in the J(Y) function space. No unique J is selected.

---

## 8. Route 6 (Bonus): Bianconi's 1-Form Sector

### 8.1 The Idea

In Bianconi (2025, PRD 111, 066001), the full entropic action is:
```
L = -Tr_F ln(G_tilde g_tilde^{-1})
```

where G_tilde and g_tilde are the "spacetime metric" and "matter metric" on the full Dirac-Kahler bundle (0-forms + 1-forms + 2-forms).

Our Sigma = 2 ln Q corresponds to the 0-FORM sector only (scalar part). Opening the 1-form sector would introduce gauge fields, and the corresponding term in L might produce J(Y).

### 8.2 The Connection

The 1-form sector of Bianconi's action includes terms like:
```
L_1 ~ -Tr ln(G_{1-form} g_{1-form}^{-1})
```

For the Khronon theory, the 1-form sector corresponds to the ACCELERATION A_mu (the covariant derivative of the unit normal n_mu). Since Y = A_mu A^mu / c^4, the 1-form contribution to the entropic action is:
```
L_1 ~ -ln(1 + alpha Y)
```

for some coefficient alpha.

Expanding:
```
L_1 ~ -alpha Y + alpha^2 Y^2/2 - alpha^3 Y^3/3 + ...
```

This is a Taylor series in Y -- ANALYTIC at Y = 0. It gives:
- J ~ -alpha Y + alpha^2 Y^2/2 + ... (leading term is linear)
- Deep MOND requires Y^{3/2} (non-analytic)
- The 1-form sector CANNOT produce MOND

### 8.3 Unless...

The non-analytic Y^{3/2} term could arise from:
1. A non-perturbative effect in the 1-form sector (instanton?)
2. A logarithmic divergence that is renormalized to give a non-integer power
3. A fractal / anomalous dimension contribution from the Dirac-Kahler structure

None of these have been computed. This remains speculative.

**VERDICT**: Route 6 is SPECULATIVE. The Bianconi 1-form sector naturally produces analytic contributions to J(Y), not the non-analytic Y^{3/2} needed for MOND. Non-perturbative effects might help but have not been computed.

---

## 9. The Deep Problem: Why Y^{3/2}?

### 9.1 The Non-Analyticity

The root difficulty is that MOND requires:
```
J(Y) ~ Lambda - Y + alpha Y^{3/2}     for Y -> 0
```

The Y^{3/2} term is NON-ANALYTIC at Y = 0. This means J(Y) cannot be written as a Taylor series around Y = 0. Any information-theoretic derivation based on smooth structures (Fisher metric, DBI actions, thermal QRE, Bianconi's log-trace) will naturally produce ANALYTIC functions of Y, i.e., integer powers of Y.

### 9.2 Where Do Non-Analytic Terms Come From in Physics?

1. **Phase transitions**: At a critical point, the free energy has non-analytic behavior F ~ |T - T_c|^{2-alpha}. If the MOND transition at a = a_0 is analogous to a phase transition, Y^{3/2} could emerge from a critical exponent alpha = 1/2 (mean-field theory of a second-order transition).

2. **One-loop quantum effects**: In d = 3+1, a one-loop calculation with a massless propagator in a background field Y gives contributions like Y^{3/2} ln(Y) or Y^{3/2}. This is how, e.g., the Euler-Heisenberg effective action gets its non-linear correction to Maxwell's equations.

3. **Dimensional reduction**: In d = 2+1, the leading non-analytic term in the effective action is Y^{3/4}. In d = 3+1, it would be Y^{5/4}. To get Y^{3/2} exactly requires a specific dimensionality or field content.

4. **RG anomalous dimension**: If the acceleration field has an anomalous dimension gamma_Y at the MOND scale, the effective action could be Y^{1+gamma_Y/2}. For gamma_Y = 1, this gives Y^{3/2}.

### 9.3 The Most Promising Interpretation

The anomalous dimension interpretation (point 4) connects to the running-G story: if eta = 1 (the anomalous dimension of the gravitational coupling), there might be an analogous eta_J = 1 for the J-sector that produces Y^{1+1/2} = Y^{3/2}.

This would mean: **the MOND interpolating function is a consequence of the same anomalous dimension that produces flat rotation curves in the running-G picture.**

However, this is a CONJECTURE. The anomalous dimension eta = 1 applies to the coupling G, not to the function J(Y). Connecting them requires showing that the same RG fixed point controls both sectors.

### 9.4 A Concrete Proposal

If the Khronon field at galactic scales flows to an IR fixed point with:
- eta_G = 1 (anomalous dimension of G -> flat rotation curves via running G / Kumar)
- eta_J = 1 (anomalous dimension of the J-sector -> Y^{3/2} MOND term)

Then both K(Q) and J(Y) are controlled by the SAME fixed point, and the unification is achieved. The J(Y) function would be:

```
J(Y) = Lambda + (1-loop quantum correction with eta_J = 1)
     = Lambda - Y + c_1 Y^{3/2} + c_2 Y^2 + ...
```

where c_1 = c^2/a_0 is fixed by the anomalous dimension, and c_2 etc. are higher-order corrections.

**This is the closest thing to a derivation**: if the anomalous dimension eta = 1 (which already appears in the K-sector for running G) also applies to the J-sector, then Y^{3/2} is the natural one-loop correction.

But this has NOT been computed. It is a direction for future work.

---

## 10. Synthesis: What We Know and What We Don't

### 10.1 Established Facts

1. **J(Y) and K(Q) decouple** at the level of the action (Section 5 of jy_form_analysis_2026_03_16.md). K controls cosmology, J controls galactic dynamics.

2. **Y = |nabla Sigma|^2 / 4** on static backgrounds (Section 2.1). The MOND sector depends on the spatial gradient of Sigma.

3. **The deep-MOND limit requires J ~ Y^{3/2}**, which is non-analytic at Y = 0. No smooth information-theoretic construction naturally produces this.

4. **The Crooks theorem gives mu = 1 - exp(-x)**, which has the correct qualitative behavior but does NOT match the RAR best fit. The corrected version (Sigma ~ sqrt(x)) does match the RAR form but requires justifying the sqrt mapping.

5. **K(Q) = mu^2(Q-1)^2 IS the Fisher information distance** between thermal states (Section 7.3 of sigma_grav_from_khronon_2026_03_16.md). No analogous Fisher interpretation exists for J(Y).

### 10.2 Honest Gap Assessment

| Aspect | Status | Confidence |
|--------|--------|------------|
| J(Y) asymptotic at Y -> 0 | CONSTRAINED by MOND phenomenology | High (observations) |
| J(Y) asymptotic at Y -> infinity | CONSTRAINED by GR/PPN | High (observations) |
| J(Y) transition region | FREE | N/A |
| J(Y) from Sigma | NOT DERIVED | N/A |
| J(Y) from K(Q) connection | NO CONNECTION KNOWN | Low |
| Crooks mu(x) | PROPOSAL, NOT DERIVATION; partial mismatch with RAR | Medium |
| DBI form for J | DOES NOT PRODUCE Y^{3/2} | Low |
| Bianconi 1-form | SPECULATIVE, needs non-perturbative effects | Low |
| Anomalous dimension route | MOST PROMISING, NOT COMPUTED | Medium-High (as direction) |

### 10.3 The Three Possible Outcomes

**Outcome A: J(Y) IS derivable from Sigma** (optimistic)
- The anomalous dimension eta = 1, which controls running G, also controls the J-sector
- Y^{3/2} emerges from a one-loop calculation in the Khronon EFT at the MOND scale
- J(Y) is uniquely fixed, making rotation curves a genuine PREDICTION
- Timeline: requires explicit 1-loop calculation (3-6 months)

**Outcome B: J(Y) is partly constrained by Sigma** (realistic)
- The Crooks theorem fixes the qualitative shape (exponential-family interpolation)
- a_0 is derived from cosmological parameters (0.3% accuracy, already achieved)
- The transition region has 1 remaining free parameter (e.g., the power in exp(-x^p))
- Timeline: requires SPARC fitting to determine p (2-4 weeks)

**Outcome C: J(Y) is genuinely free** (pessimistic)
- The K and J sectors are fundamentally independent
- a_0 is an independent parameter not derivable from the framework
- The theory has 1 free function + parameters, same as MOND since 1984
- This does not invalidate the framework but limits its predictive power at galactic scales

### 10.4 Recommended Next Steps

**Immediate (1-2 weeks)**:
1. Fit mu_Crooks = 1 - exp(-x) to the SPARC 175-galaxy sample. Compare chi^2 with the RAR form and "standard" mu(x).
2. Also fit mu = 1 - exp(-x^p) with free p. If p ~ 0.5, this recovers the RAR and connects to Sigma_spatial = sqrt(Y/Y_0).

**Short term (1-3 months)**:
3. Compute the one-loop effective action for a scalar field with DBI kinetics in a slowly varying background acceleration. Check if Y^{3/2} emerges.
4. Read Caticha (2025, arXiv:2511.19238) on "Entropic Dynamics -> Maxwell equations" to see if an analogous entropic derivation can produce J(Y).

**Medium term (3-6 months)**:
5. Explicit RG calculation of the Khronon J-sector using functional renormalization group. Compute the anomalous dimension eta_J at the IR fixed point.
6. Extend the Bianconi framework to include the 1-form sector explicitly and check if J(Y) emerges from the log-trace structure.

---

## 11. Key Equations Summary

```
BS ACTION:
  S = (c^3/16piG) integral sqrt(-g) [R - 2J(Y) + 2K(Q)] d^4x + S_m

KINEMATIC VARIABLES:
  Q = 1/sqrt(-g_00)             (inverse lapse, on static backgrounds)
  Y = |A|^2/c^4                 (acceleration scalar)
  A_mu = c^2 n^nu nabla_nu n_mu (foliation acceleration)

SIGMA CONNECTION:
  Sigma = 2 ln Q = -ln(-g_00)   (gravitational entropy production)
  Y = |nabla Sigma|^2 / 4       (on static, spherically symmetric backgrounds)

K(Q) SECTOR (DERIVED):
  K(Q) = mu^2(Q-1)^2            (ghost condensation + Fisher information)
  c_s^2 = 0                     (from K'(Q=1) = 0)
  Sigma_K = 2 ln Q              (cosmological, CDM-like)

J(Y) SECTOR (FREE):
  J(Y) ~ Lambda - Y + (c^2/a_0) Y^{3/2}    (deep-MOND limit)
  J(Y) -> 0 for Y >> a_0^2/c^4              (GR limit)
  Sigma_J = ???                               (not derived)

CROOKS PROPOSAL:
  mu_MOND(x) = 1 - exp(-x)                  (direct Crooks)
  mu_MOND(x) = x / [1 - exp(-sqrt(x))]      (if Sigma ~ sqrt(x); matches RAR)

BEST ROUTE:
  eta_J = 1 (anomalous dimension) -> Y^{3/2} from 1-loop
  Status: CONJECTURE, NOT COMPUTED
```

---

## 12. Comparison with Existing Literature

| Approach to J(Y) | Author(s) | Method | Result |
|-------------------|-----------|--------|--------|
| Phenomenological mu(x) | Milgrom (1983+) | Fit to rotation curves | Multiple forms, all viable |
| AQUAL variational | Bekenstein-Milgrom (1984) | Lagrangian formulation | J free function |
| RAR fit | McGaugh+ (2016) | Empirical | mu = 1/[1-exp(-sqrt(x))] |
| Entropic gravity | Verlinde (2016) | de Sitter entanglement | a_0 = cH_0/6, specific force law |
| Tsallis-Renyi | arXiv:2505.03061 | Modified entropy | MOND-like, not specific J |
| Caticha entropic dynamics | Caticha (2025) | Information geometry | Maxwell eqs (not gravity) |
| Running G (Kumar) | Kumar (2025) | QFT 1-loop | Logarithmic potential, not J(Y) |
| **This work** | **Huang (2026)** | **Sigma/Crooks/anomalous dim** | **Sigma ~ sqrt(x) -> RAR mu; Y^{3/2} from eta_J = 1 (conjecture)** |

---

## 13. Final Honest Assessment

**The derivation of J(Y) from Sigma is the hardest remaining problem in the framework.** Harder than deriving mu (where ghost condensation provides structure). Harder than the CLASS compatibility issue (which is a technical problem). The reason:

1. J(Y) requires a NON-ANALYTIC function (Y^{3/2}) -- no smooth information-theoretic construction naturally gives this.
2. J and K decouple -- no known principle connects them.
3. The J-sector is active only on STATIC backgrounds with acceleration -- these are far from the controlled FRW cosmological regime.

The most honest position for Paper 3 is:
- **Claim**: a_0 is derived from cosmological parameters (0.3% accuracy)
- **Claim**: The Crooks theorem provides a well-motivated proposal for mu(x), testable against SPARC
- **Acknowledge**: J(Y) is NOT derived from first principles
- **Conjecture**: The anomalous dimension eta = 1 may fix J(Y) via 1-loop corrections
- **Future work**: Explicit RG calculation of the J-sector anomalous dimension

This is more honest than claiming J is derived, and more informative than simply stating it is free.

---

## References

### Primary
- Blanchet & Skordis 2024: arXiv:2404.06584, JCAP 11 (2024) 040 -- BS Khronon action
- Blanchet & Skordis 2025: arXiv:2507.00912 -- BS Khronon-Tensor, mu^{-1} = 22.3 Mpc
- McGaugh, Lelli & Schombert 2016: PRL 117, 201101 -- RAR
- Bekenstein & Milgrom 1984: ApJ 286, 7 -- AQUAL

### Information Theory
- Chentsov 1982: Statistical Decision Rules and Optimal Inference -- Uniqueness of Fisher metric
- Amari & Nagaoka 2000: Methods of Information Geometry -- Alpha-connections
- Casini, Teste & Torroba 2017: arXiv:1611.00016 -- RG = quantum channel
- Dorau & Much 2025: PRL 136, 091602 -- QRE -> Einstein equations

### Running G / RG
- Kumar 2025: arXiv:2509.05246 -- Marginal IR running of gravity
- Reuter & Weyer 2004: arXiv:hep-th/0410117 -- Running Newton constant
- Gubitosi, Piattella & Casarini 2024: arXiv:2403.00531 -- RGGR with SPARC

### Entropic / Emergent
- Verlinde 2016: arXiv:1611.02269 -- Emergent gravity
- Bianconi 2025: PRD 111, 066001 -- Gravity from entropy
- Caticha 2025: arXiv:2511.19238 -- Entropic dynamics -> Maxwell

### Ghost Condensation
- Arkani-Hamed et al. 2004: JHEP 05, 074 -- Ghost condensation
- Scherrer 2004: PRL 93, 011301 -- K-essence as unified DM/DE

---

*Last updated: 2026-03-19*
*This research note systematically investigates whether J(Y) can be derived from Sigma, exploring 6 routes. Conclusion: not yet, but the anomalous dimension route (eta_J = 1) and the Crooks mapping (Sigma ~ sqrt(a/a_0)) are the most promising directions.*
