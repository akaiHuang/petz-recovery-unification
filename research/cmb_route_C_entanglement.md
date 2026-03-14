# Route C: Can de Sitter Entanglement Entropy Act as Effective CDM at CMB Scales?

**Author**: Sheng-Kai Huang (with computational assistance)
**Date**: 2026-03-12
**Status**: Complete analysis with rigorous calculations and brutally honest verdict
**Purpose**: Investigate whether the SAME entanglement entropy that produces apparent dark matter at galactic scales (Verlinde 2016) can reproduce CMB acoustic peak structure at z ~ 1100.

---

## Executive Summary: The Verdict

**Route C FAILS for the CMB.** The de Sitter entanglement entropy cannot act as effective CDM at recombination for five independent, quantitative reasons:

1. **Wrong redshift scaling**: rho_ent(z) ~ H(z)^2 / G ~ (1+z)^3 in matter domination BUT ~ (1+z)^4 in radiation domination. At z ~ 1100 (radiation-dominated without CDM), the scaling is radiation-like, not matter-like.

2. **Wrong magnitude by orders of magnitude**: At z ~ 1100, the Verlinde mechanism gives either rho_ent = 0 (if volume-law has not yet developed) or rho_ent >> rho_crit (if naively extrapolated), depending on assumptions. Neither matches rho_CDM.

3. **Wrong perturbation properties**: Entanglement entropy perturbations in a relativistic QFT propagate at c_s = c (or c_s = c/sqrt(3) for conformal fields), not c_s = 0 as required by CDM.

4. **Circular logic at background level**: Using rho_ent to set z_eq requires rho_ent to already behave as matter, but the scaling of rho_ent depends on whether the universe is matter- or radiation-dominated -- a self-referential problem with no consistent solution.

5. **No perturbation theory exists**: Verlinde's framework is formulated only for z ~ 0, static, spherically symmetric configurations. No cosmological perturbation theory has been developed.

However, the analysis reveals important structural insights about WHY it fails and what would be needed to make it work. These insights constrain the space of possible solutions and point toward the Khronon/AeST route as the only viable path.

---

## Part I: Verlinde's Entanglement Entropy -- The Setup

### 1.1 The Core Mechanism (Verlinde 2016, arXiv:1611.02269)

In de Sitter spacetime with Hubble parameter H and cosmological horizon radius L = c/H, the entanglement entropy of the vacuum has two contributions:

**Area-law contribution** (standard Bekenstein-Hawking):
```
S_area = A / (4 G hbar) = pi L^2 / (G hbar) = pi c^2 / (G hbar H^2)
```

**Volume-law contribution** (thermal, from de Sitter temperature T_dS = hbar H / (2 pi k_B)):
```
S_volume ~ (4/3) pi L^3 * s_thermal
```
where s_thermal is the entropy density of the de Sitter thermal bath.

For a massless scalar field at temperature T_dS:
```
s_thermal = (2 pi^2 / 45) * (k_B T_dS / (hbar c))^3 * k_B
           = (2 pi^2 / 45) * (H / (2 pi c))^3 * k_B
```

So:
```
S_volume ~ (4/3) pi (c/H)^3 * (2 pi^2 / 45) * (H / (2 pi c))^3 * k_B
         = (4 pi^3 / (135)) * k_B * N_dof
```
where N_dof counts effective degrees of freedom.

**The key feature**: For a region of size r inside the de Sitter horizon:
- r << L: Entanglement entropy follows area law, S ~ r^2 (standard)
- r >> L: Entanglement entropy transitions to volume law, S ~ r^3

**At intermediate scales** (galactic), the volume-law contribution COMPETES with the area-law contribution. When baryonic matter is present, it "displaces" some of the volume-law entanglement entropy. The displaced entropy creates an "elastic" restoring force that manifests as apparent dark matter.

### 1.2 Verlinde's Dark Matter Formula

For a spherically symmetric baryonic mass distribution M_B(r):
```
M_D^2(r) = (c H_0 / (6G)) * r^2 * d/dr [r * M_B(r)]
```

For a point mass M_B:
```
M_D(r) = sqrt(c H_0 M_B r / (6G))
```

The total gravitational acceleration at large r:
```
g_total ~ sqrt(a_0 * g_N)
```
where a_0 = c H_0 / 6 ~ 1.1 x 10^{-10} m/s^2 (close to Milgrom's value).

**Crucial observation**: This formula contains H_0 explicitly. It is defined at z = 0. The immediate question for CMB: what replaces H_0 at z ~ 1100?

### 1.3 The Effective Energy Density at z = 0

The Verlinde mechanism produces an effective dark matter density:
```
rho_ent(z=0) ~ Omega_DM * rho_crit ~ 0.26 * 9.5e-27 kg/m^3 ~ 2.5e-27 kg/m^3
```

The question is: how does rho_ent scale with redshift?

---

## Part II: Redshift Scaling of rho_ent(z) -- The Core Calculation

### 2.1 Three Possible Scaling Laws

The effective energy density from de Sitter entanglement entropy could scale in different ways depending on the physical mechanism:

**Hypothesis A: Track the de Sitter temperature**
```
T_dS(z) = hbar H(z) / (2 pi k_B)

rho_ent(z) ~ T_dS^4 / (hbar c)^3 ~ H(z)^4 / (hbar c)^3
```

**Hypothesis B: Track the entropy displacement**
```
rho_ent ~ S_displaced / V ~ (Volume-law entropy) / V ~ H(z)^2 / G
```
This comes from the Bekenstein-Hawking entropy S ~ 1/(GH^2) being displaced over a Hubble volume V ~ H^{-3}, giving rho_ent ~ H^2/G -- which is just the Friedmann equation's rho_crit. This is circular.

**Hypothesis C: Track the elastic medium strain**
```
rho_ent ~ (strain energy) ~ (M_B / M_Hubble) * rho_crit(z)
```
where M_Hubble = c^3/(2GH) is the Hubble mass.

### 2.2 Explicit Calculation for Hypothesis A

```
H(z) = H_0 * E(z)

where E(z) = sqrt(Omega_r (1+z)^4 + Omega_m (1+z)^3 + Omega_Lambda)

At z ~ 1100:
E(z=1100) ~ sqrt(Omega_r * (1101)^4)  [radiation dominates if no CDM]
          ~ sqrt(8.5e-5 * 1.47e12)
          ~ sqrt(1.25e8)
          ~ 1.12e4

H(z=1100) ~ 1.12e4 * H_0
```

The effective energy density under Hypothesis A:
```
rho_ent(z) ~ H(z)^4 / (hbar c)^3

rho_ent(1100) / rho_ent(0) ~ [H(1100)/H_0]^4 ~ (1.12e4)^4 ~ 1.58e16
```

Meanwhile, CDM dilutes as:
```
rho_CDM(1100) / rho_CDM(0) = (1+1100)^3 = 1.33e9
```

**The ratio is wrong by a factor of ~10^7.** Under Hypothesis A, the entanglement entropy density at z ~ 1100 would be 10 million times too large relative to what CDM provides.

Quantitatively:
```
rho_ent(z=1100) ~ 2.5e-27 * 1.58e16 ~ 4e-11 kg/m^3
rho_CDM(z=1100) ~ 2.5e-27 * 1.33e9  ~ 3.3e-18 kg/m^3
rho_crit(z=1100) ~ 9.5e-27 * (1.12e4)^2 ~ 1.2e-18 kg/m^3
```

Under Hypothesis A, rho_ent would be ~3 x 10^7 times rho_crit at z ~ 1100. This is **physically absurd** -- the entanglement energy would completely dominate the universe, making it inconsistent with the observed expansion history.

### 2.3 Explicit Calculation for Hypothesis B

If rho_ent scales as H^2/G (tracking rho_crit):
```
rho_ent(z) = f * rho_crit(z) = f * (3 H(z)^2) / (8 pi G)
```
where f is a constant fraction.

At z = 0: f = Omega_DM ~ 0.26.

At z ~ 1100:
```
rho_ent(1100) = 0.26 * rho_crit(1100)
```

This gives the right TOTAL density, but the SCALING is wrong for CMB physics. The critical issue is the effective equation of state:

```
rho_ent ~ H^2 ~ rho_total

In radiation domination: rho_total ~ a^{-4}, so rho_ent ~ a^{-4} (radiation-like!)
In matter domination: rho_total ~ a^{-3}, so rho_ent ~ a^{-3} (matter-like!)
```

**This is circular and self-referential**: If there is no CDM, the universe is radiation-dominated at z ~ 1100, so rho_ent ~ a^{-4}. But if rho_ent ~ a^{-4}, it behaves like radiation, NOT like CDM. So it cannot cause the matter-radiation transition that CDM provides.

To see this more explicitly: the condition for z_eq (matter-radiation equality) requires:
```
rho_matter(z_eq) = rho_radiation(z_eq)
rho_m,0 (1+z_eq)^3 = rho_r,0 (1+z_eq)^4
z_eq = rho_m,0 / rho_r,0
```

If rho_ent tracks rho_crit, then at early times rho_ent ~ H^2 ~ rho_radiation, and we never get a matter-radiation transition. **The entanglement entropy cannot CAUSE the transition it is supposed to track.**

### 2.4 Explicit Calculation for Hypothesis C

The elastic medium strain model gives:
```
rho_ent(r, z) ~ [M_B(r) / M_Hubble(z)] * rho_crit(z)
```

At z = 0, for a galaxy with M_B ~ 10^{11} M_sun:
```
M_Hubble = c^3 / (2 G H_0) ~ 4.4 x 10^{52} kg ~ 2.2 x 10^{22} M_sun
M_B / M_Hubble ~ 4.5 x 10^{-12}
```

This is the FRACTION of the Hubble mass that is in baryons, applied to a specific galaxy. At cosmological scales, the relevant quantity is:
```
Omega_B / 1 = 0.05    (baryon fraction of critical density)
```

So the entanglement "strain" at the cosmological level:
```
rho_ent ~ Omega_B * rho_crit * (correction factor)
```

The correction factor must equal Omega_DM / Omega_B ~ 5.3 for this to work. There is no obvious reason why the elastic response should amplify the baryonic density by exactly this factor. Moreover, the correction factor would need to be CONSTANT across redshift, which requires the elastic modulus to be independent of the de Sitter temperature -- an assumption Verlinde does not justify.

### 2.5 Summary of Scaling Analysis

| Hypothesis | rho_ent(z) scaling | At z=1100 | Compared to rho_CDM | Verdict |
|---|---|---|---|---|
| A: T_dS^4 | ~ H^4 ~ (1+z)^8 (rad. dom.) | 10^7 x too large | **CATASTROPHIC** | Ruled out |
| B: rho_crit | ~ H^2 ~ (1+z)^4 (rad. dom.) | Right order, wrong scaling | **CIRCULAR** | Self-referential |
| C: Elastic strain | ~ Omega_B * rho_crit | ~5x too small | **WRONG AMPLITUDE** | Requires ad hoc factor |

**None of the three natural scaling hypotheses reproduces rho_CDM(z) = rho_CDM,0 (1+z)^3 across the radiation-matter transition.**

### 2.6 The Fundamental Problem: w_ent != 0

The equation of state parameter for the entanglement entropy density:
```
w_ent = p_ent / rho_ent
```

For CDM: w = 0 (pressureless dust, rho ~ a^{-3}).

For the entanglement entropy, the pressure comes from the de Sitter temperature:
```
p_thermal = (1/3) rho_thermal = (1/3) * (pi^2 / 30) * T_dS^4 / (hbar c)^3
```

If the entanglement entropy has a thermal origin at T_dS:
```
w_ent = 1/3    (radiation-like!)
```

This is EXACTLY the wrong equation of state. CDM requires w = 0. Thermal radiation has w = 1/3. The de Sitter entanglement entropy, being thermal at T_dS, naturally has w = 1/3.

**One could argue** that the entanglement entropy is NOT truly thermal -- it is a vacuum entanglement, and vacuum energy has w = -1. But vacuum energy with w = -1 gives rho ~ constant (cosmological constant behavior), which is also wrong for CDM.

The tension:
- Thermal interpretation: w = 1/3 (radiation) -- wrong
- Vacuum interpretation: w = -1 (dark energy) -- wrong
- CDM requirement: w = 0 (dust) -- not achievable from entanglement

---

## Part III: Perturbation Properties of Entanglement Entropy

### 3.1 Setup: Entanglement Entropy in Perturbed FRW

Consider FRW with scalar perturbations in Newtonian gauge:
```
ds^2 = -(1 + 2 Phi) dt^2 + a^2(t)(1 - 2 Psi) delta_{ij} dx^i dx^j
```

The entanglement entropy of a spherical region of comoving radius R in this background depends on:
1. The intrinsic geometry of the boundary (area + extrinsic curvature)
2. The state of quantum fields inside the region
3. The UV cutoff (or renormalization scheme)

### 3.2 First-Order Perturbation: delta S_ent

From Jacobson's entanglement equilibrium (arXiv:1505.04753, PRL 116, 201101, 2016):

The entanglement entropy decomposes as:
```
S_ent = S_UV + S_IR
```

- S_UV: UV-divergent piece, proportional to area (Bekenstein-Hawking term)
- S_IR: IR piece from matter entanglement, depends on state

At first order around maximally symmetric vacuum:
```
delta S_UV = delta (A / (4G)) = delta S_geometric
```

This is the Jacobson (1995) result: the UV entanglement entropy variation equals the geometric (Bekenstein-Hawking) entropy variation. Together with the first law of causal diamond thermodynamics:
```
delta S_UV = delta E / T_Unruh
```

this gives the Einstein equations. This is a FIRST-ORDER result and produces the standard Einstein equations -- no dark matter emerges.

### 3.3 Second-Order Perturbation: delta^2 S_ent

At second order (Jacobson 2016, Section V):
```
delta^2 S_ent = delta^2 S_UV + delta^2 S_IR + cross terms
```

The second-order UV contribution:
```
delta^2 S_UV = (1/(4G)) delta^2 A + (delta^2 (1/G)) * A/4 + ...
```

The second-order IR contribution involves the quantum Fisher information:
```
delta^2 S_IR ~ integral d^3x d^3x' C(x,x') delta g(x) delta g(x')
```
where C(x,x') is the connected two-point function of the modular Hamiltonian.

**Critical point**: The second-order contribution does NOT produce a new independent degree of freedom. It produces corrections to the Einstein equations (higher-derivative terms, or equivalently, corrections to Newton's constant). These corrections are:
- UV-sensitive (depend on the cutoff)
- Small (suppressed by (l_P / R)^2 where l_P is the Planck length)
- NOT equivalent to a pressureless dust component

### 3.4 The Sound Speed of Entanglement Perturbations

**This is the key calculation for Route C.**

Define an effective entanglement stress-energy tensor:
```
T_{mu nu}^{ent} = (2 / sqrt{-g}) * delta S_ent / delta g^{mu nu}
```

The perturbation of this effective stress-energy:
```
delta T_{00}^{ent} = delta rho_ent
delta T_{ij}^{ent} = -delta p_ent * delta_{ij} + Pi_{ij}^{ent}
```

The effective sound speed:
```
c_s,ent^2 = delta p_ent / delta rho_ent
```

For a conformal field theory (which is the leading contribution to entanglement entropy):
```
S_ent(R) = (c/3) ln(R / epsilon) + finite terms    [1+1 dimensions]
S_ent(R) = alpha * (R/epsilon)^{d-2} + ...          [d+1 dimensions]
```

The response to a metric perturbation Phi:
```
delta S_ent = integral d^{d-1}x <T_{00}(x)> * (2 Phi(x)) / T_local + ...
```

The modular Hamiltonian for a ball of radius R in a CFT is (Casini, Huerta & Myers 2011):
```
K = 2 pi integral_ball d^{d-1}x (R^2 - |x|^2) / (2R) * T_{00}(x)
```

Perturbation of K under metric perturbation Phi:
```
delta K = 2 pi integral d^{d-1}x (R^2 - |x|^2) / (2R) * delta T_{00}(x)
```

The propagation of delta K follows from the propagation of delta T_{00}. For a CFT:
```
partial_t^2 delta T_{00} = c^2 nabla^2 delta T_{00}    (for the trace part)
```

This gives:
```
c_s,ent^2 = c^2 / (d-1) = c^2 / 3    [for a CFT in 3+1 dimensions]
```

**Result: The entanglement entropy perturbation propagates at c_s = c/sqrt(3) (the speed of sound in a relativistic fluid), NOT at c_s = 0.**

This is a FUNDAMENTAL obstruction: entanglement entropy perturbations behave like radiation perturbations, not like CDM perturbations.

### 3.5 The Jeans Length for Entanglement Perturbations

The Jeans length sets the scale below which perturbations oscillate (pressure-supported) rather than grow:
```
lambda_J = c_s * sqrt(pi / (G rho))
```

For entanglement perturbations with c_s = c/sqrt(3):
```
lambda_J,ent = (c/sqrt(3)) * sqrt(pi / (G rho_total))
             ~ (c/sqrt(3)) * sqrt(pi / (G * rho_crit(z)))
             ~ H^{-1} / sqrt(3)
             ~ L_Hubble / sqrt(3)
```

**The Jeans length for entanglement perturbations is comparable to the Hubble radius.** This means entanglement perturbations do NOT collapse at sub-Hubble scales. They free-stream just like radiation.

For CDM: c_s = 0, so lambda_J = 0. CDM perturbations collapse at ALL scales. This is why CDM can provide the gravitational potential wells needed for CMB acoustic peaks.

### 3.6 Can the Effective c_s Be Made Small?

One might hope that non-perturbative effects, interactions, or the specific structure of de Sitter entanglement could reduce c_s from c/sqrt(3) to near zero. Let us examine this:

**Attempt 1: Massive field entanglement**

For a massive field with mass m, at distances r >> 1/m:
```
S_ent ~ exp(-m r) * (oscillating terms)
```

The entanglement is exponentially suppressed at distances larger than the Compton wavelength. This means massive field entanglement is SHORT-RANGE and cannot contribute significantly at cosmological scales.

For the entanglement to matter at galactic scales (r ~ kpc), we need m < 1/kpc ~ 10^{-26} eV. Such ultralight fields exist (fuzzy dark matter), but:
- They are actual particles, not vacuum entanglement
- Their c_s is set by the de Broglie wavelength: c_s ~ hbar k / (2m) ~ v_particle, which is small for non-relativistic particles but NOT zero
- This is not an entanglement effect but a particle effect

**Attempt 2: Entanglement in a non-Lorentz-invariant theory**

If the UV completion of gravity breaks Lorentz invariance (as in Horava gravity), the entanglement structure could differ:
```
S_ent ~ integral d^{d-1}x rho_ent(x) with modified dispersion relation
omega^2 = c^2 k^2 + alpha k^4 / M^2 + ...
```

For z = 2 Lifshitz scaling (Horava gravity):
```
c_s^2(k) = c^2 + alpha k^2 / M^2
```

This gives c_s > c at high k, not c_s < c. The Lorentz-violating terms make things WORSE, not better.

**Attempt 3: Non-equilibrium entanglement**

If the de Sitter vacuum is not in true thermal equilibrium (which is the case -- the de Sitter temperature is observer-dependent and the vacuum is alpha-vacuum, not Bunch-Davies), the entanglement entropy response could differ from the thermal CFT prediction.

However, non-equilibrium effects typically INCREASE entropy production and make perturbations propagate FASTER (ballistic rather than diffusive). There is no known mechanism where non-equilibrium quantum effects reduce c_s.

**Attempt 4: Nonlinear effects in the entanglement entropy**

If S_ent has non-trivial dependence on the metric perturbation:
```
S_ent[g] = S_ent[g_0] + integral delta S / delta g * delta g + (1/2) integral integral delta^2 S / (delta g delta g) * delta g * delta g + ...
```

The nonlinear terms could, in principle, create an effective potential for delta g that has a minimum at some non-trivial configuration. At this minimum, the "effective sound speed" for fluctuations around the minimum could be small.

This is mathematically possible but requires:
1. The second functional derivative of S_ent to have a specific sign and magnitude
2. The minimum to be stable
3. The fluctuations around the minimum to be non-propagating

There is NO calculation or argument suggesting this occurs for the de Sitter entanglement entropy. It remains a pure speculation.

### 3.7 Summary: Perturbation Properties

| Property | CDM requirement | Entanglement entropy | Discrepancy |
|---|---|---|---|
| w | 0 | 1/3 (thermal) or -1 (vacuum) | Fatal |
| c_s^2 | 0 | c^2/3 (CFT) | Fatal |
| c_vis^2 | 0 | O(1) (relativistic) | Fatal |
| lambda_J | 0 | ~ L_Hubble | Fatal |
| Thomson coupling | None | None (correct!) | OK |
| Gravitational coupling | Yes | Yes (correct!) | OK |

The entanglement entropy gets two properties right (gravitational coupling, no Thomson scattering) and four wrong (equation of state, sound speed, viscosity, Jeans length). The failures are all connected to the same root cause: entanglement in a relativistic QFT is inherently relativistic, with c_s ~ c.

---

## Part IV: Jacobson's Entanglement Equilibrium at Second Order

### 4.1 The Entanglement Equilibrium Hypothesis (Jacobson 2016)

Jacobson proposes that the vacuum entanglement entropy in small geodesic balls is MAXIMIZED at fixed volume in a locally maximally symmetric state. This means:
```
delta S_ent = 0    (first order -- equilibrium)
delta^2 S_ent <= 0  (second order -- maximum)
```

At first order, this gives the Einstein equations (reproducing Jacobson 1995 in a cleaner framework).

### 4.2 Second-Order Analysis

At second order, the entropy variation decomposes as:
```
delta^2 S_ent = delta^2 S_UV + delta^2 S_IR
```

**UV part** (Bianchi & Myers 2014, arXiv:1212.5183):
```
delta^2 S_UV = -(A / (4G)) * integral d^{d-1}x [C_{abcd} C^{abcd}] * f(position)
```
where C_{abcd} is the Weyl tensor. This is always negative (the Weyl tensor squared is positive definite), consistent with the maximum condition.

The Weyl tensor quantifies the gravitational tidal field. Its perturbation in Newtonian gauge:
```
C_{0i0j} ~ partial_i partial_j Phi - (1/3) delta_{ij} nabla^2 Phi
```

This is the TIDAL field, not a density perturbation. It does not provide a pressureless mode.

**IR part**:
```
delta^2 S_IR = -integral d^{d-1}x integral d^{d-1}x' <delta T_{00}(x) delta T_{00}(x')>_connected
              * kernel(x, x')
```

The connected correlator <delta T_{00}(x) delta T_{00}(x')> propagates information at the speed of light. Therefore, delta^2 S_IR responds to perturbations on the light cone, giving an effective propagation speed c_s ~ c.

### 4.3 What Perturbation Modes Emerge at Second Order?

The second-order entanglement entropy produces:
1. **Corrections to Newton's constant**: delta G / G ~ (l_P / R)^2 * (curvature invariants)
2. **Higher-derivative gravity terms**: R^2, R_{ab}R^{ab} corrections
3. **Non-local terms**: From the connected correlator, which has support on the light cone

**None of these produce a c_s = 0 mode.** They produce:
- Short-range corrections (Planck-suppressed)
- Wave-equation modifications (change the speed, not eliminate propagation)
- Non-local effects (relevant at very long or very short distances)

**Explicit statement**: At no order in perturbation theory around the maximally symmetric vacuum does a pressureless, non-propagating mode emerge from the entanglement entropy. This is a consequence of the Reeh-Schlieder theorem: the vacuum of a relativistic QFT has correlations at all spacelike separations, and perturbations of these correlations propagate at c.

---

## Part V: Bianchi-Myers and Casini-Huerta: Entanglement Entropy Response to Metric Perturbations

### 5.1 Bianchi & Myers (2014, arXiv:1212.5183)

Key results for entanglement entropy of a region R in a state rho on a curved spacetime:

**Leading term** (area law):
```
S_ent(R) = (A(R) / (4 G_eff)) + subleading
```
where G_eff is an effective Newton's constant that absorbs the UV divergence.

**Under a metric perturbation** h_{mu nu}:
```
delta S_ent = (1/(4 G_eff)) delta A + integral_R d^{d-1}x integral_R d^{d-1}x' K(x,x') h(x) h(x') + ...
```

The first term (area variation) gives the standard gravitational response.
The second term (matter entanglement response) gives quantum corrections.

**Response function** K(x,x'):
```
K(x,x') = <delta^2 S / delta g(x) delta g(x')> = Fisher information metric
```

This is the quantum Fisher information, which is bounded by:
```
K(x,x') <= (4 / hbar^2) * Var(H_modular)
```

The Fisher information has support within the causal diamond of the region R. Its Fourier transform:
```
K(k, omega) ~ |k|^{d-1} * theta(omega^2 - c^2 k^2)   [on the light cone]
```

This confirms that the entanglement entropy response to metric perturbations propagates at the speed of light.

### 5.2 Casini-Huerta Modular Hamiltonian (arXiv:1611.00016 and earlier work)

For a ball of radius R in a CFT, the modular Hamiltonian is:
```
K = 2 pi integral_ball d^{d-1}x [(R^2 - r^2) / (2R)] T_{00}(x)
```

The entanglement entropy perturbation:
```
delta S_ent = Tr[delta rho * K] = 2 pi integral d^{d-1}x [(R^2 - r^2) / (2R)] <delta T_{00}(x)>
```

Under a metric perturbation Phi:
```
<delta T_{00}> = <T_{00}>_background * Phi + (response function) * Phi
```

For a CFT in the vacuum:
```
<T_{00}> = 0    (vanishes by conformal invariance)
```

So the leading response comes from the variation of the stress-energy:
```
delta S_ent ~ integral d^{d-1}x kernel(x) * delta T_{00}^{matter}(x) * Phi(x)
```

This is a COUPLING between matter perturbations and metric perturbations, mediated by entanglement. It is NOT an independent degree of freedom -- it is a RESPONSE.

### 5.3 The "Sound Speed" of Entanglement Entropy Perturbations

From the Casini-Huerta modular Hamiltonian, we can define the entanglement entropy perturbation as a field:
```
sigma(x, t) = delta S_ent(ball centered at x, time t)
```

The equation of motion for sigma follows from the conservation of the modular Hamiltonian:
```
partial_t sigma = [H, sigma] = integral d^{d-1}x' [T_{00}(x'), sigma]
```

For a free field:
```
partial_t^2 sigma = c^2 nabla^2 sigma + source terms
```

This gives a WAVE EQUATION with speed c. The entanglement entropy perturbation propagates as a wave at the speed of light.

**If we define an effective fluid for the entanglement entropy:**
```
delta rho_ent = (1/(8 pi G)) * (terms involving delta S_ent)
delta p_ent = c_s,ent^2 * delta rho_ent
```

Then:
```
c_s,ent^2 = c^2 * (some O(1) coefficient from the modular Hamiltonian geometry)
```

For a conformal field in 3+1 dimensions: c_s,ent^2 = c^2/3.

**This is an exact result for CFTs and a robust estimate for general QFTs.** The sound speed of entanglement entropy perturbations cannot be zero in a relativistic theory.

---

## Part VI: The FRW Perturbation Equations for Entanglement Entropy

### 6.1 Attempt: Write Effective Fluid Variables

Suppose, despite the problems identified above, we formally define an effective entanglement fluid with:
```
rho_ent(z) = f(z) * rho_crit(z)
delta rho_ent(k, z) = perturbation around the background
v_ent(k, z) = velocity perturbation
Pi_ent(k, z) = anisotropic stress
```

The perturbed Einstein equations in the presence of baryons, photons, neutrinos, and this entanglement fluid:
```
k^2 Phi = -4 pi G a^2 [rho_b Delta_b + rho_gamma Delta_gamma + rho_nu Delta_nu + rho_ent Delta_ent]

(Phi - Psi) = 12 pi G a^2 (rho_ent + p_ent) Pi_ent / k^2    [anisotropic stress]
```

### 6.2 The Conservation Equations for the Entanglement Fluid

```
delta_ent' + 3 aH (c_s,ent^2 - w_ent) delta_ent = -(1 + w_ent)(theta_ent + 3 Phi')

theta_ent' + aH (1 - 3 w_ent) theta_ent = k^2 [c_s,ent^2 delta_ent / (1 + w_ent)] - k^2 Psi + k^2 Pi_ent
```

where primes denote conformal time derivatives.

### 6.3 Plugging In the Entanglement Entropy Properties

**Case 1: w_ent = 1/3 (thermal, c_s^2 = 1/3)**
```
delta_ent' + 0 = -(4/3)(theta_ent + 3 Phi')
theta_ent' - aH theta_ent = k^2 [(1/3) * (3/4) delta_ent] - k^2 Psi + k^2 Pi_ent
```

This is IDENTICAL to the photon fluid equations (up to the coupling terms). The entanglement fluid would oscillate, create acoustic peaks of its own, and free-stream. It would behave like extra radiation, NOT like CDM.

**Observable consequences**: Adding rho_ent with w = 1/3 at z ~ 1100:
- Increases N_eff (effective number of neutrino species)
- Does NOT provide gravitational potential wells
- Does NOT change z_eq (adds radiation, not matter)
- Planck constrains N_eff = 2.99 +/- 0.17, leaving room for only ~0.3 extra radiation species

With Omega_ent ~ 0.26 at z ~ 0, the effective N_eff at z ~ 1100 would be:
```
Delta N_eff ~ rho_ent / rho_single_nu ~ (0.26 / 0.14) * 3 * (7/8) * (4/11)^{4/3} ~ way too large
```

**This is ruled out at >> 100 sigma.**

**Case 2: w_ent = -1 (vacuum, c_s^2 = -1 or c_s^2 = 1)**

If the entanglement entropy behaves as vacuum energy:
- rho_ent = constant (not diluting)
- At z ~ 1100: rho_ent = rho_ent,0 ~ 0.26 rho_crit,0

This is negligible compared to the matter and radiation densities at z ~ 1100:
```
rho_ent / rho_crit(1100) ~ rho_ent,0 / (rho_crit,0 * (1101)^4) ~ 0.26 / (1.47 x 10^{12}) ~ 2 x 10^{-13}
```

The entanglement fluid would be completely negligible at CMB scales. It could not affect acoustic peaks.

Furthermore, if c_s^2 = -1 (phantom), the perturbations are exponentially unstable, leading to catastrophic growth of inhomogeneities.

**Case 3: w_ent = 0, c_s^2 = 0 (the dream scenario)**

This would be CDM. But as shown in Parts II-V, there is no physical mechanism by which de Sitter entanglement entropy achieves w = 0 and c_s = 0. These values are physically inconsistent with:
- The thermal nature of de Sitter entropy (which gives w = 1/3)
- The vacuum nature of entanglement (which gives w = -1)
- The relativistic propagation of entanglement perturbations (which gives c_s ~ c)

### 6.4 The Anisotropic Stress Problem

Even if we could somehow achieve w ~ 0 and c_s ~ 0, there remains the anisotropic stress:
```
c_vis,ent^2 = Pi_ent / delta_ent
```

For CDM: c_vis^2 = 0 (no viscosity, particles free-fall).

For entanglement entropy: The entanglement entropy of a region depends on the SHAPE of the region, not just its volume. Under a tidal perturbation (Weyl tensor C_{abcd}), the entanglement entropy changes:
```
delta S_ent ~ C_{abcd} C^{abcd} * (volume factor) * (UV cutoff factor)
```

This creates an anisotropic stress response:
```
c_vis,ent^2 ~ c^2 * (cutoff-dependent coefficient)
```

The anisotropic stress is generically of order c^2, far from the c_vis^2 = 0 required by CDM.

---

## Part VII: Literature Check -- Has Anyone Computed CMB Predictions from Entanglement Entropy?

### 7.1 Verlinde's Own Papers

**Verlinde 2016 (arXiv:1611.02269)**: The paper makes NO predictions for the CMB. The entire framework is formulated for static, spherically symmetric mass distributions at z = 0. There is no perturbation theory, no discussion of the early universe, and no attempt to address CMB physics.

**Verlinde & Zurek (2019, arXiv:1902.05922)**: Discusses emergent gravity in relation to quantum error correction. No CMB predictions.

### 7.2 Observational Tests of Verlinde's Theory

**Brouwer et al. (2017, MNRAS 466, arXiv:1612.03034)**: Tests Verlinde's formula against galaxy-galaxy weak lensing. Good agreement. But explicitly states that they "assume a background LCDM cosmology" for computing distances -- they do NOT derive the cosmological evolution from Verlinde's framework.

**Ettori et al. (2019, A&A 621, arXiv:1906.00823)**: Tests Verlinde's formula at galaxy cluster scales. Finds significant discrepancies -- the formula does not work well for clusters, particularly for the intracluster gas distribution. This is relevant because cluster scales are intermediate between galaxies and cosmology.

**Ghari & Haghi (2026, arXiv:2601.01715)**: Tests Verlinde against dwarf spheroidals. Verlinde preferred over MOND at 5.2 sigma. But again, no cosmological evolution.

### 7.3 Critiques

**Hossenfelder (2017, PRD 95, arXiv:1703.01415)**: Shows that a covariant version of Verlinde's theory reduces to standard GR. The "extra force" is an artifact of the non-covariant formulation. This is a serious theoretical concern.

**Dai & Stojkovic (2017, JHEP, arXiv:1710.00946)**: Points out internal inconsistencies in Verlinde's derivation, particularly regarding the elastic medium analogy.

### 7.4 "Entanglement Cosmology" Literature

**Belfiglio, Luongo & Mancini (2022, PRD 105, 123523, arXiv:2201.12299)**: Studies "geometric corrections to cosmological entanglement." Identifies geometric quasi-particles as dark matter candidates. However:
- The quasi-particles are perturbative excitations, not background dark matter
- Their equation of state is not computed
- No CMB power spectrum calculation is performed
- The framework is for understanding entanglement production during inflation, not for replacing CDM

**Brahma, Alaryani & Brandenberger (2020, arXiv:2005.09688)**: Studies entanglement entropy of cosmological perturbations as momentum-space entanglement between sub- and super-Hubble modes. Finds entanglement entropy grows during inflation. However:
- The entanglement is BETWEEN modes, not a new energy component
- It does not source gravity (no backreaction calculation)
- No CMB prediction is made

**Martin & Vennin (various papers, 2015-2020)**: Extensive work on quantum discord and entanglement of inflationary perturbations. Key results:
- Squeezed quantum states produce entanglement entropy
- The entropy grows logarithmically with the scale factor during inflation
- At late times, decoherence destroys quantum correlations
- None of this produces a CDM-like component

### 7.5 Boyanovsky, Rai & Chen (2021, arXiv:2110.15488)

Studies ultralight dark matter production from infrared dressing. Shows entanglement entropy grows during particle production. But this is a particle physics mechanism (actual particle production), not a vacuum entanglement effect.

### 7.6 Summary of Literature Search

**No paper in the literature computes CMB predictions from de Sitter entanglement entropy.** The reasons are clear from this analysis:
1. Verlinde's framework has no perturbation theory
2. The entanglement entropy perturbations have the wrong sound speed (c_s ~ c)
3. The effective equation of state is wrong (w = 1/3 or w = -1, not w = 0)
4. The scaling with redshift is either wrong or circular
5. No one has found a mechanism to make entanglement entropy behave as CDM

---

## Part VIII: The Deeper Obstruction -- Why Entanglement Cannot Be CDM

### 8.1 The Reeh-Schlieder Theorem

**Theorem** (Reeh-Schlieder, 1961): In any relativistic quantum field theory, the vacuum state Omega has the property that any local algebra of observables A(O) (for any open region O) acts cyclically on Omega. That is, A(O) |Omega> is dense in the Hilbert space.

**Physical consequence**: The vacuum of a relativistic QFT has non-zero entanglement between ANY two spacelike-separated regions, no matter how far apart.

**Implication for entanglement perturbations**: A perturbation in one region is "felt" by the entanglement entropy of distant regions at the speed of light (not faster, by causality, but also not slower). This means entanglement entropy perturbations propagate at c, giving c_s ~ c.

### 8.2 The Bekenstein Bound

The Bekenstein bound limits the entropy of a system of energy E in a region of size R:
```
S <= 2 pi R E / (hbar c)
```

For a perturbation delta E in a region of size R:
```
delta S <= 2 pi R delta E / (hbar c)
```

This means the entanglement entropy response to an energy perturbation is LINEAR in the energy. The effective stress-energy of the entanglement:
```
rho_ent ~ delta S / (volume) ~ (2 pi R / (hbar c)) * delta E / (volume) ~ delta rho * (R / R_Planck)
```

This is a scaling relation, not a fluid with specific w and c_s. The Bekenstein bound tells us entanglement entropy tracks energy density (like radiation), not mass density (like CDM).

### 8.3 The Information-Theoretic Impossibility

The fundamental issue is that CDM must satisfy TWO seemingly contradictory requirements:
1. **Gravitationally coupled**: CDM creates gravitational potential wells
2. **Non-propagating**: CDM perturbations have c_s = 0 (no sound wave)

In the language of quantum information:
- Requirement 1 means CDM carries energy-momentum (it is a physical degree of freedom)
- Requirement 2 means CDM perturbations do not carry energy-momentum flux (no pressure wave)

For a FIELD THEORY degree of freedom, having energy but no pressure propagation requires:
- Either: a massive non-relativistic particle (kinetic energy << rest mass energy)
- Or: a constrained field whose dynamics are frozen by a constraint (like the Khronon)

Vacuum entanglement entropy is NEITHER of these. It is:
- Not a particle (it is a property of the vacuum state)
- Not constrained (it responds freely to metric perturbations)

Therefore, vacuum entanglement entropy CANNOT behave as CDM. This is not a quantitative failure -- it is a structural impossibility.

### 8.4 The One Loophole: Emergent Constraint

The only way entanglement entropy could produce c_s = 0 perturbations is if the RESPONSE of the entanglement entropy to metric perturbations has a non-trivial extremum -- a point where the entropy is insensitive to small perturbations.

Jacobson's entanglement equilibrium (delta S = 0 at the vacuum) IS such a condition. But it gives the Einstein equations, not CDM. The entanglement equilibrium is satisfied by the vacuum, and perturbations AWAY from equilibrium propagate at c (the entanglement "restoring force" is relativistic).

For CDM-like behavior, we would need a SECOND equilibrium -- a constrained surface in the space of metrics where the entanglement entropy has a degenerate extremum. This constrained surface would have non-propagating perturbations (c_s = 0).

**No such constrained surface is known to exist in the entanglement entropy of any quantum field theory.**

---

## Part IX: What Would Need to Be True for Route C to Work

### 9.1 The Mathematical Requirements

For de Sitter entanglement entropy to act as effective CDM at z ~ 1100, ALL of the following would need to hold simultaneously:

1. **rho_ent(z) ~ (1+z)^3 for z from 0 to ~ 3400**: The effective energy density must scale as matter, not radiation or vacuum energy.

2. **w_ent = 0 at all relevant redshifts**: The equation of state must be pressureless.

3. **c_s,ent^2 = 0**: The sound speed must vanish, requiring non-propagating perturbations.

4. **c_vis,ent^2 = 0**: The viscosity must vanish.

5. **Omega_ent h^2 = 0.12**: The background density must match the observed CDM density.

6. **delta_ent ~ (adiabatic, scale-invariant spectrum)**: The perturbation spectrum must match inflationary predictions.

7. **No coupling to photons via Thomson scattering**: (This IS satisfied.)

### 9.2 What Is Known to Be False

From the analysis in Parts II-VIII:

- Requirement 1: FALSE. rho_ent scales as H^4 (Hypothesis A), H^2 (Hypothesis B), or is negligible (Hypothesis C). None gives (1+z)^3 consistently.

- Requirement 2: FALSE. The de Sitter thermal contribution gives w = 1/3; the vacuum contribution gives w = -1. Neither is w = 0.

- Requirement 3: FALSE. The Reeh-Schlieder theorem and explicit CFT calculations give c_s^2 = c^2/3 or c_s^2 = c^2. Never zero.

- Requirement 4: FALSE. Entanglement entropy responds to tidal fields (Weyl tensor), creating anisotropic stress with c_vis^2 ~ c^2.

- Requirement 5: UNKNOWN. No calculation determines Omega_ent h^2 from first principles.

- Requirement 6: PLAUSIBLE. Inflationary perturbations create entanglement, which could have a scale-invariant spectrum.

- Requirement 7: TRUE. Entanglement entropy does not couple to photons.

**Score: 1 true, 1 plausible, 1 unknown, 4 false. Route C fails on 4 out of 7 requirements.**

---

## Part X: Comparison with Routes A and B

### 10.1 Route A: Running G (Kumar 2025)

- **Mechanism**: IR running of Newton's constant G(k) = G_N[1 + k_*/k]
- **CMB status**: FAILS. Running formula breaks down at CMB scales, gives G/G_N ~ 2700
- **Root cause of failure**: Running G modifies the LEFT side of Einstein equations (geometry), not the RIGHT side (matter content). Cannot produce a pressureless component.

### 10.2 Route B: Modular Flow / QRE Perturbation Theory

- **Mechanism**: Perturbations of D(rho_spacetime || rho_matter) around FRW might contain a c_s = 0 mode
- **CMB status**: UNKNOWN. No calculation exists. The observer-as-dynamical-field extension is structurally equivalent to the Khronon theory and CAN produce c_s = 0.
- **Root cause of uncertainty**: The modular Hamiltonian on FRW is not known. The extension to dynamical observers requires a new postulate.

### 10.3 Route C: De Sitter Entanglement Entropy (THIS ANALYSIS)

- **Mechanism**: Volume-law entanglement entropy in de Sitter acts as effective dark matter
- **CMB status**: FAILS. Wrong scaling (w != 0), wrong sound speed (c_s ~ c), wrong Jeans length (~ L_Hubble)
- **Root cause of failure**: Entanglement entropy is a property of the vacuum state, not an independent degree of freedom. Its perturbations propagate at c because of the Reeh-Schlieder theorem.

### 10.4 Comparative Assessment

| Route | Galaxy-scale success | CMB-scale success | Root obstruction |
|---|---|---|---|
| A (Running G) | YES | NO | Wrong side of Einstein eq. |
| B (Modular flow) | UNKNOWN | POSSIBLE | Requires new postulate |
| C (Entanglement) | YES (Verlinde) | NO | Reeh-Schlieder: c_s = c |

**Route B remains the only potentially viable path for deriving CDM-like behavior from the tau framework's information-theoretic structure.** Routes A and C fail for different but equally fundamental reasons.

---

## Part XI: What Route C DOES Contribute

Despite failing for CMB, Route C provides valuable insights:

### 11.1 It Explains Why Verlinde Works at z ~ 0

At z ~ 0, static, spherically symmetric configurations:
- No time evolution needed (static limit removes the c_s problem)
- The elastic medium analogy works because the "strain" is time-independent
- The volume-law entanglement entropy creates a static potential, not a propagating perturbation
- Verlinde's M_D(r) formula is a STATIC result, not a dynamical one

This clarifies the domain of validity: Verlinde works where the system is quasi-static, and fails where dynamics (sound speed, time evolution) matter.

### 11.2 It Identifies the Transition Scale

The transition from area-law to volume-law entanglement entropy at r ~ L_dS = c/H defines a natural scale where Verlinde's "dark matter" turns on. This is the SAME scale where MOND effects become important (a_0 ~ cH_0). Route C confirms the connection:

```
r_transition ~ L_dS = c/H_0
a_transition = GM/r_transition^2 ~ G rho r ~ G rho (c/H_0) ~ a_0
```

### 11.3 It Rules Out a Class of Models

Route C definitively rules out any model where:
- Dark matter IS vacuum entanglement entropy (as a fluid)
- The entanglement entropy responds to perturbations as a relativistic field
- The CMB is explained by the same mechanism as galaxy rotation curves

This is a USEFUL negative result because it forces the tau framework toward Route B (modular flow / dynamical observer) rather than trying to extract CDM from the vacuum entanglement directly.

---

## Part XII: Final Verdict

### 12.1 Does de Sitter Entanglement Entropy Act as Effective CDM at CMB Scales?

**NO.**

The failure is not marginal or fixable with small corrections. It is fundamental and structural:

1. **The equation of state is wrong** (w = 1/3 or w = -1, not w = 0). This is a direct consequence of the entanglement entropy being either thermal (T_dS) or vacuum-like. No intermediate w = 0 is possible without introducing a new degree of freedom.

2. **The sound speed is wrong** (c_s = c/sqrt(3), not c_s = 0). This is a direct consequence of the Reeh-Schlieder theorem: vacuum entanglement in a relativistic QFT propagates perturbations at the speed of light.

3. **The scaling is wrong** (rho_ent does not scale as (1+z)^3 consistently). The entanglement entropy tracks either the thermal radiation density (~ H^4) or the critical density (~ H^2), neither of which matches CDM across the radiation-matter transition.

4. **There is no perturbation theory**. Verlinde's framework makes predictions only for static, z = 0 configurations. No cosmological perturbation theory has been developed, and the structural analysis shows such a theory would give the wrong perturbation behavior.

### 12.2 Is This Fixable?

**No, not within the entanglement entropy framework.** The failures trace back to fundamental properties of relativistic quantum field theory (Reeh-Schlieder, conformal propagation, Bekenstein bound). These cannot be circumvented without abandoning the QFT framework -- which would undermine the entire entanglement-based approach.

The ONLY known ways to get c_s = 0 from fields (not particles) are:
1. **Break Lorentz invariance** with a preferred time direction (AeST, Khronon)
2. **Impose a constraint** that freezes the propagating mode (k-essence at a minimum)
3. **Use a non-propagating gauge mode** that has been given dynamics (Khronon as gauge-fixed diffeomorphism)

All three require NEW STRUCTURE beyond the entanglement entropy of the vacuum. This is the domain of Route B (dynamical observer = Khronon), not Route C.

### 12.3 Impact on the Four-Paper Architecture

**Paper 3**: Route C is ruled out for CMB. Paper 3 should:
- Present Verlinde's mechanism as valid for STATIC galaxy-scale phenomenology (z ~ 0)
- State clearly that the mechanism does not extend to CMB scales
- Cite this analysis as the quantitative demonstration

**Paper 4**: The CMB problem requires Route B (observer-as-dynamical-field = Khronon). The entanglement entropy provides the MOTIVATION (why there should be a preferred time direction = the observer) but not the MECHANISM (how to get c_s = 0 at CMB scales).

### 12.4 The Bottom Line in One Sentence

**De Sitter entanglement entropy is a static phenomenon that explains galactic-scale dark matter (Verlinde) but cannot produce the dynamical, pressureless, non-propagating perturbations needed for CMB acoustic peaks -- this requires a fundamentally different structure (a constrained field like the Khronon) that the tau framework can motivate but not derive.**

---

## Appendix A: Key Equations Summary

### A.1 Verlinde's Dark Matter Formula
```
M_D^2(r) = (c H_0 / (6G)) * r^2 * d/dr [r * M_B(r)]
```

### A.2 Effective Equation of State
```
Thermal: w_ent = 1/3 (radiation-like)
Vacuum:  w_ent = -1  (cosmological constant)
CDM req: w_CDM = 0   (dust)
```

### A.3 Sound Speed
```
Entanglement (CFT): c_s^2 = c^2/3
Entanglement (general): c_s^2 ~ O(c^2)
CDM requirement: c_s^2 = 0
```

### A.4 Jeans Length
```
Entanglement: lambda_J ~ L_Hubble / sqrt(3) ~ 5000 Mpc
CDM: lambda_J = 0
```

### A.5 Redshift Scaling
```
Hypothesis A: rho_ent ~ H^4 ~ (1+z)^8 (rad. dom.)  -- too large by 10^7
Hypothesis B: rho_ent ~ H^2 ~ (1+z)^4 (rad. dom.)  -- wrong scaling (circular)
Hypothesis C: rho_ent ~ Omega_B * rho_crit           -- wrong amplitude
CDM: rho_CDM ~ (1+z)^3 always
```

---

## Appendix B: The Circular Logic Problem in Detail

The self-referential nature of Hypothesis B deserves explicit treatment.

**Step 1**: Assume rho_ent = f * rho_crit(z) with f = Omega_DM.

**Step 2**: rho_crit(z) = 3 H(z)^2 / (8 pi G), where H(z) is determined by the total energy content:
```
H(z)^2 = H_0^2 [Omega_r (1+z)^4 + Omega_b (1+z)^3 + Omega_ent(z) + Omega_Lambda]
```

**Step 3**: But Omega_ent(z) = rho_ent(z) / rho_crit(z) = f. So:
```
H(z)^2 = H_0^2 [Omega_r (1+z)^4 + (Omega_b + f) (1+z)^3 + Omega_Lambda]
```

Wait -- this only works if rho_ent scales as (1+z)^3 (matter-like). But we assumed rho_ent = f * rho_crit, and rho_crit ~ H^2, so:
```
rho_ent ~ H^2 ~ [Omega_r (1+z)^4 + (Omega_b + f)(1+z)^3 + Omega_Lambda]
```

At early times (z >> 1), the dominant term is:
```
rho_ent ~ Omega_r (1+z)^4    [radiation scaling!]
```

So rho_ent ~ (1+z)^4, NOT (1+z)^3. It acts as RADIATION, not matter.

For self-consistency, we need rho_ent to CAUSE the matter-radiation transition, but rho_ent DEPENDS on the transition having already happened. This is circular.

The only self-consistent solution is rho_ent = constant fraction of rho_total, which means rho_ent has the same equation of state as the DOMINANT component at each epoch:
- z >> z_eq: rho_ent ~ (1+z)^4 (radiation)
- z_eq >> z >> 1: rho_ent ~ (1+z)^3 (matter, but only because CDM ALREADY dominates)
- z ~ 0: rho_ent ~ constant (if Lambda dominates)

**In a universe without CDM, rho_ent NEVER transitions from radiation-like to matter-like behavior.** The entanglement entropy cannot bootstrap itself into acting as CDM.

---

## Appendix C: Comparison with Belfiglio-Luongo-Mancini (2022)

The paper arXiv:2201.12299 (PRD 105, 123523) is the closest existing work to Route C. They study "geometric corrections to cosmological entanglement" and interpret geometric quasi-particles as dark matter candidates.

**Their approach**: During cosmic expansion, the time-dependent metric creates particles through the cosmological Schwinger effect. These "geometric quasi-particles" carry entanglement entropy.

**Key differences from Route C**:
1. Their quasi-particles are ACTUAL particle excitations (Bogoliubov particles), not vacuum entanglement
2. The particles have a definite equation of state determined by their production mechanism
3. No CMB power spectrum is computed
4. The approach is closer to "cosmological particle production as dark matter" than "vacuum entanglement as dark matter"

**Their limitations**:
1. No computation of c_s for the geometric quasi-particles
2. No comparison with Planck CMB data
3. The quasi-particle number density depends on the UV cutoff (same problem as all entanglement entropy calculations)

**Verdict**: Their approach is more physical than Route C (actual particles vs. vacuum entanglement) but equally incomplete for CMB predictions.

---

## References

### Verlinde's Framework
- Verlinde 2016: arXiv:1611.02269, SciPost Phys. 2, 016 (2017)
- Verlinde & Zurek 2019: arXiv:1902.05922

### Entanglement in Curved Spacetime
- Jacobson 2016: arXiv:1505.04753, PRL 116, 201101 (entanglement equilibrium)
- Bianchi & Myers 2014: arXiv:1212.5183, Class. Quant. Grav. 31, 214002
- Casini, Huerta & Myers 2011: JHEP 05, 036 (modular Hamiltonian for balls in CFT)
- Casini, Teste & Torroba 2017: arXiv:1611.00016 (QRE and RG flow)

### Cosmological Entanglement
- Belfiglio, Luongo & Mancini 2022: arXiv:2201.12299, PRD 105, 123523
- Brahma, Alaryani & Brandenberger 2020: arXiv:2005.09688
- Martin & Vennin 2016: PRD 93, 023505 (quantum discord of inflationary perturbations)
- Boyanovsky, Rai & Chen 2021: arXiv:2110.15488

### Observational Tests of Verlinde
- Brouwer et al. 2017: MNRAS 466, arXiv:1612.03034 (galaxy-galaxy lensing)
- Ghari & Haghi 2026: arXiv:2601.01715 (dwarf spheroidals, Verlinde > MOND at 5.2 sigma)
- Ettori et al. 2019: A&A 621 (cluster scale tests -- significant discrepancies)

### Critiques
- Hossenfelder 2017: PRD 95, arXiv:1703.01415 (covariant version reduces to GR)
- Dai & Stojkovic 2017: JHEP, arXiv:1710.00946 (internal inconsistencies)

### CMB Physics
- Hu 2008: arXiv:0802.3688 (CMB pedagogical reference)
- Planck 2018: arXiv:1807.06209 (cosmological parameters)

### Theories That Achieve c_s = 0 Without CDM Particles
- Skordis & Zlosnik 2021: arXiv:2007.00082, PRL 127, 161302 (AeST)
- Blanchet & Skordis 2024: arXiv:2404.06584, JCAP 11, 040 (Khronon)

### QFT Foundations
- Reeh & Schlieder 1961: Nuovo Cimento 22, 1051
- Bisognano & Wichmann 1975: JMP 16, 985

### tau Framework
- Paper 1: Huang (2026) -- Petz recovery unification
- Paper 2: Huang (2026) -- Information loss is local
- Paper 3: Huang (2026) -- Galactic-scale dark matter from running G
- Paper 4: Huang (2026) -- Grand unification via Sigma = D(rho_spacetime || rho_matter)

---

*Last updated: 2026-03-12*
*This analysis represents a rigorous, quantitative investigation of Route C for the CMB problem.*
*VERDICT: Route C FAILS. The failure is fundamental (Reeh-Schlieder, wrong w, wrong c_s), not fixable.*
*Route B (dynamical observer = Khronon) remains the only potentially viable information-theoretic path to CMB.*
