# Route 2: Modular Flow / Dorau-Much Framework

**Date**: 2026-03-07
**Agent task**: Investigate whether the Dorau-Much (2025 PRL) modular flow framework can derive Sigma_grav = r_s/r from first principles

## Executive Summary

The Dorau-Much framework (arXiv:2510.24491) derives Einstein equations from QRE on bifurcate Killing horizons, providing the most rigorous existing bridge between quantum information and gravity. While it does not directly give Sigma_grav = r_s/r, extending their mathematical toolkit -- the entanglement first law, Tolman temperature relation, and gravitational redshift -- yields Sigma_grav = r_s/r/(1 - r_s/r), which equals r_s/r in the weak-field limit. This is the most promising derivation route.

## Key Findings

### The Dorau-Much Derivation Chain

1. **States**: Vacuum omega_0 and coherent excitation omega_phi on the Weyl algebra of a free scalar field restricted to right Rindler wedge
2. **Modular operator acts geometrically**: Delta_R^{it} = D_{2pi t} (affine dilations along horizon)
3. **QRE = Energy flux**: S^rel(omega_0 || omega_phi) = -2pi integral_H U <:T_{UU}:> dU dvol_S
4. **Energy flux = Area variation** via Raychaudhuri: delta A = -integral U R_{ab} xi^a xi^b dU dvol_S
5. **Equating with Bekenstein-Hawking**: alpha <:T_{ab}:> xi^a xi^b = R_{ab} xi^a xi^b
6. **Result**: Semiclassical Einstein equations R_{ab} - (R/2)g_{ab} + Lambda g_{ab} = 8pi <:T_{ab}:>

### Three Routes to Sigma_grav = r_s/r

#### Route A: Modular Channel (MOST PROMISING)

Define N_grav as restriction from algebra A(r) to A(infinity).

Using entanglement first law + Tolman temperature T(r) = T_H / sqrt(1 - r_s/r):

```
S^rel(r) = [2pi / (kappa sqrt(1-r_s/r))] * <delta E>_phi      (at radius r)
S^rel(inf) = [2pi/kappa] * sqrt(1-r_s/r) * <delta E>_phi       (at infinity, after redshift)
```

Energy redshifts: <delta E>_inf = sqrt(1-r_s/r) * <delta E>_phi

Therefore:
```
S^rel(r) / S^rel(inf) = 1/(1 - r_s/r)
Sigma = S^rel(r) - S^rel(inf) = S^rel(inf) * r_s/r / (1 - r_s/r)
```

**Fractional QRE loss**: (S^rel(r) - S^rel(inf)) / S^rel(r) = r_s/r  [EXACT]

#### Route B: Bianconi Geometric QRE

Bianconi's action S_GfE = (1/l_P^d) integral sqrt{|g|} L d^d x with L = -Tr ln(G_tilde g_tilde^{-1})
- Local: L(r) ~ 12 beta_G r_s / r^3 (like Riemann components)
- Integrated: gives logarithmic dependence, NOT 1/r
- **Does NOT straightforwardly give r_s/r**

#### Route C: Pikovski-type (indirect)
- Total decoherence ~ integral of potential = r_s/r
- But phenomenological, not rigorous (see Route 1 report)

### The Local vs Integrated Problem

**The central technical obstacle**: r_s/r is an integrated (cumulative) quantity, but most QFT frameworks give local quantities (R_{abcd}, T_{ab}).

**Resolution**: Araki relative entropy IS naturally integrated -- it compares states on an algebra associated with a REGION, not a point. The modular Hamiltonian K_r = (2pi/kappa) H_r integrates from r to infinity. This is why Route A works.

### Critical Subtlety: Schwarzschild Does NOT Saturate the Petz Bound

From the exact result: Sigma = r_s/r / (1 - r_s/r)
- Petz bound: F >= exp(-Sigma/2) = exp(-r_s/(2r(1-r_s/r)))
- Schwarzschild fidelity: F_Schw = sqrt(1-r_s/r) < exp(-r_s/(2r))
- **Schwarzschild gives MORE information loss than the Petz minimum**
- The exponential metric saturates by construction
- This is a testable prediction: nature could sit at saturation (exponential) or beyond (Schwarzschild)

## Key Equations

```
Dorau-Much:  S^rel = -2pi integral_H U <:T_{ab}:> xi^a xi^b dU dvol_S

Modular Hamiltonian:  K_r = (2pi/kappa) H_r

Tolman temperature:  T(r) = T_H / sqrt(1 - r_s/r)

Fractional QRE loss:  [D_in - D_out] / D_in = r_s/r   [EXACT]

Entropy production:  Sigma(r) = r_s/r / (1 - r_s/r)   [exact]
                              ~ r_s/r                    [weak field]
```

## Key References

- Dorau & Much (PRL 2026): arXiv:2510.24491 -- QRE to Einstein equations
- Bianconi (PRD 2025): arXiv:2408.14391 -- Gravity from entropy
- Bianconi (Entropy 2025): arXiv:2501.09491 -- Schwarzschild QRE area law
- Modular Channels paper (2025): arXiv:2504.20457 -- Thermal filtering
- Casini, Grillo, Pontello (PRD 2019): arXiv:1903.00109 -- Relative entropy from Araki formula
- Jacobson (PRL 2016): arXiv:1505.04753 -- Entanglement equilibrium
- Wall (PRD 2012): arXiv:1105.3445 -- GSL from QRE monotonicity
- Witten (2018): arXiv:1803.04993 -- Entanglement in QFT
- UDW detector paper (PLB 2025): arXiv:2501.00229 -- QRE thermalization in Schwarzschild
- Hollands, Longo (2021): arXiv:2107.06787 -- Relative entropy in curved spacetimes
- Witten (2022): arXiv:2112.12828 -- Gravity and the crossed product
- Chandrasekaran et al. (2023): arXiv:2306.07323 -- Crossed product algebras

## Assessment

- **Likelihood of deriving Sigma_grav = r_s/r**: HIGH (fractional QRE loss = r_s/r exactly)
- **What works**: Entanglement first law + Tolman + redshift gives r_s/r as fractional distinguishability loss; the mathematics is rigorous and uses established AQFT results
- **What doesn't**: The absolute Sigma depends on the test state (not universal without normalization); making N_grav a precise CPTP map requires crossed product construction (Witten 2022)
- **The normalization question**: Defining Sigma_grav as the fractional loss (not absolute loss) is physically well-motivated but requires justification

## Next Steps

1. **For Paper 2**: Frame Sigma_grav = r_s/r as the "normalized QRE drop in the first-law regime of AQFT on Schwarzschild" -- stronger than mere ansatz
2. **Compute explicit N_grav**: Use crossed product algebra (Witten 2022, Chandrasekaran et al. 2023) to make the channel precise
3. **Test non-saturation prediction**: Compare Schwarzschild vs exponential metric predictions against UDW detector data (arXiv:2501.00229) and Pino ion trap data (2025)
4. **Upgrade the derivation**: Move from coherent state (first-law regime) to general states using full modular theory
