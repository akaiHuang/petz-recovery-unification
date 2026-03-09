# Route 1: Pikovski Gravitational Decoherence Channel

**Date**: 2026-03-07
**Agent task**: Investigate whether the Pikovski gravitational decoherence channel can DERIVE (not assume) Sigma_grav = r_s/r

## Executive Summary

The Pikovski channel (Nature Physics 11, 668, 2015) models gravitational time dilation as a dephasing channel on center-of-mass degrees of freedom. Direct computation shows that this channel **cannot** derive Sigma_grav = r_s/r because its entropy production is inherently probe-dependent (depends on Delta_E, mass, interaction time). However, r_s/r emerges naturally as the universal gravitational parameter that multiplies all probe-dependent factors, providing strong physical motivation for the ansatz.

## Key Findings

### The Pikovski Channel Structure

- **Type**: Dephasing channel in position basis
- **Hamiltonian coupling**: H_int = Phi(x) * H_0 / c^2 (gravity couples to total mass-energy)
- **Kraus operators**: K_n = sqrt(p_n) * exp(i * Phi(x) * E_n * t / (hbar * c^2))
- **Decoherence function**: Gamma(x,x',t) = Tr_int[rho_int * exp(i * Delta_Phi * (x-x') * H_0 * t / (hbar * c^2))]
- **Decoherence is Gaussian** (~exp(-t^2)), NOT exponential, because the evolution is non-Markovian

### Why Direct Derivation Fails

1. **System-dependence**: Sigma_Pikovski ~ (Delta_E * Delta_Phi * t / (hbar * c^2))^2 depends on probe properties
2. **Supremum diverges**: sup over all probes gives infinity (take Delta_E -> infinity)
3. **Normalized rate gives r_s/(2r)**, not r_s/r, and the normalization is ad hoc

### What r_s/r Actually Is in Pikovski

The dephasing rate per unit frequency per unit time equals Phi(x)/c^2 = r_s/(2r). This is the universal gravitational parameter that appears in ALL Pikovski-type calculations, regardless of probe. But it enters as a multiplicative factor in the dephasing angle, not as the entropy production itself.

## Key Equations

```
N_Pikovski(rho_cm) = sum_n p_n V_n rho_cm V_n^dag
   where V_n = exp(i Phi(x_hat) E_n t / (hbar c^2))

Sigma_Pikovski ~ (Delta_E * Delta_Phi * t / (hbar c^2))^2   [probe-dependent]

d(theta)/dt / omega = Phi(x)/c^2 = -r_s/(2r)               [universal, but not Sigma]
```

## Three Alternative Routes Identified

1. **Modular channel on local Rindler horizons** (extends Dorau-Much 2025): Define gravitational channel as restriction to causal diamond, compute QRE geometrically. Probe-independent by construction. [MOST PROMISING]

2. **Integrated proper-time deficit**: Total proper-time deficit at radius r relative to infinity ~ T * r_s/(2r). Identify entropy production with "information deficit" per unit reference time. [PHYSICAL BUT AD HOC]

3. **First law of entanglement entropy** (holographic): In AdS/CFT, delta S_rel = delta <K> - delta S_EE. For perturbation creating mass M, relative entropy on ball of radius r ~ r_s/r via JLMS result. [DEEP BUT REQUIRES AdS/CFT]

## Key References

- Pikovski et al. (2015): arXiv:1311.1095 -- Gravitational decoherence
- Dorau & Much (2025 PRL): arXiv:2510.24491 -- QRE to Einstein equations
- Moreira & Celeri (2024): arXiv:2407.21186 -- Graviton bath entropy production (no r_s/r connection)
- Basso, Maziero & Celeri (2025 PRL 134, 050406): arXiv:2405.03902 -- Quantum Crooks in curved spacetime (observer-dependent Sigma)
- Balatsky et al. (2025 PRA 111, 012411) -- Gravitational dephasing of qubits
- Santos et al. (2019) npj QI 5, 23 -- Relative entropy of coherence = dephasing entropy production

## Assessment

- **Likelihood of deriving Sigma_grav = r_s/r via Pikovski**: LOW (fundamental obstacle: probe-dependence)
- **What works**: r_s/r is the correct universal parameter; physical motivation is solid
- **What doesn't**: Cannot get a probe-independent, dimensionless entropy production from this channel alone
- **Verdict**: Pikovski is excellent MOTIVATION but not a DERIVATION route

## Next Steps

1. Pursue Route 2 (modular channel / Dorau-Much extension) as the primary derivation strategy
2. In Paper 2, cite Pikovski as physical motivation, not as derivation
3. The proper-time deficit argument may yield a semi-rigorous "physical derivation" even if not mathematically complete
4. Consider JLMS holographic route for a future paper (requires AdS/CFT machinery)
