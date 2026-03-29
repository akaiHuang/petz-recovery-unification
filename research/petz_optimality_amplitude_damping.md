# Petz Map Optimality for Amplitude Damping Channel

**Date**: 2026-03-29
**Status**: COMPLETE
**Reference**: Li et al., PRL 134, 200602 (2025)

## Summary of Findings

**The Petz recovery map with sigma = I/2 is NOT optimal for the amplitude damping channel.** In fact, it **hurts** both entanglement fidelity and state fidelity compared to doing nothing (R = identity).

## Setup

- **Channel**: Amplitude damping with transmissivity eta
  - E_0 = [[1, 0], [0, sqrt(eta)]]
  - E_1 = [[0, sqrt(1-eta)], [0, 0]]
- **Reference state**: sigma = I/2 (maximally mixed)
- **Input state**: rho = |+><+| (used in Tuna-9 experiments)
- **Tuna-9 eta values**: 0.1, 0.3, 0.5, 0.7, 0.9

## Key Results

### 1. Petz Recovery Kraus Operators (analytical)

For sigma = I/2:
```
R_0^P = [[1/sqrt(2-eta), 0], [0, 1]]
R_1^P = [[0, 0], [sqrt((1-eta)/(2-eta)), 0]]
```

### 2. Entanglement Fidelity: Petz HURTS

| eta | F_e(E) (no recovery) | F_e(Petz) | F_e(SDP optimal) | SDP = identity? |
|-----|---------------------|-----------|-----------------|----------------|
| 0.1 | 0.433114 | 0.377866 | 0.433114 | YES |
| 0.3 | 0.598861 | 0.504160 | 0.598862 | YES |
| 0.5 | 0.728553 | 0.622008 | 0.728552 | YES |
| 0.7 | 0.843330 | 0.751515 | 0.843342 | YES |
| 0.9 | 0.949342 | 0.906812 | 0.949331 | YES |

**Analytical proof**: F_e(E) = (1+sqrt(eta))^2/4, while F_e(R_P circ E) contains the term (1/sqrt(2-eta) + sqrt(eta))^2/4. Since 1/sqrt(2-eta) < 1 for all eta in (0,1), the Petz composition always has lower entanglement fidelity.

**SDP verification**: The optimal recovery for entanglement fidelity is R = identity (do nothing).

### 3. State Fidelity: Petz Also HURTS

| eta | F_s(no recovery) | F_s(Petz) | tau_Petz | tau_noR |
|-----|-----------------|-----------|----------|---------|
| 0.1 | 0.658114 | 0.614708 | 0.385292 | 0.341886 |
| 0.3 | 0.773861 | 0.710042 | 0.289958 | 0.226139 |
| 0.5 | 0.853553 | 0.788675 | 0.211325 | 0.146447 |
| 0.7 | 0.918330 | 0.866900 | 0.133100 | 0.081670 |
| 0.9 | 0.974342 | 0.952267 | 0.047733 | 0.025658 |

**Analytical**: F_s(|+>, E(|+>)) = (1+sqrt(eta))/2

### 4. Li et al. B >= 0 Condition: VIOLATED

The B matrix is NOT Hermitian (||B - B^dag|| ~ 0.4 to 1.3), so B >= 0 in the PSD sense is violated. However, all eigenvalues of B have non-negative real parts.

| eta | ||B - B^dag|| | min eig((B+B^dag)/2) | min Re(eig(B)) |
|-----|-------------|---------------------|---------------|
| 0.1 | 1.33e+00 | -0.160 | 0.000 |
| 0.3 | 1.15e+00 | -0.110 | 0.000 |
| 0.5 | 9.57e-01 | -0.072 | 0.000 |
| 0.7 | 7.43e-01 | -0.040 | 0.000 |
| 0.9 | 4.38e-01 | -0.013 | 0.000 |

### 5. Information-Theoretic Bounds

The Fawzi-Renner/Wilde bound states F(rho, R_P(E(rho))) >= exp(-Sigma), i.e., tau_Petz <= 1 - exp(-Sigma).

| eta | Sigma | tau_Petz | 1-exp(-Sigma) | Satisfied? |
|-----|-------|----------|--------------|-----------|
| 0.1 | 0.604 | 0.385 | 0.454 | YES |
| 0.3 | 0.485 | 0.290 | 0.384 | YES |
| 0.5 | 0.377 | 0.211 | 0.314 | YES |
| 0.7 | 0.260 | 0.133 | 0.229 | YES |
| 0.9 | 0.115 | 0.048 | 0.108 | YES |

The bound is always satisfied with significant margin. tau_Petz / (1-exp(-Sigma)) ranges from 0.44 to 0.85.

## Why Does Petz Hurt?

The Petz map with sigma = I/2 is designed for Bayesian retrodiction: given the output of E, infer what the input was, assuming the prior was sigma = I/2. For amplitude damping:

1. E preserves |0> perfectly
2. E damps |1> -> sqrt(eta)|1> + sqrt(1-eta)|0>
3. R_P tries to "undo" this by amplifying the |1> component
4. But this amplification also distorts |0> (which was already fine)
5. Net effect: more distortion than the original channel

The root cause is that **amplitude damping is asymmetric** (it has a fixed point |0>), while the Petz map for sigma = I/2 treats both basis states symmetrically.

## Implications for the Framework

### What tau_Petz Measures

tau_Petz = 1 - F(rho, R_P(E(rho))) is the **Bayesian retrodiction error**: how well can you infer the input state from the output, using the Petz map as the inference procedure.

This is a **valid** irrecoverability measure:
- It satisfies tau_Petz <= 1 - exp(-Sigma) (Fawzi-Renner bound)
- It is monotone in channel strength (more damping -> larger tau)
- It equals 0 for reversible channels
- It is computable in closed form

### What tau_Petz Does NOT Measure

tau_Petz does NOT measure the **optimal** irrecoverability. For the amplitude damping channel:
- tau_optimal (state fidelity) = (1 - sqrt(eta))/2 (just do nothing)
- tau_Petz > tau_optimal (Petz overestimates)

The ratio tau_Petz / tau_optimal ranges from ~1.13 (eta=0.1) to ~1.86 (eta=0.9).

### Recommendation for Tuna-9

The Tuna-9 values using tau_Petz are **conservative upper bounds** on irrecoverability. They correctly capture the qualitative behavior (more damping = more irrecoverability) but overestimate the magnitude by a factor of ~1.1-1.9.

For future work, consider:
1. **Using tau_noR = 1 - F(rho, E(rho))** as a tighter measure
2. **Optimizing sigma** to match the actual prior (not necessarily I/2)
3. **Using the rotated Petz map** (Junge et al., Wilde 2015) which is provably near-optimal

## Corrected Tuna-9 Values

| eta | tau_Petz | tau_noR | tau_opt (upper bound) | Sigma |
|-----|----------|---------|----------------------|-------|
| 0.1 | 0.385 | 0.342 | <= 0.342 | 0.604 |
| 0.3 | 0.290 | 0.226 | <= 0.226 | 0.485 |
| 0.5 | 0.211 | 0.146 | <= 0.146 | 0.377 |
| 0.7 | 0.133 | 0.082 | <= 0.082 | 0.260 |
| 0.9 | 0.048 | 0.026 | <= 0.026 | 0.115 |

## Files

- Main script: `numerical/petz_optimality_check.py`
- Debug scripts: `numerical/petz_optimality_check_v2.py`, `numerical/petz_optimality_check_v3.py`, `numerical/petz_optimality_debug.py`, `numerical/petz_optimality_final.py`
- Plot: `research/petz_optimality_final_analysis.png`
