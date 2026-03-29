# ABCD Verification Plans — Summary Index

Created: 2026-03-16

## Overview

Four parallel research directions for extending Paper 1 "The Arrow of Time from Petz Recovery":

| Option | Goal | Key Test | Timeline | Go/No-Go |
|--------|------|----------|----------|----------|
| A | Prove τ_comp ≥ Σf(Cᵢ) | α*(d=2) > 0.01 in 10⁵ random trials | 4-8 months | Week 1-2 |
| B | Generalize Saturation to σ≠I/d | Jensen gap compensation condition | 6-8 weeks | Week 3 |
| C | Continuous τ(t) via Lindbladian | τ(t)=0 for unitary < 10⁻¹² | 14 hours | Phase 1 |
| D | EDA tau-Chrono Engine | 23% improvement on 5-gate circuit | 7 weeks | Week 4 |

## Critical Discovery (Option B)

For d=2 dephasing at θ=π/4: saturation holds for ALL σ = diag(p, 1-p), even when [ρ,σ]≠0.
The σ-dependence cancels completely from both F² and exp(-ΔD).

New key quantity: **Jensen gap** = ⟨ψ|σ|ψ⟩ / exp(⟨ψ|ln σ|ψ⟩) ≥ 1
- Equality iff [ρ,σ]=0 (ρ is eigenstate of σ)
- For saturation with general σ: channel must compensate this gap exactly

## Intersection (from earlier analysis)

All four converge on: √τ_total = min Σ√τᵢ^eff(σᵢ) subject to saturation constraints

- A: Lower bound on this minimum
- B: Constraints (saturation surface)
- C: Dynamics (dτ/dt descent direction)
- D: Algorithm (numerical solver)

## Detailed Plans

Full plans are in the agent output files. Key execution files:
- Option C: 5 test cases with exact parameters
- Option D: 40+ unit tests, 4 benchmarks, 10 E2E checks
- Option A: 7 mathematical steps with numerical validation code
- Option B: 320K trial scan across 8 σ-types × 4 dimensions
