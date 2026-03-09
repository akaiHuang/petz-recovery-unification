# Path B Correction: tau_min > 0 Overclaim

## Date: 2026-03-07
## Status: Three research agents completed, overclaim identified and corrected

---

## The Overclaim (WRONG)

"Due to gravity, tau_min > 0 for ALL particles with E > 0. Perfect quantum coherence is impossible in principle."

## Why It's Wrong

### 1. Pikovski mechanism requires internal DOF (N > 0)
- H = p^2/(2m) + H_int(1 + Phi(x)/c^2)
- If N = 0 (no internal DOF), H_int is trivial → D(t) = 1 → tau = 0 identically
- Electron has N = 0 → tau_Pikovski = 0 (exact)
- Photon: not applicable (massless, no rest frame)

### 2. DSW mechanism is proportional to m^2
- DSW gravitational decoherence rate ~ G^2 m^2 d^2 / L_dS^3
- For photon (m = 0): rate = 0 exactly
- Physical reasons:
  - Massive particles produce static Coulombic gravitational fields
  - Photons produce null pp-wave fields (qualitatively different)
  - Co-propagating photons don't gravitationally interact (Tolman-Ehrenfest-Podolsky 1931)
  - DSW rate vanishes for null sources

### 3. Even for massive particles, DSW is empirically vacuous
- Electron: t_decoherence ~ 10^180 seconds (10^163 × age of universe)
- 1 kg object: t_decoherence ~ 10^98 seconds
- Measurable only for m ~ 10^46 kg (a small galaxy!)
- Not new physics — standard QM + cosmological horizon

## The Correct Statement

> For **composite objects with N > 0 internal degrees of freedom** in a spatial superposition within a gravitational potential gradient, tau_Pikovski > 0 and scales as 1 - exp(-N*Gamma*t^2). The quantum-to-classical transition is driven by **internal complexity N**, not mass alone. Elementary particles with no internal structure (N = 0) have tau_Pikovski = 0 and are NOT subject to this mechanism, consistent with matter-wave universality.

## Classification by Particle Type

| Particle | N (internal modes) | tau_Pikovski | Matter wave |
|----------|-------------------|-------------|-------------|
| Photon | 0 | Not applicable | Perfect |
| Electron | 0 | 0 (exact) | Perfect |
| Proton | ~0 | ~0 | Perfect |
| Atom (H) | Small | Very small | Perfect |
| C60 | 174 | Small but > 0 | Good |
| Nanoparticle (10^6 atoms) | ~10^6 | Significant | Difficult |
| Macroscopic | ~10^23 | ~1 | Impossible |

## Physical Insight (Corrected)

**Gravity acts as CATALYST, not CAUSE, of decoherence.**

1. Gravitational time dilation (catalyst): creates position-dependent clock rate
2. Internal degrees of freedom (fuel): provide clock that CAN run at different rates
3. Entanglement generation (process): spatial DOF become entangled with internal DOF
4. Partial trace (result): tracing out internal DOF produces mixed spatial state

Without internal structure, gravity has nothing to decohere.

## Distinction from Penrose (now cleaner)

| | tau framework | Penrose |
|---|---|---|
| Determines collapse | N (internal complexity) | m (mass alone, E_G ~ m^2) |
| Elementary particles | tau = 0 (perfect quantum) | tau > 0 (but t_c ~ 10^14 yr) |
| Same mass, different N | Different tau | Same tau |
| Modifies QM? | No | Yes |

## What's Genuinely New in Paper 1b

1. Mapping Pikovski channel through Petz recovery (new calculation)
2. tau-R bridge: irreversibility and objectivity emerge together (new link)
3. F_spatial >= exp(-Sigma_grav/2): universal recoverability bound (new application)
4. Bridge relation tau_I = 2*tau_S*(1-tau_S) (new explicit form)
5. Operational collapse = Petz recovery failure (new interpretation)

Underlying physics is Pikovski (2015). Information-theoretic wrapping is new.

## Paper 1b Corrections Needed

Line 404-407 currently says:
"tau_grav constitutes the irreducible floor of temporal asymmetry: for any massive system..."

Should say:
"For composite objects with N > 0 internal degrees of freedom, tau_grav > 0 constitutes an irreducible floor that grows with internal complexity N. Elementary particles with N = 0 have tau_Pikovski = 0, consistent with matter-wave universality."

## Key References from This Round

- Satishchandran (2025, arXiv:2508.20171): DSW review, confirms <N> ~ m^2
- Linton & Tiwari (2025, arXiv:2501.18111): confirms m^2 scaling
- Oniga & Wang (2016, PRD 93 044027): photon gravitational decoherence ~ t_P^2 omega^3
- Lagouvardos & Anastopoulos (2021, CQG 38 115012): photon grav decoherence master equation
- Xu et al. (2019, Science 366 132): Micius satellite — no photon decoherence detected
- Tolman, Ehrenfest, Podolsky (1931, Phys Rev 37 602): co-propagating light beams don't gravitationally interact
- Ran Li (2025, PRD 111 044022): precise DSW rate in de Sitter

## Lesson Learned

The logical error was: "gravity couples to all energy" → "therefore gravitational decoherence affects all particles." The second step does not follow from Pikovski. Gravity coupling universally is the equivalence principle (true). But Pikovski decoherence requires BOTH gravitational coupling AND internal structure. Without internal DOF, the channel is trivial.
