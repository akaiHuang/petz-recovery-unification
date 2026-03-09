# "Classical Gravity" Consistency Analysis

## Date: 2026-03-07
## Status: 4 research agents completed, results synthesized

---

## Core Question

If we commit to "gravity is (at least partially) classical" via Galley-Giacomini-Selby theorem, does this create logical contradictions with:
1. Penrose's own theory?
2. Our existing 4-paper framework?
3. Strong-field physics (Paper 2)?
4. Weak-field physics (Paper 3)?

---

## Executive Summary

| Comparison | Conflict Level | Resolution |
|-----------|---------------|------------|
| Penrose ↔ Classical gravity | **None** | Penrose says gravity is "neither classical nor quantum" — part of a new theory |
| Paper 1 ↔ Paper 1b | **None** | Paper 1 is gravity-agnostic; Paper 1b specializes to gravitational channel |
| Paper 2 ↔ Paper 1b | **Mild** | Paper 2 is semiclassical (QFT on classical background) — this IS Paper 1b's regime |
| Paper 3 ↔ Paper 1b | **Mild** | Kumar's QFT corrections = matter loops on classical background = perfectly fine |
| Paper 4 ↔ Paper 1b | **FUNDAMENTAL** | Σ = D(ρ_spacetime ‖ ρ_matter) requires spacetime to be a density matrix = quantum |

### Bottom Line
Papers 1-3 are safe. Paper 4 needs rethinking — either adopt Bianconi's "metric as density matrix" formalism, or use scale-dependent ontology.

---

## Part 1: Penrose Compatibility

### Penrose's actual position (NOT "gravity is classical")

Penrose does NOT say gravity is classical. His position is a **third option**:
- "Gravitization of quantum mechanics" (NOT "quantization of gravity")
- Both QM and GR are approximations of some **new, unknown theory**
- Spacetime superpositions exist **briefly** before collapse resolves them
- The instability (E_G) is what triggers collapse

### Two readings of "gravity is classical"

| | "Gravity is fundamentally classical" (Oppenheim) | "Gravity behaves classically because collapse prevents superposition" (Penrose) |
|--|---|---|
| Status of classicality | Axiom / Input | Emergent / Output |
| Collapse mechanism | Decoherence from CQ coupling | GR-QM incompatibility (E_G) |
| Is E_G meaningful? | Not directly | Yes (transient superposition) |
| Compatible with our τ? | YES | YES |

### For our framework: either version works
- Under Oppenheim: Classical gravity (input) → irreversible coupling → τ > 0. Self-consistent.
- Under Penrose: GR-QM tension → collapse (τ > 0) → effective classicality. Self-consistent.
- We need the **operational claim** ("gravity behaves classically"), not the ontological one.

### Important recent result
Grossardt et al. (2025, PRD 111, L121101): Diósi-Penrose model with classical gravity CAN generate entanglement. The BMV experiment alone cannot distinguish classical from quantum gravity. The boundary is blurrier than assumed.

---

## Part 2: Strong-Field (Paper 2) Compatibility

### Summary: Almost no conflict

| Paper 2 Element | Requires quantum gravity? | Impact of classical gravity |
|----------------|--------------------------|---------------------------|
| Exponential metric | No (classical solution) | None |
| Σ_grav = r_s/r | No (QFT on classical background) | None |
| No-horizon argument (Layer 1) | No (Petz bound is universal) | **Strengthened** (τ > 0 genuine) |
| Petz saturation (Layer 2) | No (channel property) | Pre-existing issue unchanged |
| Dorau-Much derivation | No (semiclassical Einstein eq.) | None |
| Holographic connections | Yes (full AdS/CFT) | **Weaken to analogy** |
| Page curve, ER=EPR | Yes | Paper 2 doesn't use them |

### Key insight: Dorau-Much is SEMICLASSICAL
Dorau-Much (2025 PRL) derives **semiclassical** Einstein equations from QRE. The quantum information is in the matter fields, not in the gravitational field. The spacetime geometry remains classical throughout.

This is analogous to Jacobson (1995): "It may be no more appropriate to canonically quantize the Einstein equation than it would be to quantize the wave equation for sound in air."

### The "semiclassical sweet spot"
Papers 1b and 2 both live in the regime:
1. Spacetime geometry is classical
2. Matter is quantum (full QFT)
3. Coupling is fundamentally irreversible (Galley-Giacomini-Selby)
4. Einstein equations emerge from quantum information of matter
5. Petz recovery bound constrains the geometry

**Together**: 0 < τ < 1 for all gravitational processes
- Paper 1b provides the lower bound: τ > 0 genuinely (from classical gravity)
- Paper 2 provides the upper bound: τ < 1 everywhere (from Petz bound + finite Σ)

### What changes in interpretation (not math)

| Statement | If gravity quantum | If gravity classical |
|-----------|-------------------|---------------------|
| τ > 0 | FAPP | Genuine |
| τ < 1 | Fundamental | Fundamental |
| No event horizon | Unitarity preserved | Information partially lost (bounded) |
| Holographic connection | Full AdS/CFT | Semiclassical holography only |

### One genuine tension: Saturation hypothesis
If the quantum-classical channel has fundamentally different properties from a purely quantum channel, the saturation condition F = exp(−Σ/2) might not hold. But this tension already exists independently — no known channel exactly saturates JRSWW for Σ > 0.

---

## Part 3: Weak-Field (Paper 3) Compatibility

### Summary: Safe, even strengthened

| Paper 3 Element | Requires quantum gravity? | Impact of classical gravity |
|----------------|--------------------------|---------------------------|
| Kumar's log correction | No (EFT matter loops) | **Perfect fit** |
| Running G | No (RG flow of matter) | Compatible |
| SPARC fits (Gubitosi) | No | None |
| Verlinde's mechanism | Requires entanglement | See below |
| Unified Σ equation | Depends on formulation | See below |

### Kumar's QFT correction: The QED analogy

| QED | Gravity (Kumar) |
|-----|-----------------|
| Classical Coulomb: V = −e²/r | Classical Newton: Φ = −GM/r |
| Quantum correction (Uehling): δV ~ α·ln(r) | Quantum correction: δΦ ~ G k_* ln(r) |
| EM field: classical background | Gravitational field: classical background |
| EM field need NOT be quantum | Gravity need NOT be quantum |

**Verdict**: Classical gravity + quantum matter naturally generates quantum corrections. This is standard EFT (Donoghue 1994).

### Verlinde's entanglement: WHOSE entanglement?

Verlinde's entanglement is of **microscopic degrees of freedom / vacuum fields**, not of the gravitational field itself. Gravity **emerges** from this entanglement; it is not entangled.

Analogy: Pressure is a classical, emergent force arising from quantum-mechanical molecular collisions. Similarly, gravity can emerge from entanglement while remaining classical.

**Verdict**: Compatible with classical gravity.

### The unified equation needs reinterpretation

**Problem**: Σ = D(ρ_spacetime ‖ ρ_matter) requires spacetime to be a density matrix.

**Three escape routes**:

1. **Bianconi's "metric as density matrix"** (PRD 2025): Metric has algebraic properties of density matrix without being ontologically quantum.

2. **Oppenheim's hybrid CQ formalism**: Classical probability distribution over metrics, lifted to hybrid density matrix.

3. **Reinterpret Σ**: Use Σ = D(ρ_matter(g) ‖ ρ_matter(g_flat)) — QRE of matter in curved vs flat background. This is what Dorau-Much actually compute.

**Recommendation for Paper 3**: Use Route 3. Save the full Bianconi formalism for Paper 4.

### Galley-Giacomini-Selby actually STRENGTHENS Paper 3
The theorem says: classical gravity + quantum matter = irreversibility.
Our framework says: Σ > 0 = irreversibility = time arrow.
**These are the same statement.** The theorem validates our core physical insight.

---

## Part 4: Internal Consistency of 4+1 Paper Framework

### The hierarchy of gravity assumptions

| Paper | Gravity assumption | Strength needed |
|-------|-------------------|-----------------|
| Paper 1 | None (gravity-agnostic) | Weakest |
| Paper 1b | Gravity behaves classically | Moderate |
| Paper 2 | Semiclassical (QFT on classical background) | Moderate |
| Paper 3 | EFT (classical background + matter loops) | Moderate |
| Paper 4 | Spacetime as density matrix | **Strongest** |

### The fundamental conflict: Paper 4

Paper 4 wants Σ = D(ρ_spacetime ‖ ρ_matter) — this treats spacetime as a quantum object.
Paper 1b wants gravity to be classical — spacetime is not quantum.

**These are mutually exclusive in naive formulation.**

### Resolution strategies

**Strategy A: Bianconi formalism**
- Metric plays the role of density matrix (formal analogy, not ontological claim)
- The QRE is computed on the metric's algebraic structure
- Classical gravity, formally quantum-like description

**Strategy B: Scale-dependent ontology**
- Below Planck scale: gravity may be quantum (Paper 4)
- Above Planck scale: gravity behaves classically (Papers 1b, 2, 3)
- The transition is where τ framework gives new physics

**Strategy C: "Conditional" framing**
- Papers 1-3: "If gravity is classical, THEN..."
- Paper 4: "If gravity is quantum, THEN the full unification..."
- Present both as complementary descriptions
- Let experiments decide (BMV, etc.)

**Recommended: Strategy B (scale-dependent)**
- Most physically motivated
- EFT is inherently scale-dependent
- Matches the existing narrative arc (different scales in different papers)

### The circularity question

Is this circular?
```
Quantum gravity → collapse → classical gravity → genuine τ > 0 → collapse
```

**This is a fixed-point argument, not circular.** The classicality is a self-consistent fixed point: the system "settles into" classical gravity through the collapse mechanism, and this classical behavior sustains the collapse. Like a thermostat: the temperature is what it is because the heating/cooling cycle maintains it there.

**BUT**: Galley-Giacomini-Selby requires **fundamentally** classical gravity, not merely effectively classical from decoherence. So the loop as stated doesn't quite work for GGS specifically.

**Resolution**: Use the operational version. Whether gravity is fundamentally or effectively classical, the operational consequences (τ > 0 for macroscopic systems) are the same. The distinction matters philosophically but not experimentally.

---

## Part 5: Impact on Experimental Predictions

### BMV experiment implications

| BMV Result | Impact on our framework |
|-----------|----------------------|
| Positive (entanglement found) | Gravity has quantum aspects → Paper 4 strengthened, Paper 1b weakened |
| Negative (no entanglement) | Gravity may be classical → Paper 1b strengthened |
| Ambiguous (Aziz-Howl loophole) | Classical gravity can produce entanglement → both survive |

### Key: Framework works EITHER WAY

If gravity is quantum: τ > 0 is FAPP, but operationally indistinguishable from genuine in all practical experiments. Paper 4's full unification becomes the natural endpoint.

If gravity is classical: τ > 0 is genuine (Galley-Giacomini-Selby). Paper 1b becomes the foundational result.

**Either way, Papers 1-3 are valid.** The difference is interpretive, not mathematical.

---

## Key References (New in This Round)

| Paper | Key Result |
|-------|-----------|
| Galley-Giacomini-Selby 2023 (Quantum 7, 1142) | Classical gravity + quantum matter = irreversible |
| Oppenheim 2023 (PRX 13, 041040) | Post-quantum classical gravity theory |
| Oppenheim et al. 2023 (Nat. Commun. 14, 7910) | Decoherence-diffusion tradeoff |
| Grossardt et al. 2025 (PRD 111, L121101) | Diósi-Penrose classical gravity can entangle |
| Aziz & Howl 2025 (Nature 646, 813) | Classical gravity produces entanglement (disputed) |
| Marletto et al. 2025 (arXiv:2511.00852) | Rebuttal to Aziz-Howl |
| Donoghue 1994 (PRD 50, 3874) | GR as effective field theory |
| Jacobson 1995 (PRL 75, 1260) | Einstein equations as equation of state |
| Hossenfelder 2025 (arXiv:2510.11037) | Product state constraint → collapse |

---

## Next Steps

1. **Decide Paper 1b scope**: Standalone paper or section in Paper 2?
2. **Adopt operational framing**: "Gravity behaves classically" (safe) vs "gravity IS classical" (risky)
3. **Resolve Paper 4 strategy**: Bianconi formalism vs scale-dependent ontology
4. **Calculate**: τ_grav for specific experimental scenarios (molecule interferometry)
5. **Compare predictions**: Our τ vs Diósi-Penrose vs Oppenheim — where do they differ?
