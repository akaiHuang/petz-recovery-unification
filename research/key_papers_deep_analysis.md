# Deep Analysis: Three Key Papers (2025-2026)

> Compiled: 2026-03-09
> These three papers significantly impact the τ framework and viewpoint paper

---

## 1. Liu-Bai-Scarani 2026 (PRL 136, 060203)

**"Proper and Improper Mixed States Serve as Different Prior Beliefs for Quantum State Retrodiction"**

- arXiv: [2502.10030](https://arxiv.org/abs/2502.10030)
- Journal: Physical Review Letters **136**, 060203 (2026)
- Authors: Mingxuan Liu, Ge Bai, Valerio Scarani (NUS/CQT)

### Core Result

The density matrix alone is NOT sufficient for retrodiction. Two states with **identical** density matrices — one a proper mixture (classical ignorance), the other an improper mixture (entangled subsystem) — give **completely different** retrodictions via the Petz map.

### Key Example (Qubit)

β_S = I/2 in both cases:

| Prior belief | After seeing |0⟩ | After seeing |1⟩ |
|---|---|---|
| **Proper** (|0⟩ or |1⟩ with 50/50) | Retrodict |0⟩ with certainty | Retrodict |1⟩ with certainty |
| **Improper** (half of Bell state) | Retrodict I/2 (no update) | Retrodict I/2 (no update) |

### Why This Matters for τ

**τ = 0 has TWO completely different meanings:**
1. **Closed system (unitary)**: F = 1, Σ = 0. Perfect recovery. Time genuinely absent.
2. **Improper mixture (entangled subsystem)**: "Recovery" is trivial — belief never changes. Not real reversibility, but **absence of a definite past**.

**The AND→OR transition = improper→proper transition:**
- AND regime (improper): No retrodiction possible → no definite past → superposition persists
- OR regime (proper): Retrodiction works → definite past exists → classical behavior

**Decoherence IS genuinely irreversible at the retrodictive level:** Converting improper → proper mixture changes the retrodiction structure. This is not FAPP — it's a qualitative change in what inference is possible.

### Key Equation

Prior-extended Petz map (Eq. 2):
```
R^{E,β}_ext(σ) := Tr_R[ √β · (E†(E(β_S)^{-1/2} σ E(β_S)^{-1/2}) ⊗ 1_R) · √β ]
```

Equivalence condition (Theorem 1): Two beliefs β, γ give same retrodiction iff their **second moments** match:
```
Tr_{RR'}[|√β⟩⟩⟨⟨√β|] = Tr_{RR'}[|√γ⟩⟩⟨⟨√γ|]
```

### Impact on Our Work

- Paper 1 Fix 1 (state-specific vs universal τ=0) is now MORE important
- Viewpoint paper Section 3 should incorporate this distinction
- Potentially the most important paper for the AND→OR discussion
- Related: Pinske-Mølmer (arXiv:2504.16170) — retrodiction witnesses entanglement

---

## 2. Jean-Silva-Vilasini 2025 (arXiv:2508.02463)

**"An equivalence between time-symmetry and cyclic causality in quantum theory"**

- arXiv: [2508.02463](https://arxiv.org/abs/2508.02463)
- Authors: Eliot Jean, Ralph Silva, V. Vilasini (ETH Zürich, Grenoble)
- Submitted: August 2025, 27+21 pages

### Core Result

**Bidirectional operational equivalence:**

> Multi-time states (MTS, time-symmetric QM) ⟺ P-CTC-assisted combs (cyclic causality)

Every time-symmetric quantum process can be replicated with post-selected closed timelike curves (Lloyd P-CTCs), and vice versa. The mapping is constructive and preserves all measurement statistics.

### Connection to τ = 0

The logical chain:
```
τ = 0  ⟺  time-symmetric  ⟺  P-CTC-equivalent  (by Jean et al.)
```

A closed system (τ = 0) is **operationally equivalent to having access to P-CTCs**. This is NOT saying CTCs physically exist — it's saying the operational predictions are identical.

### What P-CTCs Are

P-CTCs (Lloyd et al. 2011) = post-selected quantum teleportation. NOT Deutsch CTCs. Key differences:

| Feature | Deutsch CTCs | P-CTCs |
|---------|-------------|--------|
| Mechanism | Self-consistency | Post-selected teleportation |
| Computational power | BQP = PSPACE | BQP = PP |
| Physical basis | Unclear | Standard QM + post-selection |

### Why This Matters for τ

- τ = 0 has a much sharper physical meaning: not just "reversible" but "P-CTC-equivalent"
- The "information time machine" framing is mathematically substantiated: closed systems have bidirectional information flow structure
- The degradation from τ = 0 to τ > 0 is the loss of P-CTC structure
- The quantum eraser (τ = 0 for signal-idler pair) has P-CTC structure; coincidence counting = post-selection
- **The gap Jean et al. leave = what τ fills**: they prove a binary equivalence; τ provides the continuous interpolation with F ≥ exp(−Σ/2) as the degradation bound

### Caveats

1. Post-selection required (exponentially unlikely in practice)
2. Finite-dimensional only
3. Operational, not ontological
4. Specifically P-CTCs, not Deutsch CTCs

---

## 3. Kuang-Torii-Buscemi 2025 (arXiv:2511.20281)

**"Quantum measurement retrodiction and entropic uncertainty relations"**

- arXiv: [2511.20281](https://arxiv.org/abs/2511.20281)
- Authors: Jiaxi Kuang, Kensei Torii, Francesco Buscemi
- Submitted: November 2025, 17 pages

### Core Result

**For measurement channels (quantum-to-classical), ALL standard quantum divergences give the SAME retrodiction — the Petz map.**

This holds for: Umegaki RE, Petz-Rényi, sandwiched Rényi, geometric Rényi, Belavkin-Staszewski RE, and trace distance.

The unique optimizer is:
```
γ̂_{Q|x} = √(γ_Q) M_x √(γ_Q) / Tr{M_x γ_Q}
```
which is exactly the Petz transpose output.

### Why Uniqueness Holds for Measurements

The semiclassical (direct-sum) structure of measurement channels allows the divergence to decompose:
```
D(σ_{XQ} ‖ γ_{XQ}) = D(τ ‖ p) + Σ_x p(x) D(σ_{Q|x} ‖ γ̂_{Q|x})
```
Each conditional term is independently minimized. This decomposition works for ALL well-behaved divergences.

### When Uniqueness BREAKS

For quantum-to-quantum channels (output retains coherence), different divergences give different retrodictions. **Quantum coherence in the output makes retrodiction inherently ambiguous.**

### Why This Matters for τ

1. **τ = 1 − F_Petz is THE canonical irreversibility measure for measurements** — not one choice among many, but the unique one
2. **The equivalence chain is divergence-independent for measurements** — the bound F ≥ exp(−Σ/2) is THE bound
3. **The AND→OR transition (measurement) creates unambiguous retrodiction** — before measurement, retrodiction is ambiguous; after, it's unique
4. **Strengthens Paper 1**: should add a remark noting this divergence-independence

### Key Limitation

For general quantum channels, τ remains one distinguished choice (via Parzygnat-Buscemi axioms, DPI saturation), but not the unique minimum-change retrodiction.

---

## Synthesis: How These Three Papers Change Our Strategy

### The Emerging Picture

```
                    τ = 0                          τ → 1
                 (closed system)               (measurement)
                      │                             │
    Jean et al.  ─→  P-CTC equivalent    Kuang et al.  ─→  unique retrodiction
    Liu et al.   ─→  improper mixture     Liu et al.    ─→  proper mixture
                 ─→  no definite past                   ─→  definite past
                 ─→  AND regime                         ─→  OR regime
                 ─→  time absent                        ─→  time present
```

### What This Means for the Viewpoint Paper

1. **Section 3** (coexistence): Use Jean et al. to strengthen "τ=0 = time-symmetric = P-CTC-equivalent"
2. **Section 3** (AND→OR): Use Liu et al. to explain the improper→proper transition as the emergence of retrodictable history
3. **Section 5** (predictions): Use Kuang et al. to argue τ is THE canonical measure for measurement irreversibility
4. **Section 4** (observation paradox): The observer opening the system = converting improper to proper = enabling retrodiction but destroying time-symmetry

### The Unified Story

> In a closed system (τ=0), the mixture is improper, retrodiction is frozen (Liu et al.), and the system has P-CTC structure (Jean et al.). There is no definite past, no definite future — they coexist as equivalent descriptions.
>
> When the system opens (τ>0), the mixture becomes proper (Liu et al.), retrodiction becomes possible and unique (Kuang et al.), and the P-CTC structure degrades (Jean et al.). A definite past emerges. This IS the arrow of time.
>
> The cost of this transition is bounded by thermodynamics: F ≥ exp(−Σ/2). The continuous parameter τ interpolates between the two extremes.
