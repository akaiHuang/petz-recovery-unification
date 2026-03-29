# tau-chrono Real Hardware Validation Plan

**Goal**: Demonstrate that Bayesian noise tracking (tau-chrono) gives measurably better noise estimates than naive independent-gate models on real quantum hardware, producing publication-quality data for commercialization.

---

## Experiment Overview

### Available Hardware

| Backend | Qubits | 1Q Fidelity | 2Q (CZ) Fidelity | Connectivity |
|---------|--------|-------------|-------------------|-------------|
| **Starmon-7** | 7 | ~99.9% | ~95.5% | double-rhombus, 8 CZ pairs |
| **Tuna-5** | 5 | TBD | TBD | starfish (4 couplers to central) |
| **QX emulator** | 26 | perfect | perfect | all-to-all |

**Primary target**: Starmon-7 (most qubits, established specs)

---

## Experiment 1: Single-Gate Characterization

**Purpose**: Extract real Kraus operators for each native gate. Establish baseline noise levels.

### Gates to characterize

| Gate | Qiskit | Expected tau (from specs) | Notes |
|------|--------|--------------------------|-------|
| X | `qc.x(q)` | ~0.001 (1Q err ~0.1%) | Pauli-X |
| H | `qc.h(q)` | ~0.001 | Hadamard |
| SX (√X) | `qc.sx(q)` | ~0.001 | Native gate on many platforms |
| Rz(π/4) | `qc.rz(π/4, q)` | ~0 (virtual gate) | Often zero-error on hardware |
| S | `qc.s(q)` | ~0 (virtual) | Phase gate |
| CZ | `qc.cz(q0,q1)` | ~0.045 (2Q err ~4.5%) | **Highest noise** |

### Protocol
- 12 circuits per gate (4 inputs × 3 measurement bases)
- 4096 shots per circuit
- Repeat on qubits 0-6 to check qubit-to-qubit variation

### Expected output: Table 1

```
┌──────────┬────────┬──────────┬──────────┬──────────┬──────────┐
│ Gate     │ Qubit  │ Kraus #  │ τ_worst  │ τ_avg    │ F_proc   │
├──────────┼────────┼──────────┼──────────┼──────────┼──────────┤
│ X        │ Q0     │ 2        │ 0.0012   │ 0.0008   │ 0.9990   │
│ X        │ Q3     │ 2        │ 0.0018   │ 0.0011   │ 0.9985   │
│ H        │ Q0     │ 2        │ 0.0015   │ 0.0010   │ 0.9988   │
│ CZ       │ Q0-Q3  │ 4        │ 0.0450   │ 0.0320   │ 0.9550   │
│ Rz(π/4)  │ Q0     │ 1        │ ~0       │ ~0       │ ~1.0     │
│ ...      │ ...    │ ...      │ ...      │ ...      │ ...      │
└──────────┴────────┴──────────┴──────────┴──────────┴──────────┘
```

### Key prediction
- **1Q gates**: τ ~ 0.001, Kraus rank 2 (identity + small noise)
- **CZ gate**: τ ~ 0.03-0.05, Kraus rank 3-4 (significant noise)
- **Virtual gates (Rz, S)**: τ ~ 0 (no physical gate applied)
- **Qubit variation**: Outer qubits (Q0, Q6) may have higher τ than center (Q3)

---

## Experiment 2: Circuit-Level Bayesian vs Naive Comparison

**Purpose**: The core demo. Show that Bayesian tracking gives lower total τ.

### Circuits to test

| Circuit | Gates | Depth | Why |
|---------|-------|-------|-----|
| **A. H-Rz-X-H** | H, Rz(π/4), X, H | 4 | Simple 1Q, pure noise accumulation |
| **B. 5-gate random** | H, SX, Rz, X, H | 5 | Standard benchmark |
| **C. 10-gate deep** | 10 random 1Q gates | 10 | Stress test: deep circuit advantage |
| **D. 20-gate deep** | 20 random 1Q gates | 20 | Maximum depth test |
| **E. Bell circuit** | H, CZ | 2 (2Q) | 2-qubit with noisy CZ |
| **F. GHZ-3** | H, CZ, CZ | 3 (3Q) | 3-qubit scaling test |

### Protocol
- For each circuit: characterize every unique gate via tomography
- Run Bayesian composition with real Kraus operators
- Compare τ_Bayesian vs τ_naive (multiplicative independent model)
- Also run the actual circuit and measure output fidelity for ground-truth

### Expected output: Table 2

```
┌──────────────┬────────┬──────────────┬──────────────┬─────────────┐
│ Circuit      │ Depth  │ τ_naive      │ τ_Bayesian   │ Improvement │
├──────────────┼────────┼──────────────┼──────────────┼─────────────┤
│ H-Rz-X-H    │ 4      │ 0.0040       │ 0.0032       │ 20%         │
│ 5-gate       │ 5      │ 0.0050       │ 0.0038       │ 24%         │
│ 10-gate      │ 10     │ 0.0100       │ 0.0065       │ 35%         │
│ 20-gate      │ 20     │ 0.0199       │ 0.0110       │ 45%         │
│ Bell (2Q)    │ 2      │ 0.0460       │ 0.0350       │ 24%         │
│ GHZ-3 (3Q)  │ 3      │ 0.0880       │ 0.0590       │ 33%         │
└──────────────┴────────┴──────────────┴──────────────┴─────────────┘
```

### Key predictions
1. **Improvement scales with depth**: 4 gates ~20%, 10 gates ~35%, 20 gates ~45%
2. **2Q circuits show larger absolute τ** (CZ noise ~30× larger than 1Q)
3. **Bayesian advantage is larger when CZ gates are present** (more noise to track)
4. **Composition inequality √τ_total ≤ Σ√τ_i^eff should HOLD for all circuits**

---

## Experiment 3: Depth Scaling Curve

**Purpose**: Generate the "money plot" — improvement % vs circuit depth.

### Protocol
- Fix gate set: random 1Q gates (H, X, SX, Rz)
- Vary depth: 2, 4, 6, 8, 10, 15, 20, 30
- At each depth: compute τ_naive and τ_Bayesian
- Plot improvement % vs depth

### Expected output: Figure 1 (Depth Scaling)

```
Improvement %
    50% ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ┐
    45% │                                     ●───●
    40% │                               ●───●
    35% │                         ●───●
    30% │                   ●───●
    25% │             ●───●
    20% │       ●───●
    15% │  ●───●
    10% │ ●
     5% │●
     0% └─────┬─────┬─────┬─────┬─────┬─────┬─────
              2     5     8    10    15    20    30
                        Circuit Depth
```

**Commercial message**: "The deeper your circuit, the more you benefit from tau-chrono."

---

## Experiment 4: Qubit-to-Qubit Variation Map

**Purpose**: Show tau-chrono captures hardware non-uniformity.

### Protocol
- Run same gate (H) on all 7 qubits
- Run CZ on all 8 available pairs
- Build a "noise heatmap" of the chip

### Expected output: Figure 2 (Chip Noise Map)

```
Starmon-7 Noise Map (τ per gate)

        Q0 ─── Q1
       / \       \
      /   \       \
    Q2 ── Q3 ── Q4
           |       \
           |        \
          Q5 ── Q6

Node color: τ_H (1Q noise)    0.0005 ← green ─── red → 0.003
Edge color: τ_CZ (2Q noise)   0.02 ← green ─── red → 0.06
```

**Commercial message**: "tau-chrono tells you which qubits and connections to trust."

---

## Experiment 5: Ground-Truth Validation

**Purpose**: Verify that τ_Bayesian actually predicts real circuit fidelity better than τ_naive.

### Protocol
1. Prepare known input state |ψ_in⟩ = |0⟩
2. Run circuit on hardware, measure output ρ_out
3. Do state tomography to get full ρ_out
4. Compute actual fidelity F_actual = F(|ψ_ideal⟩, ρ_out)
5. Compare:
   - F_predicted_naive = 1 - τ_naive
   - F_predicted_Bayesian = 1 - τ_Bayesian
   - F_actual (from state tomography)

### Expected output: Table 3

```
┌──────────────┬──────────────┬──────────────┬──────────────┐
│ Circuit      │ F_actual     │ F_naive_pred │ F_Bayes_pred │
├──────────────┼──────────────┼──────────────┼──────────────┤
│ H-Rz-X-H    │ 0.9965       │ 0.9960       │ 0.9968       │
│ 10-gate      │ 0.9910       │ 0.9900       │ 0.9935       │
│ 20-gate      │ 0.9820       │ 0.9801       │ 0.9890       │
│ Bell (2Q)    │ 0.9570       │ 0.9540       │ 0.9650       │
│ GHZ-3 (3Q)  │ 0.9180       │ 0.9120       │ 0.9410       │
└──────────────┴──────────────┴──────────────┴──────────────┘
```

### Key prediction
- **F_Bayes_pred will be closer to F_actual than F_naive_pred**
- This is the "ground truth" that proves Bayesian tracking is not just mathematically better but **physically accurate**
- Residual |F_actual - F_Bayes_pred| should be smaller than |F_actual - F_naive_pred|

---

## Figures for Publication / Commercialization

### Figure 1: Depth Scaling Curve
- X-axis: circuit depth
- Y-axis: Bayesian improvement %
- Shows monotonic increase → "deeper circuits benefit more"

### Figure 2: Chip Noise Map
- Starmon-7 topology with colored nodes/edges
- Visual at a glance which qubits/connections are noisiest

### Figure 3: Prediction Accuracy
- Bar chart: |F_actual - F_predicted| for naive vs Bayesian
- Bayesian bars should be consistently shorter

### Figure 4: Per-Gate τ Breakdown
- Stacked bar chart showing τ_eff per gate in a circuit
- Shows how Bayesian tracking "compresses" later gates' τ toward 0

### Figure 5: Bayesian vs Naive Scatter
- X-axis: τ_naive for each circuit
- Y-axis: τ_Bayesian for each circuit
- All points below the y=x line → Bayesian always wins

---

## Resource Estimate

| Item | Count | Shots | Total circuits |
|------|-------|-------|---------------|
| Exp 1: 6 gates × 7 qubits | 42 tomographies | 42 × 12 × 4096 | 504 circuits |
| Exp 1: CZ × 8 pairs (2Q tomo) | 8 tomographies | 8 × 144 × 4096 | 1,152 circuits |
| Exp 2: 6 circuits | ~30 gate tomos | 30 × 12 × 4096 | 360 circuits |
| Exp 3: 8 depths | ~40 gate tomos | 40 × 12 × 4096 | 480 circuits |
| Exp 5: 5 state tomos | 5 × 3-9 circuits | 5 × 9 × 4096 | 45 circuits |
| **Total** | | | **~2,541 circuits** |

At 4096 shots each: ~10.4M total shots.
On QI free tier: should be feasible (no per-shot cost on Starmon-7).

### Time estimate
- QI emulator: ~2 minutes total
- Starmon-7: depends on queue, but ~30-60 minutes of QPU time

---

## Priority Order

1. **Exp 2 (Circuit comparison)** — fastest to run, gives the headline number
2. **Exp 3 (Depth scaling)** — the "money plot"
3. **Exp 5 (Ground truth)** — validates the approach is physically correct
4. **Exp 1 (Gate characterization)** — fills out the full picture
5. **Exp 4 (Chip map)** — visual wow factor

---

## Success Criteria

| Criterion | Threshold | Target |
|-----------|-----------|--------|
| Bayesian improvement (5-gate 1Q) | > 10% | > 20% |
| Bayesian improvement (10-gate 1Q) | > 20% | > 35% |
| Bayesian improvement (2Q circuit) | > 15% | > 25% |
| F_Bayes closer to F_actual than F_naive | All circuits | All circuits |
| Composition inequality holds | All circuits | All circuits |
| CPTP error of reconstructed channels | < 0.05 | < 0.02 |

**If all criteria pass**: tau-chrono is validated for real hardware. Ready for:
- Technical paper (PRA/PRApplied)
- Demo at quantum computing conferences
- Commercial pilot with quantum hardware vendors

---

## What Would Make This Publishable?

A paper titled something like: "Bayesian Noise Tracking Improves Quantum Circuit Error Estimation by 20-45%: Experimental Validation on a 7-Qubit Superconducting Processor"

Key claims:
1. We extract full Kraus operators via process tomography on real hardware
2. Bayesian state tracking through circuits gives 20-45% lower error estimates
3. The improvement scales with circuit depth (verified up to 20-30 gates)
4. Predictions from Bayesian model are closer to ground-truth fidelity
5. Based on Petz recovery theory (Paper 1), with rigorous composition theorem
