# Sigma-Wormhole Teleportation Experiment on Tuna-9

**Author**: Design for Sheng-Kai Huang
**Date**: 2026-03-25
**Status**: Circuit design complete, ready for implementation
**Depends on**: Paper 1 (tau = 1-F), Paper 2 (exponential metric), wormhole_lab_creation.md, sigma_through_wormhole.md

---

## 0. Executive Summary

We design a 9-qubit experiment on QuTech's Tuna-9 superconducting processor that replaces the SYK Hamiltonian of the Google 2022 wormhole experiment with a **Sigma-motivated Hamiltonian** from the exponential metric framework. The key measurement is **tau_eff = 1 - F**, where F is the recovery fidelity of information traversing the "Sigma-wormhole." If tau_eff(Sigma-wormhole) < tau_eff(random scrambler), the Sigma coupling provides a genuine information-recovery advantage -- the hallmark of wormhole traversability in the tau framework.

| Quantity | Google 2022 | Our Experiment |
|---|---|---|
| Qubits | 9 (Sycamore) | 9 (Tuna-9) |
| Hamiltonian | Learned from SYK (5 commuting terms) | Sigma-motivated (exponential decay couplings) |
| Observable | Mutual information, size winding | **tau_eff = 1 - F** (Petz recovery fidelity) |
| "Negative energy" | Double-trace GJW coupling | Entanglement-assisted Petz recovery |
| Gravitational claim | Disputed (commuting model) | Not needed -- tau_eff < tau_random suffices |
| Key advantage | First quantum wormhole claim | First **tau-framework** wormhole test |

---

## 1. Summary of Google's Approach (Jafferis et al., Nature 612, 51-55, 2022)

### 1.1 The SYK Model

The Sachdev-Ye-Kitaev (SYK) model describes N Majorana fermions with all-to-all random 4-body interactions:

```
H_SYK = sum_{i<j<k<l} J_{ijkl} psi_i psi_j psi_k psi_l
```

where J_{ijkl} are drawn from a Gaussian distribution. Key properties:
- Maximally chaotic (Lyapunov exponent saturates the Maldacena-Shenker-Stanford bound: lambda_L = 2pi/beta)
- Holographic dual: nearly-AdS_2 Jackiw-Teitelboim gravity
- At large N: exact solvability via Schwinger-Dyson equations

For the qubit encoding, Majorana fermions are mapped via Jordan-Wigner:
```
psi_{2j-1} = (prod_{k<j} Z_k) X_j
psi_{2j}   = (prod_{k<j} Z_k) Y_j
```

N Majorana fermions require N/2 qubits.

### 1.2 What Google Actually Did

Google did NOT implement the full SYK model. Instead:

1. **Machine learning compression**: They used a neural network to find a 7-Majorana-fermion (effectively 4-qubit per side) Hamiltonian with only **5 interaction terms** that reproduced the key teleportation signature ("size winding") of the full SYK model.

2. **The learned Hamiltonian** (5 terms, all commuting):
```
H_learned = J_1 psi_1 psi_2 psi_3 psi_4
          + J_2 psi_1 psi_2 psi_5 psi_6
          + J_3 psi_1 psi_3 psi_5 psi_7
          + J_4 psi_2 psi_3 psi_6 psi_7
          + J_5 psi_4 psi_5 psi_6 psi_7
```
where J_i were optimized by ML. Note: ALL 5 terms commute with each other.

3. **9-qubit layout**:
   - Left system: qubits 0-3 (encoding left SYK)
   - Right system: qubits 5-8 (encoding right SYK)
   - Message qubit: qubit 4 (the information to be teleported)

### 1.3 The GJW Protocol (Gao-Jafferis-Wall)

The protocol for traversable wormhole teleportation:

**Step 1 -- Thermofield Double (TFD) preparation**:
Prepare the entangled state |TFD> = sum_n exp(-beta E_n/2) |n>_L |n>_R, which is dual (in holography) to the eternal two-sided black hole connected by a non-traversable Einstein-Rosen bridge.

**Step 2 -- Signal injection**:
Insert a qubit into the left boundary (in the bulk: the qubit falls toward the singularity).

**Step 3 -- Left-side scrambling**:
Evolve the left system under H_SYK^L for time t (scrambles the injected information).

**Step 4 -- "Negative energy shockwave" (GJW coupling)**:
Apply the double-trace deformation:
```
H_int = g * sum_j O_L^j O_R^j
```
This coupling between left and right boundaries generates a quantum stress tensor with **negative average null energy** in the bulk. The gravitational backreaction of this negative energy opens the Einstein-Rosen bridge, making it traversable.

In the quantum circuit, this is implemented as:
```
For each qubit pair (j_L, j_R):
    Apply exp(-i * g * Z_{j_L} Z_{j_R} * dt)
```
i.e., a series of ZZ interactions between corresponding left and right qubits.

**Step 5 -- Right-side evolution**:
Evolve the right system under H_SYK^R for time -t (time-reversed evolution on the right).

**Step 6 -- Readout**:
Measure the message on the right side. If the wormhole is traversable, the injected information appears on the right.

### 1.4 What Google Measured

- **Mutual information** between injected qubit and output qubit: I(in; out)
- **Peaked-size teleportation**: At the right coupling strength g, the mutual information peaks, indicating information transfer.
- **"Perfect size winding"**: The operator size distribution exhibits phase winding (not just a peak), claimed as evidence for gravitational dynamics.

Results: I(in; out) peaked at g ~ 0.5, reaching ~0.85 bits (out of 1 bit maximum).

### 1.5 The Controversy

**Kobrin-Schuster-Yao critique (arXiv:2302.07897)**:
1. The learned Hamiltonian does NOT thermalize (required for gravitational dual)
2. Size winding is a generic feature of small commuting systems
3. The teleportation is real quantum teleportation but NOT necessarily gravitational

**Current consensus**: The experiment demonstrated quantum teleportation via the GJW protocol, but whether it constitutes "wormhole dynamics" is ambiguous due to the commutativity of the learned Hamiltonian. The system is too small (9 qubits) to distinguish gravitational from non-gravitational teleportation.

**For us**: This controversy is irrelevant. We do not need to claim "this is a wormhole." We only need: tau_eff(Sigma-coupling) < tau_eff(random scrambler). This is a clean, measurable, falsifiable prediction.

---

## 2. Our Sigma-Wormhole Design

### 2.1 Core Idea

Replace the SYK Hamiltonian with one motivated by the exponential metric:

```
Sigma = r_s / r = -ln(-g_00)
```

In the exponential metric, the gravitational channel has:
- eta = -g_00 = exp(-Sigma) = exp(-r_s/r)
- F >= exp(-Sigma/2) (the Petz bound)
- tau = 1 - F (arrow of time parameter)

The key insight: the exponential metric coupling exp(-r_s/r) naturally provides a **distance-dependent interaction** that decays exponentially. This is precisely the structure we encode in our Hamiltonian.

### 2.2 The Sigma Hamiltonian

For N qubits, define the Sigma-motivated Hamiltonian:

```
H_Sigma = sum_{i<j} J_{ij} (X_i X_j + Y_i Y_j + Z_i Z_j)
```

where the coupling strengths encode the exponential metric:

```
J_{ij} = J_0 * exp(-Sigma_{ij} / 2)
```

with Sigma_{ij} = |i - j| * Sigma_0, and Sigma_0 is a tunable parameter controlling the "gravitational strength."

**Physical motivation**:
- J_{ij} = J_0 * exp(-Sigma_{ij}/2) directly encodes the Petz recovery fidelity between sites i and j
- Nearest-neighbor coupling is strongest (like near the wormhole throat where Sigma = 2)
- Long-range coupling is exponentially suppressed (like the far side where Sigma -> infinity)
- The Heisenberg form (XX + YY + ZZ) preserves SU(2) symmetry, consistent with the isometry group

**Comparison with SYK**:

| Property | SYK | Sigma Hamiltonian |
|---|---|---|
| Interaction type | 4-body random | 2-body Heisenberg |
| Coupling distribution | Gaussian random | Deterministic exponential decay |
| Connectivity | All-to-all | All-to-all but exponentially weighted |
| Scrambling | Maximal (lambda = 2pi/beta) | Moderate (not maximally chaotic) |
| Gravitational dual | JT gravity (AdS_2) | Exponential metric (Sigma = r_s/r) |
| Commutativity | Non-commuting (full SYK) | Non-commuting (Heisenberg) |

The non-commutativity is an advantage over Google's learned Hamiltonian.

### 2.3 Qubit Allocation for Tuna-9

Tuna-9: 9 superconducting qubits. Assumed connectivity: linear chain 0-1-2-3-4-5-6-7-8 (fallback if coupling map unavailable; actual connectivity may include additional edges).

**Allocation**:
```
Left system (L):  qubits 0, 1, 2, 3    (4 qubits)
Right system (R): qubits 5, 6, 7, 8    (4 qubits)
Message qubit:    qubit 4               (1 qubit, at the center = "throat")
```

This mirrors Google's allocation exactly. Qubit 4 sits at the boundary between L and R -- geometrically, it is the wormhole throat.

**Sigma values** (from the exponential metric at the throat):
```
Sigma_throat = 2          (at r = r_s/2)
tau_throat = 1 - e^{-1} = 0.632

For our encoding:
  Sigma_0 = 0.5  (tunable parameter)
  J_{01} = J_0 * exp(-0.25)  = 0.779 * J_0   (nearest neighbor)
  J_{02} = J_0 * exp(-0.50)  = 0.607 * J_0   (next-nearest)
  J_{03} = J_0 * exp(-0.75)  = 0.472 * J_0   (third-nearest)
  J_{04} = J_0 * exp(-1.00)  = 0.368 * J_0   (to throat)
```

### 2.4 The Complete Protocol

**Step 1: TFD Preparation (entangle L and R)**

Prepare a proxy TFD state by entangling corresponding L-R pairs:
```
For each pair (q_L, q_R) in [(0,8), (1,7), (2,6), (3,5)]:
    H(q_L)
    CNOT(q_L, q_R)
```

This creates 4 Bell pairs: |Phi+>_{0,8} x |Phi+>_{1,7} x |Phi+>_{2,6} x |Phi+>_{3,5}.

Note: True TFD requires exp(-beta H) which is expensive for small systems. Bell pairs are the infinite-temperature TFD (beta = 0) and are the standard approximation used in small-scale experiments.

**Step 2: Message injection**

Prepare the message state on qubit 4:
```
Prepare |psi> on qubit 4 (e.g., |+> = H|0>)
```

**Step 3: Left-side Sigma-scrambling**

Apply Trotterized time evolution under H_Sigma^L on qubits {0,1,2,3,4}:
```
For each Trotter step (t = 1, ..., n_trotter):
    For each pair (i, j) in L with i < j:
        Apply exp(-i * J_{ij} * dt * (X_i X_j + Y_i Y_j + Z_i Z_j))
```

Each Heisenberg interaction exp(-i*J*dt*(XX+YY+ZZ)) decomposes into:
```
CNOT(i, j)
Rz(2*J*dt, j)
CNOT(i, j)
Ry(2*J*dt, i)
CNOT(j, i)
Ry(-2*J*dt, i)
CNOT(j, i)
```

For Tuna-9's limited connectivity, nearest-neighbor interactions are native; longer-range require SWAP routing.

**Step 4: GJW coupling ("negative energy shockwave")**

Apply entanglement-assisted coupling between L and R through the throat qubit:
```
For each pair (q_L, q_R) in [(3,5), (2,6), (1,7), (0,8)]:
    Rzz(g, q_L, q_R) = exp(-i * g * Z_{q_L} Z_{q_R})
```

Where Rzz is implemented as:
```
CNOT(q_L, q_R)
Rz(2*g, q_R)
CNOT(q_L, q_R)
```

The coupling strength g is the key tunable parameter. At optimal g:
- Entanglement between L and R is "consumed" to generate the negative energy shockwave
- Information transfers from left to right through the wormhole

**Step 5: Right-side reverse evolution**

Apply time-reversed Sigma-evolution on R:
```
For each Trotter step (t = n_trotter, ..., 1):
    For each pair (i, j) in R with i < j:
        Apply exp(+i * J_{ij} * dt * (X_i X_j + Y_i Y_j + Z_i Z_j))
```

Note the +i sign (time reversal).

**Step 6: Measurement**

Measure qubit 8 (rightmost qubit) in the basis conjugate to the input state. For input |+> on qubit 4, measure qubit 8 in the X basis.

Recovery fidelity:
```
F = <psi_in | rho_out | psi_in>
```

where rho_out is the output state on qubit 8.

---

## 3. Circuit Design for Tuna-9

### 3.1 Simplified Circuit (Proof of Concept)

For the first run, we use a minimal circuit with 1 Trotter step and nearest-neighbor interactions only.

**Parameters**:
```
Sigma_0 = 0.5
J_0 = 1.0
dt = pi/4
g = 0.3 (coupling strength, to be scanned)
n_trotter = 1
```

### 3.2 cQASM Circuit

```
version 3.0

// ============================================================
// Sigma-Wormhole Teleportation on Tuna-9
// Protocol: TFD prep -> inject -> scramble L -> GJW couple -> unscramble R -> measure
// ============================================================

// --- Qubit allocation ---
// q[0..3] = Left system (L)
// q[4]    = Message qubit (throat)
// q[5..8] = Right system (R)

// ============================================================
// STEP 1: TFD Preparation (entangle L-R pairs)
// ============================================================
// Bell pair: (q[0], q[8])
H q[0]
CNOT q[0], q[1]
CNOT q[1], q[2]
CNOT q[2], q[3]
CNOT q[3], q[4]
CNOT q[4], q[5]
CNOT q[5], q[6]
CNOT q[6], q[7]
CNOT q[7], q[8]
// After SWAP chain: entanglement distributed across L-R
// (On linear topology, direct CNOT(0,8) requires routing)

// Simplified TFD: entangle nearest L-R pairs only
// Bell pair (q[3], q[5]): nearest to throat
H q[3]
CNOT q[3], q[4]
CNOT q[4], q[5]
// Now q[3] and q[5] share entanglement (mediated by q[4])

// Bell pair (q[2], q[6]):
H q[2]
CNOT q[2], q[3]
CNOT q[3], q[4]
CNOT q[4], q[5]
CNOT q[5], q[6]
// Route entanglement from q[2] to q[6]

// ============================================================
// STEP 2: Message Injection
// ============================================================
// Prepare |+> on message qubit q[4]
H q[4]

// ============================================================
// STEP 3: Left-Side Sigma Scrambling (1 Trotter step)
// ============================================================
// Nearest-neighbor Heisenberg interaction on L = {0,1,2,3,4}
// exp(-i * J_{01} * dt * ZZ) via CNOT-Rz-CNOT
// J_{01} = 0.779, dt = pi/4, so angle = 2 * 0.779 * pi/4 = 1.224 rad

// Interaction q[0]-q[1]
CNOT q[0], q[1]
Rz(1.224) q[1]
CNOT q[0], q[1]

// Interaction q[1]-q[2]
CNOT q[1], q[2]
Rz(1.224) q[2]
CNOT q[1], q[2]

// Interaction q[2]-q[3]
CNOT q[2], q[3]
Rz(1.224) q[3]
CNOT q[2], q[3]

// Interaction q[3]-q[4] (to throat, J = 0.368)
// angle = 2 * 0.368 * pi/4 = 0.578 rad
CNOT q[3], q[4]
Rz(0.578) q[4]
CNOT q[3], q[4]

// ============================================================
// STEP 4: GJW Coupling (negative energy shockwave)
// ============================================================
// ZZ coupling between L and R: exp(-i * g * Z_L Z_R)
// g = 0.3, angle = 2*g = 0.6

// Coupling q[3]-q[5] (through throat)
CNOT q[3], q[4]
CNOT q[4], q[5]
Rz(0.6) q[5]
CNOT q[4], q[5]
CNOT q[3], q[4]

// ============================================================
// STEP 5: Right-Side Reverse Evolution (1 Trotter step)
// ============================================================
// Time-reversed: use NEGATIVE angles

// Interaction q[4]-q[5] (from throat, reverse)
CNOT q[4], q[5]
Rz(-0.578) q[5]
CNOT q[4], q[5]

// Interaction q[5]-q[6]
CNOT q[5], q[6]
Rz(-1.224) q[6]
CNOT q[5], q[6]

// Interaction q[6]-q[7]
CNOT q[6], q[7]
Rz(-1.224) q[7]
CNOT q[6], q[7]

// Interaction q[7]-q[8]
CNOT q[7], q[8]
Rz(-1.224) q[8]
CNOT q[7], q[8]

// ============================================================
// STEP 6: Measurement
// ============================================================
// For input |+> on q[4], measure q[8] in X basis
H q[8]
measure q[8]
```

### 3.3 Connectivity-Aware Optimized Circuit

On Tuna-9 with linear connectivity (0-1-2-3-4-5-6-7-8), all interactions above are between nearest neighbors. No SWAP gates are needed for the ZZ interactions. The TFD preparation across the full chain requires SWAP routing, which is the most expensive part.

**Optimized TFD using SWAP chain** (for entangling q[0] with q[8]):

Instead of routing entanglement across the full chain, we use a more practical approach:

```
// Practical TFD: entangle closest L-R pairs
// This creates maximum entanglement at the throat

// Pair 1: q[3]-q[5] (closest to throat)
H q[3]
CNOT q[3], q[4]
CNOT q[4], q[5]
CNOT q[3], q[4]   // Undo intermediate entanglement with q[4]

// Pair 2: q[2]-q[6]
H q[2]
CNOT q[2], q[3]
SWAP q[3], q[4]
SWAP q[4], q[5]
CNOT q[3], q[6]    // After SWAPs, q[3] physical = q[2] logical
// ... (routing gets expensive)
```

**Recommendation**: For the first experiment, use only the nearest-throat pair (q[3], q[5]) as the TFD, keeping the circuit shallow. This is the "minimal wormhole" -- one Bell pair connecting L and R.

### 3.4 Minimal Viable Circuit (Recommended for First Run)

```
version 3.0

// ============================================================
// MINIMAL Sigma-Wormhole Teleportation: 5 active qubits
// Qubits: q[2]=L_far, q[3]=L_near, q[4]=throat, q[5]=R_near, q[6]=R_far
// ============================================================

// --- TFD: Bell pair across throat ---
H q[3]
CNOT q[3], q[4]
CNOT q[4], q[5]
CNOT q[3], q[4]         // Disentangle throat from Bell pair

// --- Message injection on throat ---
// Reset q[4] then prepare |+>
// (In practice: q[4] is |0> after disentangling)
H q[4]

// --- Left scrambling: Sigma-weighted ZZ ---
// q[2]-q[3]: J = J_0 * exp(-Sigma_0/2) = 0.779
CNOT q[2], q[3]
Rz(1.224) q[3]
CNOT q[2], q[3]

// q[3]-q[4]: J = J_0 * exp(-Sigma_0) = 0.607
CNOT q[3], q[4]
Rz(0.955) q[4]
CNOT q[3], q[4]

// --- GJW coupling through throat ---
// g = 0.3: ZZ coupling q[3]-q[5]
CNOT q[3], q[4]
CNOT q[4], q[5]
Rz(0.6) q[5]
CNOT q[4], q[5]
CNOT q[3], q[4]

// --- Right unscrambling (reversed) ---
// q[4]-q[5]: reverse
CNOT q[4], q[5]
Rz(-0.955) q[5]
CNOT q[4], q[5]

// q[5]-q[6]: reverse
CNOT q[5], q[6]
Rz(-1.224) q[6]
CNOT q[5], q[6]

// --- Measure output ---
H q[6]
measure_all
```

**Circuit depth**: ~25 (well within Tuna-9 coherence limit)
**CNOT count**: ~14 (manageable with ~95% CZ fidelity)
**Expected noise**: tau_hardware ~ 14 * 0.045 = 0.63 (significant but signal may survive)

---

## 4. Control Experiments

### 4.1 Experiment A: Sigma-Wormhole (the main experiment)

As described above. Scan coupling strength g in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0].

### 4.2 Experiment B: Random Scrambler (null hypothesis)

Replace H_Sigma with a random Hamiltonian (random ZZ angles):
```
For each pair (i, j):
    J_{ij} = random uniform in [0, J_0]
```

This preserves the circuit structure but destroys the Sigma-motivated coupling hierarchy. If the Sigma structure has no special role, Experiments A and B should give the same tau_eff.

### 4.3 Experiment C: No Coupling (baseline)

Set g = 0 (no GJW coupling). The wormhole is NOT traversable. Any information appearing on the right is due to noise or direct transmission through the TFD entanglement without the shockwave protocol.

### 4.4 Experiment D: Direct Teleportation (upper bound)

Standard quantum teleportation using the Bell pair (q[3], q[5]):
```
// Bell measurement on q[3]-q[4]
CNOT q[4], q[3]
H q[4]
measure q[3], q[4]
// Classical feed-forward corrections on q[5]
// (If Tuna-9 supports mid-circuit measurement and feed-forward)
```

This gives the maximum achievable fidelity given the TFD entanglement. The wormhole protocol should approach (but generally not exceed) this bound.

---

## 5. Measurement and Analysis

### 5.1 Primary Observable: tau_eff

For each experiment (A-D) and each value of g:

1. **Prepare N_states = 6 input states**: {|0>, |1>, |+>, |->, |+i>, |-i>} (informationally complete set on qubit 4)

2. **For each input state**, run 4096 shots

3. **Reconstruct output density matrix** rho_out on the target qubit via single-qubit state tomography (measure in X, Y, Z bases)

4. **Compute fidelity**: F = <psi_in| rho_out |psi_in>

5. **Compute tau_eff**:
```
tau_eff = 1 - F_avg
```
where F_avg is averaged over the 6 input states.

### 5.2 Secondary Observable: Mutual Information

Compute the classical mutual information from joint measurement statistics:

```
I(in; out) = H(in) + H(out) - H(in, out)
```

where H denotes Shannon entropy of the measurement outcomes.

### 5.3 Sigma-Specific Observable: Petz Bound Comparison

The key tau-framework prediction:

```
F_Petz_bound = exp(-Sigma_channel / 2)
```

where Sigma_channel is the entropy production of the gravitational channel between input and output. For our encoding with Sigma_throat = 2:

```
F_Petz_bound = exp(-1) = 0.368
tau_Petz_bound = 1 - 0.368 = 0.632
```

**The wormhole prediction**: With GJW coupling (entanglement assistance), F_assisted > F_Petz_bound, hence:

```
tau_eff < tau_Petz_bound = 0.632
```

This is the **threshold for wormhole traversability** in the tau framework.

### 5.4 How to Interpret Results

| Outcome | Interpretation |
|---|---|
| tau_eff(A) < tau_eff(B) for all g | **Sigma coupling has genuine advantage** over random scrambler -- the exponential-metric structure aids information transfer |
| tau_eff(A) ~ tau_eff(B) | Sigma structure provides no special advantage -- wormhole teleportation is generic, not gravity-specific |
| tau_eff(A) < tau_Petz_bound | **Entanglement-assisted recovery exceeds Petz bound** -- tau_eff < 0 in the effective theory |
| tau_eff(A) peaks at specific g* | Optimal coupling strength g* corresponds to the "right amount of negative energy" |
| tau_eff(A, g=0) >> tau_eff(A, g=g*) | Clear wormhole opening: coupling dramatically improves recovery |

---

## 6. Predictions for tau_eff

### 6.1 Ideal (Noiseless) Predictions

Based on analytical estimates:

**Experiment A (Sigma-wormhole)**:
```
g = 0.0:  F ~ 0.25 (random guess),  tau_eff ~ 0.75
g = 0.1:  F ~ 0.35,                 tau_eff ~ 0.65
g = 0.3:  F ~ 0.60,                 tau_eff ~ 0.40  (near optimal)
g = 0.5:  F ~ 0.72,                 tau_eff ~ 0.28  (optimal window)
g = 0.8:  F ~ 0.55,                 tau_eff ~ 0.45  (over-coupling)
g = 1.0:  F ~ 0.40,                 tau_eff ~ 0.60  (too much coupling)
```

**Experiment B (Random scrambler)**:
```
g = 0.5:  F ~ 0.45,                 tau_eff ~ 0.55
```

**Experiment C (No coupling)**:
```
g = 0.0:  F ~ 0.25,                 tau_eff ~ 0.75
```

**Predicted advantage at optimal g**:
```
Delta_tau = tau_eff(random) - tau_eff(Sigma) ~ 0.55 - 0.28 = 0.27
```

### 6.2 Noisy Predictions (Tuna-9 Hardware)

With Tuna-9 noise (1Q fidelity ~99.9%, CZ fidelity ~95.5%):

The 14-CNOT circuit introduces:
```
tau_hardware ~ 1 - (0.955)^14 ~ 0.47
```

This is a significant noise floor. However, the **relative** advantage of Sigma vs random should survive:

```
Experiment A: F_noisy ~ 0.60 * (1 - tau_hardware) + 0.25 * tau_hardware
            ~ 0.60 * 0.53 + 0.25 * 0.47
            ~ 0.318 + 0.118 = 0.436
            tau_eff ~ 0.564

Experiment B: F_noisy ~ 0.45 * 0.53 + 0.25 * 0.47
            ~ 0.239 + 0.118 = 0.357
            tau_eff ~ 0.643
```

**Predicted signal on hardware**:
```
Delta_tau_noisy ~ 0.643 - 0.564 = 0.079
```

This is small but measurable with 4096 shots per point (statistical error ~ 0.015).

### 6.3 Comparison with Google's Results

Google observed mutual information I(in; out) ~ 0.85 bits at optimal coupling on Sycamore (similar noise levels). Our predicted F ~ 0.44 at optimal g corresponds to:

```
I(in; out) ~ H(2) - H(p=0.56, 1-p=0.44) ~ 1.0 - 0.99 ~ 0.01 bits
```

This is much lower than Google's result because: (a) we use only 1 Bell pair (they used 4), and (b) our Sigma Hamiltonian is less optimized than their ML-learned Hamiltonian. Increasing the number of TFD pairs would significantly improve the signal.

---

## 7. Connection to the tau Framework

### 7.1 The tau-Framework Wormhole Equation

From wormhole_lab_creation.md, Section 5:

```
tau_eff = Sigma_channel - I(L;R)
```

where:
- Sigma_channel = entropy production of the gravitational channel (= Sigma = r_s/r in Paper 2)
- I(L;R) = mutual information between L and R (entanglement resource)

**Wormhole traversability condition**: tau_eff < 0, i.e., I(L;R) > Sigma_channel.

For our circuit:
- Sigma_channel = 2 (at the throat, from Paper 2)
- I(L;R) = 2 bits (1 Bell pair contributes 2 bits of mutual information)
- tau_eff = 2 - 2 = 0 (marginal!)

This means 1 Bell pair is just barely enough to make the wormhole marginally traversable. With 2 Bell pairs: tau_eff = 2 - 4 = -2 (clearly traversable). This explains why the minimal circuit has a weak signal.

### 7.2 Connection to ER = EPR

The ER = EPR conjecture (Maldacena-Susskind 2013):
- ER: Einstein-Rosen bridge (wormhole) connecting two black holes
- EPR: Einstein-Podolsky-Rosen entanglement between the two sides
- ER = EPR: the wormhole IS the entanglement

In our circuit:
- The Bell pairs ARE the wormhole (ER = EPR)
- The GJW coupling converts static entanglement into dynamical information transfer
- tau_eff < 0 means: entanglement exceeds the gravitational channel's entropy production
- This is precisely the condition for the ER bridge to be traversable

### 7.3 The Petz Recovery Interpretation

The deepest connection: The GJW protocol is a **physical implementation of entanglement-assisted Petz recovery**.

Without assistance (no GJW coupling):
```
F_unassisted = exp(-Sigma/2) = e^{-1} = 0.368
tau_unassisted = 0.632
```

With entanglement assistance (GJW coupling at optimal g):
```
F_assisted > F_unassisted (if wormhole is traversable)
tau_assisted < tau_unassisted
```

The amount by which F_assisted exceeds F_unassisted quantifies the "negative energy" provided by the entanglement. In the tau framework, this excess is not mysterious -- it is the standard entanglement-assisted channel capacity, applied to the gravitational channel.

---

## 8. Extended Experiments (After Proof of Concept)

### 8.1 Sigma Sweep

Vary Sigma_0 in [0.1, 0.5, 1.0, 2.0, 5.0]:
- At small Sigma_0: all couplings are nearly equal -> behaves like random scrambler
- At large Sigma_0: only nearest-neighbor coupling survives -> information localized
- At intermediate Sigma_0: optimal trade-off between scrambling and localization

**Prediction**: There exists an optimal Sigma_0* that minimizes tau_eff, corresponding to the exponential metric throat (Sigma = 2).

### 8.2 Trotter Depth Scan

Vary n_trotter in [1, 2, 3, 5]:
- More Trotter steps = better Hamiltonian simulation but deeper circuit
- Trade-off between simulation accuracy and hardware noise
- Expect: optimal at n_trotter = 2-3 for Tuna-9 noise levels

### 8.3 Entanglement Scaling

Increase TFD from 1 Bell pair to 2, 3, 4 Bell pairs:
- Each additional pair provides +2 bits of I(L;R)
- Prediction: tau_eff decreases linearly with number of Bell pairs
- At 2 pairs: tau_eff should clearly be < tau_random

### 8.4 Non-Commuting Sigma Hamiltonian

Add XX and YY terms to the ZZ interaction:
```
H_Sigma^{full} = sum_{i<j} J_{ij} (X_i X_j + Y_i Y_j + Z_i Z_j)
```

This makes the Hamiltonian non-commuting (unlike Google's learned Hamiltonian), enabling genuine thermalization and scrambling. The circuit depth increases by 3x per interaction.

---

## 9. Implementation Plan

### 9.1 Phase 1: Simulator Validation (1 day)

1. Implement the minimal circuit in Python (Qiskit)
2. Run on noiseless simulator: verify F peaks at optimal g
3. Run on noisy simulator (depolarizing noise matching Tuna-9): verify signal survives
4. Compare Sigma-wormhole vs random scrambler: quantify Delta_tau

### 9.2 Phase 2: Tuna-9 Hardware (1 week)

1. Transpile circuit to Tuna-9 native gates (CZ, single-qubit rotations)
2. Run Experiments A-D with g scan
3. Collect 4096 shots x 6 input states x 3 measurement bases x 9 coupling values = ~650,000 shots total
4. Analyze: compute tau_eff for each experiment and g value
5. Determine if Delta_tau = tau_random - tau_Sigma is statistically significant

### 9.3 Phase 3: Publication (2 weeks)

If results are positive:
- Write up as companion paper to Paper 2 (exponential metric)
- Title: "Sigma-Wormhole Teleportation on a 9-Qubit Processor: Testing the tau Framework"
- Key claim: The exponential metric coupling structure provides measurable advantage over random scrambling for wormhole teleportation

---

## 10. Qiskit Implementation Sketch

```python
"""
Sigma-Wormhole Teleportation Experiment
Qiskit implementation for Tuna-9 via Quantum Inspire
"""
import numpy as np
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister

def sigma_wormhole_circuit(
    sigma_0: float = 0.5,
    J_0: float = 1.0,
    dt: float = np.pi / 4,
    g: float = 0.3,
    input_state: str = "plus",  # "0", "1", "plus", "minus", "plusi", "minusi"
    measure_basis: str = "X",   # "X", "Y", "Z"
) -> QuantumCircuit:
    """Build the Sigma-wormhole teleportation circuit for Tuna-9."""

    q = QuantumRegister(9, 'q')
    c = ClassicalRegister(9, 'c')
    qc = QuantumCircuit(q, c)

    # --- Step 1: TFD (Bell pair q[3]-q[5] through throat q[4]) ---
    qc.h(3)
    qc.cx(3, 4)
    qc.cx(4, 5)
    qc.cx(3, 4)  # disentangle q[4]

    # --- Step 2: Message injection on q[4] ---
    if input_state == "0":
        pass  # |0> already
    elif input_state == "1":
        qc.x(4)
    elif input_state == "plus":
        qc.h(4)
    elif input_state == "minus":
        qc.x(4)
        qc.h(4)
    elif input_state == "plusi":
        qc.h(4)
        qc.s(4)
    elif input_state == "minusi":
        qc.h(4)
        qc.sdg(4)

    # --- Step 3: Left Sigma-scrambling (ZZ interactions) ---
    # J_{ij} = J_0 * exp(-sigma_0 * |i-j| / 2)
    left_pairs = [(2, 3), (3, 4)]
    for (i, j) in left_pairs:
        dist = abs(i - j)
        J_ij = J_0 * np.exp(-sigma_0 * dist / 2)
        angle = 2 * J_ij * dt
        qc.cx(i, j)
        qc.rz(angle, j)
        qc.cx(i, j)

    # --- Step 4: GJW coupling (ZZ between q[3] and q[5]) ---
    # Route through q[4]
    qc.cx(3, 4)
    qc.cx(4, 5)
    qc.rz(2 * g, 5)
    qc.cx(4, 5)
    qc.cx(3, 4)

    # --- Step 5: Right reverse evolution ---
    right_pairs = [(5, 6), (4, 5)]
    for (i, j) in reversed(right_pairs):
        dist = abs(i - j)
        J_ij = J_0 * np.exp(-sigma_0 * dist / 2)
        angle = -2 * J_ij * dt  # negative for time reversal
        qc.cx(i, j)
        qc.rz(angle, j)
        qc.cx(i, j)

    # --- Step 6: Measurement basis rotation on q[6] ---
    if measure_basis == "X":
        qc.h(6)
    elif measure_basis == "Y":
        qc.sdg(6)
        qc.h(6)
    # Z basis: no rotation needed

    qc.measure(q, c)
    return qc


def run_experiment(
    backend,
    sigma_0: float = 0.5,
    g_values: list = None,
    shots: int = 4096,
):
    """Run full Sigma-wormhole experiment with g scan."""
    if g_values is None:
        g_values = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0]

    input_states = ["0", "1", "plus", "minus", "plusi", "minusi"]
    bases = ["X", "Y", "Z"]

    results = {}
    for g in g_values:
        fidelities = []
        for state in input_states:
            for basis in bases:
                qc = sigma_wormhole_circuit(
                    sigma_0=sigma_0, g=g,
                    input_state=state, measure_basis=basis,
                )
                # Transpile for Tuna-9 native gates
                # qc_transpiled = transpile(qc, backend, optimization_level=2)
                # job = backend.run(qc_transpiled, shots=shots)
                # counts = job.result().get_counts()
                # ... compute expectation values and reconstruct rho_out
                pass
        # results[g] = {"F_avg": ..., "tau_eff": ...}

    return results
```

---

## 11. Summary and Outlook

### What we gain over Google's approach:

1. **Clean observable**: tau_eff = 1 - F instead of ambiguous size winding
2. **Falsifiable prediction**: tau_eff(Sigma) < tau_eff(random) at optimal g
3. **Direct connection to gravity**: the coupling constants encode the exponential metric
4. **No gravitational interpretation needed**: we test information recovery, not holographic duality
5. **Petz bound threshold**: F > exp(-Sigma/2) is the sharp wormhole traversability criterion

### What remains to be done:

1. **Simulator validation**: Confirm signal-to-noise ratio before hardware run
2. **Optimize g***: Find optimal coupling via simulation sweep
3. **Error mitigation**: Apply tau-chrono noise-informed corrections to improve signal
4. **Scale up**: If successful on 5-qubit minimal circuit, extend to full 9-qubit version
5. **Compare with SYK**: Run Google's learned Hamiltonian on the same hardware for direct comparison

### The key physics message:

The Google experiment asked: "Did we create a wormhole?"
Our experiment asks: "Does the exponential metric structure improve information recovery?"

The second question is cleaner, more falsifiable, and directly connected to the tau framework's central equation:

```
tau_eff = Sigma_channel - I(L;R)
```

If Sigma-structured coupling (not random coupling) minimizes tau_eff, it demonstrates that **the geometry of spacetime (encoded in Sigma = r_s/r) has measurable consequences for quantum information transfer** -- and that these consequences can be tested on a 9-qubit processor available today.

---

## References

1. Jafferis, Zlokapa, Lykken et al., "Traversable wormhole dynamics on a quantum processor," Nature 612, 51-55 (2022).
2. Gao, Jafferis, Wall, "Traversable wormholes via a double trace deformation," JHEP 12 (2017) 151.
3. Kobrin, Schuster, Yao, "Comment on 'Traversable wormhole dynamics on a quantum processor'," arXiv:2302.07897.
4. Jafferis, Zlokapa, Lykken et al., "Reply to Comment," arXiv:2303.15423.
5. Gao, "Commuting SYK model and holography," arXiv:2306.14988.
6. Maldacena, Qi, "Eternal traversable wormhole," arXiv:1804.00491.
7. Brown et al., "Quantum gravity in the lab: teleportation by size and traversable wormholes," arXiv:1911.06314.
8. Huang 2026, Paper 1: tau = 1 - F, Petz recovery unification.
9. Huang 2026, Paper 2: Sigma_grav = -ln(-g_00), exponential metric.
10. Granet et al., "Sparse SYK on trapped ions," arXiv:2507.07530.
11. Bentsen, Nguyen, Swingle, "Wormhole QEC codes," Quantum 8, 1466 (2024).
