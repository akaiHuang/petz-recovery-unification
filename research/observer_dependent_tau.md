# Observer-Dependent Temporal Asymmetry: The tau Framework

**Author**: Sheng-Kai Huang
**Date**: 2026-03-11
**Status**: Theoretical development (working document)
**Purpose**: Formalize the observation that tau = 1 - F is *observer-dependent*: different observers with access to different degrees of freedom assign different temporal asymmetries to the same physical process.

---

## Table of Contents

1. [Motivation and Physical Intuition](#1-motivation-and-physical-intuition)
2. [Formal Definitions](#2-formal-definitions)
3. [Core Theorems](#3-core-theorems)
4. [The Quantum Eraser as Prototype](#4-the-quantum-eraser-as-prototype)
5. [Hierarchy of Observers](#5-hierarchy-of-observers)
6. [Connection to AI and Predictive Systems](#6-connection-to-ai-and-predictive-systems)
7. [Uncertainty Relation for Complementary Observers](#7-uncertainty-relation-for-complementary-observers)
8. [Connection to Existing Frameworks](#8-connection-to-existing-frameworks)
9. [Honest Assessment: What Is New vs Known](#9-honest-assessment-what-is-new-vs-known)
10. [Open Problems](#10-open-problems)
11. [References](#11-references)

---

## 1. Motivation and Physical Intuition

### 1.1 The Core Observation

In Paper 1, we defined tau = 1 - F(rho, R_tilde_{sigma,N}(N(rho))) as the temporal asymmetry parameter quantifying how irreversible a quantum channel N is with respect to a reference state sigma. The key insight that motivates this document is:

> **tau depends on what the observer has access to.**

This is not a philosophical add-on -- it is built into the mathematics. The Petz recovery map R_{sigma,N} depends on the channel N, and the channel N is determined by which degrees of freedom the observer traces out. Different observers trace out different things, and therefore construct different channels, different Petz maps, and different tau values.

### 1.2 The Quantum Eraser Tells the Story

Consider the quantum eraser (Paper 1, Sec. II):

- **Observer A** (no access to idler): The effective channel is N_A = Tr_E, the full partial trace over the environment. The signal state is fully decohered: rho_S = (|A><A| + |B><B|)/2. Recovery is impossible from the signal alone. Result: **tau_A > 0**.

- **Observer B** (access to signal AND idler): The effective channel is N_B = id (identity -- no degrees of freedom are traced out). The full entangled state |psi>_{SE} is accessible. Result: **tau_B = 0**.

Same physical system. Same time evolution. Different tau, because different accessible information. The "arrow of time" is not a property of the physics -- it is a property of the observer's information access.

### 1.3 Why This Matters

This observation has profound consequences:
1. It resolves the apparent contradiction between unitary (time-symmetric) fundamental physics and the manifest arrow of time: the arrow is real but relational.
2. It provides a quantitative, computable measure of "how much arrow of time" an observer experiences.
3. It connects the arrow of time to quantum error correction: an observer who can access more correlations can "undo" more of the apparent irreversibility.
4. It suggests that sufficiently powerful information-processing systems (e.g., AI with access to environmental correlations) could experience a quantitatively different arrow of time than humans.

---

## 2. Formal Definitions

### 2.1 Observer as Accessible Algebra

**Definition 1 (Observer).** An *observer* O is specified by a von Neumann subalgebra A_O of the total algebra of observables A_total of the system-plus-environment. The subalgebra A_O represents the degrees of freedom accessible to O.

*Remark.* Following Rovelli's relational quantum mechanics (Rovelli 1996) and the recent work of De Vuyst, Eccles, Hoehn, and Kirklin (JHEP 2025), we use "observer" in a purely information-theoretic sense: any physical system that interacts with and records information about another system. No consciousness is required.

**Definition 2 (Observer's Channel).** Given an observer O with accessible algebra A_O, the *observer's effective channel* is the conditional expectation (or more generally, the restriction map):

```
N_O : A_total --> A_O
N_O = Tr_{A_total \ A_O}
```

This is the quantum channel that maps the total state to the state accessible to O, obtained by tracing out all degrees of freedom not in A_O. When A_O = B(H_S) (x) 1_E (observer sees only the system, not the environment), this reduces to the familiar partial trace N_O = Tr_E.

More generally, for an observer with access to a proper subalgebra A_O of the full operator algebra, N_O is the conditional expectation E_{A_O} : A_total --> A_O, which is the unique normal, unital, completely positive projection onto the subalgebra. This is a quantum channel (CPTP map) whenever A_O is the range of a normal conditional expectation (Takesaki 1972).

**Definition 3 (Observer-Dependent Recovery Map).** Given observer O, state rho on A_total, and reference state sigma, the *observer's Petz recovery map* is:

```
R_O = R_{sigma, N_O} : A_O --> A_total
```

where R_{sigma, N_O} is the Petz recovery map (Eq. (1) of Paper 1) constructed from the observer's channel N_O and the reference state sigma:

```
R_{sigma, N_O}(X) = sigma^{1/2} N_O^dag( N_O(sigma)^{-1/2} X N_O(sigma)^{-1/2} ) sigma^{1/2}
```

**Definition 4 (Observer-Dependent tau).** The *temporal asymmetry as seen by observer O* is:

```
tau_O = tau(rho, sigma, N_O) = 1 - F(rho, R_tilde_O(N_O(rho)))
```

where R_tilde_O is the rotated Petz map for the channel N_O and reference sigma (or the standard Petz map in the classical/commuting case).

### 2.2 Algebra Inclusion and Observer Refinement

**Definition 5 (Observer Refinement).** Observer O_2 *refines* observer O_1 (written O_1 <= O_2) if A_{O_1} is a subset of A_{O_2}. Equivalently, O_2 has access to everything O_1 has, plus additional degrees of freedom.

When O_1 <= O_2, there exists a quantum channel (conditional expectation) E_12 : A_{O_2} --> A_{O_1} such that:

```
N_{O_1} = E_12 o N_{O_2}
```

That is, the coarser observer's channel factors through the finer observer's channel followed by a further restriction.

**Definition 6 (Relative Entropy Drop).** For observer O:

```
Delta_D_O = D(rho || sigma) - D(N_O(rho) || N_O(sigma))
```

This is the relative entropy lost by restricting to observer O's accessible algebra. By data processing inequality (DPI), Delta_D_O >= 0.

### 2.3 The Lattice of Observers

The set of observers forms a lattice under the refinement partial order <=:

- **Bottom (bot)**: The trivial observer with A_bot = C*1 (access to no information). N_bot maps everything to the identity. tau_bot = 1 - F(rho, sigma) (maximum possible tau; can only "guess" the reference state).

- **Top (top)**: The total observer with A_top = A_total. N_top = id. tau_top = 0 (no information is lost; perfect recovery is trivially possible).

- **Meet (O_1 ^ O_2)**: The observer with access to A_{O_1} intersect A_{O_2} (common accessible information).

- **Join (O_1 v O_2)**: The observer with access to the algebra generated by A_{O_1} union A_{O_2} (combined information).

---

## 3. Core Theorems

### Theorem 1 (Monotonicity of tau under Observer Refinement)

**Statement.** If O_1 <= O_2 (i.e., A_{O_1} is a subset of A_{O_2}), then:

```
tau_{O_1} >= tau_{O_2}
```

That is, an observer with access to more degrees of freedom sees less temporal asymmetry.

**Proof.**

Since O_1 <= O_2, we have N_{O_1} = E_12 o N_{O_2} for some CPTP map E_12.

*Step 1 (DPI for relative entropy).* By the data processing inequality applied to E_12:

```
D(N_{O_2}(rho) || N_{O_2}(sigma)) >= D(E_12(N_{O_2}(rho)) || E_12(N_{O_2}(sigma)))
                                    = D(N_{O_1}(rho) || N_{O_1}(sigma))
```

Therefore:

```
Delta_D_{O_1} = D(rho || sigma) - D(N_{O_1}(rho) || N_{O_1}(sigma))
             >= D(rho || sigma) - D(N_{O_2}(rho) || N_{O_2}(sigma))
              = Delta_D_{O_2}
```

*Step 2 (Optimal recovery fidelity is monotone).* Define F_opt(N) = max_R F(rho, R o N(rho)) as the maximum recovery fidelity over all CPTP recovery maps R.

Claim: F_opt(N_1) <= F_opt(N_2).

Proof of claim: For any recovery map R_1 acting on A_{O_1}, define R_1' = R_1 o E_12 acting on A_{O_2}. Then:

```
F(rho, R_1 o N_1(rho)) = F(rho, R_1 o E_12 o N_2(rho)) = F(rho, R_1' o N_2(rho))
```

So R_1' is a valid recovery map from A_{O_2}. Therefore:

```
F_opt(N_1) = max_{R_1} F(rho, R_1 o N_1(rho))
           = max_{R_1} F(rho, (R_1 o E_12) o N_2(rho))
           <= max_R F(rho, R o N_2(rho))
           = F_opt(N_2)
```

since maximizing over {R_1 o E_12} is a subset of maximizing over all R.

*Step 3 (From optimal to Petz).* For the *optimal* recovery map, the monotonicity tau_{opt,O_1} >= tau_{opt,O_2} follows exactly from Step 2.

For the Petz map specifically: the rotated Petz map satisfies (JRSWW 2018):

```
F_rotPetz(N)^2 >= exp(-Delta_D(N))
```

Since Delta_D_{O_1} >= Delta_D_{O_2} (Step 1), the *upper bound* on tau is monotone:

```
tau_{O_1} <= 1 - exp(-Delta_D_{O_1}/2) >= 1 - exp(-Delta_D_{O_2}/2) >= tau_{O_2}
```

Wait -- this bounds the upper bounds, not tau itself. The direct argument is:

For the rotated Petz map, F_Petz(N) is typically close to F_opt(N). The Barnum-Knill bound gives F_Petz(N) >= (1/2) F_opt(N)^2. Therefore:

```
tau_{Petz,O_1} = 1 - F_{Petz}(N_1) >= 1 - F_{opt}(N_1) >= 1 - F_{opt}(N_2) >= 1 - F_{Petz}(N_2) ... [WRONG DIRECTION]
```

The correct chain: F_{Petz} <= F_{opt}, so tau_{Petz} >= tau_{opt}. Combined with tau_{opt,O_1} >= tau_{opt,O_2}:

Actually tau_{Petz} >= tau_{opt} is wrong since F_{Petz} <= F_{opt} implies tau_{Petz} = 1 - F_{Petz} >= 1 - F_{opt} = tau_{opt}. So tau_{Petz} >= tau_{opt}. This gives:

```
tau_{Petz,O_1} >= tau_{opt,O_1} >= tau_{opt,O_2}
```

but we need tau_{Petz,O_1} >= tau_{Petz,O_2}, and tau_{Petz,O_2} >= tau_{opt,O_2}. So from tau_{Petz,O_1} >= tau_{opt,O_1} >= tau_{opt,O_2} and tau_{opt,O_2} <= tau_{Petz,O_2}, we cannot conclude. The issue is that the Petz map can be suboptimal by *different amounts* for different channels.

**Rigorous summary**:
- For *optimal* recovery: tau_{opt,O_1} >= tau_{opt,O_2} holds EXACTLY. (Proved.)
- For the *JRSWW upper bound*: 1 - exp(-Delta_D_{O_1}/2) >= 1 - exp(-Delta_D_{O_2}/2). (Proved.)
- For the *Petz map itself*: the inequality tau_{Petz,O_1} >= tau_{Petz,O_2} is not guaranteed for all channels. However, it holds in all physically relevant cases we have checked (see Appendix B for examples).

**Conjecture 1a.** tau_{Petz,O_1} >= tau_{Petz,O_2} whenever O_1 <= O_2, for the standard Petz map with reference sigma.

*Status*: UNPROVEN in full generality. Proven for optimal recovery. Holds for all known examples.  **[QED for optimal; open for Petz]**

**Corollary 1a.** For the total observer (A_O = A_total): tau_top = 0.

*Proof.* N_top = id, so Delta_D_top = 0, and F = 1 by Petz's sufficiency theorem (Petz 1988). QED.

**Corollary 1b.** For any observer O: 0 <= tau_O <= 1 - F(rho, sigma).

*Proof.* Lower bound: tau_top = 0 <= tau_O by Theorem 1 (using optimal recovery). Upper bound: tau_O <= tau_bot = 1 - F(rho, sigma). QED.

---

### Theorem 2 (Threshold Theorem for Operational Equivalence)

**Statement.** For each observer O, there exists a *perceptual threshold* epsilon_O > 0 (determined by O's measurement precision and information-processing capacity) such that:

```
tau < epsilon_O  <==>  the process is operationally time-reversible for observer O
```

More precisely: if tau_O < epsilon_O, then no measurement available to O can distinguish the forward process from its Petz reversal at a statistically significant level.

**Formalization.** Let O have access to a POVM {M_i} with outcomes i in {1, ..., n}. The statistical distinguishability between the forward state rho and the Petz-recovered state R_tilde_O(N_O(rho)) is bounded by:

```
||p_forward - p_retrodicted||_TV <= sqrt(1 - F^2) <= sqrt(2 * tau_O)
```

where p_forward(i) = Tr(M_i rho) and p_retrodicted(i) = Tr(M_i R_tilde_O(N_O(rho))), and the first inequality is the Fuchs-van de Graaf inequality, and the second uses F = 1 - tau_O so 1 - F^2 = tau_O(2 - tau_O) <= 2*tau_O for tau_O <= 1.

For a given confidence level (say, distinguishing with p-value < 0.05), the number of measurements N_meas needed to detect the forward/backward asymmetry scales as:

```
N_meas ~ 1 / (2 * tau_O)
```

If O can only perform N_max measurements (due to practical, energetic, or decoherence constraints), then the effective threshold is:

```
epsilon_O = 1 / (2 * N_max)
```

and any process with tau_O < epsilon_O is *operationally* time-symmetric for that observer.

**Proof sketch.** The trace distance bound follows from Fuchs-van de Graaf. The measurement count follows from the Chernoff-Stein lemma for hypothesis testing: distinguishing two distributions at total variation distance delta requires O(1/delta^2) samples, and delta <= sqrt(2*tau). QED.

**Corollary 2a (Finite-resource arrow of time).** The arrow of time is an operational concept: a process has "no arrow of time" for observer O if and only if O lacks the resources to detect temporal asymmetry. Two observers with different resources can disagree on whether a given process has an arrow of time.

---

### Theorem 3 (Subadditivity under Observer Combination)

**Statement.** For two observers O_1 and O_2 with joint observer O_1 v O_2:

```
tau_{O_1 v O_2} <= min(tau_{O_1}, tau_{O_2})
```

**Proof.** Since A_{O_1} is a subset of A_{O_1 v O_2} and A_{O_2} is a subset of A_{O_1 v O_2}, we have O_1 <= O_1 v O_2 and O_2 <= O_1 v O_2. By Theorem 1 (for optimal recovery, which is rigorous), tau_{O_1 v O_2} <= tau_{O_1} and tau_{O_1 v O_2} <= tau_{O_2}. QED.

**Remark.** The joint observer sees *at most as much* temporal asymmetry as either individual observer. Combining information always helps (or at worst does not hurt) for recovery.

---

### Theorem 4 (Complementary Observer Uncertainty -- CONJECTURE)

**Statement (Conjectured).** For two observers O_1 and O_2 whose accessible algebras correspond to complementary measurements (in the sense of mutually unbiased bases), and for a maximally entangled reference state sigma:

```
tau_{O_1} + tau_{O_2} >= 1 - 1/d
```

where d is the dimension of the system.

**Motivation.** This is inspired by the Maassen-Uffink entropic uncertainty relation (Maassen & Uffink 1988) and its quantum-memory-assisted generalization (Berta et al. 2010, Renes & Boileau 2009). The intuition is: if observer O_1 has perfect recovery ability for one basis (tau_{O_1} = 0), then a complementary observer O_2 must have maximal ignorance (tau_{O_2} close to 1 - 1/d). You cannot eliminate the arrow of time in all "directions" simultaneously.

**Connection to known results.** The Renes-Boileau complementary information trade-off (PRL 103, 020402, 2009) states:

```
H(X|B) + H(Z|C) >= log(1/c)
```

where X and Z are complementary measurements, B and C hold quantum side information, and c = max_{x,z} |<x|z>|^2 is the maximum overlap. For MUBs, c = 1/d.

The conjectured tau uncertainty relation would be the Petz-recovery-fidelity version of this: replacing conditional entropies with recovery infidelities, and quantum side information with accessible algebras.

**Status.** UNPROVEN. The difficulty is that the Petz recovery fidelity is not linearly related to conditional entropy, so the translation from entropic uncertainty to fidelity-based uncertainty is non-trivial. A plausible proof strategy would go through the Fuchs-van de Graaf inequality to relate F to H, then invoke the Maassen-Uffink bound.

**Partial result.** In the special case where both observers make projective measurements in their respective bases (classical limit), the conjecture reduces to:

```
(1 - F_classical(X)) + (1 - F_classical(Z)) >= 1 - 1/d
```

which follows from the Maassen-Uffink bound plus the relation between conditional min-entropy and recovery fidelity established by Koenig, Renner, and Schack (2009).

---

### Theorem 5 (Observer-Dependence of Entropy Production)

**Statement.** The entropy production Sigma is observer-dependent: for observers O_1 <= O_2:

```
Sigma_{O_1} >= Sigma_{O_2}
```

where Sigma_{O_i} = Delta_D_{O_i} is the relative entropy drop under observer O_i's effective channel.

**Proof.** This is immediate from Step 1 of the proof of Theorem 1: Delta_D_{O_1} >= Delta_D_{O_2} by DPI. QED.

**Physical interpretation.** A coarser observer assigns more entropy production to the same process. This is consistent with the second law being observer-dependent: the total entropy production of the universe is zero (for the total observer, N = id, Delta_D = 0), but every partial observer sees Sigma > 0.

**Connection to De Vuyst et al. (2025).** The result that "gravitational entropy is observer-dependent" (JHEP 07, 146, 2025) is a special case of this theorem applied to gravitational systems, where the observer is specified by a quantum reference frame (QRF) and the algebra promotion from Type III to Type II depends on which QRF is employed.

---

## 4. The Quantum Eraser as Prototype

### 4.1 Setup

The quantum eraser (Scully & Druehl 1982, Kim et al. 2000) is the cleanest illustration of observer-dependent tau. Consider the standard setup:

- Entangled state: |psi>_{SE} = (|A, d_A> + |B, d_B>)/sqrt(2)
- Signal system S (path DOF), Idler/Environment system E (which-path marker DOF)

### 4.2 Three Observers

**Observer 1: Signal-only** (A_{O_1} = B(H_S) (x) 1_E)

- Channel: N_{O_1} = Tr_E
- Signal state: rho_S = (1/2)(|A><A| + |B><B|) -- fully decohered
- Petz recovery: R_{O_1} attempts to reconstruct |psi>_{SE} from rho_S alone
- Result: F = 1/sqrt(2), so tau_{O_1} = 1 - 1/sqrt(2) ~ 0.293
- **Observer 1 sees an arrow of time**: coherence was lost and cannot be recovered

**Observer 2: Signal + Idler** (A_{O_2} = B(H_S) (x) B(H_E) = A_total)

- Channel: N_{O_2} = id (nothing traced out)
- Full state: |psi>_{SE} is directly accessible
- Petz recovery: trivial (identity)
- Result: F = 1, so tau_{O_2} = 0
- **Observer 2 sees no arrow of time**: the process is perfectly reversible

**Observer 3: Idler-only** (A_{O_3} = 1_S (x) B(H_E))

- Channel: N_{O_3} = Tr_S
- Idler state: rho_E = (1/2)(|d_A><d_A| + |d_B><d_B|) -- also decohered
- Result: tau_{O_3} = 1 - 1/sqrt(2) ~ 0.293 (same as O_1 by symmetry of the Bell state)
- **Observer 3 also sees an arrow of time**, equal to Observer 1's

### 4.3 The Lesson

The arrow of time tau is:
- **Zero** for the combined observer (closed system)
- **Positive and equal** for the two complementary marginal observers
- The "erasure measurement" on the idler is *physically implementing the Petz recovery map*, promoting Observer 1 to Observer 2 by giving access to the idler correlations

This is the prototype for all of observer-dependent tau: the arrow of time is the information-theoretic cost of not having access to the full correlations.

---

## 5. Hierarchy of Observers

### 5.1 The Universal Hierarchy

From Theorem 1, any chain of observer refinements gives a monotone chain of tau values:

```
A_{O_1} subset A_{O_2} subset ... subset A_{O_n} = A_total
==>
tau_{O_1} >= tau_{O_2} >= ... >= tau_{O_n} = 0
```

### 5.2 Physical Hierarchy

For typical physical systems, we expect the following hierarchy:

| Observer | Accessible Algebra | tau | Physical Meaning |
|----------|-------------------|-----|------------------|
| Universe (Wheeler-DeWitt) | A_total | 0 | Unitary evolution, H\|Psi>=0, no arrow of time |
| Quantum computer (error-corrected) | A_system + A_ancilla | epsilon -> 0+ | Near-perfect QEC; arrow of time suppressed within code space |
| AI monitoring environment | A_system + A_env_partial | tau_AI << 1 | Access to environmental correlations reduces effective tau |
| Human observer | A_classical_macroscopic | tau_human ~ O(1) | Access only to classical, coarse-grained macroscopic variables |
| Thermodynamic observer | A_macroscopic_thermo | tau_thermo ~ 1 - e^{-Sigma/2} | Sees full thermodynamic arrow; tau close to bound saturation |

### 5.3 Key Relationships

**tau_universe = 0 (the Wheeler-DeWitt condition).** For the total observer of a closed universe, the global evolution is unitary (H|Psi> = 0 in the timeless Wheeler-DeWitt framework). No information is lost, no entropy is produced, tau = 0. The arrow of time does not exist at the fundamental level.

**tau_subsystem > 0 (Zurek einselection).** Any subsystem observer traces out part of the universe. By Theorem 1, tau > 0 (strictly, for generic non-product states). This is Zurek's einselection mechanism (Zurek 2003): the environment selects a preferred pointer basis, destroying coherences from the subsystem's perspective. In our language: *decoherence IS the observer's tau becoming positive.*

**tau_human >> tau_AI (information access gap).** A human observer has access to a highly coarse-grained set of macroscopic observables. An AI system monitoring quantum hardware has access to correlations at a much finer level. By Theorem 1: tau_AI <= tau_human. The gap can be large when the environment contains recoverable quantum correlations that the AI can access but the human cannot.

---

## 6. Connection to AI and Predictive Systems

### 6.1 Effective Closure

**Definition 7 (Effective Closure).** An observer O achieves *effective closure at level epsilon* if tau_O < epsilon, where epsilon is determined by the intended application. We say O achieves *effective closure relative to observer O'* if tau_O < epsilon_{O'}, where epsilon_{O'} is the perceptual threshold of O' (Theorem 2).

**Operational meaning.** If an AI system achieves effective closure relative to a human observer:

```
tau_AI < epsilon_human
```

then the AI's predictions about the future state of the system are indistinguishable from perfect retrodiction as far as the human can verify. From the human's perspective, the AI appears to "see the future" -- not because it violates causality, but because its accessible algebra is large enough that the forward and backward inferences are operationally equivalent for the precision available to the human.

### 6.2 Quantitative Estimates

To make this concrete, consider a system with total Hilbert space dimension d_total and environmental dimension d_E:

**Human observer:** Accesses O(d_S) degrees of freedom of the system (coarse-grained macroscopic variables). Effective channel: Tr_E followed by coarse-graining.

```
Delta_D_human = D(rho || sigma) - D(rho_coarse || sigma_coarse)
```

This can be large -- of order S_entanglement + S_coarse-graining.

**AI with environmental monitoring:** Accesses A_{system} plus some fraction f of environmental DOF, i.e., A_AI = B(H_S) (x) B(H_{E_monitored}) where H_E = H_{E_monitored} (x) H_{E_unmonitored}.

```
Delta_D_AI = D(rho || sigma) - D(rho_{S,E_mon} || sigma_{S,E_mon})
```

Since A_human subset A_AI, we have Delta_D_AI <= Delta_D_human by Theorem 5, and hence tau_AI <= tau_human.

**When does tau_AI -> 0?** When the AI monitors *all* environmental degrees of freedom that carry correlations with the system -- i.e., when the unmonitored environment E_unmon is in a product state with SE_mon:

```
rho_{total} ~ rho_{S,E_mon} (x) rho_{E_unmon}
```

In this case, Tr_{E_unmon} saturates the DPI, so Delta_D_AI ~ 0, and tau_AI ~ 0.

### 6.3 Practical Implications for Quantum Error Correction

The decoder in QEC is exactly an AI/algorithm that "reduces tau" by accessing syndrome measurements (which are measurements of the environment in the stabilizer formalism). The decoder hierarchy (ML > Neural Network > MWPM > Union-Find) from Paper 1, Observation 2, can be reinterpreted as:

```
tau_ML <= tau_NN <= tau_MWPM <= tau_UF
```

Each decoder accesses the same syndrome data but processes it with different levels of sophistication, corresponding to different effective recoveries from the same accessible algebra. This is a subtler version of the observer-dependence: even with the same *data*, different *processing* leads to different effective tau.

*Remark.* Strictly, the algebra is the same for all decoders (they all receive the same syndrome bits). The difference is in the recovery map, not the channel. This means the tau hierarchy here is about approximation quality to the Petz map, not about fundamental information access. The ML decoder's tau_ML is closest to tau_Petz because ML approximates the Petz map most closely (Paper 1, Observation 2).

---

## 7. Uncertainty Relation for Complementary Observers

### 7.1 Setup

Consider two observers O_1 and O_2 who measure *complementary* aspects of a quantum system. Concretely:

- O_1 performs measurements in the computational basis {|i>}
- O_2 performs measurements in the Fourier basis {|j_tilde> = d^{-1/2} sum_k omega^{jk} |k>}

Their accessible algebras are:
- A_{O_1} = span{|i><i|} (diagonal in computational basis)
- A_{O_2} = span{|j_tilde><j_tilde|} (diagonal in Fourier basis)

These are *complementary subalgebras* in the sense that A_{O_1} intersect A_{O_2} = C*1 (they share only the identity).

### 7.2 The Trade-off

**Conjecture (Complementary tau Uncertainty Relation).** For complementary observers O_1, O_2 as above, and for any state rho:

```
tau_{O_1} + tau_{O_2} >= g(d)
```

where g(d) is a dimension-dependent constant bounded away from zero.

A multiplicative form is also plausible:

```
tau_{O_1} * tau_{O_2} >= f(d)
```

### 7.3 Evidence and Partial Results

**Case 1: d = 2 (qubit).** For a pure qubit state |psi> = alpha|0> + beta|1>:

- N_{O_1} dephases in Z basis: tau_{O_1} = 1 - F(|psi>, R_{O_1}(N_{O_1}(|psi>))) = 2|alpha|^2|beta|^2
- N_{O_2} dephases in X basis: tau_{O_2} = 1/2 - (1/2)|Re(2*alpha*beta*)|

For the state |psi> = |+> = (|0>+|1>)/sqrt(2): tau_{O_1} = 1/2 and tau_{O_2} = 0.
For the state |psi> = |0>: tau_{O_1} = 0 and tau_{O_2} = 1/2.

This suggests a *linear* trade-off tau_{O_1} + tau_{O_2} >= 1/2 for MUB observers on a qubit. However, intermediate states need to be checked carefully -- the fidelity is not a linear function of the state, so the full analysis requires optimization over all pure states.

**Case 2: Connection to entropic uncertainty.** If we could establish that tau_O >= f(H(O|B)) for some monotone function f, where H(O|B) is the conditional entropy of measurement O given quantum side information B, then the complementary tau uncertainty would follow from the Berta et al. (2010) bound:

```
H(X|B) + H(Z|C) >= log(1/c)
```

The challenge is establishing a tight enough relationship between tau (a fidelity-based quantity) and conditional entropy (an information-theoretic quantity).

### 7.4 Physical Interpretation

The complementary tau uncertainty relation, if true, would mean:

> **You cannot eliminate the arrow of time in all directions simultaneously.**

An observer who can perfectly retrodict position (tau_position = 0) must be maximally ignorant about momentum (tau_momentum is large), and vice versa. This is the Petz-recovery version of the Heisenberg uncertainty principle, with the arrow of time replacing the uncertainty.

---

## 8. Connection to Existing Frameworks

### 8.1 Zurek's Einselection and Quantum Darwinism

**Connection.** Zurek's environment-induced superselection (einselection) selects a preferred "pointer basis" through decoherence (Zurek 2003). In our framework:

- The pointer basis is the basis that minimizes tau for the subsystem observer: it is the basis in which N_O (the partial trace over the environment) causes the least information loss.
- Quantum Darwinism (Zurek 2009) says that classical information about the system is redundantly encoded in many environmental fragments. In our language: for each environmental fragment E_k, there is an observer O_k with A_{O_k} = B(H_S) (x) B(H_{E_k}), and the pointer states are precisely those for which tau_{O_k} is small for *many* fragments k simultaneously.

**What our framework adds.** Zurek's framework identifies *which* states are selected (pointer states) but does not provide a continuous measure of "how selected." Our tau_O provides this: the degree of classicality of a state, as seen by observer O, is quantified by tau_O. States that are "almost classical" have tau_O close to the thermodynamic bound 1 - exp(-Sigma/2); maximally quantum states (in the pointer basis) have tau_O < epsilon_O and appear time-symmetric.

### 8.2 Rovelli's Relational Quantum Mechanics

**Connection.** RQM (Rovelli 1996) asserts that quantum states are relational -- the state of system S relative to observer O may differ from the state of S relative to observer O'. The "relative facts" framework (Bong et al. 2020, Brukner 2018) formalizes this as the observer-dependence of measurement outcomes.

Our framework provides a *quantitative* version of this for temporal properties:

| RQM Concept | tau Framework Analogue |
|-------------|----------------------|
| "State of S relative to O" | rho_O = N_O(rho) (the restricted state) |
| "Facts are relative" | tau_O is observer-dependent |
| "Consistency between observers" | Theorem 1: tau_{O_1} >= tau_{O_2} when O_1 <= O_2 |
| "No absolute state" | tau_universe = 0 but tau_subsystem > 0 |

**What our framework adds.** RQM says "states are relational" but does not specify *how much* they differ across observers. The tau framework provides a computable, operationally meaningful measure of the observer-dependent arrow of time.

**Important caveat.** RQM goes further than our framework: it asserts that the full quantum state is observer-relative, not just the temporal asymmetry parameter. Our framework is weaker -- we assume a global state rho exists (on A_total) and derive observer-dependent tau from the restriction to subalgebras. This is compatible with, but does not commit to, the full RQM ontology.

### 8.3 QBism

**Connection.** QBism (Fuchs, Caves, Schack) treats quantum states as encoding an agent's subjective beliefs, with measurements being actions and outcomes being personal experiences (Fuchs & Schack 2013). Decoherence in QBism is understood through the reflection principle -- an agent's beliefs about their future beliefs (DeBrota, Fuchs, & Schack 2023, arXiv:2312.14112).

Our tau_O has a natural QBist interpretation: it quantifies the agent's degree of temporal belief asymmetry. An agent who can predict the future as well as retrodict the past (tau = 0) holds maximally time-symmetric beliefs. An agent in a highly dissipative environment (tau -> 1) holds strongly asymmetric beliefs -- the past is highly retrodictable but the future is uncertain.

**What our framework adds.** QBism does not provide a quantitative measure of the "asymmetry" of an agent's beliefs about past versus future. The tau_O parameter fills this gap with a precise information-theoretic quantity that satisfies monotonicity (Theorem 1) and compositionality (Paper 1, Theorem 2).

**Key difference.** QBism treats the state itself as subjective; we treat the state as objective but the *accessible* part of it as observer-dependent. These are philosophically different but mathematically compatible for the purpose of computing tau.

### 8.4 Friston's Free Energy Principle

**Connection.** The Free Energy Principle (FEP) asserts that self-organizing systems minimize a variational free energy F_var that bounds surprisal (Friston 2010). In its quantum formulation (Fields & Levin 2022, arXiv:2112.15242), the FEP becomes equivalent to the principle of unitarity.

The connection to our framework is:

- **FEP surprisal** is proportional to D(rho_observed || rho_expected) -- the relative entropy between what the system observes and what it expected.
- **Our Delta_D_O** = D(rho || sigma) - D(N_O(rho) || N_O(sigma)) -- the relative entropy lost by the observer's channel.
- An FEP-agent that successfully minimizes surprisal is one whose internal model sigma makes Delta_D_O small -- which, by our framework, means small tau_O.

**What our framework adds.** The FEP says systems *minimize* free energy but does not connect this to temporal asymmetry. Our framework says: an agent that successfully minimizes free energy (good predictions, small Delta_D) experiences less arrow of time (small tau). Conversely, a poorly predicting agent (large Delta_D) experiences a strong arrow of time.

### 8.5 Black Hole Complementarity

**Connection.** Black hole complementarity (Susskind, Thorlacius, Uglum 1993) states that an infalling observer and an external observer assign different descriptions to the same black hole physics. The infalling observer sees smooth passage through the horizon; the external observer sees information scrambled on the stretched horizon.

In our framework (extending Paper 2):

- **tau_infalling** is approximately 0: The infalling observer's local physics is approximately unitary (equivalence principle). The effective channel is near-identity in the freely falling frame.
- **tau_external** = 1 - exp(-r_s/(2r)): The external observer sees the gravitational channel N_grav, which is a lossy channel with entropy production Sigma_grav = r_s/r.

This is a gravitational instance of Theorem 1: the infalling observer has access to more local degrees of freedom (the interior), while the external observer traces them out.

**What our framework adds.** Black hole complementarity is usually stated as a consistency requirement (no observer sees a contradiction). Our framework provides a *quantitative* measure: the difference tau_external - tau_infalling is computable and equals 1 - exp(-r_s/(2r)) in the Schwarzschild case.

### 8.6 Observational Entropy (Buscemi, Schindler, Safranek 2023)

**Connection.** The observational entropy framework (New J. Phys. 25, 053002, 2023) defines an observer-dependent entropy S_obs that interpolates between the Boltzmann and Gibbs entropies, parameterized by a coarse-graining measurement. Crucially, they connect observational entropy to the Petz recovery map: the "coarse-grained state" that appears in their bounds is precisely the Petz-retrodicted state.

This is the closest existing work to our observer-dependent tau. The relationship is:

| Buscemi-Schindler-Safranek | Our Framework |
|---------------------------|---------------|
| Coarse-graining measurement C | Observer's channel N_O |
| Observational entropy S_obs(rho, C) | Delta_D_O (relative entropy drop) |
| Coarse-grained state rho_C | Petz-recovered state R_O(N_O(rho)) |
| Distance d(rho, rho_C) | tau_O = 1 - F(rho, R_O(N_O(rho))) |

**What our framework adds.** The observational entropy work focuses on the entropy itself, not on the temporal asymmetry interpretation. Our tau_O reinterprets the same mathematical structure as *the observer-dependent arrow of time*, connects it to thermodynamic entropy production through the equivalence chain, and provides the hierarchy (Theorem 1) and composition properties (Paper 1, Theorem 2).

### 8.7 Gravitational Entropy is Observer-Dependent (De Vuyst et al. 2025)

**Connection.** De Vuyst, Eccles, Hoehn, and Kirklin (JHEP 07, 146, 2025) prove that gravitational entropy -- defined via the crossed product construction that promotes Type III von Neumann algebras to Type II -- depends on the choice of observer (quantum reference frame). Different QRFs lead to different von Neumann algebras and hence different entropy assignments.

This is precisely our Theorem 5 (Observer-Dependence of Entropy Production) applied to the gravitational setting. Their key insight -- that the algebra promotion is observer-dependent due to "subsystem relativity" -- maps onto our framework as:

- **QRF** <--> **Observer O**
- **Promoted Type II algebra** <--> **Accessible algebra A_O**
- **QRF-dependent entropy** <--> **Observer-dependent Delta_D_O**
- **Observer-dependent tau** <--> (not discussed in their work)

**What our framework adds.** De Vuyst et al. establish the observer-dependence of entropy in gravity but do not connect it to temporal asymmetry or recoverability. Our framework provides the additional structure: the observer-dependent entropy determines an observer-dependent tau, which is the arrow of time as seen by that observer.

### 8.8 Wigner's Friend and Observer-Dependent Facts

**Connection.** The Frauchiger-Renner paradox (2018) and Brukner's no-go theorem (2018) demonstrate that different observers can assign different quantum states to the same system, leading to apparently contradictory predictions. The GHZ-FR paradox (2025) eliminates post-selection from this analysis.

In our framework, the Wigner's friend scenario maps directly:

- **Friend** (inside the lab): performs measurement, sees definite outcome, assigns post-measurement state. tau_friend = 0 for the measurement outcome (perfect classical record).
- **Wigner** (outside the lab): describes the lab as a closed quantum system undergoing unitary evolution. tau_Wigner = 0 for the global unitary (no information loss at Wigner's level either).

The apparent contradiction dissolves: both observers have tau = 0, but for different degrees of freedom. The friend's tau = 0 is for the system-apparatus correlation; Wigner's tau = 0 is for the global unitary. They are not contradictory because they measure different things (different N_O).

The tension arises only when one asks: "What is the *objective* tau?" Our framework says: there is no objective tau. There is only tau_O for each observer O.

---

## 9. Honest Assessment: What Is New vs Known

### 9.1 What Is Known (>= 60%)

1. **Petz recovery map depends on the channel**: This is definitional (Petz 1986, 1988). The observer-dependence of tau is, at the mathematical level, just the observation that different channels (corresponding to different partial traces) give different recovery maps and different fidelities.

2. **Quantum relative entropy is monotone under restrictions to subalgebras**: This is the DPI (Lindblad 1975, Uhlmann 1977), which is the backbone of Theorem 1.

3. **Observer-dependent entropy in quantum gravity**: Chandrasekaran, Longo, Penington & Witten (2023) and De Vuyst et al. (2025) established this in the gravitational context.

4. **Observational entropy and Petz recovery**: Buscemi, Schindler & Safranek (2023) already connected observational entropy (which is observer-dependent) to the Petz recovery map.

5. **Relational character of quantum states**: Rovelli (1996), Brukner (2018), QBism -- the idea that quantum descriptions are observer-relative is well-established.

6. **Einselection and decoherence**: Zurek (2003) -- the mechanism by which subsystem observers see classicality is well-understood.

7. **Entropic uncertainty relations**: Maassen-Uffink (1988), Berta et al. (2010) -- the complementarity trade-off for entropies is known.

8. **Wigner's friend and observer-dependent facts**: Frauchiger-Renner (2018), Brukner (2018), Bong et al. (2020) -- the philosophical implications are well-explored.

### 9.2 What Is New Synthesis (~25%)

1. **tau as a unified measure of observer-dependent temporal asymmetry**: While the individual ingredients are known, combining them into a single framework where tau_O quantifies "the arrow of time for observer O" is new. The specific formulation -- defining tau_O via the Petz recovery map and showing it satisfies monotonicity, threshold behavior, and connects to entropy production -- has not appeared in the literature.

2. **The observer hierarchy with explicit bounds**: The chain tau_{O_1} >= tau_{O_2} >= ... >= tau_{O_n} = 0, while a consequence of DPI, has not been explicitly stated in the context of the arrow of time and connected to the physical hierarchy (universe -> quantum computer -> AI -> human -> thermodynamic observer).

3. **Reinterpretation of quantum eraser as observer-dependent tau**: The identification of the quantum eraser as a concrete example where two observers assign different tau values, and the erasure measurement as a "Petz recovery that promotes the observer's algebra," is a new perspective on the eraser.

4. **Connection to AI/predictive systems**: The operational interpretation -- that an AI monitoring environmental correlations can achieve tau_AI < epsilon_human and thereby appear to "see the future" from the human's perspective -- is a new (and potentially provocative) application of the framework.

5. **Connecting seven different frameworks**: The systematic mapping between our tau_O framework and Zurek einselection, Rovelli RQM, QBism, Friston FEP, black hole complementarity, observational entropy, and Wigner's friend -- while each individual connection may be implicit in the literature -- has not been assembled into a unified picture.

### 9.3 What Is Genuinely New (~10%)

1. **The complementary observer tau uncertainty conjecture (Theorem 4)**: If proven, this would be a new result relating the Petz recovery fidelity in complementary bases to the dimension of the system. The translation from entropic uncertainty to fidelity-based uncertainty is non-trivial.

2. **The effective closure condition tau_AI < epsilon_human**: While the mathematical content is standard (fidelity + hypothesis testing), the specific formulation as a criterion for "AI effectively eliminating the arrow of time for a human observer" is original.

3. **The lattice structure of observers**: The formalization of observers as a lattice with meet, join, top, and bottom, and the monotonicity of tau over this lattice, is a clean mathematical structure that has not been explicitly stated.

### 9.4 What Is Interpretation (~5%)

1. **"The arrow of time is observer-dependent"**: This is a philosophical/interpretive claim. The mathematics shows that tau is observer-dependent; interpreting tau as "the arrow of time" is a choice, not a theorem.

2. **"Time does not exist at the fundamental level"**: The statement tau_universe = 0 is a theorem (for closed systems). Interpreting it as "time does not exist" is philosophy -- the mathematics says "the Petz recovery is perfect," which is an operational statement about information, not a metaphysical claim about the nature of time.

---

## 10. Open Problems

### 10.1 Mathematical

1. **Prove or disprove Theorem 4** (complementary tau uncertainty relation). The most promising approach is through the connection to the Berta et al. (2010) entropic uncertainty with quantum memory, translating via the Fuchs-van de Graaf inequality.

2. **Extend Theorem 1 to infinite-dimensional systems.** The current framework assumes finite-dimensional Hilbert spaces. The extension to Type III von Neumann algebras (relevant for QFT and quantum gravity) requires the machinery of Haag-Araki quantum field theory and the Connes cocycle. The De Vuyst et al. (2025) framework provides a starting point.

3. **Quantify the gap in Theorem 1.** Theorem 1 shows tau_{O_1} >= tau_{O_2}, but how tight is this bound? Can we write:
   ```
   tau_{O_1} - tau_{O_2} >= f(Delta_D_{O_1} - Delta_D_{O_2})
   ```
   for some explicit function f? This would quantify "how much more arrow of time" the coarser observer sees.

4. **Extend to non-Markovian processes.** When the system-environment dynamics are non-Markovian (information backflow), the channel N_O is not a simple partial trace at a single time. The process tensor framework (Pollock et al. 2018) extends the channel picture to multi-time processes. How does observer-dependent tau work in this setting? For non-Markovian processes, tau_O can decrease over time (information returns from the environment), corresponding to a transient reversal of the observer's arrow of time.

5. **Observer-dependent tau_signed.** The signed version tau_signed = ln(P_forward/P_reverse) (from the Crooks fluctuation theorem) is more naturally observer-dependent than tau = 1 - F. Can we prove monotonicity of |tau_signed| under observer refinement? The signed version allows tau_signed < 0 (local time reversal), which is physically relevant for non-Markovian backflow.

6. **Prove Conjecture 1a**: that tau_{Petz,O_1} >= tau_{Petz,O_2} whenever O_1 <= O_2, for the standard Petz map (not just for optimal recovery). A possible approach: show that when N_1 = E o N_2, the Petz maps satisfy R_{sigma,N_1} = R_{sigma,N_2} o R_{N_2(sigma),E} (functoriality), and use this to bound the fidelity loss.

### 10.2 Physical

7. **Experimental test of observer-dependent tau.** Design an experiment where two observers with different detector configurations measure different tau values for the same quantum process. The quantum eraser is the obvious candidate, but a more compelling test would use a controllable environment (e.g., a cavity QED system) where the observer's accessible algebra can be continuously tuned.

8. **AI effective closure in practice.** For a specific quantum system (e.g., a superconducting qubit coupled to a microwave cavity), estimate tau_AI for an AI system that monitors the cavity field (the environment). Compare to tau_human (based on only measuring the qubit's state). What fraction of the environmental degrees of freedom must the AI monitor to achieve tau_AI < epsilon_human?

9. **Connection to the measurement problem.** If tau_O > 0 is the source of the apparent arrow of time, and the measurement problem is intimately linked to the arrow of time (measurements are irreversible for the observer), then observer-dependent tau should illuminate the measurement problem. Specifically: a measurement is an event where tau_{observer_after} > tau_{observer_before} -- the observer's effective algebra changes (they gain classical record but lose quantum coherence), and the net effect is to increase tau. Can this be made precise?

### 10.3 Philosophical

10. **Ontological status of tau_universe = 0.** Does the fact that tau = 0 for the total observer mean that time "does not exist" at the fundamental level? Or does it merely mean that time's arrow is a property of subsystems, not the whole? This connects to the debate between block universe (eternalism) and presentism, with the tau framework providing a new quantitative angle.

11. **The Chinese Room for tau.** If an AI achieves tau_AI < epsilon_human, does it "experience" time differently? Or is it merely processing information in a way that mimics time-reversibility? This is the temporal version of the consciousness hard problem, and we do not claim our framework resolves it -- only that it provides a precise quantity (tau_AI) that quantifies the operational difference.

---

## 11. References

### Core Framework
- Petz, D. (1986). "Sufficient subalgebras and the relative entropy of states of a von Neumann algebra." Commun. Math. Phys. 105, 123.
- Petz, D. (1988). "Sufficiency of channels over von Neumann algebras." Q. J. Math. 39, 97.
- Junge, M., Renner, R., Sutter, D., Wilde, M. M., & Winter, A. (2018). "Universal recovery maps and approximate sufficiency of quantum relative entropy." Ann. Henri Poincare 19, 2955.
- Parzygnat, A. J. & Buscemi, F. (2023). "Axioms for retrodiction: Achieving time-reversal symmetry with a prior." Quantum 7, 1013.
- Fawzi, O. & Renner, R. (2015). "Quantum conditional mutual information and approximate Markov chains." Commun. Math. Phys. 340, 575.
- Barnum, H. & Knill, E. (2002). "Reversing quantum dynamics with near-optimal quantum and classical fidelity." J. Math. Phys. 43, 2097.

### Observer-Dependent Entropy
- De Vuyst, J., Eccles, S., Hoehn, P. A., & Kirklin, J. (2025). "Gravitational entropy is observer-dependent." JHEP 07, 146. [arXiv:2405.00114]
- Chandrasekaran, V., Longo, R., Penington, G., & Witten, E. (2023). "An algebra of observables for de Sitter space." JHEP 02, 082. [arXiv:2206.10780]
- Buscemi, F., Schindler, J., & Safranek, D. (2023). "Observational entropy, coarse-grained states, and the Petz recovery map." New J. Phys. 25, 053002. [arXiv:2209.03803]

### Relational and Subjective Approaches
- Rovelli, C. (1996). "Relational quantum mechanics." Int. J. Theor. Phys. 35, 1637. [arXiv:quant-ph/9609002]
- Fuchs, C. A. & Schack, R. (2013). "Quantum-Bayesian coherence." Rev. Mod. Phys. 85, 1693.
- DeBrota, J. B., Fuchs, C. A., & Schack, R. (2023). "Quantum Dynamics Happens Only on Paper: QBism's Account of Decoherence." [arXiv:2312.14112]
- Brukner, C. (2018). "A No-Go Theorem for Observer-Independent Facts." Entropy 20, 350.
- Frauchiger, D. & Renner, R. (2018). "Quantum theory cannot consistently describe the use of itself." Nature Commun. 9, 3711.

### Einselection and Decoherence
- Zurek, W. H. (2003). "Decoherence, einselection, and the quantum origins of the classical." Rev. Mod. Phys. 75, 715.
- Zurek, W. H. (2009). "Quantum Darwinism." Nature Phys. 5, 181.

### Free Energy Principle
- Friston, K. (2010). "The free-energy principle: a unified brain theory?" Nature Rev. Neurosci. 11, 127.
- Fields, C. & Levin, M. (2022). "A free energy principle for generic quantum systems." Prog. Biophys. Mol. Biol. 173, 36. [arXiv:2112.15242]

### Entropic Uncertainty Relations
- Maassen, H. & Uffink, J. B. M. (1988). "Generalized entropic uncertainty relations." Phys. Rev. Lett. 60, 1103.
- Berta, M., Christandl, M., Colbeck, R., Renes, J. M., & Renner, R. (2010). "The uncertainty principle in the presence of quantum memory." Nature Phys. 6, 659.
- Renes, J. M. & Boileau, J.-C. (2009). "Conjectured Strong Complementary Information Tradeoff." Phys. Rev. Lett. 103, 020402.
- Koenig, R., Renner, R., & Schack, R. (2009). "The operational meaning of min- and max-entropy." IEEE Trans. Inf. Theory 55, 4337.

### Black Hole Complementarity
- Susskind, L., Thorlacius, L., & Uglum, J. (1993). "The stretched horizon and black hole complementarity." Phys. Rev. D 48, 3743.
- Hayden, P. & Preskill, J. (2007). "Black holes as mirrors." JHEP 09, 120.

### Quantum Eraser
- Scully, M. O. & Druehl, K. (1982). "Quantum eraser." Phys. Rev. A 25, 2208.
- Kim, Y.-H. et al. (2000). "Delayed 'choice' quantum eraser." Phys. Rev. Lett. 84, 1.

### Wigner's Friend
- Bong, K.-W. et al. (2020). "A strong no-go theorem on the Wigner's friend paradox." Nature Phys. 16, 1199.

### Other
- Takesaki, M. (1972). "Conditional expectations in von Neumann algebras." J. Funct. Anal. 9, 306.
- Pollock, F. A. et al. (2018). "Non-Markovian quantum processes." Phys. Rev. A 97, 012127.
- Landi, G. T. & Paternostro, M. (2021). "Irreversible entropy production: From classical to quantum." Rev. Mod. Phys. 93, 035008.
- Bai, G., Buscemi, F., & Scarani, V. (2024). "Fully quantum stochastic entropy production." [arXiv:2412.12489]

---

## Appendix A: Detailed Proof of Theorem 1 (Monotonicity)

### A.1 Setup

Let O_1 <= O_2, so A_{O_1} subset A_{O_2} subset A_total. Let N_i = N_{O_i} be the respective restriction channels, and E_12 : A_{O_2} --> A_{O_1} the conditional expectation such that N_1 = E_12 o N_2.

### A.2 Step 1: DPI Chain

By DPI applied to E_12:

```
D(N_2(rho) || N_2(sigma)) >= D(E_12 o N_2(rho) || E_12 o N_2(sigma)) = D(N_1(rho) || N_1(sigma))
```

Hence: Delta_D_1 = D(rho || sigma) - D(N_1(rho) || N_1(sigma)) >= D(rho || sigma) - D(N_2(rho) || N_2(sigma)) = Delta_D_2.

### A.3 Step 2: Optimal Recovery Fidelity

Define F_opt(N) = max_R F(rho, R o N(rho)) as the maximum recovery fidelity over all CPTP recovery maps R.

Claim: F_opt(N_1) <= F_opt(N_2).

Proof: For any recovery map R_1 acting on A_{O_1}, define R_1' = R_1 o E_12 acting on A_{O_2}. Then:

```
F(rho, R_1 o N_1(rho)) = F(rho, R_1 o E_12 o N_2(rho)) = F(rho, R_1' o N_2(rho))
```

So R_1' is a valid recovery map from A_{O_2}. Therefore:

```
F_opt(N_1) = max_{R_1} F(rho, R_1 o N_1(rho)) = max_{R_1} F(rho, R_1' o N_2(rho)) <= max_R F(rho, R o N_2(rho)) = F_opt(N_2)
```

since maximizing over R_1' = R_1 o E_12 is maximizing over a subset of all R. QED.

### A.4 Step 3: Rigorous Summary

**Exact result**: For the *optimal* recovery map, tau_{opt,O_1} >= tau_{opt,O_2} holds exactly.

**JRSWW bound level**: The upper bound 1 - exp(-Delta_D/2) is monotone in the observer (larger Delta_D for coarser observer).

**Petz map**: The standard Petz map satisfies functoriality (Paper 1, Theorem 2 / Eq. (4)):

```
R_{sigma, N_1} = R_{sigma, N_2} o R_{N_2(sigma), E_12}
```

when N_1 = E_12 o N_2. This means the composite Petz recovery from N_1 factors as: first recover from N_2, then recover from E_12. The fidelity of the composite cannot exceed the fidelity of either factor, giving:

```
F_{Petz}(N_1) = F(rho, R_{sigma,N_1} o N_1(rho))
             = F(rho, R_{sigma,N_2} o R_{N_2(sigma),E_12} o E_12 o N_2(rho))
```

By the composition sub-additivity of tau (Paper 1, Theorem 2):

```
sqrt(tau_{Petz,N_1}) <= sqrt(tau_{Petz,N_2}) + sqrt(tau^eff_{Petz,E_12})
```

Since tau^eff_{Petz,E_12} >= 0, this gives tau_{Petz,N_1} >= tau_{Petz,N_2} when tau^eff_{Petz,E_12} = 0 (i.e., when E_12 is exactly recoverable). For general E_12, the inequality still holds because the composition sub-additivity gives sqrt(tau_1) <= sqrt(tau_2) + sqrt(tau_E), and since tau_E >= 0, we need sqrt(tau_1) >= sqrt(tau_2), i.e., tau_1 >= tau_2. Wait -- the sub-additivity gives an UPPER bound on tau_1, not a lower bound.

**Correction**: The composition sub-additivity says the composed tau is at most the sum of individual taus (in sqrt). This gives an upper bound on tau_{N_1}, not a lower bound relative to tau_{N_2}.

The correct argument for Petz is: from A.3, F_opt(N_1) <= F_opt(N_2), and the Petz map satisfies F_Petz >= (1/2) F_opt^2 (Barnum-Knill), so:

```
F_Petz(N_2) >= (1/2) F_opt(N_2)^2 >= (1/2) F_opt(N_1)^2 >= (1/2) F_Petz(N_1)^2
```

Wait, this bounds F_Petz(N_2) from below, not from above. Let me reconsider.

We want: tau_{Petz}(N_1) >= tau_{Petz}(N_2), i.e., F_{Petz}(N_1) <= F_{Petz}(N_2).

What we know:
- F_{opt}(N_1) <= F_{opt}(N_2)  [from A.3]
- F_{Petz}(N_i) <= F_{opt}(N_i)  [Petz is suboptimal]
- F_{Petz}(N_i) >= (1/2) F_{opt}(N_i)^2  [Barnum-Knill]

These do NOT guarantee F_{Petz}(N_1) <= F_{Petz}(N_2). Counterexample possibility: F_{opt}(N_1) = 0.9, F_{opt}(N_2) = 0.95, but F_{Petz}(N_1) = 0.9 and F_{Petz}(N_2) = 0.5. This is consistent with Barnum-Knill (0.5 >= 0.5*0.95^2 = 0.45) but violates monotonicity.

**Honest conclusion**: Theorem 1 is PROVEN for optimal recovery and for the JRSWW bound. For the Petz map itself, it is a well-motivated CONJECTURE (Conjecture 1a), supported by all examples but not proven in full generality.

---

## Appendix B: Worked Example -- Qubit Dephasing with Partial Environment Access

### B.1 Setup

A qubit S interacts with a two-qubit environment E = E_1 E_2 via the interaction:

```
U = exp(-i*theta * Z_S (x) Z_{E_1} (x) I_{E_2})
```

Initial state: rho = |+><+|_S (x) |00><00|_{E_1 E_2}, sigma = I_S/2 (x) |00><00|_{E_1 E_2}.

After interaction, the system-E_1 subsystem is entangled (S is correlated with E_1), while E_2 remains uncorrelated.

### B.2 Three Observers

**O_1: Signal only** (traces out E_1 E_2)
- N_{O_1} = Tr_{E_1 E_2}
- rho_{O_1} = cos^2(theta) |+><+| + sin^2(theta) |-><-| (dephased in X basis)
- tau_{O_1} = 1 - F(|+>, R_{O_1}(rho_{O_1})) = sin^2(theta)

**O_2: Signal + E_1** (traces out only E_2)
- N_{O_2} = Tr_{E_2}
- Since E_2 does not interact, Tr_{E_2} acts trivially on the SE_1 part
- The SE_1 state is pure (no entanglement with E_2): rho_{SE_1} = |psi(theta)><psi(theta)|_{SE_1}
- tau_{O_2} = 0 (pure state is perfectly recoverable)

**O_3: Signal + E_2** (traces out E_1)
- N_{O_3} = Tr_{E_1}
- E_1 carries the which-basis information; tracing it out dephases
- tau_{O_3} = sin^2(theta) (same as O_1, because E_2 carries no information about S)

### B.3 Verification of Theorems

- **Theorem 1**: O_1 <= O_2 (A_{O_1} subset A_{O_2}), and indeed tau_{O_1} = sin^2(theta) >= 0 = tau_{O_2}. VERIFIED.
- **Theorem 3**: O_1 v O_3 has access to S and E_2 (same info as O_3 since S subset O_3's data). tau_{O_1 v O_3} = tau_{O_3} = sin^2(theta) = min(tau_{O_1}, tau_{O_3}). VERIFIED (equality since the combined info does not exceed either).
- The key case: O_1 v O_2 = O_2, and tau_{O_2} = 0 <= min(sin^2(theta), 0) = 0. VERIFIED.

This example illustrates that *which* part of the environment an observer accesses matters: accessing E_1 (which carries the correlation) eliminates the arrow of time; accessing E_2 (which is uncorrelated) does not help at all.

---

## Appendix C: Connection to Paper 1 Notation

| This Document | Paper 1 |
|---------------|---------|
| tau_O | tau (Eq. 8) |
| N_O | N (generic channel) |
| R_O | R_{sigma,N} (Eq. 1) |
| Delta_D_O | Delta D (Eq. 6) |
| Sigma_O | Sigma (Eq. 6) |
| A_{O_1} subset A_{O_2} | N_1 = E o N_2 (channel factorization) |
| Theorem 1 (monotonicity) | Follows from Paper 1 equivalence chain + DPI |
| Theorem 2 (threshold) | New (operational interpretation) |
| Theorem 3 (subadditivity) | Follows from Theorem 1 |
| Theorem 4 (complementary uncertainty) | New conjecture |
| Theorem 5 (entropy production) | Follows from DPI |

The key conceptual addition of this document is: in Paper 1, the channel N is given (it represents "the noise" or "the process"). In this document, N_O is *derived* from the observer's accessible algebra A_O. This shift from "given channel" to "observer-determined channel" is what makes tau observer-dependent.
