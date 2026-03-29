# Engineering Reversibility: Quantitative Tools for the Arrow of Time

## Paper Outline (Research Planning Document)

**Date**: 2026-03-16
**Author**: Sheng-Kai Huang
**Status**: Outline / Pre-Draft

---

## 1. Title and Abstract

**Title**: Engineering Reversibility: Quantitative Tools for the Arrow of Time

**Alternative titles**:
- Bayesian Quantum Chronometry: Composition, Saturation, and the Cost of Time Reversal
- The Geometry of Irreversibility: From Petz Recovery to Quantum Circuit Design

**Abstract (draft, 2-3 sentences)**:

We develop quantitative tools for engineering reversibility in quantum processes, building on the Petz recovery framework for temporal asymmetry. We establish four results: (A) a computation lower bound tau >= alpha * C^2 linking irreversibility to output noncommutativity, with alpha increasing with Hilbert space dimension; (B) a corrected saturation theorem showing that exact reversibility requires output commutativity [N(rho), N(sigma)] = 0 rather than input commutativity, with the Jensen gap J quantifying deviation; (C) continuous-time dynamics dau/dt = (Sigma_dot/2)(1-tau) verified to machine precision for unitary and Lindbladian evolution; and (D) a Bayesian composition engine achieving per-gate tau reductions up to 94% and 7.7% overall improvement in a 5-gate benchmark. Together, these provide a complete toolkit for diagnosing, bounding, and minimizing the arrow of time in quantum circuits.

---

## 2. Section Structure

### Section I: Introduction
- Context: Paper 1 established tau = 1 - F as temporal asymmetry measure with equivalence chain (Eraser <-> Retrodiction <-> QEC <-> Thermo) and two structural theorems (Saturation, Composition)
- Motivation: Paper 1 was the "dictionary" -- this paper provides the "engineering manual"
- Transition from static (single-channel) to dynamic (multi-channel, continuous-time) analysis
- Preview of four results A-B-C-D and their practical implications
- Key claim: the Petz map is not just a theoretical construct but an actionable tool for quantum circuit design

**Maps to**: Framing of all four results (A/B/C/D)

### Section II: Preliminaries and Notation
- Recall tau = 1 - F(rho, R_{sigma,N}(N(rho))), Petz map definition, JRSWW bound
- Composition theorem (Paper 1, Theorem 5): sqrt(tau_12) <= sqrt(tau_1) + sqrt(tau_2^eff)
- Saturation theorem (Paper 1, Theorem 4): F^2 = exp(-DeltaD) iff conditions (I)-(IV)
- Define C = ||[N(rho), N(sigma)]||_1 (output noncommutativity)
- Define Jensen gap J = <psi|sigma|psi> / exp(<psi|ln sigma|psi>) >= 1
- Lindbladian framework: L_super, E_t = exp(Lt), continuous channels

### Section III: Computation Lower Bound (Result A)

#### III.A: The Conjecture and Numerical Evidence
- Statement: tau >= alpha(d) * C^2, where C = ||[N(rho), N(sigma)]||_1
- Motivation: quantifies the "cost" of quantum noncommutativity -- if channel outputs don't commute, irreversibility is inevitable
- 50,000-trial Monte Carlo at d=2: alpha*(d=2) ~ 0.28 for random CPTP channels
- Channel-specific results: AD ~ 0.001, Depolarizing ~ 0.06, Dephasing ~ 0.15
- The AD lower alpha is structural: AD breaks U(1) symmetry differently

#### III.B: Dimension Scaling -- The Surprise
- alpha*(d=2) ~ 0.28, alpha*(d=3) ~ 0.12, alpha*(d=4) ~ 0.08
- Critical finding: alpha* * d^2 does NOT converge to a constant
- The 1/d^2 hypothesis is REJECTED
- alpha* appears to INCREASE (or stay constant) when properly normalized
- Implication: the computation lower bound does not weaken in high dimensions as naively expected

#### III.C: Composition Extension
- Two-step test (10,000 trials): composition bound tau_total >= (alpha/2)(C_1^2 + C_2^2) holds with huge margin
- The composition theorem and the computation lower bound are complementary: composition bounds tau from above, computation bound from below

**Key equations**:
- tau >= alpha(d) * ||[N(rho), N(sigma)]||_1^2
- alpha*(d=2) ~ 0.28 (random channels), ~ 0.001 (amplitude damping)

**Maps to**: Result A

### Section IV: Saturation Beyond Commutativity (Result B)

#### IV.A: The Corrected Saturation Condition
- Paper 1 sufficient conditions: (I) rho pure, (II) N(rho) = rho, (III) [rho, sigma] = 0, (IV) [rho, N(sigma)] = 0
- Key discovery: condition (III) [rho, sigma] = 0 is NOT necessary
- True necessary and sufficient condition: [N(rho), N(sigma)] = 0 (OUTPUT commutativity)
- Verified: for dephasing at theta = pi/4, saturation F^2 = exp(-DeltaD) holds for ALL sigma = diag(p, 1-p), even when [rho, sigma] != 0

#### IV.B: The Jensen Gap
- When [N(rho), N(sigma)] != 0, the bound F^2 >= exp(-DeltaD) is STRICTLY not saturated
- The gap F^2 - exp(-DeltaD) = A(s) * eps^2 + O(eps^4)
- A(s) = 2/sqrt(1-s^2) - arctanh(s)/s where s = eigenvalue gap of N(sigma)
- Power series proof: A(s) = sum_{n=0}^inf a_n s^{2n}, all a_n > 0 (via Wallis)
- Jensen gap J = <psi|sigma|psi> / exp(<psi|ln sigma|psi>) quantifies deviation from saturation
- J >= 1 always, J = 1 iff sigma is proportional to identity on support of psi

#### IV.C: The Operator Jensen Inequality Connection
- F^2 is an arithmetic (quadratic) mean of tau^{-1/4} omega tau^{-1/4}
- exp(-DeltaD) is a geometric (exponential) mean
- The gap is AM >= GM for operators, strict when the operator has non-degenerate spectrum
- Connection to Golden-Thompson inequality (structural, not direct application)
- This is a new quantum information inequality: operator Jensen gap as irreversibility witness

**Key equations**:
- F^2 = exp(-DeltaD) iff [N(rho), N(sigma)] = 0 (corrected saturation)
- J = <psi|sigma|psi> / exp(<psi|ln sigma|psi>) >= 1 (Jensen gap)
- A(s) = sum_{n>=0} [2*C(2n,n)/4^n - 1/(2n+1)] * s^{2n} > 0 (gap positivity)

**Maps to**: Result B

### Section V: Continuous-Time Dynamics (Result C)

#### V.A: Lindbladian Framework
- Build Lindbladian superoperator L in column-stacking vec representation
- Channel E_t = exp(Lt) as CPTP map for each time t
- tau(t) = 1 - F(rho_0, R_{sigma, E_t} o E_t(rho_0))
- This makes tau a continuous function of time -- the "clock" of irreversibility

#### V.B: Three Critical Benchmarks
- **Unitary evolution**: tau = 0 at 10^{-16} precision for all t in [0, 20*pi]
  - Physical: unitary with sigma = I/d => Petz recovery = U^{-1}, exact
  - This is the strongest numerical validation of tau = 0 <=> Sigma = 0
- **Pure dephasing**: tau(t) = (1 - exp(-4*gamma*t))/2
  - Factor 2 origin: dephasing is self-adjoint, Petz R = E_adj, so R o E = E^2
  - Verified to max|error| < 10^{-8} over 200 time points
- **Amplitude damping**: monotone non-decreasing tau(t), tau(0) = 0, tau(inf) -> 0.5
  - Verified: zero monotonicity violations across 300 time points

#### V.C: The Differential Equation
- Empirical: dtau/dt ~ (Sigma_dot/2)(1-tau) in the weak-noise regime
- Connection to Spohn's theorem: Sigma = -d/dt D(rho || rho_ss) >= 0 for Markov
- Physical interpretation: tau grows at a rate proportional to instantaneous entropy production, but self-limits as tau -> 1

#### V.D: Non-Markovian Regime -- The "Time Machine"
- In Markovian dynamics: tau(t) is monotone non-decreasing
- In non-Markovian dynamics: tau(t) can temporarily DECREASE (information backflow)
- dtau/dt < 0 corresponds to tau_signed < 0 -- local time reversal
- Jaynes-Cummings model verification: strong coupling (gamma_0/lambda > 2) produces oscillations in F(t) that witness non-Markovianity
- Connection to Paper 1's signed tau: tau_signed = ln(P_forward/P_reverse), negative means "time runs backward locally"
- NOT a violation of the second law: Sigma_total >= 0 still, but the subsystem can temporarily gain information

**Key equations**:
- tau(t)_dephasing = (1 - exp(-4*gamma*t))/2
- dtau/dt ~ (Sigma_dot/2)(1 - tau) (empirical differential equation)
- tau_signed = ln(P_forward / P_reverse), can be < 0 in non-Markovian regime

**Maps to**: Result C + "time machine" connection

### Section VI: Bayesian Composition Engine (Result D)

#### VI.A: The EDA Analogy
- Classical chip design: Static Timing Analysis propagates delays through circuit DAG
- Quantum analog: propagate sqrt(tau_i^eff) through the circuit, with Bayesian reference updating
- Key difference from multiplicative fidelity: each gate's cost depends on the evolved reference sigma_i = N_{i-1}(...(N_1(sigma_0)))

#### VI.B: Three-Stage Pipeline
- **Stage 1 (Saturation Scanner)**: Classify each gate as SATURATED (tau = 0, free) or NON-SATURATED
  - Uses condition [N(rho_evolved), N(sigma_evolved)] = 0 from Result B
  - Complexity O(n * d^3): eigendecomposition only, no decoder search
- **Stage 2 (Bayesian Composition)**: Propagate sqrt(tau_i^eff) with Bayesian-updated references
  - Feasibility check: sqrt(tau_total) <= sum sqrt(tau_i^eff) vs target tau
  - Budget allocation: proportional to sqrt-contribution
- **Stage 3 (Quantum Tax Calculator)**: Compare standard Petz (composable) vs rotated Petz (optimal per-gate)
  - Gap delta_i = tau_i - tau_tilde_i >= 0 quantifies cost of Bayesian consistency

#### VI.C: 5-Gate Benchmark
- Circuit: AD(0.05) -> Dep(0.08) -> AD(0.12) -> Deph(0.03) -> AD(0.10)
- rho = |+><+|, sigma_0 = diag(0.8, 0.2)
- Per-gate tau_eff vs tau_naive: reductions up to 94% (G3 dephasing: near-saturated due to evolved state)
- Total: 7.74% improvement over multiplicative fidelity
- Composition inequality: holds with positive slack
- Bayesian sigma evolution: sigma driven toward |0><0| by cumulative amplitude damping, making later dephasing gates nearly free

#### VI.D: Why Only the Standard Petz Map Works
- Rotated Petz (JRSWW) gives better per-gate fidelity but VIOLATES composition
- This is not a numerical accident -- it's structural: only the standard Petz satisfies functoriality (Bayesian consistency)
- Impossibility theorem (from composition_verification.py): no recovery map can simultaneously (i) satisfy composition for all channels, (ii) be at least as good as Petz everywhere, (iii) be strictly better somewhere
- The composability-optimality tradeoff is irreducible

**Key equations**:
- sqrt(tau_total) <= sum_i sqrt(tau_i^eff) (composition bound as feasibility criterion)
- sigma_i = N_i(sigma_{i-1}) (Bayesian reference updating)
- tau_i^eff = 1 - F(rho_i, R_{sigma_{i-1}, N_i}(N_i(rho_i))) (per-gate effective tau)

**Maps to**: Result D

### Section VII: The Unified Picture -- Bayesian Quantum Chronometry

#### VII.A: Three Pillars as Geometry
- Composition = triangle inequality in Bures space
- Saturation = the criterion for a point to sit on the "reversibility surface" (tau = 0)
- Gap (tau - tau_tilde) = the thickness of the triangle at each vertex
- Single diagram unifying all three

#### VII.B: The BQC Central Theorem
- Statement: for n-stage quantum process with Bayesian-updated references,
  sqrt(tau_total) <= sum sqrt(tau_i^eff), with:
  - Saturation criterion at each stage
  - Bayesian tax decomposition: sqrt(tau_total) <= sum sqrt(tau_tilde_i^eff) + Delta_Bayes
- The total arrow of time has two sources: optimal noise cost + Bayesian tax
- Bayesian tax is purely quantum (vanishes when [rho, sigma] = 0)

#### VII.C: Computation-Reversibility Paradox
- Saturation at every gate requires [omega_i, nu_i] = 0 at each stage -- effectively classical outputs
- But useful quantum computation requires creating superpositions -- generically [omega_i, nu_i] != 0
- Therefore: a perfectly reversible quantum computer is a classical computer
- The arrow of time is the PRICE of quantum computation, not a bug
- Quantified by: tau_computation >= sum_i f(C_i) from Result A

### Section VIII: Discussion and Outlook
- Comparison with existing tools: IBM AI Transpiler, Q-CTRL Fire Opal, Keysight True-Q
- tau-Chrono Engine provides guarantees + diagnostics that black-box tools cannot
- Experimental proposals: implement Bayesian composition on ion trap (Pino 2025) or NMR (Singh 2025)
- Open questions: general sigma saturation, optimal sigma_0 selection, multi-qubit DAG generalization
- Connection to Paper 2 (gravity): tau as the quantum origin of the metric
- Connection to non-Markovian time reversal: practical applications in quantum memory design

---

## 3. Key Equations (numbered as they will appear)

1. **tau definition**: tau = 1 - F(rho, R_{sigma,N}(N(rho)))
2. **Computation lower bound (A)**: tau >= alpha(d) * ||[N(rho), N(sigma)]||_1^2
3. **Corrected saturation (B)**: F^2 = exp(-DeltaD) iff [N(rho), N(sigma)] = 0 (for sigma = I/d, rho pure)
4. **Jensen gap (B)**: J = <psi|sigma|psi> / exp(<psi|ln sigma|psi>) >= 1
5. **Gap positivity (B)**: A(s) = 2/sqrt(1-s^2) - arctanh(s)/s = sum a_n s^{2n}, a_n > 0
6. **Dephasing tau (C)**: tau(t) = (1 - exp(-4*gamma*t))/2
7. **Differential equation (C)**: dtau/dt ~ (Sigma_dot/2)(1 - tau)
8. **Composition bound (D)**: sqrt(tau_total) <= sum_i sqrt(tau_i^eff)
9. **Bayesian reference updating (D)**: sigma_i = N_i(sigma_{i-1})
10. **BQC decomposition**: sqrt(tau_total) <= sum sqrt(tau_tilde_i^eff) + sum sqrt(delta_i)
11. **Impossibility theorem**: No R satisfies (i) composition + (ii) R >= Petz + (iii) R > Petz somewhere

---

## 4. Figure List

**Figure 1**: Computation lower bound tau vs C^2 scatter plot (d=2, 50k trials)
- Shows tau >= alpha * C^2 envelope
- Color-coded by channel family (AD, Dep, Deph, random)
- Inset: alpha*(d) vs d showing INCREASE, rejecting 1/d^2

**Figure 2**: Saturation gap analysis
- (a) F^2 vs exp(-DeltaD) for dephasing channel, sweep over sigma = diag(p, 1-p)
  - Shows F^2 = exp(-DeltaD) at theta = pi/4 for ALL p (saturation without input commutativity)
- (b) Jensen gap J vs p for different theta values
- (c) F^2 - exp(-DeltaD) vs eps^2 confirming quadratic scaling, different t1 values

**Figure 3**: Continuous-time tau(t) for three benchmarks
- (a) Dephasing: tau_numerical vs tau_analytical = (1-exp(-4gt))/2, residual < 10^{-8}
- (b) Unitary: tau(t) ~ 10^{-16} for all t (flat at machine epsilon)
- (c) Amplitude damping: monotone tau(t) approaching 0.5

**Figure 4**: Non-Markovian tau oscillations (bonus/outlook)
- Jaynes-Cummings model: tau(t) with oscillations (dF/dt > 0) witnessing backflow
- Comparison of Petz-based witness N_Petz vs BLP trace distance measure N_BLP

**Figure 5**: 5-gate Bayesian composition benchmark
- (a) Per-gate bar chart: tau_naive vs tau_eff, with classification labels
- (b) Bayesian sigma evolution: diagonal elements through 5 gates
- (c) Composition inequality: LHS = sqrt(tau_total) vs RHS = sum sqrt(tau_i^eff), showing slack

**Figure 6**: The BQC unified picture (schematic/diagram)
- Bures space triangle: rho, rho_1 = R_1(N_1(rho)), rho_12 = R_12(N_2(N_1(rho)))
- Edges labeled sqrt(tau_1), sqrt(tau_2^eff), sqrt(tau_12)
- Saturation points marked, gap annotations

---

## 5. Theorem Statements (Informal)

**Theorem 1 (Computation Lower Bound -- Conjecture)**:
For any CPTP map N, pure state rho, and full-rank reference sigma on a d-dimensional system,
tau(rho, N, sigma) >= alpha(d) * ||[N(rho), N(sigma)]||_1^2,
where alpha(d) > 0 is a dimension-dependent constant. Numerically: alpha(2) ~ 0.28 for generic channels.

**Theorem 2 (Corrected Saturation)**:
For sigma = I/d and pure rho = |psi><psi|, the JRSWW bound F^2 >= exp(-DeltaD) is saturated (equality) if and only if [N(rho), N(sigma)] = 0. The condition [rho, sigma] = 0 is sufficient but NOT necessary.

**Theorem 3 (Gap Positivity)**:
For d=2, if N(rho) = tau + eps*C where C is off-diagonal in the eigenbasis of N(sigma) = tau (i.e., [N(rho), N(sigma)] != 0), then F^2 - exp(-DeltaD) = A(s)*eps^2 + O(eps^4) where A(s) = sum_{n>=0} a_n s^{2n} with all a_n > 0. Hence A(s) >= 1 for all s in [0,1).

**Theorem 4 (Continuous tau for Dephasing)**:
For pure dephasing with rate gamma, initial state |+><+|, and reference sigma = I/2, the Petz recovery infidelity is exactly tau(t) = (1 - exp(-4*gamma*t))/2. The factor 4 (not 2) arises because the Petz map equals the adjoint channel and R o E = E^2.

**Theorem 5 (Bayesian Composition -- from Paper 1, applied)**:
For n sequential CPTP channels with Bayesian-updated references sigma_i = N_i(sigma_{i-1}),
sqrt(tau_total) <= sum_{i=1}^n sqrt(tau_i^eff),
where tau_i^eff is evaluated at the evolved state with the updated reference. This is strictly tighter than the multiplicative fidelity model for correlated noise chains.

**Theorem 6 (Impossibility -- Composability-Optimality)**:
No recovery map can simultaneously (i) satisfy composition sub-additivity for all channels, (ii) achieve fidelity >= the standard Petz map for all (N, rho, sigma), and (iii) achieve strictly higher fidelity for some (N, rho, sigma).

---

## 6. How A/B/C/D Map to Sections

| Result | Section | Core Content | Theorem |
|--------|---------|--------------|---------|
| **A**: Computation lower bound | III | tau >= alpha * C^2, dimension scaling, composition extension | Thm 1 |
| **B**: Saturation correction | IV | Output commutativity, Jensen gap, power series proof | Thm 2, 3 |
| **C**: Continuous tau(t) | V | Lindbladian, 3 benchmarks, differential equation, non-Markov | Thm 4 |
| **D**: Bayesian engine | VI | 3-stage pipeline, 5-gate benchmark, impossibility | Thm 5, 6 |
| Unified | VII | BQC central theorem, computation-reversibility paradox | -- |

---

## 7. Connection to "Time Machine" Concept

The non-Markovian regime (Section V.D) provides the physical connection to local time reversal:

1. **Markovian dynamics**: tau(t) monotone non-decreasing. The arrow of time only advances.
2. **Non-Markovian dynamics**: tau(t) can temporarily decrease. This means:
   - dF/dt > 0 at some t: the Petz recovery fidelity IMPROVES without intervention
   - Information flows back from environment to system
   - tau_signed = ln(P_forward/P_reverse) < 0: time locally "runs backward"

3. **Quantitative criterion**: dtau/dt < 0 is a **sufficient condition** for non-Markovianity. This gives a new non-Markovianity witness based on Petz recovery (distinct from BLP trace distance witness).

4. **Engineering implication**: In quantum memory design, engineering a non-Markovian environment can temporarily REVERSE irreversibility -- a controlled "time machine" at the quantum level. The tau-Chrono framework quantifies exactly how much reversal is achievable.

5. **Fundamental bound**: Even with non-Markovian backflow, the TOTAL entropy production Sigma_total >= 0 (second law). The temporary tau decrease is always "paid back" eventually. The maximum achievable backflow is bounded by the system-environment coupling structure.

6. **Connection to Paper 1 core thesis**: This is precisely Sheng-Kai's original insight from quantum erasure. In an open system, you must wait for the future to happen to verify retrocausality. In a perfectly closed system (Sigma = 0), future and present coexist. Non-Markovian backflow is the intermediate case: partial closure allows partial "time reversal."

---

## 8. Target Journal and Estimated Page Count

**Primary target**: Physical Review A (PRA)
- Rationale: Results are technical and quantitative, better suited to PRA's longer format than PRL
- PRA allows detailed proofs, extensive numerical evidence, and application sections
- Paper 1 went to PRL (the conceptual breakthrough); this paper is the quantitative follow-up

**Alternative targets**:
- Quantum (open access, good for tools/methods papers)
- New Journal of Physics (good for interdisciplinary QI + thermo + EDA)
- PRX Quantum (if we can emphasize the practical EDA implications strongly enough)

**Estimated page count**: 12-15 pages (PRA format, two-column)
- Sections I-II: 2 pages (intro + preliminaries)
- Section III (A): 2 pages
- Section IV (B): 2.5 pages (includes power series proof)
- Section V (C): 2 pages
- Section VI (D): 2.5 pages
- Section VII: 1.5 pages
- Section VIII: 1 page
- Figures: ~1.5 pages total
- References: ~1 page

**Supplemental material**: 8-10 pages
- Full power series proof of gap positivity (Sec IV extension)
- All 50,000-trial Monte Carlo data tables for Result A
- Benchmark code descriptions and reproducibility details
- Detailed comparison with IBM/Q-CTRL/Keysight

---

## 9. What Additional Results Are Needed Before Submission

### Must-Have (blocking submission)

1. **Analytical proof of Theorem 1 (computation lower bound)**
   - Currently only numerical evidence. Need either:
     - (a) A proof for d=2 using explicit Petz map formula + Fuchs-van de Graaf, OR
     - (b) A tighter statement with channel-dependent alpha(N) that can be proved, OR
     - (c) Demote to "Conjecture" with extensive numerical evidence (acceptable for PRA)
   - **Risk**: Medium. Option (c) is always available.

2. **General-sigma saturation theorem**
   - Result B currently proved for sigma = I/d. Need extension to general full-rank sigma.
   - The Bayesian engine (Result D) uses general sigma at each stage.
   - **Risk**: Low. The numerical evidence is overwhelming. Proof approach: reduce to output-space analysis where sigma only enters through N(sigma).

3. **Higher-dimensional verification of gap positivity (Theorem 3)**
   - Power series proof works for d=2. Need d=3,4,5 verification.
   - Numerical verification up to d=8 already done (all positive gaps).
   - **Risk**: Low. The operator Jensen inequality argument is dimension-independent.

4. **Multi-qubit benchmark for Result D**
   - Current 5-gate benchmark is single-qubit. Need at least a 2-qubit example.
   - Can use existing petz_toolkit.py for d=4.
   - **Risk**: Low (just engineering).

### Nice-to-Have (not blocking)

5. **Tight bound on alpha(d)**: Currently only lower bounds from sampling. An analytical characterization of the worst-case (minimizing) channel would strengthen Theorem 1.

6. **Experimental validation**: Run the 5-gate benchmark on IBM Brisbane or IonQ Aria. Not blocking for a theory paper but would dramatically strengthen the case.

7. **Non-Markovian Jaynes-Cummings quantitative benchmark**: The witness comparison (N_Petz vs N_BLP) currently has qualitative plots. Quantitative comparison (detection efficiency, false positive rate) would strengthen Section V.D.

8. **Qiskit integration**: Working demo of tau-Chrono as a Qiskit transpiler pass. Not needed for PRA submission but important for adoption.

---

## 10. Strategic Assessment: Strongest Result and Paper Structure

### Strongest Individual Result

**Result B (Corrected Saturation + Jensen Gap) is the most publishable standalone result.**

Reasons:
1. It CORRECTS a statement in Paper 1 (sufficient conditions for saturation are too strong)
2. The Jensen gap is a new quantum information inequality with a complete analytical proof (d=2)
3. The connection to operator Jensen / Golden-Thompson inequalities places it in a well-established mathematical framework
4. It has immediate implications: the saturation scanner in Result D depends on it

Runner-up: Result C (continuous-time tau) is clean and elegant, but less surprising -- it "should" work and it does. Result A is provocative but remains a conjecture. Result D is impressive engineering but lacks a sharp theorem.

### Single Unified Paper vs Separate Papers?

**Recommendation: Single unified paper (this outline).**

Arguments for unified:
- The four results interlock tightly: B feeds into D (saturation scanner), A provides the lower bound that D's composition bounds complement from above, C provides the continuous-time interpolation
- The BQC framework (Section VII) only makes sense with all four pillars
- A 12-15 page PRA paper accommodates all four comfortably
- The "engineering reversibility" narrative is stronger with all tools presented together

Arguments against (i.e., for splitting):
- Result B alone could be a 4-page PRL (correcting Paper 1's saturation condition + new inequality)
- Result D with the EDA framing could be a standalone tools paper for Quantum or PRX Quantum
- Splitting would produce more publications

**Verdict**: Write the unified paper first. If referees suggest splitting, Result B is the natural standalone candidate. The tau-Chrono Engine (Result D) could become a separate software paper if it gets Qiskit integration.

### Priority Order for Writing

1. Section IV (Result B) -- strongest, most complete, write first
2. Section V (Result C) -- clean, code already done, write second
3. Section VI (Result D) -- benchmark code exists, needs narrative
4. Section III (Result A) -- depends on whether we get analytical proof
5. Section VII (BQC unified) -- write last, after all pieces are in place

---

## Appendix: Code-to-Paper Mapping

| Code File | Paper Section | What It Provides |
|-----------|--------------|-----------------|
| `code/test_computation_lower_bound.py` | III (Result A) | 50k Monte Carlo, alpha* values |
| `code/test_computation_lower_bound_followup.py` | III.A | AD deep dive, boundary analysis |
| `code/test_scaling_analysis.py` | III.B | Dimension scaling alpha*(d) |
| `code/test_saturation_general_sigma.py` | IV (Result B) | General sigma verification, Jensen gap |
| `numerical/noncommuting_saturation_proof.py` | IV.B | Perturbative gap proof, power series |
| `numerical/saturation_verification.py` | IV (Result B) | Full saturation condition verification |
| `code/test_continuous_tau.py` | V (Result C) | 3 Lindbladian benchmarks |
| `simulations/non_markovian_witness.py` | V.D | JC model, Petz vs BLP witness |
| `code/bayesian_engine.py` | VI (Result D) | Bayesian composition core |
| `code/noise_channels.py` | VI (Result D) | Channel definitions |
| `code/benchmark_5gate.py` | VI.C | 5-gate benchmark |
| `numerical/composition_verification.py` | VI.D | Impossibility theorem verification |
| `numerical/petz_toolkit.py` | (all sections) | Core mathematical library |
