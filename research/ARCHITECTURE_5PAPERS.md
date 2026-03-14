# MASTER ARCHITECTURE: The 5-Paper tau Framework

**Author**: Sheng-Kai Huang
**Compiled**: 2026-03-11
**Status**: Master reference document for the entire research program

---

## THE ONE-PAGE PITCH

A single quantity -- tau = 1 - F, where F is the Petz recovery fidelity -- measures how irreversible a quantum channel is. When tau = 0, the past can be perfectly reconstructed from the present; when tau = 1, all information is lost. Paper 1 proves that this parameter unifies the arrow of time, quantum error correction, and thermodynamic entropy production through a four-way equivalence chain (tau=0 iff Sigma=0 iff quantum Markov chain iff perfect QEC). Paper 2 identifies the gravitational redshift as a Petz recovery bound -- the coordinate speed of light IS the recovery fidelity, c_eff/c = exp(-Sigma_grav/2) -- yielding an exponential metric with no event horizons and testable deviations from general relativity. Paper 1b recasts wavefunction collapse as Petz recovery failure, giving an explicit tau(t) curve that subsumes Penrose's gravitational collapse timescale. Paper 3 shows that the same Sigma, evaluated with the quantum-corrected gravitational channel at galactic scales (where Newton's G runs logarithmically), produces flat rotation curves without dark matter particles. Paper 5 formalizes the observation that tau is observer-dependent: different observers with access to different degrees of freedom see different arrows of time, connecting gravitational complementarity, quantum reference frames, and (speculatively) the structure of intelligence. The entire program converges on Paper 4's unified equation Sigma = D(rho_spacetime || rho_matter) -- one quantum relative entropy, solved under different boundary conditions, generating all known gravitational phenomenology from black holes to cosmic expansion.

---

## TABLE OF CONTENTS

1. [The Five Papers: Summary and Status](#1-the-five-papers)
2. [Dependency Map](#2-dependency-map)
3. [Optimal Publication Order](#3-optimal-publication-order)
4. [Complete Results Inventory](#4-complete-results-inventory)
5. [Open Problems by Severity](#5-open-problems-by-severity)
6. [Redundancies and Conflicts](#6-redundancies-and-conflicts)
7. [Master Timeline](#7-master-timeline)
8. [Critical Path Analysis](#8-critical-path-analysis)

---

## 1. THE FIVE PAPERS

### Paper 1: The Arrow of Time from Petz Recovery
**Status**: COMPLETE. Awaiting arXiv endorsement (Wilde contacted).
**Pages**: 6 + 22 (supplemental material)
**Target**: PRL
**GitHub**: https://github.com/akaiHuang/petz-recovery-unification
**File**: `paper/petz_recovery_unification.tex`

**Core results**:
- Definition: tau = 1 - F(rho, R_Petz(N(rho))) in [0,1]
- Equivalence chain: tau=0 iff Sigma=0 iff I(A;E|B)=0 iff QMC iff perfect QEC
- Bound: tau <= 1 - exp(-Sigma/2) [from JRSWW]
- Theorem 3: Saturation iff quantum sufficient statistic [NEW]
- Theorem 4: sqrt(tau) triangle inequality + composition sub-additivity [NEW]
- Quantum eraser = Petz recovery (explicit F=1 calculation) [NEW synthesis]
- Retrodiction Landauer Principle [NEW]
- 1620/1620 numerical verification points passed

**Validation**: Buscemi-Bai-Scarani (arXiv:2412.12489) independently constructed similar structure.

---

### Paper 1b: Collapse as Retrodiction Failure
**Status**: COMPLETE DRAFT. Needs final polish.
**Pages**: 7
**Target**: PRA
**File**: `paper/paper1b_collapse_retrodiction.tex`

**Core results**:
- Collapse := tau -> 1 (operational, within standard QM, no new physics)
- Gravitational dephasing channel: explicit Kraus operators, Petz map (R_Petz = N for dephasing)
- tau(t) = (1 - exp(-2Gamma*t^2))/2 [Pikovski mechanism]
- tau(t) = 1 - exp(-E_G*t/2hbar) [under Disi-Penrose assumption]
- Internal-spatial bridge: tau_internal = 2*tau_spatial*(1 - tau_spatial) [NEW]
- Penrose's tau_P = hbar/E_G is a special point on our continuous tau(t) curve
- Entanglement complementarity: Sigma_int = Sigma_spatial = S_entanglement

**Key limitation**: E_G is INPUT, not output. Cannot derive Disi-Penrose rate from first principles. Pikovski mechanism != Disi-Penrose mechanism in general.

---

### Paper 2: Gravity from tau
**Status**: DRAFT (9 pages). GRF Essay (4 pages) COMPLETE. Full version needs compression to 6 pages.
**Target**: PRD (full) + GRF Essay Competition (deadline 3/31/2026)
**Files**: `paper/paper2_gravity_tau.tex`, `paper/paper2_grf_essay.tex`

**Core results**:
- Sigma_grav = -ln(-g_00) [common mathematical core of 3 derivation routes]
- Triple identification: c_eff/c = exp(-Sigma_grav/2) = F_bound [NEW synthesis]
- For exponential metric: g_00 = -exp(-r_s/r), Sigma = r_s/r exactly
- Layer 1 (robust): tau < 1 everywhere -> no event horizons (model-independent)
- Layer 2 (conditional): exponential metric from approximate Petz saturation
- 5 testable predictions: shadow +4.6%, ISCO +5.6%, QNM ~4.4%, echo at ~4.2 r_s/c, NS redshift 19% different
- 3 first-principles derivation routes for Sigma = r_s/r:
  1. Modular flow / Dorau-Much extension: fractional QRE loss = r_s/r [EXACT]
  2. Gravitational Landauer (Herrera 2020): Tolman-modified erasure -> Sigma = r_s/r
  3. Quantum channel (pure loss bosonic): -ln(transmissivity) = r_s/r

**Key limitations**:
- No canonical gravitational CPTP map N_grav from first principles (the "channel problem")
- Petz bound saturation impossible for Sigma > 0; approximate saturation with gap O((r_s/r)^2)
- g_rr = exp(+r_s/r) requires isometry assumption g_00*g_rr = -1 (not derived from QI)

---

### Paper 3: Galactic Dynamics from Quantum Information Recovery
**Status**: RESEARCH STAGE. Outline exists, key resources identified, CMB problem acknowledged.
**Target**: PRL (letter) + Supplemental
**File**: `paper/paper3_rotation_curves_outline.tex`

**Core thesis**: The SAME Sigma = D(rho||sigma) - D(N(rho)||N(sigma)), evaluated with the full quantum-corrected gravitational channel N_grav(k) where G runs with scale, naturally gives:
- Small scales (r << r_c): Standard Newton (Sigma ~ r_s/r, as in Paper 2)
- Large scales (r >> r_c): Logarithmic correction giving flat rotation curves

**Key ingredients**:
- Kumar (2025, PLB 871): eta=1 marginal IR running -> Phi = -GM/r + (2GMk_*/pi)*ln(r/r_0)
- Gubitosi et al. (2024): SPARC 100-galaxy fits with running G, competitive with NFW
- Ghari & Haghi (2026): Verlinde favored over MOND at 5.2sigma in dwarf spheroidals
- Casini-Huerta theorem: QRE monotonicity under RG = DPI for coarse-graining
- Verlinde (2016): Volume-law entanglement -> apparent dark matter (complementary picture)

**New contributions planned**:
1. Sigma(r) = r_s/r + alpha*ln(r/r_c) as unified QI quantity bridging strong/weak field
2. DPI-motivated bound on anomalous dimension: eta <= 1 - exp(-Sigma/2) [CONJECTURE]
3. Reinterpretation: "dark matter" = extra information loss from IR running of gravitational channel
4. Connection to SPARC/RAR data

**FATAL PROBLEM**: CMB acoustic peaks CANNOT be explained by running G alone. Every CMB observable changes by 20-50% without CDM. The framework must either accept CDM particles, derive CDM-like fields from QRE, or connect to AeST-like theories. Paper 3 should adopt "honest agnosticism" (Option D) and focus on galactic-scale successes.

---

### Paper 4: Grand Unification -- Sigma = D(rho_spacetime || rho_matter)
**Status**: CONCEPTUAL. Key frameworks identified. No calculations.
**Target**: After Papers 1-3 accepted
**File**: `paper/paper4_grand_unification_outline.tex`

**Core equation**: Sigma = D(rho_spacetime || rho_matter)

One QRE, different boundary conditions -> different physics:
- Near a black hole: exponential metric (Paper 2)
- At galactic scales: flat rotation curves (Paper 3)
- At cosmological scales: accelerating expansion + Friedmann equations

**Key resources**:
- Bianconi (2025, PRD 111, 066001): Gravity action = QRE between spacetime and matter metrics
- Dorau-Much (2025, PRL): QRE -> Einstein equations
- Kodama vector for cosmological tau in FRW spacetimes
- Padmanabhan CosMIn: finite cosmic information requires Lambda > 0

**The single most important calculation**: Solve full non-linearized GfE field equations for static spherically symmetric vacuum. Does f(r) = exp(-r_s/r) emerge?

**Risk**: HIGH. Most ambitious paper. Depends on Papers 1-3 being accepted and validated.

---

### Paper 5: Observer-Dependent Temporal Asymmetry
**Status**: RESEARCH NOTES. Extensive theoretical development in `observer_dependent_tau.md` and `ai_tau_minimizer.md`.
**Target**: Foundations of Physics or similar
**Files**: `research/observer_dependent_tau.md`, `research/ai_tau_minimizer.md`

**Core insight**: tau = 1 - F is observer-dependent. Different observers with access to different subsystems assign different temporal asymmetries to the same physical process.

**Key results (developed but not yet in paper form)**:
- Observer hierarchy: tau_god = 0 <= tau_informed <= tau_ignorant <= 1
- Monotonicity: If O1's access is a subset of O2's, then tau_O1 >= tau_O2
- Quantum eraser re-derived: "erasing" = changing observer from partial to full access
- Gravitational complementarity: infalling vs. distant observer see different tau [connects to De Vuyst et al. 2025 JHEP]
- Complementary-observer uncertainty relation: tau_A + tau_B >= tau_min > 0 for complementary subsystems [CONJECTURED]
- Intelligence as tau-minimization: an agent that builds a predictive model effectively constructs the Petz recovery map, reducing tau toward zero [connects to Friston's Free Energy Principle]
- Three fundamental barriers to tau -> 0: chaos (Lyapunov), quantum (Heisenberg), computational irreducibility (Wolfram)

**Assessment**: The observer-dependence of tau is mathematically sound and supported by independent literature (Basso-Celeri PRL 2025, De Vuyst et al. JHEP 2025). The intelligence/AI connection is more speculative but provides a striking physical interpretation. This could be a standalone paper or the philosophical capstone of the series.

---

## 2. DEPENDENCY MAP

```
Paper 1 (tau definition, equivalence chain)
  |
  |--- Paper 1b (tau applied to collapse)
  |       |
  |       |--- [Penrose timescale subsumed]
  |       |--- [Pikovski decoherence = Petz problem]
  |
  |--- Paper 2 (tau applied to gravity)
  |       |
  |       |--- Requires: Sigma_grav = r_s/r (3 routes)
  |       |--- Requires: Approximate Petz saturation
  |       |--- Does NOT require: Paper 1b
  |       |
  |       |--- Paper 3 (tau at galactic scales)
  |       |       |
  |       |       |--- Requires: Paper 2's Sigma = -ln(g_00)
  |       |       |--- Requires: Kumar (2025) running G
  |       |       |--- Requires: Casini-Huerta RG = DPI
  |       |       |--- Does NOT require: Paper 1b
  |       |       |--- BLOCKED by: CMB problem (acknowledged, not fatal)
  |       |       |
  |       |       |--- Paper 4 (grand unification)
  |       |               |
  |       |               |--- Requires: Papers 2+3 validated
  |       |               |--- Requires: Bianconi GfE or equivalent
  |       |               |--- Requires: Cosmological Sigma via Kodama
  |       |               |--- BLOCKED by: GfE spherical solution (not computed)
  |       |
  |       |--- Paper 5 (observer-dependence)
  |               |
  |               |--- Requires: Paper 1's tau definition
  |               |--- Uses: Paper 2's gravitational complementarity
  |               |--- Does NOT require: Papers 3 or 4
  |               |--- Independent enough to publish in parallel
  |
  |--- Viewpoint Paper ("Why You Must Wait")
          |--- Requires: Paper 1
          |--- Target: Foundations of Physics
          |--- Status: Complete draft (23 pages)
```

### Key dependencies:
- **Paper 1 is foundational**: ALL other papers depend on it
- **Paper 1b is a sidebar**: Important but not on the critical path for Papers 2-4
- **Paper 2 depends only on Paper 1**: Can proceed independently
- **Paper 3 depends on Paper 2**: Needs Sigma = -ln(g_00) established
- **Paper 4 depends on Papers 2+3**: Cannot start calculations until earlier results validated
- **Paper 5 depends on Paper 1 (minimally Paper 2)**: Can be parallelized

---

## 3. OPTIMAL PUBLICATION ORDER

### Current numbering: 1 -> 1b -> 2 -> 3 -> 4 -> 5
### Recommended numbering: KEEP AS IS, with parallelization

**Rationale**: The logical chain (define tau -> apply to gravity -> extend to galaxies -> unify) is correct. Paper 5 (observer-dependence) should be written in parallel with Papers 3-4, not sequentially after them.

**Recommended publication timeline**:

```
Phase 1 (NOW - Q2 2026):
  [1] Paper 1 on arXiv (awaiting endorsement)
  [2] Paper 2 GRF Essay submitted (deadline 3/31/2026)
  [3] Paper 1b polished and submitted to PRA

Phase 2 (Q2-Q3 2026):
  [4] Paper 2 full version (6 pages) submitted to PRD
  [5] Viewpoint paper submitted to Foundations of Physics
  [6] Paper 5 draft started (can parallel with everything)

Phase 3 (Q3-Q4 2026):
  [7] Paper 3 calculations (SPARC fitting, DPI bound on eta)
  [8] Paper 3 submitted to PRL
  [9] Paper 5 submitted to Foundations of Physics

Phase 4 (2027+):
  [10] Paper 4 calculations (GfE solutions, cosmological Sigma)
  [11] Paper 4 submitted
```

### What can be parallelized:
- Paper 1b and Paper 2 GRF Essay: INDEPENDENT, do simultaneously
- Paper 5 and Paper 3: INDEPENDENT, can overlap
- Paper 2 full and Viewpoint: INDEPENDENT

### What CANNOT be parallelized:
- Paper 3 needs Paper 2 accepted (or at least on arXiv) for credibility
- Paper 4 needs Papers 2+3 validated

---

## 4. COMPLETE RESULTS INVENTORY

### 4.1 Proven Mathematical Theorems

| # | Result | Paper | Status |
|---|--------|-------|--------|
| T1 | tau=0 iff Sigma=0 iff I(A;E|B)=0 iff QMC iff perfect QEC | P1 | PROVED (new synthesis) |
| T2 | tau <= 1 - exp(-Sigma/2) | P1 | PROVED (restated JRSWW) |
| T3 | Saturation iff quantum sufficient statistic | P1 SM | PROVED (new theorem) |
| T4 | sqrt(tau) triangle inequality + tau_12 <= tau_1 + tau_2 - tau_1*tau_2 | P1 SM | PROVED (new theorem) |
| T5 | Unitary channel -> F=1, tau=0, sigma-independent | P1 | PROVED (corollary) |
| T6 | Quantum eraser F=1 from Petz | P1 | CALCULATED (new) |
| T7 | Gravitational dephasing: R_Petz = N (self-adjoint) | P1b | CALCULATED |
| T8 | tau(t) = (1-exp(-2Gamma*t^2))/2 for Pikovski | P1b | DERIVED |
| T9 | tau_internal = 2*tau_spatial*(1-tau_spatial) | P1b | DERIVED (new) |
| T10 | NEC violation => dSigma_Raych/ds < 0 (tau can decrease) | P2 | PROVED with caveats |
| T11 | Approximate Petz saturation: gap ~ O((r_s/r)^2) in weak field | P2 | PROVED |

### 4.2 New Identifications / Syntheses

| # | Result | Paper | Status |
|---|--------|-------|--------|
| S1 | Retrodiction Landauer Principle: W_retro >= k_B*T*Sigma/2 | P1 | NEW principle |
| S2 | c_eff/c = exp(-Sigma_grav/2) = F_bound | P2 | NEW triple identification |
| S3 | Sigma_grav = -ln(-g_00) = r_s/r [3 routes] | P2 | NEW identification (motivated, not proved) |
| S4 | n(r) = 1/F_bound = exp(+Sigma/2) (refractive index) | P2 | NEW identification |
| S5 | tau=1 at apparent horizon (not event horizon) | P2 | NEW (connects to QES) |
| S6 | Penrose's E_G -> Sigma mapping; tau(t_Penrose) = 0.39 | P1b | NEW subsumption |
| S7 | Collapse = AND->OR via Petz recovery failure | P1b | NEW operational definition |
| S8 | Sigma(r) = r_s/r + alpha*ln(r/r_c) bridging strong/weak field | P3 | PLANNED |
| S9 | "Dark matter" = extra information loss from IR running | P3 | PLANNED interpretation |
| S10 | tau is observer-dependent: tau_O depends on O's access set | P5 | DEVELOPED (in notes) |
| S11 | Intelligence = tau-minimization capacity | P5 | SPECULATIVE |
| S12 | Sigma_eff = Sigma_grav - I(L;R) for entanglement-assisted channels | P2 | HEURISTIC |

### 4.3 Testable Predictions

| # | Prediction | Paper | Distinguishes from | Testable with |
|---|------------|-------|--------------------|---------------|
| P1 | F >= exp(-Sigma/2) for all quantum channels | P1 | Nothing (bound) | Png/Singh data reanalysis |
| P2 | Collapse is continuous tau(t), not step-function | P1b | Penrose (step) | Mesoscopic superposition experiments |
| P3 | tau(t_Penrose) = 0.39 (not 0 or 1) | P1b | Penrose (tau=1) | Same as P2 |
| P4 | No event horizons (tau < 1 everywhere) | P2 L1 | GR (Schwarzschild) | Theoretical distinction |
| P5 | Black hole shadow 4.6% larger than GR | P2 L2 | Schwarzschild | EHT (future precision) |
| P6 | ISCO 5.6% larger | P2 L2 | Schwarzschild | X-ray timing |
| P7 | QNM frequency shifted ~4.4% | P2 L2 | Schwarzschild | LIGO/Virgo/KAGRA |
| P8 | GW echoes at Delta_t ~ 4.2 r_s/c | P2 L2 | Planck-wall models (89x longer) | LIGO O5+ |
| P9 | Neutron star redshift 19% different | P2 L2 | Schwarzschild | NICER |
| P10 | Flat rotation curves from running G | P3 | CDM halos | SPARC, BIG-SPARC |
| P11 | Universal k_* ~ 2.7e-2 kpc^-1 | P3 | Per-galaxy CDM fits | SPARC statistical analysis |

### 4.4 Numerical Verifications

| Verification | Result | Paper |
|-------------|--------|-------|
| tau-R bridge (multiple channels) | PASSED | P1 |
| SBS formation dynamics | PASSED | P1 |
| N-scaling bound | PASSED | P1 |
| 5 genuinely new predictions | ALL PASSED | P1 |
| Saturation theorem (thousands of random channels) | PASSED | P1 |
| Table I: Earth -> event horizon values | COMPUTED | P2 GRF |

---

## 5. OPEN PROBLEMS BY SEVERITY

### FATAL (would invalidate core claims if proven wrong)

| Problem | Affects | Status | Mitigation |
|---------|---------|--------|------------|
| **None identified** | - | - | The framework is built on proven QI theorems; the extensions are conjectural but the foundation is solid |

### SERIOUS (major obstacles to specific papers)

| # | Problem | Affects | Severity | Status |
|---|---------|---------|----------|--------|
| G1 | **No canonical gravitational CPTP map** N_grav with Delta_D = -ln(-g_00) | P2 | SERIOUS | OPEN. 5+ candidate channels surveyed, none exact. Best: Trejo-Calderon modular channel |
| G2 | **CMB acoustic peaks impossible without CDM** from running G alone | P3 | SERIOUS | ACKNOWLEDGED. 20-50% discrepancy on every observable. Must adopt "honest agnosticism" |
| G3 | **Petz saturation structurally impossible** for Sigma > 0 (Golden-Thompson inequality) | P2 | SERIOUS | MITIGATED. Approximate saturation with O((r_s/r)^2) gap suffices for weak-field predictions |
| G4 | **Born rule / single outcome not explained** by tau | P1b | SERIOUS | ACKNOWLEDGED honestly. tau explains PROCESS (AND->OR) but not SELECTION (which outcome) |
| G5 | **g_rr not derivable from QI** -- requires isometry assumption | P2 | SERIOUS | OPEN. Bisognano-Wichmann spatial modular flow is most promising unexplored route |
| G6 | **GfE spherically symmetric solution not computed** | P4 | SERIOUS | BLOCKING Paper 4. Feasible calculation (numerical ODE) but nobody has done it |

### MODERATE (limit scope but don't invalidate)

| # | Problem | Affects | Status |
|---|---------|---------|--------|
| M1 | Pikovski mechanism != Disi-Penrose mechanism | P1b | ACKNOWLEDGED. Match only under specific conditions |
| M2 | Exact Sigma in strong field: r_s/r vs -ln(1-r_s/r)? | P2 | TESTABLE prediction (exponential vs Schwarzschild) |
| M3 | Bullet Cluster has no mechanism in tau framework | P3 | OPEN. Requires non-linear structure formation |
| M4 | Sigma_eff = Sigma_grav - I(L;R) is heuristic, not derived | P2 | OPEN. Needs verification in SYK or 2d CFT |
| M5 | eta=1 for IR running not derived from QRE first principles | P3 | OPEN. Most promising: DPI bound on eta via Casini-Huerta |
| M6 | tau_signed < 0 qualitative prediction already known from HP/ER=EPR | P2 | The QUANTITATIVE formula is new; qualitative is not |
| M7 | Observer-dependent tau: complementary-observer uncertainty relation unproven | P5 | CONJECTURED. tau_A + tau_B >= tau_min |
| M8 | Intelligence = tau-minimization: three fundamental barriers not fully quantified | P5 | DEVELOPED but formulas approximate |

### MINOR (cosmetic or presentation issues)

| # | Problem | Affects | Status |
|---|---------|---------|--------|
| m1 | Paper 2 needs compression from 9 to 6 pages | P2 | PLANNED (detailed revision plan exists) |
| m2 | Exponential metric known since Yilmaz (1958) -- must be clear about what is new vs old | P2 | NOTED in novelty analysis |
| m3 | tau = 1-F definition may be dismissed as "just notation" | P1 | MITIGATED by emphasizing equivalence chain + theorems |
| m4 | Viewpoint paper (23pp) may be too long for target journal | Viewpoint | Can be shortened |

---

## 6. REDUNDANCIES AND CONFLICTS

### 6.1 Internal Redundancies (consolidation needed)

| Redundancy | Files | Action |
|------------|-------|--------|
| Paper 2 revision plan exists in BOTH revision_plan.md AND the tex file comments | paper2_revision_plan.md, paper2_gravity_tau.tex | Use revision_plan.md as canonical |
| Three separate "derivation route" files overlap with paper2_first_principles_integration.md | derivation_route{1,2,3}.md, paper2_first_principles_integration.md | Integration report supersedes individual routes |
| Paper 3 rotation curves discussed in 4+ files | paper3_rotation_curves.md, paper3_rotation_curves_deep_research.md, paper3_scale_dependent_sigma.md, paper2_to_paper3_bridge.md | Deep research file is canonical |
| Observer-dependent tau discussed in observer_dependent_tau.md AND unified_intersection.md AND ai_tau_minimizer.md | Three files | observer_dependent_tau.md is canonical for Paper 5 |

### 6.2 Internal Conflicts (resolved or flagged)

| Conflict | Resolution |
|----------|------------|
| Paper 2 line 225-226 conflated Petz sufficiency with bound saturation | FIXED (2026-03-07). Approximate saturation acknowledged |
| Sigma_grav = r_s/r (exponential) vs -ln(1-r_s/r) (Schwarzschild) | NOT a conflict: two different metrics, agree in weak field, differ at O((r_s/r)^2). Nature decides -- TESTABLE |
| Pikovski gives tau ~ 10^-40 at galactic scales, too small for anything | NOT a conflict with Paper 3: Paper 3 uses running G mechanism, not Pikovski decoherence |
| Paper 1b claims Penrose "subsumed" but E_G is input not output | HONEST framing: "Penrose's timescale is reproduced under DP assumption; the framework does not derive E_G from first principles" |
| Bianconi GQRE uses R_s/r^3 (curvature), our Sigma uses R_s/r (potential) | REAL difference. GQRE and Sigma_grav are mathematically distinct. Relation = integration over path, but not exact identification. Paper 4 must address this |

### 6.3 External Conflicts

**No external contradictions found across 33+ papers surveyed (as of 2026-03-09).**

Closest potential tensions:
- Bai-Buscemi-Scarani (2024): parallel construction, but they do not connect to gravity or QEC
- Donoghue critique of asymptotic safety (2020): challenges the IR running of G used in Paper 3. Mitigated by using Kumar's EFT argument (not AS per se)
- Gran Sasso (2021): rules out parameter-free Disi model. Penrose version with adjusted smearing still viable

---

## 7. MASTER TIMELINE

```
2026-03 -------- NOW --------
  [x] Paper 1 on GitHub
  [x] Paper 1b draft complete (7pp)
  [x] Paper 2 GRF Essay draft complete (4pp)
  [x] Viewpoint paper draft complete (23pp)
  [x] Paper 2 full draft exists (9pp)
  [ ] Paper 1 arXiv endorsement        <- PRIORITY 1
  [ ] Paper 2 GRF Essay submission      <- PRIORITY 2 (deadline 3/31)

2026-04
  [ ] Paper 1b final polish + submission to PRA
  [ ] Paper 1 fix branch merge + push

2026-05-06
  [ ] Paper 2 full version: compress 9->6pp
  [ ] Paper 2 full submission to PRD
  [ ] Viewpoint submission to Foundations of Physics
  [ ] Begin Paper 5 draft (observer-dependent tau)
  [ ] Begin Paper 3 calculations (SPARC fitting with Sigma(r))

2026-07-09
  [ ] Png/Singh experimental data reanalysis proposal
  [ ] Paper 3 core calculations complete
  [ ] Paper 3 submission to PRL
  [ ] Paper 5 submission

2026-10+
  [ ] Paper 4 calculations begin (GfE solutions)
  [ ] Respond to referee reports on Papers 1, 1b, 2

2027+
  [ ] Paper 4 submission
  [ ] Experimental verification results (Png, Singh, or others)
```

---

## 8. CRITICAL PATH ANALYSIS

### The bottleneck: Paper 1 on arXiv

Everything downstream depends on Paper 1 being publicly available. Without arXiv posting:
- Paper 2 cannot cite it properly
- Experimentalists (Png, Singh) cannot see it to reanalyze their data
- The community cannot verify independent construction (Buscemi et al.)
- All subsequent papers lack their foundation reference

**Action**: Follow up with Wilde on endorsement. Consider alternative endorsers (Buscemi, Scarani, Parzygnat).

### The second bottleneck: Paper 2 GRF Essay deadline (3/31/2026)

This is a hard deadline. The essay is complete but needs:
- Final compilation check
- Verify all equations
- Submit

### The third bottleneck: The channel problem (G1)

The absence of a canonical gravitational CPTP map is Paper 2's Achilles heel. Without it, the triple identification c_eff/c = F_bound remains "motivated conjecture" rather than theorem. The three derivation routes provide strong circumstantial evidence, but a skeptical referee will press hard.

**Most promising resolution**: Trejo-Calderon modular thermal filter channel. This would need an explicit construction showing Delta_D = -ln(-g_00) for a well-defined Hilbert space and reference state.

### The fourth bottleneck: CMB problem (G2)

This blocks Paper 3 from being a complete "dark matter replacement" paper. The recommended strategy is honest agnosticism: Paper 3 focuses on galactic-scale success, acknowledges CMB limitation, points to AeST as proof that MOND-like theories CAN fit CMB with the right field content, and defers the full cosmological question to Paper 4.

---

## APPENDIX A: NOVELTY SCORECARD (from novelty_analysis.md)

| Claim | Category | Novelty | Significance |
|-------|----------|---------|-------------|
| Four-way equivalence chain | New synthesis (C) | **HIGH** | **HIGH** |
| Saturation theorem | New math (E) | **HIGH** | Moderate |
| sqrt(tau) triangle inequality | New math (E) | **HIGH** | Moderate |
| c_eff/c = F_bound | New synthesis (C/D) | **HIGH** | **HIGH** (if channel problem solved) |
| GW echo prediction at 4.2 r_s/c | New prediction (F) | Moderate | **HIGH** (testable) |
| tau = 1 at apparent horizon | New synthesis (C/D) | Moderate-High | Moderate-High |
| tau(t) collapse dynamics | New formula (E) | Moderate | Moderate |
| Retrodiction Landauer Principle | New synthesis (C/D) | Moderate | Moderate |
| Sigma_grav = r_s/r | New notation (B/D) | Low | Low-Moderate |
| Exponential metric from tau | New interpretation (D) | Low | Moderate |

**Overall**: ~10-15% genuinely new content, ~25% new synthesis, ~25% new interpretation, ~35-40% known results in new notation. Comparable to successful conceptual physics papers.

---

## APPENDIX B: COMPLETE FILE INDEX

### Paper manuscripts
| File | Content | Status |
|------|---------|--------|
| `paper/petz_recovery_unification.tex` | Paper 1 main | Complete |
| `paper/petz_recovery_supplemental.tex` | Paper 1 supplemental | Complete |
| `paper/paper1b_collapse_retrodiction.tex` | Paper 1b | Draft |
| `paper/paper2_gravity_tau.tex` | Paper 2 full | Draft (9pp, needs compression) |
| `paper/paper2_grf_essay.tex` | Paper 2 GRF essay | Complete (4pp) |
| `paper/paper2_supplement.tex` | Paper 2 supplemental | Draft |
| `paper/viewpoint_temporal_coexistence.tex` | Viewpoint paper | Complete draft |
| `paper/paper3_rotation_curves_outline.tex` | Paper 3 outline | Outline only |
| `paper/paper4_grand_unification_outline.tex` | Paper 4 outline | Outline only |

### Key research files (canonical versions)
| File | Content | Canonical for |
|------|---------|---------------|
| `research/paper2_first_principles_integration.md` | Sigma = r_s/r derivations | Paper 2 foundations |
| `research/paper1b_penrose_collapse_from_tau.md` | Collapse = Petz failure | Paper 1b foundations |
| `research/paper3_rotation_curves_deep_research.md` | Running G + SPARC + DPI | Paper 3 foundations |
| `research/paper3_CMB_without_DM.md` | CMB problem analysis | Paper 3 limitations |
| `research/paper4_bianconi_gfe_research.md` | Bianconi GfE analysis | Paper 4 foundations |
| `research/observer_dependent_tau.md` | Observer-dependent tau | Paper 5 foundations |
| `research/ai_tau_minimizer.md` | Intelligence = tau minimization | Paper 5 (speculative) |
| `research/paper2_canonical_channel.md` | Channel problem survey | Paper 2 gap G1 |
| `research/paper2_tau_negative_synthesis.md` | tau < 0 analysis | Paper 2 discussion |
| `research/paper2_to_paper3_bridge.md` | Bridge P2->P3 | Papers 2-3 connection |
| `research/novelty_analysis.md` | Honest novelty assessment | All papers |
| `research/mathematical_proofs.md` | Formal proof verification | All papers |
| `research/MAXWELL_COMPARISON.md` | Framework comparison | Motivation/context |
| `research/paper2_petz_gravity_latest.md` | Latest literature (2024-2026) | Paper 2 references |
| `research/unified_intersection.md` | Where all paths meet | Papers 4-5 |

---

## APPENDIX C: KEY EXTERNAL REFERENCES (MUST-CITE)

### Foundational (Papers 1-2)
1. Petz (1986, 1988): Recovery map, sufficiency theorem
2. Fawzi-Renner (2015): F >= exp(-Sigma/2)
3. JRSWW (2018): Universal recovery
4. Parzygnat-Buscemi (2023): Petz = unique retrodiction functor
5. Crooks (1999): Fluctuation theorem
6. Pikovski et al. (2015): Gravitational decoherence
7. Penrose (1996): Gravitational collapse timescale
8. Jacobson (1995): Thermodynamic Einstein equations

### Critical for Paper 2
9. Dorau-Much (2025, PRL): QRE -> Einstein equations
10. Herrera (2020): Gravitational Landauer principle
11. Basso-Celeri (2025, PRL): Observer-dependent Sigma in curved spacetime
12. Bianconi (2025, PRD): Gravity action = QRE
13. Yilmaz (1958) / Papapetrou: Exponential metric

### Critical for Paper 3
14. Kumar (2025, PLB): eta=1 marginal running -> rotation curves
15. Gubitosi et al. (2024): SPARC 100-galaxy running G fits
16. Casini-Huerta (2017, JHEP): QRE monotonicity = RG irreversibility
17. Verlinde (2016): Emergent gravity and dark matter
18. Ghari-Haghi (2026): Verlinde > MOND at 5.2sigma
19. Skordis-Zlosnik (2021, PRL): AeST fits CMB without CDM particles

### Independent verification
20. Bai-Buscemi-Scarani (2024, arXiv:2412.12489): Parallel construction

---

## APPENDIX D: THE PHILOSOPHY

> "We have not discovered new things. We have organized the nature of known things, so everyone can stand on this foundation to explore the boundary."
> -- Sheng-Kai Huang

The framework's genuine contribution is SYNTHESIS, not discovery. It provides a unified language (tau) that connects quantum information recovery, thermodynamic irreversibility, quantum error correction, gravitational physics, and (tentatively) cosmology through a single parameter. This is the physics equivalent of a Rosetta Stone -- translating between communities that speak different mathematical dialects about the same underlying structure.

The risk is that a referee dismisses it as "just notation." The defense is threefold:
1. The mathematical theorems (T3, T4) are provably new
2. The testable predictions (P4-P9) are quantitatively specific
3. The four-way equivalence chain (T1) has never been assembled, despite each link being known

---

*This document supersedes VISION_4PAPERS.md (which described 4 papers) and PROGRESS_REPORT_2026_03_09.md (which was a snapshot). It should be updated whenever a paper's status changes or a new problem/result is identified.*
