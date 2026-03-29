# Outreach Emails — tau-chrono

## 1. Q-CTRL (Michael Biercuk / André Carvalho)

**To**: michael.biercuk@q-ctrl.com or support@q-ctrl.com
**Subject**: Complementary tool for Fire Opal — physics-based circuit fidelity prediction

Hi Michael,

I'm Sheng-Kai Huang, an independent quantum physics researcher in Taiwan working on Petz recovery maps applied to circuit-level noise prediction.

I built tau-chrono, an open-source tool that predicts circuit fidelity *before* execution using Bayesian composition of the Petz recovery map. On QuTech's Tuna hardware, it predicts noise 46% more accurately than independent gate models at depth 50.

I think tau-chrono and Fire Opal are naturally complementary:

- Fire Opal suppresses errors during execution (the best treatment)
- tau-chrono predicts which circuits are worth running before execution (the best diagnosis)
- Together: diagnose first, then treat — fewer wasted QPU cycles, better outcomes

The key difference from ML-based predictors: tau-chrono requires zero training data (just per-gate tomography), works across hardware platforms, and provides mathematically provable bounds via the Petz recovery inequality.

One line to try:
  pip install tau-chrono
  python -c "from tau_chrono.demo import quick_demo; quick_demo()"

GitHub: https://github.com/akaiHuang/petz-recovery-unification

I'm not selling anything — the tool is MIT-licensed. I'd love to explore whether integrating pre-execution prediction with Fire Opal's runtime optimization could benefit your customers.

Best,
Sheng-Kai Huang
akai@fawstudio.com
Taiwan

---

## 2. QuTech / Quantum Inspire Team

**To**: support@quantum-inspire.com or info@qutech.nl
**Subject**: 46% noise prediction improvement on your Tuna hardware

Hi Quantum Inspire team,

I'm Sheng-Kai Huang, an independent researcher using your platform. I wanted to share some results from experiments on your Tuna superconducting processor.

Using a tool I built called tau-chrono (Bayesian noise tracking via Petz recovery maps), I measured the following improvements over standard independent gate noise models:

  Depth 10: +13% more accurate prediction
  Depth 20: +30%
  Depth 30: +38%
  Depth 50: +46%

The composition inequality (a mathematical guarantee from Petz recovery theory) holds in all cases. All data was collected on your hardware.

tau-chrono is open-source (MIT license):
  GitHub: https://github.com/akaiHuang/petz-recovery-unification

I'd be happy to share the full dataset or collaborate on a benchmark report. These results suggest your hardware is more capable than standard noise models indicate — which could be valuable for your users and for showcasing the platform.

Best regards,
Sheng-Kai Huang
akai@fawstudio.com

---

## 3. Nathan Shammah (Unitary Fund CTO)

**To**: nathan@unitary.fund
**Subject**: Microgrant application — tau-chrono: Petz recovery for circuit noise prediction

Hi Nathan,

I'm Sheng-Kai Huang, an independent researcher in Taiwan. I'd like to apply for a Unitary Fund microgrant for tau-chrono, an open-source quantum noise prediction tool.

What it does: Uses Petz recovery maps + Bayesian composition to predict circuit fidelity before execution. Unlike ML-based approaches (which need thousands of training circuits), tau-chrono requires only per-gate process tomography and provides mathematically provable bounds.

Results on real hardware (QuTech Tuna superconducting):
  - 46% more accurate noise prediction at depth 50
  - Composition inequality holds at all depths
  - All depths show positive Bayesian improvement

What makes it different from pyGSTi/Mitiq/Fire Opal:
  - Predicts BEFORE execution (not post-hoc mitigation)
  - Zero training data (not ML)
  - Physics-based with provable bounds (Petz recovery)
  - 169 tests passing, pip-installable

Current state: v0.2.0 on GitHub, validated on real hardware, 7,600+ lines of code.

  pip install tau-chrono
  GitHub: https://github.com/akaiHuang/petz-recovery-unification

The grant would support: IBM Quantum validation, documentation, and community building.

Best,
Sheng-Kai Huang
akai@fawstudio.com

---

## 4. Naoki Kanazawa (IBM Quantum, Qiskit Experiments)

**To**: knzwnao@jp.ibm.com (or via GitHub/LinkedIn)
**Subject**: Bayesian noise tracking complements Qiskit noise learning

Hi Naoki,

I've been following your work on Qiskit Experiments and noise characterization. I built tau-chrono, an open-source tool that uses Petz recovery maps for circuit-level noise prediction — a complementary approach to IBM's sparse Pauli-Lindblad noise learning.

Key difference: tau-chrono propagates a Bayesian reference state (sigma) through the circuit gate-by-gate, giving each gate's Petz recovery map the correct prior. This captures the noise "saturation" effect that independent models miss.

On QuTech's Tuna hardware:
  - Depth 30: 38% more accurate than independent model
  - Depth 50: 46% more accurate
  - Composition inequality (provable bound) holds at all depths

I think this could complement Qiskit's noise learning pipeline — the Pauli-Lindblad model captures the per-layer noise structure, while Bayesian Petz tracking captures how noise compounds differently at different circuit depths.

  pip install tau-chrono
  GitHub: https://github.com/akaiHuang/petz-recovery-unification

Would love to discuss or validate on IBM hardware.

Best,
Sheng-Kai Huang
akai@fawstudio.com

---

## 5. Kenneth Brown Lab (Duke University)

**To**: devansh.bhardwaj@duke.edu (or kenneth.r.brown@duke.edu)
**Subject**: Petz recovery meets drifting noise estimation

Hi Devansh / Prof. Brown,

I read your work on adaptive estimation of drifting noise in QEC — really interesting approach to a critical problem.

I've been working on a related problem from the Petz recovery side: tracking noise evolution through circuits using Bayesian composition of the reference state. The result is tau-chrono, an open-source tool that predicts circuit fidelity before execution.

The connection to your work: your approach tracks noise drift over *time* (between calibrations). tau-chrono tracks noise accumulation over *depth* (through the circuit). Combining both could give a complete picture: time-varying noise parameters (your method) fed into depth-aware Bayesian composition (our method).

Hardware results (QuTech Tuna, superconducting):
  - 46% more accurate at depth 50
  - Provable composition inequality from Petz recovery bound
  - All depths positive improvement

  pip install tau-chrono
  GitHub: https://github.com/akaiHuang/petz-recovery-unification

Would love your thoughts on potential synergies.

Best,
Sheng-Kai Huang
akai@fawstudio.com

---

## 6. Guillaume Verdon (Extropic AI)

**To**: guillaume@extropic.ai (or via X/Twitter DM)
**Subject**: Petz recovery bounds for thermodynamic computing efficiency

Hi Guillaume,

I'm an independent researcher working on Petz recovery maps and quantum relative entropy. I noticed that the core math directly applies to thermodynamic computing.

The Petz recovery bound F >= exp(-Sigma/2) quantifies the irreversibility cost of any stochastic process — including Langevin dynamics on your TSU hardware. In the classical limit, it reduces to:

  Sigma_irr = D_KL(p||q) - D_KL(N(p)||N(q)) >= -2 log BC(p, R_Bayes∘N(p))

This gives a tight efficiency bound for thermodynamic sampling:
  - eta = Delta_D_useful / Sigma_total ∈ [0,1]
  - Convergence certificate: F² >= exp(-D₀(1 - e^{-2αt}))

I've built tau-chrono, a tool that computes these bounds for quantum circuits (46% accuracy improvement on real hardware). Adapting it for THRML would give Extropic developers:
  - Quantify energy waste per computation step
  - Certify sampler convergence with provable bounds
  - Optimize algorithm efficiency before running on hardware

  GitHub: https://github.com/akaiHuang/petz-recovery-unification

Happy to build a proof-of-concept on top of THRML if useful.

Best,
Sheng-Kai Huang
akai@fawstudio.com
