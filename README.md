# Universal Recovery: The Petz Map as Retrodiction Functor, Error Corrector, and Entropy Bound

**Unifying Quantum Erasers, QEC, and Thermodynamics via the Petz Map.**

[![arXiv](https://img.shields.io/badge/arXiv-TBD-b31b1b.svg)](https://arxiv.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## Overview

The Petz recovery map appears independently in quantum sufficiency theory, quantum error correction, categorical retrodiction, and quantum thermodynamics. This work proves that these appearances are not analogies but manifestations of a **single mathematical structure**.

We define a **temporal asymmetry parameter**

$$\tau = 1 - F\!\big(\rho,\, \widetilde{\mathcal{R}}_{\sigma,\mathcal{N}}(\mathcal{N}(\rho))\big)$$

that quantifies the degree to which a quantum process distinguishes past from future:
- $\tau = 0$: closed, unitary evolution (no arrow of time)
- $\tau > 0$: open system coupled to environment (arrow of time emerges)

## Core Results

### Theorem 1: Uniqueness of Retrodiction
The Petz map is the **unique** retrodiction functor satisfying Bayesian consistency, normalization, and the classical limit ([Parzygnat-Buscemi, 2023](https://doi.org/10.22331/q-2023-05-15-1013)).

### Theorem 2: Retrodiction Implies Near-Optimal QEC
The categorical uniqueness of the Petz map implies its near-optimality as a quantum error correction decoder, combining [Barnum-Knill (2002)](https://doi.org/10.1063/1.1459754) and [Junge et al. (2018)](https://doi.org/10.1007/s00023-018-0716-0).

### Theorem 3: Master Inequality Chain
For any quantum channel $\mathcal{N}$, reference state $\sigma$, and input $\rho$:

$$-\log F^2 \;\leq\; I(A;E|B) \;\leq\; \Sigma \;\leq\; \Delta D$$

linking recovery fidelity, conditional mutual information, entropy production, and relative entropy decrease into a single chain.

### Equivalence Chain
For faithful $\sigma$:

$$\tau = 0 \;\Longleftrightarrow\; \Sigma = 0 \;\Longleftrightarrow\; I(A;E|B) = 0 \;\Longleftrightarrow\; \text{Quantum Markov Chain}$$

The arrow of time ($\tau > 0$) is equivalent to environmental coupling, entropy production, and imperfect recoverability.

## Physical Predictions

1. **Post-selection as thermodynamic filtering**: For stabilizer codes under noise, rejecting high-entropy-production syndromes yields improvement $\ln R(\varepsilon, d) = \alpha(p)\, d + \beta(\varepsilon)$, where $\alpha$ is decoder-independent.

2. **Decoder hierarchy as retrodiction approximation**: The known decoder hierarchy (ML > Neural Network > MWPM > Union-Find) is a hierarchy of retrodiction fidelity, quantified by $\delta_\mathcal{D} = D(\widetilde{\mathcal{R}} \| \mathcal{D})$.

## Repository Structure

```
petz-recovery-unification/
├── paper/
│   ├── petz_recovery_unification.tex      # Main text (PRL format, 5 pages)
│   ├── petz_recovery_unification.pdf      # Compiled PDF
│   ├── petz_recovery_supplemental.tex     # Supplemental material (9 pages)
│   └── petz_recovery_supplemental.pdf     # Compiled PDF
├── simulations/
│   ├── master_inequality_chain.py         # Verify Theorem 3 numerically
│   ├── quantum_eraser_petz.py             # Quantum eraser = Petz recovery
│   └── tau_vs_entropy_production.py       # tau bound visualization
├── README.md
└── LICENSE
```

## Simulations

The `simulations/` directory contains Python scripts that numerically verify the paper's core results.

### Requirements

```bash
pip install numpy scipy qutip matplotlib
```

### Running

```bash
# Verify the master inequality chain (Theorem 3)
python simulations/master_inequality_chain.py

# Demonstrate quantum eraser as Petz recovery
python simulations/quantum_eraser_petz.py

# Visualize tau vs entropy production bound
python simulations/tau_vs_entropy_production.py
```

Each script generates plots and prints numerical verification of the corresponding theorem.

## Key References

- D. Petz, "Sufficiency of channels over von Neumann algebras," Q. J. Math. **39**, 97 (1988).
- A. J. Parzygnat and F. Buscemi, "Axioms for retrodiction," Quantum **7**, 1013 (2023).
- M. Junge, R. Renner, D. Sutter, M. M. Wilde, and A. Winter, "Universal recovery maps," Ann. Henri Poincare **19**, 2955 (2018).
- O. Fawzi and R. Renner, "Quantum conditional mutual information and approximate Markov chains," Commun. Math. Phys. **340**, 575 (2015).
- G. E. Crooks, "Entropy production fluctuation theorem," Phys. Rev. E **60**, 2721 (1999).

## Citation

```bibtex
@article{Huang2026petz,
  author  = {Sheng-Kai Huang},
  title   = {Universal Recovery: The Petz Map as Retrodiction Functor, Error Corrector, and Entropy Bound},
  year    = {2026},
  note    = {arXiv:TBD}
}
```

## Author

**Sheng-Kai Huang** — Independent Researcher
- Email: akai@fawstudio.com

## License

This project is licensed under the MIT License - see [LICENSE](LICENSE) for details.
