# tau-chrono

**Bayesian Noise Tracker for Quantum Circuits.** Track noise evolution through your quantum circuit with per-gate tau estimates **50--100% more accurate** than independent noise models, backed by provable composition inequalities from the Petz recovery map.

## Why tau-chrono?

| Feature | pyGSTi | Mitiq | Q-CTRL | True-Q | **tau-chrono** |
|---|---|---|---|---|---|
| Bayesian noise tracking | - | - | - | - | **Yes** |
| Provable composition bounds | - | - | - | - | **Yes** |
| Petz recovery metric | - | - | - | - | **Yes** |
| Process tomography | Yes | - | - | Yes | **Yes** |
| Error mitigation | - | Yes | Yes | - | **Yes (tau-informed)** |

## Installation

```bash
pip install tau-chrono                # core (numpy only)
pip install tau-chrono[qiskit]        # + IBM Quantum support
pip install tau-chrono[viz]           # + matplotlib plots
pip install tau-chrono[all]           # everything
```

## Quick start (no hardware needed)

```python
import numpy as np
from tau_chrono import amplitude_damping, depolarizing, bayesian_compose

# Define a 3-gate noise sequence
channels = [amplitude_damping(0.05), depolarizing(0.08), amplitude_damping(0.12)]
names    = ["AD(0.05)", "Dep(0.08)", "AD(0.12)"]

# Input and reference states
rho   = 0.5 * np.array([[1, 1], [1, 1]], dtype=complex)   # |+>
sigma = np.diag([0.8, 0.2]).astype(complex)

# Run Bayesian analysis
result = bayesian_compose(channels, sigma, rho, channel_names=names)

print(f"Improvement over naive: {result.improvement_percent:.1f}%")
print(f"Composition inequality: {'PASS' if result.composition_holds else 'FAIL'}")
for gr in result.gate_results:
    print(f"  {gr.channel_name}: tau_eff={gr.tau_eff:.4f} ({gr.classification})")
```

## CLI

```bash
tau-chrono analyze "ad:0.05,dep:0.08,ad:0.12"   # quick analysis
tau-chrono benchmark --suite 1q                   # run benchmarks
```

## Hardware validation (Tuna-9)

Validated on QuTech Quantum Inspire Tuna-9 superconducting backend:

| Depth | Naive tau | Bayesian tau | Improvement |
|:-----:|-----------|-------------|:-----------:|
| 8     | 0.323     | 0.290       | 10.2%       |
| 15    | 0.509     | 0.401       | 21.3%       |
| 30    | 0.759     | 0.477       | **37.1%**   |

## Free vs Pro

| | Free | Pro |
|---|---|---|
| Core tau analysis | Yes | Yes |
| Up to 20 gates | Yes | Yes |
| Unlimited gates | - | Yes |
| 2-qubit tomography | - | Yes |
| Non-Markovian correction | - | Yes |
| Error mitigation (ZNE/PEC) | - | Yes |
| Noise drift monitoring | - | Yes |

## License

MIT
