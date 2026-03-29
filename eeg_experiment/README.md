# EEG Consciousness-Retrodiction Experiment

## Hypothesis

**Testing τ ≤ 1 − exp(−Σ/2) across consciousness states**

From the Petz recovery framework (Paper 1, Huang 2025):
- τ (retrodiction error) measures how poorly we can reconstruct a prior brain state from the current one
- Σ (relative entropy production) quantifies irreversibility
- The Petz recovery bound: τ ≤ 1 − exp(−Σ/2) sets a fundamental limit

**Prediction**: Different consciousness states (waking, sleep stages, anesthesia)
correspond to different positions along the τ–Σ curve:

| State | Expected Σ | Expected τ | Interpretation |
|-------|-----------|-----------|----------------|
| Deep sleep (N3) | Low | Low | High recoverability — brain state is "simple", nearly reversible |
| Light sleep (N1/N2) | Medium | Medium | Intermediate |
| REM | High | High | Dream states: complex, irreversible dynamics |
| Waking (eyes open) | Highest | Highest | Maximum information processing, maximum entropy production |
| Waking (eyes closed) | High | High | Slightly less than eyes-open |
| Anesthesia | Very low | Very low | "Frozen" dynamics, minimal entropy production |

The key test: **Does the Petz bound τ ≤ 1 − exp(−Σ/2) hold empirically,
and does the gap between τ and the bound vary systematically with consciousness level?**

## Connection to Core Theory

From Sheng-Kai Huang's framework:
- Σ = D(ρ_spacetime || ρ_matter) — a single QRE with different boundary conditions
- In EEG context: Σ measures the irreversibility of brain-state transitions
- τ = 1 − F where F is the fidelity between actual and retrodicted states
- The Petz recovery map R_σ,N provides the optimal retrodiction channel

**Physical intuition**: In deep sleep (N3), the brain enters a highly ordered,
low-entropy state with slow-wave oscillations. This is analogous to a
"more closed" system — closer to the zero-entropy ideal where retrodiction
becomes perfect (τ → 0). Conversely, waking consciousness with rich sensory
input is maximally "open" — retrodiction fails most (τ → 1).

## Data Sources

### Primary: Sleep-EDF Database Expanded (PhysioNet)
- 197 whole-night polysomnographic recordings
- Hypnogram annotations: Wake, N1, N2, N3, REM, Movement
- EEG channels: Fpz-Cz, Pz-Oz (100 Hz)
- URL: https://physionet.org/content/sleep-edfx/1.0.0/
- License: Open Data Commons Attribution License v1.0

### Future extensions:
- OpenNeuro: meditation/mindfulness EEG datasets
- CAP Sleep Database (PhysioNet): cyclic alternating pattern
- Temple University Hospital EEG Corpus: pathological states
- Anesthesia EEG datasets (if available publicly)

## Project Structure

```
eeg_experiment/
├── README.md              — This file
├── data/
│   └── sleep-edf/         — Raw .edf files from PhysioNet
├── scripts/
│   ├── download_data.py   — Download Sleep-EDF subset
│   ├── load_and_preview.py — Load EDF, extract stages, plot samples
│   ├── compute_sigma_tau.py — Compute Σ and τ per epoch per stage
│   └── analyze_bound.py   — Test τ ≤ 1-exp(-Σ/2) across states
├── results/
│   ├── figures/            — Generated plots
│   └── statistics/         — CSV summaries
└── .venv/                 — Python virtual environment
```

## Quick Start

```bash
# Activate virtual environment
source .venv/bin/activate

# Download data (5 subjects, ~500 MB)
python scripts/download_data.py

# Load and preview data
python scripts/load_and_preview.py

# Compute retrodiction metrics
python scripts/compute_sigma_tau.py

# Test the Petz bound
python scripts/analyze_bound.py
```

## Methods

### Step 1: EEG Preprocessing
- Load EDF files with MNE-Python
- Extract EEG channels (Fpz-Cz, Pz-Oz)
- Band-pass filter: 0.5–45 Hz
- Epoch into 30-second windows (standard sleep scoring)

### Step 2: State Extraction
- Parse hypnogram annotations for sleep stage labels
- Group epochs by stage: W, N1, N2, N3, REM

### Step 3: Compute Σ (Relative Entropy Production)
For each 30s epoch, estimate:
- Power spectral density → probability distribution over frequencies
- Transition matrices between consecutive epochs
- Σ = D(P_forward || P_backward) via KL divergence of transition probabilities
- Alternative: spectral entropy, sample entropy, permutation entropy

### Step 4: Compute τ (Retrodiction Error)
- τ = 1 − F where F = fidelity between:
  - The actual previous epoch's state
  - The Petz-retrodicted state from the current epoch
- Operationally: fit forward model, apply Petz map, measure reconstruction error

### Step 5: Test the Bound
- Plot τ vs Σ for all epochs, colored by sleep stage
- Overlay the theoretical curve τ = 1 − exp(−Σ/2)
- Statistical tests: does the bound hold? How tight is it per stage?

## References

1. Huang, S.-K. (2025). "Retrodiction, Recovery, and the Flow of Time."
   arXiv preprint. DOI: 10.5281/zenodo.14538923
2. Kemp et al. (2000). "Analysis of a sleep-dependent neuronal feedback loop."
   IEEE-BME 47(9):1185-1194.
3. Goldberger et al. (2000). "PhysioBank, PhysioToolkit, and PhysioNet."
   Circulation 101(23):e215-e220.
4. Petz, D. (1986). "Sufficient subalgebras and the relative entropy of states."
   Commun. Math. Phys. 105, 123-131.

## License

Data: Open Data Commons Attribution License v1.0 (PhysioNet)
Code: MIT License
