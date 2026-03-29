# Experimental Protocol: Measuring the Petz Recovery Saturation Gap in EEG Data

**Author**: Sheng-Kai Huang
**Date**: 2026-03-24
**Status**: Protocol design (pre-implementation)
**Novelty**: This measurement has NEVER been performed in neuroscience.

---

## 0. Executive Summary

We propose the first measurement of the **Petz recovery saturation gap**

$$\Delta = F_{\text{actual}} - e^{-\Sigma/2}$$

in human EEG data across consciousness states. Here $F_{\text{actual}}$ is the fidelity achieved by the explicit Petz recovery map (Bayesian retrodiction) applied to the neural channel, $e^{-\Sigma/2}$ is the JRSWW lower bound (Junge et al. 2018), and $\Sigma$ is the entropy production of the neural channel (Schnakenberg formula).

**Key theoretical insight**: The JRSWW bound $F \geq e^{-\Sigma/2}$ is known to be non-saturated for $\Sigma > 0$ (Gao-Junge-LaRacuente 2023). The gap $\Delta$ quantifies how much *better* the brain performs at retrodiction than the fundamental limit. We hypothesize that $\Delta$ varies systematically with consciousness level and may connect to integrated information ($\Phi$).

**Data**: 23-subject Sleep-EDF dataset already downloaded (8,620 epochs across Wake, N1, N2, N3, REM).

---

## 1. Theoretical Foundation

### 1.1 The Petz Bound for Classical Markov Chains

For a discrete-time Markov chain with transition matrix $T$ and stationary distribution $\pi$, the JRSWW bound specializes to (see `layer4_consciousness_conditions.md`, Section 2.4):

$$F(p,\; R_{\pi,T} \circ T(p)) \;\geq\; e^{-\Sigma_{\text{path}}/2}$$

where:
- $p$ = initial distribution over neural states
- $T(p)$ = evolved distribution after one time step
- $R_{\pi,T}$ = Petz recovery map (Bayesian retrodiction)
- $\Sigma_{\text{path}} = D_{\text{KL}}(P_{\text{fwd}} \| P_{\text{bwd}})$ = path-level entropy production

At stationarity (NESS), $\Sigma_{\text{path}} = \Sigma_{\text{Schnakenberg}}$.

### 1.2 The Saturation Gap

Define:

$$\Delta = F_{\text{actual}} - e^{-\Sigma/2}$$

Properties:
- $\Delta \geq 0$ always (by the JRSWW bound)
- $\Delta = 0$ only when $\Sigma = 0$ (perfect recovery = no entropy production)
- $\Delta > 0$ for all $\Sigma > 0$ (proven non-saturation for positive EP)

The **normalized gap** is more informative:

$$\bar{\Delta} = \frac{F_{\text{actual}} - e^{-\Sigma/2}}{1 - e^{-\Sigma/2}}$$

This ranges from 0 (at the bound) to 1 (perfect recovery despite nonzero EP). Physically, $\bar{\Delta}$ measures the fraction of "allowed recovery room" that the brain actually uses.

### 1.3 Connection to Previous Results

From our existing 23-subject analysis (`multi_subject_stats.json`):
- Wake: $\Sigma \approx 0.595$, $F \approx 0.965$, bound $= e^{-0.595/2} \approx 0.743$
- N3: $\Sigma \approx 0.114$, $F \approx 0.987$, bound $= e^{-0.114/2} \approx 0.945$

These imply $\Delta_{\text{Wake}} \approx 0.222$ and $\Delta_{\text{N3}} \approx 0.042$. However, the previous analysis used a **simplified Petz map** (PSD-based Bayesian retrodiction with uniform prior). This protocol designs the **exact Petz recovery map** for the empirically estimated transition matrix.

### 1.4 Why This Matters

1. **First empirical measurement** of how close neural dynamics come to a fundamental information-theoretic bound
2. **Consciousness signature**: If $\Delta$ correlates with consciousness level, it provides a new mathematically grounded consciousness metric
3. **Bridge to IIT**: The gap $\Delta$ is structurally analogous to IIT's $\Phi$ — both measure "more than the sum of parts" in information-theoretic terms
4. **Falsifiable prediction**: We predict $\Delta_{\text{Wake}} > \Delta_{\text{N3}} > 0$ (conscious states recover better than the bound predicts)

---

## 2. Data and Prerequisites

### 2.1 Dataset

**Sleep-EDF Database Expanded** (PhysioNet, already downloaded):
- 23 subjects, ~49 EDF files in `/eeg_experiment/data/sleep-edf/`
- Channels: Fpz-Cz, Pz-Oz (100 Hz sampling)
- Annotations: Wake (W), N1, N2, N3, REM
- Total epochs: 8,620 (30-second epochs)
- Distribution: W=1840, N1=1451, N2=1840, N3=1661, REM=1828

### 2.2 Python Dependencies

```
mne>=1.6.0           # EEG loading and preprocessing
numpy>=1.24          # Core numerics
scipy>=1.11          # Welch PSD, linear algebra, statistics
scikit-learn>=1.3    # PCA for state space reduction
matplotlib>=3.8      # Visualization
seaborn>=0.13        # Statistical plotting
pandas>=2.1          # Data management
tqdm>=4.66           # Progress bars
statsmodels>=0.14    # Bootstrapping and statistical tests
```

Install:
```bash
cd /Users/akaihuangm1/Desktop/github/petz-recovery-unification/eeg_experiment
source .venv/bin/activate
pip install mne numpy scipy scikit-learn matplotlib seaborn pandas tqdm statsmodels
```

### 2.3 Computational Requirements

- RAM: ~8 GB (transition matrices are sparse but can be large)
- Time estimate: ~2-4 hours for full 23-subject analysis (dominated by Method 3)
- Storage: ~500 MB for results

---

## 3. Preprocessing Pipeline

### 3.1 EEG Loading and Filtering

Reuse the existing `improved_sigma.py` infrastructure:

```
INPUT:  Raw EDF files (SC40xxE0-PSG.edf + SC40xxEx-Hypnogram.edf)
OUTPUT: epochs[n_epochs, n_channels, n_samples], stages[n_epochs], sfreq

Steps:
  1. mne.io.read_raw_edf(psg_path, preload=True)
  2. Pick EEG channels only (Fpz-Cz, Pz-Oz)
  3. Bandpass filter: 0.5-45 Hz (4th order Butterworth, zero-phase)
  4. Segment into 30-second non-overlapping epochs
  5. Parse hypnogram annotations → stage labels per epoch
  6. Reject epochs with amplitude > 500 μV (artifact rejection)
```

### 3.2 State Space Construction

The critical design choice. We need a discrete state space to build transition matrices.

**Method A: Frequency-band binarization (primary, from improved_sigma.py)**

```
INPUT:  epoch[n_channels, 3000]  (30s at 100 Hz)
OUTPUT: state_sequence[3000], n_states

Steps:
  1. Expand 2 channels to 8 virtual channels via bandpass filtering:
     - delta (0.5-4 Hz), theta (4-8 Hz), alpha (8-13 Hz), beta (13-30 Hz)
     - per physical channel → 2 × 4 = 8 virtual channels
  2. Binarize each virtual channel using running median (window=500 samples)
  3. PCA to max_bits=6 effective dimensions → 2^6 = 64 states
  4. Encode: state_t = Σ_k binary_k(t) × 2^k
```

**Method B: Voltage amplitude binning (secondary, simpler)**

```
INPUT:  epoch[n_channels, 3000]
OUTPUT: state_sequence[3000], n_states

Steps:
  1. Average across channels → signal[3000]
  2. Compute amplitude histogram with n_bins (default: 32)
  3. Assign each sample to a bin → state_sequence
  4. n_states = n_bins
```

**Method C: Time-frequency binning (tertiary, highest resolution)**

```
INPUT:  epoch[n_channels, 3000]
OUTPUT: state_sequence[n_windows], n_states

Steps:
  1. Short-time Fourier transform (window=256, overlap=128)
  2. Extract power in 5 canonical bands (delta, theta, alpha, beta, gamma)
  3. Discretize each band power into 3 levels (low/medium/high)
     using tercile thresholds computed per-subject
  4. Combined state: 3^5 = 243 possible states
  5. Merge rare states (< 5 occurrences) into nearest neighbor
```

**Recommendation**: Use Method A as primary (consistent with existing code). Use Method B for robustness checks. Method C for a higher-dimensional verification.

### 3.3 Sub-epoch Windowing

For measuring the gap WITHIN a single 30-second epoch (intra-epoch dynamics):

```
INPUT:  epoch[n_channels, 3000]
PARAM:  window_size = 1s (100 samples)
        step_size = 0.5s (50 samples, 50% overlap)
OUTPUT: List of sub-epoch windows

This gives ~59 windows per epoch → 58 transitions for the transition matrix.
```

For measuring the gap ACROSS epochs (inter-epoch dynamics):

```
INPUT:  consecutive same-stage epochs
OUTPUT: One transition per epoch pair

This gives fewer data points but captures slow dynamics.
```

**Recommendation**: Use intra-epoch analysis as primary (more data, better statistics). Use inter-epoch as validation.

---

## 4. Core Algorithms

### 4.1 Algorithm 1: Compute Transition Matrix $T$

```python
def compute_transition_matrix(states, n_states, lag=1):
    """
    Estimate T[i,j] = P(state_{t+lag} = j | state_t = i)

    INPUT:
      states: array of integer state indices, length L
      n_states: total number of possible states
      lag: time lag (default 1)

    OUTPUT:
      T: row-stochastic transition matrix [n_states × n_states]
      counts: raw count matrix (for statistical tests)
    """
    counts = np.zeros((n_states, n_states), dtype=np.float64)
    for t in range(len(states) - lag):
        i, j = int(states[t]), int(states[t + lag])
        counts[i, j] += 1

    # Row-normalize with Laplace smoothing (avoid zero rows)
    alpha = 1e-6  # Laplace smoothing parameter
    T = (counts + alpha) / (counts.sum(axis=1, keepdims=True) + alpha * n_states)

    return T, counts
```

### 4.2 Algorithm 2: Compute Stationary Distribution $\pi$

```python
def stationary_distribution(T):
    """
    Compute π such that π^T T = π^T (left eigenvector of eigenvalue 1).

    INPUT:  T: row-stochastic matrix [n × n]
    OUTPUT: pi: stationary distribution [n]

    Method: eigenvalue decomposition of T^T, extract eigenvector for λ=1.
    Fallback: power iteration if eigendecomposition fails.
    """
    n = T.shape[0]
    eigenvalues, eigenvectors = scipy.linalg.eig(T.T, left=False, right=True)
    idx = np.argmin(np.abs(eigenvalues - 1.0))
    pi = np.real(eigenvectors[:, idx])
    pi = np.abs(pi)

    if pi.sum() > 0:
        pi /= pi.sum()
    else:
        # Fallback: power iteration
        pi = np.ones(n) / n
        for _ in range(1000):
            pi_new = pi @ T
            if np.allclose(pi_new, pi, atol=1e-12):
                break
            pi = pi_new
        pi = np.abs(pi)
        pi /= pi.sum()

    # Ensure strict positivity (needed for Petz map)
    eps = 1e-15
    pi = np.maximum(pi, eps)
    pi /= pi.sum()

    return pi
```

### 4.3 Algorithm 3: Compute Schnakenberg $\Sigma$

```python
def schnakenberg_entropy_production(T, pi, epsilon=1e-20):
    """
    Schnakenberg entropy production rate:

    Σ = Σ_{i,j} π_i T_{ij} ln(π_i T_{ij} / (π_j T_{ji}))

    This equals D_KL(P_fwd || P_bwd) where:
      P_fwd(i,j) = π_i T_{ij}    (forward joint distribution)
      P_bwd(i,j) = π_j T_{ji}    (backward joint distribution)

    INPUT:
      T: transition matrix [n × n]
      pi: stationary distribution [n]

    OUTPUT:
      sigma: Schnakenberg EP rate (nats per step), guaranteed >= 0
    """
    n = T.shape[0]

    # Joint distributions
    P_fwd = T * pi[:, np.newaxis]   # P_fwd[i,j] = pi[i] * T[i,j]
    P_bwd = T.T * pi[:, np.newaxis] # P_bwd[i,j] = pi[j] * T[j,i]
    # Note: P_bwd[i,j] = pi[j] T[j,i] after transposition
    # Actually need: P_bwd(i→j) = π_j T_{ji}
    # So P_bwd[i,j] = pi[j] * T[j,i]
    P_bwd_correct = np.zeros_like(T)
    for i in range(n):
        for j in range(n):
            P_bwd_correct[i, j] = pi[j] * T[j, i]

    # Vectorized version:
    P_fwd = T * pi[:, np.newaxis]           # P_fwd[i,j] = pi_i T_ij
    P_bwd = (T * pi[:, np.newaxis]).T       # P_bwd[j,i] = pi_i T_ij → P_bwd[i,j] = pi_j T_ji

    # KL divergence: Σ = Σ_{ij} P_fwd[i,j] * ln(P_fwd[i,j] / P_bwd[i,j])
    mask = (P_fwd > epsilon) & (P_bwd > epsilon)

    P_f = np.where(mask, np.maximum(P_fwd, epsilon), epsilon)
    P_b = np.where(mask, np.maximum(P_bwd, epsilon), epsilon)

    sigma = np.sum(np.where(mask, P_fwd * np.log(P_f / P_b), 0.0))

    return max(0.0, sigma)
```

### 4.4 Algorithm 4: Construct Explicit Petz Recovery Map $R_{\pi,T}$

**This is the key new algorithm not in the existing codebase.**

```python
def petz_recovery_map(T, pi):
    """
    Construct the Petz recovery map R for transition matrix T
    with reference state π.

    For classical Markov chains, the Petz map is Bayesian retrodiction:

      R[j → i] = π_i T_{ij} / Σ_k π_k T_{kj}
                = π_i T_{ij} / π_j^{(out)}

    where π_j^{(out)} = Σ_k π_k T_{kj} = (T^T π)_j

    This is EXACTLY Bayes' theorem:
      P(X_t = i | X_{t+1} = j) = P(X_{t+1}=j | X_t=i) P(X_t=i) / P(X_{t+1}=j)

    INPUT:
      T: forward transition matrix [n × n], T[i,j] = P(j|i)
      pi: stationary distribution [n]

    OUTPUT:
      R: recovery (retrodiction) matrix [n × n], R[j,i] = P(i|j)
         i.e., R[j,i] is the probability of recovering state i
         given that we observed state j
    """
    n = T.shape[0]

    # Output distribution under stationarity: π_out = T^T π
    pi_out = T.T @ pi  # pi_out[j] = Σ_i π_i T_{ij}

    # Petz recovery matrix
    R = np.zeros((n, n), dtype=np.float64)
    for j in range(n):
        if pi_out[j] > 1e-15:
            for i in range(n):
                R[j, i] = pi[i] * T[i, j] / pi_out[j]
        else:
            # If state j is never reached, use uniform recovery
            R[j, :] = 1.0 / n

    # Vectorized version:
    # R[j, i] = pi[i] * T[i, j] / pi_out[j]
    # R = (T * pi[:, None]).T / pi_out[:, None]
    # But be careful with zero pi_out entries

    return R


def apply_petz_recovery(R, q):
    """
    Apply the Petz recovery map to distribution q.

    INPUT:
      R: recovery matrix [n × n], R[j,i] = P(retrodicted state = i | observed = j)
      q: distribution over observed states [n]

    OUTPUT:
      p_recovered: retrodicted distribution [n]

    p_recovered[i] = Σ_j R[j,i] * q[j]
                   = Σ_j [π_i T_{ij} / π_j^{out}] * q[j]
    """
    p_recovered = R.T @ q  # Sum over j: Σ_j R[j,i] q[j]

    # Normalize (should already sum to 1, but numerical safety)
    total = p_recovered.sum()
    if total > 0:
        p_recovered /= total
    else:
        p_recovered = np.ones(len(q)) / len(q)

    return p_recovered
```

### 4.5 Algorithm 5: Compute Actual Fidelity $F_{\text{actual}}$

```python
def classical_fidelity(p, q):
    """
    Classical fidelity (squared Bhattacharyya coefficient):

    F(p, q) = (Σ_i √(p_i q_i))^2

    For classical distributions, this is the Uhlmann fidelity
    restricted to diagonal density matrices.

    INPUT:  p, q: probability distributions [n]
    OUTPUT: F: fidelity in [0, 1]
    """
    bc = np.sum(np.sqrt(np.maximum(p, 0) * np.maximum(q, 0)))
    return bc ** 2


def compute_petz_fidelity(T, pi, p_initial):
    """
    Compute the actual Petz recovery fidelity for a specific initial state.

    1. Apply forward channel: p_initial → T(p_initial)
    2. Construct Petz map: R = R_{π,T}
    3. Apply recovery: T(p_initial) → R(T(p_initial))
    4. Compute fidelity: F = F(p_initial, R ∘ T(p_initial))

    INPUT:
      T: transition matrix [n × n]
      pi: stationary distribution [n]
      p_initial: initial distribution [n]

    OUTPUT:
      F_actual: fidelity of Petz recovery
      p_recovered: the retrodicted distribution
    """
    # Step 1: Forward evolution
    p_evolved = T.T @ p_initial  # p_evolved[j] = Σ_i T[i,j] p_initial[i]
    # (Note: if T is row-stochastic, T[i,j] = P(j|i), then
    #  the evolved distribution is p_evolved[j] = Σ_i p[i] T[i,j] = (p^T T)[j])
    p_evolved = p_initial @ T
    p_evolved = np.maximum(p_evolved, 1e-15)
    p_evolved /= p_evolved.sum()

    # Step 2: Construct Petz recovery map
    R = petz_recovery_map(T, pi)

    # Step 3: Apply recovery
    p_recovered = apply_petz_recovery(R, p_evolved)

    # Step 4: Compute fidelity
    F_actual = classical_fidelity(p_initial, p_recovered)

    return F_actual, p_recovered
```

### 4.6 Algorithm 6: Compute the Saturation Gap $\Delta$

```python
def compute_saturation_gap(T, pi, p_initial):
    """
    Compute the Petz recovery saturation gap.

    Δ = F_actual - e^{-Σ/2}
    Δ_normalized = (F_actual - e^{-Σ/2}) / (1 - e^{-Σ/2})

    INPUT:
      T: transition matrix
      pi: stationary distribution
      p_initial: initial state distribution

    OUTPUT:
      delta: raw gap (≥ 0)
      delta_norm: normalized gap (0 = at bound, 1 = perfect recovery)
      F_actual: actual fidelity
      sigma: Schnakenberg EP
      bound: e^{-Σ/2}
    """
    # Compute Σ
    sigma = schnakenberg_entropy_production(T, pi)

    # Compute JRSWW bound
    bound = np.exp(-sigma / 2.0)

    # Compute F_actual via explicit Petz map
    F_actual, p_recovered = compute_petz_fidelity(T, pi, p_initial)

    # Raw gap
    delta = F_actual - bound

    # Normalized gap
    if bound < 1.0 - 1e-10:
        delta_norm = delta / (1.0 - bound)
    else:
        delta_norm = 0.0  # Σ ≈ 0, bound ≈ 1

    return {
        'delta': max(0.0, delta),
        'delta_norm': np.clip(delta_norm, 0.0, 1.0),
        'F_actual': F_actual,
        'sigma': sigma,
        'bound': bound,
    }
```

### 4.7 Algorithm 7: Epoch-Level Gap Estimation

```python
def compute_epoch_gap(epoch_data, sfreq, method='binarize',
                      max_bits=6, lag=1, n_bins=32):
    """
    Full pipeline for one epoch.

    INPUT:
      epoch_data: [n_channels, n_samples] EEG data for one 30s epoch
      sfreq: sampling frequency (Hz)
      method: 'binarize' (Method A) or 'voltage_bin' (Method B)
      max_bits: bits for state space reduction (Method A)
      n_bins: number of bins (Method B)
      lag: transition lag in samples

    OUTPUT:
      result dict with delta, delta_norm, F_actual, sigma, bound, diagnostics
    """
    if method == 'binarize':
        # Method A: frequency-band binarization
        expanded = expand_to_frequency_bands(epoch_data, sfreq)
        binary = binarize_channels(expanded, method='running_median')
        states, n_states = reduce_state_space(binary, max_bits=max_bits)
    elif method == 'voltage_bin':
        # Method B: voltage amplitude binning
        signal = epoch_data.mean(axis=0)  # average channels
        bin_edges = np.linspace(
            np.percentile(signal, 1),
            np.percentile(signal, 99),
            n_bins + 1
        )
        states = np.digitize(signal, bin_edges[1:-1])
        n_states = n_bins
    else:
        raise ValueError(f"Unknown method: {method}")

    # Build transition matrix
    T, counts = compute_transition_matrix(states, n_states, lag=lag)

    # Stationary distribution
    pi = stationary_distribution(T)

    # Empirical initial distribution (first quarter of epoch)
    quarter = len(states) // 4
    state_counts = np.bincount(states[:quarter], minlength=n_states).astype(float)
    p_initial = state_counts / state_counts.sum()
    p_initial = np.maximum(p_initial, 1e-15)
    p_initial /= p_initial.sum()

    # Compute gap
    result = compute_saturation_gap(T, pi, p_initial)

    # Diagnostics
    result['n_states_occupied'] = np.sum(pi > 1e-10)
    result['n_transitions'] = int(counts.sum())
    result['markov_order'] = _test_markov_order(states, n_states)

    return result


def _test_markov_order(states, n_states, max_order=3):
    """
    Test whether Markov(1) is adequate by comparing with Markov(2).
    Returns the BIC difference: positive = Markov(1) preferred.
    """
    n = len(states)

    # Markov(1): log-likelihood
    T1, counts1 = compute_transition_matrix(states, n_states, lag=1)
    ll1 = 0.0
    for t in range(n - 1):
        i, j = int(states[t]), int(states[t + 1])
        if T1[i, j] > 1e-15:
            ll1 += np.log(T1[i, j])
    k1 = n_states * (n_states - 1)  # free parameters
    bic1 = -2 * ll1 + k1 * np.log(n)

    # Markov(2): log-likelihood (simplified — check 2-step memory)
    # Build T2[i,j,k] = P(state_{t+1}=k | state_{t-1}=i, state_t=j)
    counts2 = np.zeros((n_states, n_states, n_states))
    for t in range(1, n - 1):
        i, j, k = int(states[t-1]), int(states[t]), int(states[t+1])
        counts2[i, j, k] += 1

    ll2 = 0.0
    for i in range(n_states):
        for j in range(n_states):
            row_sum = counts2[i, j, :].sum()
            if row_sum > 0:
                T2_ijk = counts2[i, j, :] / row_sum
                for k in range(n_states):
                    if counts2[i, j, k] > 0 and T2_ijk[k] > 1e-15:
                        ll2 += counts2[i, j, k] * np.log(T2_ijk[k])
    k2 = n_states ** 2 * (n_states - 1)
    bic2 = -2 * ll2 + k2 * np.log(n)

    return bic1 - bic2  # positive = Markov(1) preferred
```

---

## 5. Full Analysis Pipeline

### 5.1 Main Loop: All Subjects, All Stages

```python
def run_full_analysis(data_dir, results_dir, method='binarize',
                      max_bits=6, lag=1):
    """
    Main analysis loop over all subjects and consciousness stages.

    PSEUDOCODE:

    FOR each subject in data_dir:
        Load EEG data + hypnogram
        FOR each stage in [Wake, N1, N2, N3, REM]:
            Collect all epochs with this stage label
            Group into runs of consecutive same-stage epochs
            FOR each run with >= 3 consecutive epochs:
                FOR each epoch in the run:
                    result = compute_epoch_gap(epoch, sfreq, method, ...)
                    Store: (subject, stage, epoch_idx, result)

    Aggregate results by stage
    Compute statistics (mean, std, median, IQR, bootstrap CI)
    Generate figures
    Run statistical tests
    """
    all_results = []

    pairs = find_subject_pairs(data_dir)

    for psg_path, hyp_path, subject_id in tqdm(pairs, desc="Subjects"):
        try:
            epochs, ch_names, sfreq, stages = load_eeg_data(
                str(psg_path), str(hyp_path)
            )
        except Exception as e:
            logger.warning(f"Skipping {subject_id}: {e}")
            continue

        for stage in ['W', 'N1', 'N2', 'N3', 'R']:
            stage_mask = (stages == stage)
            stage_indices = np.where(stage_mask)[0]

            if len(stage_indices) < 3:
                continue

            # Find runs of consecutive same-stage epochs
            runs = _find_consecutive_runs(stage_indices, min_length=3)

            for run in runs:
                for epoch_idx in run:
                    epoch_data = epochs[epoch_idx]

                    # Artifact rejection
                    if np.max(np.abs(epoch_data)) > 500e-6:
                        continue

                    result = compute_epoch_gap(
                        epoch_data, sfreq,
                        method=method,
                        max_bits=max_bits,
                        lag=lag
                    )
                    result['subject'] = subject_id
                    result['stage'] = stage
                    result['epoch_idx'] = int(epoch_idx)
                    all_results.append(result)

    return pd.DataFrame(all_results)
```

### 5.2 Multi-Scale Gap Analysis

```python
def compute_multiscale_gap(epoch_data, sfreq, max_bits=6):
    """
    Compute the saturation gap at multiple time scales.

    Lags (for 100 Hz):
      lag=1    → 10 ms   (gamma band dynamics)
      lag=10   → 100 ms  (beta band dynamics)
      lag=50   → 500 ms  (alpha band dynamics)
      lag=100  → 1 s     (theta band dynamics)
      lag=500  → 5 s     (delta band dynamics)

    This reveals WHICH timescale contributes most to the gap.
    """
    lags = [1, 10, 50, 100, 500]
    results = {}

    for lag in lags:
        timescale_ms = lag * 1000.0 / sfreq
        result = compute_epoch_gap(
            epoch_data, sfreq,
            method='binarize',
            max_bits=max_bits,
            lag=lag
        )
        results[f'lag_{lag}'] = {
            'timescale_ms': timescale_ms,
            **result
        }

    return results
```

### 5.3 Cross-Channel Gap Analysis

```python
def compute_cross_channel_gap(epoch_data, sfreq, max_bits=6):
    """
    Compute the gap separately for different channel configurations:

    1. Fpz-Cz only (frontal)
    2. Pz-Oz only (occipital)
    3. Both channels jointly

    This tests whether the gap is region-specific.
    """
    results = {}

    n_channels = epoch_data.shape[0]
    configs = {
        'frontal': epoch_data[:1, :],
        'occipital': epoch_data[1:2, :] if n_channels > 1 else epoch_data[:1, :],
        'joint': epoch_data,
    }

    for name, data in configs.items():
        result = compute_epoch_gap(
            data, sfreq, method='binarize', max_bits=max_bits
        )
        results[name] = result

    return results
```

---

## 6. Statistical Analysis Plan

### 6.1 Primary Analyses

**Analysis 1: Gap by consciousness state**
```
H0: Δ does not differ across consciousness states
H1: Δ_Wake > Δ_N1 > Δ_N2 > Δ_N3 (ordered alternative)

Test: Jonckheere-Terpstra test for ordered alternatives
Post-hoc: Mann-Whitney U with Bonferroni correction for all pairs
Effect size: rank-biserial correlation r

Report: median Δ, IQR, and bootstrap 95% CI per stage
```

**Analysis 2: Normalized gap by consciousness state**
```
Same as Analysis 1, but using Δ̄ (normalized gap).
This controls for the effect of different Σ magnitudes across stages.
```

**Analysis 3: Bound compliance**
```
For each epoch, test: F_actual >= e^{-Σ/2}?
Report: compliance rate (%) per stage
Expected: ~100% (violations indicate estimation errors, not physics)
Non-compliance analysis: which epochs violate, and why?
```

### 6.2 Secondary Analyses

**Analysis 4: Within-subject consistency**
```
For each subject, compute median Δ per stage.
ICC(2,1) intraclass correlation across subjects.
This tests whether Δ is a reliable individual-level measure.
```

**Analysis 5: Timescale dependence**
```
From multiscale gap analysis:
Plot Δ(timescale) for each consciousness state.
Does the gap peak at a specific timescale?
Hypothesis: alpha-band timescale (~100 ms) shows largest Δ during Wake.
```

**Analysis 6: Correlation with traditional EEG markers**
```
Compute per epoch:
  - Spectral entropy (Shannon entropy of PSD)
  - Alpha power (8-13 Hz relative power)
  - Theta/alpha ratio
  - Permutation entropy

Correlate each with Δ.
Report: Spearman rank correlations with bootstrap CI.
```

### 6.3 Robustness Checks

**Check 1: State space size sensitivity**
```
Repeat primary analysis with max_bits = {4, 5, 6, 7, 8}.
Report: how Δ changes with state space dimensionality.
If Δ is robust across max_bits, the result is trustworthy.
```

**Check 2: Method sensitivity**
```
Compare Δ from Method A (binarize) vs Method B (voltage bins) vs Method C (TF bins).
Kruskal-Wallis test across methods per stage.
```

**Check 3: Lag sensitivity**
```
Compute Δ at lag = {1, 2, 5, 10, 20}.
Report: how Δ changes with lag.
```

**Check 4: Markov order adequacy**
```
Report BIC difference (Markov(1) vs Markov(2)) per epoch.
Flag epochs where Markov(1) is inadequate.
Rerun analysis excluding those epochs.
```

**Check 5: Bootstrap confidence intervals**
```
For each stage:
  1. Resample epochs with replacement (1000 iterations)
  2. Compute median Δ for each bootstrap sample
  3. Report 95% CI from percentiles
```

---

## 7. Visualization Plan

### 7.1 Figure 1: The Saturation Gap Landscape

**Panel A**: Scatter plot of F_actual vs Σ for all epochs.
- x-axis: Σ (Schnakenberg entropy production)
- y-axis: F (actual Petz fidelity)
- Color: consciousness state (Wake=red, N1=orange, N2=green, N3=blue, REM=purple)
- Overlay: $F = e^{-\Sigma/2}$ curve (JRSWW bound, black dashed)
- Overlay: F = 1 line (perfect recovery)
- The vertical distance between each point and the dashed curve IS the gap Δ.

**Panel B**: Same as Panel A but with error bars (mean ± SEM per stage).

### 7.2 Figure 2: Gap Distribution by State

**Panel A**: Violin plots of Δ by consciousness state.
- x-axis: consciousness state (ordered W, N1, N2, N3, REM)
- y-axis: Δ (raw gap)
- Show median, IQR, individual points

**Panel B**: Same for normalized gap Δ̄.

### 7.3 Figure 3: Multiscale Gap Profile

**Panel A**: Δ vs timescale for each consciousness state.
- x-axis: timescale (log scale, ms)
- y-axis: median Δ
- Lines colored by consciousness state
- Error bands: IQR

**Panel B**: Δ̄ vs timescale (normalized).

### 7.4 Figure 4: Individual Subject Trajectories

**Panel**: For 4 representative subjects, plot Δ as a function of time through the night.
- x-axis: time (hours)
- y-axis: Δ
- Background shading: hypnogram (sleep stage)
- Shows how Δ tracks consciousness transitions in real time.

### 7.5 Figure 5: Correlation Matrix

Heatmap of Spearman correlations between:
- Δ, Δ̄, Σ, F_actual, spectral entropy, alpha power, theta/alpha ratio, permutation entropy

### 7.6 Figure 6: Bound Compliance

**Panel A**: Histogram of F_actual / bound per stage.
- Values > 1: bound satisfied
- Values < 1: bound violated (should be rare)

**Panel B**: Compliance rate vs max_bits (robustness check).

---

## 8. Expected Results and Interpretation

### 8.1 Primary Predictions

Based on existing data (Section 1.3) and theoretical reasoning:

| Stage | Expected Σ | Expected F_actual | Expected bound | Expected Δ | Expected Δ̄ |
|-------|-----------|-------------------|---------------|-----------|------------|
| Wake  | ~0.5-0.7  | ~0.95-0.97        | ~0.72-0.78    | ~0.18-0.25 | ~0.70-0.85 |
| N1    | ~0.15-0.20| ~0.94-0.96        | ~0.90-0.93    | ~0.02-0.06 | ~0.30-0.50 |
| N2    | ~0.10-0.15| ~0.95-0.97        | ~0.93-0.95    | ~0.01-0.04 | ~0.20-0.40 |
| N3    | ~0.08-0.12| ~0.98-0.99        | ~0.94-0.96    | ~0.02-0.05 | ~0.40-0.60 |
| REM   | ~0.12-0.18| ~0.95-0.97        | ~0.91-0.94    | ~0.02-0.06 | ~0.30-0.50 |

**Key prediction**: $\Delta_{\text{Wake}} \gg \Delta_{\text{N2}}$ (waking brain recovers far above the bound).

**Surprise prediction**: $\Delta_{\text{N3}}$ may be unexpectedly large (slow-wave sleep has highly structured dynamics that enable efficient retrodiction despite low EP).

### 8.2 Interpretation Framework

**Scenario A: Δ monotonically increases with consciousness level**
- Wake > REM > N1 > N2 > N3
- Interpretation: Conscious brains are "better recoverers" — they generate more entropy but compensate with richer internal models.
- Connection: $\Delta$ measures the brain's "excess retrodiction capacity" = a form of integrated information.

**Scenario B: Δ is large for BOTH Wake and N3, small for N1/N2**
- U-shaped: Wake ≈ N3 > REM > N1 > N2
- Interpretation: Both highly conscious (Wake) and highly organized (N3) states deviate from the bound, but for different reasons:
  - Wake: complex dynamics + strong causal structure → good recovery
  - N3: low EP (Σ small, bound near 1) + very structured slow waves → also good recovery
  - N1/N2: moderate EP + disorganized dynamics → closer to the bound

**Scenario C: Normalized gap Δ̄ reveals the true ordering**
- Since Δ depends on Σ (larger Σ → more room for gap), the normalized gap Δ̄ = Δ/(1 - bound) may show a different ordering than raw Δ.
- If Δ̄_Wake ≫ Δ̄_N3, then waking consciousness is genuinely better at recovery (not just because Σ is larger).

### 8.3 Connection to IIT

If Δ correlates with consciousness level:
- $\Delta \sim \Phi$ (integrated information) up to monotone transformation
- The advantage: $\Delta$ is computable from EEG data in polynomial time, while $\Phi$ requires exponential time
- This would provide a **tractable proxy for $\Phi$** grounded in rigorous information theory

### 8.4 What Would Invalidate This Approach

1. **Bound violations > 5%**: Would indicate systematic errors in Σ or F estimation
2. **No stage dependence of Δ**: Would mean the gap is noise, not signal
3. **Δ depends strongly on state-space parameters**: Would mean the result is an artifact of discretization
4. **Markov(1) universally inadequate**: Would require higher-order Markov models (computationally expensive but doable)

---

## 9. Implementation Plan

### 9.1 Phase 1: Core Implementation (Day 1)

```
File: eeg_experiment/scripts/saturation_gap.py

Implement:
  - compute_transition_matrix() [exists, adapt]
  - stationary_distribution() [exists, adapt]
  - schnakenberg_entropy_production() [exists, adapt]
  - petz_recovery_map() [NEW]
  - apply_petz_recovery() [NEW]
  - compute_petz_fidelity() [NEW]
  - compute_saturation_gap() [NEW]
  - compute_epoch_gap() [NEW]

Test on 1 subject first.
```

### 9.2 Phase 2: Full Analysis (Day 2)

```
File: eeg_experiment/scripts/saturation_gap_analysis.py

Implement:
  - run_full_analysis() — all 23 subjects
  - Statistical tests (Kruskal-Wallis, Jonckheere-Terpstra, Mann-Whitney)
  - Bootstrap confidence intervals
  - All 6 figures
  - Export to CSV + JSON
```

### 9.3 Phase 3: Robustness (Day 3)

```
File: eeg_experiment/scripts/saturation_gap_robustness.py

Implement:
  - State-space size sweep (max_bits = 4-8)
  - Method comparison (A vs B vs C)
  - Lag sensitivity
  - Markov order test
  - Leave-one-subject-out cross-validation
```

### 9.4 Phase 4: Multiscale + Cross-Channel (Day 4)

```
File: eeg_experiment/scripts/saturation_gap_multiscale.py

Implement:
  - compute_multiscale_gap() — 5 timescales per epoch
  - compute_cross_channel_gap() — frontal vs occipital vs joint
  - Figure 3 (multiscale profile)
  - Figure 4 (individual trajectories)
```

---

## 10. Validation Strategy

### 10.1 Synthetic Data Validation

Before running on real EEG, validate on synthetic data where the true gap is known:

```python
def validate_on_synthetic():
    """
    Test 1: Detailed balance chain (Σ = 0 exactly)
    → F_actual should be 1.0, Δ should be 0.0

    Test 2: Known irreversible chain (analytical Σ)
    → Verify Σ_estimated matches Σ_analytical
    → Verify F_actual > exp(-Σ/2)

    Test 3: Ornstein-Uhlenbeck process (known closed form)
    → Σ = (drift/diffusion)^2 analytically
    → Compare numerical Δ with analytical prediction
    """
    # Test 1: Reversible chain (detailed balance)
    n = 10
    pi = np.random.dirichlet(np.ones(n))
    # Construct T satisfying detailed balance: pi_i T_ij = pi_j T_ji
    S = np.random.rand(n, n)
    S = (S + S.T) / 2  # symmetric
    T = S * pi[np.newaxis, :] / pi[:, np.newaxis]
    # Normalize rows
    T = T / T.sum(axis=1, keepdims=True)

    sigma = schnakenberg_entropy_production(T, pi)
    assert sigma < 1e-10, f"Reversible chain should have Σ ≈ 0, got {sigma}"

    p0 = np.random.dirichlet(np.ones(n))
    result = compute_saturation_gap(T, pi, p0)
    assert result['F_actual'] > 0.999, f"Reversible chain: F should ≈ 1"
    print("Test 1 PASSED: Reversible chain, Σ=0, F≈1, Δ≈0")

    # Test 2: Irreversible chain
    T_irr = np.random.dirichlet(np.ones(n), size=n)  # random stochastic matrix
    pi_irr = stationary_distribution(T_irr)
    sigma_irr = schnakenberg_entropy_production(T_irr, pi_irr)
    result_irr = compute_saturation_gap(T_irr, pi_irr, pi_irr)

    assert result_irr['F_actual'] >= result_irr['bound'] - 1e-10, \
        f"JRSWW bound violated! F={result_irr['F_actual']}, bound={result_irr['bound']}"
    print(f"Test 2 PASSED: Irreversible chain, Σ={sigma_irr:.4f}, "
          f"F={result_irr['F_actual']:.4f}, bound={result_irr['bound']:.4f}, "
          f"Δ={result_irr['delta']:.4f}")

    # Test 3: Scale test — Δ should increase with irreversibility
    deltas = []
    for strength in np.linspace(0, 2, 20):
        T_test = np.eye(n) * (1 - strength * 0.1) + strength * 0.1 * T_irr
        T_test = T_test / T_test.sum(axis=1, keepdims=True)
        pi_test = stationary_distribution(T_test)
        r = compute_saturation_gap(T_test, pi_test, pi_test)
        deltas.append(r['delta'])

    print(f"Test 3: Δ range = [{min(deltas):.4f}, {max(deltas):.4f}]")
```

### 10.2 Consistency Checks on Real Data

```
Check 1: F_actual ∈ [0, 1] for all epochs
Check 2: Σ ≥ 0 for all epochs
Check 3: Δ ≥ -ε for all epochs (small negative values = numerical error)
Check 4: Δ_normalized ∈ [0, 1] for all epochs
Check 5: π is a valid probability distribution (sums to 1, all ≥ 0)
Check 6: R is a valid stochastic matrix (rows sum to 1)
Check 7: Bound compliance ≥ 95% (small violations from estimation noise)
```

---

## 11. Reporting Template

### 11.1 Results Summary Table

```
| Stage | N_epochs | Σ (median [IQR]) | F_actual (median [IQR]) |
|-------|----------|-------------------|-------------------------|
| Wake  |          |                   |                         |
| N1    |          |                   |                         |
| N2    |          |                   |                         |
| N3    |          |                   |                         |
| REM   |          |                   |                         |

(continued)
| Bound (median) | Δ (median [IQR]) | Δ̄ (median [IQR]) | Compliance (%) |
|----------------|-------------------|---------------------|----------------|
|                |                   |                     |                |
```

### 11.2 Statistical Tests Table

```
| Test | Statistic | p-value | Interpretation |
|------|-----------|---------|----------------|
| Kruskal-Wallis (Δ across stages) | | | |
| Jonckheere-Terpstra (ordered: W>N1>N2>N3) | | | |
| Mann-Whitney (W vs N3) | | | |
| Mann-Whitney (REM vs NREM) | | | |
| Spearman (Δ vs spectral entropy) | | | |
```

---

## 12. Potential Publication Angle

**Title**: "The Petz Recovery Saturation Gap as a Consciousness Signature: First Measurement from Human EEG"

**Abstract** (draft): We present the first empirical measurement of the Petz recovery saturation gap — the difference between the actual neural retrodiction fidelity and the JRSWW information-theoretic lower bound — across human consciousness states. Using 23-subject whole-night polysomnographic EEG, we model neural dynamics as a classical Markov chain and construct the explicit Petz recovery map (Bayesian retrodiction channel). We find that the saturation gap Δ = F_actual - exp(-Σ/2) varies systematically with consciousness level, with waking states showing the largest gap (best recovery relative to the bound) and N2 sleep showing the smallest. The normalized gap Δ̄ provides a tractable, mathematically grounded consciousness metric that may serve as a polynomial-time proxy for integrated information Φ.

**Target journals**: Physical Review E (biophysics section), PLOS Computational Biology, or Journal of Neuroscience Methods.

---

## 13. Known Limitations and Mitigations

| Limitation | Impact | Mitigation |
|-----------|--------|-----------|
| Markov(1) assumption may be too simple | Underestimates memory → biases Σ | Test Markov order; use lag > 1 |
| Only 2 EEG channels | Limited spatial resolution | Frequency-band expansion → 8 virtual channels; cross-channel analysis |
| 100 Hz sampling limits gamma band | Miss high-frequency dynamics | Focus on delta-beta range; note limitation |
| State space discretization is arbitrary | Δ depends on binning | Robustness sweep over max_bits and n_bins |
| Classical (not quantum) treatment | Ignores quantum coherence | Classical approximation is standard for EEG; note as limitation |
| No anesthesia data in Sleep-EDF | Cannot test strongest prediction | Propose as future work with separate dataset |
| Laplace smoothing introduces bias | Biases T toward uniform | Use minimal smoothing (α = 1e-6); test sensitivity |

---

## 14. References

1. Huang, S.-K. (2025). "Retrodiction, Recovery, and the Flow of Time." arXiv preprint. DOI: 10.5281/zenodo.14538923
2. Junge, R., Renner, R., Sutter, D., Wilde, M. M., Winter, A. (2018). "Universal recovery maps and approximate sufficiency of quantum relative entropy." Ann. Henri Poincare 19, 2955-2978.
3. Gao, L., Junge, M., LaRacuente, N. (2023). "Relative entropy decay and complete positivity mixing time." arXiv:2309.02074.
4. Sanz Perl, Y., et al. (2021). "Nonequilibrium brain dynamics as a signature of consciousness." Phys. Rev. E 104, 014411.
5. Petz, D. (1986). "Sufficient subalgebras and the relative entropy of states." Commun. Math. Phys. 105, 123-131.
6. Schnakenberg, J. (1976). "Network theory of microscopic and macroscopic behavior of master equation systems." Rev. Mod. Phys. 48, 571.
7. Kemp, B., et al. (2000). "Analysis of a sleep-dependent neuronal feedback loop." IEEE-BME 47(9):1185-1194.
8. Tononi, G. (2004). "An information integration theory of consciousness." BMC Neuroscience 5, 42.
9. Albantakis, L., et al. (2023). "Integrated information theory (IIT) 4.0." PLOS Comput. Biol.
10. Fawzi, O., Renner, R. (2015). "Quantum conditional mutual information and approximate Markov chains." Commun. Math. Phys. 340, 575-611.

---

## Appendix A: Complete Pseudocode Summary

```
ALGORITHM: Measure Petz Recovery Saturation Gap from EEG

INPUT:
  - EEG recordings (EDF format) with hypnogram annotations
  - Parameters: max_bits=6, lag=1, n_bootstrap=1000

OUTPUT:
  - Table of (subject, stage, Δ, Δ̄, Σ, F_actual, bound) per epoch
  - Statistical tests and figures

PROCEDURE:

1. PREPROCESSING
   FOR each subject:
     a. Load EDF file (mne.io.read_raw_edf)
     b. Pick EEG channels (Fpz-Cz, Pz-Oz)
     c. Bandpass filter 0.5-45 Hz
     d. Segment into 30-second epochs
     e. Parse hypnogram → stage label per epoch
     f. Reject epochs with max |amplitude| > 500 μV

2. STATE SPACE CONSTRUCTION (per epoch)
   a. Expand 2 channels × 4 frequency bands = 8 virtual channels
   b. Binarize each virtual channel (running median threshold)
   c. PCA → max_bits effective bits
   d. Encode state: s(t) = Σ_k bit_k(t) × 2^k
   → state sequence s[0..2999], n_states = 2^max_bits

3. TRANSITION MATRIX ESTIMATION (per epoch)
   a. Count transitions: C[i,j] = #{t : s(t)=i, s(t+lag)=j}
   b. Laplace smoothing: C'[i,j] = C[i,j] + α
   c. Normalize: T[i,j] = C'[i,j] / Σ_j C'[i,j]

4. STATIONARY DISTRIBUTION (per epoch)
   a. Solve π^T T = π^T via eigendecomposition of T^T
   b. Take eigenvector for eigenvalue closest to 1
   c. Normalize: π ← |π| / Σ|π_i|
   d. Enforce positivity: π_i ← max(π_i, 10^{-15})

5. ENTROPY PRODUCTION (per epoch)
   a. P_fwd[i,j] = π[i] T[i,j]
   b. P_bwd[i,j] = π[j] T[j,i]
   c. Σ = Σ_{i,j} P_fwd[i,j] × ln(P_fwd[i,j] / P_bwd[i,j])
   d. Σ ← max(0, Σ)

6. PETZ RECOVERY MAP (per epoch)  ← KEY NEW STEP
   a. π_out[j] = Σ_i π[i] T[i,j]
   b. R[j,i] = π[i] T[i,j] / π_out[j]    (Bayes' theorem)
   c. Verify: R is row-stochastic

7. PETZ FIDELITY (per epoch)
   a. p_initial = empirical distribution from first quarter of epoch
   b. p_evolved = p_initial × T                (forward channel)
   c. p_recovered = R^T × p_evolved             (Petz recovery)
   d. F_actual = (Σ_i √(p_initial[i] × p_recovered[i]))^2

8. SATURATION GAP (per epoch)
   a. bound = exp(-Σ/2)
   b. Δ = F_actual - bound
   c. Δ̄ = Δ / (1 - bound)    if bound < 1

9. AGGREGATION
   a. Group results by (subject, stage)
   b. Compute median, IQR, bootstrap 95% CI per stage
   c. Statistical tests (Kruskal-Wallis, Jonckheere-Terpstra)
   d. Generate figures 1-6

10. VALIDATION
    a. Run synthetic data tests (reversible chain, irreversible chain)
    b. Check bound compliance ≥ 95%
    c. Robustness sweep over max_bits = {4,5,6,7,8}
```

---

## Appendix B: Key Mathematical Identities Used

**Identity 1** (Petz map = Bayes' theorem for classical channels):
$$R_{\pi,T}(j \to i) = \frac{\pi_i \, T_{ij}}{\sum_k \pi_k \, T_{kj}} = P(X_t = i \mid X_{t+1} = j)$$

**Identity 2** (Schnakenberg EP = KL divergence of path measures):
$$\Sigma_{\text{Sch}} = D_{\text{KL}}(P_{\text{fwd}} \| P_{\text{bwd}}) = \sum_{i,j} \pi_i T_{ij} \ln\frac{\pi_i T_{ij}}{\pi_j T_{ji}}$$

**Identity 3** (JRSWW bound, classical specialization):
$$F(p, R_{\pi,T} \circ T(p)) \geq e^{-\Sigma_{\text{path}}/2}$$

**Identity 4** (Classical fidelity):
$$F(p, q) = \left(\sum_i \sqrt{p_i \, q_i}\right)^2$$

**Identity 5** (Saturation gap):
$$\Delta = F_{\text{actual}} - e^{-\Sigma/2} \geq 0$$

**Identity 6** (Non-saturation theorem, Gao-Junge-LaRacuente 2023):
$$\Sigma > 0 \implies \Delta > 0$$

The bound is saturated if and only if $\Sigma = 0$ (Petz 1986), which corresponds to a reversible channel (detailed balance).
