#!/usr/bin/env python3
"""
compute_sigma.py -- Entropy production rate from EEG time series
================================================================

Implements the Schnakenberg entropy production rate (Sigma) and related
irreversibility measures for multichannel EEG, following:

    Sanz Perl et al. (2021) Physical Review E 104, 014112
    "Nonequilibrium brain dynamics as a signature of consciousness"

Theory
------
The brain is modeled as a nonequilibrium Markov process. From binarized
multichannel EEG we construct a discrete state space and estimate the
transition matrix P(i -> j). The Schnakenberg entropy production rate:

    Sigma = (1/2) sum_{i,j} (P_ij pi_i - P_ji pi_j) ln(P_ij pi_i / P_ji pi_j)

where pi is the stationary distribution (left eigenvector of P with eigenvalue 1).
Sigma >= 0 by construction, with equality iff detailed balance holds (equilibrium).

Connection to Petz recovery (Paper 1):
    Sigma quantifies the failure of time-reversal symmetry in neural dynamics.
    In the Petz framework, exp(-Sigma) bounds the fidelity of retrodiction.
    Higher Sigma = more irreversible dynamics = harder to retrodict past states.

Expected ordering (Sanz Perl et al.):
    Sigma: Wakefulness > N1 > N2 > N3
    Consciousness tracks irreversibility / broken detailed balance.

Additional measures
-------------------
- Lempel-Ziv Complexity (LZC): Algorithmic complexity of binarized signal
- Sample Entropy (SampEn): Regularity measure (Richman & Moorman 2000)
- Permutation Entropy (PE): Ordinal pattern statistics (Bandt & Pompe 2002)
- Time-Reversal Classification: ML-based irreversibility test

Usage
-----
    # Synthetic data test
    python compute_sigma.py --synthetic

    # Real EEG from Sleep-EDF
    python compute_sigma.py --edf /path/to/PSG.edf --hypno /path/to/Hypnogram.edf

    # Batch processing
    python compute_sigma.py --data-dir /path/to/sleep-edf/

Author: Sheng-Kai Huang (akai@fawstudio.com)
Date: 2026-03-19
"""

import argparse
import json
import logging
import os
import sys
import warnings
from collections import Counter
from itertools import product as iter_product
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from numpy.typing import NDArray
from scipy import signal as sp_signal
from scipy.linalg import eig
from scipy.spatial.distance import pdist, squareform

# Suppress MNE verbose output
os.environ["MNE_LOGGING_LEVEL"] = "WARNING"

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


# ===========================================================================
# 1. SCHNAKENBERG ENTROPY PRODUCTION RATE
# ===========================================================================

def binarize_channels(
    data: NDArray[np.float64],
    method: str = "median",
) -> NDArray[np.int8]:
    """
    Binarize multichannel EEG: each sample is 0 or 1 per channel.

    Parameters
    ----------
    data : array, shape (n_channels, n_samples)
    method : 'median' (default) or 'mean'

    Returns
    -------
    binary : array, shape (n_channels, n_samples), dtype int8
    """
    if method == "median":
        threshold = np.median(data, axis=1, keepdims=True)
    elif method == "mean":
        threshold = np.mean(data, axis=1, keepdims=True)
    else:
        raise ValueError(f"Unknown binarization method: {method}")
    return (data >= threshold).astype(np.int8)


def reduce_state_space(
    binary: NDArray[np.int8],
    max_bits: int = 8,
    method: str = "pca",
) -> NDArray[np.int64]:
    """
    Reduce the binary multichannel state to a tractable number of dimensions.

    With N channels, the full state space has 2^N states. For N > ~12 this
    becomes intractable. We reduce via PCA on the binary matrix, then re-binarize
    the top components.

    Parameters
    ----------
    binary : array, shape (n_channels, n_samples)
    max_bits : maximum number of bits for state encoding (default 8 -> 256 states)
    method : 'pca' or 'direct' (use first max_bits channels)

    Returns
    -------
    states : array of integer state labels, shape (n_samples,)
    n_states : number of distinct states
    """
    n_channels, n_samples = binary.shape

    if n_channels <= max_bits:
        # Direct encoding -- state space is already tractable
        effective = binary
    elif method == "pca":
        # PCA on the transposed binary matrix
        from sklearn.decomposition import PCA
        pca = PCA(n_components=max_bits)
        projected = pca.fit_transform(binary.T)  # (n_samples, max_bits)
        # Re-binarize the projections around their medians
        med = np.median(projected, axis=0, keepdims=True)
        effective = (projected >= med).astype(np.int8).T  # (max_bits, n_samples)
    elif method == "direct":
        effective = binary[:max_bits, :]
    else:
        raise ValueError(f"Unknown reduction method: {method}")

    n_bits = effective.shape[0]
    # Encode as integer: state = sum_k bit_k * 2^k
    powers = 2 ** np.arange(n_bits, dtype=np.int64)
    states = effective.T @ powers  # (n_samples,)
    n_states = 2 ** n_bits

    return states, n_states


def compute_transition_matrix(
    states: NDArray[np.int64],
    n_states: int,
    lag: int = 1,
) -> NDArray[np.float64]:
    """
    Estimate the transition probability matrix P(i -> j) from a state sequence.

    P[i, j] = Prob(state_{t+lag} = j | state_t = i)

    Parameters
    ----------
    states : array of integer state labels
    n_states : total number of possible states
    lag : time lag for transitions (default 1)

    Returns
    -------
    P : transition matrix, shape (n_states, n_states)
        Row-stochastic: P[i, :].sum() == 1
    """
    counts = np.zeros((n_states, n_states), dtype=np.float64)
    for t in range(len(states) - lag):
        i, j = states[t], states[t + lag]
        counts[i, j] += 1

    # Normalize rows, handling zero-count rows
    row_sums = counts.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1  # avoid division by zero
    P = counts / row_sums

    return P


def stationary_distribution(P: NDArray[np.float64]) -> NDArray[np.float64]:
    """
    Compute the stationary distribution pi of transition matrix P.

    Solves pi @ P = pi, i.e., pi is the left eigenvector with eigenvalue 1.

    Parameters
    ----------
    P : transition matrix, shape (n_states, n_states), row-stochastic

    Returns
    -------
    pi : stationary distribution, shape (n_states,), sums to 1
    """
    n = P.shape[0]

    # Left eigenvector: pi @ P = pi  <=>  P^T @ pi = pi
    eigenvalues, eigenvectors = eig(P.T, left=False, right=True)

    # Find eigenvector closest to eigenvalue 1
    idx = np.argmin(np.abs(eigenvalues - 1.0))
    pi = np.real(eigenvectors[:, idx])

    # Ensure non-negative and normalize
    pi = np.abs(pi)
    pi_sum = pi.sum()
    if pi_sum > 0:
        pi /= pi_sum
    else:
        # Fallback: uniform
        pi = np.ones(n) / n

    return pi


def schnakenberg_entropy_production(
    P: NDArray[np.float64],
    pi: NDArray[np.float64],
    epsilon: float = 1e-20,
) -> float:
    """
    Compute the Schnakenberg entropy production rate.

    Sigma = (1/2) sum_{i,j} (P_ij pi_i - P_ji pi_j) ln(P_ij pi_i / P_ji pi_j)

    This equals zero iff detailed balance holds (P_ij pi_i = P_ji pi_j for all i,j),
    and is strictly positive for any nonequilibrium steady state.

    Parameters
    ----------
    P : transition matrix, shape (n, n)
    pi : stationary distribution, shape (n,)
    epsilon : regularization to avoid log(0)

    Returns
    -------
    sigma : Schnakenberg entropy production rate (nats per transition)
    """
    n = P.shape[0]
    sigma = 0.0

    for i in range(n):
        for j in range(i + 1, n):
            # Forward and backward probability currents
            J_ij = P[i, j] * pi[i]
            J_ji = P[j, i] * pi[j]

            if J_ij < epsilon and J_ji < epsilon:
                continue

            # Regularize
            J_ij_reg = max(J_ij, epsilon)
            J_ji_reg = max(J_ji, epsilon)

            # Schnakenberg formula: sum over pairs (i<j)
            sigma += (J_ij - J_ji) * np.log(J_ij_reg / J_ji_reg)

    return sigma  # the factor 1/2 with sum over all i,j = sum over i<j without 1/2


def schnakenberg_entropy_production_vectorized(
    P: NDArray[np.float64],
    pi: NDArray[np.float64],
    epsilon: float = 1e-20,
) -> float:
    """
    Vectorized version of Schnakenberg Sigma. Faster for large state spaces.

    Sigma = (1/2) sum_{i,j} (J_ij - J_ji) ln(J_ij / J_ji)
    where J_ij = P_ij * pi_i (probability current from i to j).
    """
    # Probability current matrix: J[i,j] = P[i,j] * pi[i]
    J = P * pi[:, np.newaxis]

    # Only upper triangle (avoid double counting)
    J_forward = np.triu(J, k=1)
    J_backward = np.triu(J.T, k=1)

    # Mask out negligible entries
    mask = (J_forward > epsilon) | (J_backward > epsilon)

    J_f = np.where(mask, np.maximum(J_forward, epsilon), epsilon)
    J_b = np.where(mask, np.maximum(J_backward, epsilon), epsilon)

    # Schnakenberg formula
    sigma_matrix = (J_forward - J_backward) * np.log(J_f / J_b)
    sigma_matrix = np.where(mask, sigma_matrix, 0.0)

    return float(sigma_matrix.sum())


# ===========================================================================
# 2. LEMPEL-ZIV COMPLEXITY (LZC)
# ===========================================================================

def lempel_ziv_complexity(binary_sequence: NDArray) -> float:
    """
    Compute the Lempel-Ziv complexity of a binary sequence.

    Counts the number of distinct subsequences encountered when scanning
    left to right. Normalized by n/log2(n) to give a value in [0, 1].

    Reference: Lempel & Ziv (1976), IEEE Trans. Inform. Theory.

    Parameters
    ----------
    binary_sequence : 1D array of 0s and 1s

    Returns
    -------
    lzc : normalized Lempel-Ziv complexity in [0, 1]
    """
    s = binary_sequence.flatten().astype(int)
    n = len(s)
    if n <= 1:
        return 0.0

    # Convert to string for substring operations
    s_str = "".join(map(str, s))

    # Count distinct words
    complexity = 1
    start = 0
    current = 1

    while current < n:
        # Check if substring s[start:current+1] has been seen in s[0:current]
        substr = s_str[start : current + 1]
        if substr not in s_str[:current]:
            complexity += 1
            start = current + 1
            current = start
        else:
            current += 1

    # Normalize
    if n > 0:
        lzc = complexity / (n / max(np.log2(n), 1))
    else:
        lzc = 0.0

    return lzc


def lempel_ziv_complexity_multichannel(
    binary: NDArray[np.int8],
) -> float:
    """
    Compute LZC for multichannel binary EEG.

    Concatenates all channels into a single binary string.

    Parameters
    ----------
    binary : array, shape (n_channels, n_samples)

    Returns
    -------
    lzc : normalized Lempel-Ziv complexity
    """
    # Concatenate channels
    concat = binary.flatten()
    return lempel_ziv_complexity(concat)


# ===========================================================================
# 3. SAMPLE ENTROPY
# ===========================================================================

def sample_entropy(
    x: NDArray[np.float64],
    m: int = 2,
    r: Optional[float] = None,
    r_factor: float = 0.2,
) -> float:
    """
    Compute Sample Entropy of a 1D time series.

    SampEn(m, r) = -ln(A/B) where
    - B = number of template matches of length m
    - A = number of template matches of length m+1
    A "match" means max|x[i+k] - x[j+k]| < r for k = 0, ..., m-1

    Reference: Richman & Moorman (2000), Am. J. Physiol.

    Parameters
    ----------
    x : 1D time series
    m : embedding dimension (default 2)
    r : tolerance threshold (default 0.2 * std(x))
    r_factor : used if r is None: r = r_factor * std(x)

    Returns
    -------
    sampen : sample entropy (nats). Higher = more complex/irregular.
    """
    x = np.asarray(x, dtype=np.float64).flatten()
    N = len(x)

    if r is None:
        r = r_factor * np.std(x)
        if r == 0:
            return 0.0

    def count_templates(m_val):
        """Count matching template pairs of length m_val."""
        templates = np.array([x[i : i + m_val] for i in range(N - m_val)])
        n_templates = len(templates)
        count = 0
        for i in range(n_templates):
            for j in range(i + 1, n_templates):
                if np.max(np.abs(templates[i] - templates[j])) < r:
                    count += 1
        return count

    B = count_templates(m)
    A = count_templates(m + 1)

    if B == 0 or A == 0:
        return float("inf") if A == 0 and B > 0 else 0.0

    return -np.log(A / B)


def sample_entropy_fast(
    x: NDArray[np.float64],
    m: int = 2,
    r: Optional[float] = None,
    r_factor: float = 0.2,
    max_samples: int = 2000,
) -> float:
    """
    Fast approximate Sample Entropy using random subsampling for long series.

    Parameters
    ----------
    x : 1D time series
    m : embedding dimension
    r : tolerance (default 0.2 * std)
    r_factor : if r is None, r = r_factor * std(x)
    max_samples : subsample if series is longer than this

    Returns
    -------
    sampen : approximate sample entropy
    """
    x = np.asarray(x, dtype=np.float64).flatten()
    N = len(x)

    if r is None:
        r = r_factor * np.std(x)
        if r == 0:
            return 0.0

    def count_matches_vectorized(m_val):
        templates = np.array([x[i : i + m_val] for i in range(N - m_val)])
        n = len(templates)
        if n > max_samples:
            idx = np.random.choice(n, max_samples, replace=False)
            templates_sub = templates[idx]
        else:
            templates_sub = templates
            idx = np.arange(n)

        count = 0
        for i in range(len(templates_sub)):
            # Compare template i against all templates after it
            diffs = np.max(np.abs(templates - templates_sub[i]), axis=1)
            count += np.sum(diffs < r) - 1  # exclude self-match
        return max(count, 0)

    B = count_matches_vectorized(m)
    A = count_matches_vectorized(m + 1)

    if B == 0:
        return 0.0
    if A == 0:
        return float("inf")

    return -np.log(A / B)


# ===========================================================================
# 4. PERMUTATION ENTROPY
# ===========================================================================

def permutation_entropy(
    x: NDArray[np.float64],
    order: int = 3,
    delay: int = 1,
    normalize: bool = True,
) -> float:
    """
    Compute Permutation Entropy of a 1D time series.

    For each time point, the ordinal pattern of (x[t], x[t+d], ..., x[t+(m-1)d])
    is determined (the rank ordering). PE is the Shannon entropy of the
    distribution of ordinal patterns.

    Reference: Bandt & Pompe (2002), Phys. Rev. Lett. 88, 174102.

    Parameters
    ----------
    x : 1D time series
    order : embedding dimension m (default 3, giving 3! = 6 patterns)
    delay : time delay d (default 1)
    normalize : if True, divide by log(m!) to get PE in [0, 1]

    Returns
    -------
    pe : permutation entropy (normalized if requested)
    """
    x = np.asarray(x, dtype=np.float64).flatten()
    N = len(x)
    n_patterns = N - (order - 1) * delay

    if n_patterns <= 0:
        return 0.0

    # Extract ordinal patterns
    patterns = []
    for i in range(n_patterns):
        window = x[i : i + order * delay : delay]
        pattern = tuple(np.argsort(window))
        patterns.append(pattern)

    # Count pattern frequencies
    counter = Counter(patterns)
    total = sum(counter.values())

    # Shannon entropy
    pe = 0.0
    for count in counter.values():
        p = count / total
        if p > 0:
            pe -= p * np.log(p)

    if normalize:
        from math import factorial
        max_entropy = np.log(factorial(order))
        if max_entropy > 0:
            pe /= max_entropy

    return pe


# ===========================================================================
# 5. TIME-REVERSAL CLASSIFICATION TEST
# ===========================================================================

def extract_features_for_reversal(
    data: NDArray[np.float64],
    n_segments: int = 50,
    segment_length: int = 500,
) -> Tuple[NDArray, NDArray]:
    """
    Create forward and time-reversed segments for classification.

    Features extracted per segment:
    - Autocovariance at lags 1..5
    - Skewness of increments
    - Third-order statistics (time-asymmetric by nature)

    Parameters
    ----------
    data : array, shape (n_channels, n_samples)
    n_segments : number of forward/backward segment pairs
    segment_length : samples per segment

    Returns
    -------
    X : feature matrix, shape (2 * n_segments, n_features)
    y : labels, 0 = forward, 1 = reversed
    """
    n_channels, n_samples = data.shape
    max_start = n_samples - segment_length

    if max_start <= 0:
        raise ValueError("Data too short for the requested segment length")

    starts = np.random.choice(max_start, size=min(n_segments, max_start), replace=False)

    features_list = []
    labels = []

    for start in starts:
        segment = data[:, start : start + segment_length]

        for direction, label in [("forward", 0), ("backward", 1)]:
            if direction == "backward":
                seg = segment[:, ::-1].copy()
            else:
                seg = segment.copy()

            feats = _extract_temporal_features(seg)
            features_list.append(feats)
            labels.append(label)

    X = np.array(features_list)
    y = np.array(labels)

    return X, y


def _extract_temporal_features(segment: NDArray[np.float64]) -> NDArray[np.float64]:
    """
    Extract temporal features that break time-reversal symmetry.

    These features are asymmetric under t -> -t, making them useful
    for distinguishing forward from backward time series.
    """
    n_ch, n_t = segment.shape
    features = []

    for ch in range(n_ch):
        x = segment[ch]
        dx = np.diff(x)

        # Autocovariance at several lags
        x_centered = x - x.mean()
        var = np.var(x)
        if var > 0:
            for lag in [1, 2, 3, 5, 10]:
                if lag < n_t:
                    acov = np.mean(x_centered[:-lag] * x_centered[lag:]) / var
                    features.append(acov)
                else:
                    features.append(0.0)
        else:
            features.extend([0.0] * 5)

        # Skewness of increments (time-asymmetric)
        if np.std(dx) > 0:
            skew_dx = np.mean((dx - dx.mean()) ** 3) / (np.std(dx) ** 3)
        else:
            skew_dx = 0.0
        features.append(skew_dx)

        # Third-order temporal correlation: <x(t)^2 * x(t+1)> - <x(t) * x(t+1)^2>
        # This is strictly zero for time-reversible Gaussian processes
        if n_t > 1:
            toc = np.mean(x[:-1] ** 2 * x[1:]) - np.mean(x[:-1] * x[1:] ** 2)
            features.append(toc)
        else:
            features.append(0.0)

        # Delayed mutual information asymmetry proxy
        # I(x_t; x_{t+lag}) vs I(x_t; x_{t-lag}) — approximated by correlation
        for lag in [1, 3]:
            if lag < n_t:
                fwd = np.corrcoef(x[:-lag], x[lag:])[0, 1] if np.std(x) > 0 else 0
                features.append(fwd)
            else:
                features.append(0.0)

    return np.array(features, dtype=np.float64)


def time_reversal_classification(
    data: NDArray[np.float64],
    n_segments: int = 100,
    segment_length: int = 500,
    n_splits: int = 5,
) -> Dict[str, float]:
    """
    Train a classifier to distinguish forward vs time-reversed EEG.

    Classification accuracy above 50% indicates irreversibility.
    The further from 50%, the larger the entropy production.

    Parameters
    ----------
    data : array, shape (n_channels, n_samples)
    n_segments : segments per direction
    segment_length : samples per segment
    n_splits : cross-validation folds

    Returns
    -------
    results : dict with 'accuracy_mean', 'accuracy_std', 'irreversibility_index'
    """
    from sklearn.ensemble import GradientBoostingClassifier
    from sklearn.model_selection import StratifiedKFold, cross_val_score
    from sklearn.preprocessing import StandardScaler

    X, y = extract_features_for_reversal(data, n_segments, segment_length)

    # Remove NaN/Inf
    mask = np.all(np.isfinite(X), axis=1)
    X, y = X[mask], y[mask]

    if len(y) < 2 * n_splits:
        return {
            "accuracy_mean": 0.5,
            "accuracy_std": 0.0,
            "irreversibility_index": 0.0,
        }

    # Standardize
    scaler = StandardScaler()
    X = scaler.fit_transform(X)

    clf = GradientBoostingClassifier(
        n_estimators=100,
        max_depth=3,
        random_state=42,
    )

    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)
    scores = cross_val_score(clf, X, y, cv=cv, scoring="accuracy")

    acc_mean = scores.mean()
    acc_std = scores.std()

    # Irreversibility index: 0 = reversible, 1 = maximally irreversible
    irr_index = 2 * abs(acc_mean - 0.5)

    return {
        "accuracy_mean": float(acc_mean),
        "accuracy_std": float(acc_std),
        "irreversibility_index": float(irr_index),
    }


# ===========================================================================
# 6. EPOCH-WISE COMPUTATION PIPELINE
# ===========================================================================

def compute_epoch_measures(
    data: NDArray[np.float64],
    sfreq: float,
    max_bits: int = 6,
    sampen_channel: int = 0,
    pe_order: int = 3,
    use_fast_sampen: bool = True,
) -> Dict[str, float]:
    """
    Compute all entropy/irreversibility measures for a single epoch.

    Parameters
    ----------
    data : array, shape (n_channels, n_samples) -- one epoch
    sfreq : sampling frequency in Hz
    max_bits : max bits for state space reduction
    sampen_channel : which channel to use for SampEn (default 0 = first)
    pe_order : permutation entropy embedding order
    use_fast_sampen : use fast approximate SampEn

    Returns
    -------
    measures : dict with keys 'sigma', 'lzc', 'sampen', 'pe'
    """
    n_channels, n_samples = data.shape

    # 1. Binarize
    binary = binarize_channels(data, method="median")

    # 2. Reduce state space and get state sequence
    states, n_states = reduce_state_space(binary, max_bits=max_bits, method="pca")

    # 3. Transition matrix
    P = compute_transition_matrix(states, n_states, lag=1)

    # 4. Stationary distribution
    pi = stationary_distribution(P)

    # 5. Schnakenberg Sigma
    sigma = schnakenberg_entropy_production_vectorized(P, pi)

    # 6. Lempel-Ziv Complexity
    lzc = lempel_ziv_complexity_multichannel(binary)

    # 7. Sample Entropy (single channel, for speed)
    ch_data = data[min(sampen_channel, n_channels - 1)]
    if use_fast_sampen:
        sampen = sample_entropy_fast(ch_data, m=2, r_factor=0.2, max_samples=1000)
    else:
        sampen = sample_entropy(ch_data, m=2, r_factor=0.2)

    # 8. Permutation Entropy (single channel)
    pe = permutation_entropy(ch_data, order=pe_order, delay=1, normalize=True)

    return {
        "sigma": float(sigma),
        "lzc": float(lzc),
        "sampen": float(sampen) if np.isfinite(sampen) else -1.0,
        "pe": float(pe),
        "n_visited_states": int(len(np.unique(states))),
        "n_total_states": int(n_states),
    }


# ===========================================================================
# 7. EEG DATA LOADING AND PROCESSING
# ===========================================================================

def load_eeg_data(
    edf_path: str,
    hypno_path: Optional[str] = None,
    epoch_duration: float = 30.0,
    pick_channels: Optional[List[str]] = None,
    bandpass: Tuple[float, float] = (0.5, 45.0),
) -> Tuple[NDArray, List[str], float, Optional[NDArray]]:
    """
    Load EEG from EDF file and optionally parse hypnogram.

    Parameters
    ----------
    edf_path : path to PSG .edf file
    hypno_path : path to hypnogram .edf file (optional)
    epoch_duration : epoch length in seconds (default 30)
    pick_channels : list of channel names to select (default: all EEG)
    bandpass : (low_freq, high_freq) for filtering

    Returns
    -------
    epochs : array, shape (n_epochs, n_channels, n_samples_per_epoch)
    ch_names : list of channel names
    sfreq : sampling frequency
    stages : array of sleep stage labels per epoch, or None
    """
    import mne

    logger.info(f"Loading EDF: {edf_path}")
    raw = mne.io.read_raw_edf(edf_path, preload=True, verbose=False)

    # Pick EEG channels
    if pick_channels:
        available = [ch for ch in pick_channels if ch in raw.ch_names]
        if available:
            raw.pick_channels(available)
        else:
            logger.warning(f"Requested channels not found. Using all: {raw.ch_names}")
    else:
        # Try to pick EEG-like channels
        eeg_channels = [ch for ch in raw.ch_names if "EEG" in ch.upper() or
                        "Fpz" in ch or "Pz" in ch or "Cz" in ch or "Fz" in ch or
                        "C3" in ch or "C4" in ch or "O1" in ch or "O2" in ch or
                        "F3" in ch or "F4" in ch]
        if eeg_channels:
            raw.pick_channels(eeg_channels)
        else:
            logger.info(f"No obvious EEG channels found, using all: {raw.ch_names}")

    # Bandpass filter
    raw.filter(bandpass[0], bandpass[1], verbose=False)

    sfreq = raw.info["sfreq"]
    ch_names = raw.ch_names
    data = raw.get_data()  # (n_channels, n_total_samples)

    n_channels, n_total = data.shape
    samples_per_epoch = int(epoch_duration * sfreq)
    n_epochs = n_total // samples_per_epoch

    logger.info(f"  Channels: {n_channels}, sfreq: {sfreq} Hz, "
                f"total: {n_total} samples, epochs: {n_epochs}")

    # Segment into epochs
    epochs = np.zeros((n_epochs, n_channels, samples_per_epoch))
    for i in range(n_epochs):
        start = i * samples_per_epoch
        end = start + samples_per_epoch
        epochs[i] = data[:, start:end]

    # Parse hypnogram if provided
    stages = None
    if hypno_path:
        stages = parse_hypnogram(hypno_path, n_epochs, epoch_duration)

    return epochs, ch_names, sfreq, stages


def parse_hypnogram(
    hypno_path: str,
    n_epochs: int,
    epoch_duration: float = 30.0,
) -> NDArray:
    """
    Parse sleep stage annotations from a hypnogram EDF file.

    Returns array of stage labels: 'W', 'N1', 'N2', 'N3', 'R' (REM).
    """
    import mne

    logger.info(f"Parsing hypnogram: {hypno_path}")
    annot = mne.read_annotations(hypno_path)

    # Map annotations to standard labels
    stage_map = {
        "Sleep stage W": "W",
        "Sleep stage 1": "N1",
        "Sleep stage 2": "N2",
        "Sleep stage 3": "N3",
        "Sleep stage 4": "N3",  # Stage 4 merged into N3
        "Sleep stage R": "R",
        "Sleep stage ?": "?",
        "Movement time": "MT",
    }

    stages = np.array(["?"] * n_epochs, dtype="U3")

    for ann in annot:
        onset = ann["onset"]
        duration = ann["duration"]
        desc = ann["description"]

        stage = stage_map.get(desc, "?")
        start_epoch = int(onset / epoch_duration)
        n_ann_epochs = max(1, int(duration / epoch_duration))

        for e in range(start_epoch, min(start_epoch + n_ann_epochs, n_epochs)):
            stages[e] = stage

    # Count stages
    unique, counts = np.unique(stages, return_counts=True)
    for s, c in zip(unique, counts):
        logger.info(f"  Stage {s}: {c} epochs ({c * epoch_duration / 60:.1f} min)")

    return stages


# ===========================================================================
# 8. SYNTHETIC DATA GENERATION
# ===========================================================================

def generate_synthetic_eeg(
    n_channels: int = 4,
    sfreq: float = 100.0,
    epoch_duration: float = 30.0,
    n_epochs_per_stage: int = 10,
    seed: int = 42,
) -> Tuple[NDArray, float, NDArray]:
    """
    Generate synthetic multi-channel EEG-like data with known properties.

    Different "sleep stages" are simulated with different spectral profiles:
    - W (wake): dominant alpha (8-12 Hz) + beta (12-30 Hz), high complexity
    - N1: reduced alpha, some theta (4-8 Hz)
    - N2: theta dominant, sleep spindles (12-14 Hz bursts)
    - N3: delta (0.5-4 Hz) dominant, low complexity
    - R (REM): mixed frequencies, moderate complexity

    Key: Wake has the most irreversible dynamics (asymmetric waveforms),
    N3 has the least (nearly symmetric slow waves).

    Parameters
    ----------
    n_channels : number of simulated EEG channels
    sfreq : sampling frequency
    epoch_duration : epoch length in seconds
    n_epochs_per_stage : epochs per stage
    seed : random seed

    Returns
    -------
    epochs : array, shape (n_total_epochs, n_channels, samples_per_epoch)
    sfreq : sampling frequency
    stages : array of stage labels
    """
    rng = np.random.RandomState(seed)
    samples_per_epoch = int(epoch_duration * sfreq)
    t = np.arange(samples_per_epoch) / sfreq

    stage_list = ["W", "N1", "N2", "N3", "R"]
    n_total = n_epochs_per_stage * len(stage_list)

    epochs = np.zeros((n_total, n_channels, samples_per_epoch))
    stages = np.empty(n_total, dtype="U3")

    # Stage-specific spectral parameters
    # (alpha_power, beta_power, theta_power, delta_power, noise_level, nonlinearity)
    stage_params = {
        "W":  {"alpha": 1.0, "beta": 0.8, "theta": 0.2, "delta": 0.1,
               "noise": 0.5, "nonlin": 0.4},
        "N1": {"alpha": 0.5, "beta": 0.3, "theta": 0.5, "delta": 0.3,
               "noise": 0.4, "nonlin": 0.25},
        "N2": {"alpha": 0.2, "beta": 0.1, "theta": 0.7, "delta": 0.5,
               "noise": 0.3, "nonlin": 0.15},
        "N3": {"alpha": 0.05, "beta": 0.02, "theta": 0.2, "delta": 1.0,
               "noise": 0.2, "nonlin": 0.05},
        "R":  {"alpha": 0.3, "beta": 0.5, "theta": 0.6, "delta": 0.2,
               "noise": 0.5, "nonlin": 0.3},
    }

    idx = 0
    for stage in stage_list:
        params = stage_params[stage]
        for ep in range(n_epochs_per_stage):
            for ch in range(n_channels):
                # Phase randomization per channel
                phase = rng.uniform(0, 2 * np.pi, 10)

                # Oscillatory components
                alpha = params["alpha"] * np.sin(2 * np.pi * 10 * t + phase[0])
                beta = params["beta"] * np.sin(2 * np.pi * 20 * t + phase[1])
                theta = params["theta"] * np.sin(2 * np.pi * 6 * t + phase[2])
                delta = params["delta"] * np.sin(2 * np.pi * 1.5 * t + phase[3])

                # White noise
                noise = params["noise"] * rng.randn(samples_per_epoch)

                # Combine
                x = alpha + beta + theta + delta + noise

                # Add nonlinearity (asymmetric sawtooth-like distortion)
                # This introduces time-reversal asymmetry
                nl = params["nonlin"]
                x = x + nl * x ** 2 / (1 + np.abs(x))

                # Add inter-channel correlation (spatial mixing)
                if ch > 0:
                    x += 0.3 * epochs[idx, 0, :] * rng.uniform(0.5, 1.0)

                epochs[idx, ch, :] = x

            stages[idx] = stage
            idx += 1

    # Shuffle to avoid order effects in classification
    perm = rng.permutation(n_total)
    epochs = epochs[perm]
    stages = stages[perm]

    logger.info(f"Generated synthetic EEG: {n_total} epochs, "
                f"{n_channels} channels, {sfreq} Hz")

    return epochs, sfreq, stages


# ===========================================================================
# 9. MAIN ANALYSIS PIPELINE
# ===========================================================================

def analyze_epochs(
    epochs: NDArray,
    sfreq: float,
    stages: Optional[NDArray] = None,
    max_bits: int = 6,
    run_time_reversal: bool = True,
    time_reversal_n_segments: int = 50,
) -> Dict:
    """
    Full analysis pipeline: compute all measures for all epochs.

    Parameters
    ----------
    epochs : array, shape (n_epochs, n_channels, samples_per_epoch)
    sfreq : sampling frequency
    stages : optional array of stage labels
    max_bits : max bits for state space reduction
    run_time_reversal : whether to run time-reversal classification
    time_reversal_n_segments : segments per epoch for reversal test

    Returns
    -------
    results : dict with per-epoch measures and stage-wise summaries
    """
    n_epochs, n_channels, n_samples = epochs.shape
    logger.info(f"Analyzing {n_epochs} epochs ({n_channels} ch, "
                f"{n_samples} samples each)...")

    # Per-epoch measures
    epoch_results = []
    for i in range(n_epochs):
        if (i + 1) % 10 == 0 or i == 0:
            logger.info(f"  Epoch {i + 1}/{n_epochs}...")

        measures = compute_epoch_measures(
            epochs[i], sfreq, max_bits=max_bits, use_fast_sampen=True,
        )

        if stages is not None:
            measures["stage"] = str(stages[i])

        epoch_results.append(measures)

    # Time-reversal classification per stage (if requested)
    reversal_results = {}
    if run_time_reversal and stages is not None:
        unique_stages = [s for s in ["W", "N1", "N2", "N3", "R"]
                         if s in stages]
        for stage in unique_stages:
            stage_mask = stages == stage
            stage_epochs = epochs[stage_mask]

            if len(stage_epochs) < 3:
                continue

            # Concatenate epochs for this stage
            stage_data = stage_epochs.reshape(
                stage_epochs.shape[0] * stage_epochs.shape[1] // n_channels,
                n_channels, -1
            )
            # Actually just concatenate samples
            concat_data = np.concatenate(
                [stage_epochs[j] for j in range(len(stage_epochs))],
                axis=1,
            )

            logger.info(f"  Time-reversal test for stage {stage} "
                        f"({concat_data.shape[1]} samples)...")

            seg_len = min(int(5 * sfreq), concat_data.shape[1] // 4)
            n_seg = min(time_reversal_n_segments,
                        concat_data.shape[1] // seg_len - 1)

            if n_seg >= 10 and seg_len >= 100:
                try:
                    rev = time_reversal_classification(
                        concat_data,
                        n_segments=n_seg,
                        segment_length=seg_len,
                    )
                    reversal_results[stage] = rev
                    logger.info(f"    Accuracy: {rev['accuracy_mean']:.3f} "
                                f"+/- {rev['accuracy_std']:.3f}, "
                                f"Irr. index: {rev['irreversibility_index']:.3f}")
                except Exception as e:
                    logger.warning(f"    Time-reversal test failed: {e}")
            else:
                logger.info(f"    Skipping (not enough data: {n_seg} segments "
                            f"of {seg_len} samples)")

    # Stage-wise summary
    stage_summary = {}
    if stages is not None:
        unique_stages = np.unique(stages)
        for stage in unique_stages:
            stage_epochs = [r for r in epoch_results if r.get("stage") == stage]
            if not stage_epochs:
                continue

            summary = {}
            for key in ["sigma", "lzc", "sampen", "pe"]:
                values = [r[key] for r in stage_epochs if r[key] >= 0]
                if values:
                    summary[key] = {
                        "mean": float(np.mean(values)),
                        "std": float(np.std(values)),
                        "median": float(np.median(values)),
                        "n": len(values),
                    }

            if stage in reversal_results:
                summary["time_reversal"] = reversal_results[stage]

            stage_summary[stage] = summary

    return {
        "epoch_results": epoch_results,
        "stage_summary": stage_summary,
        "metadata": {
            "n_epochs": n_epochs,
            "n_channels": n_channels,
            "sfreq": sfreq,
            "max_bits": max_bits,
            "samples_per_epoch": n_samples,
        },
    }


def print_results(results: Dict) -> None:
    """Pretty-print the analysis results."""
    meta = results["metadata"]
    print("\n" + "=" * 70)
    print("EEG ENTROPY PRODUCTION ANALYSIS")
    print("=" * 70)
    print(f"Epochs: {meta['n_epochs']}, Channels: {meta['n_channels']}, "
          f"sfreq: {meta['sfreq']} Hz")
    print(f"State space: {2 ** meta['max_bits']} states "
          f"({meta['max_bits']} bits)")
    print(f"Epoch length: {meta['samples_per_epoch']} samples "
          f"({meta['samples_per_epoch'] / meta['sfreq']:.1f} s)")

    summary = results.get("stage_summary", {})
    if summary:
        print("\n" + "-" * 70)
        print("STAGE-WISE SUMMARY")
        print("-" * 70)

        # Desired ordering
        order = ["W", "R", "N1", "N2", "N3"]
        stages_present = [s for s in order if s in summary]
        stages_other = [s for s in summary if s not in order]
        stages_to_print = stages_present + stages_other

        # Header
        print(f"{'Stage':<8} {'Sigma':>12} {'LZC':>12} {'SampEn':>12} "
              f"{'PermEn':>12} {'TR Acc':>12}")
        print(f"{'':8} {'mean+/-std':>12} {'mean+/-std':>12} {'mean+/-std':>12} "
              f"{'mean+/-std':>12} {'mean+/-std':>12}")
        print("-" * 70)

        for stage in stages_to_print:
            s = summary[stage]
            parts = [f"{stage:<8}"]

            for key in ["sigma", "lzc", "sampen", "pe"]:
                if key in s:
                    m, sd = s[key]["mean"], s[key]["std"]
                    parts.append(f"{m:>6.4f}+/-{sd:<4.4f}")
                else:
                    parts.append(f"{'N/A':>12}")

            if "time_reversal" in s:
                tr = s["time_reversal"]
                parts.append(f"{tr['accuracy_mean']:>5.3f}+/-{tr['accuracy_std']:<4.3f}")
            else:
                parts.append(f"{'N/A':>12}")

            print(" ".join(parts))

        # Check expected ordering: Sigma(W) > Sigma(N1) > Sigma(N2) > Sigma(N3)
        print("\n" + "-" * 70)
        print("IRREVERSIBILITY ORDERING CHECK (Sanz Perl et al. prediction)")
        print("-" * 70)

        expected_order = ["W", "N1", "N2", "N3"]
        present = [s for s in expected_order if s in summary and "sigma" in summary[s]]

        if len(present) >= 2:
            sigmas = [(s, summary[s]["sigma"]["mean"]) for s in present]
            print("Expected: Sigma(W) > Sigma(N1) > Sigma(N2) > Sigma(N3)")
            print("Observed: " + " > ".join(
                [f"Sigma({s})={v:.4f}" for s, v in sigmas]
            ))

            # Check if ordering is correct
            ordered = all(sigmas[i][1] >= sigmas[i + 1][1]
                          for i in range(len(sigmas) - 1))
            if ordered:
                print("RESULT: CONFIRMED -- Irreversibility tracks consciousness level")
            else:
                print("RESULT: ORDERING VIOLATED -- Check data quality or parameters")
        else:
            print("Not enough stages to check ordering.")

    # Connection to Petz recovery
    print("\n" + "-" * 70)
    print("CONNECTION TO PETZ RECOVERY MAP")
    print("-" * 70)
    print("Schnakenberg Sigma = entropy production rate = broken detailed balance")
    print("In the Petz framework (Paper 1):")
    print("  F(rho, Petz[sigma]) >= exp(-Sigma/2)")
    print("  Higher Sigma => lower recovery fidelity => more irreversible")
    print("  Wake: high Sigma => strong time arrow => hard to retrodict")
    print("  N3:   low Sigma  => weak time arrow  => easier retrodiction")
    print("  Consciousness ~ degree of non-equilibrium ~ Sigma")


def save_results(results: Dict, output_path: str) -> None:
    """Save results to JSON file."""

    # Convert numpy types for JSON serialization
    def convert(obj):
        if isinstance(obj, (np.integer,)):
            return int(obj)
        elif isinstance(obj, (np.floating,)):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        return obj

    class NumpyEncoder(json.JSONEncoder):
        def default(self, obj):
            result = convert(obj)
            if result is not obj:
                return result
            return super().default(obj)

    with open(output_path, "w") as f:
        json.dump(results, f, indent=2, cls=NumpyEncoder)

    logger.info(f"Results saved to: {output_path}")


# ===========================================================================
# 10. CLI ENTRY POINT
# ===========================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Compute entropy production rate (Sigma) from EEG time series",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Test with synthetic data
  python compute_sigma.py --synthetic

  # Process a single EDF file
  python compute_sigma.py --edf recording.edf --hypno hypnogram.edf

  # Batch process Sleep-EDF dataset
  python compute_sigma.py --data-dir ./data/sleep-edf/

  # Adjust state space resolution
  python compute_sigma.py --synthetic --max-bits 4
        """,
    )

    parser.add_argument("--synthetic", action="store_true",
                        help="Run with synthetic EEG data")
    parser.add_argument("--edf", type=str, default=None,
                        help="Path to PSG EDF file")
    parser.add_argument("--hypno", type=str, default=None,
                        help="Path to hypnogram EDF file")
    parser.add_argument("--data-dir", type=str, default=None,
                        help="Directory containing Sleep-EDF files")
    parser.add_argument("--output", type=str, default=None,
                        help="Output JSON path (default: auto)")
    parser.add_argument("--max-bits", type=int, default=6,
                        help="Max bits for state space (default 6 -> 64 states)")
    parser.add_argument("--epoch-duration", type=float, default=30.0,
                        help="Epoch duration in seconds (default 30)")
    parser.add_argument("--no-time-reversal", action="store_true",
                        help="Skip time-reversal classification")
    parser.add_argument("--channels", type=str, nargs="+", default=None,
                        help="Channel names to pick (e.g., EEG Fpz-Cz)")

    args = parser.parse_args()

    # Determine output path
    base_dir = Path(__file__).resolve().parent.parent
    results_dir = base_dir / "results"
    results_dir.mkdir(exist_ok=True)

    if args.synthetic:
        # ============== SYNTHETIC TEST ==============
        logger.info("=" * 50)
        logger.info("SYNTHETIC EEG TEST")
        logger.info("=" * 50)

        epochs, sfreq, stages = generate_synthetic_eeg(
            n_channels=4,
            sfreq=100.0,
            epoch_duration=args.epoch_duration,
            n_epochs_per_stage=10,
        )

        results = analyze_epochs(
            epochs, sfreq, stages,
            max_bits=args.max_bits,
            run_time_reversal=not args.no_time_reversal,
        )

        print_results(results)

        output = args.output or str(results_dir / "synthetic_sigma_results.json")
        save_results(results, output)

    elif args.edf:
        # ============== SINGLE EDF FILE ==============
        logger.info("=" * 50)
        logger.info(f"PROCESSING: {args.edf}")
        logger.info("=" * 50)

        epochs, ch_names, sfreq, stages = load_eeg_data(
            args.edf,
            hypno_path=args.hypno,
            epoch_duration=args.epoch_duration,
            pick_channels=args.channels,
        )

        results = analyze_epochs(
            epochs, sfreq, stages,
            max_bits=args.max_bits,
            run_time_reversal=not args.no_time_reversal,
        )

        print_results(results)

        if args.output:
            output = args.output
        else:
            edf_name = Path(args.edf).stem
            output = str(results_dir / f"{edf_name}_sigma_results.json")
        save_results(results, output)

    elif args.data_dir:
        # ============== BATCH PROCESSING ==============
        data_dir = Path(args.data_dir)
        logger.info("=" * 50)
        logger.info(f"BATCH PROCESSING: {data_dir}")
        logger.info("=" * 50)

        # Find PSG files (Sleep-EDF naming convention)
        psg_files = sorted(data_dir.glob("*PSG.edf")) + sorted(data_dir.glob("*-PSG.edf"))
        if not psg_files:
            psg_files = sorted(data_dir.glob("*.edf"))
            psg_files = [f for f in psg_files if "Hypnogram" not in f.name]

        if not psg_files:
            logger.error(f"No EDF files found in {data_dir}")
            sys.exit(1)

        logger.info(f"Found {len(psg_files)} PSG files")

        all_results = {}
        for psg in psg_files:
            # Find matching hypnogram
            stem = psg.stem.replace("-PSG", "").replace("PSG", "")
            hypno_candidates = list(data_dir.glob(f"*{stem}*Hypnogram*.edf"))
            hypno = hypno_candidates[0] if hypno_candidates else None

            logger.info(f"\nProcessing: {psg.name}")
            if hypno:
                logger.info(f"  Hypnogram: {hypno.name}")

            try:
                epochs, ch_names, sfreq, stages = load_eeg_data(
                    str(psg),
                    hypno_path=str(hypno) if hypno else None,
                    epoch_duration=args.epoch_duration,
                    pick_channels=args.channels,
                )

                results = analyze_epochs(
                    epochs, sfreq, stages,
                    max_bits=args.max_bits,
                    run_time_reversal=not args.no_time_reversal,
                )

                print_results(results)
                all_results[psg.name] = results

            except Exception as e:
                logger.error(f"  Failed: {e}")
                all_results[psg.name] = {"error": str(e)}

        output = args.output or str(results_dir / "batch_sigma_results.json")
        save_results(all_results, output)

    else:
        parser.print_help()
        print("\nRun with --synthetic for a quick test!")
        sys.exit(1)

    logger.info("\nDone.")


if __name__ == "__main__":
    main()
