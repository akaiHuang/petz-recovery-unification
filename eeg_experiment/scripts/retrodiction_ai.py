#!/usr/bin/env python3
"""
retrodiction_ai.py -- Neural Petz Recovery Map for EEG Retrodiction
====================================================================

Trains a Transformer-based model to retrodict past EEG states from
present states, implementing an *approximate Petz recovery map*.

Concept
-------
Given multichannel EEG at time t, predict the EEG at time t - Delta_t.
This is "retrodiction" -- inferring the past from the present.

The retrodiction fidelity is defined as:

    F = 1 - MSE_test / Var_test

Paper 1 (Huang 2026) predicts the Petz bound:

    F >= exp(-Sigma / 2)

where Sigma is the Schnakenberg entropy production rate. The bound
states that irreversible dynamics (high Sigma) make retrodiction
harder (lower F), and the Petz recovery map is the *optimal*
retrodiction channel.

Architecture
------------
- Transformer encoder (4 heads, 2 layers, dim=128)
  with positional encoding for temporal structure
- Input:  EEG(t),      shape (n_channels, n_samples)
- Output: EEG(t - dt), shape (n_channels, n_samples)

Baselines
---------
- Linear regression (scikit-learn Ridge)
- Autoregressive model AR(p) via least-squares

The key plot overlays F vs exp(-Sigma/2) for each consciousness
state (Wake, N1, N2, N3, REM), checking whether the Petz bound
holds and how tight it is.

Usage
-----
    # Synthetic validation (random walk with known Sigma)
    python retrodiction_ai.py --synthetic

    # Real EEG from compute_sigma.py results
    python retrodiction_ai.py --edf /path/to/PSG.edf --hypno /path/to/Hypno.edf

    # Load pre-computed sigma results for the plot
    python retrodiction_ai.py --synthetic --sigma-json /path/to/sigma_results.json

Author: Sheng-Kai Huang (akai@fawstudio.com)
Date: 2026-03-19
"""

import argparse
import json
import logging
import math
import os
import sys
import warnings
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from numpy.typing import NDArray

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# PyTorch imports
# ---------------------------------------------------------------------------
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, Dataset, TensorDataset

DEVICE = torch.device("cuda" if torch.cuda.is_available() else
                       "mps" if torch.backends.mps.is_available() else "cpu")
logger.info(f"Using device: {DEVICE}")


# ===========================================================================
# 1. SYNTHETIC DATA: RANDOM WALK WITH KNOWN SIGMA
# ===========================================================================

def generate_ornstein_uhlenbeck(
    n_channels: int = 8,
    n_samples: int = 256,
    n_segments: int = 200,
    dt: float = 1.0 / 256,
    theta: float = 1.0,
    epsilon_coupling: float = 0.3,
    sigma_noise: float = 1.0,
    mu: NDArray | None = None,
    seed: int = 42,
) -> Tuple[NDArray, NDArray, float]:
    """
    Generate multivariate Ornstein-Uhlenbeck process with broken detailed balance.

    The linear Langevin equation is:

        dx = -A (x - mu) dt + sigma * dW

    where A = A_sym + A_anti, with A_sym = theta * I (restoring force)
    and A_anti an antisymmetric nearest-neighbor ring coupling.

    A_sym ensures stability (all eigenvalues of A have positive real part).
    A_anti breaks detailed balance and produces nonzero Sigma.

    Analytical entropy production rate (Seifert 2012, Spinney & Ford 2012):
        Sigma_rate = 2 * Tr(A_anti @ D^{-1} @ C_ss @ A_anti^T)
    where D = sigma^2 * I is the diffusion matrix and C_ss solves
    the Lyapunov equation A @ C_ss + C_ss @ A^T = D.

    Since A_sym = theta * I >> A_anti, C_ss is always positive definite,
    and Sigma >= 0 by construction.

    Parameters
    ----------
    n_channels : number of channels
    n_samples : samples per segment (1 segment = 1 "epoch")
    n_segments : number of consecutive segments
    dt : time step
    theta : mean reversion rate (must be > 0)
    epsilon_coupling : antisymmetric coupling strength (breaks detailed balance)
    sigma_noise : noise amplitude
    mu : mean (default zeros)
    seed : random seed

    Returns
    -------
    segments_input : array (n_segments-1, n_channels, n_samples) -- "present"
    segments_target : array (n_segments-1, n_channels, n_samples) -- "past"
    sigma_true : analytical entropy production rate per segment (nats)
    """
    from scipy.linalg import solve_continuous_lyapunov

    rng = np.random.default_rng(seed)
    if mu is None:
        mu = np.zeros(n_channels)

    total_samples = n_samples * n_segments

    # Build drift matrix: A = A_sym + A_anti
    # A_sym = theta * I  (positive definite, ensures stability)
    A_sym = theta * np.eye(n_channels)

    # A_anti: antisymmetric nearest-neighbor ring coupling
    A_anti = np.zeros((n_channels, n_channels))
    for i in range(n_channels):
        j = (i + 1) % n_channels
        A_anti[i, j] = epsilon_coupling
        A_anti[j, i] = -epsilon_coupling

    A = A_sym + A_anti

    # Verify stability: all eigenvalues of A have positive real part
    eigvals = np.linalg.eigvals(A)
    assert np.all(np.real(eigvals) > 0), (
        f"Drift matrix A is not stable! Eigenvalues: {eigvals}"
    )

    # Diffusion matrix D = sigma^2 * I
    D = sigma_noise**2 * np.eye(n_channels)

    # Steady-state covariance: A @ C_ss + C_ss @ A^T = D
    C_ss = solve_continuous_lyapunov(A, D)

    # Verify C_ss is positive definite
    eigvals_C = np.linalg.eigvalsh(C_ss)
    if np.min(eigvals_C) < 0:
        # Symmetrize (numerical noise)
        C_ss = (C_ss + C_ss.T) / 2
        eigvals_C = np.linalg.eigvalsh(C_ss)
        assert np.min(eigvals_C) > -1e-10, (
            f"C_ss not PSD! Min eigenvalue: {np.min(eigvals_C)}"
        )
        C_ss = C_ss + max(0, -np.min(eigvals_C) + 1e-12) * np.eye(n_channels)

    # Analytical entropy production rate (per unit time):
    # Sigma_rate = 2 * Tr(A_anti @ D^{-1} @ C_ss @ A_anti^T)
    D_inv = np.linalg.inv(D)
    sigma_rate = 2.0 * np.trace(A_anti @ D_inv @ C_ss @ A_anti.T)
    assert sigma_rate >= -1e-10, f"Sigma_rate should be >= 0, got {sigma_rate}"
    sigma_rate = max(0.0, sigma_rate)

    # Convert to per-segment entropy production
    segment_duration = n_samples * dt
    sigma_per_segment = sigma_rate * segment_duration

    # Analytical optimal retrodiction fidelity for the OU process.
    # For a multivariate OU process, the optimal predictor of x(t - T) given
    # the full trajectory x(s) for s in [t, t+T] is the conditional expectation
    # under the posterior. For segment-to-segment retrodiction (predicting
    # the full past segment from the full present segment), the optimal
    # fidelity can be bounded using the propagator exp(-A*T).
    # A simple upper bound: F_opt <= ||exp(-A*T)||_F^2 / n_channels
    # (normalized Frobenius norm of the propagator).
    from scipy.linalg import expm
    propagator = expm(-A * segment_duration)
    # For channel-averaged MSE: the residual variance fraction is
    # approximately 1 - trace(propagator @ propagator.T) / n_channels
    prop_overlap = np.trace(propagator @ propagator.T) / n_channels
    f_optimal_approx = prop_overlap  # approximate upper bound

    logger.info(
        f"    OU process: theta={theta:.2f}, eps={epsilon_coupling:.3f}, "
        f"sigma_rate={sigma_rate:.4f}/s, Sigma/segment={sigma_per_segment:.4f}, "
        f"F_opt_approx={f_optimal_approx:.4f}"
    )

    # Simulate using Euler-Maruyama
    # Generate a long continuous time series, then extract overlapping segments.
    # Total length: enough for n_segments overlapping windows with stride n_samples.
    # Each window has length n_samples; stride = n_samples (non-overlapping pairs).
    total_sim = n_samples * (n_segments + 1)  # +1 for the last pair
    x = np.zeros((n_channels, total_sim))
    x[:, 0] = rng.multivariate_normal(mu, C_ss)
    sqrt_dt = np.sqrt(dt)

    for t in range(1, total_sim):
        drift = -A @ (x[:, t - 1] - mu) * dt
        diffusion = sigma_noise * sqrt_dt * rng.standard_normal(n_channels)
        x[:, t] = x[:, t - 1] + drift + diffusion

    # Extract segment pairs: present and past separated by n_samples steps.
    # segment_i_present = x[:, i*n_samples : (i+1)*n_samples]
    # segment_i_past = x[:, (i-1)*n_samples : i*n_samples]
    segments_present = []
    segments_past = []
    for i in range(1, n_segments + 1):
        start_present = i * n_samples
        end_present = (i + 1) * n_samples
        start_past = (i - 1) * n_samples
        end_past = i * n_samples
        if end_present <= total_sim:
            segments_present.append(x[:, start_present:end_present])
            segments_past.append(x[:, start_past:end_past])

    segments_input = np.array(segments_present)   # (n_pairs, n_channels, n_samples)
    segments_target = np.array(segments_past)

    return segments_input, segments_target, sigma_per_segment, f_optimal_approx


def generate_synthetic_consciousness_states(
    n_channels: int = 8,
    n_samples: int = 256,
    n_segments_per_state: int = 200,
    seed: int = 42,
) -> Dict[str, Dict[str, Any]]:
    """
    Generate synthetic data mimicking different consciousness states.

    Each state has a different level of irreversibility (Sigma),
    controlled by the antisymmetric coupling epsilon:
    - Wake:  high coupling -> high Sigma (most irreversible)
    - N1:    moderate coupling
    - N2:    lower coupling
    - N3:    very low coupling -> near equilibrium (most reversible)
    - REM:   moderate-high coupling (dream state)

    The ordering Sigma: W > R > N1 > N2 > N3 matches empirical
    findings of Sanz Perl et al. (2021 PRE).

    Returns
    -------
    states : dict mapping state_name -> {
        'input': array (n_pairs, n_channels, n_samples),
        'target': array (n_pairs, n_channels, n_samples),
        'sigma_true': float (>= 0),
    }
    """
    # Coupling epsilon controls Sigma; theta controls autocorrelation time.
    #
    # Key physics: for the Petz bound to be tight, the retrodiction must be
    # nontrivial (high F) and Sigma must be modest. With OU process:
    # - Autocorrelation time ~ 1/theta. Segments of duration T = n_samples*dt.
    # - If T << 1/theta, consecutive segments are highly correlated -> high F.
    # - Sigma = sigma_rate * T = 2*Tr(A_anti D^-1 C_ss A_anti^T) * T.
    #
    # We use theta = 0.2 (autocorrelation time = 5s) with T = 1s segments,
    # so exp(-theta*T) = exp(-0.2) ~ 0.82 -> theoretical max F ~ 0.82.
    # Then epsilon controls how much of that is irreversible (Sigma).
    state_params = {
        "W":  {"theta": 0.2, "epsilon_coupling": 0.18, "sigma_noise": 1.0},
        "R":  {"theta": 0.2, "epsilon_coupling": 0.14, "sigma_noise": 1.0},
        "N1": {"theta": 0.2, "epsilon_coupling": 0.10, "sigma_noise": 1.0},
        "N2": {"theta": 0.2, "epsilon_coupling": 0.06, "sigma_noise": 1.0},
        "N3": {"theta": 0.2, "epsilon_coupling": 0.02, "sigma_noise": 1.0},
    }

    states = {}
    for i, (name, params) in enumerate(state_params.items()):
        logger.info(f"  Generating state {name}...")
        inp, tgt, sigma_true, f_opt = generate_ornstein_uhlenbeck(
            n_channels=n_channels,
            n_samples=n_samples,
            n_segments=n_segments_per_state,
            dt=1.0 / n_samples,
            theta=params["theta"],
            epsilon_coupling=params["epsilon_coupling"],
            sigma_noise=params["sigma_noise"],
            seed=seed + i * 1000,
        )

        states[name] = {
            "input": inp,
            "target": tgt,
            "sigma_true": sigma_true,
            "f_optimal": f_opt,
        }

        logger.info(
            f"  State {name}: Sigma = {sigma_true:.4f}, F_opt = {f_opt:.4f}, "
            f"  eps={params['epsilon_coupling']:.2f}, "
            f"  {inp.shape[0]} pairs"
        )

    return states


# ===========================================================================
# 2. DATASET AND DATA LOADING
# ===========================================================================

class EEGRetrodictionDataset(Dataset):
    """
    Dataset of (present_EEG, past_EEG) pairs for retrodiction training.
    """

    def __init__(
        self,
        inputs: NDArray,
        targets: NDArray,
        normalize: bool = True,
    ):
        """
        Parameters
        ----------
        inputs : array (n_pairs, n_channels, n_samples) -- present EEG
        targets : array (n_pairs, n_channels, n_samples) -- past EEG
        normalize : if True, z-score normalize per channel
        """
        assert inputs.shape == targets.shape
        self.n_pairs, self.n_channels, self.n_samples = inputs.shape

        if normalize:
            # Compute stats from inputs (present), apply to both
            all_data = np.concatenate([inputs, targets], axis=0)
            self.channel_mean = all_data.mean(axis=(0, 2), keepdims=True)  # (1, C, 1)
            self.channel_std = all_data.std(axis=(0, 2), keepdims=True) + 1e-8
            inputs = (inputs - self.channel_mean) / self.channel_std
            targets = (targets - self.channel_mean) / self.channel_std
        else:
            self.channel_mean = np.zeros((1, self.n_channels, 1))
            self.channel_std = np.ones((1, self.n_channels, 1))

        # Shape for Transformer: (n_pairs, n_samples, n_channels)
        # Transformer expects (batch, seq_len, features)
        self.inputs = torch.FloatTensor(inputs.transpose(0, 2, 1))
        self.targets = torch.FloatTensor(targets.transpose(0, 2, 1))

        # Store variance of targets for fidelity calculation
        self.target_var = self.targets.var().item()

    def __len__(self):
        return self.n_pairs

    def __getitem__(self, idx):
        return self.inputs[idx], self.targets[idx]


def train_test_split(
    dataset: EEGRetrodictionDataset,
    test_fraction: float = 0.2,
    seed: int = 42,
) -> Tuple[DataLoader, DataLoader, float]:
    """
    Split dataset into train/test and return DataLoaders.

    Returns
    -------
    train_loader, test_loader, target_variance
    """
    n = len(dataset)
    n_test = max(1, int(n * test_fraction))
    n_train = n - n_test

    generator = torch.Generator().manual_seed(seed)
    train_set, test_set = torch.utils.data.random_split(
        dataset, [n_train, n_test], generator=generator
    )

    train_loader = DataLoader(train_set, batch_size=32, shuffle=True)
    test_loader = DataLoader(test_set, batch_size=64, shuffle=False)

    # Compute test target variance
    test_targets = []
    for idx in test_set.indices:
        _, tgt = dataset[idx]
        test_targets.append(tgt)
    test_targets = torch.stack(test_targets)
    test_var = test_targets.var().item()

    return train_loader, test_loader, test_var


# ===========================================================================
# 3. TRANSFORMER MODEL (Approximate Petz Recovery Map)
# ===========================================================================

class PositionalEncoding(nn.Module):
    """Standard sinusoidal positional encoding."""

    def __init__(self, d_model: int, max_len: int = 512, dropout: float = 0.1):
        super().__init__()
        self.dropout = nn.Dropout(p=dropout)

        pe = torch.zeros(max_len, d_model)
        position = torch.arange(0, max_len, dtype=torch.float).unsqueeze(1)
        div_term = torch.exp(
            torch.arange(0, d_model, 2).float() * (-math.log(10000.0) / d_model)
        )
        pe[:, 0::2] = torch.sin(position * div_term)
        pe[:, 1::2] = torch.cos(position * div_term)
        pe = pe.unsqueeze(0)  # (1, max_len, d_model)
        self.register_buffer("pe", pe)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """x: (batch, seq_len, d_model)"""
        x = x + self.pe[:, : x.size(1)]
        return self.dropout(x)


class PetzRecoveryTransformer(nn.Module):
    """
    Transformer-based approximate Petz recovery map.

    Maps EEG(t) -> EEG(t - dt): retrodiction of past brain states.

    Architecture:
        Input projection: Linear(n_channels, d_model)
        Positional encoding
        Transformer encoder (n_layers, n_heads)
        Output projection: Linear(d_model, n_channels)
    """

    def __init__(
        self,
        n_channels: int = 8,
        d_model: int = 128,
        n_heads: int = 4,
        n_layers: int = 2,
        dim_feedforward: int = 256,
        dropout: float = 0.1,
        max_len: int = 512,
    ):
        super().__init__()
        self.n_channels = n_channels
        self.d_model = d_model

        # Input/output projections
        self.input_proj = nn.Linear(n_channels, d_model)
        self.pos_encoder = PositionalEncoding(d_model, max_len, dropout)

        # Transformer encoder layers
        encoder_layer = nn.TransformerEncoderLayer(
            d_model=d_model,
            nhead=n_heads,
            dim_feedforward=dim_feedforward,
            dropout=dropout,
            batch_first=True,
            activation="gelu",
        )
        self.transformer_encoder = nn.TransformerEncoder(
            encoder_layer, num_layers=n_layers
        )

        # Output projection
        self.output_proj = nn.Sequential(
            nn.LayerNorm(d_model),
            nn.Linear(d_model, n_channels),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        x : (batch, seq_len, n_channels) -- present EEG
        returns : (batch, seq_len, n_channels) -- retrodicted past EEG

        Uses a skip connection: output = x + delta, where delta is the
        learned correction. This is natural because the past and present
        are highly correlated -- the model learns the *change*.
        """
        # Project to model dimension
        h = self.input_proj(x)          # (B, T, d_model)
        h = self.pos_encoder(h)         # (B, T, d_model)

        # Transformer encoder (self-attention across time)
        h = self.transformer_encoder(h)  # (B, T, d_model)

        # Project back to channel space (learns the correction)
        delta = self.output_proj(h)      # (B, T, n_channels)

        # Skip connection: past ~ present + learned correction
        return x + delta

    def count_parameters(self) -> int:
        return sum(p.numel() for p in self.parameters() if p.requires_grad)


# ===========================================================================
# 4. BASELINE MODELS
# ===========================================================================

class LinearRetrodiction:
    """
    Linear regression baseline: Ridge regression from present to past.
    Flattens (n_samples, n_channels) -> single vector.
    """

    def __init__(self, alpha: float = 1.0):
        from sklearn.linear_model import Ridge
        self.model = Ridge(alpha=alpha)
        self.fitted = False

    def fit(self, train_loader: DataLoader) -> None:
        inputs, targets = [], []
        for x, y in train_loader:
            inputs.append(x.reshape(x.size(0), -1).numpy())
            targets.append(y.reshape(y.size(0), -1).numpy())
        X = np.concatenate(inputs)
        Y = np.concatenate(targets)
        self.model.fit(X, Y)
        self.fitted = True
        self.output_shape = (x.size(1), x.size(2))

    def predict(self, test_loader: DataLoader) -> Tuple[NDArray, NDArray]:
        """Returns (predictions, targets) both as numpy arrays."""
        assert self.fitted
        preds, trues = [], []
        for x, y in test_loader:
            x_flat = x.reshape(x.size(0), -1).numpy()
            pred = self.model.predict(x_flat)
            preds.append(pred)
            trues.append(y.reshape(y.size(0), -1).numpy())
        return np.concatenate(preds), np.concatenate(trues)

    def evaluate(self, test_loader: DataLoader) -> Tuple[float, float]:
        """Returns (MSE, fidelity)."""
        preds, trues = self.predict(test_loader)
        mse = np.mean((preds - trues) ** 2)
        var = np.var(trues)
        fidelity = 1.0 - mse / var if var > 0 else 0.0
        return mse, fidelity


class ARRetrodiction:
    """
    Autoregressive baseline: fit AR(p) model per channel via least-squares.
    Uses the *present* segment to predict the *past* segment by fitting
    a channel-wise linear map with temporal lags.
    """

    def __init__(self, order: int = 5):
        self.order = order
        self.weights = None  # per-channel weights

    def fit(self, train_loader: DataLoader) -> None:
        # Collect all data
        inputs_all, targets_all = [], []
        for x, y in train_loader:
            inputs_all.append(x.numpy())   # (B, T, C)
            targets_all.append(y.numpy())
        inputs_all = np.concatenate(inputs_all)   # (N, T, C)
        targets_all = np.concatenate(targets_all)

        n_pairs, n_samples, n_channels = inputs_all.shape
        p = self.order
        self.weights = []

        for ch in range(n_channels):
            # Build design matrix from present signal for this channel
            # Use lagged features from the present to predict the past
            X_rows, Y_rows = [], []
            for i in range(n_pairs):
                for t in range(p, n_samples):
                    # Features: present channel values at lags 0..p-1
                    features = inputs_all[i, t - p : t, ch][::-1].copy()
                    X_rows.append(features)
                    Y_rows.append(targets_all[i, t, ch])

            X = np.array(X_rows)
            Y = np.array(Y_rows)

            # Add bias column
            X_bias = np.column_stack([X, np.ones(len(X))])

            # Least-squares solve
            w, _, _, _ = np.linalg.lstsq(X_bias, Y, rcond=None)
            self.weights.append(w)

    def predict(self, test_loader: DataLoader) -> Tuple[NDArray, NDArray]:
        inputs_all, targets_all = [], []
        for x, y in test_loader:
            inputs_all.append(x.numpy())
            targets_all.append(y.numpy())
        inputs_all = np.concatenate(inputs_all)
        targets_all = np.concatenate(targets_all)

        n_pairs, n_samples, n_channels = inputs_all.shape
        p = self.order
        preds = np.zeros_like(targets_all)

        for ch in range(n_channels):
            w = self.weights[ch]
            for i in range(n_pairs):
                for t in range(p, n_samples):
                    features = inputs_all[i, t - p : t, ch][::-1].copy()
                    features_bias = np.append(features, 1.0)
                    preds[i, t, ch] = features_bias @ w
                # Fill first p steps with simple copy
                preds[i, :p, ch] = inputs_all[i, :p, ch]

        preds_flat = preds.reshape(-1)
        trues_flat = targets_all.reshape(-1)
        return preds_flat, trues_flat

    def evaluate(self, test_loader: DataLoader) -> Tuple[float, float]:
        preds, trues = self.predict(test_loader)
        mse = np.mean((preds - trues) ** 2)
        var = np.var(trues)
        fidelity = 1.0 - mse / var if var > 0 else 0.0
        return mse, fidelity


# ===========================================================================
# 5. TRAINING LOOP
# ===========================================================================

def train_transformer(
    model: PetzRecoveryTransformer,
    train_loader: DataLoader,
    test_loader: DataLoader,
    n_epochs: int = 50,
    lr: float = 1e-3,
    weight_decay: float = 1e-4,
    patience: int = 10,
    verbose: bool = True,
) -> Dict[str, List[float]]:
    """
    Train the Transformer retrodiction model.

    Returns
    -------
    history : dict with 'train_loss', 'test_loss' per epoch
    """
    model = model.to(DEVICE)
    optimizer = optim.AdamW(model.parameters(), lr=lr, weight_decay=weight_decay)
    scheduler = optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=n_epochs)
    criterion = nn.MSELoss()

    history = {"train_loss": [], "test_loss": []}
    best_test_loss = float("inf")
    best_state = None
    no_improve = 0

    for epoch in range(1, n_epochs + 1):
        # --- Train ---
        model.train()
        train_loss_sum = 0.0
        n_batches = 0
        for x, y in train_loader:
            x, y = x.to(DEVICE), y.to(DEVICE)
            optimizer.zero_grad()
            pred = model(x)
            loss = criterion(pred, y)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            optimizer.step()
            train_loss_sum += loss.item()
            n_batches += 1

        train_loss = train_loss_sum / n_batches
        history["train_loss"].append(train_loss)

        # --- Evaluate ---
        model.eval()
        test_loss_sum = 0.0
        n_test = 0
        with torch.no_grad():
            for x, y in test_loader:
                x, y = x.to(DEVICE), y.to(DEVICE)
                pred = model(x)
                loss = criterion(pred, y)
                test_loss_sum += loss.item() * x.size(0)
                n_test += x.size(0)

        test_loss = test_loss_sum / n_test
        history["test_loss"].append(test_loss)

        scheduler.step()

        # Early stopping
        if test_loss < best_test_loss:
            best_test_loss = test_loss
            best_state = {k: v.cpu().clone() for k, v in model.state_dict().items()}
            no_improve = 0
        else:
            no_improve += 1

        if verbose and (epoch % 5 == 0 or epoch == 1):
            logger.info(
                f"  Epoch {epoch:3d}/{n_epochs}: "
                f"train_loss={train_loss:.6f}, test_loss={test_loss:.6f}"
            )

        if no_improve >= patience:
            logger.info(f"  Early stopping at epoch {epoch} (patience={patience})")
            break

    # Load best model
    if best_state is not None:
        model.load_state_dict(best_state)
    model = model.to(DEVICE)

    return history


def evaluate_transformer(
    model: PetzRecoveryTransformer,
    test_loader: DataLoader,
    test_var: float,
) -> Tuple[float, float]:
    """
    Evaluate retrodiction fidelity.

    Returns
    -------
    mse : mean squared error on test set
    fidelity : F = 1 - MSE / Var
    """
    model.eval()
    model = model.to(DEVICE)
    criterion = nn.MSELoss(reduction="sum")

    total_mse = 0.0
    total_elements = 0

    with torch.no_grad():
        for x, y in test_loader:
            x, y = x.to(DEVICE), y.to(DEVICE)
            pred = model(x)
            total_mse += criterion(pred, y).item()
            total_elements += y.numel()

    mse = total_mse / total_elements
    fidelity = 1.0 - mse / test_var if test_var > 0 else 0.0
    return mse, fidelity


# ===========================================================================
# 6. FULL PIPELINE: TRAIN + EVALUATE PER STATE
# ===========================================================================

def run_retrodiction_pipeline(
    states: Dict[str, Dict[str, Any]],
    n_epochs: int = 50,
    lr: float = 1e-3,
    d_model: int = 128,
    n_heads: int = 4,
    n_layers: int = 2,
    verbose: bool = True,
) -> Dict[str, Dict[str, float]]:
    """
    For each consciousness state:
    1. Train Transformer retrodiction model
    2. Train Linear and AR baselines
    3. Compute retrodiction fidelity F
    4. Compare with Petz bound exp(-Sigma/2)

    Returns
    -------
    results : dict mapping state_name -> {
        'sigma': float,
        'petz_bound': float,           # exp(-sigma/2)
        'fidelity_transformer': float,
        'fidelity_linear': float,
        'fidelity_ar': float,
        'mse_transformer': float,
        'mse_linear': float,
        'mse_ar': float,
        'bound_satisfied_transformer': bool,
        'bound_satisfied_linear': bool,
        'bound_satisfied_ar': bool,
    }
    """
    results = {}

    for state_name, state_data in states.items():
        logger.info(f"\n{'='*60}")
        logger.info(f"Processing state: {state_name}")
        logger.info(f"{'='*60}")

        inputs = state_data["input"]
        targets = state_data["target"]
        sigma = state_data["sigma_true"]
        f_optimal = state_data.get("f_optimal", None)
        petz_bound = math.exp(-sigma / 2.0)

        n_channels = inputs.shape[1]
        n_samples = inputs.shape[2]

        logger.info(
            f"  Data: {inputs.shape[0]} pairs, "
            f"{n_channels} channels, {n_samples} samples"
        )
        logger.info(
            f"  Sigma = {sigma:.4f}, Petz bound = {petz_bound:.4f}"
            + (f", F_opt = {f_optimal:.4f}" if f_optimal else "")
        )

        # Create dataset
        dataset = EEGRetrodictionDataset(inputs, targets, normalize=True)
        train_loader, test_loader, test_var = train_test_split(dataset, test_fraction=0.2)

        logger.info(f"  Test variance = {test_var:.4f}")

        # --- Transformer ---
        logger.info(f"\n  [Transformer] Training...")
        model = PetzRecoveryTransformer(
            n_channels=n_channels,
            d_model=d_model,
            n_heads=n_heads,
            n_layers=n_layers,
            dim_feedforward=d_model * 2,
            max_len=n_samples,
        )
        logger.info(f"  Parameters: {model.count_parameters():,}")

        history = train_transformer(
            model, train_loader, test_loader,
            n_epochs=n_epochs, lr=lr, verbose=verbose
        )
        mse_tf, fid_tf = evaluate_transformer(model, test_loader, test_var)
        logger.info(
            f"  [Transformer] MSE={mse_tf:.6f}, "
            f"F={fid_tf:.4f}, bound={petz_bound:.4f}, "
            f"satisfied={'YES' if fid_tf >= petz_bound else 'NO'}"
        )

        # --- Linear baseline ---
        logger.info(f"\n  [Linear] Training...")
        linear_model = LinearRetrodiction(alpha=1.0)
        linear_model.fit(train_loader)
        mse_lin, fid_lin = linear_model.evaluate(test_loader)
        logger.info(
            f"  [Linear] MSE={mse_lin:.6f}, "
            f"F={fid_lin:.4f}, bound={petz_bound:.4f}, "
            f"satisfied={'YES' if fid_lin >= petz_bound else 'NO'}"
        )

        # --- AR baseline ---
        logger.info(f"\n  [AR(5)] Training...")
        ar_model = ARRetrodiction(order=5)
        ar_model.fit(train_loader)
        mse_ar, fid_ar = ar_model.evaluate(test_loader)
        logger.info(
            f"  [AR(5)] MSE={mse_ar:.6f}, "
            f"F={fid_ar:.4f}, bound={petz_bound:.4f}, "
            f"satisfied={'YES' if fid_ar >= petz_bound else 'NO'}"
        )

        results[state_name] = {
            "sigma": sigma,
            "petz_bound": petz_bound,
            "f_optimal": f_optimal if f_optimal is not None else -1.0,
            "fidelity_transformer": fid_tf,
            "fidelity_linear": fid_lin,
            "fidelity_ar": fid_ar,
            "mse_transformer": mse_tf,
            "mse_linear": mse_lin,
            "mse_ar": mse_ar,
            "bound_satisfied_transformer": bool(fid_tf >= petz_bound),
            "bound_satisfied_linear": bool(fid_lin >= petz_bound),
            "bound_satisfied_ar": bool(fid_ar >= petz_bound),
        }

    return results


# ===========================================================================
# 7. REAL EEG DATA LOADING (from Sleep-EDF)
# ===========================================================================

def load_real_eeg_states(
    edf_path: str,
    hypno_path: str,
    n_samples: int = 256,
    max_channels: int = 8,
) -> Dict[str, Dict[str, Any]]:
    """
    Load real EEG from Sleep-EDF and organize by sleep stage.

    Parameters
    ----------
    edf_path : path to PSG .edf file
    hypno_path : path to Hypnogram .edf file
    n_samples : samples per segment (epoch)
    max_channels : max EEG channels to use

    Returns
    -------
    states : dict in same format as generate_synthetic_consciousness_states
    """
    try:
        import mne
    except ImportError:
        logger.error("MNE-Python required for real EEG. Install with: pip install mne")
        sys.exit(1)

    mne.set_log_level("WARNING")

    # Load raw EEG
    raw = mne.io.read_raw_edf(edf_path, preload=True, verbose=False)

    # Pick EEG channels only
    eeg_picks = mne.pick_types(raw.info, eeg=True)
    if len(eeg_picks) > max_channels:
        eeg_picks = eeg_picks[:max_channels]
    raw.pick(eeg_picks)

    sfreq = raw.info["sfreq"]
    # Resample if needed to get desired n_samples per second
    if sfreq != n_samples:
        raw.resample(n_samples, verbose=False)

    data = raw.get_data()  # (n_channels, total_samples)
    n_channels = data.shape[0]

    # Parse hypnogram
    annot = mne.read_annotations(hypno_path)
    stage_map = {
        "Sleep stage W": "W",
        "Sleep stage 1": "N1",
        "Sleep stage 2": "N2",
        "Sleep stage 3": "N3",
        "Sleep stage 4": "N3",
        "Sleep stage R": "R",
    }

    # Build epoch -> stage mapping (30-second epochs)
    epoch_duration = 30  # seconds
    n_total_samples = data.shape[1]
    n_epochs_total = n_total_samples // (n_samples * epoch_duration)

    # Segment data into 1-second windows
    n_windows = n_total_samples // n_samples
    windows = data[:, : n_windows * n_samples].reshape(
        n_channels, n_windows, n_samples
    ).transpose(1, 0, 2)  # (n_windows, n_channels, n_samples)

    # Map each 1-second window to a sleep stage
    window_stages = []
    for w_idx in range(n_windows):
        t_sec = w_idx * (n_samples / n_samples)  # = w_idx seconds
        # Find which annotation covers this time
        stage = "?"
        for ann in annot:
            if ann["onset"] <= t_sec < ann["onset"] + ann["duration"]:
                stage = stage_map.get(ann["description"], "?")
                break
        window_stages.append(stage)

    window_stages = np.array(window_stages)

    # Organize by state
    states = {}
    for state_name in ["W", "N1", "N2", "N3", "R"]:
        mask = window_stages == state_name
        state_windows = windows[mask]

        if len(state_windows) < 10:
            logger.warning(
                f"  State {state_name}: only {len(state_windows)} windows, skipping"
            )
            continue

        # Create consecutive pairs (input=present, target=past)
        # Find consecutive windows within the same state
        state_indices = np.where(mask)[0]
        consecutive_pairs = []
        for i in range(len(state_indices) - 1):
            if state_indices[i + 1] == state_indices[i] + 1:
                consecutive_pairs.append(
                    (state_indices[i + 1], state_indices[i])
                )

        if len(consecutive_pairs) < 5:
            logger.warning(
                f"  State {state_name}: only {len(consecutive_pairs)} "
                f"consecutive pairs, skipping"
            )
            continue

        input_idx = [p[0] for p in consecutive_pairs]
        target_idx = [p[1] for p in consecutive_pairs]
        inputs = windows[input_idx]
        targets = windows[target_idx]

        # Sigma will be loaded from compute_sigma.py results
        states[state_name] = {
            "input": inputs,
            "target": targets,
            "sigma_true": 0.0,  # placeholder -- to be filled from sigma results
        }
        logger.info(
            f"  State {state_name}: {len(consecutive_pairs)} pairs, "
            f"{n_channels} channels"
        )

    return states


# ===========================================================================
# 8. VISUALIZATION
# ===========================================================================

def plot_petz_bound(
    results: Dict[str, Dict[str, float]],
    output_path: str = "petz_bound_verification.png",
    title: str = "Petz Recovery Bound: F vs exp(-Sigma/2)",
) -> None:
    """
    The key plot: retrodiction fidelity vs entropy production.

    x-axis: Sigma (entropy production)
    y-axis: F (retrodiction fidelity)
    Overlay: F = exp(-Sigma/2) curve (Petz bound)

    Points above the curve -> bound satisfied.
    Points near the curve -> brain is near-optimal retrodiction machine.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(16, 7))

    # ---- Left panel: F vs Sigma ----
    ax = axes[0]

    # Petz bound curve
    sigma_max = max(r["sigma"] for r in results.values())
    sigma_range = np.linspace(0, sigma_max * 1.3, 200)
    petz_curve = np.exp(-sigma_range / 2)
    ax.plot(sigma_range, petz_curve, "k-", linewidth=2.5,
            label=r"$F = e^{-\Sigma/2}$ (Petz bound)")
    ax.fill_between(sigma_range, 0, petz_curve, alpha=0.08, color="red",
                    label="Bound violation region")

    # State markers
    colors = {"W": "#e74c3c", "N1": "#e67e22", "N2": "#f1c40f",
              "N3": "#2ecc71", "R": "#3498db"}
    markers = {"W": "o", "N1": "s", "N2": "D", "N3": "^", "R": "v"}

    state_order = ["W", "N1", "N2", "N3", "R"]

    # Plot each state
    for state in state_order:
        if state not in results:
            continue
        r = results[state]
        sigma = r["sigma"]
        c = colors.get(state, "gray")
        m = markers.get(state, "o")

        # Transformer (large filled)
        ax.scatter(sigma, r["fidelity_transformer"], c=c, marker=m, s=180,
                   edgecolors="black", linewidth=1.5, zorder=5,
                   label=f"{state} Transformer: F={r['fidelity_transformer']:.3f}")
        # Linear (small filled, semi-transparent)
        ax.scatter(sigma, r["fidelity_linear"], c=c, marker=m, s=70,
                   edgecolors="black", linewidth=0.5, alpha=0.5, zorder=4)
        # AR (hollow markers)
        ax.scatter(sigma, r["fidelity_ar"], marker=m, s=70,
                   edgecolors=c, linewidth=1.5, alpha=0.6, zorder=3,
                   facecolors="none")

        # Annotate sigma value
        ax.annotate(f"{state}", (sigma, r["fidelity_transformer"]),
                    textcoords="offset points", xytext=(8, 8), fontsize=9,
                    fontweight="bold", color=c)

    ax.set_xlabel(r"$\Sigma$ (entropy production per segment)", fontsize=14)
    ax.set_ylabel(r"$F$ (retrodiction fidelity)", fontsize=14)
    ax.set_title("Petz Recovery Bound Verification", fontsize=15)
    ax.legend(fontsize=8, loc="lower left")
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlim(left=-0.05)
    ax.grid(True, alpha=0.3)

    # Annotate bound status
    n_satisfied = sum(1 for r in results.values() if r["bound_satisfied_transformer"])
    n_total = len(results)
    status = "ALL PASS" if n_satisfied == n_total else f"{n_satisfied}/{n_total} PASS"
    ax.text(0.98, 0.98,
            f"Bound check: {status}\n"
            f"Filled = Transformer\n"
            f"Small = Linear / Hollow = AR",
            transform=ax.transAxes, fontsize=9,
            verticalalignment="top", horizontalalignment="right",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.8))

    # ---- Right panel: Fidelity comparison bar chart ----
    ax2 = axes[1]

    states_present = [s for s in state_order if s in results]
    n_states = len(states_present)
    x_pos = np.arange(n_states)

    has_fopt = any(results[s].get("f_optimal", -1) >= 0 for s in states_present)
    n_bars = 5 if has_fopt else 4
    width = 0.8 / n_bars

    fid_tf = [results[s]["fidelity_transformer"] for s in states_present]
    fid_lin = [results[s]["fidelity_linear"] for s in states_present]
    fid_ar = [results[s]["fidelity_ar"] for s in states_present]
    petz_b = [results[s]["petz_bound"] for s in states_present]

    offset = 0
    ax2.bar(x_pos + offset * width - (n_bars - 1) * width / 2, fid_tf, width,
            label="Transformer", color="#3498db", edgecolor="black")
    offset += 1
    ax2.bar(x_pos + offset * width - (n_bars - 1) * width / 2, fid_lin, width,
            label="Linear", color="#e67e22", edgecolor="black")
    offset += 1
    ax2.bar(x_pos + offset * width - (n_bars - 1) * width / 2, fid_ar, width,
            label="AR(5)", color="#2ecc71", edgecolor="black")
    offset += 1
    ax2.bar(x_pos + offset * width - (n_bars - 1) * width / 2, petz_b, width,
            label=r"$e^{-\Sigma/2}$ (bound)", color="#e74c3c", edgecolor="black", alpha=0.7)
    offset += 1

    if has_fopt:
        fid_opt = [max(0, results[s].get("f_optimal", 0)) for s in states_present]
        ax2.bar(x_pos + offset * width - (n_bars - 1) * width / 2, fid_opt, width,
                label=r"$F_{\mathrm{opt}}$ (analytical)", color="#9b59b6",
                edgecolor="black", alpha=0.7)

    ax2.set_xlabel("Consciousness State", fontsize=14)
    ax2.set_ylabel("Retrodiction Fidelity F", fontsize=14)
    ax2.set_title("Model Comparison by State", fontsize=15)
    ax2.set_xticks(x_pos)
    ax2.set_xticklabels(states_present, fontsize=12)
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.3, axis="y")
    ax2.set_ylim(0, 1.05)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    logger.info(f"Plot saved to: {output_path}")
    plt.close()


def plot_training_curves(
    all_histories: Dict[str, Dict[str, List[float]]],
    output_path: str = "training_curves.png",
) -> None:
    """Plot training loss curves for all states."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = {"W": "#e74c3c", "N1": "#e67e22", "N2": "#f1c40f",
              "N3": "#2ecc71", "R": "#3498db"}

    fig, ax = plt.subplots(figsize=(10, 6))

    for state, history in all_histories.items():
        c = colors.get(state, "gray")
        epochs = range(1, len(history["train_loss"]) + 1)
        ax.plot(epochs, history["train_loss"], "-", color=c, alpha=0.5,
                label=f"{state} train")
        ax.plot(epochs, history["test_loss"], "--", color=c, linewidth=2,
                label=f"{state} test")

    ax.set_xlabel("Epoch", fontsize=14)
    ax.set_ylabel("MSE Loss", fontsize=14)
    ax.set_title("Retrodiction Training Curves", fontsize=15)
    ax.legend(fontsize=9, ncol=2)
    ax.set_yscale("log")
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    logger.info(f"Training curves saved to: {output_path}")
    plt.close()


def plot_retrodiction_examples(
    model: PetzRecoveryTransformer,
    test_loader: DataLoader,
    state_name: str,
    output_path: str = "retrodiction_example.png",
    n_examples: int = 3,
) -> None:
    """Plot example retrodictions: actual past vs predicted past."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    model.eval()
    model = model.to(DEVICE)

    # Get a batch
    x_batch, y_batch = next(iter(test_loader))
    x_batch = x_batch[:n_examples].to(DEVICE)
    y_batch = y_batch[:n_examples]

    with torch.no_grad():
        pred_batch = model(x_batch).cpu()

    fig, axes = plt.subplots(n_examples, 2, figsize=(14, 4 * n_examples))
    if n_examples == 1:
        axes = axes.reshape(1, -1)

    for i in range(n_examples):
        # Channel 0
        ch = 0
        t = np.arange(y_batch.shape[1])

        ax_left = axes[i, 0]
        ax_left.plot(t, y_batch[i, :, ch].numpy(), "b-", alpha=0.7, label="Actual past")
        ax_left.plot(t, pred_batch[i, :, ch].numpy(), "r--", alpha=0.7, label="Predicted past")
        ax_left.set_title(f"{state_name} Example {i+1}: Channel {ch}", fontsize=12)
        ax_left.legend(fontsize=9)
        ax_left.set_xlabel("Time step")
        ax_left.set_ylabel("Amplitude (normalized)")
        ax_left.grid(True, alpha=0.3)

        # All channels heatmap
        ax_right = axes[i, 1]
        error = (pred_batch[i] - y_batch[i]).abs().numpy().T  # (C, T)
        im = ax_right.imshow(error, aspect="auto", cmap="hot", interpolation="nearest")
        ax_right.set_title(f"Prediction Error |pred - true|", fontsize=12)
        ax_right.set_xlabel("Time step")
        ax_right.set_ylabel("Channel")
        plt.colorbar(im, ax=ax_right)

    plt.suptitle(f"Retrodiction Examples: {state_name}", fontsize=14, y=1.01)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    logger.info(f"Example plot saved to: {output_path}")
    plt.close()


# ===========================================================================
# 9. SUMMARY TABLE
# ===========================================================================

def print_summary_table(results: Dict[str, Dict[str, float]]) -> None:
    """Print a formatted summary table."""
    has_fopt = any(r.get("f_optimal", -1) >= 0 for r in results.values())

    header = (
        f"{'State':<6} {'Sigma':>7} {'Petz':>7} "
        + (f"{'F_opt':>7} " if has_fopt else "")
        + f"{'F_TF':>7} {'F_Lin':>7} {'F_AR':>7} "
        f"{'Gap':>7} {'Bound?':>7}"
    )
    width = len(header)
    logger.info("\n" + "=" * width)
    logger.info("RETRODICTION FIDELITY SUMMARY")
    logger.info("=" * width)
    logger.info(header)
    logger.info("-" * width)

    order = ["W", "N1", "N2", "N3", "R"]
    for state in order:
        if state not in results:
            continue
        r = results[state]
        gap = r["fidelity_transformer"] - r["petz_bound"]
        satisfied = "YES" if r["bound_satisfied_transformer"] else "NO"
        fopt_str = f"{r['f_optimal']:>7.4f} " if has_fopt and r.get("f_optimal", -1) >= 0 else ""
        logger.info(
            f"{state:<6} {r['sigma']:>7.4f} {r['petz_bound']:>7.4f} "
            + fopt_str
            + f"{r['fidelity_transformer']:>7.4f} {r['fidelity_linear']:>7.4f} "
            f"{r['fidelity_ar']:>7.4f} {gap:>+7.4f} {satisfied:>7}"
        )

    logger.info("-" * width)

    # Overall summary
    n_satisfied = sum(
        1 for r in results.values() if r["bound_satisfied_transformer"]
    )
    n_total = len(results)
    logger.info(
        f"\nPetz bound satisfied: {n_satisfied}/{n_total} states "
        f"({'ALL PASS' if n_satisfied == n_total else 'SOME FAIL'})"
    )

    # Check ordering: higher Sigma -> lower F?
    sigmas = [(s, results[s]["sigma"]) for s in order if s in results]
    fids = [(s, results[s]["fidelity_transformer"]) for s in order if s in results]
    sigma_order = sorted(sigmas, key=lambda x: x[1])
    fid_order = sorted(fids, key=lambda x: -x[1])  # descending F

    if [s[0] for s in sigma_order] == [f[0] for f in fid_order]:
        logger.info("Ordering check: Sigma ordering MATCHES inverse F ordering")
    else:
        logger.info(
            f"Ordering check: Sigma order = {[s[0] for s in sigma_order]}, "
            f"F order (desc) = {[f[0] for f in fid_order]}"
        )

    # Nonlinearity advantage
    for state in order:
        if state not in results:
            continue
        r = results[state]
        tf_advantage = r["fidelity_transformer"] - r["fidelity_linear"]
        logger.info(
            f"  {state}: Transformer advantage over Linear = {tf_advantage:+.4f}"
        )


# ===========================================================================
# 10. MAIN
# ===========================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Neural Petz Recovery Map for EEG Retrodiction"
    )
    parser.add_argument(
        "--synthetic", action="store_true",
        help="Run with synthetic OU-process data (known Sigma)"
    )
    parser.add_argument(
        "--edf", type=str, default=None,
        help="Path to PSG .edf file"
    )
    parser.add_argument(
        "--hypno", type=str, default=None,
        help="Path to Hypnogram .edf file"
    )
    parser.add_argument(
        "--sigma-json", type=str, default=None,
        help="Path to pre-computed sigma results JSON (from compute_sigma.py)"
    )
    parser.add_argument(
        "--n-channels", type=int, default=8,
        help="Number of EEG channels (default: 8)"
    )
    parser.add_argument(
        "--n-samples", type=int, default=256,
        help="Samples per epoch/segment (default: 256 = 1 second at 256 Hz)"
    )
    parser.add_argument(
        "--n-epochs-train", type=int, default=80,
        help="Training epochs (default: 80)"
    )
    parser.add_argument(
        "--d-model", type=int, default=128,
        help="Transformer model dimension (default: 128)"
    )
    parser.add_argument(
        "--n-heads", type=int, default=4,
        help="Number of attention heads (default: 4)"
    )
    parser.add_argument(
        "--n-layers", type=int, default=2,
        help="Number of Transformer layers (default: 2)"
    )
    parser.add_argument(
        "--lr", type=float, default=1e-3,
        help="Learning rate (default: 1e-3)"
    )
    parser.add_argument(
        "--output-dir", type=str, default=None,
        help="Output directory for plots and results"
    )
    args = parser.parse_args()

    # Determine output directory
    script_dir = Path(__file__).resolve().parent
    if args.output_dir:
        output_dir = Path(args.output_dir)
    else:
        output_dir = script_dir.parent / "results"
    output_dir.mkdir(parents=True, exist_ok=True)

    # -----------------------------------------------------------------------
    # MODE 1: Synthetic data
    # -----------------------------------------------------------------------
    if args.synthetic or (args.edf is None):
        logger.info("=" * 60)
        logger.info("SYNTHETIC DATA MODE (Ornstein-Uhlenbeck with known Sigma)")
        logger.info("=" * 60)

        states = generate_synthetic_consciousness_states(
            n_channels=args.n_channels,
            n_samples=args.n_samples,
            n_segments_per_state=500,
            seed=42,
        )

        results = run_retrodiction_pipeline(
            states,
            n_epochs=args.n_epochs_train,
            lr=args.lr,
            d_model=args.d_model,
            n_heads=args.n_heads,
            n_layers=args.n_layers,
        )

        print_summary_table(results)

        # Plots
        plot_petz_bound(
            results,
            output_path=str(output_dir / "petz_bound_verification.png"),
            title="Petz Recovery Bound (Synthetic OU Data)",
        )

        # Save results
        results_path = str(output_dir / "retrodiction_results.json")
        with open(results_path, "w") as f:
            json.dump(results, f, indent=2, default=float)
        logger.info(f"Results saved to: {results_path}")

    # -----------------------------------------------------------------------
    # MODE 2: Real EEG
    # -----------------------------------------------------------------------
    else:
        logger.info("=" * 60)
        logger.info("REAL EEG MODE")
        logger.info("=" * 60)

        if not args.edf or not args.hypno:
            logger.error("Both --edf and --hypno required for real EEG mode")
            sys.exit(1)

        states = load_real_eeg_states(
            args.edf, args.hypno,
            n_samples=args.n_samples,
            max_channels=args.n_channels,
        )

        # Load Sigma values from compute_sigma.py results
        if args.sigma_json:
            with open(args.sigma_json) as f:
                sigma_data = json.load(f)
            stage_summary = sigma_data.get("stage_summary", {})
            for state_name in states:
                if state_name in stage_summary:
                    sigma_val = stage_summary[state_name].get("sigma_mean", 0.0)
                    states[state_name]["sigma_true"] = sigma_val
                    logger.info(
                        f"  Loaded Sigma({state_name}) = {sigma_val:.4f} "
                        f"from {args.sigma_json}"
                    )
        else:
            logger.warning(
                "No --sigma-json provided. Sigma values will be 0. "
                "Run compute_sigma.py first to get Sigma values."
            )

        results = run_retrodiction_pipeline(
            states,
            n_epochs=args.n_epochs_train,
            lr=args.lr,
            d_model=args.d_model,
            n_heads=args.n_heads,
            n_layers=args.n_layers,
        )

        print_summary_table(results)

        plot_petz_bound(
            results,
            output_path=str(output_dir / "petz_bound_verification_real.png"),
            title="Petz Recovery Bound (Real EEG)",
        )

        results_path = str(output_dir / "retrodiction_results_real.json")
        with open(results_path, "w") as f:
            json.dump(results, f, indent=2, default=float)
        logger.info(f"Results saved to: {results_path}")

    logger.info("\nDone.")


if __name__ == "__main__":
    main()
