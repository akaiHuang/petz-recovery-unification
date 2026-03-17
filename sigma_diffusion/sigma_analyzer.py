"""
sigma_analyzer.py -- Core Sigma computation for diffusion processes.

Computes the information-theoretic divergence Sigma_t = D_KL(p_{t-1} || p_t)
between consecutive diffusion steps, for both continuous (Gaussian) and
discrete (token/logit) latent representations.

Tracks the Petz fidelity bound  F >= exp(-Sigma / 2)  at every step.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import List, Optional, Union

import numpy as np

# ---------------------------------------------------------------------------
# Lazy PyTorch import (allows numpy-only usage when torch is unavailable)
# ---------------------------------------------------------------------------
_torch = None


def _ensure_torch():
    global _torch
    if _torch is None:
        import torch

        _torch = torch
    return _torch


# ===================================================================
# Closed-form KL divergences
# ===================================================================


def kl_divergence_gaussian(
    mu1: Union[float, np.ndarray],
    var1: Union[float, np.ndarray],
    mu2: Union[float, np.ndarray],
    var2: Union[float, np.ndarray],
) -> float:
    """KL(N(mu1, diag(var1)) || N(mu2, diag(var2))) for diagonal Gaussians.

    All inputs may be scalars (single dimension) or 1-D arrays /
    tensors (independent dimensions whose contributions are summed).

    Returns
    -------
    float
        The KL divergence (>= 0).
    """
    mu1, var1, mu2, var2 = (np.asarray(x, dtype=np.float64) for x in (mu1, var1, mu2, var2))

    # Sanitise variances to avoid log(0)
    eps = 1e-12
    var1 = np.maximum(var1, eps)
    var2 = np.maximum(var2, eps)

    # D_KL = 0.5 * sum[ log(var2/var1) + (var1 + (mu1-mu2)^2) / var2 - 1 ]
    kl = 0.5 * np.sum(np.log(var2 / var1) + (var1 + (mu1 - mu2) ** 2) / var2 - 1.0)
    return float(max(kl, 0.0))


def kl_divergence_discrete(
    p: Union[np.ndarray, "torch.Tensor"],
    q: Union[np.ndarray, "torch.Tensor"],
) -> float:
    """KL(p || q) for discrete distributions (e.g. token logits after softmax).

    Parameters
    ----------
    p, q : array-like of shape (..., V)
        Probability vectors.  They need not be exactly normalised; each
        will be re-normalised along the last axis.

    Returns
    -------
    float
    """
    p = np.asarray(p, dtype=np.float64)
    q = np.asarray(q, dtype=np.float64)

    eps = 1e-12
    p = p / (p.sum(axis=-1, keepdims=True) + eps)
    q = q / (q.sum(axis=-1, keepdims=True) + eps)
    p = np.clip(p, eps, None)
    q = np.clip(q, eps, None)

    kl = np.sum(p * np.log(p / q))
    return float(max(kl, 0.0))


# ===================================================================
# Dataclasses
# ===================================================================


@dataclass
class StepDiagnostics:
    """Diagnostics for a single diffusion step."""

    t: int
    sigma_t: float
    cumulative_sigma: float
    petz_fidelity_bound: float  # exp(-cumulative_sigma / 2)
    snr_t: Optional[float] = None  # signal-to-noise ratio at this step


@dataclass
class SigmaProfile:
    """Full Sigma profile across all recorded diffusion steps."""

    steps: List[StepDiagnostics] = field(default_factory=list)
    total_sigma: float = 0.0

    # -- convenience accessors ------------------------------------------------

    def sigma_array(self) -> np.ndarray:
        """Per-step Sigma values as a 1-D numpy array."""
        return np.array([s.sigma_t for s in self.steps], dtype=np.float64)

    def cumulative_array(self) -> np.ndarray:
        """Cumulative Sigma values as a 1-D numpy array."""
        return np.array([s.cumulative_sigma for s in self.steps], dtype=np.float64)

    def fidelity_array(self) -> np.ndarray:
        """Petz fidelity bound at each step."""
        return np.array([s.petz_fidelity_bound for s in self.steps], dtype=np.float64)

    def timestep_array(self) -> np.ndarray:
        """Timestep indices as a 1-D int array."""
        return np.array([s.t for s in self.steps], dtype=np.int64)

    def __len__(self) -> int:
        return len(self.steps)


@dataclass
class SigmaResult:
    """Container returned by SigmaAnalyzer.analyze().

    Bundles the full profile together with summary statistics.
    """

    profile: SigmaProfile
    total_sigma: float
    mean_sigma: float
    max_sigma: float
    petz_fidelity_bound: float  # final cumulative bound
    num_steps: int


# ===================================================================
# SigmaTracker  -- online step-by-step recorder
# ===================================================================


class SigmaTracker:
    """Records per-step distributions and computes Sigma on the fly.

    Supports two modes, selected automatically per call:

    * **continuous** -- latent tensors are modelled as diagonal Gaussians
      (mean / variance computed per channel).
    * **discrete** -- probability vectors (e.g. post-softmax logits).
    """

    def __init__(self) -> None:
        self._steps: List[StepDiagnostics] = []
        self._cumulative: float = 0.0
        self._prev_stats: Optional[dict] = None
        self._step_counter: int = 0

    # -----------------------------------------------------------------
    # Public API
    # -----------------------------------------------------------------

    def record_step(
        self,
        latent: Union[np.ndarray, "torch.Tensor"],
        *,
        mode: str = "continuous",
        snr: Optional[float] = None,
        timestep: Optional[int] = None,
    ) -> StepDiagnostics:
        """Record a single diffusion step.

        Parameters
        ----------
        latent : array-like
            For ``mode="continuous"``: the latent tensor of shape
            ``(B, C, H, W)`` or ``(C, H, W)`` -- statistics are computed
            per channel.
            For ``mode="discrete"``: a probability vector of shape
            ``(..., V)``.
        mode : str
            ``"continuous"`` (default) or ``"discrete"``.
        snr : float, optional
            If available externally, the signal-to-noise ratio.
        timestep : int, optional
            Override the timestep index (otherwise auto-incremented).

        Returns
        -------
        StepDiagnostics
        """
        t = timestep if timestep is not None else self._step_counter

        arr = self._to_numpy(latent)

        if mode == "continuous":
            stats = self._gaussian_stats(arr)
            sigma_t = self._sigma_continuous(stats)
        elif mode == "discrete":
            stats = {"probs": arr}
            sigma_t = self._sigma_discrete(stats)
        else:
            raise ValueError(f"Unknown mode {mode!r}; use 'continuous' or 'discrete'.")

        self._prev_stats = stats
        self._cumulative += sigma_t
        self._step_counter += 1

        diag = StepDiagnostics(
            t=t,
            sigma_t=sigma_t,
            cumulative_sigma=self._cumulative,
            petz_fidelity_bound=math.exp(-self._cumulative / 2.0),
            snr_t=snr,
        )
        self._steps.append(diag)
        return diag

    def get_profile(self) -> SigmaProfile:
        """Return the full SigmaProfile accumulated so far."""
        return SigmaProfile(
            steps=list(self._steps),
            total_sigma=self._cumulative,
        )

    def reset(self) -> None:
        """Clear all recorded data."""
        self._steps.clear()
        self._cumulative = 0.0
        self._prev_stats = None
        self._step_counter = 0

    # -----------------------------------------------------------------
    # Internal helpers
    # -----------------------------------------------------------------

    @staticmethod
    def _to_numpy(x: Union[np.ndarray, "torch.Tensor"]) -> np.ndarray:
        """Convert to numpy, handling torch tensors on any device."""
        if isinstance(x, np.ndarray):
            return x.astype(np.float64)
        # Assume torch tensor
        return x.detach().float().cpu().numpy().astype(np.float64)

    @staticmethod
    def _gaussian_stats(arr: np.ndarray) -> dict:
        """Compute per-channel mean and variance.

        Accepts shapes (B, C, H, W) or (C, H, W).  Reduces over all
        axes *except* the channel axis (index -3 or 0 for 3-D).
        """
        if arr.ndim == 4:
            # (B, C, H, W) -> reduce over B, H, W
            axes = (0, 2, 3)
        elif arr.ndim == 3:
            # (C, H, W) -> reduce over H, W
            axes = (1, 2)
        elif arr.ndim == 2:
            # (C, L) -> reduce over L
            axes = (1,)
        elif arr.ndim == 1:
            axes = (0,)
        else:
            axes = tuple(range(arr.ndim))

        mu = arr.mean(axis=axes)
        var = arr.var(axis=axes) + 1e-12
        return {"mu": mu, "var": var}

    def _sigma_continuous(self, stats: dict) -> float:
        """KL between the current and previous Gaussian statistics."""
        if self._prev_stats is None or "mu" not in self._prev_stats:
            return 0.0
        return kl_divergence_gaussian(
            self._prev_stats["mu"],
            self._prev_stats["var"],
            stats["mu"],
            stats["var"],
        )

    def _sigma_discrete(self, stats: dict) -> float:
        """KL between the current and previous discrete distributions."""
        if self._prev_stats is None or "probs" not in self._prev_stats:
            return 0.0
        return kl_divergence_discrete(self._prev_stats["probs"], stats["probs"])


# ===================================================================
# SigmaAnalyzer  -- batch analysis over a full trajectory
# ===================================================================


class SigmaAnalyzer:
    """Analyse a *complete* diffusion trajectory (list of latents).

    Example
    -------
    >>> analyser = SigmaAnalyzer()
    >>> result = analyser.analyze(latents, mode="continuous")
    >>> print(result.total_sigma, result.petz_fidelity_bound)
    """

    def analyze(
        self,
        latents: list,
        *,
        mode: str = "continuous",
        snr_schedule: Optional[List[float]] = None,
    ) -> SigmaResult:
        """Run full Sigma analysis over a trajectory.

        Parameters
        ----------
        latents : list of array-like
            One entry per diffusion step (T entries total).
        mode : ``"continuous"`` or ``"discrete"``.
        snr_schedule : list of float, optional
            Per-step SNR values, if available.

        Returns
        -------
        SigmaResult
        """
        tracker = SigmaTracker()
        for i, lat in enumerate(latents):
            snr = snr_schedule[i] if snr_schedule is not None else None
            tracker.record_step(lat, mode=mode, snr=snr, timestep=i)

        profile = tracker.get_profile()
        sigma_arr = profile.sigma_array()

        return SigmaResult(
            profile=profile,
            total_sigma=profile.total_sigma,
            mean_sigma=float(sigma_arr.mean()) if len(sigma_arr) > 0 else 0.0,
            max_sigma=float(sigma_arr.max()) if len(sigma_arr) > 0 else 0.0,
            petz_fidelity_bound=math.exp(-profile.total_sigma / 2.0),
            num_steps=len(profile),
        )
