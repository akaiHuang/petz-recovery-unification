"""
adaptive_scheduler.py -- Sigma-aware adaptive scheduling for diffusion.

Uses the per-step Sigma profile to:
    1. Skip low-information steps  (compute_skip_mask)
    2. Allocate a fixed compute budget across the most informative steps
       (compute_adaptive_timesteps)
    3. Detect early-stop points where the Petz fidelity bound is already
       above a quality threshold  (compute_early_stop)

All functions are pure numpy -- no GPU dependency.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Union

import numpy as np


# ===================================================================
# Dataclass
# ===================================================================


@dataclass
class AdaptiveSchedule:
    """Result of an adaptive scheduling computation."""

    active_timesteps: np.ndarray  # int array of selected timestep indices
    step_sizes: np.ndarray  # gap between consecutive active timesteps
    speedup_ratio: float  # original_steps / active_steps
    quality_loss_estimate: float  # sum of skipped Sigma (lower is better)

    def __repr__(self) -> str:
        return (
            f"AdaptiveSchedule(active={len(self.active_timesteps)} steps, "
            f"speedup={self.speedup_ratio:.2f}x, "
            f"quality_loss_est={self.quality_loss_estimate:.4f})"
        )


# ===================================================================
# Core functions
# ===================================================================


def compute_skip_mask(
    sigma_array: Union[np.ndarray, list],
    threshold: float,
) -> np.ndarray:
    """Return a boolean mask where True = *keep* this step.

    Steps whose Sigma_t < ``threshold`` are candidates for skipping
    (mask = False).  The first and last steps are always kept.

    Parameters
    ----------
    sigma_array : 1-D array-like of per-step Sigma values.
    threshold : float
        Minimum Sigma to retain a step.

    Returns
    -------
    np.ndarray of bool, same length as *sigma_array*.
    """
    sigma = np.asarray(sigma_array, dtype=np.float64)
    mask = sigma >= threshold
    # Always keep first and last steps
    if len(mask) > 0:
        mask[0] = True
        mask[-1] = True
    return mask


def compute_adaptive_timesteps(
    sigma_array: Union[np.ndarray, list],
    budget: int,
) -> List[int]:
    """Select the *budget* most informative timestep indices.

    Strategy: rank steps by descending Sigma, pick the top ``budget``
    indices, then sort them back into temporal order.  The first and
    last steps are always included (they count toward the budget).

    Parameters
    ----------
    sigma_array : 1-D array-like of per-step Sigma values.
    budget : int
        Maximum number of steps to keep (must be >= 2 for non-trivial
        schedules).

    Returns
    -------
    List[int]
        Sorted list of selected timestep indices.
    """
    sigma = np.asarray(sigma_array, dtype=np.float64)
    T = len(sigma)

    if budget >= T:
        return list(range(T))

    budget = max(budget, 2)

    # Guarantee first and last
    guaranteed = {0, T - 1}
    remaining_budget = budget - len(guaranteed)

    if remaining_budget <= 0:
        return sorted(guaranteed)

    # Rank interior steps by Sigma (descending)
    interior = np.arange(1, T - 1)
    if len(interior) == 0:
        return sorted(guaranteed)

    order = np.argsort(-sigma[interior])
    selected = set(interior[order[:remaining_budget]])
    selected |= guaranteed

    return sorted(selected)


def compute_early_stop(
    sigma_array: Union[np.ndarray, list],
    fidelity_threshold: float,
) -> Optional[int]:
    """Find the earliest step at which the Petz fidelity bound exceeds
    *fidelity_threshold*.

    The Petz bound at step k is  F_k = exp(-cumulative_sigma_k / 2).
    We scan from the *last* step backward (reverse diffusion order),
    accumulating Sigma, until F_k drops below the threshold.

    Parameters
    ----------
    sigma_array : 1-D array of per-step Sigma values.
    fidelity_threshold : float in (0, 1]
        Desired minimum fidelity.

    Returns
    -------
    int or None
        The timestep index at which generation can stop and still
        satisfy the fidelity bound, or ``None`` if no early stop is
        possible (i.e. even all steps together do not meet the bound).
    """
    sigma = np.asarray(sigma_array, dtype=np.float64)
    T = len(sigma)

    if T == 0:
        return None

    # Cumulative Sigma from step 0 forward
    cumulative = np.cumsum(sigma)

    # Fidelity bound at each step
    fidelity = np.exp(-cumulative / 2.0)

    # Find the last step where fidelity is still >= threshold
    # (i.e. cumulative Sigma is still small enough)
    candidates = np.where(fidelity >= fidelity_threshold)[0]

    if len(candidates) == 0:
        return None

    # The latest step that still satisfies the bound -- we can stop here
    return int(candidates[-1])


# ===================================================================
# AdaptiveScheduler  -- convenience wrapper
# ===================================================================


class AdaptiveScheduler:
    """High-level scheduler that wraps the core functions.

    Example
    -------
    >>> sched = AdaptiveScheduler(budget=20, skip_threshold=0.01)
    >>> schedule = sched.schedule(sigma_profile)
    >>> print(schedule.speedup_ratio)
    """

    def __init__(
        self,
        budget: Optional[int] = None,
        skip_threshold: Optional[float] = None,
        fidelity_threshold: Optional[float] = None,
    ) -> None:
        self.budget = budget
        self.skip_threshold = skip_threshold
        self.fidelity_threshold = fidelity_threshold

    def schedule(
        self,
        sigma_array: Union[np.ndarray, list],
    ) -> AdaptiveSchedule:
        """Compute an adaptive schedule from a Sigma profile.

        The strategy is resolved in priority order:
        1. If ``budget`` is set, use top-k selection.
        2. Else if ``skip_threshold`` is set, use threshold-based masking.
        3. Else keep all steps (no speedup).

        If ``fidelity_threshold`` is also set, the schedule is further
        truncated by early stopping.

        Returns
        -------
        AdaptiveSchedule
        """
        sigma = np.asarray(sigma_array, dtype=np.float64)
        T = len(sigma)

        # -- Step selection ----------------------------------------------------
        if self.budget is not None:
            active = np.array(compute_adaptive_timesteps(sigma, self.budget), dtype=np.int64)
        elif self.skip_threshold is not None:
            mask = compute_skip_mask(sigma, self.skip_threshold)
            active = np.where(mask)[0].astype(np.int64)
        else:
            active = np.arange(T, dtype=np.int64)

        # -- Early stop --------------------------------------------------------
        if self.fidelity_threshold is not None:
            stop = compute_early_stop(sigma, self.fidelity_threshold)
            if stop is not None:
                active = active[active <= stop]
                # Ensure at least the first step
                if len(active) == 0:
                    active = np.array([0], dtype=np.int64)

        # -- Build schedule dataclass -----------------------------------------
        if len(active) < 2:
            step_sizes = np.array([1], dtype=np.int64)
        else:
            step_sizes = np.diff(active)

        skipped_mask = np.ones(T, dtype=bool)
        skipped_mask[active] = False
        quality_loss = float(sigma[skipped_mask].sum()) if skipped_mask.any() else 0.0

        speedup = T / max(len(active), 1)

        return AdaptiveSchedule(
            active_timesteps=active,
            step_sizes=step_sizes,
            speedup_ratio=speedup,
            quality_loss_estimate=quality_loss,
        )
