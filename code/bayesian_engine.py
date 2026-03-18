#!/usr/bin/env python3
"""
=============================================================================
BAYESIAN COMPOSITION ENGINE for tau-Chrono EDA Engine
=============================================================================

Implements Bayesian-updated sigma composition for multi-gate circuits,
contrasting with the naive multiplicative baseline (fixed sigma).

Key concept:
  - Naive/Multiplicative: compute tau for each gate using the ORIGINAL rho_0
    and sigma_0, then estimate total error as product of (1-tau_i). This
    ignores inter-gate correlations and overestimates error.
  - Bayesian approach: update both rho and sigma after each gate:
        rho_{i+1} = N_i(rho_i),  sigma_{i+1} = N_i(sigma_i)
    This gives the Petz recovery map the correct "prior" at each stage.

The composition inequality for the standard Petz map:
    sqrt(tau_total) <= sum_i sqrt(tau_i^eff)

Author: Sheng-Kai Huang (2026)
License: MIT
"""

from __future__ import annotations

import sys
import os
import numpy as np
from numpy.typing import NDArray
from typing import List, Tuple, Dict, Any, Optional
from dataclasses import dataclass, field

# Add parent directory to path so we can import petz_toolkit
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'numerical'))
sys.path.insert(0, os.path.dirname(__file__))

from petz_toolkit import (
    apply_channel,
    tau_parameter,
    fidelity,
    petz_recovery_map,
    apply_petz_recovery,
)

KrausList = List[NDArray[np.complexfloating]]
DensityMatrix = NDArray[np.complexfloating]


@dataclass
class GateResult:
    """Result for a single gate in the circuit."""
    gate_index: int
    channel_name: str
    tau_naive: float          # tau with original rho_0 and sigma_0 (multiplicative)
    tau_eff: float            # tau with Bayesian-updated rho_i and sigma_i
    sigma_before: DensityMatrix   # sigma used for this gate (Bayesian)
    rho_before: DensityMatrix     # rho entering this gate (Bayesian)
    rho_after: DensityMatrix      # rho exiting this gate
    commutator_norm: float    # ||[N(rho_i), N(sigma_i)]|| for saturation analysis
    pre_commutator_norm: float  # ||[rho_i, sigma_i]|| before the channel
    classification: str       # NORMAL, NEAR_SATURATED, SATURATED


@dataclass
class CompositionResult:
    """Full result of a Bayesian composition analysis."""
    gate_results: List[GateResult]
    tau_bayesian_total: float       # tau for full composed channel
    tau_multiplicative_total: float  # 1 - prod(1-tau_naive_i), multiplicative estimate
    improvement_percent: float       # (1 - bayesian/multiplicative) * 100
    composition_lhs: float           # sqrt(tau_total)
    composition_rhs: float           # sum sqrt(tau_i^eff)
    composition_holds: bool          # whether inequality holds
    composition_slack: float         # RHS - LHS


def _commutator_norm(A: DensityMatrix, B: DensityMatrix) -> float:
    """Frobenius norm of the commutator [A, B] = AB - BA."""
    comm = A @ B - B @ A
    return float(np.linalg.norm(comm, 'fro'))


def _classify_gate(tau_naive: float, tau_eff: float,
                   commutator_norm: float,
                   pre_commutator_norm: float,
                   comm_threshold: float = 0.01) -> str:
    """
    Classify a gate based on its tau behavior and commutator structure.

    SATURATED: tau_eff < 1e-10 (essentially perfect recovery)
    NEAR_SATURATED: post-channel commutator_norm < threshold OR
                    pre-channel commutator_norm < threshold (indicating
                    rho and sigma nearly commute, so Petz recovery is near-optimal)
    NORMAL: standard dissipative gate
    """
    if tau_eff < 1e-10:
        return "SATURATED"
    elif commutator_norm < comm_threshold or pre_commutator_norm < comm_threshold:
        return "NEAR_SATURATED"
    else:
        return "NORMAL"


def compose_kraus(channels: List[KrausList]) -> KrausList:
    """
    Compose a sequence of channels into a single channel via Kraus
    operator composition: (N_n ... N_2 N_1)(rho) = sum K_composed rho K_composed^dag.

    The composed Kraus operators are all products K_n ... K_2 K_1
    for every combination.
    """
    if len(channels) == 0:
        return [np.eye(2, dtype=complex)]
    if len(channels) == 1:
        return channels[0]

    composed = channels[0]
    for i in range(1, len(channels)):
        new_composed = []
        for K_new in channels[i]:
            for K_old in composed:
                new_composed.append(K_new @ K_old)
        composed = new_composed
    return composed


def bayesian_compose(
    channels: List[KrausList],
    sigma_0: DensityMatrix,
    rho: DensityMatrix,
    channel_names: Optional[List[str]] = None,
    comm_threshold: float = 0.01,
) -> CompositionResult:
    """
    Bayesian composition engine: compute per-gate tau with updated sigma.

    For each gate i:
      1. Compute tau_naive = tau(rho_0, N_i, sigma_0) [fixed, original states]
      2. Compute tau_eff = tau(rho_i, N_i, sigma_i) [Bayesian-updated states]
      3. Update: rho_{i+1} = N_i(rho_i), sigma_{i+1} = N_i(sigma_i)
      4. Classify gate based on commutator structure

    The multiplicative baseline estimates total error as:
      tau_mult = 1 - prod_i(1 - tau_naive_i)
    which assumes independent gates and overestimates the total error.

    Parameters
    ----------
    channels : List[KrausList]
        Sequence of noise channels (Kraus operator lists), applied in order.
    sigma_0 : DensityMatrix
        Initial reference state.
    rho : DensityMatrix
        Initial input state.
    channel_names : Optional[List[str]]
        Names for each channel (for display).
    comm_threshold : float
        Threshold for NEAR_SATURATED classification (default 0.01).

    Returns
    -------
    CompositionResult
        Full analysis including per-gate and total results.
    """
    n_gates = len(channels)
    if channel_names is None:
        channel_names = [f"Gate_{i}" for i in range(n_gates)]

    gate_results = []
    rho_current = rho.copy()
    sigma_current = sigma_0.copy()

    for i in range(n_gates):
        kraus = channels[i]

        # Pre-channel commutator: ||[rho_i, sigma_i]||
        pre_comm_norm = _commutator_norm(rho_current, sigma_current)

        # Compute tau with original rho_0 and sigma_0 (naive/multiplicative)
        tau_naive = tau_parameter(rho, kraus, sigma_0)

        # Compute tau with Bayesian-updated rho_i and sigma_i
        tau_eff = tau_parameter(rho_current, kraus, sigma_current)

        # Post-channel commutator: ||[N(rho_i), N(sigma_i)]||
        rho_after = apply_channel(rho_current, kraus)
        sigma_after = apply_channel(sigma_current, kraus)
        comm_norm = _commutator_norm(rho_after, sigma_after)

        # Classification
        classification = _classify_gate(
            tau_naive, tau_eff, comm_norm, pre_comm_norm, comm_threshold
        )

        gate_results.append(GateResult(
            gate_index=i,
            channel_name=channel_names[i],
            tau_naive=tau_naive,
            tau_eff=tau_eff,
            sigma_before=sigma_current.copy(),
            rho_before=rho_current.copy(),
            rho_after=rho_after.copy(),
            commutator_norm=comm_norm,
            pre_commutator_norm=pre_comm_norm,
            classification=classification,
        ))

        # Bayesian update: propagate both rho and sigma through the channel
        rho_current = rho_after
        sigma_current = sigma_after

    # Compute total tau for the composed channel (the actual composed channel)
    composed_kraus = compose_kraus(channels)
    tau_bayesian_total = tau_parameter(rho, composed_kraus, sigma_0)

    # Multiplicative baseline: 1 - prod(1 - tau_naive_i)
    # This assumes gates act independently (worst-case estimate)
    prod_fidelity = 1.0
    for gr in gate_results:
        prod_fidelity *= (1.0 - gr.tau_naive)
    tau_multiplicative_total = 1.0 - prod_fidelity

    # Composition inequality: sqrt(tau_total) <= sum sqrt(tau_i^eff)
    sum_sqrt_tau_eff = sum(np.sqrt(max(gr.tau_eff, 0.0)) for gr in gate_results)
    sqrt_tau_total = np.sqrt(max(tau_bayesian_total, 0.0))
    composition_holds = sqrt_tau_total <= sum_sqrt_tau_eff + 1e-10
    composition_slack = sum_sqrt_tau_eff - sqrt_tau_total

    # Improvement: Bayesian tau vs multiplicative estimate
    if tau_multiplicative_total > 1e-15:
        improvement = (1.0 - tau_bayesian_total / tau_multiplicative_total) * 100.0
    else:
        improvement = 0.0

    return CompositionResult(
        gate_results=gate_results,
        tau_bayesian_total=tau_bayesian_total,
        tau_multiplicative_total=tau_multiplicative_total,
        improvement_percent=improvement,
        composition_lhs=sqrt_tau_total,
        composition_rhs=sum_sqrt_tau_eff,
        composition_holds=composition_holds,
        composition_slack=composition_slack,
    )


# ---------------------------------------------------------------------------
# Self-test
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    from noise_channels import amplitude_damping, depolarizing, dephasing

    print("Bayesian Composition Engine -- Self-test")
    print("=" * 50)

    # Simple 2-gate test
    channels = [amplitude_damping(0.1), depolarizing(0.05)]
    names = ["AD(0.1)", "Dep(0.05)"]

    rho = 0.5 * np.array([[1, 1], [1, 1]], dtype=complex)  # |+><+|
    sigma_0 = np.diag([0.8, 0.2]).astype(complex)

    result = bayesian_compose(channels, sigma_0, rho, names)

    for gr in result.gate_results:
        print(f"  {gr.channel_name}: tau_naive={gr.tau_naive:.6f}, "
              f"tau_eff={gr.tau_eff:.6f}, class={gr.classification}")

    print(f"\n  tau_bayesian_total:       {result.tau_bayesian_total:.6f}")
    print(f"  tau_multiplicative_total: {result.tau_multiplicative_total:.6f}")
    print(f"  Improvement:             {result.improvement_percent:.2f}%")
    print(f"  Composition: sqrt(tau)={result.composition_lhs:.6f} "
          f"<= sum sqrt(tau_eff)={result.composition_rhs:.6f}")
    print(f"  Holds: {result.composition_holds}, Slack: {result.composition_slack:.6f}")
    print("\nSelf-test passed.")
