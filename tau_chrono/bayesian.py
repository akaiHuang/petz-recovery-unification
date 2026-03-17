"""
Bayesian composition engine for multi-gate noise tracking.

Implements Bayesian-updated sigma composition for quantum circuits.
The key insight: updating both rho and sigma after each gate gives the
Petz recovery map the correct "prior" at each stage, yielding per-gate
tau estimates that are 50-100% more accurate than independent noise models.

Public API:
    - ``bayesian_compose``: run full Bayesian composition analysis
    - ``compose_kraus``: compose a sequence of Kraus channels
    - ``GateResult``: per-gate result dataclass
    - ``CompositionResult``: full circuit result dataclass
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray
from typing import List, Optional
from dataclasses import dataclass

from .petz import apply_channel, tau_parameter
from .utils import commutator_norm

KrausList = List[NDArray[np.complexfloating]]
DensityMatrix = NDArray[np.complexfloating]


@dataclass
class GateResult:
    """Result for a single gate in the circuit."""

    gate_index: int
    channel_name: str
    tau_naive: float
    tau_eff: float
    sigma_before: DensityMatrix
    rho_before: DensityMatrix
    rho_after: DensityMatrix
    commutator_norm: float
    pre_commutator_norm: float
    classification: str  # NORMAL, NEAR_SATURATED, SATURATED


@dataclass
class CompositionResult:
    """Full result of a Bayesian composition analysis."""

    gate_results: List[GateResult]
    tau_bayesian_total: float
    tau_multiplicative_total: float
    improvement_percent: float
    composition_lhs: float   # sqrt(tau_total)
    composition_rhs: float   # sum sqrt(tau_i^eff)
    composition_holds: bool
    composition_slack: float  # RHS - LHS


def _classify_gate(
    tau_naive: float,
    tau_eff: float,
    comm_norm: float,
    pre_comm_norm: float,
    comm_threshold: float = 0.01,
) -> str:
    """Classify a gate based on tau behaviour and commutator structure."""
    if tau_eff < 1e-10:
        return "SATURATED"
    elif comm_norm < comm_threshold or pre_comm_norm < comm_threshold:
        return "NEAR_SATURATED"
    else:
        return "NORMAL"


def compose_kraus(channels: List[KrausList]) -> KrausList:
    """Compose a sequence of Kraus channels into a single channel.

    The composed Kraus operators are all products ``K_n ... K_2 K_1``
    for every combination of per-channel operators.

    Parameters
    ----------
    channels : List[KrausList]
        Ordered list of channels to compose (first applied first).

    Returns
    -------
    KrausList
        Composed Kraus operators.
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


def compose_kraus_compressed(
    channels: List[KrausList],
    max_ops: int = 64,
    compress_threshold: int = 256,
) -> KrausList:
    """Compose channels with SVD compression to limit Kraus count.

    Useful for deep 2-qubit circuits where the naive product generates
    thousands of operators.

    Parameters
    ----------
    channels : List[KrausList]
        Ordered list of channels.
    max_ops : int
        Maximum number of Kraus operators after compression.
    compress_threshold : int
        Trigger compression when count exceeds this.

    Returns
    -------
    KrausList
        Compressed composed Kraus operators.
    """
    if len(channels) == 0:
        d = 2
        return [np.eye(d, dtype=complex)]
    if len(channels) == 1:
        return channels[0]

    d = channels[0][0].shape[0]
    composed = channels[0]
    for i in range(1, len(channels)):
        new_composed = []
        for K_new in channels[i]:
            for K_old in composed:
                prod = K_new @ K_old
                if np.linalg.norm(prod) > 1e-14:
                    new_composed.append(prod)
        composed = new_composed
        if len(composed) > compress_threshold:
            composed = _compress_kraus(composed, d, max_ops)
    return composed


def _compress_kraus(
    kraus_ops: KrausList, d: int, max_ops: int = 64
) -> KrausList:
    """Compress Kraus representation using SVD."""
    n = len(kraus_ops)
    M = np.zeros((n, d * d), dtype=complex)
    for i, K in enumerate(kraus_ops):
        M[i] = K.flatten()
    U, S, Vh = np.linalg.svd(M, full_matrices=False)
    threshold = 1e-12 * S[0] if len(S) > 0 else 0
    keep = min(max_ops, int(np.sum(S > threshold)))
    keep = max(keep, 1)
    compressed: KrausList = []
    for i in range(keep):
        K_new = (S[i] * Vh[i]).reshape(d, d)
        compressed.append(K_new)
    return compressed


def bayesian_compose(
    channels: List[KrausList],
    sigma_0: DensityMatrix,
    rho: DensityMatrix,
    channel_names: Optional[List[str]] = None,
    comm_threshold: float = 0.01,
    memory_alpha: Optional[float] = None,
) -> CompositionResult:
    """Bayesian composition engine: per-gate tau with updated sigma.

    For each gate *i*:

    1. ``tau_naive = tau(rho_0, N_i, sigma_0)`` -- fixed original states
    2. ``tau_eff = tau(rho_i, N_i, sigma_i)`` -- Bayesian-updated states
    3. Update: ``rho_{i+1} = N_i(rho_i)``, ``sigma_{i+1} = N_i(sigma_i)``
    4. Classify gate from commutator structure

    The multiplicative baseline ``1 - prod(1 - tau_naive_i)`` assumes
    independent gates and overestimates total error.

    Parameters
    ----------
    channels : List[KrausList]
        Sequence of noise channels (applied in order).
    sigma_0 : DensityMatrix
        Initial reference state.
    rho : DensityMatrix
        Initial input state.
    channel_names : list[str], optional
        Human-readable names for each channel.
    comm_threshold : float
        Threshold for NEAR_SATURATED classification.
    memory_alpha : float, optional
        Non-Markovian memory parameter in ``[0, 1]``.  When provided,
        sigma retains a fraction of its previous value at each step:
        ``sigma_{i+1} = alpha * N(sigma_i) + (1-alpha) * sigma_i``.
        *  ``alpha = 1.0`` recovers fully Markovian behaviour (no memory).
        *  ``alpha = 0.0`` freezes sigma (full memory, no update).
        Typical useful range is 0.85 -- 0.95.  Default ``None`` is
        equivalent to the Markovian case (``alpha = 1``).

    Returns
    -------
    CompositionResult
        Full analysis including per-gate and total results.
    """
    n_gates = len(channels)
    if channel_names is None:
        channel_names = [f"Gate_{i}" for i in range(n_gates)]

    gate_results: List[GateResult] = []
    rho_current = rho.copy()
    sigma_current = sigma_0.copy()

    for i in range(n_gates):
        kraus = channels[i]

        pre_comm_norm = commutator_norm(rho_current, sigma_current)
        tau_naive = tau_parameter(rho, kraus, sigma_0)
        tau_eff = tau_parameter(rho_current, kraus, sigma_current)

        rho_after = apply_channel(rho_current, kraus)
        sigma_after = apply_channel(sigma_current, kraus)

        # Non-Markovian memory correction
        if memory_alpha is not None and i > 0:
            sigma_after = memory_alpha * sigma_after + (1.0 - memory_alpha) * sigma_current
            # Ensure valid density matrix
            sigma_after = sigma_after / np.trace(sigma_after).real

        comm_norm = commutator_norm(rho_after, sigma_after)

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

        rho_current = rho_after
        sigma_current = sigma_after

    # Total tau for the composed channel (use compression for deep circuits)
    if len(channels) > 8:
        composed_kraus = compose_kraus_compressed(channels)
    else:
        composed_kraus = compose_kraus(channels)
    tau_bayesian_total = tau_parameter(rho, composed_kraus, sigma_0)

    # Multiplicative baseline
    prod_fidelity = 1.0
    for gr in gate_results:
        prod_fidelity *= (1.0 - gr.tau_naive)
    tau_multiplicative_total = 1.0 - prod_fidelity

    # Composition inequality
    sum_sqrt_tau_eff = sum(np.sqrt(max(gr.tau_eff, 0.0)) for gr in gate_results)
    sqrt_tau_total = np.sqrt(max(tau_bayesian_total, 0.0))
    composition_holds = sqrt_tau_total <= sum_sqrt_tau_eff + 1e-10
    composition_slack = sum_sqrt_tau_eff - sqrt_tau_total

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
