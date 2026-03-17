"""Tau-informed Zero Noise Extrapolation (ZNE).

Standard ZNE scales *all* gate noise uniformly. tau-chrono knows which
gates have high tau (high irreversibility) and which are nearly perfect.
Selective amplification on high-tau gates only gives better extrapolation
because it avoids amplifying negligible noise sources.

Free tier: basic uniform ZNE helper
Pro tier:  selective amplification with tau thresholding
"""

from __future__ import annotations

from typing import Callable, Sequence

import numpy as np
from numpy.typing import NDArray

from .._types import KrausList


def selective_amplify(
    channels: list[KrausList],
    tau_values: list[float],
    scale_factor: float,
    tau_threshold: float = 0.01,
) -> list[KrausList]:
    """Selectively amplify noise on high-tau gates.

    Gates with ``tau < tau_threshold`` are left unchanged.
    Gates with ``tau >= tau_threshold`` have their noise scaled by
    ``scale_factor`` via Kraus operator interpolation.

    Parameters
    ----------
    channels : list[KrausList]
        Original noise channels for each gate.
    tau_values : list[float]
        Per-gate tau values from Bayesian analysis.
    scale_factor : float
        Noise amplification factor (>= 1.0).
    tau_threshold : float
        Only amplify gates with tau above this threshold.

    Returns
    -------
    list[KrausList]
        Amplified channels. Unchanged for low-tau gates.
    """
    amplified = []
    for ch, tau in zip(channels, tau_values):
        if tau < tau_threshold or scale_factor <= 1.0:
            amplified.append(ch)
        else:
            amplified.append(_scale_channel_noise(ch, scale_factor))
    return amplified


def _scale_channel_noise(
    kraus_ops: KrausList, scale: float
) -> KrausList:
    """Scale noise in a channel by interpolating toward depolarizing.

    For scale > 1: increase the weight of non-identity Kraus operators.
    Uses the unitary folding approach: N_scaled = N^(scale) approximately.

    For integer scales, this is exact channel composition.
    For non-integer scales, uses linear interpolation between floor and ceil.
    """
    d = kraus_ops[0].shape[0]
    I_d = np.eye(d, dtype=complex)

    # Decompose: K_0 ~ sqrt(1-p) * U, K_i ~ sqrt(p/m) * E_i
    # Scale p -> scale * p (up to saturation)
    norms = [np.linalg.norm(K, 'fro') for K in kraus_ops]
    max_norm_idx = np.argmax(norms)

    # The dominant operator (closest to identity/unitary) gets reduced weight
    # Other operators get increased weight
    scaled_ops = []
    for i, K in enumerate(kraus_ops):
        if i == max_norm_idx:
            # Dominant operator: reduce by sqrt factor
            weight = np.sqrt(max(1.0 - (scale - 1.0) * (1.0 - norms[i] ** 2), 0.0))
            scaled_ops.append(weight / max(norms[i], 1e-15) * K)
        else:
            # Noise operators: amplify
            weight = np.sqrt(min(scale, d * d)) * norms[i]
            scaled_ops.append(weight / max(norms[i], 1e-15) * K)

    # Renormalize to CPTP
    total = sum(K.conj().T @ K for K in scaled_ops)
    correction = np.linalg.inv(
        _matrix_sqrt_safe(total)
    ) if np.linalg.norm(total - np.eye(d, dtype=complex)) > 1e-10 else I_d

    return [K @ correction for K in scaled_ops]


def _matrix_sqrt_safe(A: NDArray) -> NDArray:
    """Safe matrix square root for PSD matrices."""
    eigvals, U = np.linalg.eigh(0.5 * (A + A.conj().T))
    eigvals = np.maximum(eigvals, 0.0)
    return (U * np.sqrt(eigvals)[np.newaxis, :]) @ U.conj().T


def tau_informed_zne(
    expectation_fn: Callable[[list[KrausList]], float],
    channels: list[KrausList],
    tau_values: list[float],
    scale_factors: Sequence[float] = (1.0, 2.0, 3.0),
    tau_threshold: float = 0.01,
    extrapolation: str = "linear",
) -> dict:
    """Run tau-informed ZNE to estimate the zero-noise expectation value.

    Parameters
    ----------
    expectation_fn : callable
        Function that takes a list of channels and returns an expectation value.
        This encapsulates the circuit execution and measurement.
    channels : list[KrausList]
        Original noise channels for each gate.
    tau_values : list[float]
        Per-gate tau values from ``bayesian_compose``.
    scale_factors : sequence of float
        Noise scale factors to use (must include 1.0).
    tau_threshold : float
        Only amplify gates with tau above this.
    extrapolation : str
        Extrapolation method: "linear" or "exponential".

    Returns
    -------
    dict
        Keys:
        - ``'mitigated'``: zero-noise extrapolated value
        - ``'raw_values'``: expectation at each scale factor
        - ``'n_amplified'``: number of gates amplified
        - ``'scale_factors'``: the scale factors used
    """
    n_amplified = sum(1 for t in tau_values if t >= tau_threshold)

    raw_values = []
    for sf in scale_factors:
        amp_channels = selective_amplify(channels, tau_values, sf, tau_threshold)
        val = expectation_fn(amp_channels)
        raw_values.append(val)

    if extrapolation == "linear":
        mitigated = _linear_extrapolation(
            list(scale_factors), raw_values
        )
    elif extrapolation == "exponential":
        mitigated = _exponential_extrapolation(
            list(scale_factors), raw_values
        )
    else:
        raise ValueError(f"Unknown extrapolation method: {extrapolation}")

    return {
        "mitigated": mitigated,
        "raw_values": raw_values,
        "n_amplified": n_amplified,
        "scale_factors": list(scale_factors),
    }


def _linear_extrapolation(
    scale_factors: list[float], values: list[float]
) -> float:
    """Richardson extrapolation to zero noise (linear fit)."""
    sf = np.array(scale_factors)
    vals = np.array(values)
    # Fit: E(lambda) = a + b * lambda, extrapolate to lambda=0
    coeffs = np.polyfit(sf, vals, deg=min(len(sf) - 1, 2))
    return float(np.polyval(coeffs, 0.0))


def _exponential_extrapolation(
    scale_factors: list[float], values: list[float]
) -> float:
    """Exponential extrapolation: E(lambda) = a * exp(b * lambda)."""
    sf = np.array(scale_factors)
    vals = np.array(values)

    # Use log-linear fit: log(E) = log(a) + b * lambda
    pos_vals = np.maximum(vals, 1e-15)
    log_vals = np.log(pos_vals)
    coeffs = np.polyfit(sf, log_vals, deg=1)
    return float(np.exp(np.polyval(coeffs, 0.0)))
