"""License management for tau-chrono Pro features.

Free tier: all core functions, up to 20 gates in bayesian_compose,
           1-qubit tomography, basic channels and metrics.

Pro tier:  unlimited gates, 2-qubit tomography, non-Markovian correction,
           error mitigation, noise drift tracking, priority support.
"""

from __future__ import annotations

import os
from pathlib import Path

_PRO_FEATURES = {
    "unlimited_gates",
    "2q_tomography",
    "non_markovian",
    "error_mitigation",
    "noise_drift",
    "advanced_viz",
}

FREE_GATE_LIMIT = 20

_cached_license: bool | None = None


def _read_license_key() -> str | None:
    """Read license key from env var or file."""
    key = os.environ.get("TAU_CHRONO_LICENSE")
    if key:
        return key.strip()

    key_path = Path.home() / ".tau-chrono" / "license.key"
    if key_path.exists():
        return key_path.read_text().strip()

    return None


def _validate_key(key: str) -> bool:
    """Validate a license key.

    For now, any non-empty key starting with 'TC-PRO-' is valid.
    In production, this would verify a cryptographic signature.
    """
    return key.startswith("TC-PRO-") and len(key) > 10


def is_pro() -> bool:
    """Check if Pro license is active."""
    global _cached_license
    if _cached_license is not None:
        return _cached_license

    key = _read_license_key()
    _cached_license = key is not None and _validate_key(key)
    return _cached_license


def require_pro(feature: str) -> None:
    """Raise if a Pro feature is used without a license.

    Parameters
    ----------
    feature : str
        Feature name for the error message.

    Raises
    ------
    RuntimeError
        If Pro license is not active.
    """
    if not is_pro():
        raise RuntimeError(
            f"tau-chrono Pro required for '{feature}'. "
            f"Get a license at https://tau-chrono.dev/pro or set "
            f"TAU_CHRONO_LICENSE env var."
        )


def check_gate_limit(n_gates: int) -> None:
    """Check if the number of gates exceeds the free tier limit.

    Parameters
    ----------
    n_gates : int
        Number of gates in the circuit.

    Raises
    ------
    RuntimeError
        If n_gates > FREE_GATE_LIMIT and no Pro license.
    """
    if n_gates > FREE_GATE_LIMIT and not is_pro():
        raise RuntimeError(
            f"Free tier supports up to {FREE_GATE_LIMIT} gates "
            f"(got {n_gates}). Upgrade to tau-chrono Pro for "
            f"unlimited circuit depth: https://tau-chrono.dev/pro"
        )
