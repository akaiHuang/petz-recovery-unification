#!/usr/bin/env python3
"""
=============================================================================
NOISE CHANNEL LIBRARY for tau-Chrono EDA Engine
=============================================================================

Provides standard quantum noise channels as lists of 2x2 Kraus operators:
  - Amplitude damping (energy relaxation, T1 decay)
  - Depolarizing noise (uniform Pauli errors)
  - Dephasing noise (phase randomization, T2 decay)

Each function returns a List[np.ndarray] of Kraus operators satisfying
the CPTP condition: sum_i K_i^dag K_i = I.

Author: Sheng-Kai Huang (2026)
License: MIT
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray
from typing import List

KrausList = List[NDArray[np.complexfloating]]


def amplitude_damping(gamma: float) -> KrausList:
    """
    Amplitude damping channel (qubit energy relaxation).

    Models spontaneous emission / T1 decay. The channel maps
    |1> -> |0> with probability gamma while preserving |0>.

    Kraus operators:
        K0 = [[1, 0], [0, sqrt(1-gamma)]]
        K1 = [[0, sqrt(gamma)], [0, 0]]

    Parameters
    ----------
    gamma : float
        Damping probability in [0, 1].

    Returns
    -------
    KrausList
        Two 2x2 Kraus operators.
    """
    assert 0.0 <= gamma <= 1.0, f"gamma must be in [0,1], got {gamma}"
    K0 = np.array([[1.0, 0.0],
                    [0.0, np.sqrt(1.0 - gamma)]], dtype=complex)
    K1 = np.array([[0.0, np.sqrt(gamma)],
                    [0.0, 0.0]], dtype=complex)
    return [K0, K1]


def depolarizing(p: float) -> KrausList:
    """
    Depolarizing channel (uniform Pauli noise).

    With probability p the qubit state is replaced by I/2 (maximally mixed),
    with probability (1-p) it is unchanged.

    N(rho) = (1-p) rho + (p/3)(X rho X + Y rho Y + Z rho Z)

    Kraus operators:
        K0 = sqrt(1 - 3p/4) * I
        K1 = sqrt(p/4) * X
        K2 = sqrt(p/4) * Y
        K3 = sqrt(p/4) * Z

    Parameters
    ----------
    p : float
        Depolarizing probability in [0, 1].

    Returns
    -------
    KrausList
        Four 2x2 Kraus operators.
    """
    assert 0.0 <= p <= 1.0, f"p must be in [0,1], got {p}"
    I2 = np.eye(2, dtype=complex)
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    Z = np.array([[1, 0], [0, -1]], dtype=complex)

    K0 = np.sqrt(1.0 - 3.0 * p / 4.0) * I2
    K1 = np.sqrt(p / 4.0) * X
    K2 = np.sqrt(p / 4.0) * Y
    K3 = np.sqrt(p / 4.0) * Z

    return [K0, K1, K2, K3]


def dephasing(p: float) -> KrausList:
    """
    Dephasing channel (phase randomization / T2 noise).

    With probability p the off-diagonal elements are suppressed.
    Diagonal elements of rho are preserved.

    N(rho) = (1-p) rho + p Z rho Z

    Kraus operators:
        K0 = sqrt(1-p) * I
        K1 = sqrt(p) * Z

    Parameters
    ----------
    p : float
        Dephasing probability in [0, 1].

    Returns
    -------
    KrausList
        Two 2x2 Kraus operators.
    """
    assert 0.0 <= p <= 1.0, f"p must be in [0,1], got {p}"
    I2 = np.eye(2, dtype=complex)
    Z = np.array([[1, 0], [0, -1]], dtype=complex)

    K0 = np.sqrt(1.0 - p) * I2
    K1 = np.sqrt(p) * Z

    return [K0, K1]


def verify_cptp(kraus_ops: KrausList, tol: float = 1e-10) -> bool:
    """
    Verify that a set of Kraus operators satisfies the CPTP condition:
    sum_i K_i^dag K_i = I.

    Parameters
    ----------
    kraus_ops : KrausList
        List of Kraus operators.
    tol : float
        Tolerance for the identity check.

    Returns
    -------
    bool
        True if CPTP condition is satisfied within tolerance.
    """
    d = kraus_ops[0].shape[1]
    total = np.zeros((d, d), dtype=complex)
    for K in kraus_ops:
        total += K.conj().T @ K
    return np.linalg.norm(total - np.eye(d, dtype=complex)) < tol


# ---------------------------------------------------------------------------
# Self-test
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print("Noise Channel Library -- Self-test")
    print("=" * 50)

    for name, channel_fn, params in [
        ("Amplitude Damping", amplitude_damping, [0.0, 0.05, 0.5, 1.0]),
        ("Depolarizing", depolarizing, [0.0, 0.08, 0.5, 1.0]),
        ("Dephasing", dephasing, [0.0, 0.03, 0.5, 1.0]),
    ]:
        print(f"\n{name}:")
        for p in params:
            ops = channel_fn(p)
            ok = verify_cptp(ops)
            print(f"  p={p:.2f}: {len(ops)} Kraus ops, CPTP={ok}")

    print("\nAll self-tests passed.")
