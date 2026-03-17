"""
Quantum noise channel library.

Provides standard 1-qubit and 2-qubit noise channels as lists of Kraus
operators satisfying the CPTP condition ``sum_i K_i^dag K_i = I``.

1-qubit channels:
    - amplitude_damping (T1 decay)
    - depolarizing (uniform Pauli errors)
    - dephasing (T2 / phase randomization)

2-qubit channels:
    - two_qubit_depolarizing
    - correlated_dephasing (ZZ noise)
    - local_noise (independent per-qubit)
    - cnot_error (ideal CNOT + depolarizing)
    - crosstalk_dephasing
    - amplitude_damping_2q
    - unitary_channel (wrap any unitary as Kraus)
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray
from typing import List

KrausList = List[NDArray[np.complexfloating]]

# -- Pauli and gate constants -----------------------------------------------

I2 = np.eye(2, dtype=complex)
X_GATE = np.array([[0, 1], [1, 0]], dtype=complex)
Y_GATE = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z_GATE = np.array([[1, 0], [0, -1]], dtype=complex)
I4 = np.eye(4, dtype=complex)

CNOT = np.array([
    [1, 0, 0, 0],
    [0, 1, 0, 0],
    [0, 0, 0, 1],
    [0, 0, 1, 0],
], dtype=complex)

CNOT_REV = np.array([
    [1, 0, 0, 0],
    [0, 0, 0, 1],
    [0, 0, 1, 0],
    [0, 1, 0, 0],
], dtype=complex)


# ---------------------------------------------------------------------------
# 1-qubit channels
# ---------------------------------------------------------------------------


def amplitude_damping(gamma: float) -> KrausList:
    """Amplitude damping channel (qubit energy relaxation / T1 decay).

    Maps ``|1>`` to ``|0>`` with probability *gamma* while preserving ``|0>``.

    Parameters
    ----------
    gamma : float
        Damping probability in ``[0, 1]``.

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
    """Depolarizing channel (uniform Pauli noise).

    ``N(rho) = (1-p) rho + (p/3)(X rho X + Y rho Y + Z rho Z)``

    Parameters
    ----------
    p : float
        Depolarizing probability in ``[0, 1]``.

    Returns
    -------
    KrausList
        Four 2x2 Kraus operators.
    """
    assert 0.0 <= p <= 1.0, f"p must be in [0,1], got {p}"
    K0 = np.sqrt(1.0 - 3.0 * p / 4.0) * I2
    K1 = np.sqrt(p / 4.0) * X_GATE
    K2 = np.sqrt(p / 4.0) * Y_GATE
    K3 = np.sqrt(p / 4.0) * Z_GATE
    return [K0, K1, K2, K3]


def dephasing(p: float) -> KrausList:
    """Dephasing channel (phase randomization / T2 noise).

    ``N(rho) = (1-p) rho + p Z rho Z``

    Parameters
    ----------
    p : float
        Dephasing probability in ``[0, 1]``.

    Returns
    -------
    KrausList
        Two 2x2 Kraus operators.
    """
    assert 0.0 <= p <= 1.0, f"p must be in [0,1], got {p}"
    K0 = np.sqrt(1.0 - p) * I2
    K1 = np.sqrt(p) * Z_GATE
    return [K0, K1]


def verify_cptp(kraus_ops: KrausList, tol: float = 1e-10) -> bool:
    """Verify that Kraus operators satisfy the CPTP condition ``sum K_i^dag K_i = I``.

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
    return bool(np.linalg.norm(total - np.eye(d, dtype=complex)) < tol)


# ---------------------------------------------------------------------------
# 2-qubit channels
# ---------------------------------------------------------------------------


def two_qubit_depolarizing(p: float) -> KrausList:
    """Two-qubit depolarizing channel: ``N(rho) = (1-p) rho + p I/4``.

    Parameters
    ----------
    p : float
        Depolarizing probability in ``[0, 1]``.

    Returns
    -------
    KrausList
        16 Kraus operators of shape 4x4.
    """
    assert 0.0 <= p <= 1.0
    paulis_1q = [I2, X_GATE, Y_GATE, Z_GATE]
    kraus: KrausList = []
    kraus.append(np.sqrt(1.0 - 15.0 * p / 16.0) * I4)
    for i, Pi in enumerate(paulis_1q):
        for j, Pj in enumerate(paulis_1q):
            if i == 0 and j == 0:
                continue
            kraus.append(np.sqrt(p / 16.0) * np.kron(Pi, Pj))
    return kraus


def correlated_dephasing(p: float) -> KrausList:
    """Correlated dephasing: ``ZZ`` noise.

    ``N(rho) = (1-p) rho + p (Z x Z) rho (Z x Z)``

    Parameters
    ----------
    p : float
        Dephasing probability in ``[0, 1]``.

    Returns
    -------
    KrausList
        Two 4x4 Kraus operators.
    """
    assert 0.0 <= p <= 1.0
    ZZ = np.kron(Z_GATE, Z_GATE)
    K0 = np.sqrt(1.0 - p) * I4
    K1 = np.sqrt(p) * ZZ
    return [K0, K1]


def local_noise(kraus_q0: KrausList, kraus_q1: KrausList) -> KrausList:
    """Independent noise on each qubit: ``N_0 x N_1``.

    Parameters
    ----------
    kraus_q0 : KrausList
        Kraus operators for qubit 0.
    kraus_q1 : KrausList
        Kraus operators for qubit 1.

    Returns
    -------
    KrausList
        Tensor-product Kraus operators of shape 4x4.
    """
    kraus: KrausList = []
    for K0 in kraus_q0:
        for K1 in kraus_q1:
            kraus.append(np.kron(K0, K1))
    return kraus


def cnot_error(p: float) -> KrausList:
    """Imperfect CNOT: ideal CNOT followed by two-qubit depolarizing noise.

    Parameters
    ----------
    p : float
        Depolarizing probability applied after the ideal CNOT.

    Returns
    -------
    KrausList
        Kraus operators incorporating both the CNOT and noise.
    """
    dep_kraus = two_qubit_depolarizing(p)
    return [K @ CNOT for K in dep_kraus]


def crosstalk_dephasing(p: float, epsilon: float) -> KrausList:
    """Crosstalk noise: Z-dephasing on qubit 0 with crosstalk leaking to qubit 1.

    Parameters
    ----------
    p : float
        Dephasing probability in ``[0, 1]``.
    epsilon : float
        Crosstalk coupling strength.

    Returns
    -------
    KrausList
        Two 4x4 Kraus operators.
    """
    assert 0.0 <= p <= 1.0
    eigvals = np.array([1 + epsilon, 1 - epsilon, -1 + epsilon, -1 - epsilon])
    V = np.diag(np.exp(1j * np.pi / 2 * eigvals))
    K0 = np.sqrt(1.0 - p) * I4
    K1 = np.sqrt(p) * V
    return [K0, K1]


def amplitude_damping_2q(gamma: float) -> KrausList:
    """Local amplitude damping on both qubits with the same rate.

    Parameters
    ----------
    gamma : float
        Damping probability in ``[0, 1]``.

    Returns
    -------
    KrausList
        Four 4x4 Kraus operators.
    """
    ad = amplitude_damping(gamma)
    return local_noise(ad, ad)


def unitary_channel(U: NDArray) -> KrausList:
    """Wrap a unitary gate as a single-Kraus-operator channel.

    Parameters
    ----------
    U : NDArray
        Unitary matrix.

    Returns
    -------
    KrausList
        Single-element list ``[U]``.
    """
    return [U.astype(complex)]


# ---------------------------------------------------------------------------
# Channel comparison metrics
# ---------------------------------------------------------------------------


def entanglement_fidelity(kraus_ops: KrausList) -> float:
    """Entanglement fidelity: how close a channel is to the identity.

    F_e(N) = (1/d^2) sum_i |Tr(K_i)|^2

    For an ideal (identity) channel, F_e = 1.
    For a fully depolarizing channel, F_e = 1/d.

    Parameters
    ----------
    kraus_ops : KrausList
        Kraus operators.

    Returns
    -------
    float
        Entanglement fidelity in [0, 1].
    """
    d = kraus_ops[0].shape[0]
    f_e = sum(abs(np.trace(K)) ** 2 for K in kraus_ops) / d ** 2
    return float(np.clip(f_e, 0.0, 1.0))


def diamond_norm_approx(kraus_N: KrausList, kraus_M: KrausList) -> float:
    """Approximate diamond norm distance between two channels.

    Uses the Choi matrix SVD upper bound:
    ||N - M||_diamond <= d * ||Choi(N) - Choi(M)||_1

    where ||.||_1 is the trace norm (sum of singular values).

    Parameters
    ----------
    kraus_N, kraus_M : KrausList
        Kraus operators of the two channels.

    Returns
    -------
    float
        Upper bound on diamond norm distance.
    """
    d = kraus_N[0].shape[0]
    choi_N = _choi_from_kraus(kraus_N)
    choi_M = _choi_from_kraus(kraus_M)
    diff = choi_N - choi_M
    sv = np.linalg.svd(diff, compute_uv=False)
    return float(d * np.sum(sv))


def _choi_from_kraus(kraus_ops: KrausList) -> NDArray:
    """Construct Choi matrix from Kraus operators.

    Uses the Choi-Jamiolkowski isomorphism:
    J = sum_{i,j} |i><j| otimes N(|i><j|)
    so that J[i*d:(i+1)*d, j*d:(j+1)*d] = N(|i><j|).
    """
    d = kraus_ops[0].shape[0]
    choi = np.zeros((d * d, d * d), dtype=complex)
    for i in range(d):
        for j in range(d):
            E_ij = np.zeros((d, d), dtype=complex)
            E_ij[i, j] = 1.0
            N_Eij = sum(K @ E_ij @ K.conj().T for K in kraus_ops)
            choi[i * d:(i + 1) * d, j * d:(j + 1) * d] = N_Eij
    return choi


def coherent_incoherent_decomposition(kraus_ops: KrausList) -> dict:
    """Decompose channel error into coherent and incoherent parts.

    Coherent error: unitary part that could be corrected by calibration.
    Incoherent error: stochastic part that requires error correction.

    Returns
    -------
    dict
        Keys:

        - ``'coherent_frac'``: fraction of error that is coherent (0 to 1)
        - ``'incoherent_frac'``: fraction that is incoherent (0 to 1)
        - ``'unitarity'``: channel unitarity (1 = unitary, 1/d = fully depol)
        - ``'avg_gate_fidelity'``: average gate fidelity
    """
    d = kraus_ops[0].shape[0]

    # Average gate fidelity via entanglement fidelity
    f_e = entanglement_fidelity(kraus_ops)
    f_avg = (d * f_e + 1) / (d + 1)

    # Unitarity: approximate from second moment of Kraus operators
    second_moment = 0.0
    for K in kraus_ops:
        K_dag_K = K.conj().T @ K
        second_moment += np.trace(K_dag_K @ K_dag_K).real
    unitarity = (d * second_moment - 1) / (d * d - 1)
    unitarity = float(np.clip(unitarity, 0.0, 1.0))

    # Coherent vs incoherent: from unitarity and fidelity
    # r = 1 - f_avg (total error)
    # unitarity = 1 means purely coherent, unitarity = 1/d means purely incoherent
    r = 1.0 - f_avg
    if r > 1e-15:
        coherent_frac = float(np.clip(
            (unitarity - 1.0 / d) / (1.0 - 1.0 / d), 0.0, 1.0
        ))
        incoherent_frac = 1.0 - coherent_frac
    else:
        coherent_frac = 0.0
        incoherent_frac = 0.0

    return {
        'coherent_frac': coherent_frac,
        'incoherent_frac': incoherent_frac,
        'unitarity': unitarity,
        'avg_gate_fidelity': f_avg,
    }
