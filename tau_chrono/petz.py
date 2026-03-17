"""
Petz recovery map, fidelity, and tau computation.

Core quantum information primitives:
    - ``apply_channel``: Kraus-representation channel application
    - ``adjoint_channel``: Heisenberg-picture adjoint
    - ``petz_recovery_map``: standard Petz recovery superoperator
    - ``apply_petz_recovery``: apply the Petz map to an operator
    - ``fidelity``: Uhlmann fidelity F(rho, sigma)
    - ``tau_parameter``: recovery failure ``tau = 1 - F(rho, R(N(rho)))``
    - ``relative_entropy``: quantum relative entropy D(rho || sigma)
    - ``delta_D``: relative entropy decrease under a channel
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray
from typing import Optional, List

from .utils import (
    clip_hermitian,
    safe_eigendecomposition,
    matrix_sqrt,
    matrix_inv_sqrt,
    _EPS,
)

DensityMatrix = NDArray[np.complexfloating]
KrausList = List[NDArray[np.complexfloating]]
SuperOp = NDArray[np.complexfloating]


# ---------------------------------------------------------------------------
# Channel application
# ---------------------------------------------------------------------------


def apply_channel(rho: DensityMatrix, kraus_ops: KrausList) -> DensityMatrix:
    """Apply a quantum channel in Kraus representation.

    ``N(rho) = sum_i K_i rho K_i^dag``

    Parameters
    ----------
    rho : DensityMatrix
        Input density matrix of shape ``(d_in, d_in)``.
    kraus_ops : KrausList
        Kraus operators, each of shape ``(d_out, d_in)``.

    Returns
    -------
    DensityMatrix
        Output ``N(rho)`` of shape ``(d_out, d_out)``.
    """
    result = np.zeros((kraus_ops[0].shape[0], kraus_ops[0].shape[0]), dtype=complex)
    for K in kraus_ops:
        result += K @ rho @ K.conj().T
    return result


def adjoint_channel(X: NDArray, kraus_ops: KrausList) -> NDArray:
    """Apply the adjoint (Heisenberg picture) channel.

    ``N^dag(X) = sum_i K_i^dag X K_i``

    Parameters
    ----------
    X : NDArray
        Operator in the output space, shape ``(d_out, d_out)``.
    kraus_ops : KrausList
        Kraus operators, each of shape ``(d_out, d_in)``.

    Returns
    -------
    NDArray
        ``N^dag(X)`` of shape ``(d_in, d_in)``.
    """
    d_in = kraus_ops[0].shape[1]
    result = np.zeros((d_in, d_in), dtype=complex)
    for K in kraus_ops:
        result += K.conj().T @ X @ K
    return result


# ---------------------------------------------------------------------------
# Petz recovery map
# ---------------------------------------------------------------------------


def petz_recovery_map(kraus_ops: KrausList, sigma: DensityMatrix) -> SuperOp:
    """Construct the standard Petz recovery map as a superoperator matrix.

    ``R_{sigma,N}(X) = sigma^{1/2} N^dag(N(sigma)^{-1/2} X N(sigma)^{-1/2}) sigma^{1/2}``

    Parameters
    ----------
    kraus_ops : KrausList
        Kraus operators of the channel N.
    sigma : DensityMatrix
        Reference state (PSD).

    Returns
    -------
    SuperOp
        Column-vectorised superoperator of shape ``(d_in^2, d_out^2)``.
    """
    d_in = sigma.shape[0]
    d_out = kraus_ops[0].shape[0]

    sigma_sqrt = matrix_sqrt(sigma)
    tau_N = apply_channel(sigma, kraus_ops)
    tau_N_inv_sqrt = matrix_inv_sqrt(tau_N)

    S_R = np.zeros((d_in * d_in, d_out * d_out), dtype=complex)
    for i in range(d_out):
        for j in range(d_out):
            E_ij = np.zeros((d_out, d_out), dtype=complex)
            E_ij[i, j] = 1.0
            sandwiched = tau_N_inv_sqrt @ E_ij @ tau_N_inv_sqrt
            adj_result = adjoint_channel(sandwiched, kraus_ops)
            R_Eij = sigma_sqrt @ adj_result @ sigma_sqrt
            S_R[:, j * d_out + i] = R_Eij.flatten(order="F")
    return S_R


def apply_petz_recovery(
    X: DensityMatrix,
    kraus_ops: KrausList,
    sigma: DensityMatrix,
    S_R: Optional[SuperOp] = None,
) -> DensityMatrix:
    """Apply the standard Petz recovery map to an operator.

    Parameters
    ----------
    X : DensityMatrix
        Operator in the output space of N.
    kraus_ops : KrausList
        Kraus operators of the channel N.
    sigma : DensityMatrix
        Reference state.
    S_R : SuperOp, optional
        Precomputed superoperator (if None it is computed on the fly).

    Returns
    -------
    DensityMatrix
        ``R_{sigma,N}(X)``.
    """
    d_in = sigma.shape[0]
    if S_R is None:
        S_R = petz_recovery_map(kraus_ops, sigma)
    vec_X = X.flatten(order="F")
    vec_RX = S_R @ vec_X
    return clip_hermitian(vec_RX.reshape(d_in, d_in, order="F"))


# ---------------------------------------------------------------------------
# Quantum information quantities
# ---------------------------------------------------------------------------


def fidelity(rho: DensityMatrix, sigma: DensityMatrix) -> float:
    """Uhlmann fidelity (squared convention).

    ``F(rho, sigma) = (Tr sqrt(sqrt(rho) sigma sqrt(rho)))^2``

    Parameters
    ----------
    rho, sigma : DensityMatrix
        Density matrices.

    Returns
    -------
    float
        Fidelity in ``[0, 1]``.
    """
    sqrt_rho = matrix_sqrt(rho)
    M = clip_hermitian(sqrt_rho @ sigma @ sqrt_rho)
    eigvals = np.linalg.eigvalsh(M)
    eigvals = np.maximum(eigvals, 0.0)
    F = np.sum(np.sqrt(eigvals)) ** 2
    return float(np.clip(F.real, 0.0, 1.0))


def tau_parameter(
    rho: DensityMatrix,
    kraus_ops: KrausList,
    sigma: DensityMatrix,
) -> float:
    """Recovery failure parameter for the standard Petz map.

    ``tau = 1 - F(rho, R_{sigma,N}(N(rho)))``

    ``tau = 0`` means perfect recovery (no time arrow).
    ``tau = 1`` means complete failure (maximal time arrow).

    Parameters
    ----------
    rho : DensityMatrix
        Input state.
    kraus_ops : KrausList
        Kraus operators of the channel N.
    sigma : DensityMatrix
        Reference state.

    Returns
    -------
    float
        ``tau`` in ``[0, 1]``.
    """
    omega = apply_channel(rho, kraus_ops)
    S_R = petz_recovery_map(kraus_ops, sigma)
    recovered = apply_petz_recovery(omega, kraus_ops, sigma, S_R)
    F = fidelity(rho, recovered)
    return float(np.clip(1.0 - F, 0.0, 1.0))


def relative_entropy(rho: DensityMatrix, sigma: DensityMatrix) -> float:
    """Quantum relative entropy ``D(rho || sigma) = Tr[rho (log rho - log sigma)]``.

    Returns ``+inf`` when ``supp(rho)`` is not contained in ``supp(sigma)``.

    Parameters
    ----------
    rho, sigma : DensityMatrix
        Density matrices.

    Returns
    -------
    float
        ``D(rho || sigma)``.
    """
    eigvals_rho, U_rho = safe_eigendecomposition(rho)
    eigvals_sigma, U_sigma = safe_eigendecomposition(sigma)

    overlap = U_sigma.conj().T @ U_rho
    T = np.abs(overlap) ** 2

    for k in range(len(eigvals_rho)):
        if eigvals_rho[k] > _EPS:
            weight_outside = np.sum(T[eigvals_sigma <= _EPS, k])
            if weight_outside > 1e-8:
                return float("inf")

    mask_rho = eigvals_rho > _EPS
    term1 = np.sum(eigvals_rho[mask_rho] * np.log(eigvals_rho[mask_rho]))

    mask_sigma = eigvals_sigma > _EPS
    log_q = np.zeros_like(eigvals_sigma)
    log_q[mask_sigma] = np.log(eigvals_sigma[mask_sigma])
    term2 = 0.0
    for k in range(len(eigvals_rho)):
        if eigvals_rho[k] > _EPS:
            term2 += eigvals_rho[k] * np.sum(T[mask_sigma, k] * log_q[mask_sigma])

    result = float(term1 - term2)
    if result < 0 and result > -1e-10:
        result = 0.0
    return result


def delta_D(
    rho: DensityMatrix,
    sigma: DensityMatrix,
    kraus_ops: KrausList,
) -> float:
    """Relative entropy decrease under channel N.

    ``DeltaD = D(rho || sigma) - D(N(rho) || N(sigma)) >= 0``

    Parameters
    ----------
    rho, sigma : DensityMatrix
        Input states.
    kraus_ops : KrausList
        Kraus operators.

    Returns
    -------
    float
        ``DeltaD >= 0``.
    """
    D_before = relative_entropy(rho, sigma)
    omega = apply_channel(rho, kraus_ops)
    tau_N = apply_channel(sigma, kraus_ops)
    D_after = relative_entropy(omega, tau_N)

    if np.isinf(D_before) or np.isinf(D_after):
        return float("inf")

    result = D_before - D_after
    if result < 0 and result > -1e-10:
        result = 0.0
    return float(result)
