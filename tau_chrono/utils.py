"""
Matrix utilities for quantum information computations.

Provides numerically stable implementations of matrix square root,
inverse square root, logarithm, and trace norm for positive semidefinite
matrices via eigendecomposition.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

_EPS = 1e-12  # threshold for treating eigenvalues as zero


def clip_hermitian(A: NDArray) -> NDArray:
    """Force a matrix to be exactly Hermitian by averaging with its adjoint.

    Parameters
    ----------
    A : NDArray
        Input square matrix.

    Returns
    -------
    NDArray
        Hermitian matrix ``0.5 * (A + A^dag)``.
    """
    return 0.5 * (A + A.conj().T)


def safe_eigendecomposition(A: NDArray) -> tuple[NDArray, NDArray]:
    """Eigendecompose a Hermitian matrix with numerical safeguards.

    Clips eigenvalues to >= 0 (for PSD matrices).

    Parameters
    ----------
    A : NDArray
        Hermitian matrix.

    Returns
    -------
    tuple[NDArray, NDArray]
        ``(eigenvalues, eigenvectors)`` where eigenvectors are column vectors.
    """
    A = clip_hermitian(A)
    eigvals, eigvecs = np.linalg.eigh(A)
    eigvals = np.maximum(eigvals, 0.0)
    return eigvals, eigvecs


def matrix_power_psd(A: NDArray, power: float) -> NDArray:
    """Compute ``A^power`` for a positive semidefinite Hermitian matrix.

    Uses eigendecomposition for numerical stability.  Zero eigenvalues
    are handled correctly: ``0^p = 0`` for ``p > 0``, pseudoinverse for ``p < 0``.

    Parameters
    ----------
    A : NDArray
        PSD Hermitian matrix.
    power : float
        Exponent.

    Returns
    -------
    NDArray
        Matrix ``A^power``.
    """
    eigvals, U = safe_eigendecomposition(A)
    powered = np.zeros_like(eigvals)
    mask = eigvals > _EPS
    powered[mask] = eigvals[mask] ** power
    return (U * powered[np.newaxis, :]) @ U.conj().T


def matrix_sqrt(A: NDArray) -> NDArray:
    """Matrix square root of a PSD matrix.

    Parameters
    ----------
    A : NDArray
        PSD Hermitian matrix.

    Returns
    -------
    NDArray
        ``A^{1/2}``.
    """
    return matrix_power_psd(A, 0.5)


def matrix_inv_sqrt(A: NDArray) -> NDArray:
    """Pseudoinverse square root of a PSD matrix.

    Returns ``A^{-1/2}`` on the support of *A* and zero on the kernel.

    Parameters
    ----------
    A : NDArray
        PSD Hermitian matrix.

    Returns
    -------
    NDArray
        ``A^{-1/2}`` (pseudoinverse).
    """
    return matrix_power_psd(A, -0.5)


def matrix_log(A: NDArray) -> NDArray:
    """Matrix logarithm of a PSD matrix via eigendecomposition.

    Zero eigenvalues produce ``-inf`` entries.

    Parameters
    ----------
    A : NDArray
        PSD Hermitian matrix.

    Returns
    -------
    NDArray
        ``log(A)``.
    """
    eigvals, U = safe_eigendecomposition(A)
    with np.errstate(divide="ignore"):
        log_eigvals = np.where(eigvals > _EPS, np.log(eigvals), -np.inf)
    return (U * log_eigvals[np.newaxis, :]) @ U.conj().T


def trace_norm(A: NDArray) -> float:
    """Trace norm (Schatten 1-norm) of a matrix.

    Parameters
    ----------
    A : NDArray
        Input matrix.

    Returns
    -------
    float
        ``||A||_1 = Tr(sqrt(A^dag A))``.
    """
    svs = np.linalg.svd(A, compute_uv=False)
    return float(np.sum(svs))


def commutator_norm(A: NDArray, B: NDArray) -> float:
    """Frobenius norm of the commutator ``[A, B] = AB - BA``.

    Parameters
    ----------
    A, B : NDArray
        Square matrices of the same shape.

    Returns
    -------
    float
        ``||[A, B]||_F``.
    """
    comm = A @ B - B @ A
    return float(np.linalg.norm(comm, "fro"))
