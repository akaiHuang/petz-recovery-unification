#!/usr/bin/env python3
"""
=============================================================================
PETZ RECOVERY MAP TOOLKIT
First open-source implementation of the Petz recovery map

Provides a complete numerical toolkit for:
  - Constructing and applying the standard Petz recovery map
  - Constructing and applying the rotated Petz map (JRSWW universal recovery)
  - Computing quantum information quantities (fidelity, relative entropy, tau)
  - Verifying the JRSWW bound F^2 >= exp(-DeltaD) [uses rotated Petz]
  - Checking saturation conditions
  - Verifying composition sub-additivity sqrt(tau_12) <= sqrt(tau_1) + sqrt(tau_2)
    [unique to the standard Petz map]

Two recovery maps are implemented:

  1. Standard Petz map (Petz 1986, 1988):
     R_{sigma,N}(X) = sigma^{1/2} N^dag(N(sigma)^{-1/2} X N(sigma)^{-1/2}) sigma^{1/2}
     - Unique Bayesian-consistent (retrodiction functor) recovery map
     - Satisfies composition sub-additivity: sqrt(tau_12) <= sqrt(tau_1) + sqrt(tau_2)
     - Satisfies the weaker Barnum-Knill bound

  2. Rotated Petz map (Junge-Renner-Sutter-Wilde-Winter 2018):
     D_{sigma,t}(X) = sigma^{(1+it)/2} N^dag(N(sigma)^{(-1-it)/2} X N(sigma)^{(-1+it)/2}) sigma^{(1-it)/2}
     D_sigma(X) = int dt beta(t) D_{sigma,t}(X)
     where beta(t) = pi / (2(1 + cosh(pi*t)))
     - Satisfies the strong JRSWW bound: F^2 >= exp(-DeltaD)
     - Does NOT satisfy composition sub-additivity in general

Dependencies: NumPy, SciPy only. No QuTiP or other quantum libraries required.

Notation follows:
  Huang (2026), "Petz Recovery Map as Retrodiction Functor:
  Equivalence of Quantum Retrodiction, Recovery, and the Crooks Theorem"

Author: Sheng-Kai Huang (2026)
License: MIT
=============================================================================
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray
from typing import Optional

# Type aliases
DensityMatrix = NDArray[np.complexfloating]
KrausList = list[NDArray[np.complexfloating]]
SuperOp = NDArray[np.complexfloating]

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

_EPS = 1e-12  # threshold for treating eigenvalues as zero


def _clip_hermitian(A: NDArray) -> NDArray:
    """Force matrix to be exactly Hermitian by averaging with adjoint."""
    return 0.5 * (A + A.conj().T)


def _safe_eigendecomposition(A: NDArray) -> tuple[NDArray, NDArray]:
    """
    Eigendecomposition of a Hermitian matrix with numerical safeguards.

    Returns eigenvalues (clipped to >= 0 for positive-semidefinite matrices)
    and the unitary matrix of eigenvectors (columns).
    """
    A = _clip_hermitian(A)
    eigvals, eigvecs = np.linalg.eigh(A)
    eigvals = np.maximum(eigvals, 0.0)
    return eigvals, eigvecs


def _matrix_power_psd(A: NDArray, power: float) -> NDArray:
    """
    Compute A^power for a positive semidefinite Hermitian matrix A.

    Uses eigendecomposition for numerical stability.  Zero eigenvalues
    are handled correctly: 0^p = 0 for p > 0, and pseudoinverse for p < 0.
    """
    eigvals, U = _safe_eigendecomposition(A)
    powered = np.zeros_like(eigvals)
    mask = eigvals > _EPS
    powered[mask] = eigvals[mask] ** power
    return (U * powered[np.newaxis, :]) @ U.conj().T


def _matrix_power_hermitian(A: NDArray, power: complex) -> NDArray:
    """
    Compute A^power for a Hermitian matrix A, where power can be complex.

    Used for modular automorphisms sigma^{it} where power = i*t.
    """
    A = _clip_hermitian(A)
    eigvals, U = np.linalg.eigh(A)
    eigvals = np.maximum(eigvals, 0.0)
    powered = np.zeros(len(eigvals), dtype=complex)
    mask = eigvals > _EPS
    powered[mask] = eigvals[mask] ** power
    return (U * powered[np.newaxis, :]) @ U.conj().T


def _matrix_sqrt_psd(A: NDArray) -> NDArray:
    """Matrix square root of a PSD matrix via eigendecomposition."""
    return _matrix_power_psd(A, 0.5)


def _matrix_log_psd(A: NDArray) -> NDArray:
    """
    Matrix logarithm of a PSD matrix.

    Returns log(A) computed via eigendecomposition.  Zero eigenvalues
    produce -inf entries (handled by np.log with a warning suppressed).
    """
    eigvals, U = _safe_eigendecomposition(A)
    with np.errstate(divide="ignore"):
        log_eigvals = np.where(eigvals > _EPS, np.log(eigvals), -np.inf)
    return (U * log_eigvals[np.newaxis, :]) @ U.conj().T


def _pseudoinverse_sqrt(A: NDArray) -> NDArray:
    """
    Compute the pseudoinverse square root of a PSD matrix:
        A^{-1/2} on supp(A),  0 on ker(A).
    """
    return _matrix_power_psd(A, -0.5)


# ---------------------------------------------------------------------------
# Core quantum channel functions
# ---------------------------------------------------------------------------


def apply_channel(rho: DensityMatrix, kraus_ops: KrausList) -> DensityMatrix:
    """
    Apply a quantum channel in Kraus representation.

    N(rho) = sum_i K_i @ rho @ K_i^dag

    Parameters
    ----------
    rho : DensityMatrix
        Input density matrix (or general operator), shape (d_in, d_in).
    kraus_ops : KrausList
        List of Kraus operators K_i, each of shape (d_out, d_in).

    Returns
    -------
    DensityMatrix
        Output N(rho), shape (d_out, d_out).
    """
    result = np.zeros((kraus_ops[0].shape[0], kraus_ops[0].shape[0]), dtype=complex)
    for K in kraus_ops:
        result += K @ rho @ K.conj().T
    return result


def adjoint_channel(X: NDArray, kraus_ops: KrausList) -> NDArray:
    """
    Apply the adjoint (Heisenberg picture) channel.

    N^dag(X) = sum_i K_i^dag @ X @ K_i

    Parameters
    ----------
    X : NDArray
        Input operator, shape (d_out, d_out).
    kraus_ops : KrausList
        List of Kraus operators K_i, each of shape (d_out, d_in).

    Returns
    -------
    NDArray
        Output operator N^dag(X), shape (d_in, d_in).
    """
    d_in = kraus_ops[0].shape[1]
    result = np.zeros((d_in, d_in), dtype=complex)
    for K in kraus_ops:
        result += K.conj().T @ X @ K
    return result


# ---------------------------------------------------------------------------
# Petz recovery map (standard)
# ---------------------------------------------------------------------------


def petz_recovery_map(
    kraus_ops: KrausList, sigma: DensityMatrix
) -> SuperOp:
    """
    Construct the standard Petz recovery map as a superoperator matrix.

    The Petz recovery map R_{sigma, N} is defined by:

        R_{sigma, N}(X) = sigma^{1/2} N^dag( N(sigma)^{-1/2} X N(sigma)^{-1/2} ) sigma^{1/2}

    where N^dag is the adjoint channel and the inverse square root is taken
    as a pseudoinverse on the support of N(sigma).

    Parameters
    ----------
    kraus_ops : KrausList
        Kraus operators {K_i} of the channel N, each shape (d_out, d_in).
    sigma : DensityMatrix
        Reference state, shape (d_in, d_in).  Must be positive semidefinite.

    Returns
    -------
    SuperOp
        Superoperator matrix of shape (d_in^2, d_out^2) in column-vectorized
        form, such that vec(R(X)) = S_R @ vec(X).

    Notes
    -----
    This is the FIRST public open-source implementation of the Petz recovery
    map.  The map was introduced by D. Petz (1986, 1988) and plays a central
    role in quantum information theory as the canonical recovery channel.

    When sigma is full-rank and N is the identity channel, R_{sigma, N} is
    the identity.  When N is a measurement channel, R reconstructs the
    Bayesian retrodiction.
    """
    d_in = sigma.shape[0]
    d_out = kraus_ops[0].shape[0]

    sigma_sqrt = _matrix_sqrt_psd(sigma)
    tau_N = apply_channel(sigma, kraus_ops)       # N(sigma)
    tau_N_inv_sqrt = _pseudoinverse_sqrt(tau_N)    # N(sigma)^{-1/2}

    # Build superoperator by action on basis: R(|i><j|) for all i,j
    S_R = np.zeros((d_in * d_in, d_out * d_out), dtype=complex)

    for i in range(d_out):
        for j in range(d_out):
            E_ij = np.zeros((d_out, d_out), dtype=complex)
            E_ij[i, j] = 1.0

            # R(E_ij) = sigma^{1/2} N^dag(tau_N^{-1/2} E_ij tau_N^{-1/2}) sigma^{1/2}
            sandwiched = tau_N_inv_sqrt @ E_ij @ tau_N_inv_sqrt
            adj_result = adjoint_channel(sandwiched, kraus_ops)
            R_Eij = sigma_sqrt @ adj_result @ sigma_sqrt

            # Column-vectorized: |i><j| -> e_{j*d_out + i}
            S_R[:, j * d_out + i] = R_Eij.flatten(order="F")

    return S_R


def apply_petz_recovery(
    X: DensityMatrix,
    kraus_ops: KrausList,
    sigma: DensityMatrix,
    S_R: Optional[SuperOp] = None,
) -> DensityMatrix:
    """
    Apply the standard Petz recovery map to an operator X.

    Parameters
    ----------
    X : DensityMatrix
        Input operator in the output space of N, shape (d_out, d_out).
    kraus_ops : KrausList
        Kraus operators of the channel N.
    sigma : DensityMatrix
        Reference state, shape (d_in, d_in).
    S_R : SuperOp, optional
        Precomputed Petz recovery superoperator.  If None, it is computed.

    Returns
    -------
    DensityMatrix
        R_{sigma,N}(X), shape (d_in, d_in).
    """
    d_in = sigma.shape[0]
    if S_R is None:
        S_R = petz_recovery_map(kraus_ops, sigma)
    vec_X = X.flatten(order="F")
    vec_RX = S_R @ vec_X
    return _clip_hermitian(vec_RX.reshape(d_in, d_in, order="F"))


# ---------------------------------------------------------------------------
# Rotated Petz recovery map (JRSWW universal recovery)
# ---------------------------------------------------------------------------


def rotated_petz_recovery(
    X: DensityMatrix,
    kraus_ops: KrausList,
    sigma: DensityMatrix,
    n_points: int = 201,
    t_max: float = 8.0,
) -> DensityMatrix:
    """
    Apply the rotated Petz recovery map (JRSWW universal recovery) to X.

    The rotated Petz map D_sigma is defined by:

        D_sigma(X) = integral dt beta(t) D_{sigma,t}(X)

    where:
        D_{sigma,t}(X) = sigma^{(1+it)/2} N^dag(N(sigma)^{(-1-it)/2} X N(sigma)^{(-1+it)/2}) sigma^{(1-it)/2}

    and beta(t) = pi / (2(1 + cosh(pi*t))) is the JRSWW weight function.

    At t=0, D_{sigma,0} reduces to the standard Petz recovery map.
    The integration over modular automorphisms handles quantum
    noncommutativity and guarantees the JRSWW bound:

        F^2(rho, D_sigma(N(rho))) >= exp(-DeltaD)

    Reference: Junge, Renner, Sutter, Wilde, Winter (2018),
    "Universal recovery maps and approximate sufficiency of quantum
    relative entropy," Ann. Henri Poincare 19, 2955.

    Parameters
    ----------
    X : DensityMatrix
        Input operator in the output space of N, shape (d_out, d_out).
    kraus_ops : KrausList
        Kraus operators {K_i} of the channel N.
    sigma : DensityMatrix
        Reference state, shape (d_in, d_in).  Must be positive definite.
    n_points : int
        Number of quadrature points for the integral over t.
        Default 201 (sufficient for qubit systems).
    t_max : float
        Integration range [-t_max, t_max].  Default 8.0.
        The weight beta(t) decays as exp(-pi|t|), so t_max=8
        captures > 99.9999% of the weight.

    Returns
    -------
    DensityMatrix
        D_sigma(X), shape (d_in, d_in).
    """
    d_in = sigma.shape[0]

    tau_N = apply_channel(sigma, kraus_ops)

    # Eigendecompose sigma and N(sigma) once for efficient complex powers
    eigvals_s, U_s = _safe_eigendecomposition(sigma)
    eigvals_t, U_t = _safe_eigendecomposition(tau_N)

    # Trapezoidal quadrature
    t_values = np.linspace(-t_max, t_max, n_points)
    dt = t_values[1] - t_values[0]
    weights = np.full(n_points, dt)
    weights[0] = dt / 2
    weights[-1] = dt / 2

    # JRSWW weight: beta(t) = pi / (2(1 + cosh(pi*t)))
    beta_values = np.pi / (2.0 * (1.0 + np.cosh(np.pi * t_values)))

    result = np.zeros((d_in, d_in), dtype=complex)

    for k, t in enumerate(t_values):
        w = weights[k] * beta_values[k]
        if abs(w) < 1e-30:
            continue

        # sigma^{(1+it)/2} and sigma^{(1-it)/2}
        sig_L = (U_s * (eigvals_s ** ((1 + 1j * t) / 2))[np.newaxis, :]) @ U_s.conj().T
        sig_R = (U_s * (eigvals_s ** ((1 - 1j * t) / 2))[np.newaxis, :]) @ U_s.conj().T

        # N(sigma)^{(-1-it)/2} and N(sigma)^{(-1+it)/2}
        # Use pseudoinverse for zero eigenvalues
        pow_L = np.where(eigvals_t > _EPS, eigvals_t ** ((-1 - 1j * t) / 2), 0.0)
        pow_R = np.where(eigvals_t > _EPS, eigvals_t ** ((-1 + 1j * t) / 2), 0.0)
        tau_L = (U_t * pow_L[np.newaxis, :]) @ U_t.conj().T
        tau_R = (U_t * pow_R[np.newaxis, :]) @ U_t.conj().T

        # D_{sigma,t}(X) = sig_L @ N^dag(tau_L @ X @ tau_R) @ sig_R
        inner = tau_L @ X @ tau_R
        adj_result = adjoint_channel(inner, kraus_ops)
        result += w * (sig_L @ adj_result @ sig_R)

    return _clip_hermitian(result)


# ---------------------------------------------------------------------------
# Quantum information quantities
# ---------------------------------------------------------------------------


def fidelity(rho: DensityMatrix, sigma: DensityMatrix) -> float:
    """
    Uhlmann fidelity between two density matrices.

    F(rho, sigma) = (Tr sqrt(sqrt(rho) sigma sqrt(rho)))^2

    Parameters
    ----------
    rho : DensityMatrix
        First density matrix.
    sigma : DensityMatrix
        Second density matrix.

    Returns
    -------
    float
        Fidelity F(rho, sigma) in [0, 1].

    Notes
    -----
    Uses eigendecomposition of sqrt(rho) @ sigma @ sqrt(rho) for
    numerical stability.  Clips result to [0, 1].
    """
    sqrt_rho = _matrix_sqrt_psd(rho)
    M = sqrt_rho @ sigma @ sqrt_rho
    M = _clip_hermitian(M)

    eigvals = np.linalg.eigvalsh(M)
    eigvals = np.maximum(eigvals, 0.0)

    F = np.sum(np.sqrt(eigvals)) ** 2
    return float(np.clip(F.real, 0.0, 1.0))


def relative_entropy(rho: DensityMatrix, sigma: DensityMatrix) -> float:
    """
    Quantum relative entropy D(rho || sigma) = Tr[rho (log rho - log sigma)].

    Computed via eigendecompositions to avoid inf * 0 issues that arise
    when rho has zero eigenvalues (e.g., pure states).

    Uses the formula:
        D(rho || sigma) = sum_i p_i log(p_i) - sum_{i,j} p_i |<u_i|v_j>|^2 log(q_j)

    where (p_i, |u_i>) are eigenvalues/vectors of rho and (q_j, |v_j>) of sigma.

    Parameters
    ----------
    rho : DensityMatrix
        State rho (density matrix).
    sigma : DensityMatrix
        State sigma (density matrix).

    Returns
    -------
    float
        D(rho || sigma).  Returns +inf if supp(rho) is not contained
        in supp(sigma).
    """
    eigvals_rho, U_rho = _safe_eigendecomposition(rho)
    eigvals_sigma, U_sigma = _safe_eigendecomposition(sigma)

    # Overlap matrix: T_{jk} = |<v_j|u_k>|^2
    overlap = U_sigma.conj().T @ U_rho
    T = np.abs(overlap) ** 2  # T[j, k] = |<v_j | u_k>|^2

    # Check support condition: supp(rho) must be in supp(sigma)
    for k in range(len(eigvals_rho)):
        if eigvals_rho[k] > _EPS:
            weight_outside = np.sum(T[eigvals_sigma <= _EPS, k])
            if weight_outside > 1e-8:
                return float("inf")

    # Tr[rho log rho] = sum_i p_i log(p_i)  (0 log 0 = 0)
    mask_rho = eigvals_rho > _EPS
    term1 = np.sum(eigvals_rho[mask_rho] * np.log(eigvals_rho[mask_rho]))

    # Tr[rho log sigma] = sum_{i,j} p_i |<u_i|v_j>|^2 log(q_j)
    mask_sigma = eigvals_sigma > _EPS
    log_q = np.zeros_like(eigvals_sigma)
    log_q[mask_sigma] = np.log(eigvals_sigma[mask_sigma])
    # sum_j T[j, k] * log(q_j) for each k, then weight by p_k
    term2 = 0.0
    for k in range(len(eigvals_rho)):
        if eigvals_rho[k] > _EPS:
            term2 += eigvals_rho[k] * np.sum(T[mask_sigma, k] * log_q[mask_sigma])

    result = float(term1 - term2)

    if result < 0 and result > -1e-10:
        result = 0.0

    return result


def tau_parameter(
    rho: DensityMatrix,
    kraus_ops: KrausList,
    sigma: DensityMatrix,
) -> float:
    """
    Compute the tau parameter measuring recovery failure (standard Petz).

    tau = 1 - F(rho, R_{sigma,N}(N(rho)))

    tau = 0: perfect recovery (no time arrow).
    tau = 1: complete failure (maximal time arrow).

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
        tau in [0, 1].
    """
    omega = apply_channel(rho, kraus_ops)
    S_R = petz_recovery_map(kraus_ops, sigma)
    recovered = apply_petz_recovery(omega, kraus_ops, sigma, S_R)
    F = fidelity(rho, recovered)
    return float(np.clip(1.0 - F, 0.0, 1.0))


def tau_parameter_rotated(
    rho: DensityMatrix,
    kraus_ops: KrausList,
    sigma: DensityMatrix,
    n_points: int = 201,
) -> float:
    """
    Compute the tau parameter using the rotated Petz map.

    tau_rot = 1 - F(rho, D_sigma(N(rho)))

    The JRSWW bound guarantees: F^2 >= exp(-DeltaD), i.e.,
    tau_rot <= 1 - sqrt(exp(-DeltaD)).

    Parameters
    ----------
    rho : DensityMatrix
        Input state.
    kraus_ops : KrausList
        Kraus operators of the channel N.
    sigma : DensityMatrix
        Reference state.
    n_points : int
        Number of quadrature points for the rotated Petz integral.

    Returns
    -------
    float
        tau_rot in [0, 1].
    """
    omega = apply_channel(rho, kraus_ops)
    recovered = rotated_petz_recovery(omega, kraus_ops, sigma, n_points=n_points)
    F = fidelity(rho, recovered)
    return float(np.clip(1.0 - F, 0.0, 1.0))


def delta_D(
    rho: DensityMatrix,
    sigma: DensityMatrix,
    kraus_ops: KrausList,
) -> float:
    """
    Relative entropy decrease under channel N.

    DeltaD = D(rho || sigma) - D(N(rho) || N(sigma))

    By the data processing inequality, DeltaD >= 0.

    Parameters
    ----------
    rho : DensityMatrix
        Input state.
    sigma : DensityMatrix
        Reference state.
    kraus_ops : KrausList
        Kraus operators of the channel N.

    Returns
    -------
    float
        DeltaD >= 0.
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


# ---------------------------------------------------------------------------
# Verification functions
# ---------------------------------------------------------------------------


def verify_jrsww_bound(
    rho: DensityMatrix,
    sigma: DensityMatrix,
    kraus_ops: KrausList,
    n_points: int = 201,
) -> dict:
    """
    Verify the JRSWW bound using the rotated Petz map:

        F(rho, R_tilde(N(rho))) >= exp(-DeltaD)

    where F is the Uhlmann fidelity (squared convention):
        F(rho, sigma) = (Tr sqrt(sqrt(rho) sigma sqrt(rho)))^2.

    The bound holds for the rotated (twirled) Petz map, which averages
    the standard Petz map over modular automorphisms.  See
    Junge-Renner-Sutter-Wilde-Winter (2018), Corollary 5.1.

    Note: In the convention where F_root = Tr sqrt(sqrt(rho) sigma sqrt(rho)),
    the bound is F_root^2 >= exp(-DeltaD).  Since our fidelity() already
    returns F_root^2, we compare fidelity() directly with exp(-DeltaD).

    Also reports the standard Petz map fidelity for comparison.

    Parameters
    ----------
    rho : DensityMatrix
        Input state.
    sigma : DensityMatrix
        Reference state (should be full-rank for finite DeltaD).
    kraus_ops : KrausList
        Kraus operators of the channel N.
    n_points : int
        Quadrature points for the rotated Petz integral.

    Returns
    -------
    dict
        Keys:
        - 'fidelity_rotated': F(rho, R_tilde(N(rho))) [squared convention]
        - 'fidelity_standard': F(rho, R_Petz(N(rho))) [squared convention]
        - 'exp_neg_delta_D': exp(-DeltaD)
        - 'delta_D': DeltaD
        - 'gap_rotated': F_rot - exp(-DeltaD) (should be >= 0)
        - 'gap_standard': F_std - exp(-DeltaD) (may be < 0)
        - 'bound_holds': whether JRSWW bound holds for rotated Petz
        - 'tau_standard': 1 - F_std (standard Petz tau)
        - 'tau_rotated': 1 - F_rot (rotated Petz tau)
    """
    # Standard Petz: fidelity() returns F = (Tr sqrt(M))^2
    tau_std = tau_parameter(rho, kraus_ops, sigma)
    F_std = 1.0 - tau_std  # This is the squared fidelity

    # Rotated Petz
    tau_rot = tau_parameter_rotated(rho, kraus_ops, sigma, n_points=n_points)
    F_rot = 1.0 - tau_rot  # This is the squared fidelity

    # DeltaD
    dD = delta_D(rho, sigma, kraus_ops)
    exp_neg_dD = np.exp(-dD) if not np.isinf(dD) else 0.0

    # JRSWW bound: F(rho, R(N(rho))) >= exp(-DeltaD)
    # where F is the squared fidelity (our convention)
    gap_rot = F_rot - exp_neg_dD
    gap_std = F_std - exp_neg_dD
    bound_holds = gap_rot >= -1e-8

    return {
        "fidelity_rotated": F_rot,
        "fidelity_standard": F_std,
        "exp_neg_delta_D": exp_neg_dD,
        "delta_D": dD,
        "gap_rotated": gap_rot,
        "gap_standard": gap_std,
        "bound_holds": bool(bound_holds),
        "tau_standard": tau_std,
        "tau_rotated": tau_rot,
    }


def check_saturation_conditions(
    rho: DensityMatrix,
    sigma: DensityMatrix,
    kraus_ops: KrausList,
    n_points: int = 201,
) -> dict:
    """
    Check conditions for saturation of the JRSWW bound.

    When saturation occurs (DeltaD = 0, i.e., the channel is sufficient
    for distinguishing rho from sigma), both the standard and rotated Petz
    maps achieve perfect recovery.

    For sigma = I/d and pure rho, near-saturation is related to:
      - [N(rho), N(sigma)] ~ 0 (commutativity of outputs)

    Parameters
    ----------
    rho : DensityMatrix
        Input state (ideally pure for saturation analysis).
    sigma : DensityMatrix
        Reference state (ideally I/d).
    kraus_ops : KrausList
        Kraus operators of the channel N.
    n_points : int
        Quadrature points for the rotated Petz integral.

    Returns
    -------
    dict
        Diagnostic information about saturation conditions.
    """
    d = rho.shape[0]
    omega = apply_channel(rho, kraus_ops)
    tau_N = apply_channel(sigma, kraus_ops)

    purity = np.trace(rho @ rho).real
    is_pure = abs(purity - 1.0) < 1e-8
    sigma_is_mm = np.linalg.norm(sigma - np.eye(d) / d) < 1e-8

    # Commutator [N(rho), N(sigma)]
    commutator = omega @ tau_N - tau_N @ omega
    comm_norm = np.linalg.norm(commutator, "fro")

    # Eigenvalue ratio analysis
    eigvals_omega, U_omega = _safe_eigendecomposition(omega)
    tau_N_in_basis = U_omega.conj().T @ tau_N @ U_omega
    support_mask = eigvals_omega > _EPS
    ratios = []
    for idx in np.where(support_mask)[0]:
        ratios.append(tau_N_in_basis[idx, idx].real / eigvals_omega[idx])
    ratios_constant = (
        (max(ratios) - min(ratios)) < 1e-8 if len(ratios) > 1 else True
    )

    result = verify_jrsww_bound(rho, sigma, kraus_ops, n_points)

    return {
        "is_pure": is_pure,
        "sigma_is_maximally_mixed": sigma_is_mm,
        "omega": omega,
        "tau_N": tau_N,
        "commutator_norm": comm_norm,
        "eigenvalue_ratios_constant": ratios_constant,
        "eigenvalue_ratios": ratios,
        "gap_rotated": result["gap_rotated"],
        "gap_standard": result["gap_standard"],
        "fidelity_rotated": result["fidelity_rotated"],
        "fidelity_standard": result["fidelity_standard"],
        "exp_neg_delta_D": result["exp_neg_delta_D"],
    }


def verify_composition(
    kraus1: KrausList,
    kraus2: KrausList,
    rho: DensityMatrix,
    sigma: DensityMatrix,
) -> dict:
    """
    Verify the composition sub-additivity for the standard Petz map:

        sqrt(tau_12) <= sqrt(tau_1) + sqrt(tau_2_eff)

    This property is unique to the standard Petz map and follows from its
    functoriality (Bayesian consistency).  The rotated Petz map does NOT
    satisfy this in general.

    Parameters
    ----------
    kraus1 : KrausList
        Kraus operators of the first channel N_1.
    kraus2 : KrausList
        Kraus operators of the second channel N_2.
    rho : DensityMatrix
        Input state.
    sigma : DensityMatrix
        Reference state.

    Returns
    -------
    dict
        Keys:
        - 'tau_composed': tau for N_2 o N_1
        - 'tau_1': tau for N_1
        - 'tau_2_eff': tau for N_2 with evolved states
        - 'sqrt_tau_composed': sqrt(tau_12)
        - 'sqrt_tau_1_plus_sqrt_tau_2': sqrt(tau_1) + sqrt(tau_2_eff)
        - 'slack': RHS - LHS (should be >= 0)
        - 'bound_holds': bool
    """
    kraus_composed = [K2 @ K1 for K2 in kraus2 for K1 in kraus1]

    tau_12 = tau_parameter(rho, kraus_composed, sigma)
    tau_1 = tau_parameter(rho, kraus1, sigma)

    rho_after_1 = apply_channel(rho, kraus1)
    sigma_after_1 = apply_channel(sigma, kraus1)
    tau_2_eff = tau_parameter(rho_after_1, kraus2, sigma_after_1)

    sqrt_tau_12 = np.sqrt(max(tau_12, 0.0))
    sqrt_tau_1 = np.sqrt(max(tau_1, 0.0))
    sqrt_tau_2 = np.sqrt(max(tau_2_eff, 0.0))

    rhs = sqrt_tau_1 + sqrt_tau_2
    slack = rhs - sqrt_tau_12
    bound_holds = slack >= -1e-10

    return {
        "tau_composed": tau_12,
        "tau_1": tau_1,
        "tau_2_eff": tau_2_eff,
        "sqrt_tau_composed": sqrt_tau_12,
        "sqrt_tau_1_plus_sqrt_tau_2": rhs,
        "slack": slack,
        "bound_holds": bool(bound_holds),
    }


# ---------------------------------------------------------------------------
# Random generation
# ---------------------------------------------------------------------------


def random_cptp_map(
    d_in: int,
    d_out: int,
    rank: Optional[int] = None,
    rng: Optional[np.random.Generator] = None,
) -> KrausList:
    """
    Generate a random CPTP map via Stinespring dilation.

    Samples a random isometry V: C^{d_in} -> C^{d_out} (x) C^{d_env}
    and extracts Kraus operators by partial-tracing the environment.

    Parameters
    ----------
    d_in : int
        Input dimension.
    d_out : int
        Output dimension.
    rank : int, optional
        Number of Kraus operators (environment dimension).
        Defaults to d_in * d_out.
    rng : np.random.Generator, optional
        Random number generator for reproducibility.

    Returns
    -------
    KrausList
        List of Kraus operators, each shape (d_out, d_in).
    """
    if rng is None:
        rng = np.random.default_rng()
    if rank is None:
        rank = d_in * d_out

    G = rng.standard_normal((d_out * rank, d_in)) + 1j * rng.standard_normal(
        (d_out * rank, d_in)
    )
    Q, _ = np.linalg.qr(G, mode="reduced")

    return [Q[k * d_out : (k + 1) * d_out, :] for k in range(rank)]


def random_pure_state(
    d: int, rng: Optional[np.random.Generator] = None
) -> DensityMatrix:
    """
    Generate a random pure state |psi><psi| (Haar-distributed).

    Parameters
    ----------
    d : int
        Hilbert space dimension.
    rng : np.random.Generator, optional
        Random number generator.

    Returns
    -------
    DensityMatrix
        Rank-1 density matrix |psi><psi|, shape (d, d).
    """
    if rng is None:
        rng = np.random.default_rng()
    psi = rng.standard_normal(d) + 1j * rng.standard_normal(d)
    psi /= np.linalg.norm(psi)
    return np.outer(psi, psi.conj())


def random_density_matrix(
    d: int,
    rank: Optional[int] = None,
    rng: Optional[np.random.Generator] = None,
) -> DensityMatrix:
    """
    Generate a random density matrix of given rank.

    Uses the partial trace method: rho = Tr_env(|psi><psi|).

    Parameters
    ----------
    d : int
        Hilbert space dimension.
    rank : int, optional
        Rank of the density matrix.  Defaults to d (full rank).
    rng : np.random.Generator, optional
        Random number generator.

    Returns
    -------
    DensityMatrix
        Density matrix of shape (d, d) with the specified rank.
    """
    if rng is None:
        rng = np.random.default_rng()
    if rank is None:
        rank = d

    G = rng.standard_normal((d, rank)) + 1j * rng.standard_normal((d, rank))
    rho = G @ G.conj().T
    rho /= np.trace(rho)
    return _clip_hermitian(rho)


# ---------------------------------------------------------------------------
# Batch verification
# ---------------------------------------------------------------------------


def batch_verify_jrsww(
    n_samples: int = 1000,
    d: int = 2,
    sigma_type: str = "maximally_mixed",
    seed: int = 42,
    n_points: int = 201,
) -> dict:
    """
    Batch verification of the JRSWW bound F^2 >= exp(-DeltaD).

    Uses the rotated Petz map, for which the bound is a theorem
    (Junge et al. 2018, Corollary 5.1).

    Parameters
    ----------
    n_samples : int
        Number of random samples.
    d : int
        Hilbert space dimension.
    sigma_type : str
        'maximally_mixed' for sigma = I/d, 'random' for random full-rank sigma.
    seed : int
        Random seed for reproducibility.
    n_points : int
        Quadrature points for the rotated Petz integral.

    Returns
    -------
    dict
        Verification results including gap statistics.
    """
    rng = np.random.default_rng(seed)
    gaps_rot = []
    gaps_std = []
    violations = 0

    for _ in range(n_samples):
        rho = random_density_matrix(d, rng=rng)
        if sigma_type == "maximally_mixed":
            sigma = np.eye(d, dtype=complex) / d
        else:
            sigma = random_density_matrix(d, rng=rng)

        kraus_ops = random_cptp_map(d, d, rng=rng)
        result = verify_jrsww_bound(rho, sigma, kraus_ops, n_points)

        gaps_rot.append(result["gap_rotated"])
        gaps_std.append(result["gap_standard"])
        if not result["bound_holds"]:
            violations += 1

    gaps_rot = np.array(gaps_rot)
    gaps_std = np.array(gaps_std)

    return {
        "all_pass": violations == 0,
        "n_violations": violations,
        "min_gap_rotated": float(np.min(gaps_rot)),
        "max_gap_rotated": float(np.max(gaps_rot)),
        "mean_gap_rotated": float(np.mean(gaps_rot)),
        "min_gap_standard": float(np.min(gaps_std)),
        "max_gap_standard": float(np.max(gaps_std)),
        "mean_gap_standard": float(np.mean(gaps_std)),
        "n_standard_violations": int(np.sum(gaps_std < -1e-8)),
        "n_samples": n_samples,
    }


def batch_verify_saturation(
    n_samples: int = 1000,
    d: int = 2,
    seed: int = 42,
    n_points: int = 201,
) -> dict:
    """
    Batch verification of saturation conditions.

    For sigma = I/d and pure rho, checks the relationship between
    the commutator norm ||[N(rho), N(sigma)]|| and the gap
    F^2 - exp(-DeltaD).

    Parameters
    ----------
    n_samples : int
        Number of random samples.
    d : int
        Hilbert space dimension.
    seed : int
        Random seed for reproducibility.
    n_points : int
        Quadrature points for the rotated Petz integral.

    Returns
    -------
    dict
        Saturation analysis results.
    """
    rng = np.random.default_rng(seed)
    sigma = np.eye(d, dtype=complex) / d

    commuting_gaps = []
    noncommuting_gaps = []
    comm_norms = []
    all_gaps = []

    comm_threshold = 1e-8

    for _ in range(n_samples):
        rho = random_pure_state(d, rng=rng)
        kraus_ops = random_cptp_map(d, d, rng=rng)

        result = check_saturation_conditions(rho, sigma, kraus_ops, n_points)
        comm_norm = result["commutator_norm"]
        gap = result["gap_rotated"]

        comm_norms.append(comm_norm)
        all_gaps.append(gap)

        if comm_norm < comm_threshold:
            commuting_gaps.append(gap)
        else:
            noncommuting_gaps.append(gap)

    comm_norms = np.array(comm_norms)
    all_gaps = np.array(all_gaps)
    if len(comm_norms) > 1 and np.std(comm_norms) > 0 and np.std(all_gaps) > 0:
        correlation = float(np.corrcoef(comm_norms, all_gaps)[0, 1])
    else:
        correlation = float("nan")

    return {
        "n_samples": n_samples,
        "n_commuting": len(commuting_gaps),
        "n_noncommuting": len(noncommuting_gaps),
        "commuting_gaps": commuting_gaps,
        "noncommuting_gaps": noncommuting_gaps,
        "commuting_mean_gap": float(np.mean(commuting_gaps)) if commuting_gaps else float("nan"),
        "noncommuting_mean_gap": float(np.mean(noncommuting_gaps)) if noncommuting_gaps else float("nan"),
        "correlation": correlation,
    }


def batch_verify_composition(
    n_samples: int = 1000,
    d: int = 2,
    seed: int = 42,
) -> dict:
    """
    Batch verification of composition sub-additivity (standard Petz only):

        sqrt(tau_12) <= sqrt(tau_1) + sqrt(tau_2_eff)

    Parameters
    ----------
    n_samples : int
        Number of random samples.
    d : int
        Hilbert space dimension.
    seed : int
        Random seed for reproducibility.

    Returns
    -------
    dict
        Composition bound verification results.
    """
    rng = np.random.default_rng(seed)
    slacks = []
    violations = 0

    sigma = np.eye(d, dtype=complex) / d

    for _ in range(n_samples):
        rho = random_density_matrix(d, rng=rng)
        kraus1 = random_cptp_map(d, d, rng=rng)
        kraus2 = random_cptp_map(d, d, rng=rng)

        result = verify_composition(kraus1, kraus2, rho, sigma)
        slacks.append(result["slack"])
        if not result["bound_holds"]:
            violations += 1

    slacks = np.array(slacks)
    return {
        "all_pass": violations == 0,
        "n_violations": violations,
        "min_slack": float(np.min(slacks)),
        "max_slack": float(np.max(slacks)),
        "mean_slack": float(np.mean(slacks)),
        "n_samples": n_samples,
    }


# ---------------------------------------------------------------------------
# Main: run all batch verifications
# ---------------------------------------------------------------------------


def main() -> None:
    """Run all three batch verifications and print results."""
    N = 1000
    d = 2

    print("=" * 72)
    print("PETZ RECOVERY MAP TOOLKIT -- Batch Verification Suite")
    print("First open-source implementation of the Petz recovery map")
    print(f"Testing with {N} samples in d = {d}")
    print("=" * 72)

    # --- JRSWW bound ---
    print("\n" + "-" * 72)
    print("TEST 1: JRSWW Bound  F^2(rho, R_tilde(N(rho))) >= exp(-DeltaD)")
    print("  Uses the rotated Petz map (JRSWW universal recovery).")
    print("  Also reports standard Petz for comparison.")
    print("-" * 72)

    print("\n  (a) sigma = I/d (maximally mixed):")
    result_mm = batch_verify_jrsww(N, d, sigma_type="maximally_mixed")
    print(f"      JRSWW bound (rotated Petz):")
    print(f"        All pass: {result_mm['all_pass']}")
    print(f"        Violations: {result_mm['n_violations']}/{N}")
    print(f"        Gap range: [{result_mm['min_gap_rotated']:.2e}, {result_mm['max_gap_rotated']:.2e}]")
    print(f"        Mean gap:  {result_mm['mean_gap_rotated']:.6f}")
    print(f"      Standard Petz (for comparison):")
    print(f"        Violations of exp(-DeltaD): {result_mm['n_standard_violations']}/{N}")
    print(f"        Gap range: [{result_mm['min_gap_standard']:.2e}, {result_mm['max_gap_standard']:.2e}]")
    print(f"        Mean gap:  {result_mm['mean_gap_standard']:.6f}")

    print("\n  (b) sigma = random full-rank state:")
    result_rnd = batch_verify_jrsww(N, d, sigma_type="random")
    print(f"      JRSWW bound (rotated Petz):")
    print(f"        All pass: {result_rnd['all_pass']}")
    print(f"        Violations: {result_rnd['n_violations']}/{N}")
    print(f"        Gap range: [{result_rnd['min_gap_rotated']:.2e}, {result_rnd['max_gap_rotated']:.2e}]")
    print(f"        Mean gap:  {result_rnd['mean_gap_rotated']:.6f}")
    print(f"      Standard Petz (for comparison):")
    print(f"        Violations of exp(-DeltaD): {result_rnd['n_standard_violations']}/{N}")
    print(f"        Gap range: [{result_rnd['min_gap_standard']:.2e}, {result_rnd['max_gap_standard']:.2e}]")

    # --- Saturation conditions ---
    print("\n" + "-" * 72)
    print("TEST 2: Saturation Analysis (sigma = I/d, pure rho)")
    print("-" * 72)

    result_sat = batch_verify_saturation(N, d)
    print(f"  Commuting cases:    {result_sat['n_commuting']}/{N}")
    print(f"  Noncommuting cases: {result_sat['n_noncommuting']}/{N}")
    if result_sat["n_commuting"] > 0:
        print(f"  Mean gap (commuting):    {result_sat['commuting_mean_gap']:.2e}")
    if result_sat["n_noncommuting"] > 0:
        print(f"  Mean gap (noncommuting): {result_sat['noncommuting_mean_gap']:.2e}")
    print(f"  Correlation(||[N(rho),N(sigma)]||, gap): {result_sat['correlation']:.4f}")

    if result_sat["n_commuting"] > 0:
        comm_gaps = np.array(result_sat["commuting_gaps"])
        near_saturated = np.sum(np.abs(comm_gaps) < 1e-6)
        print(f"  Commuting near saturation (gap < 1e-6): "
              f"{near_saturated}/{result_sat['n_commuting']}")

    # --- Composition bound ---
    print("\n" + "-" * 72)
    print("TEST 3: Composition  sqrt(tau_12) <= sqrt(tau_1) + sqrt(tau_2_eff)")
    print("  Uses the standard Petz map (unique to Petz via functoriality).")
    print("-" * 72)

    result_comp = batch_verify_composition(N, d)
    print(f"  All pass:   {result_comp['all_pass']}")
    print(f"  Violations: {result_comp['n_violations']}/{N}")
    print(f"  Slack range: [{result_comp['min_slack']:.2e}, {result_comp['max_slack']:.2e}]")
    print(f"  Mean slack:  {result_comp['mean_slack']:.6f}")

    # --- Summary ---
    print("\n" + "=" * 72)
    print("SUMMARY")
    print("=" * 72)
    all_pass = (
        result_mm["all_pass"]
        and result_rnd["all_pass"]
        and result_comp["all_pass"]
    )
    if all_pass:
        print("  ALL BOUNDS VERIFIED SUCCESSFULLY")
        total = 3 * N
        print(f"  Total tests: {total}")
        print("  Zero violations across all tests.")
    else:
        print("  WARNING: Some bounds were violated!")
        total_v = (
            result_mm["n_violations"]
            + result_rnd["n_violations"]
            + result_comp["n_violations"]
        )
        print(f"  Total violations: {total_v}")
    print("=" * 72)


if __name__ == "__main__":
    main()
