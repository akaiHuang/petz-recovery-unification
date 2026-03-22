"""
solver_magnus.py -- Magnus expansion solver for CMB Boltzmann equations on GPU (MLX).

MATHEMATICAL BREAKTHROUGH: Instead of ~600 RK4 or implicit steps per k-mode,
use the matrix exponential of the linear system to take O(10-50) large steps.

The Boltzmann equations are LINEAR:  y'(tau) = A(tau) * y(tau)

where y = [eta, delta_c, delta_b, theta_b, F_g0...F_g{lg}, F_n0...F_n{ln}]
and A(tau) depends ONLY on background quantities (calH, a, R, kappa_dot).

The EXACT solution over [tau_1, tau_2] is:
    y(tau_2) = expm(Omega) * y(tau_1)

where Omega is the Magnus expansion:
    1st order: Omega_1 = integral A(tau) dtau  (trapezoidal: (dt/2)(A1+A2))
    2nd order: Omega_2 = (dt^2/12) [A1, A2]  (commutator correction)

STRATEGY FOR STIFFNESS:
    Thomson scattering makes |kappa_dot| ~ 10^7 at early times, creating
    eigenvalues of A up to ~10^11. The matrix exponential handles this EXACTLY
    (no CFL condition!) -- the stiff modes just decay exponentially. But we
    need sufficient scaling-and-squaring iterations.

    The key insight: the stiff eigenvalues are all NEGATIVE REAL (damping).
    exp(lambda * dt) -> 0 for lambda << 0, so the stiff modes self-annihilate.
    We just need enough Taylor terms and squaring to capture this correctly.

    Two-phase approach:
    Phase 1 (tight coupling, tau < tau_switch): Use TCA-0 reduced system
        with only 5 + ln_max + 1 variables. No photon hierarchy stiffness.
    Phase 2 (free streaming, tau > tau_switch): Full hierarchy with Magnus.
        After recombination, |kappa_dot| drops to ~0 and eigenvalues are O(k).

Pipeline:
  1. Background (numpy, ~10ms)
  2. Phase 1: TCA-0 propagation on GPU (small system, fast)
  3. Phase 2: Full Magnus propagation on GPU (full hierarchy)
  4. Gauge transform + LOS integration -> C_l

Usage:
  python -m mlx_class.solver_magnus [--N_k 500] [--N_steps 30] [--order 2]

Author: Sheng-Kai Huang, 2026
"""
import os
os.environ['MPLCONFIGDIR'] = '/tmp/mpl_config'

import numpy as np
import mlx.core as mx
import time
import sys

from .background import (
    Background, A_s, n_s, k_pivot, T_CMB,
    Omega_r as _OMEGA_R,
    Omega_b as _OMEGA_B,
    Omega_c as _OMEGA_C,
    Omega_m as _OMEGA_M,
    H0_Mpc as _H0_MPC,
)
from .perturbations_neutrino import _f_nu, _OMEGA_GAMMA, _OMEGA_NU
from .perturbations_sync import (
    adiabatic_ic_sync, n_var_sync, idx_fn_start, idx_fg,
    gauge_transform, diagnose_metric,
    IDX_ETA, IDX_DELTA_C, IDX_DELTA_B, IDX_THETA_B, IDX_FG_START,
    L_GAMMA_MAX, L_NU_MAX_SYNC,
)

_H02 = float(_H0_MPC ** 2)
_Og = float(_OMEGA_GAMMA)
_On = float(_OMEGA_NU)
_Ob = float(_OMEGA_B)
_Oc = float(_OMEGA_C)


# ============================================================================
# Build the TCA-0 coefficient matrix A_tca(tau) for reduced system
# ============================================================================
# TCA-0 state: [eta, delta_c, delta_b, theta_b, delta_g, F_n0...F_n{ln}]
# nvar_tca = 5 + ln_max + 1
#
# In tight coupling: theta_gamma = theta_b, F_g,l = 0 for l >= 2
# delta_g evolved: delta_g' = -(4/3) theta_b - (2/3) h'
# theta_b evolved: combined baryon-photon momentum equation:
#   theta_b' = [-R*calH*theta_b + (k^2/4)*delta_g] / (1 + R)

_T_E = 0   # eta
_T_DC = 1  # delta_c
_T_DB = 2  # delta_b
_T_TB = 3  # theta_b
_T_DG = 4  # delta_g
_T_FN = 5  # neutrino hierarchy starts here


def _nv_tca(ln_max):
    """Number of variables in TCA-0 system."""
    return 5 + ln_max + 1


def build_A_tca(k, calH, a, R, kappa_dot, ln_max):
    """
    Build coefficient matrix for TCA-0 reduced system.

    In tight coupling, the photon hierarchy collapses:
    - theta_gamma = theta_b (coupled)
    - F_g,1 = 4*theta_b/(3*k), F_g,l = 0 for l >= 2
    - Combined momentum: theta_b' = [-R*calH*theta_b + k^2/4*delta_g] / (1+R)

    Returns A_tca: (nv, nv) matrix.
    """
    nv = _nv_tca(ln_max)
    k2 = k * k
    ia = 1.0 / a
    ia2 = ia * ia

    A = np.zeros((nv, nv), dtype=np.float64)

    # h' coefficients in TCA-0 state:
    # h' = (2/calH)[k^2 eta + (3/2)H0^2(Og/a^2 dg + On/a^2 dn + Ob/a db + Oc/a dc)]
    h_c = np.zeros(nv)
    h_c[_T_E] = (2.0 / calH) * k2
    h_c[_T_DC] = (2.0 / calH) * 1.5 * _H02 * _Oc * ia
    h_c[_T_DB] = (2.0 / calH) * 1.5 * _H02 * _Ob * ia
    h_c[_T_DG] = (2.0 / calH) * 1.5 * _H02 * _Og * ia2
    h_c[_T_FN] = (2.0 / calH) * 1.5 * _H02 * _On * ia2

    # eta' coefficients:
    # theta_g = theta_b in TCA, theta_n = (3/4)k F_n,1
    eta_c = np.zeros(nv)
    eta_c[_T_TB] = 1.5 * _H02 / k2 * (_Ob * ia + (4.0/3.0) * _Og * ia2)
    if ln_max >= 1:
        eta_c[_T_FN + 1] = 1.5 * _H02 / k2 * (4.0/3.0) * _On * ia2 * 0.75 * k

    # Row 0: eta' = eta_c . y
    A[_T_E, :] = eta_c

    # Row 1: delta_c' = -h'/2
    A[_T_DC, :] = -0.5 * h_c

    # Row 2: delta_b' = -theta_b - h'/2
    A[_T_DB, :] = -0.5 * h_c
    A[_T_DB, _T_TB] += -1.0

    # Row 3: theta_b' = [-R*calH*theta_b + k^2/4*delta_g] / (1+R)
    inv_1pR = 1.0 / (1.0 + R)
    A[_T_TB, _T_TB] = -R * calH * inv_1pR
    A[_T_TB, _T_DG] = k2 / 4.0 * inv_1pR

    # Row 4: delta_g' = -(4/3)*theta_b - (2/3)*h'
    A[_T_DG, _T_TB] = -4.0 / 3.0
    A[_T_DG, :] += -(2.0 / 3.0) * h_c

    # Neutrino hierarchy (same as full system, no collisions)
    fn = _T_FN

    # l=0: F_n,0' = -k F_n,1 - (2/3) h'
    A[fn, :] += -(2.0 / 3.0) * h_c
    if ln_max >= 1:
        A[fn, fn + 1] += -k

    # l=1: F_n,1' = (k/3)(F_n,0 - 2F_n,2)
    if ln_max >= 1:
        A[fn + 1, fn] = k / 3.0
        if ln_max >= 2:
            A[fn + 1, fn + 2] = -2.0 * k / 3.0

    # l=2: F_n,2' = (2k/5)F_n,1 - (3k/5)F_n,3 + (4/15)h' + (8/15)eta'
    if ln_max >= 2:
        A[fn + 2, fn + 1] = 2.0 * k / 5.0
        if ln_max >= 3:
            A[fn + 2, fn + 3] = -3.0 * k / 5.0
        A[fn + 2, :] += (4.0 / 15.0) * h_c
        A[fn + 2, :] += (8.0 / 15.0) * eta_c

    # l=3..ln_max-1: free streaming
    for ell in range(3, ln_max):
        idx_n = fn + ell
        A[idx_n, idx_n - 1] = k * ell / (2.0 * ell + 1.0)
        A[idx_n, idx_n + 1] = -k * (ell + 1.0) / (2.0 * ell + 1.0)

    # l=ln_max: truncation
    if ln_max >= 1:
        idx_nlm = fn + ln_max
        A[idx_nlm, idx_nlm - 1] = k * ln_max / (2.0 * ln_max + 1.0)
        A[idx_nlm, idx_nlm] = -0.01 * k

    return A


def tca_ic(k, tau_init, ln_max):
    """
    Build TCA-0 initial conditions from standard adiabatic ICs.

    Maps the full sync gauge ICs to the TCA-0 reduced state.
    """
    from .perturbations_sync import adiabatic_ic_sync as _aic
    y_full = _aic(k, tau_init, L_GAMMA_MAX, ln_max)

    nv = _nv_tca(ln_max)
    fn_s_full = idx_fn_start(L_GAMMA_MAX)

    y_tca = np.zeros(nv, dtype=np.float64)
    y_tca[_T_E] = y_full[IDX_ETA]
    y_tca[_T_DC] = y_full[IDX_DELTA_C]
    y_tca[_T_DB] = y_full[IDX_DELTA_B]
    y_tca[_T_TB] = y_full[IDX_THETA_B]
    y_tca[_T_DG] = y_full[IDX_FG_START]  # delta_gamma = F_g,0
    # Neutrinos
    y_tca[_T_FN:_T_FN + ln_max + 1] = y_full[fn_s_full:fn_s_full + ln_max + 1]

    return y_tca


def tca_to_full(y_tca, k, calH, a, R, kappa_dot, lg_max, ln_max):
    """
    Promote TCA-0 state to full hierarchy state at the switching time.

    Photon monopole and dipole are set from TCA state.
    Higher photon multipoles initialized from their quasi-static equilibrium:
    F_g,l ~ 0 for l >= 2 (Thomson damped).
    """
    nvar = n_var_sync(lg_max, ln_max)
    fn_s = idx_fn_start(lg_max)

    y = np.zeros(nvar, dtype=np.float64)
    y[IDX_ETA] = y_tca[_T_E]
    y[IDX_DELTA_C] = y_tca[_T_DC]
    y[IDX_DELTA_B] = y_tca[_T_DB]
    y[IDX_THETA_B] = y_tca[_T_TB]

    # Photon hierarchy from TCA
    y[IDX_FG_START] = y_tca[_T_DG]           # F_g,0 = delta_gamma
    y[IDX_FG_START + 1] = 4.0 * y_tca[_T_TB] / (3.0 * k)  # F_g,1 = 4*theta_b/(3*k)
    # F_g,l = 0 for l >= 2 (Thomson damped -- quasi-static equilibrium)
    # Optionally set F_g,2 from its quasi-static value:
    if lg_max >= 2:
        abs_kd = abs(kappa_dot)
        if abs_kd > 1e-10:
            # From dF_g,2/dt = 0: F_g,2 ~ (2k/5)F_g,1 / ((9/10)|kd|)
            y[IDX_FG_START + 2] = (2.0 * k / 5.0) * y[IDX_FG_START + 1] / (
                (9.0 / 10.0) * abs_kd + 0.01 * k)

    # Neutrinos: copy directly
    y[fn_s:fn_s + ln_max + 1] = y_tca[_T_FN:_T_FN + ln_max + 1]

    return y


# ============================================================================
# Vectorized A-matrix builders (all k at once for one tau)
# ============================================================================

def build_A_matrices_vectorized(k_arr, calH, a, R, kappa_dot, lg_max, ln_max):
    """
    Build A matrices for ALL k-modes at a given tau, vectorized.

    Returns: (N_k, nvar, nvar) array.
    Much faster than looping build_A_matrix over k.
    """
    N_k = len(k_arr)
    nvar = n_var_sync(lg_max, ln_max)
    fn_s = idx_fn_start(lg_max)
    ig0 = IDX_FG_START

    k = k_arr
    k2 = k * k
    ia = 1.0 / a
    ia2 = ia * ia
    abs_kd = abs(kappa_dot)

    A = np.zeros((N_k, nvar, nvar), dtype=np.float32)

    # h' coefficients: (N_k, nvar)
    h_c = np.zeros((N_k, nvar), dtype=np.float32)
    h_c[:, IDX_ETA] = (2.0 / calH) * k2
    h_c[:, IDX_DELTA_C] = (2.0 / calH) * 1.5 * _H02 * _Oc * ia
    h_c[:, IDX_DELTA_B] = (2.0 / calH) * 1.5 * _H02 * _Ob * ia
    h_c[:, ig0] = (2.0 / calH) * 1.5 * _H02 * _Og * ia2
    h_c[:, fn_s] = (2.0 / calH) * 1.5 * _H02 * _On * ia2

    # eta' coefficients: (N_k, nvar)
    eta_c = np.zeros((N_k, nvar), dtype=np.float32)
    eta_c[:, IDX_THETA_B] = 1.5 * _H02 / k2 * _Ob * ia
    eta_c[:, ig0 + 1] = 1.5 * _H02 / k2 * (4.0/3.0) * _Og * ia2 * 0.75 * k
    if ln_max >= 1:
        eta_c[:, fn_s + 1] = 1.5 * _H02 / k2 * (4.0/3.0) * _On * ia2 * 0.75 * k

    # Row 0: eta'
    A[:, IDX_ETA, :] = eta_c

    # Row 1: delta_c' = -h'/2
    A[:, IDX_DELTA_C, :] = -0.5 * h_c

    # Row 2: delta_b' = -theta_b - h'/2
    A[:, IDX_DELTA_B, :] = -0.5 * h_c
    A[:, IDX_DELTA_B, IDX_THETA_B] += -1.0

    # Row 3: theta_b'
    R_safe = max(R, 1e-30)
    A[:, IDX_THETA_B, IDX_THETA_B] = -calH - abs_kd / R_safe
    A[:, IDX_THETA_B, ig0 + 1] = abs_kd / R_safe * 0.75 * k

    # Photon l=0
    A[:, ig0, :] += -(2.0 / 3.0) * h_c
    A[:, ig0, ig0 + 1] += -k

    # Photon l=1
    A[:, ig0 + 1, ig0] = k / 3.0
    if lg_max >= 2:
        A[:, ig0 + 1, ig0 + 2] = -2.0 * k / 3.0
    A[:, ig0 + 1, ig0 + 1] += -abs_kd
    A[:, ig0 + 1, IDX_THETA_B] += abs_kd * 4.0 / (3.0 * k)

    # Photon l=2
    if lg_max >= 2:
        A[:, ig0 + 2, ig0 + 1] = 2.0 * k / 5.0
        if lg_max >= 3:
            A[:, ig0 + 2, ig0 + 3] = -3.0 * k / 5.0
        A[:, ig0 + 2, :] += (4.0 / 15.0) * h_c
        A[:, ig0 + 2, :] += (8.0 / 15.0) * eta_c
        A[:, ig0 + 2, ig0 + 2] += -(9.0 / 10.0) * abs_kd

    # Photon l=3..lg_max-1
    for ell in range(3, lg_max):
        idx = ig0 + ell
        A[:, idx, idx - 1] = k * ell / (2.0 * ell + 1.0)
        A[:, idx, idx + 1] = -k * (ell + 1.0) / (2.0 * ell + 1.0)
        A[:, idx, idx] += -abs_kd

    # Photon l=lg_max
    if lg_max >= 3:
        idx_lm = ig0 + lg_max
        A[:, idx_lm, idx_lm - 1] = k * lg_max / (2.0 * lg_max + 1.0)
        A[:, idx_lm, idx_lm] = -abs_kd - 0.01 * k

    # Neutrino l=0
    A[:, fn_s, :] += -(2.0 / 3.0) * h_c
    if ln_max >= 1:
        A[:, fn_s, fn_s + 1] += -k

    # Neutrino l=1
    if ln_max >= 1:
        A[:, fn_s + 1, fn_s] = k / 3.0
        if ln_max >= 2:
            A[:, fn_s + 1, fn_s + 2] = -2.0 * k / 3.0

    # Neutrino l=2
    if ln_max >= 2:
        A[:, fn_s + 2, fn_s + 1] = 2.0 * k / 5.0
        if ln_max >= 3:
            A[:, fn_s + 2, fn_s + 3] = -3.0 * k / 5.0
        A[:, fn_s + 2, :] += (4.0 / 15.0) * h_c
        A[:, fn_s + 2, :] += (8.0 / 15.0) * eta_c

    # Neutrino l=3..ln_max-1
    for ell in range(3, ln_max):
        idx_n = fn_s + ell
        A[:, idx_n, idx_n - 1] = k * ell / (2.0 * ell + 1.0)
        A[:, idx_n, idx_n + 1] = -k * (ell + 1.0) / (2.0 * ell + 1.0)

    # Neutrino l=ln_max
    if ln_max >= 1:
        idx_nlm = fn_s + ln_max
        A[:, idx_nlm, idx_nlm - 1] = k * ln_max / (2.0 * ln_max + 1.0)
        A[:, idx_nlm, idx_nlm] = -0.01 * k

    return A


def build_A_tca_vectorized(k_arr, calH, a, R, kappa_dot, ln_max):
    """
    Build TCA-0 A matrices for ALL k-modes at a given tau, vectorized.

    Returns: (N_k, nv, nv) array.
    """
    N_k = len(k_arr)
    nv = _nv_tca(ln_max)
    fn = _T_FN

    k = k_arr
    k2 = k * k
    ia = 1.0 / a
    ia2 = ia * ia

    A = np.zeros((N_k, nv, nv), dtype=np.float32)

    # h' coefficients
    h_c = np.zeros((N_k, nv), dtype=np.float32)
    h_c[:, _T_E] = (2.0 / calH) * k2
    h_c[:, _T_DC] = (2.0 / calH) * 1.5 * _H02 * _Oc * ia
    h_c[:, _T_DB] = (2.0 / calH) * 1.5 * _H02 * _Ob * ia
    h_c[:, _T_DG] = (2.0 / calH) * 1.5 * _H02 * _Og * ia2
    h_c[:, fn] = (2.0 / calH) * 1.5 * _H02 * _On * ia2

    # eta' coefficients
    eta_c = np.zeros((N_k, nv), dtype=np.float32)
    eta_c[:, _T_TB] = 1.5 * _H02 / k2 * (_Ob * ia + (4.0/3.0) * _Og * ia2)
    if ln_max >= 1:
        eta_c[:, fn + 1] = 1.5 * _H02 / k2 * (4.0/3.0) * _On * ia2 * 0.75 * k

    inv_1pR = 1.0 / (1.0 + R)

    A[:, _T_E, :] = eta_c
    A[:, _T_DC, :] = -0.5 * h_c
    A[:, _T_DB, :] = -0.5 * h_c
    A[:, _T_DB, _T_TB] += -1.0
    A[:, _T_TB, _T_TB] = -R * calH * inv_1pR
    A[:, _T_TB, _T_DG] = k2 / 4.0 * inv_1pR
    A[:, _T_DG, _T_TB] = -4.0 / 3.0
    A[:, _T_DG, :] += -(2.0 / 3.0) * h_c

    # Neutrinos
    A[:, fn, :] += -(2.0 / 3.0) * h_c
    if ln_max >= 1:
        A[:, fn, fn + 1] += -k
        A[:, fn + 1, fn] = k / 3.0
        if ln_max >= 2:
            A[:, fn + 1, fn + 2] = -2.0 * k / 3.0

    if ln_max >= 2:
        A[:, fn + 2, fn + 1] = 2.0 * k / 5.0
        if ln_max >= 3:
            A[:, fn + 2, fn + 3] = -3.0 * k / 5.0
        A[:, fn + 2, :] += (4.0 / 15.0) * h_c
        A[:, fn + 2, :] += (8.0 / 15.0) * eta_c

    for ell in range(3, ln_max):
        idx_n = fn + ell
        A[:, idx_n, idx_n - 1] = k * ell / (2.0 * ell + 1.0)
        A[:, idx_n, idx_n + 1] = -k * (ell + 1.0) / (2.0 * ell + 1.0)

    if ln_max >= 1:
        idx_nlm = fn + ln_max
        A[:, idx_nlm, idx_nlm - 1] = k * ln_max / (2.0 * ln_max + 1.0)
        A[:, idx_nlm, idx_nlm] = -0.01 * k

    return A


# ============================================================================
# Build the FULL coefficient matrix A(tau) for sync gauge (single k)
# ============================================================================

def build_A_matrix(k, calH, a, R, kappa_dot, lg_max, ln_max):
    """
    Build the coefficient matrix A for the full sync gauge Boltzmann system.

    y' = A(tau) * y  where y = [eta, dc, db, tb, Fg0..Fg{lg}, Fn0..Fn{ln}]

    The matrix A encodes all linear couplings including Einstein constraints
    (h' and eta' are linear in y and absorbed into A).

    Parameters
    ----------
    k : float, wavenumber (Mpc^-1)
    calH : float, conformal Hubble (Mpc^-1)
    a : float, scale factor
    R : float, baryon loading 3*rho_b/(4*rho_gamma)
    kappa_dot : float, Thomson rate (negative by convention)
    lg_max, ln_max : int, hierarchy truncation

    Returns
    -------
    A : ndarray (nvar, nvar)
    """
    nvar = n_var_sync(lg_max, ln_max)
    fn_s = idx_fn_start(lg_max)
    k2 = k * k
    ia = 1.0 / a
    ia2 = ia * ia
    abs_kd = abs(kappa_dot)

    A = np.zeros((nvar, nvar), dtype=np.float64)

    # h' coefficients
    h_c = np.zeros(nvar)
    h_c[IDX_ETA] = (2.0 / calH) * k2
    h_c[IDX_DELTA_C] = (2.0 / calH) * 1.5 * _H02 * _Oc * ia
    h_c[IDX_DELTA_B] = (2.0 / calH) * 1.5 * _H02 * _Ob * ia
    h_c[IDX_FG_START] = (2.0 / calH) * 1.5 * _H02 * _Og * ia2
    h_c[fn_s] = (2.0 / calH) * 1.5 * _H02 * _On * ia2

    # eta' coefficients
    eta_c = np.zeros(nvar)
    eta_c[IDX_THETA_B] = 1.5 * _H02 / k2 * _Ob * ia
    eta_c[IDX_FG_START + 1] = 1.5 * _H02 / k2 * (4.0/3.0) * _Og * ia2 * 0.75 * k
    if ln_max >= 1:
        eta_c[fn_s + 1] = 1.5 * _H02 / k2 * (4.0/3.0) * _On * ia2 * 0.75 * k

    # Row 0: eta'
    A[IDX_ETA, :] = eta_c

    # Row 1: delta_c' = -h'/2
    A[IDX_DELTA_C, :] = -0.5 * h_c

    # Row 2: delta_b' = -theta_b - h'/2
    A[IDX_DELTA_B, :] = -0.5 * h_c
    A[IDX_DELTA_B, IDX_THETA_B] += -1.0

    # Row 3: theta_b' = -calH*theta_b + |kd|/R*(theta_g - theta_b)
    A[IDX_THETA_B, IDX_THETA_B] = -calH - abs_kd / max(R, 1e-30)
    A[IDX_THETA_B, IDX_FG_START + 1] = abs_kd / max(R, 1e-30) * 0.75 * k

    # Photon hierarchy
    ig0 = IDX_FG_START

    # l=0: F_g,0' = -k F_g,1 - (2/3) h'
    A[ig0, :] += -(2.0 / 3.0) * h_c
    A[ig0, ig0 + 1] += -k

    # l=1: F_g,1' = (k/3)(F_g,0 - 2F_g,2) + |kd|(-F_g,1 + 4 theta_b/(3k))
    ig1 = ig0 + 1
    A[ig1, ig0] = k / 3.0
    if lg_max >= 2:
        A[ig1, ig0 + 2] = -2.0 * k / 3.0
    A[ig1, ig1] += -abs_kd
    A[ig1, IDX_THETA_B] += abs_kd * 4.0 / (3.0 * k)

    # l=2: F_g,2' = (2k/5)F_g,1 - (3k/5)F_g,3 + (4/15)h' + (8/15)eta' - (9/10)|kd|F_g,2
    if lg_max >= 2:
        ig2 = ig0 + 2
        A[ig2, ig1] = 2.0 * k / 5.0
        if lg_max >= 3:
            A[ig2, ig0 + 3] = -3.0 * k / 5.0
        A[ig2, :] += (4.0 / 15.0) * h_c
        A[ig2, :] += (8.0 / 15.0) * eta_c
        A[ig2, ig2] += -(9.0 / 10.0) * abs_kd

    # l=3..lg_max-1: streaming + Thomson
    for ell in range(3, lg_max):
        idx = ig0 + ell
        A[idx, idx - 1] = k * ell / (2.0 * ell + 1.0)
        A[idx, idx + 1] = -k * (ell + 1.0) / (2.0 * ell + 1.0)
        A[idx, idx] += -abs_kd

    # l=lg_max: truncation
    if lg_max >= 3:
        idx_lm = ig0 + lg_max
        A[idx_lm, idx_lm - 1] = k * lg_max / (2.0 * lg_max + 1.0)
        A[idx_lm, idx_lm] = -abs_kd - 0.01 * k

    # Neutrino hierarchy (no collisions)
    # l=0: F_n,0' = -k F_n,1 - (2/3) h'
    A[fn_s, :] += -(2.0 / 3.0) * h_c
    if ln_max >= 1:
        A[fn_s, fn_s + 1] += -k

    # l=1
    if ln_max >= 1:
        A[fn_s + 1, fn_s] = k / 3.0
        if ln_max >= 2:
            A[fn_s + 1, fn_s + 2] = -2.0 * k / 3.0

    # l=2
    if ln_max >= 2:
        A[fn_s + 2, fn_s + 1] = 2.0 * k / 5.0
        if ln_max >= 3:
            A[fn_s + 2, fn_s + 3] = -3.0 * k / 5.0
        A[fn_s + 2, :] += (4.0 / 15.0) * h_c
        A[fn_s + 2, :] += (8.0 / 15.0) * eta_c

    # l=3..ln_max-1
    for ell in range(3, ln_max):
        idx_n = fn_s + ell
        A[idx_n, idx_n - 1] = k * ell / (2.0 * ell + 1.0)
        A[idx_n, idx_n + 1] = -k * (ell + 1.0) / (2.0 * ell + 1.0)

    # l=ln_max: truncation
    if ln_max >= 1:
        idx_nlm = fn_s + ln_max
        A[idx_nlm, idx_nlm - 1] = k * ln_max / (2.0 * ln_max + 1.0)
        A[idx_nlm, idx_nlm] = -0.01 * k

    return A


# ============================================================================
# Build adaptive time grid
# ============================================================================

def build_magnus_time_grid(bg, N_steps=300, tau_start=None, tau_end=None):
    """
    Build a non-uniform time grid for full-hierarchy Magnus steps.

    Three regions:
    1. tau_start -> tau_vis_lo:  dense linear (captures early ISW + streaming)
    2. tau_vis_lo -> tau_vis_hi: VERY dense (~1 Mpc, captures visibility peak)
    3. tau_vis_hi -> tau_end:    sparse (late ISW, reionization)

    Since the full hierarchy captures ALL physics (Silk damping, streaming,
    Thomson scattering), we need dense sampling where A(tau) changes most:
    the decoupling epoch around tau_rec +/- 50 Mpc.
    """
    if tau_start is None:
        tau_start = bg.tau_grid[1]
    if tau_end is None:
        tau_end = bg.tau_0 * 0.98

    peak_idx = np.argmax(bg.visibility_grid)
    tau_rec = bg.tau_grid[peak_idx]

    g_peak = bg.visibility_grid[peak_idx]
    half_max = g_peak / 2.0
    above = bg.visibility_grid > half_max
    if np.any(above):
        tau_above = bg.tau_grid[above]
        fwhm = tau_above[-1] - tau_above[0]
    else:
        fwhm = 20.0
    sigma = fwhm / 2.355

    # Visibility region: generous +/- 5 sigma around peak
    tau_vis_lo = max(tau_rec - 5.0 * sigma, tau_start + 1.0)
    tau_vis_hi = min(tau_rec + 5.0 * sigma, tau_end * 0.5)

    # Allocate: 30% pre-visibility, 60% visibility, 10% late
    N_pre = max(int(0.30 * N_steps), 20)
    N_vis = max(int(0.60 * N_steps), 60)
    N_late = max(N_steps - N_pre - N_vis, 10)

    parts = []

    # Pre-visibility: tau_start -> tau_vis_lo (linear, dense enough for ISW)
    if tau_vis_lo > tau_start + 1.0:
        tau_pre = np.linspace(tau_start, tau_vis_lo, N_pre + 1)
        parts.append(tau_pre)

    # Visibility region: VERY DENSE linear grid (~0.5-1 Mpc spacing)
    tau_vis_grid = np.linspace(tau_vis_lo, tau_vis_hi, N_vis + 1)
    parts.append(tau_vis_grid)

    # Late: tau_vis_hi -> tau_end (sparse)
    tau_late = np.linspace(tau_vis_hi, tau_end, N_late + 1)
    parts.append(tau_late)

    tau_all = np.sort(np.unique(np.concatenate(parts)))
    return tau_all


def build_tca_time_grid(bg, N_steps_tca=30):
    """
    Build time grid for TCA-0 phase: from tau_init to tau_switch.

    Switch EARLY (|kd| ~ 5) so that:
    1. TCA-0 is still very accurate (|kd| >> k_max = 0.35)
    2. The full-hierarchy Magnus can take over with dt ~ 0.7 Mpc steps

    The TCA phase only covers the very early radiation era where
    the photon-baryon fluid is tightly coupled.
    """
    tau_init = bg.tau_grid[1]

    # Switch at |kd| ~ 5: tight coupling is still excellent here
    # (|kd|/k = 14 at k=0.35, so TCA-0 error is ~ (k/|kd|)^2 ~ 0.5%)
    kd_grid = np.abs(bg.kappa_dot_grid)
    threshold = 5.0  # Mpc^-1

    # Find where kd drops below threshold
    idx_switch = np.searchsorted(-kd_grid, -threshold)
    if idx_switch >= len(bg.tau_grid) - 1:
        idx_switch = len(bg.tau_grid) // 2
    tau_switch = bg.tau_grid[min(idx_switch, len(bg.tau_grid) - 1)]

    # Safety bounds
    peak_idx = np.argmax(bg.visibility_grid)
    tau_rec = bg.tau_grid[peak_idx]
    tau_switch = min(tau_switch, tau_rec - 100.0)  # well before recombination
    tau_switch = max(tau_switch, tau_init * 10.0)

    # Build log-spaced grid for TCA phase
    tau_tca = np.geomspace(tau_init, tau_switch, N_steps_tca + 1)

    return tau_tca, tau_switch


# ============================================================================
# Batched matrix exponential on GPU
# ============================================================================

def batched_matrix_exp(M, n_terms=16):
    """
    Compute exp(M) for a batch of matrices M: shape (N_k, n, n).

    Uses scaling-and-squaring with Taylor series.

    For stability with large negative eigenvalues (Thomson damping),
    we use aggressive scaling to ensure ||M_scaled|| << 1.
    """
    N_k, n, _ = M.shape

    # Estimate norm (max absolute element as proxy for spectral radius)
    norms = mx.max(mx.abs(M.reshape(N_k, -1)), axis=1)  # (N_k,)

    # Use uniform scaling (max over all k) for GPU efficiency
    max_norm = float(mx.max(norms))

    # Scale so that ||M_scaled|| < 0.5 (conservative for convergence)
    s_max = max(int(np.ceil(np.log2(max_norm + 1e-10))) + 1, 0)
    s_max = min(s_max, 30)  # cap at 2^30 to avoid infinite loops

    scale = 1.0 / (2.0 ** s_max) if s_max > 0 else 1.0
    M_scaled = M * scale

    # Taylor series: exp(M_s) = I + M_s + M_s^2/2! + ...
    I = mx.broadcast_to(mx.eye(n, dtype=M.dtype), (N_k, n, n))
    result = I + M_scaled
    M_power = M_scaled

    for k_ord in range(2, n_terms + 1):
        M_power = mx.matmul(M_power, M_scaled) * (1.0 / k_ord)
        result = result + M_power

    # Repeated squaring: exp(M) = (exp(M/2^s))^{2^s}
    for _ in range(s_max):
        result = mx.matmul(result, result)

    return result


# ============================================================================
# Full Magnus solver (two-phase)
# ============================================================================

def solve_magnus(N_k=500, k_min=3e-4, k_max=0.35,
                 N_steps_tca=20, N_steps_full=300,
                 lg_max=L_GAMMA_MAX, ln_max=L_NU_MAX_SYNC,
                 order=2, n_terms=16,
                 verbose=True):
    """
    Solve Boltzmann equations using Magnus expansion on GPU.

    Two-phase approach:
      Phase 1 (TCA-0): tau_init -> tau_switch, reduced system (no photon hierarchy stiffness)
      Phase 2 (full):  tau_switch -> tau_end, full hierarchy with Magnus

    Parameters
    ----------
    N_k : int
        Number of k-modes.
    N_steps_tca : int
        Number of Magnus steps in TCA phase.
    N_steps_full : int
        Number of Magnus steps in full phase.
    order : int (1 or 2)
        Magnus expansion order.
    n_terms : int
        Taylor series terms for matrix exponential.
    verbose : bool
        Print progress.

    Returns
    -------
    dict with C_l, D_l, timing, etc.
    """
    t_total = time.time()

    # ==================================================================
    # Step 1: Background
    # ==================================================================
    if verbose:
        print("=" * 70)
        print("Magnus Expansion Boltzmann Solver (GPU, Two-Phase)")
        print("=" * 70)
        print(f"\n--- Step 1: Background ---")

    t0 = time.time()
    bg = Background(khronon=False, recombination='recfast')
    bg.solve()
    t_bg = time.time() - t0
    if verbose:
        print(f"Background: {t_bg:.2f}s")

    # ==================================================================
    # Step 2: Setup
    # ==================================================================
    if verbose:
        print(f"\n--- Step 2: Setup ---")
        sys.stdout.flush()

    t0 = time.time()
    k_arr = np.geomspace(k_min, k_max, N_k).astype(np.float64)
    nvar_full = n_var_sync(lg_max, ln_max)
    nvar_tca = _nv_tca(ln_max)

    # Build TCA time grid
    tau_tca, tau_switch = build_tca_time_grid(bg, N_steps_tca)
    N_tau_tca = len(tau_tca)

    # Build full Magnus time grid (from tau_switch to tau_end)
    tau_end = bg.tau_0 * 0.98
    tau_full = build_magnus_time_grid(bg, N_steps_full,
                                       tau_start=tau_switch, tau_end=tau_end)
    N_tau_full = len(tau_full)

    # Combined tau grid for snapshots
    tau_all = np.sort(np.unique(np.concatenate([tau_tca, tau_full])))
    N_tau_all = len(tau_all)

    if verbose:
        print(f"N_k={N_k}, order={order}, n_terms={n_terms}")
        print(f"TCA phase: {N_tau_tca} points, tau=[{tau_tca[0]:.4f}, {tau_tca[-1]:.1f}] Mpc")
        print(f"  tau_switch = {tau_switch:.1f} Mpc (|kd|={abs(float(bg.kappa_dot_at_tau(np.array([tau_switch]))[0])):.2f})")
        print(f"Full phase: {N_tau_full} points, tau=[{tau_full[0]:.1f}, {tau_full[-1]:.1f}] Mpc")
        print(f"State: TCA={nvar_tca} vars, Full={nvar_full} vars")
        print(f"Total snapshots: {N_tau_all}")

    # Evaluate background on TCA grid
    calH_tca = bg.calH_at_tau(tau_tca)
    a_tca = bg.a_at_tau(tau_tca)
    R_tca = bg.R_at_tau(tau_tca)
    kd_tca = bg.kappa_dot_at_tau(tau_tca)

    # Evaluate background on full grid
    calH_full = bg.calH_at_tau(tau_full)
    a_full = bg.a_at_tau(tau_full)
    R_full = bg.R_at_tau(tau_full)
    kd_full = bg.kappa_dot_at_tau(tau_full)

    t_setup = time.time() - t0
    if verbose:
        print(f"Setup: {t_setup:.2f}s")

    # ==================================================================
    # Step 3: Phase 1 -- TCA-0 Magnus propagation
    # ==================================================================
    if verbose:
        print(f"\n--- Step 3: Phase 1 (TCA-0 Magnus) ---")
        sys.stdout.flush()

    t0 = time.time()

    # Build A_tca matrices vectorized: all k at once per tau
    A_tca_all = np.zeros((N_k, N_tau_tca, nvar_tca, nvar_tca), dtype=np.float32)
    for it in range(N_tau_tca):
        R_val = max(float(R_tca[it]), 1e-10)
        A_tca_all[:, it] = build_A_tca_vectorized(
            k_arr, float(calH_tca[it]), float(a_tca[it]),
            R_val, float(kd_tca[it]), ln_max
        )

    t_build_tca = time.time() - t0

    # Initial conditions (TCA)
    y0_tca = np.zeros((N_k, nvar_tca), dtype=np.float64)
    for ik in range(N_k):
        y0_tca[ik] = tca_ic(k_arr[ik], tau_tca[0], ln_max)

    # GPU propagation
    t0_gpu = time.time()
    A_tca_mx = mx.array(A_tca_all)
    y_tca_mx = mx.array(y0_tca.astype(np.float32))

    y_tca_snapshots = [y_tca_mx]

    for i in range(N_tau_tca - 1):
        dt = float(tau_tca[i + 1] - tau_tca[i])
        A1 = A_tca_mx[:, i]
        A2 = A_tca_mx[:, i + 1]

        Omega = (dt / 2.0) * (A1 + A2)
        if order >= 2:
            comm = mx.matmul(A1, A2) - mx.matmul(A2, A1)
            Omega = Omega + (dt * dt / 12.0) * comm

        expOmega = batched_matrix_exp(Omega, n_terms=n_terms)
        y_tca_mx = mx.einsum('kji,ki->kj', expOmega, y_tca_mx)
        y_tca_snapshots.append(y_tca_mx)

    mx.eval(y_tca_mx)
    t_tca_gpu = time.time() - t0_gpu

    if verbose:
        mem_tca = A_tca_all.nbytes / 1e6
        print(f"A matrix build: {t_build_tca:.2f}s ({mem_tca:.0f} MB)")
        print(f"TCA GPU propagation: {t_tca_gpu:.4f}s")

        # Check for NaN/Inf
        y_tca_np = np.array(y_tca_mx)
        n_bad = np.sum(~np.isfinite(y_tca_np))
        if n_bad > 0:
            print(f"  WARNING: {n_bad} non-finite values in TCA output")
        else:
            print(f"  TCA output: all finite, max|y|={np.max(np.abs(y_tca_np)):.4e}")

    # ==================================================================
    # Step 4: Promote TCA -> full hierarchy at tau_switch
    # ==================================================================
    if verbose:
        print(f"\n--- Step 4: TCA -> Full promotion ---")
        sys.stdout.flush()

    y_tca_final = np.array(y_tca_mx).astype(np.float64)
    y_full_init = np.zeros((N_k, nvar_full), dtype=np.float64)

    calH_sw = float(bg.calH_at_tau(np.array([tau_switch]))[0])
    a_sw = float(bg.a_at_tau(np.array([tau_switch]))[0])
    R_sw = max(float(bg.R_at_tau(np.array([tau_switch]))[0]), 1e-10)
    kd_sw = float(bg.kappa_dot_at_tau(np.array([tau_switch]))[0])

    for ik in range(N_k):
        y_full_init[ik] = tca_to_full(
            y_tca_final[ik], k_arr[ik], calH_sw, a_sw, R_sw, kd_sw,
            lg_max, ln_max
        )

    if verbose:
        print(f"Promoted {N_k} modes: TCA({nvar_tca}) -> Full({nvar_full})")

    # ==================================================================
    # Step 5: Phase 2 -- Full Magnus propagation
    # ==================================================================
    if verbose:
        print(f"\n--- Step 5: Phase 2 (Full Magnus) ---")
        sys.stdout.flush()

    t0 = time.time()

    # Build full A matrices vectorized: all k at once per tau
    A_full_all = np.zeros((N_k, N_tau_full, nvar_full, nvar_full), dtype=np.float32)
    for it in range(N_tau_full):
        R_val = max(float(R_full[it]), 1e-10)
        A_full_all[:, it] = build_A_matrices_vectorized(
            k_arr, float(calH_full[it]), float(a_full[it]),
            R_val, float(kd_full[it]), lg_max, ln_max
        )

    t_build_full = time.time() - t0

    # GPU propagation with adaptive order:
    # Use order=1 when |kappa_dot| is still significant (Thomson damping active).
    # The 2nd-order commutator correction is harmful when A changes rapidly.
    # Use order=2 only after |kd| drops below 0.01 Mpc^-1.
    t0_gpu = time.time()
    A_full_mx = mx.array(A_full_all)
    y_full_mx = mx.array(y_full_init.astype(np.float32))

    y_full_snapshots = [y_full_mx]

    for i in range(N_tau_full - 1):
        dt = float(tau_full[i + 1] - tau_full[i])
        A1 = A_full_mx[:, i]
        A2 = A_full_mx[:, i + 1]

        # 1st order: midpoint trapezoidal rule
        Omega = (dt / 2.0) * (A1 + A2)

        # 2nd order commutator: only when Thomson damping is negligible
        # and the step size is moderate (dt * ||A|| < 10)
        kd_mid = 0.5 * (abs(float(kd_full[i])) + abs(float(kd_full[min(i+1, N_tau_full-1)])))
        if order >= 2 and kd_mid < 0.01 and dt < 100.0:
            comm = mx.matmul(A1, A2) - mx.matmul(A2, A1)
            Omega = Omega + (dt * dt / 12.0) * comm

        expOmega = batched_matrix_exp(Omega, n_terms=n_terms)
        y_full_mx = mx.einsum('kji,ki->kj', expOmega, y_full_mx)
        y_full_snapshots.append(y_full_mx)

    mx.eval(y_full_mx)
    t_full_gpu = time.time() - t0_gpu

    if verbose:
        mem_full = A_full_all.nbytes / 1e6
        print(f"A matrix build: {t_build_full:.2f}s ({mem_full:.0f} MB)")
        print(f"Full GPU propagation: {t_full_gpu:.4f}s")

        y_full_np = np.array(y_full_mx)
        n_bad = np.sum(~np.isfinite(y_full_np))
        if n_bad > 0:
            print(f"  WARNING: {n_bad} non-finite values in full output")
        else:
            print(f"  Full output: all finite, max|y|={np.max(np.abs(y_full_np)):.4e}")

    # ==================================================================
    # Step 6: Gauge transform at all snapshots
    # ==================================================================
    if verbose:
        print(f"\n--- Step 6: Gauge transform ---")
        sys.stdout.flush()

    t0 = time.time()

    # Build LOS snapshot grid from the combined tau grid
    # We need the full-hierarchy state at each snapshot
    # For TCA snapshots, promote TCA state to full hierarchy
    # For full snapshots, use directly

    # Identify which snapshots are TCA vs full
    tau_snap = tau_all
    N_snap = len(tau_snap)
    a_snap = bg.a_at_tau(tau_snap)
    calH_snap = bg.calH_at_tau(tau_snap)
    R_snap = bg.R_at_tau(tau_snap)
    kd_snap = bg.kappa_dot_at_tau(tau_snap)

    Phi_N = np.zeros((N_k, N_snap), dtype=np.float64)
    Psi_N = np.zeros((N_k, N_snap), dtype=np.float64)
    Theta0_N = np.zeros((N_k, N_snap), dtype=np.float64)
    vb_N = np.zeros((N_k, N_snap), dtype=np.float64)
    Theta2_arr = np.zeros((N_k, N_snap), dtype=np.float64)

    fn_s = idx_fn_start(lg_max)

    # Process TCA snapshots
    for snap_idx in range(N_snap):
        tau_val = tau_snap[snap_idx]
        a_val = float(a_snap[snap_idx])
        calH_val = float(calH_snap[snap_idx])
        R_val = max(float(R_snap[snap_idx]), 1e-10)
        kd_val = float(kd_snap[snap_idx])

        if tau_val <= tau_switch:
            # Find corresponding TCA snapshot index
            tca_idx = np.searchsorted(tau_tca, tau_val)
            tca_idx = min(tca_idx, len(y_tca_snapshots) - 1)
            y_tca_snap = np.array(y_tca_snapshots[tca_idx]).astype(np.float64)

            for ik in range(N_k):
                # Promote to full hierarchy for gauge transform
                y_f = tca_to_full(y_tca_snap[ik], k_arr[ik], calH_val,
                                   a_val, R_val, kd_val, lg_max, ln_max)

                h_prime, eta_prime = diagnose_metric(
                    y_f, k_arr[ik], calH_val, a_val, lg_max, ln_max)
                gt = gauge_transform(
                    y_f, k_arr[ik], calH_val, a_val, lg_max, ln_max,
                    h_prime, eta_prime)

                Phi_N[ik, snap_idx] = gt['Phi_N']
                Psi_N[ik, snap_idx] = gt['Psi_N']
                eta_val = y_f[IDX_ETA]
                Theta0_N[ik, snap_idx] = y_f[IDX_FG_START] / 4.0 - eta_val + gt['Phi_N']
                vb_N[ik, snap_idx] = (y_f[IDX_THETA_B] / k_arr[ik]
                                       + (h_prime + 6 * eta_prime) / (2 * k_arr[ik]))
                Theta2_arr[ik, snap_idx] = y_f[IDX_FG_START + 2] / 4.0
        else:
            # Find corresponding full snapshot index
            full_idx = np.searchsorted(tau_full, tau_val)
            full_idx = min(full_idx, len(y_full_snapshots) - 1)
            y_full_snap = np.array(y_full_snapshots[full_idx]).astype(np.float64)

            for ik in range(N_k):
                y = y_full_snap[ik]

                h_prime, eta_prime = diagnose_metric(
                    y, k_arr[ik], calH_val, a_val, lg_max, ln_max)
                gt = gauge_transform(
                    y, k_arr[ik], calH_val, a_val, lg_max, ln_max,
                    h_prime, eta_prime)

                Phi_N[ik, snap_idx] = gt['Phi_N']
                Psi_N[ik, snap_idx] = gt['Psi_N']
                eta_val = y[IDX_ETA]
                Theta0_N[ik, snap_idx] = y[IDX_FG_START] / 4.0 - eta_val + gt['Phi_N']
                vb_N[ik, snap_idx] = (y[IDX_THETA_B] / k_arr[ik]
                                       + (h_prime + 6 * eta_prime) / (2 * k_arr[ik]))
                Theta2_arr[ik, snap_idx] = y[IDX_FG_START + 2] / 4.0

    t_gauge = time.time() - t0

    # Check for NaN in gauge-transformed quantities
    n_nan = np.sum(~np.isfinite(Phi_N)) + np.sum(~np.isfinite(Psi_N))
    if verbose:
        if n_nan > 0:
            print(f"  WARNING: {n_nan} NaN/Inf in gauge-transformed output")
            # Replace NaN with 0 to allow LOS integration to proceed
            Phi_N = np.nan_to_num(Phi_N, nan=0.0, posinf=0.0, neginf=0.0)
            Psi_N = np.nan_to_num(Psi_N, nan=0.0, posinf=0.0, neginf=0.0)
            Theta0_N = np.nan_to_num(Theta0_N, nan=0.0, posinf=0.0, neginf=0.0)
            vb_N = np.nan_to_num(vb_N, nan=0.0, posinf=0.0, neginf=0.0)
            Theta2_arr = np.nan_to_num(Theta2_arr, nan=0.0, posinf=0.0, neginf=0.0)
        else:
            print(f"  All gauge-transformed quantities finite")
        print(f"Gauge transform: {t_gauge:.2f}s")

    # ==================================================================
    # Step 7: Phi'+Psi' by cubic spline
    # ==================================================================
    from scipy.interpolate import CubicSpline

    PhiPsi = Phi_N + Psi_N
    PhiPsi_prime = np.zeros_like(PhiPsi)
    for ik in range(N_k):
        if np.all(np.isfinite(PhiPsi[ik])):
            cs = CubicSpline(tau_snap, PhiPsi[ik])
            PhiPsi_prime[ik] = cs(tau_snap, 1)

    # ==================================================================
    # Step 8: LOS integration -> C_l
    # ==================================================================
    if verbose:
        print(f"\n--- Step 7: LOS integration ---")
        sys.stdout.flush()

    t0 = time.time()

    ell_values = np.unique(np.concatenate([
        np.arange(2, 30, 1),
        np.arange(30, 100, 2),
        np.arange(100, 500, 4),
        np.arange(500, 1500, 8),
        np.arange(1500, 2501, 12),
    ])).astype(int)
    N_ell = len(ell_values)

    P_R = A_s * (k_arr / k_pivot) ** (n_s - 1.0)
    g_snap = bg.visibility_at_tau(tau_snap)
    kappa_snap_vals = bg.kappa_at_tau(tau_snap)
    exp_neg_kappa = np.exp(-kappa_snap_vals)
    chi_snap = bg.tau_0 - tau_snap
    dtau_snap = np.diff(tau_snap)

    # Bessel table
    try:
        from .bessel_cache import CachedBesselTable
        x_max_bessel = float(np.max(k_arr) * np.max(chi_snap)) * 1.05 + 50.0
        bessel_table = CachedBesselTable(ell_values, x_max=x_max_bessel,
                                          N_cheb=64, seg_width=80,
                                          verbose=verbose)
        use_gpu_bessel = True
    except Exception as e:
        if verbose:
            print(f"[WARNING] Bessel table failed ({e}), using scipy")
        use_gpu_bessel = False

    Delta_l = np.zeros((N_ell, N_k), dtype=np.float64)
    Delta_l_E = np.zeros((N_ell, N_k), dtype=np.float64)

    alpha_P = 1.7
    pol_prefactor = 0.75 * alpha_P

    eps_prefactors = np.zeros(N_ell)
    for il, ell in enumerate(ell_values):
        l_int = int(ell)
        if l_int >= 2:
            eps_prefactors[il] = np.sqrt(
                float((l_int - 1) * l_int * (l_int + 1) * (l_int + 2)))

    if verbose:
        print(f"Computing LOS transfer functions (TT + EE)...")
        print(f"  {N_ell} ell values, {N_snap} snapshots, {N_k} k-modes")
        sys.stdout.flush()

    for it in range(N_snap):
        chi = chi_snap[it]
        if chi <= 0:
            continue

        if it == 0:
            w = 0.5 * dtau_snap[0]
        elif it == N_snap - 1:
            w = 0.5 * dtau_snap[-1]
        else:
            w = 0.5 * (dtau_snap[max(0, it - 1)] +
                        dtau_snap[min(it, len(dtau_snap) - 1)])

        S_SW = g_snap[it] * (Theta0_N[:, it] + Psi_N[:, it])
        S_Dop = g_snap[it] * vb_N[:, it]
        S_ISW = exp_neg_kappa[it] * PhiPsi_prime[:, it]
        S_E = g_snap[it] * pol_prefactor * Theta2_arr[:, it]

        if use_gpu_bessel:
            x_arr = k_arr * chi
            jl_mx, jlp_mx = bessel_table.eval_jl_jlp(x_arr.astype(np.float32))
            mx.eval(jl_mx, jlp_mx)
            jl = np.array(jl_mx)
            jlp = np.array(jlp_mx)

            x_safe = np.where(np.abs(x_arr) < 1e-10, 1e-10, x_arr)
            inv_x2 = 1.0 / (x_safe ** 2)

            for il in range(N_ell):
                Delta_l[il, :] += w * (
                    S_SW * jl[il, :] +
                    S_Dop * jlp[il, :] +
                    S_ISW * jl[il, :]
                )
                if eps_prefactors[il] > 0:
                    eps_l = eps_prefactors[il] * jl[il, :] * inv_x2
                    eps_l = np.where(np.abs(x_arr) < 1e-10, 0.0, eps_l)
                    Delta_l_E[il, :] += w * S_E * eps_l
        else:
            from scipy.special import spherical_jn
            x = k_arr * chi
            x_safe = np.where(np.abs(x) < 1e-10, 1e-10, x)
            inv_x2 = 1.0 / (x_safe ** 2)
            for il, ell in enumerate(ell_values):
                jl = spherical_jn(int(ell), x)
                jlp = spherical_jn(int(ell), x, derivative=True)
                Delta_l[il, :] += w * (S_SW * jl + S_Dop * jlp + S_ISW * jl)
                if eps_prefactors[il] > 0:
                    eps_l = eps_prefactors[il] * jl * inv_x2
                    eps_l = np.where(np.abs(x) < 1e-10, 0.0, eps_l)
                    Delta_l_E[il, :] += w * S_E * eps_l

    # No Silk damping correction needed: the full photon hierarchy in
    # Phase 2 naturally captures Thomson scattering and diffusion damping.

    t_los = time.time() - t0
    if verbose:
        print(f"LOS integration: {t_los:.2f}s")

    # ==================================================================
    # Step 9: C_l
    # ==================================================================
    t0 = time.time()

    lnk = np.log(k_arr)
    dlnk = np.diff(lnk)
    ell_f = ell_values.astype(float)
    uK2 = (T_CMB * 1e6) ** 2
    norm = (2.0 / 3.0) ** 2

    # TT
    integrand_TT = P_R[None, :] * Delta_l ** 2
    mid_TT = 0.5 * (integrand_TT[:, :-1] + integrand_TT[:, 1:])
    Cl_TT = 4.0 * np.pi * np.sum(mid_TT * dlnk[None, :], axis=1)
    Cl_TT = np.maximum(Cl_TT, 0.0)
    Dl_TT = norm * ell_f * (ell_f + 1.0) * Cl_TT / (2.0 * np.pi) * uK2

    # EE
    integrand_EE = P_R[None, :] * Delta_l_E ** 2
    mid_EE = 0.5 * (integrand_EE[:, :-1] + integrand_EE[:, 1:])
    Cl_EE = 4.0 * np.pi * np.sum(mid_EE * dlnk[None, :], axis=1)
    Cl_EE = np.maximum(Cl_EE, 0.0)
    Dl_EE = norm * ell_f * (ell_f + 1.0) * Cl_EE / (2.0 * np.pi) * uK2

    # TE
    integrand_TE = P_R[None, :] * Delta_l * Delta_l_E
    mid_TE = 0.5 * (integrand_TE[:, :-1] + integrand_TE[:, 1:])
    Cl_TE = 4.0 * np.pi * np.sum(mid_TE * dlnk[None, :], axis=1)
    Dl_TE = norm * ell_f * (ell_f + 1.0) * Cl_TE / (2.0 * np.pi) * uK2

    t_cl = time.time() - t0
    t_elapsed = time.time() - t_total

    t_magnus_total = t_tca_gpu + t_full_gpu

    if verbose:
        print(f"\n{'='*70}")
        print(f"TIMING SUMMARY")
        print(f"{'='*70}")
        print(f"  Background:          {t_bg:.2f}s")
        print(f"  Setup:               {t_setup:.2f}s")
        print(f"  A build (TCA):       {t_build_tca:.2f}s")
        print(f"  A build (Full):      {t_build_full:.2f}s")
        print(f"  Magnus GPU (TCA):    {t_tca_gpu:.4f}s")
        print(f"  Magnus GPU (Full):   {t_full_gpu:.4f}s")
        print(f"  Magnus GPU TOTAL:    {t_magnus_total:.4f}s  <-- ODE replacement")
        print(f"  Gauge transform:     {t_gauge:.2f}s")
        print(f"  LOS integration:     {t_los:.2f}s")
        print(f"  C_l sum:             {t_cl:.3f}s")
        print(f"  TOTAL:               {t_elapsed:.2f}s")
        print(f"{'='*70}")

        # Report peak heights
        idx_220 = np.argmin(np.abs(ell_values - 220))
        idx_550 = np.argmin(np.abs(ell_values - 540))
        idx_800 = np.argmin(np.abs(ell_values - 810))
        print(f"\nPeak heights:")
        print(f"  Dl(l~220) = {Dl_TT[idx_220]:.0f} uK^2  (CLASS: ~5700)")
        print(f"  Dl(l~540) = {Dl_TT[idx_550]:.0f} uK^2  (CLASS: ~3000)")
        print(f"  Dl(l~810) = {Dl_TT[idx_800]:.0f} uK^2  (CLASS: ~2500)")

    return {
        'ell': ell_values,
        'Cl_TT': Cl_TT,
        'Dl_TT': Dl_TT,
        'Cl_EE': Cl_EE,
        'Dl_EE': Dl_EE,
        'Cl_TE': Cl_TE,
        'Dl_TE': Dl_TE,
        'k_arr': k_arr,
        'bg': bg,
        'tau_grid': tau_all,
        'tau_tca': tau_tca,
        'tau_full': tau_full,
        'tau_switch': tau_switch,
        'Delta_l': Delta_l,
        'Delta_l_E': Delta_l_E,
        'Phi_N': Phi_N,
        'Psi_N': Psi_N,
        'Theta0_N': Theta0_N,
        'vb_N': vb_N,
        'Theta2_arr': Theta2_arr,
        'timing': {
            'background': t_bg,
            'setup': t_setup,
            'build_tca': t_build_tca,
            'build_full': t_build_full,
            'magnus_tca': t_tca_gpu,
            'magnus_full': t_full_gpu,
            'magnus_total': t_magnus_total,
            'gauge': t_gauge,
            'los': t_los,
            'total': t_elapsed,
        },
    }


# ============================================================================
# Compare with reference
# ============================================================================

def compare_with_reference(result, class_file=None, verbose=True):
    """Compare Magnus solver with CLASS reference data."""
    import os
    base_dir = os.path.dirname(__file__)

    ref_files = [
        os.path.join(base_dir, 'class_comparison.dat'),
        os.path.join(base_dir, 'cl_sync_gauge.dat'),
        os.path.join(base_dir, 'cl_lcdm.dat'),
    ]
    if class_file:
        ref_files.insert(0, class_file)

    ref_data = None
    for f in ref_files:
        if os.path.exists(f):
            try:
                ref_data = np.loadtxt(f)
                if verbose:
                    print(f"Reference: {f}")
                break
            except Exception:
                continue

    if ref_data is None:
        if verbose:
            print("No reference data found.")
        return None

    ell_ref = ref_data[:, 0]
    Dl_ref = ref_data[:, 1]

    from scipy.interpolate import interp1d
    ell_us = result['ell'].astype(float)
    Dl_us = result['Dl_TT']

    mask = (ell_ref >= ell_us[0]) & (ell_ref <= ell_us[-1])
    ell_cmp = ell_ref[mask]
    Dl_ref_cmp = Dl_ref[mask]

    f_interp = interp1d(ell_us, Dl_us, kind='cubic')
    Dl_us_cmp = f_interp(ell_cmp)

    mask2 = ell_cmp > 10
    frac_err = np.abs(Dl_us_cmp[mask2] - Dl_ref_cmp[mask2]) / (np.abs(Dl_ref_cmp[mask2]) + 1e-30)
    rms = np.sqrt(np.mean(frac_err ** 2))
    max_err = np.max(frac_err)

    if verbose:
        print(f"\nRMS fractional error (ell > 10): {rms:.4f} ({rms*100:.1f}%)")
        print(f"Max fractional error: {max_err:.4f} ({max_err*100:.1f}%)")

    return {
        'rms': rms,
        'max_err': max_err,
        'ell_cmp': ell_cmp,
        'Dl_us': Dl_us_cmp,
        'Dl_ref': Dl_ref_cmp,
    }


# ============================================================================
# Convergence test
# ============================================================================

def convergence_test(N_k=200, steps_list=None, verbose=True):
    """Test convergence with increasing N_steps."""
    if steps_list is None:
        steps_list = [10, 20, 30, 50, 80]

    results = {}
    for N_steps in steps_list:
        if verbose:
            print(f"\n{'#'*70}")
            print(f"# N_steps_full = {N_steps}")
            print(f"{'#'*70}")
        res = solve_magnus(N_k=N_k, N_steps_full=N_steps, verbose=verbose)
        results[N_steps] = res

    if verbose:
        print(f"\n{'='*70}")
        print("CONVERGENCE SUMMARY")
        print(f"{'='*70}")
        print(f"{'N_steps':>8} {'Time(s)':>8} {'Dl(220)':>10} {'Dl(1000)':>10}")
        print("-" * 40)

        for N_steps in steps_list:
            res = results[N_steps]
            ell = res['ell']
            Dl = res['Dl_TT']
            idx_220 = np.argmin(np.abs(ell - 220))
            idx_1000 = np.argmin(np.abs(ell - 1000))
            t = res['timing']['total']
            print(f"{N_steps:8d} {t:8.2f} {Dl[idx_220]:10.1f} {Dl[idx_1000]:10.1f}")

    return results


# ============================================================================
# Main
# ============================================================================

def main():
    """Run the Magnus solver."""
    import argparse
    parser = argparse.ArgumentParser(description='Magnus expansion Boltzmann solver')
    parser.add_argument('--N_k', type=int, default=500, help='Number of k-modes')
    parser.add_argument('--N_steps', type=int, default=300, help='Full-phase Magnus steps')
    parser.add_argument('--N_steps_tca', type=int, default=20, help='TCA-phase Magnus steps')
    parser.add_argument('--order', type=int, default=2, choices=[1, 2],
                        help='Magnus expansion order')
    parser.add_argument('--convergence', action='store_true',
                        help='Run convergence test')
    parser.add_argument('--plot', action='store_true', help='Plot results')
    args = parser.parse_args()

    if args.convergence:
        results = convergence_test(N_k=200, verbose=True)
        if args.plot:
            _plot_convergence(results)
        return

    result = solve_magnus(
        N_k=args.N_k,
        N_steps_tca=args.N_steps_tca,
        N_steps_full=args.N_steps,
        order=args.order,
        verbose=True,
    )

    cmp = compare_with_reference(result, verbose=True)

    outfile = os.path.join(os.path.dirname(__file__), 'cl_magnus.dat')
    np.savetxt(outfile,
               np.column_stack([result['ell'], result['Dl_TT'],
                                result['Dl_EE'], result['Dl_TE']]),
               header='ell  Dl_TT  Dl_EE  Dl_TE',
               fmt='%8d %15.6e %15.6e %15.6e')
    print(f"\nResults saved to {outfile}")

    if args.plot:
        _plot_result(result, cmp)


def _plot_result(result, cmp=None):
    """Plot D_l^TT and comparison."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 1, figsize=(12, 8),
                              gridspec_kw={'height_ratios': [3, 1]})

    ell = result['ell']
    Dl = result['Dl_TT']

    axes[0].plot(ell, Dl, 'b-', lw=1.5, label='Magnus solver')
    if cmp is not None:
        axes[0].plot(cmp['ell_cmp'], cmp['Dl_ref'], 'r--', lw=1,
                     label=f"Reference (RMS={cmp['rms']*100:.1f}%)")
    axes[0].set_ylabel(r'$D_\ell^{TT}$ [$\mu K^2$]')
    axes[0].set_xlim(2, 2500)
    axes[0].set_ylim(0, None)
    axes[0].legend()
    axes[0].set_title('Magnus Expansion Boltzmann Solver (Two-Phase)')

    if cmp is not None:
        frac = (cmp['Dl_us'] - cmp['Dl_ref']) / (cmp['Dl_ref'] + 1e-30)
        axes[1].plot(cmp['ell_cmp'], frac * 100, 'k-', lw=0.5)
        axes[1].axhline(0, color='gray', ls='--')
        axes[1].set_ylabel('Residual [%]')
        axes[1].set_xlim(2, 2500)
        axes[1].set_ylim(-50, 50)

    axes[1].set_xlabel(r'Multipole $\ell$')
    plt.tight_layout()
    outfile = os.path.join(os.path.dirname(__file__), 'cl_magnus.png')
    plt.savefig(outfile, dpi=150)
    print(f"Plot saved to {outfile}")
    plt.close()


def _plot_convergence(results):
    """Plot convergence test."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(10, 6))
    for N_steps, res in sorted(results.items()):
        ell = res['ell']
        Dl = res['Dl_TT']
        t = res['timing']['total']
        ax.plot(ell, Dl, label=f'N={N_steps} ({t:.1f}s)')

    ax.set_xlabel(r'$\ell$')
    ax.set_ylabel(r'$D_\ell^{TT}$ [$\mu K^2$]')
    ax.set_xlim(2, 2500)
    ax.set_ylim(0, None)
    ax.legend(fontsize=8)
    ax.set_title('Magnus Solver Convergence')
    plt.tight_layout()
    outfile = os.path.join(os.path.dirname(__file__), 'cl_magnus_convergence.png')
    plt.savefig(outfile, dpi=150)
    print(f"Plot saved to {outfile}")
    plt.close()


if __name__ == '__main__':
    main()
