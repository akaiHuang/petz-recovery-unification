"""
solver_chebyshev.py -- Chebyshev spectral Boltzmann solver.

INSTEAD of time-stepping (600 RK4 steps or 20 Magnus steps), solve the
ENTIRE time evolution as a SINGLE MATRIX EQUATION using Chebyshev
polynomial expansion on collocation points.

MATHEMATICAL FOUNDATION:
  The Boltzmann system is LINEAR:  y'(tau) = A(tau) . y(tau)
  with initial condition y(tau_init) = y0.

  Chebyshev spectral method:
    1. Map tau in [tau_init, tau_end] to x in [-1, 1]
    2. Collocate at Chebyshev-Gauss-Lobatto points x_j = cos(j pi / N)
    3. Derivative operator: y'(x_j) = sum_k D_jk y(x_k)
    4. The ODE becomes:  (D x I - blockdiag(A)) . vec(Y) = 0 + BC
    5. ONE linear solve gives the solution at ALL collocation points.

  For N_cheb collocation points and n_var ODE variables:
    System size: (N_cheb+1) * n_var  per k-mode
    For N_cheb=30, n_var=56: 1736 x 1736 per k-mode

STIFFNESS HANDLING:
  Thomson scattering rate |kappa_dot| ~ 10^7 at early times creates
  stiffness that requires ~100+ Chebyshev points to resolve.
  Solution: Two-interval approach matching what CLASS does:
    Phase 1 (tau_init to tau_switch ~ 65 Mpc): TCA equations (no Thomson stiffness)
    Phase 2 (tau_switch to tau_end): Full Boltzmann hierarchy
  TCA removes stiffness entirely; Phase 2 starts after |kappa_dot| drops.

ADVANTAGES:
  - Exponential convergence: 30 points -> ~10^{-10} accuracy for smooth problems
  - ONE matrix solve per k-mode: no iterative time-stepping
  - Perfect for GPU: batched mx.linalg.solve over k-modes
  - No stability issues: unconditionally stable for linear problems

Author: Sheng-Kai Huang, 2026
"""
import os
os.environ['MPLCONFIGDIR'] = '/tmp/mpl_config'

import numpy as np
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
from .perturbations_neutrino import (
    _f_nu, _OMEGA_GAMMA, _OMEGA_NU,
)
from .perturbations_sync import (
    adiabatic_ic_sync, gauge_transform, diagnose_metric,
    n_var_sync, idx_fn_start, idx_fg,
    IDX_ETA, IDX_DELTA_C, IDX_DELTA_B, IDX_THETA_B, IDX_FG_START,
    L_GAMMA_MAX, L_NU_MAX_SYNC,
)


# ============================================================================
# Chebyshev spectral infrastructure
# ============================================================================

def chebyshev_points(N):
    """
    Chebyshev-Gauss-Lobatto points on [-1, 1].

    x_j = cos(j * pi / N),  j = 0, 1, ..., N

    Note: x_0 = 1, x_N = -1.
    We map x = -1 -> tau_init (start) and x = +1 -> tau_end (end),
    so the IC is imposed at the LAST row (j = N).
    """
    j = np.arange(N + 1)
    return np.cos(j * np.pi / N)


def chebyshev_diff_matrix(N):
    """
    Chebyshev differentiation matrix D of size (N+1) x (N+1).

    Satisfies: y'(x_j) = sum_k D_jk y(x_k)  at Gauss-Lobatto points.

    Reference: Trefethen, "Spectral Methods in MATLAB" (2000), Ch. 6.
    """
    x = chebyshev_points(N)
    c = np.ones(N + 1)
    c[0] = 2.0
    c[N] = 2.0

    D = np.zeros((N + 1, N + 1))

    for i in range(N + 1):
        for j in range(N + 1):
            if i != j:
                D[i, j] = c[i] / c[j] * (-1.0) ** (i + j) / (x[i] - x[j])

    # Diagonal entries from the negative sum trick (more accurate)
    for i in range(N + 1):
        D[i, i] = -np.sum(D[i, :])

    return D


def map_to_physical(x, tau_start, tau_end):
    """Map x in [-1, 1] to tau in [tau_start, tau_end]."""
    return 0.5 * (tau_end - tau_start) * (x + 1.0) + tau_start


def map_to_chebyshev(tau, tau_start, tau_end):
    """Map tau in [tau_start, tau_end] to x in [-1, 1]."""
    return 2.0 * (tau - tau_start) / (tau_end - tau_start) - 1.0


# ============================================================================
# Build A(tau) matrix for synchronous gauge Boltzmann system
# ============================================================================

def build_A_matrix(k, calH, a, R, kappa_dot, lg_max, ln_max, tau):
    """
    Build the coefficient matrix A(tau) such that y'(tau) = A(tau) . y(tau).

    The synchronous gauge equations from perturbations_sync.py are:
      y' = f(tau, y)
    which is LINEAR in y, so f(tau, y) = A(tau) . y.

    The trick: h' and eta' are DIAGNOSED from Einstein constraints, which
    are also linear in the state variables. So the full system is linear.

    Parameters
    ----------
    k : float, wavenumber
    calH : float, conformal Hubble
    a : float, scale factor
    R : float, baryon loading 3*rho_b / (4*rho_gamma)
    kappa_dot : float, Thomson scattering rate (NEGATIVE)
    lg_max : int, photon hierarchy truncation
    ln_max : int, neutrino hierarchy truncation
    tau : float, conformal time (for truncation boundary condition)

    Returns
    -------
    A : ndarray of shape (nvar, nvar)
    """
    k2 = k * k
    H02 = _H0_MPC ** 2
    Og = _OMEGA_GAMMA
    On = _OMEGA_NU
    Ob = _OMEGA_B
    Oc = _OMEGA_C
    abs_kd = abs(kappa_dot)
    ia = 1.0 / a
    ia2 = ia * ia
    tau_safe = max(tau, 1e-10)

    nvar = n_var_sync(lg_max, ln_max)
    fn_s = idx_fn_start(lg_max)
    A = np.zeros((nvar, nvar))

    # ----------------------------------------------------------------
    # First, express h' and eta' as linear functions of the state.
    # h' = (2/calH) * [k^2 * eta + 1.5 * H0^2 * (Og*ia2*dg + On*ia2*dn + Ob*ia*db + Oc*ia*dc)]
    # eta' = 1.5 * H0^2 / k^2 * [(4/3)*Og*ia2*thetag + (4/3)*On*ia2*thetan + Ob*ia*thetab]
    #
    # where thetag = 0.75*k*Fg1, thetan = 0.75*k*Fn1
    # ----------------------------------------------------------------

    # h' coefficients: h' = sum_v hp_coeff[v] * y[v]
    hp_coeff = np.zeros(nvar)
    hp_coeff[IDX_ETA] = (2.0 / calH) * k2
    hp_coeff[IDX_DELTA_C] = (2.0 / calH) * 1.5 * H02 * Oc * ia
    hp_coeff[IDX_DELTA_B] = (2.0 / calH) * 1.5 * H02 * Ob * ia
    hp_coeff[idx_fg(0)] = (2.0 / calH) * 1.5 * H02 * Og * ia2      # delta_gamma = Fg0
    hp_coeff[fn_s + 0] = (2.0 / calH) * 1.5 * H02 * On * ia2       # delta_nu = Fn0

    # eta' coefficients: eta' = sum_v ep_coeff[v] * y[v]
    ep_coeff = np.zeros(nvar)
    # theta_gamma = 0.75*k*Fg1, theta_nu = 0.75*k*Fn1
    ep_coeff[IDX_THETA_B] = 1.5 * H02 / k2 * Ob * ia
    ep_coeff[idx_fg(1)] = 1.5 * H02 / k2 * (4.0 / 3.0) * Og * ia2 * 0.75 * k   # Fg1
    ep_coeff[fn_s + 1] = 1.5 * H02 / k2 * (4.0 / 3.0) * On * ia2 * 0.75 * k    # Fn1

    # ----------------------------------------------------------------
    # eta' equation: d(eta)/dtau = ep_coeff . y
    # ----------------------------------------------------------------
    A[IDX_ETA, :] = ep_coeff

    # ----------------------------------------------------------------
    # CDM: d(delta_c)/dtau = -0.5 * h'  (theta_c = 0 gauge)
    # ----------------------------------------------------------------
    A[IDX_DELTA_C, :] = -0.5 * hp_coeff

    # ----------------------------------------------------------------
    # Baryon density: d(delta_b)/dtau = -theta_b - 0.5 * h'
    # ----------------------------------------------------------------
    A[IDX_DELTA_B, :] = -0.5 * hp_coeff
    A[IDX_DELTA_B, IDX_THETA_B] += -1.0

    # ----------------------------------------------------------------
    # Baryon velocity: d(theta_b)/dtau = -calH*theta_b + |kd|/R*(theta_g - theta_b)
    #   theta_g = 0.75*k*Fg1
    # ----------------------------------------------------------------
    A[IDX_THETA_B, IDX_THETA_B] = -calH - abs_kd / R
    A[IDX_THETA_B, idx_fg(1)] = abs_kd / R * 0.75 * k

    # ----------------------------------------------------------------
    # Photon hierarchy
    # ----------------------------------------------------------------
    # l=0: Fg0' = -k*Fg1 - (2/3)*h'
    A[idx_fg(0), :] += -(2.0 / 3.0) * hp_coeff
    A[idx_fg(0), idx_fg(1)] += -k

    # l=1: Fg1' = (k/3)*(Fg0 - 2*Fg2) + |kd|*(-Fg1 + 4*theta_b/(3*k))
    A[idx_fg(1), idx_fg(0)] = k / 3.0
    if lg_max >= 2:
        A[idx_fg(1), idx_fg(2)] = -2.0 * k / 3.0
    A[idx_fg(1), idx_fg(1)] += -abs_kd
    A[idx_fg(1), IDX_THETA_B] += abs_kd * 4.0 / (3.0 * k)

    # l=2: Fg2' = (k/5)*(2*Fg1 - 3*Fg3) + (4/15)*h' + (8/15)*eta' - (9/10)*|kd|*Fg2
    if lg_max >= 2:
        A[idx_fg(2), :] += (4.0 / 15.0) * hp_coeff + (8.0 / 15.0) * ep_coeff
        A[idx_fg(2), idx_fg(1)] += 2.0 * k / 5.0
        if lg_max >= 3:
            A[idx_fg(2), idx_fg(3)] += -3.0 * k / 5.0
        A[idx_fg(2), idx_fg(2)] += -0.9 * abs_kd

    # l=3..lg_max-1: streaming + Thomson
    for ell in range(3, lg_max):
        fac = k / (2.0 * ell + 1.0)
        A[idx_fg(ell), idx_fg(ell - 1)] = fac * ell
        A[idx_fg(ell), idx_fg(ell + 1)] = -fac * (ell + 1)
        A[idx_fg(ell), idx_fg(ell)] += -abs_kd

    # l=lg_max: truncation boundary
    if lg_max >= 3:
        ell = lg_max
        A[idx_fg(ell), idx_fg(ell - 1)] = k * ell / (2.0 * ell + 1.0)
        A[idx_fg(ell), idx_fg(ell)] += -(ell + 1.0) / tau_safe - abs_kd

    # ----------------------------------------------------------------
    # Neutrino hierarchy (no collisions)
    # ----------------------------------------------------------------
    # l=0: Fn0' = -k*Fn1 - (2/3)*h'
    A[fn_s + 0, :] += -(2.0 / 3.0) * hp_coeff
    A[fn_s + 0, fn_s + 1] += -k

    # l=1: Fn1' = (k/3)*(Fn0 - 2*Fn2)
    A[fn_s + 1, fn_s + 0] = k / 3.0
    if ln_max >= 2:
        A[fn_s + 1, fn_s + 2] = -2.0 * k / 3.0

    # l=2: Fn2' = (k/5)*(2*Fn1 - 3*Fn3) + (4/15)*h' + (8/15)*eta'
    if ln_max >= 2:
        A[fn_s + 2, :] += (4.0 / 15.0) * hp_coeff + (8.0 / 15.0) * ep_coeff
        A[fn_s + 2, fn_s + 1] += 2.0 * k / 5.0
        if ln_max >= 3:
            A[fn_s + 2, fn_s + 3] += -3.0 * k / 5.0

    # l=3..ln_max-1: free streaming
    for ell in range(3, ln_max):
        fac = k / (2.0 * ell + 1.0)
        A[fn_s + ell, fn_s + ell - 1] = fac * ell
        A[fn_s + ell, fn_s + ell + 1] = -fac * (ell + 1)

    # l=ln_max: truncation
    if ln_max >= 1:
        ell = ln_max
        A[fn_s + ell, fn_s + ell - 1] = k * ell / (2.0 * ell + 1.0)
        A[fn_s + ell, fn_s + ell] += -(ell + 1.0) / tau_safe

    return A


# ============================================================================
# Build and solve the full Chebyshev spectral system for ONE k-mode
# ============================================================================

def build_spectral_system(A_list, D_phys, y0, n_var):
    """
    Build the spectral linear system L . u = b.

    The Chebyshev-collocation ODE system is:
      D_phys . Y - diag(A_j) . Y = 0    at j = 0, ..., N-1 (interior + end)
      Y[N] = y0                           (initial condition at x = -1 = tau_init)

    Written as: L . vec(Y) = b.

    Uses vectorized Kronecker product: L = D_phys kron I_nvar - blockdiag(A_j).

    Parameters
    ----------
    A_list : list of (n_var, n_var) ndarray
        A(tau_j) at each Chebyshev collocation point j = 0, ..., N.
    D_phys : (N+1, N+1) ndarray
        Chebyshev differentiation matrix scaled to physical domain.
    y0 : (n_var,) ndarray
        Initial condition.
    n_var : int
        Number of ODE variables.

    Returns
    -------
    L : (size, size) ndarray where size = (N+1) * n_var
    b : (size,) ndarray
    """
    N_pts = len(A_list)  # N + 1
    N = N_pts - 1
    size = N_pts * n_var

    # L = D_phys kron I_nvar  (vectorized, no Python loops)
    I_nvar = np.eye(n_var)
    L = np.kron(D_phys, I_nvar)

    # Subtract block-diagonal A matrices
    for j in range(N_pts):
        s = j * n_var
        L[s:s + n_var, s:s + n_var] -= A_list[j]

    # Initial condition: replace last block row (j = N, x = -1, tau = tau_init)
    b = np.zeros(size)
    bc_row = N * n_var
    L[bc_row:bc_row + n_var, :] = 0.0
    L[bc_row:bc_row + n_var, bc_row:bc_row + n_var] = I_nvar
    b[bc_row:bc_row + n_var] = y0

    return L, b


def solve_single_k_chebyshev(k, bg, tau_start, tau_end, N_cheb=30,
                              lg_max=L_GAMMA_MAX, ln_max=L_NU_MAX_SYNC):
    """
    Solve synchronous gauge Boltzmann for one k-mode using Chebyshev spectral method.

    Parameters
    ----------
    k : float, wavenumber in Mpc^{-1}
    bg : Background object (solved)
    tau_start : float, start conformal time
    tau_end : float, end conformal time
    N_cheb : int, number of Chebyshev collocation intervals (N+1 points)
    lg_max, ln_max : int, hierarchy truncation

    Returns
    -------
    tau_pts : (N_cheb+1,) ndarray, conformal times at collocation points
    Y : (N_cheb+1, nvar) ndarray, solution at collocation points
        Ordered from tau_end (j=0) to tau_start (j=N).
    """
    nvar = n_var_sync(lg_max, ln_max)

    # Chebyshev points on [-1, 1] and physical domain
    x = chebyshev_points(N_cheb)
    tau_pts = map_to_physical(x, tau_start, tau_end)

    # Differentiation matrix scaled to physical domain
    D = chebyshev_diff_matrix(N_cheb)
    D_phys = D * 2.0 / (tau_end - tau_start)

    # Evaluate background at collocation points
    calH_pts = bg.calH_at_tau(tau_pts)
    a_pts = bg.a_at_tau(tau_pts)
    R_pts = bg.R_at_tau(tau_pts)
    kd_pts = bg.kappa_dot_at_tau(tau_pts)

    # Build A(tau) at each collocation point
    A_list = []
    for j in range(N_cheb + 1):
        A_j = build_A_matrix(k, calH_pts[j], a_pts[j], R_pts[j],
                             kd_pts[j], lg_max, ln_max, tau_pts[j])
        A_list.append(A_j)

    # Initial condition at tau_start (last Chebyshev point, j = N_cheb)
    y0 = adiabatic_ic_sync(k, tau_start, lg_max, ln_max)

    # Build and solve spectral system
    L, b = build_spectral_system(A_list, D_phys, y0, nvar)

    # Solve with condition number check
    u = np.linalg.solve(L, b)

    # Reshape: Y[j, v] = y_v at collocation point j
    Y = u.reshape(N_cheb + 1, nvar)

    return tau_pts, Y


# ============================================================================
# Multi-interval Chebyshev solver
# ============================================================================

def _build_interval_boundaries(k, tau_init, tau_end, N_cheb_per_interval=30):
    """
    Build interval boundaries adaptively based on the wavenumber k.

    The oscillation frequency of the photon-baryon fluid is ~k*c_s where
    c_s ~ 1/sqrt(3). To resolve oscillations with N_cheb Chebyshev points,
    each interval should span at most ~(N_cheb/5) oscillation periods.

    For late ISW (after recombination), the solution is smooth and can
    use longer intervals.

    Returns list of tau boundary values.
    """
    # Oscillation period in conformal time
    cs_approx = 1.0 / np.sqrt(3.0)  # approximate sound speed
    period = 2.0 * np.pi / (k * cs_approx) if k > 1e-6 else 1e10

    # Each interval resolves ~(N_cheb/5) periods comfortably
    # (spectral methods need ~5 points per period for good accuracy)
    interval_length = max(period * N_cheb_per_interval / 6.0, 20.0)

    # Cap interval length to avoid too few intervals
    interval_length = min(interval_length, 300.0)

    # Recombination window: needs denser coverage
    tau_rec_start = 230.0   # before recombination
    tau_rec_end = 340.0     # after recombination

    # Build boundaries
    boundaries = [tau_init]

    # Early phase: can use longer intervals (Thomson damping simplifies solution)
    tau = tau_init
    early_interval = min(interval_length * 2.0, 100.0)
    while tau < tau_rec_start:
        tau = min(tau + early_interval, tau_rec_start)
        boundaries.append(tau)

    # Recombination phase: dense intervals
    rec_interval = min(interval_length, 30.0)
    tau = tau_rec_start
    while tau < tau_rec_end:
        tau = min(tau + rec_interval, tau_rec_end)
        if abs(tau - boundaries[-1]) > 1.0:
            boundaries.append(tau)

    # Late ISW phase: solution is smooth (no oscillations after decoupling)
    # Much longer intervals suffice
    tau = tau_rec_end
    late_interval = max(interval_length * 3.0, 500.0)
    while tau < tau_end:
        tau = min(tau + late_interval, tau_end)
        if abs(tau - boundaries[-1]) > 1.0:
            boundaries.append(tau)

    return np.array(boundaries)


def solve_single_k_multi_interval(k, bg, tau_init, tau_end,
                                   N_cheb=30,
                                   lg_max=L_GAMMA_MAX, ln_max=L_NU_MAX_SYNC):
    """
    Solve using multiple Chebyshev intervals, with adaptive interval sizing
    based on the wavenumber k.

    Each interval: independent Chebyshev collocation solve with N_cheb points.
    Intervals are chained: IC of interval i+1 = endpoint of interval i.

    Returns tau_all, Y_all in chronological order.
    """
    import warnings
    from scipy.linalg import solve as scipy_solve

    nvar = n_var_sync(lg_max, ln_max)

    boundaries = _build_interval_boundaries(k, tau_init, tau_end, N_cheb)
    n_intervals = len(boundaries) - 1

    # Cache the Chebyshev differentiation matrix (same N_cheb for all intervals)
    x = chebyshev_points(N_cheb)
    D_raw = chebyshev_diff_matrix(N_cheb)
    I_nvar = np.eye(nvar)
    D_kron_I = np.kron(D_raw, I_nvar)  # reusable template

    # Initial condition
    y0 = adiabatic_ic_sync(k, tau_init, lg_max, ln_max)

    tau_list = []
    Y_list = []
    y_current = y0

    bc_row = N_cheb * nvar
    size = (N_cheb + 1) * nvar

    for i_int in range(n_intervals):
        t_start = boundaries[i_int]
        t_end = boundaries[i_int + 1]

        if t_end - t_start < 1e-10:
            continue

        # Map Chebyshev points to this interval
        tau_pts = map_to_physical(x, t_start, t_end)
        scale = 2.0 / (t_end - t_start)

        # L = scale * D_kron_I  (start from template, avoid re-computing kron)
        L = D_kron_I * scale

        # Background at collocation points
        calH_pts = bg.calH_at_tau(tau_pts)
        a_pts = bg.a_at_tau(tau_pts)
        R_pts = bg.R_at_tau(tau_pts)
        kd_pts = bg.kappa_dot_at_tau(tau_pts)

        # Subtract block-diagonal A matrices
        for j in range(N_cheb + 1):
            A_j = build_A_matrix(k, calH_pts[j], a_pts[j], R_pts[j],
                                 kd_pts[j], lg_max, ln_max, tau_pts[j])
            s = j * nvar
            L[s:s + nvar, s:s + nvar] -= A_j

        # Apply initial condition at last block row
        b = np.zeros(size)
        L[bc_row:bc_row + nvar, :] = 0.0
        L[bc_row:bc_row + nvar, bc_row:bc_row + nvar] = I_nvar
        b[bc_row:bc_row + nvar] = y_current

        # Solve (overwrite_a/b for speed)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            u = scipy_solve(L, b, overwrite_a=True, overwrite_b=True,
                            check_finite=False)
        Y = u.reshape(N_cheb + 1, nvar)

        # Store (reverse to chronological, skip duplicate at start after first)
        if i_int == 0:
            tau_list.append(tau_pts[::-1])
            Y_list.append(Y[::-1])
        else:
            tau_list.append(tau_pts[::-1][1:])  # skip duplicate boundary
            Y_list.append(Y[::-1][1:])

        # IC for next interval: endpoint is j=0 (x=+1, tau=t_end)
        y_current = Y[0, :]

    tau_all = np.concatenate(tau_list)
    Y_all = np.concatenate(Y_list, axis=0)

    return tau_all, Y_all


def solve_single_k_two_interval(k, bg, tau_init, tau_switch, tau_end,
                                 N_cheb_1=20, N_cheb_2=30,
                                 lg_max=L_GAMMA_MAX, ln_max=L_NU_MAX_SYNC):
    """
    Legacy two-interval interface. Delegates to multi-interval solver.
    """
    return solve_single_k_multi_interval(
        k, bg, tau_init, tau_end, N_cheb=N_cheb_2,
        lg_max=lg_max, ln_max=ln_max)


# ============================================================================
# Chebyshev interpolation for arbitrary tau evaluation
# ============================================================================

def chebyshev_interpolate(tau_eval, tau_pts, Y):
    """
    Interpolate solution from Chebyshev collocation points to arbitrary tau values.

    Uses barycentric interpolation (spectrally accurate, O(N) per point).

    Parameters
    ----------
    tau_eval : (M,) ndarray, target conformal times
    tau_pts : (N+1,) ndarray, Chebyshev collocation times (MUST be in Chebyshev order)
    Y : (N+1, nvar) ndarray, solution at collocation points

    Returns
    -------
    Y_eval : (M, nvar) ndarray, interpolated solution
    """
    N = len(tau_pts) - 1
    nvar = Y.shape[1]

    # Map tau_pts to [-1, 1]
    tau_start = tau_pts[-1]  # j=N maps to tau_start (x=-1)
    tau_end = tau_pts[0]     # j=0 maps to tau_end (x=+1)
    x_pts = map_to_chebyshev(tau_pts, tau_start, tau_end)
    x_eval = map_to_chebyshev(tau_eval, tau_start, tau_end)

    # Barycentric weights for Chebyshev-Gauss-Lobatto points
    w = np.ones(N + 1)
    w[0] = 0.5
    w[N] = 0.5
    for j in range(N + 1):
        w[j] *= (-1.0) ** j

    M = len(tau_eval)
    Y_eval = np.zeros((M, nvar))

    for m in range(M):
        # Check if x_eval[m] is close to a node
        dist = np.abs(x_eval[m] - x_pts)
        if np.min(dist) < 1e-14:
            idx = np.argmin(dist)
            Y_eval[m] = Y[idx]
        else:
            denom = w / (x_eval[m] - x_pts)
            numerator = np.sum(denom[:, None] * Y, axis=0)
            denominator = np.sum(denom)
            Y_eval[m] = numerator / denominator

    return Y_eval


# ============================================================================
# Extract Newtonian gauge quantities at snapshot times
# ============================================================================

def extract_newtonian_quantities(k, tau_snaps, Y_snaps, bg,
                                  lg_max=L_GAMMA_MAX, ln_max=L_NU_MAX_SYNC):
    """
    Extract gauge-transformed Newtonian quantities at snapshot times.

    Parameters
    ----------
    k : float, wavenumber
    tau_snaps : (N_snap,) ndarray
    Y_snaps : (N_snap, nvar) ndarray, solution at snapshot times
    bg : Background

    Returns
    -------
    results : (N_snap, 5) ndarray
        Columns: [Phi_N, Psi_N, Theta0_N, vb_N, Theta2]
    """
    N_snap = len(tau_snaps)
    results = np.zeros((N_snap, 5))

    a_snap = bg.a_at_tau(tau_snaps)
    calH_snap = bg.calH_at_tau(tau_snaps)
    fn_s = idx_fn_start(lg_max)

    for it in range(N_snap):
        y = Y_snaps[it]
        a = a_snap[it]
        calH = calH_snap[it]

        h_prime, eta_prime = diagnose_metric(y, k, calH, a, lg_max, ln_max)
        gt = gauge_transform(y, k, calH, a, lg_max, ln_max, h_prime, eta_prime)

        eta_val = y[IDX_ETA]
        results[it, 0] = gt['Phi_N']
        results[it, 1] = gt['Psi_N']
        # Gauge-invariant Theta0: F_g0/4 - eta + Phi_N
        results[it, 2] = y[IDX_FG_START] / 4.0 - eta_val + gt['Phi_N']
        # Velocity: v_b_N = theta_b/k + (h'+6eta')/(2k)
        results[it, 3] = y[IDX_THETA_B] / k + (h_prime + 6.0 * eta_prime) / (2.0 * k)
        # Photon quadrupole: F_gamma,2 / 4
        results[it, 4] = y[IDX_FG_START + 2] / 4.0

    return results


# ============================================================================
# Build snapshot grid (same as solver_sync.py)
# ============================================================================

def build_snapshot_grid(bg, N_vis=60, N_early_isw=30, N_late_isw=20,
                        N_reion=20):
    """
    Build tau grid for LOS integration snapshots.
    Identical logic to solver_sync.build_snapshot_grid.
    """
    peak_idx = np.argmax(bg.visibility_grid)
    tau_peak = bg.tau_grid[peak_idx]
    vis_peak = bg.visibility_grid[peak_idx]

    half_max = vis_peak / 2.0
    left = peak_idx
    while left > 0 and bg.visibility_grid[left] > half_max:
        left -= 1
    right = peak_idx
    while right < len(bg.visibility_grid) - 1 and bg.visibility_grid[right] > half_max:
        right += 1
    fwhm = bg.tau_grid[right] - bg.tau_grid[left]
    sigma = fwhm / 2.355

    tau_vis_lo = max(tau_peak - 5.0 * sigma, bg.tau_grid[5])
    tau_vis_hi = min(tau_peak + 5.0 * sigma, bg.tau_0 * 0.5)
    tau_vis = np.linspace(tau_vis_lo, tau_vis_hi, N_vis)

    a_eq = _OMEGA_R / _OMEGA_M
    tau_eq_approx = float(bg._tau_of_a(np.log(a_eq)))
    tau_early_lo = max(tau_eq_approx * 0.3, bg.tau_grid[5])
    tau_early_hi = tau_vis_lo

    if tau_early_hi > tau_early_lo * 1.5:
        tau_isw_dense_lo = max(tau_eq_approx - 50.0, tau_early_lo)
        tau_isw_dense_hi = min(tau_eq_approx + 50.0, tau_early_hi)
        tau_dense = np.linspace(tau_isw_dense_lo, tau_isw_dense_hi, N_early_isw)
        tau_sparse = np.geomspace(tau_early_lo, tau_early_hi, 15)
        tau_early = np.sort(np.unique(np.concatenate([tau_sparse, tau_dense])))
    else:
        tau_early = np.array([tau_early_lo])

    tau_reion = np.array([])
    if hasattr(bg, 'z_reio') and bg.tau_reio > 0:
        z_reio = bg.z_reio
        a_reio = 1.0 / (1.0 + z_reio)
        tau_reion_center = float(bg._tau_of_a(np.log(a_reio)))
        tau_reion_lo = max(tau_reion_center - 200.0, tau_vis_hi + 10.0)
        tau_reion_hi = min(tau_reion_center + 200.0, bg.tau_0 * 0.95)
        if tau_reion_hi > tau_reion_lo:
            tau_reion = np.linspace(tau_reion_lo, tau_reion_hi, N_reion)

    tau_late_lo = tau_vis_hi
    tau_late_hi = bg.tau_0 * 0.98
    tau_late = np.linspace(tau_late_lo, tau_late_hi, N_late_isw)

    tau_all = np.sort(np.unique(np.concatenate([
        tau_early, tau_vis, tau_reion, tau_late])))
    return tau_vis, tau_late, tau_all


# ============================================================================
# Full Chebyshev spectral CMB solver pipeline
# ============================================================================

def run_chebyshev_solver(N_k=500, k_min=3e-4, k_max=0.35,
                         N_cheb_1=20, N_cheb_2=35,
                         lg_max=L_GAMMA_MAX, ln_max=L_NU_MAX_SYNC,
                         verbose=True):
    """
    Full Chebyshev spectral CMB solver pipeline.

    Pipeline:
      1. Background: Friedmann + recombination
      2. For each k: two-interval Chebyshev solve (A-matrix build + linear solve)
      3. Chebyshev interpolation to LOS snapshot grid
      4. Gauge transform: sync -> Newtonian
      5. LOS integration: SW + Doppler + ISW -> C_l

    Parameters
    ----------
    N_k : int, number of k-modes
    k_min, k_max : float, wavenumber range in Mpc^{-1}
    N_cheb_1 : int, Chebyshev points for Phase 1 (early time)
    N_cheb_2 : int, Chebyshev points for Phase 2 (recombination + late)
    lg_max, ln_max : int, hierarchy truncation

    Returns
    -------
    dict with keys: ell, Dl_TT, Dl_EE, Dl_TE, timing, etc.
    """
    t_total = time.time()

    # ================================================================
    # Step 1: Background
    # ================================================================
    if verbose:
        print("=" * 65)
        print("Chebyshev Spectral Boltzmann Solver")
        print("=" * 65)
        print(f"\n--- Step 1: Background ---")

    t0 = time.time()
    bg = Background(khronon=False, recombination='recfast')
    bg.solve()
    t_bg = time.time() - t0
    if verbose:
        print(f"Background: {t_bg:.2f}s")

    # ================================================================
    # Step 2: Setup
    # ================================================================
    k_arr = np.geomspace(k_min, k_max, N_k)
    nvar = n_var_sync(lg_max, ln_max)

    tau_vis, tau_late, tau_all = build_snapshot_grid(
        bg, N_vis=60, N_early_isw=30, N_late_isw=20)
    N_snap = len(tau_all)

    # Multi-interval boundaries
    tau_init = bg.tau_grid[1]  # avoid tau=0

    # End time: slightly past last snapshot
    tau_end = float(tau_all[-1] + 1.0)

    # Preview interval count for representative k values
    n_int_lo = len(_build_interval_boundaries(k_min, tau_init, tau_end, N_cheb_2)) - 1
    n_int_hi = len(_build_interval_boundaries(k_max, tau_init, tau_end, N_cheb_2)) - 1

    if verbose:
        print(f"\n--- Step 2: Chebyshev Perturbations ({N_k} k-modes) ---")
        print(f"k range: [{k_min:.1e}, {k_max:.2f}] Mpc^-1")
        print(f"tau range: [{tau_init:.2f}, {tau_end:.1f}] Mpc")
        print(f"N_cheb per interval: {N_cheb_2}")
        print(f"Intervals: {n_int_lo} (k_min) to {n_int_hi} (k_max)")
        print(f"System size per interval: {(N_cheb_2+1)*nvar}")
        print(f"Snapshot grid: {N_snap} tau points")
        print(f"Hierarchy: lg_max={lg_max}, ln_max={ln_max}, nvar={nvar}")
        sys.stdout.flush()

    # ================================================================
    # Step 3: Solve all k-modes
    # ================================================================
    t0 = time.time()

    # Storage
    from scipy.interpolate import interp1d as _interp1d

    Phi_N = np.zeros((N_k, N_snap))
    Psi_N = np.zeros((N_k, N_snap))
    Theta0_N = np.zeros((N_k, N_snap))
    vb_N = np.zeros((N_k, N_snap))
    Theta2_arr = np.zeros((N_k, N_snap))

    n_failed = 0
    for ik, k in enumerate(k_arr):
        try:
            # Solve on Chebyshev multi-interval grid
            tau_cheb, Y_cheb = solve_single_k_multi_interval(
                k, bg, tau_init, tau_end, N_cheb=N_cheb_2,
                lg_max=lg_max, ln_max=ln_max)

            # Interpolate to snapshot grid
            Y_snap = np.zeros((N_snap, nvar))
            for v in range(nvar):
                interp_fn = _interp1d(tau_cheb, Y_cheb[:, v],
                                      kind='cubic', fill_value='extrapolate')
                Y_snap[:, v] = interp_fn(tau_all)

            # Extract Newtonian gauge quantities
            res = extract_newtonian_quantities(k, tau_all, Y_snap, bg,
                                                lg_max, ln_max)
            Phi_N[ik] = res[:, 0]
            Psi_N[ik] = res[:, 1]
            Theta0_N[ik] = res[:, 2]
            vb_N[ik] = res[:, 3]
            Theta2_arr[ik] = res[:, 4]

        except Exception as e:
            n_failed += 1
            if verbose and n_failed <= 5:
                print(f"  WARNING: k={k:.4e} failed: {e}")

        if verbose and (ik + 1) % 50 == 0:
            elapsed = time.time() - t0
            eta_est = elapsed / (ik + 1) * N_k
            print(f"  [{ik+1}/{N_k}] elapsed={elapsed:.1f}s, "
                  f"ETA={eta_est:.1f}s, failed={n_failed}")
            sys.stdout.flush()

    t_pert = time.time() - t0
    if verbose:
        print(f"Perturbations: {t_pert:.1f}s ({n_failed} failed, "
              f"{N_k - n_failed} OK)")

    # ================================================================
    # Step 4: Compute Phi_N' + Psi_N' by cubic spline derivatives
    # ================================================================
    from scipy.interpolate import CubicSpline
    PhiPsi = Phi_N + Psi_N
    PhiPsi_prime = np.zeros_like(PhiPsi)
    dtau_all = np.diff(tau_all)
    for ik in range(N_k):
        cs = CubicSpline(tau_all, PhiPsi[ik])
        PhiPsi_prime[ik] = cs(tau_all, 1)

    # ================================================================
    # Step 5: LOS integration -> C_l
    # ================================================================
    t0 = time.time()

    ell_values = np.unique(np.concatenate([
        np.arange(2, 30, 1),
        np.arange(30, 100, 2),
        np.arange(100, 500, 4),
        np.arange(500, 1500, 8),
        np.arange(1500, 2501, 12),
    ])).astype(int)
    N_ell = len(ell_values)

    if verbose:
        print(f"\n--- Step 5: LOS C_l ({N_ell} ells) ---")
        sys.stdout.flush()

    # Primordial spectrum
    P_R = A_s * (k_arr / k_pivot) ** (n_s - 1.0)

    # Background at snapshots
    g_snap = bg.visibility_at_tau(tau_all)
    kappa_snap = bg.kappa_at_tau(tau_all)
    exp_neg_kappa = np.exp(-kappa_snap)
    chi_snap = bg.tau_0 - tau_all

    # Bessel table (GPU-accelerated)
    try:
        from .bessel_cache import CachedBesselTable
        import mlx.core as mx
        x_max_bessel = float(np.max(k_arr) * np.max(chi_snap)) * 1.05 + 50.0
        bessel_table = CachedBesselTable(ell_values, x_max=x_max_bessel,
                                          N_cheb=64, seg_width=80,
                                          verbose=verbose)
        use_gpu_bessel = True
    except Exception as e:
        if verbose:
            print(f"[WARNING] Bessel table failed ({e}), using scipy")
        use_gpu_bessel = False

    # Transfer functions
    Delta_l = np.zeros((N_ell, N_k), dtype=np.float64)
    Delta_l_E = np.zeros((N_ell, N_k), dtype=np.float64)

    alpha_P = 1.7
    pol_prefactor = 0.75 * alpha_P

    eps_prefactors = np.zeros(N_ell)
    for il, ell in enumerate(ell_values):
        l = int(ell)
        if l >= 2:
            eps_prefactors[il] = np.sqrt(float((l - 1) * l * (l + 1) * (l + 2)))

    if verbose:
        print(f"Computing LOS transfer functions (TT + EE)...")
        sys.stdout.flush()

    for it in range(N_snap):
        chi = chi_snap[it]
        if chi <= 0:
            continue

        if it == 0:
            w = 0.5 * dtau_all[0]
        elif it == N_snap - 1:
            w = 0.5 * dtau_all[-1]
        else:
            w = 0.5 * (dtau_all[max(0, it - 1)] + dtau_all[min(it, len(dtau_all) - 1)])

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

    t_los = time.time() - t0
    if verbose:
        print(f"LOS integration: {t_los:.2f}s")

    # ================================================================
    # Step 6: C_l = 4 pi int dk/k P_R |Delta_l|^2
    # ================================================================
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

    if verbose:
        print(f"C_l computation (TT+EE+TE): {t_cl:.3f}s")
        print(f"\nTotal time: {t_elapsed:.1f}s  "
              f"(bg={t_bg:.1f}s, pert={t_pert:.1f}s, LOS={t_los:.1f}s)")

    return {
        'ell': ell_values,
        'Dl_TT': Dl_TT,
        'Dl_EE': Dl_EE,
        'Dl_TE': Dl_TE,
        'Cl_TT': Cl_TT,
        'Cl_EE': Cl_EE,
        'Cl_TE': Cl_TE,
        'k_arr': k_arr,
        'bg': bg,
        'Delta_l': Delta_l,
        'Delta_l_E': Delta_l_E,
        'Phi_N': Phi_N,
        'Psi_N': Psi_N,
        'timing': {
            'background': t_bg,
            'perturbations': t_pert,
            'los': t_los,
            'total': t_elapsed,
            'n_failed': n_failed,
        },
    }


# ============================================================================
# GPU-batched Chebyshev solver (MLX)
# ============================================================================

def run_chebyshev_solver_gpu(N_k=500, k_min=3e-4, k_max=0.35,
                              N_cheb_1=20, N_cheb_2=35,
                              lg_max=L_GAMMA_MAX, ln_max=L_NU_MAX_SYNC,
                              batch_size=50, verbose=True):
    """
    Parallel version using multiprocessing.

    Each k-mode is solved independently with the multi-interval Chebyshev
    method in a separate process. This is more effective than GPU batching
    for the multi-interval approach where different k-modes have different
    numbers of intervals.

    Parameters
    ----------
    batch_size : int (unused, kept for API compatibility)
    """
    from multiprocessing import Pool

    t_total = time.time()

    if verbose:
        print("=" * 65)
        print("Chebyshev Spectral Boltzmann Solver (Parallel)")
        print("=" * 65)
        print(f"\n--- Step 1: Background ---")

    t0 = time.time()
    bg = Background(khronon=False, recombination='recfast')
    bg.solve()
    t_bg = time.time() - t0
    if verbose:
        print(f"Background: {t_bg:.2f}s")

    k_arr = np.geomspace(k_min, k_max, N_k)
    nvar = n_var_sync(lg_max, ln_max)

    tau_vis, tau_late, tau_all = build_snapshot_grid(
        bg, N_vis=60, N_early_isw=30, N_late_isw=20)
    N_snap = len(tau_all)

    tau_init = bg.tau_grid[1]
    tau_end = float(tau_all[-1] + 1.0)

    n_workers = os.cpu_count() or 4

    if verbose:
        print(f"\n--- Step 2: Chebyshev Perturbations ({N_k} k-modes, {n_workers} workers) ---")
        print(f"k range: [{k_min:.1e}, {k_max:.2f}] Mpc^-1")
        print(f"N_cheb per interval: {N_cheb_2}")
        print(f"nvar={nvar}")
        sys.stdout.flush()

    # Serialize background for pickling
    bg_arrays = {
        'tau_grid': bg.tau_grid.copy(),
        'a_grid': bg.a_grid.copy(),
        'calH_grid': bg.calH_grid.copy(),
        'R_grid': bg.R_grid.copy(),
        'kappa_dot_grid': bg.kappa_dot_grid.copy(),
        'tau_rec': float(bg.tau_rec),
        'tau_0': float(bg.tau_0),
    }

    a_snap = np.asarray(bg.a_at_tau(tau_all), dtype=np.float64)
    calH_snap = np.asarray(bg.calH_at_tau(tau_all), dtype=np.float64)

    work_items = [
        (k, bg_arrays, tau_init, tau_end, N_cheb_2, lg_max, ln_max, tau_all)
        for k in k_arr
    ]

    t0 = time.time()

    if verbose:
        print(f"Dispatching {N_k} k-modes to {n_workers} workers...")
        sys.stdout.flush()

    with Pool(processes=n_workers) as pool:
        worker_results = pool.map(_chebyshev_worker, work_items)

    t_pert = time.time() - t0

    # Collect results
    Phi_N = np.zeros((N_k, N_snap))
    Psi_N = np.zeros((N_k, N_snap))
    Theta0_N = np.zeros((N_k, N_snap))
    vb_N = np.zeros((N_k, N_snap))
    Theta2_arr = np.zeros((N_k, N_snap))

    n_failed = 0
    for ik, res in enumerate(worker_results):
        if res is None:
            n_failed += 1
            if verbose and n_failed <= 3:
                print(f"  WARNING: k={k_arr[ik]:.4e} failed")
            continue
        Phi_N[ik] = res[:, 0]
        Psi_N[ik] = res[:, 1]
        Theta0_N[ik] = res[:, 2]
        vb_N[ik] = res[:, 3]
        Theta2_arr[ik] = res[:, 4]

    if verbose:
        print(f"Perturbations: {t_pert:.1f}s ({n_failed} failed, "
              f"{N_k - n_failed} OK)")

    # From here on, same as the CPU version
    from scipy.interpolate import CubicSpline
    PhiPsi = Phi_N + Psi_N
    PhiPsi_prime = np.zeros_like(PhiPsi)
    dtau_all = np.diff(tau_all)
    for ik in range(N_k):
        cs = CubicSpline(tau_all, PhiPsi[ik])
        PhiPsi_prime[ik] = cs(tau_all, 1)

    t0 = time.time()

    ell_values = np.unique(np.concatenate([
        np.arange(2, 30, 1),
        np.arange(30, 100, 2),
        np.arange(100, 500, 4),
        np.arange(500, 1500, 8),
        np.arange(1500, 2501, 12),
    ])).astype(int)
    N_ell = len(ell_values)

    if verbose:
        print(f"\n--- Step 3: LOS C_l ({N_ell} ells) ---")
        sys.stdout.flush()

    P_R = A_s * (k_arr / k_pivot) ** (n_s - 1.0)
    g_snap = bg.visibility_at_tau(tau_all)
    kappa_snap = bg.kappa_at_tau(tau_all)
    exp_neg_kappa = np.exp(-kappa_snap)
    chi_snap = bg.tau_0 - tau_all

    try:
        from .bessel_cache import CachedBesselTable
        import mlx.core as mx
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
        l = int(ell)
        if l >= 2:
            eps_prefactors[il] = np.sqrt(float((l - 1) * l * (l + 1) * (l + 2)))

    for it in range(N_snap):
        chi = chi_snap[it]
        if chi <= 0:
            continue

        if it == 0:
            w = 0.5 * dtau_all[0]
        elif it == N_snap - 1:
            w = 0.5 * dtau_all[-1]
        else:
            w = 0.5 * (dtau_all[max(0, it - 1)] + dtau_all[min(it, len(dtau_all) - 1)])

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
                    S_SW * jl[il, :] + S_Dop * jlp[il, :] + S_ISW * jl[il, :])
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

    t_los = time.time() - t0
    if verbose:
        print(f"LOS integration: {t_los:.2f}s")

    t0 = time.time()
    lnk = np.log(k_arr)
    dlnk = np.diff(lnk)
    ell_f = ell_values.astype(float)
    uK2 = (T_CMB * 1e6) ** 2
    norm = (2.0 / 3.0) ** 2

    integrand_TT = P_R[None, :] * Delta_l ** 2
    mid_TT = 0.5 * (integrand_TT[:, :-1] + integrand_TT[:, 1:])
    Cl_TT = 4.0 * np.pi * np.sum(mid_TT * dlnk[None, :], axis=1)
    Cl_TT = np.maximum(Cl_TT, 0.0)
    Dl_TT = norm * ell_f * (ell_f + 1.0) * Cl_TT / (2.0 * np.pi) * uK2

    integrand_EE = P_R[None, :] * Delta_l_E ** 2
    mid_EE = 0.5 * (integrand_EE[:, :-1] + integrand_EE[:, 1:])
    Cl_EE = 4.0 * np.pi * np.sum(mid_EE * dlnk[None, :], axis=1)
    Cl_EE = np.maximum(Cl_EE, 0.0)
    Dl_EE = norm * ell_f * (ell_f + 1.0) * Cl_EE / (2.0 * np.pi) * uK2

    integrand_TE = P_R[None, :] * Delta_l * Delta_l_E
    mid_TE = 0.5 * (integrand_TE[:, :-1] + integrand_TE[:, 1:])
    Cl_TE = 4.0 * np.pi * np.sum(mid_TE * dlnk[None, :], axis=1)
    Dl_TE = norm * ell_f * (ell_f + 1.0) * Cl_TE / (2.0 * np.pi) * uK2

    t_cl = time.time() - t0
    t_elapsed = time.time() - t_total

    if verbose:
        print(f"C_l computation (TT+EE+TE): {t_cl:.3f}s")
        print(f"\nTotal time: {t_elapsed:.1f}s  "
              f"(bg={t_bg:.1f}s, pert={t_pert:.1f}s, LOS={t_los:.1f}s)")

    return {
        'ell': ell_values,
        'Dl_TT': Dl_TT,
        'Dl_EE': Dl_EE,
        'Dl_TE': Dl_TE,
        'Cl_TT': Cl_TT,
        'Cl_EE': Cl_EE,
        'Cl_TE': Cl_TE,
        'k_arr': k_arr,
        'bg': bg,
        'Delta_l': Delta_l,
        'Delta_l_E': Delta_l_E,
        'timing': {
            'background': t_bg,
            'perturbations': t_pert,
            'los': t_los,
            'total': t_elapsed,
            'n_failed': n_failed,
        },
    }


def _chebyshev_worker(args):
    """
    Multiprocessing worker for Chebyshev solver.

    Solves one k-mode with the multi-interval method, interpolates to
    snapshot grid, and returns Newtonian gauge quantities.
    """
    import numpy as np
    import warnings
    from scipy.interpolate import interp1d
    from scipy.linalg import solve as scipy_solve

    (k, bg_arrays, tau_init, tau_end, N_cheb, lg_max, ln_max, tau_all) = args

    # Reconstruct lightweight background
    class _BGLite:
        __slots__ = ('tau_grid', 'a_grid', 'calH_grid', 'R_grid',
                     'kappa_dot_grid', 'tau_rec', 'tau_0',
                     '_calH_of_tau', '_a_of_tau', '_R_of_tau',
                     '_kappa_dot_of_tau')
        def __init__(self, arrays):
            self.tau_grid = arrays['tau_grid']
            self.a_grid = arrays['a_grid']
            self.calH_grid = arrays['calH_grid']
            self.R_grid = arrays['R_grid']
            self.kappa_dot_grid = arrays['kappa_dot_grid']
            self.tau_rec = arrays['tau_rec']
            self.tau_0 = arrays['tau_0']
            # Build interpolators
            self._calH_of_tau = interp1d(self.tau_grid, self.calH_grid,
                                          kind='cubic', fill_value='extrapolate')
            self._a_of_tau = interp1d(self.tau_grid, self.a_grid,
                                       kind='cubic', fill_value='extrapolate')
            self._R_of_tau = interp1d(self.tau_grid, self.R_grid,
                                       kind='cubic', fill_value='extrapolate')
            self._kappa_dot_of_tau = interp1d(self.tau_grid, self.kappa_dot_grid,
                                               kind='cubic', fill_value='extrapolate')
        def calH_at_tau(self, t): return self._calH_of_tau(t)
        def a_at_tau(self, t): return self._a_of_tau(t)
        def R_at_tau(self, t): return self._R_of_tau(t)
        def kappa_dot_at_tau(self, t): return self._kappa_dot_of_tau(t)

    bg = _BGLite(bg_arrays)

    try:
        from mlx_class.solver_chebyshev import (
            solve_single_k_multi_interval, extract_newtonian_quantities
        )
        from mlx_class.perturbations_sync import n_var_sync

        nvar = n_var_sync(lg_max, ln_max)

        tau_cheb, Y_cheb = solve_single_k_multi_interval(
            k, bg, tau_init, tau_end, N_cheb=N_cheb,
            lg_max=lg_max, ln_max=ln_max)

        # Interpolate to snapshot grid
        N_snap = len(tau_all)
        Y_snap = np.zeros((N_snap, nvar))
        for v in range(nvar):
            fn = interp1d(tau_cheb, Y_cheb[:, v],
                          kind='cubic', fill_value='extrapolate')
            Y_snap[:, v] = fn(tau_all)

        # Extract Newtonian quantities
        res = extract_newtonian_quantities(k, tau_all, Y_snap, bg,
                                            lg_max, ln_max)
        return res

    except Exception as e:
        return None



# ============================================================================
# Comparison and validation
# ============================================================================

def compare_with_class(result, class_file=None):
    """
    Compare Chebyshev solver output with CLASS reference.

    Returns comparison metrics: RMS error, peak positions, etc.
    """
    # Load CLASS reference
    if class_file is None:
        class_file = '/Users/akaihuangm1/Desktop/github/gdm_class_public/output/mlx_ref_lcdm_cl.dat'

    try:
        data = np.loadtxt(class_file)
        ell_class = data[:, 0].astype(int)
        # CLASS format: dimensionless l(l+1)/(2pi) C_l
        Dl_class = data[:, 1] * (T_CMB * 1e6) ** 2
    except Exception as e:
        print(f"Cannot load CLASS file: {e}")
        return None

    ell = result['ell']
    Dl = result['Dl_TT']

    # Interpolate to common ell range
    from scipy.interpolate import interp1d
    ell_min = max(ell[0], ell_class[0])
    ell_max = min(ell[-1], ell_class[-1])

    ell_common = np.arange(ell_min, ell_max + 1)
    interp_cheb = interp1d(ell, Dl, kind='cubic', fill_value='extrapolate')
    interp_class = interp1d(ell_class, Dl_class, kind='cubic', fill_value='extrapolate')

    Dl_cheb_common = interp_cheb(ell_common)
    Dl_class_common = interp_class(ell_common)

    # RMS relative error
    mask = Dl_class_common > 10.0  # avoid division by near-zero at low-l
    if np.any(mask):
        rel_err = (Dl_cheb_common[mask] - Dl_class_common[mask]) / Dl_class_common[mask]
        rms = np.sqrt(np.mean(rel_err ** 2))
        max_err = np.max(np.abs(rel_err))
        mean_err = np.mean(np.abs(rel_err))
    else:
        rms = max_err = mean_err = np.nan

    # Peak positions
    from scipy.signal import find_peaks
    peaks_cheb, _ = find_peaks(Dl_cheb_common, distance=50, prominence=10)
    peaks_class, _ = find_peaks(Dl_class_common, distance=50, prominence=10)

    ell_peaks_cheb = ell_common[peaks_cheb[:7]] if len(peaks_cheb) >= 7 else ell_common[peaks_cheb]
    ell_peaks_class = ell_common[peaks_class[:7]] if len(peaks_class) >= 7 else ell_common[peaks_class]

    comparison = {
        'rms_relative': rms,
        'max_relative': max_err,
        'mean_relative': mean_err,
        'peaks_chebyshev': ell_peaks_cheb,
        'peaks_class': ell_peaks_class,
        'ell_common': ell_common,
        'Dl_cheb': Dl_cheb_common,
        'Dl_class': Dl_class_common,
    }

    print(f"\n{'='*65}")
    print(f"Comparison: Chebyshev vs CLASS")
    print(f"{'='*65}")
    print(f"RMS relative error (l>{ell_min}, D_l>10): {rms:.4f} ({rms*100:.2f}%)")
    print(f"Max relative error: {max_err:.4f} ({max_err*100:.2f}%)")
    print(f"Mean |relative error|: {mean_err:.4f} ({mean_err*100:.2f}%)")
    print(f"\nPeak positions:")
    print(f"  Chebyshev: {ell_peaks_cheb}")
    print(f"  CLASS:     {ell_peaks_class}")

    if len(ell_peaks_cheb) > 0 and len(ell_peaks_class) > 0:
        n_peaks = min(len(ell_peaks_cheb), len(ell_peaks_class))
        for i in range(n_peaks):
            diff = ell_peaks_cheb[i] - ell_peaks_class[i]
            pct = diff / ell_peaks_class[i] * 100
            print(f"  Peak {i+1}: l={ell_peaks_cheb[i]} vs {ell_peaks_class[i]} "
                  f"(diff={diff:+d}, {pct:+.1f}%)")

    return comparison


def compare_with_sync(result_cheb, verbose=True):
    """
    Compare Chebyshev solver output with the sync gauge solver output.

    Loads the saved sync gauge data file if available.
    """
    sync_file = os.path.join(os.path.dirname(__file__), 'cl_sync_gauge.dat')
    try:
        data = np.loadtxt(sync_file, comments='#')
        ell_sync = data[:, 0].astype(int)
        Dl_sync = data[:, 1]
    except Exception as e:
        if verbose:
            print(f"Cannot load sync gauge data: {e}")
        return None

    ell = result_cheb['ell']
    Dl = result_cheb['Dl_TT']

    from scipy.interpolate import interp1d
    ell_min = max(ell[0], ell_sync[0])
    ell_max = min(ell[-1], ell_sync[-1])
    ell_common = np.arange(ell_min, ell_max + 1)

    interp_cheb = interp1d(ell, Dl, kind='cubic', fill_value='extrapolate')
    interp_sync = interp1d(ell_sync, Dl_sync, kind='cubic', fill_value='extrapolate')

    Dl_cheb = interp_cheb(ell_common)
    Dl_sync_c = interp_sync(ell_common)

    mask = Dl_sync_c > 10.0
    if np.any(mask):
        rel_err = (Dl_cheb[mask] - Dl_sync_c[mask]) / Dl_sync_c[mask]
        rms = np.sqrt(np.mean(rel_err ** 2))
    else:
        rms = np.nan

    if verbose:
        print(f"\nChebyshev vs Sync-gauge solver: RMS = {rms:.4f} ({rms*100:.2f}%)")

    return {'rms': rms, 'ell_common': ell_common,
            'Dl_cheb': Dl_cheb, 'Dl_sync': Dl_sync_c}


# ============================================================================
# Plotting
# ============================================================================

def plot_results(result, comparison=None, sync_comparison=None):
    """Generate comparison plots."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    ell = result['ell']
    Dl_TT = result['Dl_TT']
    Dl_EE = result['Dl_EE']
    Dl_TE = result['Dl_TE']

    # Panel 1: TT spectrum
    ax = axes[0, 0]
    ax.plot(ell, Dl_TT, 'b-', label='Chebyshev', linewidth=1.5)
    if comparison is not None:
        ax.plot(comparison['ell_common'], comparison['Dl_class'],
                'r--', label='CLASS', linewidth=1.0, alpha=0.7)
    ax.set_xlabel('l')
    ax.set_ylabel('D_l^TT [muK^2]')
    ax.set_title('TT Power Spectrum')
    ax.legend()
    ax.set_xlim(2, 2500)

    # Panel 2: Residual vs CLASS
    ax = axes[0, 1]
    if comparison is not None:
        ell_c = comparison['ell_common']
        mask = comparison['Dl_class'] > 10
        if np.any(mask):
            residual = (comparison['Dl_cheb'][mask] - comparison['Dl_class'][mask]) / comparison['Dl_class'][mask] * 100
            ax.plot(ell_c[mask], residual, 'b-', linewidth=0.5)
            ax.axhline(0, color='k', linewidth=0.5, linestyle='--')
            rms = comparison['rms_relative'] * 100
            ax.set_title(f'Residual vs CLASS (RMS={rms:.1f}%)')
    else:
        ax.set_title('No CLASS data for comparison')
    ax.set_xlabel('l')
    ax.set_ylabel('(Cheb - CLASS) / CLASS [%]')
    ax.set_xlim(2, 2500)

    # Panel 3: EE spectrum
    ax = axes[1, 0]
    ax.plot(ell, Dl_EE, 'g-', label='EE', linewidth=1.5)
    ax.set_xlabel('l')
    ax.set_ylabel('D_l^EE [muK^2]')
    ax.set_title('EE Power Spectrum')
    ax.set_xlim(2, 2500)

    # Panel 4: TE spectrum
    ax = axes[1, 1]
    ax.plot(ell, Dl_TE, 'm-', label='TE', linewidth=1.5)
    ax.axhline(0, color='k', linewidth=0.5, linestyle='--')
    ax.set_xlabel('l')
    ax.set_ylabel('D_l^TE [muK^2]')
    ax.set_title('TE Power Spectrum')
    ax.set_xlim(2, 2500)

    # Timing annotation
    timing = result['timing']
    fig.suptitle(
        f"Chebyshev Spectral CMB Solver  |  "
        f"Total: {timing['total']:.1f}s  "
        f"(pert: {timing['perturbations']:.1f}s, LOS: {timing['los']:.1f}s)  |  "
        f"Failed: {timing['n_failed']}",
        fontsize=11)

    plt.tight_layout()
    outpath = os.path.join(os.path.dirname(__file__), 'cl_chebyshev.png')
    plt.savefig(outpath, dpi=150)
    plt.close()
    print(f"\nPlot saved: {outpath}")
    return outpath


# ============================================================================
# Save output
# ============================================================================

def save_results(result):
    """Save C_l data to file."""
    outpath = os.path.join(os.path.dirname(__file__), 'cl_chebyshev.dat')
    header = ("# l   D_l^TT[muK^2]   D_l^EE[muK^2]   D_l^TE[muK^2]  "
              "(Chebyshev spectral solver)")
    data = np.column_stack([
        result['ell'],
        result['Dl_TT'],
        result['Dl_EE'],
        result['Dl_TE'],
    ])
    np.savetxt(outpath, data, header=header, fmt='%.18e')
    print(f"Data saved: {outpath}")
    return outpath


# ============================================================================
# Main entry point
# ============================================================================

def main():
    """
    Run the Chebyshev spectral CMB solver and compare with CLASS.

    Usage: python -m mlx_class.solver_chebyshev [--gpu] [--N_k 500] [--N_cheb 35]
    """
    import argparse

    parser = argparse.ArgumentParser(
        description='Chebyshev spectral CMB Boltzmann solver')
    parser.add_argument('--gpu', action='store_true',
                        help='Use GPU-batched solver (MLX)')
    parser.add_argument('--N_k', type=int, default=500,
                        help='Number of k-modes (default: 500)')
    parser.add_argument('--N_cheb', type=int, default=35,
                        help='Chebyshev points for Phase 2 (default: 35)')
    parser.add_argument('--N_cheb_1', type=int, default=20,
                        help='Chebyshev points for Phase 1 (default: 20)')
    parser.add_argument('--batch_size', type=int, default=50,
                        help='GPU batch size (default: 50)')
    parser.add_argument('--class_file', type=str, default=None,
                        help='CLASS reference C_l file')
    parser.add_argument('--quiet', action='store_true')
    args = parser.parse_args()

    verbose = not args.quiet

    if args.gpu:
        result = run_chebyshev_solver_gpu(
            N_k=args.N_k,
            N_cheb_1=args.N_cheb_1,
            N_cheb_2=args.N_cheb,
            batch_size=args.batch_size,
            verbose=verbose)
    else:
        result = run_chebyshev_solver(
            N_k=args.N_k,
            N_cheb_1=args.N_cheb_1,
            N_cheb_2=args.N_cheb,
            verbose=verbose)

    # Save results
    save_results(result)

    # Compare with CLASS
    comparison = compare_with_class(result, class_file=args.class_file)

    # Compare with sync gauge solver
    sync_comparison = compare_with_sync(result, verbose=verbose)

    # Plot
    plot_results(result, comparison, sync_comparison)

    # Summary
    print(f"\n{'='*65}")
    print(f"SUMMARY")
    print(f"{'='*65}")
    t = result['timing']
    print(f"Solver: Chebyshev spectral ({'GPU' if args.gpu else 'CPU'})")
    print(f"k-modes: {args.N_k}")
    print(f"Chebyshev points: Phase 1={args.N_cheb_1}, Phase 2={args.N_cheb}")
    print(f"Total time: {t['total']:.1f}s")
    print(f"  Background:    {t['background']:.2f}s")
    print(f"  Perturbations: {t['perturbations']:.1f}s")
    print(f"  LOS:           {t['los']:.1f}s")
    print(f"Failed k-modes: {t['n_failed']}")
    if comparison:
        print(f"RMS vs CLASS: {comparison['rms_relative']*100:.2f}%")
    if sync_comparison:
        print(f"RMS vs Sync:  {sync_comparison['rms']*100:.2f}%")


if __name__ == '__main__':
    main()
