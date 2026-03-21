"""
solver_fast.py — Optimized ODE solver for the Boltzmann hierarchy.

PERFORMANCE (500 k-modes, l_gamma_max=25, l_nu_max=20, Apple M1):
  perturbations_hires.py:  8.6s  (1099 steps, Python loops over l)
  solver_fast.py (fast):   0.5s  ( 599 steps, vectorized + compiled)   17x
  solver_fast.py (hires):  0.9s  (1099 steps, vectorized + compiled)    9x

ACCURACY:
  Same-grid:    corr=1.000000, RMS=0.006%  (numerically identical)
  Reduced-grid: corr=0.99994,  RMS=1.15%   (from coarser time steps)

OPTIMIZATIONS:
  1. Vectorized hierarchy: Replace Python for-loops over multipoles with
     matrix operations (tridiagonal coupling via slicing).
  2. mx.compile: Compile the step functions so MLX traces the graph once.
  3. Reduced eval frequency: Only sync GPU 3-4 times during integration.
  4. Reduced grid: Coarser time steps in TCA phase where solution is smooth.
  5. Precomputed scalars: Background quantities packed into arrays, no per-step
     Python float conversions.

PHYSICS: Identical to perturbations_hires.py:
  - Strang splitting: collision (analytic) + streaming (RK4)
  - IMEX exponential integrator for the stiff Phi equation
  - TCA -> full hierarchy transition
  - Neutrino hierarchy with l_nu_max

Author: Sheng-Kai Huang, 2026
"""
import numpy as np
import mlx.core as mx
import time

from .background import (
    Omega_r as _OMEGA_R,
    Omega_b as _OMEGA_B,
    Omega_c as _OMEGA_C,
    H0_Mpc as _H0_MPC,
)

from .perturbations_neutrino import (
    N_EFF, _f_nu, _OMEGA_GAMMA, _OMEGA_NU, L_NU_MAX,
)


# ============================================================================
# Variable layout (same as perturbations_hires.py)
# ============================================================================

IDX_PHI = 0
IDX_DELTA_B = 1
IDX_V_B = 2
IDX_DELTA_C = 3
IDX_V_C = 4
IDX_THETA_START = 5


def idx_theta(l):
    return IDX_THETA_START + l


def idx_n_start(l_gamma_max):
    return IDX_THETA_START + l_gamma_max + 1


def n_var_total(l_gamma_max, l_nu_max):
    return IDX_THETA_START + (l_gamma_max + 1) + (l_nu_max + 1)


# ============================================================================
# Phi IMEX helpers
# ============================================================================

def _phi_lambda(k_arr, calH):
    return -(k_arr * k_arr / (3.0 * calH) + calH)


def _phi_int_factor(lam, dt):
    lam_dt = lam * dt
    exact = (mx.exp(lam_dt) - 1.0) / lam
    taylor = dt * (1.0 + 0.5 * lam_dt + lam_dt * lam_dt / 6.0)
    return mx.where(mx.abs(lam_dt) < 1e-4, taylor, exact)


def _set_phi(y, Phi_new):
    return mx.concatenate([Phi_new[:, None], y[:, 1:]], axis=1)


# ============================================================================
# Vectorized streaming derivative (NO Python loops over multipoles)
# ============================================================================

def _deriv_streaming_full(y, k_arr, calH, a, R, tau_val,
                          Og, On, Ob, Oc, H0, lgmax, lnmax, _n_start):
    """
    Vectorized streaming + gravity RHS for the full hierarchy phase.

    The photon hierarchy l=3..l_max-1 is computed in one vectorized operation
    instead of a Python loop. Same for neutrinos l=2..l_nu_max-1.
    """
    Phi = y[:, IDX_PHI]
    delta_b = y[:, IDX_DELTA_B]
    v_b = y[:, IDX_V_B]
    delta_c = y[:, IDX_DELTA_C]
    v_c = y[:, IDX_V_C]
    Theta_0 = y[:, idx_theta(0)]
    Theta_1 = y[:, idx_theta(1)]

    N_0 = y[:, _n_start]
    N_1 = y[:, _n_start + 1]
    N_2 = y[:, _n_start + 2]

    H02 = H0 * H0
    inv_a2 = 1.0 / (a * a)

    # Psi with anisotropic stress
    Theta_2 = y[:, idx_theta(2)]
    k2 = k_arr * k_arr
    Psi = Phi - 12.0 * H02 * inv_a2 / k2 * (On * N_2 + Og * Theta_2)

    # Phi_dot (Psi=Phi for stability in Phi evolution)
    S = (Og * inv_a2 * 4.0 * Theta_0
         + On * inv_a2 * 4.0 * N_0
         + Ob / a * delta_b
         + Oc / a * delta_c)
    Phi_dot = -calH * Phi - k2 * Phi / (3.0 * calH) - 0.5 * H02 / calH * S

    # ---- Photon hierarchy (vectorized) ----
    # l=0
    dTheta_0 = -k_arr * Theta_1 - Phi_dot
    # l=1 (streaming only)
    dTheta_1 = k_arr / 3.0 * (Theta_0 + Psi)

    # l=2
    Theta_3 = y[:, idx_theta(3)]
    dTheta_2 = k_arr / 5.0 * (2.0 * Theta_1 - 3.0 * Theta_3)

    # l=3 to l_max-1: VECTORIZED
    # dTheta_l = k/(2l+1) * [l * Theta_{l-1} - (l+1) * Theta_{l+1}]
    # Shape: y[:, idx_theta(3):idx_theta(lgmax)] has lgmax-3 columns
    if lgmax > 3:
        l_vals = mx.arange(3, lgmax).astype(mx.float32)  # shape (lgmax-3,)
        inv_2lp1 = 1.0 / (2.0 * l_vals + 1.0)  # (lgmax-3,)
        # Theta_{l-1} for l=3..lgmax-1: Theta_2..Theta_{lgmax-2}
        theta_lm1 = y[:, idx_theta(2):idx_theta(lgmax - 1)]  # (N_k, lgmax-3)
        # Theta_{l+1} for l=3..lgmax-1: Theta_4..Theta_{lgmax}
        theta_lp1 = y[:, idx_theta(4):idx_theta(lgmax + 1)]  # (N_k, lgmax-3)
        dTheta_mid = k_arr[:, None] * inv_2lp1[None, :] * (
            l_vals[None, :] * theta_lm1 - (l_vals[None, :] + 1.0) * theta_lp1
        )
    else:
        dTheta_mid = None

    # l=l_max: truncation
    tau_safe = mx.maximum(tau_val, mx.array(1e-10))
    Theta_lm1_last = y[:, idx_theta(lgmax - 1)]
    Theta_lmax = y[:, idx_theta(lgmax)]
    dTheta_lmax = (k_arr * Theta_lm1_last * lgmax / (2.0 * lgmax + 1.0)
                   - (lgmax + 1.0) / tau_safe * Theta_lmax)

    # ---- Baryons ----
    d_v_b = -calH * v_b + k_arr * Psi
    d_delta_b = -k_arr * v_b - 3.0 * Phi_dot

    # ---- CDM ----
    d_delta_c = -(k_arr * v_c + 3.0 * Phi_dot)
    d_v_c = -calH * v_c + k_arr * Psi

    # ---- Neutrino hierarchy (vectorized) ----
    dN_0 = -k_arr * N_1 - Phi_dot
    dN_1 = k_arr / 3.0 * (N_0 - 2.0 * N_2 + Psi)

    # l=2 to l_nu_max-1: VECTORIZED
    if lnmax > 2:
        l_nu_vals = mx.arange(2, lnmax).astype(mx.float32)
        inv_2lp1_nu = 1.0 / (2.0 * l_nu_vals + 1.0)
        n_lm1 = y[:, _n_start + 1:_n_start + lnmax - 1]
        n_lp1 = y[:, _n_start + 3:_n_start + lnmax + 1]
        dN_mid = k_arr[:, None] * inv_2lp1_nu[None, :] * (
            l_nu_vals[None, :] * n_lm1 - (l_nu_vals[None, :] + 1.0) * n_lp1
        )
    else:
        dN_mid = None

    # l=l_nu_max: truncation
    N_prev = y[:, _n_start + lnmax - 1]
    N_last = y[:, _n_start + lnmax]
    dN_lmax = (k_arr * N_prev * lnmax / (2.0 * lnmax + 1.0)
               - (lnmax + 1.0) / tau_safe * N_last)

    # ---- Stack result ----
    parts = [
        Phi_dot[:, None],
        d_delta_b[:, None],
        d_v_b[:, None],
        d_delta_c[:, None],
        d_v_c[:, None],
        dTheta_0[:, None],
        dTheta_1[:, None],
        dTheta_2[:, None],
    ]
    if dTheta_mid is not None:
        parts.append(dTheta_mid)
    parts.append(dTheta_lmax[:, None])

    parts.append(dN_0[:, None])
    parts.append(dN_1[:, None])
    if dN_mid is not None:
        parts.append(dN_mid)
    parts.append(dN_lmax[:, None])

    return mx.concatenate(parts, axis=1)


def _deriv_streaming_tca(y, k_arr, calH, a, R, tau_val,
                         Og, On, Ob, Oc, H0, lgmax, lnmax, _n_start):
    """
    Vectorized streaming + gravity RHS for the TCA phase.
    Theta_l = 0 for l >= 2 (tight coupling).
    """
    Phi = y[:, IDX_PHI]
    delta_b = y[:, IDX_DELTA_B]
    v_b = y[:, IDX_V_B]
    delta_c = y[:, IDX_DELTA_C]
    v_c = y[:, IDX_V_C]
    Theta_0 = y[:, idx_theta(0)]
    Theta_1 = y[:, idx_theta(1)]

    N_0 = y[:, _n_start]
    N_1 = y[:, _n_start + 1]
    N_2 = y[:, _n_start + 2]

    H02 = H0 * H0
    inv_a2 = 1.0 / (a * a)
    k2 = k_arr * k_arr

    # Psi (no photon anisotropic stress in TCA)
    Psi = Phi - 12.0 * H02 * On * N_2 * inv_a2 / k2

    # Phi_dot
    S = (Og * inv_a2 * 4.0 * Theta_0
         + On * inv_a2 * 4.0 * N_0
         + Ob / a * delta_b
         + Oc / a * delta_c)
    Phi_dot = -calH * Phi - k2 * Phi / (3.0 * calH) - 0.5 * H02 / calH * S

    # Photon TCA
    dTheta_0 = -k_arr * Theta_1 - Phi_dot
    dTheta_1 = (k_arr / 3.0 * (Theta_0 + Psi) - calH * R * Theta_1) / (1.0 + R)

    # Baryons in TCA
    d_v_b = 3.0 * dTheta_1
    d_delta_b = -k_arr * v_b - 3.0 * Phi_dot

    # CDM
    d_delta_c = -(k_arr * v_c + 3.0 * Phi_dot)
    d_v_c = -calH * v_c + k_arr * Psi

    # Neutrino hierarchy (vectorized)
    dN_0 = -k_arr * N_1 - Phi_dot
    dN_1 = k_arr / 3.0 * (N_0 - 2.0 * N_2 + Psi)

    tau_safe = mx.maximum(tau_val, mx.array(1e-10))

    if lnmax > 2:
        l_nu_vals = mx.arange(2, lnmax).astype(mx.float32)
        inv_2lp1_nu = 1.0 / (2.0 * l_nu_vals + 1.0)
        n_lm1 = y[:, _n_start + 1:_n_start + lnmax - 1]
        n_lp1 = y[:, _n_start + 3:_n_start + lnmax + 1]
        dN_mid = k_arr[:, None] * inv_2lp1_nu[None, :] * (
            l_nu_vals[None, :] * n_lm1 - (l_nu_vals[None, :] + 1.0) * n_lp1
        )
    else:
        dN_mid = None

    N_prev = y[:, _n_start + lnmax - 1]
    N_last = y[:, _n_start + lnmax]
    dN_lmax = (k_arr * N_prev * lnmax / (2.0 * lnmax + 1.0)
               - (lnmax + 1.0) / tau_safe * N_last)

    # Stack — zeros for Theta_2..Theta_lgmax
    nvar = n_var_total(lgmax, lnmax)
    n_theta_zero = lgmax - 1  # Theta_2 through Theta_lgmax

    parts = [
        Phi_dot[:, None],
        d_delta_b[:, None],
        d_v_b[:, None],
        d_delta_c[:, None],
        d_v_c[:, None],
        dTheta_0[:, None],
        dTheta_1[:, None],
        mx.zeros((y.shape[0], n_theta_zero)),  # Theta_2..Theta_lgmax
        dN_0[:, None],
        dN_1[:, None],
    ]
    if dN_mid is not None:
        parts.append(dN_mid)
    parts.append(dN_lmax[:, None])

    return mx.concatenate(parts, axis=1)


# ============================================================================
# Vectorized collision operator
# ============================================================================

def _apply_collision_vectorized(y, k_arr, R, kd, lgmax, _n_start, dt):
    """
    Thomson collision operator applied analytically (no Python scalar conversions).

    Parameters are all MLX scalars, not Python floats, so this can be compiled.
    """
    # Theta_l damping for l >= 2
    exp_kd_dt_t2 = mx.exp(0.9 * kd * dt)   # Theta_2: (9/10)*|kd| damping
    exp_kd_dt_hi = mx.exp(kd * dt)          # Theta_l (l>=3): |kd| damping

    # Build damping mask: multiply each column by its factor
    # Columns: [Phi, db, vb, dc, vc, T0, T1, T2, T3...Tlmax, N0...Nlmax]
    # No damping for Phi..T1 (first 7 columns)
    # T2: exp(0.9*kd*dt)
    # T3..Tlmax: exp(kd*dt)
    # N0..Nlmax: no damping

    parts = [
        y[:, :idx_theta(2)],                                          # Phi..T1: no change
        y[:, idx_theta(2):idx_theta(2)+1] * exp_kd_dt_t2,            # T2
    ]
    if lgmax >= 3:
        parts.append(y[:, idx_theta(3):idx_theta(lgmax)+1] * exp_kd_dt_hi)  # T3..Tlmax
    parts.append(y[:, _n_start:])                                     # neutrinos: no change
    y = mx.concatenate(parts, axis=1)

    # (Theta_1, v_b) coupled system: exact eigenmode solution
    R_safe = mx.maximum(R, mx.array(1e-10))
    lambda_2 = kd * (1.0 + 1.0 / R_safe)  # negative (decaying)
    exp_lam2_dt = mx.exp(lambda_2 * dt)

    Theta_1 = y[:, idx_theta(1)]
    v_b = y[:, IDX_V_B]

    slip = 3.0 * Theta_1 - v_b
    c2 = R_safe / (3.0 * (R_safe + 1.0)) * slip
    c1 = Theta_1 - c2

    Theta_1_new = c1 + c2 * exp_lam2_dt
    v_b_new = 3.0 * c1 - 3.0 / R_safe * c2 * exp_lam2_dt

    # Replace Theta_1 and v_b
    y = mx.concatenate([
        y[:, :IDX_V_B],
        v_b_new[:, None],
        y[:, IDX_V_B+1:idx_theta(1)],
        Theta_1_new[:, None],
        y[:, idx_theta(1)+1:]
    ], axis=1)

    return y


# ============================================================================
# IMEX RK4 streaming step (compilable — no Python branching on data)
# ============================================================================

def _imex_rk4_step_full(y, k_arr, calH, a, R, tau_val,
                        Og, On, Ob, Oc, H0, lgmax, lnmax, _n_start, dtau):
    """IMEX RK4 for the streaming part (full hierarchy). Compilable."""
    lam = _phi_lambda(k_arr, calH)
    exp_lam_dt = mx.exp(lam * dtau)
    exp_lam_half = mx.exp(lam * 0.5 * dtau)

    def get_F(state):
        return _phi_source_fast(state, k_arr, calH, a, Og, On, Ob, Oc, H0, _n_start)

    def rhs(state):
        return _deriv_streaming_full(state, k_arr, calH, a, R, tau_val,
                                     Og, On, Ob, Oc, H0, lgmax, lnmax, _n_start)

    Phi_1 = y[:, IDX_PHI]
    F1 = get_F(y)
    k1 = rhs(y)

    Phi_2 = exp_lam_half * Phi_1 + _phi_int_factor(lam, 0.5 * dtau) * F1
    y2 = y + 0.5 * dtau * k1
    y2 = _set_phi(y2, Phi_2)
    F2 = get_F(y2)
    k2 = rhs(y2)

    Phi_3 = exp_lam_half * Phi_1 + _phi_int_factor(lam, 0.5 * dtau) * F2
    y3 = y + 0.5 * dtau * k2
    y3 = _set_phi(y3, Phi_3)
    F3 = get_F(y3)
    k3 = rhs(y3)

    Phi_4 = exp_lam_dt * Phi_1 + _phi_int_factor(lam, dtau) * F3
    y4 = y + dtau * k3
    y4 = _set_phi(y4, Phi_4)
    F4 = get_F(y4)
    k4 = rhs(y4)

    y_new = y + (dtau / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)

    F_mid = 0.5 * (F2 + F3)
    Phi_new = (exp_lam_dt * Phi_1
               + (dtau / 6.0) * (exp_lam_dt * F1
                                  + 4.0 * exp_lam_half * F_mid
                                  + F4))
    y_new = _set_phi(y_new, Phi_new)
    return y_new


def _imex_rk4_step_tca(y, k_arr, calH, a, R, tau_val,
                       Og, On, Ob, Oc, H0, lgmax, lnmax, _n_start, dtau):
    """IMEX RK4 for the streaming part (TCA). Compilable."""
    lam = _phi_lambda(k_arr, calH)
    exp_lam_dt = mx.exp(lam * dtau)
    exp_lam_half = mx.exp(lam * 0.5 * dtau)

    def get_F(state):
        return _phi_source_fast(state, k_arr, calH, a, Og, On, Ob, Oc, H0, _n_start)

    def rhs(state):
        return _deriv_streaming_tca(state, k_arr, calH, a, R, tau_val,
                                    Og, On, Ob, Oc, H0, lgmax, lnmax, _n_start)

    Phi_1 = y[:, IDX_PHI]
    F1 = get_F(y)
    k1 = rhs(y)

    Phi_2 = exp_lam_half * Phi_1 + _phi_int_factor(lam, 0.5 * dtau) * F1
    y2 = y + 0.5 * dtau * k1
    y2 = _set_phi(y2, Phi_2)
    F2 = get_F(y2)
    k2 = rhs(y2)

    Phi_3 = exp_lam_half * Phi_1 + _phi_int_factor(lam, 0.5 * dtau) * F2
    y3 = y + 0.5 * dtau * k2
    y3 = _set_phi(y3, Phi_3)
    F3 = get_F(y3)
    k3 = rhs(y3)

    Phi_4 = exp_lam_dt * Phi_1 + _phi_int_factor(lam, dtau) * F3
    y4 = y + dtau * k3
    y4 = _set_phi(y4, Phi_4)
    F4 = get_F(y4)
    k4 = rhs(y4)

    y_new = y + (dtau / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)

    F_mid = 0.5 * (F2 + F3)
    Phi_new = (exp_lam_dt * Phi_1
               + (dtau / 6.0) * (exp_lam_dt * F1
                                  + 4.0 * exp_lam_half * F_mid
                                  + F4))
    y_new = _set_phi(y_new, Phi_new)

    # Enforce v_b = 3 * Theta_1 in TCA
    Theta_1_new = y_new[:, idx_theta(1)]
    v_b_locked = 3.0 * Theta_1_new
    y_new = mx.concatenate([
        y_new[:, :IDX_V_B],
        v_b_locked[:, None],
        y_new[:, IDX_V_B + 1:]
    ], axis=1)

    return y_new


def _phi_source_fast(state, k_arr, calH, a, Og, On, Ob, Oc, H0, _n_start):
    """Phi source function — no branching, compilable."""
    H02 = H0 * H0
    inv_a2 = 1.0 / (a * a)
    S = (Og * inv_a2 * 4.0 * state[:, idx_theta(0)]
         + On * inv_a2 * 4.0 * state[:, _n_start]
         + Ob / a * state[:, IDX_DELTA_B]
         + Oc / a * state[:, IDX_DELTA_C])
    return -0.5 * H02 / calH * S


# ============================================================================
# Full Strang step for full hierarchy (compilable)
# ============================================================================

def _strang_step_full(y, k_arr, calH, a, R, tau_val, kd,
                      Og, On, Ob, Oc, H0, lgmax, lnmax, _n_start, dtau):
    """Strang splitting step for full hierarchy. All args are MLX arrays."""
    # Half collision
    y = _apply_collision_vectorized(y, k_arr, R, kd, lgmax, _n_start, 0.5 * dtau)
    # Full streaming
    y = _imex_rk4_step_full(y, k_arr, calH, a, R, tau_val,
                            Og, On, Ob, Oc, H0, lgmax, lnmax, _n_start, dtau)
    # Half collision
    y = _apply_collision_vectorized(y, k_arr, R, kd, lgmax, _n_start, 0.5 * dtau)
    return y


# ============================================================================
# Time grid — reduced step count
# ============================================================================

def build_tau_grid_fast(bg, k_max, N_early=200, N_late=400):
    """
    Build a coarser grid than the HiRes solver.

    The HiRes solver uses N_early=300 + N_late=800 ~ 1100 steps.
    We use N_early=200 + N_late=400 ~ 600 steps.

    For RK4 with Strang splitting, the CFL is governed by:
      dtau < C / (k_max * c_s) where C ~ 2.8 for RK4
    We use 15 steps/period instead of 25 (still well within stability).
    """
    tau_rec = bg.tau_rec
    tau_init = bg.tau_grid[1]
    tau_early_end = min(10.0, 0.3 * tau_rec)

    tau_early = np.geomspace(tau_init, tau_early_end, N_early)

    cs_approx = 1.0 / np.sqrt(3.0)
    dtau_acoustic = 2.0 * np.pi / (k_max * cs_approx) / 15.0
    N_needed = int((tau_rec - tau_early_end) / dtau_acoustic) + 10
    N_late_actual = max(N_late, N_needed)

    tau_late = np.linspace(tau_early_end, tau_rec, N_late_actual)
    tau_grid = np.unique(np.concatenate([tau_early, tau_late]))
    return tau_grid


# ============================================================================
# Initial conditions
# ============================================================================

def adiabatic_ic_fast(k_arr_np, bg, l_gamma_max, l_nu_max=L_NU_MAX):
    """Same as perturbations_hires.adiabatic_ic_hires."""
    N_k = len(k_arr_np)
    tau_init = bg.tau_grid[1]
    nvar = n_var_total(l_gamma_max, l_nu_max)
    _n_start = idx_n_start(l_gamma_max)

    y0 = np.zeros((N_k, nvar), dtype=np.float32)
    y0[:, IDX_PHI] = 1.0
    y0[:, IDX_DELTA_B] = -1.5
    y0[:, IDX_V_B] = k_arr_np * tau_init / 6.0
    y0[:, IDX_DELTA_C] = -1.5
    y0[:, IDX_V_C] = k_arr_np * tau_init / 6.0

    y0[:, idx_theta(0)] = -0.5
    y0[:, idx_theta(1)] = k_arr_np * tau_init / 18.0

    y0[:, _n_start + 0] = -0.5
    y0[:, _n_start + 1] = k_arr_np * tau_init / 18.0
    if l_nu_max >= 2:
        y0[:, _n_start + 2] = (k_arr_np * tau_init) ** 2 / 60.0
    for l in range(3, min(l_nu_max + 1, 6)):
        prod_val = 1.0
        for j in range(1, l + 1):
            prod_val *= (2 * j + 1)
        y0[:, _n_start + l] = (k_arr_np * tau_init) ** l / prod_val

    return mx.array(y0)


# ============================================================================
# Fast solver class
# ============================================================================

class FastBoltzmannSolver:
    """
    Optimized Boltzmann solver. Same physics as HiResBoltzmannSolver,
    but 5-10x faster through:
      1. Vectorized hierarchy (no Python loops over l)
      2. mx.compile for step functions
      3. Reduced eval frequency
      4. Reduced grid density
      5. All background quantities pre-packed as MLX arrays (no per-step conversion)
    """

    def __init__(self, bg, k_arr_Mpc, l_gamma_max=25, l_nu_max=L_NU_MAX,
                 tca_threshold=50.0, grid_density='fast'):
        self.bg = bg
        self.k_arr_np = np.asarray(k_arr_Mpc, dtype=np.float32)
        self.N_k = len(self.k_arr_np)
        self.k_arr = mx.array(self.k_arr_np)
        self.l_gamma_max = l_gamma_max
        self.l_nu_max = l_nu_max
        self.tca_threshold = tca_threshold
        self.n_var = n_var_total(l_gamma_max, l_nu_max)
        self._n_start = idx_n_start(l_gamma_max)

        k_max = float(self.k_arr_np[-1])

        if grid_density == 'fast':
            self.tau_grid = build_tau_grid_fast(bg, k_max)
        elif grid_density == 'hires':
            from .perturbations_hires import build_tau_grid_hires
            self.tau_grid = build_tau_grid_hires(bg, k_max)
        else:
            self.tau_grid = build_tau_grid_fast(bg, k_max,
                                                N_early=grid_density[0],
                                                N_late=grid_density[1])

        self.Omega_gamma = _OMEGA_GAMMA
        self.Omega_nu = _OMEGA_NU
        self.Omega_b = _OMEGA_B
        self.Omega_c = float(bg.Omega_cdm)

        print(f"[Fast] N_k={self.N_k}, l_gamma_max={l_gamma_max}, "
              f"l_nu_max={l_nu_max}, N_var={self.n_var}")
        print(f"[Fast] Grid: {len(self.tau_grid)} steps "
              f"[{self.tau_grid[0]:.2e}, {self.tau_grid[-1]:.1f}] Mpc")

    def solve(self):
        t0 = time.time()
        bg = self.bg
        k_arr = self.k_arr
        tau_grid = self.tau_grid
        lgmax = self.l_gamma_max
        lnmax = self.l_nu_max
        _ns = self._n_start

        y = adiabatic_ic_fast(self.k_arr_np, bg, lgmax, lnmax)

        # Precompute ALL background quantities as MLX arrays
        calH_np = bg.calH_at_tau(tau_grid).astype(np.float32)
        R_np = bg.R_at_tau(tau_grid).astype(np.float32)
        a_np = bg.a_at_tau(tau_grid).astype(np.float32)
        kappa_dot_np = bg.kappa_dot_at_tau(tau_grid).astype(np.float32)

        # Pack into MLX arrays — one per time step (avoid per-step conversion)
        calH_all = mx.array(calH_np)
        R_all = mx.array(R_np)
        a_all = mx.array(a_np)
        kd_all = mx.array(kappa_dot_np)
        tau_all = mx.array(tau_grid.astype(np.float32))

        # Precompute dtau for each step
        dtau_np = np.diff(tau_grid).astype(np.float32)
        dtau_all = mx.array(dtau_np)

        _Og = self.Omega_gamma
        _On = self.Omega_nu
        _Ob = self.Omega_b
        _Oc = self.Omega_c
        _H0 = _H0_MPC
        k_max = float(self.k_arr_np[-1])

        # Determine TCA switch point
        abs_kd = np.abs(kappa_dot_np)
        tca_switch_idx = len(tau_grid) - 1
        for i in range(len(tau_grid)):
            if abs_kd[i] / k_max < self.tca_threshold:
                tca_switch_idx = i
                break

        tau_switch = tau_grid[tca_switch_idx]
        z_switch = 1.0 / float(bg.a_at_tau(np.array([tau_switch]))[0]) - 1.0
        print(f"[Fast] TCA -> full at tau={tau_switch:.1f} (z={z_switch:.0f}, "
              f"step {tca_switch_idx}/{len(tau_grid)})")

        # ---- Compile step functions ----
        t_compile_start = time.time()

        # TCA step: just the streaming RK4 (no collision splitting in TCA)
        def tca_step(y, calH_i, a_i, R_i, tau_i, dtau_i):
            return _imex_rk4_step_tca(y, k_arr, calH_i, a_i, R_i, tau_i,
                                      _Og, _On, _Ob, _Oc, _H0,
                                      lgmax, lnmax, _ns, dtau_i)

        # Full hierarchy step with Strang splitting
        def full_step(y, calH_i, a_i, R_i, tau_i, kd_i, dtau_i):
            return _strang_step_full(y, k_arr, calH_i, a_i, R_i, tau_i, kd_i,
                                     _Og, _On, _Ob, _Oc, _H0,
                                     lgmax, lnmax, _ns, dtau_i)

        compiled_tca_step = mx.compile(tca_step)
        compiled_full_step = mx.compile(full_step)

        # Warmup compile with a dummy step
        y_dummy = mx.zeros_like(y)
        _ = compiled_tca_step(y_dummy, calH_all[0], a_all[0], R_all[0],
                              tau_all[0], dtau_all[0])
        mx.eval(_)
        if tca_switch_idx < len(tau_grid) - 1:
            _ = compiled_full_step(y_dummy, calH_all[0], a_all[0], R_all[0],
                                   tau_all[0], kd_all[0], dtau_all[0])
            mx.eval(_)
        t_compile = time.time() - t_compile_start
        print(f"[Fast] Compile warmup: {t_compile:.2f}s")

        # ---- TCA phase ----
        t_tca_start = time.time()
        # Determine eval points: divide TCA into ~2 chunks
        tca_eval_interval = max(tca_switch_idx // 2, 1)

        for i in range(tca_switch_idx):
            y = compiled_tca_step(y, calH_all[i], a_all[i], R_all[i],
                                  tau_all[i], dtau_all[i])
            if (i + 1) % tca_eval_interval == 0:
                mx.eval(y)

        mx.eval(y)
        t_tca = time.time() - t_tca_start
        print(f"[Fast] TCA phase: {t_tca:.2f}s ({tca_switch_idx} steps)")

        # Seed Theta_2 at transition
        kd_switch = abs_kd[tca_switch_idx]
        if kd_switch > 1e-10:
            Theta_1_sw = y[:, idx_theta(1)]
            Theta_2_seed = (4.0 / 9.0) * k_arr * Theta_1_sw / kd_switch
            _idx_t2 = idx_theta(2)
            y = mx.concatenate([
                y[:, :_idx_t2],
                Theta_2_seed[:, None],
                y[:, _idx_t2 + 1:]
            ], axis=1)
            mx.eval(y)

        # ---- Full hierarchy phase ----
        t_full_start = time.time()
        n_full = len(tau_grid) - 1 - tca_switch_idx
        full_eval_interval = max(n_full // 3, 1)

        for i in range(tca_switch_idx, len(tau_grid) - 1):
            y = compiled_full_step(y, calH_all[i], a_all[i], R_all[i],
                                   tau_all[i], kd_all[i], dtau_all[i])
            step_in_phase = i - tca_switch_idx
            if (step_in_phase + 1) % full_eval_interval == 0:
                mx.eval(y)

        mx.eval(y)
        t_full = time.time() - t_full_start
        t_total = time.time() - t0
        print(f"[Fast] Full hierarchy: {t_full:.2f}s ({n_full} steps)")
        print(f"[Fast] Total: {t_total:.2f}s")

        return FastResult(y=y, k_arr=self.k_arr_np, bg=bg,
                          l_gamma_max=lgmax, l_nu_max=lnmax,
                          Omega_gamma=self.Omega_gamma,
                          Omega_nu=self.Omega_nu, H0=_H0)


# ============================================================================
# Result class (compatible with HiResResult)
# ============================================================================

class FastResult:
    """Result container, compatible with HiResResult."""

    def __init__(self, y, k_arr, bg, l_gamma_max, l_nu_max,
                 Omega_gamma, Omega_nu, H0):
        self.y = y
        self.k_arr = k_arr
        self.bg = bg
        self.l_gamma_max = l_gamma_max
        self.l_nu_max = l_nu_max
        self.Omega_gamma = Omega_gamma
        self.Omega_nu = Omega_nu
        self.H0 = H0
        self.N_k = len(k_arr)

    def source_at_recombination(self):
        """Same as HiResResult.source_at_recombination."""
        y_np = np.array(self.y)
        _n_start = idx_n_start(self.l_gamma_max)

        Phi = y_np[:, IDX_PHI]
        Theta_0 = y_np[:, idx_theta(0)]
        v_b = y_np[:, IDX_V_B]
        Theta_2 = y_np[:, idx_theta(2)] if self.l_gamma_max >= 2 else np.zeros_like(Phi)
        N_2 = y_np[:, _n_start + 2] if self.l_nu_max >= 2 else np.zeros_like(Phi)

        H02 = self.H0 ** 2
        a_rec = float(self.bg.a_at_tau(np.array([self.bg.tau_rec]))[0])
        Psi = (Phi
               - 12.0 * H02 * self.Omega_nu * N_2 / (a_rec ** 2 * self.k_arr ** 2)
               - 12.0 * H02 * self.Omega_gamma * Theta_2 / (a_rec ** 2 * self.k_arr ** 2))

        silk = np.exp(-(self.k_arr / self.bg.k_D) ** 2)
        return (Theta_0 * silk, Phi * silk, Psi * silk,
                v_b * silk, Theta_2 * silk)

    def sw_source(self):
        """Sachs-Wolfe source: Theta_0 + Psi (silk-damped)."""
        Theta_0, Phi, Psi, v_b, Theta_2 = self.source_at_recombination()
        return Theta_0 + Psi

    def photon_multipoles(self):
        """Return all photon multipoles at tau_rec (no damping)."""
        y_np = np.array(self.y)
        multipoles = {}
        for l in range(self.l_gamma_max + 1):
            multipoles[l] = y_np[:, idx_theta(l)]
        return multipoles


# ============================================================================
# Benchmark & validation
# ============================================================================

def benchmark_fast_solver():
    """
    Benchmark the fast solver against HiRes and validate accuracy.

    NOTE on compile warmup:
      The first mx.compile() call in a fresh Python process takes ~3-4s
      because MLX compiles Metal shaders. This is a one-time cost:
      subsequent calls (including in MCMC chains) take ~0.04s.
      Integration time (the actual ODE solve) is ~0.5s regardless.
    """
    from .background import Background
    from .perturbations_hires import HiResBoltzmannSolver

    print("=" * 70)
    print("BENCHMARK: FastBoltzmannSolver vs HiResBoltzmannSolver")
    print("=" * 70)

    bg = Background(khronon=False, recombination='peebles')
    bg.solve()

    k_arr = np.geomspace(5e-4, 0.35, 500).astype(np.float32)

    # ---- Fast solver (1st call — includes MLX shader compilation) ----
    print("\n--- Fast Solver (1st call) ---")
    solver_fast = FastBoltzmannSolver(bg, k_arr)
    t0 = time.time()
    result_fast = solver_fast.solve()
    t_fast_1st = time.time() - t0

    # ---- Fast solver (2nd call — compile cached) ----
    print("\n--- Fast Solver (2nd call, cached) ---")
    solver_fast2 = FastBoltzmannSolver(bg, k_arr)
    t0 = time.time()
    result_fast2 = solver_fast2.solve()
    t_fast_2nd = time.time() - t0

    # ---- HiRes solver (reference) ----
    print("\n--- HiRes Solver (reference) ---")
    solver_hires = HiResBoltzmannSolver(bg, k_arr)
    t0 = time.time()
    result_hires = solver_hires.solve()
    t_hires = time.time() - t0

    # ---- Compare ----
    print("\n--- Comparison ---")
    Theta_0_f, Phi_f, Psi_f, v_b_f, T2_f = result_fast.source_at_recombination()
    Theta_0_h, Phi_h, Psi_h, v_b_h, T2_h = result_hires.source_at_recombination()

    sw_fast = Theta_0_f + Psi_f
    sw_hires = Theta_0_h + Psi_h

    corr = np.corrcoef(sw_fast, sw_hires)[0, 1]
    max_diff = np.max(np.abs(sw_fast - sw_hires))
    rms_diff = np.sqrt(np.mean((sw_fast - sw_hires)**2))
    rms_ref = np.sqrt(np.mean(sw_hires**2))

    print(f"  Correlation:    {corr:.6f}")
    print(f"  Max abs diff:   {max_diff:.6f}")
    print(f"  RMS diff:       {rms_diff:.6f}")
    print(f"  RMS relative:   {rms_diff/rms_ref*100:.2f}%")

    print(f"\n--- Timing ---")
    print(f"  HiRes:               {t_hires:.2f}s")
    print(f"  Fast (1st call):     {t_fast_1st:.2f}s  (includes MLX shader warmup)")
    print(f"  Fast (2nd call):     {t_fast_2nd:.2f}s  (compile cached)")
    print(f"  Speedup (cached):    {t_hires/t_fast_2nd:.1f}x")

    # Check pass criteria using the cached (amortized) time
    if corr > 0.999 and t_fast_2nd < 1.0:
        print(f"\n  TARGET MET: <1s and corr>0.999")
    elif corr > 0.999:
        print(f"\n  ACCURACY OK (corr>0.999) but speed target not met ({t_fast_2nd:.2f}s > 1s)")
    else:
        print(f"\n  WARNING: Accuracy below threshold (corr={corr:.6f})")

    print("=" * 70)
    return result_fast, result_hires


if __name__ == '__main__':
    benchmark_fast_solver()
