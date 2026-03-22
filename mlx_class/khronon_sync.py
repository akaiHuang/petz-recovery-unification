"""
khronon_sync.py -- Khronon dark matter CMB solver using synchronous gauge backend.

Implements the Khronon DBI modification to the standard synchronous gauge
Boltzmann equations. The key physics:

  - Standard CDM: delta_c' = -h'/2, theta_c = 0 (gauge choice)
  - Khronon DBI:  delta_K' = -theta_K - h'/2 - 3*calH*c_s^2*delta_K
                   theta_K' = -calH*theta_K + c_s^2*k^2*delta_K + k^2*alpha_sync
  - where alpha_sync = (h' + 6*eta') / (2*k^2) is the sync gauge shift

  - c_s^2(k, a) = alpha_DBI * (k/k_J)^2 / (1 + (k/k_J)^2)
  - k_J(a) = sqrt(1.5 * Omega_K * H0^2 / a)
  - alpha_DBI = Q_0/2 = 0.5 (ghost condensation)

NOTE: In synchronous gauge, the CDM rest frame sets theta_c = 0. For Khronon,
we break this gauge choice since c_s^2 != 0 gives the Khronon field its own
velocity perturbation. The state vector is extended to include theta_K.

Usage:
  python -m mlx_class.khronon_sync [--N_k 100] [--alpha_dbi 0.5]

Author: Sheng-Kai Huang, 2026
"""
import os
os.environ['MPLCONFIGDIR'] = '/tmp/mpl_config'

import numpy as np
import time
import sys

from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d, CubicSpline

from .background import (
    Background, A_s, n_s, k_pivot, T_CMB,
    Omega_r as _OMEGA_R,
    Omega_b as _OMEGA_B,
    Omega_c as _OMEGA_C,
    Omega_m as _OMEGA_M,
    H0_Mpc as _H0_MPC,
    z_rec, a_rec,
)
from .perturbations_neutrino import (
    _f_nu, _OMEGA_GAMMA, _OMEGA_NU,
)
from .perturbations_sync import (
    L_GAMMA_MAX, L_NU_MAX_SYNC,
    IDX_ETA, IDX_DELTA_C, IDX_DELTA_B, IDX_THETA_B, IDX_FG_START,
    idx_fg, idx_fn_start, idx_fn, n_var_sync,
    adiabatic_ic_sync, gauge_transform, diagnose_metric,
)
from .solver_sync import build_snapshot_grid


# ============================================================================
# Khronon DBI parameters
# ============================================================================

Q_0 = 1.0
ALPHA_DBI_DEFAULT = Q_0 / 2.0  # = 0.5, ghost condensation bare value


def khronon_cs2(k, a, Omega_K, H0_Mpc, alpha_dbi=ALPHA_DBI_DEFAULT):
    """
    DBI effective sound speed squared for the Khronon field.

    c_s^2(k, a) = alpha_DBI * x^2 / (1 + x^2)
    where x = k / k_J(a),  k_J(a) = sqrt(1.5 * Omega_K * H0^2 / a)
    """
    k_J_sq = 1.5 * Omega_K * H0_Mpc**2 / a
    x_sq = k**2 / k_J_sq
    return alpha_dbi * x_sq / (1.0 + x_sq)


# ============================================================================
# Khronon sync gauge state vector
# ============================================================================
# Extended state: standard sync variables + theta_K appended at the end.
#
# Standard: [eta, delta_c, delta_b, theta_b, F_g0..F_g{lg}, F_n0..F_n{ln}]
# Khronon:  [eta, delta_K, delta_b, theta_b, F_g0..F_g{lg}, F_n0..F_n{ln}, theta_K]
#
# Note: IDX_DELTA_C (=1) now holds delta_K (Khronon density contrast).
# theta_K is at position n_var_sync(lg_max, ln_max) (appended).

IDX_DELTA_K = IDX_DELTA_C  # reuse slot, but with different evolution equation
IDX_THETA_K_OFFSET = 0     # set dynamically as n_var_sync(lg, ln)


def n_var_khronon(lg_max, ln_max):
    """Total variables: standard sync + theta_K."""
    return n_var_sync(lg_max, ln_max) + 1


def idx_theta_K(lg_max, ln_max):
    """Index of theta_K in the extended state vector."""
    return n_var_sync(lg_max, ln_max)


# ============================================================================
# Khronon sync gauge RHS
# ============================================================================

def make_khronon_sync_rhs(k, bg, lg_max=L_GAMMA_MAX, ln_max=L_NU_MAX_SYNC,
                          alpha_dbi=ALPHA_DBI_DEFAULT):
    """
    Build the RHS function for Khronon DBI in synchronous gauge.

    Modifications from standard sync gauge:
    1. delta_K evolves with pressure terms (not just -h'/2)
    2. theta_K is evolved (not zero as in CDM gauge choice)
    3. Poisson equation includes Khronon velocity contribution

    Parameters
    ----------
    k : float, wavenumber in Mpc^-1
    bg : Background object (solved)
    lg_max, ln_max : int, hierarchy truncation
    alpha_dbi : float, DBI strength parameter

    Returns
    -------
    rhs : callable f(tau, y) -> dy/dtau
    nvar : int, length of extended state vector
    """
    k2 = k * k
    H02 = _H0_MPC ** 2
    Og, On, Ob = _OMEGA_GAMMA, _OMEGA_NU, _OMEGA_B
    Ok = _OMEGA_C  # Khronon replaces CDM
    nvar_std = n_var_sync(lg_max, ln_max)
    nvar = nvar_std + 1  # +1 for theta_K
    fn_s = idx_fn_start(lg_max)
    i_thetaK = nvar_std  # index of theta_K

    # Pre-tabulate background
    _tau = bg.tau_grid.copy()
    _calH = bg.calH_grid.copy()
    _a = bg.a_grid.copy()
    _R = bg.R_grid.copy()
    _kd = bg.kappa_dot_grid.copy()

    def rhs(tau, y):
        calH = np.interp(tau, _tau, _calH)
        a = np.interp(tau, _tau, _a)
        R = np.interp(tau, _tau, _R)
        kappa_dot = np.interp(tau, _tau, _kd)
        abs_kd = abs(kappa_dot)

        ia = 1.0 / a
        ia2 = ia * ia
        tau_safe = max(tau, 1e-10)

        # --- Extract state ---
        eta = y[IDX_ETA]
        delta_K = y[IDX_DELTA_K]
        delta_b = y[IDX_DELTA_B]
        theta_b = y[IDX_THETA_B]
        theta_K = y[i_thetaK]

        Fg = y[IDX_FG_START: IDX_FG_START + lg_max + 1]
        Fn = y[fn_s: fn_s + ln_max + 1]

        delta_g = Fg[0]
        delta_n = Fn[0]
        theta_g = 0.75 * k * Fg[1]
        theta_n = 0.75 * k * Fn[1]

        # --- Khronon sound speed ---
        cs2 = khronon_cs2(k, a, Ok, _H0_MPC, alpha_dbi)

        # --- Diagnose h' from 00-constraint ---
        # Uses delta_K in place of delta_c (same slot, same weight in Poisson)
        src00 = (Og * ia2 * delta_g + On * ia2 * delta_n
                 + Ob * ia * delta_b + Ok * ia * delta_K)
        h_prime = (2.0 / calH) * (k2 * eta + 1.5 * H02 * src00)

        # --- Diagnose eta' from 0i-constraint ---
        # Khronon contributes theta_K (CDM had theta_c = 0)
        src0i = ((4.0 / 3.0) * Og * ia2 * theta_g
                 + (4.0 / 3.0) * On * ia2 * theta_n
                 + Ob * ia * theta_b
                 + Ok * ia * theta_K)
        eta_prime = 1.5 * H02 / k2 * src0i

        # --- Khronon evolution (GDM equations) ---
        # Continuity: delta_K' = -theta_K - h'/2 - 3*calH*c_s^2*delta_K
        d_delta_K = -theta_K - 0.5 * h_prime - 3.0 * calH * cs2 * delta_K

        # Euler: theta_K' = -calH*theta_K + c_s^2*k^2*delta_K
        d_theta_K = -calH * theta_K + cs2 * k2 * delta_K

        # --- Baryons ---
        d_delta_b = -theta_b - 0.5 * h_prime
        d_theta_b = -calH * theta_b + abs_kd / R * (theta_g - theta_b)

        # --- Photon hierarchy (same as standard sync) ---
        dFg = np.zeros(lg_max + 1)
        dFg[0] = -k * Fg[1] - (2.0 / 3.0) * h_prime
        F2g = Fg[2] if lg_max >= 2 else 0.0
        dFg[1] = ((k / 3.0) * (Fg[0] - 2.0 * F2g)
                  + abs_kd * (-Fg[1] + 4.0 * theta_b / (3.0 * k)))
        if lg_max >= 2:
            F3g = Fg[3] if lg_max >= 3 else 0.0
            dFg[2] = ((k / 5.0) * (2.0 * Fg[1] - 3.0 * F3g)
                      + (4.0 / 15.0) * h_prime + (8.0 / 15.0) * eta_prime
                      - (9.0 / 10.0) * abs_kd * Fg[2])
        for ell in range(3, lg_max):
            dFg[ell] = (k / (2.0 * ell + 1.0)
                        * (ell * Fg[ell - 1] - (ell + 1) * Fg[ell + 1])
                        - abs_kd * Fg[ell])
        if lg_max >= 3:
            dFg[lg_max] = (k * Fg[lg_max - 1] * lg_max / (2.0 * lg_max + 1.0)
                           - (lg_max + 1.0) / tau_safe * Fg[lg_max]
                           - abs_kd * Fg[lg_max])

        # --- Neutrino hierarchy (same as standard sync) ---
        dFn = np.zeros(ln_max + 1)
        dFn[0] = -k * Fn[1] - (2.0 / 3.0) * h_prime
        F2n = Fn[2] if ln_max >= 2 else 0.0
        dFn[1] = (k / 3.0) * (Fn[0] - 2.0 * F2n)
        if ln_max >= 2:
            F3n = Fn[3] if ln_max >= 3 else 0.0
            dFn[2] = ((k / 5.0) * (2.0 * Fn[1] - 3.0 * F3n)
                      + (4.0 / 15.0) * h_prime + (8.0 / 15.0) * eta_prime)
        for ell in range(3, ln_max):
            dFn[ell] = (k / (2.0 * ell + 1.0)
                        * (ell * Fn[ell - 1] - (ell + 1) * Fn[ell + 1]))
        if ln_max >= 1:
            dFn[ln_max] = (k * Fn[ln_max - 1] * ln_max / (2.0 * ln_max + 1.0)
                           - (ln_max + 1.0) / tau_safe * Fn[ln_max])

        # --- Assemble ---
        dydt = np.zeros(nvar)
        dydt[IDX_ETA] = eta_prime
        dydt[IDX_DELTA_K] = d_delta_K
        dydt[IDX_DELTA_B] = d_delta_b
        dydt[IDX_THETA_B] = d_theta_b
        dydt[IDX_FG_START: IDX_FG_START + lg_max + 1] = dFg
        dydt[fn_s: fn_s + ln_max + 1] = dFn
        dydt[i_thetaK] = d_theta_K

        return dydt

    return rhs, nvar


# ============================================================================
# Khronon initial conditions (adiabatic, with theta_K)
# ============================================================================

def adiabatic_ic_khronon(k, tau_init, lg_max=L_GAMMA_MAX,
                         ln_max=L_NU_MAX_SYNC):
    """
    Adiabatic ICs for Khronon in sync gauge.

    At early times, c_s^2 -> 0 for all k (since a -> 0 makes k_J -> inf),
    so Khronon ICs are identical to CDM ICs with theta_K = 0.
    The extended state vector has theta_K appended.
    """
    y0_std = adiabatic_ic_sync(k, tau_init, lg_max, ln_max)
    # Append theta_K = 0 (adiabatic, CDM-like at early times)
    y0 = np.zeros(len(y0_std) + 1, dtype=np.float64)
    y0[:len(y0_std)] = y0_std
    y0[-1] = 0.0  # theta_K
    return y0


# ============================================================================
# Solve a single k-mode with Khronon DBI
# ============================================================================

def solve_single_k_khronon(k, bg, tau_end, lg_max=L_GAMMA_MAX,
                           ln_max=L_NU_MAX_SYNC,
                           alpha_dbi=ALPHA_DBI_DEFAULT,
                           method='Radau', rtol=1e-6, atol=1e-9):
    """
    Solve the Khronon sync gauge Boltzmann equations for a single k-mode.

    Returns the OdeSolution object. The state vector has theta_K appended
    at position n_var_sync(lg_max, ln_max).
    """
    tau_init = bg.tau_grid[1]
    rhs_fn, nvar = make_khronon_sync_rhs(k, bg, lg_max, ln_max, alpha_dbi)
    y0 = adiabatic_ic_khronon(k, tau_init, lg_max, ln_max)

    sol = solve_ivp(rhs_fn, [tau_init, tau_end], y0,
                    method=method, dense_output=True,
                    rtol=rtol, atol=atol)
    return sol


# ============================================================================
# Diagnose metric for Khronon state (handles theta_K in constraints)
# ============================================================================

def diagnose_metric_khronon(y, k, calH, a, lg_max, ln_max):
    """
    Compute h' and eta' from Einstein constraints with Khronon theta_K.

    The 0i-constraint now includes the Khronon velocity:
      eta' = (3/2) H0^2/k^2 [...  + Omega_K/a * theta_K]
    """
    k2 = k * k
    H02 = _H0_MPC ** 2
    Og, On, Ob, Ok = _OMEGA_GAMMA, _OMEGA_NU, _OMEGA_B, _OMEGA_C
    fn_s = idx_fn_start(lg_max)
    i_thetaK = n_var_sync(lg_max, ln_max)
    ia = 1.0 / a
    ia2 = ia * ia

    eta = y[IDX_ETA]
    delta_K = y[IDX_DELTA_K]
    delta_b = y[IDX_DELTA_B]
    theta_b = y[IDX_THETA_B]
    theta_K = y[i_thetaK]
    delta_g = y[idx_fg(0)]
    delta_n = y[fn_s]
    theta_g = 0.75 * k * y[idx_fg(1)]
    theta_n = 0.75 * k * y[fn_s + 1]

    # 00-constraint (same structure, delta_K in CDM slot)
    src00 = (Og * ia2 * delta_g + On * ia2 * delta_n
             + Ob * ia * delta_b + Ok * ia * delta_K)
    h_prime = (2.0 / calH) * (k2 * eta + 1.5 * H02 * src00)

    # 0i-constraint (now includes theta_K)
    src0i = ((4.0 / 3.0) * Og * ia2 * theta_g
             + (4.0 / 3.0) * On * ia2 * theta_n
             + Ob * ia * theta_b
             + Ok * ia * theta_K)
    eta_prime = 1.5 * H02 / k2 * src0i

    return h_prime, eta_prime


# ============================================================================
# Full Khronon sync gauge pipeline
# ============================================================================

def run_khronon_sync_solver(N_k=100, k_min=3e-4, k_max=0.35,
                            alpha_dbi=ALPHA_DBI_DEFAULT,
                            method='Radau', lg_max=L_GAMMA_MAX,
                            ln_max=L_NU_MAX_SYNC,
                            rtol=1e-6, atol=1e-9, verbose=True):
    """
    Full Khronon sync gauge pipeline: background -> perturbations -> C_l.

    Mirrors solver_sync.run_sync_solver but with Khronon DBI equations.
    """
    t_total = time.time()

    # ================================================================
    # Step 1: Background
    # ================================================================
    if verbose:
        print("=" * 65)
        print("Khronon Synchronous Gauge Boltzmann Solver")
        print(f"  alpha_DBI = {alpha_dbi:.4f}")
        print("=" * 65)
        print("\n--- Step 1: Background ---")

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

    tau_vis, tau_late, tau_all = build_snapshot_grid(
        bg, N_vis=60, N_early_isw=30, N_late_isw=20)
    N_snap = len(tau_all)

    k_J_rec = np.sqrt(1.5 * _OMEGA_C * _H0_MPC**2 / a_rec)

    if verbose:
        nvar = n_var_khronon(lg_max, ln_max)
        print(f"\n--- Step 2: Perturbations ({N_k} k-modes, {method}) ---")
        print(f"k range: [{k_min:.1e}, {k_max:.2f}] Mpc^-1")
        print(f"Snapshot grid: {N_snap} tau points")
        print(f"Hierarchy: lg_max={lg_max}, ln_max={ln_max}, nvar={nvar}")
        print(f"k_J(z_rec) = {k_J_rec:.4f} Mpc^-1, l_J ~ {k_J_rec * bg.D_A:.0f}")
        sys.stdout.flush()

    # ================================================================
    # Step 3: Solve all k-modes
    # ================================================================
    t0 = time.time()
    nvar_std = n_var_sync(lg_max, ln_max)
    fn_s = idx_fn_start(lg_max)

    # Storage: Newtonian gauge quantities at each (k, tau) snapshot
    Phi_N = np.zeros((N_k, N_snap))
    Psi_N = np.zeros((N_k, N_snap))
    Theta0_N = np.zeros((N_k, N_snap))
    vb_N = np.zeros((N_k, N_snap))
    Theta2_arr = np.zeros((N_k, N_snap))

    a_snap = bg.a_at_tau(tau_all)
    calH_snap = bg.calH_at_tau(tau_all)

    n_failed = 0
    for ik, k_val in enumerate(k_arr):
        sol = solve_single_k_khronon(
            k_val, bg, tau_all[-1] + 1.0,
            lg_max, ln_max, alpha_dbi, method, rtol, atol)

        if not sol.success:
            n_failed += 1
            if verbose and n_failed <= 3:
                print(f"  WARNING: k={k_val:.4e} failed: {sol.message}")
            continue

        for it in range(N_snap):
            tau = tau_all[it]
            y = sol.sol(tau)
            calH = calH_snap[it]
            a = a_snap[it]

            # Use Khronon-aware metric diagnosis
            h_prime, eta_prime = diagnose_metric_khronon(
                y, k_val, calH, a, lg_max, ln_max)

            # Gauge transform uses the standard sync variables (first nvar_std)
            # The gauge_transform function only reads eta, Fg, Fn, theta_b
            # which are in the same positions.
            gt = gauge_transform(y[:nvar_std], k_val, calH, a,
                                 lg_max, ln_max, h_prime, eta_prime)

            Phi_N[ik, it] = gt['Phi_N']
            Psi_N[ik, it] = gt['Psi_N']

            eta_val = y[IDX_ETA]
            Theta0_N[ik, it] = y[IDX_FG_START] / 4.0 - eta_val + gt['Phi_N']
            vb_N[ik, it] = (y[IDX_THETA_B] / k_val
                            + (h_prime + 6 * eta_prime) / (2 * k_val))
            Theta2_arr[ik, it] = y[IDX_FG_START + 2]

        if verbose and (ik + 1) % 20 == 0:
            elapsed = time.time() - t0
            eta_est = elapsed / (ik + 1) * N_k
            print(f"  [{ik+1}/{N_k}] elapsed={elapsed:.1f}s, ETA={eta_est:.1f}s")
            sys.stdout.flush()

    t_pert = time.time() - t0
    if verbose:
        print(f"Perturbations: {t_pert:.1f}s ({n_failed} failed)")

    # ================================================================
    # Step 4: Phi_N' + Psi_N' by finite differences
    # ================================================================
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

    P_R = A_s * (k_arr / k_pivot) ** (n_s - 1.0)

    g_snap = bg.visibility_at_tau(tau_all)
    kappa_snap = bg.kappa_at_tau(tau_all)
    exp_neg_kappa = np.exp(-kappa_snap)
    chi_snap = bg.tau_0 - tau_all

    # Bessel table
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

    if verbose:
        print(f"Computing LOS transfer functions...")
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
            w = 0.5 * (dtau_all[max(0, it - 1)]
                       + dtau_all[min(it, len(dtau_all) - 1)])

        S_SW = g_snap[it] * (Theta0_N[:, it] + Psi_N[:, it])
        S_Dop = g_snap[it] * vb_N[:, it]
        S_ISW = exp_neg_kappa[it] * PhiPsi_prime[:, it]

        if use_gpu_bessel:
            x_arr = k_arr * chi
            jl_mx, jlp_mx = bessel_table.eval_jl_jlp(x_arr.astype(np.float32))
            mx.eval(jl_mx, jlp_mx)
            jl = np.array(jl_mx)
            jlp = np.array(jlp_mx)

            for il in range(N_ell):
                Delta_l[il, :] += w * (
                    S_SW * jl[il, :] +
                    S_Dop * jlp[il, :] +
                    S_ISW * jl[il, :]
                )
        else:
            from scipy.special import spherical_jn
            x = k_arr * chi
            for il, ell in enumerate(ell_values):
                jl = spherical_jn(int(ell), x)
                jlp = spherical_jn(int(ell), x, derivative=True)
                Delta_l[il, :] += w * (S_SW * jl + S_Dop * jlp + S_ISW * jl)

    t_los = time.time() - t0
    if verbose:
        print(f"LOS integration: {t_los:.2f}s")

    # ================================================================
    # Step 6: C_l
    # ================================================================
    t0 = time.time()

    lnk = np.log(k_arr)
    dlnk = np.diff(lnk)

    integrand = P_R[None, :] * Delta_l ** 2
    mid = 0.5 * (integrand[:, :-1] + integrand[:, 1:])
    Cl = 4.0 * np.pi * np.sum(mid * dlnk[None, :], axis=1)
    Cl = np.maximum(Cl, 0.0)

    ell_f = ell_values.astype(float)
    norm = (2.0 / 3.0) ** 2
    Dl = norm * ell_f * (ell_f + 1.0) * Cl / (2.0 * np.pi) * (T_CMB * 1e6) ** 2

    t_cl = time.time() - t0
    t_elapsed = time.time() - t_total

    if verbose:
        print(f"C_l computation: {t_cl:.3f}s")
        print(f"\nTotal time: {t_elapsed:.1f}s")

    return {
        'ell': ell_values,
        'Cl': Cl,
        'Dl': Dl,
        'k_arr': k_arr,
        'bg': bg,
        'Delta_l': Delta_l,
        'Phi_N': Phi_N,
        'Psi_N': Psi_N,
        'alpha_dbi': alpha_dbi,
        'timing': {
            'background': t_bg,
            'perturbations': t_pert,
            'los': t_los,
            'total': t_elapsed,
        },
    }


# ============================================================================
# Comparison: LCDM (sync) vs Khronon (sync)
# ============================================================================

def sync_khronon_comparison(N_k=100, alpha_dbi=ALPHA_DBI_DEFAULT,
                            save_dir=None, verbose=True):
    """
    Run LCDM and Khronon solvers both in sync gauge and compare.

    This provides an apples-to-apples comparison where the ONLY difference
    is the Khronon DBI c_s^2 in the dark matter perturbation equations.

    Parameters
    ----------
    N_k : int
        Number of k-modes (default: 100 for manageable runtime).
    alpha_dbi : float
        DBI strength parameter (default: 0.5, ghost condensation).
    save_dir : str or None
        Directory for output files.
    verbose : bool
        Print progress.

    Returns
    -------
    dict with keys: ell, Dl_lcdm, Dl_khronon, residual, timing
    """
    if save_dir is None:
        save_dir = os.path.dirname(os.path.abspath(__file__))

    if verbose:
        print("=" * 65)
        print("KHRONON vs LCDM: Synchronous Gauge Comparison")
        print(f"  N_k = {N_k}, alpha_DBI = {alpha_dbi:.4f}")
        print("=" * 65)

    t_start = time.time()

    # --- LCDM (standard sync gauge) ---
    if verbose:
        print("\n>>> Running LCDM (sync gauge)...")

    from .solver_sync import run_sync_solver
    result_lcdm = run_sync_solver(
        N_k=N_k, verbose=verbose)
    ell = result_lcdm['ell']
    Dl_lcdm = result_lcdm['Dl']

    # --- Khronon DBI (sync gauge) ---
    if verbose:
        print("\n>>> Running Khronon DBI (sync gauge)...")

    result_khronon = run_khronon_sync_solver(
        N_k=N_k, alpha_dbi=alpha_dbi, verbose=verbose)
    Dl_khronon = result_khronon['Dl']

    # --- Residual ---
    Dl_ref = np.maximum(Dl_lcdm, 1e-10)
    residual = (Dl_khronon - Dl_lcdm) / Dl_ref

    t_total = time.time() - t_start

    # --- Summary ---
    if verbose:
        print("\n" + "=" * 65)
        print("SYNC GAUGE COMPARISON SUMMARY")
        print("=" * 65)
        print(f"Total time: {t_total:.1f}s")
        print(f"  LCDM:    {result_lcdm['timing']['total']:.1f}s")
        print(f"  Khronon: {result_khronon['timing']['total']:.1f}s")

        for region, mask in [
            ("l < 500", ell < 500),
            ("500-2000", (ell >= 500) & (ell < 2000)),
            ("l >= 2000", ell >= 2000),
        ]:
            if np.any(mask):
                r = residual[mask]
                print(f"  {region:10s}: max |res| = {np.max(np.abs(r)):.4e}, "
                      f"mean = {np.mean(np.abs(r)):.4e}")

        k_J_rec = np.sqrt(1.5 * _OMEGA_C * _H0_MPC**2 / a_rec)
        print(f"\n  k_J(z_rec) = {k_J_rec:.4f} Mpc^-1")

    # --- Plot ---
    _plot_sync_comparison(ell, Dl_lcdm, Dl_khronon, residual,
                          alpha_dbi, t_total, save_dir)

    # --- Save data ---
    data_path = os.path.join(save_dir, 'khronon_sync_comparison.dat')
    np.savetxt(data_path,
               np.column_stack([ell, Dl_lcdm, Dl_khronon, residual]),
               header='ell  Dl_LCDM  Dl_Khronon  residual',
               fmt='%6d' + '  %.8e' * 3)
    if verbose:
        print(f"\nData: {data_path}")

    return {
        'ell': ell,
        'Dl_lcdm': Dl_lcdm,
        'Dl_khronon': Dl_khronon,
        'residual': residual,
        'result_lcdm': result_lcdm,
        'result_khronon': result_khronon,
        'timing': {
            'lcdm': result_lcdm['timing']['total'],
            'khronon': result_khronon['timing']['total'],
            'total': t_total,
        },
    }


def _plot_sync_comparison(ell, Dl_lcdm, Dl_khronon, residual,
                          alpha_dbi, t_total, save_dir):
    """Generate the 2-panel sync gauge comparison figure."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 9),
                                        gridspec_kw={'height_ratios': [3, 1.5]})

        # Panel 1: D_l overlay
        ax1.plot(ell, Dl_lcdm, 'b-', lw=1.5, alpha=0.8,
                 label=r'$\Lambda$CDM (sync gauge)')
        ax1.plot(ell, Dl_khronon, 'r--', lw=1.5, alpha=0.7,
                 label=r'Khronon DBI ($\alpha$='
                       + f'{alpha_dbi:.2f}, sync gauge)')
        ax1.set_ylabel(r'$D_\ell$ [$\mu K^2$]', fontsize=13)
        ax1.set_title(
            f'Sync Gauge: Khronon vs LCDM  '
            r'($\alpha_{\rm DBI}$' + f' = {alpha_dbi:.2f}, '
            f'{t_total:.0f}s)',
            fontsize=14)
        ax1.legend(fontsize=11, loc='upper right')
        ax1.set_xlim(2, ell[-1])
        if np.max(Dl_lcdm) > 0:
            ax1.set_ylim(0, np.max(Dl_lcdm) * 1.3)
        ax1.tick_params(labelsize=11)

        # Panel 2: Residual
        ax2.plot(ell, residual * 100, 'r-', lw=1.0, alpha=0.8)
        ax2.axhline(0, color='gray', ls='-', alpha=0.3)
        ax2.axhline(1, color='red', ls=':', alpha=0.3, label='1% level')
        ax2.axhline(-1, color='red', ls=':', alpha=0.3)
        ax2.set_ylabel(r'Residual [%]', fontsize=12)
        ax2.set_xlabel(r'$\ell$', fontsize=13)
        ax2.legend(fontsize=9, loc='best')
        ax2.set_xlim(2, ell[-1])

        valid = np.isfinite(residual)
        if np.any(valid):
            ymax = min(max(np.max(np.abs(residual[valid])) * 130, 2.0), 200.0)
            ax2.set_ylim(-ymax, ymax)
        ax2.tick_params(labelsize=11)

        plt.tight_layout()
        out_path = os.path.join(save_dir, 'khronon_sync_vs_lcdm.png')
        plt.savefig(out_path, dpi=150)
        if True:  # verbose
            print(f"Plot: {out_path}")
        plt.close()

    except ImportError:
        print("[khronon_sync] matplotlib not available, skipping plot")


# ============================================================================
# CLI
# ============================================================================

def main():
    import argparse
    parser = argparse.ArgumentParser(
        description='Khronon DBI vs LCDM in synchronous gauge')
    parser.add_argument('--N_k', type=int, default=100,
                        help='Number of k-modes (default: 100)')
    parser.add_argument('--alpha_dbi', type=float, default=ALPHA_DBI_DEFAULT,
                        help=f'DBI strength (default: {ALPHA_DBI_DEFAULT})')
    parser.add_argument('--khronon_only', action='store_true',
                        help='Run only the Khronon solver (skip comparison)')
    args = parser.parse_args()

    if args.khronon_only:
        result = run_khronon_sync_solver(
            N_k=args.N_k, alpha_dbi=args.alpha_dbi, verbose=True)
        ell = result['ell']
        Dl = result['Dl']
        print("\n--- D_l values at key multipoles ---")
        for target_ell in [2, 10, 50, 100, 220, 500, 800, 1000, 1500, 2000]:
            idx = np.argmin(np.abs(ell - target_ell))
            print(f"  l={ell[idx]:5d}: D_l = {Dl[idx]:10.2f} uK^2")
    else:
        sync_khronon_comparison(
            N_k=args.N_k, alpha_dbi=args.alpha_dbi, verbose=True)


if __name__ == '__main__':
    main()
