"""
solver_sync_class_los.py -- CLASS exact "efficient" source function LOS integrator.

Implements CLASS's NEWTONIAN GAUGE efficient source functions (perturbations.c
lines 7518-7528), computed from synchronous gauge ODE solutions via gauge
transformation to Bardeen potentials.

CLASS EFFICIENT SOURCES (Newtonian gauge, perturbations.c lines 7518-7528):

  S_t0 = [SW]   g * (delta_g/4 + Psi)
       + [ISW]  g * (Phi - Psi)  +  exp(-kappa) * 2 * Phi'
       + [Dop]  (1/k^2) * [g * theta_b_N_dot  +  g' * theta_b_N]

  S_t1 = exp(-kappa) * k * (Psi - Phi)
       [TINY: proportional to anisotropic stress only]

  S_t2 = g * P

The full temperature transfer function:
  Delta_T_l(k) = int dtau [S_t0 * j_l  +  S_t1 * j_l'  +  S_t2 * R_l^(2)]

where R_l^(2) = (3*j_l'' + j_l)/2 in flat space.

KEY INSIGHT (CLASS, perturbations.c line 7518):
  The "simple" Newtonian form puts g*theta_b_N/k in S_t1 (via j_l'),
  alongside exp(-kappa)*k*psi. These two terms nearly cancel at recombination,
  causing numerical precision loss. The "efficient" form uses integration by
  parts to move the Doppler contribution into S_t0 (via j_l), leaving only
  the tiny anisotropic stress piece in S_t1.

KEY FEATURES:
  - ODE solved in synchronous gauge (stable, no gauge mode)
  - Gauge transform to Bardeen potentials Phi, Psi for source assembly
  - Phi' computed analytically from the Poisson equation
  - theta_b_N and theta_b_N_dot from gauge transformation + Boltzmann EOM
  - Dense tau grid (~500+ points in visibility) to resolve g'(tau)
  - Full polarization support via l_pol_max parameter
  - S_t2 quadrupole/polarization via specialized radial function

Usage:
  python -m mlx_class.solver_sync_class_los [--N_k 500] [--parallel]

Author: Sheng-Kai Huang, 2026
"""
import os
os.environ['MPLCONFIGDIR'] = '/tmp/mpl_config'

import numpy as np
import time
import sys
from multiprocessing import Pool

from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d, CubicSpline

from .background import (
    Background, A_s, n_s, k_pivot, T_CMB,
    Omega_r as _OMEGA_R,
    Omega_b as _OMEGA_B,
    Omega_c as _OMEGA_C,
    Omega_m as _OMEGA_M,
    Omega_L as _OMEGA_L,
    H0_Mpc as _H0_MPC,
)
from .perturbations_neutrino import (
    _f_nu, _OMEGA_GAMMA, _OMEGA_NU,
)
from .perturbations_sync import (
    make_sync_rhs, adiabatic_ic_sync, diagnose_metric, gauge_transform,
    n_var_sync, idx_fn_start, idx_fg, idx_e_start, idx_e,
    IDX_ETA, IDX_DELTA_C, IDX_DELTA_B, IDX_THETA_B, IDX_FG_START,
    L_GAMMA_MAX, L_NU_MAX_SYNC, L_POL_MAX,
)
from .solver_sync import build_snapshot_grid


# ============================================================================
# Build a DENSE tau grid suitable for efficient source functions
# ============================================================================

def build_dense_snapshot_grid(bg, N_vis=500, N_early_isw=60, N_late_isw=40,
                              N_reion=30):
    """
    Build a dense tau grid for LOS integration with efficient source functions.

    The efficient source functions involve g'(tau) (the visibility derivative),
    which is a sharp, sign-changing function around the visibility peak.
    CLASS uses O(1000-2000) adaptively-spaced tau points. We use ~500 in the
    visibility region plus denser coverage elsewhere.

    Returns
    -------
    tau_vis, tau_late, tau_all
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

    tau_vis_lo = max(tau_peak - 6.0 * sigma, bg.tau_grid[5])
    tau_vis_hi = min(tau_peak + 6.0 * sigma, bg.tau_0 * 0.5)
    tau_vis = np.linspace(tau_vis_lo, tau_vis_hi, N_vis)

    a_eq = _OMEGA_R / _OMEGA_M
    tau_eq_approx = float(bg._tau_of_a(np.log(a_eq)))
    tau_early_lo = max(tau_eq_approx * 0.3, bg.tau_grid[5])
    tau_early_hi = tau_vis_lo

    if tau_early_hi > tau_early_lo * 1.5:
        tau_isw_dense_lo = max(tau_eq_approx - 60.0, tau_early_lo)
        tau_isw_dense_hi = min(tau_eq_approx + 60.0, tau_early_hi)
        tau_dense = np.linspace(tau_isw_dense_lo, tau_isw_dense_hi, N_early_isw)
        tau_sparse = np.geomspace(tau_early_lo, tau_early_hi, 20)
        tau_early = np.sort(np.unique(np.concatenate([tau_sparse, tau_dense])))
    else:
        tau_early = np.array([tau_early_lo])

    tau_reion = np.array([])
    if hasattr(bg, 'z_reio') and bg.tau_reio > 0:
        z_reio = bg.z_reio
        a_reio = 1.0 / (1.0 + z_reio)
        tau_reio_center = float(bg._tau_of_a(np.log(a_reio)))
        tau_reion_lo = max(tau_reio_center - 200.0, tau_vis_hi + 10.0)
        tau_reion_hi = min(tau_reio_center + 200.0, bg.tau_0 * 0.95)
        if tau_reion_hi > tau_reion_lo:
            tau_reion = np.linspace(tau_reion_lo, tau_reion_hi, N_reion)

    tau_late_lo = tau_vis_hi
    tau_late_hi = bg.tau_0 * 0.98
    tau_late = np.linspace(tau_late_lo, tau_late_hi, N_late_isw)

    tau_all = np.sort(np.unique(np.concatenate([
        tau_early, tau_vis, tau_reion, tau_late])))
    return tau_vis, tau_late, tau_all


# ============================================================================
# Compute visibility derivative g'(tau) via cubic spline
# ============================================================================

def compute_g_prime(bg, tau_arr):
    """Compute g'(tau) = dg/dtau via cubic spline on the background grid."""
    cs = CubicSpline(bg.tau_grid, bg.visibility_grid)
    return cs(tau_arr, 1)


# ============================================================================
# Worker: solve one k-mode, return Newtonian gauge quantities
# ============================================================================

def _solve_single_k_class_worker(args):
    """
    Worker for multiprocessing. Solves one k-mode and returns quantities
    needed for the Newtonian gauge efficient source functions.

    Returns numpy array of shape (N_snap, N_col) with columns:
      0: Phi_N       -- Bardeen potential Phi
      1: Psi_N       -- Bardeen potential Psi
      2: Theta0_N    -- Newtonian gauge temperature monopole (delta_g_N/4)
      3: theta_b_N   -- Newtonian gauge baryon velocity divergence
      4: theta_b_N_dot -- d(theta_b_N)/dtau (analytical from Boltzmann + gauge)
      5: Phi_N_prime -- d(Phi_N)/dtau from Poisson + sync variables
      6: Theta2      -- F_gamma,2/4 (gauge-invariant)
      7: P_pol       -- polarization combination P = (E_0+E_2+2*F_g2)/8

    Returns None if the ODE solve fails.
    """
    import numpy as np
    from scipy.integrate import solve_ivp

    (k, bg_arrays, tau_end, lg_max, ln_max, l_pol_max, method, rtol, atol,
     tau_all, a_snap, calH_snap) = args

    class _BGLite:
        __slots__ = ('tau_grid', 'a_grid', 'calH_grid', 'R_grid',
                     'kappa_dot_grid', 'tau_rec', 'tau_0')
        def __init__(self, arrays):
            self.tau_grid = arrays['tau_grid']
            self.a_grid = arrays['a_grid']
            self.calH_grid = arrays['calH_grid']
            self.R_grid = arrays['R_grid']
            self.kappa_dot_grid = arrays['kappa_dot_grid']
            self.tau_rec = arrays['tau_rec']
            self.tau_0 = arrays['tau_0']

    bg_lite = _BGLite(bg_arrays)

    from mlx_class.perturbations_sync import (
        make_sync_rhs, adiabatic_ic_sync, diagnose_metric, gauge_transform,
        IDX_ETA, IDX_DELTA_C, IDX_DELTA_B, IDX_THETA_B, IDX_FG_START,
        idx_fn_start as _idx_fn_start, idx_e_start as _idx_e_start,
    )
    from mlx_class.background import H0_Mpc as _H0_val
    from mlx_class.perturbations_neutrino import _OMEGA_GAMMA as Og, _OMEGA_NU as On

    H02 = _H0_val ** 2
    k2 = k * k

    _has_pol = l_pol_max > 0
    fn_s = _idx_fn_start(lg_max, l_pol_max)
    if _has_pol:
        _e_s = _idx_e_start(lg_max)

    tau_init = bg_lite.tau_grid[1]
    rhs_fn, nvar = make_sync_rhs(k, bg_lite, lg_max, ln_max, l_pol_max)
    y0 = adiabatic_ic_sync(k, tau_init, lg_max, ln_max, l_pol_max)

    sol = solve_ivp(rhs_fn, [tau_init, tau_end], y0,
                    method=method, dense_output=True,
                    rtol=rtol, atol=atol)

    if not sol.success:
        return None

    _tau = bg_lite.tau_grid
    _calH = bg_lite.calH_grid
    _a = bg_lite.a_grid
    _R = bg_lite.R_grid
    _kd = bg_lite.kappa_dot_grid

    N_snap = len(tau_all)
    N_COL = 8
    results = np.zeros((N_snap, N_COL), dtype=np.float64)

    for it in range(N_snap):
        tau = tau_all[it]
        y = sol.sol(tau)
        a = a_snap[it]
        calH = calH_snap[it]

        R = float(np.interp(tau, _tau, _R))
        kappa_dot = float(np.interp(tau, _tau, _kd))
        abs_kd = abs(kappa_dot)

        eta_val = y[IDX_ETA]
        theta_b_S = y[IDX_THETA_B]
        Fg = y[IDX_FG_START: IDX_FG_START + lg_max + 1]
        Fn = y[fn_s: fn_s + ln_max + 1]

        theta_g_S = 0.75 * k * Fg[1]
        sigma_g = 0.5 * Fg[2] if lg_max >= 2 else 0.0
        sigma_n = 0.5 * Fn[2] if ln_max >= 2 else 0.0

        # Diagnose h' and eta' from Einstein constraints
        h_prime, eta_prime = diagnose_metric(y, k, calH, a, lg_max, ln_max,
                                             l_pol_max)

        # Gauge transformation to Newtonian potentials
        gt = gauge_transform(y, k, calH, a, lg_max, ln_max, h_prime,
                             eta_prime, l_pol_max)
        Phi_N = gt['Phi_N']
        Psi_N = gt['Psi_N']
        alpha_val = gt['alpha']

        # Theta0_N: gauge-invariant temperature monopole
        # = F_g0/4 - eta + Phi_N (avoids catastrophic cancellation)
        Theta0_N = Fg[0] / 4.0 - eta_val + Phi_N

        # theta_b_N = theta_b_S + k^2 * alpha
        theta_b_N = theta_b_S + k2 * alpha_val

        # alpha' from 4th Einstein equation:
        # alpha' = -2*calH*alpha + eta - (9/2)*(a^2/k^2)*(rho+p)*sigma
        aniso_stress = 6.0 * H02 / (k2 * a * a) * (Og * sigma_g + On * sigma_n)
        alpha_prime = -2.0 * calH * alpha_val + eta_val - aniso_stress

        # theta_b_S_dot from Boltzmann EOM:
        theta_b_S_dot = -calH * theta_b_S + abs_kd / R * (theta_g_S - theta_b_S)

        # theta_b_N_dot = theta_b_S_dot + k^2 * alpha'
        theta_b_N_dot = theta_b_S_dot + k2 * alpha_prime

        # Phi_N_prime from sync gauge:
        # Phi_N = eta - calH*alpha
        # Phi_N' = eta' - calH'*alpha - calH*alpha'
        # where calH' = d(calH)/dtau
        # But we can also use the Newtonian Poisson equation:
        # Phi' = -calH*Psi + (3/2)*(calH^2/k^2)*sum (1+w)*Omega_i*theta_i/a^{1+3w}
        # This is what CLASS uses (perturbations.c line 6494):
        # Phi' = -calH*Psi + (3/2)*(a^2/k^2)*(rho+p)*theta_tot
        # where (rho+p)*theta_tot = H0^2*[(4/3)*Og/a^4*theta_g + (4/3)*On/a^4*theta_n + Ob/a^3*theta_b]_N
        # In Newtonian gauge.
        #
        # Actually, simpler: from sync variables:
        # eta' is already computed from the 0i constraint (momentum)
        # alpha is already computed
        # We need calH' -- but to avoid computing it per-snapshot, use the
        # identity: Phi_N' = eta' - (a''/a)*alpha - calH*alpha'
        # Wait no, Phi = eta - calH*alpha, so
        # Phi' = eta' - calH'*alpha - calH*alpha'
        # This requires calH'. But calH' is:
        # calH' = calH^2 - 0.5*H0^2*(4*Omega_r/a^2 + 3*Omega_m/a)
        from mlx_class.background import Omega_r as Or, Omega_m as Om
        calH_prime = calH**2 - 0.5 * H02 * (4.0 * Or / a**2 + 3.0 * Om / a)
        Phi_N_prime = eta_prime - calH_prime * alpha_val - calH * alpha_prime

        # Polarization combination P
        if _has_pol:
            E0 = y[_e_s]
            E2 = y[_e_s + 2] if l_pol_max >= 2 else 0.0
            P_pol = (E0 + E2 + 2.0 * Fg[2]) / 8.0
        else:
            P_pol = Fg[2] / 4.0 if lg_max >= 2 else 0.0

        results[it, 0] = Phi_N
        results[it, 1] = Psi_N
        results[it, 2] = Theta0_N
        results[it, 3] = theta_b_N
        results[it, 4] = theta_b_N_dot
        results[it, 5] = Phi_N_prime
        results[it, 6] = Fg[2] / 4.0 if lg_max >= 2 else 0.0  # Theta_2
        results[it, 7] = P_pol

    return results


# ============================================================================
# Full pipeline with CLASS efficient source functions
# ============================================================================

def run_sync_class_los(N_k=500, k_min=3e-4, k_max=0.35, method='BDF',
                       lg_max=L_GAMMA_MAX, ln_max=L_NU_MAX_SYNC,
                       l_pol_max=L_POL_MAX,
                       rtol=1e-6, atol=1e-9, n_workers=None,
                       N_vis=500, N_early_isw=60, N_late_isw=40,
                       verbose=True):
    """
    Parallel synchronous gauge pipeline with CLASS exact efficient source functions.

    The ODE is solved in synchronous gauge (stable), then gauge-transformed to
    Bardeen potentials. The LOS integration uses CLASS's Newtonian-gauge efficient
    source formulas which eliminate the Doppler-gravitational cancellation.

    Parameters
    ----------
    N_k : int
        Number of k-modes (default 500).
    k_min, k_max : float
        k range in Mpc^-1.
    method : str
        ODE solver method ('BDF' or 'Radau').
    lg_max, ln_max : int
        Hierarchy truncation.
    l_pol_max : int
        E-mode polarization hierarchy truncation. 8 = recommended.
    rtol, atol : float
        ODE solver tolerances.
    n_workers : int or None
        Number of worker processes.
    N_vis : int
        Tau points in visibility region (500 for efficient).
    verbose : bool
        Print progress.

    Returns
    -------
    dict with ell, Dl_TT, Dl_EE, Dl_TE, etc.
    """
    t_total = time.time()

    if n_workers is None:
        n_workers = min(os.cpu_count() or 4, 8)

    # ================================================================
    # Step 1: Background
    # ================================================================
    if verbose:
        print("=" * 70)
        print("Synchronous Gauge Solver -- CLASS Efficient Source Functions")
        print("=" * 70)
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

    tau_vis, tau_late, tau_all = build_dense_snapshot_grid(
        bg, N_vis=N_vis, N_early_isw=N_early_isw, N_late_isw=N_late_isw)
    N_snap = len(tau_all)

    if verbose:
        print(f"\n--- Step 2: Perturbations ({N_k} k-modes, {method}, "
              f"{n_workers} workers) ---")
        print(f"k range: [{k_min:.1e}, {k_max:.2f}] Mpc^-1")
        print(f"Snapshot grid: {N_snap} tau points (DENSE for efficient sources)")
        print(f"  Visibility: {len(tau_vis)} pts, "
              f"tau=[{tau_vis[0]:.1f}, {tau_vis[-1]:.1f}]")
        print(f"  Late ISW:   {len(tau_late)} pts, "
              f"tau=[{tau_late[0]:.1f}, {tau_late[-1]:.1f}]")
        if hasattr(bg, 'z_reio') and bg.tau_reio > 0:
            print(f"  Reionization: z_reio={bg.z_reio:.2f}, "
                  f"tau_reio={bg.tau_reio:.4f}")
        print(f"Hierarchy: lg_max={lg_max}, ln_max={ln_max}, "
              f"l_pol_max={l_pol_max}, "
              f"nvar={n_var_sync(lg_max, ln_max, l_pol_max)}")
        if l_pol_max > 0:
            print(f"  Polarization: ENABLED ({l_pol_max+1} E-mode multipoles)")
        sys.stdout.flush()

    # ================================================================
    # Step 3: Precompute background at snapshots
    # ================================================================
    a_snap = np.asarray(bg.a_at_tau(tau_all), dtype=np.float64)
    calH_snap = np.asarray(bg.calH_at_tau(tau_all), dtype=np.float64)
    g_snap = np.asarray(bg.visibility_at_tau(tau_all), dtype=np.float64)
    g_prime_snap = compute_g_prime(bg, tau_all)
    kappa_snap = np.asarray(bg.kappa_at_tau(tau_all), dtype=np.float64)
    exp_neg_kappa = np.exp(-kappa_snap)
    chi_snap = bg.tau_0 - tau_all

    # ================================================================
    # Step 4: Solve all k-modes in parallel
    # ================================================================
    t0 = time.time()

    bg_arrays = {
        'tau_grid': bg.tau_grid.copy(),
        'a_grid': bg.a_grid.copy(),
        'calH_grid': bg.calH_grid.copy(),
        'R_grid': bg.R_grid.copy(),
        'kappa_dot_grid': bg.kappa_dot_grid.copy(),
        'tau_rec': float(bg.tau_rec),
        'tau_0': float(bg.tau_0),
    }

    tau_end = float(tau_all[-1] + 1.0)

    work_items = [
        (k, bg_arrays, tau_end, lg_max, ln_max, l_pol_max, method, rtol, atol,
         tau_all, a_snap, calH_snap)
        for k in k_arr
    ]

    if verbose:
        print(f"Dispatching {N_k} k-modes to {n_workers} workers...")
        sys.stdout.flush()

    with Pool(processes=n_workers) as pool:
        worker_results = pool.map(_solve_single_k_class_worker, work_items)

    t_pert = time.time() - t0

    # ================================================================
    # Collect results
    # ================================================================
    Phi_N = np.zeros((N_k, N_snap))
    Psi_N = np.zeros((N_k, N_snap))
    Theta0_N = np.zeros((N_k, N_snap))
    theta_b_N = np.zeros((N_k, N_snap))
    theta_b_N_dot = np.zeros((N_k, N_snap))
    Phi_N_prime = np.zeros((N_k, N_snap))
    Theta2_arr = np.zeros((N_k, N_snap))
    P_pol_arr = np.zeros((N_k, N_snap))

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
        theta_b_N[ik] = res[:, 3]
        theta_b_N_dot[ik] = res[:, 4]
        Phi_N_prime[ik] = res[:, 5]
        Theta2_arr[ik] = res[:, 6]
        P_pol_arr[ik] = res[:, 7]

    if verbose:
        print(f"Perturbations: {t_pert:.1f}s ({n_failed} failed, "
              f"{N_k - n_failed} OK)")

    # ================================================================
    # Step 5: Smooth Phi_N_prime via cubic spline (optional refinement)
    # ================================================================
    # The analytical Phi_N_prime can be noisy for superhorizon modes.
    # Use cubic spline as a cross-check / smoothing.
    Phi_N_prime_spline = np.zeros_like(Phi_N)
    for ik in range(N_k):
        cs = CubicSpline(tau_all, Phi_N[ik])
        Phi_N_prime_spline[ik] = cs(tau_all, 1)

    # Use spline-smoothed derivative for the ISW term for stability
    # The analytical derivative is used as a cross-check
    # For modes where analytical and spline agree well, prefer analytical
    # For noisy superhorizon modes, spline is more stable
    for ik in range(N_k):
        # Relative difference between analytical and spline
        mask = np.abs(Phi_N_prime[ik]) > 1e-15
        if np.sum(mask) > 10:
            rel_diff = np.abs(Phi_N_prime[ik, mask] - Phi_N_prime_spline[ik, mask])
            max_diff = np.max(rel_diff)
            if max_diff > 0.1 * np.max(np.abs(Phi_N_prime[ik, mask])):
                # Spline is more stable for this mode
                Phi_N_prime[ik] = Phi_N_prime_spline[ik]

    # Also smooth theta_b_N_dot via spline as a cross-check
    theta_b_N_dot_spline = np.zeros_like(theta_b_N)
    for ik in range(N_k):
        cs = CubicSpline(tau_all, theta_b_N[ik])
        theta_b_N_dot_spline[ik] = cs(tau_all, 1)

    # Use analytical theta_b_N_dot (from Boltzmann EOM) by default,
    # fall back to spline if they differ too much
    for ik in range(N_k):
        mask = np.abs(theta_b_N_dot[ik]) > 1e-15
        if np.sum(mask) > 10:
            rel_diff = np.abs(theta_b_N_dot[ik, mask] - theta_b_N_dot_spline[ik, mask])
            max_diff = np.max(rel_diff)
            if max_diff > 0.3 * np.max(np.abs(theta_b_N_dot[ik, mask])):
                theta_b_N_dot[ik] = theta_b_N_dot_spline[ik]

    # ================================================================
    # Step 6: Assemble CLASS Newtonian efficient source functions
    # ================================================================
    if verbose:
        print(f"\n--- Step 6: Assembling CLASS efficient sources ---")
        sys.stdout.flush()

    # Pre-broadcast for vectorized computation
    g_2d = g_snap[np.newaxis, :]               # (1, N_snap)
    gp_2d = g_prime_snap[np.newaxis, :]        # (1, N_snap)
    emk_2d = exp_neg_kappa[np.newaxis, :]      # (1, N_snap)
    k2_col = (k_arr ** 2)[:, np.newaxis]       # (N_k, 1)
    k_col = k_arr[:, np.newaxis]               # (N_k, 1)

    # -----------------------------------------------------------------
    # S_t0: CLASS Newtonian gauge efficient monopole source
    # (perturbations.c lines 7518-7528, all switches = 1)
    #
    # S_t0 = g * (delta_g_N/4 + Psi)                    [SW]
    #       + g * (Phi - Psi) + exp(-kappa) * 2 * Phi'  [ISW]
    #       + (1/k^2) * (g * theta_b_N_dot + g' * theta_b_N)  [Doppler IBP]
    #
    # Combining g-terms:
    #   g * (delta_g_N/4 + Psi + Phi - Psi) = g * (Theta0_N + Phi)
    # But note: Theta0_N = delta_g_N/4 = F_g0/4 - eta + Phi_N (our storage)
    # So Theta0_N + Psi is the standard SW source.
    # Adding g*(Phi-Psi) gives: g*(Theta0_N + Phi).
    # -----------------------------------------------------------------

    # SW + ISW (g-part): g * (Theta0_N + Phi_N)
    # ISW (exp-part): exp(-kappa) * 2 * Phi'
    # Doppler IBP: (g * theta_b_N_dot + g' * theta_b_N) / k^2
    S_t0 = (
        g_2d * (Theta0_N + Phi_N)                    # SW + ISW g-part
        + emk_2d * 2.0 * Phi_N_prime                  # ISW exp-part
        + (g_2d * theta_b_N_dot + gp_2d * theta_b_N) / k2_col  # Doppler IBP
    )  # (N_k, N_snap)

    # S_t1: tiny dipole source (anisotropic stress only)
    # S_t1 = exp(-kappa) * k * (Psi - Phi)
    S_t1 = emk_2d * k_col * (Psi_N - Phi_N)  # (N_k, N_snap)

    # S_t2: polarization quadrupole source
    # S_t2 = g * P
    S_t2 = g_2d * P_pol_arr  # (N_k, N_snap)

    # S_p: E-mode polarization source
    # S_p = sqrt(6) * g * P
    S_p = np.sqrt(6.0) * g_2d * P_pol_arr  # (N_k, N_snap)

    if verbose:
        max_S_t0 = np.max(np.abs(S_t0))
        max_S_t1 = np.max(np.abs(S_t1))
        max_S_t2 = np.max(np.abs(S_t2))
        print(f"  max|S_t0| = {max_S_t0:.4e}")
        print(f"  max|S_t1| = {max_S_t1:.4e}  "
              f"(ratio to S_t0: {max_S_t1/max(max_S_t0,1e-30):.4e})")
        print(f"  max|S_t2| = {max_S_t2:.4e}  "
              f"(ratio to S_t0: {max_S_t2/max(max_S_t0,1e-30):.4e})")
        sys.stdout.flush()

    # ================================================================
    # Step 7: LOS integration
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
        print(f"\n--- Step 7: LOS C_l ({N_ell} ells, {N_snap} tau pts) ---")
        sys.stdout.flush()

    P_R = A_s * (k_arr / k_pivot) ** (n_s - 1.0)
    dtau_all = np.diff(tau_all)

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
    Delta_l_E = np.zeros((N_ell, N_k), dtype=np.float64)

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
            w = 0.5 * (dtau_all[max(0, it-1)]
                       + dtau_all[min(it, len(dtau_all)-1)])

        x_arr = k_arr * chi  # (N_k,)

        if use_gpu_bessel:
            jl_mx, jlp_mx = bessel_table.eval_jl_jlp(x_arr.astype(np.float32))
            mx.eval(jl_mx, jlp_mx)
            jl_all = np.array(jl_mx)     # (N_ell, N_k)
            jlp_all = np.array(jlp_mx)   # (N_ell, N_k)

            x_safe = np.where(np.abs(x_arr) < 1e-10, 1e-10, x_arr)
            inv_x = 1.0 / x_safe
            inv_x2 = inv_x ** 2

            for il in range(N_ell):
                jl = jl_all[il, :]
                jlp = jlp_all[il, :]
                ell = int(ell_values[il])

                # S_t2 radial function: R_l^(2) = (3*j_l'' + j_l)/2
                # j_l'' from Bessel ODE: j_l'' = (l(l+1)/x^2 - 1)*j_l - 2*j_l'/x
                if ell >= 2:
                    jlpp = (ell * (ell + 1) * inv_x2 - 1.0) * jl - 2.0 * inv_x * jlp
                    R_l2 = (3.0 * jlpp + jl) / 2.0
                else:
                    R_l2 = 0.0

                # TT transfer: S_t0*j_l + S_t1*j_l' + S_t2*R_l^(2)
                Delta_l[il, :] += w * (
                    S_t0[:, it] * jl
                    + S_t1[:, it] * jlp
                    + S_t2[:, it] * R_l2
                )

                # EE transfer: S_p * epsilon_l(x)
                if eps_prefactors[il] > 0:
                    eps_l = eps_prefactors[il] * jl * inv_x2
                    eps_l = np.where(np.abs(x_arr) < 1e-10, 0.0, eps_l)
                    Delta_l_E[il, :] += w * S_p[:, it] * eps_l
        else:
            from scipy.special import spherical_jn
            x_safe = np.where(np.abs(x_arr) < 1e-10, 1e-10, x_arr)
            inv_x = 1.0 / x_safe
            inv_x2 = inv_x ** 2

            for il, ell in enumerate(ell_values):
                ell_int = int(ell)
                jl = spherical_jn(ell_int, x_arr)
                jlp = spherical_jn(ell_int, x_arr, derivative=True)

                if ell_int >= 2:
                    jlpp = (ell_int*(ell_int+1)*inv_x2 - 1.0)*jl - 2.0*inv_x*jlp
                    R_l2 = (3.0 * jlpp + jl) / 2.0
                else:
                    R_l2 = 0.0

                Delta_l[il, :] += w * (
                    S_t0[:, it] * jl
                    + S_t1[:, it] * jlp
                    + S_t2[:, it] * R_l2
                )

                if eps_prefactors[il] > 0:
                    eps_l = eps_prefactors[il] * jl * inv_x2
                    eps_l = np.where(np.abs(x_arr) < 1e-10, 0.0, eps_l)
                    Delta_l_E[il, :] += w * S_p[:, it] * eps_l

    t_los = time.time() - t0
    if verbose:
        print(f"LOS integration (CLASS efficient): {t_los:.2f}s")

    # ================================================================
    # Step 8: C_l = 4 pi int dk/k P_R |Delta_l|^2
    # ================================================================
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
        'Cl': Cl_TT,
        'Dl': Dl_TT,
        'Cl_TT': Cl_TT,
        'Dl_TT': Dl_TT,
        'Cl_EE': Cl_EE,
        'Dl_EE': Dl_EE,
        'Cl_TE': Cl_TE,
        'Dl_TE': Dl_TE,
        'k_arr': k_arr,
        'bg': bg,
        'Delta_l': Delta_l,
        'Delta_l_E': Delta_l_E,
        'Theta2_arr': Theta2_arr,
        'P_pol_arr': P_pol_arr,
        'S_t0': S_t0,
        'S_t1': S_t1,
        'S_t2': S_t2,
        'tau_all': tau_all,
        'l_pol_max': l_pol_max,
        'timing': {
            'background': t_bg,
            'perturbations': t_pert,
            'los': t_los,
            'total': t_elapsed,
            'n_workers': n_workers,
        },
    }


# ============================================================================
# CLI entry point
# ============================================================================

def main():
    import argparse

    parser = argparse.ArgumentParser(
        description='Sync gauge solver with CLASS exact efficient source functions')
    parser.add_argument('--N_k', type=int, default=500)
    parser.add_argument('--method', type=str, default='BDF')
    parser.add_argument('--rtol', type=float, default=1e-6)
    parser.add_argument('--atol', type=float, default=1e-9)
    parser.add_argument('--lg_max', type=int, default=L_GAMMA_MAX)
    parser.add_argument('--ln_max', type=int, default=L_NU_MAX_SYNC)
    parser.add_argument('--l_pol_max', type=int, default=L_POL_MAX)
    parser.add_argument('--N_vis', type=int, default=500)
    parser.add_argument('--n_workers', type=int, default=None)
    args = parser.parse_args()

    result = run_sync_class_los(
        N_k=args.N_k, method=args.method,
        lg_max=args.lg_max, ln_max=args.ln_max,
        l_pol_max=args.l_pol_max,
        rtol=args.rtol, atol=args.atol,
        N_vis=args.N_vis,
        n_workers=args.n_workers,
        verbose=True)

    ell = result['ell']
    Dl_TT = result['Dl_TT']
    Dl_EE = result['Dl_EE']
    Dl_TE = result['Dl_TE']

    print("\n--- D_l values at key multipoles (CLASS EFFICIENT) ---")
    print(f"{'l':>6s}  {'D_l^TT':>12s}  {'D_l^EE':>12s}  {'D_l^TE':>12s}")
    print("-" * 50)
    for target_ell in [2, 10, 50, 100, 220, 537, 800, 1000, 1500, 2000]:
        idx = np.argmin(np.abs(ell - target_ell))
        print(f"  {ell[idx]:5d}  {Dl_TT[idx]:12.2f}  "
              f"{Dl_EE[idx]:12.4f}  {Dl_TE[idx]:12.2f}")

    # ================================================================
    # Compare with CLASS reference
    # ================================================================
    CLASS_FILE_TT = ('/Users/akaihuangm1/Desktop/github/gdm_class_public/'
                     'output/mlx_ref_lcdm_cl.dat')
    try:
        data = np.loadtxt(CLASS_FILE_TT)
        class_ell = data[:, 0].astype(int)
        class_Dl = data[:, 1] * (T_CMB * 1e6) ** 2

        mask = (class_ell >= 2) & (class_ell <= 2500)
        ell_use = class_ell[mask]
        Dl_ref = class_Dl[mask]
        f_test = interp1d(ell, Dl_TT, kind='cubic', fill_value='extrapolate')
        Dl_test_interp = f_test(ell_use)
        Dl_max = np.max(np.abs(Dl_ref))
        rel_err = (Dl_test_interp - Dl_ref) / Dl_max
        rms = np.sqrt(np.mean(rel_err ** 2))

        print(f"\n{'='*60}")
        print(f"  COMPARISON WITH CLASS REFERENCE (TT)")
        print(f"{'='*60}")
        print(f"  RMS error (l=2-2500): {rms*100:.3f}%")

        for target_l, name in [(220, "1st peak"), (537, "2nd peak"),
                                (810, "3rd peak"), (1120, "4th peak")]:
            idx_c = np.argmin(np.abs(class_ell - target_l))
            idx_o = np.argmin(np.abs(ell - target_l))
            ratio = Dl_TT[idx_o] / class_Dl[idx_c]
            print(f"  D_l({target_l:4d}) ratio (ours/CLASS): {ratio:.4f}"
                  f"  [{name}]")

        for l_lo, l_hi in [(2, 30), (30, 100), (100, 500),
                            (500, 1000), (1000, 2500)]:
            m = (ell_use >= l_lo) & (ell_use < l_hi)
            if np.sum(m) > 0:
                rms_bin = np.sqrt(np.mean(rel_err[m] ** 2))
                print(f"  RMS l=[{l_lo:4d},{l_hi:4d}): {rms_bin*100:.3f}%")

        print(f"{'='*60}")

    except Exception as e:
        print(f"\n[INFO] Could not compare with CLASS reference: {e}")

    outfile = 'cl_sync_class_los.dat'
    np.savetxt(outfile,
               np.column_stack([ell, Dl_TT, Dl_EE, Dl_TE]),
               header='ell  Dl_TT  Dl_EE  Dl_TE  [uK^2]',
               fmt=['%5d', '%15.6f', '%15.6f', '%15.6f'])
    print(f"\nSaved: {outfile}")


if __name__ == '__main__':
    main()
