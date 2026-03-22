"""
solver_sync_tca.py -- Sync gauge Boltzmann solver with CLASS TCA shear/slip.

This implements CLASS's Tight Coupling Approximation (TCA) properly:
  1. Two-phase integration: TCA phase (reduced, fast) + full hierarchy phase
  2. At the TCA->full switch:
     - F_gamma,2 seeded from TCA analytical shear (includes metric_shear)
     - F_gamma,3 seeded from second-order TCA
     - Polarization E_0..E_3 seeded from TCA equilibrium
  3. Full photon polarization hierarchy (E_0..E_10) with Thomson coupling
  4. Baryon sound speed c_s^2 k^2 delta_b

Includes a diagnostic mode that compares F_gamma,2 from ODE vs TCA analytical.

Target: TT RMS < 3% (down from 6.57%).

Usage:
  python -m mlx_class.solver_sync_tca [--N_k 500] [--parallel] [--diagnostic]

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
    H0_Mpc as _H0_MPC,
)
from .perturbations_neutrino import (
    _f_nu, _OMEGA_GAMMA, _OMEGA_NU,
)
from .perturbations_sync_tca import (
    make_sync_rhs_tca, adiabatic_ic_sync_tca, gauge_transform_tca,
    diagnose_metric_tca, seed_tca_at_switch, compute_tca_shear,
    compute_tca_slip, tca_diagnostic,
    n_var_sync_tca, idx_fn_start, idx_fg, idx_epol_start,
    IDX_ETA, IDX_DELTA_C, IDX_DELTA_B, IDX_THETA_B, IDX_FG_START,
    L_GAMMA_MAX, L_POL_MAX, L_NU_MAX_SYNC,
)
from .solver_sync import build_snapshot_grid


# ============================================================================
# TCA switch criterion (same as CLASS)
# ============================================================================

def find_tca_switch_time(bg, k, tau_c_over_tau_h_thr=0.015,
                         tau_c_over_tau_k_thr=0.01):
    """
    Find the TCA switch time using CLASS's criterion:
    TCA switches OFF when EITHER:
      tau_c / tau_h >= 0.015  (photon mean free path ~ Hubble time)
      tau_c / tau_k >= 0.01   (photon mean free path ~ mode wavelength)

    where tau_c = 1/|kappa_dot|, tau_h = 1/(aH) = 1/calH, tau_k = 1/k.
    """
    abs_kd = np.abs(bg.kappa_dot_grid)
    calH = bg.calH_grid
    tau_arr = bg.tau_grid

    # tau_c / tau_h = calH / |kd|
    # tau_c / tau_k = k / |kd|
    ratio_h = calH / np.maximum(abs_kd, 1e-30)
    ratio_k = k / np.maximum(abs_kd, 1e-30)

    # TCA is on when BOTH ratios are below threshold
    # TCA switches off at the FIRST time either exceeds threshold
    tca_off = (ratio_h >= tau_c_over_tau_h_thr) | (ratio_k >= tau_c_over_tau_k_thr)
    switch_idx = np.argmax(tca_off)

    if switch_idx == 0 and not tca_off[0]:
        # TCA never switches off (very small k)
        return float(tau_arr[-1])

    tau_switch = float(tau_arr[switch_idx])

    # Safety: don't switch too early
    tau_init = float(tau_arr[1])
    tau_switch = max(tau_switch, tau_init * 5)

    return tau_switch


# ============================================================================
# Solve a single k-mode with TCA seeding
# ============================================================================

def solve_single_k_tca(k, bg, tau_end, lg_max=L_GAMMA_MAX, lp_max=L_POL_MAX,
                       ln_max=L_NU_MAX_SYNC, method='Radau', rtol=1e-6,
                       atol=1e-9, use_tca_seeding=True, pol_approx='equilibrium'):
    """
    Solve the synchronous gauge Boltzmann equations for a single k-mode
    with CLASS-style TCA seeding at the switch point.

    Two-phase integration (when use_tca_seeding=True):
      Phase 1: tau_init -> tau_switch (full hierarchy, F_gamma,2 ~ 0 early)
      Phase 2: tau_switch -> tau_end  (full hierarchy with seeded F_gamma,2)

    When use_tca_seeding=False: single-shot integration, let Radau find the
    stiff equilibrium naturally. The polarization hierarchy still provides
    the correct Pi feedback to reduce effective Silk damping.
    """
    tau_init = float(bg.tau_grid[1])

    rhs_fn, nvar = make_sync_rhs_tca(k, bg, lg_max, lp_max, ln_max,
                                     pol_approx=pol_approx)
    y0 = adiabatic_ic_sync_tca(k, tau_init, lg_max, lp_max, ln_max)

    if not use_tca_seeding:
        # Single-shot: let Radau handle the stiff dynamics
        sol = solve_ivp(rhs_fn, [tau_init, tau_end], y0,
                        method=method, dense_output=True,
                        rtol=rtol, atol=atol)
        return sol

    # Two-phase with TCA seeding
    tau_switch = find_tca_switch_time(bg, k)
    tau_switch = min(tau_switch, tau_end * 0.9)

    if tau_switch > tau_init * 3 and tau_switch < tau_end * 0.8:
        # Phase 1: evolve to TCA switch point
        sol1 = solve_ivp(rhs_fn, [tau_init, tau_switch], y0,
                         method=method, dense_output=True,
                         rtol=rtol, atol=atol)
        if not sol1.success:
            # Fallback: no TCA seeding
            sol = solve_ivp(rhs_fn, [tau_init, tau_end], y0,
                            method=method, dense_output=True,
                            rtol=rtol, atol=atol)
            return sol

        # Get state at switch point
        y_switch = sol1.sol(tau_switch).copy()

        # Compute background at switch
        abs_kd_switch = float(np.interp(tau_switch, bg.tau_grid,
                                         np.abs(bg.kappa_dot_grid)))
        calH_switch = float(np.interp(tau_switch, bg.tau_grid, bg.calH_grid))
        a_switch = float(np.interp(tau_switch, bg.tau_grid, bg.a_grid))
        R_switch = float(np.interp(tau_switch, bg.tau_grid, bg.R_grid))

        # Diagnose metric at switch
        h_prime_sw, eta_prime_sw = diagnose_metric_tca(
            y_switch, k, calH_switch, a_switch, lg_max, lp_max, ln_max)

        # Seed TCA values (includes metric_shear, l=3, polarization)
        seed_tca_at_switch(y_switch, k, abs_kd_switch, calH_switch,
                           a_switch, lg_max, lp_max,
                           h_prime_sw, eta_prime_sw, R_switch)

        # Phase 2: evolve from switch to end
        sol2 = solve_ivp(rhs_fn, [tau_switch, tau_end], y_switch,
                         method=method, dense_output=True,
                         rtol=rtol, atol=atol)

        # Combine into a single dense output
        class CombinedSolution:
            def __init__(self, sa, sb, ts):
                self.success = sb.success
                self.message = sb.message
                self.t_switch = ts
                self._sa = sa
                self._sb = sb

            def sol(self, t):
                if t <= self.t_switch:
                    return self._sa.sol(t)
                else:
                    return self._sb.sol(t)

        combined = CombinedSolution(sol1, sol2, tau_switch)
        return combined
    else:
        # No TCA phase or switch too late; solve in one shot
        sol = solve_ivp(rhs_fn, [tau_init, tau_end], y0,
                        method=method, dense_output=True,
                        rtol=rtol, atol=atol)
        return sol


# ============================================================================
# Multiprocessing worker
# ============================================================================

def _solve_single_k_tca_worker(args):
    """Worker for multiprocessing. Returns (N_snap, 6) array or None."""
    import numpy as np
    from scipy.integrate import solve_ivp

    (k, bg_arrays, tau_end, lg_max, lp_max, ln_max, method, rtol, atol,
     use_tca_seeding, pol_approx, tau_all, a_snap, calH_snap) = args

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

    from mlx_class.perturbations_sync_tca import (
        make_sync_rhs_tca, adiabatic_ic_sync_tca, gauge_transform_tca,
        diagnose_metric_tca, seed_tca_at_switch,
        n_var_sync_tca, idx_fn_start, idx_fg, idx_epol_start,
        IDX_ETA, IDX_THETA_B, IDX_FG_START,
    )
    from mlx_class.solver_sync_tca import find_tca_switch_time

    tau_init = float(bg_lite.tau_grid[1])

    # Find TCA switch
    tau_switch = find_tca_switch_time(bg_lite, k)
    tau_switch = min(tau_switch, tau_end * 0.9)

    rhs_fn, nvar = make_sync_rhs_tca(k, bg_lite, lg_max, lp_max, ln_max,
                                     pol_approx=pol_approx)
    y0 = adiabatic_ic_sync_tca(k, tau_init, lg_max, lp_max, ln_max)

    # Integration: either two-phase with TCA seeding or single-shot
    if not use_tca_seeding:
        # Single-shot: let Radau handle stiff equilibrium
        sol = solve_ivp(rhs_fn, [tau_init, tau_end], y0,
                        method=method, dense_output=True,
                        rtol=rtol, atol=atol)
        if not sol.success:
            return None
        sol_combined = sol
    elif tau_switch > tau_init * 3 and tau_switch < tau_end * 0.8:
        sol1 = solve_ivp(rhs_fn, [tau_init, tau_switch], y0,
                         method=method, dense_output=True,
                         rtol=rtol, atol=atol)
        if not sol1.success:
            sol = solve_ivp(rhs_fn, [tau_init, tau_end], y0,
                            method=method, dense_output=True,
                            rtol=rtol, atol=atol)
            if not sol.success:
                return None
            sol_combined = sol
        else:
            y_switch = sol1.sol(tau_switch).copy()

            abs_kd_sw = float(np.interp(tau_switch, bg_lite.tau_grid,
                                         np.abs(bg_lite.kappa_dot_grid)))
            calH_sw = float(np.interp(tau_switch, bg_lite.tau_grid,
                                       bg_lite.calH_grid))
            a_sw = float(np.interp(tau_switch, bg_lite.tau_grid,
                                    bg_lite.a_grid))
            R_sw = float(np.interp(tau_switch, bg_lite.tau_grid,
                                    bg_lite.R_grid))

            h_prime_sw, eta_prime_sw = diagnose_metric_tca(
                y_switch, k, calH_sw, a_sw, lg_max, lp_max, ln_max)

            seed_tca_at_switch(y_switch, k, abs_kd_sw, calH_sw, a_sw,
                               lg_max, lp_max, h_prime_sw, eta_prime_sw, R_sw)

            sol2 = solve_ivp(rhs_fn, [tau_switch, tau_end], y_switch,
                             method=method, dense_output=True,
                             rtol=rtol, atol=atol)
            if not sol2.success:
                return None

            class _Combined:
                def __init__(self, s1, s2, ts):
                    self.success = s2.success
                    self.t_switch = ts
                    self._s1 = s1
                    self._s2 = s2
                def sol(self, t):
                    return self._s1.sol(t) if t <= self.t_switch else self._s2.sol(t)

            sol_combined = _Combined(sol1, sol2, tau_switch)
    else:
        sol = solve_ivp(rhs_fn, [tau_init, tau_end], y0,
                        method=method, dense_output=True,
                        rtol=rtol, atol=atol)
        if not sol.success:
            return None
        sol_combined = sol

    # Evaluate at snapshots
    N_snap = len(tau_all)
    results = np.zeros((N_snap, 6), dtype=np.float64)
    ep_s = idx_epol_start(lg_max)

    for it in range(N_snap):
        tau = tau_all[it]
        y = sol_combined.sol(tau)
        a = a_snap[it]
        calH = calH_snap[it]

        h_prime, eta_prime = diagnose_metric_tca(
            y, k, calH, a, lg_max, lp_max, ln_max)
        gt = gauge_transform_tca(
            y, k, calH, a, lg_max, lp_max, ln_max, h_prime, eta_prime,
            pol_approx=pol_approx)

        eta_val = y[IDX_ETA]
        results[it, 0] = gt['Phi_N']
        results[it, 1] = gt['Psi_N']
        results[it, 2] = y[IDX_FG_START] / 4.0 - eta_val + gt['Phi_N']
        results[it, 3] = y[IDX_THETA_B] / k + (h_prime + 6*eta_prime) / (2*k)
        results[it, 4] = y[IDX_FG_START + 2] / 4.0
        results[it, 5] = gt['Pi'] / 4.0

    return results


# ============================================================================
# Diagnostic: compare F_gamma,2 from ODE vs TCA analytical
# ============================================================================

def run_tca_diagnostic(bg, k_values=None, verbose=True):
    """
    Diagnostic: for several k-modes, compare F_gamma,2(ODE) vs F_gamma,2(TCA)
    at several tau values around recombination.

    This tells us whether the Radau solver correctly finds the stiff equilibrium
    for F_gamma,2 during tight coupling, or whether explicit TCA seeding is needed.
    """
    if k_values is None:
        k_values = [0.01, 0.02, 0.05, 0.1, 0.2]

    lg_max, lp_max, ln_max = L_GAMMA_MAX, L_POL_MAX, L_NU_MAX_SYNC

    # Tau sampling: before recombination through visibility peak
    tau_test_values = np.linspace(50, 350, 30)

    if verbose:
        print("\n" + "=" * 75)
        print("TCA DIAGNOSTIC: F_gamma,2 (ODE) vs F_gamma,2 (TCA analytical)")
        print("=" * 75)
        print(f"{'k':>8s}  {'tau':>8s}  {'|kd|/k':>10s}  "
              f"{'Fg2_ODE':>12s}  {'Fg2_TCA':>12s}  {'ratio':>8s}  {'status':>8s}")
        print("-" * 75)

    all_results = {}

    for k in k_values:
        tau_end = float(tau_test_values[-1]) + 10.0
        sol = solve_single_k_tca(k, bg, tau_end, lg_max, lp_max, ln_max,
                                 method='Radau', rtol=1e-7, atol=1e-10)

        if not sol.success:
            if verbose:
                print(f"  k={k:.3f}: ODE FAILED")
            continue

        k_results = []
        for tau in tau_test_values:
            y = sol.sol(float(tau))
            a = float(bg.a_at_tau(np.array([tau]))[0])
            calH = float(bg.calH_at_tau(np.array([tau]))[0])
            abs_kd = float(np.interp(tau, bg.tau_grid,
                                      np.abs(bg.kappa_dot_grid)))

            diag = tca_diagnostic(y, k, calH, a, lg_max, lp_max, ln_max, abs_kd)
            k_results.append(diag)

            ratio_kd_k = abs_kd / k
            in_tca = ratio_kd_k > 50  # well within TCA

            if verbose and (abs(tau - 100) < 10 or abs(tau - 200) < 10
                           or abs(tau - 280) < 10 or abs(tau - 300) < 10):
                status = "TCA" if in_tca else "TRANS" if ratio_kd_k > 5 else "FULL"
                print(f"  {k:8.4f}  {tau:8.1f}  {ratio_kd_k:10.1f}  "
                      f"{diag['Fg2_ode']:12.4e}  {diag['Fg2_tca']:12.4e}  "
                      f"{diag['ratio']:8.3f}  {status:>8s}")

        all_results[k] = k_results

    if verbose:
        print("-" * 75)
        print("\nInterpretation:")
        print("  ratio ~ 1.0: ODE resolves TCA equilibrium correctly (seeding NOT needed)")
        print("  ratio != 1.0: ODE misses TCA equilibrium (seeding IS needed)")
        print("  During TCA (|kd|/k > 50), F_gamma,2 should match TCA analytical.")
        print("  At transition (|kd|/k ~ 5-50), some deviation is expected.")
        print()

    return all_results


# ============================================================================
# Full pipeline (parallel)
# ============================================================================

def run_sync_solver_tca_parallel(N_k=500, k_min=3e-4, k_max=0.35,
                                 method='Radau',
                                 lg_max=L_GAMMA_MAX, lp_max=L_POL_MAX,
                                 ln_max=L_NU_MAX_SYNC,
                                 rtol=1e-6, atol=1e-9, n_workers=None,
                                 use_tca_seeding=True, pol_approx='equilibrium',
                                 verbose=True):
    """
    Parallel synchronous gauge pipeline with CLASS TCA seeding + polarization.
    """
    t_total = time.time()

    if n_workers is None:
        n_workers = os.cpu_count() or 4

    # ================================================================
    # Step 1: Background
    # ================================================================
    if verbose:
        print("=" * 65)
        print("Synchronous Gauge Boltzmann Solver (CLASS TCA)")
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

    tau_vis, tau_late, tau_all = build_snapshot_grid(
        bg, N_vis=60, N_early_isw=30, N_late_isw=20)
    N_snap = len(tau_all)

    nvar = n_var_sync_tca(lg_max, lp_max, ln_max)

    if verbose:
        print(f"\n--- Step 2: Perturbations ({N_k} k-modes, {method}, "
              f"{n_workers} workers) ---")
        print(f"k range: [{k_min:.1e}, {k_max:.2f}] Mpc^-1")
        print(f"Snapshot grid: {N_snap} tau points")
        print(f"  Visibility: {len(tau_vis)} pts, "
              f"tau=[{tau_vis[0]:.1f}, {tau_vis[-1]:.1f}]")
        print(f"  Late ISW:   {len(tau_late)} pts, "
              f"tau=[{tau_late[0]:.1f}, {tau_late[-1]:.1f}]")
        if hasattr(bg, 'z_reio') and bg.tau_reio > 0:
            print(f"  Reionization: z_reio={bg.z_reio:.2f}, "
                  f"tau_reio={bg.tau_reio:.4f}")
        print(f"Hierarchy: lg_max={lg_max}, lp_max={lp_max}, "
              f"ln_max={ln_max}, nvar={nvar}")
        print(f"TCA seeding: {'ON' if use_tca_seeding else 'OFF'}")
        print(f"Polarization: {pol_approx}")

        # Show TCA switch times for a few k-modes
        if use_tca_seeding:
            pass  # print below
        for ktest in [0.01, 0.05, 0.1, 0.3]:
            tau_sw = find_tca_switch_time(bg, ktest)
            abs_kd_sw = float(np.interp(tau_sw, bg.tau_grid,
                                         np.abs(bg.kappa_dot_grid)))
            print(f"  TCA switch: k={ktest:.3f} -> tau={tau_sw:.1f} "
                  f"(|kd|={abs_kd_sw:.0f})")
        sys.stdout.flush()

    # ================================================================
    # Step 3: Solve all k-modes in parallel
    # ================================================================
    t0 = time.time()

    a_snap = np.asarray(bg.a_at_tau(tau_all), dtype=np.float64)
    calH_snap = np.asarray(bg.calH_at_tau(tau_all), dtype=np.float64)

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
        (k, bg_arrays, tau_end, lg_max, lp_max, ln_max, method, rtol, atol,
         use_tca_seeding, pol_approx, tau_all, a_snap, calH_snap)
        for k in k_arr
    ]

    if verbose:
        print(f"Dispatching {N_k} k-modes to {n_workers} workers...")
        sys.stdout.flush()

    with Pool(processes=n_workers) as pool:
        worker_results = pool.map(_solve_single_k_tca_worker, work_items)

    t_pert = time.time() - t0

    # ================================================================
    # Collect results
    # ================================================================
    Phi_N = np.zeros((N_k, N_snap))
    Psi_N = np.zeros((N_k, N_snap))
    Theta0_N = np.zeros((N_k, N_snap))
    vb_N = np.zeros((N_k, N_snap))
    Theta2_arr = np.zeros((N_k, N_snap))
    Pi4_arr = np.zeros((N_k, N_snap))

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
        Pi4_arr[ik] = res[:, 5]

    if verbose:
        print(f"Perturbations: {t_pert:.1f}s ({n_failed} failed, "
              f"{N_k - n_failed} OK)")

    # ================================================================
    # Step 4: Phi' + Psi' via cubic spline
    # ================================================================
    PhiPsi = Phi_N + Psi_N
    PhiPsi_prime = np.zeros_like(PhiPsi)
    dtau_all = np.diff(tau_all)
    for ik in range(N_k):
        cs = CubicSpline(tau_all, PhiPsi[ik])
        PhiPsi_prime[ik] = cs(tau_all, 1)

    # ================================================================
    # Step 5: LOS integration -> C_l (TT + EE + TE)
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

    # Transfer functions
    Delta_l = np.zeros((N_ell, N_k), dtype=np.float64)
    Delta_l_E = np.zeros((N_ell, N_k), dtype=np.float64)

    # E-mode source uses the full Pi from polarization hierarchy
    # S_E = g(tau) * (3/4) * Pi_Theta where Pi_Theta = Pi4_arr
    eps_prefactors = np.zeros(N_ell)
    for il, ell in enumerate(ell_values):
        l = int(ell)
        if l >= 2:
            eps_prefactors[il] = np.sqrt(float((l - 1) * l * (l + 1) * (l + 2)))

    if verbose:
        print(f"Computing LOS transfer functions (TT + EE with TCA + pol)...")
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

        # TT sources
        S_SW = g_snap[it] * (Theta0_N[:, it] + Psi_N[:, it])
        S_Dop = g_snap[it] * vb_N[:, it]
        S_ISW = exp_neg_kappa[it] * PhiPsi_prime[:, it]

        # EE source using full Pi
        S_E = g_snap[it] * 0.75 * Pi4_arr[:, it]

        if use_gpu_bessel:
            x_arr = k_arr * chi
            jl_mx, jlp_mx = bessel_table.eval_jl_jlp(
                x_arr.astype(np.float32))
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
                Delta_l[il, :] += w * (
                    S_SW * jl + S_Dop * jlp + S_ISW * jl)
                if eps_prefactors[il] > 0:
                    eps_l = eps_prefactors[il] * jl * inv_x2
                    eps_l = np.where(np.abs(x) < 1e-10, 0.0, eps_l)
                    Delta_l_E[il, :] += w * S_E * eps_l

    t_los = time.time() - t0
    if verbose:
        print(f"LOS integration: {t_los:.2f}s")

    # ================================================================
    # Step 6: C_l
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
        'Phi_N': Phi_N,
        'Psi_N': Psi_N,
        'Theta2_arr': Theta2_arr,
        'Pi4_arr': Pi4_arr,
        'timing': {
            'background': t_bg,
            'perturbations': t_pert,
            'los': t_los,
            'total': t_elapsed,
            'n_workers': n_workers,
        },
    }


# ============================================================================
# Main
# ============================================================================

def main():
    import argparse

    parser = argparse.ArgumentParser(description='Sync gauge solver with CLASS TCA')
    parser.add_argument('--N_k', type=int, default=500)
    parser.add_argument('--parallel', action='store_true', default=True)
    parser.add_argument('--diagnostic', action='store_true',
                        help='Run TCA diagnostic before full solve')
    parser.add_argument('--no-plot', action='store_true')
    parser.add_argument('--no-tca-seed', action='store_true',
                        help='Disable TCA seeding (no phase switch)')
    parser.add_argument('--pol-approx', default='equilibrium',
                        choices=['full', 'equilibrium', 'none'],
                        help='Polarization approximation mode')
    args = parser.parse_args()

    # ================================================================
    # Run TCA diagnostic first (if requested)
    # ================================================================
    if args.diagnostic:
        print("\n*** Running TCA diagnostic ***\n")
        bg = Background(khronon=False, recombination='recfast')
        bg.solve()
        run_tca_diagnostic(bg, k_values=[0.01, 0.02, 0.05, 0.1, 0.2])
        print()

    # ================================================================
    # Full pipeline
    # ================================================================
    use_tca_seeding = not args.no_tca_seed
    pol_approx = args.pol_approx
    print(f"*** TCA seeding: {'ON' if use_tca_seeding else 'OFF'}, "
          f"Polarization: {pol_approx} ***\n")

    result = run_sync_solver_tca_parallel(
        N_k=args.N_k, use_tca_seeding=use_tca_seeding,
        pol_approx=pol_approx, verbose=True)

    ell = result['ell']
    Dl_TT = result['Dl_TT']
    Dl_EE = result['Dl_EE']
    Dl_TE = result['Dl_TE']

    print("\n--- D_l values at key multipoles ---")
    print(f"{'l':>6s}  {'D_l^TT':>12s}  {'D_l^EE':>12s}  {'D_l^TE':>12s}")
    print("-" * 50)
    for target_ell in [2, 10, 50, 100, 220, 500, 800, 1000, 1500, 2000]:
        idx = np.argmin(np.abs(ell - target_ell))
        print(f"  {ell[idx]:5d}  {Dl_TT[idx]:12.2f}  "
              f"{Dl_EE[idx]:12.4f}  {Dl_TE[idx]:12.2f}")

    # --- Save output ---
    outdir = os.path.dirname(os.path.abspath(__file__))
    outfile = os.path.join(outdir, 'cl_sync_tca.dat')
    np.savetxt(outfile,
               np.column_stack([ell, Dl_TT, Dl_EE, Dl_TE]),
               header='ell  Dl_TT  Dl_EE  Dl_TE',
               fmt='%8d %15.6e %15.6e %15.6e')
    print(f"\nSaved: {outfile}")

    # --- Compare with CLASS ---
    class_file = os.path.join(outdir, 'class_comparison.dat')
    if os.path.exists(class_file):
        class_data = np.loadtxt(class_file)
        ell_class = class_data[:, 0]
        Dl_class = class_data[:, 1]

        from scipy.interpolate import interp1d
        our_interp = interp1d(ell.astype(float), Dl_TT, kind='cubic',
                              fill_value='extrapolate')
        Dl_ours = our_interp(ell_class)

        # RMS in 50 < l < 2000
        mask = (ell_class >= 50) & (ell_class <= 2000) & (Dl_class > 0)
        if np.sum(mask) > 10:
            ratio = Dl_ours[mask] / Dl_class[mask]
            rms = np.sqrt(np.mean((ratio - 1.0)**2)) * 100
            print(f"\n--- Comparison with CLASS (50 < l < 2000) ---")
            print(f"TT RMS error: {rms:.2f}%")

            # Per-bin analysis
            bins = [(50, 200, "low-l"), (200, 500, "peaks 1-2"),
                    (500, 1000, "peaks 3-5"), (1000, 2000, "damping")]
            for lo, hi, name in bins:
                bmask = mask & (ell_class >= lo) & (ell_class <= hi)
                if np.sum(bmask) > 5:
                    bratio = Dl_ours[bmask] / Dl_class[bmask]
                    brms = np.sqrt(np.mean((bratio - 1.0)**2)) * 100
                    bmean = (np.mean(bratio) - 1) * 100
                    print(f"  {name:>12s} (l={lo}-{hi}): "
                          f"RMS={brms:.2f}%, mean bias={bmean:+.2f}%")

            # Peak analysis
            from scipy.signal import argrelmax
            peak_mask = (ell_class > 100) & (ell_class < 1200)
            ell_pr = ell_class[peak_mask]
            Dl_pr = Dl_class[peak_mask]
            ours_pr = Dl_ours[peak_mask]

            peak_idx = argrelmax(Dl_pr, order=20)[0]
            if len(peak_idx) >= 3:
                print(f"\nPeak analysis:")
                for ip in range(min(len(peak_idx), 4)):
                    pi = peak_idx[ip]
                    l_p = ell_pr[pi]
                    class_val = Dl_pr[pi]
                    ours_val = ours_pr[pi]
                    ratio_p = ours_val / class_val
                    print(f"  Peak {ip+1}: l={l_p:.0f}, "
                          f"CLASS={class_val:.1f}, "
                          f"ours={ours_val:.1f}, "
                          f"ratio={ratio_p:.3f}")

        # Also compare with v1 (old solver without TCA)
        v1_file = os.path.join(outdir, 'cl_sync_gauge.dat')
        if os.path.exists(v1_file):
            v1_data = np.loadtxt(v1_file)
            if v1_data.shape[1] >= 2:
                ell_v1 = v1_data[:, 0]
                Dl_v1 = v1_data[:, 1]
                v1_interp = interp1d(ell_v1, Dl_v1, kind='cubic',
                                     fill_value='extrapolate')
                Dl_v1_at_class = v1_interp(ell_class[mask])
                v1_ratio = Dl_v1_at_class / Dl_class[mask]
                v1_rms = np.sqrt(np.mean((v1_ratio - 1.0)**2)) * 100
                print(f"\nImprovement:")
                print(f"  v1 (no TCA):       {v1_rms:.2f}%")
                print(f"  TCA (this solver): {rms:.2f}%")

        # Compare with v2 physics (if exists)
        v2_file = os.path.join(outdir, 'cl_sync_v2_physics.dat')
        if os.path.exists(v2_file):
            v2_data = np.loadtxt(v2_file)
            if v2_data.shape[1] >= 2:
                ell_v2 = v2_data[:, 0]
                Dl_v2 = v2_data[:, 1]
                v2_interp = interp1d(ell_v2, Dl_v2, kind='cubic',
                                     fill_value='extrapolate')
                Dl_v2_at_class = v2_interp(ell_class[mask])
                v2_ratio = Dl_v2_at_class / Dl_class[mask]
                v2_rms = np.sqrt(np.mean((v2_ratio - 1.0)**2)) * 100
                print(f"  v2 (pol only):     {v2_rms:.2f}%")

    # --- Plot ---
    if not args.no_plot:
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt

            fig, axes = plt.subplots(2, 1, figsize=(12, 10),
                                     gridspec_kw={'height_ratios': [3, 1]})

            # Top panel: D_l^TT
            ax1 = axes[0]
            ax1.plot(ell, Dl_TT, 'b-', linewidth=1.5, label='TCA solver')

            if os.path.exists(class_file):
                ax1.plot(ell_class, Dl_class, 'k--', linewidth=1,
                         label='CLASS', alpha=0.7)

            v1_file = os.path.join(outdir, 'cl_sync_gauge.dat')
            if os.path.exists(v1_file):
                v1_data = np.loadtxt(v1_file)
                ax1.plot(v1_data[:, 0], v1_data[:, 1], 'r:', linewidth=1,
                         label='v1 (no TCA)', alpha=0.5)

            ax1.set_xlabel('l')
            ax1.set_ylabel(r'$D_l^{TT}$ [$\mu K^2$]')
            ax1.set_title('CMB TT Power Spectrum: CLASS TCA Shear + Polarization')
            ax1.legend()
            ax1.set_xlim(2, 2500)
            ax1.set_ylim(0, None)

            # Bottom panel: ratio to CLASS
            ax2 = axes[1]
            if os.path.exists(class_file):
                our_at_class = our_interp(ell_class)
                mask_plot = (ell_class >= 2) & (ell_class <= 2500) & (Dl_class > 0)
                ax2.plot(ell_class[mask_plot],
                         our_at_class[mask_plot] / Dl_class[mask_plot],
                         'b-', linewidth=1, label='TCA/CLASS')

                if os.path.exists(v1_file):
                    v1_data = np.loadtxt(v1_file)
                    v1_interp = interp1d(v1_data[:, 0], v1_data[:, 1],
                                         kind='cubic', fill_value='extrapolate')
                    v1_at_class = v1_interp(ell_class)
                    ax2.plot(ell_class[mask_plot],
                             v1_at_class[mask_plot] / Dl_class[mask_plot],
                             'r:', linewidth=1, label='v1/CLASS', alpha=0.5)

                ax2.axhline(y=1.0, color='k', linestyle='-', alpha=0.3)
                ax2.axhline(y=0.97, color='g', linestyle='--', alpha=0.3)
                ax2.axhline(y=1.03, color='g', linestyle='--', alpha=0.3,
                            label='+/-3%')
                ax2.set_ylim(0.85, 1.15)
                ax2.set_xlabel('l')
                ax2.set_ylabel('ratio')
                ax2.legend()

            plt.tight_layout()
            figfile = os.path.join(outdir, 'sync_tca_vs_class.png')
            plt.savefig(figfile, dpi=150)
            plt.close()
            print(f"\nPlot saved: {figfile}")

        except Exception as e:
            print(f"Plotting failed: {e}")


if __name__ == '__main__':
    main()

