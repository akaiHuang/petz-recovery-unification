"""
solver_sync_v2_physics.py -- Synchronous gauge Boltzmann solver with polarization.

Uses perturbations_sync_v2.py which adds THREE missing physics pieces:
  1. Photon E-mode polarization hierarchy (E_0..E_10)
  2. TCA shear seeding at the switch point
  3. Baryon sound speed c_s^2 k^2 delta_b

Target: TT RMS < 3% (down from 6.57%), second peak ratio > 0.85.

Two pipelines:
  run_sync_solver_v2()          -- sequential (for debugging)
  run_sync_solver_v2_parallel() -- multiprocessing (production)

Usage:
  python -m mlx_class.solver_sync_v2_physics [--parallel] [--N_k 500]

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
from .perturbations_sync_v2 import (
    make_sync_rhs_v2, adiabatic_ic_sync_v2, gauge_transform_v2,
    diagnose_metric_v2, seed_tca_shear,
    n_var_sync_v2, idx_fn_start, idx_fg, idx_epol_start,
    IDX_ETA, IDX_DELTA_C, IDX_DELTA_B, IDX_THETA_B, IDX_FG_START,
    L_GAMMA_MAX, L_POL_MAX, L_NU_MAX_SYNC,
)

# Also import the old solver's snapshot grid builder
from .solver_sync import build_snapshot_grid


# ============================================================================
# Solve a single k-mode with TCA seeding
# ============================================================================

def solve_single_k_v2(k, bg, tau_end, lg_max=L_GAMMA_MAX, lp_max=L_POL_MAX,
                      ln_max=L_NU_MAX_SYNC, method='Radau', rtol=1e-6,
                      atol=1e-9, tca_switch_kd=30.0):
    """
    Solve the synchronous gauge Boltzmann equations for a single k-mode
    WITH polarization hierarchy and TCA shear seeding.

    Two-phase integration:
      Phase 1: tau_init -> tau_switch (TCA: F_gamma,2 = E_l = 0, evolving normally)
      Phase 2: tau_switch -> tau_end  (full hierarchy with seeded shear)

    The switch happens when |kappa_dot| / k drops below tca_switch_kd,
    i.e., when the Thomson scattering rate is only ~30x the mode frequency.

    Parameters
    ----------
    k : float, wavenumber in Mpc^{-1}
    bg : Background object (solved)
    tau_end : float, end conformal time
    lg_max, lp_max, ln_max : hierarchy truncations
    method : str, ODE method
    rtol, atol : tolerances
    tca_switch_kd : float, |kd|/k threshold for switching

    Returns
    -------
    sol : OdeSolution object with dense_output
    """
    tau_init = bg.tau_grid[1]

    # Find TCA switch time: when |kappa_dot| / k < tca_switch_kd
    kd_grid = np.abs(bg.kappa_dot_grid)
    ratio_grid = kd_grid / k
    # Find the first time ratio drops below threshold
    below = np.where(ratio_grid < tca_switch_kd)[0]
    if len(below) > 0 and bg.tau_grid[below[0]] > tau_init * 2:
        tau_switch = float(bg.tau_grid[below[0]])
        # Clamp: don't switch too early or too late
        tau_switch = max(tau_switch, tau_init * 5)
        tau_switch = min(tau_switch, tau_end * 0.5)
    else:
        # Mode is always in the perturbative regime; no TCA phase
        tau_switch = tau_init

    rhs_fn, nvar = make_sync_rhs_v2(k, bg, lg_max, lp_max, ln_max)
    y0 = adiabatic_ic_sync_v2(k, tau_init, lg_max, lp_max, ln_max)

    if tau_switch > tau_init * 3:
        # Phase 1: evolve to switch point
        sol1 = solve_ivp(rhs_fn, [tau_init, tau_switch], y0,
                         method=method, dense_output=True,
                         rtol=rtol, atol=atol)
        if not sol1.success:
            # Fall back to no TCA seeding
            sol = solve_ivp(rhs_fn, [tau_init, tau_end], y0,
                            method=method, dense_output=True,
                            rtol=rtol, atol=atol)
            return sol

        # Get state at switch point
        y_switch = sol1.sol(tau_switch)

        # Seed TCA shear values
        abs_kd_switch = float(np.interp(tau_switch, bg.tau_grid,
                                         np.abs(bg.kappa_dot_grid)))
        seed_tca_shear(y_switch, k, abs_kd_switch, lg_max, lp_max)

        # Phase 2: evolve from switch to end
        sol2 = solve_ivp(rhs_fn, [tau_switch, tau_end], y_switch,
                         method=method, dense_output=True,
                         rtol=rtol, atol=atol)

        # Combine the two solutions into a single dense_output wrapper
        class CombinedSolution:
            """Wraps two OdeSolution objects into one with dense_output."""
            def __init__(self, sol_a, sol_b, t_switch):
                self.success = sol_b.success
                self.message = sol_b.message
                self.t_switch = t_switch
                self.sol_a = sol_a
                self.sol_b = sol_b

            def sol(self, t):
                if t <= self.t_switch:
                    return self.sol_a.sol(t)
                else:
                    return self.sol_b.sol(t)

        combined = CombinedSolution(sol1, sol2, tau_switch)
        combined.success = sol2.success
        return combined
    else:
        # No TCA phase needed; solve in one shot
        sol = solve_ivp(rhs_fn, [tau_init, tau_end], y0,
                        method=method, dense_output=True,
                        rtol=rtol, atol=atol)
        return sol


# ============================================================================
# Multiprocessing worker
# ============================================================================

def _solve_single_k_v2_worker(args):
    """
    Worker function for multiprocessing. Solves one k-mode with polarization.

    Returns a numpy array of shape (N_snap, 6) with columns:
      [Phi_N, Psi_N, Theta0_N, vb_N, Theta2, Pi/4]
    or None if the ODE solve fails.
    """
    import numpy as np
    from scipy.integrate import solve_ivp

    (k, bg_arrays, tau_end, lg_max, lp_max, ln_max, method, rtol, atol,
     tca_switch_kd, tau_all, a_snap, calH_snap) = args

    # Reconstruct lightweight background
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

    from mlx_class.perturbations_sync_v2 import (
        make_sync_rhs_v2, adiabatic_ic_sync_v2, gauge_transform_v2,
        diagnose_metric_v2, seed_tca_shear,
        n_var_sync_v2, idx_fn_start, idx_fg, idx_epol_start,
        IDX_ETA, IDX_THETA_B, IDX_FG_START,
    )

    tau_init = bg_lite.tau_grid[1]

    # Find TCA switch time
    kd_grid = np.abs(bg_lite.kappa_dot_grid)
    ratio_grid = kd_grid / k
    below = np.where(ratio_grid < tca_switch_kd)[0]
    if len(below) > 0 and bg_lite.tau_grid[below[0]] > tau_init * 2:
        tau_switch = float(bg_lite.tau_grid[below[0]])
        tau_switch = max(tau_switch, tau_init * 5)
        tau_switch = min(tau_switch, tau_end * 0.5)
    else:
        tau_switch = tau_init

    rhs_fn, nvar = make_sync_rhs_v2(k, bg_lite, lg_max, lp_max, ln_max)
    y0 = adiabatic_ic_sync_v2(k, tau_init, lg_max, lp_max, ln_max)

    # Two-phase integration with TCA seeding
    if tau_switch > tau_init * 3:
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
            y_switch = sol1.sol(tau_switch)
            abs_kd_switch = float(np.interp(
                tau_switch, bg_lite.tau_grid, np.abs(bg_lite.kappa_dot_grid)))
            seed_tca_shear(y_switch, k, abs_kd_switch, lg_max, lp_max)

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
                    if t <= self.t_switch:
                        return self._s1.sol(t)
                    else:
                        return self._s2.sol(t)

            sol_combined = _Combined(sol1, sol2, tau_switch)
    else:
        sol = solve_ivp(rhs_fn, [tau_init, tau_end], y0,
                        method=method, dense_output=True,
                        rtol=rtol, atol=atol)
        if not sol.success:
            return None
        sol_combined = sol

    # Evaluate at snapshots and gauge transform
    N_snap = len(tau_all)
    results = np.zeros((N_snap, 6), dtype=np.float64)

    for it in range(N_snap):
        tau = tau_all[it]
        y = sol_combined.sol(tau)
        a = a_snap[it]
        calH = calH_snap[it]

        h_prime, eta_prime = diagnose_metric_v2(
            y, k, calH, a, lg_max, lp_max, ln_max)
        gt = gauge_transform_v2(
            y, k, calH, a, lg_max, lp_max, ln_max, h_prime, eta_prime)

        eta_val = y[IDX_ETA]
        results[it, 0] = gt['Phi_N']
        results[it, 1] = gt['Psi_N']
        # Gauge-invariant Theta0_N
        results[it, 2] = y[IDX_FG_START] / 4.0 - eta_val + gt['Phi_N']
        # Newtonian gauge baryon velocity
        results[it, 3] = y[IDX_THETA_B] / k + (h_prime + 6 * eta_prime) / (2 * k)
        # Photon quadrupole: F_gamma,2 / 4
        results[it, 4] = y[IDX_FG_START + 2] / 4.0
        # Polarization anisotropy Pi/4 for E-mode source
        results[it, 5] = gt['Pi'] / 4.0

    return results


# ============================================================================
# Full pipeline (parallel version)
# ============================================================================

def run_sync_solver_v2_parallel(N_k=500, k_min=3e-4, k_max=0.35, method='Radau',
                                lg_max=L_GAMMA_MAX, lp_max=L_POL_MAX,
                                ln_max=L_NU_MAX_SYNC,
                                rtol=1e-6, atol=1e-9, n_workers=None,
                                tca_switch_kd=30.0, verbose=True):
    """
    Parallel synchronous gauge pipeline with polarization.
    """
    t_total = time.time()

    if n_workers is None:
        n_workers = os.cpu_count() or 4

    # ================================================================
    # Step 1: Background
    # ================================================================
    if verbose:
        print("=" * 65)
        print("Synchronous Gauge Boltzmann Solver v2 (WITH POLARIZATION)")
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

    nvar = n_var_sync_v2(lg_max, lp_max, ln_max)

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
        print(f"TCA switch: |kd|/k < {tca_switch_kd}")
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
         tca_switch_kd, tau_all, a_snap, calH_snap)
        for k in k_arr
    ]

    if verbose:
        print(f"Dispatching {N_k} k-modes to {n_workers} workers...")
        sys.stdout.flush()

    with Pool(processes=n_workers) as pool:
        worker_results = pool.map(_solve_single_k_v2_worker, work_items)

    t_pert = time.time() - t0

    # ================================================================
    # Collect results (6 columns: Phi_N, Psi_N, Theta0_N, vb_N, Theta2, Pi/4)
    # ================================================================
    Phi_N = np.zeros((N_k, N_snap))
    Psi_N = np.zeros((N_k, N_snap))
    Theta0_N = np.zeros((N_k, N_snap))
    vb_N = np.zeros((N_k, N_snap))
    Theta2_arr = np.zeros((N_k, N_snap))
    Pi4_arr = np.zeros((N_k, N_snap))   # Pi/4 for E-mode source

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

    # E-mode source: g(tau) * (3/4) * Pi/4 * 4 = g(tau) * (3/4) * Pi
    # Pi = F_gamma,2 + E_0 + E_2 is the full polarization anisotropy.
    # In the old code, we approximated Pi ~ alpha_P * F_gamma,2 (with alpha_P = 1.7).
    # Now we have the actual Pi from the polarization hierarchy.
    # S_E = g(tau) * (3/4) * Pi = g(tau) * (3/4) * 4 * Pi4_arr = 3 * g(tau) * Pi4_arr
    # Wait: Pi4_arr = Pi/4, so (3/4)*Pi = (3/4)*4*Pi4_arr = 3*Pi4_arr
    # Actually the E-mode LOS source is:
    # S_E(k,tau) = g(tau) * (3/4) * Pi(k,tau) where Pi = Theta_2 + Theta_P0 + Theta_P2
    # In our conventions: Theta_2 = F_g2/4, Theta_P0 = E_0/4, Theta_P2 = E_2/4
    # So Pi_Theta = Theta_2 + Theta_P0 + Theta_P2 = (F_g2 + E_0 + E_2)/4 = Pi4_arr
    # And S_E = g(tau) * (3/4) * Pi_Theta = g(tau) * (3/4) * Pi4_arr
    # The extra factor of 4 from F -> Theta is already absorbed in Pi4_arr.

    eps_prefactors = np.zeros(N_ell)
    for il, ell in enumerate(ell_values):
        l = int(ell)
        if l >= 2:
            eps_prefactors[il] = np.sqrt(float((l - 1) * l * (l + 1) * (l + 2)))

    if verbose:
        print(f"Computing LOS transfer functions (TT + EE with polarization)...")
        sys.stdout.flush()

    for it in range(N_snap):
        chi = chi_snap[it]
        if chi <= 0:
            continue

        # Trapezoidal weight
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

        # EE source: uses the full Pi from polarization hierarchy
        # S_E = g(tau) * (3/4) * Pi_Theta where Pi_Theta = Pi4_arr
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
# Full pipeline (sequential version, for debugging)
# ============================================================================

def run_sync_solver_v2(N_k=300, k_min=3e-4, k_max=0.35, method='Radau',
                       lg_max=L_GAMMA_MAX, lp_max=L_POL_MAX,
                       ln_max=L_NU_MAX_SYNC,
                       rtol=1e-6, atol=1e-9, tca_switch_kd=30.0,
                       verbose=True):
    """
    Sequential synchronous gauge pipeline with polarization.
    """
    t_total = time.time()

    # ================================================================
    # Step 1: Background
    # ================================================================
    if verbose:
        print("=" * 65)
        print("Synchronous Gauge Boltzmann Solver v2 (WITH POLARIZATION)")
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

    nvar = n_var_sync_v2(lg_max, lp_max, ln_max)

    if verbose:
        print(f"\n--- Step 2: Perturbations ({N_k} k-modes, {method}) ---")
        print(f"k range: [{k_min:.1e}, {k_max:.2f}] Mpc^-1")
        print(f"Snapshot grid: {N_snap} tau points")
        print(f"Hierarchy: lg_max={lg_max}, lp_max={lp_max}, "
              f"ln_max={ln_max}, nvar={nvar}")
        print(f"TCA switch: |kd|/k < {tca_switch_kd}")
        sys.stdout.flush()

    # ================================================================
    # Step 3: Solve all k-modes
    # ================================================================
    t0 = time.time()

    Phi_N = np.zeros((N_k, N_snap))
    Psi_N = np.zeros((N_k, N_snap))
    Theta0_N = np.zeros((N_k, N_snap))
    vb_N = np.zeros((N_k, N_snap))
    Theta2_arr = np.zeros((N_k, N_snap))
    Pi4_arr = np.zeros((N_k, N_snap))

    a_snap = bg.a_at_tau(tau_all)
    calH_snap = bg.calH_at_tau(tau_all)

    n_failed = 0
    for ik, k in enumerate(k_arr):
        sol = solve_single_k_v2(k, bg, tau_all[-1] + 1.0,
                                lg_max, lp_max, ln_max,
                                method, rtol, atol, tca_switch_kd)
        if not sol.success:
            n_failed += 1
            if verbose and n_failed <= 3:
                print(f"  WARNING: k={k:.4e} failed")
            continue

        for it in range(N_snap):
            tau = tau_all[it]
            y = sol.sol(tau)
            calH = calH_snap[it]
            a = a_snap[it]

            h_prime, eta_prime = diagnose_metric_v2(
                y, k, calH, a, lg_max, lp_max, ln_max)
            gt = gauge_transform_v2(
                y, k, calH, a, lg_max, lp_max, ln_max,
                h_prime, eta_prime)

            eta_val = y[IDX_ETA]
            Phi_N[ik, it] = gt['Phi_N']
            Psi_N[ik, it] = gt['Psi_N']
            Theta0_N[ik, it] = y[IDX_FG_START] / 4.0 - eta_val + gt['Phi_N']
            vb_N[ik, it] = (y[IDX_THETA_B] / k
                            + (h_prime + 6 * eta_prime) / (2 * k))
            Theta2_arr[ik, it] = y[IDX_FG_START + 2] / 4.0
            Pi4_arr[ik, it] = gt['Pi'] / 4.0

        if verbose and (ik + 1) % 20 == 0:
            elapsed = time.time() - t0
            eta_est = elapsed / (ik + 1) * N_k
            print(f"  [{ik+1}/{N_k}] elapsed={elapsed:.1f}s, "
                  f"ETA={eta_est:.1f}s")
            sys.stdout.flush()

    t_pert = time.time() - t0
    if verbose:
        print(f"Perturbations: {t_pert:.1f}s ({n_failed} failed)")

    # ================================================================
    # Step 4-6: Same as parallel version
    # ================================================================
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
        print(f"\n--- Step 5: LOS C_l ({N_ell} ells) ---")
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

    eps_prefactors = np.zeros(N_ell)
    for il, ell in enumerate(ell_values):
        l = int(ell)
        if l >= 2:
            eps_prefactors[il] = np.sqrt(float((l - 1) * l * (l + 1) * (l + 2)))

    if verbose:
        print(f"Computing LOS transfer functions (TT + EE with polarization)...")
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

        S_SW = g_snap[it] * (Theta0_N[:, it] + Psi_N[:, it])
        S_Dop = g_snap[it] * vb_N[:, it]
        S_ISW = exp_neg_kappa[it] * PhiPsi_prime[:, it]
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

    # C_l computation
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
        print(f"\nTotal time: {t_elapsed:.1f}s")

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
        },
    }


# ============================================================================
# CLI with comparison to CLASS and v1
# ============================================================================

def main():
    import argparse

    parser = argparse.ArgumentParser(
        description='Synchronous gauge Boltzmann solver v2 (with polarization)')
    parser.add_argument('--N_k', type=int, default=500)
    parser.add_argument('--method', type=str, default='Radau')
    parser.add_argument('--rtol', type=float, default=1e-6)
    parser.add_argument('--atol', type=float, default=1e-9)
    parser.add_argument('--lg_max', type=int, default=L_GAMMA_MAX)
    parser.add_argument('--lp_max', type=int, default=L_POL_MAX)
    parser.add_argument('--ln_max', type=int, default=L_NU_MAX_SYNC)
    parser.add_argument('--parallel', action='store_true', default=True,
                        help='Use multiprocessing (default: True)')
    parser.add_argument('--sequential', action='store_true',
                        help='Force sequential mode')
    parser.add_argument('--n_workers', type=int, default=None)
    parser.add_argument('--tca_switch', type=float, default=30.0,
                        help='TCA switch threshold |kd|/k (default: 30)')
    parser.add_argument('--no-plot', action='store_true',
                        help='Skip plotting')
    args = parser.parse_args()

    if args.sequential:
        result = run_sync_solver_v2(
            N_k=args.N_k, method=args.method,
            lg_max=args.lg_max, lp_max=args.lp_max, ln_max=args.ln_max,
            rtol=args.rtol, atol=args.atol,
            tca_switch_kd=args.tca_switch,
            verbose=True)
    else:
        result = run_sync_solver_v2_parallel(
            N_k=args.N_k, method=args.method,
            lg_max=args.lg_max, lp_max=args.lp_max, ln_max=args.ln_max,
            rtol=args.rtol, atol=args.atol,
            n_workers=args.n_workers,
            tca_switch_kd=args.tca_switch,
            verbose=True)

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
    outfile = os.path.join(outdir, 'cl_sync_v2_physics.dat')
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

        # Interpolate our result to CLASS ell values
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

            # Peak analysis
            from scipy.signal import argrelmax
            peak_mask = (ell_class > 100) & (ell_class < 1200)
            ell_peak_range = ell_class[peak_mask]
            Dl_peak_range = Dl_class[peak_mask]
            ours_peak_range = Dl_ours[peak_mask]

            # Find CLASS peaks
            peak_idx = argrelmax(Dl_peak_range, order=20)[0]
            if len(peak_idx) >= 3:
                print(f"\nPeak analysis:")
                for ip in range(min(len(peak_idx), 4)):
                    pi = peak_idx[ip]
                    l_p = ell_peak_range[pi]
                    class_val = Dl_peak_range[pi]
                    ours_val = ours_peak_range[pi]
                    ratio_p = ours_val / class_val
                    print(f"  Peak {ip+1}: l={l_p:.0f}, "
                          f"CLASS={class_val:.1f}, "
                          f"ours={ours_val:.1f}, "
                          f"ratio={ratio_p:.3f}")

    # --- Compare with v1 (old solver) ---
    v1_file = os.path.join(outdir, 'cl_sync_gauge.dat')
    if os.path.exists(v1_file):
        v1_data = np.loadtxt(v1_file)
        if v1_data.shape[1] >= 2:
            ell_v1 = v1_data[:, 0]
            Dl_v1 = v1_data[:, 1]
            print(f"\n--- v1 vs v2 comparison ---")
            v1_interp = interp1d(ell_v1, Dl_v1, kind='cubic',
                                 fill_value='extrapolate')
            Dl_v1_at_class = v1_interp(ell_class[mask])
            v1_ratio = Dl_v1_at_class / Dl_class[mask]
            v1_rms = np.sqrt(np.mean((v1_ratio - 1.0)**2)) * 100
            print(f"v1 TT RMS error: {v1_rms:.2f}%")
            print(f"v2 TT RMS error: {rms:.2f}%")
            print(f"Improvement: {v1_rms:.2f}% -> {rms:.2f}%")

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
            ax1.plot(ell, Dl_TT, 'b-', linewidth=1.5, label='v2 (polarization)')

            if os.path.exists(class_file):
                ax1.plot(ell_class, Dl_class, 'k--', linewidth=1,
                         label='CLASS', alpha=0.7)

            if os.path.exists(v1_file):
                ax1.plot(ell_v1, Dl_v1, 'r:', linewidth=1,
                         label='v1 (no polarization)', alpha=0.5)

            ax1.set_xlabel('l')
            ax1.set_ylabel(r'$D_l^{TT}$ [$\mu K^2$]')
            ax1.set_title('CMB TT Power Spectrum: v2 with Polarization Feedback')
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
                         'b-', linewidth=1, label='v2/CLASS')

                if os.path.exists(v1_file):
                    v1_at_class = v1_interp(ell_class)
                    ax2.plot(ell_class[mask_plot],
                             v1_at_class[mask_plot] / Dl_class[mask_plot],
                             'r:', linewidth=1, label='v1/CLASS', alpha=0.5)

                ax2.axhline(1.0, color='k', linestyle='--', alpha=0.3)
                ax2.axhline(0.97, color='gray', linestyle=':', alpha=0.3)
                ax2.axhline(1.03, color='gray', linestyle=':', alpha=0.3)
                ax2.set_xlim(2, 2500)
                ax2.set_ylim(0.85, 1.15)
                ax2.set_xlabel('l')
                ax2.set_ylabel('ratio to CLASS')
                ax2.legend()

            plt.tight_layout()
            plot_file = os.path.join(outdir, 'cl_sync_v2_physics.png')
            plt.savefig(plot_file, dpi=150)
            print(f"Plot saved: {plot_file}")
            plt.close()

        except ImportError:
            print("[WARNING] matplotlib not available, skipping plot")


if __name__ == '__main__':
    main()
