"""
solver_hybrid.py -- Two-phase optimized Boltzmann solver for CMB C_l.

KEY INSIGHT: In the full Boltzmann solve, ~90% of scipy's time is spent
on the late-time free-streaming phase (tau > 350 Mpc, ~14000 Radau steps
per high-k mode). But this late evolution is only needed for the late ISW
effect, which is significant only at low k (< 0.01 Mpc^{-1}) and low ell
(< ~30). By splitting the solve and only evolving low-k modes through the
late phase, we skip ~90% of the expensive high-k late-time steps.

ARCHITECTURE:

  Phase 1: scipy Radau (per k-mode, parallel)  tau_init -> tau_switch
    - Covers tight coupling + recombination + visibility peak (ALL k-modes).
    - Returns gauge-transformed LOS source functions at visibility snapshots.
    - Also returns state at tau_switch for Phase 2 continuation.
    - tau_switch = visibility FWHM right edge + 50 Mpc (~337 Mpc).
    - Cost: ~25-45s with 10 workers (500 modes, ~0.5-1s avg per mode).

  Phase 2: scipy Radau (per k-mode, parallel)  tau_switch -> tau_end
    - ONLY for low-k modes (k < k_isw_max = 0.01 Mpc^-1) where late ISW
      matters for ell < 30.
    - High-k modes skip this: j_l(k*chi) oscillates rapidly when k*chi >> ell,
      so the ISW integral averages to zero. Potentials are extrapolated as
      constant from the last Phase 1 value to avoid CubicSpline derivative
      artifacts.
    - Cost: ~3-10s with 10 workers (~250 low-k modes, ~0.2s each).

  Phase 3: LOS integration with GPU Bessel functions -> C_l.
    - Same LOS formula as solver_sync.py (SW + Doppler + ISW).
    - Cost: ~1-2s.

PERFORMANCE (500 k-modes, M1 Max, 10 workers):
  Hybrid:  35-60s total   (7% RMS vs full scipy)
  scipy:   260-310s total (reference)
  Speedup: ~5-9x

ACCURACY: 7% RMS vs full scipy (same physics, same sync gauge).
The only approximation is constant-potential extrapolation for high-k modes
at late times, which is physically exact to O(Phi_dot * Delta_tau_ISW).

GAUGE: Synchronous gauge, CDM rest frame (Ma & Bertschinger 1995).

Usage:
  python -m mlx_class.solver_hybrid [--N_k 500] [--n_workers 10] [--plot]

Author: Sheng-Kai Huang, 2026
"""
import os
os.environ['MPLCONFIGDIR'] = '/tmp/mpl_config'

import numpy as np
import time
import sys
from multiprocessing import Pool

from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

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
    make_sync_rhs, adiabatic_ic_sync, gauge_transform, diagnose_metric,
    n_var_sync, idx_fn_start, idx_fg,
    IDX_ETA, IDX_DELTA_C, IDX_DELTA_B, IDX_THETA_B, IDX_FG_START,
    L_GAMMA_MAX, L_NU_MAX_SYNC,
)
from .solver_sync import build_snapshot_grid
from .solver_magnus import compare_with_reference


# ============================================================================
# Phase 1 worker: scipy through visibility peak (module-level for pickling)
# ============================================================================

def _phase1_worker(args):
    """
    Phase 1 worker: solve one k-mode from tau_init to tau_switch.
    Returns state at tau_switch + gauge-transformed quantities at early snapshots.
    """
    import numpy as np
    from scipy.integrate import solve_ivp

    (k, bg_arrays, tau_init, tau_switch, lg_max, ln_max,
     rtol, atol, snap_tau, snap_a, snap_calH) = args

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
        make_sync_rhs as _make_rhs,
        adiabatic_ic_sync as _aic,
        gauge_transform as _gt,
        diagnose_metric as _dm,
        IDX_ETA as _ETA, IDX_FG_START as _FG, IDX_THETA_B as _TB,
    )

    rhs_fn, nvar = _make_rhs(k, bg_lite, lg_max, ln_max)
    y0 = _aic(k, tau_init, lg_max, ln_max)

    sol = solve_ivp(rhs_fn, [tau_init, tau_switch + 1.0], y0,
                    method='Radau', dense_output=True,
                    rtol=rtol, atol=atol)

    if not sol.success:
        return None

    y_switch = sol.sol(tau_switch).copy()

    N_snap = len(snap_tau)
    results = np.zeros((N_snap, 5), dtype=np.float64)

    for it in range(N_snap):
        tau = snap_tau[it]
        if tau > tau_switch:
            break
        y = sol.sol(tau)
        a = snap_a[it]
        calH = snap_calH[it]

        h_prime, eta_prime = _dm(y, k, calH, a, lg_max, ln_max)
        gt = _gt(y, k, calH, a, lg_max, ln_max, h_prime, eta_prime)

        eta_val = y[_ETA]
        results[it, 0] = gt['Phi_N']
        results[it, 1] = gt['Psi_N']
        results[it, 2] = y[_FG] / 4.0 - eta_val + gt['Phi_N']
        results[it, 3] = y[_TB] / k + (h_prime + 6 * eta_prime) / (2 * k)
        results[it, 4] = y[_FG + 2] / 4.0

    return (y_switch, results)


# ============================================================================
# Phase 2 worker: late ISW for low-k modes (module-level for pickling)
# ============================================================================

def _phase2_worker(args):
    """
    Phase 2 worker: solve one low-k mode from tau_switch to tau_end.
    Returns gauge-transformed quantities at late snapshot points.
    """
    import numpy as np
    from scipy.integrate import solve_ivp

    (k, bg_arrays, y_switch, tau_switch, tau_end, lg_max, ln_max,
     rtol, atol, late_snap_tau, late_snap_a, late_snap_calH) = args

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
        make_sync_rhs as _make_rhs,
        gauge_transform as _gt,
        diagnose_metric as _dm,
        IDX_ETA as _ETA, IDX_FG_START as _FG, IDX_THETA_B as _TB,
    )

    rhs_fn, nvar = _make_rhs(k, bg_lite, lg_max, ln_max)

    sol = solve_ivp(rhs_fn, [tau_switch, tau_end + 1.0], y_switch,
                    method='Radau', dense_output=True,
                    rtol=rtol, atol=atol)

    if not sol.success:
        return None

    N_late = len(late_snap_tau)
    results = np.zeros((N_late, 5), dtype=np.float64)

    for it in range(N_late):
        tau = late_snap_tau[it]
        y = sol.sol(tau)
        a = late_snap_a[it]
        calH = late_snap_calH[it]

        h_prime, eta_prime = _dm(y, k, calH, a, lg_max, ln_max)
        gt = _gt(y, k, calH, a, lg_max, ln_max, h_prime, eta_prime)

        eta_val = y[_ETA]
        results[it, 0] = gt['Phi_N']
        results[it, 1] = gt['Psi_N']
        results[it, 2] = y[_FG] / 4.0 - eta_val + gt['Phi_N']
        results[it, 3] = y[_TB] / k + (h_prime + 6 * eta_prime) / (2 * k)
        results[it, 4] = y[_FG + 2] / 4.0

    return results


def find_tca_switch_tau(bg, tau_offset=50.0):
    """
    Find the scipy Phase 1 -> Phase 2 switch time.
    Switch well after the visibility peak so Phase 1 captures all acoustic physics.
    """
    peak_idx = np.argmax(bg.visibility_grid)
    tau_rec = bg.tau_grid[peak_idx]
    g_peak = bg.visibility_grid[peak_idx]

    # Find visibility FWHM right edge
    half_max = g_peak / 2.0
    right = peak_idx
    while right < len(bg.visibility_grid) - 1 and bg.visibility_grid[right] > half_max:
        right += 1
    tau_vis_hi = bg.tau_grid[right]

    tau_switch = tau_vis_hi + tau_offset
    tau_switch = min(tau_switch, bg.tau_0 * 0.5)
    return float(tau_switch)


# ============================================================================
# Main hybrid solver
# ============================================================================

def solve_hybrid(N_k=500, k_min=3e-4, k_max=0.35,
                 lg_max=L_GAMMA_MAX, ln_max=L_NU_MAX_SYNC,
                 tau_offset=50.0, k_isw_max=0.01,
                 rtol=1e-6, atol=1e-9,
                 n_workers=None,
                 verbose=True):
    """
    Hybrid solver for CMB power spectrum.

    Phase 1: scipy Radau for all k-modes through visibility peak.
    Phase 2: scipy Radau for low-k modes only through late ISW.
    Phase 3: LOS integration -> C_l.

    Parameters
    ----------
    N_k : int
        Number of k-modes.
    tau_offset : float
        Mpc after visibility FWHM right edge to place the Phase 1/2 switch.
    k_isw_max : float
        Max wavenumber for Phase 2 late ISW solve (Mpc^-1).
    rtol, atol : float
        scipy tolerances.
    n_workers : int or None
        Number of parallel workers. Default: cpu_count.
    verbose : bool
        Print progress.
    """
    t_total = time.time()

    if n_workers is None:
        n_workers = os.cpu_count() or 4

    # ==================================================================
    # Step 1: Background
    # ==================================================================
    if verbose:
        print("=" * 70)
        print("Hybrid Boltzmann Solver (scipy Phase 1+2)")
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
    nvar = n_var_sync(lg_max, ln_max)

    tau_init = bg.tau_grid[1]
    tau_switch = find_tca_switch_tau(bg, tau_offset=tau_offset)
    tau_end = bg.tau_0 * 0.98

    # Build LOS snapshot grid
    tau_vis, tau_late, tau_snap = build_snapshot_grid(
        bg, N_vis=60, N_early_isw=30, N_late_isw=20)
    N_snap = len(tau_snap)

    # Split snapshots into early (Phase 1) and late (Phase 2)
    mask_early = tau_snap <= tau_switch
    mask_late = tau_snap > tau_switch
    n_early_snap = int(np.sum(mask_early))
    n_late_snap = int(np.sum(mask_late))
    early_indices = np.where(mask_early)[0]
    late_indices = np.where(mask_late)[0]

    late_snap_tau = tau_snap[mask_late]
    late_snap_a = np.asarray(bg.a_at_tau(late_snap_tau), dtype=np.float64)
    late_snap_calH = np.asarray(bg.calH_at_tau(late_snap_tau), dtype=np.float64)

    # Identify low-k modes for Phase 2
    isw_mask = k_arr < k_isw_max
    n_isw_modes = int(np.sum(isw_mask))
    isw_k_indices = np.where(isw_mask)[0]

    a_snap = np.asarray(bg.a_at_tau(tau_snap), dtype=np.float64)
    calH_snap = np.asarray(bg.calH_at_tau(tau_snap), dtype=np.float64)

    t_setup = time.time() - t0

    if verbose:
        peak_idx = np.argmax(bg.visibility_grid)
        tau_rec = bg.tau_grid[peak_idx]
        print(f"N_k={N_k}, hierarchy: lg={lg_max}, ln={ln_max}, nvar={nvar}")
        print(f"Visibility peak: tau_rec = {tau_rec:.1f} Mpc")
        print(f"Phase 1: tau=[{tau_init:.4f}, {tau_switch:.1f}] Mpc "
              f"(ALL {N_k} modes)")
        print(f"Phase 2: tau=[{tau_switch:.1f}, {tau_end:.1f}] Mpc "
              f"(low-k only: {n_isw_modes} modes with k < {k_isw_max})")
        print(f"LOS snapshots: {N_snap} ({n_early_snap} early + {n_late_snap} late)")
        print(f"Setup: {t_setup:.2f}s")

    # ==================================================================
    # Step 3: Phase 1 -- scipy through visibility peak (all k)
    # ==================================================================
    if verbose:
        print(f"\n--- Step 3: Phase 1 ({N_k} modes, {n_workers} workers) ---")
        sys.stdout.flush()

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

    work_items_p1 = [
        (k, bg_arrays, tau_init, tau_switch, lg_max, ln_max,
         rtol, atol, tau_snap[:n_early_snap], a_snap[:n_early_snap],
         calH_snap[:n_early_snap])
        for k in k_arr
    ]

    with Pool(processes=n_workers) as pool:
        p1_results = pool.map(_phase1_worker, work_items_p1)

    t_phase1 = time.time() - t0

    # Collect Phase 1 results
    Phi_N = np.zeros((N_k, N_snap), dtype=np.float64)
    Psi_N = np.zeros((N_k, N_snap), dtype=np.float64)
    Theta0_N = np.zeros((N_k, N_snap), dtype=np.float64)
    vb_N = np.zeros((N_k, N_snap), dtype=np.float64)
    Theta2_arr = np.zeros((N_k, N_snap), dtype=np.float64)

    y_at_switch = np.zeros((N_k, nvar), dtype=np.float64)

    n_failed_p1 = 0
    for ik, res in enumerate(p1_results):
        if res is None:
            n_failed_p1 += 1
            y_at_switch[ik] = adiabatic_ic_sync(
                k_arr[ik], tau_switch, lg_max, ln_max)
            if verbose and n_failed_p1 <= 3:
                print(f"  WARNING: k={k_arr[ik]:.4e} Phase 1 failed")
        else:
            y_switch, snap_results = res
            y_at_switch[ik] = y_switch
            for it_e in range(n_early_snap):
                snap_idx = early_indices[it_e]
                Phi_N[ik, snap_idx] = snap_results[it_e, 0]
                Psi_N[ik, snap_idx] = snap_results[it_e, 1]
                Theta0_N[ik, snap_idx] = snap_results[it_e, 2]
                vb_N[ik, snap_idx] = snap_results[it_e, 3]
                Theta2_arr[ik, snap_idx] = snap_results[it_e, 4]

    if verbose:
        print(f"Phase 1: {t_phase1:.2f}s ({n_failed_p1} failed)")

    # ==================================================================
    # Step 4: Phase 2 -- late ISW for low-k modes only
    # ==================================================================
    if verbose:
        print(f"\n--- Step 4: Phase 2 ({n_isw_modes} low-k modes, "
              f"{n_workers} workers) ---")
        sys.stdout.flush()

    t0 = time.time()

    if n_isw_modes > 0 and n_late_snap > 0:
        work_items_p2 = [
            (k_arr[ik], bg_arrays, y_at_switch[ik], tau_switch, tau_end,
             lg_max, ln_max, rtol, atol,
             late_snap_tau, late_snap_a, late_snap_calH)
            for ik in isw_k_indices
        ]

        with Pool(processes=n_workers) as pool:
            p2_results = pool.map(_phase2_worker, work_items_p2)

        n_failed_p2 = 0
        for i_isw, res in enumerate(p2_results):
            ik = isw_k_indices[i_isw]
            if res is None:
                n_failed_p2 += 1
                if verbose and n_failed_p2 <= 3:
                    print(f"  WARNING: k={k_arr[ik]:.4e} Phase 2 failed")
            else:
                for it_l in range(n_late_snap):
                    snap_idx = late_indices[it_l]
                    Phi_N[ik, snap_idx] = res[it_l, 0]
                    Psi_N[ik, snap_idx] = res[it_l, 1]
                    Theta0_N[ik, snap_idx] = res[it_l, 2]
                    vb_N[ik, snap_idx] = res[it_l, 3]
                    Theta2_arr[ik, snap_idx] = res[it_l, 4]
    else:
        n_failed_p2 = 0

    t_phase2 = time.time() - t0

    if verbose:
        print(f"Phase 2: {t_phase2:.2f}s ({n_failed_p2} failed)")

    # For high-k modes (k >= k_isw_max) that skipped Phase 2, the late
    # snapshots still have Phi_N = Psi_N = 0. We need smooth values to
    # avoid CubicSpline derivative spikes at the early/late boundary.
    # Extrapolate the last early-snapshot values to late snapshots.
    # After recombination, potentials are ~constant for subhorizon modes,
    # and the visibility-weighted ISW contribution is negligible anyway.
    if n_early_snap > 0 and n_late_snap > 0:
        last_early_idx = early_indices[-1]
        for ik in range(N_k):
            if not isw_mask[ik]:
                # High-k mode: extrapolate Phi_N, Psi_N from last early value
                Phi_N[ik, late_indices] = Phi_N[ik, last_early_idx]
                Psi_N[ik, late_indices] = Psi_N[ik, last_early_idx]
                # Theta0_N and vb_N don't matter (g(tau) ~ 0 at late times)
                # but keep them smooth too
                Theta0_N[ik, late_indices] = Theta0_N[ik, last_early_idx]
                vb_N[ik, late_indices] = vb_N[ik, last_early_idx]
                Theta2_arr[ik, late_indices] = Theta2_arr[ik, last_early_idx]

    # Clean up NaN/Inf
    for arr in [Phi_N, Psi_N, Theta0_N, vb_N, Theta2_arr]:
        np.nan_to_num(arr, copy=False, nan=0.0, posinf=0.0, neginf=0.0)

    # ==================================================================
    # Step 5: Phi' + Psi' by cubic spline
    # ==================================================================
    PhiPsi = Phi_N + Psi_N
    PhiPsi_prime = np.zeros_like(PhiPsi)
    for ik in range(N_k):
        if np.all(np.isfinite(PhiPsi[ik])):
            cs = CubicSpline(tau_snap, PhiPsi[ik])
            PhiPsi_prime[ik] = cs(tau_snap, 1)

    # ==================================================================
    # Step 6: LOS integration -> C_l
    # ==================================================================
    if verbose:
        print(f"\n--- Step 5: LOS integration ---")
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
            for il, ell_val in enumerate(ell_values):
                jl = spherical_jn(int(ell_val), x)
                jlp = spherical_jn(int(ell_val), x, derivative=True)
                Delta_l[il, :] += w * (S_SW * jl + S_Dop * jlp + S_ISW * jl)
                if eps_prefactors[il] > 0:
                    eps_l = eps_prefactors[il] * jl * inv_x2
                    eps_l = np.where(np.abs(x) < 1e-10, 0.0, eps_l)
                    Delta_l_E[il, :] += w * S_E * eps_l

    t_los = time.time() - t0
    if verbose:
        print(f"LOS integration: {t_los:.2f}s")

    # ==================================================================
    # Step 7: C_l
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

    if verbose:
        print(f"\n{'='*70}")
        print(f"TIMING SUMMARY")
        print(f"{'='*70}")
        print(f"  Background:          {t_bg:.2f}s")
        print(f"  Setup:               {t_setup:.2f}s")
        print(f"  Phase 1 (all k):     {t_phase1:.2f}s  "
              f"({N_k} modes, {n_workers} workers)")
        print(f"  Phase 2 (low-k ISW): {t_phase2:.2f}s  "
              f"({n_isw_modes} modes, {n_workers} workers)")
        print(f"  LOS integration:     {t_los:.2f}s")
        print(f"  C_l sum:             {t_cl:.3f}s")
        print(f"  TOTAL:               {t_elapsed:.2f}s")
        print(f"{'='*70}")

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
            'phase1': t_phase1,
            'phase2': t_phase2,
            'los': t_los,
            'total': t_elapsed,
            'n_workers': n_workers,
            'n_isw_modes': n_isw_modes,
        },
    }


# ============================================================================
# Plot
# ============================================================================

def _plot_result(result, cmp=None):
    """Plot D_l^TT and comparison."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 1, figsize=(12, 8),
                              gridspec_kw={'height_ratios': [3, 1]})

    ell = result['ell']
    Dl = result['Dl_TT']

    axes[0].plot(ell, Dl, 'b-', lw=1.5, label='Hybrid solver')
    if cmp is not None:
        axes[0].plot(cmp['ell_cmp'], cmp['Dl_ref'], 'r--', lw=1,
                     label=f"Reference (RMS={cmp['rms']*100:.1f}%)")
    axes[0].set_ylabel(r'$D_\ell^{TT}$ [$\mu K^2$]')
    axes[0].set_xlim(2, 2500)
    axes[0].set_ylim(0, None)
    axes[0].legend()
    axes[0].set_title('Hybrid Boltzmann Solver (scipy Phase 1+2)')

    if cmp is not None:
        frac = (cmp['Dl_us'] - cmp['Dl_ref']) / (cmp['Dl_ref'] + 1e-30)
        axes[1].plot(cmp['ell_cmp'], frac * 100, 'k-', lw=0.5)
        axes[1].axhline(0, color='gray', ls='--')
        axes[1].set_ylabel('Residual [%]')
        axes[1].set_xlim(2, 2500)
        axes[1].set_ylim(-50, 50)

    axes[1].set_xlabel(r'Multipole $\ell$')
    plt.tight_layout()
    outfile = os.path.join(os.path.dirname(__file__), 'cl_hybrid.png')
    plt.savefig(outfile, dpi=150)
    print(f"Plot saved to {outfile}")
    plt.close()


# ============================================================================
# CLI main
# ============================================================================

def main():
    """Run the hybrid solver."""
    import argparse
    parser = argparse.ArgumentParser(
        description='Hybrid Boltzmann solver')
    parser.add_argument('--N_k', type=int, default=500,
                        help='Number of k-modes')
    parser.add_argument('--tau_offset', type=float, default=50.0,
                        help='Mpc after visibility FWHM to place switch')
    parser.add_argument('--k_isw_max', type=float, default=0.01,
                        help='Max k for late ISW Phase 2 (Mpc^-1)')
    parser.add_argument('--n_workers', type=int, default=None,
                        help='Number of parallel workers')
    parser.add_argument('--plot', action='store_true',
                        help='Plot results')
    args = parser.parse_args()

    result = solve_hybrid(
        N_k=args.N_k,
        tau_offset=args.tau_offset,
        k_isw_max=args.k_isw_max,
        n_workers=args.n_workers,
        verbose=True,
    )

    cmp = compare_with_reference(result, verbose=True)

    outfile = os.path.join(os.path.dirname(__file__), 'cl_hybrid.dat')
    np.savetxt(outfile,
               np.column_stack([result['ell'], result['Dl_TT'],
                                result['Dl_EE'], result['Dl_TE']]),
               header='ell  Dl_TT  Dl_EE  Dl_TE',
               fmt='%8d %15.6e %15.6e %15.6e')
    print(f"\nResults saved to {outfile}")

    if args.plot:
        _plot_result(result, cmp)


if __name__ == '__main__':
    main()
