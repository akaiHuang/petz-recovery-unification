"""
solver_sync_efficient.py -- Sync gauge solver with CLASS-style efficient source functions.

The key difference from solver_sync.py is the LOS integration:
This uses integration by parts on the Doppler term to eliminate the
j_l' oscillatory integral, following CLASS (perturbations.c lines 7518-7563).

solver_sync.py already did one IBP (on exp(-kappa)*k*Psi), giving:
  S_t0 = g*(Theta0+Psi) + exp(-kappa)*(Phi'+Psi')  [SW + ISW]
  S_t1 = g*v_b = g*theta_b_N/k                      [Doppler in j_l']
  Delta_l = int [S_t0 * j_l + S_t1 * j_l'] dtau

EFFICIENT form (this file -- IBP the Doppler too):
  S_t0_eff = g*(Theta0+Psi) + exp(-kappa)*(Phi'+Psi')  [SW + ISW, unchanged]
           + (g*theta_b_N' + g'*theta_b_N) / k^2       [Doppler IBP'd into j_l]
  S_t1_eff = 0                                          [nothing in j_l'!]
  Delta_l = int S_t0_eff * j_l dtau

IMPORTANT NOTE (2026-03-22):
  The Doppler IBP form requires a MUCH denser tau grid than the simple form
  because g'(tau) (the visibility derivative) is a sharp, sign-changing function
  around the visibility peak. CLASS uses O(1000-2000) adaptively-spaced tau points
  to resolve this. With our current 200-point grid, the IBP form gives WORSE
  results (9-12% RMS) than the simple form (6.57% RMS) because the g'*theta_b_N
  term is not well-resolved.

  To make this work at <1% RMS, the tau grid must be densified to ~500-1000 points
  in the visibility region, or the source function must be evaluated on a much
  finer grid internally (e.g., by re-interpolating the ODE solutions).

All quantities are in NEWTONIAN gauge (gauge-invariant combinations).

Usage:
  python -m mlx_class.solver_sync_efficient [--N_k 500] [--parallel]

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
from .perturbations_sync import (
    make_sync_rhs, adiabatic_ic_sync, gauge_transform, diagnose_metric,
    n_var_sync, idx_fn_start, idx_fg,
    IDX_ETA, IDX_DELTA_C, IDX_DELTA_B, IDX_THETA_B, IDX_FG_START,
    L_GAMMA_MAX, L_NU_MAX_SYNC,
)
from .solver_sync import build_snapshot_grid


# ============================================================================
# Worker: solve one k-mode and return quantities needed for efficient LOS
# ============================================================================

def _solve_single_k_efficient_worker(args):
    """
    Worker for multiprocessing. Solves one k-mode and returns quantities
    needed for the CLASS-style efficient source functions.

    Returns numpy array of shape (N_snap, 9) with columns:
      [Phi_N, Psi_N, Theta0_N, theta_b_N, theta_b_N_prime,
       Phi_N_prime, Theta2, alpha, alpha_prime]
    or None if the ODE solve fails.

    Key difference from solver_sync._solve_single_k_worker:
    - Returns Newtonian gauge theta_b (not v_b = theta_b/k)
    - Returns theta_b_N' (time derivative of Newtonian baryon velocity)
    - Returns Phi_N' separately (not Phi_N' + Psi_N')
    - Returns alpha, alpha_prime for potential use
    """
    import numpy as np
    from scipy.integrate import solve_ivp

    (k, bg_arrays, tau_end, lg_max, ln_max, method, rtol, atol,
     tau_all, a_snap, calH_snap) = args

    # Reconstruct lightweight background
    class _BGLite:
        __slots__ = ('tau_grid', 'a_grid', 'calH_grid', 'R_grid', 'kappa_dot_grid',
                     'tau_rec', 'tau_0')
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
        make_sync_rhs, adiabatic_ic_sync, gauge_transform, diagnose_metric,
        IDX_ETA, IDX_DELTA_C, IDX_DELTA_B, IDX_THETA_B, IDX_FG_START,
    )
    from mlx_class.background import H0_Mpc as _H0_val
    from mlx_class.perturbations_neutrino import _OMEGA_GAMMA as Og, _OMEGA_NU as On

    H02 = _H0_val ** 2
    k2 = k * k

    # Build RHS and initial conditions
    tau_init = bg_lite.tau_grid[1]
    rhs_fn, nvar = make_sync_rhs(k, bg_lite, lg_max, ln_max)
    y0 = adiabatic_ic_sync(k, tau_init, lg_max, ln_max)

    # Solve ODE
    sol = solve_ivp(rhs_fn, [tau_init, tau_end], y0,
                    method=method, dense_output=True,
                    rtol=rtol, atol=atol)

    if not sol.success:
        return None

    # Precompute background interpolation arrays
    _tau = bg_lite.tau_grid
    _calH = bg_lite.calH_grid
    _a = bg_lite.a_grid
    _R = bg_lite.R_grid
    _kd = bg_lite.kappa_dot_grid  # negative

    # Evaluate at snapshot points
    N_snap = len(tau_all)
    # Columns: Phi_N, Psi_N, Theta0_N, theta_b_N, theta_b_N_prime,
    #          Phi_N_prime, Theta2, alpha_val, alpha_prime_val
    results = np.zeros((N_snap, 9), dtype=np.float64)

    fn_s = lg_max + 1 + IDX_FG_START  # start of neutrino hierarchy

    for it in range(N_snap):
        tau = tau_all[it]
        y = sol.sol(tau)
        a = a_snap[it]
        calH = calH_snap[it]

        # Background quantities at this tau
        R = float(np.interp(tau, _tau, _R))
        kappa_dot = float(np.interp(tau, _tau, _kd))  # negative
        abs_kd = abs(kappa_dot)

        # Diagnose metric perturbations
        h_prime, eta_prime = diagnose_metric(y, k, calH, a, lg_max, ln_max)

        # Gauge transformation
        gt = gauge_transform(y, k, calH, a, lg_max, ln_max, h_prime, eta_prime)

        eta_val = y[IDX_ETA]
        theta_b_S = y[IDX_THETA_B]  # sync gauge
        Fg = y[IDX_FG_START: IDX_FG_START + lg_max + 1]
        Fn = y[fn_s: fn_s + ln_max + 1]

        # Newtonian gauge potentials
        Phi_N = gt['Phi_N']
        Psi_N = gt['Psi_N']

        # Use alpha from gauge_transform (already capped for superhorizon modes)
        alpha_val = gt['alpha']

        # alpha' = -2*calH*alpha + eta - (9/2)*(a^2/k^2)*pi_tot
        # where pi_tot = (rho+p)*sigma for all species
        sigma_g = 0.5 * Fg[2] if lg_max >= 2 else 0.0
        sigma_n = 0.5 * Fn[2] if ln_max >= 2 else 0.0
        # pi_tot in CLASS units: (rho+p)*sigma -> H0^2 * [...]
        # The anisotropic stress term: 12*H0^2/(a^2*k^2) * (Og*sigma_g + On*sigma_n) = Phi_N - Psi_N
        # So: (9/2)*(a^2/k^2)*pi_tot = (3/4)*(Phi_N - Psi_N) ... no, let me be more careful.
        #
        # From CLASS (line 6583-6586):
        # alpha' = -2*calH*alpha + eta - (9/2)*(a^2/k^2)*(rho+p)_shear
        # In our units: (rho+p)_shear = H0^2 * [(4/3)*Og/a^4 * sigma_g + (4/3)*On/a^4 * sigma_n]
        # So: (9/2)*(a^2/k^2)*(rho+p)_shear = (9/2)*(a^2/k^2)*H0^2*[(4/3)*Og/a^4*sigma_g + (4/3)*On/a^4*sigma_n]
        #   = 6*H0^2/(k^2*a^2) * [Og*sigma_g + On*sigma_n]
        # This is exactly = (1/2)*(Phi_N - Psi_N) since
        #   Psi_N = Phi_N - 12*H0^2/(a^2*k^2)*(Og*sigma_g + On*sigma_n)
        #   => Phi_N - Psi_N = 12*H0^2/(a^2*k^2)*(Og*sigma_g + On*sigma_n)
        # So: (9/2)*(a^2/k^2)*pi = 6*H0^2/(k^2*a^2)*(Og*sigma_g+On*sigma_n) = (1/2)*(Phi_N-Psi_N)
        aniso_stress_term = 6.0 * H02 / (k2 * a * a) * (Og * sigma_g + On * sigma_n)
        alpha_prime_val = -2.0 * calH * alpha_val + eta_val - aniso_stress_term

        # Newtonian gauge baryon velocity: theta_b_N = theta_b_S + k^2 * alpha
        theta_b_N = theta_b_S + k2 * alpha_val

        # Sync gauge theta_b' from EOM:
        # theta_b_S' = -calH * theta_b_S + |kappa_dot| / R * (theta_gamma_S - theta_b_S)
        theta_gamma_S = 0.75 * k * Fg[1]
        theta_b_S_prime = -calH * theta_b_S + abs_kd / R * (theta_gamma_S - theta_b_S)

        # Newtonian gauge theta_b': theta_b_N' = theta_b_S' + k^2 * alpha'
        theta_b_N_prime = theta_b_S_prime + k2 * alpha_prime_val

        # Gauge-invariant temperature monopole:
        # Theta0_N = F_g0/4 - eta + Phi_N (avoids large calH*alpha subtraction)
        Theta0_N = Fg[0] / 4.0 - eta_val + Phi_N

        # Photon quadrupole (gauge-invariant for l>=2)
        Theta2 = Fg[2] / 4.0 if lg_max >= 2 else 0.0

        results[it, 0] = Phi_N
        results[it, 1] = Psi_N
        results[it, 2] = Theta0_N
        results[it, 3] = theta_b_N
        results[it, 4] = theta_b_N_prime
        # Phi_N_prime: compute later via cubic spline over tau
        results[it, 5] = 0.0  # placeholder
        results[it, 6] = Theta2
        results[it, 7] = alpha_val
        results[it, 8] = alpha_prime_val

    return results


# ============================================================================
# Compute visibility derivative g'(tau)
# ============================================================================

def compute_g_prime(bg, tau_all):
    """
    Compute g'(tau) = dg/dtau analytically.

    g(tau) = |kappa_dot| * exp(-kappa)

    g'(tau) = d/dtau[|kappa_dot| * exp(-kappa)]
            = |kappa_dot|' * exp(-kappa) + |kappa_dot| * (-kappa_dot) * exp(-kappa)
            = (|kappa_dot|' + |kappa_dot|^2) * exp(-kappa)

    where kappa_dot < 0, so |kappa_dot| = -kappa_dot, and
    d(|kappa_dot|)/dtau = -kappa_ddot.

    In practice, we use a cubic spline on the precomputed visibility to get
    a smooth derivative, which is more robust than the analytical formula
    (which requires kappa_ddot that we don't store).
    """
    # Use cubic spline on the full background grid for smoothness
    cs = CubicSpline(bg.tau_grid, bg.visibility_grid)
    g_prime = cs(tau_all, 1)  # first derivative
    return g_prime


# ============================================================================
# Full pipeline (parallel version with efficient source functions)
# ============================================================================

def run_sync_efficient(N_k=500, k_min=3e-4, k_max=0.35, method='Radau',
                       lg_max=L_GAMMA_MAX, ln_max=L_NU_MAX_SYNC,
                       rtol=1e-6, atol=1e-9, n_workers=None,
                       verbose=True):
    """
    Parallel synchronous gauge pipeline with CLASS-style efficient source functions.

    The ODE solving is identical to solver_sync.py. Only the LOS integration
    changes: instead of the simple (SW + Doppler*j_l' + ISW) form, we use
    the efficient form where Doppler is moved into S_t0 via integration by parts.

    Parameters
    ----------
    N_k : int
        Number of k-modes.
    n_workers : int or None
        Number of worker processes. Defaults to os.cpu_count().
    """
    t_total = time.time()

    if n_workers is None:
        n_workers = os.cpu_count() or 4

    # ================================================================
    # Step 1: Background
    # ================================================================
    if verbose:
        print("=" * 65)
        print("Synchronous Gauge Boltzmann Solver (EFFICIENT SOURCE FUNCTIONS)")
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

    # The efficient source functions need a DENSER tau grid than the simple form
    # because they involve g'(tau) (visibility derivative) which is a sharp,
    # sign-changing function around the visibility peak. The simple form
    # uses j_l' which is smooth and requires fewer tau points.
    # Use 200 points in the visibility region (vs 60 in solver_sync.py).
    tau_vis, tau_late, tau_all = build_snapshot_grid(
        bg, N_vis=200, N_early_isw=40, N_late_isw=30)
    N_snap = len(tau_all)

    if verbose:
        print(f"\n--- Step 2: Perturbations ({N_k} k-modes, {method}, "
              f"{n_workers} workers) ---")
        print(f"k range: [{k_min:.1e}, {k_max:.2f}] Mpc^-1")
        print(f"Snapshot grid: {N_snap} tau points")
        print(f"  Visibility: {len(tau_vis)} pts, tau=[{tau_vis[0]:.1f}, {tau_vis[-1]:.1f}]")
        print(f"  Late ISW:   {len(tau_late)} pts, tau=[{tau_late[0]:.1f}, {tau_late[-1]:.1f}]")
        if hasattr(bg, 'z_reio') and bg.tau_reio > 0:
            print(f"  Reionization: z_reio={bg.z_reio:.2f}, "
                  f"tau_reio={bg.tau_reio:.4f}")
        print(f"Hierarchy: lg_max={lg_max}, ln_max={ln_max}, "
              f"nvar={n_var_sync(lg_max, ln_max)}")
        sys.stdout.flush()

    # ================================================================
    # Step 3: Solve all k-modes in parallel
    # ================================================================
    t0 = time.time()

    # Precompute background at snapshots
    a_snap = np.asarray(bg.a_at_tau(tau_all), dtype=np.float64)
    calH_snap = np.asarray(bg.calH_at_tau(tau_all), dtype=np.float64)

    # Serialize background as raw numpy arrays (picklable)
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

    # Build work items
    work_items = [
        (k, bg_arrays, tau_end, lg_max, ln_max, method, rtol, atol,
         tau_all, a_snap, calH_snap)
        for k in k_arr
    ]

    if verbose:
        print(f"Dispatching {N_k} k-modes to {n_workers} workers...")
        sys.stdout.flush()

    with Pool(processes=n_workers) as pool:
        worker_results = pool.map(_solve_single_k_efficient_worker, work_items)

    t_pert = time.time() - t0

    # ================================================================
    # Collect results
    # ================================================================
    Phi_N = np.zeros((N_k, N_snap))
    Psi_N = np.zeros((N_k, N_snap))
    Theta0_N = np.zeros((N_k, N_snap))
    theta_b_N = np.zeros((N_k, N_snap))       # Newtonian gauge theta_b (NOT v_b)
    theta_b_N_prime = np.zeros((N_k, N_snap))  # d(theta_b_N)/dtau
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
        theta_b_N[ik] = res[:, 3]
        theta_b_N_prime[ik] = res[:, 4]
        # Phi_N_prime placeholder (column 5) -- will compute via spline
        Theta2_arr[ik] = res[:, 6]

    if verbose:
        print(f"Perturbations: {t_pert:.1f}s ({n_failed} failed, "
              f"{N_k - n_failed} OK)")

    # ================================================================
    # Step 4: Compute derivatives by cubic spline
    # ================================================================
    # Note: We use the SAME ISW formulation as solver_sync.py (which already
    # did one IBP on exp(-kappa)*k*psi), but replace the Doppler j_l' integral
    # with the IBP'd Doppler j_l form.

    # (Phi_N + Psi_N)' for ISW
    PhiPsi = Phi_N + Psi_N
    PhiPsi_prime = np.zeros_like(PhiPsi)
    for ik in range(N_k):
        cs = CubicSpline(tau_all, PhiPsi[ik])
        PhiPsi_prime[ik] = cs(tau_all, 1)

    # theta_b_N' via cubic spline (more robust than the analytical formula
    # which involves alpha_prime and is sensitive to the alpha capping)
    theta_b_N_prime_spline = np.zeros_like(theta_b_N)
    for ik in range(N_k):
        cs = CubicSpline(tau_all, theta_b_N[ik])
        theta_b_N_prime_spline[ik] = cs(tau_all, 1)
    # Use spline derivative instead of analytical
    theta_b_N_prime = theta_b_N_prime_spline

    # ================================================================
    # Step 5: LOS integration with EFFICIENT source functions
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
        print(f"\n--- Step 5: LOS C_l with EFFICIENT sources ({N_ell} ells) ---")
        sys.stdout.flush()

    # Primordial spectrum
    P_R = A_s * (k_arr / k_pivot) ** (n_s - 1.0)

    # Background at snapshots
    g_snap = bg.visibility_at_tau(tau_all)
    kappa_snap = bg.kappa_at_tau(tau_all)
    exp_neg_kappa = np.exp(-kappa_snap)
    chi_snap = bg.tau_0 - tau_all
    dtau_all = np.diff(tau_all)

    # g'(tau) -- visibility derivative
    g_prime_snap = compute_g_prime(bg, tau_all)

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

    # Polarization
    alpha_P = 1.7
    pol_prefactor = 0.75 * alpha_P

    eps_prefactors = np.zeros(N_ell)
    for il, ell in enumerate(ell_values):
        l = int(ell)
        if l >= 2:
            eps_prefactors[il] = np.sqrt(float((l - 1) * l * (l + 1) * (l + 2)))

    if verbose:
        print(f"Computing LOS transfer functions with efficient sources (TT + EE)...")
        sys.stdout.flush()

    # ----------------------------------------------------------------
    # EFFICIENT LOS FORMULA:
    #
    # Our old solver_sync.py already did one IBP (moving exp(-kappa)*k*Psi
    # from S_t1 into S_t0), giving:
    #   S_t0 = g*(Theta0+Psi) + exp(-kappa)*(Phi'+Psi')  [j_l]
    #   S_t1 = g*v_b = g*theta_b_N/k                     [j_l']
    #
    # The remaining cancellation is in S_t1: the Doppler j_l' integral
    # oscillates and loses precision due to the finite tau grid.
    #
    # We now IBP the Doppler too:
    #   int g*theta_b_N/k * j_l'(x) dtau
    #   = (1/k^2) int [g'*theta_b_N + g*theta_b_N'] * j_l(x) dtau
    #
    # This gives:
    #   S_t0_eff = g*(Theta0+Psi) + exp(-kappa)*(Phi'+Psi')
    #            + (g*theta_b_N' + g'*theta_b_N)/k^2         [Doppler IBP]
    #   S_t1_eff = 0  (no j_l' integral needed!)
    #
    # All source terms are now multiplied by j_l only, eliminating
    # the problematic j_l' oscillatory integral entirely.
    # ----------------------------------------------------------------

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
            w = 0.5 * (dtau_all[max(0, it-1)] + dtau_all[min(it, len(dtau_all)-1)])

        # --- SW: g * (Theta0_N + Psi_N) ---
        S_SW = g_snap[it] * (Theta0_N[:, it] + Psi_N[:, it])  # (N_k,)

        # --- ISW: exp(-kappa) * (Phi_N' + Psi_N') ---
        # (Same as solver_sync.py -- already has the psi IBP baked in)
        S_ISW = exp_neg_kappa[it] * PhiPsi_prime[:, it]  # (N_k,)

        # --- Doppler IBP: (g * theta_b_N' + g' * theta_b_N) / k^2 ---
        # This replaces the old g*v_b * j_l' with a j_l integral,
        # eliminating the oscillatory j_l' integration entirely.
        S_Dop_IBP = (g_snap[it] * theta_b_N_prime[:, it]
                     + g_prime_snap[it] * theta_b_N[:, it]) / (k_arr ** 2)  # (N_k,)

        # --- Total efficient source (all multiplied by j_l) ---
        S_t0_eff = S_SW + S_ISW + S_Dop_IBP  # (N_k,)

        # --- EE source (unchanged): g * pol_prefactor * Theta_2 ---
        S_E = g_snap[it] * pol_prefactor * Theta2_arr[:, it]  # (N_k,)

        if use_gpu_bessel:
            x_arr = k_arr * chi
            jl_mx, jlp_mx = bessel_table.eval_jl_jlp(x_arr.astype(np.float32))
            mx.eval(jl_mx, jlp_mx)
            jl = np.array(jl_mx)    # (N_ell, N_k)
            jlp = np.array(jlp_mx)  # (N_ell, N_k)

            x_safe = np.where(np.abs(x_arr) < 1e-10, 1e-10, x_arr)
            inv_x2 = 1.0 / (x_safe ** 2)

            for il in range(N_ell):
                # TT: Delta_l += w * S_t0_eff * j_l  (NO j_l' term!)
                Delta_l[il, :] += w * S_t0_eff * jl[il, :]
                # EE: same as before
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
                Delta_l[il, :] += w * S_t0_eff * jl
                if eps_prefactors[il] > 0:
                    eps_l = eps_prefactors[il] * jl * inv_x2
                    eps_l = np.where(np.abs(x) < 1e-10, 0.0, eps_l)
                    Delta_l_E[il, :] += w * S_E * eps_l

    t_los = time.time() - t0
    if verbose:
        print(f"LOS integration (efficient): {t_los:.2f}s")

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
        description='Sync gauge Boltzmann solver with CLASS efficient source functions')
    parser.add_argument('--N_k', type=int, default=500)
    parser.add_argument('--method', type=str, default='Radau')
    parser.add_argument('--rtol', type=float, default=1e-6)
    parser.add_argument('--atol', type=float, default=1e-9)
    parser.add_argument('--lg_max', type=int, default=L_GAMMA_MAX)
    parser.add_argument('--ln_max', type=int, default=L_NU_MAX_SYNC)
    parser.add_argument('--n_workers', type=int, default=None,
                        help='Number of worker processes')
    args = parser.parse_args()

    result = run_sync_efficient(
        N_k=args.N_k, method=args.method,
        lg_max=args.lg_max, ln_max=args.ln_max,
        rtol=args.rtol, atol=args.atol,
        n_workers=args.n_workers,
        verbose=True)

    ell = result['ell']
    Dl_TT = result['Dl_TT']
    Dl_EE = result['Dl_EE']
    Dl_TE = result['Dl_TE']

    print("\n--- D_l values at key multipoles (EFFICIENT) ---")
    print(f"{'l':>6s}  {'D_l^TT':>12s}  {'D_l^EE':>12s}  {'D_l^TE':>12s}")
    print("-" * 50)
    for target_ell in [2, 10, 50, 100, 220, 537, 800, 1000, 1500, 2000]:
        idx = np.argmin(np.abs(ell - target_ell))
        print(f"  {ell[idx]:5d}  {Dl_TT[idx]:12.2f}  {Dl_EE[idx]:12.4f}  {Dl_TE[idx]:12.2f}")

    # Compare with CLASS reference if available
    CLASS_FILE_TT = '/Users/akaihuangm1/Desktop/github/gdm_class_public/output/mlx_ref_lcdm_cl.dat'
    CLASS_FILE_POL = '/Users/akaihuangm1/Desktop/github/gdm_class_public/output/lcdm_pol_ref_cl.dat'
    try:
        data = np.loadtxt(CLASS_FILE_TT)
        class_ell = data[:, 0].astype(int)
        class_Dl = data[:, 1] * (T_CMB * 1e6) ** 2

        # RMS error
        from scipy.interpolate import interp1d as interp1d_local
        mask = (class_ell >= 2) & (class_ell <= 2500)
        ell_use = class_ell[mask]
        Dl_ref = class_Dl[mask]
        f_test = interp1d_local(ell, Dl_TT, kind='cubic', fill_value='extrapolate')
        Dl_test_interp = f_test(ell_use)
        Dl_max = np.max(np.abs(Dl_ref))
        rel_err = (Dl_test_interp - Dl_ref) / Dl_max
        rms = np.sqrt(np.mean(rel_err ** 2))

        print(f"\n--- Comparison with CLASS ---")
        print(f"RMS error (l=2-2500): {rms*100:.2f}%")

        # Second peak test
        idx_537_class = np.argmin(np.abs(class_ell - 537))
        idx_537_ours = np.argmin(np.abs(ell - 537))
        ratio_537 = Dl_TT[idx_537_ours] / class_Dl[idx_537_class]
        print(f"D_l(537) ratio (ours/CLASS): {ratio_537:.3f}")
        print(f"  Ours:  {Dl_TT[idx_537_ours]:.2f} muK^2")
        print(f"  CLASS: {class_Dl[idx_537_class]:.2f} muK^2")

        # First peak
        idx_220_class = np.argmin(np.abs(class_ell - 220))
        idx_220_ours = np.argmin(np.abs(ell - 220))
        ratio_220 = Dl_TT[idx_220_ours] / class_Dl[idx_220_class]
        print(f"D_l(220) ratio (ours/CLASS): {ratio_220:.3f}")

        # Third peak
        idx_800_class = np.argmin(np.abs(class_ell - 800))
        idx_800_ours = np.argmin(np.abs(ell - 800))
        ratio_800 = Dl_TT[idx_800_ours] / class_Dl[idx_800_class]
        print(f"D_l(800) ratio (ours/CLASS): {ratio_800:.3f}")

    except Exception as e:
        print(f"\n[INFO] Could not compare with CLASS reference: {e}")


if __name__ == '__main__':
    main()
