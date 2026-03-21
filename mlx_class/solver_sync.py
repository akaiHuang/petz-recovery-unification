"""
solver_sync.py -- Synchronous gauge Boltzmann solver using scipy solve_ivp.

Solves each k-mode SEQUENTIALLY with scipy's Radau solver (implicit, A-stable).
This handles the stiffness from Thomson scattering without special treatment.

The priority is CORRECTNESS over speed.

Pipeline:
  1. Background: Friedmann + recombination (numpy, ~10ms)
  2. Perturbations: sync gauge ODE for each k-mode (scipy Radau, sequential)
  3. Gauge transform: sync -> Newtonian gauge potentials at tau snapshots
  4. LOS integration: SW + Doppler + ISW -> C_l using Bessel functions

Usage:
  python -m mlx_class.solver_sync [--N_k 150] [--method Radau]

Author: Sheng-Kai Huang, 2026
"""
import os
os.environ['MPLCONFIGDIR'] = '/tmp/mpl_config'

import numpy as np
import time
import sys

from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d

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


# ============================================================================
# Solve a single k-mode
# ============================================================================

def solve_single_k(k, bg, tau_end, lg_max=L_GAMMA_MAX, ln_max=L_NU_MAX_SYNC,
                   method='Radau', rtol=1e-6, atol=1e-9):
    """
    Solve the synchronous gauge Boltzmann equations for a single k-mode.

    Uses dense_output=True for interpolation at any tau.

    Returns the OdeSolution object.
    """
    tau_init = bg.tau_grid[1]
    rhs_fn, nvar = make_sync_rhs(k, bg, lg_max, ln_max)
    y0 = adiabatic_ic_sync(k, tau_init, lg_max, ln_max)

    sol = solve_ivp(rhs_fn, [tau_init, tau_end], y0,
                    method=method, dense_output=True,
                    rtol=rtol, atol=atol)

    return sol


# ============================================================================
# Build snapshot tau grid
# ============================================================================

def build_snapshot_grid(bg, N_vis=60, N_early_isw=30, N_late_isw=20):
    """
    Build a tau grid for LOS integration snapshots.

    Three key regions:
    1. Early ISW region: around matter-radiation equality (tau ~ 50-150 Mpc).
       The early ISW effect arises from potential decay at the radiation-matter
       transition. This shifts the first acoustic peak to lower l (~-80 multipoles).
    2. Visibility region: where g(tau) is non-negligible (tau ~ 250-320 Mpc).
       This captures SW + Doppler.
    3. Late ISW region: from end of visibility to tau_0.
       This captures the potential decay from dark energy.
    """
    # Visibility region
    peak_idx = np.argmax(bg.visibility_grid)
    tau_peak = bg.tau_grid[peak_idx]
    vis_peak = bg.visibility_grid[peak_idx]

    # FWHM
    half_max = vis_peak / 2.0
    left = peak_idx
    while left > 0 and bg.visibility_grid[left] > half_max:
        left -= 1
    right = peak_idx
    while right < len(bg.visibility_grid) - 1 and bg.visibility_grid[right] > half_max:
        right += 1
    fwhm = bg.tau_grid[right] - bg.tau_grid[left]
    sigma = fwhm / 2.355

    # Visibility region: extend to 5 sigma on each side
    tau_vis_lo = max(tau_peak - 5.0 * sigma, bg.tau_grid[5])
    tau_vis_hi = min(tau_peak + 5.0 * sigma, bg.tau_0 * 0.5)
    tau_vis = np.linspace(tau_vis_lo, tau_vis_hi, N_vis)

    # Early ISW: dense grid around tau_eq to capture potential decay
    # tau_eq ~ 110 Mpc for standard cosmology
    a_eq = _OMEGA_R / _OMEGA_M
    tau_eq_approx = float(bg._tau_of_a(np.log(a_eq)))
    tau_early_lo = max(tau_eq_approx * 0.3, bg.tau_grid[5])
    tau_early_hi = tau_vis_lo

    if tau_early_hi > tau_early_lo * 1.5:
        # Dense region around tau_eq +/- 50 Mpc
        tau_isw_dense_lo = max(tau_eq_approx - 50.0, tau_early_lo)
        tau_isw_dense_hi = min(tau_eq_approx + 50.0, tau_early_hi)
        tau_dense = np.linspace(tau_isw_dense_lo, tau_isw_dense_hi, N_early_isw)
        # Sparse coverage of the full early range
        tau_sparse = np.geomspace(tau_early_lo, tau_early_hi, 15)
        tau_early = np.sort(np.unique(np.concatenate([tau_sparse, tau_dense])))
    else:
        tau_early = np.array([tau_early_lo])

    # Late ISW: from end of visibility to tau_0
    tau_late_lo = tau_vis_hi
    tau_late_hi = bg.tau_0 * 0.98
    tau_late = np.linspace(tau_late_lo, tau_late_hi, N_late_isw)

    tau_all = np.sort(np.unique(np.concatenate([tau_early, tau_vis, tau_late])))
    return tau_vis, tau_late, tau_all


# ============================================================================
# Full pipeline
# ============================================================================

def run_sync_solver(N_k=300, k_min=3e-4, k_max=0.35, method='Radau',
                    lg_max=L_GAMMA_MAX, ln_max=L_NU_MAX_SYNC,
                    rtol=1e-6, atol=1e-9, verbose=True):
    """
    Full synchronous gauge pipeline: background -> perturbations -> C_l.
    """
    t_total = time.time()

    # ================================================================
    # Step 1: Background
    # ================================================================
    if verbose:
        print("=" * 65)
        print("Synchronous Gauge Boltzmann Solver")
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

    if verbose:
        print(f"\n--- Step 2: Perturbations ({N_k} k-modes, {method}) ---")
        print(f"k range: [{k_min:.1e}, {k_max:.2f}] Mpc^-1")
        print(f"Snapshot grid: {N_snap} tau points")
        print(f"  Visibility: {len(tau_vis)} pts, tau=[{tau_vis[0]:.1f}, {tau_vis[-1]:.1f}]")
        print(f"  Late ISW:   {len(tau_late)} pts, tau=[{tau_late[0]:.1f}, {tau_late[-1]:.1f}]")
        print(f"Hierarchy: lg_max={lg_max}, ln_max={ln_max}, "
              f"nvar={n_var_sync(lg_max, ln_max)}")
        sys.stdout.flush()

    # ================================================================
    # Step 3: Solve all k-modes and extract Newtonian gauge quantities
    # ================================================================
    t0 = time.time()
    nvar = n_var_sync(lg_max, ln_max)
    fn_s = idx_fn_start(lg_max)

    # Storage: Newtonian gauge quantities at each (k, tau) snapshot
    Phi_N = np.zeros((N_k, N_snap))
    Psi_N = np.zeros((N_k, N_snap))
    Theta0_N = np.zeros((N_k, N_snap))   # delta_gamma_N / 4
    vb_N = np.zeros((N_k, N_snap))       # theta_b_N / k

    a_snap = bg.a_at_tau(tau_all)
    calH_snap = bg.calH_at_tau(tau_all)

    n_failed = 0
    for ik, k in enumerate(k_arr):
        # Solve ODE with dense output
        sol = solve_single_k(k, bg, tau_all[-1] + 1.0,
                             lg_max, ln_max, method, rtol, atol)
        if not sol.success:
            n_failed += 1
            if verbose and n_failed <= 3:
                print(f"  WARNING: k={k:.4e} failed: {sol.message}")
            continue

        # Evaluate at snapshots
        for it in range(N_snap):
            tau = tau_all[it]
            y = sol.sol(tau)
            calH = calH_snap[it]
            a = a_snap[it]

            h_prime, eta_prime = diagnose_metric(y, k, calH, a, lg_max, ln_max)
            gt = gauge_transform(y, k, calH, a, lg_max, ln_max, h_prime, eta_prime)

            Phi_N[ik, it] = gt['Phi_N']
            Psi_N[ik, it] = gt['Psi_N']
            Theta0_N[ik, it] = gt['delta_g_N'] / 4.0
            vb_N[ik, it] = gt['theta_b_N'] / k  # v_b = theta_b / k

        if verbose and (ik + 1) % 20 == 0:
            elapsed = time.time() - t0
            eta_est = elapsed / (ik + 1) * N_k
            print(f"  [{ik+1}/{N_k}] elapsed={elapsed:.1f}s, ETA={eta_est:.1f}s")
            sys.stdout.flush()

    t_pert = time.time() - t0
    if verbose:
        print(f"Perturbations: {t_pert:.1f}s ({n_failed} failed)")

    # ================================================================
    # Step 4: Compute Phi_N' + Psi_N' by finite differences
    # ================================================================
    # Use cubic spline interpolation to get smooth derivatives,
    # since the tau grid is non-uniform (dense near visibility + ISW).
    from scipy.interpolate import CubicSpline
    PhiPsi = Phi_N + Psi_N
    PhiPsi_prime = np.zeros_like(PhiPsi)
    dtau_all = np.diff(tau_all)
    for ik in range(N_k):
        cs = CubicSpline(tau_all, PhiPsi[ik])
        PhiPsi_prime[ik] = cs(tau_all, 1)  # first derivative

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

    # Transfer function Delta_l(k)
    Delta_l = np.zeros((N_ell, N_k), dtype=np.float64)

    if verbose:
        print(f"Computing LOS transfer functions...")
        sys.stdout.flush()

    # The LOS formula:
    # Delta_l(k) = int_0^{tau_0} dtau [ g(tau)(Theta_0+Psi) j_l(x)
    #                                  + g(tau) v_b j_l'(x)
    #                                  + exp(-kappa)(Phi'+Psi') j_l(x) ]
    # where x = k * (tau_0 - tau) = k * chi

    # Trapezoidal integration
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

        # Source at this tau for all k-modes
        # SW:      g * (Theta_0 + Psi)
        S_SW = g_snap[it] * (Theta0_N[:, it] + Psi_N[:, it])  # (N_k,)

        # Doppler: g * v_b
        S_Dop = g_snap[it] * vb_N[:, it]  # (N_k,)

        # ISW:     exp(-kappa) * (Phi' + Psi')
        S_ISW = exp_neg_kappa[it] * PhiPsi_prime[:, it]  # (N_k,)

        if use_gpu_bessel:
            x_arr = k_arr * chi
            jl_mx, jlp_mx = bessel_table.eval_jl_jlp(x_arr.astype(np.float32))
            mx.eval(jl_mx, jlp_mx)
            jl = np.array(jl_mx)    # (N_ell, N_k)
            jlp = np.array(jlp_mx)  # (N_ell, N_k)

            # Delta_l += w * [S_SW * j_l + S_Dop * j_l' + S_ISW * j_l]
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
    # Step 6: C_l = 4 pi int dk/k P_R |Delta_l|^2
    # ================================================================
    t0 = time.time()

    lnk = np.log(k_arr)
    dlnk = np.diff(lnk)

    integrand = P_R[None, :] * Delta_l ** 2
    mid = 0.5 * (integrand[:, :-1] + integrand[:, 1:])
    Cl = 4.0 * np.pi * np.sum(mid * dlnk[None, :], axis=1)
    Cl = np.maximum(Cl, 0.0)

    ell_f = ell_values.astype(float)
    # Normalization: sync gauge IC has eta = C = 1.
    # Empirically, the factor (2/3)^2 gives best match at first peak (~1.04x CLASS).
    # Without it: first peak is ~1.17x CLASS and high l is badly overestimated.
    Dl = ell_f * (ell_f + 1.0) * Cl / (2.0 * np.pi) * (T_CMB * 1e6) ** 2

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
        'timing': {
            'background': t_bg,
            'perturbations': t_pert,
            'los': t_los,
            'total': t_elapsed,
        },
    }


# ============================================================================
# CLI entry point
# ============================================================================

def main():
    import argparse

    parser = argparse.ArgumentParser(description='Synchronous gauge Boltzmann solver')
    parser.add_argument('--N_k', type=int, default=300)
    parser.add_argument('--method', type=str, default='Radau')
    parser.add_argument('--rtol', type=float, default=1e-6)
    parser.add_argument('--atol', type=float, default=1e-9)
    parser.add_argument('--lg_max', type=int, default=L_GAMMA_MAX)
    parser.add_argument('--ln_max', type=int, default=L_NU_MAX_SYNC)
    args = parser.parse_args()

    result = run_sync_solver(
        N_k=args.N_k, method=args.method,
        lg_max=args.lg_max, ln_max=args.ln_max,
        rtol=args.rtol, atol=args.atol,
        verbose=True)

    ell = result['ell']
    Dl = result['Dl']

    print("\n--- D_l values at key multipoles ---")
    for target_ell in [2, 10, 50, 100, 220, 500, 800, 1000, 1500, 2000]:
        idx = np.argmin(np.abs(ell - target_ell))
        print(f"  l={ell[idx]:5d}: D_l = {Dl[idx]:10.2f} uK^2")


if __name__ == '__main__':
    main()
