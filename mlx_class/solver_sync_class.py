"""
solver_sync_class.py -- CLASS-style three-phase Boltzmann solver.

Implements the full CLASS approximation strategy:
  Phase 1 (TCA):  Tight Coupling Approximation with analytical shear/slip
  Phase 2 (Full): Full Boltzmann hierarchy (photon + polarization + neutrino)
  Phase 3 (RSA):  Radiation Streaming Approximation (algebraic photons/neutrinos)

Key improvements over solver_sync.py:
  1. TCA phase: smaller ODE system (~22 vars vs ~56), faster early integration
  2. Analytical shear/slip: correct sigma_g to first order in tau_c
  3. Proper seeding: F_g2, F_g3, E_0..E_3 seeded from TCA at the switch
  4. RSA phase: only 4 variables after k*tau > 45, huge speed gain for late times
  5. Baryon sound speed c_s^2 k^2 delta_b included in all phases

Target: ~0.3% RMS accuracy vs CLASS (improvement from ~3%).

Usage:
  python -m mlx_class.solver_sync_class [--N_k 500] [--parallel]

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
from .perturbations_sync_class import (
    solve_three_phase,
    L_GAMMA_MAX, L_POL_MAX, L_NU_MAX,
    n_var_full,
)
from .solver_sync import build_snapshot_grid


# ============================================================================
# Multiprocessing worker
# ============================================================================

def _solve_worker_class(args):
    """
    Worker function for multiprocessing. Solves one k-mode with three-phase.

    Returns a numpy array of shape (N_snap, 6) with columns:
      [Phi_N, Psi_N, Theta0_N, vb_N, Theta2, Pi_over_4]
    Returns None if the solve fails.
    """
    import numpy as np

    (k, bg_arrays, tau_end, lg_max, lp_max, ln_max, method, rtol, atol,
     tau_all, a_snap, calH_snap) = args

    # Reconstruct lightweight background
    class _BGLite:
        __slots__ = ('tau_grid', 'a_grid', 'calH_grid', 'R_grid',
                     'kappa_dot_grid', 'kappa_grid', 'tau_rec', 'tau_0')
        def __init__(self, arrays):
            self.tau_grid = arrays['tau_grid']
            self.a_grid = arrays['a_grid']
            self.calH_grid = arrays['calH_grid']
            self.R_grid = arrays['R_grid']
            self.kappa_dot_grid = arrays['kappa_dot_grid']
            self.kappa_grid = arrays['kappa_grid']
            self.tau_rec = arrays['tau_rec']
            self.tau_0 = arrays['tau_0']

    bg_lite = _BGLite(bg_arrays)

    # Import three-phase solver
    from mlx_class.perturbations_sync_class import solve_three_phase

    try:
        combined = solve_three_phase(k, bg_lite, tau_end,
                                     lg_max=lg_max, lp_max=lp_max,
                                     ln_max=ln_max,
                                     method=method, rtol=rtol, atol=atol)
    except Exception as e:
        return None

    if not combined.success:
        return None

    # Evaluate at snapshot points
    N_snap = len(tau_all)
    results = np.zeros((N_snap, 6), dtype=np.float64)

    for it in range(N_snap):
        tau = tau_all[it]
        try:
            gt = combined.eval(tau)
        except Exception:
            continue

        results[it, 0] = gt['Phi_N']
        results[it, 1] = gt['Psi_N']
        results[it, 2] = gt['Theta0_N']
        results[it, 3] = gt['vb_N']
        results[it, 4] = gt['Theta2']
        results[it, 5] = gt['Pi_over_4']

    return results


# ============================================================================
# Full pipeline (parallel)
# ============================================================================

def run_sync_solver_class(N_k=500, k_min=3e-4, k_max=0.35, method='BDF',
                          lg_max=L_GAMMA_MAX, lp_max=L_POL_MAX,
                          ln_max=L_NU_MAX,
                          rtol=1e-6, atol=1e-9,
                          n_workers=None, parallel=True,
                          verbose=True):
    """
    CLASS-style three-phase synchronous gauge pipeline.

    Parameters
    ----------
    N_k : int
        Number of k-modes.
    lg_max : int
        Photon temperature hierarchy truncation (CLASS default: 12).
    lp_max : int
        Photon polarization hierarchy truncation (CLASS default: 10).
    ln_max : int
        Neutrino hierarchy truncation (CLASS default: 17).
    parallel : bool
        Use multiprocessing (default True).
    """
    t_total = time.time()

    if n_workers is None:
        n_workers = min(os.cpu_count() or 4, 8)

    # ================================================================
    # Step 1: Background
    # ================================================================
    if verbose:
        print("=" * 65)
        print("CLASS-Style Three-Phase Boltzmann Solver")
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

    if verbose:
        print(f"\n--- Step 2: Three-Phase Perturbations ({N_k} k-modes, "
              f"{method}) ---")
        print(f"k range: [{k_min:.1e}, {k_max:.2f}] Mpc^-1")
        print(f"Snapshot grid: {N_snap} tau points")
        print(f"  Visibility: {len(tau_vis)} pts, "
              f"tau=[{tau_vis[0]:.1f}, {tau_vis[-1]:.1f}]")
        print(f"  Late ISW:   {len(tau_late)} pts, "
              f"tau=[{tau_late[0]:.1f}, {tau_late[-1]:.1f}]")
        if hasattr(bg, 'z_reio') and bg.tau_reio > 0:
            print(f"  Reionization: z_reio={bg.z_reio:.2f}, "
                  f"tau_reio={bg.tau_reio:.4f}")
        nvar_f = n_var_full(lg_max, lp_max, ln_max)
        print(f"Hierarchy: lg_max={lg_max}, lp_max={lp_max}, ln_max={ln_max}")
        print(f"  Full phase: {nvar_f} vars")
        from .perturbations_sync_class import n_var_tca, N_VAR_RSA
        print(f"  TCA phase:  {n_var_tca(ln_max)} vars")
        print(f"  RSA phase:  {N_VAR_RSA} vars")
        print(f"  Phases: TCA -> Full Hierarchy -> RSA")
        sys.stdout.flush()

    # ================================================================
    # Step 3: Solve all k-modes
    # ================================================================
    t0 = time.time()

    a_snap = np.asarray(bg.a_at_tau(tau_all), dtype=np.float64)
    calH_snap = np.asarray(bg.calH_at_tau(tau_all), dtype=np.float64)

    # Serialize background
    bg_arrays = {
        'tau_grid': bg.tau_grid.copy(),
        'a_grid': bg.a_grid.copy(),
        'calH_grid': bg.calH_grid.copy(),
        'R_grid': bg.R_grid.copy(),
        'kappa_dot_grid': bg.kappa_dot_grid.copy(),
        'kappa_grid': bg.kappa_grid.copy(),
        'tau_rec': float(bg.tau_rec),
        'tau_0': float(bg.tau_0),
    }

    tau_end = float(tau_all[-1] + 1.0)

    if parallel and n_workers > 1:
        # Build work items
        work_items = [
            (k, bg_arrays, tau_end, lg_max, lp_max, ln_max,
             method, rtol, atol, tau_all, a_snap, calH_snap)
            for k in k_arr
        ]

        if verbose:
            print(f"Dispatching {N_k} k-modes to {n_workers} workers...")
            sys.stdout.flush()

        with Pool(processes=n_workers) as pool:
            worker_results = pool.map(_solve_worker_class, work_items)
    else:
        # Sequential
        worker_results = []
        for ik, k in enumerate(k_arr):
            combined = solve_three_phase(k, bg, tau_end,
                                         lg_max=lg_max, lp_max=lp_max,
                                         ln_max=ln_max,
                                         method=method, rtol=rtol, atol=atol)

            if not combined.success:
                worker_results.append(None)
                continue

            results = np.zeros((N_snap, 6), dtype=np.float64)
            for it in range(N_snap):
                tau = tau_all[it]
                try:
                    gt = combined.eval(tau)
                except Exception:
                    continue
                results[it, 0] = gt['Phi_N']
                results[it, 1] = gt['Psi_N']
                results[it, 2] = gt['Theta0_N']
                results[it, 3] = gt['vb_N']
                results[it, 4] = gt['Theta2']
                results[it, 5] = gt['Pi_over_4']

            worker_results.append(results)

            if verbose and (ik + 1) % 50 == 0:
                elapsed = time.time() - t0
                eta_est = elapsed / (ik + 1) * N_k
                print(f"  [{ik+1}/{N_k}] elapsed={elapsed:.1f}s, "
                      f"ETA={eta_est:.1f}s")
                sys.stdout.flush()

    t_pert = time.time() - t0

    # ================================================================
    # Collect results
    # ================================================================
    Phi_N = np.zeros((N_k, N_snap))
    Psi_N = np.zeros((N_k, N_snap))
    Theta0_N = np.zeros((N_k, N_snap))
    vb_N = np.zeros((N_k, N_snap))
    Theta2_arr = np.zeros((N_k, N_snap))
    Pi_arr = np.zeros((N_k, N_snap))

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
        Pi_arr[ik] = res[:, 5]

    if verbose:
        print(f"Perturbations: {t_pert:.1f}s ({n_failed} failed, "
              f"{N_k - n_failed} OK)")

    # ================================================================
    # Step 4: Compute Phi_N' + Psi_N' by cubic spline derivatives
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

    # EE source: full Pi from polarization hierarchy
    pol_prefactor = 0.75  # (3/4) * Pi/4, Pi_arr has Pi/4

    eps_prefactors = np.zeros(N_ell)
    for il, ell in enumerate(ell_values):
        l = int(ell)
        if l >= 2:
            eps_prefactors[il] = np.sqrt(float((l - 1) * l * (l + 1) * (l + 2)))

    if verbose:
        print(f"Computing LOS transfer functions (TT + EE)...")
        print(f"  EE source: full Pi from polarization hierarchy")
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
        S_E = g_snap[it] * pol_prefactor * Pi_arr[:, it]

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
                Delta_l[il, :] += w * (S_SW * jl + S_Dop * jlp
                                       + S_ISW * jl)
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
        'Pi_arr': Pi_arr,
        'l_pol_max': lp_max,
        'timing': {
            'background': t_bg,
            'perturbations': t_pert,
            'los': t_los,
            'total': t_elapsed,
            'n_workers': n_workers if parallel else 1,
        },
    }


# ============================================================================
# CLI entry point
# ============================================================================

def main():
    """CLI entry point: python -m mlx_class.solver_sync_class"""
    import argparse

    parser = argparse.ArgumentParser(
        description='CLASS-style three-phase sync gauge Boltzmann solver')
    parser.add_argument('--N_k', type=int, default=500)
    parser.add_argument('--k_min', type=float, default=3e-4)
    parser.add_argument('--k_max', type=float, default=0.35)
    parser.add_argument('--method', default='BDF')
    parser.add_argument('--lg_max', type=int, default=L_GAMMA_MAX)
    parser.add_argument('--lp_max', type=int, default=L_POL_MAX)
    parser.add_argument('--ln_max', type=int, default=L_NU_MAX)
    parser.add_argument('--rtol', type=float, default=1e-6)
    parser.add_argument('--atol', type=float, default=1e-9)
    parser.add_argument('--parallel', action='store_true', default=True)
    parser.add_argument('--sequential', action='store_true')
    parser.add_argument('--n_workers', type=int, default=None)
    parser.add_argument('--no-plot', action='store_true')
    args = parser.parse_args()

    parallel = not args.sequential

    result = run_sync_solver_class(
        N_k=args.N_k, k_min=args.k_min, k_max=args.k_max,
        method=args.method,
        lg_max=args.lg_max, lp_max=args.lp_max, ln_max=args.ln_max,
        rtol=args.rtol, atol=args.atol,
        n_workers=args.n_workers, parallel=parallel,
        verbose=True)

    # ================================================================
    # Print peak locations and heights
    # ================================================================
    ell = result['ell']
    Dl_TT = result['Dl_TT']

    print("\n--- TT Peak Analysis ---")
    # Find peaks: local maxima
    peaks = []
    for i in range(1, len(Dl_TT) - 1):
        if Dl_TT[i] > Dl_TT[i-1] and Dl_TT[i] > Dl_TT[i+1]:
            peaks.append((int(ell[i]), float(Dl_TT[i])))

    for ip, (l_peak, dl_peak) in enumerate(peaks[:7]):
        print(f"  Peak {ip+1}: l={l_peak}, D_l={dl_peak:.1f} uK^2")

    # ================================================================
    # Compare with CLASS (if available)
    # ================================================================
    try:
        from .class_comparison import load_class_reference
        ref = load_class_reference()
        if ref is not None:
            from scipy.interpolate import interp1d
            ref_interp = interp1d(ref['ell'], ref['Dl_TT'],
                                  kind='cubic', fill_value='extrapolate')
            mask = (ell >= 2) & (ell <= 2500)
            ref_vals = ref_interp(ell[mask])
            our_vals = Dl_TT[mask]
            valid = np.abs(ref_vals) > 1.0
            if np.any(valid):
                rms = np.sqrt(np.mean(
                    ((our_vals[valid] - ref_vals[valid]) / ref_vals[valid])**2))
                print(f"\nRMS error vs CLASS (l=2-2500): {rms*100:.2f}%")
    except Exception:
        pass

    # ================================================================
    # Plot (if matplotlib available and not suppressed)
    # ================================================================
    if not args.no_plot:
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt

            fig, axes = plt.subplots(1, 2, figsize=(14, 5))

            # TT
            ax = axes[0]
            ax.plot(ell, Dl_TT, 'b-', label='CLASS-style 3-phase')
            try:
                from .class_comparison import load_class_reference
                ref = load_class_reference()
                if ref is not None:
                    ax.plot(ref['ell'], ref['Dl_TT'], 'r--',
                            label='CLASS reference', alpha=0.7)
            except Exception:
                pass
            ax.set_xlabel('l')
            ax.set_ylabel(r'$D_l^{TT}$ [$\mu K^2$]')
            ax.set_title('TT Power Spectrum')
            ax.legend()
            ax.set_xlim(2, 2500)

            # EE
            ax = axes[1]
            Dl_EE = result['Dl_EE']
            ax.plot(ell, Dl_EE, 'b-', label='CLASS-style 3-phase')
            try:
                if ref is not None and 'Dl_EE' in ref:
                    ax.plot(ref['ell'], ref['Dl_EE'], 'r--',
                            label='CLASS reference', alpha=0.7)
            except Exception:
                pass
            ax.set_xlabel('l')
            ax.set_ylabel(r'$D_l^{EE}$ [$\mu K^2$]')
            ax.set_title('EE Power Spectrum')
            ax.legend()
            ax.set_xlim(2, 2500)

            plt.tight_layout()
            out_path = 'cl_class_style.png'
            plt.savefig(out_path, dpi=150)
            print(f"\nPlot saved to {out_path}")
            plt.close()
        except ImportError:
            print("\n[NOTE] matplotlib not available, skipping plot")


if __name__ == '__main__':
    main()
