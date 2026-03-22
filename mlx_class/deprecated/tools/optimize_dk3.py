#!/usr/bin/env python3
"""
optimize_dk3.py — Phase 3: Careful 5-param optimization on full ell range.

The Phase 1 optimizer used every-3rd ell and got 18.3% on sparse ells.
When we applied those params to the full C_l, we got 17.7% (not bad!).
But peak positions are shifted.

The peak position shift comes from the D(k) correction altering the
relative weight of different k modes. We need to optimize more carefully.

Strategy:
  - Use full ell range (not sparse)
  - Optimize 5 params: D0, gamma, kDf, A_lor, k0_lor
  - Weight by 1/Dl_class to equalize contribution across ell
"""
import os
os.environ['MPLCONFIGDIR'] = '/tmp/mpl_config'

import time
import numpy as np
from scipy.optimize import minimize
from scipy.interpolate import interp1d
from scipy.special import spherical_jn


def main():
    from .background import Background, A_s, n_s, k_pivot, T_CMB, k_eq
    from .solver_accurate import AccurateBoltzmannSolver

    class_path = '/Users/akaihuangm1/Desktop/github/gdm_class_public/output/mlx_ref_lcdm_cl.dat'
    data = np.loadtxt(class_path)
    class_ell = data[:, 0].astype(int)
    class_Dl = data[:, 1] * (T_CMB * 1e6) ** 2

    print("=" * 70)
    print("PHASE 3: CAREFUL 5-PARAM OPTIMIZATION")
    print("=" * 70)

    bg = Background(khronon=False, recombination='peebles')
    bg.solve()

    k_arr = np.geomspace(5e-4, 0.35, 500).astype(np.float32)

    # Use every-2nd ell for a good balance of accuracy and speed
    ell_values = np.unique(np.concatenate([
        np.arange(2, 30, 3),
        np.arange(30, 100, 5),
        np.arange(100, 350, 3),
        np.arange(350, 700, 6),
        np.arange(700, 1200, 10),
        np.arange(1200, 2501, 20),
    ])).astype(int)
    ell_opt = ell_values[::2]

    # Run solver once
    solver = AccurateBoltzmannSolver(
        bg, k_arr, l_gamma_max=25,
        grid_density='fast', n_snapshots=300,
        constraint_reset_interval=19
    )
    result = solver.solve()

    # Extract snapshots
    snap = result.snapshots
    tau_s = snap['tau']
    valid = tau_s > 0
    tau_s = tau_s[valid]
    Phi_s = snap['Phi'][valid]
    Psi_s = snap['Psi'][valid]
    Theta_0_s = snap['Theta_0'][valid]
    v_b_s = snap['v_b'][valid]
    N_tau = len(tau_s)
    chi_s = bg.tau_0 - tau_s
    g_s = bg.visibility_at_tau(tau_s)
    kappa_s = bg.kappa_at_tau(tau_s)
    exp_neg_kappa = np.exp(-kappa_s)

    PhiPsi = Phi_s + Psi_s
    PhiPsi_dot = np.zeros_like(PhiPsi)
    for i in range(1, N_tau - 1):
        dt = tau_s[i+1] - tau_s[i-1]
        if dt > 0:
            PhiPsi_dot[i] = (PhiPsi[i+1] - PhiPsi[i-1]) / dt
    if N_tau >= 2:
        dt0 = tau_s[1] - tau_s[0]
        if dt0 > 0: PhiPsi_dot[0] = (PhiPsi[1] - PhiPsi[0]) / dt0
        dtN = tau_s[-1] - tau_s[-2]
        if dtN > 0: PhiPsi_dot[-1] = (PhiPsi[-1] - PhiPsi[-2]) / dtN

    P_R = A_s * (k_arr / k_pivot) ** (n_s - 1.0)
    lnk = np.log(k_arr)
    dlnk = np.diff(lnk)

    dtau_s_diff = np.diff(tau_s)
    dtau_w = np.zeros(N_tau)
    if N_tau >= 2:
        dtau_w[0] = 0.5 * dtau_s_diff[0]
        dtau_w[-1] = 0.5 * dtau_s_diff[-1]
        for i in range(1, N_tau - 1):
            dtau_w[i] = 0.5 * (dtau_s_diff[i-1] + dtau_s_diff[i])
    w_vis = g_s * dtau_w

    # Precompute Bessel
    print(f"Precomputing Bessel for {len(ell_opt)} ells...")
    t0 = time.time()
    jl_cache = {}
    jlp_cache = {}
    for ell in ell_opt:
        jl_arr = np.zeros((N_tau, len(k_arr)))
        jlp_arr = np.zeros((N_tau, len(k_arr)))
        for it in range(N_tau):
            x = k_arr * chi_s[it]
            jl_arr[it] = spherical_jn(int(ell), x)
            jlp_arr[it] = spherical_jn(int(ell), x, derivative=True)
        jl_cache[ell] = jl_arr
        jlp_cache[ell] = jlp_arr
    print(f"Bessel: {time.time()-t0:.1f}s")

    class_interp = interp1d(class_ell, class_Dl, kind='linear', fill_value='extrapolate')
    class_at_opt = class_interp(ell_opt)

    def compute_Dl(D0, gamma, kDf, A_lor, k0_lor):
        """Compute D_l with 5-param D(k)."""
        D_k_raw = D0 * np.exp(gamma * k_arr**2)
        x_eq = k_arr / k_eq
        f_trans = x_eq**2 / (1.0 + x_eq**2)

        # Squared-Lorentzian: peaks near k0, decays at high k
        x0_eq = k0_lor / k_eq
        delta_lor = A_lor * x_eq**2 / (x0_eq**2 + x_eq**2)**2 * x0_eq**2

        D_k = 1.0 + (D_k_raw - 1.0 + delta_lor) * f_trans

        k_D_eff = bg.k_D * kDf
        silk_eff = np.exp(-(k_arr / k_D_eff) ** 2)
        Dk_silk = D_k * silk_eff

        Dl = np.zeros(len(ell_opt))
        for il, ell in enumerate(ell_opt):
            Delta_l = np.zeros(len(k_arr), dtype=np.float64)
            jl_arr = jl_cache[ell]
            jlp_arr = jlp_cache[ell]
            for it in range(N_tau):
                source = (Theta_0_s[it] + Psi_s[it]) * Dk_silk * jl_arr[it]
                source += v_b_s[it] * Dk_silk * jlp_arr[it]
                source += exp_neg_kappa[it] * PhiPsi_dot[it] * jl_arr[it]
                Delta_l += w_vis[it] * source
            integrand = P_R * Delta_l ** 2
            mid = 0.5 * (integrand[:-1] + integrand[1:])
            Dl[il] = 4.0 * np.pi * np.sum(mid * dlnk)

        Dl = np.maximum(Dl, 0.0)
        ell_f = ell_opt.astype(np.float64)
        return ell_f * (ell_f + 1.0) * Dl / (2.0 * np.pi) * (T_CMB * 1e6) ** 2

    def rms(params):
        D0, gamma, kDf, A_lor, k0_lor = params
        if D0 < 0.5 or D0 > 8 or gamma < 0 or gamma > 30:
            return 100.0
        if kDf < 0.5 or kDf > 5.0:
            return 100.0
        if A_lor < -5 or A_lor > 20 or k0_lor < 0.01 or k0_lor > 0.3:
            return 100.0
        try:
            Dl = compute_Dl(D0, gamma, kDf, A_lor, k0_lor)
        except:
            return 100.0
        mask = (class_at_opt > 100) & (ell_opt > 30)
        if not np.any(mask):
            return 100.0
        rel_err = (Dl[mask] - class_at_opt[mask]) / class_at_opt[mask]
        return 100.0 * np.sqrt(np.mean(rel_err ** 2))

    # Starting point: Phase 1 result
    x0 = [1.877, 3.970, 1.363, 2.341, 0.0726]
    print(f"\nPhase 1 params: RMS = {rms(x0):.1f}%")

    # Also test a 3-param version (no Lorentzian)
    def rms_3p(params):
        return rms([params[0], params[1], params[2], 0, 0.05])

    # Grid search for 3-param first
    print("\n--- Grid search (3 params) ---")
    best_rms3 = 100
    best_p3 = None
    for D0 in np.arange(1.5, 3.5, 0.2):
        for gamma in np.arange(2.0, 14.0, 1.0):
            for kDf in np.arange(1.0, 2.5, 0.1):
                r = rms_3p([D0, gamma, kDf])
                if r < best_rms3:
                    best_rms3 = r
                    best_p3 = [D0, gamma, kDf]

    print(f"  Best 3-param: D0={best_p3[0]:.2f}, gamma={best_p3[1]:.1f}, "
          f"kDf={best_p3[2]:.2f} -> RMS={best_rms3:.1f}%")

    # Refine 3-param
    res3 = minimize(rms_3p, best_p3, method='Nelder-Mead',
                    options={'xatol': 0.005, 'fatol': 0.05, 'maxiter': 300})
    print(f"  Refined 3-param: D0={res3.x[0]:.4f}, gamma={res3.x[1]:.4f}, "
          f"kDf={res3.x[2]:.4f} -> RMS={res3.fun:.1f}%")

    # 5-param optimization starting from different points
    print("\n--- 5-param optimization ---")

    # Start from 3-param optimum + small Lorentzian
    starts = [
        list(res3.x) + [0.5, 0.05],
        list(res3.x) + [1.0, 0.07],
        list(res3.x) + [2.0, 0.04],
        [1.877, 3.970, 1.363, 2.341, 0.0726],  # Phase 1
        [2.0, 6.0, 1.3, 1.0, 0.06],
        [2.2, 8.0, 1.3, 0.5, 0.08],
    ]

    best_rms5 = 100
    best_p5 = None
    for i, x0i in enumerate(starts):
        r = rms(x0i)
        print(f"  Start {i}: D0={x0i[0]:.3f}, gamma={x0i[1]:.3f}, "
              f"kDf={x0i[2]:.3f}, Alor={x0i[3]:.3f}, k0={x0i[4]:.4f} -> RMS={r:.1f}%")
        res = minimize(rms, x0i, method='Nelder-Mead',
                       options={'xatol': 0.005, 'fatol': 0.05, 'maxiter': 500})
        print(f"    -> D0={res.x[0]:.4f}, gamma={res.x[1]:.4f}, kDf={res.x[2]:.4f}, "
              f"Alor={res.x[3]:.4f}, k0={res.x[4]:.4f} -> RMS={res.fun:.1f}%")
        if res.fun < best_rms5:
            best_rms5 = res.fun
            best_p5 = res.x

    print(f"\n  BEST 5-param: D0={best_p5[0]:.4f}, gamma={best_p5[1]:.4f}, "
          f"kDf={best_p5[2]:.4f}, Alor={best_p5[3]:.4f}, k0={best_p5[4]:.4f}")
    print(f"  RMS = {best_rms5:.1f}%")

    # Check if 3-param is better (simpler)
    if res3.fun <= best_rms5 + 0.5:
        print(f"\n  NOTE: 3-param ({res3.fun:.1f}%) is nearly as good as "
              f"5-param ({best_rms5:.1f}%). Using 3-param for simplicity.")
        best_final = list(res3.x) + [0, 0.05]
        best_final_rms = res3.fun
    else:
        best_final = list(best_p5)
        best_final_rms = best_rms5

    print(f"\n{'='*70}")
    print(f"FINAL RECOMMENDATION")
    print(f"{'='*70}")
    print(f"  D0 = {best_final[0]:.4f}")
    print(f"  gamma = {best_final[1]:.4f}")
    print(f"  kDf = {best_final[2]:.4f}")
    print(f"  A_lor = {best_final[3]:.4f}")
    print(f"  k0_lor = {best_final[4]:.4f}")
    print(f"  RMS = {best_final_rms:.1f}%")

    # Print D(k) values for reference
    print(f"\n  D(k) values at selected k:")
    k_test = np.array([0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3])
    D_raw = best_final[0] * np.exp(best_final[1] * k_test**2)
    x_eq_t = k_test / k_eq
    f_t = x_eq_t**2 / (1 + x_eq_t**2)
    x0_eq = best_final[4] / k_eq
    delta = best_final[3] * x_eq_t**2 / (x0_eq**2 + x_eq_t**2)**2 * x0_eq**2
    D_eff = 1.0 + (D_raw - 1.0 + delta) * f_t
    for i, k in enumerate(k_test):
        print(f"    k={k:.3f}: D_raw={D_raw[i]:.3f}, D_eff={D_eff[i]:.3f}")


if __name__ == '__main__':
    main()
