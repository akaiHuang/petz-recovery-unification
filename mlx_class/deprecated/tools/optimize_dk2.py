#!/usr/bin/env python3
"""
optimize_dk2.py — Phase 2: Deep optimization of D(k).

KEY FINDING from Phase 1:
  l_max=25 and l_max=50 give IDENTICAL amplitudes.
  => The deficit is NOT from hierarchy truncation.
  => It's from the IMEX time-stepping (Strang splitting) numerical damping.

APPROACH:
  1. Try denser time grids (2x, 4x) to see if amplitude improves
  2. Try a per-peak D(k) table (no functional form constraint)
  3. Optimize D(k) as a spline to get the absolute floor RMS
  4. Apply best result to solver_accurate.py

Author: Sheng-Kai Huang, 2026
"""
import os
os.environ['MPLCONFIGDIR'] = '/tmp/mpl_config'

import time
import numpy as np
from scipy.optimize import minimize
from scipy.interpolate import interp1d
from scipy.special import spherical_jn


def run_phase2():
    from .background import Background, A_s, n_s, k_pivot, T_CMB, k_eq
    from .solver_accurate import AccurateBoltzmannSolver

    class_path = '/Users/akaihuangm1/Desktop/github/gdm_class_public/output/mlx_ref_lcdm_cl.dat'
    data = np.loadtxt(class_path)
    class_ell = data[:, 0].astype(int)
    class_Dl = data[:, 1] * (T_CMB * 1e6) ** 2

    print("=" * 70)
    print("PHASE 2: DEEP D(k) OPTIMIZATION")
    print("=" * 70)

    bg = Background(khronon=False, recombination='peebles')
    bg.solve()

    k_arr = np.geomspace(5e-4, 0.35, 500).astype(np.float32)

    # ================================================================
    # TEST 1: Dense time grid (2x)
    # ================================================================
    print(f"\n{'='*70}")
    print("TEST 1: Dense time grid (2x)")
    print(f"{'='*70}")

    t0 = time.time()
    solver_2x = AccurateBoltzmannSolver(
        bg, k_arr, l_gamma_max=25,
        grid_density=(400, 800),  # 2x denser
        n_snapshots=300,
        constraint_reset_interval=19
    )
    result_2x = solver_2x.solve()
    t_2x = time.time() - t0

    # Compare amplitudes with standard grid
    t0 = time.time()
    solver_1x = AccurateBoltzmannSolver(
        bg, k_arr, l_gamma_max=25,
        grid_density='fast',
        n_snapshots=300,
        constraint_reset_interval=19
    )
    result_1x = solver_1x.solve()
    t_1x = time.time() - t0

    # Extract at visibility peak
    def get_ThPsi_at_vis_peak(result):
        snap = result.snapshots
        tau_s = snap['tau']
        valid = tau_s > 0
        g_s = bg.visibility_at_tau(tau_s[valid])
        i_peak = np.argmax(g_s)
        return (snap['Theta_0'][valid][i_peak] + snap['Psi'][valid][i_peak])

    ThPsi_1x = get_ThPsi_at_vis_peak(result_1x)
    ThPsi_2x = get_ThPsi_at_vis_peak(result_2x)

    print(f"\n  1x grid: {t_1x:.2f}s, 2x grid: {t_2x:.2f}s")
    print(f"\n  Amplitude comparison (2x/1x):")
    for k_test in [0.02, 0.05, 0.08, 0.10, 0.15, 0.20, 0.25]:
        idx = np.argmin(np.abs(k_arr - k_test))
        ratio = ThPsi_2x[idx] / ThPsi_1x[idx] if abs(ThPsi_1x[idx]) > 1e-10 else 0
        print(f"    k={k_test:.2f}: 1x={ThPsi_1x[idx]:.5f}, "
              f"2x={ThPsi_2x[idx]:.5f}, ratio={ratio:.4f}")

    # ================================================================
    # TEST 2: 4x dense grid
    # ================================================================
    print(f"\n{'='*70}")
    print("TEST 2: Very dense time grid (4x)")
    print(f"{'='*70}")

    t0 = time.time()
    solver_4x = AccurateBoltzmannSolver(
        bg, k_arr, l_gamma_max=25,
        grid_density=(800, 1600),  # 4x denser
        n_snapshots=300,
        constraint_reset_interval=19
    )
    result_4x = solver_4x.solve()
    t_4x = time.time() - t0

    ThPsi_4x = get_ThPsi_at_vis_peak(result_4x)

    print(f"\n  4x grid: {t_4x:.2f}s")
    print(f"\n  Amplitude comparison (Nx/1x):")
    for k_test in [0.02, 0.05, 0.08, 0.10, 0.15, 0.20, 0.25]:
        idx = np.argmin(np.abs(k_arr - k_test))
        r2 = ThPsi_2x[idx] / ThPsi_1x[idx] if abs(ThPsi_1x[idx]) > 1e-10 else 0
        r4 = ThPsi_4x[idx] / ThPsi_1x[idx] if abs(ThPsi_1x[idx]) > 1e-10 else 0
        print(f"    k={k_test:.2f}: 1x={ThPsi_1x[idx]:.5f}, "
              f"2x/1x={r2:.4f}, 4x/1x={r4:.4f}")

    # ================================================================
    # TEST 3: Per-k spline D(k) optimization
    # ================================================================
    print(f"\n{'='*70}")
    print("TEST 3: Spline D(k) optimization against CLASS")
    print(f"{'='*70}")

    # Use the 1x result for speed
    result = result_1x
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

    # Ell grid — use a good sampling
    ell_values = np.unique(np.concatenate([
        np.arange(2, 30, 3),
        np.arange(30, 100, 5),
        np.arange(100, 350, 3),
        np.arange(350, 700, 6),
        np.arange(700, 1200, 10),
        np.arange(1200, 2501, 20),
    ])).astype(int)

    # Use every 2nd ell for optimization speed
    ell_opt = ell_values[::2]

    # Precompute Bessel
    print("Precomputing Bessel functions...")
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
    print(f"Bessel precompute: {time.time()-t0:.1f}s ({len(ell_opt)} ells)")

    class_interp = interp1d(class_ell, class_Dl, kind='linear', fill_value='extrapolate')
    class_at_opt = class_interp(ell_opt)

    def compute_Dl_general(Dk_silk):
        """Compute D_l with arbitrary per-k correction array."""
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

    def rms_general(Dk_silk):
        Dl = compute_Dl_general(Dk_silk)
        mask = (class_at_opt > 100) & (ell_opt > 30)
        if not np.any(mask):
            return 100.0
        rel_err = (Dl[mask] - class_at_opt[mask]) / class_at_opt[mask]
        return 100.0 * np.sqrt(np.mean(rel_err ** 2))

    # Spline-based D(k): define D at N_nodes evenly-spaced k-nodes
    # interpolate to all k with cubic spline
    N_nodes = 15
    k_nodes = np.linspace(0.001, 0.30, N_nodes)

    def params_to_Dk_silk(params):
        """Convert spline node values to per-k D(k)*silk array."""
        D_nodes = params[:N_nodes]
        k_D_factor = params[N_nodes] if len(params) > N_nodes else 1.4

        # Interpolate to all k
        D_interp = interp1d(k_nodes, D_nodes, kind='cubic',
                            fill_value=(D_nodes[0], D_nodes[-1]),
                            bounds_error=False)
        D_k = D_interp(k_arr.astype(np.float64))
        D_k = np.maximum(D_k, 0.1)

        # Transition for superhorizon
        x_eq = k_arr / k_eq
        f_trans = x_eq**2 / (1.0 + x_eq**2)
        D_k_eff = 1.0 + (D_k - 1.0) * f_trans

        k_D_eff = bg.k_D * k_D_factor
        silk_eff = np.exp(-(k_arr / k_D_eff) ** 2)
        return D_k_eff * silk_eff

    def rms_spline(params):
        try:
            Dk_s = params_to_Dk_silk(params)
            return rms_general(Dk_s)
        except:
            return 100.0

    # Initialize with best 3-param result from phase 1
    D0_init = 2.20
    gamma_init = 8.14
    D_init_nodes = D0_init * np.exp(gamma_init * k_nodes**2)
    x0_spline = list(D_init_nodes) + [1.27]  # k_D_factor

    print(f"\n  Initial spline RMS: {rms_spline(x0_spline):.1f}%")

    # Optimize
    print("  Optimizing spline D(k) with Nelder-Mead...")
    res_spline = minimize(rms_spline, x0_spline, method='Nelder-Mead',
                          options={'maxiter': 5000, 'xatol': 0.005, 'fatol': 0.05,
                                   'adaptive': True})
    print(f"  Spline optimized: RMS = {res_spline.fun:.1f}%")

    # Print the spline nodes
    D_opt_nodes = res_spline.x[:N_nodes]
    kDf_opt = res_spline.x[N_nodes]
    print(f"  k_D_factor = {kDf_opt:.4f}")
    print(f"  D(k) spline nodes:")
    for i in range(N_nodes):
        print(f"    k={k_nodes[i]:.4f}: D={D_opt_nodes[i]:.4f}")

    # ================================================================
    # TEST 4: ISW-only term analysis
    # ================================================================
    print(f"\n{'='*70}")
    print("TEST 4: Contribution analysis (SW, Doppler, ISW separately)")
    print(f"{'='*70}")

    # With the optimized spline D(k)
    Dk_silk_opt = params_to_Dk_silk(res_spline.x)

    # SW only
    Dl_sw = np.zeros(len(ell_opt))
    Dl_dop = np.zeros(len(ell_opt))
    Dl_isw = np.zeros(len(ell_opt))
    Dl_all = np.zeros(len(ell_opt))

    for il, ell in enumerate(ell_opt):
        Delta_sw = np.zeros(len(k_arr), dtype=np.float64)
        Delta_dop = np.zeros(len(k_arr), dtype=np.float64)
        Delta_isw = np.zeros(len(k_arr), dtype=np.float64)
        Delta_all = np.zeros(len(k_arr), dtype=np.float64)
        jl_arr = jl_cache[ell]
        jlp_arr = jlp_cache[ell]
        for it in range(N_tau):
            sw = (Theta_0_s[it] + Psi_s[it]) * Dk_silk_opt * jl_arr[it]
            dop = v_b_s[it] * Dk_silk_opt * jlp_arr[it]
            isw = exp_neg_kappa[it] * PhiPsi_dot[it] * jl_arr[it]
            Delta_sw += w_vis[it] * sw
            Delta_dop += w_vis[it] * dop
            Delta_isw += w_vis[it] * isw
            Delta_all += w_vis[it] * (sw + dop + isw)

        for Delta, Dl_out in [(Delta_sw, Dl_sw), (Delta_dop, Dl_dop),
                               (Delta_isw, Dl_isw), (Delta_all, Dl_all)]:
            integrand = P_R * Delta ** 2
            mid = 0.5 * (integrand[:-1] + integrand[1:])
            idx_Dl = il
            val = 4.0 * np.pi * np.sum(mid * dlnk)
            ell_f = float(ell)
            if Delta is Delta_sw:
                Dl_sw[il] = max(val, 0) * ell_f * (ell_f + 1) / (2*np.pi) * (T_CMB*1e6)**2
            elif Delta is Delta_dop:
                Dl_dop[il] = max(val, 0) * ell_f * (ell_f + 1) / (2*np.pi) * (T_CMB*1e6)**2
            elif Delta is Delta_isw:
                Dl_isw[il] = max(val, 0) * ell_f * (ell_f + 1) / (2*np.pi) * (T_CMB*1e6)**2
            else:
                Dl_all[il] = max(val, 0) * ell_f * (ell_f + 1) / (2*np.pi) * (T_CMB*1e6)**2

    # Print at peaks
    from scipy.signal import find_peaks
    from scipy.ndimage import gaussian_filter1d
    Dl_s = gaussian_filter1d(Dl_all, sigma=5)
    mask_pk = ell_opt > 100
    offset_pk = np.argmax(mask_pk)
    pks, _ = find_peaks(Dl_s[mask_pk], distance=20, prominence=10)
    pks = pks + offset_pk

    print(f"\n  Peak analysis:")
    for i, p in enumerate(pks[:7]):
        ell_p = ell_opt[p]
        class_v = class_interp(ell_p)
        print(f"    Peak {i+1}: l={ell_p}, All={Dl_all[p]:.0f}, "
              f"SW={Dl_sw[p]:.0f}, Dop={Dl_dop[p]:.0f}, ISW={Dl_isw[p]:.0f}, "
              f"CLASS={class_v:.0f}, ratio={Dl_all[p]/class_v:.3f}")

    # ================================================================
    # TEST 5: What if we use a completely free per-k multiplicative correction?
    # ================================================================
    print(f"\n{'='*70}")
    print("TEST 5: Per-peak RMS — what if each peak's D were perfect?")
    print(f"{'='*70}")

    # At each CLASS peak, find what D factor would match
    class_smooth = gaussian_filter1d(class_Dl, sigma=15)
    mask_c = class_ell > 100
    offset_c = np.argmax(mask_c)
    pks_c, _ = find_peaks(class_smooth[mask_c], distance=150, prominence=50)
    pks_c = pks_c + offset_c

    print(f"  CLASS peaks:")
    for i, p in enumerate(pks_c[:7]):
        print(f"    Peak {i+1}: l={class_ell[p]}, D_l={class_Dl[p]:.0f}")

    print(f"\n  Our peaks (spline-optimized):")
    for i, p in enumerate(pks[:7]):
        ell_p = ell_opt[p]
        class_v = class_interp(ell_p)
        ratio = np.sqrt(Dl_all[p] / class_v) if class_v > 0 else 0
        print(f"    Peak {i+1}: l={ell_p}, D_l={Dl_all[p]:.0f}, "
              f"CLASS={class_v:.0f}, sqrt(ratio)={ratio:.3f}")

    # ================================================================
    # TEST 6: Apply best spline to solver_accurate and compute
    #          full C_l with DENSE ell sampling for fair RMS
    # ================================================================
    print(f"\n{'='*70}")
    print("TEST 6: Full C_l with optimized spline D(k)")
    print(f"{'='*70}")

    # Full ell grid
    ell_full = ell_values
    class_at_full = class_interp(ell_full)

    # Need full Bessel for this
    print("Computing full C_l (this may take a few minutes)...")
    t0 = time.time()

    Dl_full = np.zeros(len(ell_full))
    for il, ell in enumerate(ell_full):
        Delta_l = np.zeros(len(k_arr), dtype=np.float64)
        for it in range(N_tau):
            x = k_arr * chi_s[it]
            jl = spherical_jn(int(ell), x)
            jlp = spherical_jn(int(ell), x, derivative=True)
            sw = (Theta_0_s[it] + Psi_s[it]) * Dk_silk_opt * jl
            dop = v_b_s[it] * Dk_silk_opt * jlp
            isw = exp_neg_kappa[it] * PhiPsi_dot[it] * jl
            Delta_l += w_vis[it] * (sw + dop + isw)
        integrand = P_R * Delta_l ** 2
        mid = 0.5 * (integrand[:-1] + integrand[1:])
        Dl_full[il] = 4.0 * np.pi * np.sum(mid * dlnk)
        if (il+1) % 50 == 0:
            print(f"  {il+1}/{len(ell_full)} ells")

    Dl_full = np.maximum(Dl_full, 0.0)
    ell_f = ell_full.astype(np.float64)
    Dl_full = ell_f * (ell_f + 1.0) * Dl_full / (2.0 * np.pi) * (T_CMB * 1e6) ** 2

    mask_sig = (class_at_full > 100) & (ell_full > 30)
    rel_err = (Dl_full[mask_sig] - class_at_full[mask_sig]) / class_at_full[mask_sig]
    rms_full = 100.0 * np.sqrt(np.mean(rel_err ** 2))
    max_err = 100.0 * np.max(np.abs(rel_err))

    print(f"\n  Full C_l time: {time.time()-t0:.1f}s")
    print(f"  RMS: {rms_full:.1f}%")
    print(f"  Max: {max_err:.1f}%")

    # Per-band RMS
    for lo, hi in [(100, 500), (500, 1000), (1000, 2000), (2000, 2500)]:
        m = mask_sig & (ell_full >= lo) & (ell_full < hi)
        if np.any(m):
            print(f"  l=[{lo},{hi}): RMS = {100*np.sqrt(np.mean(((Dl_full[m]-class_at_full[m])/class_at_full[m])**2)):.1f}%")

    # ================================================================
    # SUMMARY
    # ================================================================
    print(f"\n{'='*70}")
    print("PHASE 2 FINAL SUMMARY")
    print(f"{'='*70}")
    print(f"\n  Grid density test (amplitude ratio Nx/1x at k=0.1):")
    idx_01 = np.argmin(np.abs(k_arr - 0.1))
    print(f"    1x: {ThPsi_1x[idx_01]:.5f}")
    print(f"    2x: {ThPsi_2x[idx_01]:.5f} ({ThPsi_2x[idx_01]/ThPsi_1x[idx_01]:.4f}x)")
    print(f"    4x: {ThPsi_4x[idx_01]:.5f} ({ThPsi_4x[idx_01]/ThPsi_1x[idx_01]:.4f}x)")
    print(f"\n  RMS vs CLASS:")
    print(f"    Current D(k) [Phase 1]:     24.6%")
    print(f"    3-param optimized:          21.9%")
    print(f"    Spline optimized (full):    {rms_full:.1f}%")
    print(f"\n  Speed: {t_1x:.2f}s (solver) + C_l integration")
    print(f"\n  Spline D(k) nodes (for implementation):")
    for i in range(N_nodes):
        print(f"    k={k_nodes[i]:.4f}: D={D_opt_nodes[i]:.4f}")
    print(f"    k_D_factor = {kDf_opt:.4f}")

    # Plot
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(2, 1, figsize=(14, 10),
                                  gridspec_kw={'height_ratios': [3, 1]}, sharex=True)
        fig.subplots_adjust(hspace=0.06)

        ax = axes[0]
        ax.plot(class_ell, class_Dl, 'k-', lw=1.2, alpha=0.7, label='CLASS')
        ax.plot(ell_full, Dl_full, 'r-', lw=1.5, alpha=0.9,
                label=f'Spline D(k) optimized (RMS={rms_full:.1f}%)')
        ax.set_ylabel(r'$\mathcal{D}_\ell^{TT}$ [$\mu$K$^2$]', fontsize=14)
        ax.set_title('Spline-optimized D(k) vs CLASS', fontsize=13)
        ax.legend(fontsize=10)
        ax.set_xlim(30, 2500)
        ax.grid(alpha=0.15)

        ax2 = axes[1]
        residual = np.where(mask_sig, (Dl_full - class_at_full) / class_at_full * 100, 0)
        ax2.plot(ell_full, residual, 'r-', lw=0.8)
        residual_smooth = gaussian_filter1d(residual, sigma=5)
        ax2.plot(ell_full, residual_smooth, 'r-', lw=1.5, alpha=0.8)
        ax2.axhline(0, color='k', ls='--', lw=0.8)
        ax2.fill_between(ell_full, -10, 10, color='green', alpha=0.08)
        ax2.set_ylabel('Residual [%]', fontsize=11)
        ax2.set_xlabel(r'$\ell$', fontsize=14)
        ax2.set_ylim(-60, 30)
        ax2.grid(alpha=0.15)
        ax2.text(0.02, 0.90, f'RMS = {rms_full:.1f}%',
                 transform=ax2.transAxes, fontsize=10, va='top',
                 bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

        out_dir = os.path.dirname(os.path.abspath(__file__))
        plot_path = os.path.join(out_dir, 'dk_optimization_phase2.png')
        plt.savefig(plot_path, dpi=150, bbox_inches='tight')
        print(f"\n  Plot saved: {plot_path}")
        plt.close()
    except Exception as e:
        print(f"  Plot error: {e}")

    return {
        'rms_spline': rms_full,
        'D_nodes': D_opt_nodes,
        'k_nodes': k_nodes,
        'kDf': kDf_opt,
    }


if __name__ == '__main__':
    run_phase2()
