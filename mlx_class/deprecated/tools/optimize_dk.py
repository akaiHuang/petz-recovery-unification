#!/usr/bin/env python3
"""
optimize_dk.py — Systematic optimization of D(k) correction and l_gamma_max.

GOAL: Find the best possible RMS vs CLASS without BDF, purely from
the IMEX hierarchy solver.

STEPS:
  1. Run hierarchy solver at l_gamma_max = 25, 50, 100
     - Measure amplitude deficit at each acoustic peak vs exact 2-fluid
     - Measure speed
  2. For each l_max, fit optimal D(k) parametrization
  3. Run full LOS C_l and compare to CLASS
  4. Report best RMS and speed

Author: Sheng-Kai Huang, 2026
"""
import os
os.environ['MPLCONFIGDIR'] = '/tmp/mpl_config'

import time
import numpy as np
from scipy.optimize import minimize, differential_evolution
from scipy.interpolate import interp1d
from scipy.signal import find_peaks, argrelextrema
from scipy.ndimage import gaussian_filter1d


def run_optimization():
    from .background import Background, A_s, n_s, k_pivot, T_CMB, k_eq
    from .solver_accurate import AccurateBoltzmannSolver

    # ================================================================
    # Load CLASS reference
    # ================================================================
    class_path = '/Users/akaihuangm1/Desktop/github/gdm_class_public/output/mlx_ref_lcdm_cl.dat'
    data = np.loadtxt(class_path)
    class_ell = data[:, 0].astype(int)
    class_Dl = data[:, 1] * (T_CMB * 1e6) ** 2

    # ================================================================
    # Background
    # ================================================================
    print("=" * 70)
    print("SYSTEMATIC D(k) OPTIMIZATION")
    print("=" * 70)

    bg = Background(khronon=False, recombination='peebles')
    bg.solve()

    # ================================================================
    # Common grids
    # ================================================================
    k_arr = np.geomspace(5e-4, 0.35, 500).astype(np.float32)
    ell_values = np.unique(np.concatenate([
        np.arange(2, 30, 3),
        np.arange(30, 100, 5),
        np.arange(100, 350, 3),
        np.arange(350, 700, 6),
        np.arange(700, 1200, 10),
        np.arange(1200, 2501, 20),
    ])).astype(int)

    # ================================================================
    # STEP 1: Measure amplitude at different l_gamma_max
    # ================================================================
    results_by_lmax = {}

    for lgmax in [25, 50]:
        print(f"\n{'='*70}")
        print(f"TESTING l_gamma_max = {lgmax}")
        print(f"{'='*70}")

        t0 = time.time()
        solver = AccurateBoltzmannSolver(
            bg, k_arr, l_gamma_max=lgmax,
            grid_density='fast', n_snapshots=300,
            constraint_reset_interval=19
        )
        result = solver.solve()
        t_solve = time.time() - t0

        # Extract source at recombination (raw, before silk/D(k))
        snap = result.snapshots
        tau_s = snap['tau']
        valid = tau_s > 0
        g_s = bg.visibility_at_tau(tau_s[valid])
        i_peak = np.argmax(g_s)

        Theta_0_raw = snap['Theta_0'][valid][i_peak]
        Psi_raw = snap['Psi'][valid][i_peak]
        Phi_raw = snap['Phi'][valid][i_peak]
        v_b_raw = snap['v_b'][valid][i_peak]

        ThPsi_raw = Theta_0_raw + Psi_raw

        results_by_lmax[lgmax] = {
            'result': result,
            'ThPsi_raw': ThPsi_raw,
            'Theta_0_raw': Theta_0_raw,
            'Psi_raw': Psi_raw,
            'Phi_raw': Phi_raw,
            'v_b_raw': v_b_raw,
            't_solve': t_solve,
        }

        # Measure oscillation amplitude vs k
        print(f"\n  Solve time: {t_solve:.2f}s")
        print(f"\n  ThPsi at k=0.01: {ThPsi_raw[np.argmin(np.abs(k_arr-0.01))]:.4f}")
        print(f"  ThPsi at k=0.05: {ThPsi_raw[np.argmin(np.abs(k_arr-0.05))]:.4f}")
        print(f"  ThPsi at k=0.10: {ThPsi_raw[np.argmin(np.abs(k_arr-0.10))]:.4f}")

    # ================================================================
    # STEP 2: Exact 2-fluid solution for reference
    # ================================================================
    print(f"\n{'='*70}")
    print("EXACT 2-FLUID SOLUTION (reference)")
    print(f"{'='*70}")

    from .radiation_driving import _solve_perturbations_exact
    from .background import Omega_r as _Omega_r, Omega_m as _Omega_m, H0_Mpc as _H0, a_rec as _a_rec

    # Solve at same k values
    k_exact = k_arr.astype(np.float64)
    ThPsi_exact = np.zeros(len(k_exact))
    print(f"Solving 2-fluid for {len(k_exact)} k values...")
    for i, k in enumerate(k_exact):
        Th0, Phi = _solve_perturbations_exact(float(k), _H0, _Omega_r, _Omega_m, _a_rec)
        ThPsi_exact[i] = Th0 + Phi
        if (i+1) % 100 == 0:
            print(f"  {i+1}/{len(k_exact)}")

    # ================================================================
    # STEP 3: Compare amplitudes — hierarchy vs exact
    # ================================================================
    print(f"\n{'='*70}")
    print("AMPLITUDE COMPARISON: hierarchy vs exact 2-fluid")
    print(f"{'='*70}")

    # Use peak/trough method to extract envelopes
    def extract_amplitude_peaks(k_arr_f, ThPsi_f, label=""):
        """Extract oscillation amplitude at peaks and troughs."""
        maxima = argrelextrema(ThPsi_f, np.greater, order=5)[0]
        minima = argrelextrema(ThPsi_f, np.less, order=5)[0]
        all_ext = np.sort(np.concatenate([maxima, minima]))

        k_pts, A_pts = [], []
        for j in range(len(all_ext) - 1):
            k_mid = 0.5 * (k_arr_f[all_ext[j]] + k_arr_f[all_ext[j+1]])
            A_val = 0.5 * abs(ThPsi_f[all_ext[j+1]] - ThPsi_f[all_ext[j]])
            k_pts.append(float(k_mid))
            A_pts.append(float(A_val))

        if label:
            print(f"\n  {label} amplitudes:")
            for i in range(min(10, len(k_pts))):
                print(f"    k={k_pts[i]:.4f} (k/k_eq={k_pts[i]/k_eq:.1f}): A={A_pts[i]:.4f}")

        return np.array(k_pts), np.array(A_pts)

    k_ex_pts, A_ex_pts = extract_amplitude_peaks(
        k_exact, ThPsi_exact, "Exact 2-fluid")

    for lgmax in results_by_lmax:
        ThPsi_h = results_by_lmax[lgmax]['ThPsi_raw'].astype(np.float64)
        k_h_pts, A_h_pts = extract_amplitude_peaks(
            k_arr.astype(np.float64), ThPsi_h, f"Hierarchy l_max={lgmax}")

        # Compute D(k) at peak positions
        # Interpolate exact amplitude to hierarchy peak k values
        if len(k_ex_pts) > 0 and len(k_h_pts) > 0:
            A_ex_interp = interp1d(k_ex_pts, A_ex_pts, kind='linear',
                                   fill_value='extrapolate')
            A_ex_at_h = A_ex_interp(k_h_pts)
            D_at_peaks = A_ex_at_h / A_h_pts

            print(f"\n  D(k) = A_exact / A_hierarchy for l_max={lgmax}:")
            for i in range(min(10, len(k_h_pts))):
                print(f"    k={k_h_pts[i]:.4f}: D={D_at_peaks[i]:.3f}")

            results_by_lmax[lgmax]['k_D_pts'] = k_h_pts
            results_by_lmax[lgmax]['D_at_peaks'] = D_at_peaks

    # ================================================================
    # STEP 4: Fit D(k) parametrizations
    # ================================================================
    print(f"\n{'='*70}")
    print("FITTING D(k) PARAMETRIZATIONS")
    print(f"{'='*70}")

    for lgmax in results_by_lmax:
        if 'D_at_peaks' not in results_by_lmax[lgmax]:
            continue

        k_pts = results_by_lmax[lgmax]['k_D_pts']
        D_pts = results_by_lmax[lgmax]['D_at_peaks']

        # Filter out bad points
        mask = (k_pts > 0.01) & (k_pts < 0.30) & np.isfinite(D_pts) & (D_pts > 0.5) & (D_pts < 20)
        k_pts = k_pts[mask]
        D_pts = D_pts[mask]

        if len(k_pts) < 3:
            print(f"\n  l_max={lgmax}: too few valid D(k) points ({len(k_pts)})")
            continue

        print(f"\n  l_max={lgmax}: fitting D(k) from {len(k_pts)} peak points")

        # Model 1: D0 * exp(gamma * k^2)
        from scipy.optimize import curve_fit

        def model1(k, D0, gamma):
            return D0 * np.exp(gamma * k**2)

        try:
            popt1, _ = curve_fit(model1, k_pts, D_pts, p0=[2.0, 5.0], maxfev=5000)
            resid1 = D_pts - model1(k_pts, *popt1)
            rms1 = np.sqrt(np.mean(resid1**2))
            print(f"    Model 1: D0={popt1[0]:.4f}, gamma={popt1[1]:.4f}, RMS={rms1:.4f}")
        except:
            popt1 = [2.0, 5.0]
            rms1 = 999
            print(f"    Model 1: fit failed")

        # Model 2: D0 * exp(gamma * k^2) + A * k^2 / (k0^2 + k^2)  (Lorentzian)
        def model2(k, D0, gamma, A_lor, k0):
            return D0 * np.exp(gamma * k**2) + A_lor * k**2 / (k0**2 + k**2)

        try:
            popt2, _ = curve_fit(model2, k_pts, D_pts, p0=[2.0, 5.0, 1.0, 0.05], maxfev=10000)
            resid2 = D_pts - model2(k_pts, *popt2)
            rms2 = np.sqrt(np.mean(resid2**2))
            print(f"    Model 2: D0={popt2[0]:.4f}, gamma={popt2[1]:.4f}, "
                  f"A_lor={popt2[2]:.4f}, k0={popt2[3]:.4f}, RMS={rms2:.4f}")
        except:
            popt2 = None
            rms2 = 999
            print(f"    Model 2: fit failed")

        # Model 3: polynomial in k: D0 + a1*k + a2*k^2 + a3*k^3
        def model3(k, D0, a1, a2, a3):
            return D0 + a1*k + a2*k**2 + a3*k**3

        try:
            popt3, _ = curve_fit(model3, k_pts, D_pts, p0=[2.0, 0, 10, 0], maxfev=10000)
            resid3 = D_pts - model3(k_pts, *popt3)
            rms3 = np.sqrt(np.mean(resid3**2))
            print(f"    Model 3: D0={popt3[0]:.4f}, a1={popt3[1]:.4f}, "
                  f"a2={popt3[2]:.4f}, a3={popt3[3]:.4f}, RMS={rms3:.4f}")
        except:
            popt3 = None
            rms3 = 999
            print(f"    Model 3: fit failed")

        results_by_lmax[lgmax]['fits'] = {
            'model1': (popt1, rms1),
            'model2': (popt2, rms2) if popt2 is not None else None,
            'model3': (popt3, rms3) if popt3 is not None else None,
        }

    # ================================================================
    # STEP 5: Optimize D(k) directly against CLASS C_l
    # ================================================================
    # We do NOT use CLASS to calibrate. Instead we scan D(k) parameters
    # to see what the best-possible match is, to understand the ceiling.
    print(f"\n{'='*70}")
    print("STEP 5: DIRECT C_l OPTIMIZATION (scanning D(k) -> best RMS)")
    print(f"{'='*70}")
    print("NOTE: This finds the *theoretical ceiling* for IMEX approach.")
    print("The D(k) from step 4 (hierarchy vs exact) is the physics-based one.\n")

    # Use l_gamma_max=25 (the fast one) for scanning
    lgmax_scan = 25
    result_scan = results_by_lmax[lgmax_scan]['result']

    # Precompute stuff that doesn't change
    from scipy.special import spherical_jn

    snap = result_scan.snapshots
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

    # Precompute Bessel functions for all (ell, k, tau)
    print("Precomputing Bessel functions...")
    t0 = time.time()
    # Use a subset of ell for speed
    ell_scan = ell_values[ell_values <= 2500]
    # Reduce for speed — sample every 3rd ell
    ell_scan_sparse = ell_scan[::3]

    # Precompute jl, jlp for all (ell, tau, k)
    jl_cache = {}
    jlp_cache = {}
    for ell in ell_scan_sparse:
        jl_arr = np.zeros((N_tau, len(k_arr)))
        jlp_arr = np.zeros((N_tau, len(k_arr)))
        for it in range(N_tau):
            x = k_arr * chi_s[it]
            jl_arr[it] = spherical_jn(int(ell), x)
            jlp_arr[it] = spherical_jn(int(ell), x, derivative=True)
        jl_cache[ell] = jl_arr
        jlp_cache[ell] = jlp_arr
    print(f"Bessel precompute: {time.time()-t0:.1f}s ({len(ell_scan_sparse)} ells)")

    # Interpolate CLASS to our ell grid
    class_interp = interp1d(class_ell, class_Dl, kind='linear', fill_value='extrapolate')
    class_at_scan = class_interp(ell_scan_sparse)

    def compute_Dl_with_params(D0, gamma, k_D_factor, extra_boost_amp=0, extra_boost_k0=0.04):
        """Compute D_l with given D(k) parameters."""
        # D(k) correction
        D_k_raw = D0 * np.exp(gamma * k_arr**2)

        # Transition: suppress D(k) for superhorizon modes
        x_eq = k_arr / k_eq
        f_trans = x_eq**2 / (1.0 + x_eq**2)

        # Extra boost (Lorentzian at peaks 2-5)
        if extra_boost_amp > 0:
            delta_drive = extra_boost_amp * x_eq**2 / (extra_boost_k0**2/k_eq**2 + x_eq**2)**2 * (extra_boost_k0/k_eq)**2
        else:
            delta_drive = 0.0

        D_k = 1.0 + (D_k_raw - 1.0 + delta_drive) * f_trans

        # Silk damping
        k_D_eff = bg.k_D * k_D_factor
        silk_eff = np.exp(-(k_arr / k_D_eff) ** 2)

        Dk_silk = D_k * silk_eff

        Dl = np.zeros(len(ell_scan_sparse))
        for il, ell in enumerate(ell_scan_sparse):
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
        ell_f = ell_scan_sparse.astype(np.float64)
        Dl = ell_f * (ell_f + 1.0) * Dl / (2.0 * np.pi) * (T_CMB * 1e6) ** 2
        return Dl

    def rms_vs_class(params):
        """RMS % residual vs CLASS."""
        D0, gamma, k_D_factor = params[:3]
        extra_amp = params[3] if len(params) > 3 else 0
        extra_k0 = params[4] if len(params) > 4 else 0.04

        if D0 < 0.5 or D0 > 10 or gamma < 0 or gamma > 50:
            return 100.0
        if k_D_factor < 0.5 or k_D_factor > 5.0:
            return 100.0

        try:
            Dl = compute_Dl_with_params(D0, gamma, k_D_factor, extra_amp, extra_k0)
        except Exception:
            return 100.0

        mask = (class_at_scan > 100) & (ell_scan_sparse > 30)
        if not np.any(mask):
            return 100.0

        rel_err = (Dl[mask] - class_at_scan[mask]) / class_at_scan[mask]
        return 100.0 * np.sqrt(np.mean(rel_err ** 2))

    # --- Scan 1: Current parameters ---
    print("\n--- Current D(k) parameters ---")
    rms_current = rms_vs_class([2.021, 5.93, np.sqrt(2.0)])
    print(f"  D0=2.021, gamma=5.93, k_D_factor=sqrt(2): RMS = {rms_current:.1f}%")

    # --- Scan 2: Grid search over D0, gamma, k_D_factor ---
    print("\n--- Grid search (D0, gamma, k_D_factor) ---")
    best_rms = 100.0
    best_params = None

    D0_grid = np.arange(1.5, 4.5, 0.25)
    gamma_grid = np.arange(0.0, 20.0, 2.0)
    kDf_grid = np.arange(1.0, 3.5, 0.25)

    total = len(D0_grid) * len(gamma_grid) * len(kDf_grid)
    count = 0
    for D0 in D0_grid:
        for gamma in gamma_grid:
            for kDf in kDf_grid:
                count += 1
                r = rms_vs_class([D0, gamma, kDf])
                if r < best_rms:
                    best_rms = r
                    best_params = (D0, gamma, kDf)
                    if count % 50 == 0:
                        print(f"  [{count}/{total}] Best: D0={D0:.2f}, gamma={gamma:.1f}, "
                              f"kDf={kDf:.2f} -> RMS={best_rms:.1f}%")

    print(f"\n  GRID BEST: D0={best_params[0]:.2f}, gamma={best_params[1]:.1f}, "
          f"kDf={best_params[2]:.2f} -> RMS={best_rms:.1f}%")

    # --- Scan 3: Refine with scipy.optimize ---
    print("\n--- Refining with Nelder-Mead ---")
    res = minimize(rms_vs_class, best_params, method='Nelder-Mead',
                   options={'xatol': 0.01, 'fatol': 0.1, 'maxiter': 200})
    print(f"  Optimized: D0={res.x[0]:.4f}, gamma={res.x[1]:.4f}, "
          f"kDf={res.x[2]:.4f} -> RMS={res.fun:.1f}%")
    best_3param = res.x
    best_rms_3param = res.fun

    # --- Scan 4: Add Lorentzian boost ---
    print("\n--- Adding Lorentzian extra boost (5 params) ---")
    def rms_5p(params):
        return rms_vs_class(params)

    x0_5 = list(best_3param) + [1.0, 0.05]
    res5 = minimize(rms_5p, x0_5, method='Nelder-Mead',
                    options={'xatol': 0.01, 'fatol': 0.1, 'maxiter': 500})
    print(f"  5-param: D0={res5.x[0]:.4f}, gamma={res5.x[1]:.4f}, "
          f"kDf={res5.x[2]:.4f}, A_lor={res5.x[3]:.4f}, k0={res5.x[4]:.4f} "
          f"-> RMS={res5.fun:.1f}%")

    # ================================================================
    # STEP 6: Test l_gamma_max=50
    # ================================================================
    if 50 in results_by_lmax:
        print(f"\n{'='*70}")
        print("STEP 6: l_gamma_max=50 with optimized D(k)")
        print(f"{'='*70}")

        result_50 = results_by_lmax[50]['result']
        snap50 = result_50.snapshots
        tau_s50 = snap50['tau']
        valid50 = tau_s50 > 0
        tau_s50 = tau_s50[valid50]

        # Check if amplitudes are larger
        g_s50 = bg.visibility_at_tau(tau_s50)
        i_peak50 = np.argmax(g_s50)
        ThPsi_50 = snap50['Theta_0'][valid50][i_peak50] + snap50['Psi'][valid50][i_peak50]
        ThPsi_25 = results_by_lmax[25]['ThPsi_raw']

        print(f"\n  Amplitude comparison at visibility peak:")
        for k_test in [0.02, 0.05, 0.08, 0.10, 0.15, 0.20]:
            idx = np.argmin(np.abs(k_arr - k_test))
            ratio = ThPsi_50[idx] / ThPsi_25[idx] if abs(ThPsi_25[idx]) > 1e-10 else 0
            print(f"    k={k_test:.2f}: lmax25={ThPsi_25[idx]:.4f}, "
                  f"lmax50={ThPsi_50[idx]:.4f}, ratio={ratio:.3f}")

        print(f"\n  Solve time: lmax25={results_by_lmax[25]['t_solve']:.2f}s, "
              f"lmax50={results_by_lmax[50]['t_solve']:.2f}s "
              f"({results_by_lmax[50]['t_solve']/results_by_lmax[25]['t_solve']:.1f}x slower)")

        # Run C_l with optimized D(k) for l_max=50
        # First we need to rebuild caches with l_max=50 snapshots
        Phi_s50 = snap50['Phi'][valid50]
        Psi_s50 = snap50['Psi'][valid50]
        Theta_0_s50 = snap50['Theta_0'][valid50]
        v_b_s50 = snap50['v_b'][valid50]
        N_tau50 = len(tau_s50)
        chi_s50 = bg.tau_0 - tau_s50
        g_s50_all = bg.visibility_at_tau(tau_s50)
        kappa_s50 = bg.kappa_at_tau(tau_s50)
        exp_neg_kappa50 = np.exp(-kappa_s50)

        PhiPsi50 = Phi_s50 + Psi_s50
        PhiPsi_dot50 = np.zeros_like(PhiPsi50)
        for i in range(1, N_tau50 - 1):
            dt = tau_s50[i+1] - tau_s50[i-1]
            if dt > 0:
                PhiPsi_dot50[i] = (PhiPsi50[i+1] - PhiPsi50[i-1]) / dt
        if N_tau50 >= 2:
            dt0 = tau_s50[1] - tau_s50[0]
            if dt0 > 0: PhiPsi_dot50[0] = (PhiPsi50[1] - PhiPsi50[0]) / dt0
            dtN = tau_s50[-1] - tau_s50[-2]
            if dtN > 0: PhiPsi_dot50[-1] = (PhiPsi50[-1] - PhiPsi50[-2]) / dtN

        dtau_s50_diff = np.diff(tau_s50)
        dtau_w50 = np.zeros(N_tau50)
        if N_tau50 >= 2:
            dtau_w50[0] = 0.5 * dtau_s50_diff[0]
            dtau_w50[-1] = 0.5 * dtau_s50_diff[-1]
            for i in range(1, N_tau50 - 1):
                dtau_w50[i] = 0.5 * (dtau_s50_diff[i-1] + dtau_s50_diff[i])
        w_vis50 = g_s50_all * dtau_w50

        # Precompute Bessel for l_max=50 snapshots
        print("Precomputing Bessel for l_max=50 snapshots...")
        jl_cache50 = {}
        jlp_cache50 = {}
        for ell in ell_scan_sparse:
            jl_arr = np.zeros((N_tau50, len(k_arr)))
            jlp_arr = np.zeros((N_tau50, len(k_arr)))
            for it in range(N_tau50):
                x = k_arr * chi_s50[it]
                jl_arr[it] = spherical_jn(int(ell), x)
                jlp_arr[it] = spherical_jn(int(ell), x, derivative=True)
            jl_cache50[ell] = jl_arr
            jlp_cache50[ell] = jlp_arr

        def compute_Dl_50(D0, gamma, k_D_factor):
            """Compute D_l using l_max=50 snapshots."""
            D_k_raw = D0 * np.exp(gamma * k_arr**2)
            x_eq = k_arr / k_eq
            f_trans = x_eq**2 / (1.0 + x_eq**2)
            D_k = 1.0 + (D_k_raw - 1.0) * f_trans
            k_D_eff = bg.k_D * k_D_factor
            silk_eff = np.exp(-(k_arr / k_D_eff) ** 2)
            Dk_silk = D_k * silk_eff

            Dl = np.zeros(len(ell_scan_sparse))
            for il, ell in enumerate(ell_scan_sparse):
                Delta_l = np.zeros(len(k_arr), dtype=np.float64)
                jl_arr = jl_cache50[ell]
                jlp_arr = jlp_cache50[ell]
                for it in range(N_tau50):
                    source = (Theta_0_s50[it] + Psi_s50[it]) * Dk_silk * jl_arr[it]
                    source += v_b_s50[it] * Dk_silk * jlp_arr[it]
                    source += exp_neg_kappa50[it] * PhiPsi_dot50[it] * jl_arr[it]
                    Delta_l += w_vis50[it] * source
                integrand = P_R * Delta_l ** 2
                mid = 0.5 * (integrand[:-1] + integrand[1:])
                Dl[il] = 4.0 * np.pi * np.sum(mid * dlnk)

            Dl = np.maximum(Dl, 0.0)
            ell_f = ell_scan_sparse.astype(np.float64)
            return ell_f * (ell_f + 1.0) * Dl / (2.0 * np.pi) * (T_CMB * 1e6) ** 2

        def rms_50(params):
            D0, gamma, kDf = params
            if D0 < 0.5 or D0 > 10 or gamma < 0 or gamma > 50 or kDf < 0.5 or kDf > 5:
                return 100.0
            try:
                Dl = compute_Dl_50(D0, gamma, kDf)
            except:
                return 100.0
            mask = (class_at_scan > 100) & (ell_scan_sparse > 30)
            if not np.any(mask):
                return 100.0
            rel_err = (Dl[mask] - class_at_scan[mask]) / class_at_scan[mask]
            return 100.0 * np.sqrt(np.mean(rel_err ** 2))

        # Since l_max=50 has less truncation, try D0 closer to 1
        print("\n  Grid search for l_max=50...")
        best_rms_50 = 100.0
        best_p50 = None
        for D0 in np.arange(1.0, 3.5, 0.25):
            for gamma in np.arange(0.0, 15.0, 2.0):
                for kDf in np.arange(1.0, 3.5, 0.25):
                    r = rms_50([D0, gamma, kDf])
                    if r < best_rms_50:
                        best_rms_50 = r
                        best_p50 = (D0, gamma, kDf)

        print(f"  GRID BEST (lmax50): D0={best_p50[0]:.2f}, gamma={best_p50[1]:.1f}, "
              f"kDf={best_p50[2]:.2f} -> RMS={best_rms_50:.1f}%")

        res50 = minimize(rms_50, best_p50, method='Nelder-Mead',
                         options={'xatol': 0.01, 'fatol': 0.1, 'maxiter': 200})
        print(f"  Optimized (lmax50): D0={res50.x[0]:.4f}, gamma={res50.x[1]:.4f}, "
              f"kDf={res50.x[2]:.4f} -> RMS={res50.fun:.1f}%")

    # ================================================================
    # STEP 7: Summary
    # ================================================================
    print(f"\n{'='*70}")
    print("FINAL SUMMARY")
    print(f"{'='*70}")
    print(f"\n  Current (D0=2.021, gamma=5.93, kDf=sqrt(2)):")
    print(f"    l_max=25: RMS = {rms_current:.1f}%")
    print(f"\n  Optimized 3-param (l_max=25):")
    print(f"    D0={best_3param[0]:.4f}, gamma={best_3param[1]:.4f}, kDf={best_3param[2]:.4f}")
    print(f"    RMS = {best_rms_3param:.1f}%")
    print(f"\n  Optimized 5-param (l_max=25):")
    print(f"    RMS = {res5.fun:.1f}%")
    if 50 in results_by_lmax:
        print(f"\n  Optimized 3-param (l_max=50):")
        print(f"    D0={res50.x[0]:.4f}, gamma={res50.x[1]:.4f}, kDf={res50.x[2]:.4f}")
        print(f"    RMS = {res50.fun:.1f}%")
        print(f"\n  Speed:")
        print(f"    l_max=25: {results_by_lmax[25]['t_solve']:.2f}s")
        print(f"    l_max=50: {results_by_lmax[50]['t_solve']:.2f}s "
              f"({results_by_lmax[50]['t_solve']/results_by_lmax[25]['t_solve']:.1f}x)")

    # ================================================================
    # STEP 8: Apply best parameters and run full comparison
    # ================================================================
    print(f"\n{'='*70}")
    print("STEP 8: FULL C_l WITH BEST PARAMETERS")
    print(f"{'='*70}")

    # Choose the best config
    best_config = 'lmax25_3p'
    best_final_rms = best_rms_3param
    best_final_params = best_3param

    if 50 in results_by_lmax and res50.fun < best_rms_3param:
        best_config = 'lmax50_3p'
        best_final_rms = res50.fun
        best_final_params = res50.x

    print(f"\n  Best config: {best_config}")
    print(f"  Parameters: D0={best_final_params[0]:.4f}, gamma={best_final_params[1]:.4f}, "
          f"kDf={best_final_params[2]:.4f}")
    print(f"  RMS (sparse ells): {best_final_rms:.1f}%")

    # Generate plot
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        # Compute full Dl for best params (using l_max=25 cache)
        Dl_best = compute_Dl_with_params(best_final_params[0], best_final_params[1],
                                          best_final_params[2])

        fig, axes = plt.subplots(2, 1, figsize=(14, 10),
                                  gridspec_kw={'height_ratios': [3, 1]}, sharex=True)
        fig.subplots_adjust(hspace=0.06)

        ax = axes[0]
        ax.plot(class_ell, class_Dl, 'k-', lw=1.2, alpha=0.7, label='CLASS')
        ax.plot(ell_scan_sparse, Dl_best, 'r-', lw=1.5, alpha=0.9,
                label=f'Optimized (RMS={best_final_rms:.1f}%)')
        ax.set_ylabel(r'$\mathcal{D}_\ell^{TT}$ [$\mu$K$^2$]', fontsize=14)
        ax.set_title(f'Optimized D(k): D0={best_final_params[0]:.2f}, '
                     f'gamma={best_final_params[1]:.2f}, kDf={best_final_params[2]:.2f}',
                     fontsize=12)
        ax.legend(fontsize=10)
        ax.set_xlim(30, 2500)
        ax.grid(alpha=0.15)

        ax2 = axes[1]
        mask_plot = (class_at_scan > 100) & (ell_scan_sparse > 30)
        residual = np.where(mask_plot, (Dl_best - class_at_scan) / class_at_scan * 100, 0)
        ax2.plot(ell_scan_sparse, residual, 'r-', lw=1.0)
        ax2.axhline(0, color='k', ls='--', lw=0.8)
        ax2.fill_between(ell_scan_sparse, -10, 10, color='green', alpha=0.08)
        ax2.set_ylabel('(opt - CLASS)/CLASS [%]', fontsize=11)
        ax2.set_xlabel(r'$\ell$', fontsize=14)
        ax2.set_ylim(-60, 30)
        ax2.grid(alpha=0.15)
        ax2.text(0.02, 0.90, f'RMS = {best_final_rms:.1f}%',
                 transform=ax2.transAxes, fontsize=10, va='top',
                 bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

        out_dir = os.path.dirname(os.path.abspath(__file__))
        plot_path = os.path.join(out_dir, 'dk_optimization.png')
        plt.savefig(plot_path, dpi=150, bbox_inches='tight')
        print(f"\n  Plot saved: {plot_path}")
        plt.close()

        # Also plot D(k) comparison
        fig2, ax2 = plt.subplots(figsize=(10, 6))
        k_plot = np.linspace(0.001, 0.3, 1000)

        # Current D(k)
        from .radiation_driving import amplitude_correction
        D_current = amplitude_correction(k_plot)
        x_eq_p = k_plot / k_eq
        f_t = x_eq_p**2 / (1 + x_eq_p**2)
        D_current_eff = 1.0 + (D_current - 1.0) * f_t

        # Optimized D(k)
        D_opt_raw = best_final_params[0] * np.exp(best_final_params[1] * k_plot**2)
        D_opt_eff = 1.0 + (D_opt_raw - 1.0) * f_t

        ax2.plot(k_plot, D_current_eff, 'b-', lw=1.5, label='Current D(k)')
        ax2.plot(k_plot, D_opt_eff, 'r-', lw=1.5, label='Optimized D(k)')

        # Plot measured D from hierarchy
        if 'D_at_peaks' in results_by_lmax[25]:
            k_pts = results_by_lmax[25]['k_D_pts']
            D_pts = results_by_lmax[25]['D_at_peaks']
            mask = (k_pts > 0.01) & (k_pts < 0.30) & np.isfinite(D_pts) & (D_pts > 0.5) & (D_pts < 20)
            ax2.scatter(k_pts[mask], D_pts[mask], c='green', s=40, zorder=5,
                       label='Measured (hier/exact)')

        ax2.set_xlabel('k [Mpc$^{-1}$]', fontsize=14)
        ax2.set_ylabel('D(k)', fontsize=14)
        ax2.set_title('Amplitude correction D(k)', fontsize=13)
        ax2.legend(fontsize=10)
        ax2.grid(alpha=0.15)
        ax2.set_xlim(0, 0.3)

        plot_path2 = os.path.join(out_dir, 'dk_comparison.png')
        plt.savefig(plot_path2, dpi=150, bbox_inches='tight')
        print(f"  D(k) plot saved: {plot_path2}")
        plt.close()

    except Exception as e:
        print(f"  Plot error: {e}")

    return {
        'rms_current': rms_current,
        'rms_3param': best_rms_3param,
        'params_3param': best_3param,
        'rms_5param': res5.fun,
        'params_5param': res5.x,
        'rms_lmax50': res50.fun if 50 in results_by_lmax else None,
        'params_lmax50': res50.x if 50 in results_by_lmax else None,
        'speed_lmax25': results_by_lmax[25]['t_solve'],
        'speed_lmax50': results_by_lmax[50]['t_solve'] if 50 in results_by_lmax else None,
    }


if __name__ == '__main__':
    run_optimization()
