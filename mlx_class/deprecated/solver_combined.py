# DEPRECATED: Use solver_production.py instead.
# This file is kept for backward compatibility only.
import warnings as _warnings
_warnings.warn(
    "mlx_class.solver_combined is deprecated. Use mlx_class.solver_production instead.",
    DeprecationWarning,
    stacklevel=2,
)

#!/usr/bin/env python3
"""
solver_combined.py -- Ultimate combined CMB solver pipeline.

Merges ALL improvements into one pipeline:
  1. Peebles recombination     -> correct visibility peak at z~1090 (not z~1043)
  2. Neutrino IMEX solver      -> Psi != Phi, N_eff=3.046, l_nu_max=20
  3. ODE source extraction     -> correct amplitude/phase from solving equations
  4. LOS projection            -> visibility-weighted Bessel with Doppler
  5. Late ISW                  -> dark energy contribution at low l

Key innovation:
  ODE source (correct amplitude/phase from solving Boltzmann equations)
  + LOS projection (correct peak shift from visibility weighting)
  + Peebles (correct visibility shape)
  + Neutrinos (correct radiation content, Psi != Phi)

STRATEGY for the Bessel integral:
  The visibility function is narrow (sigma ~ 16 Mpc << tau_rec ~ 281 Mpc).
  The source S(k,tau) varies slowly across this window.
  We use a hybrid approach:
    - Pre-compute visibility-weighted Bessel: <g*j_l>(k) and <g*j_l'>(k)
    - The source S_SW(k) and S_Dop(k) are evaluated at the visibility peak
      (where g*S is maximized), then projected with the pre-computed Bessel weights.
    - This is equivalent to the "frozen source" LOS integral used in spectra_los.py
      but with the ODE-derived source instead of analytic Hu-Sugiyama.

Author: Sheng-Kai Huang, 2026
"""
import os
import sys
import time
import numpy as np

# Ensure package import works
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import mlx.core as mx
from scipy.special import spherical_jn
from scipy.interpolate import interp1d
from scipy.signal import find_peaks
from scipy.ndimage import gaussian_filter1d

from mlx_class.background import (
    Background, A_s, n_s, k_pivot, T_CMB, k_eq,
    Omega_r as _OMEGA_R, Omega_b as _OMEGA_B,
    Omega_c as _OMEGA_C, H0_Mpc as _H0_MPC,
    a_eq as _A_EQ, Omega_m as _OMEGA_M,
    Omega_L as _OMEGA_L,
)
from mlx_class.perturbations_neutrino import (
    NeutrinoBoltzmannSolver, NeutrinoResult, L_NU_MAX,
    IDX_PHI, IDX_THETA_0, IDX_THETA_1, IDX_N_START,
    n_var_total, adiabatic_ic, imex_rk4_step,
    _OMEGA_GAMMA, _OMEGA_NU, _f_nu,
)


# ============================================================================
# Planck reference peaks (TT power spectrum)
# ============================================================================
PLANCK_PEAKS = np.array([220, 540, 810, 1120, 1420])
PLANCK_HEIGHTS = np.array([5780, 2530, 2560, 1490, 840])  # approximate D_l in uK^2
PLANCK_RATIO_12 = 2.5  # approximate 1st/2nd peak height ratio


# ============================================================================
# Visibility function quadrature
# ============================================================================

def _build_visibility_quadrature(bg, N_tau=25):
    """Build quadrature points around the visibility peak."""
    peak_idx = np.argmax(bg.visibility_grid)
    tau_peak = bg.tau_grid[peak_idx]
    vis_peak = bg.visibility_grid[peak_idx]

    # Find FWHM
    half_max = vis_peak / 2.0
    left_idx = peak_idx
    while left_idx > 0 and bg.visibility_grid[left_idx] > half_max:
        left_idx -= 1
    right_idx = peak_idx
    while right_idx < len(bg.visibility_grid) - 1 and bg.visibility_grid[right_idx] > half_max:
        right_idx += 1

    fwhm = bg.tau_grid[right_idx] - bg.tau_grid[left_idx]
    sigma = fwhm / 2.355

    tau_lo = max(tau_peak - 4.0 * sigma, bg.tau_grid[1])
    tau_hi = min(tau_peak + 4.0 * sigma, bg.tau_grid[-2])
    tau_vis = np.linspace(tau_lo, tau_hi, N_tau)
    g_vis = bg.visibility_at_tau(tau_vis)

    return tau_vis, g_vis, tau_peak, sigma


# ============================================================================
# Late ISW
# ============================================================================

def _eisenstein_hu_transfer(k_arr):
    """Eisenstein-Hu (1998) zero-baryon transfer function."""
    h = 0.6736
    Om = _OMEGA_M
    Ob = _OMEGA_B
    Gamma = Om * h * np.exp(-Ob * (1.0 + np.sqrt(2.0 * h) / Om))
    q = k_arr / (Gamma * h)
    L = np.log(2.0 * np.e + 1.8 * q)
    C = 14.2 + 731.0 / (1.0 + 62.5 * q)
    return L / (L + C * q ** 2)


def _late_isw_phi_dot(k_arr, tau_arr, bg):
    """Compute Phi_dot(k,tau) in the dark energy era."""
    N_tau = len(tau_arr)
    N_k = len(k_arr)

    a_arr = bg.a_at_tau(tau_arr)
    calH_arr = bg.calH_at_tau(tau_arr)
    E_arr = calH_arr / (a_arr * _H0_MPC)

    T_k = _eisenstein_hu_transfer(k_arr)
    Phi_plateau = 0.9 * T_k

    Phi_dot = np.zeros((N_tau, N_k), dtype=np.float64)
    for it in range(N_tau):
        a = a_arr[it]
        calH = calH_arr[it]
        E = E_arr[it]
        Omega_m_a = _OMEGA_M / (a ** 3 * E ** 2)
        f_growth = Omega_m_a ** 0.55
        decay_rate = calH * (f_growth - 1.0)
        Phi_dot[it] = Phi_plateau * decay_rate

    return Phi_dot


# ============================================================================
# Neutrino solver to tau_rec
# ============================================================================

def solve_neutrino_to_rec(bg, k_arr_np, l_nu_max=L_NU_MAX):
    """
    Run the standard NeutrinoBoltzmannSolver to tau_rec and extract source.

    Returns (S_SW, S_Dop) as 1D arrays (N_k,) at tau_rec with Silk damping.
    """
    solver = NeutrinoBoltzmannSolver(bg, k_arr_np, l_nu_max=l_nu_max)
    result = solver.solve()
    Theta_0, Phi, Psi, v_b, N_0, N_2 = result.source_at_recombination()

    # SW source: Theta_0 + Psi (with neutrino anisotropic stress)
    S_SW = Theta_0 + Psi
    # Doppler source: v_b = 3*Theta_1 (already computed in source_at_recombination)
    S_Dop = v_b

    return S_SW, S_Dop, result


# ============================================================================
# Upsample source from coarse ODE k-grid to fine k-grid
# ============================================================================

def upsample_source(k_coarse, source_coarse, D_A, ell_max=2500):
    """
    Upsample the smooth source function S(k) from the coarse ODE k-grid
    to a fine uniform k-grid that resolves j_l(k*chi) oscillations.

    j_l oscillation period in k: pi/chi ~ pi/D_A.
    We need dk < pi/(2*D_A) for Nyquist, use ~5 points per oscillation.
    """
    dk_nyquist = np.pi / (2.0 * D_A)
    dk_fine = dk_nyquist / 1.5  # ~3 points per Bessel oscillation

    k_min = float(k_coarse[0])
    k_max_needed = (ell_max + 500.0) / D_A
    k_max = min(float(k_coarse[-1]), k_max_needed)

    k_fine = np.arange(k_min, k_max, dk_fine).astype(np.float64)

    f_source = interp1d(k_coarse, source_coarse, kind='cubic',
                        fill_value=0.0, bounds_error=False)
    source_fine = f_source(k_fine)

    return k_fine, source_fine


# ============================================================================
# Main combined solver pipeline
# ============================================================================

def run_combined_pipeline(N_k=500, ell_max=2500, N_tau_vis=25, N_tau_late_isw=10,
                          l_nu_max=L_NU_MAX):
    """
    Ultimate combined CMB C_l pipeline.

    Steps:
      1. Background with Peebles recombination
      2. Neutrino IMEX solver to tau_rec: get S_SW(k), S_Dop(k)
      3. Upsample smooth ODE source to fine k-grid
      4. Compute visibility-weighted Bessel tables (once)
      5. LOS projection: Delta_l(k) = S_SW*<g*j_l> + S_Dop*<g*j_l'>
      6. Add late ISW
      7. C_l integration

    The key insight: the source S(k) is FROZEN at tau_rec (like spectra_los.py
    analytic mode), but the source comes from the ODE solver rather than an
    analytic formula. The visibility weighting is in the Bessel tables.

    Returns: dict with ell, Cl, Dl, peaks, timing info
    """
    t_start = time.time()

    print("=" * 72)
    print("  ULTIMATE COMBINED SOLVER")
    print("  Peebles + Neutrinos + ODE source + LOS projection")
    print("=" * 72)

    # ================================================================
    # Step 1: Background with Peebles recombination
    # ================================================================
    print("\n--- Step 1: Background (Peebles recombination) ---")
    t0 = time.time()
    bg = Background(khronon=False, recombination='peebles')
    bg.solve()
    t_bg = time.time() - t0

    # Report visibility peak
    peak_idx_bg = np.argmax(bg.visibility_grid)
    z_vis_peak = 1.0 / bg.a_grid[peak_idx_bg] - 1
    print(f"[Combined] Background: {t_bg:.2f}s")
    print(f"[Combined] Visibility peak: z = {z_vis_peak:.0f} "
          f"(Peebles: should be ~1090)")
    print(f"[Combined] tau_rec = {bg.tau_rec:.1f} Mpc, "
          f"D_A = {bg.D_A:.1f} Mpc, r_s = {bg.r_s:.1f} Mpc")

    # ================================================================
    # Step 2: Solve neutrino IMEX to tau_rec
    # ================================================================
    print("\n--- Step 2: Neutrino IMEX solver ---")
    k_arr = np.geomspace(5e-4, 0.35, N_k).astype(np.float32)
    N_k_actual = len(k_arr)

    t0 = time.time()
    S_SW_coarse, S_Dop_coarse, nu_result = solve_neutrino_to_rec(
        bg, k_arr, l_nu_max=l_nu_max)
    t_solve = time.time() - t0

    print(f"[Combined] SW source range: "
          f"[{np.min(S_SW_coarse):.4f}, {np.max(S_SW_coarse):.4f}]")
    print(f"[Combined] Doppler source range: "
          f"[{np.min(S_Dop_coarse):.4f}, {np.max(S_Dop_coarse):.4f}]")
    print(f"[Combined] Psi/Phi (superhorizon): "
          f"{np.mean(nu_result.phi_psi_ratio()[:10]):.4f}")

    # ================================================================
    # Step 3: Upsample to fine k-grid
    # ================================================================
    print("\n--- Step 3: Upsample source ---")
    t0 = time.time()
    k_fine, S_SW_fine = upsample_source(k_arr, S_SW_coarse, bg.D_A, ell_max)
    _, S_Dop_fine = upsample_source(k_arr, S_Dop_coarse, bg.D_A, ell_max)
    N_kf = len(k_fine)
    t_upsample = time.time() - t0
    print(f"[Combined] Upsampled: {N_k_actual} -> {N_kf} k-points "
          f"(dk = {k_fine[1]-k_fine[0]:.2e} Mpc^-1)")

    # Primordial power spectrum on fine grid
    P_R_fine = A_s * (k_fine / k_pivot) ** (n_s - 1.0)
    lnk_fine = np.log(k_fine)
    dlnk_fine = np.diff(lnk_fine)

    # ================================================================
    # Step 4: Visibility quadrature + Bessel tables
    # ================================================================
    print("\n--- Step 4: Visibility-weighted Bessel tables ---")
    tau_vis, g_vis, tau_peak, sigma_vis = _build_visibility_quadrature(
        bg, N_tau_vis)
    chi_vis = bg.tau_0 - tau_vis
    dtau_vis = tau_vis[1] - tau_vis[0]

    print(f"[Combined] Visibility: {N_tau_vis} points around "
          f"tau = {tau_peak:.1f} Mpc (sigma = {sigma_vis:.1f} Mpc)")

    ell_values = np.unique(np.concatenate([
        np.arange(2, 30, 1),
        np.arange(30, 100, 2),
        np.arange(100, 500, 4),
        np.arange(500, 1500, 8),
        np.arange(1500, min(2501, ell_max + 1), 12),
    ])).astype(int)
    ell_values = ell_values[ell_values <= ell_max]
    N_ell = len(ell_values)

    # Pre-compute visibility-weighted Bessel: <g*j_l>(k) and <g*j_l'>(k)
    t0 = time.time()
    w_all = g_vis * dtau_vis  # (N_tau_vis,)

    gjl = np.zeros((N_ell, N_kf), dtype=np.float64)
    gjlp = np.zeros((N_ell, N_kf), dtype=np.float64)

    print(f"[Combined] Computing <g*j_l> ({N_ell} ells x {N_kf} k "
          f"x {N_tau_vis} tau)...")

    for il, ell in enumerate(ell_values):
        ell_int = int(ell)
        for it in range(N_tau_vis):
            x = k_fine * chi_vis[it]
            jl_val = spherical_jn(ell_int, x)
            jlp_val = spherical_jn(ell_int, x, derivative=True)
            gjl[il] += w_all[it] * jl_val
            gjlp[il] += w_all[it] * jlp_val

        if il % 50 == 0 and il > 0:
            print(f"  ... ell = {ell} ({il}/{N_ell})")

    t_bessel = time.time() - t0
    print(f"[Combined] Bessel tables: {t_bessel:.2f}s")

    # ================================================================
    # Step 5: Transfer function Delta_l(k) = S_SW * <g*j_l> + S_Dop * <g*j_l'>
    # ================================================================
    t0 = time.time()

    Delta_l_sw = np.zeros((N_ell, N_kf), dtype=np.float64)
    Delta_l_dop = np.zeros((N_ell, N_kf), dtype=np.float64)

    for il in range(N_ell):
        Delta_l_sw[il] = S_SW_fine * gjl[il]
        Delta_l_dop[il] = S_Dop_fine * gjlp[il]

    Delta_l = Delta_l_sw + Delta_l_dop

    t_transfer = time.time() - t0
    print(f"[Combined] Transfer: {t_transfer:.4f}s")

    # ================================================================
    # Step 6: Late ISW contribution (on coarser k-grid for speed)
    # ================================================================
    t0 = time.time()
    tau_late_lo = tau_vis[-1]
    tau_late_hi = bg.tau_0 * 0.98
    tau_late = np.linspace(tau_late_lo, tau_late_hi, N_tau_late_isw)
    chi_late = bg.tau_0 - tau_late
    kappa_late = bg.kappa_at_tau(tau_late)
    exp_neg_kappa_late = np.exp(-kappa_late)

    # Late ISW uses coarser k-grid (ISW is a low-l effect, doesn't need fine k)
    k_late = np.geomspace(float(k_fine[0]), float(k_fine[-1]),
                          min(500, N_kf)).astype(np.float64)
    late_Phi_dot = _late_isw_phi_dot(k_late, tau_late, bg)

    Delta_l_lisw_coarse = np.zeros((N_ell, len(k_late)), dtype=np.float64)
    if N_tau_late_isw > 1:
        dtau_late = np.diff(tau_late)
        for it in range(N_tau_late_isw - 1):
            integrand_lo = exp_neg_kappa_late[it] * 2.0 * late_Phi_dot[it]
            integrand_hi = exp_neg_kappa_late[it + 1] * 2.0 * late_Phi_dot[it + 1]
            integrand_avg = 0.5 * (integrand_lo + integrand_hi)
            dt = dtau_late[it]
            for il, ell in enumerate(ell_values):
                ell_int = int(ell)
                jl_lo = spherical_jn(ell_int, k_late * chi_late[it])
                jl_hi = spherical_jn(ell_int, k_late * chi_late[it + 1])
                jl_avg = 0.5 * (jl_lo + jl_hi)
                Delta_l_lisw_coarse[il] += integrand_avg * jl_avg * dt

    # Interpolate late ISW Delta_l onto fine k-grid
    Delta_l_lisw = np.zeros((N_ell, N_kf), dtype=np.float64)
    for il in range(N_ell):
        f_lisw = interp1d(k_late, Delta_l_lisw_coarse[il], kind='cubic',
                          fill_value=0.0, bounds_error=False)
        Delta_l_lisw[il] = f_lisw(k_fine)

    Delta_l_total = Delta_l + Delta_l_lisw
    t_lisw = time.time() - t0
    print(f"[Combined] Late ISW: {t_lisw:.2f}s")

    # ================================================================
    # Step 7: C_l = 4*pi * int dk/k P_R |Delta_l|^2
    # ================================================================
    t0 = time.time()
    dlnk_f32 = dlnk_fine.astype(np.float32)

    def _compute_Dl(Delta):
        integrand = (P_R_fine[None, :] * Delta ** 2).astype(np.float32)
        integrand_gpu = mx.array(integrand)
        dlnk_gpu = mx.array(dlnk_f32)
        mid = 0.5 * (integrand_gpu[:, :-1] + integrand_gpu[:, 1:])
        Cl_gpu = 4.0 * np.pi * mx.sum(mid * dlnk_gpu[None, :], axis=1)
        mx.eval(Cl_gpu)
        Cl_np = np.maximum(np.array(Cl_gpu), 0.0)
        ell_f = ell_values.astype(float)
        return Cl_np, ell_f * (ell_f + 1.0) * Cl_np / (2.0 * np.pi) * (T_CMB * 1e6) ** 2

    Cl, Dl = _compute_Dl(Delta_l_total)
    _, Dl_sw_only = _compute_Dl(Delta_l_sw)
    _, Dl_sw_dop = _compute_Dl(Delta_l)  # SW + Doppler, no late ISW

    t_cl = time.time() - t0
    print(f"[Combined] C_l integration: {t_cl:.4f}s")

    # ================================================================
    # Step 8: Find peaks
    # ================================================================
    t_total = time.time() - t_start

    print(f"\n[Combined] Total time: {t_total:.2f}s")
    print(f"  Background:  {t_bg:.2f}s")
    print(f"  ODE solver:  {t_solve:.2f}s")
    print(f"  Upsample:    {t_upsample:.2f}s")
    print(f"  Bessel proj: {t_bessel:.2f}s")
    print(f"  Late ISW:    {t_lisw:.2f}s")

    # Find peaks in smoothed D_l
    # Only search above l=100 to skip ISW bump at low ell
    l_full = np.arange(2, ell_max + 1)
    f_interp = interp1d(ell_values, Dl, kind='cubic', fill_value='extrapolate')
    Dl_full = np.maximum(f_interp(l_full), 0.0)

    Dl_smooth = gaussian_filter1d(Dl_full, sigma=15)

    # Search for acoustic peaks only above l=100
    l_search_mask = l_full >= 100
    l_search = l_full[l_search_mask]
    Dl_search = Dl_smooth[l_search_mask]
    peaks_idx_rel, props = find_peaks(Dl_search, distance=120, prominence=20)

    peak_ells = l_search[peaks_idx_rel]
    peak_heights = Dl_search[peaks_idx_rel]

    return {
        'ell_values': ell_values,
        'Cl': Cl,
        'Dl': Dl,
        'Dl_sw_only': Dl_sw_only,
        'Dl_sw_dop': Dl_sw_dop,
        'l_full': l_full,
        'Dl_full': Dl_full,
        'Dl_smooth': Dl_smooth,
        'peak_ells': peak_ells,
        'peak_heights': peak_heights,
        'bg': bg,
        'k_arr': k_arr,
        'k_fine': k_fine,
        'tau_vis': tau_vis,
        'g_vis': g_vis,
        'tau_peak': tau_peak,
        'sigma_vis': sigma_vis,
        'z_vis_peak': z_vis_peak,
        't_total': t_total,
        't_bg': t_bg,
        't_solve': t_solve,
        't_upsample': t_upsample,
        't_bessel': t_bessel,
        't_lisw': t_lisw,
    }


# ============================================================================
# Peak analysis and reporting
# ============================================================================

def print_peak_report(result):
    """Print detailed peak comparison with Planck."""
    peak_ells = result['peak_ells']
    peak_heights = result['peak_heights']

    print("\n" + "=" * 72)
    print("  PEAK POSITIONS vs PLANCK")
    print("=" * 72)
    print(f"  {'Peak':<6} {'Planck':>8} {'Found':>8} {'Error':>10} "
          f"{'Height':>12} {'Planck H':>12}")
    print(f"  {'-'*58}")

    n = min(len(peak_ells), len(PLANCK_PEAKS))
    total_err = 0.0
    for i in range(n):
        err = (peak_ells[i] - PLANCK_PEAKS[i]) / PLANCK_PEAKS[i] * 100
        total_err += abs(err)
        print(f"  {i+1:<6} {PLANCK_PEAKS[i]:>8} {peak_ells[i]:>8} "
              f"{err:>+9.1f}% {peak_heights[i]:>10.0f} "
              f"{PLANCK_HEIGHTS[i]:>10.0f} uK^2")

    print(f"  {'-'*58}")
    if n > 0:
        print(f"  Mean absolute error: {total_err/n:.1f}%")
    if n >= 2:
        ratio = peak_heights[0] / peak_heights[1]
        print(f"  1st/2nd peak ratio: {ratio:.3f}  (Planck ~ {PLANCK_RATIO_12})")
    if n >= 3:
        ratio_13 = peak_heights[0] / peak_heights[2]
        print(f"  1st/3rd peak ratio: {ratio_13:.3f}  (Planck ~ {5780/2560:.2f})")

    print(f"\n  Visibility peak: z = {result['z_vis_peak']:.0f}")
    print(f"  Total time: {result['t_total']:.2f}s")
    print("=" * 72)


# ============================================================================
# Plot generation
# ============================================================================

def generate_plot(result, out_path):
    """Generate the combined_all.png plot."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    ell_values = result['ell_values']
    Dl = result['Dl']
    l_full = result['l_full']
    Dl_full = result['Dl_full']
    Dl_smooth = result['Dl_smooth']
    Dl_sw_only = result['Dl_sw_only']
    Dl_sw_dop = result['Dl_sw_dop']
    peak_ells = result['peak_ells']
    peak_heights = result['peak_heights']
    bg = result['bg']
    tau_vis = result['tau_vis']
    g_vis = result['g_vis']
    t_total = result['t_total']
    z_vis_peak = result['z_vis_peak']

    fig, axes = plt.subplots(2, 2, figsize=(16, 12))

    # --- Panel 1: Full D_l spectrum ---
    ax = axes[0, 0]
    ax.plot(l_full, Dl_full, 'b-', lw=0.6, alpha=0.4, label='Raw $D_\\ell$')
    ax.plot(l_full, Dl_smooth, 'b-', lw=1.8, label='Smoothed $D_\\ell$')

    # Mark found peaks
    colors = ['red', 'orange', 'green', 'purple', 'brown', 'pink', 'cyan']
    n_peaks_show = min(7, len(peak_ells))
    for i in range(n_peaks_show):
        c = colors[i % len(colors)]
        ax.plot(peak_ells[i], peak_heights[i], 'o', color=c, ms=9, zorder=5)
        label_text = f'l={peak_ells[i]}'
        if i < len(PLANCK_PEAKS):
            err = (peak_ells[i] - PLANCK_PEAKS[i]) / PLANCK_PEAKS[i] * 100
            label_text += f' ({err:+.0f}%)'
        ax.annotate(label_text, (peak_ells[i], peak_heights[i] * 1.06),
                    fontsize=8, ha='center', color=c, fontweight='bold')

    # Planck reference lines
    for lp in PLANCK_PEAKS:
        ax.axvline(lp, color='green', ls=':', alpha=0.3, lw=1)
    ax.axvline(PLANCK_PEAKS[0], color='green', ls=':', alpha=0.3, lw=1,
               label='Planck peaks')

    ax.set_xlabel('Multipole $\\ell$', fontsize=12)
    ax.set_ylabel('$D_\\ell$ [$\\mu K^2$]', fontsize=12)
    ax.set_title(f'CMB TT (Peebles + Neutrinos + ODE + LOS)', fontsize=12)
    ax.legend(fontsize=9, loc='upper right')
    ax.set_xlim(2, 2500)
    if np.max(Dl_smooth) > 0:
        ax.set_ylim(0, np.max(Dl_smooth) * 1.35)

    # --- Panel 2: Component decomposition ---
    ax = axes[0, 1]
    f_sw = interp1d(ell_values, Dl_sw_only, kind='cubic', fill_value='extrapolate')
    f_swdop = interp1d(ell_values, Dl_sw_dop, kind='cubic', fill_value='extrapolate')
    Dl_sw_full = np.maximum(f_sw(l_full), 0.0)
    Dl_swdop_full = np.maximum(f_swdop(l_full), 0.0)

    Dl_sw_sm = gaussian_filter1d(Dl_sw_full, sigma=15)
    Dl_swdop_sm = gaussian_filter1d(Dl_swdop_full, sigma=15)

    ax.plot(l_full, Dl_smooth, 'b-', lw=1.8, label='Total (SW+Dop+ISW)')
    ax.plot(l_full, Dl_swdop_sm, 'r--', lw=1.2, alpha=0.7, label='SW + Doppler')
    ax.plot(l_full, Dl_sw_sm, 'g:', lw=1.5, alpha=0.8, label='SW only')

    ax.set_xlabel('Multipole $\\ell$', fontsize=12)
    ax.set_ylabel('$D_\\ell$ [$\\mu K^2$]', fontsize=12)
    ax.set_title('Component Decomposition', fontsize=12)
    ax.legend(fontsize=9)
    ax.set_xlim(2, 1500)
    if np.max(Dl_smooth[:1498]) > 0:
        ax.set_ylim(0, np.max(Dl_smooth[:1498]) * 1.35)

    # --- Panel 3: Visibility function ---
    ax = axes[1, 0]
    z_grid = 1.0 / bg.a_grid - 1.0
    mask_vis = (z_grid > 800) & (z_grid < 1500)
    ax.plot(z_grid[mask_vis], bg.visibility_grid[mask_vis], 'b-', lw=1.5,
            label='$g(\\tau)$ (Peebles)')
    z_vis_tau = 1.0 / bg.a_at_tau(tau_vis) - 1.0
    ax.plot(z_vis_tau, g_vis, 'ro', ms=4, alpha=0.6,
            label=f'Quadrature ({len(tau_vis)} pts)')
    ax.axvline(z_vis_peak, color='orange', ls='--', alpha=0.7,
               label=f'Peak z={z_vis_peak:.0f}')
    ax.axvline(1090, color='green', ls=':', alpha=0.5,
               label='Planck $z_*$=1090')
    ax.set_xlabel('Redshift $z$', fontsize=12)
    ax.set_ylabel('Visibility $g(\\tau)$', fontsize=12)
    ax.set_title('Visibility Function (Peebles Recombination)', fontsize=12)
    ax.legend(fontsize=9)
    ax.invert_xaxis()

    # --- Panel 4: Peak comparison bar chart ---
    ax = axes[1, 1]
    n_plot = min(5, len(peak_ells))
    if n_plot >= 1:
        x = np.arange(n_plot)
        width = 0.35
        bars_planck = ax.bar(x - width/2, PLANCK_PEAKS[:n_plot], width,
                             color='green', alpha=0.6, label='Planck')
        bars_found = ax.bar(x + width/2, peak_ells[:n_plot], width,
                            color='blue', alpha=0.6, label='Combined solver')

        for i in range(n_plot):
            err = (peak_ells[i] - PLANCK_PEAKS[i]) / PLANCK_PEAKS[i] * 100
            ax.annotate(f'{err:+.1f}%',
                        (x[i] + width/2, peak_ells[i] + 15),
                        fontsize=9, ha='center', color='blue', fontweight='bold')

        ax.set_xticks(x)
        ax.set_xticklabels([f'Peak {i+1}' for i in range(n_plot)])
        ax.set_ylabel('Multipole $\\ell$', fontsize=12)
        ax.set_title('Peak Positions vs Planck', fontsize=12)
        ax.legend(fontsize=10)

    plt.suptitle(
        f'Ultimate Combined Solver: '
        f'Peebles + Neutrinos ($N_{{eff}}$=3.046) + ODE + LOS  '
        f'({t_total:.1f}s)',
        fontsize=14, y=1.01)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\nPlot saved: {out_path}")


# ============================================================================
# Main
# ============================================================================

def main():
    t_start = time.time()

    # Run combined pipeline
    result = run_combined_pipeline(
        N_k=500,
        ell_max=2500,
        N_tau_vis=15,
        N_tau_late_isw=10,
        l_nu_max=L_NU_MAX,
    )

    # Print peak report
    print_peak_report(result)

    # Generate plot
    out_dir = os.path.dirname(os.path.abspath(__file__))
    out_path = os.path.join(out_dir, 'combined_all.png')
    generate_plot(result, out_path)

    # Final summary
    peak_ells = result['peak_ells']
    peak_heights = result['peak_heights']
    t_total = result['t_total']

    print("\n" + "=" * 72)
    print("  FINAL SUMMARY")
    print("=" * 72)
    print(f"  Physics included:")
    print(f"    [x] Peebles 3-level recombination (vis peak z={result['z_vis_peak']:.0f})")
    print(f"    [x] Neutrinos N_eff=3.046, l_nu_max={L_NU_MAX} (Psi != Phi)")
    print(f"    [x] ODE source (no WKB/analytic approximation)")
    print(f"    [x] LOS visibility-weighted Bessel projection")
    print(f"    [x] Doppler term (j_l' projection)")
    print(f"    [x] Late ISW (dark energy)")
    print(f"    [x] Silk damping")
    print(f"    [x] k-grid upsampled: {len(result['k_arr'])} -> {len(result['k_fine'])} "
          f"(dk = {result['k_fine'][1]-result['k_fine'][0]:.2e})")

    n = min(len(peak_ells), len(PLANCK_PEAKS))
    if n > 0:
        errors = [(peak_ells[i] - PLANCK_PEAKS[i]) / PLANCK_PEAKS[i] * 100
                  for i in range(n)]
        mean_err = np.mean(np.abs(errors))
        print(f"\n  Peak positions (vs Planck):")
        for i in range(n):
            status = "OK" if abs(errors[i]) < 15 else "WARN"
            print(f"    Peak {i+1}: l={peak_ells[i]:>5} "
                  f"(Planck {PLANCK_PEAKS[i]:>5}, err={errors[i]:+.1f}%) [{status}]")
        print(f"    Mean |error|: {mean_err:.1f}%")

    if n >= 2:
        ratio = peak_heights[0] / peak_heights[1]
        print(f"\n  Height ratios:")
        print(f"    1st/2nd: {ratio:.3f}  (Planck ~ {PLANCK_RATIO_12})")
    if n >= 3:
        ratio_13 = peak_heights[0] / peak_heights[2]
        print(f"    1st/3rd: {ratio_13:.3f}  (Planck ~ {5780/2560:.2f})")

    print(f"\n  Total time: {t_total:.2f}s")
    print("=" * 72)

    return result


if __name__ == '__main__':
    main()
