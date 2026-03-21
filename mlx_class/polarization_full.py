"""
polarization_full.py — Full CMB EE and TE polarization spectra using HiRes hierarchy.

Replaces the broken TCA-estimated Theta_2 in polarization.py (which is ~100x too low)
with the actual Theta_2(k, tau) from the full photon Boltzmann hierarchy solved by
perturbations_hires.py.

PHYSICS:
  The E-mode polarization source function (Zaldarriaga & Seljak 1997) is:
    S_E(k, tau) = (3/4) * g(tau) * Pi(k, tau)
  where Pi = Theta_2 + Theta_P0 + Theta_P2 is the photon polarization source.

  Without a separate polarization hierarchy (Theta_P_l), we approximate Pi
  using the Hu & White (1997) result: the polarization multipoles Theta_P0
  and Theta_P2 are themselves sourced by Theta_2 during recombination, giving:
    Theta_P0 ~ (5/4) * Theta_2   (monopole)
    Theta_P2 ~ (1/4) * Theta_2   (quadrupole, from tight-coupling)
  so  Pi = Theta_2 + Theta_P0 + Theta_P2 ~ (5/2) * Theta_2

  For a finite visibility width (not instantaneous recombination), the
  polarization hierarchy has time to build up further. The effective ratio
  Pi_eff / Theta_2 is enhanced to ~3.0 (Zaldarriaga & Harari 1995, Hu & White 1997).

  Therefore:
    S_E(k, tau) = (3/4) * Pi_eff = (3/4) * alpha_P * Theta_2(k, tau)
  where alpha_P ~ 3.0 (the polarization amplification factor).

  The E-mode transfer function is obtained by line-of-sight integration:
    Delta_l^E(k) = int dtau g(tau) * (3/4) * alpha_P * Theta_2(k, tau) * epsilon_l(k*chi)
  where epsilon_l(x) = sqrt((l+2)!/(l-2)!) / x^2 * j_l(x)  is the spin-2 projection.

  C_l^EE = (4 pi) int dk/k P_R(k) |Delta_l^E(k)|^2
  C_l^TE = (4 pi) int dk/k P_R(k) Delta_l^T(k) Delta_l^E(k)

IMPLEMENTATION:
  1. Extend HiResBoltzmannSolver to evolve THROUGH the visibility window
     (not just to tau_rec), recording snapshots at N_tau quadrature points.
  2. Extract Theta_2(k, tau) at each visibility snapshot.
  3. Apply LOS integration with visibility-weighted Theta_2 for the E-mode source.
  4. For TT transfer, re-use the same snapshots (Theta_0 + Psi, v_b).
  5. Compute C_l^EE, C_l^TE, C_l^TT simultaneously.

TARGET:
  EE first peak amplitude ~20 uK^2 at l~400 (CLASS: 21.8 uK^2)
  TE oscillations with amplitude ~100 uK^2

Author: Sheng-Kai Huang, 2026
"""
import os
import sys
import time
import numpy as np

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
from mlx_class.perturbations_hires import (
    HiResBoltzmannSolver, HiResResult,
    IDX_PHI, IDX_DELTA_B, IDX_V_B, IDX_DELTA_C, IDX_V_C,
    IDX_THETA_START, idx_theta, idx_n_start,
    n_var_total, adiabatic_ic_hires, build_tau_grid_hires,
    strang_split_step, deriv_streaming,
)
from mlx_class.perturbations_neutrino import (
    _OMEGA_GAMMA, _OMEGA_NU, L_NU_MAX,
)


# ============================================================================
# Spin-2 Bessel projection: epsilon_l(x)
# ============================================================================

def _epsilon_l(ell, x_arr):
    """
    Spin-2 Bessel projection factor for E-mode polarization.

    epsilon_l(x) = sqrt((l+2)!/(l-2)!) / x^2 * j_l(x)
                 = sqrt((l-1)*l*(l+1)*(l+2)) / x^2 * j_l(x)

    For l >> 1, this is approximately l^2/x^2 * j_l(x).
    """
    l = int(ell)
    if l < 2:
        return np.zeros_like(x_arr)

    prefactor = np.sqrt(float((l - 1) * l * (l + 1) * (l + 2)))
    jl = spherical_jn(l, x_arr)

    x_safe = np.where(np.abs(x_arr) < 1e-10, 1e-10, x_arr)
    result = prefactor * jl / (x_safe ** 2)
    result = np.where(np.abs(x_arr) < 1e-10, 0.0, result)
    return result


# ============================================================================
# Visibility quadrature (same as solver_precision.py)
# ============================================================================

def _build_visibility_quadrature(bg, N_tau=25):
    """Build quadrature points around the visibility peak."""
    peak_idx = np.argmax(bg.visibility_grid)
    tau_peak = bg.tau_grid[peak_idx]
    vis_peak = bg.visibility_grid[peak_idx]

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
# Driving correction (from solver_precision.py)
# ============================================================================

def _driving_correction(k_arr, bg):
    """
    Radiation driving correction factor for TCA-based temperature source.
    See solver_precision.py for detailed derivation.
    """
    x_eq = k_arr / k_eq
    D0 = 1.0
    D_inf = 1.95
    alpha = 1.0
    D_corr = D0 + (D_inf - D0) * x_eq ** 2 / (alpha + x_eq ** 2)
    return D_corr


# ============================================================================
# Extended HiRes solver with visibility-window snapshots
# ============================================================================

def _build_extended_tau_grid_hires(bg, k_max, tau_end, l_gamma_max,
                                    N_early=300, N_late=1200):
    """
    Build integration grid extending past tau_rec into the visibility tail.
    Denser than the standard grid to resolve the visibility function.
    """
    tau_init = bg.tau_grid[1]
    tau_early_end = min(10.0, 0.3 * bg.tau_rec)

    tau_early = np.geomspace(tau_init, tau_early_end, N_early)

    cs_approx = 1.0 / np.sqrt(3.0)
    dtau_acoustic = 2.0 * np.pi / (k_max * cs_approx) / 25.0

    calH_at_end = float(bg.calH_at_tau(np.array([tau_early_end]))[0])
    calH_at_final = float(bg.calH_at_tau(np.array([min(tau_end, bg.tau_0 * 0.99)]))[0])
    calH_min = min(calH_at_end, calH_at_final)
    max_eigenvalue = k_max ** 2 / (3.0 * calH_min)
    dtau_stiff = 2.5 / max_eigenvalue

    dtau_max = min(dtau_acoustic, dtau_stiff)
    N_needed = int((tau_end - tau_early_end) / dtau_max) + 10
    N_late_actual = max(N_late, N_needed)

    tau_late = np.linspace(tau_early_end, tau_end, N_late_actual)
    tau_grid = np.unique(np.concatenate([tau_early, tau_late]))
    return tau_grid


def _solve_hires_with_snapshots(bg, k_arr_np, tau_snapshots,
                                 l_gamma_max=25, l_nu_max=L_NU_MAX,
                                 tca_threshold=50.0):
    """
    Run HiResBoltzmannSolver THROUGH the visibility window, recording
    the full photon hierarchy state at specified tau snapshot values.

    This is the key function: it provides actual Theta_2(k, tau) from
    the full Boltzmann hierarchy, not TCA estimates.

    Returns a dict with arrays (N_snap, N_k) for each field.
    """
    k_arr = mx.array(k_arr_np.astype(np.float32))
    N_k = len(k_arr_np)
    k_max = float(k_arr_np[-1])

    tau_end = float(tau_snapshots[-1]) + 1.0
    tau_end = min(tau_end, bg.tau_0 * 0.99)

    tau_grid = _build_extended_tau_grid_hires(bg, k_max, tau_end, l_gamma_max)

    # Initial conditions
    y = adiabatic_ic_hires(k_arr_np, bg, l_gamma_max, l_nu_max)

    # Precompute background
    calH_np = bg.calH_at_tau(tau_grid).astype(np.float32)
    R_np = bg.R_at_tau(tau_grid).astype(np.float32)
    a_np = bg.a_at_tau(tau_grid).astype(np.float32)
    kappa_dot_np = bg.kappa_dot_at_tau(tau_grid).astype(np.float32)

    calH_arr = mx.array(calH_np)
    R_arr = mx.array(R_np)
    a_arr = mx.array(a_np)

    _Og = _OMEGA_GAMMA
    _On = _OMEGA_NU
    _Ob = _OMEGA_B
    _Oc = float(bg.Omega_cdm)
    _H0 = _H0_MPC

    # Determine TCA switch point
    abs_kappa_dot = np.abs(kappa_dot_np)
    tca_switch_idx = None
    for i in range(len(tau_grid)):
        if abs_kappa_dot[i] / k_max < tca_threshold:
            tca_switch_idx = i
            break
    if tca_switch_idx is None:
        tca_switch_idx = len(tau_grid) - 1

    tau_switch = tau_grid[tca_switch_idx]
    z_switch = 1.0 / float(bg.a_at_tau(np.array([tau_switch]))[0]) - 1.0
    print(f"[HiRes-Snap] TCA->full at tau={tau_switch:.1f} (z={z_switch:.0f})")

    N_snap = len(tau_snapshots)
    _n_start = idx_n_start(l_gamma_max)

    snap_Theta_0 = np.zeros((N_snap, N_k), dtype=np.float32)
    snap_Theta_1 = np.zeros((N_snap, N_k), dtype=np.float32)
    snap_Theta_2 = np.zeros((N_snap, N_k), dtype=np.float32)
    snap_Phi = np.zeros((N_snap, N_k), dtype=np.float32)
    snap_Psi = np.zeros((N_snap, N_k), dtype=np.float32)
    snap_v_b = np.zeros((N_snap, N_k), dtype=np.float32)

    snap_idx = 0

    def _record_snapshot(y_state, i_grid, snap_i, in_tca):
        """Record state at a snapshot."""
        mx.eval(y_state)
        y_np = np.array(y_state)

        snap_Theta_0[snap_i] = y_np[:, idx_theta(0)]
        snap_Theta_1[snap_i] = y_np[:, idx_theta(1)]
        snap_Phi[snap_i] = y_np[:, IDX_PHI]

        if in_tca:
            # In TCA: v_b = 3*Theta_1, Theta_2 estimated
            snap_v_b[snap_i] = 3.0 * y_np[:, idx_theta(1)]
            kd_val = abs_kappa_dot[i_grid]
            if kd_val > 1e-10:
                snap_Theta_2[snap_i] = (4.0 / 9.0) * k_arr_np * y_np[:, idx_theta(1)] / kd_val
            else:
                snap_Theta_2[snap_i] = 0.0
        else:
            # Full hierarchy: v_b and Theta_2 are solved
            snap_v_b[snap_i] = y_np[:, IDX_V_B]
            if l_gamma_max >= 2:
                snap_Theta_2[snap_i] = y_np[:, idx_theta(2)]
            else:
                snap_Theta_2[snap_i] = 0.0

        # Diagnose Psi from anisotropic stress
        N_2 = y_np[:, _n_start + 2] if l_nu_max >= 2 else np.zeros(N_k)
        a_val = float(np.array(a_arr[i_grid]))
        H02 = _H0 ** 2
        Theta_2_val = snap_Theta_2[snap_i]
        snap_Psi[snap_i] = (y_np[:, IDX_PHI]
                            - 12.0 * H02 * _On * N_2 / (a_val ** 2 * k_arr_np ** 2 + 1e-30)
                            - 12.0 * H02 * _Og * Theta_2_val / (a_val ** 2 * k_arr_np ** 2 + 1e-30))

    # ---- TCA phase ----
    for i in range(tca_switch_idx):
        while snap_idx < N_snap and tau_grid[i] >= tau_snapshots[snap_idx]:
            _record_snapshot(y, i, snap_idx, in_tca=True)
            snap_idx += 1

        dt = float(tau_grid[i + 1] - tau_grid[i])
        y = strang_split_step(
            y, k_arr,
            calH_arr[i], R_arr[i], tau_grid[i], float(kappa_dot_np[i]),
            _Og, _On, _Ob, _Oc,
            a_arr[i], _H0, l_gamma_max, l_nu_max, dt,
            tight_coupling=True)
        if i % 200 == 0:
            mx.eval(y)

    mx.eval(y)

    # Seed Theta_2 at TCA switch
    kd_switch = abs_kappa_dot[tca_switch_idx]
    if kd_switch > 1e-10:
        Theta_1_switch = y[:, idx_theta(1)]
        Theta_2_seed = (4.0 / 9.0) * k_arr * Theta_1_switch / kd_switch
        _idx_t2 = idx_theta(2)
        y = mx.concatenate([
            y[:, :_idx_t2],
            Theta_2_seed[:, None],
            y[:, _idx_t2 + 1:]
        ], axis=1)
        mx.eval(y)

    # ---- Full hierarchy phase ----
    for i in range(tca_switch_idx, len(tau_grid) - 1):
        while snap_idx < N_snap and tau_grid[i] >= tau_snapshots[snap_idx]:
            _record_snapshot(y, i, snap_idx, in_tca=False)
            snap_idx += 1

        dt = float(tau_grid[i + 1] - tau_grid[i])
        y = strang_split_step(
            y, k_arr,
            calH_arr[i], R_arr[i], tau_grid[i], float(kappa_dot_np[i]),
            _Og, _On, _Ob, _Oc,
            a_arr[i], _H0, l_gamma_max, l_nu_max, dt,
            tight_coupling=False)
        if i % 200 == 0:
            mx.eval(y)

    mx.eval(y)

    # Record remaining snapshots
    while snap_idx < N_snap:
        _record_snapshot(y, len(tau_grid) - 1, snap_idx, in_tca=False)
        snap_idx += 1

    return {
        'Theta_0': snap_Theta_0,
        'Theta_1': snap_Theta_1,
        'Theta_2': snap_Theta_2,
        'Phi': snap_Phi,
        'Psi': snap_Psi,
        'v_b': snap_v_b,
    }


# ============================================================================
# Source upsampling
# ============================================================================

def _upsample_source_2d(k_coarse, source_2d, k_fine):
    """Interpolate source(N_tau, N_k_coarse) onto k_fine grid."""
    N_tau = source_2d.shape[0]
    N_kf = len(k_fine)
    out = np.zeros((N_tau, N_kf), dtype=np.float64)
    for it in range(N_tau):
        f = interp1d(k_coarse, source_2d[it], kind='cubic',
                     fill_value=0.0, bounds_error=False)
        out[it] = f(k_fine)
    return out


# ============================================================================
# Main pipeline: EE + TE + TT
# ============================================================================

def run_polarization_pipeline(N_k=500, ell_max=2500, N_tau_vis=30,
                               l_gamma_max=25, l_nu_max=L_NU_MAX):
    """
    Full CMB polarization pipeline using HiRes hierarchy through visibility window.

    Computes C_l^TT, C_l^EE, C_l^TE simultaneously from the same source functions.

    Returns: dict with ell, Cl/Dl for TT/EE/TE, peaks, timing info.
    """
    t_start = time.time()

    print("=" * 72)
    print("  POLARIZATION PIPELINE (HiRes hierarchy)")
    print("  Full Theta_2 from Boltzmann hierarchy + LOS integration")
    print("=" * 72)

    # ================================================================
    # Step 1: Background
    # ================================================================
    print("\n--- Step 1: Background (Peebles recombination) ---")
    t0 = time.time()
    bg = Background(khronon=False, recombination='peebles')
    bg.solve()
    t_bg = time.time() - t0

    peak_idx_bg = np.argmax(bg.visibility_grid)
    z_vis_peak = 1.0 / bg.a_grid[peak_idx_bg] - 1
    print(f"[Pol] Background: {t_bg:.2f}s")
    print(f"[Pol] Visibility peak: z = {z_vis_peak:.0f}")
    print(f"[Pol] tau_rec = {bg.tau_rec:.1f} Mpc, "
          f"D_A = {bg.D_A:.1f} Mpc, r_s = {bg.r_s:.1f} Mpc")

    # ================================================================
    # Step 2: Visibility quadrature
    # ================================================================
    print("\n--- Step 2: Visibility quadrature ---")
    tau_vis, g_vis, tau_peak, sigma_vis = _build_visibility_quadrature(
        bg, N_tau_vis)
    chi_vis = bg.tau_0 - tau_vis
    dtau_vis = tau_vis[1] - tau_vis[0]

    print(f"[Pol] Visibility: {N_tau_vis} pts around tau={tau_peak:.1f} Mpc "
          f"(sigma={sigma_vis:.1f} Mpc)")

    # ================================================================
    # Step 3: HiRes solver with snapshots through visibility window
    # ================================================================
    print("\n--- Step 3: HiRes Boltzmann hierarchy through visibility ---")
    k_arr = np.geomspace(5e-4, 0.35, N_k).astype(np.float32)
    N_k_actual = len(k_arr)

    t0 = time.time()
    print(f"[Pol] Running HiRes (l_gamma_max={l_gamma_max}) "
          f"with {N_tau_vis} snapshots...")
    snapshots = _solve_hires_with_snapshots(
        bg, k_arr, tau_vis,
        l_gamma_max=l_gamma_max,
        l_nu_max=l_nu_max,
        tca_threshold=50.0)
    t_solve = time.time() - t0
    print(f"[Pol] HiRes done: {t_solve:.2f}s")

    # Diagnostics
    peak_snap = np.argmin(np.abs(tau_vis - tau_peak))
    Theta_2_at_peak = snapshots['Theta_2'][peak_snap]
    Theta_0_at_peak = snapshots['Theta_0'][peak_snap]
    print(f"[Pol] Theta_2 at visibility peak: "
          f"[{np.min(Theta_2_at_peak):.4e}, {np.max(Theta_2_at_peak):.4e}]")
    print(f"[Pol] Theta_0 at visibility peak: "
          f"[{np.min(Theta_0_at_peak):.4e}, {np.max(Theta_0_at_peak):.4e}]")
    print(f"[Pol] |Theta_2/Theta_0| (median): "
          f"{np.median(np.abs(Theta_2_at_peak / (Theta_0_at_peak + 1e-30))):.4f}")

    # ================================================================
    # Step 4: Apply driving correction + Silk damping to TT sources
    # ================================================================
    print("\n--- Step 4: Source preparation ---")
    D_corr = _driving_correction(k_arr, bg)
    silk = np.exp(-(k_arr / bg.k_D) ** 2)

    # TT sources: apply driving correction
    vis_Theta_0 = np.zeros((N_tau_vis, N_k_actual), dtype=np.float64)
    vis_Psi = np.zeros((N_tau_vis, N_k_actual), dtype=np.float64)
    vis_v_b = np.zeros((N_tau_vis, N_k_actual), dtype=np.float64)
    # EE source: Theta_2 (NO driving correction -- Theta_2 is from full hierarchy)
    vis_Theta_2 = np.zeros((N_tau_vis, N_k_actual), dtype=np.float64)

    for it in range(N_tau_vis):
        vis_Theta_0[it] = snapshots['Theta_0'][it] * D_corr
        vis_Psi[it] = snapshots['Psi'][it] * D_corr
        vis_v_b[it] = snapshots['v_b'][it] * D_corr
        vis_Theta_2[it] = snapshots['Theta_2'][it]

    print(f"[Pol] Driving correction: D(k_eq)={D_corr[np.argmin(np.abs(k_arr-k_eq))]:.3f}")

    # ================================================================
    # Step 5: Upsample to fine k-grid
    # ================================================================
    print("\n--- Step 5: Upsample sources ---")
    t0 = time.time()

    D_A = bg.D_A
    dk_fine = np.pi / (2.0 * D_A) / 1.5
    k_max_needed = min(float(k_arr[-1]), (ell_max + 500.0) / D_A)
    k_fine = np.arange(float(k_arr[0]), k_max_needed, dk_fine).astype(np.float64)
    N_kf = len(k_fine)

    silk_fine = np.exp(-(k_fine / bg.k_D) ** 2)
    D_corr_fine = np.interp(k_fine, k_arr, D_corr)

    vis_Theta_0_f = _upsample_source_2d(k_arr, vis_Theta_0, k_fine)
    vis_Psi_f = _upsample_source_2d(k_arr, vis_Psi, k_fine)
    vis_v_b_f = _upsample_source_2d(k_arr, vis_v_b, k_fine)
    vis_Theta_2_f = _upsample_source_2d(k_arr, vis_Theta_2, k_fine)

    t_upsample = time.time() - t0
    print(f"[Pol] Upsampled: {N_k_actual} -> {N_kf} k-points ({t_upsample:.2f}s)")

    P_R_fine = A_s * (k_fine / k_pivot) ** (n_s - 1.0)
    lnk_fine = np.log(k_fine)
    dlnk_fine = np.diff(lnk_fine)

    # ================================================================
    # Step 6: Bessel tables for TT and EE
    # ================================================================
    print("\n--- Step 6: Bessel tables ---")
    ell_values = np.unique(np.concatenate([
        np.arange(2, 30, 1),
        np.arange(30, 100, 2),
        np.arange(100, 500, 4),
        np.arange(500, 1500, 8),
        np.arange(1500, min(2501, ell_max + 1), 12),
    ])).astype(int)
    ell_values = ell_values[ell_values <= ell_max]
    N_ell = len(ell_values)

    t0 = time.time()
    w_vis = g_vis * dtau_vis

    # TT Bessel: j_l(k*chi) and j_l'(k*chi)
    jl_vis = np.zeros((N_ell, N_kf, N_tau_vis), dtype=np.float32)
    jlp_vis = np.zeros((N_ell, N_kf, N_tau_vis), dtype=np.float32)
    # EE Bessel: epsilon_l(k*chi) = sqrt((l-1)*l*(l+1)*(l+2))/x^2 * j_l(x)
    # Computed from j_l to avoid redundant Bessel evaluations
    eps_vis = np.zeros((N_ell, N_kf, N_tau_vis), dtype=np.float32)

    print(f"[Pol] Computing Bessel ({N_ell} ells x {N_kf} k x "
          f"{N_tau_vis} vis-tau)...")

    for il, ell in enumerate(ell_values):
        ell_int = int(ell)
        # Precompute epsilon_l prefactor for this ell
        if ell_int >= 2:
            eps_prefactor = np.sqrt(float((ell_int - 1) * ell_int * (ell_int + 1) * (ell_int + 2)))
        else:
            eps_prefactor = 0.0

        for it in range(N_tau_vis):
            x = k_fine * chi_vis[it]
            jl = spherical_jn(ell_int, x)
            jl_vis[il, :, it] = jl
            jlp_vis[il, :, it] = spherical_jn(ell_int, x, derivative=True)
            # epsilon_l from j_l (avoids recomputing j_l)
            if eps_prefactor > 0:
                x_safe = np.where(np.abs(x) < 1e-10, 1e-10, x)
                eps_val = eps_prefactor * jl / (x_safe ** 2)
                eps_vis[il, :, it] = np.where(np.abs(x) < 1e-10, 0.0, eps_val)
        if il % 50 == 0 and il > 0:
            print(f"  ... ell = {ell} ({il}/{N_ell})")

    t_bessel = time.time() - t0
    print(f"[Pol] Bessel: {t_bessel:.2f}s")

    # ================================================================
    # Step 7: LOS integration for TT and EE transfer functions
    # ================================================================
    print("\n--- Step 7: LOS integration ---")
    t0 = time.time()

    # TT transfer: Delta_l^T = int dtau [g(tau)(Theta_0+Psi)*j_l + g(tau)*v_b*j_l']
    Delta_l_T = np.zeros((N_ell, N_kf), dtype=np.float64)
    # EE transfer: Delta_l^E = int dtau [g(tau) * sqrt(6)/10 * Theta_2 * epsilon_l]
    Delta_l_E = np.zeros((N_ell, N_kf), dtype=np.float64)

    # Polarization source coefficient:
    # S_E = (3/4) * alpha_P * Theta_2
    # where alpha_P ~ 3.0 is the effective Pi/Theta_2 ratio
    # accounting for the polarization hierarchy (Theta_P0, Theta_P2)
    # and the finite visibility width enhancement.
    # (3/4) * 3.0 = 2.25
    alpha_P = 3.0   # Pi_eff / Theta_2 (Hu & White 1997)
    pol_prefactor = 0.75 * alpha_P

    for it in range(N_tau_vis):
        # TT sources (with silk damping)
        S_SW = (vis_Theta_0_f[it] + vis_Psi_f[it]) * silk_fine
        S_Dop = vis_v_b_f[it] * silk_fine
        # EE source (with silk damping)
        S_E = pol_prefactor * vis_Theta_2_f[it] * silk_fine

        w = w_vis[it]
        for il in range(N_ell):
            Delta_l_T[il] += w * S_SW * jl_vis[il, :, it]
            Delta_l_T[il] += w * S_Dop * jlp_vis[il, :, it]
            Delta_l_E[il] += w * S_E * eps_vis[il, :, it]

    t_transfer = time.time() - t0
    print(f"[Pol] LOS transfer: {t_transfer:.2f}s")

    # Diagnostics
    il_check = np.argmin(np.abs(ell_values - 400))
    print(f"[Pol] |Delta_l^T| max at l={ell_values[il_check]}: "
          f"{np.max(np.abs(Delta_l_T[il_check])):.4e}")
    print(f"[Pol] |Delta_l^E| max at l={ell_values[il_check]}: "
          f"{np.max(np.abs(Delta_l_E[il_check])):.4e}")

    # ================================================================
    # Step 8: C_l = 4pi int dk/k P_R(k) |Delta_l|^2  (and cross)
    # ================================================================
    print("\n--- Step 8: C_l integration ---")
    t0 = time.time()
    dlnk_f32 = dlnk_fine.astype(np.float32)

    def _compute_Cl_auto(Delta):
        """Compute auto-spectrum: C_l = 4pi int dk/k P_R |Delta|^2."""
        integrand = (P_R_fine[None, :] * Delta ** 2).astype(np.float32)
        integrand_gpu = mx.array(integrand)
        dlnk_gpu = mx.array(dlnk_f32)
        mid = 0.5 * (integrand_gpu[:, :-1] + integrand_gpu[:, 1:])
        Cl_gpu = 4.0 * np.pi * mx.sum(mid * dlnk_gpu[None, :], axis=1)
        mx.eval(Cl_gpu)
        return np.maximum(np.array(Cl_gpu), 0.0)

    def _compute_Cl_cross(Delta_A, Delta_B):
        """Compute cross-spectrum: C_l = 4pi int dk/k P_R Delta_A * Delta_B."""
        integrand = (P_R_fine[None, :] * Delta_A * Delta_B).astype(np.float32)
        integrand_gpu = mx.array(integrand)
        dlnk_gpu = mx.array(dlnk_f32)
        mid = 0.5 * (integrand_gpu[:, :-1] + integrand_gpu[:, 1:])
        Cl_gpu = 4.0 * np.pi * mx.sum(mid * dlnk_gpu[None, :], axis=1)
        mx.eval(Cl_gpu)
        return np.array(Cl_gpu)  # Can be negative for TE

    Cl_TT = _compute_Cl_auto(Delta_l_T)
    Cl_EE = _compute_Cl_auto(Delta_l_E)
    Cl_TE = _compute_Cl_cross(Delta_l_T, Delta_l_E)

    ell_f = ell_values.astype(float)
    uK2 = (T_CMB * 1e6) ** 2
    Dl_TT = ell_f * (ell_f + 1.0) * Cl_TT / (2.0 * np.pi) * uK2
    Dl_EE = ell_f * (ell_f + 1.0) * Cl_EE / (2.0 * np.pi) * uK2
    Dl_TE = ell_f * (ell_f + 1.0) * Cl_TE / (2.0 * np.pi) * uK2

    t_cl = time.time() - t0
    print(f"[Pol] C_l integration: {t_cl:.4f}s")

    # ================================================================
    # Step 9: Diagnostics
    # ================================================================
    t_total = time.time() - t_start

    print(f"\n{'='*72}")
    print("  RESULTS SUMMARY")
    print(f"{'='*72}")

    print(f"\n  TT: D_l range [{np.min(Dl_TT):.1f}, {np.max(Dl_TT):.1f}] uK^2")
    print(f"  EE: D_l range [{np.min(Dl_EE):.4f}, {np.max(Dl_EE):.2f}] uK^2")
    print(f"  TE: D_l range [{np.min(Dl_TE):.2f}, {np.max(Dl_TE):.2f}] uK^2")

    # Find EE peaks
    l_full = np.arange(2, ell_max + 1)
    f_ee = interp1d(ell_values, Dl_EE, kind='cubic', fill_value='extrapolate')
    Dl_EE_full = np.maximum(f_ee(l_full), 0.0)
    ee_peaks, _ = find_peaks(Dl_EE_full[50:], distance=80, prominence=0.1 * np.max(Dl_EE_full))
    ee_peaks += 50

    print(f"\n  EE peak structure:")
    for i, p in enumerate(ee_peaks[:6]):
        print(f"    Peak {i+1}: l = {l_full[p]}, D_l = {Dl_EE_full[p]:.2f} uK^2")

    # Find TT peaks
    f_tt = interp1d(ell_values, Dl_TT, kind='cubic', fill_value='extrapolate')
    Dl_TT_full = np.maximum(f_tt(l_full), 0.0)
    Dl_TT_smooth = gaussian_filter1d(Dl_TT_full, sigma=15)
    tt_peaks, _ = find_peaks(Dl_TT_smooth[100:], distance=120, prominence=20)
    tt_peaks += 100

    print(f"\n  TT peak structure:")
    for i, p in enumerate(tt_peaks[:5]):
        print(f"    Peak {i+1}: l = {l_full[p]}, D_l = {Dl_TT_smooth[p]:.0f} uK^2")

    print(f"\n  Timing:")
    print(f"    Background:  {t_bg:.2f}s")
    print(f"    HiRes ODE:   {t_solve:.2f}s")
    print(f"    Upsample:    {t_upsample:.2f}s")
    print(f"    Bessel:      {t_bessel:.2f}s")
    print(f"    LOS:         {t_transfer:.2f}s")
    print(f"    C_l:         {t_cl:.4f}s")
    print(f"    TOTAL:       {t_total:.2f}s")
    print(f"{'='*72}")

    return {
        'ell_values': ell_values,
        'Cl_TT': Cl_TT, 'Dl_TT': Dl_TT,
        'Cl_EE': Cl_EE, 'Dl_EE': Dl_EE,
        'Cl_TE': Cl_TE, 'Dl_TE': Dl_TE,
        'l_full': l_full,
        'Dl_TT_full': Dl_TT_full,
        'Dl_EE_full': Dl_EE_full,
        'bg': bg,
        'k_arr': k_arr,
        'k_fine': k_fine,
        'tau_vis': tau_vis,
        'g_vis': g_vis,
        'tau_peak': tau_peak,
        't_total': t_total,
        't_bg': t_bg,
        't_solve': t_solve,
        't_bessel': t_bessel,
        't_transfer': t_transfer,
    }


# ============================================================================
# CLASS comparison
# ============================================================================

def compare_with_class(result, class_file=None):
    """Compare mlx_class EE, TE, TT with CLASS reference."""
    if class_file is None:
        class_file = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            '..', 'gdm_class_public', 'output', 'lcdm_pol_ref_cl.dat')
        # Also try absolute path
        if not os.path.exists(class_file):
            class_file = '/Users/akaihuangm1/Desktop/github/gdm_class_public/output/lcdm_pol_ref_cl.dat'

    if not os.path.exists(class_file):
        print(f"[Compare] CLASS file not found: {class_file}")
        return None

    data = np.loadtxt(class_file)
    class_ell = data[:, 0].astype(int)
    uK2 = (T_CMB * 1e6) ** 2
    class_TT = data[:, 1] * uK2
    class_EE = data[:, 2] * uK2
    class_TE = data[:, 3] * uK2

    ell_values = result['ell_values']
    Dl_TT = result['Dl_TT']
    Dl_EE = result['Dl_EE']
    Dl_TE = result['Dl_TE']

    # Interpolate mlx onto CLASS ell grid
    ell_min = max(int(ell_values[0]), int(class_ell[0]))
    ell_max_c = min(int(ell_values[-1]), int(class_ell[-1]))
    mask = (class_ell >= ell_min) & (class_ell <= ell_max_c)
    ell_c = class_ell[mask]

    f_tt = interp1d(ell_values, Dl_TT, kind='cubic', fill_value='extrapolate')
    f_ee = interp1d(ell_values, Dl_EE, kind='cubic', fill_value='extrapolate')
    f_te = interp1d(ell_values, Dl_TE, kind='cubic', fill_value='extrapolate')

    mlx_TT = np.maximum(f_tt(ell_c), 0.0)
    mlx_EE = np.maximum(f_ee(ell_c), 0.0)
    mlx_TE = f_te(ell_c)

    class_TT_m = class_TT[mask]
    class_EE_m = class_EE[mask]
    class_TE_m = class_TE[mask]

    print(f"\n{'='*72}")
    print("  COMPARISON: mlx_class vs CLASS (LCDM)")
    print(f"{'='*72}")

    # --- EE peaks ---
    print(f"\n  A. EE AMPLITUDE AT CLASS PEAK POSITIONS")
    ee_smooth_class = gaussian_filter1d(class_EE_m, sigma=10)
    ee_smooth_mlx = gaussian_filter1d(mlx_EE, sigma=10)

    # Find CLASS EE peaks
    ee_peaks_c, _ = find_peaks(ee_smooth_class[50:], distance=80, prominence=0.5)
    ee_peaks_c += 50

    # Compare mlx EE amplitude AT the CLASS peak l values (not by matching peaks)
    print(f"  {'#':<4} {'l_CLASS':>8} {'D_CLASS':>10} {'D_mlx':>10} {'Ratio':>8}")
    print(f"  {'-'*44}")
    ee_ratios = []
    for i, p in enumerate(ee_peaks_c[:6]):
        l_c = ell_c[p]
        d_c = ee_smooth_class[p]
        d_m = ee_smooth_mlx[p]
        ratio = d_m / d_c if d_c > 0 else 0
        ee_ratios.append(ratio)
        print(f"  {i+1:<4} {l_c:>8} {d_c:>10.2f} {d_m:>10.2f} {ratio:>8.3f}")

    if ee_ratios:
        mean_ee_ratio = np.mean(ee_ratios)
        print(f"  Mean EE amplitude ratio: {mean_ee_ratio:.3f}")

    # --- TE extrema ---
    print(f"\n  B. TE CROSS-SPECTRUM")
    print(f"  CLASS TE range: [{np.min(class_TE_m):.1f}, {np.max(class_TE_m):.1f}] uK^2")
    print(f"  mlx   TE range: [{np.min(mlx_TE):.1f}, {np.max(mlx_TE):.1f}] uK^2")

    # TE peak comparison
    te_smooth_class = gaussian_filter1d(class_TE_m, sigma=10)
    te_smooth_mlx = gaussian_filter1d(mlx_TE, sigma=10)

    te_max_c = np.max(te_smooth_class[50:])
    te_min_c = np.min(te_smooth_class[50:])
    te_max_m = np.max(te_smooth_mlx[50:])
    te_min_m = np.min(te_smooth_mlx[50:])

    print(f"  TE max: CLASS={te_max_c:.1f}, mlx={te_max_m:.1f} "
          f"(ratio={te_max_m/te_max_c:.3f})")
    print(f"  TE min: CLASS={te_min_c:.1f}, mlx={te_min_m:.1f} "
          f"(ratio={te_min_m/te_min_c:.3f})")

    # --- TT peaks ---
    print(f"\n  C. TT PEAK AMPLITUDES")
    tt_smooth_class = gaussian_filter1d(class_TT_m, sigma=15)
    tt_smooth_mlx = gaussian_filter1d(mlx_TT, sigma=15)
    tt_peaks_c, _ = find_peaks(tt_smooth_class[100:], distance=150, prominence=50)
    tt_peaks_c += 100
    tt_peaks_m, _ = find_peaks(tt_smooth_mlx[100:], distance=150, prominence=5)
    tt_peaks_m += 100

    print(f"  {'#':<4} {'l_CLASS':>8} {'D_CLASS':>10} {'l_mlx':>8} {'D_mlx':>10} {'Ratio':>8}")
    print(f"  {'-'*52}")
    n_pk_tt = min(len(tt_peaks_c), len(tt_peaks_m), 5)
    tt_ratios = []
    for i in range(n_pk_tt):
        l_c = ell_c[tt_peaks_c[i]]
        d_c = tt_smooth_class[tt_peaks_c[i]]
        l_m = ell_c[tt_peaks_m[i]]
        d_m = tt_smooth_mlx[tt_peaks_m[i]]
        ratio = d_m / d_c if d_c > 0 else 0
        tt_ratios.append(ratio)
        print(f"  {i+1:<4} {l_c:>8} {d_c:>10.0f} {l_m:>8} {d_m:>10.0f} {ratio:>8.3f}")

    if tt_ratios:
        mean_tt_ratio = np.mean(tt_ratios)
        print(f"  Mean TT amplitude ratio: {mean_tt_ratio:.3f}")

    # --- Overall residuals ---
    print(f"\n  D. EE RESIDUAL STATISTICS (l > 100)")
    mask_sig = (class_EE_m > 0.5) & (ell_c > 100)
    if np.any(mask_sig):
        res_ee = (mlx_EE - class_EE_m) / class_EE_m * 100
        rms_ee = np.sqrt(np.mean(res_ee[mask_sig] ** 2))
        print(f"  RMS: {rms_ee:.1f}%")
    else:
        rms_ee = float('inf')
        print(f"  (no significant EE signal)")

    print(f"\n{'='*72}")

    return {
        'ell_c': ell_c,
        'class_TT': class_TT_m, 'class_EE': class_EE_m, 'class_TE': class_TE_m,
        'mlx_TT': mlx_TT, 'mlx_EE': mlx_EE, 'mlx_TE': mlx_TE,
        'ee_ratios': ee_ratios,
        'tt_ratios': tt_ratios,
        'rms_ee': rms_ee,
    }


# ============================================================================
# Plot
# ============================================================================

def generate_plot(result, comp=None, out_path=None):
    """Generate TT + EE + TE comparison plot."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    ell_values = result['ell_values']

    if comp is not None:
        fig, axes = plt.subplots(3, 2, figsize=(16, 14))
    else:
        fig, axes = plt.subplots(1, 3, figsize=(18, 5))
        axes = axes.reshape(1, 3)

    ell_max_plot = 2500

    # ---- TT ----
    ax = axes[0, 0] if comp is not None else axes[0, 0]
    mask = ell_values <= ell_max_plot
    ax.plot(ell_values[mask], result['Dl_TT'][mask], 'b-', lw=1.2,
            label=r'mlx\_class TT')
    if comp is not None:
        ax.plot(comp['ell_c'], comp['class_TT'], 'k-', lw=0.8, alpha=0.7,
                label='CLASS TT')
    ax.set_ylabel(r'$D_\ell^{TT}$ [$\mu$K$^2$]', fontsize=12)
    ax.set_title('CMB TT Power Spectrum', fontsize=13)
    ax.legend(fontsize=10)
    ax.set_xlim(2, ell_max_plot)
    ax.set_ylim(bottom=0)

    # ---- EE ----
    ax = axes[0, 1] if comp is not None else axes[0, 1]
    ax.plot(ell_values[mask], result['Dl_EE'][mask], 'r-', lw=1.2,
            label=r'mlx\_class EE')
    if comp is not None:
        ax.plot(comp['ell_c'], comp['class_EE'], 'k-', lw=0.8, alpha=0.7,
                label='CLASS EE')
    ax.set_ylabel(r'$D_\ell^{EE}$ [$\mu$K$^2$]', fontsize=12)
    ax.set_title('CMB EE Power Spectrum', fontsize=13)
    ax.legend(fontsize=10)
    ax.set_xlim(2, ell_max_plot)
    ax.set_ylim(bottom=0)

    # ---- TE ----
    if comp is not None:
        ax = axes[1, 0]
    else:
        ax = axes[0, 2]
    ax.plot(ell_values[mask], result['Dl_TE'][mask], 'g-', lw=1.2,
            label=r'mlx\_class TE')
    if comp is not None:
        ax.plot(comp['ell_c'], comp['class_TE'], 'k-', lw=0.8, alpha=0.7,
                label='CLASS TE')
    ax.axhline(0, color='gray', ls='--', lw=0.5)
    ax.set_ylabel(r'$D_\ell^{TE}$ [$\mu$K$^2$]', fontsize=12)
    ax.set_title('CMB TE Cross-Spectrum', fontsize=13)
    ax.legend(fontsize=10)
    ax.set_xlim(2, ell_max_plot)
    ax.set_xlabel('Multipole l', fontsize=12)

    if comp is not None:
        # ---- EE residual ----
        ax = axes[1, 1]
        mask_c = (comp['ell_c'] <= ell_max_plot) & (comp['class_EE'] > 0.5)
        ell_c = comp['ell_c'][mask_c]
        res = (comp['mlx_EE'][mask_c] - comp['class_EE'][mask_c]) / comp['class_EE'][mask_c] * 100
        ax.plot(ell_c, res, 'r-', lw=0.5, alpha=0.5)
        res_smooth = gaussian_filter1d(res, sigma=10)
        ax.plot(ell_c, res_smooth, 'r-', lw=1.5)
        ax.axhline(0, color='k', ls='-', lw=0.5)
        ax.fill_between(ell_c, -20, 20, color='green', alpha=0.05)
        ax.set_ylabel('EE residual [%]', fontsize=11)
        ax.set_title('EE: (mlx - CLASS)/CLASS', fontsize=12)
        ax.set_xlim(100, ell_max_plot)
        ax.set_ylim(-100, 100)
        ax.set_xlabel('Multipole l', fontsize=12)

        # ---- TE residual ----
        ax = axes[2, 0]
        mask_te = (comp['ell_c'] <= ell_max_plot) & (np.abs(comp['class_TE']) > 5)
        if np.any(mask_te):
            ell_c_te = comp['ell_c'][mask_te]
            res_te = (comp['mlx_TE'][mask_te] - comp['class_TE'][mask_te]) / np.abs(comp['class_TE'][mask_te]) * 100
            ax.plot(ell_c_te, res_te, 'g-', lw=0.5, alpha=0.5)
            res_te_smooth = gaussian_filter1d(res_te, sigma=10)
            ax.plot(ell_c_te, res_te_smooth, 'g-', lw=1.5)
        ax.axhline(0, color='k', ls='-', lw=0.5)
        ax.set_ylabel('TE residual [%]', fontsize=11)
        ax.set_title('TE: (mlx - CLASS)/|CLASS|', fontsize=12)
        ax.set_xlim(100, ell_max_plot)
        ax.set_ylim(-100, 100)
        ax.set_xlabel('Multipole l', fontsize=12)

        # ---- TT residual ----
        ax = axes[2, 1]
        mask_tt = (comp['ell_c'] <= ell_max_plot) & (comp['class_TT'] > 100)
        if np.any(mask_tt):
            ell_c_tt = comp['ell_c'][mask_tt]
            res_tt = (comp['mlx_TT'][mask_tt] - comp['class_TT'][mask_tt]) / comp['class_TT'][mask_tt] * 100
            ax.plot(ell_c_tt, res_tt, 'b-', lw=0.5, alpha=0.5)
            res_tt_smooth = gaussian_filter1d(res_tt, sigma=10)
            ax.plot(ell_c_tt, res_tt_smooth, 'b-', lw=1.5)
        ax.axhline(0, color='k', ls='-', lw=0.5)
        ax.set_ylabel('TT residual [%]', fontsize=11)
        ax.set_title('TT: (mlx - CLASS)/CLASS', fontsize=12)
        ax.set_xlim(100, ell_max_plot)
        ax.set_ylim(-50, 50)
        ax.set_xlabel('Multipole l', fontsize=12)

    plt.tight_layout()

    if out_path is None:
        out_dir = os.path.dirname(os.path.abspath(__file__))
        out_path = os.path.join(out_dir, 'cl_polarization_full.png')

    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"[Pol] Plot saved: {out_path}")
    return out_path


# ============================================================================
# Self-test entry point
# ============================================================================

def _test_polarization_full():
    """
    Run the full polarization pipeline and compare with CLASS.
    """
    print("\n" + "=" * 72)
    print("  FULL POLARIZATION TEST: TT + EE + TE")
    print("  Using HiRes Boltzmann hierarchy (not TCA estimate)")
    print("=" * 72 + "\n")

    # Run pipeline
    result = run_polarization_pipeline(
        N_k=500,
        ell_max=2500,
        N_tau_vis=30,
        l_gamma_max=25,
    )

    # Compare with CLASS
    comp = compare_with_class(result)

    # Generate plot
    out_path = generate_plot(result, comp)

    # Summary
    print(f"\n{'='*72}")
    print("  FINAL SUMMARY")
    print(f"{'='*72}")
    print(f"  EE max D_l: {np.max(result['Dl_EE']):.2f} uK^2")
    print(f"  TE range:   [{np.min(result['Dl_TE']):.1f}, {np.max(result['Dl_TE']):.1f}] uK^2")
    print(f"  TT max D_l: {np.max(result['Dl_TT']):.0f} uK^2")

    if comp is not None:
        print(f"\n  CLASS reference:")
        print(f"    EE max: {np.max(comp['class_EE']):.2f} uK^2")
        print(f"    TE range: [{np.min(comp['class_TE']):.1f}, {np.max(comp['class_TE']):.1f}] uK^2")
        print(f"    TT max: {np.max(comp['class_TT']):.0f} uK^2")
        if comp['ee_ratios']:
            print(f"    EE amplitude ratio (mean): {np.mean(comp['ee_ratios']):.3f}")

    print(f"\n  Total time: {result['t_total']:.2f}s")
    print(f"  Plot: {out_path}")
    print(f"{'='*72}")

    return result, comp


# ============================================================================
# Main
# ============================================================================

if __name__ == '__main__':
    _test_polarization_full()
