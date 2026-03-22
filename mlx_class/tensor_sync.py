"""
tensor_sync.py -- BB power spectrum from gravitational waves (synchronous gauge framework).

Tensor perturbations (gravitational waves) are gauge-invariant, so the physics
is identical to tensor.py. The difference is purely computational: we use
scipy's stiff ODE solvers (Radau/BDF) instead of GPU-batched RK4 + Strang
splitting, matching the solver_sync.py pattern.

PHYSICS:
  Tensor mode h_T(k, tau) satisfies the source-free wave equation:
    h_T'' + 2 calH h_T' + k^2 h_T = 0

  Inflationary ICs: h_T(tau_init) = 1, h_T'(tau_init) = 0.

  Two coupled photon hierarchies (Zaldarriaga & Seljak 1997):
    Temperature: Theta_T_l  (sourced by h_T'/6 at l=0)
    Polarization: Theta_P_l (sourced through Pi = Theta_T_2 + Theta_P_0 + Theta_P_2)
    Thomson collision: kappa_dot terms couple the hierarchies.

  BB source (line-of-sight):
    Delta_l^B(k) = int dtau g(tau) * S_BB(k, tau) * epsilon_l(k*chi)
    where S_BB = alpha_BB * h_T'(k, tau)  and
          epsilon_l(x) = sqrt((l-1)*l*(l+1)*(l+2)) * j_l(x) / x^2

  Power spectrum:
    C_l^BB = 4 pi int dk/k  P_T(k) |Delta_l^B(k)|^2
    P_T(k) = r * A_s * (k/k_pivot)^{n_T}   with n_T = -r/8 (consistency relation)

Usage:
    python -m mlx_class.tensor_sync --r 0.01
    python -m mlx_class.tensor_sync --r 0.01 0.1 --plot

Author: Sheng-Kai Huang, 2026
"""
import os
os.environ['MPLCONFIGDIR'] = '/tmp/mpl_config'

import sys
import time
import numpy as np
from multiprocessing import Pool

from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
from scipy.special import spherical_jn
from scipy.signal import find_peaks
from scipy.ndimage import gaussian_filter1d

from .background import (
    Background, A_s, n_s, k_pivot, T_CMB,
    Omega_r as _OMEGA_R,
    Omega_b as _OMEGA_B,
    Omega_m as _OMEGA_M,
    H0_Mpc as _H0_MPC,
)


# ============================================================================
# Constants for the tensor photon hierarchy
# ============================================================================

L_TENSOR_MAX = 10  # Default truncation for tensor photon hierarchy


def _n_tensor_var(l_max):
    """Total number of variables per k-mode: h, h', 2*(l_max+1) photon moments."""
    # State: [h, h', Theta_T_0..l_max, Theta_P_0..l_max]
    return 2 + 2 * (l_max + 1)


# ============================================================================
# Tensor ODE right-hand side (full Boltzmann, scipy-compatible)
# ============================================================================

def make_tensor_rhs(k, bg, l_max=L_TENSOR_MAX):
    """
    Build the RHS function for the combined tensor wave + photon hierarchy ODE.

    State vector layout:
      y[0] = h_T          (tensor metric perturbation)
      y[1] = h_T'         (time derivative)
      y[2..2+l_max]       = Theta_T_0 .. Theta_T_{l_max}  (temperature)
      y[2+l_max+1..end]   = Theta_P_0 .. Theta_P_{l_max}  (polarization)

    Parameters
    ----------
    k : float
        Wavenumber in Mpc^{-1}.
    bg : Background
        Solved background object.
    l_max : int
        Truncation multipole for photon hierarchy.

    Returns
    -------
    rhs : callable(tau, y) -> dy/dtau
    nvar : int
    """
    # Precompute background interpolators as numpy arrays for speed
    tau_grid = bg.tau_grid.copy()
    calH_grid = bg.calH_grid.copy()
    kappa_dot_grid = bg.kappa_dot_grid.copy()

    calH_interp = interp1d(tau_grid, calH_grid, kind='cubic',
                           fill_value='extrapolate')
    kd_interp = interp1d(tau_grid, kappa_dot_grid, kind='cubic',
                         fill_value='extrapolate')

    k2 = k * k
    nvar = _n_tensor_var(l_max)

    # Index helpers
    I_H = 0
    I_HP = 1
    I_T = lambda l: 2 + l                   # Theta_T_l
    I_P = lambda l: 2 + (l_max + 1) + l     # Theta_P_l

    def rhs(tau, y):
        calH = float(calH_interp(tau))
        kd = float(kd_interp(tau))  # kd < 0

        dy = np.zeros(nvar)

        h = y[I_H]
        hp = y[I_HP]

        # --- Tensor wave equation: h'' + 2 calH h' + k^2 h = 0 ---
        dy[I_H] = hp
        dy[I_HP] = -2.0 * calH * hp - k2 * h

        # --- Temperature hierarchy ---
        T = np.zeros(l_max + 2)  # T[l_max+1] = 0 (truncation)
        for l in range(l_max + 1):
            T[l] = y[I_T(l)]

        P = np.zeros(l_max + 2)  # P[l_max+1] = 0 (truncation)
        for l in range(l_max + 1):
            P[l] = y[I_P(l)]

        # Polarization source: Pi = T_2 + P_0 + P_2
        Pi = T[2] + P[0] + P[2] if l_max >= 2 else 0.0

        # Theta_T_0: -k T_1 - h'/6
        dy[I_T(0)] = -k * T[1] - hp / 6.0

        # Theta_T_1: k/3 (T_0 - 2 T_2) + kd T_1
        if l_max >= 1:
            dy[I_T(1)] = k / 3.0 * (T[0] - 2.0 * T[2]) + kd * T[1]

        # Theta_T_2: k/5 (2 T_1 - 3 T_3) + kd (T_2 - Pi/10)
        if l_max >= 2:
            dy[I_T(2)] = (k / 5.0 * (2.0 * T[1] - 3.0 * T[3])
                          + kd * (T[2] - Pi / 10.0))

        # Theta_T_l, l=3..l_max
        for l in range(3, l_max + 1):
            dy[I_T(l)] = (k / (2.0 * l + 1.0) * (l * T[l-1] - (l + 1.0) * T[l+1])
                          + kd * T[l])

        # --- Polarization hierarchy ---
        # Theta_P_0: -k P_1 + kd (P_0 - Pi/2)
        dy[I_P(0)] = -k * P[1] + kd * (P[0] - Pi / 2.0)

        # Theta_P_1: k/3 (P_0 - 2 P_2) + kd P_1
        if l_max >= 1:
            dy[I_P(1)] = k / 3.0 * (P[0] - 2.0 * P[2]) + kd * P[1]

        # Theta_P_2: k/5 (2 P_1 - 3 P_3) + kd (P_2 - Pi/10)
        if l_max >= 2:
            dy[I_P(2)] = (k / 5.0 * (2.0 * P[1] - 3.0 * P[3])
                          + kd * (P[2] - Pi / 10.0))

        # Theta_P_l, l=3..l_max
        for l in range(3, l_max + 1):
            dy[I_P(l)] = (k / (2.0 * l + 1.0) * (l * P[l-1] - (l + 1.0) * P[l+1])
                          + kd * P[l])

        return dy

    return rhs, nvar


# ============================================================================
# Initial conditions
# ============================================================================

def tensor_ic(l_max=L_TENSOR_MAX):
    """
    Adiabatic tensor initial conditions.

    h_T = 1 (unit amplitude, primordial spectrum provides power)
    h_T' = 0
    All photon moments = 0 (tight coupling drives them from h')
    """
    nvar = _n_tensor_var(l_max)
    y0 = np.zeros(nvar)
    y0[0] = 1.0   # h_T = 1
    y0[1] = 0.0   # h_T' = 0
    return y0


# ============================================================================
# Solve single tensor k-mode
# ============================================================================

def solve_tensor_k(k, bg, tau_end, l_max=L_TENSOR_MAX,
                   method='Radau', rtol=1e-8, atol=1e-10):
    """
    Solve the tensor Boltzmann system for a single k using scipy.

    Returns the OdeSolution object with dense_output=True.
    """
    tau_init = bg.tau_grid[1]
    rhs, nvar = make_tensor_rhs(k, bg, l_max)
    y0 = tensor_ic(l_max)

    sol = solve_ivp(rhs, [tau_init, tau_end], y0,
                    method=method, dense_output=True,
                    rtol=rtol, atol=atol)
    return sol


# ============================================================================
# Visibility quadrature (same as tensor.py)
# ============================================================================

def _build_visibility_quadrature(bg, N_tau=40):
    """Build quadrature points around the visibility peak for tensor LOS."""
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

    # Extend to 5 sigma for thorough coverage
    tau_lo = max(tau_peak - 5.0 * sigma, bg.tau_grid[1])
    tau_hi = min(tau_peak + 5.0 * sigma, bg.tau_grid[-2])
    tau_vis = np.linspace(tau_lo, tau_hi, N_tau)
    g_vis = bg.visibility_at_tau(tau_vis)

    return tau_vis, g_vis, tau_peak, sigma


# ============================================================================
# Spin-2 Bessel projection
# ============================================================================

def _epsilon_l(ell, x_arr):
    """
    Spin-2 Bessel projection factor for B-mode polarization.
    epsilon_l(x) = sqrt((l-1)*l*(l+1)*(l+2)) / x^2 * j_l(x)
    """
    l = int(ell)
    if l < 2:
        return np.zeros_like(x_arr)
    prefactor = np.sqrt(float((l - 1) * l * (l + 1) * (l + 2)))
    jl = spherical_jn(l, x_arr)
    x_safe = np.where(np.abs(x_arr) < 1e-10, 1e-10, x_arr)
    result = prefactor * jl / (x_safe ** 2)
    return np.where(np.abs(x_arr) < 1e-10, 0.0, result)


# ============================================================================
# Multiprocessing worker for tensor k-modes
# ============================================================================

def _tensor_worker(args):
    """
    Worker: solve one tensor k-mode and extract h_T' at snapshot taus.

    Returns numpy array of shape (N_snap, 3): [h, h', Pi] at each tau snapshot.
    Returns None on failure.
    """
    import numpy as np
    from scipy.integrate import solve_ivp
    from scipy.interpolate import interp1d

    (k, bg_arrays, tau_end, l_max, method, rtol, atol, tau_snaps) = args

    # Reconstruct lightweight background
    class _BGLite:
        __slots__ = ('tau_grid', 'a_grid', 'calH_grid', 'kappa_dot_grid',
                     'tau_rec', 'tau_0', 'R_grid')
        def __init__(self, arr):
            self.tau_grid = arr['tau_grid']
            self.a_grid = arr['a_grid']
            self.calH_grid = arr['calH_grid']
            self.kappa_dot_grid = arr['kappa_dot_grid']
            self.tau_rec = arr['tau_rec']
            self.tau_0 = arr['tau_0']
            self.R_grid = arr.get('R_grid', np.zeros_like(self.a_grid))

    bg_lite = _BGLite(bg_arrays)

    # Import and build RHS
    from mlx_class.tensor_sync import make_tensor_rhs, tensor_ic

    rhs_fn, nvar = make_tensor_rhs(k, bg_lite, l_max)
    y0 = tensor_ic(l_max)

    tau_init = bg_lite.tau_grid[1]

    sol = solve_ivp(rhs_fn, [tau_init, tau_end], y0,
                    method=method, dense_output=True,
                    rtol=rtol, atol=atol)

    if not sol.success:
        return None

    N_snap = len(tau_snaps)
    # Columns: h, h', Pi
    results = np.zeros((N_snap, 3), dtype=np.float64)

    I_T = lambda l: 2 + l
    I_P = lambda l: 2 + (l_max + 1) + l

    for it in range(N_snap):
        tau = tau_snaps[it]
        if tau < tau_init:
            # Before ODE start: ICs
            results[it, 0] = 1.0
            results[it, 1] = 0.0
            results[it, 2] = 0.0
            continue

        y = sol.sol(tau)
        results[it, 0] = y[0]   # h_T
        results[it, 1] = y[1]   # h_T'

        if l_max >= 2:
            T2 = y[I_T(2)]
            P0 = y[I_P(0)]
            P2 = y[I_P(2)]
            results[it, 2] = T2 + P0 + P2   # Pi
        else:
            results[it, 2] = 0.0

    return results


# ============================================================================
# Main tensor BB pipeline
# ============================================================================

def run_tensor_sync(r=0.01, n_t=None, N_k=300, ell_max=500,
                    N_tau_vis=40, l_max=L_TENSOR_MAX,
                    method='Radau', rtol=1e-8, atol=1e-10,
                    n_workers=None, verbose=True):
    """
    Full tensor BB pipeline using synchronous-gauge-framework scipy solvers.

    Tensor perturbations are gauge-invariant, so the ODE is the same in all
    gauges. The "sync" label here refers to the solver framework (scipy Radau,
    multiprocessing), not the gauge.

    Parameters
    ----------
    r : float
        Tensor-to-scalar ratio.
    n_t : float or None
        Tensor spectral tilt. None => -r/8 (consistency relation).
    N_k : int
        Number of k-modes.
    ell_max : int
        Maximum multipole for BB spectrum.
    N_tau_vis : int
        Number of visibility quadrature points.
    l_max : int
        Truncation multipole for tensor photon hierarchy.
    method : str
        ODE solver method ('Radau' or 'BDF').
    n_workers : int or None
        Number of parallel workers. None => auto.
    verbose : bool
        Print progress.

    Returns
    -------
    dict with keys: ell, Dl_BB, Dl_TT_tensor, Dl_EE_tensor, r, n_t, bg, ...
    """
    if n_t is None:
        n_t = -r / 8.0

    if n_workers is None:
        n_workers = min(os.cpu_count() or 4, 8)

    t_start = time.time()

    if verbose:
        print("=" * 72)
        print("  TENSOR BB PIPELINE (sync framework, scipy Radau)")
        print(f"  r = {r}, n_t = {n_t:.6f}")
        print(f"  N_k = {N_k}, ell_max = {ell_max}, l_max = {l_max}")
        print(f"  Workers = {n_workers}, method = {method}")
        print("=" * 72)

    # ================================================================
    # Step 1: Background
    # ================================================================
    if verbose:
        print("\n--- Step 1: Background ---")
    t0 = time.time()
    bg = Background(khronon=False, recombination='recfast')
    bg.solve()
    t_bg = time.time() - t0
    if verbose:
        print(f"Background: {t_bg:.2f}s")

    # ================================================================
    # Step 2: Visibility quadrature
    # ================================================================
    if verbose:
        print("\n--- Step 2: Visibility quadrature ---")
    tau_vis, g_vis, tau_peak, sigma_vis = _build_visibility_quadrature(bg, N_tau_vis)
    chi_vis = bg.tau_0 - tau_vis
    dtau_vis = tau_vis[1] - tau_vis[0]

    if verbose:
        print(f"Visibility: {N_tau_vis} pts around tau={tau_peak:.1f} Mpc "
              f"(sigma={sigma_vis:.1f} Mpc)")

    # ================================================================
    # Step 3: k-grid and primordial spectrum
    # ================================================================
    D_A = bg.D_A
    k_min = max(2e-4, 2.0 / D_A)
    k_max = min(0.25, max(0.2, (ell_max + 100.0) / D_A))
    k_arr = np.geomspace(k_min, k_max, N_k)

    P_T = r * A_s * (k_arr / k_pivot) ** n_t

    if verbose:
        print(f"\n--- Step 3: k-grid ---")
        print(f"k range: [{k_arr[0]:.2e}, {k_arr[-1]:.2e}] Mpc^-1, N_k={N_k}")

    # ================================================================
    # Step 4: Solve tensor ODE for all k-modes (parallel)
    # ================================================================
    if verbose:
        print(f"\n--- Step 4: Tensor ODE ({N_k} modes, {n_workers} workers) ---")
        sys.stdout.flush()

    t0 = time.time()

    tau_end = float(tau_vis[-1] + 10.0)
    tau_end = min(tau_end, bg.tau_0 * 0.99)

    # Serialize background
    bg_arrays = {
        'tau_grid': bg.tau_grid.copy(),
        'a_grid': bg.a_grid.copy(),
        'calH_grid': bg.calH_grid.copy(),
        'kappa_dot_grid': bg.kappa_dot_grid.copy(),
        'tau_rec': float(bg.tau_rec),
        'tau_0': float(bg.tau_0),
    }

    work_items = [
        (k, bg_arrays, tau_end, l_max, method, rtol, atol, tau_vis)
        for k in k_arr
    ]

    with Pool(processes=n_workers) as pool:
        worker_results = pool.map(_tensor_worker, work_items)

    t_ode = time.time() - t0

    # Collect results
    snap_h = np.zeros((N_tau_vis, N_k), dtype=np.float64)
    snap_hdot = np.zeros((N_tau_vis, N_k), dtype=np.float64)
    snap_Pi = np.zeros((N_tau_vis, N_k), dtype=np.float64)

    n_failed = 0
    for ik, res in enumerate(worker_results):
        if res is None:
            n_failed += 1
            continue
        snap_h[:, ik] = res[:, 0]
        snap_hdot[:, ik] = res[:, 1]
        snap_Pi[:, ik] = res[:, 2]

    if verbose:
        print(f"Tensor ODE: {t_ode:.2f}s ({n_failed} failed, {N_k - n_failed} OK)")

        # Diagnostics at visibility peak
        peak_snap = np.argmin(np.abs(tau_vis - tau_peak))
        print(f"  h  at vis peak: [{np.min(snap_h[peak_snap]):.4e}, "
              f"{np.max(snap_h[peak_snap]):.4e}]")
        print(f"  h' at vis peak: [{np.min(snap_hdot[peak_snap]):.4e}, "
              f"{np.max(snap_hdot[peak_snap]):.4e}]")
        print(f"  Pi at vis peak: [{np.min(snap_Pi[peak_snap]):.4e}, "
              f"{np.max(snap_Pi[peak_snap]):.4e}]")

    # ================================================================
    # Step 5: Source upsampling to fine k-grid
    # ================================================================
    if verbose:
        print(f"\n--- Step 5: Upsample sources ---")
    t0 = time.time()

    dk_fine = np.pi / (2.0 * D_A) / 1.5
    k_max_fine = min(float(k_arr[-1]), (ell_max + 200.0) / D_A)
    k_fine = np.arange(float(k_arr[0]), k_max_fine, dk_fine).astype(np.float64)
    N_kf = len(k_fine)

    # Silk damping for tensor modes
    silk_fine = np.exp(-(k_fine / bg.k_D) ** 2)

    # Upsample h' onto fine grid
    vis_hdot_f = np.zeros((N_tau_vis, N_kf), dtype=np.float64)
    for it in range(N_tau_vis):
        f = interp1d(k_arr, snap_hdot[it], kind='cubic',
                     fill_value=0.0, bounds_error=False)
        vis_hdot_f[it] = f(k_fine)

    P_T_fine = r * A_s * (k_fine / k_pivot) ** n_t
    lnk_fine = np.log(k_fine)
    dlnk_fine = np.diff(lnk_fine)

    t_upsample = time.time() - t0
    if verbose:
        print(f"Upsampled: {N_k} -> {N_kf} k-points ({t_upsample:.2f}s)")

    # ================================================================
    # Step 6: Bessel tables
    # ================================================================
    if verbose:
        print(f"\n--- Step 6: Bessel tables ---")

    ell_values = np.unique(np.concatenate([
        np.arange(2, 20, 1),
        np.arange(20, 50, 2),
        np.arange(50, 150, 3),
        np.arange(150, min(501, ell_max + 1), 5),
    ])).astype(int)
    ell_values = ell_values[ell_values <= ell_max]
    N_ell = len(ell_values)

    t0 = time.time()
    w_vis = g_vis * dtau_vis  # quadrature weights

    # Precompute Bessel j_l and epsilon_l at all (ell, k, tau_vis) combinations
    jl_vis = np.zeros((N_ell, N_kf, N_tau_vis), dtype=np.float64)
    eps_vis = np.zeros((N_ell, N_kf, N_tau_vis), dtype=np.float64)

    if verbose:
        print(f"Computing Bessel ({N_ell} ells x {N_kf} k x {N_tau_vis} vis-tau)...")
        sys.stdout.flush()

    for il, ell in enumerate(ell_values):
        ell_int = int(ell)
        if ell_int >= 2:
            eps_prefactor = np.sqrt(float((ell_int - 1) * ell_int *
                                          (ell_int + 1) * (ell_int + 2)))
        else:
            eps_prefactor = 0.0

        for it in range(N_tau_vis):
            x = k_fine * chi_vis[it]
            jl = spherical_jn(ell_int, x)
            jl_vis[il, :, it] = jl

            if eps_prefactor > 0:
                x_safe = np.where(np.abs(x) < 1e-10, 1e-10, x)
                eps_val = eps_prefactor * jl / (x_safe ** 2)
                eps_vis[il, :, it] = np.where(np.abs(x) < 1e-10, 0.0, eps_val)

        if il % 30 == 0 and il > 0 and verbose:
            print(f"  ... ell = {ell} ({il}/{N_ell})")

    t_bessel = time.time() - t0
    if verbose:
        print(f"Bessel: {t_bessel:.2f}s")

    # ================================================================
    # Step 7: LOS integration
    # ================================================================
    if verbose:
        print(f"\n--- Step 7: LOS integration ---")
    t0 = time.time()

    # BB source coefficient (calibrated to match CLASS).
    # alpha_BB = sqrt(6) * alpha_P_tensor / 4
    # where alpha_P_tensor is tuned so that for r=0.1, D_l^BB_peak ~ 0.01 uK^2.
    alpha_P_tensor = 1.87
    alpha_BB = np.sqrt(6.0) * alpha_P_tensor / 4.0   # ~ 1.15

    Delta_l_BB = np.zeros((N_ell, N_kf), dtype=np.float64)
    Delta_l_TT_tensor = np.zeros((N_ell, N_kf), dtype=np.float64)
    Delta_l_EE_tensor = np.zeros((N_ell, N_kf), dtype=np.float64)

    for it in range(N_tau_vis):
        S_BB = alpha_BB * vis_hdot_f[it] * silk_fine
        S_TT = vis_hdot_f[it] / 2.0 * silk_fine
        S_EE = S_BB  # tensor EE has same source as BB

        w = w_vis[it]
        for il in range(N_ell):
            Delta_l_BB[il] += w * S_BB * eps_vis[il, :, it]
            Delta_l_EE_tensor[il] += w * S_EE * eps_vis[il, :, it]
            Delta_l_TT_tensor[il] += w * S_TT * jl_vis[il, :, it]

    t_los = time.time() - t0
    if verbose:
        print(f"LOS integration: {t_los:.2f}s")

        il_check = np.argmin(np.abs(ell_values - 80))
        print(f"|Delta_l^BB| max at l={ell_values[il_check]}: "
              f"{np.max(np.abs(Delta_l_BB[il_check])):.4e}")

    # ================================================================
    # Step 8: C_l integration
    # ================================================================
    if verbose:
        print(f"\n--- Step 8: C_l integration ---")
    t0 = time.time()

    def _compute_Cl(Delta):
        integrand = P_T_fine[None, :] * Delta ** 2
        mid = 0.5 * (integrand[:, :-1] + integrand[:, 1:])
        Cl = 4.0 * np.pi * np.sum(mid * dlnk_fine[None, :], axis=1)
        return np.maximum(Cl, 0.0)

    Cl_BB = _compute_Cl(Delta_l_BB)
    Cl_TT_tensor = _compute_Cl(Delta_l_TT_tensor)
    Cl_EE_tensor = _compute_Cl(Delta_l_EE_tensor)

    ell_f = ell_values.astype(float)
    uK2 = (T_CMB * 1e6) ** 2

    Dl_BB = ell_f * (ell_f + 1.0) * Cl_BB / (2.0 * np.pi) * uK2
    Dl_TT_tensor = ell_f * (ell_f + 1.0) * Cl_TT_tensor / (2.0 * np.pi) * uK2
    Dl_EE_tensor = ell_f * (ell_f + 1.0) * Cl_EE_tensor / (2.0 * np.pi) * uK2

    t_cl = time.time() - t0
    if verbose:
        print(f"C_l integration: {t_cl:.4f}s")

    # ================================================================
    # Step 9: Interpolate to full ell range
    # ================================================================
    l_full = np.arange(2, ell_max + 1)

    f_bb = interp1d(ell_values, Dl_BB, kind='cubic', fill_value='extrapolate')
    Dl_BB_full = np.maximum(f_bb(l_full), 0.0)

    f_tt_t = interp1d(ell_values, Dl_TT_tensor, kind='cubic',
                      fill_value='extrapolate')
    Dl_TT_tensor_full = np.maximum(f_tt_t(l_full), 0.0)

    f_ee_t = interp1d(ell_values, Dl_EE_tensor, kind='cubic',
                      fill_value='extrapolate')
    Dl_EE_tensor_full = np.maximum(f_ee_t(l_full), 0.0)

    # ================================================================
    # Step 10: Summary
    # ================================================================
    t_total = time.time() - t_start

    bb_smooth = gaussian_filter1d(Dl_BB_full, sigma=3)
    bb_peaks, _ = find_peaks(bb_smooth, distance=30,
                             prominence=0.001 * np.max(bb_smooth + 1e-30))

    if verbose:
        print(f"\n{'='*72}")
        print("  TENSOR BB RESULTS (sync framework)")
        print(f"{'='*72}")
        print(f"\n  r = {r}, n_t = {n_t:.6f}")
        print(f"\n  BB:  D_l range [{np.min(Dl_BB):.4e}, {np.max(Dl_BB):.4e}] uK^2")
        print(f"  TT (tensor): [{np.min(Dl_TT_tensor):.4e}, "
              f"{np.max(Dl_TT_tensor):.4e}] uK^2")
        print(f"  EE (tensor): [{np.min(Dl_EE_tensor):.4e}, "
              f"{np.max(Dl_EE_tensor):.4e}] uK^2")

        if len(bb_peaks) > 0:
            print(f"\n  BB peak structure:")
            for i, p in enumerate(bb_peaks[:5]):
                print(f"    Peak {i+1}: l = {l_full[p]}, "
                      f"D_l = {Dl_BB_full[p]:.4e} uK^2")
        else:
            peak_l = l_full[np.argmax(Dl_BB_full)]
            peak_val = np.max(Dl_BB_full)
            print(f"\n  BB maximum: l = {peak_l}, D_l = {peak_val:.4e} uK^2")

        print(f"\n  Reference: CLASS r=0.1 -> D_l^BB ~ 0.01 uK^2 at l ~ 80-100")
        print(f"  => for r={r}, expect D_l^BB ~ {r * 0.1:.4e} uK^2")

        print(f"\n  Timing:")
        print(f"    Background:  {t_bg:.2f}s")
        print(f"    Tensor ODE:  {t_ode:.2f}s")
        print(f"    Upsample:    {t_upsample:.2f}s")
        print(f"    Bessel:      {t_bessel:.2f}s")
        print(f"    LOS:         {t_los:.2f}s")
        print(f"    C_l:         {t_cl:.4f}s")
        print(f"    TOTAL:       {t_total:.2f}s")
        print(f"{'='*72}")

    return {
        'r': r,
        'n_t': n_t,
        'ell_values': ell_values,
        'Cl_BB': Cl_BB,
        'Dl_BB': Dl_BB,
        'Cl_TT_tensor': Cl_TT_tensor,
        'Dl_TT_tensor': Dl_TT_tensor,
        'Cl_EE_tensor': Cl_EE_tensor,
        'Dl_EE_tensor': Dl_EE_tensor,
        'l_full': l_full,
        'Dl_BB_full': Dl_BB_full,
        'Dl_TT_tensor_full': Dl_TT_tensor_full,
        'Dl_EE_tensor_full': Dl_EE_tensor_full,
        'bg': bg,
        'k_arr': k_arr,
        'k_fine': k_fine,
        'tau_vis': tau_vis,
        'g_vis': g_vis,
        'tau_peak': tau_peak,
        'snap_h': snap_h,
        'snap_hdot': snap_hdot,
        'snap_Pi': snap_Pi,
        'timing': {
            'background': t_bg,
            'tensor_ode': t_ode,
            'upsample': t_upsample,
            'bessel': t_bessel,
            'los': t_los,
            'cl': t_cl,
            'total': t_total,
        },
    }


# ============================================================================
# Compare multiple r values
# ============================================================================

def compare_r_values(r_values=None, **kwargs):
    """Run tensor pipeline for multiple r values and compare."""
    if r_values is None:
        r_values = [0.01, 0.1]

    results = {}
    for r_val in r_values:
        print(f"\n{'#'*72}")
        print(f"  Running r = {r_val}")
        print(f"{'#'*72}")
        results[r_val] = run_tensor_sync(r=r_val, **kwargs)

    print(f"\n{'='*72}")
    print("  COMPARISON ACROSS r VALUES")
    print(f"{'='*72}")
    print(f"\n  {'r':>8} {'BB max (uK^2)':>15} {'BB peak l':>10} "
          f"{'TT_tensor max':>15} {'EE_tensor max':>15}")
    print(f"  {'-'*68}")

    for r_val in r_values:
        res = results[r_val]
        bb_max = np.max(res['Dl_BB_full'])
        bb_peak_l = res['l_full'][np.argmax(res['Dl_BB_full'])]
        tt_max = np.max(res['Dl_TT_tensor_full'])
        ee_max = np.max(res['Dl_EE_tensor_full'])
        print(f"  {r_val:>8.3f} {bb_max:>15.4e} {bb_peak_l:>10} "
              f"{tt_max:>15.4e} {ee_max:>15.4e}")

    if len(r_values) >= 2:
        r1, r2 = r_values[0], r_values[1]
        bb1 = np.max(results[r1]['Dl_BB_full'])
        bb2 = np.max(results[r2]['Dl_BB_full'])
        if bb1 > 0:
            actual_ratio = bb2 / bb1
            expected_ratio = r2 / r1
            print(f"\n  BB scaling check:")
            print(f"    BB(r={r2})/BB(r={r1}) = {actual_ratio:.3f}")
            print(f"    Expected (r2/r1):       {expected_ratio:.3f}")
            print(f"    Agreement:              {actual_ratio/expected_ratio:.3f}")

    print(f"{'='*72}")
    return results


# ============================================================================
# Plot
# ============================================================================

def generate_plot(results, out_path=None):
    """Generate tensor BB (+ TT, EE) comparison plot."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    if out_path is None:
        out_path = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            'cl_tensor_bb_sync.png')

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    colors = {0.001: 'green', 0.01: 'blue', 0.03: 'cyan',
              0.05: 'orange', 0.1: 'red', 0.2: 'darkred'}

    # BB
    ax = axes[0]
    for r_val, res in sorted(results.items()):
        c = colors.get(r_val, 'purple')
        ax.plot(res['l_full'], res['Dl_BB_full'], color=c, lw=1.5,
                label=f'r = {r_val}')
    ax.set_xlabel(r'$\ell$', fontsize=12)
    ax.set_ylabel(r'$D_\ell^{BB}$ [$\mu$K$^2$]', fontsize=12)
    ax.set_title('CMB BB Power Spectrum (Tensor, sync)', fontsize=13)
    ax.legend(fontsize=10)
    ax.set_xlim(2, 500)
    ax.set_ylim(bottom=0)

    # TT tensor
    ax = axes[1]
    for r_val, res in sorted(results.items()):
        c = colors.get(r_val, 'purple')
        ax.plot(res['l_full'], res['Dl_TT_tensor_full'], color=c, lw=1.5,
                label=f'r = {r_val}')
    ax.set_xlabel(r'$\ell$', fontsize=12)
    ax.set_ylabel(r'$D_\ell^{TT,\mathrm{tensor}}$ [$\mu$K$^2$]', fontsize=12)
    ax.set_title('CMB TT (Tensor contribution)', fontsize=13)
    ax.legend(fontsize=10)
    ax.set_xlim(2, 500)

    # EE tensor
    ax = axes[2]
    for r_val, res in sorted(results.items()):
        c = colors.get(r_val, 'purple')
        ax.plot(res['l_full'], res['Dl_EE_tensor_full'], color=c, lw=1.5,
                label=f'r = {r_val}')
    ax.set_xlabel(r'$\ell$', fontsize=12)
    ax.set_ylabel(r'$D_\ell^{EE,\mathrm{tensor}}$ [$\mu$K$^2$]', fontsize=12)
    ax.set_title('CMB EE (Tensor contribution)', fontsize=13)
    ax.legend(fontsize=10)
    ax.set_xlim(2, 500)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\n[Tensor sync] Plot saved: {out_path}")


# ============================================================================
# CLI
# ============================================================================

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        description='Tensor BB power spectrum (sync framework)')
    parser.add_argument('--r', type=float, nargs='+', default=[0.01],
                        help='Tensor-to-scalar ratio(s). Default: 0.01')
    parser.add_argument('--ell-max', type=int, default=500,
                        help='Maximum multipole')
    parser.add_argument('--N-k', type=int, default=300,
                        help='Number of k-modes')
    parser.add_argument('--l-max', type=int, default=L_TENSOR_MAX,
                        help='Photon hierarchy truncation')
    parser.add_argument('--method', default='Radau',
                        choices=['Radau', 'BDF', 'RK45'],
                        help='ODE solver method')
    parser.add_argument('--workers', type=int, default=None,
                        help='Number of parallel workers')
    parser.add_argument('--plot', action='store_true', default=True,
                        help='Generate plot (default: True)')
    parser.add_argument('--no-plot', dest='plot', action='store_false')
    args = parser.parse_args()

    common_kwargs = dict(
        N_k=args.N_k, ell_max=args.ell_max, l_max=args.l_max,
        method=args.method, n_workers=args.workers,
    )

    if len(args.r) == 1:
        result = run_tensor_sync(r=args.r[0], **common_kwargs)
        results = {args.r[0]: result}

        # Print key D_l^BB values
        ell_ref = [10, 50, 80, 100, 150, 200, 300, 500]
        print(f"\nD_l^BB at key multipoles (r = {args.r[0]}):")
        print(f"  {'ell':>6}  {'D_l^BB (uK^2)':>15}")
        print(f"  {'-'*24}")
        for ell in ell_ref:
            if ell <= args.ell_max:
                idx = ell - 2  # l_full starts at l=2
                if idx < len(result['Dl_BB_full']):
                    print(f"  {ell:>6}  {result['Dl_BB_full'][idx]:>15.4e}")
    else:
        results = compare_r_values(r_values=args.r, **common_kwargs)

    if args.plot:
        generate_plot(results)
