"""
spectra_los.py -- Proper line-of-sight C_l integration.

PHYSICS:
  Delta_l(k) = int_0^{tau_0} dtau S(k,tau) j_l(k * chi)

  S(k,tau) = g(tau) * [Theta_0 + Psi](k,tau) * j_l(k*chi)            [SW]
           + g(tau) * v_b(k,tau) * j_l'(k*chi)                        [Doppler]
           + e^{-kappa} * [Phi' + Psi'](k,tau) * j_l(k*chi)          [ISW]

  where chi = tau_0 - tau, g(tau) = |kappa_dot| * exp(-kappa).

TWO SOURCE MODES:

  (A) Analytic: Hu-Sugiyama gravitational driving phase shift.
      S_SW(k)  = [A(k)*cos(k*r_s_eff + phi(k)) + delta_DC(k)] * silk
      S_Dop(k) = [-3*c_s/(1+R)*A(k)*sin(k*r_s_eff + phi(k))] * silk
      phi(k) = phi_inf * (1 - exp(-alpha*k/k_eq)), fast-saturating.
      r_s_eff = r_s - 5.0 Mpc (visibility width correction).
      Source FROZEN at tau_*, projected with visibility-weighted Bessel.

  (B) ODE: Run TCA Boltzmann solver THROUGH recombination.
      ODE snapshots at each visibility quadrature point give actual
      Theta_0(k,tau_i), Theta_1(k,tau_i), Phi(k,tau_i), Phi'(k,tau_i).
      No WKB extrapolation. No phase shift calibration.

Both modes include:
  - Visibility function g(tau) integration (~25 points near tau_rec)
  - Early ISW: Phi' at tau_eq from ODE or analytic (l ~ 100-300)
  - Late ISW: analytic Phi' from dark energy decay (l < 20)

PERFORMANCE:
  Bessel dominates: N_ell * N_k * N_tau evaluations of j_l.
  Analytic mode: ~4s. ODE mode: ~15s (includes integration through tau_rec).

Author: Sheng-Kai Huang, 2026
"""
import numpy as np
import mlx.core as mx
import time
from scipy.special import spherical_jn
from scipy.interpolate import interp1d
from scipy.signal import find_peaks

from .background import (Background, A_s, n_s, k_pivot, T_CMB, k_eq,
                          Omega_r as _OMEGA_R, Omega_b as _OMEGA_B,
                          Omega_c as _OMEGA_C, H0_Mpc as _H0_MPC,
                          a_eq as _A_EQ, Omega_m as _OMEGA_M,
                          Omega_L as _OMEGA_L)
from .perturbations import (deriv_tca, rk4_step, adiabatic_ic_tca,
                             TCA_PHI, TCA_DELTA_B,
                             TCA_DELTA_C, TCA_V_C, TCA_THETA_0,
                             TCA_THETA_1, TCA_N_VAR,
                             _poisson_phi_dot)


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
# Analytic source with gravitational driving phase shift
# ============================================================================

def _analytic_source(k_arr, bg):
    """
    Compute the analytic source function at recombination.

    Includes the Hu-Sugiyama gravitational driving phase shift from
    potential decay during radiation domination.

    Physics of the phase shift model:
      - During radiation domination, the gravitational potential Phi decays
        inside the horizon. This drives the photon-baryon oscillation,
        adding a phase shift phi_grav(k) to the acoustic oscillation.
      - The asymptotic phase shift is phi_inf ~ 0.85, calibrated to the
        Planck 2018 peak positions. The standard Hu & Sugiyama (1996)
        value is 0.227*pi ~ 0.71; the larger value here compensates for
        the Bessel projection in the C_l integral.
      - The transition from phi=0 (matter era modes) to phi_inf (radiation
        era modes) uses a fast-saturating exponential:
          phi(k) = phi_inf * (1 - exp(-alpha * k/k_eq))
        By k/k_eq ~ 2, the phase is >99% of phi_inf. This matches the
        physics: modes entering the horizon during radiation domination
        experience the full potential decay.
      - The effective sound horizon r_s_eff = r_s - 5.0 Mpc accounts for
        the finite duration of recombination: the visibility function
        has width sigma ~ 10 Mpc, and modes effectively freeze out at
        slightly different times.

    Calibration (Planck 2018 peaks):
      Peak 1: l=220 (err < 1%), Peak 3: l=819 (err ~ 1%),
      Peak 4: l=1127 (err ~ 1%), Peak 5: l=1424 (err ~ 1.5%).
      Peak 2: l=503 (err ~ 6%) -- this is a known limitation of the
      analytic approximation; the correct l=537.5 requires full Boltzmann
      evolution with baryon drag, photon diffusion, and anisotropic stress.

    Parameters
    ----------
    k_arr : ndarray (N_k,)
    bg : Background

    Returns
    -------
    S_SW : ndarray (N_k,) -- Sachs-Wolfe source (Theta_0 + Phi)
    S_Dop : ndarray (N_k,) -- Doppler source (v_b)
    """
    R_rec = bg.R_rec
    k_D = bg.k_D

    x_eq = k_arr / k_eq

    # Baryon damping: oscillation amplitude reduced by (1+R)^{-1/4}
    baryon_damp = (1.0 + R_rec) ** (-0.25)

    # Sound speed at recombination
    cs_rec = 1.0 / np.sqrt(3.0 * (1.0 + R_rec))

    # Driving amplitude from radiation-era potential decay
    # D(k) = 1 for k << k_eq (matter domination: Phi constant, no driving)
    # D(k) -> 5 for k >> k_eq (radiation: Phi decays, drives oscillation ~5x)
    D_drive = 1.0 + 4.0 * x_eq ** 2 / (1.0 + x_eq ** 2)

    # Oscillation amplitude
    A_osc = (1.0 / 3.0) * D_drive * baryon_damp

    # ----------------------------------------------------------------
    # Gravitational driving phase shift (Hu & Sugiyama 1996)
    # ----------------------------------------------------------------
    # phi -> 0 for k << k_eq (no driving in matter era)
    # phi -> phi_inf for k >> k_eq (radiation-era driving)
    #
    # Fast-saturating exponential transition: by x_eq ~ 2 the phase is
    # already >99% of phi_inf. This is physically correct because modes
    # entering the horizon during radiation domination (x_eq > 1)
    # experience the full potential decay.
    #
    # phi_inf = 0.85 calibrated to Planck peak positions.
    # The standard Hu-Sugiyama value is 0.227*pi ~ 0.71, but the
    # Bessel projection in the C_l integral and early ISW contribution
    # shift the effective peak positions, requiring a larger input value.
    phi_inf = 0.85
    alpha_phi = 3.0
    phi_grav = phi_inf * (1.0 - np.exp(-alpha_phi * x_eq))

    # ----------------------------------------------------------------
    # Effective sound horizon
    # ----------------------------------------------------------------
    # The finite duration of recombination means different k-modes
    # effectively "freeze" at slightly different times, giving an
    # effective r_s that differs from the instantaneous value.
    # r_s_eff ~ r_s - delta_r_s, where delta_r_s ~ cs * sigma_vis
    # sigma_vis ~ 10 Mpc, cs ~ 0.57, so delta_r_s ~ 5 Mpc
    # Calibrated: r_s_eff = bg.r_s - 5.0 Mpc
    r_s_eff = bg.r_s - 5.0

    # ----------------------------------------------------------------
    # Baryon zero-shift from baryon loading
    # ----------------------------------------------------------------
    # This shifts the zero-point of oscillation, creating odd-even asymmetry
    # zero_shift = R/(3*(1+R)) * Phi(k)/Phi_0 at recombination
    # For k << k_eq: Phi ~ 0.9*Phi_0, zero_shift ~ 0.128 * 0.9 = 0.115
    # For k >> k_eq: Phi ~ 0, zero_shift ~ 0
    # Use Eisenstein-Hu transfer function for Phi(k)
    _T_k = _eisenstein_hu_transfer(k_arr)
    zero_shift = R_rec / (3.0 * (1.0 + R_rec)) * 0.9 * _T_k

    # Silk damping
    silk = np.exp(-(k_arr / k_D) ** 2)

    # Sachs-Wolfe source: Theta_0 + Phi
    S_SW = (A_osc * np.cos(k_arr * r_s_eff + phi_grav) + zero_shift) * silk

    # Doppler source: v_b = 3*Theta_1
    # In TCA: Theta_1 = -(c_s/(1+R)) * A * sin(kr_s + phi)
    S_Dop = -3.0 * cs_rec / (1.0 + R_rec) * A_osc * np.sin(
        k_arr * r_s_eff + phi_grav
    ) * silk

    return S_SW, S_Dop


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


# ============================================================================
# Extended ODE solver: integrates THROUGH recombination with snapshots
# ============================================================================

def _build_extended_tau_grid(bg, k_max, tau_end, N_early=300, N_late=800):
    """Build integration grid extending to tau_end (can be > tau_rec)."""
    tau_init = bg.tau_grid[1]
    tau_early_end = min(10.0, 0.3 * tau_end)

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


def _solve_with_snapshots_extended(bg, k_arr_np, tau_snapshots):
    """
    Run the TCA ODE solver THROUGH recombination and record state at
    specified tau snapshots. Returns actual Theta_0, Theta_1, Phi, Phi_dot.
    """
    k_arr = mx.array(k_arr_np.astype(np.float32))
    N_k = len(k_arr_np)

    tau_end = float(tau_snapshots[-1]) + 1.0
    tau_end = min(tau_end, bg.tau_0 * 0.99)
    tau_grid = _build_extended_tau_grid(bg, float(k_arr_np[-1]), tau_end)

    y = adiabatic_ic_tca(k_arr_np, bg)

    calH_all = mx.array(bg.calH_at_tau(tau_grid).astype(np.float32))
    R_all = mx.array(bg.R_at_tau(tau_grid).astype(np.float32))
    a_all = mx.array(bg.a_at_tau(tau_grid).astype(np.float32))
    _Or, _Ob, _Oc, _H0 = _OMEGA_R, _OMEGA_B, float(bg.Omega_cdm), _H0_MPC

    N_snap = len(tau_snapshots)
    snap_Theta_0 = np.zeros((N_snap, N_k), dtype=np.float32)
    snap_Theta_1 = np.zeros((N_snap, N_k), dtype=np.float32)
    snap_Phi = np.zeros((N_snap, N_k), dtype=np.float32)
    snap_Phi_dot = np.zeros((N_snap, N_k), dtype=np.float32)

    snap_idx = 0
    for i in range(len(tau_grid) - 1):
        while snap_idx < N_snap and tau_grid[i] >= tau_snapshots[snap_idx]:
            mx.eval(y)
            y_np = np.array(y)
            snap_Theta_0[snap_idx] = y_np[:, TCA_THETA_0]
            snap_Theta_1[snap_idx] = y_np[:, TCA_THETA_1]
            snap_Phi[snap_idx] = y_np[:, TCA_PHI]
            ci, ri, ai = calH_all[i], R_all[i], a_all[i]
            rhs = deriv_tca(y, k_arr, ci, ri, _Or, _Ob, _Oc, ai, _H0)
            mx.eval(rhs)
            snap_Phi_dot[snap_idx] = np.array(rhs)[:, TCA_PHI]
            snap_idx += 1

        dt = float(tau_grid[i + 1] - tau_grid[i])
        ci, ri, ai = calH_all[i], R_all[i], a_all[i]
        def rhs_fn(y_in, _c=ci, _R=ri, _a=ai):
            return deriv_tca(y_in, k_arr, _c, _R, _Or, _Ob, _Oc, _a, _H0)
        y = rk4_step(y, rhs_fn, dt)
        if i % 200 == 0:
            mx.eval(y)

    mx.eval(y)
    while snap_idx < N_snap:
        y_np = np.array(y)
        snap_Theta_0[snap_idx] = y_np[:, TCA_THETA_0]
        snap_Theta_1[snap_idx] = y_np[:, TCA_THETA_1]
        snap_Phi[snap_idx] = y_np[:, TCA_PHI]
        rhs = deriv_tca(y, k_arr, calH_all[-1], R_all[-1],
                        _Or, _Ob, _Oc, a_all[-1], _H0)
        mx.eval(rhs)
        snap_Phi_dot[snap_idx] = np.array(rhs)[:, TCA_PHI]
        snap_idx += 1

    return {
        'Theta_0': snap_Theta_0, 'Theta_1': snap_Theta_1,
        'Phi': snap_Phi, 'Phi_dot': snap_Phi_dot,
    }


# ============================================================================
# Early ISW: analytic Phi_dot during radiation-matter transition
# ============================================================================

def _early_isw_phi_dot_analytic(k_arr, tau_arr, bg):
    """
    Compute Phi_dot(k,tau) during the radiation-matter transition.
    Sub-horizon Phi ~ T(x) = 3(sin x - x cos x)/x^3, x = k*tau/sqrt(3).
    """
    N_tau = len(tau_arr)
    a_arr = bg.a_at_tau(tau_arr)
    Phi_dot = np.zeros((N_tau, len(k_arr)), dtype=np.float64)

    for it in range(N_tau):
        tau = tau_arr[it]
        a = a_arr[it]
        f_rad = _OMEGA_R / (a * (_OMEGA_R / a + _OMEGA_M))

        x = k_arr * tau / np.sqrt(3.0)
        x_safe = np.where(x > 1e-6, x, 1.0)
        sin_x = np.sin(x_safe)
        cos_x = np.cos(x_safe)

        dT_dx = np.where(
            x > 1e-6,
            3.0 * ((x_safe**2 - 3.0) * sin_x + 3.0 * x_safe * cos_x) / x_safe**4,
            0.0
        )
        Phi_dot[it] = f_rad * dT_dx * k_arr / np.sqrt(3.0)

    return Phi_dot


# ============================================================================
# Late ISW: Phi_dot from dark energy
# ============================================================================

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
# Main LOS integration
# ============================================================================

def compute_cl_los(bg, k_arr_Mpc, ell_values, use_ode=False,
                   N_tau_vis=25, N_tau_early_isw=12, N_tau_late_isw=10):
    """
    Compute C_l using proper line-of-sight integration.

    Parameters
    ----------
    bg : Background
        Solved background cosmology.
    k_arr_Mpc : ndarray (N_k,)
        Wavenumber array in Mpc^-1.
    ell_values : ndarray (N_ell,)
        Multipole values.
    use_ode : bool
        If True, run ODE solver through recombination for source functions.
        If False (default), use analytic source with Hu-Sugiyama phase shift.
    N_tau_vis : int
        Number of quadrature points for visibility integral.
    N_tau_early_isw : int
        Number of tau points for early ISW integral.
    N_tau_late_isw : int
        Number of tau points for late ISW integral.

    Returns
    -------
    ell_values, Cl, Dl, info_dict
    """
    t_total = time.time()

    k_arr = np.asarray(k_arr_Mpc, dtype=np.float64)
    N_k = len(k_arr)
    N_ell = len(ell_values)

    # Primordial power spectrum
    P_R = A_s * (k_arr / k_pivot) ** (n_s - 1.0)
    lnk = np.log(k_arr)
    dlnk = np.diff(lnk).astype(np.float32)

    # Silk damping
    k_D = bg.k_D
    silk = np.exp(-(k_arr / k_D) ** 2)

    # ================================================================
    # Build tau quadrature grids
    # ================================================================

    # 1. Visibility region
    tau_vis, g_vis, tau_peak, sigma_vis = _build_visibility_quadrature(
        bg, N_tau_vis
    )
    chi_vis = bg.tau_0 - tau_vis
    dtau_vis = tau_vis[1] - tau_vis[0]

    # 2. Early ISW: from tau_eq/3 to start of visibility region
    tau_eq = float(bg._tau_of_a(np.log(_A_EQ)))
    tau_early_lo = max(tau_eq * 0.3, bg.tau_grid[10])
    tau_early_hi = tau_vis[0]
    if tau_early_hi > tau_early_lo * 1.5:
        tau_early_isw = np.geomspace(tau_early_lo, tau_early_hi, N_tau_early_isw)
    else:
        tau_early_isw = np.array([tau_early_lo])
    chi_early_isw = bg.tau_0 - tau_early_isw
    kappa_early = bg.kappa_at_tau(tau_early_isw)
    exp_neg_kappa_early = np.exp(-kappa_early)
    N_early_pts = len(tau_early_isw)

    # 3. Late ISW: from end of visibility region to near tau_0
    tau_late_lo = tau_vis[-1]
    tau_late_hi = bg.tau_0 * 0.98
    tau_late = np.linspace(tau_late_lo, tau_late_hi, N_tau_late_isw)
    chi_late = bg.tau_0 - tau_late
    kappa_late = bg.kappa_at_tau(tau_late)
    exp_neg_kappa_late = np.exp(-kappa_late)
    N_late_pts = len(tau_late)

    N_total_tau = N_tau_vis + N_early_pts + N_late_pts
    print(f"[LOS] Visibility: {N_tau_vis} pts around tau={tau_peak:.1f} Mpc "
          f"(sigma={sigma_vis:.1f} Mpc)")
    print(f"[LOS] Early ISW: {N_early_pts} pts "
          f"[{tau_early_isw[0]:.1f}, {tau_early_isw[-1]:.1f}] Mpc (tau_eq={tau_eq:.1f})")
    print(f"[LOS] Late ISW: {N_late_pts} pts "
          f"[{tau_late[0]:.0f}, {tau_late[-1]:.0f}] Mpc")
    print(f"[LOS] Total quadrature: {N_total_tau} tau points")

    # ================================================================
    # Step 1: Source functions
    # ================================================================
    t0 = time.time()

    if use_ode:
        # Run ODE through recombination: snapshots at early ISW + visibility
        all_ode_tau = np.sort(np.unique(np.concatenate([
            tau_early_isw, tau_vis
        ])))
        print(f"[LOS] Running ODE with {len(all_ode_tau)} snapshots "
              f"through tau={all_ode_tau[-1]:.1f} Mpc...")
        snapshots = _solve_with_snapshots_extended(
            bg, k_arr.astype(np.float32), all_ode_tau
        )
        t_ode = time.time() - t0
        print(f"[LOS] ODE done: {t_ode:.2f}s")

        def _find_snap_idx(tau_val):
            return int(np.argmin(np.abs(all_ode_tau - tau_val)))

        # Extract visibility region sources
        vis_Theta_0 = np.zeros((N_tau_vis, N_k), dtype=np.float64)
        vis_Theta_1 = np.zeros((N_tau_vis, N_k), dtype=np.float64)
        vis_Phi = np.zeros((N_tau_vis, N_k), dtype=np.float64)
        vis_Phi_dot = np.zeros((N_tau_vis, N_k), dtype=np.float64)
        for it in range(N_tau_vis):
            idx = _find_snap_idx(tau_vis[it])
            vis_Theta_0[it] = snapshots['Theta_0'][idx]
            vis_Theta_1[it] = snapshots['Theta_1'][idx]
            vis_Phi[it] = snapshots['Phi'][idx]
            vis_Phi_dot[it] = snapshots['Phi_dot'][idx]

        # Early ISW from ODE
        early_Phi_dot = np.zeros((N_early_pts, N_k), dtype=np.float64)
        for it in range(N_early_pts):
            idx = _find_snap_idx(tau_early_isw[it])
            early_Phi_dot[it] = snapshots['Phi_dot'][idx]

        source_mode = 'ode'
        print(f"[LOS] Source: ODE through recombination")
    else:
        # Analytic source with Hu-Sugiyama phase shift
        S_SW, S_Dop = _analytic_source(k_arr, bg)
        # Early ISW from analytic model
        early_Phi_dot = _early_isw_phi_dot_analytic(k_arr, tau_early_isw, bg)
        source_mode = 'analytic'
        print(f"[LOS] Source: analytic with Hu-Sugiyama driving phase")
        t_ode = time.time() - t0

    # Late ISW (always analytic)
    late_Phi_dot = _late_isw_phi_dot(k_arr, tau_late, bg)

    # ================================================================
    # Step 2: Bessel tables
    # ================================================================
    t0 = time.time()
    print(f"[LOS] Computing Bessel ({N_ell} ells x {N_k} k x {N_total_tau} tau)...")

    # Visibility Bessels: accumulate <g*j_l> or store full 3D table for ODE
    w_all = g_vis * dtau_vis

    if source_mode == 'analytic':
        # For analytic: accumulate visibility-weighted Bessel directly
        gjl = np.zeros((N_ell, N_k), dtype=np.float64)
        gjlp = np.zeros((N_ell, N_k), dtype=np.float64)
        for il, ell in enumerate(ell_values):
            ell_int = int(ell)
            for it in range(N_tau_vis):
                x = k_arr * chi_vis[it]
                jl_val = spherical_jn(ell_int, x)
                jlp_val = spherical_jn(ell_int, x, derivative=True)
                gjl[il] += w_all[it] * jl_val
                gjlp[il] += w_all[it] * jlp_val
    else:
        # For ODE: need per-tau Bessel values (source varies with tau)
        jl_vis = np.zeros((N_ell, N_k, N_tau_vis), dtype=np.float32)
        jlp_vis = np.zeros((N_ell, N_k, N_tau_vis), dtype=np.float32)
        for il, ell in enumerate(ell_values):
            ell_int = int(ell)
            for it in range(N_tau_vis):
                x = k_arr * chi_vis[it]
                jl_vis[il, :, it] = spherical_jn(ell_int, x)
                jlp_vis[il, :, it] = spherical_jn(ell_int, x, derivative=True)

    # Early ISW Bessels
    jl_early = np.zeros((N_ell, N_k, N_early_pts), dtype=np.float32)
    for il, ell in enumerate(ell_values):
        ell_int = int(ell)
        for it in range(N_early_pts):
            jl_early[il, :, it] = spherical_jn(ell_int, k_arr * chi_early_isw[it])

    # Late ISW Bessels
    jl_late = np.zeros((N_ell, N_k, N_late_pts), dtype=np.float32)
    for il, ell in enumerate(ell_values):
        ell_int = int(ell)
        for it in range(N_late_pts):
            jl_late[il, :, it] = spherical_jn(ell_int, k_arr * chi_late[it])

    t_bessel = time.time() - t0
    print(f"[LOS] Bessel: {t_bessel:.2f}s")

    # ================================================================
    # Step 3: Transfer function Delta_l(k)
    # ================================================================
    t0 = time.time()

    Delta_l_sw = np.zeros((N_ell, N_k), dtype=np.float64)
    Delta_l_dop = np.zeros((N_ell, N_k), dtype=np.float64)
    Delta_l_eisw = np.zeros((N_ell, N_k), dtype=np.float64)
    Delta_l_lisw = np.zeros((N_ell, N_k), dtype=np.float64)

    # --- SW + Doppler from visibility integral ---
    if source_mode == 'analytic':
        for il in range(N_ell):
            Delta_l_sw[il] = S_SW * gjl[il]
            Delta_l_dop[il] = S_Dop * gjlp[il]
    else:
        # ODE: source varies with tau
        for it in range(N_tau_vis):
            S0 = (vis_Theta_0[it] + vis_Phi[it]) * silk
            S1 = 3.0 * vis_Theta_1[it] * silk
            w = w_all[it]
            for il in range(N_ell):
                Delta_l_sw[il] += w * S0 * jl_vis[il, :, it]
                Delta_l_dop[il] += w * S1 * jlp_vis[il, :, it]

        # Also add ISW from within visibility region (Phi_dot nonzero near tau_rec)
        kappa_vis = bg.kappa_at_tau(tau_vis)
        exp_neg_kappa_vis = np.exp(-kappa_vis)
        for it in range(N_tau_vis):
            S_ISW_vis = exp_neg_kappa_vis[it] * 2.0 * vis_Phi_dot[it]
            for il in range(N_ell):
                Delta_l_eisw[il] += dtau_vis * S_ISW_vis * jl_vis[il, :, it]

    Delta_l = Delta_l_sw + Delta_l_dop

    # --- Early ISW integral ---
    if N_early_pts > 1:
        dtau_early = np.diff(tau_early_isw)
        for it in range(N_early_pts - 1):
            integrand_lo = exp_neg_kappa_early[it] * 2.0 * early_Phi_dot[it]
            integrand_hi = exp_neg_kappa_early[it + 1] * 2.0 * early_Phi_dot[it + 1]
            integrand_avg = 0.5 * (integrand_lo + integrand_hi)
            dt = dtau_early[it]
            for il in range(N_ell):
                jl_avg = 0.5 * (jl_early[il, :, it] + jl_early[il, :, it + 1])
                Delta_l_eisw[il] += integrand_avg * jl_avg * dt

    Delta_l += Delta_l_eisw

    # --- Late ISW integral ---
    if N_late_pts > 1:
        dtau_late = np.diff(tau_late)
        for it in range(N_late_pts - 1):
            integrand_lo = exp_neg_kappa_late[it] * 2.0 * late_Phi_dot[it]
            integrand_hi = exp_neg_kappa_late[it + 1] * 2.0 * late_Phi_dot[it + 1]
            integrand_avg = 0.5 * (integrand_lo + integrand_hi)
            dt = dtau_late[it]
            for il in range(N_ell):
                jl_avg = 0.5 * (jl_late[il, :, it] + jl_late[il, :, it + 1])
                Delta_l_lisw[il] += integrand_avg * jl_avg * dt

    Delta_l += Delta_l_lisw

    t_transfer = time.time() - t0
    print(f"[LOS] Transfer: {t_transfer:.2f}s")

    # ================================================================
    # Step 4: C_l = 4pi int dk/k P_R |Delta_l|^2
    # ================================================================
    t0 = time.time()
    integrand = (P_R[None, :] * Delta_l ** 2).astype(np.float32)

    integrand_gpu = mx.array(integrand)
    dlnk_gpu = mx.array(dlnk)
    mid = 0.5 * (integrand_gpu[:, :-1] + integrand_gpu[:, 1:])
    Cl_gpu = 4.0 * np.pi * mx.sum(mid * dlnk_gpu[None, :], axis=1)
    mx.eval(Cl_gpu)

    Cl = np.maximum(np.array(Cl_gpu), 0.0)
    ell_f = ell_values.astype(float)
    Dl = ell_f * (ell_f + 1.0) * Cl / (2.0 * np.pi) * (T_CMB * 1e6) ** 2

    # Diagnostics
    def _cl_from_delta(Delta):
        integ = (P_R[None, :] * Delta ** 2).astype(np.float32)
        integ_gpu = mx.array(integ)
        m = 0.5 * (integ_gpu[:, :-1] + integ_gpu[:, 1:])
        cl = 4.0 * np.pi * mx.sum(m * dlnk_gpu[None, :], axis=1)
        mx.eval(cl)
        cl_np = np.maximum(np.array(cl), 0.0)
        return ell_f * (ell_f + 1.0) * cl_np / (2.0 * np.pi) * (T_CMB * 1e6) ** 2

    Dl_sw_only = _cl_from_delta(Delta_l_sw)
    Dl_sw_dop = _cl_from_delta(Delta_l_sw + Delta_l_dop)

    t_cl = time.time() - t0
    print(f"[LOS] C_l integration: {t_cl:.4f}s")

    t_elapsed = time.time() - t_total
    print(f"[LOS] Total: {t_elapsed:.2f}s")

    info = {
        'Dl_sw_only': Dl_sw_only,
        'Dl_sw_dop': Dl_sw_dop,
        'Delta_l_sw': Delta_l_sw,
        'Delta_l_dop': Delta_l_dop,
        'Delta_l_eisw': Delta_l_eisw,
        'Delta_l_lisw': Delta_l_lisw,
        'tau_vis': tau_vis,
        'g_vis': g_vis,
        'tau_peak': tau_peak,
        'sigma_vis': sigma_vis,
        'source_mode': source_mode,
    }

    return ell_values, Cl, Dl, info


# ============================================================================
# Peak finder
# ============================================================================

def find_peaks_in_spectrum(ell_values, Dl, n_peaks=7, smooth_sigma=15):
    """Find acoustic peak positions in D_l spectrum.

    Uses Gaussian smoothing to average over Bessel oscillations
    (period ~50 in l) and find the envelope acoustic peaks.
    Default smooth_sigma=15 is chosen to smooth Bessel oscillations
    while preserving acoustic peaks (spacing ~300 in l).
    """
    from scipy.ndimage import gaussian_filter1d

    ell_fine = np.arange(int(ell_values[0]), int(ell_values[-1]) + 1)
    f = interp1d(ell_values, Dl, kind='cubic', fill_value='extrapolate')
    Dl_fine = np.maximum(f(ell_fine), 0.0)

    Dl_smooth = gaussian_filter1d(Dl_fine, sigma=smooth_sigma)
    peaks, props = find_peaks(Dl_smooth, distance=80, prominence=50)

    result = []
    for p in peaks[:n_peaks]:
        result.append((int(ell_fine[p]), float(Dl_smooth[p])))

    return result, ell_fine, Dl_fine, Dl_smooth
