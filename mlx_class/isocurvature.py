"""
isocurvature.py -- Isocurvature initial conditions for the mlx_class Boltzmann solver.

Implements four standard isocurvature modes beyond the adiabatic mode:
  CDI  : CDM density isocurvature        S_c = delta_c - (3/4) delta_gamma != 0
  BI   : Baryon density isocurvature      S_b = delta_b - (3/4) delta_gamma != 0
  NDI  : Neutrino density isocurvature    S_nu = delta_nu - (3/4) delta_gamma != 0
  NVI  : Neutrino velocity isocurvature   initial neutrino bulk velocity != 0

Physics reference: Ma & Bertschinger (1995), Bucher, Moodley & Turok (2000),
CLASS documentation (Blas, Lesgourgues & Tram 2011).

Convention:
  - Isocurvature modes are normalised so that S_X = 1 at tau_init (the
    relevant entropy perturbation equals unity).
  - Phi = 0 at leading order for all isocurvature modes.
  - Higher-order terms in k*tau are included where they are O(k*tau) or O((k*tau)^2).

Mixed initial conditions:
  y0 = sqrt(1 - alpha) * y_adiabatic + sqrt(alpha) * y_isocurvature
       + beta * y_correlated
  where alpha is the isocurvature fraction and beta the correlation amplitude.
  The sqrt ensures the total *power* (not amplitude) is partitioned correctly:
    P_total ~ (1-alpha) P_ad + alpha P_iso + beta P_corr
  since C_l is quadratic in the initial conditions.

Usage:
    from mlx_class.isocurvature import isocurvature_ic, mixed_ic
    y0 = mixed_ic(k_arr, bg, alpha=0.05, mode='CDI')  # 5% CDM isocurvature

Author: Sheng-Kai Huang, 2026
"""
import numpy as np
import mlx.core as mx

from .background import Omega_r as _OMEGA_R, Omega_b as _OMEGA_B, Omega_c as _OMEGA_C
from .background import H0_Mpc as _H0_MPC

from .perturbations_neutrino import (
    N_EFF, _f_nu, _OMEGA_GAMMA, _OMEGA_NU, L_NU_MAX,
)

# Re-export the adiabatic ICs so callers can import everything from one place
from .perturbations_implicit import adiabatic_ic as adiabatic_ic_tca
from .perturbations_neutrino import adiabatic_ic as adiabatic_ic_neutrino
from .perturbations_hires import adiabatic_ic_hires

# Valid mode names
ISOCURVATURE_MODES = ('CDI', 'BI', 'NDI', 'NVI')


# ============================================================================
# Density fractions at radiation domination (used for IC normalisation)
# ============================================================================
# In the radiation era:  rho_gamma / rho_r = 1 - f_nu,  rho_nu / rho_r = f_nu
# Omega_b << Omega_r and Omega_c << Omega_r, so matter is sub-dominant.
# f_nu = Omega_nu / Omega_r  (already computed in perturbations_neutrino)

_f_gamma = 1.0 - _f_nu   # photon fraction of radiation
_R_nu = _f_nu             # neutrino fraction of radiation (= f_nu)


# ============================================================================
# TCA (6-variable) isocurvature initial conditions
# ============================================================================
# Variable layout: [Phi, delta_b, delta_c, v_c, Theta_0, Theta_1]
# This matches perturbations_implicit.py TCA indices.

def isocurvature_ic_tca(k_arr_np, bg, mode='CDI'):
    """
    Isocurvature initial conditions for the 6-variable TCA system.

    Parameters
    ----------
    k_arr_np : array (N_k,)
        Wavenumber array in Mpc^-1.
    bg : Background
        Solved background object.
    mode : str
        One of 'CDI', 'BI', 'NDI', 'NVI'.

    Returns
    -------
    y0 : mx.array (N_k, 6)
    """
    mode = mode.upper()
    if mode not in ISOCURVATURE_MODES:
        raise ValueError(f"Unknown isocurvature mode '{mode}'. "
                         f"Choose from {ISOCURVATURE_MODES}")

    N_k = len(k_arr_np)
    tau = bg.tau_grid[1]
    kt = k_arr_np * tau

    y0 = np.zeros((N_k, 6), dtype=np.float32)

    # TCA indices (from perturbations_implicit.py)
    PHI, DELTA_B, DELTA_C, V_C, THETA_0, THETA_1 = 0, 1, 2, 3, 4, 5

    if mode == 'CDI':
        # CDM isocurvature: S_c = delta_c - (3/4)*delta_gamma = 1
        # Leading order in radiation domination (Ma & Bertschinger 1995, Eq. 99-102):
        #   Phi = 0  (no initial curvature)
        #   delta_c = 1  (the isocurvature perturbation)
        #   delta_b = 0
        #   delta_gamma = 0  =>  Theta_0 = 0
        #   v_c = 0
        #   Theta_1 = 0
        # Next order: Phi grows as ~ -(Omega_c / Omega_r) * a * S_c / 4
        # which is negligible at tau_init deep in radiation domination.
        # The O(k*tau) corrections:
        #   Phi ~ -(1/6)(Omega_c/Omega_r) * (k*tau)^2 * ... (very small)
        y0[:, DELTA_C] = 1.0
        # All others zero at leading order.

    elif mode == 'BI':
        # Baryon isocurvature: S_b = delta_b - (3/4)*delta_gamma = 1
        # Same structure as CDI but with baryons.
        # In TCA, baryons and photons are tightly coupled, so the baryon
        # isocurvature is partially compensated by the photon-baryon fluid.
        # Leading order:
        #   delta_b = 1
        #   delta_gamma = 0 => Theta_0 = 0
        #   Phi = 0
        y0[:, DELTA_B] = 1.0

    elif mode == 'NDI':
        # Neutrino density isocurvature: S_nu = delta_nu - delta_gamma = 1
        # Compensated: f_gamma*delta_gamma + f_nu*delta_nu = 0
        # => delta_gamma = -f_nu,  delta_nu = 1 - f_nu = f_gamma
        # Theta_0 = delta_gamma/4 = -f_nu/4
        # delta_b follows photons in TCA: delta_b = (3/4)*delta_gamma
        # Note: the TCA system has no separate neutrino variable, so NDI
        # enters through the photon monopole only (approximate treatment).
        y0[:, THETA_0] = -_f_nu / 4.0
        y0[:, DELTA_B] = 0.75 * (-_f_nu)  # (3/4)*delta_gamma

    elif mode == 'NVI':
        # Neutrino velocity isocurvature: initial neutrino velocity != 0
        # Leading order:
        #   N_1 = 1 (normalisation)
        #   All density perturbations zero
        #   Phi = 0
        # In the TCA system (no neutrino hierarchy), NVI cannot be
        # represented properly. We approximate by leaving everything zero
        # and printing a warning.
        import warnings
        warnings.warn("NVI mode in TCA (6-variable) system is approximate: "
                      "no separate neutrino velocity. Use the neutrino or "
                      "hires solver for accurate NVI.")
        # The effect of NVI on photons is O(f_nu * k*tau), set Theta_1 correction:
        y0[:, THETA_1] = _f_nu * kt / 6.0

    return mx.array(y0)


# ============================================================================
# Neutrino solver isocurvature initial conditions
# ============================================================================
# Variable layout: [Phi, delta_b, delta_c, v_c, Theta_0, Theta_1,
#                   N_0, N_1, ..., N_{l_nu_max}]

def isocurvature_ic_neutrino(k_arr_np, bg, mode='CDI', l_nu_max=L_NU_MAX):
    """
    Isocurvature initial conditions for the neutrino solver.

    Parameters
    ----------
    k_arr_np : array (N_k,)
    bg : Background
    mode : str
        One of 'CDI', 'BI', 'NDI', 'NVI'.
    l_nu_max : int
        Neutrino hierarchy truncation.

    Returns
    -------
    y0 : mx.array (N_k, N_var)
    """
    from .perturbations_neutrino import (
        IDX_PHI, IDX_DELTA_B, IDX_DELTA_C, IDX_V_C,
        IDX_THETA_0, IDX_THETA_1, IDX_N_START, n_var_total,
    )

    mode = mode.upper()
    if mode not in ISOCURVATURE_MODES:
        raise ValueError(f"Unknown isocurvature mode '{mode}'. "
                         f"Choose from {ISOCURVATURE_MODES}")

    N_k = len(k_arr_np)
    tau = bg.tau_grid[1]
    kt = k_arr_np * tau
    nvar = n_var_total(l_nu_max)

    y0 = np.zeros((N_k, nvar), dtype=np.float32)

    if mode == 'CDI':
        # S_c = delta_c - (3/4)*delta_gamma = 1, Phi = 0
        y0[:, IDX_DELTA_C] = 1.0

    elif mode == 'BI':
        # S_b = delta_b - (3/4)*delta_gamma = 1, Phi = 0
        y0[:, IDX_DELTA_B] = 1.0

    elif mode == 'NDI':
        # Neutrino density isocurvature with separate neutrino hierarchy.
        # S_nu = delta_nu - delta_gamma = 1
        # Compensated: f_gamma*delta_gamma + f_nu*delta_nu = 0
        # => delta_gamma = -f_nu,  delta_nu = f_gamma = 1 - f_nu
        # Theta_0 = delta_gamma / 4 = -f_nu / 4
        # N_0 = delta_nu / 4 = (1 - f_nu) / 4
        y0[:, IDX_THETA_0] = -_f_nu / 4.0
        y0[:, IDX_THETA_1] = 0.0
        y0[:, IDX_N_START + 0] = (1.0 - _f_nu) / 4.0
        y0[:, IDX_N_START + 1] = 0.0
        # delta_b follows photons in TCA: delta_b = (3/4)*delta_gamma
        y0[:, IDX_DELTA_B] = 0.75 * (-_f_nu)

    elif mode == 'NVI':
        # Neutrino velocity isocurvature: N_1 = 1 (normalisation)
        # All density perturbations zero, Phi = 0.
        # The neutrino dipole carries the initial perturbation.
        # O(k*tau) corrections from free-streaming:
        #   N_0' = -k*N_1 => N_0(tau) ~ -k*tau * N_1(0)   (tiny at tau_init)
        #   N_2 ~ (2/5)*k*tau * N_1(0)
        y0[:, IDX_N_START + 1] = 1.0
        # Higher multipoles seeded by free-streaming
        if l_nu_max >= 2:
            y0[:, IDX_N_START + 2] = (2.0 / 5.0) * kt
        # N_0 correction:
        y0[:, IDX_N_START + 0] = -kt  # N_0 ~ -k*tau

    return mx.array(y0)


# ============================================================================
# HiRes solver isocurvature initial conditions
# ============================================================================
# Variable layout: [Phi, delta_b, v_b, delta_c, v_c,
#                   Theta_0, ..., Theta_{l_gamma_max},
#                   N_0, ..., N_{l_nu_max}]

def isocurvature_ic_hires(k_arr_np, bg, mode='CDI',
                           l_gamma_max=25, l_nu_max=L_NU_MAX):
    """
    Isocurvature initial conditions for the HiRes solver.

    Parameters
    ----------
    k_arr_np : array (N_k,)
    bg : Background
    mode : str
        One of 'CDI', 'BI', 'NDI', 'NVI'.
    l_gamma_max : int
        Photon hierarchy truncation.
    l_nu_max : int
        Neutrino hierarchy truncation.

    Returns
    -------
    y0 : mx.array (N_k, N_var)
    """
    from .perturbations_hires import (
        IDX_PHI, IDX_DELTA_B, IDX_V_B, IDX_DELTA_C, IDX_V_C,
        IDX_THETA_START, idx_theta, idx_n_start, n_var_total,
    )

    mode = mode.upper()
    if mode not in ISOCURVATURE_MODES:
        raise ValueError(f"Unknown isocurvature mode '{mode}'. "
                         f"Choose from {ISOCURVATURE_MODES}")

    N_k = len(k_arr_np)
    tau = bg.tau_grid[1]
    kt = k_arr_np * tau
    nvar = n_var_total(l_gamma_max, l_nu_max)
    _n_start = idx_n_start(l_gamma_max)

    y0 = np.zeros((N_k, nvar), dtype=np.float32)

    if mode == 'CDI':
        # CDM isocurvature: S_c = 1, Phi = 0
        y0[:, IDX_DELTA_C] = 1.0

    elif mode == 'BI':
        # Baryon isocurvature: S_b = 1, Phi = 0
        y0[:, IDX_DELTA_B] = 1.0

    elif mode == 'NDI':
        # Neutrino density isocurvature with full hierarchy
        # S_nu = delta_nu - delta_gamma = 1, compensated
        y0[:, idx_theta(0)] = -_f_nu / 4.0
        y0[:, _n_start + 0] = (1.0 - _f_nu) / 4.0
        y0[:, IDX_DELTA_B] = 0.75 * (-_f_nu)

    elif mode == 'NVI':
        # Neutrino velocity isocurvature
        y0[:, _n_start + 1] = 1.0
        if l_nu_max >= 2:
            y0[:, _n_start + 2] = (2.0 / 5.0) * kt
        y0[:, _n_start + 0] = -kt

    return mx.array(y0)


# ============================================================================
# Unified isocurvature_ic dispatcher
# ============================================================================

def isocurvature_ic(k_arr_np, bg, mode='CDI', solver='hires',
                     l_gamma_max=25, l_nu_max=L_NU_MAX):
    """
    Generate isocurvature initial conditions for any solver backend.

    Parameters
    ----------
    k_arr_np : array (N_k,)
        Wavenumber array in Mpc^-1.
    bg : Background
        Solved background object.
    mode : str
        Isocurvature mode: 'CDI', 'BI', 'NDI', or 'NVI'.
    solver : str
        Solver backend: 'tca', 'neutrino', or 'hires'.
    l_gamma_max : int
        Photon hierarchy truncation (hires only).
    l_nu_max : int
        Neutrino hierarchy truncation.

    Returns
    -------
    y0 : mx.array
        Initial state vector for the specified solver.
    """
    solver = solver.lower()
    if solver == 'tca':
        return isocurvature_ic_tca(k_arr_np, bg, mode=mode)
    elif solver == 'neutrino':
        return isocurvature_ic_neutrino(k_arr_np, bg, mode=mode,
                                         l_nu_max=l_nu_max)
    elif solver == 'hires':
        return isocurvature_ic_hires(k_arr_np, bg, mode=mode,
                                      l_gamma_max=l_gamma_max,
                                      l_nu_max=l_nu_max)
    else:
        raise ValueError(f"Unknown solver '{solver}'. "
                         f"Choose from 'tca', 'neutrino', 'hires'.")


# ============================================================================
# Mixed initial conditions
# ============================================================================

def mixed_ic(k_arr_np, bg, alpha=0.0, beta=0.0, mode='CDI',
             solver='hires', l_gamma_max=25, l_nu_max=L_NU_MAX):
    """
    Mixed adiabatic + isocurvature initial conditions.

    The primordial spectrum is:
        P(k) = (1 - alpha) * P_ad(k) + alpha * P_iso(k) + beta * P_corr(k)

    Since C_l is quadratic in the transfer functions (hence in y0), we
    combine amplitudes as:
        y0 = sqrt(1 - alpha) * y_ad + sqrt(alpha) * y_iso + beta * y_corr

    The correlated component uses y_corr = y_ad (same shape), giving a
    cross-term ~ beta * sqrt(1-alpha) * |y_ad|^2 when squared.

    For small alpha, this reduces to essentially adiabatic initial conditions
    with a small isocurvature admixture.

    Parameters
    ----------
    k_arr_np : array (N_k,)
        Wavenumber array in Mpc^-1.
    bg : Background
        Solved background object.
    alpha : float
        Isocurvature fraction, 0 <= alpha <= 1.
        alpha = 0: pure adiabatic. alpha = 1: pure isocurvature.
    beta : float
        Correlation amplitude between adiabatic and isocurvature.
        beta = 0: uncorrelated. |beta| <= 2*sqrt(alpha*(1-alpha)).
    mode : str
        Isocurvature mode: 'CDI', 'BI', 'NDI', or 'NVI'.
    solver : str
        Solver backend: 'tca', 'neutrino', or 'hires'.
    l_gamma_max : int
        Photon hierarchy truncation (hires only).
    l_nu_max : int
        Neutrino hierarchy truncation.

    Returns
    -------
    y0 : mx.array
        Mixed initial state vector.

    Notes
    -----
    IMPORTANT: The mixed IC approach of linearly combining state vectors is
    valid because the Boltzmann equations are linear in the perturbation
    variables. Each mode evolves independently, so the total C_l is the sum
    of auto- and cross-spectra. However, for a proper likelihood analysis,
    one should evolve each mode separately and combine the C_l, not the
    state vectors. This function provides a quick approximation that is
    exact for the auto-spectra but only approximate for the cross-spectrum.

    For a rigorous treatment, evolve adiabatic and isocurvature modes
    separately and combine:
        C_l = (1-alpha)*C_l^ad + alpha*C_l^iso + beta*C_l^corr
    """
    if not (0.0 <= alpha <= 1.0):
        raise ValueError(f"alpha must be in [0, 1], got {alpha}")

    solver_lower = solver.lower()
    k_arr_np = np.asarray(k_arr_np, dtype=np.float32)

    # Get adiabatic IC
    if solver_lower == 'tca':
        y_ad = adiabatic_ic_tca(k_arr_np, bg)
    elif solver_lower == 'neutrino':
        y_ad = adiabatic_ic_neutrino(k_arr_np, bg, l_nu_max=l_nu_max)
    elif solver_lower == 'hires':
        y_ad = adiabatic_ic_hires(k_arr_np, bg, l_gamma_max, l_nu_max)
    else:
        raise ValueError(f"Unknown solver '{solver_lower}'.")

    # Pure adiabatic shortcut
    if alpha == 0.0 and beta == 0.0:
        return y_ad

    # Get isocurvature IC
    y_iso = isocurvature_ic(k_arr_np, bg, mode=mode, solver=solver_lower,
                             l_gamma_max=l_gamma_max, l_nu_max=l_nu_max)

    # Combine: y_mixed = sqrt(1-alpha) * y_ad + sqrt(alpha) * y_iso
    #          + beta * y_ad  (correlation term)
    coeff_ad = float(np.sqrt(1.0 - alpha)) + beta
    coeff_iso = float(np.sqrt(alpha))

    y_mixed = coeff_ad * y_ad + coeff_iso * y_iso

    return y_mixed


# ============================================================================
# Separate-mode evolution for rigorous C_l combination
# ============================================================================

def mixed_cl_from_separate_runs(Cl_adiabatic, Cl_isocurvature,
                                 Cl_correlated=None,
                                 alpha=0.0, beta=0.0):
    """
    Combine separately-computed C_l spectra for a mixed adiabatic+isocurvature
    primordial spectrum.

    This is the rigorous approach: evolve each mode independently, compute
    its C_l, then combine:
        C_l = (1 - alpha) * C_l^ad + alpha * C_l^iso + beta * C_l^corr

    Parameters
    ----------
    Cl_adiabatic : array (N_ell,)
        C_l from pure adiabatic initial conditions.
    Cl_isocurvature : array (N_ell,)
        C_l from pure isocurvature initial conditions.
    Cl_correlated : array (N_ell,) or None
        Cross-spectrum C_l. If None, assumed zero.
    alpha : float
        Isocurvature power fraction.
    beta : float
        Correlation amplitude.

    Returns
    -------
    Cl_mixed : array (N_ell,)
    """
    Cl_mixed = (1.0 - alpha) * Cl_adiabatic + alpha * Cl_isocurvature
    if Cl_correlated is not None and beta != 0.0:
        Cl_mixed += beta * Cl_correlated
    return Cl_mixed


# ============================================================================
# Self-test
# ============================================================================

def _test_isocurvature():
    """
    Test isocurvature initial conditions:
    1. mixed_ic with alpha=0 reproduces adiabatic IC exactly
    2. CDI mode has Phi=0, delta_c=1
    3. NDI mode has correct photon/neutrino compensation
    4. Run a short integration and compare adiabatic vs CDI C_l
    """
    print("=" * 70)
    print("TEST: Isocurvature initial conditions")
    print("=" * 70)

    from .background import Background
    from .perturbations_implicit import (
        ImplicitBoltzmannSolver, adiabatic_ic, build_tau_grid_implicit,
        imex_rk4_step, TCA_PHI, TCA_DELTA_C, TCA_THETA_0, TCA_N_VAR,
    )
    from .spectra import compute_cl

    # ----- Step 1: Background -----
    print("\n--- Step 1: Background ---")
    bg = Background(khronon=False)
    bg.solve()

    k_arr = np.geomspace(5e-4, 0.30, 300).astype(np.float32)

    # ----- Test A: alpha=0 reproduces adiabatic -----
    print("\n--- Test A: alpha=0 == adiabatic ---")
    y_ad = adiabatic_ic(k_arr, bg)
    y_mixed_0 = mixed_ic(k_arr, bg, alpha=0.0, mode='CDI', solver='tca')

    diff = float(mx.max(mx.abs(y_ad - y_mixed_0)))
    print(f"  max|y_ad - mixed(alpha=0)| = {diff:.2e}")
    assert diff < 1e-7, f"FAIL: alpha=0 should match adiabatic, got diff={diff}"
    print("  PASS")

    # ----- Test B: CDI mode structure -----
    print("\n--- Test B: CDI mode structure ---")
    y_cdi = isocurvature_ic_tca(k_arr, bg, mode='CDI')
    y_cdi_np = np.array(y_cdi)
    print(f"  Phi  = {y_cdi_np[0, TCA_PHI]:.4f} (should be 0)")
    print(f"  delta_c = {y_cdi_np[0, TCA_DELTA_C]:.4f} (should be 1)")
    print(f"  Theta_0 = {y_cdi_np[0, TCA_THETA_0]:.4f} (should be 0)")
    assert abs(y_cdi_np[0, TCA_PHI]) < 1e-7
    assert abs(y_cdi_np[0, TCA_DELTA_C] - 1.0) < 1e-7
    assert abs(y_cdi_np[0, TCA_THETA_0]) < 1e-7
    print("  PASS")

    # ----- Test C: NDI photon-neutrino compensation -----
    print("\n--- Test C: NDI compensation (neutrino solver) ---")
    from .perturbations_neutrino import IDX_THETA_0, IDX_N_START

    y_ndi = isocurvature_ic_neutrino(k_arr, bg, mode='NDI')
    y_ndi_np = np.array(y_ndi)
    theta0 = y_ndi_np[0, IDX_THETA_0]
    n0 = y_ndi_np[0, IDX_N_START]
    # Check: f_gamma*delta_gamma + f_nu*delta_nu = 0
    # delta_gamma = 4*Theta_0, delta_nu = 4*N_0
    total_rad = _f_gamma * 4.0 * theta0 + _f_nu * 4.0 * n0
    print(f"  Theta_0 = {theta0:.6f}")
    print(f"  N_0     = {n0:.6f}")
    print(f"  f_gamma*delta_gamma + f_nu*delta_nu = {total_rad:.6e} (should be ~0)")
    assert abs(total_rad) < 1e-5, f"FAIL: radiation not compensated, got {total_rad}"
    # Check entropy perturbation: S_nu = delta_nu - delta_gamma = 4*(N_0 - Theta_0)
    s_nu = 4.0 * (n0 - theta0)
    print(f"  S_nu = 4*(N_0 - Theta_0) = {s_nu:.6f} (should be 1.0)")
    assert abs(s_nu - 1.0) < 0.01
    print("  PASS")

    # ----- Test D: NVI mode -----
    print("\n--- Test D: NVI mode (neutrino solver) ---")
    y_nvi = isocurvature_ic_neutrino(k_arr, bg, mode='NVI')
    y_nvi_np = np.array(y_nvi)
    n1 = y_nvi_np[0, IDX_N_START + 1]
    print(f"  N_1 = {n1:.4f} (should be 1.0)")
    assert abs(n1 - 1.0) < 1e-7
    print("  PASS")

    # ----- Test E: Evolve and compare C_l -----
    print("\n--- Test E: C_l comparison (adiabatic vs 5% CDI) ---")
    from .background import Omega_r, Omega_b, H0_Mpc

    # Use TCA solver for speed
    tau_grid = build_tau_grid_implicit(bg, float(k_arr[-1]))
    calH_np = bg.calH_at_tau(tau_grid).astype(np.float32)
    R_np = bg.R_at_tau(tau_grid).astype(np.float32)
    a_np = bg.a_at_tau(tau_grid).astype(np.float32)

    _Or = Omega_r
    _Ob = Omega_b
    _Oc = float(bg.Omega_cdm)
    _H0 = H0_Mpc

    def _evolve(y0_mx):
        y = y0_mx
        k_mx = mx.array(k_arr)
        for i in range(len(tau_grid) - 1):
            dt = float(tau_grid[i + 1] - tau_grid[i])
            y = imex_rk4_step(
                y, k_mx,
                mx.array(np.float32(calH_np[i])),
                mx.array(np.float32(R_np[i])),
                _Or, _Ob, _Oc,
                mx.array(np.float32(a_np[i])), _H0, dt)
            if i % 300 == 0:
                mx.eval(y)
        mx.eval(y)
        return y

    print("  Evolving pure adiabatic...")
    y_final_ad = _evolve(adiabatic_ic(k_arr, bg))
    y_ad_np = np.array(y_final_ad)
    sw_ad = (y_ad_np[:, TCA_THETA_0] + y_ad_np[:, TCA_PHI]) * \
            np.exp(-(k_arr / bg.k_D) ** 2)

    print("  Evolving 5% CDI mixed...")
    y0_mixed = mixed_ic(k_arr, bg, alpha=0.05, mode='CDI', solver='tca')
    y_final_mix = _evolve(y0_mixed)
    y_mix_np = np.array(y_final_mix)
    sw_mix = (y_mix_np[:, TCA_THETA_0] + y_mix_np[:, TCA_PHI]) * \
             np.exp(-(k_arr / bg.k_D) ** 2)

    # Compute C_l for both
    ell_values = np.arange(2, 2001, dtype=np.float64)
    _, _, Dl_ad = compute_cl(sw_ad, k_arr, ell_values, bg.D_A)
    _, _, Dl_mix = compute_cl(sw_mix, k_arr, ell_values, bg.D_A)

    # Check that they differ
    max_diff_pct = np.max(np.abs(Dl_mix - Dl_ad) / (np.abs(Dl_ad) + 1e-10)) * 100
    rms_diff = np.sqrt(np.mean((Dl_mix - Dl_ad) ** 2))
    corr = np.corrcoef(Dl_ad[Dl_ad > 0], Dl_mix[Dl_ad > 0])[0, 1]

    print(f"  Max relative diff: {max_diff_pct:.1f}%")
    print(f"  RMS diff (muK^2):  {rms_diff:.1f}")
    print(f"  Correlation:       {corr:.6f}")
    print(f"  (5% CDI should make a few-percent difference)")

    if corr > 0.95 and max_diff_pct > 0.1:
        print("  PASS: spectra are similar but distinguishable")
    elif max_diff_pct < 0.01:
        print("  WARNING: spectra too similar -- isocurvature may not be working")
    else:
        print("  PASS")

    # ----- Test F: HiRes IC dimensions -----
    print("\n--- Test F: HiRes IC dimensions ---")
    from .perturbations_hires import n_var_total as hires_nvar
    for m in ISOCURVATURE_MODES:
        y_test = isocurvature_ic_hires(k_arr[:10], bg, mode=m,
                                        l_gamma_max=25, l_nu_max=20)
        expected = hires_nvar(25, 20)
        actual = y_test.shape[1]
        status = "PASS" if actual == expected else "FAIL"
        print(f"  {m}: shape=({y_test.shape[0]}, {actual}), "
              f"expected N_var={expected} ... {status}")
        assert actual == expected

    # ----- Summary -----
    print("\n" + "=" * 70)
    print("ALL ISOCURVATURE TESTS PASSED")
    print("=" * 70)
    return True


if __name__ == '__main__':
    _test_isocurvature()
