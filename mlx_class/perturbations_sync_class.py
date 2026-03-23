"""
perturbations_sync_class.py -- CLASS-style three-phase Boltzmann solver.

Implements the three approximation regimes from CLASS (perturbations.c):

Phase 1 -- TCA (Tight Coupling Approximation):
    During tight coupling (|kappa_dot| >> k, calH), the photon-baryon plasma
    is nearly perfect fluid.  CLASS does NOT evolve the full photon hierarchy.
    Instead it evolves only {delta_g, theta_b, delta_b, eta, delta_c} plus
    the neutrino hierarchy, using analytical formulas for photon shear (sigma_g)
    and the photon-baryon slip (theta_g - theta_b).

Phase 2 -- Full Boltzmann hierarchy:
    After the TCA switch, the full photon temperature + polarization +
    neutrino hierarchies are evolved.  This is the same as perturbations_sync.py
    but with polarization feedback (Pi = F_g2 + E_0 + E_2).

Phase 3 -- RSA (Radiation Streaming Approximation):
    When k*tau > 45 and tau > tau_free_streaming, the photon multipoles
    oscillate too rapidly to track.  CLASS replaces them with algebraic
    approximations from the metric:
        delta_g = (4/k^2) * (calH*h' - k^2*eta)     [sync gauge]
        theta_g = -(1/2) * h'                         [sync gauge]
    Neutrinos get the same treatment (no Thomson terms).
    Only {eta, delta_c, delta_b, theta_b} are evolved.

References:
    Ma & Bertschinger (1995) ApJ 455, 7
    Blas, Lesgourgues & Tram (2011) JCAP 07, 034 -- compromise_CLASS TCA
    CLASS perturbations.c (lines 6064-6166, 8919-9060, 10181-10294)

Author: Sheng-Kai Huang, 2026
"""

import numpy as np
from scipy.integrate import solve_ivp

from .background import (
    Omega_r as _OMEGA_R,
    Omega_b as _OMEGA_B,
    Omega_c as _OMEGA_C,
    Omega_m as _OMEGA_M,
    H0_Mpc as _H0_MPC,
    T_CMB as _T_CMB,
    m_H_SI as _m_H_SI,
    Y_He as _Y_He,
    c_SI as _c_SI,
)
from .perturbations_neutrino import (
    _f_nu, _OMEGA_GAMMA, _OMEGA_NU,
)


# ============================================================================
# Physical constants
# ============================================================================
_k_B_SI = 1.380649e-23       # J/K
_mu_mol = 1.0 / (1.0 - 0.75 * _Y_He)   # mean molecular weight


# ============================================================================
# Default hierarchy truncation (CLASS defaults)
# ============================================================================
L_GAMMA_MAX = 12    # photon temperature (CLASS default l_max_g)
L_POL_MAX = 10      # photon polarization (CLASS default l_max_pol_g)
L_NU_MAX = 17       # neutrino (CLASS default l_max_ur)


# ============================================================================
# Baryon sound speed
# ============================================================================

def baryon_cs2(a):
    """
    Baryon sound speed squared c_s^2 = k_B T_b / (mu m_H c^2).
    During tight coupling T_b ~ T_CMB / a.
    """
    T_b = _T_CMB / a
    cs2 = _k_B_SI * T_b / (_mu_mol * _m_H_SI * _c_SI**2)
    return cs2


# ============================================================================
# Phase 1: TCA state vector layout
# ============================================================================
# During TCA, photon hierarchy is truncated to {delta_g, theta_g} (no shear
# or higher multipoles).  Following CLASS (Section 4.1-4.3 of the strategy doc):
#   [eta, delta_c, delta_b, theta_b, delta_g, theta_g,
#    F_n0, F_n1, ..., F_n{ln_max}]
#
# CLASS evolves BOTH delta_g AND theta_g during TCA.
# sigma_g = (16/45)*tau_c*(theta_g + metric_shear) is analytical.
# The slip theta_g - theta_b is maintained by the coupled theta_b/theta_g eqs.

TCA_IDX_ETA = 0
TCA_IDX_DELTA_C = 1
TCA_IDX_DELTA_B = 2
TCA_IDX_THETA_B = 3
TCA_IDX_DELTA_G = 4
TCA_IDX_THETA_G = 5
TCA_IDX_FN_START = 6


def n_var_tca(ln_max):
    """Number of TCA-phase variables: 6 + (ln_max + 1) neutrino multipoles."""
    return TCA_IDX_FN_START + ln_max + 1


# ============================================================================
# Phase 2: Full hierarchy state vector layout
# ============================================================================
FULL_IDX_ETA = 0
FULL_IDX_DELTA_C = 1
FULL_IDX_DELTA_B = 2
FULL_IDX_THETA_B = 3
FULL_IDX_FG_START = 4


def full_idx_fg(l):
    """Index of photon multipole F_gamma,l in full state vector."""
    return FULL_IDX_FG_START + l


def full_idx_epol_start(lg_max):
    """Starting index of E-mode polarization hierarchy."""
    return FULL_IDX_FG_START + lg_max + 1


def full_idx_epol(l, lg_max):
    """Index of polarization multipole E_l."""
    return full_idx_epol_start(lg_max) + l


def full_idx_fn_start(lg_max, lp_max):
    """Starting index of neutrino hierarchy in full state vector."""
    return full_idx_epol_start(lg_max) + lp_max + 1


def full_idx_fn(l, lg_max, lp_max):
    """Index of neutrino multipole F_nu,l."""
    return full_idx_fn_start(lg_max, lp_max) + l


def n_var_full(lg_max, lp_max, ln_max):
    """Total number of full-phase variables."""
    return FULL_IDX_FG_START + (lg_max + 1) + (lp_max + 1) + (ln_max + 1)


# ============================================================================
# Phase 3: RSA state vector layout
# ============================================================================
# Only evolve {eta, delta_c, delta_b, theta_b}.
# Photon and neutrino quantities are diagnosed algebraically.

RSA_IDX_ETA = 0
RSA_IDX_DELTA_C = 1
RSA_IDX_DELTA_B = 2
RSA_IDX_THETA_B = 3
N_VAR_RSA = 4


# ============================================================================
# Switch criteria (CLASS defaults)
# ============================================================================

def find_tca_switch_time(bg, k,
                         tau_c_over_tau_h_thr=0.015,
                         tau_c_over_tau_k_thr=0.01):
    """
    Find the conformal time at which TCA switches off.

    CLASS criterion: TCA is ON when BOTH
        tau_c / tau_h < 0.015   AND   tau_c / tau_k < 0.01
    where tau_c = 1/|kd|, tau_h = 1/calH, tau_k = 1/k.

    Equivalently: calH/|kd| < 0.015  AND  k/|kd| < 0.01.
    """
    abs_kd = np.abs(bg.kappa_dot_grid)
    calH = bg.calH_grid
    tau_arr = bg.tau_grid

    ratio_h = calH / np.maximum(abs_kd, 1e-30)
    ratio_k = k / np.maximum(abs_kd, 1e-30)

    tca_off = (ratio_h >= tau_c_over_tau_h_thr) | (ratio_k >= tau_c_over_tau_k_thr)
    switch_idx = np.argmax(tca_off)

    if switch_idx == 0 and not tca_off[0]:
        return float(tau_arr[-1])

    tau_switch = float(tau_arr[switch_idx])

    # Safety: don't switch too early
    tau_init = float(tau_arr[1])
    tau_switch = max(tau_switch, tau_init * 5)

    return tau_switch


def find_rsa_switch_time(bg, k, k_tau_threshold=45.0):
    """
    Find the conformal time at which RSA switches on.

    CLASS criterion: RSA is ON when ALL of:
        k * tau > 45                (mode deep inside horizon)
        tau > tau_free_streaming    (after recombination)

    tau_free_streaming ~ tau at which visibility function has fallen
    to a small fraction of peak.  We approximate it as the time when
    kappa < 1 (optically thin).
    """
    tau_arr = bg.tau_grid
    kappa = bg.kappa_grid

    # Condition 1: k * tau > threshold
    k_tau = k * tau_arr
    cond_ktau = k_tau > k_tau_threshold

    # Condition 2: optically thin (past free-streaming surface)
    # Use kappa < 1 as a proxy for tau_free_streaming
    cond_fs = kappa < 1.0

    both = cond_ktau & cond_fs
    switch_idx = np.argmax(both)

    if switch_idx == 0 and not both[0]:
        return float(tau_arr[-1])  # RSA never activates

    return float(tau_arr[switch_idx])


# ============================================================================
# TCA analytical shear and slip
# ============================================================================

def compute_tca_shear(theta_g, h_prime, eta_prime, abs_kd):
    """
    First-order TCA photon shear (CLASS line 10042).

    sigma_g = (16/45) * tau_c * (theta_g + metric_shear)
    metric_shear = (h' + 6*eta') / 2   [sync gauge]

    Returns (sigma_g, F_gamma_2 = 2*sigma_g).
    """
    tau_c = 1.0 / abs_kd if abs_kd > 0 else 0.0
    metric_shear = 0.5 * (h_prime + 6.0 * eta_prime)
    sigma_g = (16.0 / 45.0) * tau_c * (theta_g + metric_shear)
    Fg2 = 2.0 * sigma_g
    return sigma_g, Fg2


def compute_tca_slip(k, calH, R, delta_g, theta_b, theta_g,
                     sigma_g, abs_kd, h_prime, eta_prime, a, cs2_b):
    """
    First-order photon-baryon slip (compromise_CLASS).

    slip = tau_c / (1+R) * [
        -2*calH*R/(1+R)*theta_b
        + k^2*(delta_g/4 - sigma_g)
        + k^2*cs2_b*delta_b   -- small, dropped for simplicity
    ]

    For sync gauge, metric_euler = 0. The alpha-dependent correction is a
    second-order effect in the compromise_CLASS formula; we include the
    leading first-order slip.
    """
    if abs_kd <= 0:
        return 0.0

    tau_c = 1.0 / abs_kd
    k2 = k * k

    # First-order slip (CLASS "first_order_CLASS" variant)
    term1 = -2.0 * calH * R / (1.0 + R) * theta_b
    term2 = k2 * (delta_g / 4.0 - sigma_g)

    slip = tau_c / (1.0 + R) * (term1 + term2)

    return slip


# ============================================================================
# Phase 1: TCA RHS
# ============================================================================

def make_tca_rhs(k, bg, ln_max=L_NU_MAX):
    """
    Build the RHS function for Phase 1 (TCA).

    Evolved variables: [eta, delta_c, delta_b, theta_b, delta_g, theta_g,
                        F_n0, ..., F_n{ln_max}]

    sigma_g is computed analytically at each step.

    Baryon velocity (CLASS line 8956):
        theta_b' = [-calH*theta_b + k^2*cs2_b*delta_b
                     + k^2*R*(delta_g/4 - sigma_g)] / (1+R)

    Photon velocity (CLASS line 9050):
        theta_g' = -(theta_b' + calH*theta_b - k^2*cs2_b*delta_b) / R
                   + k^2*(delta_g/4 - sigma_g)

    In synchronous gauge, metric_euler = 0.
    """
    k2 = k * k
    H02 = _H0_MPC ** 2
    Og, On, Ob, Oc = _OMEGA_GAMMA, _OMEGA_NU, _OMEGA_B, _OMEGA_C
    nvar = n_var_tca(ln_max)
    fn_s = TCA_IDX_FN_START

    # Pre-tabulate background
    _tau = bg.tau_grid.copy()
    _calH = bg.calH_grid.copy()
    _a = bg.a_grid.copy()
    _R = bg.R_grid.copy()
    _kd = bg.kappa_dot_grid.copy()

    def rhs(tau, y):
        calH = np.interp(tau, _tau, _calH)
        a = np.interp(tau, _tau, _a)
        R = np.interp(tau, _tau, _R)
        kappa_dot = np.interp(tau, _tau, _kd)
        abs_kd = abs(kappa_dot)
        tau_safe = max(tau, 1e-10)

        ia = 1.0 / a
        ia2 = ia * ia

        # --- Extract state ---
        eta = y[TCA_IDX_ETA]
        delta_c = y[TCA_IDX_DELTA_C]
        delta_b = y[TCA_IDX_DELTA_B]
        theta_b = y[TCA_IDX_THETA_B]
        delta_g = y[TCA_IDX_DELTA_G]
        theta_g = y[TCA_IDX_THETA_G]
        Fn = y[fn_s: fn_s + ln_max + 1]

        delta_n = Fn[0]
        theta_n = 0.75 * k * Fn[1]

        # h' from 00-constraint
        src00 = (Og * ia2 * delta_g + On * ia2 * delta_n
                 + Ob * ia * delta_b + Oc * ia * delta_c)
        h_prime = (2.0 / calH) * (k2 * eta + 1.5 * H02 * src00)

        # eta' from 0i-constraint
        src0i = ((4.0 / 3.0) * Og * ia2 * theta_g
                 + (4.0 / 3.0) * On * ia2 * theta_n
                 + Ob * ia * theta_b)
        eta_prime = 1.5 * H02 / k2 * src0i

        # TCA analytical shear
        sigma_g, Fg2 = compute_tca_shear(theta_g, h_prime, eta_prime, abs_kd)

        # Baryon sound speed
        cs2_b = baryon_cs2(a)

        # --- CDM ---
        d_delta_c = -0.5 * h_prime

        # --- Baryons ---
        d_delta_b = -theta_b - 0.5 * h_prime

        # Baryon velocity: same as full hierarchy, with Thomson coupling
        # theta_b' = -calH*theta_b + cs2_b*k^2*delta_b + |kd|/R*(theta_g - theta_b)
        d_theta_b = (-calH * theta_b
                     + cs2_b * k2 * delta_b
                     + abs_kd / R * (theta_g - theta_b))

        # --- Photon density ---
        # delta_g' = -(4/3)*theta_g - (2/3)*h'
        d_delta_g = -(4.0 / 3.0) * theta_g - (2.0 / 3.0) * h_prime

        # --- Photon velocity ---
        # Same as full hierarchy l=1, but with sigma_g analytical:
        # theta_g' = k^2*(delta_g/4 - sigma_g) + |kd|*(theta_b - theta_g)
        # (metric_euler = 0 in sync gauge)
        d_theta_g = (k2 * (delta_g / 4.0 - sigma_g)
                     + abs_kd * (theta_b - theta_g))

        # --- Neutrino hierarchy (no collisions) ---
        dFn = np.zeros(ln_max + 1)
        dFn[0] = -k * Fn[1] - (2.0 / 3.0) * h_prime

        F2n = Fn[2] if ln_max >= 2 else 0.0
        dFn[1] = (k / 3.0) * (Fn[0] - 2.0 * F2n)

        if ln_max >= 2:
            F3n = Fn[3] if ln_max >= 3 else 0.0
            dFn[2] = ((k / 5.0) * (2.0 * Fn[1] - 3.0 * F3n)
                      + (4.0 / 15.0) * h_prime + (8.0 / 5.0) * eta_prime)

        for ell in range(3, ln_max):
            dFn[ell] = (k / (2.0 * ell + 1.0)
                        * (ell * Fn[ell - 1] - (ell + 1) * Fn[ell + 1]))

        if ln_max >= 1:
            dFn[ln_max] = (k * Fn[ln_max - 1] * ln_max / (2.0 * ln_max + 1.0)
                           - (ln_max + 1.0) / tau_safe * Fn[ln_max])

        # --- Assemble ---
        dydt = np.zeros(nvar)
        dydt[TCA_IDX_ETA] = eta_prime
        dydt[TCA_IDX_DELTA_C] = d_delta_c
        dydt[TCA_IDX_DELTA_B] = d_delta_b
        dydt[TCA_IDX_THETA_B] = d_theta_b
        dydt[TCA_IDX_DELTA_G] = d_delta_g
        dydt[TCA_IDX_THETA_G] = d_theta_g
        dydt[fn_s: fn_s + ln_max + 1] = dFn

        return dydt

    return rhs, nvar


# ============================================================================
# Phase 2: Full hierarchy RHS
# ============================================================================

def make_full_rhs(k, bg, lg_max=L_GAMMA_MAX, lp_max=L_POL_MAX,
                  ln_max=L_NU_MAX):
    """
    Build the RHS for Phase 2 (full Boltzmann hierarchy with polarization).

    This is essentially the same as perturbations_sync_tca.make_sync_rhs_tca
    with pol_approx='full', but with the corrected eta' coefficient (8/5).
    """
    k2 = k * k
    H02 = _H0_MPC ** 2
    Og, On, Ob, Oc = _OMEGA_GAMMA, _OMEGA_NU, _OMEGA_B, _OMEGA_C
    nvar = n_var_full(lg_max, lp_max, ln_max)
    ep_s = full_idx_epol_start(lg_max)
    fn_s = full_idx_fn_start(lg_max, lp_max)

    _tau = bg.tau_grid.copy()
    _calH = bg.calH_grid.copy()
    _a = bg.a_grid.copy()
    _R = bg.R_grid.copy()
    _kd = bg.kappa_dot_grid.copy()

    def rhs(tau, y):
        calH = np.interp(tau, _tau, _calH)
        a = np.interp(tau, _tau, _a)
        R = np.interp(tau, _tau, _R)
        kappa_dot = np.interp(tau, _tau, _kd)
        abs_kd = abs(kappa_dot)
        tau_safe = max(tau, 1e-10)

        ia = 1.0 / a
        ia2 = ia * ia

        # --- Extract state ---
        eta = y[FULL_IDX_ETA]
        delta_c = y[FULL_IDX_DELTA_C]
        delta_b = y[FULL_IDX_DELTA_B]
        theta_b = y[FULL_IDX_THETA_B]

        Fg = y[FULL_IDX_FG_START: FULL_IDX_FG_START + lg_max + 1]
        E = y[ep_s: ep_s + lp_max + 1]
        Fn = y[fn_s: fn_s + ln_max + 1]

        delta_g = Fg[0]
        delta_n = Fn[0]
        theta_g = 0.75 * k * Fg[1]
        theta_n = 0.75 * k * Fn[1]

        # --- Diagnose h' ---
        src00 = (Og * ia2 * delta_g + On * ia2 * delta_n
                 + Ob * ia * delta_b + Oc * ia * delta_c)
        h_prime = (2.0 / calH) * (k2 * eta + 1.5 * H02 * src00)

        # --- Diagnose eta' ---
        src0i = ((4.0 / 3.0) * Og * ia2 * theta_g
                 + (4.0 / 3.0) * On * ia2 * theta_n
                 + Ob * ia * theta_b)
        eta_prime = 1.5 * H02 / k2 * src0i

        # --- CDM ---
        d_delta_c = -0.5 * h_prime

        # --- Baryons ---
        d_delta_b = -theta_b - 0.5 * h_prime
        cs2_b = baryon_cs2(a)
        d_theta_b = (-calH * theta_b
                     + cs2_b * k2 * delta_b
                     + abs_kd / R * (theta_g - theta_b))

        # --- Polarization source Pi ---
        Fg2 = Fg[2] if lg_max >= 2 else 0.0
        E0 = E[0] if lp_max >= 0 else 0.0
        E2 = E[2] if lp_max >= 2 else 0.0
        Pi = Fg2 + E0 + E2

        # --- Photon temperature hierarchy ---
        dFg = np.zeros(lg_max + 1)

        # l=0
        dFg[0] = -k * Fg[1] - (2.0 / 3.0) * h_prime

        # l=1
        F2g = Fg[2] if lg_max >= 2 else 0.0
        dFg[1] = ((k / 3.0) * (Fg[0] - 2.0 * F2g)
                  + abs_kd * (-Fg[1] + 4.0 * theta_b / (3.0 * k)))

        # l=2 with polarization feedback
        if lg_max >= 2:
            F3g = Fg[3] if lg_max >= 3 else 0.0
            dFg[2] = ((k / 5.0) * (2.0 * Fg[1] - 3.0 * F3g)
                      + (4.0 / 15.0) * h_prime + (8.0 / 5.0) * eta_prime
                      - abs_kd * (Fg[2] - Pi / 10.0))

        # l=3..lg_max-1
        for ell in range(3, lg_max):
            dFg[ell] = (k / (2.0 * ell + 1.0)
                        * (ell * Fg[ell - 1] - (ell + 1) * Fg[ell + 1])
                        - abs_kd * Fg[ell])

        # l=lg_max: truncation
        if lg_max >= 3:
            dFg[lg_max] = (k * Fg[lg_max - 1] * lg_max / (2.0 * lg_max + 1.0)
                           - (lg_max + 1.0) / tau_safe * Fg[lg_max]
                           - abs_kd * Fg[lg_max])

        # --- E-mode polarization hierarchy ---
        dE = np.zeros(lp_max + 1)

        if lp_max >= 0:
            E1 = E[1] if lp_max >= 1 else 0.0
            dE[0] = -k * E1 - abs_kd * E[0] + 0.5 * abs_kd * Pi

        if lp_max >= 1:
            E0_val = E[0]
            E2_val = E[2] if lp_max >= 2 else 0.0
            dE[1] = (k / 3.0) * (E0_val - 2.0 * E2_val) - abs_kd * E[1]

        if lp_max >= 2:
            E1_val = E[1]
            E3_val = E[3] if lp_max >= 3 else 0.0
            dE[2] = ((k / 5.0) * (2.0 * E1_val - 3.0 * E3_val)
                     - abs_kd * E[2] + 0.1 * abs_kd * Pi)

        for ell in range(3, lp_max):
            dE[ell] = (k / (2.0 * ell + 1.0)
                       * (ell * E[ell - 1] - (ell + 1) * E[ell + 1])
                       - abs_kd * E[ell])

        if lp_max >= 3:
            dE[lp_max] = (k * E[lp_max - 1] * lp_max / (2.0 * lp_max + 1.0)
                          - (lp_max + 1.0) / tau_safe * E[lp_max]
                          - abs_kd * E[lp_max])

        # --- Neutrino hierarchy ---
        dFn = np.zeros(ln_max + 1)
        dFn[0] = -k * Fn[1] - (2.0 / 3.0) * h_prime

        F2n = Fn[2] if ln_max >= 2 else 0.0
        dFn[1] = (k / 3.0) * (Fn[0] - 2.0 * F2n)

        if ln_max >= 2:
            F3n = Fn[3] if ln_max >= 3 else 0.0
            dFn[2] = ((k / 5.0) * (2.0 * Fn[1] - 3.0 * F3n)
                      + (4.0 / 15.0) * h_prime + (8.0 / 5.0) * eta_prime)

        for ell in range(3, ln_max):
            dFn[ell] = (k / (2.0 * ell + 1.0)
                        * (ell * Fn[ell - 1] - (ell + 1) * Fn[ell + 1]))

        if ln_max >= 1:
            dFn[ln_max] = (k * Fn[ln_max - 1] * ln_max / (2.0 * ln_max + 1.0)
                           - (ln_max + 1.0) / tau_safe * Fn[ln_max])

        # --- Assemble ---
        dydt = np.zeros(nvar)
        dydt[FULL_IDX_ETA] = eta_prime
        dydt[FULL_IDX_DELTA_C] = d_delta_c
        dydt[FULL_IDX_DELTA_B] = d_delta_b
        dydt[FULL_IDX_THETA_B] = d_theta_b
        dydt[FULL_IDX_FG_START: FULL_IDX_FG_START + lg_max + 1] = dFg
        dydt[ep_s: ep_s + lp_max + 1] = dE
        dydt[fn_s: fn_s + ln_max + 1] = dFn

        return dydt

    return rhs, nvar


# ============================================================================
# Phase 3: RSA RHS
# ============================================================================

def make_rsa_rhs(k, bg):
    """
    Build the RHS for Phase 3 (Radiation Streaming Approximation).

    Only evolves {eta, delta_c, delta_b, theta_b}.
    Photon/neutrino contributions to the Einstein equations are computed
    algebraically from the metric.
    """
    k2 = k * k
    H02 = _H0_MPC ** 2
    Og, On, Ob, Oc = _OMEGA_GAMMA, _OMEGA_NU, _OMEGA_B, _OMEGA_C

    _tau = bg.tau_grid.copy()
    _calH = bg.calH_grid.copy()
    _a = bg.a_grid.copy()
    _R = bg.R_grid.copy()
    _kd = bg.kappa_dot_grid.copy()

    def rhs(tau, y):
        calH = np.interp(tau, _tau, _calH)
        a = np.interp(tau, _tau, _a)
        R = np.interp(tau, _tau, _R)
        kappa_dot = np.interp(tau, _tau, _kd)
        abs_kd = abs(kappa_dot)

        ia = 1.0 / a
        ia2 = ia * ia

        eta = y[RSA_IDX_ETA]
        delta_c = y[RSA_IDX_DELTA_C]
        delta_b = y[RSA_IDX_DELTA_B]
        theta_b = y[RSA_IDX_THETA_B]

        # ---- RSA metric computation ----
        # We need h' and eta' but they depend on delta_g, delta_n, theta_g, theta_n
        # which are algebraic functions of the metric itself.
        #
        # RSA (synchronous gauge, CLASS line 10252-10253):
        #   delta_g = (4/k^2) * (calH*h' - k^2*eta)
        #   theta_g = -(1/2) * h'
        #   delta_n = delta_g  (same RSA formula)
        #   theta_n = theta_g
        #
        # Substituting into the Einstein constraint for h':
        #   h' = (2/calH) * [k^2*eta + (3/2)*H0^2*(
        #           Og/a^2 * delta_g + On/a^2 * delta_n + Ob/a * delta_b + Oc/a * delta_c)]
        #
        # Let rho_rad_frac = (Og + On)/a^2 and substitute delta_g = (4/k^2)*(calH*h' - k^2*eta):
        #   h' = (2/calH) * [k^2*eta + (3/2)*H0^2*(
        #           rho_rad_frac * (4/k^2)*(calH*h' - k^2*eta) + Ob/a * delta_b + Oc/a * delta_c)]
        #
        # Solve for h':
        #   h' * [1 - (2/calH)*(3/2)*H0^2*rho_rad_frac*(4*calH/k^2)]
        #      = (2/calH)*[k^2*eta + (3/2)*H0^2*(
        #           -rho_rad_frac*4*eta + Ob/a*delta_b + Oc/a*delta_c)]
        #
        #   h' * [1 - 12*H0^2*rho_rad_frac/k^2]
        #      = (2/calH)*[k^2*eta - 6*H0^2*rho_rad_frac*eta
        #           + (3/2)*H0^2*(Ob/a*delta_b + Oc/a*delta_c)]

        rho_rad = (Og + On) * ia2
        coeff = 12.0 * H02 * rho_rad / k2

        rhs_h = ((2.0 / calH) * (k2 * eta - 6.0 * H02 * rho_rad * eta
                  + 1.5 * H02 * (Ob * ia * delta_b + Oc * ia * delta_c)))

        h_prime = rhs_h / (1.0 - coeff) if abs(1.0 - coeff) > 1e-20 else rhs_h

        # RSA algebraic values
        rsa_delta_g = (4.0 / k2) * (calH * h_prime - k2 * eta)
        rsa_theta_g = -0.5 * h_prime
        rsa_delta_n = rsa_delta_g
        rsa_theta_n = rsa_theta_g

        # Reionization correction (CLASS rsa_MD_with_reio, sync gauge):
        if abs_kd > 1e-10:
            rsa_delta_g += -(4.0 / k2) * (-abs_kd) * (theta_b + h_prime / 2.0)

        # eta' from 0i-constraint with RSA theta_g, theta_n
        src0i = ((4.0 / 3.0) * Og * ia2 * rsa_theta_g
                 + (4.0 / 3.0) * On * ia2 * rsa_theta_n
                 + Ob * ia * theta_b)
        eta_prime = 1.5 * H02 / k2 * src0i

        # --- CDM ---
        d_delta_c = -0.5 * h_prime

        # --- Baryons ---
        d_delta_b = -theta_b - 0.5 * h_prime
        cs2_b = baryon_cs2(a)
        # After decoupling, the Thomson coupling is negligible (kd -> 0)
        # but we keep the term for reionization
        theta_g_eff = rsa_theta_g
        d_theta_b = (-calH * theta_b
                     + cs2_b * k2 * delta_b
                     + abs_kd / R * (theta_g_eff - theta_b))

        # --- Assemble ---
        dydt = np.zeros(N_VAR_RSA)
        dydt[RSA_IDX_ETA] = eta_prime
        dydt[RSA_IDX_DELTA_C] = d_delta_c
        dydt[RSA_IDX_DELTA_B] = d_delta_b
        dydt[RSA_IDX_THETA_B] = d_theta_b

        return dydt

    return rhs, N_VAR_RSA


# ============================================================================
# RSA algebraic evaluation (for source function computation)
# ============================================================================

def rsa_evaluate(y_rsa, k, calH, a, abs_kd, theta_b_for_reio=None):
    """
    Compute RSA algebraic values for photon/neutrino quantities.

    Given the RSA state vector [eta, delta_c, delta_b, theta_b],
    returns dict with delta_g, theta_g, etc.
    """
    k2 = k * k
    H02 = _H0_MPC ** 2
    Og, On, Ob, Oc = _OMEGA_GAMMA, _OMEGA_NU, _OMEGA_B, _OMEGA_C

    eta = y_rsa[RSA_IDX_ETA]
    delta_c = y_rsa[RSA_IDX_DELTA_C]
    delta_b = y_rsa[RSA_IDX_DELTA_B]
    theta_b = y_rsa[RSA_IDX_THETA_B]

    ia = 1.0 / a
    ia2 = ia * ia

    rho_rad = (Og + On) * ia2
    coeff = 12.0 * H02 * rho_rad / k2

    rhs_h = ((2.0 / calH) * (k2 * eta - 6.0 * H02 * rho_rad * eta
              + 1.5 * H02 * (Ob * ia * delta_b + Oc * ia * delta_c)))
    h_prime = rhs_h / (1.0 - coeff) if abs(1.0 - coeff) > 1e-20 else rhs_h

    rsa_delta_g = (4.0 / k2) * (calH * h_prime - k2 * eta)
    rsa_theta_g = -0.5 * h_prime
    rsa_delta_n = rsa_delta_g
    rsa_theta_n = rsa_theta_g

    # Reionization correction
    if abs_kd > 1e-10:
        rsa_delta_g += -(4.0 / k2) * (-abs_kd) * (theta_b + h_prime / 2.0)

    src0i = ((4.0 / 3.0) * Og * ia2 * rsa_theta_g
             + (4.0 / 3.0) * On * ia2 * rsa_theta_n
             + Ob * ia * theta_b)
    eta_prime = 1.5 * H02 / k2 * src0i

    return {
        'h_prime': h_prime,
        'eta_prime': eta_prime,
        'delta_g': rsa_delta_g,
        'theta_g': rsa_theta_g,
        'delta_n': rsa_delta_n,
        'theta_n': rsa_theta_n,
        'sigma_g': 0.0,  # no shear in RSA
        'sigma_n': 0.0,
        'Pi': 0.0,  # no polarization in RSA
    }


# ============================================================================
# Initial conditions (TCA phase)
# ============================================================================

def adiabatic_ic_tca(k, tau_init, ln_max=L_NU_MAX):
    """
    Adiabatic growing-mode ICs for the TCA phase.

    Same physics as Ma & Bertschinger (1995) Sec. V, but with the
    reduced TCA state vector including theta_g.
    """
    R_nu = _OMEGA_NU / (_OMEGA_GAMMA + _OMEGA_NU)
    C = 1.0
    x = k * tau_init
    x2 = x * x

    nvar = n_var_tca(ln_max)
    fn_s = TCA_IDX_FN_START
    y0 = np.zeros(nvar, dtype=np.float64)

    # eta
    fac_eta = (5.0 + 4.0 * R_nu) / (12.0 * (15.0 + 4.0 * R_nu))
    y0[TCA_IDX_ETA] = C * (1.0 - fac_eta * x2)

    # CDM and baryons
    y0[TCA_IDX_DELTA_C] = -0.5 * C * x2
    y0[TCA_IDX_DELTA_B] = -0.5 * C * x2

    # Baryon velocity (tight coupled, ~ 0 at early times)
    y0[TCA_IDX_THETA_B] = 0.0

    # Photon monopole
    y0[TCA_IDX_DELTA_G] = -(2.0 / 3.0) * C * x2

    # Photon velocity (tight coupled to baryons, ~ 0 at early times)
    y0[TCA_IDX_THETA_G] = 0.0

    # Neutrino monopole
    y0[fn_s + 0] = -(2.0 / 3.0) * C * x2
    y0[fn_s + 1] = 0.0

    # Neutrino quadrupole
    if ln_max >= 2:
        fac_n2 = 32.0 * (5.0 + R_nu) / (45.0 * (15.0 + 4.0 * R_nu))
        y0[fn_s + 2] = fac_n2 * C * x2

    # Higher neutrino multipoles
    for ell in range(3, min(ln_max + 1, 6)):
        prod = 1.0
        for j in range(1, ell + 1):
            prod *= (2 * j + 1)
        y0[fn_s + ell] = C * x ** ell / prod

    return y0


# ============================================================================
# Seeding: TCA -> Full transition
# ============================================================================

def seed_full_from_tca(y_tca, k, bg, tau_switch,
                       lg_max=L_GAMMA_MAX, lp_max=L_POL_MAX,
                       ln_max=L_NU_MAX):
    """
    Create the full-phase state vector from the TCA state at the switch point.

    Seeding follows CLASS perturbations.c lines 4302-4398:
      - delta_g, theta_b, delta_b, delta_c, eta: direct copy
      - theta_g: from theta_b + slip
      - F_gamma,1: from theta_g via theta_g = (3/4)*k*F_g1
      - F_gamma,2: from TCA analytical shear (2*sigma_g)
      - F_gamma,3: second-order TCA: (6/7)*k/|kd|*sigma_g
      - E_0 = (5/2)*sigma_g,  E_1 = k/(2*|kd|)*sigma_g
      - E_2 = (1/2)*sigma_g,  E_3 = (3/14)*k/|kd|*sigma_g
      - All higher multipoles = 0
    """
    nvar = n_var_full(lg_max, lp_max, ln_max)
    y_full = np.zeros(nvar, dtype=np.float64)

    # Interpolate background at switch time
    abs_kd = float(np.interp(tau_switch, bg.tau_grid, np.abs(bg.kappa_dot_grid)))
    calH = float(np.interp(tau_switch, bg.tau_grid, bg.calH_grid))
    a = float(np.interp(tau_switch, bg.tau_grid, bg.a_grid))
    R = float(np.interp(tau_switch, bg.tau_grid, bg.R_grid))
    k2 = k * k
    H02 = _H0_MPC ** 2
    Og, On, Ob, Oc = _OMEGA_GAMMA, _OMEGA_NU, _OMEGA_B, _OMEGA_C
    ia = 1.0 / a
    ia2 = ia * ia

    # Direct copies
    y_full[FULL_IDX_ETA] = y_tca[TCA_IDX_ETA]
    y_full[FULL_IDX_DELTA_C] = y_tca[TCA_IDX_DELTA_C]
    y_full[FULL_IDX_DELTA_B] = y_tca[TCA_IDX_DELTA_B]
    y_full[FULL_IDX_THETA_B] = y_tca[TCA_IDX_THETA_B]

    delta_g = y_tca[TCA_IDX_DELTA_G]
    theta_g = y_tca[TCA_IDX_THETA_G]
    theta_b = y_tca[TCA_IDX_THETA_B]
    delta_b = y_tca[TCA_IDX_DELTA_B]
    delta_c = y_tca[TCA_IDX_DELTA_C]
    eta = y_tca[TCA_IDX_ETA]

    fn_s_tca = TCA_IDX_FN_START

    # Neutrino state: direct copy
    fn_s_full = full_idx_fn_start(lg_max, lp_max)
    n_nu = ln_max + 1
    y_full[fn_s_full: fn_s_full + n_nu] = y_tca[fn_s_tca: fn_s_tca + n_nu]

    delta_n = y_tca[fn_s_tca]
    theta_n = 0.75 * k * y_tca[fn_s_tca + 1]

    # Compute h', eta' at switch point
    src00 = (Og * ia2 * delta_g + On * ia2 * delta_n
             + Ob * ia * delta_b + Oc * ia * delta_c)
    h_prime = (2.0 / calH) * (k2 * eta + 1.5 * H02 * src00)

    src0i = ((4.0 / 3.0) * Og * ia2 * theta_g
             + (4.0 / 3.0) * On * ia2 * theta_n
             + Ob * ia * theta_b)
    eta_prime = 1.5 * H02 / k2 * src0i

    # TCA analytical shear (using theta_g from TCA state)
    sigma_g, Fg2 = compute_tca_shear(theta_g, h_prime, eta_prime, abs_kd)

    # Photon hierarchy
    y_full[full_idx_fg(0)] = delta_g           # F_g,0 = delta_g
    y_full[full_idx_fg(1)] = (4.0 / (3.0 * k)) * theta_g   # theta_g = (3/4)*k*F_g,1

    if lg_max >= 2:
        y_full[full_idx_fg(2)] = Fg2           # F_g,2 = 2*sigma_g

    if lg_max >= 3 and abs_kd > 0:
        y_full[full_idx_fg(3)] = (6.0 / 7.0) * (k / abs_kd) * sigma_g

    # Polarization seeding (CLASS lines 4368-4388)
    ep_s = full_idx_epol_start(lg_max)
    if lp_max >= 0:
        y_full[ep_s + 0] = (5.0 / 2.0) * sigma_g     # E_0
    if lp_max >= 1 and abs_kd > 0:
        y_full[ep_s + 1] = (k / abs_kd) * (3.0 / 6.0) * sigma_g  # E_1
    if lp_max >= 2:
        y_full[ep_s + 2] = 0.5 * sigma_g              # E_2
    if lp_max >= 3 and abs_kd > 0:
        y_full[ep_s + 3] = (k / abs_kd) * (3.0 / 14.0) * sigma_g  # E_3

    # All higher l = 0 (already from np.zeros)

    return y_full


# ============================================================================
# Seeding: Full -> RSA transition
# ============================================================================

def seed_rsa_from_full(y_full):
    """
    Create RSA state vector from the full-phase state.

    Simply copies {eta, delta_c, delta_b, theta_b}.
    """
    y_rsa = np.zeros(N_VAR_RSA, dtype=np.float64)
    y_rsa[RSA_IDX_ETA] = y_full[FULL_IDX_ETA]
    y_rsa[RSA_IDX_DELTA_C] = y_full[FULL_IDX_DELTA_C]
    y_rsa[RSA_IDX_DELTA_B] = y_full[FULL_IDX_DELTA_B]
    y_rsa[RSA_IDX_THETA_B] = y_full[FULL_IDX_THETA_B]
    return y_rsa


# ============================================================================
# Gauge transformation for each phase
# ============================================================================

def gauge_transform_full(y, k, calH, a, lg_max, lp_max, ln_max,
                         h_prime, eta_prime):
    """
    Gauge transform from sync to Newtonian for full-phase state vector.
    """
    k2 = k * k
    H02 = _H0_MPC ** 2
    Og, On = _OMEGA_GAMMA, _OMEGA_NU
    fn_s = full_idx_fn_start(lg_max, lp_max)
    ep_s = full_idx_epol_start(lg_max)

    eta = y[FULL_IDX_ETA]
    Fg = y[FULL_IDX_FG_START: FULL_IDX_FG_START + lg_max + 1]
    E = y[ep_s: ep_s + lp_max + 1]
    Fn = y[fn_s: fn_s + ln_max + 1]
    theta_b = y[FULL_IDX_THETA_B]

    alpha = (h_prime + 6.0 * eta_prime) / (2.0 * k2)
    if calH > 0:
        max_alpha = 5.0 * abs(eta) / calH
        alpha = np.clip(alpha, -max_alpha, max_alpha)

    Phi_N = eta - calH * alpha

    sigma_g = 0.5 * Fg[2] if lg_max >= 2 else 0.0
    sigma_n = 0.5 * Fn[2] if ln_max >= 2 else 0.0
    aniso = 6.0 * H02 / (a * a * k2) * (Og * sigma_g + On * sigma_n)
    Psi_N = Phi_N - aniso

    delta_g_N = Fg[0] - 4.0 * calH * alpha
    theta_b_N = theta_b + k2 * alpha

    # Pi / 4
    Fg2 = Fg[2] if lg_max >= 2 else 0.0
    E0 = E[0] if lp_max >= 0 else 0.0
    E2 = E[2] if lp_max >= 2 else 0.0
    Pi = Fg2 + E0 + E2

    return {
        'Phi_N': Phi_N,
        'Psi_N': Psi_N,
        'delta_g_N': delta_g_N,
        'theta_b_N': theta_b_N,
        'alpha': alpha,
        'Theta0_N': Fg[0] / 4.0 - eta + Phi_N,
        'vb_N': theta_b / k + (h_prime + 6.0 * eta_prime) / (2.0 * k),
        'Theta2': Fg[2] / 4.0 if lg_max >= 2 else 0.0,
        'Pi_over_4': Pi / 4.0,
    }


def gauge_transform_tca(y_tca, k, calH, a, bg, tau, ln_max=L_NU_MAX):
    """
    Gauge transform for TCA-phase state vector.
    """
    k2 = k * k
    H02 = _H0_MPC ** 2
    Og, On, Ob, Oc = _OMEGA_GAMMA, _OMEGA_NU, _OMEGA_B, _OMEGA_C

    eta = y_tca[TCA_IDX_ETA]
    delta_c = y_tca[TCA_IDX_DELTA_C]
    delta_b = y_tca[TCA_IDX_DELTA_B]
    theta_b = y_tca[TCA_IDX_THETA_B]
    delta_g = y_tca[TCA_IDX_DELTA_G]
    theta_g = y_tca[TCA_IDX_THETA_G]
    fn_s = TCA_IDX_FN_START
    Fn = y_tca[fn_s: fn_s + ln_max + 1]

    ia = 1.0 / a
    ia2 = ia * ia
    delta_n = Fn[0]
    theta_n = 0.75 * k * Fn[1]

    abs_kd = float(np.interp(tau, bg.tau_grid, np.abs(bg.kappa_dot_grid)))
    R = float(np.interp(tau, bg.tau_grid, bg.R_grid))

    # h', eta'
    src00 = (Og * ia2 * delta_g + On * ia2 * delta_n
             + Ob * ia * delta_b + Oc * ia * delta_c)
    h_prime = (2.0 / calH) * (k2 * eta + 1.5 * H02 * src00)

    src0i = ((4.0 / 3.0) * Og * ia2 * theta_g
             + (4.0 / 3.0) * On * ia2 * theta_n
             + Ob * ia * theta_b)
    eta_prime = 1.5 * H02 / k2 * src0i

    # TCA analytical quantities
    sigma_g, Fg2 = compute_tca_shear(theta_g, h_prime, eta_prime, abs_kd)

    alpha = (h_prime + 6.0 * eta_prime) / (2.0 * k2)
    if calH > 0:
        max_alpha = 5.0 * abs(eta) / calH
        alpha = np.clip(alpha, -max_alpha, max_alpha)

    Phi_N = eta - calH * alpha
    sigma_n = 0.5 * Fn[2] if ln_max >= 2 else 0.0
    aniso = 6.0 * H02 / (a * a * k2) * (Og * sigma_g + On * sigma_n)
    Psi_N = Phi_N - aniso

    # TCA polarization: P = (5/8)*sigma_g (CLASS line from source functions doc)
    # Pi/4 = (5/8)*sigma_g / ... wait, Pi = F_g2 + E_0 + E_2 = 5/2*Fg2 = 5*sigma_g
    # Pi_over_4 = 5*sigma_g/4
    Pi = 5.0 * sigma_g   # = F_g2 + E_0 + E_2 in TCA equilibrium

    return {
        'Phi_N': Phi_N,
        'Psi_N': Psi_N,
        'delta_g_N': delta_g - 4.0 * calH * alpha,
        'theta_b_N': theta_b + k2 * alpha,
        'alpha': alpha,
        'Theta0_N': delta_g / 4.0 - eta + Phi_N,
        'vb_N': theta_b / k + (h_prime + 6.0 * eta_prime) / (2.0 * k),
        'Theta2': Fg2 / 4.0,
        'Pi_over_4': Pi / 4.0,
        'h_prime': h_prime,
        'eta_prime': eta_prime,
    }


def gauge_transform_rsa(y_rsa, k, calH, a, abs_kd):
    """
    Gauge transform for RSA-phase state vector.
    """
    k2 = k * k
    H02 = _H0_MPC ** 2
    Og, On, Ob, Oc = _OMEGA_GAMMA, _OMEGA_NU, _OMEGA_B, _OMEGA_C

    rsa_vals = rsa_evaluate(y_rsa, k, calH, a, abs_kd)
    h_prime = rsa_vals['h_prime']
    eta_prime = rsa_vals['eta_prime']
    delta_g = rsa_vals['delta_g']

    eta = y_rsa[RSA_IDX_ETA]
    theta_b = y_rsa[RSA_IDX_THETA_B]

    alpha = (h_prime + 6.0 * eta_prime) / (2.0 * k2)
    if calH > 0:
        max_alpha = 5.0 * abs(eta) / calH
        alpha = np.clip(alpha, -max_alpha, max_alpha)

    Phi_N = eta - calH * alpha
    # sigma_g = sigma_n = 0 in RSA
    Psi_N = Phi_N

    return {
        'Phi_N': Phi_N,
        'Psi_N': Psi_N,
        'delta_g_N': delta_g - 4.0 * calH * alpha,
        'theta_b_N': theta_b + k2 * alpha,
        'alpha': alpha,
        'Theta0_N': delta_g / 4.0 - eta + Phi_N,
        'vb_N': theta_b / k + (h_prime + 6.0 * eta_prime) / (2.0 * k),
        'Theta2': 0.0,
        'Pi_over_4': 0.0,
        'h_prime': h_prime,
        'eta_prime': eta_prime,
    }


# ============================================================================
# Diagnose metric from full-phase state
# ============================================================================

def diagnose_metric_full(y, k, calH, a, lg_max, lp_max, ln_max):
    """Compute h' and eta' from Einstein constraints for full-phase state."""
    k2 = k * k
    H02 = _H0_MPC ** 2
    Og, On, Ob, Oc = _OMEGA_GAMMA, _OMEGA_NU, _OMEGA_B, _OMEGA_C
    fn_s = full_idx_fn_start(lg_max, lp_max)
    ia = 1.0 / a
    ia2 = ia * ia

    eta = y[FULL_IDX_ETA]
    delta_c = y[FULL_IDX_DELTA_C]
    delta_b = y[FULL_IDX_DELTA_B]
    theta_b = y[FULL_IDX_THETA_B]
    delta_g = y[full_idx_fg(0)]
    delta_n = y[fn_s]
    theta_g = 0.75 * k * y[full_idx_fg(1)]
    theta_n = 0.75 * k * y[fn_s + 1]

    src00 = (Og * ia2 * delta_g + On * ia2 * delta_n
             + Ob * ia * delta_b + Oc * ia * delta_c)
    h_prime = (2.0 / calH) * (k2 * eta + 1.5 * H02 * src00)

    src0i = ((4.0 / 3.0) * Og * ia2 * theta_g
             + (4.0 / 3.0) * On * ia2 * theta_n
             + Ob * ia * theta_b)
    eta_prime = 1.5 * H02 / k2 * src0i

    return h_prime, eta_prime


# ============================================================================
# Three-phase solver for a single k-mode
# ============================================================================

def solve_three_phase(k, bg, tau_end,
                      lg_max=L_GAMMA_MAX, lp_max=L_POL_MAX,
                      ln_max=L_NU_MAX,
                      method='BDF', rtol=1e-6, atol=1e-9,
                      verbose=False):
    """
    Three-phase integration for a single k-mode.

    Phase 1: tau_init -> tau_tca_end  (TCA, reduced system)
    Phase 2: tau_tca_end -> tau_rsa_start  (full Boltzmann hierarchy)
    Phase 3: tau_rsa_start -> tau_end  (RSA, algebraic radiation)

    Returns a CombinedSolution object with .sol(tau) dense output.
    """
    tau_init = float(bg.tau_grid[1])

    # --- Find switch times ---
    tau_tca_end = find_tca_switch_time(bg, k)
    tau_rsa_start = find_rsa_switch_time(bg, k)

    # Ensure proper ordering
    tau_tca_end = max(tau_tca_end, tau_init * 5)
    tau_tca_end = min(tau_tca_end, tau_end * 0.8)

    # If RSA starts before TCA ends, just do full hierarchy (no TCA/RSA)
    if tau_rsa_start <= tau_tca_end * 1.1:
        tau_rsa_start = tau_end + 1.0  # disable RSA

    # Ensure RSA doesn't start too early
    tau_rsa_start = max(tau_rsa_start, tau_tca_end * 1.5)

    if verbose:
        print(f"    k={k:.4e}: TCA=[{tau_init:.1f},{tau_tca_end:.1f}], "
              f"Full=[{tau_tca_end:.1f},{min(tau_rsa_start, tau_end):.1f}]"
              + (f", RSA=[{tau_rsa_start:.1f},{tau_end:.1f}]"
                 if tau_rsa_start < tau_end else ""))

    # ====== Phase 1: TCA ======
    rhs_tca, nvar_tca = make_tca_rhs(k, bg, ln_max)
    y0_tca = adiabatic_ic_tca(k, tau_init, ln_max)

    sol1 = solve_ivp(rhs_tca, [tau_init, tau_tca_end], y0_tca,
                     method=method, dense_output=True,
                     rtol=rtol, atol=atol)

    if not sol1.success:
        if verbose:
            print(f"    WARNING: TCA phase failed for k={k:.4e}: {sol1.message}")
        # Fallback: skip TCA, do full hierarchy from start
        return _fallback_full(k, bg, tau_init, tau_end, tau_rsa_start,
                              lg_max, lp_max, ln_max, method, rtol, atol)

    # ====== Seed Full from TCA ======
    y_tca_end = sol1.sol(tau_tca_end)
    y_full_start = seed_full_from_tca(y_tca_end, k, bg, tau_tca_end,
                                      lg_max, lp_max, ln_max)

    # ====== Phase 2: Full hierarchy ======
    rhs_full, nvar_full_ = make_full_rhs(k, bg, lg_max, lp_max, ln_max)

    tau_phase2_end = min(tau_rsa_start, tau_end)

    sol2 = solve_ivp(rhs_full, [tau_tca_end, tau_phase2_end], y_full_start,
                     method=method, dense_output=True,
                     rtol=rtol, atol=atol)

    if not sol2.success:
        if verbose:
            print(f"    WARNING: Full phase failed for k={k:.4e}: {sol2.message}")
        return _fallback_full(k, bg, tau_init, tau_end, tau_rsa_start,
                              lg_max, lp_max, ln_max, method, rtol, atol)

    # ====== Phase 3: RSA (if needed) ======
    sol3 = None
    if tau_rsa_start < tau_end:
        y_full_end = sol2.sol(tau_phase2_end)
        y_rsa_start = seed_rsa_from_full(y_full_end)

        rhs_rsa, nvar_rsa = make_rsa_rhs(k, bg)
        sol3 = solve_ivp(rhs_rsa, [tau_rsa_start, tau_end], y_rsa_start,
                         method=method, dense_output=True,
                         rtol=rtol, atol=atol)

        if not sol3.success and verbose:
            print(f"    WARNING: RSA phase failed for k={k:.4e}: {sol3.message}")
            sol3 = None  # will fall back to Phase 2 extrapolation

    # ====== Combine solutions ======
    return ThreePhaseSolution(sol1, sol2, sol3,
                              tau_tca_end, tau_rsa_start, tau_end,
                              k, bg, lg_max, lp_max, ln_max)


def _fallback_full(k, bg, tau_init, tau_end, tau_rsa_start,
                   lg_max, lp_max, ln_max, method, rtol, atol):
    """Fallback: full hierarchy from tau_init to tau_end (no TCA)."""
    from .perturbations_sync_tca import (
        make_sync_rhs_tca, adiabatic_ic_sync_tca,
    )
    rhs_fn, nvar = make_sync_rhs_tca(k, bg, lg_max, lp_max, ln_max,
                                     pol_approx='full')
    y0 = adiabatic_ic_sync_tca(k, tau_init, lg_max, lp_max, ln_max)
    sol = solve_ivp(rhs_fn, [tau_init, tau_end], y0,
                    method=method, dense_output=True,
                    rtol=rtol, atol=atol)
    # Wrap as a ThreePhaseSolution with only Phase 2
    return FallbackSolution(sol, k, bg, lg_max, lp_max, ln_max)


# ============================================================================
# Combined solution classes
# ============================================================================

class ThreePhaseSolution:
    """
    Combined solution from the three-phase integration.

    Provides .eval(tau) which returns a dict of Newtonian-gauge quantities
    for any tau in [tau_init, tau_end], dispatching to the correct phase.
    """
    def __init__(self, sol1, sol2, sol3,
                 tau_tca_end, tau_rsa_start, tau_end,
                 k, bg, lg_max, lp_max, ln_max):
        self.sol1 = sol1  # TCA
        self.sol2 = sol2  # Full
        self.sol3 = sol3  # RSA (may be None)
        self.tau_tca_end = tau_tca_end
        self.tau_rsa_start = tau_rsa_start
        self.tau_end = tau_end
        self.k = k
        self.bg = bg
        self.lg_max = lg_max
        self.lp_max = lp_max
        self.ln_max = ln_max
        self.success = sol2.success

    def eval(self, tau):
        """
        Evaluate at conformal time tau. Returns a dict with keys:
        Phi_N, Psi_N, Theta0_N, vb_N, Theta2, Pi_over_4
        """
        k = self.k
        a = float(np.interp(tau, self.bg.tau_grid, self.bg.a_grid))
        calH = float(np.interp(tau, self.bg.tau_grid, self.bg.calH_grid))
        abs_kd = float(np.interp(tau, self.bg.tau_grid,
                                  np.abs(self.bg.kappa_dot_grid)))

        if tau <= self.tau_tca_end:
            # Phase 1: TCA
            y_tca = self.sol1.sol(tau)
            return gauge_transform_tca(y_tca, k, calH, a, self.bg, tau,
                                       self.ln_max)

        elif tau <= self.tau_rsa_start or self.sol3 is None:
            # Phase 2: Full hierarchy
            # Clamp to the end of Phase 2 if we overshoot
            tau_eval = min(tau, self.sol2.t[-1]) if hasattr(self.sol2, 't') else tau
            y_full = self.sol2.sol(tau_eval)
            h_prime, eta_prime = diagnose_metric_full(
                y_full, k, calH, a, self.lg_max, self.lp_max, self.ln_max)
            gt = gauge_transform_full(
                y_full, k, calH, a, self.lg_max, self.lp_max, self.ln_max,
                h_prime, eta_prime)
            return gt

        else:
            # Phase 3: RSA
            y_rsa = self.sol3.sol(tau)
            return gauge_transform_rsa(y_rsa, k, calH, a, abs_kd)


class FallbackSolution:
    """Wrapper for single-phase full-hierarchy fallback."""

    def __init__(self, sol, k, bg, lg_max, lp_max, ln_max):
        self.sol = sol
        self.k = k
        self.bg = bg
        self.lg_max = lg_max
        self.lp_max = lp_max
        self.ln_max = ln_max
        self.success = sol.success

    def eval(self, tau):
        k = self.k
        a = float(np.interp(tau, self.bg.tau_grid, self.bg.a_grid))
        calH = float(np.interp(tau, self.bg.tau_grid, self.bg.calH_grid))

        y = self.sol.sol(tau)

        # Use the perturbations_sync_tca layout (same as full-phase layout)
        from .perturbations_sync_tca import (
            diagnose_metric_tca as _diag_met,
            gauge_transform_tca as _gt,
            IDX_ETA as _ETA, IDX_FG_START as _FG_S,
            IDX_THETA_B as _TB,
            idx_epol_start as _ep_s,
        )
        h_prime, eta_prime = _diag_met(y, k, calH, a,
                                       self.lg_max, self.lp_max, self.ln_max)
        gt = _gt(y, k, calH, a, self.lg_max, self.lp_max, self.ln_max,
                 h_prime, eta_prime, pol_approx='full')

        eta_val = y[_ETA]
        ep_s_val = _ep_s(self.lg_max)
        Fg2 = y[_FG_S + 2] if self.lg_max >= 2 else 0.0
        E0 = y[ep_s_val] if self.lp_max >= 0 else 0.0
        E2 = y[ep_s_val + 2] if self.lp_max >= 2 else 0.0

        return {
            'Phi_N': gt['Phi_N'],
            'Psi_N': gt['Psi_N'],
            'Theta0_N': y[_FG_S] / 4.0 - eta_val + gt['Phi_N'],
            'vb_N': y[_TB] / k + (h_prime + 6.0 * eta_prime) / (2.0 * k),
            'Theta2': Fg2 / 4.0,
            'Pi_over_4': (Fg2 + E0 + E2) / 4.0,
        }
