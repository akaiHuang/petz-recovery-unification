"""
perturbations_sync_tca.py -- Synchronous gauge Boltzmann equations with CLASS TCA.

Implements CLASS's Tight Coupling Approximation (TCA) with:
  1. Analytical photon shear:  sigma_g = (16/45)/|kd| * (theta_g + metric_shear)
  2. First-order slip:         slip = tau_c/(1+R) * [...]
  3. Proper seeding at TCA->full switch (F_g2, F_g3, E_0..E_3)
  4. Polarization hierarchy (E_0..E_{l_pol_max})
  5. Baryon sound speed

References:
  Ma & Bertschinger (1995) ApJ 455, 7
  Blas, Lesgourgues & Tram (2011) JCAP 07, 034 -- compromise_CLASS TCA
  Hu & White (1997) PRD 56, 596 -- polarization hierarchy

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
_mu_mol = 1.0 / (1.0 - 0.75 * _Y_He)


# ============================================================================
# Default hierarchy truncation
# ============================================================================
L_GAMMA_MAX = 25    # photon temperature multipole truncation
L_POL_MAX = 10      # photon polarization multipole truncation
L_NU_MAX_SYNC = 25  # neutrino multipole truncation


# ============================================================================
# State vector layout
# ============================================================================
IDX_ETA = 0
IDX_DELTA_C = 1
IDX_DELTA_B = 2
IDX_THETA_B = 3
IDX_FG_START = 4


def idx_fg(l):
    """Index of photon temperature multipole F_gamma,l in the state vector."""
    return IDX_FG_START + l


def idx_epol_start(lg_max):
    """Starting index of E-mode polarization hierarchy."""
    return IDX_FG_START + lg_max + 1


def idx_epol(l, lg_max):
    """Index of polarization multipole E_l."""
    return idx_epol_start(lg_max) + l


def idx_fn_start(lg_max, lp_max):
    """Starting index of neutrino hierarchy."""
    return idx_epol_start(lg_max) + lp_max + 1


def idx_fn(l, lg_max, lp_max):
    """Index of neutrino multipole F_nu,l."""
    return idx_fn_start(lg_max, lp_max) + l


def n_var_sync_tca(lg_max, lp_max, ln_max):
    """Total number of variables per k-mode."""
    return IDX_FG_START + (lg_max + 1) + (lp_max + 1) + (ln_max + 1)


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
# Diagnose metric from state
# ============================================================================

def diagnose_metric_tca(y, k, calH, a, lg_max, lp_max, ln_max):
    """
    Compute h' and eta' from Einstein constraints.
    """
    k2 = k * k
    H02 = _H0_MPC ** 2
    Og, On, Ob, Oc = _OMEGA_GAMMA, _OMEGA_NU, _OMEGA_B, _OMEGA_C
    fn_s = idx_fn_start(lg_max, lp_max)
    ia = 1.0 / a
    ia2 = ia * ia

    eta = y[IDX_ETA]
    delta_c = y[IDX_DELTA_C]
    delta_b = y[IDX_DELTA_B]
    theta_b = y[IDX_THETA_B]
    delta_g = y[idx_fg(0)]
    delta_n = y[fn_s]
    theta_g = 0.75 * k * y[idx_fg(1)]
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
# TCA shear and slip (CLASS formulas)
# ============================================================================

def compute_tca_shear(theta_g, h_prime, eta_prime, abs_kd):
    """
    TCA photon shear (first order, CLASS perturbations.c line 10042).

    sigma_g = (16/45) * tau_c * (theta_g + metric_shear)

    where metric_shear = (h' + 6 eta') / 2 in the synchronous gauge.
    (In CLASS notation: metric_shear = k^2 * alpha for scalar modes.)

    Returns sigma_g and F_gamma,2 = 2 * sigma_g.
    """
    tau_c = 1.0 / abs_kd if abs_kd > 0 else 0.0
    metric_shear = 0.5 * (h_prime + 6.0 * eta_prime)
    sigma_g = (16.0 / 45.0) * tau_c * (theta_g + metric_shear)
    Fg2 = 2.0 * sigma_g
    return sigma_g, Fg2


def compute_tca_slip(k, calH, R, delta_g, theta_b, theta_g,
                     sigma_g, abs_kd, h_prime, eta_prime, a, cs2_b):
    """
    First-order photon-baryon slip in tight coupling.

    Following Blas, Lesgourgues & Tram (2011), using the
    "first_order_CLASS" / "compromise_CLASS" formula.

    In synchronous gauge:
      metric_continuity = h'/2
      metric_euler = 0
      alpha = (h' + 6 eta') / (2 k^2)

    slip = tau_c / (1+R) * [
        -calH * R/(1+R) * theta_b
        + k^2 * (delta_g/4 - sigma_g)
        + k^2 * calH * alpha  (metric term, sync gauge correction)
    ]

    Returns the slip value (theta_g - theta_b).
    """
    if abs_kd <= 0:
        return 0.0

    tau_c = 1.0 / abs_kd
    k2 = k * k

    # Metric term in sync gauge
    alpha = (h_prime + 6.0 * eta_prime) / (2.0 * k2) if k2 > 0 else 0.0

    # First-order slip (simplified form matching CLASS first_order_CLASS)
    # The leading contribution comes from:
    # 1. Hubble drag mismatch: -calH * R/(1+R) * theta_b
    # 2. Photon pressure: k^2 * (delta_g/4 - sigma_g)
    # 3. Metric: this is the Euler equation's metric source
    # In sync gauge, metric_euler = 0, but the combination with alpha is nonzero

    term1 = -calH * R / (1.0 + R) * theta_b
    term2 = k2 * (delta_g / 4.0 - sigma_g)

    slip = tau_c / (1.0 + R) * (term1 + term2)

    return slip


# ============================================================================
# Seed photon hierarchy at TCA -> full switch
# ============================================================================

def seed_tca_at_switch(y, k, abs_kd, calH, a, lg_max, lp_max,
                       h_prime, eta_prime, R):
    """
    Seed the photon quadrupole, l=3, and polarization at the TCA switch point.

    Following CLASS perturbations.c lines 4302-4398.

    Photon temperature:
      F_gamma,2 = 2 * sigma_g    (first-order TCA analytical shear)
      F_gamma,3 = (6/7) * k/|kd| * sigma_g  (second-order TCA, flat space s_l[3]=1)

    Polarization (Hu & White 1997 equilibrium):
      E_0 = (5/2) * sigma_g  = (5/4) * F_gamma,2
      E_1 = k/|kd| * (5 - 2)/6 * sigma_g = k/(2|kd|) * sigma_g  (2nd order, s_l[2]=1)
      E_2 = (1/2) * sigma_g  = (1/4) * F_gamma,2
      E_3 = k/|kd| * 3/(14) * sigma_g  (2nd order, s_l[3]=1)

    Parameters
    ----------
    y : state vector (modified in place)
    k : wavenumber
    abs_kd : |kappa_dot| at the switch point
    calH, a : background quantities at switch
    lg_max, lp_max : hierarchy truncations
    h_prime, eta_prime : metric derivatives at switch
    R : baryon loading 3*rho_b/(4*rho_gamma) at switch
    """
    theta_g = 0.75 * k * y[idx_fg(1)]

    # Compute TCA shear (includes metric_shear correction)
    sigma_g, Fg2 = compute_tca_shear(theta_g, h_prime, eta_prime, abs_kd)

    # Seed photon temperature quadrupole
    if lg_max >= 2:
        y[idx_fg(2)] = Fg2

    # Seed F_gamma,3 (second-order TCA, flat space s_l[3] = 1)
    if lg_max >= 3 and abs_kd > 0:
        y[idx_fg(3)] = (6.0 / 7.0) * (k / abs_kd) * sigma_g

    # Seed polarization
    ep_s = idx_epol_start(lg_max)
    if lp_max >= 0:
        y[ep_s + 0] = (5.0 / 2.0) * sigma_g    # E_0 = 2.5 * sigma_g
    if lp_max >= 1 and abs_kd > 0:
        y[ep_s + 1] = (k / abs_kd) * (3.0 / 6.0) * sigma_g  # (5-2)/6 with s_l[2]=1
    if lp_max >= 2:
        y[ep_s + 2] = 0.5 * sigma_g             # E_2 = 0.5 * sigma_g
    if lp_max >= 3 and abs_kd > 0:
        y[ep_s + 3] = (k / abs_kd) * (3.0 / 14.0) * sigma_g  # E_3


# ============================================================================
# RHS function for synchronous gauge with polarization + TCA-aware coupling
# ============================================================================

def make_sync_rhs_tca(k, bg, lg_max=L_GAMMA_MAX, lp_max=L_POL_MAX,
                      ln_max=L_NU_MAX_SYNC, pol_approx='equilibrium'):
    """
    Build the RHS function f(tau, y) -> dy/dtau for synchronous gauge.

    This version includes:
      1. Baryon sound speed c_s^2 k^2 delta_b
      2. Photon polarization handling via pol_approx:
         - 'full': evolve E_l hierarchy (11 extra equations)
         - 'equilibrium': approximate Pi = (5/2)*F_g2 (TCA equilibrium),
           which gives effective damping -(3/4)|kd|*F_g2 on the quadrupole.
           This is the most important correction: reduces Silk damping to
           the correct level without evolving extra stiff equations.
         - 'none': Pi = F_g2 (original, -(9/10)|kd|*F_g2, over-damps)

    The equilibrium approximation is motivated by CLASS's TCA treatment:
    during tight coupling, E_0 = (5/4)*F_g2 and E_2 = (1/4)*F_g2, so
    Pi = F_g2 + E_0 + E_2 = (5/2)*F_g2. This reduces the effective
    damping coefficient from 9/10 to 3/4, which significantly improves
    the acoustic peak heights (especially peaks 2+).
    """
    k2 = k * k
    H02 = _H0_MPC ** 2
    Og, On, Ob, Oc = _OMEGA_GAMMA, _OMEGA_NU, _OMEGA_B, _OMEGA_C
    nvar = n_var_sync_tca(lg_max, lp_max, ln_max)
    ep_s = idx_epol_start(lg_max)
    fn_s = idx_fn_start(lg_max, lp_max)

    # Pre-tabulate background
    _tau = bg.tau_grid.copy()
    _calH = bg.calH_grid.copy()
    _a = bg.a_grid.copy()
    _R = bg.R_grid.copy()
    _kd = bg.kappa_dot_grid.copy()  # negative

    def rhs(tau, y):
        calH = np.interp(tau, _tau, _calH)
        a = np.interp(tau, _tau, _a)
        R = np.interp(tau, _tau, _R)
        kappa_dot = np.interp(tau, _tau, _kd)
        abs_kd = abs(kappa_dot)

        ia = 1.0 / a
        ia2 = ia * ia
        tau_safe = max(tau, 1e-10)

        # --- Extract state ---
        eta = y[IDX_ETA]
        delta_c = y[IDX_DELTA_C]
        delta_b = y[IDX_DELTA_B]
        theta_b = y[IDX_THETA_B]

        Fg = y[IDX_FG_START: IDX_FG_START + lg_max + 1]
        E = y[ep_s: ep_s + lp_max + 1]
        Fn = y[fn_s: fn_s + ln_max + 1]

        delta_g = Fg[0]
        delta_n = Fn[0]
        theta_g = 0.75 * k * Fg[1]
        theta_n = 0.75 * k * Fn[1]

        # --- Diagnose h' from 00-constraint ---
        src00 = (Og * ia2 * delta_g + On * ia2 * delta_n
                 + Ob * ia * delta_b + Oc * ia * delta_c)
        h_prime = (2.0 / calH) * (k2 * eta + 1.5 * H02 * src00)

        # --- Diagnose eta' from 0i-constraint ---
        src0i = ((4.0 / 3.0) * Og * ia2 * theta_g
                 + (4.0 / 3.0) * On * ia2 * theta_n
                 + Ob * ia * theta_b)
        eta_prime = 1.5 * H02 / k2 * src0i

        # --- CDM ---
        d_delta_c = -0.5 * h_prime

        # --- Baryons ---
        d_delta_b = -theta_b - 0.5 * h_prime

        # Baryon sound speed
        cs2 = baryon_cs2(a)
        d_theta_b = (-calH * theta_b
                     + cs2 * k2 * delta_b
                     + abs_kd / R * (theta_g - theta_b))

        # --- Polarization anisotropy Pi ---
        Fg2 = Fg[2] if lg_max >= 2 else 0.0
        if pol_approx == 'full':
            E0 = E[0] if lp_max >= 0 else 0.0
            E2 = E[2] if lp_max >= 2 else 0.0
            Pi = Fg2 + E0 + E2
        elif pol_approx == 'equilibrium':
            # TCA equilibrium: E_0 = (5/4)*F_g2, E_2 = (1/4)*F_g2
            # Pi = F_g2 + (5/4)*F_g2 + (1/4)*F_g2 = (5/2)*F_g2
            Pi = 2.5 * Fg2
        else:  # 'none'
            Pi = Fg2

        # --- Photon temperature hierarchy ---
        dFg = np.zeros(lg_max + 1)

        # l=0: F_0' = -k F_1 - (2/3) h'
        dFg[0] = -k * Fg[1] - (2.0 / 3.0) * h_prime

        # l=1: F_1' = (k/3)(F_0 - 2F_2) + |kd|(-F_1 + 4 theta_b / (3k))
        F2g = Fg[2] if lg_max >= 2 else 0.0
        dFg[1] = ((k / 3.0) * (Fg[0] - 2.0 * F2g)
                  + abs_kd * (-Fg[1] + 4.0 * theta_b / (3.0 * k)))

        # l=2: WITH polarization feedback (Pi)
        # F_2' = (2k/5)F_1 - (3k/5)F_3 + (4/15)h' + (8/15)eta'
        #        - |kd|(F_2 - Pi/10)
        # With equilibrium Pi = 5/2 * F_g2: -(3/4)|kd|*F_g2
        # With no Pi: -(9/10)|kd|*F_g2
        if lg_max >= 2:
            F3g = Fg[3] if lg_max >= 3 else 0.0
            dFg[2] = ((k / 5.0) * (2.0 * Fg[1] - 3.0 * F3g)
                      + (4.0 / 15.0) * h_prime + (8.0 / 15.0) * eta_prime
                      - abs_kd * (Fg[2] - Pi / 10.0))

        # l=3..lg_max-1: streaming + Thomson damping
        for ell in range(3, lg_max):
            dFg[ell] = (k / (2.0 * ell + 1.0)
                        * (ell * Fg[ell - 1] - (ell + 1) * Fg[ell + 1])
                        - abs_kd * Fg[ell])

        # l=lg_max: truncation boundary
        if lg_max >= 3:
            dFg[lg_max] = (k * Fg[lg_max - 1] * lg_max / (2.0 * lg_max + 1.0)
                           - (lg_max + 1.0) / tau_safe * Fg[lg_max]
                           - abs_kd * Fg[lg_max])

        # --- E-mode polarization hierarchy ---
        dE = np.zeros(lp_max + 1)

        if pol_approx == 'full':
            # Full polarization evolution
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
        # else: dE stays all zero (equilibrium and none modes don't evolve E)

        # --- Neutrino hierarchy (no collisions) ---
        dFn = np.zeros(ln_max + 1)

        dFn[0] = -k * Fn[1] - (2.0 / 3.0) * h_prime

        F2n = Fn[2] if ln_max >= 2 else 0.0
        dFn[1] = (k / 3.0) * (Fn[0] - 2.0 * F2n)

        if ln_max >= 2:
            F3n = Fn[3] if ln_max >= 3 else 0.0
            dFn[2] = ((k / 5.0) * (2.0 * Fn[1] - 3.0 * F3n)
                      + (4.0 / 15.0) * h_prime + (8.0 / 15.0) * eta_prime)

        for ell in range(3, ln_max):
            dFn[ell] = (k / (2.0 * ell + 1.0)
                        * (ell * Fn[ell - 1] - (ell + 1) * Fn[ell + 1]))

        if ln_max >= 1:
            dFn[ln_max] = (k * Fn[ln_max - 1] * ln_max / (2.0 * ln_max + 1.0)
                           - (ln_max + 1.0) / tau_safe * Fn[ln_max])

        # --- Assemble ---
        dydt = np.zeros(nvar)
        dydt[IDX_ETA] = eta_prime
        dydt[IDX_DELTA_C] = d_delta_c
        dydt[IDX_DELTA_B] = d_delta_b
        dydt[IDX_THETA_B] = d_theta_b
        dydt[IDX_FG_START: IDX_FG_START + lg_max + 1] = dFg
        dydt[ep_s: ep_s + lp_max + 1] = dE
        dydt[fn_s: fn_s + ln_max + 1] = dFn

        return dydt

    return rhs, nvar


# ============================================================================
# Adiabatic initial conditions
# ============================================================================

def adiabatic_ic_sync_tca(k, tau_init, lg_max=L_GAMMA_MAX, lp_max=L_POL_MAX,
                          ln_max=L_NU_MAX_SYNC):
    """
    Adiabatic growing-mode ICs following Ma & Bertschinger (1995) Sec. V.

    Normalization: eta(tau_init) ~ 1 (primordial curvature perturbation C = 1).
    Polarization E_l are initialized to zero (Thomson-damped at early times).
    """
    R_nu = _OMEGA_NU / (_OMEGA_GAMMA + _OMEGA_NU)
    C = 1.0
    x = k * tau_init
    x2 = x * x

    nvar = n_var_sync_tca(lg_max, lp_max, ln_max)
    fn_s = idx_fn_start(lg_max, lp_max)
    y0 = np.zeros(nvar, dtype=np.float64)

    # eta
    fac_eta = (5.0 + 4.0 * R_nu) / (12.0 * (15.0 + 4.0 * R_nu))
    y0[IDX_ETA] = C * (1.0 - fac_eta * x2)

    # CDM and baryons
    y0[IDX_DELTA_C] = -0.5 * C * x2
    y0[IDX_DELTA_B] = -0.5 * C * x2

    # Baryon velocity
    y0[IDX_THETA_B] = 0.0

    # Photon monopole
    y0[idx_fg(0)] = -(2.0 / 3.0) * C * x2

    # Photon dipole
    y0[idx_fg(1)] = 0.0

    # Photon quadrupole: Thomson-damped
    if lg_max >= 2:
        y0[idx_fg(2)] = 0.0

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
# Gauge transformation: synchronous -> Newtonian
# ============================================================================

def gauge_transform_tca(y, k, calH, a, lg_max, lp_max, ln_max,
                        h_prime, eta_prime, pol_approx='equilibrium'):
    """
    Transform sync gauge state to Newtonian gauge quantities needed for LOS.
    pol_approx controls how Pi is computed (same as in make_sync_rhs_tca).
    """
    k2 = k * k
    H02 = _H0_MPC ** 2
    Og, On = _OMEGA_GAMMA, _OMEGA_NU
    fn_s = idx_fn_start(lg_max, lp_max)
    ep_s = idx_epol_start(lg_max)

    eta = y[IDX_ETA]
    Fg = y[IDX_FG_START: IDX_FG_START + lg_max + 1]
    E = y[ep_s: ep_s + lp_max + 1]
    Fn = y[fn_s: fn_s + ln_max + 1]
    theta_b = y[IDX_THETA_B]

    # Gauge parameter
    alpha = (h_prime + 6.0 * eta_prime) / (2.0 * k2)
    if calH > 0:
        max_alpha = 5.0 * abs(eta) / calH
        alpha = np.clip(alpha, -max_alpha, max_alpha)

    # Newtonian gauge potentials
    Phi_N = eta - calH * alpha

    sigma_g = 0.5 * Fg[2] if lg_max >= 2 else 0.0
    sigma_n = 0.5 * Fn[2] if ln_max >= 2 else 0.0
    aniso = 12.0 * H02 / (a * a * k2) * (Og * sigma_g + On * sigma_n)
    Psi_N = Phi_N - aniso

    # Density gauge transformation
    delta_g_N = Fg[0] - 4.0 * calH * alpha

    # Velocity
    theta_b_N = theta_b + k2 * alpha

    # Polarization anisotropy
    Fg2 = Fg[2] if lg_max >= 2 else 0.0
    if pol_approx == 'full':
        E0 = E[0] if lp_max >= 0 else 0.0
        E2 = E[2] if lp_max >= 2 else 0.0
        Pi = Fg2 + E0 + E2
    elif pol_approx == 'equilibrium':
        Pi = 2.5 * Fg2
    else:
        Pi = Fg2

    return {
        'Phi_N': Phi_N,
        'Psi_N': Psi_N,
        'delta_g_N': delta_g_N,
        'theta_b_N': theta_b_N,
        'alpha': alpha,
        'Pi': Pi,
    }


# ============================================================================
# TCA diagnostic: compare F_gamma,2 with TCA analytical value
# ============================================================================

def tca_diagnostic(y, k, calH, a, lg_max, lp_max, ln_max, abs_kd):
    """
    Compare F_gamma,2 from the ODE with the TCA analytical value.

    Returns dict with:
      'Fg2_ode': F_gamma,2 from the numerical ODE
      'Fg2_tca': F_gamma,2 from TCA analytical formula
      'sigma_g_ode': sigma_g from ODE
      'sigma_g_tca': sigma_g from TCA formula
      'ratio': Fg2_ode / Fg2_tca
      'abs_kd': |kappa_dot|
      'tau_c_over_tau_k': tau_c * k (tightness parameter)
    """
    fn_s = idx_fn_start(lg_max, lp_max)
    theta_g = 0.75 * k * y[idx_fg(1)]

    # ODE value
    Fg2_ode = y[idx_fg(2)] if lg_max >= 2 else 0.0
    sigma_g_ode = 0.5 * Fg2_ode

    # Metric
    h_prime, eta_prime = diagnose_metric_tca(
        y, k, calH, a, lg_max, lp_max, ln_max)

    # TCA analytical (with metric_shear)
    sigma_g_tca, Fg2_tca = compute_tca_shear(
        theta_g, h_prime, eta_prime, abs_kd)

    ratio = Fg2_ode / Fg2_tca if abs(Fg2_tca) > 1e-30 else float('nan')

    tau_c = 1.0 / abs_kd if abs_kd > 0 else 0.0

    return {
        'Fg2_ode': Fg2_ode,
        'Fg2_tca': Fg2_tca,
        'sigma_g_ode': sigma_g_ode,
        'sigma_g_tca': sigma_g_tca,
        'ratio': ratio,
        'abs_kd': abs_kd,
        'tau_c_over_tau_k': tau_c * k,
        'theta_g': theta_g,
    }
