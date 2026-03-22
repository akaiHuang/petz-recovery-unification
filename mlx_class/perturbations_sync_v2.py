"""
perturbations_sync_v2.py -- Synchronous gauge Boltzmann equations WITH polarization.

Extends perturbations_sync.py with three missing physics pieces:
  1. Photon polarization hierarchy (E_0, E_1, ..., E_{l_pol_max})
     - Thomson scattering couples temperature quadrupole F_γ,2 to E_0, E_2
     - Polarization feedback: (1/10)|κ̇|(E_0 + E_2) in F_γ,2 equation
  2. TCA shear seeding at the TCA->full hierarchy switch point
     - F_γ,2(switch) = (32/45) θ_γ / |κ̇|  (not zero)
     - E_0(switch) = (5/4) F_γ,2           (equilibrium)
     - E_2(switch) = F_γ,2 / 4             (equilibrium)
  3. Baryon sound speed c_s² k² δ_b term (small but included for completeness)

References:
  Ma & Bertschinger (1995) ApJ 455, 7 [astro-ph/9506072]
  Hu & White (1997) PRD 56, 596 [astro-ph/9702170] -- polarization hierarchy
  Kosowsky (1996) Ann. Phys. 246, 49 -- Thomson scattering

STATE VECTOR per k-mode:
  [eta, delta_c, delta_b, theta_b,
   F_g0, F_g1, ..., F_g{lg_max},
   E_0, E_1, ..., E_{l_pol_max},
   F_n0, F_n1, ..., F_n{ln_max}]

  where:
    eta         = conformal metric perturbation (evolved)
    delta_c     = CDM density contrast
    delta_b     = baryon density contrast
    theta_b     = baryon velocity divergence
    F_g,l       = photon temperature hierarchy multipoles
    E_l         = photon E-mode polarization hierarchy multipoles
    F_n,l       = neutrino Boltzmann hierarchy multipoles

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
# Physical constants for baryon sound speed
# ============================================================================
_k_B_SI = 1.380649e-23       # J/K
_mu_mol = 1.0 / (1.0 - 0.75 * _Y_He)  # mean molecular weight ~ 1/(1-0.75*Y_He) ~ 1.22
# Actually for neutral hydrogen-helium gas: mu = 4/(4 - 3*Y_He) ~ 1.22
# But a better approximation: mu = m_mean / m_H where m_mean accounts for He
# Standard: mu = 1 / (X_H + Y_He/4) where X_H = 1-Y_He
# So mu = 1 / ((1-Y_He) + Y_He/4) = 1 / (1 - 3*Y_He/4) ~ 1.22
_mu_mol = 1.0 / (1.0 - 0.75 * _Y_He)


# ============================================================================
# Default hierarchy truncation
# ============================================================================
L_GAMMA_MAX = 25    # photon temperature multipole truncation
L_POL_MAX = 10      # photon polarization multipole truncation (more damped)
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


def n_var_sync_v2(lg_max, lp_max, ln_max):
    """Total number of variables per k-mode."""
    return IDX_FG_START + (lg_max + 1) + (lp_max + 1) + (ln_max + 1)


# ============================================================================
# Baryon sound speed
# ============================================================================

def baryon_cs2(a):
    """
    Baryon sound speed squared c_s^2 = k_B T_b / (mu m_H c^2).

    During tight coupling T_b ≈ T_CMB / a (Compton heating keeps baryons
    at the photon temperature). After decoupling T_b ~ 1/a^2 (adiabatic cooling).

    For simplicity, we use T_b = T_CMB / a (valid during recombination where
    this term matters most).

    Returns c_s^2 in natural units (dimensionless, c=1).
    """
    T_b = _T_CMB / a  # baryon temperature in Kelvin
    cs2 = _k_B_SI * T_b / (_mu_mol * _m_H_SI * _c_SI**2)
    return cs2


# ============================================================================
# RHS function for synchronous gauge with polarization
# ============================================================================

def make_sync_rhs_v2(k, bg, lg_max=L_GAMMA_MAX, lp_max=L_POL_MAX,
                     ln_max=L_NU_MAX_SYNC):
    """
    Build the RHS function f(tau, y) -> dy/dtau for synchronous gauge
    WITH photon polarization hierarchy and baryon sound speed.

    Parameters
    ----------
    k : float, wavenumber in Mpc^{-1}
    bg : Background object (solved)
    lg_max : int, photon temperature hierarchy truncation
    lp_max : int, photon polarization hierarchy truncation
    ln_max : int, neutrino hierarchy truncation

    Returns
    -------
    rhs : callable f(tau, y) -> dy/dtau
    nvar : int, length of state vector
    """
    k2 = k * k
    H02 = _H0_MPC ** 2
    Og, On, Ob, Oc = _OMEGA_GAMMA, _OMEGA_NU, _OMEGA_B, _OMEGA_C
    nvar = n_var_sync_v2(lg_max, lp_max, ln_max)
    ep_s = idx_epol_start(lg_max)
    fn_s = idx_fn_start(lg_max, lp_max)

    # Pre-tabulate background for fast np.interp
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
        theta_g = 0.75 * k * Fg[1]   # theta_gamma = (3/4) k F_1
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

        # --- CDM (theta_c = 0, gauge choice) ---
        d_delta_c = -0.5 * h_prime

        # --- Baryons ---
        d_delta_b = -theta_b - 0.5 * h_prime

        # Baryon sound speed: c_s^2 k^2 delta_b
        cs2 = baryon_cs2(a)
        d_theta_b = (-calH * theta_b
                     + cs2 * k2 * delta_b
                     + abs_kd / R * (theta_g - theta_b))

        # --- Polarization anisotropy Pi ---
        # Pi = F_gamma,2 + E_0 + E_2  (the source for Thomson scattering l=2)
        E0 = E[0] if lp_max >= 0 else 0.0
        E2 = E[2] if lp_max >= 2 else 0.0
        Fg2 = Fg[2] if lg_max >= 2 else 0.0
        Pi = Fg2 + E0 + E2

        # --- Photon temperature hierarchy ---
        dFg = np.zeros(lg_max + 1)

        # l=0: F_0' = -k F_1 - (2/3) h'
        dFg[0] = -k * Fg[1] - (2.0 / 3.0) * h_prime

        # l=1: F_1' = (k/3)(F_0 - 2F_2) + |kd|(-F_1 + 4 theta_b / (3k))
        F2g = Fg[2] if lg_max >= 2 else 0.0
        dFg[1] = ((k / 3.0) * (Fg[0] - 2.0 * F2g)
                  + abs_kd * (-Fg[1] + 4.0 * theta_b / (3.0 * k)))

        # l=2: WITH polarization feedback
        # F_2' = (2k/5)F_1 - (3k/5)F_3 + (4/15)h' + (8/15)eta'
        #        - |kd|(F_2 - Pi/10)
        # Expanding: -|kd|F_2 + |kd|Pi/10 = -(9/10)|kd|F_2 + (1/10)|kd|(E_0+E_2)
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
        # Hu & White (1997), Kosowsky (1996)
        #
        # E_0' = -k E_1 - |kd| E_0 + (1/2)|kd| Pi
        # E_1' = (k/3)(E_0 - 2 E_2) - |kd| E_1
        # E_2' = (k/5)(2 E_1 - 3 E_3) - |kd| E_2 + (1/10)|kd| Pi
        # E_l' = k/(2l+1)(l E_{l-1} - (l+1) E_{l+1}) - |kd| E_l   for l >= 3
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

        # l=3..lp_max-1: streaming + Thomson damping
        for ell in range(3, lp_max):
            dE[ell] = (k / (2.0 * ell + 1.0)
                       * (ell * E[ell - 1] - (ell + 1) * E[ell + 1])
                       - abs_kd * E[ell])

        # l=lp_max: truncation boundary
        if lp_max >= 3:
            dE[lp_max] = (k * E[lp_max - 1] * lp_max / (2.0 * lp_max + 1.0)
                          - (lp_max + 1.0) / tau_safe * E[lp_max]
                          - abs_kd * E[lp_max])

        # --- Neutrino hierarchy (no collisions) ---
        dFn = np.zeros(ln_max + 1)

        # l=0
        dFn[0] = -k * Fn[1] - (2.0 / 3.0) * h_prime

        # l=1
        F2n = Fn[2] if ln_max >= 2 else 0.0
        dFn[1] = (k / 3.0) * (Fn[0] - 2.0 * F2n)

        # l=2
        if ln_max >= 2:
            F3n = Fn[3] if ln_max >= 3 else 0.0
            dFn[2] = ((k / 5.0) * (2.0 * Fn[1] - 3.0 * F3n)
                      + (4.0 / 15.0) * h_prime + (8.0 / 15.0) * eta_prime)

        # l=3..ln_max-1: free streaming
        for ell in range(3, ln_max):
            dFn[ell] = (k / (2.0 * ell + 1.0)
                        * (ell * Fn[ell - 1] - (ell + 1) * Fn[ell + 1]))

        # l=ln_max: truncation
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
# Adiabatic initial conditions (synchronous gauge, CDM rest frame)
# ============================================================================

def adiabatic_ic_sync_v2(k, tau_init, lg_max=L_GAMMA_MAX, lp_max=L_POL_MAX,
                         ln_max=L_NU_MAX_SYNC):
    """
    Adiabatic growing-mode ICs following Ma & Bertschinger (1995) Sec. V.

    Normalization: eta(tau_init) ~ 1 (primordial curvature perturbation C = 1).

    Polarization E_l are initialized to zero (Thomson-damped at early times).
    They are seeded at the TCA switch point if TCA seeding is used.
    """
    R_nu = _OMEGA_NU / (_OMEGA_GAMMA + _OMEGA_NU)
    C = 1.0
    x = k * tau_init
    x2 = x * x

    nvar = n_var_sync_v2(lg_max, lp_max, ln_max)
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

    # Photon monopole: delta_gamma = -(2/3) C x^2
    y0[idx_fg(0)] = -(2.0 / 3.0) * C * x2

    # Photon dipole: O(x^3), negligible
    y0[idx_fg(1)] = 0.0

    # Photon quadrupole: Thomson-damped, set to 0
    if lg_max >= 2:
        y0[idx_fg(2)] = 0.0

    # Polarization: all zero initially (Thomson-damped)
    # E_l are at indices ep_s to ep_s + lp_max

    # Neutrino monopole
    y0[fn_s + 0] = -(2.0 / 3.0) * C * x2

    # Neutrino dipole
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
# TCA shear seeding
# ============================================================================

def seed_tca_shear(y, k, abs_kd, lg_max, lp_max):
    """
    Seed the photon quadrupole and polarization at the TCA switch point.

    In tight coupling, the photon quadrupole has a small but important
    equilibrium value from the steady-state of the l=2 equation:
        sigma_gamma = (16/45) theta_gamma / |kappa_dot|
        F_gamma,2   = 2 sigma_gamma = (32/45) theta_gamma / |kappa_dot|

    The polarization equilibrium values (Hu & White 1997):
        E_0 = (5/4) F_gamma,2
        E_2 = (1/4) F_gamma,2

    Parameters
    ----------
    y : state vector (modified in place)
    k : wavenumber
    abs_kd : |kappa_dot| at the switch point
    lg_max : photon temperature hierarchy truncation
    lp_max : polarization hierarchy truncation
    """
    theta_g = 0.75 * k * y[idx_fg(1)]

    # TCA equilibrium shear
    if abs_kd > 0:
        Fg2_tca = (32.0 / 45.0) * theta_g / abs_kd
    else:
        Fg2_tca = 0.0

    # Seed photon quadrupole
    if lg_max >= 2:
        y[idx_fg(2)] = Fg2_tca

    # Seed polarization
    ep_s = idx_epol_start(lg_max)
    if lp_max >= 0:
        y[ep_s + 0] = (5.0 / 4.0) * Fg2_tca  # E_0
    if lp_max >= 2:
        y[ep_s + 2] = (1.0 / 4.0) * Fg2_tca  # E_2


# ============================================================================
# Gauge transformation: synchronous -> Newtonian
# ============================================================================

def gauge_transform_v2(y, k, calH, a, lg_max, lp_max, ln_max,
                       h_prime, eta_prime):
    """
    Transform sync gauge state to Newtonian gauge quantities needed for LOS.

    Same as perturbations_sync.gauge_transform but aware of the polarization
    hierarchy layout.

    Returns
    -------
    dict with keys: 'Phi_N', 'Psi_N', 'delta_g_N', 'theta_b_N', 'Pi'
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

    # Polarization anisotropy Pi = F_gamma,2 + E_0 + E_2
    E0 = E[0] if lp_max >= 0 else 0.0
    E2 = E[2] if lp_max >= 2 else 0.0
    Fg2 = Fg[2] if lg_max >= 2 else 0.0
    Pi = Fg2 + E0 + E2

    return {
        'Phi_N': Phi_N,
        'Psi_N': Psi_N,
        'delta_g_N': delta_g_N,
        'theta_b_N': theta_b_N,
        'alpha': alpha,
        'Pi': Pi,
    }


def diagnose_metric_v2(y, k, calH, a, lg_max, lp_max, ln_max):
    """
    Compute h' and eta' from Einstein constraints for a given state.

    Same physics as perturbations_sync.diagnose_metric but aware of the
    polarization hierarchy layout.
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
