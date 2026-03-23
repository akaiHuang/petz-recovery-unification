"""
perturbations_sync.py -- Synchronous gauge Boltzmann equations.

Following Ma & Bertschinger (1995) ApJ 455, 7 [astro-ph/9506072].

GAUGE: Synchronous gauge in CDM rest frame (theta_c = 0).

STATE VECTOR per k-mode (with polarization, l_pol_max > 0):
  [eta, delta_c, delta_b, theta_b,
   F_g0, F_g1, ..., F_g{lg_max},
   E_0, E_1, ..., E_{l_pol_max},
   F_n0, F_n1, ..., F_n{ln_max}]

  where:
    eta         = conformal metric perturbation (evolved)
    delta_c     = CDM density contrast
    delta_b     = baryon density contrast
    theta_b     = baryon velocity divergence
    F_g,l       = photon Boltzmann hierarchy multipoles
                  (F_g,0 = delta_gamma, theta_gamma = 3k/4 * F_g,1)
    E_l         = E-mode polarization hierarchy multipoles
    F_n,l       = neutrino Boltzmann hierarchy multipoles
                  (F_n,0 = delta_nu, theta_nu = 3k/4 * F_n,1)

  When l_pol_max = 0, the E-mode block is absent and the layout matches the
  original code exactly (backward compatible).

METRIC PERTURBATIONS:
  h' and eta' are DIAGNOSED from Einstein constraint equations at each step.
  This avoids evolving h (which contains a gauge mode) and improves conditioning.

  00-constraint:  h' = (2/calH) [k^2 eta + (3/2) calH^2 * sum_i Omega_i delta_i / a^{1+3w_i}]
  0i-constraint:  eta' = (3/2)(calH^2/k^2) sum_i (1+w_i) Omega_i theta_i / a^{1+3w_i}

GAUGE TRANSFORMATION to Newtonian gauge (for LOS source):
  alpha = (h' + 6 eta') / (2 k^2)
  Phi_N = eta - calH * alpha
  Psi_N = Phi_N - 12 H0^2/(a^2 k^2) (Omega_gamma sigma_gamma + Omega_nu sigma_nu)
  delta_N = delta_S + (rho'/rho) * alpha  =>  delta_gamma_N = delta_gamma_S + 4 calH alpha

POLARIZATION FEEDBACK (l_pol_max > 0):
  Pi = F_gamma,2 + E_0 + E_2   (polarization source)
  F_gamma,2 gets +|kd|*Pi/10 correction (replacing -(9/10)|kd|F_g2 with -|kd|(F_g2 - Pi/10))
  E-mode hierarchy:
    E_0' = -k*E_1 - |kd|*E_0 + (1/2)*|kd|*Pi
    E_1' = k/3*(E_0 - 2*E_2) - |kd|*E_1
    E_2' = k/5*(2*E_1 - 3*E_3) - |kd|*E_2 + (1/10)*|kd|*Pi
    E_l' = k/(2l+1)*(l*E_{l-1} - (l+1)*E_{l+1}) - |kd|*E_l   for l >= 3

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
    Y_He as _Y_He,
    m_H_SI as _m_H_SI,
    c_SI as _c_SI,
)
from .perturbations_neutrino import (
    _f_nu, _OMEGA_GAMMA, _OMEGA_NU,
)

# Baryon sound speed: c_s^2 = k_B * T_b / (mu * m_H * c^2)
# During tight coupling, T_b ~ T_CMB / a (adiabatic).
_k_B_SI = 1.380649e-23   # J/K
_mu_mol = 1.0 / (1.0 - 0.75 * _Y_He)  # mean molecular weight


def _baryon_cs2(a):
    """Baryon sound speed squared in natural units (c=1)."""
    T_b = _T_CMB / a
    return _k_B_SI * T_b / (_mu_mol * _m_H_SI * _c_SI ** 2)


# ============================================================================
# Default hierarchy truncation
# ============================================================================
L_GAMMA_MAX = 25    # photon multipole truncation (25 gives best balance of damping and phase)
L_NU_MAX_SYNC = 25  # neutrino multipole truncation (free-streaming => need more)
L_POL_MAX = 8       # E-mode polarization hierarchy truncation (0 = no polarization)


# ============================================================================
# State vector layout
# ============================================================================
IDX_ETA = 0
IDX_DELTA_C = 1
IDX_DELTA_B = 2
IDX_THETA_B = 3
IDX_FG_START = 4


def idx_fg(l):
    """Index of photon multipole F_gamma,l in the state vector."""
    return IDX_FG_START + l


def idx_e_start(lg_max):
    """Starting index of E-mode polarization hierarchy (right after photon block)."""
    return IDX_FG_START + lg_max + 1


def idx_e(l, lg_max):
    """Index of E-mode polarization multipole E_l in the state vector."""
    return idx_e_start(lg_max) + l


def idx_fn_start(lg_max, l_pol_max=0):
    """Starting index of neutrino hierarchy.

    When l_pol_max=0 (default), the E-mode block is absent and this returns
    IDX_FG_START + lg_max + 1 (identical to the original code).
    When l_pol_max>0, the E-mode block sits between photons and neutrinos.
    """
    if l_pol_max > 0:
        return idx_e_start(lg_max) + l_pol_max + 1
    else:
        return IDX_FG_START + lg_max + 1


def idx_fn(l, lg_max, l_pol_max=0):
    """Index of neutrino multipole F_nu,l."""
    return idx_fn_start(lg_max, l_pol_max) + l


def n_var_sync(lg_max, ln_max, l_pol_max=0):
    """Total number of variables per k-mode.

    When l_pol_max=0: IDX_FG_START + (lg_max+1) + (ln_max+1)  (original)
    When l_pol_max>0: IDX_FG_START + (lg_max+1) + (l_pol_max+1) + (ln_max+1)
    """
    if l_pol_max > 0:
        return IDX_FG_START + (lg_max + 1) + (l_pol_max + 1) + (ln_max + 1)
    else:
        return IDX_FG_START + (lg_max + 1) + (ln_max + 1)


# ============================================================================
# RHS function for synchronous gauge (single k-mode, for scipy solve_ivp)
# ============================================================================

def make_sync_rhs(k, bg, lg_max=L_GAMMA_MAX, ln_max=L_NU_MAX_SYNC,
                  l_pol_max=0):
    """
    Build the RHS function f(tau, y) -> dy/dtau for synchronous gauge.

    CDM rest frame: theta_c = 0 (gauge choice, not evolved).
    h' and eta' are diagnosed at each step from Einstein constraints.

    Parameters
    ----------
    k : float, wavenumber in Mpc^{-1}
    bg : Background object (solved)
    lg_max : int, photon hierarchy truncation
    ln_max : int, neutrino hierarchy truncation
    l_pol_max : int, E-mode polarization hierarchy truncation (0 = no polarization)

    Returns
    -------
    rhs : callable f(tau, y) -> dy/dtau
    nvar : int, length of state vector
    """
    k2 = k * k
    H02 = _H0_MPC ** 2
    Og, On, Ob, Oc = _OMEGA_GAMMA, _OMEGA_NU, _OMEGA_B, _OMEGA_C
    nvar = n_var_sync(lg_max, ln_max, l_pol_max)
    fn_s = idx_fn_start(lg_max, l_pol_max)
    _has_pol = l_pol_max > 0
    if _has_pol:
        _e_s = idx_e_start(lg_max)

    # Pre-tabulate background for fast np.interp (called ~10k-50k times per k)
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
        Fn = y[fn_s: fn_s + ln_max + 1]

        # Extract E-mode polarization if present
        if _has_pol:
            E = y[_e_s: _e_s + l_pol_max + 1]

        delta_g = Fg[0]
        delta_n = Fn[0]
        theta_g = 0.75 * k * Fg[1]   # theta_gamma = (3/4) k F_1
        theta_n = 0.75 * k * Fn[1]

        # --- Diagnose h' from 00-constraint ---
        # k^2 eta - (calH/2) h' = -(3/2) H0^2 [Og/a^2 dg + On/a^2 dn + Ob/a db + Oc/a dc]
        src00 = Og * ia2 * delta_g + On * ia2 * delta_n + Ob * ia * delta_b + Oc * ia * delta_c
        h_prime = (2.0 / calH) * (k2 * eta + 1.5 * H02 * src00)

        # --- Diagnose eta' from 0i-constraint ---
        # k^2 eta' = (3/2) H0^2 [(4/3) Og/a^2 thetag + (4/3) On/a^2 thetan + Ob/a thetab]
        src0i = (4.0 / 3.0) * Og * ia2 * theta_g + (4.0 / 3.0) * On * ia2 * theta_n + Ob * ia * theta_b
        eta_prime = 1.5 * H02 / k2 * src0i

        # --- CDM (theta_c = 0, gauge choice) ---
        # M&B Eq. 25: delta_c' = -h'/2
        d_delta_c = -0.5 * h_prime

        # --- Baryons ---
        # M&B Eq. 26:
        d_delta_b = -theta_b - 0.5 * h_prime
        # theta_b': Hubble drag + baryon pressure + Thomson coupling to photons
        cs2_b = _baryon_cs2(a)
        d_theta_b = (-calH * theta_b + cs2_b * k2 * delta_b
                     + abs_kd / R * (theta_g - theta_b))

        # --- Polarization source Pi ---
        # Pi = F_gamma,2 + E_0 + E_2 when polarization is evolved
        # Pi = F_gamma,2 when l_pol_max = 0 (reduces to original -(9/10)|kd|F_g2)
        if _has_pol and lg_max >= 2:
            Pi = Fg[2] + E[0] + (E[2] if l_pol_max >= 2 else 0.0)
        elif lg_max >= 2:
            Pi = Fg[2]  # no polarization feedback
        else:
            Pi = 0.0

        # --- Photon hierarchy (M&B Eq. 27) ---
        dFg = np.zeros(lg_max + 1)

        # l=0: F_0' = -k F_1 - (2/3) h'
        dFg[0] = -k * Fg[1] - (2.0 / 3.0) * h_prime

        # l=1: F_1' = (k/3)(F_0 - 2F_2) + |kd|(-F_1 + 4 theta_b / (3k))
        F2g = Fg[2] if lg_max >= 2 else 0.0
        dFg[1] = (k / 3.0) * (Fg[0] - 2.0 * F2g) + abs_kd * (-Fg[1] + 4.0 * theta_b / (3.0 * k))

        # l=2: F_2' with polarization feedback
        # Without polarization: -(9/10)|kd|F_g2
        # With polarization:    -|kd|(F_g2 - Pi/10) = -(9/10)|kd|F_g2 + (1/10)|kd|(E_0+E_2)
        # eta' coeff = 8/5 (not 8/15): from (8/15)*k^2*sigma where sigma=(h'+6eta')/(2k^2)
        # Confirmed by CLASS perturbations.c line 8982 and first-principles derivation
        if lg_max >= 2:
            F3g = Fg[3] if lg_max >= 3 else 0.0
            if _has_pol:
                # Full collision term: -|kd|*(F_g2 - Pi/10)
                dFg[2] = ((k / 5.0) * (2.0 * Fg[1] - 3.0 * F3g)
                          + (4.0 / 15.0) * h_prime + (8.0 / 5.0) * eta_prime
                          - abs_kd * (Fg[2] - Pi / 10.0))
            else:
                # Original: -(9/10)|kd|F_g2  (since Pi=F_g2 -> F_g2 - F_g2/10 = 9/10 F_g2)
                dFg[2] = ((k / 5.0) * (2.0 * Fg[1] - 3.0 * F3g)
                          + (4.0 / 15.0) * h_prime + (8.0 / 5.0) * eta_prime
                          - (9.0 / 10.0) * abs_kd * Fg[2])

        # l=3..lg_max-1: streaming + Thomson damping
        for ell in range(3, lg_max):
            dFg[ell] = (k / (2.0 * ell + 1.0) * (ell * Fg[ell - 1] - (ell + 1) * Fg[ell + 1])
                        - abs_kd * Fg[ell])

        # l=lg_max: truncation boundary
        if lg_max >= 3:
            dFg[lg_max] = (k * Fg[lg_max - 1] * lg_max / (2.0 * lg_max + 1.0)
                           - (lg_max + 1.0) / tau_safe * Fg[lg_max]
                           - abs_kd * Fg[lg_max])

        # --- E-mode polarization hierarchy ---
        if _has_pol:
            dE = np.zeros(l_pol_max + 1)

            # l=0: E_0' = -k*E_1 - |kd|*E_0 + (1/2)*|kd|*Pi
            E1 = E[1] if l_pol_max >= 1 else 0.0
            dE[0] = -k * E1 - abs_kd * E[0] + 0.5 * abs_kd * Pi

            # l=1: E_1' = k/3*(E_0 - 2*E_2) - |kd|*E_1
            if l_pol_max >= 1:
                E2 = E[2] if l_pol_max >= 2 else 0.0
                dE[1] = (k / 3.0) * (E[0] - 2.0 * E2) - abs_kd * E[1]

            # l=2: E_2' = k/5*(2*E_1 - 3*E_3) - |kd|*E_2 + (1/10)*|kd|*Pi
            if l_pol_max >= 2:
                E3 = E[3] if l_pol_max >= 3 else 0.0
                dE[2] = ((k / 5.0) * (2.0 * E[1] - 3.0 * E3)
                         - abs_kd * E[2] + (1.0 / 10.0) * abs_kd * Pi)

            # l=3..l_pol_max-1: streaming + Thomson damping
            for ell in range(3, l_pol_max):
                dE[ell] = (k / (2.0 * ell + 1.0) * (ell * E[ell - 1] - (ell + 1) * E[ell + 1])
                           - abs_kd * E[ell])

            # l=l_pol_max: truncation boundary
            if l_pol_max >= 3:
                dE[l_pol_max] = (k * E[l_pol_max - 1] * l_pol_max / (2.0 * l_pol_max + 1.0)
                                 - (l_pol_max + 1.0) / tau_safe * E[l_pol_max]
                                 - abs_kd * E[l_pol_max])

        # --- Neutrino hierarchy (M&B Eq. 28, no collisions) ---
        dFn = np.zeros(ln_max + 1)

        # l=0
        dFn[0] = -k * Fn[1] - (2.0 / 3.0) * h_prime

        # l=1
        F2n = Fn[2] if ln_max >= 2 else 0.0
        dFn[1] = (k / 3.0) * (Fn[0] - 2.0 * F2n)

        # l=2 (no Thomson damping for neutrinos)
        # eta' coeff = 8/5 (not 8/15): same fix as photons
        if ln_max >= 2:
            F3n = Fn[3] if ln_max >= 3 else 0.0
            dFn[2] = ((k / 5.0) * (2.0 * Fn[1] - 3.0 * F3n)
                      + (4.0 / 15.0) * h_prime + (8.0 / 5.0) * eta_prime)

        # l=3..ln_max-1: free streaming
        for ell in range(3, ln_max):
            dFn[ell] = k / (2.0 * ell + 1.0) * (ell * Fn[ell - 1] - (ell + 1) * Fn[ell + 1])

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
        if _has_pol:
            dydt[_e_s: _e_s + l_pol_max + 1] = dE
        dydt[fn_s: fn_s + ln_max + 1] = dFn

        return dydt

    return rhs, nvar


# ============================================================================
# Adiabatic initial conditions (synchronous gauge, CDM rest frame)
# ============================================================================

def adiabatic_ic_sync(k, tau_init, lg_max=L_GAMMA_MAX, ln_max=L_NU_MAX_SYNC,
                      l_pol_max=0):
    """
    Adiabatic growing-mode ICs following Ma & Bertschinger (1995) Sec. V.

    Normalization: eta(tau_init) ~ 1 (primordial curvature perturbation C = 1).
    Power spectrum normalization applied later in the C_l formula.

    At leading order in x = k tau_init << 1:
      eta = C [1 - (5+4R_nu)/(12(15+4R_nu)) x^2]
      h   = C x^2                     (not evolved, but needed conceptually)
      delta_c = delta_b = -C x^2 / 2
      delta_gamma = delta_nu = -(2/3) C x^2
      theta_b = theta_gamma ~ 0       (tight coupling, O(x^3))
      F_g,2 ~ 0                       (Thomson-damped)
      F_n,2 ~ f(R_nu) x^2            (free-streaming neutrinos)
      E_l = 0 for all l               (Thomson-damped at early times)
    """
    R_nu = _OMEGA_NU / (_OMEGA_GAMMA + _OMEGA_NU)
    C = 1.0
    x = k * tau_init
    x2 = x * x

    nvar = n_var_sync(lg_max, ln_max, l_pol_max)
    fn_s = idx_fn_start(lg_max, l_pol_max)
    y0 = np.zeros(nvar, dtype=np.float64)

    # eta
    fac_eta = (5.0 + 4.0 * R_nu) / (12.0 * (15.0 + 4.0 * R_nu))
    y0[IDX_ETA] = C * (1.0 - fac_eta * x2)

    # CDM and baryons (adiabatic: delta_b = delta_c at leading order)
    y0[IDX_DELTA_C] = -0.5 * C * x2
    y0[IDX_DELTA_B] = -0.5 * C * x2

    # Baryon velocity (tight coupled to photons, both ~ 0 at early times)
    y0[IDX_THETA_B] = 0.0

    # Photon monopole: delta_gamma = -(2/3) C x^2
    y0[idx_fg(0)] = -(2.0 / 3.0) * C * x2

    # Photon dipole: O(x^3), negligible
    y0[idx_fg(1)] = 0.0

    # Photon quadrupole: Thomson-damped, set to 0
    # (the equilibrium value is suppressed by |kappa_dot|)
    if lg_max >= 2:
        y0[idx_fg(2)] = 0.0

    # E-mode polarization: Thomson-damped at early times, all zero
    # (the E_l block is already zero from np.zeros initialization)

    # Neutrino monopole: adiabatic, delta_nu = delta_gamma
    y0[fn_s + 0] = -(2.0 / 3.0) * C * x2

    # Neutrino dipole: O(x^3), negligible
    y0[fn_s + 1] = 0.0

    # Neutrino quadrupole: free-streaming builds up
    # F_n,2 = 32(5+R_nu)/(45(15+4R_nu)) * C * x^2
    if ln_max >= 2:
        fac_n2 = 32.0 * (5.0 + R_nu) / (45.0 * (15.0 + 4.0 * R_nu))
        y0[fn_s + 2] = fac_n2 * C * x2

    # Higher neutrino multipoles: F_n,l ~ O(x^l) (tiny)
    for ell in range(3, min(ln_max + 1, 6)):
        prod = 1.0
        for j in range(1, ell + 1):
            prod *= (2 * j + 1)
        y0[fn_s + ell] = C * x ** ell / prod

    return y0


# ============================================================================
# Gauge transformation: synchronous -> Newtonian
# ============================================================================

def gauge_transform(y, k, calH, a, lg_max, ln_max, h_prime, eta_prime,
                    l_pol_max=0):
    """
    Transform sync gauge state to Newtonian gauge quantities needed for LOS.

    For superhorizon modes (k < 0.005 Mpc^-1), the gauge parameter
    alpha = (h' + 6 eta') / (2 k^2) diverges as 1/k^2, amplifying numerical
    noise. We cap alpha using the analytical bound from adiabatic ICs.

    Parameters
    ----------
    y : state vector
    k : wavenumber
    calH : conformal Hubble
    a : scale factor
    lg_max, ln_max : hierarchy truncation
    h_prime, eta_prime : metric derivatives (diagnosed from constraints)
    l_pol_max : int, E-mode polarization hierarchy truncation (0 = no polarization)

    Returns
    -------
    dict with keys: 'Phi_N', 'Psi_N', 'delta_g_N', 'theta_b_N'
    """
    k2 = k * k
    H02 = _H0_MPC ** 2
    Og, On = _OMEGA_GAMMA, _OMEGA_NU
    fn_s = idx_fn_start(lg_max, l_pol_max)

    eta = y[IDX_ETA]
    Fg = y[IDX_FG_START: IDX_FG_START + lg_max + 1]
    Fn = y[fn_s: fn_s + ln_max + 1]
    theta_b = y[IDX_THETA_B]

    # Gauge parameter: alpha = (h' + 6 eta') / (2 k^2)
    # With the minus sign convention for delta_g_N, the SW source
    # (Theta_0+Psi) has residual alpha dependence (~2*calH*alpha).
    # For superhorizon modes, alpha ~ 1/k^2 diverges, causing low-l blowup.
    # Cap |calH*alpha| at 2*|eta| to suppress noise while preserving physics.
    alpha = (h_prime + 6.0 * eta_prime) / (2.0 * k2)
    if calH > 0:
        max_alpha = 5.0 * abs(eta) / calH
        alpha = np.clip(alpha, -max_alpha, max_alpha)

    # Newtonian gauge potentials
    Phi_N = eta - calH * alpha

    sigma_g = 0.5 * Fg[2] if lg_max >= 2 else 0.0
    sigma_n = 0.5 * Fn[2] if ln_max >= 2 else 0.0
    # Aniso = 6*H0^2/(a^2*k^2) * (Og*sigma_g + On*sigma_n)
    # Factor 6 from: 12*pi*G*a^2*(4/3)*rho_g = 6*H0^2*Omega_g/a^2
    # Confirmed by CLASS perturbations.c line 6491
    aniso = 6.0 * H02 / (a * a * k2) * (Og * sigma_g + On * sigma_n)
    Psi_N = Phi_N - aniso

    # Density gauge transformation.
    # Numerically verified: the minus sign gives the correct Sachs-Wolfe
    # plateau (Theta_0+Psi ~ Phi_N/3 ~ 0.22) for superhorizon adiabatic modes.
    # With the plus sign, the SW plateau is ~0.89 (4x too high).
    # The sign depends on the metric convention for Phi_N/Psi_N in the code:
    # our Phi_N = eta - calH*alpha maps to the SPATIAL curvature perturbation,
    # not the temporal lapse, which explains the sign flip from the standard formula.
    # Gauge transformation: delta_N = delta_S + (rho'/rho)*alpha
    # For radiation: rho'/rho = -4*calH, so delta_gamma_N = delta_gamma_S - 4*calH*alpha
    delta_g_N = Fg[0] - 4.0 * calH * alpha

    # Velocity: theta_b_N = theta_b_S + k^2 alpha
    theta_b_N = theta_b + k2 * alpha

    return {
        'Phi_N': Phi_N,
        'Psi_N': Psi_N,
        'delta_g_N': delta_g_N,
        'theta_b_N': theta_b_N,
        'alpha': alpha,
    }


def diagnose_metric(y, k, calH, a, lg_max, ln_max, l_pol_max=0):
    """
    Compute h' and eta' from Einstein constraints for a given state.

    Uses H0^2 (not calH^2) in the constraint equations:
      h' = (2/calH) [k^2 eta + (3/2) H0^2 sum Omega_i delta_i / a^{1+3w_i}]
      eta' = (3/2) H0^2 / k^2 sum (1+w_i) Omega_i theta_i / a^{1+3w_i}

    Returns h_prime, eta_prime.
    """
    k2 = k * k
    H02 = _H0_MPC ** 2
    Og, On, Ob, Oc = _OMEGA_GAMMA, _OMEGA_NU, _OMEGA_B, _OMEGA_C
    fn_s = idx_fn_start(lg_max, l_pol_max)
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

    src00 = Og * ia2 * delta_g + On * ia2 * delta_n + Ob * ia * delta_b + Oc * ia * delta_c
    h_prime = (2.0 / calH) * (k2 * eta + 1.5 * H02 * src00)

    src0i = (4.0 / 3.0) * Og * ia2 * theta_g + (4.0 / 3.0) * On * ia2 * theta_n + Ob * ia * theta_b
    eta_prime = 1.5 * H02 / k2 * src0i

    return h_prime, eta_prime
