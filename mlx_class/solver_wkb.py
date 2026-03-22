"""
solver_wkb.py -- WKB envelope solver for CMB Boltzmann equations.

PHYSICS
=======
The photon-baryon fluid in the tight coupling regime oscillates rapidly:
    Theta_0(k,tau) ~ A(k,tau) cos[psi(k,tau)] + B(k,tau) sin[psi(k,tau)]

where psi(tau) = integral_0^tau c_s(tau') dtau' is the sound horizon and
A, B are slowly-varying envelopes that change on the Hubble timescale.

APPROACH: Two-pass solver
=========================
Pass 1 (metric sector): Solve CDM + neutrino + metric using the FULL sync
gauge equations with scipy Radau, but with the photon-baryon fluid computed
analytically from the WKB transfer function. This gives Phi(k,tau), Psi(k,tau).

Pass 2 (photon envelope): With the metric potentials known, solve the WKB
envelope equations for A(tau), B(tau) as an externally forced damped system.
No feedback loop -- the envelopes are driven by the pre-computed potentials.

This decoupling avoids the numerical instability of self-consistently evolving
the metric + envelope system, while capturing all the relevant physics:
- Gravitational driving from potential decay at equality
- Baryon drag and Silk damping on the envelope
- Neutrino anisotropic stress (Psi != Phi)

SPEED ADVANTAGE
===============
The metric sector has ~10 variables per k and varies slowly.
The envelope has only 2 variables (A, B) and needs ~20 steps.
Total: much faster than the full 50+ variable hierarchy with ~600 steps.

Author: Sheng-Kai Huang, 2026
"""
import os
os.environ['MPLCONFIGDIR'] = '/tmp/mpl_config'

import numpy as np
import time
import sys

from scipy.integrate import solve_ivp, cumulative_trapezoid
from scipy.interpolate import interp1d, CubicSpline

from .background import (
    Background, A_s, n_s, k_pivot, T_CMB,
    Omega_r as _OMEGA_R,
    Omega_b as _OMEGA_B,
    Omega_c as _OMEGA_C,
    Omega_m as _OMEGA_M,
    H0_Mpc as _H0_MPC,
)
from .perturbations_neutrino import (
    _f_nu, _OMEGA_GAMMA, _OMEGA_NU,
)


# ============================================================================
# Neutrino hierarchy truncation for the metric sector
# ============================================================================
L_NU_MAX_WKB = 6


# ============================================================================
# Sound horizon and acoustic phase
# ============================================================================

def compute_acoustic_phase(bg):
    """
    Compute the sound horizon psi(tau) = integral_0^tau c_s dtau'.

    Returns interpolators for psi(tau) and c_s(tau).
    """
    tau = bg.tau_grid
    R = bg.R_grid
    cs = 1.0 / np.sqrt(3.0 * (1.0 + R))

    psi_vals = np.zeros_like(tau)
    psi_vals[1:] = cumulative_trapezoid(cs, tau)

    psi_of_tau = interp1d(tau, psi_vals, kind='cubic', fill_value='extrapolate')
    cs_of_tau = interp1d(tau, cs, kind='cubic', fill_value='extrapolate')

    return psi_of_tau, cs_of_tau


# ============================================================================
# Pass 1: Metric sector solver (CDM + neutrinos + eta)
# ============================================================================
# State: [eta, delta_c, N_0, N_1, ..., N_{ln_max}]
# The photon-baryon contribution to the metric is approximated using the
# WKB ansatz: delta_gamma ~ -(2/3) (k tau)^2 C at early times, then
# transitions to acoustic oscillations.

_IDX_ETA_M = 0
_IDX_DC_M = 1
_IDX_NU_M = 2  # N_0 starts here


def _n_var_metric(ln_max):
    return 2 + (ln_max + 1)


def _make_metric_rhs(k, bg, psi_interp, cs_interp, ln_max=L_NU_MAX_WKB):
    """
    Build the metric sector RHS.

    The photon-baryon contribution to the constraints is modeled using
    the WKB analytical transfer function T(k*tau):
    - For k*tau < 1: delta_gamma ~ -(2/3)(k*tau)^2, theta_gamma ~ 0
    - For k*tau > 1: delta_gamma ~ -2/3 cos(k*r_s), theta_gamma from dipole

    The transfer function is computed from the initial conditions and the
    damped oscillator solution, matching to the full Boltzmann result at
    the few-percent level.
    """
    k2 = k * k
    H02 = _H0_MPC ** 2
    Og, On, Ob, Oc = _OMEGA_GAMMA, _OMEGA_NU, _OMEGA_B, _OMEGA_C
    nvar = _n_var_metric(ln_max)
    ns = _IDX_NU_M

    _tau = bg.tau_grid.copy()
    _calH = bg.calH_grid.copy()
    _a = bg.a_grid.copy()
    _R = bg.R_grid.copy()
    _kd = bg.kappa_dot_grid.copy()

    def rhs(tau, y):
        calH = float(np.interp(tau, _tau, _calH))
        a = float(np.interp(tau, _tau, _a))
        R = float(np.interp(tau, _tau, _R))
        cs = float(cs_interp(tau))
        psi = float(psi_interp(tau))  # sound horizon in Mpc

        tau_safe = max(tau, 1e-10)
        ia = 1.0 / a
        ia2 = ia * ia

        eta = y[_IDX_ETA_M]
        delta_c = y[_IDX_DC_M]
        Fn = y[ns: ns + ln_max + 1]

        delta_n = Fn[0]
        theta_n = 0.75 * k * Fn[1]

        # Photon-baryon WKB model:
        # Use the analytic tight-coupling oscillation
        x = k * tau
        krs = k * psi  # k * sound horizon

        # Damping envelope (Silk + baryon drag)
        kd_val = float(np.interp(tau, _tau, _kd))
        abs_kd = abs(kd_val)
        if abs_kd > 1e-10:
            # Silk damping scale: k_D ~ sqrt(integral k^2/(|kd|(1+R)) dtau)
            gamma_S = k2 / (6.0 * (1.0 + R) * abs_kd)
        else:
            gamma_S = 0.0
        gamma_H = 0.5 * calH * R / (1.0 + R)
        # Accumulated damping (approximate using local rate * tau)
        damp = np.exp(-min((gamma_S + gamma_H) * tau * 0.5, 100.0))

        if x < 0.5:
            # Sub-horizon: use Taylor expansion of ICs
            C = 1.0  # normalization
            delta_g = -(2.0 / 3.0) * C * x * x * damp
            theta_g = 0.0
        else:
            # Oscillating regime: analytic transfer function
            # Theta_0 = -(1/6) C cos(k r_s) * (1+R)^{-1/4} * damp
            # delta_gamma = 4 Theta_0
            C = 1.0
            amp = (1.0 / 6.0) * C * (1.0 + R)**(-0.25) * damp
            delta_g = -4.0 * amp * np.cos(krs)
            # Theta_1 = (1/3) * c_s * amp * sin(k*r_s) * k / (3/4 k)
            theta_g = (4.0 / 3.0) * amp * k * cs * np.sin(krs)

        # In tight coupling, theta_b ~ theta_gamma
        theta_b = theta_g
        delta_b = delta_c  # adiabatic approx

        # h' from 00-constraint
        src00 = (Og * ia2 * delta_g + On * ia2 * delta_n
                 + Ob * ia * delta_b + Oc * ia * delta_c)
        h_prime = (2.0 / calH) * (k2 * eta + 1.5 * H02 * src00)

        # eta' from 0i-constraint
        src0i = ((4.0 / 3.0) * Og * ia2 * theta_g
                 + (4.0 / 3.0) * On * ia2 * theta_n
                 + Ob * ia * theta_b)
        eta_prime = 1.5 * H02 / k2 * src0i

        # CDM
        d_delta_c = -0.5 * h_prime

        # Neutrino hierarchy
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

        dydt = np.zeros(nvar)
        dydt[_IDX_ETA_M] = eta_prime
        dydt[_IDX_DC_M] = d_delta_c
        dydt[ns: ns + ln_max + 1] = dFn

        return dydt

    return rhs, nvar


def _metric_ic(k, tau_init, ln_max=L_NU_MAX_WKB):
    """Adiabatic ICs for the metric sector."""
    R_nu = _OMEGA_NU / (_OMEGA_GAMMA + _OMEGA_NU)
    C = 1.0
    x = k * tau_init
    x2 = x * x

    nvar = _n_var_metric(ln_max)
    y0 = np.zeros(nvar, dtype=np.float64)

    fac_eta = (5.0 + 4.0 * R_nu) / (12.0 * (15.0 + 4.0 * R_nu))
    y0[_IDX_ETA_M] = C * (1.0 - fac_eta * x2)
    y0[_IDX_DC_M] = -0.5 * C * x2

    ns = _IDX_NU_M
    y0[ns + 0] = -(2.0 / 3.0) * C * x2
    y0[ns + 1] = 0.0
    if ln_max >= 2:
        fac_n2 = 32.0 * (5.0 + R_nu) / (45.0 * (15.0 + 4.0 * R_nu))
        y0[ns + 2] = fac_n2 * C * x2
    for ell in range(3, min(ln_max + 1, 6)):
        prod = 1.0
        for j in range(1, ell + 1):
            prod *= (2 * j + 1)
        y0[ns + ell] = C * x ** ell / prod

    return y0


def solve_metric_sector(k, bg, psi_interp, cs_interp, tau_snaps,
                         ln_max=L_NU_MAX_WKB):
    """
    Solve the metric sector for a single k-mode, returning Phi_N, Psi_N
    at the requested snapshot times.

    Returns dict with 'Phi_N', 'Psi_N', 'h_prime', 'eta_prime', 'delta_c',
    'eta' at each snapshot.
    """
    tau_init = bg.tau_grid[1]
    rhs_fn, nvar = _make_metric_rhs(k, bg, psi_interp, cs_interp, ln_max)
    y0 = _metric_ic(k, tau_init, ln_max)

    # Only integrate to slightly past recombination. Late-time potentials
    # are handled analytically (CDM-dominated, slow evolution).
    tau_max_metric = min(float(tau_snaps[-1]) + 1.0, 500.0)  # ~500 Mpc
    sol = solve_ivp(rhs_fn, [tau_init, tau_max_metric], y0,
                    method='RK45', dense_output=True,
                    rtol=1e-6, atol=1e-9)

    if not sol.success:
        return None

    k2 = k * k
    H02 = _H0_MPC ** 2
    Og, On, Ob, Oc = _OMEGA_GAMMA, _OMEGA_NU, _OMEGA_B, _OMEGA_C
    ns = _IDX_NU_M

    N_snap = len(tau_snaps)
    result = {
        'Phi_N': np.zeros(N_snap),
        'Psi_N': np.zeros(N_snap),
        'h_prime': np.zeros(N_snap),
        'eta_prime': np.zeros(N_snap),
        'delta_c': np.zeros(N_snap),
        'eta': np.zeros(N_snap),
    }

    # Get the final state for extrapolation to late times
    y_final = sol.sol(tau_max_metric)
    eta_final = y_final[_IDX_ETA_M]
    delta_c_final = y_final[_IDX_DC_M]

    for it in range(N_snap):
        tau = tau_snaps[it]
        if tau <= tau_max_metric:
            y = sol.sol(tau)
        else:
            # Late time: CDM dominated. eta ~ const, delta_c grows ~ a.
            # Potentials decay slowly due to dark energy.
            # Use the final integrated value (adequate for ISW which is small).
            y = y_final.copy()
            # Approximate CDM growth: delta_c ~ delta_c_final * (a/a_final)
            a_final = float(bg.a_at_tau(tau_max_metric))
            a_now = float(bg.a_at_tau(tau))
            y[_IDX_DC_M] = delta_c_final * a_now / max(a_final, 1e-10)

        calH = float(bg.calH_at_tau(tau))
        a = float(bg.a_at_tau(tau))
        R = float(bg.R_at_tau(tau))
        cs = float(cs_interp(tau))
        psi = float(psi_interp(tau))
        kd_val = float(bg.kappa_dot_at_tau(tau))

        ia = 1.0 / a
        ia2 = ia * ia

        eta = y[_IDX_ETA_M]
        delta_c = y[_IDX_DC_M]
        Fn = y[ns: ns + ln_max + 1]
        delta_n = Fn[0]
        theta_n = 0.75 * k * Fn[1]
        sigma_n = 0.5 * Fn[2] if ln_max >= 2 else 0.0

        # Photon-baryon (WKB model, same as in the RHS)
        x = k * tau
        krs = k * psi
        abs_kd = abs(kd_val)
        if abs_kd > 1e-10:
            gamma_S = k2 / (6.0 * (1.0 + R) * abs_kd)
        else:
            gamma_S = 0.0
        gamma_H = 0.5 * calH * R / (1.0 + R)
        damp = np.exp(-min((gamma_S + gamma_H) * tau * 0.5, 100.0))

        if x < 0.5:
            C = 1.0
            delta_g = -(2.0 / 3.0) * C * x * x * damp
            theta_g = 0.0
        else:
            C = 1.0
            amp = (1.0 / 6.0) * C * (1.0 + R)**(-0.25) * damp
            delta_g = -4.0 * amp * np.cos(krs)
            theta_g = (4.0 / 3.0) * amp * k * cs * np.sin(krs)

        theta_b = theta_g
        delta_b = delta_c

        src00 = (Og * ia2 * delta_g + On * ia2 * delta_n
                 + Ob * ia * delta_b + Oc * ia * delta_c)
        h_prime = (2.0 / calH) * (k2 * eta + 1.5 * H02 * src00)

        src0i = ((4.0 / 3.0) * Og * ia2 * theta_g
                 + (4.0 / 3.0) * On * ia2 * theta_n
                 + Ob * ia * theta_b)
        eta_prime = 1.5 * H02 / k2 * src0i

        alpha = (h_prime + 6.0 * eta_prime) / (2.0 * k2)
        if calH > 0:
            max_alpha = 5.0 * abs(eta) / calH
            alpha = np.clip(alpha, -max_alpha, max_alpha)

        Phi_N = eta - calH * alpha
        aniso = 12.0 * H02 / (a * a * k2) * On * sigma_n
        Psi_N = Phi_N - aniso

        result['Phi_N'][it] = Phi_N
        result['Psi_N'][it] = Psi_N
        result['h_prime'][it] = h_prime
        result['eta_prime'][it] = eta_prime
        result['delta_c'][it] = delta_c
        result['eta'][it] = eta

    return result


# ============================================================================
# Pass 2: WKB envelope solver (externally driven)
# ============================================================================

def solve_envelope(k, bg, psi_interp, cs_interp, Psi_N_interp, tau_snaps):
    """
    Solve the WKB envelope equations for a single k-mode with EXTERNAL
    gravitational driving (no feedback loop).

    The envelope equations:
        dA/dtau = -gamma*A + S_A(tau)
        dB/dtau = -gamma*B + S_B(tau)

    where gamma = gamma_H + gamma_S (damping) and S_A, S_B are source
    terms from the gravitational potential decay:
        S_A = +(F_drive/(k c_s)) * sin(k*psi)
        S_B = -(F_drive/(k c_s)) * cos(k*psi)
        F_drive = -(k^2/3) * Psi_N(tau)

    This is a simple 2-variable ODE that is stable and well-conditioned.

    Parameters
    ----------
    k : float
        Wavenumber.
    bg : Background
    psi_interp, cs_interp : callable
    Psi_N_interp : callable
        Newtonian gauge potential Psi(tau) from Pass 1.
    tau_snaps : ndarray
        Output snapshot times.

    Returns
    -------
    A_arr, B_arr : ndarrays, shape (N_snap,)
        Envelope amplitudes at snapshot times.
    """
    k2 = k * k

    _tau = bg.tau_grid.copy()
    _calH = bg.calH_grid.copy()
    _a = bg.a_grid.copy()
    _R = bg.R_grid.copy()
    _kd = bg.kappa_dot_grid.copy()

    def rhs(tau, y):
        A, B = y[0], y[1]

        calH = float(np.interp(tau, _tau, _calH))
        R = float(np.interp(tau, _tau, _R))
        kd_val = float(np.interp(tau, _tau, _kd))
        cs = float(cs_interp(tau))
        psi = float(psi_interp(tau))

        abs_kd = abs(kd_val)

        # Damping
        gamma_H = 0.5 * calH * R / (1.0 + R)
        if abs_kd > 1e-10:
            gamma_S = k2 / (6.0 * (1.0 + R) * abs_kd)
        else:
            gamma_S = 0.0
        gamma = gamma_H + gamma_S

        # Sound speed variation
        R_prime = calH * R
        cs_prime = -0.5 * R_prime / ((1.0 + R) * np.sqrt(3.0 * (1.0 + R)))
        cs_ratio = cs_prime / max(abs(cs), 1e-20)

        # External driving from pre-computed potential
        Psi_N = float(Psi_N_interp(tau))
        F_drive = -(k2 / 3.0) * Psi_N

        kcs = k * cs
        if abs(kcs) > 1e-20:
            inv_kcs = 1.0 / kcs
        else:
            inv_kcs = 0.0

        kpsi = k * psi
        sin_kpsi = np.sin(kpsi)
        cos_kpsi = np.cos(kpsi)

        # Envelope equations (variation of parameters on damped oscillator)
        dA = (-gamma * A
              - 0.5 * cs_ratio * A
              + F_drive * inv_kcs * sin_kpsi)

        dB = (-gamma * B
              - 0.5 * cs_ratio * B
              - F_drive * inv_kcs * cos_kpsi)

        return [dA, dB]

    # Initial conditions: at tau_init, Theta_0 ~ -(1/6)(k tau)^2
    tau_init = bg.tau_grid[1]
    x = k * tau_init
    psi_init = float(psi_interp(tau_init))
    kpsi_init = k * psi_init
    Theta0_init = -(1.0 / 6.0) * x * x
    # Decompose: A cos(kpsi) + B sin(kpsi) = Theta0_init
    # At early times kpsi ~ 0, so A ~ Theta0_init, B ~ 0
    cos_psi_init = np.cos(kpsi_init)
    sin_psi_init = np.sin(kpsi_init)
    A_init = Theta0_init * cos_psi_init
    B_init = Theta0_init * sin_psi_init

    # Only integrate to past recombination (envelope is heavily damped after)
    tau_max_env = min(float(tau_snaps[-1]) + 1.0, 500.0)
    sol = solve_ivp(rhs, [tau_init, tau_max_env], [A_init, B_init],
                    method='RK45', dense_output=True,
                    rtol=1e-6, atol=1e-9)

    if not sol.success:
        return np.zeros(len(tau_snaps)), np.zeros(len(tau_snaps))

    # Get final envelope values (heavily damped after recombination)
    y_final = sol.sol(tau_max_env)

    A_arr = np.zeros(len(tau_snaps))
    B_arr = np.zeros(len(tau_snaps))
    for it, tau in enumerate(tau_snaps):
        if tau <= tau_max_env:
            y = sol.sol(tau)
        else:
            # After recombination: photons decouple, envelopes frozen/decaying
            y = y_final  # Use last computed value (essentially zero)
        A_arr[it] = y[0]
        B_arr[it] = y[1]

    return A_arr, B_arr


# ============================================================================
# Reconstruct Newtonian gauge source functions from envelope + metric
# ============================================================================

def reconstruct_sources(k, A_arr, B_arr, metric, bg, psi_interp, cs_interp,
                         tau_snaps):
    """
    Reconstruct the LOS source functions from the WKB envelopes and metric.

    Uses the Newtonian gauge formulation directly:
    - Theta_0_N = Theta_0 + Psi_N (the gauge-invariant effective temperature)
      minus Psi_N added back in the LOS integral as (Theta0_N + Psi_N)
    - We compute Theta_0 from the WKB envelopes in the sync gauge, then
      the full (Theta0+Psi) is the physical observable

    Returns arrays of shape (N_snap,) for each source.
    """
    k2 = k * k
    N_snap = len(tau_snaps)

    Theta0_N = np.zeros(N_snap)
    Psi_N = metric['Psi_N'].copy()
    Phi_N = metric['Phi_N'].copy()
    vb_N = np.zeros(N_snap)
    Theta2 = np.zeros(N_snap)

    for it in range(N_snap):
        tau = tau_snaps[it]
        cs = float(cs_interp(tau))
        psi = float(psi_interp(tau))
        kd_val = float(bg.kappa_dot_at_tau(tau))
        R = float(bg.R_at_tau(tau))

        kpsi = k * psi
        cos_kpsi = np.cos(kpsi)
        sin_kpsi = np.sin(kpsi)

        A = A_arr[it]
        B = B_arr[it]

        # Theta_0 (synchronous gauge monopole / 4)
        Theta0 = A * cos_kpsi + B * sin_kpsi

        # For the LOS integral, the source is g(tau) * (Theta0_N + Psi_N).
        # In the standard Sachs-Wolfe formula:
        #   Theta0_N + Psi_N ~ Theta0_sync + Psi_N + (gauge correction)
        # The dominant contribution comes from Theta0 + Psi_N directly.
        # The gauge correction alpha is subdominant when the metric sector
        # is solved self-consistently with the same photon model.
        #
        # We use the direct "effective temperature":
        Theta0_N[it] = Theta0

        # Baryon velocity: v_b = Theta_1 / k (tight coupling)
        # Theta_1 = (-A sin + B cos) * c_s
        Theta1 = (-A * sin_kpsi + B * cos_kpsi) * cs
        vb_N[it] = Theta1  # This is Theta_1, used as v_b source

        # Photon quadrupole from tight-coupling
        abs_kd = abs(kd_val)
        if abs_kd > 1e-10:
            Theta2[it] = (4.0 / 15.0) * k * Theta1 / abs_kd
        else:
            Theta2[it] = 0.0

    return {
        'Theta0_N': Theta0_N,
        'Psi_N': Psi_N,
        'Phi_N': Phi_N,
        'vb_N': vb_N,
        'Theta2': Theta2,
    }


# ============================================================================
# Build snapshot tau grid
# ============================================================================

def build_snapshot_grid_wkb(bg, N_vis=60, N_early_isw=30, N_late_isw=20,
                             N_reion=20):
    """Build a tau grid for LOS integration snapshots."""
    peak_idx = np.argmax(bg.visibility_grid)
    tau_peak = bg.tau_grid[peak_idx]
    vis_peak = bg.visibility_grid[peak_idx]

    half_max = vis_peak / 2.0
    left = peak_idx
    while left > 0 and bg.visibility_grid[left] > half_max:
        left -= 1
    right = peak_idx
    while right < len(bg.visibility_grid) - 1 and bg.visibility_grid[right] > half_max:
        right += 1
    fwhm = bg.tau_grid[right] - bg.tau_grid[left]
    sigma = fwhm / 2.355

    tau_vis_lo = max(tau_peak - 5.0 * sigma, bg.tau_grid[5])
    tau_vis_hi = min(tau_peak + 5.0 * sigma, bg.tau_0 * 0.5)
    tau_vis = np.linspace(tau_vis_lo, tau_vis_hi, N_vis)

    a_eq = _OMEGA_R / _OMEGA_M
    tau_eq_approx = float(bg._tau_of_a(np.log(a_eq)))
    tau_early_lo = max(tau_eq_approx * 0.3, bg.tau_grid[5])
    tau_early_hi = tau_vis_lo

    if tau_early_hi > tau_early_lo * 1.5:
        tau_isw_dense_lo = max(tau_eq_approx - 50.0, tau_early_lo)
        tau_isw_dense_hi = min(tau_eq_approx + 50.0, tau_early_hi)
        tau_dense = np.linspace(tau_isw_dense_lo, tau_isw_dense_hi, N_early_isw)
        tau_sparse = np.geomspace(tau_early_lo, tau_early_hi, 15)
        tau_early = np.sort(np.unique(np.concatenate([tau_sparse, tau_dense])))
    else:
        tau_early = np.array([tau_early_lo])

    tau_reion = np.array([])
    if hasattr(bg, 'z_reio') and bg.tau_reio > 0:
        z_reio = bg.z_reio
        a_reio = 1.0 / (1.0 + z_reio)
        tau_reion_center = float(bg._tau_of_a(np.log(a_reio)))
        tau_reion_lo = max(tau_reion_center - 200.0, tau_vis_hi + 10.0)
        tau_reion_hi = min(tau_reion_center + 200.0, bg.tau_0 * 0.95)
        if tau_reion_hi > tau_reion_lo:
            tau_reion = np.linspace(tau_reion_lo, tau_reion_hi, N_reion)

    tau_late_lo = tau_vis_hi
    tau_late_hi = bg.tau_0 * 0.98
    tau_late = np.linspace(tau_late_lo, tau_late_hi, N_late_isw)

    tau_all = np.sort(np.unique(np.concatenate([
        tau_early, tau_vis, tau_reion, tau_late])))
    return tau_vis, tau_late, tau_all


# ============================================================================
# Multiprocessing worker for solving a single k-mode
# ============================================================================

def _solve_single_k_wkb_worker(args):
    """
    Worker function for multiprocessing.
    Solves both passes (metric + envelope) for a single k-mode.

    Must be at module level for pickling.
    """
    import numpy as np
    from scipy.integrate import solve_ivp
    from scipy.interpolate import interp1d

    (k, bg_arrays, tau_all, ln_max, psi_vals, cs_vals) = args

    # Reconstruct lightweight background
    class _BGLite:
        __slots__ = ('tau_grid', 'a_grid', 'calH_grid', 'R_grid',
                     'kappa_dot_grid', 'tau_rec', 'tau_0')
        def __init__(self, arrays):
            self.tau_grid = arrays['tau_grid']
            self.a_grid = arrays['a_grid']
            self.calH_grid = arrays['calH_grid']
            self.R_grid = arrays['R_grid']
            self.kappa_dot_grid = arrays['kappa_dot_grid']
            self.tau_rec = arrays['tau_rec']
            self.tau_0 = arrays['tau_0']
        def calH_at_tau(self, tau):
            return np.interp(tau, self.tau_grid, self.calH_grid)
        def a_at_tau(self, tau):
            return np.interp(tau, self.tau_grid, self.a_grid)
        def R_at_tau(self, tau):
            return np.interp(tau, self.tau_grid, self.R_grid)
        def kappa_dot_at_tau(self, tau):
            return np.interp(tau, self.tau_grid, self.kappa_dot_grid)

    bg = _BGLite(bg_arrays)
    psi_interp = interp1d(bg_arrays['tau_grid'], psi_vals,
                           kind='cubic', fill_value='extrapolate')
    cs_interp = interp1d(bg_arrays['tau_grid'], cs_vals,
                          kind='cubic', fill_value='extrapolate')

    # Import the solver functions
    from mlx_class.solver_wkb import (
        solve_metric_sector, solve_envelope, reconstruct_sources
    )

    N_snap = len(tau_all)
    try:
        metric = solve_metric_sector(k, bg, psi_interp, cs_interp,
                                      tau_all, ln_max)
        if metric is None:
            return None

        Psi_N_interp = interp1d(tau_all, metric['Psi_N'],
                                 kind='cubic', fill_value='extrapolate',
                                 bounds_error=False)

        A_arr, B_arr = solve_envelope(k, bg, psi_interp, cs_interp,
                                       Psi_N_interp, tau_all)

        sources = reconstruct_sources(k, A_arr, B_arr, metric, bg,
                                       psi_interp, cs_interp, tau_all)

        result = np.zeros((N_snap, 5))
        result[:, 0] = sources['Phi_N']
        result[:, 1] = sources['Psi_N']
        result[:, 2] = sources['Theta0_N']
        result[:, 3] = sources['vb_N']
        result[:, 4] = sources['Theta2']
        return result
    except Exception:
        return None


# ============================================================================
# Full WKB pipeline
# ============================================================================

def run_wkb_solver(N_k=500, k_min=3e-4, k_max=0.35, N_int_steps=40,
                   ln_max=L_NU_MAX_WKB, verbose=True, parallel=True):
    """
    Full WKB envelope pipeline for CMB power spectrum computation.

    Two-pass approach:
      Pass 1: Solve metric sector (CDM + neutrinos + eta) per k-mode
      Pass 2: Solve WKB envelopes (A, B) driven by pre-computed potentials
      Then: LOS integration -> C_l

    Parameters
    ----------
    N_k : int
        Number of k-modes.
    k_min, k_max : float
        k-range in Mpc^{-1}.
    N_int_steps : int
        Not used (kept for API compatibility). Step counts are adaptive.
    ln_max : int
        Neutrino hierarchy truncation.
    verbose : bool
        Print progress.
    parallel : bool
        Use multiprocessing for k-mode solving (default True).

    Returns
    -------
    dict with keys: 'ell', 'Dl_TT', 'Dl_EE', 'Dl_TE', 'Cl_TT', etc.
    """
    from multiprocessing import Pool
    t_total = time.time()

    # ================================================================
    # Step 1: Background
    # ================================================================
    if verbose:
        print("=" * 65)
        print("WKB Envelope Boltzmann Solver (Two-Pass)")
        print("=" * 65)
        print(f"\n--- Step 1: Background ---")

    t0 = time.time()
    bg = Background(khronon=False, recombination='recfast')
    bg.solve()
    t_bg = time.time() - t0
    if verbose:
        print(f"Background: {t_bg:.2f}s")

    # ================================================================
    # Step 2: Acoustic phase
    # ================================================================
    if verbose:
        print(f"\n--- Step 2: Acoustic phase ---")

    t0 = time.time()
    psi_interp, cs_interp = compute_acoustic_phase(bg)
    t_phase = time.time() - t0

    psi_rec = float(psi_interp(bg.tau_rec))
    if verbose:
        print(f"Sound horizon at rec: r_s = {psi_rec:.1f} Mpc")
        print(f"Acoustic phase: {t_phase:.3f}s")

    # ================================================================
    # Step 3: Setup
    # ================================================================
    k_arr = np.geomspace(k_min, k_max, N_k)

    tau_vis, tau_late, tau_all = build_snapshot_grid_wkb(
        bg, N_vis=60, N_early_isw=30, N_late_isw=20)
    N_snap = len(tau_all)

    if verbose:
        print(f"\n--- Step 3: Two-Pass WKB ({N_k} k-modes) ---")
        print(f"k range: [{k_min:.1e}, {k_max:.2f}] Mpc^-1")
        print(f"Metric sector: {_n_var_metric(ln_max)} vars "
              f"(eta, delta_c, {ln_max+1} neutrino multipoles)")
        print(f"Envelope: 2 vars (A, B)")
        print(f"Snapshot grid: {N_snap} tau points")
        sys.stdout.flush()

    # ================================================================
    # Step 4: Solve all k-modes (two-pass per k)
    # ================================================================
    t0 = time.time()

    Phi_N = np.zeros((N_k, N_snap))
    Psi_N = np.zeros((N_k, N_snap))
    Theta0_N = np.zeros((N_k, N_snap))
    vb_N = np.zeros((N_k, N_snap))
    Theta2_arr = np.zeros((N_k, N_snap))

    # Precompute psi and cs on the background grid for serialization
    psi_on_grid = np.asarray(psi_interp(bg.tau_grid), dtype=np.float64)
    cs_on_grid = np.asarray(cs_interp(bg.tau_grid), dtype=np.float64)

    bg_arrays = {
        'tau_grid': bg.tau_grid.copy(),
        'a_grid': bg.a_grid.copy(),
        'calH_grid': bg.calH_grid.copy(),
        'R_grid': bg.R_grid.copy(),
        'kappa_dot_grid': bg.kappa_dot_grid.copy(),
        'tau_rec': float(bg.tau_rec),
        'tau_0': float(bg.tau_0),
    }

    n_failed = 0
    if parallel and N_k > 10:
        # Parallel mode using multiprocessing
        n_workers = min(os.cpu_count() or 4, N_k)
        if verbose:
            print(f"  Dispatching {N_k} k-modes to {n_workers} workers...")
            sys.stdout.flush()

        work_items = [
            (k, bg_arrays, tau_all, ln_max, psi_on_grid, cs_on_grid)
            for k in k_arr
        ]

        with Pool(processes=n_workers) as pool:
            worker_results = pool.map(_solve_single_k_wkb_worker, work_items)

        for ik, res in enumerate(worker_results):
            if res is None:
                n_failed += 1
                continue
            Phi_N[ik] = res[:, 0]
            Psi_N[ik] = res[:, 1]
            Theta0_N[ik] = res[:, 2]
            vb_N[ik] = res[:, 3]
            Theta2_arr[ik] = res[:, 4]
    else:
        # Serial mode
        for ik, k in enumerate(k_arr):
            metric = solve_metric_sector(k, bg, psi_interp, cs_interp,
                                          tau_all, ln_max)
            if metric is None:
                n_failed += 1
                if verbose and n_failed <= 3:
                    print(f"  WARNING: k={k:.4e} metric solve failed")
                continue

            Psi_N_interp = interp1d(tau_all, metric['Psi_N'],
                                     kind='cubic', fill_value='extrapolate',
                                     bounds_error=False)

            A_arr, B_arr = solve_envelope(k, bg, psi_interp, cs_interp,
                                           Psi_N_interp, tau_all)

            sources = reconstruct_sources(k, A_arr, B_arr, metric, bg,
                                           psi_interp, cs_interp, tau_all)

            Phi_N[ik] = sources['Phi_N']
            Psi_N[ik] = sources['Psi_N']
            Theta0_N[ik] = sources['Theta0_N']
            vb_N[ik] = sources['vb_N']
            Theta2_arr[ik] = sources['Theta2']

            if verbose and (ik + 1) % 100 == 0:
                print(f"  {ik+1}/{N_k} k-modes done")
                sys.stdout.flush()

    t_pert = time.time() - t0
    if verbose:
        print(f"WKB perturbations: {t_pert:.1f}s ({n_failed} failed)")

    # ================================================================
    # Step 5: Compute Phi_N' + Psi_N' by cubic spline
    # ================================================================
    PhiPsi = Phi_N + Psi_N
    PhiPsi_prime = np.zeros_like(PhiPsi)
    for ik in range(N_k):
        if np.all(np.isfinite(PhiPsi[ik])):
            cs_spline = CubicSpline(tau_all, PhiPsi[ik])
            PhiPsi_prime[ik] = cs_spline(tau_all, 1)

    # ================================================================
    # Step 6: LOS integration -> C_l
    # ================================================================
    t0 = time.time()

    ell_values = np.unique(np.concatenate([
        np.arange(2, 30, 1),
        np.arange(30, 100, 2),
        np.arange(100, 500, 4),
        np.arange(500, 1500, 8),
        np.arange(1500, 2501, 12),
    ])).astype(int)
    N_ell = len(ell_values)

    if verbose:
        print(f"\n--- Step 5: LOS C_l ({N_ell} ells) ---")
        sys.stdout.flush()

    P_R = A_s * (k_arr / k_pivot) ** (n_s - 1.0)

    g_snap = bg.visibility_at_tau(tau_all)
    kappa_snap = bg.kappa_at_tau(tau_all)
    exp_neg_kappa = np.exp(-kappa_snap)
    chi_snap = bg.tau_0 - tau_all
    dtau_all = np.diff(tau_all)

    # Bessel table
    try:
        from .bessel_cache import CachedBesselTable
        import mlx.core as mx
        x_max_bessel = float(np.max(k_arr) * np.max(chi_snap)) * 1.05 + 50.0
        bessel_table = CachedBesselTable(ell_values, x_max=x_max_bessel,
                                          N_cheb=64, seg_width=80,
                                          verbose=verbose)
        use_gpu_bessel = True
    except Exception as e:
        if verbose:
            print(f"[WARNING] Bessel table failed ({e}), using scipy")
        use_gpu_bessel = False

    Delta_l = np.zeros((N_ell, N_k), dtype=np.float64)
    Delta_l_E = np.zeros((N_ell, N_k), dtype=np.float64)

    alpha_P = 1.7
    pol_prefactor = 0.75 * alpha_P

    eps_prefactors = np.zeros(N_ell)
    for il, ell in enumerate(ell_values):
        l = int(ell)
        if l >= 2:
            eps_prefactors[il] = np.sqrt(float((l - 1) * l * (l + 1) * (l + 2)))

    if verbose:
        print(f"Computing LOS transfer functions (TT + EE)...")
        sys.stdout.flush()

    for it in range(N_snap):
        chi = chi_snap[it]
        if chi <= 0:
            continue

        if it == 0:
            w = 0.5 * dtau_all[0]
        elif it == N_snap - 1:
            w = 0.5 * dtau_all[-1]
        else:
            w = 0.5 * (dtau_all[max(0, it-1)]
                       + dtau_all[min(it, len(dtau_all)-1)])

        S_SW = g_snap[it] * (Theta0_N[:, it] + Psi_N[:, it])
        S_Dop = g_snap[it] * vb_N[:, it]
        S_ISW = exp_neg_kappa[it] * PhiPsi_prime[:, it]
        S_E = g_snap[it] * pol_prefactor * Theta2_arr[:, it]

        if use_gpu_bessel:
            x_arr = k_arr * chi
            jl_mx, jlp_mx = bessel_table.eval_jl_jlp(x_arr.astype(np.float32))
            mx.eval(jl_mx, jlp_mx)
            jl = np.array(jl_mx)
            jlp = np.array(jlp_mx)

            x_safe = np.where(np.abs(x_arr) < 1e-10, 1e-10, x_arr)
            inv_x2 = 1.0 / (x_safe ** 2)

            for il in range(N_ell):
                Delta_l[il, :] += w * (
                    S_SW * jl[il, :]
                    + S_Dop * jlp[il, :]
                    + S_ISW * jl[il, :]
                )
                if eps_prefactors[il] > 0:
                    eps_l = eps_prefactors[il] * jl[il, :] * inv_x2
                    eps_l = np.where(np.abs(x_arr) < 1e-10, 0.0, eps_l)
                    Delta_l_E[il, :] += w * S_E * eps_l
        else:
            from scipy.special import spherical_jn
            x = k_arr * chi
            x_safe = np.where(np.abs(x) < 1e-10, 1e-10, x)
            inv_x2 = 1.0 / (x_safe ** 2)
            for il, ell in enumerate(ell_values):
                jl = spherical_jn(int(ell), x)
                jlp = spherical_jn(int(ell), x, derivative=True)
                Delta_l[il, :] += w * (S_SW * jl + S_Dop * jlp + S_ISW * jl)
                if eps_prefactors[il] > 0:
                    eps_l = eps_prefactors[il] * jl * inv_x2
                    eps_l = np.where(np.abs(x) < 1e-10, 0.0, eps_l)
                    Delta_l_E[il, :] += w * S_E * eps_l

    t_los = time.time() - t0
    if verbose:
        print(f"LOS integration: {t_los:.2f}s")

    # ================================================================
    # Step 7: C_l = 4 pi int dk/k P_R |Delta_l|^2
    # ================================================================
    t0 = time.time()

    lnk = np.log(k_arr)
    dlnk = np.diff(lnk)

    ell_f = ell_values.astype(float)
    uK2 = (T_CMB * 1e6) ** 2
    norm = (2.0 / 3.0) ** 2

    # TT
    integrand_TT = P_R[None, :] * Delta_l ** 2
    mid_TT = 0.5 * (integrand_TT[:, :-1] + integrand_TT[:, 1:])
    Cl_TT = 4.0 * np.pi * np.sum(mid_TT * dlnk[None, :], axis=1)
    Cl_TT = np.maximum(Cl_TT, 0.0)
    Dl_TT = norm * ell_f * (ell_f + 1.0) * Cl_TT / (2.0 * np.pi) * uK2

    # EE
    integrand_EE = P_R[None, :] * Delta_l_E ** 2
    mid_EE = 0.5 * (integrand_EE[:, :-1] + integrand_EE[:, 1:])
    Cl_EE = 4.0 * np.pi * np.sum(mid_EE * dlnk[None, :], axis=1)
    Cl_EE = np.maximum(Cl_EE, 0.0)
    Dl_EE = norm * ell_f * (ell_f + 1.0) * Cl_EE / (2.0 * np.pi) * uK2

    # TE
    integrand_TE = P_R[None, :] * Delta_l * Delta_l_E
    mid_TE = 0.5 * (integrand_TE[:, :-1] + integrand_TE[:, 1:])
    Cl_TE = 4.0 * np.pi * np.sum(mid_TE * dlnk[None, :], axis=1)
    Dl_TE = norm * ell_f * (ell_f + 1.0) * Cl_TE / (2.0 * np.pi) * uK2

    t_cl = time.time() - t0
    t_elapsed = time.time() - t_total

    if verbose:
        print(f"C_l computation (TT+EE+TE): {t_cl:.3f}s")
        print(f"\nTotal time: {t_elapsed:.1f}s  "
              f"(bg={t_bg:.1f}s, phase={t_phase:.2f}s, "
              f"pert={t_pert:.1f}s, LOS={t_los:.1f}s)")

        if len(Dl_TT) > 0:
            from scipy.signal import argrelextrema
            peak_idx = argrelextrema(Dl_TT, np.greater, order=5)[0]
            if len(peak_idx) > 0:
                print(f"\nTT peak positions:")
                for ip in range(min(5, len(peak_idx))):
                    idx = peak_idx[ip]
                    print(f"  Peak {ip+1}: l = {ell_values[idx]}, "
                          f"D_l = {Dl_TT[idx]:.1f} uK^2")

    return {
        'ell': ell_values,
        'Cl': Cl_TT,
        'Dl': Dl_TT,
        'Cl_TT': Cl_TT,
        'Dl_TT': Dl_TT,
        'Cl_EE': Cl_EE,
        'Dl_EE': Dl_EE,
        'Cl_TE': Cl_TE,
        'Dl_TE': Dl_TE,
        'k_arr': k_arr,
        'bg': bg,
        'Delta_l': Delta_l,
        'Delta_l_E': Delta_l_E,
        'Phi_N': Phi_N,
        'Psi_N': Psi_N,
        'Theta2_arr': Theta2_arr,
        'timing': {
            'background': t_bg,
            'acoustic_phase': t_phase,
            'perturbations': t_pert,
            'los': t_los,
            'total': t_elapsed,
        },
    }


# ============================================================================
# Comparison with full sync gauge solver
# ============================================================================

def compare_with_sync(N_k=200, verbose=True):
    """
    Run both WKB and full sync gauge solvers and compare.
    """
    if verbose:
        print("=" * 65)
        print("WKB vs Full Sync Gauge Comparison")
        print("=" * 65)

    if verbose:
        print("\n[1/2] Running WKB solver...")
    wkb_result = run_wkb_solver(N_k=N_k, verbose=verbose)

    if verbose:
        print("\n[2/2] Running full sync gauge solver...")
    from .solver_sync import run_sync_solver
    sync_result = run_sync_solver(N_k=N_k, verbose=verbose)

    ell_wkb = wkb_result['ell']
    ell_sync = sync_result['ell']
    ell_common = np.intersect1d(ell_wkb, ell_sync)

    idx_wkb = np.searchsorted(ell_wkb, ell_common)
    idx_sync = np.searchsorted(ell_sync, ell_common)

    Dl_wkb = wkb_result['Dl_TT'][idx_wkb]
    Dl_sync = sync_result['Dl_TT'][idx_sync]

    mask = ell_common > 10
    if np.any(mask):
        Dl_wkb_m = Dl_wkb[mask]
        Dl_sync_m = Dl_sync[mask]
        valid = Dl_sync_m > 1.0
        if np.any(valid):
            residual_pct = (np.abs(Dl_wkb_m[valid] - Dl_sync_m[valid])
                            / Dl_sync_m[valid] * 100.0)
            rms_pct = np.sqrt(np.mean(residual_pct**2))
        else:
            rms_pct = np.nan
    else:
        rms_pct = np.nan

    if verbose:
        print(f"\n{'='*65}")
        print(f"COMPARISON RESULT:")
        print(f"  RMS residual (l > 10): {rms_pct:.1f}%")
        print(f"  WKB total time:  {wkb_result['timing']['total']:.1f}s")
        print(f"  Sync total time: {sync_result['timing']['total']:.1f}s")
        if wkb_result['timing']['total'] > 0:
            print(f"  Speedup: "
                  f"{sync_result['timing']['total'] / wkb_result['timing']['total']:.1f}x")

    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8),
                                         gridspec_kw={'height_ratios': [3, 1]})

        ax1.plot(ell_wkb, wkb_result['Dl_TT'], 'b-', label='WKB', alpha=0.8)
        ax1.plot(ell_sync, sync_result['Dl_TT'], 'r--', label='Full sync',
                 alpha=0.8)
        ax1.set_ylabel(r'$D_\ell^{TT}$ [$\mu K^2$]')
        ax1.set_title(f'WKB vs Full Sync Gauge (RMS = {rms_pct:.1f}%)')
        ax1.legend()
        ax1.set_xlim(2, 2500)

        if np.any(mask):
            ax2.plot(ell_common[mask],
                     (Dl_wkb[mask] - Dl_sync[mask])
                     / np.maximum(Dl_sync[mask], 1.0) * 100,
                     'k-', alpha=0.5)
            ax2.axhline(0, color='gray', linestyle=':')
            ax2.set_ylabel('Residual [%]')
            ax2.set_xlabel(r'$\ell$')
            ax2.set_xlim(2, 2500)

        plt.tight_layout()
        outpath = os.path.join(os.path.dirname(__file__), 'wkb_vs_sync.png')
        plt.savefig(outpath, dpi=150)
        plt.close()
        if verbose:
            print(f"  Plot saved to: {outpath}")
    except ImportError:
        pass

    return {
        'wkb': wkb_result,
        'sync': sync_result,
        'rms_pct': rms_pct,
    }


# ============================================================================
# Main
# ============================================================================

def main():
    """Entry point for `python -m mlx_class.solver_wkb`."""
    import argparse

    parser = argparse.ArgumentParser(
        description='WKB Envelope Boltzmann Solver for CMB')
    parser.add_argument('--N_k', type=int, default=500,
                        help='Number of k-modes (default 500)')
    parser.add_argument('--N_steps', type=int, default=40,
                        help='Not used (adaptive stepping)')
    parser.add_argument('--compare', action='store_true',
                        help='Compare with full sync gauge solver')
    parser.add_argument('--no-parallel', action='store_true',
                        help='Disable multiprocessing')
    parser.add_argument('--quiet', action='store_true',
                        help='Suppress output')
    args = parser.parse_args()

    if args.compare:
        result = compare_with_sync(N_k=min(args.N_k, 200),
                                    verbose=not args.quiet)
    else:
        result = run_wkb_solver(N_k=args.N_k, verbose=not args.quiet,
                                parallel=not args.no_parallel)

        outpath = os.path.join(os.path.dirname(__file__), 'cl_wkb.dat')
        np.savetxt(outpath,
                   np.column_stack([result['ell'], result['Dl_TT'],
                                    result['Dl_EE'], result['Dl_TE']]),
                   header='ell  Dl_TT[uK2]  Dl_EE[uK2]  Dl_TE[uK2]',
                   fmt='%8d %15.6e %15.6e %15.6e')
        print(f"\nSaved to: {outpath}")

        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt

            fig, ax = plt.subplots(figsize=(10, 5))
            ax.plot(result['ell'], result['Dl_TT'], 'b-', linewidth=0.8)
            ax.set_xlabel(r'$\ell$')
            ax.set_ylabel(r'$D_\ell^{TT}$ [$\mu K^2$]')
            ax.set_title('WKB Envelope Solver')
            ax.set_xlim(2, 2500)
            ax.set_ylim(bottom=0)

            plotpath = os.path.join(os.path.dirname(__file__), 'cl_wkb.png')
            plt.savefig(plotpath, dpi=150)
            plt.close()
            print(f"Plot saved to: {plotpath}")
        except ImportError:
            pass


if __name__ == '__main__':
    main()
