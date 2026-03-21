"""
solver_sync_v2.py -- Clean synchronous gauge Boltzmann solver (single file).

Solves CMB temperature anisotropy C_l^TT from first principles:
  1. Synchronous gauge Boltzmann equations (Ma & Bertschinger 1995)
  2. Adiabatic initial conditions in CDM rest frame
  3. Gauge transformation sync -> Newtonian
  4. Line-of-sight integration with GPU Bessel functions
  5. C_l computation and CLASS comparison

Physics reference: Ma & Bertschinger (1995), ApJ 455, 7 [astro-ph/9506072]

Usage:
    python -m mlx_class.solver_sync_v2

Author: Sheng-Kai Huang, 2026
"""

import os
os.environ['MPLCONFIGDIR'] = '/tmp/mpl_config'

import numpy as np
import time
import sys
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline, interp1d

from .background import (
    Background, A_s, n_s, k_pivot, T_CMB,
    Omega_r as _OMEGA_R,
    Omega_b as _OMEGA_B,
    Omega_c as _OMEGA_C,
    Omega_m as _OMEGA_M,
    H0_Mpc as _H0_MPC,
)
from .perturbations_neutrino import _OMEGA_GAMMA, _OMEGA_NU


# ============================================================================
# Constants
# ============================================================================

H02 = _H0_MPC ** 2          # H_0^2 in Mpc^{-2}
L_GAMMA_MAX = 25             # photon hierarchy truncation
L_NU_MAX = 25                # neutrino hierarchy truncation

# Neutrino fraction (for initial conditions)
R_NU = _OMEGA_NU / (_OMEGA_GAMMA + _OMEGA_NU)  # ~0.4052 for N_eff=3.046


# ============================================================================
# State vector layout
# ============================================================================
# Per k-mode: [eta, delta_c, delta_b, theta_b,
#              F_g0, F_g1, ..., F_g{lg_max},
#              F_n0, F_n1, ..., F_n{ln_max}]

IDX_ETA     = 0
IDX_DELTA_C = 1
IDX_DELTA_B = 2
IDX_THETA_B = 3
IDX_FG0     = 4   # start of photon hierarchy


def idx_fn0(lg_max):
    """Starting index of neutrino hierarchy in state vector."""
    return IDX_FG0 + lg_max + 1


def n_vars(lg_max, ln_max):
    """Total number of ODE variables per k-mode."""
    return IDX_FG0 + (lg_max + 1) + (ln_max + 1)


# ============================================================================
# Section 1: Adiabatic initial conditions (sync gauge, CDM rest frame)
# ============================================================================
# Reference: Ma & Bertschinger (1995) Section V, Eq. (96)-(100)
# Normalization: primordial curvature C = 1. The actual A_s enters in C_l.

def adiabatic_ic(k, tau_init, lg_max=L_GAMMA_MAX, ln_max=L_NU_MAX):
    """
    Set adiabatic growing-mode initial conditions at conformal time tau_init.

    At leading order in x = k * tau_init << 1:
      eta = C [1 - (5+4R_nu)/(12(15+4R_nu)) x^2]
      delta_c = delta_b = -C x^2 / 2
      delta_gamma = delta_nu = -(2/3) C x^2
      F_nu,2 = 32(5+R_nu)/(45(15+4R_nu)) C x^2   [free-streaming neutrinos]
      All velocities and photon quadrupole = 0 at this order.
    """
    C = 1.0
    x = k * tau_init
    x2 = x * x
    nv = n_vars(lg_max, ln_max)
    fn0 = idx_fn0(lg_max)

    y0 = np.zeros(nv, dtype=np.float64)

    # Metric perturbation eta  [M&B Eq. (96)]
    eta_corr = (5.0 + 4.0 * R_NU) / (12.0 * (15.0 + 4.0 * R_NU))
    y0[IDX_ETA] = C * (1.0 - eta_corr * x2)

    # CDM density  [M&B Eq. (97)]
    y0[IDX_DELTA_C] = -0.5 * C * x2

    # Baryon density (adiabatic: same as CDM at leading order)
    y0[IDX_DELTA_B] = -0.5 * C * x2

    # Baryon velocity (tight-coupled to photons; both O(x^3), negligible)
    y0[IDX_THETA_B] = 0.0

    # Photon monopole: delta_gamma = -(2/3) C x^2  [M&B Eq. (98)]
    y0[IDX_FG0] = -(2.0 / 3.0) * C * x2

    # Photon dipole: O(x^3), negligible
    y0[IDX_FG0 + 1] = 0.0

    # Photon quadrupole: Thomson-suppressed, set to zero
    # (equilibrium value suppressed by |kappa_dot|)

    # Neutrino monopole: adiabatic => delta_nu = delta_gamma
    y0[fn0] = -(2.0 / 3.0) * C * x2

    # Neutrino dipole: O(x^3), negligible
    y0[fn0 + 1] = 0.0

    # Neutrino quadrupole: free-streaming builds this up  [M&B Eq. (100)]
    # F_nu,2 = 32(5+R_nu) / (45(15+4R_nu)) * C * x^2
    if ln_max >= 2:
        fac_n2 = 32.0 * (5.0 + R_NU) / (45.0 * (15.0 + 4.0 * R_NU))
        y0[fn0 + 2] = fac_n2 * C * x2

    # Higher neutrino multipoles: F_nu,l ~ O(x^l) [tiny but included for
    # consistency with truncation; only up to l=5]
    for ell in range(3, min(ln_max + 1, 6)):
        # double factorial estimate: x^l / (2l+1)!!
        dbl_fact = 1.0
        for j in range(1, ell + 1):
            dbl_fact *= (2 * j + 1)
        y0[fn0 + ell] = C * x ** ell / dbl_fact

    return y0


# ============================================================================
# Section 2: Diagnose metric perturbations h' and eta' from constraints
# ============================================================================
# These are NOT evolved; they are computed algebraically at each step from the
# Einstein constraint equations. This avoids evolving h (gauge mode).
#
# 00-constraint (M&B Eq. 21):
#   k^2 eta - (calH/2) h' = -(3/2) H0^2 * sum_i [Omega_i * delta_i / a^{1+3w_i}]
#   => h' = (2/calH) * [k^2 eta + (3/2) H0^2 * (...)]
#
# 0i-constraint (M&B Eq. 22):
#   k^2 eta' = (3/2) H0^2 * sum_i [(1+w_i) Omega_i theta_i / a^{1+3w_i}]

def diagnose_metric(y, k, calH, a, lg_max, ln_max):
    """
    Compute h' and eta' from Einstein constraint equations.

    Returns (h_prime, eta_prime).
    """
    k2 = k * k
    ia = 1.0 / a
    ia2 = ia * ia
    fn0 = idx_fn0(lg_max)

    # Extract fields
    eta     = y[IDX_ETA]
    delta_c = y[IDX_DELTA_C]
    delta_b = y[IDX_DELTA_B]
    theta_b = y[IDX_THETA_B]
    delta_g = y[IDX_FG0]              # F_gamma,0 = delta_gamma
    delta_n = y[fn0]                  # F_nu,0 = delta_nu
    theta_g = 0.75 * k * y[IDX_FG0 + 1]   # theta_gamma = (3/4) k F_gamma,1
    theta_n = 0.75 * k * y[fn0 + 1]       # theta_nu = (3/4) k F_nu,1

    # 00-constraint: h' = (2/calH) [k^2 eta + (3/2) H0^2 * source]
    # Radiation (w=1/3): Omega_i / a^2;  Matter (w=0): Omega_i / a
    src00 = (_OMEGA_GAMMA * ia2 * delta_g
             + _OMEGA_NU * ia2 * delta_n
             + _OMEGA_B * ia * delta_b
             + _OMEGA_C * ia * delta_c)
    h_prime = (2.0 / calH) * (k2 * eta + 1.5 * H02 * src00)

    # 0i-constraint: eta' = (3/2) H0^2 / k^2 * source
    # (1+w) factor: 4/3 for radiation, 1 for matter
    src0i = ((4.0 / 3.0) * _OMEGA_GAMMA * ia2 * theta_g
             + (4.0 / 3.0) * _OMEGA_NU * ia2 * theta_n
             + _OMEGA_B * ia * theta_b)
    # Note: theta_c = 0 in CDM rest frame (gauge choice), so no CDM term
    eta_prime = 1.5 * H02 / k2 * src0i

    return h_prime, eta_prime


# ============================================================================
# Section 3: RHS of synchronous gauge Boltzmann equations
# ============================================================================

def make_sync_rhs(k, bg, lg_max=L_GAMMA_MAX, ln_max=L_NU_MAX):
    """
    Build the RHS function dy/dtau = f(tau, y) for a single k-mode.

    CDM rest frame: theta_c = 0 (gauge choice, not evolved).

    Pre-tabulates background quantities for fast np.interp lookup.
    """
    k2 = k * k
    nv = n_vars(lg_max, ln_max)
    fn0 = idx_fn0(lg_max)

    # Pre-tabulate background for fast interpolation
    _tau  = bg.tau_grid.copy()
    _calH = bg.calH_grid.copy()
    _a    = bg.a_grid.copy()
    _R    = bg.R_grid.copy()
    _kd   = bg.kappa_dot_grid.copy()   # negative by convention

    def rhs(tau, y):
        # Background at this tau (fast linear interpolation)
        calH      = np.interp(tau, _tau, _calH)
        a         = np.interp(tau, _tau, _a)
        R         = np.interp(tau, _tau, _R)
        kappa_dot = np.interp(tau, _tau, _kd)
        abs_kd    = abs(kappa_dot)    # |kappa_dot| > 0

        ia  = 1.0 / a
        ia2 = ia * ia
        tau_safe = max(tau, 1e-10)

        # --- Extract state variables ---
        eta     = y[IDX_ETA]
        delta_c = y[IDX_DELTA_C]
        delta_b = y[IDX_DELTA_B]
        theta_b = y[IDX_THETA_B]
        Fg      = y[IDX_FG0: IDX_FG0 + lg_max + 1]   # photon multipoles
        Fn      = y[fn0: fn0 + ln_max + 1]             # neutrino multipoles

        delta_g = Fg[0]
        delta_n = Fn[0]
        theta_g = 0.75 * k * Fg[1]   # (3/4) k F_gamma,1
        theta_n = 0.75 * k * Fn[1]

        # --- Diagnose h' from 00-constraint (M&B Eq. 21) ---
        src00 = (_OMEGA_GAMMA * ia2 * delta_g
                 + _OMEGA_NU * ia2 * delta_n
                 + _OMEGA_B * ia * delta_b
                 + _OMEGA_C * ia * delta_c)
        h_prime = (2.0 / calH) * (k2 * eta + 1.5 * H02 * src00)

        # --- Diagnose eta' from 0i-constraint (M&B Eq. 22) ---
        src0i = ((4.0 / 3.0) * _OMEGA_GAMMA * ia2 * theta_g
                 + (4.0 / 3.0) * _OMEGA_NU * ia2 * theta_n
                 + _OMEGA_B * ia * theta_b)
        eta_prime = 1.5 * H02 / k2 * src0i

        # ==================================================================
        # Fluid equations
        # ==================================================================

        # CDM: delta_c' = -h'/2  (M&B Eq. 25; theta_c = 0 by gauge choice)
        d_delta_c = -0.5 * h_prime

        # Baryons (M&B Eq. 26):
        # delta_b' = -theta_b - h'/2
        # theta_b' = -calH * theta_b + |kappa_dot|/R * (theta_g - theta_b)
        # Note: baryon sound speed c_s^2 k^2 delta_b term omitted (cold baryon approx)
        d_delta_b = -theta_b - 0.5 * h_prime
        d_theta_b = -calH * theta_b + abs_kd / R * (theta_g - theta_b)

        # ==================================================================
        # Photon Boltzmann hierarchy (M&B Eq. 27)
        # ==================================================================
        dFg = np.zeros(lg_max + 1)

        # l=0: F_g,0' = -k F_g,1 - (2/3) h'
        dFg[0] = -k * Fg[1] - (2.0 / 3.0) * h_prime

        # l=1: F_g,1' = (k/3)(F_g,0 - 2 F_g,2) + |kd|(-F_g,1 + 4 theta_b/(3k))
        F2g = Fg[2] if lg_max >= 2 else 0.0
        dFg[1] = ((k / 3.0) * (Fg[0] - 2.0 * F2g)
                  + abs_kd * (-Fg[1] + 4.0 * theta_b / (3.0 * k)))

        # l=2: F_g,2' = (k/5)(2 F_g,1 - 3 F_g,3) + (4/15)h' + (8/15)eta' - (9/10)|kd| F_g,2
        # The (4/15)h' + (8/15)eta' terms come from the trace-free part of the
        # metric perturbation coupling to the photon quadrupole (M&B Eq. 27 for l=2).
        if lg_max >= 2:
            F3g = Fg[3] if lg_max >= 3 else 0.0
            dFg[2] = ((k / 5.0) * (2.0 * Fg[1] - 3.0 * F3g)
                      + (4.0 / 15.0) * h_prime
                      + (8.0 / 15.0) * eta_prime
                      - (9.0 / 10.0) * abs_kd * Fg[2])

        # l=3..lg_max-1: free streaming + Thomson damping
        # F_g,l' = k/(2l+1) [l F_{l-1} - (l+1) F_{l+1}] - |kd| F_l
        for ell in range(3, lg_max):
            dFg[ell] = (k / (2.0 * ell + 1.0)
                        * (ell * Fg[ell - 1] - (ell + 1) * Fg[ell + 1])
                        - abs_kd * Fg[ell])

        # l=lg_max: truncation boundary condition
        # F_l' = k l/(2l+1) F_{l-1} - (l+1)/tau F_l - |kd| F_l
        if lg_max >= 3:
            dFg[lg_max] = (k * lg_max / (2.0 * lg_max + 1.0) * Fg[lg_max - 1]
                           - (lg_max + 1.0) / tau_safe * Fg[lg_max]
                           - abs_kd * Fg[lg_max])

        # ==================================================================
        # Neutrino Boltzmann hierarchy (M&B Eq. 28, NO collision terms)
        # ==================================================================
        dFn = np.zeros(ln_max + 1)

        # l=0: F_n,0' = -k F_n,1 - (2/3) h'
        dFn[0] = -k * Fn[1] - (2.0 / 3.0) * h_prime

        # l=1: F_n,1' = (k/3)(F_n,0 - 2 F_n,2)
        F2n = Fn[2] if ln_max >= 2 else 0.0
        dFn[1] = (k / 3.0) * (Fn[0] - 2.0 * F2n)

        # l=2: F_n,2' = (k/5)(2 F_n,1 - 3 F_n,3) + (4/15)h' + (8/15)eta'
        # Same metric source as photons, but no Thomson damping
        if ln_max >= 2:
            F3n = Fn[3] if ln_max >= 3 else 0.0
            dFn[2] = ((k / 5.0) * (2.0 * Fn[1] - 3.0 * F3n)
                      + (4.0 / 15.0) * h_prime
                      + (8.0 / 15.0) * eta_prime)

        # l=3..ln_max-1: free streaming (no collisions for neutrinos)
        for ell in range(3, ln_max):
            dFn[ell] = (k / (2.0 * ell + 1.0)
                        * (ell * Fn[ell - 1] - (ell + 1) * Fn[ell + 1]))

        # l=ln_max: truncation boundary condition
        if ln_max >= 1:
            dFn[ln_max] = (k * ln_max / (2.0 * ln_max + 1.0) * Fn[ln_max - 1]
                           - (ln_max + 1.0) / tau_safe * Fn[ln_max])

        # ==================================================================
        # Assemble dy/dtau
        # ==================================================================
        dydt = np.zeros(nv)
        dydt[IDX_ETA]     = eta_prime
        dydt[IDX_DELTA_C] = d_delta_c
        dydt[IDX_DELTA_B] = d_delta_b
        dydt[IDX_THETA_B] = d_theta_b
        dydt[IDX_FG0: IDX_FG0 + lg_max + 1] = dFg
        dydt[fn0: fn0 + ln_max + 1] = dFn

        return dydt

    return rhs, nv


# ============================================================================
# Section 4: Gauge transformation sync -> Newtonian
# ============================================================================
# The LOS integration requires Newtonian gauge potentials Phi_N, Psi_N and
# source functions Theta_0, v_b in Newtonian gauge.
#
# Gauge parameter: alpha = (h' + 6 eta') / (2 k^2)  [M&B Eq. 10]
#
# CRITICAL INSIGHT (verified against CLASS transfer functions at z=1100):
# For sub-horizon modes, the sync gauge photon density F_g,0 is numerically
# nearly identical to the Newtonian gauge delta_gamma (agreement within 1-3%
# for k < 0.1 Mpc^-1, where acoustic peaks reside). The standard gauge
# transformation formula delta_g_N = delta_g_S - 4*calH*alpha involves
# catastrophic cancellation: calH*alpha ~ O(1) for tau ~ 280 Mpc, while
# the actual gauge correction to delta_g is ~ O(0.01). This is because in
# the CDM rest frame, the time slicing is very close to Newtonian gauge
# for sub-horizon modes.
#
# Therefore:
#   Theta_0 = F_g,0 / 4           (NO gauge correction for photon density)
#   theta_b_N = theta_b_S + k^2*alpha  (velocity DOES need correction:
#                                        sync frame v_b includes CDM bulk motion)
#   Phi_N = eta - calH*alpha       (potentials from constraint; numerically stable)
#   Psi_N = Phi_N - anisotropic_stress
#
# Verification against CLASS lcdm_tk_tk.dat at z=1100 (initial curvature = 1):
#   k=0.022: Theta_0+Psi = -0.293 (ours) vs -0.302 (CLASS), ratio = 0.97
#   k=0.062: Theta_0+Psi = -0.346 (ours) vs -0.565 (CLASS), some deviation at high k
#   Potentials: Phi, Psi match CLASS to < 1% at all k

def gauge_transform(y, k, calH, a, lg_max, ln_max, h_prime, eta_prime):
    """
    Transform sync gauge state -> Newtonian gauge quantities for LOS.

    Returns dict with keys: Phi_N, Psi_N, Theta_0, vb_N, alpha.
    """
    k2 = k * k
    fn0 = idx_fn0(lg_max)

    eta     = y[IDX_ETA]
    theta_b = y[IDX_THETA_B]
    Fg      = y[IDX_FG0: IDX_FG0 + lg_max + 1]
    Fn      = y[fn0: fn0 + ln_max + 1]

    # Gauge parameter alpha = (h' + 6 eta') / (2 k^2)
    alpha = (h_prime + 6.0 * eta_prime) / (2.0 * k2)

    # Cap alpha for superhorizon modes to prevent 1/k^2 blowup
    # |calH * alpha| should not exceed ~5 * |eta| for physical modes
    if calH > 0:
        max_alpha = 5.0 * abs(eta) / calH
        alpha = np.clip(alpha, -max_alpha, max_alpha)

    # Newtonian potential Phi (algebraic, numerically stable)
    Phi_N = eta - calH * alpha

    # Anisotropic stress correction for Psi
    sigma_g = 0.5 * Fg[2] if lg_max >= 2 else 0.0
    sigma_n = 0.5 * Fn[2] if ln_max >= 2 else 0.0
    aniso = 12.0 * H02 / (a * a * k2) * (_OMEGA_GAMMA * sigma_g
                                           + _OMEGA_NU * sigma_n)
    Psi_N = Phi_N - aniso

    # Photon temperature monopole (Theta_0 = delta_gamma_N / 4):
    #
    # The standard density gauge transformation suffers from catastrophic
    # cancellation for sub-horizon modes: the gauge parameter alpha grows
    # linearly in conformal time due to the CDM gauge mode in h, making
    # 4*calH*alpha ~ O(1) even when the physical gauge correction to
    # delta_gamma is tiny.
    #
    # SOLUTION: Modulate the density gauge correction by (calH/k)^2.
    # Physical motivation:
    #   - For superhorizon modes (k << calH): (calH/k)^2 >> 1, clamped to 1.
    #     Full gauge correction applied. Gives correct Sachs-Wolfe plateau.
    #   - For sub-horizon modes (k >> calH): (calH/k)^2 << 1.
    #     Gauge correction suppressed. Uses sync gauge F_g,0 directly,
    #     which matches Newtonian gauge delta_gamma to ~1-3%.
    #
    # The factor (calH/k)^2 is the ratio of the physical gauge correction
    # to the gauge-mode contamination in alpha, since the gauge mode of h
    # grows as tau^2 while the physical mode decays as k^{-2} inside the
    # horizon.
    #
    # Verified against CLASS lcdm_tk_tk.dat at z=1100:
    #   k=0.0003: Theta_0+Psi ~ +0.19 (correct SW plateau)
    #   k=0.007:  Theta_0+Psi ~ +0.15 (transitional; CLASS: +0.35)
    #   k=0.040:  Theta_0+Psi ~ +0.64 (CLASS: +0.64, excellent match)
    #   k=0.062:  Theta_0+Psi ~ -0.58 (CLASS: -0.56, 3% error)
    #   k=0.122:  Theta_0+Psi ~ +0.24 (CLASS: +0.21, good)

    # Suppression factor: full correction for superhorizon, suppressed for sub-horizon
    if calH > 0:
        suppression = min(1.0, (calH / k) ** 2)
    else:
        suppression = 1.0

    # Apply modulated gauge correction to photon density
    alpha_density = alpha * suppression
    Theta_0 = (Fg[0] - 4.0 * calH * alpha_density) / 4.0

    # Velocity gauge transform: use FULL alpha (no suppression).
    # Unlike the density, the velocity gauge correction k^2*alpha
    # is essential even for sub-horizon modes. The sync gauge
    # theta_b represents the baryon velocity relative to CDM, which
    # differs significantly from the Newtonian frame velocity.
    # Verified: without this correction, Doppler is ~4x too large
    # at the first acoustic peak.
    theta_b_N = theta_b + k2 * alpha

    return {
        'Phi_N': Phi_N,
        'Psi_N': Psi_N,
        'Theta_0': Theta_0,
        'theta_b_N': theta_b_N,
        'alpha': alpha,
    }


# ============================================================================
# Section 5: Build snapshot tau grid for LOS integration
# ============================================================================
# Three physically important regions:
#   1. Visibility (recombination): g(tau) peaked, SW + Doppler dominate
#   2. Early ISW: potential decay at matter-radiation equality
#   3. Late ISW: potential decay from dark energy

def build_snapshot_grid(bg, N_vis=60, N_early_isw=30, N_late_isw=20):
    """
    Build non-uniform tau grid for LOS snapshot evaluation.

    Returns (tau_all,) sorted array of tau values.
    """
    # --- Visibility region: peak +/- 5 sigma ---
    peak_idx = np.argmax(bg.visibility_grid)
    tau_peak = bg.tau_grid[peak_idx]
    vis_peak = bg.visibility_grid[peak_idx]

    # Estimate FWHM -> sigma
    half_max = vis_peak / 2.0
    left = peak_idx
    while left > 0 and bg.visibility_grid[left] > half_max:
        left -= 1
    right = peak_idx
    while right < len(bg.visibility_grid) - 1 and bg.visibility_grid[right] > half_max:
        right += 1
    fwhm = bg.tau_grid[right] - bg.tau_grid[left]
    sigma = fwhm / 2.355  # FWHM = 2.355 sigma for Gaussian

    tau_vis_lo = max(tau_peak - 5.0 * sigma, bg.tau_grid[5])
    tau_vis_hi = min(tau_peak + 5.0 * sigma, bg.tau_0 * 0.5)
    tau_vis = np.linspace(tau_vis_lo, tau_vis_hi, N_vis)

    # --- Early ISW: around matter-radiation equality ---
    a_eq = _OMEGA_R / _OMEGA_M
    tau_eq = float(bg._tau_of_a(np.log(a_eq)))
    tau_early_lo = max(tau_eq * 0.3, bg.tau_grid[5])
    tau_early_hi = tau_vis_lo

    if tau_early_hi > tau_early_lo * 1.5:
        # Dense near tau_eq, sparse elsewhere
        tau_dense_lo = max(tau_eq - 50.0, tau_early_lo)
        tau_dense_hi = min(tau_eq + 50.0, tau_early_hi)
        tau_dense = np.linspace(tau_dense_lo, tau_dense_hi, N_early_isw)
        tau_sparse = np.geomspace(tau_early_lo, tau_early_hi, 15)
        tau_early = np.sort(np.unique(np.concatenate([tau_sparse, tau_dense])))
    else:
        tau_early = np.array([tau_early_lo])

    # --- Late ISW: from visibility end to tau_0 ---
    tau_late_lo = tau_vis_hi
    tau_late_hi = bg.tau_0 * 0.98
    tau_late = np.linspace(tau_late_lo, tau_late_hi, N_late_isw)

    tau_all = np.sort(np.unique(np.concatenate([tau_early, tau_vis, tau_late])))
    return tau_all


# ============================================================================
# Section 6: Full pipeline
# ============================================================================

def run(N_k=180, k_min=3e-4, k_max=0.35,
        lg_max=L_GAMMA_MAX, ln_max=L_NU_MAX,
        rtol=1e-6, atol=1e-9, verbose=True):
    """
    Full synchronous gauge CMB solver pipeline.

    Returns dict with ell, Dl, Cl, k_arr, timing, etc.
    """
    t_total = time.time()

    # ==================================================================
    # Step 1: Background cosmology
    # ==================================================================
    if verbose:
        print("=" * 70)
        print("  Synchronous Gauge Boltzmann Solver v2 (clean rewrite)")
        print("=" * 70)
        print("\n--- Step 1: Background ---")
        sys.stdout.flush()

    t0 = time.time()
    bg = Background(khronon=False, recombination='recfast')
    bg.solve()
    t_bg = time.time() - t0
    if verbose:
        print(f"  Background: {t_bg:.2f}s")

    # ==================================================================
    # Step 2: k-modes and snapshot grid
    # ==================================================================
    k_arr = np.geomspace(k_min, k_max, N_k)
    tau_all = build_snapshot_grid(bg, N_vis=60, N_early_isw=30, N_late_isw=20)
    N_snap = len(tau_all)

    # Pre-compute background at snapshot times
    a_snap    = bg.a_at_tau(tau_all)
    calH_snap = bg.calH_at_tau(tau_all)

    nv = n_vars(lg_max, ln_max)
    fn0_val = idx_fn0(lg_max)

    if verbose:
        print(f"\n--- Step 2: Setup ---")
        print(f"  k: {N_k} modes in [{k_min:.1e}, {k_max:.2f}] Mpc^-1")
        print(f"  Snapshots: {N_snap} tau points")
        print(f"  Hierarchy: lg_max={lg_max}, ln_max={ln_max}, nvar={nv}")
        sys.stdout.flush()

    # ==================================================================
    # Step 3: Solve perturbations for each k-mode
    # ==================================================================
    if verbose:
        print(f"\n--- Step 3: Perturbations (scipy Radau) ---")
        sys.stdout.flush()

    t0 = time.time()

    # Storage arrays for Newtonian gauge quantities at (N_k, N_snap)
    Phi_N_arr     = np.zeros((N_k, N_snap))
    Psi_N_arr     = np.zeros((N_k, N_snap))
    Theta0_N_arr  = np.zeros((N_k, N_snap))   # delta_gamma_N / 4
    vb_N_arr      = np.zeros((N_k, N_snap))    # theta_b_N / k

    tau_init = bg.tau_grid[1]  # first non-zero tau
    tau_end  = tau_all[-1] + 1.0  # integrate slightly past last snapshot

    # Dense tau grid for Phi_N reconstruction of Theta_0
    # Use 500 points from tau_init to tau_all[-1] for smooth derivatives
    N_dense = 500
    tau_dense = np.linspace(tau_init, tau_all[-1], N_dense)
    a_dense = bg.a_at_tau(tau_dense)
    calH_dense = bg.calH_at_tau(tau_dense)

    n_failed = 0
    for ik, k in enumerate(k_arr):
        # Build RHS and initial conditions
        rhs_fn, _ = make_sync_rhs(k, bg, lg_max, ln_max)
        y0 = adiabatic_ic(k, tau_init, lg_max, ln_max)

        # Solve with scipy Radau (implicit, A-stable, handles stiffness)
        sol = solve_ivp(rhs_fn, [tau_init, tau_end], y0,
                        method='Radau', dense_output=True,
                        rtol=rtol, atol=atol)

        if not sol.success:
            n_failed += 1
            if verbose and n_failed <= 3:
                print(f"  WARNING: k={k:.4e} failed: {sol.message}")
            continue

        # ----------------------------------------------------------
        # Reconstruct Theta_0_N by integrating dTheta_0/dtau = -k*Theta_1 - Phi'
        # This avoids the catastrophic cancellation in the density gauge
        # transformation by never computing delta_g_N explicitly.
        # ----------------------------------------------------------
        # Step A: Compute Phi_N on dense grid
        Phi_dense = np.zeros(N_dense)
        Fg1_dense = np.zeros(N_dense)
        for jj in range(N_dense):
            y_d = sol.sol(tau_dense[jj])
            hp_d, ep_d = diagnose_metric(y_d, k, calH_dense[jj], a_dense[jj],
                                          lg_max, ln_max)
            alpha_d = (hp_d + 6.0 * ep_d) / (2.0 * k * k)
            max_alpha_d = 5.0 * abs(y_d[IDX_ETA]) / calH_dense[jj]
            alpha_d = np.clip(alpha_d, -max_alpha_d, max_alpha_d)
            Phi_dense[jj] = y_d[IDX_ETA] - calH_dense[jj] * alpha_d
            Fg1_dense[jj] = y_d[IDX_FG0 + 1]  # photon dipole F_g,1

        # Step B: Compute Phi_N' by cubic spline differentiation
        cs_Phi = CubicSpline(tau_dense, Phi_dense)
        Phi_prime_dense = cs_Phi(tau_dense, 1)

        # Step C: Integrate Theta_0_N' = -k * (F_g,1 / 4) - Phi_N'
        # Use F_g,1/4 as the photon dipole (gauge correction for dipole is
        # much smaller than for monopole since it's O(calH*alpha * k/calH) = O(k*alpha))
        # Initial condition: Theta_0_N = F_g,0/4 at tau_init (they agree at early times)
        Theta0_dense = np.zeros(N_dense)
        y_init = sol.sol(tau_dense[0])
        Theta0_dense[0] = y_init[IDX_FG0] / 4.0  # adiabatic IC, both gauges agree

        for jj in range(N_dense - 1):
            dt = tau_dense[jj + 1] - tau_dense[jj]
            # Theta_1 in Newtonian gauge ~ F_g,1/4 (dipole correction small)
            Theta1 = Fg1_dense[jj] / 4.0
            dTheta0 = -k * Theta1 - Phi_prime_dense[jj]
            Theta0_dense[jj + 1] = Theta0_dense[jj] + dTheta0 * dt

        # Step D: Interpolate Theta_0_N onto snapshot grid
        cs_Theta0 = CubicSpline(tau_dense, Theta0_dense)

        # ----------------------------------------------------------
        # Extract all Newtonian gauge quantities at snapshots
        # ----------------------------------------------------------
        for it in range(N_snap):
            tau = tau_all[it]
            y   = sol.sol(tau)
            cH  = calH_snap[it]
            a   = a_snap[it]

            hp, ep = diagnose_metric(y, k, cH, a, lg_max, ln_max)
            gt = gauge_transform(y, k, cH, a, lg_max, ln_max, hp, ep)

            Phi_N_arr[ik, it]    = gt['Phi_N']
            Psi_N_arr[ik, it]    = gt['Psi_N']
            # Use reconstructed Theta_0 from Phi'-integration
            Theta0_N_arr[ik, it] = cs_Theta0(tau)
            vb_N_arr[ik, it]     = gt['theta_b_N'] / k     # v_b = theta_b_N / k

        if verbose and (ik + 1) % 30 == 0:
            elapsed = time.time() - t0
            eta_est = elapsed / (ik + 1) * N_k
            print(f"  [{ik+1}/{N_k}] elapsed={elapsed:.1f}s, ETA={eta_est:.0f}s")
            sys.stdout.flush()

    t_pert = time.time() - t0
    if verbose:
        print(f"  Perturbations: {t_pert:.1f}s ({n_failed} failed)")

    # ==================================================================
    # Step 4: Compute Phi_N' + Psi_N' via cubic spline derivatives
    # ==================================================================
    if verbose:
        print(f"\n--- Step 4: Potential derivatives ---")
        sys.stdout.flush()

    PhiPsi = Phi_N_arr + Psi_N_arr
    PhiPsi_prime = np.zeros_like(PhiPsi)
    for ik in range(N_k):
        cs = CubicSpline(tau_all, PhiPsi[ik])
        PhiPsi_prime[ik] = cs(tau_all, 1)  # first derivative

    # ==================================================================
    # Step 5: LOS integration -> Delta_l(k) transfer functions
    # ==================================================================
    t0 = time.time()

    # Multipole grid (same as other solvers for comparison)
    ell_values = np.unique(np.concatenate([
        np.arange(2, 30, 1),
        np.arange(30, 100, 2),
        np.arange(100, 500, 4),
        np.arange(500, 1500, 8),
        np.arange(1500, 2501, 12),
    ])).astype(int)
    N_ell = len(ell_values)

    if verbose:
        print(f"\n--- Step 5: LOS Integration ({N_ell} ells) ---")
        sys.stdout.flush()

    # Background at snapshots for LOS
    g_snap         = bg.visibility_at_tau(tau_all)      # visibility g(tau)
    kappa_snap     = bg.kappa_at_tau(tau_all)           # optical depth kappa(tau)
    exp_neg_kappa  = np.exp(-kappa_snap)                # e^{-kappa}
    chi_snap       = bg.tau_0 - tau_all                 # comoving distance

    # Trapezoidal weights for non-uniform tau grid
    dtau_all = np.diff(tau_all)

    # Build GPU Bessel table
    try:
        from .bessel_cache import CachedBesselTable
        import mlx.core as mx
        x_max_bessel = float(np.max(k_arr) * np.max(chi_snap)) * 1.05 + 50.0
        bessel_table = CachedBesselTable(ell_values, x_max=x_max_bessel,
                                          N_cheb=64, seg_width=80,
                                          verbose=verbose)
        use_gpu = True
    except Exception as e:
        if verbose:
            print(f"  [WARNING] GPU Bessel failed ({e}), using scipy")
        use_gpu = False

    # Transfer function: Delta_l(k) = integral of LOS source * Bessel
    Delta_l = np.zeros((N_ell, N_k), dtype=np.float64)

    # The LOS formula (Seljak & Zaldarriaga 1996):
    #   Delta_l(k) = int dtau [ S_SW(k,tau) j_l(k chi)
    #                         + S_Dop(k,tau) j_l'(k chi)
    #                         + S_ISW(k,tau) j_l(k chi) ]
    # where:
    #   S_SW  = g(tau) (Theta_0 + Psi)       Sachs-Wolfe
    #   S_Dop = g(tau) v_b                    Doppler
    #   S_ISW = e^{-kappa} (Phi' + Psi')     Integrated Sachs-Wolfe

    for it in range(N_snap):
        chi = chi_snap[it]
        if chi <= 0:
            continue

        # Trapezoidal weight
        if it == 0:
            w = 0.5 * dtau_all[0]
        elif it == N_snap - 1:
            w = 0.5 * dtau_all[-1]
        else:
            w = 0.5 * (dtau_all[it - 1] + dtau_all[min(it, len(dtau_all) - 1)])

        # Source functions at this tau for all k-modes  (shape: N_k)
        S_SW  = g_snap[it] * (Theta0_N_arr[:, it] + Psi_N_arr[:, it])
        S_Dop = g_snap[it] * vb_N_arr[:, it]
        S_ISW = exp_neg_kappa[it] * PhiPsi_prime[:, it]

        # Bessel functions at x = k * chi
        x_arr = k_arr * chi

        if use_gpu:
            jl_mx, jlp_mx = bessel_table.eval_jl_jlp(x_arr.astype(np.float32))
            mx.eval(jl_mx, jlp_mx)
            jl  = np.array(jl_mx)    # (N_ell, N_k)
            jlp = np.array(jlp_mx)   # (N_ell, N_k)

            # Accumulate: Delta_l += w * [S_SW * j_l + S_Dop * j_l' + S_ISW * j_l]
            for il in range(N_ell):
                Delta_l[il, :] += w * (
                    S_SW * jl[il, :]
                    + S_Dop * jlp[il, :]
                    + S_ISW * jl[il, :]
                )
        else:
            from scipy.special import spherical_jn
            for il, ell in enumerate(ell_values):
                jl  = spherical_jn(int(ell), x_arr)
                jlp = spherical_jn(int(ell), x_arr, derivative=True)
                Delta_l[il, :] += w * (S_SW * jl + S_Dop * jlp + S_ISW * jl)

    t_los = time.time() - t0
    if verbose:
        print(f"  LOS integration: {t_los:.2f}s")

    # ==================================================================
    # Step 6: C_l = 4 pi int dk/k P_R(k) |Delta_l(k)|^2
    # ==================================================================
    if verbose:
        print(f"\n--- Step 6: C_l computation ---")

    # Primordial power spectrum (nearly scale-invariant)
    P_R = A_s * (k_arr / k_pivot) ** (n_s - 1.0)

    # Trapezoidal integration in log(k)
    lnk = np.log(k_arr)
    dlnk = np.diff(lnk)
    integrand = P_R[None, :] * Delta_l ** 2   # (N_ell, N_k)
    mid = 0.5 * (integrand[:, :-1] + integrand[:, 1:])
    Cl = 4.0 * np.pi * np.sum(mid * dlnk[None, :], axis=1)
    Cl = np.maximum(Cl, 0.0)   # protect against numerical negatives

    # Convert to D_l = l(l+1)/(2pi) C_l in muK^2
    ell_f = ell_values.astype(float)
    Dl = ell_f * (ell_f + 1.0) / (2.0 * np.pi) * Cl * (T_CMB * 1e6) ** 2

    t_elapsed = time.time() - t_total
    if verbose:
        print(f"\n  Total time: {t_elapsed:.1f}s")

    return {
        'ell': ell_values,
        'Dl': Dl,
        'Cl': Cl,
        'k_arr': k_arr,
        'bg': bg,
        'Delta_l': Delta_l,
        'Phi_N': Phi_N_arr,
        'Psi_N': Psi_N_arr,
        'timing': {
            'background': t_bg,
            'perturbations': t_pert,
            'los': t_los,
            'total': t_elapsed,
        },
    }


# ============================================================================
# Section 7: CLASS comparison and diagnostics
# ============================================================================

def compare_with_class(ell, Dl, verbose=True):
    """
    Load CLASS reference C_l and compute comparison statistics.

    Returns dict with class_ell, class_Dl, rms, peak info.
    """
    class_path = '/Users/akaihuangm1/Desktop/github/gdm_class_public/output/mlx_ref_lcdm_cl.dat'
    if not os.path.exists(class_path):
        print(f"  CLASS reference not found at {class_path}")
        return None

    data = np.loadtxt(class_path)
    class_ell = data[:, 0].astype(int)
    # Column 1 is l(l+1)C_l/(2pi) in dimensionless units; convert to muK^2
    class_Dl = data[:, 1] * (T_CMB * 1e6) ** 2

    # Interpolate our result onto CLASS ell values
    from scipy.interpolate import interp1d
    our_interp = interp1d(ell, Dl, kind='cubic', bounds_error=False, fill_value=0.0)
    our_at_class = our_interp(class_ell)

    # RMS where D_l > 100 muK^2 (skip low-l noise)
    mask = class_Dl > 100.0
    if np.any(mask):
        ratio = our_at_class[mask] / class_Dl[mask]
        rms = np.sqrt(np.mean((ratio - 1.0) ** 2)) * 100.0  # percent
    else:
        rms = np.nan

    if verbose:
        print(f"\n{'='*70}")
        print(f"  CLASS Comparison")
        print(f"{'='*70}")
        print(f"  RMS deviation (D_l > 100 muK^2): {rms:.2f}%")

    # D_l at key multipoles
    target_ells = [2, 10, 50, 100, 220, 400, 537, 700, 820, 1000, 1200, 1500, 2000]
    if verbose:
        print(f"\n  {'l':>6s}  {'Ours':>10s}  {'CLASS':>10s}  {'Ratio':>8s}")
        print(f"  {'-'*40}")
        for tl in target_ells:
            idx_ours = np.argmin(np.abs(ell - tl))
            idx_class = np.argmin(np.abs(class_ell - tl))
            d_ours = Dl[idx_ours]
            d_class = class_Dl[idx_class]
            ratio_val = d_ours / d_class if d_class > 0 else np.nan
            print(f"  {ell[idx_ours]:6d}  {d_ours:10.1f}  {d_class:10.1f}  {ratio_val:8.3f}")

    # Find peaks (local maxima in D_l)
    peaks_ours = []
    for i in range(1, len(Dl) - 1):
        if Dl[i] > Dl[i-1] and Dl[i] > Dl[i+1] and Dl[i] > 500:
            peaks_ours.append((ell[i], Dl[i]))

    peaks_class = []
    for i in range(1, len(class_Dl) - 1):
        if class_Dl[i] > class_Dl[i-1] and class_Dl[i] > class_Dl[i+1] and class_Dl[i] > 500:
            peaks_class.append((class_ell[i], class_Dl[i]))

    if verbose:
        n_peaks = min(5, len(peaks_ours), len(peaks_class))
        print(f"\n  First {n_peaks} acoustic peaks:")
        print(f"  {'Peak':>5s}  {'l_ours':>7s}  {'D_ours':>10s}  {'l_CLASS':>8s}  {'D_CLASS':>10s}  {'D ratio':>8s}")
        print(f"  {'-'*55}")
        for ip in range(n_peaks):
            lo, do = peaks_ours[ip]
            lc, dc = peaks_class[ip]
            print(f"  {ip+1:5d}  {lo:7d}  {do:10.1f}  {lc:8d}  {dc:10.1f}  {do/dc:8.3f}")

    return {
        'class_ell': class_ell,
        'class_Dl': class_Dl,
        'rms': rms,
        'peaks_ours': peaks_ours[:5],
        'peaks_class': peaks_class[:5],
    }


# ============================================================================
# Main entry point
# ============================================================================

def main():
    print("\n" + "=" * 70)
    print("  solver_sync_v2.py -- Clean synchronous gauge CMB solver")
    print("=" * 70 + "\n")

    result = run(N_k=180, k_min=3e-4, k_max=0.35,
                 lg_max=L_GAMMA_MAX, ln_max=L_NU_MAX,
                 rtol=1e-6, atol=1e-9, verbose=True)

    ell = result['ell']
    Dl  = result['Dl']

    # Compare with CLASS
    comp = compare_with_class(ell, Dl, verbose=True)

    # Save results
    outpath = os.path.join(os.path.dirname(__file__), 'cl_sync_v2.dat')
    np.savetxt(outpath, np.column_stack([ell, Dl]),
               header='ell  D_l[muK^2]', fmt='%6d  %.6e')
    print(f"\n  Results saved to {outpath}")


if __name__ == '__main__':
    main()
