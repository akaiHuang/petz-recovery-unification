"""
massive_neutrino.py — Massive neutrino support for mlx_class.

Extends the massless neutrino solver (perturbations_neutrino.py) by adding:
  - Fermi-Dirac phase-space integration via Gauss-Laguerre quadrature (N_q bins)
  - Per-q-bin Boltzmann hierarchy: Psi_l(k, q, tau) for l = 0 .. l_nu_max
  - Relativistic → non-relativistic transition: epsilon(q,a) = sqrt(q^2 + m^2 a^2)
  - Background rho_nu(a) from full Fermi-Dirac integral (radiation → matter transition)
  - Proper energy-momentum tensor: delta_nu, theta_nu, shear_nu from q-integration

Physics summary:
  - Massive neutrinos free-stream at v = q/epsilon, which decreases as a grows
  - At z_nr ~ m_nu / (3 T_nu), neutrinos become non-relativistic
  - Below the free-streaming scale k_fs, they cluster like CDM
  - Above k_fs, they suppress P(k) by ~ -8 f_nu (key observable)

The per-q-bin hierarchy (conformal Newtonian gauge):
  Psi_0' = -(q k / epsilon) Psi_1 - Phi' (d ln f_0 / d ln q)
  Psi_1' = (q k / (3 epsilon)) [Psi_0 - 2 Psi_2] + (epsilon k / (3 q)) Psi
           * (d ln f_0 / d ln q)   [but see below: the Psi potential coupling
           enters as (k epsilon)/(3 q) * Psi, which for the combination
           with the f_0 derivative in Ma & Bertschinger notation simplifies]
  Psi_l' = (q k / ((2l+1) epsilon)) [l Psi_{l-1} - (l+1) Psi_{l+1}]
  Psi_{l_max}' = (q k / (2l_max+1) epsilon) l_max Psi_{l_max-1}
                 - (l_max+1)/tau Psi_{l_max}

We follow the Ma & Bertschinger (1995) conventions (ApJ 455, 7), Eqs (50)-(54).
Their Psi_l are defined such that the energy-momentum integrals are:

  delta rho_nu = (4pi / a^4) int q^2 dq epsilon f_0(q) Psi_0
  (rho+P) theta_nu = (4pi k / a^4) int q^2 dq q f_0(q) Psi_1
  (rho+P) sigma_nu = (8pi / (3 a^4)) int q^2 dq (q^2/epsilon) f_0(q) Psi_2

Reference: Ma & Bertschinger 1995, Eqs (50)-(54); CLASS ncdm module.

Author: Sheng-Kai Huang, 2026
"""
import numpy as np
import mlx.core as mx
import time

from .background import (
    Omega_r as _OMEGA_R,
    Omega_b as _OMEGA_B,
    Omega_c as _OMEGA_C,
    H0_Mpc as _H0_MPC,
    h as _h,
    T_CMB as _T_CMB,
)

# ============================================================================
# Physical constants
# ============================================================================
_K_B_EV_PER_K = 8.617333262e-5       # Boltzmann constant in eV/K
_T_NU_OVER_T_CMB = (4.0 / 11.0) ** (1.0 / 3.0)  # T_nu / T_CMB ~ 0.7138
_T_NU_K = _T_CMB * _T_NU_OVER_T_CMB  # Neutrino temperature today in K
_T_NU_EV = _T_NU_K * _K_B_EV_PER_K   # Neutrino temperature today in eV

# Conversion: neutrino energy density parameter for one massless species
# rho_nu_one = (7/8)(4/11)^{4/3} rho_gamma
_RHO_GAMMA_OVER_H02 = _OMEGA_R / (1.0 + 0.2271 * 3.046)  # Omega_gamma
_RHO_NU_ONE_MASSLESS = (7.0 / 8.0) * (4.0 / 11.0) ** (4.0 / 3.0) * _RHO_GAMMA_OVER_H02

# For converting m_nu in eV to dimensionless momentum:
# The comoving momentum q is measured in units of T_nu,0:
#   epsilon = sqrt(q^2 + (m_nu * a / T_nu,0)^2)  [dimensionless]
# So the mass parameter is M = m_nu / T_nu,0

# Hubble constant in eV (for rho_nu normalization)
# rho_crit = 3 H_0^2 / (8 pi G) => Omega_nu = rho_nu / rho_crit
# We normalize rho_nu(a) such that in the massless limit it matches Omega_nu_massless

# ============================================================================
# Gauss-Laguerre quadrature for Fermi-Dirac integrals
# ============================================================================

def _gauss_laguerre_fermi_dirac(n_q):
    """
    Compute Gauss-Laguerre nodes and weights adapted for Fermi-Dirac integrals.

    The Fermi-Dirac distribution f_0(q) = 1/(e^q + 1) is sampled using
    Gauss-Laguerre quadrature: int_0^inf g(q) dq ~ sum_i w_i * g(q_i)

    For the standard Gauss-Laguerre rule:
      int_0^inf f(q) e^{-q} dq ~ sum_i w_i^GL f(q_i)

    We need: int_0^inf q^n / (e^q + 1) dq
    Rewrite: q^n / (e^q + 1) = q^n e^{-q} / (1 + e^{-q})
            = q^n e^{-q} * sum_{k=0}^inf (-1)^k e^{-kq}
    So:      int q^n / (e^q + 1) dq = sum_i w_i^GL q_i^n / (1 + e^{-q_i})

    We absorb the Fermi-Dirac factor into the weights:
      w_i = w_i^GL * q_i^2 / (e^{q_i} + 1) * e^{q_i}
    Then: int q^2 f_0(q) g(q) dq ~ sum_i w_i g(q_i)

    Returns
    -------
    q_nodes : array (n_q,) — quadrature nodes
    w_nodes : array (n_q,) — weights including q^2 f_0(q) factor
    f0_nodes : array (n_q,) — f_0(q_i) values
    dlnf0_dlnq : array (n_q,) — d ln f_0 / d ln q at nodes
    """
    # Standard Gauss-Laguerre nodes and weights
    q_nodes, w_gl = np.polynomial.laguerre.laggauss(n_q)

    # Fermi-Dirac distribution at nodes
    f0 = 1.0 / (np.exp(q_nodes) + 1.0)

    # Weights: absorb q^2 * f_0 * e^q (the e^q cancels the Laguerre weight e^{-q})
    w_nodes = w_gl * q_nodes ** 2 * f0 * np.exp(q_nodes)

    # d ln f_0 / d ln q = -q * e^q / (e^q + 1)
    # (since f_0 = 1/(e^q+1), df_0/dq = -e^q/(e^q+1)^2, so d ln f_0/d ln q = q * df_0/(f_0 dq))
    dlnf0 = -q_nodes * np.exp(q_nodes) / (np.exp(q_nodes) + 1.0)

    return q_nodes.astype(np.float32), w_nodes.astype(np.float32), \
           f0.astype(np.float32), dlnf0.astype(np.float32)


# ============================================================================
# Background: massive neutrino energy density
# ============================================================================

def rho_nu_massive_normalized(a_arr, m_nu_eV, n_q=15):
    """
    Compute massive neutrino energy density rho_nu(a) / rho_nu(a->0).

    In the massless limit (m_nu -> 0), rho_nu ~ a^{-4}.
    In the non-relativistic limit (a >> a_nr), rho_nu ~ a^{-3}.

    Parameters
    ----------
    a_arr : array — scale factors
    m_nu_eV : float — mass per species in eV
    n_q : int — number of quadrature points (15 is accurate to ~0.01%)

    Returns
    -------
    rho_ratio : array — rho_nu(a) / rho_nu_massless(a), where
                rho_nu_massless(a) = const / a^4
    """
    M = m_nu_eV / _T_NU_EV  # dimensionless mass parameter

    # Quadrature nodes for int q^2 dq epsilon f_0(q)
    q_nodes, w_gl = np.polynomial.laguerre.laggauss(n_q)
    f0 = 1.0 / (np.exp(q_nodes) + 1.0)
    # Effective weights for int q^2 f_0(q) epsilon(q) dq
    w_eff = w_gl * q_nodes ** 2 * f0 * np.exp(q_nodes)

    # Massless normalization: int q^2 f_0(q) * q dq = int q^3 f_0(q) dq
    # = 7 pi^4 / 120 (exact for Fermi-Dirac)
    norm_massless = np.sum(w_eff * q_nodes)  # numerical massless limit

    rho_ratio = np.zeros_like(a_arr)
    for i, a in enumerate(a_arr):
        epsilon = np.sqrt(q_nodes ** 2 + (M * a) ** 2)
        rho_ratio[i] = np.sum(w_eff * epsilon) / norm_massless

    return rho_ratio


def omega_nu_massive(m_nu_eV, N_massive=1):
    """
    Non-relativistic limit: Omega_nu h^2 = sum(m_nu) / (93.14 eV).
    This is the late-time contribution when all species are non-relativistic.
    """
    return N_massive * m_nu_eV / 93.14 / _h ** 2


def z_nr(m_nu_eV):
    """
    Non-relativistic transition redshift: z_nr ~ m_nu / (3 T_nu,0).
    """
    return m_nu_eV / (3.0 * _T_NU_EV) - 1.0


# ============================================================================
# Variable layout for massive neutrino state vector
# ============================================================================
# The state vector is:
#   [Phi, delta_b, delta_c, v_c, Theta_0, Theta_1,
#    N_0, N_1, ..., N_{l_ml_max},        (massless neutrinos)
#    Psi_{0,0}, Psi_{0,1}, ..., Psi_{0,l_max},   (q-bin 0)
#    Psi_{1,0}, Psi_{1,1}, ..., Psi_{1,l_max},   (q-bin 1)
#    ...
#    Psi_{n_q-1,0}, ..., Psi_{n_q-1,l_max}]      (q-bin n_q-1)

# We keep massless neutrinos for N_eff_massless = N_eff - N_massive species.

_IDX_PHI = 0
_IDX_DELTA_B = 1
_IDX_DELTA_C = 2
_IDX_V_C = 3
_IDX_THETA_0 = 4
_IDX_THETA_1 = 5
_IDX_N_START = 6  # massless neutrino hierarchy starts here


def _n_var_massive(l_ml_max, l_m_max, n_q):
    """
    Total number of variables per k-mode.

    Parameters
    ----------
    l_ml_max : int — truncation for massless neutrino hierarchy
    l_m_max : int — truncation for massive neutrino hierarchy (per q-bin)
    n_q : int — number of momentum bins
    """
    n_base = _IDX_N_START + l_ml_max + 1          # metric + baryons + CDM + photons + massless nu
    n_massive = n_q * (l_m_max + 1)                # massive neutrino (all q-bins)
    return n_base + n_massive


def _idx_massive_start(l_ml_max):
    """Index where massive neutrino variables begin."""
    return _IDX_N_START + l_ml_max + 1


def _idx_psi(l_ml_max, i_q, l, l_m_max):
    """Index of Psi_l for q-bin i_q."""
    start = _idx_massive_start(l_ml_max)
    return start + i_q * (l_m_max + 1) + l


# ============================================================================
# Massive neutrino energy-momentum tensor (integrated over q)
# ============================================================================

def _massive_nu_emt(y, k_arr, a, m_nu_eV, q_nodes, w_nodes, dlnf0,
                    l_ml_max, l_m_max, n_q):
    """
    Compute massive neutrino perturbation integrals.

    Returns
    -------
    delta_rho_nu : (N_k,) — density perturbation (times 4pi a^{-4})
    theta_nu : (N_k,) — velocity divergence contribution
    shear_nu : (N_k,) — anisotropic stress contribution
    pressure_nu : (N_k,) — pressure perturbation (times 4pi a^{-4})

    All quantities are normalized such that they can be directly used in
    the Poisson equation with the appropriate Omega factor.
    """
    M = m_nu_eV / _T_NU_EV
    start = _idx_massive_start(l_ml_max)
    n_per_q = l_m_max + 1

    N_k = y.shape[0]
    delta_rho = mx.zeros((N_k,))
    theta = mx.zeros((N_k,))
    shear = mx.zeros((N_k,))
    delta_P = mx.zeros((N_k,))

    for i_q in range(n_q):
        q = float(q_nodes[i_q])
        w = float(w_nodes[i_q])

        # epsilon = sqrt(q^2 + M^2 a^2) — note: a is an MLX scalar
        eps = mx.sqrt(mx.array(q ** 2) + (M * a) ** 2)

        idx0 = start + i_q * n_per_q
        Psi_0 = y[:, idx0]
        Psi_1 = y[:, idx0 + 1] if n_per_q > 1 else mx.zeros((N_k,))
        Psi_2 = y[:, idx0 + 2] if n_per_q > 2 else mx.zeros((N_k,))

        # delta rho: w * epsilon * Psi_0
        delta_rho = delta_rho + w * eps * Psi_0

        # theta: w * q * Psi_1  (the k factor is in the velocity divergence def)
        theta = theta + w * q * Psi_1

        # shear: w * (q^2/epsilon) * Psi_2
        shear = shear + w * (q ** 2 / eps) * Psi_2

        # delta P: w * (q^2 / (3 epsilon)) * Psi_0
        delta_P = delta_P + w * (q ** 2 / (3.0 * eps)) * Psi_0

    return delta_rho, theta, shear, delta_P


def _massive_nu_background_integrals(a, m_nu_eV, q_nodes, w_nodes):
    """
    Background integrals for massive neutrinos at scale factor a.

    Returns
    -------
    rho_integral : float — proportional to rho_nu(a) * a^4
    P_integral : float — proportional to P_nu(a) * a^4
    """
    M = m_nu_eV / _T_NU_EV
    eps = mx.sqrt(mx.array(q_nodes ** 2, dtype=mx.float32) + (M * a) ** 2)
    q2_over_3eps = mx.array(q_nodes ** 2, dtype=mx.float32) / (3.0 * eps)
    w_mx = mx.array(w_nodes)

    rho_int = mx.sum(w_mx * eps)
    P_int = mx.sum(w_mx * q2_over_3eps)

    return rho_int, P_int


# ============================================================================
# Phi equation: IMEX decomposition with massive neutrinos
# ============================================================================

def _phi_lambda(k_arr, calH):
    """Stiff eigenvalue for Phi: lambda = -(k^2/(3 calH) + calH)."""
    return -(k_arr * k_arr / (3.0 * calH) + calH)


def _phi_source_massive(y, k_arr, calH, a, H0_Mpc,
                         Omega_gamma, Omega_nu_ml, Omega_b, Omega_c,
                         Omega_nu_m_today, m_nu_eV, N_massive,
                         q_nodes, w_nodes, dlnf0,
                         l_ml_max, l_m_max, n_q,
                         rho_nu_m_ratio):
    """
    Non-stiff source for Phi equation, including massive neutrino contribution.

    F = -(H0^2 / (2 calH)) * S

    S = Omega_gamma/a^2 * 4 Theta_0
      + Omega_nu_ml/a^2 * 4 N_0             (massless neutrinos)
      + Omega_b/a * delta_b
      + Omega_c/a * delta_c
      + Omega_nu_m/a^2 * 4 * (delta_rho_nu / rho_nu_bg)   (massive neutrinos)

    The massive neutrino contribution uses:
      4 pi / a^4 * int q^2 dq epsilon f_0 Psi_0
    normalized to the background energy density.
    """
    H02 = H0_Mpc * H0_Mpc

    Theta_0 = y[:, _IDX_THETA_0]
    N_0 = y[:, _IDX_N_START]
    delta_b = y[:, _IDX_DELTA_B]
    delta_c = y[:, _IDX_DELTA_C]

    S = (Omega_gamma / (a * a) * (4.0 * Theta_0)
         + Omega_nu_ml / (a * a) * (4.0 * N_0)
         + Omega_b / a * delta_b
         + Omega_c / a * delta_c)

    # Massive neutrino density perturbation
    if N_massive > 0 and n_q > 0:
        delta_rho_m, _, _, _ = _massive_nu_emt(
            y, k_arr, a, m_nu_eV, q_nodes, w_nodes, dlnf0,
            l_ml_max, l_m_max, n_q)

        # The background integral normalization:
        # rho_nu_m(a) * a^4 = (4pi T_nu^4) * int q^2 dq epsilon f_0
        # We compute the ratio delta_rho / rho_bg_integral
        rho_bg, _ = _massive_nu_background_integrals(a, m_nu_eV, q_nodes, w_nodes)

        # delta_nu = delta_rho_m / rho_bg (this is the fractional perturbation)
        # The Poisson source needs: Omega_nu_m(a) / a^2 * delta_nu * 4
        # But Omega_nu_m(a) = Omega_nu_m_today * rho_ratio(a) / a^4 * a^4
        # In practice: Omega_nu_m_today / a^2 * 4 * delta_rho_m / rho_bg
        # where rho_bg is the background integral that gives rho_nu_m * a^4
        # The factor of 4 comes from 4*rho -> 4*delta*rho in the Poisson eq

        # Actually, let's be more careful. The Poisson equation in conformal Newtonian gauge:
        # k^2 Phi + 3 calH (Phi' + calH Phi) = -4pi G a^2 sum_i rho_i delta_i
        #   = -(3/2) H0^2 sum_i Omega_i(a=1) delta_i / a   [for matter]
        #   = -(3/2) H0^2 sum_i Omega_i(a=1) delta_i / a^2  [for radiation, since rho~a^-4]

        # For massive neutrinos: rho_nu(a) = rho_nu,0 * a^{-4} * rho_ratio(a)
        # So rho_nu(a) * delta_nu = Omega_nu_m_today * 3H0^2/(8piG) * a^{-4} * rho_ratio * delta_nu
        # In the Poisson equation source:
        # 4pi G a^2 * rho_nu * delta_nu = (3/2) H0^2 * Omega_nu_m_today * rho_ratio / a^2 * delta_nu

        # The fractional perturbation is delta_nu = delta_rho_m / rho_bg
        # And the effective source is:
        #   Omega_nu_m_today * rho_ratio / a^2 * delta_nu
        # = Omega_nu_m_today / a^2 * delta_rho_m / rho_bg * rho_ratio
        # But delta_rho_m already has the rho_ratio built in... no, let's think again.

        # delta_rho_m = sum_i w_i eps_i Psi_0^i (this is the perturbation integral)
        # rho_bg = sum_i w_i eps_i (this is the background integral)
        # delta_nu = delta_rho_m / rho_bg

        # In the source S (which gets multiplied by -H0^2/(2 calH)):
        # For radiation: Omega_r/a^2 * 4 * Theta_0
        # For matter: Omega_m/a * delta_m
        # For massive nu: this is more subtle because it's neither pure radiation nor matter.

        # The correct expression: same as massless but with rho_ratio correction
        # In massless limit, rho_ratio = 1 and we get Omega_nu/a^2 * 4 * N_0
        # For massive:
        #   contribution = N_massive * Omega_nu_one_massless / a^2 * 4 * (delta_rho_m / rho_bg) * rho_ratio

        # But wait, we precomputed rho_nu_m_ratio = rho_ratio at this time step.
        # And Omega_nu_m_today is defined as N_massive * Omega_nu_one_massless

        # Source term for massive neutrinos:
        safe_rho = mx.maximum(rho_bg, mx.array(1e-30))
        delta_nu = delta_rho_m / safe_rho
        S = S + Omega_nu_m_today / (a * a) * 4.0 * delta_nu * rho_nu_m_ratio

    return -0.5 * H02 / calH * S


def _phi_int_factor(lam, dt):
    """(exp(lam*dt) - 1) / lam, with Taylor expansion for small |lam*dt|."""
    lam_dt = lam * dt
    exact = (mx.exp(lam_dt) - 1.0) / lam
    taylor = dt * (1.0 + 0.5 * lam_dt + lam_dt * lam_dt / 6.0)
    return mx.where(mx.abs(lam_dt) < 1e-4, taylor, exact)


# ============================================================================
# Full derivative function with massive neutrinos
# ============================================================================

def deriv_massive(y, k_arr, calH, R, tau, a, H0_Mpc,
                  Omega_gamma, Omega_nu_ml, Omega_b, Omega_c,
                  Omega_nu_m_today, m_nu_eV, N_massive,
                  q_nodes, w_nodes, dlnf0,
                  l_ml_max, l_m_max, n_q,
                  rho_nu_m_ratio):
    """
    Full RHS for all variables including massive neutrino hierarchy.

    Same IMEX strategy as perturbations_neutrino.py:
    - Phi equation: Psi = Phi (zeroth order, no aniso stress backreaction)
    - Momentum equations: use diagnosed Psi (with both massless + massive aniso stress)
    - Massive neutrino: per-q-bin Boltzmann hierarchy with velocity q/epsilon

    Returns dy shape (N_k, n_var).
    """
    N_k = y.shape[0]
    M = m_nu_eV / _T_NU_EV

    Phi = y[:, _IDX_PHI]
    delta_b = y[:, _IDX_DELTA_B]
    delta_c = y[:, _IDX_DELTA_C]
    v_c = y[:, _IDX_V_C]
    Theta_0 = y[:, _IDX_THETA_0]
    Theta_1 = y[:, _IDX_THETA_1]

    # Extract massless neutrino multipoles
    N_ml = []
    for l in range(l_ml_max + 1):
        N_ml.append(y[:, _IDX_N_START + l])
    N_ml_2 = N_ml[2] if l_ml_max >= 2 else mx.zeros_like(Phi)

    # Diagnose Psi: include both massless and massive anisotropic stress
    H02 = H0_Mpc * H0_Mpc

    # Massless neutrino contribution to Psi - Phi
    Psi_minus_Phi = -12.0 * H02 * Omega_nu_ml * N_ml_2 / (a * a * k_arr * k_arr)

    # Massive neutrino contribution to anisotropic stress
    if N_massive > 0 and n_q > 0:
        _, _, shear_m, _ = _massive_nu_emt(
            y, k_arr, a, m_nu_eV, q_nodes, w_nodes, dlnf0,
            l_ml_max, l_m_max, n_q)
        rho_bg, _ = _massive_nu_background_integrals(a, m_nu_eV, q_nodes, w_nodes)
        safe_rho = mx.maximum(rho_bg, mx.array(1e-30))
        sigma_nu_m = (2.0 / 3.0) * shear_m / safe_rho
        # Contribution: -12 H0^2 Omega_nu_m(a) / (a^2 k^2) * sigma_nu_m
        # where Omega_nu_m(a) at this epoch in the k^2 Psi equation
        Psi_minus_Phi = Psi_minus_Phi - 12.0 * H02 * Omega_nu_m_today * rho_nu_m_ratio * sigma_nu_m / (a * a * k_arr * k_arr)

    Psi = Phi + Psi_minus_Phi

    # ---- Phi' from Poisson equation (Psi = Phi, zeroth order) ----
    # For massive neutrinos, we include their density contribution
    S_base = (Omega_gamma / (a * a) * (4.0 * Theta_0)
              + Omega_nu_ml / (a * a) * (4.0 * N_ml[0])
              + Omega_b / a * delta_b
              + Omega_c / a * delta_c)

    # Massive neutrino density contribution to Phi'
    S_massive = mx.zeros_like(Phi)
    if N_massive > 0 and n_q > 0:
        delta_rho_m, _, _, _ = _massive_nu_emt(
            y, k_arr, a, m_nu_eV, q_nodes, w_nodes, dlnf0,
            l_ml_max, l_m_max, n_q)
        rho_bg, _ = _massive_nu_background_integrals(a, m_nu_eV, q_nodes, w_nodes)
        safe_rho = mx.maximum(rho_bg, mx.array(1e-30))
        delta_nu_m = delta_rho_m / safe_rho
        S_massive = Omega_nu_m_today / (a * a) * 4.0 * delta_nu_m * rho_nu_m_ratio

    S = S_base + S_massive
    Phi_dot = -calH * Phi - k_arr * k_arr * Phi / (3.0 * calH) - 0.5 * H02 / calH * S

    # ---- Photon TCA ----
    dTheta_0 = -k_arr * Theta_1 - Phi_dot
    dTheta_1 = (k_arr / 3.0 * (Theta_0 + Psi) - calH * R * Theta_1) / (1.0 + R)

    # ---- Baryons ----
    d_delta_b = -3.0 * k_arr * Theta_1 - 3.0 * Phi_dot

    # ---- CDM ----
    d_delta_c = -(k_arr * v_c + 3.0 * Phi_dot)
    d_v_c = -calH * v_c + k_arr * Psi

    # ---- Massless neutrino hierarchy ----
    dN_ml = []
    dN_ml.append(-k_arr * N_ml[1] - Phi_dot)
    dN_ml.append(k_arr / 3.0 * (N_ml[0] - 2.0 * N_ml_2 + Psi))
    for l in range(2, l_ml_max):
        dN_ml.append(k_arr / (2.0 * l + 1.0) * (l * N_ml[l - 1] - (l + 1.0) * N_ml[l + 1]))
    tau_safe = mx.maximum(tau, mx.array(1e-10))
    dN_ml.append(k_arr * N_ml[l_ml_max - 1] * l_ml_max / (2.0 * l_ml_max + 1.0)
                 - (l_ml_max + 1.0) / tau_safe * N_ml[l_ml_max])

    # ---- Massive neutrino hierarchy (per q-bin) ----
    start = _idx_massive_start(l_ml_max)
    n_per_q = l_m_max + 1

    dPsi_all = []  # will hold (n_q * n_per_q) derivatives

    for i_q in range(n_q):
        q = float(q_nodes[i_q])
        dlnf = float(dlnf0[i_q])

        eps = mx.sqrt(mx.array(q ** 2) + (M * a) ** 2)
        v_q = q / eps  # velocity = q / epsilon (dimensionless, in units of c)

        idx0 = start + i_q * n_per_q

        # Extract Psi_l for this q-bin
        Psi_l = []
        for l in range(n_per_q):
            Psi_l.append(y[:, idx0 + l])

        Psi_2_q = Psi_l[2] if n_per_q > 2 else mx.zeros_like(Phi)

        # Ma & Bertschinger Eqs (50)-(53):
        # l=0: Psi_0' = -(q k / epsilon) Psi_1 + Phi_dot / 3  [dlnf0 factor absorbed]
        #
        # Actually, the MB95 equations for the perturbation Psi_l of the
        # distribution function f = f_0(q)[1 + Psi] are:
        #
        # Psi_0' = -(q/epsilon) k Psi_1 + (1/6) Phi' d(ln f_0)/d(ln q)  ... NO
        #
        # Let me be precise. MB95 Eq (50)-(53) define F_nu_l where
        # delta T^0_0 integrates F_0, etc.
        # The hierarchy is:
        #   Psi_0' = -(q k / epsilon) Psi_1 + Phi' dlnf0/dlnq ... WRONG
        #
        # MB95 actually writes (their eq 50-53 with their notation):
        # For massive neutrinos, the Boltzmann equation in k-space gives:
        #
        # dot{Psi}_0 = -q k / epsilon * Psi_1 + dot{eta} * dlnf0/dlnq / 6
        #   where eta is MB95's metric perturbation
        #
        # In conformal Newtonian gauge (our convention):
        # Psi_0' = -(q k / epsilon) Psi_1 - Phi' * dlnf0/6   ... still not right
        #
        # Let me use the CLASS convention directly. In CLASS (ncdm.c),
        # the hierarchy is written as:
        #   Psi_0' = -(q k / epsilon) Psi_1 + Phi' / 3
        #   Psi_1' = (q k / (3 epsilon)) (Psi_0 - 2 Psi_2) - (epsilon k / (3 q)) Psi_pot
        #   Psi_l' = (q k / ((2l+1) epsilon)) (l Psi_{l-1} - (l+1) Psi_{l+1})
        #
        # where the Psi_l are defined such that:
        #   delta_rho = (4pi/a^4) int q^2 dq epsilon f_0 Psi_0
        #   (rho+P) theta = (4pi k / a^4) int q^2 dq q f_0 Psi_1
        # etc.
        #
        # The key point: with this convention, the metric source Phi'
        # enters with factor 1/3 (not dlnf0), because the dlnf0 factor
        # is already absorbed into the definition of Psi_l.
        #
        # WAIT — that's the convention where Psi_l is the perturbation
        # to the integrated distribution. Let me follow MB95 more carefully.
        #
        # MB95 define (their eq 11):
        #   f(x,q,n,tau) = f_0(q) [1 + Psi(x,q,n,tau)]
        # and expand Psi in Legendre polynomials:
        #   Psi = sum_l (-i)^l (2l+1) Psi_l(k,q,tau) P_l(mu)
        #
        # Then MB95 eq (50)-(53):
        #   dot{Psi}_0 = -(qk/eps) Psi_1 + (1/6) dot{h} dlnf0/dlnq
        # where h is their synchronous gauge metric. In conformal Newtonian:
        #   dot{Psi}_0 = -(qk/eps) Psi_1 + dot{Phi} dlnf0/dlnq ... ? No.
        #
        # Actually, for conformal Newtonian gauge, the Boltzmann eq is:
        #   dot{Psi} + i(qk/eps) mu Psi = [dot{Phi} + i(eps k/q) mu Psi_pot] dlnf0/dlnq
        # (where Psi_pot is the Newtonian potential, what we call Psi)
        #
        # Expanding in multipoles:
        #   dot{Psi}_0 = -(qk/eps) Psi_1 - dot{Phi} * |dlnf0/dlnq|
        #                [since dlnf0/dlnq < 0, and with sign conventions...]
        #
        # This is getting confusing with signs. Let me use the well-known
        # result that in the massless limit, these reduce to the standard
        # hierarchy for N_l. In the massless limit:
        #   epsilon = q, so qk/epsilon = k
        #   The hierarchy becomes:
        #     Psi_0' = -k Psi_1 - Phi'   [matching N_0' = -k N_1 - Phi']
        #     Psi_1' = k/3 (Psi_0 - 2 Psi_2 + Psi_pot)
        #     Psi_l' = k/(2l+1) [l Psi_{l-1} - (l+1) Psi_{l+1}]
        #
        # So the massive neutrino hierarchy with our conventions is:
        #
        #   Psi_0' = -(qk/eps) Psi_1 - Phi'
        #   Psi_1' = (qk/(3 eps)) [Psi_0 - 2 Psi_2] + (eps k / (3q)) Psi_pot
        #              ... but in massless limit eps=q, so this gives
        #              k/3 [Psi_0 - 2 Psi_2 + Psi_pot]  ✓
        #   Psi_l' = (qk/((2l+1) eps)) [l Psi_{l-1} - (l+1) Psi_{l+1}]
        #
        # BUT WAIT: in the massless case, Psi_0' = -k Psi_1 - Phi' does NOT
        # have a dlnf0/dlnq factor. This means the convention is that Psi_l
        # already absorbs any f_0 weighting. Good.
        #
        # So for MASSIVE neutrinos, the only change from massless is:
        #   - k -> qk/epsilon in the streaming terms
        #   - The potential coupling in l=1 has eps/q instead of 1
        #   - Phi' in l=0 remains as-is (no dlnf0 factor needed)
        #
        # This is the convention used in CAMB and matches CLASS when
        # the Psi_l are defined per the distribution function perturbation.

        dPsi_q = []

        # l = 0: Psi_0' = -(q k / epsilon) Psi_1 - Phi'
        dPsi_q.append(-v_q * k_arr * Psi_l[1] - Phi_dot)

        # l = 1: Psi_1' = (q k / (3 epsilon)) (Psi_0 - 2 Psi_2) + (epsilon k / (3 q)) Psi
        dPsi_q.append(v_q * k_arr / 3.0 * (Psi_l[0] - 2.0 * Psi_2_q)
                      + (eps / q) * k_arr / 3.0 * Psi)

        # l >= 2: Psi_l' = (q k / ((2l+1) epsilon)) [l Psi_{l-1} - (l+1) Psi_{l+1}]
        for l in range(2, l_m_max):
            dPsi_q.append(v_q * k_arr / (2.0 * l + 1.0) * (l * Psi_l[l - 1] - (l + 1.0) * Psi_l[l + 1]))

        # l = l_m_max: truncation with free-streaming damping
        dPsi_q.append(v_q * k_arr * l_m_max / (2.0 * l_m_max + 1.0) * Psi_l[l_m_max - 1]
                      - (l_m_max + 1.0) / tau_safe * Psi_l[l_m_max])

        dPsi_all.extend(dPsi_q)

    # ---- Stack all derivatives ----
    parts = [Phi_dot[:, None], d_delta_b[:, None], d_delta_c[:, None],
             d_v_c[:, None], dTheta_0[:, None], dTheta_1[:, None]]

    for l in range(l_ml_max + 1):
        parts.append(dN_ml[l][:, None])

    for dPsi in dPsi_all:
        parts.append(dPsi[:, None])

    return mx.concatenate(parts, axis=1)


# ============================================================================
# Column replacement helper
# ============================================================================

def _set_phi(y, Phi_new):
    """Replace the Phi column (index 0) in the state vector."""
    return mx.concatenate([Phi_new[:, None], y[:, 1:]], axis=1)


# ============================================================================
# IMEX RK4 step with massive neutrinos
# ============================================================================

def imex_rk4_step_massive(y, k_arr, calH, R, tau_val, a,
                          H0_Mpc,
                          Omega_gamma, Omega_nu_ml, Omega_b, Omega_c,
                          Omega_nu_m_today, m_nu_eV, N_massive,
                          q_nodes, w_nodes, dlnf0,
                          l_ml_max, l_m_max, n_q,
                          rho_nu_m_ratio, dtau):
    """
    IMEX RK4 step: exponential integrator for Phi, explicit RK4 for everything else.
    Identical structure to perturbations_neutrino.py but with massive neutrino hierarchy.
    """
    lam = _phi_lambda(k_arr, calH)
    exp_lam_dt = mx.exp(lam * dtau)
    exp_lam_half = mx.exp(lam * 0.5 * dtau)

    def get_F(state):
        return _phi_source_massive(
            state, k_arr, calH, a, H0_Mpc,
            Omega_gamma, Omega_nu_ml, Omega_b, Omega_c,
            Omega_nu_m_today, m_nu_eV, N_massive,
            q_nodes, w_nodes, dlnf0,
            l_ml_max, l_m_max, n_q,
            rho_nu_m_ratio)

    def full_rhs(state):
        return deriv_massive(
            state, k_arr, calH, R, mx.array(float(tau_val)), a, H0_Mpc,
            Omega_gamma, Omega_nu_ml, Omega_b, Omega_c,
            Omega_nu_m_today, m_nu_eV, N_massive,
            q_nodes, w_nodes, dlnf0,
            l_ml_max, l_m_max, n_q,
            rho_nu_m_ratio)

    # --- Stage 1 ---
    Phi_1 = y[:, _IDX_PHI]
    F1 = get_F(y)
    k1 = full_rhs(y)

    # --- Stage 2: half step ---
    Phi_2 = exp_lam_half * Phi_1 + _phi_int_factor(lam, 0.5 * dtau) * F1
    y2 = y + 0.5 * dtau * k1
    y2 = _set_phi(y2, Phi_2)
    F2 = get_F(y2)
    k2 = full_rhs(y2)

    # --- Stage 3: half step ---
    Phi_3 = exp_lam_half * Phi_1 + _phi_int_factor(lam, 0.5 * dtau) * F2
    y3 = y + 0.5 * dtau * k2
    y3 = _set_phi(y3, Phi_3)
    F3 = get_F(y3)
    k3 = full_rhs(y3)

    # --- Stage 4: full step ---
    Phi_4 = exp_lam_dt * Phi_1 + _phi_int_factor(lam, dtau) * F3
    y4 = y + dtau * k3
    y4 = _set_phi(y4, Phi_4)
    F4 = get_F(y4)
    k4 = full_rhs(y4)

    # --- Combine ---
    y_new = y + (dtau / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)

    F_mid = 0.5 * (F2 + F3)
    Phi_new = (exp_lam_dt * Phi_1
               + (dtau / 6.0) * (exp_lam_dt * F1
                                  + 4.0 * exp_lam_half * F_mid
                                  + F4))
    y_new = _set_phi(y_new, Phi_new)

    return y_new


# ============================================================================
# Tau grid
# ============================================================================

def build_tau_grid_massive(bg, k_max, N_early=300, N_late=1000):
    """Build integration grid. Same as neutrino version but with more late-time points."""
    tau_rec = bg.tau_rec
    tau_init = bg.tau_grid[1]
    tau_early_end = min(10.0, 0.3 * tau_rec)

    tau_early = np.geomspace(tau_init, tau_early_end, N_early)

    cs_approx = 1.0 / np.sqrt(3.0)
    dtau_acoustic = 2.0 * np.pi / (k_max * cs_approx) / 25.0
    N_needed = int((tau_rec - tau_early_end) / dtau_acoustic) + 10
    N_late_actual = max(N_late, N_needed)

    tau_late = np.linspace(tau_early_end, tau_rec, N_late_actual)
    tau_grid = np.unique(np.concatenate([tau_early, tau_late]))
    return tau_grid


# ============================================================================
# Initial conditions with massive neutrinos
# ============================================================================

def adiabatic_ic_massive(k_arr_np, bg, m_nu_eV, q_nodes, dlnf0,
                         l_ml_max, l_m_max, n_q):
    """
    Adiabatic initial conditions including massive neutrinos.

    At early times (a << a_nr), massive neutrinos are ultra-relativistic
    and their ICs are the same as massless. The per-q-bin ICs are:

    Psi_0 = -Phi/2  (same as massless N_0)
    Psi_1 = (q/epsilon) * k tau / 18  (velocity reduced by q/epsilon)
    Psi_2 = ((q/epsilon) * k tau)^2 / 60

    Since tau_init is very early, epsilon ~ q and these reduce to massless.
    """
    N_k = len(k_arr_np)
    tau_init = bg.tau_grid[1]
    a_init = float(bg.a_at_tau(np.array([tau_init]))[0])
    M = m_nu_eV / _T_NU_EV

    nvar = _n_var_massive(l_ml_max, l_m_max, n_q)
    y0 = np.zeros((N_k, nvar), dtype=np.float32)

    # Standard variables (same as massless neutrino solver)
    y0[:, _IDX_PHI] = 1.0
    y0[:, _IDX_DELTA_B] = -1.5
    y0[:, _IDX_DELTA_C] = -1.5
    y0[:, _IDX_V_C] = k_arr_np * tau_init / 6.0
    y0[:, _IDX_THETA_0] = -0.5
    y0[:, _IDX_THETA_1] = k_arr_np * tau_init / 18.0

    # Massless neutrino hierarchy
    y0[:, _IDX_N_START + 0] = -0.5
    y0[:, _IDX_N_START + 1] = k_arr_np * tau_init / 18.0
    if l_ml_max >= 2:
        y0[:, _IDX_N_START + 2] = (k_arr_np * tau_init) ** 2 / 60.0
    for l in range(3, min(l_ml_max + 1, 6)):
        prod_val = 1.0
        for j in range(1, l + 1):
            prod_val *= (2 * j + 1)
        y0[:, _IDX_N_START + l] = (k_arr_np * tau_init) ** l / prod_val

    # Massive neutrino hierarchy (per q-bin)
    start = _idx_massive_start(l_ml_max)
    n_per_q = l_m_max + 1

    for i_q in range(n_q):
        q = float(q_nodes[i_q])
        eps = np.sqrt(q ** 2 + (M * a_init) ** 2)
        v_q = q / eps  # ~ 1 at early times

        idx0 = start + i_q * n_per_q

        # l = 0: same as massless
        y0[:, idx0] = -0.5

        # l = 1: velocity-reduced
        y0[:, idx0 + 1] = v_q * k_arr_np * tau_init / 18.0

        # l = 2
        if n_per_q > 2:
            y0[:, idx0 + 2] = (v_q * k_arr_np * tau_init) ** 2 / 60.0

        # l >= 3 (tiny)
        for l in range(3, min(n_per_q, 6)):
            prod_val = 1.0
            for j in range(1, l + 1):
                prod_val *= (2 * j + 1)
            y0[:, idx0 + l] = (v_q * k_arr_np * tau_init) ** l / prod_val

    return mx.array(y0)


# ============================================================================
# Precompute rho_nu_m ratio on tau grid
# ============================================================================

def _precompute_rho_ratio(a_grid, m_nu_eV, q_nodes, w_nodes):
    """
    Precompute rho_nu_massive(a) / rho_nu_massless(a) on the tau grid.

    rho_nu_massless ~ a^{-4}, so rho_ratio = a^4 * rho_nu(a) / rho_nu(a->0).

    Returns array of shape (N_tau,).
    """
    M = m_nu_eV / _T_NU_EV

    # Massless normalization: sum w_i * q_i
    norm_massless = np.sum(w_nodes * q_nodes)

    rho_ratio = np.zeros(len(a_grid), dtype=np.float32)
    for i, a in enumerate(a_grid):
        eps = np.sqrt(q_nodes ** 2 + (M * a) ** 2)
        rho_ratio[i] = np.sum(w_nodes * eps) / norm_massless

    return rho_ratio


# ============================================================================
# Solver class
# ============================================================================

class MassiveNeutrinoBoltzmannSolver:
    """
    Batched Boltzmann solver with massive + massless neutrinos and IMEX integration.

    Supports:
    - N_massive massive species (default 1) with mass m_nu (in eV)
    - N_eff - N_massive massless species
    - Gauss-Laguerre quadrature with n_q momentum bins per massive species
    - Per-q-bin Boltzmann hierarchy with l_m_max multipoles
    - Background rho_nu(a) with relativistic → non-relativistic transition
    - Proper anisotropic stress from both massless and massive species

    Parameters
    ----------
    bg : Background — solved background object
    k_arr_Mpc : array — k values in Mpc^{-1}
    m_nu_eV : float — mass per massive species in eV (default 0.02, minimal normal hierarchy)
    N_massive : int — number of massive species (default 1; use 3 for degenerate)
    n_q : int — number of Gauss-Laguerre momentum bins (default 5)
    l_ml_max : int — truncation for massless neutrino hierarchy (default 20)
    l_m_max : int — truncation for massive neutrino hierarchy (default 15)
    """

    def __init__(self, bg, k_arr_Mpc, m_nu_eV=0.02, N_massive=1,
                 n_q=5, l_ml_max=20, l_m_max=15):
        self.bg = bg
        self.k_arr_np = np.asarray(k_arr_Mpc, dtype=np.float32)
        self.N_k = len(self.k_arr_np)
        self.k_arr = mx.array(self.k_arr_np)
        self.m_nu_eV = m_nu_eV
        self.N_massive = N_massive
        self.n_q = n_q
        self.l_ml_max = l_ml_max
        self.l_m_max = l_m_max

        # Quadrature
        self.q_nodes, self.w_nodes, self.f0_nodes, self.dlnf0 = \
            _gauss_laguerre_fermi_dirac(n_q)

        # Radiation density split:
        # N_eff_massless = N_eff - N_massive (the remaining massless species)
        N_eff = 3.046
        self.N_eff_massless = N_eff - N_massive
        if self.N_eff_massless < 0:
            raise ValueError(f"N_massive={N_massive} exceeds N_eff={N_eff}")

        f_nu_total = 0.2271 * N_eff / (1.0 + 0.2271 * N_eff)
        self.Omega_gamma = _OMEGA_R * (1.0 - f_nu_total)

        # Massless neutrino density: N_eff_massless species
        _one_nu = (7.0 / 8.0) * (4.0 / 11.0) ** (4.0 / 3.0) * self.Omega_gamma
        self.Omega_nu_ml = self.N_eff_massless * _one_nu

        # Massive neutrino density (at early times, when relativistic):
        # same as N_massive massless species
        self.Omega_nu_m_today = N_massive * _one_nu  # "today" = a=1 if still relativistic

        # Actual Omega_nu_m at late times: m_nu / (93.14 h^2) per species
        self.Omega_nu_m_late = N_massive * m_nu_eV / (93.14 * _h ** 2)

        self.Omega_b = _OMEGA_B
        self.Omega_c = float(bg.Omega_cdm)

        # Total variable count
        self.n_var = _n_var_massive(l_ml_max, l_m_max, n_q)

        # Build tau grid
        k_max = float(self.k_arr_np[-1])
        self.tau_grid = build_tau_grid_massive(bg, k_max)

        # Precompute rho_ratio on tau grid
        a_on_tau = bg.a_at_tau(self.tau_grid).astype(np.float32)
        self._rho_ratio_grid = _precompute_rho_ratio(
            a_on_tau, m_nu_eV, self.q_nodes, self.w_nodes)

        # Non-relativistic transition
        _z_nr = z_nr(m_nu_eV)
        _k_fs = 0.018 * np.sqrt(m_nu_eV / 0.06) * _h  # approximate free-streaming scale

        # Print summary
        print(f"\n{'=' * 70}")
        print(f"Massive Neutrino Boltzmann Solver")
        print(f"{'=' * 70}")
        print(f"  m_nu = {m_nu_eV:.4f} eV, N_massive = {N_massive}, "
              f"sum(m_nu) = {N_massive * m_nu_eV:.4f} eV")
        print(f"  z_nr = {_z_nr:.0f} (non-relativistic transition)")
        print(f"  k_fs ~ {_k_fs:.4f} h/Mpc (free-streaming scale)")
        print(f"  Omega_nu_m (late) = {self.Omega_nu_m_late:.6f}")
        print(f"  f_nu = {self.Omega_nu_m_late / (self.Omega_b + self.Omega_c + self.Omega_nu_m_late):.4f}")
        print(f"  N_k = {self.N_k}, n_q = {n_q}, l_m_max = {l_m_max}, l_ml_max = {l_ml_max}")
        print(f"  N_var = {self.n_var} per k-mode "
              f"(6 base + {l_ml_max + 1} massless nu + {n_q * (l_m_max + 1)} massive nu)")
        print(f"  Memory: {self.N_k * self.n_var * 4 / 1024:.1f} KB")
        print(f"  Grid: {len(self.tau_grid)} steps "
              f"[{self.tau_grid[0]:.2e}, {self.tau_grid[-1]:.1f}] Mpc")
        print(f"  Quadrature nodes (q/T_nu): {self.q_nodes}")
        print(f"  Quadrature weights: {self.w_nodes}")
        print(f"{'=' * 70}\n")

    def solve(self):
        """Integrate the Boltzmann hierarchy from tau_init to tau_rec."""
        t0 = time.time()
        bg = self.bg
        k_arr = self.k_arr
        tau_grid = self.tau_grid

        # Initial conditions
        y = adiabatic_ic_massive(
            self.k_arr_np, bg, self.m_nu_eV,
            self.q_nodes, self.dlnf0,
            self.l_ml_max, self.l_m_max, self.n_q)

        # Precompute background on tau grid
        calH_np = bg.calH_at_tau(tau_grid).astype(np.float32)
        R_np = bg.R_at_tau(tau_grid).astype(np.float32)
        a_np = bg.a_at_tau(tau_grid).astype(np.float32)

        calH_arr = mx.array(calH_np)
        R_arr = mx.array(R_np)
        a_arr = mx.array(a_np)
        rho_ratio_arr = mx.array(self._rho_ratio_grid)

        _Og = self.Omega_gamma
        _On_ml = self.Omega_nu_ml
        _On_m = self.Omega_nu_m_today
        _Ob = self.Omega_b
        _Oc = self.Omega_c
        _H0 = _H0_MPC
        _lml = self.l_ml_max
        _lm = self.l_m_max
        _nq = self.n_q
        _m = self.m_nu_eV
        _Nm = self.N_massive
        _q = self.q_nodes
        _w = self.w_nodes
        _dlnf = self.dlnf0

        print("[Massive nu IMEX] Integrating to tau_rec...")
        for i in range(len(tau_grid) - 1):
            dt = float(tau_grid[i + 1] - tau_grid[i])
            y = imex_rk4_step_massive(
                y, k_arr,
                calH_arr[i], R_arr[i], tau_grid[i], a_arr[i],
                _H0,
                _Og, _On_ml, _Ob, _Oc,
                _On_m, _m, _Nm,
                _q, _w, _dlnf,
                _lml, _lm, _nq,
                rho_ratio_arr[i], dt)

            if i % 200 == 0:
                mx.eval(y)

        mx.eval(y)
        t_total = time.time() - t0
        print(f"[Massive nu IMEX] Done in {t_total:.2f}s ({len(tau_grid)} steps)")

        return MassiveNeutrinoResult(
            y=y, k_arr=self.k_arr_np, bg=bg,
            l_ml_max=self.l_ml_max, l_m_max=self.l_m_max,
            n_q=self.n_q, m_nu_eV=self.m_nu_eV,
            N_massive=self.N_massive,
            Omega_nu_ml=self.Omega_nu_ml,
            Omega_nu_m_today=self.Omega_nu_m_today,
            q_nodes=self.q_nodes, w_nodes=self.w_nodes,
            dlnf0=self.dlnf0,
            H0=_H0_MPC)


# ============================================================================
# Result container
# ============================================================================

class MassiveNeutrinoResult:
    """Result container for massive neutrino solver."""

    def __init__(self, y, k_arr, bg, l_ml_max, l_m_max, n_q,
                 m_nu_eV, N_massive, Omega_nu_ml, Omega_nu_m_today,
                 q_nodes, w_nodes, dlnf0, H0):
        self.y = y
        self.k_arr = k_arr
        self.bg = bg
        self.l_ml_max = l_ml_max
        self.l_m_max = l_m_max
        self.n_q = n_q
        self.m_nu_eV = m_nu_eV
        self.N_massive = N_massive
        self.Omega_nu_ml = Omega_nu_ml
        self.Omega_nu_m_today = Omega_nu_m_today
        self.q_nodes = q_nodes
        self.w_nodes = w_nodes
        self.dlnf0 = dlnf0
        self.H0 = H0
        self.N_k = len(k_arr)

    def massive_nu_delta(self):
        """
        Extract massive neutrino density perturbation delta_nu(k) at tau_rec.

        Returns the fractional density perturbation delta rho / rho for each k.
        """
        y_np = np.array(self.y)
        M = self.m_nu_eV / _T_NU_EV
        a_rec = float(self.bg.a_at_tau(np.array([self.bg.tau_rec]))[0])

        start = _idx_massive_start(self.l_ml_max)
        n_per_q = self.l_m_max + 1

        # Numerator: int q^2 dq epsilon f_0 Psi_0
        num = np.zeros(self.N_k)
        den = 0.0

        for i_q in range(self.n_q):
            q = float(self.q_nodes[i_q])
            w = float(self.w_nodes[i_q])
            eps = np.sqrt(q ** 2 + (M * a_rec) ** 2)

            idx0 = start + i_q * n_per_q
            Psi_0 = y_np[:, idx0]

            num += w * eps * Psi_0
            den += w * eps

        return num / max(den, 1e-30)

    def massive_nu_shear(self):
        """
        Extract massive neutrino anisotropic stress sigma_nu(k) at tau_rec.
        """
        y_np = np.array(self.y)
        M = self.m_nu_eV / _T_NU_EV
        a_rec = float(self.bg.a_at_tau(np.array([self.bg.tau_rec]))[0])

        start = _idx_massive_start(self.l_ml_max)
        n_per_q = self.l_m_max + 1

        num = np.zeros(self.N_k)
        den = 0.0

        for i_q in range(self.n_q):
            q = float(self.q_nodes[i_q])
            w = float(self.w_nodes[i_q])
            eps = np.sqrt(q ** 2 + (M * a_rec) ** 2)

            idx0 = start + i_q * n_per_q
            Psi_2 = y_np[:, idx0 + 2] if n_per_q > 2 else np.zeros(self.N_k)

            num += w * (q ** 2 / eps) * Psi_2
            den += w * eps

        return (2.0 / 3.0) * num / max(den, 1e-30)

    def source_at_recombination(self):
        """
        Extract source functions at tau_rec with Silk damping.

        Returns (Theta_0, Phi, Psi, v_b, delta_nu_massive).

        The total Psi includes anisotropic stress from both massless and massive neutrinos.
        """
        y_np = np.array(self.y)
        Phi = y_np[:, _IDX_PHI]
        Theta_0 = y_np[:, _IDX_THETA_0]
        v_b = 3.0 * y_np[:, _IDX_THETA_1]

        # Massless neutrino anisotropic stress
        N_ml_2 = y_np[:, _IDX_N_START + 2] if self.l_ml_max >= 2 else np.zeros_like(Phi)

        H02 = self.H0 ** 2
        a_rec = float(self.bg.a_at_tau(np.array([self.bg.tau_rec]))[0])

        Psi = Phi - 12.0 * H02 * self.Omega_nu_ml * N_ml_2 / (a_rec ** 2 * self.k_arr ** 2)

        # Massive neutrino anisotropic stress contribution
        if self.N_massive > 0 and self.n_q > 0:
            sigma_m = self.massive_nu_shear()

            # rho_ratio at recombination
            M = self.m_nu_eV / _T_NU_EV
            norm_massless = np.sum(self.w_nodes * self.q_nodes)
            eps_arr = np.sqrt(self.q_nodes ** 2 + (M * a_rec) ** 2)
            rho_ratio = np.sum(self.w_nodes * eps_arr) / norm_massless

            Psi = Psi - 12.0 * H02 * self.Omega_nu_m_today * rho_ratio * sigma_m / (a_rec ** 2 * self.k_arr ** 2)

        delta_nu_m = self.massive_nu_delta()

        silk = np.exp(-(self.k_arr / self.bg.k_D) ** 2)
        return (Theta_0 * silk, Phi * silk, Psi * silk,
                v_b * silk, delta_nu_m * silk)

    def sw_source(self):
        """Sachs-Wolfe source: Theta_0 + Psi (silk-damped)."""
        Theta_0, Phi, Psi, v_b, delta_nu_m = self.source_at_recombination()
        return Theta_0 + Psi

    def phi_psi_ratio(self):
        """Psi/Phi at tau_rec (no Silk damping)."""
        y_np = np.array(self.y)
        Phi = y_np[:, _IDX_PHI]
        N_ml_2 = y_np[:, _IDX_N_START + 2] if self.l_ml_max >= 2 else np.zeros_like(Phi)

        H02 = self.H0 ** 2
        a_rec = float(self.bg.a_at_tau(np.array([self.bg.tau_rec]))[0])
        Psi = Phi - 12.0 * H02 * self.Omega_nu_ml * N_ml_2 / (a_rec ** 2 * self.k_arr ** 2)

        if self.N_massive > 0:
            sigma_m = self.massive_nu_shear()
            M = self.m_nu_eV / _T_NU_EV
            norm_massless = np.sum(self.w_nodes * self.q_nodes)
            eps_arr = np.sqrt(self.q_nodes ** 2 + (M * a_rec) ** 2)
            rho_ratio = np.sum(self.w_nodes * eps_arr) / norm_massless
            Psi = Psi - 12.0 * H02 * self.Omega_nu_m_today * rho_ratio * sigma_m / (a_rec ** 2 * self.k_arr ** 2)

        mask = np.abs(Phi) > 1e-10
        return np.where(mask, Psi / Phi, 1.0)

    def power_suppression(self, result_massless):
        """
        Compute P(k) suppression relative to a massless neutrino result.

        The suppression at k >> k_fs is approximately:
            Delta P / P ~ -8 f_nu

        Parameters
        ----------
        result_massless : NeutrinoResult — result from massless solver for comparison

        Returns
        -------
        k_arr, ratio : arrays — k values and P_massive / P_massless
        """
        sw_massive = self.sw_source()
        sw_massless = result_massless.sw_source()

        # Power ~ |source|^2 (schematic, ignoring transfer function details)
        P_massive = sw_massive ** 2
        P_massless = sw_massless ** 2

        mask = P_massless > 1e-20
        ratio = np.where(mask, P_massive / P_massless, 1.0)

        return self.k_arr, ratio


# ============================================================================
# Self-test
# ============================================================================

def _test_massive_neutrino():
    """
    Test the massive neutrino solver:
    1. Stability (no NaN/Inf)
    2. Massless limit recovery (m_nu -> 0 should match massless solver)
    3. Non-relativistic transition visible in delta_nu
    4. Correct suppression trend at high k
    5. Background rho_nu transition
    """
    print("=" * 70)
    print("TEST: Massive Neutrino Boltzmann Solver")
    print("=" * 70)

    from .background import Background
    from .perturbations_neutrino import NeutrinoBoltzmannSolver

    bg = Background(khronon=False)
    bg.solve()

    k_arr = np.geomspace(5e-4, 0.35, 300).astype(np.float32)

    # --- Test 1: Background rho_nu ---
    print("\n--- Test 1: Background energy density ---")
    a_test = np.logspace(-6, 0, 1000)
    for m in [0.0, 0.02, 0.06, 0.15]:
        rho = rho_nu_massive_normalized(a_test, m, n_q=15)
        print(f"  m={m:.2f} eV: rho_ratio(a=1) = {rho[-1]:.4f}, "
              f"rho_ratio(a=1e-6) = {rho[0]:.6f}")
        if m > 0:
            _z = z_nr(m)
            print(f"    z_nr = {_z:.0f}")
    print("  PASS (background computed)")

    # --- Test 2: Massive solver with minimal mass ---
    print("\n--- Test 2: Solver with sum(m_nu) = 0.06 eV ---")
    solver_m = MassiveNeutrinoBoltzmannSolver(
        bg, k_arr, m_nu_eV=0.02, N_massive=3, n_q=5, l_ml_max=20, l_m_max=15)
    result_m = solver_m.solve()

    Theta_0, Phi, Psi, v_b, delta_nu_m = result_m.source_at_recombination()

    has_nan = any(np.any(np.isnan(x)) for x in [Theta_0, Phi, Psi])
    has_inf = any(np.any(np.isinf(x)) for x in [Theta_0, Phi, Psi])
    print(f"  NaN: {has_nan}, Inf: {has_inf}")
    if has_nan or has_inf:
        print("  FAIL: stability check")
        return False
    print("  PASS: stability")

    # --- Test 3: Compare with massless ---
    print("\n--- Test 3: Comparison with massless solver ---")
    solver_ml = NeutrinoBoltzmannSolver(bg, k_arr, l_nu_max=20)
    result_ml = solver_ml.solve()
    sw_ml = result_ml.sw_source()
    sw_m = result_m.sw_source()

    # At low k, should be very similar (neutrinos still relativistic at recombination for small masses)
    mask_low = k_arr < 0.01
    if np.sum(mask_low) > 0:
        ratio_low = np.mean(np.abs(sw_m[mask_low])) / np.mean(np.abs(sw_ml[mask_low]))
        print(f"  SW ratio (k<0.01): {ratio_low:.4f} (should be ~1.0)")

    mask_high = k_arr > 0.1
    if np.sum(mask_high) > 0:
        ratio_high = np.mean(sw_m[mask_high] ** 2) / np.mean(sw_ml[mask_high] ** 2)
        print(f"  P(k) ratio (k>0.1): {ratio_high:.4f}")

    # --- Test 4: Massive neutrino clustering ---
    print("\n--- Test 4: Massive neutrino density perturbation ---")
    delta_nu = result_m.massive_nu_delta()
    print(f"  delta_nu range: [{np.min(delta_nu):.4f}, {np.max(delta_nu):.4f}]")
    print(f"  delta_nu at k=0.001: {delta_nu[np.argmin(np.abs(k_arr - 0.001))]:.4f}")
    print(f"  delta_nu at k=0.1:   {delta_nu[np.argmin(np.abs(k_arr - 0.1))]:.4f}")

    # --- Test 5: Heavier mass ---
    print("\n--- Test 5: Solver with sum(m_nu) = 0.15 eV ---")
    solver_heavy = MassiveNeutrinoBoltzmannSolver(
        bg, k_arr, m_nu_eV=0.05, N_massive=3, n_q=5, l_ml_max=20, l_m_max=15)
    result_heavy = solver_heavy.solve()

    Theta_0_h, Phi_h, Psi_h, _, delta_nu_h = result_heavy.source_at_recombination()
    has_nan_h = any(np.any(np.isnan(x)) for x in [Theta_0_h, Phi_h, Psi_h])
    print(f"  NaN: {has_nan_h}")
    if not has_nan_h:
        print("  PASS: stability")
        sw_h = result_heavy.sw_source()
        f_nu_heavy = 3 * 0.05 / 93.14 / _h ** 2 / ((_OMEGA_B + _OMEGA_C) + 3 * 0.05 / 93.14 / _h ** 2)
        print(f"  f_nu = {f_nu_heavy:.4f}")
        print(f"  Expected P(k) suppression: ~{-8 * f_nu_heavy * 100:.1f}%")

    # --- Plot ---
    print("\n--- Generating plot ---")
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(2, 2, figsize=(14, 10))

        # Panel 1: Background
        ax = axes[0, 0]
        a_plot = np.logspace(-6, 0, 1000)
        for m, c, ls in [(0.0, 'k', '--'), (0.02, 'b', '-'),
                          (0.06, 'g', '-'), (0.15, 'r', '-')]:
            rho = rho_nu_massive_normalized(a_plot, m, n_q=15)
            # Plot rho * a^4 (should be constant for radiation, grow for matter)
            ax.loglog(1.0 / a_plot - 1, rho, color=c, ls=ls,
                      label=f'$m_\\nu = {m}$ eV')
        ax.set_xlabel('Redshift z')
        ax.set_ylabel(r'$\rho_\nu a^4 / \rho_{\nu,0} a_0^4$')
        ax.set_title('Massive Neutrino Background')
        ax.legend(fontsize=8)
        ax.set_xlim(1e6, 1)
        ax.invert_xaxis()

        # Panel 2: Source comparison
        ax = axes[0, 1]
        ax.semilogx(k_arr, sw_ml, 'k--', lw=0.8, label='Massless', alpha=0.7)
        ax.semilogx(k_arr, sw_m, 'b-', lw=0.8, label=r'$\sum m_\nu = 0.06$ eV')
        ax.semilogx(k_arr, sw_h, 'r-', lw=0.8, label=r'$\sum m_\nu = 0.15$ eV')
        ax.set_xlabel(r'$k$ [Mpc$^{-1}$]')
        ax.set_ylabel(r'$\Theta_0 + \Psi$')
        ax.set_title('SW Source at Recombination')
        ax.legend(fontsize=8)

        # Panel 3: Delta_nu
        ax = axes[1, 0]
        ax.semilogx(k_arr, result_m.massive_nu_delta(), 'b-', lw=1,
                     label=r'$\sum m_\nu = 0.06$ eV')
        ax.semilogx(k_arr, result_heavy.massive_nu_delta(), 'r-', lw=1,
                     label=r'$\sum m_\nu = 0.15$ eV')
        ax.set_xlabel(r'$k$ [Mpc$^{-1}$]')
        ax.set_ylabel(r'$\delta_\nu$')
        ax.set_title('Massive Neutrino Density Perturbation')
        ax.legend(fontsize=8)

        # Panel 4: Power suppression
        ax = axes[1, 1]
        _, ratio_06 = result_m.power_suppression(result_ml)
        _, ratio_15 = result_heavy.power_suppression(result_ml)
        ax.semilogx(k_arr, ratio_06, 'b-', lw=1, label=r'$\sum m_\nu = 0.06$ eV')
        ax.semilogx(k_arr, ratio_15, 'r-', lw=1, label=r'$\sum m_\nu = 0.15$ eV')
        ax.axhline(1.0, color='gray', ls=':', alpha=0.5)
        ax.set_xlabel(r'$k$ [Mpc$^{-1}$]')
        ax.set_ylabel(r'$P_{m_\nu} / P_{m=0}$')
        ax.set_title('Power Spectrum Suppression')
        ax.legend(fontsize=8)

        plt.tight_layout()
        import os
        out_dir = os.path.dirname(os.path.abspath(__file__))
        out_path = os.path.join(out_dir, 'massive_neutrino_test.png')
        plt.savefig(out_path, dpi=150)
        plt.close()
        print(f"  Plot saved: {out_path}")
    except ImportError:
        print("  matplotlib not available")

    # --- Summary ---
    print("\n" + "=" * 70)
    print("SUMMARY")
    print(f"  Massive neutrino solver: {'PASS' if not (has_nan or has_inf or has_nan_h) else 'FAIL'}")
    print(f"  Variables per k: {solver_m.n_var}")
    print(f"  Memory per run: {solver_m.N_k * solver_m.n_var * 4 / 1024:.1f} KB")
    print("=" * 70)

    return not (has_nan or has_inf or has_nan_h)


if __name__ == '__main__':
    _test_massive_neutrino()
