"""
halofit.py — Non-linear matter power spectrum via revised Halofit.

Implements the Takahashi et al. (2012) revision of Halofit (Smith et al. 2003)
to compute the non-linear P(k) from the linear P(k).

Algorithm:
  1. From linear P(k,z), compute sigma^2(R) = int dk k^2 P_L(k) W^2(kR) / (2 pi^2)
     where W(x) = 3(sin x - x cos x) / x^3 is the real-space top-hat window.
  2. Find R_nl (= R_sigma) where sigma(R_nl) = 1 (non-linear scale).
  3. k_sigma = 1/R_sigma,  n_eff = -2 d(ln sigma)/d(ln R) - 3,
     C = -d^2(ln sigma)/d(ln R)^2.
  4. Halofit decomposition:
       Delta^2_nl(k) = Delta^2_Q(k) + Delta^2_H(k)
     where Delta^2_Q is the quasi-linear (two-halo) term and Delta^2_H is the
     one-halo term.

References:
  Smith et al. (2003), MNRAS 341, 1311  [arXiv:astro-ph/0207664]
  Takahashi et al. (2012), ApJ 761, 152  [arXiv:1208.2190]

Input:  linear P(k) from matter_pk.py (or any array of k, P_L(k))
Output: non-linear P(k) at the same k values

Author: Sheng-Kai Huang, 2026
"""
import numpy as np
from scipy.optimize import brentq
from scipy.interpolate import interp1d

from .background import (
    Omega_b, Omega_m, Omega_r, Omega_L,
    H0_Mpc, h as h_param, A_s, n_s, k_pivot,
)
from .matter_pk import growth_factor_integral


# ============================================================================
# Omega_m(z) and Omega_DE(z) for the fitting formulae
# ============================================================================

def _Omega_m_z(z, Omega_m_val=None, Omega_L_val=None):
    """Omega_m(z) = Omega_m (1+z)^3 / E^2(z)."""
    if Omega_m_val is None:
        Omega_m_val = Omega_m
    if Omega_L_val is None:
        Omega_L_val = Omega_L
    a = 1.0 / (1.0 + z)
    E2 = Omega_r / a**4 + Omega_m_val / a**3 + Omega_L_val
    return Omega_m_val / a**3 / E2


def _Omega_DE_z(z, Omega_m_val=None, Omega_L_val=None):
    """Omega_DE(z) = Omega_Lambda / E^2(z). (w = -1 assumed.)"""
    if Omega_m_val is None:
        Omega_m_val = Omega_m
    if Omega_L_val is None:
        Omega_L_val = Omega_L
    a = 1.0 / (1.0 + z)
    E2 = Omega_r / a**4 + Omega_m_val / a**3 + Omega_L_val
    return Omega_L_val / E2


# ============================================================================
# Top-hat window and sigma^2(R) computation
# ============================================================================

def _tophat_window(x):
    """
    Fourier-space real-space top-hat window function:
      W(x) = 3 (sin x - x cos x) / x^3
    Handles x -> 0 limit.  Works for both scalar and array input.
    """
    scalar = np.ndim(x) == 0
    x = np.atleast_1d(np.asarray(x, dtype=np.float64))
    W = np.ones_like(x)
    mask = np.abs(x) > 1e-8
    xm = x[mask]
    W[mask] = 3.0 * (np.sin(xm) - xm * np.cos(xm)) / xm**3
    mask2 = ~mask & (np.abs(x) > 0)
    if np.any(mask2):
        x2 = x[mask2]
        W[mask2] = 1.0 - x2**2 / 10.0
    if scalar:
        return float(W[0])
    return W


def sigma_squared(R, k_arr, Pk_arr):
    """
    Compute sigma^2(R) = (1/2pi^2) int_0^inf dk k^2 P(k) W^2(kR).

    Uses trapezoidal integration in log(k) space.

    Parameters
    ----------
    R : float
        Smoothing radius (same units as 1/k, typically Mpc).
    k_arr : array
        Wavenumber grid (must be sorted, positive).
    Pk_arr : array
        Linear power spectrum P(k).

    Returns
    -------
    sigma2 : float
        Variance of the smoothed density field.
    """
    lnk = np.log(k_arr)
    x = k_arr * R
    W = _tophat_window(x)
    # Integrand: d(sigma^2)/d(ln k) = k^3 P(k) W^2(kR) / (2 pi^2)
    integrand = k_arr**3 * Pk_arr * W**2 / (2.0 * np.pi**2)
    # Trapezoidal in ln(k)
    return float(np.trapezoid(integrand, lnk))


# ============================================================================
# Generate linear P(k) using BBKS + Sugiyama shape parameter
# ============================================================================

def generate_linear_pk(k_Mpc, z=0.0, Omega_m_val=None, Omega_b_val=None,
                       Omega_L_val=None, h_val=None, sigma8_target=None):
    """
    Generate the linear matter power spectrum using the BBKS transfer function
    with the Sugiyama (1995) shape parameter correction.

    The full formula:
      P(k,z) = (2pi^2/k^3) * (4/25) * (k/H0)^4 / Omega_m^2
               * A_s * (k/k_pivot)^(n_s-1)
               * T^2(k) * D^2(z)

    This is used both for testing Halofit and for extending solver P(k)
    to higher k where the Boltzmann solver may not have data.

    Parameters
    ----------
    k_Mpc : array
        Wavenumber grid in Mpc^-1.
    z : float
        Output redshift.
    Omega_m_val, Omega_b_val, Omega_L_val : float or None
        Cosmological parameters (default: from background).
    h_val : float or None
        Hubble parameter (default: from background).
    sigma8_target : float or None
        If given, rescale P(k) to match this sigma_8 at z=0.

    Returns
    -------
    Pk : array
        Linear P(k) in Mpc^3.
    """
    if Omega_m_val is None:
        Omega_m_val = Omega_m
    if Omega_b_val is None:
        Omega_b_val = Omega_b
    if Omega_L_val is None:
        Omega_L_val = Omega_L
    if h_val is None:
        h_val = h_param

    # Sugiyama (1995) shape parameter
    Gamma = Omega_m_val * h_val * np.exp(
        -Omega_b_val * (1.0 + np.sqrt(2.0 * h_val) / Omega_m_val))

    # BBKS transfer function
    k_hMpc = k_Mpc / h_val
    q = k_hMpc / Gamma
    T = (np.log(1.0 + 2.34 * q) / (2.34 * q)
         * (1.0 + 3.89 * q + (16.1 * q)**2 + (5.46 * q)**3
            + (6.71 * q)**4)**(-0.25))

    # Growth factor
    a_out = 1.0 / (1.0 + z)
    D_out = growth_factor_integral(a_out, Omega_m_val, Omega_L_val)

    # Primordial power spectrum
    P_R = A_s * (k_Mpc / k_pivot)**(n_s - 1.0)

    # Full P(k): includes Poisson equation factor (4/25)(k/H0)^4/Omega_m^2
    H0 = h_val * 100.0 / 299792.458  # H0 in Mpc^-1
    fac = (4.0 / 25.0) * (k_Mpc / H0)**4 / Omega_m_val**2

    Pk = (2.0 * np.pi**2 / k_Mpc**3) * fac * P_R * T**2 * D_out**2

    # Optional sigma_8 rescaling
    if sigma8_target is not None:
        R8 = 8.0 / h_val
        s8_sq = sigma_squared(R8, k_Mpc, Pk)
        Pk *= sigma8_target**2 / max(s8_sq, 1e-200)

    return Pk


def _extend_pk(k_arr, Pk_arr, k_max_ext=200.0, N_ext=500):
    """
    Extend P(k) to higher k using the BBKS transfer function,
    normalised to match the input data at the junction.

    Essential because Halofit needs P(k) up to k ~ 100 Mpc^-1
    for the sigma(R) integration at the non-linear scale.

    Parameters
    ----------
    k_arr, Pk_arr : arrays
        Input P(k) in Mpc^-1 and Mpc^3.
    k_max_ext : float
        Extend to this maximum k (Mpc^-1).
    N_ext : int
        Number of extra points.

    Returns
    -------
    k_ext, Pk_ext : arrays
        Extended P(k) arrays.
    """
    k_last = k_arr[-1]
    if k_last >= k_max_ext:
        return k_arr, Pk_arr

    # Generate BBKS P(k) over the combined range
    k_extra = np.geomspace(k_last * 1.005, k_max_ext, N_ext)
    k_all = np.concatenate([k_arr, k_extra])

    # BBKS P(k) shape at the full range (z=0, not normalised)
    Pk_bbks_all = generate_linear_pk(k_all, z=0.0)

    # Also generate BBKS at the data k range for normalisation
    Pk_bbks_data = generate_linear_pk(k_arr, z=0.0)

    # Find normalisation by matching at mid-to-high k of the data
    # (avoid the very ends where numerical noise may be larger)
    N_match = min(100, len(k_arr) // 3)
    if N_match < 5:
        N_match = 5
    # Use the middle-to-end portion of the data where the shape is reliable
    i_start = max(0, len(k_arr) - 2 * N_match)
    i_end = len(k_arr) - N_match // 4  # avoid very last points
    if i_end <= i_start:
        i_start = 0
        i_end = len(k_arr)

    ratios = Pk_arr[i_start:i_end] / np.maximum(Pk_bbks_data[i_start:i_end], 1e-200)
    # Use median for robustness against outliers
    norm = np.median(ratios)

    # Apply normalisation to extended portion
    Pk_ext_part = Pk_bbks_all[len(k_arr):] * norm

    k_ext = np.concatenate([k_arr, k_extra])
    Pk_ext = np.concatenate([Pk_arr, Pk_ext_part])

    return k_ext, Pk_ext


# ============================================================================
# Find the non-linear scale R_sigma where sigma(R_sigma) = 1
# ============================================================================

def _find_R_sigma(k_arr, Pk_arr, R_min=0.05, R_max=500.0):
    """
    Find R_sigma such that sigma(R_sigma) = 1.

    Automatically extends P(k) to high k if needed.

    Parameters
    ----------
    k_arr, Pk_arr : arrays
        Linear P(k).
    R_min, R_max : float
        Bracket for root finding (Mpc).

    Returns
    -------
    R_sigma : float
        Non-linear scale.
    """
    # Extend P(k) to high k for accurate small-R integration
    k_ext, Pk_ext = _extend_pk(k_arr, Pk_arr, k_max_ext=200.0)

    def f(lnR):
        R = np.exp(lnR)
        return sigma_squared(R, k_ext, Pk_ext) - 1.0

    # Adaptive bracketing
    f_min = f(np.log(R_min))
    f_max = f(np.log(R_max))

    while f_min < 0 and R_min > 1e-4:
        R_min *= 0.1
        f_min = f(np.log(R_min))

    while f_max > 0 and R_max < 1e5:
        R_max *= 10.0
        f_max = f(np.log(R_max))

    if f_min * f_max > 0:
        raise RuntimeError(
            f"Cannot bracket R_sigma: "
            f"sigma^2(R={R_min:.6f})={f_min+1:.6f}, "
            f"sigma^2(R={R_max:.4f})={f_max+1:.6f}. "
            f"k_range=[{k_ext[0]:.2e}, {k_ext[-1]:.2e}]")

    lnR_sigma = brentq(f, np.log(R_min), np.log(R_max), xtol=1e-8, rtol=1e-8)
    return np.exp(lnR_sigma)


# ============================================================================
# Effective spectral index n_eff and curvature C
# ============================================================================

def _compute_neff_C(R_sigma, k_arr, Pk_arr, dlnR=0.01):
    """
    Compute n_eff and C at R_sigma via finite differences.

    n_eff = -d(ln sigma^2)/d(ln R)|_{R_sigma} - 3
    C     = -d^2(ln sigma^2)/d(ln R)^2|_{R_sigma}
    """
    k_ext, Pk_ext = _extend_pk(k_arr, Pk_arr, k_max_ext=200.0)

    lnR0 = np.log(R_sigma)
    R_vals = np.exp([lnR0 - dlnR, lnR0, lnR0 + dlnR])

    s2 = np.array([sigma_squared(R, k_ext, Pk_ext) for R in R_vals])
    lns2 = np.log(np.maximum(s2, 1e-200))

    dlns2_dlnR = (lns2[2] - lns2[0]) / (2.0 * dlnR)
    d2lns2_dlnR2 = (lns2[2] - 2.0 * lns2[1] + lns2[0]) / dlnR**2

    n_eff = -dlns2_dlnR - 3.0
    C = -d2lns2_dlnR2

    return n_eff, C


# ============================================================================
# Takahashi et al. (2012) Halofit fitting formulae
# ============================================================================

def _halofit_params(n_eff, C, Om_z, ODE_z):
    """
    Compute all Halofit fitting parameters from Takahashi et al. (2012).

    Equations (A5)-(A15) of arXiv:1208.2190.
    """
    n = n_eff
    Cv = C

    # Eq. (A5)
    a_n = 10.0**(
        1.5222 + 2.8553 * n + 2.3706 * n**2 + 0.9903 * n**3
        + 0.2250 * n**4 - 0.6038 * Cv
        + 0.1749 * ODE_z * (1.0 - 0.5562 * ODE_z)
    )

    # Eq. (A6)
    b_n = 10.0**(
        -0.5642 + 0.5864 * n + 0.5716 * n**2
        - 1.5474 * Cv
        + 0.2279 * ODE_z * (1.0 - 0.4236 * ODE_z)
    )

    # Eq. (A7)
    c_n = 10.0**(
        0.3698 + 2.0404 * n + 0.8161 * n**2
        + 0.5869 * Cv
    )

    # Eq. (A8)
    gamma_n = 0.1971 - 0.0843 * n + 0.8460 * Cv

    # Eq. (A9)
    alpha_n = np.abs(
        6.0835 + 1.3373 * n - 0.1959 * n**2 - 5.5274 * Cv
    )

    # Eq. (A10)
    beta_n = (
        2.0379 - 0.7354 * n + 0.3157 * n**2
        + 1.2490 * n**3 + 0.3980 * n**4 - 0.1682 * Cv
    )

    # Eq. (A11)
    mu_n = 0.0

    # Eq. (A12)
    nu_n = 10.0**(5.2105 + 3.6902 * n)

    # Eq. (A13)-(A15): f1, f2, f3 for w=-1 (LCDM)
    f1 = Om_z**(-0.0732)
    f2 = Om_z**(-0.1423)
    f3 = Om_z**(0.0725)

    return dict(
        a_n=a_n, b_n=b_n, c_n=c_n, gamma_n=gamma_n,
        alpha_n=alpha_n, beta_n=beta_n, mu_n=mu_n, nu_n=nu_n,
        f1=f1, f2=f2, f3=f3,
    )


def _halofit_Delta2(k, k_sigma, Delta2_L, params):
    """
    Compute non-linear dimensionless power spectrum Delta^2_nl(k).

    Takahashi et al. (2012) Eq. (A1)-(A4):
      Delta^2_nl = Delta^2_Q + Delta^2_H
    """
    a_n = params['a_n']
    b_n = params['b_n']
    c_n = params['c_n']
    gamma_n = params['gamma_n']
    alpha_n = params['alpha_n']
    beta_n = params['beta_n']
    mu_n = params['mu_n']
    nu_n = params['nu_n']
    f1 = params['f1']
    f2 = params['f2']
    f3 = params['f3']

    y = k / k_sigma

    # Quasi-linear (two-halo) term — Eq. (A2)
    Delta2_Q = Delta2_L * (
        (1.0 + Delta2_L)**beta_n
        / (1.0 + alpha_n * Delta2_L)
    ) * np.exp(-y / 4.0 - y**2 / 8.0)

    # One-halo term — Eq. (A3)
    Delta2_H = (
        a_n * y**(3.0 * f1)
        / (1.0 + b_n * y**f2 + (c_n * f3 * y)**(3.0 - gamma_n))
    )

    # mu_n / nu_n damping — Eq. (A4)
    Delta2_H *= 1.0 / (1.0 + mu_n / y + nu_n / y**2)

    return Delta2_Q + Delta2_H


# ============================================================================
# Main interface: compute non-linear P(k)
# ============================================================================

def compute_nonlinear_pk(k_arr, Pk_linear, z=0.0,
                         Omega_m_val=None, Omega_L_val=None,
                         verbose=True):
    """
    Compute the non-linear matter power spectrum P_nl(k) from P_L(k)
    using the Takahashi et al. (2012) revised Halofit.

    Parameters
    ----------
    k_arr : array
        Wavenumber in Mpc^-1 (code units).
    Pk_linear : array
        Linear P_L(k) in Mpc^3.  Must already be at the correct redshift z.
    z : float
        Redshift.
    Omega_m_val : float or None
        Matter density parameter (default: from background).
    Omega_L_val : float or None
        Dark energy density parameter (default: from background).
    verbose : bool
        Print diagnostics.

    Returns
    -------
    Pk_nonlinear : array
        Non-linear P(k) in Mpc^3, same shape as k_arr.
    """
    if Omega_m_val is None:
        Omega_m_val = Omega_m
    if Omega_L_val is None:
        Omega_L_val = Omega_L

    # Step 1: Find R_sigma
    R_sigma = _find_R_sigma(k_arr, Pk_linear)
    k_sigma = 1.0 / R_sigma

    if verbose:
        print(f"[Halofit] z = {z:.2f}")
        print(f"[Halofit] R_sigma = {R_sigma:.4f} Mpc  =>  "
              f"k_sigma = {k_sigma:.4f} Mpc^-1  "
              f"({k_sigma/h_param:.4f} h/Mpc)")

    # Step 2: n_eff and C
    n_eff, C = _compute_neff_C(R_sigma, k_arr, Pk_linear)

    if verbose:
        print(f"[Halofit] n_eff = {n_eff:.4f},  C = {C:.4f}")

    # Step 3: Omega_m(z) and Omega_DE(z)
    Om_z = _Omega_m_z(z, Omega_m_val, Omega_L_val)
    ODE_z = _Omega_DE_z(z, Omega_m_val, Omega_L_val)

    if verbose:
        print(f"[Halofit] Omega_m(z={z:.1f}) = {Om_z:.4f},  "
              f"Omega_DE(z={z:.1f}) = {ODE_z:.4f}")

    # Step 4: Fitting parameters
    params = _halofit_params(n_eff, C, Om_z, ODE_z)

    if verbose:
        print(f"[Halofit] a_n={params['a_n']:.4f}, b_n={params['b_n']:.4f}, "
              f"c_n={params['c_n']:.4f}, gamma_n={params['gamma_n']:.4f}")
        print(f"[Halofit] alpha_n={params['alpha_n']:.4f}, "
              f"beta_n={params['beta_n']:.4f}, "
              f"nu_n={params['nu_n']:.2e}")

    # Step 5: Dimensionless linear power spectrum
    Delta2_L = k_arr**3 * Pk_linear / (2.0 * np.pi**2)

    # Step 6: Non-linear dimensionless power spectrum
    Delta2_nl = _halofit_Delta2(k_arr, k_sigma, Delta2_L, params)

    # Step 7: Convert back to P(k)
    Pk_nonlinear = Delta2_nl * (2.0 * np.pi**2) / k_arr**3

    if verbose:
        for k_check in [0.1, 0.3, 1.0, 5.0]:
            idx = np.searchsorted(k_arr, k_check)
            if idx < len(k_arr):
                boost = Pk_nonlinear[idx] / max(Pk_linear[idx], 1e-200)
                print(f"[Halofit] Boost at k={k_arr[idx]:.3f} Mpc^-1: "
                      f"P_nl/P_L = {boost:.2f}")

    return Pk_nonlinear


def compute_nonlinear_pk_hMpc(k_hMpc, Pk_hMpc3, z=0.0,
                               Omega_m_val=None, Omega_L_val=None,
                               verbose=True):
    """
    Convenience wrapper: input/output in h/Mpc and (Mpc/h)^3 units.
    """
    k_Mpc = k_hMpc * h_param
    Pk_Mpc3 = Pk_hMpc3 / h_param**3

    Pk_nl_Mpc3 = compute_nonlinear_pk(k_Mpc, Pk_Mpc3, z=z,
                                       Omega_m_val=Omega_m_val,
                                       Omega_L_val=Omega_L_val,
                                       verbose=verbose)

    return Pk_nl_Mpc3 * h_param**3


# ============================================================================
# Full pipeline: Boltzmann solver result -> non-linear P(k)
# ============================================================================

def nonlinear_pk_from_result(result, solver_type='implicit', z_out=0.0,
                              k_fine=None, N_k_fine=2000, k_max=50.0,
                              verbose=True):
    """
    End-to-end: Boltzmann solver result -> linear P(k) -> non-linear P(k).

    Parameters
    ----------
    result : solver result
        Output from any mlx_class Boltzmann solver.
    solver_type : str
        'implicit', 'neutrino', 'hires', 'explicit'.
    z_out : float
        Output redshift.
    k_fine : array or None
        Fine k-grid (Mpc^-1). If None, uses geomspace up to k_max.
    N_k_fine : int
        Number of k points if k_fine is None.
    k_max : float
        Maximum k in Mpc^-1 for the fine grid.
    verbose : bool
        Print diagnostics.

    Returns
    -------
    k_hMpc : array
        k in h/Mpc.
    Pk_L_hMpc3 : array
        Linear P(k) in (Mpc/h)^3.
    Pk_NL_hMpc3 : array
        Non-linear P(k) in (Mpc/h)^3.
    k_Mpc : array
        k in Mpc^-1.
    Pk_L_Mpc3 : array
        Linear P(k) in Mpc^3.
    Pk_NL_Mpc3 : array
        Non-linear P(k) in Mpc^3.
    """
    from .matter_pk import compute_matter_pk

    if k_fine is None:
        k_min = result.k_arr[0]
        k_fine = np.geomspace(k_min, k_max, N_k_fine).astype(np.float64)

    k_h, Pk_L_h, k_Mpc, Pk_L_Mpc = compute_matter_pk(
        result, solver_type=solver_type, z_out=z_out,
        k_fine=k_fine, N_k_fine=N_k_fine)

    bg = result.bg
    Omega_m_eff = Omega_b + bg.Omega_cdm
    Omega_L_eff = 1.0 - Omega_m_eff - Omega_r

    Pk_NL_Mpc = compute_nonlinear_pk(
        k_Mpc, Pk_L_Mpc, z=z_out,
        Omega_m_val=Omega_m_eff, Omega_L_val=Omega_L_eff,
        verbose=verbose)

    Pk_NL_h = Pk_NL_Mpc * h_param**3

    return k_h, Pk_L_h, Pk_NL_h, k_Mpc, Pk_L_Mpc, Pk_NL_Mpc


# ============================================================================
# Plotting
# ============================================================================

def plot_nonlinear_pk(k_hMpc, Pk_L, Pk_NL, z_out=0.0,
                      k_ref=None, Pk_ref_NL=None, ref_label='CLASS (Halofit)',
                      out_path=None):
    """
    Plot linear vs non-linear P(k), with optional reference comparison.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import os

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Left: P(k)
    ax = axes[0]
    ax.loglog(k_hMpc, Pk_L, 'b--', lw=1.2, alpha=0.7, label='Linear')
    ax.loglog(k_hMpc, Pk_NL, 'r-', lw=1.5, label='Halofit (non-linear)')
    if k_ref is not None and Pk_ref_NL is not None:
        ax.loglog(k_ref, Pk_ref_NL, 'k:', lw=1.2, label=ref_label)
    ax.set_xlabel(r'$k$ [$h$ Mpc$^{-1}$]', fontsize=13)
    ax.set_ylabel(r'$P(k)$ [(Mpc/$h$)$^3$]', fontsize=13)
    ax.set_title(f'Matter Power Spectrum (z={z_out:.0f})', fontsize=14)
    ax.legend(fontsize=11)
    ax.set_xlim(k_hMpc[0], k_hMpc[-1])

    # Right: Boost factor
    ax = axes[1]
    boost = Pk_NL / np.maximum(Pk_L, 1e-200)
    ax.semilogx(k_hMpc, boost, 'r-', lw=1.5, label='Halofit boost')
    ax.axhline(1.0, color='gray', ls='--', lw=0.8)
    if k_ref is not None and Pk_ref_NL is not None:
        Pk_L_interp = interp1d(np.log(k_hMpc), np.log(np.maximum(Pk_L, 1e-200)),
                                kind='cubic', fill_value='extrapolate')
        Pk_L_at_ref = np.exp(Pk_L_interp(np.log(k_ref)))
        boost_ref = Pk_ref_NL / np.maximum(Pk_L_at_ref, 1e-200)
        ax.semilogx(k_ref, boost_ref, 'k:', lw=1.2, label=ref_label)
    ax.set_xlabel(r'$k$ [$h$ Mpc$^{-1}$]', fontsize=13)
    ax.set_ylabel(r'$P_{\rm NL}(k) / P_{\rm L}(k)$', fontsize=13)
    ax.set_title(f'Non-linear Boost (z={z_out:.0f})', fontsize=14)
    ax.legend(fontsize=11)
    ax.set_xlim(k_hMpc[0], k_hMpc[-1])

    plt.tight_layout()

    if out_path is None:
        out_dir = os.path.dirname(os.path.abspath(__file__))
        out_path = os.path.join(out_dir, f'halofit_pk_z{z_out:.0f}.png')

    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"[Halofit] Plot saved: {out_path}")
    return out_path


# ============================================================================
# Self-test
# ============================================================================

def _test_halofit():
    """
    Test the Halofit implementation at z=0 and z=1.

    Uses BBKS linear P(k) normalised to sigma_8 = 0.811 (Planck 2018).
    This is independent of the Boltzmann solver.
    """
    import os
    print("=" * 70)
    print("TEST: Halofit Non-Linear Power Spectrum")
    print("=" * 70)

    # --- Generate BBKS linear P(k) at z=0 ---
    k_Mpc = np.geomspace(1e-4, 200.0, 10000)
    Pk_z0 = generate_linear_pk(k_Mpc, z=0.0, sigma8_target=0.811)

    k_hMpc = k_Mpc / h_param
    Pk_hMpc3 = Pk_z0 * h_param**3

    print(f"\n  Generated BBKS P(k, z=0): {len(k_Mpc)} points, "
          f"k=[{k_hMpc[0]:.2e}, {k_hMpc[-1]:.2e}] h/Mpc")

    # --- Test 1: sigma_8 ---
    print("\n--- Test 1: sigma(R=8/h Mpc) ---")
    R_8 = 8.0 / h_param
    sigma8 = np.sqrt(sigma_squared(R_8, k_Mpc, Pk_z0))
    print(f"  sigma_8 = {sigma8:.4f} (target: 0.811)")
    if abs(sigma8 - 0.811) < 0.01:
        print("  PASS")
    else:
        print("  WARNING: sigma_8 does not match target")

    # sigma at other scales
    for Rh in [12, 5, 3, 1]:
        s = np.sqrt(sigma_squared(Rh / h_param, k_Mpc, Pk_z0))
        print(f"  sigma(R={Rh} Mpc/h) = {s:.4f}")

    # --- Test 2: R_sigma, n_eff, C ---
    print("\n--- Test 2: Non-linear scale ---")
    R_sigma = _find_R_sigma(k_Mpc, Pk_z0)
    k_sigma = 1.0 / R_sigma
    print(f"  R_sigma = {R_sigma:.4f} Mpc = {R_sigma * h_param:.4f} Mpc/h")
    print(f"  k_sigma = {k_sigma:.4f} Mpc^-1 = {k_sigma / h_param:.4f} h/Mpc")

    n_eff, C = _compute_neff_C(R_sigma, k_Mpc, Pk_z0)
    print(f"  n_eff = {n_eff:.4f} (expect ~ -1.5 to -1.7)")
    print(f"  C = {C:.4f} (expect ~ -0.1 to 0.3)")

    R_sigma_h = R_sigma * h_param
    if 2.0 < R_sigma_h < 10.0:
        print("  PASS: R_sigma in expected range")
    else:
        print(f"  WARNING: R_sigma = {R_sigma_h:.4f} Mpc/h outside expected range")

    # --- Test 3: Non-linear P(k) at z=0 ---
    print("\n--- Test 3: Non-linear P(k) at z=0 ---")
    Pk_NL_z0 = compute_nonlinear_pk(k_Mpc, Pk_z0, z=0.0, verbose=True)
    Pk_NL_hMpc3 = Pk_NL_z0 * h_param**3

    # Check boost at representative scales
    for kh_check in [0.01, 0.05, 0.1, 0.3, 1.0, 5.0, 10.0]:
        idx = np.argmin(np.abs(k_hMpc - kh_check))
        boost = Pk_NL_hMpc3[idx] / max(Pk_hMpc3[idx], 1e-200)
        print(f"  k={k_hMpc[idx]:.4f} h/Mpc: P_NL/P_L = {boost:.4f}")

    # Linear regime: boost ~ 1 at k < 0.05
    mask_lin = k_hMpc < 0.05
    if np.sum(mask_lin) > 0:
        avg_boost = np.mean(Pk_NL_hMpc3[mask_lin] / Pk_hMpc3[mask_lin])
        print(f"  Average boost at k < 0.05 h/Mpc: {avg_boost:.4f} (expect ~1)")
        if abs(avg_boost - 1.0) < 0.10:
            print("  PASS: Linear regime preserved")
        else:
            print(f"  WARNING: Linear regime boost = {avg_boost:.4f}")

    # Non-linear regime: boost >> 1 at k > 1
    mask_nl = k_hMpc > 1.0
    if np.sum(mask_nl) > 0:
        avg_boost_nl = np.mean(Pk_NL_hMpc3[mask_nl] / np.maximum(Pk_hMpc3[mask_nl], 1e-200))
        print(f"  Average boost at k > 1 h/Mpc: {avg_boost_nl:.1f} (expect >> 1)")
        if avg_boost_nl > 2.0:
            print("  PASS: Non-linear boost present")
        else:
            print("  WARNING: Non-linear boost too small")

    # --- Test 4: z=1 ---
    print("\n--- Test 4: Non-linear P(k) at z=1 ---")
    D_0 = growth_factor_integral(1.0, Omega_m, Omega_L)
    D_1 = growth_factor_integral(0.5, Omega_m, Omega_L)  # a=0.5 => z=1
    growth_ratio = D_1 / D_0
    print(f"  D(z=1)/D(z=0) = {growth_ratio:.4f}")

    Pk_z1 = Pk_z0 * growth_ratio**2
    Pk_NL_z1 = compute_nonlinear_pk(k_Mpc, Pk_z1, z=1.0, verbose=True)

    # At z=1, non-linear boost should be smaller (less evolved)
    idx_1 = np.argmin(np.abs(k_hMpc - 1.0))
    boost_z0 = Pk_NL_z0[idx_1] / max(Pk_z0[idx_1], 1e-200)
    boost_z1 = Pk_NL_z1[idx_1] / max(Pk_z1[idx_1], 1e-200)
    print(f"  Boost at k=1 h/Mpc: z=0 -> {boost_z0:.2f}, z=1 -> {boost_z1:.2f}")
    if boost_z1 < boost_z0:
        print("  PASS: Boost decreases with redshift (as expected)")
    else:
        print("  WARNING: Boost does not decrease with redshift")

    # --- Test 5: Stability ---
    print("\n--- Test 5: Stability ---")
    has_nan = np.any(np.isnan(Pk_NL_hMpc3))
    has_inf = np.any(np.isinf(Pk_NL_hMpc3))
    has_neg = np.any(Pk_NL_hMpc3 < 0)
    print(f"  NaN: {has_nan}, Inf: {has_inf}, Negative: {has_neg}")
    if not (has_nan or has_inf or has_neg):
        print("  PASS: All P_NL values are positive and finite")
    else:
        print("  FAIL!")

    # --- Test 6: Compare with CLASS reference if available ---
    data_dir = os.path.dirname(os.path.abspath(__file__))
    pk_file_z0 = os.path.join(data_dir, 'pk_lcdm_z0.dat')
    if os.path.exists(pk_file_z0):
        print("\n--- Test 6: Compare with CLASS reference P(k) ---")
        data_class = np.loadtxt(pk_file_z0)
        k_class = data_class[:, 0]  # h/Mpc
        Pk_class = data_class[:, 1]  # (Mpc/h)^3

        # Interpolate our BBKS P(k) to CLASS k grid
        log_Pk_interp = interp1d(np.log(k_hMpc), np.log(np.maximum(Pk_hMpc3, 1e-200)),
                                  kind='cubic', fill_value='extrapolate')
        Pk_bbks_at_class = np.exp(log_Pk_interp(np.log(k_class)))

        # Compare shape (ratio) at linear scales
        mask_compare = (k_class > 0.01) & (k_class < 0.1)
        if np.sum(mask_compare) > 5:
            ratios = Pk_class[mask_compare] / np.maximum(Pk_bbks_at_class[mask_compare], 1e-200)
            print(f"  CLASS/BBKS ratio at 0.01 < k < 0.1: "
                  f"mean={np.mean(ratios):.3f}, std={np.std(ratios):.3f}")
            print("  (Ratio should be ~1 if shapes match; normalisation may differ)")

    # --- Plots ---
    print("\n--- Generating plots ---")
    plot_nonlinear_pk(k_hMpc, Pk_hMpc3, Pk_NL_hMpc3, z_out=0.0)

    Pk_z1_h = Pk_z1 * h_param**3
    Pk_NL_z1_h = Pk_NL_z1 * h_param**3
    plot_nonlinear_pk(k_hMpc, Pk_z1_h, Pk_NL_z1_h, z_out=1.0)

    print("\n" + "=" * 70)
    print("HALOFIT TEST COMPLETE")
    print("=" * 70)

    return k_hMpc, Pk_hMpc3, Pk_NL_hMpc3


if __name__ == '__main__':
    _test_halofit()
