"""
boltzmann.py — Photon Boltzmann hierarchy and Thomson scattering.

Implements the multipole expansion of the photon distribution:
    Theta_l(k, tau)  for l = 0, 1, 2, ..., l_max

Hierarchy (conformal Newtonian gauge):
    Theta_0' = -k Theta_1 - Phi'
    Theta_1' = (k/3)(Theta_0 + Psi) - kappa_dot (Theta_1 - v_b/3)
    Theta_l' = k/(2l+1) [l Theta_{l-1} - (l+1) Theta_{l+1}]
               - kappa_dot Theta_l    (l >= 2, except l=2 has polarization)

For l = 2: includes anisotropic stress (Pi = Theta_2) which couples to Phi - Psi.

Tight-coupling approximation (TCA) for tau << tau_rec:
    When |kappa_dot| >> k, expand in (k/kappa_dot) to avoid stiff ODE.
"""
import numpy as np


# ---------------------------------------------------------------------------
# Truncation: l_max and the closure relation
# ---------------------------------------------------------------------------
L_MAX_DEFAULT = 12   # Boltzmann hierarchy truncation
# For high accuracy one needs l_max ~ 25+, but for the simplified solver
# l_max = 12 gives reasonable results for l < 1500 in C_l.


def n_photon_vars(l_max):
    """Number of photon multipole variables: Theta_0 through Theta_{l_max}."""
    return l_max + 1


# ---------------------------------------------------------------------------
# Tight-coupling check
# ---------------------------------------------------------------------------
def is_tight_coupling(kappa_dot, k, threshold=10.0):
    """
    Check whether tight-coupling approximation should be used.
    Criterion: |kappa_dot| / k > threshold
    """
    if k == 0:
        return True
    return abs(kappa_dot) / k > threshold


# ---------------------------------------------------------------------------
# Photon hierarchy RHS
# ---------------------------------------------------------------------------
def photon_hierarchy_rhs(Theta, k, Phi_dot, Psi, kappa_dot, v_b, l_max):
    """
    Compute time derivatives of photon multipoles Theta_0 ... Theta_{l_max}.

    Parameters
    ----------
    Theta : array of length (l_max+1)
        Photon multipoles.
    k : float
        Wavenumber.
    Phi_dot : float
        Time derivative of metric perturbation Phi.
    Psi : float
        Metric perturbation Psi.
    kappa_dot : float
        Thomson scattering rate d(kappa)/d(tau) (negative).
    v_b : float
        Baryon velocity (divergence).
    l_max : int
        Maximum multipole.

    Returns
    -------
    dTheta : array of length (l_max+1)
    """
    dTheta = np.zeros(l_max + 1)

    # l = 0: monopole
    dTheta[0] = -k * Theta[1] - Phi_dot

    # l = 1: dipole
    dTheta[1] = (k / 3.0) * (Theta[0] + Psi) - kappa_dot * (Theta[1] - v_b / 3.0)

    # l = 2: quadrupole (includes polarization source Pi/10, simplified here)
    if l_max >= 2:
        Theta_3 = Theta[3] if l_max >= 3 else 0.0
        dTheta[2] = (k / 5.0) * (2 * Theta[1] - 3 * Theta_3) \
                     - kappa_dot * (Theta[2] - Theta[2] / 10.0)
        # The -kappa_dot * (9/10) Theta_2 term: scattering damps quadrupole
        # Full: -kappa_dot [Theta_2 - Pi/10]  where Pi = Theta_2 + Theta_P0 + Theta_P2
        # Simplified (no polarization): -kappa_dot * 9/10 * Theta_2
        dTheta[2] = (k / 5.0) * (2 * Theta[1] - 3 * Theta_3) \
                     - kappa_dot * 0.9 * Theta[2]

    # l >= 3: higher multipoles
    for l in range(3, l_max):
        dTheta[l] = k / (2*l + 1) * (l * Theta[l-1] - (l+1) * Theta[l+1]) \
                     - kappa_dot * Theta[l]

    # l = l_max: truncation (use closure relation)
    if l_max >= 3:
        # Closure: Theta_{l_max+1} ~ (2*l_max+1)/k/tau * Theta_{l_max} - Theta_{l_max-1}
        # Simplified: just set Theta_{l_max+1} = 0 (absorbing boundary)
        # Better: use l/(k*tau) approximation but need tau — use simple closure
        dTheta[l_max] = k * Theta[l_max - 1] - (l_max + 1) * Theta[l_max] * 0  # simplified
        dTheta[l_max] = k / (2*l_max + 1) * (l_max * Theta[l_max-1]) \
                         - kappa_dot * Theta[l_max]
        # Drop Theta_{l_max+1} term (= 0 truncation)

    return dTheta


# ---------------------------------------------------------------------------
# Tight-coupling equations
# ---------------------------------------------------------------------------
def tight_coupling_rhs(Theta0, Theta1, k, Phi_dot, Psi, calH, R, kappa_dot):
    """
    Tight-coupling approximation: only evolve Theta_0 and Theta_1.
    Higher multipoles are suppressed by powers of k/|kappa_dot|.

    In TCA:
        Theta_0' = -k Theta_1 - Phi'
        Theta_1' = [k(Theta_0 + Psi)/3 - calH * R * Theta_1 / (1+R)] / (1 + R/(1+R))
                  (combined baryon-photon fluid, slip = 0 to leading order)

    Parameters
    ----------
    R : float
        R = 3 rho_b / (4 rho_gamma) = R_factor * a

    Returns
    -------
    dTheta0, dTheta1 : floats
    """
    dTheta0 = -k * Theta1 - Phi_dot

    # Combined momentum equation (photon-baryon slip = 0):
    # (1+R) Theta_1' = k/3 (Theta_0 + Psi) - R calH Theta_1
    dTheta1 = (k / 3.0 * (Theta0 + Psi) - calH * R * Theta1) / (1.0 + R)

    return dTheta0, dTheta1


# ---------------------------------------------------------------------------
# Transfer function extraction
# ---------------------------------------------------------------------------
def transfer_Theta_l(Theta_grid, l):
    """
    Extract Theta_l(k, tau_0) from the integrated hierarchy.

    For C_l computation we need Theta_l at today, but more precisely
    we want the *source-integrated* transfer:
        Theta_l(k) = integral_0^{tau_0} S(k,tau) j_l(k(tau_0-tau)) dtau

    For the simplified solver, we use the "instantaneous recombination"
    line-of-sight approach (see power_spectrum.py).
    """
    return Theta_grid[l]
