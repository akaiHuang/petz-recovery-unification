"""
gdm_class_ode_extraction.py
============================
Complete perturbation ODE system extracted from gdm_class_public/source/perturbations.c.
Translated into Python pseudocode with exact line references.

Source: /Users/akaihuangm1/Desktop/github/gdm_class_public/source/perturbations.c
Function: perturb_derivs() [lines 8658-9556]
Metric eqs: perturb_einstein() [lines 6414-6588]
Stress-energy sums: perturb_total_stress_energy() [lines 6688-7124]
GDM functions: cs2_gdm_of_a_and_k() [line 10442], cv2_gdm_of_a_and_k() [line 10494]
Opacity: thermodynamics.c [line 126]

Convention (CLASS internal units):
  - Conformal time tau [Mpc]
  - Comoving wavenumber k [1/Mpc]
  - Densities rho in CLASS units where H^2 = sum(rho_i) [i.e. 8*pi*G/3 absorbed]
  - Primes (') = d/dtau
  - a_prime_over_a = aH = calH (conformal Hubble)
  - R = (4/3) * rho_gamma / rho_b  (baryon loading factor, line 8788)
  - s_l[l] = sqrt(1 - K*l*(l-1)/k^2) curvature factors (= 1 for flat K=0)
  - s2_squared = 1 - 3K/k^2 (line 8809)
  - k2 = k*k

Two gauges are supported:
  Newtonian gauge: metric_continuity = -3*Phi', metric_euler = k^2*Psi, metric_shear = 0
  Synchronous gauge: metric_continuity = h'/2, metric_euler = 0, metric_shear = k^2*alpha

All equations below are in a gauge-covariant form using the three "metric_" sources.
"""

import numpy as np


# ===========================================================================
# 0. NOTATION AND COEFFICIENT DEFINITIONS
# ===========================================================================
"""
Background quantities (from pvecback[]):
  a              = scale factor
  calH           = a' / a = aH (conformal Hubble parameter)
  rho_g          = photon energy density  (CLASS units)
  rho_b          = baryon energy density
  rho_cdm        = CDM energy density
  rho_ur         = massless neutrino energy density
  rho_gdm        = GDM (generalized dark matter) energy density
  w_gdm          = GDM equation of state p_gdm/rho_gdm
  ca2_gdm        = GDM adiabatic sound speed squared = w - w'/(3*(1+w)*calH)
  H              = Hubble rate H(a)

Thermodynamic quantities (from pvecthermo[]):
  dkappa = kappa_dot = d(kappa)/d(tau)   [Thomson scattering rate, 1/Mpc]
         = (1+z)^2 * n_e_today * x_e * sigma_T * (Mpc/m)
         (line 126 of thermodynamics.c)
         where n_e_today = pth->n_e, x_e = ionization fraction, sigma_T = Thomson cross section

  cb2    = baryon adiabatic sound speed squared = c_b^2

  g      = visibility function = dkappa * exp(-kappa)
  exp_m_kappa = exp(-kappa)

User-defined functions:
  cs2_gdm = cs2_gdm_of_a_and_k(a, k)  [sound speed in GDM rest frame]
  cv2_gdm = cv2_gdm_of_a_and_k(a, k)  [viscosity parameter for GDM]

Perturbation variables (y[]):
  delta_g   = photon density contrast (= 4 * Theta_0 in Dodelson convention)
  theta_g   = photon velocity divergence (= k * F_1 * 3/4 ... see below)
  shear_g   = F_2/2 = photon quadrupole (sigma_gamma)
  pol0_g .. polN_g = photon E-mode polarization multipoles G_l
  delta_b, theta_b = baryon density and velocity
  delta_cdm, theta_cdm = CDM density and velocity
  delta_gdm, theta_gdm, shear_gdm = GDM density, velocity, shear
  delta_ur, theta_ur, shear_ur, l3_ur, ... = massless neutrino hierarchy (N_l)
  phi (Newtonian gauge) or eta (synchronous gauge) = metric variable evolved as ODE

NOTE on conventions:
  CLASS uses delta_g = F_0 (the monopole of the distribution function),
  theta_X = (rho+p) * divergence of velocity / (rho+p), but for photons
  theta_g is defined so that the standard Boltzmann hierarchy reads:
    delta_g' = -4/3 * (theta_g + metric_continuity)
    theta_g' = k^2 * (delta_g/4 - s2^2 * shear_g) + metric_euler + ...
  This is the Ma & Bertschinger (1995) convention where theta = ik^j v_j.
"""


# ===========================================================================
# 1. METRIC EQUATIONS (Newtonian gauge)
# ===========================================================================
# Source: perturb_einstein(), lines 6460-6494
"""
Step 1: Compute total stress-energy sums (see Section 8 below):
  delta_rho         = sum_i rho_i * delta_i
  rho_plus_p_theta  = sum_i (rho_i + p_i) * theta_i
  rho_plus_p_shear  = sum_i (rho_i + p_i) * sigma_i
  delta_p           = sum_i delta_p_i

Step 2: Einstein equations in Newtonian gauge
"""

def metric_newtonian(phi, a, k, calH, delta_rho, rho_plus_p_theta, rho_plus_p_shear):
    """
    Newtonian gauge metric equations.
    phi is evolved as a dynamical variable (more stable than constraint).

    Lines 6491, 6494:
    """
    a2 = a**2
    k2 = k**2

    # Psi from anisotropic stress (line 6491)
    # Psi = Phi - 4.5 * (a^2/k^2) * sum_i (rho_i + p_i) * sigma_i
    Psi = phi - 4.5 * (a2 / k2) * rho_plus_p_shear

    # Phi' from the (0,i) Einstein equation (line 6494)
    # Phi' = -calH * Psi + 1.5 * (a^2/k^2) * sum_i (rho_i + p_i) * theta_i
    phi_prime = -calH * Psi + 1.5 * (a2 / k2) * rho_plus_p_theta

    return Psi, phi_prime


def metric_synchronous(eta, a, k, calH, delta_rho, rho_plus_p_theta,
                       rho_plus_p_shear, delta_p, K=0.0):
    """
    Synchronous gauge metric equations.
    eta is evolved as a dynamical variable.

    Lines 6520, 6545, 6548-6551, 6554, 6583-6586:
    """
    a2 = a**2
    k2 = k**2
    s2_squared = 1.0 - 3.0 * K / k2

    # h' from (0,0) Einstein eq (line 6520)
    h_prime = (k2 * s2_squared * eta + 1.5 * a2 * delta_rho) / (0.5 * calH)

    # eta' from (0,i) Einstein eq (line 6545)
    eta_prime = (1.5 * a2 * rho_plus_p_theta
                 + 0.5 * K * h_prime) / k2 / s2_squared

    # h'' from pressure equation (line 6548-6551)
    h_prime_prime = (-2.0 * calH * h_prime
                     + 2.0 * k2 * s2_squared * eta
                     - 9.0 * a2 * delta_p)

    # alpha = (h' + 6*eta') / (2*k^2)  (line 6554)
    alpha = (h_prime + 6.0 * eta_prime) / 2.0 / k2

    # alpha' from shear equation (line 6583-6586)
    alpha_prime = (-2.0 * calH * alpha
                   + eta
                   - 4.5 * (a2 / k2) * rho_plus_p_shear)

    return h_prime, eta_prime, h_prime_prime, alpha, alpha_prime


def metric_source_terms(gauge, phi_prime=0.0, Psi=0.0,
                        h_prime=0.0, alpha=0.0, alpha_prime=0.0):
    """
    Compute the three universal metric source terms.
    Lines 8887-8903.
    """
    if gauge == 'newtonian':
        metric_continuity = -3.0 * phi_prime          # line 8898
        metric_euler = k2 * Psi                       # line 8899
        metric_shear = 0.0                            # line 8900
        metric_ufa_class = -6.0 * phi_prime           # line 8902
    elif gauge == 'synchronous':
        metric_continuity = h_prime / 2.0             # line 8889
        metric_euler = 0.0                            # line 8890
        metric_shear = k2 * alpha                     # line 8891
        metric_ufa_class = h_prime / 2.0              # line 8893
    return metric_continuity, metric_euler, metric_shear, metric_ufa_class


# ===========================================================================
# 2. PHOTON TEMPERATURE HIERARCHY (Theta_l)
# ===========================================================================
# Source: perturb_derivs(), lines 8921-9006
#
# CLASS stores: delta_g = F_0, theta_g ~ F_1, shear_g = F_2/2, then l3_g, l4_g, ...
# The index scheme is: y[index_pt_delta_g + l] for l >= 3 (line 8997)
"""
Photon polarization tensor P0 (line 8972):
  Pi = P0 = (G_0 + G_2 + 2*s_l[2]*F_2) / 8
  where G_l = pol_l_g are polarization multipoles, F_2 = shear_g
"""

def photon_hierarchy_rhs(delta_g, theta_g, shear_g, theta_l,
                         pol_l, k, dkappa,
                         metric_continuity, metric_euler, metric_shear,
                         s_l, l_max_g, l_max_pol_g, theta_b, cotKgen):
    """
    Full photon Boltzmann hierarchy.

    Parameters
    ----------
    delta_g : float, F_0
    theta_g : float, photon velocity divergence
    shear_g : float, F_2/2
    theta_l : array [F_3, F_4, ..., F_{l_max}]
    pol_l   : array [G_0, G_1, G_2, G_3, ..., G_{l_max_pol}]
    k       : wavenumber
    dkappa  : Thomson scattering rate d(kappa)/d(tau) [positive quantity]
    metric_continuity, metric_euler, metric_shear : gauge-dependent source terms
    s_l     : curvature factors s_l[l] (= 1 for flat)
    l_max_g : maximum photon temperature multipole
    l_max_pol_g : maximum polarization multipole
    theta_b : baryon velocity divergence
    cotKgen : generalized cotK for truncation

    Returns
    -------
    d_delta_g, d_theta_g, d_shear_g, d_theta_l, d_pol_l
    """
    k2 = k * k
    s2 = s_l[2]  # sqrt(1 - 2K/k^2), = 1 for flat
    s2_squared = s2 * s2

    # Polarization source P0 (line 8972)
    # Pi = (G_0 + G_2 + 2*s2*shear_g) / 8
    P0 = (pol_l[0] + pol_l[2] + 2.0 * s2 * shear_g) / 8.0

    # --- Temperature monopole (line 8925) ---
    d_delta_g = -4.0 / 3.0 * (theta_g + metric_continuity)

    # --- Temperature dipole (line 8976-8979) ---
    d_theta_g = (k2 * (delta_g / 4.0 - s2_squared * shear_g)
                 + metric_euler
                 + dkappa * (theta_b - theta_g))

    # --- Temperature quadrupole / shear (line 8982-8985) ---
    d_shear_g = 0.5 * (
        8.0 / 15.0 * (theta_g + metric_shear)
        - 3.0 / 5.0 * k * s_l[3] / s_l[2] * theta_l[0]  # theta_l[0] = F_3
        - dkappa * (2.0 * shear_g - 4.0 / 5.0 / s_l[2] * P0)
    )

    # --- Temperature l=3 (line 8989-8992) ---
    l = 3
    d_l3_g = (k / (2.0 * l + 1.0)
              * (l * s_l[l] * 2.0 * s_l[2] * shear_g
                 - (l + 1) * s_l[l + 1] * theta_l[1])  # theta_l[1] = F_4
              - dkappa * theta_l[0])

    # --- Temperature l=4..l_max-1 (line 8995-8999) ---
    d_theta_l = np.zeros(l_max_g - 2)  # for l=3..l_max
    d_theta_l[0] = d_l3_g
    # theta_l[j] corresponds to F_{j+3}, y[index_pt_delta_g + j + 3]
    for l in range(4, l_max_g):
        j = l - 3  # index into theta_l array
        d_theta_l[j] = (k / (2.0 * l + 1.0)
                        * (l * s_l[l] * theta_l[j - 1]
                           - (l + 1) * s_l[l + 1] * theta_l[j + 1])
                        - dkappa * theta_l[j])

    # --- Temperature l=l_max (truncation, line 9003-9006) ---
    l = l_max_g
    j = l - 3
    d_theta_l[j] = (k * (s_l[l] * theta_l[j - 1]
                         - (1.0 + l) * cotKgen * theta_l[j])
                    - dkappa * theta_l[j])

    # --- Polarization G_0 (line 9010-9012) ---
    d_pol = np.zeros(l_max_pol_g + 1)
    d_pol[0] = (-k * pol_l[1]
                - dkappa * (pol_l[0] - 4.0 * P0))

    # --- Polarization G_1 (line 9016-9018) ---
    d_pol[1] = (k / 3.0 * (pol_l[0] - 2.0 * s_l[2] * pol_l[2])
                - dkappa * pol_l[1])

    # --- Polarization G_2 (line 9022-9024) ---
    d_pol[2] = (k / 5.0 * (2.0 * s_l[2] * pol_l[1] - 3.0 * s_l[3] * pol_l[3])
                - dkappa * (pol_l[2] - 4.0 / 5.0 * P0))

    # --- Polarization G_l, l=3..l_max_pol-1 (line 9028-9031) ---
    for l in range(3, l_max_pol_g):
        d_pol[l] = (k / (2.0 * l + 1.0)
                    * (l * s_l[l] * pol_l[l - 1]
                       - (l + 1.0) * s_l[l + 1] * pol_l[l + 1])
                    - dkappa * pol_l[l])

    # --- Polarization l=l_max_pol (truncation, line 9035-9038) ---
    l = l_max_pol_g
    d_pol[l] = (k * (s_l[l] * pol_l[l - 1]
                     - (l + 1) * cotKgen * pol_l[l])
                - dkappa * pol_l[l])

    return d_delta_g, d_theta_g, d_shear_g, d_theta_l, d_pol


# ===========================================================================
# 3. PHOTON TIGHT-COUPLING APPROXIMATION (TCA)
# ===========================================================================
# Source: lines 8948-9053
"""
When TCA is on (ppw->approx[index_ap_tca] == tca_on):

  Photon shear sigma_g in Newtonian gauge (line 6785):
    shear_g = (16/45) * theta_g / dkappa

  Photon shear sigma_g in synchronous gauge (line 6569):
    shear_g = (16/45) / dkappa * (theta_g + k^2 * alpha)

  Baryon velocity in TCA (line 8956-8960):
    theta_b' = 1/(1+R) * [-calH*theta_b
                           + k^2*(delta_p_b/rho_b + R*(delta_g/4 - s2^2*shear_g_tca))
                           + R * slip]
               + metric_euler

  where slip = tca_slip is computed in perturb_tca_slip_and_shear() [line 8951]

  Photon velocity in TCA (line 9050-9052):
    theta_g' = -(theta_b' + calH*theta_b - k^2*delta_p_b/rho_b) / R
               + k^2*(delta_g/4 - s2^2*shear_g_tca)
               + (1+R)/R * metric_euler
"""


def photon_baryon_tca(delta_g, theta_g, theta_b, delta_p_b_over_rho_b,
                      k, calH, R, dkappa, metric_continuity, metric_euler,
                      shear_g_tca, tca_slip, s2_squared=1.0):
    """
    Tight-coupling approximation for photon-baryon system.
    Lines 8925, 8956-8960, 9050-9052.

    R = (4/3)*rho_gamma/rho_b  (line 8788)
    """
    k2 = k * k

    # Photon monopole (same as exact, line 8925)
    d_delta_g = -4.0 / 3.0 * (theta_g + metric_continuity)

    # Baryon velocity in TCA (line 8956-8960)
    d_theta_b = ((-calH * theta_b
                  + k2 * (delta_p_b_over_rho_b
                          + R * (delta_g / 4.0 - s2_squared * shear_g_tca))
                  + R * tca_slip) / (1.0 + R)
                 + metric_euler)

    # Photon velocity in TCA (line 9050-9052)
    d_theta_g = (-(d_theta_b + calH * theta_b - k2 * delta_p_b_over_rho_b) / R
                 + k2 * (0.25 * delta_g - s2_squared * shear_g_tca)
                 + (1.0 + R) / R * metric_euler)

    return d_delta_g, d_theta_g, d_theta_b


# ===========================================================================
# 4. BARYONS (delta_b, theta_b)
# ===========================================================================
# Source: lines 8929-8962

def baryon_rhs(delta_b, theta_b, theta_g, delta_g,
               k, calH, R, dkappa, cb2,
               metric_continuity, metric_euler,
               tca=False, shear_g_tca=0.0, tca_slip=0.0,
               delta_p_b_over_rho_b=None, s2_squared=1.0):
    """
    Baryon density and velocity equations.

    Parameters
    ----------
    cb2 : baryon adiabatic sound speed squared (from thermodynamics)
    dkappa : Thomson scattering rate
    R : (4/3)*rho_gamma/rho_b
    """
    k2 = k * k

    if delta_p_b_over_rho_b is None:
        delta_p_b_over_rho_b = cb2 * delta_b  # Ma & Bertschinger approx (line 8830)

    # Baryon density (line 8931)
    d_delta_b = -(theta_b + metric_continuity)

    if not tca:
        # Without TCA (line 8940-8944)
        d_theta_b = (-calH * theta_b
                     + metric_euler
                     + k2 * delta_p_b_over_rho_b
                     + R * dkappa * (theta_g - theta_b))
    else:
        # With TCA (line 8956-8960)
        d_theta_b = ((-calH * theta_b
                      + k2 * (delta_p_b_over_rho_b
                              + R * (delta_g / 4.0 - s2_squared * shear_g_tca))
                      + R * tca_slip) / (1.0 + R)
                     + metric_euler)

    return d_delta_b, d_theta_b


# ===========================================================================
# 5. CDM (delta_cdm, theta_cdm)
# ===========================================================================
# Source: lines 9056-9073

def cdm_rhs(delta_cdm, theta_cdm, k, calH,
            metric_continuity, metric_euler, gauge='newtonian'):
    """
    CDM perturbation equations.

    Newtonian gauge (lines 9063-9065):
      delta_cdm' = -(theta_cdm + metric_continuity)
      theta_cdm' = -calH * theta_cdm + metric_euler

    Synchronous gauge (line 9071):
      delta_cdm' = -metric_continuity   [theta_cdm = 0 by gauge choice]
    """
    if gauge == 'newtonian':
        d_delta_cdm = -(theta_cdm + metric_continuity)
        d_theta_cdm = -calH * theta_cdm + metric_euler
    elif gauge == 'synchronous':
        d_delta_cdm = -metric_continuity
        d_theta_cdm = 0.0
    return d_delta_cdm, d_theta_cdm


# ===========================================================================
# 6. GDM / KHRONON (delta_gdm, theta_gdm, [shear_gdm])
# ===========================================================================
# Source: lines 9173-9216, with auxiliary defs at 10442-10506
"""
GDM is the key new species in gdm_class_public.
It has three free functions of time (and optionally k):
  w_gdm   = equation of state  [from background]
  cs2_gdm = rest-frame sound speed squared [cs2_gdm_of_a_and_k, line 10442]
  cv2_gdm = viscosity speed squared        [cv2_gdm_of_a_and_k, line 10494]
  ca2_gdm = adiabatic sound speed          [from background]

Non-adiabatic pressure (line 9185):
  Pi_nad = (cs2 - ca2) * [delta_gdm + 3*calH*(1+w)*theta_gdm/k^2]

  If k2_Pinad_gdm == True (line 9187-9188):
    Pi_nad *= k^2   (makes it k-dependent)

Shear can be algebraic or dynamic (controlled by dynamic_shear_gdm flag).

Algebraic shear in Newtonian gauge (line 6473):
  sigma_gdm = 8*cv2 / (15*(1+w)*calH) * theta_gdm

Algebraic shear in synchronous gauge (line 6560):
  sigma_gdm = 8*cv2 / (15*(1+w)*calH) * (theta_gdm + k^2*alpha)
"""


def gdm_rhs(delta_gdm, theta_gdm, shear_gdm,
            k, calH, w_gdm, ca2_gdm, cs2_gdm, cv2_gdm,
            metric_continuity, metric_euler, metric_shear,
            k2_Pinad=False, dynamic_shear=True, s2_squared=1.0):
    """
    Generalized Dark Matter (GDM / Khronon) perturbation equations.
    Lines 9173-9216.

    Parameters
    ----------
    delta_gdm, theta_gdm, shear_gdm : perturbation variables
    w_gdm, ca2_gdm, cs2_gdm, cv2_gdm : EOS and sound speed parameters
    k2_Pinad : if True, Pi_nad is multiplied by k^2
    dynamic_shear : if True, shear_gdm is evolved dynamically
    """
    k2 = k * k

    # Non-adiabatic pressure (line 9185)
    pinad_gdm = ((cs2_gdm - ca2_gdm)
                 * (delta_gdm + 3.0 * calH * (1.0 + w_gdm) * theta_gdm / k2))

    if k2_Pinad:  # line 9187-9188
        pinad_gdm *= k2

    # GDM density (line 9192-9194)
    d_delta_gdm = (-(1.0 + w_gdm) * (theta_gdm + metric_continuity)
                   + 3.0 * calH * ((w_gdm - ca2_gdm) * delta_gdm - pinad_gdm))

    # GDM velocity (line 9204-9207)
    d_theta_gdm = (-(1.0 - 3.0 * ca2_gdm) * calH * theta_gdm
                   + k2 / (1.0 + w_gdm) * (ca2_gdm * delta_gdm + pinad_gdm)
                   + metric_euler
                   - s2_squared * k2 * shear_gdm)

    # GDM shear (dynamic case, line 9210-9213)
    d_shear_gdm = 0.0
    if dynamic_shear:
        d_shear_gdm = (-3.0 * calH * shear_gdm
                       + 8.0 / 3.0 * cv2_gdm / (1.0 + w_gdm)
                       * (theta_gdm + metric_shear))

    return d_delta_gdm, d_theta_gdm, d_shear_gdm


def gdm_algebraic_shear_newtonian(theta_gdm, cv2_gdm, w_gdm, calH):
    """
    Algebraic (non-dynamic) GDM shear in Newtonian gauge.
    Line 6473.
    """
    return 8.0 * cv2_gdm / (15.0 * (1.0 + w_gdm) * calH) * theta_gdm


def gdm_algebraic_shear_synchronous(theta_gdm, cv2_gdm, w_gdm, calH, k2, alpha):
    """
    Algebraic (non-dynamic) GDM shear in synchronous gauge.
    Line 6560.
    """
    return 8.0 * cv2_gdm / (15.0 * (1.0 + w_gdm) * calH) * (theta_gdm + alpha * k2)


def cs2_gdm_khronon_dbi(ca2_gdm, k, a, rho_gdm, w_gdm):
    """
    Khronon DBI k-dependent sound speed (Blanchet & Skordis 2025, eq.4.38).
    Line 10461-10473.

    cs2(t,k) = ca2 / (1 + ca2 * k^2 / (1.5 * rho_gdm * a^2 * (1+w)))

    At large k: cs2 -> 0   (CDM-like on small scales)
    At small k: cs2 -> ca2  (adiabatic on large scales)

    Note: In CLASS units, 4*pi*G*rho_phys = (3/2)*rho_class
    """
    denom_factor = 1.5 * rho_gdm * a**2 * (1.0 + w_gdm)
    if denom_factor > 0.0 and ca2_gdm >= 0.0:
        cs2 = ca2_gdm / (1.0 + ca2_gdm * k**2 / denom_factor)
    else:
        cs2 = 0.0
    return min(max(cs2, 0.0), 1.0)


# ===========================================================================
# 7. MASSLESS NEUTRINOS / ULTRA-RELATIVISTIC RELICS (ur)
# ===========================================================================
# Source: lines 9322-9402
"""
Massless neutrino hierarchy: delta_ur, theta_ur, shear_ur, l3_ur, ...
Identical structure to photons but WITHOUT scattering (dkappa=0).

CLASS supports a generalized effective sound speed via three_ceff2_ur and
three_cvis2_ur (= 3*c_{eff}^2 and 3*c_{vis}^2), normally both = 1.

The "ultra-fluid approximation" (ufa) replaces the full hierarchy with
a truncated 3-variable system when safe to do so.
"""


def ur_hierarchy_rhs(delta_ur, theta_ur, shear_ur, ur_l, k, calH,
                     metric_continuity, metric_euler, metric_shear,
                     s_l, l_max_ur, cotKgen,
                     three_ceff2_ur=1.0, three_cvis2_ur=1.0,
                     metric_ufa_class=0.0, tau=0.0,
                     ufa='off', ufa_type='CLASS'):
    """
    Massless neutrino (ur) Boltzmann hierarchy.
    Lines 9331-9402.

    Parameters
    ----------
    delta_ur, theta_ur, shear_ur : monopole, dipole, quadrupole
    ur_l : array [N_3, N_4, ..., N_{l_max}]  (higher multipoles)
    three_ceff2_ur : 3 * c_{eff}^2  (default 1.0 for standard neutrinos)
    three_cvis2_ur : 3 * c_{vis}^2  (default 1.0)
    ufa : 'off' for exact hierarchy, 'on' for fluid approximation
    ufa_type : 'mb' (Ma & Bertschinger), 'hu', or 'CLASS'
    """
    k2 = k * k

    # --- ur density (line 9331-9335) ---
    d_delta_ur = (-4.0 / 3.0 * (theta_ur + metric_continuity)
                  + (1.0 - three_ceff2_ur) * calH
                  * (delta_ur + 4.0 * calH * theta_ur / k2))

    # --- ur velocity (line 9338-9342) ---
    d_theta_ur = (k2 * (three_ceff2_ur * delta_ur / 4.0
                        - s_l[2]**2 * shear_ur)
                  + metric_euler
                  - (1.0 - three_ceff2_ur) * calH * theta_ur)

    if ufa == 'off':
        # --- Exact ur shear (line 9347-9352) ---
        d_shear_ur = 0.5 * (
            # standard term, possibly modified by cvis2
            three_cvis2_ur * 8.0 / 15.0 * (theta_ur + metric_shear)
            - 3.0 / 5.0 * k * s_l[3] / s_l[2] * ur_l[0]  # ur_l[0] = N_3
        )
        # NOTE: the actual code writes this equivalently as:
        # d_shear = 0.5*(8./15.*(theta+ms) - 3./5.*k*s3/s2*N3
        #               - (1-3cvis2)*(8./15.*(theta+ms)))
        # which simplifies to the above when collected.

        # --- l=3 (line 9355-9357) ---
        d_ur_l = np.zeros(l_max_ur - 2)
        l = 3
        d_ur_l[0] = (k / (2.0 * l + 1.0)
                     * (l * 2.0 * s_l[l] * s_l[2] * shear_ur
                        - (l + 1) * s_l[l + 1] * ur_l[1]))

        # --- l=4..l_max-1 (line 9360-9362) ---
        for l in range(4, l_max_ur):
            j = l - 3
            d_ur_l[j] = (k / (2.0 * l + 1.0)
                         * (l * s_l[l] * ur_l[j - 1]
                            - (l + 1) * s_l[l + 1] * ur_l[j + 1]))

        # --- l=l_max (truncation, line 9366-9368) ---
        l = l_max_ur
        j = l - 3
        d_ur_l[j] = k * (s_l[l] * ur_l[j - 1]
                         - (1.0 + l) * cotKgen * ur_l[j])

    else:
        # Ultra-fluid approximation (ufa): only evolve shear (line 9374-9401)
        d_ur_l = np.zeros(0)

        if ufa_type == 'mb':  # Ma & Bertschinger (line 9377-9381)
            d_shear_ur = (-3.0 / tau * shear_ur
                          + 2.0 / 3.0 * (theta_ur + metric_shear))

        elif ufa_type == 'hu':  # Hu (line 9386-9390)
            d_shear_ur = (-3.0 * calH * shear_ur
                          + 2.0 / 3.0 * (theta_ur + metric_shear))

        elif ufa_type == 'CLASS':  # CLASS (line 9395-9399)
            d_shear_ur = (-3.0 / tau * shear_ur
                          + 2.0 / 3.0 * (theta_ur + metric_ufa_class))

    return d_delta_ur, d_theta_ur, d_shear_ur, d_ur_l


# ===========================================================================
# 8. TOTAL STRESS-ENERGY SUMS (for Einstein equations)
# ===========================================================================
# Source: perturb_total_stress_energy(), lines 6860-7124
"""
These sums feed into the metric equations (Section 1).
CLASS computes them incrementally; here we show the structure.
"""


def compute_total_stress_energy(
    # photon
    rho_g, delta_g, theta_g, shear_g,
    # baryon
    rho_b, delta_b, theta_b, delta_p_b_over_rho_b,
    # CDM
    rho_cdm=0, delta_cdm=0, theta_cdm=0,
    # massless neutrinos
    rho_ur=0, delta_ur=0, theta_ur=0, shear_ur=0,
    # GDM
    rho_gdm=0, delta_gdm=0, theta_gdm=0, shear_gdm=0,
    w_gdm=0, cs2_gdm=0, ca2_gdm=0,
    k=1.0, calH=1.0,
    gauge='newtonian'
):
    """
    Compute total perturbed stress-energy sums.
    Lines 6860-7124.

    Returns dict with: delta_rho, rho_plus_p_theta, rho_plus_p_shear, delta_p
    """
    k2 = k * k

    # Photon + Baryon base (lines 6860-6866)
    delta_rho = rho_g * delta_g + rho_b * delta_b
    rho_plus_p_theta = 4.0 / 3.0 * rho_g * theta_g + rho_b * theta_b
    rho_plus_p_shear = 4.0 / 3.0 * rho_g * shear_g
    delta_p = 1.0 / 3.0 * rho_g * delta_g + rho_b * delta_p_b_over_rho_b

    # CDM (lines 6880-6882)
    delta_rho += rho_cdm * delta_cdm
    if gauge == 'newtonian':
        rho_plus_p_theta += rho_cdm * theta_cdm
    # CDM has no pressure or shear

    # Massless neutrinos (lines 6941-6944)
    delta_rho += rho_ur * delta_ur
    rho_plus_p_theta += 4.0 / 3.0 * rho_ur * theta_ur
    rho_plus_p_shear += 4.0 / 3.0 * rho_ur * shear_ur
    delta_p += 1.0 / 3.0 * rho_ur * delta_ur

    # GDM (lines 7104-7111)
    delta_rho += rho_gdm * delta_gdm
    rho_plus_p_theta += (1.0 + w_gdm) * rho_gdm * theta_gdm
    delta_p += (cs2_gdm * rho_gdm * delta_gdm
                + 3.0 / k2 * calH * (1.0 + w_gdm) * (cs2_gdm - ca2_gdm)
                * rho_gdm * theta_gdm)
    # GDM shear is added in perturb_einstein, not here for algebraic case
    # For dynamic shear (line 7109-7111):
    rho_plus_p_shear += (1.0 + w_gdm) * rho_gdm * shear_gdm

    return {
        'delta_rho': delta_rho,
        'rho_plus_p_theta': rho_plus_p_theta,
        'rho_plus_p_shear': rho_plus_p_shear,
        'delta_p': delta_p,
    }


# ===========================================================================
# 9. THOMSON SCATTERING AND RECOMBINATION
# ===========================================================================
# Source: thermodynamics.c, lines 125-139, 166-185
"""
Opacity (Thomson scattering rate):
  kappa_dot = d(kappa)/d(tau)  [1/Mpc]
            = (1+z)^2 * n_e_today * x_e(z) * sigma_T * (Mpc_over_m)

  where:
    n_e_today = 3*H0^2*Omega_b*(1-Y_He) / (8*pi*G*m_H)  [total hydrogen number density today]
    x_e(z) = free electron fraction from recombination history (Recfast/HyRec)
    sigma_T = Thomson cross section = 6.6524e-29 m^2
    Mpc_over_m = 3.0857e22

Derived quantities (stored in pvecthermo[]):
  index_th_dkappa    : kappa_dot [1/Mpc]                            (line 166)
  index_th_ddkappa   : kappa_dot_dot = d^2(kappa)/d(tau)^2          (line 139)
  index_th_exp_m_kappa : exp(-kappa)                                (line 170)
  index_th_g         : visibility g = kappa_dot * exp(-kappa)       (line 171)
  index_th_cb2       : baryon sound speed squared c_b^2             (line 185)

The scattering appears in Boltzmann equations as:
  - kappa_dot * (Theta_l)  for l >= 2  [damping of anisotropies]
  - kappa_dot * (theta_b - theta_g)    [photon-baryon momentum exchange]
  - kappa_dot * (pol_l - source)       [polarization generation]

Perturbed recombination (optional, lines 8832-8868):
  When ppt->has_perturbed_recombination == True, two extra variables are evolved:
    delta_temp : perturbation to baryon temperature
    delta_chi  : perturbation to ionization fraction x_e

  delta_chi' = -alpha_rec * a * chi * n_H * (delta_alpha + delta_chi + delta_b) * Mpc/c
      (line 9111)

  delta_temp' = 2/3 * delta_b' - a * Compton_CR * (T_cmb/a)^4 * chi/(1+chi+fHe)
                * [(1 - T_cmb*a_today/(a*T_b)) * (delta_g + delta_chi*(1+fHe)/(1+chi+fHe))
                   + T_cmb*a_today/(a*T_b) * (delta_temp - delta_g/4)]
      (line 9114-9116)

  where:
    alpha_rec = case-B recombination coefficient [m^3/s]
              = 1.14 * 4.309e-19 * (T_b/1e4)^(-0.6166) / (1 + 0.6703*(T_b/1e4)^0.53)
    Compton_CR = 8/3 * sigma_T * a_rad / (m_e * c^2) * Mpc_over_m  [1/Mpc / K^4]
    n_H = (a_today/a)^3 * N_now
    fHe = Y_He / (4*(1-Y_He))
"""


# ===========================================================================
# 10. METRIC EVOLUTION EQUATIONS (final ODE)
# ===========================================================================
# Source: lines 9540-9554

def metric_evolution_rhs(gauge, phi_prime=0.0, eta_prime=0.0):
    """
    The metric perturbation evolved as an ODE variable.

    Synchronous gauge (line 9546):
      d(eta)/d(tau) = eta_prime  [computed from Einstein equations in Section 1]

    Newtonian gauge (line 9552):
      d(phi)/d(tau) = phi_prime  [computed from Einstein equations in Section 1]
    """
    if gauge == 'newtonian':
        return phi_prime     # line 9552
    elif gauge == 'synchronous':
        return eta_prime     # line 9546


# ===========================================================================
# 11. COMPLETE ODE SYSTEM (putting it all together)
# ===========================================================================
"""
The full state vector for Newtonian gauge with standard species + GDM is:

  y = [phi,                                    # 1 metric
       delta_b, theta_b,                       # 2 baryon
       delta_cdm, theta_cdm,                   # 2 CDM (theta_cdm=0 in sync gauge)
       delta_gdm, theta_gdm, (shear_gdm),     # 2-3 GDM
       delta_g, theta_g, shear_g,              # 3+  photon temperature
       F_3, ..., F_{l_max_g},                  #     higher photon T multipoles
       G_0, G_1, ..., G_{l_max_pol},           #     photon polarization
       delta_ur, theta_ur, shear_ur,           # 3+  massless neutrino
       N_3, ..., N_{l_max_ur}]                 #     higher nu multipoles

The RHS is dy/dtau = f(y, tau) computed by:

  1. Interpolate background: a(tau), H(a), rho_i(a), w_gdm(a), ca2_gdm(a)
  2. Interpolate thermodynamics: kappa_dot(tau), cb2(tau)
  3. Compute total stress-energy: delta_rho, (rho+p)*theta, (rho+p)*sigma, delta_p
  4. Compute metric from Einstein equations -> Psi, Phi', (or h', eta', alpha, alpha')
  5. Compute metric source terms: metric_continuity, metric_euler, metric_shear
  6. Compute each species' evolution from Sections 2-6 above

Below is the complete assembly pseudocode:
"""


def perturb_derivs(tau, y, params):
    """
    Complete RHS of the perturbation ODE system.
    This is the Python equivalent of perturb_derivs() in perturbations.c:8658.

    Parameters
    ----------
    tau : conformal time [Mpc]
    y : state vector (see above)
    params : dict containing all background, thermo, and config data

    Returns
    -------
    dy : time derivative of state vector
    """
    k = params['k']
    k2 = k * k
    gauge = params['gauge']

    # ---- Step 1: Background ---- (line 8751-8809)
    a = params['a_of_tau'](tau)
    a2 = a * a
    calH = params['calH_of_tau'](tau)   # a' / a
    rho_g = params['rho_g'](a)
    rho_b = params['rho_b'](a)
    rho_cdm = params['rho_cdm'](a)
    rho_ur = params['rho_ur'](a)
    rho_gdm = params['rho_gdm'](a)
    w_gdm = params['w_gdm'](a)
    ca2_gdm = params['ca2_gdm'](a)
    R = 4.0 / 3.0 * rho_g / rho_b       # line 8788

    # ---- Step 2: Thermodynamics ---- (line 8760-8768)
    dkappa = params['dkappa_of_tau'](tau)   # Thomson scattering rate
    cb2 = params['cb2_of_tau'](tau)         # baryon sound speed squared

    # ---- Step 3: Read state variables ---- (line 8814-8829)
    idx = params['idx']  # index map
    phi = y[idx['phi']]
    delta_b = y[idx['delta_b']]
    theta_b = y[idx['theta_b']]
    delta_cdm = y[idx['delta_cdm']]
    theta_cdm = y[idx['theta_cdm']]
    delta_gdm = y[idx['delta_gdm']]
    theta_gdm = y[idx['theta_gdm']]
    delta_g = y[idx['delta_g']]
    theta_g = y[idx['theta_g']]
    shear_g = y[idx['shear_g']]
    delta_ur = y[idx['delta_ur']]
    theta_ur = y[idx['theta_ur']]
    shear_ur = y[idx['shear_ur']]

    # ---- Step 4: GDM sound speeds ---- (line 9181-9182)
    cs2_gdm = params['cs2_gdm'](a, k)     # cs2_gdm_of_a_and_k
    cv2_gdm = params['cv2_gdm'](a, k)     # cv2_gdm_of_a_and_k
    dynamic_shear = params.get('dynamic_shear_gdm', True)

    if dynamic_shear:
        shear_gdm = y[idx['shear_gdm']]
    else:
        # Algebraic shear (line 6473 / 6560)
        if gauge == 'newtonian':
            shear_gdm = gdm_algebraic_shear_newtonian(
                theta_gdm, cv2_gdm, w_gdm, calH)
        else:
            shear_gdm = 0.0  # computed after alpha is known

    # ---- Step 5: Total stress-energy ---- (lines 6860-7124)
    delta_p_b_over_rho_b = cb2 * delta_b
    sums = compute_total_stress_energy(
        rho_g, delta_g, theta_g, shear_g,
        rho_b, delta_b, theta_b, delta_p_b_over_rho_b,
        rho_cdm, delta_cdm, theta_cdm,
        rho_ur, delta_ur, theta_ur, shear_ur,
        rho_gdm, delta_gdm, theta_gdm, shear_gdm,
        w_gdm, cs2_gdm, ca2_gdm, k, calH, gauge)

    # ---- Step 6: Metric Einstein equations ---- (lines 6460-6588)
    if gauge == 'newtonian':
        Psi, phi_prime = metric_newtonian(
            phi, a, k, calH,
            sums['delta_rho'], sums['rho_plus_p_theta'],
            sums['rho_plus_p_shear'])
        mc = -3.0 * phi_prime
        me = k2 * Psi
        ms = 0.0

    elif gauge == 'synchronous':
        eta = phi  # in sync gauge, the evolved variable is eta, not phi
        h_p, eta_p, h_pp, alpha, alpha_p = metric_synchronous(
            eta, a, k, calH,
            sums['delta_rho'], sums['rho_plus_p_theta'],
            sums['rho_plus_p_shear'], sums['delta_p'])
        mc = h_p / 2.0
        me = 0.0
        ms = k2 * alpha

        # Now compute algebraic shear for GDM in sync gauge
        if not dynamic_shear:
            shear_gdm = gdm_algebraic_shear_synchronous(
                theta_gdm, cv2_gdm, w_gdm, calH, k2, alpha)

    # ---- Step 7: Species equations ---- (lines 8919-9216)
    dy = np.zeros_like(y)

    # Metric (line 9546 / 9552)
    if gauge == 'newtonian':
        dy[idx['phi']] = phi_prime
    else:
        dy[idx['phi']] = eta_p  # eta' for sync gauge

    # Photon temperature hierarchy (Section 2, lines 8925-9006)
    dy[idx['delta_g']] = -4.0 / 3.0 * (theta_g + mc)

    dy[idx['theta_g']] = (k2 * (delta_g / 4.0 - shear_g)
                          + me + dkappa * (theta_b - theta_g))

    # ... (full hierarchy requires looping over multipoles, see photon_hierarchy_rhs)

    # Baryon (lines 8931, 8940-8944)
    dy[idx['delta_b']] = -(theta_b + mc)

    dy[idx['theta_b']] = (-calH * theta_b + me
                          + k2 * delta_p_b_over_rho_b
                          + R * dkappa * (theta_g - theta_b))

    # CDM (lines 9063-9065)
    if gauge == 'newtonian':
        dy[idx['delta_cdm']] = -(theta_cdm + mc)
        dy[idx['theta_cdm']] = -calH * theta_cdm + me
    else:
        dy[idx['delta_cdm']] = -mc

    # GDM (lines 9192-9213)
    d_dg, d_tg, d_sg = gdm_rhs(
        delta_gdm, theta_gdm, shear_gdm,
        k, calH, w_gdm, ca2_gdm, cs2_gdm, cv2_gdm,
        mc, me, ms,
        k2_Pinad=params.get('k2_Pinad_gdm', False),
        dynamic_shear=dynamic_shear)
    dy[idx['delta_gdm']] = d_dg
    dy[idx['theta_gdm']] = d_tg
    if dynamic_shear:
        dy[idx['shear_gdm']] = d_sg

    # Massless neutrinos (lines 9331-9368)
    dy[idx['delta_ur']] = -4.0 / 3.0 * (theta_ur + mc)
    dy[idx['theta_ur']] = k2 * (delta_ur / 4.0 - shear_ur) + me
    # ... (full hierarchy requires looping, see ur_hierarchy_rhs)

    return dy


# ===========================================================================
# 12. SUMMARY TABLE: All ODE Equations with Line Numbers
# ===========================================================================
"""
NEWTONIAN GAUGE EQUATIONS (flat space, s_l=1, K=0):

  dkappa = Thomson scattering rate (thermodynamics.c:126)
  R      = (4/3)*rho_g/rho_b                              (line 8788)
  Psi    = Phi - 4.5*(a^2/k^2)*sum[(rho+p)*sigma]         (line 6491)
  Phi'   = -calH*Psi + 1.5*(a^2/k^2)*sum[(rho+p)*theta]   (line 6494)
  mc     = -3*Phi'                                         (line 8898)
  me     = k^2*Psi                                         (line 8899)

  PHOTON TEMPERATURE (lines 8925-9006):
  -------
  delta_g'  = -(4/3)*(theta_g + mc)
  theta_g'  = k^2*(delta_g/4 - shear_g) + me + dkappa*(theta_b - theta_g)
  shear_g'  = (1/2)*[8/15*(theta_g) - 3/5*k*F_3
                      - dkappa*(2*shear_g - 4/5*P0)]
  F_l'      = k/(2l+1)*[l*F_{l-1} - (l+1)*F_{l+1}] - dkappa*F_l   (l >= 3)
  F_lmax'   = k*[F_{lmax-1} - (lmax+1)/(k*tau)*F_lmax] - dkappa*F_lmax

  P0 = (G_0 + G_2 + 2*shear_g) / 8

  PHOTON POLARIZATION (lines 9010-9038):
  -------
  G_0' = -k*G_1 - dkappa*(G_0 - 4*P0)
  G_1' = k/3*(G_0 - 2*G_2) - dkappa*G_1
  G_2' = k/5*(2*G_1 - 3*G_3) - dkappa*(G_2 - 4/5*P0)
  G_l' = k/(2l+1)*(l*G_{l-1} - (l+1)*G_{l+1}) - dkappa*G_l   (l >= 3)
  G_lmax' = k*(G_{lmax-1} - (lmax+1)/(k*tau)*G_lmax) - dkappa*G_lmax

  BARYON (lines 8931, 8940-8944):
  ------
  delta_b' = -(theta_b + mc)
  theta_b' = -calH*theta_b + me + k^2*cb2*delta_b + R*dkappa*(theta_g - theta_b)

  CDM (lines 9063-9065):
  ---
  delta_c' = -(theta_c + mc)
  theta_c' = -calH*theta_c + me

  GDM/KHRONON (lines 9192-9213):
  -----------
  pinad  = (cs2 - ca2)*[delta_gdm + 3*calH*(1+w)*theta_gdm/k^2]
  delta_gdm' = -(1+w)*(theta_gdm + mc) + 3*calH*[(w-ca2)*delta_gdm - pinad]
  theta_gdm' = -(1-3*ca2)*calH*theta_gdm + k^2/(1+w)*(ca2*delta_gdm + pinad) + me - k^2*sigma_gdm
  sigma_gdm' = -3*calH*sigma_gdm + 8/3*cv2/(1+w)*(theta_gdm)      [if dynamic shear]
  sigma_gdm  = 8*cv2/(15*(1+w)*calH)*theta_gdm                     [if algebraic shear]

  MASSLESS NEUTRINOS (lines 9331-9368):
  ------------------
  delta_ur' = -(4/3)*(theta_ur + mc)
  theta_ur' = k^2*(delta_ur/4 - shear_ur) + me
  shear_ur' = (1/2)*[8/15*(theta_ur) - 3/5*k*N_3]
  N_l'      = k/(2l+1)*[l*N_{l-1} - (l+1)*N_{l+1}]   (l >= 3)
  N_lmax'   = k*[N_{lmax-1} - (lmax+1)/(k*tau)*N_lmax]

  METRIC (lines 9550-9552):
  ------
  Phi' = [from Einstein equations above]   (evolved as ODE, line 9552)

=== KEY: KHRONON DBI SOUND SPEED (line 10461-10473) ===

  cs2_gdm(a,k) = ca2 / (1 + ca2*k^2 / (1.5*rho_gdm*a^2*(1+w)))

  This gives:  cs2 -> 0      at large k  (CDM-like on small scales)
               cs2 -> ca2    at small k  (adiabatic on large scales)

  When ca2 -> 0 (as in ghost condensation / Khronon):
    cs2 = 0 identically at all scales
    => pinad = 0, and GDM reduces to pressureless dust (CDM)

=== KEY: ALGEBRAIC vs DYNAMIC SHEAR ===

  If dynamic_shear_gdm = False (default for simple GDM):
    sigma_gdm is NOT evolved; instead computed algebraically from theta_gdm
    sigma_gdm = 8*cv2 / (15*(1+w)*calH) * theta_gdm   [Newtonian gauge]

  If dynamic_shear_gdm = True:
    sigma_gdm is an additional ODE variable with:
    sigma_gdm' = -3*calH*sigma_gdm + (8/3)*cv2/(1+w)*(theta_gdm + metric_shear)
"""
