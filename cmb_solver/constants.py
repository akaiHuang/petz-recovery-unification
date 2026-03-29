"""
Physical constants and cosmological parameters.

All units: natural units where c = hbar = k_B = 1,
except where explicitly noted. Distances in Mpc, time in Mpc/c.
"""
import numpy as np

# ---------------------------------------------------------------------------
# Fundamental constants (SI)
# ---------------------------------------------------------------------------
c_SI = 2.99792458e8          # m/s
hbar_SI = 1.054571817e-34    # J s
kB_SI = 1.380649e-23         # J/K
G_SI = 6.67430e-11           # m^3 kg^-1 s^-2
sigma_T_SI = 6.6524587e-29   # Thomson cross-section, m^2
m_e_SI = 9.1093837e-31       # electron mass, kg
m_p_SI = 1.6726219e-27       # proton mass, kg
m_H_SI = 1.6735575e-27       # hydrogen atom mass, kg

# ---------------------------------------------------------------------------
# Unit conversions
# ---------------------------------------------------------------------------
Mpc_SI = 3.0856775814913673e22   # 1 Mpc in metres
km_per_s_per_Mpc = 1e3 / Mpc_SI  # (km/s/Mpc) -> 1/s

# ---------------------------------------------------------------------------
# Cosmological parameters (Planck 2018 best-fit, flat LCDM baseline)
# ---------------------------------------------------------------------------
H0_km_s_Mpc = 67.36           # km/s/Mpc
H0 = H0_km_s_Mpc * km_per_s_per_Mpc   # 1/s
H0_Mpc = H0_km_s_Mpc * 1e3 / c_SI / Mpc_SI * Mpc_SI  # actually we want H0 in 1/Mpc units
# H0 in inverse-Mpc (natural time unit = Mpc/c):
H0_inv_Mpc = H0_km_s_Mpc / (c_SI * 1e-3)  # 1/Mpc  (H0/(c) in 1/Mpc)
# More carefully: H0 = 67.36 km/s/Mpc = 67.36e3 m/s / (3.0857e22 m) = 2.184e-18 /s
# In conformal-time units (tau in Mpc/c): H_conf = a*H, we work with conformal Hubble aH.
# For the ODE we'll use tau in Mpc, so H0 in 1/Mpc:
H0_per_Mpc = H0_km_s_Mpc * 1e3 / c_SI / Mpc_SI * Mpc_SI
# Simpler: H0 [1/s] * (Mpc [m]) / c [m/s] = H0 in 1/Mpc... let me just compute:
_H0_si = H0_km_s_Mpc * 1e3 / Mpc_SI   # H0 in 1/s
H0_invMpc = _H0_si * Mpc_SI / c_SI     # H0 in units where time = Mpc/c
# H0_invMpc ~ 67.36 / 299792.458 ~ 2.247e-4 /Mpc

T_CMB = 2.7255                 # K
Omega_b = 0.0493
Omega_K_khronon = 0.265        # Khronon sector (replaces CDM)
Omega_gamma = 5.38e-5
Omega_Lambda = 1.0 - Omega_b - Omega_K_khronon - Omega_gamma
# Neutrinos (massless, N_eff = 3.046):
N_eff = 3.046
Omega_nu = N_eff * (7.0 / 8.0) * (4.0 / 11.0) ** (4.0 / 3.0) * Omega_gamma
# Omega_r = Omega_gamma + Omega_nu (total radiation)
Omega_r = Omega_gamma + Omega_nu

# Recombination parameters (Saha approx)
Y_He = 0.2454                 # Helium mass fraction
n_H0 = (1 - Y_He) * 3 * H0_invMpc**2 * Omega_b / (8 * np.pi * G_SI) \
    if False else None  # We'll compute n_H from rho_b

# Primordial power spectrum
A_s = 2.1e-9                  # amplitude at pivot
n_s = 0.9649                  # spectral index
k_pivot = 0.05                # 1/Mpc

# ---------------------------------------------------------------------------
# Khronon / Sigma = 2 ln Q  parameters
# ---------------------------------------------------------------------------
# mu parameter: BS2025 value or H0/c
# mu = 1 / 22.3   # 1/Mpc  (BS2025 value)
mu_H0_over_c = H0_invMpc      # H0/c in 1/Mpc  (~2.25e-4 /Mpc)
mu = mu_H0_over_c             # default: mu = H0/c
lambda_DBI = 1.0              # DBI parameter

# ---------------------------------------------------------------------------
# Derived quantities
# ---------------------------------------------------------------------------
rho_crit_over_3H2 = 1.0 / (8 * np.pi)  # rho_crit = 3H0^2/(8piG) in G=1 units
# We'll work in units where 8piG/3 = H0^2 / rho_crit, i.e. Friedmann is
# H^2 = H0^2 [Omega_r a^-4 + Omega_b a^-3 + Omega_K rho_K/rho_K0 + Omega_Lambda]

# Baryon-photon ratio
R_factor = 3 * Omega_b / (4 * Omega_gamma)  # R = R_factor * a

# Recombination redshift (approximate)
z_rec = 1089.92
a_rec = 1.0 / (1.0 + z_rec)

# Today
a_today = 1.0

# Initial scale factor for integration
a_init = 1e-6   # deep in radiation domination

# ---------------------------------------------------------------------------
# Convenient print
# ---------------------------------------------------------------------------
def print_params():
    """Print all cosmological parameters."""
    print("=" * 60)
    print("CMB Solver: Sigma = 2 ln Q  (Khronon cosmology)")
    print("=" * 60)
    print(f"H0           = {H0_km_s_Mpc} km/s/Mpc")
    print(f"H0 (1/Mpc)   = {H0_invMpc:.6e}")
    print(f"Omega_b      = {Omega_b}")
    print(f"Omega_K(khr) = {Omega_K_khronon}")
    print(f"Omega_gamma  = {Omega_gamma}")
    print(f"Omega_nu     = {Omega_nu:.6e}")
    print(f"Omega_Lambda = {Omega_Lambda:.6f}")
    print(f"T_CMB        = {T_CMB} K")
    print(f"mu           = {mu:.6e} /Mpc")
    print(f"A_s          = {A_s}")
    print(f"n_s          = {n_s}")
    print(f"z_rec        = {z_rec}")
    print("=" * 60)
