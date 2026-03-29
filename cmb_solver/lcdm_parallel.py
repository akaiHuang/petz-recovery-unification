#!/usr/bin/env python3
"""
LCDM CMB solver -- parallelized with multiprocessing
Uses analytic transfer function (fast), parallelizes C_l over ell.
"""
import os
os.environ['MPLCONFIGDIR'] = '/tmp/mpl_config'

import numpy as np
from scipy.special import spherical_jn
from scipy.interpolate import interp1d
from multiprocessing import Pool, cpu_count
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import time

print(f"CPUs available: {cpu_count()}")
t_start = time.time()

# ===== Cosmological Parameters =====
h = 0.6736
H0 = h * 100  # km/s/Mpc
H0_inv_Mpc = 2997.9 / h  # c/H0 in Mpc
T_CMB = 2.7255  # K
omega_b = 0.02237
omega_c = 0.12
omega_r = 4.15e-5 * h**(-2)  # photons + 3 massless nu
Omega_b = omega_b / h**2
Omega_c = omega_c / h**2
Omega_r = omega_r / h**2
Omega_L = 1.0 - Omega_b - Omega_c - Omega_r
A_s = 2.1e-9
n_s = 0.9649

# ===== Background =====
def E(a):
    """H(a)/H0"""
    return np.sqrt(Omega_r/a**4 + (Omega_b + Omega_c)/a**3 + Omega_L)

# Conformal time: tau = integral da/(a^2 H)
a_arr = np.logspace(-7, 0, 10000)
dtau = np.diff(a_arr) / (a_arr[:-1]**2 * E(a_arr[:-1]) * H0 / 2997.9)  # in Mpc
tau_arr = np.concatenate([[0], np.cumsum(dtau)])

# Recombination
a_rec = 1.0/1101.0
idx_rec = np.searchsorted(a_arr, a_rec)
tau_rec = tau_arr[idx_rec]
tau_0 = tau_arr[-1]
D_A = tau_0 - tau_rec  # Angular diameter distance to LSS (comoving, Mpc)

# Sound horizon
R_arr = 3.0 * Omega_b * a_arr / (4.0 * Omega_r)  # baryon loading
cs_arr = 1.0 / np.sqrt(3.0 * (1.0 + R_arr))
dtau_cs = cs_arr[:-1] * dtau
r_s = np.sum(dtau_cs[:idx_rec])  # sound horizon at recombination

# Equality
a_eq = Omega_r / (Omega_b + Omega_c)
k_eq = 1.0 / (a_eq * H0_inv_Mpc * E(a_eq))

# Silk damping scale
k_D = 0.15 * (omega_b / 0.022)**0.25  # approximate, Mpc^-1

R_rec = R_arr[idx_rec]
print(f"D_A = {D_A:.1f} Mpc")
print(f"r_s = {r_s:.1f} Mpc")
print(f"l_1 ~ pi*D_A/r_s = {np.pi*D_A/r_s:.0f}")
print(f"a_eq = {a_eq:.5f}, k_eq = {k_eq:.5f} Mpc^-1")
print(f"k_D = {k_D:.4f} Mpc^-1")
print(f"R(rec) = {R_rec:.3f}")

# ===== Source Function (analytic) =====
N_k = 200
k_min = 1e-4
k_max = 0.3
k_arr = np.geomspace(k_min, k_max, N_k)
Phi_prim = 1.0

def compute_source(k):
    kr_s = k * r_s
    x_eq = k / k_eq

    # Baryon damping
    baryon_damp = (1.0 + R_rec)**(-0.25)
    zero_shift = R_rec * Phi_prim / (3.0 * (1.0 + R_rec))

    # Driving boost
    driving_boost = 1.0 + 2.0 * x_eq**2 / (1.0 + x_eq**2)

    # Sachs-Wolfe + oscillation
    A_osc = (Phi_prim / 3.0) * driving_boost * baryon_damp
    D_silk = np.exp(-(k/k_D)**2)

    S = (A_osc * np.cos(kr_s) + zero_shift) * D_silk
    return S

source_arr = np.array([compute_source(k) for k in k_arr])
print(f"Source computed in {time.time()-t_start:.1f}s")

# ===== C_l Computation (parallelized) =====
print(f"\nComputing C_l with {cpu_count()} cores...")
t0 = time.time()

k_pivot = 0.05
P_prim = A_s * (k_arr / k_pivot)**(n_s - 1.0)
lnk = np.log(k_arr)
dlnk = np.diff(lnk)
integrand_base = P_prim * source_arr**2

# l grid
l_sparse = np.unique(np.concatenate([
    np.arange(2, 30, 1),
    np.arange(30, 100, 3),
    np.arange(100, 500, 5),
    np.arange(500, 1500, 10),
    np.arange(1500, 2501, 15),
])).astype(int)

def compute_Cl(l):
    """Compute C_l for a single multipole."""
    x = k_arr * D_A
    jl = spherical_jn(l, x)
    integrand = integrand_base * jl**2
    mid = 0.5 * (integrand[:-1] + integrand[1:])
    return 4.0 * np.pi * np.sum(mid * dlnk)

# Parallel execution
with Pool(cpu_count()) as pool:
    Cl_sparse = np.array(pool.map(compute_Cl, l_sparse))

print(f"C_l done in {time.time()-t0:.1f}s")

# Interpolate
l_full = np.arange(2, 2501)
Cl_interp = interp1d(l_sparse, Cl_sparse, kind='cubic', fill_value='extrapolate')
Cl_full = np.maximum(Cl_interp(l_full), 0)
Dl_full = l_full * (l_full + 1.0) * Cl_full / (2.0*np.pi) * (T_CMB*1e6)**2

# ===== Find Peaks =====
from scipy.signal import find_peaks
peaks, _ = find_peaks(Dl_full, distance=80, prominence=100)
print(f"\nPeaks found at l = {l_full[peaks][:5]}")
print(f"Peak heights: {Dl_full[peaks][:5]:.0f}")
if len(peaks) >= 2:
    print(f"Peak ratio (1st/2nd): {Dl_full[peaks[0]]/Dl_full[peaks[1]]:.2f}")

# ===== Plot =====
fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(l_full, Dl_full, 'b-', lw=1.5, label='Our solver (LCDM)')
ax.axvline(220, color='gray', ls='--', alpha=0.5, label='l=220 (expected 1st peak)')
ax.axvline(540, color='gray', ls=':', alpha=0.5, label='l=540 (expected 2nd peak)')
for p in peaks[:5]:
    ax.axvline(l_full[p], color='red', ls='--', alpha=0.3)
    ax.annotate(f'l={l_full[p]}', (l_full[p], Dl_full[p]), fontsize=8)
ax.set_xlabel('Multipole l')
ax.set_ylabel(r'$D_l = l(l+1)C_l/2\pi$ [$\mu K^2$]')
ax.set_title('CMB TT Power Spectrum (LCDM, analytic solver)')
ax.legend()
ax.set_xlim(2, 2500)
ax.set_ylim(0, max(Dl_full)*1.2)
plt.tight_layout()
plt.savefig('/Users/akaihuangm1/Desktop/github/petz-recovery-unification/cmb_solver/lcdm_parallel_result.png', dpi=150)
print(f"\nPlot saved. Total time: {time.time()-t_start:.1f}s")
