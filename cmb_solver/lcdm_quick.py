#!/usr/bin/env python3
"""
Minimal CMB TT Power Spectrum Calculator -- LCDM baseline
=========================================================

Physics:
  - FRW background with radiation, matter, Lambda
  - Photon Boltzmann hierarchy (Theta_0, Theta_1, Theta_2) + CDM + baryons
  - Tight-coupling approximation before recombination
  - Newtonian gauge, Phi = Psi (no anisotropic stress)
  - Adiabatic initial conditions
  - C_l via source projection: integrate |S(k)|^2 * j_l^2 over k

Target: ~10-20% accuracy on first 3 acoustic peak positions.

Author: Sheng-Kai Huang / Claude
Date: 2026-03-19
"""

import os
os.environ['MPLCONFIGDIR'] = '/tmp/mpl_config'

import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.special import spherical_jn
from scipy.interpolate import interp1d
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import time

print("="*65)
print(" Minimal LCDM CMB TT Spectrum Calculator")
print("="*65)

# ===========================================================================
# 1. Cosmological parameters (Planck 2018 best-fit)
# ===========================================================================
h      = 0.6736
Ob_h2  = 0.02237
Oc_h2  = 0.12
T_CMB  = 2.7255
n_s    = 0.9649
A_s    = 2.1e-9

Omega_b = Ob_h2 / h**2
Omega_c = Oc_h2 / h**2
Og_h2   = 2.469e-5 * (T_CMB / 2.7255)**4
Omega_g = Og_h2 / h**2
N_eff   = 3.044
Omega_nu = N_eff * (7.0/8.0) * (4.0/11.0)**(4.0/3.0) * Omega_g
Omega_r = Omega_g + Omega_nu
Omega_m = Omega_b + Omega_c
Omega_L = 1.0 - Omega_m - Omega_r
f_nu    = Omega_nu / Omega_r

# Natural units: c=1, distances in Mpc
H0 = h * 100.0 / 299792.458  # [1/Mpc]

print(f"\nParameters:")
print(f"  h={h}, Ob={Omega_b:.4f}, Oc={Omega_c:.4f}, Or={Omega_r:.6f}, OL={Omega_L:.4f}")

# ===========================================================================
# 2. Background cosmology
# ===========================================================================
def Hubble(a):
    return H0 * np.sqrt(Omega_r/a**4 + Omega_m/a**3 + Omega_L)

def aHubble(a):
    return a * Hubble(a)

def dtau_da(a):
    return 1.0 / (a**2 * Hubble(a))

print("\nComputing background...")
t0 = time.time()

a_ini = 1e-7
N_bg = 3000
a_arr = np.geomspace(a_ini, 1.0, N_bg)
tau_arr = np.zeros(N_bg)
tau_arr[0] = a_ini / (H0 * np.sqrt(Omega_r))  # rad-dom analytic

for i in range(1, N_bg):
    val, _ = quad(dtau_da, a_arr[i-1], a_arr[i], limit=80)
    tau_arr[i] = tau_arr[i-1] + val

tau_of_a = interp1d(np.log(a_arr), tau_arr, kind='cubic', fill_value='extrapolate')
a_of_tau = interp1d(tau_arr, np.log(a_arr), kind='cubic', fill_value='extrapolate')

def get_a(tau):
    return np.exp(float(a_of_tau(tau)))

def get_tau(a):
    return float(tau_of_a(np.log(a)))

a_rec   = 1.0 / 1101.0
tau_rec = get_tau(a_rec)
tau_0   = tau_arr[-1]
D_A     = tau_0 - tau_rec

# Sound horizon
R_of_a = lambda a: 3.0*Omega_b*a / (4.0*Omega_g)
cs_of_a = lambda a: 1.0/np.sqrt(3.0*(1.0 + R_of_a(a)))

N_sh = 1000
a_sh = np.geomspace(a_ini, a_rec, N_sh)
tau_sh = np.array([get_tau(a) for a in a_sh])
cs_sh  = np.array([cs_of_a(a) for a in a_sh])
r_s = np.trapezoid(cs_sh, tau_sh)

l_first_peak = np.pi * D_A / r_s
print(f"  tau_rec = {tau_rec:.1f} Mpc, tau_0 = {tau_0:.1f} Mpc")
print(f"  D_A = {D_A:.0f} Mpc, r_s = {r_s:.1f} Mpc")
print(f"  Predicted l_1 = pi*D_A/r_s = {l_first_peak:.0f}")
print(f"  R(a_rec) = {R_of_a(a_rec):.3f}")
print(f"  Background done in {time.time()-t0:.1f}s")

# ===========================================================================
# 3. Perturbation solver
# ===========================================================================
# Tight-coupling + Boltzmann hierarchy approach
# Variables: [Theta_0, Theta_1, Theta_2, delta_c, v_c, Phi]
# In tight coupling: baryon velocity v_b = 3*Theta_1

# Pre-compute background on a fine tau grid for fast lookup
N_fine = 10000
tau_fine = np.linspace(tau_arr[0], tau_rec * 1.3, N_fine)
a_fine = np.array([get_a(t) for t in tau_fine])
aH_fine = np.array([aHubble(a) for a in a_fine])
R_fine  = np.array([R_of_a(a) for a in a_fine])

# Scattering rate kappa_dot
def kappa_dot_of_a(a):
    z = 1.0/a - 1.0
    z_rec_c = 1100.0
    dz = 80.0
    X_e = 0.5 * (1.0 - np.tanh((z_rec_c - z) / dz))
    return 0.0691 * H0 * (Ob_h2/0.022) * (1+z)**2 * X_e

kd_fine = np.array([kappa_dot_of_a(a) for a in a_fine])

# Fast interpolators
_aH_interp = interp1d(tau_fine, aH_fine, kind='linear', fill_value='extrapolate')
_R_interp  = interp1d(tau_fine, R_fine,  kind='linear', fill_value='extrapolate')
_a_interp  = interp1d(tau_fine, a_fine,  kind='linear', fill_value='extrapolate')
_kd_interp = interp1d(tau_fine, kd_fine, kind='linear', fill_value='extrapolate')


def solve_mode(k, tau_start, tau_end):
    """Solve perturbation equations for one k-mode."""

    def derivs(tau, y):
        Th0, Th1, Th2, dc, vc, Phi = y

        a  = float(_a_interp(tau))
        aH = float(_aH_interp(tau))
        R  = float(_R_interp(tau))
        kd = float(_kd_interp(tau))

        Psi = Phi  # no anisotropic stress

        # ---- Einstein: Phi_dot from momentum constraint ----
        # Phi' + aH*Psi = -(4piG a^2 / k) sum rho_i (1+w_i) v_i
        # = -(1.5*H0^2/k) * [Omega_g/a^2 * (4/3)*3*Th1 + Omega_nu/a^2*(4/3)*3*Th1
        #                     + Omega_c/a * vc + Omega_b/a * 3*Th1]
        v_b = 3.0 * Th1
        vel_src = (1.5*H0**2/k * (
            (Omega_g + Omega_nu)/a**2 * 4.0 * Th1   # (4/3)*3*Th1 = 4*Th1
            + Omega_c/a * vc
            + Omega_b/a * v_b
        ))
        Phi_dot = -aH * Psi - vel_src

        # ---- Photon hierarchy ----
        dTh0 = -k * Th1 - Phi_dot

        if kd > 5.0 * k:
            # Tight coupling
            dTh1 = (k/(3.0*(1.0+R))) * (Th0 - 2.0*Th2 + Psi) - R*aH*Th1/(1.0+R)
        else:
            dTh1 = (k/3.0) * (Th0 - 2.0*Th2 + Psi) - kd*(Th1 - v_b/3.0)

        # Theta_2: truncated hierarchy, with damping
        dTh2 = (2.0*k/5.0)*Th1 - (3.0/5.0)*k*0.0 - (9.0/10.0)*kd*Th2  # Theta_3 ~ 0
        # Add truncation damping
        dTh2 -= 3.0/(tau + 1.0) * Th2

        # ---- CDM ----
        ddc = -k * vc - 3.0 * Phi_dot
        dvc = -aH * vc + k * Psi

        # ---- Phi ----
        dPhi = Phi_dot

        return [dTh0, dTh1, dTh2, ddc, dvc, dPhi]

    # Initial conditions (adiabatic, radiation dominated, super-horizon)
    Phi0 = 1.0
    # Neutrino correction
    Phi0 *= (10.0 + 4.0*f_nu) / (15.0 + 4.0*f_nu)

    y0 = [
        -Phi0/2.0,                      # Theta_0
        k * tau_start * Phi0 / 6.0,     # Theta_1
        0.0,                             # Theta_2
        -3.0*Phi0/2.0,                  # delta_c
        k * tau_start * Phi0 / 2.0,     # v_c
        Phi0,                            # Phi
    ]

    sol = solve_ivp(derivs, (tau_start, tau_end), y0,
                    method='RK45', rtol=1e-5, atol=1e-8,
                    max_step=(tau_end-tau_start)/150,
                    dense_output=True)

    if not sol.success:
        print(f"  WARNING k={k:.5f}: {sol.message}")

    return sol


# ===========================================================================
# 4. Compute transfer function
# ===========================================================================
print("\nComputing transfer function T(k)...")
t0 = time.time()

N_k = 30
k_min = 2.0 / D_A
k_max = 2800.0 / D_A
k_arr = np.geomspace(k_min, k_max, N_k)

tau_start = tau_arr[0]
tau_end   = tau_rec * 1.2

source_arr = np.zeros(N_k)  # Theta_0 + Psi at recombination
Phi_rec_arr = np.zeros(N_k)

for i, k in enumerate(k_arr):
    if (i+1) % 50 == 0 or i == 0:
        elapsed = time.time() - t0
        print(f"  [{i+1}/{N_k}] k={k:.5f} 1/Mpc (l~{k*D_A:.0f})  [{elapsed:.1f}s]")

    sol = solve_mode(k, tau_start, tau_end)

    # Evaluate at tau_rec using dense output
    y_rec = sol.sol(tau_rec)
    Th0 = y_rec[0]
    Phi = y_rec[5]
    Psi = Phi

    source_arr[i] = Th0 + Psi  # Sachs-Wolfe source
    Phi_rec_arr[i] = Phi

elapsed = time.time() - t0
print(f"  Transfer function done in {elapsed:.1f}s")

# Quick check: plot source
print(f"  Source range: [{source_arr.min():.4f}, {source_arr.max():.4f}]")
print(f"  Source at k*D_A~220: {source_arr[np.argmin(np.abs(k_arr*D_A - 220))]:.4f}")

# ===========================================================================
# 5. Compute C_l via projection
# ===========================================================================
print("\nComputing C_l...")
t0 = time.time()

k_pivot = 0.05  # Mpc^-1
P_prim = A_s * (k_arr / k_pivot)**(n_s - 1.0)

# Precompute log-k weights for trapezoidal integration in log-space
lnk = np.log(k_arr)
dlnk = np.diff(lnk)

# l values -- use sparser grid then interpolate
l_sparse = np.unique(np.concatenate([
    np.arange(2, 50, 1),
    np.arange(50, 200, 3),
    np.arange(200, 800, 5),
    np.arange(800, 2501, 8),
])).astype(int)

Cl_sparse = np.zeros(len(l_sparse))

# The integrand for each l: P_prim(k) * source(k)^2 * j_l(k*D_A)^2
# Vectorized over k for each l
integrand_base = P_prim * source_arr**2  # shape (N_k,)

for il, l in enumerate(l_sparse):
    if il % 100 == 0 and il > 0:
        print(f"    l = {l} [{time.time()-t0:.1f}s]")

    # Vectorized spherical Bessel
    x = k_arr * D_A
    jl = spherical_jn(l, x)

    integrand = integrand_base * jl**2

    # Trapezoidal in ln(k)
    mid = 0.5 * (integrand[:-1] + integrand[1:])
    Cl_sparse[il] = 4.0 * np.pi * np.sum(mid * dlnk)

# Interpolate to full l range
l_full = np.arange(2, 2501)
Cl_interp = interp1d(l_sparse, Cl_sparse, kind='cubic', fill_value='extrapolate')
Cl_full = Cl_interp(l_full)
Cl_full = np.maximum(Cl_full, 0)  # avoid negative from interpolation

# D_l = l(l+1)C_l/(2pi) in muK^2
Dl_full = l_full * (l_full + 1.0) * Cl_full / (2.0*np.pi) * (T_CMB * 1e6)**2

elapsed = time.time() - t0
print(f"  C_l done in {elapsed:.1f}s")

# ===========================================================================
# 6. Peak analysis
# ===========================================================================
print("\n" + "="*65)
print(" Peak Analysis")
print("="*65)

from scipy.signal import savgol_filter, find_peaks

Dl_smooth = savgol_filter(Dl_full, min(101, len(Dl_full)//2*2-1), 3)
Dl_smooth = np.maximum(Dl_smooth, 0)

peaks_idx, _ = find_peaks(Dl_smooth, distance=100, prominence=100)
peaks_idx = peaks_idx[l_full[peaks_idx] > 100]  # skip SW plateau

expected = {1: (220, 5800), 2: (540, 2500), 3: (810, 4000)}

print(f"\n  {'#':<4} {'l_found':<10} {'l_expect':<10} {'D_l found':<12} {'D_l exp':<10} {'l_err':<8}")
print(f"  {'-'*54}")
for i, pk in enumerate(peaks_idx[:5]):
    l_f = l_full[pk]
    D_f = Dl_smooth[pk]
    if i+1 in expected:
        l_e, D_e = expected[i+1]
        err = abs(l_f - l_e)/l_e * 100
        print(f"  {i+1:<4} {l_f:<10} {l_e:<10} {D_f:<12.0f} {D_e:<10} {err:<8.1f}%")
    else:
        print(f"  {i+1:<4} {l_f:<10} {'--':<10} {D_f:<12.0f} {'--':<10}")

if len(peaks_idx) >= 2:
    r12 = Dl_smooth[peaks_idx[0]] / Dl_smooth[peaks_idx[1]]
    print(f"\n  Peak1/Peak2 ratio: {r12:.2f} (expected ~2.3)")
if len(peaks_idx) >= 3:
    r13 = Dl_smooth[peaks_idx[0]] / Dl_smooth[peaks_idx[2]]
    print(f"  Peak1/Peak3 ratio: {r13:.2f} (expected ~1.5)")

# ===========================================================================
# 7. Plot
# ===========================================================================
print("\nGenerating plot...")

fig, axes = plt.subplots(2, 1, figsize=(12, 10), gridspec_kw={'height_ratios': [3, 1]})

ax = axes[0]
ax.plot(l_full, Dl_full, 'b-', lw=0.5, alpha=0.4, label='Raw $C_\\ell$')
ax.plot(l_full, Dl_smooth, 'b-', lw=1.8, label='Smoothed')

# Mark found peaks
for i, pk in enumerate(peaks_idx[:5]):
    ax.plot(l_full[pk], Dl_smooth[pk], 'ro', ms=8)
    ax.annotate(f'$\\ell$={l_full[pk]}', (l_full[pk], Dl_smooth[pk]),
                textcoords="offset points", xytext=(8, 8), fontsize=10, color='r')

# Mark expected peaks
for i, (l_e, D_e) in expected.items():
    ax.axvline(l_e, color='green', ls=':', alpha=0.5, lw=1)
ax.plot([220, 540, 810], [5800, 2500, 4000], 'g^', ms=10, zorder=5,
        label='Planck peaks (approx)')

ax.set_xlabel(r'Multipole $\ell$', fontsize=14)
ax.set_ylabel(r'$\mathcal{D}_\ell = \ell(\ell+1)C_\ell/2\pi$ [$\mu$K$^2$]', fontsize=14)
ax.set_title(r'CMB TT Power Spectrum --- Minimal $\Lambda$CDM Solver', fontsize=15)
ax.set_xlim(2, 2500)
ax.set_ylim(bottom=0)
ax.legend(fontsize=11, loc='upper right')
ax.grid(True, alpha=0.3)

# Bottom panel: transfer function
ax2 = axes[1]
ax2.plot(k_arr * D_A, source_arr, 'b-', lw=1.2)
ax2.axhline(0, color='k', lw=0.5)
ax2.set_xlabel(r'$k \cdot D_A$ (effective $\ell$)', fontsize=13)
ax2.set_ylabel(r'$\Theta_0 + \Psi$', fontsize=13)
ax2.set_title('Source function at recombination', fontsize=13)
ax2.set_xlim(2, 2500)
ax2.grid(True, alpha=0.3)

plt.tight_layout()
outpath = '/Users/akaihuangm1/Desktop/github/petz-recovery-unification/cmb_solver/lcdm_tt_spectrum.png'
plt.savefig(outpath, dpi=150, bbox_inches='tight')
print(f"  Plot saved: {outpath}")

# Save data
np.savez('/Users/akaihuangm1/Desktop/github/petz-recovery-unification/cmb_solver/lcdm_results.npz',
         l=l_full, Dl_muK2=Dl_full, Dl_smooth=Dl_smooth,
         k=k_arr, source=source_arr, D_A=D_A, tau_rec=tau_rec, tau_0=tau_0, r_s=r_s)
print("  Data saved: lcdm_results.npz")

print("\n" + "="*65)
print(" SUMMARY")
print("="*65)
print(f"  Sound horizon r_s = {r_s:.1f} Mpc")
print(f"  Angular distance D_A = {D_A:.0f} Mpc")
print(f"  Predicted l_1 = pi*D_A/r_s = {l_first_peak:.0f}")
print(f"  Framework ready for CDM -> Khronon replacement.")
print("="*65)
