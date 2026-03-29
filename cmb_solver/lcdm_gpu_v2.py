#!/usr/bin/env python3
"""
LCDM CMB TT solver — Hybrid: scipy Bessel (exact) + MLX GPU (fast integration)
Total time: < 1 second on Apple Silicon
"""
import os
os.environ['MPLCONFIGDIR'] = '/tmp/mpl_config'

import numpy as np
import mlx.core as mx
from scipy.special import spherical_jn
from scipy.interpolate import interp1d
from scipy.signal import find_peaks
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import time

print(f"MLX: {mx.default_device()}")
t_start = time.time()

# ===== Parameters =====
h = 0.6736; T_CMB = 2.7255; A_s = 2.1e-9; n_s = 0.9649
omega_b = 0.02237; omega_c = 0.12
c_km = 299792.458
H0_km = h * 100; H0_Mpc = H0_km / c_km
Omega_b = omega_b/h**2; Omega_c = omega_c/h**2
Omega_r = 2.469e-5 * (1 + 0.2271 * 3.046) / h**2  # photons + neutrinos
Omega_m = Omega_b + Omega_c; Omega_L = 1 - Omega_m - Omega_r

# ===== Background =====
a = np.logspace(-7, 0, 50000)
E = np.sqrt(Omega_r/a**4 + Omega_m/a**3 + Omega_L)
dtau = np.diff(a) / (a[:-1]**2 * E[:-1] * H0_Mpc)
tau = np.concatenate([[0], np.cumsum(dtau)])

a_rec = 1/1101.0
idx_rec = np.searchsorted(a, a_rec)
tau_rec = tau[idx_rec]; tau_0 = tau[-1]
D_A = tau_0 - tau_rec

R = 3*Omega_b*a/(4*Omega_r)
cs = 1/np.sqrt(3*(1+R))
r_s = np.sum(cs[:idx_rec-1]*dtau[:idx_rec-1])
R_rec = R[idx_rec]

a_eq = Omega_r / Omega_m
k_eq = np.sqrt(2*Omega_m*H0_Mpc**2/a_eq)
k_D = 0.15 * (omega_b/0.022)**0.25

print(f"D_A={D_A:.0f} Mpc, r_s={r_s:.1f} Mpc, l₁~{np.pi*D_A/r_s:.0f}")
print(f"R_rec={R_rec:.3f}, k_eq={k_eq:.5f}, k_D={k_D:.4f}")

# ===== Source function (analytic, all k at once) =====
N_k = 500
k_arr = np.geomspace(5e-5, 0.35, N_k)

kr_s = k_arr * r_s
x_eq = k_arr / k_eq
baryon_damp = (1 + R_rec)**(-0.25)
zero_shift = R_rec / (3*(1+R_rec))

# Driving boost from potential decay in radiation era
driving = 1.0 + 2.5 * x_eq**2 / (1 + x_eq**2)

# Sachs-Wolfe + acoustic oscillation
A_osc = (1.0/3.0) * driving * baryon_damp

# Silk damping
silk = np.exp(-(k_arr/k_D)**2)

# Source = [amplitude * cos(k*r_s) + zero_point] * damping
source = (A_osc * np.cos(kr_s) + zero_shift) * silk

# Primordial power spectrum
P_prim = A_s * (k_arr/0.05)**(n_s-1)
integrand_base = P_prim * source**2

print(f"Source: {time.time()-t_start:.2f}s")

# ===== C_l: exact Bessel + GPU integration =====
print("Computing C_l (hybrid scipy+GPU)...")
t0 = time.time()

x_arr = k_arr * D_A
lnk = np.log(k_arr)
dlnk = np.diff(lnk)

l_sparse = np.unique(np.concatenate([
    np.arange(2, 30, 1), np.arange(30, 100, 2),
    np.arange(100, 500, 4), np.arange(500, 1500, 8),
    np.arange(1500, 2501, 12)
])).astype(int)

# Step 1: Bessel table (scipy, exact)
jl2_table = np.zeros((len(l_sparse), N_k))
for il, l in enumerate(l_sparse):
    jl2_table[il] = spherical_jn(l, x_arr)**2
t_bessel = time.time() - t0
print(f"  Bessel: {t_bessel:.2f}s")

# Step 2: GPU integration (all l at once)
t1 = time.time()
jl2_gpu = mx.array(jl2_table)
base_gpu = mx.array(integrand_base)
dlnk_gpu = mx.array(dlnk)

integrand_all = jl2_gpu * base_gpu  # (N_l, N_k) broadcast
mid = 0.5 * (integrand_all[:, :-1] + integrand_all[:, 1:])
Cl_sparse = np.array(4*np.pi * mx.sum(mid * dlnk_gpu, axis=1))
t_gpu = time.time() - t1
print(f"  GPU integration: {t_gpu:.3f}s")
print(f"  Total C_l: {time.time()-t0:.2f}s")

# Interpolate to full l
l_full = np.arange(2, 2501)
Cl_interp = interp1d(l_sparse, Cl_sparse, kind='cubic', fill_value='extrapolate')
Cl_full = np.maximum(Cl_interp(l_full), 0)
Dl_full = l_full*(l_full+1)*Cl_full/(2*np.pi) * (T_CMB*1e6)**2

# ===== Results =====
peaks, _ = find_peaks(Dl_full, distance=60, prominence=20)
print(f"\n{'='*50}")
print(f"RESULTS")
print(f"{'='*50}")
print(f"Peaks at l = {l_full[peaks][:5]}")
if len(peaks) >= 3:
    print(f"Heights: {[f'{Dl_full[p]:.0f} μK²' for p in peaks[:3]]}")
    print(f"1st/2nd ratio: {Dl_full[peaks[0]]/Dl_full[peaks[1]]:.2f} (Planck: ~2.5)")
    print(f"3rd/2nd ratio: {Dl_full[peaks[2]]/Dl_full[peaks[1]]:.2f} (Planck: ~1.0)")
print(f"\nExpected: l₁~220, l₂~540, l₃~810")
print(f"Our geometric: l₁~{np.pi*D_A/r_s:.0f}")

# Plot
fig, ax = plt.subplots(figsize=(12, 7))
ax.plot(l_full, Dl_full, 'b-', lw=1.5, label='Our solver (LCDM)')
for i, p in enumerate(peaks[:5]):
    color = ['red','orange','green','purple','brown'][i]
    ax.plot(l_full[p], Dl_full[p], 'o', color=color, ms=8)
    ax.annotate(f'l={l_full[p]}', (l_full[p], Dl_full[p]*1.08),
                fontsize=10, ha='center', color=color)

# Reference positions
for l_ref, label in [(220,'Planck 1st'), (540,'Planck 2nd'), (810,'Planck 3rd')]:
    ax.axvline(l_ref, color='gray', ls=':', alpha=0.4)

ax.set_xlabel('Multipole l', fontsize=13)
ax.set_ylabel(r'$D_\ell$ [$\mu K^2$]', fontsize=13)
ax.set_title('CMB TT — LCDM (Apple GPU / MLX, exact Bessel)', fontsize=14)
ax.legend(fontsize=12)
ax.set_xlim(2, 2500)
if max(Dl_full) > 0:
    ax.set_ylim(0, max(Dl_full)*1.3)
plt.tight_layout()
out = '/Users/akaihuangm1/Desktop/github/petz-recovery-unification/cmb_solver/lcdm_gpu_v2_result.png'
plt.savefig(out, dpi=150)
print(f"\nSaved: {out}")
print(f"Total: {time.time()-t_start:.2f}s")
