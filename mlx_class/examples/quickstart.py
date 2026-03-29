#!/usr/bin/env python3
"""
quickstart.py -- mlx_class quick-start demo.

Demonstrates the main capabilities of the mlx_class GPU-accelerated
CMB Boltzmann solver in a single, self-contained script.

What this script does:
  1. Solves the background cosmology (Friedmann + recombination)
  2. Runs the IMEX Boltzmann solver for 500 k-modes on the GPU
  3. Computes the CMB TT power spectrum C_l
  4. Finds and prints acoustic peak positions
  5. Computes the matter power spectrum P(k)
  6. Runs the Khronon dark matter model for comparison
  7. Displays total runtime

Requirements:
  pip install mlx numpy scipy matplotlib

Usage:
  python quickstart.py

Author: Sheng-Kai Huang, 2026
"""
import sys
import os
import time

# ---------------------------------------------------------------------------
# 0. Setup -- make sure mlx_class is importable
# ---------------------------------------------------------------------------
# If running from the examples/ directory, add the parent so that
# "import mlx_class" works without a pip install.
_parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _parent not in sys.path:
    sys.path.insert(0, os.path.dirname(_parent))

import numpy as np

# Check that MLX is available (Apple Silicon only)
try:
    import mlx.core as mx
    print(f"MLX device: {mx.default_device()}")
except ImportError:
    print("ERROR: Apple MLX is required.  Install with:  pip install mlx")
    print("       MLX only runs on macOS with Apple Silicon (M1/M2/M3/M4).")
    sys.exit(1)

t_total_start = time.time()

# ===========================================================================
# 1. Background cosmology
# ===========================================================================
print("\n" + "=" * 65)
print("  STEP 1: Background cosmology")
print("=" * 65)

from mlx_class.background import Background

t0 = time.time()
bg = Background(khronon=False, recombination='tanh')   # standard LCDM
bg.solve()
t_bg = time.time() - t0

print(f"\nBackground solved in {t_bg:.3f}s")
print(f"  Conformal time today   tau_0   = {bg.tau_0:.1f} Mpc")
print(f"  Conformal time at rec  tau_rec = {bg.tau_rec:.1f} Mpc")
print(f"  Comoving distance      D_A     = {bg.D_A:.1f} Mpc")
print(f"  Sound horizon          r_s     = {bg.r_s:.1f} Mpc")
print(f"  First peak estimate    l_1     ~ {np.pi * bg.D_A / bg.r_s:.0f}")

# ===========================================================================
# 2. IMEX Boltzmann solver
# ===========================================================================
print("\n" + "=" * 65)
print("  STEP 2: IMEX Boltzmann solver  (500 k-modes on GPU)")
print("=" * 65)

from mlx_class.perturbations_implicit import ImplicitBoltzmannSolver

# Logarithmically-spaced k-modes from 5e-5 to 0.35 Mpc^{-1}
k_arr = np.geomspace(5e-5, 0.35, 500).astype(np.float32)

t0 = time.time()
solver = ImplicitBoltzmannSolver(bg, k_arr)
result = solver.solve()
t_solve = time.time() - t0

# Extract source functions at recombination (with Silk damping)
Theta_0, Phi, v_b = result.source_at_recombination()
source_SW = Theta_0 + Phi    # Sachs-Wolfe monopole source

print(f"\nSolver finished in {t_solve:.2f}s")
print(f"  Source (Theta_0 + Phi) range: [{np.min(source_SW):.4f}, {np.max(source_SW):.4f}]")

# ===========================================================================
# 3. CMB TT power spectrum
# ===========================================================================
print("\n" + "=" * 65)
print("  STEP 3: CMB TT angular power spectrum C_l")
print("=" * 65)

from mlx_class.spectra import compute_cl
from mlx_class.main import upsample_source, early_isw_template

# The source S(k) is smooth, but j_l(k*D_A) oscillates fast.
# Upsample the source to a fine k-grid that resolves the Bessel oscillations.
k_fine, source_fine = upsample_source(k_arr, source_SW, bg.D_A, ell_max=2500)
print(f"Upsampled: {len(k_arr)} -> {len(k_fine)} k-points")

# Multipole values (sparse grid to keep runtime short)
ell_values = np.unique(np.concatenate([
    np.arange(2, 30, 1),
    np.arange(30, 100, 2),
    np.arange(100, 500, 4),
    np.arange(500, 1500, 8),
    np.arange(1500, 2501, 12),
])).astype(int)

t0 = time.time()
ell_out, Cl, Dl = compute_cl(source_fine, k_fine, ell_values, bg.D_A)
t_cl = time.time() - t0

# Add the phenomenological early-ISW template (improves first-peak position)
Dl_ISW = early_isw_template(ell_out, Dl, bg)
Dl = Dl + Dl_ISW

print(f"C_l computed in {t_cl:.2f}s  ({len(ell_values)} multipoles)")

# ===========================================================================
# 4. Find acoustic peaks
# ===========================================================================
print("\n" + "=" * 65)
print("  STEP 4: Acoustic peak positions")
print("=" * 65)

from scipy.signal import find_peaks
from scipy.interpolate import interp1d

# Interpolate to a dense ell grid and find peaks
l_dense = np.arange(2, 2501)
f_interp = interp1d(ell_out, Dl, kind='cubic', fill_value='extrapolate')
Dl_dense = np.maximum(f_interp(l_dense), 0.0)

peaks_idx, _ = find_peaks(Dl_dense, distance=150, prominence=100)

# Planck 2018 reference peak positions
planck_peaks = [220, 540, 810, 1120, 1420]

print(f"\n  {'Peak':<6} {'Found':>8} {'Planck':>8} {'Error':>10}")
print(f"  {'-' * 36}")
n_show = min(5, len(peaks_idx))
for i in range(n_show):
    l_found = l_dense[peaks_idx[i]]
    l_planck = planck_peaks[i]
    err_pct = (l_found - l_planck) / l_planck * 100
    print(f"  {i+1:<6} {l_found:>8} {l_planck:>8} {err_pct:>+9.1f}%")

if n_show >= 2:
    ratio_12 = Dl_dense[peaks_idx[0]] / Dl_dense[peaks_idx[1]]
    print(f"\n  1st/2nd peak height ratio: {ratio_12:.2f}  (Planck ~ 2.5)")

# ===========================================================================
# 5. Matter power spectrum P(k)
# ===========================================================================
print("\n" + "=" * 65)
print("  STEP 5: Matter power spectrum P(k) at z=0")
print("=" * 65)

from mlx_class.matter_pk import compute_matter_pk

t0 = time.time()
k_h, Pk_h, k_Mpc, Pk_Mpc = compute_matter_pk(
    result, solver_type='implicit', z_out=0.0)
t_pk = time.time() - t0

print(f"\nP(k) computed in {t_pk:.2f}s")
print(f"  k range:  [{k_h[0]:.2e}, {k_h[-1]:.2e}] h/Mpc")
print(f"  P(k) range: [{np.min(Pk_h):.2e}, {np.max(Pk_h):.2e}] (Mpc/h)^3")

# ===========================================================================
# 6. Khronon dark matter comparison
# ===========================================================================
print("\n" + "=" * 65)
print("  STEP 6: Khronon dark matter model comparison")
print("=" * 65)

from mlx_class.khronon import KhrononBoltzmannSolver

# Khronon background: Omega_K = 0.265 (replaces CDM)
bg_khr = Background(khronon=True)
bg_khr.solve()

# Solve with DBI sound speed
t0 = time.time()
solver_khr = KhrononBoltzmannSolver(bg_khr, k_arr, alpha_dbi=0.5)
result_khr = solver_khr.solve()
t_khr = time.time() - t0

Theta_0_k, Phi_k, v_b_k = result_khr.source_at_recombination()
source_khr = Theta_0_k + Phi_k

# Compute Khronon C_l
k_khr_fine, source_khr_fine = upsample_source(k_arr, source_khr, bg_khr.D_A)
ell_khr, Cl_khr, Dl_khr = compute_cl(source_khr_fine, k_khr_fine, ell_values, bg_khr.D_A)

print(f"\nKhronon solver: {t_khr:.2f}s")

# Compute residual
Dl_ref_safe = np.maximum(Dl, 1e-10)
residual_pct = (Dl_khr - Dl) / Dl_ref_safe * 100
mask_200_2000 = (ell_out >= 200) & (ell_out <= 2000)
if np.any(mask_200_2000):
    max_res = np.max(np.abs(residual_pct[mask_200_2000]))
    mean_res = np.mean(np.abs(residual_pct[mask_200_2000]))
    print(f"  Residual (200 < l < 2000):  max = {max_res:.2f}%,  mean = {mean_res:.2f}%")

# ===========================================================================
# 7. Plot everything
# ===========================================================================
print("\n" + "=" * 65)
print("  STEP 7: Generating plots")
print("=" * 65)

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # --- Panel (a): CMB TT spectrum ---
    ax = axes[0, 0]
    ax.plot(l_dense, Dl_dense, 'b-', lw=1.2, label=r'$\Lambda$CDM (IMEX)')
    for i in range(n_show):
        ax.plot(l_dense[peaks_idx[i]], Dl_dense[peaks_idx[i]],
                'ro', ms=7)
        ax.annotate(f'l={l_dense[peaks_idx[i]]}',
                    (l_dense[peaks_idx[i]], Dl_dense[peaks_idx[i]] * 1.06),
                    fontsize=8, ha='center', color='red')
    ax.set_xlabel('Multipole l')
    ax.set_ylabel(r'$D_\ell$ [$\mu K^2$]')
    ax.set_title('(a) CMB TT Power Spectrum')
    ax.legend(fontsize=10)
    ax.set_xlim(2, 2500)

    # --- Panel (b): Matter power spectrum ---
    ax = axes[0, 1]
    ax.loglog(k_h, Pk_h, 'b-', lw=1.5)
    ax.set_xlabel(r'$k$ [$h$ Mpc$^{-1}$]')
    ax.set_ylabel(r'$P(k)$ [(Mpc/$h$)$^3$]')
    ax.set_title('(b) Matter Power Spectrum P(k), z=0')

    # --- Panel (c): Khronon vs LCDM ---
    ax = axes[1, 0]
    f_khr = interp1d(ell_khr, Dl_khr, kind='cubic', fill_value='extrapolate')
    Dl_khr_dense = np.maximum(f_khr(l_dense), 0.0)
    ax.plot(l_dense, Dl_dense, 'b-', lw=1.2, label=r'$\Lambda$CDM')
    ax.plot(l_dense, Dl_khr_dense, 'r--', lw=1.2, label='Khronon (DBI)')
    ax.set_xlabel('Multipole l')
    ax.set_ylabel(r'$D_\ell$ [$\mu K^2$]')
    ax.set_title(r'(c) Khronon vs $\Lambda$CDM')
    ax.legend(fontsize=10)
    ax.set_xlim(2, 2500)

    # --- Panel (d): Residual ---
    ax = axes[1, 1]
    ax.plot(ell_out, residual_pct, 'k-', lw=0.8)
    ax.axhline(0, color='gray', ls='-', alpha=0.3)
    ax.axhline(1, color='red', ls=':', alpha=0.4)
    ax.axhline(-1, color='red', ls=':', alpha=0.4)
    ax.set_xlabel('Multipole l')
    ax.set_ylabel('Residual [%]')
    ax.set_title('(d) Khronon Residual (C-A)/A')
    ax.set_xlim(2, 2500)

    plt.tight_layout()
    out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            'quickstart_output.png')
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"  Plot saved: {out_path}")

except ImportError:
    print("  matplotlib not installed -- skipping plots.")
    print("  Install with:  pip install matplotlib")

# ===========================================================================
# Summary
# ===========================================================================
t_total = time.time() - t_total_start

print("\n" + "=" * 65)
print("  SUMMARY")
print("=" * 65)
print(f"  Background:        {t_bg:.3f}s")
print(f"  IMEX solver:       {t_solve:.2f}s")
print(f"  C_l integration:   {t_cl:.2f}s")
print(f"  P(k):              {t_pk:.2f}s")
print(f"  Khronon solver:    {t_khr:.2f}s")
print(f"  ---------------------------------")
print(f"  TOTAL:             {t_total:.2f}s")
print(f"\n  GPU:  {mx.default_device()}")
print("=" * 65)
