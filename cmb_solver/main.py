#!/usr/bin/env python3
"""
main.py — CMB TT Power Spectrum from Sigma = 2 ln Q (Khronon cosmology).

Main driver:
1. Background cosmology (Friedmann + ghost condensation Khronon)
2. Recombination (tanh fit calibrated to RECFAST)
3. Perturbation evolution (two-phase: tight coupling then free streaming)
4. Line-of-sight integration for C_l
5. Comparison with Planck 2018 TT data

Usage:
    python -m cmb_solver.main [--fast] [--ell_max 1500] [--n_k 50] [--output DIR]

Theory: Sigma = 2 ln Q, ghost-condensation Khronon, c_s^2 = 0.
Background is EXACTLY CDM; differences from LCDM appear in perturbation sector.

Author: Sheng-Kai Huang, 2026
"""
import argparse
import os
import sys
import time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from cmb_solver import constants as C
from cmb_solver.background import BackgroundSolver
from cmb_solver.power_spectrum import PowerSpectrumCalculator, compute_cl_fast
from cmb_solver.compare import (
    plot_comparison, plot_residuals, compute_chi2, load_planck_data
)


def main():
    parser = argparse.ArgumentParser(
        description='CMB TT power spectrum from Sigma = 2 ln Q')
    parser.add_argument('--fast', action='store_true',
                        help='Fast mode: fewer k-modes, lower resolution')
    parser.add_argument('--ell_max', type=int, default=1500,
                        help='Maximum multipole (default: 1500)')
    parser.add_argument('--n_k', type=int, default=50,
                        help='Number of k-modes (default: 50)')
    parser.add_argument('--n_tau', type=int, default=300,
                        help='Number of tau points per k (default: 300)')
    parser.add_argument('--l_max', type=int, default=8,
                        help='Boltzmann hierarchy l_max (default: 8)')
    parser.add_argument('--output', type=str, default='output',
                        help='Output directory (default: output)')
    parser.add_argument('--no_plot', action='store_true',
                        help='Skip plotting')
    args = parser.parse_args()

    if args.fast:
        args.n_k = 30
        args.n_tau = 200
        args.l_max = 6

    output_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), args.output)
    os.makedirs(output_dir, exist_ok=True)

    # ===================================================================
    C.print_params()
    t_total = time.time()

    # Step 1: Background
    print("\n" + "=" * 60)
    print("STEP 1: Background cosmology")
    print("=" * 60)
    t0 = time.time()
    bg = BackgroundSolver(N_points=15000)
    bg.solve()
    print(f"  Time: {time.time()-t0:.1f} s")

    # Step 2: Recombination
    print("\n" + "=" * 60)
    print("STEP 2: Recombination")
    print("=" * 60)
    t0 = time.time()
    bg.compute_recombination()
    print(f"  Time: {time.time()-t0:.1f} s")

    # Step 3: Background plots
    if not args.no_plot:
        plot_background(bg, output_dir)

    # Step 4: C_l
    print("\n" + "=" * 60)
    print("STEP 3: CMB TT power spectrum")
    print("=" * 60)
    t0 = time.time()

    # Multipoles
    ell_values = np.unique(np.concatenate([
        np.arange(2, min(30, args.ell_max + 1), 2),
        np.arange(30, min(100, args.ell_max + 1), 5),
        np.arange(100, min(500, args.ell_max + 1), 10),
        np.arange(500, args.ell_max + 1, 20),
    ])).astype(int)

    print(f"[C_l] {len(ell_values)} multipoles: {ell_values[0]} to {ell_values[-1]}")

    ell_out, Cl, Dl = compute_cl_fast(
        bg, ell_values,
        n_k=args.n_k,
        n_tau=args.n_tau,
        l_max_hierarchy=args.l_max,
        progress=True
    )

    print(f"\n  Time: {time.time()-t0:.0f} s")

    # Save
    result_file = os.path.join(output_dir, 'cl_khronon.dat')
    np.savetxt(result_file,
               np.column_stack([ell_out, Cl, Dl]),
               header='ell  C_l  D_l[muK^2]',
               fmt='%6d  %.6e  %.6e')
    print(f"[Output] Saved to {result_file}")

    # Step 5: Compare
    print("\n" + "=" * 60)
    print("STEP 4: Comparison with Planck 2018")
    print("=" * 60)

    # Only compare ell >= 30 (low ell ISW needs more careful treatment)
    chi2, ndof, chi2_dof = compute_chi2(ell_out, Dl, ell_range=(30, 2500))
    print(f"  chi^2 = {chi2:.1f} (ell >= 30)")
    print(f"  N_dof = {ndof}")
    print(f"  chi^2/dof = {chi2_dof:.2f}")

    if not args.no_plot:
        fig_path = os.path.join(output_dir, 'cl_comparison.png')
        plot_comparison(ell_out, Dl, save_path=fig_path, show=False)
        res_path = os.path.join(output_dir, 'cl_residuals.png')
        plot_residuals(ell_out, Dl, save_path=res_path, show=False)

    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"  Theory: Sigma = 2 ln Q (Khronon, c_s^2 = 0)")
    print(f"  Background: identical to LCDM (ghost condensation)")
    print(f"  Omega_K = {C.Omega_K_khronon}, Q_0 = {bg.Q0:.4f}")
    print(f"  Sigma_0 = {bg.Sigma_bg:.4f}")
    print(f"  ell range: {ell_out[0]} -- {ell_out[-1]}")
    print(f"  chi^2/dof (ell>=30) = {chi2_dof:.2f}")
    print(f"  Total time: {time.time()-t_total:.0f} s")
    print(f"  Output: {output_dir}/")
    print("=" * 60)


def plot_background(bg, output_dir):
    """Plot background quantities."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    z_grid = 1.0 / bg.a_grid - 1

    # Hubble
    ax = axes[0, 0]
    ax.loglog(1 + z_grid, bg.H_grid, 'b-')
    ax.set_xlabel('1+z')
    ax.set_ylabel('H(z) / H0')
    ax.set_title('Hubble parameter')
    ax.grid(True, alpha=0.3)

    # Density fractions
    ax = axes[0, 1]
    H2 = bg.H_grid**2
    Omega_r_z = C.Omega_r * bg.a_grid**(-4) / H2
    Omega_b_z = C.Omega_b * bg.a_grid**(-3) / H2
    Omega_K_z = C.Omega_K_khronon * bg.a_grid**(-3) / H2
    Omega_L_z = C.Omega_Lambda / H2

    ax.semilogx(1 + z_grid, Omega_r_z, 'y-', label=r'$\Omega_r$')
    ax.semilogx(1 + z_grid, Omega_b_z, 'g-', label=r'$\Omega_b$')
    ax.semilogx(1 + z_grid, Omega_K_z, 'r-', label=r'$\Omega_K$ (Khronon)')
    ax.semilogx(1 + z_grid, Omega_L_z, 'b-', label=r'$\Omega_\Lambda$')
    ax.set_xlabel('1+z')
    ax.set_ylabel(r'$\Omega_i(z)$')
    ax.set_title('Density fractions')
    ax.set_ylim(0, 1.1)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # Ionization fraction
    ax = axes[1, 0]
    ax.semilogx(1 + z_grid, bg.x_e_grid, 'k-')
    ax.axvline(1 + C.z_rec, color='red', ls='--', alpha=0.5, label=f'z={C.z_rec:.0f}')
    ax.set_xlabel('1+z')
    ax.set_ylabel(r'$x_e$')
    ax.set_title('Ionization fraction')
    ax.set_xlim(100, 3000)
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Visibility
    ax = axes[1, 1]
    vis = bg.visibility / np.max(bg.visibility) if np.max(bg.visibility) > 0 else bg.visibility
    ax.semilogx(1 + z_grid, vis, 'k-')
    ax.axvline(1 + C.z_rec, color='red', ls='--', alpha=0.5)
    ax.set_xlabel('1+z')
    ax.set_ylabel('Visibility (normalized)')
    ax.set_title('Visibility function')
    ax.set_xlim(800, 2000)
    ax.grid(True, alpha=0.3)

    plt.suptitle(r'Khronon Background: $\Sigma = 2\ln Q$, Ghost Condensation', fontsize=14)
    plt.tight_layout()
    fig_path = os.path.join(output_dir, 'background.png')
    plt.savefig(fig_path, dpi=150, bbox_inches='tight')
    print(f"[Plot] Background saved to {fig_path}")
    plt.close()


if __name__ == '__main__':
    main()
