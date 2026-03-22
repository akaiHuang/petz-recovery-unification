# DEPRECATED: Use tests/test_basic.py instead.
# This file is kept for backward compatibility only.
import warnings as _warnings
_warnings.warn(
    "mlx_class.test_combined is deprecated. Use mlx_class.tests.test_basic instead.",
    DeprecationWarning,
    stacklevel=2,
)

#!/usr/bin/env python3
"""
test_combined.py — Combined neutrino + ISW accuracy test.

Pipeline:
  1. Run NeutrinoBoltzmannSolver (with N_eff=3.046 massless neutrinos)
  2. Extract source at recombination (Theta_0 + Psi with Silk damping)
  3. Upsample to fine k-grid (resolves j_l oscillations)
  4. Compute C_l with Bessel j_l^2 projection
  5. Add early ISW template (matter-radiation transition)
  6. Find peaks and compare with Planck (220, 540, 810, 1120, 1420)
  7. Compare with no-neutrino IMEX result to quantify improvement

Author: Sheng-Kai Huang, 2026
"""
import os
import sys
import time
import numpy as np

# Ensure package import works
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from mlx_class.background import Background, T_CMB
from mlx_class.perturbations_neutrino import NeutrinoBoltzmannSolver, L_NU_MAX
from mlx_class.perturbations_implicit import ImplicitBoltzmannSolver
from mlx_class.spectra import compute_cl
from mlx_class.main import upsample_source, early_isw_template


# ============================================================================
# Planck reference peaks (TT power spectrum)
# ============================================================================
PLANCK_PEAKS = np.array([220, 540, 810, 1120, 1420])
PLANCK_RATIO_12 = 2.5  # approximate 1st/2nd peak height ratio


def find_peaks_in_dl(ell_values, Dl, ell_max=2500, l_min=100):
    """
    Find acoustic peaks in D_l spectrum by interpolating to dense ell grid.

    Only searches for peaks at l > l_min to avoid the ISW bump at low ell.
    The CMB acoustic peaks are at l ~ 220, 540, 810, 1120, 1420.
    """
    from scipy.signal import find_peaks
    from scipy.interpolate import interp1d

    l_full = np.arange(2, ell_max + 1)
    f = interp1d(ell_values, Dl, kind='cubic', fill_value='extrapolate')
    Dl_full = np.maximum(f(l_full), 0.0)

    # Only search for peaks above l_min (skip ISW bump)
    mask = l_full >= l_min
    l_search = l_full[mask]
    Dl_search = Dl_full[mask]

    # Find peaks with constraints appropriate for CMB TT acoustic peaks
    peaks_idx, props = find_peaks(Dl_search, distance=120, prominence=20)

    peak_ells = l_search[peaks_idx]
    peak_heights = Dl_search[peaks_idx]

    return peak_ells, peak_heights, l_full, Dl_full


def print_peak_table(label, peak_ells, peak_heights, planck=PLANCK_PEAKS):
    """Print a formatted peak comparison table."""
    n = min(len(peak_ells), len(planck))
    print(f"\n{'='*70}")
    print(f"  {label}")
    print(f"{'='*70}")
    print(f"  {'Peak':<6} {'Planck':>8} {'Found':>8} {'Error':>10} {'Height':>12}")
    print(f"  {'-'*50}")
    for i in range(n):
        err = (peak_ells[i] - planck[i]) / planck[i] * 100
        print(f"  {i+1:<6} {planck[i]:>8} {peak_ells[i]:>8} {err:>+9.1f}% "
              f"{peak_heights[i]:>10.0f} uK^2")
    if n >= 2:
        ratio = peak_heights[0] / peak_heights[1]
        print(f"  {'-'*50}")
        print(f"  1st/2nd ratio: {ratio:.3f}  (Planck ~ {PLANCK_RATIO_12})")
    print(f"{'='*70}")
    return n


def run_solver(label, solver_class, bg, k_arr, use_psi=False, **kwargs):
    """Run a solver and return (ell, Dl_with_ISW, Dl_no_ISW, time_elapsed)."""
    print(f"\n{'#'*70}")
    print(f"# {label}")
    print(f"{'#'*70}")

    t0 = time.time()

    # Create and run solver
    solver = solver_class(bg, k_arr, **kwargs)
    result = solver.solve()
    t_solve = time.time() - t0

    # Extract source
    if use_psi:
        # Neutrino solver: Theta_0 + Psi (proper SW source with anisotropic stress)
        Theta_0, Phi, Psi, v_b, N_0, N_2 = result.source_at_recombination()
        source_SW = Theta_0 + Psi
        print(f"[{label}] Using Theta_0 + Psi (neutrino anisotropic stress)")
        print(f"[{label}] Psi/Phi superhorizon: {np.mean(result.phi_psi_ratio()[:10]):.4f}")
    else:
        # No-neutrino solver: Theta_0 + Phi
        Theta_0, Phi, v_b = result.source_at_recombination()
        source_SW = Theta_0 + Phi
        print(f"[{label}] Using Theta_0 + Phi (no anisotropic stress)")

    print(f"[{label}] Source range: [{np.min(source_SW):.4f}, {np.max(source_SW):.4f}]")

    # Upsample to fine k-grid
    k_fine, source_fine = upsample_source(k_arr, source_SW, bg.D_A, ell_max=2500)
    print(f"[{label}] Upsampled: {len(k_arr)} -> {len(k_fine)} k-points "
          f"(dk={k_fine[1]-k_fine[0]:.2e})")

    # Build ell grid
    ell_values = np.unique(np.concatenate([
        np.arange(2, 30, 1),
        np.arange(30, 100, 2),
        np.arange(100, 500, 4),
        np.arange(500, 1500, 8),
        np.arange(1500, 2501, 12),
    ])).astype(int)

    # Compute C_l (Bessel j_l^2 projection)
    ell_out, Cl, Dl_acoustic = compute_cl(source_fine, k_fine, ell_values, bg.D_A)

    # Add early ISW template
    Dl_ISW = early_isw_template(ell_out, Dl_acoustic, bg)
    Dl_total = Dl_acoustic + Dl_ISW
    print(f"[{label}] Early ISW peak: {np.max(Dl_ISW):.0f} uK^2")

    t_total = time.time() - t0
    print(f"[{label}] Total time: {t_total:.2f}s (solve: {t_solve:.2f}s)")

    return ell_out, Dl_total, Dl_acoustic, t_total, result


def main():
    t_start = time.time()

    print("=" * 70)
    print("COMBINED TEST: Neutrino IMEX + ISW + Bessel Projection")
    print("=" * 70)

    # ---- Step 1: Background ----
    print("\n--- Step 1: Background cosmology ---")
    t0 = time.time()
    bg = Background(khronon=False)
    bg.solve()
    t_bg = time.time() - t0
    print(f"Background: {t_bg:.2f}s")

    # ---- k-grid ----
    k_arr = np.geomspace(5e-4, 0.35, 500).astype(np.float32)
    print(f"k-grid: [{k_arr[0]:.4f}, {k_arr[-1]:.4f}] Mpc^-1, N_k={len(k_arr)}")

    # ---- Step 2: Run WITH neutrinos ----
    ell_nu, Dl_nu, Dl_nu_acoustic, t_nu, result_nu = run_solver(
        "Neutrino IMEX (N_eff=3.046)",
        NeutrinoBoltzmannSolver, bg, k_arr,
        use_psi=True, l_nu_max=L_NU_MAX)

    # ---- Step 3: Run WITHOUT neutrinos (reference) ----
    ell_ref, Dl_ref, Dl_ref_acoustic, t_ref, result_ref = run_solver(
        "No-neutrino IMEX (reference)",
        ImplicitBoltzmannSolver, bg, k_arr,
        use_psi=False)

    # ---- Step 4: Find peaks ----
    peaks_nu, heights_nu, l_full, Dl_nu_full = find_peaks_in_dl(ell_nu, Dl_nu)
    peaks_ref, heights_ref, _, Dl_ref_full = find_peaks_in_dl(ell_ref, Dl_ref)

    n_nu = print_peak_table("WITH Neutrinos + ISW", peaks_nu, heights_nu)
    n_ref = print_peak_table("WITHOUT Neutrinos + ISW (reference)", peaks_ref, heights_ref)

    # ---- Step 5: Quantify neutrino improvement ----
    print(f"\n{'='*70}")
    print("  NEUTRINO IMPROVEMENT ANALYSIS")
    print(f"{'='*70}")

    n_compare = min(len(peaks_nu), len(peaks_ref), len(PLANCK_PEAKS))
    if n_compare >= 1:
        err_nu = np.abs(peaks_nu[:n_compare] - PLANCK_PEAKS[:n_compare]) / PLANCK_PEAKS[:n_compare] * 100
        err_ref = np.abs(peaks_ref[:n_compare] - PLANCK_PEAKS[:n_compare]) / PLANCK_PEAKS[:n_compare] * 100

        print(f"\n  Mean peak position error:")
        print(f"    With neutrinos:    {np.mean(err_nu):.1f}%")
        print(f"    Without neutrinos: {np.mean(err_ref):.1f}%")
        if np.mean(err_ref) > 0:
            improvement = (np.mean(err_ref) - np.mean(err_nu)) / np.mean(err_ref) * 100
            print(f"    Improvement:       {improvement:+.1f}%")

    if n_compare >= 2:
        ratio_nu = heights_nu[0] / heights_nu[1]
        ratio_ref = heights_ref[0] / heights_ref[1]
        print(f"\n  1st/2nd peak height ratio:")
        print(f"    With neutrinos:    {ratio_nu:.3f}  (Planck ~ {PLANCK_RATIO_12})")
        print(f"    Without neutrinos: {ratio_ref:.3f}")
        print(f"    Planck target:     {PLANCK_RATIO_12:.3f}")

    # Neutrino effect on SW plateau (low ell)
    mask_sw = ell_nu < 30
    if np.sum(mask_sw) > 0:
        sw_ratio = np.mean(Dl_nu[mask_sw]) / np.mean(Dl_ref[mask_sw]) if np.mean(Dl_ref[mask_sw]) > 0 else 1.0
        print(f"\n  SW plateau (l<30): neutrino/reference = {sw_ratio:.4f} ({(sw_ratio-1)*100:+.1f}%)")

    # ---- Step 6: Timing summary ----
    t_total = time.time() - t_start
    print(f"\n{'='*70}")
    print("  TIMING SUMMARY")
    print(f"{'='*70}")
    print(f"  Background:       {t_bg:.2f}s")
    print(f"  Neutrino solver:  {t_nu:.2f}s")
    print(f"  Reference solver: {t_ref:.2f}s")
    print(f"  Total:            {t_total:.2f}s")
    print(f"{'='*70}")

    # ---- Step 7: Generate plot ----
    print("\n--- Generating combined plot ---")
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from scipy.interpolate import interp1d

        fig, axes = plt.subplots(2, 2, figsize=(16, 12))

        # --- Panel 1: Full D_l spectra ---
        ax = axes[0, 0]
        ax.plot(l_full, Dl_nu_full, 'b-', lw=1.5,
                label=r'With $\nu$ ($N_{eff}=3.046$) + ISW')
        ax.plot(l_full, Dl_ref_full, 'r--', lw=1.2, alpha=0.7,
                label=r'Without $\nu$ + ISW')

        # Mark neutrino peaks
        colors = ['blue', 'dodgerblue', 'steelblue', 'slateblue', 'navy']
        for i, (lp, hp) in enumerate(zip(peaks_nu[:5], heights_nu[:5])):
            ax.plot(lp, hp, 'o', color=colors[min(i, len(colors)-1)], ms=9, zorder=5)
            ax.annotate(f'l={lp}', (lp, hp * 1.06),
                        fontsize=9, ha='center', color=colors[min(i, len(colors)-1)],
                        fontweight='bold')

        # Planck reference lines
        for lp in PLANCK_PEAKS:
            ax.axvline(lp, color='green', ls=':', alpha=0.35, lw=1)
        ax.axvline(PLANCK_PEAKS[0], color='green', ls=':', alpha=0.35, lw=1,
                   label='Planck peaks')

        ax.set_xlabel('Multipole l', fontsize=12)
        ax.set_ylabel(r'$D_\ell$ [$\mu K^2$]', fontsize=12)
        ax.set_title('CMB TT Power Spectrum', fontsize=13)
        ax.legend(fontsize=10, loc='upper right')
        ax.set_xlim(2, 2500)
        if np.max(Dl_nu_full) > 0:
            ax.set_ylim(0, np.max(Dl_nu_full) * 1.3)

        # --- Panel 2: Neutrino effect ratio ---
        ax = axes[0, 1]
        mask_pos = Dl_ref_full > 1.0
        ratio_cl = np.ones_like(Dl_nu_full)
        ratio_cl[mask_pos] = Dl_nu_full[mask_pos] / Dl_ref_full[mask_pos]
        ax.plot(l_full, ratio_cl, 'k-', lw=0.8, alpha=0.7)
        ax.axhline(1.0, color='gray', ls='-', alpha=0.4)
        ax.fill_between(l_full, 0.95, 1.05, color='green', alpha=0.1)
        ax.set_xlabel('Multipole l', fontsize=12)
        ax.set_ylabel(r'$D_\ell^{\nu} / D_\ell^{no\,\nu}$', fontsize=12)
        ax.set_title('Neutrino Effect on Power Spectrum', fontsize=13)
        ax.set_xlim(2, 2500)
        ax.set_ylim(0.5, 1.5)

        # --- Panel 3: ISW decomposition ---
        ax = axes[1, 0]
        # Recompute ISW for plotting
        Dl_ISW_nu = early_isw_template(ell_nu, Dl_nu_acoustic, bg)
        f_ac = interp1d(ell_nu, Dl_nu_acoustic, kind='cubic', fill_value='extrapolate')
        f_isw = interp1d(ell_nu, Dl_ISW_nu, kind='cubic', fill_value='extrapolate')
        Dl_ac_full = np.maximum(f_ac(l_full), 0.0)
        Dl_isw_full = np.maximum(f_isw(l_full), 0.0)

        ax.plot(l_full, Dl_nu_full, 'b-', lw=1.5, label='Total (acoustic + ISW)')
        ax.plot(l_full, Dl_ac_full, 'r--', lw=1.2, alpha=0.7, label='Acoustic only')
        ax.plot(l_full, Dl_isw_full, 'g:', lw=1.5, alpha=0.8, label='Early ISW template')
        ax.set_xlabel('Multipole l', fontsize=12)
        ax.set_ylabel(r'$D_\ell$ [$\mu K^2$]', fontsize=12)
        ax.set_title('ISW Decomposition (with neutrinos)', fontsize=13)
        ax.legend(fontsize=10)
        ax.set_xlim(2, 1500)
        if np.max(Dl_nu_full[:1498]) > 0:
            ax.set_ylim(0, np.max(Dl_nu_full[:1498]) * 1.3)

        # --- Panel 4: Peak comparison bar chart ---
        ax = axes[1, 1]
        n_plot = min(5, len(peaks_nu), len(peaks_ref))
        if n_plot >= 1:
            x = np.arange(n_plot)
            width = 0.25
            bars_planck = ax.bar(x - width, PLANCK_PEAKS[:n_plot], width,
                                  color='green', alpha=0.6, label='Planck')
            bars_nu = ax.bar(x, peaks_nu[:n_plot], width,
                              color='blue', alpha=0.6, label='With neutrinos')
            bars_ref = ax.bar(x + width, peaks_ref[:n_plot], width,
                               color='red', alpha=0.6, label='Without neutrinos')

            # Add error annotations
            for i in range(n_plot):
                err = (peaks_nu[i] - PLANCK_PEAKS[i]) / PLANCK_PEAKS[i] * 100
                ax.annotate(f'{err:+.1f}%',
                            (x[i], peaks_nu[i] + 20),
                            fontsize=8, ha='center', color='blue')

            ax.set_xticks(x)
            ax.set_xticklabels([f'Peak {i+1}' for i in range(n_plot)])
            ax.set_ylabel('Multipole l', fontsize=12)
            ax.set_title('Peak Positions vs Planck', fontsize=13)
            ax.legend(fontsize=9)

        plt.suptitle(f'Combined Neutrino + ISW Test  (total: {t_total:.1f}s)',
                     fontsize=15, y=1.01)
        plt.tight_layout()

        out_dir = os.path.dirname(os.path.abspath(__file__))
        out_path = os.path.join(out_dir, 'combined_neutrino_isw.png')
        plt.savefig(out_path, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"Plot saved: {out_path}")

    except ImportError as e:
        print(f"matplotlib not available: {e}")

    # ---- Final summary ----
    print(f"\n{'='*70}")
    print("  FINAL SUMMARY")
    print(f"{'='*70}")

    if len(peaks_nu) >= 5:
        for i in range(min(5, len(peaks_nu))):
            err = (peaks_nu[i] - PLANCK_PEAKS[i]) / PLANCK_PEAKS[i] * 100
            status = "OK" if abs(err) < 15 else "WARN"
            print(f"  Peak {i+1}: l={peaks_nu[i]:>5} (Planck {PLANCK_PEAKS[i]:>5}, "
                  f"err={err:+.1f}%) [{status}]")
    elif len(peaks_nu) >= 1:
        for i in range(len(peaks_nu)):
            if i < len(PLANCK_PEAKS):
                err = (peaks_nu[i] - PLANCK_PEAKS[i]) / PLANCK_PEAKS[i] * 100
                print(f"  Peak {i+1}: l={peaks_nu[i]:>5} (Planck {PLANCK_PEAKS[i]:>5}, "
                      f"err={err:+.1f}%)")
        print(f"  WARNING: Only {len(peaks_nu)} peaks found (expected 5)")

    if len(peaks_nu) >= 2:
        ratio = heights_nu[0] / heights_nu[1]
        print(f"  1st/2nd ratio: {ratio:.3f} (Planck ~ {PLANCK_RATIO_12})")

    print(f"  Total time: {t_total:.1f}s")
    print(f"{'='*70}")

    return 0


if __name__ == '__main__':
    sys.exit(main())
