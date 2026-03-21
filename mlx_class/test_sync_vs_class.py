"""
test_sync_vs_class.py -- Compare synchronous gauge solver against CLASS reference.

Loads CLASS reference C_l from disk, runs the sync gauge solver, and prints
a detailed comparison: RMS error, peak positions, D_l at key multipoles.

Usage:
  cd /Users/akaihuangm1/Desktop/github/petz-recovery-unification
  python -m mlx_class.test_sync_vs_class

Author: Sheng-Kai Huang, 2026
"""
import os
os.environ['MPLCONFIGDIR'] = '/tmp/mpl_config'

import numpy as np
import time
import sys

from .background import T_CMB
from .solver_sync import run_sync_solver


# ============================================================================
# Load CLASS reference
# ============================================================================

CLASS_FILE = '/Users/akaihuangm1/Desktop/github/gdm_class_public/output/mlx_ref_lcdm_cl.dat'


def load_class_cl(path=CLASS_FILE):
    """Load CLASS cl.dat: dimensionless l(l+1)/(2pi) C_l -> D_l [muK^2]."""
    data = np.loadtxt(path)
    ell = data[:, 0].astype(int)
    # CLASS outputs l(l+1)/(2pi) C_l in dimensionless units
    Dl_muK2 = data[:, 1] * (T_CMB * 1e6) ** 2
    return ell, Dl_muK2


# ============================================================================
# Comparison metrics
# ============================================================================

def compute_rms(ell_ref, Dl_ref, ell_test, Dl_test, ell_min=2, ell_max=2500):
    """Compute RMS relative error between test and reference D_l."""
    from scipy.interpolate import interp1d

    # Interpolate test onto reference ell grid
    mask = (ell_ref >= ell_min) & (ell_ref <= ell_max)
    ell_use = ell_ref[mask]
    Dl_ref_use = Dl_ref[mask]

    f_test = interp1d(ell_test, Dl_test, kind='cubic', fill_value='extrapolate')
    Dl_test_interp = f_test(ell_use)

    # Relative error (avoid division by near-zero)
    Dl_max = np.max(np.abs(Dl_ref_use))
    rel_err = (Dl_test_interp - Dl_ref_use) / Dl_max

    rms = np.sqrt(np.mean(rel_err ** 2))
    return rms, ell_use, Dl_ref_use, Dl_test_interp


def find_peaks(ell, Dl, n_peaks=7, smooth_sigma=15):
    """Find acoustic peak positions and heights in D_l spectrum."""
    from scipy.ndimage import gaussian_filter1d
    from scipy.signal import find_peaks as sp_find_peaks
    from scipy.interpolate import interp1d

    ell_fine = np.arange(int(ell[0]), int(ell[-1]) + 1)
    f = interp1d(ell, Dl, kind='cubic', fill_value='extrapolate')
    Dl_fine = np.maximum(f(ell_fine), 0.0)

    Dl_smooth = gaussian_filter1d(Dl_fine, sigma=smooth_sigma)
    peaks, _ = sp_find_peaks(Dl_smooth, distance=80, prominence=50)

    result = []
    for p in peaks[:n_peaks]:
        result.append((int(ell_fine[p]), float(Dl_smooth[p])))
    return result


# ============================================================================
# Main comparison
# ============================================================================

def main():
    print("=" * 70)
    print("SYNCHRONOUS GAUGE SOLVER vs CLASS COMPARISON")
    print("=" * 70)
    print()

    # Load CLASS reference
    print("--- Loading CLASS reference ---")
    if not os.path.exists(CLASS_FILE):
        print(f"ERROR: CLASS file not found: {CLASS_FILE}")
        sys.exit(1)
    class_ell, class_Dl = load_class_cl()
    print(f"  Loaded {len(class_ell)} multipoles (l={class_ell[0]}..{class_ell[-1]})")
    class_peaks = find_peaks(class_ell, class_Dl)
    print(f"  CLASS peaks: {[p[0] for p in class_peaks[:5]]}")
    print()

    # Run synchronous gauge solver
    # Use default lg_max and ln_max from perturbations_sync.py
    result = run_sync_solver(
        N_k=180, k_min=3e-4, k_max=0.35,
        method='Radau',
        rtol=1e-6, atol=1e-9,
        verbose=True)

    ell = result['ell']
    Dl = result['Dl']

    # ================================================================
    # Comparison
    # ================================================================
    print("\n" + "=" * 70)
    print("COMPARISON RESULTS")
    print("=" * 70)

    # RMS
    rms_full, ell_use, Dl_ref_use, Dl_test_interp = compute_rms(
        class_ell, class_Dl, ell, Dl, ell_min=2, ell_max=2500)

    rms_low, _, _, _ = compute_rms(
        class_ell, class_Dl, ell, Dl, ell_min=2, ell_max=100)

    rms_peaks, _, _, _ = compute_rms(
        class_ell, class_Dl, ell, Dl, ell_min=100, ell_max=1500)

    rms_tail, _, _, _ = compute_rms(
        class_ell, class_Dl, ell, Dl, ell_min=1500, ell_max=2500)

    print(f"\nRMS relative error (normalized to peak D_l):")
    print(f"  Full range (2-2500):  {rms_full*100:.2f}%")
    print(f"  Low ell (2-100):      {rms_low*100:.2f}%")
    print(f"  Peak region (100-1500): {rms_peaks*100:.2f}%")
    print(f"  Damping tail (1500-2500): {rms_tail*100:.2f}%")

    # Peak comparison
    print(f"\n--- Peak positions ---")
    print(f"{'Peak':>6s}  {'CLASS l':>8s}  {'Sync l':>8s}  {'Err':>7s}  "
          f"{'CLASS D_l':>10s}  {'Sync D_l':>10s}  {'Amp Ratio':>10s}")
    print("-" * 70)

    sync_peaks = find_peaks(ell, Dl)
    n_compare = min(len(class_peaks), len(sync_peaks), 5)
    for i in range(n_compare):
        c_l, c_d = class_peaks[i]
        s_l, s_d = sync_peaks[i]
        err = s_l - c_l
        ratio = s_d / c_d if c_d > 0 else 0.0
        print(f"  {i+1:4d}  {c_l:8d}  {s_l:8d}  {err:+5d}  "
              f"{c_d:10.1f}  {s_d:10.1f}  {ratio:10.3f}")

    # D_l at key multipoles
    print(f"\n--- D_l at key multipoles ---")
    print(f"{'ell':>6s}  {'CLASS':>12s}  {'Sync':>12s}  {'Ratio':>8s}  {'Err%':>8s}")
    print("-" * 55)

    from scipy.interpolate import interp1d
    f_sync = interp1d(ell, Dl, kind='cubic', fill_value='extrapolate')
    f_class = interp1d(class_ell, class_Dl, kind='cubic', fill_value='extrapolate')

    for target in [2, 10, 50, 100, 220, 400, 537, 700, 820, 1000, 1200, 1500, 2000]:
        d_c = float(f_class(target))
        d_s = float(f_sync(target))
        ratio = d_s / d_c if abs(d_c) > 1e-6 else 0.0
        err = (d_s - d_c) / d_c * 100 if abs(d_c) > 1e-6 else 0.0
        print(f"  {target:5d}  {d_c:12.2f}  {d_s:12.2f}  {ratio:8.3f}  {err:+7.1f}%")

    # Timing summary
    print(f"\n--- Timing ---")
    t = result['timing']
    print(f"  Background:    {t['background']:.2f}s")
    print(f"  Perturbations: {t['perturbations']:.1f}s")
    if 'gauge' in t:
        print(f"  Gauge xform:   {t['gauge']:.2f}s")
    print(f"  LOS integral:  {t['los']:.2f}s")
    print(f"  TOTAL:         {t['total']:.1f}s")

    # Save output
    out_dir = os.path.dirname(os.path.abspath(__file__))
    out_file = os.path.join(out_dir, 'cl_sync_gauge.dat')
    np.savetxt(out_file, np.column_stack([ell, Dl]),
               header='l   D_l[muK^2]  (sync gauge solver)')
    print(f"\n  Saved: {out_file}")

    # Try to save a comparison plot
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(2, 1, figsize=(12, 8), height_ratios=[3, 1],
                                  sharex=True)

        # Top: spectra
        ax = axes[0]
        ax.plot(class_ell, class_Dl, 'k-', lw=1.5, label='CLASS', alpha=0.8)
        ax.plot(ell, Dl, 'r-', lw=1.0, label='Sync gauge', alpha=0.8)
        ax.set_ylabel(r'$D_\ell$ [$\mu$K$^2$]')
        ax.set_title(f'Synchronous Gauge Solver vs CLASS  (RMS={rms_full*100:.1f}%)')
        ax.legend()
        ax.set_xlim(2, 2500)

        # Bottom: residual
        ax = axes[1]
        Dl_class_interp = f_class(ell)
        residual = (Dl - Dl_class_interp) / np.max(np.abs(class_Dl)) * 100
        ax.plot(ell, residual, 'b-', lw=0.8)
        ax.axhline(0, color='gray', ls='--', lw=0.5)
        ax.set_ylabel('Residual (%)')
        ax.set_xlabel(r'Multipole $\ell$')
        ax.set_ylim(-50, 50)

        plt.tight_layout()
        fig_file = os.path.join(out_dir, 'sync_vs_class.png')
        fig.savefig(fig_file, dpi=150)
        print(f"  Plot: {fig_file}")
        plt.close()
    except Exception as e:
        print(f"  (Plot skipped: {e})")

    print(f"\n{'='*70}")
    print(f"SUMMARY: RMS = {rms_full*100:.2f}% vs CLASS")
    print(f"{'='*70}")


if __name__ == '__main__':
    main()
