#!/usr/bin/env python3
"""
test_dk_optimized.py — Quick validation: run solver_accurate with
optimized D(k) parameters and compare to CLASS.
"""
import os
os.environ['MPLCONFIGDIR'] = '/tmp/mpl_config'

import time
import numpy as np
from scipy.interpolate import interp1d
from scipy.signal import find_peaks
from scipy.ndimage import gaussian_filter1d


def main():
    from .background import Background, T_CMB
    from .solver_accurate import AccurateBoltzmannSolver

    # Load CLASS
    class_path = '/Users/akaihuangm1/Desktop/github/gdm_class_public/output/mlx_ref_lcdm_cl.dat'
    data = np.loadtxt(class_path)
    class_ell = data[:, 0].astype(int)
    class_Dl = data[:, 1] * (T_CMB * 1e6) ** 2

    print("=" * 70)
    print("OPTIMIZED D(k) VALIDATION")
    print("=" * 70)

    bg = Background(khronon=False, recombination='peebles')
    bg.solve()

    k_arr = np.geomspace(5e-4, 0.35, 500).astype(np.float32)
    ell_values = np.unique(np.concatenate([
        np.arange(2, 30, 3),
        np.arange(30, 100, 5),
        np.arange(100, 350, 3),
        np.arange(350, 700, 6),
        np.arange(700, 1200, 10),
        np.arange(1200, 2501, 20),
    ])).astype(int)

    # Run solver
    t0 = time.time()
    solver = AccurateBoltzmannSolver(
        bg, k_arr, l_gamma_max=25,
        grid_density='fast', n_snapshots=300,
        constraint_reset_interval=19
    )
    result = solver.solve()
    t_solve = time.time() - t0

    # Compute C_l
    t0 = time.time()
    ell_out, Cl_out, Dl_out = result.compute_cl_los(
        ell_values, include_doppler=True, include_isw=True)
    t_cl = time.time() - t0

    print(f"\nSolve time: {t_solve:.2f}s, C_l time: {t_cl:.1f}s")
    print(f"Total: {t_solve + t_cl:.1f}s")

    # Compare with CLASS
    f_mlx = interp1d(ell_values, Dl_out, kind='cubic', fill_value='extrapolate')
    ell_min = max(int(ell_values[0]), int(class_ell[0]))
    ell_max_c = min(int(ell_values[-1]), int(class_ell[-1]))
    mask = (class_ell >= ell_min) & (class_ell <= ell_max_c)
    ell_c = class_ell[mask]
    Dl_class = class_Dl[mask]
    Dl_mlx = np.maximum(f_mlx(ell_c), 0.0)

    # RMS
    mask_sig = (Dl_class > 100.0) & (ell_c > 30)
    residual = np.where(mask_sig, (Dl_mlx - Dl_class) / Dl_class * 100.0, 0.0)
    rms_all = np.sqrt(np.mean(residual[mask_sig] ** 2))
    max_all = np.max(np.abs(residual[mask_sig]))

    print(f"\n{'='*70}")
    print(f"RESULTS")
    print(f"{'='*70}")
    print(f"  RMS:  {rms_all:.1f}%")
    print(f"  Max:  {max_all:.1f}%")

    for lo, hi in [(100, 500), (500, 1000), (1000, 2000), (2000, 2500)]:
        m = mask_sig & (ell_c >= lo) & (ell_c < hi)
        if np.any(m):
            print(f"  l=[{lo},{hi}): RMS = {np.sqrt(np.mean(residual[m]**2)):.1f}%")

    # Peak analysis
    Dl_cs = gaussian_filter1d(Dl_class, sigma=15)
    Dl_ms = gaussian_filter1d(Dl_mlx, sigma=15)

    mask_peak = ell_c > 100
    offset = np.argmax(mask_peak)
    peaks_c, _ = find_peaks(Dl_cs[mask_peak], distance=150, prominence=50)
    peaks_m, _ = find_peaks(Dl_ms[mask_peak], distance=150, prominence=5)
    peaks_c = peaks_c + offset
    peaks_m = peaks_m + offset

    print(f"\n  PEAK POSITIONS")
    print(f"  {'#':<4} {'l_CLASS':>8} {'l_mlx':>8} {'Err':>10}")
    print(f"  {'-'*34}")
    n_pk = min(len(peaks_c), len(peaks_m), 7)
    for i in range(n_pk):
        l_c = ell_c[peaks_c[i]]
        l_m = ell_c[peaks_m[i]]
        err = (l_m - l_c) / l_c * 100
        print(f"    {i+1:<2} {l_c:>8} {l_m:>8} {err:>+9.1f}%")

    print(f"\n  PEAK AMPLITUDES")
    print(f"  {'#':<4} {'l':>6} {'CLASS':>10} {'mlx':>10} {'Ratio':>8}")
    print(f"  {'-'*44}")
    ratios = []
    for i in range(min(len(peaks_c), 7)):
        p = peaks_c[i]
        d_c = Dl_class[p]
        d_m = Dl_mlx[p]
        r = d_m / d_c if d_c > 0 else 0
        ratios.append(r)
        print(f"    {i+1:<2} {ell_c[p]:>6} {d_c:>10.0f} {d_m:>10.0f} {r:>8.3f}")
    if ratios:
        print(f"    Mean ratio: {np.mean(ratios):.3f}")

    # Plot
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(2, 1, figsize=(14, 10),
                                  gridspec_kw={'height_ratios': [3, 1]}, sharex=True)
        fig.subplots_adjust(hspace=0.06)

        ax = axes[0]
        ax.plot(ell_c, Dl_class, 'k-', lw=1.2, alpha=0.7, label='CLASS')
        ax.plot(ell_values, Dl_out, 'r-', lw=1.5, alpha=0.9,
                label=f'solver_accurate (RMS={rms_all:.1f}%)')
        ax.set_ylabel(r'$\mathcal{D}_\ell^{TT}$ [$\mu$K$^2$]', fontsize=14)
        ax.set_title(f'Optimized solver_accurate vs CLASS (RMS={rms_all:.1f}%)', fontsize=13)
        ax.legend(fontsize=10)
        ax.set_xlim(30, 2500)
        ax.grid(alpha=0.15)

        ax2 = axes[1]
        residual_smooth = gaussian_filter1d(residual, sigma=10)
        ax2.plot(ell_c, residual, 'b-', lw=0.5, alpha=0.3)
        ax2.plot(ell_c, residual_smooth, 'b-', lw=1.5, alpha=0.8)
        ax2.axhline(0, color='k', ls='--', lw=0.8)
        ax2.fill_between(ell_c, -10, 10, color='green', alpha=0.08)
        ax2.set_ylabel('Residual [%]', fontsize=11)
        ax2.set_xlabel(r'$\ell$', fontsize=14)
        ax2.set_ylim(-60, 30)
        ax2.grid(alpha=0.15)
        ax2.text(0.02, 0.90, f'RMS = {rms_all:.1f}%',
                 transform=ax2.transAxes, fontsize=10, va='top',
                 bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

        out_dir = os.path.dirname(os.path.abspath(__file__))
        plot_path = os.path.join(out_dir, 'dk_optimized_validation.png')
        plt.savefig(plot_path, dpi=150, bbox_inches='tight')
        print(f"\n  Plot saved: {plot_path}")
        plt.close()
    except Exception as e:
        print(f"  Plot error: {e}")


if __name__ == '__main__':
    main()
