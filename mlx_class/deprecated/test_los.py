# DEPRECATED: Use tests/test_basic.py instead.
# This file is kept for backward compatibility only.
import warnings as _warnings
_warnings.warn(
    "mlx_class.test_los is deprecated. Use mlx_class.tests.test_basic instead.",
    DeprecationWarning,
    stacklevel=2,
)

#!/usr/bin/env python3
"""
test_los.py -- Test the line-of-sight C_l integration.

Compares:
  1. Instantaneous recombination (old method: spectra.py)
  2. LOS integration with visibility function + Hu-Sugiyama phase shift

Reports peak positions and validates against Planck observed peaks.

Usage:
    python -m mlx_class.test_los [--N_k 500] [--ell_max 2500] [--no_plot]
"""
import os
os.environ['MPLCONFIGDIR'] = '/tmp/mpl_config'

import argparse
import time
import numpy as np
import mlx.core as mx

from .background import Background
from .perturbations import AnalyticTransfer
from .spectra import compute_cl
from .spectra_los import compute_cl_los, find_peaks_in_spectrum


# Planck 2018 observed peak positions
PLANCK_PEAKS = {
    1: {'l': 220.0, 'Dl': 5720},
    2: {'l': 537.5, 'Dl': 2582},
    3: {'l': 810.8, 'Dl': 2509},
    4: {'l': 1120,  'Dl': 1130},
    5: {'l': 1445,  'Dl': 802},
}


def main():
    parser = argparse.ArgumentParser(description='Test LOS C_l integration')
    parser.add_argument('--N_k', type=int, default=500)
    parser.add_argument('--ell_max', type=int, default=2500)
    parser.add_argument('--no_plot', action='store_true')
    args = parser.parse_args()

    print("=" * 70)
    print("Line-of-Sight C_l Integration Test")
    print(f"GPU: {mx.default_device()}")
    print("=" * 70)

    # Background
    print("\n--- Background ---")
    bg = Background()
    bg.solve()

    # Grids
    k_arr = np.geomspace(5e-5, 0.35, args.N_k).astype(np.float32)
    ell_values = np.unique(np.concatenate([
        np.arange(2, 30, 1), np.arange(30, 100, 2),
        np.arange(100, 500, 4), np.arange(500, 1500, 8),
        np.arange(1500, min(2501, args.ell_max + 1), 12),
    ])).astype(int)
    ell_values = ell_values[ell_values <= args.ell_max]

    # ================================================================
    # Method 1: Instantaneous recombination (existing code)
    # ================================================================
    print("\n" + "=" * 70)
    print("Method 1: Instantaneous recombination (analytic, no phase shift)")
    print("=" * 70)
    t0 = time.time()

    transfer = AnalyticTransfer(bg)
    source_SW = transfer.compute_source(k_arr)
    source_Dop = transfer.compute_doppler(k_arr)
    ell_inst, Cl_inst, Dl_inst = compute_cl(
        source_SW, k_arr, ell_values, bg.D_A, source_Dop=source_Dop
    )
    t_inst = time.time() - t0

    peaks_inst, _, _, _ = find_peaks_in_spectrum(ell_inst, Dl_inst)
    print(f"Time: {t_inst:.2f}s")
    print("Peaks (instantaneous + Doppler):")
    for i, (l_peak, dl_peak) in enumerate(peaks_inst[:5]):
        print(f"  Peak {i + 1}: l = {l_peak}, D_l = {dl_peak:.0f} muK^2")

    # ================================================================
    # Method 2: LOS integration with Hu-Sugiyama driving phase
    # ================================================================
    print("\n" + "=" * 70)
    print("Method 2: LOS integration (analytic + driving phase + visibility)")
    print("=" * 70)
    t0 = time.time()

    ell_los, Cl_los, Dl_los, info_los = compute_cl_los(
        bg, k_arr, ell_values
    )
    t_los = time.time() - t0

    peaks_los, ell_fine, Dl_fine, Dl_smooth = find_peaks_in_spectrum(
        ell_los, Dl_los
    )
    print(f"\nTime: {t_los:.2f}s")
    print("Peaks (LOS):")
    for i, (l_peak, dl_peak) in enumerate(peaks_los[:5]):
        print(f"  Peak {i + 1}: l = {l_peak}, D_l = {dl_peak:.0f} muK^2")

    # ================================================================
    # Comparison
    # ================================================================
    print("\n" + "=" * 70)
    print("COMPARISON: LOS vs Instantaneous vs Planck")
    print("=" * 70)

    header = (f"{'Peak':>6} {'Planck':>8} {'Instant.':>10} "
              f"{'LOS':>10} {'Inst err':>10} {'LOS err':>10}")
    print(f"\n{header}")
    print("-" * len(header))

    for i in range(1, min(6, len(peaks_los) + 1)):
        l_planck = PLANCK_PEAKS[i]['l']
        l_inst = peaks_inst[i - 1][0] if i - 1 < len(peaks_inst) else None
        l_los = peaks_los[i - 1][0] if i - 1 < len(peaks_los) else None

        inst_err = f"{abs(l_inst - l_planck) / l_planck * 100:.1f}%" if l_inst else "N/A"
        los_err = f"{abs(l_los - l_planck) / l_planck * 100:.1f}%" if l_los else "N/A"

        l_inst_s = str(l_inst) if l_inst else "N/A"
        l_los_s = str(l_los) if l_los else "N/A"

        print(f"  {i:>4}  {l_planck:>8.1f} {l_inst_s:>10} "
              f"{l_los_s:>10} {inst_err:>10} {los_err:>10}")

    # Peak height ratios
    if len(peaks_los) >= 2:
        print(f"\nPeak height ratios (LOS):")
        print(f"  1st/2nd = {peaks_los[0][1] / peaks_los[1][1]:.2f} "
              f"(Planck: ~2.2)")
    if len(peaks_los) >= 3:
        print(f"  1st/3rd = {peaks_los[0][1] / peaks_los[2][1]:.2f} "
              f"(Planck: ~2.3)")

    # SW decomposition
    if 'Dl_sw_only' in info_los:
        Dl_sw = info_los['Dl_sw_only']
        Dl_swd = info_los['Dl_sw_dop']
        mask_low = ell_values < 30
        if np.sum(mask_low) > 0:
            print(f"\nLow-l (l<30) mean D_l:")
            print(f"  SW only:      {np.mean(Dl_sw[mask_low]):.0f} muK^2")
            print(f"  SW+Doppler:   {np.mean(Dl_swd[mask_low]):.0f} muK^2")
            print(f"  Full (+ISW):  {np.mean(Dl_los[mask_low]):.0f} muK^2")

    # Timing
    print(f"\nTiming: instantaneous {t_inst:.2f}s, LOS {t_los:.2f}s")

    # ================================================================
    # Plot
    # ================================================================
    if not args.no_plot:
        _plot(ell_inst, Dl_inst, peaks_inst,
              ell_los, Dl_los, peaks_los, info_los,
              ell_fine, Dl_smooth, args)

    # Save
    out_dir = os.path.dirname(os.path.abspath(__file__))
    data_path = os.path.join(out_dir, 'cl_lcdm_los.dat')
    np.savetxt(data_path,
               np.column_stack([ell_los, Cl_los, Dl_los]),
               header='ell  C_l  D_l[muK^2]', fmt='%6d  %.8e  %.8e')
    print(f"\nData saved: {data_path}")


def _plot(ell_inst, Dl_inst, peaks_inst,
          ell_los, Dl_los, peaks_los, info_los,
          ell_fine, Dl_smooth, args):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 12),
                                        gridspec_kw={'height_ratios': [3, 1]})

        # Top: spectra
        ax1.plot(ell_inst, Dl_inst, 'b--', alpha=0.6, lw=1,
                 label='Instantaneous (no phase shift)')
        ax1.plot(ell_los, Dl_los, 'r-', lw=1.5,
                 label='LOS (Hu-Sugiyama phase + visibility)')
        ax1.plot(ell_fine, Dl_smooth, 'r-', alpha=0.25, lw=3)

        # Planck peak positions
        for i in range(1, 4):
            ax1.axvline(PLANCK_PEAKS[i]['l'], color='green', ls=':', alpha=0.5)

        # Mark LOS peaks
        colors = ['red', 'orange', 'darkgreen', 'purple', 'brown']
        for i, (l_p, dl_p) in enumerate(peaks_los[:5]):
            ax1.plot(l_p, dl_p, 'v', color=colors[i], ms=10, zorder=5)
            ax1.annotate(f'l={l_p}', (l_p, dl_p * 1.08),
                         fontsize=10, ha='center', color=colors[i],
                         fontweight='bold')

        # Mark instantaneous peaks
        for i, (l_p, dl_p) in enumerate(peaks_inst[:5]):
            ax1.plot(l_p, dl_p, '^', color='blue', ms=7, alpha=0.5, zorder=4)

        ax1.set_xlabel('Multipole l', fontsize=13)
        ax1.set_ylabel(r'$D_\ell$ [$\mu$K$^2$]', fontsize=13)
        ax1.set_title('CMB TT: LOS Integration vs Instantaneous', fontsize=14)
        ax1.legend(fontsize=11)
        ax1.set_xlim(2, args.ell_max)
        max_Dl = max(np.max(Dl_inst), np.max(Dl_los))
        if max_Dl > 0:
            ax1.set_ylim(0, max_Dl * 1.4)

        # Bottom: SW decomposition
        if 'Dl_sw_only' in info_los:
            ax2.plot(ell_los, info_los['Dl_sw_only'], 'b-', lw=1, alpha=0.7,
                     label='SW only')
            ax2.plot(ell_los, info_los['Dl_sw_dop'], 'g-', lw=1, alpha=0.7,
                     label='SW + Doppler')
            ax2.plot(ell_los, Dl_los, 'r-', lw=1.5,
                     label='Full (SW+Dop+ISW)')
            ax2.legend(fontsize=10)

        ax2.set_xlabel('Multipole l', fontsize=13)
        ax2.set_ylabel(r'$D_\ell$ [$\mu$K$^2$]', fontsize=13)
        ax2.set_title('Source Decomposition', fontsize=12)
        ax2.set_xlim(2, args.ell_max)

        plt.tight_layout()
        out_dir = os.path.dirname(os.path.abspath(__file__))
        out_path = os.path.join(out_dir, 'test_los_comparison.png')
        plt.savefig(out_path, dpi=150)
        print(f"Plot: {out_path}")
        plt.close()

    except ImportError:
        print("matplotlib not available, skipping plot.")


if __name__ == '__main__':
    main()
