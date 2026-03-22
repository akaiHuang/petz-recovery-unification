"""
test_lmax_scan.py -- Find optimal l_gamma_max for sync gauge solver.

Scans l_gamma_max in [20, 25, 30, 35, 40] using N_k=100 k-modes for speed.
Compares TT power spectrum against CLASS reference and reports:
  - RMS error
  - D_l ratios at key multipoles
  - Peak positions

Usage:
  python -m mlx_class.test_lmax_scan
"""
import numpy as np
import sys
import time

from .solver_sync import run_sync_solver
from .background import T_CMB
from scipy.interpolate import interp1d
from scipy.signal import find_peaks
from scipy.ndimage import gaussian_filter1d


def load_class_reference():
    """Load CLASS TT reference spectrum."""
    class_path = '/Users/akaihuangm1/Desktop/github/gdm_class_public/output/mlx_ref_lcdm_cl.dat'
    data = np.loadtxt(class_path)
    ell_class = data[:, 0].astype(int)
    # Column 2 is dimensionless D_l = l(l+1)C_l/(2pi). Convert to uK^2.
    Dl_class = data[:, 1] * (T_CMB * 1e6) ** 2
    return ell_class, Dl_class


def find_peaks_in_spectrum(ell, Dl, n_peaks=5):
    """Find peak positions in D_l spectrum using smoothing + scipy find_peaks."""
    Dl_smooth = gaussian_filter1d(Dl, sigma=8)
    mask = ell > 100
    offset = np.argmax(mask)
    pks, props = find_peaks(Dl_smooth[mask], distance=60, prominence=10)
    peak_ells = [int(ell[p + offset]) for p in pks[:n_peaks]]
    return peak_ells


def run_scan():
    """Run l_gamma_max scan and print comparison table."""
    # Load CLASS reference
    ell_class, Dl_class = load_class_reference()
    class_interp = interp1d(ell_class, Dl_class, kind='linear',
                            fill_value='extrapolate')

    # CLASS peak positions (for reference)
    class_peaks = find_peaks_in_spectrum(ell_class, Dl_class)
    print(f"CLASS peak positions: {class_peaks}")

    lg_max_values = [20, 25, 30, 35, 40]
    key_ells = [100, 220, 537, 820, 1200]
    results = {}

    for lg_max in lg_max_values:
        print(f"\n{'=' * 65}")
        print(f"  l_gamma_max = {lg_max}")
        print(f"{'=' * 65}")
        sys.stdout.flush()

        t0 = time.time()
        try:
            result = run_sync_solver(N_k=100, lg_max=lg_max, verbose=False)
        except Exception as e:
            print(f"  FAILED: {e}")
            results[lg_max] = None
            continue
        elapsed = time.time() - t0

        ell = result['ell']
        Dl = result['Dl_TT']

        # Compute RMS error (only where CLASS > 100 uK^2 to avoid noise floor)
        class_at_ell = class_interp(ell)
        mask = class_at_ell > 100
        if np.sum(mask) > 0:
            rms = 100.0 * np.sqrt(np.mean(
                ((Dl[mask] - class_at_ell[mask]) / class_at_ell[mask]) ** 2
            ))
        else:
            rms = np.nan

        # D_l ratios at key multipoles
        ratios = {}
        for l_target in key_ells:
            idx = np.argmin(np.abs(ell - l_target))
            class_val = class_interp(l_target)
            if class_val > 0:
                ratios[l_target] = Dl[idx] / class_val
            else:
                ratios[l_target] = np.nan

        # Find peaks
        peaks = find_peaks_in_spectrum(ell, Dl)

        results[lg_max] = {
            'rms': rms,
            'time': elapsed,
            'ratios': ratios,
            'peaks': peaks,
            'ell': ell,
            'Dl': Dl,
        }

        print(f"  RMS:   {rms:.2f}%")
        print(f"  Time:  {elapsed:.0f}s")
        print(f"  Peaks: {peaks}")
        for l_target, r in ratios.items():
            print(f"  D_l({l_target})/CLASS: {r:.3f}")
        sys.stdout.flush()

    # ================================================================
    # Summary table
    # ================================================================
    print(f"\n{'=' * 90}")
    print(f"  SUMMARY TABLE")
    print(f"{'=' * 90}")

    header = (f"  {'lg_max':>8s} {'RMS%':>8s} {'Time':>7s} "
              f"{'l=220':>8s} {'l=537':>8s} {'l=820':>8s} {'l=1200':>8s} "
              f"{'Peak1':>7s} {'Peak2':>7s} {'Peak3':>7s}")
    print(header)
    print(f"  {'-' * 85}")

    for lg_max in lg_max_values:
        r = results.get(lg_max)
        if r is None:
            print(f"  {lg_max:>8d}   FAILED")
            continue

        p1 = str(r['peaks'][0]) if len(r['peaks']) > 0 else 'N/A'
        p2 = str(r['peaks'][1]) if len(r['peaks']) > 1 else 'N/A'
        p3 = str(r['peaks'][2]) if len(r['peaks']) > 2 else 'N/A'

        print(f"  {lg_max:>8d} {r['rms']:>8.2f} {r['time']:>6.0f}s "
              f"{r['ratios'][220]:>8.3f} {r['ratios'][537]:>8.3f} "
              f"{r['ratios'][820]:>8.3f} {r['ratios'][1200]:>8.3f} "
              f"{p1:>7s} {p2:>7s} {p3:>7s}")

    # CLASS reference row
    cp1 = str(class_peaks[0]) if len(class_peaks) > 0 else 'N/A'
    cp2 = str(class_peaks[1]) if len(class_peaks) > 1 else 'N/A'
    cp3 = str(class_peaks[2]) if len(class_peaks) > 2 else 'N/A'
    print(f"  {'CLASS':>8s} {'0.00':>8s} {'---':>7s} "
          f"{'1.000':>8s} {'1.000':>8s} {'1.000':>8s} {'1.000':>8s} "
          f"{cp1:>7s} {cp2:>7s} {cp3:>7s}")

    # ================================================================
    # Recommendation
    # ================================================================
    print(f"\n{'=' * 90}")
    print(f"  RECOMMENDATION")
    print(f"{'=' * 90}")

    # Find the lg_max with lowest RMS
    valid = {k: v for k, v in results.items() if v is not None}
    if valid:
        best_rms = min(valid.items(), key=lambda x: x[1]['rms'])
        print(f"  Best RMS:  lg_max={best_rms[0]} -> {best_rms[1]['rms']:.2f}%")

        # Find best peak 3 accuracy (closest to CLASS peak 3)
        if len(class_peaks) > 2:
            class_p3 = class_peaks[2]
            best_p3 = min(
                ((k, v) for k, v in valid.items() if len(v['peaks']) > 2),
                key=lambda x: abs(x[1]['peaks'][2] - class_p3),
                default=None
            )
            if best_p3:
                p3_err = best_p3[1]['peaks'][2] - class_p3
                print(f"  Best Peak3: lg_max={best_p3[0]} -> "
                      f"l={best_p3[1]['peaks'][2]} (CLASS={class_p3}, "
                      f"shift={p3_err:+d})")

        # Composite score: weight RMS 70%, peak3 accuracy 30%
        print(f"\n  Composite score (70% RMS + 30% peak3 accuracy):")
        scores = {}
        for k, v in valid.items():
            rms_score = v['rms']
            if len(class_peaks) > 2 and len(v['peaks']) > 2:
                p3_score = abs(v['peaks'][2] - class_peaks[2]) / class_peaks[2] * 100
            else:
                p3_score = 50.0  # penalty
            composite = 0.7 * rms_score + 0.3 * p3_score
            scores[k] = composite
            print(f"    lg_max={k}: RMS={rms_score:.2f}%, "
                  f"Peak3_err={p3_score:.1f}%, composite={composite:.2f}")

        best_composite = min(scores.items(), key=lambda x: x[1])
        print(f"\n  >>> RECOMMENDED: lg_max = {best_composite[0]} "
              f"(composite score = {best_composite[1]:.2f}) <<<")
    else:
        print("  No valid results to compare.")


if __name__ == '__main__':
    run_scan()
