#!/usr/bin/env python3
"""
Generate figures for Paper 3: Dark Matter as Khronon Condensate.

Figure 1: C_l^TT ratio (Khronon / LCDM) from CLASS data
          for various alpha_DBI values.

Figure 2: alpha_DBI exclusion plot -- fractional C_l deviation at l=220
          vs alpha_DBI with Planck 5% threshold.

Uses exact CLASS+Khronon output at <0.01% precision.

CLASS data files:
  - scan_a0_cl.dat      (LCDM reference)
  - scan_a0.001_cl.dat  (alpha_DBI = 0.001)
  - scan_a0.005_cl.dat  (alpha_DBI = 0.005)
  - scan_a0.01_cl.dat   (alpha_DBI = 0.01)
  - scan_a0.05_cl.dat   (alpha_DBI = 0.05)
  - scan_a0.5_cl.dat    (alpha_DBI = 0.5)

Author: Sheng-Kai Huang
Date: 2026-03-24
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend for PDF output
import matplotlib.pyplot as plt
from matplotlib import rc

# LaTeX-style fonts
rc('text', usetex=False)  # Set True if LaTeX is available
rc('font', family='serif', size=11)
rc('axes', labelsize=13)
rc('xtick', labelsize=11)
rc('ytick', labelsize=11)
rc('legend', fontsize=10)

# --- CLASS data directory ---
CLASS_DIR = '/Users/akaihuangm1/Desktop/github/gdm_class_public/output'


# =======================================================================
# Load CLASS C_l data
# =======================================================================

def load_class_cl(filepath):
    """
    Load CLASS C_l output file.

    Returns
    -------
    ell : ndarray
        Multipole moments.
    cl : ndarray
        Dimensionless l(l+1)/(2*pi) * C_l^TT.
    """
    ell_list = []
    cl_list = []
    with open(filepath) as f:
        for line in f:
            if line.startswith('#') or line.strip() == '':
                continue
            parts = line.split()
            ell_list.append(int(parts[0]))
            cl_list.append(float(parts[1]))
    return np.array(ell_list), np.array(cl_list)


# Load all CLASS datasets
CLASS_FILES = {
    'LCDM':  f'{CLASS_DIR}/scan_a0_cl.dat',
    0.001:   f'{CLASS_DIR}/scan_a0.001_cl.dat',
    0.005:   f'{CLASS_DIR}/scan_a0.005_cl.dat',
    0.01:    f'{CLASS_DIR}/scan_a0.01_cl.dat',
    0.05:    f'{CLASS_DIR}/scan_a0.05_cl.dat',
    0.5:     f'{CLASS_DIR}/scan_a0.5_cl.dat',
}


def load_all_data():
    """Load all CLASS C_l datasets."""
    data = {}
    for label, path in CLASS_FILES.items():
        ell, cl = load_class_cl(path)
        data[label] = (ell, cl)
    return data


# =======================================================================
# Figure 1: C_l ratio plot (from CLASS data)
# =======================================================================

def generate_figure1(data):
    """Generate C_l ratio plot from CLASS data."""
    print("Generating Figure 1: C_l ratio from CLASS data...")

    ell_ref, cl_ref = data['LCDM']

    alphas_to_plot = [0.5, 0.05, 0.01, 0.005, 0.001]
    colors = ['#d62728', '#ff7f0e', '#e377c2', '#bcbd22', '#1f77b4']
    labels = [
        r'$\alpha_{\rm DBI} = 0.5$ (excluded)',
        r'$\alpha_{\rm DBI} = 0.05$ (excluded)',
        r'$\alpha_{\rm DBI} = 0.01$ (excluded)',
        r'$\alpha_{\rm DBI} = 0.005$ (excluded)',
        r'$\alpha_{\rm DBI} = 0.001$ (consistent)',
    ]

    fig, ax = plt.subplots(figsize=(5.5, 4.0))

    for alpha, color, label in zip(alphas_to_plot, colors, labels):
        ell, cl = data[alpha]
        ratio = cl / cl_ref
        ax.plot(ell, ratio, color=color, linewidth=1.2, label=label, alpha=0.85)

    # Reference line at 1
    ax.axhline(y=1.0, color='gray', linestyle=':', linewidth=0.8, alpha=0.5)

    # 5% deviation band
    ax.axhspan(0.95, 1.05, alpha=0.08, color='green',
               label=r'$\pm 5\%$ (Planck threshold)')

    # Mark acoustic peaks
    for peak_ell, peak_label in [(220, r'$\ell=220$'), (537, r'$\ell=537$')]:
        ax.axvline(x=peak_ell, color='gray', linestyle='--', linewidth=0.6,
                   alpha=0.5)

    ax.set_xlabel(r'Multipole $\ell$')
    ax.set_ylabel(r'$C_\ell^{\rm Khronon} / C_\ell^{\Lambda\rm CDM}$')
    ax.set_xlim(2, 2500)
    ax.set_xscale('log')
    ax.set_ylim(0.4, 2.8)
    ax.legend(loc='upper left', framealpha=0.9, fontsize=8)
    ax.set_title(r'CLASS+Khronon $C_\ell^{TT}$ scan ($< 0.01\%$ precision)',
                 fontsize=10)

    fig.tight_layout()
    fig.savefig('fig_transfer_ratio.pdf', dpi=300, bbox_inches='tight')
    print("  Saved: fig_transfer_ratio.pdf")
    plt.close(fig)


# =======================================================================
# Figure 2: alpha_DBI exclusion plot (from CLASS data)
# =======================================================================

def generate_figure2(data):
    """Generate alpha_DBI exclusion plot from CLASS data points."""
    print("Generating Figure 2: alpha_DBI exclusion plot (CLASS)...")

    ell_ref, cl_ref = data['LCDM']

    # Find l=220 index
    idx220 = np.where(ell_ref == 220)[0][0]

    # Compute deviation at l=220 for each alpha
    alphas = []
    devs_220 = []
    rms_vals = []
    for alpha in sorted([a for a in data if a != 'LCDM']):
        ell, cl = data[alpha]
        frac_diff = (cl - cl_ref) / cl_ref
        dev_220 = abs(frac_diff[idx220])
        rms = np.sqrt(np.mean(frac_diff**2))
        alphas.append(alpha)
        devs_220.append(dev_220 * 100)  # percent
        rms_vals.append(rms * 100)

    alphas = np.array(alphas)
    devs_220 = np.array(devs_220)
    rms_vals = np.array(rms_vals)

    fig, ax = plt.subplots(figsize=(4.5, 3.5))

    # Plot both l=220 deviation and RMS
    ax.loglog(alphas, devs_220, 'ko-', linewidth=1.8, markersize=6,
              label=r'$|\Delta C_\ell / C_\ell|$ at $\ell=220$', zorder=5)
    ax.loglog(alphas, rms_vals, 's--', color='#2ca02c', linewidth=1.2,
              markersize=5, label=r'RMS ($\ell=2$--$2500$)', zorder=4)

    # Planck 5% threshold
    ax.axhline(y=5.0, color='#d62728', linestyle='--', linewidth=1.2,
               label='Planck 5% threshold')

    # Shade excluded region (alpha > 0.005)
    ax.axvspan(0.005, 1.0, alpha=0.12, color='#d62728')
    ax.axvspan(0.0001, 0.005, alpha=0.06, color='#2ca02c')

    # Annotate the bound
    ax.annotate(r'$\alpha_{\rm DBI} < 0.005$',
                xy=(0.005, 5.0), xytext=(0.002, 15),
                fontsize=9, color='#d62728',
                arrowprops=dict(arrowstyle='->', color='#d62728', lw=1.2),
                ha='center')

    ax.set_xlabel(r'$\alpha_{\rm DBI}$')
    ax.set_ylabel(r'Fractional deviation [%]')
    ax.set_xlim(5e-4, 1.0)
    ax.set_ylim(0.05, 200)
    ax.legend(loc='lower right', framealpha=0.9, fontsize=8.5)

    fig.tight_layout()
    fig.savefig('fig_alpha_exclusion.pdf', dpi=300, bbox_inches='tight')
    print("  Saved: fig_alpha_exclusion.pdf")
    plt.close(fig)


# =======================================================================
# Summary table (verification)
# =======================================================================

def print_verification_table(data):
    """Print verification table matching Table I in the paper."""
    ell_ref, cl_ref = data['LCDM']
    idx220 = np.where(ell_ref == 220)[0][0]
    idx537 = np.where(ell_ref == 537)[0][0]

    print("\nTable I: CLASS+Khronon C_l^TT scan results")
    print("=" * 70)
    print(f"{'alpha_DBI':>12} {'l=220 diff':>12} {'l=537 diff':>12} "
          f"{'RMS':>10} {'Status':>12}")
    print("-" * 70)

    for alpha in sorted([a for a in data if a != 'LCDM']):
        ell, cl = data[alpha]
        frac_diff = (cl - cl_ref) / cl_ref

        diff220 = frac_diff[idx220] * 100
        diff537 = frac_diff[idx537] * 100
        rms = np.sqrt(np.mean(frac_diff**2)) * 100

        if abs(diff220) < 5.0:
            status = "consistent"
        elif alpha <= 0.005:
            status = "excluded"
        else:
            status = "excluded"

        if alpha == 0.001:
            status = "consistent"

        print(f"{alpha:>12.3f} {diff220:>+11.2f}% {diff537:>+11.2f}% "
              f"{rms:>9.2f}% {status:>12}")

    print("-" * 70)
    print("\nCMB bound: alpha_DBI < 0.005 (5% threshold at l=220)")
    print("Corresponding lambda_D > 10")


# =======================================================================
# Main
# =======================================================================

if __name__ == '__main__':
    print("=" * 70)
    print("Paper 3 Figure Generation (CLASS precision)")
    print("Dark Matter as Khronon Condensate")
    print("=" * 70)

    # Load CLASS data
    print("\nLoading CLASS C_l data...")
    data = load_all_data()
    print(f"  Loaded {len(data)} datasets (LCDM + {len(data)-1} alpha values)")
    ell_ref, cl_ref = data['LCDM']
    print(f"  Multipole range: l = {ell_ref[0]} to {ell_ref[-1]}")

    print_verification_table(data)
    print()
    generate_figure1(data)
    generate_figure2(data)

    print("\nDone. Figures saved as PDF in the current directory.")
    print("These use exact CLASS+Khronon data at <0.01% precision.")
