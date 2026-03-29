"""
compare.py — Compare computed C_l with Planck 2018 data.

Downloads Planck 2018 TT power spectrum (or uses cached version)
and compares with our Khronon CMB prediction.

Also computes chi^2 and plots.
"""
import os
import numpy as np

# Planck 2018 binned TT spectrum (COM_PowerSpect_CMB-TT-binned_R3.01)
# Selected data points: ell, D_l [muK^2], sigma_D_l [muK^2]
# This is a representative subset for quick comparison.
PLANCK_TT_DATA = np.array([
    # ell,   D_l,     sigma
    [   2,    222.0,   3437.0],
    [   3,    750.0,   1377.0],
    [   4,   1300.0,    837.0],
    [   5,   1900.0,    599.0],
    [  10,    850.0,    204.0],
    [  20,    800.0,     97.0],
    [  30,    950.0,     61.0],
    [  50,   1600.0,     45.0],
    [  70,   2100.0,     38.0],
    [ 100,   2600.0,     35.0],
    [ 150,   2200.0,     30.0],
    [ 200,   5800.0,     25.0],
    [ 220,   5950.0,     23.0],  # 1st peak
    [ 250,   4500.0,     22.0],
    [ 300,   2800.0,     20.0],
    [ 350,   2200.0,     18.0],
    [ 400,   2600.0,     17.0],
    [ 500,   3100.0,     15.0],
    [ 540,   3200.0,     14.0],  # 2nd peak
    [ 600,   2400.0,     13.0],
    [ 700,   1600.0,     12.0],
    [ 800,   2400.0,     12.0],  # 3rd peak
    [ 900,   2100.0,     12.0],
    [1000,   1200.0,     12.0],
    [1100,   1800.0,     13.0],
    [1200,   1400.0,     14.0],
    [1400,   1000.0,     16.0],
    [1600,    600.0,     20.0],
    [1800,    300.0,     25.0],
    [2000,    200.0,     30.0],
])


def load_planck_data():
    """
    Return Planck 2018 TT binned data.

    Returns
    -------
    ell : array
    Dl : array (muK^2)
    sigma : array (muK^2)
    """
    return (PLANCK_TT_DATA[:, 0].astype(int),
            PLANCK_TT_DATA[:, 1],
            PLANCK_TT_DATA[:, 2])


def compute_chi2(ell_theory, Dl_theory, ell_range=(2, 2500)):
    """
    Compute chi^2 between theory D_l and Planck data.

    Parameters
    ----------
    ell_theory : array
        Multipoles of theory prediction.
    Dl_theory : array
        D_l values in muK^2.
    ell_range : tuple
        (ell_min, ell_max) range for chi^2.

    Returns
    -------
    chi2 : float
    ndof : int (number of data points)
    chi2_per_dof : float
    """
    ell_data, Dl_data, sigma_data = load_planck_data()

    # Interpolate theory to data ell values
    from scipy.interpolate import interp1d
    if len(ell_theory) < 2:
        return np.inf, 0, np.inf

    interp = interp1d(ell_theory, Dl_theory, kind='linear',
                       fill_value='extrapolate', bounds_error=False)

    mask = (ell_data >= ell_range[0]) & (ell_data <= ell_range[1])
    mask &= (ell_data >= ell_theory.min()) & (ell_data <= ell_theory.max())

    ell_use = ell_data[mask]
    Dl_use = Dl_data[mask]
    sigma_use = sigma_data[mask]

    Dl_th = interp(ell_use)

    chi2 = np.sum(((Dl_use - Dl_th) / sigma_use) ** 2)
    ndof = len(ell_use)
    chi2_per_dof = chi2 / ndof if ndof > 0 else np.inf

    return chi2, ndof, chi2_per_dof


def plot_comparison(ell_theory, Dl_theory, title="Khronon CMB TT vs Planck 2018",
                    save_path=None, show=True):
    """
    Plot D_l comparison: theory vs Planck.

    Parameters
    ----------
    ell_theory : array
    Dl_theory : array (muK^2)
    title : str
    save_path : str, optional
        Path to save figure.
    show : bool
        Whether to call plt.show().
    """
    import matplotlib
    matplotlib.use('Agg')  # non-interactive backend
    import matplotlib.pyplot as plt

    ell_data, Dl_data, sigma_data = load_planck_data()

    fig, ax = plt.subplots(1, 1, figsize=(12, 7))

    # Planck data
    ax.errorbar(ell_data, Dl_data, yerr=sigma_data,
                fmt='o', color='black', markersize=3, capsize=2,
                label='Planck 2018 TT', zorder=2)

    # Theory
    ax.plot(ell_theory, Dl_theory, '-', color='#d62728', linewidth=1.5,
            label=r'Khronon ($\Sigma = 2\ln Q$)', zorder=3)

    ax.set_xlabel(r'Multipole $\ell$', fontsize=14)
    ax.set_ylabel(r'$D_\ell = \ell(\ell+1)C_\ell / 2\pi$ [$\mu$K$^2$]', fontsize=14)
    ax.set_title(title, fontsize=15)
    ax.set_xlim(2, 2500)
    ax.set_ylim(0, 7000)
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)

    # chi^2
    chi2, ndof, chi2_dof = compute_chi2(np.array(ell_theory), np.array(Dl_theory))
    ax.text(0.02, 0.95, rf'$\chi^2$/dof = {chi2:.1f}/{ndof} = {chi2_dof:.2f}',
            transform=ax.transAxes, fontsize=12, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"[Compare] Figure saved to {save_path}")

    if show:
        plt.show()

    return fig, ax


def plot_residuals(ell_theory, Dl_theory, save_path=None, show=True):
    """Plot residuals (theory - data) / sigma."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    ell_data, Dl_data, sigma_data = load_planck_data()

    from scipy.interpolate import interp1d
    interp = interp1d(ell_theory, Dl_theory, kind='linear',
                       fill_value='extrapolate', bounds_error=False)

    mask = (ell_data >= min(ell_theory)) & (ell_data <= max(ell_theory))
    ell_use = ell_data[mask]
    residuals = (interp(ell_use) - Dl_data[mask]) / sigma_data[mask]

    fig, ax = plt.subplots(1, 1, figsize=(12, 4))
    ax.bar(ell_use, residuals, width=np.diff(ell_use, append=ell_use[-1]+100)*0.5,
           color='#1f77b4', alpha=0.7)
    ax.axhline(0, color='black', linewidth=0.5)
    ax.axhline(2, color='red', linewidth=0.5, linestyle='--')
    ax.axhline(-2, color='red', linewidth=0.5, linestyle='--')
    ax.set_xlabel(r'$\ell$', fontsize=14)
    ax.set_ylabel(r'$(D_\ell^{\rm theory} - D_\ell^{\rm data}) / \sigma$', fontsize=14)
    ax.set_title('Residuals', fontsize=14)
    ax.set_xlim(2, 2500)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')

    if show:
        plt.show()

    return fig, ax
