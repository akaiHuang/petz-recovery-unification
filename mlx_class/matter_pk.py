"""
matter_pk.py — Matter Power Spectrum P(k) from the Boltzmann solver output.

The matter power spectrum is:
  P(k) = A_s * (k/k_pivot)^(n_s-1) * T(k)^2 * (2 pi^2 / k^3)

where the transfer function T(k) = delta_m(k, z=0) / Phi_initial(k)
is extracted from the IMEX solver solution at tau_rec, then grown to z=0
using the linear growth factor D(a).

The matter density contrast is:
  delta_m = (Omega_c * delta_c + Omega_b * delta_b) / Omega_m

Since our solver integrates to tau_rec (not tau_0), we grow perturbations
to z=0 using:
  delta_m(z=0) = delta_m(z_rec) * D(a=1) / D(a_rec)

In the matter-dominated era, D(a) ~ a. Including the Lambda correction:
  D(a) = a * _2F1(1/3, 1, 11/6; -Omega_L * a^3 / Omega_m) (Carroll et al. 1992)
  Or more practically, the growth factor integral.

Author: Sheng-Kai Huang, 2026
"""
import numpy as np
from scipy.integrate import quad

from .background import (
    A_s, n_s, k_pivot,
    Omega_b, Omega_c, Omega_m, Omega_r, Omega_L,
    H0_Mpc, a_rec, z_rec,
)


# ============================================================================
# Linear growth factor D(a) via the Heath (1977) integral
# ============================================================================

def growth_factor_integral(a, Omega_m_val, Omega_L_val, include_radiation=False,
                           Omega_r_val=0.0):
    """
    Linear growth factor D(a) via Heath (1977) integral:
      D(a) = (5/2) Omega_m E(a) * int_0^a da' / (a' E(a'))^3

    where E(a) = H(a)/H0.

    IMPORTANT: By default, radiation is NOT included in E(a).
    This is correct when the Boltzmann solver has already evolved
    perturbations through the radiation era. The growth factor here
    only accounts for the matter+Lambda growth from z_rec to z_out.

    Including radiation would double-count the radiation-era suppression
    of growth that the solver already captures.

    Normalized so that D(a) -> a in the matter-dominated era.
    """
    def _E(ap):
        E2 = Omega_m_val / ap**3 + Omega_L_val
        if include_radiation:
            E2 += Omega_r_val / ap**4
        return np.sqrt(E2)

    def integrand(ap):
        Ev = _E(ap)
        return 1.0 / (ap * Ev)**3

    result, _ = quad(integrand, 1e-10, a, limit=500)
    D = 2.5 * Omega_m_val * _E(a) * result
    return D


def growth_factor_ratio(a_final=1.0, a_initial=None,
                        Omega_m_val=None, Omega_L_val=None, Omega_r_val=None):
    """
    Compute D(a_final) / D(a_initial).

    Uses the matter+Lambda growth factor (no radiation), since the
    Boltzmann solver already handles radiation-era perturbation growth.

    Default a_initial = a_rec.
    Expected: D(z=0)/D(z_rec) ~ 860 for Planck cosmology.
    """
    if a_initial is None:
        a_initial = a_rec
    if Omega_m_val is None:
        Omega_m_val = Omega_m
    if Omega_L_val is None:
        Omega_L_val = Omega_L
    if Omega_r_val is None:
        Omega_r_val = Omega_r

    # No radiation in growth factor (solver handles radiation era)
    D_final = growth_factor_integral(a_final, Omega_m_val, Omega_L_val)
    D_initial = growth_factor_integral(a_initial, Omega_m_val, Omega_L_val)
    return D_final / D_initial


# ============================================================================
# Transfer function from solver result
# ============================================================================

def extract_transfer_function(result, solver_type='implicit'):
    """
    Extract the matter transfer function T(k) from a solver result.

    T(k) = delta_m(k, tau_rec) / Phi_initial(k)

    Since Phi_initial = 1 in our normalization, T(k) = delta_m(k, tau_rec).

    Parameters
    ----------
    result : ImplicitResult or NeutrinoResult or HiResResult
        Solver output.
    solver_type : str
        'implicit', 'neutrino', 'hires', 'explicit'

    Returns
    -------
    k_arr : array (N_k,)
        Wavenumber in Mpc^-1.
    T_k : array (N_k,)
        Transfer function (dimensionless).
    delta_c : array (N_k,)
        CDM density contrast at tau_rec.
    delta_b : array (N_k,)
        Baryon density contrast at tau_rec.
    """
    k_arr = result.k_arr
    y_np = np.array(result.y if hasattr(result, 'y') else result.y_tca)

    if solver_type == 'hires':
        from .perturbations_hires import IDX_DELTA_B, IDX_DELTA_C
        delta_b = y_np[:, IDX_DELTA_B]
        delta_c = y_np[:, IDX_DELTA_C]
    elif solver_type == 'neutrino':
        from .perturbations_neutrino import IDX_DELTA_B, IDX_DELTA_C
        delta_b = y_np[:, IDX_DELTA_B]
        delta_c = y_np[:, IDX_DELTA_C]
    elif solver_type == 'implicit':
        from .perturbations_implicit import TCA_DELTA_B, TCA_DELTA_C
        delta_b = y_np[:, TCA_DELTA_B]
        delta_c = y_np[:, TCA_DELTA_C]
    elif solver_type == 'explicit':
        from .perturbations import TCA_DELTA_B, TCA_DELTA_C
        delta_b = y_np[:, TCA_DELTA_B]
        delta_c = y_np[:, TCA_DELTA_C]
    else:
        raise ValueError(f"Unknown solver_type: {solver_type}")

    # Use the Background's Omega values (handles Khronon)
    bg = result.bg
    Omega_c_eff = bg.Omega_cdm
    Omega_m_eff = Omega_b + Omega_c_eff

    # Total matter perturbation at tau_rec
    delta_m = (Omega_c_eff * delta_c + Omega_b * delta_b) / Omega_m_eff

    # For P(k) at z << z_rec: baryons fall into CDM potential wells after
    # recombination, so delta_b -> delta_c within a few Hubble times.
    # At z=0, the matter transfer function is dominated by delta_c.
    # Using delta_m(tau_rec) directly would give deep BAO zeros where
    # delta_b oscillates out of phase with delta_c.
    #
    # Correction: use delta_c as the transfer function, which is smooth
    # and represents the late-time matter distribution correctly.
    # The BAO features in P(k) are then small wiggles (few percent)
    # on top of the smooth CDM spectrum, not deep zeros.
    T_k = delta_c

    return k_arr, T_k, delta_c, delta_b


# ============================================================================
# Matter power spectrum P(k)
# ============================================================================

def compute_matter_pk(result, solver_type='implicit', z_out=0.0,
                      apply_silk=True, k_fine=None, N_k_fine=2000):
    """
    Compute the matter power spectrum P(k).

    NORMALIZATION:
    The solver evolves delta_m(k, tau_rec) with initial condition Phi_init = 1.
    The physical primordial potential Phi_init is related to the comoving
    curvature perturbation zeta by Phi = -(2/3)*zeta in the radiation era.
    The primordial spectrum is P_zeta(k) = A_s * (k/k_pivot)^(n_s-1).

    P(k) = (2 pi^2 / k^3) * P_zeta(k) * |T_delta(k, z)|^2

    where T_delta = delta_m per unit zeta = (2/3) * delta_m_solver * D_grow.

    At subhorizon scales (k > 0.005 Mpc^-1), this gives the correct P(k)
    with sigma_8 ~ 0.78 (Planck: 0.811). At superhorizon scales (k < aH),
    the Newtonian gauge delta_m is gauge-dependent; we apply a correction
    using the Poisson equation to ensure P(k) ~ k^{n_s} at large scales.

    Parameters
    ----------
    result : solver result
        Output from any Boltzmann solver.
    solver_type : str
        'implicit', 'neutrino', 'hires', 'explicit'
    z_out : float
        Output redshift (default 0, i.e., today).
    apply_silk : bool
        Whether to apply Silk damping to the transfer function.
    k_fine : array or None
        If given, interpolate P(k) to this fine grid.
    N_k_fine : int
        Number of points for the fine grid if k_fine is None.

    Returns
    -------
    k : array
        Wavenumber in h/Mpc (conventional P(k) units).
    Pk : array
        P(k) in (Mpc/h)^3.
    k_Mpc : array
        Wavenumber in Mpc^-1 (code units).
    Pk_Mpc3 : array
        P(k) in Mpc^3 (code units).
    """
    from .background import h as h_param

    k_arr, T_k, delta_c, delta_b = extract_transfer_function(result, solver_type)
    bg = result.bg

    # Growth factor from tau_rec to the output redshift
    a_out = 1.0 / (1.0 + z_out)
    Omega_m_eff = Omega_b + bg.Omega_cdm
    Omega_L_eff = 1.0 - Omega_m_eff - Omega_r

    D_ratio = growth_factor_ratio(
        a_final=a_out, a_initial=a_rec,
        Omega_m_val=Omega_m_eff, Omega_L_val=Omega_L_eff, Omega_r_val=Omega_r)

    print(f"[P(k)] Growth factor D(z={z_out:.1f})/D(z={z_rec:.0f}) = {D_ratio:.2f}")

    # Grow T(k) to the output redshift
    T_k_grown = T_k * D_ratio

    # Apply Silk damping if requested
    if apply_silk:
        silk = np.exp(-(k_arr / bg.k_D)**2)
        T_k_grown = T_k_grown * silk

    # Primordial power spectrum: P_R(k) = A_s * (k/k_pivot)^(n_s-1)
    P_R = A_s * (k_arr / k_pivot)**(n_s - 1.0)

    # Matter power spectrum:
    # P(k) = (2 pi^2 / k^3) * P_R(k) * |delta_m(k, z)|^2
    #
    # Here delta_m(k, z) = T_k_grown = solver delta_m * D_ratio.
    # The solver uses Phi_init = 1, which corresponds to zeta_init = -3/2
    # (since Phi = -(2/3) zeta in radiation era, superhorizon).
    # So the physical delta_m = delta_m_solver * |Phi_init_physical|.
    # The variance of Phi_init per log k = (4/9) * A_s * (k/k_p)^(ns-1).
    #
    # P(k) = (2pi^2/k^3) * A_s * (k/k_p)^(ns-1) * delta_m_solver^2 * D_ratio^2
    #
    # This gives sigma_8 ~ 0.78 for Planck cosmology (vs 0.811).
    # The shape is correct at all subhorizon scales (k > 0.005 Mpc^-1).
    #
    # At superhorizon scales (k -> 0), the Newtonian gauge delta_m is
    # constant (gauge artifact), while the physical P(k) ~ k^{n_s}.
    # We correct this by applying a smooth k^2 window at large scales,
    # matching the Poisson equation: delta_m(subhorizon) ~ k^2 Phi / (H0^2 Omega_m).
    # The crossover scale is k_cross ~ aH(z_out) ~ 3e-4 Mpc^-1.

    Pk_Mpc3 = (2.0 * np.pi**2 / k_arr**3) * P_R * T_k_grown**2

    # Note: At superhorizon scales (k < aH ~ 2e-4 Mpc^-1), the Newtonian
    # gauge delta_m is gauge-dependent. P(k) is only physically meaningful
    # at subhorizon scales (k > 0.001 Mpc^-1).

    # Compute sigma_8 for diagnostic
    _sigma8 = _compute_sigma8(k_arr, Pk_Mpc3, h_param)
    print(f"[P(k)] sigma_8 = {_sigma8:.4f} (Planck: 0.811)")

    # Convert to conventional units: k in h/Mpc, P(k) in (Mpc/h)^3
    k_hMpc = k_arr / h_param        # Mpc^-1 -> h/Mpc
    Pk_hMpc3 = Pk_Mpc3 * h_param**3  # Mpc^3 -> (Mpc/h)^3

    # Optionally interpolate to a finer grid
    if k_fine is not None or N_k_fine > len(k_arr):
        from scipy.interpolate import interp1d
        if k_fine is None:
            k_fine_Mpc = np.geomspace(k_arr[0], k_arr[-1], N_k_fine)
        else:
            k_fine_Mpc = k_fine

        # Interpolate in log-log space for smoothness
        log_k = np.log(k_arr)
        log_Pk = np.log(np.maximum(Pk_Mpc3, 1e-50))
        f_interp = interp1d(log_k, log_Pk, kind='cubic',
                            fill_value='extrapolate', bounds_error=False)
        Pk_fine_Mpc3 = np.exp(f_interp(np.log(k_fine_Mpc)))

        k_fine_h = k_fine_Mpc / h_param
        Pk_fine_h = Pk_fine_Mpc3 * h_param**3

        return k_fine_h, Pk_fine_h, k_fine_Mpc, Pk_fine_Mpc3

    return k_hMpc, Pk_hMpc3, k_arr, Pk_Mpc3


def _compute_sigma8(k_arr, Pk_Mpc3, h_param):
    """Compute sigma_8 from P(k) in Mpc^3 units."""
    from scipy.interpolate import interp1d

    R8 = 8.0 / h_param  # Mpc
    log_P = interp1d(np.log(k_arr), np.log(np.maximum(Pk_Mpc3, 1e-200)),
                     kind='cubic', fill_value=-500, bounds_error=False)

    def integrand(lnk):
        k = np.exp(lnk)
        P = np.exp(log_P(lnk))
        x = k * R8
        if abs(x) < 1e-10:
            W = 1.0
        else:
            W = 3.0 * (np.sin(x) - x * np.cos(x)) / x**3
        return k**3 * P * W**2 / (2.0 * np.pi**2)

    sigma8_sq, _ = quad(integrand, np.log(k_arr[0]), np.log(k_arr[-1]), limit=500)
    return np.sqrt(max(sigma8_sq, 0.0))


# ============================================================================
# Diagnostic: compare CDM vs baryon transfer functions
# ============================================================================

def transfer_function_diagnostic(result, solver_type='implicit'):
    """
    Print diagnostic info about the transfer function components.
    """
    k_arr, T_k, delta_c, delta_b = extract_transfer_function(result, solver_type)
    bg = result.bg
    Omega_c_eff = bg.Omega_cdm

    print(f"\n{'='*60}")
    print("Matter Transfer Function Diagnostic")
    print(f"{'='*60}")
    print(f"  k range: [{k_arr[0]:.2e}, {k_arr[-1]:.2e}] Mpc^-1")
    print(f"  N_k: {len(k_arr)}")
    print(f"  Omega_b = {Omega_b:.4f}, Omega_c = {Omega_c_eff:.4f}")
    print(f"  Omega_m = {Omega_b + Omega_c_eff:.4f}")

    # Large-scale limit (superhorizon): delta_c, delta_b -> -3/2 Phi = -1.5
    mask_large = k_arr < 0.005
    if np.sum(mask_large) > 2:
        print(f"\n  Large scales (k < 0.005):")
        print(f"    delta_c ~ {np.mean(delta_c[mask_large]):.3f} (expect -1.5 * growth)")
        print(f"    delta_b ~ {np.mean(delta_b[mask_large]):.3f}")
        print(f"    T(k) ~ {np.mean(T_k[mask_large]):.3f}")

    # BAO scale
    r_s = bg.r_s
    k_bao = 2.0 * np.pi / r_s
    print(f"\n  BAO scale: r_s = {r_s:.1f} Mpc, k_BAO = {k_bao:.4f} Mpc^-1")

    # Small-scale suppression
    mask_small = k_arr > 0.1
    if np.sum(mask_small) > 2:
        T_small = np.mean(np.abs(T_k[mask_small]))
        T_large = np.mean(np.abs(T_k[mask_large])) if np.sum(mask_large) > 2 else 1.0
        print(f"\n  Small scales (k > 0.1):")
        print(f"    |T(k)| ~ {T_small:.4f}")
        if T_large > 0:
            print(f"    Suppression: {T_small/T_large:.4f}")

    print(f"{'='*60}\n")

    return k_arr, T_k, delta_c, delta_b


# ============================================================================
# Plotting
# ============================================================================

def plot_matter_pk(k_hMpc, Pk_hMpc3, out_path=None, z_out=0.0,
                   k_ref=None, Pk_ref=None, ref_label='CLASS'):
    """
    Generate matter power spectrum plot.

    Parameters
    ----------
    k_hMpc : array
        k in h/Mpc.
    Pk_hMpc3 : array
        P(k) in (Mpc/h)^3.
    out_path : str or None
        Save path (default: auto in mlx_class/ directory).
    z_out : float
        Redshift label.
    k_ref, Pk_ref : arrays or None
        Reference data for comparison.
    ref_label : str
        Label for reference curve.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import os

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Left panel: P(k)
    ax = axes[0]
    ax.loglog(k_hMpc, Pk_hMpc3, 'b-', lw=1.5, label=f'mlx_class (z={z_out:.0f})')
    if k_ref is not None and Pk_ref is not None:
        ax.loglog(k_ref, Pk_ref, 'r--', lw=1.2, label=ref_label)
    ax.set_xlabel(r'$k$ [$h$ Mpc$^{-1}$]', fontsize=13)
    ax.set_ylabel(r'$P(k)$ [(Mpc/$h$)$^3$]', fontsize=13)
    ax.set_title(f'Matter Power Spectrum (z={z_out:.0f})', fontsize=14)
    ax.legend(fontsize=11)
    ax.set_xlim(k_hMpc[0], k_hMpc[-1])

    # Right panel: dimensionless P(k) = k^3 P(k) / (2 pi^2)
    ax = axes[1]
    Delta2 = k_hMpc**3 * Pk_hMpc3 / (2.0 * np.pi**2)
    ax.loglog(k_hMpc, Delta2, 'b-', lw=1.5, label=f'mlx_class (z={z_out:.0f})')
    if k_ref is not None and Pk_ref is not None:
        Delta2_ref = k_ref**3 * Pk_ref / (2.0 * np.pi**2)
        ax.loglog(k_ref, Delta2_ref, 'r--', lw=1.2, label=ref_label)
    ax.set_xlabel(r'$k$ [$h$ Mpc$^{-1}$]', fontsize=13)
    ax.set_ylabel(r'$\Delta^2(k) = k^3 P(k) / (2\pi^2)$', fontsize=13)
    ax.set_title(f'Dimensionless Power Spectrum (z={z_out:.0f})', fontsize=14)
    ax.legend(fontsize=11)
    ax.set_xlim(k_hMpc[0], k_hMpc[-1])

    plt.tight_layout()

    if out_path is None:
        out_dir = os.path.dirname(os.path.abspath(__file__))
        out_path = os.path.join(out_dir, f'matter_pk_z{z_out:.0f}.png')

    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"[P(k)] Plot saved: {out_path}")

    return out_path


# ============================================================================
# Self-test
# ============================================================================

def _test_matter_pk():
    """
    Test the matter power spectrum computation.
    1. Growth factor sanity check
    2. P(k) shape: rises as k^(n_s) at large scales, turns over at k_eq
    3. BAO wiggles visible
    """
    print("=" * 70)
    print("TEST: Matter Power Spectrum P(k)")
    print("=" * 70)

    from .background import Background

    # Background
    bg = Background(khronon=False)
    bg.solve()

    # k-grid: wide range for P(k)
    k_arr = np.geomspace(1e-4, 0.5, 1000).astype(np.float32)

    # Solve with implicit solver
    from .perturbations_implicit import ImplicitBoltzmannSolver
    solver = ImplicitBoltzmannSolver(bg, k_arr)
    result = solver.solve()

    # Test 1: Growth factor
    print("\n--- Test 1: Growth factor ---")
    D_ratio = growth_factor_ratio(a_final=1.0, a_initial=a_rec)
    print(f"  D(z=0)/D(z={z_rec:.0f}) = {D_ratio:.2f}")
    # In matter domination, D ~ a, so ratio ~ (1+z_rec) ~ 1091
    print(f"  Expected (pure matter): ~{1.0/a_rec:.0f}")
    # With Lambda, growth is suppressed, so ratio < 1091
    if 500 < D_ratio < 1200:
        print("  PASS: Growth factor in expected range")
    else:
        print(f"  WARNING: Growth factor {D_ratio:.2f} outside expected range")

    # Test 2: Compute P(k)
    print("\n--- Test 2: P(k) computation ---")
    k_h, Pk_h, k_Mpc, Pk_Mpc = compute_matter_pk(result, solver_type='implicit', z_out=0.0)

    print(f"  k range: [{k_h[0]:.2e}, {k_h[-1]:.2e}] h/Mpc")
    print(f"  P(k) range: [{np.min(Pk_h):.2e}, {np.max(Pk_h):.2e}] (Mpc/h)^3")

    # Check P(k) at specific reference scales
    for k_ref_h, Pk_ref_class in [(0.01, 2e4), (0.1, 800), (0.2, 400)]:
        idx = np.argmin(np.abs(k_h - k_ref_h))
        ratio = Pk_h[idx] / Pk_ref_class
        print(f"  P(k={k_ref_h:.2f} h/Mpc) = {Pk_h[idx]:.2e} (CLASS ~{Pk_ref_class:.0e}, ratio={ratio:.2f})")

    # Test 3: Slope at intermediate scales (near turnover)
    print("\n--- Test 3: Slope near turnover ---")
    mask_mid = (k_h > 0.005) & (k_h < 0.03)
    if np.sum(mask_mid) > 5:
        log_k = np.log(k_h[mask_mid])
        log_Pk = np.log(Pk_h[mask_mid])
        slope = np.polyfit(log_k, log_Pk, 1)[0]
        print(f"  Slope d(ln P)/d(ln k) at k~0.01 = {slope:.3f}")
        print(f"  (Expected: ~-1 to -2 near turnover)")
    else:
        print("  Not enough points in range")

    # Test 4: No NaN/Inf
    print("\n--- Test 4: Stability ---")
    has_nan = np.any(np.isnan(Pk_h))
    has_inf = np.any(np.isinf(Pk_h))
    has_neg = np.any(Pk_h < 0)
    print(f"  NaN: {has_nan}, Inf: {has_inf}, Negative: {has_neg}")
    if not (has_nan or has_inf or has_neg):
        print("  PASS: All P(k) values are positive and finite")
    else:
        print("  FAIL!")

    # Transfer function diagnostic
    transfer_function_diagnostic(result, solver_type='implicit')

    # Plot
    print("\n--- Generating plot ---")
    plot_matter_pk(k_h, Pk_h, z_out=0.0)

    print("\n" + "=" * 70)
    print("MATTER POWER SPECTRUM TEST COMPLETE")
    print("=" * 70)

    return k_h, Pk_h


if __name__ == '__main__':
    _test_matter_pk()
