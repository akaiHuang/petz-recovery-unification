"""
pk_lensing_sync.py -- Matter power spectrum P(k) and CMB lensing
                      from the synchronous gauge Boltzmann solver.

P(k) from sync gauge:
  1. Solve each k-mode with solve_single_k to get delta_c, delta_b at tau_rec
  2. Grow to z=0 using the Heath (1977) growth factor integral
  3. Compute P(k) = (2 pi^2 / k^3) * P_R(k) * |T(k)|^2

CMB lensing from sync gauge:
  1. Use Phi_N + Psi_N from run_sync_solver output (already on snapshot grid)
  2. Compute the Weyl potential power spectrum P_{Phi+Psi}(k, z) at each chi
  3. Limber integral -> C_l^{phi phi}
  4. Flat-sky convolution -> lensed C_l^TT

Author: Sheng-Kai Huang, 2026
"""
import numpy as np
import time
import sys
from scipy.integrate import quad
from scipy.interpolate import interp1d, CubicSpline

from .background import (
    Background, A_s, n_s, k_pivot, T_CMB,
    Omega_b, Omega_r, Omega_L,
    H0_Mpc, h as h_param, z_rec, a_rec,
)
from .perturbations_sync import (
    IDX_ETA, IDX_DELTA_C, IDX_DELTA_B, IDX_THETA_B, IDX_FG_START,
    L_GAMMA_MAX, L_NU_MAX_SYNC,
    diagnose_metric, gauge_transform,
)
from .solver_sync import solve_single_k
from .matter_pk import growth_factor_integral, growth_factor_ratio, _compute_sigma8


# ============================================================================
# P(k) from synchronous gauge solver
# ============================================================================

def compute_pk_from_sync(bg, k_arr=None, N_k=200, k_min=1e-4, k_max=0.5,
                         z_out=0.0, apply_silk=True,
                         method='Radau', rtol=1e-6, atol=1e-9,
                         verbose=True):
    """
    Compute the matter power spectrum P(k) using the synchronous gauge solver.

    Solves perturbation equations for each k-mode up to tau_rec, then extracts
    delta_c and delta_b to build the matter transfer function.  The transfer
    function is grown to z_out using the linear growth factor.

    Parameters
    ----------
    bg : Background (solved)
    k_arr : array or None
        If given, use these k values.  Otherwise, generate log-spaced grid.
    N_k : int
        Number of k-modes (used if k_arr is None).
    k_min, k_max : float
        k range in Mpc^-1 (used if k_arr is None).
    z_out : float
        Output redshift (default 0).
    apply_silk : bool
        Apply Silk damping envelope.
    method : str
        ODE method for solve_ivp.
    rtol, atol : float
        ODE tolerances.
    verbose : bool

    Returns
    -------
    dict with keys:
        'k_hMpc'    : k in h/Mpc
        'Pk_hMpc3'  : P(k) in (Mpc/h)^3
        'k_Mpc'     : k in Mpc^-1
        'Pk_Mpc3'   : P(k) in Mpc^3
        'sigma8'    : sigma_8 value
        'T_k'       : transfer function (unnormalized)
        'delta_c'   : CDM density contrast at tau_rec
        'delta_b'   : baryon density contrast at tau_rec
        'D_ratio'   : growth factor D(z_out)/D(z_rec)
    """
    t0 = time.time()

    if k_arr is None:
        k_arr = np.geomspace(k_min, k_max, N_k)
    else:
        k_arr = np.asarray(k_arr, dtype=np.float64)

    N_k_actual = len(k_arr)
    Omega_c_eff = bg.Omega_cdm
    Omega_m_eff = Omega_b + Omega_c_eff

    if verbose:
        print("=" * 60)
        print("P(k) from Synchronous Gauge Solver")
        print("=" * 60)
        print(f"  N_k = {N_k_actual}, k = [{k_arr[0]:.1e}, {k_arr[-1]:.1e}] Mpc^-1")
        print(f"  Omega_b = {Omega_b:.4f}, Omega_c = {Omega_c_eff:.4f}")
        sys.stdout.flush()

    # We only need to integrate to slightly past tau_rec
    tau_end = bg.tau_rec * 1.05

    # Background quantities at tau_rec for gauge transformation
    calH_rec = float(bg.calH_at_tau(np.array([bg.tau_rec]))[0])
    a_rec_val = float(bg.a_at_tau(np.array([bg.tau_rec]))[0])

    # Solve each k-mode and extract NEWTONIAN GAUGE delta_c, delta_b at tau_rec.
    # The sync gauge delta_c/delta_b are gauge-dependent; we must transform them
    # to the Newtonian gauge for the physical matter power spectrum.
    #
    # Gauge transformation (Ma & Bertschinger 1995):
    #   delta_c_N = delta_c_S + 3 calH alpha   (CDM: w=0, rho'/rho = -3calH)
    #   delta_b_N = delta_b_S + 3 calH alpha   (baryons: w=0 at recombination)
    # where alpha = (h' + 6 eta') / (2 k^2).
    delta_c_N = np.zeros(N_k_actual)
    delta_b_N = np.zeros(N_k_actual)
    delta_c_S = np.zeros(N_k_actual)
    delta_b_S = np.zeros(N_k_actual)
    n_failed = 0

    for ik, k in enumerate(k_arr):
        sol = solve_single_k(k, bg, tau_end,
                             L_GAMMA_MAX, L_NU_MAX_SYNC,
                             method, rtol, atol)
        if not sol.success:
            n_failed += 1
            if verbose and n_failed <= 3:
                print(f"  WARNING: k={k:.4e} failed: {sol.message}")
            continue

        # Evaluate at tau_rec using dense output
        y_rec = sol.sol(bg.tau_rec)
        delta_c_S[ik] = y_rec[IDX_DELTA_C]
        delta_b_S[ik] = y_rec[IDX_DELTA_B]

        # Gauge transformation: sync -> Newtonian
        h_prime, eta_prime = diagnose_metric(
            y_rec, k, calH_rec, a_rec_val, L_GAMMA_MAX, L_NU_MAX_SYNC)
        alpha = (h_prime + 6.0 * eta_prime) / (2.0 * k * k)

        # CDM and baryons are pressureless: rho'/rho = -3 calH
        # delta_N = delta_S - (rho'/rho) alpha = delta_S + 3 calH alpha
        delta_c_N[ik] = y_rec[IDX_DELTA_C] + 3.0 * calH_rec * alpha
        delta_b_N[ik] = y_rec[IDX_DELTA_B] + 3.0 * calH_rec * alpha

        if verbose and (ik + 1) % 50 == 0:
            elapsed = time.time() - t0
            eta_est = elapsed / (ik + 1) * N_k_actual
            print(f"  [{ik+1}/{N_k_actual}] elapsed={elapsed:.1f}s, ETA={eta_est:.0f}s")
            sys.stdout.flush()

    t_pert = time.time() - t0
    if verbose:
        print(f"  Perturbations: {t_pert:.1f}s ({n_failed} failed)")

    # Matter transfer function.
    # At subhorizon scales (k >> aH), the gauge transformation is negligible
    # (alpha ~ 0) and delta_c_N ≈ delta_c_S.  Use delta_c_S as the transfer
    # function: it represents the CDM density contrast in the CDM rest frame,
    # which is the physical late-time matter distribution.
    # Baryons fall into CDM wells after recombination (delta_b -> delta_c at
    # z << z_rec), avoiding deep BAO zeros from baryon oscillations.
    T_k = delta_c_S.copy()

    # Growth factor from z_rec to z_out
    a_out = 1.0 / (1.0 + z_out)
    Omega_L_eff = 1.0 - Omega_m_eff - Omega_r
    D_ratio = growth_factor_ratio(
        a_final=a_out, a_initial=a_rec,
        Omega_m_val=Omega_m_eff, Omega_L_val=Omega_L_eff, Omega_r_val=Omega_r)

    if verbose:
        print(f"  D(z={z_out:.1f})/D(z={z_rec:.0f}) = {D_ratio:.2f}")

    T_k_grown = T_k * D_ratio

    # Silk damping
    if apply_silk:
        silk = np.exp(-(k_arr / bg.k_D)**2)
        T_k_grown = T_k_grown * silk

    # Primordial power spectrum: P_R(k) = A_s * (k/k_pivot)^(n_s-1)
    P_R = A_s * (k_arr / k_pivot) ** (n_s - 1.0)

    # P(k) = (2 pi^2 / k^3) * P_R(k) * (C/R)^2 * |T(k)|^2
    #
    # NORMALIZATION:
    # The sync gauge uses C = 1 as the initial curvature perturbation.
    # The primordial comoving curvature perturbation R = -(3/2) C, so
    # |C/R|^2 = (2/3)^2.  The primordial power spectrum P_R(k) = A_s * ...
    # is per unit R^2, so we need (C/R)^2 to convert delta_c (per unit C)
    # to delta_c per unit R.
    #
    # This is the same normalization used in the C_l calculation
    # (solver_sync.py line 379).
    norm_sync = (2.0 / 3.0) ** 2   # (C/R)^2
    Pk_Mpc3 = norm_sync * (2.0 * np.pi**2 / k_arr**3) * P_R * T_k_grown**2

    # Sigma_8
    sigma8 = _compute_sigma8(k_arr, Pk_Mpc3, h_param)

    if verbose:
        print(f"  sigma_8 = {sigma8:.4f} (Planck: 0.811)")

    # Convert to h/Mpc units
    k_hMpc = k_arr / h_param
    Pk_hMpc3 = Pk_Mpc3 * h_param**3

    t_total = time.time() - t0
    if verbose:
        print(f"  Total P(k) time: {t_total:.1f}s")

    return {
        'k_hMpc': k_hMpc,
        'Pk_hMpc3': Pk_hMpc3,
        'k_Mpc': k_arr,
        'Pk_Mpc3': Pk_Mpc3,
        'sigma8': sigma8,
        'T_k': T_k,
        'delta_c_N': delta_c_N,
        'delta_b_N': delta_b_N,
        'delta_c_S': delta_c_S,
        'delta_b_S': delta_b_S,
        'D_ratio': D_ratio,
    }


# ============================================================================
# Lensing potential C_l^{phi phi} from sync gauge Phi_N, Psi_N
# ============================================================================

def compute_cl_phiphi_from_sync(sync_result, ell_max=3000, verbose=True):
    """
    Compute C_l^{phi phi} from the synchronous gauge solver output,
    using the actual Phi_N + Psi_N Weyl potential (no Eisenstein-Hu fit).

    The lensing potential power spectrum via Limber approximation:

      C_l^{pp} = int_0^{chi_*} dchi  W(chi)^2  P_{Weyl}(k = (l+0.5)/chi, chi)

    where W(chi) = 2 (chi_* - chi) / (chi_* chi) is the lensing kernel
    and P_{Weyl}(k, chi) is the power spectrum of (Phi_N + Psi_N)/2.

    Instead of a continuous P_{Weyl}, we build it from the discrete k-grid
    of the sync solver by interpolating |Phi_N + Psi_N|^2 in (k, tau) space.

    Parameters
    ----------
    sync_result : dict
        Output from run_sync_solver, with keys:
        'k_arr', 'bg', 'Phi_N', 'Psi_N', and tau_all from snapshot grid.
    ell_max : int
    verbose : bool

    Returns
    -------
    ell_phi : array (int, 2..ell_max)
    Cl_phiphi : array
    """
    t0 = time.time()

    bg = sync_result['bg']
    k_arr = sync_result['k_arr']
    Phi_N = sync_result['Phi_N']   # (N_k, N_snap)
    Psi_N = sync_result['Psi_N']   # (N_k, N_snap)
    N_k = len(k_arr)

    Omega_c_eff = bg.Omega_cdm
    Omega_m_eff = Omega_b + Omega_c_eff
    Omega_L_eff = 1.0 - Omega_m_eff - Omega_r

    # Weyl potential at each snapshot: (Phi + Psi) / 2
    Weyl = 0.5 * (Phi_N + Psi_N)   # (N_k, N_snap)

    # We need the snapshot tau grid.  Reconstruct from the solver.
    # The sync solver stores tau_all implicitly via build_snapshot_grid.
    # Re-build it here.
    from .solver_sync import build_snapshot_grid
    tau_vis, tau_late, tau_all = build_snapshot_grid(
        bg, N_vis=60, N_early_isw=30, N_late_isw=20)
    N_snap = len(tau_all)

    if Phi_N.shape[1] != N_snap:
        raise ValueError(
            f"Snapshot grid mismatch: Phi_N has {Phi_N.shape[1]} columns "
            f"but rebuild gives {N_snap} snapshots. "
            f"Make sure sync_result comes from a standard run_sync_solver call.")

    chi_star = bg.D_A
    chi_all = bg.tau_0 - tau_all       # comoving distance at each snapshot
    a_all = bg.a_at_tau(tau_all)

    if verbose:
        print(f"\n[Lensing-sync] Computing C_l^{{phi phi}} from Weyl potential")
        print(f"  N_k = {N_k}, N_snap = {N_snap}, ell_max = {ell_max}")
        sys.stdout.flush()

    # Build the Weyl power spectrum P_Weyl(k) at each snapshot.
    # For each snapshot i, P_Weyl(k) ~ |Weyl(k, tau_i)|^2 is what the solver
    # gives for unit initial amplitude.  We multiply by the primordial P_R(k)
    # and the correct normalization.
    #
    # The sync gauge has C = 1, R = -(3/2)C, so Phi_N per unit R is
    # Phi_N_solver * (2/3).  The Weyl potential per unit R is:
    # W(k,tau) = (2/3) * Weyl_solver(k,tau)
    #
    # P_{Weyl}(k, chi) = (2pi^2 / k^3) * P_R(k) * |W_per_R|^2
    #                   = (2pi^2 / k^3) * P_R(k) * (4/9) * |Weyl_solver|^2
    #
    # The Limber integral:
    # C_l^{pp} = int dchi [2(chi*-chi)/(chi* chi)]^2 * P_Weyl(k=(l+0.5)/chi, chi)

    P_R = A_s * (k_arr / k_pivot) ** (n_s - 1.0)
    norm_R2 = (2.0 / 3.0) ** 2  # (C/R)^2

    # For each snapshot, compute (2pi^2/k^3) * P_R * (4/9) * Weyl^2 as a function of k.
    # Then interpolate to the Limber k = (l+0.5)/chi for each ell.

    ell_phi = np.arange(2, ell_max + 1)
    ell_f = ell_phi.astype(float)
    Cl_pp = np.zeros(len(ell_phi))

    # Only use snapshots with chi > 0 and chi < chi_star
    valid_snap = (chi_all > 10.0) & (chi_all < chi_star - 10.0)
    snap_indices = np.where(valid_snap)[0]

    if len(snap_indices) < 5:
        print("[Lensing-sync] WARNING: Too few valid snapshots for lensing")
        return ell_phi, Cl_pp

    # Lensing kernel at each snapshot
    chi_valid = chi_all[snap_indices]
    W_lens = 2.0 * (chi_star - chi_valid) / (chi_star * chi_valid)

    # Trapezoidal weights for the chi integral (note: chi decreases with snapshot index
    # since tau increases -> chi = tau_0 - tau decreases).  Sort by chi for integration.
    sort_idx = np.argsort(chi_valid)
    chi_sorted = chi_valid[sort_idx]
    snap_sorted = snap_indices[sort_idx]
    W_sorted = W_lens[sort_idx]

    dchi = np.zeros(len(chi_sorted))
    dchi[0] = 0.5 * (chi_sorted[1] - chi_sorted[0])
    dchi[-1] = 0.5 * (chi_sorted[-1] - chi_sorted[-2])
    for i in range(1, len(chi_sorted) - 1):
        dchi[i] = 0.5 * (chi_sorted[i + 1] - chi_sorted[i - 1])

    # Build log-log interpolators for Weyl^2 at each snapshot
    log_k = np.log(k_arr)

    for ii in range(len(chi_sorted)):
        chi = chi_sorted[ii]
        isnap = snap_sorted[ii]
        w_chi = W_sorted[ii]

        # Weyl potential squared at this snapshot for each k-mode
        weyl_sq = Weyl[:, isnap] ** 2

        # P_Weyl(k) for unit-R perturbation
        P_weyl_k = norm_R2 * (2.0 * np.pi**2 / k_arr**3) * P_R * weyl_sq

        # Interpolate in log-log space
        log_P = np.log(np.maximum(P_weyl_k, 1e-200))
        f_pk = interp1d(log_k, log_P, kind='linear',
                        fill_value=-460.0, bounds_error=False)

        # Limber: k = (l + 0.5) / chi
        k_ell = (ell_f + 0.5) / chi

        # Only evaluate within the solver k-range
        in_range = (k_ell >= k_arr[0]) & (k_ell <= k_arr[-1])
        P_ell = np.zeros_like(k_ell)
        if np.any(in_range):
            P_ell[in_range] = np.exp(f_pk(np.log(k_ell[in_range])))

        # Accumulate: dchi * W^2 * P_Weyl / chi^2
        # (the chi^2 cancels part of the kernel but we keep the standard form)
        Cl_pp += dchi[ii] * w_chi**2 * P_ell

    t_elapsed = time.time() - t0
    if verbose:
        print(f"[Lensing-sync] C_l^{{phi phi}} done in {t_elapsed:.2f}s")
        for l_diag in [50, 100, 500, 1000]:
            idx = l_diag - 2
            if idx < len(Cl_pp) and Cl_pp[idx] > 0:
                print(f"  l={l_diag}: l^4 C_l^{{pp}} = {l_diag**4 * Cl_pp[idx]:.4e}")

    return ell_phi, Cl_pp


# ============================================================================
# Lensed C_l^TT from sync gauge
# ============================================================================

def apply_lensing_sync(sync_result, ell_max_phi=None, verbose=True):
    """
    Apply CMB lensing to the sync gauge C_l^TT spectrum.

    Uses the Weyl potential (Phi_N + Psi_N) from the sync solver to compute
    the lensing potential, then applies the flat-sky lensing convolution
    from lensing.py.

    Parameters
    ----------
    sync_result : dict
        Output from run_sync_solver.
    ell_max_phi : int or None
        Maximum ell for the lensing potential (default: max ell + 500).
    verbose : bool

    Returns
    -------
    dict with keys:
        'ell_lensed'  : multipoles (dense)
        'Dl_lensed'   : lensed D_l in muK^2 (dense)
        'ell_unlensed': original ell
        'Dl_unlensed' : original D_l in muK^2
        'ell_phi'     : lensing potential ell
        'Cl_phiphi'   : lensing potential C_l
        'sigma_lens'  : (sigma^2, d_rms_arcmin)
    """
    from .lensing import lens_cl_tt, compute_sigma_lens

    ell_in = sync_result['ell']
    Dl_in = sync_result['Dl']
    ell_f = ell_in.astype(float)

    if ell_max_phi is None:
        ell_max_phi = int(ell_in[-1]) + 500

    if verbose:
        print(f"\n{'='*60}")
        print("CMB LENSING FROM SYNC GAUGE (Weyl potential)")
        print(f"{'='*60}")

    # Step 1: Lensing potential from Phi_N + Psi_N
    ell_phi, Cl_pp = compute_cl_phiphi_from_sync(
        sync_result, ell_max=ell_max_phi, verbose=verbose)

    # Step 2: Diagnostics
    sigma2, d_rms = compute_sigma_lens(ell_phi, Cl_pp)
    if verbose:
        print(f"[Lensing-sync] sigma_lens^2 = {sigma2:.6e}")
        print(f"[Lensing-sync] RMS deflection = {d_rms:.2f} arcmin (Planck: ~2.7)")

    # Step 3: Convert D_l to C_l
    Cl_in = Dl_in / (ell_f * (ell_f + 1.0) / (2.0 * np.pi) * (T_CMB * 1e6)**2)

    # Step 4: Apply flat-sky lensing convolution
    ell_lensed, Cl_lensed = lens_cl_tt(ell_in, Cl_in, ell_phi, Cl_pp)

    # Step 5: Convert back to D_l
    ell_l_f = ell_lensed.astype(float)
    Dl_lensed = ell_l_f * (ell_l_f + 1.0) / (2.0 * np.pi) * Cl_lensed * (T_CMB * 1e6)**2

    if verbose:
        print(f"[Lensing-sync] Lensed D_l range: "
              f"[{np.min(Dl_lensed):.0f}, {np.max(Dl_lensed):.0f}] muK^2")

    return {
        'ell_lensed': ell_lensed,
        'Dl_lensed': Dl_lensed,
        'ell_unlensed': ell_in,
        'Dl_unlensed': Dl_in,
        'ell_phi': ell_phi,
        'Cl_phiphi': Cl_pp,
        'sigma_lens': (sigma2, d_rms),
    }


# ============================================================================
# Plotting utilities
# ============================================================================

def plot_pk_comparison(pk_sync, pk_ref=None, ref_label='Newtonian gauge',
                       out_path=None):
    """
    Plot P(k) from sync gauge, optionally compared with another result.

    Parameters
    ----------
    pk_sync : dict from compute_pk_from_sync
    pk_ref : dict or None
        If dict, must have 'k_hMpc' and 'Pk_hMpc3' keys.
    ref_label : str
    out_path : str or None
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import os

    k = pk_sync['k_hMpc']
    Pk = pk_sync['Pk_hMpc3']

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Left: P(k)
    ax = axes[0]
    ax.loglog(k, Pk, 'b-', lw=1.5, label=f"Sync gauge ($\\sigma_8$={pk_sync['sigma8']:.3f})")
    if pk_ref is not None:
        ax.loglog(pk_ref['k_hMpc'], pk_ref['Pk_hMpc3'], 'r--', lw=1.2,
                  label=ref_label)
    ax.set_xlabel(r'$k$ [$h$ Mpc$^{-1}$]', fontsize=13)
    ax.set_ylabel(r'$P(k)$ [(Mpc/$h$)$^3$]', fontsize=13)
    ax.set_title('Matter Power Spectrum', fontsize=14)
    ax.legend(fontsize=11)

    # Right: ratio (if reference available)
    ax = axes[1]
    if pk_ref is not None:
        k_common = np.geomspace(
            max(k[0], pk_ref['k_hMpc'][0]),
            min(k[-1], pk_ref['k_hMpc'][-1]),
            500)
        f_sync = interp1d(np.log(k), np.log(np.maximum(Pk, 1e-50)),
                          kind='cubic', fill_value='extrapolate')
        f_ref = interp1d(np.log(pk_ref['k_hMpc']),
                         np.log(np.maximum(pk_ref['Pk_hMpc3'], 1e-50)),
                         kind='cubic', fill_value='extrapolate')
        ratio = np.exp(f_sync(np.log(k_common)) - f_ref(np.log(k_common)))
        ax.semilogx(k_common, ratio, 'b-', lw=1.5)
        ax.axhline(1.0, color='k', ls='--', lw=0.5)
        ax.fill_between(k_common, 0.9, 1.1, color='green', alpha=0.1)
        ax.set_ylabel('Sync / Reference', fontsize=13)
        ax.set_title('P(k) Ratio', fontsize=14)
    else:
        Delta2 = k**3 * Pk / (2.0 * np.pi**2)
        ax.loglog(k, Delta2, 'b-', lw=1.5)
        ax.set_ylabel(r'$\Delta^2(k)$', fontsize=13)
        ax.set_title('Dimensionless Power Spectrum', fontsize=14)
    ax.set_xlabel(r'$k$ [$h$ Mpc$^{-1}$]', fontsize=13)

    plt.tight_layout()

    if out_path is None:
        out_dir = os.path.dirname(os.path.abspath(__file__))
        out_path = os.path.join(out_dir, 'pk_sync_gauge.png')

    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"[P(k)-sync] Plot saved: {out_path}")
    return out_path


def plot_lensing_sync(lens_result, class_data=None, out_path=None, ell_max=2500):
    """
    Plot lensed vs unlensed D_l from sync gauge, with optional CLASS comparison.

    Parameters
    ----------
    lens_result : dict from apply_lensing_sync
    class_data : dict or None, from lensing.load_class_lensed
    out_path : str or None
    ell_max : int
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import os

    ell_u = lens_result['ell_unlensed']
    Dl_u = lens_result['Dl_unlensed']
    ell_l = lens_result['ell_lensed']
    Dl_l = lens_result['Dl_lensed']
    ell_phi = lens_result['ell_phi']
    Cl_pp = lens_result['Cl_phiphi']

    fig, axes = plt.subplots(3, 1, figsize=(14, 14),
                             gridspec_kw={'height_ratios': [3, 1, 2]})
    fig.subplots_adjust(hspace=0.25)

    # Panel 1: D_l spectra
    ax = axes[0]
    m_u = ell_u <= ell_max
    m_l = ell_l <= ell_max
    ax.plot(ell_u[m_u], Dl_u[m_u], 'b-', lw=1.0, alpha=0.5, label='Unlensed (sync)')
    ax.plot(ell_l[m_l], Dl_l[m_l], 'r-', lw=1.5, label='Lensed (sync)')
    if class_data is not None:
        m_c = class_data['ell'] <= ell_max
        ax.plot(class_data['ell'][m_c], class_data['TT'][m_c],
                'k--', lw=1.2, alpha=0.7, label='CLASS lensed')
    ax.set_ylabel(r'$\mathcal{D}_\ell^{TT}$ [$\mu$K$^2$]', fontsize=13)
    ax.set_title('CMB Lensing from Sync Gauge Weyl Potential', fontsize=14)
    ax.legend(fontsize=11)
    ax.set_xlim(2, ell_max)
    ax.set_ylim(0, None)

    # Panel 2: Lensing residual
    ax2 = axes[1]
    ell_common = np.arange(30, min(ell_max, int(ell_u[-1]), int(ell_l[-1])) + 1)
    f_u = interp1d(ell_u.astype(float), Dl_u, kind='cubic', fill_value='extrapolate')
    f_l = interp1d(ell_l.astype(float), Dl_l, kind='cubic', fill_value='extrapolate')
    Du = f_u(ell_common)
    Dl = f_l(ell_common)
    residual = np.where(Du > 50, (Dl - Du) / Du * 100, 0.0)
    ax2.plot(ell_common, residual, 'r-', lw=0.8)
    ax2.axhline(0, color='k', ls='--', lw=0.5)
    ax2.fill_between(ell_common, -10, 10, color='green', alpha=0.08)
    ax2.set_ylabel('(Lensed-Unlensed)/Unlensed [%]', fontsize=11)
    ax2.set_ylim(-15, 15)
    ax2.set_xlabel(r'Multipole $\ell$', fontsize=13)

    # Panel 3: C_l^{phi phi}
    ax3 = axes[2]
    ell_pp = ell_phi.astype(float)
    Dl_pp = ell_pp**2 * (ell_pp + 1)**2 * Cl_pp / (2.0 * np.pi)
    m_pp = ell_pp <= ell_max
    ax3.semilogy(ell_pp[m_pp], Dl_pp[m_pp], 'b-', lw=1.5,
                 label=r'Sync gauge $[\ell(\ell+1)]^2 C_\ell^{\phi\phi}/(2\pi)$')
    if class_data is not None and 'phiphi' in class_data:
        ell_c = class_data['ell'].astype(float)
        Cpp_c = class_data['phiphi']
        Dl_pp_c = ell_c**2 * (ell_c + 1)**2 * Cpp_c / (2.0 * np.pi)
        m_cc = ell_c <= ell_max
        ax3.semilogy(ell_c[m_cc], Dl_pp_c[m_cc], 'k--', lw=1.2,
                     alpha=0.7, label='CLASS')
    ax3.set_xlabel(r'Multipole $\ell$', fontsize=13)
    ax3.set_ylabel(r'$[\ell(\ell+1)]^2 C_\ell^{\phi\phi}/(2\pi)$', fontsize=12)
    ax3.set_title('Lensing Potential from Sync Gauge Weyl Potential', fontsize=13)
    ax3.legend(fontsize=10)
    ax3.set_xlim(2, ell_max)

    plt.tight_layout()

    if out_path is None:
        out_dir = os.path.dirname(os.path.abspath(__file__))
        out_path = os.path.join(out_dir, 'lensing_sync_gauge.png')

    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"[Lensing-sync] Plot saved: {out_path}")
    return out_path


# ============================================================================
# Self-test / CLI entry point
# ============================================================================

def main():
    """
    Run P(k) and lensing from the sync gauge solver.

    Usage: python -m mlx_class.pk_lensing_sync [--N_k_pk 200] [--N_k_cl 50] [--skip_pk] [--skip_lensing]
    """
    import argparse
    import os

    parser = argparse.ArgumentParser(
        description='P(k) and CMB lensing from sync gauge solver')
    parser.add_argument('--N_k_pk', type=int, default=200,
                        help='Number of k-modes for P(k)')
    parser.add_argument('--N_k_cl', type=int, default=50,
                        help='Number of k-modes for C_l (lensing source)')
    parser.add_argument('--skip_pk', action='store_true',
                        help='Skip P(k) computation')
    parser.add_argument('--skip_lensing', action='store_true',
                        help='Skip lensing computation')
    parser.add_argument('--method', type=str, default='Radau')
    args = parser.parse_args()

    print("=" * 70)
    print("P(k) AND CMB LENSING FROM SYNCHRONOUS GAUGE SOLVER")
    print("=" * 70)

    # ---------------------------------------------------------------
    # Background
    # ---------------------------------------------------------------
    bg = Background(khronon=False, recombination='recfast')
    bg.solve()

    pk_sync = None
    lens_result = None

    # ---------------------------------------------------------------
    # Task 1: Matter power spectrum P(k)
    # ---------------------------------------------------------------
    if not args.skip_pk:
        print(f"\n{'='*60}")
        print("TASK 1: Matter Power Spectrum P(k)")
        print(f"{'='*60}")

        pk_sync = compute_pk_from_sync(
            bg, N_k=args.N_k_pk,
            k_min=1e-4, k_max=0.5,
            z_out=0.0, apply_silk=True,
            method=args.method, verbose=True)

        # Comparison with Newtonian gauge (if implicit solver is available)
        pk_ref = None
        try:
            from .perturbations_implicit import ImplicitBoltzmannSolver
            from .matter_pk import compute_matter_pk as compute_pk_newton

            k_newton = np.geomspace(1e-4, 0.5, 1000).astype(np.float32)
            solver = ImplicitBoltzmannSolver(bg, k_newton)
            result = solver.solve()
            k_h, Pk_h, k_Mpc, Pk_Mpc = compute_pk_newton(
                result, solver_type='implicit', z_out=0.0)
            pk_ref = {'k_hMpc': k_h, 'Pk_hMpc3': Pk_h}
            print(f"\n  Newtonian gauge sigma_8 for comparison (from implicit solver)")
        except Exception as e:
            print(f"\n  [INFO] Could not run Newtonian gauge comparison: {e}")

        plot_pk_comparison(pk_sync, pk_ref)

        # Print comparison at reference scales
        print(f"\n  P(k) at reference scales:")
        k_h = pk_sync['k_hMpc']
        Pk_h = pk_sync['Pk_hMpc3']
        for k_ref_val, Pk_class_approx in [(0.01, 2e4), (0.1, 800), (0.2, 350)]:
            idx = np.argmin(np.abs(k_h - k_ref_val))
            print(f"    k={k_ref_val:.2f} h/Mpc: P(k)={Pk_h[idx]:.2e} "
                  f"(CLASS ~{Pk_class_approx:.0e})")

    # ---------------------------------------------------------------
    # Task 2: CMB lensing
    # ---------------------------------------------------------------
    if not args.skip_lensing:
        print(f"\n{'='*60}")
        print("TASK 2: CMB Lensing from Sync Gauge")
        print(f"{'='*60}")

        # Run sync solver for C_l + lensing data
        from .solver_sync import run_sync_solver

        print(f"\nRunning sync gauge solver for lensing (N_k={args.N_k_cl})...")
        sync_result = run_sync_solver(
            N_k=args.N_k_cl, k_min=3e-4, k_max=0.35,
            method=args.method, verbose=True)

        # Apply lensing
        lens_result = apply_lensing_sync(sync_result, verbose=True)

        # Try loading CLASS reference
        class_data = None
        class_path = ('/Users/akaihuangm1/Desktop/github/'
                      'gdm_class_public/output/running_mu_lcdm_ref_cl_lensed.dat')
        if os.path.exists(class_path):
            from .lensing import load_class_lensed
            class_data = load_class_lensed(class_path)
            print(f"\n  Loaded CLASS lensed reference: {len(class_data['ell'])} multipoles")

        plot_lensing_sync(lens_result, class_data=class_data)

        # Lensing diagnostics
        sigma2, d_rms = lens_result['sigma_lens']
        print(f"\n  Summary:")
        print(f"    RMS deflection: {d_rms:.2f} arcmin (Planck: ~2.7)")
        print(f"    Unlensed D_l peak: {np.max(lens_result['Dl_unlensed']):.0f} muK^2")
        print(f"    Lensed D_l peak: {np.max(lens_result['Dl_lensed']):.0f} muK^2")

    # ---------------------------------------------------------------
    # Summary
    # ---------------------------------------------------------------
    print(f"\n{'='*70}")
    print("COMPLETE")
    if pk_sync is not None:
        print(f"  P(k): sigma_8 = {pk_sync['sigma8']:.4f}")
    if lens_result is not None:
        sigma2, d_rms = lens_result['sigma_lens']
        print(f"  Lensing: RMS deflection = {d_rms:.2f} arcmin")
    print(f"{'='*70}")


if __name__ == '__main__':
    main()
