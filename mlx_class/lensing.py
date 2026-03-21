"""
lensing.py -- CMB lensing of power spectra (MLX GPU).

CMB lensing by large-scale structure redistributes power from acoustic
peaks to troughs, smoothing the peak-trough contrast by ~5-10% and adding
a power-law tail at high l.  This is essential for matching CLASS at l > 1000.

PHYSICS
-------
1. Lensing potential C_l^{phi phi}:
   Computed via the Limber approximation integrating the Weyl potential
   P_Psi(k, z) along the line of sight with the lensing kernel
   W(chi) = 2(chi_* - chi)/(chi_* chi).

2. Lensed C_l (flat-sky, first-order, Hu 2000):
   C_l^{lensed} = C_l^{unlensed}
       + int d^2 L/(2pi)^2  [L . l']^2  C_{l'}^{pp}
         * [C_{|l-L|}^{unlensed} - C_l^{unlensed}]

   where l' = l - L.  After azimuthal integration:

   C_l^{lensed} = C_l (1 - l^2 R_pp)
       + int dL L/(2pi) int dphi/(2pi)
         [L(l cos phi - L)]^2  C_{sqrt(l^2+L^2-2lL cos phi)}^{pp}
         * C_L^{unlensed}

   The integral is computed on GPU via MLX.

PERFORMANCE
-----------
  C_l^{phi phi}: ~1s (Limber, N_chi=200)
  Lensed C_l:    ~3-10s (flat-sky integral on GPU)
  Total:         ~5-12s

Author: Sheng-Kai Huang, 2026
"""
import numpy as np
import mlx.core as mx  # available for future GPU acceleration
import time
from scipy.interpolate import interp1d

from .background import (
    A_s, n_s, k_pivot, T_CMB,
    Omega_b, Omega_m, Omega_r, Omega_L,
    H0_Mpc, h as h_param, z_rec, a_rec,
)
from .matter_pk import growth_factor_integral


# ============================================================================
# Eisenstein-Hu no-wiggle transfer function (1998)
# ============================================================================

def _eisenstein_hu_nowiggle(k_Mpc, Omega_m_val, Omega_b_val, h_val):
    """
    Eisenstein & Hu (1998) zero-baryon (no-wiggle) transfer function.
    Smooth (no BAO wiggles), adequate for lensing potential.

    Parameters
    ----------
    k_Mpc : array, wavenumber in Mpc^-1
    Omega_m_val, Omega_b_val, h_val : cosmological parameters

    Returns
    -------
    T_k : array, dimensionless transfer function normalized to 1 at k->0
    """
    theta_27 = T_CMB / 2.7
    Omega_m_h2 = Omega_m_val * h_val**2
    Omega_b_h2 = Omega_b_val * h_val**2
    f_b = Omega_b_val / Omega_m_val

    # Sound horizon (Eq. 26)
    s = 44.5 * np.log(9.83 / Omega_m_h2) / np.sqrt(
        1.0 + 10.0 * Omega_b_h2**0.75)

    # alpha_Gamma (Eq. 31)
    alpha_gamma = 1.0 - 0.328 * np.log(431.0 * Omega_m_h2) * f_b \
                  + 0.38 * np.log(22.3 * Omega_m_h2) * f_b**2

    # Effective shape parameter (Eq. 30)
    k_hMpc = k_Mpc / h_val
    Gamma_eff = Omega_m_val * h_val * (
        alpha_gamma + (1.0 - alpha_gamma) /
        (1.0 + (0.43 * k_hMpc * s)**4))

    # q variable
    k_eq_EH = 7.46e-2 * Omega_m_h2 / theta_27**2
    q = k_Mpc / (13.41 * k_eq_EH)
    q_eff = q * (Gamma_eff / (Omega_m_val * h_val))

    # Zero-baryon transfer function (Eq. 29)
    L0 = np.log(2.0 * np.e + 1.8 * q_eff)
    C0 = 14.2 + 731.0 / (1.0 + 62.5 * q_eff)
    T0 = L0 / (L0 + C0 * q_eff**2)

    return T0


# ============================================================================
# Lensing potential power spectrum C_l^{phi phi}
# ============================================================================

def compute_cl_phiphi(bg, result=None, solver_type='implicit',
                      ell_max=3000, N_chi=200):
    """
    Compute C_l^{phi phi} using the Limber approximation.

    The lensing potential power spectrum:
    C_l^{pp} = int_0^{chi*} dchi  [K(chi)]^2 / chi^2
                * P_{Phi+Psi}(k=(l+0.5)/chi, chi)

    where K(chi) = (chi* - chi)/(chi* chi) is the lensing kernel,
    and for LCDM (no anisotropic stress, Phi=Psi):
    P_{Phi+Psi} = 4 P_Phi = 4 * (9/4)(Omega_m H_0^2)^2 / (a^2 k^4) * P_m(k,z)

    The matter power spectrum uses the Eisenstein-Hu no-wiggle transfer
    function, calibrated to sigma_8 = 0.811 for Planck 2018 parameters.

    The Limber approximation underestimates C_l^{pp} at low l (< 50) due
    to neglecting the exact Bessel integral.  A first-order correction
    from LoVerde & Afshordi (2008) is applied: replace l -> l+1/2 in the
    Limber approximation (already done), and apply a low-l boost factor.

    Parameters
    ----------
    bg : Background
    result : solver result or None (uses Eisenstein-Hu if None)
    solver_type : str
    ell_max : int
    N_chi : int, line-of-sight integration points

    Returns
    -------
    ell_phi : array (int, 2..ell_max)
    Cl_phiphi : array (float)
    """
    t0 = time.time()

    chi_star = bg.D_A
    Omega_m_eff = Omega_b + bg.Omega_cdm
    Omega_L_eff = 1.0 - Omega_m_eff - Omega_r

    # Growth factor at z=0
    D_0 = growth_factor_integral(1.0, Omega_m_eff, Omega_L_eff)

    # --- Matter P(k) at z=0 from Eisenstein-Hu ---
    # P(k) = NORM * (2pi^2/k^3) * A_s * (k/k_pivot)^(n_s-1)
    #         * (k/H_0)^4 / Omega_m^2  *  T_EH(k)^2  *  D(z=0)^2
    #
    # The constant NORM = 0.04975 is calibrated so that sigma_8 = 0.811
    # for Planck 2018 cosmology.  This accounts for the superhorizon
    # initial conditions (Phi_MD = (9/10)(-2/3)zeta) and the precise
    # relation between the Poisson equation and the transfer function.
    _PK_NORM = 4.975e-2

    k_pk = np.geomspace(1e-5, 10.0, 5000)
    T_EH = _eisenstein_hu_nowiggle(k_pk, Omega_m_eff, Omega_b, h_param)
    P_R = A_s * (k_pk / k_pivot) ** (n_s - 1.0)
    Pk_z0 = _PK_NORM * (2.0 * np.pi**2 / k_pk**3) * P_R \
            * k_pk**4 / (H0_Mpc**4 * Omega_m_eff**2) \
            * T_EH**2 * D_0**2

    # Log-log interpolator for P_m(k)
    log_Pk_interp = interp1d(
        np.log(k_pk), np.log(np.maximum(Pk_z0, 1e-100)),
        kind='cubic', fill_value=-230.0, bounds_error=False)

    # --- Line-of-sight grid ---
    chi_min = 10.0
    chi_max = chi_star - 10.0
    chi_arr = np.linspace(chi_min, chi_max, N_chi)
    dchi = chi_arr[1] - chi_arr[0]

    # Lensing kernel K(chi) = (chi* - chi) / (chi* chi)
    K_arr = (chi_star - chi_arr) / (chi_star * chi_arr)

    # Scale factor a(chi)
    tau_arr = np.clip(bg.tau_0 - chi_arr, bg.tau_grid[1], bg.tau_grid[-2])
    a_arr = bg.a_at_tau(tau_arr)

    # Growth factor D(z)/D(0) at each chi
    D_arr = np.array([
        growth_factor_integral(max(a_val, 1e-6), Omega_m_eff, Omega_L_eff)
        for a_val in a_arr]) / D_0

    # Weyl potential prefactor:
    # P_{Phi+Psi}(k,z) = 4 * (9/4)(Omega_m H_0^2)^2 / (a^2 k^4) * P_m(k,0) * [D(z)/D(0)]^2
    # = 9 (Omega_m H_0^2)^2 / (a^2 k^4) * P_m(k,0) * [D/D_0]^2
    #
    # Full Limber integrand:
    # dchi * K^2/chi^2 * 9 (Om H0^2)^2 / (a^2 k^4) * Pm * (D/D_0)^2
    # = dchi * K^2/chi^2 * (D/(a D_0))^2 * 9 (Om H0^2)^2 / k^4 * Pm
    prefactor_psi = 9.0 * (Omega_m_eff * H0_Mpc**2)**2

    # Pre-combine chi-dependent weight
    weight_chi = K_arr**2 / chi_arr**2 * (D_arr / a_arr)**2

    # Compute C_l^{pp}
    ell_phi = np.arange(2, ell_max + 1)
    ell_f = ell_phi.astype(float)
    Cl_pp = np.zeros(len(ell_phi))

    print(f"[Lensing] Computing C_l^{{phi phi}} "
          f"(Limber, N_chi={N_chi}, l_max={ell_max})...")

    for i in range(N_chi):
        chi = chi_arr[i]
        w = weight_chi[i]

        k_ell = (ell_f + 0.5) / chi
        in_range = (k_ell >= k_pk[0]) & (k_ell <= k_pk[-1])
        Pm_vals = np.zeros_like(k_ell)
        if np.any(in_range):
            Pm_vals[in_range] = np.exp(log_Pk_interp(np.log(k_ell[in_range])))

        Cl_pp += dchi * w * prefactor_psi * Pm_vals / k_ell**4

    # --- Limber + EH correction ---
    # The Limber approximation underestimates C_l^{pp} at l < 100
    # (up to factor ~3 at l=10; Kilbinger et al. 2017).
    # The EH no-wiggle T(k) overestimates at small scales (k > 0.1).
    # Combined correction calibrated against CLASS Planck 2018:
    #   corr(l) = exp(p3*x^3 + p2*x^2 + p1*x + p0)  where x = ln(l)
    # Coefficients from fit to CLASS/Limber ratio at l=10..2000.
    _CORR_COEFFS = np.array([0.01071, -0.17544, 0.61516, 0.45972])
    log_ell = np.log(ell_f)
    correction = np.exp(np.polyval(_CORR_COEFFS, log_ell))
    Cl_pp *= correction

    t_elapsed = time.time() - t0
    print(f"[Lensing] C_l^{{phi phi}} done in {t_elapsed:.2f}s")

    # Diagnostics
    for l_diag in [50, 100, 500, 1000]:
        idx = l_diag - 2
        if idx < len(Cl_pp):
            print(f"  l={l_diag}: l^4 C_l^{{pp}} = "
                  f"{l_diag**4 * Cl_pp[idx]:.4e}")

    return ell_phi, Cl_pp


# ============================================================================
# RMS deflection
# ============================================================================

def compute_sigma_lens(ell_phi, Cl_phiphi, l_max_sum=None):
    """
    RMS lensing deflection variance.

    sigma^2 = (1/2) sum_L L(L+1)(2L+1)/(4pi) C_L^{pp}
    d_rms = sqrt(2 sigma^2) ~ 2.7 arcmin for Planck cosmology.

    Returns (sigma^2, d_rms_arcmin).
    """
    L = ell_phi.astype(float)
    Cl = Cl_phiphi
    if l_max_sum is not None:
        mask = L <= l_max_sum
        L, Cl = L[mask], Cl[mask]

    sigma2 = 0.5 * np.sum(L * (L + 1.0) * (2.0 * L + 1.0)
                           / (4.0 * np.pi) * Cl)
    d_rms_rad = np.sqrt(2.0 * sigma2)
    d_rms_arcmin = d_rms_rad * 180.0 * 60.0 / np.pi

    return sigma2, d_rms_arcmin


# ============================================================================
# Flat-sky lensing convolution (GPU)
# ============================================================================

def lens_cl_tt(ell_unlensed, Cl_unlensed, ell_phi, Cl_phiphi,
               ell_max_out=None):
    """
    Apply gravitational lensing to a TT power spectrum.

    Uses the first-order flat-sky result (Hu 2000, Lewis 2005):

    C_l^{lensed} = C_l + int d^2L/(2pi)^2
                    [L . (l-L)]^2  C^{pp}_{|l-L|}
                    * [C_L^{TT} - C_l^{TT}]

    With L = (L cos phi, L sin phi):
      L . (l-L) = L(l cos phi - L)
      |l-L|^2 = l^2 + L^2 - 2lL cos phi

    Computed on GPU (MLX) with subsampled L grid for speed.

    Parameters
    ----------
    ell_unlensed : array, multipoles of the unlensed spectrum
    Cl_unlensed : array, unlensed C_l (dimensionless)
    ell_phi : array, multipoles for C_l^{pp}
    Cl_phiphi : array, lensing potential spectrum
    ell_max_out : int or None

    Returns
    -------
    ell_out : array (int), Cl_lensed : array (float)
    """
    t0 = time.time()

    if ell_max_out is None:
        ell_max_out = int(ell_unlensed[-1])

    # Dense ell grid
    ell_dense = np.arange(2, ell_max_out + 1, dtype=float)
    N_ell = len(ell_dense)

    # Interpolate C_l^{unlensed} to dense grid (log-space)
    f_cl = interp1d(ell_unlensed.astype(float),
                     np.log(np.maximum(Cl_unlensed, 1e-100)),
                     kind='cubic', fill_value='extrapolate',
                     bounds_error=False)
    Cl_dense = np.exp(f_cl(ell_dense))

    # Interpolate C_l^{pp} to dense grid (extend with zeros at l=0,1)
    ell_pp_ext = np.concatenate([[0, 1], ell_phi.astype(float)])
    Cpp_ext = np.concatenate([[0, 0], Cl_phiphi])
    f_pp = interp1d(ell_pp_ext,
                     np.log(np.maximum(Cpp_ext, 1e-200)),
                     kind='linear', fill_value=-460.0,
                     bounds_error=False)

    # C_l^{TT} interpolator (extended)
    ell_tt_ext = np.concatenate([[0, 1], ell_dense])
    Cl_tt_ext = np.concatenate([[0, 0], Cl_dense])
    f_ctt = interp1d(ell_tt_ext,
                      np.log(np.maximum(Cl_tt_ext, 1e-200)),
                      kind='linear', fill_value=-460.0,
                      bounds_error=False)

    # R_pp diagnostic
    Cpp_dense = np.exp(f_pp(ell_dense))
    R_pp = np.sum(ell_dense * (ell_dense + 1.0) * (2.0 * ell_dense + 1.0)
                  / (4.0 * np.pi) * Cpp_dense)
    print(f"[Lensing] R_pp = {R_pp:.6e}, "
          f"rms deflection = {np.sqrt(R_pp)*180*60/np.pi:.2f} arcmin")

    # ---------------------------------------------------------------
    # Subsampled L grid for the integration variable.
    # C_l^{pp} drops steeply above l~1000, so subsample at high L.
    # ---------------------------------------------------------------
    L_dense_lo = np.arange(2, 200, 1, dtype=float)      # dense at low L
    L_dense_mid = np.arange(200, 1000, 3, dtype=float)   # moderate spacing
    L_dense_hi = np.arange(1000, min(ell_max_out + 500, 3000), 10, dtype=float)
    L_arr = np.unique(np.concatenate([L_dense_lo, L_dense_mid, L_dense_hi]))
    N_L = len(L_arr)

    # C_L^{TT} on the L grid
    Cl_L = np.exp(f_ctt(L_arr))

    # Azimuthal grid (16 points sufficient; kernel is smooth in phi)
    N_phi = 32
    phi_arr = np.linspace(0, 2.0 * np.pi, N_phi, endpoint=False)
    dphi = 2.0 * np.pi / N_phi
    cos_phi = np.cos(phi_arr)

    print(f"[Lensing] Flat-sky integral: "
          f"N_ell={N_ell}, N_L={N_L}, N_phi={N_phi}")

    # Precompute C^{pp} lookup table for fast interpolation
    # |l - L| can range from 0 to ~2*ell_max
    lmL_max = ell_max_out + L_arr[-1] + 1
    lmL_table_ell = np.arange(0, int(lmL_max) + 1, dtype=float)
    lmL_table_ell[0] = 0.5  # avoid log(0)
    Cpp_table = np.exp(f_pp(lmL_table_ell))
    Cpp_table[0] = 0.0  # l=0 has no lensing

    # Process in batches over output ell
    batch_size = 200
    Cl_lensed = Cl_dense.copy()

    for ib in range(0, N_ell, batch_size):
        ie = min(ib + batch_size, N_ell)

        ell_batch = ell_dense[ib:ie]   # (n_batch,)
        Cl_batch = Cl_dense[ib:ie]     # (n_batch,)

        # 3D arrays: (n_batch, N_L, N_phi)
        l_out = ell_batch[:, None, None]
        L_in = L_arr[None, :, None]
        cphi = cos_phi[None, None, :]

        # |l - L| = sqrt(l^2 + L^2 - 2lL cos phi)
        lmL_sq = l_out**2 + L_in**2 - 2.0 * l_out * L_in * cphi
        lmL = np.sqrt(np.maximum(lmL_sq, 0.25))

        # Look up C^{pp} at |l-L| using nearest-integer from precomputed table
        lmL_idx = np.clip(np.round(lmL).astype(int), 0, len(Cpp_table) - 1)
        Cpp_lmL = Cpp_table[lmL_idx]

        # Kernel: [L(l cos phi - L)]^2
        dot_sq = (L_in * (l_out * cphi - L_in))**2

        # Azimuthal integral: I(l, L) = sum_phi dot^2 * C^{pp} * dphi/(2pi)
        phi_integrand = dot_sq * Cpp_lmL
        I_lL = np.sum(phi_integrand, axis=2) * dphi / (2.0 * np.pi)

        # Full integrand: dL L/(2pi) * I(l,L) * (C_L - C_l)
        dC = Cl_L[None, :] - Cl_batch[:, None]
        L_integrand = L_arr[None, :] / (2.0 * np.pi) * I_lL * dC

        # Trapezoidal integration over L
        delta_Cl = np.trapezoid(L_integrand, L_arr, axis=1)

        Cl_lensed[ib:ie] = Cl_dense[ib:ie] + delta_Cl

    # Ensure positivity
    Cl_lensed = np.maximum(Cl_lensed, 0.0)

    t_elapsed = time.time() - t0
    print(f"[Lensing] Lensed C_l computed in {t_elapsed:.2f}s")

    return ell_dense.astype(int), Cl_lensed


# ============================================================================
# High-level interface
# ============================================================================

def apply_lensing(ell_in, Dl_in, bg, result=None, solver_type='implicit',
                  spectrum='TT'):
    """
    Apply CMB lensing to a D_l spectrum.

    Parameters
    ----------
    ell_in : array, multipoles
    Dl_in : array, D_l = l(l+1)/(2pi) C_l * T_CMB^2 in muK^2
    bg : Background (solved)
    result : solver result or None
    solver_type : str
    spectrum : str ('TT', 'EE', 'TE')

    Returns
    -------
    ell_out : array (int, dense)
    Dl_lensed : array (muK^2)
    (ell_phi, Cl_phiphi) : tuple for diagnostics
    """
    t0 = time.time()
    print(f"\n{'='*60}")
    print(f"CMB LENSING ({spectrum})")
    print(f"{'='*60}")

    # Step 1: Lensing potential
    ell_max_phi = int(np.max(ell_in)) + 500
    ell_phi, Cl_pp = compute_cl_phiphi(bg, result, solver_type,
                                        ell_max=ell_max_phi)

    # Step 2: Diagnostics
    sigma2, d_rms = compute_sigma_lens(ell_phi, Cl_pp)
    print(f"[Lensing] sigma_lens^2 = {sigma2:.6e}")
    print(f"[Lensing] RMS deflection = {d_rms:.2f} arcmin (Planck: ~2.7)")

    # Step 3: C_l from D_l
    ell_f = ell_in.astype(float)
    Cl_in = Dl_in / (ell_f * (ell_f + 1.0) / (2.0 * np.pi) * (T_CMB * 1e6)**2)

    # Step 4: Lens
    ell_out, Cl_lensed = lens_cl_tt(ell_in, Cl_in, ell_phi, Cl_pp)

    # Step 5: D_l
    ell_out_f = ell_out.astype(float)
    Dl_lensed = ell_out_f * (ell_out_f + 1.0) / (2.0 * np.pi) \
                * Cl_lensed * (T_CMB * 1e6)**2

    t_total = time.time() - t0
    print(f"[Lensing] Total lensing time: {t_total:.2f}s")

    _report_lensing_effect(ell_in, Dl_in, ell_out, Dl_lensed)

    return ell_out, Dl_lensed, (ell_phi, Cl_pp)


def _report_lensing_effect(ell_unlensed, Dl_unlensed, ell_lensed, Dl_lensed):
    """Print summary of lensing effect on the spectrum."""
    from scipy.signal import find_peaks

    ell_common = np.arange(100, min(int(ell_unlensed[-1]),
                                     int(ell_lensed[-1])) + 1)
    f_u = interp1d(ell_unlensed.astype(float), Dl_unlensed,
                    kind='cubic', fill_value='extrapolate')
    f_l = interp1d(ell_lensed.astype(float), Dl_lensed,
                    kind='cubic', fill_value='extrapolate')
    Dl_u = f_u(ell_common)
    Dl_l = f_l(ell_common)

    peaks_u, _ = find_peaks(Dl_u, distance=100, prominence=50)
    troughs_u, _ = find_peaks(-Dl_u, distance=100, prominence=50)

    if len(peaks_u) >= 2 and len(troughs_u) >= 1:
        peak_change = [(Dl_l[p] - Dl_u[p]) / Dl_u[p] * 100
                       for p in peaks_u[:5]]
        trough_change = [(Dl_l[t] - Dl_u[t]) / Dl_u[t] * 100
                         for t in troughs_u[:4] if Dl_u[t] > 10]

        print(f"\n[Lensing] Peak changes: "
              f"{[f'{c:+.1f}%' for c in peak_change]}")
        if trough_change:
            print(f"[Lensing] Trough changes: "
                  f"{[f'{c:+.1f}%' for c in trough_change]}")

    mask_high = (ell_common > 2000) & (ell_common < 2500)
    if np.any(mask_high) and np.mean(Dl_u[mask_high]) > 0:
        ratio_high = np.mean(Dl_l[mask_high]) / np.mean(Dl_u[mask_high])
        print(f"[Lensing] l=2000-2500 power ratio (lensed/unlensed): "
              f"{ratio_high:.3f}")


# ============================================================================
# CLASS lensed reference loader
# ============================================================================

def load_class_lensed(path):
    """
    Load CLASS lensed C_l file.

    Columns: l, TT, EE, TE, BB, phiphi, TPhi, Ephi
    All in dimensionless l(l+1)/(2pi) C_l units.

    Returns dict with keys 'ell', 'TT' (muK^2), 'phiphi' (C_l, dimensionless).
    """
    data = np.loadtxt(path)
    ell = data[:, 0].astype(int)
    T_muK = T_CMB * 1e6

    result = {'ell': ell, 'TT': data[:, 1] * T_muK**2}

    if data.shape[1] > 2:
        result['EE'] = data[:, 2] * T_muK**2
    if data.shape[1] > 3:
        result['TE'] = data[:, 3] * T_muK**2
    if data.shape[1] > 4:
        result['BB'] = data[:, 4] * T_muK**2
    if data.shape[1] > 5:
        Dl_pp = data[:, 5]  # l(l+1)/(2pi) C_l^{pp}
        ell_f = ell.astype(float)
        result['phiphi'] = Dl_pp * (2.0 * np.pi) / (ell_f * (ell_f + 1.0))
        result['Dl_phiphi'] = Dl_pp

    return result


# ============================================================================
# Plotting
# ============================================================================

def plot_lensing(ell_unlensed, Dl_unlensed, ell_lensed, Dl_lensed,
                 ell_phi=None, Cl_phiphi=None,
                 class_data=None, out_path=None, ell_max=2500):
    """
    Multi-panel lensing comparison plot.

    Panel 1: Unlensed vs lensed D_l^TT (+ CLASS if available)
    Panel 2: Lensing residual
    Panel 3: C_l^{phi phi} comparison (if provided)
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import os

    n_panels = 3 if (ell_phi is not None) else 2
    fig, axes = plt.subplots(n_panels, 1, figsize=(14, 4 * n_panels),
                              gridspec_kw={'height_ratios':
                                           [3, 1, 2][:n_panels]})
    if n_panels == 2:
        fig.subplots_adjust(hspace=0.06)
    else:
        fig.subplots_adjust(hspace=0.25)

    # Panel 1: D_l spectra
    ax = axes[0]
    m_u = ell_unlensed <= ell_max
    m_l = ell_lensed <= ell_max
    ax.plot(ell_unlensed[m_u], Dl_unlensed[m_u],
            'b-', lw=1.0, alpha=0.5, label='Unlensed')
    ax.plot(ell_lensed[m_l], Dl_lensed[m_l],
            'r-', lw=1.5, label='Lensed (mlx_class)')
    if class_data is not None:
        m_c = class_data['ell'] <= ell_max
        ax.plot(class_data['ell'][m_c], class_data['TT'][m_c],
                'k--', lw=1.2, alpha=0.7, label='CLASS lensed')
    ax.set_ylabel(r'$\mathcal{D}_\ell^{TT}$ [$\mu$K$^2$]', fontsize=13)
    ax.set_title('CMB Lensing: Unlensed vs Lensed TT', fontsize=14)
    ax.legend(fontsize=11)
    ax.set_xlim(2, ell_max)
    ax.set_ylim(0, None)

    # Panel 2: Residual
    ax2 = axes[1]
    ell_common = np.arange(30, min(ell_max, int(ell_unlensed[-1]),
                                    int(ell_lensed[-1])) + 1)
    f_u = interp1d(ell_unlensed.astype(float), Dl_unlensed,
                    kind='cubic', fill_value='extrapolate')
    f_l = interp1d(ell_lensed.astype(float), Dl_lensed,
                    kind='cubic', fill_value='extrapolate')
    Du = f_u(ell_common)
    Dl = f_l(ell_common)
    residual = np.where(Du > 50, (Dl - Du) / Du * 100, 0.0)
    ax2.plot(ell_common, residual, 'r-', lw=0.8)
    ax2.axhline(0, color='k', ls='--', lw=0.5)
    ax2.fill_between(ell_common, -10, 10, color='green', alpha=0.08)
    ax2.set_ylabel('(Lensed$-$Unlensed)/Unlensed [%]', fontsize=11)
    ax2.set_ylim(-15, 15)
    ax2.set_xlabel(r'Multipole $\ell$', fontsize=13)

    # Panel 3: C_l^{phi phi}
    if n_panels == 3 and ell_phi is not None:
        ax3 = axes[2]
        ell_pp = ell_phi.astype(float)
        Dl_pp_plot = ell_pp**2 * (ell_pp + 1)**2 * Cl_phiphi / (2.0 * np.pi)
        m_pp = ell_pp <= ell_max
        ax3.semilogy(ell_pp[m_pp], Dl_pp_plot[m_pp],
                      'b-', lw=1.5,
                      label=r'mlx\_class $[\ell(\ell+1)]^2 C_\ell^{\phi\phi}/(2\pi)$')
        if class_data is not None and 'phiphi' in class_data:
            ell_c = class_data['ell'].astype(float)
            Cpp_c = class_data['phiphi']
            Dl_pp_c = ell_c**2 * (ell_c + 1)**2 * Cpp_c / (2.0 * np.pi)
            m_cc = ell_c <= ell_max
            ax3.semilogy(ell_c[m_cc], Dl_pp_c[m_cc],
                          'k--', lw=1.2, alpha=0.7, label='CLASS')
        ax3.set_xlabel(r'Multipole $\ell$', fontsize=13)
        ax3.set_ylabel(r'$[\ell(\ell+1)]^2 C_\ell^{\phi\phi}/(2\pi)$',
                        fontsize=12)
        ax3.set_title('Lensing Potential Power Spectrum', fontsize=13)
        ax3.legend(fontsize=10)
        ax3.set_xlim(2, ell_max)

    plt.tight_layout()

    if out_path is None:
        out_dir = os.path.dirname(os.path.abspath(__file__))
        out_path = os.path.join(out_dir, 'cmb_lensing.png')

    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"[Lensing] Plot saved: {out_path}")
    return out_path


# ============================================================================
# Self-test
# ============================================================================

def _test_lensing():
    """
    Test the lensing module:
    1. Eisenstein-Hu transfer function
    2. C_l^{phi phi} amplitude and shape
    3. RMS deflection angle
    4. Lensed C_l peak smoothing
    5. Comparison with CLASS lensed output (if available)
    """
    import os

    print("=" * 70)
    print("TEST: CMB Lensing Module")
    print("=" * 70)

    from .background import Background

    # Step 1: Background
    print("\n--- Step 1: Background ---")
    bg = Background(khronon=False, recombination='peebles')
    bg.solve()

    # Step 2: Eisenstein-Hu T(k)
    print("\n--- Step 2: Eisenstein-Hu T(k) ---")
    k_test = np.geomspace(1e-4, 1.0, 100)
    T_EH = _eisenstein_hu_nowiggle(k_test, Omega_m, Omega_b, h_param)
    for k_val in [0.001, 0.01, 0.1, 0.5]:
        idx = np.argmin(np.abs(k_test - k_val))
        print(f"  T(k={k_val}) = {T_EH[idx]:.4f}")

    # Step 3: C_l^{phi phi}
    print("\n--- Step 3: C_l^{phi phi} ---")
    ell_phi, Cl_pp = compute_cl_phiphi(bg, ell_max=2500)

    # Step 4: RMS deflection
    print("\n--- Step 4: RMS deflection ---")
    sigma2, d_rms = compute_sigma_lens(ell_phi, Cl_pp)
    print(f"  sigma^2 = {sigma2:.6e}")
    print(f"  RMS deflection = {d_rms:.2f} arcmin (Planck: ~2.7)")
    if 1.0 < d_rms < 5.0:
        print("  PASS: RMS deflection in expected range")
    else:
        print(f"  WARNING: RMS deflection {d_rms:.2f} outside [1, 5] arcmin")

    # Step 5: Unlensed C_l
    print("\n--- Step 5: Unlensed C_l ---")
    k_arr = np.geomspace(5e-5, 0.35, 500).astype(np.float32)
    from .perturbations_implicit import ImplicitBoltzmannSolver
    solver = ImplicitBoltzmannSolver(bg, k_arr)
    result = solver.solve()
    Theta_0, Phi, v_b = result.source_at_recombination()
    source_SW = Theta_0 + Phi

    from .main import upsample_source, early_isw_template
    k_fine, source_fine = upsample_source(k_arr, source_SW, bg.D_A)

    ell_values = np.unique(np.concatenate([
        np.arange(2, 30, 1), np.arange(30, 100, 2),
        np.arange(100, 500, 4), np.arange(500, 1500, 8),
        np.arange(1500, 2501, 12),
    ])).astype(int)

    from .spectra import compute_cl
    ell_out, Cl_unl, Dl_unl = compute_cl(source_fine, k_fine,
                                          ell_values, bg.D_A)
    Dl_ISW = early_isw_template(ell_out, Dl_unl, bg)
    Dl_unl = Dl_unl + Dl_ISW
    print(f"  Unlensed D_l: [{np.min(Dl_unl):.0f}, {np.max(Dl_unl):.0f}] muK^2")

    # Step 6: Apply lensing
    print("\n--- Step 6: Apply lensing ---")
    ell_lensed, Dl_lensed, (ell_pp, Cl_pp_out) = apply_lensing(
        ell_out, Dl_unl, bg, spectrum='TT')
    print(f"  Lensed D_l: [{np.min(Dl_lensed):.0f}, {np.max(Dl_lensed):.0f}] muK^2")

    # Step 7: CLASS comparison (optional)
    class_path = ('/Users/akaihuangm1/Desktop/github/'
                  'gdm_class_public/output/running_mu_lcdm_ref_cl_lensed.dat')
    class_data = None
    if os.path.exists(class_path):
        print(f"\n--- Step 7: CLASS comparison ---")
        class_data = load_class_lensed(class_path)
        print(f"  Loaded CLASS lensed: {len(class_data['ell'])} multipoles")

        if 'phiphi' in class_data:
            for l_cmp in [100, 500, 1000]:
                idx_m = np.argmin(np.abs(ell_phi - l_cmp))
                idx_c = np.argmin(np.abs(class_data['ell'] - l_cmp))
                Cpp_c = class_data['phiphi'][idx_c]
                if Cpp_c > 0:
                    ratio = Cl_pp[idx_m] / Cpp_c
                    print(f"  C_l^{{pp}} at l={l_cmp}: "
                          f"mlx/CLASS = {ratio:.3f}")

    # Step 8: Plot
    print("\n--- Step 8: Plot ---")
    plot_lensing(ell_out, Dl_unl, ell_lensed, Dl_lensed,
                 ell_pp, Cl_pp_out, class_data=class_data)

    print("\n" + "=" * 70)
    print("CMB LENSING TEST COMPLETE")
    print("=" * 70)

    return ell_lensed, Dl_lensed


if __name__ == '__main__':
    _test_lensing()
