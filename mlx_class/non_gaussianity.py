"""
non_gaussianity.py — Primordial non-Gaussianity (f_NL) support for mlx_class.

Local-type non-Gaussianity:
    Phi(x) = phi(x) + f_NL * (phi^2(x) - <phi^2>)

where phi is the Gaussian primordial potential. The primary observable effect
is a scale-dependent correction to the galaxy bias at large scales:

    Delta_b(k,z) = (b_1 - 1) * f_NL * delta_c * 3 * Omega_m * H0^2
                   / (k^2 * T(k) * D(z) * c^2)

This is the Dalal et al. (2008) formula. At k ~ 0.001 h/Mpc, even modest
f_NL ~ 5 produces ~10% corrections to clustering, making it a powerful
probe of primordial physics.

For the CMB, the leading effect is the primordial bispectrum:
    B(l1,l2,l3) ~ f_NL * (C_l1 * C_l2 + C_l1 * C_l3 + C_l2 * C_l3)

References:
    Dalal et al. (2008), PRD 77, 123514 [arXiv:0710.4560]
    Matarrese & Verde (2008), ApJL 677, L77 [arXiv:0801.4826]
    Komatsu & Spergel (2001), PRD 63, 063002 [arXiv:astro-ph/0005036]

Input:  Background object (for Omega_m, H0, growth factor, transfer function)
Output: scale-dependent bias b(k,z), CMB bispectrum template B(l1,l2,l3)

Author: Sheng-Kai Huang, 2026
"""
import numpy as np
from scipy.integrate import quad

from .background import (
    Omega_m, Omega_r, Omega_L,
    H0_Mpc, h as h_param, A_s, n_s, k_pivot,
)
from .matter_pk import growth_factor_integral


# ============================================================================
# Constants
# ============================================================================

# Critical overdensity for spherical collapse (EdS approximation)
DELTA_C = 1.686

# Speed of light in Mpc * H0 units: c = 1 / H0_Mpc (since H0_Mpc is in 1/Mpc)
# Actually c in km/s / (H0 in km/s/Mpc) gives c/H0 in Mpc.
# More directly: c / H0 = c_km_s / H0_km_s_Mpc  [Mpc]
# Here we work in natural units where k is in Mpc^-1, so the Poisson
# equation factor is just 3/2 * Omega_m * H0^2 / k^2.
# Since we use conformal-Hubble units (H0_Mpc in 1/Mpc), no c factor needed.


# ============================================================================
# BBKS transfer function (Bardeen, Bond, Kaiser, Szalay 1986)
# ============================================================================

def bbks_transfer_function(k_Mpc, Omega_m_val=None, h_val=None, Omega_b_val=None):
    """
    BBKS transfer function T(k) for CDM, with Sugiyama (1995) shape parameter.

    T(q) = ln(1 + 2.34 q) / (2.34 q)
           * [1 + 3.89 q + (16.1 q)^2 + (5.46 q)^3 + (6.71 q)^4]^{-1/4}

    where q = k / (Gamma * h) and Gamma is the shape parameter:
        Gamma = Omega_m * h * exp(-Omega_b - sqrt(2h) * Omega_b / Omega_m)

    Parameters
    ----------
    k_Mpc : array
        Wavenumber in Mpc^-1.
    Omega_m_val : float or None
        Total matter density parameter. Default: module-level Omega_m.
    h_val : float or None
        Hubble parameter. Default: module-level h.
    Omega_b_val : float or None
        Baryon density parameter. Default: uses Omega_b from background.

    Returns
    -------
    T_k : array
        Transfer function T(k), normalized to 1 at k -> 0.
    """
    from .background import Omega_b as Omega_b_default

    if Omega_m_val is None:
        Omega_m_val = Omega_m
    if h_val is None:
        h_val = h_param
    if Omega_b_val is None:
        Omega_b_val = Omega_b_default

    # Sugiyama shape parameter
    Gamma = Omega_m_val * h_val * np.exp(
        -Omega_b_val * (1.0 + np.sqrt(2.0 * h_val) / Omega_m_val)
    )

    # q = k / (Gamma * h) in h/Mpc units, but k is in Mpc^-1
    # k [Mpc^-1] = k_h [h/Mpc] * h  =>  k_h = k / h
    # q = k_h / Gamma = k / (h * Gamma)
    q = k_Mpc / (h_val * Gamma + 1e-30)

    # BBKS formula
    T_k = np.log(1.0 + 2.34 * q) / (2.34 * q + 1e-30)
    T_k *= (1.0 + 3.89 * q + (16.1 * q)**2
            + (5.46 * q)**3 + (6.71 * q)**4)**(-0.25)

    # Fix k -> 0 limit
    T_k = np.where(q < 1e-10, 1.0, T_k)

    return T_k


# ============================================================================
# Linear growth factor D(z) (normalized to D(0) = 1)
# ============================================================================

def growth_factor_at_z(z, Omega_m_val=None, Omega_L_val=None):
    """
    Linear growth factor D(z), normalized so that D(z=0) = 1.

    Uses the Heath (1977) integral from matter_pk.py.

    Parameters
    ----------
    z : float or array
        Redshift(s).
    Omega_m_val : float or None
        Matter density parameter.
    Omega_L_val : float or None
        Dark energy density parameter.

    Returns
    -------
    D_z : float or array
        Growth factor D(z) with D(0) = 1.
    """
    if Omega_m_val is None:
        Omega_m_val = Omega_m
    if Omega_L_val is None:
        Omega_L_val = Omega_L

    D_0 = growth_factor_integral(1.0, Omega_m_val, Omega_L_val)

    scalar_input = np.ndim(z) == 0
    z_arr = np.atleast_1d(np.asarray(z, dtype=np.float64))
    a_arr = 1.0 / (1.0 + z_arr)

    D_z = np.array([
        growth_factor_integral(a, Omega_m_val, Omega_L_val)
        for a in a_arr
    ])

    result = D_z / D_0
    if scalar_input:
        return float(result[0])
    return result


# ============================================================================
# Scale-dependent bias: Dalal et al. (2008)
# ============================================================================

def scale_dependent_bias(k_Mpc, z, f_NL, b0, bg=None,
                         Omega_m_val=None, Omega_L_val=None,
                         T_k=None):
    """
    Galaxy bias with local-type non-Gaussianity correction.

    b(k, z) = b_0 + Delta_b(k, z)

    where
        Delta_b(k, z) = (b_0 - 1) * delta_c * f_NL
                         * 3 * Omega_m * H0^2 / (k^2 * T(k) * D(z))

    This is the Dalal et al. (2008) formula for scale-dependent bias
    induced by local f_NL. The 1/k^2 scaling makes the effect dominant
    at large scales (small k).

    Parameters
    ----------
    k_Mpc : array
        Wavenumber in Mpc^-1.
    z : float
        Redshift.
    f_NL : float
        Local non-Gaussianity parameter. Planck (2018): f_NL = -0.9 +/- 5.1.
    b0 : float
        Gaussian (scale-independent) linear bias at this redshift.
    bg : Background instance or None
        If provided, uses bg.Omega_cdm for Khronon cosmology.
        If None, uses module-level Planck parameters.
    Omega_m_val : float or None
        Override matter density parameter.
    Omega_L_val : float or None
        Override dark energy density parameter.
    T_k : array or None
        Pre-computed transfer function T(k). If None, uses BBKS approximation.

    Returns
    -------
    b_k : array
        Scale-dependent bias b(k, z). Shape matches k_Mpc.
    """
    k_Mpc = np.atleast_1d(np.asarray(k_Mpc, dtype=np.float64))

    # Cosmological parameters
    if Omega_m_val is None:
        if bg is not None:
            from .background import Omega_b
            Omega_m_val = Omega_b + bg.Omega_cdm
        else:
            Omega_m_val = Omega_m
    if Omega_L_val is None:
        if bg is not None:
            Omega_L_val = 1.0 - Omega_m_val - Omega_r
        else:
            Omega_L_val = Omega_L

    # f_NL = 0 => constant bias
    if f_NL == 0.0:
        return np.full_like(k_Mpc, b0)

    # Transfer function
    if T_k is None:
        T_k = bbks_transfer_function(k_Mpc, Omega_m_val=Omega_m_val)
    else:
        T_k = np.atleast_1d(np.asarray(T_k, dtype=np.float64))

    # Growth factor D(z) normalized to D(0) = 1
    D_z = growth_factor_at_z(z, Omega_m_val=Omega_m_val, Omega_L_val=Omega_L_val)

    # Dalal et al. (2008) scale-dependent bias correction
    # Delta_b = (b0 - 1) * delta_c * f_NL * 3 * Omega_m * H0^2 / (k^2 * T(k) * D(z))
    #
    # Units: H0_Mpc is in Mpc^-1, k is in Mpc^-1, so H0^2 / k^2 is dimensionless.
    # T(k) and D(z) are dimensionless. Result is dimensionless correction to bias.
    numerator = (b0 - 1.0) * DELTA_C * f_NL * 3.0 * Omega_m_val * H0_Mpc**2
    denominator = k_Mpc**2 * T_k * D_z

    # Guard against division by zero at k -> 0 (unphysical)
    safe_denom = np.where(np.abs(denominator) < 1e-30, 1e-30, denominator)
    Delta_b = numerator / safe_denom

    b_k = b0 + Delta_b

    return b_k


# ============================================================================
# Galaxy angular power spectrum with f_NL
# ============================================================================

def galaxy_cl_with_fnl(ell_values, z_mean, f_NL, b0, Pk_interp,
                       dNdz=None, bg=None,
                       Omega_m_val=None, Omega_L_val=None,
                       k_min=1e-4, k_max=1.0, N_k=500):
    """
    Galaxy angular power spectrum C_l^gg with f_NL scale-dependent bias.

    Uses the Limber approximation:
        C_l^gg = int dz (dN/dz)^2 * H(z) / (chi(z))^2
                 * b(k=l/chi, z)^2 * P(k=l/chi, z)

    Parameters
    ----------
    ell_values : array
        Multipole values.
    z_mean : float
        Mean redshift of the galaxy sample.
    f_NL : float
        Local non-Gaussianity parameter.
    b0 : float
        Gaussian linear bias.
    Pk_interp : callable
        P(k) interpolator, P(k_Mpc) in Mpc^3 at z=0.
    dNdz : callable or None
        Normalized redshift distribution dN/dz. If None, uses a delta
        function at z_mean (thin-shell approximation).
    bg : Background or None
        Background cosmology object.
    Omega_m_val : float or None
        Matter density override.
    Omega_L_val : float or None
        Dark energy density override.
    k_min, k_max : float
        k range in Mpc^-1.
    N_k : int
        Number of k points.

    Returns
    -------
    ell_values : array
        Multipole values.
    Cl_gg : array
        Galaxy angular power spectrum.
    """
    ell_values = np.atleast_1d(np.asarray(ell_values, dtype=np.float64))

    if Omega_m_val is None:
        if bg is not None:
            from .background import Omega_b
            Omega_m_val = Omega_b + bg.Omega_cdm
        else:
            Omega_m_val = Omega_m
    if Omega_L_val is None:
        Omega_L_val = 1.0 - Omega_m_val - Omega_r

    # Comoving distance chi(z) using Hubble parameter
    # chi(z) = c * int_0^z dz' / H(z')
    # In Mpc units with H0_Mpc: chi = int_0^z dz' / E(z') / H0_Mpc
    def _E(z_val):
        a = 1.0 / (1.0 + z_val)
        return np.sqrt(Omega_r / a**4 + Omega_m_val / a**3 + Omega_L_val)

    # Thin-shell (Limber) at z_mean
    chi_mean, _ = quad(lambda zp: 1.0 / (_E(zp) * H0_Mpc), 0.0, z_mean,
                       limit=200)

    # Growth factor
    D_z = growth_factor_at_z(z_mean, Omega_m_val, Omega_L_val)

    # For each ell, k = (ell + 0.5) / chi (Limber)
    k_ell = (ell_values + 0.5) / chi_mean  # Mpc^-1

    # Transfer function for bias correction
    T_k_ell = bbks_transfer_function(k_ell, Omega_m_val=Omega_m_val)

    # Scale-dependent bias
    b_k = scale_dependent_bias(
        k_ell, z_mean, f_NL, b0, bg=bg,
        Omega_m_val=Omega_m_val, Omega_L_val=Omega_L_val,
        T_k=T_k_ell
    )

    # P(k, z) = D(z)^2 * P(k, z=0)
    Pk_ell = np.array([Pk_interp(k) for k in k_ell]) * D_z**2

    # Limber C_l (thin shell):
    # C_l = b(k)^2 * P(k) * H(z_mean) / chi(z_mean)^2
    H_z_mean = _E(z_mean) * H0_Mpc  # in Mpc^-1
    Cl_gg = b_k**2 * Pk_ell * H_z_mean / chi_mean**2

    return ell_values, Cl_gg


# ============================================================================
# CMB bispectrum template for local f_NL
# ============================================================================

def compute_fnl_bispectrum(Cl, ell_values, f_NL):
    """
    Compute the reduced CMB bispectrum template for local-type non-Gaussianity.

    For the local ansatz Phi = phi + f_NL (phi^2 - <phi^2>), the angle-averaged
    reduced bispectrum is:

        b_{l1 l2 l3} = 2 * f_NL * (C_{l1} * C_{l2} + C_{l1} * C_{l3}
                                    + C_{l2} * C_{l3})

    This returns the reduced bispectrum for all valid triangles
    (l1 <= l2 <= l3, l1+l2 >= l3) from the input ell_values.

    Parameters
    ----------
    Cl : array (N_ell,)
        CMB angular power spectrum C_l (NOT D_l). Must correspond to ell_values.
    ell_values : array (N_ell,)
        Multipole values.
    f_NL : float
        Local non-Gaussianity parameter.

    Returns
    -------
    triangles : list of tuples
        List of (l1, l2, l3) triangles.
    B_l1l2l3 : array
        Reduced bispectrum for each triangle.
    """
    Cl = np.atleast_1d(np.asarray(Cl, dtype=np.float64))
    ell_values = np.atleast_1d(np.asarray(ell_values, dtype=np.int64))

    # Build ell -> Cl lookup
    ell_min = int(ell_values[0])
    ell_max = int(ell_values[-1])
    Cl_dict = {}
    for i, ell in enumerate(ell_values):
        Cl_dict[int(ell)] = Cl[i]

    # Enumerate valid triangles (l1 <= l2 <= l3, triangle inequality)
    triangles = []
    B_values = []

    for i, l1 in enumerate(ell_values):
        l1 = int(l1)
        for j in range(i, len(ell_values)):
            l2 = int(ell_values[j])
            # l3 range from triangle inequality: |l1-l2| <= l3 <= l1+l2
            l3_min = max(l2, abs(l1 - l2))  # also l3 >= l2
            l3_max = min(l1 + l2, ell_max)
            for l3 in range(l3_min, l3_max + 1):
                if l3 not in Cl_dict:
                    continue
                # Parity: l1 + l2 + l3 must be even for scalar CMB
                if (l1 + l2 + l3) % 2 != 0:
                    continue

                C1 = Cl_dict[l1]
                C2 = Cl_dict[l2]
                C3 = Cl_dict[l3]

                b_123 = 2.0 * f_NL * (C1 * C2 + C1 * C3 + C2 * C3)
                triangles.append((l1, l2, l3))
                B_values.append(b_123)

    return triangles, np.array(B_values, dtype=np.float64)


# ============================================================================
# Cumulative S/N for f_NL from bispectrum
# ============================================================================

def bispectrum_snr(Cl, ell_values, f_NL, f_sky=1.0):
    """
    Estimate the signal-to-noise ratio (S/N) for detecting f_NL from
    the CMB temperature bispectrum.

    (S/N)^2 = sum_{l1<=l2<=l3} B_{l1l2l3}^2 / Var(B)

    where Var(B) ~ C_{l1} * C_{l2} * C_{l3} * (symmetry factor).

    Parameters
    ----------
    Cl : array
        CMB power spectrum C_l.
    ell_values : array
        Multipole values.
    f_NL : float
        Non-Gaussianity amplitude.
    f_sky : float
        Sky fraction (default 1.0 for full sky).

    Returns
    -------
    snr : float
        Total signal-to-noise ratio.
    """
    triangles, B_vals = compute_fnl_bispectrum(Cl, ell_values, f_NL)

    if len(B_vals) == 0:
        return 0.0

    Cl = np.atleast_1d(np.asarray(Cl, dtype=np.float64))
    ell_values = np.atleast_1d(np.asarray(ell_values, dtype=np.int64))
    Cl_dict = {int(ell_values[i]): Cl[i] for i in range(len(ell_values))}

    snr_sq = 0.0
    for idx, (l1, l2, l3) in enumerate(triangles):
        C1 = Cl_dict[l1]
        C2 = Cl_dict[l2]
        C3 = Cl_dict[l3]

        # Symmetry factor
        if l1 == l2 == l3:
            Delta_123 = 6
        elif l1 == l2 or l2 == l3:
            Delta_123 = 2
        else:
            Delta_123 = 1

        # Variance of the bispectrum estimator (Gaussian cosmic variance)
        var_B = C1 * C2 * C3 * Delta_123

        if var_B > 0:
            snr_sq += B_vals[idx]**2 / var_B

    snr_sq *= f_sky
    return np.sqrt(max(snr_sq, 0.0))


# ============================================================================
# Self-test
# ============================================================================

def _test_non_gaussianity():
    """
    Test non-Gaussianity module.

    1. f_NL = 0 should give constant bias.
    2. f_NL != 0 should give k-dependent bias, diverging at small k.
    3. Bispectrum should scale linearly with f_NL.
    4. Bispectrum should vanish when f_NL = 0.
    """
    print("=" * 70)
    print("TEST: Primordial Non-Gaussianity (f_NL)")
    print("=" * 70)

    # k grid
    k_arr = np.geomspace(1e-4, 1.0, 500)  # Mpc^-1

    # --- Test 1: f_NL = 0 gives constant bias ---
    print("\n--- Test 1: f_NL = 0 => constant bias ---")
    b0 = 1.5
    b_k_zero = scale_dependent_bias(k_arr, z=1.0, f_NL=0.0, b0=b0)
    assert np.allclose(b_k_zero, b0), "FAIL: f_NL=0 should give constant bias"
    print(f"  b(k) = {b_k_zero[0]:.6f} ... {b_k_zero[-1]:.6f} (expect {b0:.1f})")
    print("  PASS: constant bias for f_NL = 0")

    # --- Test 2: f_NL = 10, scale-dependent bias ---
    print("\n--- Test 2: f_NL = 10, scale-dependent bias ---")
    f_NL_test = 10.0
    b_k_fnl = scale_dependent_bias(k_arr, z=1.0, f_NL=f_NL_test, b0=b0)
    print(f"  b(k=1e-4) = {b_k_fnl[0]:.4f}")
    print(f"  b(k=1e-2) = {b_k_fnl[np.argmin(np.abs(k_arr - 0.01))]:.4f}")
    print(f"  b(k=1e-1) = {b_k_fnl[np.argmin(np.abs(k_arr - 0.1))]:.4f}")
    print(f"  b(k=1.0)  = {b_k_fnl[-1]:.4f}")

    # At large k, bias should approach b0
    assert np.abs(b_k_fnl[-1] - b0) < 0.01, "FAIL: bias should approach b0 at high k"
    # At small k, bias should deviate significantly
    assert np.abs(b_k_fnl[0] - b0) > 0.1, "FAIL: bias should deviate at low k"
    # Bias should decrease with k (for positive f_NL and b0 > 1)
    assert b_k_fnl[0] > b_k_fnl[-1], "FAIL: bias should be larger at low k for f_NL > 0"
    print("  PASS: scale-dependent bias with correct k-dependence")

    # --- Test 3: f_NL linearity ---
    print("\n--- Test 3: f_NL linearity ---")
    b_k_fnl5 = scale_dependent_bias(k_arr, z=1.0, f_NL=5.0, b0=b0)
    b_k_fnl20 = scale_dependent_bias(k_arr, z=1.0, f_NL=20.0, b0=b0)
    Delta_b_5 = b_k_fnl5 - b0
    Delta_b_20 = b_k_fnl20 - b0
    # Delta_b should scale linearly with f_NL
    ratio = Delta_b_20 / (Delta_b_5 + 1e-30)
    mask_valid = np.abs(Delta_b_5) > 1e-10
    ratio_valid = ratio[mask_valid]
    assert np.allclose(ratio_valid, 4.0, atol=0.01), (
        f"FAIL: Delta_b should scale as f_NL (expected ratio 4, got {ratio_valid[0]:.4f})"
    )
    print(f"  Delta_b(f_NL=20) / Delta_b(f_NL=5) = {ratio_valid[0]:.4f} (expect 4.0)")
    print("  PASS: linear scaling with f_NL")

    # --- Test 4: redshift dependence ---
    print("\n--- Test 4: redshift dependence ---")
    b_z0 = scale_dependent_bias(k_arr, z=0.0, f_NL=f_NL_test, b0=b0)
    b_z1 = scale_dependent_bias(k_arr, z=1.0, f_NL=f_NL_test, b0=b0)
    b_z2 = scale_dependent_bias(k_arr, z=2.0, f_NL=f_NL_test, b0=b0)
    # At higher z, D(z) is smaller => Delta_b is larger
    Delta_z0 = np.abs(b_z0[0] - b0)
    Delta_z1 = np.abs(b_z1[0] - b0)
    Delta_z2 = np.abs(b_z2[0] - b0)
    print(f"  |Delta_b(z=0, k=1e-4)| = {Delta_z0:.4f}")
    print(f"  |Delta_b(z=1, k=1e-4)| = {Delta_z1:.4f}")
    print(f"  |Delta_b(z=2, k=1e-4)| = {Delta_z2:.4f}")
    assert Delta_z2 > Delta_z1 > Delta_z0, "FAIL: Delta_b should grow with z"
    print("  PASS: Delta_b increases with redshift (as 1/D(z))")

    # --- Test 5: CMB bispectrum ---
    print("\n--- Test 5: CMB bispectrum ---")
    ell_test = np.array([2, 3, 4, 5, 6, 7, 8, 9, 10])
    # Mock C_l ~ 1/ell^2 (rough CMB shape)
    Cl_test = 1e-10 / ell_test**2

    # f_NL = 0
    tri_0, B_0 = compute_fnl_bispectrum(Cl_test, ell_test, f_NL=0.0)
    assert np.allclose(B_0, 0.0), "FAIL: B = 0 for f_NL = 0"
    print(f"  f_NL = 0: {len(tri_0)} triangles, all B = 0")
    print("  PASS: bispectrum vanishes for f_NL = 0")

    # f_NL = 10
    tri_10, B_10 = compute_fnl_bispectrum(Cl_test, ell_test, f_NL=10.0)
    print(f"  f_NL = 10: {len(tri_10)} triangles")
    if len(B_10) > 0:
        print(f"  B range: [{np.min(B_10):.4e}, {np.max(B_10):.4e}]")
        # Check linearity: f_NL = 20 should give 2x the bispectrum
        _, B_20 = compute_fnl_bispectrum(Cl_test, ell_test, f_NL=20.0)
        ratio_B = B_20 / (B_10 + 1e-50)
        assert np.allclose(ratio_B[np.abs(B_10) > 1e-30], 2.0, atol=1e-10), \
            "FAIL: bispectrum should scale linearly with f_NL"
        print("  PASS: bispectrum scales linearly with f_NL")

    # --- Test 6: BBKS transfer function ---
    print("\n--- Test 6: BBKS transfer function ---")
    T_k = bbks_transfer_function(k_arr)
    print(f"  T(k=1e-4) = {T_k[0]:.6f} (expect ~1.0)")
    print(f"  T(k=0.1)  = {T_k[np.argmin(np.abs(k_arr - 0.1))]:.6f}")
    print(f"  T(k=1.0)  = {T_k[-1]:.6f}")
    assert np.abs(T_k[0] - 1.0) < 0.01, "FAIL: T(k) -> 1 at small k"
    assert T_k[-1] < T_k[0], "FAIL: T(k) should decrease at high k"
    print("  PASS: BBKS transfer function has correct shape")

    # --- Test 7: S/N estimate ---
    print("\n--- Test 7: Bispectrum S/N ---")
    snr = bispectrum_snr(Cl_test, ell_test, f_NL=10.0)
    print(f"  S/N(f_NL=10, ell=2-10) = {snr:.4f}")
    snr_zero = bispectrum_snr(Cl_test, ell_test, f_NL=0.0)
    assert snr_zero == 0.0, "FAIL: S/N should be 0 for f_NL = 0"
    assert snr > 0.0, "FAIL: S/N should be > 0 for f_NL != 0"
    print("  PASS: S/N is positive for f_NL != 0 and zero for f_NL = 0")

    # --- Test 8: No NaN/Inf ---
    print("\n--- Test 8: Numerical stability ---")
    has_nan = np.any(np.isnan(b_k_fnl)) or np.any(np.isnan(B_10))
    has_inf = np.any(np.isinf(b_k_fnl)) or np.any(np.isinf(B_10))
    print(f"  NaN: {has_nan}, Inf: {has_inf}")
    assert not has_nan, "FAIL: NaN detected"
    assert not has_inf, "FAIL: Inf detected"
    print("  PASS: all values finite")

    print("\n" + "=" * 70)
    print("NON-GAUSSIANITY TEST COMPLETE — ALL PASSED")
    print("=" * 70)


if __name__ == '__main__':
    _test_non_gaussianity()
