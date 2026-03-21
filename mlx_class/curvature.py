"""
curvature.py -- Non-flat spatial geometry (Omega_k != 0) for mlx_class.

Extends Background to support open (Omega_k > 0) and closed (Omega_k < 0)
universes. For |Omega_k| < 0.01 (current Planck + BAO constraints), the
effect on perturbation evolution is negligible; only background distances
are modified.

Friedmann equation with curvature:
    H^2(a) = H0^2 [Omega_r/a^4 + Omega_m/a^3 + Omega_k/a^2 + Omega_L]

where:
    Omega_k = 1 - Omega_r - Omega_m - Omega_L

Comoving angular diameter distance:
    chi(z) = integral_0^z c dz' / H(z')           (coordinate distance)
    D_M(z) = S_k(chi)                              (transverse comoving)

    S_k(chi) = (c/H0) / sqrt(|Omega_k|) * sinh(sqrt(Omega_k)  * H0*chi/c)   [open,   Omega_k > 0]
             = chi                                                             [flat,   Omega_k = 0]
             = (c/H0) / sqrt(|Omega_k|) * sin(sqrt(|Omega_k|) * H0*chi/c)    [closed, Omega_k < 0]

Angular diameter distance:
    D_A(z) = D_M(z) / (1+z)

Luminosity distance:
    D_L(z) = D_M(z) * (1+z) = D_A(z) * (1+z)^2

Comoving volume element:
    dV/dz/dOmega = D_M(z)^2 * c / H(z)

Author: Sheng-Kai Huang, 2026
"""
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import interp1d

from .background import (
    Background,
    h, H0_km_s_Mpc, c_km_s, H0_Mpc,
    omega_b, omega_c,
    Omega_b, Omega_c, Omega_r, Omega_m, Omega_L,
    T_CMB, Y_He,
    a_eq, k_eq,
    A_s, n_s, k_pivot,
    z_rec, a_rec,
    sigma_T_SI, c_SI, Mpc_SI, m_H_SI, G_SI, H0_SI,
)


# ============================================================================
# Curvature transverse comoving distance
# ============================================================================

def transverse_comoving_distance(chi, Omega_k, H0=H0_Mpc, c=1.0):
    """
    Convert coordinate (radial) comoving distance chi to transverse
    comoving distance D_M = S_k(chi), accounting for spatial curvature.

    Parameters
    ----------
    chi : array_like
        Radial comoving distance (in Mpc, same units as c/H0).
    Omega_k : float
        Curvature density parameter. >0 open, <0 closed, 0 flat.
    H0 : float
        Hubble constant in same units as 1/chi (default: 1/Mpc).
    c : float
        Speed of light in distance_units * H0_units (default: 1 for
        conformal-time Mpc where chi is already in Mpc).

    Returns
    -------
    D_M : ndarray
        Transverse comoving distance in Mpc.
    """
    chi = np.asarray(chi, dtype=np.float64)

    if abs(Omega_k) < 1e-12:
        # Flat
        return chi.copy()

    # sqrt(|Omega_k|) * H0 / c  has units 1/Mpc
    sqrtK = np.sqrt(abs(Omega_k)) * H0 / c
    R_curv = c / (H0 * np.sqrt(abs(Omega_k)))   # curvature radius in Mpc

    if Omega_k > 0:
        # Open: S_k = R_curv * sinh(chi / R_curv)
        return R_curv * np.sinh(chi / R_curv)
    else:
        # Closed: S_k = R_curv * sin(chi / R_curv)
        return R_curv * np.sin(chi / R_curv)


# ============================================================================
# BackgroundCurved class
# ============================================================================

class BackgroundCurved(Background):
    """
    Background cosmology with non-flat spatial geometry.

    Extends Background to include curvature Omega_k in the Friedmann equation
    and to compute the correct transverse comoving distances.

    For |Omega_k| < 0.01, perturbation evolution is unaffected; only the
    geometric projection (D_A, D_L, volume element) changes.

    Parameters
    ----------
    Omega_k : float
        Spatial curvature parameter. Default 0.0 (flat).
        Positive = open (hyperbolic), negative = closed (spherical).
        Omega_L is derived as 1 - Omega_m - Omega_r - Omega_k.
    **kwargs
        Passed to Background.__init__ (N_points, a_min, a_max, khronon,
        recombination).

    Examples
    --------
    # Flat (should match Background exactly)
    >>> bg = BackgroundCurved(Omega_k=0.0).solve()

    # Slightly open universe
    >>> bg = BackgroundCurved(Omega_k=+0.01).solve()

    # Slightly closed universe
    >>> bg = BackgroundCurved(Omega_k=-0.01).solve()
    """

    def __init__(self, Omega_k=0.0, **kwargs):
        # Initialize parent (sets Omega_m, Omega_L, grids, etc.)
        super().__init__(**kwargs)

        self.Omega_k = Omega_k

        # Re-derive Omega_L to satisfy the constraint:
        #   Omega_r + Omega_m + Omega_k + Omega_L = 1
        self.Omega_L = 1.0 - self.Omega_m - Omega_r - self.Omega_k

        if self.Omega_L < 0:
            raise ValueError(
                f"Omega_L = {self.Omega_L:.6f} < 0. The combination "
                f"Omega_m={self.Omega_m:.4f}, Omega_r={Omega_r:.6f}, "
                f"Omega_k={self.Omega_k:.4f} is unphysical."
            )

    def solve(self):
        """
        Compute background quantities with spatial curvature.

        Overrides Background.solve() to use the modified Friedmann equation:
            E^2(a) = Omega_r/a^4 + Omega_m/a^3 + Omega_k/a^2 + Omega_L

        The conformal time, sound horizon, and recombination are computed
        as in the flat case. After that, the angular diameter distance to
        the last scattering surface is corrected for curvature.
        """
        a = np.logspace(np.log10(self.a_min), np.log10(self.a_max), self.N)

        # --- Modified Friedmann equation with curvature ---
        E = np.sqrt(
            Omega_r / a**4
            + self.Omega_m / a**3
            + self.Omega_k / a**2
            + self.Omega_L
        )

        # Conformal time: dtau = da / (a^2 H) = da / (a^2 E H0)
        dtau = np.diff(a) / (a[:-1]**2 * E[:-1] * H0_Mpc)
        tau = np.concatenate([[0.0], np.cumsum(dtau)])

        self.a_grid = a
        self.E_grid = E
        self.tau_grid = tau
        self.dtau = dtau

        # Recombination index
        self.idx_rec = np.searchsorted(a, a_rec)
        self.tau_rec = tau[self.idx_rec]
        self.tau_0 = tau[-1]

        # Radial comoving distance to last scattering (conformal time difference)
        chi_rec = self.tau_0 - self.tau_rec

        # Apply curvature correction to get the transverse comoving distance
        # In conformal coordinates, chi is in Mpc, and H0_Mpc = H0 in 1/Mpc,
        # c = 1 in these units (since dtau = dt/a and chi = integral c dtau).
        self.chi_rec = chi_rec
        self.D_M_rec = transverse_comoving_distance(
            chi_rec, self.Omega_k, H0=H0_Mpc, c=1.0
        )

        # D_A is the transverse comoving distance (for CMB peak positions,
        # what matters is D_A_comoving = D_M; the angular diameter distance
        # proper is D_M / (1+z_rec), but for l_1 = pi * D_M / r_s we use D_M)
        self.D_A = self.D_M_rec

        # Baryon loading R = 3 rho_b / (4 rho_gamma)
        Omega_gamma = Omega_r / (1.0 + 0.2271 * 3.046)
        self.R_grid = 3.0 * Omega_b * a / (4.0 * Omega_gamma)
        self.R_rec = self.R_grid[self.idx_rec]

        # Sound speed and sound horizon
        cs = 1.0 / np.sqrt(3.0 * (1.0 + self.R_grid))
        self.cs_grid = cs
        self.r_s = np.sum(cs[:self.idx_rec - 1] * dtau[:self.idx_rec - 1])

        # Conformal Hubble: calH = a * H = a * E * H0_Mpc
        self.calH_grid = a * E * H0_Mpc

        # Silk damping scale (approximate)
        self.k_D = 0.15 * (omega_b / 0.022)**0.25

        # Interpolators
        self._build_interpolators()

        # Recombination
        if self.recombination_method == 'peebles':
            self._compute_recombination_peebles()
        else:
            self._compute_recombination()

        return self

    # === Distance measures with curvature ===

    def H_of_a(self, a):
        """
        Hubble parameter H(a) in km/s/Mpc.

        Parameters
        ----------
        a : array_like
            Scale factor.

        Returns
        -------
        H : ndarray
            H(a) in km/s/Mpc.
        """
        a = np.asarray(a, dtype=np.float64)
        E = np.sqrt(
            Omega_r / a**4
            + self.Omega_m / a**3
            + self.Omega_k / a**2
            + self.Omega_L
        )
        return E * H0_km_s_Mpc

    def comoving_distance(self, z):
        """
        Radial comoving distance chi(z) in Mpc.

        chi(z) = integral_0^z c dz' / H(z')

        Parameters
        ----------
        z : array_like
            Redshift(s).

        Returns
        -------
        chi : ndarray
            Radial comoving distance in Mpc.
        """
        z = np.atleast_1d(np.asarray(z, dtype=np.float64))
        z_max = np.max(z)
        z_fine = np.linspace(0, z_max, max(10000, int(z_max * 1000)))
        a_fine = 1.0 / (1.0 + z_fine)
        H_fine = self.H_of_a(a_fine)  # km/s/Mpc
        integrand = c_km_s / H_fine  # Mpc
        chi_cum = cumulative_trapezoid(integrand, z_fine, initial=0.0)
        chi_interp = interp1d(z_fine, chi_cum, kind='cubic')
        return chi_interp(z)

    def transverse_distance(self, z):
        """
        Transverse comoving distance D_M(z) in Mpc.

        D_M = S_k(chi), where S_k accounts for spatial curvature.

        Parameters
        ----------
        z : array_like
            Redshift(s).

        Returns
        -------
        D_M : ndarray
            Transverse comoving distance in Mpc.
        """
        chi = self.comoving_distance(z)
        # Here chi is in Mpc and H0 in km/s/Mpc, so R_curv = c/H0 in Mpc
        return transverse_comoving_distance(
            chi, self.Omega_k, H0=H0_km_s_Mpc / c_km_s, c=1.0
        )

    def angular_diameter_distance(self, z):
        """
        Angular diameter distance D_A(z) in Mpc.

        D_A = D_M(z) / (1+z)

        Parameters
        ----------
        z : array_like
            Redshift(s).

        Returns
        -------
        D_A : ndarray
            Angular diameter distance in Mpc.
        """
        z = np.atleast_1d(np.asarray(z, dtype=np.float64))
        D_M = self.transverse_distance(z)
        return D_M / (1.0 + z)

    def luminosity_distance(self, z):
        """
        Luminosity distance D_L(z) in Mpc.

        D_L = D_M(z) * (1+z)

        Parameters
        ----------
        z : array_like
            Redshift(s).

        Returns
        -------
        D_L : ndarray
            Luminosity distance in Mpc.
        """
        z = np.atleast_1d(np.asarray(z, dtype=np.float64))
        D_M = self.transverse_distance(z)
        return D_M * (1.0 + z)

    def comoving_volume_element(self, z):
        """
        Comoving volume element dV/dz/dOmega in Mpc^3/sr.

        dV/dz/dOmega = D_M(z)^2 * c / H(z)

        Parameters
        ----------
        z : array_like
            Redshift(s).

        Returns
        -------
        dVdzdOmega : ndarray
            Volume element in Mpc^3/sr.
        """
        z = np.atleast_1d(np.asarray(z, dtype=np.float64))
        a = 1.0 / (1.0 + z)
        D_M = self.transverse_distance(z)
        H = self.H_of_a(a)  # km/s/Mpc
        return D_M**2 * c_km_s / H

    def age_of_universe(self):
        """
        Age of the universe t_0 in Gyr.

        t_0 = integral_0^1 da / (a * H(a))
        """
        a = np.logspace(-10, 0, 100000)
        E = np.sqrt(
            Omega_r / a**4
            + self.Omega_m / a**3
            + self.Omega_k / a**2
            + self.Omega_L
        )
        H = E * H0_SI  # 1/s
        integrand = 1.0 / (a * H)  # seconds
        try:
            t0_s = np.trapezoid(integrand, a)
        except AttributeError:
            t0_s = np.trapz(integrand, a)
        t0_Gyr = t0_s / (3600 * 24 * 365.25 * 1e9)
        return t0_Gyr

    def deceleration_parameter(self, a):
        """
        Deceleration parameter q(a).

        q = 0.5 * [Omega_r(a) + Omega_m(a) - Omega_k(a) + (1+3w_L)*Omega_L(a)]

        For w_L = -1 (cosmological constant):
        q = 0.5 * [Omega_r(a) + Omega_m(a) - Omega_k(a) - 2*Omega_L(a)]
        """
        a = np.asarray(a, dtype=np.float64)
        E2 = (Omega_r / a**4 + self.Omega_m / a**3
              + self.Omega_k / a**2 + self.Omega_L)

        Omega_r_a = Omega_r / a**4 / E2
        Omega_m_a = self.Omega_m / a**3 / E2
        Omega_k_a = self.Omega_k / a**2 / E2
        Omega_L_a = self.Omega_L / E2

        # q = -a*ddot{a}/dot{a}^2
        # = 0.5 * (Omega_r(a) + Omega_m(a)) - Omega_L(a)
        # curvature contributes 0 to the deceleration (it's in the Friedmann
        # equation but doesn't contribute to pressure)
        return 0.5 * (Omega_r_a + Omega_m_a) - Omega_L_a


# ============================================================================
# Self-tests
# ============================================================================

def _test_flat_consistency():
    """Verify BackgroundCurved(Omega_k=0) matches Background exactly."""
    print("=" * 65)
    print("TEST 1: Flat consistency (BackgroundCurved(Omega_k=0) vs Background)")
    print("=" * 65)

    bg_std = Background(khronon=False).solve()
    bg_curv = BackgroundCurved(Omega_k=0.0, khronon=False).solve()

    checks = [
        ('tau_0', bg_std.tau_0, bg_curv.tau_0),
        ('tau_rec', bg_std.tau_rec, bg_curv.tau_rec),
        ('D_A', bg_std.D_A, bg_curv.D_A),
        ('r_s', bg_std.r_s, bg_curv.r_s),
        ('R_rec', bg_std.R_rec, bg_curv.R_rec),
        ('Omega_L', bg_std.Omega_L, bg_curv.Omega_L),
    ]

    all_pass = True
    for name, val_std, val_curv in checks:
        rel_err = abs(val_std - val_curv) / abs(val_std) if val_std != 0 else abs(val_curv)
        status = "PASS" if rel_err < 1e-10 else "FAIL"
        if status == "FAIL":
            all_pass = False
        print(f"  {name:>10s}: std={val_std:.6f}  curv={val_curv:.6f}"
              f"  err={rel_err:.2e}  [{status}]")

    # E(a) should match everywhere
    E_err = np.max(np.abs(bg_std.E_grid - bg_curv.E_grid) / bg_std.E_grid)
    status = "PASS" if E_err < 1e-10 else "FAIL"
    if status == "FAIL":
        all_pass = False
    print(f"  {'E(a) max':>10s}: err={E_err:.2e}  [{status}]")

    print(f"\n  Result: {'ALL PASS' if all_pass else 'SOME FAILED'}\n")
    return all_pass


def _test_open_universe():
    """Test open universe (Omega_k > 0)."""
    print("=" * 65)
    print("TEST 2: Open universe (Omega_k = +0.01)")
    print("=" * 65)

    bg_flat = BackgroundCurved(Omega_k=0.0).solve()
    bg_open = BackgroundCurved(Omega_k=+0.01).solve()

    print(f"\n  Omega_k = {bg_open.Omega_k:+.4f}")
    print(f"  Omega_L = {bg_open.Omega_L:.6f}  (flat: {bg_flat.Omega_L:.6f})")
    print(f"  tau_0   = {bg_open.tau_0:.1f} Mpc  (flat: {bg_flat.tau_0:.1f})")
    print(f"  chi_rec = {bg_open.chi_rec:.1f} Mpc  (flat D_A: {bg_flat.D_A:.1f})")
    print(f"  D_M_rec = {bg_open.D_M_rec:.1f} Mpc  (D_A: {bg_open.D_A:.1f})")
    print(f"  r_s     = {bg_open.r_s:.2f} Mpc  (flat: {bg_flat.r_s:.2f})")

    # Open universe: sinh > identity => D_M > chi
    assert bg_open.D_M_rec > bg_open.chi_rec, \
        "Open: D_M should be > chi (sinh > identity)"
    print(f"  D_M > chi (sinh): PASS")

    # For small positive Omega_k, with Omega_L reduced, the universe
    # expands slightly differently. D_M is typically larger than flat D_A.
    l1_flat = np.pi * bg_flat.D_A / bg_flat.r_s
    l1_open = np.pi * bg_open.D_A / bg_open.r_s
    print(f"\n  l_1 (flat) = {l1_flat:.1f}")
    print(f"  l_1 (open) = {l1_open:.1f}")
    print(f"  l_1 shift:   {(l1_open - l1_flat)/l1_flat*100:+.2f}%")

    # r_s should be nearly the same (curvature negligible at recombination)
    r_s_err = abs(bg_open.r_s - bg_flat.r_s) / bg_flat.r_s
    print(f"\n  r_s change: {r_s_err*100:.4f}%")
    assert r_s_err < 0.01, f"r_s changed too much: {r_s_err*100:.2f}%"
    print(f"  r_s stability < 1%: PASS")

    # Age
    age_flat = bg_flat.age_of_universe()
    age_open = bg_open.age_of_universe()
    print(f"\n  Age (flat) = {age_flat:.2f} Gyr")
    print(f"  Age (open) = {age_open:.2f} Gyr")

    print()
    return True


def _test_closed_universe():
    """Test closed universe (Omega_k < 0)."""
    print("=" * 65)
    print("TEST 3: Closed universe (Omega_k = -0.01)")
    print("=" * 65)

    bg_flat = BackgroundCurved(Omega_k=0.0).solve()
    bg_closed = BackgroundCurved(Omega_k=-0.01).solve()

    print(f"\n  Omega_k  = {bg_closed.Omega_k:+.4f}")
    print(f"  Omega_L  = {bg_closed.Omega_L:.6f}  (flat: {bg_flat.Omega_L:.6f})")
    print(f"  tau_0    = {bg_closed.tau_0:.1f} Mpc  (flat: {bg_flat.tau_0:.1f})")
    print(f"  chi_rec  = {bg_closed.chi_rec:.1f} Mpc")
    print(f"  D_M_rec  = {bg_closed.D_M_rec:.1f} Mpc  (D_A: {bg_closed.D_A:.1f})")
    print(f"  r_s      = {bg_closed.r_s:.2f} Mpc  (flat: {bg_flat.r_s:.2f})")

    # Closed universe: sin < identity => D_M < chi
    assert bg_closed.D_M_rec < bg_closed.chi_rec, \
        "Closed: D_M should be < chi (sin < identity)"
    print(f"  D_M < chi (sin): PASS")

    l1_flat = np.pi * bg_flat.D_A / bg_flat.r_s
    l1_closed = np.pi * bg_closed.D_A / bg_closed.r_s
    print(f"\n  l_1 (flat)   = {l1_flat:.1f}")
    print(f"  l_1 (closed) = {l1_closed:.1f}")
    print(f"  l_1 shift:     {(l1_closed - l1_flat)/l1_flat*100:+.2f}%")

    # r_s should be nearly the same
    r_s_err = abs(bg_closed.r_s - bg_flat.r_s) / bg_flat.r_s
    print(f"\n  r_s change: {r_s_err*100:.4f}%")
    assert r_s_err < 0.01
    print(f"  r_s stability < 1%: PASS")

    # Age
    age_flat = bg_flat.age_of_universe()
    age_closed = bg_closed.age_of_universe()
    print(f"\n  Age (flat)   = {age_flat:.2f} Gyr")
    print(f"  Age (closed) = {age_closed:.2f} Gyr")

    print()
    return True


def _test_distances():
    """Test distance measures for curved universes."""
    print("=" * 65)
    print("TEST 4: Distance measures (Etherington reciprocity)")
    print("=" * 65)

    all_pass = True

    for Ok_val, label in [(0.0, 'flat'), (+0.01, 'open'), (-0.01, 'closed')]:
        bg = BackgroundCurved(Omega_k=Ok_val).solve()

        z_test = np.array([0.5, 1.0, 2.0, 5.0])
        chi = bg.comoving_distance(z_test)
        D_M = bg.transverse_distance(z_test)
        D_A = bg.angular_diameter_distance(z_test)
        D_L = bg.luminosity_distance(z_test)

        # Etherington reciprocity: D_L = D_A * (1+z)^2
        D_L_from_DA = D_A * (1.0 + z_test)**2
        err = np.max(np.abs(D_L - D_L_from_DA) / D_L)
        status = "PASS" if err < 1e-10 else "FAIL"
        if status == "FAIL":
            all_pass = False
        print(f"  [{label:>6s}] Etherington D_L = D_A*(1+z)^2: err={err:.2e}  [{status}]")

        # D_A = D_M / (1+z)
        D_A_from_DM = D_M / (1.0 + z_test)
        err2 = np.max(np.abs(D_A - D_A_from_DM) / D_A)
        status2 = "PASS" if err2 < 1e-10 else "FAIL"
        if status2 == "FAIL":
            all_pass = False
        print(f"  [{label:>6s}] D_A = D_M / (1+z):             err={err2:.2e}  [{status2}]")

        # For flat case, D_M = chi
        if abs(Ok_val) < 1e-12:
            err3 = np.max(np.abs(D_M - chi) / chi)
            status3 = "PASS" if err3 < 1e-10 else "FAIL"
            if status3 == "FAIL":
                all_pass = False
            print(f"  [{label:>6s}] D_M = chi (flat):              err={err3:.2e}  [{status3}]")

    print(f"\n  Result: {'ALL PASS' if all_pass else 'SOME FAILED'}\n")
    return all_pass


def _test_symmetry():
    """Test that Omega_k = +epsilon and -epsilon give opposite D_A shifts."""
    print("=" * 65)
    print("TEST 5: Symmetry (open/closed give opposite D_A shifts)")
    print("=" * 65)

    bg_flat = BackgroundCurved(Omega_k=0.0).solve()
    bg_open = BackgroundCurved(Omega_k=+0.005).solve()
    bg_closed = BackgroundCurved(Omega_k=-0.005).solve()

    D_A_flat = bg_flat.D_A
    delta_open = bg_open.D_A - D_A_flat
    delta_closed = bg_closed.D_A - D_A_flat

    print(f"  D_A (flat)   = {D_A_flat:.2f} Mpc")
    print(f"  D_A (open)   = {bg_open.D_A:.2f} Mpc  (delta = {delta_open:+.2f})")
    print(f"  D_A (closed) = {bg_closed.D_A:.2f} Mpc  (delta = {delta_closed:+.2f})")

    # The shifts should be in opposite directions
    # Open (Omega_k > 0) => less Omega_L => different expansion history
    # Plus sinh > identity, so D_M is boosted
    # Closed (Omega_k < 0) => more Omega_L => different expansion history
    # Plus sin < identity, so D_M is reduced
    assert delta_open * delta_closed < 0, \
        "Open and closed should shift D_A in opposite directions"
    print(f"  Opposite sign shifts: PASS")

    # To leading order in Omega_k, shifts should be approximately equal
    # in magnitude (not exact due to Omega_L adjustment)
    ratio = abs(delta_open / delta_closed)
    print(f"  |delta_open / delta_closed| = {ratio:.3f}  (expect ~1)")
    assert 0.5 < ratio < 2.0, f"Asymmetry too large: ratio = {ratio:.3f}"
    print(f"  Approximate symmetry: PASS")

    print()
    return True


def _test_volume_element():
    """Test comoving volume element."""
    print("=" * 65)
    print("TEST 6: Comoving volume element")
    print("=" * 65)

    bg = BackgroundCurved(Omega_k=+0.01).solve()

    z_test = np.array([0.5, 1.0, 2.0])
    dV = bg.comoving_volume_element(z_test)

    print(f"  {'z':>5s}  {'dV/dz/dOmega (Mpc^3/sr)':>28s}")
    print(f"  " + "-" * 36)
    for z, v in zip(z_test, dV):
        print(f"  {z:5.1f}  {v:28.2f}")

    # Volume element should be positive and finite
    assert np.all(dV > 0), "Volume element must be positive"
    assert np.all(np.isfinite(dV)), "Volume element must be finite"
    print(f"\n  Positive and finite: PASS")

    print()
    return True


def run_all_tests():
    """Run all curvature self-tests."""
    print("\n" + "=" * 65)
    print("CURVATURE MODULE SELF-TESTS")
    print("=" * 65 + "\n")

    results = {}
    results['Flat consistency'] = _test_flat_consistency()
    results['Open universe'] = _test_open_universe()
    results['Closed universe'] = _test_closed_universe()
    results['Distance measures'] = _test_distances()
    results['Symmetry'] = _test_symmetry()
    results['Volume element'] = _test_volume_element()

    print("\n" + "=" * 65)
    print("SUMMARY")
    print("=" * 65)
    for name, passed in results.items():
        print(f"  {name:30s}: {'PASS' if passed else 'FAIL'}")

    n_pass = sum(results.values())
    n_total = len(results)
    print(f"\n  {n_pass}/{n_total} tests passed.")
    print("=" * 65)

    return all(results.values())


# ============================================================================
# CLI entry point
# ============================================================================

if __name__ == '__main__':
    import sys
    success = run_all_tests()
    sys.exit(0 if success else 1)
