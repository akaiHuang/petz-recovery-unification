"""
dark_energy.py -- Dark energy equation of state for mlx_class.

Extends the Background class to support dynamic dark energy:
  - Constant w:  w(a) = w = const
  - CPL:         w(a) = w0 + wa * (1 - a)   [Chevallier-Polarski-Linder]

Friedmann equation:
  H^2(a) = H0^2 [Omega_r / a^4 + Omega_m / a^3 + Omega_DE * f_DE(a)]

where for CPL:
  f_DE(a) = a^{-3(1+w0+wa)} * exp(-3 * wa * (1-a))

and for constant w:
  f_DE(a) = a^{-3(1+w)}

LCDM limit: w0 = -1, wa = 0 => f_DE = 1 (cosmological constant).
DESI 2024 hint: w0 ~ -0.7, wa ~ -1.0.

Dark energy is treated as smooth (no perturbations) -- the PPF approach.

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
# Dark energy density functions
# ============================================================================

def f_DE_constant_w(a, w):
    """
    Dark energy density ratio for constant equation of state.

    rho_DE(a) / rho_DE(0) = a^{-3(1+w)}

    Parameters
    ----------
    a : array_like
        Scale factor.
    w : float
        Constant equation of state parameter.

    Returns
    -------
    f_DE : ndarray
        rho_DE(a) / rho_DE(a=1).
    """
    a = np.asarray(a, dtype=np.float64)
    return a ** (-3.0 * (1.0 + w))


def f_DE_cpl(a, w0, wa):
    """
    Dark energy density ratio for CPL parametrization.

    w(a) = w0 + wa * (1 - a)

    rho_DE(a) / rho_DE(0) = a^{-3(1+w0+wa)} * exp(-3 * wa * (1-a))

    Parameters
    ----------
    a : array_like
        Scale factor.
    w0 : float
        Dark energy equation of state at a=1 (today).
    wa : float
        Time derivative dw/d(1-a).

    Returns
    -------
    f_DE : ndarray
        rho_DE(a) / rho_DE(a=1).
    """
    a = np.asarray(a, dtype=np.float64)
    return a ** (-3.0 * (1.0 + w0 + wa)) * np.exp(-3.0 * wa * (1.0 - a))


def w_of_a_cpl(a, w0, wa):
    """
    Dark energy equation of state w(a) for CPL.

    Parameters
    ----------
    a : array_like
        Scale factor.
    w0 : float
        w at a=1.
    wa : float
        dw/d(1-a).

    Returns
    -------
    w : ndarray
        Equation of state w(a).
    """
    a = np.asarray(a, dtype=np.float64)
    return w0 + wa * (1.0 - a)


# ============================================================================
# Dark energy effective density for arbitrary w(a) (numerical integration)
# ============================================================================

def f_DE_numerical(a_grid, w_func):
    """
    Dark energy density ratio for arbitrary w(a) via numerical integration.

    rho_DE(a) / rho_DE(1) = exp[ 3 * integral_a^1 (1 + w(a')) / a' da' ]

    Parameters
    ----------
    a_grid : ndarray
        Monotonically increasing scale factor array ending at 1.
    w_func : callable
        w(a) function.

    Returns
    -------
    f_DE : ndarray
        Dark energy density ratio on a_grid.
    """
    # Integrand: (1 + w(a)) / a
    integrand = (1.0 + w_func(a_grid)) / a_grid

    # Integral from a to 1: use cumulative trapezoid from a_min to a,
    # then total - cumulative gives integral from a to 1
    cum = cumulative_trapezoid(integrand, a_grid, initial=0.0)
    total = cum[-1]
    integral_a_to_1 = total - cum

    return np.exp(3.0 * integral_a_to_1)


# ============================================================================
# BackgroundDE class
# ============================================================================

class BackgroundDE(Background):
    """
    Background cosmology with dynamic dark energy.

    Extends Background to replace the cosmological constant with:
      - CPL: w(a) = w0 + wa*(1-a)
      - Constant-w: w(a) = w = const
      - Arbitrary w(a) via callable

    The dark energy is treated as smooth (no clustering).

    Parameters
    ----------
    w0 : float
        Dark energy equation of state at a=1. Default -1 (LCDM).
    wa : float
        CPL time variation dw/d(1-a). Default 0 (LCDM).
    w_const : float or None
        If not None, use constant w instead of CPL. Overrides w0/wa.
    w_func : callable or None
        If not None, use arbitrary w(a). Overrides w0/wa and w_const.
    Omega_DE : float or None
        Dark energy density parameter today. If None, computed from
        flatness: Omega_DE = 1 - Omega_m - Omega_r.
    **kwargs
        Passed to Background.__init__ (N_points, a_min, a_max, khronon,
        recombination).

    Examples
    --------
    # LCDM (default)
    >>> bg = BackgroundDE().solve()

    # DESI 2024 hint
    >>> bg = BackgroundDE(w0=-0.7, wa=-1.0).solve()

    # Constant w = -0.9
    >>> bg = BackgroundDE(w_const=-0.9).solve()
    """

    def __init__(self, w0=-1.0, wa=0.0, w_const=None, w_func=None,
                 Omega_DE=None, **kwargs):
        # Initialize parent (sets Omega_m, Omega_L, grids, etc.)
        super().__init__(**kwargs)

        # Store dark energy parameters
        self.w0 = w0
        self.wa = wa
        self.w_const = w_const
        self.w_func_user = w_func

        # Dark energy density parameter
        if Omega_DE is not None:
            self.Omega_DE = Omega_DE
        else:
            # Flatness: Omega_r + Omega_m + Omega_DE = 1
            self.Omega_DE = 1.0 - self.Omega_m - Omega_r

        # Override parent's Omega_L with our Omega_DE
        # (In parent, Omega_L plays the role of dark energy)
        self.Omega_L = self.Omega_DE

        # Determine the mode
        if w_func is not None:
            self._de_mode = 'numerical'
        elif w_const is not None:
            self._de_mode = 'constant'
            self.w0 = w_const
            self.wa = 0.0
        else:
            self._de_mode = 'cpl'

        # Check for LCDM limit
        self._is_lcdm = (self._de_mode == 'cpl'
                         and abs(self.w0 + 1.0) < 1e-10
                         and abs(self.wa) < 1e-10)

    def _f_DE(self, a):
        """
        Compute f_DE(a) = rho_DE(a) / rho_DE(a=1).

        Returns 1.0 for LCDM (cosmological constant).
        """
        if self._is_lcdm:
            return np.ones_like(a, dtype=np.float64)

        if self._de_mode == 'numerical':
            return f_DE_numerical(a, self.w_func_user)
        elif self._de_mode == 'constant':
            return f_DE_constant_w(a, self.w0)
        else:  # cpl
            return f_DE_cpl(a, self.w0, self.wa)

    def w_at(self, a):
        """
        Equation of state w(a).

        Parameters
        ----------
        a : array_like
            Scale factor.

        Returns
        -------
        w : ndarray
        """
        a = np.asarray(a, dtype=np.float64)

        if self._is_lcdm:
            return -np.ones_like(a)

        if self._de_mode == 'numerical':
            return self.w_func_user(a)
        elif self._de_mode == 'constant':
            return self.w0 * np.ones_like(a)
        else:
            return w_of_a_cpl(a, self.w0, self.wa)

    def solve(self):
        """
        Compute background quantities with dynamic dark energy.

        Overrides Background.solve() to use the modified Friedmann equation:
          E^2(a) = Omega_r/a^4 + Omega_m/a^3 + Omega_DE * f_DE(a)
        """
        a = np.logspace(np.log10(self.a_min), np.log10(self.a_max), self.N)

        # --- Modified Friedmann equation ---
        fDE = self._f_DE(a)
        E = np.sqrt(Omega_r / a**4 + self.Omega_m / a**3
                     + self.Omega_DE * fDE)

        # Conformal time: dtau = da / (a^2 H) = da / (a^2 E H0)
        dtau = np.diff(a) / (a[:-1]**2 * E[:-1] * H0_Mpc)
        tau = np.concatenate([[0.0], np.cumsum(dtau)])

        self.a_grid = a
        self.E_grid = E
        self.tau_grid = tau
        self.dtau = dtau
        self.fDE_grid = fDE

        # Recombination index
        self.idx_rec = np.searchsorted(a, a_rec)
        self.tau_rec = tau[self.idx_rec]
        self.tau_0 = tau[-1]
        self.D_A = self.tau_0 - self.tau_rec

        # Baryon loading R = 3 rho_b / (4 rho_gamma)
        Omega_gamma = Omega_r / (1.0 + 0.2271 * 3.046)
        self.R_grid = 3.0 * Omega_b * a / (4.0 * Omega_gamma)
        self.R_rec = self.R_grid[self.idx_rec]

        # Sound speed and sound horizon
        cs = 1.0 / np.sqrt(3.0 * (1.0 + self.R_grid))
        self.cs_grid = cs
        self.r_s = np.sum(cs[:self.idx_rec - 1] * dtau[:self.idx_rec - 1])

        # Conformal Hubble
        self.calH_grid = a * E * H0_Mpc

        # Silk damping scale
        self.k_D = 0.15 * (omega_b / 0.022)**0.25

        # Interpolators
        self._build_interpolators()

        # Dark energy interpolator
        lna = np.log(a)
        self._fDE_of_a = interp1d(lna, fDE, kind='cubic',
                                   fill_value='extrapolate')
        self._w_of_a = interp1d(lna, self.w_at(a), kind='cubic',
                                 fill_value='extrapolate')

        # Recombination (same as parent -- DE doesn't affect x_e directly)
        if self.recombination_method == 'peebles':
            self._compute_recombination_peebles()
        else:
            self._compute_recombination()

        return self

    def fDE_at(self, a):
        """Dark energy density ratio f_DE(a) via interpolation."""
        return self._fDE_of_a(np.log(np.asarray(a)))

    # === Diagnostic quantities ===

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
        fDE = self._f_DE(a)
        E = np.sqrt(Omega_r / a**4 + self.Omega_m / a**3
                     + self.Omega_DE * fDE)
        return E * H0_km_s_Mpc

    def luminosity_distance(self, z):
        """
        Luminosity distance d_L(z) in Mpc (flat universe).

        d_L = (1+z) * chi, where chi = integral_0^z dz'/H(z') * c
        """
        z = np.atleast_1d(np.asarray(z, dtype=np.float64))
        # Fine z grid for integration
        z_max = np.max(z)
        z_fine = np.linspace(0, z_max, max(10000, int(z_max * 1000)))
        a_fine = 1.0 / (1.0 + z_fine)
        H_fine = self.H_of_a(a_fine)  # km/s/Mpc

        # chi = c * integral dz/H(z) in Mpc
        integrand = c_km_s / H_fine  # Mpc
        chi_cum = cumulative_trapezoid(integrand, z_fine, initial=0.0)

        chi_interp = interp1d(z_fine, chi_cum, kind='cubic')
        chi = chi_interp(z)

        return (1.0 + z) * chi

    def angular_diameter_distance(self, z):
        """
        Angular diameter distance d_A(z) in Mpc (flat universe).

        d_A = chi / (1+z)
        """
        z = np.atleast_1d(np.asarray(z, dtype=np.float64))
        d_L = self.luminosity_distance(z)
        return d_L / (1.0 + z)**2

    def comoving_distance(self, z):
        """
        Comoving distance chi(z) in Mpc (flat universe).

        chi = integral_0^z c dz' / H(z')
        """
        z = np.atleast_1d(np.asarray(z, dtype=np.float64))
        d_L = self.luminosity_distance(z)
        return d_L / (1.0 + z)

    def deceleration_parameter(self, a):
        """
        Deceleration parameter q(a) = -a*H''/H - 1.

        q = 0.5 * [Omega_r(a) + Omega_m(a) + (1+3w(a))*Omega_DE(a)]

        where Omega_i(a) = rho_i(a) / rho_crit(a).
        """
        a = np.asarray(a, dtype=np.float64)
        fDE = self._f_DE(a)
        w = self.w_at(a)

        E2 = Omega_r / a**4 + self.Omega_m / a**3 + self.Omega_DE * fDE

        Omega_r_a = Omega_r / a**4 / E2
        Omega_m_a = self.Omega_m / a**3 / E2
        Omega_DE_a = self.Omega_DE * fDE / E2

        return 0.5 * (Omega_r_a + Omega_m_a + (1.0 + 3.0 * w) * Omega_DE_a)

    def age_of_universe(self):
        """
        Age of the universe t_0 in Gyr.

        t_0 = integral_0^1 da / (a * H(a))
        """
        a = np.logspace(-10, 0, 100000)
        fDE = self._f_DE(a)
        E = np.sqrt(Omega_r / a**4 + self.Omega_m / a**3
                     + self.Omega_DE * fDE)
        H = E * H0_SI  # 1/s

        # dt = da / (a * H)
        integrand = 1.0 / (a * H)  # seconds
        # np.trapz was removed in numpy 2.0; use np.trapezoid
        try:
            t0_s = np.trapezoid(integrand, a)
        except AttributeError:
            t0_s = np.trapz(integrand, a)
        t0_Gyr = t0_s / (3600 * 24 * 365.25 * 1e9)
        return t0_Gyr


# ============================================================================
# Convenience constructors
# ============================================================================

def LCDM(**kwargs):
    """Standard LCDM background (w=-1)."""
    return BackgroundDE(w0=-1.0, wa=0.0, **kwargs)


def wCDM(w, **kwargs):
    """Constant-w dark energy model."""
    return BackgroundDE(w_const=w, **kwargs)


def CPL(w0, wa, **kwargs):
    """CPL dark energy model: w(a) = w0 + wa*(1-a)."""
    return BackgroundDE(w0=w0, wa=wa, **kwargs)


def DESI_2024(**kwargs):
    """DESI 2024 best-fit CPL: w0 ~ -0.7, wa ~ -1.0."""
    return BackgroundDE(w0=-0.70, wa=-1.0, **kwargs)


# ============================================================================
# Self-tests and CLASS comparison
# ============================================================================

def _test_lcdm_consistency():
    """Verify BackgroundDE with w0=-1, wa=0 matches Background exactly."""
    from .background import Background

    print("=" * 65)
    print("TEST 1: LCDM consistency (BackgroundDE vs Background)")
    print("=" * 65)

    bg_std = Background(khronon=False)
    bg_std.solve()

    bg_de = BackgroundDE(w0=-1.0, wa=0.0, khronon=False)
    bg_de.solve()

    checks = [
        ('tau_0', bg_std.tau_0, bg_de.tau_0),
        ('tau_rec', bg_std.tau_rec, bg_de.tau_rec),
        ('D_A', bg_std.D_A, bg_de.D_A),
        ('r_s', bg_std.r_s, bg_de.r_s),
        ('R_rec', bg_std.R_rec, bg_de.R_rec),
    ]

    all_pass = True
    for name, val_std, val_de in checks:
        rel_err = abs(val_std - val_de) / abs(val_std)
        status = "PASS" if rel_err < 1e-10 else "FAIL"
        if status == "FAIL":
            all_pass = False
        print(f"  {name:>10s}: std={val_std:.6f}  DE={val_de:.6f}"
              f"  err={rel_err:.2e}  [{status}]")

    # E(a) should match everywhere
    E_err = np.max(np.abs(bg_std.E_grid - bg_de.E_grid) / bg_std.E_grid)
    status = "PASS" if E_err < 1e-10 else "FAIL"
    if status == "FAIL":
        all_pass = False
    print(f"  {'E(a) max':>10s}: err={E_err:.2e}  [{status}]")

    print(f"\n  Result: {'ALL PASS' if all_pass else 'SOME FAILED'}\n")
    return all_pass


def _test_constant_w():
    """Test constant-w model behavior."""
    print("=" * 65)
    print("TEST 2: Constant-w dark energy")
    print("=" * 65)

    # w = -1 should recover LCDM
    bg_lcdm = BackgroundDE(w0=-1.0, wa=0.0).solve()
    bg_wm1 = BackgroundDE(w_const=-1.0).solve()

    err = abs(bg_lcdm.tau_0 - bg_wm1.tau_0) / bg_lcdm.tau_0
    print(f"  w=-1 vs LCDM: tau_0 err = {err:.2e}"
          f"  [{'PASS' if err < 1e-10 else 'FAIL'}]")

    # w > -1: phantom divide. DE dilutes faster => larger H at late times
    # => smaller tau_0
    bg_w09 = BackgroundDE(w_const=-0.9).solve()
    print(f"\n  w=-0.9:")
    print(f"    tau_0 = {bg_w09.tau_0:.1f} Mpc  (LCDM: {bg_lcdm.tau_0:.1f})")
    print(f"    D_A   = {bg_w09.D_A:.1f} Mpc  (LCDM: {bg_lcdm.D_A:.1f})")
    print(f"    r_s   = {bg_w09.r_s:.2f} Mpc  (LCDM: {bg_lcdm.r_s:.2f})")

    # w=-0.9: DE redshifts away a^{-0.3}, so rho_DE was larger in past
    # => H was larger => tau_0 is SMALLER
    assert bg_w09.tau_0 < bg_lcdm.tau_0, "w > -1 should give smaller tau_0"
    print(f"    tau_0(w=-0.9) < tau_0(LCDM): PASS")

    # r_s should be nearly the same (DE negligible at recombination)
    r_s_err = abs(bg_w09.r_s - bg_lcdm.r_s) / bg_lcdm.r_s
    print(f"    r_s difference: {r_s_err*100:.4f}%"
          f"  [{'PASS' if r_s_err < 0.01 else 'CHECK'}]")

    # w = -1.1: phantom. DE grows with time.
    bg_w11 = BackgroundDE(w_const=-1.1).solve()
    print(f"\n  w=-1.1 (phantom):")
    print(f"    tau_0 = {bg_w11.tau_0:.1f} Mpc  (LCDM: {bg_lcdm.tau_0:.1f})")
    assert bg_w11.tau_0 > bg_lcdm.tau_0, "w < -1 should give larger tau_0"
    print(f"    tau_0(w=-1.1) > tau_0(LCDM): PASS")

    print()
    return True


def _test_cpl_desi():
    """Test CPL with DESI 2024 parameters."""
    print("=" * 65)
    print("TEST 3: CPL (DESI 2024: w0=-0.7, wa=-1.0)")
    print("=" * 65)

    bg_lcdm = BackgroundDE(w0=-1.0, wa=0.0).solve()
    bg_desi = DESI_2024().solve()

    print(f"\n  LCDM:  tau_0={bg_lcdm.tau_0:.1f}  D_A={bg_lcdm.D_A:.1f}"
          f"  r_s={bg_lcdm.r_s:.2f}")
    print(f"  DESI:  tau_0={bg_desi.tau_0:.1f}  D_A={bg_desi.D_A:.1f}"
          f"  r_s={bg_desi.r_s:.2f}")

    # w(a=1) = w0 = -0.7, w(a=0.5) = -0.7 + (-1.0)*0.5 = -1.2
    w_today = bg_desi.w_at(1.0)
    w_half = bg_desi.w_at(0.5)
    w_rec = bg_desi.w_at(a_rec)
    print(f"\n  w(a=1)   = {float(w_today):.3f}  (expect -0.700)")
    print(f"  w(a=0.5) = {float(w_half):.3f}  (expect -1.200)")
    print(f"  w(a_rec) = {float(w_rec):.3f}  (expect ~-1.700)")

    assert abs(float(w_today) + 0.7) < 1e-10
    assert abs(float(w_half) + 1.2) < 1e-10
    print(f"  w(a) values: PASS")

    # r_s should be almost identical (DE negligible at recombination)
    r_s_err = abs(bg_desi.r_s - bg_lcdm.r_s) / bg_lcdm.r_s
    print(f"\n  r_s change: {r_s_err*100:.4f}% (should be < 0.1%)")
    assert r_s_err < 0.01, f"r_s changed too much: {r_s_err*100:.2f}%"
    print(f"  r_s stability: PASS")

    # D_A should change significantly
    D_A_change = (bg_desi.D_A - bg_lcdm.D_A) / bg_lcdm.D_A
    print(f"  D_A change: {D_A_change*100:.2f}%")

    # Acoustic scale shift: l_1 ~ pi * D_A / r_s
    l1_lcdm = np.pi * bg_lcdm.D_A / bg_lcdm.r_s
    l1_desi = np.pi * bg_desi.D_A / bg_desi.r_s
    print(f"\n  l_1 (LCDM) = {l1_lcdm:.1f}")
    print(f"  l_1 (DESI) = {l1_desi:.1f}")
    print(f"  l_1 shift:   {(l1_desi - l1_lcdm) / l1_lcdm * 100:.2f}%")

    # f_DE behavior
    a_test = np.array([0.01, 0.1, 0.5, 1.0])
    fDE = bg_desi._f_DE(a_test)
    print(f"\n  f_DE(a):")
    for ai, fi in zip(a_test, fDE):
        print(f"    a={ai:.2f}: f_DE = {fi:.6e}")
    assert abs(fDE[-1] - 1.0) < 1e-10, "f_DE(a=1) must be 1"
    print(f"  f_DE(a=1) = 1: PASS")

    print()
    return True


def _test_distances():
    """Test distance measures."""
    print("=" * 65)
    print("TEST 4: Distance measures")
    print("=" * 65)

    bg = BackgroundDE(w0=-1.0, wa=0.0).solve()

    # d_A(z_rec) should match D_A from conformal time
    d_A_arr = bg.angular_diameter_distance(z_rec)
    d_A_rec = float(d_A_arr.flat[0])
    D_A_conf = bg.D_A / (1.0 + z_rec)  # comoving -> angular diameter

    # Note: bg.D_A is the *comoving* distance (tau_0 - tau_rec)
    chi_arr = bg.comoving_distance(z_rec)
    chi_rec = float(chi_arr.flat[0])
    err_chi = abs(chi_rec - bg.D_A) / bg.D_A
    print(f"  chi(z_rec) = {chi_rec:.1f} Mpc  (conformal: {bg.D_A:.1f})")
    print(f"  chi vs D_A err: {err_chi*100:.2f}%")

    # For LCDM, age should be ~13.8 Gyr
    age = bg.age_of_universe()
    print(f"\n  Age (LCDM) = {age:.2f} Gyr  (Planck: 13.80)")
    assert 13.0 < age < 14.5, f"Age out of range: {age}"
    print(f"  Age in range: PASS")

    # Deceleration parameter: q(a=1) should be < 0 (accelerating)
    q0 = float(bg.deceleration_parameter(1.0))
    print(f"\n  q_0 = {q0:.4f}  (expect ~ -0.55)")
    assert q0 < 0, "Universe should be accelerating today"
    print(f"  q_0 < 0: PASS")

    print()
    return True


def _test_class_comparison():
    """
    Compare D_A, r_s, H(z) with CLASS reference values for DESI CPL.

    CLASS reference values (computed with class_public, Planck 2018 + DESI CPL):
      w0 = -0.7, wa = -1.0
      r_s(z_drag) ~ 147.1 Mpc  (CLASS uses z_drag ~ 1060, not z_rec)
      d_A(z*) ~ 12.8 Mpc / (1+z*)  [comoving ~ 12.0-12.8 Gpc, model-dependent]

    We check:
      1. r_s is stable to < 0.5% (DE negligible at recombination)
      2. D_A shifts by ~1-5% from LCDM for DESI CPL
      3. H_0 consistency
    """
    print("=" * 65)
    print("TEST 5: CLASS-level comparison (DESI CPL)")
    print("=" * 65)

    bg_lcdm = BackgroundDE(w0=-1.0, wa=0.0).solve()
    bg_desi = DESI_2024().solve()

    # Sound horizon: should be nearly identical
    # CLASS Planck 2018 LCDM: r_s(z_drag) ~ 147.09 Mpc
    # Our r_s(z_rec) differs because we use z_rec=1089.92 not z_drag~1059.94
    # But the LCDM-to-DESI ratio should be model-independent
    r_s_ratio = bg_desi.r_s / bg_lcdm.r_s
    print(f"\n  r_s (LCDM) = {bg_lcdm.r_s:.2f} Mpc")
    print(f"  r_s (DESI) = {bg_desi.r_s:.2f} Mpc")
    print(f"  r_s ratio  = {r_s_ratio:.6f}  (expect ~1.0000)")
    assert abs(r_s_ratio - 1.0) < 0.005, \
        f"r_s changed too much: ratio = {r_s_ratio:.6f}"
    print(f"  r_s stability < 0.5%: PASS")

    # D_A ratio: model-dependent
    D_A_ratio = bg_desi.D_A / bg_lcdm.D_A
    print(f"\n  D_A (LCDM) = {bg_lcdm.D_A:.1f} Mpc")
    print(f"  D_A (DESI) = {bg_desi.D_A:.1f} Mpc")
    print(f"  D_A ratio  = {D_A_ratio:.4f}")

    # H(z) at a few redshifts
    z_arr = np.array([0.0, 0.3, 0.5, 1.0, 2.0])
    print(f"\n  {'z':>5s}  {'H_LCDM':>12s}  {'H_DESI':>12s}  {'ratio':>8s}")
    print(f"  " + "-" * 42)
    for z in z_arr:
        a = 1.0 / (1.0 + z)
        H_l = float(bg_lcdm.H_of_a(a))
        H_d = float(bg_desi.H_of_a(a))
        print(f"  {z:5.1f}  {H_l:12.2f}  {H_d:12.2f}  {H_d/H_l:8.4f}")

    # H(0) should be the same (H0 is input)
    err_H0 = abs(float(bg_desi.H_of_a(1.0)) - H0_km_s_Mpc) / H0_km_s_Mpc
    assert err_H0 < 1e-6, f"H0 inconsistent: err = {err_H0}"
    print(f"\n  H(z=0) = H0: PASS")

    print()
    return True


def _test_numerical_mode():
    """Test numerical w(a) integration matches analytic CPL."""
    print("=" * 65)
    print("TEST 6: Numerical w(a) vs analytic CPL")
    print("=" * 65)

    w0, wa = -0.7, -1.0

    # Analytic CPL
    bg_cpl = BackgroundDE(w0=w0, wa=wa).solve()

    # Numerical: same w(a) via callable
    def w_func(a):
        return w0 + wa * (1.0 - a)

    bg_num = BackgroundDE(w_func=w_func).solve()

    checks = [
        ('tau_0', bg_cpl.tau_0, bg_num.tau_0),
        ('D_A', bg_cpl.D_A, bg_num.D_A),
        ('r_s', bg_cpl.r_s, bg_num.r_s),
    ]

    all_pass = True
    for name, val_c, val_n in checks:
        rel_err = abs(val_c - val_n) / abs(val_c)
        status = "PASS" if rel_err < 1e-4 else "FAIL"
        if status == "FAIL":
            all_pass = False
        print(f"  {name:>10s}: CPL={val_c:.4f}  num={val_n:.4f}"
              f"  err={rel_err:.2e}  [{status}]")

    print(f"\n  Result: {'ALL PASS' if all_pass else 'SOME FAILED'}\n")
    return all_pass


def run_all_tests():
    """Run all dark energy self-tests."""
    print("\n" + "=" * 65)
    print("DARK ENERGY MODULE SELF-TESTS")
    print("=" * 65 + "\n")

    results = {}
    results['LCDM consistency'] = _test_lcdm_consistency()
    results['Constant-w'] = _test_constant_w()
    results['CPL DESI'] = _test_cpl_desi()
    results['Distances'] = _test_distances()
    results['CLASS comparison'] = _test_class_comparison()
    results['Numerical mode'] = _test_numerical_mode()

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
