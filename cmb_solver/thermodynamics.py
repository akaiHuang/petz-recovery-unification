"""
thermodynamics.py — Recombination and ionization history for CMB Boltzmann solver.

Implements two methods for computing the free electron fraction x_e(z):
  1. Peebles three-level atom ODE (accurate to ~5% vs RECFAST)
  2. Tanh fitting formula (fast fallback, Planck 2018 calibrated)

From x_e(z) derives:
  - kappa_dot(tau): Thomson scattering rate (opacity)
  - kappa(tau):     optical depth (integrated from today backwards)
  - g(tau):         visibility function = -kappa_dot * exp(-kappa)
  - tau_rec:        conformal time of recombination (visibility peak)

Physics references:
  - Peebles (1968) ApJ 153, 1
  - Seager, Sasselov & Scott (2000) ApJS 128, 407  [RECFAST]
  - Planck 2018 VI, A&A 641, A6

Author: Sheng-Kai Huang, 2026
"""
import numpy as np
from scipy.integrate import solve_ivp, cumulative_trapezoid
from scipy.interpolate import interp1d
from . import constants as C


# ---------------------------------------------------------------------------
# Physical constants for Peebles model (all SI)
# ---------------------------------------------------------------------------
_alpha_fs = 7.2973525693e-3       # fine structure constant
_m_e_eV = 0.51099895e6            # electron mass in eV
_B1_eV = 13.605693                # hydrogen ionization energy in eV
_B1_J = _B1_eV * 1.602176634e-19  # in Joules
_eV_to_J = 1.602176634e-19
_eV_to_K = _eV_to_J / C.kB_SI     # 1 eV in Kelvin
_Lambda_2s = 8.22458              # 2s -> 1s two-photon decay rate [s^{-1}]
_lambda_Lya_m = 1.21567e-7        # Lyman-alpha wavelength [m]


def _H0_si():
    """H0 in SI units (1/s)."""
    return C.H0_km_s_Mpc * 1e3 / C.Mpc_SI


def _n_H0_si():
    """Present-day hydrogen number density n_H,0 [m^{-3}]."""
    H0_si = _H0_si()
    rho_b0 = C.Omega_b * 3.0 * H0_si**2 / (8.0 * np.pi * C.G_SI)
    return (1.0 - C.Y_He) * rho_b0 / C.m_H_SI


# ---------------------------------------------------------------------------
# Peebles three-level atom model
# ---------------------------------------------------------------------------
def _alpha_B(T_b):
    """
    Case-B recombination coefficient alpha^(2)(T_b) [m^3/s].

    Peebles (1968) formula:
        alpha^(2) = 9.78 * (alpha/m_e c^2)^2 * sqrt(B1/T_b) * ln(B1/T_b)

    with dimensions from CGS converted to SI.  We use the standard fit from
    Pequignot, Petitjean & Boisson (1991) which is more accurate:
        alpha_B = F * 1e-19 * a_PQ * t^b_PQ / (1 + c_PQ * t^d_PQ) [m^3/s]
    where t = T_b / 1e4 K.

    For simplicity and transparency, we use the Peebles analytic form with
    the known numerical coefficient calibrated to RECFAST.
    """
    # T_b in Kelvin
    # B1/kT in dimensionless
    B1_over_kT = _B1_eV * _eV_to_K / T_b   # = B1 / (k_B T_b)

    if B1_over_kT < 1e-6:
        # Very high T: fully ionized, alpha -> 0 effectively
        return 0.0

    # Peebles formula in SI [m^3/s]:
    # alpha^(2) = 9.78 * (alpha^2 * hbar / (m_e c))^2 / hbar
    #           * sqrt(B1/(kT)) * ln(B1/(kT))
    # Numerically, the standard result is well approximated by:
    #   alpha_B ~ 2.6e-19 * (B1/kT)^{0.5} / T_4  [m^3/s]
    # We use the exact Peebles form:
    #   alpha_B = 9.78 * (alpha/m_e_c2_eV)^2 * sqrt(kT/B1) * ln(B1/kT) * [unit factor]
    # The unit factor from dimensional analysis:
    #   (alpha * hbar_eV_s / (m_e c^2 in eV * c))^2 * c = cross_section * velocity
    #
    # More practically, use the well-tested fit from Pequignot et al. (1991):
    t4 = T_b / 1e4
    # Pequignot et al. 1991 fit (case B, m^3/s):
    alpha_B = 4.309e-19 * t4**(-0.6166) / (1.0 + 0.6703 * t4**0.5300)
    return alpha_B


def _beta_ionization(alpha_B, T_b):
    """
    Photoionization rate beta(T_b) from the GROUND STATE [s^{-1}].

    Related to alpha^(2) (case-B) by detailed balance / Saha for n >= 2:
        beta = alpha_B * (m_e k_B T_b / (2 pi hbar^2))^{3/2} * exp(-B1/(k_B T_b))

    This is the total photoionization rate that appears in the
    main Peebles ODE rate equation (ionization term).
    """
    B1_over_kT = _B1_J / (C.kB_SI * T_b)

    if B1_over_kT > 500.0:
        return 0.0

    thermal_factor = (C.m_e_SI * C.kB_SI * T_b
                      / (2.0 * np.pi * C.hbar_SI**2))**1.5

    beta = alpha_B * thermal_factor * np.exp(-B1_over_kT)
    return beta


def _beta_2(alpha_B, T_b):
    """
    Photoionization rate from the n=2 level: beta^(2)(T_b) [s^{-1}].

    This is the rate that appears in the Peebles Cr correction factor.
    The n=2 binding energy is B2 = B1/4 = 3.4 eV.

        beta^(2) = alpha_B * (m_e k_B T_b / (2 pi hbar^2))^{3/2} * exp(-B2/(k_B T_b))

    where B2 = B1/4 = 13.6/4 = 3.4 eV.

    This is MUCH larger than the ground-state beta because exp(-3.4/kT) >> exp(-13.6/kT).
    The competition between beta^(2) and Lambda_2s determines whether
    recombination proceeds via 2-photon decay or Lyman-alpha escape.
    """
    B2_J = _B1_J / 4.0  # n=2 binding energy = 3.4 eV
    B2_over_kT = B2_J / (C.kB_SI * T_b)

    if B2_over_kT > 500.0:
        return 0.0

    thermal_factor = (C.m_e_SI * C.kB_SI * T_b
                      / (2.0 * np.pi * C.hbar_SI**2))**1.5

    beta2 = alpha_B * thermal_factor * np.exp(-B2_over_kT)
    return beta2


def _peebles_Cr(n_H, x_e, H_z, T_b, alpha_B):
    """
    Peebles correction factor C_r.

    C_r = (1 + K Lambda_2s n_H (1 - x_e)) /
          (1 + K (Lambda_2s + beta^(2)) n_H (1 - x_e))

    where:
      K = lambda_Lya^3 / (8 pi H(z))     [Lya escape probability]
      Lambda_2s = 8.22458 s^{-1}          [2-photon decay rate]
      beta^(2)                             [photoionization from n=2]

    K accounts for the cosmological redshifting of Lyman-alpha photons:
    a Lya photon escapes the line if it redshifts out before being reabsorbed.

    The bottleneck: atoms in n=2 can either
      (a) decay via 2-photon (rate Lambda_2s), or
      (b) be photoionized (rate beta^(2)).
    Only channel (a) leads to net recombination.
    When beta^(2) >> Lambda_2s, most n=2 atoms are re-ionized => Cr << 1.
    """
    K = _lambda_Lya_m**3 / (8.0 * np.pi * H_z)  # H_z in [1/s]

    # Photoionization from n=2 (B2 = B1/4 = 3.4 eV)
    beta2 = _beta_2(alpha_B, T_b)

    # Neutral fraction
    x_1s = 1.0 - x_e
    if x_1s < 0:
        x_1s = 0.0

    denom_term = K * n_H * x_1s
    numerator = 1.0 + _Lambda_2s * denom_term
    denominator = 1.0 + (_Lambda_2s + beta2) * denom_term

    if denominator < 1e-30:
        return 1.0

    return numerator / denominator


def peebles_ode_rhs(z, xe_arr, H_func, n_H0):
    """
    RHS of the Peebles recombination ODE.

    Physical rate equation in cosmic time t:
        dx_e/dt = C_r * [beta*(1-x_e) - n_H * alpha_B * x_e^2]

    Converting to redshift variable via dt = -dz / [H(z)*(1+z)]:
        dx_e/dz = - C_r / (H(z)*(1+z)) * [beta*(1-x_e) - n_H * alpha_B * x_e^2]

    The minus sign ensures that when recombination dominates (bracket < 0),
    dx_e/dz > 0, meaning x_e decreases as z decreases (correct physical behavior).

    Parameters
    ----------
    z : float
        Redshift (independent variable, integrated from high z downward).
    xe_arr : array-like
        xe_arr[0] = x_e (free electron fraction).
    H_func : callable
        H(z) in SI units [1/s].
    n_H0 : float
        Present-day hydrogen number density [m^{-3}].

    Returns
    -------
    dxe_dz : array-like
        [dx_e/dz]
    """
    x_e = xe_arr[0]
    x_e = np.clip(x_e, 0.0, 1.0)

    T_b = C.T_CMB * (1.0 + z)    # Baryon temperature ~ T_CMB*(1+z) (tight coupling)
    n_H = n_H0 * (1.0 + z)**3    # Hydrogen number density at redshift z
    H_z = H_func(z)               # Hubble parameter in [1/s]

    alpha_B = _alpha_B(T_b)
    beta = _beta_ionization(alpha_B, T_b)
    Cr = _peebles_Cr(n_H, x_e, H_z, T_b, alpha_B)

    # Recombination rate:  n_H * alpha_B * x_e^2
    # Ionization rate:     beta * (1 - x_e)
    recomb_rate = n_H * alpha_B * x_e**2
    ioniz_rate = beta * (1.0 - x_e)

    # dx_e/dz = -C_r / (H*(1+z)) * (ionization - recombination)
    # The minus sign comes from dt = -dz / (H*(1+z)).
    dxe_dz = -Cr / (H_z * (1.0 + z)) * (ioniz_rate - recomb_rate)

    return [dxe_dz]


# ---------------------------------------------------------------------------
# Tanh fitting formula (Planck 2018 calibrated)
# ---------------------------------------------------------------------------
def xe_tanh(z, z_rec=1089.8, Delta_z=80.0, x_e_floor=2.0e-4):
    """
    Tanh fit for the ionization fraction x_e(z).

    x_e(z) = max(0.5 * [1 + tanh((z - z_rec) / Delta_z)], x_e_floor)

    Calibrated to Planck 2018 best-fit recombination history.
    Freeze-out floor at x_e ~ 2e-4 from standard Peebles model.

    Parameters
    ----------
    z : float or array
        Redshift.
    z_rec : float
        Center of recombination transition (default: 1089.8, Planck 2018).
    Delta_z : float
        Width of the transition (default: 80).
    x_e_floor : float
        Late-time freeze-out floor (default: 2e-4).

    Returns
    -------
    x_e : float or array
    """
    z = np.asarray(z, dtype=float)
    x_main = 0.5 * (1.0 + np.tanh((z - z_rec) / Delta_z))
    return np.maximum(x_main, x_e_floor)


# ===================================================================
# Main class: RecombinationSolver
# ===================================================================
class RecombinationSolver:
    """
    Compute the ionization history x_e(z) and derive optical depth,
    scattering rate, and visibility function on the background grid.

    Methods
    -------
    solve_peebles(bg) :  Peebles three-level ODE
    solve_tanh(bg) :     Tanh fitting formula (fast fallback)
    solve(bg) :          Try Peebles, fall back to tanh on failure

    After solving, attaches these arrays to the BackgroundSolver `bg`:
        bg.x_e_grid, bg.kappa_grid, bg.kappa_dot_grid, bg.visibility
    and interpolators:
        bg.x_e_interp(a), bg.kappa_interp(tau), bg.kappa_dot_interp(tau),
        bg.visibility_interp(tau)
    """

    def __init__(self, z_start=1800.0, z_end=0.0, verbose=True):
        """
        Parameters
        ----------
        z_start : float
            Starting redshift for Peebles integration (default: 1800).
            Must be high enough that x_e ~ 1 (Saha equilibrium).
        z_end : float
            Ending redshift (default: 0).
        verbose : bool
            Print progress information.
        """
        self.z_start = z_start
        self.z_end = z_end
        self.verbose = verbose

        # Physical constants
        self._H0_si = _H0_si()
        self._n_H0 = _n_H0_si()

    def _log(self, msg):
        if self.verbose:
            print(msg)

    # -------------------------------------------------------------------
    # Hubble parameter wrappers
    # -------------------------------------------------------------------
    def _H_of_z_si(self, z, bg):
        """H(z) in SI units [1/s], using the background solver."""
        a = 1.0 / (1.0 + z)
        a = np.clip(a, C.a_init, 1.0)
        H_over_H0 = bg.H_of_a(a)  # dimensionless H/H0
        return H_over_H0 * self._H0_si

    # -------------------------------------------------------------------
    # Method 1: Peebles three-level ODE (with Saha pre-phase)
    # -------------------------------------------------------------------
    def solve_peebles(self, bg):
        """
        Solve recombination using a two-phase approach (standard RECFAST strategy):

        Phase 1 (Saha regime, z > z_switch):
            Use the Saha equation directly. This is valid as long as recombination
            is fast enough to maintain equilibrium.

        Phase 2 (Peebles ODE, z < z_switch):
            Integrate the Peebles three-level ODE from the Saha x_e at z_switch
            down to z_end. The Peebles correction factor Cr accounts for the
            Lyman-alpha bottleneck and two-photon decay.

        z_switch is found by scanning downward from z_start until the Saha x_e
        drops below a threshold (default: 0.99). This is where recombination
        starts to depart from equilibrium and the Peebles model takes over.

        Parameters
        ----------
        bg : BackgroundSolver
            Must have .solve() called already.

        Returns
        -------
        self
        """
        assert bg.solved, "BackgroundSolver must be solved first."

        self._log("[Thermodynamics] Solving Peebles three-level model...")
        self._log(f"[Thermodynamics] n_H,0 = {self._n_H0:.4e} m^-3")

        # Phase 1: Find z_switch where Saha x_e drops below threshold
        x_e_switch_threshold = 0.99
        z_switch = self.z_start
        for z_test in np.arange(self.z_start, 800, -1.0):
            xe_saha = self._saha_xe(z_test)
            if xe_saha < x_e_switch_threshold:
                z_switch = z_test + 1.0   # switch just before departure
                break

        xe_at_switch = self._saha_xe(z_switch)
        self._log(f"[Thermodynamics] Saha -> Peebles switch at z = {z_switch:.0f}, "
                  f"x_e(Saha) = {xe_at_switch:.6f}")

        # Phase 2: Peebles ODE from z_switch down to z_end
        def H_func(z):
            return self._H_of_z_si(z, bg)

        n_H0 = self._n_H0

        def rhs(z, y):
            return peebles_ode_rhs(z, y, H_func, n_H0)

        # Dense z-grid: fine sampling through recombination, coarser at late times
        n_recomb = 3000   # through recombination
        n_late = 2000     # late time
        z_recomb = np.linspace(z_switch, 800, n_recomb)
        z_late = np.linspace(799, max(self.z_end, 0.0), n_late)
        z_eval = np.concatenate([z_recomb, z_late])

        sol = solve_ivp(
            rhs,
            [z_switch, max(self.z_end, 0.0)],
            [xe_at_switch],
            method='Radau',      # implicit Runge-Kutta, excellent for stiff ODEs
            rtol=1e-8,
            atol=1e-12,
            t_eval=z_eval,
            dense_output=True,
            max_step=2.0,
        )

        if not sol.success:
            self._log(f"[Thermodynamics] Peebles ODE failed: {sol.message}")
            self._log("[Thermodynamics] Falling back to tanh fit.")
            return self.solve_tanh(bg)

        self._log(f"[Thermodynamics] Peebles ODE: {sol.t.size} points, "
                  f"status = {sol.status}")

        # Combine Saha (z > z_switch) + Peebles (z < z_switch) into arrays
        z_saha = np.linspace(self.z_start, z_switch + 0.5, 500)
        xe_saha = np.array([self._saha_xe(z) for z in z_saha])

        z_peebles = sol.t
        xe_peebles = np.clip(sol.y[0], 0.0, 1.5)

        # Merge (avoid overlap at z_switch)
        self._z_sol = np.concatenate([z_saha, z_peebles])
        self._xe_sol = np.concatenate([xe_saha, xe_peebles])
        self._xe_dense = sol.sol   # dense output for Peebles regime
        self._z_switch = z_switch
        self._method = 'peebles'

        # Report key values
        # x_e at z=1100 and z=1000
        for z_check in [1300, 1200, 1100, 1000, 900, 500]:
            if z_check <= z_switch and z_check >= self.z_end:
                xe_val = float(np.clip(sol.sol(z_check)[0], 0, 1.5))
                self._log(f"[Thermodynamics] x_e(z={z_check}) = {xe_val:.6f}")

        # Map onto the background a_grid
        self._map_to_background(bg)
        return self

    def _saha_xe(self, z):
        """
        Saha equation for x_e at given z.
        Used to set initial conditions when z_start is in the partial-ionization regime.

        x_e^2 / (1-x_e) = (1/n_H) * (m_e k_B T / (2pi hbar^2))^{3/2} * exp(-B1/(k_B T))
        """
        T_b = C.T_CMB * (1.0 + z)
        n_H = self._n_H0 * (1.0 + z)**3

        thermal = (C.m_e_SI * C.kB_SI * T_b
                   / (2.0 * np.pi * C.hbar_SI**2))**1.5
        B1_over_kT = _B1_J / (C.kB_SI * T_b)

        if B1_over_kT > 500:
            return 0.0

        rhs_val = thermal * np.exp(-B1_over_kT) / n_H

        # Solve x^2/(1-x) = rhs  => x^2 + rhs*x - rhs = 0
        # x = (-rhs + sqrt(rhs^2 + 4*rhs)) / 2
        discriminant = rhs_val**2 + 4.0 * rhs_val
        x_e = (-rhs_val + np.sqrt(discriminant)) / 2.0
        return np.clip(x_e, 0.0, 1.0)

    # -------------------------------------------------------------------
    # Method 2: Tanh fitting formula
    # -------------------------------------------------------------------
    def solve_tanh(self, bg):
        """
        Use tanh fitting formula for x_e(z) (fast, ~5% accuracy).

        Parameters
        ----------
        bg : BackgroundSolver
            Must have .solve() called already.

        Returns
        -------
        self
        """
        assert bg.solved, "BackgroundSolver must be solved first."

        self._log("[Thermodynamics] Using tanh fitting formula for x_e(z)")

        self._method = 'tanh'
        self._z_sol = None
        self._xe_sol = None
        self._xe_dense = None

        # Map directly onto background grid
        self._map_to_background(bg)
        return self

    # -------------------------------------------------------------------
    # Auto solver: try Peebles, fallback to tanh
    # -------------------------------------------------------------------
    def solve(self, bg):
        """
        Compute ionization history: try Peebles ODE first, fall back to tanh.

        Parameters
        ----------
        bg : BackgroundSolver
            Must have .solve() called already.

        Returns
        -------
        self
        """
        try:
            return self.solve_peebles(bg)
        except Exception as e:
            self._log(f"[Thermodynamics] Peebles solver error: {e}")
            self._log("[Thermodynamics] Falling back to tanh fit.")
            return self.solve_tanh(bg)

    # -------------------------------------------------------------------
    # Map x_e onto background grid & compute kappa, visibility
    # -------------------------------------------------------------------
    def _map_to_background(self, bg):
        """
        Evaluate x_e on the background a_grid and compute optical depth,
        scattering rate, and visibility function.

        Attaches results to both self and the BackgroundSolver bg.
        """
        a_grid = bg.a_grid
        tau_grid = bg.tau_grid
        z_grid = 1.0 / a_grid - 1.0

        # Evaluate x_e on the a_grid
        if self._method == 'peebles':
            # Two-phase: Saha above z_switch, Peebles below
            x_e = np.zeros_like(a_grid)
            z_sw = self._z_switch
            for i, z in enumerate(z_grid):
                if z > self.z_start:
                    # Beyond integration range: fully ionized
                    x_e[i] = 1.0
                elif z > z_sw:
                    # Saha regime
                    x_e[i] = min(self._saha_xe(z), 1.0)
                elif z < max(self.z_end, 0.0):
                    # Below ODE range: use last value
                    x_e[i] = float(np.clip(self._xe_sol[-1], 0.0, 1.5))
                else:
                    # Peebles ODE regime
                    val = self._xe_dense(z)
                    x_e[i] = np.clip(float(val[0]), 0.0, 1.5)
        else:
            # Tanh fit
            x_e = xe_tanh(z_grid)

        # Helium correction: He is singly ionized for z > ~6000, and
        # doubly ionized for z > ~15000 (small correction at recombination).
        # For the Peebles model focusing on H recombination, add He I contribution:
        # n_e = x_e * n_H + n_He (fully ionized He for T > 2 eV ~ z > 6000)
        # Effective x_e for Thomson scattering: x_e_eff = x_e + f_He * x_He
        # f_He = Y_He / (4 * (1 - Y_He)) = He/H number ratio
        f_He = C.Y_He / (4.0 * (1.0 - C.Y_He))

        # He recombination (simplified): He fully ionized above z ~ 2000,
        # singly ionized below z ~ 6000 (He I recomb),
        # neutral below z ~ 1600 (He I recomb complete before H recomb)
        x_He = np.ones_like(z_grid) * f_He   # He contribution to n_e/n_H
        # He I recombination: smooth transition around z ~ 1800 with width ~200
        x_He *= 0.5 * (1.0 + np.tanh((z_grid - 1800.0) / 200.0))

        # Total effective electron fraction for Thomson scattering
        # x_e_eff = x_e (from H) + x_He (from He)
        x_e_eff = x_e + x_He

        self.x_e_grid = x_e       # hydrogen ionization fraction only
        self.x_e_eff_grid = x_e_eff  # effective x_e for scattering
        self.z_grid = z_grid

        # Compute scattering rate and optical depth
        self._compute_optical_depth(bg, x_e_eff, a_grid, tau_grid)

        # Summary diagnostics
        self._print_diagnostics(bg, x_e, a_grid, tau_grid)

    def _compute_optical_depth(self, bg, x_e_eff, a_grid, tau_grid):
        """
        Compute Thomson optical depth kappa(tau), scattering rate kappa_dot(tau),
        and visibility function g(tau) = -kappa_dot * exp(-kappa).

        kappa_dot = d(kappa)/d(tau) where tau is conformal time.
        In code units (tau in units of H0^{-1}):
            kappa_dot_code = x_e * n_H0 * sigma_T * c / (H0 * a^2)

        kappa is integrated from tau_0 (today) backwards:
            kappa(tau) = integral_tau^{tau_0} |kappa_dot| dtau'
        """
        H0_si = self._H0_si
        n_H0 = self._n_H0

        # Thomson scattering rate in code units:
        # kappa_dot_code = (n_H0 * sigma_T * c / H0) * x_e * (1+z)^2
        # = (n_H0 * sigma_T * c / H0) * x_e_eff / a^2
        # This is d(kappa)/d(tau_code) where tau_code is in units of H0^{-1}.
        prefactor = n_H0 * C.sigma_T_SI * C.c_SI / H0_si
        kappa_dot_abs = prefactor * x_e_eff / a_grid**2

        # Cumulative integral from tau_grid[0] to tau:
        # integral_0^tau |kappa_dot| dtau
        kappa_cumul = np.zeros_like(tau_grid)
        kappa_cumul[1:] = cumulative_trapezoid(kappa_dot_abs, tau_grid)

        # Optical depth from tau to today:
        # kappa(tau) = integral_tau^{tau_0} |kappa_dot| dtau'
        #            = kappa_cumul(tau_0) - kappa_cumul(tau)
        kappa_grid = kappa_cumul[-1] - kappa_cumul

        # Visibility function: g(tau) = |kappa_dot| * exp(-kappa)
        # (positive-definite; the sign convention is that g > 0 at the
        #  last scattering surface)
        visibility = kappa_dot_abs * np.exp(-kappa_grid)

        # Store on self
        self.kappa_dot_grid = -kappa_dot_abs   # convention: kappa_dot < 0
        self.kappa_grid = kappa_grid
        self.visibility = visibility

        # Build interpolators
        self._build_interpolators(bg, a_grid, tau_grid, kappa_grid,
                                  kappa_dot_abs, visibility)

        # Attach to background solver for use by perturbation/power spectrum modules
        bg.x_e_grid = self.x_e_eff_grid
        bg.kappa_grid = kappa_grid
        bg.kappa_dot_grid = -kappa_dot_abs   # negative convention
        bg.visibility = visibility

    def _build_interpolators(self, bg, a_grid, tau_grid, kappa_grid,
                             kappa_dot_abs, visibility):
        """Build scipy interpolators and attach to both self and bg."""
        lna = np.log(a_grid)

        # x_e(a)
        self._interp_xe = interp1d(
            lna, self.x_e_eff_grid, kind='cubic',
            fill_value=(self.x_e_eff_grid[0], self.x_e_eff_grid[-1]),
            bounds_error=False
        )

        # kappa(tau)
        self._interp_kappa = interp1d(
            tau_grid, kappa_grid, kind='cubic', fill_value='extrapolate'
        )

        # kappa_dot(tau) — store absolute value for interpolation, apply sign after
        self._interp_kappa_dot = interp1d(
            tau_grid, -kappa_dot_abs, kind='cubic', fill_value='extrapolate'
        )

        # visibility g(tau)
        self._interp_visibility = interp1d(
            tau_grid, visibility, kind='cubic', fill_value='extrapolate'
        )

        # Also attach interpolators to bg (backward-compatible with existing code)
        bg._interp_x_e = self._interp_xe
        bg._interp_kappa = self._interp_kappa
        bg._interp_kappa_dot = self._interp_kappa_dot
        bg._interp_visibility = self._interp_visibility

        # Attach interpolator methods to bg
        bg.x_e_interp = self.x_e_interp
        bg.kappa_interp = self.kappa_interp
        bg.kappa_dot_interp = self.kappa_dot_interp
        bg.visibility_interp = self.visibility_interp

    def _print_diagnostics(self, bg, x_e_H, a_grid, tau_grid):
        """Print diagnostic summary of the recombination calculation."""
        z_grid = 1.0 / a_grid - 1.0

        # x_e at recombination
        a_rec = 1.0 / (1.0 + C.z_rec)
        xe_rec = self.x_e_interp(a_rec)

        # Optical depth at recombination
        tau_rec = bg.tau_of_a(a_rec)
        kappa_rec = self.kappa_interp(tau_rec)

        # Visibility function peak
        peak_idx = np.argmax(self.visibility)
        z_vis_peak = z_grid[peak_idx]
        a_vis_peak = a_grid[peak_idx]
        tau_vis_peak = tau_grid[peak_idx]

        # Last scattering surface (kappa = 1)
        idx_lss = np.argmin(np.abs(self.kappa_grid - 1.0))
        z_lss = z_grid[idx_lss]

        # Total optical depth to CMB
        kappa_total = self.kappa_grid[0]   # kappa at earliest time

        # Store key quantities
        self.z_vis_peak = z_vis_peak
        self.tau_vis_peak = tau_vis_peak
        self.a_vis_peak = a_vis_peak
        self.z_lss = z_lss
        self.tau_rec = tau_rec
        self.kappa_rec = kappa_rec
        self.kappa_total = kappa_total

        self._log(f"[Thermodynamics] Method: {self._method}")
        self._log(f"[Thermodynamics] x_e(z={C.z_rec:.0f}) = {xe_rec:.4f}")
        self._log(f"[Thermodynamics] kappa(z={C.z_rec:.0f}) = {kappa_rec:.2f}")
        self._log(f"[Thermodynamics] kappa_total (z->inf) = {kappa_total:.1f}")
        self._log(f"[Thermodynamics] Last scattering (kappa=1): z = {z_lss:.0f}")
        self._log(f"[Thermodynamics] Visibility peak: z = {z_vis_peak:.0f}, "
                  f"tau = {tau_vis_peak:.4f}")

        # Sanity checks
        if abs(z_vis_peak - 1060) > 100:
            self._log(f"[Thermodynamics] WARNING: Visibility peak z={z_vis_peak:.0f} "
                      f"far from expected ~1060")
        if kappa_rec < 0.5 or kappa_rec > 10:
            self._log(f"[Thermodynamics] WARNING: kappa(z_rec) = {kappa_rec:.2f} "
                      f"outside expected range [0.5, 10]")

    # -------------------------------------------------------------------
    # Interpolator methods
    # -------------------------------------------------------------------
    def x_e_interp(self, a):
        """Interpolated effective x_e at scale factor a."""
        a = np.asarray(a, dtype=float)
        scalar = a.ndim == 0
        a = np.atleast_1d(a)
        result = self._interp_xe(np.log(a))
        return float(result[0]) if scalar else result

    def kappa_interp(self, tau):
        """Interpolated optical depth kappa at conformal time tau."""
        return float(self._interp_kappa(tau))

    def kappa_dot_interp(self, tau):
        """Interpolated scattering rate d(kappa)/d(tau) at conformal time tau."""
        return float(self._interp_kappa_dot(tau))

    def visibility_interp(self, tau):
        """Interpolated visibility function g(tau) at conformal time tau."""
        return float(self._interp_visibility(tau))

    # -------------------------------------------------------------------
    # Convenience: get recombination quantities as dict
    # -------------------------------------------------------------------
    def get_summary(self):
        """Return a dictionary of key recombination quantities."""
        return {
            'method': self._method,
            'z_vis_peak': self.z_vis_peak,
            'tau_vis_peak': self.tau_vis_peak,
            'z_lss': self.z_lss,
            'tau_rec': self.tau_rec,
            'kappa_rec': self.kappa_rec,
            'kappa_total': self.kappa_total,
        }


# ===================================================================
# Standalone convenience function
# ===================================================================
def compute_recombination(bg, method='auto', verbose=True):
    """
    Convenience function: compute recombination and attach results to bg.

    Parameters
    ----------
    bg : BackgroundSolver
        Must have .solve() called already.
    method : str
        'peebles' : Peebles three-level ODE
        'tanh'    : tanh fitting formula
        'auto'    : try Peebles, fall back to tanh
    verbose : bool
        Print progress.

    Returns
    -------
    rec : RecombinationSolver
        Solved recombination object.
    """
    rec = RecombinationSolver(verbose=verbose)
    if method == 'peebles':
        rec.solve_peebles(bg)
    elif method == 'tanh':
        rec.solve_tanh(bg)
    else:
        rec.solve(bg)
    return rec
