"""
background.py — Background cosmology solver for Khronon (Sigma = 2 ln Q).

Physics of ghost condensation (K'(Q_0) = 0):
  - The Khronon field sits at the minimum of K(Q), so Q_0 = const.
  - delta_bg = Q_0 - 1 = const.
  - rho_K(a) = rho_K,0 / a^3  — EXACTLY CDM at background level.
  - All nontrivial physics enters at perturbation level (c_s^2 = 0, w_K).

Friedmann equation:
    H^2(a) = H0^2 [Omega_r a^{-4} + Omega_m a^{-3} + Omega_Lambda]
    where Omega_m = Omega_b + Omega_K (identical to LCDM background).
"""
import numpy as np
from scipy.integrate import cumulative_trapezoid, solve_ivp
from scipy.interpolate import interp1d
from . import constants as C


def _H0_si():
    """H0 in SI (1/s)."""
    return C.H0_km_s_Mpc * 1e3 / C.Mpc_SI


class BackgroundSolver:
    """
    Solve the background Friedmann equation for Khronon cosmology.
    At the background level, the Khronon is exactly CDM (ghost condensation).
    """

    def __init__(self, a_init=None, a_final=None, N_points=10000):
        self.a_init = a_init or C.a_init
        self.a_final = a_final or C.a_today
        self.N_points = N_points
        self._compute_Q0()
        self.Omega_m = C.Omega_b + C.Omega_K_khronon
        self.solved = False

    def _compute_Q0(self):
        """Compute background Q_0 from Omega_K and mu."""
        R = 3 * C.H0_invMpc**2 * C.Omega_K_khronon / C.mu**2
        self.delta_bg = -1.0 + np.sqrt(1.0 + R)
        self.Q0 = 1.0 + self.delta_bg
        self.Sigma_bg = 2.0 * np.log(self.Q0)
        self.w_K_bg = self.delta_bg / (2.0 + self.delta_bg)

        print(f"[Background] Ghost condensation: Q_0 = {self.Q0:.6f}")
        print(f"[Background] delta_bg = {self.delta_bg:.6f}, "
              f"Sigma = {self.Sigma_bg:.6f}, w_K = {self.w_K_bg:.6f}")

    def H_of_a(self, a):
        """Hubble parameter H(a)/H0."""
        a = np.asarray(a, dtype=float)
        h2 = C.Omega_r * a**(-4) + self.Omega_m * a**(-3) + C.Omega_Lambda
        return np.sqrt(h2)

    def conformal_hubble(self, a):
        """calH = aH(a)/H0."""
        return a * self.H_of_a(a)

    def rho_K_normalized(self, a):
        """rho_K(a)/rho_crit,0 = Omega_K/a^3."""
        return C.Omega_K_khronon / a**3

    def w_K_of_a(self, a):
        return self.w_K_bg

    def delta_K_interp(self, a):
        return self.delta_bg

    def solve(self):
        """Integrate conformal time tau(a)."""
        a_grid = np.geomspace(self.a_init, self.a_final, self.N_points)
        H_grid = self.H_of_a(a_grid)

        integrand = 1.0 / (a_grid**2 * H_grid)
        tau_grid = np.zeros_like(a_grid)
        tau_grid[1:] = cumulative_trapezoid(integrand, a_grid)

        self.a_grid = a_grid
        self.tau_grid = tau_grid
        self.H_grid = H_grid

        self._build_interpolators()
        self.solved = True

        tau_0_Mpc = self.tau_grid[-1] / C.H0_invMpc
        tau_rec_Mpc = self.tau_of_a(C.a_rec) / C.H0_invMpc
        a_eq = C.Omega_r / self.Omega_m

        print(f"[Background] tau_0 = {tau_0_Mpc:.1f} Mpc, "
              f"tau_rec = {tau_rec_Mpc:.1f} Mpc")
        print(f"[Background] z_eq = {1/a_eq - 1:.0f}")
        return self

    def _build_interpolators(self):
        lna = np.log(self.a_grid)
        self._interp_H = interp1d(lna, self.H_grid, kind='cubic',
                                   fill_value='extrapolate')
        self._interp_tau = interp1d(lna, self.tau_grid, kind='cubic',
                                     fill_value='extrapolate')
        self._interp_a_of_tau = interp1d(self.tau_grid, self.a_grid,
                                          kind='cubic', fill_value='extrapolate')

    def H_interp(self, a):
        return float(self._interp_H(np.log(a)))

    def tau_of_a(self, a):
        return float(self._interp_tau(np.log(a)))

    def a_of_tau(self, tau):
        return float(self._interp_a_of_tau(tau))

    # -----------------------------------------------------------------------
    # Recombination
    # -----------------------------------------------------------------------
    def compute_recombination(self):
        """
        Compute ionization history x_e(z) and optical depth kappa(tau).

        Uses a tanh fitting function calibrated to match RECFAST/Peebles results:
          x_e(z) = 0.5 * [1 + tanh((z - z_rec) / Delta_z)]
        with z_rec = 1089, Delta_z = 80, and a freeze-out floor at x_e ~ 2e-4.

        This gives kappa(z=1090) ~ 2.5, visibility peak at z ~ 1043,
        consistent with standard recombination codes.
        """
        return self._recombination_tanh_fit()

    def _recombination_tanh_fit(self):
        """
        Fallback: use a tanh fitting function for x_e(z).
        Calibrated to match Peebles/RECFAST results.

        x_e(z) = 0.5 * [1 + tanh((z - z_rec) / Delta_z)]

        where z_rec ~ 1089 and Delta_z ~ 80.
        At low z, freeze-out at x_e ~ 2e-4.
        """
        print("[Recombination] Using tanh fit for x_e(z)")

        H0_si = _H0_si()
        rho_b0_SI = C.Omega_b * 3 * H0_si**2 / (8 * np.pi * C.G_SI)
        n_H0_SI = (1 - C.Y_He) * rho_b0_SI / C.m_H_SI

        z_rec = 1089.0
        Delta_z = 80.0
        x_e_floor = 2e-4   # freeze-out

        x_e = np.zeros_like(self.a_grid)
        for i, a in enumerate(self.a_grid):
            z = 1.0 / a - 1.0
            x_main = 0.5 * (1.0 + np.tanh((z - z_rec) / Delta_z))
            x_e[i] = max(x_main, x_e_floor)

        self.x_e_grid = x_e
        self._compute_optical_depth(n_H0_SI, H0_si)
        return self

    def _compute_optical_depth(self, n_H0_SI, H0_si):
        """Compute optical depth kappa(tau) and visibility g(tau) from x_e."""
        # d(kappa)/d(tau_code) = x_e * n_H0 * sigma_T * c / (H0 * a^2)
        kappa_rate = self.x_e_grid * n_H0_SI * C.sigma_T_SI * C.c_SI / (
            H0_si * self.a_grid**2)

        kappa_cumul = cumulative_trapezoid(kappa_rate, self.tau_grid, initial=0)
        kappa_grid = kappa_cumul[-1] - kappa_cumul

        self.kappa_grid = kappa_grid
        self.kappa_dot_grid = -kappa_rate
        self.visibility = kappa_rate * np.exp(-kappa_grid)

        # Interpolators
        lna = np.log(self.a_grid)
        self._interp_x_e = interp1d(lna, self.x_e_grid, kind='cubic',
                                      fill_value=(1.0, self.x_e_grid[-1]),
                                      bounds_error=False)
        self._interp_kappa = interp1d(self.tau_grid, kappa_grid, kind='cubic',
                                       fill_value='extrapolate')
        self._interp_kappa_dot = interp1d(self.tau_grid, self.kappa_dot_grid,
                                           kind='cubic', fill_value='extrapolate')
        self._interp_visibility = interp1d(self.tau_grid, self.visibility,
                                            kind='cubic', fill_value='extrapolate')

        # Summary
        tau_rec = self.tau_of_a(C.a_rec)
        kappa_rec = self.kappa_interp(tau_rec)
        peak_idx = np.argmax(self.visibility)
        z_peak = 1.0 / self.a_grid[peak_idx] - 1
        idx_lss = np.argmin(np.abs(kappa_grid - 1.0))
        z_lss = 1.0 / self.a_grid[idx_lss] - 1

        print(f"[Recombination] x_e(z=1090) = {self.x_e_interp(C.a_rec):.4f}")
        print(f"[Recombination] kappa(z=1090) = {kappa_rec:.2f}")
        print(f"[Recombination] Last scattering (kappa=1): z = {z_lss:.0f}")
        print(f"[Recombination] Visibility peak: z = {z_peak:.0f}")

    def x_e_interp(self, a):
        return float(self._interp_x_e(np.log(a)))

    def kappa_interp(self, tau):
        return float(self._interp_kappa(tau))

    def kappa_dot_interp(self, tau):
        return float(self._interp_kappa_dot(tau))

    def visibility_interp(self, tau):
        return float(self._interp_visibility(tau))
