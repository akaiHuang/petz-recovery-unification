"""
background.py — Background cosmology for Khronon (Sigma = 2 ln Q).

Ghost condensation: K'(Q_0) = 0 => rho_K ~ a^{-3} exactly (CDM-like background).
Friedmann: H^2(a) = H0^2 [Omega_r/a^4 + Omega_m/a^3 + Omega_Lambda]

This is fast numpy code (~10 ms). No GPU needed.
"""
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import interp1d


# ============================================================================
# Cosmological parameters (Planck 2018)
# ============================================================================
h = 0.6736
H0_km_s_Mpc = h * 100.0
c_km_s = 299792.458
H0_Mpc = H0_km_s_Mpc / c_km_s          # H0 in 1/Mpc

omega_b = 0.02237
omega_c = 0.12
Omega_b = omega_b / h**2
Omega_c = omega_c / h**2                 # Standard CDM (or Khronon)
Omega_r = 2.469e-5 * (1 + 0.2271 * 3.046) / h**2   # photons + 3.046 neutrinos
Omega_m = Omega_b + Omega_c
Omega_L = 1.0 - Omega_m - Omega_r

T_CMB = 2.7255   # K
Y_He = 0.2454

# Derived
a_eq = Omega_r / Omega_m                 # matter-radiation equality
k_eq = np.sqrt(2 * Omega_m * H0_Mpc**2 / a_eq)

# Primordial spectrum
A_s = 2.1e-9
n_s = 0.9649
k_pivot = 0.05   # Mpc^{-1}

# Recombination
z_rec = 1089.92
a_rec = 1.0 / (1.0 + z_rec)

# Thomson scattering
sigma_T_SI = 6.6524587e-29   # m^2
c_SI = 2.99792458e8          # m/s
Mpc_SI = 3.0856775814913673e22
m_H_SI = 1.6735575e-27       # kg
G_SI = 6.67430e-11
H0_SI = H0_km_s_Mpc * 1e3 / Mpc_SI  # 1/s


class Background:
    """Solve background Friedmann + recombination. All in Mpc units."""

    def __init__(self, N_points=50000, a_min=1e-7, a_max=1.0, khronon=False,
                 recombination='tanh'):
        """
        Parameters
        ----------
        khronon : bool
            If True, use Omega_K_khronon = 0.265 replacing CDM.
            If False, standard LCDM.
        recombination : str
            'tanh'    : fast tanh fit (~5% accuracy in x_e). Default.
            'peebles' : Peebles three-level atom ODE (~0.1% accuracy in x_e).
                        Better visibility function, peak heights, and damping tail.
            'recfast' : RECFAST-like solver with helium, T_b evolution, and
                        enhanced fudge factors. Visibility FWHM ~ 20-25 Mpc.
                        Recommended for CMB precision work.
        """
        self.N = N_points
        self.a_min = a_min
        self.a_max = a_max
        self.khronon = khronon
        self.recombination_method = recombination

        if khronon:
            self.Omega_cdm = 0.265        # Khronon sector
        else:
            self.Omega_cdm = Omega_c

        self.Omega_m = Omega_b + self.Omega_cdm
        self.Omega_L = 1.0 - self.Omega_m - Omega_r

    def solve(self):
        """Compute conformal time tau(a), sound horizon, recombination."""
        a = np.logspace(np.log10(self.a_min), np.log10(self.a_max), self.N)
        E = np.sqrt(Omega_r / a**4 + self.Omega_m / a**3 + self.Omega_L)

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
        self.D_A = self.tau_0 - self.tau_rec  # comoving distance to LSS

        # Baryon loading R = 3 rho_b / (4 rho_gamma) = 3 Omega_b a / (4 Omega_gamma)
        # With neutrinos: Omega_gamma = Omega_r / (1 + 0.2271 * N_eff)
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
        if self.recombination_method == 'recfast':
            self._compute_recombination_recfast()
        elif self.recombination_method == 'peebles':
            self._compute_recombination_peebles()
        else:
            self._compute_recombination()

        return self

    def _build_interpolators(self):
        lna = np.log(self.a_grid)
        self._a_of_tau = interp1d(self.tau_grid, self.a_grid,
                                   kind='cubic', fill_value='extrapolate')
        self._tau_of_a = interp1d(lna, self.tau_grid,
                                   kind='cubic', fill_value='extrapolate')
        self._E_of_a = interp1d(lna, self.E_grid,
                                 kind='cubic', fill_value='extrapolate')
        self._calH_of_tau = interp1d(self.tau_grid, self.calH_grid,
                                      kind='cubic', fill_value='extrapolate')
        self._R_of_tau = interp1d(self.tau_grid,
                                   self.R_grid,
                                   kind='cubic', fill_value='extrapolate')

    def _compute_recombination(self):
        """Tanh fit for x_e(z), calibrated to RECFAST."""
        z_center = 1089.0
        Delta_z = 80.0
        x_e_floor = 2e-4

        z_grid = 1.0 / self.a_grid - 1.0
        x_e = np.maximum(0.5 * (1.0 + np.tanh((z_grid - z_center) / Delta_z)), x_e_floor)
        self.x_e_grid = x_e

        # Thomson scattering rate: d(kappa)/d(tau) in Mpc^{-1}
        rho_b0_SI = Omega_b * 3 * H0_SI**2 / (8 * np.pi * G_SI)
        n_H0 = (1 - Y_He) * rho_b0_SI / m_H_SI

        # kappa_dot = - x_e * n_H * sigma_T * c / (H0 * a^2)  [in code units: tau in Mpc]
        # More precisely: d(kappa)/d(tau_Mpc) = x_e * n_H0 * sigma_T * (c / H0_Mpc_to_si) / a^2
        # = x_e * n_H0 * sigma_T * Mpc_SI / a^2
        kappa_rate = x_e * n_H0 * sigma_T_SI * Mpc_SI / self.a_grid**2

        # Optical depth: kappa(tau) = integral from tau to tau_0 of kappa_rate dtau
        kappa_cumul = cumulative_trapezoid(kappa_rate, self.tau_grid, initial=0.0)
        kappa = kappa_cumul[-1] - kappa_cumul

        self.kappa_grid = kappa
        self.kappa_dot_grid = -kappa_rate   # negative by convention
        self.visibility_grid = kappa_rate * np.exp(-kappa)

        # Interpolators
        self._kappa_dot_of_tau = interp1d(self.tau_grid, self.kappa_dot_grid,
                                           kind='cubic', fill_value='extrapolate')
        self._visibility_of_tau = interp1d(self.tau_grid, self.visibility_grid,
                                            kind='cubic', fill_value='extrapolate')
        self._kappa_of_tau = interp1d(self.tau_grid, kappa,
                                       kind='cubic', fill_value='extrapolate')

        # Summary
        peak_idx = np.argmax(self.visibility_grid)
        z_peak = 1.0 / self.a_grid[peak_idx] - 1
        print(f"[Background] tau_0 = {self.tau_0:.1f} Mpc, "
              f"tau_rec = {self.tau_rec:.1f} Mpc, "
              f"D_A = {self.D_A:.1f} Mpc")
        print(f"[Background] r_s = {self.r_s:.1f} Mpc, "
              f"l_1 ~ {np.pi * self.D_A / self.r_s:.0f}")
        print(f"[Background] R_rec = {self.R_rec:.3f}, "
              f"k_D = {self.k_D:.4f} Mpc^-1")
        print(f"[Background] Visibility peak: z = {z_peak:.0f}")

    def _compute_recombination_peebles(self):
        """
        Peebles three-level atom recombination.

        Uses the ODE solver from recombination.py, then computes kappa, visibility
        on the same grids as the tanh method for full compatibility.
        """
        from .recombination import PeeblesRecombination

        solver = PeeblesRecombination(self, verbose=True)
        x_e_eff = solver.solve()

        if x_e_eff is None:
            # Peebles failed, fall back to tanh
            print("[Background] Peebles recombination failed, falling back to tanh.")
            self._compute_recombination()
            return

        self.x_e_grid = x_e_eff

        # Thomson scattering rate: same formula as tanh method
        rho_b0_SI = Omega_b * 3 * H0_SI**2 / (8 * np.pi * G_SI)
        n_H0 = (1 - Y_He) * rho_b0_SI / m_H_SI
        kappa_rate = x_e_eff * n_H0 * sigma_T_SI * Mpc_SI / self.a_grid**2

        # Optical depth
        kappa_cumul = cumulative_trapezoid(kappa_rate, self.tau_grid, initial=0.0)
        kappa = kappa_cumul[-1] - kappa_cumul

        self.kappa_grid = kappa
        self.kappa_dot_grid = -kappa_rate
        self.visibility_grid = kappa_rate * np.exp(-kappa)

        # Interpolators (same interface as tanh)
        self._kappa_dot_of_tau = interp1d(self.tau_grid, self.kappa_dot_grid,
                                           kind='cubic', fill_value='extrapolate')
        self._visibility_of_tau = interp1d(self.tau_grid, self.visibility_grid,
                                            kind='cubic', fill_value='extrapolate')
        self._kappa_of_tau = interp1d(self.tau_grid, kappa,
                                       kind='cubic', fill_value='extrapolate')

        # Summary
        peak_idx = np.argmax(self.visibility_grid)
        z_peak = 1.0 / self.a_grid[peak_idx] - 1
        print(f"[Background] tau_0 = {self.tau_0:.1f} Mpc, "
              f"tau_rec = {self.tau_rec:.1f} Mpc, "
              f"D_A = {self.D_A:.1f} Mpc")
        print(f"[Background] r_s = {self.r_s:.1f} Mpc, "
              f"l_1 ~ {np.pi * self.D_A / self.r_s:.0f}")
        print(f"[Background] R_rec = {self.R_rec:.3f}, "
              f"k_D = {self.k_D:.4f} Mpc^-1")
        print(f"[Background] Visibility peak (Peebles): z = {z_peak:.0f}")

    def _compute_recombination_recfast(self):
        """
        RECFAST-like recombination with helium, T_b evolution, and
        enhanced fudge factors.

        Uses the RecfastRecombination solver from recombination.py, then
        computes kappa, visibility on the same grids for full compatibility.
        """
        from .recombination import RecfastRecombination

        solver = RecfastRecombination(self, verbose=True)
        x_e_eff = solver.solve()

        if x_e_eff is None:
            print("[Background] RECFAST recombination failed, "
                  "falling back to Peebles.")
            self._compute_recombination_peebles()
            return

        self.x_e_grid = x_e_eff

        # Thomson scattering rate: same formula as other methods
        rho_b0_SI = Omega_b * 3 * H0_SI**2 / (8 * np.pi * G_SI)
        n_H0 = (1 - Y_He) * rho_b0_SI / m_H_SI
        kappa_rate = x_e_eff * n_H0 * sigma_T_SI * Mpc_SI / self.a_grid**2

        # Optical depth
        kappa_cumul = cumulative_trapezoid(kappa_rate, self.tau_grid, initial=0.0)
        kappa = kappa_cumul[-1] - kappa_cumul

        self.kappa_grid = kappa
        self.kappa_dot_grid = -kappa_rate
        self.visibility_grid = kappa_rate * np.exp(-kappa)

        # Interpolators (same interface)
        self._kappa_dot_of_tau = interp1d(self.tau_grid, self.kappa_dot_grid,
                                           kind='cubic', fill_value='extrapolate')
        self._visibility_of_tau = interp1d(self.tau_grid, self.visibility_grid,
                                            kind='cubic', fill_value='extrapolate')
        self._kappa_of_tau = interp1d(self.tau_grid, kappa,
                                       kind='cubic', fill_value='extrapolate')

        # Summary with visibility FWHM diagnostic
        peak_idx = np.argmax(self.visibility_grid)
        z_peak = 1.0 / self.a_grid[peak_idx] - 1
        g_peak = self.visibility_grid[peak_idx]

        # Compute FWHM in Mpc
        half_max = g_peak / 2.0
        above_half = self.visibility_grid > half_max
        if np.any(above_half):
            tau_above = self.tau_grid[above_half]
            fwhm = tau_above[-1] - tau_above[0]
        else:
            fwhm = 0.0

        # Integral of visibility
        g_integral = np.trapezoid(self.visibility_grid, self.tau_grid)

        print(f"[Background] tau_0 = {self.tau_0:.1f} Mpc, "
              f"tau_rec = {self.tau_rec:.1f} Mpc, "
              f"D_A = {self.D_A:.1f} Mpc")
        print(f"[Background] r_s = {self.r_s:.1f} Mpc, "
              f"l_1 ~ {np.pi * self.D_A / self.r_s:.0f}")
        print(f"[Background] R_rec = {self.R_rec:.3f}, "
              f"k_D = {self.k_D:.4f} Mpc^-1")
        print(f"[Background] Visibility peak (RECFAST): z = {z_peak:.0f}, "
              f"g_peak = {g_peak:.4f}")
        print(f"[Background] Visibility FWHM = {fwhm:.1f} Mpc, "
              f"integral = {g_integral:.4f}")

    # === Evaluation on arbitrary tau arrays (for MLX precomputation) ===

    def calH_at_tau(self, tau_arr):
        """Conformal Hubble at tau array. Returns numpy array."""
        return self._calH_of_tau(tau_arr)

    def R_at_tau(self, tau_arr):
        """Baryon loading R at tau array."""
        return self._R_of_tau(tau_arr)

    def kappa_dot_at_tau(self, tau_arr):
        """Thomson scattering rate at tau array."""
        return self._kappa_dot_of_tau(tau_arr)

    def visibility_at_tau(self, tau_arr):
        """Visibility function at tau array."""
        return self._visibility_of_tau(tau_arr)

    def kappa_at_tau(self, tau_arr):
        """Optical depth at tau array."""
        return self._kappa_of_tau(tau_arr)

    def a_at_tau(self, tau_arr):
        """Scale factor at tau array."""
        return self._a_of_tau(tau_arr)
