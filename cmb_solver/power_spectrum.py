"""
power_spectrum.py — CMB TT angular power spectrum C_l calculator.

Uses line-of-sight integration (Seljak & Zaldarriaga 1996):
    C_l = (4pi) int dk/k P_R(k) |Delta_l(k)|^2
    Delta_l(k) = int S(k,tau) j_l[k(tau_0-tau)] dtau

Source function (dominant terms):
    S = g * (Theta_0 + Psi)  +  g * v_b / k  +  exp(-kappa) * (Phi' + Psi')
      = Sachs-Wolfe         +   Doppler     +   ISW
"""
import numpy as np
from scipy.special import spherical_jn
from scipy.integrate import trapezoid
from . import constants as C
from . import perturbations as PT
from . import boltzmann as BZ


class PowerSpectrumCalculator:
    """Compute CMB TT power spectrum D_l = l(l+1)C_l/(2pi)."""

    def __init__(self, bg, l_max_hierarchy=BZ.L_MAX_DEFAULT):
        self.bg = bg
        self.l_max_hierarchy = l_max_hierarchy

    def primordial_spectrum(self, k_phys):
        """Dimensionless primordial power Delta^2_R(k) where k in 1/Mpc."""
        return C.A_s * (k_phys / C.k_pivot) ** (C.n_s - 1)

    def _solve_k_mode(self, k, tau_grid):
        """
        Solve perturbation equations for wavenumber k and evaluate on tau_grid.
        Returns Phi, Theta_0, v_b, Phi_dot arrays.
        Phi_dot computed from finite differences of the ODE solution.
        """
        bg = self.bg
        pert = PT.PerturbationSolver(k, bg, l_max=self.l_max_hierarchy)
        sol = pert.solve(rtol=1e-4, atol=1e-6)

        n = len(tau_grid)
        Phi = np.zeros(n)
        Theta_0 = np.zeros(n)
        v_b = np.zeros(n)

        tau_min = sol.tau_init
        tau_max = sol.tau_final

        for i, tau in enumerate(tau_grid):
            tau_c = np.clip(tau, tau_min, tau_max)
            y = sol.evaluate(tau_c)
            Phi[i] = y[PT.IDX_PHI]
            Theta_0[i] = y[PT.IDX_THETA_START]
            v_b[i] = y[PT.IDX_V_B]

        # Phi_dot from finite differences (reflects actual ODE evolution)
        Phi_dot = np.gradient(Phi, tau_grid)

        return Phi, Theta_0, v_b, Phi_dot, sol.success

    def compute_source(self, k, tau_grid, Phi, Theta_0, v_b, Phi_dot):
        """Compute source function S(k, tau) on the tau grid."""
        bg = self.bg
        n = len(tau_grid)
        source = np.zeros(n)

        for i, tau in enumerate(tau_grid):
            g_vis = bg.visibility_interp(tau)
            kappa = bg.kappa_interp(tau)

            # SW
            sw = g_vis * (Theta_0[i] + Phi[i])  # Psi = Phi

            # Doppler
            doppler = g_vis * v_b[i] / max(k, 1e-10) if k > 0 else 0.0

            # ISW (Phi' + Psi' = 2 Phi')
            isw = np.exp(-kappa) * 2.0 * Phi_dot[i]

            source[i] = sw + doppler + isw

        return source

    def compute_transfer(self, k, ell_values, n_tau=400):
        """
        Compute Delta_l(k) via line-of-sight integration.
        """
        bg = self.bg
        tau_0 = bg.tau_grid[-1]
        tau_rec = bg.tau_of_a(C.a_rec)

        # Tau grid: concentrated around recombination + some late-time for ISW
        tau_min = max(bg.tau_grid[1], tau_rec * 0.5)
        n_rec = int(0.8 * n_tau)
        n_late = n_tau - n_rec

        tau_rec_width = 0.2 * tau_rec
        tau_lo = max(tau_min, tau_rec - 4 * tau_rec_width)
        tau_hi = min(tau_0, tau_rec + 4 * tau_rec_width)

        tau_near = np.linspace(tau_lo, tau_hi, n_rec)
        tau_late = np.linspace(tau_hi * 1.01, tau_0, n_late + 1)[1:]
        tau_grid = np.sort(np.unique(np.concatenate([tau_near, tau_late])))

        # Solve perturbation
        Phi, Theta_0, v_b, Phi_dot, success = self._solve_k_mode(k, tau_grid)
        if not success:
            return {ell: 0.0 for ell in ell_values}

        # Source function
        source = self.compute_source(k, tau_grid, Phi, Theta_0, v_b, Phi_dot)

        # Line-of-sight integration
        chi_grid = tau_0 - tau_grid
        Delta_l = {}
        for ell in ell_values:
            x = k * chi_grid
            jl = spherical_jn(int(ell), x)
            Delta_l[ell] = trapezoid(source * jl, tau_grid)

        return Delta_l

    def compute_cl(self, ell_values, k_grid=None, n_k=60, n_tau=400,
                   progress=True):
        """
        Compute C_l for given multipoles.

        C_l = (4pi) int dk/k P_R(k) |Delta_l(k)|^2
        """
        ell_values = np.array(ell_values, dtype=int)

        if k_grid is None:
            tau_0 = self.bg.tau_grid[-1]
            tau_rec = self.bg.tau_of_a(C.a_rec)
            chi_rec = tau_0 - tau_rec

            k_min = max(1.0 / tau_0, 0.3)
            k_max = float(max(ell_values)) * 2.0 / chi_rec
            k_max = max(k_max, 10.0 * k_min)

            k_grid = np.geomspace(k_min, k_max, n_k)

        if progress:
            print(f"[C_l] {len(ell_values)} multipoles, {len(k_grid)} k-modes")
            print(f"[C_l] k range: [{k_grid[0]:.2f}, {k_grid[-1]:.1f}] (H0 units)")

        # Compute Delta_l(k) for each k
        Delta_l_arr = np.zeros((len(k_grid), len(ell_values)))

        for i_k, k in enumerate(k_grid):
            if progress and (i_k % max(1, len(k_grid) // 10) == 0):
                print(f"  [{100*(i_k+1)/len(k_grid):5.1f}%] k = {k:.2f}")

            Delta_l = self.compute_transfer(k, ell_values, n_tau=n_tau)
            for i_ell, ell in enumerate(ell_values):
                Delta_l_arr[i_k, i_ell] = Delta_l[ell]

        # k in H0 units -> physical k in 1/Mpc
        k_phys = k_grid * C.H0_invMpc

        # Integrate: C_l = 4pi int dk/k P_R(k) |Delta_l(k)|^2
        Cl = np.zeros(len(ell_values))
        for i_ell in range(len(ell_values)):
            P_R = self.primordial_spectrum(k_phys)
            integrand = P_R * Delta_l_arr[:, i_ell]**2 / k_grid
            Cl[i_ell] = 4.0 * np.pi * trapezoid(integrand, k_grid)

        # D_l in muK^2
        T_muK = C.T_CMB * 1e6
        Dl = ell_values * (ell_values + 1) * Cl / (2.0 * np.pi) * T_muK**2

        return ell_values, Cl, Dl


def compute_cl_fast(bg, ell_values, n_k=50, n_tau=300, l_max_hierarchy=8,
                    progress=True):
    """Convenience wrapper."""
    calc = PowerSpectrumCalculator(bg, l_max_hierarchy=l_max_hierarchy)
    return calc.compute_cl(ell_values, n_k=n_k, n_tau=n_tau, progress=progress)
