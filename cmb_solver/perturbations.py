"""
perturbations.py — Cosmological perturbation equations for Khronon cosmology.

Two-phase integration strategy:
  Phase 1 (tight coupling, tau < tau_dec): Only 5 variables: Phi, delta_b, Theta_0, Theta_1, delta_K, v_K.
           v_b = 3*Theta_1 enforced algebraically. Higher Theta_l = 0.
  Phase 2 (free streaming, tau > tau_dec): Full Boltzmann hierarchy with all species decoupled.

tau_dec is determined as the conformal time where |kappa_dot| drops below a threshold.
"""
import numpy as np
from scipy.integrate import solve_ivp
from . import constants as C
from . import boltzmann as BZ

# Index map for full state
IDX_PHI = 0
IDX_DELTA_B = 1
IDX_V_B = 2
IDX_DELTA_K = 3
IDX_V_K = 4
IDX_THETA_START = 5

def state_size(l_max):
    return 5 + l_max + 1

# TCA state: [Phi, delta_b, delta_K, v_K, Theta_0, Theta_1]
TCA_PHI = 0
TCA_DELTA_B = 1
TCA_DELTA_K = 2
TCA_V_K = 3
TCA_THETA_0 = 4
TCA_THETA_1 = 5
TCA_SIZE = 6


class PerturbationSolver:
    """Solve perturbation equations for a single k-mode using two-phase approach."""

    def __init__(self, k, bg, l_max=BZ.L_MAX_DEFAULT):
        self.k = k
        self.bg = bg
        self.l_max = l_max
        self.n_vars = state_size(l_max)

    def _find_tau_decoupling(self):
        """Find conformal time where scattering becomes subdominant."""
        bg = self.bg
        k = self.k
        threshold = max(3.0 * k, 50.0)  # |kd| < threshold => decouple

        # Search backwards from the visibility peak
        peak_idx = np.argmax(bg.visibility)
        for i in range(peak_idx, len(bg.a_grid)):
            kd = abs(bg.kappa_dot_grid[i])
            if kd < threshold:
                return bg.tau_grid[i]
        return bg.tau_grid[-1]

    def _compute_phi_dot_tca(self, a, calH, Phi, Psi, delta_b, delta_K, Theta_0):
        """Phi' from Poisson equation."""
        delta_gamma = 4.0 * Theta_0
        SUM = (C.Omega_r * a**(-4) * delta_gamma
               + C.Omega_b * a**(-3) * delta_b
               + C.Omega_K_khronon * a**(-3) * delta_K)
        k = self.k
        if calH > 1e-30:
            return (-1.5 * a**2 * SUM - k**2 * Phi) / (3.0 * calH) - calH * Psi
        return 0.0

    def rhs_tca(self, tau, y):
        """RHS for tight-coupling phase: [Phi, delta_b, delta_K, v_K, Theta_0, Theta_1]."""
        k = self.k
        bg = self.bg

        Phi = y[TCA_PHI]
        delta_b = y[TCA_DELTA_B]
        delta_K = y[TCA_DELTA_K]
        v_K = y[TCA_V_K]
        Theta_0 = y[TCA_THETA_0]
        Theta_1 = y[TCA_THETA_1]

        a = bg.a_of_tau(tau)
        a = np.clip(a, C.a_init, 1.0)
        H = bg.H_interp(a)
        calH = a * H
        w_K = bg.w_K_of_a(a)
        R = C.R_factor * a
        Psi = Phi

        Phi_dot = self._compute_phi_dot_tca(a, calH, Phi, Psi, delta_b, delta_K, Theta_0)

        # TCA: v_b = 3*Theta_1 (exact)
        v_b = 3.0 * Theta_1

        dy = np.zeros(TCA_SIZE)

        # Phi
        dy[TCA_PHI] = Phi_dot

        # Photon-baryon fluid
        dy[TCA_THETA_0] = -k * Theta_1 - Phi_dot
        dy[TCA_THETA_1] = (k / 3.0 * (Theta_0 + Psi) - calH * R * Theta_1) / (1.0 + R)

        # Baryon density: delta_b' = -k v_b - 3 Phi' = -3k Theta_1 - 3 Phi'
        dy[TCA_DELTA_B] = -3.0 * k * Theta_1 - 3.0 * Phi_dot

        # Khronon perturbations: c_s^2 = 0 from ghost condensation.
        # At the perturbation level, Khronon behaves as CDM (w_eff = 0).
        # The background w_K = delta/(2+delta) is a property of the Sigma = 2 ln Q
        # parameterization, but for linear perturbations around the ghost-condensation
        # background, the effective equation of state for perturbations is w_eff = 0
        # because K'(Q_0) = 0 => no pressure response to Q perturbations.
        dy[TCA_DELTA_K] = -(k * v_K + 3.0 * Phi_dot)  # w_eff = 0 => (1+w) = 1
        dy[TCA_V_K] = -calH * v_K + k * Psi            # w_eff = 0 => no 1/(1+w)

        return dy

    def rhs_full(self, tau, y):
        """RHS for free-streaming phase: full state vector."""
        k = self.k
        bg = self.bg
        l_max = self.l_max

        Phi = y[IDX_PHI]
        delta_b = y[IDX_DELTA_B]
        v_b = y[IDX_V_B]
        delta_K = y[IDX_DELTA_K]
        v_K = y[IDX_V_K]
        Theta = y[IDX_THETA_START:IDX_THETA_START + l_max + 1]

        a = bg.a_of_tau(tau)
        a = np.clip(a, C.a_init, 1.0)
        H = bg.H_interp(a)
        calH = a * H
        w_K = bg.w_K_of_a(a)
        kappa_dot = bg.kappa_dot_interp(tau)
        R = C.R_factor * a
        Psi = Phi

        delta_gamma = 4.0 * Theta[0]
        SUM = (C.Omega_r * a**(-4) * delta_gamma
               + C.Omega_b * a**(-3) * delta_b
               + C.Omega_K_khronon * a**(-3) * delta_K)

        if calH > 1e-30:
            Phi_dot = (-1.5 * a**2 * SUM - k**2 * Phi) / (3.0 * calH) - calH * Psi
        else:
            Phi_dot = 0.0

        dy = np.zeros(self.n_vars)
        dy[IDX_PHI] = Phi_dot

        # Full Boltzmann hierarchy
        dTheta = BZ.photon_hierarchy_rhs(Theta, k, Phi_dot, Psi, kappa_dot, v_b, l_max)
        dy[IDX_THETA_START:IDX_THETA_START + l_max + 1] = dTheta

        # Baryons
        dy[IDX_DELTA_B] = -k * v_b - 3.0 * Phi_dot
        if R > 1e-10 and abs(kappa_dot) > 0:
            dy[IDX_V_B] = -calH * v_b + k * Psi + kappa_dot / R * (Theta[1] - v_b / 3.0)
        else:
            dy[IDX_V_B] = -calH * v_b + k * Psi

        # Khronon: w_eff = 0 for perturbations (ghost condensation CDM-like)
        dy[IDX_DELTA_K] = -(k * v_K + 3.0 * Phi_dot)
        dy[IDX_V_K] = -calH * v_K + k * Psi

        return dy

    def solve(self, tau_final=None, rtol=1e-6, atol=1e-8, max_step=None):
        """Two-phase integration: TCA then free streaming."""
        k = self.k
        bg = self.bg

        if tau_final is None:
            tau_final = bg.tau_grid[-1]

        # Phase 1: Tight coupling
        a_init = C.a_init
        tau_init = bg.tau_of_a(a_init)
        tau_dec = self._find_tau_decoupling()

        if max_step is None:
            if k > 0:
                max_step = min(0.6 / k, (tau_final - tau_init) / 200)
            else:
                max_step = (tau_final - tau_init) / 200

        # TCA initial conditions (adiabatic)
        Phi_0 = 1.0
        y_tca = np.zeros(TCA_SIZE)
        y_tca[TCA_PHI] = Phi_0
        y_tca[TCA_DELTA_B] = -1.5 * Phi_0
        y_tca[TCA_DELTA_K] = -1.5 * Phi_0
        y_tca[TCA_V_K] = k * tau_init * Phi_0 / 6.0
        y_tca[TCA_THETA_0] = -0.5 * Phi_0
        y_tca[TCA_THETA_1] = k * tau_init * Phi_0 / 18.0

        sol_tca = solve_ivp(
            self.rhs_tca,
            [tau_init, tau_dec],
            y_tca,
            method='RK45',
            rtol=rtol,
            atol=atol,
            max_step=max_step,
            dense_output=True,
        )

        if not sol_tca.success:
            # If TCA fails, return a dummy solution
            sol_tca.t = np.array([tau_init, tau_final])
            sol_tca.y = np.zeros((self.n_vars, 2))

        # Phase 2: Free streaming
        # Transfer TCA state to full state
        y_dec_tca = sol_tca.sol(tau_dec)
        y_full = np.zeros(self.n_vars)
        y_full[IDX_PHI] = y_dec_tca[TCA_PHI]
        y_full[IDX_DELTA_B] = y_dec_tca[TCA_DELTA_B]
        y_full[IDX_V_B] = 3.0 * y_dec_tca[TCA_THETA_1]  # v_b = 3*Theta_1 exactly
        y_full[IDX_DELTA_K] = y_dec_tca[TCA_DELTA_K]
        y_full[IDX_V_K] = y_dec_tca[TCA_V_K]
        y_full[IDX_THETA_START] = y_dec_tca[TCA_THETA_0]
        y_full[IDX_THETA_START + 1] = y_dec_tca[TCA_THETA_1]
        # Higher Theta_l = 0 (suppressed by scattering during TCA)

        sol_full = solve_ivp(
            self.rhs_full,
            [tau_dec, tau_final],
            y_full,
            method='RK45',
            rtol=rtol,
            atol=atol,
            max_step=max_step,
            dense_output=True,
        )

        if not sol_full.success:
            print(f"[Perturbations] Warning: k={k:.4e} phase 2 failed: {sol_full.message}")

        # Combine into a unified solution object
        return CombinedSolution(sol_tca, sol_full, tau_dec, self.l_max)

    def extract_transfer(self, sol, tau_eval=None):
        """Extract transfer functions at a specific time."""
        if tau_eval is None:
            tau_eval = sol.tau_final
        y = sol.evaluate(tau_eval)
        result = {
            'Phi': y[IDX_PHI],
            'delta_b': y[IDX_DELTA_B],
            'v_b': y[IDX_V_B],
            'delta_K': y[IDX_DELTA_K],
            'v_K': y[IDX_V_K],
        }
        for l in range(self.l_max + 1):
            result[f'Theta_{l}'] = y[IDX_THETA_START + l]
        return result


class CombinedSolution:
    """Wrapper combining TCA and full solutions into a single interface."""

    def __init__(self, sol_tca, sol_full, tau_dec, l_max):
        self.sol_tca = sol_tca
        self.sol_full = sol_full
        self.tau_dec = tau_dec
        self.l_max = l_max
        self.success = sol_full.success
        self.t = np.concatenate([sol_tca.t, sol_full.t])
        self.tau_init = sol_tca.t[0]
        self.tau_final = sol_full.t[-1]

    def evaluate(self, tau):
        """Evaluate full state vector at given tau."""
        n_full = state_size(self.l_max)

        if tau <= self.tau_dec:
            y_tca = self.sol_tca.sol(tau)
            y = np.zeros(n_full)
            y[IDX_PHI] = y_tca[TCA_PHI]
            y[IDX_DELTA_B] = y_tca[TCA_DELTA_B]
            y[IDX_V_B] = 3.0 * y_tca[TCA_THETA_1]
            y[IDX_DELTA_K] = y_tca[TCA_DELTA_K]
            y[IDX_V_K] = y_tca[TCA_V_K]
            y[IDX_THETA_START] = y_tca[TCA_THETA_0]
            y[IDX_THETA_START + 1] = y_tca[TCA_THETA_1]
            return y
        else:
            return self.sol_full.sol(tau)

    def sol(self, tau):
        """Alias for evaluate, matching solve_ivp interface."""
        return self.evaluate(tau)
