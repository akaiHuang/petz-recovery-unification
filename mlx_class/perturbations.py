"""
perturbations.py — Batched Boltzmann ODE solver on Apple MLX GPU.

All k-modes integrated simultaneously via batched RK4.
The stiff Poisson equation eigenvalue k^2/(3 calH) determines the step size.

Author: Sheng-Kai Huang, 2026
"""
import numpy as np
import mlx.core as mx
import time

from .background import Omega_r as _OMEGA_R, Omega_b as _OMEGA_B, Omega_c as _OMEGA_C
from .background import H0_Mpc as _H0_MPC

# TCA indices
TCA_PHI = 0
TCA_DELTA_B = 1
TCA_DELTA_C = 2
TCA_V_C = 3
TCA_THETA_0 = 4
TCA_THETA_1 = 5
TCA_N_VAR = 6


def _poisson_phi_dot(Phi, delta_gamma, delta_b, delta_c, k_arr, calH, a,
                     Omega_r, Omega_b, Omega_c, H0_Mpc):
    H02 = H0_Mpc * H0_Mpc
    density = (Omega_r / (a * a) * delta_gamma
               + Omega_b / a * delta_b
               + Omega_c / a * delta_c)
    return (-1.5 * H02 * density - k_arr * k_arr * Phi) / (3.0 * calH) - calH * Phi


def deriv_tca(y, k_arr, calH, R, Omega_r, Omega_b, Omega_c, a, H0_Mpc):
    Phi = y[:, TCA_PHI]
    delta_b = y[:, TCA_DELTA_B]
    delta_c = y[:, TCA_DELTA_C]
    v_c = y[:, TCA_V_C]
    Theta_0 = y[:, TCA_THETA_0]
    Theta_1 = y[:, TCA_THETA_1]

    Phi_dot = _poisson_phi_dot(Phi, 4.0*Theta_0, delta_b, delta_c,
                                k_arr, calH, a, Omega_r, Omega_b, Omega_c, H0_Mpc)
    dTheta_0 = -k_arr * Theta_1 - Phi_dot
    dTheta_1 = (k_arr / 3.0 * (Theta_0 + Phi) - calH * R * Theta_1) / (1.0 + R)
    d_delta_b = -3.0 * k_arr * Theta_1 - 3.0 * Phi_dot
    d_delta_c = -(k_arr * v_c + 3.0 * Phi_dot)
    d_v_c = -calH * v_c + k_arr * Phi

    return mx.stack([Phi_dot, d_delta_b, d_delta_c, d_v_c, dTheta_0, dTheta_1], axis=1)


def rk4_step(y, deriv_fn, dtau):
    k1 = deriv_fn(y)
    k2 = deriv_fn(y + 0.5 * dtau * k1)
    k3 = deriv_fn(y + 0.5 * dtau * k2)
    k4 = deriv_fn(y + dtau * k3)
    return y + (dtau / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


def build_tau_grid(bg, k_max, N_early=300):
    """
    Build TCA grid with step size determined by the Poisson stiffness.
    The stiff eigenvalue is k^2 / (3 calH). RK4 stability: dtau < 2.8 / eigenvalue.
    """
    tau_rec = bg.tau_rec
    tau_switch = 0.85 * tau_rec
    tau_early_end = min(10.0, 0.3 * tau_switch)
    tau_init = bg.tau_grid[1]

    # Early phase: geomspace
    tau_early = np.geomspace(tau_init, tau_early_end, N_early)

    # Late phase: step size limited by max eigenvalue k^2/(3 calH_min)
    # calH(tau) increases monotonically during radiation+matter era in this range
    # calH at tau_early_end:
    calH_at_end = float(bg.calH_at_tau(np.array([tau_early_end]))[0])
    calH_at_switch = float(bg.calH_at_tau(np.array([tau_switch]))[0])
    calH_min = min(calH_at_end, calH_at_switch)

    # Maximum eigenvalue
    max_eigenvalue = k_max**2 / (3.0 * calH_min)
    dtau_max = 2.5 / max_eigenvalue  # safety factor (2.5 < 2.8)

    N_late = max(200, int((tau_switch - tau_early_end) / dtau_max) + 10)
    tau_late = np.linspace(tau_early_end, tau_switch, N_late)

    tau_tca = np.unique(np.concatenate([tau_early, tau_late]))
    return tau_tca


def adiabatic_ic_tca(k_arr_np, bg):
    N_k = len(k_arr_np)
    tau_init = bg.tau_grid[1]
    y0 = np.zeros((N_k, TCA_N_VAR), dtype=np.float32)
    y0[:, TCA_PHI] = 1.0
    y0[:, TCA_DELTA_B] = -1.5
    y0[:, TCA_DELTA_C] = -1.5
    y0[:, TCA_V_C] = k_arr_np * tau_init / 6.0
    y0[:, TCA_THETA_0] = -0.5
    y0[:, TCA_THETA_1] = k_arr_np * tau_init / 18.0
    return mx.array(y0)


class BatchedBoltzmannSolver:
    def __init__(self, bg, k_arr_Mpc, l_max=12, mode='fast'):
        self.bg = bg
        self.k_arr_np = np.asarray(k_arr_Mpc, dtype=np.float32)
        self.N_k = len(self.k_arr_np)
        self.k_arr = mx.array(self.k_arr_np)
        k_max = float(self.k_arr_np[-1])
        self.tau_tca = build_tau_grid(bg, k_max)
        print(f"[Perturbations] N_k={self.N_k}, TCA: {len(self.tau_tca)} steps "
              f"[{self.tau_tca[0]:.2e}, {self.tau_tca[-1]:.1f}] Mpc")

    def solve(self):
        t0 = time.time()
        bg = self.bg
        k_arr = self.k_arr
        tau_grid = self.tau_tca

        y = adiabatic_ic_tca(self.k_arr_np, bg)
        calH = mx.array(bg.calH_at_tau(tau_grid).astype(np.float32))
        R = mx.array(bg.R_at_tau(tau_grid).astype(np.float32))
        a = mx.array(bg.a_at_tau(tau_grid).astype(np.float32))
        _Or, _Ob, _Oc, _H0 = _OMEGA_R, _OMEGA_B, float(bg.Omega_cdm), _H0_MPC

        print("[Perturbations] Integrating TCA...")
        for i in range(len(tau_grid) - 1):
            dt = float(tau_grid[i+1] - tau_grid[i])
            ci, ri, ai = calH[i], R[i], a[i]
            def rhs(y_in, _c=ci, _R=ri, _a=ai):
                return deriv_tca(y_in, k_arr, _c, _R, _Or, _Ob, _Oc, _a, _H0)
            y = rk4_step(y, rhs, dt)
            if i % 200 == 0:
                mx.eval(y)

        mx.eval(y)
        t_total = time.time() - t0
        print(f"[Perturbations] Done in {t_total:.2f}s ({len(tau_grid)} steps)")
        return PerturbationResult(y_tca=y, k_arr=self.k_arr_np, bg=bg)


class PerturbationResult:
    def __init__(self, y_tca, k_arr, bg):
        self.y_tca = y_tca
        self.k_arr = k_arr
        self.bg = bg
        self.N_k = len(k_arr)

    def source_at_recombination(self):
        """
        Extract source at TCA endpoint (0.85*tau_rec).
        Silk damping applied analytically.
        Note: evaluated at r_s(0.85*tau_rec) ~ 125 Mpc, not r_s(tau_rec) ~ 144 Mpc.
        This gives correct peak SPACING but ~15% offset in absolute positions.
        """
        y_np = np.array(self.y_tca)
        Theta_0 = y_np[:, TCA_THETA_0]
        Phi = y_np[:, TCA_PHI]
        v_b = 3.0 * y_np[:, TCA_THETA_1]

        silk = np.exp(-(self.k_arr / self.bg.k_D)**2)
        return Theta_0 * silk, Phi * silk, v_b * silk


class AnalyticTransfer:
    """Analytic transfer function (no ODE)."""
    def __init__(self, bg, khronon_correction=0.0):
        self.bg = bg
        self.khr_corr = khronon_correction
        from .background import k_eq
        self.k_eq = k_eq
        self.r_s = bg.r_s
        self.R_rec = bg.R_rec
        self.k_D = bg.k_D

    def compute_source(self, k_arr):
        """Compute monopole source S0 = Theta_0 + Phi (Sachs-Wolfe)."""
        kr_s = k_arr * self.r_s
        x_eq = k_arr / self.k_eq
        baryon_damp = (1.0 + self.R_rec)**(-0.25)
        zero_shift = self.R_rec / (3.0 * (1.0 + self.R_rec))
        driving = 1.0 + 2.5 * x_eq**2 / (1.0 + x_eq**2)
        A_osc = (1.0 / 3.0) * driving * baryon_damp
        silk = np.exp(-(k_arr / self.k_D)**2)
        source = (A_osc * np.cos(kr_s) + zero_shift) * silk
        if self.khr_corr != 0.0:
            source *= (1.0 + self.khr_corr)
        return source

    def compute_doppler(self, k_arr):
        """
        Compute Doppler source S1 = v_b from WKB approximation.
        v_b = 3 Theta_1 = -3 c_s A_osc sin(k r_s) / (1+R) * silk
        The Doppler term shifts acoustic peaks by ~10-15% in l.
        """
        kr_s = k_arr * self.r_s
        x_eq = k_arr / self.k_eq
        baryon_damp = (1.0 + self.R_rec)**(-0.25)
        driving = 1.0 + 2.5 * x_eq**2 / (1.0 + x_eq**2)
        A_osc = (1.0 / 3.0) * driving * baryon_damp

        cs_rec = 1.0 / np.sqrt(3.0 * (1.0 + self.R_rec))
        silk = np.exp(-(k_arr / self.k_D)**2)

        # v_b = -3 cs/(1+R) * A_osc * sin(kr_s) * silk
        v_b = -3.0 * cs_rec / (1.0 + self.R_rec) * A_osc * np.sin(kr_s) * silk
        if self.khr_corr != 0.0:
            v_b *= (1.0 + self.khr_corr)
        return v_b
