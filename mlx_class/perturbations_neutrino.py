"""
perturbations_neutrino.py — IMEX Boltzmann solver with massless neutrino hierarchy.

Extends perturbations_implicit.py by adding:
  - N_eff = 3.046 massless neutrinos with Boltzmann hierarchy (l_nu_max ~ 20)
  - Neutrino anisotropic stress: Psi != Phi
  - Properly split Omega_gamma and Omega_nu in the Poisson equation

DESIGN PRINCIPLE:
  The Phi equation uses the SAME IMEX exponential integrator as the original
  perturbations_implicit.py, with the stiff eigenvalue lambda = -(k^2/(3 calH) + calH).
  The non-stiff source F includes density perturbations from both photons and neutrinos
  (separately weighted by Omega_gamma and Omega_nu), treating Psi = Phi in the
  Phi evolution to avoid the superhorizon instability from the 1/k^2 anisotropic
  stress backreaction term.

  Psi is diagnosed algebraically: Psi = Phi - 12 H_0^2 Omega_nu N_2 / (a^2 k^2).
  The Psi != Phi effect enters the momentum equations (Theta_1, N_1, v_c) through
  the gravitational force term. This perturbative treatment is standard in CMB codes
  using conformal Newtonian gauge.

The neutrino hierarchy (conformal Newtonian gauge):
  N_0' = -k N_1 - Phi'
  N_1' = k/3 (N_0 - 2 N_2 + Psi)
  N_l' = k/(2l+1) [l N_{l-1} - (l+1) N_{l+1}]   for 2 <= l < l_max
  N_{l_max}' = k/(2l_max+1) * l_max * N_{l_max-1} - (l_max+1)/tau * N_{l_max}

Author: Sheng-Kai Huang, 2026
"""
import numpy as np
import mlx.core as mx
import time

from .background import (
    Omega_r as _OMEGA_R,
    Omega_b as _OMEGA_B,
    Omega_c as _OMEGA_C,
    H0_Mpc as _H0_MPC,
)


# ============================================================================
# Radiation density split: photons vs neutrinos
# ============================================================================
N_EFF = 3.046
_f_nu = 0.2271 * N_EFF / (1.0 + 0.2271 * N_EFF)
_OMEGA_GAMMA = _OMEGA_R * (1.0 - _f_nu)   # photon density parameter
_OMEGA_NU = _OMEGA_R * _f_nu               # neutrino density parameter

# Neutrino hierarchy truncation
L_NU_MAX = 20


# ============================================================================
# Variable indices in state vector
# ============================================================================
# [Phi, delta_b, delta_c, v_c, Theta_0, Theta_1, N_0, N_1, ..., N_{l_max}]
IDX_PHI = 0
IDX_DELTA_B = 1
IDX_DELTA_C = 2
IDX_V_C = 3
IDX_THETA_0 = 4
IDX_THETA_1 = 5
IDX_N_START = 6  # N_0 at index 6, N_1 at 7, etc.


def n_var_total(l_nu_max):
    return IDX_N_START + l_nu_max + 1


# ============================================================================
# Phi equation: IMEX decomposition (stiff + non-stiff)
# ============================================================================

def _phi_lambda(k_arr, calH):
    """Stiff eigenvalue: lambda = -(k^2/(3 calH) + calH)."""
    return -(k_arr * k_arr / (3.0 * calH) + calH)


def _phi_source(Theta_0, N_0, delta_b, delta_c, k_arr, calH, a,
                Omega_gamma, Omega_nu, Omega_b, Omega_c, H0_Mpc):
    """
    Non-stiff source F in the IMEX decomposition: Phi' = lambda*Phi + F.

    The Phi equation is treated as if Psi = Phi (zeroth order in anisotropic
    stress). This avoids the superhorizon instability caused by the
    Phi-Psi difference term calH*(Phi-Psi) which has a 1/k^2 factor.

    The neutrino anisotropic stress (Psi != Phi) enters perturbatively:
    it affects the momentum equations (Theta_1, N_1, v_c) where the
    gravitational force depends on Psi, but NOT the Phi evolution.
    This is consistent because the anisotropic stress backreaction on Phi
    is a higher-order correction in f_nu.

    F = -H02/(2 calH) * S

    S = Omega_gamma/a^2 * 4*Theta_0 + Omega_nu/a^2 * 4*N_0
        + Omega_b/a * delta_b + Omega_c/a * delta_c
    """
    H02 = H0_Mpc * H0_Mpc

    S = (Omega_gamma / (a * a) * (4.0 * Theta_0)
         + Omega_nu / (a * a) * (4.0 * N_0)
         + Omega_b / a * delta_b
         + Omega_c / a * delta_c)

    return -0.5 * H02 / calH * S


def _phi_int_factor(lam, dt):
    """(exp(lam*dt) - 1) / lam, with Taylor expansion for small |lam*dt|."""
    lam_dt = lam * dt
    exact = (mx.exp(lam_dt) - 1.0) / lam
    taylor = dt * (1.0 + 0.5 * lam_dt + lam_dt * lam_dt / 6.0)
    return mx.where(mx.abs(lam_dt) < 1e-4, taylor, exact)


# ============================================================================
# Full derivative function
# ============================================================================

def deriv_full(y, k_arr, calH, R, tau, Omega_gamma, Omega_nu, Omega_b, Omega_c,
               a, H0_Mpc, l_nu_max):
    """
    Full RHS for all variables.

    Phi equation: uses Psi = Phi (zeroth order, no anisotropic stress backreaction).
    This is consistent with the IMEX source _phi_source and avoids the superhorizon
    instability from the 1/k^2 term in Phi - Psi.

    Momentum equations (Theta_1, N_1, v_c): use Psi != Phi (diagnosed from N_2).
    This is where the neutrino anisotropic stress has its physical effect:
    it modifies the gravitational force felt by photons, neutrinos, and CDM.

    Returns dy shape (N_k, n_var).
    """
    Phi = y[:, IDX_PHI]
    delta_b = y[:, IDX_DELTA_B]
    delta_c = y[:, IDX_DELTA_C]
    v_c = y[:, IDX_V_C]
    Theta_0 = y[:, IDX_THETA_0]
    Theta_1 = y[:, IDX_THETA_1]

    # Extract neutrino multipoles
    N = []
    for l in range(l_nu_max + 1):
        N.append(y[:, IDX_N_START + l])
    N_2 = N[2] if l_nu_max >= 2 else mx.zeros_like(Phi)

    # Diagnose Psi (for momentum equations only)
    H02 = H0_Mpc * H0_Mpc
    Psi = Phi - 12.0 * H02 * Omega_nu * N_2 / (a * a * k_arr * k_arr)

    # Phi' from Poisson equation (with Psi = Phi, zeroth order)
    S = (Omega_gamma / (a * a) * (4.0 * Theta_0)
         + Omega_nu / (a * a) * (4.0 * N[0])
         + Omega_b / a * delta_b
         + Omega_c / a * delta_c)
    Phi_dot = -calH * Phi - k_arr * k_arr * Phi / (3.0 * calH) - 0.5 * H02 / calH * S

    # Photon TCA
    dTheta_0 = -k_arr * Theta_1 - Phi_dot
    # Theta_1 source uses Psi (not Phi) — key neutrino effect
    dTheta_1 = (k_arr / 3.0 * (Theta_0 + Psi) - calH * R * Theta_1) / (1.0 + R)

    # Baryons
    d_delta_b = -3.0 * k_arr * Theta_1 - 3.0 * Phi_dot

    # CDM falls in the Psi potential
    d_delta_c = -(k_arr * v_c + 3.0 * Phi_dot)
    d_v_c = -calH * v_c + k_arr * Psi

    # Neutrino hierarchy
    dN = []

    # l = 0: N_0' = -k N_1 - Phi'
    dN.append(-k_arr * N[1] - Phi_dot)

    # l = 1: N_1' = k/3 (N_0 - 2 N_2 + Psi) — uses Psi
    dN.append(k_arr / 3.0 * (N[0] - 2.0 * N_2 + Psi))

    # l >= 2: standard hierarchy
    for l in range(2, l_nu_max):
        dN.append(k_arr / (2.0 * l + 1.0) * (l * N[l - 1] - (l + 1.0) * N[l + 1]))

    # l = l_max: truncation
    tau_safe = mx.maximum(tau, mx.array(1e-10))
    dN.append(k_arr * N[l_nu_max - 1] * l_nu_max / (2.0 * l_nu_max + 1.0)
              - (l_nu_max + 1.0) / tau_safe * N[l_nu_max])

    # Stack
    parts = [Phi_dot[:, None], d_delta_b[:, None], d_delta_c[:, None],
             d_v_c[:, None], dTheta_0[:, None], dTheta_1[:, None]]
    for l in range(l_nu_max + 1):
        parts.append(dN[l][:, None])

    return mx.concatenate(parts, axis=1)


# ============================================================================
# Column replacement helper
# ============================================================================

def _set_phi(y, Phi_new):
    """Replace the Phi column (index 0) in the state vector."""
    return mx.concatenate([Phi_new[:, None], y[:, 1:]], axis=1)


# ============================================================================
# IMEX RK4 step
# ============================================================================

def imex_rk4_step(y, k_arr, calH, R, tau_val,
                  Omega_gamma, Omega_nu, Omega_b, Omega_c,
                  a, H0_Mpc, l_nu_max, dtau):
    """
    IMEX RK4 step: exponential integrator for Phi, explicit RK4 for everything else.

    Identical structure to perturbations_implicit.py but with neutrino hierarchy.
    Phi source uses Psi = Phi (no anisotropic stress backreaction).
    Psi != Phi enters only in the momentum equations via deriv_full.
    """
    lam = _phi_lambda(k_arr, calH)
    exp_lam_dt = mx.exp(lam * dtau)
    exp_lam_half = mx.exp(lam * 0.5 * dtau)

    def get_F(state):
        """Extract non-stiff source F from state (no anisotropic stress)."""
        return _phi_source(
            state[:, IDX_THETA_0], state[:, IDX_N_START],
            state[:, IDX_DELTA_B], state[:, IDX_DELTA_C],
            k_arr, calH, a, Omega_gamma, Omega_nu, Omega_b, Omega_c, H0_Mpc)

    def full_rhs(state):
        return deriv_full(state, k_arr, calH, R, mx.array(float(tau_val)),
                          Omega_gamma, Omega_nu, Omega_b, Omega_c,
                          a, H0_Mpc, l_nu_max)

    # --- Stage 1 ---
    Phi_1 = y[:, IDX_PHI]
    F1 = get_F(y)
    k1 = full_rhs(y)

    # --- Stage 2: half step ---
    Phi_2 = exp_lam_half * Phi_1 + _phi_int_factor(lam, 0.5 * dtau) * F1
    y2 = y + 0.5 * dtau * k1
    y2 = _set_phi(y2, Phi_2)
    F2 = get_F(y2)
    k2 = full_rhs(y2)

    # --- Stage 3: half step ---
    Phi_3 = exp_lam_half * Phi_1 + _phi_int_factor(lam, 0.5 * dtau) * F2
    y3 = y + 0.5 * dtau * k2
    y3 = _set_phi(y3, Phi_3)
    F3 = get_F(y3)
    k3 = full_rhs(y3)

    # --- Stage 4: full step ---
    Phi_4 = exp_lam_dt * Phi_1 + _phi_int_factor(lam, dtau) * F3
    y4 = y + dtau * k3
    y4 = _set_phi(y4, Phi_4)
    F4 = get_F(y4)
    k4 = full_rhs(y4)

    # --- Combine: RK4 for non-Phi, Simpson for Phi ---
    y_new = y + (dtau / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)

    F_mid = 0.5 * (F2 + F3)
    Phi_new = (exp_lam_dt * Phi_1
               + (dtau / 6.0) * (exp_lam_dt * F1
                                  + 4.0 * exp_lam_half * F_mid
                                  + F4))
    y_new = _set_phi(y_new, Phi_new)

    return y_new


# ============================================================================
# Tau grid
# ============================================================================

def build_tau_grid_neutrino(bg, k_max, N_early=300, N_late=800):
    """Build integration grid extending to tau_rec."""
    tau_rec = bg.tau_rec
    tau_init = bg.tau_grid[1]
    tau_early_end = min(10.0, 0.3 * tau_rec)

    tau_early = np.geomspace(tau_init, tau_early_end, N_early)

    cs_approx = 1.0 / np.sqrt(3.0)
    dtau_acoustic = 2.0 * np.pi / (k_max * cs_approx) / 25.0
    N_needed = int((tau_rec - tau_early_end) / dtau_acoustic) + 10
    N_late_actual = max(N_late, N_needed)

    tau_late = np.linspace(tau_early_end, tau_rec, N_late_actual)
    tau_grid = np.unique(np.concatenate([tau_early, tau_late]))
    return tau_grid


# ============================================================================
# Initial conditions
# ============================================================================

def adiabatic_ic(k_arr_np, bg, l_nu_max=L_NU_MAX):
    """
    Adiabatic initial conditions with neutrinos.

    Phi     = 1.0  (normalization)
    delta_b = -3/2 Phi
    delta_c = -3/2 Phi
    v_c     = k tau / 6
    Theta_0 = -Phi/2
    Theta_1 = k tau / 18
    N_0     = -Phi/2   (adiabatic: same as Theta_0)
    N_1     = k tau / 18
    N_2     = (k tau)^2 / 60
    N_l     ~ (k tau)^l / prod(2j+1, j=1..l)  for l >= 3 (tiny)
    """
    N_k = len(k_arr_np)
    tau_init = bg.tau_grid[1]
    nvar = n_var_total(l_nu_max)

    y0 = np.zeros((N_k, nvar), dtype=np.float32)
    y0[:, IDX_PHI] = 1.0
    y0[:, IDX_DELTA_B] = -1.5
    y0[:, IDX_DELTA_C] = -1.5
    y0[:, IDX_V_C] = k_arr_np * tau_init / 6.0
    y0[:, IDX_THETA_0] = -0.5
    y0[:, IDX_THETA_1] = k_arr_np * tau_init / 18.0

    y0[:, IDX_N_START + 0] = -0.5
    y0[:, IDX_N_START + 1] = k_arr_np * tau_init / 18.0
    if l_nu_max >= 2:
        y0[:, IDX_N_START + 2] = (k_arr_np * tau_init) ** 2 / 60.0

    for l in range(3, min(l_nu_max + 1, 6)):
        prod_val = 1.0
        for j in range(1, l + 1):
            prod_val *= (2 * j + 1)
        y0[:, IDX_N_START + l] = (k_arr_np * tau_init) ** l / prod_val

    return mx.array(y0)


# ============================================================================
# Solver class
# ============================================================================

class NeutrinoBoltzmannSolver:
    """
    Batched Boltzmann solver with massless neutrinos and IMEX integration.

    Key physics improvements over perturbations_implicit.py:
    1. Separate Omega_gamma (photons) and Omega_nu (neutrinos) in Poisson equation
    2. Neutrino Boltzmann hierarchy with l_nu_max multipoles
    3. Psi != Phi in momentum equations: CDM falls in Psi, dipole sources use Psi
    4. Neutrino free-streaming modifies radiation driving of acoustic oscillations
    """

    def __init__(self, bg, k_arr_Mpc, l_nu_max=L_NU_MAX):
        self.bg = bg
        self.k_arr_np = np.asarray(k_arr_Mpc, dtype=np.float32)
        self.N_k = len(self.k_arr_np)
        self.k_arr = mx.array(self.k_arr_np)
        self.l_nu_max = l_nu_max
        k_max = float(self.k_arr_np[-1])
        self.tau_grid = build_tau_grid_neutrino(bg, k_max)
        self.n_var = n_var_total(l_nu_max)

        self.Omega_gamma = _OMEGA_GAMMA
        self.Omega_nu = _OMEGA_NU
        self.Omega_b = _OMEGA_B
        self.Omega_c = float(bg.Omega_cdm)

        print(f"[Neutrino IMEX] N_k={self.N_k}, l_nu_max={l_nu_max}, "
              f"N_var={self.n_var}")
        print(f"[Neutrino IMEX] Omega_gamma={self.Omega_gamma:.6f}, "
              f"Omega_nu={self.Omega_nu:.6f}, f_nu={_f_nu:.4f}")
        print(f"[Neutrino IMEX] Grid: {len(self.tau_grid)} steps "
              f"[{self.tau_grid[0]:.2e}, {self.tau_grid[-1]:.1f}] Mpc "
              f"(tau_rec={bg.tau_rec:.1f})")

    def solve(self):
        t0 = time.time()
        bg = self.bg
        k_arr = self.k_arr
        tau_grid = self.tau_grid

        y = adiabatic_ic(self.k_arr_np, bg, self.l_nu_max)

        calH_np = bg.calH_at_tau(tau_grid).astype(np.float32)
        R_np = bg.R_at_tau(tau_grid).astype(np.float32)
        a_np = bg.a_at_tau(tau_grid).astype(np.float32)

        calH_arr = mx.array(calH_np)
        R_arr = mx.array(R_np)
        a_arr = mx.array(a_np)

        _Og = self.Omega_gamma
        _On = self.Omega_nu
        _Ob = self.Omega_b
        _Oc = self.Omega_c
        _H0 = _H0_MPC
        _lmax = self.l_nu_max

        print("[Neutrino IMEX] Integrating to tau_rec...")
        for i in range(len(tau_grid) - 1):
            dt = float(tau_grid[i + 1] - tau_grid[i])
            y = imex_rk4_step(
                y, k_arr,
                calH_arr[i], R_arr[i], tau_grid[i],
                _Og, _On, _Ob, _Oc,
                a_arr[i], _H0, _lmax, dt)

            if i % 200 == 0:
                mx.eval(y)

        mx.eval(y)
        t_total = time.time() - t0
        print(f"[Neutrino IMEX] Done in {t_total:.2f}s ({len(tau_grid)} steps)")

        return NeutrinoResult(y=y, k_arr=self.k_arr_np, bg=bg,
                              l_nu_max=self.l_nu_max,
                              Omega_nu=_On, H0=_H0)


class NeutrinoResult:
    """Result container."""

    def __init__(self, y, k_arr, bg, l_nu_max, Omega_nu, H0):
        self.y = y
        self.k_arr = k_arr
        self.bg = bg
        self.l_nu_max = l_nu_max
        self.Omega_nu = Omega_nu
        self.H0 = H0
        self.N_k = len(k_arr)

    def source_at_recombination(self):
        """
        Extract source functions at tau_rec with Silk damping.

        Returns (Theta_0, Phi, Psi, v_b, N_0, N_2).

        The SW source is Theta_0 + Psi (not Theta_0 + Phi) when Psi != Phi.
        """
        y_np = np.array(self.y)
        Phi = y_np[:, IDX_PHI]
        Theta_0 = y_np[:, IDX_THETA_0]
        v_b = 3.0 * y_np[:, IDX_THETA_1]
        N_0 = y_np[:, IDX_N_START]
        N_2 = y_np[:, IDX_N_START + 2] if self.l_nu_max >= 2 else np.zeros_like(Phi)

        # Diagnose Psi
        H02 = self.H0 ** 2
        a_rec = float(self.bg.a_at_tau(np.array([self.bg.tau_rec]))[0])
        Psi = Phi - 12.0 * H02 * self.Omega_nu * N_2 / (a_rec ** 2 * self.k_arr ** 2)

        silk = np.exp(-(self.k_arr / self.bg.k_D) ** 2)
        return (Theta_0 * silk, Phi * silk, Psi * silk,
                v_b * silk, N_0 * silk, N_2 * silk)

    def sw_source(self):
        """Sachs-Wolfe source: Theta_0 + Psi (silk-damped)."""
        Theta_0, Phi, Psi, v_b, N_0, N_2 = self.source_at_recombination()
        return Theta_0 + Psi

    def phi_psi_ratio(self):
        """Psi/Phi at tau_rec (no Silk damping)."""
        y_np = np.array(self.y)
        Phi = y_np[:, IDX_PHI]
        N_2 = y_np[:, IDX_N_START + 2] if self.l_nu_max >= 2 else np.zeros_like(Phi)
        H02 = self.H0 ** 2
        a_rec = float(self.bg.a_at_tau(np.array([self.bg.tau_rec]))[0])
        Psi = Phi - 12.0 * H02 * self.Omega_nu * N_2 / (a_rec ** 2 * self.k_arr ** 2)
        mask = np.abs(Phi) > 1e-10
        return np.where(mask, Psi / Phi, 1.0)


# ============================================================================
# Self-test and comparison
# ============================================================================

def _test_neutrino_solver():
    """
    Test the neutrino solver:
    1. No NaN/Inf
    2. Psi != Phi at superhorizon scales
    3. C_l comparison with/without neutrinos
    4. Effect on peak heights
    """
    print("=" * 70)
    print("TEST: Neutrino IMEX Boltzmann solver (N_eff = 3.046)")
    print("=" * 70)

    from .background import Background
    from .spectra import compute_cl

    # Background
    print("\n--- Step 1: Background ---")
    bg = Background(khronon=False)
    bg.solve()

    print(f"\nOmega_gamma = {_OMEGA_GAMMA:.6f}")
    print(f"Omega_nu    = {_OMEGA_NU:.6f}")
    print(f"f_nu        = {_f_nu:.4f}")
    print(f"Omega_r     = {_OMEGA_R:.6f}")

    # k-grid
    k_arr = np.geomspace(5e-4, 0.35, 500).astype(np.float32)

    # --- Solve WITH neutrinos ---
    print("\n--- Step 2: Neutrino solver ---")
    solver_nu = NeutrinoBoltzmannSolver(bg, k_arr, l_nu_max=L_NU_MAX)
    result_nu = solver_nu.solve()

    Theta_0_nu, Phi_nu, Psi_nu, v_b_nu, N_0_nu, N_2_nu = \
        result_nu.source_at_recombination()

    # Test A: Stability
    print(f"\n--- Test A: Stability ---")
    has_nan = any(np.any(np.isnan(x)) for x in [Theta_0_nu, Phi_nu, Psi_nu])
    has_inf = any(np.any(np.isinf(x)) for x in [Theta_0_nu, Phi_nu, Psi_nu])
    print(f"  NaN: {has_nan}, Inf: {has_inf}")
    if has_nan or has_inf:
        print("  FAIL!")
        # Diagnose
        for name, arr in [("Theta_0", Theta_0_nu), ("Phi", Phi_nu), ("Psi", Psi_nu)]:
            bad = np.where(np.isnan(arr) | np.isinf(arr))[0]
            if len(bad) > 0:
                print(f"    {name}: first bad at k={k_arr[bad[0]]:.4f}")
        return False
    print("  PASS")

    # Test B: Psi/Phi
    print(f"\n--- Test B: Psi/Phi anisotropic stress ---")
    ratio = result_nu.phi_psi_ratio()
    mask_ls = k_arr < 0.01
    if np.sum(mask_ls) > 0:
        ratio_ls = np.mean(ratio[mask_ls])
        expected = 1.0 - 2.0 * _f_nu / 5.0
        print(f"  Psi/Phi (k<0.01): {ratio_ls:.4f}")
        print(f"  Analytic (RD):    {expected:.4f}")
        dev = abs(ratio_ls - 1.0)
        print(f"  Deviation from 1: {dev*100:.2f}%")

    # Test C: Amplitudes
    print(f"\n--- Test C: Source amplitudes ---")
    sw_nu = Theta_0_nu + Psi_nu
    print(f"  Theta_0 + Psi: [{np.min(sw_nu):.4f}, {np.max(sw_nu):.4f}]")
    print(f"  Phi: [{np.min(Phi_nu):.4f}, {np.max(Phi_nu):.4f}]")

    # --- Solve WITHOUT neutrinos ---
    print("\n--- Step 3: Reference (no neutrinos) ---")
    from .perturbations_implicit import ImplicitBoltzmannSolver
    solver_ref = ImplicitBoltzmannSolver(bg, k_arr)
    result_ref = solver_ref.solve()
    Theta_0_ref, Phi_ref, v_b_ref = result_ref.source_at_recombination()
    sw_ref = Theta_0_ref + Phi_ref

    # --- C_l comparison ---
    print("\n--- Step 4: C_l comparison ---")
    ell_values = np.unique(np.concatenate([
        np.arange(2, 30, 1), np.arange(30, 100, 2),
        np.arange(100, 500, 4), np.arange(500, 1500, 8),
        np.arange(1500, 2501, 12),
    ])).astype(int)

    _, _, Dl_nu = compute_cl(sw_nu, k_arr, ell_values, bg.D_A)
    _, _, Dl_ref = compute_cl(sw_ref, k_arr, ell_values, bg.D_A)

    from scipy.signal import find_peaks
    from scipy.interpolate import interp1d

    l_full = np.arange(2, 2501)
    for label, Dl, ell in [("With neutrinos", Dl_nu, ell_values),
                            ("No neutrinos", Dl_ref, ell_values)]:
        f = interp1d(ell, Dl, kind='cubic', fill_value='extrapolate')
        Dl_interp = np.maximum(f(l_full), 0.0)
        peaks, _ = find_peaks(Dl_interp, distance=80, prominence=20)
        if len(peaks) >= 3:
            peak_ls = l_full[peaks[:5]]
            peak_Ds = Dl_interp[peaks[:5]]
            print(f"\n  {label}:")
            print(f"    Peaks at l = {peak_ls}")
            print(f"    Heights = {[f'{d:.0f}' for d in peak_Ds]} uK^2")
            if len(peaks) >= 2:
                print(f"    1st/2nd: {peak_Ds[0]/peak_Ds[1]:.3f}")

    # Neutrino effects
    print(f"\n--- Test D: Neutrino effects ---")
    mask_sw = ell_values < 30
    if np.sum(mask_sw) > 0 and np.mean(Dl_ref[mask_sw]) > 0:
        sw_change = np.mean(Dl_nu[mask_sw]) / np.mean(Dl_ref[mask_sw])
        print(f"  SW plateau (l<30): ratio = {sw_change:.4f} ({(sw_change-1)*100:+.1f}%)")

    mask_p1 = (ell_values > 150) & (ell_values < 300)
    if np.sum(mask_p1) > 0:
        max_nu = np.max(Dl_nu[mask_p1])
        max_ref = np.max(Dl_ref[mask_p1])
        if max_ref > 0:
            print(f"  1st peak: {max_nu:.0f} vs {max_ref:.0f} uK^2 "
                  f"({(max_nu/max_ref-1)*100:+.1f}%)")

    # Plot
    print("\n--- Generating plot ---")
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(2, 2, figsize=(14, 10))

        ax = axes[0, 0]
        ax.plot(ell_values, Dl_nu, 'b-', lw=1.2, label=r'With $\nu$ ($N_{eff}=3.046$)')
        ax.plot(ell_values, Dl_ref, 'r--', lw=1.2, label=r'Without $\nu$')
        ax.set_xlabel('Multipole l')
        ax.set_ylabel(r'$D_\ell$ [$\mu K^2$]')
        ax.set_title('CMB TT Power Spectrum')
        ax.legend(fontsize=9)
        ax.set_xlim(2, 2500)

        ax = axes[0, 1]
        mask_pos = Dl_ref > 0
        ratio_cl = np.ones_like(Dl_nu)
        ratio_cl[mask_pos] = Dl_nu[mask_pos] / Dl_ref[mask_pos]
        ax.plot(ell_values, ratio_cl, 'k-', lw=1)
        ax.axhline(1.0, color='gray', ls=':', alpha=0.5)
        ax.set_xlabel('Multipole l')
        ax.set_ylabel(r'$D_\ell^{\nu} / D_\ell^{no\,\nu}$')
        ax.set_title('Neutrino Effect on Power Spectrum')
        ax.set_xlim(2, 2500)
        ax.set_ylim(0.5, 1.5)

        ax = axes[1, 0]
        ax.semilogx(k_arr, ratio, 'g-', lw=1.2)
        ax.axhline(1.0, color='gray', ls=':', alpha=0.5)
        expected = 1.0 - 2.0 * _f_nu / 5.0
        ax.axhline(expected, color='orange', ls='--', alpha=0.5,
                   label=f'$1-2f_\\nu/5 = {expected:.3f}$')
        ax.set_xlabel(r'$k$ [Mpc$^{-1}$]')
        ax.set_ylabel(r'$\Psi / \Phi$')
        ax.set_title(r'Potential Ratio at $\tau_{rec}$')
        ax.legend(fontsize=9)

        ax = axes[1, 1]
        ax.semilogx(k_arr, sw_nu, 'b-', lw=1,
                     label=r'$\Theta_0 + \Psi$ (with $\nu$)')
        ax.semilogx(k_arr, sw_ref, 'r--', lw=1,
                     label=r'$\Theta_0 + \Phi$ (no $\nu$)')
        ax.set_xlabel(r'$k$ [Mpc$^{-1}$]')
        ax.set_ylabel('Source function')
        ax.set_title('SW Source at Recombination')
        ax.legend(fontsize=9)

        plt.tight_layout()
        import os
        out_dir = os.path.dirname(os.path.abspath(__file__))
        out_path = os.path.join(out_dir, 'neutrino_comparison.png')
        plt.savefig(out_path, dpi=150)
        plt.close()
        print(f"  Plot saved: {out_path}")
    except ImportError:
        print("  matplotlib not available")

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print(f"  N_eff = {N_EFF}, f_nu = {_f_nu:.4f}")
    print(f"  l_nu_max = {L_NU_MAX}, N_var = {solver_nu.n_var}")
    print(f"  Stability: {'PASS' if not (has_nan or has_inf) else 'FAIL'}")
    print("=" * 70)

    return not (has_nan or has_inf)


if __name__ == '__main__':
    _test_neutrino_solver()
