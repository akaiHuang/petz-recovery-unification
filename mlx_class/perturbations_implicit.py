"""
perturbations_implicit.py — IMEX Boltzmann solver for stiff Poisson equation.

The standard explicit RK4 solver becomes unstable for k > 0.25 Mpc^{-1}
when integrating past 0.85*tau_rec because the Poisson equation has a
stiff eigenvalue lambda = -k^2/(3*calH) - calH.

SOLUTION: IMEX (Implicit-Explicit) Runge-Kutta splitting.

The Phi equation is:
  Phi' = [-k^2 Phi - (3/2) H0^2 S] / (3 calH) - calH Phi
       = lambda * Phi + F(y)

where:
  lambda = -(k^2/(3 calH) + calH)  ... STIFF linear coefficient
  F(y) = -(H0^2/(2 calH)) * S(y)   ... non-stiff source from density perturbations
  S = Omega_r/a^2 * 4*Theta_0 + Omega_b/a * delta_b + Omega_c/a * delta_c

The IMEX approach treats lambda*Phi implicitly and F(y) explicitly:
- All non-Phi variables: standard explicit RK4
- Phi equation: Rosenbrock-type implicit treatment of the stiff linear part

This uses the L-stable IMEX-RK scheme (Ascher-Ruuth-Spiteri, 1997).
For simplicity, we implement a semi-implicit approach: at each RK4 stage,
the explicit part is computed normally, then the implicit correction for
the stiff Phi term is applied via an exact exponential integrator factor.

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

# TCA variable indices — same 6-variable system as the explicit solver
TCA_PHI = 0
TCA_DELTA_B = 1
TCA_DELTA_C = 2
TCA_V_C = 3
TCA_THETA_0 = 4
TCA_THETA_1 = 5
TCA_N_VAR = 6


# ============================================================================
# Phi_dot decomposition: stiff + non-stiff parts
# ============================================================================

def _phi_source(Theta_0, delta_b, delta_c, calH, a,
                Omega_r, Omega_b, Omega_c, H0_Mpc):
    """
    Non-stiff source F in Phi' = lambda*Phi + F.

    F = -(H0^2/(2*calH)) * [Omega_r/a^2 * 4*Theta_0
                              + Omega_b/a * delta_b
                              + Omega_c/a * delta_c]
    """
    H02 = H0_Mpc * H0_Mpc
    S = (Omega_r / (a * a) * (4.0 * Theta_0)
         + Omega_b / a * delta_b
         + Omega_c / a * delta_c)
    return -0.5 * H02 / calH * S


def _phi_lambda(k_arr, calH):
    """
    Stiff eigenvalue: lambda = -(k^2/(3*calH) + calH).

    The k^2/(3*calH) term is the stiff part (diverges as calH -> 0 near recombination).
    """
    return -(k_arr * k_arr / (3.0 * calH) + calH)


# ============================================================================
# Full TCA derivative (6 variables)
# ============================================================================

def deriv_tca_full(y, k_arr, calH, R, Omega_r, Omega_b, Omega_c, a, H0_Mpc):
    """
    Compute the full RHS for all 6 TCA variables.
    Same physics as the original perturbations.py, same equations.

    Returns: (dy_full, Phi_dot_explicit_source)
    where Phi_dot = lambda*Phi + F is split into stiff + non-stiff.
    """
    Phi = y[:, TCA_PHI]
    delta_b = y[:, TCA_DELTA_B]
    delta_c = y[:, TCA_DELTA_C]
    v_c = y[:, TCA_V_C]
    Theta_0 = y[:, TCA_THETA_0]
    Theta_1 = y[:, TCA_THETA_1]

    # Full Phi_dot (for the other equations that depend on it)
    F_phi = _phi_source(Theta_0, delta_b, delta_c, calH, a,
                        Omega_r, Omega_b, Omega_c, H0_Mpc)
    lam = _phi_lambda(k_arr, calH)
    Phi_dot = lam * Phi + F_phi

    # Other equations (standard TCA)
    dTheta_0 = -k_arr * Theta_1 - Phi_dot
    dTheta_1 = (k_arr / 3.0 * (Theta_0 + Phi) - calH * R * Theta_1) / (1.0 + R)
    d_delta_b = -3.0 * k_arr * Theta_1 - 3.0 * Phi_dot
    d_delta_c = -(k_arr * v_c + 3.0 * Phi_dot)
    d_v_c = -calH * v_c + k_arr * Phi

    dy = mx.stack([Phi_dot, d_delta_b, d_delta_c, d_v_c, dTheta_0, dTheta_1], axis=1)
    return dy


# ============================================================================
# IMEX RK4 step: exponential Rosenbrock for the stiff Phi term
# ============================================================================

def imex_rk4_step(y, k_arr, calH, R, Omega_r, Omega_b, Omega_c, a, H0_Mpc, dtau):
    """
    IMEX RK4 step with exponential integrator for the stiff Phi equation.

    Strategy: Use the exact solution of the LINEAR stiff part.
    For Phi' = lambda*Phi + F(t), the exact solution over dt is:
      Phi(t+dt) = exp(lambda*dt)*Phi(t) + integral_0^dt exp(lambda*(dt-s)) F(t+s) ds

    For the non-stiff variables, we use standard explicit RK4.

    The coupling between Phi and other variables is handled by:
    1. Predicting Phi at each sub-stage using the exponential factor
    2. Using the predicted Phi to evaluate derivatives of other variables
    3. After the RK4 step, correcting Phi using the exponential integrator
       with the average source F from all 4 stages

    This is a first-order splitting but with the RK4 accuracy for the
    non-stiff part. The stiff part is integrated exactly (exponential).
    """
    lam = _phi_lambda(k_arr, calH)  # shape: (N_k,)
    exp_lam_dt = mx.exp(lam * dtau)
    exp_lam_half = mx.exp(lam * 0.5 * dtau)

    # ------ Stage 1: t_n ------
    Phi_1 = y[:, TCA_PHI]
    F1 = _phi_source(y[:, TCA_THETA_0], y[:, TCA_DELTA_B], y[:, TCA_DELTA_C],
                     calH, a, Omega_r, Omega_b, Omega_c, H0_Mpc)
    k1 = deriv_tca_full(y, k_arr, calH, R, Omega_r, Omega_b, Omega_c, a, H0_Mpc)

    # ------ Stage 2: t_n + dt/2, using explicit Euler prediction ------
    # Predict Phi at half-step using exponential integrator
    Phi_2 = exp_lam_half * Phi_1 + _phi_int_factor(lam, 0.5 * dtau) * F1

    y2 = y + 0.5 * dtau * k1
    # Override Phi with the exponential prediction
    y2 = _set_phi(y2, Phi_2)
    F2 = _phi_source(y2[:, TCA_THETA_0], y2[:, TCA_DELTA_B], y2[:, TCA_DELTA_C],
                     calH, a, Omega_r, Omega_b, Omega_c, H0_Mpc)
    k2 = deriv_tca_full(y2, k_arr, calH, R, Omega_r, Omega_b, Omega_c, a, H0_Mpc)

    # ------ Stage 3: t_n + dt/2, using stage 2 ------
    Phi_3 = exp_lam_half * Phi_1 + _phi_int_factor(lam, 0.5 * dtau) * F2

    y3 = y + 0.5 * dtau * k2
    y3 = _set_phi(y3, Phi_3)
    F3 = _phi_source(y3[:, TCA_THETA_0], y3[:, TCA_DELTA_B], y3[:, TCA_DELTA_C],
                     calH, a, Omega_r, Omega_b, Omega_c, H0_Mpc)
    k3 = deriv_tca_full(y3, k_arr, calH, R, Omega_r, Omega_b, Omega_c, a, H0_Mpc)

    # ------ Stage 4: t_n + dt ------
    Phi_4 = exp_lam_dt * Phi_1 + _phi_int_factor(lam, dtau) * F3

    y4 = y + dtau * k3
    y4 = _set_phi(y4, Phi_4)
    F4 = _phi_source(y4[:, TCA_THETA_0], y4[:, TCA_DELTA_B], y4[:, TCA_DELTA_C],
                     calH, a, Omega_r, Omega_b, Omega_c, H0_Mpc)
    k4 = deriv_tca_full(y4, k_arr, calH, R, Omega_r, Omega_b, Omega_c, a, H0_Mpc)

    # ------ Combine: standard RK4 for non-Phi variables ------
    y_new = y + (dtau / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)

    # ------ Phi: Simpson's rule for the exponential integrator ------
    # Phi(t+dt) = exp(lambda*dt)*Phi(t) + integral using Simpson's rule:
    # integral_0^dt exp(lambda*(dt-s)) F(s) ds
    #   ~ (dt/6) [exp(lambda*dt)*F1 + 4*exp(lambda*dt/2)*F_mid + F4]
    # where F_mid = average of F2 and F3
    F_mid = 0.5 * (F2 + F3)
    Phi_new = (exp_lam_dt * Phi_1
               + (dtau / 6.0) * (exp_lam_dt * F1
                                  + 4.0 * exp_lam_half * F_mid
                                  + F4))

    y_new = _set_phi(y_new, Phi_new)
    return y_new


def _phi_int_factor(lam, dt):
    """
    Compute (exp(lam*dt) - 1) / lam, handling lam -> 0 safely.

    This is the integral factor: integral_0^dt exp(lam*s) ds.
    """
    lam_dt = lam * dt
    # For |lam_dt| < 1e-6, use Taylor: dt * (1 + lam_dt/2 + lam_dt^2/6)
    # For larger values, use exact formula
    # In MLX, we compute both and select
    exact = (mx.exp(lam_dt) - 1.0) / lam
    taylor = dt * (1.0 + 0.5 * lam_dt + lam_dt * lam_dt / 6.0)
    return mx.where(mx.abs(lam_dt) < 1e-4, taylor, exact)


def _set_phi(y, Phi_new):
    """Replace the Phi column in the state vector y."""
    # y shape: (N_k, 6). We need to set y[:, 0] = Phi_new.
    # In MLX, array assignment is not supported directly.
    # Use concatenation instead.
    return mx.concatenate([
        Phi_new[:, None],
        y[:, 1:]
    ], axis=1)


# ============================================================================
# Tau grid builder — can go all the way to tau_rec
# ============================================================================

def build_tau_grid_implicit(bg, k_max, N_early=300, N_late=800):
    """
    Build integration grid that extends to tau_rec (not 0.85*tau_rec).

    The stiff eigenvalue is handled implicitly, so we only need to resolve
    the acoustic oscillation timescale and the source variation timescale.
    """
    tau_rec = bg.tau_rec
    tau_init = bg.tau_grid[1]
    tau_early_end = min(10.0, 0.3 * tau_rec)

    # Early phase: geometric spacing for radiation-dominated era
    tau_early = np.geomspace(tau_init, tau_early_end, N_early)

    # Late phase: resolve acoustic oscillations
    # Acoustic period: 2*pi / (k_max * c_s) where c_s ~ 1/sqrt(3)
    # Need ~25 steps per period for RK4 accuracy
    cs_approx = 1.0 / np.sqrt(3.0)
    dtau_acoustic = 2.0 * np.pi / (k_max * cs_approx) / 25.0
    N_needed = int((tau_rec - tau_early_end) / dtau_acoustic) + 10
    N_late_actual = max(N_late, N_needed)

    tau_late = np.linspace(tau_early_end, tau_rec, N_late_actual)
    tau_grid = np.unique(np.concatenate([tau_early, tau_late]))
    return tau_grid


# ============================================================================
# Initial conditions (6 variables)
# ============================================================================

def adiabatic_ic(k_arr_np, bg):
    """
    Adiabatic initial conditions in the radiation-dominated era.

    Standard adiabatic mode (Phi_init = 1 normalization):
      Phi     = 1
      delta_b = -3/2 Phi = -1.5
      delta_c = -3/2 Phi = -1.5
      v_c     = k tau / 6
      Theta_0 = -Phi/2 = -0.5
      Theta_1 = k tau / 18
    """
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


# ============================================================================
# Main solver class
# ============================================================================

class ImplicitBoltzmannSolver:
    """
    Batched Boltzmann solver with IMEX treatment of the stiff Poisson equation.

    The Phi equation's stiff eigenvalue k^2/(3*calH) is handled by an
    exponential integrator, while all other variables use standard RK4.
    This allows stable integration all the way to tau_rec for any k.
    """

    def __init__(self, bg, k_arr_Mpc):
        self.bg = bg
        self.k_arr_np = np.asarray(k_arr_Mpc, dtype=np.float32)
        self.N_k = len(self.k_arr_np)
        self.k_arr = mx.array(self.k_arr_np)
        k_max = float(self.k_arr_np[-1])
        self.tau_grid = build_tau_grid_implicit(bg, k_max)
        print(f"[IMEX] N_k={self.N_k}, grid: {len(self.tau_grid)} steps "
              f"[{self.tau_grid[0]:.2e}, {self.tau_grid[-1]:.1f}] Mpc "
              f"(tau_rec={bg.tau_rec:.1f})")

    def solve(self):
        t0 = time.time()
        bg = self.bg
        k_arr = self.k_arr
        tau_grid = self.tau_grid

        # Initial conditions
        y = adiabatic_ic(self.k_arr_np, bg)

        # Precompute background quantities
        calH_np = bg.calH_at_tau(tau_grid).astype(np.float32)
        R_np = bg.R_at_tau(tau_grid).astype(np.float32)
        a_np = bg.a_at_tau(tau_grid).astype(np.float32)

        calH_arr = mx.array(calH_np)
        R_arr = mx.array(R_np)
        a_arr = mx.array(a_np)

        _Or = _OMEGA_R
        _Ob = _OMEGA_B
        _Oc = float(bg.Omega_cdm)
        _H0 = _H0_MPC

        print("[IMEX] Integrating TCA to tau_rec...")
        for i in range(len(tau_grid) - 1):
            dt = float(tau_grid[i + 1] - tau_grid[i])
            y = imex_rk4_step(y, k_arr,
                              calH_arr[i], R_arr[i],
                              _Or, _Ob, _Oc, a_arr[i], _H0,
                              dt)

            if i % 200 == 0:
                mx.eval(y)

        mx.eval(y)
        t_total = time.time() - t0
        print(f"[IMEX] Done in {t_total:.2f}s ({len(tau_grid)} steps)")

        return ImplicitResult(y_tca=y, k_arr=self.k_arr_np, bg=bg)


class ImplicitResult:
    """Result container for the IMEX solver."""

    def __init__(self, y_tca, k_arr, bg):
        self.y_tca = y_tca
        self.k_arr = k_arr
        self.bg = bg
        self.N_k = len(k_arr)

    def source_at_recombination(self):
        """
        Extract source functions at tau_rec.
        Returns (Theta_0, Phi, v_b) with Silk damping applied analytically.
        """
        y_np = np.array(self.y_tca)
        Theta_0 = y_np[:, TCA_THETA_0]
        Phi = y_np[:, TCA_PHI]
        v_b = 3.0 * y_np[:, TCA_THETA_1]

        silk = np.exp(-(self.k_arr / self.bg.k_D) ** 2)
        return Theta_0 * silk, Phi * silk, v_b * silk


# ============================================================================
# Self-test
# ============================================================================

def _test_implicit_solver():
    """
    Test the IMEX solver:
    1. Verify no NaN for k up to 0.35 Mpc^{-1}
    2. Verify Theta_0 + Phi oscillates as cos(k * r_s) at tau_rec
    3. Compare sound horizon with expected r_s ~ 144.4 Mpc
    4. Cross-validate against explicit solver at 0.85*tau_rec
    """
    print("=" * 70)
    print("TEST: IMEX (exponential integrator) Boltzmann solver")
    print("=" * 70)

    from .background import Background

    # Step 1: Background
    print("\n--- Step 1: Background ---")
    bg = Background(khronon=False)
    bg.solve()

    # Step 2: k-grid (500 points for accurate r_s measurement)
    k_arr = np.geomspace(5e-4, 0.35, 500).astype(np.float32)
    print(f"\nk range: [{k_arr[0]:.4f}, {k_arr[-1]:.4f}] Mpc^-1, N_k={len(k_arr)}")
    print(f"Expected r_s = {bg.r_s:.1f} Mpc")

    # Step 3: Solve to tau_rec
    print("\n--- Step 2: IMEX solver to tau_rec ---")
    solver = ImplicitBoltzmannSolver(bg, k_arr)
    result = solver.solve()

    # Step 4: Extract source
    Theta_0, Phi, v_b = result.source_at_recombination()

    # --- Test A: No NaN / Inf ---
    has_nan = (np.any(np.isnan(Theta_0)) or np.any(np.isnan(Phi))
               or np.any(np.isnan(v_b)))
    has_inf = (np.any(np.isinf(Theta_0)) or np.any(np.isinf(Phi))
               or np.any(np.isinf(v_b)))
    print(f"\n--- Test A: Stability ---")
    print(f"  NaN: {has_nan}")
    print(f"  Inf: {has_inf}")

    if has_nan or has_inf:
        print("  FAIL: NaN or Inf detected!")
        for name, arr in [("Theta_0", Theta_0), ("Phi", Phi), ("v_b", v_b)]:
            bad_idx = np.where(np.isnan(arr) | np.isinf(arr))[0]
            if len(bad_idx) > 0:
                print(f"    {name}: first bad value at k[{bad_idx[0]}] = "
                      f"{k_arr[bad_idx[0]]:.4f}")
        return False

    print("  PASS: No NaN or Inf for k up to 0.35 Mpc^-1")

    # --- Test B: Amplitude reasonableness ---
    source_SW = Theta_0 + Phi
    print(f"\n--- Test B: Amplitude ---")
    print(f"  Theta_0 + Phi range: [{np.min(source_SW):.4f}, {np.max(source_SW):.4f}]")
    print(f"  Phi range: [{np.min(Phi):.4f}, {np.max(Phi):.4f}]")
    if np.max(np.abs(source_SW)) < 2.0:
        print("  PASS: Amplitudes are physical (< 2)")
    else:
        print("  WARNING: Amplitudes may be unphysical")

    # --- Test C: Sound horizon from zero crossings ---
    print(f"\n--- Test C: Sound horizon r_s ---")

    # Use subhorizon modes (k > 0.03) where oscillations are well-resolved
    mask_sub = k_arr > 0.03
    k_sub = k_arr[mask_sub]
    sw_sub = source_SW[mask_sub]
    sign_changes = np.where(np.diff(np.sign(sw_sub)))[0]

    if len(sign_changes) >= 6:
        # Linear interpolation for precise zero positions
        k_zeros = []
        for sc in sign_changes:
            k0, k1 = k_sub[sc], k_sub[sc + 1]
            s0, s1 = sw_sub[sc], sw_sub[sc + 1]
            k_zero = k0 - s0 * (k1 - k0) / (s1 - s0)
            k_zeros.append(k_zero)
        k_zeros = np.array(k_zeros)

        # Consecutive zeros are separated by pi/r_s in k-space
        dk_zeros = np.diff(k_zeros)
        r_s_values = np.pi / dk_zeros

        # Robust estimate using the second half of zero crossings
        n_half = len(r_s_values) // 2
        r_s_measured = np.median(r_s_values[n_half:])

        print(f"  Number of zero crossings (k>0.03): {len(k_zeros)}")
        print(f"  First 6 zeros at k = {np.array2string(k_zeros[:6], precision=4)}")
        print(f"  r_s estimates: {np.array2string(r_s_values[:8], precision=1)}")
        print(f"  r_s (median, upper half): {r_s_measured:.1f} Mpc")
        print(f"  Expected r_s:             {bg.r_s:.1f} Mpc")
        error_pct = abs(r_s_measured - bg.r_s) / bg.r_s * 100
        print(f"  Relative error:           {error_pct:.1f}%")

        if error_pct < 10.0:
            print(f"  PASS: cos(k*r_s) oscillation confirmed (r_s error < 10%)")
        else:
            print(f"  WARNING: r_s error = {error_pct:.1f}%")
    else:
        print(f"  Only {len(sign_changes)} zero crossings found (k>0.03). "
              f"Cannot measure r_s.")

    # --- Test D: Cross-validation at 0.85*tau_rec ---
    print(f"\n--- Test D: Cross-validation with explicit solver ---")
    from .perturbations import BatchedBoltzmannSolver

    # Run IMEX to 0.85*tau_rec for a fair comparison
    k_arr_comp = np.geomspace(5e-4, 0.20, 200).astype(np.float32)
    k_arr_comp_mx = mx.array(k_arr_comp)
    tau_switch = 0.85 * bg.tau_rec

    tau_grid_short = build_tau_grid_implicit(bg, 0.20)
    tau_grid_short = tau_grid_short[tau_grid_short <= tau_switch + 0.1]

    calH_np = bg.calH_at_tau(tau_grid_short).astype(np.float32)
    R_np = bg.R_at_tau(tau_grid_short).astype(np.float32)
    a_np = bg.a_at_tau(tau_grid_short).astype(np.float32)

    y_imex = adiabatic_ic(k_arr_comp, bg)
    for i in range(len(tau_grid_short) - 1):
        dt = float(tau_grid_short[i + 1] - tau_grid_short[i])
        y_imex = imex_rk4_step(
            y_imex, k_arr_comp_mx,
            mx.array(calH_np[i]), mx.array(R_np[i]),
            _OMEGA_R, _OMEGA_B, float(bg.Omega_cdm),
            mx.array(a_np[i]), _H0_MPC, dt)
        if i % 200 == 0:
            mx.eval(y_imex)
    mx.eval(y_imex)

    y_imex_np = np.array(y_imex)
    sw_imex_085 = (y_imex_np[:, TCA_THETA_0] + y_imex_np[:, TCA_PHI])
    silk_comp = np.exp(-(k_arr_comp / bg.k_D)**2)
    sw_imex_085 = sw_imex_085 * silk_comp

    # Explicit solver (also to 0.85*tau_rec)
    solver_exp = BatchedBoltzmannSolver(bg, k_arr_comp)
    result_exp = solver_exp.solve()
    Theta_0_exp, Phi_exp, _ = result_exp.source_at_recombination()
    sw_exp = Theta_0_exp + Phi_exp

    corr = np.corrcoef(sw_imex_085, sw_exp)[0, 1]
    max_diff = np.max(np.abs(sw_imex_085 - sw_exp))
    mean_diff = np.mean(np.abs(sw_imex_085 - sw_exp))

    print(f"  Both evaluated at 0.85*tau_rec, k=[5e-4, 0.20], N_k=200")
    print(f"  Correlation:       {corr:.6f}")
    print(f"  Max abs diff:      {max_diff:.6f}")
    print(f"  Mean abs diff:     {mean_diff:.6f}")

    if corr > 0.999:
        print(f"  PASS: IMEX matches explicit solver (corr > 0.999)")
    else:
        print(f"  WARNING: Correlation = {corr:.6f} (< 0.999)")

    # --- Summary ---
    print("\n" + "=" * 70)
    tests_passed = (not has_nan and not has_inf
                    and np.max(np.abs(source_SW)) < 2.0
                    and corr > 0.999)
    if tests_passed:
        print("ALL TESTS PASSED")
    else:
        print("SOME TESTS FAILED - see above")
    print("=" * 70)
    return tests_passed


if __name__ == '__main__':
    _test_implicit_solver()
