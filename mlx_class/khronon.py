"""
khronon.py — Khronon dark matter CMB solver with DBI correction.

Theory (Sigma = 2 ln Q, ghost condensation):
  - Background: K'(Q_0) = 0 => rho_K ~ a^{-3} exactly => IDENTICAL to LCDM
  - Perturbations: c_s^2 = 0 at ghost condensation => CDM-like
  - DBI correction: the full DBI kinetic structure K(Q) generates a
    k-dependent effective sound speed at next-to-leading order:

      c_s^2(k, a) = alpha_DBI * (k / k_J(a))^2 / (1 + (k / k_J(a))^2)

    where:
      k_J(a) = sqrt(1.5 * Omega_K * H0^2 / a) is the Khronon Jeans scale
      alpha_DBI = DBI NLO parameter controlling perturbation sound speed

    At ghost condensation: alpha_DBI = Q_0/2 = 0.5 (bare value)
    This is the UPPER BOUND. Physical values may be smaller due to
    loop corrections and renormalization of K(Q).

  The effect enters through the full GDM (Generalized Dark Matter) treatment:
    1. CDM Euler: pressure gradient  v_c' += c_s^2 * k * delta_c
    2. CDM continuity: pressure damping  delta_c' += -3*H*c_s^2*delta_c
    3. Poisson: GDM velocity correction to effective density

  Three comparison modes:
    Mode A: "Pure DBI" — same Omega, different c_s^2 only
    Mode B: "Full Khronon" — Omega_K=0.265 + DBI (both bg + perturbation)
    Mode C: "Background only" — Omega_K=0.265 but c_s^2=0 (background shift)

  This is the KEY differentiator: no other GPU solver implements this physics.

Author: Sheng-Kai Huang, 2026
"""
import os
os.environ['MPLCONFIGDIR'] = '/tmp/mpl_config'

import numpy as np
import mlx.core as mx
import time

from .background import (
    Background,
    Omega_r as _OMEGA_R,
    Omega_b as _OMEGA_B,
    Omega_c as _OMEGA_C,
    H0_Mpc as _H0_MPC,
    A_s, n_s, k_pivot, T_CMB,
    z_rec, a_rec,
)
from .perturbations_implicit import (
    ImplicitBoltzmannSolver,
    ImplicitResult,
    TCA_PHI, TCA_DELTA_B, TCA_DELTA_C, TCA_V_C,
    TCA_THETA_0, TCA_THETA_1, TCA_N_VAR,
    _phi_source, _phi_lambda, _phi_int_factor, _set_phi,
    build_tau_grid_implicit,
    adiabatic_ic,
)
from .spectra import compute_cl


# ============================================================================
# Khronon DBI parameters
# ============================================================================

Q_0 = 1.0
# Ghost condensation bare value: alpha_DBI = Q_0/2 = 0.5
# This is an upper bound; physical value depends on DBI normalization
ALPHA_DBI_DEFAULT = Q_0 / 2.0


def khronon_cs2(k_arr, a, Omega_K, H0_Mpc, alpha_dbi=ALPHA_DBI_DEFAULT):
    """
    DBI effective sound speed squared for the Khronon field.

    c_s^2(k, a) = alpha_DBI * x^2 / (1 + x^2)
    where x = k / k_J(a),  k_J(a) = sqrt(1.5 * Omega_K * H0^2 / a)

    Properties:
      - k << k_J:  c_s^2 ~ alpha * k^2 / k_J^2 -> 0  (CDM-like)
      - k >> k_J:  c_s^2 -> alpha_DBI               (DBI saturation)
      - k_J(a_rec) ~ 0.005 Mpc^{-1}
    """
    k_J_sq = 1.5 * Omega_K * H0_Mpc * H0_Mpc / a
    x_sq = k_arr * k_arr / k_J_sq
    cs2 = alpha_dbi * x_sq / (1.0 + x_sq)
    return cs2


# ============================================================================
# Modified Phi source (GDM Poisson equation)
# ============================================================================

def _phi_source_khronon(Theta_0, delta_b, delta_c, v_c, calH, a, k_arr,
                         Omega_r, Omega_b, Omega_K, H0_Mpc, cs2_k):
    """
    Poisson equation source with Khronon GDM pressure correction.

    Standard: S = Omega_r/a^2 * 4*Theta_0 + Omega_b/a * delta_b + Omega_K/a * delta_c
    GDM adds: + Omega_K/a * 3 * cs2 * calH * v_c / k
    (velocity correction from comoving-gauge effective density)
    """
    H02 = H0_Mpc * H0_Mpc
    S = (Omega_r / (a * a) * (4.0 * Theta_0)
         + Omega_b / a * delta_b
         + Omega_K / a * delta_c)
    # GDM pressure correction
    S = S + Omega_K / a * 3.0 * cs2_k * calH * v_c / (k_arr + 1e-30)
    return -0.5 * H02 / calH * S


# ============================================================================
# Modified TCA derivative (full GDM)
# ============================================================================

def deriv_tca_khronon(y, k_arr, calH, R, Omega_r, Omega_b, Omega_K,
                      a, H0_Mpc, cs2_k):
    """
    Full RHS for 6 TCA variables with Khronon GDM treatment.

    Three modifications from CDM (all proportional to cs2_k):
    1. Euler:      v_c' += cs2 * k * delta_c   (pressure gradient)
    2. Continuity: delta_c' += -3*calH*cs2*delta_c  (pressure damping)
    3. Poisson:    GDM velocity correction in source

    For cs2 -> 0: reduces exactly to standard CDM equations.
    """
    Phi = y[:, TCA_PHI]
    delta_b = y[:, TCA_DELTA_B]
    delta_c = y[:, TCA_DELTA_C]
    v_c = y[:, TCA_V_C]
    Theta_0 = y[:, TCA_THETA_0]
    Theta_1 = y[:, TCA_THETA_1]

    # Modified Poisson source
    F_phi = _phi_source_khronon(Theta_0, delta_b, delta_c, v_c, calH, a,
                                 k_arr, Omega_r, Omega_b, Omega_K, H0_Mpc,
                                 cs2_k)
    lam = _phi_lambda(k_arr, calH)
    Phi_dot = lam * Phi + F_phi

    # Photon-baryon (unchanged)
    dTheta_0 = -k_arr * Theta_1 - Phi_dot
    dTheta_1 = (k_arr / 3.0 * (Theta_0 + Phi) - calH * R * Theta_1) / (1.0 + R)
    d_delta_b = -3.0 * k_arr * Theta_1 - 3.0 * Phi_dot

    # Khronon GDM equations
    d_delta_c = -(k_arr * v_c + 3.0 * Phi_dot) - 3.0 * calH * cs2_k * delta_c
    d_v_c = -calH * v_c + k_arr * Phi + cs2_k * k_arr * delta_c

    dy = mx.stack([Phi_dot, d_delta_b, d_delta_c, d_v_c, dTheta_0, dTheta_1],
                  axis=1)
    return dy


# ============================================================================
# IMEX step with Khronon GDM
# ============================================================================

def imex_rk4_step_khronon(y, k_arr, calH, R, Omega_r, Omega_b, Omega_K,
                           a, H0_Mpc, dtau, cs2_k):
    """
    IMEX RK4 step: exponential integrator for stiff Phi + GDM corrections.
    """
    lam = _phi_lambda(k_arr, calH)
    exp_lam_dt = mx.exp(lam * dtau)
    exp_lam_half = mx.exp(lam * 0.5 * dtau)

    def _F(y_in):
        return _phi_source_khronon(
            y_in[:, TCA_THETA_0], y_in[:, TCA_DELTA_B],
            y_in[:, TCA_DELTA_C], y_in[:, TCA_V_C],
            calH, a, k_arr, Omega_r, Omega_b, Omega_K, H0_Mpc, cs2_k)

    def _K(y_in):
        return deriv_tca_khronon(y_in, k_arr, calH, R, Omega_r, Omega_b,
                                  Omega_K, a, H0_Mpc, cs2_k)

    Phi_1 = y[:, TCA_PHI]
    F1 = _F(y)
    k1 = _K(y)

    Phi_2 = exp_lam_half * Phi_1 + _phi_int_factor(lam, 0.5 * dtau) * F1
    y2 = y + 0.5 * dtau * k1
    y2 = _set_phi(y2, Phi_2)
    F2 = _F(y2)
    k2 = _K(y2)

    Phi_3 = exp_lam_half * Phi_1 + _phi_int_factor(lam, 0.5 * dtau) * F2
    y3 = y + 0.5 * dtau * k2
    y3 = _set_phi(y3, Phi_3)
    F3 = _F(y3)
    k3 = _K(y3)

    Phi_4 = exp_lam_dt * Phi_1 + _phi_int_factor(lam, dtau) * F3
    y4 = y + dtau * k3
    y4 = _set_phi(y4, Phi_4)
    F4 = _F(y4)
    k4 = _K(y4)

    y_new = y + (dtau / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)

    F_mid = 0.5 * (F2 + F3)
    Phi_new = (exp_lam_dt * Phi_1
               + (dtau / 6.0) * (exp_lam_dt * F1
                                  + 4.0 * exp_lam_half * F_mid
                                  + F4))
    y_new = _set_phi(y_new, Phi_new)
    return y_new


# ============================================================================
# Khronon Boltzmann Solver
# ============================================================================

class KhrononBoltzmannSolver:
    """
    Batched Boltzmann solver for the Khronon dark matter model.

    Implements the full GDM treatment: modified Euler, continuity, and
    Poisson equations with k-dependent DBI sound speed.
    """

    def __init__(self, bg, k_arr_Mpc, alpha_dbi=ALPHA_DBI_DEFAULT):
        self.bg = bg
        self.k_arr_np = np.asarray(k_arr_Mpc, dtype=np.float32)
        self.N_k = len(self.k_arr_np)
        self.k_arr = mx.array(self.k_arr_np)
        self.alpha_dbi = alpha_dbi

        k_max = float(self.k_arr_np[-1])
        self.tau_grid = build_tau_grid_implicit(bg, k_max)

        k_J_rec = np.sqrt(1.5 * bg.Omega_cdm * _H0_MPC**2 / a_rec)
        print(f"[Khronon IMEX] N_k={self.N_k}, alpha_DBI={alpha_dbi:.4f}, "
              f"grid: {len(self.tau_grid)} steps "
              f"[{self.tau_grid[0]:.2e}, {self.tau_grid[-1]:.1f}] Mpc")
        print(f"  k_J(z_rec) = {k_J_rec:.4f} Mpc^-1, l_J ~ {k_J_rec * bg.D_A:.0f}")

    def solve(self):
        t0 = time.time()
        bg = self.bg
        k_arr = self.k_arr
        tau_grid = self.tau_grid

        y = adiabatic_ic(self.k_arr_np, bg)

        calH_np = bg.calH_at_tau(tau_grid).astype(np.float32)
        R_np = bg.R_at_tau(tau_grid).astype(np.float32)
        a_np = bg.a_at_tau(tau_grid).astype(np.float32)

        calH_arr = mx.array(calH_np)
        R_arr = mx.array(R_np)
        a_arr = mx.array(a_np)

        _Or = _OMEGA_R
        _Ob = _OMEGA_B
        _Ok = float(bg.Omega_cdm)
        _H0 = _H0_MPC
        alpha = self.alpha_dbi

        print("[Khronon IMEX] Integrating with GDM equations...")
        for i in range(len(tau_grid) - 1):
            dt = float(tau_grid[i + 1] - tau_grid[i])
            cs2_k = khronon_cs2(k_arr, a_arr[i], _Ok, _H0, alpha)
            y = imex_rk4_step_khronon(
                y, k_arr, calH_arr[i], R_arr[i],
                _Or, _Ob, _Ok, a_arr[i], _H0, dt, cs2_k)
            if i % 200 == 0:
                mx.eval(y)

        mx.eval(y)
        t_total = time.time() - t0
        print(f"[Khronon IMEX] Done in {t_total:.2f}s ({len(tau_grid)} steps)")
        return KhrononResult(y_tca=y, k_arr=self.k_arr_np, bg=bg,
                             alpha_dbi=alpha)


class KhrononResult(ImplicitResult):
    """Result container for the Khronon solver."""

    def __init__(self, y_tca, k_arr, bg, alpha_dbi):
        super().__init__(y_tca, k_arr, bg)
        self.alpha_dbi = alpha_dbi

    def effective_cs2_at_rec(self):
        Omega_K = self.bg.Omega_cdm
        k_J_sq = 1.5 * Omega_K * _H0_MPC**2 / a_rec
        x_sq = self.k_arr**2 / k_J_sq
        return self.alpha_dbi * x_sq / (1.0 + x_sq)


# ============================================================================
# Helper: run a single model and return D_l
# ============================================================================

def _run_model(label, bg, k_arr, ell_values, alpha_dbi=None):
    """Run IMEX (or Khronon IMEX) and return D_l."""
    print(f"\n{'=' * 40}")
    print(f"  {label}")
    print(f"{'=' * 40}")

    if alpha_dbi is not None and alpha_dbi > 0:
        solver = KhrononBoltzmannSolver(bg, k_arr, alpha_dbi=alpha_dbi)
    else:
        solver = ImplicitBoltzmannSolver(bg, k_arr)

    result = solver.solve()
    Theta0, Phi, vb = result.source_at_recombination()
    source_SW = Theta0 + Phi
    ell_out, Cl, Dl = compute_cl(source_SW, k_arr, ell_values, bg.D_A)
    return Dl, source_SW, result


# ============================================================================
# Main comparison: 3-way (LCDM, Khronon pure DBI, Khronon full)
# ============================================================================

def run_comparison(N_k=500, ell_max=3000, alpha_dbi=ALPHA_DBI_DEFAULT,
                   save_dir=None):
    """
    Run LCDM and Khronon solvers with 3-way comparison:

    A. LCDM:          Omega_c standard, c_s^2 = 0
    B. Khronon (DBI): Omega_c standard, c_s^2 = DBI  (isolate perturbation effect)
    C. Khronon (full): Omega_K = 0.265, c_s^2 = DBI  (full prediction)

    This separates the two effects:
    - Residual (B - A): pure DBI perturbation effect
    - Residual (C - A): combined background + perturbation effect
    """
    if save_dir is None:
        save_dir = os.path.dirname(os.path.abspath(__file__))

    print("=" * 70)
    print("KHRONON vs LCDM: 3-WAY CMB COMPARISON")
    print(f"  N_k = {N_k}, ell_max = {ell_max}")
    print(f"  alpha_DBI = {alpha_dbi:.4f} (ghost condensation: Q_0/2)")
    print(f"  GPU: {mx.default_device()}")
    print("=" * 70)

    t_start = time.time()

    k_arr = np.geomspace(5e-5, 0.35, N_k).astype(np.float32)

    ell_values = np.unique(np.concatenate([
        np.arange(2, 30, 1),
        np.arange(30, 100, 2),
        np.arange(100, 500, 4),
        np.arange(500, 1500, 6),
        np.arange(1500, min(3001, ell_max + 1), 8),
    ])).astype(int)
    ell_values = ell_values[ell_values <= ell_max]

    # --- Model A: LCDM ---
    bg_lcdm = Background(khronon=False)
    bg_lcdm.solve()
    Dl_A, sw_A, res_A = _run_model(
        "A: LCDM (Omega_c=0.264, c_s^2=0)", bg_lcdm, k_arr, ell_values)

    # --- Model B: Khronon pure DBI (same background) ---
    bg_dbi = Background(khronon=False)  # SAME background
    bg_dbi.solve()
    Dl_B, sw_B, res_B = _run_model(
        f"B: Khronon DBI (Omega_c=0.264, alpha={alpha_dbi:.4f})",
        bg_dbi, k_arr, ell_values, alpha_dbi=alpha_dbi)

    # --- Model C: Khronon full (Omega_K=0.265 + DBI) ---
    bg_full = Background(khronon=True)  # Omega_K = 0.265
    bg_full.solve()
    Dl_C, sw_C, res_C = _run_model(
        f"C: Khronon full (Omega_K=0.265, alpha={alpha_dbi:.4f})",
        bg_full, k_arr, ell_values, alpha_dbi=alpha_dbi)

    # --- Residuals ---
    Dl_ref = np.maximum(Dl_A, 1e-10)
    res_BA = (Dl_B - Dl_A) / Dl_ref   # Pure DBI effect
    res_CA = (Dl_C - Dl_A) / Dl_ref   # Full Khronon effect
    res_CB = (Dl_C - Dl_B) / Dl_ref   # Background-only effect

    # c_s^2(k) at recombination
    cs2_k = res_B.effective_cs2_at_rec() if hasattr(res_B, 'effective_cs2_at_rec') \
        else np.zeros_like(k_arr)

    # --- Summary ---
    t_total = time.time() - t_start
    print("\n" + "=" * 70)
    print("3-WAY COMPARISON SUMMARY")
    print("=" * 70)
    print(f"Total time: {t_total:.1f}s")

    for label, res in [("B-A (pure DBI)", res_BA),
                        ("C-A (full Khronon)", res_CA),
                        ("C-B (bg only)", res_CB)]:
        print(f"\n  [{label}]:")
        for region, mask in [
            ("l < 500", ell_values < 500),
            ("500-2000", (ell_values >= 500) & (ell_values < 2000)),
            ("l >= 2000", ell_values >= 2000),
        ]:
            if np.any(mask):
                r = res[mask]
                print(f"    {region:10s}: max |res| = {np.max(np.abs(r)):.4e}, "
                      f"mean = {np.mean(np.abs(r)):.4e}")

    k_J_rec = np.sqrt(1.5 * _OMEGA_C * _H0_MPC**2 / a_rec)
    print(f"\n  Khronon Jeans: k_J(z_rec) = {k_J_rec:.4f} Mpc^-1, "
          f"l_J ~ {k_J_rec * bg_lcdm.D_A:.0f}")
    print(f"  Omega_K - Omega_c = "
          f"{bg_full.Omega_cdm - bg_lcdm.Omega_cdm:.4f} "
          f"({(bg_full.Omega_cdm - bg_lcdm.Omega_cdm)/bg_lcdm.Omega_cdm*100:.2f}%)")

    # --- Plot ---
    _plot_3way(ell_values, Dl_A, Dl_B, Dl_C, res_BA, res_CA, res_CB,
               k_arr, cs2_k, sw_A, sw_B, sw_C,
               alpha_dbi, t_total, save_dir)

    # --- Save ---
    data_path = os.path.join(save_dir, 'khronon_comparison.dat')
    np.savetxt(data_path,
               np.column_stack([ell_values, Dl_A, Dl_B, Dl_C,
                                res_BA, res_CA]),
               header=('ell  Dl_LCDM  Dl_Khr_DBI  Dl_Khr_full  '
                        'res_pureDBI  res_full'),
               fmt='%6d' + '  %.8e' * 5)
    print(f"\nData: {data_path}")

    return {
        'ell': ell_values,
        'Dl_lcdm': Dl_A,
        'Dl_khronon_dbi': Dl_B,
        'Dl_khronon_full': Dl_C,
        'residual_dbi': res_BA,
        'residual_full': res_CA,
        'residual_bg': res_CB,
        'cs2_k': cs2_k,
        'k_arr': k_arr,
        't_total': t_total,
    }


def _plot_3way(ell, Dl_A, Dl_B, Dl_C, res_BA, res_CA, res_CB,
               k_arr, cs2_k, sw_A, sw_B, sw_C,
               alpha_dbi, t_total, save_dir):
    """Generate the 4-panel 3-way comparison figure."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(4, 1, figsize=(13, 18),
                                  gridspec_kw={'height_ratios': [3, 1.5, 1.2, 1]})

        # ---- Panel 1: D_l overlay ----
        ax1 = axes[0]
        ax1.plot(ell, Dl_A, 'b-', lw=1.5, alpha=0.8,
                 label=r'A: $\Lambda$CDM ($\Omega_c=0.264$, $c_s^2=0$)')
        ax1.plot(ell, Dl_B, 'r--', lw=1.5, alpha=0.7,
                 label=r'B: Khronon DBI ($\Omega_c=0.264$, '
                       r'$\alpha_{\rm DBI}$=' + f'{alpha_dbi:.2f})')
        ax1.plot(ell, Dl_C, 'g:', lw=1.5, alpha=0.7,
                 label=r'C: Khronon full ($\Omega_K=0.265$ + DBI)')

        ax1.axvline(2000, color='gray', ls=':', alpha=0.5)
        ax1.set_ylabel(r'$D_\ell$ [$\mu K^2$]', fontsize=13)
        ax1.set_title(
            f'CMB TT: 3-Way Comparison  '
            r'($\alpha_{\rm DBI}$' + f' = {alpha_dbi:.2f}, '
            f'GPU {t_total:.1f}s)',
            fontsize=14)
        ax1.legend(fontsize=9, loc='upper right')
        ax1.set_xlim(2, ell[-1])
        if np.max(Dl_A) > 0:
            ax1.set_ylim(0, np.max(Dl_A) * 1.3)
        ax1.tick_params(labelsize=11)

        # ---- Panel 2: Residuals ----
        ax2 = axes[1]
        ax2.plot(ell, res_BA * 100, 'r-', lw=1.0, alpha=0.8,
                 label='B-A: pure DBI effect')
        ax2.plot(ell, res_CA * 100, 'g--', lw=1.0, alpha=0.8,
                 label='C-A: full Khronon')
        ax2.plot(ell, res_CB * 100, 'k:', lw=0.8, alpha=0.6,
                 label=r'C-B: $\Delta\Omega$ background only')

        ax2.axhline(0, color='gray', ls='-', alpha=0.3)
        ax2.axhline(1, color='red', ls=':', alpha=0.3, label='1% level')
        ax2.axhline(-1, color='red', ls=':', alpha=0.3)
        ax2.axvline(2000, color='gray', ls=':', alpha=0.5)

        ax2.set_ylabel(r'Residual [%]', fontsize=12)
        ax2.legend(fontsize=8, loc='best', ncol=2)
        ax2.set_xlim(2, ell[-1])

        all_res = np.concatenate([res_BA, res_CA, res_CB]) * 100
        valid = np.isfinite(all_res)
        if np.any(valid):
            ymax = min(max(np.max(np.abs(all_res[valid])) * 1.3, 2.0), 200.0)
            ax2.set_ylim(-ymax, ymax)
        ax2.tick_params(labelsize=11)

        # ---- Panel 3: Transfer function comparison ----
        ax3 = axes[2]
        ax3.plot(k_arr, sw_A, 'b-', lw=0.7, alpha=0.7, label=r'$\Lambda$CDM')
        ax3.plot(k_arr, sw_B, 'r--', lw=0.7, alpha=0.6, label='Khronon DBI')
        ax3.plot(k_arr, sw_C, 'g:', lw=0.7, alpha=0.6, label='Khronon full')
        ax3.set_xlabel(r'$k$ [Mpc$^{-1}$]', fontsize=12)
        ax3.set_ylabel(r'$(\Theta_0 + \Phi) \times$ Silk', fontsize=11)
        ax3.set_title('Transfer function at recombination', fontsize=11)
        ax3.legend(fontsize=9)
        ax3.set_xlim(k_arr[0], k_arr[-1])
        ax3.tick_params(labelsize=11)

        # ---- Panel 4: c_s^2(k) ----
        ax4 = axes[3]
        ax4.loglog(k_arr, cs2_k, 'g-', lw=1.5)
        ax4.axhline(alpha_dbi, color='orange', ls='--', alpha=0.6,
                     label=r'$\alpha_{\rm DBI}$ = ' + f'{alpha_dbi:.2f}')
        k_J_rec = np.sqrt(1.5 * _OMEGA_C * _H0_MPC**2 / a_rec)
        ax4.axvline(k_J_rec, color='purple', ls=':', alpha=0.6,
                     label=f'$k_J$ = {k_J_rec:.4f}')
        ax4.set_xlabel(r'$k$ [Mpc$^{-1}$]', fontsize=13)
        ax4.set_ylabel(r'$c_s^2(k)$ at $z_{\rm rec}$', fontsize=11)
        ax4.set_title('Khronon DBI sound speed', fontsize=11)
        ax4.legend(fontsize=9)
        ax4.set_xlim(k_arr[0], k_arr[-1])
        ax4.tick_params(labelsize=11)

        plt.tight_layout()
        out_path = os.path.join(save_dir, 'khronon_vs_lcdm.png')
        plt.savefig(out_path, dpi=150)
        print(f"Plot: {out_path}")
        plt.close()

    except ImportError:
        print("[Khronon] matplotlib not available, skipping plot")


# ============================================================================
# CLI
# ============================================================================

def main():
    import argparse
    parser = argparse.ArgumentParser(
        description='Khronon vs LCDM CMB comparison (MLX GPU)')
    parser.add_argument('--N_k', type=int, default=500)
    parser.add_argument('--ell_max', type=int, default=3000)
    parser.add_argument('--alpha_dbi', type=float, default=ALPHA_DBI_DEFAULT,
                        help=f'DBI strength (default: {ALPHA_DBI_DEFAULT})')
    args = parser.parse_args()
    run_comparison(N_k=args.N_k, ell_max=args.ell_max,
                   alpha_dbi=args.alpha_dbi)


if __name__ == '__main__':
    main()
