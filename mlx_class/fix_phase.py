"""
fix_phase.py — Diagnose and fix the photon oscillation phase error.

Findings from the 2026-03-22 diagnostic session:

ROOT CAUSES IDENTIFIED:
1. Omega_gamma rounding: 2.469e-5 vs exact 2.4730e-5 (0.16% error) -- FIXED
2. Conformal time integration: left-endpoint vs trapezoidal -- FIXED
3. Neutrino-to-photon ratio: 0.2271 approximation vs exact 7/8*(4/11)^{4/3} -- FIXED
4. Anisotropic stress formula (Psi computation): factor needs verification
5. Remaining ~5% CDM deficit at k=0.1, growing with k -- UNDER INVESTIGATION

FIXES APPLIED:
- background.py: Omega_gamma now computed from T_CMB via Stefan-Boltzmann (exact)
- background.py: Conformal time uses trapezoidal rule (2nd order)
- background.py: R = 3*Omega_b*a/(4*Omega_gamma) uses exact neutrino ratio
- perturbations_neutrino.py: f_nu uses exact 7/8*(4/11)^{4/3} ratio

IMPACT:
- Sound horizon: 144.39 Mpc (was 144.43), CLASS gives 144.33 -- improved to 0.04%
- Conformal time at z=0: matches CLASS to 0.0003% (was 0.012%)
- Transfer functions: CDM deficit unchanged (~5% at k=0.1)
  - This 5% deficit is NOT from Omega_gamma or integration order
  - It's a cumulative effect from ~0.1-0.2% background discrepancies
  - At z=5000 it's already 4% at k=0.12, growing to 6% by z=1100

REMAINING WORK:
- The 5% CDM deficit likely affects the C_l through the metric source in the
  LOS integral, but the DIRECT photon oscillation (delta_gamma) matches CLASS
  to ~3% for the first few oscillation cycles
- The second peak C_l deficit (0.63x) may be dominated by:
  (a) The Psi_N formula (anisotropic stress factor)
  (b) The LOS integration grid resolution
  (c) The photon phase shift accumulation (~0.5-0.9% in zero crossings)

Usage:
    python -m mlx_class.fix_phase
"""
import os
os.environ['MPLCONFIGDIR'] = '/tmp/mpl_config'

import numpy as np
import time

from .background import Background, H0_Mpc as _H0_MPC, h, omega_b, T_CMB
from .perturbations_neutrino import _OMEGA_GAMMA, _OMEGA_NU
from .perturbations_sync import (
    make_sync_rhs, adiabatic_ic_sync,
    IDX_ETA, IDX_DELTA_C, IDX_DELTA_B, IDX_THETA_B, IDX_FG_START,
    idx_fn_start, L_GAMMA_MAX, L_NU_MAX_SYNC,
)


def main():
    from scipy.integrate import solve_ivp

    t0 = time.time()
    print("=" * 70)
    print("PHOTON PHASE ERROR — DIAGNOSTIC SUMMARY")
    print("=" * 70)

    # Background
    bg = Background(khronon=False, recombination='recfast')
    bg.solve()

    print(f"\nOmega_gamma = {_OMEGA_GAMMA:.10e}")
    print(f"Omega_nu    = {_OMEGA_NU:.10e}")
    print(f"r_s = {bg.r_s:.4f} Mpc")
    print(f"R_rec = {bg.R_rec:.6f}")

    # Load CLASS reference
    class_path = '/Users/akaihuangm1/Desktop/github/gdm_class_public/output/mlx_match_tk.dat'
    try:
        data = np.loadtxt(class_path)
    except FileNotFoundError:
        print(f"\n[WARNING] CLASS reference not found at {class_path}")
        print("Run CLASS with mlx_match.ini to generate it.")
        return

    k_class = data[:, 0] * h
    d_cdm_class = data[:, 3]
    d_g_class = data[:, 1]

    # Transfer function comparison
    tau_init = bg.tau_grid[1]
    lg_max, ln_max = L_GAMMA_MAX, L_NU_MAX_SYNC
    a_eval = 1.0 / 1101.0
    tau_eval = float(bg._tau_of_a(np.log(a_eval)))

    mask = (k_class > 0.005) & (k_class < 0.15)
    k_test = k_class[mask]

    print(f"\nTransfer function comparison at z=1100 ({len(k_test)} modes):")
    print(f"{'k':>10} | {'d_c/CLASS':>10} {'d_g/CLASS':>10}")
    print("-" * 38)

    ratios_dc = []
    ratios_dg = []
    for ik in range(0, len(k_test), 3):  # every 3rd for speed
        k = float(k_test[ik])
        rhs_fn, nvar = make_sync_rhs(k, bg, lg_max, ln_max)
        y0 = adiabatic_ic_sync(k, tau_init, lg_max, ln_max)
        sol = solve_ivp(rhs_fn, [tau_init, tau_eval + 5], y0,
                        method='Radau', dense_output=True, rtol=1e-8, atol=1e-11)
        if not sol.success:
            continue
        y = sol.sol(tau_eval)

        idx = np.argmin(np.abs(k_class - k))
        r_dc = y[IDX_DELTA_C] / d_cdm_class[idx]
        r_dg = y[IDX_FG_START] / d_g_class[idx] if abs(d_g_class[idx]) > 1e-10 else float('nan')
        ratios_dc.append(r_dc)
        if np.isfinite(r_dg) and 0.5 < abs(r_dg) < 2.0:
            ratios_dg.append(r_dg)
        print(f"{k:10.6f} | {r_dc:10.4f} {r_dg:10.4f}")

    ratios_dc = np.array(ratios_dc)
    ratios_dg = np.array(ratios_dg)

    print(f"\nSUMMARY:")
    print(f"  delta_cdm/CLASS: mean={np.mean(ratios_dc):.4f}, range=[{min(ratios_dc):.4f}, {max(ratios_dc):.4f}]")
    if len(ratios_dg) > 0:
        print(f"  delta_gamma/CLASS (away from zeros): mean={np.mean(ratios_dg):.4f}")

    print(f"\nDiagnostic time: {time.time()-t0:.1f}s")


if __name__ == '__main__':
    main()
