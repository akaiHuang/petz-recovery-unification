#!/usr/bin/env python3
"""
MLX GPU-Accelerated Black Hole Shadow Ray Tracer
=================================================
Compares Schwarzschild vs Exponential metric predictions for black hole shadows.

Physics
-------
Both metrics describe the same mass M with r_s = 2GM/c^2:

  Schwarzschild:  ds^2 = -(1 - r_s/r) dt^2 + dr^2/(1 - r_s/r) + r^2 dOmega^2
  Exponential:    ds^2 = -exp(-r_s/r) dt^2 + exp(+r_s/r) dr^2 + r^2 dOmega^2

The exponential metric arises from Sigma = -ln(-g_00) = r_s/r, where Sigma is
the quantum relative entropy between spacetime and matter channels (Paper 2).
Key differences from Schwarzschild:
  - No event horizon (exp(-r_s/r) > 0 for all r > 0)
  - Smaller photon sphere: r_ph = r_s/2  (vs 3r_s/2)
  - Smaller shadow:        b_c = e*r_s/2 (vs 3*sqrt(3)*r_s/2)
  - Smaller ISCO:          r_ISCO = r_s   (vs 3*r_s)
  - Shadow 47.7% smaller for the same mass

Orbit Equation
--------------
For null geodesics with u = 1/r:

  (du/dphi)^2 = h(u)

  Schwarzschild:  h(u) = 1/b^2 - u^2 + r_s u^3
  Exponential:    h(u) = 1/b^2 - u^2 exp(-r_s u)

Rewritten as second-order ODE:  d^2u/dphi^2 = (1/2) dh/du

  Schwarzschild:  d^2u/dphi^2 = -u + (3/2) r_s u^2
  Exponential:    d^2u/dphi^2 = (1/2) u exp(-r_s u)(r_s u - 2)

Exact analytic results:
  Schwarzschild: r_ph = 3r_s/2, b_c = 3*sqrt(3)*r_s/2 = 2.5981 r_s
  Exponential:   r_ph = r_s/2,  b_c = e*r_s/2           = 1.3591 r_s

Deflection angle for a ray starting at finite r_obs with impact parameter b:
  delta = 2*phi_turn - 2*arctan(sqrt(r_obs^2/b^2 - 1))

Reference
---------
S.-K. Huang, "Exponential metric from quantum relative entropy" (Paper 2).
Sigma = D(rho_spacetime || rho_matter) => g_00 = -exp(-Sigma).

Usage
-----
    python black_hole_shadow.py

Author: Sheng-Kai Huang / Generated with Claude Code
"""

import os
import time
import math
import mlx.core as mx
import numpy as np

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

R_S = 1.0       # Schwarzschild radius r_s = 2GM/c^2 = 1 (sets the scale)
R_OBS = 500.0   # Observer distance [r_s]

# ---------------------------------------------------------------------------
# 1. Analytic results
# ---------------------------------------------------------------------------

def analytic_results():
    """Return exact analytic values for both metrics."""
    results = {}

    # --- Schwarzschild ---
    rph = 1.5 * R_S
    bc = 1.5 * math.sqrt(3) * R_S       # 3*sqrt(3)/2 * r_s
    risco = 3.0 * R_S
    f_ph = 1.0 - R_S / rph              # 1/3
    z_ph = 1.0 / math.sqrt(f_ph) - 1.0  # sqrt(3) - 1
    sigma_ph = -math.log(f_ph)           # ln(3)
    results["schwarzschild"] = dict(
        r_ph=rph, b_crit=bc, r_isco=risco, f_at_rph=f_ph,
        redshift_at_rph=z_ph, sigma_at_rph=sigma_ph,
        has_horizon=True, r_horizon=R_S,
    )

    # --- Exponential ---
    rph = 0.5 * R_S
    bc = 0.5 * math.e * R_S             # e/2 * r_s
    risco = 1.0 * R_S
    f_ph = math.exp(-R_S / rph)          # exp(-2)
    z_ph = 1.0 / math.sqrt(f_ph) - 1.0  # e - 1
    sigma_ph = R_S / rph                 # 2
    results["exponential"] = dict(
        r_ph=rph, b_crit=bc, r_isco=risco, f_at_rph=f_ph,
        redshift_at_rph=z_ph, sigma_at_rph=sigma_ph,
        has_horizon=False, r_horizon=None,
    )

    results["shadow_ratio"] = results["exponential"]["b_crit"] / results["schwarzschild"]["b_crit"]
    results["shadow_diff_pct"] = (results["shadow_ratio"] - 1.0) * 100.0

    return results


# ---------------------------------------------------------------------------
# 2. MLX GPU-accelerated RK4 orbit integrator
# ---------------------------------------------------------------------------

def _integrate_batch(accel_fn, b_arr, u_capture, n_steps, phi_max):
    """
    Core RK4 integrator for d^2u/dphi^2 = accel(u).

    All N rays are integrated in parallel on the GPU.
    A ray is classified as:
      - CAPTURED if u > u_capture (fell into the compact object)
      - SCATTERED if it turns around and returns to large r (u < u_escape)
      - INDETERMINATE if neither condition met by end of integration
    """
    N = b_arr.shape[0]
    dphi = phi_max / n_steps

    u = mx.full((N,), 1.0 / R_OBS)
    b_inv2 = 1.0 / (b_arr * b_arr)
    v0_sq = b_inv2 - u * u  # At r_obs >> r_s, both metrics give h ~ 1/b^2 - u^2
    v = mx.sqrt(mx.maximum(v0_sq, mx.zeros((N,))))

    u_max = mx.array(u)
    has_turned = mx.zeros((N,), dtype=mx.bool_)
    captured = mx.zeros((N,), dtype=mx.bool_)
    escaped = mx.zeros((N,), dtype=mx.bool_)
    phi_at_turn = mx.full((N,), phi_max)

    u_escape = mx.array(1.0 / (0.8 * R_OBS))  # Ray has effectively escaped

    eval_every = 80

    for step in range(n_steps):
        # RK4
        a1 = accel_fn(u)
        u1h = u + 0.5 * dphi * v
        v1h = v + 0.5 * dphi * a1

        a2 = accel_fn(u1h)
        u2h = u + 0.5 * dphi * v1h
        v2h = v + 0.5 * dphi * a2

        a3 = accel_fn(u2h)
        u3 = u + dphi * v2h
        v3 = v + dphi * a3

        a4 = accel_fn(u3)

        u_new = u + (dphi / 6.0) * (v + 2.0 * v1h + 2.0 * v2h + v3)
        v_new = v + (dphi / 6.0) * (a1 + 2.0 * a2 + 2.0 * a3 + a4)

        # Only update active rays (not captured or escaped)
        active = ~(captured | escaped)
        u = mx.where(active, u_new, u)
        v = mx.where(active, v_new, v)
        u_max = mx.maximum(u_max, mx.where(active, u, u_max))

        # Turning point: v crosses from positive to negative
        newly_turned = active & (~has_turned) & (v < 0)
        phi_now = (step + 1) * dphi
        phi_at_turn = mx.where(newly_turned, mx.full((N,), phi_now), phi_at_turn)
        has_turned = has_turned | newly_turned

        # Periodic graph flush + capture/escape check
        if (step + 1) % eval_every == 0:
            mx.eval(u, v, u_max, has_turned, phi_at_turn, captured, escaped)
            # Capture: u exceeds threshold
            newly_captured = active & (u > u_capture)
            captured = captured | newly_captured
            # Escape: ray has turned and returned to large r
            newly_escaped = active & has_turned & (u < u_escape)
            escaped = escaped | newly_escaped
            mx.eval(captured, escaped)

    mx.eval(u, v, u_max, has_turned, captured, escaped, phi_at_turn)
    return phi_at_turn, captured, escaped, u_max


def _compute_deflection(b_arr, phi_at_turn, captured):
    """Compute deflection angle corrected for finite observer distance."""
    b2 = b_arr * b_arr
    robs2 = R_OBS * R_OBS
    ratio = mx.sqrt(mx.maximum(robs2 / b2 - 1.0, mx.zeros(b_arr.shape)))
    phi_flat = 2.0 * mx.arctan(ratio)
    deflection = 2.0 * phi_at_turn - phi_flat
    deflection = mx.where(captured, mx.full(b_arr.shape, float('nan')), deflection)
    return deflection


def integrate_schwarzschild(b_arr, n_steps=5000, phi_max=6.0 * np.pi):
    """Integrate Schwarzschild null geodesics."""
    def accel(uu):
        return -uu + 1.5 * R_S * uu * uu

    u_capture = 1.0 / (0.05 * R_S)  # r < 0.05 r_s
    phi_turn, captured, escaped, u_max = _integrate_batch(
        accel, b_arr, u_capture, n_steps, phi_max)

    return dict(
        phi_turn=phi_turn,
        deflection=_compute_deflection(b_arr, phi_turn, captured),
        captured=captured,
        u_max=u_max,
        r_min=1.0 / u_max,
    )


def integrate_exponential(b_arr, n_steps=8000, phi_max=10.0 * np.pi):
    """Integrate exponential metric null geodesics.

    Needs more steps/phi than Schwarzschild because the photon sphere is
    at r = r_s/2 (deeper), and rays near b_c spiral many times.
    """
    def accel(uu):
        return 0.5 * uu * mx.exp(-R_S * uu) * (R_S * uu - 2.0)

    u_capture = 10.0  # r < 0.1 r_s (well inside photon sphere at r_s/2)
    phi_turn, captured, escaped, u_max = _integrate_batch(
        accel, b_arr, u_capture, n_steps, phi_max)

    return dict(
        phi_turn=phi_turn,
        deflection=_compute_deflection(b_arr, phi_turn, captured),
        captured=captured,
        u_max=u_max,
        r_min=1.0 / u_max,
    )


# ---------------------------------------------------------------------------
# 3. Shadow radius finder (bisection)
# ---------------------------------------------------------------------------

def find_shadow_bisection(integrate_fn, b_lo, b_hi, n_bisect=35, **kwargs):
    """Bisection: b < b_c → captured, b > b_c → scattered."""
    lo, hi = b_lo, b_hi
    for _ in range(n_bisect):
        mid = 0.5 * (lo + hi)
        res = integrate_fn(mx.array([mid]), **kwargs)
        if bool(res["captured"][0]):
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


# ---------------------------------------------------------------------------
# 4. Plotting
# ---------------------------------------------------------------------------

def make_plots(b_s, rs_s, b_e, rs_e, ana, save_path):
    """Generate four-panel comparison plot."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec

    defl_s, cap_s, rmin_s = rs_s["deflection"], rs_s["captured"], rs_s["r_min"]
    defl_e, cap_e, rmin_e = rs_e["deflection"], rs_e["captured"], rs_e["r_min"]

    fig = plt.figure(figsize=(18, 14))
    gs = GridSpec(2, 2, hspace=0.32, wspace=0.28)

    bc_s = ana["schwarzschild"]["b_crit"]
    bc_e = ana["exponential"]["b_crit"]

    # --- Panel 1: Deflection angle ---
    ax1 = fig.add_subplot(gs[0, 0])
    ms = ~cap_s & np.isfinite(defl_s)
    me = ~cap_e & np.isfinite(defl_e)
    ax1.plot(b_s[ms], np.degrees(defl_s[ms]), 'b-', lw=1.0,
             label="Schwarzschild", alpha=0.8)
    ax1.plot(b_e[me], np.degrees(defl_e[me]), 'r-', lw=1.0,
             label="Exponential", alpha=0.8)
    ax1.axvline(bc_s, color='b', ls='--', lw=0.8,
                label=f'$b_c^{{Sch}}$ = {bc_s:.4f} $r_s$')
    ax1.axvline(bc_e, color='r', ls='--', lw=0.8,
                label=f'$b_c^{{Exp}}$ = {bc_e:.4f} $r_s$')
    ax1.set_xlabel("Impact parameter $b$ [$r_s$]", fontsize=12)
    ax1.set_ylabel("Deflection angle [deg]", fontsize=12)
    ax1.set_title("Photon Deflection Angle", fontsize=13, fontweight='bold')
    ax1.set_ylim(-5, 200)
    ax1.legend(fontsize=9, loc='upper right')
    ax1.grid(True, alpha=0.3)

    # --- Panel 2: Closest approach ---
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.plot(b_s, rmin_s, 'b-', lw=1.2, label="Schwarzschild", alpha=0.8)
    ax2.plot(b_e, rmin_e, 'r-', lw=1.2, label="Exponential", alpha=0.8)
    ax2.axhline(ana["schwarzschild"]["r_ph"], color='b', ls=':', lw=0.8,
                label=f'$r_{{ph}}^{{Sch}}$ = {ana["schwarzschild"]["r_ph"]:.2f} $r_s$')
    ax2.axhline(ana["exponential"]["r_ph"], color='r', ls=':', lw=0.8,
                label=f'$r_{{ph}}^{{Exp}}$ = {ana["exponential"]["r_ph"]:.2f} $r_s$')
    ax2.axhline(R_S, color='k', ls='--', lw=0.8, alpha=0.5, label='$r_s$ (Sch. horizon)')
    ax2.set_xlabel("Impact parameter $b$ [$r_s$]", fontsize=12)
    ax2.set_ylabel("Closest approach $r_{min}$ [$r_s$]", fontsize=12)
    ax2.set_title("Closest Approach Distance", fontsize=13, fontweight='bold')
    ax2.set_xlim(0, 6)
    ax2.set_ylim(0, 6)
    ax2.legend(fontsize=9, loc='upper left')
    ax2.grid(True, alpha=0.3)

    # --- Panel 3: Shadow image ---
    ax3 = fig.add_subplot(gs[1, 0])
    theta = np.linspace(0, 2 * np.pi, 500)
    ax3.fill(bc_s * np.cos(theta), bc_s * np.sin(theta),
             color='navy', alpha=0.25, label=f"Schwarzschild ($b_c$={bc_s:.3f})")
    ax3.fill(bc_e * np.cos(theta), bc_e * np.sin(theta),
             color='darkred', alpha=0.35, label=f"Exponential ($b_c$={bc_e:.3f})")
    ax3.plot(bc_s * np.cos(theta), bc_s * np.sin(theta), 'b-', lw=2)
    ax3.plot(bc_e * np.cos(theta), bc_e * np.sin(theta), 'r-', lw=2)
    ax3.set_xlim(-3.5, 3.5)
    ax3.set_ylim(-3.5, 3.5)
    ax3.set_aspect('equal')
    ax3.set_xlabel(r"$\alpha$ [$r_s$]", fontsize=12)
    ax3.set_ylabel(r"$\beta$ [$r_s$]", fontsize=12)
    ax3.set_title("Black Hole Shadow (Observer's Sky)", fontsize=13, fontweight='bold')
    ax3.legend(fontsize=9, loc='upper right')
    ax3.grid(True, alpha=0.3)

    # --- Panel 4: Effective potential ---
    ax4 = fig.add_subplot(gs[1, 1])
    r_arr = np.linspace(0.15, 8.0, 2000)
    f_sch = np.where(r_arr > R_S, 1.0 - R_S / r_arr, np.nan)
    V_sch = f_sch / r_arr**2
    V_exp = np.exp(-R_S / r_arr) / r_arr**2

    ax4.plot(r_arr, V_sch, 'b-', lw=1.5, label=r"$V_{eff}^{Sch}$")
    ax4.plot(r_arr, V_exp, 'r-', lw=1.5, label=r"$V_{eff}^{Exp}$")

    V_sch_pk = (1.0 / 3.0) / (1.5)**2
    V_exp_pk = np.exp(-2.0) / (0.5)**2
    ax4.plot(1.5, V_sch_pk, 'bo', ms=8, zorder=5)
    ax4.plot(0.5, V_exp_pk, 'ro', ms=8, zorder=5)
    ax4.annotate(f"$r_{{ph}}$ = 1.5 $r_s$", xy=(1.5, V_sch_pk),
                 xytext=(3.0, V_sch_pk + 0.01), fontsize=9, color='blue',
                 arrowprops=dict(arrowstyle='->', color='blue', lw=0.8))
    ax4.annotate(f"$r_{{ph}}$ = 0.5 $r_s$", xy=(0.5, V_exp_pk),
                 xytext=(2.0, V_exp_pk + 0.05), fontsize=9, color='red',
                 arrowprops=dict(arrowstyle='->', color='red', lw=0.8))

    ax4.axhline(1 / bc_s**2, color='b', ls=':', lw=0.6, alpha=0.5,
                label=f'$1/b_c^{{2,Sch}}$ = {1/bc_s**2:.4f}')
    ax4.axhline(1 / bc_e**2, color='r', ls=':', lw=0.6, alpha=0.5,
                label=f'$1/b_c^{{2,Exp}}$ = {1/bc_e**2:.4f}')

    ax4.axvline(R_S, color='k', ls='--', lw=0.8, alpha=0.4, label='$r_s$')
    ax4.set_xlabel("$r$ [$r_s$]", fontsize=12)
    ax4.set_ylabel(r"$V_{eff}(r) = f(r)/r^2$ [$r_s^{-2}$]", fontsize=12)
    ax4.set_title("Photon Effective Potential", fontsize=13, fontweight='bold')
    ax4.set_ylim(-0.02, 0.65)
    ax4.legend(fontsize=8, loc='upper right')
    ax4.grid(True, alpha=0.3)

    fig.suptitle(
        "Black Hole Shadow: Schwarzschild vs Exponential Metric\n"
        r"$\Sigma = -\ln(-g_{00}) = r_s/r$"
        f"  |  Shadow difference: {ana['shadow_diff_pct']:.1f}%",
        fontsize=15, fontweight='bold', y=0.98
    )

    plt.savefig(save_path, dpi=180, bbox_inches='tight', facecolor='white')
    print(f"\nPlot saved to: {save_path}")
    plt.close()


# ---------------------------------------------------------------------------
# 5. Main
# ---------------------------------------------------------------------------

def main():
    t_start = time.time()

    print("=" * 72)
    print("  BLACK HOLE SHADOW RAY TRACER")
    print("  Schwarzschild vs Exponential Metric (MLX GPU-accelerated)")
    print("=" * 72)
    print()

    ana = analytic_results()
    s = ana["schwarzschild"]
    e = ana["exponential"]

    # ---- Analytic table ----
    print("ANALYTIC RESULTS (exact)")
    print("-" * 72)
    print(f"{'Quantity':<35} {'Schwarzschild':>16} {'Exponential':>16}")
    print("-" * 72)
    rows = [
        ("Photon sphere r_ph [r_s]",     f"{s['r_ph']:.4f}",           f"{e['r_ph']:.4f}"),
        ("Critical impact b_c [r_s]",    f"{s['b_crit']:.4f}",         f"{e['b_crit']:.4f}"),
        ("  b_c formula",                "3*sqrt(3)/2",                "e/2"),
        ("ISCO radius r_ISCO [r_s]",     f"{s['r_isco']:.4f}",         f"{e['r_isco']:.4f}"),
        ("Event horizon",               "r_H = r_s",                  "None"),
        ("f(r_ph) = -g_00(r_ph)",        f"{s['f_at_rph']:.6f}",       f"{e['f_at_rph']:.6f}"),
        ("  f formula at r_ph",          "1/3",                        "exp(-2)"),
        ("Redshift z at r_ph",           f"{s['redshift_at_rph']:.4f}", f"{e['redshift_at_rph']:.4f}"),
        ("Sigma at r_ph",               f"{s['sigma_at_rph']:.4f}",    f"{e['sigma_at_rph']:.4f}"),
        ("  Sigma formula",             "ln(3)",                       "2"),
    ]
    for label, v1, v2 in rows:
        print(f"{label:<35} {v1:>16} {v2:>16}")
    print("-" * 72)
    print(f"{'Shadow ratio b_c(Exp)/b_c(Sch)':<35} {ana['shadow_ratio']:>16.4f}")
    print(f"{'Shadow difference':<35} {ana['shadow_diff_pct']:>15.1f}%")
    print()

    # ---- Sigma table ----
    print("SIGMA = -ln(-g_00) AT KEY RADII")
    print("-" * 72)
    print(f"{'Radius':<25} {'Schwarzschild':>20} {'Exponential':>20}")
    print("-" * 72)
    for label, rv in [("r = 10 r_s (weak)", 10.0),
                      ("r = 3 r_s (Sch ISCO)", 3.0),
                      ("r = 1.5 r_s (Sch r_ph)", 1.5),
                      ("r = r_s (Sch horizon)", 1.0),
                      ("r = 0.5 r_s (Exp r_ph)", 0.5)]:
        f_val = 1.0 - R_S / rv
        if f_val > 0:
            sig_s = f"{-np.log(f_val):.6f}"
        elif f_val == 0:
            sig_s = "infinity"
        else:
            sig_s = "undefined (r<r_H)"
        print(f"{label:<25} {sig_s:>20} {R_S / rv:>20.6f}")
    print("-" * 72)
    print()

    # ---- GPU ray tracing ----
    print("MLX GPU RAY TRACING")
    print(f"Device: {mx.default_device()}")
    print(f"Observer: r_obs = {R_OBS:.0f} r_s")
    print("-" * 72)

    # ---- Step 1: Shadow radius (bisection) ----
    print("\n[1] Shadow radius by bisection (35 iterations)...")

    t0 = time.time()
    bc_s_num = find_shadow_bisection(
        integrate_schwarzschild, 0.5, 5.0, n_bisect=35, n_steps=8000)
    dt = time.time() - t0
    err_s = abs(bc_s_num - s['b_crit']) / s['b_crit'] * 100
    print(f"  Schwarzschild: b_c = {bc_s_num:.6f}  "
          f"(exact: {s['b_crit']:.6f})  err = {err_s:.3f}%  [{dt:.1f}s]")

    t0 = time.time()
    bc_e_num = find_shadow_bisection(
        integrate_exponential, 0.05, 3.0, n_bisect=35,
        n_steps=15000, phi_max=12.0 * np.pi)
    dt = time.time() - t0
    err_e = abs(bc_e_num - e['b_crit']) / e['b_crit'] * 100
    print(f"  Exponential:   b_c = {bc_e_num:.6f}  "
          f"(exact: {e['b_crit']:.6f})  err = {err_e:.3f}%  [{dt:.1f}s]")

    # ---- Step 2: Full deflection profiles ----
    N_RAYS = 10000
    print(f"\n[2] Deflection profiles ({N_RAYS} rays each)...")

    t0 = time.time()
    b_s_mx = mx.linspace(0.1, 10.0, N_RAYS)
    res_s = integrate_schwarzschild(b_s_mx, n_steps=6000, phi_max=8.0 * np.pi)
    dt = time.time() - t0
    n_cap_s = int(mx.sum(res_s["captured"]).item())
    print(f"  Schwarzschild: {n_cap_s} captured, "
          f"{N_RAYS - n_cap_s} scattered  [{dt:.1f}s]")

    t0 = time.time()
    b_e_mx = mx.linspace(0.02, 10.0, N_RAYS)
    res_e = integrate_exponential(b_e_mx, n_steps=10000, phi_max=12.0 * np.pi)
    dt = time.time() - t0
    n_cap_e = int(mx.sum(res_e["captured"]).item())
    print(f"  Exponential:   {n_cap_e} captured, "
          f"{N_RAYS - n_cap_e} scattered  [{dt:.1f}s]")

    # ---- Step 3: Weak-field check ----
    print("\n[3] Weak-field deflection (b = 50 r_s)...")
    b_wk = mx.array([50.0])
    d_s = float(integrate_schwarzschild(b_wk, n_steps=8000, phi_max=4*np.pi)["deflection"][0])
    d_e = float(integrate_exponential(b_wk, n_steps=8000, phi_max=4*np.pi)["deflection"][0])
    d_th = 2.0 * R_S / 50.0  # Leading-order GR: delta = 2r_s/b
    print(f"  Schwarzschild: {np.degrees(d_s):.4f} deg")
    print(f"  Exponential:   {np.degrees(d_e):.4f} deg")
    print(f"  GR 1st order:  {np.degrees(d_th):.4f} deg")
    # 2nd-order correction for Schwarzschild: delta = 2r_s/b + (15pi/16)(r_s/b)^2
    d_th2_sch = 2.0 * R_S / 50.0 + (15 * np.pi / 16) * (R_S / 50.0)**2
    print(f"  GR 2nd order (Sch): {np.degrees(d_th2_sch):.4f} deg")

    # ---- Convert to numpy ----
    b_s_np = np.array(b_s_mx.tolist())
    b_e_np = np.array(b_e_mx.tolist())

    def to_np(res):
        return {
            "deflection": np.array(res["deflection"].tolist()),
            "captured": np.array(res["captured"].tolist()).astype(bool),
            "r_min": np.array(res["r_min"].tolist()),
        }

    rs_s_np = to_np(res_s)
    rs_e_np = to_np(res_e)

    # ---- Numerical shadow from ray count ----
    # The shadow radius can also be estimated from the capture fraction
    b_s_cap = b_s_np[rs_s_np["captured"]]
    b_e_cap = b_e_np[rs_e_np["captured"]]
    bc_s_from_rays = np.max(b_s_cap) if len(b_s_cap) > 0 else 0
    bc_e_from_rays = np.max(b_e_cap) if len(b_e_cap) > 0 else 0

    # ---- Summary ----
    print()
    print("=" * 72)
    print("  SUMMARY OF RESULTS")
    print("=" * 72)
    print()
    print(f"{'Quantity':<42} {'Schwarzschild':>13} {'Exponential':>13}")
    print("-" * 72)
    print(f"{'Photon sphere r_ph / r_s':<42} {s['r_ph']:>13.4f} {e['r_ph']:>13.4f}")
    print(f"{'Shadow b_c / r_s (analytic)':<42} {s['b_crit']:>13.4f} {e['b_crit']:>13.4f}")
    print(f"{'Shadow b_c / r_s (bisection)':<42} {bc_s_num:>13.6f} {bc_e_num:>13.6f}")
    print(f"{'Shadow b_c / r_s (from rays)':<42} {bc_s_from_rays:>13.4f} {bc_e_from_rays:>13.4f}")
    print(f"{'ISCO r_ISCO / r_s':<42} {s['r_isco']:>13.4f} {e['r_isco']:>13.4f}")
    print(f"{'Event horizon':<42} {'r_s':>13} {'None':>13}")
    print(f"{'Redshift at r_ph':<42} {s['redshift_at_rph']:>13.4f} {e['redshift_at_rph']:>13.4f}")
    print(f"{'Sigma at r_ph':<42} {s['sigma_at_rph']:>13.4f} {e['sigma_at_rph']:>13.4f}")
    cap_label = f"Captured rays (of {N_RAYS})"
    print(f"{cap_label:<42} {n_cap_s:>13d} {n_cap_e:>13d}")
    print("-" * 72)
    print(f"{'Shadow ratio (Exp/Sch)':<42} {ana['shadow_ratio']:>13.4f}")
    print(f"{'Shadow difference':<42} {ana['shadow_diff_pct']:>12.1f}%")
    print()
    print("PHYSICAL INTERPRETATION:")
    print(f"  The exponential metric predicts a shadow {abs(ana['shadow_diff_pct']):.1f}%")
    print(f"  SMALLER than Schwarzschild for the same mass M.")
    print()
    print(f"  Schwarzschild has an event horizon at r = r_s where Sigma diverges.")
    print(f"  The exponential metric has NO horizon: Sigma = r_s/r is finite")
    print(f"  everywhere. Photons can orbit at r = r_s/2, deep inside the")
    print(f"  would-be Schwarzschild horizon.")
    print()
    print(f"  At the photon sphere:")
    print(f"    Schwarzschild: Sigma = ln(3) = {s['sigma_at_rph']:.4f}")
    print(f"    Exponential:   Sigma = 2")
    print()
    print(f"  For M87* / Sgr A* (EHT): angular shadow size = b_c / D_L")
    print(f"  The {abs(ana['shadow_diff_pct']):.0f}% difference should be resolvable")
    print(f"  by ngEHT (targeting few-% shadow size precision).")
    print()

    # ---- Plot ----
    script_dir = os.path.dirname(os.path.abspath(__file__))
    plot_path = os.path.join(script_dir, "black_hole_shadow.png")
    print("Generating comparison plot...")
    make_plots(b_s_np, rs_s_np, b_e_np, rs_e_np, ana, plot_path)

    print(f"\nTotal runtime: {time.time() - t_start:.1f}s")
    print("Done.")


if __name__ == "__main__":
    main()
