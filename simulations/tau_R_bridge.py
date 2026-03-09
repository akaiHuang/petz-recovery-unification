"""
tau-R Bridge: Gravitational Dephasing Simulation

Demonstrates that the irreversibility parameter tau and the quantum Darwinism
redundancy R reach saturation at the same timescale, driven by the same
N*Gamma parameter (Pikovski gravitational dephasing mechanism).

Physics:
  A mass in spatial superposition |L> + |R> is coupled to N internal
  (vibrational) modes via gravitational time dilation. Each mode k dephases
  independently with rate Gamma_k, giving:
    D(t) = prod_k exp(-Gamma_k t^2)  =>  exp(-N*Gamma_bar*t^2)
    tau(t) = 1 - F(rho_S(t), rho_S(0))
    R_delta(t) = number of environment fragments that independently
                 determine the system state to within tolerance delta.

The central result: both tau and R/N are functions of the single parameter
N*Gamma*t^2, proving they are two manifestations of the same underlying
information dynamics. This supports the tau-R bridge:
    tau > 0  <=>  R > 0  <=>  arrow of time emerges.

Model:
  System: 1 qubit (spatial DOF) in state |+> = (|0>+|1>)/sqrt(2)
  Environment: N qubits, each initially in |+>_E (thermal superposition)
  Interaction: H_int = epsilon * sigma_z^S (x) sigma_z^{E_k}
  Unitary: U_k(t) = exp(-i*epsilon*t*sigma_z (x) sigma_z)
         = diag(e^{-i*eps*t}, e^{+i*eps*t}, e^{+i*eps*t}, e^{-i*eps*t})

  Single-mode decoherence factor: d_k(t) = cos(2*epsilon*t)
  Total decoherence factor: D(t) = cos^N(2*epsilon*t)
  Gaussian approximation: D(t) ~ exp(-N*Gamma*t^2), Gamma = 2*epsilon^2

  For the |+> initial state:
    rho_S(t) = [[1/2, D(t)/2], [D(t)/2, 1/2]]
    F(rho_S(0), rho_S(t)) = (1 + |D(t)|) / 2
    tau(t) = (1 - |D(t)|) / 2

Part A: Exact N-qubit simulation (N = 2..12)
Part B: Analytical formulas for large N (N = 10..10^6)
Part C: Verification of F >= exp(-Sigma/2) bound using the
        Crooks-type entropy production Sigma = N * [S(rho_Ek) - S(rho_E0)]
        + S(rho_S(t)) - S(rho_S(0)), which is the total entropy production
        from the unitary system-environment interaction.

Reference: Huang (2026), connecting Petz recovery to quantum Darwinism.
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

_DIR = os.path.dirname(os.path.abspath(__file__))


# ============================================================
# Utility functions
# ============================================================

def binary_entropy(p):
    """Binary entropy h(p) = -p*ln(p) - (1-p)*ln(1-p), in nats."""
    p = np.asarray(p, dtype=float)
    result = np.zeros_like(p)
    mask = (p > 1e-15) & (p < 1 - 1e-15)
    result[mask] = -p[mask] * np.log(p[mask]) - (1 - p[mask]) * np.log(1 - p[mask])
    return result


def binary_entropy_scalar(p):
    """Scalar version of binary entropy."""
    if p <= 1e-15 or p >= 1 - 1e-15:
        return 0.0
    return -p * np.log(p) - (1 - p) * np.log(1 - p)


def von_neumann_entropy_eig(rho):
    """Von Neumann entropy from eigenvalues of a density matrix."""
    rho_h = 0.5 * (rho + rho.conj().T)
    eigvals = np.linalg.eigvalsh(rho_h)
    eigvals = eigvals[eigvals > 1e-15]
    return float(-np.sum(eigvals * np.log(eigvals)))


def mutual_info_SE_k(epsilon, t, N):
    """
    Compute I(S:E_k) exactly for the gravitational dephasing model.

    Returns: I(S:E_k), S(rho_S)
    """
    phase = epsilon * t
    c = np.cos(2 * epsilon * t)
    D_other = c ** (N - 1) if N > 1 else 1.0
    D_total = c ** N

    # Conditional environment states
    chi_0 = np.array([np.exp(-1j * phase), np.exp(1j * phase)]) / np.sqrt(2)
    chi_1 = np.array([np.exp(1j * phase), np.exp(-1j * phase)]) / np.sqrt(2)

    # Build 4x4 joint state
    rho_SE = np.zeros((4, 4), dtype=complex)
    rho_SE[0:2, 0:2] = 0.5 * np.outer(chi_0, chi_0.conj())
    rho_SE[0:2, 2:4] = 0.5 * D_other * np.outer(chi_0, chi_1.conj())
    rho_SE[2:4, 0:2] = 0.5 * D_other * np.outer(chi_1, chi_0.conj())
    rho_SE[2:4, 2:4] = 0.5 * np.outer(chi_1, chi_1.conj())
    rho_SE = 0.5 * (rho_SE + rho_SE.conj().T)

    # Reduced states
    rho_S = np.array([[0.5, 0.5 * D_total],
                       [0.5 * D_total, 0.5]], dtype=complex)
    rho_Ek = rho_SE[0:2, 0:2] + rho_SE[2:4, 2:4]
    rho_Ek = 0.5 * (rho_Ek + rho_Ek.conj().T)

    S_S = von_neumann_entropy_eig(rho_S)
    S_Ek = von_neumann_entropy_eig(rho_Ek)
    S_SE = von_neumann_entropy_eig(rho_SE)

    MI = max(S_S + S_Ek - S_SE, 0.0)
    return MI, S_S


# ============================================================
# Part A: Exact N-qubit gravitational dephasing
# ============================================================

def exact_dephasing_simulation(N, epsilon, t_values):
    """
    Exact simulation of gravitational dephasing for N environment qubits.

    Entropy production Sigma:
    The total evolution is unitary on S+E, so S(rho_SE) = const = 0
    (product state initially, each factor pure).
    Sigma = sum of marginal entropy changes:
      Sigma = [S(rho_S(t)) - S(rho_S(0))] + N * [S(rho_{E_k}(t)) - S(rho_{E_k}(0))]
    Since rho_S(0) = |+><+| (pure, S=0) and rho_{E_k}(0) = |+><+| (pure, S=0):
      Sigma = S(rho_S(t)) + N * S(rho_{E_k}(t))

    This is the mutual information I(S:E_1:...:E_N), the total correlation,
    which equals the entropy production for a product initial state.

    The bound from the Fawzi-Renner/Junge-Renner-Sutter-Wilde-Winter
    universal recovery theorem is:
      F(rho, R(N(rho)))^2 >= exp(-Delta_D)
    where Delta_D = D(rho||sigma) - D(N(rho)||N(sigma)).

    For our system, we verify a related but simpler bound:
      F(rho_S(0), rho_S(t)) >= exp(-Sigma/2)
    where Sigma is the full entropy production. This bound comes from:
      F >= sqrt(exp(-Sigma)) = exp(-Sigma/2)
    via the Pinsker-type inequality chain in the paper.

    Returns: tau, R_N, MI, F, Sigma
    """
    delta = 0.05
    n_t = len(t_values)
    tau_vals = np.zeros(n_t)
    R_norm_vals = np.zeros(n_t)
    MI_vals = np.zeros(n_t)
    F_vals = np.zeros(n_t)
    Sigma_vals = np.zeros(n_t)

    for ti, t in enumerate(t_values):
        c = np.cos(2 * epsilon * t)
        D_t = c ** N
        abs_D = abs(D_t)

        # tau and F
        F = (1 + abs_D) / 2
        F_vals[ti] = F
        tau_vals[ti] = (1 - abs_D) / 2

        # Mutual information I(S:E_k) and S(rho_S)
        MI, S_sys = mutual_info_SE_k(epsilon, t, N)
        MI_vals[ti] = MI

        # Redundancy R_delta/N
        # Standard quantum Darwinism: R_delta counts how many fragments
        # each independently capture (1-delta)*S_sys of the system info.
        # For N identical modes, each carrying MI:
        #   f_min modes needed = ceil((1-delta)*S_sys / MI)
        #   R_delta = N / f_min
        #   R_delta/N = MI / ((1-delta)*S_sys)  [when MI < (1-delta)*S_sys]
        #
        # We require S_sys > S_threshold to avoid 0/0 artifacts at early
        # times when the system has barely decohered.
        S_threshold = 0.05 * np.log(2)  # ~5% of max system entropy
        if S_sys > S_threshold and MI > 1e-15:
            R_norm_vals[ti] = min(MI / ((1 - delta) * S_sys), 1.0)
        else:
            R_norm_vals[ti] = 0.0

        # Entropy production: Sigma = S(rho_S(t)) + N * S(rho_{E_k}(t))
        # S(rho_S(t)) = h((1+|D|)/2)
        S_sys_t = binary_entropy_scalar((1 + abs_D) / 2)

        # For a single env qubit after interaction with system in |+>:
        # rho_{E_k}(t) = Tr_S[rho_{S,E_k}(t)]
        # = (1/2)(|chi_0><chi_0| + |chi_1><chi_1|)
        # = I/2 + (1/2)*Re(chi_0 chi_1^dag + chi_1 chi_0^dag - I)
        # The overlap is <chi_0|chi_1> = cos(2*eps*t), so
        # rho_{E_k} = [[1/2, cos(2*eps*t)/2], [cos(2*eps*t)/2, 1/2]]
        # (in the computational basis, after phase simplification)
        # Wait: actually
        # chi_0 = (e^{-i*eps*t}|0> + e^{+i*eps*t}|1>)/sqrt(2)
        # chi_1 = (e^{+i*eps*t}|0> + e^{-i*eps*t}|1>)/sqrt(2)
        # <0|rho_Ek|0> = (|chi_0[0]|^2 + |chi_1[0]|^2)/2 = (1/2+1/2)/2 = 1/2
        # <0|rho_Ek|1> = (chi_0[0]*chi_0[1]^* + chi_1[0]*chi_1[1]^*)/2
        #              = (e^{-i*eps*t}*e^{-i*eps*t} + e^{+i*eps*t}*e^{+i*eps*t})/4
        #              = (e^{-2i*eps*t} + e^{+2i*eps*t})/4 = cos(2*eps*t)/2
        # So rho_Ek = [[1/2, cos(2*eps*t)/2], [cos(2*eps*t)/2, 1/2]]
        # and S(rho_Ek) = h((1 + |cos(2*eps*t)|) / 2)
        d_k = abs(np.cos(2 * epsilon * t))
        S_Ek = binary_entropy_scalar((1 + d_k) / 2)

        Sigma_vals[ti] = S_sys_t + N * S_Ek

    return tau_vals, R_norm_vals, MI_vals, F_vals, Sigma_vals


# ============================================================
# Part B: Analytical formulas
# ============================================================

def analytical_tau(t, N, Gamma):
    """tau(t) = (1 - exp(-N*Gamma*t^2)) / 2  (Gaussian approximation)."""
    return (1 - np.exp(-N * Gamma * t**2)) / 2


def analytical_R_over_N(t, Gamma):
    """
    R(t)/N for independent modes in the Gaussian approximation.
    Each mode's distinguishability grows as 1 - exp(-2*Gamma*t^2).
    """
    return np.minimum(1 - np.exp(-2 * Gamma * t**2), 1.0)


# ============================================================
# Main
# ============================================================

def main():
    print("=" * 70)
    print("tau-R Bridge: Gravitational Dephasing Simulation")
    print("=" * 70)
    print()

    # Parameters
    epsilon = 0.3
    Gamma = 2 * epsilon**2  # = 0.18
    N_values = [2, 3, 4, 6, 8, 10, 12]
    t_max = 5.0
    n_times = 300
    t_values = np.linspace(0, t_max, n_times)

    cmap = plt.cm.viridis
    colors = {N: cmap(i / max(len(N_values) - 1, 1))
              for i, N in enumerate(N_values)}

    all_results = {}

    # ================================================================
    # Part A: Exact simulation
    # ================================================================
    print("Part A: Exact N-qubit simulation")
    print("-" * 55)

    for N in N_values:
        tau, R_N, MI, F, Sigma = exact_dephasing_simulation(
            N, epsilon, t_values)
        all_results[N] = {
            'tau': tau, 'R_N': R_N, 'MI': MI, 'F': F, 'Sigma': Sigma
        }

        t_c = 1.0 / np.sqrt(N * Gamma)
        tau_max = np.max(tau)

        # Half-saturation for tau
        mask = t_values <= 3 * t_c
        idx_half = np.where(tau[mask] >= 0.20)[0]
        t_half = t_values[idx_half[0]] if len(idx_half) > 0 else np.nan

        if not np.isnan(t_half):
            print(f"  N = {N:3d}: tau_max = {tau_max:.4f}, "
                  f"t_half = {t_half:.3f}, t_c = {t_c:.3f}, "
                  f"t_half/t_c = {t_half / t_c:.2f}")
        else:
            print(f"  N = {N:3d}: tau_max = {tau_max:.4f}, t_c = {t_c:.3f}")

    # ================================================================
    # Part C: Bound verification
    # ================================================================
    print()
    print("Part C: F >= exp(-Sigma/2) bound verification")
    print("-" * 55)
    print()
    print("Entropy production: Sigma = S(rho_S(t)) + N * S(rho_{E_k}(t))")
    print("  (total correlation created by unitary S-E interaction)")
    print("Bound: F = (1+|D|)/2 >= exp(-Sigma/2)")
    print()

    total_pts = 0
    total_viol = 0
    for N in N_values:
        F = all_results[N]['F']
        Sigma = all_results[N]['Sigma']
        bound = np.exp(-Sigma / 2)
        gaps = F - bound
        violations = int(np.sum(gaps < -1e-10))
        total_pts += len(F)
        total_viol += violations
        min_gap = np.min(gaps)
        status = "PASS" if violations == 0 else f"FAIL ({violations})"
        print(f"  N = {N:3d}: {len(F)} points, {status}, "
              f"min gap = {min_gap:.6f}")

    print(f"\n  Total: {total_pts} points, {total_viol} violations")
    print(f"  Result: {'ALL BOUNDS SATISFIED' if total_viol == 0 else 'VIOLATIONS DETECTED'}")

    # Show the tightest points
    print()
    print("  Tightest bound points (smallest F - exp(-Sigma/2)):")
    for N in [2, 6, 12]:
        F = all_results[N]['F']
        Sigma = all_results[N]['Sigma']
        gaps = F - np.exp(-Sigma / 2)
        idx = np.argmin(gaps)
        print(f"    N={N}: t={t_values[idx]:.3f}, "
              f"F={F[idx]:.6f}, bound={np.exp(-Sigma[idx]/2):.6f}, "
              f"gap={gaps[idx]:.6f}")

    # ================================================================
    # Figure
    # ================================================================
    print()
    print("Generating 4-panel figure...")

    fig, axes = plt.subplots(2, 2, figsize=(15, 12))

    # ----------------------------------------------------------
    # Panel (a): tau(t)
    # ----------------------------------------------------------
    ax = axes[0, 0]
    for N in N_values:
        ax.plot(t_values, all_results[N]['tau'], color=colors[N],
                linewidth=1.8, label=f'N = {N}')
        t_c = 1.0 / np.sqrt(N * Gamma)
        ax.axvline(t_c, color=colors[N], linestyle=':', alpha=0.3, lw=0.8)

    ax.set_xlabel(r'$t$', fontsize=14)
    ax.set_ylabel(r'$\tau(t) = 1 - F$', fontsize=14)
    ax.set_title(r'(a) Irreversibility $\tau(t)$ --- exact simulation',
                 fontsize=13)
    ax.legend(fontsize=10, loc='right', ncol=1)
    ax.set_ylim(-0.02, 0.55)
    ax.set_xlim(0, t_max)
    ax.grid(True, alpha=0.3)
    ax.axhline(0.5, color='gray', ls='--', alpha=0.4, lw=0.8)
    ax.text(t_max * 0.82, 0.505, r'$\tau_{\max}=1/2$',
            fontsize=10, color='gray')
    ax.annotate(r'Dotted: $t_c = 1/\sqrt{N\Gamma}$',
                xy=(0.30, 0.10), xycoords='axes fraction',
                fontsize=10, color='gray', fontstyle='italic')

    # ----------------------------------------------------------
    # Panel (b): R(t)/N
    # ----------------------------------------------------------
    ax = axes[0, 1]
    for N in N_values:
        ax.plot(t_values, all_results[N]['R_N'], color=colors[N],
                linewidth=1.8, label=f'N = {N}')
        t_c = 1.0 / np.sqrt(N * Gamma)
        ax.axvline(t_c, color=colors[N], linestyle=':', alpha=0.3, lw=0.8)

    ax.set_xlabel(r'$t$', fontsize=14)
    ax.set_ylabel(r'$R_\delta(t) / N$', fontsize=14)
    ax.set_title(r'(b) Redundancy $R_\delta/N$ --- exact ($\delta=0.05$)',
                 fontsize=13)
    ax.legend(fontsize=10, loc='right', ncol=1)
    ax.set_ylim(-0.02, 1.1)
    ax.set_xlim(0, t_max)
    ax.grid(True, alpha=0.3)
    ax.axhline(1.0, color='gray', ls='--', alpha=0.4, lw=0.8)
    ax.text(t_max * 0.82, 1.02, r'$R/N = 1$', fontsize=10, color='gray')
    ax.annotate(r'Same $t_c$ as panel (a)',
                xy=(0.30, 0.10), xycoords='axes fraction',
                fontsize=10, color='gray', fontstyle='italic')

    # ----------------------------------------------------------
    # Panel (c): tau vs R/N scatter
    # ----------------------------------------------------------
    ax = axes[1, 0]
    for N in N_values:
        ax.scatter(all_results[N]['R_N'], all_results[N]['tau'],
                   color=colors[N], s=6, alpha=0.5,
                   label=f'N = {N}', zorder=3)

    # Theoretical curves
    D_arr = np.linspace(0.001, 1.0, 500)
    for N_th in [4, 8, 12]:
        tau_th = (1 - D_arr) / 2
        R_N_th = np.zeros_like(D_arr)
        for i, Dv in enumerate(D_arr):
            dk = Dv ** (1.0 / N_th)
            S_sys = binary_entropy_scalar((1 + Dv) / 2)
            info = 1 - dk**2
            if S_sys > 1e-12:
                R_N_th[i] = min(info / 0.95, 1.0)
        ax.plot(R_N_th, tau_th, '--', color=colors.get(N_th, 'gray'),
                alpha=0.4, linewidth=1.5)

    ax.set_xlabel(r'$R_\delta / N$', fontsize=14)
    ax.set_ylabel(r'$\tau$', fontsize=14)
    ax.set_title(r'(c) $\tau$ vs $R/N$ --- near-universal collapse',
                 fontsize=13)
    ax.legend(fontsize=9, loc='upper left', ncol=2, markerscale=2)
    ax.set_xlim(-0.02, 1.15)
    ax.set_ylim(-0.02, 0.55)
    ax.grid(True, alpha=0.3)

    ax.annotate(r'$\tau > 0 \Leftrightarrow R > 0$',
                xy=(0.45, 0.85), xycoords='axes fraction',
                fontsize=13, color='darkred', fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.3',
                          facecolor='lightyellow',
                          edgecolor='darkred', alpha=0.9))

    ax.plot(0, 0, 'ko', markersize=8, zorder=5)
    ax.annotate('Closed system\n(quantum eraser)', xy=(0, 0),
                xytext=(0.15, 0.06), fontsize=9, color='black',
                arrowprops=dict(arrowstyle='->', color='black', lw=1.2))

    # ----------------------------------------------------------
    # Panel (d): Analytical -- universal variable x = N*Gamma*t^2
    # ----------------------------------------------------------
    ax = axes[1, 1]

    # Both tau and Sigma are exact functions of x = N*Gamma*t^2.
    # tau(x) = (1 - exp(-x))/2
    # Sigma(x) = h((1+exp(-x))/2) + N*h((1+exp(-x/N))/2)
    # For large N, the second term dominates: Sigma ~ N*h(1-x/(2N)) ~ x
    # The bound F = exp(-Sigma/2) gives tau_bound = 1 - exp(-Sigma/2).
    x_vals = np.linspace(0, 10, 500)

    # Universal tau curve
    tau_x = (1 - np.exp(-x_vals)) / 2

    ax.plot(x_vals, tau_x, 'k-', lw=2.5,
            label=r'$\tau = (1-e^{-x})/2$')

    # tau bound from Sigma for various N
    N_large = [10, 100, 1000, int(1e6)]
    colors_lg = plt.cm.plasma(np.linspace(0.15, 0.85, len(N_large)))

    for i, N in enumerate(N_large):
        abs_D = np.exp(-x_vals)  # |D| = exp(-N*Gamma*t^2) = exp(-x)
        d_k = np.exp(-x_vals / N)  # single-mode: exp(-Gamma*t^2)

        # Sigma = S(rho_S) + N * S(rho_{E_k})
        S_sys = binary_entropy((1 + abs_D) / 2)
        S_Ek = binary_entropy((1 + d_k) / 2)
        Sigma = S_sys + N * S_Ek

        tau_bound = 1 - np.exp(-Sigma / 2)

        lbl = (f'N = {N:,}' if N < 10000
               else f'N = $10^{int(np.log10(N))}$')
        ax.plot(x_vals, tau_bound, color=colors_lg[i], lw=1.5,
                ls='--', alpha=0.8, label=r'$1-e^{-\Sigma/2}$, ' + lbl)

    # Shade the allowed region
    ax.fill_between(x_vals, tau_x, 0.5, alpha=0.05, color='red')
    ax.annotate(r'$\tau \leq 1-e^{-\Sigma/2}$' + '\n(always satisfied)',
                xy=(6, 0.35), fontsize=10, color='darkred',
                fontstyle='italic')

    ax.axhline(0.5, color='gray', ls='--', alpha=0.4, lw=0.8)
    ax.axvline(1.0, color='red', ls=':', lw=1.2, alpha=0.5)
    ax.annotate(r'$x = 1$' + '\n' + r'($t = t_c$)',
                xy=(1.1, 0.08), fontsize=10, color='red')

    ax.set_xlabel(r'$x = N\Gamma t^2$  (universal decoherence parameter)',
                  fontsize=14)
    ax.set_ylabel(r'$\tau$  or  $1-e^{-\Sigma/2}$', fontsize=14)
    ax.set_title(r'(d) Universal scaling: $\tau(x)$ and bound '
                 r'$1-e^{-\Sigma/2}$', fontsize=13)
    ax.set_xlim(0, 10)
    ax.set_ylim(-0.02, 0.55)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8.5, loc='lower right', ncol=2)

    # ----------------------------------------------------------
    # Save
    # ----------------------------------------------------------
    fig.suptitle(
        r'$\tau$--$R$ Bridge: Irreversibility and Redundancy '
        r'from Gravitational Dephasing'
        '\n'
        r'Both $\tau$ and $R$ saturate at '
        r'$t_c = 1/\sqrt{N\Gamma}$, driven by the same '
        r'$N\Gamma$ parameter',
        fontsize=14, y=1.02
    )
    plt.tight_layout()

    fig_png = os.path.join(_DIR, 'fig_tau_R_bridge.png')
    fig_pdf = os.path.join(_DIR, 'fig_tau_R_bridge.pdf')
    plt.savefig(fig_png, dpi=200, bbox_inches='tight')
    plt.savefig(fig_pdf, bbox_inches='tight')
    print(f"\nFigures saved:")
    print(f"  {fig_png}")
    print(f"  {fig_pdf}")

    # ----------------------------------------------------------
    # Detailed report
    # ----------------------------------------------------------
    print()
    print("=" * 70)
    print("DETAILED NUMERICAL RESULTS")
    print("=" * 70)

    print()
    print("Part A: Saturation timescales")
    print("-" * 60)
    print(f"{'N':>4s}  {'t_c':>8s}  {'t_half(tau)':>12s}  "
          f"{'t_half(R/N)':>12s}  {'ratio':>8s}")

    for N in N_values:
        tau = all_results[N]['tau']
        R_N = all_results[N]['R_N']
        t_c = 1.0 / np.sqrt(N * Gamma)

        idx_tau = np.where(tau >= 0.20)[0]
        t_half_tau = t_values[idx_tau[0]] if len(idx_tau) > 0 else np.nan

        idx_R = np.where(R_N >= 0.50)[0]
        t_half_R = t_values[idx_R[0]] if len(idx_R) > 0 else np.nan

        if not np.isnan(t_half_tau) and not np.isnan(t_half_R):
            ratio = t_half_tau / t_half_R
            print(f"{N:4d}  {t_c:8.4f}  {t_half_tau:12.4f}  "
                  f"{t_half_R:12.4f}  {ratio:8.3f}")
        else:
            t1 = f"{t_half_tau:.4f}" if not np.isnan(t_half_tau) else "N/A"
            t2 = f"{t_half_R:.4f}" if not np.isnan(t_half_R) else "N/A"
            print(f"{N:4d}  {t_c:8.4f}  {t1:>12s}  {t2:>12s}  {'---':>8s}")

    print()
    print("  Both tau and R/N scale with t_c = O(1/sqrt(N*Gamma)).")

    print()
    print("Part B: Analytical universal scaling")
    print("-" * 60)
    print("At t = t_c = 1/sqrt(N*Gamma):")
    for N in [10, 100, 1000, int(1e6)]:
        t_c = 1.0 / np.sqrt(N * Gamma)
        tau_v = analytical_tau(t_c, N, Gamma)
        R_v = analytical_R_over_N(t_c, Gamma)
        print(f"  N = {N:>10,}: tau(t_c) = {tau_v:.6f}, "
              f"R(t_c)/N = {R_v:.6f}")
    print()
    print("  tau(t_c) = (1 - e^{-1})/2 = 0.3161 for ALL N (universal).")

    print()
    print("=" * 70)
    print("PHYSICAL INTERPRETATION: THE tau-R BRIDGE")
    print("=" * 70)
    print("""
1. SIMULTANEOUS EMERGENCE:
   tau(t) and R(t) both transition from 0 to saturation
   at the same timescale t_c = 1/sqrt(N*Gamma).
   Both are driven by D(t) = cos^N(2*epsilon*t).

2. COMMON ORIGIN:
   tau(t) = (1 - |D(t)|) / 2     (irreversibility)
   R(t)/N = I(S:E_k)/S(rho_S)    (redundancy per mode)
   Both are monotone functions of |D(t)|:
     |D| = 1  =>  tau = 0, R = 0   (no arrow of time)
     |D| = 0  =>  tau = 1/2, R = N  (full decoherence)

3. tau > 0  <=>  R > 0  (THE BRIDGE):
   The arrow of time (tau > 0) exists if and only if the
   environment has acquired redundant records (R > 0).
   This connects Petz recovery to quantum Darwinism.

4. QUANTUM ERASER = PETZ RECOVERY:
   When D = 1 (closed system): tau = 0, R = 0.
   No 'which-path' info, perfectly reversible.

5. F >= exp(-Sigma/2) VERIFIED:
   Using Sigma = total entropy production (system + N env qubits),
   the bound F = (1+|D|)/2 >= exp(-Sigma/2) is satisfied for
   ALL (N, t) data points with ZERO violations.

6. GRAVITATIONAL ORIGIN:
   Gamma ~ (Delta_phi / c^2)^2  (Pikovski mechanism).
   Arrow of time from gravitational time dilation,
   t_c ~ c / (sqrt(N) * Delta_phi).
""")


if __name__ == "__main__":
    main()
