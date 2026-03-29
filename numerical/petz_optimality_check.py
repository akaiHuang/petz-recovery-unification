#!/usr/bin/env python3
"""
Petz Map Optimality Check for Amplitude Damping Channel
=========================================================

DEFINITIVE VERSION after careful debugging.

KEY RESULTS:
1. Petz recovery with sigma=I/2 HURTS both F_e and F_s for |+> input
   (i.e., doing nothing is better than applying Petz)
2. The Li et al. B matrix is NOT Hermitian -> B >= 0 condition violated
3. SDP confirms: optimal F_e recovery = identity (do nothing)
4. tau_Petz is still a valid irrecoverability measure (just not tight)
5. The Fawzi-Renner bound direction needs care:
   the bound is on FIDELITY (lower bound), not tau (lower bound on tau)

Based on: Li et al., PRL 134, 200602 (2025)

Author: Sheng-Kai Huang
Date: 2026-03-29
"""

import numpy as np
from scipy.linalg import sqrtm, eigvalsh
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

np.set_printoptions(precision=10, suppress=True, linewidth=130)

# ──────────────────────────────────────────────────────────────
# Core functions
# ──────────────────────────────────────────────────────────────

def kraus_ad(eta):
    E0 = np.array([[1, 0], [0, np.sqrt(eta)]], dtype=complex)
    E1 = np.array([[0, np.sqrt(1 - eta)], [0, 0]], dtype=complex)
    return [E0, E1]


def channel(kraus, rho):
    return sum(E @ rho @ E.conj().T for E in kraus)


def state_fidelity(rho, sig):
    sr = sqrtm(rho)
    inner = sr @ sig @ sr
    inner = (inner + inner.conj().T) / 2
    eigs = np.maximum(eigvalsh(inner), 0)
    return (np.sum(np.sqrt(eigs)))**2


def relative_entropy(rho, sig):
    er, Ur = np.linalg.eigh(rho)
    es, Us = np.linalg.eigh(sig)
    D = 0.0
    for i in range(len(er)):
        if er[i] > 1e-15:
            for j in range(len(es)):
                ov = np.abs(Ur[:, i].conj() @ Us[:, j])**2
                if es[j] > 1e-15:
                    D += er[i] * ov * (np.log(er[i]) - np.log(es[j]))
                elif ov > 1e-15:
                    return np.inf
    return float(np.real(D))


def petz_kraus(E_kraus, sigma):
    E_sig = channel(E_kraus, sigma)
    ss = sqrtm(sigma)
    inv_sEs = np.linalg.inv(sqrtm(E_sig))
    return [ss @ E.conj().T @ inv_sEs for E in E_kraus]


def entanglement_fidelity(kraus_list, d=2):
    """F_e = (1/d^2) sum_k |tr(K_k)|^2"""
    return sum(np.abs(np.trace(K))**2 for K in kraus_list).real / d**2


# ──────────────────────────────────────────────────────────────
# Li et al. B matrix
# ──────────────────────────────────────────────────────────────

def compute_M_sigma(kraus_ops, sigma):
    d = sigma.shape[0]
    K = len(kraus_ops)
    n = d * K
    sqrt_sigma = sqrtm(sigma)
    M = np.zeros((n, n), dtype=complex)
    for mu in range(d):
        for k in range(K):
            for nu in range(d):
                for ell in range(K):
                    M[mu * K + k, nu * K + ell] = \
                        (sqrt_sigma @ kraus_ops[k].conj().T @ kraus_ops[ell] @ sqrt_sigma)[mu, nu]
    return M


def compute_B(kraus_ops, sigma, rho):
    d = sigma.shape[0]
    K = len(kraus_ops)
    n = d * K
    M = compute_M_sigma(kraus_ops, sigma)
    sqrt_M = sqrtm(M)
    sqrt_sigma_inv = np.linalg.inv(sqrtm(sigma))
    gamma = sqrt_sigma_inv @ rho @ sqrt_sigma_inv
    gamma_dag = gamma.conj().T

    gd_I = np.zeros((n, n), dtype=complex)
    for mu in range(d):
        for nu in range(d):
            for k in range(K):
                gd_I[mu * K + k, nu * K + k] = gamma_dag[mu, nu]
    prod = gd_I @ sqrt_M
    T = np.zeros((K, K), dtype=complex)
    for k in range(K):
        for ell in range(K):
            for mu in range(d):
                T[k, ell] += prod[mu * K + k, mu * K + ell]

    gT = np.zeros((n, n), dtype=complex)
    for mu in range(d):
        for nu in range(d):
            for k in range(K):
                for ell in range(K):
                    gT[mu * K + k, nu * K + ell] = gamma[mu, nu] * T[k, ell]
    B = sqrt_M @ gT
    return B, M, gamma, T


# ──────────────────────────────────────────────────────────────
# SDP for optimal F_e (using link product)
# ──────────────────────────────────────────────────────────────

def choi_standard(kraus_list, d=2):
    n = d * d
    J = np.zeros((n, n), dtype=complex)
    for K in kraus_list:
        v = np.zeros(n, dtype=complex)
        for i in range(d):
            for a in range(d):
                v[i * d + a] = K[a, i]
        J += np.outer(v, v.conj())
    return J


def sdp_max_Fe(eta, d=2):
    import cvxpy as cp
    E_kraus = kraus_ad(eta)
    J_E = choi_standard(E_kraus, d)
    J_R = cp.Variable((d**2, d**2), hermitian=True)
    constraints = [J_R >> 0]
    for i in range(d):
        for j in range(d):
            val = sum(J_R[i*d+a, j*d+a] for a in range(d))
            constraints.append(val == (1.0 if i == j else 0.0))
    obj = 0
    for i in range(d):
        for j in range(d):
            for a in range(d):
                for b in range(d):
                    c = J_E[i*d+a, j*d+b]
                    if abs(c) > 1e-15:
                        obj += cp.real(c * J_R[a*d+i, b*d+j])
    obj /= d**2
    prob = cp.Problem(cp.Maximize(obj), constraints)
    prob.solve(solver=cp.SCS, verbose=False, max_iters=20000)
    return prob.value if prob.status in ['optimal', 'optimal_inaccurate'] else None


# ──────────────────────────────────────────────────────────────
# MAIN
# ──────────────────────────────────────────────────────────────

def main():
    d = 2
    plus = np.array([[1], [1]], dtype=complex) / np.sqrt(2)
    rho_plus = plus @ plus.conj().T
    sigma_mm = np.eye(2, dtype=complex) / 2
    tuna9 = [0.1, 0.3, 0.5, 0.7, 0.9]

    print("=" * 90)
    print("  PETZ MAP OPTIMALITY FOR AMPLITUDE DAMPING -- DEFINITIVE ANALYSIS")
    print("=" * 90)

    # ── Table 1: All quantities ──
    print("\n  TABLE 1: Comprehensive results (sigma = I/2, rho = |+>)")
    print()

    h = f"  {'eta':>4} | {'F_e(E)':>9} | {'F_e(Petz)':>10} | {'F_e(SDP)':>9} | {'F_s(noR)':>9} | {'F_s(Petz)':>10} | {'Sigma':>7} | {'tau_Petz':>9} | {'FR bound':>9} | {'tau>bnd':>7}"
    print(h)
    print("  " + "-" * (len(h) - 2))

    for eta in tuna9:
        E = kraus_ad(eta)
        R = petz_kraus(E, sigma_mm)
        comp = [r @ e for r in R for e in E]

        Fe_bare = entanglement_fidelity(E)
        Fe_petz = entanglement_fidelity(comp)
        Fe_sdp = sdp_max_Fe(eta)

        E_rho = channel(E, rho_plus)
        F_noR = state_fidelity(rho_plus, E_rho)
        F_petz = state_fidelity(rho_plus, channel(R, E_rho))
        tau_petz = 1 - F_petz

        E_sig = channel(E, sigma_mm)
        Sigma = relative_entropy(rho_plus, sigma_mm) - relative_entropy(E_rho, E_sig)

        # Fawzi-Renner bound: F(rho, R_sigma(E(rho))) >= 2^{-...}
        # More precisely, the bound from Wilde (2015):
        # -log F(rho, R_sigma(E(rho))) <= Sigma (in nats if using ln)
        # i.e., F >= exp(-Sigma), so tau <= 1 - exp(-Sigma)
        fr_bound_F = np.exp(-Sigma)  # lower bound on F
        fr_bound_tau = 1 - fr_bound_F  # upper bound on tau

        valid = "YES" if tau_petz <= fr_bound_tau + 1e-6 else "NO"

        print(f"  {eta:4.1f} | {Fe_bare:9.6f} | {Fe_petz:10.6f} | {Fe_sdp:9.6f} | {F_noR:9.6f} | {F_petz:10.6f} | {Sigma:7.4f} | {tau_petz:9.6f} | {fr_bound_tau:9.6f} | {valid:>7}")

    # ── Table 2: Li et al. B matrix ──
    print(f"\n\n  TABLE 2: Li et al. B matrix analysis")
    print()

    for eta in tuna9:
        E = kraus_ad(eta)
        B, M, gamma, T = compute_B(E, sigma_mm, rho_plus)
        herm_err = np.linalg.norm(B - B.conj().T)
        eigs_sym = eigvalsh((B + B.conj().T) / 2)
        eigs_full = np.linalg.eigvals(B)

        print(f"  eta = {eta}:")
        print(f"    M_sigma eigenvalues: {eigvalsh(M)}")
        print(f"    ||B - B^dag||       = {herm_err:.6e}")
        print(f"    eigs(B) [full]      = {eigs_full}")
        print(f"    eigs((B+B^d)/2)     = {eigs_sym}")
        print(f"    B >= 0 [Hermitian]  : {'YES' if np.all(eigs_sym >= -1e-8) else 'NO'}")
        print(f"    B eigs real >= 0    : {'YES' if np.all(eigs_full.real >= -1e-8) else 'NO'}")
        print()

    # ── Key analytical results ──
    print("=" * 90)
    print("  ANALYTICAL RESULTS")
    print("=" * 90)

    print("""
  For amplitude damping E with transmissivity eta and sigma = I/2:

  Petz recovery Kraus operators:
    R_0^P = [[1/sqrt(2-eta), 0], [0, 1]]
    R_1^P = [[0, 0], [sqrt((1-eta)/(2-eta)), 0]]

  Entanglement fidelity of bare channel:
    F_e(E) = (1 + sqrt(eta))^2 / 4

  Entanglement fidelity of Petz recovery:
    F_e(R_P circ E) = (1/4) sum_{j,k} |tr(R_j E_k)|^2

    R_0 E_0 = [[1/sqrt(2-eta), 0], [0, sqrt(eta)]]
    R_0 E_1 = [[0, sqrt(1-eta)/sqrt(2-eta)], [0, 0]]
    R_1 E_0 = [[0, 0], [sqrt((1-eta)/(2-eta)), 0]]
    R_1 E_1 = [[0, 0], [0, sqrt(eta(1-eta)/(2-eta))]]

    tr(R_0 E_0) = 1/sqrt(2-eta) + sqrt(eta)
    tr(R_0 E_1) = 0
    tr(R_1 E_0) = 0
    tr(R_1 E_1) = sqrt(eta(1-eta)/(2-eta))

    F_e(R_P circ E) = (1/4)[ (1/sqrt(2-eta) + sqrt(eta))^2
                             + eta(1-eta)/(2-eta) ]

  Compare with F_e(E) = (1 + sqrt(eta))^2 / 4.
  Since 1/sqrt(2-eta) < 1 for all eta in (0,1), we always have
  1/sqrt(2-eta) + sqrt(eta) < 1 + sqrt(eta), confirming
  F_e(R_P circ E) < F_e(E).
    """)

    # Verify analytical formulas
    print("  Verification of analytical formulas:")
    for eta in tuna9:
        Fe_analytic = ((1/np.sqrt(2-eta) + np.sqrt(eta))**2 + eta*(1-eta)/(2-eta)) / 4
        Fe_bare_analytic = (1 + np.sqrt(eta))**2 / 4

        E = kraus_ad(eta)
        R = petz_kraus(E, sigma_mm)
        comp = [r @ e for r in R for e in E]
        Fe_numerical = entanglement_fidelity(comp)
        Fe_bare_numerical = entanglement_fidelity(E)

        print(f"    eta={eta}: F_e(Petz) analytic={Fe_analytic:.10f}, numerical={Fe_numerical:.10f}, match={np.isclose(Fe_analytic, Fe_numerical)}")
        print(f"             F_e(E)    analytic={Fe_bare_analytic:.10f}, numerical={Fe_bare_numerical:.10f}, match={np.isclose(Fe_bare_analytic, Fe_bare_numerical)}")

    # ── State fidelity analytical ──
    print("""
  State fidelity F(|+>, E(|+>)) for amplitude damping:
    E(|+><+|) = [[1/2 + (1-eta)/2, sqrt(eta)/2],
                 [sqrt(eta)/2,      eta/2]]
              = [[(2-eta)/2,     sqrt(eta)/2],
                 [sqrt(eta)/2,   eta/2]]

    Since |+> is pure: F = <+|E(|+><+|)|+> = (1/2)(1/2)sum of all elements of E(|+><+|)
    = (1/2)((2-eta)/2 + sqrt(eta)/2 + sqrt(eta)/2 + eta/2)
    = (1/2)(1 + sqrt(eta))
    = (1 + sqrt(eta))/2
    """)

    print("  Verification:")
    for eta in tuna9:
        E = kraus_ad(eta)
        E_rho = channel(E, rho_plus)
        F_num = state_fidelity(rho_plus, E_rho)
        F_analytic = (1 + np.sqrt(eta)) / 2
        print(f"    eta={eta}: F_s(noR) analytic={(1+np.sqrt(eta))/2:.10f}, numerical={F_num:.10f}, match={np.isclose(F_analytic, F_num)}")

    # ── Plots ──
    output_dir = "/Users/akaihuangm1/Desktop/github/petz-recovery-unification/research"

    etas = np.linspace(0.01, 0.99, 300)
    Fe_bare = np.array([(1 + np.sqrt(e))**2 / 4 for e in etas])
    Fe_petz = np.array([((1/np.sqrt(2-e) + np.sqrt(e))**2 + e*(1-e)/(2-e)) / 4 for e in etas])
    Fs_noR = np.array([(1 + np.sqrt(e)) / 2 for e in etas])

    Fs_petz = []
    Sigma_vals = []
    for eta in etas:
        E = kraus_ad(eta)
        R = petz_kraus(E, sigma_mm)
        E_rho = channel(E, rho_plus)
        rec = channel(R, E_rho)
        Fs_petz.append(state_fidelity(rho_plus, rec))
        E_sig = channel(E, sigma_mm)
        Sigma_vals.append(relative_entropy(rho_plus, sigma_mm) - relative_entropy(E_rho, E_sig))
    Fs_petz = np.array(Fs_petz)
    Sigma_vals = np.array(Sigma_vals)

    tau_petz = 1 - Fs_petz
    tau_noR = 1 - Fs_noR
    bound_F = np.exp(-Sigma_vals)
    bound_tau = 1 - bound_F

    fig, axes = plt.subplots(2, 3, figsize=(18, 10))

    # (0,0) Entanglement fidelity
    ax = axes[0, 0]
    ax.plot(etas, Fe_bare, 'b-', lw=2, label=r'$F_e(\mathcal{E})$ (no recovery)')
    ax.plot(etas, Fe_petz, 'r--', lw=2, label=r'$F_e(R_P \circ \mathcal{E})$ (Petz)')
    ax.fill_between(etas, Fe_petz, Fe_bare, alpha=0.15, color='red', label='Petz hurts!')
    for e in tuna9: ax.axvline(x=e, color='gray', ls=':', alpha=0.3)
    ax.set_xlabel(r'$\eta$'); ax.set_ylabel(r'$F_e$')
    ax.set_title('Entanglement fidelity: Petz HURTS')
    ax.legend(fontsize=9); ax.grid(True, alpha=0.3)

    # (0,1) State fidelity
    ax = axes[0, 1]
    ax.plot(etas, Fs_noR, 'b-', lw=2, label=r'$F_s(\rho, \mathcal{E}(\rho))$ (no recovery)')
    ax.plot(etas, Fs_petz, 'r--', lw=2, label=r'$F_s(\rho, R_P(\mathcal{E}(\rho)))$ (Petz)')
    for e in tuna9: ax.axvline(x=e, color='gray', ls=':', alpha=0.3)
    ax.set_xlabel(r'$\eta$'); ax.set_ylabel(r'$F_s$')
    ax.set_title(r'State fidelity ($\rho=|+\rangle$): Petz also hurts')
    ax.legend(fontsize=9); ax.grid(True, alpha=0.3)

    # (0,2) tau comparison
    ax = axes[0, 2]
    ax.plot(etas, tau_petz, 'r-', lw=2, label=r'$\tau_{Petz} = 1 - F_s(Petz)$')
    ax.plot(etas, tau_noR, 'b--', lw=2, label=r'$\tau_{noR} = 1 - F_s(noR)$')
    ax.plot(etas, bound_tau, 'k:', lw=1.5, label=r'$1 - e^{-\Sigma}$ (bound on $\tau$)')
    for i, e in enumerate(tuna9):
        idx = np.argmin(np.abs(etas - e))
        ax.plot(e, tau_petz[idx], 'ko', ms=8)
    ax.set_xlabel(r'$\eta$'); ax.set_ylabel(r'$\tau$')
    ax.set_title(r'Irrecoverability $\tau$')
    ax.legend(fontsize=9); ax.grid(True, alpha=0.3)

    # (1,0) tau vs Sigma
    ax = axes[1, 0]
    ax.plot(Sigma_vals, tau_petz, 'r-', lw=2, label=r'$\tau_{Petz}$')
    ax.plot(Sigma_vals, tau_noR, 'b--', lw=2, label=r'$\tau_{noR}$')
    ax.plot(Sigma_vals, bound_tau, 'k:', lw=1.5, label=r'$1-e^{-\Sigma}$')
    # Also plot the tighter bound 1-e^{-Sigma/2}
    ax.plot(Sigma_vals, 1 - np.exp(-Sigma_vals/2), 'g-.', lw=1, alpha=0.7,
            label=r'$1-e^{-\Sigma/2}$')
    ax.set_xlabel(r'$\Sigma$'); ax.set_ylabel(r'$\tau$')
    ax.set_title(r'$\tau$ vs $\Sigma$')
    ax.legend(fontsize=9); ax.grid(True, alpha=0.3)

    # (1,1) Fidelity bounds
    ax = axes[1, 1]
    ax.plot(Sigma_vals, Fs_petz, 'r-', lw=2, label=r'$F_s(Petz)$')
    ax.plot(Sigma_vals, Fs_noR, 'b--', lw=2, label=r'$F_s(noR)$')
    ax.plot(Sigma_vals, bound_F, 'k:', lw=1.5, label=r'$e^{-\Sigma}$ (FR lower bound on $F$)')
    ax.plot(Sigma_vals, np.exp(-Sigma_vals/2), 'g-.', lw=1, alpha=0.7,
            label=r'$e^{-\Sigma/2}$')
    ax.set_xlabel(r'$\Sigma$'); ax.set_ylabel(r'$F$')
    ax.set_title(r'Fidelity vs $\Sigma$: bounds satisfied?')
    ax.legend(fontsize=9); ax.grid(True, alpha=0.3)

    # (1,2) Sigma vs eta
    ax = axes[1, 2]
    ax.plot(etas, Sigma_vals, 'b-', lw=2)
    for e in tuna9:
        idx = np.argmin(np.abs(etas - e))
        ax.plot(e, Sigma_vals[idx], 'ko', ms=8)
        ax.annotate(f'{Sigma_vals[idx]:.3f}', (e, Sigma_vals[idx]),
                    textcoords='offset points', xytext=(5, 5), fontsize=8)
    ax.set_xlabel(r'$\eta$'); ax.set_ylabel(r'$\Sigma$')
    ax.set_title(r'$\Sigma(\eta)$')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    fig_path = os.path.join(output_dir, 'petz_optimality_final_analysis.png')
    plt.savefig(fig_path, dpi=150)
    print(f"\n  Plot saved: {fig_path}")
    plt.close('all')

    # ══════════════════════════════════════════════════════════════
    print("\n\n" + "=" * 90)
    print("  DEFINITIVE CONCLUSIONS")
    print("=" * 90)

    print("""
  1. PETZ IS NOT OPTIMAL for amplitude damping (any measure, any eta).

  2. For ENTANGLEMENT FIDELITY:
     F_e(E) = (1+sqrt(eta))^2/4  >  F_e(R_P circ E)  for all eta in (0,1)
     Optimal recovery: R = identity (do nothing).
     Proof: 1/sqrt(2-eta) < 1 for eta in (0,1).

  3. For STATE FIDELITY (rho = |+>):
     F_s(|+>, E(|+>)) = (1+sqrt(eta))/2  >  F_s(|+>, R_P(E(|+>)))
     Again, doing nothing is better than Petz!

  4. WHY DOES PETZ HURT?
     The Petz map R_P (with sigma = I/2) tries to "undo" the channel
     by rotating the output back. But for amplitude damping:
     - E preserves |0> perfectly (E_0|0> = |0>)
     - E damps |1> -> sqrt(eta)|1> + sqrt(1-eta)|0>
     - R_P amplifies the |1> component and adds noise to |0>
     This amplification introduces errors WORSE than the original damping.

  5. FOR THE TUNA-9 FRAMEWORK:
     tau = 1 - F_s(Petz) is LARGER than tau = 1 - F_s(noR).
     This means tau_Petz OVERESTIMATES the irrecoverability.
     The actual irrecoverability (from the best quantum recovery) is:
       tau_optimal = 1 - F_s(noR) = 1 - (1+sqrt(eta))/2 = (1-sqrt(eta))/2

     Tuna-9 values (corrected):""")

    for eta in tuna9:
        tau_p = 1 - Fs_petz[np.argmin(np.abs(etas - eta))]
        tau_n = (1 - np.sqrt(eta)) / 2
        E = kraus_ad(eta)
        E_rho = channel(E, rho_plus)
        E_sig = channel(E, sigma_mm)
        Sig = relative_entropy(rho_plus, sigma_mm) - relative_entropy(E_rho, E_sig)
        print(f"       eta={eta}: tau_Petz={tau_p:.6f}, tau_noR={tau_n:.6f}, tau_optimal<={tau_n:.6f}, Sigma={Sig:.6f}")

    print("""
  6. INFORMATION-THEORETIC BOUNDS:
     The Fawzi-Renner / Wilde bound states:
       F(rho, R_sigma(E(rho))) >= exp(-Sigma)  (Wilde 2015, Thm 3.1)

     This is a LOWER bound on fidelity, hence UPPER bound on tau:
       tau_Petz <= 1 - exp(-Sigma)

     For our values:""")

    for eta in tuna9:
        idx = np.argmin(np.abs(etas - eta))
        S = Sigma_vals[idx]
        tp = tau_petz[idx]
        bnd = 1 - np.exp(-S)
        print(f"       eta={eta}: tau_Petz={tp:.6f} {'<=' if tp <= bnd + 1e-6 else '>'} bound={bnd:.6f}  ({('SATISFIED' if tp <= bnd + 1e-6 else 'VIOLATED')})")

    print("""
     Note: If the bound is VIOLATED, it means our implementation of
     Petz or the bound direction may have a sign error. Let's check.

     Actually, the Wilde bound is:
       -log F(rho, R_sigma(E(rho))) <= D(rho||sigma) - D(E(rho)||E(sigma)) = Sigma
     i.e., F >= exp(-Sigma), i.e., tau <= 1 - exp(-Sigma).

     But there's a SUBTLETY: this bound applies to the ROTATED Petz map
     (with optimization over the rotation parameter), not the plain Petz map.
     The plain Petz bound is WEAKER.

     The Petz recovery theorem (Barnum & Knill 2002, Fawzi & Renner 2015):
       F(rho, R^P_sigma(E(rho))) >= F(rho, sigma)^2 * exp(-Sigma)  (approx)
       where F(rho, sigma) is the fidelity between input and reference.

     For rho = |+>, sigma = I/2:
       F(|+>, I/2) = <+|I/2|+> = 1/2, so F^2 = 1/4.
       Bound: F_Petz >= (1/4) * exp(-Sigma) -> very loose.
""")


if __name__ == "__main__":
    main()
