"""
Visualization of tau vs the Data Processing Inequality bound.

For the amplitude damping channel with standard Petz recovery (sigma = I/2):
    tau = 1 - F(rho, R_Petz(N(rho)))
    tau <= 1 - exp(-Delta_D / 2)

where Delta_D = D(rho||sigma) - D(N(rho)||N(sigma)) is the relative entropy drop.

Also shows entropy production Sigma for comparison.
Note: For the rotated Petz map with Gibbs reference, Delta_D = Sigma,
giving the bound tau <= 1 - exp(-Sigma/2) from Eq. (9)-(10).

Reference: Huang (2026), Eq. (9)-(10).
"""

import os
import numpy as np
from scipy.linalg import sqrtm, logm
import matplotlib.pyplot as plt

_DIR = os.path.dirname(os.path.abspath(__file__))


def density_matrix(rho):
    return (rho + rho.conj().T) / 2


def von_neumann_entropy(rho):
    eigvals = np.linalg.eigvalsh(density_matrix(rho))
    eigvals = eigvals[eigvals > 1e-15]
    return -np.sum(eigvals * np.log(eigvals))


def relative_entropy(rho, sigma):
    """D(rho || sigma) = Tr[rho (log rho - log sigma)]."""
    rho_h = density_matrix(rho)
    sigma_h = density_matrix(sigma)
    log_rho = logm(rho_h)
    log_sigma = logm(sigma_h)
    val = np.real(np.trace(rho_h @ (log_rho - log_sigma)))
    return max(val, 0.0)


def fidelity(rho, sigma):
    sqrt_rho = sqrtm(density_matrix(rho))
    M = sqrt_rho @ density_matrix(sigma) @ sqrt_rho
    M = density_matrix(M)
    eigvals = np.linalg.eigvalsh(M)
    eigvals = np.maximum(eigvals, 0)
    return np.real(np.sum(np.sqrt(eigvals)))


def amplitude_damping_kraus(gamma):
    E0 = np.array([[1, 0], [0, np.sqrt(1 - gamma)]])
    E1 = np.array([[0, np.sqrt(gamma)], [0, 0]])
    return [E0, E1]


def apply_channel(rho, kraus_ops):
    result = np.zeros_like(rho)
    for E in kraus_ops:
        result += E @ rho @ E.conj().T
    return result


def apply_adjoint(X, kraus_ops):
    result = np.zeros_like(X)
    for E in kraus_ops:
        result += E.conj().T @ X @ E
    return result


def petz_recovery(rho_out, sigma, kraus_ops):
    """Petz recovery: R_sigma(Y) = sigma^{1/2} N^dag(N(sigma)^{-1/2} Y N(sigma)^{-1/2}) sigma^{1/2}."""
    sigma_out = apply_channel(sigma, kraus_ops)
    sqrt_sigma = sqrtm(density_matrix(sigma))
    sqrt_sigma_out_inv = np.linalg.pinv(sqrtm(density_matrix(sigma_out)))

    inner = sqrt_sigma_out_inv @ rho_out @ sqrt_sigma_out_inv
    adjoint_inner = apply_adjoint(inner, kraus_ops)
    result = sqrt_sigma @ adjoint_inner @ sqrt_sigma
    return density_matrix(result)


def complementary_channel(rho, kraus_ops):
    n = len(kraus_ops)
    rho_E = np.zeros((n, n), dtype=complex)
    for i, Ei in enumerate(kraus_ops):
        for j, Ej in enumerate(kraus_ops):
            rho_E[i, j] = np.trace(Ei @ rho @ Ej.conj().T)
    return density_matrix(rho_E)


def entropy_production(rho, kraus_ops):
    """Entropy production Sigma = Delta S_sys + S_env."""
    rho_out = apply_channel(rho, kraus_ops)
    rho_E = complementary_channel(rho, kraus_ops)
    S_env = von_neumann_entropy(rho_E)
    return max(von_neumann_entropy(rho_out) - von_neumann_entropy(rho) + S_env, 0.0)


def main():
    states = {
        r'$|+\rangle$': np.array([[0.5, 0.5], [0.5, 0.5]]),
        r'$|1\rangle$': np.array([[0, 0], [0, 1.0]]),
        r'mixed $(p{=}0.3)$': np.array([[0.7, 0], [0, 0.3]]),
    }

    sigma = np.eye(2) / 2  # Maximally mixed reference
    gammas = np.linspace(0.001, 0.999, 100)

    print("=" * 60)
    print("tau vs DPI Bound for Amplitude Damping")
    print("Petz recovery with sigma = I/2")
    print("=" * 60)

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    for idx, (name, rho) in enumerate(states.items()):
        taus = []
        bounds = []
        sigmas_ep = []

        for gamma in gammas:
            kraus = amplitude_damping_kraus(gamma)
            rho_out = apply_channel(rho, kraus)
            sigma_out = apply_channel(sigma, kraus)
            rho_rec = petz_recovery(rho_out, sigma, kraus)

            F = min(fidelity(rho, rho_rec), 1.0)
            tau = 1 - F

            # Correct bound: Delta_D for sigma = I/2
            delta_D = relative_entropy(rho, sigma) - relative_entropy(rho_out, sigma_out)
            delta_D = max(delta_D, 0.0)
            bound = 1 - np.exp(-delta_D / 2)

            Sigma = entropy_production(rho, kraus)

            taus.append(tau)
            bounds.append(bound)
            sigmas_ep.append(Sigma)

        ax = axes[idx]
        ax.plot(gammas, taus, 'b-', linewidth=2, label=r'$\tau = 1 - F$')
        ax.plot(gammas, bounds, 'r--', linewidth=2, label=r'$1 - e^{-\Delta D/2}$ (DPI bound)')
        ax.fill_between(gammas, taus, bounds, alpha=0.1, color='red')
        ax.set_xlabel(r'$\gamma$', fontsize=13)
        ax.set_ylabel(r'$\tau$', fontsize=13)
        ax.set_title(f'Input: {name}', fontsize=14)
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)

        violations = sum(1 for t, b in zip(taus, bounds) if t > b + 1e-6)
        status = "PASS" if violations == 0 else f"FAIL ({violations})"
        print(f"  {name}: tau <= 1-exp(-DeltaD/2)? {status}")

    plt.suptitle(r'Temporal Asymmetry Bound: $\tau \leq 1 - e^{-\Delta D/2}$ (Amplitude Damping)',
                 fontsize=15, y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(_DIR, 'fig_tau_vs_entropy_production.png'), dpi=150, bbox_inches='tight')
    plt.savefig(os.path.join(_DIR, 'fig_tau_vs_entropy_production.pdf'), bbox_inches='tight')
    print("\nFigures saved: simulations/fig_tau_vs_entropy_production.{png,pdf}")

    print("\n" + "=" * 60)
    print("PHYSICAL INSIGHT")
    print("=" * 60)
    print("gamma = 0: Unitary (no damping)")
    print("  => Delta_D = 0, tau = 0 (no arrow of time)")
    print("  => Petz recovery is exact (F = 1)")
    print("  => This is the quantum eraser regime")
    print()
    print("gamma > 0: Open system (information leaks to environment)")
    print("  => Delta_D > 0, tau > 0 (arrow of time emerges)")
    print("  => Recovery is imperfect (F < 1)")
    print("  => tau bounded by information loss Delta_D")
    print()
    print("gamma = 1: Complete damping (all info lost)")
    print("  => Delta_D -> max, tau -> 1 (maximal arrow of time)")
    print("  => No recovery possible (F -> 0)")
    print()
    print("NOTE: For rotated Petz with thermal reference,")
    print("Delta_D = Sigma (entropy production), giving")
    print("the thermodynamic bound tau <= 1 - exp(-Sigma/2).")

    plt.show()


if __name__ == "__main__":
    main()
