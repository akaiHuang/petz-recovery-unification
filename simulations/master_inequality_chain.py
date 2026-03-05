"""
Numerical verification of the Petz Recovery Bound (Theorem 3).

For the amplitude damping channel with parameter gamma:
    -log F^2  <=  Delta_D  (Data Processing Inequality)

where:
    F = fidelity of Petz recovery
    Delta_D = D(rho||sigma) - D(N(rho)||N(sigma))  (relative entropy drop)

Also computes entropy production Sigma for comparison.

Reference: Huang (2026), Theorem 3.
"""

import numpy as np
from scipy.linalg import sqrtm, logm
import matplotlib.pyplot as plt


def density_matrix(rho):
    """Ensure Hermiticity."""
    return (rho + rho.conj().T) / 2


def von_neumann_entropy(rho):
    """S(rho) = -Tr[rho log rho]."""
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
    """F(rho, sigma) = Tr[sqrt(sqrt(rho) sigma sqrt(rho))]."""
    sqrt_rho = sqrtm(density_matrix(rho))
    M = sqrt_rho @ density_matrix(sigma) @ sqrt_rho
    M = density_matrix(M)
    eigvals = np.linalg.eigvalsh(M)
    eigvals = np.maximum(eigvals, 0)
    return np.real(np.sum(np.sqrt(eigvals)))


def amplitude_damping_kraus(gamma):
    """Return Kraus operators for amplitude damping channel."""
    E0 = np.array([[1, 0], [0, np.sqrt(1 - gamma)]])
    E1 = np.array([[0, np.sqrt(gamma)], [0, 0]])
    return [E0, E1]


def apply_channel(rho, kraus_ops):
    """Apply channel N(rho) = sum_k E_k rho E_k^dag."""
    result = np.zeros_like(rho)
    for E in kraus_ops:
        result += E @ rho @ E.conj().T
    return result


def apply_adjoint(X, kraus_ops):
    """Apply adjoint N^dag(X) = sum_k E_k^dag X E_k."""
    result = np.zeros_like(X)
    for E in kraus_ops:
        result += E.conj().T @ X @ E
    return result


def petz_recovery(rho_out, sigma, kraus_ops):
    """Petz recovery map: R_sigma(Y) = sigma^{1/2} N^dag(N(sigma)^{-1/2} Y N(sigma)^{-1/2}) sigma^{1/2}."""
    sigma_out = apply_channel(sigma, kraus_ops)

    sqrt_sigma = sqrtm(density_matrix(sigma))
    sqrt_sigma_out_inv = np.linalg.pinv(sqrtm(density_matrix(sigma_out)))

    # Inner: N(sigma)^{-1/2} rho_out N(sigma)^{-1/2}
    inner = sqrt_sigma_out_inv @ rho_out @ sqrt_sigma_out_inv

    # Apply channel adjoint: N^dag(inner)
    adjoint_inner = apply_adjoint(inner, kraus_ops)

    # Outer: sigma^{1/2} ... sigma^{1/2}
    result = sqrt_sigma @ adjoint_inner @ sqrt_sigma
    return density_matrix(result)


def complementary_channel(rho, kraus_ops):
    """Complementary channel output (environment state)."""
    n = len(kraus_ops)
    rho_E = np.zeros((n, n), dtype=complex)
    for i, Ei in enumerate(kraus_ops):
        for j, Ej in enumerate(kraus_ops):
            rho_E[i, j] = np.trace(Ei @ rho @ Ej.conj().T)
    return density_matrix(rho_E)


def compute_quantities(gamma, rho, sigma):
    """Compute recovery fidelity, relative entropy drop, and entropy production."""
    kraus = amplitude_damping_kraus(gamma)
    rho_out = apply_channel(rho, kraus)
    sigma_out = apply_channel(sigma, kraus)

    # Delta D = D(rho||sigma) - D(N(rho)||N(sigma))  [information loss]
    D_before = relative_entropy(rho, sigma)
    D_after = relative_entropy(rho_out, sigma_out)
    delta_D = max(D_before - D_after, 0.0)

    # Entropy production Sigma = Delta S_sys + S_env
    rho_E = complementary_channel(rho, kraus)
    S_env = von_neumann_entropy(rho_E)
    Sigma = max(von_neumann_entropy(rho_out) - von_neumann_entropy(rho) + S_env, 0.0)

    # Petz recovery fidelity (with correct adjoint)
    rho_recovered = petz_recovery(rho_out, sigma, kraus)
    F = min(fidelity(rho, rho_recovered), 1.0)
    neg_log_F2 = -2 * np.log(max(F, 1e-15))

    return neg_log_F2, delta_D, Sigma, F


def main():
    states = {
        r'$|+\rangle$': np.array([[0.5, 0.5], [0.5, 0.5]]),
        r'$|1\rangle$': np.array([[0, 0], [0, 1.0]]),
        r'mixed $(p{=}0.3)$': np.array([[0.7, 0], [0, 0.3]]),
    }

    sigma = np.eye(2) / 2  # Maximally mixed reference
    gammas = np.linspace(0.01, 0.99, 50)

    print("=" * 70)
    print("Verification: -log F^2 <= Delta_D  (Petz Recovery Bound)")
    print("Reference state: sigma = I/2 (maximally mixed)")
    print("=" * 70)

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    all_pass = True
    for idx, (name, rho) in enumerate(states.items()):
        neg_log_F2_vals = []
        delta_D_vals = []
        Sigma_vals = []
        tau_vals = []
        tau_bound_vals = []

        print(f"\nInput state: {name}")
        print(f"{'gamma':>8} {'−logF²':>10} {'ΔD':>10} {'Σ':>10} {'τ':>8} {'bound':>8} {'OK?':>5}")
        print("-" * 60)

        state_ok = True
        for gamma in gammas:
            neg_log_F2, delta_D, Sigma, F = compute_quantities(gamma, rho, sigma)
            tau = 1 - F
            tau_bound = 1 - np.exp(-delta_D / 2)

            neg_log_F2_vals.append(neg_log_F2)
            delta_D_vals.append(delta_D)
            Sigma_vals.append(Sigma)
            tau_vals.append(tau)
            tau_bound_vals.append(tau_bound)

            tol = 1e-6
            ok = neg_log_F2 <= delta_D + tol and tau <= tau_bound + tol
            if not ok:
                state_ok = False
                all_pass = False

            if abs(gamma - 0.01) < 0.01 or abs(gamma - 0.25) < 0.01 or \
               abs(gamma - 0.5) < 0.01 or abs(gamma - 0.75) < 0.01 or \
               abs(gamma - 0.99) < 0.01:
                print(f"{gamma:8.2f} {neg_log_F2:10.4f} {delta_D:10.4f} "
                      f"{Sigma:10.4f} {tau:8.4f} {tau_bound:8.4f} "
                      f"{'YES' if ok else 'NO':>5}")

        print(f"  => -log F^2 <= Delta_D: {'PASS' if state_ok else 'FAIL'}")

        ax = axes[idx]
        ax.plot(gammas, neg_log_F2_vals, 'b-', linewidth=2, label=r'$-\log F^2$')
        ax.plot(gammas, delta_D_vals, 'r--', linewidth=2, label=r'$\Delta D$ (bound)')
        ax.plot(gammas, Sigma_vals, 'g:', linewidth=2, label=r'$\Sigma$ (entropy prod.)')
        ax.fill_between(gammas, neg_log_F2_vals, delta_D_vals, alpha=0.1, color='blue',
                        label='Gap (recovery slack)')
        ax.set_xlabel(r'Damping parameter $\gamma$', fontsize=13)
        ax.set_ylabel('Value (nats)', fontsize=13)
        ax.set_title(f'Input: {name}', fontsize=14)
        ax.legend(fontsize=10, loc='upper left')
        ax.grid(True, alpha=0.3)

    plt.suptitle(r'Petz Recovery Bound: $-\log F^2 \leq \Delta D$ (Amplitude Damping)',
                 fontsize=15, y=1.02)
    plt.tight_layout()

    print("\n" + "=" * 70)
    print(f"All bounds satisfied: {'YES' if all_pass else 'NO'}")
    print("=" * 70)
    print()
    print("PHYSICAL INSIGHT")
    print("-" * 40)
    print("Delta_D = information lost to environment (relative entropy drop)")
    print("-log F^2 = recovery difficulty (Petz map quality)")
    print("Sigma = total entropy production (thermodynamic cost)")
    print()
    print("The bound -log F^2 <= Delta_D is guaranteed by")
    print("the Data Processing Inequality (Petz 1988).")
    print("For the rotated Petz with Gibbs reference, Delta_D = Sigma.")

    plt.savefig('simulations/fig_master_inequality_chain.png', dpi=150, bbox_inches='tight')
    plt.savefig('simulations/fig_master_inequality_chain.pdf', bbox_inches='tight')
    print("\nFigures saved: simulations/fig_master_inequality_chain.{png,pdf}")
    plt.show()


if __name__ == "__main__":
    main()
