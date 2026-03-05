"""
Quantum Eraser as Petz Recovery: Explicit Calculation.

Demonstrates that the quantum eraser is precisely the Petz recovery map
applied to the partial trace channel, achieving F = 1 and tau = 0.

Reference: Huang (2026), Section III and Supplemental Material S1.
"""

import numpy as np
from scipy.linalg import sqrtm
import matplotlib.pyplot as plt


def fidelity(rho, sigma):
    """Uhlmann fidelity F(rho, sigma)."""
    sqrt_rho = sqrtm(rho)
    M = sqrt_rho @ sigma @ sqrt_rho
    M = (M + M.conj().T) / 2
    eigvals = np.linalg.eigvalsh(M)
    eigvals = np.maximum(eigvals, 0)
    return np.real(np.sum(np.sqrt(eigvals)))


def partial_trace_E(rho_SE, dim_S, dim_E):
    """Partial trace over environment E: Tr_E[rho_SE]."""
    rho_SE = rho_SE.reshape(dim_S, dim_E, dim_S, dim_E)
    return np.trace(rho_SE, axis1=1, axis2=3)


def petz_recovery_partial_trace(rho_S, sigma_SE, dim_S, dim_E):
    """
    Petz recovery map for channel N = Tr_E.
    R(rho_S) = sigma_SE^{1/2} (N^dag(N(sigma)^{-1/2} rho_S N(sigma)^{-1/2})) sigma_SE^{1/2}

    For N = Tr_E: N^dag(X_S) = X_S tensor I_E
    """
    d = dim_S * dim_E

    # Step 1: N(sigma) = Tr_E[sigma_SE]
    sigma_S = partial_trace_E(sigma_SE, dim_S, dim_E)

    # Step 2: N(sigma)^{-1/2}
    sigma_S_inv_sqrt = np.linalg.pinv(sqrtm(sigma_S))

    # Step 3: sigma_SE^{1/2}
    sigma_SE_sqrt = sqrtm(sigma_SE)

    # Step 4: Construct the map
    # N^dag(sigma_S_inv_sqrt @ rho_S @ sigma_S_inv_sqrt) = (sigma_S_inv_sqrt @ rho_S @ sigma_S_inv_sqrt) tensor I_E
    inner = sigma_S_inv_sqrt @ rho_S @ sigma_S_inv_sqrt
    inner_lifted = np.kron(inner, np.eye(dim_E))

    # Step 5: Apply sigma_SE^{1/2} on both sides
    result = sigma_SE_sqrt @ inner_lifted @ sigma_SE_sqrt

    return (result + result.conj().T) / 2


def main():
    dim_S = 2  # Signal: {|A>, |B>}
    dim_E = 2  # Idler/Environment: {|d_A>, |d_B>}

    print("=" * 60)
    print("Quantum Eraser as Petz Recovery")
    print("=" * 60)

    # Sweep over partial distinguishability
    # |psi> = (|A>|d_A> + |B>|d_B(theta)>)/sqrt(2)
    # where <d_A|d_B(theta)> = cos(theta)
    # theta=0: no which-path info, theta=pi/2: full which-path info

    thetas = np.linspace(0, np.pi / 2, 50)
    fidelities = []
    taus = []
    visibilities = []

    for theta in thetas:
        # Detector states
        d_A = np.array([1, 0])
        d_B = np.array([np.cos(theta), np.sin(theta)])

        # Entangled state |psi_SE>
        A = np.array([1, 0])  # path A
        B = np.array([0, 1])  # path B

        psi_SE = (np.kron(A, d_A) + np.kron(B, d_B)) / np.sqrt(2)
        sigma_SE = np.outer(psi_SE, psi_SE.conj())

        # Signal state after tracing out idler
        rho_S = partial_trace_E(sigma_SE, dim_S, dim_E)

        # Apply Petz recovery
        rho_recovered = petz_recovery_partial_trace(rho_S, sigma_SE, dim_S, dim_E)

        # Fidelity with original
        F = fidelity(sigma_SE, rho_recovered)
        F = min(F, 1.0)
        tau = 1 - F

        # Interference visibility
        V = np.abs(rho_S[0, 1]) * 2  # off-diagonal coherence

        fidelities.append(F)
        taus.append(tau)
        visibilities.append(V)

    # Print key values
    print(f"\n{'theta/pi':>10} {'<dA|dB>':>10} {'F':>10} {'tau':>10} {'Visibility':>10}")
    print("-" * 55)
    for i, theta in enumerate(thetas):
        if theta in [0] or i == len(thetas) - 1 or abs(theta - np.pi / 4) < 0.04:
            overlap = np.cos(theta)
            print(f"{theta/np.pi:10.3f} {overlap:10.4f} {fidelities[i]:10.6f} {taus[i]:10.2e} {visibilities[i]:10.4f}")

    # Full which-path case
    print(f"\n{'theta=pi/2 (orthogonal detectors, full which-path info):':}")
    print(f"  F = {fidelities[-1]:.6f}")
    print(f"  tau = {taus[-1]:.2e}")
    print(f"  => Petz recovery is EXACT (F=1, tau=0)")
    print(f"     because signal-idler is a CLOSED system (Sigma=0)")

    # theta=0 case
    print(f"\n{'theta=0 (identical detectors, no which-path info):':}")
    print(f"  F = {fidelities[0]:.6f}")
    print(f"  tau = {taus[0]:.2e}")
    print(f"  Visibility = {visibilities[0]:.4f} (full interference)")

    # Plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    ax1.plot(thetas / np.pi, fidelities, 'b-', linewidth=2)
    ax1.axhline(y=1.0, color='r', linestyle='--', alpha=0.5, label='F = 1 (exact recovery)')
    ax1.set_xlabel(r'$\theta / \pi$ (detector distinguishability)', fontsize=13)
    ax1.set_ylabel(r'Recovery fidelity $F$', fontsize=13)
    ax1.set_title('Petz Recovery Fidelity for Quantum Eraser', fontsize=14)
    ax1.set_ylim(0.99, 1.005)
    ax1.legend(fontsize=12)
    ax1.grid(True, alpha=0.3)

    ax2.plot(thetas / np.pi, taus, 'b-', linewidth=2, label=r'$\tau = 1 - F$')
    ax2.set_xlabel(r'$\theta / \pi$ (detector distinguishability)', fontsize=13)
    ax2.set_ylabel(r'$\tau$ (temporal asymmetry)', fontsize=13)
    ax2.set_title(r'Quantum Eraser: $\tau = 0$ for all $\theta$', fontsize=14)
    ax2.legend(fontsize=12)
    ax2.grid(True, alpha=0.3)
    ax2.ticklabel_format(style='scientific', axis='y', scilimits=(0,0))

    plt.tight_layout()
    plt.savefig('simulations/fig_quantum_eraser_petz.png', dpi=150, bbox_inches='tight')
    plt.savefig('simulations/fig_quantum_eraser_petz.pdf', bbox_inches='tight')
    print("\nFigures saved: simulations/fig_quantum_eraser_petz.{png,pdf}")
    plt.show()


if __name__ == "__main__":
    main()
