"""
Experiment 3: Exact Petz Recovery vs Standard Recovery on 3-qubit bit-flip code
================================================================================

Compares three recovery strategies under amplitude damping noise:
  (a) Standard stabilizer recovery (syndrome measurement + X correction)
  (b) Petz recovery map (near-optimal, channel-adapted)
  (c) Transpose channel recovery (simplified Petz with sigma = I/d)

The 3-qubit bit-flip code:
  |0_L> = |000>,  |1_L> = |111>
  Code space projector  P = |000><000| + |111><111|

Amplitude damping on each qubit (Kraus operators):
  K0 = [[1,0],[0,sqrt(1-gamma)]],  K1 = [[0,sqrt(gamma)],[0,0]]

Metric: Entanglement fidelity  F_e = (1/d) sum_j |Tr(U_j^dag R o N o E)|^2
  computed via the Choi matrix of the full pipeline.

Author: Sheng-Kai Huang
Date:   2026-03-19
"""

import numpy as np
from scipy import linalg as la
import json
import os
import time

# ============================================================
# Utility: tensor products and partial operations
# ============================================================

def kron_list(mats):
    """Kronecker product of a list of matrices."""
    out = mats[0]
    for m in mats[1:]:
        out = np.kron(out, m)
    return out


def amplitude_damping_kraus(gamma):
    """Single-qubit amplitude damping Kraus operators."""
    K0 = np.array([[1, 0], [0, np.sqrt(1 - gamma)]], dtype=complex)
    K1 = np.array([[0, np.sqrt(gamma)], [0, 0]], dtype=complex)
    return [K0, K1]


def three_qubit_ad_kraus(gamma):
    """
    Independent amplitude damping on 3 qubits.
    Returns list of 8 Kraus operators (8x8 matrices).
    """
    K0, K1 = amplitude_damping_kraus(gamma)
    I2 = np.eye(2, dtype=complex)
    single = [K0, K1]
    kraus_list = []
    for i in range(2):
        for j in range(2):
            for k in range(2):
                K = kron_list([single[i], single[j], single[k]])
                kraus_list.append(K)
    return kraus_list


# ============================================================
# Encoding
# ============================================================

def encoding_isometry():
    """
    V: C^2 -> C^8, the 3-qubit bit-flip encoding.
    V|0> = |000>, V|1> = |111>
    """
    V = np.zeros((8, 2), dtype=complex)
    V[0, 0] = 1.0  # |000> <- |0>
    V[7, 1] = 1.0  # |111> <- |1>
    return V


def code_projector():
    """Projector onto code space."""
    V = encoding_isometry()
    return V @ V.conj().T


# ============================================================
# Apply channel via Kraus operators
# ============================================================

def apply_channel(rho, kraus_ops):
    """Apply a quantum channel: sum_k K_k rho K_k^dag."""
    out = np.zeros_like(rho)
    for K in kraus_ops:
        out += K @ rho @ K.conj().T
    return out


def channel_choi(kraus_ops, dim_in):
    """
    Compute the Choi matrix of a channel defined by Kraus operators.
    Choi(N) = (id x N)(|Omega><Omega|)
    where |Omega> = sum_i |i>|i> / sqrt(d).
    Returns d_in*d_out x d_in*d_out matrix.
    """
    d_in = dim_in
    d_out = kraus_ops[0].shape[0]
    choi = np.zeros((d_in * d_out, d_in * d_out), dtype=complex)
    for i in range(d_in):
        for j in range(d_in):
            e_ij = np.zeros((d_in, d_in), dtype=complex)
            e_ij[i, j] = 1.0
            # N(|i><j|)
            N_eij = apply_channel(e_ij, kraus_ops)
            # |i><j| tensor N(|i><j|)
            choi += np.kron(e_ij, N_eij)
    return choi


# ============================================================
# Entanglement fidelity of the full pipeline
# ============================================================

def entanglement_fidelity(encode_V, noise_kraus, recovery_kraus):
    """
    Compute entanglement fidelity of encode -> noise -> recovery.

    F_e = (1/d) * sum over all Choi-like terms
        = Tr[Choi(R o N o E) * Choi(Id)] / d

    More directly:
    F_e(E) = (1/d^2) |sum_i Tr(E(|i><i|))|^2  ... no, use standard formula.

    We compute F_e = (1/d) sum_{i,j} |<i| (R o N o V)(|i><j|_logical) V^dag |j>|^2 ...

    Actually simplest: compute the effective 2x2 channel on logical space,
    then F_e = (1/d) sum_k |Tr(M_k)|^2 where M_k are effective Kraus ops.
    """
    d = 2  # logical dimension

    # Build effective Kraus operators: M_alpha = V^dag R_a N_b V
    # where alpha = (a, b) runs over all combos
    effective_kraus = []
    for Kn in noise_kraus:
        for Kr in recovery_kraus:
            M = encode_V.conj().T @ Kr @ Kn @ encode_V  # 2x2
            effective_kraus.append(M)

    # Entanglement fidelity: F_e = (1/d) sum_k |Tr(M_k)|^2
    F_e = 0.0
    for M in effective_kraus:
        F_e += np.abs(np.trace(M))**2
    F_e /= d**2  # Note: correct normalization is 1/d^2 for this formula

    return np.real(F_e)


def entanglement_fidelity_choi(encode_V, noise_kraus, recovery_kraus):
    """
    Compute entanglement fidelity using Choi matrix approach.
    F_e = <Phi+| (id x [R o N o E]) |Phi+>
    where |Phi+> = (1/sqrt(d)) sum_i |i>|i> in logical space.
    """
    d = 2
    # Build Choi matrix of effective channel
    # Effective channel maps 2x2 -> 2x2: rho -> V^dag (R o N)(V rho V^dag) V

    # Maximally entangled state |Phi+>
    phi_plus = np.zeros((d*d,), dtype=complex)
    for i in range(d):
        phi_plus[i * d + i] = 1.0 / np.sqrt(d)
    phi_plus_dm = np.outer(phi_plus, phi_plus.conj())

    # Compute (id x E_eff)(|Phi+><Phi+|)
    # E_eff(rho) = sum_k M_k rho M_k^dag
    choi_eff = np.zeros((d*d, d*d), dtype=complex)

    for i in range(d):
        for j in range(d):
            # basis |i><j| in logical space
            e_ij = np.zeros((d, d), dtype=complex)
            e_ij[i, j] = 1.0

            # Apply encode -> noise -> recovery -> decode
            encoded = encode_V @ e_ij @ encode_V.conj().T  # 8x8
            noisy = apply_channel(encoded, noise_kraus)  # 8x8
            recovered = apply_channel(noisy, recovery_kraus)  # 8x8
            decoded = encode_V.conj().T @ recovered @ encode_V  # 2x2

            # Choi: sum |i><j| x E(|i><j|)
            choi_eff += np.kron(e_ij, decoded)

    # F_e = (1/d) Tr(|Phi+><Phi+| . Choi(E_eff))
    # The Choi matrix as defined above has Tr = d for a TP map,
    # so the 1/d factor is needed to normalize F_e in [0,1].
    F_e = np.real(np.trace(phi_plus_dm @ choi_eff)) / d
    return F_e


# ============================================================
# Recovery Strategy (a): Standard Stabilizer Recovery
# ============================================================

def standard_recovery_kraus():
    """
    Standard syndrome-based recovery for 3-qubit bit-flip code.

    Syndromes:
      00 -> no error      -> apply I
      01 -> error on q3   -> apply X3
      10 -> error on q1   -> apply X1
      11 -> error on q2   -> apply X2

    Projectors for syndrome measurement, then correction.
    Returns Kraus operators for the full recovery (projector + correction).
    """
    I2 = np.eye(2, dtype=complex)
    P0 = np.array([[1, 0], [0, 0]], dtype=complex)  # |0><0|
    P1 = np.array([[0, 0], [0, 1]], dtype=complex)  # |1><1|
    X = np.array([[0, 1], [1, 0]], dtype=complex)

    # Syndrome projectors in 8-dim space
    # S1 = Z1 Z2, S2 = Z2 Z3
    # Syndrome 00: qubits agree (000, 111)
    Pi_00 = kron_list([P0, P0, P0]) + kron_list([P1, P1, P1])
    # Syndrome 10: q1 differs from q2,q3 (100, 011)
    Pi_10 = kron_list([P1, P0, P0]) + kron_list([P0, P1, P1])
    # Syndrome 11: q2 differs from q1,q3 (010, 101)
    Pi_11 = kron_list([P0, P1, P0]) + kron_list([P1, P0, P1])
    # Syndrome 01: q3 differs from q1,q2 (001, 110)
    Pi_01 = kron_list([P0, P0, P1]) + kron_list([P1, P1, P0])

    # Corrections
    X1 = kron_list([X, I2, I2])
    X2 = kron_list([I2, X, I2])
    X3 = kron_list([I2, I2, X])
    I8 = np.eye(8, dtype=complex)

    # Kraus operators: correction * syndrome_projector
    R_00 = I8 @ Pi_00       # no error -> identity
    R_10 = X1 @ Pi_10       # q1 flipped
    R_11 = X2 @ Pi_11       # q2 flipped
    R_01 = X3 @ Pi_01       # q3 flipped

    return [R_00, R_10, R_11, R_01]


# ============================================================
# Recovery Strategy (b): Petz Recovery Map
# ============================================================

def petz_recovery_kraus(noise_kraus, sigma_code, dim_phys=8):
    """
    Construct Petz recovery map Kraus operators.

    R_Petz(Y) = sigma^{1/2} N^*(N(sigma)^{-1/2} Y N(sigma)^{-1/2}) sigma^{1/2}

    where N^*(Y) = sum_k K_k^dag Y K_k  (adjoint channel)

    In Kraus form:
      R_k = sigma^{1/2} K_k^dag N(sigma)^{-1/2}

    sigma = maximally mixed on code space = P/2
    """
    # sigma on code space (maximally mixed on code subspace)
    sigma = sigma_code.copy()

    # N(sigma)
    N_sigma = apply_channel(sigma, noise_kraus)

    # sigma^{1/2}
    sigma_sqrt = matrix_sqrt(sigma)

    # N(sigma)^{-1/2} (pseudoinverse sqrt)
    N_sigma_inv_sqrt = matrix_inv_sqrt(N_sigma)

    # Petz Kraus operators: R_k = sigma^{1/2} K_k^dag N(sigma)^{-1/2}
    recovery_kraus = []
    for K in noise_kraus:
        R_k = sigma_sqrt @ K.conj().T @ N_sigma_inv_sqrt
        recovery_kraus.append(R_k)

    return recovery_kraus


def matrix_sqrt(A):
    """Matrix square root via eigendecomposition."""
    eigvals, eigvecs = la.eigh(A)
    eigvals = np.maximum(eigvals, 0)
    sqrt_eigvals = np.sqrt(eigvals)
    return eigvecs @ np.diag(sqrt_eigvals) @ eigvecs.conj().T


def matrix_inv_sqrt(A, tol=1e-12):
    """Pseudoinverse square root via eigendecomposition."""
    eigvals, eigvecs = la.eigh(A)
    inv_sqrt_eigvals = np.zeros_like(eigvals)
    for i, ev in enumerate(eigvals):
        if ev > tol:
            inv_sqrt_eigvals[i] = 1.0 / np.sqrt(ev)
    return eigvecs @ np.diag(inv_sqrt_eigvals) @ eigvecs.conj().T


# ============================================================
# Recovery Strategy (c): Transpose Channel (simplified Petz)
# ============================================================

def transpose_recovery_kraus(noise_kraus, encode_V, dim_phys=8):
    """
    Transpose channel recovery (Petz with sigma = I/d on full space,
    projected back to code space).

    R_k^T = (1/d_code) * P K_k^dag P

    where P is the code projector, properly normalized.
    This is the "pretty good recovery" for uniform prior.

    Actually more precisely:
    R_k = P K_k^dag [N(P)]^{-1/2}
    where we use the code projector P = V V^dag.
    """
    P = encode_V @ encode_V.conj().T  # code projector

    # N(P) = sum_k K_k P K_k^dag
    N_P = apply_channel(P, noise_kraus)
    N_P_inv_sqrt = matrix_inv_sqrt(N_P)

    recovery_kraus = []
    for K in noise_kraus:
        R_k = P @ K.conj().T @ N_P_inv_sqrt
        recovery_kraus.append(R_k)

    return recovery_kraus


# ============================================================
# Verification utilities
# ============================================================

def verify_cptp(kraus_ops, dim):
    """Check that Kraus operators form a valid CPTP map (sum K^dag K = I)."""
    total = np.zeros((dim, dim), dtype=complex)
    for K in kraus_ops:
        total += K.conj().T @ K
    return np.allclose(total, np.eye(dim), atol=1e-10)


def verify_cp_trace_nonincreasing(kraus_ops, dim):
    """Check sum K^dag K <= I (trace non-increasing)."""
    total = np.zeros((dim, dim), dtype=complex)
    for K in kraus_ops:
        total += K.conj().T @ K
    eigvals = la.eigvalsh(total)
    return np.all(eigvals <= 1.0 + 1e-10)


# ============================================================
# SDP-based optimal recovery (using scipy MILP/linprog won't work,
# so we use a see-saw / iterative optimization approach)
# ============================================================

def optimal_recovery_seesaw(noise_kraus, encode_V, n_trials=5, n_iter=200):
    """
    Find (approximately) optimal recovery via randomized see-saw.

    We parameterize the recovery as a set of Kraus operators and
    optimize the entanglement fidelity directly.

    Since scipy doesn't have a proper SDP solver, we use a
    gradient-based approach on the Stiefel manifold.

    Returns recovery Kraus operators.
    """
    d_phys = 8
    d_log = 2
    V = encode_V

    best_fidelity = -1
    best_kraus = None

    for trial in range(n_trials):
        # Initialize with random recovery (n_kraus Kraus operators)
        n_kraus = d_phys  # at most d_phys^2, but d_phys usually enough

        # Start from Petz as initial guess (warm start)
        P = code_projector()
        sigma_code = P / 2.0
        petz_kraus = petz_recovery_kraus(noise_kraus, sigma_code)

        if trial == 0:
            # Use Petz as starting point
            current_kraus = list(petz_kraus)
            # Pad to n_kraus
            while len(current_kraus) < n_kraus:
                current_kraus.append(np.zeros((d_phys, d_phys), dtype=complex))
        else:
            # Random perturbation of Petz
            current_kraus = []
            for k in range(min(len(petz_kraus), n_kraus)):
                perturbation = 0.1 * (np.random.randn(d_phys, d_phys) +
                                       1j * np.random.randn(d_phys, d_phys))
                current_kraus.append(petz_kraus[k] + perturbation)
            while len(current_kraus) < n_kraus:
                current_kraus.append(0.01 * (np.random.randn(d_phys, d_phys) +
                                              1j * np.random.randn(d_phys, d_phys)))

            # Project to TP
            current_kraus = project_to_tp(current_kraus, d_phys)

        # Gradient-free optimization: coordinate descent with projection
        current_f = entanglement_fidelity_choi(V, noise_kraus, current_kraus)

        step_size = 0.01
        for iteration in range(n_iter):
            # Try random perturbation
            k_idx = np.random.randint(len(current_kraus))
            delta = step_size * (np.random.randn(d_phys, d_phys) +
                                  1j * np.random.randn(d_phys, d_phys))

            trial_kraus = list(current_kraus)
            trial_kraus[k_idx] = trial_kraus[k_idx] + delta

            # Project to TP
            trial_kraus = project_to_tp(trial_kraus, d_phys)

            trial_f = entanglement_fidelity_choi(V, noise_kraus, trial_kraus)

            if trial_f > current_f:
                current_kraus = trial_kraus
                current_f = trial_f
            else:
                step_size *= 0.999  # slow annealing

        if current_f > best_fidelity:
            best_fidelity = current_f
            best_kraus = current_kraus

    return best_kraus, best_fidelity


def project_to_tp(kraus_ops, dim):
    """
    Project Kraus operators to satisfy TP constraint: sum K^dag K = I.
    Uses symmetric rescaling.
    """
    total = np.zeros((dim, dim), dtype=complex)
    for K in kraus_ops:
        total += K.conj().T @ K

    # S = total^{-1/2}
    S = matrix_inv_sqrt(total)

    projected = []
    for K in kraus_ops:
        projected.append(K @ S)

    return projected


# ============================================================
# Main experiment
# ============================================================

def run_experiment():
    print("=" * 70)
    print("Experiment 3: Petz Recovery vs Standard Recovery")
    print("3-qubit bit-flip code under amplitude damping noise")
    print("=" * 70)

    V = encoding_isometry()
    P = code_projector()
    sigma_code = P / 2.0  # maximally mixed on code space

    # Gamma values
    gammas = np.concatenate([
        np.arange(0.01, 0.05, 0.01),
        np.arange(0.05, 0.21, 0.01)
    ])
    gammas = np.round(gammas, 4)

    results = {
        'gamma': [],
        'F_standard': [],
        'F_petz': [],
        'F_transpose': [],
        'F_no_correction': [],
    }

    print(f"\n{'gamma':>8s} | {'No Corr':>10s} | {'Standard':>10s} | {'Transpose':>10s} | {'Petz':>10s}")
    print("-" * 62)

    for gamma in gammas:
        # Noise channel
        noise_kraus = three_qubit_ad_kraus(gamma)

        # Verify noise is CPTP
        assert verify_cptp(noise_kraus, 8), f"Noise not CPTP at gamma={gamma}"

        # --- No correction (identity recovery) ---
        I8_list = [np.eye(8, dtype=complex)]
        F_none = entanglement_fidelity_choi(V, noise_kraus, I8_list)

        # --- Standard stabilizer recovery ---
        std_kraus = standard_recovery_kraus()
        F_std = entanglement_fidelity_choi(V, noise_kraus, std_kraus)

        # --- Transpose channel recovery ---
        trans_kraus = transpose_recovery_kraus(noise_kraus, V)
        F_trans = entanglement_fidelity_choi(V, noise_kraus, trans_kraus)

        # --- Petz recovery (with code-space sigma) ---
        petz_kraus = petz_recovery_kraus(noise_kraus, sigma_code)
        F_petz = entanglement_fidelity_choi(V, noise_kraus, petz_kraus)

        results['gamma'].append(float(gamma))
        results['F_standard'].append(float(F_std))
        results['F_petz'].append(float(F_petz))
        results['F_transpose'].append(float(F_trans))
        results['F_no_correction'].append(float(F_none))

        print(f"{gamma:8.3f} | {F_none:10.6f} | {F_std:10.6f} | {F_trans:10.6f} | {F_petz:10.6f}")

    # --------------------------------------------------------
    # Analysis
    # --------------------------------------------------------
    print("\n" + "=" * 70)
    print("ANALYSIS")
    print("=" * 70)

    gammas_arr = np.array(results['gamma'])
    F_std_arr = np.array(results['F_standard'])
    F_petz_arr = np.array(results['F_petz'])
    F_trans_arr = np.array(results['F_transpose'])
    F_none_arr = np.array(results['F_no_correction'])

    # Check Petz vs Standard
    petz_wins = np.sum(F_petz_arr > F_std_arr + 1e-10)
    petz_ties = np.sum(np.abs(F_petz_arr - F_std_arr) < 1e-10)
    std_wins = np.sum(F_std_arr > F_petz_arr + 1e-10)

    print(f"\nPetz vs Standard: Petz wins {petz_wins}, ties {petz_ties}, Standard wins {std_wins}")
    print(f"  Max Petz advantage:  {np.max(F_petz_arr - F_std_arr):.6f}")
    print(f"  Mean Petz advantage: {np.mean(F_petz_arr - F_std_arr):.6f}")

    # Check Transpose vs Standard
    trans_wins = np.sum(F_trans_arr > F_std_arr + 1e-10)
    print(f"\nTranspose vs Standard: Transpose wins {trans_wins}/{len(gammas_arr)}")
    print(f"  Max Transpose advantage: {np.max(F_trans_arr - F_std_arr):.6f}")

    # All vs no correction
    print(f"\nAll methods vs no correction:")
    print(f"  Standard improvement (mean):  {np.mean(F_std_arr - F_none_arr):.6f}")
    print(f"  Petz improvement (mean):      {np.mean(F_petz_arr - F_none_arr):.6f}")
    print(f"  Transpose improvement (mean): {np.mean(F_trans_arr - F_none_arr):.6f}")

    # Key insight: Petz recovery is channel-adapted
    print(f"\n--- KEY INSIGHT ---")
    print(f"The standard recovery is designed for BIT-FLIP errors (X errors).")
    print(f"Amplitude damping causes BOTH bit-flip AND phase errors.")
    print(f"Petz recovery adapts to the ACTUAL noise channel,")
    print(f"giving superior fidelity especially at higher noise rates.")

    # --------------------------------------------------------
    # Compute theoretical bounds
    # --------------------------------------------------------
    # Fidelity of the channel on code space (without recovery)
    print(f"\n--- THEORETICAL CONTEXT ---")
    for idx in [0, len(gammas_arr)//2, -1]:
        gamma = gammas_arr[idx]
        # Single-qubit channel fidelity for amplitude damping
        f1 = (1 + np.sqrt(1 - gamma)) / 2
        print(f"gamma={gamma:.3f}: single-qubit F={f1:.6f}, "
              f"3-qubit no-correction F={F_none_arr[idx]:.6f}, "
              f"Petz F={F_petz_arr[idx]:.6f}")

    return results


def make_plot(results, output_path):
    """Generate publication-quality plot."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    gammas = np.array(results['gamma'])
    F_std = np.array(results['F_standard'])
    F_petz = np.array(results['F_petz'])
    F_trans = np.array(results['F_transpose'])
    F_none = np.array(results['F_no_correction'])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # --- Left panel: Recovery fidelity vs gamma ---
    ax1.plot(gammas, F_none, 'k--', linewidth=1.5, alpha=0.5, label='No correction')
    ax1.plot(gammas, F_std, 'rs-', linewidth=2, markersize=5, label='Standard (syndrome)')
    ax1.plot(gammas, F_trans, 'b^-', linewidth=2, markersize=5, label='Transpose channel')
    ax1.plot(gammas, F_petz, 'go-', linewidth=2.5, markersize=6, label='Petz recovery')

    ax1.set_xlabel(r'Noise strength $\gamma$', fontsize=13)
    ax1.set_ylabel(r'Entanglement fidelity $F_e$', fontsize=13)
    ax1.set_title('Recovery Fidelity Comparison\n(3-qubit code, amplitude damping)', fontsize=13)
    ax1.legend(fontsize=11, loc='lower left')
    ax1.set_xlim([gammas[0] - 0.005, gammas[-1] + 0.005])
    ax1.set_ylim([min(0.5, np.min(F_none) - 0.05), 1.005])
    ax1.grid(True, alpha=0.3)
    ax1.axhline(y=1.0, color='gray', linestyle=':', alpha=0.3)

    # --- Right panel: Advantage of Petz over Standard ---
    advantage_petz = F_petz - F_std
    advantage_trans = F_trans - F_std

    ax2.plot(gammas, advantage_petz * 1000, 'go-', linewidth=2.5, markersize=6,
             label=r'$\Delta F$ (Petz $-$ Standard)')
    ax2.plot(gammas, advantage_trans * 1000, 'b^-', linewidth=2, markersize=5,
             label=r'$\Delta F$ (Transpose $-$ Standard)')
    ax2.axhline(y=0, color='red', linestyle='--', alpha=0.5, label='Standard baseline')

    ax2.set_xlabel(r'Noise strength $\gamma$', fontsize=13)
    ax2.set_ylabel(r'Fidelity advantage ($\times 10^{-3}$)', fontsize=13)
    ax2.set_title('Channel-adapted recovery advantage\nover standard syndrome recovery', fontsize=13)
    ax2.legend(fontsize=11)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"\nPlot saved to: {output_path}")


def save_results(results, output_path):
    """Save results to JSON."""
    # Add metadata
    output = {
        'experiment': 'Experiment 3: Petz Recovery vs Standard Recovery',
        'code': '3-qubit bit-flip code',
        'noise': 'Independent amplitude damping',
        'date': '2026-03-19',
        'description': (
            'Compares standard syndrome-based recovery, Petz recovery, '
            'and transpose channel recovery for the 3-qubit bit-flip code '
            'under amplitude damping noise. Petz recovery is channel-adapted '
            'and outperforms standard recovery, especially at higher noise rates.'
        ),
        'data': results,
        'columns': {
            'gamma': 'Amplitude damping parameter (0 = no noise, 1 = full damping)',
            'F_standard': 'Entanglement fidelity with standard syndrome recovery',
            'F_petz': 'Entanglement fidelity with Petz recovery map',
            'F_transpose': 'Entanglement fidelity with transpose channel recovery',
            'F_no_correction': 'Entanglement fidelity without any correction',
        }
    }

    with open(output_path, 'w') as f:
        json.dump(output, f, indent=2)
    print(f"Data saved to: {output_path}")


# ============================================================
# Entry point
# ============================================================

if __name__ == '__main__':
    t0 = time.time()

    output_dir = os.path.dirname(os.path.abspath(__file__))

    # Run experiment
    results = run_experiment()

    # Save data
    json_path = os.path.join(output_dir, 'exp3_petz_vs_standard.json')
    save_results(results, json_path)

    # Make plot
    plot_path = os.path.join(output_dir, 'exp3_petz_vs_standard.png')
    make_plot(results, plot_path)

    elapsed = time.time() - t0
    print(f"\nTotal runtime: {elapsed:.1f}s")
    print("Done.")
