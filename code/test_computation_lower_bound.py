#!/usr/bin/env python3
"""
Numerical test of the Computation Lower Bound conjecture:

    tau >= alpha * C^2

where:
    tau   = 1 - F(rho, R_sigma,N(N(rho)))   (Petz recovery infidelity)
    C     = || [N(rho), N(sigma)] ||_1        (trace-norm of commutator of outputs)
    alpha > 0  is a conjectured universal constant

Author: Sheng-Kai Huang
Date: 2026-03-16
"""

import numpy as np
from scipy.linalg import sqrtm, svdvals, polar
import time
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

np.random.seed(42)

# ============================================================
# Utility functions
# ============================================================

def random_unitary(d):
    """Haar-random unitary via QR decomposition."""
    Z = (np.random.randn(d, d) + 1j * np.random.randn(d, d)) / np.sqrt(2)
    Q, R = np.linalg.qr(Z)
    diag = np.diag(R)
    phase = diag / np.abs(diag)
    return Q @ np.diag(phase)


def random_pure_state(d):
    """Random pure state density matrix."""
    psi = np.random.randn(d) + 1j * np.random.randn(d)
    psi /= np.linalg.norm(psi)
    return np.outer(psi, psi.conj())


def random_density_matrix(d, full_rank=True):
    """Random density matrix. If full_rank, ensure all eigenvalues > 0."""
    if full_rank:
        # Wishart-type: A A^dagger with A random
        A = np.random.randn(d, d) + 1j * np.random.randn(d, d)
        rho = A @ A.conj().T
        # Add small identity to guarantee full rank
        rho += 0.01 * np.eye(d)
        rho /= np.trace(rho)
    else:
        rho = random_pure_state(d)
    return rho


def random_kraus_operators(d, num_kraus=4):
    """
    Generate random CPTP map via Stinespring dilation.
    Returns list of Kraus operators K_i such that sum K_i^dag K_i = I.
    """
    # Random isometry V: d -> d * num_kraus
    d_out = d * num_kraus
    # Generate random d_out x d matrix, then make it an isometry
    A = np.random.randn(d_out, d) + 1j * np.random.randn(d_out, d)
    # QR to get isometry
    Q, _ = np.linalg.qr(A)
    V = Q[:, :d]  # d_out x d isometry: V^dag V = I_d

    # Extract Kraus operators
    kraus = []
    for i in range(num_kraus):
        K_i = V[i*d:(i+1)*d, :]  # d x d block
        kraus.append(K_i)

    # Verify CPTP: sum K_i^dag K_i should = I
    check = sum(K.conj().T @ K for K in kraus)
    # Normalize if needed (should be close to I already)
    # Use Cholesky-like correction
    S = sqrtm(check)
    S_inv = np.linalg.inv(S)
    kraus = [K @ S_inv for K in kraus]

    return kraus


def apply_channel(kraus, rho):
    """Apply quantum channel N(rho) = sum_i K_i rho K_i^dag."""
    d = rho.shape[0]
    result = np.zeros((d, d), dtype=complex)
    for K in kraus:
        result += K @ rho @ K.conj().T
    return result


def adjoint_channel(kraus, X):
    """Apply adjoint channel N^dag(X) = sum_i K_i^dag X K_i."""
    d = X.shape[0]
    result = np.zeros((d, d), dtype=complex)
    for K in kraus:
        result += K.conj().T @ X @ K
    return result


def matrix_sqrt_safe(A):
    """Compute matrix square root safely."""
    # Hermitianize
    A = (A + A.conj().T) / 2
    vals, vecs = np.linalg.eigh(A)
    vals = np.maximum(vals, 0)
    return vecs @ np.diag(np.sqrt(vals)) @ vecs.conj().T


def matrix_inv_sqrt_safe(A, eps=1e-12):
    """Compute matrix inverse square root safely."""
    A = (A + A.conj().T) / 2
    vals, vecs = np.linalg.eigh(A)
    vals = np.maximum(vals, eps)
    return vecs @ np.diag(1.0 / np.sqrt(vals)) @ vecs.conj().T


def petz_recovery(kraus, sigma, X):
    """
    Petz recovery map: R_{sigma,N}(X) = sigma^{1/2} N^dag( N(sigma)^{-1/2} X N(sigma)^{-1/2} ) sigma^{1/2}

    Input X is in the output space of N.
    Returns recovered state in input space.
    """
    N_sigma = apply_channel(kraus, sigma)
    N_sigma_inv_sqrt = matrix_inv_sqrt_safe(N_sigma)
    sigma_sqrt = matrix_sqrt_safe(sigma)

    # N(sigma)^{-1/2} X N(sigma)^{-1/2}
    sandwiched = N_sigma_inv_sqrt @ X @ N_sigma_inv_sqrt

    # N^dag of sandwiched
    adj_result = adjoint_channel(kraus, sandwiched)

    # sigma^{1/2} ... sigma^{1/2}
    recovered = sigma_sqrt @ adj_result @ sigma_sqrt

    # Hermitianize and normalize
    recovered = (recovered + recovered.conj().T) / 2
    tr = np.trace(recovered).real
    if tr > 1e-15:
        recovered /= tr
    return recovered


def fidelity(rho, sigma):
    """
    Quantum fidelity F(rho, sigma) = (Tr sqrt(sqrt(rho) sigma sqrt(rho)))^2
    For numerical stability, use SVD-based computation.
    """
    sqrt_rho = matrix_sqrt_safe(rho)
    M = sqrt_rho @ sigma @ sqrt_rho
    M = (M + M.conj().T) / 2
    vals = np.linalg.eigvalsh(M)
    vals = np.maximum(vals, 0)
    F = (np.sum(np.sqrt(vals)))**2
    return min(F.real, 1.0)


def trace_norm(A):
    """Trace norm ||A||_1 = sum of singular values."""
    return np.sum(svdvals(A))


def commutator_trace_norm(A, B):
    """||[A, B]||_1 = ||AB - BA||_1"""
    comm = A @ B - B @ A
    return trace_norm(comm)


def compute_tau_and_C(kraus, rho, sigma):
    """
    Compute tau and C for a given (N, rho, sigma) triple.

    tau = 1 - F(rho, R_{sigma,N}(N(rho)))
    C = ||[N(rho), N(sigma)]||_1
    """
    N_rho = apply_channel(kraus, rho)
    N_sigma = apply_channel(kraus, sigma)

    # C = trace norm of commutator
    C = commutator_trace_norm(N_rho, N_sigma)

    # Petz recovery of N(rho)
    recovered = petz_recovery(kraus, sigma, N_rho)

    # tau = 1 - F
    F = fidelity(rho, recovered)
    tau = 1.0 - F

    return tau, C


# ============================================================
# Specific channel families
# ============================================================

def amplitude_damping_kraus(gamma):
    """Amplitude damping channel with parameter gamma in [0,1]."""
    K0 = np.array([[1, 0], [0, np.sqrt(1 - gamma)]], dtype=complex)
    K1 = np.array([[0, np.sqrt(gamma)], [0, 0]], dtype=complex)
    return [K0, K1]


def depolarizing_kraus(p):
    """Depolarizing channel: rho -> (1-p)*rho + p*I/2, for p in [0,1]."""
    # Kraus operators: sqrt(1 - 3p/4)*I, sqrt(p/4)*X, sqrt(p/4)*Y, sqrt(p/4)*Z
    I = np.eye(2, dtype=complex)
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    Z = np.array([[1, 0], [0, -1]], dtype=complex)
    K0 = np.sqrt(1 - 3*p/4) * I
    K1 = np.sqrt(p/4) * X
    K2 = np.sqrt(p/4) * Y
    K3 = np.sqrt(p/4) * Z
    return [K0, K1, K2, K3]


def dephasing_kraus(lam):
    """Dephasing channel: rho -> (1-lam)*rho + lam*Z*rho*Z, for lam in [0,0.5]."""
    # Reparametrize: lam in [0,1] -> effective lambda = lam/2
    effective_lam = lam / 2
    I = np.eye(2, dtype=complex)
    Z = np.array([[1, 0], [0, -1]], dtype=complex)
    K0 = np.sqrt(1 - effective_lam) * I
    K1 = np.sqrt(effective_lam) * Z
    return [K0, K1]


# ============================================================
# Task 1: Single-step bound for d=2
# ============================================================

def task1_qubit(n_trials=50000):
    print("=" * 70)
    print("TASK 1: Single-step bound for d=2 (qubit)")
    print(f"  Trials: {n_trials}")
    print("=" * 70)

    d = 2
    ratios = []
    taus = []
    Cs = []
    skipped = 0

    t0 = time.time()
    for i in range(n_trials):
        rho = random_pure_state(d)
        sigma = random_density_matrix(d, full_rank=True)
        kraus = random_kraus_operators(d, num_kraus=4)

        try:
            tau, C = compute_tau_and_C(kraus, rho, sigma)
        except Exception:
            skipped += 1
            continue

        if C < 1e-10:
            skipped += 1
            continue
        if tau < 0:
            tau = 0.0  # numerical noise

        ratio = tau / (C**2)
        ratios.append(ratio)
        taus.append(tau)
        Cs.append(C)

        if (i + 1) % 10000 == 0:
            elapsed = time.time() - t0
            print(f"  ... {i+1}/{n_trials} done ({elapsed:.1f}s)")

    elapsed = time.time() - t0
    ratios = np.array(ratios)
    taus = np.array(taus)
    Cs = np.array(Cs)

    alpha_star = np.min(ratios)
    alpha_5th = np.percentile(ratios, 5)
    alpha_median = np.median(ratios)

    print(f"\n  Results (d=2):")
    print(f"    Valid trials:  {len(ratios)} (skipped {skipped})")
    print(f"    Time:          {elapsed:.1f}s")
    print(f"    alpha*(d=2) = min(tau/C^2) = {alpha_star:.8f}")
    print(f"    5th percentile:              {alpha_5th:.6f}")
    print(f"    Median:                      {alpha_median:.6f}")
    print(f"    Mean:                        {np.mean(ratios):.6f}")
    print(f"    tau range:     [{np.min(taus):.2e}, {np.max(taus):.6f}]")
    print(f"    C range:       [{np.min(Cs):.2e}, {np.max(Cs):.6f}]")

    violations = np.sum(ratios < 0)
    print(f"    Violations (tau/C^2 < 0): {violations}")

    # Check the 10 smallest ratios
    idx_sorted = np.argsort(ratios)[:10]
    print(f"\n    10 smallest tau/C^2 ratios:")
    for rank, idx in enumerate(idx_sorted):
        print(f"      #{rank+1}: tau/C^2 = {ratios[idx]:.8f}  (tau={taus[idx]:.2e}, C={Cs[idx]:.4f})")

    print()
    return alpha_star, ratios, taus, Cs


# ============================================================
# Task 2: Specific channel families
# ============================================================

def task2_channel_families(n_pairs=500):
    print("=" * 70)
    print("TASK 2: Specific channel families (d=2)")
    print(f"  Pairs per parameter value: {n_pairs}")
    print("=" * 70)

    d = 2
    param_values = np.linspace(0.01, 0.99, 100)

    channel_families = {
        "Amplitude Damping": amplitude_damping_kraus,
        "Depolarizing":      depolarizing_kraus,
        "Dephasing":         dephasing_kraus,
    }

    family_results = {}

    for name, kraus_func in channel_families.items():
        print(f"\n  --- {name} ---")
        t0 = time.time()
        all_ratios = []
        all_taus = []
        all_Cs = []
        skipped = 0

        for param in param_values:
            kraus = kraus_func(param)
            for _ in range(n_pairs):
                rho = random_pure_state(d)
                sigma = random_density_matrix(d, full_rank=True)
                try:
                    tau, C = compute_tau_and_C(kraus, rho, sigma)
                except Exception:
                    skipped += 1
                    continue

                if C < 1e-10:
                    skipped += 1
                    continue
                if tau < 0:
                    tau = 0.0

                all_ratios.append(tau / (C**2))
                all_taus.append(tau)
                all_Cs.append(C)

        elapsed = time.time() - t0
        all_ratios = np.array(all_ratios)

        alpha_star = np.min(all_ratios) if len(all_ratios) > 0 else float('inf')
        family_results[name] = alpha_star

        print(f"    Valid trials:  {len(all_ratios)} (skipped {skipped})")
        print(f"    Time:          {elapsed:.1f}s")
        print(f"    alpha*({name}) = min(tau/C^2) = {alpha_star:.8f}")
        print(f"    5th percentile:  {np.percentile(all_ratios, 5):.6f}" if len(all_ratios) > 0 else "")
        print(f"    Median:          {np.median(all_ratios):.6f}" if len(all_ratios) > 0 else "")

    print()
    return family_results


# ============================================================
# Task 3: Two-step composition test
# ============================================================

def task3_composition(n_trials=10000, alpha_global=None):
    print("=" * 70)
    print("TASK 3: Two-step composition test (d=2)")
    print(f"  Trials: {n_trials}")
    print(f"  Using alpha = {alpha_global:.8f}" if alpha_global else "  alpha not set")
    print("=" * 70)

    d = 2
    violations = 0
    valid = 0
    skipped = 0
    margin_ratios = []  # tau_total / (alpha/2 * (C1^2 + C2^2))

    t0 = time.time()
    for i in range(n_trials):
        rho = random_pure_state(d)
        sigma = random_density_matrix(d, full_rank=True)
        kraus1 = random_kraus_operators(d, num_kraus=4)
        kraus2 = random_kraus_operators(d, num_kraus=4)

        try:
            # Step 1 outputs
            N1_rho = apply_channel(kraus1, rho)
            N1_sigma = apply_channel(kraus1, sigma)
            C1 = commutator_trace_norm(N1_rho, N1_sigma)

            # Step 2 outputs (composed)
            N2N1_rho = apply_channel(kraus2, N1_rho)
            N2N1_sigma = apply_channel(kraus2, N1_sigma)
            C2 = commutator_trace_norm(N2N1_rho, N2N1_sigma)

            # Composed channel Kraus operators
            composed_kraus = []
            for K2 in kraus2:
                for K1 in kraus1:
                    composed_kraus.append(K2 @ K1)

            # Total tau for composed channel
            recovered = petz_recovery(composed_kraus, sigma, N2N1_rho)
            F = fidelity(rho, recovered)
            tau_total = max(1.0 - F, 0.0)

        except Exception:
            skipped += 1
            continue

        bound_value = C1**2 + C2**2
        if bound_value < 1e-20:
            skipped += 1
            continue

        valid += 1
        rhs = (alpha_global / 2) * bound_value
        if rhs > 0:
            margin = tau_total / rhs
            margin_ratios.append(margin)
            if tau_total < rhs - 1e-12:  # small tolerance
                violations += 1

        if (i + 1) % 5000 == 0:
            elapsed = time.time() - t0
            print(f"  ... {i+1}/{n_trials} done ({elapsed:.1f}s)")

    elapsed = time.time() - t0
    margin_ratios = np.array(margin_ratios)

    print(f"\n  Results:")
    print(f"    Valid trials:  {valid} (skipped {skipped})")
    print(f"    Time:          {elapsed:.1f}s")
    print(f"    Violations:    {violations} / {valid}")
    print(f"    Violation rate: {violations/max(valid,1)*100:.2f}%")
    if len(margin_ratios) > 0:
        print(f"    Min margin (tau/(alpha/2*(C1^2+C2^2))): {np.min(margin_ratios):.6f}")
        print(f"    5th percentile margin:                   {np.percentile(margin_ratios, 5):.6f}")
        print(f"    Median margin:                           {np.median(margin_ratios):.6f}")

    print()
    return violations, valid


# ============================================================
# Task 4: Dimension scaling
# ============================================================

def task4_dimension_scaling(dims_trials):
    print("=" * 70)
    print("TASK 4: Dimension scaling")
    print("=" * 70)

    results = {}

    for d, n_trials in dims_trials:
        print(f"\n  --- d={d}, {n_trials} trials ---")
        t0 = time.time()
        ratios = []
        skipped = 0

        for i in range(n_trials):
            rho = random_pure_state(d)
            sigma = random_density_matrix(d, full_rank=True)
            kraus = random_kraus_operators(d, num_kraus=4)

            try:
                tau, C = compute_tau_and_C(kraus, rho, sigma)
            except Exception:
                skipped += 1
                continue

            if C < 1e-10:
                skipped += 1
                continue
            if tau < 0:
                tau = 0.0

            ratios.append(tau / (C**2))

            if (i + 1) % 5000 == 0:
                elapsed = time.time() - t0
                print(f"    ... {i+1}/{n_trials} done ({elapsed:.1f}s)")

        elapsed = time.time() - t0
        ratios = np.array(ratios)

        alpha_star = np.min(ratios) if len(ratios) > 0 else float('inf')
        results[d] = alpha_star

        print(f"    Valid trials:  {len(ratios)} (skipped {skipped})")
        print(f"    Time:          {elapsed:.1f}s")
        print(f"    alpha*(d={d}) = min(tau/C^2) = {alpha_star:.8f}")
        print(f"    5th percentile:  {np.percentile(ratios, 5):.6f}" if len(ratios) > 0 else "")
        print(f"    Median:          {np.median(ratios):.6f}" if len(ratios) > 0 else "")

    # Scaling analysis
    print(f"\n  --- Scaling Analysis ---")
    dims = sorted(results.keys())
    for d in dims:
        print(f"    d={d}: alpha* = {results[d]:.8f}")

    if len(dims) >= 2:
        print(f"\n    Testing alpha*(d) ~ 1/d^2 hypothesis:")
        for d in dims:
            print(f"      d={d}: alpha*(d) * d^2 = {results[d] * d**2:.6f}")

        # Ratio test
        if len(dims) >= 2:
            d1, d2 = dims[0], dims[1]
            ratio = results[d1] / results[d2]
            expected = (d2/d1)**2
            print(f"\n      alpha*(d={d1})/alpha*(d={d2}) = {ratio:.4f}")
            print(f"      Expected if ~1/d^2: (d2/d1)^2 = {expected:.4f}")
            if len(dims) >= 3:
                d3 = dims[2]
                ratio2 = results[d1] / results[d3]
                expected2 = (d3/d1)**2
                print(f"      alpha*(d={d1})/alpha*(d={d3}) = {ratio2:.4f}")
                print(f"      Expected if ~1/d^2: (d3/d1)^2 = {expected2:.4f}")

    print()
    return results


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":
    print()
    print("*" * 70)
    print("*  COMPUTATION LOWER BOUND CONJECTURE: tau >= alpha * C^2          *")
    print("*  Numerical Test Suite                                            *")
    print("*" * 70)
    print()

    # Task 1
    alpha_d2, ratios_d2, taus_d2, Cs_d2 = task1_qubit(n_trials=50000)

    # Task 2
    family_results = task2_channel_families(n_pairs=500)

    # Task 3 — use alpha from Task 1
    # Use the minimum across all d=2 results
    alpha_all_d2 = min(alpha_d2, min(family_results.values()))
    print(f"  Global alpha for composition test: {alpha_all_d2:.8f}")
    violations, valid = task3_composition(n_trials=10000, alpha_global=alpha_all_d2)

    # Task 4
    dim_results = task4_dimension_scaling([(2, 10000), (3, 10000), (4, 5000)])

    # ============================================================
    # SUMMARY
    # ============================================================
    print()
    print("=" * 70)
    print("FINAL SUMMARY")
    print("=" * 70)
    print()
    print(f"  Task 1 (d=2, 50k trials):")
    print(f"    alpha*(d=2) = {alpha_d2:.8f}")
    print(f"    CONJECTURE {'SUPPORTED' if alpha_d2 > 0 else 'VIOLATED'}: alpha* > 0")
    print()
    print(f"  Task 2 (specific channels):")
    for name, alpha in family_results.items():
        print(f"    alpha*({name}) = {alpha:.8f}")
    print()
    print(f"  Task 3 (two-step composition, 10k trials):")
    print(f"    Violations: {violations}/{valid}")
    print(f"    Composition bound {'HOLDS' if violations == 0 else 'VIOLATED'}")
    print()
    print(f"  Task 4 (dimension scaling):")
    for d, alpha in sorted(dim_results.items()):
        print(f"    alpha*(d={d}) = {alpha:.8f}")
    if len(dim_results) >= 2:
        dims = sorted(dim_results.keys())
        for d in dims:
            print(f"    alpha*(d={d}) * d^2 = {dim_results[d] * d**2:.6f}")
    print()
    print("*" * 70)
    print("*  END OF NUMERICAL TEST                                           *")
    print("*" * 70)
