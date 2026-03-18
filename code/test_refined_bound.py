#!/usr/bin/env python3
"""
=============================================================================
TEST: Refined Petz Recovery Lower Bound
=============================================================================

Conjecture:  tau >= C^2 * Sigma / (2 * D(rho||sigma))

where:
  tau   = 1 - F(rho, R_{sigma,N}(N(rho)))     (Petz recovery infidelity)
  C     = ||[N(rho), N(sigma)]||_1             (output commutator trace norm)
  Sigma = D(rho||sigma) - D(N(rho)||N(sigma))  (DPI deficit)
  D     = D(rho||sigma)                        (quantum relative entropy)

The original conjecture tau >= beta * C^2 * Sigma with constant beta > 0
was DISPROVEN (beta* = 0). This refined bound normalizes by 2*D(rho||sigma).

Also tests alternative bounds:
  (A) tau >= C^2 * Sigma   / (2*D)   [proposed]
  (B) tau >= C^2 * Sigma^2 / (2*D)   [quadratic Sigma]
  (C) tau >= C^2           / (2*D)   [without Sigma]

Author: Sheng-Kai Huang (2026)
Seed: 42
"""

import numpy as np
from scipy.linalg import svdvals, sqrtm
import time
import warnings
import sys

warnings.filterwarnings("ignore")
np.random.seed(42)

# ===========================================================================
# Core helper functions
# ===========================================================================

def random_pure_state(d):
    """Random pure state as density matrix."""
    psi = np.random.randn(d) + 1j * np.random.randn(d)
    psi /= np.linalg.norm(psi)
    return np.outer(psi, psi.conj())

def random_full_rank_state(d, rank_boost=0.05):
    """Random full-rank density matrix."""
    A = np.random.randn(d, d) + 1j * np.random.randn(d, d)
    rho = A @ A.conj().T + rank_boost * np.eye(d)
    return rho / np.trace(rho).real

def random_kraus_operators(d, num_kraus=4):
    """Generate random CPTP channel via Stinespring dilation."""
    d_env = num_kraus
    # Random isometry V: d -> d * d_env
    A = np.random.randn(d * d_env, d) + 1j * np.random.randn(d * d_env, d)
    Q, _ = np.linalg.qr(A)
    V = Q[:, :d]
    kraus = [V[i*d:(i+1)*d, :] for i in range(d_env)]
    # Enforce exact CPTP: renormalize
    check = sum(K.conj().T @ K for K in kraus)
    S_inv = np.linalg.inv(sqrtm(check))
    kraus = [K @ S_inv for K in kraus]
    return kraus

def random_unitary_channel(d):
    """Random unitary channel (single Kraus operator)."""
    A = np.random.randn(d, d) + 1j * np.random.randn(d, d)
    U, _ = np.linalg.qr(A)
    return [U]

def amplitude_damping_kraus(gamma):
    """Amplitude damping channel for qubit."""
    K0 = np.array([[1.0, 0.0], [0.0, np.sqrt(1.0 - gamma)]], dtype=complex)
    K1 = np.array([[0.0, np.sqrt(gamma)], [0.0, 0.0]], dtype=complex)
    return [K0, K1]

def depolarizing_kraus(p):
    """Depolarizing channel for qubit."""
    I2 = np.eye(2, dtype=complex)
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    Z = np.array([[1, 0], [0, -1]], dtype=complex)
    K0 = np.sqrt(1.0 - 3.0 * p / 4.0) * I2
    K1 = np.sqrt(p / 4.0) * X
    K2 = np.sqrt(p / 4.0) * Y
    K3 = np.sqrt(p / 4.0) * Z
    return [K0, K1, K2, K3]

def dephasing_kraus(p):
    """Dephasing channel for qubit."""
    I2 = np.eye(2, dtype=complex)
    Z = np.array([[1, 0], [0, -1]], dtype=complex)
    K0 = np.sqrt(1.0 - p) * I2
    K1 = np.sqrt(p) * Z
    return [K0, K1]

# ===========================================================================
# Quantum information functions
# ===========================================================================

def apply_channel(kraus, rho):
    """Apply CPTP channel: N(rho) = sum_i K_i rho K_i^dag."""
    return sum(K @ rho @ K.conj().T for K in kraus)

def adjoint_channel(kraus, X):
    """Adjoint (Heisenberg picture): N*(X) = sum_i K_i^dag X K_i."""
    return sum(K.conj().T @ X @ K for K in kraus)

def matrix_sqrt_safe(A):
    """Hermitian positive-semidefinite square root via eigendecomposition."""
    A = (A + A.conj().T) / 2
    vals, vecs = np.linalg.eigh(A)
    vals = np.maximum(vals, 0)
    return vecs @ np.diag(np.sqrt(vals)) @ vecs.conj().T

def matrix_inv_sqrt_safe(A, eps=1e-12):
    """Regularized inverse square root."""
    A = (A + A.conj().T) / 2
    vals, vecs = np.linalg.eigh(A)
    vals = np.maximum(vals, eps)
    return vecs @ np.diag(1.0 / np.sqrt(vals)) @ vecs.conj().T

def matrix_log_safe(A, eps=1e-15):
    """Matrix logarithm via eigendecomposition (Hermitian)."""
    A = (A + A.conj().T) / 2
    vals, vecs = np.linalg.eigh(A)
    vals = np.maximum(vals, eps)
    return vecs @ np.diag(np.log(vals)) @ vecs.conj().T

def petz_recovery(kraus, sigma, X):
    """
    Petz recovery map: R_{sigma,N}(X) = sigma^{1/2} N*(N(sigma)^{-1/2} X N(sigma)^{-1/2}) sigma^{1/2}
    normalized to be a density matrix.
    """
    N_sigma = apply_channel(kraus, sigma)
    N_sigma_inv_sqrt = matrix_inv_sqrt_safe(N_sigma)
    sigma_sqrt = matrix_sqrt_safe(sigma)
    sandwiched = N_sigma_inv_sqrt @ X @ N_sigma_inv_sqrt
    adj_result = adjoint_channel(kraus, sandwiched)
    recovered = sigma_sqrt @ adj_result @ sigma_sqrt
    recovered = (recovered + recovered.conj().T) / 2
    tr = np.trace(recovered).real
    if tr > 1e-15:
        recovered /= tr
    return recovered

def fidelity(rho, sigma):
    """Uhlmann fidelity F(rho, sigma) = (Tr sqrt(sqrt(rho) sigma sqrt(rho)))^2."""
    sqrt_rho = matrix_sqrt_safe(rho)
    M = sqrt_rho @ sigma @ sqrt_rho
    M = (M + M.conj().T) / 2
    vals = np.linalg.eigvalsh(M)
    vals = np.maximum(vals, 0)
    return min((np.sum(np.sqrt(vals)))**2, 1.0)

def trace_norm(A):
    """Schatten 1-norm."""
    return np.sum(svdvals(A))

def commutator_trace_norm(A, B):
    """||[A,B]||_1."""
    return trace_norm(A @ B - B @ A)

def relative_entropy(rho, sigma, eps=1e-15):
    """D(rho||sigma) = Tr(rho (ln rho - ln sigma))."""
    log_rho = matrix_log_safe(rho, eps)
    log_sigma = matrix_log_safe(sigma, eps)
    return np.trace(rho @ (log_rho - log_sigma)).real

# ===========================================================================
# Compute all quantities for a trial
# ===========================================================================

def compute_quantities(kraus, rho, sigma):
    """
    Returns (tau, C, Sigma, D_rho_sigma, valid) where valid=False if
    the trial should be skipped (e.g., near-zero denominators).
    """
    N_rho = apply_channel(kraus, rho)
    N_sigma = apply_channel(kraus, sigma)

    # D(rho||sigma) and D(N(rho)||N(sigma))
    D_pre = relative_entropy(rho, sigma)
    D_post = relative_entropy(N_rho, N_sigma)

    if D_pre < 1e-12:
        return None  # rho ~ sigma, trivial

    Sigma = D_pre - D_post
    if Sigma < -1e-10:
        # DPI violated => numerical issue
        return None

    Sigma = max(Sigma, 0.0)

    # Commutator
    C = commutator_trace_norm(N_rho, N_sigma)

    if C < 1e-14:
        return None  # bound is trivially 0; skip

    # Petz recovery
    recovered = petz_recovery(kraus, sigma, N_rho)
    F = fidelity(rho, recovered)
    tau = 1.0 - F

    # Safety: tau should be >= 0
    tau = max(tau, 0.0)

    return (tau, C, Sigma, D_pre)

# ===========================================================================
# Parametric qubit state from angles
# ===========================================================================

def pure_state_theta(theta):
    """Pure qubit state |psi> = cos(theta/2)|0> + sin(theta/2)|1>."""
    c = np.cos(theta / 2)
    s = np.sin(theta / 2)
    psi = np.array([c, s], dtype=complex)
    return np.outer(psi, psi.conj())

def sigma_from_p(p):
    """sigma = diag(p, 1-p) full-rank for 0<p<1."""
    return np.diag(np.array([p, 1.0 - p], dtype=complex))

# ===========================================================================
# Test runner for a batch
# ===========================================================================

def test_batch(name, generate_trial, n_trials, verbose=True):
    """
    Run n_trials. For each, generate_trial() -> (kraus, rho, sigma, d).
    Returns dict with statistics.
    """
    violations_A = 0   # tau >= C^2 Sigma / (2D)
    violations_B = 0   # tau >= C^2 Sigma^2 / (2D)
    violations_Bh = 0  # tau >= C^2 Sigma^2 / (4D)  [half-strength B]
    violations_C = 0   # tau >= C^2 / (2D)
    min_ratio_A = float('inf')
    min_ratio_B = float('inf')
    min_ratio_Bh = float('inf')
    min_ratio_C = float('inf')
    valid_count = 0

    t0 = time.time()

    for i in range(n_trials):
        trial = generate_trial()
        if trial is None:
            continue
        kraus, rho, sigma = trial

        result = compute_quantities(kraus, rho, sigma)
        if result is None:
            continue

        tau, C, Sigma, D_pre = result
        valid_count += 1

        C2 = C ** 2

        # Bound A: tau >= C^2 * Sigma / (2*D)
        bound_A = C2 * Sigma / (2.0 * D_pre)
        if bound_A > 1e-15:
            ratio_A = tau / bound_A
            min_ratio_A = min(min_ratio_A, ratio_A)
            if ratio_A < 1.0 - 1e-9:
                violations_A += 1

        # Bound B: tau >= C^2 * Sigma^2 / (2*D)
        bound_B = C2 * Sigma**2 / (2.0 * D_pre)
        if bound_B > 1e-15:
            ratio_B = tau / bound_B
            min_ratio_B = min(min_ratio_B, ratio_B)
            if ratio_B < 1.0 - 1e-9:
                violations_B += 1

        # Bound Bh: tau >= C^2 * Sigma^2 / (4*D)  [half-strength]
        bound_Bh = C2 * Sigma**2 / (4.0 * D_pre)
        if bound_Bh > 1e-15:
            ratio_Bh = tau / bound_Bh
            min_ratio_Bh = min(min_ratio_Bh, ratio_Bh)
            if ratio_Bh < 1.0 - 1e-9:
                violations_Bh += 1

        # Bound C: tau >= C^2 / (2*D)
        bound_C = C2 / (2.0 * D_pre)
        if bound_C > 1e-15:
            ratio_C = tau / bound_C
            min_ratio_C = min(min_ratio_C, ratio_C)
            if ratio_C < 1.0 - 1e-9:
                violations_C += 1

    elapsed = time.time() - t0

    if verbose:
        print(f"\n{'='*70}")
        print(f"  {name}")
        print(f"{'='*70}")
        print(f"  Trials: {n_trials}  |  Valid: {valid_count}  |  Time: {elapsed:.1f}s")
        print(f"  ---------------------------------------------------------------")
        print(f"  Bound A: tau >= C^2*Sigma/(2D)")
        print(f"    Violations: {violations_A}/{valid_count}")
        print(f"    Min ratio:  {min_ratio_A:.6f}" if min_ratio_A < float('inf') else "    Min ratio:  N/A")
        print(f"    HOLDS: {'YES' if violations_A == 0 else 'NO'}")
        print(f"  Bound B: tau >= C^2*Sigma^2/(2D)")
        print(f"    Violations: {violations_B}/{valid_count}")
        print(f"    Min ratio:  {min_ratio_B:.6f}" if min_ratio_B < float('inf') else "    Min ratio:  N/A")
        print(f"    HOLDS: {'YES' if violations_B == 0 else 'NO'}")
        print(f"  Bound Bh: tau >= C^2*Sigma^2/(4D)  [half-strength]")
        print(f"    Violations: {violations_Bh}/{valid_count}")
        print(f"    Min ratio:  {min_ratio_Bh:.6f}" if min_ratio_Bh < float('inf') else "    Min ratio:  N/A")
        print(f"    HOLDS: {'YES' if violations_Bh == 0 else 'NO'}")
        print(f"  Bound C: tau >= C^2/(2D)")
        print(f"    Violations: {violations_C}/{valid_count}")
        print(f"    Min ratio:  {min_ratio_C:.6f}" if min_ratio_C < float('inf') else "    Min ratio:  N/A")
        print(f"    HOLDS: {'YES' if violations_C == 0 else 'NO'}")

    return {
        'name': name,
        'n_trials': n_trials,
        'valid': valid_count,
        'viol_A': violations_A, 'min_A': min_ratio_A,
        'viol_B': violations_B, 'min_B': min_ratio_B,
        'viol_Bh': violations_Bh, 'min_Bh': min_ratio_Bh,
        'viol_C': violations_C, 'min_C': min_ratio_C,
    }

# ===========================================================================
# Trial generators
# ===========================================================================

def gen_random_d2():
    kraus = random_kraus_operators(2, num_kraus=np.random.randint(2, 5))
    rho = random_pure_state(2)
    sigma = random_full_rank_state(2)
    return (kraus, rho, sigma)

def gen_amplitude_damping():
    gamma = np.random.uniform(0.001, 0.999)
    kraus = amplitude_damping_kraus(gamma)
    theta = np.random.uniform(0.01, np.pi - 0.01)
    rho = pure_state_theta(theta)
    p = np.random.uniform(0.05, 0.95)
    sigma = sigma_from_p(p)
    return (kraus, rho, sigma)

def gen_depolarizing():
    p = np.random.uniform(0.001, 0.999)
    kraus = depolarizing_kraus(p)
    rho = random_pure_state(2)
    sigma = random_full_rank_state(2)
    return (kraus, rho, sigma)

def gen_dephasing():
    p = np.random.uniform(0.001, 0.999)
    kraus = dephasing_kraus(p)
    rho = random_pure_state(2)
    sigma = random_full_rank_state(2)
    return (kraus, rho, sigma)

def gen_unitary_d2():
    kraus = random_unitary_channel(2)
    rho = random_pure_state(2)
    sigma = random_full_rank_state(2)
    return (kraus, rho, sigma)

def gen_random_d(d):
    def _gen():
        kraus = random_kraus_operators(d, num_kraus=np.random.randint(2, d+2))
        rho = random_pure_state(d)
        sigma = random_full_rank_state(d)
        return (kraus, rho, sigma)
    return _gen

# ===========================================================================
# Adversarial test
# ===========================================================================

def adversarial_test():
    """
    Worst-case trajectory: gamma=eps, theta=pi/2-eps, p=1-eps
    and variant p=1-eps^2.
    """
    print(f"\n{'='*70}")
    print(f"  ADVERSARIAL TEST")
    print(f"{'='*70}")

    epsilons = [0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001, 0.00005, 0.00001, 0.000001]

    print(f"\n  Trajectory 1: gamma=eps, theta=pi/2-eps, p=1-eps")
    print(f"  {'eps':>12s} {'tau':>12s} {'C':>12s} {'Sigma':>12s} {'D':>12s} {'ratio_A':>12s} {'status':>8s}")
    print(f"  {'-'*80}")

    violations_1 = 0
    for eps in epsilons:
        gamma = eps
        theta = np.pi / 2 - eps
        p = 1.0 - eps

        kraus = amplitude_damping_kraus(gamma)
        rho = pure_state_theta(theta)
        sigma = sigma_from_p(p)

        result = compute_quantities(kraus, rho, sigma)
        if result is None:
            print(f"  {eps:>12.1e} {'SKIP':>12s}")
            continue

        tau, C, Sigma, D_pre = result
        bound_A = C**2 * Sigma / (2.0 * D_pre) if D_pre > 1e-15 else 0.0
        ratio_A = tau / bound_A if bound_A > 1e-15 else float('inf')
        status = "OK" if ratio_A >= 1.0 - 1e-9 else "FAIL"
        if status == "FAIL":
            violations_1 += 1
        print(f"  {eps:>12.1e} {tau:>12.6e} {C:>12.6e} {Sigma:>12.6e} {D_pre:>12.6e} {ratio_A:>12.6f} {status:>8s}")

    print(f"\n  Trajectory 2: gamma=eps, theta=pi/2-eps, p=1-eps^2")
    print(f"  {'eps':>12s} {'tau':>12s} {'C':>12s} {'Sigma':>12s} {'D':>12s} {'ratio_A':>12s} {'status':>8s}")
    print(f"  {'-'*80}")

    violations_2 = 0
    for eps in epsilons:
        gamma = eps
        theta = np.pi / 2 - eps
        p_val = 1.0 - eps**2

        if p_val >= 1.0:
            p_val = 1.0 - 1e-15

        kraus = amplitude_damping_kraus(gamma)
        rho = pure_state_theta(theta)
        sigma = sigma_from_p(p_val)

        result = compute_quantities(kraus, rho, sigma)
        if result is None:
            print(f"  {eps:>12.1e} {'SKIP':>12s}")
            continue

        tau, C, Sigma, D_pre = result
        bound_A = C**2 * Sigma / (2.0 * D_pre) if D_pre > 1e-15 else 0.0
        ratio_A = tau / bound_A if bound_A > 1e-15 else float('inf')
        status = "OK" if ratio_A >= 1.0 - 1e-9 else "FAIL"
        if status == "FAIL":
            violations_2 += 1
        print(f"  {eps:>12.1e} {tau:>12.6e} {C:>12.6e} {Sigma:>12.6e} {D_pre:>12.6e} {ratio_A:>12.6f} {status:>8s}")

    total_viol = violations_1 + violations_2
    print(f"\n  Adversarial violations: {total_viol} (traj1: {violations_1}, traj2: {violations_2})")
    return total_viol

# ===========================================================================
# MAIN
# ===========================================================================

def main():
    print("=" * 70)
    print("  REFINED PETZ RECOVERY LOWER BOUND TEST")
    print("  Conjecture: tau >= C^2 * Sigma / (2 * D(rho||sigma))")
    print("=" * 70)
    print(f"  Seed: 42  |  Date: 2026-03-17")
    print()

    all_results = []

    # ----- 1. Main test: 50k random, d=2 -----
    r = test_batch("1. MAIN TEST: Random CPTP, d=2", gen_random_d2, 50000)
    all_results.append(r)

    # ----- 2. Channel family tests (10k each) -----
    r = test_batch("2a. Amplitude Damping (d=2)", gen_amplitude_damping, 10000)
    all_results.append(r)

    r = test_batch("2b. Depolarizing (d=2)", gen_depolarizing, 10000)
    all_results.append(r)

    r = test_batch("2c. Dephasing (d=2)", gen_dephasing, 10000)
    all_results.append(r)

    r = test_batch("2d. Random Unitary (d=2)", gen_unitary_d2, 10000)
    all_results.append(r)

    # ----- 3. Adversarial -----
    adv_viol = adversarial_test()

    # ----- 4. Higher dimensions -----
    for d in [3, 4, 5]:
        r = test_batch(f"4. Higher dim d={d}", gen_random_d(d), 5000)
        all_results.append(r)

    # ===========================================================================
    # FINAL SUMMARY
    # ===========================================================================
    print(f"\n\n{'#'*70}")
    print(f"  FINAL SUMMARY")
    print(f"{'#'*70}")

    print(f"\n  {'Test':<35s} {'Bnd A':>8s} {'Bnd B':>8s} {'Bnd Bh':>8s} {'Bnd C':>8s}")
    print(f"  {'':35s} {'C2S/2D':>8s} {'C2S2/2D':>8s} {'C2S2/4D':>8s} {'C2/2D':>8s}")
    print(f"  {'-'*75}")

    total_viol_A = 0
    total_viol_B = 0
    total_viol_Bh = 0
    total_viol_C = 0
    global_min_A = float('inf')
    global_min_B = float('inf')
    global_min_Bh = float('inf')
    global_min_C = float('inf')

    for r in all_results:
        total_viol_A += r['viol_A']
        total_viol_B += r['viol_B']
        total_viol_Bh += r['viol_Bh']
        total_viol_C += r['viol_C']
        global_min_A = min(global_min_A, r['min_A'])
        global_min_B = min(global_min_B, r['min_B'])
        global_min_Bh = min(global_min_Bh, r['min_Bh'])
        global_min_C = min(global_min_C, r['min_C'])

        def _s(v):
            return f"{v}" if v > 0 else "0"

        print(f"  {r['name']:<35s} {_s(r['viol_A']):>8s} {_s(r['viol_B']):>8s} {_s(r['viol_Bh']):>8s} {_s(r['viol_C']):>8s}")

    print(f"  {'-'*75}")
    print(f"  {'TOTAL VIOLATIONS':<35s} {total_viol_A:>8d} {total_viol_B:>8d} {total_viol_Bh:>8d} {total_viol_C:>8d}")

    def _mr(v):
        return f"{v:.4f}" if v < float('inf') else "inf"
    print(f"  {'GLOBAL MIN RATIO':<35s} {_mr(global_min_A):>8s} {_mr(global_min_B):>8s} {_mr(global_min_Bh):>8s} {_mr(global_min_C):>8s}")

    print(f"\n  Adversarial test violations: {adv_viol}")

    print(f"\n  VERDICT:")
    print(f"    Bound A  [tau >= C^2*Sigma/(2D)]   : {'HOLDS' if total_viol_A == 0 and adv_viol == 0 else 'VIOLATED'}")
    print(f"    Bound B  [tau >= C^2*Sigma^2/(2D)]  : {'HOLDS' if total_viol_B == 0 else 'VIOLATED'}")
    print(f"    Bound Bh [tau >= C^2*Sigma^2/(4D)]  : {'HOLDS' if total_viol_Bh == 0 else 'VIOLATED'}")
    print(f"    Bound C  [tau >= C^2/(2D)]          : {'HOLDS' if total_viol_C == 0 else 'VIOLATED'}")

    # Tightest bound analysis
    print(f"\n  TIGHTNESS (lower min ratio = tighter bound):")
    bounds = [
        ("A:  C^2*Sigma/(2D)", global_min_A, total_viol_A + adv_viol),
        ("B:  C^2*Sigma^2/(2D)", global_min_B, total_viol_B),
        ("Bh: C^2*Sigma^2/(4D)", global_min_Bh, total_viol_Bh),
        ("C:  C^2/(2D)", global_min_C, total_viol_C),
    ]
    for name, minr, viol in sorted(bounds, key=lambda x: x[1]):
        status = "HOLDS" if viol == 0 else "VIOLATED"
        minr_str = f"{minr:.6f}" if minr < float('inf') else "inf"
        print(f"    {name:<25s}  min_ratio={minr_str:<12s}  {status}")

    holding = [(n, m) for n, m, v in bounds if v == 0]
    if holding:
        tightest = min(holding, key=lambda x: x[1])
        print(f"\n  >> TIGHTEST VALID BOUND: {tightest[0]} (min ratio = {tightest[1]:.6f})")
    else:
        print(f"\n  >> ALL BOUNDS VIOLATED in some regime.")

    print(f"\n{'#'*70}")
    print(f"  END OF REPORT")
    print(f"{'#'*70}")


if __name__ == "__main__":
    main()
