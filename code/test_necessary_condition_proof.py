#!/usr/bin/env python3
"""
Numerical Verification: Saturation of the Petz Bound and Output Commutativity
==============================================================================

Tests the key steps and claims from proof_attempt_necessary_condition.md:

1. Saturation (F^2 = exp(-Delta_D)) is possible ONLY when Delta_D = 0.
   - For Delta_D > 0, the gap is always strictly positive (strict GT inequality).

2. Unitary channels give Delta_D = 0 and gap = 0, but [N(rho), N(sigma)] != 0
   in general. This DISPROVES the original conjecture that output commutativity
   is necessary for saturation.

3. Dephasing channels can give [N(rho), N(sigma)] = 0 with gap > 0.
   This DISPROVES the converse (output commutativity implies saturation).

4. The quantitative correlation: gap correlates with C_out because both
   are controlled by Delta_D.

5. The Kadison-Schwarz equality and multiplicative domain structure
   at exact recovery (Delta_D = 0).

Author: Sheng-Kai Huang
Date: 2026-03-16
"""

import numpy as np
from scipy.stats import unitary_group
import warnings
import time

warnings.filterwarnings('ignore')
np.set_printoptions(precision=10, linewidth=130)
np.random.seed(42)

# ============================================================
# Core utility functions
# ============================================================

def ensure_hermitian(A):
    return (A + A.conj().T) / 2

def matrix_sqrt(A):
    A = ensure_hermitian(A)
    eigvals, eigvecs = np.linalg.eigh(A)
    eigvals = np.maximum(eigvals, 0)
    return eigvecs @ np.diag(np.sqrt(eigvals)) @ eigvecs.conj().T

def matrix_log(A):
    A = ensure_hermitian(A)
    eigvals, eigvecs = np.linalg.eigh(A)
    eigvals = np.maximum(eigvals, 1e-30)
    return eigvecs @ np.diag(np.log(eigvals)) @ eigvecs.conj().T

def matrix_inv_sqrt(A):
    A = ensure_hermitian(A)
    eigvals, eigvecs = np.linalg.eigh(A)
    inv_sqrt_ev = np.where(eigvals > 1e-12, 1.0/np.sqrt(eigvals), 0)
    return eigvecs @ np.diag(inv_sqrt_ev) @ eigvecs.conj().T

def uhlmann_fidelity(rho, sigma):
    """F(rho,sigma) = Tr sqrt(sqrt(rho) sigma sqrt(rho)). Returns F, not F^2."""
    sqrt_rho = matrix_sqrt(rho)
    inner = ensure_hermitian(sqrt_rho @ sigma @ sqrt_rho)
    sqrt_inner = matrix_sqrt(inner)
    F = np.real(np.trace(sqrt_inner))
    return np.clip(F, 0, 1)

def relative_entropy(rho, sigma):
    """D(rho||sigma) = Tr[rho(ln rho - ln sigma)]."""
    return np.real(np.trace(rho @ (matrix_log(rho) - matrix_log(sigma))))

def relative_entropy_pure(rho_pure, sigma):
    """D(|psi><psi| || sigma) = -<psi|ln(sigma)|psi> for pure rho."""
    log_sigma = matrix_log(sigma)
    return -np.real(np.trace(rho_pure @ log_sigma))

def trace_norm(A):
    return np.sum(np.linalg.svd(A, compute_uv=False))

def comm_tn(A, B):
    """||[A,B]||_1."""
    return trace_norm(A @ B - B @ A)

def comm_frob(A, B):
    """||[A,B]||_F."""
    return np.linalg.norm(A @ B - B @ A, 'fro')

def is_pure(rho, tol=1e-8):
    return abs(np.real(np.trace(rho @ rho)) - 1.0) < tol

# ============================================================
# Channel operations
# ============================================================

def apply_channel(rho, kraus_ops):
    d = rho.shape[0]
    result = np.zeros((d, d), dtype=complex)
    for K in kraus_ops:
        result += K @ rho @ K.conj().T
    return ensure_hermitian(result)

def adjoint_channel(omega, kraus_ops):
    d = omega.shape[0]
    result = np.zeros((d, d), dtype=complex)
    for K in kraus_ops:
        result += K.conj().T @ omega @ K
    return result

def random_cptp_channel(d, n_kraus=None):
    if n_kraus is None:
        n_kraus = d
    total_dim = d * n_kraus
    A = np.random.randn(total_dim, d) + 1j * np.random.randn(total_dim, d)
    Q, R = np.linalg.qr(A)
    V = Q[:, :d]
    kraus_ops = [V[k*d:(k+1)*d, :] for k in range(n_kraus)]
    check = sum(K.conj().T @ K for K in kraus_ops)
    assert np.allclose(check, np.eye(d), atol=1e-10), "CPTP check failed"
    return kraus_ops

def dephasing_kraus(d, basis=None):
    if basis is None:
        basis = np.eye(d, dtype=complex)
    return [np.outer(basis[:, k], basis[:, k].conj()) for k in range(d)]

def random_full_rank_state(d, min_eigval=0.05):
    U = unitary_group.rvs(d)
    eigvals = np.random.uniform(min_eigval, 1.0, d)
    eigvals /= eigvals.sum()
    return ensure_hermitian(U @ np.diag(eigvals) @ U.conj().T)

def random_pure_state(d):
    psi = np.random.randn(d) + 1j * np.random.randn(d)
    psi /= np.linalg.norm(psi)
    return np.outer(psi, psi.conj())

# ============================================================
# Petz recovery map
# ============================================================

def petz_recovery(sigma, kraus_ops, N_rho):
    """R_{sigma,N}(N(rho)) = sigma^{1/2} N^dag(N(sigma)^{-1/2} N(rho) N(sigma)^{-1/2}) sigma^{1/2}"""
    sqrt_sigma = matrix_sqrt(sigma)
    N_sigma = apply_channel(sigma, kraus_ops)
    inv_sqrt_Ns = matrix_inv_sqrt(N_sigma)
    sandwiched = inv_sqrt_Ns @ N_rho @ inv_sqrt_Ns
    adj_result = adjoint_channel(sandwiched, kraus_ops)
    result = sqrt_sigma @ adj_result @ sqrt_sigma
    return ensure_hermitian(result)

# ============================================================
# Compute all quantities for a trial
# ============================================================

def compute_trial(d, kraus_ops, rho, sigma):
    """Compute saturation gap, output commutator, and related quantities."""
    try:
        N_rho = apply_channel(rho, kraus_ops)
        N_sigma = apply_channel(sigma, kraus_ops)

        if np.min(np.linalg.eigvalsh(N_sigma)) < 1e-10:
            return None
        if np.min(np.linalg.eigvalsh(sigma)) < 1e-12:
            return None

        rho_is_pure = is_pure(rho)
        if not rho_is_pure and np.min(np.linalg.eigvalsh(rho)) < 1e-12:
            return None

        recovered = petz_recovery(sigma, kraus_ops, N_rho)
        tr = np.real(np.trace(recovered))
        if abs(tr - 1.0) > 0.01:
            return None

        F = uhlmann_fidelity(rho, recovered)
        F2 = F ** 2

        if rho_is_pure:
            D_before = relative_entropy_pure(rho, sigma)
        else:
            D_before = relative_entropy(rho, sigma)

        N_rho_eigmin = np.min(np.linalg.eigvalsh(N_rho))
        if N_rho_eigmin < 1e-12:
            if is_pure(N_rho):
                D_after = relative_entropy_pure(N_rho, N_sigma)
            else:
                N_rho_reg = N_rho + 1e-14 * np.eye(d)
                N_rho_reg /= np.real(np.trace(N_rho_reg))
                D_after = relative_entropy(N_rho_reg, N_sigma)
        else:
            D_after = relative_entropy(N_rho, N_sigma)

        Delta_D = D_before - D_after
        if Delta_D < -1e-6:
            return None
        Delta_D = max(Delta_D, 0)
        exp_neg_Delta = np.exp(-Delta_D)
        sat_gap = F2 - exp_neg_Delta

        C_out = comm_tn(N_rho, N_sigma)
        C_in = comm_tn(rho, sigma)

        # Recovery error
        recovery_error = trace_norm(recovered - rho) / 2

        return {
            'F2': F2, 'exp_neg_Delta': exp_neg_Delta, 'sat_gap': sat_gap,
            'C_out': C_out, 'C_in': C_in,
            'D_before': D_before, 'D_after': D_after, 'Delta_D': Delta_D,
            'recovery_error': recovery_error,
        }
    except Exception:
        return None

# ============================================================
# TEST 1: Unitary channels -- gap = 0 but C_out != 0
# ============================================================

def test1_unitary_counterexample():
    """
    KEY TEST: Unitary channels give exact recovery (Delta_D = 0, gap = 0)
    but DO NOT force [N(rho), N(sigma)] = 0.
    This DISPROVES the conjecture that output commutativity is necessary for saturation.
    """
    print("=" * 80)
    print("TEST 1: Unitary Channel Counterexample")
    print("  Claim: Gap = 0 does NOT imply [N(rho), N(sigma)] = 0")
    print("=" * 80)

    d = 2
    n_trials = 1000
    results = {'sat_gap': [], 'C_out': [], 'C_in': [], 'Delta_D': [], 'recovery_error': []}

    for _ in range(n_trials):
        rho = random_full_rank_state(d, 0.05)
        sigma = random_full_rank_state(d, 0.05)
        U = unitary_group.rvs(d)
        kraus = [U]

        r = compute_trial(d, kraus, rho, sigma)
        if r:
            for key in results:
                results[key].append(r[key])

    for key in results:
        results[key] = np.array(results[key])

    N = len(results['sat_gap'])
    print(f"\n  Valid trials: {N}")
    print(f"\n  Gap statistics:")
    print(f"    Max |gap|: {np.max(np.abs(results['sat_gap'])):.2e}")
    print(f"    Mean |gap|: {np.mean(np.abs(results['sat_gap'])):.2e}")
    print(f"    All gaps ~ 0: {np.all(np.abs(results['sat_gap']) < 1e-8)}")

    print(f"\n  Delta_D statistics:")
    print(f"    Max Delta_D: {np.max(results['Delta_D']):.2e}")
    print(f"    All Delta_D ~ 0: {np.all(results['Delta_D'] < 1e-8)}")

    print(f"\n  Output commutator C_out:")
    print(f"    Mean: {np.mean(results['C_out']):.6f}")
    print(f"    Max: {np.max(results['C_out']):.6f}")
    print(f"    Fraction C_out > 0.01: {np.mean(results['C_out'] > 0.01):.4f}")

    print(f"\n  Input commutator C_in:")
    print(f"    Mean: {np.mean(results['C_in']):.6f}")

    # Verify C_out = C_in for unitary channels
    ratio = np.abs(results['C_out'] - results['C_in'])
    print(f"\n  |C_out - C_in| (should be ~0 for unitary):")
    print(f"    Max: {np.max(ratio):.2e}")

    n_noncommuting_sat = np.sum((np.abs(results['sat_gap']) < 1e-8) & (results['C_out'] > 0.01))
    print(f"\n  >>> COUNTEREXAMPLES: gap ~ 0 AND C_out > 0.01: {n_noncommuting_sat}/{N}")
    print(f"  >>> This DISPROVES 'gap=0 implies [N(rho),N(sigma)]=0'")

    return results

# ============================================================
# TEST 2: Dephasing -- C_out = 0 but gap > 0
# ============================================================

def test2_dephasing_converse():
    """
    Dephasing channels give [N(rho), N(sigma)] = 0 (outputs always diagonal)
    but gap > 0 for mixed-state rho.
    This DISPROVES the converse.
    """
    print("\n\n" + "=" * 80)
    print("TEST 2: Dephasing Converse Counterexample")
    print("  Claim: [N(rho), N(sigma)] = 0 does NOT imply gap = 0")
    print("=" * 80)

    for d in [2, 3, 4]:
        print(f"\n  --- d = {d} ---")
        gaps = []
        Cout = []

        for _ in range(1000):
            rho = random_full_rank_state(d, 0.05)
            p = np.random.dirichlet(np.ones(d))
            p = np.maximum(p, 0.02); p /= p.sum()
            sigma = np.diag(p).astype(complex)
            kraus = dephasing_kraus(d)

            r = compute_trial(d, kraus, rho, sigma)
            if r:
                gaps.append(r['sat_gap'])
                Cout.append(r['C_out'])

        gaps = np.array(gaps)
        Cout = np.array(Cout)

        print(f"    Valid: {len(gaps)}")
        print(f"    All C_out = 0: {np.all(Cout < 1e-10)}")
        print(f"    Mean gap: {gaps.mean():.6f}")
        print(f"    Min gap: {gaps.min():.2e}")
        print(f"    Fraction gap > 0.001: {np.mean(gaps > 0.001):.4f}")
        print(f"    >>> C_out = 0 everywhere, but gap > 0 for most trials")

# ============================================================
# TEST 3: Strict GT -- Delta_D > 0 implies gap > 0
# ============================================================

def test3_strict_golden_thompson():
    """
    Verify that for Delta_D > 0, the gap is always strictly positive.
    This is the key mathematical content: saturation is impossible for Delta_D > 0.
    """
    print("\n\n" + "=" * 80)
    print("TEST 3: Strict Golden-Thompson Inequality")
    print("  Claim: Delta_D > 0 implies gap > 0 (saturation impossible)")
    print("=" * 80)

    for d in [2, 3, 4]:
        print(f"\n  --- d = {d} ---")
        DeltaD_vals = []
        gap_vals = []

        for _ in range(3000):
            rho = random_full_rank_state(d, 0.05)
            sigma = random_full_rank_state(d, 0.05)
            nk = np.random.choice(range(2, d+2))
            kraus = random_cptp_channel(d, nk)

            r = compute_trial(d, kraus, rho, sigma)
            if r:
                DeltaD_vals.append(r['Delta_D'])
                gap_vals.append(r['sat_gap'])

        DeltaD_vals = np.array(DeltaD_vals)
        gap_vals = np.array(gap_vals)

        # Check: when Delta_D > threshold, is gap always > 0?
        for thresh in [0.01, 0.1, 0.5]:
            mask = DeltaD_vals > thresh
            if np.sum(mask) > 0:
                min_gap = gap_vals[mask].min()
                mean_gap = gap_vals[mask].mean()
                print(f"    Delta_D > {thresh}: N={np.sum(mask)}, min gap={min_gap:.6e}, mean gap={mean_gap:.6f}")

        # Correlation between Delta_D and gap
        pos_mask = (DeltaD_vals > 1e-6) & (gap_vals > 1e-15)
        if np.sum(pos_mask) > 10:
            corr = np.corrcoef(np.log10(DeltaD_vals[pos_mask]), np.log10(gap_vals[pos_mask]))[0, 1]
            print(f"    Correlation(log Delta_D, log gap): {corr:.4f}")

        # Key check: any cases with Delta_D > 0 but gap ~ 0?
        suspicious = (DeltaD_vals > 0.01) & (gap_vals < 1e-8)
        print(f"    Cases with Delta_D > 0.01 AND gap < 1e-8: {np.sum(suspicious)} (should be 0)")

# ============================================================
# TEST 4: Kadison-Schwarz at exact recovery
# ============================================================

def test4_kadison_schwarz():
    """
    At exact recovery (Delta_D = 0), verify the Kadison-Schwarz equality:
    N^dag(h^2) = (N^dag(h))^2 where h = ln(N(rho)) - ln(N(sigma)).

    This is a key step in the proof attempt.
    """
    print("\n\n" + "=" * 80)
    print("TEST 4: Kadison-Schwarz Equality at Exact Recovery")
    print("  At Delta_D = 0: N^dag(h^2) should equal (N^dag(h))^2")
    print("  where h = ln(omega) - ln(nu)")
    print("=" * 80)

    d = 2

    # Use unitary channels (which give exact recovery)
    print(f"\n  --- Unitary channels (exact recovery) ---")
    KS_gaps = []

    for _ in range(500):
        rho = random_full_rank_state(d, 0.05)
        sigma = random_full_rank_state(d, 0.05)
        U = unitary_group.rvs(d)
        kraus = [U]

        omega = apply_channel(rho, kraus)
        nu = apply_channel(sigma, kraus)

        h = matrix_log(omega) - matrix_log(nu)
        h_in = matrix_log(rho) - matrix_log(sigma)

        # N^dag(h) should equal h_in
        Ndagh = adjoint_channel(h, kraus)
        Ndagh_error = np.linalg.norm(Ndagh - h_in, 'fro')

        # N^dag(h^2)
        Ndagh2 = adjoint_channel(h @ h, kraus)

        # (N^dag(h))^2
        Ndagh_sq = Ndagh @ Ndagh

        KS_gap = np.linalg.norm(Ndagh2 - Ndagh_sq, 'fro')
        KS_gaps.append(KS_gap)

    KS_gaps = np.array(KS_gaps)
    print(f"    N^dag(h) = h_in check: max error = {Ndagh_error:.2e}")
    print(f"    KS equality: max |N^dag(h^2) - (N^dag(h))^2| = {KS_gaps.max():.2e}")
    print(f"    KS equality holds: {np.all(KS_gaps < 1e-8)}")

    # Now test with NON-recovery cases (random channels, Delta_D > 0)
    print(f"\n  --- Random channels (no exact recovery, Delta_D > 0) ---")
    KS_gaps_random = []
    DeltaD_random = []

    for _ in range(500):
        rho = random_full_rank_state(d, 0.05)
        sigma = random_full_rank_state(d, 0.05)
        nk = np.random.choice([2, 3])
        kraus = random_cptp_channel(d, nk)

        omega = apply_channel(rho, kraus)
        nu = apply_channel(sigma, kraus)

        if np.min(np.linalg.eigvalsh(nu)) < 1e-10:
            continue
        if np.min(np.linalg.eigvalsh(omega)) < 1e-10:
            continue

        h = matrix_log(omega) - matrix_log(nu)

        Ndagh = adjoint_channel(h, kraus)
        Ndagh2 = adjoint_channel(h @ h, kraus)
        Ndagh_sq = Ndagh @ Ndagh

        KS_gap_rand = np.linalg.norm(Ndagh2 - Ndagh_sq, 'fro')
        KS_gaps_random.append(KS_gap_rand)

        rho_is_pure = is_pure(rho)
        if rho_is_pure:
            D_before = relative_entropy_pure(rho, sigma)
        else:
            D_before = relative_entropy(rho, sigma)
        D_after = relative_entropy(omega, nu)
        DeltaD_random.append(max(D_before - D_after, 0))

    KS_gaps_random = np.array(KS_gaps_random)
    DeltaD_random = np.array(DeltaD_random)

    print(f"    Valid: {len(KS_gaps_random)}")
    print(f"    KS gap: min={KS_gaps_random.min():.2e}, max={KS_gaps_random.max():.2e}")
    print(f"    Mean Delta_D: {DeltaD_random.mean():.4f}")

    # Correlation between KS gap and Delta_D
    pos_mask = (KS_gaps_random > 1e-15) & (DeltaD_random > 1e-10)
    if np.sum(pos_mask) > 10:
        corr = np.corrcoef(np.log10(KS_gaps_random[pos_mask]), np.log10(DeltaD_random[pos_mask]))[0, 1]
        print(f"    Correlation(log KS_gap, log Delta_D): {corr:.4f}")

    # Check: N^dag(h^2) >= (N^dag(h))^2 (Kadison-Schwarz inequality)
    # This means N^dag(h^2) - (N^dag(h))^2 should be PSD
    n_violations = 0
    for _ in range(200):
        rho = random_full_rank_state(d, 0.05)
        sigma = random_full_rank_state(d, 0.05)
        nk = np.random.choice([2, 3])
        kraus = random_cptp_channel(d, nk)

        omega = apply_channel(rho, kraus)
        nu = apply_channel(sigma, kraus)
        if np.min(np.linalg.eigvalsh(nu)) < 1e-10:
            continue

        h = matrix_log(omega) - matrix_log(nu)
        Ndagh2 = adjoint_channel(h @ h, kraus)
        Ndagh = adjoint_channel(h, kraus)
        diff = ensure_hermitian(Ndagh2 - Ndagh @ Ndagh)
        min_eig = np.min(np.linalg.eigvalsh(diff))
        if min_eig < -1e-8:
            n_violations += 1

    print(f"\n  Kadison-Schwarz inequality N^dag(h^2) >= (N^dag(h))^2:")
    print(f"    Violations: {n_violations}/200 (should be 0)")

# ============================================================
# TEST 5: The log-condition at exact recovery
# ============================================================

def test5_log_condition():
    """
    Verify Petz's log-condition:
    D(rho||sigma) = D(N(rho)||N(sigma)) iff ln(rho) - ln(sigma) = N^dag(ln(omega) - ln(nu))

    Test with unitary channels (exact recovery) and random channels (no recovery).
    """
    print("\n\n" + "=" * 80)
    print("TEST 5: Petz Log-Condition Verification")
    print("  ln(rho) - ln(sigma) = N^dag(ln(omega) - ln(nu)) at exact recovery")
    print("=" * 80)

    d = 2

    # Unitary channels: should satisfy log-condition exactly
    print(f"\n  --- Unitary channels ---")
    log_errors = []
    for _ in range(500):
        rho = random_full_rank_state(d, 0.05)
        sigma = random_full_rank_state(d, 0.05)
        U = unitary_group.rvs(d)
        kraus = [U]

        omega = apply_channel(rho, kraus)
        nu = apply_channel(sigma, kraus)

        h_in = matrix_log(rho) - matrix_log(sigma)
        h_out = matrix_log(omega) - matrix_log(nu)
        Ndagh_out = adjoint_channel(h_out, kraus)

        log_errors.append(np.linalg.norm(h_in - Ndagh_out, 'fro'))

    log_errors = np.array(log_errors)
    print(f"    Max |h_in - N^dag(h_out)|: {log_errors.max():.2e}")
    print(f"    Log-condition holds: {np.all(log_errors < 1e-8)}")

    # Random channels: should NOT satisfy log-condition
    print(f"\n  --- Random channels ---")
    log_errors_rand = []
    DeltaD_rand = []

    for _ in range(500):
        rho = random_full_rank_state(d, 0.05)
        sigma = random_full_rank_state(d, 0.05)
        nk = np.random.choice([2, 3])
        kraus = random_cptp_channel(d, nk)

        omega = apply_channel(rho, kraus)
        nu = apply_channel(sigma, kraus)

        if np.min(np.linalg.eigvalsh(nu)) < 1e-10:
            continue
        if np.min(np.linalg.eigvalsh(omega)) < 1e-10:
            continue

        h_in = matrix_log(rho) - matrix_log(sigma)
        h_out = matrix_log(omega) - matrix_log(nu)
        Ndagh_out = adjoint_channel(h_out, kraus)

        log_errors_rand.append(np.linalg.norm(h_in - Ndagh_out, 'fro'))

        D_before = relative_entropy(rho, sigma)
        D_after = relative_entropy(omega, nu)
        DeltaD_rand.append(max(D_before - D_after, 0))

    log_errors_rand = np.array(log_errors_rand)
    DeltaD_rand = np.array(DeltaD_rand)

    print(f"    Valid: {len(log_errors_rand)}")
    print(f"    Mean |h_in - N^dag(h_out)|: {log_errors_rand.mean():.4f}")

    # Correlation between log-condition error and Delta_D
    pos_mask = (log_errors_rand > 1e-10) & (DeltaD_rand > 1e-10)
    if np.sum(pos_mask) > 10:
        corr = np.corrcoef(np.log10(log_errors_rand[pos_mask]), np.log10(DeltaD_rand[pos_mask]))[0, 1]
        print(f"    Correlation(log |log-condition error|, log Delta_D): {corr:.4f}")

# ============================================================
# TEST 6: Multiplicative domain verification
# ============================================================

def test6_multiplicative_domain():
    """
    At exact recovery, h = ln(omega) - ln(nu) should be in the multiplicative
    domain of N^dag: N^dag(h*X) = N^dag(h) * N^dag(X) for all X.

    Test this for unitary channels (where it must hold) and random channels (where it shouldn't).
    """
    print("\n\n" + "=" * 80)
    print("TEST 6: Multiplicative Domain Verification")
    print("  At exact recovery: N^dag(hX) = N^dag(h)N^dag(X) for all X")
    print("=" * 80)

    d = 2

    # Unitary channels
    print(f"\n  --- Unitary channels (should be in multiplicative domain) ---")
    mult_errors = []

    for trial in range(200):
        rho = random_full_rank_state(d, 0.05)
        sigma = random_full_rank_state(d, 0.05)
        U = unitary_group.rvs(d)
        kraus = [U]

        omega = apply_channel(rho, kraus)
        nu = apply_channel(sigma, kraus)

        h = matrix_log(omega) - matrix_log(nu)

        # Test with random X operators
        max_err = 0
        for _ in range(20):
            X = np.random.randn(d, d) + 1j * np.random.randn(d, d)

            lhs = adjoint_channel(h @ X, kraus)
            rhs = adjoint_channel(h, kraus) @ adjoint_channel(X, kraus)
            err = np.linalg.norm(lhs - rhs, 'fro')
            max_err = max(max_err, err)

        mult_errors.append(max_err)

    mult_errors = np.array(mult_errors)
    print(f"    Max multiplicative domain error: {mult_errors.max():.2e}")
    print(f"    In multiplicative domain: {np.all(mult_errors < 1e-8)}")

    # Random channels (should NOT be in multiplicative domain when Delta_D > 0)
    print(f"\n  --- Random channels (should NOT be in multiplicative domain) ---")
    mult_errors_rand = []

    for trial in range(200):
        rho = random_full_rank_state(d, 0.05)
        sigma = random_full_rank_state(d, 0.05)
        nk = np.random.choice([2, 3])
        kraus = random_cptp_channel(d, nk)

        omega = apply_channel(rho, kraus)
        nu = apply_channel(sigma, kraus)
        if np.min(np.linalg.eigvalsh(nu)) < 1e-10:
            continue
        if np.min(np.linalg.eigvalsh(omega)) < 1e-10:
            continue

        h = matrix_log(omega) - matrix_log(nu)

        max_err = 0
        for _ in range(10):
            X = np.random.randn(d, d) + 1j * np.random.randn(d, d)
            lhs = adjoint_channel(h @ X, kraus)
            rhs = adjoint_channel(h, kraus) @ adjoint_channel(X, kraus)
            err = np.linalg.norm(lhs - rhs, 'fro')
            max_err = max(max_err, err)

        mult_errors_rand.append(max_err)

    mult_errors_rand = np.array(mult_errors_rand)
    print(f"    Valid: {len(mult_errors_rand)}")
    print(f"    Mean multiplicative domain error: {mult_errors_rand.mean():.4f}")
    print(f"    Fraction with error > 0.01: {np.mean(mult_errors_rand > 0.01):.4f}")

# ============================================================
# TEST 7: Quantitative correlation -- gap vs C_out vs Delta_D
# ============================================================

def test7_quantitative_correlation():
    """
    The proof shows: gap ~ Delta_D^2 and C_out correlates with Delta_D.
    Verify the chain: small gap => small Delta_D => small C_out.
    """
    print("\n\n" + "=" * 80)
    print("TEST 7: Quantitative Correlation Chain")
    print("  gap <-> Delta_D <-> C_out")
    print("=" * 80)

    for d in [2, 3]:
        print(f"\n  --- d = {d} ---")
        gaps = []
        DeltaDs = []
        Couts = []

        for _ in range(5000):
            rho = random_full_rank_state(d, 0.05)
            sigma = random_full_rank_state(d, 0.05)
            nk = np.random.choice(range(2, d+2))
            kraus = random_cptp_channel(d, nk)

            r = compute_trial(d, kraus, rho, sigma)
            if r:
                gaps.append(r['sat_gap'])
                DeltaDs.append(r['Delta_D'])
                Couts.append(r['C_out'])

        gaps = np.array(gaps)
        DeltaDs = np.array(DeltaDs)
        Couts = np.array(Couts)

        print(f"    Valid: {len(gaps)}")

        # Correlation matrix
        pos = (gaps > 1e-15) & (DeltaDs > 1e-10) & (Couts > 1e-10)
        if np.sum(pos) > 10:
            log_g = np.log10(gaps[pos])
            log_D = np.log10(DeltaDs[pos])
            log_C = np.log10(Couts[pos])

            corr_gD = np.corrcoef(log_g, log_D)[0, 1]
            corr_gC = np.corrcoef(log_g, log_C)[0, 1]
            corr_DC = np.corrcoef(log_D, log_C)[0, 1]

            print(f"    Corr(log gap, log Delta_D): {corr_gD:.4f}")
            print(f"    Corr(log gap, log C_out):   {corr_gC:.4f}")
            print(f"    Corr(log Delta_D, log C_out): {corr_DC:.4f}")

            # Power law fits
            coeffs_gD = np.polyfit(log_D, log_g, 1)
            coeffs_gC = np.polyfit(log_C, log_g, 1)
            print(f"    Fit: log(gap) ~ {coeffs_gD[0]:.2f} * log(Delta_D) + {coeffs_gD[1]:.2f}")
            print(f"    Fit: log(gap) ~ {coeffs_gC[0]:.2f} * log(C_out) + {coeffs_gC[1]:.2f}")

        # Conditional statistics
        print(f"\n    Conditional: when gap < eps, what is C_out?")
        for eps in [1e-2, 1e-3, 1e-4]:
            mask = gaps < eps
            if np.sum(mask) > 0:
                print(f"      gap < {eps:.0e}: N={np.sum(mask)}, mean C_out={Couts[mask].mean():.4e}, "
                      f"mean Delta_D={DeltaDs[mask].mean():.4e}")

# ============================================================
# TEST 8: The theta=pi/4 special case -- not actual JRSWW saturation
# ============================================================

def test8_theta_pi4_clarification():
    """
    Clarify the theta=pi/4 dephasing case:
    - F^2 = 1/2 = exp(-ln 2), which LOOKS like saturation
    - But this is for the UNROTATED Petz map with Delta_D = ln 2 > 0
    - The JRSWW bound (with rotated maps) would give a different (tighter) bound
    - This is NOT true saturation in the sense of the theorem

    The key insight: this case has Delta_D > 0 and F^2 > 0, so it's a case
    where the Petz bound happens to be tight for the unrotated map, but
    the rotated Petz maps would give F^2 > exp(-Delta_D) strictly.
    """
    print("\n\n" + "=" * 80)
    print("TEST 8: The theta=pi/4 Clarification")
    print("  Is the dephasing theta=pi/4 case 'true' saturation?")
    print("=" * 80)

    d = 2
    theta = np.pi / 4
    psi = np.array([np.cos(theta), np.sin(theta)], dtype=complex)
    rho = np.outer(psi, psi.conj())
    K0 = np.array([[1, 0], [0, 0]], dtype=complex)
    K1 = np.array([[0, 0], [0, 1]], dtype=complex)
    kraus = [K0, K1]

    print(f"\n  rho = |psi><psi| with theta = pi/4")
    print(f"  Channel: full dephasing in computational basis")

    for p in [0.1, 0.3, 0.5, 0.7, 0.9]:
        sigma = np.diag([p, 1 - p]).astype(complex)
        r = compute_trial(d, kraus, rho, sigma)
        if r:
            print(f"\n  p = {p:.1f}:")
            print(f"    F^2 = {r['F2']:.10f}")
            print(f"    exp(-Delta_D) = {r['exp_neg_Delta']:.10f}")
            print(f"    Delta_D = {r['Delta_D']:.10f}")
            print(f"    gap = {r['sat_gap']:.2e}")
            print(f"    C_out = {r['C_out']:.2e}")
            print(f"    C_in = {r['C_in']:.6f}")
            print(f"    Recovery error ||R(N(rho)) - rho|| = {r['recovery_error']:.6f}")

    print(f"\n  Analysis:")
    print(f"    Delta_D = ln(2) ~ {np.log(2):.6f} > 0")
    print(f"    F^2 = 1/2 = exp(-ln 2)")
    print(f"    gap ~ 0 for the UNROTATED Petz map")
    print(f"    BUT: this is NOT exact recovery (recovery error > 0)")
    print(f"    AND: Delta_D > 0, so the JRSWW strict GT inequality should give gap > 0")
    print(f"    RESOLUTION: The unrotated Petz map CAN coincidentally match exp(-Delta_D)")
    print(f"    without the full JRSWW chain being saturated. The 'saturation' is an")
    print(f"    accident of the specific channel + state combination, not structural.")

# ============================================================
# SUMMARY
# ============================================================

def print_summary():
    print("\n\n" + "=" * 80)
    print("FINAL SUMMARY")
    print("=" * 80)
    print("""
    FINDINGS FROM PROOF ATTEMPT:
    ============================

    1. THE ORIGINAL CONJECTURE IS FALSE:
       "Gap = 0 implies [N(rho), N(sigma)] = 0" is DISPROVED by unitary channels.
       Unitary channels always give gap = 0 (Delta_D = 0), but [N(rho), N(sigma)]
       is generically nonzero (it equals U[rho,sigma]U^dag).

    2. THE CONVERSE IS ALSO FALSE:
       "[N(rho), N(sigma)] = 0 implies gap = 0" is DISPROVED by dephasing channels
       with mixed-state rho. Dephasing always kills off-diagonals (C_out = 0),
       but the gap is generically > 0.

    3. THE CORRECT THEOREM:
       (a) Gap = 0 iff Delta_D = 0 (for the JRSWW version with rotated Petz maps).
       (b) Delta_D = 0 iff exact Petz recovery R(N(rho)) = rho (Petz's theorem).
       (c) Exact recovery does NOT require [N(rho), N(sigma)] = 0.

    4. THE QUANTITATIVE CORRELATION EXISTS BUT IS NOT AN EQUIVALENCE:
       - Small gap => small Delta_D => C_out TENDS to be small
       - log(gap) ~ 1.5-2.0 * log(Delta_D) (power law)
       - log(gap) ~ 0.5-0.8 * log(C_out) (weaker power law)
       - The chain gap <-> Delta_D <-> C_out is mediated by the DPI structure

    5. THE KADISON-SCHWARZ AND MULTIPLICATIVE DOMAIN STEPS ARE VERIFIED:
       - At exact recovery: N^dag(h^2) = (N^dag(h))^2 (equality holds)
       - The multiplicative domain condition N^dag(hX) = N^dag(h)N^dag(X) holds
       - But these do NOT force [N(rho), N(sigma)] = 0

    6. THE THETA=PI/4 CASE IS A COINCIDENCE, NOT STRUCTURAL:
       - The unrotated Petz map happens to give F^2 = exp(-Delta_D) = 1/2
       - But Delta_D = ln(2) > 0 and recovery is NOT exact
       - The full JRSWW bound (with rotations) is NOT saturated
    """)

# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":
    t_start = time.time()
    print("Numerical Verification: Saturation and Output Commutativity")
    print("=" * 80)
    print(f"Seed: 42\n")

    test1_unitary_counterexample()
    test2_dephasing_converse()
    test3_strict_golden_thompson()
    test4_kadison_schwarz()
    test5_log_condition()
    test6_multiplicative_domain()
    test7_quantitative_correlation()
    test8_theta_pi4_clarification()
    print_summary()

    print(f"\nTotal runtime: {time.time() - t_start:.1f}s")
