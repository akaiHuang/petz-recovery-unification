#!/usr/bin/env python3
"""
Critical investigation: Does alpha*(d, gamma) -> 0 as gamma -> 0?

If so, then the true infimum is 0 and the conjecture tau >= alpha*C^2
with a UNIVERSAL alpha > 0 would be false. Instead, alpha would be
channel-dependent.

We test very small gamma values and see how alpha*(gamma) behaves.
"""

import numpy as np
from scipy.linalg import sqrtm, svdvals
import time
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

np.random.seed(789)

# Utility functions
def random_pure_state(d):
    psi = np.random.randn(d) + 1j * np.random.randn(d)
    psi /= np.linalg.norm(psi)
    return np.outer(psi, psi.conj())

def random_density_matrix(d, full_rank=True):
    if full_rank:
        A = np.random.randn(d, d) + 1j * np.random.randn(d, d)
        rho = A @ A.conj().T
        rho += 0.01 * np.eye(d)
        rho /= np.trace(rho)
    else:
        rho = random_pure_state(d)
    return rho

def matrix_sqrt_safe(A):
    A = (A + A.conj().T) / 2
    vals, vecs = np.linalg.eigh(A)
    vals = np.maximum(vals, 0)
    return vecs @ np.diag(np.sqrt(vals)) @ vecs.conj().T

def matrix_inv_sqrt_safe(A, eps=1e-12):
    A = (A + A.conj().T) / 2
    vals, vecs = np.linalg.eigh(A)
    vals = np.maximum(vals, eps)
    return vecs @ np.diag(1.0 / np.sqrt(vals)) @ vecs.conj().T

def apply_channel(kraus, rho):
    d = rho.shape[0]
    result = np.zeros((d, d), dtype=complex)
    for K in kraus:
        result += K @ rho @ K.conj().T
    return result

def adjoint_channel(kraus, X):
    d = X.shape[0]
    result = np.zeros((d, d), dtype=complex)
    for K in kraus:
        result += K.conj().T @ X @ K
    return result

def petz_recovery(kraus, sigma, X):
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
    sqrt_rho = matrix_sqrt_safe(rho)
    M = sqrt_rho @ sigma @ sqrt_rho
    M = (M + M.conj().T) / 2
    vals = np.linalg.eigvalsh(M)
    vals = np.maximum(vals, 0)
    F = (np.sum(np.sqrt(vals)))**2
    return min(F.real, 1.0)

def trace_norm(A):
    return np.sum(svdvals(A))

def commutator_trace_norm(A, B):
    return trace_norm(A @ B - B @ A)

def compute_tau_and_C(kraus, rho, sigma):
    N_rho = apply_channel(kraus, rho)
    N_sigma = apply_channel(kraus, sigma)
    C = commutator_trace_norm(N_rho, N_sigma)
    recovered = petz_recovery(kraus, sigma, N_rho)
    F = fidelity(rho, recovered)
    tau = 1.0 - F
    return tau, C

def gen_amplitude_damping_kraus(d, gamma):
    K0 = np.diag([1.0] + [np.sqrt(1 - gamma)] * (d - 1)).astype(complex)
    kraus = [K0]
    for j in range(1, d):
        Kj = np.zeros((d, d), dtype=complex)
        Kj[0, j] = np.sqrt(gamma)
        kraus.append(Kj)
    return kraus

def random_kraus_operators(d, num_kraus=4):
    d_out = d * num_kraus
    A = np.random.randn(d_out, d) + 1j * np.random.randn(d_out, d)
    Q, _ = np.linalg.qr(A)
    V = Q[:, :d]
    kraus = []
    for i in range(num_kraus):
        K_i = V[i*d:(i+1)*d, :]
        kraus.append(K_i)
    check = sum(K.conj().T @ K for K in kraus)
    S = sqrtm(check)
    S_inv = np.linalg.inv(S)
    kraus = [K @ S_inv for K in kraus]
    return kraus


if __name__ == "__main__":
    print()
    print("*" * 74)
    print("*  GAMMA -> 0 LIMIT ANALYSIS                                          *")
    print("*  Does min(tau/C^2) -> 0 as gamma -> 0?                              *")
    print("*" * 74)
    print()

    # Test with geometrically decreasing gamma
    gamma_values = [0.5, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001, 0.00005, 0.00001]
    n_pairs = 2000

    for d in [2, 3, 4, 5]:
        print(f"\n  === d = {d} ===")
        print(f"  {'gamma':>12s}  {'min(tau/C^2)':>14s}  {'min(tau)':>12s}  {'min(C)':>12s}  {'5th pct':>12s}  {'median':>12s}")
        print(f"  {'--------':>12s}  {'----------':>14s}  {'--------':>12s}  {'--------':>12s}  {'--------':>12s}  {'--------':>12s}")

        for gamma in gamma_values:
            kraus = gen_amplitude_damping_kraus(d, gamma)
            ratios = []
            taus = []
            Cs = []
            for _ in range(n_pairs):
                rho = random_pure_state(d)
                sigma = random_density_matrix(d, full_rank=True)
                try:
                    tau, C = compute_tau_and_C(kraus, rho, sigma)
                except:
                    continue
                if C < 1e-14:
                    continue
                if tau < 0:
                    tau = 0.0
                ratios.append(tau / C**2)
                taus.append(tau)
                Cs.append(C)

            ratios = np.array(ratios)
            taus = np.array(taus)
            Cs = np.array(Cs)

            if len(ratios) > 0:
                print(f"  {gamma:12.6f}  {np.min(ratios):14.8f}  {np.min(taus):12.2e}  "
                      f"{np.min(Cs):12.2e}  {np.percentile(ratios,5):12.6f}  {np.median(ratios):12.6f}")
            else:
                print(f"  {gamma:12.6f}  {'N/A':>14s}")

    # KEY ANALYSIS: Perturbative regime
    print()
    print("=" * 74)
    print("  PERTURBATIVE ANALYSIS: tau and C scaling with gamma")
    print("=" * 74)
    print()
    print("  For amplitude damping with small gamma:")
    print("    N(rho) ~ rho + gamma * correction")
    print("    tau ~ O(gamma^a), C ~ O(gamma^b)")
    print("    => tau/C^2 ~ O(gamma^(a-2b))")
    print("    If a >= 2b, then tau/C^2 stays bounded as gamma -> 0")
    print("    If a < 2b, then tau/C^2 -> 0 as gamma -> 0")
    print()

    # Numerical check: measure the scaling exponents
    for d in [2, 3, 4]:
        print(f"\n  --- d = {d}: Measuring tau(gamma) and C(gamma) scaling ---")
        gamma_test = np.array([0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001])

        # Use FIXED states to isolate gamma dependence
        np.random.seed(100 + d)
        rho = random_pure_state(d)
        sigma = random_density_matrix(d, full_rank=True)

        taus_fixed = []
        Cs_fixed = []
        for g in gamma_test:
            kraus = gen_amplitude_damping_kraus(d, g)
            try:
                tau, C = compute_tau_and_C(kraus, rho, sigma)
            except:
                taus_fixed.append(np.nan)
                Cs_fixed.append(np.nan)
                continue
            taus_fixed.append(max(tau, 1e-30))
            Cs_fixed.append(max(C, 1e-30))

        taus_fixed = np.array(taus_fixed)
        Cs_fixed = np.array(Cs_fixed)
        ratios_fixed = taus_fixed / Cs_fixed**2

        # Log-log slopes
        mask = (gamma_test <= 0.01) & ~np.isnan(taus_fixed) & (taus_fixed > 1e-25) & (Cs_fixed > 1e-25)
        if np.sum(mask) >= 2:
            log_g = np.log(gamma_test[mask])
            log_tau = np.log(taus_fixed[mask])
            log_C = np.log(Cs_fixed[mask])

            slope_tau = np.polyfit(log_g, log_tau, 1)[0]
            slope_C = np.polyfit(log_g, log_C, 1)[0]

            print(f"    tau ~ gamma^{slope_tau:.3f}")
            print(f"    C   ~ gamma^{slope_C:.3f}")
            print(f"    => tau/C^2 ~ gamma^({slope_tau:.3f} - 2*{slope_C:.3f}) = gamma^{slope_tau - 2*slope_C:.3f}")
            if slope_tau - 2*slope_C >= -0.1:
                print(f"    => tau/C^2 BOUNDED (stays finite) as gamma -> 0")
            else:
                print(f"    => tau/C^2 -> {'0' if slope_tau - 2*slope_C < 0 else 'inf'} as gamma -> 0")
        else:
            print(f"    Insufficient data for slope fit")

        print(f"\n    gamma        tau              C               tau/C^2")
        for g, t, c, r in zip(gamma_test, taus_fixed, Cs_fixed, ratios_fixed):
            print(f"    {g:10.5f}  {t:14.2e}  {c:14.2e}  {r:14.6f}")

    # Now: what is the TRUE minimum across ALL channel types (not just amp damp)?
    print()
    print("=" * 74)
    print("  CROSS-CHECK: Random channels at small C values")
    print("=" * 74)
    print("  (Channels CLOSE to identity have small C but also small tau)")
    print()

    for d in [2, 3, 4, 5, 6, 7, 8]:
        # Generate channels that are perturbatively close to identity
        n_trials = 5000
        min_ratio = float('inf')
        count_valid = 0

        for _ in range(n_trials):
            rho = random_pure_state(d)
            sigma = random_density_matrix(d, full_rank=True)

            # Near-identity channel: K0 = sqrt(1-eps)*I, K1 = sqrt(eps)*random
            eps = np.random.uniform(0.0001, 0.01)
            K0 = np.sqrt(1 - eps) * np.eye(d, dtype=complex)
            # Random traceless K1
            A = np.random.randn(d, d) + 1j * np.random.randn(d, d)
            A = A / np.linalg.norm(A) * np.sqrt(eps)  # scale
            # Make CPTP: need K0^dag K0 + K1^dag K1 = I
            # K0^dag K0 = (1-eps) I
            # Need K1^dag K1 = eps I => K1 = sqrt(eps) U
            from scipy.stats import unitary_group
            if d >= 2:
                try:
                    U = unitary_group.rvs(d)
                except:
                    U = np.eye(d)
            K1 = np.sqrt(eps) * U
            kraus = [K0, K1]

            try:
                tau, C = compute_tau_and_C(kraus, rho, sigma)
            except:
                continue
            if C < 1e-14:
                continue
            if tau < 0:
                tau = 0.0
            ratio = tau / C**2
            count_valid += 1
            if ratio < min_ratio:
                min_ratio = ratio

        print(f"  d={d}: min(tau/C^2) = {min_ratio:.8f}  (valid={count_valid})")

    print()
    print("*" * 74)
    print("*  FINAL VERDICT                                                       *")
    print("*" * 74)
    print()
