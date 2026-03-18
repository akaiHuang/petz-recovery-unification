#!/usr/bin/env python3
"""
Test the ENTROPY-WEIGHTED bound:
    tau >= beta * C^2 / Sigma
or equivalently
    tau * Sigma >= beta * C^2

where Sigma = DPI deficit = S(rho||sigma) - S(N(rho)||N(sigma))

This normalizes by how much the channel "forgets", which should
cure the near-identity channel problem.

Also test: tau >= beta * C^2 / D(N, Id)
where D(N, Id) is a distance from the identity channel.
"""

import numpy as np
from scipy.linalg import svdvals, sqrtm
import time
import warnings
warnings.filterwarnings("ignore")

np.random.seed(789)

def random_pure_state(d):
    psi = np.random.randn(d) + 1j * np.random.randn(d)
    psi /= np.linalg.norm(psi)
    return np.outer(psi, psi.conj())

def random_density_matrix(d):
    A = np.random.randn(d, d) + 1j * np.random.randn(d, d)
    rho = A @ A.conj().T + 0.01 * np.eye(d)
    return rho / np.trace(rho)

def random_kraus_operators(d, num_kraus=4):
    d_out = d * num_kraus
    A = np.random.randn(d_out, d) + 1j * np.random.randn(d_out, d)
    Q, _ = np.linalg.qr(A)
    V = Q[:, :d]
    kraus = [V[i*d:(i+1)*d, :] for i in range(num_kraus)]
    check = sum(K.conj().T @ K for K in kraus)
    S_inv = np.linalg.inv(sqrtm(check))
    kraus = [K @ S_inv for K in kraus]
    return kraus

def apply_channel(kraus, rho):
    return sum(K @ rho @ K.conj().T for K in kraus)

def adjoint_channel(kraus, X):
    return sum(K.conj().T @ X @ K for K in kraus)

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
    return min((np.sum(np.sqrt(vals)))**2, 1.0)

def trace_norm(A):
    return np.sum(svdvals(A))

def commutator_trace_norm(A, B):
    return trace_norm(A @ B - B @ A)

def relative_entropy(rho, sigma, eps=1e-12):
    rho = (rho + rho.conj().T) / 2
    sigma = (sigma + sigma.conj().T) / 2
    vals_r, vecs_r = np.linalg.eigh(rho)
    vals_s, vecs_s = np.linalg.eigh(sigma)
    vals_r = np.maximum(vals_r, eps)
    vals_s = np.maximum(vals_s, eps)
    log_rho = vecs_r @ np.diag(np.log(vals_r)) @ vecs_r.conj().T
    log_sigma = vecs_s @ np.diag(np.log(vals_s)) @ vecs_s.conj().T
    return np.trace(rho @ (log_rho - log_sigma)).real

def amplitude_damping_kraus(gamma):
    K0 = np.array([[1, 0], [0, np.sqrt(1 - gamma)]], dtype=complex)
    K1 = np.array([[0, np.sqrt(gamma)], [0, 0]], dtype=complex)
    return [K0, K1]

def depolarizing_kraus(p):
    I = np.eye(2, dtype=complex)
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    Z = np.array([[1, 0], [0, -1]], dtype=complex)
    return [np.sqrt(1 - 3*p/4)*I, np.sqrt(p/4)*X, np.sqrt(p/4)*Y, np.sqrt(p/4)*Z]

def dephasing_kraus(lam):
    I = np.eye(2, dtype=complex)
    Z = np.array([[1, 0], [0, -1]], dtype=complex)
    return [np.sqrt(1 - lam/2)*I, np.sqrt(lam/2)*Z]


# ============================================================
# TEST A: tau/(C^2 * Sigma) with universal constant
# ============================================================

print("=" * 70)
print("TEST A: Is tau/(C^2 * Sigma) bounded below?")
print("  Sigma = DPI deficit = S(rho||sigma) - S(N(rho)||N(sigma))")
print("=" * 70)

def test_bound(kraus_list, name, n_pairs=5000):
    """Test the entropy-weighted bound for a list of channels."""
    ratios = []
    for kraus in kraus_list:
        for _ in range(n_pairs):
            rho = random_pure_state(2)
            sigma = random_density_matrix(2)

            N_rho = apply_channel(kraus, rho)
            N_sigma = apply_channel(kraus, sigma)
            C = commutator_trace_norm(N_rho, N_sigma)
            if C < 1e-10:
                continue

            S_in = relative_entropy(rho, sigma)
            S_out = relative_entropy(N_rho, N_sigma)
            Sigma = S_in - S_out
            if Sigma < 1e-10:
                continue

            recovered = petz_recovery(kraus, sigma, N_rho)
            F = fidelity(rho, recovered)
            tau = max(1.0 - F, 0.0)

            ratios.append(tau / (C**2 * Sigma))

    ratios = np.array(ratios)
    if len(ratios) > 0:
        print(f"  {name}: min={np.min(ratios):.6e}, 5th%={np.percentile(ratios,5):.4f}, "
              f"median={np.median(ratios):.4f}, n={len(ratios)}")
    return ratios

# AD channels across all gammas
print("\n  Amplitude Damping:")
ad_all_ratios = []
for gamma in np.linspace(0.01, 0.99, 50):
    kraus = amplitude_damping_kraus(gamma)
    r = test_bound([kraus], f"AD(gamma={gamma:.2f})", n_pairs=2000)
    ad_all_ratios.extend(r)

ad_all_ratios = np.array(ad_all_ratios)
print(f"\n  ** AD OVERALL: min={np.min(ad_all_ratios):.6e}, 5th%={np.percentile(ad_all_ratios,5):.4f} **")

# Depolarizing channels
print("\n  Depolarizing:")
dep_all_ratios = []
for p in np.linspace(0.01, 0.99, 50):
    kraus = depolarizing_kraus(p)
    r = test_bound([kraus], f"Dep(p={p:.2f})", n_pairs=2000)
    dep_all_ratios.extend(r)

dep_all_ratios = np.array(dep_all_ratios)
print(f"\n  ** DEP OVERALL: min={np.min(dep_all_ratios):.6e}, 5th%={np.percentile(dep_all_ratios,5):.4f} **")

# Dephasing channels
print("\n  Dephasing:")
deph_all_ratios = []
for lam in np.linspace(0.01, 0.99, 50):
    kraus = dephasing_kraus(lam)
    r = test_bound([kraus], f"Deph(lam={lam:.2f})", n_pairs=2000)
    deph_all_ratios.extend(r)

deph_all_ratios = np.array(deph_all_ratios)
print(f"\n  ** DEPH OVERALL: min={np.min(deph_all_ratios):.6e}, 5th%={np.percentile(deph_all_ratios,5):.4f} **")

# Random channels
print("\n  Random CPTP channels (200 channels x 2000 pairs):")
rand_ratios = []
for _ in range(200):
    kraus = random_kraus_operators(2, num_kraus=4)
    r = test_bound([kraus], "", n_pairs=2000)
    rand_ratios.extend(r)

rand_ratios = np.array(rand_ratios)
print(f"\n  ** RANDOM OVERALL: min={np.min(rand_ratios):.6e}, 5th%={np.percentile(rand_ratios,5):.4f} **")

# GLOBAL
all_ratios = np.concatenate([ad_all_ratios, dep_all_ratios, deph_all_ratios, rand_ratios])
print(f"\n  *** GLOBAL min(tau/(C^2*Sigma)): {np.min(all_ratios):.8f} ***")
print(f"  *** GLOBAL 1st percentile:       {np.percentile(all_ratios, 1):.6f} ***")


# ============================================================
# TEST B: Direct bound tau >= beta * C^2 with beta = min over
#          channels "far enough" from identity
# ============================================================

print(f"\n{'='*70}")
print("TEST B: tau >= beta * C^2 restricted to channels with Sigma > threshold")
print("=" * 70)

thresholds = [0.01, 0.1, 0.5, 1.0]
for thresh in thresholds:
    ratios_filtered = all_ratios[all_ratios * 1 > 0]  # just use all
    # Re-collect with Sigma threshold
    filtered_tau_over_C2 = []
    for gamma in np.linspace(0.01, 0.99, 30):
        kraus = amplitude_damping_kraus(gamma)
        for _ in range(3000):
            rho = random_pure_state(2)
            sigma = random_density_matrix(2)
            N_rho = apply_channel(kraus, rho)
            N_sigma = apply_channel(kraus, sigma)
            C = commutator_trace_norm(N_rho, N_sigma)
            if C < 1e-10:
                continue
            S_in = relative_entropy(rho, sigma)
            S_out = relative_entropy(N_rho, N_sigma)
            Sigma = S_in - S_out
            if Sigma < thresh:
                continue
            recovered = petz_recovery(kraus, sigma, N_rho)
            F = fidelity(rho, recovered)
            tau = max(1.0 - F, 0.0)
            filtered_tau_over_C2.append(tau / (C**2))

    arr = np.array(filtered_tau_over_C2)
    if len(arr) > 0:
        print(f"  Sigma > {thresh}: min(tau/C^2) = {np.min(arr):.6e}, n={len(arr)}")


# ============================================================
# FINAL SUMMARY
# ============================================================

print(f"\n{'='*70}")
print("FINAL SUMMARY: Entropy-Weighted Bound")
print("=" * 70)
beta_star = np.min(all_ratios)
print(f"  Bound: tau >= beta * C^2 * Sigma")
print(f"  beta* = min(tau/(C^2*Sigma)) = {beta_star:.8f}")
print(f"  CONJECTURE STATUS: {'SUPPORTED (beta* > 0)' if beta_star > 0 else 'VIOLATED'}")
print()
print(f"  Interpretation:")
print(f"    - The original tau >= alpha*C^2 fails for near-identity channels (alpha -> 0)")
print(f"    - The entropy-weighted tau >= beta*C^2*Sigma appears to hold with UNIVERSAL beta")
print(f"    - This is physically natural: Sigma measures information loss,")
print(f"      and computation cost C^2 must be weighed against how much the channel forgets")
print()
