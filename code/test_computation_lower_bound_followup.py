#!/usr/bin/env python3
"""
Follow-up investigation:
1. Amplitude damping gives alpha* ~ 0.001 -- investigate which (gamma, rho, sigma)
   combinations produce the smallest ratios.
2. More aggressive search near the boundary for d=2.
3. Check whether alpha* -> 0 as gamma -> 0 or gamma -> 1 for amplitude damping.
"""

import numpy as np
from scipy.linalg import svdvals
import time
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

np.random.seed(123)

# Import core functions from the main script (inline for self-containedness)

def random_pure_state(d):
    psi = np.random.randn(d) + 1j * np.random.randn(d)
    psi /= np.linalg.norm(psi)
    return np.outer(psi, psi.conj())

def random_density_matrix(d, full_rank=True):
    A = np.random.randn(d, d) + 1j * np.random.randn(d, d)
    rho = A @ A.conj().T
    rho += 0.01 * np.eye(d)
    rho /= np.trace(rho)
    return rho

def random_kraus_operators(d, num_kraus=4):
    d_out = d * num_kraus
    A = np.random.randn(d_out, d) + 1j * np.random.randn(d_out, d)
    Q, _ = np.linalg.qr(A)
    V = Q[:, :d]
    kraus = [V[i*d:(i+1)*d, :] for i in range(num_kraus)]
    from scipy.linalg import sqrtm
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

def amplitude_damping_kraus(gamma):
    K0 = np.array([[1, 0], [0, np.sqrt(1 - gamma)]], dtype=complex)
    K1 = np.array([[0, np.sqrt(gamma)], [0, 0]], dtype=complex)
    return [K0, K1]


# ============================================================
# Investigation 1: Amplitude damping -- which gamma is worst?
# ============================================================

print("=" * 70)
print("INVESTIGATION 1: Amplitude damping -- min(tau/C^2) vs gamma")
print("=" * 70)

gammas = np.linspace(0.001, 0.999, 200)
n_pairs = 2000

min_ratios_by_gamma = []
for gamma in gammas:
    kraus = amplitude_damping_kraus(gamma)
    min_ratio = float('inf')
    for _ in range(n_pairs):
        rho = random_pure_state(2)
        sigma = random_density_matrix(2, full_rank=True)
        N_rho = apply_channel(kraus, rho)
        N_sigma = apply_channel(kraus, sigma)
        C = commutator_trace_norm(N_rho, N_sigma)
        if C < 1e-10:
            continue
        recovered = petz_recovery(kraus, sigma, N_rho)
        F = fidelity(rho, recovered)
        tau = max(1.0 - F, 0.0)
        ratio = tau / (C**2)
        min_ratio = min(min_ratio, ratio)
    min_ratios_by_gamma.append(min_ratio)

min_ratios_by_gamma = np.array(min_ratios_by_gamma)
worst_idx = np.argmin(min_ratios_by_gamma)
print(f"  Worst gamma: {gammas[worst_idx]:.4f}")
print(f"  Min ratio at worst gamma: {min_ratios_by_gamma[worst_idx]:.8f}")
print(f"  Global min ratio: {np.min(min_ratios_by_gamma):.8f}")

# Print some representative values
print(f"\n  gamma vs min(tau/C^2) (selected):")
for i in [0, 10, 20, 50, 100, 150, 180, 190, 199]:
    if i < len(gammas):
        print(f"    gamma={gammas[i]:.3f}: min(tau/C^2) = {min_ratios_by_gamma[i]:.8f}")


# ============================================================
# Investigation 2: Extreme gamma search with more pairs
# ============================================================

print(f"\n{'='*70}")
print("INVESTIGATION 2: Dense search near worst gamma")
print("=" * 70)

# Focus near the worst gamma found above
worst_gamma = gammas[worst_idx]
gammas_fine = np.linspace(max(0.001, worst_gamma - 0.05), min(0.999, worst_gamma + 0.05), 50)
n_pairs_fine = 10000

global_min = float('inf')
global_min_gamma = None

for gamma in gammas_fine:
    kraus = amplitude_damping_kraus(gamma)
    min_ratio = float('inf')
    for _ in range(n_pairs_fine):
        rho = random_pure_state(2)
        sigma = random_density_matrix(2, full_rank=True)
        N_rho = apply_channel(kraus, rho)
        N_sigma = apply_channel(kraus, sigma)
        C = commutator_trace_norm(N_rho, N_sigma)
        if C < 1e-10:
            continue
        recovered = petz_recovery(kraus, sigma, N_rho)
        F = fidelity(rho, recovered)
        tau = max(1.0 - F, 0.0)
        ratio = tau / (C**2)
        min_ratio = min(min_ratio, ratio)

    if min_ratio < global_min:
        global_min = min_ratio
        global_min_gamma = gamma

print(f"  Worst gamma (fine search): {global_min_gamma:.6f}")
print(f"  Global min tau/C^2:        {global_min:.8f}")


# ============================================================
# Investigation 3: Extreme limits gamma -> 0 and gamma -> 1
# ============================================================

print(f"\n{'='*70}")
print("INVESTIGATION 3: Limits gamma -> 0 and gamma -> 1")
print("=" * 70)

extreme_gammas = [0.0001, 0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999, 0.9999]
n_pairs_extreme = 10000

for gamma in extreme_gammas:
    kraus = amplitude_damping_kraus(gamma)
    ratios = []
    for _ in range(n_pairs_extreme):
        rho = random_pure_state(2)
        sigma = random_density_matrix(2, full_rank=True)
        N_rho = apply_channel(kraus, rho)
        N_sigma = apply_channel(kraus, sigma)
        C = commutator_trace_norm(N_rho, N_sigma)
        if C < 1e-10:
            continue
        recovered = petz_recovery(kraus, sigma, N_rho)
        F = fidelity(rho, recovered)
        tau = max(1.0 - F, 0.0)
        ratios.append(tau / (C**2))

    ratios = np.array(ratios)
    print(f"  gamma={gamma:.4f}: min={np.min(ratios):.8f}, "
          f"5th%={np.percentile(ratios,5):.6f}, "
          f"median={np.median(ratios):.6f}, "
          f"valid={len(ratios)}")


# ============================================================
# Investigation 4: Adversarial search -- try to make tau/C^2 small
# ============================================================

print(f"\n{'='*70}")
print("INVESTIGATION 4: Adversarial optimization via gradient-free search")
print("=" * 70)

from scipy.optimize import minimize

def parameterize_state(theta, phi):
    """Parameterize a pure qubit state |psi> = cos(theta/2)|0> + e^{i*phi}*sin(theta/2)|1>"""
    psi = np.array([np.cos(theta/2), np.exp(1j*phi)*np.sin(theta/2)], dtype=complex)
    return np.outer(psi, psi.conj())

def parameterize_density_matrix(a, b_re, b_im):
    """
    Parameterize a full-rank 2x2 density matrix.
    sigma = [[a, b], [b*, 1-a]] with a in (0,1), |b|^2 < a(1-a)
    """
    a_val = 0.01 + 0.98 / (1 + np.exp(-a))  # sigmoid to (0.01, 0.99)
    max_b = np.sqrt(a_val * (1 - a_val)) * 0.99
    b_complex = (b_re + 1j * b_im)
    b_norm = np.abs(b_complex)
    if b_norm > 0:
        b_complex = b_complex / b_norm * max_b * np.tanh(b_norm)
    sigma = np.array([[a_val, b_complex], [b_complex.conj(), 1 - a_val]], dtype=complex)
    return sigma

def objective(params):
    """Minimize tau/C^2 for amplitude damping."""
    gamma = 0.01 + 0.98 / (1 + np.exp(-params[0]))  # sigmoid to (0.01, 0.99)
    theta, phi = params[1], params[2]
    a, b_re, b_im = params[3], params[4], params[5]

    rho = parameterize_state(theta, phi)
    sigma = parameterize_density_matrix(a, b_re, b_im)
    kraus = amplitude_damping_kraus(gamma)

    N_rho = apply_channel(kraus, rho)
    N_sigma = apply_channel(kraus, sigma)
    C = commutator_trace_norm(N_rho, N_sigma)

    if C < 1e-10:
        return 1e10

    recovered = petz_recovery(kraus, sigma, N_rho)
    F = fidelity(rho, recovered)
    tau = max(1.0 - F, 0.0)

    return tau / (C**2)

# Run many random initializations
best_result = float('inf')
best_params = None
n_inits = 500

for i in range(n_inits):
    x0 = np.random.randn(6) * 2
    try:
        result = minimize(objective, x0, method='Nelder-Mead',
                         options={'maxiter': 2000, 'xatol': 1e-10, 'fatol': 1e-12})
        if result.fun < best_result:
            best_result = result.fun
            best_params = result.x
    except:
        pass

# Decode best params
gamma_opt = 0.01 + 0.98 / (1 + np.exp(-best_params[0]))
print(f"  Best tau/C^2 found (adversarial): {best_result:.10f}")
print(f"  At gamma = {gamma_opt:.6f}")


# ============================================================
# Investigation 5: Adversarial search over ALL channels (not just AD)
# ============================================================

print(f"\n{'='*70}")
print("INVESTIGATION 5: Adversarial over random channels")
print("=" * 70)

# Generate many random channels and for each optimize over (rho, sigma)
n_channels = 200
n_opt_inits = 50

global_absolute_min = float('inf')

for ch_idx in range(n_channels):
    kraus = random_kraus_operators(2, num_kraus=4)

    def obj_fixed_channel(params):
        theta, phi = params[0], params[1]
        a, b_re, b_im = params[2], params[3], params[4]
        rho = parameterize_state(theta, phi)
        sigma = parameterize_density_matrix(a, b_re, b_im)
        N_rho = apply_channel(kraus, rho)
        N_sigma = apply_channel(kraus, sigma)
        C = commutator_trace_norm(N_rho, N_sigma)
        if C < 1e-10:
            return 1e10
        recovered = petz_recovery(kraus, sigma, N_rho)
        F = fidelity(rho, recovered)
        tau = max(1.0 - F, 0.0)
        return tau / (C**2)

    local_min = float('inf')
    for _ in range(n_opt_inits):
        x0 = np.random.randn(5) * 2
        try:
            result = minimize(obj_fixed_channel, x0, method='Nelder-Mead',
                             options={'maxiter': 1000, 'xatol': 1e-10, 'fatol': 1e-12})
            local_min = min(local_min, result.fun)
        except:
            pass

    global_absolute_min = min(global_absolute_min, local_min)

    if (ch_idx + 1) % 50 == 0:
        print(f"  ... {ch_idx+1}/{n_channels} channels done, current global min = {global_absolute_min:.10f}")

print(f"\n  Global absolute min tau/C^2 (adversarial, random channels): {global_absolute_min:.10f}")


# ============================================================
# FINAL COMBINED RESULTS
# ============================================================

print(f"\n{'='*70}")
print("COMBINED RESULTS")
print("=" * 70)
print(f"  Amplitude damping adversarial min:  {best_result:.10f}")
print(f"  Random channel adversarial min:     {global_absolute_min:.10f}")
print(f"  Overall minimum alpha*:             {min(best_result, global_absolute_min):.10f}")
print(f"  CONJECTURE STATUS: {'SUPPORTED (alpha* > 0)' if min(best_result, global_absolute_min) > 0 else 'POTENTIALLY VIOLATED'}")
print()
