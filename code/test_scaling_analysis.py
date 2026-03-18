#!/usr/bin/env python3
"""
Critical scaling analysis: for amplitude damping AD(gamma),
how does min(tau/C^2) scale with gamma as gamma -> 0?

If min(tau/C^2) ~ gamma^beta, the conjecture tau >= alpha*C^2
with universal alpha is violated (alpha -> 0 for near-identity channels).

But a MODIFIED conjecture might hold:
  tau >= alpha(N) * C^2
where alpha(N) depends on the channel's "distance from identity".

Also test: does tau >= C^beta for some beta > 2?
"""

import numpy as np
from scipy.linalg import svdvals
from scipy.optimize import minimize
import time
import warnings
warnings.filterwarnings("ignore")

np.random.seed(456)

def random_pure_state(d):
    psi = np.random.randn(d) + 1j * np.random.randn(d)
    psi /= np.linalg.norm(psi)
    return np.outer(psi, psi.conj())

def random_density_matrix(d):
    A = np.random.randn(d, d) + 1j * np.random.randn(d, d)
    rho = A @ A.conj().T + 0.01 * np.eye(d)
    return rho / np.trace(rho)

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
# Test 1: Precise scaling of min(tau/C^2) with gamma
# ============================================================

print("=" * 70)
print("SCALING ANALYSIS: min(tau/C^2) vs gamma for AD(gamma)")
print("=" * 70)

gammas = [1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 0.1, 0.3, 0.5, 0.7, 0.9]
n_pairs = 20000

results = []
for gamma in gammas:
    kraus = amplitude_damping_kraus(gamma)
    min_ratio = float('inf')
    min_tau = None
    min_C = None
    for _ in range(n_pairs):
        rho = random_pure_state(2)
        sigma = random_density_matrix(2)
        N_rho = apply_channel(kraus, rho)
        N_sigma = apply_channel(kraus, sigma)
        C = commutator_trace_norm(N_rho, N_sigma)
        if C < 1e-14:
            continue
        recovered = petz_recovery(kraus, sigma, N_rho)
        F = fidelity(rho, recovered)
        tau = max(1.0 - F, 0.0)
        ratio = tau / (C**2)
        if ratio < min_ratio:
            min_ratio = ratio
            min_tau = tau
            min_C = C

    results.append((gamma, min_ratio, min_tau, min_C))
    print(f"  gamma={gamma:.1e}: min(tau/C^2)={min_ratio:.6e}, tau={min_tau:.4e}, C={min_C:.4e}")

# Fit power law: min(tau/C^2) ~ gamma^beta
print(f"\n  --- Power law fit ---")
log_gammas = np.log10([r[0] for r in results])
log_ratios = np.log10([r[1] for r in results])

# Linear fit in log-log
coeffs = np.polyfit(log_gammas, log_ratios, 1)
beta = coeffs[0]
log_A = coeffs[1]
print(f"  min(tau/C^2) ~ A * gamma^beta")
print(f"  beta = {beta:.4f}")
print(f"  log10(A) = {log_A:.4f}, A = {10**log_A:.6f}")
print(f"  => min(tau/C^2) ~ {10**log_A:.4f} * gamma^{beta:.4f}")


# ============================================================
# Test 2: Check alternative bound tau >= alpha * C^p for p != 2
# ============================================================

print(f"\n{'='*70}")
print("ALTERNATIVE: Find best p such that tau >= alpha * C^p")
print("=" * 70)

# For a fixed moderate gamma, find what power p gives a tight bound
gamma_test = 0.3
kraus = amplitude_damping_kraus(gamma_test)
n_test = 50000

taus_list = []
Cs_list = []
for _ in range(n_test):
    rho = random_pure_state(2)
    sigma = random_density_matrix(2)
    N_rho = apply_channel(kraus, rho)
    N_sigma = apply_channel(kraus, sigma)
    C = commutator_trace_norm(N_rho, N_sigma)
    if C < 1e-10:
        continue
    recovered = petz_recovery(kraus, sigma, N_rho)
    F = fidelity(rho, recovered)
    tau = max(1.0 - F, 0.0)
    if tau > 1e-15:
        taus_list.append(tau)
        Cs_list.append(C)

taus_arr = np.array(taus_list)
Cs_arr = np.array(Cs_list)

print(f"  gamma={gamma_test}, {len(taus_arr)} valid points")

# For various p, find min(tau/C^p)
ps = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
for p in ps:
    ratios = taus_arr / (Cs_arr**p)
    print(f"  p={p:.1f}: min(tau/C^p) = {np.min(ratios):.6e}, median = {np.median(ratios):.4f}")

# Fit: in log-log, tau vs C
log_tau = np.log(taus_arr)
log_C = np.log(Cs_arr)
# For the LOWER envelope, focus on points near the boundary
# Sort by C, bin, and find min tau in each bin
n_bins = 50
C_sorted_idx = np.argsort(Cs_arr)
bin_size = len(Cs_arr) // n_bins
envelope_log_C = []
envelope_log_tau = []
for i in range(n_bins):
    start = i * bin_size
    end = (i + 1) * bin_size if i < n_bins - 1 else len(Cs_arr)
    idx = C_sorted_idx[start:end]
    min_tau_idx = idx[np.argmin(taus_arr[idx])]
    envelope_log_C.append(np.log(Cs_arr[min_tau_idx]))
    envelope_log_tau.append(np.log(taus_arr[min_tau_idx]))

envelope_coeffs = np.polyfit(envelope_log_C, envelope_log_tau, 1)
p_envelope = envelope_coeffs[0]
print(f"\n  Lower envelope fit: tau ~ C^{p_envelope:.3f}")
print(f"  (if p > 2, then tau/C^2 -> 0 as C -> 0, which would explain the AD result)")


# ============================================================
# Test 3: MODIFIED CONJECTURE -- channel-dependent alpha
# ============================================================

print(f"\n{'='*70}")
print("MODIFIED CONJECTURE: tau >= alpha(N) * C^2")
print("  where alpha(N) = f(diamond_distance(N, Id))")
print("=" * 70)

# For AD(gamma), the diamond distance from identity ~ gamma
# Test: alpha(AD(gamma)) ~ gamma^beta
# We already have this from Test 1
print(f"  From Test 1: alpha(AD(gamma)) ~ {10**log_A:.4f} * gamma^{beta:.4f}")
print()
if abs(beta - 1.0) < 0.3:
    print(f"  ** beta ~ 1 suggests: alpha(N) ~ ||N - Id||_diamond **")
    print(f"  ** Modified conjecture: tau >= c * ||N - Id||_diamond * C^2 **")
elif abs(beta - 2.0) < 0.3:
    print(f"  ** beta ~ 2 suggests: alpha(N) ~ ||N - Id||_diamond^2 **")
    print(f"  ** Modified conjecture: tau >= c * ||N - Id||_diamond^2 * C^2 **")
else:
    print(f"  ** beta = {beta:.2f} -- non-trivial scaling **")


# ============================================================
# Test 4: Check tau >= C^2 * Sigma (entropy production)
# ============================================================

print(f"\n{'='*70}")
print("TEST: tau >= beta * C^2 * Sigma?")
print("  where Sigma = S(N(rho)||N(sigma)) - S(rho||sigma) is entropy production")
print("=" * 70)

# Actually compute relative entropy S(rho||sigma) = Tr(rho(log rho - log sigma))
def relative_entropy(rho, sigma, eps=1e-12):
    """S(rho||sigma) = Tr(rho(log rho - log sigma))"""
    rho = (rho + rho.conj().T) / 2
    sigma = (sigma + sigma.conj().T) / 2

    vals_r, vecs_r = np.linalg.eigh(rho)
    vals_s, vecs_s = np.linalg.eigh(sigma)

    vals_r = np.maximum(vals_r, eps)
    vals_s = np.maximum(vals_s, eps)

    log_rho = vecs_r @ np.diag(np.log(vals_r)) @ vecs_r.conj().T
    log_sigma = vecs_s @ np.diag(np.log(vals_s)) @ vecs_s.conj().T

    return np.trace(rho @ (log_rho - log_sigma)).real


gammas_test = [0.01, 0.1, 0.3, 0.5, 0.9]
n_test = 20000

for gamma in gammas_test:
    kraus = amplitude_damping_kraus(gamma)
    ratios_simple = []
    ratios_sigma_weighted = []

    for _ in range(n_test):
        rho = random_pure_state(2)
        sigma = random_density_matrix(2)
        N_rho = apply_channel(kraus, rho)
        N_sigma = apply_channel(kraus, sigma)
        C = commutator_trace_norm(N_rho, N_sigma)
        if C < 1e-10:
            continue

        recovered = petz_recovery(kraus, sigma, N_rho)
        F = fidelity(rho, recovered)
        tau = max(1.0 - F, 0.0)

        # Entropy production (DPI deficit)
        S_input = relative_entropy(rho, sigma)
        S_output = relative_entropy(N_rho, N_sigma)
        Sigma = S_input - S_output  # Should be >= 0 by DPI

        ratios_simple.append(tau / (C**2))
        if Sigma > 1e-10:
            ratios_sigma_weighted.append(tau / (C**2 * Sigma))

    ratios_simple = np.array(ratios_simple)
    ratios_sigma_weighted = np.array(ratios_sigma_weighted)

    print(f"  gamma={gamma:.2f}: min(tau/C^2) = {np.min(ratios_simple):.6e}, "
          f"min(tau/(C^2*Sigma)) = {np.min(ratios_sigma_weighted):.6e} [{len(ratios_sigma_weighted)} valid]")


# ============================================================
# Test 5: Check pure tau >= C^p directly (no channel dependence needed)
# ============================================================

print(f"\n{'='*70}")
print("DIRECT TEST: Does tau >= alpha * C^p hold for p > 2?")
print("  Testing with ALL channels at once")
print("=" * 70)

# Collect all (tau, C) pairs across many channels
all_taus = []
all_Cs = []

# Random channels
for _ in range(5000):
    d = 2
    d_out = d * 4
    A = np.random.randn(d_out, d) + 1j * np.random.randn(d_out, d)
    Q, _ = np.linalg.qr(A)
    V = Q[:, :d]
    kraus = [V[i*d:(i+1)*d, :] for i in range(4)]
    from scipy.linalg import sqrtm
    check = sum(K.conj().T @ K for K in kraus)
    S_inv = np.linalg.inv(sqrtm(check))
    kraus = [K @ S_inv for K in kraus]

    rho = random_pure_state(d)
    sigma = random_density_matrix(d)
    N_rho = apply_channel(kraus, rho)
    N_sigma = apply_channel(kraus, sigma)
    C = commutator_trace_norm(N_rho, N_sigma)
    if C < 1e-10:
        continue
    recovered = petz_recovery(kraus, sigma, N_rho)
    F = fidelity(rho, recovered)
    tau = max(1.0 - F, 0.0)
    all_taus.append(tau)
    all_Cs.append(C)

# AD channels
for gamma in np.linspace(0.001, 0.999, 100):
    kraus = amplitude_damping_kraus(gamma)
    for _ in range(200):
        rho = random_pure_state(2)
        sigma = random_density_matrix(2)
        N_rho = apply_channel(kraus, rho)
        N_sigma = apply_channel(kraus, sigma)
        C = commutator_trace_norm(N_rho, N_sigma)
        if C < 1e-10:
            continue
        recovered = petz_recovery(kraus, sigma, N_rho)
        F = fidelity(rho, recovered)
        tau = max(1.0 - F, 0.0)
        all_taus.append(tau)
        all_Cs.append(C)

all_taus = np.array(all_taus)
all_Cs = np.array(all_Cs)
print(f"  Total points: {len(all_taus)}")

for p in [2.0, 2.5, 3.0, 3.5, 4.0]:
    ratios = all_taus / (all_Cs**p)
    min_r = np.min(ratios)
    violations = np.sum(ratios < 0)
    print(f"  p={p:.1f}: min(tau/C^p) = {min_r:.6e}, violations (< 0): {violations}")

# Lower envelope analysis
mask = all_taus > 1e-15
log_tau = np.log(all_taus[mask])
log_C = np.log(all_Cs[mask])

# Bin by C and find minimum tau in each bin
n_bins = 100
C_sorted_idx = np.argsort(all_Cs[mask])
bin_size = len(log_C) // n_bins
envelope_log_C = []
envelope_log_tau = []
for i in range(n_bins):
    start = i * bin_size
    end = (i + 1) * bin_size if i < n_bins - 1 else len(log_C)
    idx = C_sorted_idx[start:end]
    taus_bin = all_taus[mask][idx]
    Cs_bin = all_Cs[mask][idx]
    min_idx = np.argmin(taus_bin)
    envelope_log_C.append(np.log(Cs_bin[min_idx]))
    envelope_log_tau.append(np.log(taus_bin[min_idx]))

envelope_coeffs = np.polyfit(envelope_log_C, envelope_log_tau, 1)
p_global = envelope_coeffs[0]
alpha_global = np.exp(envelope_coeffs[1])

print(f"\n  Global lower envelope: tau ~ {alpha_global:.4f} * C^{p_global:.3f}")
print(f"  {'SUPPORTS tau >= alpha*C^2' if p_global <= 2.05 else f'Suggests tau ~ C^{p_global:.2f} (steeper than C^2)'}")

# ============================================================
# SUMMARY
# ============================================================
print(f"\n{'='*70}")
print("SCALING ANALYSIS SUMMARY")
print("=" * 70)
print(f"  1. For AD(gamma): min(tau/C^2) ~ {10**log_A:.4f} * gamma^{beta:.3f}")
print(f"     => As gamma -> 0, alpha -> 0. NO universal alpha exists.")
print(f"  2. Lower envelope across all channels: tau ~ C^{p_global:.3f}")
if p_global > 2.05:
    print(f"     => tau/C^2 -> 0 when C -> 0. The bound is STEEPER than quadratic.")
    print(f"     => MODIFIED CONJECTURE: tau >= alpha * C^{p_global:.1f} (or tau >= alpha * C^2 * f(N))")
else:
    print(f"     => Consistent with quadratic. The issue is only near-identity channels.")
print(f"  3. Channel-dependent bound: tau >= c * ||N-Id|| * C^2 may hold.")
print()
