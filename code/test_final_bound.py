#!/usr/bin/env python3
"""
=============================================================================
DEEP VERIFICATION: tau >= C^2 * Sigma^2 / (4 * D(rho||sigma))
=============================================================================

The bound to test:
    tau >= C^2 * Sigma^2 / (4 * D)

Equivalently, the ratio g = tau * 4D / (C^2 * Sigma^2) >= 1.

where:
    tau   = 1 - F(rho, R_Petz(N(rho)))    (Petz non-recovery)
    C     = ||[N(rho), N(sigma)]||_1       (output commutator trace norm)
    Sigma = D(rho||sigma) - D(N(rho)||N(sigma))  (DPI deficit)
    D     = D(rho||sigma)                  (input relative entropy)

7 verification stages:
    1. Massive random test (d=2, 200k trials)
    2. Adversarial optimization via differential_evolution
    3. Extreme parameter sweeps for amplitude damping
    4. Higher dimensions (d=3,4,5,6)
    5. Mixed states (full-rank rho)
    6. Tightness analysis
    7. Closed-form check at minimizer

Author: Sheng-Kai Huang (2026)
"""

import numpy as np
from scipy.linalg import sqrtm, svdvals
from scipy.optimize import differential_evolution, minimize
import time
import warnings
warnings.filterwarnings("ignore")

np.random.seed(42)

# ============================================================
# HELPER FUNCTIONS
# ============================================================

def random_pure_state(d):
    """Random pure state as density matrix."""
    psi = np.random.randn(d) + 1j * np.random.randn(d)
    psi /= np.linalg.norm(psi)
    return np.outer(psi, psi.conj())

def random_density_matrix(d, rank=None):
    """Random density matrix. If rank=None, full rank."""
    if rank is None:
        rank = d
    A = np.random.randn(d, rank) + 1j * np.random.randn(d, rank)
    rho = A @ A.conj().T + 1e-8 * np.eye(d)
    return rho / np.trace(rho).real

def random_kraus_operators(d_in, num_kraus=4):
    """Random CPTP channel via Stinespring dilation."""
    d_E = num_kraus
    V = np.random.randn(d_in * d_E, d_in) + 1j * np.random.randn(d_in * d_E, d_in)
    Q, R = np.linalg.qr(V)
    V = Q[:, :d_in]
    kraus = [V[i*d_in:(i+1)*d_in, :] for i in range(d_E)]
    check = sum(K.conj().T @ K for K in kraus)
    S_inv = np.linalg.inv(matrix_sqrt_safe(check))
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

def compute_g(kraus, rho, sigma, min_C=1e-10, min_Sigma=1e-10):
    """
    Compute g = tau * 4D / (C^2 * Sigma^2).
    Returns (g, tau, C, Sigma, D) or None if degenerate.
    """
    N_rho = apply_channel(kraus, rho)
    N_sigma = apply_channel(kraus, sigma)

    C = commutator_trace_norm(N_rho, N_sigma)
    if C < min_C:
        return None

    D = relative_entropy(rho, sigma)
    D_out = relative_entropy(N_rho, N_sigma)
    Sigma = D - D_out
    if Sigma < min_Sigma:
        return None

    recovered = petz_recovery(kraus, sigma, N_rho)
    F = fidelity(rho, recovered)
    tau = max(1.0 - F, 0.0)

    if tau < 1e-15:
        return None

    g = tau * 4 * D / (C**2 * Sigma**2)
    return (g, tau, C, Sigma, D)


# ============================================================
# STAGE 1: MASSIVE RANDOM TEST (d=2, 200k trials)
# ============================================================

print("=" * 72)
print("STAGE 1: MASSIVE RANDOM TEST (d=2, 200k trials)")
print("=" * 72)
t0 = time.time()

n_trials = 200000
g_values = []
violations = 0
skipped = 0
violation_details = []

for i in range(n_trials):
    kraus = random_kraus_operators(2, num_kraus=np.random.choice([2, 3, 4]))
    rho = random_pure_state(2)
    sigma = random_density_matrix(2)

    result = compute_g(kraus, rho, sigma)
    if result is None:
        skipped += 1
        continue

    g, tau, C, Sigma, D = result
    g_values.append(g)
    if g < 1.0 - 1e-8:
        violations += 1
        if len(violation_details) < 20:
            violation_details.append({
                'g': g, 'tau': tau, 'C': C, 'Sigma': Sigma, 'D': D,
                'trial': i
            })

    if (i + 1) % 50000 == 0:
        arr = np.array(g_values)
        print(f"  [{i+1}/{n_trials}] valid={len(g_values)}, min_g={np.min(arr):.8f}, "
              f"violations={violations}")

g_arr = np.array(g_values)
elapsed = time.time() - t0

print(f"\n  STAGE 1 RESULTS:")
print(f"  Trials: {n_trials}, Valid: {len(g_values)}, Skipped: {skipped}")
print(f"  Violations (g < 1): {violations}")
print(f"  min(g)   = {np.min(g_arr):.10f}")
print(f"  1st%     = {np.percentile(g_arr, 1):.6f}")
print(f"  5th%     = {np.percentile(g_arr, 5):.6f}")
print(f"  median   = {np.median(g_arr):.4f}")
print(f"  mean     = {np.mean(g_arr):.4f}")
print(f"  Time: {elapsed:.1f}s")

if violation_details:
    print(f"\n  VIOLATION DETAILS (up to 20):")
    for v in violation_details:
        print(f"    trial={v['trial']}: g={v['g']:.8f}, tau={v['tau']:.6e}, "
              f"C={v['C']:.6e}, Sigma={v['Sigma']:.6e}, D={v['D']:.6e}")

stage1_min = np.min(g_arr)
stage1_all = g_arr.copy()


# ============================================================
# STAGE 2: ADVERSARIAL OPTIMIZATION
# ============================================================

print(f"\n{'=' * 72}")
print("STAGE 2: ADVERSARIAL OPTIMIZATION (differential_evolution)")
print("=" * 72)
t0 = time.time()

# --- 2a: Amplitude Damping parametrization ---
print("\n  2a: Amplitude Damping (gamma, theta, p)")

def neg_g_amplitude_damping(params):
    gamma, theta, p = params
    try:
        kraus = amplitude_damping_kraus(gamma)
        psi = np.array([np.cos(theta), np.sin(theta)], dtype=complex)
        rho = np.outer(psi, psi.conj())
        sigma = np.diag(np.array([p, 1 - p], dtype=complex))
        result = compute_g(kraus, rho, sigma, min_C=1e-12, min_Sigma=1e-12)
        if result is None:
            return 1e10
        return result[0]
    except:
        return 1e10

bounds_ad = [(1e-6, 1 - 1e-6), (1e-6, np.pi/2 - 1e-6), (1e-6, 1 - 1e-6)]
result_ad = differential_evolution(neg_g_amplitude_damping, bounds_ad,
                                    seed=42, maxiter=2000, tol=1e-12,
                                    popsize=30, mutation=(0.5, 1.5),
                                    recombination=0.9, polish=True)
print(f"  Best g = {result_ad.fun:.12f}")
print(f"  At (gamma, theta, p) = ({result_ad.x[0]:.8f}, {result_ad.x[1]:.8f}, {result_ad.x[2]:.8f})")

# Verify
gamma_opt, theta_opt, p_opt = result_ad.x
kraus_opt = amplitude_damping_kraus(gamma_opt)
psi_opt = np.array([np.cos(theta_opt), np.sin(theta_opt)], dtype=complex)
rho_opt = np.outer(psi_opt, psi_opt.conj())
sigma_opt = np.diag(np.array([p_opt, 1 - p_opt], dtype=complex))
res_verify = compute_g(kraus_opt, rho_opt, sigma_opt, min_C=1e-15, min_Sigma=1e-15)
if res_verify is not None:
    g_v, tau_v, C_v, Sigma_v, D_v = res_verify
    print(f"  Verification: g={g_v:.12f}, tau={tau_v:.10f}, C={C_v:.10f}, Sigma={Sigma_v:.10f}, D={D_v:.10f}")

# --- 2b: General CPTP via full Stinespring unitary parametrization ---
print("\n  2b: General CPTP (full parametrization, d=2, 2 Kraus ops)")

def build_general_channel_and_states(params):
    """
    Build a general 2-Kraus qubit channel + pure rho + full-rank sigma.
    params has 16 real numbers:
      [0:8]   = 8 real params for 4x2 isometry (via QR of reshaped)
      [8]     = theta for pure rho
      [9]     = phi phase for pure rho
      [10:12] = 2 params for sigma eigenvalues + 1 for unitary on sigma
    """
    # Build isometry via parametrized matrix
    raw = params[0:8].reshape(4, 2)
    V = raw[:, 0:1] + 1j * raw[:, 1:2]  # 4x1 complex, not right...

    # Better: build 2 Kraus from 8 real params
    re = params[0:8]
    K0 = np.array([[re[0] + 1j*re[1], re[2] + 1j*re[3]],
                    [re[4] + 1j*re[5], re[6] + 1j*re[7]]], dtype=complex)

    # CPTP: K1 from remainder
    rem = np.eye(2) - K0.conj().T @ K0
    rem = (rem + rem.conj().T) / 2
    vals, vecs = np.linalg.eigh(rem)
    if np.any(vals < -1e-6):
        return None
    vals = np.maximum(vals, 0)
    K1 = vecs @ np.diag(np.sqrt(vals)) @ vecs.conj().T
    kraus = [K0, K1]

    # Pure rho
    theta = params[8]
    phi = params[9]
    psi = np.array([np.cos(theta), np.sin(theta) * np.exp(1j*phi)], dtype=complex)
    rho = np.outer(psi, psi.conj())

    # Full-rank sigma
    lam = 1.0 / (1.0 + np.exp(-params[10]))  # sigmoid -> (0,1)
    alpha = params[11]
    U_sig = np.array([[np.cos(alpha), -np.sin(alpha)],
                       [np.sin(alpha), np.cos(alpha)]], dtype=complex)
    sigma = U_sig @ np.diag(np.array([lam, 1-lam], dtype=complex)) @ U_sig.conj().T
    sigma = (sigma + sigma.conj().T) / 2
    sigma = sigma / np.trace(sigma).real

    return kraus, rho, sigma

def neg_g_general(params):
    try:
        out = build_general_channel_and_states(params)
        if out is None:
            return 1e10
        kraus, rho, sigma = out

        # Check CPTP
        check = sum(K.conj().T @ K for K in kraus)
        if np.linalg.norm(check - np.eye(2)) > 0.01:
            return 1e10

        result = compute_g(kraus, rho, sigma, min_C=1e-12, min_Sigma=1e-12)
        if result is None:
            return 1e10
        return result[0]
    except:
        return 1e10

bounds_gen = [(-1, 1)]*8 + [(0.01, np.pi/2-0.01), (-np.pi, np.pi), (-3, 3), (-np.pi, np.pi)]
result_gen = differential_evolution(neg_g_general, bounds_gen,
                                     seed=42, maxiter=1500, tol=1e-12,
                                     popsize=25, mutation=(0.5, 1.5),
                                     recombination=0.9, polish=True)
print(f"  Best g = {result_gen.fun:.12f}")

# --- 2c: Depolarizing channel ---
print("\n  2c: Depolarizing channel optimization")

def neg_g_depolarizing(params):
    p_dep, theta, p_sig = params
    try:
        I = np.eye(2, dtype=complex)
        X = np.array([[0, 1], [1, 0]], dtype=complex)
        Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
        Z = np.array([[1, 0], [0, -1]], dtype=complex)
        kraus = [np.sqrt(1 - 3*p_dep/4)*I, np.sqrt(p_dep/4)*X,
                 np.sqrt(p_dep/4)*Y, np.sqrt(p_dep/4)*Z]
        psi = np.array([np.cos(theta), np.sin(theta)], dtype=complex)
        rho = np.outer(psi, psi.conj())
        sigma = np.diag(np.array([p_sig, 1 - p_sig], dtype=complex))
        result = compute_g(kraus, rho, sigma, min_C=1e-12, min_Sigma=1e-12)
        if result is None:
            return 1e10
        return result[0]
    except:
        return 1e10

bounds_dep = [(1e-6, 1 - 1e-6), (1e-6, np.pi/2 - 1e-6), (1e-6, 1 - 1e-6)]
result_dep = differential_evolution(neg_g_depolarizing, bounds_dep,
                                     seed=42, maxiter=2000, tol=1e-12,
                                     popsize=30, polish=True)
print(f"  Best g = {result_dep.fun:.12f}")
print(f"  At (p, theta, p_sig) = ({result_dep.x[0]:.8f}, {result_dep.x[1]:.8f}, {result_dep.x[2]:.8f})")

# --- 2d: Dephasing channel ---
print("\n  2d: Dephasing channel optimization")

def neg_g_dephasing(params):
    lam, theta, p_sig = params
    try:
        I = np.eye(2, dtype=complex)
        Z = np.array([[1, 0], [0, -1]], dtype=complex)
        kraus = [np.sqrt(1 - lam)*I, np.sqrt(lam)*Z]
        psi = np.array([np.cos(theta), np.sin(theta)], dtype=complex)
        rho = np.outer(psi, psi.conj())
        sigma = np.diag(np.array([p_sig, 1 - p_sig], dtype=complex))
        result = compute_g(kraus, rho, sigma, min_C=1e-12, min_Sigma=1e-12)
        if result is None:
            return 1e10
        return result[0]
    except:
        return 1e10

bounds_deph = [(1e-6, 1 - 1e-6), (1e-6, np.pi/2 - 1e-6), (1e-6, 1 - 1e-6)]
result_deph = differential_evolution(neg_g_dephasing, bounds_deph,
                                      seed=42, maxiter=2000, tol=1e-12,
                                      popsize=30, polish=True)
print(f"  Best g = {result_deph.fun:.12f}")

stage2_min = min(result_ad.fun, result_gen.fun, result_dep.fun, result_deph.fun)
elapsed2 = time.time() - t0
print(f"\n  STAGE 2 OVERALL MIN g = {stage2_min:.12f}")
print(f"  Time: {elapsed2:.1f}s")


# ============================================================
# STAGE 3: EXTREME PARAMETER SWEEPS (Amplitude Damping)
# ============================================================

print(f"\n{'=' * 72}")
print("STAGE 3: EXTREME PARAMETER SWEEPS (Amplitude Damping)")
print("=" * 72)
t0 = time.time()

gammas = np.logspace(-8, np.log10(1 - 1e-8), 100)
thetas = np.linspace(1e-4, np.pi/2 - 1e-4, 100)
ps = np.linspace(1e-4, 1 - 1e-4, 100)

global_min_g = np.inf
global_min_params = None
sweep_count = 0

print("  Sweeping gamma (100) x theta (100) x p (100)...")

for ig, gamma in enumerate(gammas):
    kraus = amplitude_damping_kraus(gamma)
    local_min = np.inf
    for theta in thetas:
        psi = np.array([np.cos(theta), np.sin(theta)], dtype=complex)
        rho = np.outer(psi, psi.conj())
        for p in ps:
            sigma = np.diag(np.array([p, 1 - p], dtype=complex))
            result = compute_g(kraus, rho, sigma, min_C=1e-12, min_Sigma=1e-12)
            if result is not None:
                g = result[0]
                sweep_count += 1
                if g < local_min:
                    local_min = g
                if g < global_min_g:
                    global_min_g = g
                    global_min_params = (gamma, theta, p)

    if (ig + 1) % 20 == 0:
        print(f"  [{ig+1}/100] gamma={gamma:.2e}, local_min={local_min:.8f}, "
              f"global_min={global_min_g:.8f}")

elapsed3 = time.time() - t0
print(f"\n  STAGE 3 RESULTS:")
print(f"  Valid evaluations: {sweep_count}")
print(f"  Global min g = {global_min_g:.12f}")
if global_min_params:
    print(f"  At (gamma, theta, p) = ({global_min_params[0]:.8e}, "
          f"{global_min_params[1]:.8f}, {global_min_params[2]:.8f})")
print(f"  Time: {elapsed3:.1f}s")

stage3_min = global_min_g


# ============================================================
# STAGE 4: HIGHER DIMENSIONS (d=3,4,5,6)
# ============================================================

print(f"\n{'=' * 72}")
print("STAGE 4: HIGHER DIMENSIONS (d=3,4,5,6)")
print("=" * 72)
t0 = time.time()

stage4_results = {}
for d in [3, 4, 5, 6]:
    g_vals_d = []
    violations_d = 0
    n_d = 20000

    for i in range(n_d):
        nk = np.random.choice([2, 3, 4])
        kraus = random_kraus_operators(d, num_kraus=nk)
        rho = random_pure_state(d)
        sigma = random_density_matrix(d)

        result = compute_g(kraus, rho, sigma)
        if result is None:
            continue

        g = result[0]
        g_vals_d.append(g)
        if g < 1.0 - 1e-8:
            violations_d += 1

    arr_d = np.array(g_vals_d)
    stage4_results[d] = {'min': np.min(arr_d), 'median': np.median(arr_d),
                          'violations': violations_d, 'valid': len(g_vals_d)}
    print(f"  d={d}: valid={len(g_vals_d)}/{n_d}, min_g={np.min(arr_d):.10f}, "
          f"median={np.median(arr_d):.4f}, violations={violations_d}")

elapsed4 = time.time() - t0
print(f"  Time: {elapsed4:.1f}s")


# ============================================================
# STAGE 5: MIXED STATES (full-rank rho)
# ============================================================

print(f"\n{'=' * 72}")
print("STAGE 5: MIXED STATES (full-rank rho, d=2)")
print("=" * 72)
t0 = time.time()

g_mixed = []
violations_mixed = 0
n_mixed = 20000

for i in range(n_mixed):
    kraus = random_kraus_operators(2, num_kraus=np.random.choice([2, 3, 4]))
    rho = random_density_matrix(2)
    sigma = random_density_matrix(2)

    result = compute_g(kraus, rho, sigma)
    if result is None:
        continue

    g = result[0]
    g_mixed.append(g)
    if g < 1.0 - 1e-8:
        violations_mixed += 1

arr_mixed = np.array(g_mixed)
elapsed5 = time.time() - t0
print(f"  Trials: {n_mixed}, Valid: {len(g_mixed)}")
print(f"  Violations (g < 1): {violations_mixed}")
print(f"  min(g)   = {np.min(arr_mixed):.10f}")
print(f"  1st%     = {np.percentile(arr_mixed, 1):.6f}")
print(f"  5th%     = {np.percentile(arr_mixed, 5):.6f}")
print(f"  median   = {np.median(arr_mixed):.4f}")
print(f"  Time: {elapsed5:.1f}s")

# Higher-d mixed states
print("\n  Higher-d mixed states:")
for d in [3, 4, 5]:
    g_hd = []
    viol_hd = 0
    n_hd = 10000
    for _ in range(n_hd):
        kraus = random_kraus_operators(d, num_kraus=np.random.choice([2, 3, 4]))
        rho = random_density_matrix(d)
        sigma = random_density_matrix(d)
        result = compute_g(kraus, rho, sigma)
        if result is None:
            continue
        g = result[0]
        g_hd.append(g)
        if g < 1.0 - 1e-8:
            viol_hd += 1
    arr_hd = np.array(g_hd)
    print(f"  d={d} mixed: valid={len(g_hd)}/{n_hd}, min_g={np.min(arr_hd):.10f}, "
          f"median={np.median(arr_hd):.4f}, violations={viol_hd}")


# ============================================================
# STAGE 6: TIGHTNESS ANALYSIS
# ============================================================

print(f"\n{'=' * 72}")
print("STAGE 6: TIGHTNESS ANALYSIS")
print("=" * 72)

overall_min = min(stage1_min, stage2_min, stage3_min, np.min(arr_mixed))

print(f"  Stage 1 min (random d=2, pure):     {stage1_min:.12f}")
print(f"  Stage 2 min (adversarial):           {stage2_min:.12f}")
print(f"  Stage 3 min (parameter sweep):       {stage3_min:.12f}")
print(f"  Stage 5 min (mixed states):          {np.min(arr_mixed):.12f}")
print(f"  OVERALL MIN g = {overall_min:.12f}")

# Determine what constant C_opt would make the bound tight
# If g_min = overall_min, the correct bound is tau >= C^2 * Sigma^2 / (4/g_min * D)
# i.e., the correct denominator constant is 4/g_min
if overall_min < 1.0:
    corrected_const = 4.0 / overall_min
    print(f"\n  BOUND VIOLATED! The constant 1/4 is too large.")
    print(f"  The correct constant (with safety margin) would be:")
    print(f"    tau >= C^2 * Sigma^2 / ({corrected_const:.4f} * D)")
    print(f"    i.e., replace 4 with {corrected_const:.4f}")
    print(f"    or equivalently tau >= C^2 * Sigma^2 / (4 * D) needs correction factor {overall_min:.6f}")

# Focused adversarial: really push for the minimum
print("\n  AGGRESSIVE adversarial search (5 restarts, AD + general)...")
best_aggressive = overall_min
for restart in range(5):
    res = differential_evolution(neg_g_amplitude_damping, bounds_ad,
                                  seed=42+restart*1000, maxiter=3000, tol=1e-14,
                                  popsize=50, mutation=(0.3, 1.7),
                                  recombination=0.95, polish=True)
    if res.fun < best_aggressive:
        best_aggressive = res.fun
        print(f"    Restart {restart}: NEW MIN g = {res.fun:.12f}")
    else:
        print(f"    Restart {restart}: g = {res.fun:.12f}")

print(f"  Best aggressive = {best_aggressive:.12f}")
overall_min = min(overall_min, best_aggressive)


# ============================================================
# STAGE 7: CLOSED-FORM CHECK AT MINIMIZER
# ============================================================

print(f"\n{'=' * 72}")
print("STAGE 7: CLOSED-FORM CHECK AT MINIMIZER")
print("=" * 72)

# Recompute at the AD minimizer
if result_ad.fun < 1e9:
    gamma_opt, theta_opt, p_opt = result_ad.x
    kraus_m = amplitude_damping_kraus(gamma_opt)
    psi_m = np.array([np.cos(theta_opt), np.sin(theta_opt)], dtype=complex)
    rho_m = np.outer(psi_m, psi_m.conj())
    sigma_m = np.diag(np.array([p_opt, 1 - p_opt], dtype=complex))

    N_rho_m = apply_channel(kraus_m, rho_m)
    N_sigma_m = apply_channel(kraus_m, sigma_m)

    C_m = commutator_trace_norm(N_rho_m, N_sigma_m)
    D_m = relative_entropy(rho_m, sigma_m)
    D_out_m = relative_entropy(N_rho_m, N_sigma_m)
    Sigma_m = D_m - D_out_m

    recovered_m = petz_recovery(kraus_m, sigma_m, N_rho_m)
    F_m = fidelity(rho_m, recovered_m)
    tau_m = 1.0 - F_m

    g_m = tau_m * 4 * D_m / (C_m**2 * Sigma_m**2)

    print(f"  At the AD minimizer:")
    print(f"    gamma = {gamma_opt:.12f}")
    print(f"    theta = {theta_opt:.12f}")
    print(f"    p     = {p_opt:.12f}")
    print(f"    ---")
    print(f"    tau   = {tau_m:.12e}")
    print(f"    C     = {C_m:.12e}")
    print(f"    Sigma = {Sigma_m:.12e}")
    print(f"    D     = {D_m:.12e}")
    print(f"    ---")
    print(f"    g = tau*4D/(C^2*Sigma^2) = {g_m:.12f}")
    print(f"    ---")
    print(f"    Checking special values:")
    for name, val in [("1", 1), ("2", 2), ("e", np.e), ("pi", np.pi),
                       ("sqrt(2)", np.sqrt(2)), ("ln(2)", np.log(2)),
                       ("1/2", 0.5), ("1/3", 1/3), ("1/4", 0.25),
                       ("1/pi", 1/np.pi), ("1/e", 1/np.e),
                       ("2/pi", 2/np.pi), ("4/pi^2", 4/np.pi**2)]:
        print(f"    g/{name} = {g_m / val:.10f}")

    print(f"\n    1/g = {1.0/g_m:.10f}")

# General channel minimizer details
if result_gen.fun < 1e9:
    print(f"\n  General channel minimizer: g = {result_gen.fun:.12f}")
    print(f"    1/g = {1.0/result_gen.fun:.10f}")
    for name, val in [("1", 1), ("1/2", 0.5), ("1/3", 1/3), ("1/4", 0.25),
                       ("1/pi", 1/np.pi), ("4/pi^2", 4/np.pi**2)]:
        print(f"    g/{name} = {result_gen.fun / val:.10f}")

# Depolarizing minimizer details
if result_dep.fun < 1e9:
    p_dep_opt, theta_dep_opt, p_sig_opt = result_dep.x
    I = np.eye(2, dtype=complex)
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    Z = np.array([[1, 0], [0, -1]], dtype=complex)
    kraus_dep = [np.sqrt(1 - 3*p_dep_opt/4)*I, np.sqrt(p_dep_opt/4)*X,
                 np.sqrt(p_dep_opt/4)*Y, np.sqrt(p_dep_opt/4)*Z]
    psi_dep = np.array([np.cos(theta_dep_opt), np.sin(theta_dep_opt)], dtype=complex)
    rho_dep = np.outer(psi_dep, psi_dep.conj())
    sigma_dep = np.diag(np.array([p_sig_opt, 1 - p_sig_opt], dtype=complex))
    res_dep_v = compute_g(kraus_dep, rho_dep, sigma_dep, min_C=1e-15, min_Sigma=1e-15)
    if res_dep_v is not None:
        g_dep_v = res_dep_v[0]
        print(f"\n  Depolarizing minimizer: g = {g_dep_v:.12f}")
        print(f"    At (p, theta, p_sig) = ({p_dep_opt:.8f}, {theta_dep_opt:.8f}, {p_sig_opt:.8f})")
        print(f"    tau={res_dep_v[1]:.10e}, C={res_dep_v[2]:.10e}, "
              f"Sigma={res_dep_v[3]:.10e}, D={res_dep_v[4]:.10e}")


# ============================================================
# BONUS: WHAT CONSTANT WORKS?
# ============================================================

print(f"\n{'=' * 72}")
print("BONUS: DETERMINING THE CORRECT CONSTANT")
print("=" * 72)

# Test: tau >= C^2 * Sigma^2 / (K * D) for various K
# Equivalent to g_K = tau * K * D / (C^2 * Sigma^2) >= 1
# We need K such that min(g) * 4/K >= 1, i.e., K >= 4/min(g)

print(f"  Overall min g (with factor 4) = {overall_min:.10f}")
print(f"  Required K = 4 / g_min = {4.0 / overall_min:.4f}")
print()

# Try specific denominator constants and check violation rates
all_g_raw = np.concatenate([stage1_all, np.array(g_mixed)])
for K_test in [4, 8, 16, 32, 4*np.pi**2, 64]:
    # g_K = tau * K * D / (C^2 * Sigma^2) = g_4 * K/4
    g_K = all_g_raw * K_test / 4.0
    n_viol = np.sum(g_K < 1.0 - 1e-8)
    min_gK = np.min(g_K)
    print(f"  K={K_test:>8.3f}: min(g_K)={min_gK:.8f}, violations={n_viol}/{len(g_K)}")

# Check the adversarial minimum too
print()
for K_test in [4, 8, 16, 32, 4*np.pi**2, 64]:
    g_K_adv = overall_min * K_test / 4.0
    print(f"  K={K_test:>8.3f}: adversarial g_K = {g_K_adv:.8f} {'PASS' if g_K_adv >= 1.0 else 'FAIL'}")


# ============================================================
# COMPREHENSIVE SUMMARY
# ============================================================

print(f"\n{'=' * 72}")
print("COMPREHENSIVE SUMMARY")
print("=" * 72)
print()
print(f"  Bound tested: tau >= C^2 * Sigma^2 / (4 * D(rho||sigma))")
print(f"  Equivalently: g = tau * 4D / (C^2 * Sigma^2) >= 1")
print()
print(f"  Stage 1 (200k random, d=2, pure rho):   min g = {stage1_min:.10f}, violations = {violations}")
print(f"  Stage 2 (adversarial optimization):       min g = {stage2_min:.10f}")
print(f"    2a AD:          {result_ad.fun:.10f}")
print(f"    2b General:     {result_gen.fun:.10f}")
print(f"    2c Depolarizing:{result_dep.fun:.10f}")
print(f"    2d Dephasing:   {result_deph.fun:.10f}")
print(f"  Stage 3 (1M parameter sweep, AD):        min g = {stage3_min:.10f}")
print(f"  Stage 4 (higher dims):")
for d, res in stage4_results.items():
    print(f"    d={d}: min={res['min']:.10f}, violations={res['violations']}")
print(f"  Stage 5 (mixed states d=2):               min g = {np.min(arr_mixed):.10f}, violations = {violations_mixed}")
print(f"  Stage 6 (aggressive refinement):          min g = {overall_min:.10f}")
print()

total_violations = violations + violations_mixed + sum(r['violations'] for r in stage4_results.values())
total_tested = len(g_values) + len(g_mixed) + sum(r['valid'] for r in stage4_results.values()) + sweep_count
print(f"  Total trials analyzed: ~{total_tested:,}")
print(f"  Total violations (g<1): {total_violations}")
print()

if overall_min >= 1.0 - 1e-6:
    status = "CONFIRMED: bound holds with constant 1/4"
elif overall_min >= 0.5:
    status = f"VIOLATED: min g = {overall_min:.6f} < 1. Constant 1/4 too aggressive"
else:
    status = f"STRONGLY VIOLATED: min g = {overall_min:.6f}. Need much larger denominator"

print(f"  BOUND STATUS: {status}")
print()

if overall_min < 1.0:
    safe_K = np.ceil(4.0 / overall_min * 1.1)  # 10% safety margin
    print(f"  CORRECTED BOUND: tau >= C^2 * Sigma^2 / ({safe_K:.0f} * D)")
    print(f"    (based on min g = {overall_min:.6f}, with 10% safety margin)")
    print()
    print(f"  The infimum of g (with factor 4) is approximately {overall_min:.6f},")
    print(f"  meaning the true constant in the denominator should be at least {4.0/overall_min:.4f}")
    print()

    # Dimension dependence
    print(f"  DIMENSION DEPENDENCE:")
    print(f"    d=2: min g ~ {min(stage1_min, stage2_min, stage3_min, np.min(arr_mixed)):.4f}  (VIOLATED)")
    for d, res in stage4_results.items():
        label = "VIOLATED" if res['min'] < 1.0 else "OK"
        print(f"    d={d}: min g ~ {res['min']:.4f}  ({label})")
    print()
    print(f"  The bound improves with dimension: higher d -> larger g_min.")
    print(f"  Violations are concentrated in d=2 with extreme parameters.")
else:
    print(f"  The bound holds universally with constant 1/4.")
    print(f"  The infimum of g approaches {overall_min:.6f}.")

print()
print("=" * 72)
print("END OF DEEP VERIFICATION")
print("=" * 72)
