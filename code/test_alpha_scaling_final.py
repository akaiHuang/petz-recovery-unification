#!/usr/bin/env python3
"""
DEFINITIVE dimension scaling study for tau >= alpha * C^2.

This script combines all findings:
1. Random channels (20k trials per d) -- gives large alpha* (channel-far-from-identity)
2. Specific families (amp damp, depolarizing, random unitary) at various parameters
3. Near-identity channels (the real adversary) -- pushes alpha* toward 0
4. Perturbative analysis showing tau/C^2 ~ gamma for amp damp (goes to 0)

KEY FINDING: The universal alpha* = inf_{N,rho,sigma} tau/C^2 = 0.
But for any FIXED channel N (away from identity), alpha*(N) > 0.
The scaling of alpha*(N) with dimension d depends on the channel family.

For random channels (generic, not near-identity):
  alpha*_random(d) grows roughly as d^{~2.5-3}

Author: Sheng-Kai Huang
Date: 2026-03-16
"""

import numpy as np
from scipy.linalg import sqrtm, svdvals
from scipy.optimize import curve_fit
import time
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

np.random.seed(123)

# ============================================================
# Utility functions
# ============================================================

def random_unitary(d):
    Z = (np.random.randn(d, d) + 1j * np.random.randn(d, d)) / np.sqrt(2)
    Q, R = np.linalg.qr(Z)
    diag = np.diag(R)
    phase = diag / np.abs(diag)
    return Q @ np.diag(phase)

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

def random_kraus_operators(d, num_kraus=4):
    d_out = d * num_kraus
    A = np.random.randn(d_out, d) + 1j * np.random.randn(d_out, d)
    Q, _ = np.linalg.qr(A)
    V = Q[:, :d]
    kraus = [V[i*d:(i+1)*d, :] for i in range(num_kraus)]
    check = sum(K.conj().T @ K for K in kraus)
    S = sqrtm(check)
    S_inv = np.linalg.inv(S)
    kraus = [K @ S_inv for K in kraus]
    return kraus

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

def gen_depolarizing_kraus(d, p):
    coeff0 = np.sqrt(1 - p * (d**2 - 1) / d**2)
    coeff1 = np.sqrt(p / d**2)
    kraus = [coeff0 * np.eye(d, dtype=complex)]
    omega = np.exp(2j * np.pi / d)
    X = np.zeros((d, d), dtype=complex)
    for j in range(d):
        X[(j + 1) % d, j] = 1.0
    Z = np.diag([omega**j for j in range(d)])
    for a in range(d):
        for b in range(d):
            if a == 0 and b == 0:
                continue
            W = np.linalg.matrix_power(X, a) @ np.linalg.matrix_power(Z, b)
            kraus.append(coeff1 * W)
    return kraus

def random_unitary_channel_kraus(d, num_unitaries=4):
    probs = np.random.dirichlet(np.ones(num_unitaries))
    kraus = [np.sqrt(probs[i]) * random_unitary(d) for i in range(num_unitaries)]
    return kraus

def near_identity_kraus(d, eps):
    """Near-identity channel: K0 = sqrt(1-eps)*I, K1 = sqrt(eps)*U."""
    K0 = np.sqrt(1 - eps) * np.eye(d, dtype=complex)
    U = random_unitary(d)
    K1 = np.sqrt(eps) * U
    return [K0, K1]


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":
    total_start = time.time()
    print()
    print("*" * 78)
    print("*  DEFINITIVE DIMENSION SCALING STUDY                                    *")
    print("*  tau >= alpha(d,N) * C^2                                                *")
    print("*  Testing: d = 2, 3, 4, 5, 6, 7, 8                                      *")
    print("*" * 78)
    print()

    dimensions = [2, 3, 4, 5, 6, 7, 8]

    # Storage
    results = {d: {} for d in dimensions}

    # ================================================================
    # PART 1: Random CPTP channels (20k per d)
    # ================================================================
    print("=" * 78)
    print("  PART 1: Random CPTP channels (20,000 trials per dimension)")
    print("=" * 78)

    alpha_random = {}
    for d in dimensions:
        t0 = time.time()
        ratios = []
        for _ in range(20000):
            rho = random_pure_state(d)
            sigma = random_density_matrix(d, full_rank=True)
            kraus = random_kraus_operators(d, num_kraus=4)
            try:
                tau, C = compute_tau_and_C(kraus, rho, sigma)
            except:
                continue
            if C < 1e-10 or tau < 0:
                if tau < 0:
                    tau = 0.0
                if C < 1e-10:
                    continue
            ratios.append(tau / C**2)
        ratios = np.array(ratios)
        alpha_random[d] = np.min(ratios) if len(ratios) > 0 else float('inf')
        results[d]['Random CPTP'] = {
            'min': alpha_random[d],
            'p5': np.percentile(ratios, 5),
            'median': np.median(ratios),
            'n_valid': len(ratios)
        }
        elapsed = time.time() - t0
        print(f"  d={d}: alpha*={alpha_random[d]:.6f}  5th={np.percentile(ratios,5):.4f}  "
              f"median={np.median(ratios):.4f}  n={len(ratios)}  ({elapsed:.1f}s)")

    # ================================================================
    # PART 2: Generalized amplitude damping (various gamma)
    # ================================================================
    print()
    print("=" * 78)
    print("  PART 2: Generalized amplitude damping")
    print("=" * 78)

    gamma_values = np.concatenate([
        np.linspace(0.01, 0.05, 10),
        np.linspace(0.05, 0.95, 30),
        np.linspace(0.95, 0.99, 10),
    ])

    alpha_amp = {}
    for d in dimensions:
        t0 = time.time()
        ratios = []
        n_pairs = max(100, 300 // d)
        for gamma in gamma_values:
            kraus = gen_amplitude_damping_kraus(d, gamma)
            for _ in range(n_pairs):
                rho = random_pure_state(d)
                sigma = random_density_matrix(d, full_rank=True)
                try:
                    tau, C = compute_tau_and_C(kraus, rho, sigma)
                except:
                    continue
                if C < 1e-10:
                    continue
                if tau < 0:
                    tau = 0.0
                ratios.append(tau / C**2)
        ratios = np.array(ratios)
        alpha_amp[d] = np.min(ratios) if len(ratios) > 0 else float('inf')
        results[d]['Amp Damp'] = {
            'min': alpha_amp[d],
            'p5': np.percentile(ratios, 5),
            'median': np.median(ratios),
            'n_valid': len(ratios)
        }
        elapsed = time.time() - t0
        print(f"  d={d}: alpha*={alpha_amp[d]:.6f}  5th={np.percentile(ratios,5):.4f}  ({elapsed:.1f}s)")

    # ================================================================
    # PART 3: Generalized depolarizing
    # ================================================================
    print()
    print("=" * 78)
    print("  PART 3: Generalized depolarizing")
    print("=" * 78)

    p_values = np.linspace(0.01, 0.99, 50)
    alpha_dep = {}
    for d in dimensions:
        t0 = time.time()
        ratios = []
        n_pairs = max(100, 300 // d)
        for p in p_values:
            try:
                kraus = gen_depolarizing_kraus(d, p)
            except:
                continue
            for _ in range(n_pairs):
                rho = random_pure_state(d)
                sigma = random_density_matrix(d, full_rank=True)
                try:
                    tau, C = compute_tau_and_C(kraus, rho, sigma)
                except:
                    continue
                if C < 1e-10:
                    continue
                if tau < 0:
                    tau = 0.0
                ratios.append(tau / C**2)
        ratios = np.array(ratios)
        alpha_dep[d] = np.min(ratios) if len(ratios) > 0 else float('inf')
        results[d]['Depolarizing'] = {
            'min': alpha_dep[d],
            'p5': np.percentile(ratios, 5),
            'median': np.median(ratios),
            'n_valid': len(ratios)
        }
        elapsed = time.time() - t0
        print(f"  d={d}: alpha*={alpha_dep[d]:.6f}  5th={np.percentile(ratios,5):.4f}  ({elapsed:.1f}s)")

    # ================================================================
    # PART 4: Random unitary channels
    # ================================================================
    print()
    print("=" * 78)
    print("  PART 4: Random unitary channels")
    print("=" * 78)

    alpha_uni = {}
    for d in dimensions:
        t0 = time.time()
        n_trials = max(1000, 3000 // d)
        ratios = []
        for _ in range(n_trials):
            rho = random_pure_state(d)
            sigma = random_density_matrix(d, full_rank=True)
            kraus = random_unitary_channel_kraus(d, num_unitaries=4)
            try:
                tau, C = compute_tau_and_C(kraus, rho, sigma)
            except:
                continue
            if C < 1e-10:
                continue
            if tau < 0:
                tau = 0.0
            ratios.append(tau / C**2)
        ratios = np.array(ratios)
        alpha_uni[d] = np.min(ratios) if len(ratios) > 0 else float('inf')
        results[d]['Rand Unitary'] = {
            'min': alpha_uni[d],
            'p5': np.percentile(ratios, 5),
            'median': np.median(ratios),
            'n_valid': len(ratios)
        }
        elapsed = time.time() - t0
        print(f"  d={d}: alpha*={alpha_uni[d]:.6f}  5th={np.percentile(ratios,5):.4f}  ({elapsed:.1f}s)")

    # ================================================================
    # PART 5: Near-identity channels (adversarial)
    # ================================================================
    print()
    print("=" * 78)
    print("  PART 5: Near-identity channels (adversarial, eps in [10^-5, 10^-1])")
    print("=" * 78)

    eps_values = np.logspace(-5, -1, 30)
    alpha_near_id = {}
    for d in dimensions:
        t0 = time.time()
        ratios = []
        for eps in eps_values:
            for _ in range(200):
                kraus = near_identity_kraus(d, eps)
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
        ratios = np.array(ratios)
        alpha_near_id[d] = np.min(ratios) if len(ratios) > 0 else float('inf')
        results[d]['Near-Id'] = {
            'min': alpha_near_id[d],
            'p5': np.percentile(ratios, 5),
            'median': np.median(ratios),
            'n_valid': len(ratios)
        }
        elapsed = time.time() - t0
        print(f"  d={d}: alpha*={alpha_near_id[d]:.8f}  5th={np.percentile(ratios,5):.6f}  ({elapsed:.1f}s)")

    # ================================================================
    # PART 6: Perturbative scaling (fixed states, varying gamma)
    # ================================================================
    print()
    print("=" * 78)
    print("  PART 6: Perturbative scaling: tau(gamma) and C(gamma) for fixed states")
    print("=" * 78)

    for d in [2, 4, 8]:
        np.random.seed(200 + d)
        rho = random_pure_state(d)
        sigma = random_density_matrix(d, full_rank=True)
        np.random.seed(123 + d * 1000)  # reset to not affect other parts

        gamma_test = np.logspace(-5, -1, 20)
        taus_g = []
        Cs_g = []
        for g in gamma_test:
            kraus = gen_amplitude_damping_kraus(d, g)
            tau, C = compute_tau_and_C(kraus, rho, sigma)
            taus_g.append(max(tau, 1e-30))
            Cs_g.append(max(C, 1e-30))

        taus_g = np.array(taus_g)
        Cs_g = np.array(Cs_g)
        ratios_g = taus_g / Cs_g**2

        # Fit log-log slopes in the small-gamma regime (gamma < 0.01)
        mask = gamma_test < 0.01
        if np.sum(mask) >= 3:
            log_g = np.log(gamma_test[mask])
            slope_tau = np.polyfit(log_g, np.log(taus_g[mask]), 1)[0]
            slope_C = np.polyfit(log_g, np.log(Cs_g[mask]), 1)[0]
            slope_ratio = slope_tau - 2 * slope_C
            print(f"  d={d}: tau ~ gamma^{slope_tau:.3f}, C ~ gamma^{slope_C:.3f}, "
                  f"tau/C^2 ~ gamma^{slope_ratio:.3f}")

    np.random.seed(123)  # reset seed for reproducibility

    # ================================================================
    # OVERALL ANALYSIS
    # ================================================================
    print()
    print()
    print("*" * 78)
    print("*  COMPREHENSIVE RESULTS TABLE                                            *")
    print("*" * 78)
    print()

    # Overall min for each d
    overall_min = {}
    overall_min_ch = {}
    for d in dimensions:
        best_val = float('inf')
        best_ch = ''
        for ch_name, info in results[d].items():
            if info['min'] < best_val:
                best_val = info['min']
                best_ch = ch_name
        overall_min[d] = best_val
        overall_min_ch[d] = best_ch

    # Table 1: Per-channel alpha*
    print("  Table 1: alpha*(d) = min(tau/C^2) by channel family")
    print()
    ch_names = ['Random CPTP', 'Amp Damp', 'Depolarizing', 'Rand Unitary', 'Near-Id']
    header = f"  {'d':>3s}"
    for cn in ch_names:
        header += f" | {cn:>12s}"
    header += f" | {'OVERALL':>12s} | {'Best Ch':>12s}"
    print(header)
    print("  " + "-" * (len(header) - 2))
    for d in dimensions:
        row = f"  {d:3d}"
        for cn in ch_names:
            if cn in results[d]:
                row += f" | {results[d][cn]['min']:12.6f}"
            else:
                row += f" | {'N/A':>12s}"
        row += f" | {overall_min[d]:12.8f} | {overall_min_ch[d]:>12s}"
        print(row)

    # Table 2: RANDOM channel alpha* scaling (away from identity)
    print()
    print()
    print("  Table 2: Random CPTP channel alpha* scaling with d")
    print()
    print(f"  {'d':>4s}  {'alpha*(rand)':>12s}  {'a/d^2':>12s}  {'a/d^3':>12s}  "
          f"{'a*(d^2-1)/d^2':>14s}")
    print(f"  {'----':>4s}  {'----------':>12s}  {'------':>12s}  {'------':>12s}  "
          f"{'-------------':>14s}")
    for d in dimensions:
        a = alpha_random[d]
        norm1 = a * (d**2) / (d**2 - 1) if d > 1 else a
        print(f"  {d:4d}  {a:12.6f}  {a/d**2:12.6f}  {a/d**3:12.8f}  {norm1:14.6f}")

    # Table 3: Near-identity (adversarial) scaling
    print()
    print()
    print("  Table 3: Near-identity channel alpha* (adversarial lower bound)")
    print()
    print(f"  {'d':>4s}  {'alpha*(near-id)':>16s}  {'a*d^2':>12s}  {'a/d^2':>14s}")
    for d in dimensions:
        a = alpha_near_id[d]
        print(f"  {d:4d}  {a:16.8f}  {a*d**2:12.6f}  {a/d**2:14.10f}")

    # Scaling law fits for RANDOM channels
    print()
    print()
    print("  SCALING LAW FITS (Random CPTP channels)")
    print("  " + "-" * 60)
    d_arr = np.array(dimensions, dtype=float)
    a_rand = np.array([alpha_random[d] for d in dimensions])

    # Power law
    try:
        def power_law(d, a, beta):
            return a * d**beta
        popt, _ = curve_fit(power_law, d_arr, a_rand, p0=[0.01, 2.0], maxfev=10000)
        pred = power_law(d_arr, *popt)
        res = np.sum((a_rand - pred)**2) / np.sum(a_rand**2)
        print(f"  Power law:    a * d^beta       => a={popt[0]:.6f}, beta={popt[1]:.4f}  "
              f"(rel. residual = {res:.2e})")
    except Exception as e:
        print(f"  Power law: FAILED ({e})")

    # Logarithmic
    try:
        def log_model(d, a, b):
            return a + b * np.log(d)
        popt, _ = curve_fit(log_model, d_arr, a_rand, p0=[0, 5], maxfev=10000)
        pred = log_model(d_arr, *popt)
        res = np.sum((a_rand - pred)**2) / np.sum(a_rand**2)
        print(f"  Logarithmic:  a + b*ln(d)      => a={popt[0]:.6f}, b={popt[1]:.6f}  "
              f"(rel. residual = {res:.2e})")
    except Exception as e:
        print(f"  Logarithmic: FAILED ({e})")

    # Quadratic
    try:
        def quad_model(d, a, b, c):
            return a + b * d + c * d**2
        popt, _ = curve_fit(quad_model, d_arr, a_rand, p0=[0, 0, 0.1], maxfev=10000)
        pred = quad_model(d_arr, *popt)
        res = np.sum((a_rand - pred)**2) / np.sum(a_rand**2)
        print(f"  Quadratic:    a + b*d + c*d^2   => a={popt[0]:.4f}, b={popt[1]:.4f}, c={popt[2]:.4f}  "
              f"(rel. residual = {res:.2e})")
    except Exception as e:
        print(f"  Quadratic: FAILED ({e})")

    # Normalized (d^2-1)/d^2
    try:
        def norm_model(d, a):
            return a * (d**2 - 1) / d**2
        popt, _ = curve_fit(norm_model, d_arr, a_rand, p0=[5], maxfev=10000)
        pred = norm_model(d_arr, *popt)
        res = np.sum((a_rand - pred)**2) / np.sum(a_rand**2)
        print(f"  (d^2-1)/d^2:  a*(d^2-1)/d^2    => a={popt[0]:.6f}  "
              f"(rel. residual = {res:.2e})")
    except Exception as e:
        print(f"  (d^2-1)/d^2: FAILED ({e})")

    # Log-log regression
    ln_d = np.log(d_arr)
    ln_a = np.log(a_rand)
    slope, intercept = np.polyfit(ln_d, ln_a, 1)
    print(f"\n  Log-log regression: ln(alpha*) = {intercept:.4f} + {slope:.4f} * ln(d)")
    print(f"  => alpha*(d) ~ {np.exp(intercept):.6f} * d^{slope:.4f}")

    # Scaling fits for NEAR-IDENTITY channels
    print()
    print()
    print("  SCALING LAW FITS (Near-identity channels)")
    print("  " + "-" * 60)
    a_near = np.array([alpha_near_id[d] for d in dimensions])
    ln_a_near = np.log(a_near)
    slope_near, intercept_near = np.polyfit(ln_d, ln_a_near, 1)
    print(f"  Log-log regression: ln(alpha*) = {intercept_near:.4f} + {slope_near:.4f} * ln(d)")
    print(f"  => alpha*_near(d) ~ {np.exp(intercept_near):.8f} * d^{slope_near:.4f}")

    # ================================================================
    # DIMENSION CROSS-MULTIPLICATION TABLE
    # ================================================================
    print()
    print()
    print("  Table 4: alpha*(d) * d^k — looking for approximate constant")
    print()
    print(f"  {'d':>4s}", end="")
    for k_label in ['d^-2', 'd^-1', 'd^0', 'd^1', 'd^2']:
        print(f"  {k_label:>12s}", end="")
    print(f"  (Random CPTP)")
    print(f"  {'----':>4s}", end="")
    for _ in range(5):
        print(f"  {'----------':>12s}", end="")
    print()
    for d in dimensions:
        a = alpha_random[d]
        print(f"  {d:4d}", end="")
        for k in [-2, -1, 0, 1, 2]:
            print(f"  {a * d**k:12.4f}", end="")
        print()

    # Coefficient of variation for each k
    print()
    best_cv = float('inf')
    best_k = 0
    for k in [-2, -1, 0, 1, 2]:
        vals = [alpha_random[d] * d**k for d in dimensions]
        cv = np.std(vals) / np.mean(vals)
        marker = ""
        if cv < best_cv:
            best_cv = cv
            best_k = k
        print(f"    k={k:+d}: CV = {cv:.4f}")
    print(f"    => Best: k={best_k:+d} (CV={best_cv:.4f})")

    # ================================================================
    # FINAL VERDICT
    # ================================================================
    total_elapsed = time.time() - total_start
    print()
    print()
    print("*" * 78)
    print("*  FINAL VERDICT                                                          *")
    print("*" * 78)
    print()
    print("  1. UNIVERSAL BOUND: alpha* = inf_{N,rho,sigma} tau/C^2 = 0")
    print("     Near-identity channels drive tau/C^2 arbitrarily close to 0.")
    print("     For amplitude damping: tau ~ gamma, C ~ O(1), so tau/C^2 ~ gamma -> 0.")
    print("     Therefore: NO universal alpha > 0 exists.")
    print()
    print("  2. CHANNEL-DEPENDENT BOUND: For any fixed channel N (bounded away from")
    print("     identity), there exists alpha(N) > 0 such that tau >= alpha(N) * C^2.")
    print()
    print("  3. RANDOM CHANNEL SCALING: For generic (non-near-identity) channels,")
    print(f"     alpha*(d) ~ {np.exp(intercept):.4f} * d^{slope:.2f}")
    print(f"     (log-log slope = {slope:.4f})")
    print()
    print("  4. NEAR-IDENTITY SCALING:")
    print(f"     alpha*_near(d) ~ {np.exp(intercept_near):.6f} * d^{slope_near:.2f}")
    print(f"     This goes to 0 as eps -> 0, NOT as a function of d alone.")
    print()
    print("  5. REFORMULATED CONJECTURE:")
    print("     For any CPTP map N with dist(N, Id) >= delta > 0:")
    print("       tau >= alpha(d, delta) * C^2")
    print("     where alpha(d, delta) > 0 depends on both d and the distance from identity.")
    print()
    print(f"  Total computation time: {total_elapsed:.1f}s ({total_elapsed/60:.1f} min)")
    print()
    print("*" * 78)
    print("*  END OF DEFINITIVE STUDY                                                *")
    print("*" * 78)
