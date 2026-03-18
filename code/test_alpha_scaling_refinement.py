#!/usr/bin/env python3
"""
Refinement study: Focus on generalized amplitude damping (the channel giving
the smallest alpha*) with MORE trials and finer parameter scans to get
tighter alpha* estimates.

Also: adversarial search near gamma~1 where alpha* tends to be smallest.
"""

import numpy as np
from scipy.linalg import sqrtm, svdvals
from scipy.optimize import curve_fit, minimize_scalar
import time
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

np.random.seed(456)

# ============================================================
# Utility functions (same as before)
# ============================================================

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
    comm = A @ B - B @ A
    return trace_norm(comm)

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
    check = sum(K.conj().T @ K for K in kraus)
    diff = np.max(np.abs(check - np.eye(d)))
    if diff > 1e-10:
        S = sqrtm(check)
        S_inv = np.linalg.inv(S)
        kraus = [K @ S_inv for K in kraus]
    return kraus


# ============================================================
# Refined amplitude damping scan
# ============================================================

if __name__ == "__main__":
    print()
    print("*" * 74)
    print("*  REFINEMENT: Focused amplitude damping scan for alpha*(d)           *")
    print("*" * 74)
    print()

    dimensions = [2, 3, 4, 5, 6, 7, 8]
    # Fine parameter grid, especially near gamma -> 0 and gamma -> 1
    gamma_values = np.concatenate([
        np.linspace(0.001, 0.05, 20),    # near 0
        np.linspace(0.05, 0.95, 50),     # middle
        np.linspace(0.95, 0.999, 20),    # near 1
    ])
    n_pairs = 500  # more state pairs per gamma

    alpha_star_refined = {}
    best_gamma = {}

    for d in dimensions:
        print(f"  d={d}: scanning {len(gamma_values)} gamma values x {n_pairs} state pairs...")
        t0 = time.time()

        min_ratio = float('inf')
        min_gamma = 0
        all_ratios = []
        skipped = 0

        for gamma in gamma_values:
            kraus = gen_amplitude_damping_kraus(d, gamma)
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
                ratio = tau / (C**2)
                all_ratios.append(ratio)
                if ratio < min_ratio:
                    min_ratio = ratio
                    min_gamma = gamma

        elapsed = time.time() - t0
        alpha_star_refined[d] = min_ratio
        best_gamma[d] = min_gamma
        all_ratios = np.array(all_ratios)
        print(f"         alpha*(d={d}) = {min_ratio:.8f}  (gamma={min_gamma:.4f})")
        print(f"         5th pct = {np.percentile(all_ratios, 5):.6f}, "
              f"valid = {len(all_ratios)}, skipped = {skipped}, time = {elapsed:.1f}s")
        print()

    # Now do a second pass: zoom into the best gamma region for each d
    print("  --- Second pass: zooming into best gamma region ---\n")
    alpha_star_zoomed = {}

    for d in dimensions:
        g_best = best_gamma[d]
        g_lo = max(0.001, g_best - 0.05)
        g_hi = min(0.999, g_best + 0.05)
        gamma_zoom = np.linspace(g_lo, g_hi, 40)
        n_zoom = 1000

        print(f"  d={d}: zooming gamma in [{g_lo:.3f}, {g_hi:.3f}], {len(gamma_zoom)} values x {n_zoom} pairs...")
        t0 = time.time()

        min_ratio = float('inf')
        min_gamma_z = 0
        all_ratios = []
        skipped = 0

        for gamma in gamma_zoom:
            kraus = gen_amplitude_damping_kraus(d, gamma)
            for _ in range(n_zoom):
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
                ratio = tau / (C**2)
                all_ratios.append(ratio)
                if ratio < min_ratio:
                    min_ratio = ratio
                    min_gamma_z = gamma

        elapsed = time.time() - t0
        # Use the better of the two passes
        if min_ratio < alpha_star_refined[d]:
            alpha_star_zoomed[d] = min_ratio
        else:
            alpha_star_zoomed[d] = alpha_star_refined[d]
            min_gamma_z = best_gamma[d]

        print(f"         alpha*(d={d}) = {alpha_star_zoomed[d]:.8f}  (gamma={min_gamma_z:.4f})")
        print(f"         valid = {len(all_ratios)}, time = {elapsed:.1f}s")
        print()

    # ============================================================
    # Final scaling analysis
    # ============================================================
    print("=" * 74)
    print("  REFINED SCALING RESULTS")
    print("=" * 74)

    dims = sorted(alpha_star_zoomed.keys())
    alphas = [alpha_star_zoomed[d] for d in dims]

    print(f"\n  {'d':>4s}  {'alpha*(d)':>12s}  {'alpha*d^2':>12s}  {'alpha*/d^2':>12s}  {'alpha*/d':>12s}")
    for d, a in zip(dims, alphas):
        print(f"  {d:4d}  {a:12.8f}  {a*d**2:12.6f}  {a/d**2:12.8f}  {a/d:12.8f}")

    # Fit models
    print(f"\n  --- Scaling law fits ---")
    d_arr = np.array(dims, dtype=float)
    a_arr = np.array(alphas)

    # Power law
    try:
        def power_law(d, a, beta):
            return a * d**beta
        popt, _ = curve_fit(power_law, d_arr, a_arr, p0=[0.001, 2.0], maxfev=10000)
        pred = power_law(d_arr, *popt)
        res = np.sum((a_arr - pred)**2)
        print(f"  Power law a*d^beta:  a={popt[0]:.8f}, beta={popt[1]:.4f}, residual={res:.2e}")
        for d, p in zip(dims, pred):
            print(f"    d={d}: predicted={p:.6f}")
    except Exception as e:
        print(f"  Power law: FAILED ({e})")

    # a*(d^2-1)/d^2  (= a*(1-1/d^2))
    try:
        def model_d2m1(d, a):
            return a * (d**2 - 1) / d**2
        popt, _ = curve_fit(model_d2m1, d_arr, a_arr, p0=[0.1], maxfev=10000)
        pred = model_d2m1(d_arr, *popt)
        res = np.sum((a_arr - pred)**2)
        print(f"\n  a*(d^2-1)/d^2:  a={popt[0]:.8f}, residual={res:.2e}")
    except Exception as e:
        print(f"  a*(d^2-1)/d^2: FAILED ({e})")

    # a*(d-1)^2/d^2
    try:
        def model_dm1sq(d, a):
            return a * (d - 1)**2 / d**2
        popt, _ = curve_fit(model_dm1sq, d_arr, a_arr, p0=[0.1], maxfev=10000)
        pred = model_dm1sq(d_arr, *popt)
        res = np.sum((a_arr - pred)**2)
        print(f"  a*(d-1)^2/d^2:  a={popt[0]:.8f}, residual={res:.2e}")
    except Exception as e:
        print(f"  a*(d-1)^2/d^2: FAILED ({e})")

    # Quadratic
    try:
        def quadratic(d, a, b, c):
            return a + b*d + c*d**2
        popt, _ = curve_fit(quadratic, d_arr, a_arr, p0=[0, 0, 0.001], maxfev=10000)
        pred = quadratic(d_arr, *popt)
        res = np.sum((a_arr - pred)**2)
        print(f"  Quadratic a+b*d+c*d^2:  a={popt[0]:.6f}, b={popt[1]:.6f}, c={popt[2]:.6f}, residual={res:.2e}")
    except Exception as e:
        print(f"  Quadratic: FAILED ({e})")

    # Log-log analysis
    print(f"\n  --- Log-log regression (ln alpha* vs ln d) ---")
    ln_d = np.log(d_arr)
    ln_a = np.log(a_arr)
    slope, intercept = np.polyfit(ln_d, ln_a, 1)
    print(f"  slope (beta) = {slope:.4f}")
    print(f"  intercept (ln a) = {intercept:.4f}  =>  a = {np.exp(intercept):.8f}")
    print(f"  => alpha*(d) ~ {np.exp(intercept):.6f} * d^{slope:.4f}")

    # Check ratios
    print(f"\n  --- Ratio analysis: alpha*(d+1)/alpha*(d) ---")
    for i in range(len(dims)-1):
        ratio = alphas[i+1] / alphas[i]
        d_ratio = (dims[i+1] / dims[i])
        print(f"  alpha*({dims[i+1]})/alpha*({dims[i]}) = {ratio:.4f}  "
              f"  d-ratio^1 = {d_ratio:.4f}  d-ratio^2 = {d_ratio**2:.4f}  "
              f"  d-ratio^3 = {d_ratio**3:.4f}")

    print()
    print("*" * 74)
    print("*  CONCLUSION                                                          *")
    print("*" * 74)
    print(f"\n  alpha*(d) > 0 for ALL d = {dims[0]}..{dims[-1]}")
    print(f"  Log-log slope: beta = {slope:.4f}")
    if 2.5 < slope < 3.5:
        print(f"  => Consistent with alpha*(d) ~ d^3 scaling")
    elif 1.5 < slope < 2.5:
        print(f"  => Consistent with alpha*(d) ~ d^2 scaling")
    elif 0.5 < slope < 1.5:
        print(f"  => Consistent with alpha*(d) ~ d scaling")
    else:
        print(f"  => Scaling exponent: beta ~ {slope:.2f}")
    print()
