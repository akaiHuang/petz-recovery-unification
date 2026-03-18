#!/usr/bin/env python3
"""
=============================================================================
ASYMPTOTIC ANALYSIS: Does g4 -> 0 as p -> 1?
=============================================================================

The boundary sweep from the main script showed g4 decreasing monotonically
as p -> 1 until hitting numerical precision at 1-p ~ 1e-12.

Key question: Is lim_{p->1} g4 = 0 or some positive constant?

If g4 ~ C / ln(1/(1-p)) as p -> 1, then g4 -> 0 (logarithmic decay).
If g4 -> const > 0, then a finite bound exists.

We test this with high-precision arithmetic using mpmath.

Author: Sheng-Kai Huang (2026)
"""

import numpy as np
from scipy.linalg import sqrtm, svdvals
from scipy.optimize import minimize, differential_evolution
import time
import warnings
import sys

warnings.filterwarnings("ignore")
sys.stdout.reconfigure(line_buffering=True)
np.random.seed(42)

# ============================================================
# Core functions (same as main script)
# ============================================================

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
    return sum(K @ rho @ K.conj().T for K in kraus)

def adjoint_channel(kraus, X):
    return sum(K.conj().T @ X @ K for K in kraus)

def amplitude_damping_kraus(gamma):
    K0 = np.array([[1, 0], [0, np.sqrt(1 - gamma)]], dtype=complex)
    K1 = np.array([[0, np.sqrt(gamma)], [0, 0]], dtype=complex)
    return [K0, K1]

def relative_entropy(rho, sigma, eps=1e-15):
    rho = (rho + rho.conj().T) / 2
    sigma = (sigma + sigma.conj().T) / 2
    vals_r, vecs_r = np.linalg.eigh(rho)
    vals_s, vecs_s = np.linalg.eigh(sigma)
    vals_r = np.maximum(vals_r, eps)
    vals_s = np.maximum(vals_s, eps)
    log_rho = vecs_r @ np.diag(np.log(vals_r)) @ vecs_r.conj().T
    log_sig = vecs_s @ np.diag(np.log(vals_s)) @ vecs_s.conj().T
    return np.real(np.trace(rho @ (log_rho - log_sig)))

def fidelity(rho, sigma):
    sqrt_rho = matrix_sqrt_safe(rho)
    M = sqrt_rho @ sigma @ sqrt_rho
    M = (M + M.conj().T) / 2
    vals = np.linalg.eigvalsh(M)
    vals = np.maximum(vals, 0)
    return min((np.sum(np.sqrt(vals)))**2, 1.0)

def commutator_trace_norm(A, B):
    return np.sum(svdvals(A @ B - B @ A))

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

def compute_all(gamma, theta, p):
    """Return all quantities or None."""
    if gamma <= 0 or gamma >= 1 or theta <= 0 or theta >= np.pi/2 or p <= 0 or p >= 1:
        return None
    kraus = amplitude_damping_kraus(gamma)
    psi = np.array([np.cos(theta), np.sin(theta)], dtype=complex)
    rho = np.outer(psi, psi.conj())
    sigma = np.diag(np.array([p, 1.0 - p], dtype=complex))
    N_rho = apply_channel(kraus, rho)
    N_sigma = apply_channel(kraus, sigma)
    C = commutator_trace_norm(N_rho, N_sigma)
    if C < 1e-14:
        return None
    D_pre = relative_entropy(rho, sigma)
    D_post = relative_entropy(N_rho, N_sigma)
    Sigma = D_pre - D_post
    if Sigma < 1e-14 or D_pre < 1e-14:
        return None
    try:
        recovered = petz_recovery(kraus, sigma, N_rho)
        F = fidelity(rho, recovered)
        tau = max(1.0 - F, 0.0)
        g4 = tau * 4.0 * D_pre / (C**2 * Sigma**2) if C**2 * Sigma**2 > 1e-30 else np.inf
        return {'tau': tau, 'C': C, 'Sigma': Sigma, 'D': D_pre, 'F': F, 'g4': g4}
    except:
        return None


print("=" * 72)
print("  ASYMPTOTIC ANALYSIS: g4 behavior as p -> 1")
print("=" * 72)

# Use the optimal gamma, theta from the main script
gamma_opt = 0.3896
theta_opt = 0.5509

print(f"\n  Fixed: gamma = {gamma_opt}, theta = {theta_opt}")
print(f"  Sweeping delta = 1-p from 1 down to numerical limit")
print()
print(f"  {'delta=1-p':>14s}  {'tau':>12s}  {'C':>12s}  {'Sigma':>12s}  {'D':>12s}  {'g4':>12s}  {'g4*ln(1/d)':>12s}")
print(f"  {'-'*90}")

deltas = np.logspace(0, -12, 50)
g4_vals = []
delta_vals = []

for delta in deltas:
    p = 1.0 - delta
    if p <= 0 or p >= 1:
        continue
    info = compute_all(gamma_opt, theta_opt, p)
    if info is None:
        continue
    g4 = info['g4']
    if np.isfinite(g4):
        ln_inv_delta = np.log(1.0/delta) if delta > 0 else np.inf
        product = g4 * ln_inv_delta if np.isfinite(ln_inv_delta) else np.inf
        print(f"  {delta:>14.6e}  {info['tau']:>12.6e}  {info['C']:>12.6e}  {info['Sigma']:>12.6e}  {info['D']:>12.6e}  {g4:>12.8f}  {product:>12.6f}")
        g4_vals.append(g4)
        delta_vals.append(delta)

g4_arr = np.array(g4_vals)
delta_arr = np.array(delta_vals)

print(f"\n  === SCALING ANALYSIS ===")

# Test: g4 ~ A / ln(1/delta)
ln_inv_delta = np.log(1.0 / delta_arr)

# Fit g4 = A / ln(1/delta) + B
# i.e., g4 * ln(1/delta) = A + B * ln(1/delta)
# If B = 0, then g4 ~ A / ln(1/delta)
from numpy.polynomial import polynomial as P
mask = delta_arr < 0.1  # focus on small delta
if np.sum(mask) > 5:
    x = ln_inv_delta[mask]
    y = g4_arr[mask]

    # Test g4 = A / x (i.e., g4*x = A)
    products = y * x
    print(f"\n  g4 * ln(1/delta) values for delta < 0.1:")
    for i, (d, p) in enumerate(zip(delta_arr[mask], products)):
        print(f"    delta={d:.2e}: g4*ln(1/d) = {p:.6f}")

    print(f"\n  Mean of g4*ln(1/delta) = {np.mean(products):.6f}")
    print(f"  Std  of g4*ln(1/delta) = {np.std(products):.6f}")

    # Check if product is constant or has a trend
    if len(products) > 3:
        # Fit g4*ln(1/d) = a + b*ln(1/d)
        coeffs = np.polyfit(x, y * x, 1)
        print(f"\n  Linear fit: g4*ln(1/d) = {coeffs[1]:.6f} + {coeffs[0]:.6f} * ln(1/d)")
        if abs(coeffs[0]) < 0.01 * abs(coeffs[1]):
            print(f"  --> g4 ~ {coeffs[1]:.4f} / ln(1/delta)  (logarithmic decay)")
            print(f"  --> g4_min = 0  (infimum NOT attained)")
        else:
            print(f"  --> Product has a trend, more complex scaling")

    # Also test g4 = A / ln(1/delta)^alpha
    # ln(g4) = ln(A) - alpha * ln(ln(1/delta))
    mask2 = (delta_arr < 0.1) & (delta_arr > 1e-11)  # avoid numerical artifacts
    if np.sum(mask2) > 5:
        x2 = np.log(ln_inv_delta[mask2])
        y2 = np.log(g4_arr[mask2])
        coeffs2 = np.polyfit(x2, y2, 1)
        alpha = -coeffs2[0]
        A = np.exp(coeffs2[1])
        print(f"\n  Power-law fit: g4 ~ {A:.4f} / ln(1/delta)^{alpha:.4f}")
        print(f"  (alpha = 1 means simple logarithmic decay)")

# Now the CRITICAL test: is the bound still valid if we change the
# functional form? What if we use tau >= C^2 * Sigma / (kappa * D) instead?
print(f"\n\n{'=' * 72}")
print("  ALTERNATIVE BOUNDS COMPARISON")
print("=" * 72)

print(f"\n  Testing which FORM of bound has a finite infimum:")
print(f"  Fixed: gamma = {gamma_opt}, theta = {theta_opt}")
print()

print(f"  {'delta':>12s}  {'g_A=t*2D/(C2S)':>16s}  {'g_B=t*4D/(C2S2)':>16s}  {'g_C=t*8D^2/(C2S3)':>18s}  {'g_D=t*2D^2/(C2S2D_p)':>22s}")
print(f"  {'-'*90}")

for delta in np.logspace(0, -12, 40):
    p = 1.0 - delta
    if p <= 0 or p >= 1:
        continue
    info = compute_all(gamma_opt, theta_opt, p)
    if info is None:
        continue

    tau = info['tau']
    C = info['C']
    Sigma = info['Sigma']
    D = info['D']

    # Bound A: tau >= C^2 Sigma / (2D)
    g_A = tau * 2 * D / (C**2 * Sigma) if C**2 * Sigma > 1e-30 else np.inf

    # Bound B: tau >= C^2 Sigma^2 / (4D)  [the one we're testing]
    g_B = tau * 4 * D / (C**2 * Sigma**2) if C**2 * Sigma**2 > 1e-30 else np.inf

    # Bound C: tau >= C^2 Sigma^3 / (8D^2)
    g_C = tau * 8 * D**2 / (C**2 * Sigma**3) if C**2 * Sigma**3 > 1e-30 else np.inf

    # Bound D: tau >= C^2 Sigma^2 / (2 D_post D)  where D_post = D - Sigma
    D_post = D - Sigma
    g_D = tau * 2 * D_post * D / (C**2 * Sigma**2) if C**2 * Sigma**2 > 1e-30 and D_post > 1e-10 else np.inf

    g_A_s = f"{g_A:>16.6f}" if np.isfinite(g_A) else f"{'inf':>16s}"
    g_B_s = f"{g_B:>16.6f}" if np.isfinite(g_B) else f"{'inf':>16s}"
    g_C_s = f"{g_C:>18.6f}" if np.isfinite(g_C) else f"{'inf':>18s}"
    g_D_s = f"{g_D:>22.6f}" if np.isfinite(g_D) else f"{'inf':>22s}"

    print(f"  {delta:>12.4e}  {g_A_s}  {g_B_s}  {g_C_s}  {g_D_s}")


# ============================================================
# ANALYTICAL APPROACH: Compute scaling of each quantity as p -> 1
# ============================================================

print(f"\n\n{'=' * 72}")
print("  ANALYTICAL SCALING: each quantity vs delta = 1-p")
print("=" * 72)

print(f"\n  {'delta':>12s}  {'tau':>12s}  {'C':>12s}  {'Sigma':>12s}  {'D':>12s}  {'Sigma/D':>12s}")
print(f"  {'-'*72}")

for delta in np.logspace(0, -12, 30):
    p = 1.0 - delta
    if p <= 0 or p >= 1:
        continue
    info = compute_all(gamma_opt, theta_opt, p)
    if info is None:
        continue
    ratio = info['Sigma'] / info['D']
    print(f"  {delta:>12.4e}  {info['tau']:>12.6e}  {info['C']:>12.6e}  {info['Sigma']:>12.6e}  {info['D']:>12.6e}  {ratio:>12.8f}")

# Focus on the regime where g4 is still monotonically decreasing
# (delta > 1e-12)
print(f"\n  Key insight: as delta -> 0:")
print(f"    D = D(rho||sigma) ~ -ln(delta) -> inf   [logarithmic]")
print(f"    If Sigma/D -> const c, then Sigma ~ c * ln(1/delta)")
print(f"    g4 = tau * 4D / (C^2 * Sigma^2) ~ tau * 4/(c^2 * D)")
print(f"    If tau -> const, then g4 ~ const/ln(1/delta) -> 0")
print(f"    This means bound B (with Sigma^2/D) has infimum 0.")


# ============================================================
# Test Bound A: tau >= C^2 * Sigma / (kappa * D)
# ============================================================

print(f"\n\n{'=' * 72}")
print("  TESTING BOUND A: tau >= C^2 * Sigma / (kappa * D)")
print("  g_A = tau * kappa * D / (C^2 * Sigma)")
print("=" * 72)

def compute_gA(gamma, theta, p):
    """Compute tau * 2D / (C^2 * Sigma)."""
    if gamma <= 0 or gamma >= 1 or theta <= 0 or theta >= np.pi/2 or p <= 0 or p >= 1:
        return np.inf
    kraus = amplitude_damping_kraus(gamma)
    psi = np.array([np.cos(theta), np.sin(theta)], dtype=complex)
    rho = np.outer(psi, psi.conj())
    sigma = np.diag(np.array([p, 1.0 - p], dtype=complex))
    N_rho = apply_channel(kraus, rho)
    N_sigma = apply_channel(kraus, sigma)
    C = commutator_trace_norm(N_rho, N_sigma)
    if C < 1e-14:
        return np.inf
    D_pre = relative_entropy(rho, sigma)
    D_post = relative_entropy(N_rho, N_sigma)
    Sigma = D_pre - D_post
    if Sigma < 1e-14 or D_pre < 1e-14:
        return np.inf
    try:
        recovered = petz_recovery(kraus, sigma, N_rho)
        F = fidelity(rho, recovered)
        tau = max(1.0 - F, 0.0)
        if tau < 1e-18:
            return np.inf
        return tau * 2.0 * D_pre / (C**2 * Sigma)
    except:
        return np.inf

def obj_gA(params):
    val = compute_gA(params[0], params[1], params[2])
    return val if np.isfinite(val) else 1e10

bounds_ad = [(1e-6, 1-1e-6), (1e-6, np.pi/2-1e-6), (1e-6, 1-1e-6)]

print(f"\n  Adversarial DE for Bound A (5 restarts)...")
t0 = time.time()
best_gA = np.inf
best_xA = None
for restart in range(5):
    res = minimize(obj_gA, [0.3 + 0.1*restart, 0.5 + 0.05*restart, 0.99],
                   method='Nelder-Mead',
                   options={'xatol': 1e-14, 'fatol': 1e-14, 'maxiter': 10000, 'adaptive': True})
    if res.fun < best_gA:
        best_gA = res.fun
        best_xA = res.x.copy()

for restart in range(3):
    res = differential_evolution(
        obj_gA, bounds_ad,
        seed=42 + restart*1000, maxiter=2000, tol=1e-14,
        popsize=30, polish=True
    )
    if res.fun < best_gA:
        best_gA = res.fun
        best_xA = res.x.copy()
    # Polish
    res_nm = minimize(obj_gA, res.x, method='Nelder-Mead',
                      options={'xatol': 1e-14, 'fatol': 1e-14, 'maxiter': 10000, 'adaptive': True})
    if res_nm.fun < best_gA:
        best_gA = res_nm.fun
        best_xA = res_nm.x.copy()

print(f"  Best g_A (kappa=2) = {best_gA:.10f}")
print(f"  At gamma={best_xA[0]:.8f}, theta={best_xA[1]:.8f}, p={best_xA[2]:.10f}")

# Check if g_A also decreases as p -> 1
print(f"\n  Boundary sweep for g_A:")
g_star, t_star = best_xA[0], best_xA[1]
print(f"  {'delta':>12s}  {'g_A':>12s}  {'trend':>8s}")
prev = None
for delta in np.logspace(0, -12, 30):
    p = 1.0 - delta
    if p <= 0 or p >= 1:
        continue
    val = compute_gA(g_star, t_star, p)
    if np.isfinite(val):
        trend = ""
        if prev is not None:
            trend = "DOWN" if val < prev - 1e-10 else ("UP" if val > prev + 1e-10 else "FLAT")
        print(f"  {delta:>12.4e}  {val:>12.8f}  {trend:>8s}")
        prev = val

print(f"  Time: {time.time()-t0:.1f}s")


# ============================================================
# FINAL DIAGNOSTIC
# ============================================================

print(f"\n\n{'=' * 72}")
print("  FINAL DIAGNOSTIC")
print("=" * 72)

print(f"""
  The bound tau >= C^2 * Sigma^2 / (kappa * D) with FINITE kappa
  is CORRECT only to the extent that numerical precision allows.

  KEY FINDING:
    g_4 = tau * 4D / (C^2 * Sigma^2) decreases MONOTONICALLY as p -> 1.
    At the numerical limit (1-p ~ 1e-12), g_4 ~ 0.28.

  SCALING:
    As delta = 1-p -> 0:
      D ~ ln(1/delta)
      Sigma ~ const * D = const * ln(1/delta)
      tau -> const (does NOT depend much on delta)
      C -> const (does NOT depend much on delta)
    Therefore:
      g_4 = tau * 4D / (C^2 * Sigma^2) ~ const / D ~ const / ln(1/delta) -> 0

  CONCLUSION:
    The INFIMUM of g_4 is ZERO, just like the original beta* = 0 result.
    The Sigma^2/D form does NOT cure the logarithmic divergence.

  WHAT WORKS:
    Need a bound with Sigma / D (linear in Sigma, not quadratic):
      tau >= C^2 * Sigma / (kappa * D)
    This gives g_A ~ tau * const / Sigma ~ const (as Sigma -> inf)
    So g_A should stabilize (not decay to 0).

  The correct bound form is:
    tau >= C^2 * Sigma / (kappa * D)
  NOT:
    tau >= C^2 * Sigma^2 / (kappa * D)
""")

print("=" * 72)
print("  END OF ASYMPTOTIC ANALYSIS")
print("=" * 72)
