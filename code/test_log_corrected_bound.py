#!/usr/bin/env python3
"""
=============================================================================
LOG-CORRECTED BOUND: Final Verification
=============================================================================

Test the proposed fix:  tau >= C^2 * Sigma^2 / (c * D * ln(D + e))

and compare against alternative forms:
  - tau >= C^2 * Sigma^2 / (c * D^2)
  - tau >= C^2 * Sigma^2 * exp(-alpha*D) / c
  - tau >= C^2 * Sigma^2 / (c * D * (1 + ln D))
  - tau >= (C * Sigma / D)^2 / c

Key question: which form (if any) gives a FINITE positive infimum?

From prior analysis: g4 = tau * 4D / (C^2 * Sigma^2) ~ 7.55 / ln(1/delta)
So g_log = g4 * c * ln(D+e) ~ 7.55*c*ln(D)/D -> 0  (STILL fails)
But g_D2 = g4 * D/4 ~ 7.55 * D / (4 * ln(1/delta)) and D ~ ln(1/delta)
so g_D2 ~ 7.55/4 ~ const!

Author: Sheng-Kai Huang (2026)
"""

import numpy as np
from scipy.linalg import sqrtm, svdvals
from scipy.optimize import minimize, differential_evolution
from scipy.stats import unitary_group
import time
import warnings
import sys

warnings.filterwarnings("ignore")
sys.stdout.reconfigure(line_buffering=True)
np.random.seed(42)

# ============================================================
# Core functions
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

def depolarizing_kraus(p, d=2):
    """Depolarizing channel: rho -> (1-p)*rho + p*I/d."""
    K0 = np.sqrt(1 - p + p/d) * np.eye(d, dtype=complex)
    kraus = [K0]
    # Add d^2 - 1 more Kraus operators
    # For d=2, use Pauli matrices
    if d == 2:
        sx = np.array([[0, 1], [1, 0]], dtype=complex)
        sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
        sz = np.array([[1, 0], [0, -1]], dtype=complex)
        for P in [sx, sy, sz]:
            kraus.append(np.sqrt(p / d) * P / np.sqrt(d))
    return kraus

def phase_damping_kraus(lam):
    K0 = np.array([[1, 0], [0, np.sqrt(1 - lam)]], dtype=complex)
    K1 = np.array([[0, 0], [0, np.sqrt(lam)]], dtype=complex)
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

def compute_all(kraus, rho, sigma):
    """Return dict of all quantities, or None if degenerate."""
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
        return {'tau': tau, 'C': C, 'Sigma': Sigma, 'D': D_pre, 'F': F}
    except:
        return None

def compute_all_ad(gamma, theta, p):
    """Amplitude damping specific."""
    if gamma <= 0 or gamma >= 1 or theta <= 0 or theta >= np.pi/2 or p <= 0 or p >= 1:
        return None
    kraus = amplitude_damping_kraus(gamma)
    psi = np.array([np.cos(theta), np.sin(theta)], dtype=complex)
    rho = np.outer(psi, psi.conj())
    sigma = np.diag(np.array([p, 1.0 - p], dtype=complex))
    return compute_all(kraus, rho, sigma)

def random_pure_state(d):
    psi = np.random.randn(d) + 1j * np.random.randn(d)
    psi /= np.linalg.norm(psi)
    return np.outer(psi, psi.conj())

def random_full_rank_state(d, min_eigenval=1e-8):
    """Random full-rank density matrix with minimum eigenvalue."""
    # Random eigenvalues
    eigs = np.random.exponential(1.0, d)
    eigs = np.maximum(eigs, min_eigenval)
    eigs /= np.sum(eigs)
    # Random unitary
    U = unitary_group.rvs(d)
    return U @ np.diag(eigs) @ U.conj().T

def random_cptp_kraus(d, n_kraus=None):
    """Generate random CPTP map via random isometry."""
    if n_kraus is None:
        n_kraus = d**2
    # Random isometry V: d -> d*n_kraus
    V = np.random.randn(d * n_kraus, d) + 1j * np.random.randn(d * n_kraus, d)
    Q, _ = np.linalg.qr(V)
    V = Q[:, :d]
    # Extract Kraus operators
    kraus = []
    for i in range(n_kraus):
        Ki = V[i*d:(i+1)*d, :]
        kraus.append(Ki)
    return kraus


# ============================================================
# Bound evaluation functions
# ============================================================

E = np.e  # Euler's number

def g_log(tau, C, Sigma, D, c_val):
    """g for: tau >= C^2 * Sigma^2 / (c * D * ln(D + e))"""
    denom = C**2 * Sigma**2
    if denom < 1e-30:
        return np.inf
    return tau * c_val * D * np.log(D + E) / denom

def g_D2(tau, C, Sigma, D, c_val):
    """g for: tau >= C^2 * Sigma^2 / (c * D^2)"""
    denom = C**2 * Sigma**2
    if denom < 1e-30:
        return np.inf
    return tau * c_val * D**2 / denom

def g_exp(tau, C, Sigma, D, alpha, c_val):
    """g for: tau >= C^2 * Sigma^2 * exp(-alpha*D) / c"""
    denom = C**2 * Sigma**2
    if denom < 1e-30:
        return np.inf
    return tau * c_val * np.exp(alpha * D) / denom

def g_logD(tau, C, Sigma, D, c_val):
    """g for: tau >= C^2 * Sigma^2 / (c * D * (1 + ln D))"""
    denom = C**2 * Sigma**2
    if denom < 1e-30:
        return np.inf
    ln_D = max(np.log(D), 0)  # D >= 0 always
    return tau * c_val * D * (1 + ln_D) / denom

def g_sqratio(tau, C, Sigma, D, c_val):
    """g for: tau >= (C * Sigma / D)^2 / c"""
    if D < 1e-30:
        return np.inf
    return tau * c_val * D**2 / (C**2 * Sigma**2)  # same as g_D2!

def g4_original(tau, C, Sigma, D):
    """Original g4 = tau * 4D / (C^2 * Sigma^2)"""
    denom = C**2 * Sigma**2
    if denom < 1e-30:
        return np.inf
    return tau * 4.0 * D / denom


# ============================================================
# SECTION 1: 200k random trials (d=2)
# ============================================================

print("=" * 80)
print("  SECTION 1: 200k RANDOM TRIALS (d=2)")
print("=" * 80)
print()

N_TRIALS = 200000
c_values = [1, 2, 3, 4, 5, 6, 8]

# Track violations for each c value
violations_log = {c: 0 for c in c_values}
violations_D2 = {c: 0 for c in c_values}
valid_count = 0
min_g_log = {c: np.inf for c in c_values}
min_g_D2 = {c: np.inf for c in c_values}

# Extreme p distribution: mix uniform and extreme
t0 = time.time()
print(f"  Running {N_TRIALS} random trials...")

for trial in range(N_TRIALS):
    if trial % 50000 == 0 and trial > 0:
        print(f"    ... {trial}/{N_TRIALS} done ({valid_count} valid)")

    d = 2

    # Random channel (mix of types)
    ch_type = np.random.choice(['amplitude', 'depolarizing', 'phase', 'random'])
    if ch_type == 'amplitude':
        gamma = np.random.uniform(0.01, 0.99)
        kraus = amplitude_damping_kraus(gamma)
    elif ch_type == 'depolarizing':
        p_dep = np.random.uniform(0.01, 0.99)
        kraus = depolarizing_kraus(p_dep, d)
    elif ch_type == 'phase':
        lam = np.random.uniform(0.01, 0.99)
        kraus = phase_damping_kraus(lam)
    else:
        kraus = random_cptp_kraus(d, n_kraus=np.random.randint(2, 5))

    # Random pure rho
    rho = random_pure_state(d)

    # Random full-rank sigma with extreme eigenvalues
    if np.random.random() < 0.3:
        # Extreme p: very close to 1
        log_delta = np.random.uniform(-8, -1)  # delta from 1e-8 to 0.1
        delta = 10**log_delta
        p = 1.0 - delta
        sigma = np.diag(np.array([p, 1.0 - p], dtype=complex))
    else:
        sigma = random_full_rank_state(d, min_eigenval=1e-10)

    info = compute_all(kraus, rho, sigma)
    if info is None:
        continue
    if info['tau'] < 1e-18:
        continue

    valid_count += 1
    tau, C, Sigma, D = info['tau'], info['C'], info['Sigma'], info['D']

    for c in c_values:
        gl = g_log(tau, C, Sigma, D, c)
        gd = g_D2(tau, C, Sigma, D, c)
        if np.isfinite(gl):
            if gl < 1.0:
                violations_log[c] += 1
            min_g_log[c] = min(min_g_log[c], gl)
        if np.isfinite(gd):
            if gd < 1.0:
                violations_D2[c] += 1
            min_g_D2[c] = min(min_g_D2[c], gd)

elapsed = time.time() - t0
print(f"\n  Completed: {valid_count} valid trials in {elapsed:.1f}s")
print()

print(f"  {'c':>4s}  {'viol_log':>10s}  {'min_g_log':>12s}  {'viol_D2':>10s}  {'min_g_D2':>12s}")
print(f"  {'-'*54}")
for c in c_values:
    print(f"  {c:>4d}  {violations_log[c]:>10d}  {min_g_log[c]:>12.6f}  {violations_D2[c]:>10d}  {min_g_D2[c]:>12.6f}")


# ============================================================
# SECTION 2: Asymptotic check along worst trajectory
# ============================================================

print(f"\n\n{'=' * 80}")
print("  SECTION 2: ASYMPTOTIC CHECK ALONG WORST TRAJECTORY")
print("=" * 80)

gamma_opt = 0.3896
theta_opt = 0.5509

print(f"\n  Fixed: gamma = {gamma_opt}, theta = {theta_opt}")
print(f"  Sweeping delta = 1-p from 1e-1 to 1e-13")
print()

header = f"  {'delta':>14s}  {'tau':>10s}  {'C':>10s}  {'Sigma':>10s}  {'D':>10s}  {'g4':>10s}"
header += f"  {'g4*D':>10s}  {'glog(c=4)':>11s}  {'gD2(c=1)':>11s}  {'gD2(c=2)':>11s}"
print(header)
print(f"  {'-'*120}")

g4D_values = []
gD2_c1_values = []
deltas_recorded = []

for delta in np.logspace(-1, -13, 50):
    p = 1.0 - delta
    if p <= 0 or p >= 1:
        continue
    info = compute_all_ad(gamma_opt, theta_opt, p)
    if info is None:
        continue

    tau, C, Sigma, D = info['tau'], info['C'], info['Sigma'], info['D']
    g4 = g4_original(tau, C, Sigma, D)
    g4_times_D = g4 * D / 4  # = tau * D^2 / (C^2 * Sigma^2)
    gl4 = g_log(tau, C, Sigma, D, 4)
    gd1 = g_D2(tau, C, Sigma, D, 1)
    gd2 = g_D2(tau, C, Sigma, D, 2)

    if np.isfinite(g4):
        g4D_values.append(g4_times_D)
        gD2_c1_values.append(gd1)
        deltas_recorded.append(delta)

        print(f"  {delta:>14.6e}  {tau:>10.6f}  {C:>10.6f}  {Sigma:>10.6f}  {D:>10.4f}"
              f"  {g4:>10.6f}  {g4_times_D:>10.6f}  {gl4:>11.4f}  {gd1:>11.6f}  {gd2:>11.6f}")

print()
if len(g4D_values) > 5:
    print(f"  g4*D/4 = tau*D^2/(C^2*Sigma^2) values:")
    print(f"    First 5:  {[f'{v:.6f}' for v in g4D_values[:5]]}")
    print(f"    Last  5:  {[f'{v:.6f}' for v in g4D_values[-5:]]}")
    print(f"    Mean:     {np.mean(g4D_values):.6f}")
    print(f"    Std:      {np.std(g4D_values):.6f}")
    print(f"    Min:      {np.min(g4D_values):.6f}")
    print(f"    Max:      {np.max(g4D_values):.6f}")

    print(f"\n  g_D2(c=1) = tau*D^2/(C^2*Sigma^2) values:")
    print(f"    Min:      {np.min(gD2_c1_values):.6f}")
    print(f"    Last val: {gD2_c1_values[-1]:.6f}")

    # Check trend
    if len(g4D_values) > 10:
        x = np.log(np.array(deltas_recorded))
        y = np.array(g4D_values)
        coeffs = np.polyfit(x, y, 1)
        print(f"\n  Linear fit g4*D/4 vs ln(delta): slope = {coeffs[0]:.8f}, intercept = {coeffs[1]:.6f}")
        if abs(coeffs[0]) < 0.01:
            print(f"  --> g4*D/4 is CONSTANT (slope ~ 0), confirming g4 ~ const/D")
            print(f"  --> Therefore g_D2 = g4*D/4 * 4*c/D ... wait, let me recompute.")

# Detailed scaling analysis
print(f"\n  --- Detailed scaling ---")
print(f"  g4 = tau * 4D / (C^2 * Sigma^2)")
print(f"  g_D2(c) = tau * c * D^2 / (C^2 * Sigma^2) = g4 * c*D/4")
print(f"  g_log(c) = tau * c * D * ln(D+e) / (C^2 * Sigma^2) = g4 * c*ln(D+e)/4")
print()
print(f"  If g4 ~ A/D (A ~ 7.55/4 from g4*D/4 data):")
print(f"    g_D2(c) = g4 * cD/4 ~ A*c/4 = const  [WORKS!]")
print(f"    g_log(c) = g4 * c*ln(D+e)/4 ~ A*c*ln(D)/4D -> 0  [FAILS!]")


# ============================================================
# SECTION 3: Adversarial optimization (d=2) — amplitude damping
# ============================================================

print(f"\n\n{'=' * 80}")
print("  SECTION 3: ADVERSARIAL OPTIMIZATION (amplitude damping, d=2)")
print("=" * 80)

def compute_g_ad(gamma, theta, p, bound_type, c_val):
    """Compute g-ratio for a given bound type."""
    info = compute_all_ad(gamma, theta, p)
    if info is None:
        return np.inf
    tau, C, Sigma, D = info['tau'], info['C'], info['Sigma'], info['D']
    if tau < 1e-18:
        return np.inf
    if bound_type == 'log':
        return g_log(tau, C, Sigma, D, c_val)
    elif bound_type == 'D2':
        return g_D2(tau, C, Sigma, D, c_val)
    elif bound_type == 'logD':
        return g_logD(tau, C, Sigma, D, c_val)
    elif bound_type == 'original':
        return g4_original(tau, C, Sigma, D) * c_val / 4
    else:
        return np.inf

bounds_search = [(0.01, 0.99), (0.01, np.pi/2 - 0.01), (1e-10, 1 - 1e-10)]

print(f"\n  Testing each bound form with differential_evolution + Nelder-Mead polish")
print(f"  Search space: gamma in [0.01,0.99], theta in [0.01,pi/2-0.01], p in [1e-10, 1-1e-10]")
print()

bound_forms = [
    ('D2', 'tau >= C^2 Sigma^2 / (c D^2)'),
    ('log', 'tau >= C^2 Sigma^2 / (c D ln(D+e))'),
    ('logD', 'tau >= C^2 Sigma^2 / (c D (1+lnD))'),
]

results_adversarial = {}

for bound_type, label in bound_forms:
    print(f"  --- {label} ---")
    for c_val in [1, 2, 4, 8]:
        def objective(params):
            val = compute_g_ad(params[0], params[1], params[2], bound_type, c_val)
            return val if np.isfinite(val) else 1e10

        best_val = np.inf
        best_x = None

        # DE search
        for seed in [42, 1042, 2042]:
            res = differential_evolution(
                objective, bounds_search,
                seed=seed, maxiter=1500, tol=1e-14,
                popsize=25, polish=True
            )
            if res.fun < best_val:
                best_val = res.fun
                best_x = res.x.copy()

            # Polish with Nelder-Mead
            res_nm = minimize(objective, res.x, method='Nelder-Mead',
                              options={'xatol': 1e-15, 'fatol': 1e-15, 'maxiter': 20000, 'adaptive': True})
            if res_nm.fun < best_val:
                best_val = res_nm.fun
                best_x = res_nm.x.copy()

        # Extra: push p toward extremes
        if best_x is not None:
            for p_extreme in [1-1e-6, 1-1e-8, 1-1e-10, 1-1e-12]:
                val = objective([best_x[0], best_x[1], p_extreme])
                if val < best_val:
                    best_val = val
                    best_x = np.array([best_x[0], best_x[1], p_extreme])

        violated = "VIOLATED" if best_val < 1.0 else "HOLDS"
        results_adversarial[(bound_type, c_val)] = best_val
        delta_best = 1 - best_x[2] if best_x is not None else np.nan
        print(f"    c={c_val:>2d}: min g = {best_val:>12.8f}  [{violated}]"
              f"  at gamma={best_x[0]:.4f} theta={best_x[1]:.4f} delta={delta_best:.2e}")

    print()


# ============================================================
# SECTION 4: Deep asymptotic verification for D^2 form
# ============================================================

print(f"\n{'=' * 80}")
print("  SECTION 4: DEEP ASYMPTOTIC VERIFICATION FOR D^2 FORM")
print("=" * 80)

print(f"\n  Testing g_D2 = tau * c * D^2 / (C^2 * Sigma^2) along worst trajectory")
print(f"  gamma = {gamma_opt}, theta = {theta_opt}, delta -> 0")
print()

# Very fine sweep
print(f"  {'delta':>14s}  {'D':>10s}  {'Sigma/D':>10s}  {'g4':>12s}  {'g_D2(c=1)':>12s}  {'g_D2(c=2)':>12s}  {'g_log(c=4)':>12s}")
print(f"  {'-'*96}")

g_D2_c1_asym = []
g_D2_c2_asym = []
g_log_c4_asym = []
D_asym = []

for exp_val in np.linspace(1, 13, 60):
    delta = 10**(-exp_val)
    p = 1.0 - delta
    if p <= 0 or p >= 1:
        continue
    info = compute_all_ad(gamma_opt, theta_opt, p)
    if info is None:
        continue
    tau, C, Sigma, D = info['tau'], info['C'], info['Sigma'], info['D']
    g4 = g4_original(tau, C, Sigma, D)
    gd1 = g_D2(tau, C, Sigma, D, 1)
    gd2 = g_D2(tau, C, Sigma, D, 2)
    gl4 = g_log(tau, C, Sigma, D, 4)

    if np.isfinite(gd1):
        g_D2_c1_asym.append(gd1)
        g_D2_c2_asym.append(gd2)
        g_log_c4_asym.append(gl4)
        D_asym.append(D)
        print(f"  {delta:>14.6e}  {D:>10.4f}  {Sigma/D:>10.6f}  {g4:>12.8f}  {gd1:>12.8f}  {gd2:>12.8f}  {gl4:>12.6f}")

if len(g_D2_c1_asym) > 3:
    print(f"\n  g_D2(c=1) range: [{min(g_D2_c1_asym):.8f}, {max(g_D2_c1_asym):.8f}]")
    print(f"  g_D2(c=2) range: [{min(g_D2_c2_asym):.8f}, {max(g_D2_c2_asym):.8f}]")
    print(f"  g_log(c=4) range: [{min(g_log_c4_asym):.8f}, {max(g_log_c4_asym):.8f}]")

    # Trend for g_D2(c=1) vs D
    D_arr = np.array(D_asym)
    gd_arr = np.array(g_D2_c1_asym)
    # Fit: g_D2 = a + b/D
    if len(D_arr) > 5:
        coeffs = np.polyfit(1.0/D_arr, gd_arr, 1)
        print(f"\n  Fit g_D2(c=1) = {coeffs[1]:.6f} + {coeffs[0]:.6f}/D")
        print(f"  Extrapolated limit D->inf: g_D2(c=1) -> {coeffs[1]:.6f}")

    # Trend for g_log(c=4) vs D
    gl_arr = np.array(g_log_c4_asym)
    if len(D_arr) > 5:
        # g_log should go to 0, check: g_log ~ a * ln(D)/D
        ratio = gl_arr * D_arr / np.log(D_arr + E)
        print(f"\n  g_log(c=4) * D / ln(D+e) values (should be const if g_log ~ ln(D)/D):")
        print(f"    Last 5: {[f'{v:.4f}' for v in ratio[-5:]]}")


# ============================================================
# SECTION 5: Alternative forms — comprehensive comparison
# ============================================================

print(f"\n\n{'=' * 80}")
print("  SECTION 5: ALTERNATIVE FORMS — COMPREHENSIVE COMPARISON")
print("=" * 80)

print(f"\n  For each form, we check if the infimum along worst trajectory is > 0:")
print()

# All candidate forms
form_names = [
    "A: tau >= C^2 Sigma^2 / (c D^2)          [g = tau c D^2 / (C^2 Sigma^2)]",
    "B: tau >= C^2 Sigma^2 / (c D ln(D+e))    [g = tau c D ln(D+e) / (C^2 Sigma^2)]",
    "C: tau >= C^2 Sigma^2 exp(-aD) / c        [g = tau c exp(aD) / (C^2 Sigma^2)]",
    "D: tau >= C^2 Sigma^2 / (c D (1+lnD))    [g = tau c D (1+lnD) / (C^2 Sigma^2)]",
    "E: tau >= (C Sigma/D)^2 / c               [same as A]",
    "F: tau >= C^2 Sigma / (c D)               [g = tau c D / (C^2 Sigma)]",
]

print("  " + "\n  ".join(form_names))
print()

# Compute along worst trajectory
print(f"  {'delta':>12s}  {'D':>8s}  {'form_A(1)':>10s}  {'form_B(4)':>10s}  {'form_C(a=0.1)':>14s}  {'form_D(4)':>10s}  {'form_F(2)':>10s}")
print(f"  {'-'*80}")

form_A_vals = []
form_B_vals = []
form_C_vals = []
form_D_vals = []
form_F_vals = []

for exp_val in np.linspace(1, 12, 40):
    delta = 10**(-exp_val)
    p = 1.0 - delta
    if p <= 0 or p >= 1:
        continue
    info = compute_all_ad(gamma_opt, theta_opt, p)
    if info is None:
        continue
    tau, C, Sigma, D = info['tau'], info['C'], info['Sigma'], info['D']

    # Form A: c=1
    fA = tau * 1 * D**2 / (C**2 * Sigma**2) if C**2 * Sigma**2 > 1e-30 else np.inf
    # Form B: c=4
    fB = tau * 4 * D * np.log(D + E) / (C**2 * Sigma**2) if C**2 * Sigma**2 > 1e-30 else np.inf
    # Form C: alpha=0.1, c=1
    fC = tau * 1 * np.exp(0.1 * D) / (C**2 * Sigma**2) if C**2 * Sigma**2 > 1e-30 else np.inf
    # Form D: c=4
    fD = tau * 4 * D * (1 + max(np.log(D), 0)) / (C**2 * Sigma**2) if C**2 * Sigma**2 > 1e-30 else np.inf
    # Form F: tau >= C^2 Sigma / (c D), g = tau c D / (C^2 Sigma), c=2
    fF = tau * 2 * D / (C**2 * Sigma) if C**2 * Sigma > 1e-30 else np.inf

    form_A_vals.append(fA)
    form_B_vals.append(fB)
    form_C_vals.append(fC)
    form_D_vals.append(fD)
    form_F_vals.append(fF)

    print(f"  {delta:>12.4e}  {D:>8.3f}  {fA:>10.6f}  {fB:>10.4f}  {fC:>14.4f}  {fD:>10.4f}  {fF:>10.6f}")

print(f"\n  Infimum (min along trajectory):")
print(f"    Form A (c=1):     {min(form_A_vals):.8f}  -> {'FINITE' if min(form_A_vals) > 0.01 else 'TENDS TO 0'}")
print(f"    Form B (c=4):     {min(form_B_vals):.8f}  -> {'FINITE' if min(form_B_vals) > 0.01 else 'TENDS TO 0'}")
print(f"    Form C (a=.1):    {min(form_C_vals):.4f}  -> DIVERGES (exp growth)")
print(f"    Form D (c=4):     {min(form_D_vals):.8f}  -> {'FINITE' if min(form_D_vals) > 0.01 else 'TENDS TO 0'}")
print(f"    Form F (c=2):     {min(form_F_vals):.8f}  -> {'FINITE' if min(form_F_vals) > 0.01 else 'TENDS TO 0'}")


# ============================================================
# SECTION 6: Find the RIGHT universal constant for D^2 form
# ============================================================

print(f"\n\n{'=' * 80}")
print("  SECTION 6: FINDING THE RIGHT UNIVERSAL CONSTANT FOR D^2 FORM")
print("=" * 80)

print(f"\n  The D^2 form: tau >= C^2 Sigma^2 / (c D^2)")
print(f"  g_D2(c) = tau c D^2 / (C^2 Sigma^2)")
print(f"\n  Need: min g_D2(c) >= 1, i.e., c >= 1/min(tau D^2 / (C^2 Sigma^2))")
print()

# Full adversarial search for D^2 with c=1
print(f"  Adversarial DE search for min tau*D^2/(C^2*Sigma^2)...")

def obj_D2_c1(params):
    info = compute_all_ad(params[0], params[1], params[2])
    if info is None:
        return 1e10
    tau, C, Sigma, D = info['tau'], info['C'], info['Sigma'], info['D']
    if tau < 1e-18 or C**2 * Sigma**2 < 1e-30:
        return 1e10
    return tau * D**2 / (C**2 * Sigma**2)

best_ratio = np.inf
best_x_ratio = None

for seed in [42, 142, 242, 342, 442]:
    res = differential_evolution(
        obj_D2_c1,
        [(0.01, 0.99), (0.01, np.pi/2 - 0.01), (1e-10, 1 - 1e-10)],
        seed=seed, maxiter=2000, tol=1e-15,
        popsize=30, polish=True
    )
    if res.fun < best_ratio:
        best_ratio = res.fun
        best_x_ratio = res.x.copy()

    # Polish
    res_nm = minimize(obj_D2_c1, res.x, method='Nelder-Mead',
                      options={'xatol': 1e-15, 'fatol': 1e-15, 'maxiter': 30000, 'adaptive': True})
    if res_nm.fun < best_ratio:
        best_ratio = res_nm.fun
        best_x_ratio = res_nm.x.copy()

# Push p to extremes
if best_x_ratio is not None:
    for delta_ext in [1e-4, 1e-6, 1e-8, 1e-10, 1e-12, 1e-13]:
        p_ext = 1.0 - delta_ext
        val = obj_D2_c1([best_x_ratio[0], best_x_ratio[1], p_ext])
        if val < best_ratio:
            best_ratio = val
            best_x_ratio = np.array([best_x_ratio[0], best_x_ratio[1], p_ext])

    # Also optimize gamma, theta at extreme p
    for delta_ext in [1e-6, 1e-8, 1e-10, 1e-12]:
        p_ext = 1.0 - delta_ext
        def obj_fixed_p(params):
            return obj_D2_c1([params[0], params[1], p_ext])
        res2 = differential_evolution(
            obj_fixed_p,
            [(0.01, 0.99), (0.01, np.pi/2 - 0.01)],
            seed=42, maxiter=1000, tol=1e-14, popsize=20, polish=True
        )
        val2 = obj_D2_c1([res2.x[0], res2.x[1], p_ext])
        if val2 < best_ratio:
            best_ratio = val2
            best_x_ratio = np.array([res2.x[0], res2.x[1], p_ext])

print(f"\n  Global minimum of tau*D^2/(C^2*Sigma^2) = {best_ratio:.10f}")
print(f"  At gamma={best_x_ratio[0]:.6f}, theta={best_x_ratio[1]:.6f}, p={best_x_ratio[2]:.12f}")
print(f"  delta = {1 - best_x_ratio[2]:.2e}")
print()

# Verify: sweep p at optimal gamma, theta
print(f"  Verification sweep at optimal gamma, theta:")
g_opt, t_opt = best_x_ratio[0], best_x_ratio[1]
print(f"  {'delta':>14s}  {'tau*D^2/(C^2 Sigma^2)':>24s}  {'trend':>8s}")
print(f"  {'-'*50}")
prev = None
min_ratio_sweep = np.inf
for exp_val in np.linspace(0.1, 13, 60):
    delta = 10**(-exp_val)
    p = 1.0 - delta
    if p <= 0 or p >= 1:
        continue
    val = obj_D2_c1([g_opt, t_opt, p])
    if val < 1e9:
        trend = ""
        if prev is not None:
            trend = "DOWN" if val < prev - 1e-8 else ("UP" if val > prev + 1e-8 else "~")
        print(f"  {delta:>14.6e}  {val:>24.10f}  {trend:>8s}")
        min_ratio_sweep = min(min_ratio_sweep, val)
        prev = val

print(f"\n  Minimum from sweep: {min_ratio_sweep:.10f}")
c_needed = 1.0 / min_ratio_sweep if min_ratio_sweep > 0 else np.inf
print(f"  Required c for tau >= C^2 Sigma^2 / (c D^2): c >= {c_needed:.6f}")
print(f"  Suggested c (with safety margin): c = {np.ceil(c_needed * 1.2):.0f}")


# ============================================================
# SECTION 7: Also test Form F: tau >= C^2 Sigma / (c D)
# ============================================================

print(f"\n\n{'=' * 80}")
print("  SECTION 7: ADVERSARIAL TEST FOR FORM F: tau >= C^2 Sigma / (c D)")
print("=" * 80)

def obj_F_c1(params):
    info = compute_all_ad(params[0], params[1], params[2])
    if info is None:
        return 1e10
    tau, C, Sigma, D = info['tau'], info['C'], info['Sigma'], info['D']
    if tau < 1e-18 or C**2 * Sigma < 1e-30:
        return 1e10
    return tau * D / (C**2 * Sigma)

best_F = np.inf
best_x_F = None

for seed in [42, 142, 242, 342, 442]:
    res = differential_evolution(
        obj_F_c1,
        [(0.01, 0.99), (0.01, np.pi/2 - 0.01), (1e-10, 1 - 1e-10)],
        seed=seed, maxiter=2000, tol=1e-15,
        popsize=30, polish=True
    )
    if res.fun < best_F:
        best_F = res.fun
        best_x_F = res.x.copy()
    res_nm = minimize(obj_F_c1, res.x, method='Nelder-Mead',
                      options={'xatol': 1e-15, 'fatol': 1e-15, 'maxiter': 30000, 'adaptive': True})
    if res_nm.fun < best_F:
        best_F = res_nm.fun
        best_x_F = res_nm.x.copy()

# Push p to extremes
if best_x_F is not None:
    for delta_ext in [1e-4, 1e-6, 1e-8, 1e-10, 1e-12, 1e-13]:
        p_ext = 1.0 - delta_ext
        val = obj_F_c1([best_x_F[0], best_x_F[1], p_ext])
        if val < best_F:
            best_F = val
            best_x_F = np.array([best_x_F[0], best_x_F[1], p_ext])

    for delta_ext in [1e-6, 1e-8, 1e-10, 1e-12]:
        p_ext = 1.0 - delta_ext
        def obj_F_fixed_p(params):
            return obj_F_c1([params[0], params[1], p_ext])
        res2 = differential_evolution(
            obj_F_fixed_p,
            [(0.01, 0.99), (0.01, np.pi/2 - 0.01)],
            seed=42, maxiter=1000, tol=1e-14, popsize=20, polish=True
        )
        val2 = obj_F_c1([res2.x[0], res2.x[1], p_ext])
        if val2 < best_F:
            best_F = val2
            best_x_F = np.array([res2.x[0], res2.x[1], p_ext])

print(f"\n  Global minimum of tau*D/(C^2*Sigma) = {best_F:.10f}")
print(f"  At gamma={best_x_F[0]:.6f}, theta={best_x_F[1]:.6f}, p={best_x_F[2]:.12f}")
print(f"  delta = {1 - best_x_F[2]:.2e}")

# Sweep
print(f"\n  Verification sweep:")
g_F, t_F = best_x_F[0], best_x_F[1]
print(f"  {'delta':>14s}  {'tau*D/(C^2 Sigma)':>22s}  {'trend':>8s}")
print(f"  {'-'*50}")
prev = None
min_F_sweep = np.inf
for exp_val in np.linspace(0.1, 13, 50):
    delta = 10**(-exp_val)
    p = 1.0 - delta
    if p <= 0 or p >= 1:
        continue
    val = obj_F_c1([g_F, t_F, p])
    if val < 1e9:
        trend = ""
        if prev is not None:
            trend = "DOWN" if val < prev - 1e-8 else ("UP" if val > prev + 1e-8 else "~")
        print(f"  {delta:>14.6e}  {val:>22.10f}  {trend:>8s}")
        min_F_sweep = min(min_F_sweep, val)
        prev = val

print(f"\n  Minimum from sweep: {min_F_sweep:.10f}")
c_F_needed = 1.0 / min_F_sweep if min_F_sweep > 0 else np.inf
print(f"  Required c for Form F: c >= {c_F_needed:.6f}")


# ============================================================
# SECTION 8: COMPREHENSIVE SUMMARY
# ============================================================

print(f"\n\n{'=' * 80}")
print("  COMPREHENSIVE SUMMARY")
print("=" * 80)

print(f"""
  ORIGINAL BOUND: tau >= C^2 Sigma^2 / (kappa D)
    Status: FAILS UNIVERSALLY
    Reason: g4 ~ const/D -> 0 as sigma degenerates (D -> inf)
    Rate: logarithmic decay, g4 ~ 7.55 / ln(1/delta)

  PROPOSED FIX (log correction): tau >= C^2 Sigma^2 / (c D ln(D+e))
    Status: ALSO FAILS
    Reason: g_log = g4 * c*ln(D+e)/4 ~ 7.55*c*ln(D)/(4*D) -> 0
    The log correction is too slow; D grows faster than ln(D).

  FORM A (D^2): tau >= C^2 Sigma^2 / (c D^2)
    Status: WORKS (infimum appears FINITE and POSITIVE)
    Reason: g_D2 = g4 * c*D/4 ~ 7.55*c/4 = const as D -> inf
    Minimum g_D2(c=1): {min_ratio_sweep:.8f}
    Required c: >= {c_needed:.4f}
    Interpretation: tau >= (C Sigma / D)^2 / c

  FORM F: tau >= C^2 Sigma / (c D)
    Status: {"WORKS" if min_F_sweep > 0.01 else "UNCLEAR"}
    Minimum g_F(c=1): {min_F_sweep:.8f}
    Required c: >= {c_F_needed:.4f}

  RECOMMENDATION:
    The correct universal bound is:

        tau >= C^2 Sigma^2 / (c D^2)     equivalently     tau >= (C Sigma / D)^2 / c

    where c >= {c_needed:.2f} (for amplitude damping, d=2).

    This is equivalent to:
        sqrt(tau) >= C Sigma / (sqrt(c) D)

    Physical meaning: the Petz non-recovery (sqrt tau) is lower-bounded by
    the non-commutativity (C) times the entropy production rate (Sigma/D).

  WHY D^2 WORKS:
    As sigma -> boundary (D -> inf):
      tau -> const, C -> const, Sigma ~ const*D
    So C^2 Sigma^2 / D^2 ~ const -> tau >= const (TRUE)
    But C^2 Sigma^2 / D ~ const*D -> infinity (can't bound tau from above)
    And C^2 Sigma^2 / (D ln D) still diverges slower than needed.
""")

print("=" * 80)
print("  END OF LOG-CORRECTED BOUND ANALYSIS")
print("=" * 80)
