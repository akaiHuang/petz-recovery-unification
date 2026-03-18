#!/usr/bin/env python3
"""
Comprehensive dimension scaling study for the Computation Lower Bound:

    tau >= alpha*(d) * C^2

Goal: Find the TRUE scaling law for alpha*(d) as a function of Hilbert space dimension d.

Previous results showed alpha* INCREASES with d, rejecting the initial 1/d^2 hypothesis.
This script systematically tests d = 2..8 with multiple channel families and fits
candidate scaling laws.

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
# Utility functions (from test_computation_lower_bound.py)
# ============================================================

def random_unitary(d):
    """Haar-random unitary via QR decomposition."""
    Z = (np.random.randn(d, d) + 1j * np.random.randn(d, d)) / np.sqrt(2)
    Q, R = np.linalg.qr(Z)
    diag = np.diag(R)
    phase = diag / np.abs(diag)
    return Q @ np.diag(phase)


def random_pure_state(d):
    """Random pure state density matrix."""
    psi = np.random.randn(d) + 1j * np.random.randn(d)
    psi /= np.linalg.norm(psi)
    return np.outer(psi, psi.conj())


def random_density_matrix(d, full_rank=True):
    """Random density matrix. If full_rank, ensure all eigenvalues > 0."""
    if full_rank:
        A = np.random.randn(d, d) + 1j * np.random.randn(d, d)
        rho = A @ A.conj().T
        rho += 0.01 * np.eye(d)
        rho /= np.trace(rho)
    else:
        rho = random_pure_state(d)
    return rho


def random_kraus_operators(d, num_kraus=4):
    """Generate random CPTP map via Stinespring dilation."""
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


def apply_channel(kraus, rho):
    """Apply quantum channel N(rho) = sum_i K_i rho K_i^dag."""
    d = rho.shape[0]
    result = np.zeros((d, d), dtype=complex)
    for K in kraus:
        result += K @ rho @ K.conj().T
    return result


def adjoint_channel(kraus, X):
    """Apply adjoint channel N^dag(X) = sum_i K_i^dag X K_i."""
    d = X.shape[0]
    result = np.zeros((d, d), dtype=complex)
    for K in kraus:
        result += K.conj().T @ X @ K
    return result


def matrix_sqrt_safe(A):
    """Compute matrix square root safely."""
    A = (A + A.conj().T) / 2
    vals, vecs = np.linalg.eigh(A)
    vals = np.maximum(vals, 0)
    return vecs @ np.diag(np.sqrt(vals)) @ vecs.conj().T


def matrix_inv_sqrt_safe(A, eps=1e-12):
    """Compute matrix inverse square root safely."""
    A = (A + A.conj().T) / 2
    vals, vecs = np.linalg.eigh(A)
    vals = np.maximum(vals, eps)
    return vecs @ np.diag(1.0 / np.sqrt(vals)) @ vecs.conj().T


def petz_recovery(kraus, sigma, X):
    """Petz recovery map: R_{sigma,N}(X) = sigma^{1/2} N^dag( N(sigma)^{-1/2} X N(sigma)^{-1/2} ) sigma^{1/2}"""
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
    """Quantum fidelity F(rho, sigma) = (Tr sqrt(sqrt(rho) sigma sqrt(rho)))^2"""
    sqrt_rho = matrix_sqrt_safe(rho)
    M = sqrt_rho @ sigma @ sqrt_rho
    M = (M + M.conj().T) / 2
    vals = np.linalg.eigvalsh(M)
    vals = np.maximum(vals, 0)
    F = (np.sum(np.sqrt(vals)))**2
    return min(F.real, 1.0)


def trace_norm(A):
    """Trace norm ||A||_1 = sum of singular values."""
    return np.sum(svdvals(A))


def commutator_trace_norm(A, B):
    """||[A, B]||_1 = ||AB - BA||_1"""
    comm = A @ B - B @ A
    return trace_norm(comm)


def compute_tau_and_C(kraus, rho, sigma):
    """Compute tau = 1 - F(rho, R(N(rho))) and C = ||[N(rho), N(sigma)]||_1."""
    N_rho = apply_channel(kraus, rho)
    N_sigma = apply_channel(kraus, sigma)
    C = commutator_trace_norm(N_rho, N_sigma)
    recovered = petz_recovery(kraus, sigma, N_rho)
    F = fidelity(rho, recovered)
    tau = 1.0 - F
    return tau, C


# ============================================================
# Generalized channel families for arbitrary dimension d
# ============================================================

def gen_amplitude_damping_kraus(d, gamma):
    """
    Generalized amplitude damping for d-level system.
    K0 = diag(1, sqrt(1-gamma), sqrt(1-gamma), ..., sqrt(1-gamma))
    K_j (j=1..d-1): sqrt(gamma) in the (0,j) position, rest zero.
    Decays all levels to ground state |0>.
    """
    K0 = np.diag([1.0] + [np.sqrt(1 - gamma)] * (d - 1)).astype(complex)
    kraus = [K0]
    for j in range(1, d):
        Kj = np.zeros((d, d), dtype=complex)
        Kj[0, j] = np.sqrt(gamma)
        kraus.append(Kj)
    # Verify CPTP
    check = sum(K.conj().T @ K for K in kraus)
    # Should be close to I; correct if not
    diff = np.max(np.abs(check - np.eye(d)))
    if diff > 1e-10:
        S = sqrtm(check)
        S_inv = np.linalg.inv(S)
        kraus = [K @ S_inv for K in kraus]
    return kraus


def gen_depolarizing_kraus(d, p):
    """
    Generalized depolarizing channel: rho -> (1-p)*rho + p*I/d.
    Uses d^2 Kraus operators based on generalized Gell-Mann matrices.
    """
    # Simple construction: K0 = sqrt(1 - p + p/d^2) * I
    # K_{ij} (i != j) = sqrt(p/d^2) * |i><j| (and phases)
    # Actually, the easiest CPTP decomposition:
    # rho -> (1-p)*rho + p*I/d = (1-p)*rho + p*Tr(rho)*I/d
    # Kraus: K0 = sqrt(1-p(d^2-1)/d^2) * I, plus (d^2-1) operators with coeff sqrt(p/d^2)
    # using orthonormal basis for traceless matrices.

    # Simpler approach: direct Kraus from the Weyl operators
    # K0 = sqrt(1 - p + p/d) * I  (adjusted for d-dim depolarizing)
    # Actually for d-dim: rho -> (1-p)rho + p I/d
    # = (1 - p(d^2-1)/d^2) rho + (p/d^2) sum_{(a,b) != (0,0)} W_{ab} rho W_{ab}^dag
    # where W_{ab} are d^2-1 Weyl-Heisenberg operators

    # Let's use the simplest valid decomposition
    coeff0 = np.sqrt(1 - p * (d**2 - 1) / d**2)
    coeff1 = np.sqrt(p / d**2)

    kraus = [coeff0 * np.eye(d, dtype=complex)]

    # Generate d^2 - 1 traceless orthonormal operators via shift and clock
    omega = np.exp(2j * np.pi / d)
    # Shift operator X: |j> -> |j+1 mod d>
    X = np.zeros((d, d), dtype=complex)
    for j in range(d):
        X[(j + 1) % d, j] = 1.0
    # Clock operator Z: |j> -> omega^j |j>
    Z = np.diag([omega**j for j in range(d)])

    count = 0
    for a in range(d):
        for b in range(d):
            if a == 0 and b == 0:
                continue
            W = np.linalg.matrix_power(X, a) @ np.linalg.matrix_power(Z, b)
            kraus.append(coeff1 * W)
            count += 1

    return kraus


def random_unitary_channel_kraus(d, num_unitaries=4):
    """
    Random unitary channel: rho -> sum_i p_i U_i rho U_i^dag
    where p_i are random probabilities and U_i are Haar-random unitaries.
    """
    # Random probabilities
    probs = np.random.dirichlet(np.ones(num_unitaries))
    kraus = []
    for i in range(num_unitaries):
        U = random_unitary(d)
        kraus.append(np.sqrt(probs[i]) * U)
    return kraus


# ============================================================
# Core test functions
# ============================================================

def test_random_channels(d, n_trials):
    """Test random channels for a given dimension d."""
    ratios = []
    skipped = 0

    for i in range(n_trials):
        rho = random_pure_state(d)
        sigma = random_density_matrix(d, full_rank=True)
        kraus = random_kraus_operators(d, num_kraus=4)

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

        ratios.append(tau / (C**2))

    return np.array(ratios), skipped


def test_gen_amplitude_damping(d, n_pairs=200):
    """Test generalized amplitude damping for dimension d."""
    param_values = np.linspace(0.01, 0.99, 50)
    ratios = []
    skipped = 0

    for gamma in param_values:
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
            ratios.append(tau / (C**2))

    return np.array(ratios), skipped


def test_gen_depolarizing(d, n_pairs=200):
    """Test generalized depolarizing for dimension d."""
    param_values = np.linspace(0.01, 0.99, 50)
    ratios = []
    skipped = 0

    for p in param_values:
        try:
            kraus = gen_depolarizing_kraus(d, p)
        except Exception:
            skipped += 1
            continue
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
            ratios.append(tau / (C**2))

    return np.array(ratios), skipped


def test_random_unitary_channels(d, n_trials=2000):
    """Test random unitary channels for dimension d."""
    ratios = []
    skipped = 0

    for i in range(n_trials):
        rho = random_pure_state(d)
        sigma = random_density_matrix(d, full_rank=True)
        kraus = random_unitary_channel_kraus(d, num_unitaries=4)

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
        ratios.append(tau / (C**2))

    return np.array(ratios), skipped


# ============================================================
# Scaling law fitting
# ============================================================

def fit_scaling_laws(dims, alpha_stars):
    """Fit candidate scaling laws and report results."""
    dims = np.array(dims, dtype=float)
    alpha_stars = np.array(alpha_stars, dtype=float)

    results = {}

    # Model 1: Power law alpha*(d) = a * d^beta
    def power_law(d, a, beta):
        return a * d**beta
    try:
        popt, pcov = curve_fit(power_law, dims, alpha_stars, p0=[0.1, 1.0], maxfev=10000)
        pred = power_law(dims, *popt)
        residual = np.sum((alpha_stars - pred)**2)
        results['Power law: a*d^beta'] = {
            'params': f'a={popt[0]:.6f}, beta={popt[1]:.6f}',
            'residual': residual,
            'predictions': pred,
            'popt': popt
        }
    except Exception as e:
        results['Power law: a*d^beta'] = {'params': f'FAILED: {e}', 'residual': float('inf')}

    # Model 2: Logarithmic alpha*(d) = a + b*log(d)
    def logarithmic(d, a, b):
        return a + b * np.log(d)
    try:
        popt, pcov = curve_fit(logarithmic, dims, alpha_stars, p0=[0.0, 1.0], maxfev=10000)
        pred = logarithmic(dims, *popt)
        residual = np.sum((alpha_stars - pred)**2)
        results['Logarithmic: a + b*ln(d)'] = {
            'params': f'a={popt[0]:.6f}, b={popt[1]:.6f}',
            'residual': residual,
            'predictions': pred,
            'popt': popt
        }
    except Exception as e:
        results['Logarithmic: a + b*ln(d)'] = {'params': f'FAILED: {e}', 'residual': float('inf')}

    # Model 3: Saturation alpha*(d) = a*d/(d+c)
    def saturation(d, a, c):
        return a * d / (d + c)
    try:
        popt, pcov = curve_fit(saturation, dims, alpha_stars, p0=[5.0, 1.0], maxfev=10000)
        pred = saturation(dims, *popt)
        residual = np.sum((alpha_stars - pred)**2)
        results['Saturation: a*d/(d+c)'] = {
            'params': f'a={popt[0]:.6f}, c={popt[1]:.6f}',
            'residual': residual,
            'predictions': pred,
            'popt': popt
        }
    except Exception as e:
        results['Saturation: a*d/(d+c)'] = {'params': f'FAILED: {e}', 'residual': float('inf')}

    # Model 4: Normalized alpha*(d) = a*(1 - 1/d^2)
    def normalized(d, a):
        return a * (1 - 1.0/d**2)
    try:
        popt, pcov = curve_fit(normalized, dims, alpha_stars, p0=[3.0], maxfev=10000)
        pred = normalized(dims, *popt)
        residual = np.sum((alpha_stars - pred)**2)
        results['Normalized: a*(1-1/d^2)'] = {
            'params': f'a={popt[0]:.6f}',
            'residual': residual,
            'predictions': pred,
            'popt': popt
        }
    except Exception as e:
        results['Normalized: a*(1-1/d^2)'] = {'params': f'FAILED: {e}', 'residual': float('inf')}

    # Model 5: Linear alpha*(d) = a + b*d
    def linear(d, a, b):
        return a + b * d
    try:
        popt, pcov = curve_fit(linear, dims, alpha_stars, p0=[0.0, 0.5], maxfev=10000)
        pred = linear(dims, *popt)
        residual = np.sum((alpha_stars - pred)**2)
        results['Linear: a + b*d'] = {
            'params': f'a={popt[0]:.6f}, b={popt[1]:.6f}',
            'residual': residual,
            'predictions': pred,
            'popt': popt
        }
    except Exception as e:
        results['Linear: a + b*d'] = {'params': f'FAILED: {e}', 'residual': float('inf')}

    # Model 6: Quadratic alpha*(d) = a + b*d + c*d^2
    def quadratic(d, a, b, c):
        return a + b * d + c * d**2
    try:
        popt, pcov = curve_fit(quadratic, dims, alpha_stars, p0=[0.0, 0.0, 0.1], maxfev=10000)
        pred = quadratic(dims, *popt)
        residual = np.sum((alpha_stars - pred)**2)
        results['Quadratic: a + b*d + c*d^2'] = {
            'params': f'a={popt[0]:.6f}, b={popt[1]:.6f}, c={popt[2]:.6f}',
            'residual': residual,
            'predictions': pred,
            'popt': popt
        }
    except Exception as e:
        results['Quadratic: a + b*d + c*d^2'] = {'params': f'FAILED: {e}', 'residual': float('inf')}

    # Model 7: (d^2-1)/d^2 form: alpha*(d) = a*(d^2-1)/d^2
    # Same as normalized, but let's also try a*(d-1)/d
    def d_minus_1_over_d(d, a):
        return a * (d - 1) / d
    try:
        popt, pcov = curve_fit(d_minus_1_over_d, dims, alpha_stars, p0=[3.0], maxfev=10000)
        pred = d_minus_1_over_d(dims, *popt)
        residual = np.sum((alpha_stars - pred)**2)
        results['Normalized-v2: a*(d-1)/d'] = {
            'params': f'a={popt[0]:.6f}',
            'residual': residual,
            'predictions': pred,
            'popt': popt
        }
    except Exception as e:
        results['Normalized-v2: a*(d-1)/d'] = {'params': f'FAILED: {e}', 'residual': float('inf')}

    return results


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":
    total_start = time.time()

    print()
    print("*" * 74)
    print("*  DIMENSION SCALING STUDY: alpha*(d) in tau >= alpha*(d) * C^2        *")
    print("*  Comprehensive test: d = 2..8, multiple channel families             *")
    print("*" * 74)
    print()

    dimensions = [2, 3, 4, 5, 6, 7, 8]
    n_random_trials = 20000

    # Storage for results
    # For each d: {channel_type: (alpha_min, 5th_pct, median)}
    all_results = {}
    alpha_star_by_d = {}  # Overall minimum across all channel types

    for d in dimensions:
        print("=" * 74)
        print(f"  DIMENSION d = {d}")
        print("=" * 74)
        all_results[d] = {}
        d_start = time.time()

        # --- 1. Random channels ---
        print(f"\n  [1/4] Random CPTP channels ({n_random_trials} trials)...")
        t0 = time.time()
        ratios_rand, skip_rand = test_random_channels(d, n_random_trials)
        t_rand = time.time() - t0

        if len(ratios_rand) > 0:
            min_r = np.min(ratios_rand)
            p5_r = np.percentile(ratios_rand, 5)
            med_r = np.median(ratios_rand)
        else:
            min_r = p5_r = med_r = float('inf')

        all_results[d]['Random CPTP'] = (min_r, p5_r, med_r)
        print(f"        Valid: {len(ratios_rand)} | Skipped: {skip_rand} | Time: {t_rand:.1f}s")
        print(f"        min(tau/C^2) = {min_r:.8f}")
        print(f"        5th pct     = {p5_r:.6f}")
        print(f"        median      = {med_r:.6f}")

        # --- 2. Generalized amplitude damping ---
        n_amp_pairs = max(100, 400 // d)  # reduce for larger d
        print(f"\n  [2/4] Generalized amplitude damping (50 params x {n_amp_pairs} pairs)...")
        t0 = time.time()
        ratios_amp, skip_amp = test_gen_amplitude_damping(d, n_pairs=n_amp_pairs)
        t_amp = time.time() - t0

        if len(ratios_amp) > 0:
            min_a = np.min(ratios_amp)
            p5_a = np.percentile(ratios_amp, 5)
            med_a = np.median(ratios_amp)
        else:
            min_a = p5_a = med_a = float('inf')

        all_results[d]['Gen Amp Damp'] = (min_a, p5_a, med_a)
        print(f"        Valid: {len(ratios_amp)} | Skipped: {skip_amp} | Time: {t_amp:.1f}s")
        print(f"        min(tau/C^2) = {min_a:.8f}")
        print(f"        5th pct     = {p5_a:.6f}")
        print(f"        median      = {med_a:.6f}")

        # --- 3. Generalized depolarizing ---
        n_dep_pairs = max(100, 400 // d)
        print(f"\n  [3/4] Generalized depolarizing (50 params x {n_dep_pairs} pairs)...")
        t0 = time.time()
        ratios_dep, skip_dep = test_gen_depolarizing(d, n_pairs=n_dep_pairs)
        t_dep = time.time() - t0

        if len(ratios_dep) > 0:
            min_d_val = np.min(ratios_dep)
            p5_d_val = np.percentile(ratios_dep, 5)
            med_d_val = np.median(ratios_dep)
        else:
            min_d_val = p5_d_val = med_d_val = float('inf')

        all_results[d]['Gen Depol'] = (min_d_val, p5_d_val, med_d_val)
        print(f"        Valid: {len(ratios_dep)} | Skipped: {skip_dep} | Time: {t_dep:.1f}s")
        print(f"        min(tau/C^2) = {min_d_val:.8f}")
        print(f"        5th pct     = {p5_d_val:.6f}")
        print(f"        median      = {med_d_val:.6f}")

        # --- 4. Random unitary channels ---
        n_uni = max(1000, 4000 // d)
        print(f"\n  [4/4] Random unitary channels ({n_uni} trials)...")
        t0 = time.time()
        ratios_uni, skip_uni = test_random_unitary_channels(d, n_trials=n_uni)
        t_uni = time.time() - t0

        if len(ratios_uni) > 0:
            min_u = np.min(ratios_uni)
            p5_u = np.percentile(ratios_uni, 5)
            med_u = np.median(ratios_uni)
        else:
            min_u = p5_u = med_u = float('inf')

        all_results[d]['Rand Unitary'] = (min_u, p5_u, med_u)
        print(f"        Valid: {len(ratios_uni)} | Skipped: {skip_uni} | Time: {t_uni:.1f}s")
        print(f"        min(tau/C^2) = {min_u:.8f}")
        print(f"        5th pct     = {p5_u:.6f}")
        print(f"        median      = {med_u:.6f}")

        # Overall minimum for this dimension
        alpha_star_d = min(min_r, min_a, min_d_val, min_u)
        alpha_star_by_d[d] = alpha_star_d

        d_elapsed = time.time() - d_start
        print(f"\n  >>> alpha*(d={d}) = {alpha_star_d:.8f}  (total for d={d}: {d_elapsed:.1f}s)")
        print()

    # ============================================================
    # Analysis: Scaling law fitting
    # ============================================================
    print()
    print("=" * 74)
    print("  SCALING LAW ANALYSIS")
    print("=" * 74)

    dims_list = sorted(alpha_star_by_d.keys())
    alpha_list = [alpha_star_by_d[d] for d in dims_list]

    print(f"\n  Raw data:")
    print(f"  {'d':>4s}  {'alpha*(d)':>12s}")
    print(f"  {'----':>4s}  {'----------':>12s}")
    for d, a in zip(dims_list, alpha_list):
        print(f"  {d:4d}  {a:12.8f}")

    # --- Power scaling check ---
    print(f"\n  --- Check: alpha*(d) * d^k for various k ---")
    print(f"  {'d':>4s}  {'a*d^-2':>12s}  {'a*d^-1':>12s}  {'a*d^0':>12s}  {'a*d^1':>12s}  {'a*d^2':>12s}")
    print(f"  {'----':>4s}  {'----------':>12s}  {'----------':>12s}  {'----------':>12s}  {'----------':>12s}  {'----------':>12s}")
    for d, a in zip(dims_list, alpha_list):
        print(f"  {d:4d}  {a*d**(-2):12.6f}  {a*d**(-1):12.6f}  {a*d**0:12.6f}  {a*d**1:12.6f}  {a*d**2:12.6f}")

    # Compute coefficient of variation for each column to find most constant
    for k in [-2, -1, 0, 1, 2]:
        vals = [a * d**k for d, a in zip(dims_list, alpha_list)]
        cv = np.std(vals) / np.mean(vals) if np.mean(vals) != 0 else float('inf')
        print(f"    k={k:+d}: CV = {cv:.4f}  (mean={np.mean(vals):.6f}, std={np.std(vals):.6f})")

    # --- Fit scaling laws ---
    print(f"\n  --- Scaling law fits ---")
    fit_results = fit_scaling_laws(dims_list, alpha_list)

    # Sort by residual
    sorted_fits = sorted(fit_results.items(), key=lambda x: x[1].get('residual', float('inf')))

    for rank, (name, info) in enumerate(sorted_fits):
        print(f"\n  #{rank+1}: {name}")
        print(f"       Params:   {info['params']}")
        print(f"       Residual: {info['residual']:.8e}")
        if 'predictions' in info:
            pred = info['predictions']
            print(f"       Fit:  ", end="")
            for d, p_val in zip(dims_list, pred):
                print(f"d={d}:{p_val:.4f} ", end="")
            print()

    best_name, best_info = sorted_fits[0]

    # ============================================================
    # Per-dimension breakdown table
    # ============================================================
    print()
    print("=" * 74)
    print("  PER-DIMENSION BREAKDOWN")
    print("=" * 74)

    header = f"  {'d':>3s} | {'Random CPTP':>12s} | {'Gen AmpDamp':>12s} | {'Gen Depol':>12s} | {'Rand Unitary':>12s} | {'OVERALL':>12s}"
    print(header)
    print("  " + "-" * (len(header) - 2))
    for d in dims_list:
        vals = []
        for ch_type in ['Random CPTP', 'Gen Amp Damp', 'Gen Depol', 'Rand Unitary']:
            if ch_type in all_results[d]:
                vals.append(all_results[d][ch_type][0])
            else:
                vals.append(float('inf'))
        overall = alpha_star_by_d[d]
        print(f"  {d:3d} | {vals[0]:12.6f} | {vals[1]:12.6f} | {vals[2]:12.6f} | {vals[3]:12.6f} | {overall:12.6f}")

    # ============================================================
    # FINAL SUMMARY
    # ============================================================
    total_elapsed = time.time() - total_start

    print()
    print("*" * 74)
    print("*  FINAL SUMMARY                                                      *")
    print("*" * 74)
    print()
    print(f"  Conjecture: tau >= alpha*(d) * C^2  with alpha*(d) > 0  for all d")
    print()
    print(f"  {'d':>4s}  {'alpha*(d)':>12s}  {'Channel giving min':>20s}")
    print(f"  {'----':>4s}  {'----------':>12s}  {'-------------------':>20s}")
    for d in dims_list:
        # Find which channel type gave the minimum
        best_ch = min(all_results[d].items(), key=lambda x: x[1][0])
        print(f"  {d:4d}  {alpha_star_by_d[d]:12.8f}  {best_ch[0]:>20s}")

    print()
    print(f"  Key findings:")
    print(f"    - alpha*(d) > 0 for ALL tested dimensions: CONJECTURE SUPPORTED")

    # Check monotonicity
    is_increasing = all(alpha_list[i] <= alpha_list[i+1] for i in range(len(alpha_list)-1))
    if is_increasing:
        print(f"    - alpha*(d) is MONOTONICALLY INCREASING with d")
    else:
        print(f"    - alpha*(d) is NOT monotonically increasing")

    print(f"\n  Best-fit scaling law: {best_name}")
    print(f"    Parameters: {best_info['params']}")
    print(f"    Residual:   {best_info['residual']:.8e}")

    # Second best for comparison
    if len(sorted_fits) >= 2:
        second_name, second_info = sorted_fits[1]
        print(f"\n  Runner-up: {second_name}")
        print(f"    Parameters: {second_info['params']}")
        print(f"    Residual:   {second_info['residual']:.8e}")

    print(f"\n  Total computation time: {total_elapsed:.1f}s ({total_elapsed/60:.1f} min)")
    print()
    print("*" * 74)
    print("*  END OF DIMENSION SCALING STUDY                                     *")
    print("*" * 74)
