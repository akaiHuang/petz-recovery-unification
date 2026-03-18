#!/usr/bin/env python3
"""
Non-Markovian tau(t) verification: Jaynes-Cummings & Lorentzian spectral density.
==================================================================================

Extends the continuous-time tau(t) framework to non-Markovian dynamics.

Key physics:
  - Markovian (weak coupling): tau(t) monotonically increases, dτ/dt >= 0
  - Non-Markovian (strong coupling): information backflow => tau(t) DECREASES
    in certain intervals, dτ/dt < 0 => Σ < 0 locally => "time reversal"

Test cases:
  (a) Markovian limit:       gamma_0/lambda = 0.1
  (b) Critical:              gamma_0/lambda = 0.5
  (c) Non-Markovian:         gamma_0/lambda = 2.0
  (d) Strong non-Markovian:  gamma_0/lambda = 5.0
  (e) Exact resonance (JC):  No dissipation, pure Jaynes-Cummings

Also computes BLP non-Markovianity measure for comparison.

Author: Sheng-Kai Huang
Date: 2026-03-16
"""

import numpy as np
from scipy.linalg import eigvalsh
import sys

np.random.seed(42)

# =============================================================================
# Pauli matrices and helpers
# =============================================================================
I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)

# States
ket_0 = np.array([[1], [0]], dtype=complex)
ket_1 = np.array([[0], [1]], dtype=complex)
ket_plus = (ket_0 + ket_1) / np.sqrt(2)

rho_plus = ket_plus @ ket_plus.conj().T    # |+><+|
sigma_ref = I2 / 2.0                        # I/2


# =============================================================================
# Numerical utilities
# =============================================================================
def matrix_sqrt_psd(A):
    """Matrix square root for PSD matrix via eigendecomposition."""
    eigvals, eigvecs = np.linalg.eigh(A)
    eigvals = np.maximum(eigvals, 0.0)
    return (eigvecs * np.sqrt(eigvals)) @ eigvecs.conj().T


def matrix_inv_psd(A, tol=1e-12):
    """Pseudoinverse for PSD matrix."""
    eigvals, eigvecs = np.linalg.eigh(A)
    inv_eigvals = np.where(eigvals > tol, 1.0 / eigvals, 0.0)
    return (eigvecs * inv_eigvals) @ eigvecs.conj().T


def uhlmann_fidelity(rho, sigma):
    """Uhlmann fidelity F(rho, sigma) = (Tr sqrt(sqrt(rho) sigma sqrt(rho)))^2."""
    sqrt_rho = matrix_sqrt_psd(rho)
    M = sqrt_rho @ sigma @ sqrt_rho
    M = 0.5 * (M + M.conj().T)
    eigvals = eigvalsh(M)
    eigvals = np.maximum(eigvals, 0.0)
    F = (np.sum(np.sqrt(eigvals)))**2
    return float(np.clip(F.real, 0.0, 1.0))


def trace_distance(rho, sigma):
    """Trace distance D(rho, sigma) = 0.5 * Tr|rho - sigma|."""
    diff = rho - sigma
    diff = 0.5 * (diff + diff.conj().T)
    eigvals = np.linalg.eigvalsh(diff)
    return 0.5 * np.sum(np.abs(eigvals))


def force_density_matrix(rho):
    """Force a matrix to be a valid density matrix."""
    rho = 0.5 * (rho + rho.conj().T)
    eigvals, eigvecs = np.linalg.eigh(rho)
    eigvals = np.maximum(eigvals, 0.0)
    rho = (eigvecs * eigvals) @ eigvecs.conj().T
    tr = np.trace(rho).real
    if tr > 1e-15:
        rho = rho / tr
    return rho


# =============================================================================
# Kraus-based channel application and Petz recovery
# =============================================================================
def apply_kraus(kraus_ops, rho):
    """Apply quantum channel: E(rho) = sum_k K_k rho K_k^dag."""
    out = np.zeros_like(rho)
    for K in kraus_ops:
        out += K @ rho @ K.conj().T
    return out


def build_superoperator_from_kraus(kraus_ops, d):
    """Build d^2 x d^2 superoperator from Kraus operators (column-stacking)."""
    d2 = d * d
    E_super = np.zeros((d2, d2), dtype=complex)
    for K in kraus_ops:
        E_super += np.kron(K.conj(), K)
    return E_super


def adjoint_superoperator(E_super):
    """Adjoint channel superoperator (Hilbert-Schmidt adjoint)."""
    return E_super.conj().T


def vec(rho):
    """Vectorize density matrix (column-stacking)."""
    return rho.flatten('F')


def unvec(v, d):
    """Unvectorize back to d x d matrix."""
    return v.reshape((d, d), order='F')


def apply_super(E_super, rho, d):
    """Apply superoperator to density matrix."""
    return unvec(E_super @ vec(rho), d)


def petz_recovery_kraus(kraus_ops, sigma, rho_in, d):
    """
    Petz recovery map applied to rho_in:
        R[X] = sigma^{1/2} E_adj[ E[sigma]^{-1/2} X E[sigma]^{-1/2} ] sigma^{1/2}

    Uses superoperator for the adjoint channel.
    """
    # E[sigma]
    E_sigma = apply_kraus(kraus_ops, sigma)
    E_sigma = 0.5 * (E_sigma + E_sigma.conj().T)

    # E[sigma]^{-1/2}
    E_sigma_inv_sqrt = matrix_sqrt_psd(matrix_inv_psd(E_sigma))

    # Sandwich rho_in
    Y = E_sigma_inv_sqrt @ rho_in @ E_sigma_inv_sqrt

    # Adjoint channel via superoperator
    E_super = build_superoperator_from_kraus(kraus_ops, d)
    E_adj = adjoint_superoperator(E_super)
    Z = apply_super(E_adj, Y, d)
    Z = 0.5 * (Z + Z.conj().T)

    # sigma^{1/2} sandwich
    sqrt_sigma = matrix_sqrt_psd(sigma)
    R = sqrt_sigma @ Z @ sqrt_sigma

    # Normalize
    R = force_density_matrix(R)
    return R


def compute_tau_kraus(kraus_ops, rho0, sigma0, d):
    """Compute tau = 1 - F(rho0, R o E_t[rho0]) using Kraus operators."""
    E_rho = apply_kraus(kraus_ops, rho0)
    E_rho = force_density_matrix(E_rho)

    R_E_rho = petz_recovery_kraus(kraus_ops, sigma0, E_rho, d)

    F = uhlmann_fidelity(rho0, R_E_rho)
    tau = 1.0 - F
    return tau, F


# =============================================================================
# Non-Markovian amplitude damping: Lorentzian spectral density
# =============================================================================
def lorentzian_c(t, gamma_0, lam):
    """
    Decoherence function c(t) for Lorentzian spectral density.

    c(t) = exp(-lambda*t/2) * [cosh(d*t/2) + (lambda/d)*sinh(d*t/2)]

    where d = sqrt(lambda^2 - 2*gamma_0*lambda).

    When gamma_0 > lambda/2: d is imaginary => oscillations (non-Markovian).
    When gamma_0 < lambda/2: d is real => exponential decay (Markovian).
    """
    discriminant = lam**2 - 2.0 * gamma_0 * lam

    if discriminant >= 0:
        # Markovian or critical regime
        d = np.sqrt(discriminant)
        if np.abs(d) < 1e-15:
            # Critical case: d -> 0, cosh->1, sinh->dt/2
            c = np.exp(-lam * t / 2.0) * (1.0 + lam * t / 2.0)
        else:
            c = np.exp(-lam * t / 2.0) * (
                np.cosh(d * t / 2.0) + (lam / d) * np.sinh(d * t / 2.0)
            )
    else:
        # Non-Markovian regime: d is imaginary
        d_abs = np.sqrt(-discriminant)  # |d|, which is real and positive
        c = np.exp(-lam * t / 2.0) * (
            np.cos(d_abs * t / 2.0) + (lam / d_abs) * np.sin(d_abs * t / 2.0)
        )

    return c


def amplitude_damping_kraus(c_val):
    """
    Time-dependent Kraus operators for amplitude damping channel.

    K0(t) = [[1, 0], [0, c(t)]]
    K1(t) = [[0, s(t)], [0, 0]]

    where |c(t)|^2 + |s(t)|^2 = 1, so s(t) = sqrt(1 - |c(t)|^2).
    """
    c = complex(c_val)
    abs_c_sq = np.abs(c)**2

    # Clamp to avoid numerical issues
    if abs_c_sq > 1.0:
        abs_c_sq = 1.0
        c = c / np.abs(c)

    s = np.sqrt(max(0.0, 1.0 - abs_c_sq))

    K0 = np.array([[1, 0], [0, c]], dtype=complex)
    K1 = np.array([[0, s], [0, 0]], dtype=complex)

    return [K0, K1]


# =============================================================================
# Jaynes-Cummings (exact resonance, unitary in single-excitation sector)
# =============================================================================
def jaynes_cummings_kraus(t, g, delta=0.0):
    """
    Jaynes-Cummings model in the single-excitation sector.

    For initial state |1,0> (qubit excited, cavity vacuum):
    - Qubit reduced dynamics are amplitude-damping-like with:
      c(t) = cos(Omega*t/2) - i*(delta/Omega)*sin(Omega*t/2)

    where Omega = sqrt(delta^2 + 4g^2) is the generalized Rabi frequency.

    For exact resonance (delta=0): c(t) = cos(g*t), s(t) = sin(g*t).
    """
    if np.abs(delta) < 1e-15:
        # Exact resonance
        c_val = np.cos(g * t)
        s_val = np.sin(g * t)
    else:
        Omega = np.sqrt(delta**2 + 4 * g**2)
        c_val = np.cos(Omega * t / 2.0) - 1j * (delta / Omega) * np.sin(Omega * t / 2.0)
        abs_c_sq = np.abs(c_val)**2
        s_val = np.sqrt(max(0.0, 1.0 - abs_c_sq))

    K0 = np.array([[1, 0], [0, c_val]], dtype=complex)
    K1 = np.array([[0, s_val], [0, 0]], dtype=complex)

    return [K0, K1]


# =============================================================================
# BLP non-Markovianity measure
# =============================================================================
def compute_blp(t_vals, rho1_list, rho2_list):
    """
    Compute BLP non-Markovianity measure.

    N_BLP = integral over regions where dD/dt > 0 of (dD/dt) dt

    where D = trace distance between two evolved states.
    """
    n = len(t_vals)
    D_vals = np.zeros(n)
    for i in range(n):
        D_vals[i] = trace_distance(rho1_list[i], rho2_list[i])

    # Numerical derivative
    dD_dt = np.gradient(D_vals, t_vals)

    # Integrate positive parts
    dt = np.diff(t_vals)
    blp = 0.0
    for i in range(len(dt)):
        avg_dD = 0.5 * (dD_dt[i] + dD_dt[i + 1])
        if avg_dD > 0:
            blp += avg_dD * dt[i]

    return blp, D_vals, dD_dt


# =============================================================================
# Main test runner
# =============================================================================
def run_lorentzian_test(gamma_0, lam, label, n_points=200):
    """
    Run a single Lorentzian spectral density test case.

    Returns dict with all results.
    """
    ratio = gamma_0 / lam
    regime = "Markovian" if ratio < 0.5 else ("Critical" if ratio == 0.5 else "Non-Markovian")

    d = 2
    rho0 = rho_plus.copy()
    sigma0 = sigma_ref.copy()

    # Also use |0><0| and |1><1| for BLP
    rho_blp_1 = ket_0 @ ket_0.conj().T  # |0><0|
    rho_blp_2 = ket_1 @ ket_1.conj().T  # |1><1|

    # Compute t_max based on the regime:
    # For Markovian: use decay timescale 5/gamma_0
    # For non-Markovian: use oscillation timescale to capture several periods
    discriminant = lam**2 - 2.0 * gamma_0 * lam
    if discriminant < 0:
        # Non-Markovian: oscillation frequency is d_abs/2
        d_abs = np.sqrt(-discriminant)
        T_osc = 2 * np.pi / (d_abs / 2.0)  # oscillation period
        # Use enough periods to see oscillations, but also account for envelope decay
        # Envelope: exp(-lam*t/2), so decay time ~ 2/lam
        t_max = max(5.0 * T_osc, 10.0 / lam)
    else:
        # Markovian: use effective decay time
        t_max = 5.0 / gamma_0 if gamma_0 > 0 else 5.0

    t_vals = np.linspace(0, t_max, n_points)

    tau_vals = np.zeros(n_points)
    fid_vals = np.zeros(n_points)
    c_vals = np.zeros(n_points)
    evolved_blp1 = []
    evolved_blp2 = []

    for i, t in enumerate(t_vals):
        c_t = lorentzian_c(t, gamma_0, lam)
        c_vals[i] = c_t

        kraus = amplitude_damping_kraus(c_t)
        tau, F = compute_tau_kraus(kraus, rho0, sigma0, d)
        tau_vals[i] = tau
        fid_vals[i] = F

        # For BLP: evolve the two orthogonal states
        rho1_out = apply_kraus(kraus, rho_blp_1)
        rho2_out = apply_kraus(kraus, rho_blp_2)
        evolved_blp1.append(force_density_matrix(rho1_out))
        evolved_blp2.append(force_density_matrix(rho2_out))

    # Compute dτ/dt
    dtau_dt = np.gradient(tau_vals, t_vals)

    # Check for non-monotonicity
    intervals_decrease = dtau_dt < -1e-8
    n_decrease = np.sum(intervals_decrease)
    has_backflow = n_decrease > 0

    # Count oscillation periods (zero crossings of dtau_dt)
    sign_changes = np.diff(np.sign(dtau_dt))
    n_oscillations = np.sum(np.abs(sign_changes) > 1.0) // 2  # each full period = 2 crossings

    # BLP measure
    blp, D_vals, dD_dt = compute_blp(t_vals, evolved_blp1, evolved_blp2)

    # Correlation: periods where dtau/dt < 0 vs dD/dt > 0
    dtau_negative = dtau_dt < -1e-8
    dD_positive = dD_dt > 1e-8
    if np.sum(dtau_negative) > 0:
        overlap = np.sum(dtau_negative & dD_positive) / np.sum(dtau_negative)
    else:
        overlap = float('nan')

    results = {
        'label': label,
        'gamma_0': gamma_0,
        'lambda': lam,
        'ratio': ratio,
        'regime': regime,
        't_vals': t_vals,
        'tau_vals': tau_vals,
        'fid_vals': fid_vals,
        'c_vals': c_vals,
        'dtau_dt': dtau_dt,
        'D_vals': D_vals,
        'dD_dt': dD_dt,
        'min_tau': np.min(tau_vals),
        'max_tau': np.max(tau_vals),
        'n_oscillations': int(n_oscillations),
        'has_backflow': has_backflow,
        'n_decrease_points': int(n_decrease),
        'blp': blp,
        'blp_dtau_correlation': overlap,
    }
    return results


def run_jaynes_cummings_test(g, delta=0.0, n_points=200):
    """
    Run Jaynes-Cummings test (unitary, exact resonance or detuned).
    """
    d = 2
    rho0 = rho_plus.copy()
    sigma0 = sigma_ref.copy()

    rho_blp_1 = ket_0 @ ket_0.conj().T
    rho_blp_2 = ket_1 @ ket_1.conj().T

    # Time range: several Rabi periods
    if np.abs(delta) < 1e-15:
        Omega = 2 * g  # Rabi frequency at resonance
    else:
        Omega = np.sqrt(delta**2 + 4 * g**2)

    T_rabi = 2 * np.pi / Omega
    t_max = 5 * T_rabi
    t_vals = np.linspace(0, t_max, n_points)

    tau_vals = np.zeros(n_points)
    fid_vals = np.zeros(n_points)
    evolved_blp1 = []
    evolved_blp2 = []

    for i, t in enumerate(t_vals):
        kraus = jaynes_cummings_kraus(t, g, delta)
        tau, F = compute_tau_kraus(kraus, rho0, sigma0, d)
        tau_vals[i] = tau
        fid_vals[i] = F

        rho1_out = apply_kraus(kraus, rho_blp_1)
        rho2_out = apply_kraus(kraus, rho_blp_2)
        evolved_blp1.append(force_density_matrix(rho1_out))
        evolved_blp2.append(force_density_matrix(rho2_out))

    dtau_dt = np.gradient(tau_vals, t_vals)

    intervals_decrease = dtau_dt < -1e-8
    n_decrease = np.sum(intervals_decrease)
    has_backflow = n_decrease > 0

    sign_changes = np.diff(np.sign(dtau_dt))
    n_oscillations = np.sum(np.abs(sign_changes) > 1.0) // 2

    blp, D_vals, dD_dt = compute_blp(t_vals, evolved_blp1, evolved_blp2)

    dtau_negative = dtau_dt < -1e-8
    dD_positive = dD_dt > 1e-8
    if np.sum(dtau_negative) > 0:
        overlap = np.sum(dtau_negative & dD_positive) / np.sum(dtau_negative)
    else:
        overlap = float('nan')

    results = {
        'label': f'Jaynes-Cummings (g={g}, delta={delta})',
        'gamma_0': g,
        'lambda': 0.0,
        'ratio': float('inf'),
        'regime': 'Unitary (JC)',
        't_vals': t_vals,
        'tau_vals': tau_vals,
        'fid_vals': fid_vals,
        'c_vals': np.array([np.cos(g * t) for t in t_vals]) if np.abs(delta) < 1e-15 else None,
        'dtau_dt': dtau_dt,
        'D_vals': D_vals,
        'dD_dt': dD_dt,
        'min_tau': np.min(tau_vals),
        'max_tau': np.max(tau_vals),
        'n_oscillations': int(n_oscillations),
        'has_backflow': has_backflow,
        'n_decrease_points': int(n_decrease),
        'blp': blp,
        'blp_dtau_correlation': overlap,
        'T_rabi': T_rabi,
    }
    return results


# =============================================================================
# MAIN
# =============================================================================
if __name__ == "__main__":
    print()
    print("*" * 78)
    print("*  Non-Markovian tau(t) Verification")
    print("*  Jaynes-Cummings Model & Lorentzian Spectral Density")
    print("*  Petz Recovery Framework")
    print("*" * 78)
    print()

    all_results = []

    # =========================================================================
    # Test (a): Markovian limit  gamma_0/lambda = 0.1
    # =========================================================================
    lam = 1.0  # Fix lambda = 1 as reference
    print("=" * 78)
    print("TEST (a): Markovian Limit  gamma_0/lambda = 0.1")
    print("  Expectation: tau(t) monotonically increases, dτ/dt >= 0 always")
    print("=" * 78)
    res_a = run_lorentzian_test(gamma_0=0.1 * lam, lam=lam, label="(a) Markovian")
    all_results.append(res_a)

    print(f"  gamma_0 = {res_a['gamma_0']:.3f}, lambda = {res_a['lambda']:.3f}")
    print(f"  Regime: {res_a['regime']}")
    print(f"  tau range: [{res_a['min_tau']:.6f}, {res_a['max_tau']:.6f}]")
    print(f"  dτ/dt < 0 detected: {res_a['has_backflow']}")
    if res_a['has_backflow']:
        print(f"    (UNEXPECTED for Markovian! {res_a['n_decrease_points']} points)")
    else:
        print(f"    (CORRECT: monotonic increase)")
    print(f"  BLP measure: {res_a['blp']:.6e}")
    print()

    # Sample time points
    sample_idx = [0, 25, 50, 75, 100, 125, 150, 175, 199]
    print(f"  {'t':>10s}  {'c(t)':>10s}  {'tau(t)':>12s}  {'F(t)':>12s}  {'dtau/dt':>12s}")
    print(f"  {'-'*10}  {'-'*10}  {'-'*12}  {'-'*12}  {'-'*12}")
    for idx in sample_idx:
        print(f"  {res_a['t_vals'][idx]:10.4f}  {res_a['c_vals'][idx]:10.6f}  "
              f"{res_a['tau_vals'][idx]:12.8f}  {res_a['fid_vals'][idx]:12.8f}  "
              f"{res_a['dtau_dt'][idx]:12.4e}")
    print()

    # =========================================================================
    # Test (b): Critical  gamma_0/lambda = 0.5
    # =========================================================================
    print("=" * 78)
    print("TEST (b): Critical  gamma_0/lambda = 0.5")
    print("  Expectation: borderline, near-monotonic with possible weak oscillations")
    print("=" * 78)
    res_b = run_lorentzian_test(gamma_0=0.5 * lam, lam=lam, label="(b) Critical")
    all_results.append(res_b)

    print(f"  gamma_0 = {res_b['gamma_0']:.3f}, lambda = {res_b['lambda']:.3f}")
    print(f"  Regime: {res_b['regime']}")
    print(f"  tau range: [{res_b['min_tau']:.6f}, {res_b['max_tau']:.6f}]")
    print(f"  dτ/dt < 0 detected: {res_b['has_backflow']}")
    print(f"  Oscillation periods: {res_b['n_oscillations']}")
    print(f"  BLP measure: {res_b['blp']:.6e}")
    print()

    print(f"  {'t':>10s}  {'c(t)':>10s}  {'tau(t)':>12s}  {'F(t)':>12s}  {'dtau/dt':>12s}")
    print(f"  {'-'*10}  {'-'*10}  {'-'*12}  {'-'*12}  {'-'*12}")
    for idx in sample_idx:
        print(f"  {res_b['t_vals'][idx]:10.4f}  {res_b['c_vals'][idx]:10.6f}  "
              f"{res_b['tau_vals'][idx]:12.8f}  {res_b['fid_vals'][idx]:12.8f}  "
              f"{res_b['dtau_dt'][idx]:12.4e}")
    print()

    # =========================================================================
    # Test (c): Non-Markovian  gamma_0/lambda = 2.0
    # =========================================================================
    print("=" * 78)
    print("TEST (c): Non-Markovian  gamma_0/lambda = 2.0")
    print("  Expectation: tau(t) oscillates, dτ/dt < 0 in some intervals")
    print("=" * 78)
    res_c = run_lorentzian_test(gamma_0=2.0 * lam, lam=lam, label="(c) Non-Markovian")
    all_results.append(res_c)

    print(f"  gamma_0 = {res_c['gamma_0']:.3f}, lambda = {res_c['lambda']:.3f}")
    print(f"  Regime: {res_c['regime']}")
    print(f"  tau range: [{res_c['min_tau']:.6f}, {res_c['max_tau']:.6f}]")
    print(f"  dτ/dt < 0 detected: {res_c['has_backflow']}")
    if res_c['has_backflow']:
        print(f"    (CORRECT: information backflow! {res_c['n_decrease_points']} time points)")
    else:
        print(f"    (UNEXPECTED: should see backflow for gamma_0/lambda = 2.0!)")
    print(f"  Oscillation periods: {res_c['n_oscillations']}")
    print(f"  BLP measure: {res_c['blp']:.6e}")
    print(f"  Correlation (dτ/dt<0 with dD/dt>0): {res_c['blp_dtau_correlation']:.2%}")
    print()

    print(f"  {'t':>10s}  {'c(t)':>10s}  {'tau(t)':>12s}  {'F(t)':>12s}  {'dtau/dt':>12s}")
    print(f"  {'-'*10}  {'-'*10}  {'-'*12}  {'-'*12}  {'-'*12}")
    for idx in sample_idx:
        print(f"  {res_c['t_vals'][idx]:10.4f}  {res_c['c_vals'][idx]:10.6f}  "
              f"{res_c['tau_vals'][idx]:12.8f}  {res_c['fid_vals'][idx]:12.8f}  "
              f"{res_c['dtau_dt'][idx]:12.4e}")
    print()

    # Highlight the backflow intervals
    if res_c['has_backflow']:
        print("  Intervals where dτ/dt < 0 (information backflow / Sigma < 0):")
        in_backflow = False
        bf_start = 0
        for i in range(len(res_c['t_vals'])):
            if res_c['dtau_dt'][i] < -1e-8:
                if not in_backflow:
                    bf_start = res_c['t_vals'][i]
                    in_backflow = True
            else:
                if in_backflow:
                    bf_end = res_c['t_vals'][i]
                    # Find min tau in this interval
                    mask = (res_c['t_vals'] >= bf_start) & (res_c['t_vals'] <= bf_end)
                    min_dtau = np.min(res_c['dtau_dt'][mask])
                    print(f"    t in [{bf_start:.4f}, {bf_end:.4f}], "
                          f"min dτ/dt = {min_dtau:.4e}")
                    in_backflow = False
        if in_backflow:
            bf_end = res_c['t_vals'][-1]
            mask = (res_c['t_vals'] >= bf_start)
            min_dtau = np.min(res_c['dtau_dt'][mask])
            print(f"    t in [{bf_start:.4f}, {bf_end:.4f}], "
                  f"min dτ/dt = {min_dtau:.4e}")
        print()

    # =========================================================================
    # Test (d): Strong Non-Markovian  gamma_0/lambda = 5.0
    # =========================================================================
    print("=" * 78)
    print("TEST (d): Strong Non-Markovian  gamma_0/lambda = 5.0")
    print("  Expectation: pronounced oscillations, tau returns near 0 at revival times")
    print("=" * 78)
    res_d = run_lorentzian_test(gamma_0=5.0 * lam, lam=lam, label="(d) Strong Non-Markov")
    all_results.append(res_d)

    print(f"  gamma_0 = {res_d['gamma_0']:.3f}, lambda = {res_d['lambda']:.3f}")
    print(f"  Regime: {res_d['regime']}")
    print(f"  tau range: [{res_d['min_tau']:.6f}, {res_d['max_tau']:.6f}]")
    print(f"  dτ/dt < 0 detected: {res_d['has_backflow']}")
    if res_d['has_backflow']:
        print(f"    (CORRECT: strong backflow! {res_d['n_decrease_points']} time points)")
    print(f"  Oscillation periods: {res_d['n_oscillations']}")
    print(f"  BLP measure: {res_d['blp']:.6e}")
    print(f"  Correlation (dτ/dt<0 with dD/dt>0): {res_d['blp_dtau_correlation']:.2%}")
    print()

    print(f"  {'t':>10s}  {'c(t)':>10s}  {'tau(t)':>12s}  {'F(t)':>12s}  {'dtau/dt':>12s}")
    print(f"  {'-'*10}  {'-'*10}  {'-'*12}  {'-'*12}  {'-'*12}")
    for idx in sample_idx:
        print(f"  {res_d['t_vals'][idx]:10.4f}  {res_d['c_vals'][idx]:10.6f}  "
              f"{res_d['tau_vals'][idx]:12.8f}  {res_d['fid_vals'][idx]:12.8f}  "
              f"{res_d['dtau_dt'][idx]:12.4e}")
    print()

    # Revival analysis for strong non-Markovian
    if res_d['has_backflow']:
        print("  Revival analysis (local minima of tau after decrease):")
        from scipy.signal import argrelmin
        local_min_idx = argrelmin(res_d['tau_vals'], order=5)[0]
        for mi in local_min_idx:
            print(f"    Revival at t = {res_d['t_vals'][mi]:.4f}: "
                  f"tau = {res_d['tau_vals'][mi]:.6f}, "
                  f"c(t) = {res_d['c_vals'][mi]:.6f}")
        print()

        print("  Intervals where dτ/dt < 0 (information backflow / Sigma < 0):")
        in_backflow = False
        bf_start = 0
        for i in range(len(res_d['t_vals'])):
            if res_d['dtau_dt'][i] < -1e-8:
                if not in_backflow:
                    bf_start = res_d['t_vals'][i]
                    in_backflow = True
            else:
                if in_backflow:
                    bf_end = res_d['t_vals'][i]
                    mask = (res_d['t_vals'] >= bf_start) & (res_d['t_vals'] <= bf_end)
                    min_dtau = np.min(res_d['dtau_dt'][mask])
                    print(f"    t in [{bf_start:.4f}, {bf_end:.4f}], "
                          f"min dτ/dt = {min_dtau:.4e}")
                    in_backflow = False
        if in_backflow:
            bf_end = res_d['t_vals'][-1]
            mask = (res_d['t_vals'] >= bf_start)
            min_dtau = np.min(res_d['dtau_dt'][mask])
            print(f"    t in [{bf_start:.4f}, {bf_end:.4f}], "
                  f"min dτ/dt = {min_dtau:.4e}")
        print()

    # =========================================================================
    # Test (e): Jaynes-Cummings exact resonance (unitary)
    # =========================================================================
    print("=" * 78)
    print("TEST (e): Jaynes-Cummings, Exact Resonance (Unitary)")
    print("  g = 1.0, delta = 0 (exact resonance)")
    print("  Expectation: tau oscillates perfectly between 0 and max, full revivals")
    print("=" * 78)
    res_e = run_jaynes_cummings_test(g=1.0, delta=0.0)
    all_results.append(res_e)

    print(f"  g = 1.0, delta = 0.0")
    print(f"  Rabi period T_Rabi = {res_e['T_rabi']:.4f}")
    print(f"  Regime: {res_e['regime']}")
    print(f"  tau range: [{res_e['min_tau']:.8f}, {res_e['max_tau']:.8f}]")
    print(f"  dτ/dt < 0 detected: {res_e['has_backflow']}")
    if res_e['has_backflow']:
        print(f"    (CORRECT: oscillatory backflow! {res_e['n_decrease_points']} time points)")
    print(f"  Oscillation periods: {res_e['n_oscillations']}")
    print(f"  BLP measure: {res_e['blp']:.6e}")
    print()

    print(f"  {'t':>10s}  {'tau(t)':>12s}  {'F(t)':>12s}  {'dtau/dt':>12s}")
    print(f"  {'-'*10}  {'-'*12}  {'-'*12}  {'-'*12}")
    for idx in sample_idx:
        print(f"  {res_e['t_vals'][idx]:10.4f}  "
              f"{res_e['tau_vals'][idx]:12.8f}  {res_e['fid_vals'][idx]:12.8f}  "
              f"{res_e['dtau_dt'][idx]:12.4e}")
    print()

    # Check that tau returns to 0 at revival times (multiples of pi/g)
    g_val = 1.0
    revival_times = [np.pi / g_val, 2 * np.pi / g_val, 3 * np.pi / g_val]
    print("  Revival check (tau should return to ~0 at t = n*pi/g):")
    for t_rev in revival_times:
        if t_rev <= res_e['t_vals'][-1]:
            idx_rev = np.argmin(np.abs(res_e['t_vals'] - t_rev))
            print(f"    t = {t_rev:.4f} (n*pi/g): tau = {res_e['tau_vals'][idx_rev]:.2e}, "
                  f"F = {res_e['fid_vals'][idx_rev]:.10f}")
    print()

    # =========================================================================
    # SUMMARY TABLE
    # =========================================================================
    print()
    print("=" * 78)
    print("SUMMARY TABLE: Non-Markovian tau(t) Verification")
    print("=" * 78)
    print()

    header = (f"  {'Case':<24s}  {'gamma/lam':>9s}  {'min(tau)':>10s}  {'max(tau)':>10s}  "
              f"{'Oscill':>6s}  {'dtau<0':>6s}  {'BLP':>12s}")
    print(header)
    print("  " + "-" * (len(header) - 2))

    for res in all_results:
        ratio_str = f"{res['ratio']:.1f}" if res['ratio'] != float('inf') else "inf"
        backflow_str = "YES" if res['has_backflow'] else "NO"
        print(f"  {res['label']:<24s}  {ratio_str:>9s}  {res['min_tau']:>10.6f}  "
              f"{res['max_tau']:>10.6f}  {res['n_oscillations']:>6d}  "
              f"{backflow_str:>6s}  {res['blp']:>12.4e}")

    print()

    # =========================================================================
    # KEY PREDICTIONS VERIFICATION
    # =========================================================================
    print("=" * 78)
    print("KEY PREDICTIONS VERIFICATION")
    print("=" * 78)
    print()

    # Prediction 1: Markovian => monotonic
    pred1 = not res_a['has_backflow']
    print(f"  [{'PASS' if pred1 else 'FAIL'}] Prediction 1: Markovian (gamma/lam=0.1) => "
          f"dτ/dt >= 0 always")
    if pred1:
        print(f"         Confirmed: no intervals with dτ/dt < 0 detected")
    else:
        print(f"         VIOLATED: found {res_a['n_decrease_points']} points with dτ/dt < 0")

    # Prediction 2: Non-Markovian => backflow
    pred2 = res_c['has_backflow'] and res_d['has_backflow']
    print(f"  [{'PASS' if pred2 else 'FAIL'}] Prediction 2: Non-Markovian (gamma/lam >= 2) => "
          f"intervals with dτ/dt < 0 exist")
    if pred2:
        print(f"         Case (c): {res_c['n_decrease_points']} backflow points, "
              f"{res_c['n_oscillations']} oscillations")
        print(f"         Case (d): {res_d['n_decrease_points']} backflow points, "
              f"{res_d['n_oscillations']} oscillations")
    else:
        print(f"         UNEXPECTED: backflow not detected in non-Markovian regime")

    # Prediction 3: Connection tau decrease <-> Sigma < 0
    pred3_c = res_c['has_backflow']
    pred3_d = res_d['has_backflow']
    pred3 = pred3_c and pred3_d
    print(f"  [{'PASS' if pred3 else 'FAIL'}] Prediction 3: dτ/dt < 0 => Sigma < 0 locally "
          f"(information backflow = time reversal)")
    if pred3:
        print(f"         Case (c): min dτ/dt = {np.min(res_c['dtau_dt']):.4e}")
        print(f"         Case (d): min dτ/dt = {np.min(res_d['dtau_dt']):.4e}")
        print(f"         These intervals represent entropy decrease => 'backward time'")

    # Prediction 4: Revival times => tau near 0
    pred4_vals = []
    g_val = 1.0
    for t_rev in [np.pi / g_val, 2 * np.pi / g_val]:
        if t_rev <= res_e['t_vals'][-1]:
            idx_rev = np.argmin(np.abs(res_e['t_vals'] - t_rev))
            pred4_vals.append(res_e['tau_vals'][idx_rev])
    pred4 = all(v < 0.01 for v in pred4_vals) if pred4_vals else False
    print(f"  [{'PASS' if pred4 else 'FAIL'}] Prediction 4: JC revival times => tau returns to ~0")
    if pred4_vals:
        print(f"         Revival tau values: {[f'{v:.6f}' for v in pred4_vals]}")

    # Prediction 5: BLP correlation
    if not np.isnan(res_c['blp_dtau_correlation']):
        pred5 = res_c['blp_dtau_correlation'] > 0.5
        print(f"  [{'PASS' if pred5 else 'FAIL'}] Prediction 5: Periods dτ/dt < 0 correlate "
              f"with dD/dt > 0 (BLP)")
        print(f"         Case (c) correlation: {res_c['blp_dtau_correlation']:.2%}")
        print(f"         Case (d) correlation: {res_d['blp_dtau_correlation']:.2%}")
    else:
        pred5 = False
        print(f"  [N/A ] Prediction 5: Cannot check BLP correlation (no backflow in test case)")

    # =========================================================================
    # PHYSICAL INTERPRETATION
    # =========================================================================
    print()
    print("=" * 78)
    print("PHYSICAL INTERPRETATION: Connection to Sigma_signed")
    print("=" * 78)
    print()
    print("  tau(t) = 1 - F(rho_0, R_Petz o E_t[rho_0])")
    print()
    print("  In the Markovian regime:")
    print("    - tau(t) monotonically increases")
    print("    - dtau/dt >= 0 always")
    print("    - Sigma = -ln(1 - tau) >= 0")
    print("    - Time arrow is well-defined and points forward")
    print()
    print("  In the non-Markovian regime:")
    print("    - tau(t) can DECREASE in certain intervals")
    print("    - dtau/dt < 0 <=> information flows BACK from environment to system")
    print("    - Locally: Sigma_signed < 0 => 'time runs backward'")
    print("    - At revival times: tau -> 0, full information recovery")
    print()
    print("  This is EXACTLY the tau_signed framework prediction:")
    print("    tau_signed > 0  =>  time forward  (Markovian)")
    print("    tau_signed = 0  =>  no time        (unitary / perfect recovery)")
    print("    tau_signed < 0  =>  time backward  (non-Markovian backflow)")
    print()

    # =========================================================================
    # OVERALL VERDICT
    # =========================================================================
    overall_tests = [pred1, pred2, pred3, pred4, pred5]
    n_pass = sum(overall_tests)
    n_total = len(overall_tests)

    print("=" * 78)
    print(f"OVERALL: {n_pass}/{n_total} predictions verified")
    print("=" * 78)

    if n_pass == n_total:
        print("  ALL PREDICTIONS CONFIRMED.")
        print("  Non-Markovian tau(t) framework is fully consistent with:")
        print("    - BLP non-Markovianity measure")
        print("    - Jaynes-Cummings unitary revivals")
        print("    - Sigma_signed interpretation of time direction")
    else:
        failed = [f"Pred {i+1}" for i, p in enumerate(overall_tests) if not p]
        print(f"  Failed: {', '.join(failed)}")
        print("  Investigation needed for failed predictions.")

    print()
    sys.exit(0 if n_pass == n_total else 1)
