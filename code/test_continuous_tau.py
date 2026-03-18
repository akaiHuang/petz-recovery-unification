#!/usr/bin/env python3
"""
Continuous-time tau(t) verification for Lindbladian dynamics.
=============================================================

Three critical tests:
  1. Pure dephasing (analytical benchmark)
  2. Pure unitary evolution (tau must vanish identically)
  3. Amplitude damping (Markovian monotonicity)

Framework:
  - Build Lindbladian superoperator L in column-stacking vec representation
  - Channel E_t = expm(L*t)  (superoperator matrix)
  - Petz recovery map R_{sigma, E_t}
  - tau(t) = 1 - F(rho_0, R o E_t [rho_0])
  - Uhlmann fidelity F(rho, sigma) = (Tr sqrt(sqrt(rho) sigma sqrt(rho)))^2

Author: Sheng-Kai Huang
Date: 2026-03-16
"""

import numpy as np
from scipy.linalg import expm, sqrtm, eigvalsh
import sys

# =============================================================================
# Pauli matrices and helpers
# =============================================================================
I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
sp = np.array([[0, 0], [1, 0]], dtype=complex)  # sigma_+ = |1><0|
sm = np.array([[0, 1], [0, 0]], dtype=complex)  # sigma_- = |0><1|

# States
ket_0 = np.array([[1], [0]], dtype=complex)
ket_1 = np.array([[0], [1]], dtype=complex)
ket_plus = (ket_0 + ket_1) / np.sqrt(2)

rho_plus = ket_plus @ ket_plus.conj().T   # |+><+|
rho_1 = ket_1 @ ket_1.conj().T            # |1><1|
sigma_max_mix = I2 / 2.0                    # I/2


def vec(rho):
    """Vectorize density matrix using column-stacking (Fortran order)."""
    return rho.flatten('F')


def unvec(v, d):
    """Un-vectorize back to d x d matrix (column-stacking)."""
    return v.reshape((d, d), order='F')


def build_lindbladian(d, H, L_ops, gammas):
    """
    Build the Lindbladian superoperator in column-stacking vec representation.

    L_super = -i(I x H - H^T x I)
              + sum_k gamma_k * (L_k^* x L_k
                                 - 0.5 * I x L_k^dag L_k
                                 - 0.5 * (L_k^dag L_k)^T x I)

    where x denotes Kronecker product.
    """
    d2 = d * d
    Id = np.eye(d, dtype=complex)

    # Hamiltonian part
    L_super = -1j * (np.kron(Id, H) - np.kron(H.T, Id))

    # Dissipative part
    for Lk, gk in zip(L_ops, gammas):
        LdL = Lk.conj().T @ Lk
        L_super += gk * (
            np.kron(Lk.conj(), Lk)
            - 0.5 * np.kron(Id, LdL)
            - 0.5 * np.kron(LdL.T, Id)
        )

    return L_super


def apply_channel(E_super, rho, d):
    """Apply channel superoperator: E_t[rho] = unvec(E_super @ vec(rho))."""
    return unvec(E_super @ vec(rho), d)


def adjoint_channel_super(E_super):
    """
    Adjoint channel (Hilbert-Schmidt adjoint) superoperator.
    For column-stacking: E_adj = E_super^{dagger} (conjugate transpose).
    """
    return E_super.conj().T


def matrix_sqrt_psd(A):
    """
    Compute matrix square root for a PSD matrix, forcing eigenvalues >= 0.
    More numerically stable than sqrtm for near-singular matrices.
    """
    eigvals, eigvecs = np.linalg.eigh(A)
    eigvals = np.maximum(eigvals, 0.0)
    sqrt_eigvals = np.sqrt(eigvals)
    return (eigvecs * sqrt_eigvals) @ eigvecs.conj().T


def matrix_inv_psd(A, tol=1e-12):
    """
    Pseudoinverse-style inverse for PSD matrix. Inverts eigenvalues > tol.
    """
    eigvals, eigvecs = np.linalg.eigh(A)
    inv_eigvals = np.where(eigvals > tol, 1.0 / eigvals, 0.0)
    return (eigvecs * inv_eigvals) @ eigvecs.conj().T


def uhlmann_fidelity(rho, sigma):
    """
    Uhlmann fidelity F(rho, sigma) = (Tr sqrt(sqrt(rho) sigma sqrt(rho)))^2.
    Uses eigendecomposition for numerical stability.
    """
    sqrt_rho = matrix_sqrt_psd(rho)
    M = sqrt_rho @ sigma @ sqrt_rho

    # Force Hermitian
    M = 0.5 * (M + M.conj().T)

    eigvals = eigvalsh(M)
    eigvals = np.maximum(eigvals, 0.0)

    F = (np.sum(np.sqrt(eigvals)))**2
    # Clamp to [0, 1]
    return float(np.clip(F.real, 0.0, 1.0))


def petz_recovery(E_super, sigma, rho_in, d):
    """
    Apply Petz recovery map to rho_in:
        R[X] = sigma^{1/2}  E_adj[ E[sigma]^{-1/2}  X  E[sigma]^{-1/2} ]  sigma^{1/2}

    Steps:
        1. Compute E[sigma]
        2. Compute E[sigma]^{-1/2}
        3. Apply: Y = E[sigma]^{-1/2}  rho_in  E[sigma]^{-1/2}
        4. Apply adjoint channel: Z = E_adj[Y]
        5. Sandwich: R = sigma^{1/2} Z sigma^{1/2}
    """
    # Step 1: E[sigma]
    E_sigma = apply_channel(E_super, sigma, d)
    E_sigma = 0.5 * (E_sigma + E_sigma.conj().T)  # Force Hermitian

    # Step 2: E[sigma]^{-1/2}
    E_sigma_inv_sqrt = matrix_sqrt_psd(matrix_inv_psd(E_sigma))

    # Step 3: Sandwich rho_in
    Y = E_sigma_inv_sqrt @ rho_in @ E_sigma_inv_sqrt

    # Step 4: Adjoint channel
    E_adj_super = adjoint_channel_super(E_super)
    Z = apply_channel(E_adj_super, Y, d)
    Z = 0.5 * (Z + Z.conj().T)  # Force Hermitian

    # Step 5: sigma^{1/2} sandwich
    sqrt_sigma = matrix_sqrt_psd(sigma)
    R = sqrt_sigma @ Z @ sqrt_sigma

    # Force valid density matrix
    R = 0.5 * (R + R.conj().T)
    R = R / np.trace(R).real  # Normalize

    return R


def compute_tau(E_super, rho0, sigma0, d):
    """
    Compute tau = 1 - F(rho0, R o E_t[rho0]).
    """
    # Apply channel to rho0
    E_rho = apply_channel(E_super, rho0, d)
    E_rho = 0.5 * (E_rho + E_rho.conj().T)

    # Apply Petz recovery
    R_E_rho = petz_recovery(E_super, sigma0, E_rho, d)

    # Compute fidelity
    F = uhlmann_fidelity(rho0, R_E_rho)

    tau = 1.0 - F
    return tau, F, R_E_rho


# =============================================================================
# TEST 1: Pure Dephasing
# =============================================================================
def test_pure_dephasing():
    print("=" * 72)
    print("TEST 1: Pure Dephasing")
    print("  d=2, H=0, L=[sigma_z], gamma=0.1")
    print("  rho_0 = |+><+|, sigma_0 = I/2")
    print("  Analytical: tau(t) = (1 - exp(-4*gamma*t)) / 2")
    print("  [Derivation: dephasing is self-adjoint, E[I/2]=I/2,")
    print("   so Petz R = E_adj, and R o E = E^2 => rho(2t).")
    print("   F(|+><+|, rho(2t)) = (1+exp(-4gt))/2, hence tau = (1-exp(-4gt))/2]")
    print("=" * 72)

    d = 2
    gamma = 0.1
    H = np.zeros((d, d), dtype=complex)
    L_ops = [sz]
    gammas = [gamma]

    rho0 = rho_plus.copy()
    sigma0 = sigma_max_mix.copy()

    # Build Lindbladian
    L_super = build_lindbladian(d, H, L_ops, gammas)

    # Time grid
    t_vals = np.linspace(0, 50, 200)
    tau_numerical = np.zeros(len(t_vals))
    tau_analytical = np.zeros(len(t_vals))
    fidelity_vals = np.zeros(len(t_vals))

    for i, t in enumerate(t_vals):
        E_super = expm(L_super * t)
        tau_num, F, _ = compute_tau(E_super, rho0, sigma0, d)
        tau_numerical[i] = tau_num
        fidelity_vals[i] = F
        tau_analytical[i] = (1.0 - np.exp(-4.0 * gamma * t)) / 2.0

    errors = np.abs(tau_numerical - tau_analytical)
    max_error = np.max(errors)
    mean_error = np.mean(errors)

    print(f"\n  Time range: [0, 50]  ({len(t_vals)} points)")
    print(f"  gamma = {gamma}")
    print()

    # Print sample values
    sample_indices = [0, 10, 25, 50, 100, 150, 199]
    print(f"  {'t':>8s}  {'tau_num':>14s}  {'tau_ana':>14s}  {'error':>12s}  {'F':>12s}")
    print(f"  {'-'*8}  {'-'*14}  {'-'*14}  {'-'*12}  {'-'*12}")
    for idx in sample_indices:
        print(f"  {t_vals[idx]:8.3f}  {tau_numerical[idx]:14.10f}  "
              f"{tau_analytical[idx]:14.10f}  {errors[idx]:12.2e}  {fidelity_vals[idx]:12.10f}")

    print()
    print(f"  tau(0)  = {tau_numerical[0]:.2e}  (should be 0)")
    print(f"  tau(50) = {tau_numerical[-1]:.10f}  (analytical: {tau_analytical[-1]:.10f})")
    print(f"  Max |error| = {max_error:.2e}")
    print(f"  Mean |error| = {mean_error:.2e}")

    # PASS/FAIL
    passed = max_error < 1e-8
    status = "PASS" if passed else "FAIL"
    print(f"\n  >>> Test 1 criterion: max|error| < 1e-8")
    print(f"  >>> Result: {status} (max|error| = {max_error:.2e})")

    return passed, t_vals, tau_numerical, tau_analytical


# =============================================================================
# TEST 2: Pure Unitary Evolution
# =============================================================================
def test_pure_unitary():
    print("\n" + "=" * 72)
    print("TEST 2: Pure Unitary Evolution (CRITICAL)")
    print("  d=2, H = omega*sigma_z/2, omega=1.0, NO dissipation")
    print("  rho_0 = |+><+|, sigma_0 = I/2")
    print("  tau MUST be exactly 0 for all t")
    print("=" * 72)

    d = 2
    omega = 1.0
    H = omega * sz / 2.0
    L_ops = []
    gammas = []

    rho0 = rho_plus.copy()
    sigma0 = sigma_max_mix.copy()

    # Build Lindbladian (pure Hamiltonian)
    L_super = build_lindbladian(d, H, L_ops, gammas)

    # Time grid: 10 Larmor periods = 20*pi
    t_max = 20.0 * np.pi
    t_vals = np.linspace(0, t_max, 500)
    tau_vals = np.zeros(len(t_vals))
    fidelity_vals = np.zeros(len(t_vals))

    for i, t in enumerate(t_vals):
        E_super = expm(L_super * t)
        tau_num, F, _ = compute_tau(E_super, rho0, sigma0, d)
        tau_vals[i] = tau_num
        fidelity_vals[i] = F

    max_tau = np.max(np.abs(tau_vals))
    mean_tau = np.mean(np.abs(tau_vals))

    print(f"\n  Time range: [0, 20*pi]  ({len(t_vals)} points)")
    print(f"  omega = {omega}")
    print()

    # Print sample values
    sample_indices = [0, 50, 100, 200, 300, 400, 499]
    print(f"  {'t':>10s}  {'tau(t)':>16s}  {'F(t)':>16s}")
    print(f"  {'-'*10}  {'-'*16}  {'-'*16}")
    for idx in sample_indices:
        print(f"  {t_vals[idx]:10.4f}  {tau_vals[idx]:16.2e}  {fidelity_vals[idx]:16.14f}")

    print()
    print(f"  max|tau(t)| = {max_tau:.2e}")
    print(f"  mean|tau(t)| = {mean_tau:.2e}")

    # Physical explanation
    print("\n  Physical explanation:")
    print("  For unitary evolution with sigma = I/d:")
    print("    E[sigma] = U sigma U^dag = (I/d) (since U I U^dag = I)")
    print("    Petz recovery R = sigma^{1/2} E_adj[E[sigma]^{-1/2} . E[sigma]^{-1/2}] sigma^{1/2}")
    print("    = (I/sqrt(d)) U^dag [(dI)^{1/2} . (dI)^{1/2}] U (I/sqrt(d))")
    print("    = U^dag . U  (the inverse channel)")
    print("    => R o E = identity => F = 1 => tau = 0")

    # PASS/FAIL
    passed = max_tau < 1e-10
    status = "PASS" if passed else "FAIL"
    print(f"\n  >>> Test 2 criterion: max|tau(t)| < 1e-10")
    print(f"  >>> Result: {status} (max|tau| = {max_tau:.2e})")

    return passed, t_vals, tau_vals


# =============================================================================
# TEST 3: Amplitude Damping
# =============================================================================
def test_amplitude_damping():
    print("\n" + "=" * 72)
    print("TEST 3: Amplitude Damping (Markovian Monotonicity)")
    print("  d=2, H=0, L=[sigma_-], gamma=0.05")
    print("  rho_0 = |1><1|, sigma_0 = I/2")
    print("  Checks: monotone non-decreasing, tau(0)=0, tau in [0,1]")
    print("=" * 72)

    d = 2
    gamma = 0.05
    H = np.zeros((d, d), dtype=complex)
    L_ops = [sm]
    gammas = [gamma]

    rho0 = rho_1.copy()
    sigma0 = sigma_max_mix.copy()

    # Build Lindbladian
    L_super = build_lindbladian(d, H, L_ops, gammas)

    # Time grid
    t_vals = np.linspace(0, 80, 300)
    tau_vals = np.zeros(len(t_vals))
    fidelity_vals = np.zeros(len(t_vals))

    for i, t in enumerate(t_vals):
        E_super = expm(L_super * t)
        tau_num, F, _ = compute_tau(E_super, rho0, sigma0, d)
        tau_vals[i] = tau_num
        fidelity_vals[i] = F

    # Check monotonicity
    diffs = np.diff(tau_vals)
    monotone_violations = np.sum(diffs < -1e-10)
    max_violation = np.min(diffs)

    # Check bounds
    in_bounds = np.all((tau_vals >= -1e-10) & (tau_vals <= 1.0 + 1e-10))
    tau_at_0 = tau_vals[0]
    tau_at_end = tau_vals[-1]

    print(f"\n  Time range: [0, 80]  ({len(t_vals)} points)")
    print(f"  gamma = {gamma}")
    print()

    # Print sample values
    sample_indices = [0, 15, 30, 60, 100, 150, 200, 250, 299]
    print(f"  {'t':>8s}  {'tau(t)':>14s}  {'F(t)':>14s}")
    print(f"  {'-'*8}  {'-'*14}  {'-'*14}")
    for idx in sample_indices:
        print(f"  {t_vals[idx]:8.3f}  {tau_vals[idx]:14.10f}  {fidelity_vals[idx]:14.10f}")

    print()
    print(f"  tau(0)     = {tau_at_0:.2e}  (should be 0)")
    print(f"  tau(80)    = {tau_at_end:.10f}  (should approach 0.5)")
    print(f"  All tau in [0,1]? {in_bounds}")
    print(f"  Monotonicity violations (dtau < -1e-10): {monotone_violations}")
    if monotone_violations > 0:
        worst_idx = np.argmin(diffs)
        print(f"    Worst violation: dtau = {diffs[worst_idx]:.2e} at t = {t_vals[worst_idx]:.3f}")
    else:
        print(f"    Min dtau = {max_violation:.2e} (all non-negative within tolerance)")

    # Analytical check for amplitude damping:
    # rho(t) decays |1><1| -> |0><0| with rate gamma
    # At t -> infinity, rho -> |0><0|
    # sigma = I/2
    # The steady state Petz recovery becomes imperfect
    print("\n  Analytical notes:")
    print("  Amplitude damping: rho(t) = (1-p)|1><1| + p|0><0| + off-diag terms")
    print("  where p = 1 - exp(-gamma*t)")
    print("  At t->inf: rho->|0><0|, channel loses info about initial state")
    print(f"  Expected tau(inf) -> 0.5 (for rho_0=|1><1|, sigma=I/2)")

    # PASS/FAIL
    check_mono = (monotone_violations == 0)
    check_zero = (np.abs(tau_at_0) < 1e-10)
    check_bounds = in_bounds
    passed = check_mono and check_zero and check_bounds

    status = "PASS" if passed else "FAIL"
    print(f"\n  >>> Test 3 criteria:")
    print(f"      tau(0) = 0?         {'PASS' if check_zero else 'FAIL'} (tau(0) = {tau_at_0:.2e})")
    print(f"      Monotone?           {'PASS' if check_mono else 'FAIL'} ({monotone_violations} violations)")
    print(f"      tau in [0,1]?       {'PASS' if check_bounds else 'FAIL'}")
    print(f"  >>> Result: {status}")

    return passed, t_vals, tau_vals


# =============================================================================
# Bonus: Verify superoperator construction
# =============================================================================
def verify_superoperator():
    """Quick sanity checks on the superoperator construction."""
    print("\n" + "=" * 72)
    print("SANITY CHECKS: Superoperator Construction")
    print("=" * 72)

    d = 2

    # Check 1: Unitary channel preserves trace
    H = sz / 2.0
    L_super = build_lindbladian(d, H, [], [])
    t = 1.0
    E_super = expm(L_super * t)
    rho_test = rho_plus.copy()
    rho_out = apply_channel(E_super, rho_test, d)
    tr = np.trace(rho_out).real
    print(f"  Trace preservation (unitary, t=1): Tr(E[rho]) = {tr:.15f}")

    # Check 2: Dissipative channel preserves trace
    L_super_diss = build_lindbladian(d, np.zeros((2, 2), dtype=complex), [sm], [0.1])
    E_super_diss = expm(L_super_diss * t)
    rho_out_diss = apply_channel(E_super_diss, rho_test, d)
    tr_diss = np.trace(rho_out_diss).real
    print(f"  Trace preservation (amplitude damp, t=1): Tr(E[rho]) = {tr_diss:.15f}")

    # Check 3: Steady state of dephasing is diagonal
    L_super_deph = build_lindbladian(d, np.zeros((2, 2), dtype=complex), [sz], [0.1])
    E_super_deph = expm(L_super_deph * 100.0)  # Large t
    rho_steady = apply_channel(E_super_deph, rho_plus, d)
    print(f"  Dephasing steady state (t=100):")
    print(f"    rho = [[{rho_steady[0,0].real:.6f}, {rho_steady[0,1]:.2e}],")
    print(f"           [{rho_steady[1,0]:.2e}, {rho_steady[1,1].real:.6f}]]")
    print(f"    (should be diag(0.5, 0.5))")

    # Check 4: Amplitude damping steady state is |0><0|
    E_super_amp = expm(L_super_diss * 200.0)
    rho_amp_steady = apply_channel(E_super_amp, rho_1, d)
    print(f"  Amplitude damping steady state (t=200):")
    print(f"    rho = [[{rho_amp_steady[0,0].real:.6f}, {rho_amp_steady[0,1]:.2e}],")
    print(f"           [{rho_amp_steady[1,0]:.2e}, {rho_amp_steady[1,1].real:.6f}]]")
    print(f"    (should be diag(1, 0))")

    print()


# =============================================================================
# MAIN
# =============================================================================
if __name__ == "__main__":
    print()
    print("*" * 72)
    print("*  Continuous-time tau(t) Framework Verification")
    print("*  Lindbladian Dynamics with Petz Recovery Map")
    print("*" * 72)
    print()

    # Run sanity checks first
    verify_superoperator()

    # Run the three critical tests
    pass1, t1, tau1_num, tau1_ana = test_pure_dephasing()
    pass2, t2, tau2 = test_pure_unitary()
    pass3, t3, tau3 = test_amplitude_damping()

    # ==========================================================================
    # SUMMARY
    # ==========================================================================
    print("\n" + "=" * 72)
    print("SUMMARY")
    print("=" * 72)
    results = [
        ("Test 1: Pure Dephasing (analytical)", pass1),
        ("Test 2: Pure Unitary (tau=0)", pass2),
        ("Test 3: Amplitude Damping (monotone)", pass3),
    ]
    all_pass = True
    for name, passed in results:
        status = "PASS" if passed else "FAIL"
        marker = "[OK]" if passed else "[XX]"
        print(f"  {marker} {name}: {status}")
        if not passed:
            all_pass = False

    print()
    if all_pass:
        print("  ALL 3 CRITICAL TESTS PASSED.")
        print("  The continuous-time tau(t) framework is numerically verified.")
    else:
        print("  SOME TESTS FAILED. Investigation needed.")

    print()

    # Exit with appropriate code
    sys.exit(0 if all_pass else 1)
