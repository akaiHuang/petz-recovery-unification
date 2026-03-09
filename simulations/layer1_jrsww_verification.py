"""
Layer 1 Rigorous Verification: JRSWW / Fawzi-Renner Recovery Bound

This is the ONLY layer of the tau framework that is a pure mathematical theorem.

=== PRECISE STATEMENT ===

The Junge-Renner-Sutter-Wilde-Winter (JRSWW) universal recovery theorem
(Ann. Henri Poincare 19, 2955, 2018) guarantees that for any quantum channel N,
reference state sigma, and input state rho, there exists a rotated Petz map
tilde{R} such that:

    F(rho, tilde{R}(N(rho)))^2  >=  exp(-Delta_D)           ... (*)

where F is the ROOT fidelity  F(rho,sigma) = Tr sqrt(sqrt(rho) sigma sqrt(rho)),
and Delta_D = D(rho||sigma) - D(N(rho)||N(sigma)) is the relative entropy decrease.

For sigma = I/d (maximally mixed), the modular automorphisms are trivial,
so the rotated Petz map equals the standard Petz map.  Equivalently, for
sigma = I/d, the standard Petz map already saturates the JRSWW bound.

=== FIDELITY CONVENTIONS (CRITICAL) ===

Two conventions exist in the literature:
  - Root fidelity:  f(rho,sigma) = Tr sqrt(sqrt(rho) sigma sqrt(rho))
  - Jozsa fidelity: F_J(rho,sigma) = [Tr sqrt(sqrt(rho) sigma sqrt(rho))]^2 = f^2

Our code fidelity() function returns the JOZSA convention F_J = f^2.
The paper uses root fidelity f, so the bound (*) reads:

    f^2  >=  exp(-Delta_D)

In code terms:
    fidelity(rho, R(N(rho)))  >=  exp(-Delta_D)

The temporal asymmetry parameter is:
    tau = 1 - f = 1 - sqrt(fidelity(rho, R(N(rho))))

And the upper bound on tau from (*) is:
    tau  <=  1 - exp(-Delta_D / 2)

=== WHAT THIS SCRIPT VERIFIES ===

Part A: fidelity(rho, R(N(rho))) >= exp(-Delta_D) for all (channel, p, rho)
Part B: tau = 0  iff  Delta_D = 0  (equivalence chain)
Part C: Compare tau_Petz, tau_identity, tau_bound for dephasing + |+>
Part D: Output CSV + plot + final PASS/FAIL

Reference: Huang (2026), "Petz Recovery Map as Retrodiction Functor".
"""

import os
import csv
import numpy as np
from scipy.linalg import logm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

_DIR = os.path.dirname(os.path.abspath(__file__))

np.random.seed(42)

# =====================================================================
# Quantum information utilities
# =====================================================================

def ensure_hermitian(M):
    """Force a matrix to be Hermitian."""
    return (M + M.conj().T) / 2


def ensure_density_matrix(rho):
    """Project onto valid density matrix (Hermitian, PSD, trace 1)."""
    rho = ensure_hermitian(rho)
    eigvals, eigvecs = np.linalg.eigh(rho)
    eigvals = np.maximum(eigvals, 0)
    rho = eigvecs @ np.diag(eigvals) @ eigvecs.conj().T
    tr = np.real(np.trace(rho))
    if tr > 1e-15:
        rho = rho / tr
    return rho


def matrix_sqrt_psd(M):
    """Compute PSD matrix square root via eigendecomposition."""
    M_h = ensure_hermitian(M)
    eigvals, eigvecs = np.linalg.eigh(M_h)
    eigvals = np.maximum(eigvals, 0)
    return eigvecs @ np.diag(np.sqrt(eigvals)) @ eigvecs.conj().T


def matrix_inv_sqrt(M, eps=1e-10):
    """Compute M^{-1/2} with regularization for near-zero eigenvalues."""
    M_h = ensure_hermitian(M)
    eigvals, eigvecs = np.linalg.eigh(M_h)
    safe_eigvals = np.where(eigvals > eps, eigvals, 1.0)
    eigvals_inv_sqrt = np.where(eigvals > eps, 1.0 / np.sqrt(safe_eigvals), 0.0)
    return eigvecs @ np.diag(eigvals_inv_sqrt) @ eigvecs.conj().T


def jozsa_fidelity(rho, sigma):
    """Jozsa (squared) fidelity: F_J = [Tr sqrt(sqrt(rho) sigma sqrt(rho))]^2.

    This equals f^2 where f is the root fidelity.
    For pure rho = |psi><psi|, F_J = <psi|sigma|psi>.
    """
    rho = ensure_density_matrix(rho)
    sigma = ensure_density_matrix(sigma)
    sqrt_rho = matrix_sqrt_psd(rho)
    inner = sqrt_rho @ sigma @ sqrt_rho
    inner = ensure_hermitian(inner)
    eigvals = np.linalg.eigvalsh(inner)
    eigvals = np.maximum(eigvals, 0)
    F_J = np.real(np.sum(np.sqrt(eigvals))) ** 2
    return float(min(max(F_J, 0.0), 1.0))


def relative_entropy(rho, sigma, eps=1e-12):
    """D(rho || sigma) = Tr[rho (ln rho - ln sigma)] using natural log.

    Regularizes both states to avoid log(0).
    """
    rho = ensure_density_matrix(rho)
    d = rho.shape[0]
    sigma_reg = ensure_hermitian(sigma) + eps * np.eye(d)
    sigma_reg = sigma_reg / np.real(np.trace(sigma_reg))
    rho_reg = rho + eps * np.eye(d)
    rho_reg = rho_reg / np.real(np.trace(rho_reg))
    log_rho = logm(rho_reg)
    log_sigma = logm(sigma_reg)
    val = np.real(np.trace(rho @ (log_rho - log_sigma)))
    return max(val, 0.0)


# =====================================================================
# Quantum channels (Kraus representation)
# =====================================================================

def dephasing_kraus(p):
    """Dephasing channel: K0 = sqrt((1+p)/2) I, K1 = sqrt((1-p)/2) sigma_z.

    p = 1: identity channel (no noise).
    p = 0: full dephasing (complete decoherence of off-diagonals).
    p in (0,1): partial dephasing.  Off-diagonal elements multiplied by p.
    """
    I = np.eye(2, dtype=complex)
    Z = np.array([[1, 0], [0, -1]], dtype=complex)
    K0 = np.sqrt((1 + p) / 2) * I
    K1 = np.sqrt((1 - p) / 2) * Z
    return [K0, K1]


def amplitude_damping_kraus(gamma):
    """Amplitude damping channel with damping parameter gamma.

    gamma = 0: identity.  gamma = 1: full damping to |0>.
    """
    K0 = np.array([[1.0, 0.0], [0.0, np.sqrt(1.0 - gamma)]], dtype=complex)
    K1 = np.array([[0.0, np.sqrt(gamma)], [0.0, 0.0]], dtype=complex)
    return [K0, K1]


def depolarizing_kraus(p):
    """Depolarizing channel: N(rho) = (1-p)rho + (p/3)(X rho X + Y rho Y + Z rho Z).

    p = 0: identity.  p = 3/4: fully depolarizing (output = I/2).
    p = 1: not fully depolarizing but strong noise.
    Kraus: sqrt(1-p) I, sqrt(p/3) X, sqrt(p/3) Y, sqrt(p/3) Z.
    """
    I = np.eye(2, dtype=complex)
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    Z = np.array([[1, 0], [0, -1]], dtype=complex)
    return [np.sqrt(1.0 - p) * I, np.sqrt(p / 3.0) * X,
            np.sqrt(p / 3.0) * Y, np.sqrt(p / 3.0) * Z]


def apply_channel(rho, kraus_ops):
    """N(rho) = sum_k K_k rho K_k^dag."""
    result = np.zeros_like(rho, dtype=complex)
    for K in kraus_ops:
        result += K @ rho @ K.conj().T
    return ensure_density_matrix(result)


def channel_adjoint(Y, kraus_ops):
    """N^dag(Y) = sum_k K_k^dag Y K_k."""
    result = np.zeros_like(Y, dtype=complex)
    for K in kraus_ops:
        result += K.conj().T @ Y @ K
    return result


# =====================================================================
# Recovery maps
# =====================================================================

def petz_recovery(Y, kraus_ops, sigma):
    """Standard Petz recovery map:
    R_sigma(Y) = sigma^{1/2} N^dag( N(sigma)^{-1/2} Y N(sigma)^{-1/2} ) sigma^{1/2}

    For sigma = I/d (maximally mixed), this equals the rotated Petz map
    since modular automorphisms are trivial.
    """
    N_sigma = apply_channel(sigma, kraus_ops)
    N_sigma_inv_sqrt = matrix_inv_sqrt(N_sigma)
    sigma_sqrt = matrix_sqrt_psd(sigma)
    sandwiched = N_sigma_inv_sqrt @ Y @ N_sigma_inv_sqrt
    adjoint_result = channel_adjoint(sandwiched, kraus_ops)
    result = sigma_sqrt @ adjoint_result @ sigma_sqrt
    return ensure_density_matrix(result)


# =====================================================================
# Input states
# =====================================================================

def standard_states():
    """Return 6 standard qubit states: |0>, |1>, |+>, |->, |i+>, |i->."""
    ket0 = np.array([[1], [0]], dtype=complex)
    ket1 = np.array([[0], [1]], dtype=complex)
    ketp = (ket0 + ket1) / np.sqrt(2)
    ketm = (ket0 - ket1) / np.sqrt(2)
    keti = (ket0 + 1j * ket1) / np.sqrt(2)
    ketmi = (ket0 - 1j * ket1) / np.sqrt(2)

    return {
        '|0>': ket0 @ ket0.conj().T,
        '|1>': ket1 @ ket1.conj().T,
        '|+>': ketp @ ketp.conj().T,
        '|->': ketm @ ketm.conj().T,
        '|i+>': keti @ keti.conj().T,
        '|i->': ketmi @ ketmi.conj().T,
    }


def random_pure_state():
    """Generate a Haar-random pure qubit state."""
    v = np.random.randn(2) + 1j * np.random.randn(2)
    v = v / np.linalg.norm(v)
    v = v.reshape(2, 1)
    return v @ v.conj().T


def random_mixed_state():
    """Generate a random mixed qubit state (inside the Bloch sphere)."""
    r = np.random.randn(3)
    r = r / np.linalg.norm(r) * np.random.uniform(0.1, 0.9)
    I = np.eye(2, dtype=complex)
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    Z = np.array([[1, 0], [0, -1]], dtype=complex)
    rho = (I + r[0] * X + r[1] * Y + r[2] * Z) / 2
    return ensure_density_matrix(rho)


# =====================================================================
# Core verification logic
# =====================================================================

def verify_single_case(channel_name, kraus_fn, p, rho, state_name, sigma):
    """Verify the JRSWW bound for a single (channel, p, rho, sigma) case.

    The bound is: F_J(rho, R(N(rho))) >= exp(-Delta_D)
    where F_J is the Jozsa (squared) fidelity returned by jozsa_fidelity().

    tau = 1 - sqrt(F_J)          (root fidelity convention)
    tau_bound = 1 - exp(-Delta_D/2)

    Returns a dict with all computed quantities.
    """
    kraus_ops = kraus_fn(p)

    # Apply channel
    N_rho = apply_channel(rho, kraus_ops)
    N_sigma = apply_channel(sigma, kraus_ops)

    # Petz recovery of channel output
    R_N_rho = petz_recovery(N_rho, kraus_ops, sigma)

    # Jozsa fidelity F_J = f^2 where f is root fidelity
    F_J_petz = jozsa_fidelity(rho, R_N_rho)

    # Relative entropy decrease (information loss)
    D_before = relative_entropy(rho, sigma)
    D_after = relative_entropy(N_rho, N_sigma)
    Delta_D = max(D_before - D_after, 0.0)

    # The JRSWW bound: F_J >= exp(-Delta_D)
    exp_bound = np.exp(-Delta_D)

    # Root fidelity (paper's F)
    f_petz = np.sqrt(F_J_petz)

    # tau = 1 - f (paper convention)
    tau = 1.0 - f_petz
    tau_bound = 1.0 - np.exp(-Delta_D / 2)

    # Check bound: F_J >= exp(-Delta_D)
    tol = 1e-8
    bound_holds = (F_J_petz >= exp_bound - tol)

    # Identity recovery: "do nothing" baseline
    F_J_identity = jozsa_fidelity(rho, N_rho)
    f_identity = np.sqrt(F_J_identity)
    tau_identity = 1.0 - f_identity

    return {
        'channel': channel_name,
        'p': p,
        'state': state_name,
        'D_before': D_before,
        'D_after': D_after,
        'Delta_D': Delta_D,
        'exp_bound': exp_bound,
        'F_J_petz': F_J_petz,
        'f_petz': f_petz,
        'F_J_identity': F_J_identity,
        'f_identity': f_identity,
        'tau': tau,
        'tau_bound': tau_bound,
        'tau_identity': tau_identity,
        'bound_holds': bound_holds,
        'gap_FJ': F_J_petz - exp_bound,    # F_J - exp(-Delta_D), should be >= 0
        'gap_tau': tau_bound - tau,          # tau_bound - tau, should be >= 0
    }


# =====================================================================
# Main
# =====================================================================

def main():
    print("=" * 78)
    print("  Layer 1 Rigorous Verification: JRSWW Recovery Bound")
    print("  F^2 >= exp(-Delta_D)  where F = root fidelity")
    print("  Equivalently: F_Jozsa >= exp(-Delta_D)")
    print("  Pure mathematical theorem -- no physics assumptions")
    print("=" * 78)

    # Reference state: maximally mixed
    # For sigma = I/2, rotated Petz = standard Petz (trivial modular flow)
    sigma = np.eye(2, dtype=complex) / 2

    # Channels
    channels = {
        'Dephasing': dephasing_kraus,
        'Amplitude_Damping': amplitude_damping_kraus,
        'Depolarizing': depolarizing_kraus,
    }

    # Noise parameters: 99 values from 0.01 to 0.99
    p_values = np.arange(0.01, 1.00, 0.01)

    # Input states: 6 standard + 10 random pure + 4 random mixed = 20 total
    std_states = standard_states()
    random_states = {}
    for i in range(10):
        random_states[f'random_pure_{i}'] = random_pure_state()
    for i in range(4):
        random_states[f'random_mixed_{i}'] = random_mixed_state()

    all_states = {**std_states, **random_states}
    print(f"\n  Input states: {len(all_states)} "
          f"({len(std_states)} standard + {len(random_states)} random)")
    print(f"  Noise parameters: {len(p_values)} values in [0.01, 0.99]")
    print(f"  Channels: {list(channels.keys())}")
    print(f"  Total cases: {len(channels) * len(p_values) * len(all_states)}")

    # ─── Part A: Exhaustive verification ─────────────────────────────
    print("\n" + "=" * 78)
    print("  PART A: Verify F_Jozsa(rho, R(N(rho))) >= exp(-Delta_D) for all cases")
    print("=" * 78)

    results = []
    total_cases = 0
    violations = 0
    violation_details = []

    for ch_name, ch_fn in channels.items():
        ch_count = 0
        ch_violations = 0
        ch_min_gap = np.inf
        ch_max_gap = -np.inf

        for p in p_values:
            for st_name, rho in all_states.items():
                result = verify_single_case(ch_name, ch_fn, p, rho, st_name, sigma)
                results.append(result)
                total_cases += 1
                ch_count += 1

                if not result['bound_holds']:
                    violations += 1
                    ch_violations += 1
                    violation_details.append(result)

                gap = result['gap_FJ']
                ch_min_gap = min(ch_min_gap, gap)
                ch_max_gap = max(ch_max_gap, gap)

        print(f"\n  Channel: {ch_name}")
        print(f"    Cases tested: {ch_count}")
        print(f"    Violations:   {ch_violations}")
        print(f"    F_J gap range: [{ch_min_gap:.6e}, {ch_max_gap:.6e}]")
        if ch_min_gap >= 0:
            print(f"    --> Bound satisfied for all cases (min gap >= 0)")

    print(f"\n  TOTAL cases tested: {total_cases}")
    print(f"  TOTAL violations:   {violations}")

    if violation_details:
        print("\n  VIOLATION DETAILS (first 10):")
        for v in violation_details[:10]:
            print(f"    {v['channel']} p={v['p']:.2f} state={v['state']}: "
                  f"F_J={v['F_J_petz']:.10f}, bound={v['exp_bound']:.10f}, "
                  f"gap={v['gap_FJ']:.2e}")

    # ─── Part B: Equivalence chain tau=0 iff Delta_D=0 ─────────────────
    print("\n" + "=" * 78)
    print("  PART B: Verify equivalence  tau ~ 0  <=>  Delta_D ~ 0")
    print("=" * 78)

    # For the equivalence check, we need consistent thresholds.
    # The relation tau <= 1 - exp(-DD/2) ~ DD/2 for small DD means
    # tau and DD are of comparable order.  We use a single threshold
    # that is small enough to identify genuinely "zero" cases (exact
    # sufficiency) while avoiding false positives from near-invariant states.
    thresh = 1e-8

    tau_zero_but_DD_not = []
    DD_zero_but_tau_not = []
    both_zero = 0
    neither_zero = 0

    for r in results:
        tau_is_zero = (r['tau'] < thresh)
        DD_is_zero = (r['Delta_D'] < thresh)

        if tau_is_zero and DD_is_zero:
            both_zero += 1
        elif tau_is_zero and not DD_is_zero:
            tau_zero_but_DD_not.append(r)
        elif DD_is_zero and not tau_is_zero:
            DD_zero_but_tau_not.append(r)
        else:
            neither_zero += 1

    print(f"\n  Threshold: tau < {thresh}, Delta_D < {thresh}")
    print(f"  Both tau~0 and Delta_D~0:          {both_zero}")
    print(f"  Both tau>0 and Delta_D>0:           {neither_zero}")
    print(f"  tau~0 but Delta_D>0 (VIOLATION):    {len(tau_zero_but_DD_not)}")
    print(f"  Delta_D~0 but tau>0 (VIOLATION):    {len(DD_zero_but_tau_not)}")

    if tau_zero_but_DD_not:
        print("\n  VIOLATIONS (tau~0 but Delta_D>0):")
        for v in tau_zero_but_DD_not[:5]:
            print(f"    {v['channel']} p={v['p']:.2f} state={v['state']}: "
                  f"tau={v['tau']:.2e}, Delta_D={v['Delta_D']:.2e}")

    if DD_zero_but_tau_not:
        print("\n  VIOLATIONS (Delta_D~0 but tau>0):")
        for v in DD_zero_but_tau_not[:5]:
            print(f"    {v['channel']} p={v['p']:.2f} state={v['state']}: "
                  f"tau={v['tau']:.2e}, Delta_D={v['Delta_D']:.2e}")

    equiv_pass = (len(tau_zero_but_DD_not) == 0 and
                  len(DD_zero_but_tau_not) == 0)
    print(f"\n  Equivalence tau=0 <=> Delta_D=0: {'PASS' if equiv_pass else 'FAIL'}")

    # Show which cases have tau ~ 0 (these should be channel-invariant states)
    if both_zero > 0:
        zero_examples = [(r['channel'], r['p'], r['state'])
                         for r in results if r['tau'] < thresh][:8]
        print(f"\n  Examples of tau~0 cases (channel preserves state):")
        for ch, p, st in zero_examples:
            print(f"    {ch} p={p:.2f} state={st}")

    # ─── Part C: Three tau values for dephasing + |+> ─────────────────
    print("\n" + "=" * 78)
    print("  PART C: Compare three tau values (dephasing channel, |+> input)")
    print("  tau = 1 - sqrt(F_Jozsa)  [root fidelity convention]")
    print("=" * 78)

    ket0 = np.array([[1], [0]], dtype=complex)
    ket1 = np.array([[0], [1]], dtype=complex)
    ketp = (ket0 + ket1) / np.sqrt(2)
    rho_plus = ketp @ ketp.conj().T

    p_fine = np.linspace(0.01, 0.99, 200)
    tau_petz_arr = []
    tau_identity_arr = []
    tau_bound_arr = []
    FJ_petz_arr = []
    FJ_bound_arr = []
    ordering_violations_bound = 0

    for p in p_fine:
        kraus_ops = dephasing_kraus(p)
        N_rho = apply_channel(rho_plus, kraus_ops)
        N_sigma = apply_channel(sigma, kraus_ops)
        R_N_rho = petz_recovery(N_rho, kraus_ops, sigma)

        F_J_petz = jozsa_fidelity(rho_plus, R_N_rho)
        F_J_identity = jozsa_fidelity(rho_plus, N_rho)

        D_before = relative_entropy(rho_plus, sigma)
        D_after = relative_entropy(N_rho, N_sigma)
        Delta_D = max(D_before - D_after, 0.0)

        # Root fidelities
        f_petz = np.sqrt(F_J_petz)
        f_identity = np.sqrt(F_J_identity)

        # tau values
        tau_p = 1.0 - f_petz
        tau_id = 1.0 - f_identity
        tau_b = 1.0 - np.exp(-Delta_D / 2)

        tau_petz_arr.append(tau_p)
        tau_identity_arr.append(tau_id)
        tau_bound_arr.append(tau_b)
        FJ_petz_arr.append(F_J_petz)
        FJ_bound_arr.append(np.exp(-Delta_D))

        # The JRSWW bound requires tau_Petz <= tau_bound
        tol = 1e-7
        if tau_p > tau_b + tol:
            ordering_violations_bound += 1

    tau_petz_arr = np.array(tau_petz_arr)
    tau_identity_arr = np.array(tau_identity_arr)
    tau_bound_arr = np.array(tau_bound_arr)
    FJ_petz_arr = np.array(FJ_petz_arr)
    FJ_bound_arr = np.array(FJ_bound_arr)

    # Report values at key noise levels
    for label, idx in [("p=0.10 (weak noise)", 9),
                       ("p=0.50 (mid noise)", len(p_fine)//2),
                       ("p=0.90 (strong noise)", 182)]:
        if idx < len(p_fine):
            print(f"\n  At {label} (dephasing param p={p_fine[idx]:.2f}):")
            print(f"    tau_bound    = {tau_bound_arr[idx]:.6f}  "
                  f"(upper bound from JRSWW)")
            print(f"    tau_Petz     = {tau_petz_arr[idx]:.6f}  "
                  f"(actual Petz recovery)")
            print(f"    tau_identity = {tau_identity_arr[idx]:.6f}  "
                  f"(do nothing)")

    # Ordering analysis
    petz_le_bound = np.all(tau_petz_arr <= tau_bound_arr + 1e-7)
    petz_le_identity = np.all(tau_petz_arr <= tau_identity_arr + 1e-7)
    identity_le_bound = np.all(tau_identity_arr <= tau_bound_arr + 1e-7)

    print(f"\n  Ordering analysis across all p:")
    print(f"    tau_Petz <= tau_bound (required by JRSWW):  "
          f"{'YES' if petz_le_bound else 'NO'}  "
          f"({ordering_violations_bound} violations)")
    print(f"    tau_Petz <= tau_identity (Petz is better):  "
          f"{'YES' if petz_le_identity else 'NO'}")
    print(f"    tau_identity <= tau_bound:                   "
          f"{'YES' if identity_le_bound else 'NO'}")

    # Note about dephasing: for sigma=I/2, Petz = N^dag = N (self-adjoint).
    # So R(N(rho)) = N(N(rho)) = N^2(rho), which is WORSE than N(rho).
    # Thus tau_Petz > tau_identity for dephasing! This is expected because
    # the Petz map optimizes recovery from the JRSWW perspective (F^2 >= exp(-DD)),
    # not in the sense of maximizing raw fidelity over all maps.
    print(f"\n  NOTE: For dephasing with sigma=I/2, the Petz map equals N^dag = N")
    print(f"  (self-adjoint channel), so R(N(rho)) = N^2(rho).")
    print(f"  This is worse than identity recovery N(rho).")
    print(f"  Nevertheless, tau_Petz <= tau_bound always holds (theorem).")

    # ─── Plot ─────────────────────────────────────────────────────────
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 10))

    # Panel 1: tau comparison
    ax1.plot(p_fine, tau_bound_arr, 'r-', linewidth=2.0,
             label=r'$\tau_{\mathrm{bound}} = 1 - e^{-\Delta D/2}$ (JRSWW)')
    ax1.plot(p_fine, tau_petz_arr, 'b-', linewidth=2.0,
             label=r'$\tau_{\mathrm{Petz}} = 1 - \sqrt{F_J(\rho, R(\mathcal{N}(\rho)))}$')
    ax1.plot(p_fine, tau_identity_arr, 'g--', linewidth=1.5,
             label=r'$\tau_{\mathrm{id}} = 1 - \sqrt{F_J(\rho, \mathcal{N}(\rho))}$')

    # Shade the gap between Petz and bound
    ax1.fill_between(p_fine, tau_petz_arr, tau_bound_arr,
                     where=(tau_bound_arr >= tau_petz_arr),
                     alpha=0.15, color='blue',
                     label='Margin (bound - actual)')

    ax1.set_xlabel('Dephasing parameter $p$\n'
                   r'($p \to 1$: identity, $p \to 0$: full dephasing)',
                   fontsize=11)
    ax1.set_ylabel(r'$\tau = 1 - F$  (root fidelity)', fontsize=12)
    ax1.set_title(r'Layer 1: Three $\tau$ values for dephasing channel'
                  '\n'
                  r'Input: $|+\rangle$, Reference: $\sigma = I/2$',
                  fontsize=13)
    ax1.legend(fontsize=9, loc='upper left')
    ax1.set_xlim(0, 1)
    ax1.set_ylim(-0.02, 0.55)
    ax1.grid(True, alpha=0.3)

    ax1.annotate(r'$\tau_{\mathrm{Petz}} \leq \tau_{\mathrm{bound}}$ always'
                 '\n(JRSWW theorem)',
                 xy=(0.5, 0.42), fontsize=10, color='darkred',
                 ha='center',
                 bbox=dict(boxstyle='round,pad=0.3', facecolor='lightyellow',
                           alpha=0.8))

    # Panel 2: F_Jozsa vs exp(-Delta_D)
    ax2.plot(p_fine, FJ_petz_arr, 'b-', linewidth=2.0,
             label=r'$F_J(\rho, R(\mathcal{N}(\rho)))$')
    ax2.plot(p_fine, FJ_bound_arr, 'r--', linewidth=2.0,
             label=r'$\exp(-\Delta D)$ (lower bound)')

    ax2.fill_between(p_fine, FJ_bound_arr, FJ_petz_arr,
                     where=(FJ_petz_arr >= FJ_bound_arr),
                     alpha=0.15, color='green',
                     label='Gap (verified $\\geq 0$)')

    ax2.set_xlabel('Dephasing parameter $p$', fontsize=11)
    ax2.set_ylabel(r'$F_J = f^2$ (Jozsa fidelity)', fontsize=12)
    ax2.set_title(r'Verification: $F_J \geq \exp(-\Delta D)$', fontsize=13)
    ax2.legend(fontsize=10, loc='lower right')
    ax2.set_xlim(0, 1)
    ax2.set_ylim(0, 1.05)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plot_path = os.path.join(_DIR, 'layer1_tau_comparison.png')
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.savefig(os.path.join(_DIR, 'layer1_tau_comparison.pdf'),
                bbox_inches='tight')
    print(f"\n  Plot saved to: {plot_path}")

    # ─── Part D: Save CSV and final verdict ───────────────────────────
    print("\n" + "=" * 78)
    print("  PART D: Output summary and CSV")
    print("=" * 78)

    csv_path = os.path.join(_DIR, 'layer1_results.csv')
    fieldnames = ['channel', 'p', 'state', 'D_before', 'D_after', 'Delta_D',
                  'exp_bound', 'F_J_petz', 'f_petz', 'F_J_identity', 'f_identity',
                  'tau', 'tau_bound', 'tau_identity', 'bound_holds',
                  'gap_FJ', 'gap_tau']
    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in results:
            writer.writerow(r)
    print(f"\n  CSV saved to: {csv_path}")
    print(f"  Total rows: {len(results)}")

    # Summary statistics
    print("\n  " + "-" * 74)
    print("  SUMMARY STATISTICS BY CHANNEL")
    print("  " + "-" * 74)

    for ch_name in channels:
        ch_results = [r for r in results if r['channel'] == ch_name]
        gaps_FJ = [r['gap_FJ'] for r in ch_results]
        gaps_tau = [r['gap_tau'] for r in ch_results]
        deltas = [r['Delta_D'] for r in ch_results]
        taus = [r['tau'] for r in ch_results]

        print(f"\n  {ch_name} ({len(ch_results)} cases):")
        print(f"    Delta_D range:   [{min(deltas):.6f}, {max(deltas):.6f}]")
        print(f"    tau range:       [{min(taus):.6f}, {max(taus):.6f}]")
        print(f"    F_J gap (F_J - exp(-DD)):  "
              f"[{min(gaps_FJ):.6e}, {max(gaps_FJ):.6e}]")
        print(f"    tau gap (bound - actual):  "
              f"[{min(gaps_tau):.6e}, {max(gaps_tau):.6e}]")

    # ─── Final verdict ────────────────────────────────────────────────
    print("\n" + "=" * 78)
    partA_pass = (violations == 0)
    partB_pass = equiv_pass
    partC_pass = petz_le_bound

    print(f"  Part A  F_J >= exp(-Delta_D):          "
          f"{'PASS' if partA_pass else 'FAIL'}  "
          f"[{total_cases} cases, {violations} violations]")
    print(f"  Part B  tau=0 <=> Delta_D=0:           "
          f"{'PASS' if partB_pass else 'FAIL'}")
    print(f"  Part C  tau_Petz <= tau_bound:          "
          f"{'PASS' if partC_pass else 'FAIL'}  "
          f"[{ordering_violations_bound} ordering violations]")

    overall = partA_pass and partB_pass and partC_pass
    print()
    print(f"  {'=' * 50}")
    if overall:
        print(f"  OVERALL LAYER 1 VERIFICATION:   PASS")
    else:
        print(f"  OVERALL LAYER 1 VERIFICATION:   FAIL")
    print(f"  {'=' * 50}")

    if overall:
        print(f"""
  The JRSWW bound  F^2 >= exp(-Delta_D)  is verified for ALL {total_cases}
  tested cases across three qubit channels ({', '.join(channels.keys())}),
  {len(all_states)} input states per parameter, and {len(p_values)} noise parameters.

  Convention: F = root fidelity = Tr sqrt(sqrt(rho) sigma sqrt(rho)).
  Our fidelity function returns F^2 (Jozsa convention), so the bound
  verified in code is: jozsa_fidelity(rho, R(N(rho))) >= exp(-Delta_D).

  The equivalence tau = 0 <=> Delta_D = 0 holds without exception.
  This is a pure mathematical theorem -- no physics assumptions needed.

  For sigma = I/2, the standard Petz map equals the rotated Petz map
  (trivial modular automorphisms), so the JRSWW bound applies directly.""")
    else:
        print("\n  WARNING: Verification FAILED. Check details above.")

    print("=" * 78)

    return overall


if __name__ == '__main__':
    success = main()
    raise SystemExit(0 if success else 1)
