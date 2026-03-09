"""
Test 1: tau-R Correlation with Correlated Environment (COULD FAIL)

The tau-R bridge tau = 1-exp(-N*Gamma*t^2), R = N[1-exp(-Gamma*t^2)] assumes
INDEPENDENT environment modes. But real environments have correlations.

This test introduces environment-environment coupling J and checks whether
the universal tau vs R/N curve survives.

Simulation:
  System: 1 qubit
  Environment: N = 6 qubits with CORRELATED coupling
  Hamiltonian: H = sum_k g_k sigma_z^S x sigma_z^{E_k}
                 + J sum_{k,k+1} sigma_z^{E_k} sigma_z^{E_{k+1}}
  Vary J from 0 (independent) to J >> g (strongly correlated)

What could happen:
  J = 0: tau-R curve matches our prediction (baseline)
  J > 0: tau-R curve DEVIATES -> maps the domain of validity
  Strong J: tau large but R small (irreversibility WITHOUT objectivity)

Reference: Huang (2026), tau-R bridge.
"""

import os
import numpy as np
from scipy.linalg import expm
import warnings
warnings.filterwarnings('ignore')

_DIR = os.path.dirname(os.path.abspath(__file__))


# ============================================================
# Quantum primitives
# ============================================================

# Pauli matrices
I2 = np.eye(2, dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)


def tensor_list(ops):
    """Tensor product of a list of 2x2 operators."""
    result = ops[0]
    for op in ops[1:]:
        result = np.kron(result, op)
    return result


def partial_trace(rho, dims, keep):
    """
    Partial trace over all subsystems except 'keep'.
    dims: list of subsystem dimensions
    keep: list of subsystem indices to keep
    """
    n = len(dims)
    total_dim = int(np.prod(dims))
    rho_reshaped = rho.reshape([dims[i] for i in range(n)] * 2)

    # Determine which indices to trace over
    trace_indices = [i for i in range(n) if i not in keep]

    # Trace over unwanted subsystems (from highest index to lowest)
    for idx in sorted(trace_indices, reverse=True):
        rho_reshaped = np.trace(rho_reshaped, axis1=idx, axis2=idx + n - (n - len(dims)))

    # This approach is complex for general case, use explicit method instead
    return _partial_trace_explicit(rho, dims, keep)


def _partial_trace_explicit(rho, dims, keep):
    """Explicit partial trace by iterating over basis states."""
    n = len(dims)
    keep = sorted(keep)
    trace_out = sorted([i for i in range(n) if i not in keep])

    keep_dims = [dims[i] for i in keep]
    trace_dims = [dims[i] for i in trace_out]

    d_keep = int(np.prod(keep_dims))
    d_trace = int(np.prod(trace_dims))

    result = np.zeros((d_keep, d_keep), dtype=complex)

    for t_idx in range(d_trace):
        # Convert t_idx to multi-index for traced-out subsystems
        t_multi = []
        tmp = t_idx
        for d in reversed(trace_dims):
            t_multi.insert(0, tmp % d)
            tmp //= d

        for i in range(d_keep):
            # Convert i to multi-index for kept subsystems
            i_multi = []
            tmp = i
            for d in reversed(keep_dims):
                i_multi.insert(0, tmp % d)
                tmp //= d

            for j in range(d_keep):
                j_multi = []
                tmp = j
                for d in reversed(keep_dims):
                    j_multi.insert(0, tmp % d)
                    tmp //= d

                # Build full indices
                row_full = [0] * n
                col_full = [0] * n
                ki = 0
                ti = 0
                for s in range(n):
                    if s in keep:
                        row_full[s] = i_multi[ki]
                        col_full[s] = j_multi[ki]
                        ki += 1
                    else:
                        row_full[s] = t_multi[ti]
                        col_full[s] = t_multi[ti]
                        ti += 1

                # Convert multi-index to flat index
                row_flat = 0
                col_flat = 0
                for s in range(n):
                    row_flat = row_flat * dims[s] + row_full[s]
                    col_flat = col_flat * dims[s] + col_full[s]

                result[i, j] += rho[row_flat, col_flat]

    return result


def von_neumann_entropy(rho):
    """S(rho) = -Tr[rho log rho]."""
    rho_h = (rho + rho.conj().T) / 2
    eigvals = np.linalg.eigvalsh(rho_h)
    eigvals = eigvals[eigvals > 1e-15]
    return -np.sum(eigvals * np.log(eigvals))


def fidelity_dm(rho, sigma):
    """Trace fidelity F = Tr(sqrt(sqrt(rho) sigma sqrt(rho)))."""
    from scipy.linalg import sqrtm
    sqrt_rho = sqrtm((rho + rho.conj().T) / 2)
    M = sqrt_rho @ ((sigma + sigma.conj().T) / 2) @ sqrt_rho
    M = (M + M.conj().T) / 2
    eigvals = np.linalg.eigvalsh(M)
    return np.real(np.sum(np.sqrt(np.maximum(eigvals, 0))))


# ============================================================
# Build Hamiltonian
# ============================================================

def build_hamiltonian(N_env, g, J):
    """
    Build the full system+environment Hamiltonian.

    H = sum_k g * sigma_z^S x sigma_z^{E_k}
      + J * sum_{k,k+1} sigma_z^{E_k} x sigma_z^{E_{k+1}}

    Total Hilbert space: 2^(1+N_env)
    Subsystem 0 = system, subsystems 1..N_env = environment
    """
    n_total = 1 + N_env
    dim = 2**n_total
    H = np.zeros((dim, dim), dtype=complex)

    # System-environment coupling: g * sigma_z^S x sigma_z^{E_k}
    for k in range(N_env):
        ops = [I2] * n_total
        ops[0] = sz  # system
        ops[1 + k] = sz  # env qubit k
        H += g * tensor_list(ops)

    # Environment-environment coupling: J * sigma_z^{E_k} x sigma_z^{E_{k+1}}
    for k in range(N_env - 1):
        ops = [I2] * n_total
        ops[1 + k] = sz
        ops[1 + k + 1] = sz
        H += J * tensor_list(ops)

    return H


# ============================================================
# Simulation
# ============================================================

def simulate_correlated_env(N_env, g, J, t_values):
    """
    Simulate dephasing with correlated environment.

    Initial state: |+>_S x |+>^{xN}_E  (all in |+> state)

    Returns: tau_vals, R_norm_vals, MI_per_mode
    """
    n_total = 1 + N_env
    dim = 2**n_total
    dims = [2] * n_total

    # Build Hamiltonian
    H = build_hamiltonian(N_env, g, J)

    # Initial state: |+>^{x(1+N_env)}
    plus_state = np.array([1, 1], dtype=complex) / np.sqrt(2)
    psi0 = plus_state.copy()
    for _ in range(N_env):
        psi0 = np.kron(psi0, plus_state)

    rho0_S = _partial_trace_explicit(np.outer(psi0, psi0.conj()), dims, [0])

    n_t = len(t_values)
    tau_vals = np.zeros(n_t)
    R_norm_vals = np.zeros(n_t)
    MI_per_mode_vals = np.zeros(n_t)

    for ti, t in enumerate(t_values):
        if t == 0:
            tau_vals[ti] = 0.0
            R_norm_vals[ti] = 0.0
            MI_per_mode_vals[ti] = 0.0
            continue

        # Time evolution
        U = expm(-1j * H * t)
        psi_t = U @ psi0
        rho_t = np.outer(psi_t, psi_t.conj())

        # System state
        rho_S_t = _partial_trace_explicit(rho_t, dims, [0])

        # tau = 1 - F(rho_S(0), rho_S(t))
        F = fidelity_dm(rho0_S, rho_S_t)
        tau_vals[ti] = max(1.0 - F, 0.0)

        # Mutual information and redundancy
        S_sys = von_neumann_entropy(rho_S_t)

        # Average mutual information per mode
        MI_total = 0.0
        for k in range(N_env):
            # I(S:E_k) = S(S) + S(E_k) - S(SE_k)
            rho_Ek = _partial_trace_explicit(rho_t, dims, [1 + k])
            rho_SEk = _partial_trace_explicit(rho_t, dims, [0, 1 + k])
            S_Ek = von_neumann_entropy(rho_Ek)
            S_SEk = von_neumann_entropy(rho_SEk)
            MI_k = max(S_sys + S_Ek - S_SEk, 0.0)
            MI_total += MI_k

        MI_avg = MI_total / N_env if N_env > 0 else 0.0
        MI_per_mode_vals[ti] = MI_avg

        # Redundancy: R_delta/N
        delta = 0.05
        S_threshold = 0.05 * np.log(2)
        if S_sys > S_threshold and MI_avg > 1e-15:
            R_norm_vals[ti] = min(MI_avg / ((1 - delta) * S_sys), 1.0)
        else:
            R_norm_vals[ti] = 0.0

    return tau_vals, R_norm_vals, MI_per_mode_vals


def run_test1():
    """Run Test 1: Correlated environment."""
    print("=" * 70)
    print("TEST 1: tau-R Correlation with Correlated Environment")
    print("=" * 70)
    print()

    N_env = 6
    g = 0.3
    n_times = 80
    t_max = 4.0
    t_values = np.linspace(0, t_max, n_times)

    # J values: from 0 (independent) to strongly correlated
    J_values = [0.0, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0]

    results = {}

    for J in J_values:
        print(f"  Simulating J = {J:.2f} (N_env = {N_env})...")
        tau, R_norm, MI = simulate_correlated_env(N_env, g, J, t_values)
        results[J] = {
            'tau': tau,
            'R_norm': R_norm,
            'MI': MI,
            't': t_values,
        }

        # Summary statistics
        tau_max = np.max(tau)
        R_max = np.max(R_norm)
        # Check if tau and R/N track each other
        # Use correlation coefficient on the rising part
        mask = tau > 0.01
        if np.sum(mask) > 3:
            corr = np.corrcoef(tau[mask], R_norm[mask])[0, 1]
        else:
            corr = np.nan
        print(f"    tau_max = {tau_max:.4f}, R/N_max = {R_max:.4f}, "
              f"corr(tau, R/N) = {corr:.4f}")

    # Analysis: deviation from J=0 baseline
    print()
    print("ANALYSIS: Deviation from independent-mode baseline (J=0)")
    print("-" * 60)
    print(f"{'J':>6}  {'max|tau-tau_0|':>15}  {'max|R-R_0|':>15}  "
          f"{'corr(tau,R/N)':>15}")
    print("-" * 60)

    tau_0 = results[0.0]['tau']
    R_0 = results[0.0]['R_norm']

    for J in J_values:
        tau_J = results[J]['tau']
        R_J = results[J]['R_norm']
        delta_tau = np.max(np.abs(tau_J - tau_0))
        delta_R = np.max(np.abs(R_J - R_0))
        mask = tau_J > 0.01
        if np.sum(mask) > 3:
            corr = np.corrcoef(tau_J[mask], R_J[mask])[0, 1]
        else:
            corr = np.nan
        print(f"{J:>6.2f}  {delta_tau:>15.6f}  {delta_R:>15.6f}  "
              f"{corr:>15.4f}")

    print()
    print("INTERPRETATION:")
    print("  If corr(tau, R/N) remains high for all J:")
    print("    -> tau-R bridge is ROBUST to env correlations")
    print("  If corr drops or tau >> R/N at large J:")
    print("    -> Irreversibility WITHOUT objectivity (genuinely new regime)")
    print()

    return results, t_values, J_values


if __name__ == "__main__":
    results, t_values, J_values = run_test1()
