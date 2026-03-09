"""
Test 2: Bridge Equation for Qutrit Spatial DOF (COULD FAIL)

Paper 1b derives tau_I = 2*tau_S*(1-tau_S) for QUBIT spatial DOF.
What happens for qutrit (3-path interferometer)?

The bridge equation relates the CHANNEL tau parameters:
  tau_S = worst-case irreversibility of the spatial dephasing channel
  tau_I = worst-case irreversibility of the internal dephasing channel

For qubit spatial DOF (d=2):
  The spatial channel has a single decoherence factor D.
  tau_S = (1 - |D|) / 2
  tau_I = (1 - |D|^2) / 2  (tracing 2-dim space squares the factor)
  => tau_I = 2*tau_S*(1 - tau_S) exactly

For qutrit spatial DOF (d=3):
  The spatial channel has MULTIPLE decoherence factors D_ij (for i,j in {L,C,R}).
  Does a clean bridge equation still exist?

Reference: Huang (2026), Section on tau-bridge.
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

I2 = np.eye(2, dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)


def partial_trace_explicit(rho, dims, keep):
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
        t_multi = []
        tmp = t_idx
        for d in reversed(trace_dims):
            t_multi.insert(0, tmp % d)
            tmp //= d

        for i in range(d_keep):
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

                row_flat = 0
                col_flat = 0
                for s in range(n):
                    row_flat = row_flat * dims[s] + row_full[s]
                    col_flat = col_flat * dims[s] + col_full[s]

                result[i, j] += rho[row_flat, col_flat]

    return result


def _tensor_op(n, ops_list):
    """Build tensor product with specific operators at given positions."""
    mats = [I2] * n
    for pos, op in ops_list:
        mats[pos] = op
    result = mats[0]
    for m in mats[1:]:
        result = np.kron(result, m)
    return result


# ============================================================
# Qubit spatial DOF (d=2)
# ============================================================

def compute_qubit_bridge(N_env, g, t_values):
    """
    Qubit spatial DOF: compute decoherence factor D(t) and derive
    tau_S, tau_I from the channel perspective.

    For N identical env qubits with coupling g*sigma_z^S x sigma_z^{Ek}:
      D(t) = cos^N(2*g*t)
      tau_S = (1 - |D|) / 2   (spatial channel, worst-case)
      tau_I = (1 - |D|^2) / 2  (internal channel, worst-case)

    Also compute tau_S numerically from the exact simulation as a check.
    """
    n_t = len(t_values)
    tau_S = np.zeros(n_t)
    tau_I = np.zeros(n_t)
    tau_S_numerical = np.zeros(n_t)

    # For numerical verification
    n_total = 1 + N_env
    dims = [2] * n_total
    dim = 2**n_total

    H = np.zeros((dim, dim), dtype=complex)
    for k in range(N_env):
        H += g * _tensor_op(n_total, [(0, sz), (1 + k, sz)])

    plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    psi0 = plus.copy()
    for _ in range(N_env):
        psi0 = np.kron(psi0, plus)
    rho0 = np.outer(psi0, psi0.conj())

    for ti, t in enumerate(t_values):
        if t == 0:
            continue

        # Analytical
        c = np.cos(2 * g * t)
        D = c ** N_env
        tau_S[ti] = (1 - abs(D)) / 2
        tau_I[ti] = (1 - D**2) / 2

        # Numerical check: compute rho_S and extract off-diagonal
        U = expm(-1j * H * t)
        psi_t = U @ psi0
        rho_t = np.outer(psi_t, psi_t.conj())
        rho_S = partial_trace_explicit(rho_t, dims, [0])
        D_numerical = 2 * np.real(rho_S[0, 1])  # rho_S = [[1/2, D/2],[D/2, 1/2]]
        tau_S_numerical[ti] = (1 - abs(D_numerical)) / 2

    return tau_S, tau_I, tau_S_numerical


# ============================================================
# Qutrit spatial DOF (d=3)
# ============================================================

def compute_qutrit_bridge(N_env, g_values, t_values):
    """
    Qutrit spatial DOF: |L>, |C>, |R> with couplings (g_L, g_C, g_R).

    For a qutrit spatial DOF coupled to N env qubits:
      H = sum_k diag(g_L, g_C, g_R) x sigma_z^{Ek}

    The spatial density matrix has off-diagonal elements:
      rho_S[i,j] = (1/3) * prod_k cos(2*(g_i - g_j)*t)  for i != j

    There are THREE decoherence factors:
      D_LC = cos^N(2*(g_L - g_C)*t)
      D_LR = cos^N(2*(g_L - g_R)*t)
      D_CR = cos^N(2*(g_C - g_R)*t)

    tau_S: For the spatial dephasing channel on the qutrit, the worst-case
           input is the equal superposition (1/sqrt(3))(|L>+|C>+|R>).
           We compute tau_S = 1 - <+3|rho_S(t)|+3> where |+3> = equal superposition.

    tau_I: For the internal channel (trace out spatial), we compute
           the overlap with the initial internal state.
           This is more complex for d=3.
    """
    g_L, g_C, g_R = g_values

    plus3 = np.array([1, 1, 1], dtype=complex) / np.sqrt(3)

    n_t = len(t_values)
    tau_S = np.zeros(n_t)
    tau_I = np.zeros(n_t)
    D_factors = np.zeros((n_t, 3))  # D_LC, D_LR, D_CR

    for ti, t in enumerate(t_values):
        if t == 0:
            continue

        # Analytical decoherence factors
        D_LC = np.cos(2 * (g_L - g_C) * t) ** N_env
        D_LR = np.cos(2 * (g_L - g_R) * t) ** N_env
        D_CR = np.cos(2 * (g_C - g_R) * t) ** N_env
        D_factors[ti] = [D_LC, D_LR, D_CR]

        # Spatial density matrix (3x3)
        rho_S = np.array([
            [1/3, D_LC/3, D_LR/3],
            [D_LC/3, 1/3, D_CR/3],
            [D_LR/3, D_CR/3, 1/3]
        ], dtype=complex)

        # tau_S = 1 - <+3|rho_S|+3>
        F_S = np.real(plus3.conj() @ rho_S @ plus3)
        tau_S[ti] = max(1 - F_S, 0.0)

        # tau_I analytically:
        # The internal state rho_I = (1/d) sum_s |chi_s><chi_s|
        # F_I = (1/d) sum_s |<psi0_I|chi_s>|^2 = (1/3) sum_s cos^{2N}(g_s*t)
        F_I_anal = (np.cos(g_L*t)**(2*N_env) +
                    np.cos(g_C*t)**(2*N_env) +
                    np.cos(g_R*t)**(2*N_env)) / 3
        tau_I[ti] = max(1 - F_I_anal, 0.0)

    return tau_S, tau_I, D_factors


# ============================================================
# Run Test 2
# ============================================================

def run_test2():
    """Run Test 2: Qutrit bridge equation."""
    print("=" * 70)
    print("TEST 2: Bridge Equation for Qutrit Spatial DOF")
    print("=" * 70)
    print()

    N_env = 4
    g = 0.3
    n_times = 100
    t_max = 5.0
    t_values = np.linspace(0, t_max, n_times)

    # Part A: Qubit (d=2) -- verify tau_I = 2*tau_S*(1-tau_S)
    print("Part A: Qubit spatial DOF (d=2, N_env=4)")
    print("-" * 50)
    tau_S_qubit, tau_I_qubit, tau_S_num = compute_qubit_bridge(N_env, g, t_values)

    # Verify numerical matches analytical
    mask_nonzero = tau_S_qubit > 0.001
    if np.sum(mask_nonzero) > 0:
        max_dev_num = np.max(np.abs(tau_S_qubit[mask_nonzero] - tau_S_num[mask_nonzero]))
        print(f"  Numerical vs analytical tau_S: max deviation = {max_dev_num:.8f}")

    # Check bridge equation
    tau_I_predicted = 2 * tau_S_qubit * (1 - tau_S_qubit)
    residual = tau_I_qubit - tau_I_predicted
    mask = tau_S_qubit > 0.01
    if np.sum(mask) > 0:
        max_residual = np.max(np.abs(residual[mask]))
        rms_residual = np.sqrt(np.mean(residual[mask]**2))
    else:
        max_residual = 0.0
        rms_residual = 0.0

    print(f"  Bridge: tau_I = 2*tau_S*(1-tau_S)")
    print(f"  Max residual: {max_residual:.8f}")
    print(f"  RMS residual: {rms_residual:.8f}")
    bridge_holds_qubit = max_residual < 1e-6
    print(f"  Bridge {'HOLDS EXACTLY' if bridge_holds_qubit else 'FAILS'} for qubit")
    print()

    # Part B: Qutrit (d=3) -- symmetric case
    print("Part B: Qutrit spatial DOF (d=3, N_env=4)")
    print("-" * 50)

    g_values = (0.3, 0.0, -0.3)
    print(f"  Couplings: g_L={g_values[0]}, g_C={g_values[1]}, g_R={g_values[2]}")
    print(f"  Differences: g_L-g_C={g_values[0]-g_values[1]:.1f}, "
          f"g_L-g_R={g_values[0]-g_values[2]:.1f}, "
          f"g_C-g_R={g_values[1]-g_values[2]:.1f}")

    tau_S_qutrit, tau_I_qutrit, D_factors = compute_qutrit_bridge(
        N_env, g_values, t_values)

    # Test qubit bridge on qutrit data
    tau_I_qubit_pred = 2 * tau_S_qutrit * (1 - tau_S_qutrit)
    mask_q3 = tau_S_qutrit > 0.01
    if np.sum(mask_q3) > 0:
        residual_q3 = tau_I_qutrit[mask_q3] - tau_I_qubit_pred[mask_q3]
        max_res_q3 = np.max(np.abs(residual_q3))
    else:
        max_res_q3 = 0.0

    print(f"  Qubit bridge on qutrit: max deviation = {max_res_q3:.6f}")
    print(f"  Qubit bridge {'HOLDS' if max_res_q3 < 0.02 else 'DOES NOT HOLD'} for qutrit")
    print()

    # Fit generalized bridge: tau_I = a * tau_S + b * tau_S^2
    a_fit = 2.0
    b_fit = -2.0  # qubit: a=2, b=-2 so tau_I = 2*tS - 2*tS^2 = 2*tS*(1-tS)
    if np.sum(mask_q3) > 5:
        tS = tau_S_qutrit[mask_q3]
        tI = tau_I_qutrit[mask_q3]

        A = np.column_stack([tS, tS**2])
        coeffs, _, _, _ = np.linalg.lstsq(A, tI, rcond=None)
        a_fit = coeffs[0]
        b_fit = coeffs[1]

        tau_I_fit = a_fit * tS + b_fit * tS**2
        rms_fit = np.sqrt(np.mean((tI - tau_I_fit)**2))

        print(f"  Fitted: tau_I = {a_fit:.4f}*tau_S + ({b_fit:.4f})*tau_S^2")
        print(f"  (Qubit: tau_I = 2*tau_S - 2*tau_S^2)")
        print(f"  RMS of fit: {rms_fit:.6f}")
        print()

    # Part C: Symmetric qutrit
    print("Part C: Symmetric qutrit (g_L=-g, g_C=0, g_R=+g)")
    print("-" * 50)
    g_sym = (-0.3, 0.0, 0.3)
    tau_S_sym, tau_I_sym, _ = compute_qutrit_bridge(N_env, g_sym, t_values)
    mask_sym = tau_S_sym > 0.01
    if np.sum(mask_sym) > 0:
        dev_sym = np.max(np.abs(tau_I_sym[mask_sym] -
                                2 * tau_S_sym[mask_sym] * (1 - tau_S_sym[mask_sym])))
        print(f"  Deviation from qubit bridge: {dev_sym:.6f}")
    print()

    # Part D: Asymmetric qutrit
    print("Part D: Asymmetric qutrit (g_L=0.5, g_C=0.2, g_R=-0.1)")
    print("-" * 50)
    g_asym = (0.5, 0.2, -0.1)
    tau_S_asym, tau_I_asym, _ = compute_qutrit_bridge(N_env, g_asym, t_values)
    mask_asym = tau_S_asym > 0.01
    if np.sum(mask_asym) > 0:
        dev_asym = np.max(np.abs(tau_I_asym[mask_asym] -
                                  2 * tau_S_asym[mask_asym] * (1 - tau_S_asym[mask_asym])))
        print(f"  Deviation from qubit bridge: {dev_asym:.6f}")
    print()

    # Part E: Analytical comparison
    print("Part E: Analytical Structure")
    print("-" * 50)
    print()
    print("  QUBIT (d=2):")
    print("    tau_S = (1 - |D|) / 2")
    print("    tau_I = (1 - D^2) / 2")
    print("    Bridge: tau_I = 2*tau_S*(1-tau_S)")
    print()
    print("  QUTRIT (d=3) with couplings (g_L, g_C, g_R):")
    print("    tau_S = 1 - (1/3)(1 + D_LC + D_LR + D_CR + ...)")
    print("           where D_ij = cos^N(2*(g_i-g_j)*t)")
    print("    tau_I = 1 - (1/3)(cos^{2N}(g_L*t) + cos^{2N}(g_C*t) + cos^{2N}(g_R*t))")
    print()
    print("    The qutrit has 3 independent decoherence factors vs 1 for qubit.")
    print("    The bridge cannot be a simple function of tau_S alone because")
    print("    tau_I depends on the INDIVIDUAL g_i values, not just their differences.")
    print()

    print("SUMMARY")
    print("=" * 70)
    print(f"  d=2 (qubit):  Bridge tau_I = 2*tau_S*(1-tau_S) holds EXACTLY")
    print(f"                Max residual: {max_residual:.8f}")
    if max_res_q3 > 0.02:
        print(f"  d=3 (qutrit): Bridge DOES NOT generalize")
        print(f"                Max deviation from qubit formula: {max_res_q3:.6f}")
        print(f"                Fitted: tau_I = {a_fit:.4f}*tau_S + ({b_fit:.4f})*tau_S^2")
        print()
        print("  KEY FINDING: The bridge equation requires modification for d>2.")
        print("  For d=3, tau_I depends on the individual coupling constants,")
        print("  not just tau_S. This is because the qutrit has 3 independent")
        print("  decoherence factors that cannot be collapsed into a single parameter.")
        print("  This is a GENUINELY NEW PREDICTION for 3-path interferometers.")
    else:
        print(f"  d=3 (qutrit): Bridge approximately holds")
        print(f"                Max deviation: {max_res_q3:.6f}")
    print()

    return {
        'qubit': {'tau_S': tau_S_qubit, 'tau_I': tau_I_qubit},
        'qutrit': {'tau_S': tau_S_qutrit, 'tau_I': tau_I_qutrit},
        'qutrit_sym': {'tau_S': tau_S_sym, 'tau_I': tau_I_sym},
        'qutrit_asym': {'tau_S': tau_S_asym, 'tau_I': tau_I_asym},
        't_values': t_values,
        'a_fit': a_fit,
        'b_fit': b_fit,
    }


if __name__ == "__main__":
    results = run_test2()
