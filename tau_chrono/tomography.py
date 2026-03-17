"""
Process tomography: extract Kraus operators from measurement data.

Works with any Qiskit-compatible backend (Quantum Inspire, IBM, IonQ, etc.)
by running standard input-state / measurement-basis circuits and
reconstructing the Choi matrix via linear inversion.

Public API:
    - ``build_tomography_circuits_1q``: generate 1-qubit tomography circuits
    - ``build_tomography_circuits_2q``: generate 2-qubit tomography circuits
    - ``counts_to_choi_1q``: reconstruct 4x4 Choi matrix from 12 measurements
    - ``counts_to_choi_2q``: reconstruct 16x16 Choi matrix from 144 measurements
    - ``choi_to_kraus``: eigendecompose Choi matrix into Kraus operators (CPTP-projected)
    - ``run_process_tomography_1q``: end-to-end 1Q pipeline (needs a Qiskit backend)
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray
from typing import List, Dict, Tuple, Optional, Callable

KrausList = List[NDArray[np.complexfloating]]


# ---------------------------------------------------------------------------
# Standard bases for 1-qubit QPT
# ---------------------------------------------------------------------------

# Input states: |0>, |1>, |+>, |+i>
# These span the 1-qubit state space (informationally complete)
INPUT_LABELS_1Q = ["0", "1", "+", "+i"]
INPUT_STATES_1Q = [
    np.array([[1, 0], [0, 0]], dtype=complex),            # |0><0|
    np.array([[0, 0], [0, 1]], dtype=complex),            # |1><1|
    np.array([[0.5, 0.5], [0.5, 0.5]], dtype=complex),    # |+><+|
    np.array([[0.5, -0.5j], [0.5j, 0.5]], dtype=complex), # |+i><+i|
]

# Measurement operators: Z eigenbasis is default; X and Y via rotations
MEAS_LABELS = ["Z", "X", "Y"]


# ---------------------------------------------------------------------------
# 2-qubit input states (tensor products of 1Q states)
# ---------------------------------------------------------------------------

INPUT_LABELS_2Q = [
    f"{a}{b}" for a in INPUT_LABELS_1Q for b in INPUT_LABELS_1Q
]
INPUT_STATES_2Q = [
    np.kron(s0, s1)
    for s0 in INPUT_STATES_1Q
    for s1 in INPUT_STATES_1Q
]


def _prep_circuit_1q(label: str) -> List[str]:
    """Return gate sequence to prepare the given 1-qubit input state."""
    if label == "0":
        return []
    elif label == "1":
        return ["x"]
    elif label == "+":
        return ["h"]
    elif label == "+i":
        return ["h", "s"]
    raise ValueError(f"Unknown input label: {label}")


def _meas_circuit_1q(basis: str) -> List[str]:
    """Return gate sequence before measurement for the given basis."""
    if basis == "Z":
        return []
    elif basis == "X":
        return ["h"]
    elif basis == "Y":
        return ["sdg", "h"]
    raise ValueError(f"Unknown basis: {basis}")


# ---------------------------------------------------------------------------
# Build Qiskit circuits for process tomography
# ---------------------------------------------------------------------------

def build_tomography_circuits_1q(
    gate_circuit_fn: Callable,
    shots: int = 4096,
) -> Tuple[list, List[Tuple[str, str]]]:
    """Build all process tomography circuits for a 1-qubit channel.

    Parameters
    ----------
    gate_circuit_fn : callable
        A function that takes a Qiskit QuantumCircuit and qubit index,
        and appends the gate(s) under test.
        Example: lambda qc, q: qc.h(q)  # test a Hadamard gate
    shots : int
        Number of shots per circuit.

    Returns
    -------
    circuits : list of QuantumCircuit
        The 12 tomography circuits (4 inputs x 3 bases).
    labels : list of (input_label, meas_label) tuples
        Corresponding labels for each circuit.
    """
    try:
        from qiskit import QuantumCircuit
    except ImportError:
        raise ImportError(
            "Qiskit is required for circuit construction. "
            "Install with: pip install qiskit"
        )

    circuits = []
    labels = []

    for inp in INPUT_LABELS_1Q:
        for basis in MEAS_LABELS:
            qc = QuantumCircuit(1, 1, name=f"tomo_{inp}_{basis}")

            # Prepare input state
            for gate in _prep_circuit_1q(inp):
                getattr(qc, gate)(0)

            qc.barrier()

            # Apply the gate under test
            gate_circuit_fn(qc, 0)

            qc.barrier()

            # Rotate to measurement basis
            for gate in _meas_circuit_1q(basis):
                getattr(qc, gate)(0)

            # Measure
            qc.measure(0, 0)

            circuits.append(qc)
            labels.append((inp, basis))

    return circuits, labels


def build_tomography_circuits_2q(
    gate_circuit_fn: Callable,
    shots: int = 4096,
) -> Tuple[list, List[Tuple[str, str, str]]]:
    """Build process tomography circuits for a 2-qubit channel.

    Parameters
    ----------
    gate_circuit_fn : callable
        A function that takes (qc, q0, q1) and appends the gate(s).
        Example: lambda qc, q0, q1: qc.cx(q0, q1)
    shots : int
        Number of shots per circuit.

    Returns
    -------
    circuits : list of QuantumCircuit
        The 144 tomography circuits (16 inputs x 9 bases).
    labels : list of (input_label_0, input_label_1, meas_basis_0, meas_basis_1)
    """
    try:
        from qiskit import QuantumCircuit
    except ImportError:
        raise ImportError("Qiskit required. Install: pip install qiskit")

    circuits = []
    labels = []

    for inp0 in INPUT_LABELS_1Q:
        for inp1 in INPUT_LABELS_1Q:
            for b0 in MEAS_LABELS:
                for b1 in MEAS_LABELS:
                    qc = QuantumCircuit(2, 2,
                                        name=f"tomo_{inp0}{inp1}_{b0}{b1}")

                    # Prepare input states
                    for gate in _prep_circuit_1q(inp0):
                        getattr(qc, gate)(0)
                    for gate in _prep_circuit_1q(inp1):
                        getattr(qc, gate)(1)

                    qc.barrier()
                    gate_circuit_fn(qc, 0, 1)
                    qc.barrier()

                    # Measurement rotations
                    for gate in _meas_circuit_1q(b0):
                        getattr(qc, gate)(0)
                    for gate in _meas_circuit_1q(b1):
                        getattr(qc, gate)(1)

                    qc.measure([0, 1], [0, 1])
                    circuits.append(qc)
                    labels.append((inp0, inp1, b0, b1))

    return circuits, labels


# ---------------------------------------------------------------------------
# Reconstruct Choi matrix from measurement counts
# ---------------------------------------------------------------------------

# Pauli matrices for reconstruction
_I = np.eye(2, dtype=complex)
_X = np.array([[0, 1], [1, 0]], dtype=complex)
_Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
_Z = np.array([[1, 0], [0, -1]], dtype=complex)

_PAULIS = [_I, _X, _Y, _Z]
_PAULI_LABELS = ["I", "X", "Y", "Z"]


def _expectation_from_counts(counts: Dict[str, int], qubit: int = 0) -> float:
    """Extract <Z> expectation value from measurement counts for a qubit."""
    total = sum(counts.values())
    p0, p1 = 0.0, 0.0
    for bitstring, count in counts.items():
        # Qiskit uses little-endian: rightmost bit is qubit 0
        bit = int(bitstring[-(qubit + 1)])
        if bit == 0:
            p0 += count
        else:
            p1 += count
    return (p0 - p1) / total


def counts_to_choi_1q(
    counts_list: List[Dict[str, int]],
    labels: List[Tuple[str, str]],
) -> NDArray:
    """Reconstruct the 4x4 Choi matrix for a 1-qubit channel.

    Direct reconstruction via the Choi-Jamiolkowski isomorphism:
    J = sum_{k,l} |k><l| tensor N(|k><l|)

    The output density matrices N(rho_input) are reconstructed from
    Pauli expectation values, then N(|k><l|) are extracted by
    linearity from the four informationally complete input states.

    Parameters
    ----------
    counts_list : list of dict
        Measurement counts for each tomography circuit.
    labels : list of (input_label, meas_label) tuples
        Corresponding labels.

    Returns
    -------
    NDArray
        4x4 Choi matrix (trace = d = 2, TP condition: Tr_output = I_input).
    """
    # Step 1: Build expectation value table
    exp_table = {}
    for counts, (inp, basis) in zip(counts_list, labels):
        exp_table[(inp, basis)] = _expectation_from_counts(counts, qubit=0)

    # Step 2: Reconstruct output density matrix for each input
    # rho_out = (I + <X>X + <Y>Y + <Z>Z) / 2
    output_states = {}
    for inp in INPUT_LABELS_1Q:
        ex = exp_table.get((inp, "X"), 0.0)
        ey = exp_table.get((inp, "Y"), 0.0)
        ez = exp_table.get((inp, "Z"), 0.0)
        output_states[inp] = (_I + ex * _X + ey * _Y + ez * _Z) / 2

    # Step 3: Extract N(|k><l|) by linearity
    # |0><0| = (I+Z)/2, |1><1| = (I-Z)/2
    # |+><+| = (|0>+|1>)(<0|+<1|)/2 = (|0><0|+|0><1|+|1><0|+|1><1|)/2
    # |+i><+i| = (|0>+i|1>)(<0|-i<1|)/2 = (|0><0|-i|0><1|+i|1><0|+|1><1|)/2
    N00 = output_states["0"]   # N(|0><0|)
    N11 = output_states["1"]   # N(|1><1|)
    Npp = output_states["+"]   # N(|+><+|)
    Npi = output_states["+i"]  # N(|+i><+i|)

    # From linearity:
    # Npp = (N00 + N01 + N10 + N11) / 2
    # Npi = (N00 - iN01 + iN10 + N11) / 2
    # So: N01 + N10 = 2*Npp - N00 - N11  ... (alpha)
    #     -iN01 + iN10 = 2*Npi - N00 - N11  ... (beta)
    #     => N01 = (alpha + i*beta) / 2
    #     => N10 = (alpha - i*beta) / 2
    alpha = 2 * Npp - N00 - N11
    beta = 2 * Npi - N00 - N11
    N01 = (alpha + 1j * beta) / 2
    N10 = (alpha - 1j * beta) / 2

    # Step 4: Build Choi matrix
    # J = sum_{k,l} |k><l| tensor N(|k><l|)
    # Block structure: J[block(k,l)] = N(|k><l|)
    # Index: J[k*d_out + i_out, l*d_out + j_out] = N(|k><l|)[i_out, j_out]
    choi = np.zeros((4, 4), dtype=complex)
    choi[0:2, 0:2] = N00
    choi[0:2, 2:4] = N01
    choi[2:4, 0:2] = N10
    choi[2:4, 2:4] = N11

    return choi


def _expectation_from_counts_2q(
    counts: Dict[str, int],
) -> Tuple[float, float, float]:
    """Extract <Z0>, <Z1>, and <Z0 Z1> from 2-qubit measurement counts.

    Qiskit little-endian convention: rightmost bit is qubit 0.

    Parameters
    ----------
    counts : dict
        Measurement outcome counts with 2-bit keys ("00", "01", "10", "11").

    Returns
    -------
    tuple of (exp_q0, exp_q1, exp_correlator)
        Expectation values for Z on qubit 0, Z on qubit 1, and Z⊗Z.
    """
    total = sum(counts.values())
    exp_q0 = 0.0
    exp_q1 = 0.0
    exp_corr = 0.0
    for bitstring, count in counts.items():
        # Qiskit little-endian: rightmost bit is qubit 0
        b0 = int(bitstring[-1])
        b1 = int(bitstring[-2])
        s0 = 1 - 2 * b0   # 0 -> +1, 1 -> -1
        s1 = 1 - 2 * b1
        exp_q0 += s0 * count
        exp_q1 += s1 * count
        exp_corr += s0 * s1 * count
    return exp_q0 / total, exp_q1 / total, exp_corr / total


def counts_to_choi_2q(
    counts_list: List[Dict[str, int]],
    labels: List[Tuple[str, str, str, str]],
) -> NDArray:
    """Reconstruct the 16x16 Choi matrix for a 2-qubit channel.

    Uses 144 measurement results (16 input states x 9 measurement bases)
    to reconstruct the full Choi matrix via linear inversion, analogous
    to ``counts_to_choi_1q`` but for the 2-qubit case.

    The 16 input states are tensor products of
    {|0>, |1>, |+>, |+i>} on each qubit.
    The 9 measurement bases are products of {X, Y, Z} on each qubit.

    Parameters
    ----------
    counts_list : list of dict
        Measurement counts for each of the 144 tomography circuits.
    labels : list of (inp0, inp1, basis0, basis1) tuples
        Corresponding labels for each circuit.

    Returns
    -------
    NDArray
        16x16 Choi matrix.
    """
    # ---------------------------------------------------------------
    # Step 1: Build expectation value table from counts
    # ---------------------------------------------------------------
    # exp_table[(inp0, inp1, basis0, basis1)] -> (exp_q0, exp_q1, exp_corr)
    exp_table: Dict[Tuple[str, str, str, str], Tuple[float, float, float]] = {}
    for counts, (inp0, inp1, b0, b1) in zip(counts_list, labels):
        exp_table[(inp0, inp1, b0, b1)] = _expectation_from_counts_2q(counts)

    # Pauli matrices and labels for reconstruction
    pauli_mats = {"I": _I, "X": _X, "Y": _Y, "Z": _Z}

    # ---------------------------------------------------------------
    # Step 2: Reconstruct 4x4 output density matrix for each input
    # ---------------------------------------------------------------
    # rho_out = (I⊗I + sum_{P!=II} <P> P) / 4
    # From 9 measurement bases (b0, b1) in {X,Y,Z}^2 we extract:
    #   - 9 correlators <b0⊗b1>
    #   - 3 marginals <b0⊗I> (averaged over b1)
    #   - 3 marginals <I⊗b1> (averaged over b0)
    output_states_2q: Dict[Tuple[str, str], NDArray] = {}

    for inp0 in INPUT_LABELS_1Q:
        for inp1 in INPUT_LABELS_1Q:
            # Collect marginals and correlators from the 9 basis measurements.
            # Marginals are averaged over the 3 redundant measurements.
            marginal_q0: Dict[str, float] = {}  # basis -> <basis ⊗ I>
            marginal_q1: Dict[str, float] = {}  # basis -> <I ⊗ basis>
            correlator: Dict[Tuple[str, str], float] = {}

            for b0 in MEAS_LABELS:
                q0_vals = []
                for b1 in MEAS_LABELS:
                    eq0, eq1, ecorr = exp_table[(inp0, inp1, b0, b1)]
                    correlator[(b0, b1)] = ecorr
                    q0_vals.append(eq0)
                # Average marginal for q0 over the 3 b1 measurements
                marginal_q0[b0] = float(np.mean(q0_vals))

            for b1 in MEAS_LABELS:
                q1_vals = []
                for b0 in MEAS_LABELS:
                    _, eq1, _ = exp_table[(inp0, inp1, b0, b1)]
                    q1_vals.append(eq1)
                # Average marginal for q1 over the 3 b0 measurements
                marginal_q1[b1] = float(np.mean(q1_vals))

            # Build 4x4 output density matrix
            rho_out = np.eye(4, dtype=complex) / 4.0

            # Single-qubit terms: <b⊗I> and <I⊗b>
            for b_label in MEAS_LABELS:
                # <b ⊗ I> contribution
                rho_out += marginal_q0[b_label] * np.kron(
                    pauli_mats[b_label], _I
                ) / 4.0
                # <I ⊗ b> contribution
                rho_out += marginal_q1[b_label] * np.kron(
                    _I, pauli_mats[b_label]
                ) / 4.0

            # Two-qubit correlator terms: <b0 ⊗ b1>
            for b0 in MEAS_LABELS:
                for b1 in MEAS_LABELS:
                    rho_out += correlator[(b0, b1)] * np.kron(
                        pauli_mats[b0], pauli_mats[b1]
                    ) / 4.0

            output_states_2q[(inp0, inp1)] = rho_out

    # ---------------------------------------------------------------
    # Step 3: Extract N(|i0 i1><k0 k1|) by linearity
    # ---------------------------------------------------------------
    # Strategy: factored two-step extraction.
    # Step 3a: For each fixed m1 (1Q input label for qubit 1),
    #          use 1Q linearity on qubit 0 to extract
    #          N(|i0><k0| ⊗ rho_{m1}).
    # Step 3b: For each fixed (i0, k0), use 1Q linearity on qubit 1
    #          to extract N(|i0><k0| ⊗ |i1><k1|) = N(|i0 i1><k0 k1|).
    #
    # 1Q linearity (same alpha/beta trick as counts_to_choi_1q):
    #   N(|0><0|) = output for "0"
    #   N(|1><1|) = output for "1"
    #   alpha = 2*N("+") - N("0") - N("1")
    #   beta  = 2*N("+i") - N("0") - N("1")
    #   N(|0><1|) = (alpha + i*beta) / 2
    #   N(|1><0|) = (alpha - i*beta) / 2

    def _extract_1q_linearity(
        outputs: Dict[str, NDArray],
    ) -> Dict[Tuple[int, int], NDArray]:
        """Extract N(|i><k|) from outputs keyed by 1Q input labels."""
        N00 = outputs["0"]
        N11 = outputs["1"]
        Npp = outputs["+"]
        Npi = outputs["+i"]
        alpha = 2 * Npp - N00 - N11
        beta = 2 * Npi - N00 - N11
        N01 = (alpha + 1j * beta) / 2
        N10 = (alpha - 1j * beta) / 2
        return {(0, 0): N00, (0, 1): N01, (1, 0): N10, (1, 1): N11}

    # Step 3a: Fix m1, extract qubit-0 basis elements
    # intermediate[(i0, k0, m1)] = N(|i0><k0| ⊗ rho_{m1})
    intermediate: Dict[Tuple[int, int, str], NDArray] = {}
    for m1 in INPUT_LABELS_1Q:
        q0_outputs = {m0: output_states_2q[(m0, m1)] for m0 in INPUT_LABELS_1Q}
        extracted = _extract_1q_linearity(q0_outputs)
        for (i0, k0), mat in extracted.items():
            intermediate[(i0, k0, m1)] = mat

    # Step 3b: Fix (i0, k0), extract qubit-1 basis elements
    # N_basis[(i0, k0, i1, k1)] = N(|i0 i1><k0 k1|)  [4x4 matrix]
    N_basis: Dict[Tuple[int, int, int, int], NDArray] = {}
    for i0 in range(2):
        for k0 in range(2):
            q1_outputs = {
                m1: intermediate[(i0, k0, m1)] for m1 in INPUT_LABELS_1Q
            }
            extracted = _extract_1q_linearity(q1_outputs)
            for (i1, k1), mat in extracted.items():
                N_basis[(i0, k0, i1, k1)] = mat

    # ---------------------------------------------------------------
    # Step 4: Build 16x16 Choi matrix
    # ---------------------------------------------------------------
    # J = sum_{i_in, j_in} |i_in><j_in| ⊗ N(|i_in><j_in|)
    # where i_in indexes 2-qubit computational basis: i_in = i0*2 + i1
    # J[(i_in*d_out + i_out), (j_in*d_out + j_out)]
    #   = N(|i_in><j_in|)[i_out, j_out]
    d_out_2q = 4
    choi = np.zeros((16, 16), dtype=complex)

    for i0 in range(2):
        for i1 in range(2):
            i_in = i0 * 2 + i1
            for k0 in range(2):
                for k1 in range(2):
                    j_in = k0 * 2 + k1
                    block = N_basis[(i0, k0, i1, k1)]
                    for i_out in range(d_out_2q):
                        for j_out in range(d_out_2q):
                            choi[
                                i_in * d_out_2q + i_out,
                                j_in * d_out_2q + j_out,
                            ] = block[i_out, j_out]

    return choi


def _project_psd(M: NDArray, tol: float = -1e-10) -> NDArray:
    """Project a Hermitian matrix to positive semidefinite."""
    eigvals, eigvecs = np.linalg.eigh(M)
    eigvals = np.maximum(eigvals, 0.0)
    return eigvecs @ np.diag(eigvals) @ eigvecs.conj().T


def _project_tp(choi: NDArray, d_in: int, d_out: int) -> NDArray:
    """Project Choi matrix onto the TP constraint set.

    The trace-preserving condition requires Tr_output(Choi) = I_input.
    The minimum-Frobenius-norm correction distributes the error uniformly
    over diagonal elements of each output block.

    Parameters
    ----------
    choi : NDArray
        Choi matrix of shape ``(d_in*d_out, d_in*d_out)``.
    d_in : int
        Input Hilbert space dimension.
    d_out : int
        Output Hilbert space dimension.

    Returns
    -------
    NDArray
        TP-projected Choi matrix (copy).
    """
    result = choi.copy()
    for i_in in range(d_in):
        for j_in in range(d_in):
            # Current partial trace element
            current = 0.0 + 0.0j
            for k in range(d_out):
                current += result[i_in * d_out + k, j_in * d_out + k]
            # Target: I_input[i_in, j_in]
            target = 1.0 if i_in == j_in else 0.0
            correction = (target - current) / d_out
            for k in range(d_out):
                result[i_in * d_out + k, j_in * d_out + k] += correction
    return result


def _dykstra_cptp_projection(
    choi: NDArray,
    d_in: int = 2,
    max_iter: int = 100,
    tol: float = 1e-12,
) -> NDArray:
    """Project a matrix onto the CPTP set using Dykstra's alternating projection.

    Alternates between the CP (positive semidefinite) cone and the TP
    (trace-preserving) affine subspace, using Dykstra increments to ensure
    convergence to the nearest point in the intersection.

    This implements the Projected Least-Squares (PLS) method from
    Surawy-Stepney et al. (Quantum 6, 844, 2022).

    Parameters
    ----------
    choi : NDArray
        Input Choi matrix (possibly non-physical).
    d_in : int
        Input Hilbert space dimension.
    max_iter : int
        Maximum number of alternating projection iterations.
    tol : float
        Convergence tolerance on the Frobenius norm change.

    Returns
    -------
    NDArray
        CPTP-projected Choi matrix.
    """
    d_total = choi.shape[0]
    d_out = d_total // d_in

    p_cp = np.zeros_like(choi)
    p_tp = np.zeros_like(choi)
    choi_current = choi.copy()

    for _ in range(max_iter):
        # Project onto CP (PSD cone) with Dykstra increment
        y_cp = choi_current + p_cp
        choi_cp = _project_psd(y_cp)
        p_cp = y_cp - choi_cp

        # Project onto TP with Dykstra increment
        y_tp = choi_cp + p_tp
        choi_tp = _project_tp(y_tp, d_in, d_out)
        p_tp = y_tp - choi_tp

        # Check convergence
        if np.linalg.norm(choi_tp - choi_current) < tol:
            break
        choi_current = choi_tp

    return choi_current


def choi_to_kraus(choi: NDArray, d_in: int = 2, tol: float = 1e-10) -> KrausList:
    """Convert Choi matrix to Kraus operators via eigendecomposition.

    Uses Dykstra alternating projection (PLS method) to enforce both
    complete positivity and trace preservation before extracting Kraus
    operators.

    Uses the convention J = sum_{k,l} |k><l| tensor N(|k><l|)
    so that Tr_output[J] = I_input (TP condition) and Tr[J] = d_in.

    Parameters
    ----------
    choi : NDArray
        Choi matrix of shape ``(d_in*d_out, d_in*d_out)``.
    d_in : int
        Input dimension.
    tol : float
        Eigenvalue threshold.

    Returns
    -------
    KrausList
        List of Kraus operators satisfying sum K_i^dag K_i = I.
    """
    d_total = choi.shape[0]
    d_out = d_total // d_in

    # CPTP projection via Dykstra alternating projection (PLS method)
    choi_cptp = _dykstra_cptp_projection(choi, d_in)

    eigvals, eigvecs = np.linalg.eigh(choi_cptp)

    kraus_ops = []
    for i in range(len(eigvals)):
        if eigvals[i] > tol:
            vec = eigvecs[:, i] * np.sqrt(eigvals[i])
            # Reshape: J[k*d_out + i_out, l*d_out + j_out] = N(|k><l|)[i_out,j_out]
            # Eigenvector v[k*d_out + i_out] -> K[i_out, k] (output x input)
            K = vec.reshape(d_in, d_out).T  # equivalent to reshape(d_out, d_in, order='F')
            kraus_ops.append(K)

    if len(kraus_ops) == 0:
        kraus_ops = [np.eye(d_in, dtype=complex) * 1e-10]

    return kraus_ops


# ---------------------------------------------------------------------------
# End-to-end pipeline
# ---------------------------------------------------------------------------

def run_process_tomography_1q(
    gate_circuit_fn: Callable,
    backend,
    shots: int = 4096,
    verbose: bool = True,
    timeout: int = 600,
) -> KrausList:
    """Run full 1-qubit process tomography on a backend.

    Parameters
    ----------
    gate_circuit_fn : callable
        Function that appends gate(s) to a QuantumCircuit.
        Signature: gate_circuit_fn(qc, qubit_index)
    backend
        A Qiskit-compatible backend (QI, IBM, simulator, etc.)
    shots : int
        Shots per circuit.
    verbose : bool
        Print progress.

    Returns
    -------
    KrausList
        Extracted Kraus operators.
    """
    if verbose:
        print("=== 1-Qubit Process Tomography ===")
        print(f"Building 12 tomography circuits...")

    circuits, labels = build_tomography_circuits_1q(gate_circuit_fn, shots)

    if verbose:
        print(f"Submitting {len(circuits)} circuits to {backend}...")

    # Submit all circuits
    try:
        from qiskit.compiler import transpile
        transpiled = transpile(circuits, backend, optimization_level=0)
        if verbose:
            for i, (orig, trans) in enumerate(zip(circuits, transpiled)):
                orig_ops = sum(v for k, v in orig.count_ops().items()
                               if k not in ('barrier', 'measure'))
                trans_ops = sum(v for k, v in trans.count_ops().items()
                                if k not in ('barrier', 'measure'))
                if trans_ops < orig_ops:
                    print(f"  WARNING: Circuit {i}: transpiler reduced "
                          f"{orig_ops} -> {trans_ops} ops")
        job = backend.run(transpiled, shots=shots)
        job.wait_for_final_state(timeout=timeout)
        result = job.result()
    except Exception as e:
        raise RuntimeError(f"Backend execution failed: {e}")

    if verbose:
        print("Collecting measurement results...")

    counts_list = [result.get_counts(i) for i in range(len(circuits))]

    if verbose:
        print("Reconstructing Choi matrix...")

    choi = counts_to_choi_1q(counts_list, labels)

    if verbose:
        eigvals = np.linalg.eigvalsh(choi)
        print(f"Choi eigenvalues: {np.sort(eigvals)[::-1]}")
        print(f"Choi trace: {np.trace(choi).real:.4f} (should be {2})")

    kraus_ops = choi_to_kraus(choi, d_in=2)

    if verbose:
        print(f"Extracted {len(kraus_ops)} Kraus operators")
        # Verify CPTP
        total = sum(K.conj().T @ K for K in kraus_ops)
        cptp_err = np.linalg.norm(total - np.eye(2))
        print(f"CPTP error: {cptp_err:.2e}")

    return kraus_ops


# ---------------------------------------------------------------------------
# Simulate process tomography (for testing without hardware)
# ---------------------------------------------------------------------------

def simulate_process_tomography_1q(
    kraus_ops: KrausList,
    shots: int = 100000,
    verbose: bool = True,
) -> KrausList:
    """Simulate process tomography using known Kraus operators.

    Useful for testing the reconstruction pipeline without hardware.

    Parameters
    ----------
    kraus_ops : KrausList
        The "true" Kraus operators to simulate.
    shots : int
        Simulated shots per circuit.
    verbose : bool
        Print progress.

    Returns
    -------
    KrausList
        Reconstructed Kraus operators (should match input up to gauge).
    """
    if verbose:
        print("=== Simulated Process Tomography ===")

    from .petz import apply_channel

    counts_list = []
    labels = []

    for inp in INPUT_LABELS_1Q:
        rho_in = INPUT_STATES_1Q[INPUT_LABELS_1Q.index(inp)]

        # Apply the channel
        rho_out = apply_channel(rho_in, kraus_ops)

        for basis in MEAS_LABELS:
            # Rotate to measurement basis
            if basis == "Z":
                rho_meas = rho_out
            elif basis == "X":
                H = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
                rho_meas = H @ rho_out @ H.conj().T
            elif basis == "Y":
                # Sdg then H
                Sdg = np.array([[1, 0], [0, -1j]], dtype=complex)
                H = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
                U = H @ Sdg
                rho_meas = U @ rho_out @ U.conj().T

            # Sample
            p0 = max(0, min(1, rho_meas[0, 0].real))
            n0 = np.random.binomial(shots, p0)
            n1 = shots - n0
            counts = {"0": int(n0), "1": int(n1)}

            counts_list.append(counts)
            labels.append((inp, basis))

    choi = counts_to_choi_1q(counts_list, labels)
    kraus_reconstructed = choi_to_kraus(choi, d_in=2)

    if verbose:
        print(f"Reconstructed {len(kraus_reconstructed)} Kraus operators")
        total = sum(K.conj().T @ K for K in kraus_reconstructed)
        cptp_err = np.linalg.norm(total - np.eye(2))
        print(f"CPTP error: {cptp_err:.2e}")

    return kraus_reconstructed
