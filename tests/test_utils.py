"""Tests for tau_chrono.utils -- matrix functions and numerical utilities."""

import numpy as np
import pytest
from scipy.linalg import expm

from tau_chrono.utils import (
    matrix_sqrt,
    matrix_inv_sqrt,
    matrix_log,
    matrix_power_psd,
    trace_norm,
    commutator_norm,
    clip_hermitian,
    safe_eigendecomposition,
    _EPS,
)
from .conftest import random_psd_matrix, random_hermitian_matrix


# ---------------------------------------------------------------------------
# matrix_sqrt
# ---------------------------------------------------------------------------


class TestMatrixSqrt:
    """Tests for matrix_sqrt(A)."""

    def test_sqrt_squared_equals_original_random_psd(self):
        """matrix_sqrt(A)^2 == A for random PSD matrices."""
        for seed in [0, 1, 2, 42, 100]:
            rng = np.random.default_rng(seed)
            A = random_psd_matrix(4, rng)
            sqrtA = matrix_sqrt(A)
            reconstructed = sqrtA @ sqrtA
            np.testing.assert_allclose(reconstructed, A, atol=1e-10)

    def test_sqrt_identity(self):
        """sqrt(I) == I."""
        I4 = np.eye(4, dtype=complex)
        result = matrix_sqrt(I4)
        np.testing.assert_allclose(result, I4, atol=1e-12)

    def test_sqrt_diagonal(self):
        """sqrt of diagonal matrix has sqrt eigenvalues."""
        D = np.diag([4.0, 9.0, 16.0]).astype(complex)
        result = matrix_sqrt(D)
        expected = np.diag([2.0, 3.0, 4.0]).astype(complex)
        np.testing.assert_allclose(result, expected, atol=1e-12)

    def test_sqrt_rank_deficient(self):
        """sqrt handles rank-deficient (singular) PSD matrices."""
        A = np.diag([1.0, 0.0, 4.0]).astype(complex)
        sqrtA = matrix_sqrt(A)
        reconstructed = sqrtA @ sqrtA
        np.testing.assert_allclose(reconstructed, A, atol=1e-12)

    def test_sqrt_output_is_hermitian(self):
        """sqrt(A) is Hermitian when A is Hermitian PSD."""
        A = random_psd_matrix(3, np.random.default_rng(7))
        sqrtA = matrix_sqrt(A)
        np.testing.assert_allclose(sqrtA, sqrtA.conj().T, atol=1e-12)

    def test_sqrt_output_is_psd(self):
        """sqrt(A) is PSD when A is PSD."""
        A = random_psd_matrix(3, np.random.default_rng(8))
        sqrtA = matrix_sqrt(A)
        eigvals = np.linalg.eigvalsh(sqrtA)
        assert np.all(eigvals >= -1e-12)


# ---------------------------------------------------------------------------
# matrix_inv_sqrt
# ---------------------------------------------------------------------------


class TestMatrixInvSqrt:
    """Tests for matrix_inv_sqrt(A)."""

    def test_inv_sqrt_times_sqrt_is_identity_on_support(self):
        """matrix_inv_sqrt(A) @ matrix_sqrt(A) == I on the support of A."""
        for seed in [10, 20, 30]:
            rng = np.random.default_rng(seed)
            A = random_psd_matrix(3, rng)
            sqrtA = matrix_sqrt(A)
            inv_sqrtA = matrix_inv_sqrt(A)
            product = inv_sqrtA @ sqrtA

            # Project onto support
            eigvals, U = np.linalg.eigh(A)
            support_mask = eigvals > _EPS
            P_support = U[:, support_mask] @ U[:, support_mask].conj().T
            expected = P_support
            np.testing.assert_allclose(product, expected, atol=1e-10)

    def test_inv_sqrt_full_rank(self):
        """For full-rank A, inv_sqrt(A) @ sqrt(A) == I."""
        A = np.diag([1.0, 4.0, 9.0]).astype(complex)
        product = matrix_inv_sqrt(A) @ matrix_sqrt(A)
        np.testing.assert_allclose(product, np.eye(3, dtype=complex), atol=1e-12)

    def test_inv_sqrt_rank_deficient(self):
        """Pseudoinverse handles zero eigenvalues gracefully."""
        A = np.diag([4.0, 0.0, 9.0]).astype(complex)
        inv_sqrtA = matrix_inv_sqrt(A)
        expected = np.diag([0.5, 0.0, 1.0 / 3.0]).astype(complex)
        np.testing.assert_allclose(inv_sqrtA, expected, atol=1e-12)

    def test_inv_sqrt_identity(self):
        """inv_sqrt(I) == I."""
        I3 = np.eye(3, dtype=complex)
        result = matrix_inv_sqrt(I3)
        np.testing.assert_allclose(result, I3, atol=1e-12)


# ---------------------------------------------------------------------------
# matrix_log
# ---------------------------------------------------------------------------


class TestMatrixLog:
    """Tests for matrix_log(A)."""

    def test_log_of_exp_equals_original_hermitian(self):
        """matrix_log(expm(H)) == H for Hermitian H with small norm."""
        rng = np.random.default_rng(55)
        H = random_hermitian_matrix(3, rng)
        # Scale down to keep eigenvalues reasonable
        H = H / np.linalg.norm(H) * 0.5
        expH = expm(H)
        logexpH = matrix_log(expH)
        np.testing.assert_allclose(logexpH, H, atol=1e-10)

    def test_log_identity(self):
        """log(I) == 0."""
        I3 = np.eye(3, dtype=complex)
        result = matrix_log(I3)
        np.testing.assert_allclose(result, np.zeros((3, 3), dtype=complex), atol=1e-12)

    def test_log_diagonal(self):
        """log(diag(e, e^2)) == diag(1, 2)."""
        D = np.diag([np.e, np.e**2]).astype(complex)
        result = matrix_log(D)
        expected = np.diag([1.0, 2.0]).astype(complex)
        np.testing.assert_allclose(result, expected, atol=1e-12)

    def test_log_output_is_hermitian(self):
        """log(A) is Hermitian when A is Hermitian PSD."""
        A = random_psd_matrix(3, np.random.default_rng(9))
        logA = matrix_log(A)
        # Filter out -inf entries before checking Hermiticity
        finite_mask = np.isfinite(logA)
        if np.all(finite_mask):
            np.testing.assert_allclose(logA, logA.conj().T, atol=1e-12)


# ---------------------------------------------------------------------------
# matrix_power_psd
# ---------------------------------------------------------------------------


class TestMatrixPowerPsd:
    """Tests for the general matrix_power_psd."""

    def test_power_zero_is_identity_on_support(self):
        """A^0 is the projector onto the support of A."""
        A = np.diag([3.0, 0.0, 5.0]).astype(complex)
        result = matrix_power_psd(A, 0.0)
        # 0^0 = 0 in this implementation (mask > _EPS)
        # Actually for zero eigenvalues, the mask filters them out,
        # so it's the projector onto the support
        expected = np.diag([1.0, 0.0, 1.0]).astype(complex)
        np.testing.assert_allclose(result, expected, atol=1e-12)

    def test_power_one_is_identity_transform(self):
        """A^1 == A."""
        A = random_psd_matrix(3, np.random.default_rng(77))
        result = matrix_power_psd(A, 1.0)
        np.testing.assert_allclose(result, A, atol=1e-12)

    def test_negative_power_is_pseudoinverse(self):
        """A^{-1} is the pseudoinverse on the support."""
        A = np.diag([2.0, 0.0, 5.0]).astype(complex)
        result = matrix_power_psd(A, -1.0)
        expected = np.diag([0.5, 0.0, 0.2]).astype(complex)
        np.testing.assert_allclose(result, expected, atol=1e-12)


# ---------------------------------------------------------------------------
# trace_norm
# ---------------------------------------------------------------------------


class TestTraceNorm:
    """Tests for trace_norm(A)."""

    def test_trace_norm_nonneg(self):
        """Trace norm is always >= 0."""
        rng = np.random.default_rng(3)
        for _ in range(10):
            A = rng.standard_normal((3, 3)) + 1j * rng.standard_normal((3, 3))
            assert trace_norm(A) >= 0.0

    def test_trace_norm_zero_matrix(self):
        """||0|| = 0."""
        assert trace_norm(np.zeros((2, 2), dtype=complex)) == pytest.approx(0.0)

    def test_trace_norm_identity(self):
        """||I_d|| = d."""
        for d in [2, 3, 4]:
            assert trace_norm(np.eye(d, dtype=complex)) == pytest.approx(d, abs=1e-12)

    def test_trace_norm_rank_one(self):
        """Trace norm of a rank-1 matrix equals its operator norm."""
        v = np.array([1, 2, 3], dtype=complex)
        A = np.outer(v, v.conj())
        # For rank-1 PSD: trace_norm = trace = |v|^2
        expected = np.sum(np.abs(v) ** 2)
        assert trace_norm(A) == pytest.approx(expected, abs=1e-12)

    def test_trace_norm_pauli_x(self):
        """||X|| = 2 (eigenvalues +/-1)."""
        X = np.array([[0, 1], [1, 0]], dtype=complex)
        assert trace_norm(X) == pytest.approx(2.0, abs=1e-12)


# ---------------------------------------------------------------------------
# commutator_norm
# ---------------------------------------------------------------------------


class TestCommutatorNorm:
    """Tests for commutator_norm(A, B)."""

    def test_commutator_with_self_is_zero(self):
        """[A, A] = 0."""
        A = random_psd_matrix(3, np.random.default_rng(4))
        assert commutator_norm(A, A) == pytest.approx(0.0, abs=1e-12)

    def test_commutator_with_identity_is_zero(self):
        """[A, I] = 0."""
        A = random_psd_matrix(3, np.random.default_rng(5))
        I3 = np.eye(3, dtype=complex)
        assert commutator_norm(A, I3) == pytest.approx(0.0, abs=1e-12)

    def test_commutator_with_scalar_multiple_is_zero(self):
        """[A, cA] = 0."""
        A = random_psd_matrix(2, np.random.default_rng(6))
        assert commutator_norm(A, 3.0 * A) == pytest.approx(0.0, abs=1e-12)

    def test_commutator_pauli_xy(self):
        """[X, Y] = 2iZ, ||[X, Y]||_F = 2*||Z||_F = 2*sqrt(2)."""
        X = np.array([[0, 1], [1, 0]], dtype=complex)
        Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
        # [X, Y] = XY - YX = 2iZ, ||2iZ||_F = 2 * sqrt(2)
        expected = 2.0 * np.sqrt(2.0)
        assert commutator_norm(X, Y) == pytest.approx(expected, abs=1e-12)

    def test_commutator_is_nonneg(self):
        """Commutator norm is always >= 0."""
        rng = np.random.default_rng(7)
        A = rng.standard_normal((3, 3)) + 1j * rng.standard_normal((3, 3))
        B = rng.standard_normal((3, 3)) + 1j * rng.standard_normal((3, 3))
        assert commutator_norm(A, B) >= 0.0


# ---------------------------------------------------------------------------
# clip_hermitian
# ---------------------------------------------------------------------------


class TestClipHermitian:
    """Tests for clip_hermitian(A)."""

    def test_output_is_hermitian(self):
        """clip_hermitian always produces Hermitian output."""
        rng = np.random.default_rng(11)
        for _ in range(10):
            A = rng.standard_normal((4, 4)) + 1j * rng.standard_normal((4, 4))
            H = clip_hermitian(A)
            np.testing.assert_allclose(H, H.conj().T, atol=1e-15)

    def test_hermitian_input_unchanged(self):
        """Hermitian input is returned unchanged."""
        H = random_hermitian_matrix(3, np.random.default_rng(12))
        result = clip_hermitian(H)
        np.testing.assert_allclose(result, H, atol=1e-15)

    def test_clip_identity(self):
        """clip_hermitian(I) == I."""
        I3 = np.eye(3, dtype=complex)
        np.testing.assert_allclose(clip_hermitian(I3), I3, atol=1e-15)


# ---------------------------------------------------------------------------
# safe_eigendecomposition
# ---------------------------------------------------------------------------


class TestSafeEigendecomposition:
    """Tests for safe_eigendecomposition(A)."""

    def test_eigenvalues_nonneg(self):
        """All eigenvalues are clipped to >= 0."""
        # Create a matrix with small negative eigenvalues due to numerics
        A = np.array([[1.0, 0.5], [0.5, 1e-14]], dtype=complex)
        eigvals, _ = safe_eigendecomposition(A)
        assert np.all(eigvals >= 0.0)

    def test_eigenvectors_unitary(self):
        """Eigenvector matrix is unitary."""
        A = random_psd_matrix(3, np.random.default_rng(13))
        _, U = safe_eigendecomposition(A)
        np.testing.assert_allclose(U @ U.conj().T, np.eye(3, dtype=complex), atol=1e-12)

    def test_reconstruction(self):
        """U @ diag(eigvals) @ U^dag reconstructs A."""
        A = random_psd_matrix(3, np.random.default_rng(14))
        eigvals, U = safe_eigendecomposition(A)
        reconstructed = (U * eigvals[np.newaxis, :]) @ U.conj().T
        # A is Hermitianized inside safe_eigendecomposition
        A_herm = 0.5 * (A + A.conj().T)
        np.testing.assert_allclose(reconstructed, A_herm, atol=1e-10)
