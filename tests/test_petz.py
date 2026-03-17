"""Tests for tau_chrono.petz -- Petz recovery map and quantum info quantities."""

import numpy as np
import pytest

from tau_chrono.petz import (
    apply_channel,
    adjoint_channel,
    petz_recovery_map,
    apply_petz_recovery,
    fidelity,
    tau_parameter,
    relative_entropy,
    delta_D,
)
from tau_chrono.channels import (
    amplitude_damping,
    depolarizing,
    dephasing,
    unitary_channel,
    I2,
    X_GATE,
    Z_GATE,
)
from .conftest import random_density_matrix


# ---------------------------------------------------------------------------
# apply_channel
# ---------------------------------------------------------------------------


class TestApplyChannel:
    """Tests for Kraus-representation channel application."""

    def test_identity_channel_preserves_state(self, rho_plus, identity_channel):
        """N(rho) = rho when N = identity."""
        result = apply_channel(rho_plus, identity_channel)
        np.testing.assert_allclose(result, rho_plus, atol=1e-12)

    def test_preserves_trace(self, rho_plus):
        """Channel application preserves trace for arbitrary CPTP channels."""
        for kraus in [
            depolarizing(0.3),
            amplitude_damping(0.5),
            dephasing(0.2),
        ]:
            result = apply_channel(rho_plus, kraus)
            assert np.trace(result).real == pytest.approx(1.0, abs=1e-10)

    def test_preserves_trace_random_states(self):
        """Trace preservation for random density matrices."""
        rng = np.random.default_rng(42)
        for _ in range(5):
            rho = random_density_matrix(2, rng)
            kraus = depolarizing(rng.uniform(0.01, 0.99))
            result = apply_channel(rho, kraus)
            assert np.trace(result).real == pytest.approx(1.0, abs=1e-10)

    def test_output_is_hermitian(self, rho_plus):
        """Output is Hermitian."""
        result = apply_channel(rho_plus, depolarizing(0.3))
        np.testing.assert_allclose(result, result.conj().T, atol=1e-12)

    def test_output_is_psd(self, rho_plus):
        """Output is positive semidefinite."""
        result = apply_channel(rho_plus, amplitude_damping(0.5))
        eigvals = np.linalg.eigvalsh(result)
        assert np.all(eigvals >= -1e-12)

    def test_unitary_channel_rotates(self, rho_zero):
        """X gate flips |0> to |1>."""
        kraus = unitary_channel(X_GATE)
        result = apply_channel(rho_zero, kraus)
        expected = np.array([[0, 0], [0, 1]], dtype=complex)
        np.testing.assert_allclose(result, expected, atol=1e-12)


# ---------------------------------------------------------------------------
# adjoint_channel
# ---------------------------------------------------------------------------


class TestAdjointChannel:
    """Tests for the Heisenberg-picture adjoint."""

    def test_adjoint_identity(self, identity_channel):
        """Adjoint of identity is identity."""
        X = np.array([[1, 2], [3, 4]], dtype=complex)
        result = adjoint_channel(X, identity_channel)
        np.testing.assert_allclose(result, X, atol=1e-12)

    def test_adjoint_preserves_identity_operator(self):
        """N^dag(I) = I for any CPTP channel (trace preservation)."""
        I = np.eye(2, dtype=complex)
        for kraus in [depolarizing(0.3), amplitude_damping(0.5), dephasing(0.2)]:
            result = adjoint_channel(I, kraus)
            np.testing.assert_allclose(result, I, atol=1e-10)


# ---------------------------------------------------------------------------
# fidelity
# ---------------------------------------------------------------------------


class TestFidelity:
    """Tests for Uhlmann fidelity F(rho, sigma)."""

    def test_fidelity_same_state_is_one(self, rho_plus, rho_zero, sigma_mixed):
        """F(rho, rho) = 1."""
        assert fidelity(rho_plus, rho_plus) == pytest.approx(1.0, abs=1e-10)
        assert fidelity(rho_zero, rho_zero) == pytest.approx(1.0, abs=1e-10)
        assert fidelity(sigma_mixed, sigma_mixed) == pytest.approx(1.0, abs=1e-10)

    def test_fidelity_symmetric(self, rho_plus, sigma_mixed):
        """F(rho, sigma) == F(sigma, rho)."""
        f1 = fidelity(rho_plus, sigma_mixed)
        f2 = fidelity(sigma_mixed, rho_plus)
        assert f1 == pytest.approx(f2, abs=1e-10)

    def test_fidelity_symmetric_random(self):
        """Symmetry for random density matrices."""
        rng = np.random.default_rng(99)
        for _ in range(10):
            rho = random_density_matrix(2, rng)
            sigma = random_density_matrix(2, rng)
            f1 = fidelity(rho, sigma)
            f2 = fidelity(sigma, rho)
            assert f1 == pytest.approx(f2, abs=1e-10)

    def test_fidelity_in_range(self, rho_plus, sigma_mixed):
        """0 <= F <= 1."""
        f = fidelity(rho_plus, sigma_mixed)
        assert 0.0 <= f <= 1.0 + 1e-10

    def test_fidelity_orthogonal_states_is_zero(self, rho_zero, rho_one):
        """F(|0>, |1>) = 0."""
        assert fidelity(rho_zero, rho_one) == pytest.approx(0.0, abs=1e-10)

    def test_fidelity_pure_and_mixed(self, rho_zero, sigma_maximally_mixed):
        """F(|0><0|, I/2) = 1/2."""
        f = fidelity(rho_zero, sigma_maximally_mixed)
        assert f == pytest.approx(0.5, abs=1e-10)

    def test_fidelity_random_in_range(self):
        """Fidelity is in [0, 1] for random states."""
        rng = np.random.default_rng(100)
        for _ in range(20):
            rho = random_density_matrix(2, rng)
            sigma = random_density_matrix(2, rng)
            f = fidelity(rho, sigma)
            assert 0.0 - 1e-10 <= f <= 1.0 + 1e-10


# ---------------------------------------------------------------------------
# tau_parameter
# ---------------------------------------------------------------------------


class TestTauParameter:
    """Tests for tau = 1 - F(rho, R(N(rho)))."""

    def test_tau_zero_for_identity(self, rho_plus, sigma_mixed):
        """Identity channel: perfect recovery, tau = 0."""
        kraus = unitary_channel(I2)
        tau = tau_parameter(rho_plus, kraus, sigma_mixed)
        assert tau == pytest.approx(0.0, abs=1e-8)

    def test_tau_zero_for_unitary(self, rho_plus, sigma_maximally_mixed):
        """Unitary channels are perfectly recoverable, tau = 0."""
        for U in [I2, X_GATE, Z_GATE]:
            tau = tau_parameter(rho_plus, unitary_channel(U), sigma_maximally_mixed)
            assert tau == pytest.approx(0.0, abs=1e-6)

    def test_tau_in_range(self, rho_plus, sigma_mixed):
        """tau is in [0, 1]."""
        for kraus in [depolarizing(0.3), amplitude_damping(0.5)]:
            tau = tau_parameter(rho_plus, kraus, sigma_mixed)
            assert 0.0 - 1e-10 <= tau <= 1.0 + 1e-10

    def test_tau_monotonicity_depolarizing(self, rho_plus, sigma_mixed):
        """More noise -> higher tau (monotone in p for depolarizing)."""
        taus = [
            tau_parameter(rho_plus, depolarizing(p), sigma_mixed)
            for p in [0.01, 0.1, 0.3, 0.5, 0.8]
        ]
        for i in range(len(taus) - 1):
            assert taus[i] <= taus[i + 1] + 1e-8

    def test_tau_monotonicity_amplitude_damping(self, rho_plus, sigma_mixed):
        """More noise -> higher tau (monotone in gamma for amp damping)."""
        taus = [
            tau_parameter(rho_plus, amplitude_damping(g), sigma_mixed)
            for g in [0.01, 0.1, 0.3, 0.5, 0.8]
        ]
        for i in range(len(taus) - 1):
            assert taus[i] <= taus[i + 1] + 1e-8

    def test_tau_with_maximally_mixed_sigma(self, rho_plus, sigma_maximally_mixed):
        """tau computable with sigma = I/2."""
        tau = tau_parameter(rho_plus, depolarizing(0.3), sigma_maximally_mixed)
        assert 0.0 <= tau <= 1.0


# ---------------------------------------------------------------------------
# relative_entropy
# ---------------------------------------------------------------------------


class TestRelativeEntropy:
    """Tests for D(rho || sigma)."""

    def test_relative_entropy_same_state_zero(self, rho_plus, sigma_mixed):
        """D(rho, rho) = 0."""
        assert relative_entropy(rho_plus, rho_plus) == pytest.approx(0.0, abs=1e-10)
        assert relative_entropy(sigma_mixed, sigma_mixed) == pytest.approx(
            0.0, abs=1e-10
        )

    def test_relative_entropy_nonneg(self):
        """D(rho || sigma) >= 0 (Gibbs inequality)."""
        rng = np.random.default_rng(200)
        for _ in range(20):
            rho = random_density_matrix(2, rng)
            sigma = random_density_matrix(2, rng)
            D = relative_entropy(rho, sigma)
            if np.isfinite(D):
                assert D >= -1e-10, f"D(rho||sigma) = {D} < 0, violates Gibbs"

    def test_relative_entropy_pure_vs_mixed(self, rho_zero, sigma_maximally_mixed):
        """D(|0><0| || I/2) = log(2)."""
        D = relative_entropy(rho_zero, sigma_maximally_mixed)
        assert D == pytest.approx(np.log(2), abs=1e-10)

    def test_relative_entropy_support_condition(self, rho_zero, rho_one):
        """D = +inf when supp(rho) not in supp(sigma)."""
        D = relative_entropy(rho_zero, rho_one)
        assert D == float("inf")

    def test_relative_entropy_diagonal_states(self):
        """Manual computation for diagonal states."""
        # D(diag(p) || diag(q)) = sum_i p_i log(p_i/q_i)
        p = np.array([0.7, 0.3])
        q = np.array([0.4, 0.6])
        rho = np.diag(p.astype(complex))
        sigma = np.diag(q.astype(complex))
        expected = np.sum(p * np.log(p / q))
        D = relative_entropy(rho, sigma)
        assert D == pytest.approx(expected, abs=1e-10)


# ---------------------------------------------------------------------------
# delta_D (data processing inequality)
# ---------------------------------------------------------------------------


class TestDeltaD:
    """Tests for DeltaD = D(rho||sigma) - D(N(rho)||N(sigma)) >= 0."""

    def test_delta_D_nonneg(self, rho_plus, sigma_mixed):
        """Data processing inequality: DeltaD >= 0."""
        for kraus in [depolarizing(0.3), amplitude_damping(0.5), dephasing(0.2)]:
            dd = delta_D(rho_plus, sigma_mixed, kraus)
            if np.isfinite(dd):
                assert dd >= -1e-10

    def test_delta_D_nonneg_random(self):
        """DPI for random states and channels."""
        rng = np.random.default_rng(300)
        for _ in range(10):
            rho = random_density_matrix(2, rng)
            sigma = random_density_matrix(2, rng)
            kraus = depolarizing(rng.uniform(0.01, 0.99))
            dd = delta_D(rho, sigma, kraus)
            if np.isfinite(dd):
                assert dd >= -1e-10

    def test_delta_D_identity_channel_zero(self, rho_plus, sigma_mixed):
        """Identity channel: DeltaD = 0."""
        kraus = unitary_channel(I2)
        dd = delta_D(rho_plus, sigma_mixed, kraus)
        assert dd == pytest.approx(0.0, abs=1e-8)

    def test_delta_D_same_state_zero(self, rho_plus):
        """DeltaD = 0 when rho == sigma (both D before/after are zero)."""
        dd = delta_D(rho_plus, rho_plus, depolarizing(0.3))
        assert dd == pytest.approx(0.0, abs=1e-10)


# ---------------------------------------------------------------------------
# Petz recovery map
# ---------------------------------------------------------------------------


class TestPetzRecoveryMap:
    """Tests for the Petz recovery superoperator."""

    def test_petz_identity_channel(self, rho_plus, sigma_mixed):
        """Identity channel: R(N(rho)) = rho (perfect recovery)."""
        kraus = unitary_channel(I2)
        omega = apply_channel(rho_plus, kraus)
        recovered = apply_petz_recovery(omega, kraus, sigma_mixed)
        np.testing.assert_allclose(recovered, rho_plus, atol=1e-8)

    def test_petz_output_is_hermitian(self, rho_plus, sigma_mixed):
        """Petz recovery output is Hermitian."""
        kraus = depolarizing(0.3)
        omega = apply_channel(rho_plus, kraus)
        recovered = apply_petz_recovery(omega, kraus, sigma_mixed)
        np.testing.assert_allclose(recovered, recovered.conj().T, atol=1e-10)

    def test_petz_recovery_with_precomputed_superop(self, rho_plus, sigma_mixed):
        """Precomputed superoperator gives same result."""
        kraus = depolarizing(0.3)
        omega = apply_channel(rho_plus, kraus)
        S_R = petz_recovery_map(kraus, sigma_mixed)
        r1 = apply_petz_recovery(omega, kraus, sigma_mixed, S_R=S_R)
        r2 = apply_petz_recovery(omega, kraus, sigma_mixed, S_R=None)
        np.testing.assert_allclose(r1, r2, atol=1e-10)
