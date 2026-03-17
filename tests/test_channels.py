"""Tests for tau_chrono.channels -- quantum noise channels."""

import numpy as np
import pytest

from tau_chrono.channels import (
    amplitude_damping,
    depolarizing,
    dephasing,
    two_qubit_depolarizing,
    correlated_dephasing,
    local_noise,
    cnot_error,
    crosstalk_dephasing,
    amplitude_damping_2q,
    unitary_channel,
    verify_cptp,
    entanglement_fidelity,
    diamond_norm_approx,
    coherent_incoherent_decomposition,
    I2,
    X_GATE,
    Y_GATE,
    Z_GATE,
    I4,
    CNOT,
    _choi_from_kraus,
)
from tau_chrono.petz import apply_channel


# ---------------------------------------------------------------------------
# verify_cptp for ALL channel constructors
# ---------------------------------------------------------------------------


class TestVerifyCPTP:
    """All channel constructors must produce valid CPTP maps."""

    @pytest.mark.parametrize("gamma", [0.0, 0.01, 0.1, 0.5, 0.99, 1.0])
    def test_amplitude_damping_cptp(self, gamma):
        assert verify_cptp(amplitude_damping(gamma))

    @pytest.mark.parametrize("p", [0.0, 0.01, 0.1, 0.5, 0.99, 1.0])
    def test_depolarizing_cptp(self, p):
        assert verify_cptp(depolarizing(p))

    @pytest.mark.parametrize("p", [0.0, 0.01, 0.1, 0.5, 0.99, 1.0])
    def test_dephasing_cptp(self, p):
        assert verify_cptp(dephasing(p))

    @pytest.mark.parametrize("p", [0.0, 0.1, 0.5, 1.0])
    def test_two_qubit_depolarizing_cptp(self, p):
        assert verify_cptp(two_qubit_depolarizing(p))

    @pytest.mark.parametrize("p", [0.0, 0.1, 0.5, 1.0])
    def test_correlated_dephasing_cptp(self, p):
        assert verify_cptp(correlated_dephasing(p))

    @pytest.mark.parametrize("p", [0.0, 0.1, 0.5])
    def test_cnot_error_cptp(self, p):
        assert verify_cptp(cnot_error(p))

    @pytest.mark.parametrize("p,eps", [(0.1, 0.0), (0.3, 0.05), (0.5, 0.1)])
    def test_crosstalk_dephasing_cptp(self, p, eps):
        assert verify_cptp(crosstalk_dephasing(p, eps))

    @pytest.mark.parametrize("gamma", [0.0, 0.1, 0.5, 1.0])
    def test_amplitude_damping_2q_cptp(self, gamma):
        assert verify_cptp(amplitude_damping_2q(gamma))

    def test_unitary_channel_cptp(self):
        """Unitary channels are always CPTP."""
        for U in [I2, X_GATE, Y_GATE, Z_GATE]:
            assert verify_cptp(unitary_channel(U))

    def test_local_noise_cptp(self):
        """Tensor product of CPTP channels is CPTP."""
        kraus = local_noise(depolarizing(0.1), dephasing(0.2))
        assert verify_cptp(kraus)


# ---------------------------------------------------------------------------
# Edge cases: parameter boundaries
# ---------------------------------------------------------------------------


class TestEdgeCases:
    """Boundary behaviour of channel parameters."""

    def test_amplitude_damping_gamma0_is_identity(self, rho_plus, rho_zero):
        """gamma=0 should act as identity."""
        kraus = amplitude_damping(0.0)
        np.testing.assert_allclose(apply_channel(rho_plus, kraus), rho_plus, atol=1e-12)
        np.testing.assert_allclose(apply_channel(rho_zero, kraus), rho_zero, atol=1e-12)

    def test_amplitude_damping_gamma1_full_damping(self):
        """gamma=1: any state maps to |0><0|."""
        kraus = amplitude_damping(1.0)
        zero = np.array([[1, 0], [0, 0]], dtype=complex)
        one = np.array([[0, 0], [0, 1]], dtype=complex)
        np.testing.assert_allclose(apply_channel(one, kraus), zero, atol=1e-12)
        np.testing.assert_allclose(apply_channel(zero, kraus), zero, atol=1e-12)

    def test_depolarizing_p0_is_identity(self, rho_plus):
        """p=0 depolarizing acts as identity."""
        kraus = depolarizing(0.0)
        np.testing.assert_allclose(apply_channel(rho_plus, kraus), rho_plus, atol=1e-12)

    def test_depolarizing_p1_fully_mixed(self, rho_zero):
        """p=1: full depolarizing maps to I/2."""
        kraus = depolarizing(1.0)
        result = apply_channel(rho_zero, kraus)
        expected = np.eye(2, dtype=complex) / 2.0
        np.testing.assert_allclose(result, expected, atol=1e-12)

    def test_dephasing_p0_is_identity(self, rho_plus):
        """p=0 dephasing acts as identity."""
        kraus = dephasing(0.0)
        np.testing.assert_allclose(apply_channel(rho_plus, kraus), rho_plus, atol=1e-12)

    def test_dephasing_p1_removes_coherence(self, rho_plus):
        """p=1: full dephasing kills off-diagonal elements."""
        # N(rho) = (1-p)*rho + p*Z*rho*Z.  At p=1: Z*rho*Z
        # For |+><+|: Z|+><+|Z = |-><-|
        # But actually, p=1 means K0 = 0*I, K1 = 1*Z; so N(rho) = Z*rho*Z
        kraus = dephasing(1.0)
        result = apply_channel(rho_plus, kraus)
        # Z * |+><+| * Z = |-><-| = [[0.5, -0.5], [-0.5, 0.5]]
        expected = np.array([[0.5, -0.5], [-0.5, 0.5]], dtype=complex)
        np.testing.assert_allclose(result, expected, atol=1e-12)

    def test_dephasing_half_kills_offdiagonal(self, rho_plus):
        """p=0.5: (1-p)rho + pZrhoZ averages to zero off-diagonal."""
        kraus = dephasing(0.5)
        result = apply_channel(rho_plus, kraus)
        # Off-diagonal should be zero for |+><+|:
        # (0.5)*0.5 + (0.5)*(-0.5) = 0
        np.testing.assert_allclose(result[0, 1], 0.0, atol=1e-12)
        np.testing.assert_allclose(result[1, 0], 0.0, atol=1e-12)

    def test_two_qubit_depolarizing_p0_identity(self, bell_state):
        """p=0: two-qubit depolarizing acts as identity."""
        kraus = two_qubit_depolarizing(0.0)
        np.testing.assert_allclose(
            apply_channel(bell_state, kraus), bell_state, atol=1e-12
        )

    def test_two_qubit_depolarizing_p1_fully_mixed(self, bell_state):
        """p=1: maps to I/4."""
        kraus = two_qubit_depolarizing(1.0)
        result = apply_channel(bell_state, kraus)
        expected = np.eye(4, dtype=complex) / 4.0
        np.testing.assert_allclose(result, expected, atol=1e-12)


# ---------------------------------------------------------------------------
# Channel comparison metrics
# ---------------------------------------------------------------------------


class TestEntanglementFidelity:
    """Tests for entanglement_fidelity(kraus_ops)."""

    def test_identity_channel_fidelity_one(self, identity_channel):
        """Identity channel has F_e = 1."""
        assert entanglement_fidelity(identity_channel) == pytest.approx(1.0, abs=1e-12)

    def test_unitary_channel_fidelity_one(self):
        """Unitary channels preserve entanglement fidelity 1 only for I."""
        # For X gate: F_e = |Tr(X)|^2 / 4 = 0
        assert entanglement_fidelity(unitary_channel(X_GATE)) == pytest.approx(
            0.0, abs=1e-12
        )

    def test_depolarizing_fidelity_decreases(self):
        """Entanglement fidelity decreases with noise strength."""
        f0 = entanglement_fidelity(depolarizing(0.0))
        f1 = entanglement_fidelity(depolarizing(0.1))
        f2 = entanglement_fidelity(depolarizing(0.5))
        assert f0 > f1 > f2

    def test_fidelity_in_range(self):
        """Entanglement fidelity is in [0, 1]."""
        for p in [0.0, 0.1, 0.5, 1.0]:
            f = entanglement_fidelity(depolarizing(p))
            assert 0.0 <= f <= 1.0 + 1e-12


class TestDiamondNormApprox:
    """Tests for diamond_norm_approx."""

    def test_same_channel_distance_zero(self):
        """d(N, N) = 0."""
        kraus = depolarizing(0.3)
        assert diamond_norm_approx(kraus, kraus) == pytest.approx(0.0, abs=1e-10)

    def test_identity_vs_identity(self, identity_channel):
        """d(I, I) = 0."""
        assert diamond_norm_approx(identity_channel, identity_channel) == pytest.approx(
            0.0, abs=1e-12
        )

    def test_distance_nonneg(self):
        """Diamond norm distance is >= 0."""
        k1 = depolarizing(0.1)
        k2 = depolarizing(0.5)
        assert diamond_norm_approx(k1, k2) >= -1e-12

    def test_distance_increases_with_param_diff(self):
        """Greater parameter difference => larger distance."""
        k0 = depolarizing(0.0)
        k1 = depolarizing(0.1)
        k2 = depolarizing(0.5)
        d01 = diamond_norm_approx(k0, k1)
        d02 = diamond_norm_approx(k0, k2)
        assert d02 > d01


class TestCoherentIncoherentDecomposition:
    """Tests for coherent_incoherent_decomposition."""

    def test_identity_channel(self, identity_channel):
        """Identity has no error, fractions should be zero."""
        result = coherent_incoherent_decomposition(identity_channel)
        assert result["avg_gate_fidelity"] == pytest.approx(1.0, abs=1e-12)
        assert result["coherent_frac"] == pytest.approx(0.0, abs=1e-12)
        assert result["incoherent_frac"] == pytest.approx(0.0, abs=1e-12)

    def test_fractions_sum_to_one_or_zero(self):
        """Coherent + incoherent fractions sum to 1 (or both 0 for no error)."""
        for p in [0.01, 0.1, 0.5]:
            result = coherent_incoherent_decomposition(depolarizing(p))
            total = result["coherent_frac"] + result["incoherent_frac"]
            assert total == pytest.approx(1.0, abs=1e-10)

    def test_depolarizing_decomposition_valid(self):
        """Depolarizing noise decomposition returns valid fractions."""
        result = coherent_incoherent_decomposition(depolarizing(0.1))
        assert 0.0 <= result["coherent_frac"] <= 1.0
        assert 0.0 <= result["incoherent_frac"] <= 1.0

    def test_unitarity_in_range(self):
        """Unitarity is in [0, 1]."""
        for p in [0.0, 0.1, 0.5, 1.0]:
            result = coherent_incoherent_decomposition(depolarizing(p))
            assert 0.0 <= result["unitarity"] <= 1.0 + 1e-10

    def test_unitary_channel_unitarity_one(self):
        """Unitary channel has unitarity = 1."""
        result = coherent_incoherent_decomposition(unitary_channel(X_GATE))
        assert result["unitarity"] == pytest.approx(1.0, abs=1e-10)

    def test_returns_valid_dict_keys(self):
        """Output contains required keys."""
        result = coherent_incoherent_decomposition(depolarizing(0.1))
        assert set(result.keys()) == {
            "coherent_frac",
            "incoherent_frac",
            "unitarity",
            "avg_gate_fidelity",
        }


# ---------------------------------------------------------------------------
# Choi matrix construction
# ---------------------------------------------------------------------------


class TestChoiFromKraus:
    """Tests for internal _choi_from_kraus helper."""

    def test_identity_choi_is_bell(self, identity_channel):
        """Choi of identity channel is the (unnormalized) Bell state projector."""
        choi = _choi_from_kraus(identity_channel)
        # |Phi+> = (|00> + |11>)/sqrt(2) -> projector * 2 (unnormalized)
        bell = np.zeros((4, 4), dtype=complex)
        for i in range(2):
            for j in range(2):
                bell[i * 2 + i, j * 2 + j] = 1.0
        # Choi with our convention: J[i*d+i_out, j*d+j_out] = N(|i><j|)[i_out, j_out]
        # For identity: N(|i><j|) = |i><j|, so J = sum_{i,j} |i><j| x |i><j|
        np.testing.assert_allclose(choi, bell, atol=1e-12)

    def test_choi_hermitian(self):
        """Choi matrix is Hermitian for any CPTP channel."""
        for kraus in [depolarizing(0.1), dephasing(0.3), amplitude_damping(0.2)]:
            choi = _choi_from_kraus(kraus)
            np.testing.assert_allclose(choi, choi.conj().T, atol=1e-12)

    def test_choi_trace_equals_dim(self):
        """Tr(Choi) = d for a trace-preserving channel."""
        for kraus in [depolarizing(0.1), amplitude_damping(0.5)]:
            choi = _choi_from_kraus(kraus)
            d = kraus[0].shape[0]
            assert np.trace(choi).real == pytest.approx(d, abs=1e-10)

    def test_choi_psd(self):
        """Choi matrix is PSD for a CP channel."""
        for kraus in [depolarizing(0.3), amplitude_damping(0.5)]:
            choi = _choi_from_kraus(kraus)
            eigvals = np.linalg.eigvalsh(choi)
            assert np.all(eigvals >= -1e-10)


# ---------------------------------------------------------------------------
# Assertion on invalid parameters
# ---------------------------------------------------------------------------


class TestInvalidParameters:
    """Channels reject out-of-range parameters."""

    def test_amplitude_damping_negative(self):
        with pytest.raises(AssertionError):
            amplitude_damping(-0.1)

    def test_amplitude_damping_above_one(self):
        with pytest.raises(AssertionError):
            amplitude_damping(1.1)

    def test_depolarizing_negative(self):
        with pytest.raises(AssertionError):
            depolarizing(-0.1)

    def test_depolarizing_above_one(self):
        with pytest.raises(AssertionError):
            depolarizing(1.1)

    def test_dephasing_negative(self):
        with pytest.raises(AssertionError):
            dephasing(-0.1)

    def test_dephasing_above_one(self):
        with pytest.raises(AssertionError):
            dephasing(1.1)
