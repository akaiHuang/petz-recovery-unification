"""Tests for tau_chrono.bayesian -- Bayesian composition engine."""

import numpy as np
import pytest

from tau_chrono.bayesian import (
    bayesian_compose,
    compose_kraus,
    compose_kraus_compressed,
    GateResult,
    CompositionResult,
    _classify_gate,
)
from tau_chrono.channels import (
    depolarizing,
    dephasing,
    amplitude_damping,
    verify_cptp,
    unitary_channel,
    I2,
)
from tau_chrono.petz import apply_channel
from .conftest import random_density_matrix


# ---------------------------------------------------------------------------
# compose_kraus
# ---------------------------------------------------------------------------


class TestComposeKraus:
    """Tests for Kraus channel composition."""

    def test_single_channel_returns_same(self):
        """Composing a single channel returns the same channel."""
        kraus = depolarizing(0.1)
        composed = compose_kraus([kraus])
        assert len(composed) == len(kraus)
        for K_orig, K_comp in zip(kraus, composed):
            np.testing.assert_allclose(K_comp, K_orig, atol=1e-15)

    def test_empty_list_returns_identity(self):
        """Composing no channels returns identity."""
        composed = compose_kraus([])
        assert len(composed) == 1
        np.testing.assert_allclose(composed[0], np.eye(2, dtype=complex), atol=1e-15)

    def test_composition_count(self):
        """Number of composed Kraus ops = product of per-channel counts."""
        k1 = depolarizing(0.1)  # 4 ops
        k2 = dephasing(0.2)  # 2 ops
        composed = compose_kraus([k1, k2])
        assert len(composed) == len(k1) * len(k2)

    def test_composition_is_cptp(self):
        """Composed channel is CPTP."""
        k1 = depolarizing(0.1)
        k2 = amplitude_damping(0.2)
        composed = compose_kraus([k1, k2])
        assert verify_cptp(composed, tol=1e-8)

    def test_two_identities_compose_to_identity(self):
        """I o I = I."""
        k_id = unitary_channel(I2)
        composed = compose_kraus([k_id, k_id])
        assert len(composed) == 1
        np.testing.assert_allclose(composed[0], I2, atol=1e-12)

    def test_composition_order_matters(self):
        """N2 o N1 != N1 o N2 in general."""
        k1 = amplitude_damping(0.3)
        k2 = depolarizing(0.2)
        rho = np.array([[0.5, 0.3], [0.3, 0.5]], dtype=complex)
        c12 = compose_kraus([k1, k2])
        c21 = compose_kraus([k2, k1])
        r12 = apply_channel(rho, c12)
        r21 = apply_channel(rho, c21)
        # Should differ (non-commuting channels)
        assert np.linalg.norm(r12 - r21) > 1e-6

    def test_triple_composition_cptp(self):
        """Three channels composed remain CPTP."""
        channels = [depolarizing(0.05), dephasing(0.1), amplitude_damping(0.1)]
        composed = compose_kraus(channels)
        assert verify_cptp(composed, tol=1e-8)


# ---------------------------------------------------------------------------
# compose_kraus_compressed
# ---------------------------------------------------------------------------


class TestComposeKrausCompressed:
    """Tests for SVD-compressed Kraus composition."""

    def test_compressed_preserves_cptp(self):
        """Compressed composition maintains CPTP within tolerance."""
        channels = [depolarizing(0.05)] * 5
        composed = compose_kraus_compressed(channels, max_ops=16)
        assert verify_cptp(composed, tol=1e-4)

    def test_compressed_single_channel(self):
        """Single channel returns same."""
        kraus = depolarizing(0.1)
        composed = compose_kraus_compressed([kraus])
        assert len(composed) == len(kraus)

    def test_compressed_empty_list(self):
        """Empty list returns identity."""
        composed = compose_kraus_compressed([])
        assert len(composed) == 1

    def test_compressed_limits_ops(self):
        """Compressed composition reduces operator count vs naive."""
        channels = [depolarizing(0.05)] * 10
        exact = compose_kraus(channels)
        composed = compose_kraus_compressed(
            channels, max_ops=8, compress_threshold=32
        )
        # Compressed should not exceed the exact count
        assert len(composed) <= len(exact)

    def test_compressed_vs_exact_similar_output(self):
        """Compressed composition produces similar output state to exact."""
        channels = [depolarizing(0.05), dephasing(0.05), amplitude_damping(0.05)]
        rho = np.array([[0.7, 0.2], [0.2, 0.3]], dtype=complex)
        exact = compose_kraus(channels)
        compressed = compose_kraus_compressed(channels, max_ops=32)
        r_exact = apply_channel(rho, exact)
        r_compressed = apply_channel(rho, compressed)
        np.testing.assert_allclose(r_exact, r_compressed, atol=1e-6)


# ---------------------------------------------------------------------------
# bayesian_compose
# ---------------------------------------------------------------------------


class TestBayesianCompose:
    """Tests for the main Bayesian composition engine."""

    @pytest.fixture
    def simple_circuit(self):
        """A simple 3-gate circuit."""
        channels = [depolarizing(0.05), dephasing(0.03), amplitude_damping(0.02)]
        sigma_0 = np.eye(2, dtype=complex) / 2.0
        rho = np.array([[0.5, 0.5], [0.5, 0.5]], dtype=complex)
        return channels, sigma_0, rho

    def test_composition_holds(self, simple_circuit):
        """sqrt(tau_total) <= sum(sqrt(tau_i^eff))."""
        channels, sigma_0, rho = simple_circuit
        result = bayesian_compose(channels, sigma_0, rho)
        assert result.composition_holds

    def test_composition_holds_random_circuits(self):
        """Composition inequality holds for random circuits up to depth 20."""
        rng = np.random.default_rng(42)
        for depth in [2, 5, 10, 15, 20]:
            channels = []
            for _ in range(depth):
                choice = rng.choice(3)
                p = rng.uniform(0.01, 0.15)
                if choice == 0:
                    channels.append(depolarizing(p))
                elif choice == 1:
                    channels.append(dephasing(p))
                else:
                    channels.append(amplitude_damping(p))

            rho = random_density_matrix(2, rng)
            sigma_0 = random_density_matrix(2, rng)
            result = bayesian_compose(channels, sigma_0, rho)
            assert result.composition_holds, (
                f"Composition inequality failed at depth {depth}: "
                f"LHS={result.composition_lhs:.6f} > RHS={result.composition_rhs:.6f}"
            )

    def test_tau_bayesian_le_multiplicative(self, simple_circuit):
        """tau_bayesian_total <= tau_multiplicative_total."""
        channels, sigma_0, rho = simple_circuit
        result = bayesian_compose(channels, sigma_0, rho)
        assert result.tau_bayesian_total <= result.tau_multiplicative_total + 1e-10

    def test_gate_results_count(self, simple_circuit):
        """Number of gate results matches circuit depth."""
        channels, sigma_0, rho = simple_circuit
        result = bayesian_compose(channels, sigma_0, rho)
        assert len(result.gate_results) == len(channels)

    def test_channel_names_default(self, simple_circuit):
        """Default channel names are Gate_0, Gate_1, etc."""
        channels, sigma_0, rho = simple_circuit
        result = bayesian_compose(channels, sigma_0, rho)
        for i, gr in enumerate(result.gate_results):
            assert gr.channel_name == f"Gate_{i}"

    def test_channel_names_custom(self, simple_circuit):
        """Custom channel names are propagated."""
        channels, sigma_0, rho = simple_circuit
        names = ["Depol", "Dephase", "AmpDamp"]
        result = bayesian_compose(channels, sigma_0, rho, channel_names=names)
        for gr, name in zip(result.gate_results, names):
            assert gr.channel_name == name

    def test_improvement_percent_nonneg(self):
        """Improvement percent is >= 0 when depth > 1."""
        channels = [depolarizing(0.1)] * 5
        sigma_0 = np.eye(2, dtype=complex) / 2.0
        rho = np.array([[0.5, 0.5], [0.5, 0.5]], dtype=complex)
        result = bayesian_compose(channels, sigma_0, rho)
        assert result.improvement_percent >= -1e-6

    def test_single_gate_circuit(self):
        """Single gate circuit works correctly."""
        channels = [depolarizing(0.1)]
        sigma_0 = np.eye(2, dtype=complex) / 2.0
        rho = np.array([[0.5, 0.5], [0.5, 0.5]], dtype=complex)
        result = bayesian_compose(channels, sigma_0, rho)
        assert len(result.gate_results) == 1
        assert result.composition_holds


# ---------------------------------------------------------------------------
# memory_alpha (non-Markovian memory)
# ---------------------------------------------------------------------------


class TestMemoryAlpha:
    """Tests for non-Markovian memory parameter."""

    def test_memory_alpha_1_equals_no_memory(self):
        """memory_alpha=1.0 is equivalent to no memory (fully Markovian)."""
        channels = [depolarizing(0.05)] * 3
        sigma_0 = np.eye(2, dtype=complex) / 2.0
        rho = np.array([[0.5, 0.5], [0.5, 0.5]], dtype=complex)
        r_no_mem = bayesian_compose(channels, sigma_0, rho, memory_alpha=None)
        r_alpha1 = bayesian_compose(channels, sigma_0, rho, memory_alpha=1.0)
        # The total tau should be the same (since memory_alpha=1 => no correction)
        assert r_no_mem.tau_bayesian_total == pytest.approx(
            r_alpha1.tau_bayesian_total, abs=1e-10
        )

    def test_memory_alpha_changes_tau(self):
        """Non-trivial memory_alpha changes effective tau values."""
        # Use a non-maximally-mixed sigma so memory_alpha has visible effect
        channels = [depolarizing(0.1)] * 6
        sigma_0 = np.diag([0.9, 0.1]).astype(complex)
        rho = np.array([[0.5, 0.5], [0.5, 0.5]], dtype=complex)
        r_markov = bayesian_compose(channels, sigma_0, rho, memory_alpha=None)
        r_memory = bayesian_compose(channels, sigma_0, rho, memory_alpha=0.5)
        # They should differ (memory modifies sigma evolution)
        diff = any(
            abs(g1.tau_eff - g2.tau_eff) > 1e-10
            for g1, g2 in zip(
                r_markov.gate_results[1:], r_memory.gate_results[1:]
            )
        )
        assert diff, "memory_alpha should affect per-gate tau_eff values"


# ---------------------------------------------------------------------------
# GateResult classification
# ---------------------------------------------------------------------------


class TestClassifyGate:
    """Tests for gate classification logic."""

    def test_saturated(self):
        """tau_eff < 1e-10 -> SATURATED."""
        c = _classify_gate(0.1, 0.0, 0.5, 0.5)
        assert c == "SATURATED"

    def test_near_saturated_by_comm(self):
        """Low comm_norm -> NEAR_SATURATED."""
        c = _classify_gate(0.1, 0.05, 0.005, 0.5, comm_threshold=0.01)
        assert c == "NEAR_SATURATED"

    def test_near_saturated_by_pre_comm(self):
        """Low pre_comm_norm -> NEAR_SATURATED."""
        c = _classify_gate(0.1, 0.05, 0.5, 0.005, comm_threshold=0.01)
        assert c == "NEAR_SATURATED"

    def test_normal(self):
        """High tau_eff and comm norms -> NORMAL."""
        c = _classify_gate(0.1, 0.05, 0.5, 0.5, comm_threshold=0.01)
        assert c == "NORMAL"

    def test_classification_in_result(self):
        """Classifications appear in GateResult from bayesian_compose."""
        channels = [depolarizing(0.1)] * 3
        sigma_0 = np.eye(2, dtype=complex) / 2.0
        rho = np.array([[0.5, 0.5], [0.5, 0.5]], dtype=complex)
        result = bayesian_compose(channels, sigma_0, rho)
        for gr in result.gate_results:
            assert gr.classification in {"NORMAL", "NEAR_SATURATED", "SATURATED"}


# ---------------------------------------------------------------------------
# CompositionResult fields
# ---------------------------------------------------------------------------


class TestCompositionResult:
    """Verify CompositionResult fields are populated correctly."""

    def test_all_fields_present(self):
        """All dataclass fields are populated."""
        channels = [depolarizing(0.05)] * 2
        sigma_0 = np.eye(2, dtype=complex) / 2.0
        rho = np.array([[0.5, 0.5], [0.5, 0.5]], dtype=complex)
        result = bayesian_compose(channels, sigma_0, rho)
        assert hasattr(result, "gate_results")
        assert hasattr(result, "tau_bayesian_total")
        assert hasattr(result, "tau_multiplicative_total")
        assert hasattr(result, "improvement_percent")
        assert hasattr(result, "composition_lhs")
        assert hasattr(result, "composition_rhs")
        assert hasattr(result, "composition_holds")
        assert hasattr(result, "composition_slack")

    def test_composition_slack_nonneg(self):
        """Composition slack (RHS - LHS) >= 0 when composition holds."""
        channels = [depolarizing(0.05)] * 3
        sigma_0 = np.eye(2, dtype=complex) / 2.0
        rho = np.array([[0.5, 0.5], [0.5, 0.5]], dtype=complex)
        result = bayesian_compose(channels, sigma_0, rho)
        if result.composition_holds:
            assert result.composition_slack >= -1e-10

    def test_tau_values_nonneg(self):
        """All tau values are non-negative."""
        channels = [depolarizing(0.05)] * 3
        sigma_0 = np.eye(2, dtype=complex) / 2.0
        rho = np.array([[0.5, 0.5], [0.5, 0.5]], dtype=complex)
        result = bayesian_compose(channels, sigma_0, rho)
        assert result.tau_bayesian_total >= -1e-10
        assert result.tau_multiplicative_total >= -1e-10
        for gr in result.gate_results:
            assert gr.tau_naive >= -1e-10
            assert gr.tau_eff >= -1e-10


# ---------------------------------------------------------------------------
# Deep circuits (exercising compression path)
# ---------------------------------------------------------------------------


class TestDeepCircuits:
    """Tests for circuits deeper than 8 gates (triggers compression)."""

    def test_deep_circuit_runs(self):
        """Circuit with >8 gates runs without error."""
        channels = [depolarizing(0.02)] * 12
        sigma_0 = np.eye(2, dtype=complex) / 2.0
        rho = np.array([[0.5, 0.5], [0.5, 0.5]], dtype=complex)
        result = bayesian_compose(channels, sigma_0, rho)
        assert len(result.gate_results) == 12
        assert result.composition_holds

    def test_deep_circuit_tau_in_range(self):
        """Deep circuit tau stays in [0, 1]."""
        channels = [dephasing(0.03)] * 15
        sigma_0 = np.eye(2, dtype=complex) / 2.0
        rho = np.array([[0.5, 0.5], [0.5, 0.5]], dtype=complex)
        result = bayesian_compose(channels, sigma_0, rho)
        assert 0.0 <= result.tau_bayesian_total <= 1.0 + 1e-8
