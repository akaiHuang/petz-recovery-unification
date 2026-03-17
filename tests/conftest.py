"""Common fixtures for tau-chrono test suite."""

import numpy as np
import pytest


@pytest.fixture
def rho_plus():
    """|+><+| = [[0.5, 0.5], [0.5, 0.5]]"""
    return np.array([[0.5, 0.5], [0.5, 0.5]], dtype=complex)


@pytest.fixture
def rho_zero():
    """|0><0| = [[1, 0], [0, 0]]"""
    return np.array([[1.0, 0.0], [0.0, 0.0]], dtype=complex)


@pytest.fixture
def rho_one():
    """|1><1| = [[0, 0], [0, 1]]"""
    return np.array([[0.0, 0.0], [0.0, 1.0]], dtype=complex)


@pytest.fixture
def sigma_mixed():
    """diag(0.8, 0.2) -- a non-maximally-mixed state."""
    return np.diag(np.array([0.8, 0.2], dtype=complex))


@pytest.fixture
def sigma_maximally_mixed():
    """I/2 -- the maximally mixed qubit state."""
    return np.eye(2, dtype=complex) / 2.0


@pytest.fixture
def identity_channel():
    """Identity channel: single Kraus operator = I_2."""
    return [np.eye(2, dtype=complex)]


@pytest.fixture
def rho_plus_i():
    """|+i><+i| = [[0.5, -0.5j], [0.5j, 0.5]]"""
    return np.array([[0.5, -0.5j], [0.5j, 0.5]], dtype=complex)


@pytest.fixture
def bell_state():
    """Bell state |Phi+> = (|00> + |11>)/sqrt(2), as 4x4 density matrix."""
    psi = np.array([1, 0, 0, 1], dtype=complex) / np.sqrt(2)
    return np.outer(psi, psi.conj())


# ---------------------------------------------------------------------------
# Helper for generating random density matrices
# ---------------------------------------------------------------------------

def random_density_matrix(d: int, rng: np.random.Generator = None) -> np.ndarray:
    """Generate a random density matrix of dimension d (Hilbert-Schmidt measure)."""
    if rng is None:
        rng = np.random.default_rng(42)
    G = rng.standard_normal((d, d)) + 1j * rng.standard_normal((d, d))
    rho = G @ G.conj().T
    rho /= np.trace(rho).real
    return rho


def random_psd_matrix(d: int, rng: np.random.Generator = None) -> np.ndarray:
    """Generate a random positive semidefinite matrix (not necessarily trace 1)."""
    if rng is None:
        rng = np.random.default_rng(123)
    G = rng.standard_normal((d, d)) + 1j * rng.standard_normal((d, d))
    return G @ G.conj().T


def random_hermitian_matrix(d: int, rng: np.random.Generator = None) -> np.ndarray:
    """Generate a random Hermitian matrix (may have negative eigenvalues)."""
    if rng is None:
        rng = np.random.default_rng(99)
    H = rng.standard_normal((d, d)) + 1j * rng.standard_normal((d, d))
    return 0.5 * (H + H.conj().T)
