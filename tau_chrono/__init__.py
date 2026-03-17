"""
tau-chrono: Bayesian Noise Tracker for Quantum Circuits.

Track noise evolution through your quantum circuit.
Get per-gate tau estimates 50-100% more accurate than independent noise models.
"""

__version__ = "0.2.0"

# Core API -- available at top level
from .channels import (
    amplitude_damping,
    depolarizing,
    dephasing,
    verify_cptp,
    two_qubit_depolarizing,
    correlated_dephasing,
    local_noise,
    cnot_error,
    crosstalk_dephasing,
    amplitude_damping_2q,
    unitary_channel,
    entanglement_fidelity,
    diamond_norm_approx,
    coherent_incoherent_decomposition,
)
from .petz import (
    apply_channel,
    adjoint_channel,
    petz_recovery_map,
    apply_petz_recovery,
    fidelity,
    tau_parameter,
    relative_entropy,
    delta_D,
)
from .bayesian import (
    bayesian_compose,
    compose_kraus,
    compose_kraus_compressed,
    GateResult,
    CompositionResult,
)
from .utils import (
    matrix_sqrt,
    matrix_inv_sqrt,
    matrix_log,
    trace_norm,
    commutator_norm,
)
from .tomography import (
    build_tomography_circuits_1q,
    counts_to_choi_1q,
    choi_to_kraus,
    simulate_process_tomography_1q,
)
from ._types import DensityMatrix, KrausList, SuperOp

# Visualization (optional -- requires matplotlib)
try:
    from .viz import (
        plot_tau_heatmap,
        plot_depth_scaling,
        plot_composition_inequality,
        plot_sigma_evolution,
        plot_ground_truth_comparison,
        plot_channel_diagnostics,
        plot_summary_dashboard,
    )
    _HAS_VIZ = True
except ImportError:
    _HAS_VIZ = False

__all__ = [
    # channels
    "amplitude_damping",
    "depolarizing",
    "dephasing",
    "verify_cptp",
    "two_qubit_depolarizing",
    "correlated_dephasing",
    "local_noise",
    "cnot_error",
    "crosstalk_dephasing",
    "amplitude_damping_2q",
    "unitary_channel",
    "entanglement_fidelity",
    "diamond_norm_approx",
    "coherent_incoherent_decomposition",
    # petz
    "apply_channel",
    "adjoint_channel",
    "petz_recovery_map",
    "apply_petz_recovery",
    "fidelity",
    "tau_parameter",
    "relative_entropy",
    "delta_D",
    # bayesian
    "bayesian_compose",
    "compose_kraus",
    "compose_kraus_compressed",
    "GateResult",
    "CompositionResult",
    # utils
    "matrix_sqrt",
    "matrix_inv_sqrt",
    "matrix_log",
    "trace_norm",
    "commutator_norm",
    # viz (optional)
    "plot_tau_heatmap",
    "plot_depth_scaling",
    "plot_composition_inequality",
    "plot_sigma_evolution",
    "plot_ground_truth_comparison",
    "plot_channel_diagnostics",
    "plot_summary_dashboard",
]
