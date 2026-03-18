"""
sigma_diffusion -- Core Sigma computation for diffusion processes.

Provides tools to track information-theoretic divergence (Sigma) across
diffusion steps, derive adaptive scheduling from Sigma profiles, compute
image quality metrics, and generate interactive HTML dashboards.

Theoretical basis:
    Sigma_t = D_KL(p_{t-1} || p_t)
    Petz fidelity bound:  F >= exp(-Sigma / 2)

See: Huang (2026), "Petz Recovery Map as a Universal Retrodiction Functor."
"""

from sigma_diffusion.sigma_analyzer import (
    SigmaAnalyzer,
    SigmaResult,
    SigmaProfile,
    SigmaTracker,
    StepDiagnostics,
    kl_divergence_gaussian,
    kl_divergence_discrete,
)
from sigma_diffusion.adaptive_scheduler import (
    AdaptiveScheduler,
    AdaptiveSchedule,
    compute_skip_mask,
    compute_adaptive_timesteps,
    compute_early_stop,
)

__all__ = [
    "SigmaAnalyzer",
    "SigmaResult",
    "SigmaProfile",
    "SigmaTracker",
    "StepDiagnostics",
    "AdaptiveScheduler",
    "AdaptiveSchedule",
    "kl_divergence_gaussian",
    "kl_divergence_discrete",
    "compute_skip_mask",
    "compute_adaptive_timesteps",
    "compute_early_stop",
]

__version__ = "0.1.0"
