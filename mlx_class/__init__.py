"""
mlx_class -- GPU-accelerated CMB Boltzmann solver using Apple MLX.

All k-modes integrated simultaneously on the GPU via batched ODE solvers.
Achieves ~2s total runtime on Apple M1 with 4.9% RMS accuracy vs CLASS.

Primary API:
    from mlx_class import CMBSolver
    solver = CMBSolver()                    # sync gauge (accurate, ~420s)
    solver = CMBSolver(backend='fast')      # conformal Newtonian (~2s)
    result = solver.run()
    # result.ell, result.Dl_TT, result.Dl_TE, result.Dl_EE, result.peaks, result.timing

Sync gauge solver (direct):
    from mlx_class import SyncSolver
    result = SyncSolver(N_k=500, parallel=True)
    # Returns dict: ell, Dl_TT, Dl_EE, Dl_TE, ...

Green's function solver (precomputed transfer functions):
    from mlx_class import GreensFunctionSolver

Legacy API (backward compatible):
    from mlx_class import ProductionSolver
    solver = ProductionSolver()
    result = solver.run()

Architecture:
    solver_api.py              -- CMBSolver: unified entry point (sync + fast backends)
    solver_sync.py             -- Synchronous gauge solver (scipy Radau, accurate)
    solver_greens.py           -- Green's function solver (precomputed transfers)
    solver_magnus.py           -- Magnus expansion solver (batched matrix exp)
    solver_production.py       -- ProductionSolver: fast backend (conformal Newtonian)
    solver_optimized.py        -- OptimizedSolver: production engine (4.9% RMS vs CLASS)
    solver_unified.py          -- UnifiedSolver: multi-mode pipeline (TT/TE/EE/BB/P(k)/lensing)
    solver_accurate.py         -- AccurateBoltzmannSolver: correct Psi != Phi physics
    solver_fast.py             -- FastBoltzmannSolver: 0.5s ODE via mx.compile
    background.py              -- Friedmann equation, recombination (numpy)
    perturbations_sync.py      -- Synchronous gauge ODE RHS + gauge transform
    perturbations.py           -- TCA ODE solver, batched over k (MLX GPU)
    perturbations_implicit.py  -- IMEX Boltzmann solver (stiff Poisson)
    perturbations_neutrino.py  -- + massless neutrino hierarchy
    perturbations_hires.py     -- Full photon hierarchy, Strang splitting
    massive_neutrino.py        -- + massive neutrino (Fermi-Dirac, per-q-bin)
    polarization_full.py       -- EE, TE from HiRes Theta_2
    lensing.py                 -- CMB lensing (Limber + flat-sky)
    halofit.py                 -- Nonlinear P(k) (HaloFit)
    dark_energy.py             -- w0/wa dark energy extension
    curvature.py               -- Non-flat spatial geometry (Omega_k)
    inifile.py                 -- CLASS-compatible .ini file parser
    calibration.py             -- Calibration tables and comparison tools
    bessel_gpu.py              -- GPU Bessel function table (Chebyshev)
    bessel_cache.py            -- Disk-cached Bessel table
    khronon.py                 -- Khronon DBI dark matter solver
    tensor.py                  -- Tensor (gravitational wave) perturbations
    galaxy_cl.py               -- Galaxy angular power spectrum (tomographic)
    isocurvature.py            -- Isocurvature initial conditions
    non_gaussianity.py         -- Primordial non-Gaussianity (f_NL)
    sz_effect.py               -- Sunyaev-Zeldovich effect (thermal + kinetic)
    mcmc.py                    -- MCMC parameter estimation
    radiation_driving.py       -- Radiation driving correction
    precision_fixes.py         -- Precision post-processing corrections
    recombination.py           -- Peebles three-level atom recombination
    matter_pk.py               -- Linear matter power spectrum P(k)
    spectra.py                 -- C_l projection (instantaneous recombination)
    spectra_los.py             -- Line-of-sight C_l integration
    transfer.py                -- Bessel function table (scipy fallback)
    accuracy_analysis.py       -- Accuracy analysis and diagnostics
    class_comparison.py        -- Comparison against CLASS
    main.py                    -- CLI orchestrator (legacy)
    __main__.py                -- Entry point: python -m mlx_class

Deprecated (moved to deprecated/):
    solver_combined.py         -- Use solver_production.py
    solver_precision.py        -- Use solver_production.py
    solver_scipy_bdf.py        -- Use solver_sync.py
    solver_bdf_production.py   -- Use solver_sync.py
    polarization.py            -- Use polarization_full.py
    test_combined.py           -- Use tests/test_basic.py
    test_los.py                -- Use tests/test_basic.py
    test_dk_optimized.py       -- Use test_sync_vs_class.py
    optimize_dk*.py            -- Moved to deprecated/tools/

Author: Sheng-Kai Huang, 2026
Theory: Sigma = 2 ln Q, ghost-condensation Khronon with c_s^2 = 0
"""
__version__ = "0.6.0"

from .solver_api import CMBSolver, CMBResult

# Sync gauge solver (direct function call)
from .solver_sync import run_sync_solver as SyncSolver
from .solver_sync import run_sync_solver_parallel as SyncSolverParallel

# Green's function solver
from .solver_greens import GreensFunctionSolver

# Legacy API (backward compatible)
from .solver_production import ProductionSolver, ProductionResult

__all__ = [
    'CMBSolver',
    'CMBResult',
    'SyncSolver',
    'SyncSolverParallel',
    'GreensFunctionSolver',
    'ProductionSolver',
    'ProductionResult',
    '__version__',
]
