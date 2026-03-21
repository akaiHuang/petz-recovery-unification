"""
mlx_class -- GPU-accelerated CMB Boltzmann solver using Apple MLX.

All k-modes integrated simultaneously on the GPU via batched ODE solvers.
Achieves ~2s total runtime on Apple M1 with 4.9% RMS accuracy vs CLASS.

Primary API:
    from mlx_class import ProductionSolver
    solver = ProductionSolver()
    result = solver.run()
    # result.ell, result.Dl_TT, result.Dl_TE, result.Dl_EE, result.peaks, result.timing

Architecture (44 modules):
    solver_production.py       -- ProductionSolver: recommended entry point
    solver_optimized.py        -- OptimizedSolver: production engine (4.9% RMS vs CLASS)
    solver_unified.py          -- UnifiedSolver: multi-mode pipeline (TT/TE/EE/BB/P(k)/lensing)
    solver_accurate.py         -- AccurateBoltzmannSolver: correct Psi != Phi physics
    solver_fast.py             -- FastBoltzmannSolver: 0.5s ODE via mx.compile
    background.py              -- Friedmann equation, recombination (numpy)
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

Deprecated (kept for backward compatibility):
    solver_combined.py         -- Use solver_production.py
    solver_precision.py        -- Use solver_production.py
    polarization.py            -- Use polarization_full.py
    test_combined.py           -- Use tests/test_basic.py
    test_los.py                -- Use tests/test_basic.py

Author: Sheng-Kai Huang, 2026
Theory: Sigma = 2 ln Q, ghost-condensation Khronon with c_s^2 = 0
"""
__version__ = "0.4.0"

from .solver_production import ProductionSolver, ProductionResult

__all__ = [
    'ProductionSolver',
    'ProductionResult',
    '__version__',
]
