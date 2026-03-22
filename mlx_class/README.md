# mlx-class

GPU-accelerated CMB Boltzmann solver for Apple Silicon, built on [Apple MLX](https://github.com/ml-explore/mlx).

Computes the CMB temperature (TT), polarization (EE, TE), and lensed angular power spectra, as well as the matter power spectrum P(k), by solving the coupled Einstein--Boltzmann equations with all k-modes integrated simultaneously on the GPU.

**Current status**: v0.6.0. Sync gauge solver (`solver_sync.py`) is the primary accurate backend, using scipy Radau in synchronous gauge with full photon (l=25) and neutrino (l=20) hierarchies, RECFAST recombination, and line-of-sight integration. The fast backend (`ProductionSolver`) achieves 4.9% RMS accuracy vs CLASS in ~2s on Apple M1.

## Features

- **Sync gauge solver** (`SyncSolver`): full synchronous gauge Boltzmann hierarchy with scipy Radau, accurate to ~few% vs CLASS
- **Green's function solver** (`GreensFunctionSolver`): precomputed transfer functions for instant C_l re-evaluation
- **Production solver** (`ProductionSolver`): fast conformal Newtonian backend (~2s, 4.9% RMS vs CLASS)
- **GPU-native**: all k-modes batched on Apple Silicon GPU via MLX
- **Full physics**: photon hierarchy (l=25), neutrino hierarchy (l=20), Psi != Phi
- **All spectra**: TT, TE, EE (unlensed + lensed), BB (tensor), P(k) (linear + nonlinear)
- **Line-of-sight integration**: visibility-weighted Bessel projection with SW + Doppler + ISW (early + late)
- **RECFAST recombination**: three-level atom ODE with enhanced fudge factors (~0.1% vs full RECFAST)
- **Khronon mode**: ghost-condensation Khronon dark matter with DBI sound speed
- **Dark energy**: w0/wa CPL parameterization
- **Non-flat geometry**: Omega_k support (open and closed)
- **CLASS-compatible .ini files**: drop-in parameter file support
- **MCMC**: built-in Metropolis parameter estimation
- **HaloFit**: nonlinear P(k) corrections
- **Galaxy C_l**: tomographic angular power spectra
- **SZ effect**: thermal + kinetic Sunyaev-Zeldovich
- **Tensor modes**: gravitational wave B-mode spectrum

## Installation

**Requirements**: macOS with Apple Silicon (M1/M2/M3/M4), Python >= 3.10.

```bash
# Install all dependencies
pip install mlx numpy scipy matplotlib
```

Or install the package in development mode:

```bash
# From the repository root (one level above mlx_class/)
pip install -e mlx_class/

# Or from inside the mlx_class/ directory
pip install -e .

# With plotting support
pip install -e ".[plot]"

# With development tools (pytest, etc.)
pip install -e ".[dev]"
```

## Quick start

### Recommended: CMBSolver (unified API)

```python
from mlx_class import CMBSolver

# Accurate mode: sync gauge (~420s, few% accuracy vs CLASS)
solver = CMBSolver()
result = solver.run()
print(result.ell, result.Dl_TT)  # numpy arrays
print(result.peaks)               # {'ell': [...], 'Dl': [...]}

# Fast mode: conformal Newtonian (~2s, 4.9% RMS vs CLASS)
solver = CMBSolver(backend='fast')
result = solver.run()
```

### Sync gauge solver (direct)

```python
from mlx_class import SyncSolver, SyncSolverParallel

# Sequential (N_k=300 default, ~180s)
result = SyncSolver(N_k=300, verbose=True)

# Parallel (uses all CPU cores, ~8x faster)
result = SyncSolverParallel(N_k=500, verbose=True)

# result is a dict: ell, Dl_TT, Dl_EE, Dl_TE, k_arr, bg, timing, ...
ell, Dl_TT = result['ell'], result['Dl_TT']
```

### Green's function solver (precomputed)

```python
from mlx_class import GreensFunctionSolver
from mlx_class.background import Background
import numpy as np

bg = Background(khronon=False, recombination='recfast')
bg.solve()
k_arr = np.geomspace(3e-4, 0.35, 500)

gf = GreensFunctionSolver(k_arr, bg)
gf.precompute_sources()          # one-time cost (~5 min)
gf.save('greens_fiducial.npz')   # cache to disk

result = gf.compute_cl()          # instant re-evaluation
```

### ProductionSolver (fast backend)

```python
from mlx_class import ProductionSolver

# Default mode (~1.5s cached, Peebles recombination)
solver = ProductionSolver()
result = solver.run()

# Access spectra
ell = result.ell
Dl_TT = result.Dl_TT       # TT in muK^2
Dl_TE = result.Dl_TE       # TE in muK^2
Dl_EE = result.Dl_EE       # EE in muK^2
```

### .ini file mode

```bash
# Run with CLASS-compatible .ini file
python -m mlx_class planck_bestfit.ini

# Verify .ini parsing
python -m mlx_class planck_bestfit.ini --verify
```

### Command line

```bash
# Sync gauge solver (parallel by default)
python -m mlx_class.solver_sync --N_k 500 --parallel

# Sync gauge solver (sequential)
python -m mlx_class.solver_sync --N_k 300 --sequential

# Production solver with CLASS comparison
python -m mlx_class.solver_optimized

# IMEX solver (~2.5s)
python -m mlx_class.main --implicit
```

## Performance

| Solver | Time | RMS vs CLASS | Notes |
|--------|------|-------------|-------|
| **SyncSolver** (parallel, N_k=500) | ~50s | ~few% | Full sync gauge + RECFAST |
| **SyncSolver** (sequential, N_k=300) | ~180s | ~few% | Same physics, single-threaded |
| **ProductionSolver** (default) | ~2.0s | 4.9% | Peebles + GPU Bessel + real ISW |
| IMEX | ~2.5s | ~15% | TCA only, good peaks |
| HiRes (l_gamma=25) | ~8s | ~8% | Full photon hierarchy |
| CLASS (CPU) | ~10s | -- | Reference (baseline) |

Benchmarked on Apple M1 Pro, l_max = 2500.

## Architecture

```
mlx_class/                          (v0.6.0)
  __init__.py                       -- Package init, exports CMBSolver + SyncSolver + GreensFunctionSolver
  __main__.py                       -- Entry point: python -m mlx_class [.ini or legacy]

  # --- Primary solvers ---
  solver_api.py                     -- CMBSolver: unified entry point (sync + fast backends)
  solver_sync.py                    -- SyncSolver: sync gauge + scipy Radau (accurate)
  solver_greens.py                  -- GreensFunctionSolver: precomputed transfer functions
  solver_magnus.py                  -- Magnus expansion solver (batched matrix exp)
  solver_production.py              -- ProductionSolver: fast backend (conformal Newtonian)
  solver_optimized.py               -- OptimizedSolver: production engine (4.9% RMS)

  # --- Solver internals ---
  solver_unified.py                 -- UnifiedSolver: multi-mode pipeline
  solver_accurate.py                -- AccurateBoltzmannSolver: Psi != Phi
  solver_fast.py                    -- FastBoltzmannSolver: 0.5s ODE engine

  # --- Background ---
  background.py                     -- Friedmann equation + recombination (~50ms)
  recombination.py                  -- Peebles three-level atom ODE

  # --- Perturbations (sync gauge) ---
  perturbations_sync.py             -- Sync gauge ODE RHS + gauge transform (core)
  perturbations_sync_v2.py          -- Sync gauge v2 (experimental)
  perturbations_sync_tca.py         -- Sync gauge + TCA switching

  # --- Perturbations (conformal Newtonian) ---
  perturbations.py                  -- Explicit RK4 TCA (MLX GPU)
  perturbations_implicit.py         -- IMEX solver (stiff Poisson)
  perturbations_neutrino.py         -- + massless neutrino hierarchy
  perturbations_hires.py            -- Full photon hierarchy, Strang splitting
  massive_neutrino.py               -- + massive neutrino (Fermi-Dirac)
  massive_nu_sync.py                -- Massive neutrinos in sync gauge

  # --- Spectra ---
  spectra.py                        -- C_l projection (instantaneous recombination)
  spectra_los.py                    -- Line-of-sight C_l integration
  polarization_full.py              -- EE, TE from full photon hierarchy
  matter_pk.py                      -- Linear matter power spectrum P(k)
  tensor.py                         -- Tensor (BB) perturbations (conformal Newtonian)
  tensor_sync.py                    -- Tensor perturbations (sync gauge)
  lensing.py                        -- CMB lensing (Limber + flat-sky)
  halofit.py                        -- Nonlinear P(k) (HaloFit)

  # --- Extensions ---
  dark_energy.py                    -- w0/wa dark energy (CPL)
  curvature.py                      -- Non-flat geometry (Omega_k)
  khronon.py                        -- Khronon DBI dark matter (conformal Newtonian)
  khronon_sync.py                   -- Khronon in sync gauge
  isocurvature.py                   -- Isocurvature initial conditions
  non_gaussianity.py                -- Primordial non-Gaussianity (f_NL)
  sz_effect.py                      -- Sunyaev-Zeldovich (thermal + kinetic)
  galaxy_cl.py                      -- Galaxy angular power spectrum (tomographic)

  # --- Infrastructure ---
  bessel_gpu.py                     -- GPU Bessel function table (Chebyshev)
  bessel_cache.py                   -- Disk-cached Bessel table
  ndf15.py                          -- NDF15 ODE solver
  ndf15_batched.py                  -- Batched NDF15 solver
  bdf_solver.py                     -- BDF ODE solver
  calibration.py                    -- Calibration tables
  radiation_driving.py              -- Radiation driving correction
  precision_fixes.py                -- Precision post-processing
  inifile.py                        -- CLASS-compatible .ini parser
  mcmc.py                           -- MCMC parameter estimation
  transfer.py                       -- Bessel function table (scipy fallback)
  accuracy_analysis.py              -- Accuracy analysis and diagnostics
  class_comparison.py               -- Comparison against CLASS
  main.py                           -- CLI orchestrator (legacy)

  # --- Tests ---
  tests/
    test_basic.py                   -- Comprehensive test suite (10 test classes)
    test_recombination.py           -- Recombination-specific tests
    test_sync_gauge.py              -- Sync gauge solver tests vs CLASS
  test_sync_vs_class.py             -- Sync vs CLASS comparison script

  # --- Deprecated (moved to deprecated/) ---
  deprecated/
    solver_combined.py              -- Superseded by solver_sync.py
    solver_precision.py             -- Superseded by solver_sync.py
    solver_scipy_bdf.py             -- Superseded by solver_sync.py
    solver_bdf_production.py        -- Superseded by solver_sync.py
    polarization.py                 -- Superseded by polarization_full.py
    test_combined.py                -- Superseded by test_sync_vs_class.py
    test_los.py                     -- Superseded by test_sync_vs_class.py
    test_dk_optimized.py            -- Superseded by test_sync_vs_class.py
    tools/
      optimize_dk.py                -- dk optimization (development tool)
      optimize_dk2.py               -- dk optimization v2
      optimize_dk3.py               -- dk optimization v3

  examples/
    quickstart.py                   -- Self-contained quick-start demo
  planck_bestfit.ini                -- Planck 2018 best-fit parameters

  API.md                            -- Full API reference
  ROADMAP.md                        -- Development roadmap
```

## Running tests

```bash
# Full test suite
cd mlx_class/
pytest tests/ -v

# Sync gauge tests only
pytest tests/test_sync_gauge.py -v

# Quick tests only (skip slow solver tests)
pytest tests/ -v -k "not slow"
```

Or run individual solver self-tests:

```bash
python -m mlx_class.solver_sync --N_k 500  # Sync gauge solver
python -m mlx_class.solver_optimized         # Production solver + CLASS comparison
python -m mlx_class.perturbations_implicit   # IMEX test
python -m mlx_class.perturbations_neutrino   # Neutrino test
python -m mlx_class.perturbations_hires      # HiRes convergence test
python -m mlx_class.polarization_full        # Full EE polarization test
python -m mlx_class.matter_pk                # P(k) test
```

## Theory background

This solver implements both standard LCDM and the Khronon dark matter model:

- **Sigma = 2 ln Q**: the relative entropy between spacetime and matter degrees of freedom, where Q = 1/N_phi is the inverse Khronon lapse
- **Ghost condensation**: K'(Q_0) = 0 at the Khronon background, giving c_s^2 = 0 exactly (dust-like perturbations, consistent with CMB and LSS)
- **DBI correction**: at next-to-leading order, the DBI kinetic structure generates c_s^2(k) = alpha_DBI * k^2/(k^2 + k_J^2), a k-dependent sound speed that modifies the CMB at high multipoles

For details, see: S.-K. Huang, "Petz Recovery Map as Retrodiction" (2026).

## License

MIT License. Copyright (c) 2026 Sheng-Kai Huang.

## Citation

If you use mlx-class in your research, please cite:

```bibtex
@software{huang2026mlxclass,
    author = {Huang, Sheng-Kai},
    title = {mlx-class: GPU-accelerated CMB Boltzmann solver for Apple Silicon},
    year = {2026},
    url = {https://github.com/akaiHuang/petz-recovery-unification},
}
```
