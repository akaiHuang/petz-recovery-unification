# mlx_class API Reference

Complete reference for all public classes and functions in `mlx_class` (v0.4.0, 44 modules).

All units are **Mpc** (comoving) unless otherwise noted. Wavenumbers `k` are in Mpc^{-1}. Angular power spectra `D_l` are in muK^2.

---

## Module Index (44 modules)

One-line description of every module in the package.

### Production solver (recommended entry point)

| Module | Description |
|--------|-------------|
| `solver_production.py` | `ProductionSolver` -- recommended API for all CMB computations (~2s, 4.9% RMS vs CLASS) |
| `solver_optimized.py` | `OptimizedSolver` -- production engine: FastBoltzmann ODE + GPU Bessel + real ISW + ell-remap |

### Solver internals

| Module | Description |
|--------|-------------|
| `solver_unified.py` | `UnifiedSolver` -- multi-mode pipeline dispatching to all solver tiers (TT/TE/EE/BB/P(k)/lensing) |
| `solver_accurate.py` | `AccurateBoltzmannSolver` -- correct Psi != Phi physics + extended integration beyond tau_rec |
| `solver_fast.py` | `FastBoltzmannSolver` -- 0.5s ODE engine via mx.compile + vectorized Strang splitting |

### Background cosmology

| Module | Description |
|--------|-------------|
| `background.py` | `Background` -- Friedmann equation, conformal time, sound horizon, recombination (~50ms) |
| `recombination.py` | `PeeblesRecombination` -- Peebles three-level atom ODE for x_e(z) (~0.1% vs RECFAST) |

### Perturbation solvers

| Module | Description |
|--------|-------------|
| `perturbations.py` | `AnalyticTransfer`, `BatchedBoltzmannSolver` -- analytic WKB + explicit RK4 TCA (MLX GPU) |
| `perturbations_implicit.py` | `ImplicitBoltzmannSolver` -- IMEX solver with exponential integrator for stiff Poisson |
| `perturbations_neutrino.py` | `NeutrinoBoltzmannSolver` -- IMEX + massless neutrino hierarchy (N_eff=3.046, Psi != Phi) |
| `perturbations_hires.py` | `HiResBoltzmannSolver` -- full photon hierarchy (l_gamma_max=25) + Strang splitting |
| `massive_neutrino.py` | `MassiveNeutrinoBoltzmannSolver` -- massive neutrino Fermi-Dirac per-q-bin integration |

### Spectra computation

| Module | Description |
|--------|-------------|
| `spectra.py` | `compute_cl` -- C_l projection via instantaneous recombination + Bessel j_l^2 |
| `spectra_los.py` | `compute_cl_los` -- line-of-sight C_l with visibility weighting, Doppler, ISW |
| `polarization_full.py` | `run_polarization_pipeline` -- EE and TE from full photon hierarchy Theta_2 |
| `matter_pk.py` | `compute_matter_pk` -- linear matter power spectrum P(k) with growth factor |
| `tensor.py` | `run_tensor_pipeline` -- tensor (gravitational wave) B-mode spectrum |
| `lensing.py` | `apply_lensing` -- CMB lensing C_l^phi_phi (Limber + flat-sky convolution) |
| `halofit.py` | `compute_nonlinear_pk` -- nonlinear P(k) via HaloFit (Takahashi et al. 2012) |

### Extensions

| Module | Description |
|--------|-------------|
| `dark_energy.py` | `BackgroundDE` -- w0/wa dark energy (CPL parameterization) |
| `curvature.py` | `BackgroundCurved` -- non-flat spatial geometry (Omega_k != 0) |
| `khronon.py` | `KhrononBoltzmannSolver` -- Khronon DBI dark matter with ghost condensation |
| `isocurvature.py` | Isocurvature initial conditions (CDI, NID, NIV modes) |
| `non_gaussianity.py` | Primordial non-Gaussianity (f_NL local, equilateral, orthogonal) |
| `sz_effect.py` | Sunyaev-Zeldovich effect (thermal tSZ + kinetic kSZ spectra) |
| `galaxy_cl.py` | Galaxy angular power spectrum (tomographic bins, magnification bias) |

### Infrastructure

| Module | Description |
|--------|-------------|
| `bessel_gpu.py` | `BesselTable` -- GPU spherical Bessel j_l via segmented Chebyshev interpolation |
| `bessel_cache.py` | `CachedBesselTable` -- disk-cached wrapper for BesselTable (~0.02s per load) |
| `calibration.py` | Calibration tables and comparison tools (CLASS reference spectra) |
| `radiation_driving.py` | Radiation driving correction (WKB amplification at k >> k_eq) |
| `precision_fixes.py` | Precision post-processing corrections (envelope, ell-remap) |
| `inifile.py` | `load_ini`, `create_solver`, `run_from_ini` -- CLASS-compatible .ini file parser |
| `mcmc.py` | `MetropolisMCMC` -- Metropolis-Hastings MCMC parameter estimation |
| `transfer.py` | Bessel function table via scipy (fallback for non-MLX systems) |
| `accuracy_analysis.py` | Accuracy analysis, RMS residuals, peak comparison diagnostics |
| `class_comparison.py` | Side-by-side comparison against CLASS (requires classy) |
| `main.py` | CLI orchestrator for legacy command-line interface |
| `__main__.py` | Entry point for `python -m mlx_class` (.ini or legacy mode) |
| `__init__.py` | Package init, exports `ProductionSolver`, version 0.4.0 |

### Tests and examples

| Module | Description |
|--------|-------------|
| `tests/test_basic.py` | Comprehensive test suite (10 test classes, all output types) |
| `tests/test_recombination.py` | Recombination-specific tests (tanh vs Peebles vs RECFAST) |
| `examples/quickstart.py` | Self-contained quick-start demo |
| `planck_bestfit.ini` | Planck 2018 best-fit LCDM parameters (CLASS-compatible) |

### Deprecated (kept for backward compatibility)

| Module | Replacement | Description |
|--------|-------------|-------------|
| `solver_combined.py` | `solver_production.py` | Old combined pipeline (Peebles + neutrinos + ODE + LOS) |
| `solver_precision.py` | `solver_production.py` | Old precision solver (TCA + enhanced LOS) |
| `polarization.py` | `polarization_full.py` | Basic EE from TCA-estimated Theta_2 (~100x too low) |
| `test_combined.py` | `tests/test_basic.py` | Old combined neutrino + ISW test |
| `test_los.py` | `tests/test_basic.py` | Old LOS integration test |

---

## Table of Contents (detailed API)

- [ProductionSolver](#productionsolver)
- [Background](#background)
- [Solvers](#solvers)
  - [BatchedBoltzmannSolver](#batchedboltzmannsolver)
  - [ImplicitBoltzmannSolver](#implicitboltzmannsolver)
  - [NeutrinoBoltzmannSolver](#neutrinoboltzmannsolver)
  - [HiResBoltzmannSolver](#hiresboltzmannsolver)
  - [KhrononBoltzmannSolver](#khrononboltzmannsolver)
- [Analytic Transfer](#analytictransfer)
- [Spectra](#spectra)
  - [compute_cl](#compute_cl)
  - [compute_cl_los](#compute_cl_los)
  - [find_peaks_in_spectrum](#find_peaks_in_spectrum)
- [Matter Power Spectrum](#matter-power-spectrum)
  - [compute_matter_pk](#compute_matter_pk)
- [Polarization](#polarization)
  - [compute_cl_ee](#compute_cl_ee)
- [Recombination](#recombination)
  - [PeeblesRecombination](#peeblesrecombination)
- [Constants](#constants)

---

## ProductionSolver

```python
from mlx_class import ProductionSolver
# or: from mlx_class.solver_production import ProductionSolver
```

### `ProductionSolver(h=0.6736, omega_b=0.02237, omega_cdm=0.12, ...)`

The recommended entry point for CMB TT computations. Combines accurate physics (Psi != Phi, Poisson constraint reset) with fast execution (mx.compile, cached Bessel, D(k) correction).

**Parameters**

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `h` | `float` | `0.6736` | Dimensionless Hubble parameter |
| `omega_b` | `float` | `0.02237` | Physical baryon density |
| `omega_cdm` | `float` | `0.12` | Physical CDM density |
| `l_gamma_max` | `int` | `25` | Maximum photon multipole |
| `N_k` | `int` | `500` | Number of k-modes |
| `ell_max` | `int` | `2500` | Maximum multipole |
| `recombination` | `str` | `'peebles'` | Recombination method ('peebles' or 'tanh') |
| `apply_D_correction` | `bool` | `True` | Apply D(k) radiation driving correction |
| `verbose` | `bool` | `True` | Print progress information |

**Instance methods**

| Method | Returns | Description |
|--------|---------|-------------|
| `run()` | `ProductionResult` | Run the full pipeline and return all spectra. |

**Result: `ProductionResult`**

| Attribute | Type | Description |
|-----------|------|-------------|
| `ell` | `ndarray` | Multipole array |
| `Dl_TT` | `ndarray` | TT spectrum in muK^2 |
| `Dl_TE` | `ndarray` or `None` | TE cross-spectrum (placeholder) |
| `Dl_EE` | `ndarray` or `None` | EE spectrum (placeholder) |
| `Cl_TT` | `ndarray` | Raw dimensionless C_l |
| `peaks` | `dict` | `{'ell': [...], 'Dl': [...]}` |
| `timing` | `dict` | `{'ode', 'bessel', 'los', 'total'}` |
| `k_arr` | `ndarray` | Wavenumber grid used |

**Example**

```python
from mlx_class import ProductionSolver

solver = ProductionSolver()
result = solver.run()

# TT spectrum and peaks
print(f"First peak: l={result.peaks['ell'][0]}")
print(f"TT max: {result.Dl_TT.max():.0f} muK^2")
print(f"Total time: {result.timing['total']:.1f}s")
```

For full TT + TE + EE + lensed spectra, use `OptimizedSolver` from `solver_optimized.py`.

---

## Background

```python
from mlx_class.background import Background
```

### `Background(N_points=50000, a_min=1e-7, a_max=1.0, khronon=False, recombination='tanh')`

Solve the Friedmann equation and recombination history.

**Parameters**

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `N_points` | `int` | `50000` | Number of grid points in scale factor `a`. |
| `a_min` | `float` | `1e-7` | Minimum scale factor (deep radiation era). |
| `a_max` | `float` | `1.0` | Maximum scale factor (today). |
| `khronon` | `bool` | `False` | If `True`, replace CDM with Khronon dark matter (`Omega_K = 0.265`). If `False`, standard LCDM. |
| `recombination` | `str` | `'tanh'` | Recombination model. `'tanh'`: fast tanh fit (~5% accuracy). `'peebles'`: Peebles three-level atom ODE (~0.1% accuracy). |

**Methods**

| Method | Returns | Description |
|--------|---------|-------------|
| `solve()` | `self` | Compute conformal time, sound horizon, and recombination. Must be called before using the object. |
| `calH_at_tau(tau_arr)` | `ndarray` | Conformal Hubble `aH` at given conformal times. |
| `R_at_tau(tau_arr)` | `ndarray` | Baryon loading `R = 3 rho_b / (4 rho_gamma)` at given conformal times. |
| `kappa_dot_at_tau(tau_arr)` | `ndarray` | Thomson scattering rate (negative by convention). |
| `visibility_at_tau(tau_arr)` | `ndarray` | Visibility function `g(tau) = |kappa_dot| * exp(-kappa)`. |
| `kappa_at_tau(tau_arr)` | `ndarray` | Optical depth `kappa(tau)`. |
| `a_at_tau(tau_arr)` | `ndarray` | Scale factor at given conformal times. |

**Key Attributes** (available after calling `solve()`)

| Attribute | Type | Description |
|-----------|------|-------------|
| `tau_0` | `float` | Conformal time today (Mpc). |
| `tau_rec` | `float` | Conformal time at recombination (Mpc). |
| `D_A` | `float` | Comoving angular diameter distance to LSS (Mpc). |
| `r_s` | `float` | Sound horizon at recombination (Mpc). |
| `R_rec` | `float` | Baryon loading at recombination. |
| `k_D` | `float` | Silk damping scale (Mpc^{-1}). |
| `Omega_cdm` | `float` | CDM (or Khronon) density parameter. |
| `a_grid` | `ndarray` | Scale factor grid. |
| `tau_grid` | `ndarray` | Conformal time grid. |
| `visibility_grid` | `ndarray` | Visibility function on the grid. |
| `x_e_grid` | `ndarray` | Free electron fraction on the grid. |

**Example**

```python
from mlx_class.background import Background

bg = Background(khronon=False, recombination='tanh')
bg.solve()

print(f"Sound horizon: {bg.r_s:.1f} Mpc")
print(f"First peak estimate: l ~ {3.14159 * bg.D_A / bg.r_s:.0f}")
```

---

## Solvers

All solvers share the same pattern:
1. Create with a `Background` object and a k-array.
2. Call `.solve()` to integrate the Boltzmann equations on the GPU.
3. Extract source functions from the returned result object.

### `BatchedBoltzmannSolver`

```python
from mlx_class.perturbations import BatchedBoltzmannSolver
```

#### `BatchedBoltzmannSolver(bg, k_arr_Mpc, l_max=12, mode='fast')`

Explicit RK4 Boltzmann solver in the tight-coupling approximation (TCA). Integrates to `0.85 * tau_rec`. All k-modes run simultaneously on the GPU.

**Parameters**

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `bg` | `Background` | -- | Solved background cosmology. |
| `k_arr_Mpc` | `array` | -- | Wavenumber array in Mpc^{-1}. |
| `l_max` | `int` | `12` | (Unused in current implementation; reserved for future multipole expansion.) |
| `mode` | `str` | `'fast'` | (Unused in current implementation.) |

**Methods**

| Method | Returns | Description |
|--------|---------|-------------|
| `solve()` | `PerturbationResult` | Integrate to `0.85 * tau_rec` and return result. |

**Result: `PerturbationResult`**

| Method | Returns | Description |
|--------|---------|-------------|
| `source_at_recombination()` | `(Theta_0, Phi, v_b)` | Source functions at `0.85 * tau_rec` with Silk damping applied. All arrays have shape `(N_k,)`. |

**Example**

```python
import numpy as np
from mlx_class.background import Background
from mlx_class.perturbations import BatchedBoltzmannSolver

bg = Background().solve()
k_arr = np.geomspace(5e-4, 0.2, 300).astype(np.float32)

solver = BatchedBoltzmannSolver(bg, k_arr)
result = solver.solve()

Theta_0, Phi, v_b = result.source_at_recombination()
source = Theta_0 + Phi  # Sachs-Wolfe monopole
```

---

### `ImplicitBoltzmannSolver`

```python
from mlx_class.perturbations_implicit import ImplicitBoltzmannSolver
```

#### `ImplicitBoltzmannSolver(bg, k_arr_Mpc)`

IMEX (Implicit-Explicit) Boltzmann solver. Uses an exponential integrator for the stiff Poisson equation, allowing stable integration all the way to `tau_rec` for any k.

**Parameters**

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `bg` | `Background` | -- | Solved background cosmology. |
| `k_arr_Mpc` | `array` | -- | Wavenumber array in Mpc^{-1}. |

**Methods**

| Method | Returns | Description |
|--------|---------|-------------|
| `solve()` | `ImplicitResult` | Integrate to `tau_rec` and return result. |

**Result: `ImplicitResult`**

| Method | Returns | Description |
|--------|---------|-------------|
| `source_at_recombination()` | `(Theta_0, Phi, v_b)` | Source functions at `tau_rec` with Silk damping. Arrays shape `(N_k,)`. |

| Attribute | Type | Description |
|-----------|------|-------------|
| `y_tca` | `mx.array` | Raw state vector, shape `(N_k, 6)`. |
| `k_arr` | `ndarray` | Wavenumber array. |
| `bg` | `Background` | Background object. |

**Example**

```python
import numpy as np
from mlx_class.background import Background
from mlx_class.perturbations_implicit import ImplicitBoltzmannSolver

bg = Background().solve()
k_arr = np.geomspace(5e-5, 0.35, 500).astype(np.float32)

solver = ImplicitBoltzmannSolver(bg, k_arr)
result = solver.solve()

Theta_0, Phi, v_b = result.source_at_recombination()
```

---

### `NeutrinoBoltzmannSolver`

```python
from mlx_class.perturbations_neutrino import NeutrinoBoltzmannSolver
```

#### `NeutrinoBoltzmannSolver(bg, k_arr_Mpc, l_nu_max=20)`

IMEX solver with massless neutrino Boltzmann hierarchy. Includes `N_eff = 3.046` neutrinos, separate `Omega_gamma` and `Omega_nu`, and the anisotropic stress effect `Psi != Phi`.

**Parameters**

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `bg` | `Background` | -- | Solved background cosmology. |
| `k_arr_Mpc` | `array` | -- | Wavenumber array in Mpc^{-1}. |
| `l_nu_max` | `int` | `20` | Maximum neutrino multipole. |

**Methods**

| Method | Returns | Description |
|--------|---------|-------------|
| `solve()` | `NeutrinoResult` | Integrate to `tau_rec` and return result. |

**Result: `NeutrinoResult`**

| Method | Returns | Description |
|--------|---------|-------------|
| `source_at_recombination()` | `(Theta_0, Phi, Psi, v_b, N_0, N_2)` | Source functions with Silk damping. `Psi != Phi` due to neutrino anisotropic stress. |
| `sw_source()` | `ndarray` | Sachs-Wolfe source `Theta_0 + Psi` (silk-damped). |
| `phi_psi_ratio()` | `ndarray` | `Psi/Phi` at `tau_rec` (no Silk damping). |

**Example**

```python
import numpy as np
from mlx_class.background import Background
from mlx_class.perturbations_neutrino import NeutrinoBoltzmannSolver

bg = Background().solve()
k_arr = np.geomspace(5e-4, 0.35, 500).astype(np.float32)

solver = NeutrinoBoltzmannSolver(bg, k_arr, l_nu_max=20)
result = solver.solve()

Theta_0, Phi, Psi, v_b, N_0, N_2 = result.source_at_recombination()
source = Theta_0 + Psi   # SW source with neutrino correction
```

---

### `HiResBoltzmannSolver`

```python
from mlx_class.perturbations_hires import HiResBoltzmannSolver
```

#### `HiResBoltzmannSolver(bg, k_arr_Mpc, l_gamma_max=25, l_nu_max=20, tca_threshold=50.0)`

High-resolution photon Boltzmann hierarchy solver with Strang operator splitting. Tracks photon multipoles up to `l_gamma_max` (replacing the TCA for Theta_0, Theta_1 only). Includes neutrinos, anisotropic stress from both photons and neutrinos, and a separate baryon velocity `v_b`.

Uses a two-phase strategy: TCA at early times (when Thomson scattering is fast), then switches to the full hierarchy when `|kappa_dot| / k_max` drops below the threshold.

**Parameters**

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `bg` | `Background` | -- | Solved background cosmology. |
| `k_arr_Mpc` | `array` | -- | Wavenumber array in Mpc^{-1}. |
| `l_gamma_max` | `int` | `25` | Maximum photon multipole. Higher values (50, 100) improve damping tail accuracy. |
| `l_nu_max` | `int` | `20` | Maximum neutrino multipole. |
| `tca_threshold` | `float` | `50.0` | Switch from TCA to full hierarchy when `|kappa_dot| / k_max < threshold`. |

**Methods**

| Method | Returns | Description |
|--------|---------|-------------|
| `solve()` | `HiResResult` | Integrate to `tau_rec` and return result. |

**Result: `HiResResult`**

| Method | Returns | Description |
|--------|---------|-------------|
| `source_at_recombination()` | `(Theta_0, Phi, Psi, v_b, Theta_2)` | Source functions with Silk damping. Includes photon quadrupole `Theta_2` (exact, not estimated). |
| `sw_source()` | `ndarray` | Sachs-Wolfe source `Theta_0 + Psi`. |
| `photon_multipoles()` | `dict` | All photon multipoles `{l: Theta_l}` at `tau_rec` (no damping). |
| `phi_psi_ratio()` | `ndarray` | `Psi/Phi` at `tau_rec`. |

**Example**

```python
import numpy as np
from mlx_class.background import Background
from mlx_class.perturbations_hires import HiResBoltzmannSolver

bg = Background(recombination='peebles').solve()
k_arr = np.geomspace(5e-4, 0.35, 500).astype(np.float32)

solver = HiResBoltzmannSolver(bg, k_arr, l_gamma_max=25)
result = solver.solve()

Theta_0, Phi, Psi, v_b, Theta_2 = result.source_at_recombination()
```

---

### `KhrononBoltzmannSolver`

```python
from mlx_class.khronon import KhrononBoltzmannSolver
```

#### `KhrononBoltzmannSolver(bg, k_arr_Mpc, alpha_dbi=0.5)`

IMEX Boltzmann solver for the Khronon dark matter model with DBI effective sound speed. Implements the full Generalized Dark Matter (GDM) treatment: modified Euler, continuity, and Poisson equations with k-dependent sound speed.

At ghost condensation (`K'(Q_0) = 0`), the background-level sound speed is exactly zero (`c_s^2 = 0`), but the DBI kinetic structure generates a k-dependent effective sound speed at next-to-leading order.

**Parameters**

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `bg` | `Background` | -- | Solved background (typically with `khronon=True`). |
| `k_arr_Mpc` | `array` | -- | Wavenumber array in Mpc^{-1}. |
| `alpha_dbi` | `float` | `0.5` | DBI strength parameter. Default `Q_0/2 = 0.5` is the ghost condensation bare value (upper bound). |

**Methods**

| Method | Returns | Description |
|--------|---------|-------------|
| `solve()` | `KhrononResult` | Integrate to `tau_rec` and return result. |

**Result: `KhrononResult`** (inherits from `ImplicitResult`)

| Method | Returns | Description |
|--------|---------|-------------|
| `source_at_recombination()` | `(Theta_0, Phi, v_b)` | Same interface as `ImplicitResult`. |
| `effective_cs2_at_rec()` | `ndarray` | The effective `c_s^2(k)` at recombination. |

**Example**

```python
import numpy as np
from mlx_class.background import Background
from mlx_class.khronon import KhrononBoltzmannSolver

bg = Background(khronon=True).solve()
k_arr = np.geomspace(5e-5, 0.35, 500).astype(np.float32)

solver = KhrononBoltzmannSolver(bg, k_arr, alpha_dbi=0.5)
result = solver.solve()

Theta_0, Phi, v_b = result.source_at_recombination()
cs2 = result.effective_cs2_at_rec()
```

---

## Analytic Transfer

```python
from mlx_class.perturbations import AnalyticTransfer
```

### `AnalyticTransfer(bg, khronon_correction=0.0)`

Analytic transfer function (no ODE integration). Computes the Sachs-Wolfe and Doppler sources using WKB approximations. Very fast (~0.3s total including C_l), but less accurate than ODE solvers for peak heights and phases.

**Parameters**

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `bg` | `Background` | -- | Solved background cosmology. |
| `khronon_correction` | `float` | `0.0` | Fractional correction to the source amplitude from Khronon. |

**Methods**

| Method | Returns | Description |
|--------|---------|-------------|
| `compute_source(k_arr)` | `ndarray` | Monopole source `Theta_0 + Phi` (Sachs-Wolfe) with Silk damping. |
| `compute_doppler(k_arr)` | `ndarray` | Doppler source `v_b` with Silk damping. |

**Example**

```python
import numpy as np
from mlx_class.background import Background
from mlx_class.perturbations import AnalyticTransfer

bg = Background().solve()
k_arr = np.geomspace(5e-5, 0.35, 500).astype(np.float32)

transfer = AnalyticTransfer(bg)
source_SW = transfer.compute_source(k_arr)
source_Dop = transfer.compute_doppler(k_arr)
```

---

## Spectra

### `compute_cl`

```python
from mlx_class.spectra import compute_cl
```

#### `compute_cl(source_SW, k_arr_Mpc, ell_values, D_A, source_Dop=None)`

Compute the CMB TT angular power spectrum from source functions and Bessel projection. Uses the instantaneous recombination approximation.

**Parameters**

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `source_SW` | `array (N_k,)` | -- | Sachs-Wolfe source `(Theta_0 + Phi) * silk`. |
| `k_arr_Mpc` | `array (N_k,)` | -- | Wavenumber array in Mpc^{-1}. |
| `ell_values` | `array (N_ell,)` | -- | Multipole values. |
| `D_A` | `float` | -- | Comoving angular diameter distance to LSS (Mpc). |
| `source_Dop` | `array (N_k,)` or `None` | `None` | Doppler source `v_b * silk`. If `None`, Doppler is omitted. |

**Returns**

| Name | Type | Description |
|------|------|-------------|
| `ell_values` | `array` | Multipole values (same as input). |
| `Cl` | `array` | Dimensionless `C_l`. |
| `Dl` | `array` | `l(l+1) C_l / (2 pi)` in muK^2. |

**Example**

```python
import numpy as np
from mlx_class.background import Background
from mlx_class.perturbations_implicit import ImplicitBoltzmannSolver
from mlx_class.spectra import compute_cl

bg = Background().solve()
k_arr = np.geomspace(5e-5, 0.35, 500).astype(np.float32)

solver = ImplicitBoltzmannSolver(bg, k_arr)
result = solver.solve()
Theta_0, Phi, v_b = result.source_at_recombination()

ell_values = np.arange(2, 2501, 4).astype(int)
ell, Cl, Dl = compute_cl(Theta_0 + Phi, k_arr, ell_values, bg.D_A)
```

---

### `compute_cl_los`

```python
from mlx_class.spectra_los import compute_cl_los
```

#### `compute_cl_los(bg, k_arr_Mpc, ell_values, use_ode=False, N_tau_vis=25, N_tau_early_isw=12, N_tau_late_isw=10)`

Compute C_l using proper line-of-sight (LOS) integration with visibility weighting, Doppler, early ISW, and late ISW contributions.

Two source modes:
- **Analytic** (`use_ode=False`): Hu-Sugiyama gravitational driving phase shift. Fast (~4s).
- **ODE** (`use_ode=True`): Run TCA Boltzmann solver through recombination. More accurate (~15s).

**Parameters**

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `bg` | `Background` | -- | Solved background cosmology. |
| `k_arr_Mpc` | `array (N_k,)` | -- | Wavenumber array in Mpc^{-1}. |
| `ell_values` | `array (N_ell,)` | -- | Multipole values. |
| `use_ode` | `bool` | `False` | Use ODE solver for source functions (more accurate, slower). |
| `N_tau_vis` | `int` | `25` | Number of quadrature points for the visibility integral. |
| `N_tau_early_isw` | `int` | `12` | Number of tau points for the early ISW integral. |
| `N_tau_late_isw` | `int` | `10` | Number of tau points for the late ISW integral. |

**Returns**

| Name | Type | Description |
|------|------|-------------|
| `ell_values` | `array` | Multipole values. |
| `Cl` | `array` | Dimensionless `C_l`. |
| `Dl` | `array` | `l(l+1) C_l / (2 pi)` in muK^2. |
| `info` | `dict` | Diagnostic info (decomposed D_l, visibility data, etc.). |

**Example**

```python
import numpy as np
from mlx_class.background import Background
from mlx_class.spectra_los import compute_cl_los

bg = Background(recombination='peebles').solve()
k_arr = np.geomspace(5e-4, 0.35, 500).astype(np.float64)
ell_values = np.arange(2, 2501, 4).astype(int)

ell, Cl, Dl, info = compute_cl_los(bg, k_arr, ell_values, use_ode=False)
print(f"SW-only peak: {np.max(info['Dl_sw_only']):.0f} uK^2")
```

---

### `find_peaks_in_spectrum`

```python
from mlx_class.spectra_los import find_peaks_in_spectrum
```

#### `find_peaks_in_spectrum(ell_values, Dl, n_peaks=7, smooth_sigma=15)`

Find acoustic peak positions in a D_l spectrum using Gaussian smoothing and peak finding.

**Parameters**

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `ell_values` | `array` | -- | Multipole values. |
| `Dl` | `array` | -- | `D_l` in muK^2. |
| `n_peaks` | `int` | `7` | Maximum number of peaks to return. |
| `smooth_sigma` | `int` | `15` | Gaussian smoothing width (in ell). |

**Returns**

| Name | Type | Description |
|------|------|-------------|
| `peaks` | `list` of `(l, D_l)` | Peak positions and heights. |
| `ell_fine` | `array` | Dense ell grid. |
| `Dl_fine` | `array` | Interpolated `D_l`. |
| `Dl_smooth` | `array` | Smoothed `D_l`. |

---

## Matter Power Spectrum

### `compute_matter_pk`

```python
from mlx_class.matter_pk import compute_matter_pk
```

#### `compute_matter_pk(result, solver_type='implicit', z_out=0.0, apply_silk=True, k_fine=None, N_k_fine=2000)`

Compute the matter power spectrum P(k) from a solver result. Extracts the transfer function at `tau_rec`, applies linear growth to the output redshift, and normalizes with the primordial spectrum.

**Parameters**

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `result` | solver result | -- | Output from any Boltzmann solver (Implicit, Neutrino, HiRes, or Explicit). |
| `solver_type` | `str` | `'implicit'` | `'implicit'`, `'neutrino'`, `'hires'`, or `'explicit'`. |
| `z_out` | `float` | `0.0` | Output redshift. |
| `apply_silk` | `bool` | `True` | Apply Silk damping to the transfer function. |
| `k_fine` | `array` or `None` | `None` | If given, interpolate P(k) onto this grid. |
| `N_k_fine` | `int` | `2000` | Number of fine-grid points (if `k_fine` is None and `N_k_fine > N_k`). |

**Returns**

| Name | Type | Description |
|------|------|-------------|
| `k_hMpc` | `array` | Wavenumber in h/Mpc. |
| `Pk_hMpc3` | `array` | P(k) in (Mpc/h)^3. |
| `k_Mpc` | `array` | Wavenumber in Mpc^{-1} (code units). |
| `Pk_Mpc3` | `array` | P(k) in Mpc^3 (code units). |

**Example**

```python
import numpy as np
from mlx_class.background import Background
from mlx_class.perturbations_implicit import ImplicitBoltzmannSolver
from mlx_class.matter_pk import compute_matter_pk

bg = Background().solve()
k_arr = np.geomspace(1e-4, 0.5, 1000).astype(np.float32)

solver = ImplicitBoltzmannSolver(bg, k_arr)
result = solver.solve()

k_h, Pk_h, k_Mpc, Pk_Mpc = compute_matter_pk(result, z_out=0.0)
# sigma_8 is printed automatically
```

---

## Polarization

### `compute_cl_ee`

```python
from mlx_class.polarization import compute_cl_ee
```

#### `compute_cl_ee(result, ell_values, solver_type='implicit')`

Compute the E-mode polarization power spectrum `C_l^EE`. Uses the photon quadrupole `Theta_2` at recombination (exact from HiRes solver, or estimated from the tight-coupling formula for other solvers).

**Parameters**

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `result` | solver result | -- | Output from any Boltzmann solver. |
| `ell_values` | `array` | -- | Multipole values (must be >= 2). |
| `solver_type` | `str` | `'implicit'` | `'implicit'`, `'neutrino'`, `'hires'`, or `'explicit'`. |

**Returns**

| Name | Type | Description |
|------|------|-------------|
| `ell_values` | `array` | Multipole values (>= 2). |
| `Cl_EE` | `array` | Dimensionless `C_l^EE`. |
| `Dl_EE` | `array` | `l(l+1) C_l^EE / (2 pi)` in muK^2. |

**Example**

```python
import numpy as np
from mlx_class.background import Background
from mlx_class.perturbations_implicit import ImplicitBoltzmannSolver
from mlx_class.polarization import compute_cl_ee

bg = Background().solve()
k_arr = np.geomspace(5e-4, 0.35, 500).astype(np.float32)

solver = ImplicitBoltzmannSolver(bg, k_arr)
result = solver.solve()

ell_values = np.arange(2, 2501, 4).astype(int)
ell, Cl_EE, Dl_EE = compute_cl_ee(result, ell_values, solver_type='implicit')
```

---

## Recombination

### `PeeblesRecombination`

```python
from mlx_class.recombination import PeeblesRecombination
```

#### `PeeblesRecombination(bg, z_start=1800.0, verbose=True)`

Peebles three-level atom recombination solver. Solves for `x_e(z)` using a two-phase approach: Saha equilibrium at high z, then the Peebles ODE. Provides ~0.1% accuracy vs RECFAST. Normally used internally by `Background(recombination='peebles')`.

**Parameters**

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `bg` | `Background` | -- | Solved background object. |
| `z_start` | `float` | `1800.0` | Starting redshift for integration. |
| `verbose` | `bool` | `True` | Print diagnostic information. |

**Methods**

| Method | Returns | Description |
|--------|---------|-------------|
| `solve()` | `ndarray` or `None` | Free electron fraction `x_e(a)` on the background grid (including He correction). Returns `None` if the ODE fails. |

---

## Constants

Key cosmological constants are defined in `mlx_class.background`:

```python
from mlx_class.background import (
    h,          # 0.6736 (Planck 2018)
    H0_Mpc,     # H0 in 1/Mpc
    Omega_b,    # baryon density parameter
    Omega_c,    # CDM density parameter
    Omega_r,    # radiation density parameter (photons + neutrinos)
    Omega_m,    # total matter
    Omega_L,    # dark energy
    T_CMB,      # 2.7255 K
    A_s,        # 2.1e-9 (primordial amplitude)
    n_s,        # 0.9649 (spectral index)
    k_pivot,    # 0.05 Mpc^{-1}
    z_rec,      # 1089.92
    a_rec,      # 1 / (1 + z_rec)
    a_eq,       # matter-radiation equality
    k_eq,       # equality wavenumber
    Y_He,       # 0.2454 (helium fraction)
)
```

---

## CLI Usage

The package can be run directly from the command line:

```bash
# Analytic transfer (fastest, ~0.3s)
python -m mlx_class.main

# IMEX solver (integrates to tau_rec, ~2.5s)
python -m mlx_class.main --implicit

# IMEX with matter power spectrum and polarization
python -m mlx_class.main --implicit --pk --ee

# Fast mode (N_k=200, ~0.6s)
python -m mlx_class.main --fast

# Khronon dark matter comparison
python -m mlx_class.khronon

# Combined solver (Peebles + neutrinos + ODE + LOS)
python -m mlx_class.solver_combined

# Full CLASS comparison (requires classy)
python -m mlx_class.class_comparison
```

**CLI flags for `main.py`:**

| Flag | Description |
|------|-------------|
| `--ode` | Use explicit RK4 solver |
| `--implicit` | Use IMEX solver (recommended) |
| `--khronon` | Use Khronon dark matter model |
| `--fast` | Fast mode (N_k=200, sparse grids) |
| `--pk` | Output matter power spectrum |
| `--ee` | Output E-mode polarization |
| `--compare` | Compare IMEX vs explicit solver |
| `--N_k N` | Number of k-modes (default: 500) |
| `--ell_max L` | Maximum multipole (default: 2500) |
| `--z_pk Z` | Output redshift for P(k) (default: 0) |
| `--no_plot` | Suppress plot generation |
| `--no_isw` | Disable early ISW correction |
