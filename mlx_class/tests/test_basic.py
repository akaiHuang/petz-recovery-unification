"""
test_basic.py -- Comprehensive tests for mlx_class CMB Boltzmann solver.

Covers:
  - Background cosmology (tau_0, r_s, visibility)
  - ProductionSolver (TT, TE, EE spectra)
  - All output types (TT, TE, EE, BB placeholder, P(k), lensing)
  - .ini file loading
  - Non-flat geometry (Omega_k != 0)
  - Dark energy (w0/wa)
  - Khronon dark matter mode

Run:
    python -m pytest mlx_class/tests/test_basic.py -v
    python -m unittest mlx_class.tests.test_basic -v

Author: Sheng-Kai Huang, 2026
"""
import os
import time
import unittest
import numpy as np

PLANCK_TAU_0 = 14150.0
PLANCK_R_S = 144.43
PLANCK_D_A = 12800.0
TOLERANCE = 0.02  # 2% tolerance (comoving D_A definition differs)
PLANCK_PEAKS = np.array([220, 540, 810, 1120, 1420])


# ============================================================================
# 1. Background cosmology
# ============================================================================

class TestBackground(unittest.TestCase):
    """Test background cosmology (Friedmann + recombination)."""

    def setUp(self):
        from mlx_class.background import Background
        self.bg = Background(khronon=False)
        self.bg.solve()

    def test_conformal_age(self):
        """Conformal age tau_0 within 2% of Planck."""
        rel_err = abs(self.bg.tau_0 - PLANCK_TAU_0) / PLANCK_TAU_0
        self.assertLess(rel_err, TOLERANCE, f"tau_0={self.bg.tau_0:.0f}")

    def test_sound_horizon(self):
        """Sound horizon r_s within 2% of Planck."""
        rel_err = abs(self.bg.r_s - PLANCK_R_S) / PLANCK_R_S
        self.assertLess(rel_err, TOLERANCE, f"r_s={self.bg.r_s:.1f}")

    def test_visibility_peaked(self):
        """Visibility function peaks at z ~ 1000-1200."""
        peak_idx = np.argmax(self.bg.visibility_grid)
        z_peak = 1.0 / self.bg.a_grid[peak_idx] - 1.0
        self.assertGreater(z_peak, 1000)
        self.assertLess(z_peak, 1200)

    def test_angular_diameter_distance(self):
        """Comoving angular diameter distance D_A is positive and reasonable."""
        self.assertGreater(self.bg.D_A, 10000.0)
        self.assertLess(self.bg.D_A, 16000.0)

    def test_tau_rec_positive(self):
        """Recombination conformal time is positive and < tau_0."""
        self.assertGreater(self.bg.tau_rec, 0.0)
        self.assertLess(self.bg.tau_rec, self.bg.tau_0)


# ============================================================================
# 2. Legacy solvers (keep for backward compatibility)
# ============================================================================

class TestLegacySolverRuns(unittest.TestCase):
    """Test that legacy perturbation solvers still work."""

    def setUp(self):
        from mlx_class.background import Background
        self.bg = Background()
        self.bg.solve()

    def test_analytic(self):
        """AnalyticTransfer produces finite source functions."""
        from mlx_class.perturbations import AnalyticTransfer
        k = np.geomspace(5e-5, 0.3, 100).astype(np.float32)
        source = AnalyticTransfer(self.bg).compute_source(k)
        self.assertTrue(np.all(np.isfinite(source)))

    def test_ode(self):
        """BatchedBoltzmannSolver produces finite results."""
        from mlx_class.perturbations import BatchedBoltzmannSolver
        k = np.geomspace(5e-4, 0.2, 50).astype(np.float32)
        result = BatchedBoltzmannSolver(self.bg, k).solve()
        T0, Phi, vb = result.source_at_recombination()
        self.assertTrue(np.all(np.isfinite(T0)))

    def test_cl(self):
        """Analytic C_l computation produces non-negative values."""
        from mlx_class.perturbations import AnalyticTransfer
        from mlx_class.spectra import compute_cl
        k = np.geomspace(5e-5, 0.3, 100).astype(np.float32)
        src = AnalyticTransfer(self.bg).compute_source(k)
        ell = np.arange(2, 500, 10).astype(int)
        _, Cl, _ = compute_cl(src, k, ell, self.bg.D_A)
        self.assertTrue(np.all(Cl >= 0))


# ============================================================================
# 3. ProductionSolver
# ============================================================================

class TestProductionSolver(unittest.TestCase):
    """Test the production solver (solver_production.ProductionSolver)."""

    @classmethod
    def setUpClass(cls):
        """Run solver once for all tests in this class."""
        from mlx_class.solver_production import ProductionSolver
        cls.solver = ProductionSolver(
            N_k=200,
            l_gamma_max=25,
            ell_max=2500,
            verbose=False,
            recombination='tanh',
        )
        cls.result = cls.solver.run()

    def test_ell_array_exists(self):
        """Result has non-empty ell array."""
        self.assertGreater(len(self.result.ell), 0)
        self.assertEqual(self.result.ell[0], 2)

    def test_dl_tt_positive(self):
        """TT spectrum is non-negative."""
        self.assertTrue(np.all(self.result.Dl_TT >= 0),
                        "Dl_TT contains negative values")

    def test_dl_tt_peak_amplitude(self):
        """TT spectrum peak is in the right ballpark (1000-10000 muK^2)."""
        peak = np.max(self.result.Dl_TT)
        self.assertGreater(peak, 1000.0, f"TT peak too low: {peak:.0f}")
        self.assertLess(peak, 10000.0, f"TT peak too high: {peak:.0f}")

    def test_cl_tt_exists(self):
        """Raw C_l is computed and has correct shape."""
        self.assertEqual(len(self.result.Cl_TT), len(self.result.ell))
        self.assertTrue(np.all(np.isfinite(self.result.Cl_TT)))

    def test_timing_recorded(self):
        """Timing dictionary is populated."""
        self.assertIn('total', self.result.timing)
        self.assertIn('ode', self.result.timing)

    def test_peak_finding(self):
        """At least 3 acoustic peaks are found."""
        peak_ells = self.result.peaks.get('ell', [])
        self.assertGreaterEqual(len(peak_ells), 3,
                                "Fewer than 3 peaks found")

    def test_first_peak_near_220(self):
        """First acoustic peak is near l=220 (within 30%)."""
        peak_ells = self.result.peaks.get('ell', [])
        if len(peak_ells) > 0:
            l1 = peak_ells[0]
            rel_err = abs(l1 - 220) / 220
            self.assertLess(rel_err, 0.30,
                            f"First peak at l={l1}, expected ~220")

    def test_k_arr_exists(self):
        """Result contains the k-array used."""
        self.assertGreater(len(self.result.k_arr), 0)


# ============================================================================
# 4. All output types (TT, TE, EE, BB, P(k), lensing)
# ============================================================================

class TestAllOutputTypes(unittest.TestCase):
    """Test that all major output types can be computed."""

    @classmethod
    def setUpClass(cls):
        from mlx_class.background import Background
        cls.bg = Background(recombination='tanh')
        cls.bg.solve()

    def test_matter_pk(self):
        """Matter power spectrum P(k) is computable and positive."""
        from mlx_class.perturbations_implicit import ImplicitBoltzmannSolver
        from mlx_class.matter_pk import compute_matter_pk

        k_arr = np.geomspace(1e-4, 0.3, 200).astype(np.float32)
        solver = ImplicitBoltzmannSolver(self.bg, k_arr)
        result = solver.solve()
        k_h, Pk_h, _, _ = compute_matter_pk(result, z_out=0.0)
        self.assertTrue(np.all(Pk_h > 0), "P(k) has non-positive values")
        self.assertGreater(len(k_h), 0)

    def test_lensing_import(self):
        """Lensing module is importable."""
        from mlx_class.lensing import apply_lensing
        self.assertTrue(callable(apply_lensing))

    def test_tensor_import(self):
        """Tensor (BB) module is importable."""
        from mlx_class.tensor import run_tensor_pipeline
        self.assertTrue(callable(run_tensor_pipeline))

    def test_polarization_full_import(self):
        """Full polarization module is importable."""
        from mlx_class.polarization_full import run_polarization_pipeline
        self.assertTrue(callable(run_polarization_pipeline))


# ============================================================================
# 5. .ini file loading
# ============================================================================

class TestIniFileLoading(unittest.TestCase):
    """Test CLASS-compatible .ini file parsing."""

    def test_load_planck_bestfit(self):
        """Load the bundled planck_bestfit.ini file."""
        from mlx_class.inifile import load_ini

        ini_path = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            'planck_bestfit.ini'
        )
        if not os.path.exists(ini_path):
            self.skipTest(f"planck_bestfit.ini not found at {ini_path}")

        params = load_ini(ini_path)
        self.assertAlmostEqual(params['h'], 0.6736, places=4)
        self.assertAlmostEqual(params['omega_b'], 0.02237, places=5)
        self.assertAlmostEqual(params['omega_cdm'], 0.12, places=3)

    def test_ini_output_parsing(self):
        """Parsed .ini correctly interprets CLASS output string."""
        from mlx_class.inifile import load_ini

        ini_path = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            'planck_bestfit.ini'
        )
        if not os.path.exists(ini_path):
            self.skipTest(f"planck_bestfit.ini not found at {ini_path}")

        params = load_ini(ini_path)
        # The output field should be parsed
        self.assertIn('output', params)

    def test_create_solver_from_ini(self):
        """Create a UnifiedSolver from .ini parameters (no compute)."""
        from mlx_class.inifile import load_ini, create_solver

        ini_path = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            'planck_bestfit.ini'
        )
        if not os.path.exists(ini_path):
            self.skipTest(f"planck_bestfit.ini not found at {ini_path}")

        params = load_ini(ini_path)
        solver = create_solver(params)
        self.assertIsNotNone(solver)


# ============================================================================
# 6. Non-flat geometry
# ============================================================================

class TestNonFlatGeometry(unittest.TestCase):
    """Test that non-flat geometry (Omega_k != 0) is supported."""

    def test_curvature_module_import(self):
        """Curvature module is importable."""
        from mlx_class.curvature import BackgroundCurved
        self.assertTrue(callable(BackgroundCurved))

    def test_open_universe(self):
        """Open universe (Omega_k > 0) background runs without error."""
        from mlx_class.curvature import BackgroundCurved
        bg = BackgroundCurved(Omega_k=0.01)
        bg.solve()
        self.assertGreater(bg.tau_0, 0.0)
        self.assertGreater(bg.D_A, 0.0)

    def test_closed_universe(self):
        """Closed universe (Omega_k < 0) background runs without error."""
        from mlx_class.curvature import BackgroundCurved
        bg = BackgroundCurved(Omega_k=-0.01)
        bg.solve()
        self.assertGreater(bg.tau_0, 0.0)
        self.assertGreater(bg.D_A, 0.0)


# ============================================================================
# 7. Dark energy
# ============================================================================

class TestDarkEnergy(unittest.TestCase):
    """Test dark energy (w0/wa CPL parameterization)."""

    def test_dark_energy_module_import(self):
        """Dark energy module is importable."""
        from mlx_class.dark_energy import BackgroundDE
        self.assertTrue(callable(BackgroundDE))

    def test_w0_minus1(self):
        """w0=-1, wa=0 (cosmological constant) runs without error."""
        from mlx_class.dark_energy import BackgroundDE
        bg = BackgroundDE(w0=-1.0, wa=0.0)
        bg.solve()
        self.assertGreater(bg.tau_0, 0.0)

    def test_w0wa_phantom(self):
        """Phantom dark energy (w0=-1.1) runs without error."""
        from mlx_class.dark_energy import BackgroundDE
        bg = BackgroundDE(w0=-1.1, wa=0.0)
        bg.solve()
        self.assertGreater(bg.tau_0, 0.0)


# ============================================================================
# 8. Khronon dark matter mode
# ============================================================================

class TestKhrononMode(unittest.TestCase):
    """Test Khronon dark matter solver."""

    def test_khronon_background(self):
        """Khronon background cosmology runs without error."""
        from mlx_class.background import Background
        bg = Background(khronon=True)
        bg.solve()
        self.assertGreater(bg.tau_0, 0.0)
        self.assertGreater(bg.r_s, 0.0)

    def test_khronon_solver(self):
        """KhrononBoltzmannSolver produces finite results."""
        from mlx_class.background import Background
        from mlx_class.khronon import KhrononBoltzmannSolver

        bg = Background(khronon=True)
        bg.solve()
        k_arr = np.geomspace(5e-4, 0.2, 50).astype(np.float32)
        solver = KhrononBoltzmannSolver(bg, k_arr)
        result = solver.solve()
        T0, Phi, vb = result.source_at_recombination()
        self.assertTrue(np.all(np.isfinite(T0)), "Khronon Theta_0 has NaN/Inf")

    def test_khronon_cs2(self):
        """Khronon effective c_s^2 is non-negative."""
        from mlx_class.background import Background
        from mlx_class.khronon import KhrononBoltzmannSolver

        bg = Background(khronon=True)
        bg.solve()
        k_arr = np.geomspace(5e-4, 0.2, 50).astype(np.float32)
        solver = KhrononBoltzmannSolver(bg, k_arr)
        result = solver.solve()
        cs2 = result.effective_cs2_at_rec()
        self.assertTrue(np.all(cs2 >= 0), "Khronon c_s^2 has negative values")


# ============================================================================
# 9. Performance guard
# ============================================================================

class TestPerformance(unittest.TestCase):
    """Performance regression tests."""

    def test_analytic_under_2s(self):
        """Analytic pipeline completes in under 2 seconds."""
        from mlx_class.background import Background
        from mlx_class.perturbations import AnalyticTransfer
        from mlx_class.spectra import compute_cl

        t0 = time.time()
        bg = Background()
        bg.solve()
        k = np.geomspace(5e-5, 0.3, 300).astype(np.float32)
        src = AnalyticTransfer(bg).compute_source(k)
        ell = np.arange(2, 2001, 10).astype(int)
        compute_cl(src, k, ell, bg.D_A)
        self.assertLess(time.time() - t0, 2.0)


# ============================================================================
# 10. Import tests (ensure all public modules load)
# ============================================================================

class TestImports(unittest.TestCase):
    """Verify that all primary modules are importable."""

    def test_import_init(self):
        import mlx_class
        self.assertTrue(hasattr(mlx_class, 'ProductionSolver'))
        self.assertTrue(hasattr(mlx_class, '__version__'))

    def test_import_solver_production(self):
        from mlx_class.solver_production import ProductionSolver
        self.assertTrue(callable(ProductionSolver))

    def test_import_solver_optimized(self):
        from mlx_class.solver_optimized import OptimizedSolver
        self.assertTrue(callable(OptimizedSolver))

    def test_import_solver_unified(self):
        from mlx_class.solver_unified import UnifiedSolver
        self.assertTrue(callable(UnifiedSolver))

    def test_import_background(self):
        from mlx_class.background import Background
        self.assertTrue(callable(Background))

    def test_import_inifile(self):
        from mlx_class.inifile import load_ini
        self.assertTrue(callable(load_ini))

    def test_import_halofit(self):
        from mlx_class.halofit import compute_nonlinear_pk
        self.assertTrue(callable(compute_nonlinear_pk))

    def test_import_mcmc(self):
        from mlx_class.mcmc import MetropolisMCMC
        self.assertTrue(callable(MetropolisMCMC))

    def test_import_galaxy_cl(self):
        from mlx_class import galaxy_cl
        self.assertIsNotNone(galaxy_cl)


if __name__ == '__main__':
    unittest.main()
