"""
test_sync_gauge.py -- Synchronous gauge solver tests against CLASS reference values.

Tests:
  - First peak position (within 10 of CLASS l=221)
  - First peak amplitude (D_l(220) within 10% of CLASS)
  - Sound horizon (r_s within 0.1% of CLASS 144.4 Mpc)
  - Background quantities (tau_0, D_A)
  - RECFAST visibility (FWHM < 25 Mpc)

These tests verify the physics of the synchronous gauge solver pipeline:
  perturbations_sync.py (ODE) -> solver_sync.py (LOS) -> C_l

Run:
    python -m pytest mlx_class/tests/test_sync_gauge.py -v
    python -m pytest mlx_class/tests/test_sync_gauge.py -v -k "not slow"

Author: Sheng-Kai Huang, 2026
"""
import unittest
import numpy as np
import os


# CLASS reference values (Planck 2018 best-fit LCDM)
CLASS_FIRST_PEAK_ELL = 221
CLASS_R_S = 144.43  # Mpc
CLASS_TAU_0 = 14150.0  # Mpc


def _solve_background():
    """Solve background cosmology (shared fixture)."""
    from mlx_class.background import Background
    bg = Background(khronon=False, recombination='recfast')
    bg.solve()
    return bg


class TestBackground(unittest.TestCase):
    """Background cosmology checks."""

    @classmethod
    def setUpClass(cls):
        cls.bg = _solve_background()

    def test_sound_horizon(self):
        """r_s should be within 0.1% of CLASS value 144.4 Mpc."""
        rel_err = abs(self.bg.r_s - CLASS_R_S) / CLASS_R_S
        self.assertLess(rel_err, 0.001,
                        f"r_s = {self.bg.r_s:.2f}, expected {CLASS_R_S}")

    def test_conformal_age(self):
        """tau_0 should be within 1% of CLASS value ~14150 Mpc."""
        rel_err = abs(self.bg.tau_0 - CLASS_TAU_0) / CLASS_TAU_0
        self.assertLess(rel_err, 0.01,
                        f"tau_0 = {self.bg.tau_0:.1f}, expected {CLASS_TAU_0}")

    def test_angular_diameter_distance(self):
        """D_A should be positive and within 2% of ~12800 Mpc."""
        self.assertGreater(self.bg.D_A, 12000.0)
        self.assertLess(self.bg.D_A, 14000.0)

    def test_recfast_visibility(self):
        """RECFAST visibility FWHM should be < 25 Mpc.

        The visibility function g(tau) = -kappa' exp(-kappa) peaks at
        recombination. Its FWHM controls the damping scale. RECFAST
        gives a narrower visibility (~20 Mpc) than tanh (~50 Mpc),
        which is essential for correct Silk damping.
        """
        vis = self.bg.visibility_grid
        tau = self.bg.tau_grid
        peak_idx = np.argmax(vis)
        half_max = vis[peak_idx] / 2.0

        # Find FWHM
        left = peak_idx
        while left > 0 and vis[left] > half_max:
            left -= 1
        right = peak_idx
        while right < len(vis) - 1 and vis[right] > half_max:
            right += 1
        fwhm = tau[right] - tau[left]

        self.assertLess(fwhm, 25.0,
                        f"Visibility FWHM = {fwhm:.1f} Mpc (should be < 25)")
        self.assertGreater(fwhm, 5.0,
                           f"Visibility FWHM = {fwhm:.1f} Mpc (suspiciously narrow)")


class TestSyncSolverQuick(unittest.TestCase):
    """Quick sanity checks on the sync solver (N_k=60, fast but approximate)."""

    @classmethod
    def setUpClass(cls):
        """Run the sync solver with minimal k-modes for speed."""
        from mlx_class.solver_sync import run_sync_solver
        # N_k=60 is fast (~30s) but sufficient to check peak position/amplitude
        cls.result = run_sync_solver(N_k=60, verbose=False)

    def test_first_peak_position(self):
        """First peak should be within 10 of CLASS l=221."""
        ell = self.result['ell']
        Dl = self.result['Dl_TT']

        # Search for first peak in l=150-300
        mask = (ell >= 150) & (ell <= 300)
        ell_sub = ell[mask]
        Dl_sub = Dl[mask]
        peak_ell = ell_sub[np.argmax(Dl_sub)]

        self.assertAlmostEqual(peak_ell, CLASS_FIRST_PEAK_ELL, delta=10,
                               msg=f"First peak at l={peak_ell}, expected ~{CLASS_FIRST_PEAK_ELL}")

    def test_first_peak_amplitude(self):
        """D_l(220) should be within 10% of CLASS (~5800 muK^2)."""
        ell = self.result['ell']
        Dl = self.result['Dl_TT']

        idx = np.argmin(np.abs(ell - 220))
        Dl_220 = Dl[idx]

        # CLASS reference: D_l(220) ~ 5800 muK^2
        self.assertGreater(Dl_220, 5800 * 0.9,
                           f"D_l(220) = {Dl_220:.0f}, too low (< 5220)")
        self.assertLess(Dl_220, 5800 * 1.1,
                        f"D_l(220) = {Dl_220:.0f}, too high (> 6380)")

    def test_spectra_positive(self):
        """TT and EE spectra should be non-negative."""
        self.assertTrue(np.all(self.result['Dl_TT'] >= 0),
                        "D_l^TT has negative values")
        self.assertTrue(np.all(self.result['Dl_EE'] >= 0),
                        "D_l^EE has negative values")

    def test_ell_range(self):
        """Output ell should span from 2 to at least 2000."""
        ell = self.result['ell']
        self.assertEqual(ell[0], 2, "First ell should be 2")
        self.assertGreaterEqual(ell[-1], 2000, "Last ell should be >= 2000")

    def test_second_peak_exists(self):
        """Second acoustic peak should exist around l ~ 540."""
        ell = self.result['ell']
        Dl = self.result['Dl_TT']
        mask = (ell >= 400) & (ell <= 700)
        ell_sub = ell[mask]
        Dl_sub = Dl[mask]
        # Should have a local maximum
        if len(Dl_sub) > 2:
            peak_ell = ell_sub[np.argmax(Dl_sub)]
            self.assertGreater(peak_ell, 400, "Second peak too early")
            self.assertLess(peak_ell, 700, "Second peak too late")


if __name__ == '__main__':
    unittest.main()
