"""
solver_api.py -- Unified CMB solver API for mlx_class.

Usage:
    from mlx_class import CMBSolver

    # Quick (sync gauge, scipy, ~420s but accurate)
    solver = CMBSolver()
    result = solver.run()
    print(result.Dl_TT)  # numpy array
    print(result.peaks)   # peak positions and amplitudes

    # With custom parameters
    solver = CMBSolver(h=0.6736, omega_b=0.02237, omega_cdm=0.12)
    result = solver.run()

    # Fast mode (conformal Newtonian, ~2s but less accurate)
    solver = CMBSolver(backend='fast')
    result = solver.run()

Author: Sheng-Kai Huang, 2026
"""
from dataclasses import dataclass, field
from typing import Optional, Dict, List
import numpy as np


@dataclass
class CMBResult:
    """Unified result from any CMB solver backend."""

    # Multipoles
    ell: np.ndarray

    # Power spectra in D_l = l(l+1)C_l/(2pi) [muK^2]
    Dl_TT: np.ndarray
    Dl_EE: Optional[np.ndarray] = None
    Dl_TE: Optional[np.ndarray] = None
    Dl_BB: Optional[np.ndarray] = None

    # Raw C_l (dimensionless)
    Cl_TT: Optional[np.ndarray] = None
    Cl_EE: Optional[np.ndarray] = None
    Cl_TE: Optional[np.ndarray] = None

    # Matter power spectrum
    k_pk: Optional[np.ndarray] = None
    Pk: Optional[np.ndarray] = None

    # Peaks
    peaks: Optional[Dict] = None

    # Timing
    timing: Optional[Dict] = None

    # Backend info
    backend: str = 'sync'
    rms_vs_class: Optional[float] = None


class CMBSolver:
    """
    Unified CMB power spectrum solver.

    Parameters
    ----------
    h : float
        Dimensionless Hubble parameter (default: 0.6736).
    omega_b : float
        Physical baryon density (default: 0.02237).
    omega_cdm : float
        Physical CDM density (default: 0.12).
    backend : str
        'sync' (synchronous gauge, accurate, ~420s) or
        'fast' (conformal Newtonian IMEX, ~2s, ~60% RMS).
    N_k : int or None
        Number of k-modes (default: 180 for sync, 500 for fast).
    recombination : str
        'recfast' (default for sync) or 'peebles' (default for fast).
    verbose : bool
        Print progress (default: True).
    """

    def __init__(self, h=0.6736, omega_b=0.02237, omega_cdm=0.12,
                 backend='sync', N_k=None, recombination=None,
                 verbose=True):
        self.h = h
        self.omega_b = omega_b
        self.omega_cdm = omega_cdm
        self.backend = backend
        self.N_k = N_k
        self.recombination = recombination
        self.verbose = verbose

    def run(self) -> CMBResult:
        """Run the solver and return a CMBResult."""
        if self.backend == 'sync':
            return self._run_sync()
        elif self.backend == 'fast':
            return self._run_fast()
        else:
            raise ValueError(f"Unknown backend: {self.backend!r}. "
                             f"Choose 'sync' or 'fast'.")

    def _run_sync(self) -> CMBResult:
        from .solver_sync import run_sync_solver
        N_k = self.N_k or 180
        raw = run_sync_solver(N_k=N_k, verbose=self.verbose)

        # Find peaks
        peaks = self._find_peaks(raw['ell'], raw['Dl_TT'])

        return CMBResult(
            ell=raw['ell'],
            Dl_TT=raw['Dl_TT'],
            Dl_EE=raw.get('Dl_EE'),
            Dl_TE=raw.get('Dl_TE'),
            Cl_TT=raw.get('Cl_TT'),
            Cl_EE=raw.get('Cl_EE'),
            Cl_TE=raw.get('Cl_TE'),
            peaks=peaks,
            timing=raw.get('timing'),
            backend='sync',
        )

    def _run_fast(self) -> CMBResult:
        from .solver_production import ProductionSolver
        N_k = self.N_k or 500
        recombination = self.recombination or 'peebles'
        solver = ProductionSolver(
            h=self.h, omega_b=self.omega_b, omega_cdm=self.omega_cdm,
            N_k=N_k, recombination=recombination, verbose=self.verbose)
        raw = solver.run()

        return CMBResult(
            ell=raw.ell,
            Dl_TT=raw.Dl_TT,
            Dl_EE=raw.Dl_EE,
            Dl_TE=raw.Dl_TE,
            Cl_TT=raw.Cl_TT,
            peaks=raw.peaks,
            timing=raw.timing,
            backend='fast',
        )

    @staticmethod
    def _find_peaks(ell, Dl, min_ell=100, n_peaks=7):
        """Find acoustic peaks in D_l spectrum."""
        from scipy.signal import find_peaks as _find_peaks
        from scipy.ndimage import gaussian_filter1d

        Dl_smooth = gaussian_filter1d(Dl, sigma=8)
        mask = ell > min_ell
        if not np.any(mask):
            return {'ell': [], 'Dl': []}

        offset = np.argmax(mask)
        pks, _ = _find_peaks(Dl_smooth[mask], distance=60, prominence=10)
        pks = pks + offset

        return {
            'ell': [int(ell[p]) for p in pks[:n_peaks]],
            'Dl': [float(Dl[p]) for p in pks[:n_peaks]],
        }
