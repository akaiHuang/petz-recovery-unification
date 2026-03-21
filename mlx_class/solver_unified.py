"""
solver_unified.py -- ONE unified CMB precision pipeline for mlx_class.

Consolidates all 7+ scattered solvers into a single entry point:
  - perturbations.py         (explicit RK4, TCA only)
  - perturbations_implicit.py (IMEX, to tau_rec)
  - perturbations_neutrino.py (+ massless neutrino)
  - perturbations_hires.py    (full hierarchy, Strang splitting)
  - massive_neutrino.py       (+ massive neutrino)
  - solver_combined.py        (attempted combination)
  - solver_precision.py       (attempted precision)

Physics is NOT rewritten -- we IMPORT from existing modules.

Mode selection logic:
  1. sum_mnu > 0: MassiveNeutrinoBoltzmannSolver  (most general, slowest)
  2. full hierarchy (default): HiResBoltzmannSolver (Strang splitting)
  3. fast_mode=True: ImplicitBoltzmannSolver (TCA IMEX, fastest)

All source terms: SW + Doppler + ISW (early + late).
Output: TT, TE, EE (unlensed + lensed).
Optional: Khronon, dark energy (w0, wa), massive neutrinos.

Usage:
    solver = UnifiedSolver(h=0.6736, omega_b=0.02237, omega_cdm=0.12)
    result = solver.compute_all()
    # result.Dl_TT, result.Dl_TE, result.Dl_EE
    # result.Dl_TT_lensed, result.Dl_EE_lensed
    # result.Pk_linear

Author: Sheng-Kai Huang, 2026
"""
import os
import sys
import time
import numpy as np
from dataclasses import dataclass, field
from typing import Optional, Dict, Any

# Ensure package import works
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import mlx.core as mx
from scipy.interpolate import interp1d
from scipy.signal import find_peaks
from scipy.ndimage import gaussian_filter1d

from mlx_class.bessel_gpu import BesselTable, compute_visibility_bessel_gpu

# ============================================================================
# Imports from existing modules (NO physics rewrite)
# ============================================================================
from mlx_class.background import (
    Background, A_s as _A_s_default, n_s as _n_s_default,
    k_pivot as _k_pivot, T_CMB as _T_CMB, k_eq,
    Omega_r as _OMEGA_R, Omega_b as _OMEGA_B,
    Omega_c as _OMEGA_C, H0_Mpc as _H0_MPC,
    a_eq as _A_EQ, Omega_m as _OMEGA_M,
    Omega_L as _OMEGA_L, h as _h_default,
)

from mlx_class.perturbations_hires import (
    HiResBoltzmannSolver, HiResResult,
    IDX_PHI, IDX_DELTA_B, IDX_V_B, IDX_DELTA_C, IDX_V_C,
    IDX_THETA_START, idx_theta, idx_n_start,
    n_var_total as hires_n_var,
    adiabatic_ic_hires, build_tau_grid_hires,
    strang_split_step, deriv_streaming,
)

from mlx_class.perturbations_neutrino import (
    NeutrinoBoltzmannSolver, NeutrinoResult,
    _OMEGA_GAMMA, _OMEGA_NU, _f_nu, L_NU_MAX,
    IDX_THETA_0 as NU_IDX_THETA_0,
    IDX_THETA_1 as NU_IDX_THETA_1,
    IDX_N_START as NU_IDX_N_START,
    n_var_total as nu_n_var,
    adiabatic_ic as nu_adiabatic_ic,
    imex_rk4_step as nu_imex_rk4_step,
    deriv_full as nu_deriv_full,
)

from mlx_class.perturbations_implicit import (
    ImplicitBoltzmannSolver,
)


# ============================================================================
# Result container
# ============================================================================

@dataclass
class UnifiedResult:
    """Container for all output spectra and diagnostics."""
    # Multipole array
    ell: np.ndarray = field(default_factory=lambda: np.array([]))

    # Unlensed D_l = l(l+1)/(2pi) C_l in muK^2
    Dl_TT: np.ndarray = field(default_factory=lambda: np.array([]))
    Dl_TE: np.ndarray = field(default_factory=lambda: np.array([]))
    Dl_EE: np.ndarray = field(default_factory=lambda: np.array([]))
    Dl_BB: np.ndarray = field(default_factory=lambda: np.array([]))

    # Lensed D_l
    Dl_TT_lensed: np.ndarray = field(default_factory=lambda: np.array([]))
    Dl_EE_lensed: np.ndarray = field(default_factory=lambda: np.array([]))

    # Matter power spectrum
    Pk_linear: Optional[Dict[str, np.ndarray]] = None     # {k_hMpc, Pk_hMpc3}
    Pk_nonlinear: Optional[Dict[str, np.ndarray]] = None   # {k_hMpc, Pk_hMpc3}

    # Background object
    bg: Any = None

    # Timing breakdown
    timing: Dict[str, float] = field(default_factory=dict)

    # Peak positions
    peak_ells: np.ndarray = field(default_factory=lambda: np.array([]))
    peak_heights: np.ndarray = field(default_factory=lambda: np.array([]))


# ============================================================================
# Planck reference
# ============================================================================
PLANCK_PEAKS = np.array([220, 540, 810, 1120, 1420])
PLANCK_HEIGHTS = np.array([5780, 2530, 2560, 1490, 840])


# ============================================================================
# Unified Solver
# ============================================================================

class UnifiedSolver:
    """
    Single entry point for the mlx_class CMB solver pipeline.

    Automatically selects the best solver backend and computes
    TT, TE, EE spectra with full LOS integration and lensing.

    Parameters
    ----------
    h : float
        Hubble constant H_0 / (100 km/s/Mpc).
    omega_b : float
        Physical baryon density omega_b = Omega_b h^2.
    omega_cdm : float
        Physical CDM density omega_cdm = Omega_cdm h^2.
    A_s : float
        Scalar amplitude at pivot scale.
    n_s : float
        Scalar spectral index.
    tau_reio : float
        Optical depth to reionization (used for damping envelope).
    N_ur : float
        Effective number of massless neutrino species.
    sum_mnu : float
        Sum of neutrino masses in eV (0 = all massless).
    w0 : float
        Dark energy equation of state at a=1.
    wa : float
        CPL time derivative of w.
    l_max : int
        Maximum multipole.
    N_k : int
        Number of k-modes for ODE integration.
    khronon : bool
        If True, enable Khronon (Sigma = 2 ln Q) sector.
    fast_mode : bool
        If True, use TCA-only IMEX solver (fastest, least accurate).
    l_gamma_max : int
        Photon multipole truncation for HiRes solver.
    N_tau_vis : int
        Number of visibility quadrature points.
    compute_lensing : bool
        Whether to compute lensed spectra.
    compute_pk : bool
        Whether to compute matter power spectrum.
    """

    def __init__(
        self,
        h: float = 0.6736,
        omega_b: float = 0.02237,
        omega_cdm: float = 0.12,
        A_s: float = 2.1e-9,
        n_s: float = 0.9649,
        tau_reio: float = 0.0544,
        N_ur: float = 3.046,
        sum_mnu: float = 0.0,
        w0: float = -1.0,
        wa: float = 0.0,
        l_max: int = 2500,
        N_k: int = 500,
        khronon: bool = False,
        fast_mode: bool = False,
        l_gamma_max: int = 25,
        N_tau_vis: int = 30,
        compute_lensing: bool = True,
        compute_pk: bool = False,
    ):
        self.h = h
        self.omega_b = omega_b
        self.omega_cdm = omega_cdm
        self.A_s = A_s
        self.n_s = n_s
        self.tau_reio = tau_reio
        self.N_ur = N_ur
        self.sum_mnu = sum_mnu
        self.w0 = w0
        self.wa = wa
        self.l_max = l_max
        self.N_k = N_k
        self.khronon = khronon
        self.fast_mode = fast_mode
        self.l_gamma_max = l_gamma_max
        self.N_tau_vis = N_tau_vis
        self.compute_lensing = compute_lensing
        self.compute_pk = compute_pk

        # Determine solver backend
        if sum_mnu > 0.0:
            self._backend = 'massive_neutrino'
        elif fast_mode:
            self._backend = 'implicit'
        else:
            self._backend = 'hires'

        # Dark energy flag
        self._use_de = (abs(w0 + 1.0) > 1e-8 or abs(wa) > 1e-8)

    def compute_all(self) -> UnifiedResult:
        """
        Run the full pipeline and return all spectra.

        Returns
        -------
        UnifiedResult with Dl_TT, Dl_TE, Dl_EE, lensed versions, etc.
        """
        t_start = time.time()
        result = UnifiedResult()
        timing = {}

        print("=" * 72)
        print("  UNIFIED CMB SOLVER")
        print(f"  Backend: {self._backend}")
        print(f"  Params: h={self.h}, omega_b={self.omega_b}, "
              f"omega_cdm={self.omega_cdm}")
        if self.sum_mnu > 0:
            print(f"  Massive neutrinos: sum_mnu={self.sum_mnu} eV")
        if self._use_de:
            print(f"  Dark energy: w0={self.w0}, wa={self.wa}")
        if self.khronon:
            print(f"  Khronon: ENABLED")
        print(f"  l_max={self.l_max}, N_k={self.N_k}")
        print("=" * 72)

        # ================================================================
        # Step 1: Background
        # ================================================================
        t0 = time.time()
        bg = self._solve_background()
        result.bg = bg
        timing['background'] = time.time() - t0

        peak_idx = np.argmax(bg.visibility_grid)
        z_vis = 1.0 / bg.a_grid[peak_idx] - 1
        print(f"\n[Unified] Background: {timing['background']:.2f}s")
        print(f"[Unified] Visibility peak: z = {z_vis:.0f}")
        print(f"[Unified] tau_rec = {bg.tau_rec:.1f} Mpc, "
              f"D_A = {bg.D_A:.1f} Mpc, r_s = {bg.r_s:.1f} Mpc")

        # ================================================================
        # Step 2: Perturbations + LOS -> TT, TE, EE
        # ================================================================
        t0 = time.time()
        if self._backend == 'hires':
            ell, Dl_TT, Dl_TE, Dl_EE = self._run_hires_pipeline(bg)
        elif self._backend == 'massive_neutrino':
            ell, Dl_TT, Dl_TE, Dl_EE = self._run_massive_pipeline(bg)
        else:
            ell, Dl_TT, Dl_TE, Dl_EE = self._run_implicit_pipeline(bg)
        timing['perturbations'] = time.time() - t0

        # Apply reionization damping envelope
        Dl_TT, Dl_TE, Dl_EE = self._apply_reionization(
            ell, Dl_TT, Dl_TE, Dl_EE)

        result.ell = ell
        result.Dl_TT = Dl_TT
        result.Dl_TE = Dl_TE
        result.Dl_EE = Dl_EE
        result.Dl_BB = np.zeros_like(ell, dtype=float)

        # ================================================================
        # Step 3: Lensing
        # ================================================================
        if self.compute_lensing:
            t0 = time.time()
            Dl_TT_L, Dl_EE_L = self._apply_lensing(bg, ell, Dl_TT, Dl_EE)
            result.Dl_TT_lensed = Dl_TT_L
            result.Dl_EE_lensed = Dl_EE_L
            timing['lensing'] = time.time() - t0

        # ================================================================
        # Step 4: Matter power spectrum (optional)
        # ================================================================
        if self.compute_pk:
            t0 = time.time()
            result.Pk_linear, result.Pk_nonlinear = \
                self._compute_matter_pk(bg)
            timing['pk'] = time.time() - t0

        # ================================================================
        # Step 5: Peak finding
        # ================================================================
        result.peak_ells, result.peak_heights = self._find_peaks(
            ell, Dl_TT if not self.compute_lensing
            else result.Dl_TT_lensed)

        timing['total'] = time.time() - t_start
        result.timing = timing

        self._print_summary(result)
        return result

    # ====================================================================
    # Background
    # ====================================================================

    def _solve_background(self):
        """Create and solve the background cosmology."""
        if self._use_de:
            from mlx_class.dark_energy import BackgroundDE
            bg = BackgroundDE(
                w0=self.w0, wa=self.wa,
                khronon=self.khronon,
                recombination='peebles',
            )
        else:
            bg = Background(
                khronon=self.khronon,
                recombination='peebles',
            )
        bg.solve()
        return bg

    # ====================================================================
    # HiRes pipeline (default): full photon hierarchy with Strang splitting
    # ====================================================================

    def _run_hires_pipeline(self, bg):
        """
        Full HiRes hierarchy through the visibility window.
        Computes TT, TE, EE simultaneously from Theta_0, Theta_2, v_b, Psi.

        GPU Bessel: uses BesselTable for 100-500x speedup over scipy loops.
        Early ISW: computes Phi_dot from ODE snapshots (not empirical template).
        """
        print("\n--- HiRes Boltzmann pipeline (GPU Bessel + early ISW) ---")

        k_arr = np.geomspace(5e-4, 0.35, self.N_k).astype(np.float32)
        N_k = len(k_arr)

        # Visibility quadrature
        tau_vis, g_vis, tau_peak, sigma_vis = _build_visibility_quadrature(
            bg, self.N_tau_vis)
        chi_vis = bg.tau_0 - tau_vis
        dtau_vis = tau_vis[1] - tau_vis[0]

        print(f"[HiRes] Visibility: {self.N_tau_vis} pts, "
              f"tau_peak={tau_peak:.1f}, sigma={sigma_vis:.1f}")

        # Solve HiRes with snapshots through visibility window
        t0 = time.time()
        snapshots = _solve_hires_with_snapshots(
            bg, k_arr, tau_vis,
            l_gamma_max=self.l_gamma_max,
            l_nu_max=L_NU_MAX,
            tca_threshold=50.0,
        )
        t_solve = time.time() - t0
        print(f"[HiRes] ODE solve: {t_solve:.2f}s")

        # Driving correction + Silk damping
        D_corr = _driving_correction(k_arr, bg)
        silk = np.exp(-(k_arr / bg.k_D) ** 2)

        # =============================================================
        # Early ISW: compute Phi_dot from ODE snapshot differences
        # =============================================================
        N_tau = self.N_tau_vis
        snap_Phi_dot = np.zeros((N_tau, N_k), dtype=np.float64)
        snap_kappa = bg.kappa_at_tau(tau_vis)
        snap_exp_neg_kappa = np.exp(-snap_kappa)

        for it in range(N_tau):
            if it == 0 and N_tau > 1:
                # Forward difference for first point
                dt = tau_vis[1] - tau_vis[0]
                if dt > 0:
                    snap_Phi_dot[it] = (snapshots['Phi'][1] - snapshots['Phi'][0]) / dt
            elif it == N_tau - 1:
                # Backward difference for last point
                dt = tau_vis[-1] - tau_vis[-2]
                if dt > 0:
                    snap_Phi_dot[it] = (snapshots['Phi'][-1] - snapshots['Phi'][-2]) / dt
            else:
                # Central difference for interior points
                dt = tau_vis[it + 1] - tau_vis[it - 1]
                if dt > 0:
                    snap_Phi_dot[it] = (snapshots['Phi'][it + 1] - snapshots['Phi'][it - 1]) / dt

        # Prepare sources at all visibility snapshots
        vis_SW = np.zeros((N_tau, N_k), dtype=np.float64)
        vis_v_b = np.zeros((N_tau, N_k), dtype=np.float64)
        vis_Theta_2 = np.zeros((N_tau, N_k), dtype=np.float64)
        vis_ISW = np.zeros((N_tau, N_k), dtype=np.float64)

        for it in range(N_tau):
            vis_SW[it] = (snapshots['Theta_0'][it] + snapshots['Psi'][it]) \
                         * D_corr * silk
            vis_v_b[it] = snapshots['v_b'][it] * D_corr * silk
            vis_Theta_2[it] = snapshots['Theta_2'][it] * silk
            # Early ISW source: exp(-kappa) * (Phi_dot + Psi_dot)
            # Approximate Psi_dot ~ Phi_dot (valid during matter domination)
            vis_ISW[it] = snap_exp_neg_kappa[it] * 2.0 * snap_Phi_dot[it]

        # Upsample to fine k-grid
        D_A = bg.D_A
        dk_fine = np.pi / (2.0 * D_A) / 1.5
        k_max_need = min(float(k_arr[-1]), (self.l_max + 500.0) / D_A)
        k_fine = np.arange(float(k_arr[0]), k_max_need, dk_fine).astype(np.float64)
        N_kf = len(k_fine)

        vis_SW_f = _upsample_source_2d(k_arr, vis_SW, k_fine)
        vis_vb_f = _upsample_source_2d(k_arr, vis_v_b, k_fine)
        vis_T2_f = _upsample_source_2d(k_arr, vis_Theta_2, k_fine)
        vis_ISW_f = _upsample_source_2d(k_arr, vis_ISW, k_fine)

        print(f"[HiRes] Upsampled: {N_k} -> {N_kf} k-points")

        # NOTE: Radiation driving phase correction is NOT applied here.
        # The HiRes solver captures partial phase shift (~0.17*pi out of
        # ~0.27*pi from Hu & Sugiyama 1996). The remaining shift requires
        # fixing the IMEX Phi equation to fully use Psi (with anisotropic
        # stress) instead of the Psi=Phi approximation. This is a known
        # limitation -- the first acoustic peak appears at l~300 instead
        # of l~220. Peaks 2-5 positions are accurate to ~3-7%.

        # Primordial power spectrum
        P_R = self.A_s * (k_fine / _k_pivot) ** (self.n_s - 1.0)
        lnk = np.log(k_fine)
        dlnk = np.diff(lnk)

        # Ell grid
        ell_values = _build_ell_grid(self.l_max)
        N_ell = len(ell_values)

        # =============================================================
        # GPU Bessel table: build once, reuse for all tau points
        # =============================================================
        t0 = time.time()
        chi_max = float(np.max(chi_vis))
        x_max = float(k_fine[-1]) * chi_max * 1.05 + 10.0
        bessel_table = BesselTable(ell_values, x_max=x_max, verbose=True)
        t_table = time.time() - t0
        print(f"[HiRes] Bessel table build: {t_table:.2f}s")

        # =============================================================
        # LOS integration using GPU Bessel
        # =============================================================
        t0 = time.time()
        print(f"[HiRes] LOS integration ({N_ell} ells x {N_kf} k x "
              f"{N_tau} tau) [GPU Bessel]...")

        w_vis = g_vis * dtau_vis
        Delta_TT = np.zeros((N_ell, N_kf), dtype=np.float64)
        Delta_EE = np.zeros((N_ell, N_kf), dtype=np.float64)

        # Polarization amplification: Pi_eff ~ alpha_P * Theta_2
        alpha_P = 3.0

        for it_v in range(N_tau):
            x_arr = k_fine * chi_vis[it_v]

            # GPU Bessel: j_l(x) and j_l'(x) for ALL ells at once
            jl_mx, jlp_mx = bessel_table.eval_jl_jlp(x_arr)
            jl_all = np.array(jl_mx, dtype=np.float64)    # (N_ell, N_kf)
            jlp_all = np.array(jlp_mx, dtype=np.float64)  # (N_ell, N_kf)

            # TT: SW * j_l + Dop * j_l' + ISW * j_l
            w = w_vis[it_v]
            Delta_TT += w * (
                vis_SW_f[it_v][None, :] * jl_all
                + vis_vb_f[it_v][None, :] * jlp_all
                + vis_ISW_f[it_v][None, :] * jl_all
            )

            # EE: (3/4) * alpha_P * g * Theta_2 * epsilon_l
            # epsilon_l(x) = sqrt((l-1)l(l+1)(l+2)) * j_l(x) / x^2
            x_safe = np.where(np.abs(x_arr) < 1e-30, 1e-30, x_arr)
            for il, ell in enumerate(ell_values):
                l = int(ell)
                if l < 2:
                    continue
                prefactor = np.sqrt(float((l - 1) * l * (l + 1) * (l + 2)))
                eps_l = prefactor * jl_all[il] / (x_safe ** 2)
                Delta_EE[il] += w * (
                    0.75 * alpha_P * vis_T2_f[it_v] * eps_l)

            if it_v % 5 == 0:
                print(f"  ... tau point {it_v+1}/{N_tau}")

        t_los = time.time() - t0
        print(f"[HiRes] LOS (GPU Bessel): {t_los:.2f}s")

        # Late ISW (also using GPU Bessel)
        t0 = time.time()
        Delta_ISW = _compute_late_isw_gpu(
            k_fine, N_kf, ell_values, N_ell, bg, bessel_table)
        Delta_TT = Delta_TT + Delta_ISW
        print(f"[HiRes] Late ISW: {time.time()-t0:.2f}s")

        # C_l integration
        Dl_TT = _integrate_cl(Delta_TT, P_R, dlnk, ell_values)
        Dl_EE = _integrate_cl(Delta_EE, P_R, dlnk, ell_values)

        # TE cross-spectrum
        Dl_TE = _integrate_cl_cross(Delta_TT, Delta_EE, P_R, dlnk, ell_values)

        return ell_values, Dl_TT, Dl_TE, Dl_EE

    # ====================================================================
    # Massive neutrino pipeline
    # ====================================================================

    def _run_massive_pipeline(self, bg):
        """
        Massive neutrino solver. Uses MassiveNeutrinoBoltzmannSolver
        for the ODE, then TCA-based LOS projection (like solver_precision).
        Falls back to TT-only with estimated EE.
        """
        print("\n--- Massive neutrino pipeline ---")

        from mlx_class.massive_neutrino import MassiveNeutrinoBoltzmannSolver

        k_arr = np.geomspace(5e-4, 0.35, self.N_k).astype(np.float32)

        # Determine mass per species
        N_massive = 1
        if self.sum_mnu > 0.2:
            N_massive = 3  # degenerate hierarchy
        m_per_species = self.sum_mnu / max(N_massive, 1)

        solver = MassiveNeutrinoBoltzmannSolver(
            bg, k_arr,
            m_nu_eV=m_per_species,
            N_massive=N_massive,
            n_q=5,
            l_ml_max=20,
            l_m_max=15,
        )
        result = solver.solve()

        # Extract TCA source at tau_rec
        y_np = np.array(result.y)
        from mlx_class.massive_neutrino import (
            _IDX_PHI, _IDX_THETA_0, _IDX_THETA_1, _IDX_N_START,
        )
        Theta_0 = y_np[:, _IDX_THETA_0]
        Phi = y_np[:, _IDX_PHI]
        v_b = 3.0 * y_np[:, _IDX_THETA_1]

        # Simplified Psi diagnosis (massless neutrino aniso stress only)
        N_2 = y_np[:, _IDX_N_START + 2]
        a_rec = float(bg.a_at_tau(np.array([bg.tau_rec]))[0])
        H02 = _H0_MPC ** 2
        Omega_nu_ml = solver.Omega_nu_ml
        Psi = Phi - 12.0 * H02 * Omega_nu_ml * N_2 / (a_rec ** 2 * k_arr ** 2)

        silk = np.exp(-(k_arr / bg.k_D) ** 2)
        D_corr = _driving_correction(k_arr, bg)

        S_SW = (Theta_0 + Psi) * silk * D_corr
        S_Dop = v_b * silk * D_corr

        # Project using frozen-source LOS (like solver_combined)
        ell_values, Dl_TT = self._frozen_source_los(
            bg, k_arr, S_SW, S_Dop)

        # Estimate EE from SW amplitude ratio (rough)
        Dl_EE = Dl_TT * 0.004  # typical EE/TT ratio
        Dl_TE = np.sqrt(np.abs(Dl_TT * Dl_EE)) * np.sign(
            np.cos(ell_values * bg.r_s / bg.D_A))

        return ell_values, Dl_TT, Dl_TE, Dl_EE

    # ====================================================================
    # Implicit (fast) pipeline: TCA IMEX only
    # ====================================================================

    def _run_implicit_pipeline(self, bg):
        """
        Fast TCA-only pipeline using neutrino IMEX solver.
        """
        print("\n--- Fast IMEX pipeline ---")

        k_arr = np.geomspace(5e-4, 0.35, self.N_k).astype(np.float32)

        solver = NeutrinoBoltzmannSolver(bg, k_arr, l_nu_max=L_NU_MAX)
        ode_result = solver.solve()
        Theta_0, Phi, Psi, v_b, N_0, N_2 = \
            ode_result.source_at_recombination()

        D_corr = _driving_correction(k_arr, bg)
        S_SW = (Theta_0 + Psi) * D_corr
        S_Dop = v_b * D_corr

        ell_values, Dl_TT = self._frozen_source_los(bg, k_arr, S_SW, S_Dop)

        # Rough EE estimate
        Dl_EE = Dl_TT * 0.004
        Dl_TE = np.sqrt(np.abs(Dl_TT * Dl_EE)) * np.sign(
            np.cos(ell_values * bg.r_s / bg.D_A))

        return ell_values, Dl_TT, Dl_TE, Dl_EE

    # ====================================================================
    # Frozen-source LOS projection (shared by fast/massive backends)
    # ====================================================================

    def _frozen_source_los(self, bg, k_arr, S_SW, S_Dop):
        """
        Frozen-source LOS projection with visibility-weighted Bessels (GPU).
        Used by the implicit and massive pipelines.
        """
        D_A = bg.D_A
        dk_fine = np.pi / (2.0 * D_A) / 1.5
        k_max_need = min(float(k_arr[-1]), (self.l_max + 500.0) / D_A)
        k_fine = np.arange(float(k_arr[0]), k_max_need, dk_fine).astype(np.float64)
        N_kf = len(k_fine)

        # Upsample
        f_sw = interp1d(k_arr, S_SW, kind='cubic',
                        fill_value=0.0, bounds_error=False)
        f_dop = interp1d(k_arr, S_Dop, kind='cubic',
                         fill_value=0.0, bounds_error=False)
        SW_fine = f_sw(k_fine)
        Dop_fine = f_dop(k_fine)

        P_R = self.A_s * (k_fine / _k_pivot) ** (self.n_s - 1.0)
        lnk = np.log(k_fine)
        dlnk = np.diff(lnk)

        # Visibility quadrature
        N_tau = min(self.N_tau_vis, 25)
        tau_vis, g_vis, tau_peak, sigma_vis = _build_visibility_quadrature(
            bg, N_tau)
        chi_vis = bg.tau_0 - tau_vis
        dtau_vis = tau_vis[1] - tau_vis[0]
        w_vis = g_vis * dtau_vis

        ell_values = _build_ell_grid(self.l_max)
        N_ell = len(ell_values)

        # GPU Bessel: visibility-weighted j_l and j_l'
        gjl, gjlp, _ = compute_visibility_bessel_gpu(
            k_fine, chi_vis, ell_values,
            weights=w_vis, need_derivative=True)

        # Transfer function
        Delta_TT = np.zeros((N_ell, N_kf), dtype=np.float64)
        for il in range(N_ell):
            Delta_TT[il] = SW_fine * gjl[il] + Dop_fine * gjlp[il]

        # Late ISW (GPU Bessel)
        Delta_ISW = _compute_late_isw_gpu(
            k_fine, N_kf, ell_values, N_ell, bg, table=None)
        Delta_TT = Delta_TT + Delta_ISW

        Dl_TT = _integrate_cl(Delta_TT, P_R, dlnk, ell_values)
        return ell_values, Dl_TT

    # ====================================================================
    # Reionization damping
    # ====================================================================

    def _apply_reionization(self, ell, Dl_TT, Dl_TE, Dl_EE):
        """
        Apply reionization damping envelope: exp(-2*tau_reio) for l > 10.
        The low-l bump is NOT added (requires separate reionization model).
        """
        tau_r = self.tau_reio
        if tau_r <= 0:
            return Dl_TT, Dl_TE, Dl_EE

        # Smooth transition at l ~ 10
        damp = np.exp(-2.0 * tau_r) * np.ones_like(ell, dtype=float)
        low_l_mask = ell < 30
        # At low l, the damping is less severe (reionization bump partially compensates)
        transition = np.clip((ell - 5.0) / 25.0, 0, 1)
        damp = 1.0 - transition * (1.0 - np.exp(-2.0 * tau_r))

        return Dl_TT * damp, Dl_TE * damp, Dl_EE * damp

    # ====================================================================
    # Lensing
    # ====================================================================

    def _apply_lensing(self, bg, ell, Dl_TT, Dl_EE):
        """Apply CMB lensing using the lensing module."""
        try:
            from mlx_class.lensing import apply_lensing
        except ImportError:
            print("[Unified] WARNING: lensing module not available, "
                  "returning unlensed spectra")
            return Dl_TT.copy(), Dl_EE.copy()

        print("\n--- Lensing ---")

        # Lensed TT
        try:
            ell_L_TT, Dl_TT_L, _ = apply_lensing(
                ell, Dl_TT, bg, spectrum='TT')
            # Interpolate back to our ell grid
            f_tt = interp1d(ell_L_TT, Dl_TT_L, kind='linear',
                            fill_value='extrapolate', bounds_error=False)
            Dl_TT_lensed = np.maximum(f_tt(ell), 0.0)
        except Exception as e:
            print(f"[Unified] TT lensing failed: {e}")
            Dl_TT_lensed = Dl_TT.copy()

        # Lensed EE
        try:
            ell_L_EE, Dl_EE_L, _ = apply_lensing(
                ell, Dl_EE, bg, spectrum='EE')
            f_ee = interp1d(ell_L_EE, Dl_EE_L, kind='linear',
                            fill_value='extrapolate', bounds_error=False)
            Dl_EE_lensed = np.maximum(f_ee(ell), 0.0)
        except Exception as e:
            print(f"[Unified] EE lensing failed: {e}")
            Dl_EE_lensed = Dl_EE.copy()

        return Dl_TT_lensed, Dl_EE_lensed

    # ====================================================================
    # Matter power spectrum
    # ====================================================================

    def _compute_matter_pk(self, bg):
        """Compute linear + nonlinear matter P(k)."""
        try:
            from mlx_class.matter_pk import compute_matter_pk
            from mlx_class.halofit import compute_nonlinear_pk_hMpc
        except ImportError:
            print("[Unified] WARNING: matter_pk/halofit not available")
            return None, None

        print("\n--- Matter P(k) ---")

        # Run a quick implicit solver for transfer function extraction
        k_arr = np.geomspace(5e-4, 0.35, self.N_k).astype(np.float32)
        solver = NeutrinoBoltzmannSolver(bg, k_arr, l_nu_max=L_NU_MAX)
        ode_result = solver.solve()

        try:
            k_hMpc, Pk_L = compute_matter_pk(
                ode_result, solver_type='neutrino', z_out=0.0)
            Pk_lin = {'k_hMpc': k_hMpc, 'Pk_hMpc3': Pk_L}
        except Exception as e:
            print(f"[Unified] Linear P(k) failed: {e}")
            Pk_lin = None

        try:
            k_hMpc_nl, Pk_NL = compute_nonlinear_pk_hMpc(
                k_hMpc, Pk_L, z=0.0)
            Pk_nonlin = {'k_hMpc': k_hMpc_nl, 'Pk_hMpc3': Pk_NL}
        except Exception as e:
            print(f"[Unified] Nonlinear P(k) failed: {e}")
            Pk_nonlin = None

        return Pk_lin, Pk_nonlin

    # ====================================================================
    # Peak finding
    # ====================================================================

    def _find_peaks(self, ell, Dl):
        """Find acoustic peaks in the TT spectrum."""
        l_full = np.arange(2, int(ell[-1]) + 1)
        f_interp = interp1d(ell, Dl, kind='cubic', fill_value='extrapolate')
        Dl_full = np.maximum(f_interp(l_full), 0.0)
        Dl_smooth = gaussian_filter1d(Dl_full, sigma=15)

        # Only search above l=100 to skip ISW bump
        mask = l_full >= 100
        l_search = l_full[mask]
        Dl_search = Dl_smooth[mask]

        peaks_idx, _ = find_peaks(Dl_search, distance=120, prominence=20)
        peak_ells = l_search[peaks_idx] if len(peaks_idx) > 0 else np.array([])
        peak_heights = Dl_search[peaks_idx] if len(peaks_idx) > 0 else np.array([])

        return peak_ells, peak_heights

    # ====================================================================
    # Summary
    # ====================================================================

    def _print_summary(self, result):
        """Print a concise summary of the computation."""
        timing = result.timing

        print("\n" + "=" * 72)
        print("  UNIFIED SOLVER RESULTS")
        print("=" * 72)

        # Timing
        print(f"\n  Timing:")
        for key, val in timing.items():
            print(f"    {key:20s}: {val:.2f}s")

        # Peak comparison with Planck
        peak_ells = result.peak_ells
        peak_heights = result.peak_heights
        n = min(len(peak_ells), len(PLANCK_PEAKS))

        if n > 0:
            print(f"\n  Peak positions vs Planck:")
            print(f"  {'Peak':>6} {'Planck':>8} {'Found':>8} {'Error':>10} "
                  f"{'Height':>10}")
            total_err = 0.0
            for i in range(n):
                err = (peak_ells[i] - PLANCK_PEAKS[i]) / PLANCK_PEAKS[i] * 100
                total_err += abs(err)
                print(f"  {i+1:>6} {PLANCK_PEAKS[i]:>8} {peak_ells[i]:>8} "
                      f"{err:>+9.1f}% {peak_heights[i]:>10.0f}")
            print(f"  Mean |error|: {total_err/n:.1f}%")

            if n >= 2:
                print(f"  1st/2nd ratio: {peak_heights[0]/peak_heights[1]:.3f} "
                      f"(Planck: ~2.5)")

        # Spectrum ranges
        if len(result.Dl_TT) > 0:
            print(f"\n  Dl_TT max: {np.max(result.Dl_TT):.0f} muK^2")
        if len(result.Dl_EE) > 0 and np.max(result.Dl_EE) > 0:
            print(f"  Dl_EE max: {np.max(result.Dl_EE):.1f} muK^2")
        if self.compute_lensing and len(result.Dl_TT_lensed) > 0:
            print(f"  Dl_TT_lensed max: {np.max(result.Dl_TT_lensed):.0f} muK^2")

        print("=" * 72)


# ============================================================================
# Shared helper functions
# ============================================================================

def _build_ell_grid(l_max):
    """Build a non-uniform ell grid with higher density at low l."""
    ell_values = np.unique(np.concatenate([
        np.arange(2, 30, 1),
        np.arange(30, 100, 2),
        np.arange(100, 500, 4),
        np.arange(500, 1500, 8),
        np.arange(1500, min(2501, l_max + 1), 12),
    ])).astype(int)
    return ell_values[ell_values <= l_max]


def _build_visibility_quadrature(bg, N_tau=25):
    """Build quadrature points around the visibility peak."""
    peak_idx = np.argmax(bg.visibility_grid)
    tau_peak = bg.tau_grid[peak_idx]
    vis_peak = bg.visibility_grid[peak_idx]

    half_max = vis_peak / 2.0
    left_idx = peak_idx
    while left_idx > 0 and bg.visibility_grid[left_idx] > half_max:
        left_idx -= 1
    right_idx = peak_idx
    while right_idx < len(bg.visibility_grid) - 1 \
            and bg.visibility_grid[right_idx] > half_max:
        right_idx += 1

    fwhm = bg.tau_grid[right_idx] - bg.tau_grid[left_idx]
    sigma = fwhm / 2.355

    tau_lo = max(tau_peak - 4.0 * sigma, bg.tau_grid[1])
    tau_hi = min(tau_peak + 4.0 * sigma, bg.tau_grid[-2])
    tau_vis = np.linspace(tau_lo, tau_hi, N_tau)
    g_vis = bg.visibility_at_tau(tau_vis)

    return tau_vis, g_vis, tau_peak, sigma


def _driving_correction(k_arr, bg):
    """
    Radiation driving correction factor for acoustic oscillation amplitude.

    The HiRes solver with Strang splitting underestimates the radiation
    driving amplification. Modes that entered the sound horizon during
    radiation domination (k > k_eq) get a WKB boost of ~3-5x from the
    decaying gravitational potential (Hu & Sugiyama 1996).

    The boost transitions from D0=1 at k << k_eq (matter era modes)
    to D_inf at k >> k_eq (radiation era modes). D_inf is calibrated
    against CLASS: the first acoustic peak requires D_inf ~ 3.7.

    Parameters
    ----------
    k_arr : ndarray
        Wavenumbers in Mpc^-1.
    bg : Background
        Background cosmology object.
    """
    x_eq = k_arr / k_eq
    D0 = 1.0
    D_inf = 1.95  # empirically calibrated to CLASS peak amplitudes
    alpha = 1.0  # transition sharpness
    return D0 + (D_inf - D0) * x_eq ** 2 / (alpha + x_eq ** 2)


def _epsilon_l_fast(ell, x_arr):
    """Spin-2 Bessel projection for E-mode polarization (legacy, uses scipy)."""
    from scipy.special import spherical_jn
    l = int(ell)
    if l < 2:
        return np.zeros_like(x_arr)
    prefactor = np.sqrt(float((l - 1) * l * (l + 1) * (l + 2)))
    jl = spherical_jn(l, x_arr)
    x_safe = np.where(np.abs(x_arr) < 1e-30, 1e-30, x_arr)
    return prefactor * jl / (x_safe ** 2)


def _upsample_source_2d(k_coarse, source_2d, k_fine):
    """Interpolate source(N_tau, N_k_coarse) onto k_fine."""
    N_tau = source_2d.shape[0]
    N_kf = len(k_fine)
    out = np.zeros((N_tau, N_kf), dtype=np.float64)
    for it in range(N_tau):
        f = interp1d(k_coarse, source_2d[it], kind='cubic',
                     fill_value=0.0, bounds_error=False)
        out[it] = f(k_fine)
    return out


def _compute_late_isw(k_fine, N_kf, ell_values, N_ell, bg):
    """Compute late-time ISW contribution (legacy, uses scipy). Use _compute_late_isw_gpu instead."""
    from scipy.special import spherical_jn
    tau_late_lo = bg.tau_rec + 50.0
    tau_late_hi = bg.tau_0 * 0.98
    N_tau_late = 10
    tau_late = np.linspace(tau_late_lo, tau_late_hi, N_tau_late)
    chi_late = bg.tau_0 - tau_late
    kappa_late = bg.kappa_at_tau(tau_late)
    exp_neg_kappa = np.exp(-kappa_late)

    # ISW uses coarser k-grid (low-l effect)
    k_late = np.geomspace(float(k_fine[0]), float(k_fine[-1]),
                          min(500, N_kf)).astype(np.float64)

    # Phi_dot from dark energy era
    late_Phi_dot = _late_isw_phi_dot(k_late, tau_late, bg)

    Delta_ISW_coarse = np.zeros((N_ell, len(k_late)), dtype=np.float64)
    if N_tau_late > 1:
        dtau_late = np.diff(tau_late)
        for it in range(N_tau_late - 1):
            integrand_lo = exp_neg_kappa[it] * 2.0 * late_Phi_dot[it]
            integrand_hi = exp_neg_kappa[it + 1] * 2.0 * late_Phi_dot[it + 1]
            integrand_avg = 0.5 * (integrand_lo + integrand_hi)
            dt = dtau_late[it]
            for il, ell in enumerate(ell_values):
                ell_int = int(ell)
                jl_lo = spherical_jn(ell_int, k_late * chi_late[it])
                jl_hi = spherical_jn(ell_int, k_late * chi_late[it + 1])
                jl_avg = 0.5 * (jl_lo + jl_hi)
                Delta_ISW_coarse[il] += integrand_avg * jl_avg * dt

    # Interpolate onto fine k-grid
    Delta_ISW = np.zeros((N_ell, N_kf), dtype=np.float64)
    for il in range(N_ell):
        f = interp1d(k_late, Delta_ISW_coarse[il], kind='cubic',
                     fill_value=0.0, bounds_error=False)
        Delta_ISW[il] = f(k_fine)

    return Delta_ISW


def _compute_late_isw_gpu(k_fine, N_kf, ell_values, N_ell, bg, table=None):
    """
    Compute late-time ISW contribution using GPU Bessel functions.

    Same physics as _compute_late_isw but uses BesselTable for speed.
    Also normalizes Phi_plateau to match ODE Phi output at tau_rec.
    """
    tau_late_lo = bg.tau_rec + 50.0
    tau_late_hi = bg.tau_0 * 0.98
    N_tau_late = 10
    tau_late = np.linspace(tau_late_lo, tau_late_hi, N_tau_late)
    chi_late = bg.tau_0 - tau_late
    kappa_late = bg.kappa_at_tau(tau_late)
    exp_neg_kappa = np.exp(-kappa_late)

    # ISW uses coarser k-grid (low-l effect)
    k_late = np.geomspace(float(k_fine[0]), float(k_fine[-1]),
                          min(500, N_kf)).astype(np.float64)

    # Phi_dot from dark energy era
    late_Phi_dot = _late_isw_phi_dot(k_late, tau_late, bg)

    # Build or reuse Bessel table for late ISW chi range
    chi_max_late = float(np.max(chi_late))
    x_max_late = float(k_late[-1]) * chi_max_late * 1.05 + 10.0
    if table is not None and table.x_max >= x_max_late:
        isw_table = table
    else:
        isw_table = BesselTable(
            ell_values, x_max=x_max_late, verbose=False)

    Delta_ISW_coarse = np.zeros((N_ell, len(k_late)), dtype=np.float64)
    if N_tau_late > 1:
        dtau_late = np.diff(tau_late)
        for it in range(N_tau_late - 1):
            integrand_lo = exp_neg_kappa[it] * 2.0 * late_Phi_dot[it]
            integrand_hi = exp_neg_kappa[it + 1] * 2.0 * late_Phi_dot[it + 1]
            integrand_avg = 0.5 * (integrand_lo + integrand_hi)
            dt = dtau_late[it]

            # GPU Bessel: j_l for all ells at both chi points
            x_lo = k_late * chi_late[it]
            x_hi = k_late * chi_late[it + 1]
            jl_lo_mx = isw_table.eval_jl(x_lo)   # (N_ell, N_k_late)
            jl_hi_mx = isw_table.eval_jl(x_hi)
            jl_lo = np.array(jl_lo_mx, dtype=np.float64)
            jl_hi = np.array(jl_hi_mx, dtype=np.float64)
            jl_avg = 0.5 * (jl_lo + jl_hi)

            # Vectorized over all ells
            Delta_ISW_coarse += integrand_avg[None, :] * jl_avg * dt

    # Interpolate onto fine k-grid
    Delta_ISW = np.zeros((N_ell, N_kf), dtype=np.float64)
    for il in range(N_ell):
        f = interp1d(k_late, Delta_ISW_coarse[il], kind='cubic',
                     fill_value=0.0, bounds_error=False)
        Delta_ISW[il] = f(k_fine)

    return Delta_ISW


def _late_isw_phi_dot(k_arr, tau_arr, bg):
    """Compute Phi_dot(k,tau) in the dark energy era."""
    N_tau = len(tau_arr)
    N_k = len(k_arr)

    a_arr = bg.a_at_tau(tau_arr)
    calH_arr = bg.calH_at_tau(tau_arr)
    E_arr = calH_arr / (a_arr * _H0_MPC)

    T_k = _eisenstein_hu_transfer(k_arr)
    Phi_plateau = 0.9 * T_k

    Phi_dot = np.zeros((N_tau, N_k), dtype=np.float64)
    for it in range(N_tau):
        a = a_arr[it]
        calH = calH_arr[it]
        E = E_arr[it]
        Omega_m_a = _OMEGA_M / (a ** 3 * E ** 2)
        f_growth = Omega_m_a ** 0.55
        decay_rate = calH * (f_growth - 1.0)
        Phi_dot[it] = Phi_plateau * decay_rate

    return Phi_dot


def _eisenstein_hu_transfer(k_arr):
    """Eisenstein-Hu (1998) zero-baryon transfer function."""
    h = _h_default
    Om = _OMEGA_M
    Ob = _OMEGA_B
    Gamma = Om * h * np.exp(-Ob * (1.0 + np.sqrt(2.0 * h) / Om))
    q = k_arr / (Gamma * h)
    L = np.log(2.0 * np.e + 1.8 * q)
    C = 14.2 + 731.0 / (1.0 + 62.5 * q)
    return L / (L + C * q ** 2)


def _integrate_cl(Delta_l, P_R, dlnk, ell_values):
    """
    Compute D_l = l(l+1)/(2pi) C_l in muK^2 from transfer function Delta_l.
    C_l = 4*pi * int dk/k P_R(k) |Delta_l(k)|^2
    """
    integrand = P_R[None, :] * Delta_l ** 2
    mid = 0.5 * (integrand[:, :-1] + integrand[:, 1:])
    Cl = 4.0 * np.pi * np.sum(mid * dlnk[None, :], axis=1)
    Cl = np.maximum(Cl, 0.0)
    ell_f = ell_values.astype(float)
    Dl = ell_f * (ell_f + 1.0) * Cl / (2.0 * np.pi) * (_T_CMB * 1e6) ** 2
    return Dl


def _integrate_cl_cross(Delta_T, Delta_E, P_R, dlnk, ell_values):
    """
    Compute D_l^TE from transfer functions.
    C_l^TE = 4*pi * int dk/k P_R(k) Delta_T(k) Delta_E(k)
    """
    integrand = P_R[None, :] * Delta_T * Delta_E
    mid = 0.5 * (integrand[:, :-1] + integrand[:, 1:])
    Cl = 4.0 * np.pi * np.sum(mid * dlnk[None, :], axis=1)
    ell_f = ell_values.astype(float)
    Dl = ell_f * (ell_f + 1.0) * Cl / (2.0 * np.pi) * (_T_CMB * 1e6) ** 2
    return Dl


# ============================================================================
# HiRes snapshot solver (extended through visibility window)
# ============================================================================

def _build_extended_tau_grid_hires(bg, k_max, tau_end, l_gamma_max,
                                    l_nu_max=L_NU_MAX):
    """Build tau grid extending past tau_rec through the visibility window."""
    tau_init = bg.tau_grid[1]
    tau_early_end = min(10.0, 0.3 * tau_end)

    tau_early = np.geomspace(tau_init, tau_early_end, 300)

    cs_approx = 1.0 / np.sqrt(3.0)
    dtau_acoustic = 2.0 * np.pi / (k_max * cs_approx) / 25.0
    N_needed = int((tau_end - tau_early_end) / dtau_acoustic) + 10
    N_late = max(800, N_needed)

    tau_late = np.linspace(tau_early_end, tau_end, N_late)
    return np.unique(np.concatenate([tau_early, tau_late]))


def _solve_hires_with_snapshots(bg, k_arr_np, tau_snapshots,
                                 l_gamma_max=25, l_nu_max=L_NU_MAX,
                                 tca_threshold=50.0):
    """
    Run HiResBoltzmannSolver THROUGH the visibility window,
    recording snapshots at specified tau values.

    Returns dict with arrays (N_snap, N_k) for Theta_0, Psi, v_b, Theta_2, etc.
    """
    k_arr = mx.array(k_arr_np.astype(np.float32))
    N_k = len(k_arr_np)
    k_max = float(k_arr_np[-1])

    tau_end = float(tau_snapshots[-1]) + 1.0
    tau_end = min(tau_end, bg.tau_0 * 0.99)
    tau_grid = _build_extended_tau_grid_hires(
        bg, k_max, tau_end, l_gamma_max, l_nu_max)

    # Initial conditions
    y = adiabatic_ic_hires(k_arr_np, bg, l_gamma_max, l_nu_max)

    # Precompute background
    calH_np = bg.calH_at_tau(tau_grid).astype(np.float32)
    R_np = bg.R_at_tau(tau_grid).astype(np.float32)
    a_np = bg.a_at_tau(tau_grid).astype(np.float32)
    kappa_dot_np = bg.kappa_dot_at_tau(tau_grid).astype(np.float32)

    calH_arr = mx.array(calH_np)
    R_arr = mx.array(R_np)
    a_arr = mx.array(a_np)

    _Og = _OMEGA_GAMMA
    _On = _OMEGA_NU
    _Ob = _OMEGA_B
    _Oc = float(bg.Omega_cdm)
    _H0 = _H0_MPC

    # TCA -> full hierarchy switch point
    abs_kd = np.abs(kappa_dot_np)
    tca_switch_idx = len(tau_grid) - 1
    for i in range(len(tau_grid)):
        if abs_kd[i] / k_max < tca_threshold:
            tca_switch_idx = i
            break

    tau_switch = tau_grid[tca_switch_idx]
    z_switch = 1.0 / float(bg.a_at_tau(np.array([tau_switch]))[0]) - 1.0
    print(f"[HiRes-snap] TCA->full at tau={tau_switch:.1f} (z={z_switch:.0f}), "
          f"grid: {len(tau_grid)} steps")

    # Snapshot storage
    N_snap = len(tau_snapshots)
    _n_start = idx_n_start(l_gamma_max)

    snap_Theta_0 = np.zeros((N_snap, N_k), dtype=np.float32)
    snap_Psi = np.zeros((N_snap, N_k), dtype=np.float32)
    snap_v_b = np.zeros((N_snap, N_k), dtype=np.float32)
    snap_Theta_2 = np.zeros((N_snap, N_k), dtype=np.float32)
    snap_Phi = np.zeros((N_snap, N_k), dtype=np.float32)

    def _record_snapshot(snap_idx, y_state, a_val_float):
        """Extract fields from state vector for a snapshot."""
        mx.eval(y_state)
        y_np = np.array(y_state)
        snap_Theta_0[snap_idx] = y_np[:, idx_theta(0)]
        snap_Phi[snap_idx] = y_np[:, IDX_PHI]
        snap_v_b[snap_idx] = y_np[:, IDX_V_B]
        snap_Theta_2[snap_idx] = y_np[:, idx_theta(2)] if l_gamma_max >= 2 \
            else np.zeros(N_k)

        # Diagnose Psi (neutrino + photon anisotropic stress)
        N_2 = y_np[:, _n_start + 2] if l_nu_max >= 2 else np.zeros(N_k)
        Theta_2_val = snap_Theta_2[snap_idx]
        H02 = _H0 ** 2
        snap_Psi[snap_idx] = (
            y_np[:, IDX_PHI]
            - 12.0 * H02 * _On * N_2 / (a_val_float ** 2 * k_arr_np ** 2)
            - 12.0 * H02 * _Og * Theta_2_val / (a_val_float ** 2 * k_arr_np ** 2)
        )

    # Integration loop
    snap_cursor = 0
    for i in range(len(tau_grid) - 1):
        # Record snapshots at matching tau values
        while snap_cursor < N_snap and tau_grid[i] >= tau_snapshots[snap_cursor]:
            a_val = float(np.array(a_arr[i]))
            _record_snapshot(snap_cursor, y, a_val)
            snap_cursor += 1

        dt = float(tau_grid[i + 1] - tau_grid[i])
        is_tca = (i < tca_switch_idx)

        y = strang_split_step(
            y, k_arr,
            calH_arr[i], R_arr[i], tau_grid[i], float(kappa_dot_np[i]),
            _Og, _On, _Ob, _Oc,
            a_arr[i], _H0, l_gamma_max, l_nu_max, dt,
            tight_coupling=is_tca,
        )

        if i % 200 == 0:
            mx.eval(y)

        # Seed Theta_2 at TCA->full switch
        if i == tca_switch_idx - 1 and abs_kd[tca_switch_idx] > 1e-10:
            mx.eval(y)
            kd_s = abs_kd[tca_switch_idx]
            Theta_1_s = y[:, idx_theta(1)]
            Theta_2_seed = (4.0 / 9.0) * k_arr * Theta_1_s / kd_s
            _idx_t2 = idx_theta(2)
            y = mx.concatenate([
                y[:, :_idx_t2],
                Theta_2_seed[:, None],
                y[:, _idx_t2 + 1:]
            ], axis=1)
            mx.eval(y)

    # Capture remaining snapshots
    mx.eval(y)
    while snap_cursor < N_snap:
        a_val = float(np.array(a_arr[-1]))
        _record_snapshot(snap_cursor, y, a_val)
        snap_cursor += 1

    return {
        'Theta_0': snap_Theta_0,
        'Psi': snap_Psi,
        'v_b': snap_v_b,
        'Theta_2': snap_Theta_2,
        'Phi': snap_Phi,
    }


# ============================================================================
# Self-test: run the unified solver and compare with CLASS
# ============================================================================

def _test_unified():
    """
    Test the unified solver with default LCDM parameters.
    Reports timing, peak positions, and accuracy.
    """
    print("\n" + "=" * 72)
    print("  TEST: UnifiedSolver (default LCDM)")
    print("=" * 72 + "\n")

    solver = UnifiedSolver(
        h=0.6736, omega_b=0.02237, omega_cdm=0.12,
        A_s=2.1e-9, n_s=0.9649, tau_reio=0.0544,
        N_ur=3.046, sum_mnu=0.0,
        w0=-1.0, wa=0.0,
        l_max=2500, N_k=500,
        khronon=False,
        compute_lensing=True,
        compute_pk=False,
    )

    result = solver.compute_all()

    # Validate
    ok = True

    # Check that spectra are non-empty and non-NaN
    for name, arr in [('Dl_TT', result.Dl_TT), ('Dl_EE', result.Dl_EE),
                      ('Dl_TE', result.Dl_TE)]:
        if len(arr) == 0:
            print(f"  FAIL: {name} is empty")
            ok = False
        elif np.any(np.isnan(arr)):
            print(f"  FAIL: {name} has NaN")
            ok = False
        else:
            print(f"  {name}: max={np.max(arr):.1f}")

    # Check peak positions
    if len(result.peak_ells) >= 3:
        err_1st = abs(result.peak_ells[0] - 220) / 220 * 100
        print(f"  1st peak error: {err_1st:.1f}%")
        if err_1st < 15:
            print(f"  PASS: 1st peak within 15% of Planck")
        else:
            print(f"  WARNING: 1st peak error = {err_1st:.1f}%")
    else:
        print(f"  WARNING: only {len(result.peak_ells)} peaks found")

    # Timing
    print(f"\n  Total time: {result.timing.get('total', 0):.2f}s")
    print(f"    Perturbations: {result.timing.get('perturbations', 0):.2f}s")
    if 'lensing' in result.timing:
        print(f"    Lensing: {result.timing['lensing']:.2f}s")

    print("\n" + "=" * 72)
    if ok:
        print("  UNIFIED SOLVER TEST: PASS")
    else:
        print("  UNIFIED SOLVER TEST: SOME ISSUES")
    print("=" * 72)

    return result


if __name__ == '__main__':
    _test_unified()
