"""
solver_production.py -- Production CMB Boltzmann solver.

Combines all improvements into a single, fast, accurate solver:

  FROM solver_accurate.py:
    - Correct Psi != Phi physics (F_aniso in IMEX source)
    - Poisson constraint reset (quartic blend, N=19)
    - Extended integration past tau_rec (visibility tail)
    - Snapshot storage for line-of-sight integration
    - Separate Psi computation (anisotropic stress)

  FROM solver_fast.py:
    - mx.compile for TCA and full hierarchy steps
    - Vectorized hierarchy (no Python loops over l)
    - Precomputed background quantities as MLX arrays
    - Reduced eval frequency (2-6 GPU syncs)

  FROM bessel_cache.py:
    - Disk-cached Chebyshev Bessel table
    - First run: ~4s build + save. Cached runs: ~0.02s load.
    - GPU-native evaluation via Chebyshev interpolation

  FROM radiation_driving.py:
    - D(k) amplitude correction applied per-k in C_l integral

PERFORMANCE:
  First run:  ODE ~0.5s + Bessel build ~4s + LOS ~1s  =  ~6s
  Cached run: ODE ~0.5s + Bessel load ~0.02s + LOS ~1s = ~1.5s

API:
  from mlx_class.solver_production import ProductionSolver
  solver = ProductionSolver(h=0.6736, omega_b=0.02237, omega_cdm=0.12)
  result = solver.run()
  # result.Dl_TT, result.Dl_TE, result.Dl_EE, result.peaks, result.timing

Author: Sheng-Kai Huang, 2026
"""
import numpy as np
import mlx.core as mx
import time
from dataclasses import dataclass, field
from typing import Optional, Dict, List

# ============================================================================
# Imports from existing modules
# ============================================================================

from .background import (
    Background,
    A_s, n_s, k_pivot, T_CMB,
    H0_Mpc as _H0_MPC_DEFAULT,
    Omega_r as _OMEGA_R_DEFAULT,
    Omega_b as _OMEGA_B_DEFAULT,
    Omega_c as _OMEGA_C_DEFAULT,
)

from .perturbations_neutrino import (
    N_EFF, _f_nu, _OMEGA_GAMMA, _OMEGA_NU, L_NU_MAX,
)

from .bessel_cache import CachedBesselTable

from .radiation_driving import amplitude_correction as D_k_correction

# Import physics functions from solver_accurate
from .solver_accurate import (
    IDX_PHI, IDX_DELTA_B, IDX_V_B, IDX_DELTA_C, IDX_V_C, IDX_THETA_START,
    idx_theta, idx_n_start, n_var_total,
    _phi_lambda, _phi_int_factor, _set_phi,
    _phi_source_accurate, _phi_source_accurate_tca,
    _deriv_streaming_full, _deriv_streaming_tca,
    _apply_collision_vectorized,
    _imex_rk4_step_full, _imex_rk4_step_tca,
    _strang_step_full,
    _compute_phi_constraint, _apply_constraint_reset,
    build_tau_grid, adiabatic_ic, _compute_psi,
)


# ============================================================================
# Result dataclass
# ============================================================================

@dataclass
class ProductionResult:
    """
    Production solver output.

    Attributes
    ----------
    Dl_TT : ndarray (N_ell,)
        TT power spectrum D_l = l(l+1)C_l/(2pi) in muK^2.
    Dl_TE : ndarray (N_ell,) or None
        TE cross-spectrum (placeholder, not yet implemented).
    Dl_EE : ndarray (N_ell,) or None
        EE power spectrum (placeholder, not yet implemented).
    ell : ndarray (N_ell,)
        Multipole values.
    Cl_TT : ndarray (N_ell,)
        Raw C_l (dimensionless).
    peaks : dict
        Peak positions and amplitudes: {'ell': [...], 'Dl': [...]}.
    timing : dict
        Timing breakdown: {'ode', 'bessel', 'los', 'total'}.
    k_arr : ndarray
        Wavenumber grid used.
    """
    Dl_TT: np.ndarray
    Dl_TE: Optional[np.ndarray]
    Dl_EE: Optional[np.ndarray]
    ell: np.ndarray
    Cl_TT: np.ndarray
    peaks: Dict[str, List]
    timing: Dict[str, float]
    k_arr: np.ndarray


# ============================================================================
# Production Solver
# ============================================================================

class ProductionSolver:
    """
    Production CMB power spectrum solver.

    Combines accurate physics (Psi!=Phi, constraint reset) with fast
    execution (mx.compile, cached Bessel) and D(k) correction.

    Parameters
    ----------
    h : float
        Dimensionless Hubble parameter.
    omega_b : float
        Physical baryon density.
    omega_cdm : float
        Physical CDM density.
    l_gamma_max : int
        Maximum photon multipole.
    l_nu_max : int
        Maximum neutrino multipole.
    N_k : int
        Number of k-modes.
    k_min, k_max : float
        Wavenumber range in Mpc^-1.
    ell_min, ell_max : int
        Multipole range for output.
    constraint_reset_interval : int
        Reset Poisson constraint every N steps.
    n_snapshots : int
        Number of snapshots for LOS integration.
    recombination : str
        Recombination method ('peebles' or 'tanh').
    apply_D_correction : bool
        Whether to apply D(k) radiation driving correction.
    verbose : bool
        Print progress information.
    """

    def __init__(self, h=0.6736, omega_b=0.02237, omega_cdm=0.12,
                 l_gamma_max=25, l_nu_max=L_NU_MAX,
                 N_k=500, k_min=5e-4, k_max=0.35,
                 ell_min=2, ell_max=2500,
                 constraint_reset_interval=19,
                 n_snapshots=80,
                 recombination='peebles',
                 apply_D_correction=True,
                 verbose=True):

        self.h = h
        self.omega_b = omega_b
        self.omega_cdm = omega_cdm
        self.l_gamma_max = l_gamma_max
        self.l_nu_max = l_nu_max
        self.N_k = N_k
        self.k_min = k_min
        self.k_max = k_max
        self.ell_min = ell_min
        self.ell_max = ell_max
        self.constraint_reset_interval = constraint_reset_interval
        self.n_snapshots = n_snapshots
        self.recombination = recombination
        self.apply_D_correction = apply_D_correction
        self.verbose = verbose

        # Build ell grid: dense at low-l, coarser at high-l
        self.ell_values = self._build_ell_grid(ell_min, ell_max)

        # Build k grid
        self.k_arr_np = np.geomspace(k_min, k_max, N_k).astype(np.float32)

    @staticmethod
    def _build_ell_grid(ell_min, ell_max):
        """Build non-uniform ell grid: dense near peaks, coarser at high l."""
        parts = []
        if ell_min <= 30:
            parts.append(np.arange(max(ell_min, 2), min(30, ell_max + 1), 1))
        if ell_max >= 30:
            parts.append(np.arange(30, min(100, ell_max + 1), 2))
        if ell_max >= 100:
            parts.append(np.arange(100, min(500, ell_max + 1), 4))
        if ell_max >= 500:
            parts.append(np.arange(500, min(1500, ell_max + 1), 8))
        if ell_max >= 1500:
            parts.append(np.arange(1500, ell_max + 1, 12))
        return np.unique(np.concatenate(parts)).astype(int)

    def run(self):
        """
        Execute the full pipeline: background -> ODE -> LOS -> C_l.

        Returns
        -------
        ProductionResult
            Contains Dl_TT, peaks, timing, etc.
        """
        t_total_start = time.time()
        timing = {}

        # ================================================================
        # Step 1: Background
        # ================================================================
        t0 = time.time()
        bg = Background(khronon=False, recombination=self.recombination)
        bg.solve()
        timing['background'] = time.time() - t0
        if self.verbose:
            print(f"[Production] Background: {timing['background']:.2f}s")

        # ================================================================
        # Step 2: ODE integration (from solver_accurate physics)
        # ================================================================
        t0 = time.time()
        snapshots = self._solve_ode(bg)
        timing['ode'] = time.time() - t0
        if self.verbose:
            print(f"[Production] ODE total: {timing['ode']:.2f}s")

        # ================================================================
        # Step 3: Bessel table (cached GPU)
        # ================================================================
        t0 = time.time()
        bessel_table = self._build_bessel_table(bg)
        timing['bessel'] = time.time() - t0
        if self.verbose:
            print(f"[Production] Bessel table: {timing['bessel']:.2f}s")

        # ================================================================
        # Step 4: Line-of-sight C_l with GPU Bessel
        # ================================================================
        t0 = time.time()
        ell_out, Cl, Dl = self._compute_cl_los_gpu(bg, snapshots, bessel_table)
        timing['los'] = time.time() - t0
        if self.verbose:
            print(f"[Production] LOS integration: {timing['los']:.2f}s")

        timing['total'] = time.time() - t_total_start

        # ================================================================
        # Step 5: Find peaks
        # ================================================================
        peaks = self._find_peaks(ell_out, Dl)

        if self.verbose:
            print(f"[Production] Total: {timing['total']:.2f}s")
            print(f"[Production] Peaks: {peaks['ell'][:5]}")

        return ProductionResult(
            Dl_TT=Dl,
            Dl_TE=None,
            Dl_EE=None,
            ell=ell_out,
            Cl_TT=Cl,
            peaks=peaks,
            timing=timing,
            k_arr=self.k_arr_np,
        )

    # ====================================================================
    # ODE Integration
    # ====================================================================

    def _solve_ode(self, bg):
        """
        Solve the Boltzmann hierarchy ODE with accurate physics.

        Uses solver_accurate's physics (Psi!=Phi, constraint reset)
        with solver_fast's speed (mx.compile, vectorized hierarchy).
        Returns snapshot dict for LOS integration.
        """
        k_arr_np = self.k_arr_np
        N_k = self.N_k
        k_arr = mx.array(k_arr_np)
        lgmax = self.l_gamma_max
        lnmax = self.l_nu_max
        _ns = idx_n_start(lgmax)
        n_var = n_var_total(lgmax, lnmax)

        # Time grid extending past tau_rec
        k_max_val = float(k_arr_np[-1])
        tau_grid = build_tau_grid(bg, k_max_val)

        if self.verbose:
            print(f"[Production] ODE: N_k={N_k}, l_gamma={lgmax}, "
                  f"l_nu={lnmax}, grid={len(tau_grid)} steps")

        # Initial conditions
        y = adiabatic_ic(k_arr_np, bg, lgmax, lnmax)

        # Precompute background quantities
        calH_np = bg.calH_at_tau(tau_grid).astype(np.float32)
        R_np = bg.R_at_tau(tau_grid).astype(np.float32)
        a_np = bg.a_at_tau(tau_grid).astype(np.float32)
        kappa_dot_np = bg.kappa_dot_at_tau(tau_grid).astype(np.float32)

        calH_all = mx.array(calH_np)
        R_all = mx.array(R_np)
        a_all = mx.array(a_np)
        kd_all = mx.array(kappa_dot_np)
        tau_all = mx.array(tau_grid.astype(np.float32))
        dtau_np = np.diff(tau_grid).astype(np.float32)
        dtau_all = mx.array(dtau_np)

        _Og = _OMEGA_GAMMA
        _On = _OMEGA_NU
        _Ob = _OMEGA_B_DEFAULT
        _Oc = float(bg.Omega_cdm)
        _H0 = _H0_MPC_DEFAULT

        # TCA switch point
        abs_kd = np.abs(kappa_dot_np)
        tca_threshold = 50.0
        tca_switch_idx = len(tau_grid) - 1
        for i in range(len(tau_grid)):
            if abs_kd[i] / k_max_val < tca_threshold:
                tca_switch_idx = i
                break

        # Compile step functions (closures capture k_arr, cosmological params)
        def tca_step(y, calH_i, a_i, R_i, tau_i, dtau_i):
            return _imex_rk4_step_tca(y, k_arr, calH_i, a_i, R_i, tau_i,
                                      _Og, _On, _Ob, _Oc, _H0,
                                      lgmax, lnmax, _ns, dtau_i)

        def full_step(y, calH_i, a_i, R_i, tau_i, kd_i, dtau_i):
            return _strang_step_full(y, k_arr, calH_i, a_i, R_i, tau_i, kd_i,
                                     _Og, _On, _Ob, _Oc, _H0,
                                     lgmax, lnmax, _ns, dtau_i)

        compiled_tca_step = mx.compile(tca_step)
        compiled_full_step = mx.compile(full_step)

        # Warmup compile
        y_dummy = mx.zeros_like(y)
        _ = compiled_tca_step(y_dummy, calH_all[0], a_all[0], R_all[0],
                              tau_all[0], dtau_all[0])
        mx.eval(_)
        if tca_switch_idx < len(tau_grid) - 1:
            _ = compiled_full_step(y_dummy, calH_all[0], a_all[0], R_all[0],
                                   tau_all[0], kd_all[0], dtau_all[0])
            mx.eval(_)

        # Snapshot schedule: dense around visibility peak
        N_snap = self.n_snapshots
        tau_rec = bg.tau_rec
        tau_end = tau_grid[-1]
        tau_snap_early = np.linspace(tau_grid[1], tau_rec - 50,
                                     max(10, N_snap // 6))
        tau_snap_vis = np.linspace(tau_rec - 50, tau_end,
                                   max(40, N_snap - N_snap // 6))
        tau_snapshots = np.unique(np.concatenate([tau_snap_early, tau_snap_vis]))
        N_snap = len(tau_snapshots)

        snap_grid_indices = np.searchsorted(tau_grid, tau_snapshots)
        snap_grid_indices = np.clip(snap_grid_indices, 0, len(tau_grid) - 1)

        snap_Phi = np.zeros((N_snap, N_k), dtype=np.float32)
        snap_Psi = np.zeros((N_snap, N_k), dtype=np.float32)
        snap_Theta_0 = np.zeros((N_snap, N_k), dtype=np.float32)
        snap_v_b = np.zeros((N_snap, N_k), dtype=np.float32)
        snap_tau = np.zeros(N_snap, dtype=np.float32)
        snap_a = np.zeros(N_snap, dtype=np.float32)

        grid_to_snap = {}
        for si, gi in enumerate(snap_grid_indices):
            gi_int = int(gi)
            if gi_int not in grid_to_snap:
                grid_to_snap[gi_int] = []
            grid_to_snap[gi_int].append(si)

        def take_snapshot(y_mx, grid_idx):
            if grid_idx not in grid_to_snap:
                return
            mx.eval(y_mx)
            y_np = np.array(y_mx)
            a_val = float(a_np[grid_idx])
            for si in grid_to_snap[grid_idx]:
                snap_Phi[si] = y_np[:, IDX_PHI]
                snap_Psi[si] = _compute_psi(y_np, k_arr_np, a_val,
                                             _Og, _On, _H0, _ns, lgmax)
                snap_Theta_0[si] = y_np[:, idx_theta(0)]
                snap_v_b[si] = y_np[:, IDX_V_B]
                snap_tau[si] = tau_grid[grid_idx]
                snap_a[si] = a_val

        # TCA phase with constraint resets
        N_reset = self.constraint_reset_interval
        tca_eval_interval = max(tca_switch_idx // 4, 1)
        for i in range(tca_switch_idx):
            if i in grid_to_snap:
                take_snapshot(y, i)
            y = compiled_tca_step(y, calH_all[i], a_all[i], R_all[i],
                                  tau_all[i], dtau_all[i])
            if (i + 1) % N_reset == 0:
                mx.eval(y)
                y = _apply_constraint_reset(
                    y, k_arr, calH_all[i], a_all[i],
                    _Og, _On, _Ob, _Oc, _H0, _ns, lgmax)
            if (i + 1) % tca_eval_interval == 0:
                mx.eval(y)
        mx.eval(y)

        # Seed Theta_2 at TCA -> full transition
        kd_switch = abs_kd[tca_switch_idx]
        if kd_switch > 1e-10:
            Theta_1_sw = y[:, idx_theta(1)]
            Theta_2_seed = (4.0 / 9.0) * k_arr * Theta_1_sw / kd_switch
            _idx_t2 = idx_theta(2)
            y = mx.concatenate([y[:, :_idx_t2], Theta_2_seed[:, None],
                                y[:, _idx_t2 + 1:]], axis=1)
            mx.eval(y)

        # Full hierarchy phase with constraint resets
        n_full = len(tau_grid) - 1 - tca_switch_idx
        full_eval_interval = max(n_full // 6, 1)
        for i in range(tca_switch_idx, len(tau_grid) - 1):
            if i in grid_to_snap:
                take_snapshot(y, i)
            y = compiled_full_step(y, calH_all[i], a_all[i], R_all[i],
                                   tau_all[i], kd_all[i], dtau_all[i])
            step_in_phase = i - tca_switch_idx
            if (step_in_phase + 1) % N_reset == 0:
                mx.eval(y)
                y = _apply_constraint_reset(
                    y, k_arr, calH_all[i], a_all[i],
                    _Og, _On, _Ob, _Oc, _H0, _ns, lgmax)
            if (step_in_phase + 1) % full_eval_interval == 0:
                mx.eval(y)
        mx.eval(y)

        # Final snapshot
        last_idx = len(tau_grid) - 1
        if last_idx in grid_to_snap:
            take_snapshot(y, last_idx)

        return {
            'tau': snap_tau,
            'a': snap_a,
            'Phi': snap_Phi,
            'Psi': snap_Psi,
            'Theta_0': snap_Theta_0,
            'v_b': snap_v_b,
        }

    # ====================================================================
    # Bessel Table
    # ====================================================================

    def _build_bessel_table(self, bg):
        """Build or load cached Bessel table for all ell values."""
        # x_max = k_max * (tau_0 - tau_earliest_snapshot)
        # We need x up to k_max * D_A (comoving distance to LSS)
        x_max = float(self.k_arr_np[-1]) * bg.tau_0 * 1.05 + 10.0
        return CachedBesselTable(
            self.ell_values,
            x_max=x_max,
            N_cheb=64,
            seg_width=80.0,
            verbose=self.verbose,
        )

    # ====================================================================
    # Line-of-Sight C_l with GPU Bessel + D(k) correction
    # ====================================================================

    def _compute_cl_los_gpu(self, bg, snapshots, bessel_table):
        """
        Line-of-sight C_l using GPU-cached Bessel functions.

        The source function is:
          S = D(k) * (Theta_0 + Psi) * silk * j_l(k*chi)     [Sachs-Wolfe]
            + D(k) * v_b * silk * j_l'(k*chi)                 [Doppler]
            + exp(-kappa) * (Phi'+Psi') * j_l(k*chi)          [ISW]

        D(k) is the amplitude correction from radiation_driving.py,
        applied to SW and Doppler but NOT to ISW (see that module's docs).

        The Bessel functions j_l and j_l' are evaluated on GPU via the
        cached Chebyshev table, replacing the 36s scipy loop.
        """
        k_arr = self.k_arr_np
        N_k = self.N_k
        ell_values = self.ell_values
        tau_0 = bg.tau_0

        # Extract valid snapshots
        tau_s = snapshots['tau']
        valid = tau_s > 0
        tau_s = tau_s[valid]
        Phi_s = snapshots['Phi'][valid]
        Psi_s = snapshots['Psi'][valid]
        Theta_0_s = snapshots['Theta_0'][valid]
        v_b_s = snapshots['v_b'][valid]
        N_tau = len(tau_s)

        chi_s = tau_0 - tau_s  # comoving distance for each snapshot

        # Visibility function weights
        g_s = bg.visibility_at_tau(tau_s)
        kappa_s = bg.kappa_at_tau(tau_s)
        exp_neg_kappa = np.exp(-kappa_s)

        # Phi' + Psi' from central finite differences
        PhiPsi = Phi_s + Psi_s
        PhiPsi_dot = np.zeros_like(PhiPsi)
        for i in range(1, N_tau - 1):
            dt = tau_s[i + 1] - tau_s[i - 1]
            if dt > 0:
                PhiPsi_dot[i] = (PhiPsi[i + 1] - PhiPsi[i - 1]) / dt
        if N_tau >= 2:
            dt0 = tau_s[1] - tau_s[0]
            if dt0 > 0:
                PhiPsi_dot[0] = (PhiPsi[1] - PhiPsi[0]) / dt0
            dtN = tau_s[-1] - tau_s[-2]
            if dtN > 0:
                PhiPsi_dot[-1] = (PhiPsi[-1] - PhiPsi[-2]) / dtN

        # NOTE: Silk damping is already in the ODE solution (collision terms).
        # Do NOT apply extra damping here — it was double-counting.
        silk = np.ones_like(k_arr)

        # D(k) amplitude correction
        if self.apply_D_correction:
            D_corr = D_k_correction(k_arr)
        else:
            D_corr = np.ones_like(k_arr)

        # Primordial spectrum
        P_R = A_s * (k_arr / k_pivot) ** (n_s - 1.0)
        lnk = np.log(k_arr)
        dlnk = np.diff(lnk)

        # Visibility weights: g(tau) * dtau (trapezoidal)
        dtau_s = np.diff(tau_s)
        dtau_w = np.zeros(N_tau)
        if N_tau >= 2:
            dtau_w[0] = 0.5 * dtau_s[0]
            dtau_w[-1] = 0.5 * dtau_s[-1]
            for i in range(1, N_tau - 1):
                dtau_w[i] = 0.5 * (dtau_s[i - 1] + dtau_s[i])
        w_vis = g_s * dtau_w

        if self.verbose:
            print(f"[Production] LOS: {len(ell_values)} ells, "
                  f"{N_tau} tau-steps, GPU Bessel")

        # ---------------------------------------------------------
        # Pre-evaluate ALL Bessel functions on GPU in one batch
        # ---------------------------------------------------------
        # For each tau snapshot, we need j_l(k * chi) and j_l'(k * chi)
        # x_arrays[it] = k_arr * chi_s[it], shape (N_k,)
        x_arrays = [k_arr * chi_s[it] for it in range(N_tau)]

        # Batch GPU evaluation: returns lists of mx.array, each (N_ell, N_k)
        jl_list, jlp_list = bessel_table.eval_jl_jlp_at_multiple_x(x_arrays)

        # Convert to numpy for the integration loop
        # jl_all[it] has shape (N_ell, N_k)
        jl_all = [np.array(jl) for jl in jl_list]
        jlp_all = [np.array(jlp) for jlp in jlp_list]

        # ---------------------------------------------------------
        # C_l integration: vectorized over k, loop over ell and tau
        # ---------------------------------------------------------
        N_ell = len(ell_values)
        Cl = np.zeros(N_ell, dtype=np.float64)

        # Pre-compute weighted source functions at each tau: (N_tau, N_k)
        SW_source = np.zeros((N_tau, N_k), dtype=np.float64)
        Dop_source = np.zeros((N_tau, N_k), dtype=np.float64)
        ISW_source = np.zeros((N_tau, N_k), dtype=np.float64)

        # ISW uses exp(-kappa)*dtau, NOT g(tau)*dtau
        w_isw = exp_neg_kappa * dtau_w

        for it in range(N_tau):
            # SW: D(k) * (Theta_0 + Psi) * g(tau) * dtau
            SW_source[it] = (D_corr * (Theta_0_s[it] + Psi_s[it])
                             * silk * w_vis[it])
            # Doppler: D(k) * v_b * g(tau) * dtau
            Dop_source[it] = D_corr * v_b_s[it] * silk * w_vis[it]
            # ISW: exp(-kappa) * (Phi'+Psi') * dtau  [NOT weighted by g(tau)]
            ISW_source[it] = PhiPsi_dot[it] * w_isw[it]

        for il in range(N_ell):
            Delta_l = np.zeros(N_k, dtype=np.float64)
            for it in range(N_tau):
                jl = jl_all[it][il]    # (N_k,)
                jlp = jlp_all[it][il]  # (N_k,)
                Delta_l += (SW_source[it] * jl
                            + Dop_source[it] * jlp
                            + ISW_source[it] * jl)

            # Integrate over k: int dk/k * P_R * Delta_l^2
            integrand = P_R * Delta_l ** 2
            mid = 0.5 * (integrand[:-1] + integrand[1:])
            Cl[il] = 4.0 * np.pi * np.sum(mid * dlnk)

        Cl = np.maximum(Cl, 0.0)
        ell_f = ell_values.astype(np.float64)
        Dl = ell_f * (ell_f + 1.0) * Cl / (2.0 * np.pi) * (T_CMB * 1e6) ** 2

        return ell_values, Cl, Dl

    # ====================================================================
    # Peak finding
    # ====================================================================

    @staticmethod
    def _find_peaks(ell_values, Dl, min_ell=100, n_peaks=7):
        """Find acoustic peaks in the D_l spectrum."""
        from scipy.signal import find_peaks
        from scipy.ndimage import gaussian_filter1d

        Dl_smooth = gaussian_filter1d(Dl, sigma=8)
        mask = ell_values > min_ell
        if not np.any(mask):
            return {'ell': [], 'Dl': []}

        offset = np.argmax(mask)
        pks, props = find_peaks(
            Dl_smooth[mask], distance=60,
            prominence=max(10, 0.01 * np.max(Dl_smooth[mask]))
        )
        pks = pks + offset

        peak_ells = [int(ell_values[p]) for p in pks[:n_peaks]]
        peak_Dls = [float(Dl[p]) for p in pks[:n_peaks]]

        return {'ell': peak_ells, 'Dl': peak_Dls}


# ============================================================================
# Convenience: run + CLASS comparison
# ============================================================================

def benchmark(h=0.6736, omega_b=0.02237, omega_cdm=0.12, plot=True):
    """
    Run the production solver and compare with CLASS if available.

    Returns the ProductionResult.
    """
    import os

    print("=" * 70)
    print("PRODUCTION SOLVER BENCHMARK")
    print("=" * 70)

    solver = ProductionSolver(h=h, omega_b=omega_b, omega_cdm=omega_cdm)
    result = solver.run()

    # Print timing
    print(f"\n--- Timing ---")
    for key, val in result.timing.items():
        print(f"  {key:12s}: {val:.3f}s")

    # Print peaks
    print(f"\n--- Peaks ---")
    for i, (ell_p, dl_p) in enumerate(zip(result.peaks['ell'],
                                            result.peaks['Dl'])):
        print(f"  Peak {i+1}: l={ell_p}, D_l={dl_p:.0f} muK^2")

    # CLASS comparison
    class_path = os.path.expanduser(
        '/Users/akaihuangm1/Desktop/github/gdm_class_public/output/'
        'mlx_ref_lcdm_cl.dat'
    )
    if os.path.exists(class_path):
        data = np.loadtxt(class_path)
        class_ell = data[:, 0].astype(int)
        class_Dl = data[:, 1] * (T_CMB * 1e6) ** 2

        # Interpolate CLASS to our ell grid
        from scipy.interpolate import interp1d
        class_interp = interp1d(class_ell, class_Dl, kind='linear',
                                fill_value='extrapolate')
        class_at_ell = class_interp(result.ell)

        # RMS relative error (where CLASS > 100 muK^2)
        mask = class_at_ell > 100
        if np.any(mask):
            rel_err = np.abs(result.Dl_TT[mask] - class_at_ell[mask]) / class_at_ell[mask]
            rms_pct = 100.0 * np.sqrt(np.mean(rel_err ** 2))
            max_pct = 100.0 * np.max(rel_err)
            print(f"\n--- CLASS Comparison ---")
            print(f"  RMS relative error: {rms_pct:.1f}%")
            print(f"  Max relative error: {max_pct:.1f}%")
            print(f"  (where CLASS D_l > 100 muK^2)")

        # Peak comparison
        from scipy.signal import find_peaks
        from scipy.ndimage import gaussian_filter1d
        class_smooth = gaussian_filter1d(class_Dl, sigma=15)
        mask_c = class_ell > 100
        if np.any(mask_c):
            offset_c = np.argmax(mask_c)
            pks_c, _ = find_peaks(class_smooth[mask_c], distance=150,
                                  prominence=50)
            pks_c = pks_c + offset_c
            class_peaks = class_ell[pks_c][:5]
            print(f"\n  Peak positions:")
            print(f"    {'Peak':>6s}  {'Production':>10s}  {'CLASS':>6s}  {'Delta':>6s}")
            for i in range(min(5, len(result.peaks['ell']),
                               len(class_peaks))):
                delta = result.peaks['ell'][i] - int(class_peaks[i])
                print(f"    {i+1:>6d}  {result.peaks['ell'][i]:>10d}  "
                      f"{int(class_peaks[i]):>6d}  {delta:>+6d}")

        if plot:
            try:
                import matplotlib
                matplotlib.use('Agg')
                import matplotlib.pyplot as plt

                fig, axes = plt.subplots(2, 1, figsize=(14, 10),
                                          gridspec_kw={'height_ratios': [3, 1]},
                                          sharex=True)
                fig.subplots_adjust(hspace=0.06)

                ax = axes[0]
                ax.plot(result.ell, result.Dl_TT, 'r-', lw=1.5,
                        label='Production solver')
                ax.plot(class_ell, class_Dl, 'k--', lw=1.0, alpha=0.6,
                        label='CLASS')
                ax.set_ylabel(r'$\mathcal{D}_\ell^{TT}$ [$\mu$K$^2$]',
                              fontsize=14)
                ax.set_title('Production Solver vs CLASS', fontsize=13)
                ax.legend(fontsize=10)
                ax.set_xlim(2, 2500)
                ax.grid(alpha=0.15)

                ax2 = axes[1]
                ratio = np.where(class_at_ell > 100,
                                 result.Dl_TT / class_at_ell, np.nan)
                ax2.plot(result.ell, ratio, 'r-', lw=1.0)
                ax2.axhline(1.0, color='k', ls='--', lw=0.8)
                ax2.set_ylabel('Production / CLASS', fontsize=11)
                ax2.set_xlabel(r'Multipole $\ell$', fontsize=14)
                ax2.set_ylim(0.5, 1.5)
                ax2.grid(alpha=0.15)

                out_dir = os.path.dirname(os.path.abspath(__file__))
                plot_path = os.path.join(out_dir,
                                         'production_vs_class.png')
                plt.savefig(plot_path, dpi=150, bbox_inches='tight')
                print(f"\n  Plot saved: {plot_path}")
                plt.close()
            except Exception as e:
                print(f"\n  Plot failed: {e}")
    else:
        print(f"\n[Production] CLASS reference not found at {class_path}")

    print("=" * 70)
    return result


if __name__ == '__main__':
    benchmark()
