"""
solver_greens.py -- Spectral Green's function for instant C_l computation.

The Boltzmann system y'(tau) = A(tau) * y(tau) is LINEAR in y. This linearity
means the transfer function Delta_l(k) is a LINEAR functional of the initial
conditions y0(k). Two levels of precomputation exploit this:

  Level 1 -- Source function cache:
      Solve the Boltzmann equations ONCE for normalized adiabatic ICs (C=1),
      caching the LOS source functions {Theta0+Psi, v_b, Phi'+Psi', Theta2}
      at each (k, tau_snap). The transfer function Delta_l(k) and C_l follow
      from this cache. Since the ICs are linear in the primordial amplitude,
      scanning A_s and n_s requires only rescaling the cached sources -- no
      re-solving of any ODEs.
      Cost: N_k ODE solves (one-time). C_l for any (A_s, n_s): ~50 ms.

  Level 2 -- Full propagator G(k, tau_s):
      For each k and snapshot tau_s, the propagator G maps ANY y0 -> y(tau_s).
      Precompute by solving n_var IVPs with unit-vector initial conditions.
      This enables C_l evaluation for arbitrary (non-adiabatic) initial
      conditions, including isocurvature modes.
      Cost: n_var * N_k IVPs (one-time). C_l for any y0: ~100 ms.

Both levels are implemented. Level 1 is sufficient for primordial parameter
scans (MCMC on A_s, n_s). Level 2 is needed for isocurvature studies.

Pipeline:
  1. Background (numpy, ~10 ms)
  2. Precompute source cache or full G (multiprocessing, ~30-120 s)
  3. LOS integration (Bessel functions, ~2 s)
  4. C_l evaluation for any primordial spectrum (~50 ms)

Usage:
  python -m mlx_class.solver_greens                  # Level 1 + validate
  python -m mlx_class.solver_greens --save cache.npz # save for reuse
  python -m mlx_class.solver_greens --load cache.npz # load + instant C_l
  python -m mlx_class.solver_greens --scan            # primordial scan
  python -m mlx_class.solver_greens --full-G          # Level 2 (slow)

Author: Sheng-Kai Huang, 2026
"""
import os
os.environ['MPLCONFIGDIR'] = '/tmp/mpl_config'

import numpy as np
import time
import sys
from multiprocessing import Pool
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

from .background import (
    Background, A_s, n_s, k_pivot, T_CMB,
    Omega_r as _OMEGA_R,
    Omega_b as _OMEGA_B,
    Omega_c as _OMEGA_C,
    Omega_m as _OMEGA_M,
    H0_Mpc as _H0_MPC,
)
from .perturbations_neutrino import (
    _f_nu, _OMEGA_GAMMA, _OMEGA_NU,
)
from .perturbations_sync import (
    make_sync_rhs, adiabatic_ic_sync, gauge_transform, diagnose_metric,
    n_var_sync, idx_fn_start, idx_fg,
    IDX_ETA, IDX_DELTA_C, IDX_DELTA_B, IDX_THETA_B, IDX_FG_START,
    L_GAMMA_MAX, L_NU_MAX_SYNC,
)
from .solver_sync import build_snapshot_grid


# ============================================================================
# Worker: solve one k-mode and return source functions at all snapshots
# ============================================================================

def _solve_source_worker(args):
    """
    Worker for Level 1: solve adiabatic IVP for one k-mode.

    Returns source functions at all snapshot times:
      [Phi_N, Psi_N, Theta0_N, vb_N, Theta2] at each snapshot.

    Identical to solver_sync._solve_single_k_worker but returned for
    separate archival and reuse.
    """
    import numpy as np
    from scipy.integrate import solve_ivp

    (k, bg_arrays, tau_end, lg_max, ln_max, method, rtol, atol,
     tau_all, a_snap, calH_snap) = args

    class _BGLite:
        __slots__ = ('tau_grid', 'a_grid', 'calH_grid', 'R_grid',
                     'kappa_dot_grid', 'tau_rec', 'tau_0')
        def __init__(self, arrays):
            self.tau_grid = arrays['tau_grid']
            self.a_grid = arrays['a_grid']
            self.calH_grid = arrays['calH_grid']
            self.R_grid = arrays['R_grid']
            self.kappa_dot_grid = arrays['kappa_dot_grid']
            self.tau_rec = arrays['tau_rec']
            self.tau_0 = arrays['tau_0']

    bg_lite = _BGLite(bg_arrays)

    from mlx_class.perturbations_sync import (
        make_sync_rhs, adiabatic_ic_sync, gauge_transform, diagnose_metric,
        IDX_ETA, IDX_THETA_B, IDX_FG_START,
    )

    tau_init = bg_lite.tau_grid[1]
    rhs_fn, nvar = make_sync_rhs(k, bg_lite, lg_max, ln_max)
    y0 = adiabatic_ic_sync(k, tau_init, lg_max, ln_max)

    sol = solve_ivp(rhs_fn, [tau_init, tau_end], y0,
                    method=method, dense_output=True,
                    rtol=rtol, atol=atol)

    if not sol.success:
        return None

    N_snap = len(tau_all)
    results = np.zeros((N_snap, 5), dtype=np.float64)

    for it in range(N_snap):
        tau = tau_all[it]
        y = sol.sol(tau)
        a = a_snap[it]
        calH = calH_snap[it]

        h_prime, eta_prime = diagnose_metric(y, k, calH, a, lg_max, ln_max)
        gt = gauge_transform(y, k, calH, a, lg_max, ln_max, h_prime, eta_prime)

        eta_val = y[IDX_ETA]
        results[it, 0] = gt['Phi_N']
        results[it, 1] = gt['Psi_N']
        results[it, 2] = y[IDX_FG_START] / 4.0 - eta_val + gt['Phi_N']
        results[it, 3] = y[IDX_THETA_B] / k + (h_prime + 6 * eta_prime) / (2 * k)
        results[it, 4] = y[IDX_FG_START + 2] / 4.0

    return results


# ============================================================================
# Worker: solve one column of full propagator G
# ============================================================================

def _precompute_k_worker(args):
    """
    Level 2 worker: compute ALL columns of G for a single k-mode.

    Solves n_var IVPs with unit vector initial conditions and evaluates
    at all snapshot times.

    Returns G_k: (N_snap, n_var, n_var) or None on failure.
    """
    import numpy as np
    from scipy.integrate import solve_ivp

    (k, n_var, bg_arrays, tau_end, lg_max, ln_max,
     method, rtol, atol, tau_snaps) = args

    class _BGLite:
        __slots__ = ('tau_grid', 'a_grid', 'calH_grid', 'R_grid',
                     'kappa_dot_grid', 'tau_rec', 'tau_0')
        def __init__(self, arrays):
            self.tau_grid = arrays['tau_grid']
            self.a_grid = arrays['a_grid']
            self.calH_grid = arrays['calH_grid']
            self.R_grid = arrays['R_grid']
            self.kappa_dot_grid = arrays['kappa_dot_grid']
            self.tau_rec = arrays['tau_rec']
            self.tau_0 = arrays['tau_0']

    bg_lite = _BGLite(bg_arrays)

    from mlx_class.perturbations_sync import make_sync_rhs

    tau_init = bg_lite.tau_grid[1]
    rhs_fn, _ = make_sync_rhs(k, bg_lite, lg_max, ln_max)

    N_snap = len(tau_snaps)
    G_k = np.zeros((N_snap, n_var, n_var), dtype=np.float64)

    for j in range(n_var):
        y0 = np.zeros(n_var, dtype=np.float64)
        y0[j] = 1.0

        sol = solve_ivp(rhs_fn, [tau_init, tau_end], y0,
                        method=method, dense_output=True,
                        rtol=rtol, atol=atol)

        if not sol.success:
            return None

        for it in range(N_snap):
            G_k[it, :, j] = sol.sol(tau_snaps[it])

    return G_k


# ============================================================================
# Green's Function Solver
# ============================================================================

class GreensFunctionSolver:
    """
    Precompute Boltzmann transfer functions for instant C_l evaluation.

    Two levels:
      Level 1: Source function cache (adiabatic mode only, fast precomputation).
      Level 2: Full propagator G(k, tau) (any initial conditions, slow precompute).

    Parameters
    ----------
    k_arr : ndarray
        Wavenumber grid in Mpc^-1.
    bg : Background
        Solved background cosmology object.
    lg_max : int
        Photon hierarchy truncation.
    ln_max : int
        Neutrino hierarchy truncation.
    """

    def __init__(self, k_arr, bg, lg_max=L_GAMMA_MAX, ln_max=L_NU_MAX_SYNC):
        self.k_arr = np.asarray(k_arr)
        self.bg = bg
        self.lg_max = lg_max
        self.ln_max = ln_max
        self.n_var = n_var_sync(lg_max, ln_max)
        self.N_k = len(k_arr)

        # Level 1: source function cache
        self.sources = None     # dict of (N_k, N_snap) arrays
        self.tau_snaps = None
        self.a_snap = None
        self.calH_snap = None

        # Level 1: transfer functions (after LOS)
        self.Delta_l_TT = None  # (N_ell, N_k)
        self.Delta_l_EE = None
        self.ell_values = None

        # Level 2: full propagator
        self.G = None           # (N_k, N_snap, n_var, n_var)

        # Level 2: transfer kernel
        self.T_TT = None        # (N_ell, N_k, n_var)
        self.T_EE = None

    # ================================================================
    # Level 1: Precompute source functions (adiabatic mode)
    # ================================================================

    def precompute_sources(self, method='Radau', rtol=1e-6, atol=1e-9,
                           n_workers=None, verbose=True):
        """
        Solve the Boltzmann equations for adiabatic ICs (C=1) at each k,
        caching the Newtonian gauge source functions at snapshot times.

        This is the same computational cost as the direct solver in
        solver_sync.py, but stores results for reuse.

        Cost: ~N_k * 0.5s / n_workers. With N_k=500, n_workers=10: ~25s.
        """
        if n_workers is None:
            n_workers = os.cpu_count() or 4

        t0 = time.time()

        # Build snapshot grid
        tau_vis, tau_late, tau_all = build_snapshot_grid(
            self.bg, N_vis=60, N_early_isw=30, N_late_isw=20)
        self.tau_snaps = tau_all
        N_snap = len(tau_all)

        self.a_snap = np.asarray(self.bg.a_at_tau(tau_all), dtype=np.float64)
        self.calH_snap = np.asarray(
            self.bg.calH_at_tau(tau_all), dtype=np.float64)

        if verbose:
            print(f"\n[GreensSolver] Level 1: Precomputing source functions")
            print(f"  N_k = {self.N_k}, N_snap = {N_snap}, "
                  f"n_var = {self.n_var}")
            print(f"  Workers: {n_workers}")
            sys.stdout.flush()

        # Serialize background
        bg_arrays = {
            'tau_grid': self.bg.tau_grid.copy(),
            'a_grid': self.bg.a_grid.copy(),
            'calH_grid': self.bg.calH_grid.copy(),
            'R_grid': self.bg.R_grid.copy(),
            'kappa_dot_grid': self.bg.kappa_dot_grid.copy(),
            'tau_rec': float(self.bg.tau_rec),
            'tau_0': float(self.bg.tau_0),
        }
        tau_end = float(tau_all[-1] + 1.0)

        work_items = [
            (k, bg_arrays, tau_end, self.lg_max, self.ln_max, method,
             rtol, atol, tau_all, self.a_snap, self.calH_snap)
            for k in self.k_arr
        ]

        if verbose:
            print(f"  Dispatching {self.N_k} k-modes...")
            sys.stdout.flush()

        with Pool(processes=n_workers) as pool:
            results = pool.map(_solve_source_worker, work_items)

        # Collect
        Phi_N = np.zeros((self.N_k, N_snap))
        Psi_N = np.zeros((self.N_k, N_snap))
        Theta0_N = np.zeros((self.N_k, N_snap))
        vb_N = np.zeros((self.N_k, N_snap))
        Theta2 = np.zeros((self.N_k, N_snap))

        n_failed = 0
        for ik, res in enumerate(results):
            if res is None:
                n_failed += 1
                if verbose and n_failed <= 3:
                    print(f"  WARNING: k={self.k_arr[ik]:.4e} failed")
                continue
            Phi_N[ik] = res[:, 0]
            Psi_N[ik] = res[:, 1]
            Theta0_N[ik] = res[:, 2]
            vb_N[ik] = res[:, 3]
            Theta2[ik] = res[:, 4]

        # Compute (Phi + Psi)' via cubic spline
        PhiPsi = Phi_N + Psi_N
        PhiPsi_prime = np.zeros_like(PhiPsi)
        for ik in range(self.N_k):
            cs = CubicSpline(tau_all, PhiPsi[ik])
            PhiPsi_prime[ik] = cs(tau_all, 1)

        self.sources = {
            'Phi_N': Phi_N,
            'Psi_N': Psi_N,
            'Theta0_N': Theta0_N,
            'vb_N': vb_N,
            'Theta2': Theta2,
            'PhiPsi_prime': PhiPsi_prime,
        }

        t_src = time.time() - t0
        if verbose:
            print(f"  Source precomputation: {t_src:.1f}s "
                  f"({n_failed} failed, {self.N_k - n_failed} OK)")
            mem_mb = sum(v.nbytes for v in self.sources.values()) / 1024**2
            print(f"  Source cache: {mem_mb:.1f} MB")

        return self

    # ================================================================
    # Level 1: LOS integration (builds Delta_l from cached sources)
    # ================================================================

    def build_transfer_functions(self, ell_values=None, verbose=True):
        """
        Build transfer functions Delta_l(k) from cached source functions.

        This performs the LOS integration: Bessel function projection of
        the cached sources. Same algorithm as solver_sync.py.

        After this, compute_cl() can evaluate C_l for any (A_s, n_s)
        instantaneously.
        """
        if self.sources is None:
            raise RuntimeError("Must call precompute_sources() first")

        t0 = time.time()

        if ell_values is None:
            ell_values = np.unique(np.concatenate([
                np.arange(2, 30, 1),
                np.arange(30, 100, 2),
                np.arange(100, 500, 4),
                np.arange(500, 1500, 8),
                np.arange(1500, 2501, 12),
            ])).astype(int)
        self.ell_values = ell_values
        N_ell = len(ell_values)
        N_snap = len(self.tau_snaps)
        dtau_all = np.diff(self.tau_snaps)

        if verbose:
            print(f"\n[GreensSolver] Building transfer functions")
            print(f"  N_ell = {N_ell}, N_k = {self.N_k}, "
                  f"N_snap = {N_snap}")
            sys.stdout.flush()

        # Background at snapshots
        g_snap = self.bg.visibility_at_tau(self.tau_snaps)
        kappa_snap = self.bg.kappa_at_tau(self.tau_snaps)
        exp_neg_kappa = np.exp(-kappa_snap)
        chi_snap = self.bg.tau_0 - self.tau_snaps

        # Polarization
        alpha_P = 1.7
        pol_prefactor = 0.75 * alpha_P

        eps_prefactors = np.zeros(N_ell)
        for il, ell in enumerate(ell_values):
            l = int(ell)
            if l >= 2:
                eps_prefactors[il] = np.sqrt(
                    float((l - 1) * l * (l + 1) * (l + 2)))

        # Bessel table
        try:
            from .bessel_cache import CachedBesselTable
            import mlx.core as mx
            x_max_bessel = (float(np.max(self.k_arr) * np.max(chi_snap))
                            * 1.05 + 50.0)
            bessel_table = CachedBesselTable(
                ell_values, x_max=x_max_bessel,
                N_cheb=64, seg_width=80, verbose=verbose)
            use_gpu_bessel = True
        except Exception as e:
            if verbose:
                print(f"  [WARNING] Bessel table failed ({e}), using scipy")
            use_gpu_bessel = False

        # Extract source arrays
        Phi_N = self.sources['Phi_N']
        Psi_N = self.sources['Psi_N']
        Theta0_N = self.sources['Theta0_N']
        vb_N = self.sources['vb_N']
        Theta2 = self.sources['Theta2']
        PhiPsi_prime = self.sources['PhiPsi_prime']

        # Transfer functions
        Delta_l = np.zeros((N_ell, self.N_k), dtype=np.float64)
        Delta_l_E = np.zeros((N_ell, self.N_k), dtype=np.float64)

        if verbose:
            print(f"  Computing LOS transfer functions (TT + EE)...")
            sys.stdout.flush()

        for it in range(N_snap):
            chi = chi_snap[it]
            if chi <= 0:
                continue

            if it == 0:
                w = 0.5 * dtau_all[0]
            elif it == N_snap - 1:
                w = 0.5 * dtau_all[-1]
            else:
                w = 0.5 * (dtau_all[max(0, it - 1)]
                           + dtau_all[min(it, len(dtau_all) - 1)])

            S_SW = g_snap[it] * (Theta0_N[:, it] + Psi_N[:, it])
            S_Dop = g_snap[it] * vb_N[:, it]
            S_ISW = exp_neg_kappa[it] * PhiPsi_prime[:, it]
            S_E = g_snap[it] * pol_prefactor * Theta2[:, it]

            if use_gpu_bessel:
                x_arr = self.k_arr * chi
                jl_mx, jlp_mx = bessel_table.eval_jl_jlp(
                    x_arr.astype(np.float32))
                mx.eval(jl_mx, jlp_mx)
                jl = np.array(jl_mx)
                jlp = np.array(jlp_mx)

                x_safe = np.where(np.abs(x_arr) < 1e-10, 1e-10, x_arr)
                inv_x2 = 1.0 / (x_safe ** 2)

                for il in range(N_ell):
                    Delta_l[il, :] += w * (
                        S_SW * jl[il, :] +
                        S_Dop * jlp[il, :] +
                        S_ISW * jl[il, :]
                    )
                    if eps_prefactors[il] > 0:
                        eps_l = eps_prefactors[il] * jl[il, :] * inv_x2
                        eps_l = np.where(
                            np.abs(x_arr) < 1e-10, 0.0, eps_l)
                        Delta_l_E[il, :] += w * S_E * eps_l
            else:
                from scipy.special import spherical_jn
                x = self.k_arr * chi
                x_safe = np.where(np.abs(x) < 1e-10, 1e-10, x)
                inv_x2 = 1.0 / (x_safe ** 2)
                for il, ell in enumerate(ell_values):
                    jl = spherical_jn(int(ell), x)
                    jlp = spherical_jn(int(ell), x, derivative=True)
                    Delta_l[il, :] += w * (
                        S_SW * jl + S_Dop * jlp + S_ISW * jl)
                    if eps_prefactors[il] > 0:
                        eps_l = eps_prefactors[il] * jl * inv_x2
                        eps_l = np.where(np.abs(x) < 1e-10, 0.0, eps_l)
                        Delta_l_E[il, :] += w * S_E * eps_l

        self.Delta_l_TT = Delta_l
        self.Delta_l_EE = Delta_l_E

        t_los = time.time() - t0
        if verbose:
            print(f"  LOS integration: {t_los:.2f}s")

        return self

    # ================================================================
    # Level 1: Instant C_l computation
    # ================================================================

    def compute_cl(self, A_s_val=A_s, n_s_val=n_s, verbose=True):
        """
        Compute C_l using cached transfer functions.

        Since the Boltzmann system is linear and the ICs are normalized to
        C=1, the physical transfer function is:
            Delta_l^phys(k) = Delta_l^cache(k) * scale(k)
        where scale(k) depends only on the primordial spectrum.

        For adiabatic ICs with unit curvature perturbation C=1:
            y0(k) ~ (k * tau_init)^n   (power-law in k*tau_init)
        The primordial spectrum P_R(k) = A_s * (k/k_pivot)^(n_s-1)
        enters only through the C_l integral.

        Runs in < 100 ms for any (A_s, n_s).

        Returns
        -------
        dict with keys: ell, Cl_TT, Dl_TT, Cl_EE, Dl_EE, Cl_TE, Dl_TE
        """
        if self.Delta_l_TT is None:
            raise RuntimeError(
                "Must call build_transfer_functions() first")

        t0 = time.time()

        k_arr = self.k_arr
        ell_values = self.ell_values
        N_ell = len(ell_values)

        # Primordial spectrum
        P_R = A_s_val * (k_arr / k_pivot) ** (n_s_val - 1.0)

        # C_l = 4 pi int dk/k P_R |Delta_l|^2
        lnk = np.log(k_arr)
        dlnk = np.diff(lnk)
        ell_f = ell_values.astype(float)
        uK2 = (T_CMB * 1e6) ** 2
        norm = (2.0 / 3.0) ** 2

        # TT
        integrand_TT = P_R[None, :] * self.Delta_l_TT ** 2
        mid_TT = 0.5 * (integrand_TT[:, :-1] + integrand_TT[:, 1:])
        Cl_TT = 4.0 * np.pi * np.sum(mid_TT * dlnk[None, :], axis=1)
        Cl_TT = np.maximum(Cl_TT, 0.0)
        Dl_TT = norm * ell_f * (ell_f + 1.0) * Cl_TT / (2.0 * np.pi) * uK2

        # EE
        integrand_EE = P_R[None, :] * self.Delta_l_EE ** 2
        mid_EE = 0.5 * (integrand_EE[:, :-1] + integrand_EE[:, 1:])
        Cl_EE = 4.0 * np.pi * np.sum(mid_EE * dlnk[None, :], axis=1)
        Cl_EE = np.maximum(Cl_EE, 0.0)
        Dl_EE = norm * ell_f * (ell_f + 1.0) * Cl_EE / (2.0 * np.pi) * uK2

        # TE
        integrand_TE = P_R[None, :] * self.Delta_l_TT * self.Delta_l_EE
        mid_TE = 0.5 * (integrand_TE[:, :-1] + integrand_TE[:, 1:])
        Cl_TE = 4.0 * np.pi * np.sum(mid_TE * dlnk[None, :], axis=1)
        Dl_TE = norm * ell_f * (ell_f + 1.0) * Cl_TE / (2.0 * np.pi) * uK2

        t_cl = time.time() - t0

        if verbose:
            print(f"[GreensSolver] C_l computed in {t_cl*1000:.1f} ms")

        return {
            'ell': ell_values,
            'Cl_TT': Cl_TT,
            'Dl_TT': Dl_TT,
            'Cl_EE': Cl_EE,
            'Dl_EE': Dl_EE,
            'Cl_TE': Cl_TE,
            'Dl_TE': Dl_TE,
            'Delta_l_TT': self.Delta_l_TT,
            'Delta_l_EE': self.Delta_l_EE,
            'timing': {'cl_ms': t_cl * 1000},
        }

    # ================================================================
    # Level 1: Scan primordial parameters
    # ================================================================

    def scan_primordial(self, A_s_arr, n_s_arr, verbose=True):
        """
        Compute C_l for a grid of (A_s, n_s) values.

        Since the transfer functions are independent of the primordial
        spectrum, this is essentially free after the LOS integration.

        Parameters
        ----------
        A_s_arr : array-like of float
        n_s_arr : array-like of float

        Returns
        -------
        results : list of dicts
        """
        A_s_arr = np.asarray(A_s_arr)
        n_s_arr = np.asarray(n_s_arr)
        results = []
        t0 = time.time()

        for A_s_val in A_s_arr:
            for n_s_val in n_s_arr:
                res = self.compute_cl(A_s_val, n_s_val, verbose=False)
                results.append({
                    'A_s': A_s_val,
                    'n_s': n_s_val,
                    **res,
                })

        t_total = time.time() - t0
        n_total = len(A_s_arr) * len(n_s_arr)
        if verbose:
            print(f"\n[GreensSolver] Scanned {n_total} cosmologies "
                  f"in {t_total*1000:.1f} ms "
                  f"({t_total/n_total*1000:.1f} ms/cosmology)")
        return results

    # ================================================================
    # Level 2: Full propagator precomputation
    # ================================================================

    def precompute_full_G(self, method='Radau', rtol=1e-6, atol=1e-9,
                          n_workers=None, verbose=True):
        """
        Precompute the FULL propagator G(k, tau_s) for all k and snapshots.

        This solves n_var IVPs per k-mode with unit-vector initial conditions.
        Much more expensive than Level 1, but allows C_l evaluation for
        arbitrary (non-adiabatic) initial conditions.

        Cost: ~n_var * N_k * 0.5s / n_workers. Very slow.
        """
        if n_workers is None:
            n_workers = os.cpu_count() or 4

        t0 = time.time()

        # Reuse snapshot grid from Level 1 if available
        if self.tau_snaps is None:
            tau_vis, tau_late, tau_all = build_snapshot_grid(
                self.bg, N_vis=60, N_early_isw=30, N_late_isw=20)
            self.tau_snaps = tau_all
            self.a_snap = np.asarray(
                self.bg.a_at_tau(tau_all), dtype=np.float64)
            self.calH_snap = np.asarray(
                self.bg.calH_at_tau(tau_all), dtype=np.float64)

        N_snap = len(self.tau_snaps)

        if verbose:
            print(f"\n[GreensSolver] Level 2: Precomputing full G(k, tau_s)")
            print(f"  N_k = {self.N_k}, n_var = {self.n_var}, "
                  f"N_snap = {N_snap}")
            print(f"  Total IVPs: {self.N_k * self.n_var}")
            print(f"  Workers: {n_workers}")
            sys.stdout.flush()

        bg_arrays = {
            'tau_grid': self.bg.tau_grid.copy(),
            'a_grid': self.bg.a_grid.copy(),
            'calH_grid': self.bg.calH_grid.copy(),
            'R_grid': self.bg.R_grid.copy(),
            'kappa_dot_grid': self.bg.kappa_dot_grid.copy(),
            'tau_rec': float(self.bg.tau_rec),
            'tau_0': float(self.bg.tau_0),
        }
        tau_end = float(self.tau_snaps[-1] + 1.0)

        work_items = [
            (k, self.n_var, bg_arrays, tau_end, self.lg_max, self.ln_max,
             method, rtol, atol, self.tau_snaps)
            for k in self.k_arr
        ]

        if verbose:
            print(f"  Dispatching {self.N_k} k-modes...")
            sys.stdout.flush()

        with Pool(processes=n_workers) as pool:
            results = pool.map(_precompute_k_worker, work_items)

        self.G = np.zeros(
            (self.N_k, N_snap, self.n_var, self.n_var), dtype=np.float64)
        n_failed = 0
        for ik, res in enumerate(results):
            if res is None:
                n_failed += 1
            else:
                self.G[ik] = res

        t_G = time.time() - t0
        if verbose:
            mem_mb = self.G.nbytes / 1024**2
            print(f"  Full G: {t_G:.1f}s ({n_failed} failed)")
            print(f"  G memory: {mem_mb:.1f} MB")

        return self

    # ================================================================
    # Level 2: Propagate arbitrary initial conditions
    # ================================================================

    def propagate(self, y0_all):
        """
        Apply the propagator to get the state at all snapshot times.

        Parameters
        ----------
        y0_all : ndarray of shape (N_k, n_var)

        Returns
        -------
        y_snap : ndarray of shape (N_k, N_snap, n_var)
        """
        if self.G is None:
            raise RuntimeError("Must call precompute_full_G() first")
        return np.einsum('ksnm,km->ksn', self.G, y0_all)

    def compute_cl_from_y0(self, y0_all, A_s_val=A_s, n_s_val=n_s,
                           verbose=True):
        """
        Compute C_l from arbitrary initial conditions using the full
        propagator G.

        Parameters
        ----------
        y0_all : ndarray of shape (N_k, n_var)
            Custom initial conditions for each k-mode.
        A_s_val, n_s_val : float
            Primordial spectrum parameters.

        Returns
        -------
        dict with Cl_TT, Dl_TT, etc.
        """
        if self.G is None:
            raise RuntimeError("Must call precompute_full_G() first")

        t0 = time.time()

        # Propagate
        y_snap = self.propagate(y0_all)

        # Extract source functions
        N_snap = len(self.tau_snaps)
        Phi_N = np.zeros((self.N_k, N_snap))
        Psi_N = np.zeros((self.N_k, N_snap))
        Theta0_N = np.zeros((self.N_k, N_snap))
        vb_N = np.zeros((self.N_k, N_snap))
        Theta2 = np.zeros((self.N_k, N_snap))

        for ik in range(self.N_k):
            k = self.k_arr[ik]
            for it in range(N_snap):
                y = y_snap[ik, it, :]
                a = self.a_snap[it]
                calH = self.calH_snap[it]

                h_prime, eta_prime = diagnose_metric(
                    y, k, calH, a, self.lg_max, self.ln_max)
                gt = gauge_transform(
                    y, k, calH, a, self.lg_max, self.ln_max,
                    h_prime, eta_prime)

                eta_val = y[IDX_ETA]
                Phi_N[ik, it] = gt['Phi_N']
                Psi_N[ik, it] = gt['Psi_N']
                Theta0_N[ik, it] = (y[IDX_FG_START] / 4.0 - eta_val
                                    + gt['Phi_N'])
                vb_N[ik, it] = (y[IDX_THETA_B] / k
                                + (h_prime + 6 * eta_prime) / (2 * k))
                Theta2[ik, it] = y[IDX_FG_START + 2] / 4.0

        # Phi+Psi derivative
        PhiPsi = Phi_N + Psi_N
        PhiPsi_prime = np.zeros_like(PhiPsi)
        for ik in range(self.N_k):
            cs = CubicSpline(self.tau_snaps, PhiPsi[ik])
            PhiPsi_prime[ik] = cs(self.tau_snaps, 1)

        # Store temporarily for LOS
        old_sources = self.sources
        old_Delta_TT = self.Delta_l_TT
        old_Delta_EE = self.Delta_l_EE

        self.sources = {
            'Phi_N': Phi_N, 'Psi_N': Psi_N,
            'Theta0_N': Theta0_N, 'vb_N': vb_N,
            'Theta2': Theta2, 'PhiPsi_prime': PhiPsi_prime,
        }

        # LOS integration
        self.build_transfer_functions(verbose=False)

        # Compute C_l
        result = self.compute_cl(A_s_val, n_s_val, verbose=verbose)

        # Restore original cached data
        self.sources = old_sources
        self.Delta_l_TT = old_Delta_TT
        self.Delta_l_EE = old_Delta_EE

        t_total = time.time() - t0
        if verbose:
            print(f"  Total (propagate + sources + LOS + C_l): "
                  f"{t_total*1000:.0f} ms")

        return result

    # ================================================================
    # Save / Load
    # ================================================================

    def save(self, path):
        """
        Save source cache and transfer functions to disk.

        For Level 1: saves sources + Delta_l + metadata.
        For Level 2: additionally saves G (warning: can be very large).
        """
        data = {
            'k_arr': self.k_arr,
            'tau_snaps': self.tau_snaps,
            'a_snap': self.a_snap,
            'calH_snap': self.calH_snap,
            'lg_max': np.array(self.lg_max),
            'ln_max': np.array(self.ln_max),
        }

        # Level 1 data
        if self.sources is not None:
            for key, val in self.sources.items():
                data[f'src_{key}'] = val
        if self.Delta_l_TT is not None:
            data['Delta_l_TT'] = self.Delta_l_TT
            data['Delta_l_EE'] = self.Delta_l_EE
            data['ell_values'] = self.ell_values

        # Level 2 data (optional, large)
        if self.G is not None:
            data['G'] = self.G

        np.savez_compressed(path, **data)
        size_mb = os.path.getsize(path) / 1024**2
        print(f"[GreensSolver] Saved to {path} ({size_mb:.1f} MB)")

    def load(self, path, bg=None):
        """Load precomputed data from disk."""
        data = np.load(path)

        self.k_arr = data['k_arr']
        self.N_k = len(self.k_arr)
        self.tau_snaps = data['tau_snaps']
        self.a_snap = data['a_snap']
        self.calH_snap = data['calH_snap']
        self.lg_max = int(data['lg_max'])
        self.ln_max = int(data['ln_max'])
        self.n_var = n_var_sync(self.lg_max, self.ln_max)

        if bg is not None:
            self.bg = bg

        # Level 1
        src_keys = ['Phi_N', 'Psi_N', 'Theta0_N', 'vb_N', 'Theta2',
                     'PhiPsi_prime']
        if f'src_{src_keys[0]}' in data:
            self.sources = {}
            for key in src_keys:
                self.sources[key] = data[f'src_{key}']

        if 'Delta_l_TT' in data:
            self.Delta_l_TT = data['Delta_l_TT']
            self.Delta_l_EE = data['Delta_l_EE']
            self.ell_values = data['ell_values']

        # Level 2
        if 'G' in data:
            self.G = data['G']

        has_src = self.sources is not None
        has_dl = self.Delta_l_TT is not None
        has_G = self.G is not None

        print(f"[GreensSolver] Loaded from {path}")
        print(f"  N_k={self.N_k}, n_var={self.n_var}, "
              f"N_snap={len(self.tau_snaps)}")
        print(f"  Sources: {'yes' if has_src else 'no'}, "
              f"Delta_l: {'yes' if has_dl else 'no'}, "
              f"G: {'yes' if has_G else 'no'}")

        return self


# ============================================================================
# Validation helper
# ============================================================================

def validate_against_direct(gs_result, direct_result, verbose=True):
    """
    Compare C_l from Green's function solver with direct solver.

    Parameters
    ----------
    gs_result : dict from compute_cl()
    direct_result : dict from run_sync_solver_parallel()

    Returns
    -------
    metrics : dict
    """
    ell_gs = gs_result['ell']
    ell_direct = direct_result['ell']

    common = np.intersect1d(ell_gs, ell_direct)
    idx_gs = np.searchsorted(ell_gs, common)
    idx_direct = np.searchsorted(ell_direct, common)

    Dl_TT_gs = gs_result['Dl_TT'][idx_gs]
    Dl_TT_direct = direct_result['Dl_TT'][idx_direct]

    mask = Dl_TT_direct > 1.0
    if np.any(mask):
        rel_diff = np.abs(
            Dl_TT_gs[mask] - Dl_TT_direct[mask]) / Dl_TT_direct[mask]
        max_rel = np.max(rel_diff)
        mean_rel = np.mean(rel_diff)
        median_rel = np.median(rel_diff)
    else:
        max_rel = mean_rel = median_rel = 0.0

    metrics = {
        'max_rel_diff': max_rel,
        'mean_rel_diff': mean_rel,
        'median_rel_diff': median_rel,
        'n_common_ell': len(common),
    }

    if verbose:
        print(f"\n[Validation] Green's function vs direct solver")
        print(f"  Common multipoles: {len(common)}")
        print(f"  Max relative diff (TT): {max_rel:.6f} "
              f"({max_rel*100:.4f}%)")
        print(f"  Mean relative diff: {mean_rel:.6f} "
              f"({mean_rel*100:.4f}%)")
        print(f"  Median relative diff: {median_rel:.6f} "
              f"({median_rel*100:.4f}%)")

        # Also check EE
        Dl_EE_gs = gs_result['Dl_EE'][idx_gs]
        Dl_EE_direct = direct_result['Dl_EE'][idx_direct]
        mask_ee = Dl_EE_direct > 0.01
        if np.any(mask_ee):
            rel_ee = np.abs(
                Dl_EE_gs[mask_ee] - Dl_EE_direct[mask_ee]
            ) / Dl_EE_direct[mask_ee]
            print(f"  Max relative diff (EE): {np.max(rel_ee):.6f} "
                  f"({np.max(rel_ee)*100:.4f}%)")

    return metrics


# ============================================================================
# CLI entry point
# ============================================================================

def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Spectral Green's function Boltzmann solver")
    parser.add_argument('--N_k', type=int, default=500)
    parser.add_argument('--k_min', type=float, default=3e-4)
    parser.add_argument('--k_max', type=float, default=0.35)
    parser.add_argument('--method', type=str, default='Radau')
    parser.add_argument('--rtol', type=float, default=1e-6)
    parser.add_argument('--atol', type=float, default=1e-9)
    parser.add_argument('--lg_max', type=int, default=L_GAMMA_MAX)
    parser.add_argument('--ln_max', type=int, default=L_NU_MAX_SYNC)
    parser.add_argument('--n_workers', type=int, default=None)
    parser.add_argument('--save', type=str, default=None,
                        help='Save cache to file (e.g. greens_cache.npz)')
    parser.add_argument('--load', type=str, default=None,
                        help='Load cache from file')
    parser.add_argument('--validate', action='store_true',
                        help='Run direct solver for comparison')
    parser.add_argument('--scan', action='store_true',
                        help='Scan primordial parameters')
    parser.add_argument('--full-G', action='store_true',
                        help='Precompute full propagator (Level 2)')
    args = parser.parse_args()

    t_total_start = time.time()

    # ================================================================
    # Step 1: Background
    # ================================================================
    print("=" * 65)
    print("Spectral Green's Function Boltzmann Solver")
    print("=" * 65)
    print(f"\n--- Step 1: Background ---")

    t0 = time.time()
    bg = Background(khronon=False, recombination='recfast')
    bg.solve()
    t_bg = time.time() - t0
    print(f"Background: {t_bg:.2f}s")

    # ================================================================
    # Step 2: Create or load solver
    # ================================================================
    k_arr = np.geomspace(args.k_min, args.k_max, args.N_k)
    gs = GreensFunctionSolver(k_arr, bg, args.lg_max, args.ln_max)

    if args.load:
        print(f"\n--- Loading from {args.load} ---")
        gs.load(args.load, bg=bg)
    else:
        print(f"\n--- Step 2: Precompute source functions ---")
        gs.precompute_sources(
            method=args.method, rtol=args.rtol, atol=args.atol,
            n_workers=args.n_workers, verbose=True)

        print(f"\n--- Step 3: LOS integration ---")
        gs.build_transfer_functions(verbose=True)

        if args.full_G:
            print(f"\n--- Step 3b: Full propagator G ---")
            gs.precompute_full_G(
                method=args.method, rtol=args.rtol, atol=args.atol,
                n_workers=args.n_workers, verbose=True)

        if args.save:
            gs.save(args.save)

    # ================================================================
    # Step 3: Compute C_l (instant!)
    # ================================================================
    print(f"\n--- Step 4: Compute C_l ---")
    result = gs.compute_cl(verbose=True)

    ell = result['ell']
    Dl_TT = result['Dl_TT']
    Dl_EE = result['Dl_EE']
    Dl_TE = result['Dl_TE']

    print(f"\n--- D_l values at key multipoles ---")
    print(f"{'l':>6s}  {'D_l^TT':>12s}  {'D_l^EE':>12s}  {'D_l^TE':>12s}")
    print("-" * 50)
    for target_ell in [2, 10, 50, 100, 220, 500, 800, 1000, 1500, 2000]:
        idx = np.argmin(np.abs(ell - target_ell))
        print(f"  {ell[idx]:5d}  {Dl_TT[idx]:12.2f}  "
              f"{Dl_EE[idx]:12.4f}  {Dl_TE[idx]:12.2f}")

    # ================================================================
    # Step 4 (optional): Validate against direct solver
    # ================================================================
    if args.validate:
        print(f"\n--- Validation: direct solver comparison ---")
        from .solver_sync import run_sync_solver_parallel
        direct = run_sync_solver_parallel(
            N_k=args.N_k, k_min=args.k_min, k_max=args.k_max,
            method=args.method, lg_max=args.lg_max, ln_max=args.ln_max,
            rtol=args.rtol, atol=args.atol,
            n_workers=args.n_workers, verbose=True)
        validate_against_direct(result, direct, verbose=True)

    # ================================================================
    # Step 5 (optional): Primordial parameter scan
    # ================================================================
    if args.scan:
        print(f"\n--- Primordial parameter scan ---")
        A_s_arr = np.array([1.8e-9, 2.0e-9, 2.1e-9, 2.2e-9, 2.4e-9])
        n_s_arr = np.array([0.94, 0.96, 0.9649, 0.97, 0.99])
        scan_results = gs.scan_primordial(A_s_arr, n_s_arr)

        print(f"\n{'A_s':>10s}  {'n_s':>6s}  {'D_l(220)':>10s}  "
              f"{'D_l(550)':>10s}  {'D_l(1000)':>10s}")
        print("-" * 55)
        for res in scan_results:
            idx_220 = np.argmin(np.abs(res['ell'] - 220))
            idx_550 = np.argmin(np.abs(res['ell'] - 550))
            idx_1000 = np.argmin(np.abs(res['ell'] - 1000))
            print(f"  {res['A_s']:.2e}  {res['n_s']:.4f}  "
                  f"{res['Dl_TT'][idx_220]:10.2f}  "
                  f"{res['Dl_TT'][idx_550]:10.2f}  "
                  f"{res['Dl_TT'][idx_1000]:10.2f}")

    # ================================================================
    # Step 6 (optional): Demonstrate instant re-evaluation
    # ================================================================
    print(f"\n--- Timing: 100 C_l evaluations ---")
    t0 = time.time()
    for _ in range(100):
        _ = gs.compute_cl(verbose=False)
    t_100 = time.time() - t0
    print(f"  100 evaluations in {t_100*1000:.1f} ms "
          f"({t_100/100*1000:.2f} ms each)")

    # ================================================================
    # Summary
    # ================================================================
    t_total = time.time() - t_total_start
    print(f"\n{'='*65}")
    print(f"Total time: {t_total:.1f}s")
    print(f"C_l evaluation: {result['timing']['cl_ms']:.1f} ms")
    print(f"100x re-evaluation: {t_100/100*1000:.2f} ms each")
    if args.scan:
        n_scan = len(A_s_arr) * len(n_s_arr)
        print(f"Scan: {n_scan} cosmologies in "
              f"{sum(1 for _ in scan_results)/n_scan * t_100/100*1000:.1f} ms")
    print("=" * 65)


if __name__ == '__main__':
    main()
