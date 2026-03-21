"""
solver_scipy_bdf.py -- Scipy BDF reference solver for accuracy validation.

PURPOSE:
  Use scipy.integrate.solve_ivp with method='BDF' (or 'Radau') to solve the
  Boltzmann hierarchy for each k-mode SEQUENTIALLY. This is SLOW (~2-5 min
  for 50 k-modes) but provides a gold-standard accuracy reference because:

    - BDF is an implicit multistep method, A-stable, handles stiff Thomson
      scattering without operator splitting or IMEX tricks.
    - rtol=1e-8, atol=1e-10 gives ~8 digits of precision in the ODE.
    - No Strang splitting error, no IMEX exponential integrator approximation.

  By comparing scipy BDF output against CLASS, we isolate whether accuracy
  deficits come from the ODE solver or from the physics equations themselves.

PHYSICS:
  Same Boltzmann + Einstein equations as solver_accurate.py:
    - Conformal Newtonian gauge: Phi, Psi with anisotropic stress
    - Photon hierarchy: Theta_0 ... Theta_{l_gamma_max} with Thomson collision
    - Massless neutrino hierarchy: N_0 ... N_{l_nu_max}
    - Baryons: delta_b, v_b with Thomson drag
    - CDM: delta_c, v_c
    - Psi = Phi - 12 H0^2 / (a^2 k^2) * (Omega_nu N_2 + Omega_gamma Theta_2)

  State vector per k-mode: [Phi, delta_b, v_b, delta_c, v_c,
                             Theta_0, ..., Theta_{lgmax},
                             N_0, ..., N_{lnmax}]

USAGE:
  python -m mlx_class.solver_scipy_bdf [--N_k 50] [--method BDF]

Author: Sheng-Kai Huang, 2026
"""
import os
os.environ['MPLCONFIGDIR'] = '/tmp/mpl_config'

import numpy as np
import time
import sys

from scipy.integrate import solve_ivp

from .background import (
    Background, A_s, n_s, k_pivot, T_CMB,
    Omega_r as _OMEGA_R,
    Omega_b as _OMEGA_B,
    Omega_c as _OMEGA_C,
    H0_Mpc as _H0_MPC,
    Y_He,
)

from .perturbations_neutrino import (
    N_EFF, _f_nu, _OMEGA_GAMMA, _OMEGA_NU, L_NU_MAX,
)

from .solver_accurate import (
    IDX_PHI, IDX_DELTA_B, IDX_V_B, IDX_DELTA_C, IDX_V_C,
    IDX_THETA_START, idx_theta, idx_n_start, n_var_total,
    build_tau_grid, _compute_psi,
)


# ============================================================================
# RHS function for scipy: single k-mode, pure numpy
# ============================================================================

def _make_rhs(k, bg, Og, On, Ob, Oc, H0, lgmax, lnmax, _n_start):
    """
    Return a function f(tau, y) -> dy/dtau for a single wavenumber k.

    This closure captures all cosmological parameters and background
    interpolators. The returned function is called by scipy's ODE solver.

    The equations are IDENTICAL to solver_accurate.py's
    _deriv_streaming_full + Thomson collision terms, but written for a
    single k-mode as a 1D numpy array (not batched over k).

    PERFORMANCE: Pre-tabulates background quantities as numpy arrays and
    uses np.interp (C-speed) instead of scipy interp1d objects. This is
    critical because BDF calls the RHS ~10,000-50,000 times per k-mode.
    """
    k2 = k * k
    H02 = H0 * H0

    # Pre-tabulate background on the bg.tau_grid for fast np.interp lookup.
    # np.interp is ~100x faster than scipy interp1d for scalar evaluation.
    _tau_tab = bg.tau_grid.copy()
    _calH_tab = bg.calH_grid.copy()
    _R_tab = bg.R_grid.copy()
    _a_tab = bg.a_grid.copy()
    _kd_tab = bg.kappa_dot_grid.copy()

    def rhs(tau, y):
        """Boltzmann + Einstein RHS for a single k-mode."""
        # Background quantities via fast np.interp (scalar in, scalar out)
        calH = np.interp(tau, _tau_tab, _calH_tab)
        R = np.interp(tau, _tau_tab, _R_tab)
        a = np.interp(tau, _tau_tab, _a_tab)
        kappa_dot = np.interp(tau, _tau_tab, _kd_tab)  # negative

        inv_a2 = 1.0 / (a * a)
        tau_safe = max(tau, 1e-10)

        # Extract state variables
        Phi = y[IDX_PHI]
        delta_b = y[IDX_DELTA_B]
        v_b = y[IDX_V_B]
        delta_c = y[IDX_DELTA_C]
        v_c = y[IDX_V_C]
        Theta_0 = y[idx_theta(0)]
        Theta_1 = y[idx_theta(1)]
        Theta_2 = y[idx_theta(2)] if lgmax >= 2 else 0.0
        N_0 = y[_n_start]
        N_1 = y[_n_start + 1]
        N_2 = y[_n_start + 2] if lnmax >= 2 else 0.0

        # Anisotropic stress: Psi != Phi
        Psi = Phi - 12.0 * H02 * inv_a2 / k2 * (On * N_2 + Og * Theta_2)

        # ---- Phi equation (Poisson) ----
        S = (Og * inv_a2 * 4.0 * Theta_0
             + On * inv_a2 * 4.0 * N_0
             + Ob / a * delta_b
             + Oc / a * delta_c)
        Phi_dot = -calH * Phi - k2 * Psi / (3.0 * calH) - 0.5 * H02 / calH * S

        # ---- Photon hierarchy ----
        # Thomson scattering rate (|kappa_dot| since kappa_dot is negative)
        abs_kd = abs(kappa_dot)

        dTheta_0 = -k * Theta_1 - Phi_dot

        # Theta_1: streaming + Thomson collision
        # Collision: -|kd| * (Theta_1 - v_b/3)
        dTheta_1 = (k / 3.0 * (Theta_0 - 2.0 * Theta_2 + Psi)
                     + abs_kd * (-Theta_1 + v_b / 3.0))

        # Theta_2: streaming + Thomson polarization damping
        # Collision: -0.9 * |kd| * Theta_2
        Theta_3 = y[idx_theta(3)] if lgmax >= 3 else 0.0
        dTheta_2 = (k / 5.0 * (2.0 * Theta_1 - 3.0 * Theta_3)
                     - 0.9 * abs_kd * Theta_2)

        # Higher photon multipoles: streaming + damping
        dTheta = np.zeros(lgmax + 1)
        dTheta[0] = dTheta_0
        dTheta[1] = dTheta_1
        if lgmax >= 2:
            dTheta[2] = dTheta_2
        for l in range(3, lgmax):
            Theta_lm1 = y[idx_theta(l - 1)]
            Theta_lp1 = y[idx_theta(l + 1)]
            dTheta[l] = (k / (2.0 * l + 1.0) * (l * Theta_lm1 - (l + 1) * Theta_lp1)
                          - abs_kd * y[idx_theta(l)])
        # Truncation boundary for photon hierarchy
        if lgmax >= 1:
            Theta_lgmax_m1 = y[idx_theta(lgmax - 1)]
            Theta_lgmax = y[idx_theta(lgmax)]
            dTheta[lgmax] = (k * Theta_lgmax_m1 * lgmax / (2.0 * lgmax + 1.0)
                              - (lgmax + 1.0) / tau_safe * Theta_lgmax
                              - abs_kd * Theta_lgmax)

        # ---- Baryon equations ----
        # Thomson drag on baryons: |kd|/R * (v_b - 3*Theta_1)
        d_delta_b = -k * v_b - 3.0 * Phi_dot
        d_v_b = (-calH * v_b + k * Psi
                  + abs_kd / R * (3.0 * Theta_1 - v_b)) if R > 1e-10 else (
                  -calH * v_b + k * Psi)

        # ---- CDM equations ----
        d_delta_c = -(k * v_c + 3.0 * Phi_dot)
        d_v_c = -calH * v_c + k * Psi

        # ---- Neutrino hierarchy ----
        dN = np.zeros(lnmax + 1)
        dN[0] = -k * N_1 - Phi_dot
        dN[1] = k / 3.0 * (N_0 - 2.0 * N_2 + Psi)
        for l in range(2, lnmax):
            N_lm1 = y[_n_start + l - 1]
            N_lp1 = y[_n_start + l + 1]
            dN[l] = k / (2.0 * l + 1.0) * (l * N_lm1 - (l + 1) * N_lp1)
        # Truncation boundary for neutrino hierarchy
        N_prev = y[_n_start + lnmax - 1]
        N_last = y[_n_start + lnmax]
        dN[lnmax] = (k * N_prev * lnmax / (2.0 * lnmax + 1.0)
                      - (lnmax + 1.0) / tau_safe * N_last)

        # ---- Assemble output ----
        dydt = np.zeros_like(y)
        dydt[IDX_PHI] = Phi_dot
        dydt[IDX_DELTA_B] = d_delta_b
        dydt[IDX_V_B] = d_v_b
        dydt[IDX_DELTA_C] = d_delta_c
        dydt[IDX_V_C] = d_v_c
        dydt[IDX_THETA_START:IDX_THETA_START + lgmax + 1] = dTheta
        dydt[_n_start:_n_start + lnmax + 1] = dN

        return dydt

    return rhs


# ============================================================================
# Analytical Jacobian sparsity pattern (helps BDF/Radau converge faster)
# ============================================================================

def _make_jac_sparsity(nvar, lgmax, lnmax, _n_start):
    """
    Build a boolean sparsity pattern for the Jacobian.

    Returns a (nvar, nvar) array with 1 where J[i,j] may be nonzero.
    This helps scipy's BDF/Radau use sparse LU factorization.
    """
    from scipy.sparse import lil_matrix
    J = lil_matrix((nvar, nvar), dtype=np.float64)

    # Phi equation couples to everything through the Poisson equation
    for j in range(nvar):
        J[IDX_PHI, j] = 1
    # Every equation depends on Phi (through Phi_dot and Psi)
    for i in range(nvar):
        J[i, IDX_PHI] = 1

    # delta_b, v_b couple to each other and to Theta_1
    J[IDX_DELTA_B, IDX_V_B] = 1
    J[IDX_V_B, IDX_DELTA_B] = 1
    J[IDX_V_B, idx_theta(1)] = 1
    J[idx_theta(1), IDX_V_B] = 1

    # delta_c, v_c couple to each other
    J[IDX_DELTA_C, IDX_V_C] = 1
    J[IDX_V_C, IDX_DELTA_C] = 1

    # Photon hierarchy: tridiagonal + self-coupling
    for l in range(lgmax + 1):
        idx_l = idx_theta(l)
        J[idx_l, idx_l] = 1  # self (collision)
        if l > 0:
            J[idx_l, idx_theta(l - 1)] = 1
        if l < lgmax:
            J[idx_l, idx_theta(l + 1)] = 1
    # Theta_0 couples to Theta_0, Theta_1 (and Phi through Phi_dot)
    # Theta_1 couples to Theta_0, Theta_2, Psi -> N_2
    J[idx_theta(1), _n_start + 2] = 1  # through Psi
    if lgmax >= 2:
        J[idx_theta(2), idx_theta(1)] = 1
        J[idx_theta(2), idx_theta(3)] = 1 if lgmax >= 3 else 0

    # Neutrino hierarchy: tridiagonal
    for l in range(lnmax + 1):
        idx_l = _n_start + l
        J[idx_l, idx_l] = 1
        if l > 0:
            J[idx_l, _n_start + l - 1] = 1
        if l < lnmax:
            J[idx_l, _n_start + l + 1] = 1
    # N_1 couples to Psi -> Theta_2
    J[_n_start + 1, idx_theta(2)] = 1 if lgmax >= 2 else 0

    # Phi equation couples to Theta_0, N_0, delta_b, delta_c, N_2, Theta_2
    J[IDX_PHI, idx_theta(0)] = 1
    J[IDX_PHI, _n_start] = 1
    J[IDX_PHI, IDX_DELTA_B] = 1
    J[IDX_PHI, IDX_DELTA_C] = 1
    J[IDX_PHI, _n_start + 2] = 1 if lnmax >= 2 else 0
    J[IDX_PHI, idx_theta(2)] = 1 if lgmax >= 2 else 0

    return J.tocsc()


# ============================================================================
# Initial conditions (single k-mode, numpy float64)
# ============================================================================

def _adiabatic_ic_single(k, tau_init, lgmax, lnmax, _n_start):
    """Adiabatic initial conditions for a single k-mode."""
    nvar = n_var_total(lgmax, lnmax)
    y0 = np.zeros(nvar, dtype=np.float64)

    y0[IDX_PHI] = 1.0
    y0[IDX_DELTA_B] = -1.5
    y0[IDX_V_B] = k * tau_init / 6.0
    y0[IDX_DELTA_C] = -1.5
    y0[IDX_V_C] = k * tau_init / 6.0
    y0[idx_theta(0)] = -0.5
    y0[idx_theta(1)] = k * tau_init / 18.0
    y0[_n_start + 0] = -0.5
    y0[_n_start + 1] = k * tau_init / 18.0
    if lnmax >= 2:
        y0[_n_start + 2] = (k * tau_init) ** 2 / 60.0
    for l in range(3, min(lnmax + 1, 6)):
        prod_val = 1.0
        for j in range(1, l + 1):
            prod_val *= (2 * j + 1)
        y0[_n_start + l] = (k * tau_init) ** l / prod_val

    return y0


# ============================================================================
# Core solver: solve one k-mode with scipy BDF
# ============================================================================

def solve_single_k(k, bg, tau_grid, tau_snapshots,
                   lgmax=25, lnmax=L_NU_MAX,
                   method='BDF', rtol=1e-8, atol=1e-10,
                   verbose=False):
    """
    Solve the Boltzmann hierarchy for a single k-mode using scipy.

    Parameters
    ----------
    k : float
        Wavenumber in Mpc^{-1}.
    bg : Background
        Solved background cosmology.
    tau_grid : ndarray
        Not used for stepping (scipy handles that), but defines the
        integration domain [tau_grid[0], tau_grid[-1]].
    tau_snapshots : ndarray
        Conformal times at which to record the state.
    lgmax : int
        Max photon multipole.
    lnmax : int
        Max neutrino multipole.
    method : str
        'BDF' or 'Radau' (both implicit, both handle stiffness).
    rtol, atol : float
        Tolerances for solve_ivp.
    verbose : bool

    Returns
    -------
    dict with keys:
        'Phi', 'Psi', 'Theta_0', 'v_b' : arrays of shape (N_snap,)
        'tau' : snapshot times
        'nfev' : number of RHS evaluations
        'success' : bool
    """
    _n_start = idx_n_start(lgmax)
    nvar = n_var_total(lgmax, lnmax)
    Og = _OMEGA_GAMMA
    On = _OMEGA_NU
    Ob = _OMEGA_B
    Oc = float(bg.Omega_cdm)
    H0 = _H0_MPC

    # Build RHS
    rhs = _make_rhs(k, bg, Og, On, Ob, Oc, H0, lgmax, lnmax, _n_start)

    # Initial conditions
    tau_init = float(tau_grid[0])
    tau_end = float(tau_grid[-1])
    y0 = _adiabatic_ic_single(k, tau_init, lgmax, lnmax, _n_start)

    # Solve.
    # NOTE: jac_sparsity is OMITTED -- for 52x52 the sparse LU overhead
    # is worse than dense. BDF with dense Jacobian takes ~0.1-0.5s per mode.
    t0 = time.time()
    sol = solve_ivp(
        rhs,
        t_span=(tau_init, tau_end),
        y0=y0,
        method=method,
        t_eval=tau_snapshots,
        rtol=rtol,
        atol=atol,
        dense_output=False,
        max_step=np.inf,  # let BDF choose
    )
    dt = time.time() - t0

    if verbose:
        status = "OK" if sol.success else f"FAILED: {sol.message}"
        print(f"  k={k:.4f}: {dt:.2f}s, nfev={sol.nfev}, status={status}")

    if not sol.success:
        print(f"  WARNING: k={k:.4f} failed: {sol.message}")

    # Extract snapshots
    N_snap = len(tau_snapshots)
    Phi_snap = np.zeros(N_snap)
    Psi_snap = np.zeros(N_snap)
    Theta_0_snap = np.zeros(N_snap)
    v_b_snap = np.zeros(N_snap)

    if sol.success and sol.y is not None:
        for i in range(min(N_snap, sol.y.shape[1])):
            y_i = sol.y[:, i]
            Phi_snap[i] = y_i[IDX_PHI]
            Theta_0_snap[i] = y_i[idx_theta(0)]
            v_b_snap[i] = y_i[IDX_V_B]

            # Compute Psi
            a_val = float(bg.a_at_tau(np.array([tau_snapshots[i]]))[0])
            N_2 = y_i[_n_start + 2] if lnmax >= 2 else 0.0
            Theta_2 = y_i[idx_theta(2)] if lgmax >= 2 else 0.0
            Psi_snap[i] = (Phi_snap[i]
                            - 12.0 * H0**2 / (a_val**2 * k**2)
                            * (On * N_2 + Og * Theta_2))

    return {
        'Phi': Phi_snap,
        'Psi': Psi_snap,
        'Theta_0': Theta_0_snap,
        'v_b': v_b_snap,
        'tau': tau_snapshots,
        'nfev': sol.nfev if sol.success else -1,
        'success': sol.success,
        'time': dt,
    }


# ============================================================================
# Batch solver: all k-modes (sequential, with progress)
# ============================================================================

def solve_all_k(bg, k_arr, tau_snapshots,
                lgmax=25, lnmax=L_NU_MAX,
                method='BDF', rtol=1e-8, atol=1e-10,
                verbose=True):
    """
    Solve the Boltzmann hierarchy for all k-modes sequentially.

    Parameters
    ----------
    bg : Background
        Solved background.
    k_arr : ndarray (N_k,)
        Wavenumbers in Mpc^{-1}.
    tau_snapshots : ndarray (N_snap,)
        Conformal times for snapshots.
    lgmax, lnmax : int
        Hierarchy truncation.
    method : str
        'BDF' or 'Radau'.
    rtol, atol : float

    Returns
    -------
    dict with keys:
        'Phi' : (N_snap, N_k)
        'Psi' : (N_snap, N_k)
        'Theta_0' : (N_snap, N_k)
        'v_b' : (N_snap, N_k)
        'tau' : (N_snap,)
        'k_arr' : (N_k,)
        'nfev_total' : int
        'n_failed' : int
        'total_time' : float
    """
    N_k = len(k_arr)
    N_snap = len(tau_snapshots)

    # Build tau_grid for integration domain
    k_max = float(k_arr[-1])
    tau_grid = build_tau_grid(bg, k_max, N_early=200, N_late=400)

    # Ensure tau_snapshots are within integration domain
    tau_snapshots = np.clip(tau_snapshots, tau_grid[0], tau_grid[-1])

    Phi_all = np.zeros((N_snap, N_k))
    Psi_all = np.zeros((N_snap, N_k))
    Theta_0_all = np.zeros((N_snap, N_k))
    v_b_all = np.zeros((N_snap, N_k))

    nfev_total = 0
    n_failed = 0
    t_start = time.time()

    if verbose:
        print(f"[ScipyBDF] Solving {N_k} k-modes with {method}")
        print(f"[ScipyBDF] lgmax={lgmax}, lnmax={lnmax}, "
              f"nvar={n_var_total(lgmax, lnmax)}")
        print(f"[ScipyBDF] rtol={rtol:.0e}, atol={atol:.0e}")
        print(f"[ScipyBDF] tau: [{tau_grid[0]:.2e}, {tau_grid[-1]:.1f}] Mpc, "
              f"{N_snap} snapshots")

    for ik, k in enumerate(k_arr):
        result = solve_single_k(
            k, bg, tau_grid, tau_snapshots,
            lgmax=lgmax, lnmax=lnmax,
            method=method, rtol=rtol, atol=atol,
            verbose=(verbose and (ik % max(N_k // 10, 1) == 0 or ik == N_k - 1)),
        )

        Phi_all[:, ik] = result['Phi']
        Psi_all[:, ik] = result['Psi']
        Theta_0_all[:, ik] = result['Theta_0']
        v_b_all[:, ik] = result['v_b']

        if result['success']:
            nfev_total += result['nfev']
        else:
            n_failed += 1

        # Progress
        if verbose and (ik + 1) % max(N_k // 10, 1) == 0:
            elapsed = time.time() - t_start
            eta = elapsed / (ik + 1) * (N_k - ik - 1)
            print(f"  [{ik+1}/{N_k}] elapsed={elapsed:.1f}s, ETA={eta:.1f}s, "
                  f"nfev={nfev_total}")

    total_time = time.time() - t_start

    if verbose:
        print(f"[ScipyBDF] Done: {total_time:.1f}s total, "
              f"{total_time/N_k:.2f}s/mode, "
              f"nfev={nfev_total}, failed={n_failed}/{N_k}")

    return {
        'Phi': Phi_all,
        'Psi': Psi_all,
        'Theta_0': Theta_0_all,
        'v_b': v_b_all,
        'tau': tau_snapshots,
        'k_arr': k_arr,
        'nfev_total': nfev_total,
        'n_failed': n_failed,
        'total_time': total_time,
    }


# ============================================================================
# C_l from snapshots: line-of-sight integration (reuses GPU Bessel if avail)
# ============================================================================

def compute_cl_from_snapshots(snapshots, bg, ell_values,
                              include_doppler=True, include_isw=True,
                              use_gpu_bessel=True):
    """
    Line-of-sight C_l from scipy BDF snapshots.

    This uses the SAME physics as AccurateResult.compute_cl_los but with
    NO amplitude correction D(k) -- since BDF handles the ODE accurately,
    we expect the raw source amplitudes to be correct.

    Parameters
    ----------
    snapshots : dict from solve_all_k
    bg : Background
    ell_values : ndarray of int
    include_doppler, include_isw : bool
    use_gpu_bessel : bool
        If True, try to use GPU Bessel tables for speed.

    Returns
    -------
    ell_values, Cl, Dl_muK2
    """
    from scipy.special import spherical_jn

    k_arr = snapshots['k_arr']
    N_k = len(k_arr)
    tau_s = snapshots['tau']
    Phi_s = snapshots['Phi']
    Psi_s = snapshots['Psi']
    Theta_0_s = snapshots['Theta_0']
    v_b_s = snapshots['v_b']

    valid = tau_s > 0
    tau_s = tau_s[valid]
    Phi_s = Phi_s[valid]
    Psi_s = Psi_s[valid]
    Theta_0_s = Theta_0_s[valid]
    v_b_s = v_b_s[valid]
    N_tau = len(tau_s)

    tau_0 = bg.tau_0
    chi_s = tau_0 - tau_s

    # Visibility function weights
    g_s = bg.visibility_at_tau(tau_s)
    dtau_s_diff = np.diff(tau_s)
    dtau_w = np.zeros(N_tau)
    if N_tau >= 2:
        dtau_w[0] = 0.5 * dtau_s_diff[0]
        dtau_w[-1] = 0.5 * dtau_s_diff[-1]
        for i in range(1, N_tau - 1):
            dtau_w[i] = 0.5 * (dtau_s_diff[i - 1] + dtau_s_diff[i])
    w_vis = g_s * dtau_w

    # Optical depth for ISW
    kappa_s = bg.kappa_at_tau(tau_s)
    exp_neg_kappa = np.exp(-kappa_s)

    # Phi' + Psi' from finite differences
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

    # Primordial power spectrum
    P_R = A_s * (k_arr / k_pivot) ** (n_s - 1.0)
    lnk = np.log(k_arr)
    dlnk = np.diff(lnk)

    # NO D(k) amplitude correction -- BDF should give correct amplitudes.
    # Only apply Silk damping from the background.
    k_D_eff = bg.k_D * np.sqrt(2.0)
    silk_eff = np.exp(-(k_arr / k_D_eff) ** 2)

    print(f"[ScipyBDF] C_l: {len(ell_values)} ells, {N_tau} tau-steps, "
          f"{N_k} k-modes")
    t0 = time.time()

    Cl = np.zeros(len(ell_values), dtype=np.float64)

    for il, ell in enumerate(ell_values):
        Delta_l = np.zeros(N_k, dtype=np.float64)

        for it in range(N_tau):
            x = k_arr * chi_s[it]
            jl = spherical_jn(int(ell), x)

            # Sachs-Wolfe: g * (Theta_0 + Psi) * silk * j_l
            source = (Theta_0_s[it] + Psi_s[it]) * silk_eff * jl

            if include_doppler:
                jlp = spherical_jn(int(ell), x, derivative=True)
                source += v_b_s[it] * silk_eff * jlp

            if include_isw:
                source += exp_neg_kappa[it] * PhiPsi_dot[it] * jl

            Delta_l += w_vis[it] * source

        integrand = P_R * Delta_l ** 2
        mid = 0.5 * (integrand[:-1] + integrand[1:])
        Cl[il] = 4.0 * np.pi * np.sum(mid * dlnk)

    Cl = np.maximum(Cl, 0.0)
    ell_f = ell_values.astype(np.float64)
    Dl = ell_f * (ell_f + 1.0) * Cl / (2.0 * np.pi) * (T_CMB * 1e6) ** 2

    print(f"[ScipyBDF] C_l done: {time.time() - t0:.1f}s")
    return ell_values, Cl, Dl


# ============================================================================
# Main: run the full pipeline and compare against CLASS
# ============================================================================

def main():
    import argparse

    parser = argparse.ArgumentParser(
        description='Scipy BDF Boltzmann solver -- accuracy reference')
    parser.add_argument('--N_k', type=int, default=50,
                        help='Number of k-modes (default: 50)')
    parser.add_argument('--method', type=str, default='BDF',
                        choices=['BDF', 'Radau'],
                        help='ODE solver method')
    parser.add_argument('--rtol', type=float, default=1e-8)
    parser.add_argument('--atol', type=float, default=1e-10)
    parser.add_argument('--lgmax', type=int, default=25,
                        help='Max photon multipole')
    parser.add_argument('--lnmax', type=int, default=L_NU_MAX,
                        help='Max neutrino multipole')
    parser.add_argument('--ell_max', type=int, default=2500)
    parser.add_argument('--class_file', type=str,
                        default='/Users/akaihuangm1/Desktop/github/'
                                'gdm_class_public/output/mlx_ref_lcdm_cl.dat')
    parser.add_argument('--recombination', type=str, default='peebles',
                        choices=['tanh', 'peebles'])
    parser.add_argument('--save_plot', type=str, default=None,
                        help='Save comparison plot to this path')
    args = parser.parse_args()

    print("=" * 70)
    print("SCIPY BDF REFERENCE SOLVER")
    print(f"Method: {args.method}, rtol={args.rtol:.0e}, atol={args.atol:.0e}")
    print(f"N_k={args.N_k}, lgmax={args.lgmax}, lnmax={args.lnmax}")
    print(f"Recombination: {args.recombination}")
    print("=" * 70)

    # ------------------------------------------------------------------
    # Step 1: Background
    # ------------------------------------------------------------------
    print("\n--- Step 1: Background ---")
    t0 = time.time()
    bg = Background(khronon=False, recombination=args.recombination)
    bg.solve()
    print(f"Background: {time.time() - t0:.2f}s")

    # ------------------------------------------------------------------
    # Step 2: Build k-grid and snapshot schedule
    # ------------------------------------------------------------------
    k_arr = np.geomspace(5e-4, 0.35, args.N_k).astype(np.float64)

    # Dense snapshot schedule around visibility peak
    tau_rec = bg.tau_rec
    tau_end = min(tau_rec + 100.0, bg.tau_0 * 0.5)
    tau_early = np.linspace(bg.tau_grid[1], tau_rec - 50, 50)
    tau_vis = np.linspace(tau_rec - 50, tau_end, 250)
    tau_snapshots = np.unique(np.concatenate([tau_early, tau_vis]))
    print(f"Snapshots: {len(tau_snapshots)} points, "
          f"tau=[{tau_snapshots[0]:.2e}, {tau_snapshots[-1]:.1f}] Mpc")

    # ------------------------------------------------------------------
    # Step 3: Solve with scipy BDF
    # ------------------------------------------------------------------
    print(f"\n--- Step 2: Scipy {args.method} solver ---")
    snapshots = solve_all_k(
        bg, k_arr, tau_snapshots,
        lgmax=args.lgmax, lnmax=args.lnmax,
        method=args.method, rtol=args.rtol, atol=args.atol,
    )

    # ------------------------------------------------------------------
    # Step 4: Compute C_l
    # ------------------------------------------------------------------
    print("\n--- Step 3: C_l integration ---")
    ell_values = np.unique(np.concatenate([
        np.arange(2, 30, 1),
        np.arange(30, 100, 2),
        np.arange(100, 500, 4),
        np.arange(500, 1500, 8),
        np.arange(1500, args.ell_max + 1, 12),
    ])).astype(int)

    ell_out, Cl, Dl = compute_cl_from_snapshots(
        snapshots, bg, ell_values,
        include_doppler=True, include_isw=True,
    )

    # ------------------------------------------------------------------
    # Step 5: Compare against CLASS
    # ------------------------------------------------------------------
    print("\n--- Step 4: Comparison vs CLASS ---")
    class_file = args.class_file
    if os.path.exists(class_file):
        class_data = np.loadtxt(class_file)
        class_ell = class_data[:, 0].astype(int)
        class_Dl = class_data[:, 1] * (T_CMB * 1e6) ** 2

        # Interpolate to common ell grid
        from scipy.interpolate import interp1d
        ell_common = ell_out[(ell_out >= class_ell[0]) & (ell_out <= class_ell[-1])]
        class_interp = interp1d(class_ell, class_Dl, kind='cubic')
        class_at_common = class_interp(ell_common)

        # Find our Dl at common ells
        our_interp = interp1d(ell_out, Dl, kind='cubic')
        our_at_common = our_interp(ell_common)

        # Restrict to ell >= 10 for meaningful comparison
        mask = ell_common >= 10
        ell_cmp = ell_common[mask]
        our_cmp = our_at_common[mask]
        class_cmp = class_at_common[mask]

        # Metrics
        safe = class_cmp > 0
        ratio = np.ones_like(our_cmp)
        ratio[safe] = our_cmp[safe] / class_cmp[safe]

        residual = np.zeros_like(our_cmp)
        residual[safe] = (our_cmp[safe] - class_cmp[safe]) / class_cmp[safe]

        rms = np.sqrt(np.mean(residual[safe] ** 2)) * 100
        max_err = np.max(np.abs(residual[safe])) * 100
        mean_ratio = np.mean(ratio[safe])

        print(f"\n{'='*60}")
        print(f"  SCIPY {args.method} vs CLASS (ell >= 10)")
        print(f"{'='*60}")
        print(f"  RMS residual:  {rms:.2f}%")
        print(f"  Max |residual|: {max_err:.2f}%")
        print(f"  Mean amplitude ratio: {mean_ratio:.4f}")
        print(f"  N_k = {args.N_k}, total time = {snapshots['total_time']:.1f}s")
        print(f"  Time per k-mode: {snapshots['total_time']/args.N_k:.2f}s")
        print(f"  Failed modes: {snapshots['n_failed']}/{args.N_k}")

        # Peak analysis
        from scipy.signal import find_peaks
        peaks_class, _ = find_peaks(class_Dl, height=100, distance=50)
        if len(peaks_class) >= 3:
            print(f"\n  CLASS peaks: ell = {class_ell[peaks_class[:5]]}")

        peaks_ours, _ = find_peaks(Dl, height=1, distance=20)
        if len(peaks_ours) >= 1:
            print(f"  Our peaks:   ell = {ell_out[peaks_ours[:5]]}")

        # Ratio at specific ells
        for target_ell in [220, 540, 820, 1130, 1420]:
            idx_c = np.argmin(np.abs(class_ell - target_ell))
            idx_o = np.argmin(np.abs(ell_out - target_ell))
            if class_Dl[idx_c] > 0:
                r = Dl[idx_o] / class_Dl[idx_c]
                print(f"  ell~{target_ell}: our={Dl[idx_o]:.1f}, "
                      f"CLASS={class_Dl[idx_c]:.1f}, ratio={r:.4f}")

        print(f"{'='*60}")
        print(f"\n  KEY QUESTION: Is RMS < 5%?")
        if rms < 5.0:
            print(f"  YES -> The Boltzmann equations + BDF solver achieve < 5% RMS.")
            print(f"         The IMEX solver was limiting accuracy.")
        elif rms < 20.0:
            print(f"  PARTIAL -> RMS is {rms:.1f}%, better than IMEX but not < 5%.")
            print(f"         Physics equations may need refinement.")
        else:
            print(f"  NO -> RMS is {rms:.1f}%, similar to IMEX.")
            print(f"        The problem is in the physics, not the ODE solver.")

        # ------------------------------------------------------------------
        # Optional: save plot
        # ------------------------------------------------------------------
        plot_path = args.save_plot
        if plot_path is None:
            plot_dir = os.path.dirname(os.path.abspath(__file__))
            plot_path = os.path.join(plot_dir, 'scipy_bdf_comparison.png')

        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt

            fig, axes = plt.subplots(3, 1, figsize=(10, 10),
                                     gridspec_kw={'height_ratios': [3, 1, 1]})

            # Panel 1: D_l spectra
            ax1 = axes[0]
            ax1.plot(class_ell, class_Dl, 'k-', lw=1.0, alpha=0.7, label='CLASS')
            ax1.plot(ell_out, Dl, 'r-', lw=1.0, alpha=0.8,
                     label=f'Scipy {args.method} (N_k={args.N_k})')
            ax1.set_ylabel(r'$D_\ell$ [$\mu K^2$]')
            ax1.set_xlim(2, args.ell_max)
            ax1.legend(fontsize=10)
            ax1.set_title(f'Scipy {args.method} vs CLASS  |  '
                          f'RMS={rms:.2f}%, N_k={args.N_k}, '
                          f'time={snapshots["total_time"]:.1f}s')

            # Panel 2: residual
            ax2 = axes[1]
            ax2.plot(ell_cmp, residual[safe] * 100, 'b-', lw=0.5, alpha=0.6)
            ax2.axhline(0, color='k', ls='--', lw=0.5)
            ax2.axhline(5, color='r', ls=':', lw=0.5, label='5%')
            ax2.axhline(-5, color='r', ls=':', lw=0.5)
            ax2.set_ylabel('Residual [%]')
            ax2.set_xlim(2, args.ell_max)
            ax2.set_ylim(-100, 100)
            ax2.legend(fontsize=8)

            # Panel 3: amplitude ratio
            ax3 = axes[2]
            ax3.plot(ell_cmp, ratio[safe], 'g-', lw=0.5, alpha=0.6)
            ax3.axhline(1, color='k', ls='--', lw=0.5)
            ax3.set_ylabel('Ratio (ours/CLASS)')
            ax3.set_xlabel(r'$\ell$')
            ax3.set_xlim(2, args.ell_max)
            ax3.set_ylim(0, 2)

            plt.tight_layout()
            plt.savefig(plot_path, dpi=150, bbox_inches='tight')
            print(f"\n  Plot saved: {plot_path}")
            plt.close()
        except ImportError:
            print("  (matplotlib not available, skipping plot)")

    else:
        print(f"  CLASS file not found: {class_file}")
        print("  Skipping comparison. Output D_l only:")
        for i in range(0, len(ell_out), max(len(ell_out) // 20, 1)):
            print(f"    ell={ell_out[i]:5d}  Dl={Dl[i]:.4f}")

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    print(f"\n--- Summary ---")
    print(f"  Method: {args.method}")
    print(f"  N_k: {args.N_k}")
    print(f"  Total ODE time: {snapshots['total_time']:.1f}s "
          f"({snapshots['total_time']/args.N_k:.2f}s/mode)")
    print(f"  Extrapolated 500 modes: "
          f"~{snapshots['total_time']/args.N_k * 500:.0f}s "
          f"(~{snapshots['total_time']/args.N_k * 500 / 60:.1f} min)")
    print(f"  Failed: {snapshots['n_failed']}/{args.N_k}")
    print(f"  Total RHS evaluations: {snapshots['nfev_total']}")

    return snapshots, ell_out, Cl, Dl


if __name__ == '__main__':
    main()
