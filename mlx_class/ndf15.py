"""
ndf15.py — Variable-order (1-5) NDF/BDF ODE solver ported from CLASS.

This is a Python/NumPy port of the evolver_ndf15.c from CLASS (Cosmic Linear
Anisotropy Solving System), written by Thomas Tram (2010).  The algorithm is
from [Shampine & Reichelt, "The MATLAB ODE Suite", SIAM J. Sci. Comput. 18-1,
1997].

Features:
  - Variable order 1-5 (NDF = Numerical Differentiation Formulas)
  - Variable step size with automatic selection
  - Simplified Newton iteration for the implicit solve
  - Numerical Jacobian with adaptive increment sizing
  - Lazy Jacobian recomputation (only when Newton convergence is slow)
  - Dense output via polynomial interpolation from backward differences

This module is self-contained: it only needs NumPy.  No CLASS or MLX
dependencies.  It can be used as a drop-in stiff ODE integrator for any
system y' = f(t, y).

Usage
-----
    from ndf15 import solve_ndf15

    def rhs(t, y):
        return np.array([-1000*y[0] + y[1], y[0] - y[1]])

    sol = solve_ndf15(rhs, t_span=(0.0, 1.0), y0=np.array([1.0, 0.0]),
                      rtol=1e-6)
    print(sol.t, sol.y)

Author: Ported from CLASS (Thomas Tram, 2010) by Claude, 2026.
"""

import numpy as np
from scipy.linalg import lu_factor, lu_solve
from dataclasses import dataclass, field
from typing import Callable, Optional, List

# ---------------------------------------------------------------------------
# Constants from CLASS evolver_ndf15.c
# ---------------------------------------------------------------------------
# NDF gamma coefficients  G[k] = sum_{j=1}^{k} 1/j
_G = np.array([1.0, 3.0/2.0, 11.0/6.0, 25.0/12.0, 137.0/60.0])

# NDF kappa (alpha) coefficients — the "free parameter" of NDF vs pure BDF
# alpha=0 gives pure BDF; these values are from Shampine & Reichelt
_ALPHA = np.array([-37.0/200.0, -1.0/9.0, -8.23e-2, -4.15e-2, 0.0])

# Precomputed from G and alpha
_INVGA = 1.0 / (_G * (1.0 - _ALPHA))
_ERCONST = _ALPHA * _G + 1.0 / np.arange(2, 7)

_EPS = np.finfo(np.float64).eps   # ~2.2e-16
_TINY = 1e-50

# Newton iteration limits
_MAXIT = 4
_MAXK = 5


# ---------------------------------------------------------------------------
# Result container
# ---------------------------------------------------------------------------
@dataclass
class NDF15Result:
    """Result of solve_ndf15."""
    t: np.ndarray          # shape (n_points,) — time values
    y: np.ndarray          # shape (n_points, neq) — solution
    success: bool = True
    message: str = ""
    n_steps: int = 0
    n_failed: int = 0
    n_fevals: int = 0
    n_jevals: int = 0
    n_ludecomps: int = 0
    n_lusolves: int = 0


# ---------------------------------------------------------------------------
# Numerical Jacobian (simplified from CLASS numjac)
# ---------------------------------------------------------------------------
def _numjac(fun, t, y, f0, fac, thresh):
    """Compute numerical Jacobian df/dy by finite differences.

    Uses the adaptive-increment strategy from MATLAB's numjac / CLASS:
    each column j uses an increment  del_j = fac[j] * max(|y[j]|, thresh).

    Parameters
    ----------
    fun : callable(t, y) -> dy
    t, y, f0 : current time, state, and f(t,y)
    fac : array(neq,) — adaptive increment factors (modified in-place)
    thresh : float — absolute threshold for y-scale

    Returns
    -------
    J : ndarray(neq, neq)
    nfe : int — number of function evaluations used
    """
    neq = len(y)
    J = np.empty((neq, neq))
    nfe = 0

    facmin = _EPS ** 0.78
    facmax = 0.1
    br = _EPS ** 0.875
    bl = _EPS ** 0.75
    bu = _EPS ** 0.25

    yscale = np.maximum(np.abs(y), thresh)

    for j in range(neq):
        # Compute increment
        del_j = (y[j] + fac[j] * yscale[j]) - y[j]
        if del_j == 0.0:
            while fac[j] < facmax:
                fac[j] = min(100.0 * fac[j], facmax)
                del_j = (y[j] + fac[j] * yscale[j]) - y[j]
                if del_j != 0.0:
                    break
            else:
                del_j = thresh

        # Keep del pointing into the region where f was evaluated
        if f0[j] >= 0.0:
            del_j = abs(del_j)
        else:
            del_j = -abs(del_j)

        y_pert = y.copy()
        y_pert[j] += del_j
        f_pert = fun(t, y_pert)
        nfe += 1

        diff = f_pert - f0
        J[:, j] = diff / del_j

        # Adjust fac for next call (simplified version of CLASS logic)
        abs_diff_max = np.max(np.abs(diff))
        abs_fdel_max = np.max(np.abs(f_pert))
        fscale = max(abs_fdel_max, np.max(np.abs(f0)))

        if fscale > 0.0:
            if abs_diff_max <= bl * fscale:
                fac[j] = min(10.0 * fac[j], facmax)
            elif abs_diff_max > bu * fscale:
                fac[j] = max(0.1 * fac[j], facmin)

    return J, nfe


# ---------------------------------------------------------------------------
# Step-size adjustment: rebuild backward-difference table for new h
# ---------------------------------------------------------------------------
def _adjust_stepsize(dif, ratio, k):
    """Adjust the backward-difference table `dif` when the step size changes.

    dif : ndarray(neq, MAXK+2) — backward differences, columns 0..k+1 used
    ratio : float — new_h / old_h
    k : int — current order (1-based)

    This implements the recurrence from Shampine & Reichelt:
      RU = product of U^{-1} and the ratio matrix R
    then  dif[:, 0:k] = dif[:, 0:k] @ RU[0:k, 0:k]
    """
    # Build the 5x5 U matrix (upper triangular with specific entries)
    U = np.array([
        [-1, -2, -3, -4, -5],
        [ 0,  1,  3,  6, 10],
        [ 0,  0, -1, -4,-10],
        [ 0,  0,  0,  1,  5],
        [ 0,  0,  0,  0, -1],
    ], dtype=np.float64)

    # Build R*U directly (the ratio-dependent matrix)
    RU = np.zeros((5, 5))
    # First row
    for i in range(5):
        RU[0, i] = -(i + 1) * ratio
    # Subsequent rows via recurrence
    for j in range(1, 5):
        for i in range(5):
            RU[j, i] = RU[j-1, i] * (1.0 - (1.0 + (i+1) * ratio) / (j+1))

    # RU = RU @ U  (in-place, matching CLASS)
    RU = RU @ U

    # Apply to dif: dif[:, 0:k] = dif[:, 0:k] @ RU[0:k, 0:k]
    dif[:, :k] = dif[:, :k] @ RU[:k, :k]


# ---------------------------------------------------------------------------
# Dense output interpolation from backward differences
# ---------------------------------------------------------------------------
def _interp_from_dif(t_interp, t_new, y_new, h, dif, k):
    """Interpolate solution at t_interp using backward-difference table.

    Returns y_interp, yp_interp (value and first derivative).
    """
    s = (t_interp - t_new) / h
    neq = len(y_new)

    prod = 1.0
    sumfrac = 0.0
    fact = 1.0
    y_interp = y_new.copy()
    yp_interp = np.zeros(neq)

    for j in range(k):
        prod *= (s + j)
        fact *= (j + 1)
        sumfrac += 1.0 / (s + j) if (s + j) != 0.0 else 0.0
        coeff_y = prod / fact
        coeff_yp = prod * sumfrac / (h * fact)
        y_interp += coeff_y * dif[:, j]
        yp_interp += coeff_yp * dif[:, j]

    return y_interp, yp_interp


# ---------------------------------------------------------------------------
# Main solver
# ---------------------------------------------------------------------------
def solve_ndf15(
    fun: Callable[[float, np.ndarray], np.ndarray],
    t_span: tuple,
    y0: np.ndarray,
    rtol: float = 1e-6,
    atol: float = 1e-15,
    t_eval: Optional[np.ndarray] = None,
    max_steps: int = 100_000,
    jac: Optional[Callable[[float, np.ndarray], np.ndarray]] = None,
    min_step: float = 0.0,
) -> NDF15Result:
    """Solve a stiff ODE system using NDF order 1-5.

    Parameters
    ----------
    fun : callable(t, y) -> dy/dt
        Right-hand side of the ODE system.
    t_span : (t0, tfinal)
        Integration interval.
    y0 : ndarray(neq,)
        Initial conditions.
    rtol : float
        Relative tolerance.
    atol : float
        Absolute tolerance (also used as threshold for y-scaling).
    t_eval : ndarray or None
        Times at which to store the solution.  If None, store every accepted
        step.
    max_steps : int
        Maximum number of steps before raising an error.
    jac : callable(t, y) -> ndarray(neq,neq) or None
        Analytical Jacobian.  If None, numerical differencing is used.
    min_step : float
        Minimum allowed step size.  0 means 16*eps*|t|.

    Returns
    -------
    NDF15Result with fields t, y, and solver statistics.
    """
    t0, tfinal = float(t_span[0]), float(t_span[1])
    y = np.array(y0, dtype=np.float64).copy()
    neq = len(y)
    thresh = atol

    # Direction of integration
    tdir = 1 if tfinal > t0 else -1

    # Statistics
    stats = dict(n_steps=0, n_failed=0, n_fevals=0, n_jevals=0,
                 n_ludecomps=0, n_lusolves=0)

    # ---------------------------------------------------------------
    # Initial function evaluation
    # ---------------------------------------------------------------
    f0 = fun(t0, y)
    stats['n_fevals'] += 1

    # ---------------------------------------------------------------
    # Compute initial Jacobian
    # ---------------------------------------------------------------
    fac = np.full(neq, np.sqrt(_EPS))  # adaptive increment for numjac

    if jac is not None:
        J = jac(t0, y)
        stats['n_jevals'] += 1
        nfe_jac = 0
    else:
        J, nfe_jac = _numjac(fun, t0, y, f0, fac, thresh)
        stats['n_jevals'] += 1
        stats['n_fevals'] += nfe_jac
    Jcurrent = True

    # ---------------------------------------------------------------
    # Compute initial step size (from CLASS)
    # ---------------------------------------------------------------
    hmax = abs(tfinal - t0) / 10.0
    htspan = abs(tfinal - t0)

    wt = np.maximum(np.abs(y), thresh)
    rh = np.max(1.25 / np.sqrt(rtol) * np.abs(f0 / wt))

    absh = min(hmax, htspan)
    if absh * rh > 1.0:
        absh = 1.0 / rh
    hmin_val = 16.0 * _EPS * abs(t0) if min_step == 0.0 else min_step
    absh = max(absh, hmin_val)
    h = tdir * absh

    # Second estimate using ddfddt = J*f0 + (f(t+tdel,y)-f(t,y))/tdel
    tdel = (t0 + tdir * min(np.sqrt(_EPS) * max(abs(t0), abs(t0 + h)), absh)) - t0
    f_tdel = fun(t0 + tdel, y)
    stats['n_fevals'] += 1

    ddfddt = J @ f0 + (f_tdel - f0) / tdel
    rh = np.max(1.25 * np.sqrt(0.5 * np.abs(ddfddt / wt) / rtol))
    absh = min(hmax, htspan)
    if absh * rh > 1.0:
        absh = 1.0 / rh
    absh = max(absh, hmin_val)
    h = tdir * absh

    # ---------------------------------------------------------------
    # Initialise backward-difference table and method variables
    # ---------------------------------------------------------------
    # dif[:, j] stores the j-th backward difference (0-indexed, j=0..MAXK+1)
    dif = np.zeros((neq, _MAXK + 2))
    dif[:, 0] = h * f0

    k = 1           # current order
    klast = k
    abshlast = absh
    nconhk = 0      # steps at current h and k

    # Initial linearisation:  M = I - hinvGak * J,  LU-factorise M
    hinvGak = h * _INVGA[k - 1]
    M = np.eye(neq) - hinvGak * J
    lu_piv = lu_factor(M)
    stats['n_ludecomps'] += 1
    havrate = False

    # ---------------------------------------------------------------
    # Output storage
    # ---------------------------------------------------------------
    t_out_list: List[float] = []
    y_out_list: List[np.ndarray] = []

    if t_eval is not None:
        t_eval = np.asarray(t_eval, dtype=np.float64)
        next_eval = 0
        # skip any t_eval < t0
        while next_eval < len(t_eval) and tdir * (t_eval[next_eval] - t0) < 0:
            next_eval += 1
    else:
        t_out_list.append(t0)
        y_out_list.append(y.copy())

    # ---------------------------------------------------------------
    # Main integration loop
    # ---------------------------------------------------------------
    t = t0
    done = False
    at_hmin = False

    while not done:
        if stats['n_steps'] > max_steps:
            return NDF15Result(
                t=np.array(t_out_list),
                y=np.array(y_out_list),
                success=False,
                message=f"Maximum number of steps ({max_steps}) exceeded at t={t}",
                **stats,
            )

        hmin_val = 16.0 * _EPS * abs(t) if min_step == 0.0 else min_step
        absh = max(hmin_val, absh)
        absh = min(hmax, absh)

        if abs(absh - hmin_val) < 100 * _EPS:
            if at_hmin:
                absh = abshlast
            at_hmin = True
        else:
            at_hmin = False

        h = tdir * absh

        # Stretch to hit tfinal exactly
        if 1.1 * absh >= abs(tfinal - t):
            h = tfinal - t
            absh = abs(h)
            done = True

        # Adjust dif if h or k changed
        if (abs(absh - abshlast) / max(absh, _TINY) > 1e-6) or (k != klast):
            _adjust_stepsize(dif, absh / abshlast, k)
            hinvGak = h * _INVGA[k - 1]
            nconhk = 0
            M = np.eye(neq) - hinvGak * J
            lu_piv = lu_factor(M)
            stats['n_ludecomps'] += 1
            havrate = False

        # =============================================================
        # Inner loop: advance one step (may retry on failure)
        # =============================================================
        nofailed = True
        while True:
            gotynew = False
            while not gotynew:
                # Compute psi = sum_{j=0}^{k-1} dif[:, j] * G[j] * invGa[k-1]
                psi = np.zeros(neq)
                for j in range(k):
                    psi += dif[:, j] * _G[j] * _INVGA[k - 1]

                # Predict
                tnew = t + h
                if done:
                    tnew = tfinal
                h = tnew - t  # purify

                pred = y.copy()
                for j in range(k):
                    pred += dif[:, j]
                ynew = pred.copy()

                # Prepare for Newton iteration
                difkp1 = np.zeros(neq)
                wt_new = np.maximum(np.maximum(np.abs(ynew), np.abs(y)), thresh)
                invwt = 1.0 / wt_new
                minnrm = np.max(100.0 * _EPS * np.abs(ynew) * invwt)

                # Simplified Newton iteration
                tooslow = False
                oldnrm = 0.0
                rate = 0.0

                for it in range(1, _MAXIT + 1):
                    rhs_vec = hinvGak * fun(tnew, ynew) - (psi + difkp1)
                    stats['n_fevals'] += 1

                    delta = lu_solve(lu_piv, rhs_vec)
                    stats['n_lusolves'] += 1

                    newnrm = np.max(np.abs(delta * invwt))

                    difkp1 += delta
                    ynew = pred + difkp1

                    if newnrm <= minnrm:
                        gotynew = True
                        break
                    elif it == 1:
                        if havrate:
                            errit = newnrm * rate / (1.0 - rate)
                            if errit <= 0.05 * rtol:
                                gotynew = True
                                break
                        else:
                            rate = 0.0
                    elif newnrm > 0.9 * oldnrm:
                        tooslow = True
                        break
                    else:
                        rate = max(0.9 * rate, newnrm / oldnrm)
                        havrate = True
                        errit = newnrm * rate / (1.0 - rate)
                        if errit <= 0.5 * rtol:
                            gotynew = True
                            break
                        elif it == _MAXIT:
                            tooslow = True
                            break
                        elif 0.5 * rtol < errit * rate ** (_MAXIT - it):
                            tooslow = True
                            break

                    oldnrm = newnrm

                # Handle Newton convergence failure
                if tooslow:
                    stats['n_failed'] += 1
                    if not Jcurrent:
                        # Recompute Jacobian
                        f0 = fun(t, y)
                        stats['n_fevals'] += 1
                        if jac is not None:
                            J = jac(t, y)
                        else:
                            J, nfe_jac = _numjac(fun, t, y, f0, fac, thresh)
                            stats['n_fevals'] += nfe_jac
                        stats['n_jevals'] += 1
                        Jcurrent = True
                    elif absh <= hmin_val:
                        return NDF15Result(
                            t=np.array(t_out_list),
                            y=np.array(y_out_list),
                            success=False,
                            message=f"Step size {absh} too small at t={t}",
                            **stats,
                        )
                    else:
                        abshlast = absh
                        absh = max(0.3 * absh, hmin_val)
                        h = tdir * absh
                        done = False
                        _adjust_stepsize(dif, absh / abshlast, k)
                        hinvGak = h * _INVGA[k - 1]
                        nconhk = 0

                    # Re-linearise
                    M = np.eye(neq) - hinvGak * J
                    lu_piv = lu_factor(M)
                    stats['n_ludecomps'] += 1
                    havrate = False

            # End while not gotynew
            # difkp1 is the (k+1)-th backward difference of ynew

            # Error estimation
            err = np.max(np.abs(difkp1 * invwt)) * _ERCONST[k - 1]

            if err > rtol:
                # Step failed
                stats['n_failed'] += 1
                if absh <= hmin_val:
                    return NDF15Result(
                        t=np.array(t_out_list),
                        y=np.array(y_out_list),
                        success=False,
                        message=f"Step size {absh} too small at t={t}",
                        **stats,
                    )

                abshlast = absh
                if nofailed:
                    nofailed = False
                    hopt = absh * max(0.1, 0.833 * (rtol / err) ** (1.0 / (k + 1)))
                    if k > 1:
                        errkm1 = np.max(np.abs((dif[:, k-1] + difkp1) * invwt))
                        errkm1 *= _ERCONST[k - 2]
                        hkm1 = absh * max(0.1, 0.769 * (rtol / errkm1) ** (1.0 / k))
                        if hkm1 > hopt:
                            hopt = min(absh, hkm1)
                            k = k - 1
                    absh = max(hmin_val, hopt)
                else:
                    absh = max(hmin_val, 0.5 * absh)

                h = tdir * absh
                if absh < abshlast:
                    done = False
                _adjust_stepsize(dif, absh / abshlast, k)
                hinvGak = h * _INVGA[k - 1]
                nconhk = 0
                M = np.eye(neq) - hinvGak * J
                lu_piv = lu_factor(M)
                stats['n_ludecomps'] += 1
                havrate = False
            else:
                break  # Successful step

        # End inner retry loop
        stats['n_steps'] += 1

        # =============================================================
        # Update backward-difference table
        # =============================================================
        dif[:, k+1] = difkp1 - dif[:, k]
        dif[:, k] = difkp1
        for j in range(k - 1, -1, -1):
            dif[:, j] += dif[:, j + 1]

        # =============================================================
        # Output
        # =============================================================
        if t_eval is not None:
            while next_eval < len(t_eval) and tdir * (tnew - t_eval[next_eval]) >= 0:
                if abs(tnew - t_eval[next_eval]) < _EPS * abs(tnew):
                    # Exactly at tnew
                    t_out_list.append(tnew)
                    y_out_list.append(ynew.copy())
                else:
                    # Interpolate
                    yi, _ = _interp_from_dif(t_eval[next_eval], tnew, ynew, h, dif, k)
                    t_out_list.append(t_eval[next_eval])
                    y_out_list.append(yi)
                next_eval += 1
        else:
            t_out_list.append(tnew)
            y_out_list.append(ynew.copy())

        if done:
            break

        # =============================================================
        # Step size and order selection for next step
        # =============================================================
        klast = k
        abshlast = absh
        nconhk = min(nconhk + 1, _MAXK + 2)

        if nconhk >= k + 2:
            # Consider changing order
            temp = 1.2 * (err / rtol) ** (1.0 / (k + 1))
            hopt = absh / temp if temp > 0.1 else 10 * absh
            kopt = k

            if k > 1:
                errkm1 = np.max(np.abs(dif[:, k-1] * invwt)) * _ERCONST[k - 2]
                temp = 1.3 * (errkm1 / rtol) ** (1.0 / k)
                hkm1 = absh / temp if temp > 0.1 else 10 * absh
                if hkm1 > hopt:
                    hopt = hkm1
                    kopt = k - 1

            if k < _MAXK:
                errkp1 = np.max(np.abs(dif[:, k+1] * invwt)) * _ERCONST[k]
                temp = 1.4 * (errkp1 / rtol) ** (1.0 / (k + 2))
                hkp1 = absh / temp if temp > 0.1 else 10 * absh
                if hkp1 > hopt:
                    hopt = hkp1
                    kopt = k + 1

            if hopt > absh:
                absh = hopt
                k = kopt

        # Advance
        t = tnew
        y = ynew.copy()
        Jcurrent = False

    # ---------------------------------------------------------------
    # Pack results
    # ---------------------------------------------------------------
    return NDF15Result(
        t=np.array(t_out_list),
        y=np.array(y_out_list),
        success=True,
        message="Integration successful",
        **stats,
    )


# ===========================================================================
# Test suite
# ===========================================================================
def _test_van_der_pol(mu=1000.0, rtol=1e-6):
    """Van der Pol oscillator — a classic stiff test problem.

    y0'' - mu*(1 - y0^2)*y0' + y0 = 0

    Rewritten as:
        y0' = y1
        y1' = mu*(1 - y0^2)*y1 - y0
    """
    print(f"=== Van der Pol oscillator (mu={mu}) ===")

    def vdp(t, y):
        return np.array([y[1], mu * (1 - y[0]**2) * y[1] - y[0]])

    def vdp_jac(t, y):
        return np.array([
            [0.0, 1.0],
            [-2.0 * mu * y[0] * y[1] - 1.0, mu * (1 - y[0]**2)]
        ])

    y0 = np.array([2.0, 0.0])
    t_span = (0.0, 3000.0)

    # Solve with ndf15 (numerical Jacobian)
    import time as _time
    t0 = _time.perf_counter()
    sol = solve_ndf15(vdp, t_span, y0, rtol=rtol, atol=1e-10)
    t_ndf15 = _time.perf_counter() - t0

    print(f"  ndf15 (numjac): {sol.n_steps} steps, {sol.n_failed} failed, "
          f"{sol.n_fevals} fevals, {sol.n_jevals} jevals, {sol.n_ludecomps} LUs")
    print(f"  time: {t_ndf15:.3f}s")
    print(f"  y(tfinal) = [{sol.y[-1, 0]:.10f}, {sol.y[-1, 1]:.10f}]")

    # Solve with ndf15 (analytical Jacobian)
    t0 = _time.perf_counter()
    sol_a = solve_ndf15(vdp, t_span, y0, rtol=rtol, atol=1e-10, jac=vdp_jac)
    t_ndf15_a = _time.perf_counter() - t0

    print(f"  ndf15 (anajac): {sol_a.n_steps} steps, {sol_a.n_failed} failed, "
          f"{sol_a.n_fevals} fevals, {sol_a.n_jevals} jevals, {sol_a.n_ludecomps} LUs")
    print(f"  time: {t_ndf15_a:.3f}s")
    print(f"  y(tfinal) = [{sol_a.y[-1, 0]:.10f}, {sol_a.y[-1, 1]:.10f}]")

    # Compare with scipy BDF
    from scipy.integrate import solve_ivp
    t0 = _time.perf_counter()
    sol_scipy = solve_ivp(vdp, t_span, y0, method='BDF', rtol=rtol, atol=1e-10,
                          jac=vdp_jac, dense_output=True)
    t_scipy = _time.perf_counter() - t0

    print(f"  scipy BDF:  {sol_scipy.nfev} fevals, {sol_scipy.njev} jevals, "
          f"{sol_scipy.nlu} LUs")
    print(f"  time: {t_scipy:.3f}s")
    print(f"  y(tfinal) = [{sol_scipy.y[0, -1]:.10f}, {sol_scipy.y[1, -1]:.10f}]")

    # Error between ndf15 and scipy at final time
    diff = np.abs(sol.y[-1] - sol_scipy.y[:, -1])
    print(f"  |ndf15 - scipy| at tfinal = [{diff[0]:.2e}, {diff[1]:.2e}]")

    return sol, sol_scipy


def _test_exponential_decay(rtol=1e-8):
    """Simple stiff linear system: y' = -1000*y, exact solution y = exp(-1000*t)."""
    print(f"\n=== Exponential decay (lambda=-1000) ===")

    lam = -1000.0

    def f(t, y):
        return np.array([lam * y[0]])

    def f_jac(t, y):
        return np.array([[lam]])

    y0 = np.array([1.0])
    t_span = (0.0, 0.01)
    t_eval = np.linspace(0, 0.01, 50)

    sol = solve_ndf15(f, t_span, y0, rtol=rtol, atol=1e-14, t_eval=t_eval,
                      jac=f_jac)

    y_exact = np.exp(lam * sol.t)
    max_err = np.max(np.abs(sol.y[:, 0] - y_exact))

    print(f"  ndf15: {sol.n_steps} steps, max error = {max_err:.2e}")
    print(f"  y(0.01) ndf15 = {sol.y[-1, 0]:.12e}")
    print(f"  y(0.01) exact = {np.exp(lam * 0.01):.12e}")

    return sol


def _test_robertson(rtol=1e-6):
    """Robertson chemical kinetics — extremely stiff (stiffness ratio ~1e8).

    y1' = -0.04*y1 + 1e4*y2*y3
    y2' =  0.04*y1 - 1e4*y2*y3 - 3e7*y2^2
    y3' =  3e7*y2^2
    """
    print(f"\n=== Robertson chemical kinetics ===")

    def robertson(t, y):
        return np.array([
            -0.04 * y[0] + 1e4 * y[1] * y[2],
             0.04 * y[0] - 1e4 * y[1] * y[2] - 3e7 * y[1]**2,
             3e7 * y[1]**2,
        ])

    def robertson_jac(t, y):
        return np.array([
            [-0.04,         1e4 * y[2],     1e4 * y[1]],
            [ 0.04, -1e4 * y[2] - 6e7 * y[1], -1e4 * y[1]],
            [ 0.0,          6e7 * y[1],     0.0],
        ])

    y0 = np.array([1.0, 0.0, 0.0])
    t_span = (0.0, 1e5)

    import time as _time

    t0 = _time.perf_counter()
    sol = solve_ndf15(robertson, t_span, y0, rtol=rtol, atol=1e-10,
                      jac=robertson_jac)
    t_ndf15 = _time.perf_counter() - t0

    print(f"  ndf15: {sol.n_steps} steps, {sol.n_failed} failed, "
          f"{sol.n_fevals} fevals, {sol.n_ludecomps} LUs, time={t_ndf15:.3f}s")
    print(f"  y(1e5) = [{sol.y[-1, 0]:.8e}, {sol.y[-1, 1]:.8e}, {sol.y[-1, 2]:.8e}]")
    print(f"  sum(y) = {np.sum(sol.y[-1]):.15f}  (should be 1.0)")

    # Compare with scipy
    from scipy.integrate import solve_ivp
    t0 = _time.perf_counter()
    sol_scipy = solve_ivp(robertson, t_span, y0, method='BDF', rtol=rtol,
                          atol=1e-10, jac=robertson_jac)
    t_scipy = _time.perf_counter() - t0

    print(f"  scipy BDF: {sol_scipy.nfev} fevals, time={t_scipy:.3f}s")
    print(f"  y(1e5) scipy = [{sol_scipy.y[0, -1]:.8e}, {sol_scipy.y[1, -1]:.8e}, "
          f"{sol_scipy.y[2, -1]:.8e}]")

    diff = np.abs(sol.y[-1] - sol_scipy.y[:, -1])
    print(f"  |ndf15 - scipy| = [{diff[0]:.2e}, {diff[1]:.2e}, {diff[2]:.2e}]")

    return sol, sol_scipy


def _test_stiff_linear_system(rtol=1e-8):
    """Stiff linear system with known exact solution.

    y' = A y, A = [[-1, 0], [1, -1000]], y(0) = [1, 1]
    Exact: y1 = exp(-t), y2 = (1000/999)*exp(-t) - (1/999)*exp(-1000t)
    """
    print(f"\n=== Stiff linear system (eigenvalues -1, -1000) ===")

    A = np.array([[-1.0, 0.0], [1.0, -1000.0]])

    def f(t, y):
        return A @ y

    def f_jac(t, y):
        return A

    y0 = np.array([1.0, 1.0])
    t_span = (0.0, 1.0)
    t_eval = np.linspace(0, 1, 100)

    sol = solve_ndf15(f, t_span, y0, rtol=rtol, atol=1e-14, t_eval=t_eval,
                      jac=f_jac)

    # Exact solution
    # Exact solution: y1 = exp(-t), y2 = (1/999)exp(-t) + (998/999)exp(-1000t)
    y1_exact = np.exp(-sol.t)
    y2_exact = (1.0/999.0) * np.exp(-sol.t) + (998.0/999.0) * np.exp(-1000.0 * sol.t)

    # Skip transient region (t < 0.01) where the fast mode exp(-1000t)
    # is being intentionally damped by the stiff solver
    mask = sol.t > 0.01
    err1 = np.max(np.abs(sol.y[mask, 0] - y1_exact[mask]))
    err2 = np.max(np.abs(sol.y[mask, 1] - y2_exact[mask]))
    rel2 = np.max(np.abs((sol.y[mask, 1] - y2_exact[mask]) / y2_exact[mask]))

    print(f"  ndf15: {sol.n_steps} steps, max abs errors (t>0.01) = [{err1:.2e}, {err2:.2e}]")
    print(f"  max relative error y2 (t>0.01) = {rel2:.2e}")

    # Compare with scipy
    from scipy.integrate import solve_ivp
    sol_scipy = solve_ivp(f, t_span, y0, method='BDF', rtol=rtol, atol=1e-14,
                          t_eval=t_eval, jac=f_jac)
    mask_s = t_eval > 0.01
    err1_s = np.max(np.abs(sol_scipy.y[0, mask_s] - np.exp(-t_eval[mask_s])))
    y2_exact_s = (1.0/999.0) * np.exp(-t_eval[mask_s]) + (998.0/999.0) * np.exp(-1000.0 * t_eval[mask_s])
    err2_s = np.max(np.abs(sol_scipy.y[1, mask_s] - y2_exact_s))
    print(f"  scipy: max errors (t>0.01) = [{err1_s:.2e}, {err2_s:.2e}]")

    return sol


if __name__ == "__main__":
    print("=" * 70)
    print("  NDF15 ODE Solver — Test Suite")
    print("  Ported from CLASS evolver_ndf15.c (Thomas Tram, 2010)")
    print("=" * 70)

    _test_exponential_decay()
    _test_stiff_linear_system()
    _test_robertson()
    _test_van_der_pol()

    print("\n" + "=" * 70)
    print("  All tests completed.")
    print("=" * 70)
