"""
ndf15_batched.py — Batched variable-order (1-5) NDF/BDF ODE solver on Apple Silicon GPU.

Solves N_k independent ODE systems simultaneously using MLX.  All systems
share a COMMON adaptive time grid (step size controlled by the most
restrictive mode), but each system has its own dynamics (e.g. different k
wavenumbers in a Boltzmann solver).

Algorithm: identical to ndf15.py (Shampine-Reichelt NDF, ported from CLASS
evolver_ndf15.c by Thomas Tram), but every array operation is lifted to carry
a leading batch dimension N_k and executed via MLX.

GPU / CPU split
---------------
- RHS evaluations, error norms, predictions, backward-difference updates,
  Newton residuals  -> GPU  (matmul, elementwise, reductions)
- Matrix inverse (mx.linalg.inv) -> CPU stream  (MLX linalg not yet on GPU
  as of 2026-03).  We precompute M_inv = inv(I - hinvGak * J) on CPU,
  then use GPU matmul M_inv @ rhs for each Newton iteration.

Performance
-----------
- Backward-difference table is stored as a Python list of MLX arrays to avoid
  expensive 3D slice-and-reconstruct operations.
- mx.eval() calls are minimised to ~2 per step (Newton convergence check and
  dif update).  All intermediate graph nodes are computed lazily.
- Target: 500 modes x 50 vars x ~1000 steps in seconds on Apple M-series.

Usage
-----
    from mlx_class.ndf15_batched import BatchedNDF15, solve_ndf15_batched

    solver = BatchedNDF15(fun, (t0, tf), y0, rtol=1e-6, atol=1e-8, jac=jac)
    result = solver.solve(t_eval=t_eval_array)
    # result.t : shape (n_out,)
    # result.y : shape (n_out, N_k, n_var)

Author: Claude, 2026.  Based on ndf15.py (CLASS port).
"""

import mlx.core as mx
import numpy as np
from dataclasses import dataclass
from typing import Callable, Optional

# ---------------------------------------------------------------------------
# NDF constants — Python floats to avoid numpy contamination of MLX arrays
# ---------------------------------------------------------------------------
_G_np = np.array([1.0, 3.0/2.0, 11.0/6.0, 25.0/12.0, 137.0/60.0])
_ALPHA_np = np.array([-37.0/200.0, -1.0/9.0, -8.23e-2, -4.15e-2, 0.0])
_INVGA_np = 1.0 / (_G_np * (1.0 - _ALPHA_np))
_ERCONST_np = _ALPHA_np * _G_np + 1.0 / np.arange(2, 7)

_G = [float(x) for x in _G_np]
_INVGA = [float(x) for x in _INVGA_np]
_ERCONST = [float(x) for x in _ERCONST_np]

_EPS = float(np.finfo(np.float64).eps)
_EPS32 = float(np.finfo(np.float32).eps)
_TINY = 1e-50
_MAXIT = 4
_MAXK = 5


@dataclass
class BatchedNDF15Result:
    """Result of batched NDF15 solve.

    t : mx.array (n_out,) — output time points
    y : mx.array (n_out, N_k, n_var) — solution
    """
    t: mx.array
    y: mx.array
    success: bool = True
    message: str = ""
    n_steps: int = 0
    n_failed: int = 0
    n_fevals: int = 0
    n_jevals: int = 0


# ---------------------------------------------------------------------------
# Batched numerical Jacobian
# ---------------------------------------------------------------------------
def _batched_numjac(fun, t, y, f0, thresh):
    """Compute batched numerical Jacobian by finite differences.

    Returns J (N_k, n_var, n_var) and number of function evaluations.
    """
    N_k, n_var = y.shape
    eps_sqrt = float(np.sqrt(_EPS32))

    J_cols = []
    yscale = mx.maximum(mx.abs(y), thresh)

    for j in range(n_var):
        delta_j = eps_sqrt * yscale[:, j]
        delta_j = mx.where(delta_j < _TINY, thresh, delta_j)
        y_pert = y.at[:, j].add(delta_j)
        f_pert = fun(t, y_pert)
        col_j = (f_pert - f0) / delta_j[:, None]
        J_cols.append(col_j)

    J = mx.stack(J_cols, axis=2)
    mx.eval(J)
    return J, n_var


# ---------------------------------------------------------------------------
# Step-size adjustment of backward-difference table (list-of-arrays version)
# ---------------------------------------------------------------------------
def _adjust_dif(dif_list, ratio, k):
    """Adjust backward-difference table for new step size.

    dif_list : list of mx.array, each (N_k, n_var), length MAXK+2
    ratio : float (new_h / old_h)
    k : int (current order, 1-based)

    Returns new dif_list.
    """
    RU = np.zeros((5, 5))
    for i in range(5):
        RU[0, i] = -(i + 1) * ratio
    for jr in range(1, 5):
        for i in range(5):
            RU[jr, i] = RU[jr - 1, i] * (
                1.0 - (1.0 + (i + 1) * ratio) / (jr + 1))
    U = np.array([
        [-1, -2, -3, -4, -5],
        [ 0,  1,  3,  6, 10],
        [ 0,  0, -1, -4,-10],
        [ 0,  0,  0,  1,  5],
        [ 0,  0,  0,  0, -1],
    ], dtype=np.float64)
    RU = (RU @ U)[:k, :k]  # (k, k)

    # dif_new[i] = sum_j RU[j, i] * dif_old[j]  for i, j in 0..k-1
    new_list = list(dif_list)  # shallow copy
    for i in range(k):
        val = RU[0, i] * dif_list[0]
        for j in range(1, k):
            val = val + RU[j, i] * dif_list[j]
        new_list[i] = val
    return new_list


# ---------------------------------------------------------------------------
# Dense output interpolation
# ---------------------------------------------------------------------------
def _interp_from_dif(t_interp, t_new, y_new, h, dif_list, k):
    """Interpolate batched solution at t_interp."""
    s = (t_interp - t_new) / h
    y_interp = y_new
    prod = 1.0
    fact = 1.0
    for j in range(k):
        prod *= (s + j)
        fact *= (j + 1)
        y_interp = y_interp + (prod / fact) * dif_list[j]
    return y_interp


# ---------------------------------------------------------------------------
# Main solver class
# ---------------------------------------------------------------------------
class BatchedNDF15:
    """
    Batched NDF15 (variable-order BDF) solver for N_k independent ODE systems.

    Parameters
    ----------
    fun : callable(t, y) -> dy
        t is Python float, y is mx.array (N_k, n_var). Returns (N_k, n_var).
    t_span : (float, float)
    y0 : mx.array (N_k, n_var)
    rtol, atol : float
    jac : callable(t, y) -> J (N_k, n_var, n_var), optional
    max_steps, min_step : int, float
    """

    def __init__(self, fun, t_span, y0, rtol=1e-6, atol=1e-8,
                 jac=None, max_steps=100_000, min_step=0.0):
        self.fun = fun
        self.t0, self.tfinal = float(t_span[0]), float(t_span[1])
        if not isinstance(y0, mx.array):
            y0 = mx.array(y0, dtype=mx.float32)
        self.y0 = y0
        self.N_k, self.n_var = self.y0.shape
        self.rtol = rtol
        self.atol = atol
        self.thresh = atol
        self.jac_fn = jac
        self.max_steps = max_steps
        self.min_step = min_step

    def _get_jacobian(self, t, y, f0):
        if self.jac_fn is not None:
            J = self.jac_fn(t, y)
            mx.eval(J)
            return J, 0
        return _batched_numjac(self.fun, t, y, f0, self.thresh)

    def _get_M_inv(self, hinvGak, J):
        """inv(I - hinvGak * J) on CPU, result used on GPU via matmul."""
        M = mx.eye(self.n_var) - hinvGak * J
        M_inv = mx.linalg.inv(M, stream=mx.cpu)
        mx.eval(M_inv)
        return M_inv

    def solve(self, t_eval=None):
        t0, tfinal = self.t0, self.tfinal
        N_k, n_var = self.N_k, self.n_var
        thresh = self.thresh
        rtol = self.rtol
        tdir = 1 if tfinal > t0 else -1
        eps = _EPS32

        # Stats
        nst = 0; nfl = 0; nfe = 0; nje = 0

        y = self.y0
        t = t0

        f0 = self.fun(t, y); mx.eval(f0); nfe += 1
        J, nj = self._get_jacobian(t, y, f0); nje += 1; nfe += nj
        Jcurrent = True

        # ----- initial step size -----
        hmax = abs(tfinal - t0) / 10.0
        wt = mx.maximum(mx.abs(y), thresh)
        rh = float(mx.max(1.25 / float(np.sqrt(rtol)) * mx.abs(f0) / wt))
        absh = min(hmax, abs(tfinal - t0))
        if absh * rh > 1.0:
            absh = 1.0 / rh
        hmin0 = 16.0 * eps * abs(t0) if self.min_step == 0.0 else self.min_step
        absh = max(absh, hmin0)
        h = tdir * absh

        tdel_off = tdir * min(float(np.sqrt(eps)) * max(abs(t0), abs(t0 + h)), absh)
        tdel = (t0 + tdel_off) - t0
        if abs(tdel) < _TINY:
            tdel = tdir * float(np.sqrt(eps))
        f_tdel = self.fun(t0 + tdel, y); mx.eval(f_tdel); nfe += 1
        Jf0 = (J @ f0[:, :, None])[:, :, 0]
        ddfddt = Jf0 + (f_tdel - f0) / tdel
        rh = float(mx.max(1.25 * mx.sqrt(mx.maximum(
            0.5 * mx.abs(ddfddt) / (wt * rtol), 0.0))))
        absh = min(hmax, abs(tfinal - t0))
        if absh * rh > 1.0:
            absh = 1.0 / rh
        absh = max(absh, hmin0)
        h = tdir * absh

        # ----- backward-difference table as list -----
        zero = mx.zeros((N_k, n_var))
        dif = [h * f0] + [zero] * (_MAXK + 1)  # indices 0.._MAXK+1
        mx.eval(dif[0])

        k = 1; klast = 1; abshlast = absh; nconhk = 0
        hinvGak = h * _INVGA[0]
        M_inv = self._get_M_inv(hinvGak, J)
        havrate = False

        # ----- output -----
        t_out = []; y_out = []
        if t_eval is not None:
            tev = np.asarray(t_eval, dtype=np.float64) if not isinstance(
                t_eval, mx.array) else np.array(t_eval.tolist())
            nxev = 0
            while nxev < len(tev) and tdir * (tev[nxev] - t0) < 0:
                nxev += 1
        else:
            tev = None
            t_out.append(t0); y_out.append(y)

        done = False; at_hmin = False

        # ================================================================
        #  MAIN LOOP
        # ================================================================
        while not done:
            if nst > self.max_steps:
                return self._pack(t_out, y_out, False,
                    f"Max steps ({self.max_steps}) at t={t}", nst, nfl, nfe, nje)

            hmin = 16.0 * eps * abs(t) if self.min_step == 0.0 else self.min_step
            absh = min(max(hmin, absh), hmax)
            if abs(absh - hmin) < 100 * eps:
                if at_hmin: absh = abshlast
                at_hmin = True
            else:
                at_hmin = False
            h = tdir * absh

            if 1.1 * absh >= abs(tfinal - t):
                h = tfinal - t; absh = abs(h); done = True

            # Adjust dif if h or k changed
            if (abs(absh - abshlast) / max(absh, _TINY) > 1e-6) or (k != klast):
                dif = _adjust_dif(dif, absh / abshlast, k)
                hinvGak = h * _INVGA[k - 1]
                nconhk = 0
                M_inv = self._get_M_inv(hinvGak, J)
                havrate = False
                mx.eval(*dif[:k])

            # ========== INNER: advance one step ==========
            nofailed = True
            while True:
                gotynew = False
                while not gotynew:
                    # Predictor and psi
                    invga_k = _INVGA[k - 1]
                    psi = _G[0] * invga_k * dif[0]
                    pred = y + dif[0]
                    for j in range(1, k):
                        psi = psi + _G[j] * invga_k * dif[j]
                        pred = pred + dif[j]

                    tnew = tfinal if done else t + h
                    h = tnew - t
                    ynew = pred

                    wt_new = mx.maximum(mx.maximum(mx.abs(ynew), mx.abs(y)), thresh)
                    invwt = 1.0 / wt_new

                    # Newton
                    difkp1 = mx.zeros((N_k, n_var))
                    tooslow = False; oldnrm = 0.0; rate = 0.0

                    for it in range(1, _MAXIT + 1):
                        rhs_v = hinvGak * self.fun(tnew, ynew) - (psi + difkp1)
                        delta = (M_inv @ rhs_v[:, :, None])[:, :, 0]
                        nrm_arr = mx.max(mx.abs(delta) * invwt)
                        mx.eval(delta, nrm_arr)  # single eval per Newton iter
                        nfe += 1
                        newnrm = float(nrm_arr)

                        difkp1 = difkp1 + delta
                        ynew = pred + difkp1

                        minnrm = 100.0 * eps  # simplified floor

                        if newnrm <= minnrm:
                            gotynew = True; break
                        elif it == 1:
                            if havrate and rate < 1.0:
                                if newnrm * rate / (1.0 - rate) <= 0.05 * rtol:
                                    gotynew = True; break
                            rate = 0.0
                        elif newnrm > 0.9 * oldnrm:
                            tooslow = True; break
                        else:
                            rate = max(0.9 * rate, newnrm / oldnrm)
                            havrate = True
                            errit = newnrm * rate / (1 - rate) if rate < 1 else 1e30
                            if errit <= 0.5 * rtol:
                                gotynew = True; break
                            elif it == _MAXIT:
                                tooslow = True; break
                            elif rate > 0 and errit * rate ** (_MAXIT - it) > 0.5 * rtol:
                                tooslow = True; break
                        oldnrm = newnrm

                    if tooslow:
                        nfl += 1
                        if not Jcurrent:
                            f0 = self.fun(t, y); mx.eval(f0); nfe += 1
                            J, nj = self._get_jacobian(t, y, f0); nje += 1; nfe += nj
                            Jcurrent = True
                        elif absh <= hmin:
                            return self._pack(t_out, y_out, False,
                                f"h={absh} too small at t={t}", nst, nfl, nfe, nje)
                        else:
                            abshlast = absh
                            absh = max(0.3 * absh, hmin); h = tdir * absh; done = False
                            dif = _adjust_dif(dif, absh / abshlast, k)
                            hinvGak = h * _INVGA[k - 1]; nconhk = 0
                        M_inv = self._get_M_inv(hinvGak, J)
                        havrate = False
                        mx.eval(*dif[:k])

                # gotynew is True
                mx.eval(ynew, difkp1)

                # Error
                err_raw = mx.max(mx.abs(difkp1) * invwt)
                mx.eval(err_raw)
                err = float(err_raw) * abs(_ERCONST[k - 1])

                if err > rtol:
                    nfl += 1
                    if absh <= hmin:
                        return self._pack(t_out, y_out, False,
                            f"h={absh} too small at t={t}", nst, nfl, nfe, nje)
                    abshlast = absh
                    if nofailed:
                        nofailed = False
                        hopt = absh * max(0.1, 0.833 * (rtol / err) ** (1.0 / (k + 1)))
                        if k > 1:
                            ek1 = float(mx.max(mx.abs(
                                (dif[k-1] + difkp1) * invwt))) * abs(_ERCONST[k-2])
                            hk1 = absh * max(0.1, 0.769 * (rtol / ek1) ** (1.0 / k))
                            if hk1 > hopt:
                                hopt = min(absh, hk1); k -= 1
                        absh = max(hmin, hopt)
                    else:
                        absh = max(hmin, 0.5 * absh)
                    h = tdir * absh
                    if absh < abshlast: done = False
                    dif = _adjust_dif(dif, absh / abshlast, k)
                    hinvGak = h * _INVGA[k - 1]; nconhk = 0
                    M_inv = self._get_M_inv(hinvGak, J)
                    havrate = False
                    mx.eval(*dif[:k])
                else:
                    break  # step accepted

            nst += 1

            # ========== Update backward-difference table ==========
            dif[k + 1] = difkp1 - dif[k]
            dif[k] = difkp1
            for j in range(k - 1, -1, -1):
                dif[j] = dif[j] + dif[j + 1]
            mx.eval(*dif[:k+2])

            # ========== Output ==========
            if tev is not None:
                while nxev < len(tev) and tdir * (tnew - tev[nxev]) >= 0:
                    if abs(tnew - tev[nxev]) < eps * abs(tnew):
                        t_out.append(tnew); y_out.append(ynew)
                    else:
                        yi = _interp_from_dif(tev[nxev], tnew, ynew, h, dif, k)
                        mx.eval(yi)
                        t_out.append(float(tev[nxev])); y_out.append(yi)
                    nxev += 1
            else:
                t_out.append(tnew); y_out.append(ynew)

            if done: break

            # ========== Step / order selection ==========
            klast = k; abshlast = absh
            nconhk = min(nconhk + 1, _MAXK + 2)

            if nconhk >= k + 2:
                temp = 1.2 * (err / rtol) ** (1.0 / (k + 1))
                hopt = absh / temp if temp > 0.1 else 10.0 * absh
                kopt = k

                if k > 1:
                    ek = float(mx.max(mx.abs(dif[k-1]) * invwt)) * abs(_ERCONST[k-2])
                    temp = 1.3 * (ek / rtol) ** (1.0 / k)
                    hk = absh / temp if temp > 0.1 else 10.0 * absh
                    if hk > hopt: hopt = hk; kopt = k - 1

                if k < _MAXK:
                    ek = float(mx.max(mx.abs(dif[k+1]) * invwt)) * abs(_ERCONST[k])
                    temp = 1.4 * (ek / rtol) ** (1.0 / (k + 2))
                    hk = absh / temp if temp > 0.1 else 10.0 * absh
                    if hk > hopt: hopt = hk; kopt = k + 1

                if hopt > absh: absh = hopt; k = kopt

            t = tnew; y = ynew; Jcurrent = False

        return self._pack(t_out, y_out, True, "OK", nst, nfl, nfe, nje)

    def _pack(self, t_out, y_out, ok, msg, nst, nfl, nfe, nje):
        if not t_out:
            return BatchedNDF15Result(
                mx.array([]), mx.zeros((0, self.N_k, self.n_var)),
                ok, msg, nst, nfl, nfe, nje)
        t_arr = mx.array(t_out)
        y_arr = mx.stack(y_out, axis=0)
        mx.eval(t_arr, y_arr)
        return BatchedNDF15Result(t_arr, y_arr, ok, msg, nst, nfl, nfe, nje)


# ---------------------------------------------------------------------------
# Convenience functional API
# ---------------------------------------------------------------------------
def solve_ndf15_batched(fun, t_span, y0, rtol=1e-6, atol=1e-8,
                        jac=None, t_eval=None, max_steps=100_000,
                        min_step=0.0):
    """Solve N_k independent ODE systems simultaneously.

    Parameters
    ----------
    fun : callable(t, y) -> dy
        t is float, y is mx.array (N_k, n_var). Returns (N_k, n_var).
    t_span : (t0, tfinal)
    y0 : mx.array (N_k, n_var)
    rtol, atol : float
    jac : callable(t, y) -> J (N_k, n_var, n_var), optional
    t_eval : array-like, optional
    max_steps, min_step : int, float

    Returns
    -------
    BatchedNDF15Result
    """
    solver = BatchedNDF15(fun, t_span, y0, rtol=rtol, atol=atol,
                          jac=jac, max_steps=max_steps, min_step=min_step)
    return solver.solve(t_eval=t_eval)


# ===========================================================================
# Tests
# ===========================================================================
if __name__ == '__main__':
    import time as _time

    # ==================================================================
    print("=" * 70)
    print("Test 1: Batched stiff exponential decay y' = -lambda * y")
    print("=" * 70)

    N_k = 200
    lam_arr = mx.linspace(1.0, 1000.0, N_k)

    def decay_fun(t, y):
        return -lam_arr[:, None] * y

    def decay_jac(t, y):
        J = mx.zeros((N_k, 1, 1))
        J = J.at[:, 0, 0].add(-lam_arr)
        return J

    y0_d = mx.ones((N_k, 1))

    t0 = _time.perf_counter()
    res = solve_ndf15_batched(decay_fun, (0.0, 0.01), y0_d,
                              rtol=1e-4, atol=1e-10, jac=decay_jac,
                              t_eval=[0.01])
    el = _time.perf_counter() - t0
    exact_d = mx.exp(-lam_arr * 0.01)
    mx.eval(exact_d)
    err_d = float(mx.max(mx.abs(res.y[-1, :, 0] - exact_d)))
    print(f"  N_k={N_k}, Time: {el:.3f}s, Steps: {res.n_steps}, "
          f"Err: {err_d:.2e}, OK: {res.success}")

    # ==================================================================
    print("\n" + "=" * 70)
    print("Test 2: Batched harmonic oscillators y'' + k^2 y = 0")
    print("=" * 70)

    N_k = 100
    k_arr = mx.linspace(0.1, 10.0, N_k)

    def harm_fun(t, y):
        return mx.stack([y[:, 1], -(k_arr**2) * y[:, 0]], axis=1)

    def harm_jac(t, y):
        N = y.shape[0]
        J = mx.zeros((N, 2, 2))
        J = J.at[:, 0, 1].add(mx.ones(N))
        J = J.at[:, 1, 0].add(-k_arr**2)
        return J

    y0_h = mx.zeros((N_k, 2))
    y0_h = y0_h.at[:, 0].add(1.0)

    t0 = _time.perf_counter()
    res = solve_ndf15_batched(harm_fun, (0.0, 10.0), y0_h,
                              rtol=1e-4, atol=1e-8, jac=harm_jac,
                              t_eval=[10.0])
    el = _time.perf_counter() - t0
    exact_p = mx.cos(k_arr * 10.0)
    mx.eval(exact_p)
    err_p = float(mx.max(mx.abs(res.y[-1, :, 0] - exact_p)))
    print(f"  N_k={N_k}, Time: {el:.3f}s, Steps: {res.n_steps}, "
          f"Err: {err_p:.2e}, OK: {res.success}")

    # With numerical Jacobian
    t0 = _time.perf_counter()
    res2 = solve_ndf15_batched(harm_fun, (0.0, 10.0), y0_h,
                               rtol=1e-4, atol=1e-8,
                               t_eval=[10.0])
    el2 = _time.perf_counter() - t0
    err2 = float(mx.max(mx.abs(res2.y[-1, :, 0] - exact_p)))
    print(f"  Numerical Jac: Time: {el2:.3f}s, Steps: {res2.n_steps}, Err: {err2:.2e}")

    # ==================================================================
    print("\n" + "=" * 70)
    print("Test 3: Scale test -- 500 modes x 10 vars")
    print("=" * 70)

    N_k = 500; n_var = 10
    np.random.seed(42)
    A_np = np.random.randn(N_k, n_var, n_var).astype(np.float32) * 0.1
    for i in range(N_k):
        A_np[i] -= np.diag(np.abs(A_np[i]).sum(axis=1) +
                           np.random.uniform(1, 100, n_var).astype(np.float32))
    A_mx = mx.array(A_np)

    def lin_fun(t, y):
        return (A_mx @ y[:, :, None])[:, :, 0]

    def lin_jac(t, y):
        return A_mx

    y0_l = mx.ones((N_k, n_var))

    t0 = _time.perf_counter()
    res = solve_ndf15_batched(lin_fun, (0.0, 1.0), y0_l,
                              rtol=1e-4, atol=1e-8, jac=lin_jac,
                              t_eval=[1.0])
    el = _time.perf_counter() - t0
    print(f"  N_k={N_k}, n_var={n_var}, Time: {el:.3f}s")
    print(f"  Steps: {res.n_steps}, Failed: {res.n_failed}, "
          f"Fevals: {res.n_fevals}, Jevals: {res.n_jevals}")
    print(f"  OK: {res.success}, y shape: {res.y.shape}")

    try:
        from scipy.linalg import expm
        errors = []
        for ik in [0, N_k//2, N_k-1]:
            ye = expm(A_np[ik].astype(np.float64)) @ np.ones(n_var)
            yc = np.array(res.y[-1, ik, :].tolist())
            errors.append(float(np.max(np.abs(yc - ye) / (np.abs(ye) + 1e-10))))
        print(f"  Rel errors (modes 0,{N_k//2},{N_k-1}): "
              f"{errors[0]:.2e}, {errors[1]:.2e}, {errors[2]:.2e}")
    except ImportError:
        print("  (scipy not available)")

    print("\n" + "=" * 70)
    print("All tests complete.")
    print("=" * 70)
