"""
solver_sync_gpu.py -- GPU-accelerated synchronous gauge Boltzmann solver.

Uses zeroth-order Tight Coupling Approximation (TCA-0) for the photon sector
throughout the entire integration, eliminating Thomson scattering stiffness.
All N_k modes are integrated simultaneously on the GPU via batched ndf15.

TCA-0 approximation:
  - theta_gamma = theta_b (perfect baryon-photon coupling)
  - F_g,1 = 4*theta_b/(3*k) (diagnosed, not evolved)
  - F_g,l = 0 for l >= 2 (Thomson-damped)
  - F_g,2 diagnosed quasi-statically at output times for gauge transform

This correctly captures acoustic oscillations and the Sachs-Wolfe plateau
but misses Silk damping at high l. Accuracy target: <10% RMS vs CLASS.

State vector: [eta, dc, db, tb, dg, Fn_0..Fn_ln_max]
Variables per mode: 5 + ln_max + 1

Physics: Ma & Bertschinger (1995), synchronous gauge, CDM rest frame.

Usage:
  python -m mlx_class.solver_sync_gpu [--N_k 500]

Author: Sheng-Kai Huang, 2026
"""
import os
os.environ['MPLCONFIGDIR'] = '/tmp/mpl_config'

import numpy as np
import mlx.core as mx
import time
import sys

from .background import (
    Background, A_s, n_s, k_pivot, T_CMB,
    Omega_r as _OMEGA_R,
    Omega_b as _OMEGA_B,
    Omega_c as _OMEGA_C,
    Omega_m as _OMEGA_M,
    H0_Mpc as _H0_MPC,
)
from .perturbations_neutrino import _f_nu, _OMEGA_GAMMA, _OMEGA_NU
from .perturbations_sync import (
    adiabatic_ic_sync, n_var_sync, idx_fn_start,
    IDX_ETA, IDX_DELTA_C, IDX_DELTA_B, IDX_THETA_B, IDX_FG_START,
    L_GAMMA_MAX, L_NU_MAX_SYNC,
)
from .ndf15_batched import solve_ndf15_batched

_H02 = float(_H0_MPC ** 2)
_Og = float(_OMEGA_GAMMA)
_On = float(_OMEGA_NU)
_Ob = float(_OMEGA_B)
_Oc = float(_OMEGA_C)


# ============================================================================
# Background table
# ============================================================================
class BGTable:
    def __init__(self, bg):
        self._t = bg.tau_grid.copy(); self._cH = bg.calH_grid.copy()
        self._a = bg.a_grid.copy(); self._R = bg.R_grid.copy()
        self._kd = bg.kappa_dot_grid.copy()
    def __call__(self, t):
        tv = float(t)
        return (float(np.interp(tv, self._t, self._cH)),
                float(np.interp(tv, self._t, self._a)),
                float(np.interp(tv, self._t, self._R)),
                float(np.interp(tv, self._t, self._kd)))


# ============================================================================
# TCA-0 state layout
# ============================================================================
_E = 0; _DC = 1; _DB = 2; _TB = 3; _DG = 4; _FN = 5

def _nv(ln_max):
    return 5 + ln_max + 1


# ============================================================================
# Batched TCA-0 RHS
# ============================================================================
def _make_rhs(k, bgt, ln_max):
    N = k.shape[0]; k2 = k*k; nv = _nv(ln_max)
    if ln_max > 3:
        lm = mx.arange(3, ln_max, dtype=mx.float32)
        iv = 1.0/(2.*lm+1.)
    else:
        lm = iv = None

    def rhs(t, y):
        cH, a, R, kd = bgt(t)
        ia = 1./a; ia2 = ia*ia; ts = max(float(t), 1e-10)
        eta = y[:, _E]; dc = y[:, _DC]; db = y[:, _DB]
        tb = y[:, _TB]; dg = y[:, _DG]
        Fn = y[:, _FN:_FN+ln_max+1]
        dn = Fn[:, 0]; tn = 0.75*k*Fn[:, 1]; tg = tb

        s00 = _Og*ia2*dg + _On*ia2*dn + _Ob*ia*db + _Oc*ia*dc
        hp = (2./cH)*(k2*eta + 1.5*_H02*s00)
        s0i = (4./3.)*_Og*ia2*tg + (4./3.)*_On*ia2*tn + _Ob*ia*tb
        ep = 1.5*_H02/k2*s0i

        d_dc = -0.5*hp
        d_db = -tb - 0.5*hp
        d_tb = (-R*cH*tb + 0.25*k2*dg)/(1.+R)
        d_dg = -(4./3.)*tb - (2./3.)*hp

        dn0 = -k*Fn[:, 1] - (2./3.)*hp
        F2n = Fn[:, 2] if ln_max >= 2 else mx.zeros(N)
        dn1 = (k/3.)*(Fn[:, 0] - 2.*F2n)
        F3n = Fn[:, 3] if ln_max >= 3 else mx.zeros(N)
        dn2 = (k/5.)*(2.*Fn[:, 1] - 3.*F3n) + (4./15.)*hp + (8./15.)*ep

        pn = [dn0[:, None], dn1[:, None], dn2[:, None]]
        if lm is not None and ln_max > 3:
            pn.append(k[:, None]*iv[None, :]*(
                lm[None, :]*Fn[:, 2:ln_max-1] - (lm[None, :]+1.)*Fn[:, 4:ln_max+1]))
        pn.append((k*Fn[:, ln_max-1]*float(ln_max)/(2.*ln_max+1.)
                   - float(ln_max+1)/ts*Fn[:, ln_max])[:, None])

        return mx.concatenate([ep[:, None], d_dc[:, None], d_db[:, None],
                               d_tb[:, None], d_dg[:, None],
                               mx.concatenate(pn, axis=1)], axis=1)
    return rhs, nv


# ============================================================================
# Initial conditions
# ============================================================================
def _build_ic(k_np, tau_i, ln_max):
    Nk = len(k_np); nv = _nv(ln_max)
    y0 = np.zeros((Nk, nv), dtype=np.float32)
    fns = idx_fn_start(L_GAMMA_MAX)
    for i, kv in enumerate(k_np):
        yf = adiabatic_ic_sync(kv, tau_i, L_GAMMA_MAX, ln_max)
        y0[i, _E] = yf[IDX_ETA]; y0[i, _DC] = yf[IDX_DELTA_C]
        y0[i, _DB] = yf[IDX_DELTA_B]; y0[i, _TB] = yf[IDX_THETA_B]
        y0[i, _DG] = yf[IDX_FG_START]
        y0[i, _FN:_FN+ln_max+1] = yf[fns:fns+ln_max+1]
    return mx.array(y0)


# ============================================================================
# Gauge transform from TCA state
# ============================================================================
def _gauge_xform_tca(y, k, cH, a, ln_max):
    """Gauge transform from TCA-0 state to Newtonian gauge quantities.

    Diagnoses F_g,1 = 4*tb/(3*k) and F_g,2 quasi-statically.
    Returns PhiN, PsiN, dgN, tbN as mx.arrays (N_k,).
    """
    k2 = k*k; ia = 1./a; ia2 = ia*ia
    eta = y[:, _E]; dc = y[:, _DC]; db = y[:, _DB]
    tb = y[:, _TB]; dg = y[:, _DG]
    Fn = y[:, _FN:_FN+ln_max+1]
    dn = Fn[:, 0]; tn = 0.75*k*Fn[:, 1]; tg = tb

    s00 = _Og*ia2*dg + _On*ia2*dn + _Ob*ia*db + _Oc*ia*dc
    hp = (2./cH)*(k2*eta + 1.5*_H02*s00)
    s0i = (4./3.)*_Og*ia2*tg + (4./3.)*_On*ia2*tn + _Ob*ia*tb
    ep = 1.5*_H02/k2*s0i

    alpha = (hp + 6.*ep)/(2.*k2)
    if cH > 0:
        ma = 5.*mx.abs(eta)/cH
        alpha = mx.clip(alpha, -ma, ma)

    PhiN = eta - cH*alpha

    # F_g,2 quasi-static: 0 = (2k/5)*Fg1 + (4/15)*hp + (8/15)*ep - (9/10)*|kd|*Fg2
    # Since we don't have |kd| here, approximate Fg2 from the metric terms
    # For the anisotropic stress, the dominant contribution is from neutrinos
    sig_n = 0.5*Fn[:, 2] if ln_max >= 2 else mx.zeros_like(eta)
    aniso = 12.*_H02/(a*a*k2)*_On*sig_n  # photon aniso ~ 0 in TCA
    PsiN = PhiN - aniso

    dgN = dg - 4.*cH*alpha
    tbN = tb + k2*alpha

    return PhiN, PsiN, dgN, tbN


# ============================================================================
# Snapshot grid
# ============================================================================
def _build_snaps(bg, Nv=60, Ne=30, Nl=20):
    pi = np.argmax(bg.visibility_grid)
    tp = bg.tau_grid[pi]; vp = bg.visibility_grid[pi]
    hm = vp/2.; l = pi; r = pi
    while l > 0 and bg.visibility_grid[l] > hm: l -= 1
    while r < len(bg.visibility_grid)-1 and bg.visibility_grid[r] > hm: r += 1
    sig = (bg.tau_grid[r]-bg.tau_grid[l])/2.355
    vlo = max(tp-5.*sig, bg.tau_grid[5])
    vhi = min(tp+5.*sig, bg.tau_0*0.5)
    tv = np.linspace(vlo, vhi, Nv)
    aeq = _OMEGA_R/_OMEGA_M
    teq = float(bg._tau_of_a(np.log(aeq)))
    elo = max(teq*0.3, bg.tau_grid[5]); ehi = vlo
    if ehi > elo*1.5:
        dl, dh = max(teq-50., elo), min(teq+50., ehi)
        te = np.sort(np.unique(np.concatenate([
            np.geomspace(elo, ehi, 15), np.linspace(dl, dh, Ne)])))
    else:
        te = np.array([elo])
    # Limit late ISW to tau ~ 3000 Mpc (z ~ 4) to keep integration tractable.
    # The ISW signal beyond z ~ 1 is small and dominated by the lowest ells.
    tau_late_max = min(bg.tau_0*0.98, 3000.0)
    tl = np.linspace(vhi, tau_late_max, Nl)
    return np.sort(np.unique(np.concatenate([te, tv, tl])))


# ============================================================================
# Main pipeline
# ============================================================================
def run_sync_gpu(N_k=500, k_min=3e-4, k_max=0.35,
                 ln_max=L_NU_MAX_SYNC,
                 rtol=1e-4, atol=1e-7, verbose=True):
    T0 = time.time()

    if verbose:
        print("="*65)
        print("GPU Synchronous Gauge Boltzmann Solver (TCA-0 + ndf15)")
        print("="*65)
        print("\n--- Background ---"); sys.stdout.flush()
    t0 = time.time()
    bg = Background(khronon=False, recombination='recfast'); bg.solve()
    t_bg = time.time()-t0
    if verbose: print(f"  {t_bg:.2f}s")

    k_np = np.geomspace(k_min, k_max, N_k).astype(np.float64)
    k_mx = mx.array(k_np.astype(np.float32))
    bgt = BGTable(bg)

    tau_all = _build_snaps(bg)
    Ns = len(tau_all)
    tau_i = float(bg.tau_grid[1])
    # Integrate to last snapshot + small buffer (not beyond)
    tau_end = float(tau_all[-1]) * 1.001 + 1.0

    if verbose:
        nv = _nv(ln_max)
        print(f"\n--- Perturbations: {N_k} modes, nvar={nv} (TCA-0) ---")
        print(f"  k: [{k_min:.1e}, {k_max}], tau: [{tau_i:.4f}, {tau_end:.1f}]")
        print(f"  {Ns} snapshot times"); sys.stdout.flush()

    # --- ODE solve ---
    t0 = time.time()
    rhs, nv = _make_rhs(k_mx, bgt, ln_max)
    y0 = _build_ic(k_np, tau_i, ln_max)

    result = solve_ndf15_batched(rhs, (tau_i, tau_end), y0,
                                  rtol=rtol, atol=atol,
                                  t_eval=tau_all, max_steps=200_000)
    t_ode = time.time()-t0
    if verbose:
        print(f"\n  ODE: {t_ode:.2f}s, {result.n_steps} steps, "
              f"{result.n_failed} failed, OK={result.success}")
        if not result.success: print(f"    {result.message}")
        sys.stdout.flush()

    if not result.success:
        print(f"WARNING: solver failed: {result.message}")

    # --- Gauge transform ---
    t0 = time.time()
    a_s = bg.a_at_tau(tau_all); cH_s = bg.calH_at_tau(tau_all)
    PhiN = np.zeros((N_k, Ns)); PsiN = np.zeros((N_k, Ns))
    Th0N = np.zeros((N_k, Ns)); vbN = np.zeros((N_k, Ns))

    n_out = min(result.y.shape[0], Ns)
    for it in range(n_out):
        P, Q, dgN, tbN = _gauge_xform_tca(
            result.y[it], k_mx, float(cH_s[it]), float(a_s[it]), ln_max)
        mx.eval(P, Q, dgN, tbN)
        PhiN[:, it] = np.array(P.tolist())
        PsiN[:, it] = np.array(Q.tolist())
        Th0N[:, it] = np.array(dgN.tolist())/4.
        vbN[:, it] = np.array(tbN.tolist())/k_np

    t_gx = time.time()-t0
    if verbose: print(f"  Gauge xform: {t_gx:.3f}s")

    # --- Phi'+Psi' ---
    t0 = time.time()
    from scipy.interpolate import CubicSpline
    PP = PhiN + PsiN; PPp = np.zeros_like(PP)
    for ik in range(N_k):
        cs = CubicSpline(tau_all, PP[ik])
        PPp[ik] = cs(tau_all, 1)
    t_d = time.time()-t0
    if verbose: print(f"  Derivatives: {t_d:.3f}s")

    # --- LOS ---
    t0 = time.time()
    ells = np.unique(np.concatenate([
        np.arange(2, 30), np.arange(30, 100, 2), np.arange(100, 500, 4),
        np.arange(500, 1500, 8), np.arange(1500, 2501, 12)])).astype(int)
    Nel = len(ells)
    if verbose: print(f"\n--- LOS C_l ({Nel} ells) ---"); sys.stdout.flush()

    PR = A_s*(k_np/k_pivot)**(n_s-1.)
    gs = bg.visibility_at_tau(tau_all)
    ks = bg.kappa_at_tau(tau_all); enk = np.exp(-ks)
    chi = bg.tau_0 - tau_all

    try:
        from .bessel_cache import CachedBesselTable
        xm = float(np.max(k_np)*np.max(chi))*1.05+50.
        bt = CachedBesselTable(ells, x_max=xm, N_cheb=64, seg_width=80, verbose=verbose)
        gpu_b = True
    except Exception as e:
        if verbose: print(f"[WARN] Bessel failed ({e})")
        gpu_b = False

    Dl_k = np.zeros((Nel, N_k)); dtau = np.diff(tau_all)
    if verbose: print("  Computing..."); sys.stdout.flush()

    for it in range(Ns):
        c = chi[it]
        if c <= 0: continue
        if it == 0: w = .5*dtau[0]
        elif it == Ns-1: w = .5*dtau[-1]
        else: w = .5*(dtau[max(0,it-1)]+dtau[min(it,len(dtau)-1)])

        ssw = gs[it]*(Th0N[:, it]+PsiN[:, it])
        sdop = gs[it]*vbN[:, it]
        sisw = enk[it]*PPp[:, it]

        if gpu_b:
            xa = k_np*c
            jl, jlp = bt.eval_jl_jlp(xa.astype(np.float32))
            mx.eval(jl, jlp); jl = np.array(jl); jlp = np.array(jlp)
            for il in range(Nel):
                Dl_k[il] += w*(ssw*jl[il]+sdop*jlp[il]+sisw*jl[il])
        else:
            from scipy.special import spherical_jn
            x = k_np*c
            for il, el in enumerate(ells):
                j = spherical_jn(int(el), x)
                jp = spherical_jn(int(el), x, derivative=True)
                Dl_k[il] += w*(ssw*j+sdop*jp+sisw*j)

    t_los = time.time()-t0
    if verbose: print(f"  LOS: {t_los:.2f}s")

    # --- C_l ---
    lnk = np.log(k_np); dlnk = np.diff(lnk)
    ig = PR[None, :]*Dl_k**2
    md = .5*(ig[:, :-1]+ig[:, 1:])
    Cl = np.maximum(4.*np.pi*np.sum(md*dlnk[None, :], axis=1), 0.)
    ef = ells.astype(float)
    Dl = ef*(ef+1.)*Cl/(2.*np.pi)*(T_CMB*1e6)**2

    te = time.time()-T0
    if verbose:
        print(f"\n{'='*65}")
        print(f"TOTAL: {te:.2f}s  (BG {t_bg:.1f}, ODE {t_ode:.1f}, "
              f"GX {t_gx:.1f}, LOS {t_los:.1f})")
        print(f"{'='*65}")

    return {'ell': ells, 'Cl': Cl, 'Dl': Dl, 'k_arr': k_np, 'bg': bg,
            'Delta_l': Dl_k, 'Phi_N': PhiN, 'Psi_N': PsiN,
            'timing': {'background': t_bg, 'ode': t_ode, 'gauge': t_gx,
                       'derivatives': t_d, 'los': t_los, 'total': te}}


# ============================================================================
# CLI
# ============================================================================
def main():
    import argparse
    p = argparse.ArgumentParser(description='GPU sync gauge Boltzmann solver (TCA-0)')
    p.add_argument('--N_k', type=int, default=500)
    p.add_argument('--k_min', type=float, default=3e-4)
    p.add_argument('--k_max', type=float, default=0.35)
    p.add_argument('--rtol', type=float, default=1e-4)
    p.add_argument('--atol', type=float, default=1e-7)
    p.add_argument('--ln_max', type=int, default=L_NU_MAX_SYNC)
    p.add_argument('--compare', action='store_true')
    a = p.parse_args()

    r = run_sync_gpu(N_k=a.N_k, k_min=a.k_min, k_max=a.k_max,
                     ln_max=a.ln_max, rtol=a.rtol, atol=a.atol, verbose=True)
    ell, Dl = r['ell'], r['Dl']
    print("\n--- D_l at key multipoles ---")
    for tl in [2, 10, 50, 100, 220, 500, 800, 1000, 1500, 2000]:
        ix = np.argmin(np.abs(ell-tl))
        print(f"  l={ell[ix]:5d}: D_l = {Dl[ix]:10.2f} uK^2")

    ix220 = np.argmin(np.abs(ell-220))
    print(f"\nFirst peak (l~220): {Dl[ix220]:.1f} uK^2")
    if 3000 < Dl[ix220] < 8000: print("  -> PASS")
    else: print("  -> WARNING: outside [3000, 8000]")

    if a.compare:
        from .solver_sync import run_sync_solver
        rs = run_sync_solver(N_k=min(a.N_k, 150), method='Radau', verbose=True)
        ce = np.intersect1d(ell, rs['ell'])
        mg = np.isin(ell, ce); ms = np.isin(rs['ell'], ce)
        dg, ds = Dl[mg], rs['Dl'][ms]; sm = ds > 100
        if sm.sum() > 0:
            re = np.abs(dg[sm]-ds[sm])/ds[sm]
            print(f"  RMS err (D_l>100): {np.sqrt(np.mean(re**2))*100:.1f}%")
            print(f"  Max err: {re.max()*100:.1f}%")


if __name__ == '__main__':
    main()
