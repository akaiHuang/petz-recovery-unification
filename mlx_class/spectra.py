"""
spectra.py — C_l angular power spectrum (MLX GPU).

C_l = (4pi) int d(ln k) P_R(k) |Delta_l(k)|^2

where Delta_l = SW * j_l(k D_A) + Doppler * j_l'(k D_A)

Author: Sheng-Kai Huang, 2026
"""
import numpy as np
import mlx.core as mx
from scipy.special import spherical_jn
from .background import A_s, n_s, k_pivot, T_CMB


def compute_cl(source_SW, k_arr_Mpc, ell_values, D_A, source_Dop=None):
    """
    Compute C_l from source function(s) and Bessel projection.

    Parameters
    ----------
    source_SW : array (N_k,)
        Sachs-Wolfe source: (Theta_0 + Phi) * silk_damping
    k_arr_Mpc : array (N_k,)
    ell_values : array (N_ell,)
    D_A : float
        Comoving distance to LSS (Mpc).
    source_Dop : array (N_k,) or None
        Doppler source: v_b * silk_damping. If None, Doppler term omitted.

    Returns
    -------
    ell_values, Cl, Dl (in muK^2)
    """
    import time

    N_ell = len(ell_values)
    N_k = len(k_arr_Mpc)
    x_arr = k_arr_Mpc * D_A

    # Primordial power spectrum
    P_R = A_s * (k_arr_Mpc / k_pivot)**(n_s - 1.0)

    # Bessel table
    t0 = time.time()
    jl = np.zeros((N_ell, N_k), dtype=np.float32)
    jlp = np.zeros((N_ell, N_k), dtype=np.float32) if source_Dop is not None else None
    for il, ell in enumerate(ell_values):
        jl[il] = spherical_jn(int(ell), x_arr)
        if jlp is not None:
            jlp[il] = spherical_jn(int(ell), x_arr, derivative=True)
    print(f"[Spectra] Bessel: {time.time()-t0:.2f}s ({N_ell} ells)")

    # Transfer function: Delta_l(k) = SW * jl + Dop * jl'
    Delta_l = np.array(source_SW[None, :]) * jl
    if source_Dop is not None and jlp is not None:
        Delta_l += np.array(source_Dop[None, :]) * jlp

    # C_l = 4pi int dk/k P_R |Delta_l|^2
    t0 = time.time()
    integrand = (P_R[None, :] * Delta_l**2).astype(np.float32)
    lnk = np.log(k_arr_Mpc)
    dlnk = np.diff(lnk).astype(np.float32)

    integrand_gpu = mx.array(integrand)
    dlnk_gpu = mx.array(dlnk)
    mid = 0.5 * (integrand_gpu[:, :-1] + integrand_gpu[:, 1:])
    Cl_gpu = 4.0 * np.pi * mx.sum(mid * dlnk_gpu[None, :], axis=1)
    mx.eval(Cl_gpu)
    print(f"[Spectra] GPU integration: {time.time()-t0:.4f}s")

    Cl = np.maximum(np.array(Cl_gpu), 0.0)
    ell_f = ell_values.astype(float)
    Dl = ell_f * (ell_f + 1.0) * Cl / (2.0 * np.pi) * (T_CMB * 1e6)**2

    return ell_values, Cl, Dl
