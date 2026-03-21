"""
transfer.py — Bessel function table for C_l projection.

j_l(k * D_A) computed with scipy (exact) on CPU.
This is typically the bottleneck (~0.2-0.8s for 372 ells).

Author: Sheng-Kai Huang, 2026
"""
import numpy as np
from scipy.special import spherical_jn


def compute_bessel_table(k_arr_Mpc, D_A, ell_values):
    """
    Compute j_l(k * D_A)^2 for all (ell, k) pairs.

    Returns mx.array shape (N_ell, N_k).
    """
    import mlx.core as mx
    x_arr = k_arr_Mpc * D_A
    N_ell = len(ell_values)
    N_k = len(k_arr_Mpc)

    jl2 = np.zeros((N_ell, N_k), dtype=np.float32)
    for il, ell in enumerate(ell_values):
        jl2[il] = spherical_jn(int(ell), x_arr)**2

    return mx.array(jl2)
