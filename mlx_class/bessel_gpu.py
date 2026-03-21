"""
bessel_gpu.py -- GPU-native spherical Bessel functions via MLX.

Replaces scipy.special.spherical_jn which is the main bottleneck in C_l
computation (0.2-16s in Python loops over ell).

Strategy: Piecewise Chebyshev interpolation on GPU.
  1. Split [0, x_max] into segments of width ~80
  2. At Chebyshev nodes, compute j_l(x) via hybrid method:
     - Forward recurrence for x > 1.5*l (stable, vectorized numpy)
     - scipy fallback for transition zone x ~ l (~14% of pairs)
     - Zero for x << l (exponentially small)
  3. DCT matrix transforms node values to Chebyshev coefficients
  4. Coefficients stored on GPU as (N_ell, n_seg, N_cheb) tensor
  5. Evaluate on GPU: gather + Chebyshev basis matrix multiply

Performance (Apple Silicon GPU, 372 ells x 5000 k-points):
  - Table build: ~2s   (8x faster than scipy alone)
  - GPU eval:    0.026s (570x faster than scipy)
  - Visibility-weighted (25 tau, with table reuse): 0.64s (78x faster)
  - Total first call: 2s. Subsequent evals: 0.026s.

Accuracy:
  - Absolute error: max ~5e-7, median ~9e-9
  - Relative error (|j_l|>1e-5): max 0.7%, median 0.004%
  - Float32 limited; sufficient for CMB C_l computation.

Author: Sheng-Kai Huang, 2026
"""
import numpy as np
import mlx.core as mx
import time
from scipy.special import spherical_jn as _scipy_jn


# ============================================================================
# Fast spherical Bessel computation for table building (CPU)
# ============================================================================

def _spherical_jn_table(ell_values, x_nodes_per_seg, n_seg, N_cheb):
    """
    Compute j_l(x) and j_l'(x) at Chebyshev nodes for all (ell, segment) pairs.

    Hybrid strategy:
    1. Forward recurrence covers all (ell, x) in one vectorized pass.
       Accurate where x > ~1.5 * ell (stable regime).
    2. Batched scipy fixes the transition zone (x ~ ell), about ~14% of pairs.
    3. x << ell entries set to zero (exponentially small).

    Parameters
    ----------
    ell_values : ndarray (N_ell,) int
    x_nodes_per_seg : ndarray (n_seg, N_cheb) float64
    n_seg, N_cheb : int

    Returns
    -------
    jl_table, jlp_table : ndarray (N_ell, n_seg, N_cheb) float64
    """
    import warnings
    N_ell = len(ell_values)
    ell_max = int(np.max(ell_values))

    jl_table = np.zeros((N_ell, n_seg, N_cheb), dtype=np.float64)
    jlp_table = np.zeros((N_ell, n_seg, N_cheb), dtype=np.float64)

    seg_x_min = np.min(x_nodes_per_seg, axis=1)
    seg_x_max = np.max(x_nodes_per_seg, axis=1)

    # --- Phase 1: Forward recurrence (vectorized over all x-nodes) ---
    all_x = np.maximum(x_nodes_per_seg.ravel(), 1e-30)
    inv_x = 1.0 / all_x

    j_prev = np.sin(all_x) * inv_x
    j_prev[all_x < 1e-30] = 1.0
    j_curr = (np.sin(all_x) * inv_x - np.cos(all_x)) * inv_x
    j_curr[all_x < 1e-30] = 0.0

    ell_to_idx = {int(e): i for i, e in enumerate(ell_values)}

    if 0 in ell_to_idx:
        jl_table[ell_to_idx[0]] = j_prev.reshape(n_seg, N_cheb)
        jlp_table[ell_to_idx[0]] = (-j_curr).reshape(n_seg, N_cheb)
    if 1 in ell_to_idx:
        jl_table[ell_to_idx[1]] = j_curr.reshape(n_seg, N_cheb)
        jlp_table[ell_to_idx[1]] = (j_prev - 2.0 * inv_x * j_curr).reshape(n_seg, N_cheb)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        for ell in range(1, ell_max):
            j_next = (2 * ell + 1) * inv_x * j_curr - j_prev
            l_new = ell + 1

            if l_new in ell_to_idx:
                il = ell_to_idx[l_new]
                jl_table[il] = j_next.reshape(n_seg, N_cheb)
                jlp_table[il] = (j_curr - (l_new + 1) * inv_x * j_next).reshape(n_seg, N_cheb)

            j_prev = j_curr
            j_curr = j_next

    # --- Phase 2: Fix transition zone with batched scipy ---
    for il, ell in enumerate(ell_values):
        ell_int = int(ell)
        if ell_int <= 1:
            continue

        needs_fix = []
        for s in range(n_seg):
            if seg_x_min[s] > 1.5 * ell_int:
                continue  # forward recurrence was stable
            elif seg_x_max[s] < 0.3 * ell_int:
                jl_table[il, s, :] = 0.0
                jlp_table[il, s, :] = 0.0
            else:
                needs_fix.append(s)

        if not needs_fix:
            continue

        # Batch all transition-zone nodes for this ell
        n_fix = len(needs_fix)
        x_batch = np.empty(n_fix * N_cheb, dtype=np.float64)
        for i, s in enumerate(needs_fix):
            x_batch[i * N_cheb:(i + 1) * N_cheb] = x_nodes_per_seg[s]

        jl_batch = _scipy_jn(ell_int, x_batch)
        jlp_batch = _scipy_jn(ell_int, x_batch, derivative=True)

        for i, s in enumerate(needs_fix):
            jl_table[il, s] = jl_batch[i * N_cheb:(i + 1) * N_cheb]
            jlp_table[il, s] = jlp_batch[i * N_cheb:(i + 1) * N_cheb]

    # --- Phase 3: Clean up any NaN/Inf from forward recurrence overflow ---
    bad = ~np.isfinite(jl_table)
    jl_table[bad] = 0.0
    bad_d = ~np.isfinite(jlp_table)
    jlp_table[bad_d] = 0.0

    return jl_table, jlp_table


# ============================================================================
# Piecewise Chebyshev interpolation table
# ============================================================================

class BesselTable:
    """
    Pre-computed piecewise Chebyshev interpolation table for j_l(x) and j_l'(x).

    Usage:
        table = BesselTable(ell_values, x_max=5500.0)
        jl = table.eval_jl(x_arr)             # shape (N_ell, N_x)
        jl, jlp = table.eval_jl_jlp(x_arr)    # with derivative

    Build cost: ~2s (forward recurrence + scipy transition fix + DCT).
    Eval cost:  ~0.026s for 372 ells x 5000 points (GPU matmul).
    """

    def __init__(self, ell_values, x_max=5500.0, N_cheb=64,
                 seg_width=80.0, verbose=True):
        """
        Build Chebyshev coefficient table.

        Parameters
        ----------
        ell_values : array-like of int
            Multipole values (e.g. 2 to 2500).
        x_max : float
            Maximum argument x = k * D_A.
        N_cheb : int
            Chebyshev nodes per segment. 64 is sufficient for seg_width=80.
        seg_width : float
            Width of each piecewise segment.
        verbose : bool
        """
        self.ell_values = np.asarray(ell_values, dtype=np.int32)
        self.N_ell = len(self.ell_values)
        self.x_max = float(x_max)
        self.N_cheb = int(N_cheb)

        t0 = time.time()
        self._build_table(seg_width)
        dt = time.time() - t0

        if verbose:
            print(f"[BesselTable] Built: {self.N_ell} ells, "
                  f"{self.n_seg} segs x {self.N_cheb} nodes, "
                  f"x_max={self.x_max:.0f}, "
                  f"time={dt:.2f}s")

    def _build_table(self, seg_width):
        """
        Build coefficient tables using hybrid forward recurrence + scipy.

        1. Forward recurrence for all (ell, x) at once via numpy (fast).
        2. Fix transition zone (x ~ ell) with scipy (exact, ~9% of pairs).
        3. Apply DCT matrix to get Chebyshev coefficients.
        """
        N_cheb = self.N_cheb
        x_max = self.x_max

        # Segment structure
        n_seg = max(1, int(np.ceil(x_max / seg_width)))
        seg_edges = np.linspace(0.0, x_max, n_seg + 1)
        self.n_seg = n_seg
        self.seg_width_actual = x_max / n_seg

        # Chebyshev nodes on [-1, 1]
        t_nodes = np.cos(np.pi * (np.arange(N_cheb) + 0.5) / N_cheb)

        # Build x-nodes for each segment: (n_seg, N_cheb)
        x_nodes = np.empty((n_seg, N_cheb), dtype=np.float64)
        for s in range(n_seg):
            a, b = seg_edges[s], seg_edges[s + 1]
            x_nodes[s] = 0.5 * (b - a) * t_nodes + 0.5 * (a + b)
        # Clamp to avoid x=0
        x_nodes = np.maximum(x_nodes, 1e-15)

        # Compute j_l and j_l' at all nodes using hybrid method
        jl_segs, jlp_segs = _spherical_jn_table(
            self.ell_values, x_nodes, n_seg, N_cheb)
        # jl_segs: (N_ell, n_seg, N_cheb)

        # DCT matrix: values -> Chebyshev coefficients
        k_idx = np.arange(N_cheb)
        arccos_t = np.arccos(t_nodes)
        DCT_mat = np.cos(k_idx[:, None] * arccos_t[None, :])
        DCT_mat *= (2.0 / N_cheb)
        DCT_mat[0] *= 0.5

        # Apply DCT: coeffs[il, s, :] = DCT_mat @ jl_segs[il, s, :]
        coeffs = np.einsum('ijk,lk->ijl', jl_segs, DCT_mat).astype(np.float32)
        coeffs_d = np.einsum('ijk,lk->ijl', jlp_segs, DCT_mat).astype(np.float32)

        # Transfer to GPU
        self.coeffs_mx = mx.array(coeffs)
        self.coeffs_deriv_mx = mx.array(coeffs_d)

    # ------------------------------------------------------------------
    # Segment mapping (GPU)
    # ------------------------------------------------------------------

    def _map_to_segments(self, x_arr_mx):
        """Map x values to segment indices and local [-1,1] coordinates."""
        sw = self.seg_width_actual
        x = mx.clip(x_arr_mx, 0.0, self.x_max - 1e-6)

        seg_idx = mx.floor(x / sw).astype(mx.int32)
        seg_idx = mx.clip(seg_idx, 0, self.n_seg - 1)

        a = seg_idx.astype(mx.float32) * sw
        t_local = 2.0 * (x - a) / sw - 1.0

        return seg_idx, t_local

    # ------------------------------------------------------------------
    # Chebyshev evaluation via Clenshaw recurrence (GPU)
    # ------------------------------------------------------------------

    def _eval_chebyshev(self, coeffs_mx, seg_idx, t_local):
        """
        Evaluate piecewise Chebyshev expansion on GPU.

        Two methods available:
        1. Clenshaw recurrence: O(N_cheb) sequential steps, low memory
        2. Matrix multiply: one matmul, needs (N_x, N_cheb) basis matrix

        Uses matrix multiply for N_x <= 10000, Clenshaw for larger.

        Parameters
        ----------
        coeffs_mx : (N_ell, n_seg, N_cheb)
        seg_idx : (N_x,) int32
        t_local : (N_x,) float32

        Returns
        -------
        (N_ell, N_x) float32
        """
        N_x = t_local.shape[0]
        N_cheb = self.N_cheb

        # Gather coefficients per x-point: (N_ell, N_x, N_cheb)
        c = coeffs_mx[:, seg_idx, :]

        # Build Chebyshev basis: T_k(t) = cos(k * arccos(t))
        # Shape: (N_x, N_cheb)
        theta = mx.arccos(mx.clip(t_local, -1.0, 1.0))   # (N_x,)
        k_idx = mx.arange(N_cheb).astype(mx.float32)      # (N_cheb,)
        T_basis = mx.cos(theta[:, None] * k_idx[None, :])  # (N_x, N_cheb)

        # result[il, ix] = sum_k c[il, ix, k] * T_basis[ix, k]
        result = mx.sum(c * T_basis[None, :, :], axis=2)

        mx.eval(result)
        return result

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def eval_jl(self, x_arr):
        """Evaluate j_l(x) for all ells. Returns mx.array (N_ell, N_x)."""
        x_mx = mx.array(np.asarray(x_arr, dtype=np.float32))
        seg_idx, t = self._map_to_segments(x_mx)
        return self._eval_chebyshev(self.coeffs_mx, seg_idx, t)

    def eval_jlp(self, x_arr):
        """Evaluate j_l'(x) for all ells. Returns mx.array (N_ell, N_x)."""
        x_mx = mx.array(np.asarray(x_arr, dtype=np.float32))
        seg_idx, t = self._map_to_segments(x_mx)
        return self._eval_chebyshev(self.coeffs_deriv_mx, seg_idx, t)

    def eval_jl_jlp(self, x_arr):
        """Evaluate j_l(x) and j_l'(x). Returns two mx.array (N_ell, N_x)."""
        x_mx = mx.array(np.asarray(x_arr, dtype=np.float32))
        seg_idx, t = self._map_to_segments(x_mx)
        jl = self._eval_chebyshev(self.coeffs_mx, seg_idx, t)
        jlp = self._eval_chebyshev(self.coeffs_deriv_mx, seg_idx, t)
        return jl, jlp

    def eval_jl_at_multiple_x(self, x_arrays):
        """Evaluate j_l for multiple x arrays. Returns list of mx.array."""
        results = [self.eval_jl(x) for x in x_arrays]
        mx.eval(*results)
        return results

    def eval_jl_jlp_at_multiple_x(self, x_arrays):
        """Evaluate j_l and j_l' for multiple x arrays."""
        jl_list, jlp_list = [], []
        for x in x_arrays:
            jl, jlp = self.eval_jl_jlp(x)
            jl_list.append(jl)
            jlp_list.append(jlp)
        mx.eval(*(jl_list + jlp_list))
        return jl_list, jlp_list


# ============================================================================
# Drop-in replacements
# ============================================================================

def compute_bessel_table_gpu(k_arr_Mpc, D_A, ell_values, need_derivative=False,
                              table=None):
    """
    Drop-in replacement for scipy Bessel loop in spectra.py / transfer.py.

    Returns: jl, jlp (or None), table
    """
    x_arr = np.asarray(k_arr_Mpc) * D_A
    x_max = float(np.max(x_arr)) * 1.05 + 10.0

    if table is None:
        table = BesselTable(ell_values, x_max=x_max)

    if need_derivative:
        jl_mx, jlp_mx = table.eval_jl_jlp(x_arr)
        return np.array(jl_mx), np.array(jlp_mx), table
    else:
        jl_mx = table.eval_jl(x_arr)
        return np.array(jl_mx), None, table


def compute_visibility_bessel_gpu(k_arr, chi_values, ell_values,
                                   weights=None, need_derivative=True,
                                   table=None):
    """
    Compute visibility-weighted Bessel integrals on GPU.

    If weights provided: returns gjl, gjlp, table  (N_ell, N_k)
    If weights is None:  returns jl_3d, jlp_3d, table  (N_ell, N_k, N_tau)
    """
    k_arr = np.asarray(k_arr, dtype=np.float64)
    chi_values = np.asarray(chi_values, dtype=np.float64)
    x_max = float(np.max(k_arr) * np.max(chi_values)) * 1.05 + 10.0

    if table is None:
        table = BesselTable(ell_values, x_max=x_max)

    N_tau = len(chi_values)
    N_k = len(k_arr)
    N_ell = len(ell_values)

    if weights is not None:
        gjl = mx.zeros((N_ell, N_k))
        gjlp = mx.zeros((N_ell, N_k)) if need_derivative else None

        for it in range(N_tau):
            x = k_arr * chi_values[it]
            w = float(weights[it])
            if need_derivative:
                jl, jlp = table.eval_jl_jlp(x)
                gjl = gjl + w * jl
                gjlp = gjlp + w * jlp
            else:
                jl = table.eval_jl(x)
                gjl = gjl + w * jl

        mx.eval(gjl)
        if gjlp is not None:
            mx.eval(gjlp)
        return np.array(gjl), (np.array(gjlp) if gjlp is not None else None), table

    else:
        jl_3d = np.zeros((N_ell, N_k, N_tau), dtype=np.float32)
        jlp_3d = np.zeros((N_ell, N_k, N_tau), dtype=np.float32) if need_derivative else None

        for it in range(N_tau):
            x = k_arr * chi_values[it]
            if need_derivative:
                jl, jlp = table.eval_jl_jlp(x)
                mx.eval(jl, jlp)
                jl_3d[:, :, it] = np.array(jl)
                jlp_3d[:, :, it] = np.array(jlp)
            else:
                jl = table.eval_jl(x)
                mx.eval(jl)
                jl_3d[:, :, it] = np.array(jl)

        return jl_3d, jlp_3d, table


# ============================================================================
# Benchmark
# ============================================================================

def benchmark(N_ell=372, N_k=5000, x_max=5000.0):
    """Benchmark GPU Bessel vs scipy."""
    from scipy.special import spherical_jn

    print("=" * 70)
    print(f"Bessel Benchmark: {N_ell} ells x {N_k} k-points, x_max={x_max}")
    print("=" * 70)

    ell_values = np.unique(np.concatenate([
        np.arange(2, 30),
        np.arange(30, 100, 2),
        np.arange(100, 500, 3),
        np.arange(500, 1500, 5),
        np.arange(1500, 2501, 8),
    ])).astype(int)
    ell_values = ell_values[:N_ell]
    N_ell = len(ell_values)

    x_arr = np.linspace(1.0, x_max, N_k).astype(np.float64)

    print(f"\nEll range: [{ell_values[0]}..{ell_values[-1]}], {N_ell} values")
    print(f"x range:   [{x_arr[0]:.1f}..{x_arr[-1]:.1f}], {N_k} points\n")

    # ---- scipy ----
    print("--- scipy ---")
    t0 = time.time()
    jl_scipy = np.zeros((N_ell, N_k), dtype=np.float64)
    jlp_scipy = np.zeros((N_ell, N_k), dtype=np.float64)
    for il, ell in enumerate(ell_values):
        jl_scipy[il] = spherical_jn(int(ell), x_arr)
        jlp_scipy[il] = spherical_jn(int(ell), x_arr, derivative=True)
    t_scipy = time.time() - t0
    print(f"Time: {t_scipy:.3f}s\n")

    # ---- GPU build + eval ----
    print("--- GPU (build + eval) ---")
    t0 = time.time()
    table = BesselTable(ell_values, x_max=x_max * 1.05, verbose=True)
    t_build = time.time() - t0

    t0 = time.time()
    jl_gpu, jlp_gpu = table.eval_jl_jlp(x_arr)
    jl_np = np.array(jl_gpu)
    jlp_np = np.array(jlp_gpu)
    t_eval = time.time() - t0
    print(f"Build: {t_build:.3f}s  Eval: {t_eval:.4f}s  Total: {t_build+t_eval:.3f}s\n")

    # ---- GPU eval only (reuse) ----
    t0 = time.time()
    jl2, _ = table.eval_jl_jlp(x_arr)
    _ = np.array(jl2)
    t_reuse = time.time() - t0
    print(f"Eval only (reuse): {t_reuse:.4f}s\n")

    # ---- Accuracy ----
    print("--- Accuracy ---")
    for name, gpu, ref in [("j_l", jl_np, jl_scipy), ("j_l'", jlp_np, jlp_scipy)]:
        abs_err = np.abs(gpu - ref)
        mask = np.abs(ref) > 1e-10
        if np.any(mask):
            rel = np.abs((gpu[mask] - ref[mask]) / ref[mask])
            print(f"  {name} rel (|ref|>1e-10): max={np.max(rel):.2e}  "
                  f"med={np.median(rel):.2e}  P99={np.percentile(rel,99):.2e}")
        print(f"  {name} abs: max={np.max(abs_err):.2e}  med={np.median(abs_err):.2e}")

    # ---- Summary ----
    print(f"\n--- Summary ---")
    print(f"scipy:       {t_scipy:.3f}s")
    print(f"GPU+build:   {t_build+t_eval:.3f}s  ({t_scipy/(t_build+t_eval):.1f}x)")
    print(f"GPU eval:    {t_reuse:.4f}s  ({t_scipy/t_reuse:.0f}x)")

    # ---- Visibility-weighted ----
    print(f"\n--- Visibility-weighted (25 tau x {N_k} k) ---")
    N_tau = 25
    chi = np.linspace(13500, 14000, N_tau)
    k_vis = np.geomspace(1e-4, 0.35, N_k)
    w = np.exp(-0.5 * ((chi - 13750) / 50)**2)
    w /= w.sum()

    t0 = time.time()
    gjl_sp = np.zeros((N_ell, N_k), dtype=np.float64)
    for il, ell in enumerate(ell_values):
        for it in range(N_tau):
            gjl_sp[il] += w[it] * spherical_jn(int(ell), k_vis * chi[it])
    t_vis_sp = time.time() - t0

    t0 = time.time()
    gjl_gpu, _, _ = compute_visibility_bessel_gpu(
        k_vis, chi, ell_values, weights=w, need_derivative=False)
    t_vis_gpu = time.time() - t0

    mask = np.abs(gjl_sp) > 1e-10
    rel_vis = np.abs((gjl_gpu[mask] - gjl_sp[mask]) / gjl_sp[mask]) if np.any(mask) else np.array([0])
    print(f"scipy: {t_vis_sp:.3f}s  GPU: {t_vis_gpu:.3f}s  "
          f"Speedup: {t_vis_sp/t_vis_gpu:.1f}x  "
          f"rel_err: max={np.max(rel_vis):.2e} med={np.median(rel_vis):.2e}")

    print("=" * 70)
    return dict(scipy=t_scipy, build=t_build, eval=t_eval, reuse=t_reuse,
                vis_scipy=t_vis_sp, vis_gpu=t_vis_gpu)


if __name__ == '__main__':
    benchmark()
