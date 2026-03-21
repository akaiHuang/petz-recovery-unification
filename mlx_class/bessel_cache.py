"""
bessel_cache.py -- Disk-cached Bessel table for sub-second repeat runs.

First run:  builds Chebyshev Bessel table (~2-4s), saves to ~/.mlx_class/
Cached run: loads from disk (~0.05s), skips scipy entirely.

Cache key: hash of (ell_values, x_max, N_cheb, seg_width) ensures
correctness when parameters change.

Author: Sheng-Kai Huang, 2026
"""
import os
import hashlib
import time
import json
import numpy as np
import mlx.core as mx

from mlx_class.bessel_gpu import BesselTable

# Default cache directory
CACHE_DIR = os.path.expanduser("~/.mlx_class")


def _cache_key(ell_values, x_max, N_cheb, seg_width):
    """
    Compute a deterministic hash from Bessel table parameters.

    Returns (hex_str, params_dict) where hex_str is a short hash
    and params_dict stores the full parameters for verification.
    """
    ell_arr = np.asarray(ell_values, dtype=np.int32)
    # Build a canonical byte representation
    h = hashlib.sha256()
    h.update(ell_arr.tobytes())
    h.update(np.float64(x_max).tobytes())
    h.update(np.int32(N_cheb).tobytes())
    h.update(np.float64(seg_width).tobytes())
    hex_str = h.hexdigest()[:16]

    params = {
        "ell_min": int(ell_arr[0]),
        "ell_max": int(ell_arr[-1]),
        "n_ell": len(ell_arr),
        "x_max": float(x_max),
        "N_cheb": int(N_cheb),
        "seg_width": float(seg_width),
    }
    return hex_str, params


def _cache_path(hex_key, cache_dir=CACHE_DIR):
    """Return path for cached npz file."""
    return os.path.join(cache_dir, f"bessel_{hex_key}.npz")


def _meta_path(hex_key, cache_dir=CACHE_DIR):
    """Return path for cache metadata JSON."""
    return os.path.join(cache_dir, f"bessel_{hex_key}.json")


class CachedBesselTable(BesselTable):
    """
    BesselTable with transparent disk caching.

    On first call: builds table via scipy + forward recurrence, saves
    post-DCT coefficients (float32) to ~/.mlx_class/.
    On subsequent calls with same parameters: loads from disk in ~0.05s.

    Usage (drop-in replacement for BesselTable):
        table = CachedBesselTable(ell_values, x_max=5500.0)
        jl = table.eval_jl(x_arr)
    """

    def __init__(self, ell_values, x_max=5500.0, N_cheb=64,
                 seg_width=80.0, verbose=True, cache_dir=CACHE_DIR,
                 force_rebuild=False):
        """
        Build or load cached Chebyshev coefficient table.

        Parameters
        ----------
        ell_values : array-like of int
        x_max : float
        N_cheb : int
        seg_width : float
        verbose : bool
        cache_dir : str
            Directory for cache files. Created if needed.
        force_rebuild : bool
            If True, ignore existing cache and rebuild.
        """
        self.ell_values = np.asarray(ell_values, dtype=np.int32)
        self.N_ell = len(self.ell_values)
        self.x_max = float(x_max)
        self.N_cheb = int(N_cheb)

        hex_key, params = _cache_key(ell_values, x_max, N_cheb, seg_width)
        npz_path = _cache_path(hex_key, cache_dir)
        meta_path = _meta_path(hex_key, cache_dir)

        t0 = time.time()

        if not force_rebuild and os.path.exists(npz_path) and os.path.exists(meta_path):
            # ---- CACHE HIT ----
            try:
                data = np.load(npz_path)
                coeffs = data["coeffs"]        # (N_ell, n_seg, N_cheb) float32
                coeffs_d = data["coeffs_d"]    # (N_ell, n_seg, N_cheb) float32
                n_seg = int(data["n_seg"])
                seg_width_actual = float(data["seg_width_actual"])

                # Sanity check dimensions
                assert coeffs.shape[0] == self.N_ell, \
                    f"Cache N_ell mismatch: {coeffs.shape[0]} vs {self.N_ell}"
                assert coeffs.shape[2] == self.N_cheb, \
                    f"Cache N_cheb mismatch: {coeffs.shape[2]} vs {self.N_cheb}"

                self.n_seg = n_seg
                self.seg_width_actual = seg_width_actual
                self.coeffs_mx = mx.array(coeffs)
                self.coeffs_deriv_mx = mx.array(coeffs_d)

                dt = time.time() - t0
                if verbose:
                    print(f"[BesselCache] Loaded: {self.N_ell} ells, "
                          f"{n_seg} segs x {self.N_cheb} nodes, "
                          f"x_max={self.x_max:.0f}, "
                          f"time={dt:.3f}s (cached)")
                return

            except Exception as e:
                if verbose:
                    print(f"[BesselCache] Cache load failed ({e}), rebuilding...")

        # ---- CACHE MISS: build from scratch ----
        self._build_table(seg_width)
        dt_build = time.time() - t0

        # Save to cache
        t_save = time.time()
        try:
            os.makedirs(cache_dir, exist_ok=True)
            np.savez_compressed(
                npz_path,
                coeffs=np.array(self.coeffs_mx),
                coeffs_d=np.array(self.coeffs_deriv_mx),
                n_seg=np.array(self.n_seg),
                seg_width_actual=np.array(self.seg_width_actual),
            )
            with open(meta_path, "w") as f:
                json.dump(params, f, indent=2)
            dt_save = time.time() - t_save
            if verbose:
                cache_size_mb = os.path.getsize(npz_path) / 1024 / 1024
                print(f"[BesselCache] Built: {self.N_ell} ells, "
                      f"{self.n_seg} segs x {self.N_cheb} nodes, "
                      f"x_max={self.x_max:.0f}, "
                      f"build={dt_build:.2f}s, save={dt_save:.2f}s "
                      f"({cache_size_mb:.1f} MB)")
        except Exception as e:
            if verbose:
                print(f"[BesselCache] Built: {dt_build:.2f}s (save failed: {e})")


def clear_cache(cache_dir=CACHE_DIR):
    """Remove all cached Bessel tables."""
    if not os.path.exists(cache_dir):
        print("[BesselCache] No cache directory found.")
        return
    removed = 0
    for f in os.listdir(cache_dir):
        if f.startswith("bessel_") and (f.endswith(".npz") or f.endswith(".json")):
            os.remove(os.path.join(cache_dir, f))
            removed += 1
    print(f"[BesselCache] Removed {removed} cache files from {cache_dir}")


def list_cache(cache_dir=CACHE_DIR):
    """List all cached Bessel tables with their parameters."""
    if not os.path.exists(cache_dir):
        print("[BesselCache] No cache directory found.")
        return []
    entries = []
    for f in os.listdir(cache_dir):
        if f.startswith("bessel_") and f.endswith(".json"):
            meta_path = os.path.join(cache_dir, f)
            hex_key = f.replace("bessel_", "").replace(".json", "")
            npz_path = _cache_path(hex_key, cache_dir)
            try:
                with open(meta_path) as fp:
                    params = json.load(fp)
                size_mb = os.path.getsize(npz_path) / 1024 / 1024 if os.path.exists(npz_path) else 0
                entries.append({"key": hex_key, "params": params, "size_mb": size_mb})
                print(f"  {hex_key}: {params['n_ell']} ells, "
                      f"x_max={params['x_max']:.0f}, "
                      f"N_cheb={params['N_cheb']}, "
                      f"{size_mb:.1f} MB")
            except Exception:
                pass
    if not entries:
        print("[BesselCache] No cached tables found.")
    return entries


# ============================================================================
# Benchmark: cached vs uncached
# ============================================================================

if __name__ == "__main__":
    ell_values = np.unique(np.concatenate([
        np.arange(2, 30, 1),
        np.arange(30, 100, 2),
        np.arange(100, 500, 4),
        np.arange(500, 1500, 8),
        np.arange(1500, 2501, 12),
    ])).astype(int)

    x_max = 5500.0
    print("=" * 60)
    print("Bessel Cache Benchmark")
    print("=" * 60)

    # Clear cache for clean test
    clear_cache()

    # First run (cache miss)
    print("\n--- First run (build + save) ---")
    t0 = time.time()
    table1 = CachedBesselTable(ell_values, x_max=x_max)
    t1 = time.time() - t0
    print(f"Total: {t1:.3f}s")

    # Second run (cache hit)
    print("\n--- Second run (cache load) ---")
    t0 = time.time()
    table2 = CachedBesselTable(ell_values, x_max=x_max)
    t2 = time.time() - t0
    print(f"Total: {t2:.3f}s")

    # Verify equivalence
    x_test = np.linspace(1.0, 5000.0, 1000)
    jl1 = np.array(table1.eval_jl(x_test))
    jl2 = np.array(table2.eval_jl(x_test))
    max_diff = np.max(np.abs(jl1 - jl2))
    print(f"\nMax difference between cached and fresh: {max_diff:.2e}")
    print(f"Speedup: {t1/t2:.1f}x")

    # List cache
    print("\n--- Cache contents ---")
    list_cache()
    print("=" * 60)
