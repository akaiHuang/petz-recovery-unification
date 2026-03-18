"""
metrics.py -- Image quality metrics for diffusion evaluation.

All functions accept numpy arrays or PyTorch tensors.  Lazy imports
are used so that optional dependencies (skimage, lpips) do not cause
import-time failures when they are not installed.
"""

from __future__ import annotations

import math
from typing import Union

import numpy as np

# ===================================================================
# Internal helpers
# ===================================================================


def _to_numpy_hwc(img: Union[np.ndarray, "torch.Tensor"]) -> np.ndarray:
    """Normalise *img* to a float64 numpy array with shape (H, W, C).

    Accepted input layouts:
    * (H, W)         -> (H, W, 1)
    * (H, W, C)      -> as-is
    * (C, H, W)      -> transposed to (H, W, C) when C in {1, 3, 4}
    * (B, C, H, W)   -> squeeze batch dim, then as above

    Values are expected in [0, 1] or [0, 255]; the function normalises
    to [0, 1] if max > 1.
    """
    if not isinstance(img, np.ndarray):
        # Assume torch tensor
        img = img.detach().float().cpu().numpy()

    img = img.astype(np.float64)

    # Squeeze leading batch dim
    if img.ndim == 4:
        img = img[0]

    if img.ndim == 2:
        img = img[:, :, np.newaxis]
    elif img.ndim == 3 and img.shape[0] in (1, 3, 4) and img.shape[2] not in (1, 3, 4):
        # Likely (C, H, W) -> (H, W, C)
        img = np.transpose(img, (1, 2, 0))

    if img.max() > 1.0 + 1e-3:
        img = img / 255.0

    return np.clip(img, 0.0, 1.0)


# ===================================================================
# PSNR
# ===================================================================


def compute_psnr(
    img1: Union[np.ndarray, "torch.Tensor"],
    img2: Union[np.ndarray, "torch.Tensor"],
    max_val: float = 1.0,
) -> float:
    """Peak Signal-to-Noise Ratio (PSNR) in dB.

    Parameters
    ----------
    img1, img2 : array-like
        Images with the same spatial dimensions.
    max_val : float
        Dynamic range of the images (default 1.0 for [0, 1] images).

    Returns
    -------
    float
        PSNR in dB.  Returns ``float('inf')`` if the images are identical.
    """
    a = _to_numpy_hwc(img1)
    b = _to_numpy_hwc(img2)

    mse = np.mean((a - b) ** 2)
    if mse < 1e-15:
        return float("inf")
    return float(10.0 * math.log10(max_val**2 / mse))


# ===================================================================
# SSIM (structural similarity)
# ===================================================================


def compute_ssim(
    img1: Union[np.ndarray, "torch.Tensor"],
    img2: Union[np.ndarray, "torch.Tensor"],
    win_size: int = 7,
) -> float:
    """Structural Similarity Index (SSIM).

    Tries to use ``skimage.metrics.structural_similarity`` when
    available; otherwise falls back to a lightweight manual
    implementation (uniform window, per-channel average).

    Parameters
    ----------
    img1, img2 : array-like
    win_size : int
        Side length of the sliding window (must be odd).

    Returns
    -------
    float
        Mean SSIM in [0, 1].
    """
    a = _to_numpy_hwc(img1)
    b = _to_numpy_hwc(img2)

    try:
        from skimage.metrics import structural_similarity

        # skimage expects channel_axis for colour images
        if a.shape[-1] > 1:
            return float(
                structural_similarity(a, b, win_size=win_size, channel_axis=-1, data_range=1.0)
            )
        return float(
            structural_similarity(
                a.squeeze(-1), b.squeeze(-1), win_size=win_size, data_range=1.0
            )
        )
    except ImportError:
        pass

    # --- Manual fallback (uniform-window SSIM) ----------------------------
    return float(_ssim_manual(a, b, win_size))


def _ssim_manual(a: np.ndarray, b: np.ndarray, win_size: int) -> float:
    """Simple uniform-window SSIM (no Gaussian weighting)."""
    C1 = (0.01) ** 2
    C2 = (0.03) ** 2

    # Uniform kernel
    kernel = np.ones((win_size, win_size), dtype=np.float64) / (win_size * win_size)

    ssim_channels = []
    for c in range(a.shape[-1]):
        ac = a[:, :, c]
        bc = b[:, :, c]

        mu_a = _conv2d(ac, kernel)
        mu_b = _conv2d(bc, kernel)
        mu_a_sq = mu_a * mu_a
        mu_b_sq = mu_b * mu_b
        mu_ab = mu_a * mu_b

        sigma_a_sq = _conv2d(ac * ac, kernel) - mu_a_sq
        sigma_b_sq = _conv2d(bc * bc, kernel) - mu_b_sq
        sigma_ab = _conv2d(ac * bc, kernel) - mu_ab

        # Clamp to avoid numerical negatives
        sigma_a_sq = np.maximum(sigma_a_sq, 0.0)
        sigma_b_sq = np.maximum(sigma_b_sq, 0.0)

        num = (2.0 * mu_ab + C1) * (2.0 * sigma_ab + C2)
        den = (mu_a_sq + mu_b_sq + C1) * (sigma_a_sq + sigma_b_sq + C2)

        ssim_map = num / den
        ssim_channels.append(ssim_map.mean())

    return float(np.mean(ssim_channels))


def _conv2d(img: np.ndarray, kernel: np.ndarray) -> np.ndarray:
    """2-D convolution via numpy (valid padding, no scipy dependency)."""
    from numpy.lib.stride_tricks import sliding_window_view

    kh, kw = kernel.shape
    windows = sliding_window_view(img, (kh, kw))  # (H-kh+1, W-kw+1, kh, kw)
    return np.einsum("ijkl,kl->ij", windows, kernel)


# ===================================================================
# LPIPS (Learned Perceptual Image Patch Similarity)
# ===================================================================


def compute_lpips(
    img1: Union[np.ndarray, "torch.Tensor"],
    img2: Union[np.ndarray, "torch.Tensor"],
    net: str = "alex",
) -> float:
    """LPIPS perceptual distance (lower = more similar).

    Wraps the ``lpips`` library with lazy import.  Requires both
    ``torch`` and ``lpips`` to be installed.

    Parameters
    ----------
    img1, img2 : array-like
        Images in any layout accepted by ``_to_numpy_hwc``.
    net : str
        Backbone network for LPIPS (``"alex"``, ``"vgg"``, ``"squeeze"``).

    Returns
    -------
    float
        LPIPS distance.

    Raises
    ------
    ImportError
        If ``lpips`` or ``torch`` is not installed.
    """
    import torch

    try:
        import lpips as _lpips_lib
    except ImportError:
        raise ImportError(
            "lpips is required for compute_lpips(). "
            "Install it with:  pip install lpips"
        )

    a = _to_numpy_hwc(img1)
    b = _to_numpy_hwc(img2)

    # LPIPS expects (B, C, H, W) tensors in [-1, 1]
    def _prep(x: np.ndarray) -> torch.Tensor:
        # (H, W, C) -> (1, C, H, W), [0,1] -> [-1,1]
        t = torch.from_numpy(x).permute(2, 0, 1).unsqueeze(0).float()
        t = t * 2.0 - 1.0
        # Ensure 3 channels (LPIPS expects RGB)
        if t.shape[1] == 1:
            t = t.expand(-1, 3, -1, -1)
        return t

    ta = _prep(a)
    tb = _prep(b)

    # Lazy singleton so we don't re-init the network on every call
    if not hasattr(compute_lpips, "_model_cache"):
        compute_lpips._model_cache = {}  # type: ignore[attr-defined]

    if net not in compute_lpips._model_cache:  # type: ignore[attr-defined]
        model = _lpips_lib.LPIPS(net=net, verbose=False)
        model.eval()
        compute_lpips._model_cache[net] = model  # type: ignore[attr-defined]

    model = compute_lpips._model_cache[net]  # type: ignore[attr-defined]

    with torch.no_grad():
        dist = model(ta, tb)

    return float(dist.item())
