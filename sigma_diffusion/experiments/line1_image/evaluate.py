#!/usr/bin/env python3
"""
Line 1 Evaluation: Compare baseline vs Sigma-guided images.

Computes:
    - PSNR (Peak Signal-to-Noise Ratio)
    - SSIM (Structural Similarity Index)
    - LPIPS (Learned Perceptual Image Patch Similarity) -- if available
    - Side-by-side comparison grid

Usage:
    python evaluate.py --baseline-dir ../../results/line1/baseline \
                       --guided-dir ../../results/line1/guided \
                       --output-dir ../../results/line1/evaluation
"""

from __future__ import annotations

import argparse
import json
import logging
import os
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------


def compute_psnr(img1: np.ndarray, img2: np.ndarray, max_val: float = 255.0) -> float:
    """
    Peak Signal-to-Noise Ratio between two images.
    Higher is better (identical images -> inf).
    """
    mse = np.mean((img1.astype(np.float64) - img2.astype(np.float64)) ** 2)
    if mse == 0:
        return float("inf")
    return float(20.0 * np.log10(max_val / np.sqrt(mse)))


def compute_ssim(
    img1: np.ndarray,
    img2: np.ndarray,
    window_size: int = 7,
    max_val: float = 255.0,
) -> float:
    """
    Structural Similarity Index (simplified, per-channel average).
    Range: [-1, 1], where 1 = identical.

    Uses a sliding window approach with Gaussian-like weighting.
    """
    C1 = (0.01 * max_val) ** 2
    C2 = (0.03 * max_val) ** 2

    img1 = img1.astype(np.float64)
    img2 = img2.astype(np.float64)

    # If images are RGB, compute SSIM per channel and average
    if img1.ndim == 3:
        ssim_channels = []
        for c in range(img1.shape[2]):
            ssim_channels.append(
                _ssim_single_channel(img1[:, :, c], img2[:, :, c], C1, C2, window_size)
            )
        return float(np.mean(ssim_channels))
    else:
        return _ssim_single_channel(img1, img2, C1, C2, window_size)


def _ssim_single_channel(
    img1: np.ndarray,
    img2: np.ndarray,
    C1: float,
    C2: float,
    window_size: int,
) -> float:
    """SSIM for a single channel using uniform window (fast approximation)."""
    from scipy.ndimage import uniform_filter

    mu1 = uniform_filter(img1, size=window_size)
    mu2 = uniform_filter(img2, size=window_size)

    mu1_sq = mu1 ** 2
    mu2_sq = mu2 ** 2
    mu1_mu2 = mu1 * mu2

    sigma1_sq = uniform_filter(img1 ** 2, size=window_size) - mu1_sq
    sigma2_sq = uniform_filter(img2 ** 2, size=window_size) - mu2_sq
    sigma12 = uniform_filter(img1 * img2, size=window_size) - mu1_mu2

    # Clamp to avoid negative variances from numerical issues
    sigma1_sq = np.maximum(sigma1_sq, 0)
    sigma2_sq = np.maximum(sigma2_sq, 0)

    numerator = (2 * mu1_mu2 + C1) * (2 * sigma12 + C2)
    denominator = (mu1_sq + mu2_sq + C1) * (sigma1_sq + sigma2_sq + C2)

    ssim_map = numerator / denominator
    return float(np.mean(ssim_map))


def _ssim_fallback(
    img1: np.ndarray,
    img2: np.ndarray,
    max_val: float = 255.0,
) -> float:
    """
    Ultra-simple SSIM fallback when scipy is not available.
    Uses global statistics instead of sliding window.
    """
    C1 = (0.01 * max_val) ** 2
    C2 = (0.03 * max_val) ** 2

    img1 = img1.astype(np.float64)
    img2 = img2.astype(np.float64)

    mu1 = np.mean(img1)
    mu2 = np.mean(img2)
    sigma1_sq = np.var(img1)
    sigma2_sq = np.var(img2)
    sigma12 = np.mean((img1 - mu1) * (img2 - mu2))

    numerator = (2 * mu1 * mu2 + C1) * (2 * sigma12 + C2)
    denominator = (mu1 ** 2 + mu2 ** 2 + C1) * (sigma1_sq + sigma2_sq + C2)

    return float(numerator / denominator)


def compute_ssim_safe(
    img1: np.ndarray,
    img2: np.ndarray,
    max_val: float = 255.0,
) -> float:
    """Compute SSIM with fallback if scipy is not available."""
    try:
        return compute_ssim(img1, img2, max_val=max_val)
    except ImportError:
        logger.warning("scipy not available; using global SSIM approximation.")
        return _ssim_fallback(img1, img2, max_val=max_val)


def compute_lpips(
    img1: np.ndarray,
    img2: np.ndarray,
) -> Optional[float]:
    """
    Compute LPIPS (Learned Perceptual Image Patch Similarity).
    Lower is better (0 = identical).

    Requires the `lpips` package. Returns None if not available.
    """
    try:
        import torch
        import lpips as lpips_module
    except ImportError:
        logger.warning(
            "lpips not available. Install with: pip install lpips"
        )
        return None

    try:
        loss_fn = lpips_module.LPIPS(net="alex", verbose=False)

        # Convert to torch tensors: [B, C, H, W] in [-1, 1]
        def to_tensor(img: np.ndarray) -> torch.Tensor:
            if img.ndim == 2:
                img = np.stack([img, img, img], axis=-1)
            t = torch.from_numpy(img.astype(np.float32)).permute(2, 0, 1)
            t = t.unsqueeze(0) / 127.5 - 1.0  # normalize to [-1, 1]
            return t

        t1 = to_tensor(img1)
        t2 = to_tensor(img2)

        # Resize to match if needed
        if t1.shape != t2.shape:
            import torch.nn.functional as F
            target_h = min(t1.shape[2], t2.shape[2])
            target_w = min(t1.shape[3], t2.shape[3])
            t1 = F.interpolate(t1, size=(target_h, target_w), mode="bilinear")
            t2 = F.interpolate(t2, size=(target_h, target_w), mode="bilinear")

        with torch.no_grad():
            d = loss_fn(t1, t2)
        return float(d.item())

    except Exception as e:
        logger.warning(f"LPIPS computation failed: {e}")
        return None


# ---------------------------------------------------------------------------
# Image loading
# ---------------------------------------------------------------------------


def load_image(path: str) -> Optional[np.ndarray]:
    """Load an image as a numpy array (H, W, C) in uint8."""
    try:
        from PIL import Image
        img = Image.open(path).convert("RGB")
        return np.array(img)
    except ImportError:
        logger.error("Pillow not available. Install with: pip install Pillow")
        return None
    except Exception as e:
        logger.error(f"Failed to load image {path}: {e}")
        return None


# ---------------------------------------------------------------------------
# Comparison grid
# ---------------------------------------------------------------------------


def make_comparison_grid(
    baseline_images: List[np.ndarray],
    guided_images: List[np.ndarray],
    prompts: List[str],
    metrics: List[Dict[str, Any]],
    output_path: str,
    max_display: int = 10,
) -> None:
    """
    Generate a side-by-side comparison grid:
        Left column: baseline images
        Right column: guided images
        Labels: prompt (truncated) + PSNR/SSIM
    """
    try:
        import matplotlib
        matplotlib.use("Agg")  # Non-interactive backend
        import matplotlib.pyplot as plt
        from matplotlib.gridspec import GridSpec
    except ImportError:
        logger.error(
            "matplotlib not available. Install with: pip install matplotlib"
        )
        return

    n = min(len(baseline_images), len(guided_images), max_display)
    if n == 0:
        logger.warning("No images to display.")
        return

    fig_height = max(4 * n, 8)
    fig, axes = plt.subplots(n, 2, figsize=(14, fig_height))

    if n == 1:
        axes = axes.reshape(1, 2)

    fig.suptitle(
        "Line 1: Baseline vs Sigma-Guided Diffusion",
        fontsize=16,
        fontweight="bold",
        y=0.98,
    )

    for i in range(n):
        # Baseline
        axes[i, 0].imshow(baseline_images[i])
        axes[i, 0].set_title("Baseline", fontsize=10)
        axes[i, 0].axis("off")

        # Guided
        axes[i, 1].imshow(guided_images[i])

        # Metrics label
        m = metrics[i] if i < len(metrics) else {}
        psnr = m.get("psnr", "N/A")
        ssim = m.get("ssim", "N/A")
        lpips_val = m.get("lpips", None)

        psnr_str = f"{psnr:.2f}" if isinstance(psnr, float) else str(psnr)
        ssim_str = f"{ssim:.4f}" if isinstance(ssim, float) else str(ssim)
        label = f"Guided | PSNR={psnr_str} SSIM={ssim_str}"
        if lpips_val is not None:
            label += f" LPIPS={lpips_val:.4f}"

        axes[i, 1].set_title(label, fontsize=10)
        axes[i, 1].axis("off")

        # Prompt label on the left
        prompt_short = prompts[i][:50] + ("..." if len(prompts[i]) > 50 else "")
        axes[i, 0].set_ylabel(
            prompt_short,
            fontsize=8,
            rotation=0,
            labelpad=120,
            va="center",
        )

    plt.tight_layout(rect=[0.08, 0.02, 1.0, 0.96])
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    logger.info(f"Saved comparison grid to {output_path}")


# ---------------------------------------------------------------------------
# Main evaluation
# ---------------------------------------------------------------------------


@dataclass
class ImageMetrics:
    """Metrics for a single baseline-guided image pair."""
    prompt: str
    prompt_index: int
    psnr: float
    ssim: float
    lpips: Optional[float]
    baseline_path: str
    guided_path: str


def evaluate(
    baseline_dir: Optional[str] = None,
    guided_dir: Optional[str] = None,
    output_dir: Optional[str] = None,
    compute_lpips_flag: bool = True,
    max_grid_images: int = 10,
) -> Dict[str, Any]:
    """
    Evaluate baseline vs Sigma-guided images.

    Returns a summary dict with aggregate metrics.
    """
    base_results_dir = Path(__file__).parent.parent.parent / "results" / "line1"

    if baseline_dir is None:
        baseline_dir = str(base_results_dir / "baseline")
    if guided_dir is None:
        guided_dir = str(base_results_dir / "guided")
    if output_dir is None:
        output_dir = str(base_results_dir / "evaluation")

    os.makedirs(output_dir, exist_ok=True)

    # -----------------------------------------------------------------------
    # Load metadata to find prompt count and image paths
    # -----------------------------------------------------------------------
    baseline_sigma_path = os.path.join(baseline_dir, "sigma_profile.json")
    guided_sigma_path = os.path.join(guided_dir, "sigma_profile.json")

    if os.path.exists(baseline_sigma_path):
        with open(baseline_sigma_path, "r") as f:
            baseline_profiles = json.load(f)
    else:
        logger.warning("No baseline sigma profile found; scanning images directly.")
        baseline_profiles = []

    if os.path.exists(guided_sigma_path):
        with open(guided_sigma_path, "r") as f:
            guided_profiles = json.load(f)
    else:
        logger.warning("No guided sigma profile found; scanning images directly.")
        guided_profiles = []

    # Find matching image pairs
    baseline_img_dir = os.path.join(baseline_dir, "images")
    guided_img_dir = os.path.join(guided_dir, "images")

    if not os.path.isdir(baseline_img_dir):
        raise FileNotFoundError(f"Baseline images not found at {baseline_img_dir}")
    if not os.path.isdir(guided_img_dir):
        raise FileNotFoundError(f"Guided images not found at {guided_img_dir}")

    # Find matching pairs by index
    baseline_files = sorted([
        f for f in os.listdir(baseline_img_dir) if f.endswith(".png")
    ])
    guided_files = sorted([
        f for f in os.listdir(guided_img_dir) if f.endswith(".png")
    ])

    # Build index-based matching
    def extract_index(filename: str) -> int:
        """Extract numeric index from filename like 'baseline_003.png'."""
        name = os.path.splitext(filename)[0]
        parts = name.split("_")
        for part in reversed(parts):
            try:
                return int(part)
            except ValueError:
                continue
        return -1

    baseline_by_idx = {extract_index(f): f for f in baseline_files}
    guided_by_idx = {extract_index(f): f for f in guided_files}

    common_indices = sorted(
        set(baseline_by_idx.keys()) & set(guided_by_idx.keys())
    )

    if not common_indices:
        logger.error("No matching image pairs found.")
        return {"error": "No matching image pairs"}

    logger.info(f"Found {len(common_indices)} matching image pairs.")

    # -----------------------------------------------------------------------
    # Compute metrics for each pair
    # -----------------------------------------------------------------------
    all_metrics: List[ImageMetrics] = []
    baseline_images: List[np.ndarray] = []
    guided_images: List[np.ndarray] = []
    prompts: List[str] = []

    for idx in common_indices:
        baseline_path = os.path.join(baseline_img_dir, baseline_by_idx[idx])
        guided_path = os.path.join(guided_img_dir, guided_by_idx[idx])

        img_b = load_image(baseline_path)
        img_g = load_image(guided_path)

        if img_b is None or img_g is None:
            logger.warning(f"Skipping pair {idx}: failed to load images.")
            continue

        # Resize if dimensions differ
        if img_b.shape != img_g.shape:
            from PIL import Image
            target_size = (
                min(img_b.shape[1], img_g.shape[1]),
                min(img_b.shape[0], img_g.shape[0]),
            )
            img_b_pil = Image.fromarray(img_b).resize(target_size, Image.LANCZOS)
            img_g_pil = Image.fromarray(img_g).resize(target_size, Image.LANCZOS)
            img_b = np.array(img_b_pil)
            img_g = np.array(img_g_pil)

        # PSNR
        psnr = compute_psnr(img_b, img_g)

        # SSIM
        ssim = compute_ssim_safe(img_b, img_g)

        # LPIPS
        lpips_val = None
        if compute_lpips_flag:
            lpips_val = compute_lpips(img_b, img_g)

        # Get prompt
        prompt = f"prompt_{idx}"
        if idx < len(baseline_profiles):
            prompt = baseline_profiles[idx].get("prompt", prompt)

        all_metrics.append(ImageMetrics(
            prompt=prompt,
            prompt_index=idx,
            psnr=psnr,
            ssim=ssim,
            lpips=lpips_val,
            baseline_path=baseline_path,
            guided_path=guided_path,
        ))

        baseline_images.append(img_b)
        guided_images.append(img_g)
        prompts.append(prompt)

        logger.info(
            f"  [{idx}] PSNR={psnr:.2f} SSIM={ssim:.4f}"
            + (f" LPIPS={lpips_val:.4f}" if lpips_val is not None else "")
        )

    # -----------------------------------------------------------------------
    # Aggregate metrics
    # -----------------------------------------------------------------------
    psnr_values = [m.psnr for m in all_metrics if np.isfinite(m.psnr)]
    ssim_values = [m.ssim for m in all_metrics]
    lpips_values = [m.lpips for m in all_metrics if m.lpips is not None]

    summary = {
        "num_pairs": len(all_metrics),
        "psnr": {
            "mean": float(np.mean(psnr_values)) if psnr_values else None,
            "std": float(np.std(psnr_values)) if psnr_values else None,
            "min": float(np.min(psnr_values)) if psnr_values else None,
            "max": float(np.max(psnr_values)) if psnr_values else None,
        },
        "ssim": {
            "mean": float(np.mean(ssim_values)) if ssim_values else None,
            "std": float(np.std(ssim_values)) if ssim_values else None,
            "min": float(np.min(ssim_values)) if ssim_values else None,
            "max": float(np.max(ssim_values)) if ssim_values else None,
        },
        "lpips": {
            "mean": float(np.mean(lpips_values)) if lpips_values else None,
            "std": float(np.std(lpips_values)) if lpips_values else None,
            "min": float(np.min(lpips_values)) if lpips_values else None,
            "max": float(np.max(lpips_values)) if lpips_values else None,
        } if lpips_values else None,
        "per_image": [asdict(m) for m in all_metrics],
    }

    # Load speedup metrics if available
    speedup_path = os.path.join(guided_dir, "speedup_metrics.json")
    if os.path.exists(speedup_path):
        with open(speedup_path, "r") as f:
            speedup = json.load(f)
        summary["speedup"] = {
            "step_reduction_pct": speedup.get("step_reduction_pct"),
            "wall_clock_speedup": speedup.get("wall_clock_speedup"),
            "theoretical_speedup": speedup.get("speedup_ratio_theoretical"),
        }

    # -----------------------------------------------------------------------
    # Save results
    # -----------------------------------------------------------------------
    metrics_path = os.path.join(output_dir, "metrics.json")
    with open(metrics_path, "w") as f:
        json.dump(summary, f, indent=2)
    logger.info(f"Saved metrics to {metrics_path}")

    # Comparison grid
    grid_path = os.path.join(output_dir, "comparison_grid.png")
    make_comparison_grid(
        baseline_images=baseline_images,
        guided_images=guided_images,
        prompts=prompts,
        metrics=[asdict(m) for m in all_metrics],
        output_path=grid_path,
        max_display=max_grid_images,
    )

    # -----------------------------------------------------------------------
    # Print summary
    # -----------------------------------------------------------------------
    logger.info("=" * 60)
    logger.info("EVALUATION SUMMARY")
    logger.info("=" * 60)
    logger.info(f"  Image pairs evaluated: {len(all_metrics)}")
    if psnr_values:
        logger.info(
            f"  PSNR: {np.mean(psnr_values):.2f} +/- {np.std(psnr_values):.2f} dB"
        )
    if ssim_values:
        logger.info(
            f"  SSIM: {np.mean(ssim_values):.4f} +/- {np.std(ssim_values):.4f}"
        )
    if lpips_values:
        logger.info(
            f"  LPIPS: {np.mean(lpips_values):.4f} +/- {np.std(lpips_values):.4f}"
        )
    if "speedup" in summary:
        sp = summary["speedup"]
        logger.info(
            f"  Step reduction: {sp['step_reduction_pct']:.1f}%"
        )
        logger.info(
            f"  Wall-clock speedup: {sp['wall_clock_speedup']:.2f}x"
        )

    # Quality verdict
    if psnr_values and np.mean(psnr_values) >= 25.0:
        logger.info("  VERDICT: Quality preserved (PSNR >= 25 dB)")
    elif psnr_values and np.mean(psnr_values) >= 20.0:
        logger.info("  VERDICT: Acceptable quality (PSNR >= 20 dB)")
    elif psnr_values:
        logger.info("  VERDICT: Quality degradation detected (PSNR < 20 dB)")

    logger.info("=" * 60)

    return summary


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description="Line 1 Evaluation: Compare baseline vs Sigma-guided images."
    )
    parser.add_argument(
        "--baseline-dir",
        type=str,
        default=None,
        help="Path to baseline results directory.",
    )
    parser.add_argument(
        "--guided-dir",
        type=str,
        default=None,
        help="Path to guided results directory.",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=None,
        help="Output directory for evaluation results.",
    )
    parser.add_argument(
        "--no-lpips",
        action="store_true",
        help="Skip LPIPS computation.",
    )
    parser.add_argument(
        "--max-grid-images",
        type=int,
        default=10,
        help="Maximum images in comparison grid (default: 10).",
    )
    args = parser.parse_args()

    evaluate(
        baseline_dir=args.baseline_dir,
        guided_dir=args.guided_dir,
        output_dir=args.output_dir,
        compute_lpips_flag=not args.no_lpips,
        max_grid_images=args.max_grid_images,
    )


if __name__ == "__main__":
    main()
