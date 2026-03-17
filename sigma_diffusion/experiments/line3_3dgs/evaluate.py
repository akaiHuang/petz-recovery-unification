#!/usr/bin/env python3
"""
Line 3 Experiment: Evaluation of Baseline vs Sigma-Guided 3DGS.

Compares rendered views from baseline (30k iterations) vs sigma-guided
(early-stopped) training. Computes PSNR, SSIM, and generates side-by-side
comparison images.

Note: OpenSplat renders views during training. For post-hoc rendering,
we use the final PLY files and a simple Gaussian splatting renderer,
or compare training outputs if available.

Usage:
    python evaluate.py
    python evaluate.py --baseline-dir /path/to/baseline --guided-dir /path/to/guided
    python evaluate.py --render   # Also render new views from PLY files
"""

import argparse
import json
import math
import sys
import time
from pathlib import Path
from typing import Optional

import numpy as np

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[3]
RESULTS_DIR = PROJECT_ROOT / "sigma_diffusion" / "results" / "line3"
CHECKPOINT_DIR = RESULTS_DIR / "checkpoints"
GUIDED_DIR = RESULTS_DIR / "guided"


# ---------------------------------------------------------------------------
# Image Quality Metrics
# ---------------------------------------------------------------------------
def compute_psnr(img1: np.ndarray, img2: np.ndarray) -> float:
    """Compute Peak Signal-to-Noise Ratio between two images.

    Args:
        img1, img2: numpy arrays of shape (H, W, C), values in [0, 1] or [0, 255].

    Returns:
        PSNR in dB. Higher is better.
    """
    if img1.shape != img2.shape:
        raise ValueError(f"Image shapes don't match: {img1.shape} vs {img2.shape}")

    # Ensure float
    img1 = img1.astype(np.float64)
    img2 = img2.astype(np.float64)

    # Determine value range
    max_val = 255.0 if img1.max() > 1.0 else 1.0

    mse = np.mean((img1 - img2) ** 2)
    if mse == 0:
        return float('inf')

    return 20 * math.log10(max_val) - 10 * math.log10(mse)


def compute_ssim(
    img1: np.ndarray,
    img2: np.ndarray,
    window_size: int = 11,
    C1: float = 0.01 ** 2,
    C2: float = 0.03 ** 2,
) -> float:
    """Compute Structural Similarity Index (SSIM) between two images.

    Implements the SSIM algorithm from Wang et al. (2004).

    Args:
        img1, img2: numpy arrays of shape (H, W, C), values in [0, 1] or [0, 255].
        window_size: Size of the Gaussian window.
        C1, C2: Stability constants.

    Returns:
        SSIM value in [0, 1]. Higher is better.
    """
    if img1.shape != img2.shape:
        raise ValueError(f"Image shapes don't match: {img1.shape} vs {img2.shape}")

    img1 = img1.astype(np.float64)
    img2 = img2.astype(np.float64)

    # Normalize to [0, 1]
    if img1.max() > 1.0:
        img1 = img1 / 255.0
        img2 = img2 / 255.0

    # Handle multi-channel: compute per-channel and average
    if img1.ndim == 3:
        ssim_channels = []
        for c in range(img1.shape[2]):
            ssim_channels.append(
                _compute_ssim_single_channel(img1[:, :, c], img2[:, :, c],
                                              window_size, C1, C2)
            )
        return float(np.mean(ssim_channels))
    else:
        return _compute_ssim_single_channel(img1, img2, window_size, C1, C2)


def _compute_ssim_single_channel(
    img1: np.ndarray,
    img2: np.ndarray,
    window_size: int,
    C1: float,
    C2: float,
) -> float:
    """SSIM for a single channel using uniform window (fast approximation)."""
    from scipy.ndimage import uniform_filter

    pad = window_size // 2

    # Local statistics
    mu1 = uniform_filter(img1, size=window_size)
    mu2 = uniform_filter(img2, size=window_size)

    mu1_sq = mu1 ** 2
    mu2_sq = mu2 ** 2
    mu1_mu2 = mu1 * mu2

    sigma1_sq = uniform_filter(img1 ** 2, size=window_size) - mu1_sq
    sigma2_sq = uniform_filter(img2 ** 2, size=window_size) - mu2_sq
    sigma12 = uniform_filter(img1 * img2, size=window_size) - mu1_mu2

    # Clamp negative variances (numerical artifact)
    sigma1_sq = np.maximum(sigma1_sq, 0)
    sigma2_sq = np.maximum(sigma2_sq, 0)

    # SSIM formula
    numerator = (2 * mu1_mu2 + C1) * (2 * sigma12 + C2)
    denominator = (mu1_sq + mu2_sq + C1) * (sigma1_sq + sigma2_sq + C2)

    ssim_map = numerator / denominator

    # Crop padding region and average
    return float(ssim_map[pad:-pad, pad:-pad].mean())


def _compute_ssim_fallback(img1: np.ndarray, img2: np.ndarray) -> float:
    """Fallback SSIM without scipy (simpler block-based approach)."""
    img1 = img1.astype(np.float64)
    img2 = img2.astype(np.float64)

    if img1.max() > 1.0:
        img1 /= 255.0
        img2 /= 255.0

    C1 = 0.01 ** 2
    C2 = 0.03 ** 2

    # Use 8x8 blocks
    block_size = 8
    H, W = img1.shape[:2]

    ssim_values = []
    for i in range(0, H - block_size, block_size):
        for j in range(0, W - block_size, block_size):
            block1 = img1[i:i + block_size, j:j + block_size].flatten()
            block2 = img2[i:i + block_size, j:j + block_size].flatten()

            mu1 = block1.mean()
            mu2 = block2.mean()
            sigma1_sq = block1.var()
            sigma2_sq = block2.var()
            sigma12 = np.cov(block1, block2)[0, 1]

            num = (2 * mu1 * mu2 + C1) * (2 * sigma12 + C2)
            den = (mu1 ** 2 + mu2 ** 2 + C1) * (sigma1_sq + sigma2_sq + C2)
            ssim_values.append(num / den)

    return float(np.mean(ssim_values))


# ---------------------------------------------------------------------------
# Image Loading and Comparison
# ---------------------------------------------------------------------------
def load_image(path: Path) -> Optional[np.ndarray]:
    """Load an image file as numpy array."""
    try:
        from PIL import Image
        img = Image.open(path).convert("RGB")
        return np.array(img)
    except ImportError:
        print("[WARN] Pillow not available. Install with: pip install Pillow")
        return None
    except Exception as e:
        print(f"[WARN] Could not load image {path}: {e}")
        return None


def find_rendered_images(directory: Path) -> list[Path]:
    """Find rendered output images in a directory."""
    image_extensions = {".png", ".jpg", ".jpeg", ".bmp", ".tiff"}
    images = []
    for ext in image_extensions:
        images.extend(directory.glob(f"*{ext}"))
        images.extend(directory.glob(f"renders/*{ext}"))
        images.extend(directory.glob(f"output/*{ext}"))
    return sorted(images)


def find_matching_pairs(
    baseline_images: list[Path],
    guided_images: list[Path],
) -> list[tuple[Path, Path]]:
    """Match baseline and guided images by filename."""
    baseline_names = {p.stem: p for p in baseline_images}
    guided_names = {p.stem: p for p in guided_images}

    pairs = []
    for name in sorted(set(baseline_names) & set(guided_names)):
        pairs.append((baseline_names[name], guided_names[name]))

    return pairs


def find_ground_truth_images(data_dir: Path) -> list[Path]:
    """Find ground truth images from the original dataset."""
    gt_dirs = [
        data_dir / "images",
        data_dir / "test",
        data_dir / "images_test",
    ]
    for d in gt_dirs:
        if d.is_dir():
            images = sorted(d.glob("*.png")) + sorted(d.glob("*.jpg"))
            if images:
                return images
    return []


# ---------------------------------------------------------------------------
# PLY-based Quality Estimation
# ---------------------------------------------------------------------------
def estimate_quality_from_ply(ply_path: Path) -> dict:
    """
    Estimate reconstruction quality metrics directly from PLY file properties.

    Without rendering, we can still characterize the Gaussian cloud:
    - Number of Gaussians (model complexity)
    - Opacity distribution (well-converged Gaussians have high opacity)
    - Scale distribution (very small/large scales indicate artifacts)
    - Spatial extent and density
    """
    sys.path.insert(0, str(SCRIPT_DIR))
    from sigma_analysis import parse_ply

    cloud = parse_ply(ply_path)

    metrics = {
        "num_gaussians": cloud.num_gaussians,
        "file_size_mb": ply_path.stat().st_size / (1024 * 1024),
    }

    # Position statistics
    pos = cloud.positions
    metrics["spatial_extent"] = {
        "min": pos.min(axis=0).tolist(),
        "max": pos.max(axis=0).tolist(),
        "center": pos.mean(axis=0).tolist(),
        "std": pos.std(axis=0).tolist(),
    }

    # Opacity statistics (sigmoid-activated in log space)
    if cloud.opacities is not None:
        # OpenSplat stores pre-sigmoid opacity; apply sigmoid
        activated = 1.0 / (1.0 + np.exp(-cloud.opacities))
        metrics["opacity"] = {
            "mean": float(activated.mean()),
            "median": float(np.median(activated)),
            "std": float(activated.std()),
            "high_opacity_fraction": float((activated > 0.5).mean()),
            "very_low_fraction": float((activated < 0.01).mean()),
        }

    # Scale statistics
    if cloud.scales is not None:
        # OpenSplat stores log-scale; exponentiate
        scales_exp = np.exp(cloud.scales)
        metrics["scale"] = {
            "mean": float(scales_exp.mean()),
            "median": float(np.median(scales_exp)),
            "std": float(scales_exp.std()),
            "min": float(scales_exp.min()),
            "max": float(scales_exp.max()),
        }

    return metrics


# ---------------------------------------------------------------------------
# Comparison and Visualization
# ---------------------------------------------------------------------------
def generate_comparison_image(
    pairs: list[tuple[Path, Path]],
    output_path: Path,
    max_pairs: int = 4,
):
    """Generate side-by-side comparison images."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("[WARN] matplotlib not available. Skipping comparison image.")
        return

    pairs = pairs[:max_pairs]
    if not pairs:
        print("[WARN] No image pairs to compare.")
        return

    n = len(pairs)
    fig, axes = plt.subplots(n, 2, figsize=(14, 5 * n))
    if n == 1:
        axes = [axes]

    for i, (baseline_path, guided_path) in enumerate(pairs):
        img_baseline = load_image(baseline_path)
        img_guided = load_image(guided_path)

        if img_baseline is None or img_guided is None:
            continue

        axes[i][0].imshow(img_baseline)
        axes[i][0].set_title(f"Baseline (30k iter)\n{baseline_path.name}", fontsize=10)
        axes[i][0].axis("off")

        axes[i][1].imshow(img_guided)
        axes[i][1].set_title(f"Sigma-Guided (early stop)\n{guided_path.name}", fontsize=10)
        axes[i][1].axis("off")

    plt.suptitle("Line 3: Baseline vs Sigma-Guided 3DGS", fontsize=14, fontweight="bold")
    plt.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[INFO] Comparison image saved to: {output_path}")


def generate_ply_comparison_figure(
    baseline_metrics: dict,
    guided_metrics: dict,
    output_path: Path,
):
    """Generate comparison figure based on PLY analysis."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("[WARN] matplotlib not available. Skipping PLY comparison figure.")
        return

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # Panel 1: Gaussian count comparison
    labels = ["Baseline", "Sigma-Guided"]
    counts = [
        baseline_metrics.get("num_gaussians", 0),
        guided_metrics.get("num_gaussians", 0),
    ]
    colors_bar = ["#3498db", "#2ecc71"]
    axes[0].bar(labels, counts, color=colors_bar)
    axes[0].set_title("Number of Gaussians")
    axes[0].set_ylabel("Count")
    for i, v in enumerate(counts):
        axes[0].text(i, v + v * 0.01, f"{v:,}", ha="center", va="bottom", fontsize=10)

    # Panel 2: Opacity distribution comparison
    b_opacity = baseline_metrics.get("opacity", {})
    g_opacity = guided_metrics.get("opacity", {})
    if b_opacity and g_opacity:
        metrics_names = ["mean", "median", "high_opacity_fraction"]
        x = np.arange(len(metrics_names))
        width = 0.35
        b_vals = [b_opacity.get(m, 0) for m in metrics_names]
        g_vals = [g_opacity.get(m, 0) for m in metrics_names]
        axes[1].bar(x - width / 2, b_vals, width, label="Baseline", color="#3498db")
        axes[1].bar(x + width / 2, g_vals, width, label="Guided", color="#2ecc71")
        axes[1].set_xticks(x)
        axes[1].set_xticklabels(metrics_names, rotation=15)
        axes[1].set_title("Opacity Metrics")
        axes[1].legend()
    else:
        axes[1].text(0.5, 0.5, "No opacity data", ha="center", va="center",
                    transform=axes[1].transAxes)

    # Panel 3: File size comparison
    sizes = [
        baseline_metrics.get("file_size_mb", 0),
        guided_metrics.get("file_size_mb", 0),
    ]
    axes[2].bar(labels, sizes, color=colors_bar)
    axes[2].set_title("Model Size (MB)")
    axes[2].set_ylabel("MB")
    for i, v in enumerate(sizes):
        axes[2].text(i, v + v * 0.01, f"{v:.1f}", ha="center", va="bottom", fontsize=10)
    if sizes[0] > 0:
        reduction = (1 - sizes[1] / sizes[0]) * 100
        axes[2].set_xlabel(f"Size reduction: {reduction:.1f}%")

    plt.suptitle("Line 3: PLY Model Comparison", fontsize=14, fontweight="bold")
    plt.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[INFO] PLY comparison figure saved to: {output_path}")


# ---------------------------------------------------------------------------
# Main Evaluation Pipeline
# ---------------------------------------------------------------------------
def run_evaluation(
    baseline_dir: Path,
    guided_dir: Path,
    output_dir: Path,
    data_dir: Optional[Path] = None,
):
    """Run full evaluation pipeline."""
    print(f"\n{'='*72}")
    print(f"  Line 3: Evaluation - Baseline vs Sigma-Guided 3DGS")
    print(f"{'='*72}")
    print(f"  Baseline dir:  {baseline_dir}")
    print(f"  Guided dir:    {guided_dir}")
    print(f"  Output dir:    {output_dir}")
    print(f"{'='*72}\n")

    output_dir.mkdir(parents=True, exist_ok=True)

    metrics = {
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "baseline": {},
        "guided": {},
        "image_comparisons": [],
        "ply_comparisons": {},
    }

    # --- Part 1: PLY-based comparison ---
    print("[STEP 1] PLY-based model comparison...")

    # Find final PLY files
    baseline_plys = sorted(baseline_dir.glob("*.ply"))
    guided_plys = sorted((guided_dir / "checkpoints").glob("*.ply")) if (guided_dir / "checkpoints").exists() else []

    if not guided_plys:
        guided_plys = sorted(guided_dir.glob("*.ply"))

    baseline_final = baseline_plys[-1] if baseline_plys else None
    guided_final = guided_plys[-1] if guided_plys else None

    if baseline_final:
        print(f"  Baseline PLY: {baseline_final.name}")
        baseline_ply_metrics = estimate_quality_from_ply(baseline_final)
        metrics["baseline"]["ply_metrics"] = baseline_ply_metrics
    else:
        print("  [WARN] No baseline PLY found.")
        baseline_ply_metrics = {}

    if guided_final:
        print(f"  Guided PLY:   {guided_final.name}")
        guided_ply_metrics = estimate_quality_from_ply(guided_final)
        metrics["guided"]["ply_metrics"] = guided_ply_metrics
    else:
        print("  [WARN] No guided PLY found.")
        guided_ply_metrics = {}

    if baseline_ply_metrics and guided_ply_metrics:
        metrics["ply_comparisons"] = {
            "gaussian_count_baseline": baseline_ply_metrics.get("num_gaussians", 0),
            "gaussian_count_guided": guided_ply_metrics.get("num_gaussians", 0),
            "size_mb_baseline": baseline_ply_metrics.get("file_size_mb", 0),
            "size_mb_guided": guided_ply_metrics.get("file_size_mb", 0),
        }

        # Generate PLY comparison figure
        generate_ply_comparison_figure(
            baseline_ply_metrics, guided_ply_metrics,
            output_dir / "comparison_ply.png"
        )

    # --- Part 2: Rendered image comparison (if available) ---
    print("\n[STEP 2] Rendered image comparison...")

    baseline_images = find_rendered_images(baseline_dir)
    guided_images = find_rendered_images(guided_dir)

    if baseline_images and guided_images:
        pairs = find_matching_pairs(baseline_images, guided_images)
        print(f"  Found {len(pairs)} matching image pairs.")

        # Compute PSNR and SSIM for each pair
        for baseline_path, guided_path in pairs:
            img_baseline = load_image(baseline_path)
            img_guided = load_image(guided_path)

            if img_baseline is None or img_guided is None:
                continue

            # Ensure same shape
            if img_baseline.shape != img_guided.shape:
                # Resize to match
                min_h = min(img_baseline.shape[0], img_guided.shape[0])
                min_w = min(img_baseline.shape[1], img_guided.shape[1])
                img_baseline = img_baseline[:min_h, :min_w]
                img_guided = img_guided[:min_h, :min_w]

            psnr = compute_psnr(img_baseline, img_guided)

            try:
                ssim = compute_ssim(img_baseline, img_guided)
            except ImportError:
                ssim = _compute_ssim_fallback(img_baseline, img_guided)

            comp = {
                "name": baseline_path.stem,
                "psnr_db": psnr,
                "ssim": ssim,
            }
            metrics["image_comparisons"].append(comp)
            print(f"  {baseline_path.stem}: PSNR={psnr:.2f} dB, SSIM={ssim:.4f}")

        # If we also have ground truth
        if data_dir:
            gt_images = find_ground_truth_images(data_dir)
            if gt_images:
                print(f"\n  Ground truth images found: {len(gt_images)}")
                metrics["baseline"]["vs_ground_truth"] = []
                metrics["guided"]["vs_ground_truth"] = []

                for gt_path in gt_images[:5]:  # limit to first 5
                    gt_img = load_image(gt_path)
                    if gt_img is None:
                        continue

                    # Find matching rendered images
                    for rendered_dir, key in [(baseline_dir, "baseline"), (guided_dir, "guided")]:
                        rendered_path = rendered_dir / gt_path.name
                        if not rendered_path.exists():
                            continue
                        rendered_img = load_image(rendered_path)
                        if rendered_img is None:
                            continue
                        if gt_img.shape != rendered_img.shape:
                            min_h = min(gt_img.shape[0], rendered_img.shape[0])
                            min_w = min(gt_img.shape[1], rendered_img.shape[1])
                            gt_crop = gt_img[:min_h, :min_w]
                            rendered_crop = rendered_img[:min_h, :min_w]
                        else:
                            gt_crop = gt_img
                            rendered_crop = rendered_img

                        psnr_gt = compute_psnr(gt_crop, rendered_crop)
                        try:
                            ssim_gt = compute_ssim(gt_crop, rendered_crop)
                        except ImportError:
                            ssim_gt = _compute_ssim_fallback(gt_crop, rendered_crop)

                        metrics[key]["vs_ground_truth"].append({
                            "name": gt_path.stem,
                            "psnr_db": psnr_gt,
                            "ssim": ssim_gt,
                        })
                        print(f"  {key} vs GT ({gt_path.stem}): "
                              f"PSNR={psnr_gt:.2f} dB, SSIM={ssim_gt:.4f}")

        # Generate comparison image
        if pairs:
            generate_comparison_image(pairs, output_dir / "comparison.png")
    else:
        print("  No rendered images found for comparison.")
        print("  [NOTE] OpenSplat renders are typically saved during training.")
        print("         PLY-based comparison is available above.")

    # --- Part 3: Training curve comparison ---
    print("\n[STEP 3] Training curve comparison...")

    baseline_meta_path = baseline_dir / "training_metadata.json"
    guided_info_path = guided_dir / "guided_info.json"

    if baseline_meta_path.exists() and guided_info_path.exists():
        with open(baseline_meta_path) as f:
            baseline_meta = json.load(f)
        with open(guided_info_path) as f:
            guided_info = json.load(f)

        b_losses = baseline_meta.get("losses", [])
        g_losses = guided_info.get("training_metadata", {}).get("losses", [])

        if b_losses and g_losses:
            b_final = b_losses[-1]["loss"]
            g_final = g_losses[-1]["loss"]
            loss_diff = abs(b_final - g_final)
            loss_diff_pct = (loss_diff / b_final) * 100 if b_final > 0 else 0

            metrics["loss_comparison"] = {
                "baseline_final_loss": b_final,
                "guided_final_loss": g_final,
                "absolute_difference": loss_diff,
                "relative_difference_pct": loss_diff_pct,
            }

            print(f"  Baseline final loss:  {b_final:.6f}")
            print(f"  Guided final loss:    {g_final:.6f}")
            print(f"  Difference:           {loss_diff:.6f} ({loss_diff_pct:.2f}%)")

        # Time comparison
        b_time = baseline_meta.get("elapsed_seconds", 0)
        g_time = guided_info.get("training_metadata", {}).get("elapsed_seconds", 0)
        if b_time > 0 and g_time > 0:
            speedup = b_time / g_time
            time_savings = (1 - g_time / b_time) * 100
            metrics["time_comparison"] = {
                "baseline_seconds": b_time,
                "guided_seconds": g_time,
                "speedup_factor": speedup,
                "time_savings_pct": time_savings,
            }
            print(f"\n  Baseline time:  {b_time:.1f}s")
            print(f"  Guided time:    {g_time:.1f}s")
            print(f"  Speedup:        {speedup:.2f}x ({time_savings:.1f}% savings)")

    # --- Summary ---
    print(f"\n{'='*72}")
    print(f"  Evaluation Summary")
    print(f"{'='*72}")

    if metrics["image_comparisons"]:
        avg_psnr = np.mean([c["psnr_db"] for c in metrics["image_comparisons"]
                           if not math.isinf(c["psnr_db"])])
        avg_ssim = np.mean([c["ssim"] for c in metrics["image_comparisons"]])
        metrics["summary"] = {
            "avg_psnr_db": float(avg_psnr),
            "avg_ssim": float(avg_ssim),
        }
        print(f"  Average PSNR: {avg_psnr:.2f} dB")
        print(f"  Average SSIM: {avg_ssim:.4f}")
    elif metrics["ply_comparisons"]:
        print(f"  Gaussian count: baseline={metrics['ply_comparisons']['gaussian_count_baseline']}, "
              f"guided={metrics['ply_comparisons']['gaussian_count_guided']}")
        print(f"  Model size: baseline={metrics['ply_comparisons']['size_mb_baseline']:.1f} MB, "
              f"guided={metrics['ply_comparisons']['size_mb_guided']:.1f} MB")
    else:
        print("  No quantitative comparison available yet.")
        print("  Run baseline and guided training first.")

    print(f"{'='*72}")

    # Save metrics
    metrics_path = output_dir / "metrics.json"
    with open(metrics_path, "w") as f:
        json.dump(metrics, f, indent=2)
    print(f"\n[INFO] Metrics saved to: {metrics_path}")

    return metrics


def main():
    parser = argparse.ArgumentParser(
        description="Evaluate baseline vs sigma-guided 3DGS"
    )
    parser.add_argument(
        "--baseline-dir", type=str, default=str(CHECKPOINT_DIR),
        help="Directory with baseline checkpoints/renders"
    )
    parser.add_argument(
        "--guided-dir", type=str, default=str(GUIDED_DIR),
        help="Directory with guided checkpoints/renders"
    )
    parser.add_argument(
        "--output-dir", type=str, default=str(RESULTS_DIR),
        help="Output directory for evaluation results"
    )
    parser.add_argument(
        "--data-dir", type=str, default=None,
        help="Original scene data directory (for ground truth)"
    )
    args = parser.parse_args()

    # Find data dir if not specified
    data_dir = Path(args.data_dir) if args.data_dir else None
    if data_dir is None:
        candidates = [
            SCRIPT_DIR / "data" / "minimal_test",
            SCRIPT_DIR / "data" / "garden",
            SCRIPT_DIR / "data" / "lego",
        ]
        for c in candidates:
            if c.is_dir():
                data_dir = c
                break

    run_evaluation(
        baseline_dir=Path(args.baseline_dir),
        guided_dir=Path(args.guided_dir),
        output_dir=Path(args.output_dir),
        data_dir=data_dir,
    )


if __name__ == "__main__":
    main()
