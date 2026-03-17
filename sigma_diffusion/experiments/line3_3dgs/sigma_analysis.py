#!/usr/bin/env python3
"""
Line 3 Experiment: Sigma Analysis of 3D Gaussian Splatting Checkpoints.

Loads consecutive PLY checkpoint files from OpenSplat, computes per-iteration
Sigma (parameter drift), and identifies plateau regions where optimization
has converged.

Core idea: Sigma_t = D(rho_t || rho_{t-1}) measures the information-theoretic
"cost" of one optimization step. When Sigma drops to near-zero, the model has
reached a local equilibrium -- further training yields diminishing returns.

Usage:
    python sigma_analysis.py
    python sigma_analysis.py --checkpoint-dir /path/to/checkpoints
    python sigma_analysis.py --threshold 0.001  # Custom plateau threshold
"""

import argparse
import json
import struct
import sys
import time
from dataclasses import dataclass, field
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


# ---------------------------------------------------------------------------
# PLY Parsing (handles OpenSplat output format)
# ---------------------------------------------------------------------------
@dataclass
class GaussianCloud:
    """Represents a 3D Gaussian Splatting point cloud."""
    positions: np.ndarray       # (N, 3) xyz
    colors_dc: np.ndarray       # (N, 3) DC component of SH (base color)
    colors_sh: Optional[np.ndarray] = None  # (N, K) higher-order SH coefficients
    opacities: np.ndarray = None  # (N,) sigmoid-activated opacity
    scales: np.ndarray = None     # (N, 3) log-scale
    rotations: np.ndarray = None  # (N, 4) quaternion
    num_gaussians: int = 0

    def parameter_vector(self) -> np.ndarray:
        """Concatenate all parameters into a single vector per Gaussian.

        Returns shape (N, D) where D = 3 + 3 + 1 + 3 + 4 = 14 (minimum).
        """
        components = [self.positions, self.colors_dc]
        if self.opacities is not None:
            components.append(self.opacities.reshape(-1, 1))
        if self.scales is not None:
            components.append(self.scales)
        if self.rotations is not None:
            components.append(self.rotations)
        return np.concatenate(components, axis=1)


def parse_ply(filepath: Path) -> GaussianCloud:
    """Parse a PLY file from OpenSplat into a GaussianCloud.

    Handles both ASCII and binary_little_endian PLY formats.
    OpenSplat PLY files typically contain:
      - x, y, z (positions)
      - f_dc_0, f_dc_1, f_dc_2 (DC color)
      - f_rest_0 ... f_rest_N (SH coefficients)
      - opacity
      - scale_0, scale_1, scale_2
      - rot_0, rot_1, rot_2, rot_3
    """
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"PLY file not found: {filepath}")

    file_size_mb = filepath.stat().st_size / (1024 * 1024)
    print(f"  Parsing {filepath.name} ({file_size_mb:.1f} MB)...", end=" ", flush=True)
    t0 = time.time()

    with open(filepath, "rb") as f:
        # Parse header
        header_lines = []
        while True:
            line = f.readline()
            if not line:
                raise ValueError("Unexpected end of file in PLY header")
            line_str = line.decode("ascii", errors="replace").strip()
            header_lines.append(line_str)
            if line_str == "end_header":
                break

        # Extract metadata from header
        format_line = [l for l in header_lines if l.startswith("format")]
        if not format_line:
            raise ValueError("No format line in PLY header")

        is_binary = "binary_little_endian" in format_line[0]
        is_ascii = "ascii" in format_line[0]

        # Get vertex count
        vertex_count = 0
        in_vertex_element = False
        properties = []
        for line_str in header_lines:
            if line_str.startswith("element vertex"):
                vertex_count = int(line_str.split()[-1])
                in_vertex_element = True
            elif line_str.startswith("element") and in_vertex_element:
                in_vertex_element = False
            elif line_str.startswith("property") and in_vertex_element:
                parts = line_str.split()
                prop_type = parts[1]
                prop_name = parts[2]
                properties.append((prop_name, prop_type))

        if vertex_count == 0:
            raise ValueError("No vertices found in PLY file")

        # Map property names to indices
        prop_indices = {name: i for i, (name, _) in enumerate(properties)}

        # Determine struct format for binary
        type_map = {
            "float": "f", "float32": "f",
            "double": "d", "float64": "d",
            "uchar": "B", "uint8": "B",
            "char": "b", "int8": "b",
            "short": "h", "int16": "h",
            "ushort": "H", "uint16": "H",
            "int": "i", "int32": "i",
            "uint": "I", "uint32": "I",
        }

        type_sizes = {
            "f": 4, "d": 8, "B": 1, "b": 1,
            "h": 2, "H": 2, "i": 4, "I": 4,
        }

        struct_fmt = "<" + "".join(type_map.get(t, "f") for _, t in properties)
        vertex_size = struct.calcsize(struct_fmt)

        # Read vertex data
        if is_binary:
            raw_data = f.read(vertex_count * vertex_size)
            if len(raw_data) < vertex_count * vertex_size:
                raise ValueError(
                    f"Incomplete PLY data: expected {vertex_count * vertex_size} bytes, "
                    f"got {len(raw_data)}"
                )

            # Parse all vertices at once using numpy for speed
            dt = np.dtype([(name, f"<{type_map.get(t, 'f')}") for name, t in properties])
            vertices = np.frombuffer(raw_data, dtype=dt, count=vertex_count)

        elif is_ascii:
            vertices_list = []
            for _ in range(vertex_count):
                line = f.readline().decode("ascii").strip()
                values = line.split()
                vertices_list.append([float(v) for v in values])
            vertices_arr = np.array(vertices_list, dtype=np.float32)
            # Create structured array
            dt = np.dtype([(name, np.float32) for name, _ in properties])
            vertices = np.zeros(vertex_count, dtype=dt)
            for i, (name, _) in enumerate(properties):
                vertices[name] = vertices_arr[:, i]
        else:
            raise ValueError(f"Unsupported PLY format: {format_line[0]}")

    # Extract fields
    cloud = GaussianCloud(
        positions=np.column_stack([
            vertices["x"].astype(np.float32),
            vertices["y"].astype(np.float32),
            vertices["z"].astype(np.float32),
        ]),
        colors_dc=np.zeros((vertex_count, 3), dtype=np.float32),
        num_gaussians=vertex_count,
    )

    # Colors (DC component of spherical harmonics)
    dc_names = ["f_dc_0", "f_dc_1", "f_dc_2"]
    rgb_names = ["red", "green", "blue"]

    if all(n in prop_indices for n in dc_names):
        cloud.colors_dc = np.column_stack([
            vertices[n].astype(np.float32) for n in dc_names
        ])
    elif all(n in prop_indices for n in rgb_names):
        cloud.colors_dc = np.column_stack([
            vertices[n].astype(np.float32) / 255.0 for n in rgb_names
        ])

    # Higher-order SH coefficients
    sh_names = [f"f_rest_{i}" for i in range(45)]  # up to degree 3
    sh_present = [n for n in sh_names if n in prop_indices]
    if sh_present:
        cloud.colors_sh = np.column_stack([
            vertices[n].astype(np.float32) for n in sh_present
        ])

    # Opacity
    if "opacity" in prop_indices:
        cloud.opacities = vertices["opacity"].astype(np.float32)

    # Scales
    scale_names = ["scale_0", "scale_1", "scale_2"]
    if all(n in prop_indices for n in scale_names):
        cloud.scales = np.column_stack([
            vertices[n].astype(np.float32) for n in scale_names
        ])

    # Rotations
    rot_names = ["rot_0", "rot_1", "rot_2", "rot_3"]
    if all(n in prop_indices for n in rot_names):
        cloud.rotations = np.column_stack([
            vertices[n].astype(np.float32) for n in rot_names
        ])

    elapsed = time.time() - t0
    print(f"done ({vertex_count} Gaussians, {elapsed:.2f}s)")
    return cloud


# ---------------------------------------------------------------------------
# Sigma Computation
# ---------------------------------------------------------------------------
@dataclass
class SigmaRecord:
    """Sigma measurement for a single transition between checkpoints."""
    iteration_from: int
    iteration_to: int
    sigma_mean: float           # Mean parameter drift across all Gaussians
    sigma_median: float         # Median (more robust to outliers)
    sigma_max: float            # Maximum drift (identifies fast-changing Gaussians)
    sigma_std: float            # Std dev of drift
    num_gaussians_from: int
    num_gaussians_to: int
    num_matched: int            # How many Gaussians were matched
    densification_delta: int    # Change in Gaussian count
    per_component: dict = field(default_factory=dict)  # Sigma broken down by component


def compute_sigma(
    cloud_prev: GaussianCloud,
    cloud_curr: GaussianCloud,
    iter_prev: int,
    iter_curr: int,
) -> SigmaRecord:
    """
    Compute Sigma between two consecutive Gaussian clouds.

    When Gaussian counts match, we compute element-wise drift.
    When counts differ (due to densification/pruning), we match by
    nearest-neighbor in position space and compute drift on matched pairs.
    """
    n_prev = cloud_prev.num_gaussians
    n_curr = cloud_curr.num_gaussians
    n_matched = min(n_prev, n_curr)

    # If counts match, assume 1:1 correspondence (standard case between
    # densification steps)
    if n_prev == n_curr:
        params_prev = cloud_prev.parameter_vector()
        params_curr = cloud_curr.parameter_vector()
    else:
        # Counts differ: match by nearest neighbor in position space
        # For efficiency with large clouds, use a simple truncation strategy:
        # match the first min(N_prev, N_curr) Gaussians by index
        # (OpenSplat preserves ordering for existing Gaussians)
        params_prev = cloud_prev.parameter_vector()[:n_matched]
        params_curr = cloud_curr.parameter_vector()[:n_matched]

    # Normalize each parameter dimension to [0, 1] range for fair comparison
    combined = np.concatenate([params_prev, params_curr], axis=0)
    col_min = combined.min(axis=0)
    col_max = combined.max(axis=0)
    col_range = col_max - col_min
    col_range[col_range < 1e-10] = 1.0  # avoid division by zero

    params_prev_norm = (params_prev - col_min) / col_range
    params_curr_norm = (params_curr - col_min) / col_range

    # Per-Gaussian drift: ||params_t - params_{t-1}||^2
    drift = np.sum((params_curr_norm - params_prev_norm) ** 2, axis=1)

    # Per-component Sigma breakdown
    per_component = {}
    dim_offset = 0
    component_dims = {
        "position": 3,
        "color_dc": 3,
    }
    if cloud_curr.opacities is not None:
        component_dims["opacity"] = 1
    if cloud_curr.scales is not None:
        component_dims["scale"] = 3
    if cloud_curr.rotations is not None:
        component_dims["rotation"] = 4

    for comp_name, comp_dim in component_dims.items():
        end = dim_offset + comp_dim
        if end <= params_prev_norm.shape[1]:
            comp_drift = np.sum(
                (params_curr_norm[:, dim_offset:end] - params_prev_norm[:, dim_offset:end]) ** 2,
                axis=1
            )
            per_component[comp_name] = {
                "mean": float(np.mean(comp_drift)),
                "median": float(np.median(comp_drift)),
                "max": float(np.max(comp_drift)),
            }
        dim_offset = end

    return SigmaRecord(
        iteration_from=iter_prev,
        iteration_to=iter_curr,
        sigma_mean=float(np.mean(drift)),
        sigma_median=float(np.median(drift)),
        sigma_max=float(np.max(drift)),
        sigma_std=float(np.std(drift)),
        num_gaussians_from=n_prev,
        num_gaussians_to=n_curr,
        num_matched=n_matched,
        densification_delta=n_curr - n_prev,
        per_component=per_component,
    )


# ---------------------------------------------------------------------------
# Analysis Pipeline
# ---------------------------------------------------------------------------
def find_checkpoints(checkpoint_dir: Path) -> list[tuple[int, Path]]:
    """Find and sort checkpoint PLY files by iteration number."""
    checkpoints = []
    for ply_file in sorted(checkpoint_dir.glob("*.ply")):
        # Try to extract iteration number from filename
        name = ply_file.stem
        # Patterns: checkpoint_000500.ply, iter_500.ply, 500.ply, etc.
        import re
        match = re.search(r"(\d+)", name)
        if match:
            iteration = int(match.group(1))
            checkpoints.append((iteration, ply_file))

    checkpoints.sort(key=lambda x: x[0])
    return checkpoints


def find_plateau_regions(
    sigma_records: list[SigmaRecord],
    threshold: float,
    min_consecutive: int = 3,
) -> list[dict]:
    """
    Identify plateau regions where Sigma stays below threshold
    for at least min_consecutive steps.
    """
    plateaus = []
    current_plateau = None

    for record in sigma_records:
        if record.sigma_mean < threshold:
            if current_plateau is None:
                current_plateau = {
                    "start_iter": record.iteration_from,
                    "end_iter": record.iteration_to,
                    "length": 1,
                    "mean_sigma": record.sigma_mean,
                    "records": [record],
                }
            else:
                current_plateau["end_iter"] = record.iteration_to
                current_plateau["length"] += 1
                current_plateau["records"].append(record)
                current_plateau["mean_sigma"] = np.mean(
                    [r.sigma_mean for r in current_plateau["records"]]
                )
        else:
            if current_plateau and current_plateau["length"] >= min_consecutive:
                plateaus.append(current_plateau)
            current_plateau = None

    # Don't forget the last plateau
    if current_plateau and current_plateau["length"] >= min_consecutive:
        plateaus.append(current_plateau)

    return plateaus


def generate_sigma_profile(
    sigma_records: list[SigmaRecord],
    losses: list[dict],
    plateaus: list[dict],
    output_dir: Path,
):
    """Generate sigma_profile.json with all analysis results."""
    profile = {
        "analysis_timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "num_checkpoints": len(sigma_records) + 1,
        "num_transitions": len(sigma_records),
        "sigma_records": [],
        "plateaus": [],
        "losses": losses,
        "summary": {},
    }

    for record in sigma_records:
        profile["sigma_records"].append({
            "iteration_from": record.iteration_from,
            "iteration_to": record.iteration_to,
            "sigma_mean": record.sigma_mean,
            "sigma_median": record.sigma_median,
            "sigma_max": record.sigma_max,
            "sigma_std": record.sigma_std,
            "num_gaussians_from": record.num_gaussians_from,
            "num_gaussians_to": record.num_gaussians_to,
            "num_matched": record.num_matched,
            "densification_delta": record.densification_delta,
            "per_component": record.per_component,
        })

    for plateau in plateaus:
        profile["plateaus"].append({
            "start_iter": plateau["start_iter"],
            "end_iter": plateau["end_iter"],
            "length": plateau["length"],
            "mean_sigma": float(plateau["mean_sigma"]),
        })

    # Summary statistics
    if sigma_records:
        sigmas = [r.sigma_mean for r in sigma_records]
        profile["summary"] = {
            "sigma_mean_overall": float(np.mean(sigmas)),
            "sigma_max_overall": float(np.max(sigmas)),
            "sigma_min_overall": float(np.min(sigmas)),
            "sigma_final": sigmas[-1],
            "total_gaussians_initial": sigma_records[0].num_gaussians_from,
            "total_gaussians_final": sigma_records[-1].num_gaussians_to,
            "num_plateaus": len(plateaus),
            "first_plateau_iter": plateaus[0]["start_iter"] if plateaus else None,
            "suggested_early_stop": plateaus[0]["start_iter"] if plateaus else None,
        }

    output_path = output_dir / "sigma_profile.json"
    with open(output_path, "w") as f:
        json.dump(profile, f, indent=2)

    print(f"\n[INFO] Sigma profile saved to: {output_path}")
    return profile


def generate_sigma_heatmap(
    sigma_records: list[SigmaRecord],
    losses: list[dict],
    plateaus: list[dict],
    output_dir: Path,
):
    """Generate sigma_heatmap.png visualization."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.gridspec as gridspec
    except ImportError:
        print("[WARN] matplotlib not available. Skipping heatmap generation.")
        print("       Install with: pip install matplotlib")
        return

    iterations = [(r.iteration_from + r.iteration_to) / 2 for r in sigma_records]
    sigma_means = [r.sigma_mean for r in sigma_records]
    sigma_medians = [r.sigma_median for r in sigma_records]
    sigma_maxs = [r.sigma_max for r in sigma_records]
    n_gaussians = [r.num_gaussians_to for r in sigma_records]

    fig = plt.figure(figsize=(16, 14))
    gs = gridspec.GridSpec(4, 1, height_ratios=[3, 2, 2, 1.5], hspace=0.3)

    # --- Panel 1: Sigma over iterations ---
    ax1 = fig.add_subplot(gs[0])
    ax1.semilogy(iterations, sigma_means, "b-", linewidth=1.5, label="Sigma (mean)", alpha=0.9)
    ax1.semilogy(iterations, sigma_medians, "g--", linewidth=1.0, label="Sigma (median)", alpha=0.7)
    ax1.fill_between(iterations, sigma_medians, sigma_maxs, alpha=0.15, color="blue",
                     label="median-to-max range")

    # Mark plateaus
    for plateau in plateaus:
        ax1.axvspan(plateau["start_iter"], plateau["end_iter"],
                    alpha=0.2, color="green", label="Plateau" if plateau == plateaus[0] else None)

    ax1.set_xlabel("Iteration")
    ax1.set_ylabel("Sigma (log scale)")
    ax1.set_title("Line 3: Sigma Profile for 3D Gaussian Splatting", fontsize=14, fontweight="bold")
    ax1.legend(loc="upper right")
    ax1.grid(True, alpha=0.3)

    # --- Panel 2: Per-component Sigma breakdown ---
    ax2 = fig.add_subplot(gs[1])
    component_names = ["position", "color_dc", "opacity", "scale", "rotation"]
    colors = ["#e74c3c", "#2ecc71", "#3498db", "#f39c12", "#9b59b6"]

    for comp_name, color in zip(component_names, colors):
        comp_sigmas = []
        for r in sigma_records:
            if comp_name in r.per_component:
                comp_sigmas.append(r.per_component[comp_name]["mean"])
            else:
                comp_sigmas.append(0)
        if any(s > 0 for s in comp_sigmas):
            ax2.semilogy(iterations, comp_sigmas, linewidth=1.2, color=color,
                        label=comp_name, alpha=0.8)

    ax2.set_xlabel("Iteration")
    ax2.set_ylabel("Component Sigma (log)")
    ax2.set_title("Per-Component Sigma Breakdown")
    ax2.legend(loc="upper right", ncol=3)
    ax2.grid(True, alpha=0.3)

    # --- Panel 3: Training loss (if available) ---
    ax3 = fig.add_subplot(gs[2])
    if losses:
        loss_iters = [l["iteration"] for l in losses]
        loss_vals = [l["loss"] for l in losses]
        ax3.semilogy(loss_iters, loss_vals, "r-", linewidth=1.0, alpha=0.7, label="Training Loss")
        ax3.set_xlabel("Iteration")
        ax3.set_ylabel("Loss (log scale)")
        ax3.set_title("Training Loss")
        ax3.legend()
        ax3.grid(True, alpha=0.3)
    else:
        ax3.text(0.5, 0.5, "No loss data available\n(run run_baseline.py first)",
                ha="center", va="center", transform=ax3.transAxes, fontsize=12,
                color="gray")
        ax3.set_title("Training Loss (not available)")

    # --- Panel 4: Number of Gaussians ---
    ax4 = fig.add_subplot(gs[3])
    ax4.plot(iterations, n_gaussians, "k-", linewidth=1.5)
    ax4.set_xlabel("Iteration")
    ax4.set_ylabel("# Gaussians")
    ax4.set_title("Gaussian Count (densification events)")
    ax4.grid(True, alpha=0.3)

    # Mark densification events
    for i, r in enumerate(sigma_records):
        if abs(r.densification_delta) > 0:
            ax4.axvline(x=iterations[i], color="red", alpha=0.3, linewidth=0.5)

    plt.tight_layout()
    output_path = output_dir / "sigma_heatmap.png"
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[INFO] Sigma heatmap saved to: {output_path}")


def run_analysis(
    checkpoint_dir: Path,
    output_dir: Path,
    threshold: float = 0.001,
    min_plateau_length: int = 3,
):
    """Main analysis pipeline."""
    print(f"\n{'='*72}")
    print(f"  Line 3: Sigma Analysis of 3D Gaussian Splatting")
    print(f"{'='*72}")
    print(f"  Checkpoint directory: {checkpoint_dir}")
    print(f"  Output directory:     {output_dir}")
    print(f"  Plateau threshold:    {threshold}")
    print(f"{'='*72}\n")

    # Find checkpoints
    checkpoints = find_checkpoints(checkpoint_dir)
    if len(checkpoints) < 2:
        print(f"[ERROR] Need at least 2 checkpoints, found {len(checkpoints)}")
        print(f"        Run run_baseline.py first to generate checkpoints.")
        sys.exit(1)

    print(f"[INFO] Found {len(checkpoints)} checkpoints:")
    for iteration, path in checkpoints:
        size_mb = path.stat().st_size / (1024 * 1024)
        print(f"  iter {iteration:>6d}: {path.name} ({size_mb:.1f} MB)")

    # Load training metadata (losses) if available
    meta_path = checkpoint_dir / "training_metadata.json"
    losses = []
    if meta_path.exists():
        with open(meta_path) as f:
            meta = json.load(f)
            losses = meta.get("losses", [])
        print(f"\n[INFO] Loaded {len(losses)} loss records from training metadata.")

    # Compute Sigma for consecutive pairs
    print(f"\n[INFO] Computing Sigma between consecutive checkpoints...")
    sigma_records = []
    prev_cloud = None
    prev_iter = None

    for iteration, ply_path in checkpoints:
        cloud = parse_ply(ply_path)

        if prev_cloud is not None:
            record = compute_sigma(prev_cloud, cloud, prev_iter, iteration)
            sigma_records.append(record)

            delta_str = ""
            if record.densification_delta != 0:
                delta_str = f"  [densification: {record.densification_delta:+d}]"
            print(f"  {prev_iter:>6d} -> {iteration:>6d}: "
                  f"Sigma={record.sigma_mean:.6f} "
                  f"(median={record.sigma_median:.6f}, max={record.sigma_max:.4f})"
                  f"{delta_str}")

        prev_cloud = cloud
        prev_iter = iteration

    if not sigma_records:
        print("[ERROR] No Sigma records computed. Need more checkpoints.")
        sys.exit(1)

    # Find plateaus
    plateaus = find_plateau_regions(sigma_records, threshold, min_plateau_length)

    # Print summary
    print(f"\n{'='*72}")
    print(f"  Analysis Summary")
    print(f"{'='*72}")
    sigmas = [r.sigma_mean for r in sigma_records]
    print(f"  Sigma range:       [{min(sigmas):.6f}, {max(sigmas):.6f}]")
    print(f"  Sigma mean:        {np.mean(sigmas):.6f}")
    print(f"  Initial Gaussians: {sigma_records[0].num_gaussians_from}")
    print(f"  Final Gaussians:   {sigma_records[-1].num_gaussians_to}")
    print(f"  Plateaus found:    {len(plateaus)}")

    if plateaus:
        print(f"\n  Plateau regions (Sigma < {threshold}):")
        for i, p in enumerate(plateaus):
            print(f"    [{i+1}] iter {p['start_iter']}-{p['end_iter']} "
                  f"(length={p['length']}, mean_sigma={p['mean_sigma']:.6f})")

        suggested_stop = plateaus[0]["start_iter"]
        print(f"\n  ** Suggested early stop: iteration {suggested_stop} **")
        total_iters = checkpoints[-1][0]
        savings = (1 - suggested_stop / total_iters) * 100
        print(f"  ** Potential savings: {savings:.1f}% of training iterations **")
    else:
        print(f"\n  No plateaus found with threshold={threshold}.")
        print(f"  Try a higher threshold or more checkpoints.")

    print(f"{'='*72}")

    # Generate outputs
    output_dir.mkdir(parents=True, exist_ok=True)
    profile = generate_sigma_profile(sigma_records, losses, plateaus, output_dir)
    generate_sigma_heatmap(sigma_records, losses, plateaus, output_dir)

    return profile


def main():
    parser = argparse.ArgumentParser(
        description="Sigma analysis of 3DGS checkpoints"
    )
    parser.add_argument(
        "--checkpoint-dir", type=str, default=str(CHECKPOINT_DIR),
        help="Directory containing PLY checkpoint files"
    )
    parser.add_argument(
        "--output-dir", type=str, default=str(RESULTS_DIR),
        help="Directory for analysis output"
    )
    parser.add_argument(
        "--threshold", type=float, default=0.001,
        help="Sigma threshold for plateau detection (default: 0.001)"
    )
    parser.add_argument(
        "--min-plateau-length", type=int, default=3,
        help="Minimum consecutive steps below threshold to count as plateau (default: 3)"
    )
    args = parser.parse_args()

    run_analysis(
        checkpoint_dir=Path(args.checkpoint_dir),
        output_dir=Path(args.output_dir),
        threshold=args.threshold,
        min_plateau_length=args.min_plateau_length,
    )


if __name__ == "__main__":
    main()
