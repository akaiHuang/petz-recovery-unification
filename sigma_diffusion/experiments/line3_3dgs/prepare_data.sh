#!/usr/bin/env bash
# =============================================================================
# Line 3 Experiment: Data Preparation for 3D Gaussian Splatting
# =============================================================================
# This script prepares a COLMAP-formatted dataset for OpenSplat.
#
# Usage:
#   ./prepare_data.sh                    # Download pre-made synthetic dataset
#   ./prepare_data.sh --from-images DIR  # Run COLMAP on your own images
#
# Requirements (Mac):
#   brew install colmap  (only needed for --from-images mode)
#   pip install gdown    (only needed for Google Drive downloads)
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
DATA_DIR="$SCRIPT_DIR/data"
RESULTS_DIR="$PROJECT_ROOT/sigma_diffusion/results/line3"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

info()  { echo -e "${GREEN}[INFO]${NC} $*"; }
warn()  { echo -e "${YELLOW}[WARN]${NC} $*"; }
error() { echo -e "${RED}[ERROR]${NC} $*"; }

# ---------------------------------------------------------------------------
# Option 1: Download the Mip-NeRF 360 "garden" scene (small, well-tested)
# ---------------------------------------------------------------------------
download_mipnerf360_garden() {
    info "Downloading Mip-NeRF 360 'garden' scene (COLMAP pre-processed)..."
    info "This is a standard benchmark scene for 3DGS (~1.5 GB)."

    local DEST="$DATA_DIR/garden"
    if [ -d "$DEST/images" ] && [ -d "$DEST/sparse" ]; then
        info "Garden scene already exists at $DEST. Skipping download."
        return 0
    fi

    mkdir -p "$DEST"

    # Mip-NeRF 360 dataset hosted by the original authors
    local URL="http://storage.googleapis.com/gresearch/refraw360/garden.zip"

    info "Downloading from: $URL"
    info "This may take a few minutes depending on your connection..."

    if command -v wget &>/dev/null; then
        wget -q --show-progress -O "$DATA_DIR/garden.zip" "$URL"
    elif command -v curl &>/dev/null; then
        curl -L --progress-bar -o "$DATA_DIR/garden.zip" "$URL"
    else
        error "Neither wget nor curl found. Please install one of them."
        exit 1
    fi

    info "Extracting..."
    cd "$DATA_DIR"
    unzip -q -o garden.zip -d garden_tmp

    # The zip may have a nested directory structure; normalize it
    if [ -d "garden_tmp/garden" ]; then
        mv garden_tmp/garden/* garden/ 2>/dev/null || true
        rm -rf garden_tmp
    elif [ -d "garden_tmp/images" ]; then
        mv garden_tmp/* garden/ 2>/dev/null || true
        rm -rf garden_tmp
    else
        # Move whatever we got
        mv garden_tmp/* garden/ 2>/dev/null || true
        rm -rf garden_tmp
    fi

    rm -f garden.zip
    info "Garden scene ready at: $DEST"
}

# ---------------------------------------------------------------------------
# Option 2: Download a small synthetic NeRF scene (Blender "lego")
# ---------------------------------------------------------------------------
download_nerf_synthetic_lego() {
    info "Downloading NeRF synthetic 'lego' scene (~50 MB)..."

    local DEST="$DATA_DIR/lego"
    if [ -d "$DEST/images" ] || [ -d "$DEST/train" ]; then
        info "Lego scene already exists at $DEST. Skipping download."
        return 0
    fi

    mkdir -p "$DEST"

    # NeRF synthetic dataset (Blender scenes)
    local URL="https://huggingface.co/datasets/nrtf/nerf_synthetic/resolve/main/nerf_synthetic.zip"

    info "Downloading from HuggingFace: nerf_synthetic dataset..."
    if command -v curl &>/dev/null; then
        curl -L --progress-bar -o "$DATA_DIR/nerf_synthetic.zip" "$URL"
    elif command -v wget &>/dev/null; then
        wget -q --show-progress -O "$DATA_DIR/nerf_synthetic.zip" "$URL"
    else
        error "Neither wget nor curl found."
        exit 1
    fi

    info "Extracting lego scene only..."
    cd "$DATA_DIR"
    # Extract only the lego scene from the zip
    unzip -q -o nerf_synthetic.zip "nerf_synthetic/lego/*" -d . 2>/dev/null || \
    unzip -q -o nerf_synthetic.zip "lego/*" -d . 2>/dev/null || \
    unzip -q -o nerf_synthetic.zip -d . 2>/dev/null

    # Normalize directory structure
    if [ -d "nerf_synthetic/lego" ]; then
        mv nerf_synthetic/lego/* "$DEST/" 2>/dev/null || true
        rm -rf nerf_synthetic
    fi

    rm -f nerf_synthetic.zip
    info "Lego scene ready at: $DEST"
}

# ---------------------------------------------------------------------------
# Option 3: Run COLMAP on user-provided images
# ---------------------------------------------------------------------------
run_colmap_on_images() {
    local IMG_DIR="$1"

    if [ ! -d "$IMG_DIR" ]; then
        error "Image directory not found: $IMG_DIR"
        exit 1
    fi

    local NUM_IMAGES
    NUM_IMAGES=$(find "$IMG_DIR" -maxdepth 1 \( -name "*.jpg" -o -name "*.jpeg" -o -name "*.png" \) | wc -l | tr -d ' ')

    if [ "$NUM_IMAGES" -lt 3 ]; then
        error "Need at least 3 images, found $NUM_IMAGES in $IMG_DIR"
        exit 1
    fi
    info "Found $NUM_IMAGES images in $IMG_DIR"

    # Check COLMAP
    if ! command -v colmap &>/dev/null; then
        warn "COLMAP not found. Installing via Homebrew..."
        if command -v brew &>/dev/null; then
            brew install colmap
        else
            error "Homebrew not found. Please install COLMAP manually:"
            error "  brew install colmap"
            error "  or visit: https://colmap.github.io/install.html"
            exit 1
        fi
    fi

    local DEST="$DATA_DIR/custom"
    mkdir -p "$DEST/images"
    mkdir -p "$DEST/sparse"
    mkdir -p "$DEST/database"

    # Copy images
    info "Copying images to $DEST/images/..."
    cp "$IMG_DIR"/*.{jpg,jpeg,png,JPG,JPEG,PNG} "$DEST/images/" 2>/dev/null || true

    local DB_PATH="$DEST/database/database.db"

    # Step 1: Feature extraction
    info "COLMAP Step 1/4: Feature extraction..."
    colmap feature_extractor \
        --database_path "$DB_PATH" \
        --image_path "$DEST/images" \
        --ImageReader.single_camera 1 \
        --ImageReader.camera_model OPENCV \
        --SiftExtraction.use_gpu 0

    # Step 2: Feature matching
    info "COLMAP Step 2/4: Feature matching..."
    colmap exhaustive_matcher \
        --database_path "$DB_PATH" \
        --SiftMatching.use_gpu 0

    # Step 3: Sparse reconstruction
    info "COLMAP Step 3/4: Sparse reconstruction..."
    mkdir -p "$DEST/sparse/0"
    colmap mapper \
        --database_path "$DB_PATH" \
        --image_path "$DEST/images" \
        --output_path "$DEST/sparse"

    # Step 4: Undistort images (optional, improves quality)
    info "COLMAP Step 4/4: Undistorting images..."
    colmap image_undistorter \
        --image_path "$DEST/images" \
        --input_path "$DEST/sparse/0" \
        --output_path "$DEST/undistorted" \
        --output_type COLMAP

    info "COLMAP processing complete!"
    info "Dataset ready at: $DEST"
    info "Use $DEST as the input to OpenSplat."
}

# ---------------------------------------------------------------------------
# Option 4: Create a minimal synthetic test scene (no download needed)
# ---------------------------------------------------------------------------
create_minimal_test_scene() {
    info "Creating minimal synthetic test scene (no download required)..."
    info "This uses Python to generate a simple multi-view dataset."

    local DEST="$DATA_DIR/minimal_test"
    if [ -d "$DEST/images" ] && [ -f "$DEST/sparse/0/cameras.bin" ]; then
        info "Minimal test scene already exists at $DEST. Skipping."
        return 0
    fi

    python3 << 'PYTHON_SCRIPT'
import os
import sys
import json
import struct
import numpy as np

# Check for PIL
try:
    from PIL import Image, ImageDraw
except ImportError:
    print("[WARN] Pillow not installed. Installing...")
    os.system(f"{sys.executable} -m pip install Pillow")
    from PIL import Image, ImageDraw

DEST = os.path.join(os.path.dirname(os.path.abspath("__file__")),
                     os.environ.get("DEST_DIR", "data/minimal_test"))

# Override DEST from environment
import pathlib
script_dir = pathlib.Path(__file__).parent if '__file__' in dir() else pathlib.Path.cwd()

DEST = os.environ.get("DEST_DIR", str(script_dir / "data" / "minimal_test"))
os.makedirs(os.path.join(DEST, "images"), exist_ok=True)
os.makedirs(os.path.join(DEST, "sparse", "0"), exist_ok=True)

W, H = 800, 600
NUM_VIEWS = 20
RADIUS = 4.0

# Generate a simple 3D scene: colored cubes rendered as 2D projections
def project_point(p3d, cam_R, cam_t, fx, fy, cx, cy):
    """Project 3D point to 2D given camera extrinsics and intrinsics."""
    p_cam = cam_R @ p3d + cam_t
    if p_cam[2] <= 0:
        return None
    x = fx * p_cam[0] / p_cam[2] + cx
    y = fy * p_cam[1] / p_cam[2] + cy
    return (int(x), int(y), p_cam[2])

# Simple 3D objects: colored spheres at known positions
objects = [
    {"center": np.array([0.0, 0.0, 0.0]),  "radius": 0.5, "color": (255, 50, 50)},
    {"center": np.array([1.5, 0.0, 0.0]),  "radius": 0.3, "color": (50, 255, 50)},
    {"center": np.array([-1.0, 0.5, 0.5]), "radius": 0.4, "color": (50, 50, 255)},
    {"center": np.array([0.5, -0.5, -0.5]),"radius": 0.35,"color": (255, 255, 50)},
    {"center": np.array([-0.5, 0.0, 1.0]), "radius": 0.25,"color": (255, 50, 255)},
]

fx = fy = 600.0
cx, cy = W / 2.0, H / 2.0

# Camera poses: circle around the origin
cameras_txt_lines = []
images_txt_lines = []

# Write cameras.bin in COLMAP binary format
# Camera model: PINHOLE (id=1), params: fx, fy, cx, cy
cam_id = 1
cam_model = 1  # PINHOLE

# COLMAP binary: cameras.bin
with open(os.path.join(DEST, "sparse", "0", "cameras.bin"), "wb") as f:
    f.write(struct.pack("<Q", 1))  # num_cameras
    f.write(struct.pack("<i", cam_id))
    f.write(struct.pack("<i", cam_model))
    f.write(struct.pack("<Q", W))
    f.write(struct.pack("<Q", H))
    for p in [fx, fy, cx, cy]:
        f.write(struct.pack("<d", p))

# Generate camera poses and render images
image_records = []

for i in range(NUM_VIEWS):
    angle = 2 * np.pi * i / NUM_VIEWS
    elev = 0.3 * np.sin(2 * np.pi * i / NUM_VIEWS * 2)  # slight up/down

    # Camera position (world coords)
    cam_pos = np.array([
        RADIUS * np.cos(angle),
        elev,
        RADIUS * np.sin(angle)
    ])

    # Look-at matrix (looking at origin)
    forward = -cam_pos / np.linalg.norm(cam_pos)
    right = np.cross(forward, np.array([0, 1, 0]))
    right = right / (np.linalg.norm(right) + 1e-8)
    up = np.cross(right, forward)

    # Rotation matrix (world-to-camera)
    R = np.array([right, -up, forward])  # 3x3
    t = -R @ cam_pos  # translation

    # Convert R to quaternion (COLMAP convention: qw, qx, qy, qz)
    def rot_to_quat(R):
        tr = R[0,0] + R[1,1] + R[2,2]
        if tr > 0:
            s = np.sqrt(tr + 1.0) * 2
            w = 0.25 * s
            x = (R[2,1] - R[1,2]) / s
            y = (R[0,2] - R[2,0]) / s
            z = (R[1,0] - R[0,1]) / s
        elif R[0,0] > R[1,1] and R[0,0] > R[2,2]:
            s = np.sqrt(1.0 + R[0,0] - R[1,1] - R[2,2]) * 2
            w = (R[2,1] - R[1,2]) / s
            x = 0.25 * s
            y = (R[0,1] + R[1,0]) / s
            z = (R[0,2] + R[2,0]) / s
        elif R[1,1] > R[2,2]:
            s = np.sqrt(1.0 + R[1,1] - R[0,0] - R[2,2]) * 2
            w = (R[0,2] - R[2,0]) / s
            x = (R[0,1] + R[1,0]) / s
            y = 0.25 * s
            z = (R[1,2] + R[2,1]) / s
        else:
            s = np.sqrt(1.0 + R[2,2] - R[0,0] - R[1,1]) * 2
            w = (R[1,0] - R[0,1]) / s
            x = (R[0,2] + R[2,0]) / s
            y = (R[1,2] + R[2,1]) / s
            z = 0.25 * s
        return np.array([w, x, y, z])

    quat = rot_to_quat(R)

    # Render image
    img = Image.new("RGB", (W, H), (200, 200, 220))
    draw = ImageDraw.Draw(img)

    # Sort objects by depth for painter's algorithm
    depths = []
    for obj in objects:
        p_cam = R @ obj["center"] + t
        depths.append(p_cam[2])
    order = np.argsort(depths)[::-1]  # far to near

    for idx in order:
        obj = objects[idx]
        center_proj = project_point(obj["center"], R, t, fx, fy, cx, cy)
        if center_proj is None:
            continue
        px, py, depth = center_proj
        # Apparent radius in pixels
        r_px = int(fx * obj["radius"] / depth)
        if r_px < 1:
            continue
        # Draw filled circle with shading
        for dr in range(r_px, 0, -1):
            factor = 0.3 + 0.7 * (1.0 - dr / r_px)
            color = tuple(int(c * factor) for c in obj["color"])
            draw.ellipse([px - dr, py - dr, px + dr, py + dr], fill=color)

    img_name = f"image_{i:04d}.png"
    img.save(os.path.join(DEST, "images", img_name))

    image_records.append({
        "image_id": i + 1,
        "quat": quat,
        "t": t,
        "camera_id": cam_id,
        "name": img_name,
    })

# Write images.bin
with open(os.path.join(DEST, "sparse", "0", "images.bin"), "wb") as f:
    f.write(struct.pack("<Q", len(image_records)))
    for rec in image_records:
        f.write(struct.pack("<i", rec["image_id"]))
        for q in rec["quat"]:
            f.write(struct.pack("<d", q))
        for tv in rec["t"]:
            f.write(struct.pack("<d", tv))
        f.write(struct.pack("<i", rec["camera_id"]))
        # Image name as null-terminated string
        name_bytes = rec["name"].encode("utf-8") + b"\x00"
        f.write(name_bytes)
        # Number of 2D points (0 for our synthetic scene)
        f.write(struct.pack("<Q", 0))

# Write empty points3D.bin
with open(os.path.join(DEST, "sparse", "0", "points3D.bin"), "wb") as f:
    f.write(struct.pack("<Q", 0))

print(f"[INFO] Minimal test scene created at: {DEST}")
print(f"[INFO] {NUM_VIEWS} views, {W}x{H} resolution")
print(f"[INFO] Camera intrinsics: fx=fy={fx}, cx={cx}, cy={cy}")
PYTHON_SCRIPT
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
main() {
    info "=== Line 3 Experiment: Data Preparation ==="
    info "Script directory: $SCRIPT_DIR"

    mkdir -p "$DATA_DIR"
    mkdir -p "$RESULTS_DIR"

    if [ "${1:-}" = "--from-images" ]; then
        if [ -z "${2:-}" ]; then
            error "Usage: $0 --from-images /path/to/image/folder"
            exit 1
        fi
        run_colmap_on_images "$2"
    elif [ "${1:-}" = "--minimal" ]; then
        export DEST_DIR="$DATA_DIR/minimal_test"
        create_minimal_test_scene
    elif [ "${1:-}" = "--lego" ]; then
        download_nerf_synthetic_lego
    else
        # Default: try minimal first (no download), then offer alternatives
        info ""
        info "Available data preparation modes:"
        info "  ./prepare_data.sh                     # Create minimal synthetic scene (default, no download)"
        info "  ./prepare_data.sh --minimal            # Same as above"
        info "  ./prepare_data.sh --lego               # Download NeRF synthetic lego (~50 MB)"
        info "  ./prepare_data.sh --from-images DIR    # Run COLMAP on your own images"
        info ""
        info "Defaulting to minimal synthetic scene..."
        export DEST_DIR="$DATA_DIR/minimal_test"
        create_minimal_test_scene
    fi

    # Verify output
    info ""
    info "=== Verification ==="
    local SCENE_DIR
    if [ "${1:-}" = "--from-images" ]; then
        SCENE_DIR="$DATA_DIR/custom"
    elif [ "${1:-}" = "--lego" ]; then
        SCENE_DIR="$DATA_DIR/lego"
    else
        SCENE_DIR="$DATA_DIR/minimal_test"
    fi

    if [ -d "$SCENE_DIR/images" ]; then
        local N_IMGS
        N_IMGS=$(find "$SCENE_DIR/images" -maxdepth 1 \( -name "*.jpg" -o -name "*.jpeg" -o -name "*.png" \) | wc -l | tr -d ' ')
        info "Images found: $N_IMGS"
    fi

    if [ -d "$SCENE_DIR/sparse/0" ]; then
        info "COLMAP sparse model found at: $SCENE_DIR/sparse/0"
        ls -la "$SCENE_DIR/sparse/0/"
    fi

    info ""
    info "=== Data preparation complete ==="
    info "Scene directory: $SCENE_DIR"
    info "Next step: python run_baseline.py --data-dir $SCENE_DIR"
}

main "$@"
