#!/usr/bin/env python3
"""
Download Sleep-EDF Database Expanded from PhysioNet.

Downloads a subset of subjects (PSG recordings + hypnogram annotations)
for the consciousness-retrodiction experiment.

Dataset: https://physionet.org/content/sleep-edfx/1.0.0/
License: Open Data Commons Attribution License v1.0

Usage:
    python scripts/download_data.py [--n-subjects N]
"""

import os
import sys
import urllib.request
import argparse
from pathlib import Path

# Sleep-EDF Expanded: "sleep-cassette" subset (healthy subjects)
# File naming:
#   PSG:       SC4SSNE0-PSG.edf
#   Hypnogram: SC4SSNXX-Hypnogram.edf  (XX varies per subject: EC, EH, EJ, EP, etc.)
# SS = subject number (00-99), N = night (1 or 2)
BASE_URL = "https://physionet.org/files/sleep-edfx/1.0.0/sleep-cassette"

# (PSG_id, Hypnogram_id) pairs — hypnogram suffixes vary per subject
SUBJECT_FILES = [
    ("SC4001E0-PSG.edf", "SC4001EC-Hypnogram.edf"),   # Subject 00, night 1
    ("SC4002E0-PSG.edf", "SC4002EC-Hypnogram.edf"),   # Subject 00, night 2
    ("SC4011E0-PSG.edf", "SC4011EH-Hypnogram.edf"),   # Subject 01, night 1
    ("SC4012E0-PSG.edf", "SC4012EC-Hypnogram.edf"),   # Subject 01, night 2
    ("SC4021E0-PSG.edf", "SC4021EH-Hypnogram.edf"),   # Subject 02, night 1
    ("SC4022E0-PSG.edf", "SC4022EJ-Hypnogram.edf"),   # Subject 02, night 2
    ("SC4031E0-PSG.edf", "SC4031EC-Hypnogram.edf"),   # Subject 03, night 1
    ("SC4032E0-PSG.edf", "SC4032EP-Hypnogram.edf"),   # Subject 03, night 2
    ("SC4041E0-PSG.edf", "SC4041EC-Hypnogram.edf"),   # Subject 04, night 1
    ("SC4042E0-PSG.edf", "SC4042EC-Hypnogram.edf"),   # Subject 04, night 2
]


def download_file(url: str, dest: Path, timeout: int = 300) -> bool:
    """Download a file with progress reporting."""
    if dest.exists():
        size_mb = dest.stat().st_size / (1024 * 1024)
        print(f"  Already exists: {dest.name} ({size_mb:.1f} MB), skipping.")
        return True

    print(f"  Downloading: {dest.name} ...", end="", flush=True)
    try:
        urllib.request.urlretrieve(url, str(dest))
        size_mb = dest.stat().st_size / (1024 * 1024)
        print(f" done ({size_mb:.1f} MB)")
        return True
    except Exception as e:
        print(f" FAILED: {e}")
        if dest.exists():
            dest.unlink()
        return False


def main():
    parser = argparse.ArgumentParser(description="Download Sleep-EDF data")
    parser.add_argument(
        "--n-subjects", type=int, default=5,
        help="Number of subject files to download (default: 5)"
    )
    parser.add_argument(
        "--data-dir", type=str,
        default=None,
        help="Data directory (default: auto-detect)"
    )
    args = parser.parse_args()

    # Determine data directory
    if args.data_dir:
        data_dir = Path(args.data_dir)
    else:
        script_dir = Path(__file__).resolve().parent
        data_dir = script_dir.parent / "data" / "sleep-edf"

    data_dir.mkdir(parents=True, exist_ok=True)
    print(f"Data directory: {data_dir}")
    print(f"Downloading {args.n_subjects} subject files from Sleep-EDF Expanded...")
    print(f"Source: {BASE_URL}")
    print()

    files_to_download = SUBJECT_FILES[:args.n_subjects]
    success_count = 0
    fail_count = 0

    for psg_name, hyp_name in files_to_download:
        subject_prefix = psg_name[:6]
        print(f"Subject: {subject_prefix}")

        # Download PSG (polysomnography) file
        psg_url = f"{BASE_URL}/{psg_name}"
        psg_dest = data_dir / psg_name
        if download_file(psg_url, psg_dest):
            success_count += 1
        else:
            fail_count += 1

        # Download Hypnogram (sleep stage annotations)
        hyp_url = f"{BASE_URL}/{hyp_name}"
        hyp_dest = data_dir / hyp_name
        if download_file(hyp_url, hyp_dest):
            success_count += 1
        else:
            fail_count += 1

        print()

    print(f"Download complete: {success_count} files succeeded, {fail_count} failed.")
    print(f"Files saved to: {data_dir}")

    # List downloaded files
    edf_files = sorted(data_dir.glob("*.edf"))
    if edf_files:
        print(f"\nDownloaded files ({len(edf_files)}):")
        total_size = 0
        for f in edf_files:
            size_mb = f.stat().st_size / (1024 * 1024)
            total_size += size_mb
            print(f"  {f.name:30s}  {size_mb:8.1f} MB")
        print(f"  {'TOTAL':30s}  {total_size:8.1f} MB")


if __name__ == "__main__":
    main()
