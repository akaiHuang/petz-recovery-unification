#!/usr/bin/env bash
#
# Download Sleep-EDF Expanded subjects SC4001-SC4191 (subjects 00-19, both nights)
# from PhysioNet. Uses wget with --continue and rate limiting.
# Runs 2 parallel downloads at a time.
#
# Usage: bash scripts/download_batch.sh
#

set -euo pipefail

BASE_URL="https://physionet.org/files/sleep-edfx/1.0.0/sleep-cassette"
DATA_DIR="/Users/akaihuangm1/Desktop/github/petz-recovery-unification/eeg_experiment/data/sleep-edf"

mkdir -p "$DATA_DIR"

# All file pairs: PSG + Hypnogram
# Subject 00 (nights 1-2)
# Subject 01 (nights 1-2)
# ...through Subject 19 (nights 1-2)
# Note: Subject 13 has only night 1 (SC4132 does not exist)

FILES=(
    # Subject 00
    "SC4001E0-PSG.edf"
    "SC4001EC-Hypnogram.edf"
    "SC4002E0-PSG.edf"
    "SC4002EC-Hypnogram.edf"
    # Subject 01
    "SC4011E0-PSG.edf"
    "SC4011EH-Hypnogram.edf"
    "SC4012E0-PSG.edf"
    "SC4012EC-Hypnogram.edf"
    # Subject 02
    "SC4021E0-PSG.edf"
    "SC4021EH-Hypnogram.edf"
    "SC4022E0-PSG.edf"
    "SC4022EJ-Hypnogram.edf"
    # Subject 03
    "SC4031E0-PSG.edf"
    "SC4031EC-Hypnogram.edf"
    "SC4032E0-PSG.edf"
    "SC4032EP-Hypnogram.edf"
    # Subject 04
    "SC4041E0-PSG.edf"
    "SC4041EC-Hypnogram.edf"
    "SC4042E0-PSG.edf"
    "SC4042EC-Hypnogram.edf"
    # Subject 05
    "SC4051E0-PSG.edf"
    "SC4051EC-Hypnogram.edf"
    "SC4052E0-PSG.edf"
    "SC4052EC-Hypnogram.edf"
    # Subject 06
    "SC4061E0-PSG.edf"
    "SC4061EC-Hypnogram.edf"
    "SC4062E0-PSG.edf"
    "SC4062EC-Hypnogram.edf"
    # Subject 07
    "SC4071E0-PSG.edf"
    "SC4071EC-Hypnogram.edf"
    "SC4072E0-PSG.edf"
    "SC4072EH-Hypnogram.edf"
    # Subject 08
    "SC4081E0-PSG.edf"
    "SC4081EC-Hypnogram.edf"
    "SC4082E0-PSG.edf"
    "SC4082EP-Hypnogram.edf"
    # Subject 09
    "SC4091E0-PSG.edf"
    "SC4091EC-Hypnogram.edf"
    "SC4092E0-PSG.edf"
    "SC4092EC-Hypnogram.edf"
    # Subject 10
    "SC4101E0-PSG.edf"
    "SC4101EC-Hypnogram.edf"
    "SC4102E0-PSG.edf"
    "SC4102EC-Hypnogram.edf"
    # Subject 11
    "SC4111E0-PSG.edf"
    "SC4111EC-Hypnogram.edf"
    "SC4112E0-PSG.edf"
    "SC4112EC-Hypnogram.edf"
    # Subject 12
    "SC4121E0-PSG.edf"
    "SC4121EC-Hypnogram.edf"
    "SC4122E0-PSG.edf"
    "SC4122EV-Hypnogram.edf"
    # Subject 13 (night 1 only - SC4132 does not exist)
    "SC4131E0-PSG.edf"
    "SC4131EC-Hypnogram.edf"
    # Subject 14
    "SC4141E0-PSG.edf"
    "SC4141EU-Hypnogram.edf"
    "SC4142E0-PSG.edf"
    "SC4142EU-Hypnogram.edf"
    # Subject 15
    "SC4151E0-PSG.edf"
    "SC4151EC-Hypnogram.edf"
    "SC4152E0-PSG.edf"
    "SC4152EC-Hypnogram.edf"
    # Subject 16
    "SC4161E0-PSG.edf"
    "SC4161EC-Hypnogram.edf"
    "SC4162E0-PSG.edf"
    "SC4162EC-Hypnogram.edf"
    # Subject 17
    "SC4171E0-PSG.edf"
    "SC4171EU-Hypnogram.edf"
    "SC4172E0-PSG.edf"
    "SC4172EC-Hypnogram.edf"
    # Subject 18
    "SC4181E0-PSG.edf"
    "SC4181EC-Hypnogram.edf"
    "SC4182E0-PSG.edf"
    "SC4182EC-Hypnogram.edf"
    # Subject 19
    "SC4191E0-PSG.edf"
    "SC4191EP-Hypnogram.edf"
    "SC4192E0-PSG.edf"
    "SC4192EV-Hypnogram.edf"
)

echo "============================================"
echo "Sleep-EDF Batch Download"
echo "Target: ${#FILES[@]} files (20 subjects, 39 nights)"
echo "Directory: $DATA_DIR"
echo "Rate limit: 500 KB/s per connection"
echo "Parallel downloads: 2"
echo "============================================"
echo ""

# Count how many already exist
EXISTING=0
for f in "${FILES[@]}"; do
    if [ -f "$DATA_DIR/$f" ]; then
        EXISTING=$((EXISTING + 1))
    fi
done
echo "Already downloaded: $EXISTING / ${#FILES[@]} files"
echo ""

# Create a URL list file for xargs
URLFILE=$(mktemp /tmp/sleep-edf-urls.XXXXXX)
for f in "${FILES[@]}"; do
    if [ ! -f "$DATA_DIR/$f" ]; then
        echo "$BASE_URL/$f" >> "$URLFILE"
    fi
done

REMAINING=$(wc -l < "$URLFILE" | tr -d ' ')
if [ "$REMAINING" -eq 0 ]; then
    echo "All files already downloaded!"
    rm "$URLFILE"
else
    echo "Downloading $REMAINING remaining files (2 at a time)..."
    echo ""

    # Download with xargs for parallelism
    # wget flags:
    #   --continue       : resume interrupted downloads
    #   --limit-rate=500k: polite to PhysioNet
    #   -P               : save to data directory
    #   --no-verbose     : reduce output noise
    cat "$URLFILE" | xargs -P 2 -I {} wget --continue --limit-rate=500k -P "$DATA_DIR" --no-verbose {}

    rm "$URLFILE"
fi

echo ""
echo "============================================"
echo "Download phase complete. Running verification..."
echo "============================================"
echo ""

# Verification
PASS=0
FAIL=0
MANIFEST=""

for f in "${FILES[@]}"; do
    FILEPATH="$DATA_DIR/$f"
    if [ ! -f "$FILEPATH" ]; then
        echo "MISSING: $f"
        FAIL=$((FAIL + 1))
        continue
    fi

    SIZE=$(stat -f%z "$FILEPATH" 2>/dev/null || stat -c%s "$FILEPATH" 2>/dev/null)
    SIZE_MB=$(echo "scale=2; $SIZE / 1048576" | bc)
    SIZE_KB=$(echo "scale=2; $SIZE / 1024" | bc)

    if echo "$f" | grep -q "PSG"; then
        # PSG files should be > 20 MB
        if [ "$SIZE" -lt 20971520 ]; then
            echo "WARN: $f is only ${SIZE_MB} MB (expected > 20 MB) - may be incomplete"
            FAIL=$((FAIL + 1))
        else
            PASS=$((PASS + 1))
        fi
        MANIFEST="${MANIFEST}${f}  ${SIZE_MB} MB\n"
    else
        # Hypnogram files should be > 10 KB
        if [ "$SIZE" -lt 10240 ]; then
            # Hypnogram files are typically 3-5 KB, so 10 KB is too strict
            # PhysioNet hypnograms are usually 3-6 KB, let's use 1 KB as minimum
            if [ "$SIZE" -lt 1024 ]; then
                echo "WARN: $f is only ${SIZE_KB} KB (expected > 1 KB) - may be corrupt"
                FAIL=$((FAIL + 1))
            else
                PASS=$((PASS + 1))
            fi
        else
            PASS=$((PASS + 1))
        fi
        MANIFEST="${MANIFEST}${f}  ${SIZE_KB} KB\n"
    fi
done

echo ""
echo "Verification: $PASS passed, $FAIL failed out of ${#FILES[@]} files"
echo ""

# Write manifest
MANIFEST_FILE="$DATA_DIR/MANIFEST.txt"
{
    echo "# Sleep-EDF Expanded - Downloaded Files Manifest"
    echo "# Generated: $(date -u '+%Y-%m-%d %H:%M:%S UTC')"
    echo "# Source: https://physionet.org/content/sleep-edfx/1.0.0/"
    echo "# Subjects: 00-19 (SC4001-SC4192), 20 subjects, 39 nights"
    echo "#"
    echo "# File                              Size"
    echo "# ----                              ----"
    echo ""

    TOTAL_SIZE=0
    for f in "${FILES[@]}"; do
        FILEPATH="$DATA_DIR/$f"
        if [ -f "$FILEPATH" ]; then
            SIZE=$(stat -f%z "$FILEPATH" 2>/dev/null || stat -c%s "$FILEPATH" 2>/dev/null)
            SIZE_MB=$(echo "scale=2; $SIZE / 1048576" | bc)
            TOTAL_SIZE=$(echo "$TOTAL_SIZE + $SIZE" | bc)
            printf "%-40s %8s MB\n" "$f" "$SIZE_MB"
        else
            printf "%-40s %8s\n" "$f" "MISSING"
        fi
    done

    TOTAL_MB=$(echo "scale=2; $TOTAL_SIZE / 1048576" | bc)
    echo ""
    echo "# Total: $TOTAL_MB MB"
    echo "# Files: ${#FILES[@]}"
    echo "# Subjects: 20 (00-19)"
    echo "# Nights: 39 (subject 13 has only night 1)"
} > "$MANIFEST_FILE"

echo "Manifest written to: $MANIFEST_FILE"
echo ""
echo "Done!"
