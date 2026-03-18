#!/bin/bash
# Auto-benchmark v3: polls Tuna-9 every 60s, runs benchmark when idle
# Includes: PLS CPTP projection, optimization_level=0, barriers
# Usage: nohup bash code/auto_benchmark_v3.sh &

PYTHON="/Users/akaihuangm1/.pyenv/versions/3.12.3/bin/python"
DIR="/Users/akaihuangm1/Desktop/github/petz-recovery-unification"
LOG="$DIR/results/auto_benchmark_v3.log"

export MPLCONFIGDIR="/tmp/mpl"

cd "$DIR"
echo "[$(date)] Auto-benchmark v3 started. Polling every 60s..." | tee "$LOG"

while true; do
    STATUS=$($PYTHON -c "
import warnings; warnings.filterwarnings('ignore')
from qiskit_quantuminspire.qi_provider import QIProvider
for b in QIProvider().backends():
    if b.name == 'Tuna-9': print(b.status.value)
" 2>/dev/null)

    echo "[$(date)] Tuna-9: $STATUS" | tee -a "$LOG"

    if [ "$STATUS" = "idle" ]; then
        echo "[$(date)] Tuna-9 is IDLE! Starting benchmark v3..." | tee -a "$LOG"
        $PYTHON code/hardware_benchmark.py --backend "Tuna-9" --shots 4096 --timeout 600 2>&1 | tee -a "$LOG"
        EXIT_CODE=$?
        echo "[$(date)] Benchmark v3 finished (exit=$EXIT_CODE)" | tee -a "$LOG"
        exit $EXIT_CODE
    fi

    sleep 60
done
