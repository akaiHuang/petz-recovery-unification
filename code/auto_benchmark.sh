#!/bin/bash
# Auto-benchmark: polls Tuna-9 every 60s, runs benchmark when idle
# Usage: nohup bash code/auto_benchmark.sh &

PYTHON="/Users/akaihuangm1/.pyenv/versions/3.12.3/bin/python"
DIR="/Users/akaihuangm1/Desktop/github/petz-recovery-unification"
LOG="$DIR/results/auto_benchmark.log"

cd "$DIR"
echo "[$(date)] Auto-benchmark started. Polling every 60s..." | tee "$LOG"

while true; do
    STATUS=$($PYTHON -c "
from qiskit_quantuminspire.qi_provider import QIProvider
for b in QIProvider().backends():
    if b.name == 'Tuna-9': print(b.status.value)
" 2>/dev/null)
    
    echo "[$(date)] Tuna-9: $STATUS" | tee -a "$LOG"
    
    if [ "$STATUS" = "idle" ]; then
        echo "[$(date)] Tuna-9 is IDLE! Starting benchmark..." | tee -a "$LOG"
        $PYTHON code/hardware_benchmark.py --backend "Tuna-9" --shots 4096 --timeout 600 2>&1 | tee -a "$LOG"
        echo "[$(date)] Benchmark finished (exit=$?)" | tee -a "$LOG"
        exit 0
    fi
    
    sleep 60
done
