#!/bin/bash
# Poll Tuna-9 every 60s, run benchmarks when IDLE
PYTHON=/Users/akaihuangm1/.pyenv/versions/3.12.3/bin/python
export MPLCONFIGDIR=/tmp/mpl

echo "[$(date)] Waiting for Tuna-9 IDLE..."
while true; do
    STATUS=$($PYTHON -c "
import warnings; warnings.filterwarnings('ignore')
import os; os.environ['MPLCONFIGDIR']='/tmp/mpl'
from qiskit_quantuminspire.qi_provider import QIProvider
p=QIProvider()
for b in p.backends():
    if b.name=='Tuna-9': print(b.status.name)
" 2>/dev/null)
    
    if [ "$STATUS" = "IDLE" ]; then
        echo "[$(date)] Tuna-9 is IDLE! Starting benchmarks..."
        $PYTHON code/benchmark_applications.py --backend "Tuna-9" 2>&1 | tee results/benchmarks/benchmark_log.txt
        echo "[$(date)] Benchmarks complete!"
        break
    else
        echo "[$(date)] Tuna-9: $STATUS — waiting 60s..."
        sleep 60
    fi
done
