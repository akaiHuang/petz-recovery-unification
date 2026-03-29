#!/bin/bash
# Auto-analyze benchmark results + run optimizations
# Triggered after benchmark_applications.py completes

PYTHON=/Users/akaihuangm1/.pyenv/versions/3.12.3/bin/python
PROJ=/Users/akaihuangm1/Desktop/github/petz-recovery-unification
LOG=$PROJ/results/benchmarks/benchmark_log.txt
RESULT=$PROJ/results/benchmarks/auto_analysis.txt
export MPLCONFIGDIR=/tmp/mpl

echo "[$(date)] Waiting for benchmark to complete..." >> $RESULT

# Wait for benchmark_log.txt to contain "complete" or "Done"
while true; do
    if [ -f "$LOG" ] && grep -qi "complete\|Done\|SUMMARY\|saved" "$LOG" 2>/dev/null; then
        echo "[$(date)] Benchmark completed! Starting analysis..." >> $RESULT
        break
    fi
    sleep 120
done

# Phase 1: Analyze results
$PYTHON -c "
import json, os, glob
os.environ['MPLCONFIGDIR'] = '/tmp/mpl'

results_dir = '$PROJ/results/benchmarks'
print('='*60)
print('AUTO-ANALYSIS OF BENCHMARK RESULTS')
print('='*60)

# Read all JSON results
for f in sorted(glob.glob(os.path.join(results_dir, '*.json'))):
    print(f'\n--- {os.path.basename(f)} ---')
    with open(f) as fh:
        data = json.load(fh)
    if isinstance(data, dict):
        for k, v in data.items():
            if isinstance(v, (int, float, str, bool)):
                print(f'  {k}: {v}')
            elif isinstance(v, list) and len(v) < 10:
                print(f'  {k}: {v}')
    elif isinstance(data, list) and len(data) < 20:
        for item in data[:5]:
            print(f'  {item}')

print('\n' + '='*60)
print('OPTIMIZATION SUGGESTIONS')
print('='*60)

# Check if H2 VQE shows tau-chrono enables deeper circuits
h2_file = os.path.join(results_dir, 'exp_a_h2_vqe.json')
if os.path.exists(h2_file):
    with open(h2_file) as fh:
        h2 = json.load(fh)
    print('\nH2 VQE:')
    if isinstance(h2, list):
        for r in h2:
            d = r.get('depth', '?')
            print(f'  depth={d}: naive_stop={r.get(\"naive_says_stop\")}, bayes_stop={r.get(\"bayes_says_stop\")}')
    print('  → If Bayes allows deeper → better energy → SUCCESS')

bv_file = os.path.join(results_dir, 'exp_b_bernstein_vazirani.json')  
if os.path.exists(bv_file):
    with open(bv_file) as fh:
        bv = json.load(fh)
    print('\nBernstein-Vazirani:')
    print('  → If success rate at deep circuits > naive prediction → SUCCESS')

print('\nDone.')
" >> $RESULT 2>&1

# Phase 2: Run extended depth benchmark (75, 100, 150, 200) if T9 still available
echo "[$(date)] Attempting extended depth run..." >> $RESULT
$PYTHON -c "
import warnings; warnings.filterwarnings('ignore')
import os, json, numpy as np
os.environ['MPLCONFIGDIR'] = '/tmp/mpl'

try:
    from qiskit_quantuminspire.qi_provider import QIProvider
    p = QIProvider()
    status = None
    for b in p.backends():
        if b.name == 'Tuna-9': status = b.status.name
    
    if status != 'IDLE':
        print(f'Tuna-9 is {status}, skipping extended depth run')
    else:
        print('Tuna-9 IDLE — running extended depths 75, 100, 150, 200...')
        import sys; sys.path.insert(0, '$PROJ')
        import tau_chrono as tc
        from tau_chrono.qi_bridge import get_qi_backend, extract_gate_kraus
        
        backend = get_qi_backend('Tuna-9')
        GATE_POOL = {
            'H': lambda qc, q: qc.h(q), 'X': lambda qc, q: qc.x(q),
            'SX': lambda qc, q: qc.sx(q), 'Rz': lambda qc, q: qc.rz(np.pi/4, q),
            'S': lambda qc, q: qc.s(q),
        }
        gate_kraus = {}
        for g, fn in GATE_POOL.items():
            char = extract_gate_kraus(fn, backend, g, 4096, verbose=False, timeout=600)
            gate_kraus[g] = char.kraus_ops
        
        rho = np.array([[1,0],[0,0]], dtype=complex)
        sigma = np.eye(2, dtype=complex) / 2
        gate_cycle = ['H', 'X', 'SX', 'Rz', 'S']
        
        extended = []
        for depth in [75, 100, 150, 200]:
            gate_list = [gate_cycle[i % 5] for i in range(depth)]
            channels = [gate_kraus[g] for g in gate_list]
            comp = tc.bayesian_compose(channels, sigma, rho, gate_list, min_depth=6)
            row = {'depth': depth, 'tau_naive': comp.tau_multiplicative_total,
                   'tau_bayesian': comp.tau_bayesian_total, 'improvement_pct': comp.improvement_percent}
            extended.append(row)
            print(f'  depth={depth}: naive={row[\"tau_naive\"]:.4f} bayes={row[\"tau_bayesian\"]:.4f} imp={row[\"improvement_pct\"]:+.1f}%')
        
        with open('$PROJ/results/extended_depth_75_200.json', 'w') as f:
            json.dump(extended, f, indent=2)
        print('Saved extended results.')
except Exception as e:
    print(f'Extended depth run failed: {e}')
" >> $RESULT 2>&1

echo "[$(date)] All auto tasks complete." >> $RESULT
