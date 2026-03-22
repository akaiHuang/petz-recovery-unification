# mlx_class 全面優化審計

## 當前基準（2026-03-23）
- Peak-normalized RMS: 7.4%
- Standard relative RMS: 30.6%
- 第一峰: 1.03× CLASS ✅
- 第二峰: 0.59× ❌
- 低 l: 0.40× ❌
- 速度: 285s (full scipy) / 35-60s (hybrid) / 1.1ms (Green's fn cached)

## 審計目標
1. 找出所有可優化的精度瓶頸
2. 找出所有可優化的速度瓶頸
3. 預估每項優化的效果
4. 排優先級
5. 執行

## 審計範圍
- 核心物理: perturbations_sync.py, background.py, recombination.py
- ODE 求解: solver_sync.py, solver_magnus.py, ndf15_batched.py
- LOS 積分: solver_sync.py (LOS section), bessel_cache.py
- C_l 計算: solver_sync.py (C_l section)
- 背景: background.py (Friedmann, tau, R, kappa)
- 所有 GPU 相關: MLX usage patterns

## Agent 分工
- Agent 1: 精度瓶頸（比較每個中間量 vs CLASS）
- Agent 2: 速度瓶頸（profiling + GPU 利用率）
- Agent 3: 數學/算法優化（是否有更好的公式）
