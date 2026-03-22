# mlx_class GPU 加速計劃

## M1 Max 硬體
- CPU: 10 核（8P + 2E）— 目前用 multiprocessing
- **GPU: 32 核** — 幾乎閒置（只有 Bessel 用到）
- 統一記憶體 — CPU/GPU 共享，零拷貝

## 當前瓶頸分析（按耗時排序）

| 瓶頸 | 現在 | 耗時 | GPU 可行性 | 預期加速 |
|------|------|------|-----------|---------|
| **1. Sync gauge ODE** | scipy Radau per-k (CPU) | 253s | ★★★ 高 | 253s → 2s |
| **2. Lensing flat-sky** | CPU loop | 3-10s | ★★ 中 | 10s → 1s |
| **3. Massive ν** | scipy per-q-bin (CPU) | 5-30s | ★★★ 高 | 30s → 2s |
| **4. MCMC proposals** | serial per-proposal | 2.8h/5000步 | ★★★ 高 | vmap batch |
| 5. Bessel table build | 一次性 | 2-4s | 已完成 | 570× ✅ |
| 6. LOS 積分 | GPU partial | 1.1s | ★ 低 | 已夠快 |
| 7. 背景計算 | numpy | 1.2s | ☆ 無需 | 已夠快 |

## GPU 化開發計劃

### Phase G1：Sync gauge ODE on GPU（最大瓶頸）
**目標：253s → 2s**

核心思路：把 solver_fast.py 的 IMEX+Strang GPU 架構移植到 sync gauge 方程式。

| 步驟 | 做什麼 | 改什麼 |
|------|--------|--------|
| G1.1 | **Sync gauge batched RHS on MLX** | 新檔 `perturbations_sync_gpu.py` |
| G1.2 | **Strang splitting for Thomson** | collision 解析處理（exp(-κ̇ dt)），streaming explicit RK4 |
| G1.3 | **IMEX for η equation** | η' = (3/2)H₀²/k² × src0i（無 stiffness 但需穩定） |
| G1.4 | **h' diagnosed per step** | 00-constraint，batch over k |
| G1.5 | **Gauge transform on GPU** | batch Φ_N, Ψ_N, Θ₀_N |
| G1.6 | **LOS on GPU** | 已有 bessel_gpu，改為全 MLX tensor |

**關鍵：** Strang splitting 消除了 Thomson stiffness（不需要 implicit solver），IMEX 處理 metric stiffness。整個 pipeline 在 GPU 上，無需 scipy。

solver_fast.py 已經做到 0.5s/500 modes（conformal Newtonian gauge）。Sync gauge 方程式結構類似，應該能達到相同速度。

### Phase G2：Lensing on GPU
**目標：10s → 1s**

| 步驟 | 做什麼 |
|------|--------|
| G2.1 | Limber integral vectorize on MLX |
| G2.2 | Flat-sky convolution on GPU (mx.fft) |

### Phase G3：Massive ν on GPU
**目標：30s → 2s**

| 步驟 | 做什麼 |
|------|--------|
| G3.1 | Per-q-bin hierarchy batch on GPU |
| G3.2 | 整合到 G1 的 batched solver |

### Phase G4：MCMC on GPU
**目標：2.8h → 10min**

| 步驟 | 做什麼 |
|------|--------|
| G4.1 | Batch C_l evaluation via vmap |
| G4.2 | Gradient MCMC via MLX autodiff |

## 預期最終性能

| 功能 | 現在 | GPU 後 |
|------|------|--------|
| CMB TT (500 k) | 255s | **~2s** |
| CMB TT+EE+TE | 260s | **~3s** |
| + Lensing | +10s | **+1s** |
| + Massive ν | +30s | **+2s** |
| MCMC 5000 步 | 2.8h | **~15min** |
| **全 pipeline** | **~300s** | **~6s** |

## 執行順序

```
G1 (sync ODE on GPU) ← 最大加速，最高優先
  G1.1 → G1.2 → G1.3 → G1.4 → G1.5 → G1.6
    ↓
G2 (lensing GPU) — 可與 G1.6 平行
    ↓
G3 (massive ν GPU) — 依賴 G1 架構
    ↓
G4 (MCMC GPU) — 依賴 G1 完成
```

**G1 是關鍵路徑。** 完成 G1 後其他都是增量改進。

## 不需要 GPU 的功能（保持 CPU）

| 功能 | 原因 |
|------|------|
| background.py | 10ms，無需加速 |
| recombination.py | 0.1s，ODE 太小不值得 |
| dark_energy.py | 0.05s |
| curvature.py | 0.05s |
| inifile.py | 文字解析 |
| non_gaussianity.py | 0.1s |
| sz_effect.py | template，0.2s |
| calibration.py | 後處理 |
