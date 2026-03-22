# mlx_class 開發計劃

## 當前狀態（2026-03-22）

| 指標 | 值 |
|------|-----|
| 模組數 | 60 |
| 總行數 | ~40,000 |
| CMB TT RMS vs CLASS | **8.81%** |
| CMB EE RMS | 20% |
| CMB TE RMS | 31% |
| 第一峰 D_l(220) | 5960 (CLASS: 5740, **1.04×**) |
| 速度 | 420s (scipy sequential) |

---

## Phase 1：精度 + 速度（平行進行）

### Stream A：精度改善（目標 <5% TT RMS）

| 優先 | 任務 | 現狀 | 目標 | 根因 | 預計 |
|------|------|------|------|------|------|
| A1 | 第二峰振幅 | 0.48× | >0.85× | 180 k-modes 不夠 + velocity gauge transform cancellation | 增 N_k→500，multiprocessing 平行 |
| A2 | 低 l 不足 | D_l(10)=294 | ~820 | late ISW tau 點不夠 + reionization 未加 | 加 reionization + 更多 late ISW 點 |
| A3 | 高 l 峰位偏移 | peak3 l=685 | l=815 | oscillation phase 累積 | 增 l_gamma_max 測試 30/35/40 |
| A4 | EE 精度 | 20% RMS | <10% | α_P=1.7 是近似，需完整 polarization hierarchy | 加 Θ_P0, Θ_P2 到 state vector |
| A5 | TE 精度 | 31% RMS | <15% | 依賴 TT 和 EE 的相位正確性 | 跟 A1-A4 一起改善 |

### Stream B：速度優化（目標 <30s）

| 優先 | 任務 | 現狀 | 目標 | 方法 |
|------|------|------|------|------|
| B1 | Python multiprocessing | 420s/180 modes | ~60s | 8 核平行 scipy Radau |
| B2 | 背景快取 | 每次 1.2s | 0.05s | pickle 存/讀 bg 物件 |
| B3 | GPU TCA phase | N/A | N/A | batched ndf15 跑 TCA 到 switch point |
| B4 | 分段策略 | N/A | <30s | TCA(GPU) + full hierarchy(multiprocessing) |

### Stream C：功能補齊（目標：完整 CMB pipeline）

| 優先 | 任務 | 依賴 | 預計 |
|------|------|------|------|
| C1 | Reionization τ_reio | 無 | 在 LOS 加 reion bump: g_reio × exp(-2κ_reio) |
| C2 | BB tensor modes | 無 | 移植 tensor.py 的 h(k,τ) 到 sync gauge |
| C3 | Massive neutrino | C 需改 perturbations_sync | 加 per-q-bin Fermi-Dirac hierarchy |
| C4 | Dark energy w(a) | 需改 background.py | CPL w₀/wₐ 加入 Friedmann |
| C5 | `__init__.py` 更新 | 無 | 匯出 SyncSolver API |
| C6 | 統一 API | C5 | dict → dataclass，與 ProductionSolver 相容 |

---

## Phase 2：工程完善

| 任務 | 說明 |
|------|------|
| D1 | 清理舊 solver（solver_combined, solver_precision 等） |
| D2 | Unit tests for sync gauge |
| D3 | Non-flat Ω_k |
| D4 | Isocurvature IC |
| D5 | GitHub public push |
| D6 | PyPI publish |
| D7 | JOSS 論文 |

---

## 建議執行順序（依賴關係）

```
Week 1: A1 + B1 + C1 + C5 （平行 4 路）
         │       │      │      │
         ▼       ▼      ▼      ▼
      更多k   multiproc reion  API
         │       │
         ▼       ▼
Week 2: A2 + A3 + B2   （精度微調 + 快取）
         │
         ▼
Week 3: A4 + C2 + C3   （polarization + tensor + massive ν）
         │
         ▼
Week 4: C4 + C6 + D1-D7 （dark energy + 工程完善）
```

## 可立即平行啟動的任務

| Agent | 任務 | 改動檔案 | 互不衝突 |
|-------|------|---------|---------|
| Agent 1 | **A1**: N_k→500 + multiprocessing (B1) | solver_sync.py | ✅ |
| Agent 2 | **C1**: Reionization τ_reio | solver_sync.py (LOS section only) | ⚠️ 需協調 |
| Agent 3 | **C5**: __init__.py + 統一 API | __init__.py, 新檔 | ✅ |
| Agent 4 | **A3**: l_gamma_max 測試 | perturbations_sync.py (parameter only) | ✅ |
