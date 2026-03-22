# CLASS → MLX GPU 完整移植計劃

## CLASS 原始碼分析

CLASS 總共 65,423 行 C 代碼：

| 檔案 | 行數 | 功能 | 移植優先級 |
|------|------|------|-----------|
| **perturbations.c** | **10,506** | Boltzmann ODE + 源函數 | ★★★ 最核心 |
| transfer.c | 5,127 | Bessel 投影 C_l | ★★★ |
| thermodynamics.c | 4,344 | 重組 + 光學深度 | ★★ |
| background.c | 3,730 | Friedmann + 距離 | ★★ |
| primordial.c | 3,458 | 原初功率譜 | ★ |
| lensing.c | 2,360 | CMB lensing | ★★ |
| nonlinear.c | 4,954 | Halofit P(k) | ★ |
| spectra.c | 1,661 | C_l 後處理 | ★ |
| input.c | 3,282 | 參數解析 | ☆ |
| output.c | 1,184 | 檔案輸出 | ☆ |
| tools/*.c | 9,543 | 數值工具（ODE, 陣列, Bessel） | ★★ |
| include/*.h | 8,149 | 頭文件 | — |

## 移植策略：不是逐行翻譯，是理解算法後用 GPU 重寫

### Phase 1：讀懂 perturbations.c（10,506 行）

| 任務 | 內容 | Agent |
|------|------|-------|
| P1.1 | 提取 SOURCE FUNCTION 的精確公式 | Agent A |
| P1.2 | 提取 gauge-invariant 源函數組合 | Agent A |
| P1.3 | 提取 TCA/RSA 切換邏輯 | Agent B |
| P1.4 | 提取 ndf15 在 perturbations 中的調用方式 | Agent B |
| P1.5 | 提取初始條件的完整公式（含高階修正） | Agent C |
| P1.6 | 提取 Silk damping 的精確處理 | Agent C |

### Phase 2：讀懂 transfer.c + thermodynamics.c（9,471 行）

| 任務 | 內容 | Agent |
|------|------|-------|
| P2.1 | 提取 LOS 積分的精確公式和 quadrature 方法 | Agent D |
| P2.2 | 提取 Bessel 函數的精確計算方法 | Agent D |
| P2.3 | 提取 HyRec/RECFAST 的完整重組物理 | Agent E |
| P2.4 | 提取 reionization 的精確模型 | Agent E |

### Phase 3：用 MLX GPU 重新實作

| 任務 | 內容 | 依賴 |
|------|------|------|
| P3.1 | **GPU sync gauge ODE**（IMEX+Strang on MLX） | P1.1-P1.6 |
| P3.2 | **GPU source function**（gauge-invariant，零 cancellation） | P1.1, P1.2 |
| P3.3 | **GPU LOS 積分**（全 MLX tensor） | P2.1, P2.2 |
| P3.4 | **GPU C_l quadrature** | P2.1 |
| P3.5 | 精確 recombination | P2.3, P2.4 |
| P3.6 | 精確初始條件 | P1.5 |

### Phase 4：驗證 + 優化

| 任務 | 內容 |
|------|------|
| P4.1 | 逐模組對比 CLASS 輸出（transfer function, source function, C_l） |
| P4.2 | RMS 收斂到 <0.1% |
| P4.3 | Metal GPU kernel 優化（如需要） |
| P4.4 | Autodiff pipeline（為 MCMC 準備） |

## 執行順序

```
Phase 1 (6 Agents 平行讀 CLASS):
  A: perturbations.c 源函數
  B: perturbations.c ODE/TCA 邏輯
  C: perturbations.c IC + damping
  D: transfer.c LOS + Bessel
  E: thermodynamics.c recombination
  F: tools/*.c ndf15 + 數值方法
    ↓
Phase 2 (整合理解):
  寫出完整的「CLASS 算法規格書」
    ↓
Phase 3 (3 Agents 平行實作 GPU):
  G: GPU ODE + source function
  H: GPU LOS + C_l
  I: recombination + IC
    ↓
Phase 4 (驗證):
  逐步對比到 <0.1%
```

## 預期結果

| 指標 | 目標 |
|------|------|
| TT RMS | <0.1% |
| EE RMS | <0.5% |
| TE RMS | <1% |
| 速度 | <3s（全 GPU） |
| 代碼行數 | ~5,000（Python+MLX，vs CLASS 65,000 C） |
