# mlx_class 待調整內容

## 精度：6.37% → <1%

### 唯一根因：0.5% photon phase shift

已排除的：
- ❌ Background（14000× 改善，C_l 不變）
- ❌ Gauge transform（machine precision 一致）
- ❌ LOS 源函數公式（simple = direct = efficient 都一樣）
- ❌ k-mode 數量（500→1000 不改善）
- ❌ l_gamma_max（20→40 不改善）

確認的根因：
- ✅ delta_cdm 匹配 CLASS 到 0.78%（所有 k）
- ✅ delta_gamma 振幅正確
- ✅ delta_gamma 有 ~0.5% 相位偏移

### 修法：CLASS 的 TCA slip/shear 高階修正

CLASS perturbations.c 的 `perturb_tca_slip_and_shear()`（~100 行）計算：

```
slip = (theta_gamma - theta_b) 的解析公式
     = tau_c/(1+R) × [-2*calH*R*theta_b/(1+R) + k²*(delta_g/4 - sigma_g)]

shear = sigma_g 的解析公式
      = (16/45)*tau_c*(theta_g + metric_shear)
```

「compromise_CLASS」scheme 包含 leading + selected second-order 項。
這些修正直接影響光子振蕩的相位。

### 實作方式

在 perturbations_sync.py 中：
1. TCA 期間（|κ̇|/k > 30）：用解析 slip 修正 theta_gamma
2. 在 TCA→full switch 時：用解析 shear 設 F_gamma_2 初始值
3. 不改變 ODE solver — 只改 RHS 中的耦合項

之前嘗試失敗的原因：
- 直接 forcing F_gamma_2 = TCA value 破壞了 ODE 連續性
- 正確做法是修改 COUPLING TERM（theta_b 方程中的 slip 修正）

### 預期效果
- 0.5% phase shift → <0.1%
- 6.37% C_l RMS → <1%

## 速度：已優化

| Solver | 時間 | 用途 |
|--------|------|------|
| BDF (current) | 146s | 精度最佳 |
| TCA+RSA | 21s | 日常使用 |
| Green's fn | 1.1ms | MCMC |
| Magnus GPU | 7.6s | GPU demo |

## 論文：待投稿

| 論文 | 狀態 | 下一步 |
|------|------|--------|
| Paper 1 | 等 endorsement | 催 Wilde |
| Paper 2 | 可投 | 投 PRD |
| Paper 3 | 13頁完成 | 最終校對 |
| Paper 9+ext | 6頁+延伸 | 整合 |

## 工程：待完成

- GitHub push（mlx_class 獨立 repo）
- PyPI publish
- JOSS 論文
