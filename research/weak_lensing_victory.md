# Weak Lensing Evidence: τ Framework vs NFW (2024-03-12)

## 戰果記錄

### Mistele et al. 2024 (ApJ Letters) — 關鍵勝利

**論文**: "Indefinitely Flat Circular Velocities and the Baryonic Tully–Fisher Relation from Weak Lensing"
**arXiv**: 2406.09685
**DOI**: 10.3847/2041-8213/ad54b0

**觀測結果**: 用 DESI 弱透鏡數據反推孤立星系的重力場，發現圓周速度在數百 kpc 內保持平坦，直到 ~1 Mpc 都沒有下降跡象——遠超 NFW 暈的預期截止半徑。

### 計分板

| 預測 | τ 框架 | NFW (ΛCDM) | 觀測 | 勝者 |
|------|--------|-----------|------|------|
| v(r) 在 r > r₂₀₀ | 保持平坦 (ρ∝1/r²) | 下降 ∝ 1/√r | **平坦** | **τ** |
| 到 ~1 Mpc | 不截斷 | 截斷 | **不截斷** | **τ** |
| BTFR 延伸到弱透鏡 | 自然預測 | 需要微調 | 符合 BTFR | **τ** |
| RAR 延伸 2 decades | 自然預測 (Brouwer+2021) | 需要調參 | 符合 | **τ** |

### 為什麼這是勝利

τ 框架 (Paper III):
- G(r) = G_N [1 + 2k_*r/π], k_* ≈ 0.027 kpc⁻¹
- ρ_DM ∝ 1/r², M_DM(r) ∝ r, v(r) = const
- **一個參數 k_*，自然預測，無需調整**

NFW (ΛCDM):
- ρ ∝ 1/r³ at large r → v(r) 應該下降
- 需要加 two-halo term 來事後解釋
- **每個星系需要 2 個參數 (M₂₀₀, c)**

### ΛCDM 的反駁（誠實記錄）

1. **Two-halo term**: 鄰近暈的疊加可能使 v(r) 看似平坦
2. **Comment paper** (arXiv:2407.18263): 質疑反投影方法的系統誤差
3. **Stacking bias**: 疊加上千星系可能掩蓋個別行為

### 額外支持證據

**Brouwer et al. 2021 (KiDS-1000)**:
- 用弱透鏡把 RAR 延伸 2 個數量級到低加速度端
- 與 Verlinde emergent gravity 和 MOND 吻合
- 與某些 ΛCDM 模擬有 >6σ 差異
- arXiv: 2106.11677

**NGC 6505 Einstein Ring (Euclid, O'Riordan+2025)**:
- R_Ein = 2.1 kpc, f_DM = 11.1% (+5.4/-3.5)%
- τ 框架預測: f_DM ≈ 3.5% (running G) + IMF 調整 → 相容
- 不是決定性，但不矛盾
- arXiv: 2502.06505

### 尚未贏的戰場

| 觀測 | τ 框架狀態 | ΛCDM 狀態 |
|------|-----------|-----------|
| CMB acoustic peaks | **未驗證** ⚠️ | 完美擬合 |
| BAO | 未計算 | 完美擬合 |
| 大尺度結構 | 未計算 | 完美擬合 |
| Bullet Cluster | 定性（糾纏剛性）| 定量擬合 |

### 下一步：CMB 是最後戰場

需要用 CLASS/CAMB 跑完整 C_ℓ：
- 輸入: running μ(a) = H(a)/c, η_μ = 1
- 預期: 如果 Skordis & Złośnik 2021 (AeST) 能 fit CMB，
  我們的 running Khronon 也應該能
- 關鍵項目: VERIFICATION_CHECKLIST.md #4.9
