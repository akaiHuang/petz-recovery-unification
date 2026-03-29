# Paper 2: 能做/不能做清單 + 下一步討論方向

**Date**: 2026-03-16
**對應論文**: `papers/paper2_gravity_tau.tex`

---

## 能做（有數學支持）

1. **光速 = 恢復保真度**（弱場精確）
   - c_eff/c = exp(Φ/c²) = exp(−Σ_grav/2) = F_bound
   - 弱場精確（至 2PN order），強場有 O((r_s/r)²) 差異

2. **事件視界 = τ=1 的操作性定義**
   - Σ_grav → ∞ ⟹ F → 0 ⟹ τ = 1
   - 完全資訊損失 = 事件視界

3. **Page curve 用 τ(t) 描述**
   - τ(t) ≤ 1 − exp(−S_Page(t)/2)
   - Page time = τ 從 ~1 轉向 0 的轉折點

4. **HP protocol = τ bound 的實例**
   - Cotler et al. 2019: HP decoder = Petz map（已證明）
   - 這不是新結果，是已知結果在 τ 語言中的重述

5. **糾纏楔重建 = τ bound + RT 面積**
   - F² ≥ exp(−ΔS_gen)（Chen-Penington-Salton 2020）
   - τ_holo = 1 − √(exp(−ΔS_gen)) 有直接幾何意義

6. **Pikovski 退相干 = Petz recovery 問題**
   - N_grav = Tr_internal（結構上等同量子擦除）
   - τ_gravity = 1 − F(ρ₀, R_Petz(N_grav(ρ₀)))
   - 地球上 τ ~ 10⁻¹⁹（基本底線，非巨觀來源）

7. **互補性 = observer-dependent τ**
   - 不同觀察者 = 不同 partial trace = 不同 channel = 不同 τ
   - Infalling: τ_in = 0；External: τ_ext ~ 1

8. **Scrambling rate bound（猜想）**
   - |dτ/dt| ≤ (2π/β)τ
   - 來自 MSS bound，在 τ 語言中的重述
   - 標記為 SPECULATIVE

---

## 不能做（已確認的 dead ends）

1. **暗物質**
   - τ_grav ~ 10⁻¹⁹ 在地球表面
   - 太小，無法產生星系尺度的動力學效應
   - 暗物質需要 Papers 3/4 的 Khronon 機制

2. **暗能量**
   - τ 框架不預測 Λ 的值
   - 只描述資訊損失，不產生斥力

3. **g_ab = f(τ)（反向推導不可能）**
   - g₀₀ → Σ → τ 是單向的
   - 從 τ 不能重建完整度規 g_ab
   - 需要方向性資訊（directional recovery fidelity）

4. **解決 horizon problem**
   - 需要因果連接（causal connection）
   - τ 本身不提供這個

5. **c_eff = c(1−τ) 精確成立**
   - 只是一階近似，valid for Σ ≪ 1
   - 精確關係是 c_eff/c = exp(−Σ/2)

6. **巨觀時間箭頭來自重力退相干**
   - Σ_grav ≪ Σ_thermal
   - 重力退相干是基本底線，不是主要來源
   - 巨觀時間箭頭來自熱環境

---

## 下一步討論方向

### 立即可做（本週）
1. **Pikovski channel 的 Petz recovery 明確計算**
   - 構造 N_grav = Tr_internal 的 Kraus 算子
   - 計算 F 作為 Σ_grav 的函數
   - 判斷是否飽和 Petz bound

2. **數值驗證**：exp(Φ/c²) vs F_Petz 在 amplitude damping 模型
   - 用 qubit 模型做 toy calculation
   - 驗證 saturation 假設的合理性

### 短期（2026 Q2）
3. **Gravity Research Foundation 2026 essay**（截止 3/31）
   - 可從 Paper 2 extract 一篇 1500-word essay
   - 核心賣點：c_eff/c = F_bound 的簡潔表述

4. **強場修正**：超越弱場的 τ-gravity 公式
   - 在指數度規和 Schwarzschild 之間做插值
   - 可測試預測：shadow 2.72 vs 2.60 r_s (5%)

### 中期（2026-2027）
5. **Complexity-weighted τ**
   - ER=EPR 需要區分「高 F + 低 complexity」和「高 F + 高 complexity」
   - 可能解決 firewall 問題
   - 與 computational complexity in holography 有交集

6. **Non-Markovian gravity**
   - Σ < 0 在重力場中的意義
   - 與 dark energy / 宇宙加速膨脹的可能連結
   - 極度推測性，需要具體模型

7. **動態時空**
   - 現有框架限於靜態（需要 Killing vector）
   - 擴展到 FRW / gravitational collapse
   - 可能使用 Kodama vector 或 approximate Killing vectors

---

## 與 Papers 3/4 的關係

Paper 2 討論的是 **τ 與重力的一般連結**（任何重力場）。
Papers 3/4 討論的是 **具體的 Khronon 理論**（BS action with J(Y)+K(Q)）。

兩者的關係：
- Paper 2 的 Σ_grav = -ln(-g₀₀) 是**一般結果**
- Papers 3/4 的 Khronon 提供了**具體的物質內容**（K(Q) → DM, J(Y) → MOND）
- Paper 2 不依賴 Khronon，不需要 J(Y) 或 K(Q)
- Paper 2 的唯一假設是「重力紅移 = 量子通道」
