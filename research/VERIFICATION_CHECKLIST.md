# Verification Checklist

所有結果的驗證狀態追蹤。標記：
- ✅ 已由獨立人類物理學家驗證（已發表/教科書）
- 🔢 數學可立即驗證（初等工具，任何人10分鐘內可查）
- ⚠️ AI 計算，等待人類驗證
- ❌ 已知有問題

---

## Tier 1：框架基礎（Paper 1）

| # | 結果 | 狀態 | 依據 | 驗證方法 |
|---|------|------|------|----------|
| 1.1 | τ=0 ⟺ Σ=0 ⟺ QMC ⟺ QEC 等價鏈 | ✅ | Petz 1988, Fawzi-Renner 2015, JRSWW 2018 | 已發表定理的直接組合 |
| 1.2 | F ≥ exp(−Σ/2) (JRSWW bound) | ✅ | Junge et al. 2018, Ann. Henri Poincaré | 已發表定理 |
| 1.3 | Petz map = quantum Bayes rule | ✅ | Buscemi 2024, Parzygnat 2023 | 獨立驗證 |
| 1.4 | Petz recovery 實驗實現 | ✅ | Pino 2025 (ion trap), Singh 2025 (NMR) | 實驗結果 |

---

## Tier 2：重力通道（Paper 2）

| # | 結果 | 狀態 | 依據 | 驗證方法 |
|---|------|------|------|----------|
| 2.1 | 重力紅移 = 純損通道 (η = −g₀₀) | ✅ 零件 / ⚠️ 組裝 | Ahmadi-Fuentes 2014 + 我們的推導 | 檢查 WKB → η = f(r₁)/f(r₂) 的兩個因子 |
| 2.2 | Tolman 溫度 = 通道環境 | ✅ | Tolman 1930, Birrell-Davies 1982 | 教科書 |
| 2.3 | Σ_grav = −ln(−g₀₀) (三路線收斂) | ✅ 個別 / ⚠️ 合成 | Dorau-Much 2025, Herrera 2020, Ahmadi-Fuentes 2014 | 三篇已發表論文各給出相同結果 |
| 2.4 | Greybody 修正 ≤ (1/4)(r_s/r)² | ⚠️ | 我們的計算 | 用 Regge-Wheeler Born 近似重新推導 |
| 2.5 | Mode-sum Σ 是 intensive (per d.o.f.) | ⚠️ | 我們的計算 | 檢查 AQFT relative entropy density 文獻 |
| 2.6 | Bogoliubov β 精確為零 (靜態) | 🔢 | 時間平移對稱性的直接推論 | 一行證明：靜態 → ω 守恆 → β = 0 |
| 2.7 | 指數度規 g₀₀ = −exp(−r_s/r) | ⚠️ | Petz 飽和 ansatz | 與 Schwarzschild 弱場一致，強場不同——需 GW echo 實驗 |
| 2.8 | GW echo 延遲 ~4r_s/c | ⚠️ | 我們的計算 | 等待 LIGO O5+ 數據 |

---

## Tier 3：暗物質 / 銀河尺度（Paper 3）

| # | 結果 | 狀態 | 依據 | 驗證方法 |
|---|------|------|------|----------|
| 3.1 | η_G = 1 (anomalous dimension) | ✅ | Kumar 2025, Phys. Lett. B 871 | 已發表，QFT 一圈計算 |
| 3.2 | Running G fit SPARC 100 星系 | ✅ | Gubitosi et al. 2024 | 已發表，chi² 分析 |
| 3.3 | Verlinde > MOND at 5.2σ | ✅ | Ghari & Haghi 2026 | 已發表，衛星星系數據 |
| 3.4 | η_G = 1 從 DPI + extensivity | ⚠️ | 我們的推導 | 檢查三叉排除論證的每一步 |
| 3.5 | a₀ = cH₀/(2π) 從 KMS-Crooks | ⚠️ | 我們的推導 | 檢查 de Sitter 熱週期性 → Crooks → a₀ |
| 3.6 | Σ_DM = Σ_cl − Σ_q ≥ 0 (DPI) | 🔢 | DPI 的直接應用 | 一行證明：conditional expectation 是 CPTP |
| 3.7 | Bullet Cluster 定性解釋 | ⚠️ | 糾纏剛性論證 | 需要定量 N-body 模擬 |

---

## Tier 4：μ 和 CMB（Paper 3 + 4）— 最需要驗證

| # | 結果 | 狀態 | 依據 | 驗證方法 |
|---|------|------|------|----------|
| 4.1 | μ₀ = H₀/c 是維度分析唯一解 | 🔢 | 純維度分析 | 列出 [μ] = m⁻¹ 的所有組合，驗證 d≠0 偏差 10¹²² |
| 4.2 | μ₀ = H₀/c → ρ_eff = ρ_crit/3 | 🔢 | 代入公式 | 計算 μ²c²/(8πG) = H₀²/(8πG) = ρ_crit/3 |
| 4.3 | Ω_K = (δ² + 2δ)/3 (正確公式) | 🔢 | Friedmann + K(Q) | 從 Khronon action 重新推導 |
| 4.4 | η_μ = 1 (DPI + extensivity) | ⚠️ **最關鍵** | 與 η_G 相同的論證 | 逐步檢查 extensivity 假設是否適用於 μ |
| 4.5 | η_μ = 1 (dimensional transmutation) | ⚠️ | μ̃ = μ/k 不動點 | 檢查 beta function 結構 |
| 4.6 | η_μ = 1 (modular flow) | ⚠️ | T_mod(k) = k/π | 檢查 Casini-Huerta modular Hamiltonian 正規化 |
| 4.7 | μ_bg(a) = H(a)/c → w̃₀ ~ 4×10⁻¹⁰ | ⚠️ **關鍵** | 背景 running | 從 Khronon EOM 在 FRW 上重新推導 |
| 4.8 | Khronon + c_s=0 能 fit CMB | ✅ | Skordis & Złośnik 2021, PRL | 已發表，完整 C_ℓ 計算 |
| 4.9 | τ 框架的 running μ 能 fit CMB | ⚠️ | 我們的論證 | **需要跑 CLASS/CAMB 做完整數值驗證** |
| 4.10 | T_dS → μ₀ → {a₀, Ω_DM} 統一鏈 | ⚠️ | 我們的綜合 | 檢查每一步的邏輯和數學 |

---

## Tier 5：觀察者相依（Paper 5）

| # | 結果 | 狀態 | 依據 | 驗證方法 |
|---|------|------|------|----------|
| 5.1 | τ_O 單調性 (DPI) | 🔢 | DPI 的直接應用 | 一行證明 |
| 5.2 | τ_O 閾值定理 | 🔢 | Fawzi-Renner 的直接推論 | 短證明 |
| 5.3 | τ 次可加性 | 🔢 | 標準 QI 不等式 | 短證明 |
| 5.4 | Complementary uncertainty: g(d) = (d+1)−√(2(d+1)) | ⚠️ | 2-design + Cauchy-Schwarz | **驗證三步：(i) F=√P_i (ii) ΣP_i=2 (iii) C-S** |
| 5.5 | g(2) = 3−√6 ≈ 0.5505 (tight) | ⚠️ | 上面的 d=2 特例 | 數值驗證：找到達成等號的態 |
| 5.6 | 混態不成立 (maximally mixed → τ=0) | 🔢 | 反例 | 直接計算 I/d 的 τ |

---

## 驗證優先順序

### 🔴 第一批（決定框架成敗）
1. **4.4** η_μ = 1 — 如果錯，CMB 方案崩潰
2. **4.7** w̃₀ ~ 4×10⁻¹⁰ — 如果錯，dust-like 行為不成立
3. **5.4** Complementary uncertainty — 如果錯，Paper 5 核心定理失效

### 🟠 第二批（加強框架）
4. **2.4** Greybody bound
5. **3.4** η_G = 1 的 DPI 推導
6. **4.1-4.3** μ₀ 的維度分析（應該最容易驗證）

### 🟡 第三批（完善框架）
7. **2.1** Channel theorem 完整組裝
8. **3.5** a₀ 從 KMS-Crooks
9. **4.9** 完整 C_ℓ 數值計算（最耗時）

### 🟢 第四批（未來工作）
10. **2.7-2.8** 指數度規 + echo（需要實驗數據）
11. **3.7** Bullet Cluster 定量
12. **4.10** 統一鏈完整性

---

## 驗證記錄

（每次驗證後在此記錄結果）

| 日期 | 項目 | 結果 | 備註 |
|------|------|------|------|
| — | — | — | 等待開始 |
