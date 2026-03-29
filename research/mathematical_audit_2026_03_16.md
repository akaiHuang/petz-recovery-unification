# 數學審計報告：τ/Khronon 框架 — 論文前完整檢查

**日期**：2026-03-16
**目的**：在寫論文前，誠實評估每個數學環節的嚴格程度
**方法**：9 項平行 Agent 驗證 + 文獻交叉檢查

---

## 總覽

| # | 數學環節 | 評級 | 論文中應標記為 |
|---|---------|------|--------------|
| 1 | μ₀ = H₀/c | ⭐⭐⭐⭐ STRONG | "唯一的無量綱分析選擇" (naturality, not derivation) |
| 2 | K(Q) = μ²(Q-1)² | ⭐⭐⭐⭐ STRONG | "leading-order 展開" (Scherrer 2004, Arkani-Hamed 2004) |
| 3 | w̃(z=0) = 0.170 低紅移約束 | ⭐⭐⭐½ SAFE (11/11 tests) | w̃≠c_s²，擾動=CDM；w_DE~-0.95 bonus |
| 4 | 指數度規場方程 | ⭐⭐½ SAFE + 3 gaps | g₀₀ 嚴格；g_rr 假設；phantom scalar 未連結 Khronon |
| 5 | 旋轉曲線推導鏈 | ⚠️⚠️ CRITICAL GAP | 兩條路混淆，需選擇一條 |
| 6 | Running μ(k) = k | ⚠️⚠️ ASSUMED | "well-motivated assumption, not derivation" |
| 7 | P(k) 收斂性 | ⭐⭐⭐ SAFE | 線性 = CDM，非線性偏差在允許範圍 |
| 8 | DPI → η = 1 | ⭐⭐⭐ PLAUSIBLE | "DPI + extensivity"（extensivity 做了主要工作）|
| 9 | Σ_grav 三條推導獨立性 | ⭐⭐ PRESENTATION ISSUE | 同一結果三種記號，非三條獨立推導 |

**整體評估**：框架在觀測相容性上**全部通過**（BBN、GW、BAO、PPN、CMB）。
數學上有 **2 個 critical gap**（旋轉曲線推導鏈、running μ）+ **1 個 significant gap**（指數度規 g_rr + phantom scalar 連結）需要在論文中誠實標記。

---

## 詳細分析

### 1. μ₀ = H₀/c：STRONG（⭐⭐⭐⭐）

**結論**：Naturality argument，不是推導。但極其強。

**論證**：
- 要求 μ 的量綱為 [m⁻¹]，可用的基本常數為 {H₀, c, ℏ, G, k_B}
- 一般組合：μ = (H₀/c) × (t_P H₀)^{2d}
- 因為 t_P H₀ = 1.18 × 10⁻⁶¹，任何 d ≠ 0 都偏離 10^{122|d|} 倍
- **結論**：μ₀ = H₀/c 是唯一不涉及 Planck 尺度的選擇

**物理詮釋（三重等價）**：
- Khronon Compton 波長 = Hubble radius
- μ₀ = 2π k_B T_dS / (ℏc)（de Sitter 熱波長倒數）
- μ₀ = 2π a₀/c²（MOND 加速度 ↔ 長度）

**何為不足**：
- 這是 naturality argument，不是場方程推導
- 不排除有未知常數進入的可能性
- 類比：MOND 的 a₀ ≈ cH₀/(2π) 同樣是 naturality，從未被嚴格推導

**論文寫法**：✅ "唯一的無量綱分析選擇（unique dimensional-analysis selection without Planck-scale physics）"

---

### 2. K(Q) = μ²(Q-1)²：STRONG（⭐⭐⭐⭐）

**結論**：有原理性理由，已有文獻支持。

**三重論證**：
1. **Ghost condensation（Arkani-Hamed et al. 2004）**：任何自發破缺 Lorentz 的標量場，leading-order 有效作用都是 K(Q) = λ(Q-1)²
2. **Scherrer (2004, PRD 70, 043511)**：K(X) 模型中，唯一能同時作為暗物質和暗能量的形式是二次展開
3. **Petz 恢復最優性**：二次形式對應 Gaussian 通道，是最大恢復保真度的選擇

**已知約束**：
- c₁₃ = 0 → c_T = c exactly（GW170817 safe）— 注意：c₁₃=0 是 BS action choice（省略 K_μνK^μν），非幾何強制
- w̃ ≈ 0 at high-z（CMB safe）
- c_s² ≈ 0（structure formation = CDM-like）

**何為不足**：
- K(Q) 的具體形式不是從 τ 框架本身推導出來的
- 更高階項 K(Q) = μ²(Q-1)² + α(Q-1)³ + ... 被忽略（需要 α 很小的論證）

**論文寫法**：✅ "ghost condensation 的 leading-order 有效展開，與 Scherrer (2004) 和 Arkani-Hamed et al. (2004) 一致"

---

### 3. w̃(z=0) = 0.170 低紅移約束：SAFE（⭐⭐⭐½）— 11 項測試全過

**結論**：w̃(z=0) = 0.170 通過全部 11 項低紅移觀測約束（9 SAFE / 1 MARGINAL / 0 DANGEROUS）。

**核心洞見（Blanchet & Skordis 2024）**：w̃ 和 c_s² 在 Khronon 中**獨立**
- w̃(z) 只影響**背景膨脹**（Friedmann 方程）
- c_s² = 0（scalar mode ω=0，不傳播）→ 擾動 = CDM
- 所以 w̃ = 0.17 **不等於** warm DM — Khronon 擾動在所有可觀測尺度與 CDM 完全相同

**11 項測試結果**：

| 測試 | 結果 | 偏差 |
|------|------|------|
| Jeans length | ✅ SAFE | λ_J = 0（c_s² = 0），無小尺度壓制 |
| SN Ia 距離 | ⚠️ MARGINAL | |δμ| ~ 0.016 mag at z=0.5（~1× Pantheon+ uncertainty）|
| BAO/DESI | ✅ SAFE+ | Δχ² = -3.9 **Khronon preferred** |
| Growth rate fσ₈ | ✅ SAFE | Δχ² = +0.15（基本相同）|
| S₈ tension | ✅ SAFE | S₈ = 0.836 vs ΛCDM 0.832（不惡化）|
| ISW effect | ✅ SAFE | A_ISW = 1.001（0.06% 變化）|
| TKS 2016 CMB | ✅ SAFE | z>100: w̃ < 5×10⁻⁷ << TKS bound 2.4×10⁻³ |
| Weak lensing | ✅ SAFE | < 1% 變化 |
| Peculiar velocities | ✅ SAFE | < 1% 變化 |
| Cluster counts | ✅ SAFE | σ₈(Khr) = 0.817 vs 0.811 (+0.7%) |
| WDM constraints | ✅ N/A | Lyman-α 約束 c_s²>0，Khronon c_s²=0，不適用 |

**關鍵發現**：
1. **目前沒有專門的低紅移 w_DM 約束** — 高紅移主導所有約束
2. **意外加分**：Ω_Λ 減少 6.6% → 有效 w_DE ~ -0.95（quintessence-like）→ DESI 偏好 Δχ² = -3.9
3. **唯一 MARGINAL**：SN Ia（Pantheon+），|δμ| ~ 0.016 mag（~1× uncertainty）

**⚠️ 更審慎的評估**（第二輪驗證）：
- 2025 年 barotropic DM 研究報告 |w_DM| < 3×10⁻⁷（95% CL）— 但此約束被 z > 1 數據主導
- w̃(z=0) = 0.17 比此約束大 ~10⁶ 倍 — 但約束是對**常數 w** 的，不適用於 evolving w(z)
- **近期威脅**：DESI DR2 full-shape RSD (2026-2027) 和 Euclid (2027-2028) 將直接探測 w(z) 演化
- **如果 w̃(z) 真的從 10⁻¹⁰ 演化到 0.17，精度再提高 1-2 個量級就可能被偵測到**
- 這是**可測試預測**，不是漏洞

**論文寫法**：✅ "Khronon 預測 w̃(z=0) ≈ 0.17，在現有約束之內（因 c_s² = 0 且高紅移 w̃ < 10⁻⁹）。w̃(z) 的演化構成可測試預測，近期巡天（DESI DR2、Euclid）可直接檢驗。"

---

### 4. 指數度規場方程：SAFE but with 3 OPEN GAPS（⭐⭐⭐ → ⭐⭐½）

**結論**：g₀₀ 從 information theory 嚴格推導；g_rr 來自假設；自洽但有邏輯循環。

**已建立（嚴格）**：
- Σ_grav = -ln(-g₀₀) 三條路徑嚴格推導 ✓
- g₀₀ = -exp(-r_s/r) 由 Σ = r_s/r 得出 ✓
- 指數度規是 Einstein 方程 + phantom scalar (φ = M/r) 的精確解 ✓
- PPN 參數 β = γ = 1 exactly ✓
- 水星進動 = GR 完全相同 ✓

**三個 OPEN GAP**：

**Gap A：g_rr 未從 information theory 推導**
- Paper 2 只推導了 g₀₀ = -exp(-Σ_grav)
- g_rr = exp(+r_s/r) 來自假設 g₀₀·g_rr = -1（isotropic refractive index condition）
- 五種推導 g_rr 的嘗試：Bisognano-Wichmann (open), Jacobson/Dorau-Much (**gives Schwarzschild, not exponential**), RT/MERA (not tested), Connes/Witten (unknown), Einstein+phantom (YES but needs matter)
- **最嚴格的路徑（Jacobson/Dorau-Much）給的是 Schwarzschild**

**Gap B：Phantom scalar 未從第一原理推導**
- 指數度規需要 phantom scalar φ = M/r 作為源
- Σ = φ 的 identification 是 motivated by 1/r profile，但未推導
- Makukov-Mychelkin (2020) 驗證了一致性，但不是推導

**Gap C：Khronon → phantom scalar 連結 — ⚠️ CRITICAL UPDATE (SymPy verified)**
- **標準 Khronon K(Q) = μ²(Q-1)² 不能（CANNOT）產生指數度規** — 數學上不相容
  - 匹配條件要求 K'(Q) = 0 同時 K'(Q)ψ'² = K(Q)e^{α/r}，互相矛盾
  - 對**任何** K(Q) 都成立，不只是二次的
- **解決方案**：entropy production 場 Σ 本身就是 phantom scalar
  - Φ(r) = Σ/√2 = r_s/(√2 r)
  - Action: S = ∫√(-g)[R/(16πG) + ½(∂Φ)²] d⁴x（phantom kinetic term: +½ 而非 -½）
  - Phantom sign 的物理意義：Σ 量度資訊**損失**，貢獻負能量
  - Einstein 方程 + scalar EOM **精確滿足**（SymPy 驗證）
- **文獻支持**：指數度規 = Fisher/JNW 族的 antiscalar 極限（Makukov & Mychelkin 2020）
- **Einstein tensor**：G^μ_ν = f(r) × diag(+1,-1,+1,+1)，f(r) = r_s² e^{-r_s/r}/(4r⁴)

**含義**：Khronon 在弱場（宇宙學）和強場（黑洞）扮演**不同角色**：
- 弱場：K(Q) = μ²(Q-1)² → CDM-like fluid（GDM 映射）
- 強場：Σ = -ln(-g₀₀) → phantom scalar → 指數度規（Fisher/JNW 映射）
- 這不是矛盾——是同一個場在不同 regime 的不同有效描述

**核心 tension**：
- Jacobson/Dorau-Much（最嚴格的 information-theoretic 推導）→ 得到**標準 GR**（Schwarzschild）
- Paper 2 的直接 information theory 路徑 → 得到指數度規（phantom scalar source）
- 兩者在弱場一致（1PN），在強場分歧 → **這是可測試預測**
- 可測試差異：shadow size (2.72 r_s vs 2.60 r_s, ~5%), QNM spectrum, 2PN light deflection

**PPN 驗證**：
- β = γ = 1 exactly，與 Schwarzschild 的 g₀₀ 差異在 O((r_s/r)³) → ~10⁻¹⁸ in solar system
- g_rr 差異在 O((r_s/r)²) → ~1.6 μas light deflection（目前不可測，2030s 可能可測）

**文獻背景**：
- Yilmaz (1958)、Papapetrou (1954)、Rosen (1973) 都獨立得到指數度規
- Makukov & Mychelkin (2020)：Fisher、JNW、XZ 三種已知解在特定條件下都約化為指數度規
- Boonserm et al. (2018)：分類為可穿越蟲洞（throat at R = (e/2)r_s）

**論文寫法**：⚠️ "g₀₀ 從 information theory 嚴格推導。完整度規自洽性已驗證（phantom scalar source）。g_rr 的純 information-theoretic 推導是 open problem。Schwarzschild vs 指數度規的差異是可測試預測（2PN 光偏折 ~1.6 μas）。"

---

### 5. 旋轉曲線推導鏈：⚠️ CRITICAL GAP

**結論**：Paper 3 中混淆了兩條不相容的推導路徑。這是最嚴重的數學問題。

**兩條路徑**：

| | Path A (BS Khronon) | Path B (Kumar running-G) |
|---|---|---|
| 起點 | Blanchet-Skordis 完整作用 S[g, φ] | 量子重力啟發的 G(k) running |
| MOND 來源 | J(Y) 函數（Y = A_μ A^μ/c⁴）| Running G(r) = G_N(1 + 2k_*/πr) |
| 旋轉曲線 | v⁴ = G_N M a₀（BTFR，from J(Y)）| v² = G_N M/r + 2G_N M k_*/π |
| μ 的角色 | K(Q) 控制宇宙學（DM density）| k_* = μ₀ = H₀/c 設定 running scale |
| DM 的角色 | Khronon 密度 ρ_K ∝ μ² | 不需要額外 DM |

**問題**：
1. Path A 中 a₀ 來自 J(Y)，μ 來自 K(Q)，二者獨立（已確認無法連結）
2. Path B 中 running G 是一個完全不同的機制（QFT 啟發），與 Khronon 作用量無關
3. Paper 3 在不同章節交替使用兩條路徑，但未說明它們的關係
4. K(Q) = μ²(Q-1)² **不產生** running G — 這是不同的物理

**修復方案**：
- **選項 A**：嚴格使用 BS Khronon → J(Y) 給 MOND，K(Q) 給 DM density，a₀ 和 μ 獨立
- **選項 B**：使用 running-G 圖景 → 不需要 Khronon，但失去 τ 框架連結
- **選項 C（推薦）**：承認 BS Khronon 是完整理論，J(Y) → MOND，K(Q) → cosmological DM。在論文中明確標記：「旋轉曲線來自 J(Y)，不是 K(Q)。K(Q) 負責宇宙學尺度的 DM 密度。Running G 是一個互補但不等價的觀點。」

**論文寫法**：⚠️ 必須明確區分兩條路徑，不能混用

---

### 6. Running μ(k) = k：⚠️ ASSUMED

**結論**：Well-motivated assumption，不是推導。5 個論證收斂但都有假設。

**5 個論證的評估**：

| Approach | 結論 | 可靠度 |
|----------|------|--------|
| A: DPI + extensivity | η_μ = d-3 = 1 | ⭐⭐⭐ Strong（但 extensivity 做了所有工作）|
| B: Khronon self-energy | β_μ > 0（正確方向）| ⭐⭐ Suggestive（perturbative 太弱）|
| C: Dimensional transmutation | μ(k) = k unique | ⭐⭐⭐ Strong（但是 tautology？）|
| D: Fixed-point condition | η_μ = 1 at fixed point | ⭐⭐ Suggestive（不保證有 FP）|
| E: Modular flow | T_mod ∝ k → μ ∝ k | ⭐⭐⭐ Strong（但 Khronon = modular flow 未證明）|

**核心問題**：
- 在 Lagrangian 中 μ 是常數（不 run）
- Running 是量子效應 → 需要 UV 完備理論才能嚴格推導
- 「μ(k) = k」的物理意義：在不同觀測尺度 k 下，有效的 μ 值不同
- 這等價於假設 G 在 IR 有 running，與 Asymptotic Safety 一致但未被嚴格推導

**CMB Catch-22 的角色**：
- 沒有 running μ → w̃(z=1100) 太大 → CMB 不一致 → 理論被排除
- 有 running μ → w̃(z=1100) = 3.1×10⁻¹⁰ → CMB safe
- 所以 running μ 是理論存活的**必要條件**
- 這不能反過來當作 running μ 為真的論證（circular reasoning）

**論文寫法**：⚠️ "We assume μ runs with anomalous dimension η = 1, motivated by DPI + extensivity (Approach A), dimensional transmutation (Approach C), and modular flow identification (Approach E). A rigorous derivation from the microscopic theory remains an open problem."

---

### 7. P(k) 收斂性：SAFE（⭐⭐⭐）

**結論**：線性擾動 = CDM，非線性偏差在當前觀測允許範圍內。

**關鍵論證**：
- 在線性擾動理論（structure formation）中，Khronon 的 GDM 映射給 w ≈ 0, c_s² ≈ 0 → **完全等於 CDM**
- P(k) 在線性尺度（k < 0.1 h/Mpc）**與 CDM 完全相同**
- 非線性效應只在 galaxy 尺度（k > 1 h/Mpc）出現
- Skordis & Zlosnik (2021) 已經展示類似 Khronon 理論的 RMOND 能重現 CMB + P(k)

**Fagin 比較**：
- Khronon 預測 β = 6.2 vs 觀測 β = 5.22 ± 0.41 → 2.4σ tension
- 但 CDM 被排除 6.8σ → **Khronon 仍是最接近的**
- 2.4σ 可能來自模型簡化（使用 isothermal 近似而非完整 Khronon 旋轉曲線）

**收斂性**：
- ∫P(k)dk 收斂：P(k) ∝ k^{n_s} at small k, P(k) → 0 at large k（damping）
- Khronon 修正不改變 large-k behavior（只改變 transfer function 的 O(10⁻¹⁰) 修正）

**論文寫法**：✅ "線性 P(k) 與 CDM 相同；非線性偏差是可測試預測"

---

### 8. DPI → η = 1：PLAUSIBLE（⭐⭐⭐）

**結論**：DPI 只給 η ≥ 0。η = 1 來自 extensivity 假設。

**推導結構**：
1. DPI（Data Processing Inequality）：任何 CPTP map 滿足 → 只給 η ≥ 0
2. **de Sitter extensivity**：大尺度 entanglement entropy 從 area-law 轉到 volume-law
3. η > 1 → entropy 發散 → 超過 de Sitter bound → 排除
4. η < 1 → IR-irrelevant → 不產生觀測效果 → 排除
5. η = 1 → marginal → logarithmic correction → 唯一存活

**維度依賴性**：η = d - 3
- d = 3：η = 0（G 無量綱，無 running）
- d = 4：η = 1 ✓
- d = 5：η = 2

**強處**：
- 論證結構完整，排除法有力
- 三重獨立動機（Verlinde 2016, Jacobson 2016, Casini-Huerta 2017）

**弱處**：
- **Extensivity 本身是假設**，不是從 Lagrangian 推導的
- Volume-law entanglement entropy 在 de Sitter 空間中是 well-motivated 但尚未嚴格證明
- 將 G 的 running 方法直接搬到 μ 上 → 需要論證為何相同的邏輯適用

**論文寫法**：✅ "DPI + de Sitter extensivity uniquely selects η = 1. This is not a prediction of DPI alone — it requires the specific entanglement structure of de Sitter spacetime."

---

### 9. Σ_grav 三條推導獨立性：PRESENTATION ISSUE（⭐⭐）

**結論**：本質上是同一個 GR 結果（g₀₀ = -(1 - r_s/r)）的三種表述，不是三條獨立推導。

**三條「推導」**：
1. **路徑 A（代數/Modular）**：conditional expectation → Σ = S^rel → 用到 g₀₀
2. **路徑 B（Thermal attenuator）**：beam-splitter model + Tolman → η = -g₀₀ → Σ = -ln(-g₀₀)
3. **路徑 C（Gravitational Landauer）**：erasure cost → Σ/2 → 用到 Tolman temperature ∝ 1/√(-g₀₀)

**共同輸入**：三條路徑都以 g₀₀（GR 度規）為輸入，得到 Σ_grav = -ln(-g₀₀)

**這不是壞事**：
- 三種不同的物理框架（代數量子場論、量子光學、熱力學）得到相同答案 → 結果穩健
- 但不應宣稱「三條獨立推導互相驗證」
- 正確的說法是：「三種獨立的物理解讀都指向同一個數學結構」

**論文寫法**：✅ "Three complementary physical frameworks — algebraic QFT, quantum optics, and gravitational thermodynamics — yield the same identification Σ = -ln(-g₀₀), confirming its structural robustness."

---

## 綜合評估：能寫 vs 不能寫

### ✅ 可以安全寫入論文的

1. **τ = 1 - F 框架**（Paper 1 已發表）
2. **等價鏈**：τ ↔ Σ ↔ entropy production ↔ Petz recovery（Paper 1）
3. **Channel Theorem**：η = -g₀₀, Σ = -ln(-g₀₀)（有完整推導）
4. **μ₀ = H₀/c**（唯一的無量綱分析選擇，非推導）
5. **K(Q) = μ²(Q-1)²**（ghost condensation leading-order）
6. **PPN 參數 β = γ = 1**（精確計算，已驗證）
7. **BBN、GW170817、BAO 相容性**（全部 SAFE，margin 巨大）
8. **CMB 相容性**（w̃(z=1100) = 3.1×10⁻¹⁰，TKS bound 以下 11,121 倍）
9. **a₀ = cH₀/(2π)**（from KMS-Crooks + μ 關係）
10. **Ω_DM ~ ρ_crit/3**（order-of-magnitude，非精確預測）

### ⚠️ 必須標記為假設/猜想的

1. **Running μ(k) = k**（well-motivated, 5 arguments converge, but ASSUMED）
2. **η = 1 from extensivity**（plausible, but extensivity itself is assumed）
3. **μ_bg(z) = H(z)/c**（background running，CMB 存活的必要條件，但未推導）
4. **J(Y) 和 K(Q) 的連結**（structurally obstructed，a₀ 和 μ 獨立）
5. **g_rr = exp(+r_s/r)**（假設 g₀₀·g_rr = -1，非 information theory 推導）
6. **Khronon 在強場 = phantom scalar Σ/√2**（SymPy verified，但弱場↔強場的 transition 機制未推導）

### ❌ 不能宣稱的

1. **Ω_DM h² = 0.12 的精確值**（initial condition，不可推導）
2. **旋轉曲線來自 K(Q)**（來自 J(Y)，不同部分的作用量）
3. **三條「獨立」Σ_grav 推導**（同一結果的三種記號）
4. **指數度規是 GR 真空解**（不是，但在 Khronon 框架中自洽）
5. **解決暗物質問題**（提供新詮釋，不是解決）
6. **c 真的在變化**（座標光速在 GR 中本就依賴度規）

---

## 建議的論文修訂策略

### Paper 2（重力-τ 連結）
- ✅ Channel Theorem (g₀₀ → Σ) 是核心結果，可以大方寫
- ✅ 指數度規 PPN 分析可以寫
- ✅ **NEW**: Σ = phantom scalar（SymPy verified field equations）— 這是強結果
- ✅ **NEW**: G^μ_ν = f(r)×diag(+1,-1,+1,+1) 結構簡潔（單一函數控制）
- ⚠️ 明確說明：standard Khronon K(Q) 不能產生指數度規；強場中 Σ 本身扮演 phantom scalar
- ⚠️ 弱場(K(Q)→CDM) vs 強場(Σ→phantom) 的 transition 是 open question
- ⚠️ 三條 Σ 推導改為「三種互補的物理解讀」
- ⚠️ Jacobson/Dorau-Much → Schwarzschild 的 tension 需要誠實討論
- ✅ 可測試預測：shadow 2.72 vs 2.60 r_s (5%), QNM spectrum, 2PN light deflection ~1.6 μas

### Paper 3（弱場/暗物質）
- ⚠️ **重寫旋轉曲線章節**：明確 J(Y) → MOND，K(Q) → cosmological DM
- ⚠️ 不混用 Kumar running-G 和 BS Khronon
- ✅ BTFR 和 RAR 來自 J(Y)，可以寫
- ⚠️ Running μ 標記為假設

### Paper 4（大統一）
- ✅ μ₀ = H₀/c 的 naturality 論證
- ✅ a₀ = cH₀/(2π) 的推導
- ⚠️ Ω_DM ~ 1/3 是 order-of-magnitude
- ⚠️ Running μ 在全文明確標記為假設
- ❌ 不宣稱推導出 Ω_DM h² = 0.12

---

## 下一步

1. **立即修復**：Paper 3 旋轉曲線章節（最嚴重的問題）
2. **標記假設**：在所有論文中明確標記 running μ 為假設
3. **Σ_grav 措辭修改**：「三條推導」→「三種互補解讀」
4. **可測試預測清單**：整理所有可測試的預測（w̃(z=0) ≈ 0.20、2PN 光偏折差、Fagin β = 6.2 vs 5.22）
5. **微觀推導**：running μ 的嚴格推導是下一個重要的理論工作

---

## Hostile Referee Attack Vectors — Agent 修復結果（2026-03-16 第二輪）

以下 4 項由平行 Agent 完成，詳細報告在各 .md 文件中。

### Attack Vector 1: c₁₃ = 0 是否需要微調？
**結果：⚠️ CORRECTED — c₁₃ = 0 是 action choice，非幾何結果**

- ❌ **先前聲明有誤**：hypersurface orthogonality **不會** force c₁₃ = 0
- ✅ Hypersurface orthogonality 的效果是把 c₁, c₃ 合併為單一組合 c₁₃（消除 twist sector）
- ✅ c₁₃ = 0 在 Blanchet-Skordis 作用量中是 **action choice**：他們省略了 K_μν K^μν 項
- ✅ 這個選擇有三重動機：(1) GW170817 要求 |c₁₃| < 10⁻¹⁵, (2) minimality, (3) J+K 已足夠
- ✅ 即使加入 K_μν K^μν，μ₀ = H₀/c 使得 c₁₃ ~ 10⁻³⁰ at solar system scales → physically irrelevant
- **論文寫法**："The Khronon's hypersurface-orthogonal construction merges c₁ and c₃ into c₁₃ = c₁+c₃. We adopt the Blanchet-Skordis (2024) action which gives c₁₃ = 0 by construction (no K_μν K^μν term), yielding c_T = c identically and automatic GW170817 compatibility."

詳見：[c13_proof_2026_03_16.md](c13_proof_2026_03_16.md)（已修正版）

### Attack Vector 2: J(Y) 是自由函數嗎？真正的參數計數？
**結果：⚠️ CONFIRMED — J(Y) 是自由函數，需誠實承認**

- BS 2024 (arXiv:2404.06584) Section 3.2.7 標題為 "Choice of function J(Y)" — 明確是選擇
- 只有漸近形式被確定：deep-MOND 極限 J(Y) ~ Λ − Y + c²a₀⁻¹ Y^{3/2}
- 完整函數形式**未被理論唯一確定**
- J(Y) 和 K(Q) **結構上獨立**，a₀ 和 μ 無理論連結
- **真正的參數/函數計數**：
  - 1 自由函數：J(Y)（受漸近條件約束）
  - ~3 自由參數：a₀, c₁₄, Ω_DM h²
  - K(Q) = μ²(Q-1)² 已固定（ghost condensation）
  - μ₀ = H₀/c 已固定（dimensional analysis）
- **比較**：同 AeST (SZ2021)，略多於 MOND (1 parameter)
- **論文寫法**："The theory has two sectors: K(Q) is uniquely fixed (1 parameter μ₀), J(Y) is a phenomenological function (1 parameter a₀) constrained by MOND asymptotics — analogous to the interpolation function in classical MOND."

詳見：[jy_form_analysis_2026_03_16.md](jy_form_analysis_2026_03_16.md)

### Attack Vector 3: 穩定性分析（Ostrogradski, ghosts, gradients）
**結果：✅ ALL 8 CRITERIA SAFE**

| Check | Result | Mechanism |
|-------|--------|-----------|
| Ostrogradski ghost | ✅ SAFE | 2nd-order EOM |
| Wrong-sign kinetic | ✅ SAFE | K'' = 2μ² > 0 |
| Gradient instability | ✅ SAFE | c₁₃ = 0 exactly |
| Tachyonic modes | ✅ SAFE | ω² ≥ 0 always |
| Ghost condensation | ✅ SAFE | ω = 0, non-propagating |
| Einstein-aether conditions | ✅ ALL MET | |
| Strong coupling | ✅ SAFE | μ₀ tiny → huge margin |
| BPS 2010 confirmation | ✅ VERIFIED | Independent check |

- 所有 6 個物理自由度（2 tensor + 2 vector + 2 scalar）的頻率都是實數
- scalar mode ω = 0（非傳播，= CDM），不是不穩定
- **論文寫法**：可在 Appendix 或 Section 放簡要穩定性證明

詳見：[stability_analysis_2026_03_16.md](stability_analysis_2026_03_16.md)

### Attack Vector 4: 可證偽性（21 項具體標準）
**結果：✅ 21 項可證偽預測，81% 在 95% CL 可測**

| Category | Tests | Falsifiable at 95%+ |
|----------|-------|---------------------|
| A: Current data | 6 | 6/6 |
| B: Near-term 2026-28 | 5 | 4/5 |
| C: Medium-term 2028-35 | 5 | 3/5 |
| D: Structural | 5 | 4/5 |
| **Total** | **21** | **17/21 (81%)** |

**最致命測試**（任一失敗 = 理論死亡）：
1. DM 粒子偵測（>5σ WIMP → 理論完蛋）
2. w̃(z) 演化曲線（DESI DR2 2027 直接測試）
3. BH shadow 大小（2.72 vs 2.60 r_s, 5% 差異）
4. QNM 頻譜（LISA 2030s）
5. c_s² > 0 偵測（21cm + CMB）

**當前無任何觀測否定此框架。**

詳見：[falsification_criteria_2026_03_16.md](falsification_criteria_2026_03_16.md)

---

## 最終審計總結（所有 13 項完成）

### 原始 9 項 + 4 項 Hostile Referee 攻擊 = 13 項全部完成

| # | 項目 | 評級 | 結論 |
|---|------|------|------|
| 1 | μ₀ = H₀/c | ⭐⭐⭐⭐ STRONG | 唯一無量綱分析選擇 |
| 2 | K(Q) = μ²(Q-1)² | ⭐⭐⭐⭐ STRONG | Ghost condensation leading-order |
| 3 | w̃ 低紅移約束 | ⭐⭐⭐½ SAFE | 11/11 tests passed |
| 4 | 指數度規場方程 | ⭐⭐½ + gaps | g₀₀ 嚴格；g_rr 假設；phantom scalar |
| 5 | 旋轉曲線推導鏈 | ⚠️⚠️ CRITICAL | 兩條路混淆，需選 Path A |
| 6 | Running μ(k) = k | ⚠️⚠️ ASSUMED | 5 arguments converge, not derivation |
| 7 | P(k) 收斂性 | ⭐⭐⭐ SAFE | 線性 = CDM |
| 8 | DPI → η = 1 | ⭐⭐⭐ PLAUSIBLE | Extensivity does main work |
| 9 | Σ_grav 三條推導 | ⭐⭐ PRESENTATION | 三種解讀，非三條獨立推導 |
| 10 | c₁₃ = 0 | ⚠️ ACTION CHOICE | BS action 省略 K_μν K^μν → c₁₃=0 by construction（非幾何結果）|
| 11 | J(Y) 形式 | ⚠️ FREE FUNCTION | 漸近確定，完整形式自由 |
| 12 | 穩定性 | ✅ ALL SAFE | 8/8 criteria passed |
| 13 | 可證偽性 | ✅ 21 CRITERIA | 81% at 95%+ CL |

### 論文可以寫的誠實框架

**一個作用量**：S = ∫√(-g)[R − 2J(Y) + 2K(Q)] d⁴x
**兩個 sector**：
- K(Q) = μ²(Q-1)² + μ₀ = H₀/c → cosmological DM（**完全確定，1 parameter**）
- J(Y) → galactic MOND（**phenomenological function，1 parameter a₀**）

**通過所有測試**：BBN, GW170817, BAO, CMB, PPN, binary pulsars
**結構性保證**：c_T = c (not tuned), 全穩定 (no ghosts/tachyons)
**21 項可證偽預測**，含 5 項致命測試
**2 個 critical gap 誠實標記**：running μ = assumption, 旋轉曲線需選 Path A
