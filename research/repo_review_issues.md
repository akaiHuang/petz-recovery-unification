# GitHub Repos 問題清單與修正建議

**日期**：2026-03-24
**範圍**：exponential-metric (Paper 2) + khronon-dark-matter (Paper 3)

---

## A. exponential-metric（Paper 2）— 共 14 個問題

### CRITICAL（必須修，否則會被立即駁回）

#### A1. 飽和假說無物理原理
- **位置**：paper2_gravity_tau.tex:58-60, paper2_supplement.tex:624-684
- **問題**：論文的核心宣稱 `c_eff/c = exp(-Σ/2) = F_bound` 需要 Petz bound 被 saturate，但：
  - supplement 自己證明了 Theorem (line 631)：**Σ > 0 時精確飽和不可能**
  - 然後 line 678-683 退一步說「lower bound geometry」
  - 但 abstract (line 58) 和 main text 仍用 "hypothesize saturation" 的語言
  - **自相矛盾**：你自己的定理證明飽和不可能，但 abstract 假設飽和成立
- **研究者會說**：「你自己的 supplement 就否定了你的 abstract。到底是 saturate 還是 lower bound？」
- **修正建議**：
  1. Abstract 改為：「We identify Σ_grav = -ln(-g₀₀) and show the exponential metric g₀₀ = -exp(-r_s/r) corresponds to the **lower bound geometry** — the metric where recovery fidelity equals the JRSWW bound」
  2. 不再宣稱 "saturation"，改為 "bound geometry"
  3. 把 supplement 的 gap calculation (line 687-733) 提到 main text，展示 Schwarzschild 比 bound 好 O(r_s/r)
  4. 物理論述改為：「Schwarzschild = actual metric, exponential = worst-case bound → 兩者差異 = recovery gap → 可觀測」

#### A2. 指數度規不是推導出來的
- **位置**：paper2_gravity_tau.tex:348-356
- **問題**：論文宣稱指數度規從飽和條件得到，但實際上只是**斷言** g₀₀ = -exp(-r_s/r)，沒有：
  - Variational principle（什麼 action 的解？）
  - Uniqueness 證明（為什麼不是其他 metric？）
  - 動力學機制（什麼物理過程 select 這個 metric？）
- **研究者會說**：「任何人都可以寫下一個 metric 然後說它有某種資訊詮釋。你需要告訴我它從什麼方程式出來。」
- **修正建議**：
  1. 明確改為：「The exponential metric is the **unique static spherically symmetric metric** for which Σ_grav = r_s/r exactly (not approximately)」
  2. 加一個 proposition：證明在 Σ = -ln(-g₀₀) 的定義下，Σ = r_s/r（線性）的唯一解是 g₀₀ = -exp(-r_s/r)
  3. 承認：「We do not derive this metric from an action principle; rather, we identify it as the metric that **exactly linearizes** the QRE-gravity correspondence」
  4. 把這當成一個 **characterization**（特徵化），不是 derivation

#### A3. 三條路線不獨立
- **位置**：paper2_supplement.tex Route A/B/C (lines ~103-334)
- **問題**：宣稱三個獨立路線 converge 到同一結果，但：
  - Route A（Dorau-Much modular flow）：**輸入** g₀₀ → **輸出** Σ
  - Route B（Landauer/Tolman）：**輸入** g₀₀ → **輸出** Σ
  - Route C（Bosonic channel）：**輸入** g₀₀ → **輸出** Σ
  - 全部都是 g₀₀ → Σ 的方向，不是獨立推導
- **研究者會說**：「這是同一個公式的三種寫法，不是三個獨立驗證。」
- **修正建議**：
  1. 改稱 "three consistent physical interpretations"，不稱 "three independent routes"
  2. 明確寫：「All three frameworks take the metric g₀₀ as input and produce the same information-theoretic quantity Σ. This is a **consistency check**, not independent derivation.」
  3. 加一段討論反方向（Σ → g₀₀）需要什麼額外假設

### HIGH（嚴重問題，會大幅降低可信度）

#### A4. Pikovski 章節自相矛盾
- **位置**：paper2_gravity_tau.tex:362-469, paper2_supplement.tex:1003-1035
- **問題**：
  - Main text (line 375-376)：「has the exact structure of a Petz recovery problem」
  - Supplement (line 1003)：明確標題 "Why Pikovski is NOT the gravitational channel"
  - 列出三個 mismatch（probe dependence, divergent supremum, basis）
- **研究者會說**：「你的 main text 說它是 Petz recovery，你的 supplement 說它不是。到底是哪個？」
- **修正建議**：
  1. Main text 改為：「The Pikovski mechanism has the **formal structure** of a Petz recovery problem, but is fundamentally distinct from the gravitational channel Σ_grav」
  2. 明確分兩層：
     - **形式類比**：tracing out internal DOF → τ > 0 ✓
     - **物質差異**：Σ_Pik 依賴 probe，Σ_grav 普適 ✗
  3. 把 supplement 的三個 mismatch 精簡版搬到 main text

#### A5. 黑洞章節過度宣稱
- **位置**：paper2_gravity_tau.tex Table 3 (lines ~489-611)
- **問題**：把已知結果（HP protocol, entanglement wedge）標記為 "NEW SYNTHESIS"，但實際只是把 τ 的符號貼上去
- **修正建議**：
  1. Table 標記改為三類：KNOWN（文獻已有）/ REINTERPRETATION（已知結果的 τ 語言）/ CONJECTURE（推測）
  2. HP protocol 和 EW reconstruction：改標為 KNOWN（Cotler et al. 2019 已經用 Petz map）
  3. 只有 "complementarity as observer-dependent τ" 和 "scrambling rate bound" 是真正的 NEW

#### A6. GW 預測缺誤差分析
- **位置**：gw_predictions_section.tex Tables 1-3
- **問題**：EMRI ΔΦ = 67,000 rad, ISCO shift 5.7%, shadow 4.7% — 全部沒有 error budget
- **修正建議**：
  1. 每個預測加入：source uncertainty (mass, spin) + detector sensitivity + systematic
  2. 特別是 EMRI：質量估計有 ~10% 不確定性，累積相位偏差會跟著放大
  3. 加一個 "detectability" 欄位（S/N ratio 或 sigma level）

#### A7. 缺 compiled PDF
- **位置**：repo root
- **問題**：只有 .tex 檔，沒有 PDF。閱讀者需要自己編譯
- **修正建議**：用 pdflatex 編譯所有 .tex，加入 repo（或加 GitHub Actions 自動編譯）

### MEDIUM（需要修但不致命）

#### A8. [KNOWN] / [NEW] 標記格式不專業
- **位置**：paper2_gravity_tau.tex abstract (lines 50-68)
- **問題**：方括號標記 [KNOWN], [NEW SYNTHESIS], [NEW] 看起來像草稿筆記
- **修正建議**：移除方括號標記。改在 Section 開頭用完整句子說明哪些是新的

#### A9. 自引 Paper 2 沒有 DOI
- **位置**：paper2_gravity_tau.tex bibliography
- **問題**：引用 "PaperII_exp" 但只寫 "(2026)"，沒有 arXiv ID 或 DOI
- **修正建議**：上傳 Zenodo 取得 DOI，或標明 "in preparation"

#### A10. 缺少關鍵引用
- **問題**：缺少 Caticha 2025 (arXiv:2511.19238), Blanchet-Skordis 2024 (arXiv:2404.06584) 在 main paper bibliography
- **修正建議**：加入這些引用，特別是 BS2024 因為 Khronon 是整個系列的基礎

#### A11. Wormhole 詮釋缺物理後果
- **位置**：paper2_gravity_tau.tex:351
- **問題**：提到「classified as a traversable wormhole」但沒討論：穩定性？NEC violation？可觀測後果？
- **修正建議**：加 1-2 段討論，或明確說 "the wormhole interpretation is geometrical; stability analysis is deferred to future work"

#### A12. Table 1 數值未標明近似
- **位置**：paper2_gravity_tau.tex Table 1
- **問題**：中子星 Σ = 0.33 仍用弱場公式 Σ ≈ r_s/r，但 0.33 不算 ≪ 1
- **修正建議**：Table caption 加：「Σ_grav computed using weak-field approximation; corrections are O(Σ²) ≈ 10% for neutron stars」

#### A13. notation 不一致
- **問題**：supplement 和 main text 對 η（transmissivity）的定義不同處需要統一
- **修正建議**：統一為 η = -g₀₀（intensity transmissivity），在第一次出現時明確定義

#### A14. 「能做/不能做」列表位置不對
- **位置**：paper2_gravity_tau.tex Section V (Discussion 末尾)
- **問題**：這個重要的 scope statement 藏在最後，讀者可能已經被過度宣稱的語言誤導
- **修正建議**：移到 Section II 結尾（identification 之後、應用之前），讓讀者一開始就知道邊界

---

## B. khronon-dark-matter（Paper 3）— 共 12 個問題

### CRITICAL（必須修）

#### B1. No-go 定理的覆蓋範圍
- **位置**：paper3_cs2_time_dependent.tex:328-398
- **問題**：Theorem 宣稱「there exists no function μ(a) such that...」但證明只覆蓋三種 case：
  - Case 1：δ(a_rec) ≪ 1 → 17% deviation from a⁻³
  - Case 2：δ(a_rec) ≫ 1 → c_s² = O(1)
  - Case 3：δ(a_rec) ~ O(1) → c_s² ~ 0.1-0.3
- **問題在哪**：Case 1 的結論依賴 δ 從 O(1)（今天）到 ≪ 1（recombination）的過渡造成 17% deviation。但如果 μ(a) 被特別設計來補償這個過渡呢？定理沒有排除 fine-tuned μ(a)
- **研究者會說（BS）**：「你的 Case 1 假設 δ₀ = 0.34，但如果 μ(a) 讓 δ 永遠是 O(1) 但 c_s² 仍然小呢？你只覆蓋了三個 regime，沒覆蓋 transition region。」
- **修正建議**：
  1. 選項 A（加強定理）：加入連續性論證，證明 δ(a) 在三個 case 之間的過渡不可能同時滿足所有條件。關鍵是 c_s² = δ/(2+δ) 是 δ 的單調函數，所以 c_s² < 10⁻⁶ **強制** δ < 2×10⁻⁶
  2. 選項 B（誠實標記）：改 Theorem statement 為「For the class of analytic K(Q) with ghost condensation, any μ(a) satisfying conditions (1)-(3) requires a fine-tuned transition that violates condition (3) by ≥ 17%」— 這更精確

#### B2. BS2025 arXiv 號碼可能不存在
- **位置**：paper3_weak_field.tex:1338-1341, paper3_alpha_derivation.tex:545, paper3_cs2_time_dependent.tex:536
- **問題**：arXiv:2507.00912 — 如果 BS2025 發表於 2025 年 7 月，而現在是 2026 年 3 月，這應該存在。但需要驗證
- **研究者會說**：如果號碼錯誤，這是嚴重的不專業
- **修正建議**：
  1. 立即驗證 arXiv:2507.00912 是否存在
  2. 如果不存在：查正確 ID，或改為 "L. Blanchet, E. Polito, and C. Skordis, *Khronon-Tensor theory* (in preparation, 2025)"
  3. 同時更新所有三個文件

#### B3. Crooks → μ(x) 的推導缺失
- **位置**：paper3_weak_field.tex:669-694
- **問題**：宣稱 Crooks fluctuation theorem 給出 μ(x) = 1 - e^{-√x}，但：
  - Line 674-675 只寫了 Crooks relation P(Σ)/P(-Σ) = e^Σ
  - Line 678 直接寫出 μ(x) = 1 - e^{-√x}
  - **中間的推導步驟完全缺失**：怎麼從 P(Σ)/P(-Σ) = e^Σ 得到 μ(x) = 1 - e^{-√x}？
  - 需要：Σ = √(g_bar/a₀) 的 identification + P(Σ) 的物理詮釋
- **研究者會說（BS）**：「Crooks theorem 是非平衡統計力學的結果，跟星系旋轉曲線有什麼關係？你需要解釋 Σ 在星系尺度的物理意義，以及為什麼 P(Σ) 和 MOND interpolation 有關。」
- **修正建議**：
  1. 選項 A（補推導）：加入完整推導鏈：
     - 定義 Σ_gal = 2 ln(g_obs/g_bar) 在 MOND regime
     - 用 Crooks: P(forward)/P(backward) = e^Σ
     - 定義 μ = P(forward) = 1/(1 + e^{-Σ})... 但這給 logistic，不是 1-e^{-√x}
     - **注意**：我不確定 1 - e^{-√x} 真的能從 Crooks 推出來。需要仔細檢查
  2. 選項 B（降級宣稱）：改為 "We propose the interpolating function μ(x) = 1 - e^{-√x}, which is **inspired by** (but not rigorously derived from) the Crooks fluctuation theorem"
  3. 建議選 B，因為選 A 可能推不出來

### HIGH（嚴重問題）

#### B4. μ₀ = H₀/c 的 dimensional analysis 不唯一
- **位置**：paper3_weak_field.tex:866-868 (Assumption A1)
- **問題**：論文誠實標記為 assumption，但 Section III.A (lines ~318-334) 的論證暗示這是 "natural"
  - 反例：μ₀ = √(H₀·Λ)/c 也沒有 Planck 物理，也給 O(10⁻²⁶) m⁻¹
- **修正建議**：在 dimensional analysis section 加入反例，明確說明 uniqueness 需要額外假設

#### B5. a₀–μ₀ 連結只有 13% — 過度宣傳
- **位置**：paper3_weak_field.tex:696-720, README.md line 19
- **問題**：
  - README 寫 "a₀ is derived to within 13%"，但 13% 不是 derivation
  - Main text (line 877) 誠實標記 "suggestive but unproven"，但 README 給人錯誤印象
- **修正建議**：
  1. README 改為：「a₀ = cH₀/(2π) gives 13% agreement — a suggestive but unproven coincidence」
  2. Main text 中加入其他理論也有類似 O(1) coincidence 的對照（e.g. Milgrom's a₀ ~ cH₀ 已知幾十年）

#### B6. Running μ 的三個論證都太弱
- **位置**：paper3_weak_field.tex:869-872 (Assumption A2)
- **問題**：
  - DPI + extensivity → μ ~ k：只是 consistency，不是推導
  - Dimensional transmutation：任何 running 都滿足
  - Modular flow T_mod = k/π：需要解釋為什麼 modular temperature 是對的 scale
- **修正建議**：
  1. 改為：「Three independent physical arguments **motivate** μ(k) = k, but a rigorous derivation from the microscopic theory remains open」
  2. 加一段：如果 μ(k) ≠ k（例如 μ(k) = k^{0.9}），CMB 預測如何變化？— 這會展示 robustness

#### B7. 缺 compiled PDF
- **修正建議**：同 A7

### MEDIUM（需要修但不致命）

#### B8. w(z=0) = 0.145 的物理解釋矛盾
- **位置**：paper3_weak_field.tex Table 3, lines ~735-767
- **問題**：
  - w > 0 意味著 DM 密度比 a⁻³ 稀釋更快
  - 但 DESI 偏好 w < -1/3（dark energy 方向）
  - 論文宣稱 w = 0.145 被 DESI 偏好，需要解釋這不矛盾
- **修正建議**：加一段解釋 DESI 測的是 effective dark energy w_DE，而 Khronon 的 w_K 是 dark matter equation of state，兩者混淆會導致 phantom crossing 的假象

#### B9. Fagin comparison 用 isothermal 近似
- **位置**：paper3_weak_field.tex:817-824
- **問題**：2.4σ tension 基於 isothermal 近似，但真正的 Khronon profile 可能不同
- **修正建議**：明確標記：「This prediction uses the isothermal approximation and should be revisited with full Khronon halo profiles」

#### B10. Paper 1/2 的引用格式不一致
- **問題**：Paper 1 有 Zenodo DOI，Paper 2 只有 "companion paper (2026)"
- **修正建議**：統一引用格式，Paper 2 也上 Zenodo

#### B11. Table 1 (comparison) 對 AeST 不完全公平
- **位置**：paper3_weak_field.tex lines ~253-270
- **問題**：列表偏向 Khronon（3 params vs AeST 5 params），但沒提 AeST 已有 Planck-level CMB fit
- **修正建議**：加一列「CMB precision」：AeST = Planck-compatible, Khronon = pending (requires BS2025 tensor field)

#### B12. Assumption 列表位置很好但缺一個
- **位置**：paper3_weak_field.tex:863-882
- **問題**：漏掉 A7：「The galactic sector J(Y) and cosmological sector K(Q) are decoupled — no cross-coupling terms in the action」
- **修正建議**：加 A7

---

## C. 兩個 Repo 的共同問題

| # | 問題 | 建議 |
|---|------|------|
| C1 | 沒有 PDF | 編譯後加入 repo |
| C2 | 沒有 compilation instructions | 加 Makefile 或 README 說明 |
| C3 | LICENSE 是 MIT（不適合論文） | 改為 CC-BY-4.0 |
| C4 | 沒有 .bib 檔（用 thebibliography） | 可以保留，但考慮改 .bib 方便引用管理 |
| C5 | README 沒有 Abstract | 加上論文 abstract 的英文版 |

---

## D. 修正優先級排序

### 第一輪（1-2 天，解決致命問題）

| 順序 | 問題 | 預估時間 | 效果 |
|------|------|---------|------|
| 1 | B2: 驗證/修正 BS2025 arXiv | 10 分鐘 | 消除不專業印象 |
| 2 | A1: Saturation → Lower bound geometry | 2 小時 | 消除核心矛盾 |
| 3 | A2: 加 Proposition（Σ = r_s/r 的唯一解） | 1 小時 | 正面 characterization |
| 4 | A3: 三路線改為 interpretations | 30 分鐘 | 誠實化 |
| 5 | B3: Crooks 降級為 inspired by | 30 分鐘 | 消除推導缺口 |
| 6 | C1: 編譯 PDF | 20 分鐘 | 基本專業度 |
| 7 | C3: LICENSE 改 CC-BY-4.0 | 5 分鐘 | 學術標準 |

### 第二輪（3-5 天，提升品質）

| 順序 | 問題 | 預估時間 |
|------|------|---------|
| 8 | A4: Pikovski 章節重寫 | 3 小時 |
| 9 | A5: Black hole table 重標記 | 1 小時 |
| 10 | B1: No-go 定理加強或誠實化 | 3 小時 |
| 11 | B5: README a₀ 措辭修正 | 15 分鐘 |
| 12 | B6: Running μ robustness analysis | 2 小時 |
| 13 | A6: GW error budget | 4 小時 |
| 14 | A8: 移除 [KNOWN]/[NEW] 標記 | 1 小時 |

### 第三輪（1 週，打磨）

所有 MEDIUM 問題 + PDF compilation + Zenodo DOI

---

## E. 寄送策略（修正後）

| 對象 | 需完成 | 預期反應 |
|------|--------|---------|
| **Blanchet & Skordis** (Paper 3) | B1-B3, B5, C1 | 「有趣的系統化分析，但 μ₀ 和 Crooks 需要更多工作。No-go 是有價值的。」正面 60% |
| **Dorau & Much** (Paper 2) | A1-A3, A4, C1 | 「Lower bound geometry 的概念有趣，但需要從 first principles 推導。」正面 40% |
| **Wilde** (Paper 1+2) | A1-A3, C1 | 「QI 技術正確但物理連結需要加強。」正面 50% |
| **Skordis & Zlosnik** (Paper 3) | B1-B3, B11 | 「比 AeST 更簡單但精度未驗證。」中性 50% |
