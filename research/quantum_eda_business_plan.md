# 量子 EDA 商業計劃：從 Petz Recovery 到產品

**作者：Sheng-Kai Huang**
**日期：2026-03-16**
**版本：v1.0**

---

## 一、Executive Summary

基於 Paper 1《The Arrow of Time from Petz Recovery》的數學框架，開發量子電腦的 **Noise Intelligence Engine**——一個嵌入現有量子 compiler stack 的中間層，提供基於嚴格信息論的噪聲管理。

核心技術優勢：Petz recovery map 是**唯一**滿足 Bayesian consistency 的 retrodiction functor（Parzygnat-Buscemi 2023），所有好的 quantum error correction decoder 都必須逼近 Petz map。這個唯一性定理是護城河。

**產品定位**：量子版的 Static Timing Analysis engine——不取代 compiler（Qiskit/TKET/Classiq），而是提供 noise-aware 決策層。

---

## 二、市場概況（2026）

### 2.1 市場規模

| 指標 | 數值 | 來源 |
|------|------|------|
| 量子軟體市場（2026） | $12.5 億 | Coherent MI |
| 量子軟體市場（2033） | $47.5 億（CAGR ~21%） | Coherent MI |
| 整體量子計算市場（2035） | $194 億（CAGR 29.7%） | Precedence Research |
| 量子應用開發軟體（2030） | >$100 億（CAGR >40%） | SpinQ |
| 2025 H1 量子新創融資 | $12.5 億（YoY 2x） | 產業報告 |

### 2.2 競爭者地圖

#### 量子電路 EDA（compiler 層）

| 公司 | 產品 | 融資/估值 | 定位 | 弱點 |
|------|------|----------|------|------|
| IBM | Qiskit SDK v2.2 | ~$1B 累計收入 | 效能最強（83x 快於 TKET） | 綁定 IBM 硬體 |
| Quantinuum | TKET | 估值 $100 億 | 跨平台 compiler | 效能不如 Qiskit |
| Classiq | Qmod + 自動合成 | Series C $2 億+（最大量子軟體輪） | 高階描述 → 自動電路 | 無 noise-aware 理論基礎 |
| Google | Cirq + AlphaQubit | — | 研究導向 + ML decoder | 僅支援 Google 硬體 |
| Q-CTRL | Fire Opal | 融資 $1.33 億 | 錯誤抑制中間層 | 無 QEC decoder 功能 |
| Riverlane | Deltakit + LCD | 融資 $1.25 億 | QEC 專家 | 無 compiler 功能 |
| Rigetti | pyQuil/Quilc | 上市公司 | 全棧（硬體+軟體） | 規模較小 |
| Horizon QC | Beryllium | 估值 ~$5 億（Nasdaq: HQ） | 硬體無關高階語言 | 早期階段 |

#### 量子晶片 EDA（物理層）

| 公司 | 產品 | 定位 |
|------|------|------|
| Keysight | QuantumPro EM, ADS 2026 | 超導量子晶片設計（傳統 EDA 巨頭） |
| Nanoacademic + Kothar | 首套 Quantum EDA | 量子晶片 TCAD（spin + 超導） |
| MQT（TU Munich） | 開源工具集 | 學術全棧工具 |

### 2.3 市場缺口

目前**所有**量子 compiler 都缺少一個東西：**基於嚴格信息論的統一 noise 管理層**。

現狀的三個問題：

1. **Noise 指標碎片化**：gate fidelity、diamond distance、quantum volume、T1/T2 各自為政，無統一度量
2. **誤差預算靠 heuristic**：沒有嚴格的逐層分解定理，無法做 timing budget 式的分配
3. **Decoder 設計與 compiler 脫節**：QEC decoder 獨立開發，不整合進 compilation flow

---

## 三、技術基礎：Paper 1 的五個 EDA 應用

### 3.1 τ = 1-F 作為統一 noise 指標

**來源**：Paper 1 Eq. (9), 等價鏈 Eq. (11)

**現狀**：IBM 用 layer fidelity、各家用 ad-hoc fidelity metrics、Q-fid 用 LSTM 預測。無統一標準。

**τ 的優勢**：
- 直接綁定 recoverability（不只是 error rate）
- 明確操作意義：τ = 0 → 完美恢復，τ = 1 → 完全信息損失
- 可作為 compiler 最佳化的 objective function（取代 gate count）
- Parzygnat-Buscemi 唯一性定理保證**沒有其他 Bayesian-consistent 選擇**

**EDA 應用**：compiler pass 計算每個 candidate qubit mapping 的 τ，選 minimum τ 方案。

### 3.2 飽和定理診斷 channel 品質

**來源**：Paper 1 Theorem 4（Saturation）

**核心**：$F^2 = \exp(-\Delta D)$ 當且僅當 $[\omega,\tau]=0$ 且 likelihood ratio 為常數。

**EDA 應用**：
- 一個數字判斷某 gate 是否 exactly correctable
- 不需要跑 decoder 就能從 channel statistics 判斷 decoder performance
- 硬體廠商可報告 τ per gate/per qubit pair 作為標準化品質指標

### 3.3 組合次可加性做誤差預算

**來源**：Paper 1 Theorem 5（Composition）

**核心**：$\sqrt{\tau_{12}} \leq \sqrt{\tau_1} + \sqrt{\tau_2^{\text{eff}}}$

**EDA 應用**：
- **嚴格的逐層誤差分解**——類比 VLSI Static Timing Analysis
- 分配總 τ budget 到各層，每層 compiler pass 在 budget 內最佳化
- 若累積 τ 超過閾值，提前終止並換策略
- n-channel 推廣（Corollary S19）支持任意深度電路

### 3.4 Decoder Hierarchy 指導 QEC decoder 設計

**來源**：Paper 1 Observation 2 + Retrodiction gap δ_D

**核心**：Petz map 是唯一 Bayesian-consistent decoder → 所有好的 decoder 都逼近 Petz map → retrodiction gap δ_D 排序 decoder 品質。

**EDA 應用**：
- 自動生成 Petz recovery circuit（不需要 classical decoder）
- 用 δ_D 作為 decoder 品質的統一評分標準
- 對不同架構的 decoder 做公平比較

### 3.5 等價鏈作為跨領域翻譯字典

**來源**：Paper 1 等價鏈 Eq. (11)

**核心**：τ=0 ⟺ Σ=0 ⟺ I(A;E|B)=0 ⟺ QMC

**EDA 應用**：
- 接受任何格式的 noise data（gate fidelity, T1/T2, process matrix）
- 用等價鏈轉換為統一的 τ 指標
- 報告最終電路品質時可轉回任何所需 metric
- 類比：classical EDA 的 Liberty format（統一 delay model）

---

## 四、產品架構

### 4.1 在 stack 中的位置

```
Layer 4: 應用層       [用戶的量子演算法]
            ↓
Layer 3: Compiler      [Qiskit / TKET / Classiq — 電路合成 + routing]
            ↓
Layer 2: τ-Engine ★    [Petz Noise Intelligence — 本產品]
            ↓
Layer 1: Hardware      [IBM / Google / Quantinuum / IonQ QPU]
```

τ-Engine **不取代** Qiskit 或 Classiq，而是**嵌入它們之間**，提供 noise-aware 決策層。

### 4.2 核心模組

#### 模組 A：τ-Profiler（硬體 noise 建模）

- **輸入**：硬體校準數據（T1/T2、gate fidelity、process matrix）
- **處理**：用等價鏈轉換為統一的 τ 指標
- **輸出**：每個 qubit / gate pair 的 τ noise map
- **技術基礎**：Paper 1 等價鏈 + Gibbs 特化

#### 模組 B：τ-Budget Allocator（誤差預算分配）

- **輸入**：目標總 fidelity + 電路結構
- **處理**：用 composition theorem 逐層分配 τ 預算
- **輸出**：每層的 τ 上界，作為 constraint 傳給 compiler
- **技術基礎**：Paper 1 Theorem 5 + n-channel Corollary
- **類比**：量子版的 VLSI Static Timing Analysis

#### 模組 C：τ-Optimizer（noise-aware 編譯優化）

- **輸入**：candidate 電路 + τ noise map + τ 預算
- **處理**：在所有 candidate mapping/routing 中選 minimum τ 的方案
- **輸出**：最佳化電路
- **技術基礎**：τ 作為 objective function（取代 gate count）

#### 模組 D：Decoder Synthesizer（QEC decoder 生成）

- **輸入**：error correction code + noise model
- **處理**：生成 Petz recovery circuit → 優化為目標硬體
- **輸出**：hardware-specific decoder 電路
- **技術基礎**：Petz map 的 Kraus representation + saturation theorem 判斷最優性
- **殺手功能**：不需要 classical decoder（MWPM/UF），直接用量子電路做 recovery

#### 模組 E：τ-Dashboard（診斷與報告）

- 視覺化 τ 沿電路的累積曲線
- 標記 saturation / non-saturation 的 gate
- 提供 decoder performance score（retrodiction gap δ_D）
- 匯出報告（PDF / JSON / API）

---

## 五、開發路線圖

### Phase 1：核心引擎 + PoC（Month 1-6）

| 月 | 里程碑 | 交付物 |
|----|--------|--------|
| 1-2 | τ-Profiler MVP | Python SDK，支援 IBM/Google calibration data → τ map |
| 3-4 | τ-Budget Allocator | 逐層 τ 分配器，驗證 composition theorem 在實際電路上 |
| 5-6 | Qiskit plugin PoC | 整合進 Qiskit transpiler 作為自訂 pass |

**PoC 驗證方式**：
- 用 IBM Brisbane/Sherbrooke 校準數據
- 比較 τ-optimized routing vs standard Qiskit routing 的 fidelity
- 目標：在 10-50 qubit 電路上實現 **5-15% fidelity improvement**

**所需資源**：1-2 開發者 + IBM Quantum cloud access

### Phase 2：整合與驗證（Month 7-12）

| 月 | 里程碑 | 交付物 |
|----|--------|--------|
| 7-8 | TKET 整合 | pytket 後端，支援跨平台 |
| 9-10 | Decoder Synthesizer v1 | 對 surface code d=3,5,7 生成 Petz decoder |
| 11-12 | Benchmark paper | 在 arXiv 發表 τ-Engine vs 現有工具的定量比較 |

**關鍵驗證**：
- τ-Budget 的 decoder-independent α 預測（Paper 1 Observation 1）
- Petz decoder vs MWPM/UF 的 fidelity 對比
- 在真實硬體上的 end-to-end 測試

### Phase 3：產品化 + 首批客戶（Month 13-18）

| 月 | 里程碑 | 交付物 |
|----|--------|--------|
| 13-15 | SaaS 平台 | Web API + Dashboard + 文件 |
| 16-18 | 首批企業 pilot | 與 2-3 家量子硬體/雲端公司合作 |

### Phase 4：擴展 + 規模化（Month 19-36）

- 支援所有主流 QPU（超導、離子阱、中性原子、光子）
- 整合進 Classiq 高階合成流程
- 推出 on-premise 版本給國防/金融客戶
- 考慮 Decoder Synthesizer 獨立授權

---

## 六、商業模式與定價

### 6.1 收入來源（三層）

#### SaaS 訂閱（主營收）

| Tier | 功能 | 月費 |
|------|------|------|
| Community | τ-Profiler（單 backend）、基本 Dashboard | 免費 |
| Pro | 全模組、多 backend、API access | $500/月 |
| Enterprise | 自訂 noise model、dedicated support、on-premise 選項 | $5,000-$20,000/月 |

#### IP 授權（高利潤）

- 對硬體廠商授權 Decoder Synthesizer
- 模式：$100K-$500K upfront + 5-10% royalty per QPU sold
- 參考：Riverlane 對 Deltaflow 的授權模式

#### 專業服務（早期現金流）

- Noise characterization consulting：$50K-$200K per engagement
- 客製化 Petz decoder 設計：$100K-$500K
- 學術合作 + 政府合約（DARPA / NSF / 科技部）

### 6.2 目標客戶

| 客戶類型 | 需求 | 價值主張 | 預估合約規模 |
|---------|------|---------|------------|
| 量子硬體公司（IBM, IonQ, Quantinuum） | 提升 QPU 可用性 | τ-Engine 讓硬體 fidelity 提升 5-15% | $100K-$500K/年 |
| 量子雲端平台（AWS Braket, Azure Quantum） | 差異化服務 | 整合 τ noise profiling 作為增值功能 | $200K-$1M/年 |
| 量子應用開發者 | 寫出更好的量子程式 | noise-aware compilation 一鍵完成 | $6K-$12K/年 |
| 國防/金融 | 對 QEC 有急迫需求 | Petz decoder 的理論最優性保證 | $500K-$2M/年 |
| 學術機構 | 研究工具 | Community tier 免費 + paper 引用 | 免費（品牌建設） |

### 6.3 收入預測

| 年 | ARR | 客戶數 | 主要來源 |
|----|-----|--------|---------|
| Y1 | $200K | 5-10 | 專業服務 + 早期 pilot |
| Y2 | $1M | 20-50 | Pro subscriptions + 首批 IP 授權 |
| Y3 | $5M | 100+ | Enterprise + IP royalties |
| Y5 | $20M+ | 500+ | 全面規模化 |

---

## 七、競爭優勢（Moat）

### 7.1 專利護城河

- τ-Budget Allocator（composition theorem 的 EDA 應用）具有專利潛力
- Petz Decoder Synthesizer（自動生成 Petz recovery circuit）具有專利潛力
- Paper 1 的 Theorem 4（飽和定理）和 Theorem 5（組合定理）是原創數學結果

### 7.2 理論基礎的唯一性

- Petz map 是**唯一**滿足 Bayesian consistency 的 retrodiction functor（Parzygnat-Buscemi 2023）
- 所有好的 decoder 都必須逼近 Petz map——這是定理，不是猜想
- 任何競爭者如果要做同樣的事，必須用 Petz map——進入了我們的專利範圍

### 7.3 先行者優勢

- 目前**沒有**任何量子 EDA 產品使用 Petz recovery 的信息論框架做 noise management
- Q-CTRL 做 noise suppression 但沒有信息論基礎
- Riverlane 做 QEC decoder 但沒有統一的 τ 指標

### 7.4 網路效應

- 越多硬體用 τ 指標校準 → 越多 compiler 整合 τ-Engine → τ 成為產業標準指標
- Community tier 的免費用戶建立生態系 + 鎖定效應

---

## 八、團隊需求

### 核心團隊（Phase 1-2）

| 角色 | 人數 | 技能需求 |
|------|------|---------|
| 理論物理/量子信息 | 1 | Petz map、QEC、量子熱力學（創辦人） |
| 量子軟體工程師 | 2 | Qiskit/TKET 內部架構、Python、Rust |
| 全端工程師 | 1 | Dashboard、API、SaaS 基礎設施 |

### 擴展團隊（Phase 3-4）

| 角色 | 人數 | 技能需求 |
|------|------|---------|
| 業務開發 | 1 | 量子產業人脈、B2B SaaS 銷售 |
| 應用工程師 | 2 | 客戶整合、技術支援 |
| ML 工程師 | 1 | noise model 學習、decoder 優化 |

---

## 九、融資規劃

| 輪次 | 時機 | 金額 | 用途 | 預估估值 |
|------|------|------|------|---------|
| Pre-seed | Phase 1 PoC 完成後 | $500K-$1M | 2-3 人核心團隊 + 雲端 QPU access | $3M-$5M |
| Seed | 首批付費客戶後 | $3M-$5M | 10 人團隊 + Qiskit/TKET 深度整合 | $15M-$25M |
| Series A | ARR $1M 後 | $15M-$25M | 規模化 + 國際擴展 | $75M-$150M |

**參考估值**：
- Classiq 在 Series B 時 ARR < $5M 但估值 $3 億+（量子軟體高溢價）
- Riverlane Series C $7,500 萬，估值 ~$5 億
- 量子軟體公司的估值倍數遠高於傳統 SaaS

### 潛在投資者

| 類型 | 名稱 | 理由 |
|------|------|------|
| 量子專注 VC | Entree Capital, QIC Ventures | 已投 Classiq |
| 策略投資者 | AMD Ventures, Qualcomm Ventures | 已投 Classiq，需要 EDA |
| 硬體廠商 | IonQ, IBM Ventures | 需要 noise management 層 |
| 台灣/亞太 | 國發基金, 鴻海研究院 | 台灣量子戰略 |

---

## 十、風險與對策

| 風險 | 機率 | 影響 | 對策 |
|------|------|------|------|
| 大廠自建類似功能 | 高 | 高 | 走跨平台中立路線（像 TKET vs Qiskit）；申請專利 |
| τ 改進量在實際硬體上不顯著 | 中 | 致命 | Phase 1 先 PoC 驗證，不顯著就 pivot 到純 consulting |
| QEC 時代來得比預期慢 | 中 | 中 | Decoder Synthesizer 不急推，先靠 τ-Profiler + Budget 賺錢 |
| Classiq/Riverlane 合併覆蓋此領域 | 低 | 高 | 主動尋求成為他們的 plugin 或被收購 |
| 論文未被社群接受 | 低 | 中 | 已有獨立驗證（Buscemi 2024）+ 實驗實現（Png 2025, Singh 2025） |
| 專利申請失敗 | 中 | 中 | 核心價值在 know-how + 先行者優勢，非僅靠專利 |

---

## 十一、里程碑與 Go/No-Go 決策點

### Decision Gate 1（Month 6）：PoC 結果

- **Go**：τ-optimized routing 在 IBM 硬體上實現 ≥5% fidelity improvement
- **Pivot**：improvement < 5% → 轉向 consulting + 學術工具
- **No-Go**：τ 指標與 fidelity 無相關 → 重新評估整個方向

### Decision Gate 2（Month 12）：市場驗證

- **Go**：至少 1 家付費 pilot 客戶 + benchmark paper 被引用
- **Pivot**：無付費客戶但有學術影響 → 轉向開源 + grant funding
- **No-Go**：無客戶且無影響 → 關閉

### Decision Gate 3（Month 18）：規模化準備

- **Go**：ARR > $200K + 3+ 付費客戶 → 融 Seed 輪
- **Pivot**：客戶少但單筆大 → 轉向 IP 授權模式
- **No-Go**：無法達到 product-market fit → 考慮被收購或關閉

---

## 十二、第一步行動計劃（立即可執行）

### 本週

1. 將 Paper 1 的 `numerical/` 目錄重構為 Python SDK 核心（τ-Profiler 原型）
2. 下載 IBM Brisbane 最新校準數據，計算所有 qubit pair 的 τ map

### 本月

3. 用 τ map 做 noise-aware routing，與 Qiskit default transpiler 對比 fidelity
4. 量化 improvement（目標 ≥5%）

### 下月

5. 在 arXiv 發表 "τ as a Universal Noise Metric for Quantum Circuit Optimization"
6. 開源 τ-Profiler core（建立用戶基礎 + 學術 credibility）
7. 聯繫 Classiq / Q-CTRL 探討合作可能

### Q2 2026

8. Qiskit plugin 上線 GitHub
9. 申請 USPTO 臨時專利（τ-Budget Allocator + Decoder Synthesizer）
10. 開始融 Pre-seed

---

## 附錄 A：技術細節——τ_std vs τ_rot 在 EDA 中的處理

### 問題

Paper 1 中存在兩個版本的 τ：
- τ_std = 1 - F(standard Petz)：滿足組合定理，有唯一性
- τ_rot = 1 - F(rotated Petz)：滿足更緊的 bound

### EDA 的解決方案

**全程使用 τ_std。** 理由：

1. **Composability**：EDA 的核心用例是 error budget allocation，需要組合定理 → 只有 τ_std 滿足
2. **唯一性**：Parzygnat-Buscemi 定理保證 τ_std 是唯一 Bayesian-consistent 選擇 → 標準化所需
3. **閉式計算**：τ_std 直接從 Petz map 公式算，不需要數值積分 → 工程效率
4. **差異可忽略**：在當前硬體操作區間（gate fidelity > 99%），τ_std ≈ τ_rot + O(τ²)

| 特性 | τ_std | τ_rot | EDA 需要？ |
|------|-------|-------|-----------|
| 組合定理 | ✓ | ✗ | **必須** |
| 唯一性（Bayesian） | ✓ | ✗ | **必須** |
| 閉式計算 | ✓ | ✗ | 強烈偏好 |
| 最緊 bound | ✗ | ✓ | 不需要 |

---

## 附錄 B：Paper 1 技術與 EDA 模組對應表

| Paper 1 結果 | 數學內容 | EDA 模組 | 功能 |
|-------------|---------|---------|------|
| 等價鏈 Eq. (11) | τ=0 ⟺ Σ=0 ⟺ QMC | τ-Profiler | 跨格式 noise 轉換 |
| Theorem 4（飽和） | F²=exp(-ΔD) iff 充分統計 | τ-Profiler | channel 品質一鍵診斷 |
| Theorem 5（組合） | √τ₁₂ ≤ √τ₁+√τ₂ | τ-Budget Allocator | 逐層誤差分配 |
| Corollary S19 | n-channel 推廣 | τ-Budget Allocator | 任意深度電路 |
| Observation 2 | decoder hierarchy ∝ δ_D | Decoder Synthesizer | decoder 品質評分 |
| Master chain Eq. (6) | -log F² ≤ I(A;E|B) ≤ Σ ≤ ΔD | τ-Dashboard | 多指標報告 |
| Corollary 1（unitary） | τ=0 for unitary | τ-Optimizer | 識別可省略的 recovery |
| Crooks connection | Petz = reverse channel | Decoder Synthesizer | 物理 decoder 設計 |

---

## 附錄 C：量子 EDA 產業參考資料

### 主要公司近期動態（2025-2026）

- **Classiq** Series C $2 億+（2025-2026），AMD/Qualcomm/IonQ/SoftBank 參投
- **Riverlane** LCD decoder 發表於 Nature Communications（2025.12），4x fewer qubits
- **Q-CTRL** DARPA 合約 A$3,800 萬（2025）
- **Quantinuum** Quantum Volume 33,554,432（2025.09），估值 $100 億
- **Google** AlphaQubit decoder 30% error reduction（2025）
- **Keysight** Quantum System Analysis tool in ADS 2026
- **Horizon QC** 首家量子軟體公司擁有自己的量子電腦（2025.12）
- **IBM** Kookaburra processor roadmap: 首個 qLDPC QEC 模組（2026）

### 商業模式參考

| 模式 | 代表公司 | 備註 |
|------|---------|------|
| QaaS（按次計費） | AWS Braket: $0.30/task + $0.10/shot (IonQ) | 我們不走此路 |
| SaaS 訂閱 | Classiq, Q-CTRL, Strangeworks | **主要模式** |
| IP 授權 | Riverlane ($100K-$500K + royalty) | **次要模式** |
| 硬體綁定 | IBM, Quantinuum | 不適用 |
| 政府合約 | Q-CTRL (DARPA), IonQ (DoD) | 早期現金流 |
| 上市 | Horizon QC (Nasdaq: HQ), IonQ | 長期目標 |
