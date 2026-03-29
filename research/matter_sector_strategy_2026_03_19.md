# Matter Sector 攻略：從 Spacetime 到 Standard Model
## 2026-03-19
## Sheng-Kai Huang

### 一、現狀診斷

Σ = D(ρ_spacetime ‖ ρ_matter) 已完成 spacetime sector（59% 覆蓋率）。
Matter sector（Higgs, SU(2), SU(3), fermions, 3 代）仍然空白。

已確認失敗：
- ❌ Khronon = Higgs（shift symmetry vs fixed VEV，結構不相容）
- ❌ Khronon = Relaxion（永遠滾動 vs 必須停止，動力學不相容）
- ❌ α = 1/137（管轄範圍外）
- ❌ 費米子從純玻色湧現（所有機制都失敗）
- ❌ Σ = spectral action（S_vN ≠ QRE，維度不匹配）

### 二、三重交集結構

三條路線（Connes, CW Portal, Division Algebra）交匯於：

**Scale-invariant spectral action on J₃(O)**

- Connes: A_F = C ⊕ H ⊕ M₃(C) 是 J₃(O) 的 maximal associative subalgebra
- Furey: C⊗H⊗O 的 ladder operators = fermion representations
- Szangolies: 3-qubit entanglement = O 的 information content
- Farnsworth 2025: nonassociative spectral geometry（G₂×G₂，待推廣到 SM）
- CCSvS 2019: S_vN of spectral triple = spectral action
- Chamseddine-Connes 2006: scale-invariant spectral action 存在（dilaton = Khronon）

### 三、Σ 的 Division Algebra 階梯

```
Level 0 (R): Σ ≥ 0 (DPI)           → 時間箭頭        ✅ Paper 1
Level 1 (C): Σ 需要 complex        → U(1) → EM       ✅ Paper 6
Level 2 (H): Σ 的 2-body 糾纏結構  → SU(2) → 弱力    🔨 Paper 7
Level 3 (O): Σ 的 3-body 糾纏結構  → SU(3) → 強力    🔨 Paper 8
Level 4 (J₃(O)): Σ 的 extremal 結構 → 3 代 + 質量比   🔨 Paper 10
```

每一層不是新假設，而是 Σ 在越來越複雜的系統上的自然展開。

### 四、Scale-Invariant CW + Khronon Portal

核心鏈：
1. 理論起點：classically scale-invariant（Chamseddine-Connes 2006）
2. Khronon ghost condensation → M_Pl（已有 Papers 2-4）
3. Kinetic portal λ_p|H|²(∂φ)²/Λ² 傳遞 scale breaking
4. CW one-loop → v ≠ 0
5. v ~ M_Pl exp(-c/λ_p)（指數壓制 → hierarchy 自然）

關鍵障礙：μ₀ = H₀/c ~ 10⁻³³ eV 太小。需要 M_GC ~ M_Pl（UV condensation scale）。
Running μ(k) = k 可能連接：μ(k_Pl) ~ M_Pl, μ(k_Hubble) ~ H₀/c

### 五、Paper 9 計劃（最高優先）

目標：δΣ on A_F ≈ δ(spectral action)

Toy model 第一步：S¹ × M₂(C)（circle × 2×2 matrices）
- 計算 S_vN (CCSvS formula)
- 計算 Σ = QRE between two Gibbs states
- 比較 first variations
- Go/No-Go at Week 8

如果成功：推廣到 full A_F = C ⊕ H ⊕ M₃(C) → SM + gravity from δΣ = 0

### 六、完整論文路線圖

```
已完成：
  Paper 1:  τ = 1-F                    [GitHub + Zenodo]
  Paper 1b: 坍塌 = Petz failure        [6頁, 0 errors]
  Paper 2:  Σ_grav → 指數度規          [arXiv-ready]
  Paper 3:  暗物質                      [需 DBI]
  Paper 4:  μ₀ = H₀/c                  [14頁]
  Paper 5:  Complementary Uncertainty   [6頁]
  Paper 6:  EM 統一                     [5頁, 0 errors]
  Supplement 1: Paper 1b 數學           [8頁]
  Supplement 2: Paper 6 數學            [7頁]

待完成：
  Paper 7:  SU(2) via CW + portal      [3-6 月]
  Paper 8:  SU(3) via O               [6-12 月，等 Farnsworth]
  Paper 9:  Σ on A_F = spectral action [3-6 月，★★★最高優先]
  Paper 10: 費米子質量 J₃(O) + Σ       [6-12 月]
  Paper 11: Grand finale               [12-18 月]
```

### 七、Connes 路線具體計劃

必讀論文（前 5）：
1. CCSvS 2018 (arXiv:1809.02944)
2. Dong-Khalkhali-van Suijlekom 2019 (arXiv:1903.09624)
3. Dorau-Much 2025 (arXiv:2510.24491)
4. Van Suijlekom 教科書 2nd ed (2024)
5. Chamseddine-Connes 1996 (hep-th/9606001)

數學問題（按順序）：
1. 在 A_F 上定義 Σ_F
2. δΣ_F vs δ(spectral action) 比較
3. M × F 推廣
4. SM 場方程 from δΣ = 0
5. Second variation → stability

失敗機率：~35% regularization mismatch, ~25% inner fluctuation, ~20% Lorentzian

### 八、CW Portal 具體計劃

必讀論文（前 5）：
1. Trautner 2025 (arXiv:2502.09699) Custodial Naturalness
2. de Boer 2025 (arXiv:2507.22980) Hidden Sector
3. de Boer 2025 (arXiv:2510.12882) Gravity + Hierarchy
4. arXiv:2510.00686 Kinetic Higgs Portal
5. Frasca 2025 (arXiv:2408.00093) Non-perturbative EW

計算（按順序）：
A. Portal on ghost condensation background 展開
B. Effective Higgs mass δm_H²
C. CW one-loop potential
D. v 的解析表達式
E. m_H = 125 GeV 條件

Go/No-Go checkpoints 在 Week 2, 4, 6, 8

### 九、Division Algebra 具體計劃

必讀論文（前 4）：
1. Szangolies 2025 (arXiv:2512.17328)
2. Furey 2025 (arXiv:2505.07923)
3. Boyle 2020 (arXiv:2006.16265)
4. Singh 2025 (arXiv:2508.10131)

主攻 Route B（entanglement → SU(2)）：
- 2-qubit 的 Σ = D(ρ_AB ‖ σ_A ⊗ σ_B) 的 SU(2) 結構
- quaternionic Hopf fibration 作為糾纏幾何
- 20 週路線圖

### 十、三條路線交匯的精確位置

交匯點 1：Σ 選擇 C（Paper 6 Thm 1）= Krasnov 的 complex direction = Szangolies 降維
交匯點 2：Hopf equivariance group = A_F automorphism group
交匯點 3：SO(8) triality = J₃(O) 本徵值 = 3 代
三重交集：Scale-invariant spectral action on J₃(O) with Khronon as dilaton

### 十一、已確認失敗的記錄

1. Khronon doublet = Higgs：shift symmetry vs fixed VEV 不相容
2. Khronon = Relaxion：永遠滾動 vs 停止，5 個致命矛盾
3. Hierarchy from Khronon shift symmetry：保護 Khronon 不保護 Higgs
4. Σ = spectral action（直接等式）：S_vN ≠ QRE

### 十二、關鍵文獻總清單（新增 2026-03-19）

- Chamseddine-Connes 2006: Scale-invariant spectral action
- CCSvS 2018 (arXiv:1809.02944): entropy = spectral action
- Dong-Khalkhali-van Suijlekom 2019 (arXiv:1903.09624): chemical potential
- Boyle-Farnsworth 2018 (arXiv:1604.00847): differential graded algebra
- Boyle-Farnsworth 2020 (arXiv:1910.11888): Jordan geometry
- Farnsworth 2025 (arXiv:2506.21496): nonassociative spectral geometry
- Szangolies 2025 (arXiv:2512.17328): 3-qubit → SM
- Furey 2025 (arXiv:2505.07923): division algebra superalgebra
- Singh 2025 (arXiv:2508.10131): fermion mass ratios from J₃(O)
- Krasnov 2019 (arXiv:1912.11282): SO(9) characterization
- Gourlay-Gresnigt 2026 (arXiv:2601.07857): Cl(10) + S₃
- Trautner 2025 (arXiv:2502.09699): Custodial Naturalness
- de Boer 2025 (arXiv:2510.12882): Gravity + Hierarchy
- arXiv:2510.00686: Kinetic Higgs Portal
- Frasca 2025 (arXiv:2408.00093): non-perturbative EW via portal
- Nieuviarts 2025 (arXiv:2512.15450): twisted spectral triple
- Bianconi 2025 (PRD 111, 066001): gravity from entropy
- Oda 2025 (arXiv:2509.23648): GR from ghost condensation
