# 虛數單位 = Gauge Forces：Σ 的方向結構
## 2026-03-19

## 核心觀察

Division algebra 的虛數單位數量精確對應 gauge group 的 generators：

| 數域 | 虛數單位數 | Gauge Group | Generators 數 |
|------|:---------:|------------|:-------------:|
| R | 0 | — | 0（只有時間箭頭）|
| C | 1 (i) | U(1) | 1 |
| H | 3 (i,j,k) | SU(2) | 3 |
| O | 7 (e₁...e₇) | G₂ ⊃ SU(3) | 7 → 8 |

## 在 Σ 框架中的新詮釋

### 虛數 = Σ 的相位方向

Complex Khronon φ = f·e^{iθ}：
- 實部 f（amplitude）→ 重力（Σ_grav = 2 ln Q, Q 由 f 決定）
- 虛部 θ（phase）→ 電磁（Σ_EM from phase gradient）
- **i 是重力和電磁之間的 90° 旋轉**

### 新猜想：Σ 的方向數 = 虛數數量

Σ = D(ρ‖σ) 是一個「資訊距離」。在不同的代數結構上，這個距離有不同數量的「不可約方向」：

```
1-qubit (C): Σ 有 1 個相位方向 → 1 個 gauge generator → U(1)
2-qubit (H): Σ 有 3 個糾纏方向 → 3 個 generators → SU(2)
3-qubit (O): Σ 有 7 個糾纏方向 → 7 → SU(3)?
```

### 精確地說

在 Level 2 中（今天驗證的）：
- QRE chain rule: D(ρ‖I/4) = D(ρ‖G(ρ)) + D(G(ρ)‖I/4)
- D(ρ‖G(ρ)) = SU(2)-breaking 部分（3 個方向，對應 i,j,k）
- D(G(ρ)‖I/4) = SU(2)-invariant 部分

Singlet |Ψ⁻⟩: D(ρ‖G(ρ)) = 0（所有 3 個方向都 invariant）
Triplet |Φ⁺⟩: D(ρ‖G(ρ)) = ln 3 ≈ 1.099（3 個方向都 broken）

**三個方向 = 三個虛數單位 i, j, k = SU(2) 的三個 generators。**

## 虛數的物理意義

### 為什麼 i 是必需的

Paper 6 Theorem 1 的 5 條路線都指向同一件事：Σ 需要 complex numbers。

但更深的原因（Tomita-Takesaki）：
- Modular flow: σ_t(A) = Δ^{it} A Δ^{-it}
- 這裡的 **it** 中的 i 不是裝飾
- 沒有 i → Δ^{it} 無法定義 → modular flow 不是 group（只是 semigroup）
- 沒有 modular flow → 沒有 QRE 的完整結構

**i 讓時間從半群（只能向前）變成群（可以定義逆）。但宏觀上 Σ ≥ 0 仍然保證時間箭頭。**

### 不交換性 → Non-Abelian

ℂ: ij = ji（不適用，只有一個虛數）
ℍ: ij ≠ ji（ij = k, ji = -k）→ SU(2) 是 non-Abelian
𝕆: 不交換 AND 不結合（(ab)c ≠ a(bc)）→ 更深的結構

**力的 Abelian vs non-Abelian 性質直接來自數域的交換性。**

### 非結合性 → Confinement？

八元數不結合：(e₁e₂)e₃ ≠ e₁(e₂e₃)

推測（Furey 等人方向）：
- 夸克不能單獨存在（confinement）
- 八元數的元素不能「獨立結合」（non-associativity）
- 兩者可能有深層聯繫
- **confinement = Σ 在八元數方向上的非結合性阻止了「部分恢復」？**

## 歸屬

| 觀察 | 誰 | 我們的？ |
|------|-----|---------|
| 虛數數量 = gauge generators | Dixon (1994), Furey (2012+) | 否 |
| 不交換性 → non-Abelian | 已知（教科書） | 否 |
| 非結合性 → confinement? | Furey, Boyle（推測） | 否 |
| **i = 重力-EM 的 90° 旋轉** | **無人這樣說過** | **新（我們的詮釋）** |
| **Σ 的方向數 = 虛數數量** | **無人計算過** | **新猜想（待驗證）** |
| **modular flow 需要 i → 時間從半群變群** | 部分已知（Connes），但 Σ 連結新 | **新詮釋** |

## 計算結果（2026-03-19 完成）

### Level 2：精確匹配 ✅
D_break = ln(2j+1) = ln 3 for triplet (j=1)
3 個方向 = dim V₁ = dim Im(H) = {i,j,k}
原因：su(2) ≅ Im(H)（Lie algebra 同構）

### Level 3：不匹配 ❌
D_break = ln 4 for j=3/2 symmetric subspace
4 ≠ 7（O 的虛數單位數）
原因：SU(2)_diagonal 只給 j=3/2 → dim 4。需要 G₂ 才能得到 7。

### 修正後的理解
猜想「Σ 方向數 = 虛數數量」在代數層面（QRE chain rule）只到 Level 2 成立。
Level 3 的 1,3,7 對應住在 Hopf fibration 的 fiber 維度（拓撲），不是 QRE irrep（代數）。

Level 2 是代數+拓撲都匹配的 sweet spot：su(2) ≅ Im(H) 且 Hopf fiber S³ ≅ SU(2)。

## 對論文的意義

如果「Σ 的方向數 = 虛數數量」被驗證，可以加入 Paper 7 或 Paper 9：

"The number of independent directions in which the quantum relative
entropy D(ρ‖σ) can be decomposed on an n-qubit system equals the
number of imaginary units of the corresponding division algebra.
This provides an information-theoretic explanation for why the gauge
group at each level has exactly dim(Im(A)) generators, where A is
the division algebra at that level: 1 for C (electromagnetism),
3 for H (weak force), and 7 for O (related to strong force)."
