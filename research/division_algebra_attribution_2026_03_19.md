# Division Algebra 階梯：歸屬與貢獻分析
## 2026-03-19

## 歸屬表

### 數學基礎（不是我們的）

| 發現 | 誰 | 年份 |
|------|-----|------|
| 只有 4 種 division algebra (R,C,H,O) | Hurwitz | 1898 |
| Hopf fibrations S³→S², S⁷→S⁴, S¹⁵→S⁸ | Hopf | 1931 |
| 只有這 4 種 Hopf fibration（Adams 定理）| Adams | 1960 |
| 八元數 → 夸克 SU(3) color | Günaydin & Gürsey | 1973 |
| T = R⊗C⊗H⊗O → 粒子物理 | Dixon | 1994 |
| A_F = C⊕H⊕M₃(C) → SM from spectral triple | Connes | 1996 |
| C⊗H⊗O minimal ideals → 一代 SM fermions | Furey | 2012-2025 |
| G_SM = Stab_{Spin(9)}(J) | Krasnov | 2019 |
| J₃(O) + triality → SM + 3 generations | Boyle | 2020 |
| n-qubit entanglement → Hopf → SM gauge group | Szangolies | 2025 |
| J₃(O) → parameter-free fermion mass ratios | Singh | 2025 |
| Z₂⁵-graded superalgebra → SM | Furey | 2025 |
| Nonassociative spectral geometry | Farnsworth | 2025 |
| QRE asymmetry = D(ρ‖G(ρ)) | Marvian-Spekkens | 2014 |

### Szangolies (2025) 的原始表

```
┌─────────┬──────────┬────────────┬───────────────────┐
│ Qubits  │ Algebra  │ Hopf       │ Gauge Group       │
├─────────┼──────────┼────────────┼───────────────────┤
│ 1       │ C        │ S³→S²      │ U(1)              │
│ 2       │ H        │ S⁷→S⁴     │ SU(2)×U(1)/Z₂    │
│ 3       │ O        │ S¹⁵→S⁸   │ SU(3)×SU(2)×U(1) │
└─────────┴──────────┴────────────┴───────────────────┘
```

Szangolies 沒有：Level 0, Level 4, Σ 詮釋, Petz/τ 連結。

## 我們的貢獻

### 我們加的新元素

1. **Level 0 (R → 時間箭頭)**：
   沒有人把 R（實數/經典）→ Σ ≥ 0 → 時間箭頭放在 division algebra 階梯的起點。
   Szangolies 從 1-qubit 開始。我們從 0-qubit（經典）開始。
   這對應 Paper 1 的 τ = 1-F 框架。

2. **Level 1 的 Σ 連結（Paper 6 Theorem 1）**：
   別人知道 QM 需要 complex。我們證明了 Σ = D(ρ‖σ) 的 C*-algebra 結構 +
   Petz tensor product + Poincaré + 實驗（Nature 2022）五條獨立路線收斂到 U(1)。

3. **Level 2 的 Σ_weak 定義**：
   Marvian-Spekkens 已知 QRE asymmetry = D(ρ‖G(ρ))。
   我們將此應用到 division algebra Level 2，定義：
   Σ_weak = D(G_{SU(2)}(ρ) ‖ I/4)
   並數值驗證到 10⁻¹⁶ 精度（2026-03-19）。

4. **Level 4 (J₃(O) → 3 代)**：
   Boyle 做了 J₃(O) → SM，Singh 做了質量比。
   把它放在 Level 4 作為 Level 0-3 的自然延伸是我們的組織方式。

5. **Σ 詮釋欄（整個框架的新敘事）**：
   「Σ = D(ρ_spacetime ‖ ρ_matter) 在 n-qubit 系統上的自然展開
   重現了 division algebra 階梯」— 這個統一敘事是我們的。

6. **連接到重力/暗物質/時間箭頭**：
   Furey、Szangolies、Boyle 都不做重力或暗物質。
   Connes 做重力但不做暗物質。
   把 division algebra 階梯同時連到 Papers 1-6 的全部結果
   （τ, Petz, 指數度規, Khronon DM, DBI-Σ）是我們的。

### 我們的完整表

```
┌───────┬────────┬─────────┬────────────┬──────────────────┬────────┐
│ Level │ Qubits │ Algebra │ Force      │ Σ interpretation │ Status │
├───────┼────────┼─────────┼────────────┼──────────────────┼────────┤
│ 0     │ 0      │ R       │ Time arrow │ Σ ≥ 0 (DPI)      │ ✅     │
│ 1     │ 1      │ C       │ EM U(1)    │ Σ needs complex   │ ✅     │
│ 2     │ 2      │ H       │ Weak SU(2) │ Σ chain rule      │ ✅     │
│ 3     │ 3      │ O       │ Strong SU(3)│ (to do)          │ 🔨     │
│ 4     │ —      │ J₃(O)   │ 3 gen + mass│ Σ extremal       │ 💡     │
└───────┴────────┴─────────┴────────────┴──────────────────┴────────┘

Level 0: 新（我們加的）
Level 1: Σ 連結新（Paper 6 Theorem 1）
Level 2: Σ_weak 定義新（2026-03-19 驗證）
Level 3: 數學已知（Szangolies），Σ 實現待做
Level 4: 組織方式新（Boyle/Singh 的結果放在階梯中）
```

## 論文中的正確寫法

```
"Building on Szangolies' qubit-entanglement characterization of the
SM gauge group [Szangolies 2025], we extend the correspondence in
two directions: downward to Level 0, where the real number field R
corresponds to the positivity of entropy production (Σ ≥ 0, Paper I);
and upward to Level 4, where the exceptional Jordan algebra J₃(O)
[Boyle 2020, Singh 2025] corresponds to the extremal structure of Σ.

The division algebra ladder was independently connected to the SM by
Günaydin-Gürsey (1973, octonions → color), Dixon (1994, R⊗C⊗H⊗O),
Furey (2012-2025, minimal ideals → fermions), Krasnov (2019, SO(9)
characterization), and Connes (1996, spectral triple A_F = C⊕H⊕M₃(C)).

Our contribution is the information-theoretic interpretation: at each
level, the quantum relative entropy D(ρ‖σ) evaluated on n-qubit
systems naturally exhibits the gauge structure of the corresponding
division algebra. Level 1 (C → U(1)) follows from the C*-algebraic
requirements of QRE (Theorem 1). Level 2 (H → SU(2)) follows from
the exact chain rule D(ρ‖I/d) = D(ρ‖G(ρ)) + D(G(ρ)‖I/d) under
the SU(2) twirl [cf. Marvian-Spekkens 2014]."
```

## 記錄日期
2026-03-19, created during research session on matter sector strategy.
