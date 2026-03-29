# μ 的突破：從失敗到三條收斂線索
## 2026-03-19

## 背景
μ₀ = H₀/c 被 CMB 排除 >55σ。BS2025 用 μ⁻¹ = 22.3 Mpc（常數）通過 CMB。

## 三條獨立線索

### 線索 1：μ·c² = a₀·(1+z_dec)
- μ⁻¹ = c²/[a₀·(1+z_dec)] = 22.25 Mpc（吻合 BS2025 到 0.2%）
- 物理：Khronon 質量 = MOND 加速度在光子退耦時的值
- 文獻中從未出現過

### 線索 2：μ = 2π(1+Ω_b)/r_d
- μ⁻¹ = r_d/[2π(1+Ω_b)] = 22.31 Mpc（吻合 0.04%）
- r_d = 147.09 Mpc (Planck sound horizon at baryon drag)
- 物理：Khronon 質量 = sound horizon 的因果約束 + 重子修正
- 200 = 2π × (Hubble radius)/(Sound horizon)

### 線索 3：Running μ(z=49) = H(49)/c = BS2025 的值
- BS Model 1 (223 kpc) ≈ μ(z=1077)（14% 內）
- BS Model 2 (22.3 Mpc) = μ(z=49)（精確）
- 兩個「常數模型」是同一個 running 在不同 epoch 的快照

## 翻轉的結論

μ₀ = H₀/c 作為 z=0 的值是對的。
作為所有 epoch 的常數是錯的。
Running μ(a) = H(a)/c 可能同時解決 CMB 和 MOND：
  - 早期 (z~1100): μ 大 → w~0 → CDM ✅
  - 今天 (z~0): μ 小 → MOND ✅

## 額外巧合
Ghost condensation scale M = √(μ·M_Pl) ≈ 59 meV ≈ 最小中微子質量和

## CLASS 驗證結果（2026-03-19 完成）

### Running μ(a) = H(a)/c：❌ 排除
- chi²/dof = 5524（λ_D=1, khronon cs2）
- 原因：Z = δ₀/Ω_m ≈ 1.08 = 常數 in matter era → w ≈ 0.2 at z=1100

### 常數 μ = H₀/c + DBI：❌ 排除
- chi²/dof = 178.5
- 早期：w ~ 10⁻¹⁵（CDM-like ✅ — DBI works!）
- 晚期：w ~ 0.13（太大 ❌）
- 原因：Ω_DM=0.265 + μ=H₀/c → δ₀≈0.34 → w₀≈0.13

### 根本問題
要 Ω_DM=0.265 且 w₀≈0：需要 μ >> H₀/c（使 δ₀ << 1）
μ = H₀/c 強制 δ₀ ~ O(1) → w₀ ~ O(0.1)。代數事實，無法逃脫。

### μ₀ = H₀/c 最終判決：EXCLUDED
在所有組合下（quadratic, DBI+running, DBI+constant）都被 CMB 排除。
BS2025 的 μ⁻¹ = 22.3 Mpc 通過。

### 三條線索的狀態
仍然有趣（μ·c² = a₀(1+z_dec), μ = 2π(1+Ω_b)/r_d 都吻合 22.3 Mpc）
但它們不等於 H₀/c。μ 的理論推導仍然是 open problem。

### CLASS 修改的貢獻
gdm_class_public 增加了 khronon_dbi_cs2 模式（k-dependent c_s²）
這是一個對社群有用的工具，不管我們的理論成敗。

## 歸屬
- μ·c² = a₀(1+z_dec)：新發現（2026-03-19）
- μ = 2π(1+Ω_b)/r_d：新發現（2026-03-19）
- Running μ = BS models 的快照：新觀察（2026-03-19）
- BS2025 μ⁻¹ = 22.3 Mpc：Blanchet & Skordis
- Running μ(k) = k 假設：我們的 Paper 4
