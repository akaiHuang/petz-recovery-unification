# 所需工具和論文清單
## 2026-03-20

---

## 一、論文狀態

### 已完成（可投）
| # | 論文 | 頁數 | 核心結果 | 檔案 |
|---|------|------|---------|------|
| 1 | Paper 1: τ = 1-F, Petz recovery | 6 | 等價鏈 | GitHub + Zenodo |
| 1b | 坍塌 = Petz recovery 失敗 | 6 | 本體論 + Penrose P3→P4 | paper1b_collapse_petz.tex |
| 2 | Σ_grav = 2 ln Q | 6 | Channel theorem | paper2_exponential_metric.tex |
| a₀ | MOND 加速度公式 | 4 | a₀ = 2πc²(1+Ωb)/[rd(1+zdec)], 0.3% | paper_a0_formula.tex |
| ν-MOND | Neutrino → 旋轉曲線 | 4 | 零參數鏈 + SPARC 驗證 + Crooks μ(x) | paper_neutrino_mond_chain.tex |
| 9 | Spectral Fisher = d=4 | 5 | g_QQ = 2/Q², (d-2)(d-3)=2 → d=4 | paper9_spectral_fisher.tex |
| 意識 | Consciousness as retrodiction | 8 | 五理論統一 + 23人 EEG 驗證 | paper_consciousness_sigma.tex |

### 需要寫的論文
| # | 論文 | 內容 | 依賴什麼 | 預計 |
|---|------|------|---------|------|
| 3 | 旋轉曲線 | J(Y) from Crooks + SPARC fit | SPARC 數據（已有） | 有材料，需整理 |
| 4 | μ₀ 修正版 | μ = 2π(1+Ωb)/rd（常數，不是 H₀/c） | a₀ 論文 | 需重寫（舊版 μ₀=H₀/c 被排除）|
| 5 | Complementary uncertainty | g(d) = (d+1) - √(2(d+1)) | 獨立 | 已有草稿 |
| 6 | EM 統一 | U(1), DBI-Σ, Modified Maxwell | 獨立 | 已有草稿 |
| 7 | SU(2) | Bianconi 1-form → Yang-Mills | Bianconi 1-form 計算（~7週）| 待 Bianconi 計算 |
| 8 | SU(3) | 併入 Paper 9 | Paper 9 | 不需獨立 |
| 10 | 費米子質量 | J₃(O) + Σ | NCG spectral triple | 長期 |
| 11 | Grand finale | 大統一 | 全部 | 最後 |
| JOSS | mlx_class 工具論文 | CMB solver | mlx_class 精度 0.1% | 精度達標後 |
| EEG | Petz bound in brain | 23人驗證 + retrodiction AI | 更多 subjects | 有材料 |

---

## 二、工具狀態

### mlx_class（CMB Boltzmann solver）

#### 已完成的功能
| 功能 | 模組 | 精度 |
|------|------|------|
| CMB TT | solver_unified.py | 峰位 0.9%（ISW template），振幅 -6% |
| CMB TE | polarization_full.py | 60-80% of CLASS |
| CMB EE | polarization_full.py | 86% of CLASS |
| CMB BB | tensor.py | scaling 1.9% |
| Lensing | lensing.py | C_l^φφ < 3% |
| P(k) linear | matter_pk.py | σ₈ 6% off |
| P(k) non-linear | halofit.py | Takahashi 2012 |
| Massive ν | massive_neutrino.py | 107 vars/k |
| Dark energy w(a) | dark_energy.py | CPL parametrization |
| Khronon DBI | khronon.py | 3-way comparison |
| Peebles recombination | recombination.py | visibility z=1090 |
| MCMC | mcmc.py | Metropolis-Hastings |
| GPU Bessel | bessel_gpu.py | 570x speedup |
| Fast solver | solver_fast.py | 0.5s (18x speedup) |
| 統一 solver | solver_unified.py | 全功能整合 |

#### 待修復
| 問題 | 目標 | 方法 |
|------|------|------|
| 精度 48.7% RMS vs CLASS | 0.1% | 按 accuracy_analysis.md 4 階段修 |
| 第一峰 l=303 不是 220 | l=220 | 正確 early ISW from ODE Phi_dot |
| 速度 ~10s（含 GPU Bessel） | <1s | TCA+近似 on GPU |
| EE 振幅 | 100% of CLASS | 需要完整 polarization hierarchy |

#### 待開發功能
| 功能 | 優先 | 時間 |
|------|------|------|
| 精度優化到 0.1% | ★★★ | 1-2 週 |
| 速度優化到 <1s | ★★★ | 1 週 |
| GPU Bessel 整合進統一 solver | ★★★ | 2 天 |
| Batch vmap（多參數同時跑） | ★★ | 3 天 |
| Autodiff（gradient MCMC） | ★★ | 2-4 週 |
| Non-flat Ω_k | ★ | 2 天 |
| HyRec 取代 Peebles | ★ | 1 週 |
| GitHub public push (mlx_class) | ★★ | 1 小時 |
| pip publish (PyPI) | ★★ | 1 天 |

### mlx-special（MLX 科學計算 ecosystem）

#### 已發布
| Package | PyPI | GitHub | Tests |
|---------|------|--------|-------|
| mlx-bessel | ✅ | ✅ akaiHuang/mlx-bessel | 7/7 |
| mlx-wigner | ✅ | ✅ akaiHuang/mlx-wigner | 44/44 |
| mlx-gamma | ✅ | ✅ akaiHuang/mlx-gamma | 20/20 |
| mlx-airy | ✅ | ✅ akaiHuang/mlx-airy | 16/16 |
| mlx-fisher | ✅ | ✅ akaiHuang/mlx-fisher | 30/30 |
| mlx-expm | ✅ | ✅ akaiHuang/mlx-expm | 25/25 |
| mlx-qre | ✅ | ✅ akaiHuang/mlx-qre | 35/35 |
| mlx-hyp2f1 | ✅ | ✅ akaiHuang/mlx-hyp2f1 | 21/21 |
| mlx-stft | GitHub only | ✅ akaiHuang/mlx-stft | 34/34 (PyPI 名字被佔) |

#### 網站
- https://akaihuang.github.io/mlx-special-site/

### EEG 實驗工具
| 工具 | 狀態 |
|------|------|
| Neural Σ estimator | ✅ 88.3% bound compliance |
| Retrodiction AI (Transformer) | ✅ 合成數據驗證 |
| Multi-subject analysis | ✅ 23 subjects, p < 10⁻¹⁰⁰ |
| Sleep-EDF 數據 | ✅ 25 subjects downloaded |
| Improved Σ estimator | ✅ 3 methods compared |

### 宇宙模擬
| 工具 | 狀態 | 可算什麼 |
|------|------|---------|
| anatropic 3D Euler+Poisson | ✅ 64³ | Phase 1-3（均勻→filament→星系團）|
| GR collapse code | ❌ 沒有 | Phase 4（黑洞形成）|
| QFT on curved spacetime | ❌ 沒有（全球少數人能算） | Phase 5（Hawking 蒸發）|
| 量子重力模擬 | ❌ 不存在 | Phase 6（資訊回流）|

---

## 三、理論缺口

### 致命缺口：0 個（全部消解）

### IMPORTANT 缺口
| 缺口 | 狀態 | 攻擊路線 |
|------|------|---------|
| (1+Ωb) 因子推導 | 3 候選全失敗，可能是 BS2025 μ 值不確定性 | 等 BS 更精確的 μ |
| Born rule | Petz → Gleason 路線（部分） | 長期開放問題 |
| SU(2) Level 3 | Bianconi 1-form（7 週計算） | Paper 7 |
| 3 代費米子 | J₃(O) 唯一性（沒人證明過） | 長期 |
| Vector field from Σ | BS2025 Khronon-Tensor 的 vector field 未從 Σ 推導 | Paper 4 修正版 |

### 已解決的缺口
| 缺口 | 解法 | 何時解決 |
|------|------|---------|
| Running μ(k)=k | 不需要（K 和 J 獨立） | 2026-03-20 |
| J(Y) → MOND | Crooks + Petz bound 1/2 | 2026-03-19 |
| μ 推導 | μ = 2π/r_d（5%）+ neutrino seesaw（1.7%） | 2026-03-19 |
| CLASS CMB 排除 | CMB = ΛCDM（ghost condensation c_s²=0） | 2026-03-19 |
| d=4 | g_QQ = (d-2)(d-3)/Q² = 2/Q² | 2026-03-19 |

---

## 四、Phase 4 黑洞模擬需要什麼

如果要做 Phase 4（黑洞形成）的真實模擬：

### 需要的軟體
| 軟體 | 用途 | 開源？ |
|------|------|--------|
| Einstein Toolkit | 數值廣義相對論 (BSSN formulation) | ✅ einsteintoolkit.org |
| SpEC / SpECTRE | 黑洞合併模擬 | ✅ (SXS collaboration) |
| GRChombo | AMR + GR | ✅ grchombo.org |
| Cactus | GR framework | ✅ cactuscode.org |

### 需要的計算資源
| 模擬 | 規模 | 需要 |
|------|------|------|
| 球對稱塌縮（最簡） | 1D, ~1小時 | 你的 Mac 就夠 |
| 軸對稱塌縮 | 2D, ~1天 | Mac 夠 |
| 完整 3D 塌縮 | 3D, ~1週 | 需要 cluster |

### 最實際的做法
1. 用 Einstein Toolkit 做 1D 球對稱塌縮
2. 在塌縮過程中計算 τ = 1-exp(-Σ/2)
3. 展示 τ 從 ~0 增長到 →1 的過程
4. 跟指數度規（Paper 2）對比

---

## 五、優先級總結

### 立即可做
1. mlx_class 精度 + 速度優化
2. lifecycle.html 改成 volume rendering（Phase 1-3 用真數據）
3. 整理現有論文材料

### 短期（1-2 週）
4. mlx_class push GitHub + PyPI
5. Bianconi 1-form SU(2) 計算開始

### 中期（1-2 月）
6. Paper 3 旋轉曲線完整論文
7. Paper 7 SU(2)
8. Einstein Toolkit 1D 黑洞塌縮

### 長期（3-12 月）
9. Paper 9 → SM 完整推導
10. 費米子質量
11. Grand finale
