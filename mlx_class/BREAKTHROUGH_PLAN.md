# mlx_class 數學突破計劃

## 目標：<0.1% RMS vs CLASS，不抄 CLASS 答案

## 當前狀態（2026-03-22）
- TT RMS: 6.57%
- 根因: gauge transformation cancellation（sync → Newtonian）
- 速度: ODE 1.07s (Magnus GPU) / 253s (scipy)

---

## 方向 1：無 Gauge Transform LOS（直接在 sync gauge 算 C_l）

### 原理
C_l 是 gauge-invariant。不需要轉到 Newtonian gauge。
直接用 sync gauge 變量寫 LOS 積分：

```
Δ_l(k) = ∫dτ [g(τ)·δ_γ/4·j_l + g(τ)·θ_b/k·j_l' + e^{-κ}·(η' + h'/6·P₂(cosθ))·j_l]
```

δ_γ = F_γ,0（sync gauge，精確）
θ_b = sync gauge baryon velocity（精確）
η', h' = diagnosed from constraints（精確）

**零 gauge transformation → 零 cancellation → 精度由 ODE solver 決定**

### 風險：低
sync gauge 的 LOS 公式在文獻中存在（Ma & Bertschinger 1995），只是沒人這樣 implement。

### 檔案：`solver_sync_direct.py`

---

## 方向 2：Chebyshev 頻譜法（一次矩陣求解取代 600 步 ODE）

### 原理
y'(τ) = A(τ)·y(τ) 是線性 ODE。展開 y(τ) = Σ aₙ Tₙ(τ)：

```
Σ aₙ Tₙ'(τ) = A(τ) · Σ aₙ Tₙ(τ)
```

離散化 → (D - A_cheb) · a = b （D = Chebyshev 微分矩陣）

一次矩陣求解：a = (D - A_cheb)⁻¹ · b

### 優勢
- 指數收斂：30 個 Chebyshev 點 → machine precision
- 天然 GPU（大矩陣求解）
- 所有 k-modes 可以 batch

### 風險：中
- A(τ) 隨 τ 變化（不是常數矩陣），需要 Chebyshev interpolation
- Thomson scattering 的 stiffness 可能需要特殊處理
- 系統大小：(N_cheb × n_var) × (N_cheb × n_var) ≈ (30×50) × (30×50) = 1500×1500

### 檔案：`solver_chebyshev.py`

---

## 方向 3：Autodiff HMC（可微分 Boltzmann solver）

### 原理
MLX 支持 automatic differentiation。如果整個 C_l pipeline 可微分：

```python
def Cl_from_params(params):
    bg = solve_background(params)
    y = solve_boltzmann(bg)
    Cl = compute_Cl(y)
    return Cl

grad_Cl = mx.grad(Cl_from_params)(params)  # 一次反向傳播
```

用 Hamiltonian Monte Carlo 取代 Metropolis-Hastings：
- 6 維參數空間：HMC ~200 步 vs MH ~5000 步
- 每步用 gradient 引導（不是隨機走）

### 風險：中
- scipy.integrate.solve_ivp 不可微 → 需要用 MLX 的 ODE solver
- Magnus solver 是可微的（矩陣運算全在 MLX 上）

### 檔案：`solver_autodiff.py`

---

## 方向 4：Neural Emulator（自訓練，非抄 CLASS）

### 原理
用 mlx_class 自己產生 (params, C_l) 訓練數據，訓練小型 NN：

```
訓練：1000 組 × 1s/組 = 15 分鐘
推理：<1ms/cosmology
精度上限 = mlx_class 的精度
```

### 風險：低
- 已有成熟方法（COSMOPOWER 等）
- 但精度受限於 solver 精度（如果 solver 6.57%，emulator 也是 6.57%）
- 需要先提高 solver 精度才有意義

### 檔案：`emulator.py`

---

## 方向 5：Riccati Reduction（用 4×4 矩陣取代 25 個 multipole）

### 原理
光子 Boltzmann hierarchy 的 Green's function 可以用 matrix Riccati equation 表示：

```
dΣ/dτ = A + BΣ + ΣC + ΣDΣ
```

Σ 是 4×4 矩陣（monopole, dipole, quadrupole, + metric coupling）。
等效於整個 25-multipole hierarchy。

### 風險：高
- Riccati 在 CMB 上下文中未被使用過
- 非線性（ΣDΣ 項）需要特殊處理
- 理論推導複雜

### 檔案：`solver_riccati.py`

---

## 執行計劃

### 第一輪（平行 3 路，今天開始）

```
Agent A: 方向 1（無 gauge transform LOS）  ← 風險最低，直接修精度根因
Agent B: 方向 2（Chebyshev 頻譜法）        ← 最大創新潛力
Agent C: 方向 5（Riccati reduction）       ← 最大理論突破
```

### 比較指標

| 指標 | 目標 |
|------|------|
| TT RMS vs CLASS | <1%（第一輪），<0.1%（最終） |
| D_l(220)/CLASS | 0.99-1.01 |
| D_l(537)/CLASS | >0.90（第二峰） |
| ODE 速度 | <2s (500 k-modes) |
| 總速度 | <5s |

### 第二輪（根據第一輪結果）

- 最好的方向 → 精調到 <0.1%
- 加上方向 3（autodiff）做 MCMC
- 加上方向 4（emulator）做即時推理

### 成功標準

| 等級 | RMS | 意義 |
|------|-----|------|
| Bronze | <3% | 第二峰修好，所有峰可見 |
| Silver | <1% | 可以做 Khronon vs ΛCDM 初步比較 |
| **Gold** | **<0.1%** | **可以 fit Planck 數據** |
| Platinum | <0.01% | 和 CLASS 完全等價 |
