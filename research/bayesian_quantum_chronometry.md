# Bayesian Quantum Chronometry (BQC): A Unified Framework

**Author: Sheng-Kai Huang**
**Date: 2026-03-16**

---

## 0. Motivation

Paper 1 ("The Arrow of Time from Petz Recovery") establishes three mathematical discoveries: a Composition Theorem (triangle inequality for $\sqrt{\tau}$), a Saturation Theorem (exact characterization of when $F^2 = \exp(-\Delta D)$), and the $\tau - \tilde{\tau}$ gap (the cost of Bayesian consistency). These were presented as independent results. This document argues they are three projections of a single geometric structure—**Bayesian Quantum Chronometry (BQC)**—that governs how time, noise, and irreversibility interrelate in quantum processes.

核心直覺：三個定理分別回答了三個關於量子時間箭頭的基本問題：
- **Composition**: 時間如何累積？（沿著量子電路，不可逆性如何逐級疊加）
- **Saturation**: 時間何時消失？（什麼精確條件讓 $\tau = 0$，即可逆與不可逆的分界線）
- **Gap**: 量子性的代價是什麼？（Bayesian consistency 要求在量子世界中額外付出多少代價）

BQC 的論點：這三者共同構成一個關於「量子時間」的完整幾何理論，其中 Bures 空間上的距離就是不可逆性的度量。

---

## 1. The Three Pillars

### 1.1 Composition Theorem: The Accumulation of the Arrow

For sequential quantum channels $\mathcal{N}_1, \mathcal{N}_2$ with reference state $\sigma$:

$$\sqrt{\tau_{12}} \leq \sqrt{\tau_1} + \sqrt{\tau_2^{\text{eff}}}$$

where $\tau = 1 - F(\rho, \mathcal{R}_{\sigma,\mathcal{N}}(\mathcal{N}(\rho)))$ is the retrodiction infidelity under the standard Petz map, and crucially:

$$\tau_2^{\text{eff}} = 1 - F\big(\mathcal{N}_1(\rho),\; \mathcal{R}_{\mathcal{N}_1(\sigma), \mathcal{N}_2}(\mathcal{N}_2(\mathcal{N}_1(\rho)))\big)$$

uses the **updated** reference $\mathcal{N}_1(\sigma)$, not the original $\sigma$.

**Physical content**: The total irreversibility of a multi-stage process decomposes sub-additively in the Bures metric, but each stage's contribution depends on what happened before. This is not mere additivity—it is **contextual accumulation** with Bayesian updating.

**Why $\sqrt{\tau}$ and not $\tau$?** Because $\sqrt{1-F}$ is the Bures distance (up to a factor of $\sqrt{2}$), and Bures distance satisfies the triangle inequality while fidelity does not compose additively. The square root transforms multiplicative fidelity into additive geometry.

### 1.2 Saturation Theorem: The Boundary of Reversibility

For $\sigma = I/d$ (maximally mixed reference) and pure input $\rho = |\psi\rangle\langle\psi|$, with $\omega = \mathcal{N}(\rho)$ and $\nu = \mathcal{N}(\sigma)$:

$$F^2 = \exp(-\Delta D) \quad \Longleftrightarrow \quad [\omega, \nu] = 0 \;\text{ AND }\; \omega_i/\nu_i = \text{constant}$$

**Physical content**: The bound $F^2 \geq \exp(-\Delta D)$ is tight if and only if the channel acts as a **quantum sufficient statistic**—it preserves the distinguishability ratio uniformly. This is the exact boundary between the realm where the Petz map is "perfectly matched" to the noise and where it is not.

兩個條件缺一不可：
- $[\omega, \nu] = 0$：通道輸出可以在同一個基底下對角化（經典化條件）
- $\omega_i/\nu_i = c$：似然比（likelihood ratio）在支撐上是常數（充分統計量條件）

### 1.3 The $\tau - \tilde{\tau}$ Gap: The Quantum Tax on Time Reversal

$$\tau - \tilde{\tau} \geq 0$$

where $\tau$ uses the standard (Bayesian-consistent) Petz map $\mathcal{R}_{\sigma,\mathcal{N}}$ and $\tilde{\tau}$ uses the rotated Petz map $\tilde{\mathcal{R}}_{\sigma,\mathcal{N}}$ (which is near-optimal but not Bayesian-consistent).

**Vanishing conditions**: The gap is zero when:
- $[\rho, \sigma] = 0$ (classical regime)
- $\Delta D = 0$ (exact recovery, both maps achieve $F = 1$)

**Physical content**: In the quantum regime, Bayesian consistency costs extra fidelity. You can be Bayesian (standard Petz, composes nicely, has physical meaning as retrodiction) or you can be optimal (rotated Petz, tightest bound), but not both. The gap $\tau - \tilde{\tau}$ is the **irreducible price of having a coherent theory of time reversal** in the quantum setting.

---

## 2. The Unified Story

### 2.1 BQC as a Geometric Theory of Quantum Time

在 BQC 框架中，量子過程的時間結構可以用一張圖表達：

```
State space (Bures geometry)
         ρ
        / \
  τ₁  /   \ τ₁₂
      /     \
  ρ₁ ------- ρ₁₂
     τ₂^eff

Triangle inequality: √τ₁₂ ≤ √τ₁ + √τ₂^eff
Saturation vertex:   τᵢ = 0 iff stage i is a quantum sufficient statistic
Gap at each vertex:  τᵢ ≥ τ̃ᵢ (Bayesian tax)
```

The three results then have a single geometric reading:

| Result | Geometric Role |
|--------|---------------|
| Composition | The triangle inequality in Bures space |
| Saturation | The criterion for a vertex to sit at zero (on the "reversibility surface") |
| Gap | The minimal thickness of the triangle at each vertex |

### 2.2 How They Interlock

**Question 1: What happens when you compose channels that individually satisfy saturation?**

If stage $i$ satisfies the saturation conditions with updated reference $\sigma_i = \mathcal{N}_{i-1} \circ \cdots \circ \mathcal{N}_1(\sigma)$, then $\tau_i^{\text{eff}} = 0$ at each stage. By the composition theorem:

$$\sqrt{\tau_{\text{total}}} \leq \sum_i \sqrt{\tau_i^{\text{eff}}} = 0$$

Therefore $\tau_{\text{total}} = 0$: **saturation is preserved under composition with Bayesian updating.** This is non-trivial. It says that if every stage independently acts as a quantum sufficient statistic for the Bayesian-updated reference, the entire process is perfectly reversible.

But here is the subtlety: the saturation conditions $[\omega_i, \nu_i] = 0$ must hold at each stage with the **updated** reference, not the original one. Even if $\mathcal{N}_1$ saturates for $\sigma$, it does not follow that $\mathcal{N}_2$ saturates for $\mathcal{N}_1(\sigma)$. The composition forces the conditions to propagate forward in a Bayesian chain.

**Question 2: Does the gap accumulate?**

Consider the gap at each stage: $\delta_i = \tau_i^{\text{eff}} - \tilde{\tau}_i^{\text{eff}} \geq 0$. We can derive:

$$\sqrt{\tau_{\text{total}}} \leq \sum_i \sqrt{\tau_i^{\text{eff}}} = \sum_i \sqrt{\tilde{\tau}_i^{\text{eff}} + \delta_i}$$

Since $\sqrt{a+b} \leq \sqrt{a} + \sqrt{b}$ for $a,b \geq 0$:

$$\sqrt{\tau_{\text{total}}} \leq \sum_i \sqrt{\tilde{\tau}_i^{\text{eff}}} + \sum_i \sqrt{\delta_i}$$

This yields the **Gap Composition Theorem**:

> **Theorem (Gap Composition).** For an $n$-stage quantum process,
> $$\sqrt{\tau_{\text{total}}} \leq \sqrt{\tilde{\tau}_{\text{total}}^{\text{bound}}} + \Delta_{\text{Bayes}}$$
> where $\sqrt{\tilde{\tau}_{\text{total}}^{\text{bound}}} = \sum_i \sqrt{\tilde{\tau}_i^{\text{eff}}}$ is the accumulated irreversibility under optimal (non-Bayesian) recovery, and
> $$\Delta_{\text{Bayes}} = \sum_i \sqrt{\delta_i}$$
> is the **accumulated Bayesian tax**—the total price paid for maintaining a coherent retrodiction theory across all stages.

重要觀察：$\Delta_{\text{Bayes}}$ 可以增長為 $O(\sqrt{n})$（如果每步的 gap 是常數），但 $\sqrt{\tilde{\tau}_{\text{total}}^{\text{bound}}}$ 也是 $O(\sqrt{n})$ 量級。因此 Bayesian tax 佔總不可逆性的比例在 $n \to \infty$ 時不會消失——它是一個 **persistent quantum overhead**。

**Question 3: What is the relationship between saturation and Bayesian updating?**

Saturation requires $[\omega_i, \nu_i] = 0$—the channel outputs must commute. But the Bayesian update $\sigma \mapsto \mathcal{N}(\sigma)$ generically produces a state that does not commute with $\mathcal{N}(\rho)$ for subsequent channels. Therefore:

> **Observation (Bayesian Decoherence).** The Bayesian updating process $\sigma_i = \mathcal{N}_i(\sigma_{i-1})$ acts as a form of decoherence on the reference: each channel pushes the updated reference toward the channel's fixed-point set. Saturation at stage $i+1$ requires that this "Bayesian decoherence" aligns with the next channel's action.

This gives a geometric picture: the saturation surface in state space is a kind of "classical manifold" where $[\omega, \nu] = 0$. Bayesian updating traces a trajectory through state space, and reversibility is maintained exactly when this trajectory stays on the saturation surface at every step.

---

## 3. The Central Theorem of BQC

> **Theorem (Bayesian Quantum Chronometry).** Let $\{\mathcal{N}_i\}_{i=1}^n$ be a sequence of quantum channels with initial reference $\sigma_0 = \sigma$, and define the Bayesian-updated references $\sigma_i = \mathcal{N}_i(\sigma_{i-1})$. For any input state $\rho$, the total retrodiction infidelity satisfies:
>
> $$\boxed{\sqrt{\tau_{\text{total}}} \leq \sum_{i=1}^n \sqrt{\tau_i^{\text{eff}}}}$$
>
> where each $\tau_i^{\text{eff}} = 1 - F(\rho_i, \mathcal{R}_{\sigma_{i-1}, \mathcal{N}_i}(\mathcal{N}_i(\rho_i)))$ with $\rho_i = \mathcal{N}_{i-1} \circ \cdots \circ \mathcal{N}_1(\rho)$, subject to:
>
> **(Bayesian Tax):** $\tau_i^{\text{eff}} \geq \tilde{\tau}_i^{\text{eff}}$ at each stage, with equality iff $[\rho_i, \sigma_{i-1}] = 0$.
>
> **(Saturation Criterion):** $\tau_i^{\text{eff}} = 0$ if and only if
> (i) $[\mathcal{N}_i(\rho_i), \mathcal{N}_i(\sigma_{i-1})] = 0$ and
> (ii) $\mathcal{N}_i(\rho_i)_j / \mathcal{N}_i(\sigma_{i-1})_j = c_i$ for all $j$ in the support (for $\sigma_{i-1} = I/d_i$; the general criterion requires quantum sufficiency of $\mathcal{N}_i$ for the pair $(\rho_i, \sigma_{i-1})$).
>
> **(Irreversibility Decomposition):**
> $$\sqrt{\tau_{\text{total}}} \leq \underbrace{\sum_i \sqrt{\tilde{\tau}_i^{\text{eff}}}}_{\text{optimal cost}} + \underbrace{\sum_i \sqrt{\delta_i}}_{\text{Bayesian tax}}$$
>
> The total arrow of time decomposes into an optimal (non-Bayesian) component and a Bayesian-consistency tax.

---

## 4. The Deep Paradox: Computation Requires an Arrow of Time

### 4.1 The Saturation-Computation Incompatibility

Consider a quantum computer executing a circuit of $n$ gates $\{U_i\}$, each subject to noise $\mathcal{E}_i$. The effective channel at each stage is $\mathcal{N}_i = \mathcal{E}_i \circ \mathcal{U}_i$.

For perfect error correction (i.e., $\tau_i^{\text{eff}} = 0$ at every gate), the saturation conditions require:

$$[\mathcal{N}_i(\rho_i), \mathcal{N}_i(\sigma_{i-1})] = 0 \quad \forall i$$

But useful quantum computation **requires** that the gates create superpositions and entanglement—precisely the situations where $\mathcal{N}_i(\rho_i)$ and $\mathcal{N}_i(\sigma_{i-1})$ generically do not commute (unless the noise is specifically aligned with the gate).

> **Paradox (Computation-Reversibility Tension).** A quantum computation that satisfies saturation at every gate requires $[\omega_i, \nu_i] = 0$ at every stage—meaning the channel outputs are effectively classical at each step. But a quantum computer that only produces classical outputs at each stage cannot leverage quantum parallelism.

### 4.2 Quantifying the Paradox

We can make this precise. Define the **noncommutativity measure** at stage $i$:

$$\mathcal{C}_i = \|[\mathcal{N}_i(\rho_i), \mathcal{N}_i(\sigma_{i-1})]\|_1$$

When $\mathcal{C}_i > 0$, the saturation theorem guarantees $F_i^2 > \exp(-\Delta D_i)$, meaning the Petz recovery is strictly suboptimal. More precisely, from the necessity proof (the logarithmic-mean vs geometric-mean gap):

$$\tau_i^{\text{eff}} \geq \tau_i^{\text{sat}} + f(\mathcal{C}_i)$$

where $\tau_i^{\text{sat}}$ is the irreversibility that would remain even with saturation, and $f(\mathcal{C}_i) > 0$ is a strictly positive function of the noncommutativity. For small $\mathcal{C}_i$, from the perturbation analysis in the saturation proof:

$$f(\mathcal{C}_i) \sim \sum_{j \neq k} |\omega_{jk}|^2 \left[\frac{1}{\sqrt{t_j t_k}} - \frac{1}{L(t_j, t_k)}\right]$$

where $t_j = \nu_j$ are the eigenvalues of the noise-on-reference output and $L(a,b)$ is the logarithmic mean.

This yields the **Computation Lower Bound**:

> $$\tau_{\text{computation}} \geq \sum_i f(\mathcal{C}_i)$$
>
> The total irreversibility of a quantum computation is bounded below by the total noncommutativity generated at each gate.

**Interpretation**: 時間的箭頭不是量子計算的「bug」——它是量子計算的**代價**。量子平行性（superposition、entanglement）必然產生 $\mathcal{C}_i > 0$，而 $\mathcal{C}_i > 0$ 必然產生 $\tau_i > 0$。一台完全可逆的量子電腦必然是一台無用的（經典的）電腦。

### 4.3 The Error Correction Tradeoff

This paradox resolves into a tradeoff. Quantum error correction does not eliminate the arrow of time; it **redistributes** it:

- **Without QEC**: $\tau_i^{\text{eff}}$ accumulates freely, $\sqrt{\tau_{\text{total}}} \leq \sum_i \sqrt{\tau_i}$.
- **With QEC**: The decoder $\mathcal{D}_i$ (approximating the Petz map) reduces $\tau_i^{\text{eff}}$ but introduces its own cost: the overhead of syndrome extraction, ancilla preparation, etc.
- **Fault-tolerance threshold**: The threshold theorem can be recast as: there exists a critical noise rate $p_c$ below which the per-gate effective $\tau_i^{\text{eff}}$ can be suppressed to $O(p^{d/2})$ (for code distance $d$), so that $\sqrt{\tau_{\text{total}}} = O(n \cdot p^{d/4})$, which vanishes in the limit $d \to \infty$ for $p < p_c$.

---

## 5. Physical Analogies

### 5.1 Classical EDA: Static Timing Analysis

In classical chip design, **Static Timing Analysis (STA)** computes the critical path delay as:

$$T_{\text{total}} = \sum_i d_i + \sum_j s_j$$

where $d_i$ are gate delays and $s_j$ are setup/hold margins. The delays are additive along paths.

BQC provides the quantum analog:

| Classical STA | BQC |
|--------------|-----|
| Gate delay $d_i$ | $\sqrt{\tau_i^{\text{eff}}}$ (retrodiction cost) |
| Setup/hold margin $s_j$ | $\sqrt{\delta_j}$ (Bayesian tax) |
| Critical path = $\sum d_i$ | $\sqrt{\tau_{\text{total}}} \leq \sum \sqrt{\tau_i^{\text{eff}}}$ |
| Delay is additive | $\sqrt{\tau}$ is sub-additive (triangle ineq.) |
| Context-independent | Context-dependent (Bayesian updating) |

The key difference: in classical STA, delays are fixed per gate. In BQC, each $\tau_i^{\text{eff}}$ depends on the **history** (through the updated reference $\sigma_{i-1}$). This is the quantum analog of state-dependent delay—a gate's "temporal cost" depends on the information landscape inherited from all previous gates.

**Practical implication for quantum circuit optimization**: Just as classical place-and-route minimizes critical path delay, quantum compilation should minimize $\sum_i \sqrt{\tau_i^{\text{eff}}}$ along the critical path. This provides a principled cost function for gate ordering and noise-aware compilation.

### 5.2 Thermodynamics: $\Sigma$ vs $\tau$

Entropy production $\Sigma$ and retrodiction infidelity $\tau$ are related by the master inequality chain:

$$-\log F^2 \leq I(A;E|B) \leq \Sigma \leq \Delta D$$

Both quantify irreversibility, but they have different **compositional structures**:

| Property | $\Sigma$ (entropy production) | $\sqrt{\tau}$ (Bures distance) |
|----------|------------------------------|-------------------------------|
| Composition | $\Sigma_{12} = \Sigma_1 + \Sigma_2$ (additive for Markov) | $\sqrt{\tau_{12}} \leq \sqrt{\tau_1} + \sqrt{\tau_2^{\text{eff}}}$ (sub-additive) |
| Reference dependence | Fixed (Gibbs state, equilibrium) | Updated at each stage (Bayesian) |
| Uniqueness | Many entropy production definitions | Unique (Petz is the only Bayesian retrodiction) |
| Sign | $\Sigma \geq 0$ always (Markov); can be negative (non-Markov) | $\tau \geq 0$ always; $\tau_{\text{signed}}$ can be negative |

The critical insight: $\Sigma$ is additive because it does not care about the reference updating. $\sqrt{\tau}$ is sub-additive because it incorporates Bayesian learning. The gap between additivity and sub-additivity is precisely the information that the system "learns" about itself through the process.

### 5.3 General Relativity: Proper Time Accumulation

In GR, proper time along a worldline is:

$$\Delta\tau_{\text{GR}} = \int \sqrt{-g_{\mu\nu} dx^\mu dx^\nu}$$

This is additive along the worldline (geodesic or not), with the local metric $g_{\mu\nu}$ providing the "tick rate" at each event.

In BQC, the total retrodiction cost along a quantum circuit is:

$$\sqrt{\tau_{\text{total}}} \leq \sum_i \sqrt{\tau_i^{\text{eff}}}$$

The parallel is:

| GR | BQC |
|----|-----|
| Worldline through spacetime | Quantum circuit (sequence of channels) |
| $\sqrt{-g_{\mu\nu} dx^\mu dx^\nu}$ | $\sqrt{\tau_i^{\text{eff}}}$ |
| Proper time is additive | $\sqrt{\tau}$ is sub-additive |
| Metric is local | $\tau_i^{\text{eff}}$ depends on history (non-local in time) |
| Twin paradox: different paths, different ages | Circuit paradox: different gate orderings, different $\tau_{\text{total}}$ |

BQC 的「度規」比 GR 更複雜：它不僅依賴局部（當前通道），還依賴整個歷史（通過 Bayesian updating）。但這正是量子力學的特性——歷史是不能被遺忘的，除非你付出 $\tau > 0$ 的代價。

有一個引人注目的對應：在 GR 中，自由落體（測地線）最大化固有時。在 BQC 中，理想糾錯（saturation at every stage）最小化 $\tau_{\text{total}}$。兩者都是「最優路徑」，但方向相反：GR 最大化時間的流逝，BQC 最小化時間箭頭的累積。

This reversal has a deep reason: in GR, proper time measures how much time has *elapsed* (more is "more natural"); in BQC, $\tau$ measures how much time-reversal symmetry has been *broken* (less means "more reversible"). They are complementary aspects of the same phenomenon.

---

## 6. Connection to Paper 1's Core Philosophy

Paper 1 的核心論點是：Petz recovery map 同時是 Bayesian retrodiction functor（唯一的）、近最優 QEC decoder、以及量子 Crooks 定理中的逆向通道。BQC 框架將這個「四域等價鏈」

$$\text{Eraser} \leftrightarrow \text{Retrodiction} \leftrightarrow \text{QEC} \leftrightarrow \text{Thermo}$$

從靜態的等價提升為動態的幾何理論：

1. **等價鏈告訴我們 $\tau = 0$ 時發生什麼** → BQC 告訴我們 $\tau > 0$ 時如何累積
2. **等價鏈是一個快照（單一通道）** → BQC 處理序列（多通道組合）
3. **等價鏈是四個域的字典** → BQC 是在這四個域上同時運作的幾何

特別是，BQC 揭示了一個 Paper 1 未明說的結構：**時間箭頭的不可逆性有兩個來源**——物理雜訊（$\tilde{\tau}_i^{\text{eff}}$，可以用最優解碼器消除）和量子 Bayesian tax（$\delta_i$，只要 $[\rho, \sigma] \neq 0$ 就必須付出）。第二個來源是純量子的，不可消除的，而且正是使 Petz map 區別於其他 recovery map 的關鍵特徵。

---

## 7. Open Questions and Future Directions

1. **General $\sigma$ saturation**: The Saturation Theorem (Theorem 2 of Paper 1) is proved for $\sigma = I/d$. What are the saturation conditions for general full-rank $\sigma$? In the BQC context, this is crucial because the Bayesian-updated references $\sigma_i = \mathcal{N}_i(\sigma_{i-1})$ are generically not maximally mixed.

2. **Continuous-time BQC**: The composition theorem is stated for discrete channels. In the Lindblad setting $\dot{\rho} = \mathcal{L}(\rho)$, can we define $d\tau/dt$ as an instantaneous irreversibility rate? The connection to Spohn's theorem ($\Sigma = -\frac{d}{dt}D(\rho\|\rho_{\text{ss}}) \geq 0$ for Markovian dynamics) suggests yes.

3. **Non-Markovian extension**: When $\Sigma$ can temporarily decrease (non-Markovian backflow), can $\tau_{\text{total}}$ decrease between stages? The composition theorem does not forbid this: $\tau_{12}$ can be less than $\tau_1$ if $\tau_2^{\text{eff}}$ is very small (the second stage "heals" the first).

4. **Gap quantification**: What is the exact dependence of $\delta_i = \tau_i - \tilde{\tau}_i$ on the noncommutativity $\|[\rho, \sigma]\|$? A tight bound would quantify the exact price of Bayesian consistency.

5. **Holographic interpretation**: In AdS/CFT, the Petz map appears in entanglement wedge reconstruction (Chen-Penington-Salton 2020). Does the BQC composition theorem have a gravitational dual? Does $\sqrt{\tau_{\text{total}}} \leq \sum_i \sqrt{\tau_i^{\text{eff}}}$ correspond to a sub-additivity of some bulk geometric quantity?

6. **Experimental test**: In a multi-stage quantum circuit (e.g., on ion trap or superconducting qubits), measure $\tau_i^{\text{eff}}$ at each stage via Petz recovery (already demonstrated by Png 2025, Singh 2025). Test whether $\sqrt{\tau_{\text{total}}}$ saturates the composition bound—and whether the Bayesian tax $\sum_i \sqrt{\delta_i}$ can be independently measured by comparing standard and rotated Petz recovery at each stage.

---

## 8. Summary: The BQC Dictionary

| Concept | Mathematical Expression | Physical Meaning |
|---------|------------------------|-----------------|
| Arrow of time at stage $i$ | $\tau_i^{\text{eff}} = 1 - F_i$ | Retrodiction cost with updated prior |
| Total arrow of time | $\sqrt{\tau_{\text{total}}} \leq \sum_i \sqrt{\tau_i^{\text{eff}}}$ | Sub-additive accumulation (composition) |
| Reversibility boundary | $[\omega_i, \nu_i] = 0$ AND $\omega_i/\nu_i = c$ | Quantum sufficient statistic (saturation) |
| Bayesian tax | $\delta_i = \tau_i - \tilde{\tau}_i \geq 0$ | Price of consistent retrodiction (gap) |
| Computation cost | $\tau_{\text{comp}} \geq \sum_i f(\mathcal{C}_i)$ | Noncommutativity $\to$ irreversibility |
| Zero-entropy limit | $\tau_i^{\text{eff}} = 0 \;\forall i$ | Perfectly closed system, no time arrow |
| Classical limit | $\delta_i = 0$ | Bayesian tax vanishes, $\tau = \tilde{\tau}$ |

The central message of BQC: **時間的箭頭不是一個單一的數字，而是一個沿著量子過程累積的幾何量。它的累積服從三角不等式（composition），它的消失需要精確條件（saturation），而它的量子部分是不可避免的代價（gap）。這三者共同構成了量子時間的完整理論。**
