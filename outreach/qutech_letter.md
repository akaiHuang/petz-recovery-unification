# Letter to QuTech / Quantum Inspire Team

## English Version

**To:** support@quantum-inspire.com / QuTech team
**Subject:** τ-chrono: open-source tool doubles effective circuit depth on Tuna-9 — collaboration proposal

Dear Quantum Inspire Team,

My name is Sheng-Kai Huang, an independent quantum physics researcher based in Taiwan. I have been using your Quantum Inspire platform and the Tuna-9 processor for my research on quantum noise characterization.

I developed an open-source tool called τ-chrono that applies Petz recovery map theory to predict circuit noise more accurately than the standard independent gate model. I validated it entirely on your Tuna-9 hardware, and the results are significant:

**NISQ noise prediction (validated on Tuna-9):**
- Depth 50: 46% more accurate noise prediction than independent model
- Bernstein-Vazirani: 60.4% improvement, 4 cases where the standard model was wrong
- H₂ VQE: Bayesian tracking keeps depth 4 viable (τ=0.49) where standard model says stop (τ=0.60)
- All depths from 2 to 50 show positive improvement

**QEC noise-informed decoding (validated on Tuna-9):**
- 3-qubit repetition code: 39.9% improvement with noise-informed correction
- Correctly identifies the noisiest qubit and down-weights it

In practical terms, this means Tuna-9 users can reliably run circuits roughly twice as deep as they currently do, without any hardware changes. The tool deploys in about 10 minutes using standard gate characterization data and requires zero additional training circuits.

Everything is open-source (MIT license):
- GitHub: https://github.com/akaiHuang/petz-recovery-unification
- QEC experiments: https://github.com/akaiHuang/qec-retrodiction-decoder

**I would like to propose three levels of collaboration:**

1. **Documentation reference** (minimal effort for your team): Add a mention of τ-chrono in the Quantum Inspire documentation or tutorials, so your users can benefit from improved circuit fidelity prediction. I am happy to write the tutorial myself.

2. **Platform integration** (medium effort): Integrate τ-chrono's noise prediction into the Quantum Inspire web interface, providing users with per-circuit fidelity estimates and qubit recommendations before they submit jobs.

3. **Research collaboration** (deeper engagement): I would be interested in contributing to the Quantum Inspire platform development, particularly in noise characterization, QEC decoder optimization, and hardware-software co-design. My background combines quantum information theory (Petz recovery maps, retrodiction) with practical software engineering, and I believe this combination could contribute to your team's efforts in scaling toward fault-tolerant quantum computing.

I am based in Taiwan but am open to remote collaboration or relocation. My long-term goal is to help bridge quantum information theory and hardware engineering, and Quantum Inspire is one of the platforms I most respect for making quantum computing accessible to researchers worldwide.

I would be grateful for the opportunity to discuss any of these possibilities.

Best regards,
Sheng-Kai Huang
akai@fawstudio.com
GitHub: https://github.com/akaiHuang

---

## 中文版（給台灣聯絡人參考）

**主旨：** τ-chrono 開源工具讓 Tuna-9 有效電路深度提升 50% — 合作提案

Quantum Inspire 團隊您好，

我是黃聖凱，台灣的獨立量子物理研究者。我一直在使用 Quantum Inspire 平台和 Tuna-9 處理器進行量子噪聲分析研究。

我開發了一個名為 τ-chrono 的開源工具，基於 Petz recovery map 理論，能比標準獨立閘模型更準確地預測電路噪聲。所有驗證數據都來自你們的 Tuna-9 硬體：

**NISQ 噪聲預測（Tuna-9 實測）：**
- Depth 50：比獨立模型準確 46%
- Bernstein-Vazirani：改善 60.4%
- H₂ VQE：Bayesian 追蹤在 depth 4 仍可用 (τ=0.49)，標準模型已判定失敗 (τ=0.60)

**QEC 噪聲感知解碼（Tuna-9 實測）：**
- 3-qubit repetition code：改善 39.9%

實際意義：Tuna-9 的用戶可以可靠地跑大約兩倍深的電路，不需要任何硬體改動。

工具完全開源（MIT 授權）：
- https://github.com/akaiHuang/petz-recovery-unification
- https://github.com/akaiHuang/qec-retrodiction-decoder

**我希望探討三個層次的合作：**

1. **文件推薦**：在 Quantum Inspire 文件中提及 τ-chrono，讓用戶受益。
2. **平台整合**：將噪聲預測整合進網頁介面。
3. **研究合作**：參與 Quantum Inspire 的硬體軟體協同設計。

我在台灣，但可以遠端合作或考慮前往荷蘭。

黃聖凱
akai@fawstudio.com
