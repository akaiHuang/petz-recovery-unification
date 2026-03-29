# tau-Chrono Engine: Concrete Implementation Plan

**Date**: 2026-03-16
**Author**: Sheng-Kai Huang
**Status**: Ready to Build

---

## 0. 你已經有什麼了？（Existing Assets Audit）

在開始寫任何新 code 之前，先盤點現有的 building blocks：

| 已有的 | 檔案 | 可直接用？ |
|--------|------|-----------|
| Petz recovery map (standard) | `numerical/petz_toolkit.py` | YES -- 完整實作，含 superoperator 建構 |
| Rotated Petz map (JRSWW) | `numerical/petz_toolkit.py` | YES -- 含 quadrature 積分，201 points |
| Fidelity (Uhlmann) | `numerical/petz_toolkit.py` | YES |
| Relative entropy | `numerical/petz_toolkit.py` | YES |
| tau_parameter (standard + rotated) | `numerical/petz_toolkit.py` | YES |
| delta_D computation | `numerical/petz_toolkit.py` | YES |
| Composition verification | `numerical/petz_toolkit.py` | YES -- `verify_composition()` |
| Saturation conditions check | `numerical/petz_toolkit.py` | YES -- `check_saturation_conditions()` |
| Random channel generation | `numerical/petz_toolkit.py` | YES -- `random_cptp_map()` |
| Common noise channels (AD, Dep, Deph) | 散落在 verification scripts | PARTIAL -- 需要整理成 library |

**結論**：核心數學引擎已經存在。`petz_toolkit.py` 就是 tau-Chrono 的心臟。你不需要重寫這些東西，只需要 wrap 它們成一個 EDA pipeline。

---

## 1. MVP: 第一週做什麼

### 1.1 先做哪個 Module？

**答案：Bayesian Composition Engine（Stage 2）**

理由：
- Saturation Scanner 很酷，但大多數 gate 不會恰好 saturate。它是 nice-to-have。
- Quantum Tax Calculator 需要 rotated Petz（慢，每個 gate 201 個積分點），不適合 MVP。
- Bayesian Composition Engine 是核心差異化。它直接展示「為什麼你的方法比 IBM 的 multiplicative fidelity 更好」。

### 1.2 最簡單的 Demo

**5-gate single-qubit circuit, 3 種 noise model 混合**

```
|+> --[AD(0.05)]--[Dep(0.08)]--[AD(0.12)]--[Deph(0.03)]--[AD(0.10)]--> measure
```

Demo 做三件事：
1. 計算 multiplicative fidelity (IBM 做法): `F_mult = prod(F_i)`
2. 計算 Bayesian composition bound: `sqrt(tau) <= sum(sqrt(tau_i^eff))`
3. 顯示 per-gate tau_eff < tau_naive，尤其是 dephasing gate（evolved state 接近 fixed point）

**這就是 design document Section 4 的數值例子，直接 code 化。**

### 1.3 "明天就能開始寫" 的 First Task

打開 terminal，執行：

```bash
mkdir -p /Users/akaihuangm1/Desktop/github/tau-chrono/src/tau_chrono
mkdir -p /Users/akaihuangm1/Desktop/github/tau-chrono/tests
mkdir -p /Users/akaihuangm1/Desktop/github/tau-chrono/notebooks
mkdir -p /Users/akaihuangm1/Desktop/github/tau-chrono/examples
```

然後把 `petz_toolkit.py` 複製過去作為 core engine（或 symlink），開始在上面加 EDA layer。

---

## 2. Tech Stack

```
Language:       Python 3.10+
Core deps:      numpy >= 1.24, scipy >= 1.10
Optional deps:  qiskit >= 1.0 (for circuit representation + noise models)
                matplotlib (for visualization)
                networkx (for DAG operations, optional -- 先用 list 就夠)
Dev deps:       pytest, jupyter
```

**不需要 Qiskit 才能開始。** MVP 用 raw Kraus operators，Qiskit 整合是 Phase 2。

### 安裝

```bash
cd /Users/akaihuangm1/Desktop/github/tau-chrono
python -m venv .venv
source .venv/bin/activate
pip install numpy scipy matplotlib pytest jupyter
# Qiskit is Phase 2:
# pip install qiskit qiskit-aer
```

---

## 3. File Structure

```
tau-chrono/
├── README.md                          # 用 design doc 的 Section 1 改寫
├── pyproject.toml                     # modern Python packaging
├── src/
│   └── tau_chrono/
│       ├── __init__.py                # 公開 API
│       ├── core.py                    # 從 petz_toolkit.py 搬來的核心函數
│       ├── channels.py                # 常見 noise channel 工廠函數
│       ├── composition.py             # Bayesian composition engine (Stage 2)
│       ├── saturation.py              # Saturation scanner (Stage 1)
│       ├── quantum_tax.py             # Quantum tax calculator (Stage 3)
│       ├── pipeline.py                # 三階段 pipeline 整合
│       ├── circuit.py                 # 簡易 circuit DAG 表示（不依賴 Qiskit）
│       └── visualization.py           # Chronometric path plots
├── tests/
│   ├── test_core.py                   # 測試 Petz map, fidelity, tau
│   ├── test_channels.py              # 測試 noise channel Kraus operators
│   ├── test_composition.py           # 測試 Bayesian composition bound
│   ├── test_saturation.py            # 測試 saturation detection
│   └── test_pipeline.py             # 端到端測試
├── notebooks/
│   ├── 01_quickstart.ipynb            # 5 分鐘 demo
│   ├── 02_composition_vs_multiplicative.ipynb  # 核心比較
│   └── 03_chronometric_path.ipynb     # 視覺化
├── examples/
│   ├── five_gate_demo.py              # Design doc Section 4 的完整實現
│   └── noise_budget_allocation.py     # Budget allocation demo
└── .gitignore
```

---

## 4. Core Algorithms — Implementation Priority

### Priority 1: channels.py (Day 1)

這是最容易做的，而且所有後續 demo 都需要它。

```python
"""
channels.py -- Standard quantum noise channel factories.

Each function returns a list of Kraus operators for a single-qubit channel.
"""

import numpy as np
from numpy.typing import NDArray

KrausList = list[NDArray[np.complexfloating]]


def amplitude_damping(gamma: float) -> KrausList:
    """
    Amplitude damping channel with parameter gamma in [0, 1].

    Kraus operators:
        E0 = [[1, 0], [0, sqrt(1-gamma)]]
        E1 = [[0, sqrt(gamma)], [0, 0]]

    Physical: spontaneous emission / T1 decay.
    """
    E0 = np.array([[1, 0], [0, np.sqrt(1 - gamma)]], dtype=complex)
    E1 = np.array([[0, np.sqrt(gamma)], [0, 0]], dtype=complex)
    return [E0, E1]


def depolarizing(p: float) -> KrausList:
    """
    Depolarizing channel: rho -> (1-p) rho + p/3 (X rho X + Y rho Y + Z rho Z).

    Kraus operators:
        E0 = sqrt(1 - 3p/4) I
        E1 = sqrt(p/4) X
        E2 = sqrt(p/4) Y
        E3 = sqrt(p/4) Z
    """
    I = np.eye(2, dtype=complex)
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    Z = np.array([[1, 0], [0, -1]], dtype=complex)

    E0 = np.sqrt(1 - 3 * p / 4) * I
    E1 = np.sqrt(p / 4) * X
    E2 = np.sqrt(p / 4) * Y
    E3 = np.sqrt(p / 4) * Z
    return [E0, E1, E2, E3]


def dephasing(p: float) -> KrausList:
    """
    Dephasing (phase damping) channel.

    Kraus operators:
        E0 = [[1, 0], [0, sqrt(1-p)]]
        E1 = [[0, 0], [0, sqrt(p)]]

    Physical: T2 dephasing (without energy loss).
    """
    E0 = np.array([[1, 0], [0, np.sqrt(1 - p)]], dtype=complex)
    E1 = np.array([[0, 0], [0, np.sqrt(p)]], dtype=complex)
    return [E0, E1]


def compose_kraus(kraus1: KrausList, kraus2: KrausList) -> KrausList:
    """
    Compose two channels: N2 o N1.

    Kraus operators of the composition: {K2_j @ K1_i} for all i, j.
    """
    return [K2 @ K1 for K2 in kraus2 for K1 in kraus1]


def thermal_state(beta: float, d: int = 2) -> NDArray:
    """
    Thermal state diag(p0, p1, ...) with p_k proportional to exp(-beta * k).

    For d=2: sigma = diag(e^0, e^{-beta}) / Z = diag(1, e^{-beta}) / (1 + e^{-beta}).

    Special case beta -> 0: maximally mixed.
    Special case beta -> inf: ground state |0><0|.
    """
    energies = np.arange(d, dtype=float)
    weights = np.exp(-beta * energies)
    weights /= np.sum(weights)
    return np.diag(weights.astype(complex))
```

### Priority 2: composition.py (Days 2-3)

這是 MVP 的核心。

```python
"""
composition.py -- Bayesian Composition Engine (Stage 2 of tau-Chrono).

The key insight: noise cost of gate G_i depends on the evolved reference
state from all preceding gates, not just the original reference.

Main function: bayesian_composition()
"""

import numpy as np
from dataclasses import dataclass
from typing import Optional

from .core import (
    apply_channel,
    tau_parameter,
    fidelity,
    apply_petz_recovery,
    petz_recovery_map,
)


@dataclass
class GateReport:
    """Per-gate analysis from the Bayesian composition engine."""
    gate_index: int
    channel_name: str
    tau_naive: float          # tau computed with original sigma (IBM's approach)
    tau_eff: float            # tau computed with Bayesian-updated sigma
    sqrt_tau_eff: float       # sqrt(tau_eff) -- the compositional contribution
    sigma_before: np.ndarray  # reference state entering this gate
    sigma_after: np.ndarray   # reference state leaving this gate
    saving_pct: float         # (tau_naive - tau_eff) / tau_naive * 100


@dataclass
class CompositionReport:
    """Full circuit analysis from the Bayesian composition engine."""
    gate_reports: list[GateReport]

    # Multiplicative approach (IBM baseline)
    tau_multiplicative: float          # 1 - prod(1 - tau_naive_i)
    sum_tau_naive: float               # sum(tau_naive_i) -- linear approximation

    # Bayesian approach
    sqrt_tau_bound: float              # (sum sqrt(tau_eff_i))^2
    sum_tau_eff: float                 # sum(tau_eff_i)

    # Comparison
    bayesian_saving_pct: float         # per-gate budget saving
    feasible: Optional[bool] = None    # whether tau_target is achievable
    tau_target: Optional[float] = None


def bayesian_composition(
    kraus_list: list,
    sigma_0: np.ndarray,
    rho_0: np.ndarray,
    channel_names: Optional[list[str]] = None,
    tau_target: Optional[float] = None,
) -> CompositionReport:
    """
    Run the Bayesian Composition Engine on a sequential circuit.

    Parameters
    ----------
    kraus_list : list of KrausList
        Kraus operators for each gate, in circuit order.
        kraus_list[i] is the list of Kraus operators for gate i.
    sigma_0 : ndarray
        Initial reference state.
    rho_0 : ndarray
        Initial input state.
    channel_names : list of str, optional
        Human-readable names for each gate (for reporting).
    tau_target : float, optional
        Target total tau. If provided, feasibility is checked.

    Returns
    -------
    CompositionReport
        Full analysis including per-gate budgets and comparison.
    """
    n = len(kraus_list)
    if channel_names is None:
        channel_names = [f"G{i+1}" for i in range(n)]

    # ---- Forward pass: propagate Bayesian references ----
    sigma = [None] * (n + 1)
    rho = [None] * (n + 1)
    sigma[0] = sigma_0.copy()
    rho[0] = rho_0.copy()

    for i in range(n):
        sigma[i + 1] = apply_channel(sigma[i], kraus_list[i])
        rho[i + 1] = apply_channel(rho[i], kraus_list[i])

    # ---- Compute per-gate tau ----
    gate_reports = []

    for i in range(n):
        # tau_naive: use original sigma_0 (what IBM does)
        tau_naive = tau_parameter(rho[i], kraus_list[i], sigma_0)

        # tau_eff: use Bayesian-updated sigma[i] (what we do)
        tau_eff = tau_parameter(rho[i], kraus_list[i], sigma[i])

        sqrt_tau_eff = np.sqrt(max(tau_eff, 0.0))

        saving = ((tau_naive - tau_eff) / tau_naive * 100) if tau_naive > 1e-15 else 0.0

        gate_reports.append(GateReport(
            gate_index=i,
            channel_name=channel_names[i],
            tau_naive=tau_naive,
            tau_eff=tau_eff,
            sqrt_tau_eff=sqrt_tau_eff,
            sigma_before=sigma[i],
            sigma_after=sigma[i + 1],
            saving_pct=saving,
        ))

    # ---- Aggregate metrics ----
    tau_naives = [g.tau_naive for g in gate_reports]
    tau_effs = [g.tau_eff for g in gate_reports]
    sqrt_tau_effs = [g.sqrt_tau_eff for g in gate_reports]

    # Multiplicative fidelity (IBM): F_total = prod(1 - tau_i)
    F_mult = 1.0
    for t in tau_naives:
        F_mult *= (1.0 - t)
    tau_multiplicative = 1.0 - F_mult

    # Sum of naive tau (linear approximation)
    sum_tau_naive = sum(tau_naives)

    # Bayesian bound: sqrt(tau_total) <= sum(sqrt(tau_i^eff))
    sqrt_tau_sum = sum(sqrt_tau_effs)
    sqrt_tau_bound = sqrt_tau_sum ** 2

    # Sum of effective tau
    sum_tau_eff = sum(tau_effs)

    # Saving
    bayesian_saving_pct = (
        (sum_tau_naive - sum_tau_eff) / sum_tau_naive * 100
        if sum_tau_naive > 1e-15 else 0.0
    )

    # Feasibility
    feasible = None
    if tau_target is not None:
        feasible = (sqrt_tau_bound <= tau_target)

    return CompositionReport(
        gate_reports=gate_reports,
        tau_multiplicative=tau_multiplicative,
        sum_tau_naive=sum_tau_naive,
        sqrt_tau_bound=sqrt_tau_bound,
        sum_tau_eff=sum_tau_eff,
        bayesian_saving_pct=bayesian_saving_pct,
        feasible=feasible,
        tau_target=tau_target,
    )


def print_composition_report(report: CompositionReport) -> None:
    """Pretty-print a composition report to stdout."""
    print("=" * 72)
    print("tau-Chrono Engine: Bayesian Composition Report")
    print("=" * 72)

    print(f"\n{'Gate':<12} {'tau_naive':>10} {'tau_eff':>10} "
          f"{'sqrt(tau)':>10} {'saving':>8}")
    print("-" * 55)

    for g in report.gate_reports:
        print(f"{g.channel_name:<12} {g.tau_naive:10.6f} {g.tau_eff:10.6f} "
              f"{g.sqrt_tau_eff:10.6f} {g.saving_pct:7.1f}%")

    print("-" * 55)
    print(f"\n--- Multiplicative Fidelity (IBM baseline) ---")
    print(f"  F_total = prod(1 - tau_i) => tau_mult = {report.tau_multiplicative:.6f}")
    print(f"  sum(tau_naive_i) = {report.sum_tau_naive:.6f}")

    print(f"\n--- Bayesian Composition (tau-Chrono) ---")
    print(f"  sum(tau_eff_i) = {report.sum_tau_eff:.6f}")
    print(f"  sqrt-tau bound = (sum sqrt(tau_eff_i))^2 = {report.sqrt_tau_bound:.6f}")
    print(f"  Per-gate budget saving: {report.bayesian_saving_pct:.1f}%")

    if report.feasible is not None:
        status = "FEASIBLE" if report.feasible else "INFEASIBLE"
        print(f"\n  Target tau = {report.tau_target:.6f} => {status}")

    print("=" * 72)
```

### Priority 3: saturation.py (Days 3-4)

```python
"""
saturation.py -- Saturation Scanner (Stage 1 of tau-Chrono).

Identifies gates where the Petz recovery is exact (tau_eff ~ 0),
removing them from the noise budget entirely.
"""

import numpy as np
from dataclasses import dataclass
from enum import Enum

from .core import (
    apply_channel,
    tau_parameter,
    fidelity,
)


class GateClass(Enum):
    SATURATED = "SATURATED"           # tau_eff < tol => free gate
    NEAR_SATURATED = "NEAR_SATURATED" # tau_eff < 10*tol => almost free
    NON_SATURATED = "NON_SATURATED"   # needs noise budget


@dataclass
class SaturationResult:
    gate_index: int
    classification: GateClass
    tau_eff: float
    commutator_norm: float       # ||[N(rho), N(sigma)]|| -- driver of non-saturation
    is_fixed_point: bool         # N(rho) ~ rho ?
    sigma_commutes: bool         # [rho, sigma] ~ 0 ?


def scan_saturation(
    kraus_list: list,
    sigma_0: np.ndarray,
    rho_0: np.ndarray,
    tol_saturated: float = 1e-6,
    tol_near: float = 1e-4,
    tol_commutator: float = 1e-6,
) -> list[SaturationResult]:
    """
    Scan each gate for saturation conditions.

    Saturation conditions (sufficient):
      (I)   rho is pure
      (II)  N(rho) ~ rho (fixed point)
      (III) [rho, sigma] ~ 0
      (IV)  [rho, N(sigma)] ~ 0

    When all hold: tau = 0 exactly => gate is "free" in the noise budget.

    Parameters
    ----------
    kraus_list : list of KrausList
        Kraus operators for each gate.
    sigma_0 : ndarray
        Initial reference state.
    rho_0 : ndarray
        Initial input state.
    tol_saturated : float
        Threshold for classifying as SATURATED.
    tol_near : float
        Threshold for NEAR_SATURATED.
    tol_commutator : float
        Threshold for commutator norm.

    Returns
    -------
    list of SaturationResult
    """
    n = len(kraus_list)
    results = []

    # Forward pass: propagate states
    sigma = sigma_0.copy()
    rho = rho_0.copy()

    for i in range(n):
        # Apply channel
        N_rho = apply_channel(rho, kraus_list[i])
        N_sigma = apply_channel(sigma, kraus_list[i])

        # Check fixed point: ||N(rho) - rho|| < tol
        is_fp = np.linalg.norm(N_rho - rho, 'fro') < tol_commutator

        # Check commutativity: ||[rho, sigma]||
        comm_rho_sigma = rho @ sigma - sigma @ rho
        sigma_commutes = np.linalg.norm(comm_rho_sigma, 'fro') < tol_commutator

        # Commutator of outputs: ||[N(rho), N(sigma)]||
        comm_outputs = N_rho @ N_sigma - N_sigma @ N_rho
        comm_norm = np.linalg.norm(comm_outputs, 'fro')

        # Compute actual tau_eff
        tau_eff = tau_parameter(rho, kraus_list[i], sigma)

        # Classify
        if tau_eff < tol_saturated:
            classification = GateClass.SATURATED
        elif tau_eff < tol_near:
            classification = GateClass.NEAR_SATURATED
        else:
            classification = GateClass.NON_SATURATED

        results.append(SaturationResult(
            gate_index=i,
            classification=classification,
            tau_eff=tau_eff,
            commutator_norm=comm_norm,
            is_fixed_point=is_fp,
            sigma_commutes=sigma_commutes,
        ))

        # Propagate for next gate
        sigma = N_sigma
        rho = N_rho

    return results
```

### Priority 4: quantum_tax.py (Days 4-5)

```python
"""
quantum_tax.py -- Quantum Tax Calculator (Stage 3 of tau-Chrono).

Computes the gap between standard Petz (composable) and rotated Petz
(per-gate optimal) for each gate. The gap quantifies the "cost of
Bayesian consistency".
"""

import numpy as np
from dataclasses import dataclass

from .core import (
    apply_channel,
    tau_parameter,
    tau_parameter_rotated,
)


@dataclass
class TaxReport:
    gate_index: int
    tau_petz: float           # standard Petz tau (composable)
    tau_rotated: float        # rotated Petz tau (per-gate optimal)
    gap: float                # tau_petz - tau_rotated >= 0
    recommendation: str       # "petz" or "rotated"


def compute_quantum_tax(
    kraus_list: list,
    sigma_0: np.ndarray,
    rho_0: np.ndarray,
    gap_threshold: float = 0.01,
    n_quadrature: int = 201,
) -> list[TaxReport]:
    """
    Compute the quantum tax (composability gap) for each gate.

    For each gate:
    - tau_petz: recovery quality with standard Petz (supports composition)
    - tau_rotated: recovery quality with rotated Petz (better per-gate)
    - gap = tau_petz - tau_rotated >= 0

    If gap is small: use Petz (get composition for free).
    If gap is large: consider rotated Petz at this gate (break composition).

    NOTE: This is the most expensive computation (rotated Petz requires
    numerical integration). Use n_quadrature=31 for fast estimates.
    """
    n = len(kraus_list)
    results = []

    sigma = sigma_0.copy()
    rho = rho_0.copy()

    for i in range(n):
        tau_p = tau_parameter(rho, kraus_list[i], sigma)
        tau_r = tau_parameter_rotated(
            rho, kraus_list[i], sigma, n_points=n_quadrature
        )
        gap = tau_p - tau_r

        rec = "rotated" if gap > gap_threshold else "petz"

        results.append(TaxReport(
            gate_index=i,
            tau_petz=tau_p,
            tau_rotated=tau_r,
            gap=max(gap, 0.0),
            recommendation=rec,
        ))

        # Propagate
        sigma = apply_channel(sigma, kraus_list[i])
        rho = apply_channel(rho, kraus_list[i])

    return results
```

### Priority 5: visualization.py (Day 5)

```python
"""
visualization.py -- Chronometric path visualization.

The "chronometric path" is the cumulative sqrt(tau) plot through the circuit.
It's the quantum analogue of arrival-time propagation in classical static
timing analysis.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure


def plot_chronometric_path(
    gate_names: list[str],
    sqrt_tau_effs: list[float],
    sqrt_tau_naives: list[float] = None,
    title: str = "Chronometric Path",
) -> Figure:
    """
    Plot cumulative sqrt(tau) through the circuit.

    X axis: gate index
    Y axis: cumulative sum of sqrt(tau_eff)

    Shows:
    - Blue line: Bayesian composition (tighter)
    - Red dashed: Naive/multiplicative (looser)
    - Shaded area: the Bayesian saving
    """
    n = len(gate_names)
    x = np.arange(n + 1)  # include initial point at 0

    # Cumulative sqrt-tau (Bayesian)
    cum_bayesian = np.zeros(n + 1)
    for i in range(n):
        cum_bayesian[i + 1] = cum_bayesian[i] + sqrt_tau_effs[i]

    fig, ax = plt.subplots(1, 1, figsize=(10, 5))

    ax.plot(x, cum_bayesian, 'b-o', linewidth=2, label='Bayesian sqrt(tau) path')

    if sqrt_tau_naives is not None:
        cum_naive = np.zeros(n + 1)
        for i in range(n):
            cum_naive[i + 1] = cum_naive[i] + sqrt_tau_naives[i]
        ax.plot(x, cum_naive, 'r--s', linewidth=2, label='Naive sqrt(tau) path')

        # Shade the saving
        ax.fill_between(x, cum_bayesian, cum_naive,
                        alpha=0.2, color='green', label='Bayesian saving')

    # Gate labels
    ax.set_xticks(x)
    xlabels = ['start'] + gate_names
    ax.set_xticklabels(xlabels, rotation=45, ha='right')

    ax.set_xlabel('Circuit Position')
    ax.set_ylabel('Cumulative sqrt(tau)')
    ax.set_title(title)
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    return fig


def plot_per_gate_comparison(
    gate_names: list[str],
    tau_naive: list[float],
    tau_eff: list[float],
    tau_rotated: list[float] = None,
    title: str = "Per-Gate Noise Budget Comparison",
) -> Figure:
    """
    Bar chart comparing naive vs Bayesian vs rotated tau per gate.
    """
    n = len(gate_names)
    x = np.arange(n)
    width = 0.25

    fig, ax = plt.subplots(1, 1, figsize=(10, 5))

    ax.bar(x - width, tau_naive, width, label='Naive (IBM)', color='#e74c3c', alpha=0.8)
    ax.bar(x, tau_eff, width, label='Bayesian (tau-Chrono)', color='#3498db', alpha=0.8)

    if tau_rotated is not None:
        ax.bar(x + width, tau_rotated, width, label='Rotated Petz (optimal)',
               color='#2ecc71', alpha=0.8)

    ax.set_xlabel('Gate')
    ax.set_ylabel('tau (process infidelity)')
    ax.set_title(title)
    ax.set_xticks(x)
    ax.set_xticklabels(gate_names)
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    return fig


def plot_reference_evolution(
    sigmas: list[np.ndarray],
    gate_names: list[str],
    title: str = "Reference State Evolution (Bayesian Updates)",
) -> Figure:
    """
    Plot how the reference state sigma evolves through the circuit.

    For single-qubit: shows sigma[0,0] (ground state population) at each step.
    """
    n = len(gate_names)
    populations = [sigma[0, 0].real for sigma in sigmas]

    fig, ax = plt.subplots(1, 1, figsize=(10, 4))

    x = np.arange(n + 1)
    xlabels = ['initial'] + gate_names

    ax.plot(x, populations, 'g-^', linewidth=2, markersize=8)
    ax.set_xticks(x)
    ax.set_xticklabels(xlabels, rotation=45, ha='right')

    ax.set_xlabel('Circuit Position')
    ax.set_ylabel('sigma[0,0] (ground state weight)')
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0.4, 1.0)

    plt.tight_layout()
    return fig
```

### Priority 6: pipeline.py (Day 5-6) — 把三個 Stage 串起來

```python
"""
pipeline.py -- Full tau-Chrono Engine pipeline.

Combines:
  Stage 1: Saturation Scanner
  Stage 2: Bayesian Composition Engine
  Stage 3: Quantum Tax Calculator
"""

from dataclasses import dataclass
from typing import Optional
import numpy as np

from .saturation import scan_saturation, SaturationResult, GateClass
from .composition import bayesian_composition, CompositionReport
from .quantum_tax import compute_quantum_tax, TaxReport


@dataclass
class EngineReport:
    """Complete tau-Chrono Engine output."""
    saturation: list[SaturationResult]
    composition: CompositionReport
    tax: list[TaxReport]

    # Summary statistics
    n_saturated: int
    n_near_saturated: int
    n_non_saturated: int


def run_engine(
    kraus_list: list,
    sigma_0: np.ndarray,
    rho_0: np.ndarray,
    channel_names: Optional[list[str]] = None,
    tau_target: Optional[float] = None,
    compute_tax: bool = True,
    tax_quadrature: int = 201,
) -> EngineReport:
    """
    Run the full tau-Chrono Engine pipeline.

    Parameters
    ----------
    kraus_list : list of KrausList
        Kraus operators for each gate.
    sigma_0 : ndarray
        Initial reference state.
    rho_0 : ndarray
        Initial input state.
    channel_names : list of str, optional
        Gate names for reporting.
    tau_target : float, optional
        Target tau for feasibility check.
    compute_tax : bool
        Whether to run Stage 3 (expensive). Default True.
    tax_quadrature : int
        Number of quadrature points for rotated Petz.

    Returns
    -------
    EngineReport
    """
    # Stage 1: Saturation Scanner
    sat_results = scan_saturation(kraus_list, sigma_0, rho_0)

    # Stage 2: Bayesian Composition
    comp_report = bayesian_composition(
        kraus_list, sigma_0, rho_0,
        channel_names=channel_names,
        tau_target=tau_target,
    )

    # Stage 3: Quantum Tax (optional)
    tax_results = []
    if compute_tax:
        tax_results = compute_quantum_tax(
            kraus_list, sigma_0, rho_0,
            n_quadrature=tax_quadrature,
        )

    # Summary
    n_sat = sum(1 for s in sat_results if s.classification == GateClass.SATURATED)
    n_near = sum(1 for s in sat_results if s.classification == GateClass.NEAR_SATURATED)
    n_non = sum(1 for s in sat_results if s.classification == GateClass.NON_SATURATED)

    return EngineReport(
        saturation=sat_results,
        composition=comp_report,
        tax=tax_results,
        n_saturated=n_sat,
        n_near_saturated=n_near,
        n_non_saturated=n_non,
    )
```

---

## 5. Demo Scenario: five_gate_demo.py

這就是 design document Section 4 的 code 化。把這個檔案放在 `examples/five_gate_demo.py`：

```python
#!/usr/bin/env python3
"""
five_gate_demo.py -- The canonical demo of the tau-Chrono Engine.

Reproduces the numerical example from the design document (Section 4):
  |+> --[AD(0.05)]--[Dep(0.08)]--[AD(0.12)]--[Deph(0.03)]--[AD(0.10)]-->

Shows:
  1. Multiplicative fidelity (IBM baseline)
  2. Bayesian composition bound (tau-Chrono)
  3. Saturation detection (G4 dephasing is nearly free)
  4. Chronometric path visualization
"""

import numpy as np
from tau_chrono.channels import (
    amplitude_damping, depolarizing, dephasing, thermal_state
)
from tau_chrono.pipeline import run_engine
from tau_chrono.composition import print_composition_report
from tau_chrono.visualization import (
    plot_chronometric_path,
    plot_per_gate_comparison,
    plot_reference_evolution,
)


def main():
    # ---- Circuit setup ----
    kraus_list = [
        amplitude_damping(0.05),    # G1
        depolarizing(0.08),         # G2
        amplitude_damping(0.12),    # G3
        dephasing(0.03),            # G4
        amplitude_damping(0.10),    # G5
    ]
    channel_names = [
        "AD(0.05)", "Dep(0.08)", "AD(0.12)", "Deph(0.03)", "AD(0.10)"
    ]

    # Reference: thermal state (beta=0.85 gives sigma ~ diag(0.7, 0.3))
    sigma_0 = thermal_state(beta=0.847)  # diag(0.700, 0.300)

    # Input: |+> state
    rho_0 = 0.5 * np.array([[1, 1], [1, 1]], dtype=complex)

    # ---- Run the engine ----
    report = run_engine(
        kraus_list, sigma_0, rho_0,
        channel_names=channel_names,
        tau_target=0.20,
        compute_tax=True,
        tax_quadrature=201,
    )

    # ---- Print results ----
    print_composition_report(report.composition)

    print("\n--- Saturation Analysis ---")
    for s in report.saturation:
        print(f"  {channel_names[s.gate_index]}: {s.classification.value}"
              f"  (tau_eff={s.tau_eff:.6f}, ||comm||={s.commutator_norm:.2e})")

    if report.tax:
        print("\n--- Quantum Tax Report ---")
        print(f"  {'Gate':<12} {'tau_Petz':>10} {'tau_rot':>10} {'gap':>10} {'rec':>8}")
        for t in report.tax:
            print(f"  {channel_names[t.gate_index]:<12} {t.tau_petz:10.6f} "
                  f"{t.tau_rotated:10.6f} {t.gap:10.6f} {t.recommendation:>8}")

    # ---- Visualizations ----
    gate_reports = report.composition.gate_reports

    fig1 = plot_chronometric_path(
        channel_names,
        [g.sqrt_tau_eff for g in gate_reports],
        [np.sqrt(max(g.tau_naive, 0)) for g in gate_reports],
        title="Chronometric Path: 5-Gate Mixed-Noise Circuit",
    )
    fig1.savefig("chronometric_path.png", dpi=150)

    fig2 = plot_per_gate_comparison(
        channel_names,
        [g.tau_naive for g in gate_reports],
        [g.tau_eff for g in gate_reports],
        [t.tau_rotated for t in report.tax] if report.tax else None,
        title="Per-Gate Noise Budget: Naive vs Bayesian vs Rotated",
    )
    fig2.savefig("per_gate_comparison.png", dpi=150)

    sigmas = [gate_reports[0].sigma_before] + [g.sigma_after for g in gate_reports]
    fig3 = plot_reference_evolution(sigmas, channel_names)
    fig3.savefig("reference_evolution.png", dpi=150)

    print("\nPlots saved: chronometric_path.png, per_gate_comparison.png, "
          "reference_evolution.png")


if __name__ == "__main__":
    main()
```

---

## 6. What NOT to Build (Scope Control)

| 不做 | 為什麼 |
|------|--------|
| Web UI / REST API | Library + notebooks 就夠了，先證明演算法價值 |
| Multi-qubit generalization | Composition theorem 目前只證了 sequential，tensor product 還是 open |
| Real hardware integration | IBM Brisbane 的 noise model 從 Qiskit 拿就好，Phase 2 再做 |
| Custom decoder construction | tau-Chrono 不是 decoder，它是 noise budget analyzer |
| Qiskit transpiler pass | Phase 3 的事，先把 standalone demo 做出來 |
| GPU acceleration | d=2 的矩陣運算，numpy 已經夠快了 |
| DAG with fan-in/fan-out | 先做 sequential (linear chain)，DAG 是 Phase 2 |

---

## 7. Timeline: 7-Day Sprint

假設一個人全職開發，每天 4-6 小時 coding：

### Day 1 (Mon): Project Setup + channels.py
- [ ] 建立 project structure (`mkdir`, `pyproject.toml`, `.gitignore`)
- [ ] 複製 `petz_toolkit.py` 作為 `core.py`，minor cleanup
- [ ] 實作 `channels.py` (AD, depolarizing, dephasing, thermal state)
- [ ] 寫 `test_channels.py`：驗證 Kraus operators 的 CPTP 性質
- **交付物**: `pip install -e .` 可以跑

### Day 2 (Tue): composition.py — 核心 MVP
- [ ] 實作 `bayesian_composition()` 函數
- [ ] 實作 `print_composition_report()`
- [ ] 寫 `test_composition.py`：
  - 兩個 channel 的 composition bound 成立
  - Bayesian tau_eff <= naive tau (至少在 AD → Deph 情況)
  - 5-gate demo 的數字跟 design doc Section 4 match
- **交付物**: `python -m pytest tests/test_composition.py` 全過

### Day 3 (Wed): saturation.py + integration
- [ ] 實作 `scan_saturation()`
- [ ] 寫 `test_saturation.py`：
  - Identity channel → SATURATED
  - Dephasing after AD → NEAR_SATURATED (for Z-diagonal states)
  - Strong depolarizing → NON_SATURATED
- [ ] 開始整合 `pipeline.py`
- **交付物**: Saturation scanner 能正確分類

### Day 4 (Thu): quantum_tax.py + visualization.py
- [ ] 實作 `compute_quantum_tax()`
- [ ] 實作 visualization functions (3 plots)
- [ ] 寫 `test_pipeline.py`：端到端測試
- **交付物**: 完整 pipeline 跑得通

### Day 5 (Fri): five_gate_demo.py + Notebook
- [ ] 完成 `examples/five_gate_demo.py`
- [ ] 寫 `notebooks/01_quickstart.ipynb`：5 分鐘上手
- [ ] 寫 `notebooks/02_composition_vs_multiplicative.ipynb`：深入比較
- [ ] 確保所有 plots 正確，數字跟 design doc match
- **交付物**: 3 個可執行的 demo artifacts

### Day 6 (Sat): Edge Cases + Robustness
- [ ] Test edge cases：
  - sigma 接近 singular (low rank)
  - gamma/p 接近 0 或 1
  - Identity channel (所有 gates)
  - 單一 gate 的退化情況
- [ ] 加入 input validation + error messages
- [ ] 跑一次 `random_cptp_map` 的 statistical test (1000 samples)
- **交付物**: 所有 tests pass，包括 edge cases

### Day 7 (Sun): Polish + README
- [ ] 清理 code，加 docstrings
- [ ] 寫 `notebooks/03_chronometric_path.ipynb`
- [ ] 初始化 git repo，第一個 commit
- [ ] 考慮 CLI entry point: `python -m tau_chrono run --circuit five_gate.json`
- **交付物**: 可以 demo 給別人看的 v0.1

---

## 8. Phase 2 Roadmap (After MVP, Weeks 2-4)

| Week | Task | Detail |
|------|------|--------|
| Week 2 | Qiskit Integration | `from_qiskit_circuit()` converter: 從 QuantumCircuit + NoiseModel 自動抽出 Kraus operators |
| Week 2 | IBM Noise Models | 用 `qiskit_aer.noise` 的 real device models (Brisbane, Osaka) |
| Week 3 | DAG Support | 用 networkx 表示 circuit DAG，支援 parallel gates |
| Week 3 | Budget Optimizer | 給定 tau_target，找最佳 per-gate budget allocation (convex optimization) |
| Week 4 | Benchmark Suite | 跑 QFT, QAOA, VQE 等標準 circuits，比較 multiplicative vs Bayesian |

### Phase 3 (Weeks 5-8): Qiskit Transpiler Pass

```python
# 未來的 API -- 作為 Qiskit transpiler pass
from qiskit.transpiler import PassManager
from tau_chrono.qiskit import TauChronoAnalysis, TauChronoRouting

pm = PassManager([
    TauChronoAnalysis(noise_model=backend.noise_model),
    TauChronoRouting(tau_target=0.05),
])

optimized = pm.run(circuit)
```

---

## 9. Key Technical Decisions

### Q: 為什麼不直接在 Qiskit 裡面改？
**A**: 因為 Qiskit 用 Pauli channel approximation（只存 Pauli error rates），而 Bayesian composition 需要完整的 Kraus operators。等 MVP 驗證了價值，再做 Qiskit integration。

### Q: Rotated Petz 太慢怎麼辦？
**A**: Stage 3 (Quantum Tax) 是 optional。用 `compute_tax=False` 可以跳過。而且對 single-qubit (d=2)，201 個 quadrature point 的積分大約 0.1 秒/gate，5-gate circuit 總共 0.5 秒。可以接受。如果要加速，先降到 `n_quadrature=31`。

### Q: 為什麼 MVP 只做 single-qubit？
**A**: 因為 d=2 的所有矩陣運算都是 2x2，eigendecomposition 是 O(1)。跳到 d=4 (two-qubit) 也沒問題（O(64) 的 superoperator），但 multi-qubit 的 tensor product composition 在數學上還是 open question (design doc Section 7.3.1)。

### Q: sigma_0 怎麼選？
**A**: MVP 用 thermal state。更好的選法是 `sigma_0 = maximally mixed state` (最保守) 或 `sigma_0 = estimated prior from previous calibration` (最準確)。Adaptive reference selection 是 Phase 3 的 open optimization problem。

---

## 10. Success Criteria

MVP 做完之後，你應該能展示這些東西：

1. **一張 chronometric path 圖**：藍線 (Bayesian) 明顯低於紅線 (naive)
2. **一個數字**：「Bayesian composition 在 5-gate mixed-noise circuit 上節省 X% 的 noise budget」
3. **一個 saturation detection**：G4 (dephasing) 被自動標記為 NEAR_SATURATED，budget 從 0.015 降到 ~0.002
4. **一張 bar chart**：per-gate tau 三種方法的比較
5. **一個 feasibility answer**：「Given tau_target = 0.20, this circuit is FEASIBLE/INFEASIBLE」

如果這五個東西都能跑出來，MVP 就成功了。
