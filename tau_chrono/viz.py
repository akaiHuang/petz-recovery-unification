"""
Visualization module for tau-chrono quantum noise analysis.

Provides publication-quality matplotlib figures for inspecting Bayesian
composition results, channel diagnostics, and depth-scaling behaviour.

All ``plot_*`` functions return a ``matplotlib.figure.Figure`` object
without calling ``plt.show()``, making them compatible with Jupyter
notebooks, PDF export, and headless environments.

Public API
----------
- ``plot_tau_heatmap``            : per-gate tau comparison bar chart
- ``plot_depth_scaling``          : tau vs circuit depth line plot
- ``plot_composition_inequality`` : LHS vs RHS composition check
- ``plot_sigma_evolution``        : eigenvalue trajectories of sigma
- ``plot_ground_truth_comparison``: actual vs predicted fidelity bars
- ``plot_channel_diagnostics``    : multi-panel Choi / eigenvalue / pie
"""

from __future__ import annotations

from typing import List, Optional, Sequence, TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

if TYPE_CHECKING:
    from matplotlib.figure import Figure

try:
    import matplotlib
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    from matplotlib.figure import Figure as _Figure
    from matplotlib.patches import FancyBboxPatch
    HAS_MPL = True
except ImportError:  # pragma: no cover
    HAS_MPL = False

from .bayesian import CompositionResult, GateResult
from .channels import (
    _choi_from_kraus,
    entanglement_fidelity,
    coherent_incoherent_decomposition,
)

# Type aliases
KrausList = List[NDArray[np.complexfloating]]


# ---------------------------------------------------------------------------
# Style configuration
# ---------------------------------------------------------------------------

_PALETTE = {
    "NORMAL": "#2ecc71",           # green
    "NEAR_SATURATED": "#f39c12",   # amber
    "SATURATED": "#e74c3c",        # red
}

_STYLE_DEFAULTS = {
    "font.family": "sans-serif",
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.titleweight": "bold",
    "axes.labelsize": 11,
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "figure.dpi": 120,
    "figure.facecolor": "white",
    "savefig.dpi": 200,
    "savefig.bbox": "tight",
}


def _require_mpl() -> None:
    """Raise a clear error when matplotlib is not installed."""
    if not HAS_MPL:
        raise ImportError(
            "matplotlib is required for tau_chrono.viz. "
            "Install it with:  pip install matplotlib"
        )


def _apply_style() -> None:
    """Apply the tau-chrono house style to matplotlib."""
    plt.rcParams.update(_STYLE_DEFAULTS)


def _classification_color(cls: str) -> str:
    """Map a gate classification string to its colour."""
    return _PALETTE.get(cls, "#95a5a6")


# ---------------------------------------------------------------------------
# 1. Per-gate tau heatmap (horizontal bar chart)
# ---------------------------------------------------------------------------


def plot_tau_heatmap(result: CompositionResult) -> "Figure":
    """Horizontal bar chart of per-gate ``tau_naive`` vs ``tau_eff``.

    Each gate gets a pair of bars, colour-coded by classification:
    green (NORMAL), amber (NEAR_SATURATED), red (SATURATED).
    The improvement percentage is annotated beside each pair.

    Parameters
    ----------
    result : CompositionResult
        Output of ``bayesian_compose``.

    Returns
    -------
    Figure
        Matplotlib figure.
    """
    _require_mpl()
    _apply_style()

    gates: List[GateResult] = result.gate_results
    n = len(gates)

    fig, ax = plt.subplots(figsize=(8, max(2.5, 0.7 * n + 1.2)))

    labels = [g.channel_name for g in gates]
    tau_naive = np.array([g.tau_naive for g in gates])
    tau_eff = np.array([g.tau_eff for g in gates])
    colors = [_classification_color(g.classification) for g in gates]

    y = np.arange(n)
    bar_h = 0.35

    # Naive bars (lighter, behind)
    ax.barh(
        y + bar_h / 2, tau_naive, height=bar_h,
        color=[mcolors.to_rgba(c, alpha=0.35) for c in colors],
        edgecolor=[mcolors.to_rgba(c, alpha=0.6) for c in colors],
        linewidth=0.8, label=r"$\tau_{\rm naive}$",
    )

    # Effective bars (solid, in front)
    ax.barh(
        y - bar_h / 2, tau_eff, height=bar_h,
        color=colors, edgecolor="white", linewidth=0.5,
        label=r"$\tau_{\rm eff}$ (Bayesian)",
    )

    # Improvement annotations
    x_max = max(np.max(tau_naive), np.max(tau_eff)) if n > 0 else 1.0
    for i, g in enumerate(gates):
        if g.tau_naive > 1e-15:
            improv = (1.0 - g.tau_eff / g.tau_naive) * 100.0
            sign = "+" if improv < 0 else ""  # negative = worse
            text = f"{improv:+.1f}%" if abs(improv) > 0.05 else "0.0%"
        else:
            text = "--"
        ax.text(
            x_max * 1.03, i, text,
            va="center", ha="left", fontsize=9,
            color="#2c3e50", fontweight="medium",
        )

    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlabel(r"$\tau$ (recovery failure)")
    ax.set_title("Per-gate noise: naive vs Bayesian")
    ax.legend(loc="lower right", framealpha=0.9, fontsize=9)
    ax.set_xlim(left=0, right=x_max * 1.18)
    ax.invert_yaxis()

    # Classification legend
    for cls, col in _PALETTE.items():
        ax.plot([], [], "s", color=col, markersize=7, label=cls)

    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# 2. Depth-scaling plot
# ---------------------------------------------------------------------------


def plot_depth_scaling(
    depths: Sequence[int],
    tau_naive_list: Sequence[float],
    tau_bayesian_list: Sequence[float],
    *,
    log_scale: bool = False,
    title: str = "Noise scaling with circuit depth",
) -> "Figure":
    """Line plot of total ``tau`` vs circuit depth for naive and Bayesian.

    A shaded region highlights the improvement that Bayesian tracking
    provides over the naive multiplicative estimate.

    Parameters
    ----------
    depths : sequence of int
        Circuit depths (x-axis).
    tau_naive_list : sequence of float
        Naive (multiplicative) total tau at each depth.
    tau_bayesian_list : sequence of float
        Bayesian-composed total tau at each depth.
    log_scale : bool
        If True, use logarithmic y-axis.
    title : str
        Plot title.

    Returns
    -------
    Figure
    """
    _require_mpl()
    _apply_style()

    depths = np.asarray(depths)
    tau_n = np.asarray(tau_naive_list)
    tau_b = np.asarray(tau_bayesian_list)

    fig, ax = plt.subplots(figsize=(7, 4.5))

    ax.plot(
        depths, tau_n, "o--", color="#e74c3c", linewidth=1.6,
        markersize=5, label=r"$\tau_{\rm naive}$ (multiplicative)",
    )
    ax.plot(
        depths, tau_b, "s-", color="#2980b9", linewidth=1.8,
        markersize=5, label=r"$\tau_{\rm Bayesian}$",
    )

    # Shaded improvement region
    ax.fill_between(
        depths, tau_b, tau_n,
        where=tau_n >= tau_b,
        color="#2980b9", alpha=0.12,
        label="Bayesian improvement",
    )

    if log_scale:
        ax.set_yscale("log")
        ax.set_ylabel(r"$\tau$ (log scale)")
    else:
        ax.set_ylabel(r"$\tau$ (recovery failure)")

    ax.set_xlabel("Circuit depth")
    ax.set_title(title)
    ax.legend(framealpha=0.9, fontsize=9)

    # Annotate max improvement
    if len(depths) > 0:
        diff = tau_n - tau_b
        idx_max = int(np.argmax(diff))
        if tau_n[idx_max] > 1e-15:
            pct = diff[idx_max] / tau_n[idx_max] * 100.0
            ax.annotate(
                f"max {pct:.0f}% improvement",
                xy=(depths[idx_max], tau_b[idx_max]),
                xytext=(depths[idx_max], (tau_n[idx_max] + tau_b[idx_max]) / 2),
                ha="center", fontsize=8, color="#2c3e50",
                arrowprops=dict(arrowstyle="->", color="#7f8c8d", lw=0.8),
            )

    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# 3. Composition inequality check
# ---------------------------------------------------------------------------


def plot_composition_inequality(result: CompositionResult) -> "Figure":
    """Bar chart comparing the composition inequality LHS and RHS.

    Displays ``sqrt(tau_total)`` vs ``sum sqrt(tau_i^eff)`` with a
    pass/fail indicator and the slack value.

    Parameters
    ----------
    result : CompositionResult
        Output of ``bayesian_compose``.

    Returns
    -------
    Figure
    """
    _require_mpl()
    _apply_style()

    lhs = result.composition_lhs
    rhs = result.composition_rhs
    holds = result.composition_holds
    slack = result.composition_slack

    fig, ax = plt.subplots(figsize=(5, 4))

    bar_colors = ["#2980b9", "#27ae60" if holds else "#e74c3c"]
    bars = ax.bar(
        [r"LHS: $\sqrt{\tau_{\rm total}}$",
         r"RHS: $\sum_i \sqrt{\tau_i^{\rm eff}}$"],
        [lhs, rhs],
        color=bar_colors,
        edgecolor="white", linewidth=1.2, width=0.55,
    )

    # Value labels on bars
    for bar, val in zip(bars, [lhs, rhs]):
        ax.text(
            bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.005,
            f"{val:.4f}", ha="center", va="bottom", fontsize=10,
            fontweight="bold", color="#2c3e50",
        )

    # Pass / fail banner
    status_text = "PASS" if holds else "FAIL"
    status_color = "#27ae60" if holds else "#e74c3c"
    ax.text(
        0.98, 0.95, status_text,
        transform=ax.transAxes, ha="right", va="top",
        fontsize=14, fontweight="bold", color="white",
        bbox=dict(
            boxstyle="round,pad=0.35", facecolor=status_color,
            edgecolor="none", alpha=0.9,
        ),
    )

    # Slack annotation
    ax.text(
        0.98, 0.82,
        f"slack = {slack:.4f}",
        transform=ax.transAxes, ha="right", va="top",
        fontsize=9, color="#7f8c8d",
    )

    ax.set_ylabel("Value")
    ax.set_title(r"Composition inequality: $\sqrt{\tau} \leq \sum_i \sqrt{\tau_i}$")
    ax.set_ylim(bottom=0, top=max(lhs, rhs) * 1.25 + 0.01)

    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# 4. Sigma eigenvalue evolution
# ---------------------------------------------------------------------------


def plot_sigma_evolution(result: CompositionResult) -> "Figure":
    """Plot eigenvalues of sigma at each gate step.

    Tracks how the reference state evolves through the circuit, showing
    convergence toward a fixed point.

    Parameters
    ----------
    result : CompositionResult
        Output of ``bayesian_compose``.

    Returns
    -------
    Figure
    """
    _require_mpl()
    _apply_style()

    gates: List[GateResult] = result.gate_results
    n = len(gates)

    # Collect eigenvalues at each step (sigma_before for each gate, plus
    # the sigma_after of the last gate inferred from rho_after structure).
    steps = []
    eigvals_per_step = []

    for i, g in enumerate(gates):
        eigs = np.sort(np.linalg.eigvalsh(g.sigma_before))[::-1]
        steps.append(i)
        eigvals_per_step.append(eigs)

    d = eigvals_per_step[0].shape[0] if n > 0 else 2
    eigvals_arr = np.array(eigvals_per_step)  # shape (n, d)

    fig, ax = plt.subplots(figsize=(7, 4))

    # Use a colour palette for eigenvalues
    cmap = matplotlib.colormaps.get_cmap("viridis").resampled(d)
    for k in range(d):
        ax.plot(
            steps, eigvals_arr[:, k],
            "o-", color=cmap(k), linewidth=1.5, markersize=4,
            label=rf"$\lambda_{k}$",
        )

    # Fixed-point reference: maximally mixed state
    ax.axhline(
        1.0 / d, color="#95a5a6", linestyle=":", linewidth=1.0,
        label=f"maximally mixed (1/{d})",
    )

    ax.set_xlabel("Gate index")
    ax.set_ylabel(r"Eigenvalue of $\sigma$")
    ax.set_title(r"Reference state $\sigma$ evolution through circuit")
    ax.legend(fontsize=8, framealpha=0.9, ncol=min(d + 1, 4))
    ax.set_ylim(bottom=-0.02, top=1.05)

    if n > 0:
        ax.set_xticks(steps)
        ax.set_xticklabels([g.channel_name for g in gates], rotation=45, ha="right")

    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# 5. Ground-truth comparison
# ---------------------------------------------------------------------------


def plot_ground_truth_comparison(
    circuits: Sequence[str],
    F_actual: Sequence[float],
    F_naive: Sequence[float],
    F_bayesian: Sequence[float],
    *,
    title: str = "Fidelity prediction accuracy",
) -> "Figure":
    """Grouped bar chart comparing actual vs predicted fidelities.

    Parameters
    ----------
    circuits : sequence of str
        Circuit labels (x-axis).
    F_actual : sequence of float
        Ground-truth fidelities (from full simulation or experiment).
    F_naive : sequence of float
        Naive (multiplicative) predicted fidelities.
    F_bayesian : sequence of float
        Bayesian predicted fidelities.
    title : str
        Plot title.

    Returns
    -------
    Figure
    """
    _require_mpl()
    _apply_style()

    n = len(circuits)
    x = np.arange(n)
    width = 0.25

    F_act = np.asarray(F_actual)
    F_nai = np.asarray(F_naive)
    F_bay = np.asarray(F_bayesian)

    fig, ax = plt.subplots(figsize=(max(6, 1.5 * n + 2), 4.5))

    bars_actual = ax.bar(
        x - width, F_act, width,
        color="#2c3e50", edgecolor="white", linewidth=0.5,
        label="Actual $F$",
    )
    bars_naive = ax.bar(
        x, F_nai, width,
        color="#e74c3c", edgecolor="white", linewidth=0.5,
        label=r"Naive $\hat{F}$",
    )
    bars_bayes = ax.bar(
        x + width, F_bay, width,
        color="#2980b9", edgecolor="white", linewidth=0.5,
        label=r"Bayesian $\hat{F}$",
    )

    # Difference annotations above grouped bars
    for i in range(n):
        err_naive = abs(F_nai[i] - F_act[i])
        err_bayes = abs(F_bay[i] - F_act[i])
        y_top = max(F_act[i], F_nai[i], F_bay[i]) + 0.015
        ax.text(
            x[i], y_top,
            f"$\\Delta$N={err_naive:.3f}\n$\\Delta$B={err_bayes:.3f}",
            ha="center", va="bottom", fontsize=7, color="#7f8c8d",
        )

    ax.set_xticks(x)
    ax.set_xticklabels(circuits, rotation=30, ha="right")
    ax.set_ylabel("Fidelity $F$")
    ax.set_title(title)
    ax.legend(loc="lower left", framealpha=0.9, fontsize=9)
    ax.set_ylim(bottom=0, top=min(1.15, np.max(F_act) * 1.25 + 0.05))

    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# 6. Channel diagnostics (multi-panel)
# ---------------------------------------------------------------------------


def plot_channel_diagnostics(
    kraus_ops: KrausList,
    name: str = "Channel",
) -> "Figure":
    """Multi-panel channel diagnostic plot.

    Panel layout (1 row, 3 columns):

    1. **Choi matrix** -- heatmap of ``|J_N|`` with eigenvalue annotations.
    2. **Eigenvalue spectrum** -- bar chart of Choi eigenvalues (positive
       eigenvalues in blue, negative in red for non-CP detection).
    3. **Coherent / incoherent decomposition** -- pie chart breaking down
       the channel error budget.

    Works for both 1-qubit (2x2 Kraus) and 2-qubit (4x4 Kraus) channels.

    Parameters
    ----------
    kraus_ops : KrausList
        Kraus operators defining the channel.
    name : str
        Label for the channel (used in titles).

    Returns
    -------
    Figure
    """
    _require_mpl()
    _apply_style()

    d = kraus_ops[0].shape[0]
    choi = _choi_from_kraus(kraus_ops)
    choi_eigvals = np.sort(np.linalg.eigvalsh(choi))[::-1]
    decomp = coherent_incoherent_decomposition(kraus_ops)
    f_e = entanglement_fidelity(kraus_ops)

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))

    # --- Panel 1: Choi matrix heatmap ---
    ax1 = axes[0]
    choi_abs = np.abs(choi)
    im = ax1.imshow(choi_abs, cmap="inferno", interpolation="nearest", aspect="equal")
    fig.colorbar(im, ax=ax1, fraction=0.046, pad=0.04, label="|value|")
    ax1.set_title(f"Choi matrix $|J_N|$\n({name})")
    ax1.set_xlabel("Column")
    ax1.set_ylabel("Row")

    # Tick labels for small matrices
    if d * d <= 16:
        ax1.set_xticks(range(d * d))
        ax1.set_yticks(range(d * d))
    else:
        ax1.set_xticks([0, d * d - 1])
        ax1.set_yticks([0, d * d - 1])

    # --- Panel 2: Choi eigenvalue spectrum ---
    ax2 = axes[1]
    n_eigs = len(choi_eigvals)
    eig_colors = ["#2980b9" if v >= -1e-12 else "#e74c3c" for v in choi_eigvals]
    bars = ax2.bar(
        range(n_eigs), choi_eigvals,
        color=eig_colors, edgecolor="white", linewidth=0.5,
    )
    ax2.axhline(0, color="#bdc3c7", linewidth=0.6)
    ax2.set_xlabel("Eigenvalue index")
    ax2.set_ylabel("Eigenvalue")
    ax2.set_title(f"Choi eigenvalue spectrum\n$F_e$ = {f_e:.4f}")
    ax2.set_xticks(range(n_eigs))

    # Annotate dominant eigenvalue
    if n_eigs > 0:
        idx_max = int(np.argmax(choi_eigvals))
        ax2.annotate(
            f"{choi_eigvals[idx_max]:.3f}",
            xy=(idx_max, choi_eigvals[idx_max]),
            xytext=(idx_max, choi_eigvals[idx_max] + 0.05 * max(abs(choi_eigvals[0]), 0.1)),
            ha="center", fontsize=8, color="#2c3e50",
            arrowprops=dict(arrowstyle="->", color="#7f8c8d", lw=0.7),
        )

    # --- Panel 3: Coherent / incoherent pie ---
    ax3 = axes[2]
    coh = decomp["coherent_frac"]
    inc = decomp["incoherent_frac"]
    f_avg = decomp["avg_gate_fidelity"]
    unitarity = decomp["unitarity"]

    # Handle near-identity channels where total error is negligible
    total_error = 1.0 - f_avg
    if total_error < 1e-10:
        # Near-perfect channel -- show a single "no error" wedge
        ax3.pie(
            [1.0], labels=["No error"], colors=["#27ae60"],
            startangle=90, textprops={"fontsize": 10},
        )
    else:
        sizes = [coh, inc]
        labels_pie = [
            f"Coherent\n{coh * 100:.1f}%",
            f"Incoherent\n{inc * 100:.1f}%",
        ]
        pie_colors = ["#e67e22", "#3498db"]

        wedges, texts = ax3.pie(
            sizes, labels=labels_pie, colors=pie_colors,
            startangle=90, textprops={"fontsize": 9},
            wedgeprops={"edgecolor": "white", "linewidth": 1.5},
        )

    ax3.set_title(f"Error decomposition\n$F_{{\\rm avg}}$={f_avg:.4f}  U={unitarity:.4f}")

    fig.suptitle(f"Channel diagnostics: {name}", fontsize=13, fontweight="bold", y=1.02)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Convenience: summary dashboard
# ---------------------------------------------------------------------------


def plot_summary_dashboard(result: CompositionResult) -> "Figure":
    """Four-panel summary dashboard for a composition result.

    Combines ``plot_tau_heatmap``, ``plot_composition_inequality``,
    and ``plot_sigma_evolution`` into a single figure, plus a text
    panel with key metrics.

    Parameters
    ----------
    result : CompositionResult
        Output of ``bayesian_compose``.

    Returns
    -------
    Figure
    """
    _require_mpl()
    _apply_style()

    gates = result.gate_results
    n = len(gates)

    fig = plt.figure(figsize=(14, max(8, 3 * n * 0.3 + 4)),
                      layout="constrained")
    gs = fig.add_gridspec(2, 2, hspace=0.08, wspace=0.08)

    # --- Top-left: tau heatmap ---
    ax_tau = fig.add_subplot(gs[0, 0])
    labels = [g.channel_name for g in gates]
    tau_naive = np.array([g.tau_naive for g in gates])
    tau_eff = np.array([g.tau_eff for g in gates])
    colors = [_classification_color(g.classification) for g in gates]
    y = np.arange(n)
    bar_h = 0.35

    ax_tau.barh(
        y + bar_h / 2, tau_naive, height=bar_h,
        color=[mcolors.to_rgba(c, alpha=0.35) for c in colors],
        edgecolor=[mcolors.to_rgba(c, alpha=0.6) for c in colors],
        linewidth=0.8, label=r"$\tau_{\rm naive}$",
    )
    ax_tau.barh(
        y - bar_h / 2, tau_eff, height=bar_h,
        color=colors, edgecolor="white", linewidth=0.5,
        label=r"$\tau_{\rm eff}$",
    )
    ax_tau.set_yticks(y)
    ax_tau.set_yticklabels(labels)
    ax_tau.set_xlabel(r"$\tau$")
    ax_tau.set_title("Per-gate tau")
    ax_tau.invert_yaxis()
    ax_tau.legend(fontsize=8)

    # --- Top-right: composition inequality ---
    ax_comp = fig.add_subplot(gs[0, 1])
    bar_colors = [
        "#2980b9",
        "#27ae60" if result.composition_holds else "#e74c3c",
    ]
    bars = ax_comp.bar(
        [r"$\sqrt{\tau_{\rm tot}}$", r"$\sum\sqrt{\tau_i}$"],
        [result.composition_lhs, result.composition_rhs],
        color=bar_colors, edgecolor="white", width=0.5,
    )
    for bar, val in zip(bars, [result.composition_lhs, result.composition_rhs]):
        ax_comp.text(
            bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.003,
            f"{val:.4f}", ha="center", va="bottom", fontsize=9, fontweight="bold",
        )
    status = "PASS" if result.composition_holds else "FAIL"
    sc = "#27ae60" if result.composition_holds else "#e74c3c"
    ax_comp.text(
        0.95, 0.92, status, transform=ax_comp.transAxes, ha="right", va="top",
        fontsize=12, fontweight="bold", color="white",
        bbox=dict(boxstyle="round,pad=0.3", facecolor=sc, edgecolor="none"),
    )
    ax_comp.set_title("Composition inequality")
    ax_comp.set_ylim(
        bottom=0,
        top=max(result.composition_lhs, result.composition_rhs) * 1.3 + 0.01,
    )

    # --- Bottom-left: sigma evolution ---
    ax_sig = fig.add_subplot(gs[1, 0])
    if n > 0:
        d = gates[0].sigma_before.shape[0]
        eigvals_arr = np.array([
            np.sort(np.linalg.eigvalsh(g.sigma_before))[::-1] for g in gates
        ])
        cmap = matplotlib.colormaps.get_cmap("viridis").resampled(d)
        for k in range(d):
            ax_sig.plot(
                range(n), eigvals_arr[:, k],
                "o-", color=cmap(k), linewidth=1.3, markersize=3,
                label=rf"$\lambda_{k}$",
            )
        ax_sig.axhline(1.0 / d, color="#95a5a6", ls=":", lw=0.8)
        ax_sig.set_xticks(range(n))
        ax_sig.set_xticklabels(labels, rotation=45, ha="right", fontsize=7)
    ax_sig.set_ylabel(r"$\lambda(\sigma)$")
    ax_sig.set_title(r"$\sigma$ eigenvalue evolution")

    # --- Bottom-right: summary metrics ---
    ax_text = fig.add_subplot(gs[1, 1])
    ax_text.axis("off")
    lines = [
        f"Gates:  {n}",
        f"tau_Bayesian:  {result.tau_bayesian_total:.6f}",
        f"tau_naive:     {result.tau_multiplicative_total:.6f}",
        f"Improvement:   {result.improvement_percent:.1f}%",
        f"Slack:         {result.composition_slack:.6f}",
        f"Inequality:    {'HOLDS' if result.composition_holds else 'VIOLATED'}",
    ]
    for i, g in enumerate(gates):
        cls_sym = {"NORMAL": "o", "NEAR_SATURATED": "~", "SATURATED": "!"}
        sym = cls_sym.get(g.classification, "?")
        lines.append(f"  [{sym}] {g.channel_name}: tau_eff={g.tau_eff:.4e}")

    text_block = "\n".join(lines)
    ax_text.text(
        0.05, 0.95, text_block,
        transform=ax_text.transAxes, va="top", ha="left",
        fontsize=9, fontfamily="monospace",
        bbox=dict(boxstyle="round,pad=0.5", facecolor="#ecf0f1", edgecolor="#bdc3c7"),
    )
    ax_text.set_title("Summary")

    fig.suptitle(
        "tau-chrono Bayesian Composition Report",
        fontsize=14, fontweight="bold",
    )
    return fig
