"""
dashboard.py -- Generate interactive HTML dashboards for Sigma diffusion experiments.

Uses Chart.js (loaded from CDN) for interactive charts.  Dark theme
matches the existing docs/index.html style (Apple-inspired dark UI).
Images can be embedded as base64 data URIs.
"""

from __future__ import annotations

import base64
import html as html_lib
import json
import os
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import numpy as np

# ===================================================================
# Public API
# ===================================================================


def generate_dashboard(
    results_dict: Dict[str, dict],
    output_path: Union[str, Path],
) -> str:
    """Generate a self-contained ``results.html`` dashboard.

    Parameters
    ----------
    results_dict : dict
        Mapping from experiment name (str) to a dict with keys:

        * ``sigma_profile`` : dict with ``timesteps`` (list[int]),
          ``sigma_per_step`` (list[float]), ``cumulative_sigma``
          (list[float]), ``fidelity_bound`` (list[float]).
        * ``total_sigma`` : float
        * ``speedup_ratio`` : float  (optional, default 1.0)
        * ``quality_loss_estimate`` : float  (optional)
        * ``psnr`` : float  (optional)
        * ``ssim`` : float  (optional)
        * ``lpips`` : float  (optional)
        * ``num_active_steps`` : int  (optional)
        * ``num_total_steps`` : int  (optional)
        * ``image_before_b64`` : str  (optional, base64-encoded PNG/JPEG)
        * ``image_after_b64`` : str  (optional, base64-encoded PNG/JPEG)
        * ``image_mime`` : str  (optional, e.g. ``"image/png"``)

    output_path : str or Path
        Where to write the HTML file.

    Returns
    -------
    str
        Absolute path to the generated file.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    body_sections = []
    chart_init_js: List[str] = []

    for idx, (name, data) in enumerate(results_dict.items()):
        sec_html, sec_js = _build_experiment_section(name, data, idx)
        body_sections.append(sec_html)
        chart_init_js.append(sec_js)

    full_html = _HTML_TEMPLATE.replace("{{BODY}}", "\n".join(body_sections))
    full_html = full_html.replace("{{CHART_INIT}}", "\n".join(chart_init_js))

    output_path.write_text(full_html, encoding="utf-8")
    return str(output_path.resolve())


# ===================================================================
# Internal: per-experiment section builder
# ===================================================================


def _build_experiment_section(name: str, data: dict, idx: int) -> tuple:
    """Return (html_string, js_init_string) for one experiment."""

    safe_name = html_lib.escape(name)
    canvas_sigma = f"chart_sigma_{idx}"
    canvas_cumul = f"chart_cumul_{idx}"
    canvas_fidelity = f"chart_fidelity_{idx}"

    # --- Sigma profile data -----------------------------------------------
    profile = data.get("sigma_profile", {})
    timesteps = json.dumps(profile.get("timesteps", []))
    sigma_vals = json.dumps(_sanitise_list(profile.get("sigma_per_step", [])))
    cumul_vals = json.dumps(_sanitise_list(profile.get("cumulative_sigma", [])))
    fidel_vals = json.dumps(_sanitise_list(profile.get("fidelity_bound", [])))

    # --- Summary cards ----------------------------------------------------
    cards = []
    total_sigma = data.get("total_sigma", None)
    if total_sigma is not None:
        cards.append(_card("Total Sigma", f"{total_sigma:.4f}", "blue"))

    speedup = data.get("speedup_ratio", None)
    if speedup is not None:
        cards.append(_card("Speedup", f"{speedup:.2f}x", "green"))

    ql = data.get("quality_loss_estimate", None)
    if ql is not None:
        cards.append(_card("Quality Loss Est.", f"{ql:.4f}", "orange"))

    psnr = data.get("psnr", None)
    if psnr is not None:
        psnr_str = f"{psnr:.2f} dB" if psnr != float("inf") else "inf"
        cards.append(_card("PSNR", psnr_str, "purple"))

    ssim = data.get("ssim", None)
    if ssim is not None:
        cards.append(_card("SSIM", f"{ssim:.4f}", "purple"))

    lpips_val = data.get("lpips", None)
    if lpips_val is not None:
        cards.append(_card("LPIPS", f"{lpips_val:.4f}", "purple"))

    n_active = data.get("num_active_steps", None)
    n_total = data.get("num_total_steps", None)
    if n_active is not None and n_total is not None:
        cards.append(_card("Steps", f"{n_active}/{n_total}", "text"))

    cards_html = "\n".join(cards) if cards else ""

    # --- Before / After images --------------------------------------------
    images_html = ""
    img_before = data.get("image_before_b64")
    img_after = data.get("image_after_b64")
    mime = data.get("image_mime", "image/png")
    if img_before or img_after:
        imgs = []
        if img_before:
            imgs.append(
                f'<div class="img-box"><div class="img-label">Before</div>'
                f'<img src="data:{mime};base64,{img_before}" /></div>'
            )
        if img_after:
            imgs.append(
                f'<div class="img-box"><div class="img-label">After (adaptive)</div>'
                f'<img src="data:{mime};base64,{img_after}" /></div>'
            )
        images_html = f'<div class="img-row">{"".join(imgs)}</div>'

    section_html = f"""
    <section class="experiment">
      <h2 class="exp-title">{safe_name}</h2>
      <div class="cards-row">{cards_html}</div>
      <div class="charts-row">
        <div class="chart-box"><canvas id="{canvas_sigma}"></canvas></div>
        <div class="chart-box"><canvas id="{canvas_cumul}"></canvas></div>
        <div class="chart-box"><canvas id="{canvas_fidelity}"></canvas></div>
      </div>
      {images_html}
    </section>
    """

    js = f"""
    // --- {safe_name} ---
    new Chart(document.getElementById('{canvas_sigma}'), {{
      type: 'line',
      data: {{
        labels: {timesteps},
        datasets: [{{
          label: 'Sigma_t',
          data: {sigma_vals},
          borderColor: '#2997ff',
          backgroundColor: 'rgba(41,151,255,0.1)',
          fill: true,
          tension: 0.3,
          pointRadius: 1,
        }}]
      }},
      options: {{ ...chartOpts, plugins: {{ ...chartOpts.plugins, title: {{ display:true, text:'Per-step Sigma', color:'#f5f5f7' }} }} }}
    }});
    new Chart(document.getElementById('{canvas_cumul}'), {{
      type: 'line',
      data: {{
        labels: {timesteps},
        datasets: [{{
          label: 'Cumulative Sigma',
          data: {cumul_vals},
          borderColor: '#bf5af2',
          backgroundColor: 'rgba(191,90,242,0.1)',
          fill: true,
          tension: 0.3,
          pointRadius: 1,
        }}]
      }},
      options: {{ ...chartOpts, plugins: {{ ...chartOpts.plugins, title: {{ display:true, text:'Cumulative Sigma', color:'#f5f5f7' }} }} }}
    }});
    new Chart(document.getElementById('{canvas_fidelity}'), {{
      type: 'line',
      data: {{
        labels: {timesteps},
        datasets: [{{
          label: 'Petz Fidelity Bound',
          data: {fidel_vals},
          borderColor: '#30d158',
          backgroundColor: 'rgba(48,209,88,0.1)',
          fill: true,
          tension: 0.3,
          pointRadius: 1,
        }}]
      }},
      options: {{ ...chartOpts, plugins: {{ ...chartOpts.plugins, title: {{ display:true, text:'Petz Fidelity Bound  F >= exp(-Sigma/2)', color:'#f5f5f7' }} }} }}
    }});
    """

    return section_html, js


# ===================================================================
# Helpers
# ===================================================================


def _card(label: str, value: str, colour: str) -> str:
    colour_map = {
        "blue": "#2997ff",
        "green": "#30d158",
        "orange": "#ff9f0a",
        "purple": "#bf5af2",
        "red": "#ff453a",
        "text": "#f5f5f7",
    }
    c = colour_map.get(colour, "#f5f5f7")
    return (
        f'<div class="card">'
        f'<div class="card-val" style="color:{c}">{html_lib.escape(value)}</div>'
        f'<div class="card-label">{html_lib.escape(label)}</div>'
        f"</div>"
    )


def _sanitise_list(lst: list) -> list:
    """Replace non-finite floats with None so JSON serialisation works."""
    out = []
    for v in lst:
        if isinstance(v, float) and (np.isnan(v) or np.isinf(v)):
            out.append(None)
        else:
            out.append(v)
    return out


# ===================================================================
# HTML template
# ===================================================================

_HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Sigma Diffusion Dashboard</title>
<script src="https://cdn.jsdelivr.net/npm/chart.js@4/dist/chart.umd.min.js"></script>
<style>
:root {
  --bg: #000;
  --bg-card: #111;
  --bg-card-hover: #161616;
  --text: #f5f5f7;
  --text-2: #a1a1a6;
  --text-3: #6e6e73;
  --blue: #2997ff;
  --purple: #bf5af2;
  --green: #30d158;
  --orange: #ff9f0a;
  --red: #ff453a;
  --pink: #ff375f;
  --grad: linear-gradient(90deg, #2997ff, #bf5af2, #ff375f);
}
* { margin:0; padding:0; box-sizing:border-box; }
html { scroll-behavior:smooth; -webkit-font-smoothing:antialiased; }
body {
  font-family: -apple-system, BlinkMacSystemFont, 'SF Pro Display', 'Helvetica Neue', Arial, sans-serif;
  background: var(--bg);
  color: var(--text);
  overflow-x: hidden;
  line-height: 1.5;
  padding: 0 0 80px 0;
}

/* ── Header ── */
.dash-header {
  text-align: center;
  padding: 60px 24px 40px;
  border-bottom: 1px solid rgba(255,255,255,0.06);
}
.dash-header h1 {
  font-size: clamp(28px, 4vw, 48px);
  font-weight: 700;
  letter-spacing: -0.035em;
  background: var(--grad);
  -webkit-background-clip: text;
  -webkit-text-fill-color: transparent;
  background-clip: text;
}
.dash-header p {
  font-size: 17px;
  color: var(--text-2);
  margin-top: 8px;
}

/* ── Experiment sections ── */
.experiment {
  max-width: 1200px;
  margin: 60px auto 0;
  padding: 0 24px;
}
.exp-title {
  font-size: clamp(22px, 3vw, 32px);
  font-weight: 700;
  letter-spacing: -0.02em;
  margin-bottom: 24px;
  padding-bottom: 12px;
  border-bottom: 1px solid rgba(255,255,255,0.08);
}

/* ── Cards ── */
.cards-row {
  display: flex;
  gap: 16px;
  flex-wrap: wrap;
  margin-bottom: 32px;
}
.card {
  background: var(--bg-card);
  border: 1px solid rgba(255,255,255,0.06);
  border-radius: 16px;
  padding: 20px 28px;
  min-width: 140px;
  transition: background 0.3s;
}
.card:hover { background: var(--bg-card-hover); }
.card-val {
  font-size: 28px;
  font-weight: 700;
  letter-spacing: -0.02em;
  line-height: 1.2;
}
.card-label {
  font-size: 13px;
  color: var(--text-3);
  margin-top: 4px;
  text-transform: uppercase;
  letter-spacing: 0.06em;
}

/* ── Charts ── */
.charts-row {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(320px, 1fr));
  gap: 20px;
  margin-bottom: 32px;
}
.chart-box {
  background: var(--bg-card);
  border: 1px solid rgba(255,255,255,0.06);
  border-radius: 16px;
  padding: 20px;
}
.chart-box canvas {
  width: 100% !important;
  max-height: 280px;
}

/* ── Images ── */
.img-row {
  display: flex;
  gap: 20px;
  flex-wrap: wrap;
  margin-top: 16px;
}
.img-box {
  background: var(--bg-card);
  border: 1px solid rgba(255,255,255,0.06);
  border-radius: 16px;
  padding: 16px;
  flex: 1;
  min-width: 250px;
}
.img-box img {
  width: 100%;
  border-radius: 8px;
  display: block;
}
.img-label {
  font-size: 13px;
  color: var(--text-3);
  text-transform: uppercase;
  letter-spacing: 0.06em;
  margin-bottom: 10px;
}

/* ── Footer ── */
.dash-footer {
  text-align: center;
  padding: 60px 24px 0;
  font-size: 13px;
  color: var(--text-3);
}

@media (max-width: 768px) {
  .charts-row { grid-template-columns: 1fr; }
  .cards-row { justify-content: center; }
}
</style>
</head>
<body>

<div class="dash-header">
  <h1>Sigma Diffusion Dashboard</h1>
  <p>Information-theoretic analysis &middot; Petz fidelity bound &middot; Adaptive scheduling</p>
</div>

{{BODY}}

<div class="dash-footer">
  Generated by <strong>sigma_diffusion</strong> &middot;
  Petz Recovery Unification &middot;
  Huang (2026)
</div>

<script>
const chartOpts = {
  responsive: true,
  maintainAspectRatio: false,
  scales: {
    x: {
      ticks: { color: '#6e6e73', maxTicksLimit: 12 },
      grid: { color: 'rgba(255,255,255,0.04)' },
      title: { display: true, text: 'Timestep', color: '#a1a1a6' }
    },
    y: {
      ticks: { color: '#6e6e73' },
      grid: { color: 'rgba(255,255,255,0.04)' },
    }
  },
  plugins: {
    legend: { labels: { color: '#a1a1a6' } },
    tooltip: {
      backgroundColor: '#222',
      titleColor: '#f5f5f7',
      bodyColor: '#a1a1a6',
      borderColor: 'rgba(255,255,255,0.1)',
      borderWidth: 1,
    },
  },
};

{{CHART_INIT}}
</script>
</body>
</html>
"""
