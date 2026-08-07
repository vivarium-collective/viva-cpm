"""Chemotactic-recruitment visualizations — interactive Plotly figures.

Data are baked in (from pbg_cpm_studies.chemotaxis.run) so the figures render
identically in the live and the published read-only dashboard. The headline
figure overlays the three conditions to show that recruitment requires BOTH a
secreted cue AND a competent chemotactic response.
"""
from __future__ import annotations

import json

from viva_superpowers.visualization import as_visualization

_INK = "#e5e7eb"; _MUTED = "#94a3b8"; _GRID = "#1e293b"; _BG = "#0b1220"; _CARD = "#0f172a"

# baked from workspace/chemotaxis_data/results/recruitment_summary.json (radius 15)
_MCS = [0, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375,
        400, 425, 450, 475, 500]
_SERIES = {
    "baseline":    [0.0, 0.0, 0.1667, 0.3333, 0.5, 0.5, 0.6667, 0.6667, 0.6667, 0.6667,
                    0.6667, 0.8333, 0.8333, 0.8333, 0.8333, 0.5, 0.3333, 0.5, 0.5, 0.6667, 0.6667],
    "inhibited":   [0.0] * 21,
    "adversarial": [0.0] * 21,
}
_DEFS = [
    ("baseline", "baseline · cue ✓ response ✓", "#38bdf8"),
    ("inhibited", "inhibited · cue ✓ response ✗ (intervention)", "#f43f5e"),
    ("adversarial", "adversarial · cue ✗ response ✓ (control)", "#a78bfa"),
]


def _recruitment_card(div):
    traces = []
    for key, label, color in _DEFS:
        traces.append({
            "x": _MCS, "y": _SERIES[key], "name": label, "type": "scatter",
            "mode": "lines+markers",
            "line": {"color": color, "width": 2.6, "shape": "spline", "smoothing": 0.6},
            "marker": {"color": color, "size": 6, "line": {"color": _CARD, "width": 1.5}},
            "hovertemplate": f"<b>{label}</b><br>MCS %{{x}}<br>recruited %{{y:.0%}}<extra></extra>",
        })
    layout = {
        "title": {"text": "<b>Recruitment requires both cue and response</b><br>"
                          f"<span style='font-size:12px;color:{_MUTED}'>fraction of responder cells reaching the source region (within 15 px)</span>",
                  "x": 0.02, "xanchor": "left", "font": {"size": 16, "color": _INK}},
        "paper_bgcolor": _CARD, "plot_bgcolor": _CARD,
        "font": {"family": "-apple-system, Segoe UI, Roboto, sans-serif", "color": _INK, "size": 12},
        "margin": {"l": 62, "r": 22, "t": 74, "b": 58},
        "xaxis": {"title": {"text": "Time (Monte-Carlo sweeps)"}, "gridcolor": _GRID,
                  "zeroline": False, "ticks": "outside", "tickcolor": _GRID, "linecolor": _GRID},
        "yaxis": {"title": {"text": "Recruitment index"}, "gridcolor": _GRID, "zeroline": False,
                  "range": [-0.03, 1.0], "tickformat": ".0%", "ticks": "outside",
                  "tickcolor": _GRID, "linecolor": _GRID},
        "hovermode": "x unified",
        "legend": {"orientation": "h", "y": -0.24, "x": 0, "font": {"size": 11.5}},
    }
    caption = ("<b>d=0.22 · decay=0.001 · λ=14 · rate=10</b> &nbsp;·&nbsp; Source cells secrete a "
               "diffusible cue; responders climb its gradient and accumulate at the source. "
               "Blocking the response (λ→0) <i>or</i> removing the cue (rate→0) each abolishes "
               "recruitment — the claim holds only when both are present.")
    return (
        f'<div style="background:{_BG};border:1px solid {_GRID};border-radius:14px;'
        f'padding:10px 12px 14px;max-width:780px;margin:0 auto;'
        f'font-family:-apple-system,Segoe UI,Roboto,sans-serif">'
        f'<div id="{div}" style="height:400px"></div>'
        f'<div style="color:{_MUTED};font-size:12.5px;line-height:1.55;padding:2px 8px 4px">{caption}</div>'
        f'<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>'
        f'<script>Plotly.newPlot("{div}",{json.dumps(traces)},{json.dumps(layout)},'
        f'{{responsive:true,displayModeBar:false}});</script></div>'
    )


@as_visualization(inputs={"mcs": "list[float]"}, name="ChemotaxisRecruitment",
                  demo={"mcs": _MCS})
def update_chemotaxis_recruitment(state):
    """Recruitment index over time — baseline vs inhibited vs adversarial"""
    return {"html": _recruitment_card("chemotaxis-recruitment")}
