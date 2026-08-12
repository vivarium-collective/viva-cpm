"""Receptor-level recruitment visualizations — interactive Plotly figures.

Data are baked in (from workspace/chemotaxis_data/results/receptor_baseline.json
and receptor_blocked.json) so the figures render identically in the live and
the published read-only dashboard. The headline figure overlays the two
receptor-realization conditions (baseline vs blocked) with 95% CI ribbons;
there is no receptor no-cue arm in this increment, so only two series are
shown here (unlike the three-condition ChemotaxisRecruitment figure). The
second figure exposes the activation FRONT per responder cell over time — a
readout the scalar-lambda phenomenological model has no analogue for.
"""
from __future__ import annotations

import json

from viva_superpowers.visualization import as_visualization

_INK = "#e5e7eb"; _MUTED = "#94a3b8"; _GRID = "#1e293b"; _BG = "#0b1220"; _CARD = "#0f172a"

# baked from workspace/chemotaxis_data/results/receptor_baseline.json and
# receptor_blocked.json (5 seeds [17, 29, 43, 61, 89]; radius 15)
_MCS = [0, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375,
        400, 425, 450, 475, 500]

_KD = 2.9
_CONC_SCALE = 0.02
_HILL = 2.0

_MEAN = {
    "baseline": [0.0, 0.0, 0.1, 0.2667, 0.3667, 0.4667, 0.5, 0.7333, 0.7333, 0.7667,
                 0.7333, 0.6333, 0.6667, 0.6333, 0.6667, 0.7667, 0.6, 0.6, 0.6, 0.6, 0.5667],
    "blocked": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0333, 0.0667, 0.0667, 0.0667, 0.0667, 0.0667, 0.0667],
}
_CI = {
    "baseline": [[0.0, 0.0], [0.0, 0.0], [0.02, 0.18], [0.136, 0.3973], [0.1266, 0.6067],
                 [0.2762, 0.6571], [0.269, 0.731], [0.5668, 0.8999], [0.5118, 0.9549],
                 [0.5707, 0.9627], [0.5118, 0.9549], [0.4429, 0.8238], [0.4877, 0.8456],
                 [0.4429, 0.8238], [0.4877, 0.8456], [0.5707, 0.9627], [0.3784, 0.8216],
                 [0.3784, 0.8216], [0.3784, 0.8216], [0.3784, 0.8216], [0.3451, 0.7882]],
    "blocked": [[0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0],
                [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0],
                [0.0, 0.0], [0.0, 0.0], [-0.032, 0.0987], [-0.0133, 0.1467],
                [-0.0133, 0.1467], [-0.0133, 0.1467], [-0.0133, 0.1467], [-0.0133, 0.1467],
                [-0.0133, 0.1467]],
}
_DEFS = [
    ("baseline", "receptor baseline · cue ✓ response ✓", "#38bdf8"),
    ("blocked", "receptor blocked · cue ✓ response ✗ (intervention)", "#f43f5e"),
]

# baked from receptor_baseline.json activation_by_cell (responders 2..7,
# per-tick activation fraction across the 5 seeds)
_ACTIVATION_CELLS = [2, 3, 4, 5, 6, 7]
_ACTIVATION_BY_CELL = {
    2: [0.0, 0.0, 0.6, 1.0, 1.0, 1.0, 1.0, 1.0, 0.6, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4,
        0.4, 0.4, 0.4, 0.4, 0.4],
    3: [0.0, 0.0, 0.2, 0.4, 0.6, 0.8, 0.8, 1.0, 1.0, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8,
        0.8, 0.8, 0.8, 0.8, 0.6],
    4: [0.0, 0.0, 0.0, 0.2, 0.2, 0.2, 0.2, 0.6, 0.6, 0.6, 0.6, 0.4, 0.4, 0.4, 0.4, 0.4,
        0.4, 0.4, 0.4, 0.4, 0.4],
    5: [0.0, 0.0, 0.0, 0.0, 0.2, 0.4, 0.4, 0.8, 0.8, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6,
        0.6, 0.6, 0.6, 0.6, 0.6],
    6: [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.8, 0.6, 0.6, 0.6, 0.6, 0.4, 0.4, 0.4, 0.4,
        0.2, 0.2, 0.2, 0.2, 0.2],
    7: [0.0, 0.0, 0.0, 0.0, 0.2, 0.4, 0.4, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6,
        0.6, 0.4, 0.4, 0.4, 0.4],
}


def _hex_to_rgba(hexcolor, alpha):
    h = hexcolor.lstrip("#")
    r, g, b = (int(h[i:i + 2], 16) for i in (0, 2, 4))
    return f"rgba({r},{g},{b},{alpha})"


def _recruitment_card(div):
    traces = []
    for key, label, color in _DEFS:
        y = _MEAN[key]
        lo = [c[0] for c in _CI[key]]
        hi = [c[1] for c in _CI[key]]
        rgba = _hex_to_rgba(color, 0.16)
        # invisible upper bound
        traces.append({
            "x": _MCS, "y": hi, "type": "scatter", "mode": "lines",
            "line": {"width": 0, "color": color}, "hoverinfo": "skip",
            "showlegend": False, "legendgroup": key,
        })
        # invisible lower bound, filled up to the previous (upper) trace -> CI ribbon
        traces.append({
            "x": _MCS, "y": lo, "type": "scatter", "mode": "lines",
            "line": {"width": 0, "color": color}, "fill": "tonexty",
            "fillcolor": rgba, "hoverinfo": "skip", "showlegend": False,
            "legendgroup": key,
        })
        # visible mean line on top
        traces.append({
            "x": _MCS, "y": y, "name": label, "type": "scatter",
            "mode": "lines+markers", "legendgroup": key,
            "line": {"color": color, "width": 2.6, "shape": "spline", "smoothing": 0.6},
            "marker": {"color": color, "size": 6, "line": {"color": _CARD, "width": 1.5}},
            "hovertemplate": f"<b>{label}</b><br>MCS %{{x}}<br>recruited %{{y:.0%}}<extra></extra>",
        })
    layout = {
        "title": {"text": "<b>Receptor-gated recruitment — baseline vs blocked</b><br>"
                          f"<span style='font-size:12px;color:{_MUTED}'>fraction of responder "
                          f"cells within 15 px of the source · 95% CI, 5 seeds · activation gated "
                          f"by receptor Kd = {_KD} nM (Nasser et al. 2009)</span>",
                  "x": 0.02, "xanchor": "left", "font": {"size": 16, "color": _INK}},
        "paper_bgcolor": _CARD, "plot_bgcolor": _CARD,
        "font": {"family": "-apple-system, Segoe UI, Roboto, sans-serif", "color": _INK, "size": 12},
        "margin": {"l": 62, "r": 22, "t": 84, "b": 58},
        "xaxis": {"title": {"text": "Time (Monte-Carlo sweeps)"}, "gridcolor": _GRID,
                  "zeroline": False, "ticks": "outside", "tickcolor": _GRID, "linecolor": _GRID},
        "yaxis": {"title": {"text": "Recruitment index"}, "gridcolor": _GRID, "zeroline": False,
                  "range": [-0.1, 1.02], "tickformat": ".0%", "ticks": "outside",
                  "tickcolor": _GRID, "linecolor": _GRID},
        "hovermode": "x unified",
        "legend": {"orientation": "h", "y": -0.24, "x": 0, "font": {"size": 11.5}},
    }
    caption = (f"<b>kd={_KD} nM · conc_scale={_CONC_SCALE} · hill={_HILL}</b> &nbsp;·&nbsp; "
               "Each responder carries a per-cell Hill-occupancy receptor (ReceptorSubcell) "
               "gated by a cited CXCL8-CXCR1 Kd; only once occupancy crosses threshold does it "
               "switch from naive to chemotaxing. Forcing the activated-type response off "
               "(<i>blocked</i>) abolishes recruitment even though the cue is still present and "
               "receptors still bind it — no receptor-blank / no-cue third arm exists in this "
               "increment.")
    return (
        f'<div style="background:{_BG};border:1px solid {_GRID};border-radius:14px;'
        f'padding:10px 12px 14px;max-width:780px;margin:0 auto;'
        f'font-family:-apple-system,Segoe UI,Roboto,sans-serif">'
        f'<div id="{div}" style="height:420px"></div>'
        f'<div style="color:{_MUTED};font-size:12.5px;line-height:1.55;padding:2px 8px 4px">{caption}</div>'
        f'<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>'
        f'<script>Plotly.newPlot("{div}",{json.dumps(traces)},{json.dumps(layout)},'
        f'{{responsive:true,displayModeBar:false}});</script></div>'
    )


@as_visualization(inputs={"mcs": "list[float]"}, name="ReceptorRecruitment", demo={"mcs": _MCS})
def update_receptor_recruitment(state):
    """Recruitment index over time — receptor baseline vs receptor blocked, with CI ribbons"""
    return {"html": _recruitment_card("receptor-recruitment")}


def _activation_map(div):
    z = [_ACTIVATION_BY_CELL[cid] for cid in _ACTIVATION_CELLS]
    y_labels = [f"cell {cid}" for cid in _ACTIVATION_CELLS]
    trace = {
        "type": "heatmap", "z": z, "x": _MCS, "y": y_labels,
        "colorscale": [[0.0, "#0f172a"], [0.35, "#1d4ed8"], [0.7, "#38bdf8"], [1.0, "#facc15"]],
        "zmin": 0.0, "zmax": 1.0,
        "colorbar": {"title": {"text": "activation<br>fraction", "font": {"size": 11, "color": _INK}},
                     "tickfont": {"color": _INK, "size": 10}, "outlinecolor": _GRID,
                     "thickness": 14, "len": 0.85},
        "hovertemplate": "cell %{y}<br>MCS %{x}<br>activated %{z:.0%}<extra></extra>",
        "xgap": 1, "ygap": 3,
    }
    layout = {
        "title": {"text": "<b>Receptor activation front</b><br>"
                          f"<span style='font-size:12px;color:{_MUTED}'>fraction of the 5 seeds in "
                          "which each responder's Hill-occupancy has crossed the activation "
                          f"threshold (kd = {_KD} nM, Nasser et al. 2009) · receptor baseline</span>",
                  "x": 0.02, "xanchor": "left", "font": {"size": 16, "color": _INK}},
        "paper_bgcolor": _CARD, "plot_bgcolor": _CARD,
        "font": {"family": "-apple-system, Segoe UI, Roboto, sans-serif", "color": _INK, "size": 12},
        "margin": {"l": 76, "r": 22, "t": 84, "b": 58},
        "xaxis": {"title": {"text": "Time (Monte-Carlo sweeps)"}, "gridcolor": _GRID,
                  "zeroline": False, "ticks": "outside", "tickcolor": _GRID, "linecolor": _GRID},
        "yaxis": {"title": {"text": "Responder cell (id)"}, "gridcolor": _GRID, "zeroline": False,
                  "ticks": "outside", "tickcolor": _GRID, "linecolor": _GRID, "autorange": "reversed"},
    }
    caption = ("Each row is one of the six responder cells (ids 2..7); each column a sampled MC "
               "sweep. Color is the fraction of the 5 seeds in which that responder had already "
               "crossed the Hill-occupancy activation threshold at that tick — cells nearer the "
               "source activate earlier, producing a visible activation front that a scalar-&#955; "
               "phenomenological model has no analogue for.")
    return (
        f'<div style="background:{_BG};border:1px solid {_GRID};border-radius:14px;'
        f'padding:10px 12px 14px;max-width:780px;margin:0 auto;'
        f'font-family:-apple-system,Segoe UI,Roboto,sans-serif">'
        f'<div id="{div}" style="height:360px"></div>'
        f'<div style="color:{_MUTED};font-size:12.5px;line-height:1.55;padding:2px 8px 4px">{caption}</div>'
        f'<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>'
        f'<script>Plotly.newPlot("{div}",[{json.dumps(trace)}],{json.dumps(layout)},'
        f'{{responsive:true,displayModeBar:false}});</script></div>'
    )


@as_visualization(inputs={"mcs": "list[float]"}, name="ReceptorActivationMap", demo={"mcs": _MCS})
def update_receptor_activation_map(state):
    """Per-responder receptor activation over time — the activation front (receptor baseline)"""
    return {"html": _activation_map("receptor-activation-map")}


# Plain zero-arg accessors matching the registered viz names, so the module can
# be exercised directly (dashboard preview harness / tests) without going
# through the Step lifecycle (config/core construction + accumulate/render).
def ReceptorRecruitment():
    return _recruitment_card("receptor-recruitment")


def ReceptorActivationMap():
    return _activation_map("receptor-activation-map")
