"""Receptor-level recruitment visualizations — interactive Plotly figures.

Three figures make the *spatial chemotaxis mechanism* legible and tie it to the
cited receptor Kd (CXCL8–CXCR1, 2.9 nM, Nasser et al. 2009):

* ``ReceptorRecruitmentScene`` — the headline. An animated scene over MC sweeps:
  the diffusing chemokine field as a background gradient, the (mobile) source and
  the recruitment radius, and every responder cell at its centre of mass, coloured
  by activation state and sized by receptor occupancy θ. Responders chase the
  chemokine source and are recruited as their occupancy crosses the cited Kd.
* ``ReceptorOccupancyLaw`` — evidence → mechanism. The Hill occupancy law θ(c)
  with the cited Kd marked; each responder overlaid at (local concentration,
  occupancy). With hill=2 and an activation threshold of θ=0.5, the switch sits
  exactly at c = Kd, so a cell activates precisely when its local concentration
  exceeds the cited affinity.
* ``ReceptorRecruitmentCurve`` — the outcome. Recruitment index over time,
  baseline vs blocked, with 95% CI ribbons and the pass threshold marked.

All data are baked in (recruitment summaries from receptor_baseline.json /
receptor_blocked.json; spatial frames from receptor_spatial_baseline.json via
``_receptor_spatial_data``) so the figures render identically in the live and the
published read-only dashboard.
"""
from __future__ import annotations

import json
import math

from viva_superpowers.visualization import as_visualization

from ._receptor_spatial_data import SPATIAL_BASELINE

# ---------------------------------------------------------------------------
# shared palette (reused across all three figures)
# ---------------------------------------------------------------------------
_INK = "#e5e7eb"; _MUTED = "#94a3b8"; _GRID = "#1e293b"; _BG = "#0b1220"; _CARD = "#0f172a"
_NAIVE = "#7c8ba1"      # naive responder (muted slate)
_ACTIVE = "#f59e0b"     # activated responder (vivid amber)
_KD_COLOR = "#f43f5e"   # the cited Kd reference (rose)
_RADIUS = "#38bdf8"     # recruitment radius (sky)
_BASE = "#38bdf8"       # baseline series
_BLOCK = "#f43f5e"      # blocked series
_SOURCE = "#f8fafc"     # source marker
# chemokine field: dark -> teal -> bright (low -> high toward the source)
_FIELD_SCALE = [[0.0, "#0b1220"], [0.12, "#0c2733"], [0.35, "#0e7490"],
                [0.62, "#22d3ee"], [1.0, "#ecfeff"]]

_KD = 2.9
_CONC_SCALE = 0.02
_HILL = 2.0
_ACTIVATE = 0.5

_SPATIAL = json.loads(SPATIAL_BASELINE)


def _hex_to_rgba(hexcolor, alpha):
    h = hexcolor.lstrip("#")
    r, g, b = (int(h[i:i + 2], 16) for i in (0, 2, 4))
    return f"rgba({r},{g},{b},{alpha})"


def _plot(div, traces, layout, *, height, caption, frames=None, max_width=820):
    """Emit a self-contained Plotly card (matches the dashboard convention:
    CDN plotly + Plotly.newPlot in an inline <script>)."""
    cfg = "{responsive:true,displayModeBar:false}"
    if frames:
        script = (f'Plotly.newPlot("{div}",{json.dumps(traces)},{json.dumps(layout)},{cfg})'
                  f'.then(function(){{Plotly.addFrames("{div}",{json.dumps(frames)});}});')
    else:
        script = f'Plotly.newPlot("{div}",{json.dumps(traces)},{json.dumps(layout)},{cfg});'
    return (
        f'<div style="background:{_BG};border:1px solid {_GRID};border-radius:14px;'
        f'padding:10px 12px 14px;max-width:{max_width}px;margin:0 auto;'
        f'font-family:-apple-system,Segoe UI,Roboto,sans-serif">'
        f'<div id="{div}" style="height:{height}px"></div>'
        f'<div style="color:{_MUTED};font-size:12.5px;line-height:1.55;padding:2px 8px 4px">{caption}</div>'
        f'<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>'
        f'<script>{script}</script></div>'
    )


# ===========================================================================
# Figure 1 — ReceptorRecruitmentScene (headline, animated)
# ===========================================================================
_NX, _NY = _SPATIAL["lattice_dims"]
_RADIUS_PX = _SPATIAL["recruitment_radius"]
_FIELD_X = _SPATIAL["field_x"]
_FIELD_Y = _SPATIAL["field_y"]
# saturating colour cap so the gradient shape is visible rather than a single
# hotspot (the raw peak is ~2600 a.u.); the field VALUES are never altered.
_ZMAX = 800


def _circle(cx, cy, r, n=48):
    xs = [round(cx + r * math.cos(2 * math.pi * k / n), 2) for k in range(n + 1)]
    ys = [round(cy + r * math.sin(2 * math.pi * k / n), 2) for k in range(n + 1)]
    return xs, ys


def _scene_traces(frame, *, show_colorbar):
    field = frame["field"]
    sx, sy = frame["source_com"]
    heat = {
        "type": "heatmap", "z": field, "x": _FIELD_X, "y": _FIELD_Y,
        "zmin": 0, "zmax": _ZMAX, "colorscale": _FIELD_SCALE, "zsmooth": "best",
        "showscale": bool(show_colorbar), "hoverinfo": "skip", "showlegend": False,
        "colorbar": {"title": {"text": "chemokine<br>C (a.u.)", "font": {"size": 10, "color": _INK}},
                     "tickfont": {"color": _INK, "size": 9}, "thickness": 12, "len": 0.8,
                     "outlinecolor": _GRID},
    }
    cx, cy = _circle(sx, sy, _RADIUS_PX) if sx is not None else ([], [])
    radius = {
        "type": "scatter", "x": cx, "y": cy, "mode": "lines",
        "line": {"color": _RADIUS, "width": 1.6, "dash": "dot"},
        "name": f"recruitment radius ({int(_RADIUS_PX)} px)", "legendgroup": "radius",
        "hoverinfo": "skip",
    }
    source = {
        "type": "scatter", "x": [sx], "y": [sy], "mode": "markers",
        "marker": {"symbol": "star", "size": 17, "color": _SOURCE,
                   "line": {"color": _BG, "width": 1.4}},
        "name": "source (secretes C)", "legendgroup": "source",
        "hovertemplate": "source<br>x=%{x:.0f} y=%{y:.0f}<extra></extra>",
    }

    def _grp(activated):
        xs, ys, sizes, cd = [], [], [], []
        for c in frame["cells"]:
            if not c.get("alive") or bool(c.get("activated")) != activated:
                continue
            xs.append(c["x"]); ys.append(c["y"])
            sizes.append(round(11 + 24 * c["theta"], 1))
            cd.append([c["id"], round(c["theta"], 3),
                       "activated" if c["type"] == _SPATIAL["activated_type"] else "naive",
                       round(c["conc"], 2)])
        color = _ACTIVE if activated else _NAIVE
        return {
            "type": "scatter", "x": xs, "y": ys, "mode": "markers",
            "marker": {"size": sizes if sizes else 12, "color": color,
                       "line": {"color": "#0b1220", "width": 1.4}, "opacity": 0.95},
            "customdata": cd,
            "name": "activated responder" if activated else "naive responder",
            "legendgroup": "act" if activated else "naive",
            "hovertemplate": ("cell %{customdata[0]} · %{customdata[2]}<br>"
                              "θ=%{customdata[1]:.2f}<br>c=%{customdata[3]:.1f} nM"
                              "<extra></extra>"),
        }

    return [heat, radius, source, _grp(False), _grp(True)]


def _scene_card(div):
    frames_data = _SPATIAL["frames"]
    base = _scene_traces(frames_data[0], show_colorbar=True)
    frames = [{"name": str(fr["mcs"]),
               "data": _scene_traces(fr, show_colorbar=True)}
              for fr in frames_data]
    slider_steps = [{
        "label": str(fr["mcs"]), "method": "animate",
        "args": [[str(fr["mcs"])], {"mode": "immediate",
                                     "frame": {"duration": 0, "redraw": True},
                                     "transition": {"duration": 0}}],
    } for fr in frames_data]
    layout = {
        "title": {"text": "<b>Responders climb the chemokine gradient and are recruited "
                          "as receptor occupancy crosses the cited Kd</b><br>"
                          f"<span style='font-size:12px;color:{_MUTED}'>marker colour = "
                          "activation state · size ∝ receptor occupancy θ · background = "
                          f"secreted chemokine field · Kd = {_KD} nM (Nasser 2009) · baseline, seed "
                          f"{_SPATIAL['seed']}</span>",
                  "x": 0.02, "xanchor": "left", "font": {"size": 15, "color": _INK}},
        "paper_bgcolor": _CARD, "plot_bgcolor": _BG,
        "font": {"family": "-apple-system, Segoe UI, Roboto, sans-serif", "color": _INK, "size": 12},
        "margin": {"l": 54, "r": 22, "t": 92, "b": 96},
        "xaxis": {"title": {"text": "x (lattice units, px)"}, "range": [0, _NX],
                  "gridcolor": _GRID, "zeroline": False, "constrain": "domain",
                  "ticks": "outside", "tickcolor": _GRID, "linecolor": _GRID},
        "yaxis": {"title": {"text": "y (px)"}, "range": [0, _NY], "gridcolor": _GRID,
                  "zeroline": False, "scaleanchor": "x", "scaleratio": 1,
                  "ticks": "outside", "tickcolor": _GRID, "linecolor": _GRID},
        "legend": {"orientation": "h", "y": -0.30, "x": 0, "font": {"size": 11},
                   "bgcolor": "rgba(0,0,0,0)"},
        "updatemenus": [{
            "type": "buttons", "direction": "left", "showactive": False,
            "x": 0.02, "y": 1.02, "xanchor": "left", "yanchor": "bottom",
            "pad": {"t": 0, "r": 8}, "bgcolor": _CARD, "bordercolor": _GRID,
            "font": {"color": _INK, "size": 11},
            "buttons": [
                {"label": "▶ Play", "method": "animate",
                 "args": [None, {"frame": {"duration": 550, "redraw": True},
                                  "fromcurrent": True, "transition": {"duration": 0}}]},
                {"label": "❚❚ Pause", "method": "animate",
                 "args": [[None], {"mode": "immediate", "frame": {"duration": 0, "redraw": False},
                                    "transition": {"duration": 0}}]},
            ],
        }],
        "sliders": [{
            "active": 0, "x": 0.10, "y": 0, "len": 0.88,
            "pad": {"t": 34, "b": 4}, "currentvalue": {"prefix": "MC sweeps: ",
                                                        "font": {"color": _INK, "size": 12}},
            "font": {"color": _MUTED, "size": 10}, "bgcolor": _CARD,
            "bordercolor": _GRID, "steps": slider_steps,
        }],
    }
    caption = (f"<b>Kd={_KD} nM · conc_scale={_CONC_SCALE} · hill={_HILL}</b> &nbsp;·&nbsp; "
               "Each responder carries a per-cell Hill-occupancy receptor: it reads the local "
               "chemokine C, computes θ = c<sup>h</sup>/(Kd<sup>h</sup>+c<sup>h</sup>), and switches "
               "from naive (grey) to activated (amber, chemotaxing) once θ ≥ 0.5. The source is "
               "mobile — its chemokine blob wanders and activated responders chase it, staying inside "
               "the recruitment radius, while naive or squeezed-out cells fall behind. Field colour is "
               f"shown on a saturating 0–{_ZMAX} a.u. scale to reveal the gradient (peak ≈ 2.6k a.u.).")
    return _plot(div, base, layout, height=470, caption=caption, frames=frames, max_width=860)


@as_visualization(inputs={"mcs": "list[float]"}, name="ReceptorRecruitmentScene",
                  demo={"mcs": [fr["mcs"] for fr in _SPATIAL["frames"]]})
def update_receptor_recruitment_scene(state):
    """Animated spatial scene — responders climb the chemokine gradient and are recruited (baseline)"""
    return {"html": _scene_card("receptor-scene")}


# ===========================================================================
# Figure 2 — ReceptorOccupancyLaw (evidence -> mechanism)
# ===========================================================================
def _occupancy(c):
    if c <= 0:
        return 0.0
    ch = c ** _HILL
    return ch / (_KD ** _HILL + ch)


def _occupancy_card(div):
    # concentrations seen across the run, to set a sensible x-range
    concs = [c["conc"] for fr in _SPATIAL["frames"] for c in fr["cells"] if c.get("alive")]
    cmax = max(30.0, max(concs) * 1.05) if concs else 30.0
    n = 240
    xs = [cmax * k / (n - 1) for k in range(n)]
    ys = [_occupancy(x) for x in xs]

    curve = {
        "type": "scatter", "x": [round(x, 3) for x in xs], "y": [round(y, 4) for y in ys],
        "mode": "lines", "line": {"color": "#38bdf8", "width": 3, "shape": "spline"},
        "name": "Hill occupancy law", "legendgroup": "law",
        "hovertemplate": "c=%{x:.2f} nM<br>θ=%{y:.2f}<extra></extra>",
    }

    # faint trajectory: every (conc, θ) each responder actually passed through,
    # coloured by its activation state at that frame -> cells sweep up the curve
    # from naive(low-left) to activated(up-right), crossing at c = Kd.
    traj_n = {"type": "scatter", "x": [], "y": [], "mode": "markers", "customdata": [],
              "marker": {"size": 6, "color": _NAIVE, "opacity": 0.30,
                         "line": {"width": 0}},
              "name": "responder (naive, over time)", "legendgroup": "trajn",
              "hovertemplate": "cell %{customdata[0]} · naive<br>c=%{x:.1f} nM<br>θ=%{y:.2f}<extra></extra>"}
    traj_a = {"type": "scatter", "x": [], "y": [], "mode": "markers", "customdata": [],
              "marker": {"size": 6, "color": _ACTIVE, "opacity": 0.30,
                         "line": {"width": 0}},
              "name": "responder (activated, over time)", "legendgroup": "traja",
              "hovertemplate": "cell %{customdata[0]} · activated<br>c=%{x:.1f} nM<br>θ=%{y:.2f}<extra></extra>"}
    for fr in _SPATIAL["frames"]:
        for c in fr["cells"]:
            if not c.get("alive"):
                continue
            t = traj_a if c.get("activated") else traj_n
            t["x"].append(round(c["conc"], 3)); t["y"].append(round(c["theta"], 4))
            t["customdata"].append([c["id"]])

    # bold emphasis: the FINAL baseline frame's live cells
    final = _SPATIAL["frames"][-1]
    fin_n = {"type": "scatter", "x": [], "y": [], "mode": "markers", "customdata": [],
             "marker": {"size": 15, "color": _NAIVE, "line": {"color": _INK, "width": 1.6},
                        "symbol": "circle"},
             "name": "naive @ t=final", "legendgroup": "finn",
             "hovertemplate": "cell %{customdata[0]} · naive (final)<br>c=%{x:.1f} nM<br>θ=%{y:.2f}<extra></extra>"}
    fin_a = {"type": "scatter", "x": [], "y": [], "mode": "markers", "customdata": [],
             "marker": {"size": 16, "color": _ACTIVE, "line": {"color": _INK, "width": 1.6},
                        "symbol": "circle"},
             "name": "activated @ t=final", "legendgroup": "fina",
             "hovertemplate": "cell %{customdata[0]} · activated (final)<br>c=%{x:.1f} nM<br>θ=%{y:.2f}<extra></extra>"}
    for c in final["cells"]:
        if not c.get("alive"):
            continue
        t = fin_a if c.get("activated") else fin_n
        t["x"].append(round(c["conc"], 3)); t["y"].append(round(c["theta"], 4))
        t["customdata"].append([c["id"]])

    theta_at_kd = _occupancy(_KD)
    layout = {
        "title": {"text": "<b>A cell activates precisely when its local chemokine exceeds the "
                          "cited receptor Kd</b><br>"
                          f"<span style='font-size:12px;color:{_MUTED}'>Hill occupancy θ(c) with the "
                          f"CXCL8–CXCR1 Kd marked · hill={_HILL}, threshold θ=0.5 ⇒ switch sits exactly "
                          f"at c = Kd · baseline responders overlaid</span>",
                  "x": 0.02, "xanchor": "left", "font": {"size": 15, "color": _INK}},
        "paper_bgcolor": _CARD, "plot_bgcolor": _CARD,
        "font": {"family": "-apple-system, Segoe UI, Roboto, sans-serif", "color": _INK, "size": 12},
        "margin": {"l": 64, "r": 22, "t": 88, "b": 62},
        "xaxis": {"title": {"text": "local chemokine concentration c (nM)"}, "gridcolor": _GRID,
                  "zeroline": False, "range": [0, cmax], "ticks": "outside",
                  "tickcolor": _GRID, "linecolor": _GRID},
        "yaxis": {"title": {"text": "receptor occupancy θ"}, "gridcolor": _GRID, "zeroline": False,
                  "range": [0, 1.03], "ticks": "outside", "tickcolor": _GRID, "linecolor": _GRID},
        "hovermode": "closest",
        "legend": {"orientation": "h", "y": -0.22, "x": 0, "font": {"size": 10.5},
                   "bgcolor": "rgba(0,0,0,0)"},
        "shapes": [
            {"type": "line", "x0": _KD, "x1": _KD, "y0": 0, "y1": 1.03, "yref": "y",
             "line": {"color": _KD_COLOR, "width": 2, "dash": "dash"}},
            {"type": "line", "x0": 0, "x1": cmax, "y0": _ACTIVATE, "y1": _ACTIVATE,
             "line": {"color": _MUTED, "width": 1.4, "dash": "dot"}},
        ],
        "annotations": [
            {"x": _KD, "y": 1.02, "xanchor": "left", "yanchor": "top",
             "text": f"<b>Kd = {_KD} nM</b><br>Nasser 2009", "showarrow": False,
             "font": {"color": _KD_COLOR, "size": 11}, "align": "left", "xshift": 6},
            {"x": cmax, "y": _ACTIVATE, "xanchor": "right", "yanchor": "bottom",
             "text": "activation threshold θ = 0.5", "showarrow": False,
             "font": {"color": _MUTED, "size": 10.5}, "yshift": 3},
            {"x": _KD, "y": theta_at_kd, "text": "switch point", "showarrow": True,
             "arrowhead": 2, "arrowcolor": _KD_COLOR, "ax": 42, "ay": 30,
             "font": {"color": _INK, "size": 10.5}},
        ],
    }
    caption = ("The receptor is not a free-fit knob: θ(c) uses the CXCL8–CXCR1 affinity reported by "
               "<i>Nasser et al. 2009</i> (Kd = 2.9 nM). Faint points trace every responder's "
               "(concentration, occupancy) over time; bold points are the final-frame survivors. With "
               "hill=2 and a θ=0.5 switch, the vertical Kd line and the horizontal threshold meet on "
               "the curve — cells to the right of Kd cross into the activated (amber) state, and those "
               "are the cells that chemotax and are recruited.")
    return _plot(div, [curve, traj_n, traj_a, fin_n, fin_a], layout,
                 height=430, caption=caption)


@as_visualization(inputs={"mcs": "list[float]"}, name="ReceptorOccupancyLaw",
                  demo={"mcs": [_SPATIAL["frames"][-1]["mcs"]]})
def update_receptor_occupancy_law(state):
    """Hill occupancy law with the cited Kd, responders overlaid — evidence to mechanism (baseline)"""
    return {"html": _occupancy_card("receptor-occupancy-law")}


# ===========================================================================
# Figure 3 — ReceptorRecruitmentCurve (polished baseline vs blocked)
# ===========================================================================
# baked from workspace/chemotaxis_data/results/receptor_baseline.json /
# receptor_blocked.json (5 seeds [17, 29, 43, 61, 89]; radius 15)
_MCS = [0, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375,
        400, 425, 450, 475, 500]
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
    ("baseline", "receptor baseline · cue ✓ response ✓", _BASE),
    ("blocked", "receptor blocked · cue ✓ response ✗ (intervention)", _BLOCK),
]
_PASS = 0.5  # primary-test pass threshold (recruitment index > 0.5)


def _recruitment_card(div):
    traces = []
    for key, label, color in _DEFS:
        y = _MEAN[key]
        lo = [max(0.0, c[0]) for c in _CI[key]]   # lower CI clamp is documented (index >= 0)
        hi = [c[1] for c in _CI[key]]
        rgba = _hex_to_rgba(color, 0.16)
        traces.append({"x": _MCS, "y": hi, "type": "scatter", "mode": "lines",
                       "line": {"width": 0, "color": color}, "hoverinfo": "skip",
                       "showlegend": False, "legendgroup": key})
        traces.append({"x": _MCS, "y": lo, "type": "scatter", "mode": "lines",
                       "line": {"width": 0, "color": color}, "fill": "tonexty",
                       "fillcolor": rgba, "hoverinfo": "skip", "showlegend": False,
                       "legendgroup": key})
        traces.append({
            "x": _MCS, "y": y, "name": label, "type": "scatter",
            "mode": "lines+markers", "legendgroup": key,
            "line": {"color": color, "width": 2.6, "shape": "spline", "smoothing": 0.6},
            "marker": {"color": color, "size": 6, "line": {"color": _CARD, "width": 1.5}},
            "hovertemplate": f"<b>{label}</b><br>MCS %{{x}}<br>recruited %{{y:.0%}}<extra></extra>",
        })
    layout = {
        "title": {"text": "<b>The receptor-gated response is what drives recruitment — "
                          "blocking it abolishes the effect</b><br>"
                          f"<span style='font-size:12px;color:{_MUTED}'>fraction of responders within "
                          f"{int(_RADIUS_PX)} px of the source · 95% CI over 5 seeds · activation gated "
                          f"by receptor Kd = {_KD} nM (Nasser 2009)</span>",
                  "x": 0.02, "xanchor": "left", "font": {"size": 15, "color": _INK}},
        "paper_bgcolor": _CARD, "plot_bgcolor": _CARD,
        "font": {"family": "-apple-system, Segoe UI, Roboto, sans-serif", "color": _INK, "size": 12},
        "margin": {"l": 62, "r": 22, "t": 84, "b": 60},
        "xaxis": {"title": {"text": "Time (Monte-Carlo sweeps)"}, "gridcolor": _GRID,
                  "zeroline": False, "ticks": "outside", "tickcolor": _GRID, "linecolor": _GRID},
        "yaxis": {"title": {"text": "Recruitment index"}, "gridcolor": _GRID, "zeroline": False,
                  "range": [-0.1, 1.02], "tickformat": ".0%", "ticks": "outside",
                  "tickcolor": _GRID, "linecolor": _GRID},
        "hovermode": "x unified",
        "legend": {"orientation": "h", "y": -0.24, "x": 0, "font": {"size": 11.5},
                   "bgcolor": "rgba(0,0,0,0)"},
        "shapes": [{"type": "line", "x0": 0, "x1": _MCS[-1], "y0": _PASS, "y1": _PASS,
                    "line": {"color": _MUTED, "width": 1.3, "dash": "dot"}}],
        "annotations": [{"x": 0, "y": _PASS, "xanchor": "left", "yanchor": "bottom",
                         "text": "pass threshold (0.5)", "showarrow": False,
                         "font": {"color": _MUTED, "size": 10.5}, "xshift": 4, "yshift": 2}],
    }
    caption = (f"<b>kd={_KD} nM · conc_scale={_CONC_SCALE} · hill={_HILL}</b> &nbsp;·&nbsp; "
               "Each responder carries a per-cell Hill-occupancy receptor gated by a cited CXCL8–CXCR1 "
               "Kd; only once occupancy crosses threshold does it switch from naive to chemotaxing. "
               "Forcing the activated-type response off (<i>blocked</i>) abolishes recruitment even "
               "though the cue is still present and receptors still bind it — isolating the response, "
               "not the cue, as the causal step. Ribbons are 95% CIs across 5 seeds.")
    return _plot(div, traces, layout, height=430, caption=caption)


@as_visualization(inputs={"mcs": "list[float]"}, name="ReceptorRecruitmentCurve",
                  demo={"mcs": _MCS})
def update_receptor_recruitment_curve(state):
    """Recruitment index over time — receptor baseline vs blocked, CI ribbons + pass threshold"""
    return {"html": _recruitment_card("receptor-recruitment-curve")}


# Plain zero-arg accessors matching the registered viz names, so the module can
# be exercised directly (dashboard preview harness / tests) without going
# through the Step lifecycle.
def ReceptorRecruitmentScene():
    return _scene_card("receptor-scene")


def ReceptorOccupancyLaw():
    return _occupancy_card("receptor-occupancy-law")


def ReceptorRecruitmentCurve():
    return _recruitment_card("receptor-recruitment-curve")
