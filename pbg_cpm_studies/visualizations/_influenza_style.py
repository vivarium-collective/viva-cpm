"""Shared visual system for the influenza-sego2022 study visualizations.

ONE design language applied across every Sego study (and forward-compatible with
the immune-arm increments): a validated color palette, a consistent card shell
(title + increment badge + caption), and shared Plotly layout/plot builders so
every figure reads as one system.

Colors are the data-viz reference palette (validated for CVD + contrast on the
light chart surface). Cell-state hues are semantic and fixed:

    Healthy  aqua    #1baf7a      (ships)      Macrophage  blue    #2a78d6  (future)
    Infected orange  #eb6834      (ships)      NK cell     violet  #4a3aa7  (future)
    Dead     slate   #6b7280      (future)     APC/CD8+T   magenta #e87ba4  (future)

Healthy+Infected (the states that appear in Increments 1-2) pass every CVD and
contrast gate as a pair. The full state set relies on the per-mark hover labels
that spatial maps always carry (secondary encoding) once the immune arms land.
"""
from __future__ import annotations
import json

# ── chrome / ink (light chart surface) ──────────────────────────────────────
SURFACE = "#fcfcfb"
PLANE = "#f6f6f4"
INK = "#0b0b0b"
SECONDARY = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
BASELINE = "#c3c2b7"
BORDER = "rgba(11,11,11,0.10)"
FONT = "system-ui, -apple-system, 'Segoe UI', Roboto, sans-serif"

# ── data hues (validated categorical + sequential) ──────────────────────────
CELL_STATES = {
    1: ("Healthy epithelial", "#1baf7a"),
    2: ("Infected epithelial", "#eb6834"),
    3: ("Dead epithelial", "#6b7280"),
    4: ("Macrophage", "#2a78d6"),
    5: ("NK cell", "#4a3aa7"),
    6: ("APC / CD8⁺ T", "#e87ba4"),
}
HEALTHY, INFECTED, DEAD = CELL_STATES[1][1], CELL_STATES[2][1], CELL_STATES[3][1]

# blue sequential ramp for the virus field (light -> dark = low -> high)
VIRUS_SCALE = [
    [0.0, "#fcfcfb"], [0.12, "#cde2fb"], [0.3, "#9ec5f4"], [0.5, "#6da7ec"],
    [0.7, "#3987e5"], [0.85, "#256abf"], [1.0, "#0d366b"],
]
# status hues (fixed, never themed)
GOOD, WARNING, CRITICAL = "#0ca30c", "#fab219", "#d03b3b"

# increment badge tints
_BADGE_TINT = {
    0: ("#eef2ff", "#4338ca"), 1: ("#ecfdf5", "#047857"),
    2: ("#fff7ed", "#c2410c"), 3: ("#eff6ff", "#1d4ed8"),
}


def base_layout(xtitle="", ytitle="", *, logx=False, height=380, showlegend=True,
                extra=None):
    """Shared Plotly layout: transparent plot, recessive axes, unified fonts."""
    layout = {
        "paper_bgcolor": "rgba(0,0,0,0)", "plot_bgcolor": "rgba(0,0,0,0)",
        "font": {"family": FONT, "color": INK, "size": 12},
        "margin": {"l": 62, "r": 24, "t": 14, "b": 52},
        "height": height,
        "xaxis": {"title": {"text": xtitle, "font": {"size": 12, "color": SECONDARY}},
                  "type": "log" if logx else "linear", "gridcolor": GRID,
                  "zeroline": False, "ticks": "outside", "tickcolor": GRID,
                  "linecolor": BASELINE, "tickfont": {"color": MUTED}},
        "yaxis": {"title": {"text": ytitle, "font": {"size": 12, "color": SECONDARY}},
                  "gridcolor": GRID, "zeroline": False, "ticks": "outside",
                  "tickcolor": GRID, "linecolor": BASELINE, "tickfont": {"color": MUTED}},
        "hovermode": "x unified",
        "legend": {"orientation": "h", "y": -0.22, "x": 0, "font": {"size": 12},
                   "bgcolor": "rgba(0,0,0,0)"},
        "showlegend": showlegend,
    }
    if extra:
        layout.update(extra)
    return layout


def band_trace(x, lo, hi, color, name):
    """A soft lo/hi acceptance band (filled ribbon) behind a reference line."""
    rgba = _hex_rgba(color, 0.14)
    return [
        {"x": x, "y": hi, "type": "scatter", "mode": "lines", "line": {"width": 0},
         "hoverinfo": "skip", "showlegend": False, "name": name + " hi"},
        {"x": x, "y": lo, "type": "scatter", "mode": "lines", "line": {"width": 0},
         "fill": "tonexty", "fillcolor": rgba, "hoverinfo": "skip",
         "showlegend": False, "name": name + " band"},
    ]


def line_trace(x, y, color, name, *, dash=None, width=2.4, markers=True):
    tr = {
        "x": x, "y": y, "name": name, "type": "scatter",
        "mode": "lines+markers" if markers else "lines",
        "line": {"color": color, "width": width, "shape": "spline", "smoothing": 0.5},
        "hovertemplate": f"<b>{name}</b><br>%{{x}}<br>%{{y:.3g}}<extra></extra>",
    }
    if dash:
        tr["line"]["dash"] = dash
    if markers:
        tr["marker"] = {"color": color, "size": 6, "line": {"color": SURFACE, "width": 1.4}}
    return tr


def plot(div_id, traces, layout, *, frames=None, config=None):
    """Emit a Plotly figure as a self-contained div + inline script."""
    cfg = {"responsive": True, "displayModeBar": False}
    if config:
        cfg.update(config)
    if frames is not None:
        fig = {"data": traces, "layout": layout, "frames": frames}
        script = (f'Plotly.newPlot("{div_id}",{json.dumps(fig)},{json.dumps(cfg)});')
    else:
        script = (f'Plotly.newPlot("{div_id}",{json.dumps(traces)},'
                  f'{json.dumps(layout)},{json.dumps(cfg)});')
    return (f'<div id="{div_id}"></div>'
            f'<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>'
            f'<script>{script}</script>')


def card(title, subtitle, body_html, caption, *, increment=None, kpis=None,
         max_width=820):
    """The unified card shell: header (title + increment badge), optional KPI
    row, the figure body, and a source caption. One look for every study."""
    badge = ""
    if increment is not None:
        bg, fg = _BADGE_TINT.get(increment, ("#eef2ff", "#4338ca"))
        badge = (f'<span style="background:{bg};color:{fg};font-size:11px;'
                 f'font-weight:600;padding:3px 9px;border-radius:999px;'
                 f'white-space:nowrap">Increment {increment}</span>')
    kpi_row = ""
    if kpis:
        cells = "".join(
            f'<div style="flex:1;min-width:120px"><div style="font-size:20px;'
            f'font-weight:650;color:{INK};letter-spacing:-.01em">{v}</div>'
            f'<div style="font-size:11.5px;color:{MUTED};margin-top:1px">{k}</div></div>'
            for k, v in kpis)
        kpi_row = (f'<div style="display:flex;gap:18px;flex-wrap:wrap;'
                   f'padding:2px 4px 12px">{cells}</div>')
    return (
        f'<div style="background:{SURFACE};border:1px solid {BORDER};'
        f'border-radius:14px;padding:16px 18px 14px;max-width:{max_width}px;'
        f'margin:14px auto;font-family:{FONT};'
        f'box-shadow:0 1px 2px rgba(11,11,11,.04),0 4px 16px rgba(11,11,11,.05)">'
        f'<div style="display:flex;align-items:flex-start;justify-content:space-between;'
        f'gap:12px;margin-bottom:2px">'
        f'<div><div style="font-size:16px;font-weight:680;color:{INK};'
        f'letter-spacing:-.01em;line-height:1.25">{title}</div>'
        f'<div style="font-size:12.5px;color:{SECONDARY};margin-top:2px">{subtitle}</div></div>'
        f'{badge}</div>'
        f'<div style="margin-top:10px">{kpi_row}{body_html}</div>'
        f'<div style="color:{MUTED};font-size:12px;line-height:1.55;'
        f'padding:10px 4px 2px;border-top:1px solid {GRID};margin-top:12px">{caption}</div>'
        f'</div>'
    )


def _hex_rgba(hex_color, alpha):
    h = hex_color.lstrip("#")
    r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
    return f"rgba({r},{g},{b},{alpha})"
