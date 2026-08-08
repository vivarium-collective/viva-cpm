"""Generate the animated spatial-scene figure for chemotactic recruitment.

Runs the baseline and inhibited conditions, capturing the equilibrated cue field
(a heatmap background) and the per-sweep cell centres of mass, then emits a
self-contained interactive HTML: two synchronized panels (cue+response vs cue-only)
with a play button and an MCS slider. In the baseline the responders (cyan) stream
up the gradient toward the source (gold); in the inhibited panel they stay put.

Run (from the worktree root, viva-cpm venv)::

    PYTHONPATH=. python -m pbg_cpm_studies.chemotaxis.scene
"""
from __future__ import annotations

import json
import os

from cpm import load_world

from ..composites.chemotaxis import build_spec, CUE_RATE, CHEMO_LAMBDA, NX, NY, SOURCE_TYPE

WARMUP, SWEEPS, SAMPLE = 600, 500, 25
SRC_COLOR, RESP_COLOR = "#fbbf24", "#22d3ee"


def _capture(cue_rate, chemo_lambda, want_field=False):
    w = load_world(build_spec(cue_rate=cue_rate, chemo_lambda=chemo_lambda))
    for _ in range(WARMUP):
        w.advance_fields(1)
    field = None
    if want_field:
        raw = list(w.field_conc(0))
        field = [raw[y * NX:(y + 1) * NX] for y in range(NY)]  # [ny][nx]
    frames = []

    def snap(mcs):
        coms = list(w.cell_coms()); types = list(w.cell_types())
        pts = [(round(coms[c][0], 1), round(coms[c][1], 1), int(types[c]))
               for c in range(1, len(types)) if int(types[c]) != 0]
        frames.append({"mcs": mcs, "pts": pts})

    snap(0)
    for s in range(SAMPLE, SWEEPS + 1, SAMPLE):
        w.step(SAMPLE); snap(s)
    return field, frames


def _cells_trace(pts, xaxis, yaxis, showlegend):
    xs = [p[0] for p in pts]; ys = [p[1] for p in pts]
    colors = [SRC_COLOR if p[2] == SOURCE_TYPE else RESP_COLOR for p in pts]
    sizes = [22 if p[2] == SOURCE_TYPE else 15 for p in pts]
    return {
        "x": xs, "y": ys, "mode": "markers", "type": "scatter",
        "xaxis": xaxis, "yaxis": yaxis, "showlegend": False, "hoverinfo": "skip",
        "marker": {"color": colors, "size": sizes, "line": {"color": "#0b1220", "width": 1.5},
                   "opacity": 0.96},
    }


def build_html(field, base_frames, inh_frames):
    import math
    logf = [[math.log1p(v) for v in row] for row in field]
    zmax = max(max(r) for r in logf)
    heat = lambda ax, ay: {
        "z": logf, "type": "heatmap", "xaxis": ax, "yaxis": ay,
        "colorscale": [[0, "#0b1220"], [0.15, "#0b2545"], [0.5, "#13315c"],
                       [0.8, "#3b82f6"], [1, "#93c5fd"]],
        "zmin": 0, "zmax": zmax, "showscale": False, "hoverinfo": "skip", "zsmooth": "best",
    }
    data = [heat("x", "y"), _cells_trace(base_frames[0]["pts"], "x", "y", True),
            heat("x2", "y2"), _cells_trace(inh_frames[0]["pts"], "x2", "y2", False)]
    frames = []
    for bf, inf in zip(base_frames, inh_frames):
        frames.append({
            "name": str(bf["mcs"]),
            "data": [_cells_trace(bf["pts"], "x", "y", False),
                     _cells_trace(inf["pts"], "x2", "y2", False)],
            "traces": [1, 3],
        })
    axis = lambda dom, title: dict(domain=dom, range=[0, NX], showgrid=False, zeroline=False,
                                   showticklabels=False, title={"text": title, "font": {"size": 12.5, "color": "#94a3b8"}})
    yax = dict(range=[0, NY], showgrid=False, zeroline=False, showticklabels=False, scaleanchor=None)
    layout = {
        "title": {"text": "<b>Chemotactic recruitment</b> &nbsp;&mdash;&nbsp; responders stream up the secreted gradient",
                  "x": 0.5, "xanchor": "center", "font": {"size": 17, "color": "#e5e7eb"}},
        "paper_bgcolor": "#0b1220", "plot_bgcolor": "#0b1220",
        "font": {"family": "-apple-system, Segoe UI, Roboto, sans-serif", "color": "#e5e7eb"},
        "margin": {"l": 10, "r": 10, "t": 92, "b": 70},
        "xaxis": axis([0.0, 0.485], "baseline · cue ✓ response ✓"),
        "yaxis": {**yax, "anchor": "x"},
        "xaxis2": axis([0.515, 1.0], "inhibited · cue ✓ response ✗"),
        "yaxis2": {**yax, "anchor": "x2"},
        "height": 440,
        "annotations": [
            {"xref": "paper", "yref": "paper", "x": 0.02, "y": 1.02, "showarrow": False,
             "text": "◆ source (secretes cue)", "font": {"color": SRC_COLOR, "size": 11.5}, "xanchor": "left"},
            {"xref": "paper", "yref": "paper", "x": 0.30, "y": 1.02, "showarrow": False,
             "text": "● responder", "font": {"color": RESP_COLOR, "size": 11.5}, "xanchor": "left"},
        ],
        "updatemenus": [{
            "type": "buttons", "showactive": False, "y": -0.08, "x": 0.0, "xanchor": "left",
            "bgcolor": "#1e293b", "font": {"color": "#e5e7eb"},
            "buttons": [
                {"label": "▶ Play", "method": "animate",
                 "args": [None, {"frame": {"duration": 220, "redraw": True}, "fromcurrent": True,
                                 "transition": {"duration": 140, "easing": "cubic-in-out"}}]},
                {"label": "❚❚ Pause", "method": "animate",
                 "args": [[None], {"frame": {"duration": 0, "redraw": False}, "mode": "immediate"}]},
            ],
        }],
        "sliders": [{
            "active": 0, "y": -0.06, "x": 0.12, "len": 0.86,
            "currentvalue": {"prefix": "MCS ", "font": {"color": "#94a3b8", "size": 12}},
            "font": {"color": "#94a3b8", "size": 10},
            "steps": [{"label": f["name"], "method": "animate",
                       "args": [[f["name"]], {"frame": {"duration": 0, "redraw": True},
                                              "mode": "immediate"}]} for f in frames],
        }],
    }
    return (
        '<!doctype html><html><head><meta charset="utf-8">'
        '<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>'
        '<style>body{margin:0;background:#0b1220}</style></head><body>'
        '<div id="scene" style="max-width:960px;margin:0 auto"></div><script>'
        f'Plotly.newPlot("scene",{json.dumps(data)},{json.dumps(layout)},'
        '{responsive:true,displayModeBar:false})'
        f'.then(function(gd){{Plotly.addFrames("scene",{json.dumps(frames)});}});'
        '</script></body></html>'
    )


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    ws = os.path.abspath(os.path.join(here, "..", "..", "workspace"))
    field, base_frames = _capture(CUE_RATE, CHEMO_LAMBDA, want_field=True)
    _, inh_frames = _capture(CUE_RATE, 0.0)
    html = build_html(field, base_frames, inh_frames)
    out = os.path.join(ws, "studies", "recruitment-baseline", "viz")
    os.makedirs(out, exist_ok=True)
    path = os.path.join(out, "scene.html")
    with open(path, "w") as fh:
        fh.write(html)
    print(f"wrote animated scene ({len(html)} bytes, {len(base_frames)} frames) -> {path}")
    return path


if __name__ == "__main__":
    main()
