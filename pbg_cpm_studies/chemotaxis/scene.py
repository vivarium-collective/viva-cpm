"""Generate the animated spatial-scene figure for chemotactic recruitment.

Renders the REAL Cellular Potts lattice (``world.snapshot()`` — the per-pixel
cell-label array), not just cell centroids: each cell is drawn as its actual
filled, irregular multi-pixel region, colored by type (source = gold, responder =
cyan), over the diffusing cue gradient which shows through the medium (blue ramp).
Two synchronized panels — baseline (cue+response) vs inhibited (cue only) — with a
play button and an MCS slider. In the baseline the responder cells crawl up the
gradient and pack against the source; in the inhibited panel they stay put.

Run (from the worktree root, viva-cpm venv)::

    PYTHONPATH=. python -m pbg_cpm_studies.chemotaxis.scene
"""
from __future__ import annotations

import json
import math
import os

from cpm import load_world

from ..composites.chemotaxis import build_spec, CUE_RATE, CHEMO_LAMBDA, NX, NY

WARMUP, SWEEPS, SAMPLE = 600, 500, 25

# pixel value encoding -> colorscale bands (medium shows the cue; cells are solid)
SRC_VAL, RESP_VAL, MED_MAX = 0.82, 1.0, 0.55
COLORSCALE = [
    [0.00, "#0a1020"], [0.20, "#0f2547"], [0.40, "#17539e"], [MED_MAX, "#3b82f6"],
    [MED_MAX + 0.01, "#0a1020"], [0.77, "#0a1020"],
    [0.78, "#f59e0b"], [0.90, "#fbbf24"],
    [0.905, "#22d3ee"], [1.0, "#67e8f9"],
]


def _type_by_id(world):
    types = list(world.cell_types())
    return types


def _field_scaled(world):
    raw = list(world.field_conc(0))
    logf = [math.log1p(max(0.0, v)) for v in raw]
    m = max(logf) or 1.0
    return [MED_MAX * (logf[i] / m) for i in range(len(logf))]  # flat, len nx*ny


def _z_frame(world, field_flat, types):
    """Per-pixel scalar: cue value in medium, SRC_VAL/RESP_VAL on cells. [ny][nx]."""
    lab = list(world.snapshot())
    z = []
    for y in range(NY):
        row = []
        for x in range(NX):
            i = y * NX + x
            cid = lab[i]
            t = types[cid] if cid < len(types) else 0
            if t == 1:
                row.append(SRC_VAL)
            elif t == 2:
                row.append(RESP_VAL)
            else:
                row.append(round(field_flat[i], 3))
        z.append(row)
    return z


def _capture(cue_rate, chemo_lambda):
    w = load_world(build_spec(cue_rate=cue_rate, chemo_lambda=chemo_lambda))
    for _ in range(WARMUP):
        w.advance_fields(1)
    field = _field_scaled(w)
    types = _type_by_id(w)
    zs = [{"mcs": 0, "z": _z_frame(w, field, _type_by_id(w))}]
    for s in range(SAMPLE, SWEEPS + 1, SAMPLE):
        w.step(SAMPLE)
        # re-read the cue in the medium too (it shifts slightly as cells move)
        zs.append({"mcs": s, "z": _z_frame(w, _field_scaled(w), _type_by_id(w))})
    return zs


def _heat(z, ax, ay):
    return {"z": z, "type": "heatmap", "xaxis": ax, "yaxis": ay,
            "colorscale": COLORSCALE, "zmin": 0.0, "zmax": 1.0,
            "showscale": False, "hoverinfo": "skip", "zsmooth": False}


def build_html(base, inh, right_label="inhibited · cue ✓ response ✗"):
    data = [_heat(base[0]["z"], "x", "y"), _heat(inh[0]["z"], "x2", "y2")]
    frames = [{"name": str(b["mcs"]),
               "data": [{"z": b["z"]}, {"z": i["z"]}], "traces": [0, 1]}
              for b, i in zip(base, inh)]
    axis = lambda dom, title: dict(domain=dom, range=[0, NX], showgrid=False, zeroline=False,
                                   showticklabels=False, constrain="domain",
                                   title={"text": title, "font": {"size": 12.5, "color": "#94a3b8"}})
    yax = dict(range=[0, NY], showgrid=False, zeroline=False, showticklabels=False)
    layout = {
        "title": {"text": "<b>Chemotactic recruitment</b> &nbsp;&mdash;&nbsp; a real Cellular Potts lattice; responder cells crawl up the secreted cue",
                  "x": 0.5, "xanchor": "center", "font": {"size": 16, "color": "#e5e7eb"}},
        "paper_bgcolor": "#0a1020", "plot_bgcolor": "#0a1020",
        "font": {"family": "-apple-system, Segoe UI, Roboto, sans-serif", "color": "#e5e7eb"},
        "margin": {"l": 8, "r": 8, "t": 92, "b": 74},
        "xaxis": axis([0.0, 0.485], "baseline · cue ✓ response ✓"),
        "yaxis": {**yax, "anchor": "x"},
        "xaxis2": axis([0.515, 1.0], right_label),
        "yaxis2": {**yax, "anchor": "x2"},
        "height": 430,
        "annotations": [
            {"xref": "paper", "yref": "paper", "x": 0.0, "y": 1.03, "showarrow": False, "xanchor": "left",
             "text": "■ source cells (secrete cue)", "font": {"color": "#fbbf24", "size": 11.5}},
            {"xref": "paper", "yref": "paper", "x": 0.30, "y": 1.03, "showarrow": False, "xanchor": "left",
             "text": "■ responder cells", "font": {"color": "#22d3ee", "size": 11.5}},
            {"xref": "paper", "yref": "paper", "x": 0.62, "y": 1.03, "showarrow": False, "xanchor": "left",
             "text": "cue gradient (blue → bright near source)", "font": {"color": "#3b82f6", "size": 11.5}},
        ],
        "updatemenus": [{
            "type": "buttons", "showactive": False, "y": -0.08, "x": 0.0, "xanchor": "left",
            "bgcolor": "#1e293b", "font": {"color": "#e5e7eb"},
            "buttons": [
                {"label": "▶ Play", "method": "animate",
                 "args": [None, {"frame": {"duration": 240, "redraw": True}, "fromcurrent": True,
                                 "transition": {"duration": 0}}]},
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
        '<style>body{margin:0;background:#0a1020}</style></head><body>'
        '<div id="scene" style="max-width:980px;margin:0 auto"></div><script>'
        f'Plotly.newPlot("scene",{json.dumps(data)},{json.dumps(layout)},'
        '{responsive:true,displayModeBar:false})'
        f'.then(function(gd){{Plotly.addFrames("scene",{json.dumps(frames)});}});'
        '</script></body></html>'
    )


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    ws = os.path.abspath(os.path.join(here, "..", "..", "workspace"))
    base = _capture(CUE_RATE, CHEMO_LAMBDA)
    inh = _capture(CUE_RATE, 0.0)
    adv = _capture(0.0, CHEMO_LAMBDA)
    # one animated scene per study, each contrasting the baseline with the study's
    # own condition (so inhibited/adversarial are not just the shared chart).
    scenes = {
        "recruitment-baseline":    (inh, "inhibited · cue ✓ response ✗"),
        "recruitment-inhibited":   (inh, "inhibited · cue ✓ response ✗ (this study)"),
        "recruitment-adversarial": (adv, "adversarial · cue ✗ response ✓ (this study)"),
    }
    for study, (right, label) in scenes.items():
        html = build_html(base, right, right_label=label)
        out = os.path.join(ws, "studies", study, "viz")
        os.makedirs(out, exist_ok=True)
        with open(os.path.join(out, "scene.html"), "w") as fh:
            fh.write(html)
        print(f"wrote {study} scene ({len(html)} bytes)")


if __name__ == "__main__":
    main()
