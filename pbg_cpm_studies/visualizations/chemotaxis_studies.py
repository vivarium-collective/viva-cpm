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


def _claim_layer_atlas():
    box = ("background:#0f172a;border:1px solid #1e293b;border-radius:12px;"
           "padding:12px 14px;margin:0 0 10px")
    lab = "font-size:11px;letter-spacing:.08em;text-transform:uppercase;color:#64748b"
    body = "color:#e5e7eb;font-size:13.5px;line-height:1.5;margin-top:3px"
    arrow = ('<div style="text-align:center;color:#334155;font-size:18px;'
             'margin:-4px 0 2px">&#8595;</div>')
    step = (lambda tag, color, head, txt:
            f'<div style="{box};border-left:3px solid {color}">'
            f'<div style="{lab}">{tag}</div>'
            f'<div style="{body}"><b style="color:{color}">{head}</b> &mdash; {txt}</div></div>')
    responds = "".join(
        f'<li style="margin:3px 0"><b style="color:#e5e7eb">{k}</b> '
        f'<span style="color:#94a3b8">{v}</span></li>' for k, v in [
            ("Claim is the source", "the biology, not the CPM model, is canonical"),
            ("Roles &amp; capabilities", "mechanism template names no lattice/pixels"),
            ("Semantic gap recorded", "&#955; labelled phenomenological, not hidden"),
            ("Backend-neutral contract", "recruitment_index via a per-backend extractor"),
            ("Alternatives kept", "PhysiCell / receptor-ODE realize the same claim"),
            ("Adversarial + falsifiers", "the claim can be, and is, stress-tested"),
        ])
    return (
        f'<div style="background:#0b1220;border:1px solid #1e293b;border-radius:14px;'
        f'padding:14px 16px;max-width:780px;margin:0 auto;'
        f'font-family:-apple-system,Segoe UI,Roboto,sans-serif">'
        f'<div style="color:#e5e7eb;font-size:16px;font-weight:700;margin-bottom:4px">'
        f'Claim layer &mdash; claim &rarr; mechanism &rarr; realization</div>'
        f'<div style="color:#94a3b8;font-size:12.5px;margin-bottom:12px">This investigation is '
        f'the reference vertical slice of the Biological Claim Layer: the claim is authored '
        f'backend-neutral; the CPM studies are one realization a compiler would generate.</div>'
        + step("Claim &middot; what is asserted", "#38bdf8", "A secreted cue recruits responder cells",
               "blocking the response, or removing the cue, abolishes recruitment. "
               "Context, intervention, observable, alternatives, applicability, falsifiers.") + arrow
        + step("Mechanism &middot; roles + capabilities", "#a78bfa", "Chemotactic response",
               "a motile responder, a diffusible cue, a gradient-query environment "
               "&mdash; no lattice, no pixels, no backend.") + arrow
        + step("Realization &middot; viva-cpm (fidelity + gap)", "#fbbf24", "CPM binding",
               "cue = secretion field, response = chemotaxis &#955; (a copy-attempt "
               "energy bias). <b>Semantic gap:</b> &#955; is <i>phenomenological</i>, "
               "not receptor-level &mdash; recorded, not hidden.") + arrow
        + step("Studies &middot; the realization, run", "#34d399", "baseline &middot; inhibited &middot; adversarial",
               "positive &middot; intervention on the response &middot; negative control "
               "(no cue). A future claim compiler would generate these three.")
        + f'<div style="{box};border-left:3px solid #64748b;margin-top:4px">'
        f'<div style="{lab}">Responds to the Biological-Modeling design doc</div>'
        f'<ul style="{body};padding-left:18px;margin:6px 0 2px">{responds}</ul></div>'
        f'</div>'
    )


@as_visualization(inputs={"mcs": "list[float]"}, name="ClaimLayerAtlas", demo={"mcs": []})
def update_claim_layer_atlas(state):
    """Claim layer — claim -> mechanism -> realization (with the recorded semantic gap)"""
    return {"html": _claim_layer_atlas()}
