"""Visualizations for the merks-ecm-reciprocity-2d investigation.

Every figure is built from REAL output of the bridged TST-MD engine (the
published 2D static-adhesion hybrid CPM + bead-spring ECM of Tsingos et al.
2023, driven CPU-only via viva-tstmd / Docker+HOOMD), baked into
``_tstmd_data.py``. These reproduce, at the mechanism level and in 2D, the
phenomena of Keijzer & Merks 2026 (arXiv:2609.02375): cell contraction
deforming a crosslinked ECM, fibers reorienting toward the cell (q/q0 > 1),
and crosslinking building network connectivity.
"""
from __future__ import annotations

import plotly.graph_objects as go
from plotly.subplots import make_subplots

from viva_superpowers.visualization import as_visualization

from ._tstmd_data import TSTMD_DATA

# palette (readable in light + dark)
C_AREA = "#e4572e"   # contraction
C_Q = "#1b9e77"      # reorientation q/q0
C_GC = "#3773b8"     # connectivity
C_FIBER = "rgba(130,130,130,0.45)"
C_CELL = "rgba(228,87,46,0.85)"


def _snapshot_trace(snap, name):
    fib = snap["fibers"]
    cell = snap["cell"]
    fx = [p[0] for p in fib]; fy = [p[1] for p in fib]
    cx = [p[0] for p in cell]; cy = [p[1] for p in cell]
    return (
        go.Scatter(x=fx, y=fy, mode="markers", name="ECM fibers",
                   marker=dict(size=2.5, color=C_FIBER), showlegend=False,
                   hoverinfo="skip"),
        go.Scatter(x=cx, y=cy, mode="markers", name="cell",
                   marker=dict(size=3, color=C_CELL), showlegend=False,
                   hoverinfo="skip"),
    )


def _build_contraction_reorientation() -> str:
    dyn = TSTMD_DATA["dynamics"]
    traj = dyn["trajectory"]
    mcs = [m["mcs"] for m in traj]
    area = [m["cell_area"] for m in traj]
    qratio = [m["q_ratio"] for m in traj]
    start, end = dyn["start_snapshot"], dyn["end_snapshot"]
    L = start["L"]

    fig = make_subplots(
        rows=2, cols=2,
        specs=[[{"colspan": 2, "secondary_y": True}, None], [{}, {}]],
        row_heights=[0.52, 0.48], vertical_spacing=0.14, horizontal_spacing=0.08,
        subplot_titles=(
            "Cell contracts while fibers reorient toward it",
            f"ECM + cell at start (MCS {start['mcs']})",
            f"ECM + cell at end (MCS {end['mcs']})",
        ),
    )

    fig.add_trace(go.Scatter(x=mcs, y=area, mode="lines+markers", name="cell area",
                             line=dict(color=C_AREA, width=3),
                             marker=dict(size=6)), row=1, col=1, secondary_y=False)
    fig.add_trace(go.Scatter(x=mcs, y=qratio, mode="lines+markers", name="reorientation q/q₀",
                             line=dict(color=C_Q, width=3, dash="dot"),
                             marker=dict(size=6)), row=1, col=1, secondary_y=True)
    fig.add_hline(y=1.0, line=dict(color=C_Q, width=1, dash="dash"),
                  opacity=0.4, row=1, col=1, secondary_y=True)

    for tr in _snapshot_trace(start, "start"):
        fig.add_trace(tr, row=2, col=1)
    for tr in _snapshot_trace(end, "end"):
        fig.add_trace(tr, row=2, col=2)

    fig.update_xaxes(title_text="Monte Carlo step (MCS)", row=1, col=1)
    fig.update_yaxes(title_text="cell area (lattice sites)", color=C_AREA,
                     row=1, col=1, secondary_y=False)
    fig.update_yaxes(title_text="q/q₀  (fibers toward cell)", color=C_Q,
                     row=1, col=1, secondary_y=True)
    for c in (1, 2):
        fig.update_xaxes(range=[-10, L + 10], row=2, col=c, scaleanchor=f"y{3 if c==1 else 4}",
                         showticklabels=False)
        fig.update_yaxes(range=[-10, L + 10], row=2, col=c, showticklabels=False)

    fig.update_layout(
        template="plotly_white", height=760,
        title=dict(text="<b>Mechanical reciprocity (2D):</b> a contractile cell "
                        "remodels a real crosslinked ECM", x=0.02, font=dict(size=17)),
        legend=dict(orientation="h", y=1.07, x=0.0),
        margin=dict(l=70, r=70, t=90, b=50),
    )
    return fig.to_html(full_html=True, include_plotlyjs="cdn",
                       config={"displayModeBar": False})


def _build_network_connectivity() -> str:
    conn = TSTMD_DATA["connectivity"]
    created = [c["created"] for c in conn]
    gc = [c["giant_component"] for c in conn]

    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=created, y=gc, mode="lines+markers",
        line=dict(color=C_GC, width=3), marker=dict(size=10, color=C_GC),
        hovertemplate="crosslinks: %{x}<br>giant component: %{y:.3f}<extra></extra>",
        name="giant component"))
    fig.add_annotation(x=created[0], y=gc[0], text="no crosslinks →<br>disconnected",
                       showarrow=True, arrowhead=2, ax=50, ay=-30, font=dict(size=11))
    fig.update_layout(
        template="plotly_white", height=460,
        title=dict(text="<b>Crosslinking builds ECM network connectivity</b><br>"
                        "<span style='font-size:12px;color:#666'>Giant component rises sharply "
                        "once crosslinks appear, then saturates (~0.3) at these reduced densities "
                        "— a percolation-onset signature, not full percolation.</span>",
                   x=0.02, font=dict(size=16)),
        xaxis=dict(title="crosslinks created", type="log"),
        yaxis=dict(title="giant-component fraction", range=[0, 1]),
        margin=dict(l=70, r=40, t=90, b=55),
    )
    return fig.to_html(full_html=True, include_plotlyjs="cdn",
                       config={"displayModeBar": False})


@as_visualization(inputs={"mcs": "list[float]"},
                  name="TstmdContractionReorientation", demo={"mcs": [0.0]})
def update_tstmd_contraction_reorientation(state):
    """Real TST-MD: cell contracts while fibers reorient toward it (q/q0 > 1)."""
    return {"html": _build_contraction_reorientation()}


@as_visualization(inputs={"mcs": "list[float]"},
                  name="TstmdNetworkConnectivity", demo={"mcs": [0.0]})
def update_tstmd_network_connectivity(state):
    """Real TST-MD: crosslinking builds network connectivity (percolation onset)."""
    return {"html": _build_network_connectivity()}
