"""The endosymbiosis phase diagram — enclosure (G3) plane.

The campaign's capstone deliverable is a phase diagram separating the outcome
regimes of endosymbiosis over metabolic, physical, ecological, and evolutionary
conditions. This driver computes the PHYSICAL slice that the S3 enclosure model
can resolve directly: over the two mechanical control axes — host-symbiont
adhesion ``j_hs`` and symbiont-medium surface tension ``j_sM`` — it runs the real
CPM wrapping model as a small per-cell ensemble and classifies each parameter set
into an enclosure regime:

    no enclosure  ·  partial wrapping  ·  transient internalization  ·  stable enclosure

The differential-adhesion prediction is the line ``j_hs = j_sM - j_hM``: below it
the host engulfs the symbiont. The full campaign phase diagram overlays the
metabolic (G1), ecological (G2), energetic (G4) and evolutionary (G5) gates on
top of this physical plane; here we map the enclosure axis on its own.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from .enclosure import EnclosureParams, classify_regime, run_trajectory

REGIMES = ["no enclosure", "partial wrapping", "transient internalization", "stable enclosure"]
REGIME_IDX = {r: i for i, r in enumerate(REGIMES)}


def compute(j_hs_axis=(2, 6, 10, 12, 14, 18), j_sM_axis=(10, 14, 18, 22),
            j_hM=4.0, seeds=(7, 23), nx=80, total_sweeps=1800):
    """Return the regime grid + mean wrapping grid over the (j_hs, j_sM) plane."""
    regime = np.empty((len(j_sM_axis), len(j_hs_axis)), dtype=int)
    wrap = np.empty_like(regime, dtype=float)
    for iy, jsm in enumerate(j_sM_axis):
        for ix, jhs in enumerate(j_hs_axis):
            plateaus, regimes = [], []
            for sd in seeds:
                p = EnclosureParams(nx=nx, j_hs=float(jhs), j_sM=float(jsm),
                                    j_hM=float(j_hM), seed=int(sd))
                traj = run_trajectory(p, total_sweeps=total_sweeps, block=60)
                plateaus.append(float(np.mean(traj.wrapping[-max(1, len(traj.wrapping) // 5):])))
                regimes.append(classify_regime(traj))
            # ensemble: mean plateau, modal regime (ties -> the more-enclosed one)
            wrap[iy, ix] = float(np.mean(plateaus))
            ranked = sorted(regimes, key=lambda r: REGIME_IDX[r])
            regime[iy, ix] = REGIME_IDX[max(set(ranked), key=lambda r: (regimes.count(r), REGIME_IDX[r]))]
    return {
        "j_hs_axis": list(map(float, j_hs_axis)),
        "j_sM_axis": list(map(float, j_sM_axis)),
        "j_hM": float(j_hM), "nx": nx, "total_sweeps": total_sweeps, "seeds": list(seeds),
        "regimes": REGIMES,
        "regime_grid": regime.tolist(),
        "wrapping_grid": np.round(wrap, 3).tolist(),
    }


def render_png(data, out_path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import BoundaryNorm, ListedColormap

    jhs, jsm = data["j_hs_axis"], data["j_sM_axis"]
    grid = np.array(data["regime_grid"])
    colors = ["#f4f1ea", "#e4c27a", "#c6792f", "#2f6f8f"]  # no / partial / transient / stable
    cmap = ListedColormap(colors)
    norm = BoundaryNorm(range(len(REGIMES) + 1), cmap.N)

    fig, ax = plt.subplots(figsize=(7.4, 5.4))
    ax.imshow(grid, origin="lower", aspect="auto", cmap=cmap, norm=norm,
              extent=[min(jhs), max(jhs), min(jsm), max(jsm)], interpolation="nearest")
    # the differential-adhesion prediction line  j_hs = j_sM - j_hM
    yy = np.array(jsm, dtype=float)
    ax.plot(yy - data["j_hM"], yy, color="#1b1b1b", lw=2, ls="--",
            label=r"predicted threshold  $J_{hs}=J_{sM}-J_{hM}$")
    # annotate measured mean wrapping
    for iy, y in enumerate(jsm):
        for ix, x in enumerate(jhs):
            ax.text(x, y, f"{data['wrapping_grid'][iy][ix]:.2f}", ha="center", va="center",
                    fontsize=8, color="#111")
    ax.set_xlabel(r"host–symbiont adhesion  $J_{hs}$   (lower = stronger adhesion)")
    ax.set_ylabel(r"symbiont–medium surface tension  $J_{sM}$")
    ax.set_title("Endosymbiosis phase diagram — enclosure (G3) plane\n"
                 "CPM differential-adhesion wrapping (real cpm engine)")
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in colors]
    ax.legend(handles + [ax.lines[0]], REGIMES + ["predicted threshold"],
              loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize=9, frameon=False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    return out_path


def render_html(data, out_path):
    """A self-contained interactive Plotly heatmap for the workbench (viz/*.html)."""
    jhs, jsm = data["j_hs_axis"], data["j_sM_axis"]
    colors = ["#f4f1ea", "#e4c27a", "#c6792f", "#2f6f8f"]
    # discrete colorscale over the 4 regime indices
    n = len(REGIMES)
    colorscale = []
    for i, c in enumerate(colors):
        colorscale.append([i / n, c]); colorscale.append([(i + 1) / n, c])
    text = [[f"{REGIMES[data['regime_grid'][iy][ix]]}<br>wrap {data['wrapping_grid'][iy][ix]:.2f}"
             for ix in range(len(jhs))] for iy in range(len(jsm))]
    thr_x = [y - data["j_hM"] for y in jsm]
    fig = {
        "data": [
            {"type": "heatmap", "x": jhs, "y": jsm, "z": data["regime_grid"],
             "text": text, "hoverinfo": "text", "colorscale": colorscale,
             "zmin": 0, "zmax": n, "showscale": False, "xgap": 2, "ygap": 2},
            {"type": "scatter", "x": thr_x, "y": jsm, "mode": "lines",
             "line": {"color": "#111", "width": 2.5, "dash": "dash"},
             "name": "predicted threshold  J_hs = J_sM − J_hM", "hoverinfo": "name"},
        ],
        "layout": {
            "title": {"text": "<b>Endosymbiosis phase diagram — enclosure (G3) plane</b><br>"
                      "<span style='font-size:12px;color:#64748b'>CPM differential-adhesion wrapping; "
                      "colour = regime, number = measured host coverage of the symbiont boundary</span>"},
            "xaxis": {"title": "host–symbiont adhesion  J_hs  (lower = stronger adhesion)"},
            "yaxis": {"title": "symbiont–medium surface tension  J_sM"},
            "annotations": [
                {"x": jhs[ix], "y": jsm[iy], "text": f"{data['wrapping_grid'][iy][ix]:.2f}",
                 "showarrow": False, "font": {"size": 10, "color": "#111"}}
                for iy in range(len(jsm)) for ix in range(len(jhs))],
            "margin": {"l": 70, "r": 30, "t": 70, "b": 60},
            "paper_bgcolor": "white", "plot_bgcolor": "white",
        },
    }
    legend = "".join(
        f'<span style="display:inline-block;width:12px;height:12px;background:{c};'
        f'margin:0 4px 0 12px;vertical-align:middle"></span>{r}'
        for c, r in zip(colors, REGIMES))
    html = (
        "<!doctype html><html><head><meta charset='utf-8'>"
        "<script src='https://cdn.jsdelivr.net/npm/plotly.js-dist-min@2.35.2/plotly.min.js'></script>"
        "</head><body style='font-family:-apple-system,Segoe UI,Roboto,sans-serif;margin:0'>"
        "<div id='pd' style='width:100%;max-width:820px;margin:0 auto'></div>"
        f"<div style='text-align:center;font-size:12px;color:#334155;margin:6px'>{legend}</div>"
        f"<script>Plotly.newPlot('pd',{json.dumps(fig['data'])},{json.dumps(fig['layout'])},"
        "{responsive:true,displayModeBar:false});</script></body></html>")
    Path(out_path).write_text(html)
    return out_path


def main(out_dir=None):
    data = compute()
    out_dir = Path(out_dir or ".")
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "phase_diagram.json").write_text(json.dumps(data, indent=2))
    render_png(data, out_dir / "phase_diagram.png")
    (out_dir / "viz").mkdir(exist_ok=True)
    render_html(data, out_dir / "viz" / "endosymbiosis-phase-diagram.html")
    print("regimes (rows = j_sM low->high, cols = j_hs low->high):")
    for iy, row in enumerate(data["regime_grid"]):
        print(f"  j_sM={data['j_sM_axis'][iy]:>4}:  " +
              "  ".join(REGIMES[v][:4] for v in row))
    print("wrote", out_dir / "phase_diagram.json", "and phase_diagram.png")
    return data


if __name__ == "__main__":
    import sys
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
