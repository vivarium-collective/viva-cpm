"""Influenza (Sego 2022) study visualizations — one interactive system.

A consistent, high-polish Plotly visualization per influenza-sego2022 study,
built on the shared shell in ``_influenza_style``. Data are REAL engine output
baked in ``_influenza_data`` (captured from pbg_cpm_studies.influenza at
patch_mm=0.3, seed=17 — the values the studies report), so figures render
identically in the live and the published read-only dashboard.

Studies covered: parameter-provenance (Increment 0), epithelial-sheet-baseline
(Increment 1), virus-field-infection (Increment 2). The system is N-study: the
immune-arm studies (IFN/resistance, macrophage/NK/CD8) drop in with one more
subclass each, reusing the same shell + palette.
"""
from __future__ import annotations
import json

from viva_superpowers.visualization import as_visualization

from . import _influenza_style as S
from ._influenza_data import INFLUENZA_DATA

_D = json.loads(INFLUENZA_DATA)


# ── Study A · parameter-provenance (Increment 0): Fig-3B acceptance targets ──
def _build_fig3b():
    """Digitized Fig-3B acceptance bands the Increment-9 capstone must hit"""
    fig = _D["targets"]["fig3b"]
    obs = fig["observables"]
    series = [("uninfected_cells", "Healthy epithelial", S.HEALTHY),
              ("infected_cells", "Infected epithelial", S.INFECTED),
              ("dead_cells", "Dead epithelial", S.DEAD)]
    traces = []
    for key, label, color in series:
        pts = obs.get(key) or []
        if not pts:
            continue
        t = [p["t_days"] for p in pts]
        val = [p["value"] for p in pts]
        lo = [p.get("lo", p["value"]) for p in pts]
        hi = [p.get("hi", p["value"]) for p in pts]
        traces += S.band_trace(t, lo, hi, color, label)
        traces.append(S.line_trace(t, val, color, label, markers=False))
    layout = S.base_layout("Time (days)", "Epithelial cells", height=400)
    body = S.plot("influenza-fig3b", traces, layout)
    scen = fig.get("scenario", {})
    n_disc = len(_D["params"].get("scenarios", {}).get("initial_viral_load", [])) or 0
    kpis = [("3 figure targets", "Figs 3B · 5 · 7"),
            (f'{len(obs)} observables', "digitized w/ bands"),
            ("8 discrepancies", "source vs paper, open")]
    cap = (f'<b>Targets, not results.</b> Digitized from <b>{fig.get("source","Sego 2022 Fig 3B")}</b> '
           f'(ODE reference line + shaded 50-replica band), scenario: init-infection '
           f'{scen.get("init_infection_fraction","0.05")}, {scen.get("cells","1225")} cells. '
           f'These are the acceptance bands the native reproduction is graded against at the '
           f'Increment-9 capstone — no reproduction is claimed yet.')
    return S.card("Acceptance-band targets (Sego 2022, Fig 3B)",
                           "What the capstone must reproduce — the parameter/target authority",
                           body, cap, increment=0, kpis=kpis)


# ── Study B · epithelial-sheet-baseline (Increment 1): confluent sheet ──────
def _build_sheet():
    """Confluent epithelial-sheet substrate — the cell mosaic at scale"""
    sheet = _D["sheet"]
    nx, ny = int(sheet["dims"][0]), int(sheet["dims"][1])
    labels = sheet["snapshot"]
    # Mosaic: hash each cell label into a repeating soft-hue band so adjacent
    # cells read as distinct tiles; medium (0) recedes to the surface.
    z = [[0 if labels[y * nx + x] == 0 else 1 + (labels[y * nx + x] * 2654435761 % 11)
          for x in range(nx)] for y in range(ny)]
    mosaic = [
        [0.0, S.SURFACE], [0.0001, "#d7f0e6"], [0.1, "#cfe6fb"], [0.2, "#e8dcfb"],
        [0.3, "#fde3cf"], [0.4, "#d7f0e6"], [0.5, "#fce0ec"], [0.6, "#e3eccf"],
        [0.7, "#cfe6fb"], [0.8, "#fdeecf"], [0.9, "#e8dcfb"], [1.0, "#d7f0e6"],
    ]
    heat = {"type": "heatmap", "z": z, "colorscale": mosaic, "showscale": False,
            "hoverinfo": "skip", "xgap": 0, "ygap": 0}
    layout = S.base_layout(height=460, showlegend=False, extra={
        "margin": {"l": 20, "r": 20, "t": 10, "b": 20},
        "xaxis": {"visible": False, "scaleanchor": "y", "constrain": "domain"},
        "yaxis": {"visible": False, "autorange": "reversed"},
    })
    body = S.plot("influenza-sheet", [heat], layout)
    vols = [v for v in sheet["cell_volumes"][1:]]
    mean_v = round(sum(vols) / len(vols), 1) if vols else 0
    runs = sheet["throughput_runs_mcs_per_s"]
    kpis = [(f'{len(vols):,} cells', f'{sheet["patch_mm"]} mm² confluent sheet'),
            (f'{mean_v} sites/cell', 'mean volume (target 25)'),
            (f'{min(runs):.0f}–{max(runs):.0f} MCS/s',
             f'throughput — clears {sheet["throughput_floor"]:.0f} floor')]
    cap = (f'The Increment-1 substrate: a confluent epithelial tiling on a '
           f'{nx}×{ny} lattice ({len(vols):,} five-by-five cells), each patch a distinct '
           f'tile above. Geometry holds (mean volume {mean_v} vs the target 25) and 1 mm² '
           f'throughput ({min(runs):.0f}–{max(runs):.0f} MCS/s across {len(runs)} runs) clears the '
           f'{sheet["throughput_floor"]:.0f} MCS/s feasibility floor. Substrate only — no biology yet.')
    return S.card("Confluent epithelial sheet at scale",
                           "The substrate every mechanism runs on — geometry + throughput",
                           body, cap, increment=1, kpis=kpis)


# ── Study C · virus-field-infection (Increment 2): animated spread scene ────
def _build_scene():
    """Virus field diffusing + the lesion spreading locally, over time"""
    inf = _D["infection"]
    frames_data = inf["frames"]
    fg = inf["field_grid"]
    gnx, gny = fg["nx"], fg["ny"]
    fac = fg["coarsen"]
    xs = [fac * j + fac / 2 for j in range(gnx)]
    ys = [fac * j + fac / 2 for j in range(gny)]
    zmax = max((max(max(row) for row in fr["field"]) for fr in frames_data), default=1.0) or 1.0

    def _state_color(t):
        return S.CELL_STATES.get(t, S.CELL_STATES[1])[1]

    def _frame_traces(fr):
        heat = {"type": "heatmap", "z": fr["field"], "x": xs, "y": ys,
                "colorscale": S.VIRUS_SCALE, "zmin": 0, "zmax": zmax,
                "colorbar": {"title": {"text": "virus", "side": "right",
                             "font": {"size": 11, "color": S.SECONDARY}},
                             "thickness": 12, "len": 0.85, "outlinewidth": 0,
                             "tickfont": {"size": 10, "color": S.MUTED}},
                "hovertemplate": "virus %{z:.3f}<extra></extra>"}
        cells = fr["cells"]
        cx = [c["x"] for c in cells]
        cy = [c["y"] for c in cells]
        colors = [_state_color(c["t"]) for c in cells]
        labels = [S.CELL_STATES.get(c["t"], S.CELL_STATES[1])[0] for c in cells]
        scat = {"type": "scatter", "x": cx, "y": cy, "mode": "markers",
                "marker": {"color": colors, "size": 5,
                           "line": {"color": "rgba(255,255,255,0.55)", "width": 0.5}},
                "text": labels, "hovertemplate": "%{text}<extra></extra>",
                "showlegend": False}
        return [heat, scat]

    data0 = _frame_traces(frames_data[0])
    frames = [{"name": str(fr["step"]), "data": _frame_traces(fr)} for fr in frames_data]
    steps = [{"label": str(fr["step"]), "method": "animate",
              "args": [[str(fr["step"])],
                       {"mode": "immediate", "frame": {"duration": 0, "redraw": True},
                        "transition": {"duration": 0}}]} for fr in frames_data]
    layout = S.base_layout(height=520, showlegend=False, extra={
        "margin": {"l": 20, "r": 20, "t": 10, "b": 60},
        "xaxis": {"visible": False, "scaleanchor": "y", "constrain": "domain",
                  "range": [0, fac * gnx]},
        "yaxis": {"visible": False, "autorange": "reversed", "range": [fac * gny, 0]},
        "updatemenus": [{"type": "buttons", "showactive": False, "x": 0.0, "y": -0.04,
                         "xanchor": "left", "yanchor": "top", "direction": "left",
                         "pad": {"t": 4}, "buttons": [
            {"label": "▶ Play", "method": "animate",
             "args": [None, {"mode": "immediate", "fromcurrent": True,
                             "frame": {"duration": 420, "redraw": True},
                             "transition": {"duration": 0}}]},
            {"label": "❙❙ Pause", "method": "animate",
             "args": [[None], {"mode": "immediate", "frame": {"duration": 0, "redraw": False},
                               "transition": {"duration": 0}}]}]}],
        "sliders": [{"active": 0, "x": 0.26, "len": 0.72, "y": -0.02,
                     "xanchor": "left", "yanchor": "top", "pad": {"t": 4},
                     "currentvalue": {"prefix": "update ", "xanchor": "left",
                                      "font": {"size": 12, "color": S.SECONDARY}},
                     "steps": steps}],
    })
    body = S.plot("influenza-scene", data0, layout, frames=frames)
    loc = inf["locality"]
    ser = inf["series"]
    kpis = [(f'{ser["n_I"][0]} → {ser["n_I"][-1]}', 'infected cells (60 updates)'),
            (f'{loc["null_ratio"]:.2f}', f'locality ratio (< {loc["pass_threshold"]} = local)'),
            (f'{loc["virus_diffusion_length_sites"]} sites', 'virus diffusion length')]
    cap = (f'<b>Play it.</b> The virus field (blue) diffuses from the seeded lesion; healthy cells '
           f'(<span style="color:{S.HEALTHY}">●</span>) turn infected '
           f'(<span style="color:{S.INFECTED}">●</span>) stochastically where the field is strong. '
           f'New infections land a mean {inf["actual_mean_dist"]} sites from a prior-infected cell — '
           f'{loc["null_ratio"]:.2f}× the {loc["null_mean"]}-site uniform-random null (a shuffled-field '
           f'control gives {loc["shuffled_field_control_ratio"]:.2f}), so spread is genuinely local. '
           f'Mechanism only — no resistance/IFN or death yet.')
    return S.card("Virus field + local infection spread",
                           "Increment-2 mechanism — animated over 60 updates (patch 0.3 mm, seed 17)",
                           body, cap, increment=2, kpis=kpis)


# ── Study C companion: population conservation + virus load + locality ──────
def _build_dynamics():
    """Population conservation, virus load, and the locality null-baseline"""
    inf = _D["infection"]
    ser = inf["series"]
    x = ser["steps"]
    pop = [S.line_trace(x, ser["n_H"], S.HEALTHY, "Healthy"),
           S.line_trace(x, ser["n_I"], S.INFECTED, "Infected")]
    pop_layout = S.base_layout("Update", "Cell count", height=300)
    virus = [S.line_trace(x, ser["total_virus"], "#256abf", "Total virus", markers=False)]
    virus_layout = S.base_layout("Update", "Virus field (Σ concentration)",
                                 height=260, showlegend=False)
    body = ('<div style="font-size:12px;color:%s;font-weight:600;margin:2px 4px 0">'
            'Population is conserved (n_H + n_I constant)</div>' % S.SECONDARY
            + S.plot("influenza-pop", pop, pop_layout)
            + '<div style="font-size:12px;color:%s;font-weight:600;margin:8px 4px 0">'
              'Extracellular virus accumulates</div>' % S.SECONDARY
            + S.plot("influenza-virus", virus, virus_layout))
    loc = inf["locality"]
    kpis = [(f'{ser["n_H"][0] + ser["n_I"][0]}', 'total cells (conserved)'),
            (f'{ser["total_virus"][-1]:,.0f}', 'final virus load'),
            (f'{loc["null_ratio"]:.2f} < {loc["pass_threshold"]}', 'locality: PASS')]
    cap = (f'Over 60 updates the infected count rises {ser["n_I"][0]} → {ser["n_I"][-1]} while the total '
           f'({ser["n_H"][0] + ser["n_I"][0]}) is conserved every step — no cells appear or vanish. Virus '
           f'accumulates to {ser["total_virus"][-1]:,.0f}. Locality null-baseline: observed mean spread '
           f'{inf["actual_mean_dist"]} sites vs {loc["null_mean"]} for uniform-random placement '
           f'(ratio {loc["null_ratio"]:.2f}, below the {loc["pass_threshold"]} threshold), '
           f'while a shuffled-field control fails at {loc["shuffled_field_control_ratio"]:.2f}.')
    return S.card("Infection dynamics + locality control",
                           "Increment-2 quantitative readout — conservation, load, spread locality",
                           body, cap, increment=2, kpis=kpis)


# ── registry wrappers (auto-discovered v2 Visualizations) ───────────────────
@as_visualization(inputs={"mcs": "list[float]"}, name="InfluenzaFig3BTargets", demo={"mcs": [0.0]})
def update_influenza_fig3b_targets(state):
    """Digitized Fig-3B acceptance bands the Increment-9 capstone must hit"""
    return {"html": _build_fig3b()}


@as_visualization(inputs={"mcs": "list[float]"}, name="InfluenzaEpithelialSheet", demo={"mcs": [0.0]})
def update_influenza_epithelial_sheet(state):
    """Confluent epithelial-sheet substrate — the cell mosaic at scale"""
    return {"html": _build_sheet()}


@as_visualization(inputs={"mcs": "list[float]"}, name="InfluenzaVirusFieldScene", demo={"mcs": [0.0]})
def update_influenza_virus_field_scene(state):
    """Virus field diffusing + the lesion spreading locally, over time"""
    return {"html": _build_scene()}


@as_visualization(inputs={"mcs": "list[float]"}, name="InfluenzaInfectionDynamics", demo={"mcs": [0.0]})
def update_influenza_infection_dynamics(state):
    """Population conservation, virus load, and the locality null-baseline"""
    return {"html": _build_dynamics()}


# ── zero-arg accessors (exercise outside the Step lifecycle: tests + file render) ──
def InfluenzaFig3BTargets():
    return _build_fig3b()


def InfluenzaEpithelialSheet():
    return _build_sheet()


def InfluenzaVirusFieldScene():
    return _build_scene()


def InfluenzaInfectionDynamics():
    return _build_dynamics()
