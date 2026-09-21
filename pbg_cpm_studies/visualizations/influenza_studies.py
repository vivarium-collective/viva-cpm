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
import base64
import json
import zlib

import numpy as np

from viva_superpowers.visualization import as_visualization

from . import _influenza_style as S
from ._influenza_data import INFLUENZA_DATA, INFLUENZA_STUDY_DATA
from ._influenza_spatial import INFLUENZA_SPATIAL

_D = json.loads(INFLUENZA_DATA)
# Increments 3-8 study readouts (compact REAL engine series, no field frames).
_SD = json.loads(INFLUENZA_STUDY_DATA)
# Baked CPM 2D spatial-state frames (per-site cell-TYPE mosaic) per study.
_SP = json.loads(INFLUENZA_SPATIAL)


# ── CPM 2D spatial-state videos (the primary, animated mosaic per study) ─────
def _spatial_stamp(text):
    """A recessive MCS/day timestamp pinned top-left of the mosaic frame."""
    return {"text": text, "x": 0.014, "y": 0.986, "xref": "paper", "yref": "paper",
            "xanchor": "left", "yanchor": "top", "showarrow": False,
            "font": {"size": 12.5, "color": S.SECONDARY, "family": S.FONT},
            "bgcolor": "rgba(252,252,251,0.74)", "borderpad": 3}


def _spatial_legend(codes, *, owner=False):
    """Compact colour→cell-state legend (only the states present in a study)."""
    if owner:
        items = [("#cfe6fb", "each tile = one epithelial cell"),
                 (S.MEDIUM_TINT, "open medium")]
    else:
        items = [(S.SPATIAL_STATES[v][1], S.SPATIAL_STATES[v][0]) for v in codes]
    chips = "".join(
        f'<span style="display:inline-flex;align-items:center;gap:5px;'
        f'margin:0 15px 4px 0;font-size:11.5px;color:{S.SECONDARY}">'
        f'<span style="width:11px;height:11px;border-radius:3px;background:{col};'
        f'border:1px solid rgba(11,11,11,.14);flex:none"></span>{lab}</span>'
        for col, lab in items)
    return (f'<div style="display:flex;flex-wrap:wrap;align-items:center;'
            f'padding:6px 4px 4px">{chips}</div>')


def _spatial_axes(nx, ny):
    return {
        "xaxis": {"visible": False, "scaleanchor": "y", "constrain": "domain",
                  "range": [-0.5, nx - 0.5]},
        "yaxis": {"visible": False, "autorange": "reversed",
                  "range": [ny - 0.5, -0.5]},
    }


def _spatial_controls(steps):
    return {
        "updatemenus": [{"type": "buttons", "showactive": False, "x": 0.0, "y": -0.02,
                         "xanchor": "left", "yanchor": "top", "direction": "left",
                         "pad": {"t": 4}, "buttons": [
            {"label": "▶ Play", "method": "animate",
             "args": [None, {"mode": "immediate", "fromcurrent": True,
                             "frame": {"duration": 520, "redraw": True},
                             "transition": {"duration": 0}}]},
            {"label": "❙❙ Pause", "method": "animate",
             "args": [[None], {"mode": "immediate", "frame": {"duration": 0, "redraw": False},
                               "transition": {"duration": 0}}]}]}],
        "sliders": [{"active": 0, "x": 0.26, "len": 0.72, "y": 0.0,
                     "xanchor": "left", "yanchor": "top", "pad": {"t": 4},
                     "currentvalue": {"prefix": "MCS ", "xanchor": "left",
                                      "font": {"size": 12, "color": S.SECONDARY}},
                     "steps": steps}],
    }


def _decode_grid(frame, nx, ny, dtype):
    """Inflate one baked frame (zlib+base64 row-major int lattice) to a 2D
    numpy grid — the inverse of ``_capture_spatial._encode``."""
    raw = zlib.decompress(base64.b64decode(frame["grid"]))
    return np.frombuffer(raw, dtype=np.dtype(dtype)).reshape(ny, nx)


def _type_traces(grid2d):
    # code 8 = cell-boundary sentinel (baked tessellation outline) -> extend the
    # discrete scale to 0..8 so those sites draw dark; harmless when absent.
    return [{"type": "heatmap", "z": grid2d.tolist(),
             "colorscale": S.discrete_state_colorscale(8),
             "zmin": -0.5, "zmax": 8.5, "showscale": False, "xgap": 0, "ygap": 0,
             "hoverinfo": "skip"}]


def _owner_traces(grid2d):
    # hash each cell label into a repeating soft-hue tile band; medium recedes
    g = grid2d.astype(np.int64)
    zz = np.where(g == 0, 0, 1 + (g * 2654435761 % 11)).tolist()
    return [{"type": "heatmap", "z": zz, "colorscale": S.MOSAIC_SCALE, "zmin": 0,
             "zmax": 11, "showscale": False, "xgap": 0, "ygap": 0, "hoverinfo": "skip"}]


def _build_spatial(slug, title, subtitle, caption, *, increment, kpis,
                   day=True, day_stamp="day"):
    """Render one study's baked spatial-state frames as an animated Plotly
    heatmap (▶ Play/❙❙ Pause + MCS slider, MCS/day stamp, colour→state legend,
    equal-aspect gridless mosaic). ``kind='owner'`` studies (the substrate
    sheet) render the cell-label tile mosaic; all others the cell-TYPE mosaic."""
    d = _SP[slug]
    nx, ny, owner = d["nx"], d["ny"], d["kind"] == "owner"
    dtype = d.get("dtype", "uint8")
    frames_data = d["frames"]
    decoded = [_decode_grid(fr, nx, ny, dtype) for fr in frames_data]
    div = f"influenza-spatial-{slug}"

    def _stamp_text(mcs):
        if owner:
            return f"MCS {mcs} · relaxation"
        if day:
            return f"MCS {mcs} · day {mcs / 1440.0:.2f}"
        return f"MCS {mcs}"

    def _traces(g):
        return (_owner_traces if owner else _type_traces)(g)

    frames = [{"name": str(i), "data": _traces(g),
               "layout": {"annotations": [_spatial_stamp(_stamp_text(fr["mcs"]))]}}
              for i, (fr, g) in enumerate(zip(frames_data, decoded))]
    steps = [{"label": str(fr["mcs"]), "method": "animate",
              "args": [[str(i)], {"mode": "immediate",
                                  "frame": {"duration": 0, "redraw": True},
                                  "transition": {"duration": 0}}]}
             for i, fr in enumerate(frames_data)]

    layout = S.base_layout(height=self_height(nx, ny), showlegend=False, extra={
        "margin": {"l": 14, "r": 14, "t": 12, "b": 54},
        "annotations": [_spatial_stamp(_stamp_text(frames_data[0]["mcs"]))],
        **_spatial_axes(nx, ny), **_spatial_controls(steps)})

    # Exclude the boundary sentinel (8) from the state legend — it's an outline,
    # not a cell state.
    codes = (sorted({int(v) for g in decoded for v in np.unique(g) if int(v) != 8})
             if not owner else [])
    body = (_spatial_legend(codes, owner=owner)
            + S.plot(div, _traces(decoded[0]), layout, frames=frames))
    return S.card(title, subtitle, body, caption, increment=increment, kpis=kpis)


def self_height(nx, ny):
    """Pick a card height that keeps the equal-aspect mosaic compact but legible
    (taller lattices get a little more room), clamped to a tasteful band."""
    return int(max(360, min(540, 300 + 260 * (ny / max(nx, 1)))))


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


# ── Study D · ifn-resistance (Increment 3): IFN gate cuts virus load ────────
_VIRUS_BLUE = "#256abf"


def _build_ifn():
    """IFN->resistance gate roughly halves virus load; resist plateaus ~0.55"""
    d = _SD["ifn"]
    x = d["steps"]
    wi, wo = d["with"], d["without"]
    virus = [S.line_trace(x, wo["total_virus"], S.MUTED, "no IFN (control)",
                          dash="dash", markers=False),
             S.line_trace(x, wi["total_virus"], _VIRUS_BLUE, "IFN gate active",
                          markers=False)]
    virus_layout = S.base_layout("Update", "Total virus (Σ field)", height=250)
    ni = [S.line_trace(x, wo["n_I"], S.MUTED, "no IFN (control)", dash="dash", markers=False),
          S.line_trace(x, wi["n_I"], S.INFECTED, "IFN gate active", markers=False)]
    ni_layout = S.base_layout("Update", "Infected cells", height=230)
    resist = [S.line_trace(x, wi["mean_resist"], S.GOOD, "mean resistance ρ", markers=False)]
    resist_layout = S.base_layout("Update", "mean ρ (infected cells)", height=210,
                                  showlegend=False, extra={"yaxis": {"range": [0, 1]}})
    body = ('<div style="font-size:12px;color:%s;font-weight:600;margin:2px 4px 0">'
            'Virus load — the primary, robust signal (~halved every step)</div>' % S.SECONDARY
            + S.plot("influenza-ifn-virus", virus, virus_layout)
            + '<div style="font-size:12px;color:%s;font-weight:600;margin:8px 4px 0">'
              'Infected count — weaker early, widens to ~23%% by step 40</div>' % S.SECONDARY
            + S.plot("influenza-ifn-ni", ni, ni_layout)
            + '<div style="font-size:12px;color:%s;font-weight:600;margin:8px 4px 0">'
              'Per-cell resistance ρ plateaus fast (~0.55), it does NOT ramp to 1</div>' % S.SECONDARY
            + S.plot("influenza-ifn-resist", resist, resist_layout))
    v20w, v20o = wi["total_virus"][20], wo["total_virus"][20]
    kpis = [(f'{v20w:.0f} vs {v20o:.0f}', f'virus @ step 20 (−{100*(1-v20w/v20o):.0f}%)'),
            (f'~{d["with"]["mean_resist"][40]:.2f}', 'resist ρ plateau (not →1)'),
            (f'{wi["n_I"][40]} vs {wo["n_I"][40]}', f'infected @ step 40 (−{100*(1-wi["n_I"][40]/wo["n_I"][40]):.0f}%)')]
    cap = (f'<b>Mechanism validated, reproduction PENDING.</b> Same 0.3 mm / 900-cell sheet, '
           f'seed 17, as Increment 2 — but each infected cell\'s virus secretion is throttled by '
           f'<b>(1−ρ)</b>, with ρ from that cell\'s locally-sampled IFN. The <b>primary, robust</b> '
           f'effect is on virus load (~half the control from the earliest sampled step: '
           f'{v20w:.0f} vs {v20o:.0f} at step 20). The effect on infected COUNT is real but weak '
           f'early (~4% at step 20) and only becomes sizeable (~23%) by step 40. ρ jumps to ~0.53 '
           f'after one update and plateaus in a 0.55–0.57 band — a steady brake, not a hard block. '
           f'This is NOT a claim any Sego 2022 Fig 3B/5/7 target is met (deferred to Increment 9); '
           f'note also the source\'s IFN-self-secretion feedback is omitted here, so ρ runs somewhat '
           f'HIGH vs a fully source-faithful version.')
    return S.card("IFN-driven resistance gates virus release",
                  "Increment-3 mechanism — with-IFN vs no-IFN control (seed 17, 60 updates)",
                  body, cap, increment=3, kpis=kpis)


# ── Study E · epithelial-fate (Increment 4): H→I→D lifecycle + Allee ────────
def _build_fate():
    """H->I->D lesion forms; D->H Allee recovery fires; direct Allee death rare"""
    d = _SD["fate"]
    x = d["steps"]
    comp = [S.line_trace(x, d["n_H"], S.HEALTHY, "Healthy (H)", markers=False),
            S.line_trace(x, d["n_I"], S.INFECTED, "Infected (I)", markers=False),
            S.line_trace(x, d["n_D"], S.DEAD, "Dead (D)", markers=False)]
    comp_layout = S.base_layout("Update", "Cell count", height=300)
    events = [S.line_trace(x, d["recovery_cum"], S.GOOD, "Allee recovery D→H (cumulative)",
                          markers=False),
              S.line_trace(x, d["death_cum"], S.WARNING, "Allee death H→D (cumulative)",
                          markers=False)]
    events_layout = S.base_layout("Update", "Cumulative events", height=240)
    body = ('<div style="font-size:12px;color:%s;font-weight:600;margin:2px 4px 0">'
            'Epithelial-fate composition (population conserved every update)</div>' % S.SECONDARY
            + S.plot("influenza-fate-comp", comp, comp_layout)
            + '<div style="font-size:12px;color:%s;font-weight:600;margin:8px 4px 0">'
              'Allee-driven transitions accrue organically (no hand-made geometry)</div>' % S.SECONDARY
            + S.plot("influenza-fate-events", events, events_layout))
    kpis = [(f'0 → {d["n_D"][-1]}', 'dead cells (H→I→D lesion, 200 updates)'),
            (f'{d["recovery_cum"][-1]}', 'organic D→H recovery events'),
            (f'{d["death_cum"][-1]}', 'direct H→D Allee deaths (rare, not a bug)')]
    cap = (f'<b>Mechanism validated, reproduction PENDING.</b> The lifecycle composes: infection '
           f'seeds I cells, infected-death converts them to D (n_D 0→{d["n_D"][-1]} over 200 updates, '
           f'n_H {d["n_H"][0]}→{d["n_H"][-1]}, n_I {d["n_I"][0]}→{d["n_I"][-1]}, population conserved '
           f'every step). D→H Allee recovery fires organically on local contact geometry '
           f'({d["recovery_cum"][-1]} events here; a deterministic negative control confirms a '
           f'dead cell surrounded only by dying neighbours never recovers). Direct H→D Allee death '
           f'is genuinely rare at 0.3 mm ({d["death_cum"][-1]} events) — a source-literal b_h/θ '
           f'asymmetry plus this driver\'s per-update cadence, verified NOT a wiring bug (~0.19 '
           f'expected successes over the run). Dominant death path is infection, not Allee death. '
           f'No Fig 3B/5/7 target is claimed (deferred to Increment 9).')
    return S.card("Cellularized epithelial-fate lifecycle",
                  "Increment-4 mechanism — H→I→D + D→H Allee recovery (seed 17, 200 updates)",
                  body, cap, increment=4, kpis=kpis)


# ── Study F · macrophage-response (Increment 5): chemotaxis localization ─────
def _build_macrophage():
    """Macrophages chemotax up the virus field and localize to the infection"""
    d = _SD["macrophage"]
    x = d["steps"]
    start = d["start_distance"]
    dist = [S.line_trace(x, d["off"]["dist"], S.MUTED, "λ=0 control", dash="dash", markers=False),
            S.line_trace(x, d["on"]["dist"], S.CELL_STATES[4][1], "chemotaxis on (λ=5000)",
                        markers=False)]
    dist_layout = S.base_layout("Update", "Mean distance to infection (sites)", height=290)
    com = d["on"]["com"]
    xs = [c[0] for c in com]
    ys = [c[1] for c in com]
    traj = [{"x": xs, "y": ys, "type": "scatter", "mode": "lines+markers",
             "line": {"color": S.CELL_STATES[4][1], "width": 2.2},
             "marker": {"color": S.CELL_STATES[4][1], "size": 4},
             "name": "macrophage COM", "hovertemplate": "(%{x:.1f}, %{y:.1f})<extra></extra>"},
            {"x": [xs[0]], "y": [ys[0]], "type": "scatter", "mode": "markers",
             "marker": {"color": S.GOOD, "size": 11, "symbol": "circle",
                        "line": {"color": S.SURFACE, "width": 1.5}}, "name": "start"},
            {"x": [xs[-1]], "y": [ys[-1]], "type": "scatter", "mode": "markers",
             "marker": {"color": S.INFECTED, "size": 11, "symbol": "star",
                        "line": {"color": S.SURFACE, "width": 1.5}}, "name": "end"}]
    traj_layout = S.base_layout("x (sites)", "y (sites)", height=280)
    body = ('<div style="font-size:12px;color:%s;font-weight:600;margin:2px 4px 0">'
            'Distance to infection — chemotaxis on approaches, λ=0 control drifts</div>' % S.SECONDARY
            + S.plot("influenza-mac-dist", dist, dist_layout)
            + '<div style="font-size:12px;color:%s;font-weight:600;margin:8px 4px 0">'
              'Macrophage centre-of-mass trajectory (chemotaxis on)</div>' % S.SECONDARY
            + S.plot("influenza-mac-traj", traj, traj_layout))
    on_ch = d["on"]["dist"][-1] - start
    off_ch = d["off"]["dist"][-1] - start
    kpis = [(f'{on_ch:+.1f} sites', 'distance change, seed 17 (chemotaxis on)'),
            (f'{off_ch:+.1f} sites', 'λ=0 control (same seed/geometry)'),
            ('5 / 5 seeds', 'localize (mean −13.4 vs control −0.85)')]
    cap = (f'<b>Localization validated, reproduction PENDING.</b> Both clusters start in the domain '
           f'INTERIOR at the same {start:.1f}-site separation (an unbiased control — no wall pins the '
           f'λ=0 drift). Shown: seed 17. With chemotaxis on (λ=5000), macrophages close the gap '
           f'({on_ch:+.1f} sites); the λ=0 control barely moves ({off_ch:+.1f}). Across the full '
           f'{{5,17,23,42,100}} seed set the on case localizes in 5/5 seeds (mean −13.39) vs a '
           f'near-isotropic control (mean −0.85, spread 8.44) — the on/off gap exceeds control noise '
           f'in every seed. This is a localization-mechanism claim only; no Fig 3B/5/7 target is '
           f'evaluated (Increment 9).')
    return S.card("Macrophages localize to the infection",
                  "Increment-5 mechanism — chemotaxis on vs λ=0 (seed 17, 40 updates)",
                  body, cap, increment=5, kpis=kpis)


# ── Study G · signaling-fields (Increment 6): chemokine + IL-10 fields ──────
def _build_signaling():
    """Chemokine gradient centered on the macrophage cluster; IL-10 from both sources"""
    d = _SD["signaling"]
    r = d["radial"]
    bar = [{"x": r["centers"], "y": r["means"], "type": "bar",
            "marker": {"color": S.CELL_STATES[4][1], "line": {"width": 0}},
            "hovertemplate": "%{x:.0f} sites<br>%{y:.5f}<extra></extra>", "name": "chemokine"}]
    bar_layout = S.base_layout("Distance from macrophage centroid (sites)",
                               "Mean chemokine", height=280, showlegend=False)
    x = d["steps"]
    il10 = [S.line_trace(x, d["il10_macro"], "#0f766e", "IL-10 at macrophages", markers=False),
            S.line_trace(x, d["il10_uninf"], S.WARNING, "IL-10 at uninfected (H)", markers=False)]
    il10_layout = S.base_layout("Update", "mean IL-10 at cell", height=250)
    body = ('<div style="font-size:12px;color:%s;font-weight:600;margin:2px 4px 0">'
            'Chemokine radial profile — a gradient centered on the cluster, decaying outward</div>'
            % S.SECONDARY
            + S.plot("influenza-sig-radial", bar, bar_layout)
            + '<div style="font-size:12px;color:%s;font-weight:600;margin:8px 4px 0">'
              'IL-10 accrues from BOTH regulated sources (macrophage + uninfected)</div>' % S.SECONDARY
            + S.plot("influenza-sig-il10", il10, il10_layout))
    inner, outer = r["means"][0], r["means"][-1]
    kpis = [(f'{d["chemo_near"]/d["chemo_far"]:.1f}×', 'chemokine near/far (per-cell)'),
            (f'{inner/outer:.1f}×', 'radial decay (inner/outer bin)'),
            (f'{d["il10_macro"][-1]:.1e} / {d["il10_uninf"][-1]:.1e}', 'IL-10 macro / uninfected')]
    cap = (f'<b>Field mechanism documented, reproduction PENDING.</b> Macrophage-released chemokine '
           f'forms a gradient centered on the cluster: {d["chemo_near"]:.5f} at macrophage cells vs '
           f'{d["chemo_far"]:.5f} at the uninfected patch {d["params"]["separation_sites"]} sites of '
           f'open Medium away ({d["chemo_near"]/d["chemo_far"]:.1f}×), and a monotonic '
           f'{inner/outer:.1f}× decay across 6 raw-lattice radial bins ({inner:.5f}→{outer:.5f}, '
           f'cell-position-independent). IL-10 is positive and rising from both its Hill-regulated '
           f'sources. This validates the chemokine/IL-10 FIELD mechanism, NOT a Fig 2/3A '
           f'reproduction: sig_1 here is a documented STATIC stub (its dynamic value is resolved in '
           f'Increment 8), and global boundary secretion is deferred (Increment 9).')
    return S.card("Chemokine + IL-10 signaling fields",
                  "Increment-6 mechanism — radial gradient + IL-10 sources (seed 17, 15 updates)",
                  body, cap, increment=6, kpis=kpis)


# ── Study H · cytotoxic-killing (Increment 7): NK/CD8 localize + kill ───────
def _build_cytotoxic():
    """NK/CD8 localize robustly + kill on contact; end-to-end clearance still weak"""
    d = _SD["cytotoxic"]
    x = d["steps"]
    dist = [S.line_trace(x, d["nk_dist"], S.CELL_STATES[5][1], "NK distance", markers=False),
            S.line_trace(x, d["cd8_dist"], S.CELL_STATES[6][1], "CD8⁺ distance", markers=False)]
    dist_layout = S.base_layout("Update", "Mean distance to infection (sites)", height=290)
    c = d["close"]
    cx = c["steps"]
    inf = [S.line_trace(cx, c["n_infected_nokill"], S.MUTED, "killing disabled (control)",
                       dash="dash", markers=False),
           S.line_trace(cx, c["n_infected_kill"], S.INFECTED, "killing enabled", markers=False)]
    inf_layout = S.base_layout("Update", "Infected cells", height=250,
                               extra={"yaxis": {"range": [-0.15, 1.3], "dtick": 1}})
    body = ('<div style="font-size:12px;color:%s;font-weight:600;margin:2px 4px 0">'
            'Localization at DEFAULT scale — NK/CD8 close the gap but never reach contact</div>'
            % S.SECONDARY
            + S.plot("influenza-cyto-dist", dist, dist_layout)
            + '<div style="font-size:12px;color:%s;font-weight:600;margin:8px 4px 0">'
              'Killing WORKS in a close-contact test (1→0 at update 22)</div>' % S.SECONDARY
            + S.plot("influenza-cyto-inf", inf, inf_layout))
    kill_at = next((s for s, n in zip(cx, c["n_infected_kill"]) if n == 0), None)
    kpis = [(f'{d["cd8_dist"][0]:.0f}→{d["cd8_dist"][-1]:.0f}', 'CD8⁺ distance (default scale)'),
            (f'{d["nk_dist"][0]:.0f}→{d["nk_dist"][-1]:.0f}', 'NK distance (never contacts)'),
            (f'1→0 @ update {kill_at}', 'close-contact kill (capability proven)')]
    cap = (f'<b>Capability proven, end-to-end clearance WEAK/PENDING.</b> NK/CD8 chemotax up the '
           f'chemokine field and localize robustly (5/5 seeds: CD8 mean −35.99, NK −20.37 vs a '
           f'λ=0 control near 0). Contact-killing WORKS where cells actually touch: in a close-contact '
           f'scenario the infected cell is cleared (1→0 at update {kill_at}), while the '
           f'killing-disabled control stays at 1. BUT at the default full-scale scenario the NK/CD8 '
           f'cluster closes distance ({d["cd8_dist"][0]:.0f}→{d["cd8_dist"][-1]:.0f} for CD8, '
           f'{d["nk_dist"][0]:.0f}→{d["nk_dist"][-1]:.0f} for NK) yet never reaches contact, so '
           f'n_infected never drops. End-to-end cytotoxic CLEARANCE at scale is an Increment-9 '
           f'field-magnitude calibration gap (same root cause as the 100× chemotaxis-scale flag), '
           f'not a reproduction claim.')
    return S.card("NK/CD8 cytotoxic localization + killing",
                  "Increment-7 mechanism — default-scale localization + close-contact kill (seed 17)",
                  body, cap, increment=7, kpis=kpis)


# ── Study I · global-coupling (Increment 8): hybrid Price-2015 ODE ──────────
def _build_global():
    """Hybrid 10-species global ODE coupled to the CPM patch; magnitudes uncalibrated"""
    d = _SD["global"]
    x = d["mcs"]
    ode = d["ode"]
    species = [S.line_trace(x, ode["T"], S.INFECTED, "T (TNF)", markers=False),
               S.line_trace(x, ode["X"], S.WARNING, "X (ROS)", markers=False),
               S.line_trace(x, ode["A"], S.GOOD, "A (antibody)", markers=False),
               S.line_trace(x, ode["P"], _VIRUS_BLUE, "P (APC)", markers=False)]
    ode_layout = S.base_layout("MCS", "ODE state value", height=290)
    sig = [S.line_trace(x, d["sigma1"], "#0f766e", "dynamic sig_1 = a₁₁·T + a₁₂·D",
                       markers=False)]
    sig_layout = S.base_layout("MCS", "sig_1", height=250, showlegend=False,
                               extra={"yaxis": {"exponentformat": "e"}})
    body = ('<div style="font-size:12px;color:%s;font-weight:600;margin:2px 4px 0">'
            'Systemic Price-2015 ODE species — integrated once per MCS (scipy LSODA)</div>'
            % S.SECONDARY
            + S.plot("influenza-gc-ode", species, ode_layout)
            + '<div style="font-size:12px;color:%s;font-weight:600;margin:8px 4px 0">'
              'Dynamic sig_1 rises with TNF but collapses ~1e-6 vs the retired stub (2.44)</div>'
            % S.SECONDARY
            + S.plot("influenza-gc-sig", sig, sig_layout))
    kpis = [(f'{d["sigma1"][0]:.1e}→{d["sigma1"][-1]:.1e}', 'dynamic sig_1 (stub was 2.44)'),
            ('3 / 3 stubs', 'mechanisms resolved, source-faithful'),
            ('calibration', 'magnitudes PENDING (Increment 9)')]
    cap = (f'<b>Mechanism-resolved, calibration-pending.</b> The hybrid 10-species Price-2015 global '
           f'ODE integrates once per MCS and couples bidirectionally to the spatial CPM patch (the '
           f'spatialized species are read from the world as ODE inputs, never integrated). The three '
           f'carried stubs are RESOLVED as source-faithful mechanisms: dynamic sig_1 (a₁₁·T + a₁₂·D), '
           f'ODE-driven recruitment, and NK/CD8 nearby-killing. BUT at this reduced-scale patch '
           f'(η≈{d["params"]["eta"]:.1e}) the coupling MAGNITUDES are not yet calibrated: dynamic '
           f'sig_1 collapses to {d["sigma1"][0]:.1e}–{d["sigma1"][-1]:.1e} vs the retired static stub '
           f'(2.44), and recruitment / nearby-killing stay too weak to engage naturally. The '
           f'mechanisms are wired correctly; unit-scale calibration is deferred to Increment 9.')
    return S.card("Hybrid global ODE coupling (Price 2015)",
                  "Increment-8 mechanism — 10-species ODE + dynamic sig_1 (seed 3, side 30)",
                  body, cap, increment=8, kpis=kpis)


# ── Study J · repro-fig3b (Increment 9 CAPSTONE): full model vs Fig-3B ──────
def _repro_band_traces(obs_name, target_obs, ensemble, color, label):
    """Shared by the three Increment-9 repro cards: for one observable,
    the digitized target's [lo,hi] band + ODE-reference line (if present in
    `target_obs`) PLUS the model ensemble-mean line (if present in
    `ensemble`) — the overlay each repro card's honesty rests on (model vs
    the SAME band `run.repro_fig3b`/`repro_fig5`/`repro_fig7`'s own
    `band_eval` compares against, not a second/different computation)."""
    traces = []
    pts = target_obs.get(obs_name) or []
    if pts:
        t = [p["t_days"] for p in pts]
        lo = [p.get("lo", p["value"]) for p in pts]
        hi = [p.get("hi", p["value"]) for p in pts]
        val = [p["value"] for p in pts]
        traces += S.band_trace(t, lo, hi, color, label + " target")
        traces.append(S.line_trace(t, val, color, label + " target (ODE ref)",
                                    dash="dot", width=1.4, markers=False))
    series = ensemble.get(obs_name) or []
    if series:
        mt = [p[0] for p in series]
        mv = [p[1] for p in series]
        traces.append(S.line_trace(mt, mv, color, label + " (model)", markers=False))
    return traces


def _build_repro_fig3b():
    """CAPSTONE run_full_model ensemble vs digitized Fig-3B bands — reduced scale, 0/12 in-band"""
    d = _SD["repro_fig3b"]
    ensemble = d["ensemble"]
    target = _D["targets"]["fig3b"]["observables"]

    traces_counts = []
    for key, label, color in (("uninfected_cells", "Uninfected (H)", S.HEALTHY),
                               ("infected_cells", "Infected (I)", S.INFECTED),
                               ("dead_cells", "Dead (D)", S.DEAD)):
        traces_counts += _repro_band_traces(key, target, ensemble, color, label)
    counts_layout = S.base_layout("t (days)", "cell count", height=320)

    traces_virus = _repro_band_traces("extracellular_virus", target, ensemble,
                                       _VIRUS_BLUE, "Extracellular virus")
    virus_layout = S.base_layout("t (days)", "spatial-mean concentration", height=280)

    body = ('<div style="font-size:12px;color:%s;font-weight:600;margin:2px 4px 0">'
            'Cell counts — ensemble mean (solid) vs digitized Fig-3B band (shaded)</div>'
            % S.SECONDARY
            + S.plot("influenza-repro-fig3b-counts", traces_counts, counts_layout)
            + '<div style="font-size:12px;color:%s;font-weight:600;margin:8px 4px 0">'
              'Extracellular virus — ensemble mean vs Fig-3B band</div>' % S.SECONDARY
            + S.plot("influenza-repro-fig3b-virus", traces_virus, virus_layout))

    be = d["band_eval"]
    n_pass = sum(1 for v in be["observables"].values() if v["in_band"])
    n_tot = len(be["observables"])
    sim_days = d["steps"] * 7 * 60 / 86400
    kpis = [(f'{n_pass} / {n_tot} observables', "in-band at reduced scale"),
            (f'{d["replicas"]} reps · {d["cells_per_side"]}×{d["cells_per_side"]} cells',
             "reduced-scale config"),
            ("Phase B (Mac mini)", "paper-scale ensemble still required")]
    cap = (f'<b>Reduced-scale, band_eval.passed=False — NO reproduction claimed.</b> '
           f'`run.repro_fig3b` wires every Increment 0-8 mechanism into ONE per-MCS loop '
           f'(`run_full_model`) and evaluates the ensemble mean against the digitized Fig-3B '
           f'acceptance bands (`targets/fig3b.json`, `bands.evaluate_study`). At this '
           f'reduced-scale config ({d["replicas"]} replicas, {d["cells_per_side"]}×'
           f'{d["cells_per_side"]} cells, {d["steps"]} records, ~{sim_days:.2f} simulated '
           f'days) {n_pass} of {n_tot} observables land in-band at every digitized checkpoint '
           f'— expected at this small population/short-duration scale, not a mechanism '
           f'failure. The mechanisms are wired and source-faithful; field/coupling '
           f'MAGNITUDES remain calibration-pending. Only the paper-scale 50-replica, '
           f'35×35-cell, ~3.5-day ensemble (Mac-mini Phase-B follow-up) can set a '
           f'`reproduced` verdict — this card must not be read as one.')
    return S.card("Fig-3B reproduction — reduced-scale ensemble vs acceptance band",
                  "Increment-9 CAPSTONE — run_full_model ensemble vs targets/fig3b.json",
                  body, cap, increment=9, kpis=kpis)


def _build_repro_sweep(kind):
    """Shared by `_build_repro_fig5`/`_build_repro_fig7`: the dose-response
    curve (final ensemble-mean uninfected fraction vs the swept dose) plus
    the highest-tested-dose `uninfected_cells` overlay against ITS
    scenario-filtered Fig-5/Fig-7 band (never the raw multi-scenario band
    list — mirrors `run.py`'s `_fig5_target_subset`/`_fig7_target_subset`
    scenario-grouping guard)."""
    if kind == "fig5":
        d = _SD["repro_fig5"]
        doses, by = d["loads"], d["by_load"]
        tag_key, xlabel = "viral_load_multiplier", "initial viral load multiplier"
        target_fig, fig_label = "fig5", "Fig-5"
    else:
        d = _SD["repro_fig7"]
        doses, by = d["fracs"], d["by_frac"]
        tag_key, xlabel = "initial_infection_fraction", "initial infection fraction"
        target_fig, fig_label = "fig7", "Fig-7"

    final_fracs = [by[str(dose)]["uninfected_final_frac"] for dose in doses]
    dose_traces = [S.line_trace(doses, final_fracs, S.INFECTED,
                                "final uninfected fraction (model)")]
    dose_layout = S.base_layout(xlabel, "final ensemble-mean uninfected fraction",
                                logx=True, height=300)

    rep_dose = doses[-1]  # highest tested dose — the paper's most-severe scenario
    rep = by[str(rep_dose)]
    target_all = _D["targets"][target_fig]["observables"].get("uninfected_cells", [])
    matched = [o for o in target_all if o.get(tag_key) == rep_dose]
    rep_target = {"uninfected_cells": matched} if matched else {}
    band_traces = _repro_band_traces("uninfected_cells", rep_target, rep["ensemble"],
                                     S.HEALTHY, "Uninfected (H)")
    band_layout = S.base_layout("t (days)", "uninfected_cells (ensemble mean)", height=300)

    sim_days = d["steps"] * 7 * 60 / 86400
    body = ('<div style="font-size:12px;color:%s;font-weight:600;margin:2px 4px 0">'
            f'Dose-response — final uninfected fraction vs {xlabel} (log x)</div>'
            % S.SECONDARY
            + S.plot(f"influenza-repro-{target_fig}-dose", dose_traces, dose_layout)
            + '<div style="font-size:12px;color:%s;font-weight:600;margin:8px 4px 0">'
              f'Highest tested dose ({rep_dose:g}) vs its {fig_label} acceptance band</div>'
            % S.SECONDARY
            + S.plot(f"influenza-repro-{target_fig}-band", band_traces, band_layout))

    n_fail = sum(1 for dose in doses if not by[str(dose)]["band_eval"]["passed"])
    lt = d.get("lethal_threshold")
    lt_str = f'{lt:g}' if lt is not None else "none reached"
    kpis = [(f'{len(doses) - n_fail} / {len(doses)} doses', "band_eval.passed at reduced scale"),
            (f'{d["replicas"]} reps · {d["cells_per_side"]}×{d["cells_per_side"]} cells',
             "reduced-scale config"),
            (f'model lethal_threshold={lt_str}', "this driver's own reduced-scale run")]
    cap = (f'<b>Reduced-scale, band_eval.passed=False at every tested dose — NO reproduction '
           f'claimed.</b> `run.repro_{target_fig}` sweeps `run_full_model` over {fig_label}\'s '
           f'{xlabel} scenarios ({", ".join(f"{x:g}" for x in doses)}) and evaluates EACH dose '
           f'against ONLY that dose\'s {fig_label} band subset (scenario-grouping guard). At '
           f'this reduced-scale config ({d["replicas"]} replica, {d["cells_per_side"]}×'
           f'{d["cells_per_side"]} cells, {d["steps"]} records, ~{sim_days:.2f} simulated days), '
           f'the monotone dose-response DIRECTION holds (higher dose → not-more-surviving) but '
           f'per-dose band containment fails for all {len(doses)} tested doses — expected at '
           f'this scale. Mechanisms are wired and source-faithful; magnitudes are '
           f'calibration-pending. Only the paper-scale ensemble (Mac-mini Phase-B follow-up) '
           f'can set a `reproduced` verdict.')
    return S.card(f"{fig_label} reproduction — dose-response vs acceptance band",
                  f"Increment-9 CAPSTONE — run_full_model {xlabel} sweep vs targets/{target_fig}.json",
                  body, cap, increment=9, kpis=kpis)


def _build_repro_fig5():
    """CAPSTONE viral-load sweep vs digitized Fig-5 bands — reduced scale, dose-response direction holds"""
    return _build_repro_sweep("fig5")


def _build_repro_fig7():
    """CAPSTONE infection-fraction sweep vs digitized Fig-7 bands — reduced scale, dose-response direction holds"""
    return _build_repro_sweep("fig7")


# ── per-study CPM 2D spatial-state videos (the PRIMARY viz per study) ────────
def _sp_kpis(slug, *headline):
    """A consistent KPI row: mosaic scale + frame span (from the baked frames)
    + the study's own one-line spatial headline."""
    d = _SP[slug]
    fr = d["frames"]
    span = f'MCS {fr[0]["mcs"]}→{fr[-1]["mcs"]}'
    return [(f'{d["nx"]}×{d["ny"]} lattice', 'full-resolution mosaic'),
            (f'{len(fr)} frames', span), headline]


def _build_spatial_sheet():
    """Confluent epithelial sheet relaxing under the Potts temperature (owner mosaic)"""
    cap = ('<b>The substrate, in motion.</b> The confluent 0.3&nbsp;mm / 900-cell epithelial '
           'tiling (Increment&nbsp;1) relaxes under the Cellular-Potts surface/volume energy over '
           '300 MCS — each tile is one cell, and you can watch the irregular cell boundaries '
           'settle while the sheet stays gap-free. This is geometry only: no virus, infection or '
           'immune biology runs yet (those layer on in Increments 2+). Rendered at full lattice '
           'resolution; the mean cell volume and throughput are quantified in the companion sheet card.')
    return _build_spatial("epithelial-sheet-baseline",
                          "Confluent epithelial sheet — spatial relaxation",
                          "Increment-1 substrate — the cell mosaic settling over 300 MCS (seed 17)",
                          cap, increment=1,
                          kpis=_sp_kpis("epithelial-sheet-baseline",
                                        "900 cells", "confluent, gap-free tiling"))


def _build_spatial_virus():
    """H→I lesion spreading locally across the sheet as the virus field diffuses"""
    ser = _D["infection"]["series"]
    cap = (f'<b>Play it.</b> The seeded lesion of infected cells '
           f'(<span style="color:{S.INFECTED}">■</span>) grows across the healthy sheet '
           f'(<span style="color:{S.HEALTHY}">■</span>) as the extracellular virus field diffuses '
           f'and drives stochastic H→I transitions — and new infections land NEXT TO existing '
           f'ones (locality ratio {_D["infection"]["locality"]["null_ratio"]:.2f} &lt; '
           f'{_D["infection"]["locality"]["pass_threshold"]}), not scattered uniformly. '
           f'n_I grows {ser["n_I"][0]}→{ser["n_I"][-1]} with population conserved. '
           f'Base infection mechanism only — no IFN resistance or cell death yet.')
    return _build_spatial("virus-field-infection",
                          "Virus-driven infection spread — spatial state",
                          "Increment-2 mechanism — the H→I lesion mosaic over 60 updates (seed 17)",
                          cap, increment=2,
                          kpis=_sp_kpis("virus-field-infection",
                                        f'n_I {ser["n_I"][0]}→{ser["n_I"][-1]}', "local lesion growth"))


def _build_spatial_ifn():
    """Same lesion under the IFN→resistance gate — visibly fewer infected cells"""
    d = _SD["ifn"]
    cap = (f'<b>Mechanism validated, reproduction PENDING.</b> The same 0.3&nbsp;mm / 900-cell '
           f'sheet and seed as Increment&nbsp;2, but each infected cell\'s virus secretion is '
           f'throttled by <b>(1−ρ)</b> from its locally-sampled type-I IFN. The primary effect is '
           f'on virus load (~halved every step); the mosaic shows the weaker-but-real effect on '
           f'infected COUNT — the lesion stays more contained than the no-IFN control '
           f'(n_I {d["with"]["n_I"][40]} vs {d["without"]["n_I"][40]} at step 40). ρ plateaus near '
           f'~0.55, a steady brake, not a hard block. No Fig&nbsp;3B/5/7 target is claimed here.')
    return _build_spatial("ifn-resistance",
                          "IFN-gated infection — spatial state",
                          "Increment-3 mechanism — lesion under the resistance gate (seed 17)",
                          cap, increment=3,
                          kpis=_sp_kpis("ifn-resistance",
                                        f'n_I −{100*(1-d["with"]["n_I"][40]/d["without"]["n_I"][40]):.0f}% @40',
                                        "vs no-IFN control"))


def _build_spatial_fate():
    """Full epithelial lifecycle: H→I→D lesion forms, D→H recovery fires"""
    d = _SD["fate"]
    cap = (f'<b>Mechanism validated, reproduction PENDING.</b> The lifecycle composes spatially: '
           f'infection seeds I cells (<span style="color:{S.INFECTED}">■</span>), infected-death '
           f'converts them to a dead core (<span style="color:{S.DEAD}">■</span>, n_D 0→'
           f'{d["n_D"][-1]} over 200 updates), and Allee recovery returns some dead cells to '
           f'healthy on local contact ({d["recovery_cum"][-1]} organic D→H events). Population is '
           f'conserved every step; direct H→D Allee death stays genuinely rare '
           f'({d["death_cum"][-1]} events, not a bug). Dominant death path is infection. '
           f'No Fig&nbsp;3B/5/7 target is claimed (deferred to the capstone).')
    return _build_spatial("epithelial-fate",
                          "Epithelial fate lifecycle — spatial state",
                          "Increment-4 mechanism — H→I→D + D→H recovery mosaic (seed 17, 200 updates)",
                          cap, increment=4,
                          kpis=_sp_kpis("epithelial-fate",
                                        f'n_D 0→{d["n_D"][-1]}', "dead lesion core forms"))


def _build_spatial_macrophage():
    """Macrophages chemotax up the virus field and localize to the lesion"""
    d = _SD["macrophage"]
    on_ch = d["on"]["dist"][-1] - d["start_distance"]
    cap = (f'<b>Localization validated, reproduction PENDING.</b> Macrophages '
           f'(<span style="color:{S.CELL_STATES[4][1]}">■</span>) chemotax up the virus field '
           f'(λ=5000) and migrate toward the infected patch '
           f'(<span style="color:{S.INFECTED}">■</span>), closing the gap by {on_ch:+.1f} sites '
           f'from a {d["start_distance"]:.0f}-site interior start — while a λ=0 control barely '
           f'moves ({d["off"]["dist"][-1]-d["start_distance"]:+.1f}). Across 5 seeds the on case '
           f'localizes 5/5 (mean −13.4 vs control −0.85). A localization-mechanism claim only; '
           f'no Fig&nbsp;3B/5/7 target is evaluated here.')
    return _build_spatial("macrophage-response",
                          "Macrophage recruitment — spatial state",
                          "Increment-5 mechanism — macrophages localizing to the lesion (seed 17)",
                          cap, increment=5,
                          kpis=_sp_kpis("macrophage-response",
                                        f'{on_ch:+.1f} sites', "macrophage net approach"))


def _build_spatial_signaling():
    """Chemokine + IL-10 fields around the macrophage cluster (spatial context)"""
    d = _SD["signaling"]
    cap = (f'<b>Field mechanism documented, reproduction PENDING.</b> The spatial context for the '
           f'signaling fields: the macrophage cluster '
           f'(<span style="color:{S.CELL_STATES[4][1]}">■</span>) sits a fixed '
           f'{d["params"]["separation_sites"]} sites of open medium from the uninfected epithelial '
           f'patch (<span style="color:{S.HEALTHY}">■</span>). Macrophage-released chemokine forms '
           f'a gradient centered on the cluster ({d["chemo_near"]/d["chemo_far"]:.1f}× near/far, '
           f'~{d["radial"]["means"][0]/d["radial"]["means"][-1]:.1f}× radial decay), and IL-10 '
           f'accrues from both regulated sources. This validates the FIELD mechanism, not a '
           f'Fig&nbsp;2/3A reproduction (sig_1 is a documented static stub here, resolved in '
           f'Increment&nbsp;8).')
    return _build_spatial("signaling-fields",
                          "Signaling-field scene — spatial state",
                          "Increment-6 mechanism — macrophage cluster + epithelial patch (seed 17)",
                          cap, increment=6,
                          kpis=_sp_kpis("signaling-fields",
                                        f'{d["chemo_near"]/d["chemo_far"]:.1f}×', "chemokine near/far"))


def _build_spatial_cytotoxic():
    """NK + CD8 clusters chemotax toward the lesion; contact-killing where they touch"""
    d = _SD["cytotoxic"]
    cap = (f'<b>Capability proven, end-to-end clearance WEAK/PENDING.</b> The full cytotoxic scene: '
           f'NK (<span style="color:{S.CELL_STATES[5][1]}">■</span>) and CD8⁺ '
           f'(<span style="color:{S.CELL_STATES[6][1]}">■</span>) clusters chemotax up the '
           f'chemokine field toward the infected patch '
           f'(<span style="color:{S.INFECTED}">■</span>) and localize robustly (5/5 seeds: CD8 '
           f'mean −36, NK −20). Contact-killing WORKS where cells touch (a close-contact test '
           f'clears the infected cell 1→0). BUT at this default scale the NK/CD8 clusters close '
           f'distance yet never reach contact, so n_infected does not drop — an Increment-9 '
           f'field-magnitude calibration gap, NOT a reproduction claim.')
    return _build_spatial("cytotoxic-killing",
                          "Cytotoxic response — spatial state",
                          "Increment-7 mechanism — NK/CD8 localization at default scale (seed 17)",
                          cap, increment=7,
                          kpis=_sp_kpis("cytotoxic-killing",
                                        f'CD8 {d["cd8_dist"][0]:.0f}→{d["cd8_dist"][-1]:.0f}',
                                        "localize, no contact yet"))


def _build_spatial_global():
    """The full immune scene coupled to the systemic ODE (reduced-scale patch)"""
    d = _SD["global"]
    cap = (f'<b>Mechanism-resolved, calibration-pending.</b> The full spatial immune scene at the '
           f'reduced-scale side-30 patch, coupled bidirectionally to the hybrid 10-species '
           f'Price-2015 global ODE (integrated once per MCS). Epithelial H/I cells, the '
           f'macrophage/NK/CD8 clusters, and the '
           f'<span style="color:{S.RESERVE_TINT}">▨</span> dormant recruit-reserve pool (parked in '
           f'the medium margins, activated by ODE-driven inflow) all share one lattice. At this '
           f'scale (η≈{d["params"]["eta"]:.1e}) the coupling MAGNITUDES are uncalibrated — dynamic '
           f'sig_1 collapses ~1e-6 vs the retired stub (2.44) and recruitment stays too weak to '
           f'engage naturally. Mechanisms wired correctly; unit-scale calibration deferred to '
           f'Increment&nbsp;9.')
    return _build_spatial("global-coupling",
                          "Global immune coupling — spatial state",
                          "Increment-8 mechanism — full scene + ODE-driven reserve pool (seed 3, side 30)",
                          cap, increment=8, day=False,
                          kpis=_sp_kpis("global-coupling",
                                        "10-species ODE", "coupled per MCS"))


def _build_spatial_repro(fig_label, increment_note):
    """Shared full-model money-shot mosaic for the three repro-fig* studies —
    a full-SCALE, full-DURATION single realization (the qualitative sweep),
    distinct from the reduced-scale quantitative ensemble card."""
    cap = (f'<b>Qualitative full-scale sweep — NOT a quantitative reproduction.</b> The '
           f'money shot: <b>every</b> Increment 0–10 mechanism wired into ONE per-MCS loop '
           f'(<code>run_full_model</code>), rendered as the live cell mosaic at the paper\'s '
           f'shipped <b>0.3&nbsp;mm scale</b> — a 35×35 (1,225-cell) epithelial patch with '
           f'macrophage / NK / CD8⁺ clusters, run the full <b>~3.5 simulated days</b> '
           f'(5,040 MCS) at full lattice resolution. Play it: the seeded lesion of infected '
           f'(<span style="color:{S.INFECTED}">■</span>) cells spreads across the healthy '
           f'(<span style="color:{S.HEALTHY}">■</span>) sheet, and Increment-10 ROS-driven '
           f'death then sweeps the epithelium to a dead (<span style="color:{S.DEAD}">■</span>) '
           f'core — the green→orange→grey collapse Sego 2022 Fig 2/3 shows — while macrophages '
           f'(<span style="color:{S.CELL_STATES[4][1]}">■</span>), NK '
           f'(<span style="color:{S.CELL_STATES[5][1]}">■</span>) and CD8⁺ '
           f'(<span style="color:{S.CELL_STATES[6][1]}">■</span>) act on the same lattice. '
           f'This is a <b>single-seed realization</b> for visual realism, not the calibrated '
           f'50-replica ensemble: the companion {fig_label} time-course card evaluates the '
           f'quantitative acceptance bands (still at reduced scale, band_eval.passed=False), '
           f'and full quantitative {fig_label} calibration is ongoing (Mac-mini Phase-B). '
           f'{increment_note} Read this as the mechanism\'s spatial behaviour, not a '
           f'<code>reproduced</code> verdict.')
    return _build_spatial("repro-full-model",
                          f"Full model — spatial state ({fig_label})",
                          "Increment-9/10 CAPSTONE — full run_full_model mosaic, full scale + "
                          "~3.5 days, ROS death (seed 0)",
                          cap, increment=9,
                          kpis=_sp_kpis("repro-full-model",
                                        "green→orange→grey", "ROS death sweeps the sheet"))


def _build_spatial_repro_fig3b():
    return _build_spatial_repro("Fig-3B", "The Fig-3B time-course card is the companion readout.")


def _build_spatial_repro_fig5():
    return _build_spatial_repro("Fig-5", "The Fig-5 dose-response card is the companion readout.")


def _build_spatial_repro_fig7():
    return _build_spatial_repro("Fig-7", "The Fig-7 dose-response card is the companion readout.")


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


@as_visualization(inputs={"mcs": "list[float]"}, name="InfluenzaIfnResistance", demo={"mcs": [0.0]})
def update_influenza_ifn_resistance(state):
    """IFN->resistance gate roughly halves virus load; resist plateaus ~0.55"""
    return {"html": _build_ifn()}


@as_visualization(inputs={"mcs": "list[float]"}, name="InfluenzaEpithelialFate", demo={"mcs": [0.0]})
def update_influenza_epithelial_fate(state):
    """H->I->D lesion forms; D->H Allee recovery fires; direct Allee death rare"""
    return {"html": _build_fate()}


@as_visualization(inputs={"mcs": "list[float]"}, name="InfluenzaMacrophageResponse", demo={"mcs": [0.0]})
def update_influenza_macrophage_response(state):
    """Macrophages chemotax up the virus field and localize to the infection"""
    return {"html": _build_macrophage()}


@as_visualization(inputs={"mcs": "list[float]"}, name="InfluenzaSignalingFields", demo={"mcs": [0.0]})
def update_influenza_signaling_fields(state):
    """Chemokine gradient centered on the macrophage cluster; IL-10 from both sources"""
    return {"html": _build_signaling()}


@as_visualization(inputs={"mcs": "list[float]"}, name="InfluenzaCytotoxicKilling", demo={"mcs": [0.0]})
def update_influenza_cytotoxic_killing(state):
    """NK/CD8 localize robustly + kill on contact; end-to-end clearance still weak"""
    return {"html": _build_cytotoxic()}


@as_visualization(inputs={"mcs": "list[float]"}, name="InfluenzaGlobalCoupling", demo={"mcs": [0.0]})
def update_influenza_global_coupling(state):
    """Hybrid 10-species global ODE coupled to the CPM patch; magnitudes uncalibrated"""
    return {"html": _build_global()}


@as_visualization(inputs={"mcs": "list[float]"}, name="InfluenzaReproFig3B", demo={"mcs": [0.0]})
def update_influenza_repro_fig3b(state):
    """CAPSTONE run_full_model ensemble vs digitized Fig-3B bands — reduced scale, 0/12 in-band"""
    return {"html": _build_repro_fig3b()}


@as_visualization(inputs={"mcs": "list[float]"}, name="InfluenzaReproFig5", demo={"mcs": [0.0]})
def update_influenza_repro_fig5(state):
    """CAPSTONE viral-load sweep vs digitized Fig-5 bands — reduced scale, dose-response direction holds"""
    return {"html": _build_repro_fig5()}


@as_visualization(inputs={"mcs": "list[float]"}, name="InfluenzaReproFig7", demo={"mcs": [0.0]})
def update_influenza_repro_fig7(state):
    """CAPSTONE infection-fraction sweep vs digitized Fig-7 bands — reduced scale, dose-response direction holds"""
    return {"html": _build_repro_fig7()}


# ── zero-arg accessors (exercise outside the Step lifecycle: tests + file render) ──
def InfluenzaFig3BTargets():
    return _build_fig3b()


def InfluenzaEpithelialSheet():
    return _build_sheet()


def InfluenzaVirusFieldScene():
    return _build_scene()


def InfluenzaInfectionDynamics():
    return _build_dynamics()


def InfluenzaIfnResistance():
    return _build_ifn()


def InfluenzaEpithelialFate():
    return _build_fate()


def InfluenzaMacrophageResponse():
    return _build_macrophage()


def InfluenzaSignalingFields():
    return _build_signaling()


def InfluenzaCytotoxicKilling():
    return _build_cytotoxic()


def InfluenzaGlobalCoupling():
    return _build_global()


def InfluenzaReproFig3B():
    return _build_repro_fig3b()


def InfluenzaReproFig5():
    return _build_repro_fig5()


def InfluenzaReproFig7():
    return _build_repro_fig7()


# ── CPM 2D spatial-state videos — registry wrappers (PRIMARY per study) ──────
@as_visualization(inputs={"mcs": "list[float]"}, name="InfluenzaSpatialSheet", demo={"mcs": [0.0]})
def update_influenza_spatial_sheet(state):
    """Confluent epithelial sheet relaxing under the Potts temperature (owner mosaic)"""
    return {"html": _build_spatial_sheet()}


@as_visualization(inputs={"mcs": "list[float]"}, name="InfluenzaSpatialVirusField", demo={"mcs": [0.0]})
def update_influenza_spatial_virus(state):
    """H→I lesion spreading locally across the sheet as the virus field diffuses"""
    return {"html": _build_spatial_virus()}


@as_visualization(inputs={"mcs": "list[float]"}, name="InfluenzaSpatialIfnResistance", demo={"mcs": [0.0]})
def update_influenza_spatial_ifn(state):
    """Same lesion under the IFN→resistance gate — visibly fewer infected cells"""
    return {"html": _build_spatial_ifn()}


@as_visualization(inputs={"mcs": "list[float]"}, name="InfluenzaSpatialEpithelialFate", demo={"mcs": [0.0]})
def update_influenza_spatial_fate(state):
    """Full epithelial lifecycle: H→I→D lesion forms, D→H recovery fires"""
    return {"html": _build_spatial_fate()}


@as_visualization(inputs={"mcs": "list[float]"}, name="InfluenzaSpatialMacrophage", demo={"mcs": [0.0]})
def update_influenza_spatial_macrophage(state):
    """Macrophages chemotax up the virus field and localize to the lesion"""
    return {"html": _build_spatial_macrophage()}


@as_visualization(inputs={"mcs": "list[float]"}, name="InfluenzaSpatialSignaling", demo={"mcs": [0.0]})
def update_influenza_spatial_signaling(state):
    """Chemokine + IL-10 fields around the macrophage cluster (spatial context)"""
    return {"html": _build_spatial_signaling()}


@as_visualization(inputs={"mcs": "list[float]"}, name="InfluenzaSpatialCytotoxic", demo={"mcs": [0.0]})
def update_influenza_spatial_cytotoxic(state):
    """NK + CD8 clusters chemotax toward the lesion; contact-killing where they touch"""
    return {"html": _build_spatial_cytotoxic()}


@as_visualization(inputs={"mcs": "list[float]"}, name="InfluenzaSpatialGlobalCoupling", demo={"mcs": [0.0]})
def update_influenza_spatial_global(state):
    """The full immune scene coupled to the systemic ODE (reduced-scale patch)"""
    return {"html": _build_spatial_global()}


@as_visualization(inputs={"mcs": "list[float]"}, name="InfluenzaSpatialReproFig3B", demo={"mcs": [0.0]})
def update_influenza_spatial_repro_fig3b(state):
    """Full run_full_model cell mosaic — reduced scale, calibration-pending (Fig-3B study)"""
    return {"html": _build_spatial_repro_fig3b()}


@as_visualization(inputs={"mcs": "list[float]"}, name="InfluenzaSpatialReproFig5", demo={"mcs": [0.0]})
def update_influenza_spatial_repro_fig5(state):
    """Full run_full_model cell mosaic — reduced scale, calibration-pending (Fig-5 study)"""
    return {"html": _build_spatial_repro_fig5()}


@as_visualization(inputs={"mcs": "list[float]"}, name="InfluenzaSpatialReproFig7", demo={"mcs": [0.0]})
def update_influenza_spatial_repro_fig7(state):
    """Full run_full_model cell mosaic — reduced scale, calibration-pending (Fig-7 study)"""
    return {"html": _build_spatial_repro_fig7()}


# ── spatial zero-arg accessors (tests + file render) ────────────────────────
def InfluenzaSpatialSheet():
    return _build_spatial_sheet()


def InfluenzaSpatialVirusField():
    return _build_spatial_virus()


def InfluenzaSpatialIfnResistance():
    return _build_spatial_ifn()


def InfluenzaSpatialEpithelialFate():
    return _build_spatial_fate()


def InfluenzaSpatialMacrophage():
    return _build_spatial_macrophage()


def InfluenzaSpatialSignaling():
    return _build_spatial_signaling()


def InfluenzaSpatialCytotoxic():
    return _build_spatial_cytotoxic()


def InfluenzaSpatialGlobalCoupling():
    return _build_spatial_global()


def InfluenzaSpatialReproFig3B():
    return _build_spatial_repro_fig3b()


def InfluenzaSpatialReproFig5():
    return _build_spatial_repro_fig5()


def InfluenzaSpatialReproFig7():
    return _build_spatial_repro_fig7()
