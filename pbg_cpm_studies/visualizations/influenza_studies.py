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
from ._influenza_data import INFLUENZA_DATA, INFLUENZA_STUDY_DATA

_D = json.loads(INFLUENZA_DATA)
# Increments 3-8 study readouts (compact REAL engine series, no field frames).
_SD = json.loads(INFLUENZA_STUDY_DATA)


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
