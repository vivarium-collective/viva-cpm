"""Minimal in-package visualizations for the influenza-sego2022 capability
ladder (Increments 1-5): the epithelial-sheet baseline, the virus-field
infection mechanism, the IFN/resistance mechanism, the cellularized
epithelial-fate (infection/death/Allee recovery) mechanism, and the
macrophage-localization (chemotaxis) mechanism.

The dashboard workbench isn't the delivery path for these increments (no
server running in this worktree), so these are plain matplotlib functions
(Agg, no display). This module is a MINIMAL stub -- the polished viz system
lives under `pbg_cpm_studies/visualizations/` (owned by a peer session); do
not add anything there from here.
"""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")  # headless: no display backend required

import numpy as np
from matplotlib.figure import Figure


def sheet_snapshot_figure(world) -> Figure:
    """Render a two-panel figure for a CPM ``world`` (``cpm.schema.load_world``
    result): a type-colored lattice snapshot, and a histogram of per-cell
    volumes. Returns a ``matplotlib.figure.Figure`` (not shown/saved)."""
    nx, ny, nz = world.dims()
    labels = np.asarray(world.snapshot(), dtype=np.int64).reshape(nz, ny, nx)
    types = np.asarray(world.cell_types(), dtype=np.int64)
    volumes = np.asarray(world.cell_volumes(), dtype=np.int64)

    type_lattice = types[labels[0]]  # site -> cell type, z=0 slice

    fig = Figure(figsize=(9, 4.2))
    ax_snap, ax_hist = fig.subplots(1, 2)

    ax_snap.imshow(type_lattice, origin="lower", interpolation="nearest", cmap="tab20")
    ax_snap.set_title("Epithelial sheet — cell types")
    ax_snap.set_xlabel("x (sites)")
    ax_snap.set_ylabel("y (sites)")

    cell_volumes = volumes[1:]  # exclude the medium placeholder at index 0
    ax_hist.hist(cell_volumes, bins=20)
    ax_hist.set_title("Cell volume distribution")
    ax_hist.set_xlabel("volume (sites)")
    ax_hist.set_ylabel("cell count")

    fig.tight_layout()
    return fig


def virus_infection_figure(run_result: dict) -> Figure:
    """Render a two-panel figure for a `run.run_virus_infection(...)` result
    (Increment 2, Task 2.3/2.4): (a) infected/uninfected cell counts vs
    update index, (b) total virus-field concentration vs update index.
    Returns a `matplotlib.figure.Figure` (not shown/saved)."""
    steps = run_result["steps"]
    n_H = run_result["n_H"]
    n_I = run_result["n_I"]
    total_virus = run_result["total_virus"]

    fig = Figure(figsize=(9, 4.2))
    ax_counts, ax_virus = fig.subplots(1, 2)

    ax_counts.plot(steps, n_I, label="n_I (infected)", color="firebrick")
    ax_counts.plot(steps, n_H, label="n_H (uninfected)", color="steelblue")
    ax_counts.set_title("Cell counts vs update")
    ax_counts.set_xlabel("update index")
    ax_counts.set_ylabel("cell count")
    ax_counts.legend()

    ax_virus.plot(steps, total_virus, color="darkorange")
    ax_virus.set_title("Total virus-field concentration vs update")
    ax_virus.set_xlabel("update index")
    ax_virus.set_ylabel("sum(field concentration)")

    fig.tight_layout()
    return fig


def ifn_resistance_figure(with_result: dict, without_result: dict) -> Figure:
    """Render a three-panel MINIMAL mechanism figure (Increment 3, Task 3.4)
    comparing a `run.run_virus_infection_with_ifn(...)` result against a
    `run.run_virus_infection(...)` result run with the same driver params:
    (a) n_I with-IFN vs without over update index, (b) total_virus with-IFN
    vs without, (c) the with-IFN run's mean_resist trajectory (the per-cell
    resistance signal that gates virus secretion by (1-resist); the
    no-resistance run has no such series).

    This is a minimal in-package stub, NOT the polished viz system under
    `pbg_cpm_studies/visualizations/` (owned by a peer session). Returns a
    `matplotlib.figure.Figure` (not shown/saved)."""
    steps = with_result["steps"]

    fig = Figure(figsize=(13, 4.2))
    ax_nI, ax_virus, ax_resist = fig.subplots(1, 3)

    ax_nI.plot(steps, with_result["n_I"], label="n_I (with IFN)", color="firebrick")
    ax_nI.plot(steps, without_result["n_I"], label="n_I (without IFN)",
               color="firebrick", linestyle="--")
    ax_nI.set_title("Infected cell count vs update")
    ax_nI.set_xlabel("update index")
    ax_nI.set_ylabel("cell count")
    ax_nI.legend()

    ax_virus.plot(steps, with_result["total_virus"], label="total_virus (with IFN)",
                  color="darkorange")
    ax_virus.plot(steps, without_result["total_virus"], label="total_virus (without IFN)",
                  color="darkorange", linestyle="--")
    ax_virus.set_title("Total virus-field concentration vs update")
    ax_virus.set_xlabel("update index")
    ax_virus.set_ylabel("sum(field concentration)")
    ax_virus.legend()

    ax_resist.plot(steps, with_result["mean_resist"], color="seagreen")
    ax_resist.set_title("Mean per-cell resistance (with-IFN run)")
    ax_resist.set_xlabel("update index")
    ax_resist.set_ylabel("mean resist (infected cells)")

    fig.tight_layout()
    return fig


def epithelial_fate_figure(run_result: dict) -> Figure:
    """Render a two-panel MINIMAL mechanism figure (Increment 4, Task 4.4)
    for a `run.run_epithelial_fate(...)` result: (a) cell-type composition
    n_H/n_I/n_D vs update index (the H -> I -> D lifecycle + the growing
    dead lesion), (b) cumulative Allee-driven recovery (D -> H) and death
    (H -> D) event counts vs update index -- `run_epithelial_fate` records
    these as sparse PER-UPDATE counts (`n_allee_recovery`/`n_allee_death`,
    mostly 0, occasionally 1+), so the cumulative sum is what actually shows
    the handful of organic events accruing over a run (Task 4.3: ~3-12
    recovery events/run, ~0 death events at 0.3mm -- see the study's
    calibration-flag caveat; a bare per-update line would look almost empty).

    This is a minimal in-package stub, NOT the polished viz system under
    `pbg_cpm_studies/visualizations/` (owned by a peer session). Returns a
    `matplotlib.figure.Figure` (not shown/saved)."""
    steps = run_result["steps"]
    n_allee_death = run_result["n_allee_death"]
    n_allee_recovery = run_result["n_allee_recovery"]

    death_cum = [sum(n_allee_death[: i + 1]) for i in range(len(steps))]
    recovery_cum = [sum(n_allee_recovery[: i + 1]) for i in range(len(steps))]

    fig = Figure(figsize=(9, 4.2))
    ax_counts, ax_events = fig.subplots(1, 2)

    ax_counts.plot(steps, run_result["n_H"], label="n_H (uninfected)", color="steelblue")
    ax_counts.plot(steps, run_result["n_I"], label="n_I (infected)", color="firebrick")
    ax_counts.plot(steps, run_result["n_D"], label="n_D (dead)", color="dimgray")
    ax_counts.set_title("Epithelial-fate composition vs update")
    ax_counts.set_xlabel("update index")
    ax_counts.set_ylabel("cell count")
    ax_counts.legend()

    ax_events.plot(steps, recovery_cum, label="cumulative Allee recovery (D->H)",
                    color="seagreen")
    ax_events.plot(steps, death_cum, label="cumulative Allee death (H->D)",
                    color="darkorange")
    ax_events.set_title("Cumulative Allee-driven events vs update")
    ax_events.set_xlabel("update index")
    ax_events.set_ylabel("cumulative event count")
    ax_events.legend()

    fig.tight_layout()
    return fig


def macrophage_response_figure(run_result: dict, control_result: dict | None = None) -> Figure:
    """Render a two-panel MINIMAL mechanism figure (Increment 5, Task 5.2) for
    a `run.run_macrophage_response(...)` result (Task 5.1): (a) macrophages'
    mean centre-of-mass distance to the infection centroid vs update index --
    the falsifiable localization claim (measured across 5 seeds, interior-
    placed scenario: chemotaxis on reduces the distance in 5/5 seeds, mean
    change -13.39 sites, vs a near-isotropic lambda=0 control, mean change
    -0.85 sites -- task-5.1-report.md's fix-round-1 section); (b) the
    macrophage centre-of-mass trajectory (x vs y over the run, one seed),
    with start/end markers.

    If `control_result` (a second `run_macrophage_response(...)` result, e.g.
    the lambda=0 control) is also passed, panel (a) overlays both distance
    trajectories so the on/off contrast is visible in one figure -- this is
    the study's primary evidence, not two separate plots.

    This is a minimal in-package stub, NOT the polished viz system under
    `pbg_cpm_studies/visualizations/` (owned by a peer session). Returns a
    `matplotlib.figure.Figure` (not shown/saved)."""
    steps = run_result["steps"]
    dist = run_result["mean_distance_to_infection"]
    com = run_result["macrophage_com"]
    lam = run_result.get("params", {}).get("chemotaxis_v_macro")
    run_label = f"chemotaxis on (lambda={lam:g})" if lam is not None else "chemotaxis on"

    fig = Figure(figsize=(9, 4.2))
    ax_dist, ax_traj = fig.subplots(1, 2)

    ax_dist.plot(steps, dist, label=run_label, color="firebrick")
    if control_result is not None:
        ax_dist.plot(control_result["steps"], control_result["mean_distance_to_infection"],
                     label="lambda=0 control", color="dimgray", linestyle="--")
    ax_dist.set_title("Macrophage mean distance to infection vs update")
    ax_dist.set_xlabel("update index")
    ax_dist.set_ylabel("mean distance to infection centroid (sites)")
    ax_dist.legend()

    xs = [c[0] for c in com]
    ys = [c[1] for c in com]
    ax_traj.plot(xs, ys, color="steelblue", marker="o", markersize=3, label="macrophage COM")
    ax_traj.scatter([xs[0]], [ys[0]], color="seagreen", zorder=3, label="start")
    ax_traj.scatter([xs[-1]], [ys[-1]], color="firebrick", zorder=3, label="end")
    ax_traj.set_title("Macrophage centre-of-mass trajectory")
    ax_traj.set_xlabel("x (sites)")
    ax_traj.set_ylabel("y (sites)")
    ax_traj.legend()

    fig.tight_layout()
    return fig


def global_coupling_figure(run_result: dict) -> Figure:
    """Render a two-panel MINIMAL mechanism figure (Increment 8, Task 8.6)
    for a `run.run_global_coupling(...)` result (Tasks 8.1-8.5): the hybrid
    Price-2015 global ODE (the 10 systemic/integrated species `{NB, N, T, X,
    A, B, P, W, G, O}`, `price_ode.GlobalODE`) integrated once per MCS and
    coupled BIDIRECTIONALLY to the spatial CPM patch -- discrepancy #9's
    documented hybrid correction (only these 10 states are integrated; the
    spatialized species `{H, I, M, E, K, L, C, F, V, DH}` are read from the
    CPM world each MCS as ODE *inputs*, never integrated). This task's three
    Increment-5/6/7 stubs are RESOLVED here: sig_1 is now DYNAMIC (Task 8.3,
    `a_11*T + a_12*D`, ODE TNF + spatial dead count, discrepancy #11's
    saturating IL-10-inhibited Michaelis secretion form, replacing the static
    Increment-6 `params.il10.sig_1_stub`); recruitment is ODE-driven (Task
    8.4, asymmetric chemokine/APC Hill inflows, discrepancy #12: CD8+ is
    APC(`P`)-driven with NO homeostatic baseline, unlike macrophage/NK); and
    NK/CD8 cytotoxic killing now includes a well-mixed NEARBY-population term
    (Task 8.5) alongside the Increment-7 LOCAL contact term.

    Panel (a) -- systemic-species trajectories: TNF (`T`), ROS (`X`),
    antibody (`A`), and APC (`P`) vs MCS, from `run_result["ode"]` (the
    ACTUAL scipy-LSODA-integrated ODE state each MCS). Panel (b) -- spatial
    aggregates: infected (`I`), macrophage (`M`), NK (`K`), and CD8+ (`E`)
    cell counts vs MCS, from `run_result["spatial"]`, PLUS the Task-8.3
    dynamic `sig_1` series (`run_result["sigma1"]`) on a twin y-axis (its
    scale is orders of magnitude different from the cell counts -- see the
    caveat below).

    HONEST CAVEAT (central finding, NOT a quantitative Fig-3B reproduction --
    that is Increment 9): at this reduced-scale patch (tiny `eta` ~1.4e-4 by
    default), the coupling MAGNITUDES are not yet calibrated even though the
    mechanisms are wired correctly and source-faithful. Dynamic sig_1
    collapses the Michaelis secretion-scale factor to ~1e-5-1e-7 versus the
    static Increment-6 stub's ~0.5-0.975 (task-8.3-report.md); ODE-driven
    recruitment gives no integer population growth at default scale within a
    short run (baseline-dominated only, task-8.4-report.md); the nearby
    killing term is too weak to fire via natural recruitment buildup alone
    within a feasible step budget (task-8.5-report.md). This figure renders
    whatever series `run_result` contains -- it does not itself assert
    field-vs-ODE unit-scale calibration, which is deferred to Increment 9.

    This is a minimal in-package stub, NOT the polished viz system under
    `pbg_cpm_studies/visualizations/` (owned by a peer session). Returns a
    `matplotlib.figure.Figure` (not shown/saved)."""
    mcs = run_result["mcs"]
    ode = run_result["ode"]
    spatial = run_result["spatial"]
    sigma1 = run_result["sigma1"]

    fig = Figure(figsize=(11, 4.4))
    ax_ode, ax_spatial = fig.subplots(1, 2)

    for sp, label, color in (
        ("T", "T (TNF)", "firebrick"),
        ("X", "X (ROS)", "darkorange"),
        ("A", "A (antibody)", "seagreen"),
        ("P", "P (APC)", "steelblue"),
    ):
        ax_ode.plot(mcs, ode[sp], label=label, color=color)
    ax_ode.set_title("Systemic ODE species vs MCS (Price-2015 hybrid, once/MCS)")
    ax_ode.set_xlabel("MCS")
    ax_ode.set_ylabel("ODE state value")
    ax_ode.legend(fontsize="small")

    for sp, label, color in (
        ("I", "I (infected)", "firebrick"),
        ("M", "M (macrophage)", "steelblue"),
        ("K", "K (NK)", "darkorange"),
        ("E", "E (CD8+)", "seagreen"),
    ):
        ax_spatial.plot(mcs, spatial[sp], label=label, color=color)
    ax_spatial.set_title("Spatial cell counts + dynamic sig_1 vs MCS")
    ax_spatial.set_xlabel("MCS")
    ax_spatial.set_ylabel("cell count")
    ax_spatial.legend(fontsize="small", loc="upper left")

    ax_sigma1 = ax_spatial.twinx()
    ax_sigma1.plot(mcs, sigma1, label="sigma1 = a11*T + a12*D (dynamic, Task 8.3)",
                   color="dimgray", linestyle="--")
    ax_sigma1.set_ylabel("sigma1 (macrophage secretion-scale driver)")
    ax_sigma1.legend(fontsize="small", loc="upper right")

    fig.tight_layout()
    return fig


def signaling_fields_figure(run_result: dict) -> Figure:
    """Render a two-panel MINIMAL mechanism figure (Increment 6, Task 6.2)
    for a `run.run_macrophage_signaling(...)` result (Task 6.1): (a) the
    chemokine field's radial profile around the macrophage-cluster centroid
    (`run_result["chemo_radial_profile"]` -- bar chart over the bin centers,
    peak in the innermost/macrophage-centered bin, decaying outward: Task
    6.1 measured 0.00195 -> 0.00033 across 6 bins, a 5.9x monotonic drop);
    (b) IL-10 level over time at the macrophage vs uninfected source
    populations (`il10_at_macrophages`/`il10_at_uninfected`, both positive
    and increasing -- Task 6.1 measured 2.97e-4 / 3.71e-5 at the final step).

    This validates the chemokine/IL-10 FIELD mechanism (macrophage-released,
    gradient centered on the infection) -- NOT a Figs 2/3A-style
    reproduction (verdict stays `documented`/PENDING, Increment 9): `sig_1`
    is a documented STUB (TNF-dependent dynamic value deferred to Increment
    8), `g_1`/`g_2` use the 0.3mm-scenario `eta` approximation, and global
    boundary secretion + the IL-10 uniform initial condition are deferred
    (see `signaling.py`'s and `fields.il10_hill_constants`'s docstrings, and
    task-6.1-report.md).

    This is a minimal in-package stub, NOT the polished viz system under
    `pbg_cpm_studies/visualizations/` (owned by a peer session). Returns a
    `matplotlib.figure.Figure` (not shown/saved)."""
    radial = run_result["chemo_radial_profile"]
    bin_edges = radial["bin_edges"]
    bin_means = radial["bin_means"]
    bin_centers = [(bin_edges[i] + bin_edges[i + 1]) / 2.0 for i in range(len(bin_means))]
    bin_width = (bin_edges[1] - bin_edges[0]) if len(bin_edges) > 1 else 1.0

    steps = run_result["steps"]

    fig = Figure(figsize=(9, 4.2))
    ax_radial, ax_il10 = fig.subplots(1, 2)

    ax_radial.bar(bin_centers, bin_means, width=bin_width * 0.9, color="darkorchid")
    ax_radial.set_title("Chemokine radial profile (from macrophage-cluster centroid)")
    ax_radial.set_xlabel("distance from macrophage centroid (sites)")
    ax_radial.set_ylabel("mean chemokine concentration")

    ax_il10.plot(steps, run_result["il10_at_macrophages"], label="IL-10 at macrophages",
                 color="teal")
    ax_il10.plot(steps, run_result["il10_at_uninfected"], label="IL-10 at uninfected (H)",
                 color="goldenrod")
    ax_il10.set_title("IL-10 level over time, by source population")
    ax_il10.set_xlabel("update index")
    ax_il10.set_ylabel("mean IL-10 concentration at cell")
    ax_il10.legend()

    fig.tight_layout()
    return fig


def cytotoxic_killing_figure(run_result: dict, control_result: dict | None = None,
                              no_killing_result: dict | None = None) -> Figure:
    """Render a two-panel MINIMAL mechanism figure (Increment 7, Task 7.3)
    for a `run.run_cytotoxic_response(...)` result (Task 7.1 localization +
    Task 7.2 contact-killing): (a) NK + CD8 mean centre-of-mass distance to
    the infection centroid over time -- the localization claim (measured
    across 5 seeds: NK mean change -20.37 sites vs a lambda=0 control's
    -1.02, CD8 -35.99 vs -0.54, task-7.1-report.md); (b) infected-cell count
    (`n_infected`) over time -- the killing claim. A small close-contact
    scenario shows killing works (1 -> 0 at ~update 22, task-7.2-report.md),
    but at the DEFAULT full-scale scenario NK/CD8 do not reach + contact the
    infected cells within the run budget, so `n_infected` stays unchanged --
    the killing CAPABILITY is proven, end-to-end clearance at this scale is
    NOT (an Increment-9 field-magnitude-calibration gap, same root cause as
    Task 7.1's 100x chemotaxis-scale flag). This figure does not itself
    assert which regime `run_result` came from -- it renders whatever
    `n_infected`/distance series are handed to it.

    If `control_result` (a second `run_cytotoxic_response(...)` result, e.g.
    the lambda=0 NK/CD8-chemotaxis control) is also passed, panel (a)
    overlays both NK+CD8 distance trajectories so the localization on/off
    contrast is visible in one figure. If `no_killing_result` (an
    `enable_killing=False` result, same scenario/seed otherwise) is also
    passed, panel (b) overlays both `n_infected` trajectories so the
    with-vs-without-killing contrast is visible in one figure -- together
    the study's two falsifiable claims, not two separate plots.

    This is a minimal in-package stub, NOT the polished viz system under
    `pbg_cpm_studies/visualizations/` (owned by a peer session). Returns a
    `matplotlib.figure.Figure` (not shown/saved)."""
    steps = run_result["steps"]

    fig = Figure(figsize=(9, 4.2))
    ax_dist, ax_infected = fig.subplots(1, 2)

    ax_dist.plot(steps, run_result["nk_mean_distance_to_infection"],
                 label="NK distance (chemotaxis on)", color="firebrick")
    ax_dist.plot(steps, run_result["cd8_mean_distance_to_infection"],
                 label="CD8 distance (chemotaxis on)", color="darkorange")
    if control_result is not None:
        ax_dist.plot(control_result["steps"], control_result["nk_mean_distance_to_infection"],
                     label="NK distance (lambda=0 control)", color="firebrick", linestyle="--")
        ax_dist.plot(control_result["steps"], control_result["cd8_mean_distance_to_infection"],
                     label="CD8 distance (lambda=0 control)", color="darkorange", linestyle="--")
    ax_dist.set_title("NK + CD8 mean distance to infection vs update")
    ax_dist.set_xlabel("update index")
    ax_dist.set_ylabel("mean distance to infection centroid (sites)")
    ax_dist.legend(fontsize="small")

    ax_infected.plot(steps, run_result["n_infected"], label="killing enabled",
                      color="firebrick", marker="o", markersize=3)
    if no_killing_result is not None:
        ax_infected.plot(no_killing_result["steps"], no_killing_result["n_infected"],
                          label="killing disabled (control)", color="dimgray",
                          linestyle="--", marker="o", markersize=3)
    ax_infected.set_title("Infected-cell count vs update")
    ax_infected.set_xlabel("update index")
    ax_infected.set_ylabel("n_infected")
    ax_infected.legend()

    fig.tight_layout()
    return fig
