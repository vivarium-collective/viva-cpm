import matplotlib
matplotlib.use("Agg")  # headless: no display backend required

from matplotlib.figure import Figure

from viva_cpm.schema import load_world
from viva_cpm_studies.influenza import sheet, viz
from viva_cpm_studies.influenza.run import (
    run_virus_infection, run_virus_infection_with_ifn, run_global_coupling,
    repro_fig3b, repro_fig5, repro_fig7,
)


def test_sheet_snapshot_figure_returns_figure_with_axes():
    spec = sheet.build_sheet_spec(patch_mm=0.1)  # 50x50 sites, 100 cells
    world = load_world(spec)
    world.step(1)  # a bit of dynamics so the snapshot isn't a trivial grid

    fig = viz.sheet_snapshot_figure(world)

    assert isinstance(fig, Figure)
    assert len(fig.axes) == 2
    ax_snap, ax_hist = fig.axes
    assert len(ax_snap.images) == 1        # the lattice snapshot
    assert len(ax_hist.patches) > 0        # the volume histogram bars


def test_virus_infection_figure_returns_figure_with_axes():
    result = run_virus_infection(patch_mm=0.3, steps=20)

    fig = viz.virus_infection_figure(result)

    assert isinstance(fig, Figure)
    assert len(fig.axes) == 2
    ax_counts, ax_virus = fig.axes

    # panel (a): infected/uninfected counts over update index -- non-vacuous
    assert len(ax_counts.lines) >= 1
    n_I_line = ax_counts.lines[0]
    xdata, ydata = n_I_line.get_data()
    assert len(xdata) == len(result["steps"])
    assert list(xdata) == result["steps"]
    assert list(ydata) == result["n_I"]

    # panel (b): total virus over update index -- non-vacuous
    assert len(ax_virus.lines) >= 1
    virus_line = ax_virus.lines[0]
    vx, vy = virus_line.get_data()
    assert len(vx) == len(result["steps"])
    assert list(vy) == result["total_virus"]


def _small_epithelial_fate_result():
    """Hand-built `run.run_epithelial_fate(...)`-shaped result (small,
    deterministic, no actual CPM run) -- exercises the viz against the
    driver's actual return keys without paying for a real sim in this test."""
    return {
        "steps": [0, 1, 2, 3, 4, 5],
        "n_H": [855, 850, 842, 830, 820, 812],
        "n_I": [45, 48, 52, 58, 60, 58],
        "n_D": [0, 2, 6, 12, 20, 30],
        "total_virus": [0.0, 12.3, 40.1, 88.7, 150.2, 210.9],
        "n_allee_death": [0, 0, 0, 0, 0, 0],
        "n_allee_recovery": [0, 0, 1, 0, 1, 0],
    }


def test_epithelial_fate_figure_returns_figure_with_axes():
    """Task 4.4: mechanism figure -- cell-type composition (n_H/n_I/n_D) over
    time, plus the Allee recovery/death event counts `run_epithelial_fate`
    records. Non-vacuous: every series must actually be plotted with the
    same data the driver returned."""
    result = _small_epithelial_fate_result()

    fig = viz.epithelial_fate_figure(result)

    assert isinstance(fig, Figure)
    assert len(fig.axes) == 2
    ax_counts, ax_events = fig.axes

    # panel (a): n_H/n_I/n_D vs update index -- three non-vacuous lines
    assert len(ax_counts.lines) >= 3
    for line in ax_counts.lines:
        xdata, ydata = line.get_data()
        assert len(xdata) == len(result["steps"])
        assert list(xdata) == result["steps"]
    ydata_by_line = [list(line.get_data()[1]) for line in ax_counts.lines]
    assert result["n_H"] in ydata_by_line
    assert result["n_I"] in ydata_by_line
    assert result["n_D"] in ydata_by_line

    # panel (b): cumulative Allee recovery/death event counts -- non-vacuous,
    # and reflects the (sparse) per-update event counts actually recorded.
    assert len(ax_events.lines) >= 1
    for line in ax_events.lines:
        xdata, ydata = line.get_data()
        assert len(xdata) == len(result["steps"])
    recovery_cum = [sum(result["n_allee_recovery"][: i + 1])
                     for i in range(len(result["steps"]))]
    event_ydata = [list(line.get_data()[1]) for line in ax_events.lines]
    assert recovery_cum in event_ydata
    assert max(recovery_cum) > 0  # sanity: the hand-built result has events


def _small_macrophage_result(chemotaxis_v_macro=5000.0):
    """Hand-built `run.run_macrophage_response(...)`-shaped result (Task 5.1)
    -- exercises the viz against the driver's actual return keys
    ("steps", "mean_distance_to_infection", "macrophage_com", "params")
    without paying for a real CPM+field run in this test."""
    return {
        "steps": [0, 1, 2, 3],
        "mean_distance_to_infection": [47.20, 43.0, 38.5, 33.8],
        "macrophage_com": [(25.0, 25.0), (27.2, 26.8), (29.5, 28.9), (31.8, 30.5)],
        "params": {
            "chemotaxis_v_macro": chemotaxis_v_macro, "n_macrophages": 6,
            "n_infected": 1, "epithelial_cells_per_side": 4,
            "margin_sites": 30, "separation_sites": 25, "seed": 17, "steps": 3,
            "mcs_per_update": 10, "field_warmup": 3000,
        },
    }


def _small_macrophage_control_result():
    """Same shape, lambda=0 interior control -- distance stays roughly flat
    (near-isotropic) instead of closing. Illustrative synthetic values only
    (not tied to a specific run), exercising the plotting code."""
    result = _small_macrophage_result(chemotaxis_v_macro=0.0)
    result["mean_distance_to_infection"] = [47.20, 47.5, 46.8, 47.0]
    result["macrophage_com"] = [(25.0, 25.0), (24.6, 25.4), (25.3, 24.7), (24.9, 25.2)]
    return result


def test_macrophage_response_figure_returns_figure_with_axes():
    """Task 5.2: mechanism figure -- macrophage mean-distance-to-infection
    over time (panel a) and the macrophage centre-of-mass trajectory
    (panel b). Non-vacuous: plotted data must match the driver's actual
    return series."""
    result = _small_macrophage_result()

    fig = viz.macrophage_response_figure(result)

    assert isinstance(fig, Figure)
    assert len(fig.axes) == 2
    ax_dist, ax_traj = fig.axes

    # panel (a): mean_distance_to_infection vs update index
    assert len(ax_dist.lines) >= 1
    dist_line = ax_dist.lines[0]
    xdata, ydata = dist_line.get_data()
    assert len(xdata) == len(result["steps"])
    assert list(xdata) == result["steps"]
    assert list(ydata) == result["mean_distance_to_infection"]

    # panel (b): macrophage centre-of-mass trajectory (x vs y)
    assert len(ax_traj.lines) >= 1
    tx, ty = ax_traj.lines[0].get_data()
    assert len(tx) == len(result["macrophage_com"])
    assert list(tx) == [c[0] for c in result["macrophage_com"]]
    assert list(ty) == [c[1] for c in result["macrophage_com"]]


def test_macrophage_response_figure_with_control_overlays_second_distance_line():
    """When a lambda=0 control result is also passed, panel (a) overlays
    both trajectories so the on/off contrast (Task 5.1's falsifiable claim)
    is visible in one figure."""
    result = _small_macrophage_result()
    control = _small_macrophage_control_result()

    fig = viz.macrophage_response_figure(result, control_result=control)

    ax_dist, _ax_traj = fig.axes
    assert len(ax_dist.lines) >= 2
    ydata_by_line = [list(line.get_data()[1]) for line in ax_dist.lines]
    assert result["mean_distance_to_infection"] in ydata_by_line
    assert control["mean_distance_to_infection"] in ydata_by_line


def test_ifn_resistance_figure_returns_figure_with_axes():
    """Task 3.4: mechanism figure -- n_I and total_virus, with-IFN vs
    without, plus the mean_resist trajectory (only `run_virus_infection_with_ifn`
    records `mean_resist`). Small runs (0.1mm, 5 steps) -- this test only
    checks the figure is non-vacuous, not the mechanism's magnitude (that's
    Task 3.3's report)."""
    common = dict(patch_mm=0.1, steps=5, seed=17, init_infected_frac=0.05, mcs_per_update=10)
    without_result = run_virus_infection(**common)
    with_result = run_virus_infection_with_ifn(**common)

    fig = viz.ifn_resistance_figure(with_result, without_result)

    assert isinstance(fig, Figure)
    assert len(fig.axes) == 3
    ax_nI, ax_virus, ax_resist = fig.axes

    # panel (a): n_I with-IFN vs without -- two non-vacuous lines
    assert len(ax_nI.lines) >= 2
    for line in ax_nI.lines:
        xdata, ydata = line.get_data()
        assert len(xdata) == len(with_result["steps"])

    # panel (b): total_virus with-IFN vs without -- two non-vacuous lines
    assert len(ax_virus.lines) >= 2
    for line in ax_virus.lines:
        xdata, ydata = line.get_data()
        assert len(xdata) == len(with_result["steps"])

    # panel (c): mean_resist trajectory (with-IFN run only)
    assert len(ax_resist.lines) >= 1
    rx, ry = ax_resist.lines[0].get_data()
    assert len(rx) == len(with_result["steps"])
    assert list(ry) == with_result["mean_resist"]


def _small_signaling_result():
    """Hand-built `run.run_macrophage_signaling(...)`-shaped result (Task
    6.1) -- exercises the viz against the driver's ACTUAL return keys
    ("steps", "total_chemokine", "total_il10", "chemo_at_macrophages",
    "chemo_at_uninfected", "il10_at_macrophages", "il10_at_uninfected",
    "chemo_radial_profile", "params") without paying for a real CPM+field
    run in this test. Shape/pattern (near > far, gradient decays with
    distance) matches task-6.1-report.md's measured numbers, but the exact
    values here are illustrative, not re-asserted."""
    return {
        "steps": [0, 1, 2, 3],
        "total_chemokine": [0.0, 2.1, 5.4, 9.17],
        "total_il10": [0.0, 0.15, 0.4, 0.72],
        "chemo_at_macrophages": [0.0, 0.0009, 0.0018, 0.00261],
        "chemo_at_uninfected": [0.0, 0.0002, 0.00045, 0.00062],
        "il10_at_macrophages": [0.0, 1.0e-4, 2.0e-4, 2.97e-4],
        "il10_at_uninfected": [0.0, 1.2e-5, 2.5e-5, 3.71e-5],
        "chemo_radial_profile": {
            "bin_edges": [0.0, 15.5, 31.0, 46.4, 61.9, 77.4, 92.9],
            "bin_means": [0.00195, 0.00127, 0.00098, 0.00060, 0.00038, 0.00033],
        },
        "params": {
            "chemotaxis_v_macro": 5000.0, "n_macrophages": 6, "n_infected": 1,
            "epithelial_cells_per_side": 4, "margin_sites": 30,
            "separation_sites": 25, "seed": 17, "steps": 3,
            "mcs_per_update": 10, "field_warmup": 0,
        },
    }


def test_signaling_fields_figure_returns_figure_with_axes():
    """Task 6.2: mechanism figure -- chemokine field radial profile centered
    on the macrophage cluster (panel a: peak-near-macrophages, decaying with
    distance), and IL-10 level over time at both source populations (panel
    b). Non-vacuous: plotted data must match the driver's actual return
    series/keys (`run.run_macrophage_signaling`, Task 6.1)."""
    result = _small_signaling_result()

    fig = viz.signaling_fields_figure(result)

    assert isinstance(fig, Figure)
    assert len(fig.axes) == 2
    ax_radial, ax_il10 = fig.axes

    # panel (a): chemokine radial profile -- bin_means vs bin centers,
    # non-vacuous, and actually decaying (peak near the macrophage cluster).
    bin_means = result["chemo_radial_profile"]["bin_means"]
    assert len(ax_radial.lines) >= 1 or len(ax_radial.patches) >= 1
    if ax_radial.lines:
        plotted = [list(line.get_data()[1]) for line in ax_radial.lines]
        assert bin_means in plotted
    else:
        heights = [p.get_height() for p in ax_radial.patches]
        assert list(heights) == bin_means
    assert bin_means[0] > bin_means[-1]  # peak near macrophages, decaying

    # panel (b): IL-10 over time -- at least the macrophage + uninfected
    # source trajectories, matching the driver's actual series.
    assert len(ax_il10.lines) >= 2
    ydata_by_line = [list(line.get_data()[1]) for line in ax_il10.lines]
    assert result["il10_at_macrophages"] in ydata_by_line
    assert result["il10_at_uninfected"] in ydata_by_line


def _small_cytotoxic_result(enable_killing=True, chemotaxis_on=True, n_infected_series=None):
    """Hand-built `run.run_cytotoxic_response(...)`-shaped result (Task 7.1
    localization + Task 7.2 killing) -- exercises the viz against the
    driver's ACTUAL return keys ("steps", "nk_mean_distance_to_infection",
    "cd8_mean_distance_to_infection", "total_chemokine", "n_infected",
    "params") without paying for a real CPM+field run in this test. Shape
    (chemotaxis-on closes distance, killing-enabled reduces n_infected)
    is illustrative only, matching task-7.1/7.2-report.md's qualitative
    direction, not re-asserted here."""
    if n_infected_series is None:
        n_infected_series = [1, 1, 1, 0, 0] if enable_killing else [1, 1, 1, 1, 1]
    if chemotaxis_on:
        nk_dist = [81.6, 74.2, 66.8, 61.9, 57.8]
        cd8_dist = [81.6, 68.4, 55.9, 46.1, 39.8]
        lam_nk, lam_cd8 = 500000.0, 1000000.0
    else:
        nk_dist = [81.6, 82.1, 80.9, 81.7, 81.3]
        cd8_dist = [81.6, 81.2, 82.0, 81.5, 81.9]
        lam_nk, lam_cd8 = 0.0, 0.0
    return {
        "steps": [0, 1, 2, 3, 4],
        "nk_mean_distance_to_infection": nk_dist,
        "cd8_mean_distance_to_infection": cd8_dist,
        "total_chemokine": [0.0, 1.2, 2.9, 4.7, 6.5],
        "n_infected": n_infected_series,
        "params": {
            "chemotaxis_v_macro": 5000.0, "chemotaxis_v_nk": lam_nk, "chemotaxis_v_cd8": lam_cd8,
            "enable_killing": enable_killing, "g_ik": 0.0001, "g_ie": 0.0001, "tot_ec_ODE": 250000.0,
            "n_macrophages": 6, "n_nk": 6, "n_cd8": 6, "n_infected": 1,
            "epithelial_cells_per_side": 4, "margin_sites": 20, "separation_sites": 20,
            "seed": 17, "steps": 4, "mcs_per_update": 10, "field_warmup": 0,
        },
    }


def test_cytotoxic_killing_figure_returns_figure_with_axes():
    """Task 7.3: mechanism figure -- (a) NK + CD8 mean distance-to-infection
    over time (chemotaxis on), and (b) infected-cell count over time,
    killing enabled -- non-vacuous, using `run.run_cytotoxic_response`'s
    actual return keys."""
    result = _small_cytotoxic_result(enable_killing=True, chemotaxis_on=True)

    fig = viz.cytotoxic_killing_figure(result)

    assert isinstance(fig, Figure)
    assert len(fig.axes) == 2
    ax_dist, ax_infected = fig.axes

    # panel (a): NK + CD8 distance-to-infection -- two non-vacuous lines
    # matching the driver's own series.
    assert len(ax_dist.lines) >= 2
    for line in ax_dist.lines:
        xdata, ydata = line.get_data()
        assert len(xdata) == len(result["steps"])
        assert list(xdata) == result["steps"]
    ydata_by_line = [list(line.get_data()[1]) for line in ax_dist.lines]
    assert result["nk_mean_distance_to_infection"] in ydata_by_line
    assert result["cd8_mean_distance_to_infection"] in ydata_by_line

    # panel (b): n_infected over time -- non-vacuous, matches the driver's
    # own series.
    assert len(ax_infected.lines) >= 1
    inf_ydata_by_line = [list(line.get_data()[1]) for line in ax_infected.lines]
    assert result["n_infected"] in inf_ydata_by_line


def test_cytotoxic_killing_figure_overlays_control_and_no_killing_results():
    """When a chemotaxis-off control (localization contrast, Task 7.1) and a
    killing-disabled result (killing contrast, Task 7.2) are also passed,
    both panels overlay the on/off comparison so the study's two
    falsifiable claims -- localization AND killing -- are visible in one
    figure, plus the honest weak-end-to-end-clearance gap (killing-enabled
    n_infected does not always reach 0 at default/full scale)."""
    on_result = _small_cytotoxic_result(enable_killing=True, chemotaxis_on=True)
    control_result = _small_cytotoxic_result(enable_killing=False, chemotaxis_on=False)
    no_killing_result = _small_cytotoxic_result(enable_killing=False, chemotaxis_on=True)

    fig = viz.cytotoxic_killing_figure(
        on_result, control_result=control_result, no_killing_result=no_killing_result)

    ax_dist, ax_infected = fig.axes

    # panel (a): on's NK/CD8 distance lines AND the control's, all present.
    assert len(ax_dist.lines) >= 4
    ydata_by_line = [list(line.get_data()[1]) for line in ax_dist.lines]
    assert on_result["nk_mean_distance_to_infection"] in ydata_by_line
    assert on_result["cd8_mean_distance_to_infection"] in ydata_by_line
    assert control_result["nk_mean_distance_to_infection"] in ydata_by_line
    assert control_result["cd8_mean_distance_to_infection"] in ydata_by_line

    # panel (b): killing-enabled n_infected AND killing-disabled n_infected,
    # both present -- the with-vs-without-killing contrast.
    assert len(ax_infected.lines) >= 2
    inf_ydata_by_line = [list(line.get_data()[1]) for line in ax_infected.lines]
    assert on_result["n_infected"] in inf_ydata_by_line
    assert no_killing_result["n_infected"] in inf_ydata_by_line


def test_global_coupling_figure_returns_figure_with_two_axes():
    """Task 8.6: mechanism figure for `run.run_global_coupling(...)` (Tasks
    8.1-8.5) -- the hybrid Price-2015 global ODE (10 integrated systemic
    species) coupled bidirectionally to the spatial CPM model. Panel (a):
    systemic-species trajectories (TNF `T`, ROS `X`, antibody `A`, APC `P`)
    over MCS. Panel (b): spatial cell-count aggregates (I/M/K/E) over MCS,
    plus the Task-8.3 dynamic `sig_1` series (on a twin axis, different
    scale). A small, fast real run (side=30, steps=8, seed=0) -- non-vacuous:
    plotted data must match the driver's actual return series/keys."""
    result = run_global_coupling(side=30, steps=8, seed=0)

    fig = viz.global_coupling_figure(result)

    assert isinstance(fig, Figure)
    assert len(fig.axes) >= 2
    ax_ode, ax_spatial = fig.axes[0], fig.axes[1]

    # panel (a): systemic ODE species (T/X/A/P) vs MCS -- non-vacuous,
    # matching the driver's actual "ode" series.
    assert len(ax_ode.lines) >= 4
    ode_ydata_by_line = [list(line.get_data()[1]) for line in ax_ode.lines]
    for sp in ("T", "X", "A", "P"):
        assert result["ode"][sp] in ode_ydata_by_line

    # panel (b): spatial aggregates I/M/K/E -- non-vacuous, matching the
    # driver's actual "spatial" series.
    assert len(ax_spatial.lines) >= 4
    spatial_ydata_by_line = [list(line.get_data()[1]) for line in ax_spatial.lines]
    for sp in ("I", "M", "K", "E"):
        assert result["spatial"][sp] in spatial_ydata_by_line

    # sigma1 (Task 8.3's dynamic sig_1 series) is rendered somewhere in the
    # figure too (e.g. a twin axis on the spatial panel) -- non-vacuous.
    all_ydata = [list(line.get_data()[1]) for ax in fig.axes for line in ax.lines]
    assert result["sigma1"] in all_ydata


def test_repro_fig3b_figure_returns_figure_with_axes():
    """Task 9.5: CAPSTONE viz -- `run.repro_fig3b(...)`'s ensemble overlaid
    on the digitized Fig-3B acceptance band. A tiny, fast real ensemble
    (replicas=1, cells_per_side=12, steps=8) -- non-vacuous: at least 2 axes
    (cell-count panel + virus panel), each with a plotted model line AND a
    plotted/filled band from the target."""
    result = repro_fig3b(replicas=1, cells_per_side=12, steps=8)

    fig = viz.repro_fig3b_figure(result)

    assert isinstance(fig, Figure)
    assert len(fig.axes) >= 2
    for ax in fig.axes:
        assert len(ax.lines) >= 1          # ensemble-mean model line(s)
        assert len(ax.collections) >= 1    # the fill_between acceptance band


def test_repro_sweep_figure_returns_figure_with_axes_for_fig5_and_fig7():
    """Task 9.5: CAPSTONE viz -- `run.repro_fig5`/`run.repro_fig7`'s
    per-dose sweep overlaid on the highest-dose Fig-5/Fig-7 acceptance band.
    Tiny, fast real sweeps (2 doses, replicas=1, cells_per_side=12, steps=8)
    -- non-vacuous: 2 axes (dose-response curve + band overlay), each with
    at least one plotted line."""
    fig5_result = repro_fig5(loads=(1, 10), replicas=1, cells_per_side=12, steps=8)
    fig5 = viz.repro_sweep_figure(fig5_result, "fig5")
    assert isinstance(fig5, Figure)
    assert len(fig5.axes) >= 2
    for ax in fig5.axes:
        assert len(ax.lines) >= 1

    fig7_result = repro_fig7(fracs=(0.001, 0.05), replicas=1, cells_per_side=12, steps=8)
    fig7 = viz.repro_sweep_figure(fig7_result, "fig7")
    assert isinstance(fig7, Figure)
    assert len(fig7.axes) >= 2
    for ax in fig7.axes:
        assert len(ax.lines) >= 1
