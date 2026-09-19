import matplotlib
matplotlib.use("Agg")  # headless: no display backend required

from matplotlib.figure import Figure

from cpm.schema import load_world
from pbg_cpm_studies.influenza import sheet, viz
from pbg_cpm_studies.influenza.run import run_virus_infection, run_virus_infection_with_ifn


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
