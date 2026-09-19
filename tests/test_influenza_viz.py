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
