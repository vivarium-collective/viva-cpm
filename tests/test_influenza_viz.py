import matplotlib
matplotlib.use("Agg")  # headless: no display backend required

from matplotlib.figure import Figure

from cpm.schema import load_world
from pbg_cpm_studies.influenza import sheet, viz
from pbg_cpm_studies.influenza.run import run_virus_infection


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
