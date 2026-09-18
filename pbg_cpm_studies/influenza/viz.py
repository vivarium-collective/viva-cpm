"""Minimal visualization for the epithelial-sheet baseline (Increment 1).

The dashboard workbench isn't the delivery path for this increment (no
server running in this worktree), so this is a plain in-package matplotlib
function: a cell-type lattice snapshot plus the cell-volume distribution
histogram, per the spec's minimal figure requirement (spec §5).
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
