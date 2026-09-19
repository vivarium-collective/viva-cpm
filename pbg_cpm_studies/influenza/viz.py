"""Minimal in-package visualizations for the influenza-sego2022 capability
ladder (Increments 1-3): the epithelial-sheet baseline, the virus-field
infection mechanism, and the IFN/resistance mechanism.

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
