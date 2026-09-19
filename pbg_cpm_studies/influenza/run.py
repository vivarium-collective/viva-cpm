"""Task 2.3: drive the virus field (2.1) + infection transition (2.2) over the
Increment-1 confluent epithelial sheet and report per-update population +
spatial-spread metrics.

Loop per update: `world.step(mcs_per_update)` advances Potts + the virus field
together (see `fields.py`'s module docstring -- callers must NOT call
`advance_fields` separately), then `infection_step` decides H -> I transitions
from the freshly-diffused field, and any changed cells are written back with
`world.set_cell_type`.

Behavior (lesion growth over time, staying local early) is the fidelity
criterion here, not exact counts -- this is mechanism validation, not a
Fig-3B/5/7 reproduction (Increment 9). Two independent, seed-derived
`numpy.random.Generator`s make the run fully deterministic: one for choosing
which cells start infected, one for the per-update stochastic infection
transition.
"""
from __future__ import annotations

import math

import numpy as np

from . import build, fields, sheet, transitions, types
from .params import load_params
from .resistance import cell_resistance


def run_virus_infection(patch_mm: float = 0.3, steps: int = 60, seed: int = 17,
                         init_infected_frac: float = 0.05,
                         mcs_per_update: int = 10) -> dict:
    """Run the CPM + virus field + infection transition for `steps` updates.

    Returns a dict of equal-length per-update series (index 0 = the seeded
    initial state, before any updates):
      - "steps": update index (0..steps)
      - "n_H", "n_I": uninfected / infected cell counts
      - "total_virus": sum of the virus field concentration over the lattice
      - "infected_ids": sorted list of currently-infected cell ids
      - "new_infection_dists": for each update, the list of (one per newly
        infected cell) Euclidean distances (lattice sites) from that cell's
        centroid to the nearest *already*-infected cell's centroid -- the
        locality signal (empty list if no new infections that update).
    """
    params = load_params()
    g_hv = float(params["virus"]["infection_g_hv"])

    spec = sheet.build_sheet_spec(patch_mm, seed=seed)
    n_cells = len(spec["cells"])  # ids 1..n_cells, all seeded as H

    world = build.world_from_spec(spec, finalize=False)
    field_idx = fields.add_virus_field(world)

    # Seed the initial infection BEFORE finalize (World.set_cell_type works
    # pre-finalize; doing it here, rather than post-finalize, matches the
    # brief and means the seeded I cells secrete from MCS 0).
    init_rng = np.random.default_rng(seed)
    n_init = int(round(n_cells * init_infected_frac))
    initial_infected_ids = (
        sorted(int(c) for c in init_rng.choice(np.arange(1, n_cells + 1),
                                                size=n_init, replace=False))
        if n_init > 0 else []
    )
    for cid in initial_infected_ids:
        world.set_cell_type(cid, types.I)

    world.finalize(int(spec["potts"]["seed"]))

    # Separate stream from init_rng so the choice of initial lesion and the
    # per-update stochastic transition don't share (and don't accidentally
    # correlate) random state.
    infect_rng = np.random.default_rng(seed + 1)

    coms = world.cell_coms()  # fixed centroids, index = cell id (0 = medium)
    current_types = list(world.cell_types())
    infected_ids = set(initial_infected_ids)

    result = {
        "steps": [], "n_H": [], "n_I": [], "total_virus": [],
        "infected_ids": [], "new_infection_dists": [],
    }

    def _record(step_idx, new_dists):
        result["steps"].append(step_idx)
        result["n_H"].append(sum(1 for t in current_types[1:] if t == types.H))
        result["n_I"].append(sum(1 for t in current_types[1:] if t == types.I))
        result["total_virus"].append(float(sum(world.field_conc(field_idx))))
        result["infected_ids"].append(sorted(infected_ids))
        result["new_infection_dists"].append(new_dists)

    _record(0, [])

    for step_idx in range(1, steps + 1):
        world.step(mcs_per_update)

        virus_at_cell = [world.field_mean_at_cell(field_idx, cid)
                          for cid in range(n_cells + 1)]
        new_types = transitions.infection_step(current_types, virus_at_cell, g_hv, infect_rng)

        new_dists = []
        for cid in range(1, n_cells + 1):
            if new_types[cid] == types.I and current_types[cid] != types.I:
                world.set_cell_type(cid, types.I)
                if infected_ids:
                    cx, cy, _ = coms[cid]
                    d = min(math.hypot(cx - coms[other][0], cy - coms[other][1])
                            for other in infected_ids)
                    new_dists.append(d)
                infected_ids.add(cid)

        current_types = new_types
        _record(step_idx, new_dists)

    return result


def run_virus_infection_with_ifn(patch_mm: float = 0.3, steps: int = 60, seed: int = 17,
                                  init_infected_frac: float = 0.05,
                                  mcs_per_update: int = 10) -> dict:
    """Increment-3 mechanism crux: same driver as `run_virus_infection`, but
    with the type-I IFN field (3.2) wired alongside the virus field, and a
    per-cell resistance computed from the local IFN each update that gates
    the infected cell's VIRUS secretion by `(1 - resist)` via the Task-3.1
    per-cell secretion-scale primitive.

    Each update, AFTER `world.step` advances both fields together: for every
    currently-infected cell, `f_bar = world.field_mean_at_cell(ifn_fi, cid)`,
    `resist = cell_resistance(f_bar, a_rf)`, then
    `world.set_cell_secretion_scale(virus_fi, cid, 1.0 - resist)` -- so the
    NEXT update's virus secretion (and hence the infection transition it
    feeds) is throttled by accumulated local IFN. The infection transition
    itself is unchanged from `run_virus_infection`.

    Returns the same per-update series as `run_virus_infection`, plus
    "mean_resist": the mean `resist` over currently-infected cells that
    update (0.0 when there are no infected cells yet, e.g. index 0 before
    any IFN has accumulated).
    """
    params = load_params()
    g_hv = float(params["virus"]["infection_g_hv"])
    a_rf = float(params["resistance"]["a_rf"])

    spec = sheet.build_sheet_spec(patch_mm, seed=seed)
    n_cells = len(spec["cells"])  # ids 1..n_cells, all seeded as H

    world = build.world_from_spec(spec, finalize=False)
    virus_fi = fields.add_virus_field(world)
    ifn_fi = fields.add_ifn_field(world)

    init_rng = np.random.default_rng(seed)
    n_init = int(round(n_cells * init_infected_frac))
    initial_infected_ids = (
        sorted(int(c) for c in init_rng.choice(np.arange(1, n_cells + 1),
                                                size=n_init, replace=False))
        if n_init > 0 else []
    )
    for cid in initial_infected_ids:
        world.set_cell_type(cid, types.I)

    world.finalize(int(spec["potts"]["seed"]))

    infect_rng = np.random.default_rng(seed + 1)

    coms = world.cell_coms()
    current_types = list(world.cell_types())
    infected_ids = set(initial_infected_ids)

    result = {
        "steps": [], "n_H": [], "n_I": [], "total_virus": [],
        "infected_ids": [], "new_infection_dists": [], "mean_resist": [],
    }

    def _record(step_idx, new_dists, mean_resist):
        result["steps"].append(step_idx)
        result["n_H"].append(sum(1 for t in current_types[1:] if t == types.H))
        result["n_I"].append(sum(1 for t in current_types[1:] if t == types.I))
        result["total_virus"].append(float(sum(world.field_conc(virus_fi))))
        result["infected_ids"].append(sorted(infected_ids))
        result["new_infection_dists"].append(new_dists)
        result["mean_resist"].append(mean_resist)

    _record(0, [], 0.0)

    for step_idx in range(1, steps + 1):
        world.step(mcs_per_update)

        # Per-cell resistance from the freshly-diffused local IFN, gating
        # THIS cell's virus secretion for subsequent updates. Recorded
        # BEFORE the infection transition below (which reads the virus
        # field the way it was left by the PRIOR update's secretion scale,
        # matching the source's "resist computed this MCS throttles this
        # MCS's release" ordering).
        resists = []
        for cid in sorted(infected_ids):
            f_bar = world.field_mean_at_cell(ifn_fi, cid)
            resist = cell_resistance(f_bar, a_rf)
            world.set_cell_secretion_scale(virus_fi, cid, 1.0 - resist)
            resists.append(resist)
        mean_resist = sum(resists) / len(resists) if resists else 0.0

        virus_at_cell = [world.field_mean_at_cell(virus_fi, cid)
                          for cid in range(n_cells + 1)]
        new_types = transitions.infection_step(current_types, virus_at_cell, g_hv, infect_rng)

        new_dists = []
        for cid in range(1, n_cells + 1):
            if new_types[cid] == types.I and current_types[cid] != types.I:
                world.set_cell_type(cid, types.I)
                if infected_ids:
                    cx, cy, _ = coms[cid]
                    d = min(math.hypot(cx - coms[other][0], cy - coms[other][1])
                            for other in infected_ids)
                    new_dists.append(d)
                infected_ids.add(cid)

        current_types = new_types
        _record(step_idx, new_dists, mean_resist)

    return result
