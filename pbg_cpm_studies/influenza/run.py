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
