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

from . import (allee, build, fields, immune, killing, price_ode, recruitment,
               sheet, signaling, transitions, types)
from .params import load_params
from .resistance import cell_resistance

# Epithelial cell types that participate in the Allee contact-surface
# calculation (matches the source's `ec_list = [UNINFECTED, INFECTED,
# INFECTEDRELEASING, DYING]`, i.e. H/I/D here) -- MEDIUM and any future
# immune-cell contact must be excluded from `srf_total`/`srf_uninfected`,
# not just tallied by `world.cell_contact_area_by_type`'s raw dict.
_EPITHELIAL_TYPES = (types.H, types.I, types.D)

# Task 8.4 recruitment reserve pool: the engine can only create cells
# PRE-finalize (`World.add_cell`/`seed_block` error after finalize -- see
# `crates/cpm-py/src/lib.rs`), so ODE-driven inflow cannot mint brand-new CPM
# cells mid-run the way the source's `new_immune_cell_by_type` does. Instead we
# pre-seed a pool of DORMANT reserve cells (this distinct type, which secretes
# nothing and is never counted as M/K/E) parked in the domain's Medium margins;
# an inflow event ACTIVATES one via `set_cell_type(id, target_type)` (it then
# counts + secretes like any recruited cell), and an outflow event `remove_cells`
# an active cell. Reserve cells are parked in an interior reservoir rather than
# placed at the lesion (the source seeds onto Medium near the target field's
# peak) -- an engine-imposed geometry approximation, documented in
# task-8.4-report.md; this task's observable is the POPULATION COUNT (no killing
# consumes placement until Task 8.5), for which reservoir vs lesion placement is
# immaterial.
RECRUIT_RESERVE_TYPE = 7


def _seed_recruit_pool(world, spec, *, pool_per_type, targets, target_volume,
                       lambda_volume):
    """Pre-seed (BEFORE `world.finalize`) `pool_per_type` dormant reserve cells
    (`RECRUIT_RESERVE_TYPE`) per target immune type, packed into the Medium
    strips above and below the scenario's occupied bounding box. Returns
    ``{target_type: [reserve_cell_id, ...]}`` -- FIFO queues an inflow event pops
    from to activate a cell of that type. Raises if the strips can't hold the
    requested pool (caller keeps `pool_per_type` small enough to fit)."""
    nx, ny, _nz = spec["potts"]["dims"]
    blocks = [c["seed_block"] for c in spec["cells"]]
    content_y0 = min(b[1] for b in blocks)
    content_y1 = max(b[4] for b in blocks)

    side = int(round(target_volume ** 0.5))
    pitch = side + 2  # 2-site Medium gap between reserves (matches cluster gap)

    # Candidate (x0, y0) block origins in the top strip [1, content_y0) and the
    # bottom strip (content_y1, ny), left-to-right then top-to-bottom.
    def _strip_origins(y_lo, y_hi):
        origins = []
        y = y_lo
        while y + side <= y_hi:
            x = 1
            while x + side <= nx - 1:
                origins.append((x, y))
                x += pitch
            y += pitch
        return origins

    origins = _strip_origins(1, content_y0 - 1) + _strip_origins(content_y1 + 1, ny - 1)
    need = pool_per_type * len(targets)
    if len(origins) < need:
        raise ValueError(
            f"recruit reserve pool needs {need} slots but only {len(origins)} fit "
            f"in the {nx}x{ny} domain's Medium margins; reduce recruit_pool_per_type")

    pool: dict[int, list[int]] = {t: [] for t in targets}
    slot = 0
    for t in targets:
        for _ in range(pool_per_type):
            x0, y0 = origins[slot]
            slot += 1
            cid = world.add_cell(RECRUIT_RESERVE_TYPE, float(target_volume),
                                 float(lambda_volume), 0.0, 0.0)
            world.seed_block(cid, x0, y0, 0, x0 + side, y0 + side, 1)
            pool[t].append(cid)
    return pool


def _epithelial_contact_totals(world, cid):
    """srf_uninfected/srf_D/srf_total for cell `cid`, restricted to
    epithelial neighbor types (see `_EPITHELIAL_TYPES`). `world.
    cell_contact_area_by_type(cid)` (Task 4.1) returns contact area for
    EVERY neighbor type touching `cid`, including MEDIUM (key 0) at sheet
    edges -- the source's `RecoverySteppable.step` only accumulates
    `srf_area_total`/`srf_uninfected` over neighbors whose type is in
    `ec_list`, so MEDIUM contact must NOT inflate `srf_total` here (an
    edge-of-sheet H cell touching medium is not thereby "poorly
    surrounded" in the source's sense)."""
    contact = world.cell_contact_area_by_type(cid)
    srf_uninfected = contact.get(types.H, 0)
    srf_D = contact.get(types.D, 0)
    srf_total = sum(v for t, v in contact.items() if t in _EPITHELIAL_TYPES)
    return srf_uninfected, srf_D, srf_total


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


def run_epithelial_fate(patch_mm: float = 0.3, steps: int = 200, seed: int = 17,
                         init_infected_frac: float = 0.05,
                         mcs_per_update: int = 10) -> dict:
    """Increment-4 crux driver: wire ALL epithelial fates together --
    infection (H -> I, Task 2.2), infected death (I -> D, Task 4.2), and
    the cellularized Allee effect (H -> D death / D -> H recovery, Task
    4.3, this task) -- on top of the virus + type-I IFN fields (2.1/3.2).

    Each update, AFTER `world.step` advances both fields together:
      1. per-cell resist from local IFN, for EVERY epithelial cell (H, I,
         D all carry local IFN in the source, not just I -- see
         `resistance.cell_resistance`); gate I-cell virus secretion by
         `(1 - resist)` for the NEXT update (same convention as
         `run_virus_infection_with_ifn`).
      2. infection H -> I (`transitions.infection_step`, virus field).
      3. infected death I -> D (`transitions.infected_death_step`, resist
         + `cell_death.mu_i_per_mcs`).
      4. Allee death H -> D / recovery D -> H (`allee.allee_death_rate` /
         `allee.allee_recovery_rate`), evaluated per epithelial cell from
         `world.cell_contact_area_by_type` (`_epithelial_contact_totals`).
    Steps 2-3's type changes are written to the world via `world.
    set_cell_type` BEFORE step 4 reads neighbor contact -- matching the
    source's within-MCS sequencing (`RecoverySteppable` runs after the
    infection/death steppables in the same MCS, and CC3D applies
    `cell.type` immediately, so `RecoverySteppable` sees THIS MCS's fresh
    infection/death types when tallying neighbor contact, not last
    update's stale types). Step 4's OWN death/recovery writes are still
    batched and applied together at the very end of the update (a much
    smaller residual approximation than the one just fixed -- see the
    in-loop comment and task-4.3-report.md). This otherwise follows this
    module's existing per-update-not-per-MCS driver convention (see
    `run_virus_infection`).

    Returns per-update series: "steps", "n_H"/"n_I"/"n_D", "total_virus",
    plus "n_allee_death"/"n_allee_recovery" (count of Allee-driven
    transitions THAT update -- 0 for index 0, the seeded initial state; a
    diagnostic for checking the Allee mechanism actually engaged, not part
    of the brief's required series).
    """
    params = load_params()
    g_hv = float(params["virus"]["infection_g_hv"])
    a_rf = float(params["resistance"]["a_rf"])
    mu_i = float(params["cell_death"]["mu_i_per_mcs"])
    b_h = float(params["allee"]["b_h"])
    theta = float(params["allee"]["srf_threshold"])

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

    # Independent RNG streams so the initial-lesion draw, the infection
    # transition, the infected-death transition, and the Allee draws don't
    # share (or accidentally correlate) random state.
    infect_rng = np.random.default_rng(seed + 1)
    death_rng = np.random.default_rng(seed + 2)
    allee_rng = np.random.default_rng(seed + 3)

    current_types = list(world.cell_types())

    result = {
        "steps": [], "n_H": [], "n_I": [], "n_D": [], "total_virus": [],
        "n_allee_death": [], "n_allee_recovery": [],
    }

    def _record(step_idx, n_allee_death, n_allee_recovery):
        result["steps"].append(step_idx)
        result["n_H"].append(sum(1 for t in current_types[1:] if t == types.H))
        result["n_I"].append(sum(1 for t in current_types[1:] if t == types.I))
        result["n_D"].append(sum(1 for t in current_types[1:] if t == types.D))
        result["total_virus"].append(float(sum(world.field_conc(virus_fi))))
        result["n_allee_death"].append(n_allee_death)
        result["n_allee_recovery"].append(n_allee_recovery)

    _record(0, 0, 0)

    for step_idx in range(1, steps + 1):
        world.step(mcs_per_update)

        # (1) per-cell resist from local IFN -- for every epithelial cell
        # (H, I, D), not just infected ones (H and D cells have local IFN
        # too, per the brief).
        resist_at_cell = [0.0] * (n_cells + 1)
        for cid in range(1, n_cells + 1):
            if current_types[cid] in _EPITHELIAL_TYPES:
                f_bar = world.field_mean_at_cell(ifn_fi, cid)
                resist_at_cell[cid] = cell_resistance(f_bar, a_rf)
        for cid in range(1, n_cells + 1):
            if current_types[cid] == types.I:
                world.set_cell_secretion_scale(virus_fi, cid, 1.0 - resist_at_cell[cid])

        # (2) infection H -> I.
        virus_at_cell = [world.field_mean_at_cell(virus_fi, cid)
                          for cid in range(n_cells + 1)]
        types_after_infection = transitions.infection_step(
            current_types, virus_at_cell, g_hv, infect_rng)

        # (3) infected death I -> D.
        types_after_death = transitions.infected_death_step(
            types_after_infection, resist_at_cell, mu_i, death_rng)

        # Apply the infection + infected-death writes to the world NOW,
        # BEFORE the Allee step reads neighbor contact -- matching the
        # source's within-MCS sequencing (RecoverySteppable runs AFTER the
        # infection/death steppables in the same MCS, and CC3D applies
        # cell.type changes immediately, so RecoverySteppable's
        # `get_cell_neighbor_data_list` sees THIS MCS's fresh infection/
        # death types, not last update's). Committing here (rather than
        # batching these into the final write-back below) is the fix for
        # a prior review finding: reading `cell_contact_area_by_type`
        # against stale (pre-update) neighbor types systematically biased
        # Allee death DOWN and recovery UP.
        for cid in range(1, n_cells + 1):
            if types_after_death[cid] != current_types[cid]:
                world.set_cell_type(cid, types_after_death[cid])

        # (4) Allee death H -> D / recovery D -> H, from contact geometry
        # as of the world's state AFTER this update's infection/death
        # writes (just applied above) but BEFORE any Allee writes -- i.e.
        # every cell's Allee rate this update is computed from the SAME
        # neighbor snapshot (this update's post-infection/death types).
        # Residual approximation (deliberately not fixed further, see
        # task-4.3-report.md): the source applies each cell's death/revive
        # write immediately within RecoverySteppable's own per-cell loop,
        # so a cell processed later in the source's iteration order can see
        # an EARLIER cell's Allee write this same MCS; here all Allee
        # writes are collected and applied together at the end of this
        # step, so within-Allee-pass neighbor propagation is one iteration
        # behind the source. Much smaller effect than the infection/death
        # staleness this fix addresses (Allee write density per update is
        # far lower than the infection/death write density it now reads
        # freshly).
        final_types = list(types_after_death)
        n_allee_death = 0
        n_allee_recovery = 0
        for cid in range(1, n_cells + 1):
            t = types_after_death[cid]
            if t not in (types.H, types.D):
                continue
            srf_uninfected, srf_D, srf_total = _epithelial_contact_totals(world, cid)
            resist = resist_at_cell[cid]
            if t == types.H:
                rate = allee.allee_death_rate(srf_D, srf_uninfected, srf_total,
                                               b_h, resist, theta)
                pr = 1.0 - math.exp(-rate)
                if allee_rng.random() < pr:
                    final_types[cid] = types.D
                    n_allee_death += 1
            else:  # types.D
                rate = allee.allee_recovery_rate(srf_uninfected, srf_total,
                                                  b_h, resist, theta)
                pr = 1.0 - math.exp(-rate)
                if allee_rng.random() < pr:
                    final_types[cid] = types.H
                    n_allee_recovery += 1

        for cid in range(1, n_cells + 1):
            if final_types[cid] != types_after_death[cid]:
                world.set_cell_type(cid, final_types[cid])

        current_types = final_types
        _record(step_idx, n_allee_death, n_allee_recovery)

    return result


def run_macrophage_response(*, epithelial_cells_per_side: int = 4, n_infected: int = 1,
                             n_macrophages: int = 6, margin_sites: int = 30,
                             separation_sites: int = 25,
                             steps: int = 40, seed: int = 17, mcs_per_update: int = 10,
                             field_warmup: int = 3000,
                             chemotaxis_lambda: float | None = None) -> dict:
    """Task 5.1 crux driver: build the non-confluent macrophage scenario
    (`immune.build_macrophage_scenario_spec` -- both the epithelial/infection
    patch and the macrophage cluster sit in the domain's INTERIOR, each
    ``margin_sites`` from every noflux wall; see that function's docstring
    for why interior placement matters -- a wall-pinned control biases the
    lambda=0 baseline's drift), wire the virus field (`fields.
    add_virus_field` -- I cells secrete) + macrophage chemotaxis (`immune.
    set_macrophage_chemotaxis`), warm the field up (`World.advance_fields`,
    matching `pbg_cpm_studies.chemotaxis.run`'s WARMUP convention -- builds
    an initial gradient reaching the macrophages' starting cluster before
    they start responding to it), then step the coupled CPM + field world and
    record the macrophages' localization each update.

    ``chemotaxis_lambda`` overrides ``params.yaml``'s ``macrophage.
    chemotaxis_v_macro`` (5000); pass ``0.0`` for the lambda=0 control (no
    directed chemotaxis -- see `immune.set_macrophage_chemotaxis`).

    Returns per-update series (index 0 = the seeded initial state, before any
    Potts/field updates -- but AFTER the field warmup):
      - "steps": update index (0..steps)
      - "mean_distance_to_infection": macrophages' mean centre-of-mass
        distance to the infected-cell centroid
      - "macrophage_com": macrophage centre-of-mass (x, y), for the
        companion localization-trajectory readout
    plus "params" (the scenario/run knobs, for the report).
    """
    params = load_params()
    lam = (float(params["macrophage"]["chemotaxis_v_macro"])
           if chemotaxis_lambda is None else float(chemotaxis_lambda))

    spec = immune.build_macrophage_scenario_spec(
        epithelial_cells_per_side=epithelial_cells_per_side, n_infected=n_infected,
        n_macrophages=n_macrophages, margin_sites=margin_sites,
        separation_sites=separation_sites, seed=seed)

    cell_type_by_idx = [c["type"] for c in spec["cells"]]  # spec index i -> cell id i+1
    infected_ids = [i + 1 for i, t in enumerate(cell_type_by_idx) if t == types.I]
    macrophage_ids = [i + 1 for i, t in enumerate(cell_type_by_idx) if t == types.M]

    world = build.world_from_spec(spec, finalize=False)
    virus_fi = fields.add_virus_field(world)
    world.finalize(int(spec["potts"]["seed"]))

    for _ in range(field_warmup):
        world.advance_fields(1)

    immune.set_macrophage_chemotaxis(world, virus_fi, chemotaxis_v_macro=lam)

    def _infection_centroid(coms):
        xs = [coms[cid][0] for cid in infected_ids]
        ys = [coms[cid][1] for cid in infected_ids]
        return (sum(xs) / len(xs), sum(ys) / len(ys))

    def _macrophage_mean_distance(coms, centroid):
        ds = [math.hypot(coms[cid][0] - centroid[0], coms[cid][1] - centroid[1])
              for cid in macrophage_ids]
        return sum(ds) / len(ds)

    def _macrophage_com(coms):
        xs = [coms[cid][0] for cid in macrophage_ids]
        ys = [coms[cid][1] for cid in macrophage_ids]
        return (sum(xs) / len(xs), sum(ys) / len(ys))

    result = {"steps": [], "mean_distance_to_infection": [], "macrophage_com": []}

    def _record(step_idx):
        coms = world.cell_coms()
        centroid = _infection_centroid(coms)
        result["steps"].append(step_idx)
        result["mean_distance_to_infection"].append(
            round(_macrophage_mean_distance(coms, centroid), 3))
        result["macrophage_com"].append(
            tuple(round(v, 2) for v in _macrophage_com(coms)))

    _record(0)
    for step_idx in range(1, steps + 1):
        world.step(mcs_per_update)
        _record(step_idx)

    result["params"] = {
        "chemotaxis_v_macro": lam, "n_macrophages": n_macrophages,
        "n_infected": n_infected, "epithelial_cells_per_side": epithelial_cells_per_side,
        "margin_sites": margin_sites, "separation_sites": separation_sites,
        "seed": seed, "steps": steps,
        "mcs_per_update": mcs_per_update, "field_warmup": field_warmup,
    }
    return result


def _radial_field_profile(world, field_idx, center_xy, n_bins=6):
    """Mean field concentration over `n_bins` concentric distance bins from
    `center_xy` (x, y in lattice-site units), computed over the RAW lattice
    (`world.field_conc`/`world.dims`) -- independent of where any CELL
    happens to sit, unlike `field_mean_at_cell`'s per-cell averages. The
    direct "gradient centered on a location, decaying with distance"
    evidence for `run_macrophage_signaling`'s chemokine assertion.

    `world`'s flat field array is laid out `idx = x + y*nx` (z is always 0
    in this z=1 domain) -- `crates/cpm-core/src/lattice.rs`'s
    `index`/`coords` -- so `conc.reshape(ny, nx)` (standard C-order,
    row-major) recovers `conc[y, x]`.

    Returns `(bin_edges, bin_means)`: `bin_edges` has `n_bins + 1` entries,
    0 .. the domain's max corner-to-corner distance from `center_xy`;
    `bin_means[i]` is the mean concentration over lattice sites whose
    distance to `center_xy` falls in `[bin_edges[i], bin_edges[i+1]]` (the
    last bin's upper edge is inclusive so the single farthest corner site
    isn't dropped); `nan` if no site falls in a bin.
    """
    nx, ny, _nz = world.dims()
    conc = np.asarray(world.field_conc(field_idx), dtype=float).reshape(ny, nx)
    xs, ys = np.meshgrid(np.arange(nx), np.arange(ny))
    dist = np.hypot(xs - center_xy[0], ys - center_xy[1])
    max_d = float(dist.max())
    bin_edges = np.linspace(0.0, max_d, n_bins + 1)
    bin_means = []
    for i in range(n_bins):
        lo, hi = bin_edges[i], bin_edges[i + 1]
        mask = (dist >= lo) & (dist <= hi if i == n_bins - 1 else dist < hi)
        bin_means.append(float(conc[mask].mean()) if mask.any() else float("nan"))
    return bin_edges.tolist(), bin_means


def run_macrophage_signaling(*, epithelial_cells_per_side: int = 4, n_infected: int = 1,
                              n_macrophages: int = 6, margin_sites: int = 30,
                              separation_sites: int = 25,
                              steps: int = 20, seed: int = 17, mcs_per_update: int = 10,
                              field_warmup: int = 0,
                              chemotaxis_lambda: float | None = None,
                              radial_bins: int = 6) -> dict:
    """Task 6.1 crux driver: extend `run_macrophage_response`'s Increment-5
    non-confluent macrophage scenario (`immune.build_macrophage_scenario_spec`)
    with the macrophage-released chemokine + IL-10 diffusive fields (`fields.
    add_chemokine_field`/`add_il10_field`), on top of the existing virus
    (`fields.add_virus_field`) + type-I IFN (`fields.add_ifn_field`) fields,
    and the per-cell secretion regulation this task adds
    (`signaling.macrophage_secretion_scale`/`uninfected_il10_scale`).

    Each update, AFTER `world.step` advances all four fields together
    (`fields.py`'s module docstring -- callers must not call `advance_fields`
    separately once stepping has started):
      - every macrophage (`types.M`) cell: `L_loc = world.
        field_mean_at_cell(il10_fi, cid)` (local mean IL-10), `scale =
        signaling.macrophage_secretion_scale(L_loc, *fields.
        il10_hill_constants())`, applied via `world.
        set_cell_secretion_scale` to BOTH the chemokine and the IL-10
        field's per-cell scale for that cell -- source: the SAME
        IL-10-Hill self-regulation gates both fields (`b_c`/`b_l` share
        `sig_1`/`g_1`/`g_2`/`d_2`, see `signaling.py`'s module docstring /
        the sig_1-STUB ruling there).
      - every uninfected (`types.H`) epithelial cell: `f_bar = world.
        field_mean_at_cell(ifn_fi, cid)`, `resist = resistance.
        cell_resistance(f_bar, a_rf)`, `scale = signaling.
        uninfected_il10_scale(resist)`, applied to the IL-10 field's
        per-cell scale for that cell.
    The IFN field here exists ONLY to supply `resist` for the IL-10
    uninfected-cell gate -- infection/virus-release throttling (Increment
    3/4's mechanism) is out of THIS task's scope (chemokine/IL-10 field
    fidelity) and is not modeled by this driver.

    Returns per-update series (index 0 = the seeded initial state, before
    any Potts/field updates, but AFTER `field_warmup`):
      - "steps": update index (0..steps)
      - "total_chemokine", "total_il10": field sums over the lattice
      - "chemo_at_macrophages": mean chemokine field value AT macrophage
        cells (the near-the-source reading)
      - "chemo_at_uninfected": mean chemokine field value AT uninfected (H)
        epithelial cells, `separation_sites` of open Medium away from the
        macrophage cluster (the far-from-the-source reading)
      - "il10_at_macrophages", "il10_at_uninfected": same near/far split
        for the IL-10 field (both macrophage AND uninfected cells are IL-10
        SOURCES, unlike chemokine, so this pair is a secondary diagnostic,
        not the locality signal)
    plus, computed ONCE from the FINAL state:
      - "chemo_radial_profile": `{"bin_edges": [...], "bin_means": [...]}`
        from `_radial_field_profile` around the macrophage cluster's
        centroid -- the direct, cell-position-independent "gradient
        centered on the macrophage cluster, decaying with distance"
        evidence.
      - "params" (the scenario/run knobs, for the report).
    """
    params = load_params()
    a_rf = float(params["resistance"]["a_rf"])
    sig_1, g_1, g_2, d_2 = fields.il10_hill_constants()
    lam = (float(params["macrophage"]["chemotaxis_v_macro"])
           if chemotaxis_lambda is None else float(chemotaxis_lambda))

    spec = immune.build_macrophage_scenario_spec(
        epithelial_cells_per_side=epithelial_cells_per_side, n_infected=n_infected,
        n_macrophages=n_macrophages, margin_sites=margin_sites,
        separation_sites=separation_sites, seed=seed)

    cell_type_by_idx = [c["type"] for c in spec["cells"]]  # spec index i -> cell id i+1
    macrophage_ids = [i + 1 for i, t in enumerate(cell_type_by_idx) if t == types.M]
    uninfected_ids = [i + 1 for i, t in enumerate(cell_type_by_idx) if t == types.H]

    world = build.world_from_spec(spec, finalize=False)
    virus_fi = fields.add_virus_field(world)
    ifn_fi = fields.add_ifn_field(world)
    chemo_fi = fields.add_chemokine_field(world)
    il10_fi = fields.add_il10_field(world)
    world.finalize(int(spec["potts"]["seed"]))

    for _ in range(field_warmup):
        world.advance_fields(1)

    immune.set_macrophage_chemotaxis(world, virus_fi, chemotaxis_v_macro=lam)

    def _macrophage_centroid(coms):
        xs = [coms[cid][0] for cid in macrophage_ids]
        ys = [coms[cid][1] for cid in macrophage_ids]
        return (sum(xs) / len(xs), sum(ys) / len(ys))

    def _apply_secretion_scales():
        # Freshly-diffused (post-`world.step`) field readings gate the NEXT
        # update's secretion -- same "this update's local reading throttles
        # the following update's release" convention as
        # `run_virus_infection_with_ifn`/`run_epithelial_fate` above.
        for cid in macrophage_ids:
            l_loc = world.field_mean_at_cell(il10_fi, cid)
            scale = signaling.macrophage_secretion_scale(l_loc, sig_1, g_1, g_2, d_2)
            world.set_cell_secretion_scale(chemo_fi, cid, scale)
            world.set_cell_secretion_scale(il10_fi, cid, scale)
        for cid in uninfected_ids:
            f_bar = world.field_mean_at_cell(ifn_fi, cid)
            resist = cell_resistance(f_bar, a_rf)
            scale = signaling.uninfected_il10_scale(resist)
            world.set_cell_secretion_scale(il10_fi, cid, scale)

    result = {
        "steps": [], "total_chemokine": [], "total_il10": [],
        "chemo_at_macrophages": [], "chemo_at_uninfected": [],
        "il10_at_macrophages": [], "il10_at_uninfected": [],
    }

    def _record(step_idx):
        result["steps"].append(step_idx)
        result["total_chemokine"].append(float(sum(world.field_conc(chemo_fi))))
        result["total_il10"].append(float(sum(world.field_conc(il10_fi))))

        chemo_macro = [world.field_mean_at_cell(chemo_fi, cid) for cid in macrophage_ids]
        chemo_h = [world.field_mean_at_cell(chemo_fi, cid) for cid in uninfected_ids]
        il10_macro = [world.field_mean_at_cell(il10_fi, cid) for cid in macrophage_ids]
        il10_h = [world.field_mean_at_cell(il10_fi, cid) for cid in uninfected_ids]

        result["chemo_at_macrophages"].append(sum(chemo_macro) / len(chemo_macro))
        result["chemo_at_uninfected"].append(sum(chemo_h) / len(chemo_h))
        result["il10_at_macrophages"].append(sum(il10_macro) / len(il10_macro))
        result["il10_at_uninfected"].append(sum(il10_h) / len(il10_h))

    _record(0)
    for step_idx in range(1, steps + 1):
        world.step(mcs_per_update)
        _apply_secretion_scales()
        _record(step_idx)

    macro_centroid = _macrophage_centroid(world.cell_coms())
    bin_edges, bin_means = _radial_field_profile(world, chemo_fi, macro_centroid, radial_bins)
    result["chemo_radial_profile"] = {"bin_edges": bin_edges, "bin_means": bin_means}

    result["params"] = {
        "chemotaxis_v_macro": lam, "n_macrophages": n_macrophages,
        "n_infected": n_infected, "epithelial_cells_per_side": epithelial_cells_per_side,
        "margin_sites": margin_sites, "separation_sites": separation_sites,
        "seed": seed, "steps": steps,
        "mcs_per_update": mcs_per_update, "field_warmup": field_warmup,
    }
    return result


# Task 7.1 finding (documented, NOT a per-seed tune -- see `run_cytotoxic_
# response`'s docstring): the chemokine field's steady-state concentration is
# ~2-3 orders of magnitude smaller than the virus field's at comparable
# distances (chemokine's diffusion constant, ~15.6 lat^2/MCS, is ~87x the
# virus field's ~0.18 -- `fields.py`'s `_chemokine_field_params` -- so the
# same secreted mass spreads far more thinly). `World.set_chemotaxis`'s
# engine primitive computes ΔH = -lambda*(c(dest)-c(source)) directly from
# this engine's raw field units; measured empirically (single-secretor +
# single-responder minimal worlds, see `test_influenza_nk_cd8.py`'s engine-
# sanity tests and task-7.1-report.md), the source's LITERAL
# chemotaxis_v_nk=5000/chemotaxis_v_cd8=10000 produce a chemotaxis ΔH per
# copy attempt that is ~2-3 orders of magnitude below the thermal/adhesion
# noise floor (T=10, adhesion J~10-25) at this field's concentration scale --
# not merely "weaker", but statistically undetectable within any feasible
# step budget. This is the same category as the already-documented
# linear-vs-saturating functional-form gap (`immune.set_nk_cd8_chemotaxis`'s
# docstring): the source's own `chemotaxis_v_nk / (1 + concentration)`
# formula implies an internal concentration scale of order 1, not this
# field's order 1e-3, so the literal lambda value was never calibrated
# against this engine's concentration units in the first place.
#
# `run_cytotoxic_response` compensates with a single FIXED multiplier
# (applied identically to NK and CD8, preserving their 2:1 ratio, and
# identically across every seed/run -- not re-tuned per seed) so the
# localization BEHAVIOR (this increment's fidelity criterion, not the exact
# lambda magnitude) is demonstrable within a feasible step budget. Chosen
# from a small grid scan (50/200/500x at seed=17, `task-7.1-report.md`):
# 100x gives a clear, non-saturated on/off gap (200x/500x mostly saturate --
# the NK/CD8 cluster has already closed most of the feasible distance to the
# macrophage cluster, so the gap stops growing) and was confirmed robust
# (on < off in every one of 5 seeds, comfortably beyond the off-control's own
# seed-to-seed spread) before being adopted as this driver's default -- see
# `test_influenza_nk_cd8.py::test_nk_and_cd8_localize_to_infection_across_seeds`.
NK_CD8_CHEMOTAXIS_ENGINE_SCALE = 100.0


def run_cytotoxic_response(*, epithelial_cells_per_side: int = 4, n_infected: int = 1,
                            n_macrophages: int = 6, n_nk: int = 6, n_cd8: int = 6,
                            margin_sites: int = 20, separation_sites: int = 20,
                            steps: int = 60, seed: int = 17, mcs_per_update: int = 10,
                            field_warmup: int = 0,
                            macrophage_chemotaxis_lambda: float | None = None,
                            nk_chemotaxis_lambda: float | None = None,
                            cd8_chemotaxis_lambda: float | None = None,
                            enable_killing: bool = True) -> dict:
    """Task 7.1 crux driver (localization) + Task 7.2 (this task, contact
    killing): extend `run_macrophage_signaling`'s macrophage + chemokine/
    IL-10 scenario with a THIRD interior cluster of NK (`types.K`) + CD8+
    (`types.E`) cells (`immune.build_cytotoxic_scenario_spec` -- see that
    function's docstring for the "epithelial/infection | macrophage cluster |
    NK+CD8 cluster" interior layout), wires NK/CD8 chemotaxis UP THE
    CHEMOKINE FIELD (`immune.set_nk_cd8_chemotaxis` -- NOT the virus field;
    that stays the macrophages' own chemotaxis target, unchanged from
    `run_macrophage_response`/`run_macrophage_signaling`), and (Task 7.2,
    ``enable_killing``) applies `killing.contact_kill_rate`'s LOCAL
    (surface-contact) NK/CD8 contact-kill term to every currently-infected
    cell each update.

    Field wiring/regulation (virus, IFN, chemokine, IL-10 fields; per-cell
    IL-10-Hill macrophage secretion scale gating chemokine+IL-10; uninfected
    `(1-resist)` IL-10 gate; macrophage chemotaxis up virus) is IDENTICAL to
    `run_macrophage_signaling` -- this driver only adds the NK/CD8 cells +
    their own chemotaxis, and now (Task 7.2) contact-killing, on top.

    ``nk_chemotaxis_lambda``/``cd8_chemotaxis_lambda``, when given explicitly
    (including ``0.0`` for the lambda=0 control), are passed to `immune.
    set_nk_cd8_chemotaxis` UNSCALED. When left ``None`` (the default -- the
    "chemotaxis on" case), this driver passes `params.yaml`'s
    `nk.chemotaxis_v_nk` (5000)/`cd8.chemotaxis_v_cd8` (10000) each
    multiplied by `NK_CD8_CHEMOTAXIS_ENGINE_SCALE` (see that constant's
    module-level docstring for why: the literal source values are
    statistically undetectable at this engine's chemokine-field concentration
    scale, a documented, non-per-seed-tuned engine-unit compensation, not a
    change to `set_nk_cd8_chemotaxis`'s own default, which stays the literal
    unscaled source value per the brief).
    ``macrophage_chemotaxis_lambda`` overrides `macrophage.chemotaxis_v_macro`
    (5000, unrelated to the NK/CD8 on/off comparison -- left at its default
    in every run of this driver used by this task's tests).

    Scenario knob (separate from the lambda scale above, also fixed/not
    per-seed-tuned): ``separation_sites=20`` here, vs `run_macrophage_
    response`'s ``separation_sites=25`` for the virus-field macrophage case
    -- the NK/CD8 cluster starts closer to the macrophage cluster so the
    (still much weaker than virus) chemokine gradient reaching it has a
    shorter distance to work over. Still fully interior (>= margin_sites from
    every wall), same unbiased-control rigor as `immune.
    build_cytotoxic_scenario_spec`'s docstring.

    Task 7.2 killing step (each update, AFTER the chemotaxis/secretion-scale
    updates below, for every id still in the currently-infected set): local
    IFN self-secreted by the infected cell itself (`ifn_fi`, wired the same
    as every other driver in this module) gives `resist = resistance.
    cell_resistance(f_bar, a_rf)`; `world.cell_contact_area_by_type(cid)`
    (Increment 4's contact-geometry primitive, the same one `run_epithelial_
    fate`'s Allee step uses) gives `srf_nk`/`srf_cd8` (contact area with
    `types.K`/`types.E` neighbors); `killing.contact_kill_rate` gives
    `kill_rate_nk`/`kill_rate_cd8` from `nk.g_ik`/`cd8.g_ie` and
    `scaling.ode_epithelial_population` (`tot_ec_ODE`); `Pr = 1 -
    exp(-kill_rate)`. Draw order matches the source (`ContactKillingSteppable.
    step()`, see `killing.contact_kill_rate`'s module docstring): the NK draw
    is checked FIRST, and if it kills the cell (-> `types.D`, removed from the
    infected set), the CD8 draw is SKIPPED for that cell that update (matches
    the source's `continue`-past-CD8-if-NK-already-killed ordering). A
    dedicated RNG stream (`seed + 500`, unused by any other stream in this
    module) keeps the kill draws independent of the Potts/field RNG (set at
    `world.finalize`) and of every other driver's RNG offsets.

    **DISCREPANCY #7 (source-faithful, NOT "fixed" here -- see `killing.
    contact_kill_rate`'s module docstring and params.yaml's `nk.killing.
    kill_rate_resist_direct`/`cd8.killing.kill_rate_resist_direct`):** the
    kill rate multiplies by `cell_resist` DIRECTLY (not `1 - cell_resist`) --
    unlike every other resistance consumer in this module. A cell with MORE
    local (self-secreted) IFN is thus killed FASTER by contact NK/CD8 here,
    the opposite of every other resist-gated mechanism in this codebase.
    Implemented literally per the source; flagged again in task-7.2-report.md.

    ``enable_killing=False`` builds the IDENTICAL scenario/chemotaxis (NK/CD8
    cells present and still chemotaxing) but skips the killing step entirely
    -- the brief's "NK/CD8 present but killing disabled" control, so a
    same-seed on/off comparison isolates the killing step's own effect from
    any of localization's own run-to-run variation.

    Returns per-update series (index 0 = the seeded initial state, before
    any Potts/field updates but AFTER `field_warmup`):
      - "steps"
      - "nk_mean_distance_to_infection", "cd8_mean_distance_to_infection":
        NK/CD8 cells' mean centre-of-mass distance to the infected-cell
        centroid (frozen at its last value once no infected cells remain).
      - "total_chemokine": lattice-wide chemokine sum (sanity: confirms the
        gradient actually builds over the run).
      - "n_infected": count of currently `types.I` cells (Task 7.2; constant
        at ``n_infected`` throughout when ``enable_killing=False``).
    plus "params" (the scenario/run knobs, for the report).
    """
    params = load_params()
    a_rf = float(params["resistance"]["a_rf"])
    sig_1, g_1, g_2, d_2 = fields.il10_hill_constants()
    g_ik = float(params["nk"]["g_ik"])
    g_ie = float(params["cd8"]["g_ie"])
    tot_ec = float(params["scaling"]["ode_epithelial_population"])
    macro_lam = (float(params["macrophage"]["chemotaxis_v_macro"])
                 if macrophage_chemotaxis_lambda is None else float(macrophage_chemotaxis_lambda))
    # NK_CD8_CHEMOTAXIS_ENGINE_SCALE applies ONLY to the unspecified
    # ("chemotaxis on") default -- an explicit override (including the
    # lambda=0 control) is passed through unscaled, see this function's
    # docstring.
    nk_lam = (float(params["nk"]["chemotaxis_v_nk"]) * NK_CD8_CHEMOTAXIS_ENGINE_SCALE
              if nk_chemotaxis_lambda is None else float(nk_chemotaxis_lambda))
    cd8_lam = (float(params["cd8"]["chemotaxis_v_cd8"]) * NK_CD8_CHEMOTAXIS_ENGINE_SCALE
               if cd8_chemotaxis_lambda is None else float(cd8_chemotaxis_lambda))

    spec = immune.build_cytotoxic_scenario_spec(
        epithelial_cells_per_side=epithelial_cells_per_side, n_infected=n_infected,
        n_macrophages=n_macrophages, n_nk=n_nk, n_cd8=n_cd8,
        margin_sites=margin_sites, separation_sites=separation_sites, seed=seed)

    cell_type_by_idx = [c["type"] for c in spec["cells"]]  # spec index i -> cell id i+1
    # Mutable: Task 7.2's killing step removes ids as infected cells are
    # killed (-> types.D). Every OTHER id list above is fixed for the run
    # (macrophage/uninfected/NK/CD8 cell counts don't change in this driver).
    infected_ids = set(i + 1 for i, t in enumerate(cell_type_by_idx) if t == types.I)
    macrophage_ids = [i + 1 for i, t in enumerate(cell_type_by_idx) if t == types.M]
    uninfected_ids = [i + 1 for i, t in enumerate(cell_type_by_idx) if t == types.H]
    nk_ids = [i + 1 for i, t in enumerate(cell_type_by_idx) if t == types.K]
    cd8_ids = [i + 1 for i, t in enumerate(cell_type_by_idx) if t == types.E]

    world = build.world_from_spec(spec, finalize=False)
    virus_fi = fields.add_virus_field(world)
    ifn_fi = fields.add_ifn_field(world)
    chemo_fi = fields.add_chemokine_field(world)
    il10_fi = fields.add_il10_field(world)
    world.finalize(int(spec["potts"]["seed"]))

    for _ in range(field_warmup):
        world.advance_fields(1)

    immune.set_macrophage_chemotaxis(world, virus_fi, chemotaxis_v_macro=macro_lam)
    immune.set_nk_cd8_chemotaxis(world, chemo_fi, chemotaxis_v_nk=nk_lam, chemotaxis_v_cd8=cd8_lam)

    # Task 7.2 kill-draw RNG: a dedicated stream/offset (unused by any other
    # RNG in this module -- this driver otherwise has none, since it has no
    # H->I/I->D/Allee stochastic transitions, only movement, which draws from
    # the Potts RNG set at `world.finalize` above) so kill draws don't share
    # (or accidentally correlate with) any other source of randomness.
    kill_rng = np.random.default_rng(seed + 500)

    # Cache: once every infected cell has been killed, `infected_ids` is
    # empty and the NK/CD8 distance-to-infection series has no live centroid
    # to measure against -- freeze it at its last value rather than raising.
    _last_centroid = [None]

    def _infection_centroid(coms):
        if infected_ids:
            xs = [coms[cid][0] for cid in infected_ids]
            ys = [coms[cid][1] for cid in infected_ids]
            _last_centroid[0] = (sum(xs) / len(xs), sum(ys) / len(ys))
        return _last_centroid[0]

    def _mean_distance(coms, ids, centroid):
        ds = [math.hypot(coms[cid][0] - centroid[0], coms[cid][1] - centroid[1]) for cid in ids]
        return sum(ds) / len(ds)

    def _apply_secretion_scales():
        # Same per-update convention as `run_macrophage_signaling`: freshly-
        # diffused (post-`world.step`) field readings gate the NEXT update's
        # secretion.
        for cid in macrophage_ids:
            l_loc = world.field_mean_at_cell(il10_fi, cid)
            scale = signaling.macrophage_secretion_scale(l_loc, sig_1, g_1, g_2, d_2)
            world.set_cell_secretion_scale(chemo_fi, cid, scale)
            world.set_cell_secretion_scale(il10_fi, cid, scale)
        for cid in uninfected_ids:
            f_bar = world.field_mean_at_cell(ifn_fi, cid)
            resist = cell_resistance(f_bar, a_rf)
            scale = signaling.uninfected_il10_scale(resist)
            world.set_cell_secretion_scale(il10_fi, cid, scale)

    def _apply_killing():
        # Task 7.2: local (surface-contact) NK/CD8 kill term, see this
        # function's docstring for the full derivation/source-ordering note
        # and DISCREPANCY #7 (cell_resist applied DIRECT, not 1-resist).
        if not enable_killing or not infected_ids:
            return
        cell_volumes = world.cell_volumes()
        for cid in sorted(infected_ids):
            contact = world.cell_contact_area_by_type(cid)
            srf_nk = contact.get(types.K, 0)
            srf_cd8 = contact.get(types.E, 0)
            f_bar = world.field_mean_at_cell(ifn_fi, cid)
            resist = cell_resistance(f_bar, a_rf)
            cell_volume = cell_volumes[cid]

            rate_nk = killing.contact_kill_rate(srf_nk, resist, g_ik, tot_ec, cell_volume)
            if kill_rng.random() < 1.0 - math.exp(-rate_nk):
                world.set_cell_type(cid, types.D)
                infected_ids.discard(cid)
                continue  # source: CD8 check is skipped once NK already killed the cell this MCS

            rate_cd8 = killing.contact_kill_rate(srf_cd8, resist, g_ie, tot_ec, cell_volume)
            if kill_rng.random() < 1.0 - math.exp(-rate_cd8):
                world.set_cell_type(cid, types.D)
                infected_ids.discard(cid)

    result = {
        "steps": [], "nk_mean_distance_to_infection": [],
        "cd8_mean_distance_to_infection": [], "total_chemokine": [], "n_infected": [],
    }

    def _record(step_idx):
        coms = world.cell_coms()
        centroid = _infection_centroid(coms)
        result["steps"].append(step_idx)
        result["nk_mean_distance_to_infection"].append(
            round(_mean_distance(coms, nk_ids, centroid), 3))
        result["cd8_mean_distance_to_infection"].append(
            round(_mean_distance(coms, cd8_ids, centroid), 3))
        result["total_chemokine"].append(float(sum(world.field_conc(chemo_fi))))
        result["n_infected"].append(len(infected_ids))

    _record(0)
    for step_idx in range(1, steps + 1):
        world.step(mcs_per_update)
        _apply_secretion_scales()
        _apply_killing()
        _record(step_idx)

    result["params"] = {
        "chemotaxis_v_macro": macro_lam, "chemotaxis_v_nk": nk_lam, "chemotaxis_v_cd8": cd8_lam,
        "enable_killing": enable_killing, "g_ik": g_ik, "g_ie": g_ie, "tot_ec_ODE": tot_ec,
        "n_macrophages": n_macrophages, "n_nk": n_nk, "n_cd8": n_cd8,
        "n_infected": n_infected, "epithelial_cells_per_side": epithelial_cells_per_side,
        "margin_sites": margin_sites, "separation_sites": separation_sites,
        "seed": seed, "steps": steps, "mcs_per_update": mcs_per_update,
        "field_warmup": field_warmup,
    }
    return result


def run_global_coupling(*, num_epithelial: int | None = None, side: int = 30,
                        steps: int = 20, seed: int = 17, with_immune: bool = True,
                        seed_infection_frac: float = 0.05, mcs_per_step: int = 1,
                        with_recruitment: bool = True,
                        recruit_pool_per_type: int = 40,
                        recruit_zero_signal: bool = False,
                        enable_nearby_killing: bool = True,
                        initial_K_nb: float = 0.0, initial_E_nb: float = 0.0,
                        **scenario_kw) -> dict:
    """Task 8.2 crux driver: wire the SPATIAL->ODE direction of the Sego-2022
    hybrid coupling (dossier Sec.4a). Each MCS, advance the CPM world one step,
    read the spatial aggregates (cell counts by type + z-normalized field
    integrals + the resistance-weighted infected-load terms B_ei/G_ki), feed
    them as ``inputs`` into the Task-8.1 :class:`price_ode.GlobalODE`, integrate
    one MCS (``rr.timestep()`` cadence, dt = 60 s), and record the global
    species + spatial trajectories.

    Scenario (reused verbatim from Increments 6/7): a small non-confluent
    epithelial patch with ``n_infected`` cells seeded infected (the virus +
    type-I IFN sources) and a macrophage cluster secreting chemokine + IL-10
    (`immune.build_macrophage_scenario_spec` when ``with_immune=False``);
    ``with_immune=True`` adds the Increment-7 NK + CD8 cytotoxic cluster
    (`immune.build_cytotoxic_scenario_spec`). No infection transitions (H->I)
    are wired here; Task 8.4 wires ODE-driven recruitment (inflow/outflow of
    M/K/E, see ``with_recruitment``) and Task 8.5 (this task) wires NK/CD8
    cytotoxic killing of infected cells (LOCAL contact term, reused from
    Increment 7, PLUS the well-mixed NEARBY-population term, see
    ``enable_nearby_killing`` and `_apply_cytotoxic_killing` below) on top of
    the loop this task establishes: the spatial->ODE push, AND (Task 8.3) the
    ODE->spatial sig_1 secretion feedback. The base field secretion
    (virus/IFN from I cells, chemokine/IL-10 from M/H cells) and the
    Increment-6 IL-10-Hill macrophage secretion regulation are the existing
    physics; this driver only READS the resulting aggregates and feeds in the
    dynamic ``sig_1`` (it does not reimplement field/secretion physics or the
    Michaelis functional form -- dossier §4b(i)).

    Task 8.5 killing (each MCS, AFTER ``_apply_secretion_scales``, for every
    id still in ``infected_ids``): total cytotoxic death rate = the Task-7.2
    LOCAL contact term (`killing.contact_kill_rate`, via `world.
    cell_contact_area_by_type`, summed over NK + CD8+ contact, reused exactly
    as `run_cytotoxic_response`'s killing step) PLUS the well-mixed NEARBY
    term (`killing.nearby_kill_rate`, this task) for NK and CD8+, driven by
    the Task-8.4 nearby ODE surrogates ``state["K_nb"]``/``state["E_nb"]`` (as
    of the end of the PRECEDING MCS's recruitment step -- this MCS's own
    recruitment runs later, after the ODE step, same causal-ordering
    convention as ``sig_1_dynamic``). ``Pr = 1 - exp(-total_rate)``
    (`recruitment.ul_rate_to_prob`), ONE combined draw per infected cell (a
    documented simplification of the source's four separate per-term
    Bernoulli draws -- `ContactKillingSteppable.step()`, dossier §4b(iii) --
    adopted per the brief's exact wiring spec: summing all four rates before
    a single ``1 - exp(-total_rate)`` draw). A kill sets the cell -> ``types.
    D`` and removes it from ``infected_ids``. ``enable_nearby_killing=False``
    keeps the LOCAL term only (Increment-7 behavior) -- an on/off control
    used to isolate the NEARBY term's own contribution to clearance
    (task-8.5-report.md).

    ``initial_K_nb``/``initial_E_nb`` (both default ``0.0``, a no-op):
    override the Task-8.4 nearby-surrogate ODE state's INITIAL CONDITION
    (``state["K_nb"]``/``state["E_nb"]``, ``0.0`` at the healthy IC) --
    a documented TEST SEAM, not a source parameter. At this engine's default
    small-patch scale (``eta`` tiny), the NEARBY term is empirically too weak
    to produce an observable effect within any FAST step budget via the
    natural Task-8.4 recruitment buildup alone (mirrors the already-flagged
    ``NK_CD8_CHEMOTAXIS_ENGINE_SCALE`` weak-effect-at-default-scale finding,
    task-7.1-report.md) -- this override lets a test set a large synthetic
    nearby population directly to exercise the `killing.nearby_kill_rate`
    wiring itself (task-8.5-report.md), without tuning the rate formula or
    any recruitment constant.

    **DH-input fix (load-bearing, from the Task-8.2 review):** a cell killed
    by this step is an infected cell (I -> D), i.e. dead-from-INFECTED (the
    source's ``DI``, see `price_ode.derived_inputs`'s ``DI = tot - H - I -
    DH``) -- NOT dead-from-HEALTHY (``DH``). This driver has no H -> D
    mechanism (no ROS/Allee death wired here), so every ``types.D`` cell in
    this driver is dead-from-infected; the spatial->ODE ``DH`` input is
    therefore the raw ``types.D`` count MINUS the ids this step has killed
    (tracked in ``dead_from_infected_ids``, persists across MCS), not the raw
    count itself -- else ``DH`` would over-count and (via `price_ode.
    derived_inputs`'s ``DI = tot - H - I - DH``) ``DI`` would under-count,
    suppressing the APC-production term ``dP/dt``'s ``g_pi*DI`` driver.

    ``side`` is the nominal square epithelial-patch edge in lattice sites; the
    scenario's ``epithelial_cells_per_side`` is ``round(side/CELL_SIDE_SITES)``
    and, when ``num_epithelial is None``, the ODE reference population
    ``tot_cell`` is taken as the ACTUAL epithelial (H+I) cell count in the built
    scenario -- so ``D = tot_cell - H - I = 0`` at MCS 0 (no spurious dead
    epithelium injected into Sigma1). Defaults are deliberately small/fast;
    paper-scale is Increment 9.

    HOMEOSTATIC-SEEDING CHOICE (load-bearing carry-forward from Task 8.1): the
    ODE's healthy fixed point encodes tissue-resident immune homeostasis at
    populations ``b_m, b_k`` (~= eta * homeostatic_pops_ode). This driver feeds
    the TRUE spatial M/K/E counts (it does NOT clamp them to ``max(., b_m)``).
    At these small ``side`` defaults ``eta`` is tiny so ``b_k`` is far below the
    seeded NK count -- the initial (IC) IFN-gamma ``G`` was built with ``K=b_k``,
    so feeding the larger spatial ``K`` makes ``G`` step up over MCS 1. This is
    the expected Increment-8 coupling artifact the brief flags (recruitment in
    Task 8.4 will grow the ODE populations to match); it is recorded, not
    masked. See task-8.2-report.md.

    Returns per-MCS series (one entry per MCS, index 0 = the first stepped MCS):
      - "mcs": MCS index (1..steps)
      - "ode": {sp: [...]} for each of the 10 integrated states
        (NB, N, T, X, A, B, P, W, G, O)
      - "spatial": {"H","I","M","K","E","DH","V","F","C","L": [...]} -- the fed-in
        spatial aggregates (counts and z-normalized field integrals). Task
        9.0 §4a fix: "M"/"K"/"E" are `num_immune_by_type` = the LOCAL CPM
        count PLUS the Task-8.4 NEARBY ODE surrogate (state["M_nb"/"K_nb"/
        "E_nb"], as of the end of the preceding MCS's recruitment step) --
        this IS the value fed to the ODE, not the pure CPM count.
      - "sigma1": [a_11*T + a_12*D] each MCS (T from the ODE, D from spatial
        counts) -- Task 8.3: this is the DYNAMIC sig_1 actually fed into
        `_apply_secretion_scales` (macrophage chemokine + IL-10 Michaelis
        secretion gate, `signaling.macrophage_secretion_scale`), replacing
        the static Increment-6 stub (`params.il10.sig_1_stub`). Causal
        ordering: MCS N's secretion scale uses the sig_1 computed from MCS
        (N-1)'s freshly-integrated ODE state + spatial dead count (MCS 0 uses
        D=0, T=the healthy IC's T=0.0) -- the ODE step for MCS N itself runs
        AFTER MCS N's world.step/secretion, so same-MCS sig_1 isn't causally
        available.
      - "params": the scenario/run knobs, for the report.
    """
    params = load_params()
    a_rf = float(params["resistance"]["a_rf"])
    # Task 8.3: sig_1 (the first element) is the STATIC Increment-6 stub --
    # discarded here. Only the Hill-shape constants g_1/g_2/d_2 are reused;
    # `_apply_secretion_scales` below is fed the DYNAMIC sig_1 instead
    # (a_11*T + a_12*D, computed each MCS).
    _, g_1, g_2, d_2 = fields.il10_hill_constants()
    # Task 8.5 killing coefficients: the SAME population-independent per-MCS
    # coefficients `run_cytotoxic_response` uses for the LOCAL term (`nk.
    # g_ik`/`cd8.g_ie`), also the `g_i` the NEARBY term's `killing.
    # nearby_kill_rate` divides by `eta` (dossier §4b(iii)) -- NOT `consts
    # ["g_ki"]`/`consts["b_ei"]` below, which are the DIFFERENT eta-scaled
    # coefficients driving the Task-8.4 recruitment outflow terms G_ki/B_ei.
    g_ik = float(params["nk"]["g_ik"])
    g_ie = float(params["cd8"]["g_ie"])
    tot_ec = float(params["scaling"]["ode_epithelial_population"])

    eps = max(1, int(round(side / sheet.CELL_SIDE_SITES)))
    n_epi_cells = eps * eps
    n_infected = max(1, int(round(seed_infection_frac * n_epi_cells)))
    n_infected = min(n_infected, n_epi_cells)

    if with_immune:
        spec = immune.build_cytotoxic_scenario_spec(
            epithelial_cells_per_side=eps, n_infected=n_infected, seed=seed,
            **scenario_kw)
    else:
        spec = immune.build_macrophage_scenario_spec(
            epithelial_cells_per_side=eps, n_infected=n_infected, seed=seed,
            **scenario_kw)

    cell_type_by_idx = [c["type"] for c in spec["cells"]]  # spec index i -> cell id i+1
    infected_ids = [i + 1 for i, t in enumerate(cell_type_by_idx) if t == types.I]
    # NB: no static `macrophage_ids` here (Task 9.0 §4a fix) -- the
    # secretion-scale loop below recomputes the CURRENT macrophage id list
    # every MCS, so it picks up Task-8.4-recruited macrophages too.
    uninfected_ids = [i + 1 for i, t in enumerate(cell_type_by_idx) if t == types.H]
    n_epithelial_actual = sum(1 for t in cell_type_by_idx if t in (types.H, types.I))

    # ODE reference population (tot_cell). Default: the ACTUAL epithelial count,
    # so D = tot - H - I = 0 at MCS 0 (no spurious dead epithelium in Sigma1).
    tot_cell = int(num_epithelial) if num_epithelial is not None else n_epithelial_actual
    # eta = num_epithelial / tot_ec_ODE (`scaling.eta_0p3mm`/`eta_1p0mm`'s
    # general form, dossier §4b(iii)) -- the Task-8.5 NEARBY killing term's
    # population-scale divisor; also recorded in "params" below (unchanged
    # value/definition from before this task).
    eta = tot_cell / price_ode.ODE_EPITHELIAL_POPULATION

    world = build.world_from_spec(spec, finalize=False)
    virus_fi = fields.add_virus_field(world)
    ifn_fi = fields.add_ifn_field(world)
    chemo_fi = fields.add_chemokine_field(world)
    il10_fi = fields.add_il10_field(world)

    # Task 8.4: pre-seed the dormant recruitment reserve pool BEFORE finalize
    # (the only time the engine allows new cells -- see RECRUIT_RESERVE_TYPE).
    # Targets = the immune types the scenario actually builds (macrophage-only
    # scenario has no NK/CD8 clusters, so it pools only macrophages).
    recruit_targets = [types.M, types.K, types.E] if with_immune else [types.M]
    reserve_pool: dict[int, list[int]] = {}
    if with_recruitment:
        mv = int(params["macrophage"]["volume_sites"])
        mlv = float(params["macrophage"]["lambda_volume"])
        # Reserve cells are inert loners: give them a modest Medium adhesion so
        # they hold shape and stay parked (they never chemotax/secrete while
        # dormant). One J vs each existing type + self is enough.
        for t in range(0, RECRUIT_RESERVE_TYPE):
            world.set_contact(RECRUIT_RESERVE_TYPE, t, 10.0)
        world.set_contact(RECRUIT_RESERVE_TYPE, RECRUIT_RESERVE_TYPE, 25.0)
        reserve_pool = _seed_recruit_pool(
            world, spec, pool_per_type=recruit_pool_per_type, targets=recruit_targets,
            target_volume=mv, lambda_volume=mlv)

    world.finalize(int(spec["potts"]["seed"]))

    dim_z = world.dims()[2]  # z=1 in this reproduction; §4a divides integrals by dim.z

    # Task-8.1 ODE: resolve constants at this tot_cell, seed the healthy IC.
    consts = price_ode.resolve_constants(params["price_ode"], num_epithelial=tot_cell)
    ode = price_ode.GlobalODE(consts, num_epithelial=tot_cell)
    state = price_ode.initial_state(consts, num_epithelial=tot_cell, v0=0.0)
    # Task 8.5 test seam (see docstring, "initial_K_nb/initial_E_nb"
    # paragraph): both default 0.0 (== `price_ode.initial_state`'s own
    # healthy-IC value, a no-op).
    if initial_K_nb:
        state["K_nb"] = float(initial_K_nb)
    if initial_E_nb:
        state["E_nb"] = float(initial_E_nb)
    b_ei = consts["b_ei"]
    g_ki = consts["g_ki"]
    a_11, a_12 = consts["a_11"], consts["a_12"]

    # Task 8.4 recruitment: a dedicated seeded RNG stream (offset unused by any
    # other stream in this module) keeps the inflow Poisson draws + outflow
    # Bernoulli draws deterministic and independent of the Potts RNG.
    recruit_rng = np.random.default_rng(seed + 900)
    _lr = params["coupling"]["recruitment"]["local_ratios"]
    local_ratio = {types.M: float(_lr["macro"]), types.K: float(_lr["nk"]),
                   types.E: float(_lr["cd8"])}
    # (target_type, inflow_fn(driver, consts), outflow_fn(load, consts),
    #  driver-picker, load-picker, nearby-surrogate key) -- discrepancy #12:
    # macrophage/NK are chemokine(C)-driven, CD8 is APC(P)-driven; CD8 inflow
    # has NO homeostatic baseline. Drivers are read fresh each MCS below.
    def _recruit_step(C_field, P_ode, G_ki_now, B_ei_now):
        """Run one MCS of ODE-driven recruitment (AFTER the ODE step, per the
        source's post-`rr.timestep()` `update_populations`). Mutates the world
        (activate reserves / remove cells) and the ODE nearby surrogates in
        `state`; returns the six per-type inflow/outflow rates for recording.

        `recruit_zero_signal` forces the Hill DRIVER fields (C, P) to 0 -- the
        no-signal control that isolates each law's signal-driven Hill term from
        its (signal-independent) homeostatic baseline. G_ki/B_ei (infected-load
        outflow terms) are spatial state, not the recruitment signal, so they
        are NOT zeroed."""
        C_sig = 0.0 if recruit_zero_signal else C_field
        P_sig = 0.0 if recruit_zero_signal else P_ode

        rates = {
            "macro_inflow": recruitment.macrophage_inflow(C_sig, consts),
            "nk_inflow": recruitment.nk_inflow(C_sig, consts),
            "cd8_inflow": recruitment.cd8_inflow(P_sig, consts),
            "macro_outflow": recruitment.macrophage_outflow(consts),
            "nk_outflow": recruitment.nk_outflow(G_ki_now, consts),
            "cd8_outflow": recruitment.cd8_outflow(B_ei_now, consts),
        }
        by_type = {
            types.M: (rates["macro_inflow"], rates["macro_outflow"], "M_nb"),
            types.K: (rates["nk_inflow"], rates["nk_outflow"], "K_nb"),
            types.E: (rates["cd8_inflow"], rates["cd8_outflow"], "E_nb"),
        }
        for t, (inflow, outflow, nb_key) in by_type.items():
            lr = local_ratio[t]
            # --- LOCAL fraction -> CPM cells --------------------------------
            # Outflow: remove each active cell with prob ul_rate_to_prob(lr*out)
            # (source `outflow_by_type`). Removed cells' voxels go to Medium AND
            # are relabelled to the reserve type so they stop counting/secreting.
            pr_out = recruitment.ul_rate_to_prob(lr * outflow)
            if pr_out > 0.0:
                active = [cid for cid, ct in enumerate(world.cell_types())
                          if ct == t and cid != 0]
                to_remove = [cid for cid in active if recruit_rng.random() < pr_out]
                if to_remove:
                    world.remove_cells(to_remove)
                    for cid in to_remove:
                        world.set_cell_type(cid, RECRUIT_RESERVE_TYPE)
            # Inflow: Poisson draw at lr*inflow, activate that many reserves
            # (capped by the pool; extra draws fail, matching the source's
            # `try_add` failure when no Medium/cell is available).
            n_local = recruitment.poisson_inflow_count(lr * inflow, recruit_rng)
            pool_t = reserve_pool.get(t, [])
            for _ in range(n_local):
                if not pool_t:
                    break
                world.set_cell_type(pool_t.pop(0), t)
            # --- NEARBY fraction -> ODE surrogate (M_nb/K_nb/E_nb) ----------
            # (1-lr) of the inflow accrues into the well-mixed nearby surrogate
            # rather than placing a CPM cell (dossier Sec.4b(ii) / brief); it
            # attrits at the full outflow rate. Continuous (these are ODE
            # populations). Macrophage lr=1.0 -> no nearby accrual.
            nearby_in = (1.0 - lr) * inflow
            state[nb_key] = state[nb_key] * (1.0 - recruitment.ul_rate_to_prob(outflow)) + nearby_in
        return rates

    def _apply_secretion_scales(sig_1_dynamic):
        # Increment-6 IL-10-Hill macrophage secretion regulation (existing
        # physics, same as run_macrophage_signaling), now driven by the
        # Task-8.3 DYNAMIC sig_1 = a_11*T + a_12*D (ODE TNF + spatial dead
        # count) instead of the static Increment-6 stub: local IL-10
        # self-limits the macrophage's chemokine + IL-10 release; uninfected
        # cells' IL-10 is (1-resist)-gated (unaffected by sig_1).
        # Task 9.0 §4a fix: refresh the macrophage-id list EVERY MCS (any
        # cell currently typed `types.M`, seeded OR Task-8.4-recruited),
        # not the static pre-loop `macrophage_ids` -- else a macrophage
        # activated from the reserve pool by `_recruit_step` never gets its
        # chemokine/IL-10 secretion scale set. Same "current cells of a
        # type" pattern `_recruit_step`'s own outflow draw already uses.
        current_macrophage_ids = [cid for cid, ct in enumerate(world.cell_types())
                                   if ct == types.M and cid != 0]
        for cid in current_macrophage_ids:
            l_loc = world.field_mean_at_cell(il10_fi, cid)
            scale = signaling.macrophage_secretion_scale(l_loc, sig_1_dynamic, g_1, g_2, d_2)
            world.set_cell_secretion_scale(chemo_fi, cid, scale)
            world.set_cell_secretion_scale(il10_fi, cid, scale)
        for cid in uninfected_ids:
            f_bar = world.field_mean_at_cell(ifn_fi, cid)
            resist = cell_resistance(f_bar, a_rf)
            world.set_cell_secretion_scale(il10_fi, cid, signaling.uninfected_il10_scale(resist))

    # Task 8.5 kill-draw RNG: a dedicated stream/offset (unused by any other
    # RNG in this function -- `recruit_rng` above uses seed+900) so kill
    # draws don't correlate with the recruitment inflow/outflow draws or the
    # Potts RNG (set at `world.finalize`).
    kill_rng = np.random.default_rng(seed + 800)
    # DH-input fix (see this function's docstring, "DH-input fix" paragraph):
    # ids of cells killed BY THIS STEP (dead-from-INFECTED, source `DI`), kept
    # so the spatial->ODE `DH` (dead-from-HEALTHY) input below can exclude
    # them -- this driver has no other death mechanism, so every `types.D`
    # cell is one of these.
    dead_from_infected_ids: set[int] = set()

    def _apply_cytotoxic_killing():
        """Task 8.5: for every id still in `infected_ids`, combine the
        Task-7.2 LOCAL contact term (`killing.contact_kill_rate`, NK + CD8+,
        via `world.cell_contact_area_by_type` -- reused exactly as
        `run_cytotoxic_response`'s killing step) with the well-mixed NEARBY
        term (`killing.nearby_kill_rate`, this task, using `state["K_nb"]`/
        `state["E_nb"]`, gated by `enable_nearby_killing`) into ONE total
        rate, draw death with `Pr = 1 - exp(-total_rate)`
        (`recruitment.ul_rate_to_prob`), and on a kill: `world.set_cell_type
        (cid, types.D)`, drop `cid` from `infected_ids`, record it in
        `dead_from_infected_ids`. See this function's module-level docstring
        for the full derivation/causal-ordering note."""
        if not infected_ids:
            return
        cell_volumes = world.cell_volumes()
        for cid in list(infected_ids):
            contact = world.cell_contact_area_by_type(cid)
            srf_nk = contact.get(types.K, 0)
            srf_cd8 = contact.get(types.E, 0)
            f_bar = world.field_mean_at_cell(ifn_fi, cid)
            resist = cell_resistance(f_bar, a_rf)
            cell_volume = cell_volumes[cid]

            total_rate = (
                killing.contact_kill_rate(srf_nk, resist, g_ik, tot_ec, cell_volume)
                + killing.contact_kill_rate(srf_cd8, resist, g_ie, tot_ec, cell_volume)
            )
            if enable_nearby_killing:
                total_rate += (
                    killing.nearby_kill_rate(g_ik, eta, state["K_nb"], resist)
                    + killing.nearby_kill_rate(g_ie, eta, state["E_nb"], resist)
                )

            if kill_rng.random() < 1.0 - math.exp(-total_rate):
                world.set_cell_type(cid, types.D)
                infected_ids.remove(cid)
                dead_from_infected_ids.add(cid)

    result = {
        "mcs": [],
        "ode": {sp: [] for sp in price_ode.INTEGRATED_STATES},
        "spatial": {k: [] for k in ("H", "I", "M", "K", "E", "DH", "V", "F", "C", "L")},
        "sigma1": [],
        # Task 8.4: per-MCS recruited populations (LOCAL CPM counts M/K/E, NEARBY
        # ODE surrogates M_nb/K_nb/E_nb) + the driving inflow rates, for the
        # report/trajectory. Absent when with_recruitment=False.
        "recruit": ({"M": [], "K": [], "E": [], "M_nb": [], "K_nb": [], "E_nb": [],
                     "macro_inflow_rate": [], "nk_inflow_rate": [], "cd8_inflow_rate": [],
                     "macro_outflow_rate": [], "nk_outflow_rate": [], "cd8_outflow_rate": []}
                    if with_recruitment else None),
    }

    # Task-8.3: sig_1 seed BEFORE MCS 1 -- D = tot_cell - H - I = 0 at MCS 0
    # (the same invariant the module docstring already documents for the
    # recorded sigma1 series), T = the healthy-IC ODE state's T (0.0 per
    # §2.5). This is the dynamic sig_1 the FIRST `_apply_secretion_scales`
    # call below uses; each subsequent MCS uses the sig_1 computed from the
    # PRECEDING MCS's freshly-integrated ODE state + spatial dead count (the
    # only causal ordering available: this MCS's secretion is set before this
    # MCS's own ODE step runs).
    sig_1_dynamic = a_11 * state["T"] + a_12 * 0.0

    for mcs in range(1, steps + 1):
        world.step(mcs_per_step)
        _apply_secretion_scales(sig_1_dynamic)
        # Task 8.5: NK/CD8 cytotoxic killing (LOCAL + NEARBY), BEFORE reading
        # this MCS's spatial aggregates below, so I/DH reflect any kills that
        # just happened this MCS (see this function's docstring).
        _apply_cytotoxic_killing()

        # --- read spatial aggregates (§4a) ---
        types_now = world.cell_types()
        H = sum(1 for t in types_now[1:] if t == types.H)
        I = sum(1 for t in types_now[1:] if t == types.I)
        M = sum(1 for t in types_now[1:] if t == types.M)
        K = sum(1 for t in types_now[1:] if t == types.K)
        E = sum(1 for t in types_now[1:] if t == types.E)
        # Task 8.5 DH-input fix: raw `types.D` count MINUS the ids Task 8.5's
        # killing has itself killed (dead-from-INFECTED, source `DI`) -- this
        # driver has no other death mechanism, so every `types.D` cell here IS
        # one of `dead_from_infected_ids`, and DH (dead-from-HEALTHY) stays 0.
        # See this function's docstring, "DH-input fix" paragraph.
        D_total = sum(1 for t in types_now[1:] if t == types.D)
        DH = D_total - len(dead_from_infected_ids)

        V = float(sum(world.field_conc(virus_fi))) / dim_z
        F = float(sum(world.field_conc(ifn_fi))) / dim_z
        C = float(sum(world.field_conc(chemo_fi))) / dim_z
        L = float(sum(world.field_conc(il10_fi))) / dim_z

        # resistance-weighted infected load: B_ei = Σ_infected(resist)*b_ei,
        # G_ki = Σ_infected(resist)*g_ki (§4a, ImmuneModelSteppable :1300-1306).
        sum_resist = 0.0
        for cid in infected_ids:
            f_bar = world.field_mean_at_cell(ifn_fi, cid)
            sum_resist += cell_resistance(f_bar, a_rf)
        B_ei = sum_resist * b_ei
        G_ki = sum_resist * g_ki

        # Task 9.0 §4a fix: ODE `num_immune_by_type` = LOCAL CPM count +
        # NEARBY ODE surrogate (dossier §4a), not the pure CPM count --
        # `state["M_nb"/"K_nb"/"E_nb"]` as of the end of the PRECEDING MCS's
        # `_recruit_step` (same causal-ordering convention
        # `_apply_cytotoxic_killing`'s NEARBY term already uses). Macro
        # `M_nb` stays 0 (macro `local_ratio`=1.0 -> no nearby accrual, see
        # `_recruit_step`); NK/CD8 pick up their Task-8.4 nearby surrogates.
        M_ode = M + state["M_nb"]
        K_ode = K + state["K_nb"]
        E_ode = E + state["E_nb"]

        inputs = dict(H=H, I=I, M=M_ode, K=K_ode, E=E_ode, DH=DH,
                      V=V, F=F, C=C, L=L, B_ei=B_ei, G_ki=G_ki)

        # --- integrate one MCS (dt = 60 s) ---
        state = ode.step(state, inputs, dt_seconds=60.0)

        # --- Task 8.4: ODE-driven recruitment (AFTER the ODE step, using the
        # freshly-integrated APC P + this MCS's chemokine field C and infected-
        # load terms) -- seeds/removes CPM immune cells + updates nearby
        # surrogates. Its effect on the M/K/E counts is seen by the NEXT MCS's
        # spatial->ODE push (source ordering: update_populations follows
        # rr.timestep()). ---
        if with_recruitment:
            rates = _recruit_step(C, state["P"], G_ki, B_ei)
            types_after = world.cell_types()
            result["recruit"]["M"].append(sum(1 for t in types_after[1:] if t == types.M))
            result["recruit"]["K"].append(sum(1 for t in types_after[1:] if t == types.K))
            result["recruit"]["E"].append(sum(1 for t in types_after[1:] if t == types.E))
            result["recruit"]["M_nb"].append(state["M_nb"])
            result["recruit"]["K_nb"].append(state["K_nb"])
            result["recruit"]["E_nb"].append(state["E_nb"])
            result["recruit"]["macro_inflow_rate"].append(rates["macro_inflow"])
            result["recruit"]["nk_inflow_rate"].append(rates["nk_inflow"])
            result["recruit"]["cd8_inflow_rate"].append(rates["cd8_inflow"])
            result["recruit"]["macro_outflow_rate"].append(rates["macro_outflow"])
            result["recruit"]["nk_outflow_rate"].append(rates["nk_outflow"])
            result["recruit"]["cd8_outflow_rate"].append(rates["cd8_outflow"])

        # Task 8.3: sig_1 = a_11*T + a_12*D (D from spatial counts, T from the
        # ODE state just integrated) -- recorded AND carried into next MCS's
        # `_apply_secretion_scales` call above, replacing the static
        # Increment-6 stub as the macrophage chemokine/IL-10 secretion driver.
        D = tot_cell - H - I
        sig_1_dynamic = a_11 * state["T"] + a_12 * D
        sigma1 = sig_1_dynamic

        result["mcs"].append(mcs)
        for sp in price_ode.INTEGRATED_STATES:
            result["ode"][sp].append(state[sp])
        for k, v in inputs.items():
            if k in result["spatial"]:
                result["spatial"][k].append(v)
        result["sigma1"].append(sigma1)

    result["params"] = {
        "side": side, "epithelial_cells_per_side": eps,
        "num_epithelial": tot_cell, "n_epithelial_actual": n_epithelial_actual,
        "n_infected": n_infected, "seed_infection_frac": seed_infection_frac,
        "with_immune": with_immune, "steps": steps, "seed": seed,
        "mcs_per_step": mcs_per_step, "eta": eta,
        "b_m": consts["b_m"], "b_k": consts["b_k"], "b_p": consts["b_p"],
        "with_recruitment": with_recruitment,
        "recruit_pool_per_type": recruit_pool_per_type if with_recruitment else 0,
        "recruit_zero_signal": recruit_zero_signal,
        "recruit_local_ratios": {"macro": local_ratio[types.M], "nk": local_ratio[types.K],
                                 "cd8": local_ratio[types.E]},
        "enable_nearby_killing": enable_nearby_killing,
        "n_infected_final": len(infected_ids),
        "dead_from_infected": len(dead_from_infected_ids),
    }
    return result


# Enable-flag tokens accepted by `run_full_model.enable` (subsystem isolation).
_FULL_MODEL_SUBSYSTEMS = (
    "infection", "ifn", "death", "allee", "macrophage", "chemokine",
    "nk_cd8", "killing", "recruitment", "ode",
)


def _seed_uniform_virus(world, virus_fi, target_conc, epithelial_ids):
    """Approximate the source's uniform ``Virus``-field initial condition
    (``ImmuneModelSteppable.init_fresh_immune_model``: if ``v0 > 0`` seed the
    field uniformly to ``v0`` with NO pre-infected cells) within this engine's
    API, which has NO direct field-write primitive (`crates/cpm-py/src/lib.rs`
    exposes only secretion + diffusion, no field setter) and MUST NOT be
    rebuilt this task.

    Engine-only approximation: transiently make every epithelial (H) cell a
    virus SOURCE for one field advance, sized so the per-pixel deposit equals
    ``target_conc`` (the field engine adds ``rate*dt`` once per advance, see
    `fields._per_pixel_secretion_rate`), then RESET H-cell virus secretion to 0
    so only infected cells secrete during the run proper. The result is a
    ~uniform ``target_conc`` virus concentration under the epithelial patch --
    the region where infection is decided -- with a small amount of edge
    diffusion (virus D ~ 0.18 lat^2/MCS is the smallest of any field, so the
    spread over one advance is negligible). This is a documented IC
    approximation, NOT a tuned rate: the deposit scales linearly with
    ``target_conc`` (so the fig5 viral-load dose-response is monotone) and adds
    no multiplier to any dynamics constant.
    """
    _d, _decay, dt, _substeps, _sec = fields._virus_field_params()
    # `world.set_secretion(field, TYPE, rate)` is a per-TYPE rule; every H cell
    # is the epithelial patch, so this deposits uniformly across the patch.
    world.set_secretion(virus_fi, types.H, float(target_conc) / dt)
    world.advance_fields(1)
    world.set_secretion(virus_fi, types.H, 0.0)


def run_full_model(*, cells_per_side: int, steps: int, seed: int,
                   init_infection_frac: float | None = None,
                   init_viral_load: float | None = None,
                   s_per_mcs: float = 60.0,
                   enable=_FULL_MODEL_SUBSYSTEMS,
                   mcs_per_step: int = 7,
                   n_macrophages: int = 4, n_nk: int = 4, n_cd8: int = 4,
                   margin_sites: int = 10, separation_sites: int = 8,
                   recruit_pool_per_type: int = 6,
                   **scenario_kw) -> dict:
    """Increment 9 Task 9.1 -- the CAPSTONE multiscale driver: assemble ALL
    mechanism primitives from Increments 0-8 into ONE per-MCS loop, in the
    Sego-2022 source's steppable order (dossier
    ``docs/cc3d-reference/sego2022-global-ode.md`` §4 / the plan's Global
    Constraints), and record every capstone observable.

    This is `run_global_coupling`'s ODE-coupling block (spatial->ODE push,
    dynamic ``sig_1`` macrophage secretion, ODE-driven recruitment, NK/CD8
    LOCAL+NEARBY killing, the §4a nearby-surrogate + recruited-macrophage
    fidelity fixes) reused wholesale, PLUS the two subsystems it was missing
    per the Task-8.2 review: the epithelial-fate TRANSITIONS (infection H->I,
    infected apoptosis I->D, Allee death H->D / recovery D->H) and immune-cell
    CHEMOTAXIS (macrophage up virus, NK/CD8 up chemokine). No mechanism math is
    reimplemented here -- every step calls the existing Increment-0..8 helper.

    Scenario: `immune.build_cytotoxic_scenario_spec` at
    ``epithelial_cells_per_side = cells_per_side`` (a non-confluent epithelial
    patch + interior macrophage and NK/CD8 clusters + the Task-8.4 dormant
    recruit reserve pool). ``cells_per_side`` sets the ODE population scale
    ``eta = cells_per_side**2 / 250000`` (the ODE reference population
    ``tot_cell`` is the actual epithelial cell count, ``cells_per_side**2``).

    Scenario selection (mutually exclusive -- exactly one required):
      - ``init_infection_frac``: seed ``round(frac*tot_cell)`` random-nearest-
        center epithelial cells as INFECTED (the virus/IFN sources), no initial
        virus field. ``infected[0] > 0``.
      - ``init_viral_load``: seed NO pre-infected cells; instead lay down a
        ~uniform virus-field IC at concentration ``init_viral_load`` under the
        patch (`_seed_uniform_virus`), so infection EMERGES from the field
        (H->I ∝ local virus·g_hv). ``infected[0] == 0``.

    Per-MCS SOURCE-ORDERED pipeline (Global Constraints; each step calls the
    named existing helper, gated by ``enable`` so tests can isolate a
    subsystem):
      1. ``world.step`` -- advance Potts + all four fields (virus/IFN/chemokine/
         IL-10) together (also applies the chemotaxis wired once at setup:
         macrophage↑virus, NK/CD8↑chemokine).
      2. infection H->I -- `transitions.infection_step` (local virus·``g_hv``).
      3. per-cell resistance ρ from local IFN -- `resistance.cell_resistance`.
      4. field secretion scales -- virus release by infected gated ``(1-ρ)``;
         macrophage chemokine/IL-10 via `signaling.macrophage_secretion_scale`
         driven by the DYNAMIC ``sig_1 = a_11*T + a_12*D`` (§4b(i)), macrophage
         id list refreshed each MCS (§4a fix); uninfected IL-10 gated ``(1-ρ)``.
      5. chemotaxis -- persistent engine setting wired once at setup (applied
         inside every ``world.step``); no per-MCS action.
      6. contact + nearby killing I->D -- `killing.contact_kill_rate` (LOCAL) +
         `killing.nearby_kill_rate` (well-mixed NEARBY, from ``K_nb``/``E_nb``).
      7. Allee death H->D / recovery D->H -- `allee.allee_death_rate` /
         `allee.allee_recovery_rate` over the epithelial contact geometry.
      8. infected apoptosis I->D -- `transitions.infected_death_step`
         (``mu_i·(1-ρ)``).
      9. spatial->ODE push + `price_ode.GlobalODE.step` + ODE->spatial
         recruitment -- the source's ``__update_spatial_data`` ->
         ``rr.timestep()`` -> ``update_populations`` order (§4; recruitment
         follows the integrate, using the freshly-integrated APC ``P``), reused
         verbatim from `run_global_coupling`.

    CALIBRATION HONESTY (plan Global Constraints): this task WIRES the model; it
    does NOT tune to hit figure bands. The only engine-unit scale used is the
    documented Increment-7 100× chemotaxis normalization
    (`NK_CD8_CHEMOTAXIS_ENGINE_SCALE`, see its module-level docstring -- the
    source's literal λ is statistically undetectable at this engine's
    chemokine-field concentration scale); source constants and the correct η
    are used unchanged. At the reduced η of these fast tests some subsystems
    (recruitment inflow, NEARBY killing) are near-inert -- expected, flagged for
    the 9.2-9.4 / Phase-B paper-scale calibration, and NOT tuned here.

    Returns (index 0 = the seeded INITIAL state, before any transition --
    matching `run_virus_infection`'s convention, so ``infected[0]`` reflects the
    scenario's initial seeding; indices 1..steps-1 are after each stepped MCS,
    for ``steps`` records total; ``t_days[i] = i*s_per_mcs/86400``)::

        {
          "mcs":    [0, 1, ..., steps-1],
          "t_days": [0.0, ...],
          "counts": {"uninfected","infected","dead","macrophage","nk","cd8"},
          "fields": {"virus","ifn","chemo","il10"},   # lattice field integrals
          "ode":    {sp: [...] for the 10 integrated states}
                    + {"H","I","M","K","E","DH"}       # the spatial->ODE inputs
                    (M/K/E = LOCAL CPM count + NEARBY surrogate, the §4a value
                    actually fed to the ODE),
          "params": {...},
        }
    """
    if (init_infection_frac is None) == (init_viral_load is None):
        raise ValueError(
            "run_full_model requires exactly one of init_infection_frac XOR "
            "init_viral_load (got both or neither)")

    enable = set(enable)
    params = load_params()
    a_rf = float(params["resistance"]["a_rf"])
    _, g_1, g_2, d_2 = fields.il10_hill_constants()
    g_ik = float(params["nk"]["g_ik"])
    g_ie = float(params["cd8"]["g_ie"])
    tot_ec = float(params["scaling"]["ode_epithelial_population"])
    g_hv = float(params["virus"]["infection_g_hv"])
    mu_i = float(params["cell_death"]["mu_i_per_mcs"])
    b_h = float(params["allee"]["b_h"])
    theta = float(params["allee"]["srf_threshold"])

    # --- scenario ---------------------------------------------------------
    tot_cell = cells_per_side * cells_per_side
    if init_infection_frac is not None:
        n_infected = int(round(float(init_infection_frac) * tot_cell))
        n_infected = max(0, min(n_infected, tot_cell))
    else:
        n_infected = 0  # init_viral_load: no pre-infected cells

    spec = immune.build_cytotoxic_scenario_spec(
        epithelial_cells_per_side=cells_per_side, n_infected=n_infected,
        n_macrophages=n_macrophages, n_nk=n_nk, n_cd8=n_cd8,
        margin_sites=margin_sites, separation_sites=separation_sites, seed=seed,
        **scenario_kw)

    world = build.world_from_spec(spec, finalize=False)
    virus_fi = fields.add_virus_field(world)
    ifn_fi = fields.add_ifn_field(world)
    chemo_fi = fields.add_chemokine_field(world)
    il10_fi = fields.add_il10_field(world)

    # Task-8.4 recruit reserve pool (dormant cells activated by ODE inflow).
    recruit_targets = [types.M, types.K, types.E]
    reserve_pool: dict[int, list[int]] = {}
    if "recruitment" in enable:
        mv = int(params["macrophage"]["volume_sites"])
        mlv = float(params["macrophage"]["lambda_volume"])
        for t in range(0, RECRUIT_RESERVE_TYPE):
            world.set_contact(RECRUIT_RESERVE_TYPE, t, 10.0)
        world.set_contact(RECRUIT_RESERVE_TYPE, RECRUIT_RESERVE_TYPE, 25.0)
        reserve_pool = _seed_recruit_pool(
            world, spec, pool_per_type=recruit_pool_per_type,
            targets=recruit_targets, target_volume=mv, lambda_volume=mlv)

    world.finalize(int(spec["potts"]["seed"]))

    n_cells = len(spec["cells"])
    dim_z = world.dims()[2]

    cell_type_by_idx = [c["type"] for c in spec["cells"]]
    epithelial_ids = [i + 1 for i, t in enumerate(cell_type_by_idx)
                      if t in (types.H, types.I)]

    # init_viral_load: lay down the ~uniform virus IC (engine has no field-write
    # primitive; see `_seed_uniform_virus`) BEFORE the initial-state record.
    if init_viral_load is not None:
        _seed_uniform_virus(world, virus_fi, init_viral_load, epithelial_ids)

    # --- chemotaxis (the piece run_global_coupling was missing) -----------
    if "macrophage" in enable:
        immune.set_macrophage_chemotaxis(world, virus_fi)
    if "nk_cd8" in enable:
        nk_lam = float(params["nk"]["chemotaxis_v_nk"]) * NK_CD8_CHEMOTAXIS_ENGINE_SCALE
        cd8_lam = float(params["cd8"]["chemotaxis_v_cd8"]) * NK_CD8_CHEMOTAXIS_ENGINE_SCALE
        immune.set_nk_cd8_chemotaxis(world, chemo_fi,
                                     chemotaxis_v_nk=nk_lam, chemotaxis_v_cd8=cd8_lam)

    # --- ODE (Task 8.1) ---------------------------------------------------
    consts = price_ode.resolve_constants(params["price_ode"], num_epithelial=tot_cell)
    ode = price_ode.GlobalODE(consts, num_epithelial=tot_cell)
    state = price_ode.initial_state(consts, num_epithelial=tot_cell, v0=0.0)
    b_ei = consts["b_ei"]
    g_ki = consts["g_ki"]
    a_11, a_12 = consts["a_11"], consts["a_12"]
    eta = tot_cell / price_ode.ODE_EPITHELIAL_POPULATION

    # --- independent RNG streams (offsets match the other drivers) --------
    infect_rng = np.random.default_rng(seed + 1)
    death_rng = np.random.default_rng(seed + 2)
    allee_rng = np.random.default_rng(seed + 3)
    kill_rng = np.random.default_rng(seed + 800)
    recruit_rng = np.random.default_rng(seed + 900)

    _lr = params["coupling"]["recruitment"]["local_ratios"]
    local_ratio = {types.M: float(_lr["macro"]), types.K: float(_lr["nk"]),
                   types.E: float(_lr["cd8"])}

    infected_ids = set(i + 1 for i, t in enumerate(cell_type_by_idx) if t == types.I)
    # Dead-cell provenance for the §4a DH input: Allee/ROS deaths are
    # dead-from-HEALTHY (source `DH`); killing + apoptosis are dead-from-
    # INFECTED (source `DI = tot-H-I-DH`). DH = count of currently-D cells of
    # healthy origin.
    dead_from_healthy_ids: set[int] = set()
    dead_from_infected_ids: set[int] = set()

    def _recruit_step(C_field, P_ode, G_ki_now, B_ei_now):
        """One MCS of ODE-driven recruitment (source `update_populations`,
        AFTER `rr.timestep()`): activate/remove CPM immune cells + attrit the
        nearby surrogates. Reused verbatim from `run_global_coupling`."""
        rates = {
            "macro_inflow": recruitment.macrophage_inflow(C_field, consts),
            "nk_inflow": recruitment.nk_inflow(C_field, consts),
            "cd8_inflow": recruitment.cd8_inflow(P_ode, consts),
            "macro_outflow": recruitment.macrophage_outflow(consts),
            "nk_outflow": recruitment.nk_outflow(G_ki_now, consts),
            "cd8_outflow": recruitment.cd8_outflow(B_ei_now, consts),
        }
        by_type = {
            types.M: (rates["macro_inflow"], rates["macro_outflow"], "M_nb"),
            types.K: (rates["nk_inflow"], rates["nk_outflow"], "K_nb"),
            types.E: (rates["cd8_inflow"], rates["cd8_outflow"], "E_nb"),
        }
        for t, (inflow, outflow, nb_key) in by_type.items():
            lr = local_ratio[t]
            pr_out = recruitment.ul_rate_to_prob(lr * outflow)
            if pr_out > 0.0:
                active = [cid for cid, ct in enumerate(world.cell_types())
                          if ct == t and cid != 0]
                to_remove = [cid for cid in active if recruit_rng.random() < pr_out]
                if to_remove:
                    world.remove_cells(to_remove)
                    for cid in to_remove:
                        world.set_cell_type(cid, RECRUIT_RESERVE_TYPE)
            n_local = recruitment.poisson_inflow_count(lr * inflow, recruit_rng)
            pool_t = reserve_pool.get(t, [])
            for _ in range(n_local):
                if not pool_t:
                    break
                world.set_cell_type(pool_t.pop(0), t)
            nearby_in = (1.0 - lr) * inflow
            state[nb_key] = state[nb_key] * (1.0 - recruitment.ul_rate_to_prob(outflow)) + nearby_in
        return rates

    def _counts(types_now):
        H = sum(1 for t in types_now[1:] if t == types.H)
        I = sum(1 for t in types_now[1:] if t == types.I)
        D = sum(1 for t in types_now[1:] if t == types.D)
        M = sum(1 for t in types_now[1:] if t == types.M)
        K = sum(1 for t in types_now[1:] if t == types.K)
        E = sum(1 for t in types_now[1:] if t == types.E)
        return H, I, D, M, K, E

    result = {
        "mcs": [], "t_days": [],
        "counts": {k: [] for k in ("uninfected", "infected", "dead",
                                   "macrophage", "nk", "cd8")},
        "fields": {k: [] for k in ("virus", "ifn", "chemo", "il10")},
        # `ode` carries BOTH the 10 integrated states AND the §4a spatial->ODE
        # inputs (H,I,M,K,E,DH) actually fed to the ODE -- M/K/E = LOCAL CPM
        # count + NEARBY surrogate, so e.g. r["ode"]["K"] >= the local NK count.
        "ode": {sp: [] for sp in
                tuple(price_ode.INTEGRATED_STATES) + ("H", "I", "M", "K", "E", "DH")},
    }

    def _record(i):
        types_now = world.cell_types()
        H, I, D, M, K, E = _counts(types_now)
        DH = sum(1 for cid in dead_from_healthy_ids if types_now[cid] == types.D)
        mcs_now = i * mcs_per_step  # actual MCS elapsed at record i
        result["mcs"].append(mcs_now)
        result["t_days"].append(mcs_now * s_per_mcs / 86400.0)
        result["counts"]["uninfected"].append(H)
        result["counts"]["infected"].append(I)
        result["counts"]["dead"].append(D)
        result["counts"]["macrophage"].append(M)
        result["counts"]["nk"].append(K)
        result["counts"]["cd8"].append(E)
        result["fields"]["virus"].append(float(sum(world.field_conc(virus_fi))))
        result["fields"]["ifn"].append(float(sum(world.field_conc(ifn_fi))))
        result["fields"]["chemo"].append(float(sum(world.field_conc(chemo_fi))))
        result["fields"]["il10"].append(float(sum(world.field_conc(il10_fi))))
        for sp in price_ode.INTEGRATED_STATES:
            result["ode"][sp].append(state[sp])
        result["ode"]["H"].append(H)
        result["ode"]["I"].append(I)
        result["ode"]["M"].append(M + state["M_nb"])
        result["ode"]["K"].append(K + state["K_nb"])
        result["ode"]["E"].append(E + state["E_nb"])
        result["ode"]["DH"].append(DH)

    # resist_at_cell persists across the inner-MCS loop and is read by the
    # per-record ODE push (B_ei/G_ki). Recomputed every MCS below.
    resist_at_cell = [0.0] * (n_cells + 1)

    def _one_mcs(sig_1_dynamic):
        """One source-faithful MCS of the cell-scale pipeline (steps 1-8 of the
        dossier §4 order): advance Potts + all four fields, then infection,
        resistance, field-secretion scales, contact+nearby killing, Allee
        death/recovery, and infected apoptosis. Every cell-fate transition
        fires EVERY MCS (matching the source's per-MCS steppable cadence). The
        systemic-ODE coupling (step 9) runs once per RECORDED step in the outer
        loop -- a documented coupling-cadence reduction for test speed that is
        immaterial at reduced scale (the integrated ODE species evolve slowly);
        ``sig_1`` and the nearby surrogates it feeds are refreshed there."""
        nonlocal resist_at_cell
        # (1) advance Potts + all fields (chemotaxis applied inside world.step).
        world.step(1)
        current_types = list(world.cell_types())

        # (2) infection H->I (local virus·g_hv).
        if "infection" in enable:
            virus_at_cell = [world.field_mean_at_cell(virus_fi, cid)
                             for cid in range(n_cells + 1)]
            after_inf = transitions.infection_step(current_types, virus_at_cell, g_hv, infect_rng)
            for cid in range(1, n_cells + 1):
                if after_inf[cid] == types.I and current_types[cid] != types.I:
                    world.set_cell_type(cid, types.I)
                    infected_ids.add(cid)
            current_types = after_inf

        # (3) per-cell resistance ρ from local IFN (every epithelial cell).
        resist_at_cell = [0.0] * (n_cells + 1)
        if "ifn" in enable:
            for cid in range(1, n_cells + 1):
                if current_types[cid] in _EPITHELIAL_TYPES:
                    f_bar = world.field_mean_at_cell(ifn_fi, cid)
                    resist_at_cell[cid] = cell_resistance(f_bar, a_rf)

        # (4) field-secretion scales. Virus release by infected gated (1-ρ);
        # macrophage chemokine/IL-10 via the DYNAMIC sig_1 with the §4a-
        # refreshed macrophage id list; uninfected IL-10 gated (1-ρ).
        for cid in infected_ids:
            world.set_cell_secretion_scale(virus_fi, cid, 1.0 - resist_at_cell[cid])
        if "chemokine" in enable:
            current_macrophage_ids = [cid for cid, ct in enumerate(world.cell_types())
                                      if ct == types.M and cid != 0]
            for cid in current_macrophage_ids:
                l_loc = world.field_mean_at_cell(il10_fi, cid)
                scale = signaling.macrophage_secretion_scale(l_loc, sig_1_dynamic, g_1, g_2, d_2)
                world.set_cell_secretion_scale(chemo_fi, cid, scale)
                world.set_cell_secretion_scale(il10_fi, cid, scale)
            for cid in range(1, n_cells + 1):
                if current_types[cid] == types.H:
                    world.set_cell_secretion_scale(
                        il10_fi, cid, signaling.uninfected_il10_scale(resist_at_cell[cid]))

        # (5) chemotaxis: persistent engine setting, applied inside world.step.

        # (6) contact + nearby killing I->D.
        if "killing" in enable and infected_ids:
            cell_volumes = world.cell_volumes()
            for cid in list(infected_ids):
                contact = world.cell_contact_area_by_type(cid)
                srf_nk = contact.get(types.K, 0)
                srf_cd8 = contact.get(types.E, 0)
                resist = resist_at_cell[cid]
                cvol = cell_volumes[cid]
                total_rate = (
                    killing.contact_kill_rate(srf_nk, resist, g_ik, tot_ec, cvol)
                    + killing.contact_kill_rate(srf_cd8, resist, g_ie, tot_ec, cvol)
                    + killing.nearby_kill_rate(g_ik, eta, state["K_nb"], resist)
                    + killing.nearby_kill_rate(g_ie, eta, state["E_nb"], resist)
                )
                if kill_rng.random() < 1.0 - math.exp(-total_rate):
                    world.set_cell_type(cid, types.D)
                    infected_ids.discard(cid)
                    dead_from_infected_ids.add(cid)

        # (7) Allee death H->D / recovery D->H (contact-geometry rule).
        if "allee" in enable:
            snapshot = list(world.cell_types())
            allee_writes = []
            for cid in range(1, n_cells + 1):
                t = snapshot[cid]
                if t not in (types.H, types.D):
                    continue
                srf_uninfected, srf_D, srf_total = _epithelial_contact_totals(world, cid)
                resist = resist_at_cell[cid]
                if t == types.H:
                    rate = allee.allee_death_rate(srf_D, srf_uninfected, srf_total,
                                                  b_h, resist, theta)
                    if allee_rng.random() < 1.0 - math.exp(-rate):
                        allee_writes.append((cid, types.D, "death"))
                else:  # types.D
                    rate = allee.allee_recovery_rate(srf_uninfected, srf_total,
                                                     b_h, resist, theta)
                    if allee_rng.random() < 1.0 - math.exp(-rate):
                        allee_writes.append((cid, types.H, "recover"))
            for cid, new_t, kind in allee_writes:
                world.set_cell_type(cid, new_t)
                if kind == "death":  # H->D: dead-from-HEALTHY (source DH)
                    dead_from_healthy_ids.add(cid)
                else:                # D->H: recovered, no longer dead
                    dead_from_healthy_ids.discard(cid)
                    dead_from_infected_ids.discard(cid)

        # (8) infected apoptosis I->D (mu_i·(1-ρ)).
        if "death" in enable and infected_ids:
            current_types = list(world.cell_types())
            after_apop = transitions.infected_death_step(current_types, resist_at_cell, mu_i, death_rng)
            for cid in range(1, n_cells + 1):
                if after_apop[cid] == types.D and current_types[cid] == types.I:
                    world.set_cell_type(cid, types.D)
                    infected_ids.discard(cid)
                    dead_from_infected_ids.add(cid)

    def _ode_couple(sig_1_dynamic):
        """Step 9: spatial->ODE push -> integrate -> ODE->spatial recruitment
        (source §4: __update_spatial_data -> rr.timestep() -> update_
        populations), reused wholesale from `run_global_coupling` incl. the
        §4a nearby-surrogate fix. Returns the fresh dynamic sig_1 for the next
        record's secretion. No-op (returns sig_1 unchanged) if ODE disabled."""
        if "ode" not in enable:
            return sig_1_dynamic
        types_now = world.cell_types()
        H, I, D, M, K, E = _counts(types_now)
        DH = sum(1 for cid in dead_from_healthy_ids if types_now[cid] == types.D)
        V = float(sum(world.field_conc(virus_fi))) / dim_z
        F = float(sum(world.field_conc(ifn_fi))) / dim_z
        C = float(sum(world.field_conc(chemo_fi))) / dim_z
        L = float(sum(world.field_conc(il10_fi))) / dim_z
        sum_resist = sum(resist_at_cell[cid] for cid in infected_ids)
        B_ei = sum_resist * b_ei
        G_ki = sum_resist * g_ki
        # §4a fix: num_immune_by_type = LOCAL CPM count + NEARBY surrogate.
        inputs = dict(H=H, I=I, M=M + state["M_nb"], K=K + state["K_nb"],
                      E=E + state["E_nb"], DH=DH, V=V, F=F, C=C, L=L,
                      B_ei=B_ei, G_ki=G_ki)
        # dt spans the mcs_per_step MCS advanced since the last coupling, so
        # ODE model-time stays synced to CPM time.
        new_state = ode.step(state, inputs, dt_seconds=s_per_mcs * mcs_per_step)
        state.clear()
        state.update(new_state)
        if "recruitment" in enable:
            _recruit_step(C, state["P"], G_ki, B_ei)
        # dynamic sig_1 for the next record's secretion (a_11*T + a_12*D).
        return a_11 * state["T"] + a_12 * (tot_cell - H - I)

    # sig_1 for the first secretion step: D = tot_cell - H - I = 0 at the
    # seeded initial state, T = healthy-IC T (0.0). Same causal ordering as
    # `run_global_coupling` (a record's secretion uses the PRECEDING record's
    # freshly-integrated T + spatial dead count; the initial record uses D=0).
    sig_1_dynamic = a_11 * state["T"] + a_12 * 0.0

    # Index 0 = the seeded initial state (before any transition), so
    # infected[0] reflects the scenario's initial seeding (frac -> >0, viral
    # load -> 0). Indices 1..steps-1 each advance `mcs_per_step` MCS.
    _record(0)

    for i in range(1, steps):
        for _ in range(mcs_per_step):
            _one_mcs(sig_1_dynamic)
        sig_1_dynamic = _ode_couple(sig_1_dynamic)
        _record(i)

    result["params"] = {
        "cells_per_side": cells_per_side, "tot_cell": tot_cell, "eta": eta,
        "steps": steps, "seed": seed, "s_per_mcs": s_per_mcs,
        "mcs_per_step": mcs_per_step,
        "init_infection_frac": init_infection_frac, "init_viral_load": init_viral_load,
        "n_infected_seeded": n_infected,
        "n_macrophages": n_macrophages, "n_nk": n_nk, "n_cd8": n_cd8,
        "enable": sorted(enable),
        "nk_cd8_chemotaxis_engine_scale": NK_CD8_CHEMOTAXIS_ENGINE_SCALE,
        "b_m": consts["b_m"], "b_k": consts["b_k"], "b_p": consts["b_p"],
        "n_infected_final": len(infected_ids),
        "dead_from_infected": len(dead_from_infected_ids),
        "dead_from_healthy": len(dead_from_healthy_ids),
        "recruit_pool_per_type": recruit_pool_per_type if "recruitment" in enable else 0,
        "dims": list(world.dims()),
    }
    return result
