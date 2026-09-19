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

from . import allee, build, fields, immune, killing, sheet, signaling, transitions, types
from .params import load_params
from .resistance import cell_resistance

# Epithelial cell types that participate in the Allee contact-surface
# calculation (matches the source's `ec_list = [UNINFECTED, INFECTED,
# INFECTEDRELEASING, DYING]`, i.e. H/I/D here) -- MEDIUM and any future
# immune-cell contact must be excluded from `srf_total`/`srf_uninfected`,
# not just tallied by `world.cell_contact_area_by_type`'s raw dict.
_EPITHELIAL_TYPES = (types.H, types.I, types.D)


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
