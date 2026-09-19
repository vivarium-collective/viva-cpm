"""Stochastic per-cell type transitions for the influenza-sego2022 CPM.

Each transition is a pure function: given the current per-cell types and a
relevant per-cell field, returns a NEW list of types with some cells
stochastically transitioned. The input list is never mutated.

- `infection_step` (Increment 2, Task 2.2): H -> I, driven by local virus.
- `infected_death_step` (Increment 4, Task 4.2): I -> D, driven by
  mu_i*(1-resist).
"""
from __future__ import annotations

import math

from . import types


def infection_step(types_list, virus_at_cell, g_hv, rng) -> list:
    """One MCS of the H -> I infection transition.

    `types_list`: sequence of int cell types, index = cell id. Index 0 is
    the medium sentinel (see tests/test_mitosis.py convention) and is always
    left unchanged.
    `virus_at_cell`: sequence/mapping of per-cell mean local virus
    concentration, same indexing as `types_list` (e.g. `world.field_mean_at_cell`
    per cell id).
    `g_hv`: per-MCS infection rate coefficient (params.yaml `virus.infection_g_hv`).
    `rng`: a `numpy.random.Generator` (e.g. `numpy.random.default_rng(seed)`).

    For each `types.H` cell: `rate = g_hv * v_bar`, `Pr(infect) = 1 - exp(-rate)`;
    the cell becomes `types.I` iff `rng.random() < Pr`. All other cells
    (already I/D/medium/immune types, and the index-0 sentinel) are copied
    through unchanged. Returns a new list; `types_list` is not mutated.
    """
    new_types = list(types_list)
    for cell_id in range(1, len(types_list)):
        if types_list[cell_id] != types.H:
            continue
        v_bar = virus_at_cell[cell_id]
        rate = g_hv * v_bar
        pr_infect = 1.0 - math.exp(-rate)
        if rng.random() < pr_infect:
            new_types[cell_id] = types.I
    return new_types


def infected_death_step(types_list, resist_at_cell, mu_i, rng) -> list:
    """One MCS of the I -> D infected-cell death transition (Increment 4, Task 4.2).

    `types_list`: sequence of int cell types, index = cell id. Index 0 is
    the medium sentinel and is always left unchanged.
    `resist_at_cell`: sequence/mapping of per-cell viral resistance (0..1,
    see resistance.cell_resistance), same indexing as `types_list`.
    `mu_i`: per-MCS base death rate coefficient (params.yaml
    `cell_death.mu_i_per_mcs`).
    `rng`: a `numpy.random.Generator` (e.g. `numpy.random.default_rng(seed)`).

    This is the CC3D source's SIMPLIFIED flat form (ViralCellDeathSteppable's
    own docstring: "simpified [sic] version") -- rate = mu_i * (1 - resist),
    with NO viral-load/ROS Hill saturation term. The paper's fuller printed
    Table 2 "Infected death" row is the conceptual target this omits (see
    transcription doc's discrepancy note); params.yaml mirrors what the
    source actually runs, not the paper's printed form.

    For each `types.I` cell: `rate = mu_i * (1 - resist)`,
    `Pr(die) = 1 - exp(-rate)`; the cell becomes `types.D` iff
    `rng.random() < Pr`. A fully resistant cell (resist == 1) has rate 0 and
    never dies. All other cells (H/D/medium/immune types, and the index-0
    sentinel) are copied through unchanged. Returns a new list; `types_list`
    is not mutated.
    """
    new_types = list(types_list)
    for cell_id in range(1, len(types_list)):
        if types_list[cell_id] != types.I:
            continue
        resist = resist_at_cell[cell_id]
        rate = mu_i * (1.0 - resist)
        pr_die = 1.0 - math.exp(-rate)
        if rng.random() < pr_die:
            new_types[cell_id] = types.D
    return new_types
