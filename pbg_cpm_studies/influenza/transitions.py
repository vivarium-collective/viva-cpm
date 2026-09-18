"""Stochastic virus-driven infection transition H -> I (Increment 2, Task 2.2).

Pure function: given the current per-cell types and the local virus field
seen by each cell, returns a NEW list of types with some `types.H` cells
transitioned to `types.I`, chosen stochastically. The input list is never
mutated.
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
