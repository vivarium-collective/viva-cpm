"""Task 5.1: macrophage in the world + chemotaxis up the virus field.

Increment 5 crux: introduce the macrophage cell type (`types.M`) into a CPM
world and show it chemotaxes up the extracellular-virus gradient toward the
infection (paper Fig 2B / Fig 3A -- macrophages accumulate at the lesion).
Uses existing engine primitives only (`World.add_cell`, `World.set_contact`,
`World.set_chemotaxis` on the virus field already wired by `fields.
add_virus_field`) -- no Rust change.

## Geometry: documented 2D approximation

The Increment-1..4 base sheet (`sheet.build_sheet_spec`) is a CONFLUENT z=1
patch -- every lattice site is covered by an epithelial cell, so there is no
Medium space for a macrophage to occupy or move through. The source model
instead places immune cells in a z=2 layer *above* the epithelial sheet
(a true 3rd dimension the source's CC3D lattice has and viva-cpm's z=1
increments so far do not use for the epithelium).

For this MECHANISM increment we approximate that geometry in 2D: a
NON-CONFLUENT epithelial patch (a block of H cells, with some seeded I cells
so the virus field has a source from MCS 0) sits inside a much larger domain
with open Medium around it, and macrophage (`types.M`) cells are placed in
that Medium, away from the patch. This gives macrophages literal lattice
space to migrate through in-plane while chemotaxing up the virus gradient --
the same shape as `pbg_cpm_studies/chemotaxis`'s secreting-source +
chemotaxing-responder scenario, with the VIRUS field as the attractant and
`types.I` epithelial cells as the (incidental) source.

The faithful z=2 immune layer over a CONFLUENT z=1 epithelial sheet is
deferred to Increment 9 (full reproduction increment) -- this increment only
needs to show the chemotaxis mechanism engages and localizes macrophages,
not reproduce the source's exact lattice topology.
"""
from __future__ import annotations

from .params import load_params
from .sheet import CELL_SIDE_SITES
from . import types


def _epithelial_block_cells(cells_per_side, x_off, y_off, target_volume, lambda_volume):
    cells = []
    for gy in range(cells_per_side):
        for gx in range(cells_per_side):
            x0, y0 = x_off + gx * CELL_SIDE_SITES, y_off + gy * CELL_SIDE_SITES
            cells.append({
                "type": types.H,
                "target_volume": float(target_volume),
                "lambda_volume": float(lambda_volume),
                "target_surface": 0.0, "lambda_surface": 0.0,
                "seed_block": [x0, y0, 0, x0 + CELL_SIDE_SITES, y0 + CELL_SIDE_SITES, 1],
            })
    return cells


def build_macrophage_scenario_spec(*, epithelial_cells_per_side: int = 4,
                                    n_infected: int = 1, n_macrophages: int = 6,
                                    margin_sites: int = 30, seed: int = 17) -> dict:
    """Non-confluent macrophage-response scenario (see module docstring for the
    2D-approximation rationale): a small epithelial patch (mostly H, with
    ``n_infected`` cells nearest its center seeded as I) sits centered in a
    much larger open-Medium domain; ``n_macrophages`` type-M cells are placed
    in the Medium, in a far corner of the domain -- away from the infection --
    so `run.run_macrophage_response` measures real migration, not a trivial
    already-there start.

    Adhesion (Table 3 / ``params.yaml`` ``macrophage.adhesion``) and volume
    (``params.yaml`` ``macrophage.volume_sites`` / ``lambda_volume``) are set
    from the cited macrophage params, matching every other epithelial-type row
    already used by `sheet.build_sheet_spec`.
    """
    params = load_params()
    p = params["cpm"]
    adh = params["adhesion"]
    macro = params["macrophage"]
    macro_adh = macro["adhesion"]

    epithelial_side_sites = epithelial_cells_per_side * CELL_SIDE_SITES
    n = epithelial_side_sites + 2 * margin_sites  # square domain, patch centered

    x_off = y_off = margin_sites
    cells = _epithelial_block_cells(epithelial_cells_per_side, x_off, y_off,
                                     p["cell_sites"], p["lambda_volume"])

    # Seed the `n_infected` cells nearest the patch's geometric center as
    # INFECTED (type I) so the virus field (wired via `fields.add_virus_field`,
    # which secretes from every `types.I` cell) has a source from MCS 0 --
    # same convention as `run.run_virus_infection`'s pre-finalize seeding.
    mid = (epithelial_cells_per_side - 1) / 2.0
    order = sorted(range(len(cells)), key=lambda i: (
        (i % epithelial_cells_per_side - mid) ** 2
        + (i // epithelial_cells_per_side - mid) ** 2))
    for i in order[:n_infected]:
        cells[i]["type"] = types.I

    # Macrophages: placed in the Medium along the domain's near (low-x/low-y)
    # edge -- as far from the centered epithelial patch/infection as the
    # domain allows -- spaced out so their seed blocks don't overlap.
    mv = int(macro["volume_sites"])
    mlv = float(macro["lambda_volume"])
    macro_side = int(round(mv ** 0.5))  # 5 for 25 sites
    x0, y0, gap = 2, 2, macro_side + 2
    per_row = max(1, (n - 4) // gap)
    for k in range(n_macrophages):
        row, col = divmod(k, per_row)
        cx = x0 + col * gap
        cy = y0 + row * gap
        cells.append({
            "type": types.M,
            "target_volume": float(mv),
            "lambda_volume": mlv,
            "target_surface": 0.0, "lambda_surface": 0.0,
            "seed_block": [cx, cy, 0, cx + macro_side, cy + macro_side, 1],
        })

    contact = [
        {"a": types.MEDIUM, "b": types.MEDIUM, "j": 0.0},
        {"a": types.MEDIUM, "b": types.H, "j": float(adh["epithelial_medium"])},
        {"a": types.H, "b": types.H, "j": float(adh["epithelial_epithelial"])},
        {"a": types.MEDIUM, "b": types.I, "j": float(adh["epithelial_medium"])},
        {"a": types.H, "b": types.I, "j": float(adh["epithelial_epithelial"])},
        {"a": types.I, "b": types.I, "j": float(adh["epithelial_epithelial"])},
        {"a": types.MEDIUM, "b": types.M, "j": float(macro_adh["medium_macrophage"])},
        {"a": types.H, "b": types.M, "j": float(macro_adh["uninfected_macrophage"])},
        {"a": types.I, "b": types.M, "j": float(macro_adh["infected_macrophage"])},
        {"a": types.M, "b": types.M, "j": float(macro_adh["macrophage_macrophage"])},
    ]

    return {
        "potts": {"dims": [n, n, 1], "boundary": "noflux",
                  "neighbor_order": int(p["neighbor_order"]),
                  "temperature": float(p["temperature"]), "seed": seed},
        "cells": cells,
        "contact": contact,
    }


def set_macrophage_chemotaxis(world, virus_fi, *, chemotaxis_v_macro: float | None = None) -> None:
    """Wire macrophage (`types.M`) chemotaxis up the virus field.

    ``chemotaxis_v_macro`` defaults to ``params.yaml``'s ``macrophage.
    chemotaxis_v_macro`` (5000, source `ImmuneModelInputs.py`'s
    `chemotaxis_v_macro = 5E3`). Pass ``0.0`` to build the lambda=0 control
    (chemotaxis mechanically disabled -- see `run.run_macrophage_response`).

    Note (functional-form gap, see ``params.yaml``'s ``chemotaxis.
    functional_form``): the source uses a saturating, per-cell-COM-scaled
    lambda (``chemotaxis_v_macro / (1 + concentration)``); the viva-cpm engine
    primitive (`World.set_chemotaxis`) is a single fixed (LINEAR-regime)
    lambda per (field, type). This increment uses the engine's linear form
    -- fidelity criterion is the localization BEHAVIOR, not the exact
    functional form (documented gap, deferred).
    """
    params = load_params()
    lam = (float(params["macrophage"]["chemotaxis_v_macro"])
           if chemotaxis_v_macro is None else float(chemotaxis_v_macro))
    world.set_chemotaxis(virus_fi, types.M, lam)
