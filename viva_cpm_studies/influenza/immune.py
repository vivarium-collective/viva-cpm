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
so the virus field has a source from MCS 0) and a cluster of macrophage
(`types.M`) cells both sit in the INTERIOR of a much larger domain with open
Medium around and between them -- neither cluster starts pinned against a
noflux wall (see `build_macrophage_scenario_spec`'s docstring for why: a
wall-pinned start biases the lambda=0 control's drift). This gives
macrophages literal lattice space to migrate through in-plane while
chemotaxing up the virus gradient -- the same shape as
`viva_cpm_studies/chemotaxis`'s secreting-source + chemotaxing-responder
scenario, with the VIRUS field as the attractant and `types.I` epithelial
cells as the (incidental) source.

The faithful z=2 immune layer over a CONFLUENT z=1 epithelial sheet is
deferred to Increment 9 (full reproduction increment) -- this increment only
needs to show the chemotaxis mechanism engages and localizes macrophages,
not reproduce the source's exact lattice topology.
"""
from __future__ import annotations

import math

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
                                    margin_sites: int = 30, separation_sites: int = 25,
                                    seed: int = 17) -> dict:
    """Non-confluent macrophage-response scenario (see module docstring for the
    2D-approximation rationale): a small epithelial patch (mostly H, with
    ``n_infected`` cells nearest its center seeded as I) and a cluster of
    ``n_macrophages`` type-M cells both sit in the domain's INTERIOR, each at
    least ``margin_sites`` from every noflux wall, separated horizontally by
    ``separation_sites`` of open Medium.

    Fix (Task 5.1 review finding #1 -- CONTROL-DRIFT CONFOUND): an earlier
    version placed the macrophages in a domain CORNER (pinned against two
    walls at once). A noflux wall is a hard reflector, so a corner-pinned
    cluster's thermal/wall drift is NOT an isotropic random walk -- it is
    biased to net-drift away from that corner, and because the corner was
    chosen diagonally opposite the (centered) infection, that bias happened
    to point roughly TOWARD the infection, inflating the apparent
    chemotaxis-off "control" drift and hence overstating the on/off gap
    attributable to chemotaxis. Placing BOTH clusters deep in the interior
    (``margin_sites`` default large relative to a run's random-walk
    displacement scale) removes that wall-proximity bias entirely: with no
    wall nearby, the lambda=0 control's drift has no preferred direction with
    respect to the infection, so it is a fair (isotropic) baseline for the
    chemotaxis-on comparison. `run.run_macrophage_response` validates this
    empirically across seeds (see its module docstring / task-5.1-report.md's
    fix note) rather than asserting it analytically here.

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

    mv = int(macro["volume_sites"])
    mlv = float(macro["lambda_volume"])
    macro_side = int(round(mv ** 0.5))  # 5 for 25 sites
    macro_gap = 2
    macro_cols = max(1, int(math.ceil(math.sqrt(n_macrophages))))
    macro_rows = max(1, int(math.ceil(n_macrophages / macro_cols)))
    macro_cluster_w = macro_cols * macro_side + (macro_cols - 1) * macro_gap
    macro_cluster_h = macro_rows * macro_side + (macro_rows - 1) * macro_gap

    # Domain: both clusters interior, each >= margin_sites from every wall.
    # Horizontally: margin | epithelial patch | separation | macrophage
    # cluster | margin. Vertically: both clusters are centered on the same
    # row, so each gets AT LEAST margin_sites of clearance above/below (more,
    # for whichever cluster is shorter than the tallest one).
    ny = margin_sites * 2 + max(epithelial_side_sites, macro_cluster_h)
    nx = (margin_sites * 2 + epithelial_side_sites + separation_sites
          + macro_cluster_w)

    x_off = margin_sites
    y_off = margin_sites + (max(epithelial_side_sites, macro_cluster_h) - epithelial_side_sites) // 2
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

    # Macrophages: a compact interior cluster, `separation_sites` of open
    # Medium to the right of the epithelial patch -- away from the infection,
    # but nowhere near a wall (see fix note above).
    mac_x_off = margin_sites + epithelial_side_sites + separation_sites
    mac_y_off = margin_sites + (max(epithelial_side_sites, macro_cluster_h) - macro_cluster_h) // 2
    for k in range(n_macrophages):
        row, col = divmod(k, macro_cols)
        cx = mac_x_off + col * (macro_side + macro_gap)
        cy = mac_y_off + row * (macro_side + macro_gap)
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
        "potts": {"dims": [nx, ny, 1], "boundary": "noflux",
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


def _cluster_layout(n_cells: int, side: int, gap: int) -> tuple[int, int, int, int]:
    """``(cols, rows, width, height)`` for a compact grid cluster of
    ``n_cells`` ``side`` x ``side``-site cells, ``gap`` sites of open Medium
    between adjacent cells within the cluster -- the same packing convention
    `build_macrophage_scenario_spec` uses for its macrophage cluster,
    extracted here so `build_cytotoxic_scenario_spec`'s NK/CD8 clusters below
    use the identical packing rather than a re-derived one."""
    cols = max(1, int(math.ceil(math.sqrt(n_cells))))
    rows = max(1, int(math.ceil(n_cells / cols)))
    width = cols * side + (cols - 1) * gap
    height = rows * side + (rows - 1) * gap
    return cols, rows, width, height


def _cluster_cells(cell_type: int, target_volume: float, lambda_volume: float,
                    n_cells: int, side: int, gap: int, x_off: int, y_off: int) -> list[dict]:
    """Cell specs for a compact grid cluster of ``n_cells`` cells of
    ``cell_type``, packed per `_cluster_layout` and seeded starting at
    ``(x_off, y_off)``."""
    cols, _rows, _w, _h = _cluster_layout(n_cells, side, gap)
    cells = []
    for k in range(n_cells):
        row, col = divmod(k, cols)
        cx = x_off + col * (side + gap)
        cy = y_off + row * (side + gap)
        cells.append({
            "type": cell_type,
            "target_volume": float(target_volume),
            "lambda_volume": float(lambda_volume),
            "target_surface": 0.0, "lambda_surface": 0.0,
            "seed_block": [cx, cy, 0, cx + side, cy + side, 1],
        })
    return cells


def build_cytotoxic_scenario_spec(*, epithelial_cells_per_side: int = 4,
                                   n_infected: int = 1, n_macrophages: int = 6,
                                   n_nk: int = 6, n_cd8: int = 6,
                                   margin_sites: int = 30, separation_sites: int = 25,
                                   seed: int = 17) -> dict:
    """Task 7.1 scenario: extend `build_macrophage_scenario_spec`'s
    non-confluent epithelial/infection + macrophage layout with a THIRD
    interior cluster holding NK (`types.K`) + CD8+ (`types.E`) cells,
    ``separation_sites`` further out beyond the macrophage cluster (away
    from the infection). Domain layout (every cluster INTERIOR, >=
    ``margin_sites`` from every noflux wall -- same fix-round-1 rigor as
    `build_macrophage_scenario_spec`, applied to the new cluster too: an
    unbiased lambda=0 control needs every cluster's start free of
    wall-proximity drift bias, see that function's docstring for the
    corner-control-drift finding this avoids):

        margin | epithelial/infection patch | separation | macrophage
        cluster | separation | NK+CD8 cluster | margin

    The macrophage cluster sits BETWEEN the infection and the NK/CD8
    cluster: it is both the chemokine SOURCE (Increment 6's
    `run_macrophage_signaling`) the NK/CD8 cells chemotax up, and, as it
    itself chemotaxes up the virus field toward the infection (Increment 5,
    unchanged by this task), it stays on the infection side of the NK/CD8
    cluster throughout a run -- so the chemokine gradient direction seen
    from the NK/CD8 cluster continues to point toward the infection even as
    the macrophage cluster drifts.

    NK and CD8 cells are seeded as two side-by-side sub-clusters (NK above
    CD8, separated by the same ``gap`` used within each cluster) inside the
    combined "far" zone, so each type's own localization can be measured
    independently.

    Adhesion rows for K/E vs every existing type, and K<->E, come from the
    literal Contact J values in `params.yaml`'s `nk.adhesion`/`cd8.adhesion`
    (source `ViralInfectionVTM.xml`), matching `build_macrophage_scenario_spec`'s
    convention for the macrophage rows. Volume/lambda_volume for K/E come
    from `params.yaml`'s `nk.volume_sites`/`nk.lambda_volume` and
    `cd8.volume_sites`/`cd8.lambda_volume` (both identical to macrophage's:
    25 sites, lambda_volume 9 -- see those keys' params.yaml comments).
    """
    params = load_params()
    p = params["cpm"]
    adh = params["adhesion"]
    macro = params["macrophage"]
    macro_adh = macro["adhesion"]
    nk = params["nk"]
    nk_adh = nk["adhesion"]
    cd8 = params["cd8"]
    cd8_adh = cd8["adhesion"]

    epithelial_side_sites = epithelial_cells_per_side * CELL_SIDE_SITES
    gap = 2

    mv, mlv = int(macro["volume_sites"]), float(macro["lambda_volume"])
    macro_side = int(round(mv ** 0.5))
    _mc, _mr, macro_cluster_w, macro_cluster_h = _cluster_layout(n_macrophages, macro_side, gap)

    nkv, nklv = int(nk["volume_sites"]), float(nk["lambda_volume"])
    nk_side = int(round(nkv ** 0.5))
    _nc, _nr, nk_w, nk_h = _cluster_layout(n_nk, nk_side, gap)

    cv, clv = int(cd8["volume_sites"]), float(cd8["lambda_volume"])
    cd8_side = int(round(cv ** 0.5))
    _cc, _cr, cd8_w, cd8_h = _cluster_layout(n_cd8, cd8_side, gap)

    far_cluster_w = max(nk_w, cd8_w)
    far_cluster_h = nk_h + gap + cd8_h

    ny = margin_sites * 2 + max(epithelial_side_sites, macro_cluster_h, far_cluster_h)
    nx = (margin_sites * 2 + epithelial_side_sites + separation_sites
          + macro_cluster_w + separation_sites + far_cluster_w)

    x_off = margin_sites
    y_off = margin_sites + (ny - margin_sites * 2 - epithelial_side_sites) // 2
    cells = _epithelial_block_cells(epithelial_cells_per_side, x_off, y_off,
                                     p["cell_sites"], p["lambda_volume"])

    # Seed the `n_infected` cells nearest the patch's geometric center as
    # INFECTED (type I) -- same convention as `build_macrophage_scenario_spec`.
    mid = (epithelial_cells_per_side - 1) / 2.0
    order = sorted(range(len(cells)), key=lambda i: (
        (i % epithelial_cells_per_side - mid) ** 2
        + (i // epithelial_cells_per_side - mid) ** 2))
    for i in order[:n_infected]:
        cells[i]["type"] = types.I

    mac_x_off = margin_sites + epithelial_side_sites + separation_sites
    mac_y_off = margin_sites + (ny - margin_sites * 2 - macro_cluster_h) // 2
    cells += _cluster_cells(types.M, mv, mlv, n_macrophages, macro_side, gap,
                             mac_x_off, mac_y_off)

    far_x_off = mac_x_off + macro_cluster_w + separation_sites
    far_y_center = margin_sites + (ny - margin_sites * 2) // 2
    nk_y_off = far_y_center - far_cluster_h // 2
    cd8_y_off = nk_y_off + nk_h + gap
    cells += _cluster_cells(types.K, nkv, nklv, n_nk, nk_side, gap, far_x_off, nk_y_off)
    cells += _cluster_cells(types.E, cv, clv, n_cd8, cd8_side, gap, far_x_off, cd8_y_off)

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
        {"a": types.MEDIUM, "b": types.K, "j": float(nk_adh["medium_nkcell"])},
        {"a": types.H, "b": types.K, "j": float(nk_adh["uninfected_nkcell"])},
        {"a": types.I, "b": types.K, "j": float(nk_adh["infected_nkcell"])},
        {"a": types.M, "b": types.K, "j": float(nk_adh["macrophage_nkcell"])},
        {"a": types.K, "b": types.K, "j": float(nk_adh["nkcell_nkcell"])},
        {"a": types.MEDIUM, "b": types.E, "j": float(cd8_adh["medium_cd8tcell"])},
        {"a": types.H, "b": types.E, "j": float(cd8_adh["uninfected_cd8tcell"])},
        {"a": types.I, "b": types.E, "j": float(cd8_adh["infected_cd8tcell"])},
        {"a": types.M, "b": types.E, "j": float(cd8_adh["macrophage_cd8tcell"])},
        {"a": types.E, "b": types.E, "j": float(cd8_adh["cd8tcell_cd8tcell"])},
        {"a": types.K, "b": types.E, "j": float(nk_adh["nkcell_cd8tcell"])},
    ]

    return {
        "potts": {"dims": [nx, ny, 1], "boundary": "noflux",
                  "neighbor_order": int(p["neighbor_order"]),
                  "temperature": float(p["temperature"]), "seed": seed},
        "cells": cells,
        "contact": contact,
    }


def set_nk_cd8_chemotaxis(world, chemo_fi, *, chemotaxis_v_nk: float | None = None,
                           chemotaxis_v_cd8: float | None = None) -> None:
    """Wire NK (`types.K`) + CD8+ (`types.E`) chemotaxis up the CHEMOKINE
    field (``chemo_fi``, from `fields.add_chemokine_field` -- macrophage-
    secreted, Increment 6), NOT the virus field (the macrophages' own
    chemotaxis target, see `set_macrophage_chemotaxis`) -- confirmed by
    `params.yaml`'s `nk.chemotaxis_field`/`cd8.chemotaxis_field` (both
    ``"chemokine"``).

    ``chemotaxis_v_nk``/``chemotaxis_v_cd8`` default to `params.yaml`'s
    `nk.chemotaxis_v_nk` (5000) / `cd8.chemotaxis_v_cd8` (10000) -- CD8+'s
    lambda is 2x NK's (paper Sec. 2.3: CD8+ sensitivity "twice that of NK
    cells"). Pass ``0.0`` for the lambda=0 control (either or both).

    Same linear-lambda engine-primitive gap as `set_macrophage_chemotaxis`
    (source uses a saturating, per-cell-COM-scaled lambda; `World.
    set_chemotaxis` here is a single fixed (linear-regime) lambda per
    (field, type)) -- this increment's fidelity criterion is the
    localization BEHAVIOR, not the exact functional form (documented gap,
    deferred).
    """
    params = load_params()
    lam_nk = (float(params["nk"]["chemotaxis_v_nk"])
              if chemotaxis_v_nk is None else float(chemotaxis_v_nk))
    lam_cd8 = (float(params["cd8"]["chemotaxis_v_cd8"])
               if chemotaxis_v_cd8 is None else float(chemotaxis_v_cd8))
    world.set_chemotaxis(chemo_fi, types.K, lam_nk)
    world.set_chemotaxis(chemo_fi, types.E, lam_cd8)
