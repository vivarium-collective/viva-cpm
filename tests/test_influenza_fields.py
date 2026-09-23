"""Task 2.1: extracellular virus field builder.

Verifies fields.add_virus_field(world) sets up secretion (types.I only) and
diffusion so that infecting one cell produces a local, growing, decaying-
with-distance virus field, while an all-H sheet (no secretion) stays at ~0.
"""
from __future__ import annotations

import math

from viva_cpm_studies.influenza import build, fields, sheet, types


def _dims_and_conc(world, field_idx):
    conc = world.field_conc(field_idx)
    n = int(round(math.sqrt(len(conc))))  # square sheet (dims n x n x 1)
    return n, conc


def _annulus_mean(conc, n, cx, cy, r_lo, r_hi):
    total, count = 0.0, 0
    for y in range(n):
        for x in range(n):
            r = math.hypot(x - cx, y - cy)
            if r_lo <= r < r_hi:
                total += conc[x + y * n]
                count += 1
    assert count > 0
    return total / count


def test_infected_cell_secretes_local_virus_that_decays_with_distance():
    spec = sheet.build_sheet_spec(0.3)
    # World.add_field is only valid pre-finalize (engine constraint), so build
    # unfinalized, wire the field, then finalize ourselves.
    world = build.world_from_spec(spec, finalize=False)

    # Interior cell, well away from the domain boundary (150x150 sites, 30x30
    # cell grid): pick grid cell (15, 15) -> cell id 1 + 15*30 + 15 = 466.
    cells_per_side = 150 // sheet.CELL_SIDE_SITES
    gx, gy = 15, 15
    cell_id = 1 + gy * cells_per_side + gx
    cx = gx * sheet.CELL_SIDE_SITES + sheet.CELL_SIDE_SITES / 2.0
    cy = gy * sheet.CELL_SIDE_SITES + sheet.CELL_SIDE_SITES / 2.0

    field_idx = fields.add_virus_field(world)
    world.finalize(spec["potts"]["seed"])
    world.set_cell_type(cell_id, types.I)

    world.step(30)

    n, conc = _dims_and_conc(world, field_idx)
    assert max(conc) > 0.0, "secreting cell should raise the virus field above 0"

    near = _annulus_mean(conc, n, cx, cy, 0.0, 8.0)
    far = _annulus_mean(conc, n, cx, cy, 40.0, 45.0)
    assert near > 0.0
    assert far < near, f"virus should decay with distance: near={near}, far={far}"


def test_all_uninfected_sheet_stays_virus_free():
    spec = sheet.build_sheet_spec(0.3)
    world = build.world_from_spec(spec, finalize=False)
    field_idx = fields.add_virus_field(world)
    world.finalize(spec["potts"]["seed"])

    world.step(30)

    conc = world.field_conc(field_idx)
    assert sum(conc) < 1e-6, "no I cells -> no secretion -> field should stay ~0"


def test_virus_diffusion_length_sites_constant():
    assert fields.VIRUS_DIFFUSION_LENGTH_SITES == 25


def test_infected_cell_secretes_local_ifn_that_decays_with_distance():
    spec = sheet.build_sheet_spec(0.3)
    world = build.world_from_spec(spec, finalize=False)

    cells_per_side = 150 // sheet.CELL_SIDE_SITES
    gx, gy = 15, 15
    cell_id = 1 + gy * cells_per_side + gx
    cx = gx * sheet.CELL_SIDE_SITES + sheet.CELL_SIDE_SITES / 2.0
    cy = gy * sheet.CELL_SIDE_SITES + sheet.CELL_SIDE_SITES / 2.0

    field_idx = fields.add_ifn_field(world)
    world.finalize(spec["potts"]["seed"])
    world.set_cell_type(cell_id, types.I)

    world.step(30)

    n, conc = _dims_and_conc(world, field_idx)
    assert max(conc) > 0.0, "secreting cell should raise the IFN field above 0"

    near = _annulus_mean(conc, n, cx, cy, 0.0, 8.0)
    far = _annulus_mean(conc, n, cx, cy, 40.0, 45.0)
    assert near > 0.0
    assert far < near, f"IFN should decay with distance: near={near}, far={far}"


def test_all_uninfected_sheet_stays_ifn_free():
    spec = sheet.build_sheet_spec(0.3)
    world = build.world_from_spec(spec, finalize=False)
    field_idx = fields.add_ifn_field(world)
    world.finalize(spec["potts"]["seed"])

    world.step(30)

    conc = world.field_conc(field_idx)
    assert sum(conc) < 1e-6, "no I cells -> no secretion -> field should stay ~0"


def test_ifn_diffusion_length_sites_constant():
    assert fields.IFN_DIFFUSION_LENGTH_SITES == 10


def test_ifn_spatial_extent_shorter_than_virus():
    # Sanity check (optional per brief): IFN's characteristic diffusion length
    # is shorter than the virus field's (10 sites vs 25).
    assert fields.IFN_DIFFUSION_LENGTH_SITES < fields.VIRUS_DIFFUSION_LENGTH_SITES
