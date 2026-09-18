"""Confluent epithelial sheet as a load_world spec for the Rust CPM engine.

A uniform grid of 5x5-site (10 um) cells of type H (uninfected), tiling a
square patch. No fields or immune cells yet (Increment 1). Adhesion for the
epithelial-only sheet keeps cells confluent; the immune J matrix (Table 3)
arrives in Increments 6-7.
"""
from __future__ import annotations

from .params import load_params
from . import types

CELL_SIDE_SITES = 5  # 5x5 = 25 sites = 10 um cell


def _sites_per_side(patch_mm: float, lattice_um: int) -> int:
    return int(round(patch_mm * 1000.0 / lattice_um))


def build_sheet_spec(patch_mm: float, seed: int = 17) -> dict:
    params = load_params()
    p = params["cpm"]
    adh = params["adhesion"]
    n = _sites_per_side(patch_mm, p["lattice_um"])          # sites per side
    cells_per_side = n // CELL_SIDE_SITES
    cells = []
    for gy in range(cells_per_side):
        for gx in range(cells_per_side):
            x0, y0 = gx * CELL_SIDE_SITES, gy * CELL_SIDE_SITES
            cells.append({
                "type": types.H,
                "target_volume": float(p["cell_sites"]),    # 25
                "lambda_volume": float(p["lambda_volume"]),  # 9
                "target_surface": 0.0, "lambda_surface": 0.0,
                "seed_block": [x0, y0, 0,
                               x0 + CELL_SIDE_SITES, y0 + CELL_SIDE_SITES, 1],
            })
    # Epithelial-only adhesion: cohesive H-H, higher H-medium so the sheet
    # stays confluent. Source-literal values (ViralInfectionVTM.xml), read
    # from params.yaml (single source of truth) rather than duplicated here.
    contact = [
        {"a": types.MEDIUM, "b": types.MEDIUM, "j": 0.0},
        {"a": types.MEDIUM, "b": types.H, "j": float(adh["epithelial_medium"])},   # CONTROLLER OVERRIDE: source-literal epithelial-medium (was 16.0 CC3D-typical)
        {"a": types.H, "b": types.H, "j": float(adh["epithelial_epithelial"])},    # CONTROLLER OVERRIDE: source-literal epithelial-epithelial (was 4.0 CC3D-typical)
    ]
    return {
        "potts": {"dims": [n, n, 1], "boundary": "noflux",
                  "neighbor_order": int(p["neighbor_order"]),
                  "temperature": float(p["temperature"]), "seed": seed},
        "cells": cells,
        "contact": contact,
    }
