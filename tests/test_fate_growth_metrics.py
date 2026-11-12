from viva_cpm.schema import load_world
from viva_cpm.metrics import fate_counts, viable_biomass


def _cell(type_, block):
    return {"type": type_, "target_volume": 50, "lambda_volume": 2.0,
            "target_surface": 40, "lambda_surface": 0.0, "seed_block": block}


def _world():
    # Rectangles with 10 rows; volumes 70, 30, 20 (types 1, 1, 2), detached.
    spec = {"potts": {"dims": [40, 20, 1], "boundary": "noflux", "neighbor_order": 2,
                      "temperature": 10.0, "seed": 1},
            "cells": [_cell(1, [0, 2, 0, 7, 12, 1]),
                      _cell(1, [10, 2, 0, 13, 12, 1]),
                      _cell(2, [20, 2, 0, 22, 12, 1])]}
    return load_world(spec)


def test_fate_counts_tallies_by_type():
    assert fate_counts(_world()) == {1: 2, 2: 1}


def test_viable_biomass_sums_above_floor():
    w = _world()
    vols = w.cell_volumes()[1:]
    assert sorted(vols) == [20, 30, 70]
    assert viable_biomass(w, min_volume=25) == 100.0
    assert viable_biomass(w, min_volume=1.0) == 120.0
    assert viable_biomass(w, min_volume=1000) == 0.0
    assert isinstance(viable_biomass(w), float)
