from viva_cpm.schema import load_world
from viva_cpm.metrics import aggregate_components, largest_cluster_fraction


def _cell(type_, block):
    return {"type": type_, "target_volume": 50, "lambda_volume": 2.0,
            "target_surface": 40, "lambda_surface": 0.0, "seed_block": block}


def _world(cells):
    spec = {"potts": {"dims": [40, 20, 1], "boundary": "noflux", "neighbor_order": 2,
                      "temperature": 10.0, "seed": 1},
            "cells": cells}
    return load_world(spec)


def _world_with_two_blobs(type_=1, sizes=(70, 30)):
    # Rectangles 10 rows tall; width = size // 10, so size == site count.
    # Three empty columns separate them (detached).
    a, b = sizes
    wa, wb = a // 10, b // 10
    return _world([
        _cell(type_, [0, 2, 0, wa, 12, 1]),
        _cell(type_, [wa + 3, 2, 0, wa + 3 + wb, 12, 1]),
    ])


def _world_with_two_separated_blobs(type_=1):
    return _world_with_two_blobs(type_, (50, 50))


def _world_with_single_blob(type_=1):
    return _world([_cell(type_, [5, 5, 0, 12, 12, 1])])


def test_two_detached_blobs_same_type_counts_two():
    world = _world_with_two_separated_blobs(type_=1)
    assert aggregate_components(world, [1]) == 2


def test_one_cohesive_blob_counts_one():
    world = _world_with_single_blob(type_=1)
    assert aggregate_components(world, [1]) == 1


def test_empty_mask_is_zero_not_crash():
    world = _world_with_single_blob(type_=1)
    assert aggregate_components(world, [9]) == 0
    assert largest_cluster_fraction(world, [9]) == 0.0


def test_largest_fraction_on_70_30_split():
    world = _world_with_two_blobs(type_=1, sizes=(70, 30))
    assert abs(largest_cluster_fraction(world, [1]) - 0.7) < 0.05


def test_adjacent_same_type_cells_merge_into_one_component():
    # Two distinct cells touching along x=8 -> one aggregate, not a per-cell sum.
    world = _world([_cell(1, [2, 2, 0, 8, 12, 1]), _cell(1, [8, 2, 0, 14, 12, 1])])
    assert world.n_cells() == 2
    assert aggregate_components(world, [1]) == 1
    assert largest_cluster_fraction(world, [1]) == 1.0


def test_adjacent_cells_of_different_types_merge_only_when_both_selected():
    world = _world([_cell(1, [2, 2, 0, 8, 12, 1]), _cell(2, [8, 2, 0, 14, 12, 1])])
    assert aggregate_components(world, [1, 2]) == 1
    assert aggregate_components(world, [1]) == 1
