import process_bigraph as pb

from viva_cpm.processes.cpm_process import CPMProcess


def _spec(j11=2.0, ncells=1):
    cells = [{"type": 1, "target_volume": 25, "lambda_volume": 2.0,
              "target_surface": 20, "lambda_surface": 0.0,
              "seed_block": [10, 10, 0, 15, 15, 1]}]
    return {"potts": {"dims": [40, 40, 1], "boundary": "periodic",
                      "neighbor_order": 2, "temperature": 10.0, "seed": 7},
            "cell_types": [{"name": "medium", "id": 0}, {"name": "a", "id": 1}],
            "contact": [{"a": 0, "b": 1, "j": 16.0}, {"a": 1, "b": 1, "j": j11}],
            "cells": cells}


def _run(config, updates):
    proc = CPMProcess(config, core=pb.allocate_core())
    out = {}
    for _ in range(updates):
        out = proc.update({}, 1.0)
    return proc, out


def test_mitosis_grows_cell_count():
    mit = [{"type": 1, "grow_rate": 2.0, "divide_threshold": 50.0, "reset_target": 25.0}]
    proc, _ = _run({"spec": _spec(), "mcs_per_update": 5, "mitosis": mit}, 20)
    assert proc.world.n_cells() > 1


def test_no_mitosis_config_is_noop():
    proc, out = _run({"spec": _spec(), "mcs_per_update": 5}, 5)
    assert proc.world.n_cells() == 1
    assert "n_propagules" not in out


def test_aggregate_output_present_and_nonnull():
    _, out = _run({"spec": _spec(), "mcs_per_update": 5,
                   "report_aggregate_types": [1]}, 5)
    assert out["n_propagules"] is not None
    assert out["n_propagules"] >= 1
    assert 0.0 < out["largest_cluster_fraction"] <= 1.0


def test_strong_adhesion_stays_cohesive_weak_fractures_under_growth():
    mit = [{"type": 1, "grow_rate": 2.0, "divide_threshold": 50.0, "reset_target": 25.0}]
    def go(j):
        return _run({"spec": _spec(j11=j), "mcs_per_update": 5, "mitosis": mit,
                     "report_aggregate_types": [1]}, 40)[1]
    strong, weak = go(2.0), go(14.0)
    assert strong["n_propagules"] <= weak["n_propagules"]
