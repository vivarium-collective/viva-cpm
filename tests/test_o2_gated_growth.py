import numpy as np
import process_bigraph as pb

from viva_cpm import metrics
from viva_cpm.processes.cpm_process import CPMProcess

DIMS = (40, 20, 1)


def _spec():
    def cell(x0):
        return {"type": 1, "target_volume": 25, "lambda_volume": 2.0,
                "target_surface": 20, "lambda_surface": 0.0,
                "seed_block": [x0, 8, 0, x0 + 5, 13, 1]}
    return {"potts": {"dims": list(DIMS), "boundary": "noflux",
                      "neighbor_order": 2, "temperature": 10.0, "seed": 7},
            "cell_types": [{"name": "medium", "id": 0}, {"name": "a", "id": 1},
                           {"name": "b", "id": 2}],
            "contact": [{"a": 0, "b": 1, "j": 16.0}, {"a": 0, "b": 2, "j": 16.0},
                        {"a": 1, "b": 1, "j": 2.0}, {"a": 2, "b": 2, "j": 2.0},
                        {"a": 1, "b": 2, "j": 2.0}],
            # type 1 secretes O2 locally (no diffusion); type 2 does not
            "fields": [{"name": "o2", "d": 0.0, "decay": 0.0,
                        "secretion": [{"type": 1, "rate": 1.0}]}],
            # cell 1 (type 1) is oxygenated, cell 2 (type 2) is hypoxic
            "cells": [cell(5), dict(cell(30), type=2)]}


def _proc(mitosis=(), **extra):
    cfg = {"spec": _spec(), "mcs_per_update": 5, "n_fields": 1,
           "mitosis": list(mitosis), **extra}
    p = CPMProcess(cfg, core=pb.allocate_core())
    return p


GATE = {"field": 0, "threshold": 0.5}
GATED = [{"type": t, "grow_rate": 15.0, "o2_gate": GATE} for t in (1, 2)]


def test_gated_growth_only_oxygenated():
    p = _proc(GATED)
    v0 = list(p.world.cell_volumes())
    for _ in range(6):
        p.update({}, 1.0)
    assert p.world.field_mean_at_cell(0, 1) >= 0.5 > p.world.field_mean_at_cell(0, 2)
    v1 = list(p.world.cell_volumes())
    assert v1[1] - v0[1] > 20          # oxygenated cell grew
    assert v1[2] - v0[2] <= 3          # hypoxic cell did not grow (may jitter down)
    assert (v1[1] - v0[1]) - (v1[2] - v0[2]) > 20


def test_gated_entry_no_division():
    p = _proc(GATED)
    n0 = len(p.world.cell_types())
    for _ in range(6):
        p.update({}, 1.0)
    assert len(p.world.cell_types()) == n0


def test_fate_growth_outputs_present():
    p = _proc(report_fate_growth=True)
    assert "fate_counts" in p.outputs() and "viable_biomass" in p.outputs()
    out = p.update({}, 1.0)
    assert out["fate_counts"] == {str(t): c for t, c in metrics.fate_counts(p.world).items()}
    assert all(isinstance(k, str) and isinstance(v, int) for k, v in out["fate_counts"].items())
    assert isinstance(out["viable_biomass"], float)
    assert out["viable_biomass"] == metrics.viable_biomass(p.world)
    off = _proc()
    assert "fate_counts" not in off.outputs()
    assert "fate_counts" not in off.update({}, 1.0)


def test_backward_compat_uniform_grow_divide():
    mit = [{"type": 1, "grow_rate": 10.0, "divide_threshold": 40.0, "reset_target": 25.0}]
    p = _proc(mit)
    n0 = p.world.n_cells()
    for _ in range(15):
        p.update({}, 1.0)
    assert p.world.n_cells() > n0
