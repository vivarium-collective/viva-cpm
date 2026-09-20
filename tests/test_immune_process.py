import numpy as np
import process_bigraph as pb
from pbg_cpm_studies.influenza.immune_process import ImmuneProcess
from pbg_cpm_studies.influenza import types


def _linear_field(nx, ny):
    # increasing in +x
    return [float(x) for _ in range(ny) for x in range(nx)]  # site=x+y*nx


def test_nk_agent_moves_up_chemokine_gradient():
    # NOTE: unlike the brief's illustrative snippet, `Process.__init__`
    # requires a `core` (see `epithelium_process`'s tests) -- `pb.allocate_core()`
    # supplied here, matching that convention.
    nx, ny = 30, 10
    p = ImmuneProcess({"seed": 1}, core=pb.allocate_core())
    agents = [{"id": 1, "type": int(types.K), "x": 10.0, "y": 5.0}]
    out = p.update({"chemo_field": _linear_field(nx, ny),
                    "virus_field": [0.0] * (nx * ny),
                    "dims": [nx, ny, 1], "immune_agents": agents}, 1.0)
    assert out["immune_agents"][0]["x"] > 10.0   # moved toward higher chemokine


def test_macrophage_agent_moves_up_virus_gradient():
    # Field-per-type routing: macrophages (types.M) chemotax up the VIRUS
    # field, not the chemokine field -- prove routing by giving the
    # macrophage a flat chemokine field (no signal there) and a +x-increasing
    # virus field; it should still move +x.
    nx, ny = 30, 10
    p = ImmuneProcess({"seed": 1}, core=pb.allocate_core())
    agents = [{"id": 1, "type": int(types.M), "x": 10.0, "y": 5.0}]
    out = p.update({"chemo_field": [0.0] * (nx * ny),
                    "virus_field": _linear_field(nx, ny),
                    "dims": [nx, ny, 1], "immune_agents": agents}, 1.0)
    assert out["immune_agents"][0]["x"] > 10.0   # moved toward higher virus concentration
