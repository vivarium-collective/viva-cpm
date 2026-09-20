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


def test_nk_agent_kills_adjacent_infected():
    # NOTE: like the tests above, `core=pb.allocate_core()` is required
    # (deviation from the brief's illustrative snippet, same as Task 2.1).
    p = ImmuneProcess({"seed": 2, "kill_radius": 5.0}, core=pb.allocate_core())
    agents = [{"id": 1, "type": int(types.K), "x": 10.0, "y": 5.0}]
    killed_any = False
    for _ in range(200):
        out = p.update({"chemo_field": [0.0] * 300, "virus_field": [0.0] * 300,
                        "dims": [30, 10, 1],
                        "immune_agents": agents,
                        "epithelial_positions": {"7": [11.0, 5.0]}, "infected_ids": [7],
                        "sig_1": 0.0}, 1.0)
        if 7 in out["immune_kills"]:
            killed_any = True
            break
    assert killed_any


def test_macrophage_deposits_chemokine_when_sig1_positive():
    p = ImmuneProcess({"seed": 3}, core=pb.allocate_core())
    agents = [{"id": 1, "type": int(types.M), "x": 10.0, "y": 5.0}]
    out = p.update({"chemo_field": [0.0] * 300, "virus_field": [0.0] * 300,
                    "dims": [30, 10, 1],
                    "immune_agents": agents, "epithelial_positions": {}, "infected_ids": [],
                    "sig_1": 0.02}, 1.0)
    deposits = [d for d in out["field_deposit"] if d[0] == 2]  # chemokine field idx
    assert deposits and deposits[0][4] > 0.0
