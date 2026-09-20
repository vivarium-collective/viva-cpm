import numpy as np
import process_bigraph as pb
import pytest
from pbg_cpm_studies.influenza.immune_process import ImmuneProcess
from pbg_cpm_studies.influenza import types
from pbg_cpm_studies.influenza.fields import il10_hill_constants
from pbg_cpm_studies.influenza.killing import contact_kill_rate
from pbg_cpm_studies.influenza.params import load_params
from pbg_cpm_studies.influenza.signaling import macrophage_secretion_scale


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
    #
    # Fix round 1: a non-stochastic assertion that the underlying rate is
    # actually > 0 for this scenario (decouples "rate correctly wired" from
    # "did a random draw land"), PLUS a strengthened stochastic loop (5000
    # draws, not 200) so the draw-based assertion is robust to an arbitrary
    # seed / RNG-consumption-order change, not just the specific seed=2 used
    # here. Neither the kill rate nor `kill_radius` were touched (both are
    # source-faithful per the Task 2.2 report) -- only the test's robustness.
    params = load_params()
    g_ik = float(params["nk"]["g_ik"])
    tot_ec = float(params["scaling"]["ode_epithelial_population"])
    rate = contact_kill_rate(srf_immune=1.0, cell_resist=1.0, g_i=g_ik,
                              tot_ec=tot_ec, cell_volume=1.0)
    assert rate > 0.0

    p = ImmuneProcess({"seed": 2, "kill_radius": 5.0}, core=pb.allocate_core())
    agents = [{"id": 1, "type": int(types.K), "x": 10.0, "y": 5.0}]
    killed_any = False
    for _ in range(5000):
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
    sig_1 = 0.02
    p = ImmuneProcess({"seed": 3}, core=pb.allocate_core())
    agents = [{"id": 1, "type": int(types.M), "x": 10.0, "y": 5.0}]
    out = p.update({"chemo_field": [0.0] * 300, "virus_field": [0.0] * 300,
                    "dims": [30, 10, 1],
                    "immune_agents": agents, "epithelial_positions": {}, "infected_ids": [],
                    "sig_1": sig_1}, 1.0)
    deposits = [d for d in out["field_deposit"] if d[0] == 2]  # chemokine field idx
    assert deposits and deposits[0][4] > 0.0

    # Fix round 1 (controller ruling R4): pin the WHOLE-CELL secretion factor
    # -- a macrophage agent represents one whole CPM cell, which secretes its
    # per-site rate at every one of its ~cell_sites occupied pixels on-lattice,
    # so the deposited amount must equal b_c_per_site * cell_sites * scale.
    params = load_params()
    cell_sites = int(params["cpm"]["cell_sites"])
    b_c_per_site = (float(params["chemokine"]["b_c"]) / 2.0) / cell_sites
    _, g_1, g_2, d_2 = il10_hill_constants()
    scale = macrophage_secretion_scale(0.0, sig_1, g_1, g_2, d_2)
    expected_chemo = b_c_per_site * cell_sites * scale
    assert deposits[0][4] == pytest.approx(expected_chemo)
