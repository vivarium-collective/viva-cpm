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


def test_macrophage_secretion_self_limits_with_local_il10():
    """Source-faithfulness fix (final review round): `run_full_model` samples
    the IL-10 field at each macrophage's own cell (`world.field_mean_at_cell
    (il10_fi, cid)`) and feeds it into `macrophage_secretion_scale`, which is
    monotone DECREASING in local IL-10 (self-limiting secretion). Previously
    `ImmuneProcess` hardcoded `il10_local = 0.0` (the Hill term's maximum,
    unself-limited value) -- this proves the fix: a macrophage sitting in a
    NON-ZERO local IL-10 field must secrete LESS chemokine than the exact
    same agent in an all-zero IL-10 field, all else (sig_1, position, RNG
    seed) held identical."""
    nx, ny = 30, 10
    sig_1 = 0.02
    agents = [{"id": 1, "type": int(types.M), "x": 10.0, "y": 5.0}]
    base_kwargs = {
        "chemo_field": [0.0] * (nx * ny), "virus_field": [0.0] * (nx * ny),
        "dims": [nx, ny, 1], "immune_agents": agents,
        "epithelial_positions": {}, "infected_ids": [], "sig_1": sig_1,
    }

    p_zero = ImmuneProcess({"seed": 3}, core=pb.allocate_core())
    out_zero = p_zero.update({**base_kwargs, "il10_field": [0.0] * (nx * ny)}, 1.0)

    p_nonzero = ImmuneProcess({"seed": 3}, core=pb.allocate_core())
    il10_field = [5.0] * (nx * ny)  # uniform, well above the Hill's d_2 scale
    out_nonzero = p_nonzero.update({**base_kwargs, "il10_field": il10_field}, 1.0)

    chemo_zero = [d for d in out_zero["field_deposit"] if d[0] == 2]
    chemo_nonzero = [d for d in out_nonzero["field_deposit"] if d[0] == 2]
    assert chemo_zero and chemo_nonzero
    assert chemo_nonzero[0][4] < chemo_zero[0][4], (
        f"expected non-zero local IL-10 to REDUCE macrophage chemokine "
        f"secretion (self-limiting feedback): zero-IL10={chemo_zero[0][4]!r} "
        f"nonzero-IL10={chemo_nonzero[0][4]!r}")

    # Cross-check against the shared Hill-scale helper directly, at the exact
    # sampled value (agent sits at integer site (10, 5) -> il10_local = 5.0).
    _, g_1, g_2, d_2 = il10_hill_constants()
    scale_zero = macrophage_secretion_scale(0.0, sig_1, g_1, g_2, d_2)
    scale_nonzero = macrophage_secretion_scale(5.0, sig_1, g_1, g_2, d_2)
    assert scale_nonzero < scale_zero
    params = load_params()
    cell_sites = int(params["cpm"]["cell_sites"])
    b_c_per_site = (float(params["chemokine"]["b_c"]) / 2.0) / cell_sites
    assert chemo_zero[0][4] == pytest.approx(b_c_per_site * cell_sites * scale_zero)
    assert chemo_nonzero[0][4] == pytest.approx(b_c_per_site * cell_sites * scale_nonzero)


def test_macrophage_secretion_defaults_to_zero_il10_when_field_absent():
    """Backward-compat: a caller that doesn't wire `il10_field` at all (e.g.
    an older driver/test) must fall back to the pre-fix `il10_local = 0.0`
    behavior, not raise."""
    nx, ny = 30, 10
    sig_1 = 0.02
    p = ImmuneProcess({"seed": 3}, core=pb.allocate_core())
    agents = [{"id": 1, "type": int(types.M), "x": 10.0, "y": 5.0}]
    out = p.update({"chemo_field": [0.0] * (nx * ny), "virus_field": [0.0] * (nx * ny),
                    "dims": [nx, ny, 1], "immune_agents": agents,
                    "epithelial_positions": {}, "infected_ids": [], "sig_1": sig_1}, 1.0)
    deposits = [d for d in out["field_deposit"] if d[0] == 2]
    _, g_1, g_2, d_2 = il10_hill_constants()
    scale = macrophage_secretion_scale(0.0, sig_1, g_1, g_2, d_2)
    params = load_params()
    cell_sites = int(params["cpm"]["cell_sites"])
    b_c_per_site = (float(params["chemokine"]["b_c"]) / 2.0) / cell_sites
    assert deposits[0][4] == pytest.approx(b_c_per_site * cell_sites * scale)


def test_recruitment_grows_agents_uncapped():
    # NOTE: like the tests above, `core=pb.allocate_core()` is required
    # (deviation from the brief's illustrative snippet, same as Task 2.1/2.2).
    #
    # A large macrophage inflow driver, applied every update for 50 updates,
    # must grow the macrophage agent count well beyond any fixed reserve-pool
    # cap (`run_full_model`'s pool is a handful of cells per type) -- proving
    # recruitment here mints agents UNCAPPED, the point of Task 2.3.
    p = ImmuneProcess({"seed": 4}, core=pb.allocate_core())
    agents = [{"id": 1, "type": int(types.M), "x": 5.0, "y": 5.0}]
    drivers = {"macro_inflow": 5.0, "nk_inflow": 0.0, "cd8_inflow": 0.0,
               "macro_outflow": 0.0, "nk_outflow": 0.0, "cd8_outflow": 0.0}
    for _ in range(50):
        out = p.update({"chemo_field": [0.0] * 400, "virus_field": [0.0] * 400,
                        "dims": [20, 20, 1],
                        "immune_agents": agents, "epithelial_positions": {}, "infected_ids": [],
                        "sig_1": 0.0, "recruit_drivers": drivers, "margin_box": [1, 1, 19, 19],
                        "apply_recruitment": True}, 1.0)
        agents = out["immune_agents"]
    n_macro = sum(1 for a in agents if a["type"] == int(types.M))
    assert n_macro > 100   # far beyond the old pool_per_type cap


def test_recruitment_noop_when_apply_recruitment_false():
    # Recruitment must be gated on `apply_recruitment`: a driver that would
    # otherwise clearly grow the population leaves the agent count unchanged
    # when the flag is falsy/absent (Composite-controlled cadence gate).
    p = ImmuneProcess({"seed": 5}, core=pb.allocate_core())
    agents = [{"id": 1, "type": int(types.M), "x": 5.0, "y": 5.0}]
    drivers = {"macro_inflow": 5.0, "nk_inflow": 0.0, "cd8_inflow": 0.0,
               "macro_outflow": 0.0, "nk_outflow": 0.0, "cd8_outflow": 0.0}
    out = p.update({"chemo_field": [0.0] * 400, "virus_field": [0.0] * 400,
                    "dims": [20, 20, 1],
                    "immune_agents": agents, "epithelial_positions": {}, "infected_ids": [],
                    "sig_1": 0.0, "recruit_drivers": drivers, "margin_box": [1, 1, 19, 19],
                    "apply_recruitment": False}, 1.0)
    assert len(out["immune_agents"]) == 1


def test_immune_counts_match_agent_list():
    # Task 3.2 (controller ruling R3 addendum): `immune_counts` must report
    # M/K/E counts by type from the CURRENT agent list -- the ODE's M/K/E
    # inputs. No recruitment/movement drivers here, just a static agent mix.
    p = ImmuneProcess({"seed": 7}, core=pb.allocate_core())
    agents = (
        [{"id": i, "type": int(types.M), "x": 5.0, "y": 5.0} for i in range(1, 3)]
        + [{"id": i, "type": int(types.K), "x": 5.0, "y": 5.0} for i in range(3, 7)]
        + [{"id": i, "type": int(types.E), "x": 5.0, "y": 5.0} for i in range(7, 8)]
    )
    out = p.update({"chemo_field": [0.0] * 400, "virus_field": [0.0] * 400,
                    "dims": [20, 20, 1],
                    "immune_agents": agents, "epithelial_positions": {}, "infected_ids": [],
                    "sig_1": 0.0}, 1.0)
    assert out["immune_counts"] == {"M": 2, "K": 4, "E": 1}


def test_recruitment_outflow_removes_agents():
    # A large outflow rate on an existing population of macrophages, with no
    # inflow, must shrink the agent count (proximity/killing/secretion don't
    # otherwise change agent count).
    p = ImmuneProcess({"seed": 6}, core=pb.allocate_core())
    agents = [{"id": i, "type": int(types.M), "x": 5.0, "y": 5.0} for i in range(1, 21)]
    drivers = {"macro_inflow": 0.0, "nk_inflow": 0.0, "cd8_inflow": 0.0,
               "macro_outflow": 5.0, "nk_outflow": 0.0, "cd8_outflow": 0.0}
    out = p.update({"chemo_field": [0.0] * 400, "virus_field": [0.0] * 400,
                    "dims": [20, 20, 1],
                    "immune_agents": agents, "epithelial_positions": {}, "infected_ids": [],
                    "sig_1": 0.0, "recruit_drivers": drivers, "margin_box": [1, 1, 19, 19],
                    "apply_recruitment": True}, 1.0)
    assert len(out["immune_agents"]) < 20
