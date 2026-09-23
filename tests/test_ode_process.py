import process_bigraph as pb

from viva_cpm_studies.influenza.ode_process import SystemicODEProcess
from viva_cpm_studies.influenza import price_ode
from viva_cpm_studies.influenza.run import load_params


def test_ode_process_advances_and_emits_drivers():
    # NOTE: unlike the brief's illustrative snippet, `Process.__init__`
    # requires a `core` (see `test_immune_process.py`) -- `pb.allocate_core()`
    # supplied here, matching that convention.
    p = load_params()
    consts = price_ode.resolve_constants(p["price_ode"], num_epithelial=1225)
    proc = SystemicODEProcess(
        {"consts": consts, "num_epithelial": 1225, "dt_seconds": 420.0},
        core=pb.allocate_core(),
    )
    apc_input = {"H":1000,"I":200,"M":10,"K":5,"E":2,"DH":0,
                 "V":50.0,"F":0.0,"C":0.5,"L":0.0,"B_ei":0.0,"G_ki":0.0}

    # Test-rigor fix (fix round 1): the two assertions below (ode_state keys,
    # macro_inflow > 0) are satisfiable even by a no-op step -- initial_state()
    # already has every INTEGRATED_STATES key, and macro_inflow's baseline
    # mu_m*b_m term is >0 regardless of stepping. Capture P before/after two
    # successive updates under an APC-producing input to prove the ODE
    # actually integrates AND that self.state persists across update() calls
    # (not re-initialized each call).
    p0 = proc.state["P"]
    out1 = proc.update(apc_input, 1.0)
    p1 = out1["ode_state"]["P"]
    assert p1 != p0
    assert p1 > p0, "P should rise under an APC-producing input"
    assert proc.state["P"] == p1, "self.state must persist the stepped value"

    out = proc.update(apc_input, 1.0)
    p2 = out["ode_state"]["P"]
    assert p2 != p1, "second update must continue integrating from the persisted state, not reset it"
    assert p2 > p1

    assert set(out["ode_state"]).issuperset(set(price_ode.INTEGRATED_STATES))
    assert "macro_inflow" in out["recruit_drivers"]
    assert out["recruit_drivers"]["macro_inflow"] > 0.0
