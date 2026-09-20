import process_bigraph as pb

from pbg_cpm_studies.influenza.ode_process import SystemicODEProcess
from pbg_cpm_studies.influenza import price_ode
from pbg_cpm_studies.influenza.run import load_params


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
    out = proc.update({"H":1000,"I":200,"M":10,"K":5,"E":2,"DH":0,
                       "V":50.0,"F":0.0,"C":0.5,"L":0.0,"B_ei":0.0,"G_ki":0.0}, 1.0)
    assert set(out["ode_state"]).issuperset(set(price_ode.INTEGRATED_STATES))
    assert "macro_inflow" in out["recruit_drivers"]
    assert out["recruit_drivers"]["macro_inflow"] > 0.0
