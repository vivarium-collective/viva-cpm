import process_bigraph as pb
from pbg_cpm_studies.composites import chemotaxis_adaptive as CA
from pbg_cpm_studies.chemotaxis import metrics as M


def _final_index(condition, seed=17, steps=40):
    core = pb.allocate_core()
    cfg = CA.CONDITIONS[condition]
    doc = CA.recruitment_adaptive(core=core, seed=seed, **cfg)
    comp = pb.Composite({"state": doc}, core=core)
    comp.run(steps)
    world = comp.state["cpm"]["instance"].world
    return M.recruitment_index(world, responder_type=CA.ACTIVATED_TYPE)


def test_conditions_exist():
    assert set(CA.CONDITIONS) == {"low_bg", "mid_bg", "high_bg", "high_bg_blocked", "high_bg_knockout"}


def test_low_bg_recruits_and_blocked_does_not():
    assert _final_index("low_bg") > _final_index("high_bg_blocked")
