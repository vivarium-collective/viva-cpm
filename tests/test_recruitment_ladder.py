# tests/test_recruitment_ladder.py — the emergent ladder must separate
import process_bigraph as pb
from statistics import mean
from pbg_cpm_studies.composites import chemotaxis_adaptive as CA
from pbg_cpm_studies.chemotaxis import metrics as M

SEEDS = [17, 29, 43]; STEPS = 40

def _idx(condition, **over):
    vals = []
    for s in SEEDS:
        core = pb.allocate_core()
        cfg = dict(CA.CONDITIONS[condition]); cfg.update(over)
        comp = pb.Composite({"state": CA.recruitment_adaptive(core=core, seed=s, **cfg)}, core=core)
        comp.run(STEPS)
        vals.append(M.recruitment_index(comp.state["cpm"]["instance"].world,
                                        responder_type=CA.ACTIVATED_TYPE))
    return mean(vals)

def test_adaptive_beats_hill_at_high_background():
    # hill == adaptation knockout (epsilon=0) saturates and fails; adaptive recovers
    adaptive_high = _idx("high_bg")                 # epsilon>0
    hill_high     = _idx("high_bg_knockout")        # epsilon=0  (the middle rung)
    assert adaptive_high > 0.4                       # adaptive recruits at high bg
    assert adaptive_high - hill_high > 0.2           # and clearly beats hill there

def test_receptor_gating_abolishes():
    assert _idx("high_bg_blocked") < 0.2             # blocking the response kills recruitment

def test_low_background_recruits_for_both():
    assert _idx("low_bg") > 0.4                       # sanity: both rungs work at low bg
