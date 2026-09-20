from pbg_cpm_studies.composites.influenza import full_model_composite_document
from pbg_cpm_studies.core import build_core
from process_bigraph import Composite


def test_full_model_composite_builds_and_steps():
    doc = full_model_composite_document(cells_per_side=6, seed=1, init_infection_frac=0.1)
    assert set(doc).issuperset({"epithelium", "immune", "ode"})
    core = build_core()
    comp = Composite({"state": doc}, core=core)
    comp.run(1)   # one interval, no exception
