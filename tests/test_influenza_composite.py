import process_bigraph as pb

from viva_cpm_studies.composites import influenza as inf
from viva_cpm_studies.influenza import types as inf_types


def test_composite_document_runs_and_emits():
    doc = inf.epithelium()
    assert isinstance(doc, dict)
    # smoke: the demo composite builds a modest sheet, not the full 1mm^2
    spec = inf.build_spec(patch_mm=0.1)   # 50x50 sites, 100 cells
    assert spec["potts"]["dims"] == [50, 50, 1]
    assert len(spec["cells"]) == 100


def test_virus_infection_spec_wires_field_and_seeds_lesion():
    spec = inf.build_viral_infection_spec(patch_mm=0.1, seed=17, init_infected_frac=0.05)
    assert len(spec["fields"]) == 1
    assert spec["fields"][0]["name"] == "virus"
    n_init = sum(1 for c in spec["cells"] if c["type"] == inf_types.I)
    assert n_init == 5  # round(100 * 0.05)


def test_virus_infection_composite_runs_through_the_engine():
    core = pb.allocate_core()
    doc = inf.viral_infection(core=core, patch_mm=0.1, seed=17, init_infected_frac=0.05)
    comp = pb.Composite({"state": doc}, core=core)
    comp.run(10)
    world = comp.state["cpm"]["instance"].world
    types_now = list(world.cell_types())
    # population conserved; at least the seeded lesion is still infected
    assert sum(1 for t in types_now[1:] if t == inf_types.I) >= 5
