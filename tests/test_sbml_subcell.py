import numpy as np
import process_bigraph as pb
from viva_cpm.subcellular.sbml import SBMLSubcell

# S is produced at a Wnt-gated rate and decays; Wnt is a held floating species
MODEL = """
J1: -> S; k_on*Wnt^4/(K^4 + Wnt^4);
J2: S -> ; k_off*S;
species S, Wnt;
S = 1.0; Wnt = 0.0;
k_on = 0.8; K = 0.3; k_off = 0.4;
"""

CFG = {"model": MODEL, "ligand_species": "Wnt", "state_species": "S", "ligand_scale": 1.0}


def _proc():
    return SBMLSubcell(CFG, core=pb.allocate_core())


def test_low_ligand_lets_stemness_decay():
    p = _proc()
    s_last = 1.0
    for _ in range(20):
        s_last = p.update({"ligand": 0.0}, 1.0)["state"]
    assert s_last < 0.2            # no Wnt -> S decays toward 0


def test_high_ligand_sustains_stemness():
    p = _proc()
    s_last = 1.0
    for _ in range(20):
        s_last = p.update({"ligand": 1.0}, 1.0)["state"]
    assert s_last > 0.8            # saturating Wnt -> S held high


def test_update_reads_state_from_array_species(monkeypatch):
    """viva-tellurium >= 0.1.2 bulk mode returns `species` as a numpy ARRAY ordered by
    getFloatingSpeciesIds, not a {id: value} map. SBMLSubcell must read the chosen state
    species by its POSITION — indexing an array with the string id raises IndexError (the
    CI break with libroadrunner 2.11 / viva-tellurium 0.1.2, green only on the older map-
    returning 0.1.0)."""
    p = _proc()
    p.update({"ligand": 0.0}, 1.0)                 # one real step to populate species ids
    ids = list(p._tp._species_ids)
    expected = {sid: float(i + 1) for i, sid in enumerate(ids)}   # distinct per species
    arr = np.array([expected[sid] for sid in ids], dtype=float)
    monkeypatch.setattr(p._tp, "update", lambda inp, interval: {"species": arr})
    out = p.update({"ligand": 0.0}, 1.0)
    assert out["state"] == expected[p.state_species]
