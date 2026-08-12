import process_bigraph as pb
from cpm.subcellular.receptor import ReceptorSubcell

CFG = {"kd": 10.0, "hill": 1.0, "conc_scale": 1.0, "activate_occupancy": 0.5,
       "naive_type": 2, "activated_type": 3}


def _proc(**over):
    cfg = {**CFG, **over}
    return ReceptorSubcell(cfg, core=pb.allocate_core())


def test_occupancy_half_at_kd():
    assert abs(_proc(kd=10.0).occupancy(10.0) - 0.5) < 1e-9


def test_occupancy_monotonic_and_bounded():
    p = _proc(kd=10.0)
    vals = [p.occupancy(c) for c in (0.0, 1.0, 5.0, 10.0, 50.0, 500.0)]
    assert vals[0] == 0.0
    assert all(0.0 <= v <= 1.0 for v in vals)
    assert all(b >= a for a, b in zip(vals, vals[1:]))


def test_fate_activates_above_threshold():
    p = _proc(kd=10.0, activate_occupancy=0.5)
    assert p.update({"ligand": 20.0}, 1.0)["fate"] == 3   # theta ~0.667 -> activated
    assert p.update({"ligand": 2.0}, 1.0)["fate"] == 2    # theta ~0.167 -> naive
    assert p.update({"ligand": 0.0}, 1.0)["fate"] == 2


def test_deterministic():
    p = _proc()
    assert p.update({"ligand": 12.3}, 1.0) == p.update({"ligand": 12.3}, 1.0)
