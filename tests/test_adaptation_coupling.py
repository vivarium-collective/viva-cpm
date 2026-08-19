from cpm.coupling import adaptation_coupling, ADAPTIVE_ADDR

def test_fragment_shape_and_preinit():
    cfg = {"kd": 2.9, "hill": 2.0, "conc_scale": 0.02, "activate_occupancy": 0.5,
           "epsilon": 0.1, "background": 0.0, "naive_type": 2, "activated_type": 3}
    frag = adaptation_coupling([2, 3], receptor_config=cfg)
    node = frag["receptor_2"]
    assert node["address"] == ADAPTIVE_ADDR
    assert node["inputs"]["ligand"] == ["field_at_cell", "2"]
    assert node["inputs"]["m_prev"] == ["adaptation", "2"]
    assert node["outputs"]["fate"] == ["fates", "2"]
    assert node["outputs"]["m"] == ["adaptation", "2"]
    # BOTH per-cell maps pre-initialized (map-write-to-absent-key is dropped)
    assert frag["fates"] == {"2": 2, "3": 2}
    assert frag["adaptation"] == {"2": 0.0, "3": 0.0}
