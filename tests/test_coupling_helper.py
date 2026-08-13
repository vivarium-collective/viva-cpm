from cpm.coupling import receptor_coupling

RCFG = {"kd": 10.0, "hill": 1.0, "conc_scale": 1.0, "activate_occupancy": 0.5,
        "naive_type": 2, "activated_type": 3}


def test_builds_one_subcell_per_cell_with_wiring():
    frag = receptor_coupling([4, 7], receptor_config=RCFG)
    assert frag["receptor_4"]["inputs"]["ligand"] == ["field_at_cell", "4"]
    assert frag["receptor_4"]["outputs"]["fate"] == ["fates", "4"]
    assert frag["receptor_7"]["inputs"]["ligand"] == ["field_at_cell", "7"]
    assert frag["receptor_4"]["address"].endswith("ReceptorSubcell")
    assert frag["receptor_4"]["config"]["kd"] == 10.0


def test_fates_store_preinitialized_for_every_cell():
    frag = receptor_coupling([4, 7], receptor_config=RCFG)
    assert frag["fates"] == {"4": 2, "7": 2}   # H4 pre-init fallback (naive_type)
