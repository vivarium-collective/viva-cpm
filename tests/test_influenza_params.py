from pbg_cpm_studies.influenza import params, types


def test_type_codes_are_distinct_and_ordered():
    codes = [types.MEDIUM, types.H, types.I, types.D, types.M, types.K, types.E]
    assert codes == [0, 1, 2, 3, 4, 5, 6]
    assert len(set(codes)) == 7


def test_params_load_core_cpm_values():
    p = params.load_params()
    cpm = p["cpm"]
    assert cpm["lattice_um"] == 2
    assert cpm["cell_sites"] == 25            # 5x5
    assert cpm["volume_constraint_um2"] == 100
    assert cpm["lambda_volume"] == 9
    assert cpm["temperature"] == 10
    assert cpm["neighbor_order"] == 3   # CONTROLLER OVERRIDE: source-literal (was 2)
    assert cpm["dt_min"] == 1


def test_params_carry_provenance_for_every_top_section():
    p = params.load_params()
    # every parameter section must cite a source string (paper table/eq or CC3D source)
    for section in ("cpm", "adhesion", "diffusion", "chemotaxis", "scaling"):
        assert p[section].get("source"), f"{section} missing provenance"


def test_diffusion_coefficients_match_table3():
    d = params.load_params()["diffusion"]
    assert d["virus_um2_s"] == 0.0119
    assert d["chemokine_um2_s"] == 1.04
    assert d["ifn_type1_um2_s"] == 0.520
    assert d["il10_um2_s"] == 0.327


def test_params_carry_virus_section():
    v = params.load_params()["virus"]
    for k in ("diffusion_lat2_per_mcs","decay_per_mcs","diffusion_length_cell_diam",
              "secretion_g_vi","infection_g_hv","source"):
        assert k in v, f"virus params missing {k}"
    assert v["diffusion_length_cell_diam"] == 5
    assert isinstance(v["secretion_g_vi"], (int,float)) and v["secretion_g_vi"] > 0
    assert isinstance(v["infection_g_hv"], (int,float)) and v["infection_g_hv"] > 0
