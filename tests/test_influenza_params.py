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


def test_params_carry_ifn_and_resistance_sections():
    p = params.load_params()
    ifn = p["ifn"]; res = p["resistance"]
    for k in ("diffusion_lat2_per_mcs","decay_per_mcs","diffusion_length_cell_diam","secretion_g_fp","source"):
        assert k in ifn, f"ifn missing {k}"
    assert ifn["diffusion_length_cell_diam"] == 2
    assert ifn["secretion_g_fp"] > 0
    assert "a_rf" in res and res["a_rf"] > 0 and "source" in res


def test_params_carry_death_and_allee_sections():
    p = params.load_params()
    cd = p["cell_death"]; al = p["allee"]
    assert cd.get("mu_i_per_mcs", 0) > 0 and "source" in cd
    for k in ("b_h", "srf_threshold", "source"):
        assert k in al, f"allee missing {k}"
    assert al["b_h"] > 0


def test_params_carry_macrophage_section():
    m = params.load_params()["macrophage"]
    assert m.get("chemotaxis_v_macro", 0) == 5000
    for k in ("adhesion", "source"):
        assert k in m, f"macrophage missing {k}"
    assert m["adhesion"]["macrophage_macrophage"] == 25.0
    assert m["adhesion"]["infected_macrophage"] == 20.0  # source-literal, NOT paper's collapsed 10
    assert m["volume_sites"] == 25
    assert m["lambda_volume"] == 9
    assert m["recruitment"]["chemokine_driven"] is True
    assert m["recruitment"]["hill_coefficient_h_m"] == 3
