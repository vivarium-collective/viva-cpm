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


def test_params_carry_chemokine_and_il10_sections():
    p = params.load_params(); c = p["chemokine"]; l = p["il10"]
    for k in ("diffusion_lat2_per_mcs","decay_per_mcs","diffusion_length_cell_diam","b_c","source"):
        assert k in c, f"chemokine missing {k}"
    assert c["diffusion_length_cell_diam"] == 10
    for k in ("diffusion_lat2_per_mcs","decay_per_mcs","diffusion_length_cell_diam","b_l","b_lh","mu_l","sig_1","g_1","g_2","d_2","source"):
        assert k in l, f"il10 missing {k}"
    assert l["diffusion_length_cell_diam"] == 10


def test_params_carry_nk_and_cd8_sections():
    p = params.load_params(); nk = p["nk"]; cd8 = p["cd8"]
    assert nk.get("chemotaxis_v_nk") == 5000 and "g_ik" in nk and "adhesion" in nk and "source" in nk
    assert cd8.get("chemotaxis_v_cd8") == 10000 and "g_ie" in cd8 and "adhesion" in cd8 and "source" in cd8
    # chemotaxis is up the CHEMOKINE field for both (source `chemotaxis_f_nk`/`chemotaxis_f_cd8` = "chemo")
    assert nk["chemotaxis_field"] == "chemokine" and cd8["chemotaxis_field"] == "chemokine"
    # local contact-kill form: kill_rate = g_i*tot_ec_ODE*srf_immune*cell_resist/cell.volume
    assert nk["adhesion"]["nkcell_nkcell"] == 25.0          # homotypic
    assert nk["adhesion"]["nkcell_cd8tcell"] == 25.0        # heterotypic -- source-literal, NOT paper's collapsed 10
    assert cd8["adhesion"]["cd8tcell_cd8tcell"] == 25.0     # homotypic
    assert nk["volume_sites"] == 25 and cd8["volume_sites"] == 25
    # DISCREPANCY #7: contact-killing multiplies by cell_resist DIRECTLY, not (1 - cell_resist)
    assert nk["killing"]["kill_rate_resist_direct"] is True
    assert cd8["killing"]["kill_rate_resist_direct"] is True
    assert "nearby_term_deferred" in nk["killing"] and "nearby_term_deferred" in cd8["killing"]
    assert nk["recruitment"]["chemokine_driven"] is True
    assert cd8["recruitment"]["chemokine_driven"] is False   # CD8 inflow is APC(P)-driven, not chemokine(C)-driven
