from pbg_cpm_studies.influenza import sheet, types


def test_sheet_lattice_and_cell_count_0p3mm():
    spec = sheet.build_sheet_spec(0.3)
    # 0.3 mm / 2 um = 150 sites per side; 5x5 cells -> 30 per side -> 900 cells
    assert spec["potts"]["dims"] == [150, 150, 1]
    assert len(spec["cells"]) == 900
    assert all(c["type"] == types.H for c in spec["cells"])
    assert all(c["target_volume"] == 25.0 for c in spec["cells"])
    assert all(c["lambda_volume"] == 9.0 for c in spec["cells"])


def test_sheet_lattice_and_cell_count_1mm():
    spec = sheet.build_sheet_spec(1.0)
    # 1.0 mm / 2 um = 500 sites per side; 5x5 cells -> 100 per side -> 10000 cells
    assert spec["potts"]["dims"] == [500, 500, 1]
    assert len(spec["cells"]) == 10000


def test_sheet_potts_uses_paper_cpm_constants():
    p = sheet.build_sheet_spec(0.3)["potts"]
    assert p["temperature"] == 10
    assert p["neighbor_order"] == 3   # CONTROLLER OVERRIDE: source-literal (was 2)
    assert p["boundary"] in ("noflux", "periodic")
