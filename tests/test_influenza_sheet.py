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


def test_world_builds_and_holds_volume_after_relaxation():
    from pbg_cpm_studies.influenza import build, sheet
    spec = sheet.build_sheet_spec(0.3)
    w = build.world_from_spec(spec)
    assert w.n_cells() == 900
    w.step(50)  # short relaxation
    # cell_volumes()[0] is the medium placeholder, not a real cell (see
    # tests/test_mitosis.py:25 convention). The 0.3mm sheet is exactly
    # confluent (900 * 25 = 22500 = 150x150 sites), so medium volume is
    # legitimately 0 here -- exclude it and check only real epithelial cells.
    vols = w.cell_volumes()[1:]
    mean_v = sum(vols) / len(vols)
    assert 18 <= mean_v <= 32          # ~25 sites, confluent, no collapse
    assert min(vols) > 0               # no cell vanished
