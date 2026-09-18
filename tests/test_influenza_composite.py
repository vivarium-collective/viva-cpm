from pbg_cpm_studies.composites import influenza as inf


def test_composite_document_runs_and_emits():
    doc = inf.epithelial_sheet_baseline()
    assert isinstance(doc, dict)
    # smoke: the demo composite builds a modest sheet, not the full 1mm^2
    spec = inf.build_spec(patch_mm=0.1)   # 50x50 sites, 100 cells
    assert spec["potts"]["dims"] == [50, 50, 1]
    assert len(spec["cells"]) == 100
