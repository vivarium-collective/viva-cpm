from viva_cpm_studies.influenza import bands


def test_series_in_band_basic():
    tgt = [{"t_days": 0.0, "value": 1000, "lo": 900, "hi": 1100},
           {"t_days": 1.0, "value": 300, "lo": 150, "hi": 500}]
    ok = bands.series_in_band([(0.0, 1000), (1.0, 300)], tgt)
    assert ok["in_band"] and ok["n_in"] == 2 and ok["n_checked"] == 2
    bad = bands.series_in_band([(0.0, 1000), (1.0, 50)], tgt)   # 50 < lo 150
    assert not bad["in_band"] and bad["worst_miss"] > 0


def test_soft_band_widens():
    tgt = [{"t_days": 1.0, "value": 300, "lo": 250, "hi": 350}]
    assert not bands.series_in_band([(1.0, 360)], tgt)["in_band"]
    assert bands.series_in_band([(1.0, 360)], tgt, soft=True)["in_band"]  # widened


def test_load_fig3b_has_observables():
    t = bands.load("fig3b")
    assert "uninfected_cells" in t["observables"]
