from pbg_cpm_studies.visualizations import receptor_studies as V


def test_visualizations_render_html():
    for fn in (V.ReceptorRecruitmentScene, V.ReceptorOccupancyLaw, V.ReceptorRecruitmentCurve):
        html = fn()
        assert isinstance(html, str) and "plotly" in html.lower()


def test_scene_is_animated():
    # the headline scene must ship animation frames (Play button + slider)
    html = V.ReceptorRecruitmentScene()
    assert "frames" in html.lower() and "addframes" in html.lower()
