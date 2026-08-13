from pbg_cpm_studies.visualizations import receptor_studies as V


def test_visualizations_render_html():
    for fn in (V.ReceptorRecruitmentScene, V.ReceptorOccupancyLaw, V.ReceptorRecruitmentCurve):
        html = fn()
        assert isinstance(html, str) and "plotly" in html.lower()


def test_scene_is_animated():
    # the headline scene must ship a real top-level Plotly frames array (not
    # every timepoint overlaid as static traces) driving the Play button + slider
    html = V.ReceptorRecruitmentScene()
    assert html.count('"frames"') > 0
