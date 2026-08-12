from pbg_cpm_studies.visualizations import receptor_studies as V


def test_visualizations_render_html():
    for fn in (V.ReceptorRecruitment, V.ReceptorActivationMap):
        html = fn()
        assert isinstance(html, str) and "plotly" in html.lower()
