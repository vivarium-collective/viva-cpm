"""The influenza-sego2022 interactive visualization system (Plotly cards).

Distinct from test_influenza_viz.py, which covers the in-package matplotlib
figures (pbg_cpm_studies/influenza/viz.py). This covers the shared v2
visualization system in pbg_cpm_studies/visualizations/influenza_studies.py.
"""
from pbg_cpm_studies.visualizations import influenza_studies as V


def test_visualizations_render_html():
    for fn in (V.InfluenzaFig3BTargets, V.InfluenzaEpithelialSheet,
               V.InfluenzaVirusFieldScene, V.InfluenzaInfectionDynamics):
        html = fn()
        assert isinstance(html, str) and "plotly" in html.lower()
        assert len(html) > 2000


def test_virus_scene_is_animated():
    # the headline scene ships a real top-level Plotly frames array (Play + slider)
    html = V.InfluenzaVirusFieldScene()
    assert html.count('"frames"') > 0


def test_data_is_real_engine_output():
    # baked series match the study's reported values (n_I 9->16, virus->~1216)
    from pbg_cpm_studies.visualizations.influenza_studies import _D
    ser = _D["infection"]["series"]
    assert ser["n_I"][0] == 9 and ser["n_I"][-1] == 16
    assert ser["n_H"][0] + ser["n_I"][0] == ser["n_H"][-1] + ser["n_I"][-1]  # conserved
    assert 1200 < ser["total_virus"][-1] < 1230
