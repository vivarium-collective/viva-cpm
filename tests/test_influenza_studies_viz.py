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


def test_backfilled_study_visualizations_render_html():
    # Increments 3-8: each backfilled study's zero-arg accessor returns a
    # non-empty Plotly card (mirrors the Increment 0-2 check above).
    for fn in (V.InfluenzaIfnResistance, V.InfluenzaEpithelialFate,
               V.InfluenzaMacrophageResponse, V.InfluenzaSignalingFields,
               V.InfluenzaCytotoxicKilling, V.InfluenzaGlobalCoupling):
        html = fn()
        assert isinstance(html, str) and "plotly" in html.lower()
        assert len(html) > 2000


def test_repro_study_visualizations_render_html():
    # Increment 9 (CAPSTONE, Task 9.5): the three repro-fig* dashboard cards
    # -- ensemble/sweep vs digitized acceptance band, baked from a REAL
    # reduced-scale run.repro_fig3b/fig5/fig7 capture (see _influenza_data.py).
    for fn in (V.InfluenzaReproFig3B, V.InfluenzaReproFig5, V.InfluenzaReproFig7):
        html = fn()
        assert isinstance(html, str) and "plotly" in html.lower()
        assert len(html) > 2000


def test_virus_scene_is_animated():
    # the headline scene ships a real top-level Plotly frames array (Play + slider)
    html = V.InfluenzaVirusFieldScene()
    assert html.count('"frames"') > 0


# ── CPM 2D spatial-state videos (the PRIMARY per-study visualization) ────────
_SPATIAL_ACCESSORS = (
    V.InfluenzaSpatialSheet, V.InfluenzaSpatialVirusField,
    V.InfluenzaSpatialIfnResistance, V.InfluenzaSpatialEpithelialFate,
    V.InfluenzaSpatialMacrophage, V.InfluenzaSpatialSignaling,
    V.InfluenzaSpatialCytotoxic, V.InfluenzaSpatialGlobalCoupling,
    V.InfluenzaSpatialReproFig3B, V.InfluenzaSpatialReproFig5,
    V.InfluenzaSpatialReproFig7,
)


def test_spatial_state_videos_render_animated_html():
    # Each study's primary spatial-state video renders non-empty Plotly HTML
    # with a real top-level frames array (▶ Play/❙❙ Pause + MCS slider). Fast:
    # renders from the baked _influenza_spatial blob, no live capture.
    for fn in _SPATIAL_ACCESSORS:
        html = fn()
        assert isinstance(html, str) and "plotly" in html.lower()
        assert html.count('"frames"') > 0
        assert len(html) > 2000


def test_spatial_data_is_real_multiframe():
    # the baked spatial blob has a frame series (>=2 frames) per study, with
    # grids sized to the recorded coarsened lattice (compact, side <= 60).
    from pbg_cpm_studies.visualizations.influenza_studies import _SP
    assert len(_SP) == 9  # 8 distinct scenes + the shared full-model repro scene
    for slug, d in _SP.items():
        assert len(d["frames"]) >= 2
        assert max(d["nx"], d["ny"]) <= 60
        assert len(d["frames"][0]["grid"]) == d["nx"] * d["ny"]


def test_data_is_real_engine_output():
    # baked series match the study's reported values (n_I 9->16, virus->~1216)
    from pbg_cpm_studies.visualizations.influenza_studies import _D
    ser = _D["infection"]["series"]
    assert ser["n_I"][0] == 9 and ser["n_I"][-1] == 16
    assert ser["n_H"][0] + ser["n_I"][0] == ser["n_H"][-1] + ser["n_I"][-1]  # conserved
    assert 1200 < ser["total_virus"][-1] < 1230


def test_backfilled_study_data_matches_reported_values():
    # spot-check the baked Increment 3-8 series against each study.yaml's
    # reported metrics — guards against a stale/garbled data blob.
    from pbg_cpm_studies.visualizations.influenza_studies import _SD
    # Incr 3: IFN gate roughly halves virus load at step 20 (819.16 vs 1750.69)
    assert _SD["ifn"]["with"]["total_virus"][20] < 0.6 * _SD["ifn"]["without"]["total_virus"][20]
    assert 0.5 < _SD["ifn"]["with"]["mean_resist"][40] < 0.6  # plateau ~0.55, not ~1
    # Incr 4: dead lesion forms (n_D 0 -> 14 over 200 updates), Allee death rare
    assert _SD["fate"]["n_D"][0] == 0 and _SD["fate"]["n_D"][-1] == 14
    assert _SD["fate"]["death_cum"][-1] == 0 and _SD["fate"]["recovery_cum"][-1] == 3
    # Incr 5: chemotaxis-on localizes (net approach), lambda=0 control does not
    mac = _SD["macrophage"]
    assert mac["on"]["dist"][-1] < mac["start_distance"] - 10
    # Incr 6: chemokine radial profile decays monotonically (~5.9x inner/outer)
    means = _SD["signaling"]["radial"]["means"]
    assert means == sorted(means, reverse=True) and means[0] / means[-1] > 5
    # Incr 7: close-contact killing clears the cell (1 -> 0), control stays at 1
    assert min(_SD["cytotoxic"]["close"]["n_infected_kill"]) == 0
    assert _SD["cytotoxic"]["close"]["n_infected_nokill"][-1] == 1
    # Incr 8: dynamic sig_1 collapses ~1e-6, far below the retired stub (2.44)
    assert _SD["global"]["sigma1"][-1] < 1e-3


def test_repro_study_data_matches_reported_values():
    # Increment 9 (CAPSTONE): the baked repro_fig3b/fig5/fig7 blobs must
    # match each study.yaml's own reported reduced-scale numbers exactly
    # (same seed0=0 config the studies report) -- guards against a stale/
    # garbled data blob AND against the dashboard card silently disagreeing
    # with the study's written verdict.
    from pbg_cpm_studies.visualizations.influenza_studies import _SD
    # repro-fig3b study.yaml: replicas=2, cells_per_side=15, steps=20,
    # band_eval.passed=False (0/12 observables in-band at reduced scale).
    r3 = _SD["repro_fig3b"]
    assert (r3["replicas"], r3["cells_per_side"], r3["steps"]) == (2, 15, 20)
    assert r3["band_eval"]["passed"] is False
    # repro-fig5-viral-load study.yaml: loads=(1,10000), uninfected_final_frac
    # load=1 -> 0.611, load=10000 -> 0.0, both bands fail.
    r5 = _SD["repro_fig5"]
    assert r5["loads"] == [1, 10000]
    assert abs(r5["by_load"]["1"]["uninfected_final_frac"] - 0.611111) < 1e-4
    assert r5["by_load"]["10000"]["uninfected_final_frac"] == 0.0
    assert r5["by_load"]["1"]["band_eval"]["passed"] is False
    assert r5["by_load"]["10000"]["band_eval"]["passed"] is False
    # repro-fig7-infection-fraction study.yaml: fracs=(0.001,0.05),
    # uninfected_final_frac frac=0.001 -> 1.0, frac=0.05 -> 0.9444.
    r7 = _SD["repro_fig7"]
    assert r7["fracs"] == [0.001, 0.05]
    assert r7["by_frac"]["0.001"]["uninfected_final_frac"] == 1.0
    assert abs(r7["by_frac"]["0.05"]["uninfected_final_frac"] - 0.944444) < 1e-4
