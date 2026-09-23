import pytest
from viva_cpm_studies.influenza import targets


@pytest.mark.parametrize("name", ["fig3b", "fig5", "fig7"])
def test_target_loads_and_has_bands(name):
    t = targets.load_target(name)
    assert t["figure"] == name
    assert t["observables"], "at least one observable"
    for obs, series in t["observables"].items():
        assert series, f"{obs} has points"
        for pt in series:
            assert pt["lo"] <= pt["value"] <= pt["hi"], f"{obs} band brackets value"
            assert pt["t_days"] >= 0


def test_fig3b_covers_key_observables():
    obs = targets.load_target("fig3b")["observables"]
    for key in ("uninfected_cells", "infected_cells", "extracellular_virus"):
        assert key in obs
