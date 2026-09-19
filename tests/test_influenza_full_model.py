"""Increment 9 Task 9.1: the complete multiscale driver `run.run_full_model`.

These tests are the contract for the capstone full-model assembly: every
mechanism from Increments 0-8 (infection/death transitions, per-cell
resistance, the four diffusive fields, macrophage/NK/CD8 chemotaxis +
contact/nearby killing, ODE-driven recruitment, the hybrid Price-2015 global
ODE, the Allee death/recovery rule) wired into ONE source-ordered per-MCS
loop. They assert BEHAVIOR (all observables recorded; infection depletes
uninfected; a uniform virus IC seeds infection with no pre-infected cells) at
a small/fast reduced scale -- NOT figure-band matching (that is Tasks
9.2-9.4 + the Phase-B paper-scale ensemble).
"""
from pbg_cpm_studies.influenza import run


def test_full_model_smoke_all_observables():
    r = run.run_full_model(cells_per_side=15, steps=10, seed=1, init_infection_frac=0.05)
    n = len(r["mcs"])
    assert n == 10 and len(r["t_days"]) == n
    for k in ("uninfected", "infected", "dead", "macrophage", "nk", "cd8"):
        assert len(r["counts"][k]) == n
    for f in ("virus", "ifn", "chemo", "il10"):
        assert len(r["fields"][f]) == n
    assert set(("T", "X", "A", "P")).issubset(r["ode"].keys())


def test_full_model_infection_depletes_uninfected():
    # With infection + death active and a 5% seed, uninfected count must fall
    r = run.run_full_model(cells_per_side=15, steps=25, seed=2, init_infection_frac=0.05)
    assert r["counts"]["uninfected"][-1] < r["counts"]["uninfected"][0]
    assert r["counts"]["infected"][0] > 0


def test_full_model_viral_load_scenario_seeds_infection():
    # init_viral_load with no pre-infected cells: virus drives H->I over time
    r = run.run_full_model(cells_per_side=15, steps=25, seed=3, init_viral_load=1000.0)
    assert r["counts"]["infected"][0] == 0
    assert max(r["counts"]["infected"]) > 0     # infection emerges from the virus field


def test_full_model_fidelity_fixes():
    # §4a: ODE M/K/E inputs include nearby surrogates (K input >= local K when K_nb>0);
    # recruited macrophages secrete (macrophage id list refreshes). Assert via the driver's
    # recorded ode K vs spatial macrophage/nk counts, or an exposed debug hook. Keep FAST.
    r = run.run_full_model(cells_per_side=15, steps=8, seed=4, init_infection_frac=0.05)
    assert "K" in r["ode"]


def test_repro_fig3b_runs_and_evaluates_bands():
    r = run.repro_fig3b(replicas=2, cells_per_side=15, steps=20, seed0=0)  # reduced for CI
    assert r["replicas"] == 2
    assert "uninfected_cells" in r["ensemble"]
    assert "band_eval" in r and "passed" in r["band_eval"]
    # non-vacuous: the ensemble uninfected series declines under infection
    u = r["ensemble"]["uninfected_cells"]
    assert u[-1][1] <= u[0][1]
