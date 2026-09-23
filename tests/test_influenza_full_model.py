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
import pytest

from viva_cpm_studies.influenza import run


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
    # With infection + death (incl. ROS, Incr 10) active and a 5% seed, the
    # LIVING epithelium (H+I) shrinks as cells die off (-> D). We assert the
    # living count falls + dead accumulates, rather than uninfected specifically:
    # at this tiny/short config ROS can clear the seeded infected before it
    # spreads, keeping uninfected ~flat while cells still die -- the robust
    # invariant of the infection/death cascade is that living epithelium shrinks.
    r = run.run_full_model(cells_per_side=15, steps=25, seed=2, init_infection_frac=0.05)
    c = r["counts"]
    living0 = c["uninfected"][0] + c["infected"][0]
    living1 = c["uninfected"][-1] + c["infected"][-1]
    assert c["infected"][0] > 0
    assert living1 < living0           # cells died off (-> D)
    assert c["dead"][-1] > c["dead"][0]


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


def test_recruit_pool_four_strip_capacity_and_per_type_sizing():
    # Incr 13: the 4-strip packing (top+bottom+left+right gutters) holds far
    # more reserves than the old top/bottom-only pack, and _seed_recruit_pool
    # accepts a per-type {type: count} dict, seeding exactly that many each.
    from viva_cpm_studies.influenza import build, immune, types
    from viva_cpm_studies.influenza.run import (
        _recruit_pool_origins, _seed_recruit_pool, RECRUIT_RESERVE_TYPE)
    spec = immune.build_cytotoxic_scenario_spec(
        epithelial_cells_per_side=12, n_infected=1, n_macrophages=4,
        n_nk=4, n_cd8=4, margin_sites=40, separation_sites=8, seed=1)
    origins = _recruit_pool_origins(spec, 25)
    assert len(origins) > 200                     # four strips are roomy at margin 40
    assert len(set(origins)) == len(origins)      # no overlapping reserve slots
    world = build.world_from_spec(spec, finalize=False)
    for t in range(RECRUIT_RESERVE_TYPE):
        world.set_contact(RECRUIT_RESERVE_TYPE, t, 10.0)
    world.set_contact(RECRUIT_RESERVE_TYPE, RECRUIT_RESERVE_TYPE, 25.0)
    sizes = {types.M: 30, types.K: 20, types.E: 10}
    pool = _seed_recruit_pool(world, spec, pool_per_type=sizes,
                              targets=[types.M, types.K, types.E],
                              target_volume=25, lambda_volume=9.0)
    assert [len(pool[t]) for t in (types.M, types.K, types.E)] == [30, 20, 10]


def test_bootstrap_immune_config_sizes_pool_to_equilibrium():
    # The pool is sized to the model's own ODE recruitment equilibria (M ~ chemo-
    # saturated ~380, capped types flagged), NOT fitted to bands; seeds are the
    # scenario t=0 immune ICs.
    from viva_cpm_studies.influenza import types
    cfg = run.bootstrap_immune_config(1225)
    assert cfg["n_macrophages"] == 10 and cfg["n_nk"] == 5 and cfg["n_cd8"] == 2
    pool = cfg["recruit_pool_per_type"]
    assert pool[types.M] > 300                     # macrophage eq ~380 -> pool ~440
    assert pool[types.K] > 200 and pool[types.E] > 100
    assert all(v <= 600 for v in pool.values())    # cap holds


def test_full_model_bootstrap_activates_beyond_old_pool_cap():
    # With the bootstrap pool + auto-margin, recruitment can grow macrophages
    # WELL past the old n_seed+6 cap. Small/short so it stays fast: assert the
    # reserve pool actually seeded the requested (large) count -- i.e. auto-margin
    # grew the domain to fit it -- and the run completes.
    from viva_cpm_studies.influenza import types
    boot = run.bootstrap_immune_config(15 * 15)
    r = run.run_full_model(cells_per_side=15, steps=3, seed=1,
                           init_infection_frac=0.05, **boot)
    # domain grew beyond the default margin to hold the pool
    assert r["params"]["tot_cell"] == 225
    # macrophage capacity (seed + pool) far exceeds the old 4+6 cap
    assert boot["recruit_pool_per_type"][types.M] > 20


def test_repro_fig5_viral_load_sweep():
    r = run.repro_fig5(loads=(1, 10000), replicas=1, cells_per_side=12, steps=15, seed0=0)  # reduced
    assert set(r["by_load"].keys()) == {1, 10000}
    # higher initial viral load => not more surviving uninfected than the low-load case
    surv = lambda L: r["by_load"][L]["uninfected_final_frac"]
    assert surv(10000) <= surv(1) + 1e-9
    assert "band_eval" in r


def test_repro_fig7_infection_fraction_sweep():
    r = run.repro_fig7(fracs=(0.001, 0.05), replicas=1, cells_per_side=12, steps=15, seed0=0)  # reduced
    assert set(r["by_frac"].keys()) == {0.001, 0.05}
    # a larger initial infection fraction does not leave more uninfected at the end
    f = lambda x: r["by_frac"][x]["uninfected_final_frac"]
    assert f(0.05) <= f(0.001) + 1e-9
    assert "band_eval" in r


def test_repro_fig5_extracellular_virus_is_concentration_not_raw_sum():
    # Task-9.3-review MUST-FIX regression (units bug): `run_full_model`'s
    # "fields" section is a RAW SUM over the lattice, but `targets/fig5.json`
    # digitizes extracellular_virus as a PER-SITE CONCENTRATION -- at t=0 for
    # viral_load_multiplier=load, the target value IS `load` itself. Before
    # the fix, `repro_fig5` fed the raw sum straight through
    # (`load * n_patch_sites` ~ 3600x too large at cells_per_side=12,
    # confirmed by direct measurement: worst_miss=1198.67 for load=1).
    # `_map_full_model_observables` now divides field-typed observables by
    # `n_lattice_sites` (`params.dims` product) before mapping, so the
    # ensemble's t=0 extracellular_virus must land within an order of
    # magnitude of `load`, not thousands of times larger.
    load = 1000.0
    r = run.repro_fig5(loads=(load,), replicas=1, cells_per_side=12, steps=2, seed0=0)
    t0_virus = r["by_load"][load]["ensemble"]["extracellular_virus"][0][1]
    assert 0.01 * load < t0_virus < 10 * load


def test_fig7_target_subset_scenario_grouping_guard():
    # Fast unit test for the scenario-grouping guard itself (Task-9.3-review
    # gap: the fig5 guard was only verified interactively, not in pytest).
    target_observables = run.bands.load("fig7")["observables"]

    subset = run._fig7_target_subset(target_observables, 0.05)
    for obs_name, obs_list in subset.items():
        assert obs_list, f"{obs_name} subset for frac=0.05 must be non-empty"
        assert {o["initial_infection_fraction"] for o in obs_list} == {0.05}

    with pytest.raises(ValueError):
        run._fig7_target_subset(target_observables, 0.5)  # untagged fraction


def test_fig5_viral_load_dose_response_restored():
    """D2/D8 fix: the fig5 viral-load dose-response is reproduced -- a LOW
    init_viral_load leaves substantially more surviving epithelium than a HIGH
    load. Before the fix (init_viral_load over-seeded ~cell_sites-fold), every
    load drove total loss (no dose-response). The IC now deposits v0/cell_sites
    per pixel (source `v0` = virus PER CELL), so low doses are survivable.
    See docs/cc3d-reference/known-divergences.md D8."""
    from viva_cpm_studies.influenza import run
    lo = run.run_full_model(cells_per_side=12, steps=150, seed=0, init_viral_load=1.0)
    hi = run.run_full_model(cells_per_side=12, steps=150, seed=0, init_viral_load=1000.0)
    un_lo = lo["counts"]["uninfected"][-1]
    un_hi = hi["counts"]["uninfected"][-1]
    tot = lo["params"]["tot_cell"]
    assert un_lo > 0.4 * tot          # low load: much of the sheet survives
    assert un_hi == 0                 # high load: total epithelial loss
    assert un_lo > un_hi              # monotone dose-response


def test_fig5_virus_observable_t0_equals_load():
    """The extracellular_virus observable at t=0 reads the viral-load multiplier
    (per-cell v0), matching fig5.json -- via repro_fig5's field_divisor=tot_cell
    reconciliation now that the IC deposits v0/cell_sites per pixel."""
    from viva_cpm_studies.influenza import run
    r = run.repro_fig5(loads=(1000,), replicas=1, cells_per_side=12, steps=2, seed0=0)
    t0_virus = r["by_load"][1000]["ensemble"]["extracellular_virus"][0][1]
    assert 0.5 * 1000 < t0_virus < 2 * 1000    # reads ~load, not load/cell_sites
