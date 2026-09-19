"""Task 7.2: NK/CD8 contact-mediated killing of infected cells -- the LOCAL
(surface-contact) term of `ContactKillingSteppable`, using the Increment-4
`cell_contact_area_by_type` contact-geometry primitive.

Two groups of tests:
  (a) pure: `killing.contact_kill_rate`'s algebraic properties (monotone in
      contact area + resist, 0 at zero contact area / zero cell volume).
  (b) INTEGRATION (kept FAST, per the brief -- a single deterministic run, NOT
      another multi-seed multi-minute test like `test_influenza_nk_cd8.py`'s
      `test_nk_and_cd8_localize_to_infection_across_seeds`): a small,
      close-contact variant of `run.run_cytotoxic_response`'s Task-7.1
      scenario -- infected-cell count is LOWER with NK/CD8 killing enabled
      than with it disabled, from the SAME seed (so the only difference
      between the two runs is the killing step itself, not chemotaxis or any
      other mechanism).

Scenario-size note for (b): `run_cytotoxic_response`'s Task-7.1 DEFAULTS
(6 macrophages/6 NK/6 CD8, margin_sites=20, separation_sites=20, steps=60)
were empirically confirmed (ad hoc, not part of this fast suite) to NOT
produce any NK/CD8-infected-cell contact within that step budget -- see
task-7.2-report.md's "weak killing at default scale" finding, the same
category as task-7.1-report.md's chemotaxis-magnitude gap (NK_CD8_CHEMOTAXIS_
ENGINE_SCALE). This test instead uses a SMALL, tight-margin/separation variant
of the identical scenario/driver (fewer cells, smaller domain) so the NK/CD8
cluster can physically reach the infected cell within a fast step budget --
this exercises the contact-kill MECHANISM itself (this task's scope), not the
separately-flagged chemotaxis-reach question (deferred to Incr 9 calibration).

Task 8.5 adds two more groups, exercising `run.run_global_coupling` instead
(the Increment-8 hybrid-coupling driver):
  (c) pure: `killing.nearby_kill_rate`'s algebraic properties (DIRECT in
      cell_resist -- discrepancy #7 -- and 0 at `num_nearby<=0`/`eta<=0`).
  (d) INTEGRATION (kept FAST -- a single small, deterministic
      `run_global_coupling` run, NOT a multi-seed/paper-scale run):
      `run_global_coupling`'s NEARBY-population term is empirically too weak
      at this engine's default small-patch scale (`eta` tiny) to fire within
      any FAST step budget via the natural Task-8.4 recruitment buildup alone
      -- the SAME "weak effect at default scale" category as (b)'s note
      above and task-7.1-report.md's NK_CD8_CHEMOTAXIS_ENGINE_SCALE finding
      (see task-8.5-report.md). `run_global_coupling`'s `initial_K_nb`/
      `initial_E_nb` (a documented TEST SEAM, NOT a source parameter or a
      tuned rate -- see that function's docstring) instead seed a large
      synthetic nearby population directly, so this test exercises the
      `nearby_kill_rate` WIRING itself within a fast budget.
"""
from __future__ import annotations

import math

from pbg_cpm_studies.influenza import killing, run
from pbg_cpm_studies.influenza.params import load_params

SEED = 17

# Small, close-contact scenario knobs: fixed (not per-seed-tuned), chosen so
# the NK/CD8 cluster is within physical reach of the single infected cell
# within a fast (a few seconds) step budget -- see module docstring.
_SMALL_SCENARIO = dict(
    epithelial_cells_per_side=2, n_infected=1,
    n_macrophages=2, n_nk=2, n_cd8=2,
    margin_sites=4, separation_sites=1,
    steps=45, seed=SEED,
)


# --- (a) pure contact_kill_rate ---------------------------------------------

def test_contact_kill_rate_zero_at_zero_contact_area():
    assert killing.contact_kill_rate(0.0, 0.5, 1e-6, 250000.0, 25.0) == 0.0


def test_contact_kill_rate_zero_at_zero_cell_volume():
    assert killing.contact_kill_rate(5.0, 0.5, 1e-6, 250000.0, 0.0) == 0.0


def test_contact_kill_rate_nonnegative():
    rate = killing.contact_kill_rate(3.0, 0.2, 2.14e-8, 250000.0, 25.0)
    assert rate >= 0.0


def test_contact_kill_rate_monotone_in_contact_area():
    g_i, tot_ec, cell_volume, resist = 2.14e-8, 250000.0, 25.0, 0.3
    rates = [killing.contact_kill_rate(srf, resist, g_i, tot_ec, cell_volume)
              for srf in (0.0, 1.0, 3.0, 8.0, 15.0)]
    assert rates == sorted(rates)
    assert rates[0] < rates[-1]


def test_contact_kill_rate_monotone_in_cell_resist():
    g_i, tot_ec, cell_volume, srf = 2.14e-8, 250000.0, 25.0, 5.0
    rates = [killing.contact_kill_rate(srf, resist, g_i, tot_ec, cell_volume)
              for resist in (0.0, 0.1, 0.5, 0.9, 1.0)]
    assert rates == sorted(rates)
    assert rates[0] == 0.0  # resist == 0 -> DIRECT multiplier -> rate 0 (discrepancy #7)
    assert rates[-1] > 0.0


def test_contact_kill_rate_matches_params_yaml_precomputed_coefficients():
    """Sanity cross-check against params.yaml's recorded `pr_cf_nk_loc` /
    `pr_cf_cd8_loc` (== `g_i * tot_ec_ODE`, the source's own precomputed
    LOCAL killing-probability coefficient, `ContactKillingSteppable.start()`):
    at `srf_immune=cell_resist=cell_volume=1`, `contact_kill_rate` reduces to
    exactly `g_i * tot_ec`.
    """
    params = load_params()
    g_ik = float(params["nk"]["g_ik"])
    g_ie = float(params["cd8"]["g_ie"])
    tot_ec = float(params["scaling"]["ode_epithelial_population"])
    pr_cf_nk_loc = float(params["nk"]["killing"]["pr_cf_nk_loc"])
    pr_cf_cd8_loc = float(params["cd8"]["killing"]["pr_cf_cd8_loc"])

    rate_nk = killing.contact_kill_rate(1.0, 1.0, g_ik, tot_ec, 1.0)
    rate_cd8 = killing.contact_kill_rate(1.0, 1.0, g_ie, tot_ec, 1.0)
    assert math.isclose(rate_nk, pr_cf_nk_loc, rel_tol=1e-9)
    assert math.isclose(rate_cd8, pr_cf_cd8_loc, rel_tol=1e-9)


# --- (b) INTEGRATION: infected count lower with NK/CD8 killing -------------

def test_infected_count_lower_with_nk_cd8_killing_than_without():
    """Core Task 7.2 integration assertion: from the SAME seed, the small
    close-contact scenario (module docstring), the run WITH NK/CD8 killing
    enabled ends with a LOWER (or equal, if the stochastic draw happens not
    to fire) infected-cell count than the run with killing disabled -- and,
    at this fixed scenario/seed, strictly lower (the infected cell is
    actually killed once NK/CD8 make contact).
    """
    with_killing = run.run_cytotoxic_response(**_SMALL_SCENARIO, enable_killing=True)
    without_killing = run.run_cytotoxic_response(**_SMALL_SCENARIO, enable_killing=False)

    # Same seeded initial scenario -> identical starting infected count.
    assert with_killing["n_infected"][0] == without_killing["n_infected"][0] == 1

    # Without killing: NK/CD8 present + chemotaxing, but the kill step is
    # skipped -- the infected cell survives the whole run (no other
    # mechanism in this driver can remove it).
    assert without_killing["n_infected"][-1] == 1

    # With killing: the infected-cell count at the final step is lower than
    # the no-killing control from the same seed -- NK/CD8 contact-kill fired.
    assert with_killing["n_infected"][-1] < without_killing["n_infected"][-1]
    assert with_killing["n_infected"][-1] == 0


def test_infected_count_never_increases_with_killing_enabled():
    """`n_infected` is monotonically non-increasing across the run: this
    driver has no H -> I infection mechanism (fixed, pre-seeded infected
    cells only, per the brief's "at minimum, seed infected cells for NK/CD8
    to kill" minimum bar), and killing only ever removes cells from the
    infected set."""
    result = run.run_cytotoxic_response(**_SMALL_SCENARIO, enable_killing=True)
    counts = result["n_infected"]
    assert all(counts[i] >= counts[i + 1] for i in range(len(counts) - 1))


# --- (c) pure nearby_kill_rate (Task 8.5) -----------------------------------

def test_nearby_kill_rate_form():
    r = killing.nearby_kill_rate(g_i=1e-3, eta=0.0049, num_nearby=4.0, cell_resist=0.5)
    assert math.isclose(r, (1e-3/0.0049)*4.0*0.5)
    assert killing.nearby_kill_rate(1e-3, 0.0049, 0.0, 0.5) == 0.0   # no nearby -> 0
    assert killing.nearby_kill_rate(1e-3, 0.0, 4.0, 0.5) == 0.0      # eta<=0 -> 0


def test_nearby_kill_rate_resist_direct_discrepancy7():
    # DIRECT in resist (NOT 1-resist): rate rises with resist
    assert killing.nearby_kill_rate(1e-3,0.0049,4.0,0.9) > killing.nearby_kill_rate(1e-3,0.0049,4.0,0.1)
    assert killing.nearby_kill_rate(1e-3,0.0049,4.0,0.0) == 0.0


def test_killing_disabled_control_holds_identical_scenario():
    """`enable_killing=False` keeps NK/CD8 cells present and still
    chemotaxing (same `params` scenario knobs) -- only the killing step
    itself is skipped, per the brief's "NK/CD8 present but killing disabled"
    construction for an isolated on/off comparison."""
    with_killing = run.run_cytotoxic_response(**_SMALL_SCENARIO, enable_killing=True)
    without_killing = run.run_cytotoxic_response(**_SMALL_SCENARIO, enable_killing=False)
    for key in ("chemotaxis_v_macro", "chemotaxis_v_nk", "chemotaxis_v_cd8",
                "n_macrophages", "n_nk", "n_cd8", "n_infected"):
        assert with_killing["params"][key] == without_killing["params"][key]
    assert with_killing["params"]["enable_killing"] is True
    assert without_killing["params"]["enable_killing"] is False


# --- (d) INTEGRATION: run_global_coupling NEARBY term + DH-input fix -------

# Task 8.5 test-seam scenario knobs (module docstring, group (d)): a large
# synthetic nearby NK/CD8 surrogate population (`initial_K_nb`/`initial_E_nb`
# -- NOT a tuned rate constant, see `run.run_global_coupling`'s docstring)
# large enough to make `killing.nearby_kill_rate`'s contribution dominate
# within a few MCS, at this driver's small default `side=30` patch scale.
_NEARBY_SCENARIO = dict(side=30, steps=10, seed=1, with_immune=True,
                         with_recruitment=True)


def test_global_coupling_nearby_killing_strengthens_infected_clearance():
    """Core Task 8.5 integration assertion: from the SAME seed/scenario, a
    large nearby NK/CD8 surrogate population (K_nb, E_nb > 0, via the
    `initial_K_nb`/`initial_E_nb` test seam) clears infected cells FASTER
    than the K_nb=E_nb=0 control -- the nearby term
    (`killing.nearby_kill_rate`) actually does work once wired into
    `run.run_global_coupling`."""
    with_nearby = run.run_global_coupling(**_NEARBY_SCENARIO,
                                           initial_K_nb=2000.0, initial_E_nb=100.0)
    without_nearby = run.run_global_coupling(**_NEARBY_SCENARIO,
                                              initial_K_nb=0.0, initial_E_nb=0.0)

    # Same seeded scenario -> identical seeded infected count.
    assert with_nearby["params"]["n_infected"] == without_nearby["params"]["n_infected"] > 0

    # Without a nearby population, the (already-weak-at-this-scale) LOCAL
    # term alone does not clear the infected cells within this fast budget.
    assert without_nearby["params"]["n_infected_final"] == without_nearby["params"]["n_infected"]
    assert without_nearby["params"]["dead_from_infected"] == 0

    # With a nearby population present, clearance is strictly stronger.
    assert with_nearby["params"]["n_infected_final"] < without_nearby["params"]["n_infected_final"]
    assert with_nearby["spatial"]["I"][-1] < without_nearby["spatial"]["I"][-1]
    assert with_nearby["params"]["dead_from_infected"] > 0


def test_global_coupling_dh_input_does_not_rise_when_infected_cells_are_killed():
    """Task-8.2-review DH-fix (this function's docstring, "DH-input fix"
    paragraph): a cell killed by Task 8.5 (I -> D) is dead-from-INFECTED, NOT
    dead-from-HEALTHY -- the spatial->ODE `DH` input must NOT rise as a
    result. This driver has no H -> D (ROS/Allee) mechanism, so DH must stay
    at its pre-kill value (0) throughout, even while killing is actively
    firing (`dead_from_infected > 0`, i.e. this assertion is not vacuous)."""
    result = run.run_global_coupling(**_NEARBY_SCENARIO,
                                      initial_K_nb=2000.0, initial_E_nb=100.0)
    assert result["params"]["dead_from_infected"] > 0
    assert result["params"]["n_infected_final"] < result["params"]["n_infected"]
    assert all(dh == 0 for dh in result["spatial"]["DH"])
