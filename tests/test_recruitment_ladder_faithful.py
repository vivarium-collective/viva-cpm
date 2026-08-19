"""tests/test_recruitment_ladder_faithful.py -- the FAITHFUL occupancy-space
mechanism ladder must separate dramatically (Task R2/R3).

Three responder-level chemotaxis modes on top of the real Rust engine's
occupancy-space chemotaxis (``World.set_chemotaxis_occupancy``, Task R1):

    static_lambda      raw chemotaxis -- background-invariant, recruits at
                        every background (not receptor-mediated).
    hill_occupancy     occupancy chemotaxis, FIXED kd -- saturates and
                        COLLAPSES at high uniform background.
    adaptive_receptor  occupancy chemotaxis, ADAPTIVE kd (re-centred on the
                        locally sensed mean every SAMPLE-step increment) --
                        RESCUES recruitment at high background.

See pbg_cpm_studies/model_building/mechanisms.py for the mechanism
implementation and docs/superpowers/notes/2026-08-18-faithful-ladder-
calibration.md for the empirical calibration sweep that produced the
constants asserted against here.
"""
from pbg_cpm_studies.model_building.mechanisms import simulate_condition


def test_static_recruits_at_high_background():
    # raw chemotaxis is background-invariant -- not receptor-mediated
    assert simulate_condition("static_lambda", "high_bg") > 0.4


def test_hill_occupancy_collapses_at_high_background():
    # fixed-kd occupancy chemotaxis saturates under a high uniform background:
    # theta(dest) - theta(source) -> 0, so the gradient (and recruitment) vanish
    assert simulate_condition("hill_occupancy", "high_bg") < 0.3


def test_adaptive_receptor_rescues_at_high_background():
    # A 3-seed (17,29,43) absolute threshold happened to land favorably
    # (mean 0.833) but a 15-seed resample showed adaptive@high_bg is highly
    # variable per-seed (range 0.0-1.0, mean ~0.49) -- an absolute floor near
    # 0.6 is fragile. The ROBUST result is the COMPARATIVE gap: adaptive
    # reliably and substantially beats hill (which collapses ~deterministically
    # to 0.0 at high background) by a wide margin. Widen the seed set to
    # average out per-seed variance and assert the comparative rescue plus a
    # floor the wider distribution actually supports.
    seeds = tuple(range(11, 23))  # 12 seeds
    hill_hi = simulate_condition("hill_occupancy", "high_bg", seeds=seeds)
    adaptive_hi = simulate_condition("adaptive_receptor", "high_bg", seeds=seeds)
    assert hill_hi < 0.15
    assert adaptive_hi > 0.35
    assert adaptive_hi - hill_hi > 0.3


def test_both_rungs_work_at_low_background():
    assert simulate_condition("adaptive_receptor", "low_bg") > 0.5
    assert simulate_condition("hill_occupancy", "low_bg") > 0.5


def test_receptor_gating_abolishes_hill_recruitment():
    assert simulate_condition("hill_occupancy", "high_bg_blocked") < 0.15
