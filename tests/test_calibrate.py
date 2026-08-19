"""Calibration-as-a-program: sensitivity screen + common random numbers + refine.

Uses fast synthetic objectives (no engine) so the suite stays quick — the point of
the module is that calibration becomes a cheap deterministic program, not a slow
conversation. The recruitment mechanism ladder is one real caller in production.
"""
import numpy as np

from pbg_cpm_studies.model_building import calibrate as C


def test_screen_ranks_the_influential_parameter_first():
    # y depends strongly on 'a', weakly on 'b', not at all on 'c'
    def obj(p):
        return 10.0 * p["a"] + 0.5 * p["b"] + 0.0 * p["c"]

    bounds = {"a": (0.0, 1.0), "b": (0.0, 1.0), "c": (0.0, 1.0)}
    mu = C.elementary_effects(obj, bounds, trajectories=10, seed=1)
    assert mu["a"] > mu["b"] > mu["c"]
    assert C.rank_influential(mu, rel_threshold=0.1) == ["a"]          # b,c pruned
    assert C.rank_influential(mu, top_k=2) == ["a", "b"]


def test_common_random_numbers_make_a_noisy_comparison_deterministic():
    # a "simulator" whose noise depends only on the seed (not the params): with
    # paired seeds the noise cancels and the comparison is exact + reproducible.
    def sim(p, seed):
        return 3.0 * p["x"] + np.random.default_rng(seed).normal(0, 5.0)

    seeds = (17, 29, 43, 61)
    f = C.common_random_numbers(sim, seeds)
    # deterministic across calls (same seeds every time)
    assert f({"x": 0.2}) == f({"x": 0.2})
    # the paired difference recovers the true slope*Δx exactly (noise is common)
    diff = f({"x": 1.0}) - f({"x": 0.0})
    assert abs(diff - 3.0) < 1e-9


def test_refine_emits_a_signed_margin_table_over_only_the_swept_param():
    # objective peaks at x=0.6; band [0.5,0.7] around the peak value
    def obj(p):
        return 1.0 - abs(p["x"] - 0.6)          # max 1.0 at x=0.6

    bounds = {"x": (0.0, 1.0), "y": (0.0, 1.0)}
    # band centered on the peak value (1.0) so the most-central margin IS the peak
    table = C.refine(obj, bounds, ["x"], band=(0.5, 1.5), levels=11)
    assert abs(table[0]["params"]["x"] - 0.6) < 1e-9 and table[0]["within"]  # best clears band
    assert table[0]["margin"] >= table[-1]["margin"]                  # sorted best-first
    # 'y' was not swept → held at its midpoint in every row
    assert all(r["params"]["y"] == 0.5 for r in table)


def test_calibrate_screens_then_refines_and_counts_calls():
    # 'kd' matters, 'noise_knob' does not; target band on the output
    def obj(p):
        return 0.9 * p["kd"] + 0.0 * p["noise_knob"]

    bounds = {"kd": (0.0, 1.0), "noise_knob": (0.0, 1.0)}
    rep = C.calibrate(obj, bounds, band=(0.6, 1.0), top_k=1,
                      screen_trajectories=6, refine_levels=6)
    assert rep["influential"] == ["kd"]                              # screen found it
    assert rep["within"] and rep["best"]["params"]["kd"] >= 0.66     # refine cleared the band
    # far cheaper than a full 2-D grid would be, and auditable
    assert rep["n_calls"] < 40


def test_calibrate_is_deterministic():
    def obj(p):
        return p["a"] - p["b"]
    bounds = {"a": (0.0, 1.0), "b": (0.0, 1.0)}
    r1 = C.calibrate(obj, bounds, band=(0.0, 2.0), seed=7)
    r2 = C.calibrate(obj, bounds, band=(0.0, 2.0), seed=7)
    assert r1["mu_star"] == r2["mu_star"] and r1["best"]["value"] == r2["best"]["value"]


def test_refine_parallel_matches_serial_and_runs_concurrently():
    import threading
    from pbg_cpm_studies.model_building import calibrate as C

    def obj(p):
        return 1.0 - abs(p["x"] - 0.6)

    bounds = {"x": (0.0, 1.0), "y": (0.0, 1.0)}
    serial = C.refine(obj, bounds, ["x"], band=(0.5, 1.5), levels=7)
    parallel = C.refine(obj, bounds, ["x"], band=(0.5, 1.5), levels=7, max_workers=4)
    assert parallel == serial                                     # same result, same order

    # concurrency proof: a barrier of 4 only releases if 4 evals run at once;
    # a serial map would time out on it.
    barrier = threading.Barrier(4, timeout=5)

    def barrier_obj(p):
        barrier.wait()                                           # deadlocks unless concurrent
        return p["x"]

    tbl = C.refine(barrier_obj, {"x": (0.0, 1.0)}, ["x"], band=(0.0, 1.0),
                   levels=4, max_workers=4)
    assert len(tbl) == 4                                          # all four ran together
