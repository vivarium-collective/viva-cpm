"""Task 7.1: NK (`types.K`) + CD8+ (`types.E`) cells chemotax up the
CHEMOKINE gradient (Increment 6, secreted by macrophages) toward the
infection -- Increment 7 crux (localization only, no killing -- that's
Task 7.2).

Two groups of tests, mirroring Increment 5/6's rigor:
  (a) engine sanity: a minimal world with one chemokine-secreting M cell and
      one K (or E) cell -- the immune cell's distance to the M cell
      decreases once `set_nk_cd8_chemotaxis` is wired, and (lambda=0
      control) does not net-decrease when it is not.
  (b) INTEGRATION (multi-seed, interior placement, unbiased control -- same
      fix-round-1 rigor Increment 5's review caught missing the first time,
      see `immune.build_cytotoxic_scenario_spec`'s docstring): the full
      macrophage + chemokine + NK/CD8 scenario (`run.run_cytotoxic_response`),
      asserted across >= 5 seeds -- NK and CD8 mean distance to the
      infection centroid decreases MORE with chemotaxis on than with the
      lambda=0 control, consistently across seeds; CD8 (lambda=10000)
      localizes at least as strongly as NK (lambda=5000).

Behavior (localization on, CD8 >= NK) is the fidelity criterion here, not
exact distances -- mechanism validation (Increment 7), not a Fig reproduction
(Increment 9). The chemotaxis linear-form gap (engine linear lambda vs the
source's saturating per-cell-COM lambda) is documented in `immune.
set_nk_cd8_chemotaxis`'s docstring, not re-litigated here.
"""
from __future__ import annotations

import math

from viva_cpm_studies.influenza import build, fields, immune, run, types

SEED = 17

# >= 5 seeds, fixed non-tuned knobs (same rigor as Increment 5's fix round 1 --
# `test_influenza_macrophage.py`'s MULTI_SEEDS).
MULTI_SEEDS = (5, 17, 23, 42, 100)

# Engine-sanity lambdas (deliberately NOT `params.yaml`'s literal
# chemotaxis_v_nk=5000/chemotaxis_v_cd8=10000 -- see `run.
# NK_CD8_CHEMOTAXIS_ENGINE_SCALE`'s module docstring): a single M cell's
# chemokine field is ~2-3 orders of magnitude weaker, at these lattice
# distances, than the virus field `test_influenza_macrophage.py`'s
# equivalent sanity check uses, so the literal lambda produces a
# statistically undetectable ΔH per copy attempt here. These values are
# empirically confirmed (task-7.1-report.md) to produce a clean, obvious
# single-trajectory approach in this minimal one-secretor world -- this test
# checks the `set_nk_cd8_chemotaxis` MECHANISM engages for K/E types on the
# chemokine field at all, not the literal default magnitude (that's
# `run.run_cytotoxic_response`'s job, at its own, separately-calibrated
# `NK_CD8_CHEMOTAXIS_ENGINE_SCALE`).
_NK_SANITY_LAMBDA = 5_000_000.0
_CD8_SANITY_LAMBDA = 10_000_000.0


def _minimal_chemokine_world(*, cell_type, chemotaxis_v, warmup=1500, seed=SEED, n=50):
    """A minimal 2-cell world: one M cell (chemokine source, via `fields.
    add_chemokine_field`'s standard M-cell secretion) in one corner, one
    `cell_type` (K or E) cell in the opposite corner. `World.advance_fields`
    warms the field up before any cell steps -- isolating the chemotaxis
    *response* from field build-up time, same convention as
    `test_influenza_macrophage.py`'s `_minimal_world`."""
    spec = {
        "potts": {"dims": [n, n, 1], "boundary": "noflux", "neighbor_order": 3,
                  "temperature": 10.0, "seed": seed},
        "cells": [
            {"type": types.M, "target_volume": 25.0, "lambda_volume": 9.0,
             "target_surface": 0.0, "lambda_surface": 0.0,
             "seed_block": [22, 22, 0, 27, 27, 1]},
            {"type": cell_type, "target_volume": 25.0, "lambda_volume": 9.0,
             "target_surface": 0.0, "lambda_surface": 0.0,
             "seed_block": [3, 3, 0, 8, 8, 1]},
        ],
        "contact": [
            {"a": types.MEDIUM, "b": types.MEDIUM, "j": 0.0},
            {"a": types.MEDIUM, "b": types.M, "j": 10.0},
            {"a": types.MEDIUM, "b": cell_type, "j": 10.0},
            {"a": types.M, "b": cell_type, "j": 20.0},
        ],
    }
    world = build.world_from_spec(spec, finalize=False)
    chemo_fi = fields.add_chemokine_field(world)
    world.finalize(seed)
    for _ in range(warmup):
        world.advance_fields(1)
    if cell_type == types.K:
        immune.set_nk_cd8_chemotaxis(world, chemo_fi, chemotaxis_v_nk=chemotaxis_v,
                                      chemotaxis_v_cd8=0.0)
    else:
        immune.set_nk_cd8_chemotaxis(world, chemo_fi, chemotaxis_v_nk=0.0,
                                      chemotaxis_v_cd8=chemotaxis_v)
    return world


def _distance_to_source(world):
    coms = world.cell_coms()
    mx, my, _ = coms[1]  # cell id 1 == the M (source) cell
    rx, ry, _ = coms[2]  # cell id 2 == the K/E (responder) cell
    return math.hypot(mx - rx, my - ry)


def test_nk_moves_up_chemokine_gradient_engine_sanity():
    world = _minimal_chemokine_world(cell_type=types.K, chemotaxis_v=_NK_SANITY_LAMBDA)
    d0 = _distance_to_source(world)
    for _ in range(20):
        world.step(10)
    d_final = _distance_to_source(world)
    assert d_final < d0, (
        f"NK cell should move UP the chemokine gradient (closer to the M "
        f"source): {d0:.2f} -> {d_final:.2f}"
    )


def test_cd8_moves_up_chemokine_gradient_engine_sanity():
    world = _minimal_chemokine_world(cell_type=types.E, chemotaxis_v=_CD8_SANITY_LAMBDA)
    d0 = _distance_to_source(world)
    for _ in range(20):
        world.step(10)
    d_final = _distance_to_source(world)
    assert d_final < d0, (
        f"CD8+ cell should move UP the chemokine gradient (closer to the M "
        f"source): {d0:.2f} -> {d_final:.2f}"
    )


def test_nk_and_cd8_do_not_approach_source_without_chemotaxis():
    for cell_type in (types.K, types.E):
        world = _minimal_chemokine_world(cell_type=cell_type, chemotaxis_v=0.0)
        d0 = _distance_to_source(world)
        for _ in range(20):
            world.step(10)
        d_final = _distance_to_source(world)
        assert d_final >= d0, (
            f"lambda=0 control (type={cell_type}): should NOT net-approach "
            f"the source (no directed chemotaxis): {d0:.2f} -> {d_final:.2f}"
        )


def test_nk_and_cd8_localize_to_infection_across_seeds():
    """Task 7.1 core integration assertion. Over `MULTI_SEEDS`, with the
    interior-placed macrophage + chemokine + NK/CD8 scenario at fixed,
    non-seed-tuned knobs:
      - chemotaxis-ON's NK and CD8 distance change is more negative than
        chemotaxis-OFF's (lambda=0 control) in EVERY seed;
      - the MEAN on/off gap, for both NK and CD8, robustly exceeds the OFF
        control's own seed-to-seed spread;
      - CD8 (lambda=10000) localizes AT LEAST as strongly as NK (lambda=5000)
        on average (mean CD8 on/off gap >= mean NK on/off gap).
    """
    nk_on_changes, nk_off_changes = [], []
    cd8_on_changes, cd8_off_changes = [], []

    for seed in MULTI_SEEDS:
        on = run.run_cytotoxic_response(seed=seed)
        off = run.run_cytotoxic_response(
            seed=seed, nk_chemotaxis_lambda=0.0, cd8_chemotaxis_lambda=0.0)

        nk_d0 = on["nk_mean_distance_to_infection"][0]
        cd8_d0 = on["cd8_mean_distance_to_infection"][0]
        assert nk_d0 == off["nk_mean_distance_to_infection"][0]  # same seeded scenario
        assert cd8_d0 == off["cd8_mean_distance_to_infection"][0]

        nk_on_change = on["nk_mean_distance_to_infection"][-1] - nk_d0
        nk_off_change = off["nk_mean_distance_to_infection"][-1] - nk_d0
        cd8_on_change = on["cd8_mean_distance_to_infection"][-1] - cd8_d0
        cd8_off_change = off["cd8_mean_distance_to_infection"][-1] - cd8_d0

        nk_on_changes.append(nk_on_change)
        nk_off_changes.append(nk_off_change)
        cd8_on_changes.append(cd8_on_change)
        cd8_off_changes.append(cd8_off_change)

        assert nk_on_change < nk_off_change, (
            f"seed={seed}: NK chemotaxis-on should close MORE of the gap "
            f"than the lambda=0 control: on={nk_on_change:+.2f} "
            f"off={nk_off_change:+.2f}"
        )
        assert cd8_on_change < cd8_off_change, (
            f"seed={seed}: CD8+ chemotaxis-on should close MORE of the gap "
            f"than the lambda=0 control: on={cd8_on_change:+.2f} "
            f"off={cd8_off_change:+.2f}"
        )

    mean_nk_on, mean_nk_off = sum(nk_on_changes) / len(nk_on_changes), sum(nk_off_changes) / len(nk_off_changes)
    mean_cd8_on, mean_cd8_off = sum(cd8_on_changes) / len(cd8_on_changes), sum(cd8_off_changes) / len(cd8_off_changes)
    nk_off_spread = max(nk_off_changes) - min(nk_off_changes)
    cd8_off_spread = max(cd8_off_changes) - min(cd8_off_changes)

    assert mean_nk_on < mean_nk_off - nk_off_spread, (
        f"mean NK on/off gap ({mean_nk_off - mean_nk_on:.2f}) should robustly "
        f"exceed the OFF control's own seed-to-seed spread ({nk_off_spread:.2f}): "
        f"mean_on={mean_nk_on:.2f} mean_off={mean_nk_off:.2f} "
        f"on={nk_on_changes} off={nk_off_changes}"
    )
    assert mean_cd8_on < mean_cd8_off - cd8_off_spread, (
        f"mean CD8 on/off gap ({mean_cd8_off - mean_cd8_on:.2f}) should robustly "
        f"exceed the OFF control's own seed-to-seed spread ({cd8_off_spread:.2f}): "
        f"mean_on={mean_cd8_on:.2f} mean_off={mean_cd8_off:.2f} "
        f"on={cd8_on_changes} off={cd8_off_changes}"
    )

    nk_gap = mean_nk_off - mean_nk_on
    cd8_gap = mean_cd8_off - mean_cd8_on
    assert cd8_gap >= nk_gap, (
        f"CD8+ (lambda=10000) should localize AT LEAST as strongly as NK "
        f"(lambda=5000): nk_gap={nk_gap:.2f} cd8_gap={cd8_gap:.2f}"
    )
