"""Task 5.1: macrophage in the world + chemotaxis up the virus field.

Increment-5 crux -- proves the macrophage cell type (`types.M`) chemotaxes up
the extracellular-virus gradient toward the infection (paper Fig 2B / Fig 3A
-- macrophages accumulate at the lesion), using the existing engine
chemotaxis primitive (`World.set_chemotaxis`) wired the same way
`viva_cpm_studies.chemotaxis` wires its source/responder recipe -- no Rust
change.

Two tests:
  (a) engine sanity: a minimal 2-cell world (one virus-secreting I cell, one
      chemotaxing M cell) -- the M cell's distance to the I cell decreases
      once `set_macrophage_chemotaxis` is wired, and (lambda=0 control)
      increases (thermal/wall drift) when it is not. Fully deterministic,
      single-seed -- fine as a unit check of the engine primitive itself.
  (b) integration (Task 5.1 review fix round 1): the full non-confluent
      macrophage-response scenario (`immune.build_macrophage_scenario_spec`
      + `run.run_macrophage_response`), asserted across MULTIPLE seeds --
      macrophages' mean distance to the infection centroid decreases MORE
      with chemotaxis on than with chemotaxis off (lambda=0 control),
      consistently across seeds, not just at one hand-tuned seed.

## Fix round 1 (review findings)

An earlier version of (b) ran ONLY at seed=17 with a macrophage cluster
pinned in a domain CORNER diagonally opposite the (centered) infection.
Two problems, both fixed here:

1. CONTROL-DRIFT CONFOUND -- a corner is a hard noflux-wall reflector on TWO
   sides at once, so the lambda=0 control's thermal/wall drift was not an
   isotropic random walk: it was biased to net-drift away from that corner,
   which happened to point roughly toward the (diagonally opposite)
   infection, inflating the apparent on/off gap. `immune.
   build_macrophage_scenario_spec` now places BOTH the epithelial/infection
   patch and the macrophage cluster deep in the domain's INTERIOR (each
   `margin_sites` from every wall) -- see that function's docstring.
2. SINGLE-SEED / knob-tuning -- the old run-length/warmup/margin values were
   effectively fit to seed 17's own trajectory. `run.run_macrophage_response`
   now uses fixed, physically-motivated defaults (`field_warmup=3000` is
   close to the virus field's own decay relaxation time, `sqrt(D/decay) ~ 25
   sites` / `1/decay ~ 3500 MCS`; margin/separation set the macrophages'
   start well outside that length scale) that were NOT re-tuned per seed,
   and the integration assertion below checks 5 seeds, not 1.

Behavior (localization on, no consistent localization in the lambda=0
control) is the fidelity criterion, not exact distances -- this is mechanism
validation (Increment 5), not a Fig-2B/3A reproduction (Increment 9).
"""
from __future__ import annotations

import math

from viva_cpm_studies.influenza import build, fields, immune, run, types

SEED = 17

# Task 5.1 fix round 1: >= 5 seeds, none of them the one any knob was tuned
# against (the geometry/timing knobs were fixed from physical scales -- see
# module docstring -- before this seed set was ever run).
MULTI_SEEDS = (5, 17, 23, 42, 100)


def _minimal_world(*, chemotaxis_v_macro, warmup=1500, seed=SEED, n=50):
    """A minimal 2-cell world: one I cell (virus source, via `fields.
    add_virus_field`'s standard I-cell secretion) in one corner, one M cell
    in the opposite corner. `World.advance_fields` warms the field up
    (secretion + diffusion + decay, no Potts motion) before any cell steps,
    so the gradient already reaches the M cell by step 0 -- isolating the
    chemotaxis *response* from field build-up time. (This minimal check
    deliberately keeps its own small, self-contained geometry -- it is not
    exercising `immune.build_macrophage_scenario_spec` and is not subject to
    the interior-placement fix that scenario now requires; it's a unit
    check of `World.set_chemotaxis` on type M, not a control-drift claim.)"""
    spec = {
        "potts": {"dims": [n, n, 1], "boundary": "noflux", "neighbor_order": 3,
                  "temperature": 10.0, "seed": seed},
        "cells": [
            {"type": types.I, "target_volume": 25.0, "lambda_volume": 9.0,
             "target_surface": 0.0, "lambda_surface": 0.0,
             "seed_block": [22, 22, 0, 27, 27, 1]},
            {"type": types.M, "target_volume": 25.0, "lambda_volume": 9.0,
             "target_surface": 0.0, "lambda_surface": 0.0,
             "seed_block": [3, 3, 0, 8, 8, 1]},
        ],
        "contact": [
            {"a": types.MEDIUM, "b": types.MEDIUM, "j": 0.0},
            {"a": types.MEDIUM, "b": types.I, "j": 25.0},
            {"a": types.MEDIUM, "b": types.M, "j": 10.0},
            {"a": types.I, "b": types.M, "j": 20.0},
        ],
    }
    world = build.world_from_spec(spec, finalize=False)
    virus_fi = fields.add_virus_field(world)
    world.finalize(seed)
    for _ in range(warmup):
        world.advance_fields(1)
    immune.set_macrophage_chemotaxis(world, virus_fi, chemotaxis_v_macro=chemotaxis_v_macro)
    return world


def _m_to_i_distance(world):
    coms = world.cell_coms()
    ix, iy, _ = coms[1]  # cell id 1 == the I (source) cell
    mx, my, _ = coms[2]  # cell id 2 == the M cell
    return math.hypot(ix - mx, iy - my)


def test_macrophage_moves_up_virus_gradient_engine_sanity():
    world = _minimal_world(chemotaxis_v_macro=5000.0)
    d0 = _m_to_i_distance(world)
    for _ in range(20):
        world.step(10)
    d_final = _m_to_i_distance(world)
    assert d_final < d0, (
        f"macrophage should move UP the virus gradient (closer to the I "
        f"source): {d0:.2f} -> {d_final:.2f}"
    )


def test_macrophage_does_not_approach_source_without_chemotaxis():
    world = _minimal_world(chemotaxis_v_macro=0.0)
    d0 = _m_to_i_distance(world)
    for _ in range(20):
        world.step(10)
    d_final = _m_to_i_distance(world)
    assert d_final >= d0, (
        f"lambda=0 control: macrophage should NOT net-approach the source "
        f"(no directed chemotaxis): {d0:.2f} -> {d_final:.2f}"
    )


def test_macrophages_localize_to_infection_across_seeds():
    """Task 5.1 fix round 1's core integration assertion. Over
    `MULTI_SEEDS`, with the (now interior-placed) macrophage scenario at
    fixed, non-seed-tuned knobs:
      - chemotaxis-ON's distance change is negative (net approach) in EVERY
        seed;
      - chemotaxis-ON's distance change is more negative than chemotaxis-
        OFF's in EVERY seed (the on/off gap holds seed-by-seed, not just on
        average);
      - the MEAN on/off gap is large relative to the OFF control's own
        seed-to-seed spread (i.e. not explainable by control noise alone).
    """
    on_changes = []
    off_changes = []
    for seed in MULTI_SEEDS:
        on = run.run_macrophage_response(chemotaxis_lambda=None, seed=seed)
        off = run.run_macrophage_response(chemotaxis_lambda=0.0, seed=seed)
        d0 = on["mean_distance_to_infection"][0]
        assert d0 == off["mean_distance_to_infection"][0]  # same seeded scenario

        on_change = on["mean_distance_to_infection"][-1] - d0
        off_change = off["mean_distance_to_infection"][-1] - d0
        on_changes.append(on_change)
        off_changes.append(off_change)

        assert on_change < 0, (
            f"seed={seed}: chemotaxis-on should net-approach the infection "
            f"(negative distance change), got {on_change:+.2f}"
        )
        assert on_change < off_change, (
            f"seed={seed}: chemotaxis-on should close MORE of the gap than "
            f"the lambda=0 control: on={on_change:+.2f} off={off_change:+.2f}"
        )

    mean_on = sum(on_changes) / len(on_changes)
    mean_off = sum(off_changes) / len(off_changes)
    off_spread = max(off_changes) - min(off_changes)

    assert mean_on < mean_off - off_spread, (
        f"mean on/off gap ({mean_off - mean_on:.2f}) should robustly exceed "
        f"the OFF control's own seed-to-seed spread ({off_spread:.2f}): "
        f"mean_on={mean_on:.2f} mean_off={mean_off:.2f} "
        f"on_changes={on_changes} off_changes={off_changes}"
    )
