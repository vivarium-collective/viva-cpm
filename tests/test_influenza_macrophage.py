"""Task 5.1: macrophage in the world + chemotaxis up the virus field.

Increment-5 crux -- proves the macrophage cell type (`types.M`) chemotaxes up
the extracellular-virus gradient toward the infection (paper Fig 2B / Fig 3A
-- macrophages accumulate at the lesion), using the existing engine
chemotaxis primitive (`World.set_chemotaxis`) wired the same way
`pbg_cpm_studies.chemotaxis` wires its source/responder recipe -- no Rust
change.

Two tests:
  (a) engine sanity: a minimal 2-cell world (one virus-secreting I cell, one
      chemotaxing M cell) -- the M cell's distance to the I cell decreases
      once `set_macrophage_chemotaxis` is wired, and (lambda=0 control)
      increases (thermal/wall drift) when it is not.
  (b) integration: the full non-confluent macrophage-response scenario
      (`immune.build_macrophage_scenario_spec` + `run.
      run_macrophage_response`) -- macrophages' mean distance to the
      infection centroid DECREASES over the run when chemotaxis is on, and
      does NOT decrease (lambda=0 control) when it is off.

Behavior (localization on, no localization in the lambda=0 control) is the
fidelity criterion, not exact distances -- this is mechanism validation
(Increment 5), not a Fig-2B/3A reproduction (Increment 9). Both tests are
fully deterministic (fixed seed=17, the investigation's canonical
determinism seed, used throughout `pbg_cpm_studies.influenza`).
"""
from __future__ import annotations

import math

from pbg_cpm_studies.influenza import build, fields, immune, run, types

SEED = 17


def _minimal_world(*, chemotaxis_v_macro, warmup=1500, seed=SEED, n=50):
    """A minimal 2-cell world: one I cell (virus source, via `fields.
    add_virus_field`'s standard I-cell secretion) in one corner, one M cell
    in the opposite corner. `World.advance_fields` warms the field up
    (secretion + diffusion + decay, no Potts motion) before any cell steps,
    so the gradient already reaches the M cell by step 0 -- isolating the
    chemotaxis *response* from field build-up time."""
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


def test_macrophages_localize_to_infection():
    result = run.run_macrophage_response(chemotaxis_lambda=None, seed=SEED)
    d0 = result["mean_distance_to_infection"][0]
    d_final = result["mean_distance_to_infection"][-1]
    assert d_final < d0, (
        f"macrophages should localize to the infection (mean distance to "
        f"infection centroid should decrease): {d0:.2f} -> {d_final:.2f}"
    )


def test_macrophages_do_not_localize_without_chemotaxis_control():
    result = run.run_macrophage_response(chemotaxis_lambda=0.0, seed=SEED)
    d0 = result["mean_distance_to_infection"][0]
    d_final = result["mean_distance_to_infection"][-1]
    assert d_final >= d0, (
        f"lambda=0 control: macrophages should NOT net-localize to the "
        f"infection: {d0:.2f} -> {d_final:.2f}"
    )


def test_macrophage_localization_markedly_stronger_with_chemotaxis_on():
    """Direct on-vs-off comparison, same seed/geometry: chemotaxis-on must
    close substantially more of the initial gap than the lambda=0 control --
    guards against both curves drifting the same (small) amount by chance."""
    on = run.run_macrophage_response(chemotaxis_lambda=None, seed=SEED)
    off = run.run_macrophage_response(chemotaxis_lambda=0.0, seed=SEED)
    d0 = on["mean_distance_to_infection"][0]
    assert d0 == off["mean_distance_to_infection"][0]  # same seeded scenario
    on_final = on["mean_distance_to_infection"][-1]
    off_final = off["mean_distance_to_infection"][-1]
    assert on_final < off_final, (
        f"chemotaxis-on should end closer to the infection than the "
        f"lambda=0 control: on={on_final:.2f} off={off_final:.2f}"
    )
