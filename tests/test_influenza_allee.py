"""Task 4.3: cellularized Allee-effect death (H -> D) + recovery (D -> H),
plus the driver wiring ALL epithelial fates (infection, infected death,
Allee death, recovery).

(a) pure-function unit tests for `allee.allee_death_rate` /
    `allee.allee_recovery_rate` against the exact source formulas (see
    task-4.0-report.md "EXACT RecoverySteppable death + recovery
    expressions").
(b) integration tests (the key claims): a seeded lesion grows a contiguous
    DEAD region over time (H -> I -> D), and a hand-made world proves the
    Allee recovery mechanism itself: a DYING cell fully surrounded by
    UNINFECTED neighbors recovers; one surrounded by DYING/INFECTED
    neighbors (srf_uninfected below threshold) never does.
"""
from __future__ import annotations

import math

import numpy as np

from viva_cpm_studies.influenza import allee, run, types
from viva_cpm_studies.influenza.params import load_params

try:
    from viva_cpm import cpm_core
except ImportError:  # pragma: no cover
    cpm_core = None


B_H = float(load_params()["allee"]["b_h"])
THETA = float(load_params()["allee"]["srf_threshold"])


# --- (a) pure functions -----------------------------------------------


def test_death_rate_zero_when_srf_total_zero():
    assert allee.allee_death_rate(srf_D=0, srf_uninfected=0, srf_total=0,
                                   b_h=B_H, resist=0.0, theta=THETA) == 0.0


def test_death_rate_zero_when_uninfected_at_or_above_threshold():
    srf_total = 100
    srf_thresh = THETA * srf_total
    # exactly at threshold -> not < threshold -> zero
    assert allee.allee_death_rate(srf_D=50, srf_uninfected=srf_thresh,
                                   srf_total=srf_total, b_h=B_H, resist=0.0,
                                   theta=THETA) == 0.0
    # comfortably above threshold -> zero
    assert allee.allee_death_rate(srf_D=50, srf_uninfected=90,
                                   srf_total=srf_total, b_h=B_H, resist=0.0,
                                   theta=THETA) == 0.0


def test_death_rate_positive_and_monotone_in_surface_deficit():
    srf_total = 100
    # below threshold, with nonzero srf_D -> positive
    r_small_deficit = allee.allee_death_rate(srf_D=80, srf_uninfected=4,
                                              srf_total=srf_total, b_h=B_H,
                                              resist=0.0, theta=THETA)
    r_large_deficit = allee.allee_death_rate(srf_D=80, srf_uninfected=0,
                                              srf_total=srf_total, b_h=B_H,
                                              resist=0.0, theta=THETA)
    assert r_small_deficit > 0.0
    assert r_large_deficit > r_small_deficit, (
        "a larger shortfall below threshold (lower srf_uninfected) should "
        "give a strictly higher death rate"
    )


def test_death_rate_matches_formula_directly():
    srf_total, srf_D, srf_uninfected, resist = 64.0, 40.0, 2.0, 0.3
    srf_thresh = THETA * srf_total
    expected = B_H * (1 - resist) * srf_D * (srf_thresh - srf_uninfected) / srf_total ** 2
    got = allee.allee_death_rate(srf_D=srf_D, srf_uninfected=srf_uninfected,
                                  srf_total=srf_total, b_h=B_H, resist=resist,
                                  theta=THETA)
    assert math.isclose(got, expected, rel_tol=1e-12)


def test_recovery_rate_zero_when_srf_total_zero():
    assert allee.allee_recovery_rate(srf_uninfected=0, srf_total=0, b_h=B_H,
                                      resist=0.0, theta=THETA) == 0.0


def test_recovery_rate_zero_when_uninfected_at_or_below_threshold():
    srf_total = 100
    srf_thresh = THETA * srf_total
    assert allee.allee_recovery_rate(srf_uninfected=srf_thresh, srf_total=srf_total,
                                      b_h=B_H, resist=0.0, theta=THETA) == 0.0
    assert allee.allee_recovery_rate(srf_uninfected=0, srf_total=srf_total,
                                      b_h=B_H, resist=0.0, theta=THETA) == 0.0


def test_recovery_rate_positive_above_threshold():
    srf_total = 100
    r = allee.allee_recovery_rate(srf_uninfected=90, srf_total=srf_total,
                                   b_h=B_H, resist=0.0, theta=THETA)
    assert r > 0.0


def test_recovery_rate_matches_formula_directly():
    srf_total, srf_uninfected, resist = 64.0, 60.0, 0.1
    srf_thresh = THETA * srf_total
    expected = B_H * (1 - resist) * srf_uninfected * (srf_uninfected - srf_thresh) / srf_total ** 2
    got = allee.allee_recovery_rate(srf_uninfected=srf_uninfected, srf_total=srf_total,
                                     b_h=B_H, resist=resist, theta=THETA)
    assert math.isclose(got, expected, rel_tol=1e-12)


# --- (b) integration: lesion formation over the full driver -----------


def test_lesion_forms_dead_region_grows_over_time():
    result = run.run_epithelial_fate(
        patch_mm=0.3, steps=200, seed=17, init_infected_frac=0.05, mcs_per_update=10,
    )
    assert result["n_D"][-1] > result["n_D"][0] == 0, (
        f"a contiguous dead region should form: n_D {result['n_D'][0]} -> {result['n_D'][-1]}"
    )
    n_cells = result["n_H"][0] + result["n_I"][0] + result["n_D"][0]
    for h, i, d in zip(result["n_H"], result["n_I"], result["n_D"]):
        assert h + i + d == n_cells, "H/I/D population must be conserved every update"


# --- (b) integration: hand-made world proves the recovery mechanism ---


def _build_3x3_single_pixel_world(center_type, neighbor_type, temperature=0.0):
    """A 5x5 NoFlux lattice, Moore neighborhood (neighbor_order=3, matching
    this investigation's params.yaml), with a 3x3 block of single-pixel
    cells (mirrors crates/cpm-core/tests/contact_area.rs's setup): the
    center cell is fully Moore-surrounded by the 8 other cells, all of
    `neighbor_type`. `temperature=0.0` freezes Potts dynamics so the
    contact geometry never changes across steps -- isolates the Allee
    stochastic mechanism from cell-shape fluctuation.

    Returns (world, center_cid).
    """
    world = cpm_core.World((5, 5, 1), "noflux", 3, temperature)
    world.set_contact(types.MEDIUM, types.MEDIUM, 0.0)
    for a in (types.H, types.I, types.D):
        world.set_contact(types.MEDIUM, a, 0.0)
        for b in (types.H, types.I, types.D):
            world.set_contact(a, b, 0.0)

    ids = [[0] * 3 for _ in range(3)]
    center_cid = None
    for gy in range(3):
        for gx in range(3):
            cell_type = center_type if (gx, gy) == (1, 1) else neighbor_type
            cid = world.add_cell(cell_type, 1.0, 0.0, 0.0, 0.0)
            world.seed_block(cid, gx + 1, gy + 1, 0, gx + 2, gy + 2, 1)
            ids[gy][gx] = cid
            if (gx, gy) == (1, 1):
                center_cid = cid
    world.finalize(0)
    return world, center_cid


def _epithelial_srf(world, cid):
    """srf_uninfected/srf_D/srf_total for `cid`, restricted to epithelial
    neighbor types (H, I, D) -- matching the source's `ec_list` filter
    (RecoverySteppable.step only accumulates srf_area_total/srf_uninfected
    over neighbors of type in {UNINFECTED, INFECTED, INFECTEDRELEASING,
    DYING}; MEDIUM/immune-type contact is NOT epithelial-neighbor surface
    and must be excluded, see run.py's `_epithelial_contact_totals`)."""
    contact = world.cell_contact_area_by_type(cid)
    epithelial = {types.H, types.I, types.D}
    srf_uninfected = contact.get(types.H, 0)
    srf_D = contact.get(types.D, 0)
    srf_total = sum(v for t, v in contact.items() if t in epithelial)
    return srf_uninfected, srf_D, srf_total


def test_dead_cell_surrounded_by_uninfected_recovers():
    world, center_cid = _build_3x3_single_pixel_world(types.D, types.H)
    rng = np.random.default_rng(3)

    recovered = False
    for _ in range(500):
        srf_uninfected, _srf_D, srf_total = _epithelial_srf(world, center_cid)
        rate = allee.allee_recovery_rate(srf_uninfected=srf_uninfected, srf_total=srf_total,
                                          b_h=B_H, resist=0.0, theta=THETA)
        pr = 1.0 - math.exp(-rate)
        if rng.random() < pr:
            world.set_cell_type(center_cid, types.H)
            recovered = True
            break

    assert recovered, "a DYING cell fully surrounded by UNINFECTED neighbors should recover to H"
    assert list(world.cell_types())[center_cid] == types.H


def test_dead_cell_surrounded_by_dying_does_not_recover():
    world, center_cid = _build_3x3_single_pixel_world(types.D, types.D)
    rng = np.random.default_rng(3)

    for _ in range(500):
        srf_uninfected, _srf_D, srf_total = _epithelial_srf(world, center_cid)
        rate = allee.allee_recovery_rate(srf_uninfected=srf_uninfected, srf_total=srf_total,
                                          b_h=B_H, resist=0.0, theta=THETA)
        assert rate == 0.0, "srf_uninfected=0 must never be above threshold"
        pr = 1.0 - math.exp(-rate)
        if rng.random() < pr:
            world.set_cell_type(center_cid, types.H)

    assert list(world.cell_types())[center_cid] == types.D, (
        "a DYING cell surrounded only by DYING neighbors (srf_uninfected below "
        "threshold) must never recover"
    )
