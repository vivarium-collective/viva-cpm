"""Task 2.3: virus-infection integration test -- proves a seeded lesion
spreads locally over time on top of the wired virus field (2.1) + infection
transition (2.2).

Behavior is the fidelity criterion (spread grows + stays local early), NOT
exact counts -- this is mechanism validation, not a Fig-3B/5/7 reproduction
(that's Increment 9). All RNGs are seeded for determinism.
"""
from __future__ import annotations

import math

from pbg_cpm_studies.influenza import build, fields, run, sheet, types
from pbg_cpm_studies.influenza.fields import VIRUS_DIFFUSION_LENGTH_SITES


def test_seeded_lesion_spreads_over_time():
    result = run.run_virus_infection(
        patch_mm=0.3, steps=60, seed=17, init_infected_frac=0.01, mcs_per_update=10,
    )
    assert result["n_I"][-1] > result["n_I"][0], (
        f"infected count should grow: {result['n_I'][0]} -> {result['n_I'][-1]}"
    )
    assert result["n_H"][-1] < result["n_H"][0], (
        f"uninfected count should shrink: {result['n_H'][0]} -> {result['n_H'][-1]}"
    )
    # sanity: population is conserved (no cell disappears/appears in this
    # increment -- only H -> I transitions).
    assert result["n_H"][0] + result["n_I"][0] == result["n_H"][-1] + result["n_I"][-1]


def test_new_infections_land_closer_to_prior_infection_than_random_chance():
    """Locality, made discriminating with a uniform-random null baseline.

    A loose absolute radius doesn't discriminate here: on a 150x150-site
    sheet with only ~10-20 reference infected cells, a uniform-RANDOM new
    infection would *also* land within a large radius (e.g. 75 sites) of
    some infected cell with near-certainty -- that bound can't reject "not
    local". Instead, at each update that produces new infections, compare
    the actual new-infection distances (to the nearest already-infected
    cell) against the MEAN nearest-infected distance over the full pool of
    still-uninfected cells that update (a uniform-random null baseline, with
    no sampling noise of its own since it's a population mean, not a random
    draw). Real (field-driven) spread should land much closer than chance.
    """
    result = run.run_virus_infection(
        patch_mm=0.3, steps=60, seed=17, init_infected_frac=0.01, mcs_per_update=10,
    )
    assert any(result["new_infection_dists"]), "expected at least some new infections over the run"

    # Cell centroids are fixed by the sheet geometry (independent of which
    # cells are H vs I), so rebuild the identically-seeded sheet once to read
    # them off rather than threading positions through run_virus_infection's
    # return value.
    spec = sheet.build_sheet_spec(0.3, seed=17)
    world = build.world_from_spec(spec, finalize=True)
    coms = world.cell_coms()
    n_cells = len(spec["cells"])

    def _dist(a, b):
        return math.hypot(coms[a][0] - coms[b][0], coms[a][1] - coms[b][1])

    def _nearest_dist(cid, reference_ids):
        return min(_dist(cid, other) for other in reference_ids)

    # Per-update null: the mean nearest-prior-infected distance over the
    # FULL still-uninfected population that update (not a small random
    # sample of it), so the baseline itself is low-variance -- the only
    # thing that varies event-to-event is where the ACTUAL new infections
    # land relative to that fixed population mean.
    actual_all, random_all = [], []
    early_actual = []  # first few updates that actually produce new infections

    for i in range(1, len(result["steps"])):
        new_dists = result["new_infection_dists"][i]
        if not new_dists:
            continue

        prior_infected = set(result["infected_ids"][i - 1])
        newly_infected = set(result["infected_ids"][i]) - prior_infected
        assert len(newly_infected) == len(new_dists)

        eligible_h = [cid for cid in range(1, n_cells + 1)
                      if cid not in prior_infected and cid not in newly_infected]
        population_dists = [_nearest_dist(cid, prior_infected) for cid in eligible_h]
        null_mean_this_update = sum(population_dists) / len(population_dists)
        random_all.extend([null_mean_this_update] * len(new_dists))
        actual_all.extend(new_dists)

        if len(early_actual) < 3:
            early_actual.extend(new_dists)

    assert len(actual_all) >= 5, (
        f"expected several new-infection events to compare against chance, got {len(actual_all)}"
    )

    actual_mean = sum(actual_all) / len(actual_all)
    random_mean = sum(random_all) / len(random_all)
    # Observed actual/null ratio is ~0.58-0.67 across seeds/step counts here;
    # a shuffled-field (non-local) sanity variant of this same run instead
    # gives ~1.13 (>= chance, as expected once spatial correlation is
    # destroyed) -- 0.75 sits well below the real signal and well above the
    # non-local control, so it discriminates without being seed-fragile.
    assert actual_mean < 0.75 * random_mean, (
        f"new infections should land closer to existing infection than chance "
        f"placement: actual_mean={actual_mean:.1f}, random_null_mean={random_mean:.1f}"
    )

    # Early-update locality specifically (real indices with actual new
    # infections -- the first one here doesn't happen until update ~21, so a
    # fixed [:3] update-index slice would be vacuous): tighter than the
    # 1-2x diffusion-length regime, close to the seed lesion.
    tight_radius = 2 * VIRUS_DIFFUSION_LENGTH_SITES  # 50 sites
    assert early_actual and all(d < tight_radius for d in early_actual), (
        f"early new infections should be close to the seed lesion "
        f"(<{tight_radius} sites): {early_actual}"
    )


def test_no_seed_and_no_initial_virus_means_no_infection():
    result = run.run_virus_infection(
        patch_mm=0.3, steps=60, seed=17, init_infected_frac=0.0, mcs_per_update=10,
    )
    assert all(n == 0 for n in result["n_I"])
    assert result["n_H"][0] == result["n_H"][-1]
    assert all(not dists for dists in result["new_infection_dists"])
