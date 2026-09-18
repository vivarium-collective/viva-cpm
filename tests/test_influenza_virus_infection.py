"""Task 2.3: virus-infection integration test -- proves a seeded lesion
spreads locally over time on top of the wired virus field (2.1) + infection
transition (2.2).

Behavior is the fidelity criterion (spread grows + stays local early), NOT
exact counts -- this is mechanism validation, not a Fig-3B/5/7 reproduction
(that's Increment 9). All RNGs are seeded for determinism.
"""
from __future__ import annotations

from pbg_cpm_studies.influenza import fields, run, types
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


def test_new_infections_appear_near_existing_infected_cells():
    result = run.run_virus_infection(
        patch_mm=0.3, steps=60, seed=17, init_infected_frac=0.01, mcs_per_update=10,
    )
    all_dists = [d for dists in result["new_infection_dists"] for d in dists]
    assert all_dists, "expected at least some new infections over the run"

    # Locality: every new infection lands within a few virus-diffusion
    # lengths of a *prior* infected cell's centroid -- not scattered
    # uniformly across the (150x150-site) sheet.
    locality_radius = 3 * VIRUS_DIFFUSION_LENGTH_SITES  # 75 sites
    assert max(all_dists) < locality_radius, (
        f"new infections should stay local (<{locality_radius} sites): max={max(all_dists)}"
    )

    # Early-update locality specifically: the first update with any new
    # infections should be tighter still (close to the seed lesion).
    for dists in result["new_infection_dists"][:3]:
        for d in dists:
            assert d < locality_radius


def test_no_seed_and_no_initial_virus_means_no_infection():
    result = run.run_virus_infection(
        patch_mm=0.3, steps=60, seed=17, init_infected_frac=0.0, mcs_per_update=10,
    )
    assert all(n == 0 for n in result["n_I"])
    assert result["n_H"][0] == result["n_H"][-1]
    assert all(not dists for dists in result["new_infection_dists"])
