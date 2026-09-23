"""Task 2.2: pure, stochastic virus-driven infection transition H -> I.

Also covers Task 4.2: pure, stochastic infected-cell death transition I -> D.
"""
from __future__ import annotations

import math

import numpy as np

from viva_cpm_studies.influenza import transitions, types


def test_high_virus_infects_with_near_certain_probability():
    types_list = [types.MEDIUM, types.H]
    virus_at_cell = [0.0, 1000.0]  # rate = g_hv * v_bar >> 1
    rng = np.random.default_rng(0)
    out = transitions.infection_step(types_list, virus_at_cell, g_hv=0.1, rng=rng)
    assert out[1] == types.I


def test_zero_virus_never_infects():
    types_list = [types.MEDIUM] + [types.H] * 50
    virus_at_cell = [0.0] * 51
    rng = np.random.default_rng(1)
    out = transitions.infection_step(types_list, virus_at_cell, g_hv=0.05, rng=rng)
    assert out[1:] == [types.H] * 50


def test_non_h_cells_and_medium_sentinel_are_unchanged():
    types_list = [types.H, types.I, types.D, types.M, types.K, types.E]
    # index 0 is the medium sentinel in the real convention, but here we also
    # check every non-H type is left alone regardless of virus level.
    virus_at_cell = [1000.0] * len(types_list)
    rng = np.random.default_rng(2)
    out = transitions.infection_step(types_list, virus_at_cell, g_hv=1.0, rng=rng)
    assert out[1:] == types_list[1:]


def test_medium_sentinel_at_index_0_is_always_left_unchanged():
    types_list = [types.MEDIUM, types.H]
    virus_at_cell = [1000.0, 1000.0]  # huge virus even at the medium sentinel
    rng = np.random.default_rng(3)
    out = transitions.infection_step(types_list, virus_at_cell, g_hv=1.0, rng=rng)
    assert out[0] == types.MEDIUM


def test_pure_function_does_not_mutate_input():
    types_list = [types.MEDIUM, types.H, types.H]
    virus_at_cell = [0.0, 1000.0, 1000.0]
    rng = np.random.default_rng(4)
    before = list(types_list)
    out = transitions.infection_step(types_list, virus_at_cell, g_hv=1.0, rng=rng)
    assert types_list == before
    assert out is not types_list


def test_same_seed_gives_deterministic_output():
    types_list = [types.MEDIUM] + [types.H] * 20
    virus_at_cell = [0.0] + [5.0] * 20
    out1 = transitions.infection_step(
        types_list, virus_at_cell, g_hv=0.03, rng=np.random.default_rng(42)
    )
    out2 = transitions.infection_step(
        types_list, virus_at_cell, g_hv=0.03, rng=np.random.default_rng(42)
    )
    assert out1 == out2


def test_empirical_infected_fraction_matches_probability_within_tolerance():
    n = 20000
    v_bar = 2.0
    g_hv = 0.05
    types_list = [types.MEDIUM] + [types.H] * n
    virus_at_cell = [0.0] + [v_bar] * n
    rng = np.random.default_rng(7)
    out = transitions.infection_step(types_list, virus_at_cell, g_hv=g_hv, rng=rng)

    infected = sum(1 for t in out[1:] if t == types.I)
    empirical_fraction = infected / n
    expected_fraction = 1.0 - math.exp(-g_hv * v_bar)

    assert abs(empirical_fraction - expected_fraction) < 0.01


def test_high_effective_rate_kills_unresisted_infected_cell():
    types_list = [types.MEDIUM, types.I]
    resist_at_cell = [0.0, 0.0]  # rate = mu_i * (1 - resist) >> 1
    rng = np.random.default_rng(0)
    out = transitions.infected_death_step(types_list, resist_at_cell, mu_i=10.0, rng=rng)
    assert out[1] == types.D


def test_full_resistance_never_kills_regardless_of_mu_i():
    types_list = [types.MEDIUM] + [types.I] * 50
    resist_at_cell = [0.0] + [1.0] * 50  # rate = mu_i * (1 - 1) == 0
    rng = np.random.default_rng(1)
    out = transitions.infected_death_step(types_list, resist_at_cell, mu_i=10.0, rng=rng)
    assert out[1:] == [types.I] * 50


def test_non_infected_cells_and_medium_sentinel_are_unchanged():
    types_list = [types.MEDIUM, types.H, types.I, types.D, types.M, types.K, types.E]
    resist_at_cell = [0.0] * len(types_list)  # rate = mu_i, huge for non-I cells too (ignored)
    rng = np.random.default_rng(2)
    out = transitions.infected_death_step(types_list, resist_at_cell, mu_i=10.0, rng=rng)
    assert out[0] == types.MEDIUM
    assert out[1] == types.H
    assert out[3:] == types_list[3:]


def test_infected_death_pure_function_does_not_mutate_input():
    types_list = [types.MEDIUM, types.I, types.I]
    resist_at_cell = [0.0, 0.0, 0.0]
    rng = np.random.default_rng(4)
    before = list(types_list)
    out = transitions.infected_death_step(types_list, resist_at_cell, mu_i=10.0, rng=rng)
    assert types_list == before
    assert out is not types_list


def test_infected_death_same_seed_gives_deterministic_output():
    types_list = [types.MEDIUM] + [types.I] * 20
    resist_at_cell = [0.0] + [0.3] * 20
    out1 = transitions.infected_death_step(
        types_list, resist_at_cell, mu_i=0.5, rng=np.random.default_rng(42)
    )
    out2 = transitions.infected_death_step(
        types_list, resist_at_cell, mu_i=0.5, rng=np.random.default_rng(42)
    )
    assert out1 == out2


def test_empirical_death_fraction_matches_probability_within_tolerance():
    n = 20000
    resist = 0.25
    mu_i = 0.05
    types_list = [types.MEDIUM] + [types.I] * n
    resist_at_cell = [0.0] + [resist] * n
    rng = np.random.default_rng(7)
    out = transitions.infected_death_step(types_list, resist_at_cell, mu_i=mu_i, rng=rng)

    dead = sum(1 for t in out[1:] if t == types.D)
    empirical_fraction = dead / n
    expected_fraction = 1.0 - math.exp(-mu_i * (1 - resist))

    assert abs(empirical_fraction - expected_fraction) < 0.01
