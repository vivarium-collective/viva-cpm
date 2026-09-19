"""Task 3.3: pure-function tests for per-cell viral resistance.

`cell_resistance(ifn_at_cell, a_rf) = ifn_at_cell / (a_rf + ifn_at_cell)` --
the CC3D-source `resist = f_bar/(a_rf+f_bar)` convention (params.yaml
`resistance:` block; NOT the paper's `1 - f_bar/(theta*a_f+f_bar)` form --
see that block's CONVENTION note for why).
"""
from __future__ import annotations

from pbg_cpm_studies.influenza.resistance import cell_resistance
from pbg_cpm_studies.influenza.params import load_params

A_RF = float(load_params()["resistance"]["a_rf"])


def test_zero_ifn_gives_zero_resistance():
    assert cell_resistance(0.0, A_RF) == 0.0


def test_ifn_equal_to_a_rf_gives_half_resistance():
    assert abs(cell_resistance(A_RF, A_RF) - 0.5) < 1e-12


def test_resistance_monotone_increasing_with_ifn():
    xs = [0.0, A_RF / 10, A_RF / 2, A_RF, A_RF * 2, A_RF * 10, A_RF * 1000]
    vals = [cell_resistance(x, A_RF) for x in xs]
    assert vals == sorted(vals)
    assert vals[0] == 0.0
    assert all(v1 < v2 for v1, v2 in zip(vals, vals[1:]))


def test_resistance_approaches_one_for_large_ifn():
    assert cell_resistance(1000.0, A_RF) > 0.999


def test_resistance_never_negative():
    for x in [0.0, A_RF, A_RF * 5]:
        assert cell_resistance(x, A_RF) >= 0.0


def test_negative_ifn_clamped_to_zero_resistance():
    assert cell_resistance(-1.0, A_RF) == 0.0
