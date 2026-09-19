"""Task 6.1: macrophage-released chemokine + IL-10 fields with regulated
macrophage/uninfected sources -- Increment 6 crux.

Two groups of tests:
  (a) pure-function tests for `signaling.macrophage_secretion_scale` (the
      shared IL-10-Hill self-regulation factor gating macrophage chemokine +
      IL-10 secretion) and `signaling.uninfected_il10_scale` (the uninfected
      cell IL-10 gate, `1 - resist`).
  (b) INTEGRATION: `run.run_macrophage_signaling` (extends Increment 5's
      macrophage scenario with the chemokine + IL-10 fields wired via
      `fields.add_chemokine_field`/`add_il10_field` + this task's per-cell
      regulation) -- the chemokine field is positive and HIGHEST near the
      macrophage cluster, decaying with distance (a gradient centered on the
      macrophages, per the raw-lattice radial profile AND the
      macrophage-vs-uninfected-cell per-cell readings); IL-10 is positive.

Behavior (gradient shape, IL-10 positivity) is the fidelity criterion here,
not exact magnitudes -- this is mechanism validation (Increment 6), not a
Fig-3/6/7 reproduction (Increment 9); the sig_1 STUB (see `signaling.py`'s
module docstring / `params.yaml`'s `il10.sig_1_stub`) means absolute
secretion magnitudes are not meant to be load-bearing here.
"""
from __future__ import annotations

from pbg_cpm_studies.influenza import fields, run, signaling
from pbg_cpm_studies.influenza.params import load_params

_PARAMS = load_params()
SIG_1, G_1, G_2, D_2 = fields.il10_hill_constants()


# --- (a) pure functions -----------------------------------------------

def test_hill_scale_positive_at_zero_local_il10():
    scale = signaling.macrophage_secretion_scale(0.0, SIG_1, G_1, G_2, D_2)
    assert scale > 0.0


def test_hill_scale_bounded_between_zero_and_one():
    for L in (0.0, 1e-6, 1e-3, 1.0, 1000.0):
        scale = signaling.macrophage_secretion_scale(L, SIG_1, G_1, G_2, D_2)
        assert 0.0 <= scale <= 1.0, f"L={L}: scale={scale} out of [0, 1]"


def test_hill_scale_monotone_decreasing_in_local_il10():
    # g_1, d_2 > 0 -> the ratio term (g_1*L+g_2)/(L+d_2) increases with L ->
    # the Hill factor sig_1/(sig_1+ratio) decreases with L (more local
    # IL-10 -> more self-limiting -> lower secretion scale).
    Ls = [0.0, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10.0]
    scales = [signaling.macrophage_secretion_scale(L, SIG_1, G_1, G_2, D_2) for L in Ls]
    assert all(s1 >= s2 for s1, s2 in zip(scales, scales[1:])), (
        f"expected non-increasing Hill scale with local IL-10: {list(zip(Ls, scales))}"
    )
    assert scales[0] > scales[-1], "Hill scale should genuinely vary, not be flat"


def test_hill_scale_negative_local_il10_clamped_like_zero():
    assert signaling.macrophage_secretion_scale(-1.0, SIG_1, G_1, G_2, D_2) == \
        signaling.macrophage_secretion_scale(0.0, SIG_1, G_1, G_2, D_2)


def test_hill_scale_nonpositive_sig_1_returns_zero_not_negative_or_nan():
    assert signaling.macrophage_secretion_scale(0.0, 0.0, G_1, G_2, D_2) == 0.0
    assert signaling.macrophage_secretion_scale(1.0, -1.0, G_1, G_2, D_2) == 0.0


def test_uninfected_il10_scale_is_one_minus_resist():
    for resist in (0.0, 0.25, 0.5, 0.9):
        assert abs(signaling.uninfected_il10_scale(resist) - (1.0 - resist)) < 1e-12


def test_uninfected_il10_scale_bounded_zero_one():
    for resist in (-1.0, 0.0, 0.5, 1.0, 2.0):
        scale = signaling.uninfected_il10_scale(resist)
        assert 0.0 <= scale <= 1.0, f"resist={resist}: scale={scale} out of [0, 1]"


def test_uninfected_il10_scale_full_at_zero_resist():
    assert signaling.uninfected_il10_scale(0.0) == 1.0


def test_uninfected_il10_scale_zero_at_full_resist():
    assert signaling.uninfected_il10_scale(1.0) == 0.0


# --- (b) integration -----------------------------------------------------

def test_chemokine_and_il10_fields_are_positive_after_macrophage_signaling_run():
    result = run.run_macrophage_signaling(steps=15, seed=17)

    assert result["total_chemokine"][-1] > 0.0, "chemokine field should be positive"
    assert result["total_il10"][-1] > 0.0, "IL-10 field should be positive"
    # Monotonically non-decreasing lattice-wide sums over the run (sources
    # only add; decay/diffusion redistribute but the driver's field-warmup
    # is 0 by default so there's no pre-loop head start to decay from).
    assert result["total_chemokine"][-1] > result["total_chemokine"][0]
    assert result["total_il10"][-1] > result["total_il10"][0]


def test_chemokine_field_is_highest_near_macrophages_and_decays_with_distance():
    result = run.run_macrophage_signaling(steps=15, seed=17, radial_bins=6)

    # Per-cell reading: macrophages (the source) see MORE chemokine than the
    # uninfected epithelial cells, `separation_sites` of open Medium away.
    chemo_macro = result["chemo_at_macrophages"][-1]
    chemo_far = result["chemo_at_uninfected"][-1]
    assert chemo_macro > 0.0
    assert chemo_macro > chemo_far, (
        f"chemokine at macrophages ({chemo_macro:.6f}) should exceed chemokine "
        f"at the far epithelial patch ({chemo_far:.6f})"
    )

    # Raw-lattice radial profile around the macrophage cluster's centroid,
    # independent of any cell's position: a genuine gradient CENTERED on the
    # macrophages should decay (non-increasing, with a real drop end-to-end)
    # from the innermost bin outward.
    profile = result["chemo_radial_profile"]
    means = [m for m in profile["bin_means"] if m == m]  # drop NaN (empty bins)
    assert len(means) >= 2, f"radial profile too sparse to assess a gradient: {profile}"
    assert means[0] > 0.0, "innermost (macrophage-centered) bin should be positive"
    assert means[0] > means[-1], (
        f"chemokine should decay from the macrophage-centered bin outward: {means}"
    )
    non_increasing_pairs = sum(1 for a, b in zip(means, means[1:]) if a >= b)
    assert non_increasing_pairs >= len(means) - 2, (
        f"radial profile should be predominantly non-increasing with distance "
        f"(allowing at most one local bump from noise/geometry): {means}"
    )


def test_il10_positive_from_macrophage_and_uninfected_sources():
    result = run.run_macrophage_signaling(steps=15, seed=17)

    assert result["il10_at_macrophages"][-1] > 0.0, "macrophage-sourced IL-10 should be positive"
    assert result["il10_at_uninfected"][-1] > 0.0, "uninfected-cell-sourced IL-10 should be positive"


def test_macrophage_signaling_run_is_deterministic():
    a = run.run_macrophage_signaling(steps=8, seed=17)
    b = run.run_macrophage_signaling(steps=8, seed=17)
    assert a["total_chemokine"] == b["total_chemokine"]
    assert a["total_il10"] == b["total_il10"]
    assert a["chemo_radial_profile"] == b["chemo_radial_profile"]
