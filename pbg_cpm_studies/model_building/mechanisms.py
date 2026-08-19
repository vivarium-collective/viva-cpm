"""Occupancy-space chemotaxis mechanism ladder (Task R2/R3).

Three responder-level chemotaxis MODES, applied directly to the Rust engine's
``World.set_chemotaxis`` / ``World.set_chemotaxis_occupancy`` (Task R1) on top
of the plain ``pbg_cpm_studies.composites.chemotaxis`` recruitment world (a
secreting source slab + responders + one diffusing field ``C``). No per-cell
receptor subcell is needed any more -- the mechanism *is* the chemotaxis mode
on the responder type:

  static_lambda     -- raw chemotaxis (``set_chemotaxis``, fixed lambda). The
                        Hamiltonian bias is -lambda*(c(dest)-c(source)); a
                        spatially uniform background adds the SAME constant to
                        both sides of that difference, so it cancels exactly.
                        Static recruits at every background -- it is not
                        receptor-mediated.
  hill_occupancy    -- occupancy chemotaxis (``set_chemotaxis_occupancy``)
                        with a FIXED kd. The bias is
                        -lambda*(theta(dest)-theta(source)) where
                        theta(c) = (scale*c)^hill / (kd^hill + (scale*c)^hill).
                        theta saturates at high c, so a uniform high
                        background pins both dest and source near theta~1 and
                        the gradient collapses -- recruitment fails.
  adaptive_receptor -- occupancy chemotaxis with an ADAPTIVE kd. Between
                        SAMPLE-step run increments, a slow tracked level `a`
                        is updated toward the mean concentration the
                        responders currently sense
                        (``world.field_mean_at_cell``), and kd_eff = scale*a
                        is re-applied. Re-centering theta's half-max on the
                        local mean keeps the *local* gradient on theta's
                        steep part regardless of background level -- gradient
                        sensing (and recruitment) is rescued. epsilon=0
                        disables the tracking update (the knockout condition,
                        ``high_bg_knockout``): kd_eff is initialised once at
                        t=0 and then frozen, so the mechanism degenerates to
                        plain (fixed-kd) Hill occupancy.

The engine has no uniform-field-injection API, so a real, honest background
bath is created the same way any secretion is created: EVERY cell type
(Medium, Source, and Responder alike) is given the same secretion rate for
``BACKGROUND_PRE_STEPS`` field-only ticks (``world.advance_fields`` --
diffusion+secretion, no CPM flips, so cells do not move yet) before the
run begins, then all of it is switched back off except the source's real
cue, which is restored. Secreting from every pixel (not just Medium) matters:
an earlier version only had Medium secrete, and cell-occupied pixels are
"holes" in that scheme (a cell's own pixels do not secrete), so the bath was
measurably lower under the responder cluster and higher in the open medium
beyond it -- a spurious secondary gradient competing with the real source
gradient. The rate is chosen from the field's own steady state
(``rate = level * decay``), and ``BACKGROUND_PRE_STEPS`` (>> 1/decay) is long
enough to actually reach it, so the preloaded bath is verified flat to
within noise. This is a real field value the responders actually sense (not
a bookkeeping fudge), so the adaptive rule's ``field_mean_at_cell`` readout
is honest. See ``_apply_background`` for the full details.

See docs/superpowers/notes/2026-08-18-faithful-ladder-calibration.md for the
empirical sweep that produced the constants below.
"""
from __future__ import annotations

from statistics import mean

import process_bigraph as pb

from pbg_cpm_studies.chemotaxis import metrics as M
from pbg_cpm_studies.composites import chemotaxis as CT

MECHANISMS = ("static_lambda", "hill_occupancy", "adaptive_receptor")

FIELD_IDX = 0
RESPONDER_TYPE = CT.RESPONDER_TYPE
SOURCE_TYPE = CT.SOURCE_TYPE

# --- calibrated constants (docs/superpowers/notes/2026-08-18-faithful-ladder-calibration.md) ---
CUE_RATE = CT.CUE_RATE
LAMBDA_STATIC = 14.0          # raw chemotaxis (matches chemotaxis.py's de-risked recipe)
LAMBDA_HILL = 30000.0         # hill_occupancy's chemotaxis lambda (theta differences are O(0.001-0.01), not O(10-100))
LAMBDA_ADAPTIVE = 5000.0      # adaptive_receptor's chemotaxis lambda (smaller: its slope-at-kd_eff is steeper than hill's tail)
KD = 1.0                      # hill_occupancy's FIXED receptor half-occupancy concentration
HILL = 4.0                    # Hill exponent -- sharper transition helps both adaptive's peak slope and hill's tail collapse
SCALE = 0.02                  # maps raw field units into the Hill regime (matches prior CONC_SCALE_DEFAULT)
EPSILON = 0.5                 # adaptation rate applied every SAMPLE-step increment
SAMPLE = 8                    # composite ticks between adaptation updates (mirrors run_receptor.py's SAMPLE)
BACKGROUND_PRE_STEPS = 8000   # advance_fields ticks to reach the background's true steady state (>> 1/decay=1000)

RECRUITMENT_RADIUS = CT.RECRUITMENT_RADIUS

# Background *steady-state concentration* passed to `_apply_background` (raw field units;
# since every pixel type secretes identically during the preload -- see `_apply_background`
# -- the bath converges to exactly this level everywhere, verified flat to within noise).
# scale*high_bg = 0.02*400 = 8, >> kd=1 (deep saturation, theta'~0); scale*low_bg = 0
# (theta starts at 0, and the source's own growing cue crosses kd's steep part naturally).
CONDITIONS = {
    "low_bg":           {"background": 0.0},
    "mid_bg":            {"background": 150.0},
    "high_bg":           {"background": 400.0},
    "high_bg_blocked":  {"background": 400.0, "blocked": True},
    "high_bg_knockout": {"background": 400.0, "knockout": True},
}


def _responder_ids(spec):
    """CPM cell ids (1-based, in ``cells`` list order) of the responders."""
    return [i + 1 for i, cell in enumerate(spec["cells"]) if cell["type"] == RESPONDER_TYPE]


def _responder_mean_field(world, responder_ids):
    if not responder_ids:
        return 0.0
    return mean(world.field_mean_at_cell(FIELD_IDX, cid) for cid in responder_ids)


def _apply_background(world, level, *, pre_steps=BACKGROUND_PRE_STEPS, cue_rate=CUE_RATE):
    """Pre-load a real, uniform background bath to concentration ``level``
    (raw field units) before the timed run begins.

    Two pitfalls, found empirically (see the calibration note), had to be
    designed around:

    1. If only Medium (type 0) secretes, cell-occupied pixels are "holes" in
       the bath (a cell's own pixels do not secrete), so the bath is LOWER
       wherever the responder cluster sits and HIGHER in the open medium
       beyond it -- a spurious secondary gradient that pulls responders
       toward the open side of the domain, competing with (and for the
       adaptive mechanism, sometimes swamping) the real source gradient.
       Fix: every cell type (Medium AND Source AND Responder) secretes the
       SAME background rate during the preload, so the bath is uniform
       regardless of what is sitting on each pixel.
    2. ``rate = level / pre_steps`` (a naive "budget" accumulation) never
       reaches a spatially uniform state in a physically reasonable number of
       ticks -- diffusion needs order ``domain_len^2 / D`` ticks to smooth a
       freshly-secreted profile, far more than any practical ``pre_steps``.
       Fix: use the field's own steady state (secretion balances decay:
       ``rate = level * decay``) and let ``pre_steps`` (>> 1/decay) run the
       field to that steady state, which -- because every pixel secretes
       identically -- is exactly flat at ``level`` everywhere.

    The source's own cue secretion is switched OFF for the duration of the
    preload and restored afterward -- otherwise the preload would also give
    the cue itself a free ``pre_steps``-tick equilibration head start,
    confounding "background level" with "extra cue equilibration time"
    between conditions."""
    if level <= 0:
        return
    rate = level * CT.FIELD_DECAY
    world.set_secretion(FIELD_IDX, SOURCE_TYPE, rate)
    world.set_secretion(FIELD_IDX, RESPONDER_TYPE, rate)
    world.set_secretion(FIELD_IDX, 0, rate)
    world.advance_fields(pre_steps)
    world.set_secretion(FIELD_IDX, 0, 0.0)
    world.set_secretion(FIELD_IDX, RESPONDER_TYPE, 0.0)
    world.set_secretion(FIELD_IDX, SOURCE_TYPE, cue_rate)


def _run_seed(mechanism, condition, *, seed, steps):
    if mechanism not in MECHANISMS:
        raise ValueError(f"unknown mechanism {mechanism!r}; expected one of {MECHANISMS}")
    cfg = CONDITIONS[condition]
    background = float(cfg.get("background", 0.0))
    blocked = bool(cfg.get("blocked", False))
    knockout = bool(cfg.get("knockout", False))

    core = pb.allocate_core()
    spec = CT.build_spec(cue_rate=CUE_RATE, chemo_lambda=0.0, seed=seed)
    doc = CT.composite_document(cue_rate=CUE_RATE, chemo_lambda=0.0, seed=seed)
    comp = pb.Composite({"state": doc}, core=core)
    world = comp.state["cpm"]["instance"].world
    responder_ids = _responder_ids(spec)

    _apply_background(world, background, pre_steps=BACKGROUND_PRE_STEPS, cue_rate=CUE_RATE)

    # static ignores the receptor-gating intervention (it is not receptor-mediated)
    lam_occ = LAMBDA_ADAPTIVE if mechanism == "adaptive_receptor" else LAMBDA_HILL
    lam = 0.0 if (blocked and mechanism != "static_lambda") else LAMBDA_STATIC if mechanism == "static_lambda" else lam_occ

    epsilon = 0.0 if (mechanism != "adaptive_receptor" or knockout) else EPSILON

    if mechanism == "static_lambda":
        world.set_chemotaxis(FIELD_IDX, RESPONDER_TYPE, lam)
    elif mechanism == "hill_occupancy":
        world.set_chemotaxis_occupancy(FIELD_IDX, RESPONDER_TYPE, lam, KD, HILL, SCALE)
    else:  # adaptive_receptor -- track the local mean concentration, starting from
        # whatever is already sensed at t=0 (a pre-adapted receptor for a fresh cell
        # dropped into this environment), and re-centre kd_eff on it every increment.
        a = _responder_mean_field(world, responder_ids)
        kd_eff = max(SCALE * a, 1e-6)
        world.set_chemotaxis_occupancy(FIELD_IDX, RESPONDER_TYPE, lam, kd_eff, HILL, SCALE)

    n_incr = max(1, steps // SAMPLE)
    for _ in range(n_incr):
        comp.run(SAMPLE)
        world = comp.state["cpm"]["instance"].world
        if mechanism == "adaptive_receptor" and epsilon > 0.0:
            sensed = _responder_mean_field(world, responder_ids)
            a = a + epsilon * (sensed - a)
            kd_eff = max(SCALE * a, 1e-6)
            world.set_chemotaxis_occupancy(FIELD_IDX, RESPONDER_TYPE, lam, kd_eff, HILL, SCALE)

    return M.recruitment_index(world, responder_type=RESPONDER_TYPE, radius=RECRUITMENT_RADIUS)


def simulate_condition(mechanism, condition, *, seeds=(17, 29, 43), steps=40) -> float:
    """Mean recruitment index for ``mechanism`` under ``condition``, averaged
    across ``seeds``. Builds the ``chemotaxis.recruitment`` composite world,
    sets the responder type's chemotaxis per ``mechanism``, runs it in
    SAMPLE-step increments applying the adaptation rule (adaptive_receptor
    only), and reads the recruitment index off the final live world."""
    vals = [_run_seed(mechanism, condition, seed=s, steps=steps) for s in seeds]
    return mean(vals)
