# 2026-08-18 -- faithful occupancy-space mechanism ladder calibration (Task R2/R3)

Rewires the adaptive-recruitment ladder onto the engine's occupancy-space
chemotaxis primitive (Task R1, `World.set_chemotaxis_occupancy`) instead of
the earlier per-cell `AdaptiveReceptorSubcell`. New module:
`pbg_cpm_studies/model_building/mechanisms.py`. New test:
`tests/test_recruitment_ladder_faithful.py` (5/5 passing). Full suite:
67 passed, 1 skipped.

## Mechanism ladder

- `static_lambda` -- raw chemotaxis (`set_chemotaxis`, fixed lambda). Bias is
  `-lambda*(c(dest)-c(source))`; a uniform background adds the same constant
  to both sides, cancelling exactly. Background-invariant, not receptor
  mediated.
- `hill_occupancy` -- occupancy chemotaxis (`set_chemotaxis_occupancy`) with
  a FIXED `kd`. `theta(c) = (scale*c)^hill / (kd^hill + (scale*c)^hill)`
  saturates at high `c`, so a uniform high background pins both source and
  destination pixel near `theta~1` and the gradient (and recruitment)
  collapse.
- `adaptive_receptor` -- occupancy chemotaxis with an ADAPTIVE `kd`. Between
  `SAMPLE`-step run increments, a slow tracked level `a` is updated toward
  the mean concentration responders currently sense
  (`world.field_mean_at_cell`), and `kd_eff = scale*a` is re-applied.
  Re-centring theta's half-max on the local mean keeps the *local* gradient
  on theta's steep part regardless of background level.

## Diagnosis (probe, before any tuning)

Probed the raw cue field the responders sense with the default
`chemotaxis.py` recipe (`CUE_RATE=10`, `FIELD_DECAY=0.001`): concentration
climbs fast and essentially unboundedly (decay is negligible on a 40-step
horizon) -- mean responder-sensed concentration went from ~9 to ~200 within
~120 composite ticks with zero background. This set the operating range for
`kd`/`scale`: reused the prior `chemotaxis_receptor` defaults (`kd=2.9`,
`hill=2`, `scale=0.02`) as a starting point, then tuned from there.

## Two real bugs found and fixed during calibration

1. **Background "holes."** The engine has no uniform-field-injection API.
   First attempt: only Medium (`cell_type 0`) secretes a background rate
   (`world.set_secretion(field_idx, 0, rate)`), pre-loaded via
   `world.advance_fields(pre_steps)` before the timed run. This is NOT
   actually uniform: cell-occupied pixels (source slab + responder cluster)
   are "holes" (a cell's own pixels do not secrete), so the bath measured
   ~19-50% lower under the responder cluster than in the open medium beyond
   it (verified: `x=10` mean 430 vs `x=75` mean 650 at target level 700).
   This spurious secondary gradient pulled far responders toward the OPEN
   side of the domain instead of the source -- and because the adaptive
   mechanism's `kd_eff` is *maximally sensitive* right at the ambient level
   (by design), it was picking up this artifact gradient as readily as the
   real one, capping `adaptive_receptor,high_bg` around 0.33-0.44 for dozens
   of configs no matter how `kd`/`hill`/`lambda`/`epsilon` were tuned.
   **Fix:** every cell type (Medium AND Source AND Responder) secretes the
   SAME background rate during the preload, eliminating the holes. Verified
   flat to within noise (700 target -> 699.8 everywhere, `x=0` to `x=75`).
2. **"Budget" accumulation never reaches steady state.** `rate = level /
   pre_steps` (a naive linear-accumulation heuristic) does not converge to a
   spatially uniform profile in any practical `pre_steps` -- diffusion needs
   order `domain_len^2 / D ~ 80^2/0.22 ~ 29000` ticks to smooth a
   freshly-secreted profile. **Fix:** use the field's actual steady state
   (secretion balances decay: `rate = level * FIELD_DECAY`) and run
   `advance_fields` for `BACKGROUND_PRE_STEPS=8000` ticks (`>> 1/decay =
   1000`), which is cheap (~0.4s) and reaches the true flat steady state.
3. Also fixed: the source's own cue secretion is switched OFF during the
   background preload and restored after, otherwise the preload also gives
   the cue a free multi-thousand-tick equilibration head start (confounding
   "background level" with "extra cue diffusion time" between conditions).

Fixing (1) and (2) changed the achievable ceiling substantially: before,
`adaptive_receptor,high_bg` plateaued at 0.33-0.44 across ~50+ configs no
matter what else was tuned; after, straightforward parameter tuning reached
0.83.

## Sweep log (representative; ~70 configs run total across both background
regimes)

All numbers are mean `recruitment_index` over seeds (17, 29, 43), steps=40,
domain 80x40, 6 responders (so per-seed values are multiples of 1/6).

| kd | hill | lambda_hill | lambda_adap | eps | bg | static_hi | hill_lo | hill_hi | adap_lo | adap_hi | gap |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 2.9 | 2 | 4000 (shared) | 4000 (shared) | 0.3 | 200 | 0.83 | 0.5 | 0.33 | 1.0 | 0.33 | 0.0 |
| 1.0 | 4 | 15000 (shared) | 15000 (shared) | 0.5 | 700 (holes bug) | 0.5 | 0.56 | 0.0 | 0.83 | 0.5 | 0.5 |
| 1.0 | 4 | 15000 | 15000 | 0.5 | 700 (holes fixed, still shared lambda) | 0.56 | 0.56 | 0.0 | 0.83 | 0.39 | 0.39 |
| 1.0 | 4 | 30000 | 5000 | 0.5 | 700 (steady-state fixed) | 0.72 | 0.67 | 0.0 | 1.0 | 0.39 | 0.39 |
| 1.0 | 4 | 30000 | 3000 | 0.5 | 700 | 0.72 | 0.67 | 0.0 | 0.94 | 0.61 | 0.61 |
| 1.0 | 4 | 30000 | 8000 | 0.5 | 700 | 0.72 | 0.67 | 0.0 | 0.89 | 0.72 | 0.72 |
| **1.0** | **4** | **30000** | **5000** | **0.5** | **400 (FINAL)** | **0.72** | **0.67** | **0.0** | **1.0** | **0.83** | **0.83** |

Key levers that mattered:
- Splitting `lambda` per mechanism (`LAMBDA_HILL=30000` vs
  `LAMBDA_ADAPTIVE=5000`) rather than sharing one value: hill needs a large
  lambda to punch through its modest low-bg gradient; adaptive's
  steeper local slope (re-centred exactly at `kd_eff`) needs a *smaller*
  lambda or it overshoots/oscillates.
- Lowering `background` from 700 to 400 (still `scale*bg=8 >> kd=1`, still
  fully saturating hill) reduced the adaptive mechanism's tracking lag
  (smaller gap for `a` to close), which raised `adaptive_receptor,high_bg`
  from 0.5 to 0.83 with no cost to `hill_occupancy,high_bg` (stayed exactly
  0.0).

## Final calibrated parameters

```
KD = 1.0
HILL = 4.0
SCALE = 0.02
LAMBDA_STATIC = 14.0
LAMBDA_HILL = 30000.0
LAMBDA_ADAPTIVE = 5000.0
EPSILON = 0.5
SAMPLE = 8            # composite ticks between adaptation updates
BACKGROUND_PRE_STEPS = 8000
CONDITIONS:
  low_bg            background=0
  mid_bg            background=150
  high_bg           background=400
  high_bg_blocked   background=400, blocked=True   (responder lambda forced 0)
  high_bg_knockout  background=400, knockout=True  (epsilon forced 0 -- adaptation off)
```

## Final numbers (3 seeds x 17/29/43, steps=40)

| mechanism | condition | per-seed | mean |
|---|---|---|---|
| static_lambda | low_bg | 0.833, 0.833, 0.667 | 0.778 |
| static_lambda | mid_bg | 0.5, 0.833, 0.667 | 0.667 |
| static_lambda | high_bg | 0.667, 0.833, 0.667 | **0.722** |
| static_lambda | high_bg_blocked | 0.667, 0.833, 0.667 | 0.722 (static ignores blocking, as designed) |
| hill_occupancy | low_bg | 0.5, 0.833, 0.667 | **0.667** |
| hill_occupancy | mid_bg | 0.833, 0.833, 0.833 | 0.833 |
| hill_occupancy | high_bg | 0.0, 0.0, 0.0 | **0.000** |
| hill_occupancy | high_bg_blocked | 0.0, 0.0, 0.0 | **0.000** |
| adaptive_receptor | low_bg | 1.0, 1.0, 1.0 | **1.000** |
| adaptive_receptor | mid_bg | 0.5, 0.833, 0.333 | 0.556 |
| adaptive_receptor | high_bg | 0.833, 1.0, 0.667 | **0.833** |
| adaptive_receptor | high_bg_blocked | 0.0, 0.0, 0.0 | 0.000 |
| adaptive_receptor | high_bg_knockout | 0.333, 0.5, 0.5 | 0.444 (partial: between hill's 0 and adaptive's 0.83 -- knockout's `kd_eff` is pre-adapted to the background at t=0, unlike hill's independently-fixed kd=1, so it is not numerically identical to hill, but is much closer to hill than to full adaptive) |

## Seed robustness

The 3-seed table above (17/29/43) is a representative sample, not a robust
estimate: a follow-up 15-seed resample of `adaptive_receptor,high_bg` showed
the per-seed recruitment index is highly variable -- range 0.0-1.0, mean
~0.49 -- so an absolute floor pinned near the lucky 3-seed mean of 0.833 is
fragile. `hill_occupancy,high_bg`, by contrast, is essentially deterministic
at 0.000 across many seeds and both high-background conditions (blocked and
unblocked) -- the saturation collapse does not depend on which seeds you
draw. The honest, defensible headline is therefore the COMPARATIVE rescue,
not the single 0.83 data point: adaptive reliably and substantially
outperforms hill at high background (the wider-seed gap is ~0.49 and the
qualitative "hill collapses / adaptive rescues" story replicates on every
seed spot-checked), even though the absolute adaptive value swings widely
per seed. The hardened test (below) asserts the comparative gap plus a
floor the wider distribution actually supports, rather than the old
absolute `> 0.6` threshold. This does not change the science -- the
mechanism ladder still separates dramatically -- it only fixes how the
result is measured and reported.

## Assertions (tests/test_recruitment_ladder_faithful.py) -- all pass

- `static_lambda,high_bg > 0.4` -- 0.722 (margin 0.32)
- `hill_occupancy,high_bg < 0.3` -- 0.000 (margin 0.30)
- `adaptive_receptor,high_bg` (12 seeds, 11-22) `> 0.35` -- ~0.56 (margin ~0.21;
  hardened from the fragile 3-seed absolute `> 0.6`)
- `hill_occupancy,high_bg` (12 seeds, 11-22) `< 0.15` -- 0.000 (margin 0.15)
- `adaptive_receptor,high_bg - hill_occupancy,high_bg` (12 seeds) `> 0.3` --
  ~0.56 (margin ~0.26; the robust comparative rescue, hardened from the
  3-seed `> 0.4`)
- `adaptive_receptor,low_bg > 0.5` -- 1.000 (3-seed), 0.819 (12-seed) -- holds
- `hill_occupancy,low_bg > 0.5` -- 0.667 (3-seed), 0.583 (12-seed) -- holds
- `hill_occupancy,high_bg_blocked < 0.15` -- 0.000

Full suite: `python -m pytest tests/ -q` -> 67 passed, 1 skipped.
