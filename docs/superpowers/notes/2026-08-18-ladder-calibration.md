# 2026-08-18 — Recruitment ladder calibration (Task 4)

Status: **RESOLVED** (empirically calibrated, not blocked). Kept as a note
because the diagnosis and the reasoning behind the final parameter choice
are non-obvious and worth preserving.

## Symptom

At the Task-3 placeholder params (`kd=2.9, hill=2.0, conc_scale=0.02,
activate_occupancy=0.5, epsilon=0.1`, backgrounds `0/80/200`), a 40-step
smoke run gave `recruitment_index == 0.0` at both `low_bg` and
`high_bg_blocked` — no responder ever crossed `activate_occupancy` and
flipped to `ACTIVATED_TYPE`, so nothing ever chemotaxed.

## Diagnosis (probe script, scratchpad only)

Built `recruitment_adaptive(**CONDITIONS["low_bg"])` and
`**CONDITIONS["high_bg"])`, ran 40 steps, and at each step printed, per
responder: raw `world.field_mean_at_cell(0, cid)` (ligand), `theta =
occupancy(ligand + background)`, the running `m`, and `response = theta -
m`.

Findings:

- **`low_bg`** (`background=0`): `theta` climbs slowly as the field
  diffuses out from the source (kd=2.9 needs `ligand*conc_scale ≈ 2.9`, i.e.
  `ligand ≈ 145`, to reach `theta=0.5`). By `t=40` the closest sampled
  responder only reaches `theta≈0.51`, `m≈0.26`, so `response≈0.25` — well
  short of the `0.5` threshold. **Root cause: `theta` itself barely reaches
  the activation threshold within 40 steps, let alone `theta - m`.**

- **`high_bg`** (`background=200`, `epsilon=0.1`): `theta(0) =
  occupancy(200) ≈ 0.655` immediately (background alone, before any ligand
  arrives) — already above `0.5`. But `m` starts at `0.0` (pre-init in
  `cpm.coupling.adaptation_coupling`, not pre-equilibrated to the ambient
  background), so `response(0) = theta(0) - 0 = 0.655` is a one-tick spike
  that decays geometrically (`m` chases `theta` with `epsilon=0.1`, so
  `response(t) ≈ theta(0)·0.9^t` early on). By `t≈5` `response` is already
  `≈0.39`; by `t=40` it has settled to `~0.05–0.08` per responder — again
  well under `0.5`.

**Key insight confirmed**: `response = theta - m` is fundamentally a small,
transient/fold-change quantity (observed range ≈ 0.05–0.3 over 40 steps),
never approaching the `[0,1]` range that raw Hill occupancy `theta` spans.
Reusing `theta`'s calibrated threshold (`0.5`, tuned for
`chemotaxis_receptor`'s raw-occupancy gate) for the adaptive response is a
scale mismatch — the primary bug in the Task-3 placeholder.

## A structural concern raised (and empirically resolved) during calibration

Because `epsilon=0` freezes `m` at its initial value (`0.0`), the
`high_bg_knockout` rung's `response == theta` directly (its own docstring
claim), and `theta ≥ theta - m` pointwise whenever `m ≥ 0` (true for all
`t`, since `m` is an EMA of non-negative values seeded at `0`). At `t=0`
specifically, `response_adaptive(0) == theta_knockout(0)` exactly (both
`background`-only, `m=0` for both branches at the first tick). This raised
a real worry: for *any* shared `activate_occupancy`, does knockout's
activation window (in ticks) always weakly dominate adaptive's — since
knockout's decision variable never decays back down, while adaptive's does
once `m` catches up — making `adaptive_high - hill_high > 0.2` impossible
by construction?

Empirically, no: **the two conditions are separate simulations, not a
shared trajectory**, and activation feeds back into physical position
(activated cells chemotax, changing what ligand they see next). A low
enough shared `activate_occupancy` lets adaptive responders activate
essentially as early as knockout responders (both start from the same
`background`-driven `theta(0)`), and because the adaptive responders then
keep *climbing the real gradient* (local ligand outpaces the slower-tracking
`m`, so `response` stays elevated as they approach the source, generating a
positive feedback loop of its own), the adaptive condition reliably reaches
ceiling recruitment (`1.0`) while `high_bg_knockout` plateaus around
`~0.78` (2 of 12 seed-cell observations never clear the CPM's radius-15
recruitment bar — some randomness/noise floor even though every knockout
responder saturates to `theta≈0.99` and is "activated" essentially the
whole run). This was verified, not assumed — see the sweep log below.

## Sweep log (`.venv/bin/python`, 3 seeds [17,29,43], 40 steps unless noted)

All runs measured `recruitment_index(world, responder_type=ACTIVATED_TYPE)`
(fraction of *activated* responders within radius 15 of the source
centroid) at the final tick, averaged over 3 seeds, for all 4 ladder
conditions (`low_bg`, `high_bg`, `high_bg_blocked`, `high_bg_knockout`).

| # | Config (delta from defaults) | low_bg | high_bg | blocked | knockout | gap (hb−ko) | ALL pass |
|---|---|---|---|---|---|---|---|
| 1 | baseline placeholders (act_occ=0.5) | 0.000 | 0.000 | 0.000 | 0.778 | −0.778 | no |
| 2 | act_occ=0.2 | 0.000 | 0.000 | 0.000 | 0.778 | −0.778 | no |
| 3 | act_occ=0.1 | 1.000 | 0.333 | 0.000 | 0.778 | −0.444 | no |
| 4 | act_occ=0.05 | 1.000 | 0.917 | 0.000 | 0.778 | +0.139 | no (gap<0.2) |
| 5 | act_occ=0.3 | 0.000 | 0.000 | 0.000 | 0.778 | −0.778 | no |
| 6 | act_occ=0.07 | 1.000 | 0.667 | 0.000 | 0.778 | −0.111 | no |
| 7 | act_occ=0.06 | 1.000 | 0.778 | 0.000 | 0.778 | 0.000 | no |
| 8 | act_occ=0.04 | 0.833 | 1.000 | 0.000 | 0.778 | +0.222 | **yes** |
| 9 | act_occ=0.038 | 0.833 | 1.000 | 0.000 | 0.778 | +0.222 | **yes** |
| 10 | act_occ=0.042 | 0.833 | 1.000 | 0.000 | 0.778 | +0.222 | **yes** |
| 11 | act_occ=0.044 | 0.833 | 1.000 | 0.000 | 0.778 | +0.222 | **yes** |
| 12 | act_occ=0.046 | 0.833 | 1.000 | 0.000 | 0.778 | +0.222 | **yes** |
| 13 | act_occ=0.048 | 1.000 | 0.917 | 0.000 | 0.778 | +0.139 | no |
| 14 | act_occ=0.04, steps=60 | 0.500 | 0.000 | 0.000 | 0.778 | −0.778 | no (window too long, response decays back down) |
| 15 | act_occ=0.04, epsilon=0.05 | 0.833 | 0.767 | 0.000 | 0.778 | −0.011 | no |
| 16 | act_occ=0.04, epsilon=0.15 | 0.833 | 0.333 | 0.000 | 0.778 | −0.444 | no |
| 17 | act_occ=0.04, background=150 | 0.833 | 0.917 | 0.000 | 0.778 | +0.139 | no |
| 18 | act_occ=0.04, background=300 | 0.833 | 0.889 | 0.000 | 0.778 | +0.111 | no |

Runs 8–12 identify a stable plateau `activate_occupancy ∈ [0.038, 0.046]`
(at `kd=2.9, epsilon=0.1, background=0/200, steps=40` — the Task-3
placeholders otherwise unchanged) where all three ladder assertions pass
with identical, reproducible per-seed values. `0.04` (run 8, the round
number nearest the plateau center) was chosen as the final default.

Note `high_bg_blocked` is `0.000` in every single row above — the
`receptor_gating` intervention (zeroing the field's chemotaxis list
entirely) trivially satisfies `< 0.2` regardless of `activate_occupancy`,
since blocked responders never move directionally no matter what fate they
carry.

## Final calibrated parameters

`pbg_cpm_studies/composites/chemotaxis_adaptive.py`:

- `ACTIVATE_OCCUPANCY_DEFAULT`: `0.5` → **`0.04`** (the only change; all
  other defaults — `kd=2.9`, `hill=2.0`, `conc_scale=0.02`, `epsilon=0.1`,
  and `CONDITIONS`' background levels `0 / 80 / 200` — are unchanged from
  Task 3).

## Final per-condition numbers (3 seeds [17,29,43], 40 steps, committed params)

| condition | mean recruitment_index | per-seed |
|---|---|---|
| `low_bg` | 0.833 | [1.0, 0.5, 1.0] |
| `high_bg` | 1.000 | [1.0, 1.0, 1.0] |
| `high_bg_blocked` | 0.000 | [0.0, 0.0, 0.0] |
| `high_bg_knockout` | 0.778 | [0.833, 0.833, 0.667] |
| `mid_bg` (sanity, not gated by a test) | 0.867 | [1.0, 1.0, 0.6] |

`tests/test_recruitment_ladder.py` (3 tests) and
`tests/test_chemotaxis_adaptive.py::test_low_bg_recruits_and_blocked_does_not`
all pass with these params; full suite (`tests/`) is 62 passed / 1
pre-existing skip.

## Caveats for whoever cites this in the contract (Task 6)

- The working `activate_occupancy` band is narrow (~0.008 wide out of a
  `[0,1]` domain) and the metric itself is coarse (6 responders ⇒
  increments of 1/6 ≈ 0.167), so the `+0.2` gap requirement is being met by
  a `1.0` vs `0.778` ceiling/plateau contrast, not by a wide margin of
  continuous headroom. Re-running with a different responder count, radius,
  or step count would likely require re-tuning `activate_occupancy` (run 14
  shows `steps=60` alone breaks it: `high_bg`'s response has already decayed
  back below threshold before the run ends, since responders that arrive
  near the source early keep adapting and the transient signal doesn't
  regenerate once nothing is left to climb).
- `high_bg_knockout`'s ~0.78 plateau, not literal near-zero, is what the
  ladder actually rests on — "pure Hill saturates and collapses recruitment"
  is not literally what happens (knockout still recruits most responders,
  since chemotaxis force is a fixed `lambda` once a cell is flagged
  `ACTIVATED`, independent of how *far* `theta` is past threshold). The
  separation is "adaptive reaches the ceiling, knockout plateaus below it,"
  not "knockout collapses to near-zero." Task 6's contract language should
  describe the effect this way rather than as a literal saturation-collapse,
  to avoid overclaiming.
