# Recruitment model-building — filling the studies + an agentic contract→implementation loop

**Status:** design (approved in brainstorming 2026-08-18)
**Workspace:** viva-cpm (`cpm-studies`)
**Author:** Eran Agmon (with Claude)

## 1. Motivation & goals

The `chemotactic-recruitment` investigation is a richly-authored but *static* record:
five studies that each pass, documenting one biological claim ("a secreted cue
recruits responder cells") at two levels of realization — a phenomenological
chemotaxis-λ knob and a receptor-level Hill-occupancy model calibrated to a cited
Kd. Two things are missing:

1. **The studies are not reviewer-complete.** They fill the Model/Tests/
   Simulations/Results tabs but omit the **v4 narrative-spine** sections, carry
   **thin readouts**, under-document the **simulations**, and the receptor arm
   lacks **controls / alternative hypotheses**.
2. **The investigation hides its own structure.** The three realizations are
   secretly a **mechanism-refinement ladder** — `static-λ` → `hill-occupancy` →
   (a new) *adaptive* rung — but nothing makes that ladder explicit or shows how
   a model-builder would *climb* it from evidence.

This project delivers **both**:

- **Part A** — fill all five existing studies to the canonical study-detail
  spine standard, so every tab shows and the report linter's narrative-spine
  check is satisfied.
- **Part B** — a new sibling investigation, `recruitment-model-building`, that
  runs the **agentic model-building loop** (viva-superpowers'
  contract→audit→lock→build→evaluate→navigate pipeline, the one demonstrated in
  `viva-casebook`) to build a **new, harder refinement** — *adaptive recruitment*
  (fold-change detection) — from a locked contract of tests, capturing a real,
  emergent trajectory and rendering it as a flagship pipeline report.

### Non-goals (YAGNI)

- **No live-LLM agent harness / benchmark.** The loop driver is a deterministic
  NAVIGATE policy (as in casebook's `build_loop_demo.py`). Study fields are
  authored so a live harness could be added later without rework, but it is out
  of scope here.
- **No reimplementation of the loop framework.** We reuse the installed
  `viva_superpowers.{loop_state, test_audit, test_contract, study_evaluator,
  study_verdict, module_sourcing}` — no forks.
- **No new Rust.** The adaptation mechanism is a Python-level extension of the
  existing `ReceptorSubcell` coupling; the CPM Rust core is unchanged.

## 2. The unifying idea — an emergent mechanism ladder

The existing studies are the lower rungs of a ladder of chemotactic-response
mechanisms. Crucially, the ladder is **monotone in tests-passed** only because
the contract rewards receptor realism *and* background-invariance at once — so a
cruder mechanism always leaves exactly one hard test to fail:

| Rung | Mechanism | Response signal | Passes | Fails |
|---|---|---|---|---|
| 1 | `static_lambda` | constant λ up the **raw local gradient ∇c** (not receptor-gated) | `recruits_low`, `recruits_high` (∇c is offset-invariant) | **`receptor_gating`** — blocking the receptor does nothing (there is no receptor to block) |
| 2 | `hill_occupancy` | λ ∝ Hill occupancy θ(c); gradient sensed as Δθ across the cell | `recruits_low`, `receptor_gating` (block → abolished) | **`recruits_high`** — θ saturates (→1) at high background, so Δθ → 0 and gradient sensing dies |
| 3 | `adaptive_receptor` *(new)* | λ ∝ (θ − m), slow adaptation dm/dt = ε(θ − m) | `recruits_low`, `receptor_gating`, `recruits_high` (fold-change detection) | — (knockout m=0 reverts to `hill_occupancy` → fails `recruits_high`: the discriminating control) |

Why each mechanism fails exactly one rung:
- `static_lambda` follows the **raw** cue gradient, which is unchanged by adding a
  uniform background offset — so it recruits at any background, but it is *not*
  receptor-mediated, so "block the receptor" cannot abolish recruitment. It fails
  `receptor_gating`.
- `hill_occupancy` *is* receptor-mediated (blocking abolishes recruitment), but it
  senses the gradient through **occupancy**, which saturates at high background —
  every part of the cell reads θ≈1, the intracellular contrast vanishes, and the
  cell stops climbing. It fails `recruits_high`.
- `adaptive_receptor` keeps receptor gating *and* restores background-invariance:
  a slow setpoint `m` adapts away the mean occupancy so the cell responds to the
  *local deviation* (θ − m) — fold-change detection. It passes all three; its
  adaptation-knockout (`m≡0`) collapses back to `hill_occupancy`.

Part B makes this ladder **emergent**: the loop starts from `static_lambda` and,
driven only by the locked contract, installs the mechanism the top failing hard
test pulls on — `receptor_gating` fails → install `hill_occupancy`; then
`recruits_high` fails → install `adaptive_receptor`; then all pass → DONE. The
climb `static → hill → adaptive` is a *consequence* of the policy, not a script.
The adaptation rung is genuinely new — it does not exist in the codebase today —
so the contract→implementation arc is real, not a replay. Note `receptor_gating`
is exactly the intervention the existing `recruitment-receptor-blocked` study
already documents, tying the loop to the investigation it recapitulates.

The science is standard bacterial-chemotaxis adaptation: a slow internal setpoint
`m` (methylation-like) tracks occupancy so the *steady-state* response returns to
baseline regardless of absolute cue level, while the transient/gradient response
survives — **fold-change detection** (Barkai–Leibler perfect adaptation; Tu et al.
fold-change detection). In the CPM this is an extra per-cell ODE state in the
responder's receptor coupling; the chemotaxis λ scales with r = max(0, θ − m).

## 3. Part A — fill the five existing studies

Target studies (members of `chemotactic-recruitment`):
`recruitment-baseline`, `recruitment-inhibited`, `recruitment-adversarial`,
`recruitment-receptor-baseline`, `recruitment-receptor-blocked`.

For each, bring every study-detail tab to the reviewer-complete standard:

### A1. Model tab — *show* the mechanism
- Keep the canonical `conditions.baseline` (composite + params).
- Add/extend `conditions.model_settings[]` documenting **each parameter**: name,
  value, units, role, and provenance (e.g. `chemo_lambda` = phenomenological
  copy-attempt bias; `kd = 2.9 nM` cited Nasser 2009; `hill`, `conc_scale`,
  `activate_occupancy`). The Model tab then reads as a described mechanism, not a
  bare param map.

### A2. Readouts — the full observable set
Expand the thin single-readout studies to the canonical set (each with
identifier/units/status/description):
- All: `recruitment_index`, `mean_approach`, `responder_mean_distance`.
- Receptor arm additionally: `receptor_occupancy` (θ) and `activation_fraction`
  (fraction of responders that crossed the activation threshold).
- Validate against the composite's real emitted paths (`/viva-study
  check-observables` grammar) so the readouts are canonical, not prose.

### A3. Simulations tab — detail the runs
Enrich each `runs[]` entry: `n_steps`, `seeds`, the 600-sweep `advance_fields`
field equilibration, emitter, engine settings, and `provenance`. Where a study
declares completion, ensure `runs[].outcomes` are present (they are) and fully
described so the "Ran · Tests · Verdict" strip is unambiguous.

### A4. Controls + alternative hypotheses — receptor arm
Add to `recruitment-receptor-baseline` and `recruitment-receptor-blocked`:
- `controls[]`: `recruitment-receptor-blocked` is the negative-control
  intervention (block response downstream of receptor activation → recruitment
  abolished); cross-link the adversarial (no-cue) arm.
- `alternative_hypotheses[]`: e.g. "accumulation by adhesion differences, not
  receptor-gated chemotaxis" — excluded by the unconfounded contact-J design;
  "metric artifact" — excluded by the blocked arm.

### A5. v4 narrative spine — all five
Author the missing sections the linter flags (the load-bearing ★ trio plus the
rest of the 15): `report`, `study_card`, `conclusion_verdicts` (3-track:
computational / biological / methodological verdicts), `biological_summary`,
`literature_anchors` (Nasser 2009 Kd; Glazier–Graner; chemotaxis refs), and the
remaining spine fields (`runtime`, `enforced_params`, `model_change`,
`implementation_requirements`, `design_pivot_required`) as applicable. Exact
field shapes are read from `viva_superpowers` (the narrative-spine linter /
`study_narrative`) at implementation time; the section list above is canonical.

### A6. `purpose.mechanism` tokens
Add a `purpose` block with `mechanism:` tokens (e.g. `chemotactic_recruitment`,
`receptor_occupancy`) so each study is loop-aware and satisfies the audit's
objective-coverage check.

**A verification:** `study_audit --gate` stays green (0 hard failures); the
report linter reports the narrative spine **complete** for all five; each study's
readouts pass `check-observables`.

## 4. Part B — the `recruitment-model-building` investigation

### B1. Mechanism library — real Processes
New module `pbg_cpm_studies/model_building/mechanisms.py` defining the ladder as
**real `process_bigraph` response-couplings** over the CPM engine, extending the
existing `cpm.coupling.receptor_coupling` / `ReceptorSubcell`:

- `static_lambda` — constant chemotaxis bias (the phenomenological baseline).
- `hill_occupancy` — instantaneous Hill occupancy θ (the current receptor model).
- `adaptive_receptor` — occupancy + a slow per-cell adaptation state `m`
  (dm/dt = ε(θ − m)); response signal r = max(0, θ − m); λ scales with r.
  Fold-change detection.

"Installing a mechanism" = swapping the responder's response coupling in the
composite document — a genuine model edit the loop performs. A `LIBRARY` dict maps
each mechanism → its default knobs + literature citation (mirroring casebook's
`LIBRARY`), with **no test→mechanism answer key** baked in (NAVIGATE reasons from
the failing test's `knob`/mechanism token).

`pbg_cpm_studies/model_building/navigate.py` — the deterministic NAVIGATE policy:
given the graded contract, take the most-negative **hard** test and install the
mechanism it pulls on (or take one knob-calibration step). Emergent: the climb
static→hill→adaptive is a *consequence* of this policy, never a scripted list.

### B2. The contract — `recruitment-adaptive` study
New study `workspace/studies/recruitment-adaptive/study.yaml`, schema v4, authored
to participate in the loop:

```yaml
question: |
  Does the responder recruit across a RANGE of background cue levels
  (fold-change detection), not just at one background?
purpose:
  mechanism: adaptation            # + fold_change_detection
requires: [spatial_chemotaxis, receptor_occupancy, adaptation]
sourcing:
  decision: compose                # reuse chemotaxis_receptor + build-new adaptation
  modules: [chemotaxis_receptor]
  rationale: >-
    receptor occupancy is catalogued; no module provides adaptation, so compose
    the adaptation feedback on top of the receptor coupling.
catalog:
  chemotaxis: [spatial_chemotaxis]
  chemotaxis_receptor: [spatial_chemotaxis, receptor_occupancy]
behavior_tests:                    # THE CONTRACT (locked) — 3 primary + 1 control
  - name: recruits_low                        # sanity: all rungs pass
    classification: primary
    measure: {kind: recruitment_index, condition: low_bg, stat: final}
    pass_if: {op: in_range, low: 0.5, high: 0.9}
    cites: [recruitment-baseline]
  - name: receptor_gating                     # static FAILS (no receptor); hill/adaptive pass
    classification: primary
    measure: {kind: recruitment_index, condition: high_bg_blocked, stat: final}
    pass_if: {op: "<=", value: 0.15, provenance: {kind: cited}}
    cites: [Nasser2009]
  - name: recruits_high                       # hill FAILS (θ saturates); adaptive passes
    classification: primary
    measure: {kind: recruitment_index, condition: high_bg, stat: final}
    pass_if: {op: in_range, low: 0.5, high: 0.9}
    cites: [Barkai1997, Tu2008]
  - name: adaptation_knockout_collapses       # discriminating control for the top rung
    classification: diagnostic
    control: negative
    measure: {kind: recruitment_index, condition: high_bg_knockout, stat: final}
    pass_if: {op: "<=", value: 0.15, provenance: {kind: first_principles}}
conditions:
  baseline: {composite: chemotaxis_adaptive.recruitment_adaptive, params: {kd: 2.9, epsilon: <tuned>}}
  variants: []
  model_settings: [...]            # documented adaptation knobs (ε, setpoint, kd cite)
controls:
  - {name: adaptation_knockout, kind: negative, description: "clamp m≡0 → reverts to hill_occupancy → fails recruits_high"}
readouts: [recruitment_index, adaptation_ratio, receptor_occupancy, activation_fraction]
```

`receptor_gating` is the receptor-blocked condition (recruitment must be *abolished*
when the response is blocked downstream of binding) — the test `static_lambda`
cannot pass and the driver of the static→hill step. `recruits_high` is the
fold-change test only `adaptive_receptor` passes. The **draft** carries the
one-sided form of `recruits_high` (e.g. `op: ">="`) so **AUDIT round 1** flags
discrimination drift; round 2 revises to the two-sided band above (cited). This
demonstrates the two-round audit genuinely rejecting an insufficient draft.
(`adaptation_ratio` is retained as a **readout** and a directional axis in the
report, not a hard gate, so an uncalibrated ε is not falsely failed.)

A supporting composite (`pbg_cpm_studies/composites/chemotaxis_adaptive.py`,
`@composite_generator` id `chemotaxis_adaptive.recruitment_adaptive`) exposes the
adaptive coupling under the background conditions the contract measures:
`low_bg`, `mid_bg`, `high_bg` (uniform baseline cue added on top of the source
gradient), `high_bg_blocked` (high background with the response blocked
downstream of binding — the `receptor_gating` condition), and `high_bg_knockout`
(high background with adaptation clamped, `m≡0` — the discriminating control).
The same coupling, run under different conditions; the loop grades each.

### B3. The loop driver — emergent trajectory
`scripts/build_recruitment_loop.py` (mirroring `viva-casebook/scripts/
build_loop_demo.py`) drives `recruitment-adaptive` through the full loop:

**contract → AUDIT (2 rounds) → SELECT/sourcing → prereg-LOCK → emergent
BUILD→RUN→EVALUATE→DECIDE→NAVIGATE → DONE**, using
`viva_superpowers.{loop_state, test_audit, test_contract, study_evaluator}`.

- Each iteration assembles a **real `process_bigraph.Composite`** from the
  currently-installed mechanism and runs it on a **small lattice** (for
  minutes-scale regeneration) at LOW/MID/HIGH background, grading
  `recruitment_index` and `adaptation_ratio` against the locked bands via
  `test_contract.check`.
- NAVIGATE (B1) installs `hill_occupancy` then `adaptive_receptor` as the
  high-background / adaptation tests demand → reaches **DONE** at ~iteration 3.
- **Honest GIVE_UP companion**: same loop with `adaptive_receptor` withheld from
  the library → climbs to `hill_occupancy`, improves low-bg but stalls on
  high-bg/adaptation → budget spent → **GIVE_UP** (OPEN), proving the loop
  terminates honestly instead of faking a pass.
- **Perturbation robustness**: re-grade the final model under dt/2 and ±10% knob
  perturbations.
- Writes `.pbg/loop/recruitment-adaptive.json` (real `loop_state`) +
  `workspace/investigations/recruitment-model-building/trajectory.json`
  (schema **`model_build_trajectory/v2`**: `contract`, `draft`,
  `final_composite`, `select`, `audit{round1_gate,round2_gate,axes,round1_flags}`,
  `control`, `lock{tests_hash,n_tests_locked,reopen_count,prior_hashes}`, `tests`,
  `iterations`, `result`, `giveup_companion`, `timeseries`).

### B4. Flagship report
`scripts/render_recruitment_pipeline_report.py` →
`docs/recruitment-pipeline-report.html` (viva-cpm styled; mirrors casebook's
`render_pipeline_report.py`). Reads `trajectory.json` and renders, in order:
**contract → draft → select → tests → sufficiency audit (2 rounds) → lock → the
EMERGENT build loop (a per-iteration margin matrix + NAVIGATE decisions,
illustrated with the spatial recruitment scenes at low vs high background) →
result → honest GIVE_UP companion → pipeline health.** Every value read from the
real captured trajectory; re-run the driver then the renderer to refresh.

### B5. The investigation YAML
`workspace/investigations/recruitment-model-building/investigation.yaml`, fully
populated in the style of casebook's `model-sourcing` investigation:
`question`, `hypothesis`, `lead`, `executive` (what_is_this / verdict /
verdict_status / decisions_needed), `at_a_glance` (per-rung role table), a
capability `catalog` (module→capabilities), `members: [recruitment-adaptive]`,
`acceptance_criteria`, `expert_docs`. Cross-links `chemotactic-recruitment` as the
"how these models were built" companion (and names its static/hill studies as the
lower rungs the loop recapitulates).

## 5. File plan

New:
- `pbg_cpm_studies/model_building/__init__.py`
- `pbg_cpm_studies/model_building/mechanisms.py` — the 3-rung library (real couplings)
- `pbg_cpm_studies/model_building/navigate.py` — deterministic NAVIGATE policy
- `pbg_cpm_studies/composites/chemotaxis_adaptive.py` — adaptive composite + background conditions
- `workspace/studies/recruitment-adaptive/study.yaml` — the loop's target/contract
- `workspace/investigations/recruitment-model-building/investigation.yaml`
- `scripts/build_recruitment_loop.py` — the loop driver (emergent trajectory)
- `scripts/render_recruitment_pipeline_report.py` — flagship renderer
- `docs/recruitment-pipeline-report.html` — generated flagship (committed)
- `tests/test_recruitment_mechanisms.py`, `tests/test_recruitment_loop.py`

Modified:
- the 5 existing `workspace/studies/recruitment-*/study.yaml` (Part A)
- `workspace/references/papers.bib` — add Barkai1997, Tu2008 (adaptation refs)

## 6. Testing & verification

- **Mechanism unit tests** — `static_lambda` recruits at baseline but fails
  selectivity; `hill_occupancy` passes low-bg but saturates/fails high-bg;
  `adaptive_receptor` passes across LOW/MID/HIGH; adaptation-knockout collapses
  high-bg recruitment.
- **Loop-driver tests** — reaches DONE at `adaptive_receptor`; the GIVE_UP
  companion terminates OPEN; AUDIT round 1 flags the one-sided draft and round 2
  passes; `loop_state.validate` reports no invariant violations.
- **Study/report health** — `study_audit --gate` green across the (now 17)
  studies; report-linter narrative spine complete for all five filled studies;
  the new study and investigation pass structural lint.
- **CI** — `pytest tests/` green in a fresh venv; the loop-driver test bounded to
  a small lattice so CI stays fast.

## 7. Phasing (each independently shippable)

1. **Phase 1 — Part A**: fill the five studies to the spine standard. Reviewer-
   ready; no new Python. Ship as its own PR.
2. **Phase 2 — B1–B3**: mechanism library + adaptive composite + loop driver +
   `recruitment-adaptive` study + captured `trajectory.json`. Ship as a PR.
3. **Phase 3 — B4–B5**: flagship report + investigation YAML + cross-links +
   republish the read-only workbench. Ship as a PR.

## 8. Risks & mitigations

- **Adaptation tuning** — fold-change detection needs the adaptation timescale ε
  slow relative to cell motion but fast relative to the run. Mitigation: expose ε
  as a `model_settings` knob; the loop's CALIBRATE branch can take a knob step,
  and the unit test brackets the working range.
- **Loop runtime in CI** — full recruitment lattices are slow. Mitigation: the
  driver uses a reduced lattice / step count; the flagship report is regenerated
  offline and committed, not built in CI.
- **v4 spine schema drift** — exact field shapes live in `viva_superpowers`.
  Mitigation: read the current schema at implementation time; treat the linter's
  "complete" result as the acceptance signal, not a hand-copied schema.
- **Emergence honesty** — the trajectory must be a real consequence of the
  policy, not authored. Mitigation: the driver captures grades from real Composite
  runs (as casebook does); the GIVE_UP companion + `loop_state.validate` guard
  against a faked pass.
