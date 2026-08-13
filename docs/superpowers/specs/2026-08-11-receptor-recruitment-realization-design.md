# Receptor-level realization of chemotactic recruitment

*Design spec — 2026-08-11. Workspace: `viva-cpm`. Investigation: `chemotactic-recruitment`.*

## 0. Context and motivation

This is the first **use-case-driven** increment of the "typed translator / citable-evidence"
program (see `~/ai-generated/type-guarantees-formalization.pdf` / the Domain Bridges design). A **translator** is a typed, law-bearing artifact that converts a value from
one representation/vocabulary to another — an *adapter/translator with a contract*; the Domain Bridges
design calls it a "domain bridge", and it is **not** the process-bigraph Composite bridge (structural
wiring within one vocabulary). The investigation is not merely a consumer of that formalism — it is the **forcing function to harden the
underlying substrate** (bigraph-schema / process-bigraph / viva) in the translator direction and to
**streamline** it. The discipline (§4): build only what this investigation exercises, land each piece as
real, tested framework at the right layer, and **prune the ad-hoc glue it replaces** — prune as much as
we add.

The investigation **names its own refinement and records the exact gap we close**:

- Composite already exposes the down-scale socket, wired but empty:
  `"inputs": {"fates": ["fates"]}` with `"fates": {}` (`pbg_cpm_studies/composites/chemotaxis.py`).
- Recorded semantic gap (investigation.yaml `caveats`): *"CPM chemotaxis lambda is a PHENOMENOLOGICAL
  response strength … not a receptor-level model — a recorded semantic gap, not a hidden one."*
- Named refinement (investigation.yaml `competing_frameworks`): *"Receptor-kinetics / ODE chemotaxis
  … the CPM lambda is its phenomenological coarse-graining. A claim-layer refinement contract would
  state what a receptor-level realization must preserve."*
- Open hole (claim-bundle/realization.yaml): `hole:lambda-calibration`, `validated: null`, an uncited
  `biologically_plausible: [5,20]` for λ. `papers.bib` currently has one entry (`glazier1993`, which
  contains no chemotaxis) and **zero** structured `cites`/`calibration_anchor` anywhere in the workspace.

This spec adds a **receptor-level realization** of the *same* biological claim: a per-cell
ligand–receptor occupancy model, calibrated to a real chemokine–receptor affinity, drives an
**emergent** chemotactic response — resolving the recorded gap and filling the open hole with citable
evidence — while hardening and streamlining the framework underneath it (§4).

## 1. Biological question and the insight it sharpens

Current investigation answers **on/off**: cue AND response are both necessary, with the response a
phenomenological knob (λ=14). The receptor realization answers the sharper, quantitative question:

> **Does receptor-level signaling, calibrated to a real chemokine–receptor Kd, reproduce recruitment —
> and what does making the response emergent reveal that a scalar λ hides?**

Three things it reveals that λ cannot:

1. **A concentration threshold for recruitment** set by the receptor Kd (nanomolar, cited) — a
   falsifiable quantitative prediction instead of an arbitrary λ.
2. **Spatial heterogeneity of activation** — responders nearer the source cross the occupancy threshold
   first, producing an activation front; a genuinely new readout.
3. **A mechanistic dose-response** — recruitment vs chemokine dose follows the receptor-occupancy curve,
   directly comparable to real chemotaxis dose-response assays.

## 2. Architecture (zero Rust changes)

The Rust engine stores chemotaxis λ **per cell-type** (`crates/cpm-core/src/field.rs:13`
`chemotaxis: Vec<f64>  // per cell_type … index = cell_type`; `delta_chemotaxis` looks up λ by the
moving cell's *type*). The `fates` port sets a cell's *type* (`set_cell_type`). So emergent
responsiveness couples back by switching a cell's **type**, and each responder sub-type carries its own
λ in the field spec. No engine changes are required.

### 2.1 `ReceptorSubcell` (new Process) — `cpm/subcellular/receptor.py`

A minimal per-cell receptor-occupancy model, ~40 lines, in the style of the existing
`cpm/subcellular/boolean.py` / `sbml.py`.

- **inputs:** `{"ligand": "float"}` — the local chemokine concentration, wired from the CPM output
  `field_at_cell[str(cid)]`.
- **outputs:** `{"fate": "overwrite[integer]"}` — the responder sub-type, wired to `fates[str(cid)]`.
- **config_schema:** `kd` (float; chemokine–receptor dissociation constant, cited), `hill` (float,
  default 1.0), `conc_scale` (float; maps CPM field units → concentration units, default 1.0),
  `activate_occupancy` (float in (0,1]; occupancy threshold to switch naive→activated, default 0.5),
  `naive_type` (int, default 2), `activated_type` (int, default 3).
- **update:** occupancy `θ = c^h / (Kd^h + c^h)` with `c = ligand * conc_scale`; emit
  `activated_type` if `θ ≥ activate_occupancy` else `naive_type`. (Deterministic; adaptation deferred —
  §7.) `activate_occupancy=0.5` makes the recruitment threshold coincide with `c = Kd`, so the cited Kd
  is directly visible in the dose-response.

*Rationale for a new file rather than `SBMLSubcell`:* a two-parameter Hill occupancy is the "minimal
cited ODE" chosen in brainstorming; it needs no SBML/tellurium dependency and keeps the cited Kd a
first-class config field an evidence bridge can target. `SBMLSubcell` remains the path for a fuller
published model later (§7).

### 2.2 Responder sub-types

- **naive responder** = type 2, chemotaxis λ = 0.
- **activated responder** = type 3, chemotaxis λ = 14 (the value the phenomenological realization uses).
- **Contact energies J for type 3 are set identical to type 2** so that *only* chemotaxis differs
  between naive and activated — no confound with differential adhesion. This is a required invariant of
  the composite spec (asserted in a test, §6).
- All responders initialize as **type 2 (naive)**; `ReceptorSubcell` flips a cell to type 3 when local
  occupancy crosses threshold. (Graded naive/partial/full with intermediate λ is a trivial extension for
  a smoother dose-response — §7.)

### 2.3 Composite — `pbg_cpm_studies/composites/chemotaxis_receptor.py`

New `@composite_generator` factory `recruitment_receptor`, reusing `chemotaxis.build_spec(...)` for the
CPM world (adding the type-3 contact rows) and adding, per responder cell:

- a `ReceptorSubcell` node,
- input wire `[field_at_cell, str(cid)]`,
- output wire `[fates, str(cid)]`,
- a pre-initialized `fates` store with a key per wired responder cell — required **until H5 lands**,
  because a write to an absent map key is currently dropped (`methods/apply.py:438-439`; see
  `cpm/processes/cpm_process.py` inputs docstring). Once the H5 keyed-map affordance (§4) is depended
  upon, the `fates` port is declared with it and the pre-init is deleted.

The fates-store + per-cell wiring is built via a **shared helper** (§4, H4), not hand-rolled. Conditions
are params, matching the existing slice's idiom: `baseline` (cue on, receptor competent),
`receptor_blocked` (cue on, `activated_type` λ forced to 0 — the biological intervention: receptor
antagonist / knockdown), `no_cue` (secretion rate 0 — shared negative control).

## 3. The three translator instances (in existing viva schema)

We realize translator *concepts* as concrete artifacts using schema fields that already exist and
render in the dashboard today but are currently unused workspace-wide. Each is landed as real framework
(§4), not one-off study data:

1. **Evidence translator (`calibrates`)** — add a chemokine–receptor entry to
   `workspace/references/papers.bib`, and attach `calibration_anchor{literature_target: <Kd>, cites:
   [<bibkey>], resolution: …}` to the receptor `kd` parameter in the receptor studies (and record it in
   the claim bundle's resolved hole). This fills `hole:lambda-calibration` with a cited value.
2. **Refinement translator + satisfaction condition** — the "refinement contract" the investigation
   requests, authored as a checked claim: the receptor realization **must preserve** {baseline recruits,
   response-block abolishes recruitment, no-cue abolishes recruitment} and reproduce the recruitment
   envelope of the phenomenological realization within tolerance. Encoded as (a) a `refinement:` block in
   `claim-bundle/realization.yaml` naming what must be preserved, and (b) a reusable check comparing the
   two realizations' recruitment observables (§4, H3).
3. **Scale translator (`aggregate`, cell→tissue)** — lift the recruitment-index reduction
   (`pbg_cpm_studies/chemotaxis/metrics.py: recruitment_index_from_coms`) into a **declared readout /
   typed emitted observable** with `index_by`/`units`, so the tissue-scale quantity is first-class rather
   than a post-hoc script (§4, H1).

Fields used already exist: `study.schema.json:790` (readout `cites`), `:831` (`calibration_anchor`),
`:866` (per-test `cites`), `:1324` (`literature_anchors`); dashboard renders threshold provenance
(`workspace/reports/assets/walkthrough.js`).

## 4. Infrastructure hardening & streamlining (co-equal workstream)

Each investigation need is landed as real, tested framework at the right layer, and **retires the glue
it replaces**. Net line count trends down; the new approach must not sit beside the old. Layering:
in-workspace (viva-cpm) is the proving ground; a piece is promoted upstream (viva-superpowers,
process-bigraph) only once the investigation proves it — each upstream change in its own worktree/PR.
The investigation can proceed even if an upstream promotion lags (build in-workspace first, promote when
proven; keeps the effort incremental and non-breaking).

| # | What hardens (add) | What streamlines (prune) | Layer | Test |
|---|---|---|---|---|
| **H1** | Typed cell→tissue **aggregate readout** (the first real scale bridge): a small reusable reduction that emits a declared tissue-scale observable with `index_by`/`units`, replacing post-hoc scripts. | Consolidate the two bespoke reductions (`chemotaxis/metrics.py`, `gg1993/metrics.py`) toward **one shared typed reduction module**; recruitment-index migrates onto it. | viva-cpm pkg → (promote to viva) | reduction unit tests; recruitment-index parity vs current script |
| **H2** | Make `calibration_anchor`/`cites` **load-bearing**: evaluator/rigor computes divergence vs `literature_target` and surfaces it; `bands_missing_provenance` becomes an **enforced gap check** (uncited band = flagged gap, not silent pass). | Remove reliance on hand-set thresholds where a cited anchor exists; delete the parallel "aspirational" story for anchored bands. | viva-superpowers (evaluator/rigor) — upstream | anchor-divergence unit test; gap-check on a fixture study |
| **H3** | **Refinement/satisfaction check** as reusable infra: a general two-realization equivalence check on a shared claim observable (needed by the whole claim-layer, not just here). | Replaces would-be one-off comparison scripts; single check consumed by any multi-realization claim. | viva-superpowers (claim-layer) — upstream | preservation-conditions test on receptor vs phenomenological |
| **H4** | **Robust down-scale coupling helper**: a function that builds the `fates` store + per-cell subcell wiring from a cell list, encapsulating the map-key/pre-init contract (today fragile, documented as workarounds in `cpm_process.py`). | Delete hand-rolled per-composite fates wiring; one code path for every subcellular→CPM coupling. | viva-cpm (cpm pkg) | wiring helper unit test; both composites use it |
| **H5** | **Keyed-map write affordance** in bigraph-schema (the one foundational change): an *additive* Map wrapper/flag whose `apply` **creates** absent keys instead of silently dropping them. Today `apply(Map,…)` skips `key not in state` (`methods/apply.py:438-439`), so a producer owning the full keyset must pre-seed the store or wrap in `overwrite[map]`. Secondary: harden `append_link_path` (`core.py:104-114`) only if the keyed-map emit path needs the `"node"` leaf dodge. | Delete the `fates`-store pre-init dance; drop the `overwrite[map]`-just-to-dodge-key-drop idiom. | **bigraph-schema** (own worktree/PR); viva-cpm adopts on release | new affordance creates absent keys on apply; **existing `Map` apply unchanged** (regression test) |

**Foundational scope:** the *only* bigraph-schema change is the additive keyed-map affordance (H5) — **not**
the general translator registry / `resolve`-classification (Domain Bridges design P3), which stays
speculative-free until this and later use cases inform it. H5 is additive and non-breaking: existing
`Map` semantics (silent-skip of absent keys — load-bearing for sparse per-tick updates elsewhere) are
untouched; the new affordance is opt-in per port.

**Cross-repo sequencing:** H5 develops in a bigraph-schema worktree; viva-cpm's coupling keeps the H4
pre-init as a **fallback** until H5 is merged and depended upon, so the investigation is never blocked on
the upstream release. Develop/test the coupling against the H5 worktree via editable install / PYTHONPATH.

**Not in scope now (design-doc ledger items this slice does *not* force):** dual-v4 study-vocabulary
split, `_REGISTRY` public front door, `runs_meta` DDL consolidation, `agents/0/` emit-path heuristic —
none are exercised by this CPM investigation; they stay in the general Domain Bridges rollout. This keeps
hardening honest: we harden what we touch.

### 4.1 Existing translators to build on and consolidate (reconnaissance)

The ecosystem already has translator-role constructs; the hardening **generalizes and consolidates**
these rather than adding a parallel notion:

- **Typed, already declare source→target (the parents to generalize, not sibling):**
  `bigraph_schema/methods/transform.py:17` `transform(core, src_schema, tgt_schema, state)` (schema→schema
  migration); `process_bigraph/composite_spec.py:339,395` `to_document`/`to_composite` (spec↔runtime);
  `process_bigraph/types/process.py:50` `Bridge`/`Interface` (internal↔I/O vocabulary — note this
  `Bridge(Node)` *is* the Composite bridge, typed).
- **The incumbent "adapter" your vocabulary names:** `viva-superpowers/templates/model/viva_<model>/adapters.py`
  ("Step adapters bridging mismatched ports between wrapped processes") + the pass-through-vs-adapter rule
  in `viva-superpowers/skills/viva-expert/reference.md:1005` — effectively the informal spec of the
  translator concept.
- **Strongest existing instance (typed-by-convention, validates, registry-dispatched):**
  `pbg_cpm_studies/evaluators.py:45,85` `_recruitment_index` + `register_evaluators` — a backend-neutral
  `measure` → backend observable, keyed by kind, with validation. **Increment 1's studies reuse this exact
  extractor.** H2/H3 extend it; the measure-crossing formalism should *become* this, generalized.
- **Prose-only bindings to make executable:** `claim-bundle/realization.yaml` `bindings:` /
  `intervention-binding:` declare source→target + fidelity but as human strings — give them a typed,
  resolvable form so `pbg_cpm_studies/claim_bundle_closure.py` can check a binding names a real
  extractor/param, not just id cross-refs.
- **Untyped reductions to type (H1):** `cpm/metrics.py`, `pbg_cpm_studies/chemotaxis/metrics.py`,
  `pbg_cpm_studies/gg1993/metrics.py` (state→observable). H1 also **collapses the two parallel registries**
  — `evaluators.py` (by measure-kind) and `gg1993/validate.py:35 MEASURE` (by study-slug) — into one.

**Consolidation principle:** the typed translator generalizes `transform` / `Bridge` / the evaluators
registry; it does not sit beside them. Every H-item that lands must retire its ad-hoc predecessor above.

## 5. Studies, rigor, and the claim bundle

- New studies under `workspace/investigations/chemotactic-recruitment/studies/`:
  `recruitment-receptor-baseline` (positive), `recruitment-receptor-blocked` (intervention: receptor
  blocked). The existing `recruitment-adversarial` (no cue) is realization-agnostic and is cited as the
  shared negative control.
- **Rigor (fixed in this slice):** ≥5 seeds with reported CIs on the recruitment index; the receptor
  studies carry `calibration_anchor`/`cites` (now enforced by H2); any verdict that overclaims beyond
  the tested domain is corrected. Receptor `activate_occupancy`, `kd`, `hill`, and `conc_scale` are
  recorded as `model_settings[]` with `cites`/provenance for the data-constrained ones.
- **Claim bundle:** `realization.yaml` gains a second realization (`receptor-level`) with fidelity labels
  that mark the semantic gap as **resolved for the response term** (occupancy-driven, not
  phenomenological), plus the `refinement:` preservation block. The existing phenomenological realization
  is left intact as the coarse realization.

## 6. Testing

- **Unit — `ReceptorSubcell`:** occupancy monotonic in ligand; `θ(Kd) = 0.5`; fate = activated iff
  `θ ≥ activate_occupancy`; deterministic.
- **Composite invariant:** in `recruitment_receptor`, type-3 contact-J rows equal type-2 rows (adhesion
  unconfounded); every responder cell has a pre-initialized `fates` key (via the H4 helper).
- **Behavior (multi-seed):** baseline recruits (index rises, CI excludes 0); receptor_blocked and no_cue
  do not (CI includes 0). Reuse the existing `study_evaluator` measure path.
- **Refinement/satisfaction test (H3):** receptor-baseline recruitment envelope matches phenomenological
  baseline within tolerance; the three preservation conditions hold.
- **Dose-response:** recruitment rises with cue rate and its half-max coincides with `c ≈ Kd`.
- **Framework tests:** H1 reduction parity + units; H2 anchor-divergence + gap-check; H4 wiring helper;
  H5 keyed-map affordance creates absent keys on apply **and** existing `Map` apply is unchanged
  (regression), plus a coupling test showing `fates` writes land without pre-init under H5.

## 7. Scope and increments (hardening interleaved with biology)

**Increment 1 (biology MVP + coupling hardening):** `ReceptorSubcell`; H4 coupling helper;
`recruitment_receptor` composite; `recruitment-receptor-baseline` + `recruitment-receptor-blocked`
(multi-seed, CIs); the cited Kd evidence bridge (papers.bib + calibration_anchor); Figures 1–2 (§8). **H5
develops in parallel** in a bigraph-schema worktree; viva-cpm adopts it when merged (H4 pre-init is the
fallback until then), at which point the pre-init is deleted.

**Increment 2 (evidence infra + dose-response):** H2 (calibration_anchor load-bearing + gap check);
dose-response sweep (Figure 3-data); half-max-at-Kd check.

**Increment 3 (scale bridge + refinement):** H1 (typed aggregate readout, consolidate metrics); H3
(refinement/satisfaction check) + Figure 3; claim-bundle second realization.

**Deferred (YAGNI):** receptor adaptation/desensitization; graded multi-class responsiveness; continuous
per-cell λ (would require a modest Rust addition — `delta_chemotaxis` already threads `new_owner:
CellId`); importing a full published SBML chemotaxis model via `SBMLSubcell`; promoting H1/H2/H3 from
viva-cpm into upstream packages (do once proven here).

## 8. Figures — pulled out by the need for citable evidence

1. **Receptor dose-response** — recruitment index (and mean occupancy) vs chemokine concentration, with
   the **cited Kd/EC50 marked** and the literature dose-response band shaded. Exists only because we
   needed to cite a Kd.
2. **Spatial activation map** — lattice colored by per-cell occupancy / activation class, showing the
   activation front. A readout the scalar λ cannot produce.
3. **Refinement preservation** — receptor-realization vs phenomenological-λ recruitment, overlaid with CI
   ribbons: the satisfaction condition (and the translator "round-trip") made visible.

New viz functions in `pbg_cpm_studies/visualizations/` (`@as_visualization`), following the existing
`chemotaxis_studies.py` pattern.

## 9. How this feeds the general formalism

Each instance is the concrete seed of a general translator class: the cited Kd → **evidence bridge
(`calibrates`)**; the refinement contract → **refinement bridge + satisfaction condition**; the
recruitment-index readout → **scale bridge (`aggregate`)**. Building them against a real investigation
tests the design's claims (Does `calibration_anchor` carry what an evidence bridge needs? Is the
satisfaction condition expressible as a runnable check? Is the reduction a natural scale bridge?) before
any framework is committed. H1–H4 and their prunes are the first entries of the general
delete/consolidate ledger, executed for real. Findings feed back into
`project_domain_bridges_formalization`.

## 10. Open questions

1. **Field units → concentration.** `conc_scale` maps CPM field units to the Kd's concentration units;
   its value is a modeling choice that must be recorded (provenance `theory`) rather than tuned to fit.
2. **Chemokine/receptor choice.** Neutrophil/chemokine selected; the specific ligand–receptor pair
   (e.g. fMLP–FPR1 vs CXCL8–CXCR1/2) determines the exact Kd and citation — pick during implementation
   from a well-characterized affinity.
3. **Two classes vs graded.** Start naive/activated (threshold); revisit graded if the dose-response
   needs more than a step.
4. **Promotion boundary.** H2/H3 touch viva-superpowers (upstream). Confirm the ownership split from the
   Domain Bridges design (§7.2 there): does evidence/refinement infra live in viva-superpowers or lower?
   Build in-workspace first regardless.
5. **H5 realization.** Settle in the implementation plan whether the keyed-map affordance is a new wrapper
   type (e.g. `assign[map[…]]` / an `Overwrite`-adjacent per-key variant) or a `Map` flag
   (`_open`/`_autocreate`). Prefer whichever composes with the existing grammar and needs one `apply`
   overload. Must leave bare `Map` apply byte-for-byte unchanged.
6. **Cross-repo dependency.** viva-cpm must run against the H5-fixed bigraph-schema (editable/PYTHONPATH
   during dev; a version bump once released). Track so the pre-init deletion and the dependency bump land
   together, not before.
