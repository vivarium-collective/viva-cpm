# Influenza (Sego 2022) reproduction — Increment 5 Implementation Plan (RE-SCOPED: macrophages)

> **For agentic workers:** REQUIRED SUB-SKILL: superpowers:subagent-driven-development. Steps use `- [ ]`.

## RE-SCOPE NOTE (read first)
The design spec's §5 ladder listed Increment 5 as "chemokine + IL-10 fields" and Increment 6 as
"macrophages". Reading the actual CC3D source (`ChemokineSecretionSteppable` / `IL10SecretionSteppable`),
**both chemokine and IL-10 are secreted by MACROPHAGES** (`for cell in cell_list_by_type(MACROPHAGE)`), not
by dead cells — so those fields have no source until macrophages exist. We therefore INVERT the order:
**Increment 5 = macrophages** (the field sources, and they chemotax up the already-existing virus field and
phagocytose virus/dead cells), and Increment 6 = the macrophage-released chemokine + IL-10 fields. This is a
faithful-to-source re-scope, recorded here and in the investigation.

**Goal:** Introduce the first local immune cell — the macrophage — recruited to the epithelial sheet, chemotaxing up the extracellular-virus gradient toward the infection and phagocytosing virus, so macrophages localize to the lesion (paper Fig 2B / Fig 3A immune layer).

**Architecture:** Reuse existing engine primitives — macrophage is cell type M=4 (already in `types.py`), added to the world with adhesion (Table 3 immune J) and volume; chemotaxis up the virus field via the existing `set_chemotaxis(virus_fi, M, lambda)`; virus phagocytosis via a virus-field UPTAKE (negative per-type secretion, or a per-cell uptake using the Increment-3 secretion scale). Recruitment is STUBBED (a simple fixed inflow) — real recruitment is Hill-driven by chemokines/APCs (Increment 6 / the global ODE Increment 8). NO Rust change expected (chemotaxis + negative-secretion uptake already exist; confirm in Task 5.1).

**Tech Stack:** Python 3.12, `cpm.cpm_core`, `process_bigraph`, pytest, numpy, PyYAML. Builds on merged Increments 0–4.

**Spec:** `docs/superpowers/specs/2026-09-18-influenza-sego2022-reproduction-design.md` (Increment 6 row of §5 — "Macrophage response"; note the re-scope above swaps the 5/6 order).

## Global Constraints
- CC3D source authoritative. Re-fetch to scratch per `docs/cc3d-reference/sego2022-source-notes.md`; don't vendor.
- Macrophage: type M=4 (`types.py`). Chemotaxis up virus `chemotaxis_v_macro = 5000` (Table 3; source `ImmuneModelInputs.py:chemotaxis_v_macro = 5E3`). Adhesion (Table 3): uninfected–immune 20, infected–immune 10, dead–immune 20, homotypic-immune 25 — use the SOURCE's full matrix (Task 5.0). Volume constraint as for epithelial (25 sites) unless the source differs.
- Chemotaxis functional-form gap: the paper uses the saturating/COM form `λc·c/(1+c_CM)`; the engine uses the linear CC3D form `ΔH=-λ(c_dest-c_source)`. Use the existing LINEAR form for macrophage localization (behavior: macrophages climb the virus gradient) and NOTE the gap (Increment-9 calibration). Do NOT add the saturating form unless trivial.
- Recruitment STUB: seed macrophages at a simple fixed/first-order rate at the immune layer (interspersed / on medium). Real recruitment (Hill on chemokines C / APCs P) is Increment 6/8 — clearly flagged, not faked as calibrated.
- Phagocytosis: macrophages take up extracellular virus (reduce the local virus field). Source `MacrophageInternalizationSteppable`/equivalent — extract the uptake form in Task 5.0.
- params.yaml single authority. No AI attribution. Death/type changes via set_cell_type (no removal). Work only in worktree `~/code/viva-cpm--influenza-incr5` (branch `investigation/influenza-incr5`); never touch `~/code/viva-cpm`, `viva_cpm_studies/composites/__init__.py`, or `viva_cpm_studies/visualizations/` (peer-owned). Tests: `~/code/viva-cpm--influenza-incr5/.venv/bin/python -m pytest <path> -v`.
- Rust rebuild (only if a Rust change proves necessary): `VIRTUAL_ENV="$PWD/.venv" uvx maturin develop --release`.

---

### Task 5.0: Extract & cite macrophage constants + document the re-scope
**Files:** `params.yaml` (+`macrophage:` section), `docs/cc3d-reference/sego2022-parameters.md` (macrophage rows), `tests/test_influenza_params.py` (+test).
- [ ] Re-fetch source. From `ImmuneModelInputs.py`/`ImmuneModelLib.py` + the macrophage steppables, read: `chemotaxis_v_macro` (5000), the macrophage adhesion J row (uninfected/infected/dead/homotypic-immune), macrophage volume, the virus-uptake/phagocytosis rate, the recruitment/seeding parameters (seeding fraction, local fraction macrophage 100%, and the Hill recruitment form — record it even though wired later), and how immune cells enter (the immune-layer/seeding geometry). Cite file+symbol+units + conversions.
- [ ] Add `macrophage:` section (chemotaxis_v_macro, adhesion_*, volume, virus_uptake, seeding params, source). Test asserts keys + numeric + source. Also add a short note in the transcription doc + the section documenting the 5/6 re-scope (fields are macrophage-sourced).
- [ ] Commit `feat(influenza): cite macrophage constants + record 5/6 re-scope (Incr 5)`.

### Task 5.1: Macrophage in the world + chemotaxis up virus
**Files:** `viva_cpm_studies/influenza/immune.py` (new; macrophage seeding + world wiring), `viva_cpm_studies/influenza/run.py` (extend the driver to add macrophages + chemotaxis), `tests/test_influenza_macrophage.py`.
**Interfaces:** `immune.seed_macrophages(world, spec, n_or_frac, rng)` — add macrophage (type M) cells to the world (interspersed / on medium), with the Table-3 adhesion + volume; `immune.set_macrophage_chemotaxis(world, virus_fi)` — `set_chemotaxis(virus_fi, M, chemotaxis_v_macro)`. Extend `run.run_epithelial_fate` (or a new `run_with_macrophages`) to seed macrophages and enable virus-chemotaxis each update.
- [ ] Confirm the engine's chemotaxis works for a NEW type M on the virus field (set_contact for M with all types, set_chemotaxis). If seeding macrophages into a confluent sheet needs medium space, decide the geometry (e.g. a sparser sheet, or macrophages replacing some medium at the boundary) — document it.
- [ ] Failing INTEGRATION test: seed a virus blob / infected patch, add macrophages away from it, run; assert macrophages' mean distance to the virus centroid DECREASES over time (they climb the gradient / localize to the infection). Deterministic/seeded, small/fast. (Mirror the chemotactic-recruitment study's localization readout if useful.)
- [ ] Implement, PASS. Commit `feat(influenza): macrophage cell type + chemotaxis up the virus field`.

### Task 5.2: Macrophage phagocytosis of virus
**Files:** `viva_cpm_studies/influenza/immune.py` (+uptake), `viva_cpm_studies/influenza/run.py`, `tests/test_influenza_macrophage.py` (+test).
**Interfaces:** wire macrophage virus UPTAKE — reduce the virus field where macrophages are (source-cited rate). Prefer a per-type negative secretion `set_secretion(virus_fi, M, -uptake)` if the engine supports it, else per-cell uptake via `set_cell_secretion_scale` on a negative base; confirm the mechanism in code and document.
- [ ] Failing test: with macrophages present + chemotaxing to the infection, total_virus at a comparable step is LOWER than a no-macrophage run from the same seed (macrophages clear virus). Deterministic.
- [ ] Implement, PASS. Commit `feat(influenza): macrophage phagocytosis (virus uptake)`.

### Task 5.3: macrophage-response study + viz + membership
**Files:** `workspace/studies/macrophage-response/study.yaml` (hand-author), `workspace/investigations/influenza-sego2022/investigation.yaml` (+member + re-scope note), `viva_cpm_studies/influenza/viz.py` (+figure), `tests/test_influenza_viz.py` (+test). Do NOT touch visualizations/ or composites/__init__.py.
- [ ] Add `viz.macrophage_response_figure(run_result)` (macrophage distance-to-infection over time + virus with/without macrophages, and/or a snapshot showing macrophages clustered at the lesion). Test returns a Figure (Agg).
- [ ] Hand-author `macrophage-response/study.yaml`: report HONESTLY the macrophage localization (mean distance to infection decreases) + virus reduction (cite the actual Task 5.1/5.2 numbers). Caveats: recruitment is STUBBED (real Hill recruitment on chemokines/APCs is Increment 6/8); chemotaxis uses the engine's LINEAR form not the paper's saturating/COM form (Increment-9 calibration); chemokine/IL-10 fields + NK/CD8 + global ODE are later increments; reproduction PENDING (Increment 9). verdict = documented. Record the 5/6 re-scope in the investigation overview. Wire a baseline composite (reuse or add a macrophage composite in composites/influenza.py — auto-registers via the peer __init__ import; do NOT edit __init__.py). Validate with lint-workspace.py; add investigation member.
- [ ] Run `-k influenza` (all pass). Commit `feat(influenza): macrophage-response study + viz`.

## Self-Review notes
- Spec coverage: the "Macrophage response" row (§5 Incr 6) → Tasks 5.0–5.3, brought forward per the re-scope; the chemokine+IL-10 field row (§5 Incr 5) moves to the next increment (macrophage-released).
- Stubs clearly flagged: recruitment (Hill on chemokines/APCs → Incr 6/8), chemotaxis functional form (linear vs saturating → Incr 9), no chemokine/IL-10 yet, no NK/CD8 (Incr 7), no global ODE (Incr 8).
- Reuses existing primitives (chemotaxis, secretion, contact) — confirm no Rust change needed in Task 5.1; if virus uptake needs a Rust tweak, keep it minimal + backward-compatible + Rust-tested.
