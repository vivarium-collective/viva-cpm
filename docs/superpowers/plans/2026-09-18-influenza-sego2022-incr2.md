# Influenza (Sego 2022) reproduction — Increment 2 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax.

**Goal:** Add the extracellular-virus diffusive field and the stochastic epithelial infection transition (H→I) driven by local virus, so a seeded lesion spreads across the sheet — the first field-coupled mechanism of the reproduction.

**Architecture:** Reuse the existing Rust engine primitives (`add_field`, `set_secretion`, `advance_fields`/field stepping, `field_mean_at_cell`, `set_cell_type`) — NO Rust changes. Infected epithelial cells secrete virus at a flat rate (resistance = 0 until the IFN field arrives in Increment 3); a Python transition Step reads each uninfected cell's local mean virus and stochastically flips it to infected. Wire both into a `virus-infection` composite over the Increment-1 sheet.

**Tech Stack:** Python 3.12, `cpm.cpm_core` (existing CPM+field engine), `process_bigraph`, numpy, pytest, PyYAML.

**Spec:** `docs/superpowers/specs/2026-09-18-influenza-sego2022-reproduction-design.md` (Increment 2 row of §5). Builds directly on Increment 1 (merged): `pbg_cpm_studies/influenza/{params.yaml,sheet.py,build.py,types.py}` and `pbg_cpm_studies/composites/influenza.py`.

## Global Constraints

- CC3D source (ViralInfectionVTM) is the authoritative parameter/spec; source wins over the paper. Re-fetch it to a scratch dir per `docs/cc3d-reference/sego2022-source-notes.md` when exact constants are needed; do NOT vendor it into the repo.
- Virus field units (engine field units, from `docs/cc3d-reference/sego2022-parameters.md` §3): unitless diffusion **D ≈ 0.1788 lattice²/MCS**, decay **≈ 2.861e-4 /MCS**; target diffusion length **5 cell diameters** (= 25 lattice sites, since 1 cell = 5 sites).
- Virus secretion: infected cells release at flat rate `g_vi` scaled by `(1 − resist)`; **resist = 0 in Increment 2** (no IFN field yet) → full release. The `(1−resist)` gating is deferred to Increment 3.
- Infection transition (source `ViralInternalizationSteppable`, transcription §6): per uninfected cell per MCS-update, `rate = g_hv · v̄(s)` (v̄ = local mean virus over the cell's pixels, no saturation), `Pr(infect) = 1 − exp(−rate)`; on success type H→I. Source flips straight to the releasing/infected type — we use `types.I`.
- Cell-type codes (existing `pbg_cpm_studies/influenza/types.py`): MEDIUM=0, H=1, I=2, D=3, M=4, K=5, E=6. (D and immune types are not yet populated in Increment 2.)
- Honesty: the `virus-field-infection` study validates the MECHANISM (field diffusion length, infection spread), not the full Fig-3B/5/7 reproduction — that stays PENDING (Increment 9). No fabricated reproduction verdict. No AI attribution in commits.
- Work only in the worktree `~/code/viva-cpm--influenza-incr2` (branch `investigation/influenza-incr2`); never touch `~/code/viva-cpm`. Tests: `/Users/eranagmon/code/viva-cpm--influenza-incr2/.venv/bin/python -m pytest <path> -v`.

---

### Task 2.0: Extract & cite virus/infection constants into params.yaml

**Files:**
- Modify: `pbg_cpm_studies/influenza/params.yaml` (add a `virus` section)
- Modify: `docs/cc3d-reference/sego2022-parameters.md` (fill exact `g_vi`, `g_hv`, initial-condition constants if not already present)
- Test: `tests/test_influenza_params.py` (extend)

**Interfaces:**
- Produces: `load_params()["virus"]` with keys `diffusion_lat2_per_mcs` (≈0.1788), `decay_per_mcs` (≈2.861e-4), `diffusion_length_cell_diam` (5), `secretion_g_vi` (from source), `infection_g_hv` (from source), and the initial-condition constants (initial viral load / initial infection fraction) needed to seed a run.

- [ ] **Step 1:** Re-fetch the CC3D source to a scratch dir (per source-notes). From `Simulation/ViralInfectionVTMModelInputs.py` + `ImmuneModel/ImmuneModelInputs.py`, read the literal virus **secretion rate** (`g_vi` / `secr_amount` for infected cells; source symbol e.g. `virus_secr_...`/`p_vi`), the **infection/internalization rate** (`g_hv`; source symbol e.g. `infection_...`/`p_hv` used in `ViralInternalizationSteppable`), and the unitless field D/decay (`virus_dc_im`/`virus_decay_im`) to confirm the ≈0.1788/≈2.861e-4 values in the transcription doc. Cite each with its source file+symbol. If a rate is expressed per-second/per-day, convert to per-MCS-update using `s_to_mcs=60` and the composite's `mcs_per_update` (document the conversion).
- [ ] **Step 2:** Write the failing test: `test_params_carry_virus_section` asserting `load_params()["virus"]` has all keys above with numeric values and a `source` string, and that `diffusion_length_cell_diam == 5`.
- [ ] **Step 3:** Add the `virus:` section to params.yaml (values from Step 1, each with provenance) and fill any missing exact constants into the transcription doc's virus rows.
- [ ] **Step 4:** Run the test → PASS. Commit: `feat(influenza): cite virus field + infection constants (Incr 2)`.

### Task 2.1: Virus field on the sheet (diffusion + infected-cell secretion)

**Files:**
- Create: `pbg_cpm_studies/influenza/fields.py`
- Test: `tests/test_influenza_fields.py`

**Interfaces:**
- Consumes: `build.world_from_spec` (Increment 1), `cpm_core.World.add_field/set_secretion/set_field_dynamics/advance_fields/field_conc/field_mean_at_cell`, `load_params()["virus"]`, `types`.
- Produces: `fields.add_virus_field(world) -> int` (returns the virus field index) — adds the virus field with the params' D/decay, sets `types.I` secretion to `g_vi` (resist=0 stub), and configures field dynamics for stability. A helper `fields.virus_diffusion_length_sites() -> float` returning the expected 25 (5 cell diameters × 5 sites).

- [ ] **Step 1:** Write failing tests: (a) after adding the field, seeding one infected cell, and advancing fields + stepping, the virus concentration is positive near the infected cell and decays with distance (measure the ~1/e radius is within a factor of ~2 of 25 sites — behavior, not exact); (b) an all-uninfected sheet keeps virus at 0. Use `build.world_from_spec(sheet.build_sheet_spec(0.3))`, set one cell's type to `types.I` via `set_cell_type`, `add_virus_field`, advance.
- [ ] **Step 2:** Run → FAIL (module missing).
- [ ] **Step 3:** Implement `fields.py`. Use `add_field("virus", D, decay)` with the params values, `set_secretion(fi, types.I, g_vi)`, `set_field_dynamics(fi, dt, substeps)` chosen for explicit-scheme stability (`dt*D*2*ndim < 1`; sub-step as needed — see `cpm/subcellular`/`test_fields.py` for the pattern). Read D/decay/g_vi from params.
- [ ] **Step 4:** Run → PASS. Commit: `feat(influenza): extracellular virus field + infected-cell secretion`.

### Task 2.2: Stochastic infection transition (H→I) Step

**Files:**
- Create: `pbg_cpm_studies/influenza/transitions.py`
- Test: `tests/test_influenza_transitions.py`

**Interfaces:**
- Consumes: `load_params()["virus"]["infection_g_hv"]`, `types`, numpy RNG.
- Produces: `transitions.infection_step(types_list, virus_at_cell, g_hv, rng) -> list` — a PURE function returning the new types list: for each cell of type `H`, draw `Pr = 1 − exp(−g_hv · v̄)` and flip to `I` on success; all other cells unchanged. `virus_at_cell` is a dict/list of per-cell mean virus (from `field_mean_at_cell`). Keeping it pure makes it deterministic-per-seed and unit-testable without the engine.

- [ ] **Step 1:** Write failing tests (pure function, seeded numpy RNG): (a) an H cell with high v̄ (rate ≫ 1) infects with Pr≈1 → becomes I; (b) v̄=0 → never infects; (c) already-I/other-type cells are never changed; (d) same seed → same result (determinism); (e) with a moderate rate over many cells, the empirical infection fraction ≈ `1−exp(−g_hv·v̄)` within tolerance.
- [ ] **Step 2:** Run → FAIL.
- [ ] **Step 3:** Implement `transitions.infection_step`.
- [ ] **Step 4:** Run → PASS. Commit: `feat(influenza): stochastic virus-driven infection transition H->I`.

### Task 2.3: Virus-infection composite + spread integration test

**Files:**
- Modify: `pbg_cpm_studies/composites/influenza.py` (add `virus_infection` factory)
- Test: `tests/test_influenza_virus_infection.py`

**Interfaces:**
- Consumes: `sheet.build_sheet_spec`, `build.world_from_spec`, `fields.add_virus_field`, `transitions.infection_step`, the CPMProcess.
- Produces: `composites.influenza.virus_infection(patch_mm, seed, initial_infected_frac_or_load)` returning a composite document, and a driver-style helper `run_virus_infection(patch_mm, steps, seed) -> dict` (returns per-step counts of H/I and total virus) usable by the study + integration test. Wiring: each update = CPM step (motility) → advance virus field (secretion+diffusion+decay) → read `field_mean_at_cell` per cell → `infection_step` → write back types via `set_cell_type`. Reuse the CPMProcess where it fits; if the transition/field loop is simpler as a small explicit driver over the Rust world (like `pbg_cpm_studies/gg1993` driver), do that and keep the composite as the dashboard-runnable wrapper.

- [ ] **Step 1:** Write failing integration test: seed a small sheet (0.3mm) with a few infected cells (or an initial virus blob), run `run_virus_infection` for N updates, assert (a) infected-cell count strictly increases over time (lesion spreads), (b) spread is LOCAL early (infected cells are spatially contiguous with the seed, not random across the sheet), (c) with zero initial virus AND zero initial infected, nothing infects.
- [ ] **Step 2:** Run → FAIL.
- [ ] **Step 3:** Implement the factory + `run_virus_infection` driver.
- [ ] **Step 4:** Run → PASS. Commit: `feat(influenza): virus-infection composite + lesion-spread test`.

### Task 2.4: Study + viz + investigation membership

**Files:**
- Create: `workspace/studies/virus-field-infection/study.yaml`
- Modify: `workspace/investigations/influenza-sego2022/investigation.yaml` (add member)
- Create/Modify: `pbg_cpm_studies/influenza/viz.py` (add a virus+infection figure)
- Test: `tests/test_influenza_viz.py` (extend)

**Interfaces:** hand-author the study (no dashboard server); model on `workspace/studies/epithelial-sheet-baseline/study.yaml`.

- [ ] **Step 1:** Add `viz.virus_infection_figure(run_result)` → a matplotlib Figure: virus-field heatmap at a late step + infected-cell overlay, and an infected-count-vs-time curve. Test asserts it returns a Figure with axes (Agg backend).
- [ ] **Step 2:** Hand-author `virus-field-infection/study.yaml` (schema v3, investigation influenza-sego2022, phase Simulate, status in-progress): report the measured virus diffusion length (vs the 5-cell-diameter target) and the lesion-spread behavior HONESTLY; reproduction verdict PENDING (Increment 9). Add it as an investigation member. Validate with `scripts/lint-workspace.py`.
- [ ] **Step 3:** Run the full influenza suite (`pytest tests/ -k influenza`) → all PASS. Commit: `feat(influenza): virus-field-infection study + viz`.

## Self-Review notes
- Spec coverage: Increment 2 row (virus field + infection H→I; front speed, 5-cell-diam length; Fig 6 virus target) → Tasks 2.0–2.4.
- No Rust changes: uses existing field + set_cell_type primitives (confirmed present in `crates/cpm-py/src/lib.rs`: add_field/set_secretion/set_field_dynamics/advance_fields/field_conc/field_mean_at_cell/set_cell_type).
- Resistance is 0 (full virus release) this increment; `(1−resist)` gating + the IFN field are Increment 3. The `virus` params section and secretion call are structured so Increment 3 only needs to multiply by `(1−resist)`.
- Type consistency: `add_virus_field(world)->int`, `infection_step(types_list, virus_at_cell, g_hv, rng)->list`, `run_virus_infection(...)->dict`, `virus_infection(...)` composite — used identically across 2.1–2.4.
