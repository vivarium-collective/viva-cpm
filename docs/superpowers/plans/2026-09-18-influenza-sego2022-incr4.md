# Influenza (Sego 2022) reproduction — Increment 4 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: superpowers:subagent-driven-development. Steps use `- [ ]`.

**Goal:** Complete the epithelial cell-fate lifecycle: infected cells die (`μ_i·(1−resist)`), uninfected cells die by the cellularized Allee effect when insufficiently surrounded by healthy tissue, and dead cells recover (D→H) by the symmetric Allee effect — so a lesion of dead cells forms and can recover (paper Fig 4A).

**Architecture:** One new Rust query — **contact area by cell type** (`cell_contact_area_by_type(cell_id) -> {type: area}`), the surface-fraction substrate the Allee formulas need — plus Python transition Steps for the three fates. Infected death reuses Increment 3's per-cell resist; Allee death/recovery use the new surface query. Death is `set_cell_type(→D)` (dead cells persist and can recover), NOT removal — so the Increment-3 `cell_scale`/`remove_cells` caveat stays inert.

**Tech Stack:** Rust (`crates/cpm-core`, `crates/cpm-py`), Python 3.12, `cpm.cpm_core`, pytest, numpy, PyYAML. Builds on merged Increments 0–3.

**Spec:** `docs/superpowers/specs/2026-09-18-influenza-sego2022-reproduction-design.md` (Increment 4 row of §5).

## Global Constraints
- CC3D source authoritative (source wins). Re-fetch to scratch per `docs/cc3d-reference/sego2022-source-notes.md`; don't vendor.
- **Transitions (transcription Table 2, §"death/Allee"):**
  - Infected death Î→D̂: `rate = μ_i·(1−resist)` (source `ViralCellDeathSteppable`, an explicitly SIMPLIFIED flat rate — omits the paper's printed ROS/Hill term; mirror the SIMPLIFIED source form, note the conceptual target). `Pr = 1−exp(−rate)`.
  - Uninfected Allee death Ĥ→D̂: `rate = b_h·(1−resist)·srf_area_oi·(srf_thresh − srf_uninfected)/srf_area_total²`, applied only when local uninfected surface fraction < threshold `θ_local` (source `RecoverySteppable`). srf_uninfected = area of the cell's surface contacting uninfected (H) cells; srf_area_total = total surface.
  - Recovery D̂→Ĥ (Allee a_H): structurally symmetric to the Allee death — a dead cell recovers when SUFFICIENTLY surrounded by healthy (H) tissue (swap the uninfected↔dying roles in the surface-fraction). Extract the exact source form in Task 4.0.
  - NK/CD8⁺ contact-killing γ terms are Increment 7 — DEFER.
- resist = f̄/(a_rf+f̄) from Increment 3 (`resistance.cell_resistance`); consumers use (1−resist). Cell types: MEDIUM=0,H=1,I=2,D=3 (immune 4-6 not populated yet).
- params.yaml single authority. No AI attribution. Death is set_cell_type (no removal).
- Work only in worktree `~/code/viva-cpm--influenza-incr4` (branch `investigation/influenza-incr4`); never touch `~/code/viva-cpm`, and **do NOT touch `pbg_cpm_studies/composites/__init__.py`** (a peer session owns a registration fix there). Tests: `~/code/viva-cpm--influenza-incr4/.venv/bin/python -m pytest <path> -v`.
- **Rust rebuild (VERIFIED):** after edits under `crates/`, from the worktree root run `VIRTUAL_ENV="$PWD/.venv" uvx maturin develop --release`. Rust unit tests via `cargo test -p cpm-core`.

---

### Task 4.0: Extract & cite death/Allee/recovery constants
**Files:** `params.yaml` (+`cell_death:` / `allee:` sections), `docs/cc3d-reference/sego2022-parameters.md` (fill exact `mu_i`, `b_h`, the Allee threshold `theta_local`/`srf_thresh`, and the recovery-rate form/constant), `tests/test_influenza_params.py` (+test).
- [ ] Re-fetch source. From `ImmuneModelInputs.py`/`ImmuneModelLib.py` + `ViralCellDeathSteppable`/`RecoverySteppable`, read: `mu_i` (infected death rate; per-day→per-MCS via s_t), `b_h` (Allee death coefficient), the surface threshold used in `RecoverySteppable` (`srf_thresh`/`theta_local` = θ_ODE/num_epithelial — record how it's computed), and the recovery-side rate/constant. Cite file+symbol+units + conversions.
- [ ] Add `cell_death:` (mu_i_per_mcs, source) and `allee:` (b_h, srf_threshold, recovery params, formula strings, source) sections. Test `test_params_carry_death_and_allee_sections` asserts keys + numeric + source.
- [ ] Commit `feat(influenza): cite epithelial death + Allee constants (Incr 4)`.

### Task 4.1 (RUST): contact-area-by-type query
**Files:** `crates/cpm-core/src/{field.rs or world.rs}`, `crates/cpm-py/src/lib.rs`, `crates/cpm-core/tests/contact_area.rs`, `tests/test_contact_area_by_type.py`.
**Interfaces produced:** `World.cell_contact_area_by_type(cell_id: int) -> dict[int, int]` — for the given cell, the count of face-neighbor pairs (boundary interfaces) between it and each neighbor cell type (including MEDIUM=0), i.e. per-type contact area. The total (sum over types) is the cell's contact surface.
- [ ] **Rust unit test first** (`contact_area.rs`): a small world with a known layout (e.g. one cell surrounded by cells of two types + medium), assert the returned per-type areas match the hand-counted face interfaces. Run `cargo test -p cpm-core` → FAIL.
- [ ] Implement: iterate the cell's pixels' face-neighbors (lattice.face_neighbors), and for each face crossing into a DIFFERENT owner, tally by that owner's cell_type. Return a HashMap<u16 type, i64 area> (or Vec). Reuse `lattice.face_neighbors` (crates/cpm-core/src/lattice.rs:131). Mirror `field_mean_at_cell`'s iteration style.
- [ ] `cargo test -p cpm-core` → PASS. Add the pyo3 binding. Rebuild (`VIRTUAL_ENV="$PWD/.venv" uvx maturin develop --release`). Write `tests/test_contact_area_by_type.py`: a seeded confluent sheet — a cell's per-type areas sum to its cell_surface; an interior H cell surrounded by H has ~all-H contact; injecting an I neighbor shows nonzero I-contact. Run → PASS.
- [ ] No-regression: `tests/test_influenza_*.py`, `tests/test_fields.py` → PASS. Commit `feat(engine): per-cell contact-area-by-type query (Allee surface substrate)`.

### Task 4.2: infected death transition (I→D)
**Files:** `pbg_cpm_studies/influenza/transitions.py` (+`infected_death_step`), `tests/test_influenza_transitions.py`.
**Interfaces:** `transitions.infected_death_step(types_list, resist_at_cell, mu_i, rng) -> list` — pure: for each `types.I` cell, `Pr = 1−exp(−mu_i·(1−resist))`; on success → `types.D`. Others unchanged.
- [ ] Failing tests: high mu_i·(1−resist) → I dies (→D); resist=1 (fully resistant) → never dies; non-I unchanged; determinism; empirical fraction ≈ 1−exp(−rate). Implement, PASS. Commit `feat(influenza): infected-cell death transition I->D (mu_i*(1-resist))`.

### Task 4.3: Allee death (H→D) + recovery (D→H)
**Files:** `pbg_cpm_studies/influenza/allee.py`, `pbg_cpm_studies/influenza/run.py` (extend the driver to apply all fate transitions per update), `tests/test_influenza_allee.py`.
**Interfaces:** pure functions `allee.uninfected_death_rate(srf_uninfected, srf_total, b_h, resist, srf_threshold)` and `allee.recovery_rate(...)` implementing the transcribed forms; a driver `run.run_epithelial_fate(...)` that each update reads per-cell resist (IFN) + `cell_contact_area_by_type` and applies infected death, Allee death, and recovery via set_cell_type.
- [ ] Failing tests: (a) the rate pure-functions match the transcribed formula at sample inputs (0 above threshold; positive below; monotone in the surface deficit). (b) INTEGRATION — the key claims: starting from a seeded infected patch, a contiguous DEAD region forms (H→I→D), AND a dead cell fully surrounded by H recovers to H over time (D→H), while a dead cell surrounded by dead/infected does NOT. Deterministic/seeded, small/fast.
- [ ] Implement, PASS. Commit `feat(influenza): cellularized Allee death (H->D) + recovery (D->H)`.

### Task 4.4: epithelial-fate study + viz + membership
**Files:** `workspace/studies/epithelial-fate/study.yaml` (hand-author), `workspace/investigations/influenza-sego2022/investigation.yaml` (+member), `pbg_cpm_studies/influenza/viz.py` (+figure), `tests/test_influenza_viz.py` (+test). Do NOT touch `pbg_cpm_studies/visualizations/` (peer owns).
- [ ] Add `viz.epithelial_fate_figure(run_result)` (cell-type composition H/I/D over time + a lesion snapshot showing dead region + recovered cells). Test returns a Figure (Agg).
- [ ] Hand-author `epithelial-fate/study.yaml` (schema v3, Simulate, in-progress): report HONESTLY the measured fate dynamics (lesion forms; recovery of surrounded dead cells) with the actual numbers from Task 4.3. Caveats: infected death is the SOURCE's SIMPLIFIED flat form (ROS/Hill term omitted — conceptual target noted); NK/CD8 contact-killing deferred (Increment 7); reproduction PENDING (Increment 9). verdict = documented. Wire a baseline composite (reuse virus_infection or add a fate variant in composites/influenza.py — it auto-registers via the peer's __init__ import; do NOT edit __init__.py). Validate with lint-workspace.py; add investigation member; keep status honest.
- [ ] Run `-k influenza` (all pass). Commit `feat(influenza): epithelial-fate study + viz`.

## Self-Review notes
- Spec coverage: Increment 4 row (uninfected death a_D + infected death + Allee recovery D→H; Table 2; Fig 4A) → Tasks 4.0–4.4.
- Second Rust change (4.1) — a read-only query, backward-compatible by construction (adds a method, changes no existing path). Rust unit test + Python test + no-regression required.
- Infected death mirrors the source's SIMPLIFIED flat μ_i·(1−resist); the paper's ROS/Hill term is the conceptual target, reinstated only if a later increment intends it (transcription §7 #6).
- Death = set_cell_type (dead cells persist + can recover) → the Increment-3 cell_scale/remove_cells caveat stays inert (still no cell removal).
- Deferred: NK/CD8 contact-killing (Increment 7); ROS field / global ROS (Increment 8); exact Allee threshold calibration is an Increment-9 concern if the surface-fraction magnitude diverges.
