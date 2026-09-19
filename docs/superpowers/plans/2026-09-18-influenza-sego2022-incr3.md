# Influenza (Sego 2022) reproduction — Increment 3 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: superpowers:subagent-driven-development. Steps use `- [ ]`.

**Goal:** Add the type-I IFN diffusive field and per-cell viral resistance ρ, and gate virus release by `(1−resist)` so accumulating IFN slows the infection spread — the first cell-state-regulated field source.

**Architecture:** One small, general Rust primitive — a **per-cell secretion scale** (`set_cell_secretion_scale(field_idx, cell_id, scale)`, default 1.0) applied in `advance_fields` as `rate * scale[owner] * dt`. Everything else is Python: an IFN field (like the virus field), a per-cell resistance computed from local IFN (`resist = f̄/(a_rf + f̄)`), and a driver that each step recomputes resist and sets each infected cell's virus secretion scale to `(1 − resist)`. This one primitive serves every later cell-state-regulated source (IL-10, TNF-gated chemokine, Allee death).

**Tech Stack:** Rust (`crates/cpm-core`, `crates/cpm-py`, maturin/pyo3), Python 3.12, `cpm.cpm_core`, pytest, numpy, PyYAML. Builds on merged Increments 0–2.

**Spec:** `docs/superpowers/specs/2026-09-18-influenza-sego2022-reproduction-design.md` (Increment 3 row of §5).

## Global Constraints
- CC3D source authoritative (source wins). Re-fetch to scratch per `docs/cc3d-reference/sego2022-source-notes.md`; don't vendor.
- **Resistance (source-literal, transcription §"Cellular viral resistance ρ"):** `resist = f̄/(a_rf + f̄)` — a plain Hill on local mean IFN `f̄ = amountSeenByCell/cell.volume` (NO leading `1−`, NO θ factor). Consumers apply `(1 − resist)` in a protective role — here: virus release `g_vi·(1−resist)`. At f̄=0, resist=0 → full release (consistent with Increment 2). `a_rf` is the source ODE param — extract its literal value in Task 3.0.
- **IFN field (transcription §3/§5):** unitless D ≈ 7.792 lat²/MCS, decay ≈ 0.07792/MCS, diffusion length 2 cell-diameters; infected cells release IFN at basal rate `g_fp` (APC amplification is global → Increment 8, omit here); IFN uptake by infected cells `g_fi` — include if simple, else note deferred.
- params.yaml is the single parameter authority. Per-pixel secretion uses the same z=1/per-pixel/cell_sites convention Increment 2 established (see fields.py). No AI attribution.
- Work only in worktree `~/code/viva-cpm--influenza-incr3` (branch `investigation/influenza-incr3`); never touch `~/code/viva-cpm`. Tests: `/Users/eranagmon/code/viva-cpm--influenza-incr3/.venv/bin/python -m pytest <path> -v`.
- **Rust rebuild (VERIFIED command):** after ANY edit under `crates/`, rebuild the extension into the worktree venv before running Python tests, from the worktree root:
  `VIRTUAL_ENV="$PWD/.venv" uvx maturin develop --release`
  (confirmed working: ~4s, installs editable into `.venv`). Rust unit tests run via `cargo test -p cpm-core` (cargo 1.95 present, compiles cleanly).

---

### Task 3.0: Extract & cite IFN field + resistance constants
**Files:** `pbg_cpm_studies/influenza/params.yaml` (+`ifn:` and `resistance:` sections), `docs/cc3d-reference/sego2022-parameters.md` (fill exact `g_fp`, `g_fi`, `a_rf`), `tests/test_influenza_params.py` (+test).
- [ ] Re-fetch the CC3D source. From `ImmuneModelInputs.py`/`ImmuneModelLib.py`/`ViralInfectionVTMSteppables.py` (`Type1InterferonSecretionSteppable`, `Type1InterferonModelSteppable.update_resistance`), read: IFN unitless D (`t1ifn_dc`≈7.792) + decay (`t1ifn_decay`≈0.07792), diffusion length 2 (`exp_t1ifn_dl`); basal IFN secretion `g_fp`; IFN uptake `g_fi` (note if deferring); and `a_rf` (the resistance denominator constant). Cite file+symbol+units and any per-day→per-MCS conversion.
- [ ] Add `ifn:` (diffusion_lat2_per_mcs, decay_per_mcs, diffusion_length_cell_diam=2, secretion_g_fp, uptake_g_fi?, source) and `resistance:` (a_rf, formula string `resist = f_bar/(a_rf + f_bar)`, source) sections; test `test_params_carry_ifn_and_resistance_sections` asserts keys + numeric + source + `ifn.diffusion_length_cell_diam == 2`.
- [ ] Commit `feat(influenza): cite type-I IFN field + resistance constants (Incr 3)`.

### Task 3.1 (RUST): per-cell secretion scale primitive
**Files:** `crates/cpm-core/src/field.rs`, `crates/cpm-core/src/world.rs` (plumbing if needed), `crates/cpm-py/src/lib.rs` (pyo3 binding), `crates/cpm-core/tests/field_scale.rs` (Rust unit test), `tests/test_cell_secretion_scale.py` (Python binding test).
**Interfaces produced:** `World.set_cell_secretion_scale(field_idx: int, cell_id: int, scale: float)` (pyo3) — sets a per-(field,cell) multiplier (default 1.0) applied to that cell's pixels' secretion in `advance_fields`.
- [ ] **Rust unit test first** (`crates/cpm-core/tests/field_scale.rs`): build a small world with one secreting cell, add a field with a per-type rate, set that cell's scale to 0.5, advance, and assert the accumulated concentration is half the unscaled case (and scale=0 gives no secretion). Run `cargo test -p cpm-core` → FAIL (method missing).
- [ ] Implement in `field.rs`: add per-field per-cell scale storage (e.g. `cell_scale: HashMap<CellId,f64>` or `Vec<f64>` len n_cells default 1.0) on `Field` (or a parallel structure), a `set_cell_secretion_scale(field_idx, cell_id, scale)` method, and modify the secretion loop in `advance_fields` (currently `field.rs:92-99`):
  ```rust
  let owner = self.lattice.owner(idx);            // already available
  let t = self.cells[owner as usize].cell_type as usize;
  let rate = self.fields[fi].secretion.get(t).copied().unwrap_or(0.0);
  if rate != 0.0 {
      let scale = self.fields[fi].cell_scale_for(owner);   // default 1.0
      let dt = self.fields[fi].dt;
      self.fields[fi].conc[idx] += (rate * scale * dt) as f32;
  }
  ```
  Keep the default path (no scale set) numerically identical to today (scale 1.0) so Increment-2 field tests don't regress. Reset/resize cell_scale on add_cell/remove as needed so ids stay valid.
- [ ] `cargo test -p cpm-core` → PASS. Add the pyo3 binding `set_cell_secretion_scale` in `lib.rs` mirroring `set_secretion`.
- [ ] Rebuild: `uv run maturin develop --release`. Write `tests/test_cell_secretion_scale.py` (Python): one secreting cell, scale 0.5 → half concentration vs default; scale 0.0 → ~0; unrelated cells unaffected. Run → PASS.
- [ ] Also run `tests/test_fields.py` + `tests/test_influenza_fields.py` → still PASS (no regression). Commit `feat(engine): per-cell secretion scale for cell-state-regulated field sources`.

### Task 3.2: type-I IFN field
**Files:** `pbg_cpm_studies/influenza/fields.py` (+`add_ifn_field`), `tests/test_influenza_fields.py` (+tests).
**Interfaces:** `fields.add_ifn_field(world) -> int` — adds the IFN field (params D/decay), sets `types.I` secretion to the per-pixel `g_fp` rate (same z=1/cell_sites convention as the virus field; reuse the derivation helper), returns index. `fields.IFN_DIFFUSION_LENGTH_SITES = 10` (2 cell-diam × 5 sites).
- [ ] Failing tests mirroring the virus-field tests: infected cell raises IFN locally, decays with distance; no-infected → ~0. Note IFN diffuses FARTHER-per-unit than virus is FALSE — IFN length is 2 cell-diam (10 sites) vs virus 5 cell-diam (25 sites); assert the shorter IFN range is consistent. Implement, PASS, commit `feat(influenza): type-I IFN field + infected-cell secretion`.

### Task 3.3: per-cell resistance + (1−resist)-gated virus release
**Files:** `pbg_cpm_studies/influenza/resistance.py`, `pbg_cpm_studies/influenza/run.py` (extend the driver), `tests/test_influenza_resistance.py`, `tests/test_influenza_virus_infection.py` (extend).
**Interfaces:** `resistance.cell_resistance(ifn_at_cell, a_rf) -> float` = `ifn/(a_rf+ifn)` (pure; 0 at ifn=0; →1 as ifn→∞). A driver `run.run_virus_infection_with_ifn(...)` (or extend run_virus_infection with an `ifn=True` flag) that each update: advances fields (virus+IFN via world.step), reads IFN per cell, computes resist, calls `world.set_cell_secretion_scale(virus_fi, cid, 1.0 - resist)` for infected cells, then the infection transition as before.
- [ ] Failing tests: (a) `cell_resistance` pure-function values (0 at 0; 0.5 at ifn=a_rf; monotone). (b) INTEGRATION — the key claim: with IFN+resistance ON, a seeded lesion produces LESS total virus and FEWER infected cells at a fixed late step than the Increment-2 no-resistance run from the same seed (resistance slows spread). Assert `n_I_with_ifn(t) < n_I_without_ifn(t)` and/or `total_virus_with_ifn < total_virus_without` at a comparable step, same seed. Keep small/fast.
- [ ] Implement, PASS, commit `feat(influenza): per-cell resistance gates virus release (1-resist)`.

### Task 3.4: ifn-resistance study + viz + membership
**Files:** `workspace/studies/ifn-resistance/study.yaml` (hand-author, schema v3, phase Simulate), `workspace/investigations/influenza-sego2022/investigation.yaml` (+member), `pbg_cpm_studies/influenza/viz.py` (+figure), `tests/test_influenza_viz.py` (+test).
- [ ] Add `viz.ifn_resistance_figure(run_result)` (IFN field / resistance map + the with-vs-without-IFN infected-count comparison). Test returns a Figure (Agg).
- [ ] Hand-author `ifn-resistance/study.yaml`: report HONESTLY that IFN accumulates, per-cell resistance rises, and virus release + spread are reduced vs no-resistance (cite the actual measured with/without numbers from Task 3.3). Caveats: APC-driven IFN amplification omitted (global/Increment 8); resistance's death/Allee/IL-10 roles come in later increments; reproduction verdict PENDING (Increment 9). Wire `virus_infection` (or a new ifn composite) as the baseline; validate with lint-workspace.py; add as investigation member; keep investigation status honest.
- [ ] Run `-k influenza` (all pass). Commit `feat(influenza): ifn-resistance study + viz`.

## Self-Review notes
- Spec coverage: Increment 3 row (IFN field; per-cell ρ; §2.2; Fig 6 IFN) → Tasks 3.0–3.4. The gating EFFECT (resistance slows spread) is the falsifiable mechanism claim, tested in 3.3(b).
- First Rust change (3.1) — kept minimal + general + backward-compatible (default scale 1.0 preserves Increment-2 behavior); Rust unit test + Python binding test + no-regression run required.
- Resistance's other consumers (infected death, Allee, IL-10, contact-kill) are later increments; `resistance.cell_resistance` + the per-cell secretion-scale primitive are the reusable substrate.
- Deferred/known: `a_rf`/`g_fp` exact values pending Task 3.0; APC IFN amplification (Increment 8); the source's ContactKilling resist-sign oddity (transcription §7 discrepancy) is an Increment-7 concern.
