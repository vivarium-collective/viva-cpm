# Influenza (Sego 2022) reproduction — Increment 6 Implementation Plan (chemokine + IL-10 fields)

> **For agentic workers:** REQUIRED SUB-SKILL: superpowers:subagent-driven-development. Steps use `- [ ]`.

**Goal:** Add the macrophage-released chemokine and IL-10 diffusive fields, so the chemokine field forms a gradient centered on the macrophage cluster at the infection (the attractant NK/CD8⁺ will follow in Increment 7), and IL-10 accumulates. This is the field half of the immune module (the sources — macrophages — arrived in Increment 5).

**Architecture:** Two new diffusive fields via `add_field` (existing). Per-cell secretion regulated by local IL-10 (a Hill term) and by resistance — implemented with the Increment-3 per-cell secretion-scale primitive (`set_cell_secretion_scale`): a per-type base rate × a per-cell scale computed each update from the cell's local IL-10 / resistance. NO Rust change. The uniform IL-10 initial condition is SKIPPED (no field-write API; IL-10 builds from sources — a documented Increment-9 detail).

**Tech Stack:** Python 3.12, `cpm.cpm_core`, pytest, numpy, PyYAML. Builds on merged Increments 0–5.

**Spec:** `docs/superpowers/specs/2026-09-18-influenza-sego2022-reproduction-design.md` (Increment 5 row "chemokine + IL-10 fields" — brought here after the 5/6 re-scope).

## Global Constraints (source forms — CC3D `ChemokineSecretionSteppable` / `IL10SecretionSteppable`)
- **Chemokine field** `chemo`: unitless D ≈ 15.601 lat²/MCS, decay ≈ 6.240e-3 /MCS, diffusion length 10 cell-diam (`chemo_dc`/`chemo_decay`). SOURCE = macrophages, IL-10-Hill regulated:
  `sec = b_c · sig_1 / (sig_1 + (g_1·L_loc + g_2)/(L_loc + d_2))`, where L_loc = local mean IL-10 in the cell. (b_c = source `b_c`·dim.z.) Plus a global boundary term (num_macro-based) — DEFER (global/boundary, Increment 8).
- **IL-10 field** `il10`: unitless D ≈ 4.901 lat²/MCS, decay ≈ 1.960e-3 /MCS, diffusion length 10 cell-diam (`il10_dc`/`il10_decay`). TWO local sources:
  - Macrophages (same IL-10-Hill form with `b_l`): `sec = b_l · sig_1/(sig_1 + (g_1·L_loc + g_2)/(L_loc + d_2))`.
  - **Uninfected (H) cells**, resistance-scaled: `sec = mu_l · b_lh · (1 − resist)`.
  - Global boundary term — DEFER. Uniform IC `L = b_lh·(1−mean_resist)·num_epithelial/(dim.x·dim.y)` at t=0 — SKIP (no field-write API; document).
- Same z=1/per-pixel/cell_sites secretion convention as prior fields (reuse `_per_pixel_secretion_rate`). Per-cell regulation via `set_cell_secretion_scale`.
- params.yaml single authority. No AI attribution. Work only in worktree `~/code/viva-cpm--influenza-incr6` (branch `investigation/influenza-incr6`); never touch `~/code/viva-cpm`, `composites/__init__.py`, or `visualizations/`. Tests: `~/code/viva-cpm--influenza-incr6/.venv/bin/python -m pytest <path> -v`. Prefer TARGETED test runs (the full `-k influenza` suite is ~140s+).

---

### Task 6.0: Extract & cite chemokine + IL-10 field + secretion constants
**Files:** `params.yaml` (+`chemokine:` and `il10:` sections), `docs/cc3d-reference/sego2022-parameters.md`, `tests/test_influenza_params.py` (+test).
- [ ] Re-fetch source. From `ImmuneModelInputs.py`/`ImmuneModelLib.py` + `ChemokineSecretionSteppable`/`IL10SecretionSteppable`, read: chemo/il10 unitless D+decay (`chemo_dc`/`chemo_decay`/`il10_dc`/`il10_decay`), diffusion lengths (10 each), and the secretion + regulation constants `b_c`, `b_l`, `b_lh`, `mu_l`, `sig_1` (Sigma1), `g_1`, `g_2`, `d_2`. Cite file+symbol+units + per-day→per-MCS + dim.z conventions.
- [ ] Add `chemokine:` (D, decay, length=10, b_c, source) and `il10:` (D, decay, length=10, b_l, b_lh, mu_l, sig_1, g_1, g_2, d_2, source). Test asserts keys + numeric + source + lengths==10.
- [ ] Commit `feat(influenza): cite chemokine + IL-10 field constants (Incr 6)`.

### Task 6.1: chemokine + IL-10 fields with regulated macrophage/uninfected sources
**Files:** `viva_cpm_studies/influenza/fields.py` (+`add_chemokine_field`, `add_il10_field`), `viva_cpm_studies/influenza/signaling.py` (new; the per-cell Hill/resist regulation helpers), `viva_cpm_studies/influenza/run.py` (extend the driver to add the fields + apply per-cell secretion scales each update), `tests/test_influenza_signaling.py`.
**Interfaces:**
- `fields.add_chemokine_field(world) -> int`, `fields.add_il10_field(world) -> int` — add the fields (params D/decay), set per-type base secretion (`set_secretion(chemo_fi, M, b_c_base)`; `set_secretion(il10_fi, M, b_l_base)` and `set_secretion(il10_fi, H, mu_l*b_lh_base)`), stability via `set_field_dynamics`.
- `signaling.macrophage_secretion_scale(il10_local, sig_1, g_1, g_2, d_2) -> float` = `sig_1/(sig_1 + (g_1*il10_local + g_2)/(il10_local + d_2))` (the IL-10-Hill factor; pure). `signaling.uninfected_il10_scale(resist) -> float` = `1 - resist` (pure).
- Driver: each update, after world.step advances all fields, for each macrophage set `set_cell_secretion_scale(chemo_fi, cid, hill)` and `set_cell_secretion_scale(il10_fi, cid, hill)` from its local IL-10; for each H cell set `set_cell_secretion_scale(il10_fi, cid, 1-resist)`.
- [ ] Failing tests: (a) pure regulation functions (Hill: monotone, bounded in (0,1)-ish; uninfected scale 0..1). (b) INTEGRATION: in a macrophage scenario (reuse Increment 5's), after running, the CHEMOKINE field is POSITIVE and highest near the macrophage cluster and decays with distance (a gradient centered on the macrophages); IL-10 is positive (from macrophages + uninfected cells). Deterministic/seeded, small/fast.
- [ ] Implement, PASS. Commit `feat(influenza): macrophage-released chemokine + IL-10 fields (IL-10-Hill + resist regulated)`.

### Task 6.2: signaling-fields study + viz + membership
**Files:** `workspace/studies/signaling-fields/study.yaml` (hand-author), `workspace/investigations/influenza-sego2022/investigation.yaml` (+member), `viva_cpm_studies/influenza/viz.py` (+figure), `tests/test_influenza_viz.py` (+test). Do NOT touch visualizations/ or composites/__init__.py.
- [ ] Add `viz.signaling_fields_figure(run_result)` (chemokine + IL-10 field heatmaps / radial profiles centered on the macrophage cluster). Test returns a Figure (Agg).
- [ ] Hand-author `signaling-fields/study.yaml`: report HONESTLY that the chemokine field forms a gradient centered on the macrophage cluster (cite the actual Task 6.1 numbers — peak-near-macrophages, decay-with-distance) and IL-10 accumulates. Caveats: global/boundary chemokine+IL-10 terms deferred (Incr 8); IL-10 uniform IC skipped (no field-write API; builds from sources — Incr-9 detail); recruitment still stubbed; reproduction PENDING (Incr 9). verdict = documented. Wire a baseline composite (reuse macrophage_response or add a signaling variant in composites/influenza.py — auto-registers; do NOT edit __init__.py). Validate with lint-workspace.py; add member.
- [ ] Run targeted influenza tests (all pass). Commit `feat(influenza): signaling-fields study + viz`.

## Self-Review notes
- Spec coverage: the "chemokine + IL-10 fields" row → Tasks 6.0–6.2 (macrophage-released, per the 5/6 re-scope).
- No Rust change: uses add_field + set_secretion + the Increment-3 set_cell_secretion_scale (per-cell regulation). Confirm in 6.1.
- Deferred/documented: global boundary secretion terms (Incr 8), IL-10 uniform IC (no field-write API — Incr 9 or a future small Rust field-set), the self-regulation loop's exact quantitative behavior (Incr-9 calibration).
- Sets up Increment 7: the chemokine gradient centered on the infection is the attractant for NK/CD8⁺ chemotaxis + recruitment.
