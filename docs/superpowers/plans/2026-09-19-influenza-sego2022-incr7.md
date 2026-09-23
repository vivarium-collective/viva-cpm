# Influenza (Sego 2022) reproduction — Increment 7 Implementation Plan (NK + CD8⁺ T cells)

> **For agentic workers:** REQUIRED SUB-SKILL: superpowers:subagent-driven-development. Steps use `- [ ]`.

**Goal:** Add the cytotoxic local immune cells — NK cells (K=5) and CD8⁺ T cells (E=6) — recruited (stubbed) to the sheet, chemotaxing up the chemokine gradient (from Increment 6) toward the infection, and KILLING infected cells on contact — so the infected-cell count is reduced vs a run without them (the effector arm of the immune response, Fig 2/3A).

**Architecture:** Reuse existing primitives — NK/CD8 as cell types with source-cited adhesion; chemotaxis up the chemokine field via `set_chemotaxis(chemo_fi, K/E, λ)`; contact-mediated killing via the Increment-4 `cell_contact_area_by_type` query (to get each infected cell's NK/CD8 contact area) + a Python killing Step. NO Rust change. Recruitment stubbed (chemokine-Hill, like macrophages). Haptotaxis DEFERRED (a secondary migration bias + a new engine capability — note the gap, Increment 9).

**Tech Stack:** Python 3.12, `cpm.cpm_core`, pytest, numpy, PyYAML. Builds on merged Increments 0–6.

**Spec:** `docs/superpowers/specs/2026-09-18-influenza-sego2022-reproduction-design.md` (Increment 7 row of §5).

## Global Constraints (source: ContactKillingSteppable + ChemotaxisSteppable)
- Cell types: NK = K = 5, CD8⁺ = E = 6 (`types.py`). Chemotaxis up the CHEMOKINE field: NK λ = 5000 (`chemotaxis_v_nk`), CD8 λ = 10000 (`chemotaxis_v_cd8`).
- **Local contact-killing** (source-literal, `ContactKillingSteppable.step`): for each infected cell,
  `srf_nk` = contact area with NK neighbors, `srf_cd8` = with CD8; `kill_rate_nk = (g_ik·tot_ec_ODE)·srf_nk·cell_resist/cell.volume`, `Pr = 1−exp(−kill_rate)` → cell becomes D (DYING). Same with `g_ie` for CD8.
  **DISCREPANCY #7 (source-faithful):** the source multiplies by `cell_resist` DIRECTLY (NOT `1−resist`) — unlike every other resistance consumer. Implement `cell_resist` directly per the source; FLAG the quirk in code + report (a possible intentional "resistant cells present more antigen" effect, or a source inconsistency — do not silently "fix" it).
  The "nearby population" killing term (`(g_ik/pop_scale)·num_nk·resist`) is a well-mixed/population approximation — DEFER (it needs a global nearby-count; the spatially-faithful part is the local surface-contact term).
- Differential adhesion (Table 3 / source Contact rows): homotypic-immune 25, heterotypic-immune 10, uninfected-immune 20, infected-immune 10, dead-immune 20 — use the SOURCE's literal Contact rows for NK/CD8 (Task 7.0). "Preferential attachment to infected" is realized by the J values (infected-immune 10 < homotypic 25).
- Recruitment STUBBED (source: Hill on chemokines/APCs — Increment 8) — seed NK/CD8 at a simple rate, like the macrophage scenario.
- 2D-approximation geometry (source z=2 immune layer → Increment 9), same as Increment 5. params.yaml single authority. No AI attribution. Work only in worktree `~/code/viva-cpm--influenza-incr7` (branch `investigation/influenza-incr7`); never touch `~/code/viva-cpm`, `composites/__init__.py`, or `visualizations/`. Tests: `~/code/viva-cpm--influenza-incr7/.venv/bin/python -m pytest <path> -v` (TARGETED runs; full `-k influenza` is slow).

---

### Task 7.0: Extract & cite NK + CD8 constants
**Files:** `params.yaml` (+`nk:`/`cd8:` sections), `docs/cc3d-reference/sego2022-parameters.md`, `tests/test_influenza_params.py` (+test).
- [ ] Re-fetch source. Read `chemotaxis_v_nk` (5000), `chemotaxis_v_cd8` (10000), `g_ik` (NK killing), `g_ie` (CD8 killing), `tot_ec_ODE`, the NK/CD8 adhesion Contact rows (ViralInfectionVTM.xml), NK/CD8 volume, and the recruitment Hill params (record-but-stub). Cite file+symbol+units + conversions. Note the local-vs-nearby killing split (`pr_cf_*_loc = g_i*·tot_ec`; nearby deferred) and DISCREPANCY #7 (resist-direct).
- [ ] Add `nk:`/`cd8:` sections. Test asserts chemotaxis λ (5000/10000) + g_ik/g_ie + adhesion + source.
- [ ] Commit `feat(influenza): cite NK + CD8 constants (Incr 7)`.

### Task 7.1: NK + CD8 cell types + chemotaxis up chemokines
**Files:** `viva_cpm_studies/influenza/immune.py` (+NK/CD8 seeding + `set_nk_cd8_chemotaxis`), `viva_cpm_studies/influenza/run.py` (extend the driver), `tests/test_influenza_nk_cd8.py`.
**Interfaces:** seed NK (K) + CD8 (E) cells with source-cited adhesion; `set_chemotaxis(chemo_fi, K, 5000)`, `set_chemotaxis(chemo_fi, E, 10000)`. Extend the signaling/macrophage scenario so macrophages secrete chemokine (Incr 6) and NK/CD8 chemotax up it.
- [ ] Failing INTEGRATION test (reuse the Incr-5 interior/multi-seed rigor pattern): NK and CD8 cells, seeded away from the infection, reduce their mean distance to the infection over the run (they follow the chemokine gradient), with a λ=0 control that does not — across multiple seeds (unbiased interior control, like Increment 5's fix). CD8 (λ=10000) localizes at least as strongly as NK (λ=5000).
- [ ] Implement, PASS. Commit `feat(influenza): NK + CD8 cells + chemotaxis up the chemokine gradient`.

### Task 7.2: contact-mediated killing of infected cells
**Files:** `viva_cpm_studies/influenza/killing.py` (new; the kill-rate + the killing Step logic), `viva_cpm_studies/influenza/run.py` (extend the driver to apply killing each update), `tests/test_influenza_killing.py`.
**Interfaces:** pure `killing.contact_kill_rate(srf_immune, cell_resist, g_i, tot_ec, cell_volume) -> float` = `g_i*tot_ec*srf_immune*cell_resist/cell_volume` (resist DIRECT, per source #7). A driver step: for each infected cell, get `cell_contact_area_by_type(cid)` → srf_nk = area with K, srf_cd8 = area with E; compute NK + CD8 kill rates; `Pr = 1−exp(−rate)`; on success set the infected cell to D. (Combine with the existing infection/death/Allee fate driver.)
- [ ] Failing tests: (a) pure rate function (monotone in srf_immune + resist; 0 when srf_immune=0). (b) INTEGRATION: with NK/CD8 present + chemotaxing to the infection, the infected-cell count (and/or peak infection) is LOWER at a comparable step than a no-NK/CD8 run from the same seed — cytotoxic clearance. Deterministic/seeded.
- [ ] Implement (flag discrepancy #7 in code + report), PASS. Commit `feat(influenza): NK/CD8 contact-mediated killing of infected cells`.

### Task 7.3: cytotoxic-killing study + viz + membership
**Files:** `workspace/studies/cytotoxic-killing/study.yaml` (hand-author), `workspace/investigations/influenza-sego2022/investigation.yaml` (+member), `viva_cpm_studies/influenza/viz.py` (+figure), `tests/test_influenza_viz.py` (+test). Do NOT touch visualizations/ or composites/__init__.py.
- [ ] Add `viz.cytotoxic_killing_figure(run_result)` (NK/CD8 localization + infected-count with vs without NK/CD8). Test returns a Figure (Agg).
- [ ] Hand-author `cytotoxic-killing/study.yaml`: report HONESTLY the NK/CD8 localization + the infected-cell reduction (cite actual Task 7.1/7.2 numbers). Caveats: nearby-population killing term deferred; haptotaxis DEFERRED (secondary bias + new engine capability, Incr 9); recruitment stubbed; 2D-approx geometry; DISCREPANCY #7 (resist-direct killing) flagged; chemotaxis linear-form gap; reproduction PENDING (Incr 9). verdict = documented. Wire a baseline composite (reuse/add in composites/influenza.py — auto-registers; do NOT edit __init__.py). Validate lint; add member.
- [ ] Run targeted influenza tests (all pass). Commit `feat(influenza): cytotoxic-killing study + viz`.

## Self-Review notes
- Spec coverage: the "NK & CD8⁺ T cells" row → Tasks 7.0–7.3. Contact killing uses the Increment-4 contact-area query; chemotaxis uses the Increment-6 chemokine field.
- No Rust change: chemotaxis + cell_contact_area_by_type + a killing Step all exist. Confirm in 7.1/7.2.
- Deferred/flagged: haptotaxis (Incr 9), nearby-population killing (well-mixed → Incr 8), recruitment realism (Incr 8), discrepancy #7 (resist-direct — source-faithful, flagged), differential-adhesion exactness, 2D-vs-z=2 geometry.
- Reuse Increment 5's interior/multi-seed localization-test rigor (unbiased control, multi-seed) — don't repeat the single-seed/corner-control mistake.
