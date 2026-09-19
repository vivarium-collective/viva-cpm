# Influenza (Sego 2022) — Increment 9: Capstone quantitative reproduction

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Assemble the full multiscale model into one driver, and produce the three capstone reproduction studies (`repro-fig3b`, `repro-fig5-viral-load`, `repro-fig7-infection-fraction`) that compare ensemble observables against the Increment-0 digitized acceptance bands — the "reproducibility achieved" deliverable. Also close the two remaining §4a fidelity gaps from Increment 8.

**Architecture:** All mechanism primitives already exist across Increments 0–8 (infection/death transitions, Allee recovery, per-cell resistance, the four diffusive fields, macrophage/NK/CD8 chemotaxis + contact/nearby killing, ODE-driven recruitment, the hybrid global Price-2015 ODE). This increment (A) wires them into a single `run_full_model` per-MCS loop in the source's steppable order, (B) adds a scale-agnostic acceptance-band harness (`bands.py`), (C) builds the three `repro_*` drivers + studies. Reduced-scale (0.3 mm) ensembles run locally + in CI; the paper-scale (500×500, η=0.04, 50-replica) calibration runs are executed separately on the Mac mini (Phase B, post-merge — see the compute note).

**Tech Stack:** Python, numpy, scipy (ODE), matplotlib(Agg); the existing `pbg_cpm_studies/influenza/*` package + Rust CPM engine (no Rust change expected).

**Spec:** `docs/superpowers/specs/2026-09-18-influenza-sego2022-reproduction-design.md` (Increment 9 + §6 acceptance + §7 compute).

**Source authority:** `docs/cc3d-reference/sego2022-parameters.md` and `docs/cc3d-reference/sego2022-global-ode.md` (per-MCS steppable order §4). Targets: `pbg_cpm_studies/influenza/targets/{fig3b,fig5,fig7}.json` (Increment-0 digitized bands; `value`=ODE reference, `lo`/`hi`=50-replica spatial spread; `soft:true` observables are lenient).

## Global Constraints

- **This is the FIRST increment allowed to CALIBRATE toward the numeric targets.** But calibration must stay source-anchored: prefer running at the correct scale (η) where the source constants were tuned, and use the *already-documented* engine-unit scales (e.g. the Increment-7 100× chemotaxis normalization) — do NOT introduce new ad-hoc multipliers to force a band pass. Any calibration factor MUST be documented, justified against the source, and flagged; an undocumented fudge to pass a band is a review reject.
- **Honest bands.** A study `passes` only if the observable lands within `[lo,hi]` across the ensemble; otherwise it is `documented` with the miss quantified. `soft:true` observables (Type-I IFN + downstream, per fig3b notes) use the wider band. NEVER edit a target band to make a run pass.
- **Full-model per-MCS order (dossier §4 / source steppables):** world.step → viral internalization/infection (H→I ∝ local virus·g_hv) → update resistance ρ from local IFN → field secretion (virus by infected·(1−ρ); IFN by infected + APC; chemokine/IL-10 by macrophages with dynamic sig_1) → chemotaxis (macrophage↑virus, NK/CD8↑chemokine) → contact+nearby killing (I→D) → ROS/Allee death (H→D) + recovery (D→H) → infected apoptosis (I→D, mu_i·(1−ρ)) → recruitment → global ODE step (spatial→ODE→spatial). Reuse the existing per-increment helpers; do NOT reimplement their math.
- **Discrepancies #7,#9,#10,#11,#12** all still hold (resist-direct killing; hybrid ODE; spatial-string scalings; sig_1 saturating; asymmetric recruitment). Preserve.
- **CI feasibility:** the in-repo studies + tests run at reduced scale (0.3 mm, ~1225 cells) with a REDUCED replica count (≤5) and reduced steps — fast enough for CI. The paper-scale 50-replica runs are a documented Mac-mini follow-up, NOT run in CI. Perf: keep the full test suite feasible (no multi-minute unit test).
- **Reproduction status:** reduced-scale in-repo studies are `documented` with the band comparison shown honestly; the "reproduction achieved" verdict is asserted only after the Phase-B paper-scale ensemble lands in-band (recorded post-merge). Do NOT prematurely stamp `passing`/`reproduced` on a reduced-scale run.
- No AI attribution. Work ONLY in worktree `/Users/eranagmon/code/viva-cpm--influenza-incr9`; NEVER touch canonical `/Users/eranagmon/code/viva-cpm`. Do NOT edit `composites/__init__.py` or `visualizations/`. Per-worktree venv at `.venv`.

---

## File Structure

- `pbg_cpm_studies/influenza/bands.py` — NEW. Load a target, `series_in_band`, ensemble aggregation, `evaluate_study` (Task 9.0).
- `pbg_cpm_studies/influenza/run.py` — MODIFY. The two §4a fidelity fixes (Task 9.0) + `run_full_model` (Task 9.1) + `repro_fig3b`/`repro_fig5`/`repro_fig7` drivers (Tasks 9.2–9.4).
- `pbg_cpm_studies/influenza/viz.py` — MODIFY. `repro_fig3b_figure` etc. (Task 9.5).
- `pbg_cpm_studies/composites/influenza.py` — the modular biological composites (`epithelium`/`viral_infection`/`innate_immunity`/`cytotoxic_immunity`/`systemic_ode`/`full_model`) ALREADY EXIST (merged refactor #52). The capstone studies REFERENCE `full_model`; do NOT add a new composite.
- `tests/test_influenza_bands.py` — NEW (Task 9.0).
- `tests/test_influenza_full_model.py` — NEW (Tasks 9.1–9.4).
- `tests/test_influenza_viz.py` — MODIFY (Task 9.5).
- `workspace/studies/{repro-fig3b,repro-fig5-viral-load,repro-fig7-infection-fraction}/study.yaml` — NEW (Tasks 9.2–9.4).
- `workspace/investigations/influenza-sego2022/investigation.yaml` — MODIFY (Task 9.5).

---

### Task 9.0: acceptance-band harness + the two §4a fidelity fixes

**Files:** Create `pbg_cpm_studies/influenza/bands.py`, `tests/test_influenza_bands.py`; Modify `pbg_cpm_studies/influenza/run.py` (the `run_global_coupling` helpers so `run_full_model` inherits the fixes).

**Interfaces — Produces:**
- `bands.load(fig: str) -> dict` (wraps `targets.load_target`).
- `bands.series_in_band(series: list[tuple[float,float]], target_obs: list[dict], *, soft: bool=False) -> dict` — given a model time-series `[(t_days, value)]` and the target's per-observable band list, return `{"in_band": bool, "n_checked": int, "n_in": int, "worst_miss": float}` (interpolate model to each target t; `soft` widens the band by a documented factor, e.g. ×1.5).
- `bands.aggregate_replicas(runs: list[dict], key_path) -> list[tuple]` — ensemble mean series across replicas.
- `bands.evaluate_study(ensemble: dict, fig: str) -> dict` — per-observable in-band verdicts + an overall `passed` bool.

**§4a fidelity fixes (Increment-8 carry-forwards):**
1. **Nearby-surrogate feedback:** in the spatial→ODE input build, set ODE `M = local_M + M_nb`, `K = local_K + K_nb`, `E = local_E + E_nb` (dossier §4a: `num_immune_by_type` = local CPM + nearby surrogate), instead of the pure CPM count. (Macro `M_nb≡0` since local_ratio=1.0; NK/CD8 pick up their surrogates.)
2. **Recruited-macrophage secretion:** make the macrophage-id list used for chemokine/IL-10 secretion refresh each MCS to include newly-activated (recruited) macrophages, not a static pre-loop list.

- [ ] **Step 1: Failing tests** — `tests/test_influenza_bands.py`:

```python
from pbg_cpm_studies.influenza import bands

def test_series_in_band_basic():
    tgt = [{"t_days":0.0,"value":1000,"lo":900,"hi":1100},
           {"t_days":1.0,"value":300,"lo":150,"hi":500}]
    ok = bands.series_in_band([(0.0,1000),(1.0,300)], tgt)
    assert ok["in_band"] and ok["n_in"]==2 and ok["n_checked"]==2
    bad = bands.series_in_band([(0.0,1000),(1.0,50)], tgt)   # 50 < lo 150
    assert not bad["in_band"] and bad["worst_miss"] > 0

def test_soft_band_widens():
    tgt = [{"t_days":1.0,"value":300,"lo":250,"hi":350}]
    assert not bands.series_in_band([(1.0,360)], tgt)["in_band"]
    assert bands.series_in_band([(1.0,360)], tgt, soft=True)["in_band"]  # widened

def test_load_fig3b_has_observables():
    t = bands.load("fig3b")
    assert "uninfected_cells" in t["observables"]
```

Plus a fidelity-fix test in `tests/test_influenza_full_model.py` deferred to 9.1 (needs the driver).

- [ ] **Step 2: Run, verify fail.**
- [ ] **Step 3: Implement `bands.py`** (linear interpolation of the model series to each target t; `worst_miss` = max relative distance outside `[lo,hi]`; `soft` widens `[lo,hi]` about `value` by ×1.5, documented). Apply the two §4a fixes in `run.py`'s `run_global_coupling` input-build + secretion helpers.
- [ ] **Step 4: Run `tests/test_influenza_bands.py` + `tests/test_influenza_price_ode.py` (the §4a fixes must not break 8.x); verify pass.**
- [ ] **Step 5: Commit** — `feat(influenza): acceptance-band harness + close §4a fidelity gaps (Incr 9 Task 9.0)`

---

### Task 9.1: `run_full_model` — the complete multiscale driver

**Files:** Modify `pbg_cpm_studies/influenza/run.py`; Create/extend `tests/test_influenza_full_model.py`.

**Interfaces — Produces:**
- `run.run_full_model(*, cells_per_side: int, steps: int, seed: int, init_infection_frac: float | None = None, init_viral_load: float | None = None, s_per_mcs: float = 60.0, enable=("infection","ifn","death","allee","macrophage","chemokine","nk_cd8","killing","recruitment","ode"), **kw) -> dict` returning `{"mcs":[...], "t_days":[...], "counts":{"uninfected":[...],"infected":[...],"dead":[...],"macrophage":[...],"nk":[...],"cd8":[...]}, "fields":{"virus":[...],"ifn":[...],"chemo":[...],"il10":[...]}, "ode":{sp:[...]}}` — one full-model run. Each MCS runs the source-ordered pipeline (see Global Constraints). Either `init_infection_frac` (random infected fraction) or `init_viral_load` (uniform virus IC, no pre-infected cells) sets the scenario. `cells_per_side` sets η = (cells_per_side²)/250000.

- [ ] **Step 1: Failing tests** — `tests/test_influenza_full_model.py`:

```python
from pbg_cpm_studies.influenza import run

def test_full_model_smoke_all_observables():
    r = run.run_full_model(cells_per_side=15, steps=10, seed=1, init_infection_frac=0.05)
    n = len(r["mcs"])
    assert n == 10 and len(r["t_days"]) == n
    for k in ("uninfected","infected","dead","macrophage","nk","cd8"):
        assert len(r["counts"][k]) == n
    for f in ("virus","ifn","chemo","il10"):
        assert len(r["fields"][f]) == n
    assert set(("T","X","A","P")).issubset(r["ode"].keys())

def test_full_model_infection_depletes_uninfected():
    # With infection + death active and a 5% seed, uninfected count must fall
    r = run.run_full_model(cells_per_side=15, steps=25, seed=2, init_infection_frac=0.05)
    assert r["counts"]["uninfected"][-1] < r["counts"]["uninfected"][0]
    assert r["counts"]["infected"][0] > 0

def test_full_model_viral_load_scenario_seeds_infection():
    # init_viral_load with no pre-infected cells: virus drives H->I over time
    r = run.run_full_model(cells_per_side=15, steps=25, seed=3, init_viral_load=1000.0)
    assert r["counts"]["infected"][0] == 0
    assert max(r["counts"]["infected"]) > 0     # infection emerges from the virus field

def test_full_model_fidelity_fixes():
    # §4a: ODE M/K/E inputs include nearby surrogates (K input >= local K when K_nb>0);
    # recruited macrophages secrete (macrophage id list refreshes). Assert via the driver's
    # recorded ode K vs spatial macrophage/nk counts, or an exposed debug hook. Keep FAST.
    r = run.run_full_model(cells_per_side=15, steps=8, seed=4, init_infection_frac=0.05)
    assert "K" in r["ode"]
```

- [ ] **Step 2: Run, verify fail.**
- [ ] **Step 3: Implement `run_full_model`.** Build the scenario (epithelial sheet at `cells_per_side`; seed initial infection by fraction OR uniform virus IC; seed homeostatic macrophages/NK/CD8 near-lesion + the recruit reserve pool). Each MCS execute the ordered pipeline reusing: `transitions.infection_step`, `resistance.cell_resistance`, `fields.*` secretion + `world.set_cell_secretion_scale` with dynamic `sig_1`, `immune.set_macrophage_chemotaxis`/`set_nk_cd8_chemotaxis` (apply the documented Incr-7 100× normalization at paper-relevant scale, cited), `killing.contact_kill_rate`+`nearby_kill_rate`, `allee.allee_death_rate`/`allee_recovery_rate`, `transitions.infected_death_step`, `recruitment.*`, `price_ode.GlobalODE.step`. `t_days = mcs*s_per_mcs/86400`. Record all observables. Respect `enable` flags (so tests can isolate). Reuse `run_global_coupling`'s ODE-coupling block (with the 9.0 fidelity fixes) rather than duplicating it.
- [ ] **Step 4: Run `tests/test_influenza_full_model.py`; verify pass.** Keep each test FAST (small `cells_per_side`, few steps).
- [ ] **Step 5: Commit** — `feat(influenza): run_full_model complete multiscale driver (Incr 9 Task 9.1)`

---

### Task 9.2: `repro_fig3b` driver + study

**Files:** Modify `run.py` (add `repro_fig3b`), Create `workspace/studies/repro-fig3b/study.yaml`, extend `tests/test_influenza_full_model.py`.

**Interfaces — Produces:** `run.repro_fig3b(*, replicas: int = 3, cells_per_side: int = 35, steps: int = 240, seed0: int = 0) -> dict` — runs the fig3b scenario (5% initial infection, 0.3 mm ≈ 35×35 cells → 1225) for `replicas` seeds, aggregates the ensemble mean of each observable, and calls `bands.evaluate_study(ensemble, "fig3b")`. Returns `{"ensemble":{obs:series}, "band_eval":{...}, "replicas":N, "cells_per_side":...}`. Observable names map to the fig3b target keys (uninfected_cells, infected_cells, dead_cells, virus, type1_ifn, macrophages, nk, cd8 — match the JSON's actual keys).

- [ ] **Step 1: Failing test** — in `tests/test_influenza_full_model.py`:

```python
def test_repro_fig3b_runs_and_evaluates_bands():
    r = run.repro_fig3b(replicas=2, cells_per_side=15, steps=20, seed0=0)  # reduced for CI
    assert r["replicas"] == 2
    assert "uninfected_cells" in r["ensemble"]
    assert "band_eval" in r and "passed" in r["band_eval"]
    # non-vacuous: the ensemble uninfected series declines under infection
    u = r["ensemble"]["uninfected_cells"]
    assert u[-1][1] <= u[0][1]
```

- [ ] **Step 2: Run, verify fail.**
- [ ] **Step 3: Implement `repro_fig3b`** (loop seeds → `run_full_model(init_infection_frac=0.05, ...)` → map counts/fields to the fig3b observable keys, ensemble-mean via `bands.aggregate_replicas` → `bands.evaluate_study`). The test uses a tiny reduced config; the docstring documents the paper config (`replicas=50, cells_per_side=35, steps≈2880` = 2 days at 1 min/MCS) for the Mac-mini Phase-B run.
- [ ] **Step 4: Run test; verify pass.**
- [ ] **Step 5: Hand-author `workspace/studies/repro-fig3b/study.yaml`** — schema_version 3, investigation influenza-sego2022, phase Evaluate, status in-progress, verdict `documented`, biological_validation PENDING. Report the REDUCED-scale band comparison honestly (which observables land in-band at reduced scale, which miss and by how much) and state the paper-scale 50-replica ensemble is the Mac-mini follow-up that will set the `reproduced` verdict. Target: Fig 3B. Cite the band source.
- [ ] **Step 6: Lint (`scripts/lint-workspace.py`) + commit** — `feat(influenza): repro-fig3b capstone study (Incr 9 Task 9.2)`

---

### Task 9.3: `repro_fig5` (viral-load sweep) driver + study

**Files:** Modify `run.py`, Create `workspace/studies/repro-fig5-viral-load/study.yaml`, extend the test file.

**Interfaces — Produces:** `run.repro_fig5(*, loads=(1,10,100,1000,10000), replicas=3, cells_per_side=35, steps=240, seed0=0) -> dict` — for each initial viral load, run `run_full_model(init_viral_load=load, init_infection_frac=None)` × replicas, extract the fig5 observable (e.g. final/min uninfected fraction = survival, and the lethal-threshold load), compare to `fig5.json` bands. Returns `{"by_load":{load:{...}}, "lethal_threshold":..., "band_eval":{...}}`.

- [ ] **Step 1: Failing test:**

```python
def test_repro_fig5_viral_load_sweep():
    r = run.repro_fig5(loads=(1,10000), replicas=1, cells_per_side=12, steps=15, seed0=0)  # reduced
    assert set(r["by_load"].keys()) == {1,10000}
    # higher initial viral load => not more surviving uninfected than the low-load case
    surv = lambda L: r["by_load"][L]["uninfected_final_frac"]
    assert surv(10000) <= surv(1) + 1e-9
    assert "band_eval" in r
```

- [ ] **Step 2–4:** verify fail → implement (read `fig5.json` for the exact observable + bands; the monotone dose-response is the fidelity criterion) → verify pass (reduced config; document the paper config in the docstring).
- [ ] **Step 5: Hand-author `workspace/studies/repro-fig5-viral-load/study.yaml`** (as 9.2, honest, `documented`/PENDING, paper-scale = mini follow-up). Target: Fig 5.
- [ ] **Step 6: Lint + commit** — `feat(influenza): repro-fig5 viral-load sweep study (Incr 9 Task 9.3)`

---

### Task 9.4: `repro_fig7` (infection-fraction sweep) driver + study

**Files:** Modify `run.py`, Create `workspace/studies/repro-fig7-infection-fraction/study.yaml`, extend the test file.

**Interfaces — Produces:** `run.repro_fig7(*, fracs=(0.001,0.005,0.01,0.05), replicas=3, cells_per_side=35, steps=240, seed0=0) -> dict` — for each initial infection fraction, run × replicas, extract the fig7 observable (per `fig7.json`), compare to bands. Returns `{"by_frac":{...}, "band_eval":{...}}`.

- [ ] **Step 1: Failing test:**

```python
def test_repro_fig7_infection_fraction_sweep():
    r = run.repro_fig7(fracs=(0.001,0.05), replicas=1, cells_per_side=12, steps=15, seed0=0)
    assert set(r["by_frac"].keys()) == {0.001,0.05}
    # a larger initial infection fraction does not leave more uninfected at the end
    f = lambda x: r["by_frac"][x]["uninfected_final_frac"]
    assert f(0.05) <= f(0.001) + 1e-9
    assert "band_eval" in r
```

- [ ] **Step 2–4:** verify fail → implement (read `fig7.json` for the observable + bands) → verify pass (reduced; document paper config).
- [ ] **Step 5: Hand-author `workspace/studies/repro-fig7-infection-fraction/study.yaml`** (honest, `documented`/PENDING). Target: Fig 7.
- [ ] **Step 6: Lint + commit** — `feat(influenza): repro-fig7 infection-fraction sweep study (Incr 9 Task 9.4)`

---

### Task 9.5: capstone viz + investigation summary (references the existing `full_model` composite)

**Files:** Modify `viz.py`, `workspace/investigations/influenza-sego2022/investigation.yaml`, `tests/test_influenza_viz.py`. Do NOT add a composite — `full_model` already exists (refactor #52); the three repro studies reference it in their `base_model`/condition.

**Interfaces — Produces:** `viz.repro_fig3b_figure(result)`, `viz.repro_sweep_figure(result, kind)` (Agg; observable ensemble vs target band); investigation.yaml gains the three repro studies as members + a capstone `at_a_glance`/`acceptance_criteria` entry summarizing the reduced-scale band results and the pending paper-scale verdict. The three study.yaml `base_model` refs point at `pbg_cpm_studies.composites.influenza.full_model`.

- [ ] **Step 1: Failing viz test** — build a tiny `repro_fig3b(replicas=1, cells_per_side=12, steps=8)` result, assert `repro_fig3b_figure(result)` returns a Figure with ≥2 axes (ensemble series + band overlay). Keep FAST.
- [ ] **Step 2: Run, verify fail.**
- [ ] **Step 3: Implement** the figures (overlay ensemble mean + `[lo,hi]` band from the target) and add the three studies to `investigation.yaml` (members + at_a_glance + acceptance_criteria), honestly noting reduced-scale band status + the paper-scale mini follow-up. Studies reference the existing `full_model` composite; NO new composite is added.
- [ ] **Step 4: Run `tests/test_influenza_viz.py` + `scripts/lint-workspace.py`; verify pass/OK** (new studies may be dashboard-✗ PENDING — expected).
- [ ] **Step 5: Commit** — `feat(influenza): capstone viz + investigation summary (Incr 9 Task 9.5)`

---

## Compute note (Phase B — post-merge, Mac mini)

After this branch merges, the controller syncs the repo to the Mac mini and launches the paper-scale ensembles (`repro_fig3b(replicas=50, cells_per_side=35, steps≈2880)`, `repro_fig5`/`repro_fig7` at paper config) via headless `mct` task agents (one worktree each), verifying via commits. The paper-scale band results update each study's verdict to `reproduced` (or `documented` with the quantified miss + a calibration follow-up). The in-repo tests/studies remain at reduced scale for CI. This split is per spec §7 (capstone compute budgeted to the mini).

## Self-Review notes (author)

- **Spec coverage:** Increment 9 row = 3 repro studies + full composite + paper-scale ensemble + acceptance bands → Tasks 9.1 (full driver), 9.2/9.3/9.4 (the three studies), 9.0 (band harness + §4a fidelity), 9.5 (viz/composite/summary), Phase B (paper-scale). Covered.
- **Calibration honesty is the load-bearing constraint:** no band edited, no ad-hoc multiplier; reduced-scale studies are `documented`, paper-scale verdict deferred to Phase B. A reviewer must reject any band edit or undocumented fudge factor.
- **Type consistency:** `run_full_model` result keys → `repro_*` observable maps → `bands.evaluate_study` → viz, all fixed in 9.1/9.0 and reused downstream. Observable key names come from the fig JSONs (implementers read them).
- **CI feasibility:** every in-repo test uses a tiny reduced config; the paper config lives only in docstrings + Phase B.
