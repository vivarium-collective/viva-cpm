# Receptor Recruitment — Increment 1 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a receptor-level realization of chemotactic recruitment to viva-cpm — a per-cell ligand–receptor occupancy model that drives an *emergent* chemotactic response — with a cited chemokine Kd and multi-seed rigor.

**Architecture:** A new `ReceptorSubcell` Process reads each responder's local chemokine concentration (`field_at_cell`), computes Hill receptor occupancy, and writes a responder *fate* (naive/activated cell-type) via the CPM `fates` port. A new composite wires one `ReceptorSubcell` per responder to a single `CPMProcess`; activated responders carry chemotaxis λ, naive do not. Two studies (baseline, receptor-blocked) run multi-seed; the receptor Kd is a cited `calibration_anchor`. Two figures follow.

**Tech Stack:** Python (process-bigraph, viva-superpowers), the compiled Rust CPM engine (`cpm.cpm_core`, unchanged), pytest, Plotly (baked-data viz).

## Global Constraints

- **Zero Rust / engine changes.** Reuse the per-type λ + `fates` (`set_cell_type`) mechanism as-is.
- **Terminology:** "translator" is the project term for the typed translation artifact (a.k.a. the Domain Bridges design's "domain bridge"); do not conflate with the process-bigraph **Composite bridge** (structural wiring).
- **H4 pre-init fallback.** Until the bigraph-schema keyed-map affordance (H5, separate plan) lands, the composite MUST pre-initialize the `fates` store with a key per wired responder cell — an absent map key is dropped on apply (`bigraph_schema/methods/apply.py:438-439`).
- **Contact-J invariant.** The activated responder sub-type (3) MUST share the naive responder's (2) contact energies J on every pair, so only chemotaxis differs (no adhesion confound).
- **Rigor:** ≥5 seeds with reported CIs on the recruitment index; every acceptance band on a receptor study carries `calibration_anchor`/`cites`.
- **Worktree:** All work happens in `~/code/viva-cpm--receptor-recruitment` (branch `receptor-recruitment`). Do not touch the canonical `~/code/viva-cpm`.
- **Run tests:** `pytest` from the worktree root (`[tool.pytest.ini_options] pythonpath=["."]`). The compiled engine `cpm/cpm_core*.so` is already built in the checkout.

**Scope note:** This plan is Increment 1 only. The bigraph-schema keyed-map affordance (H5), the dose-response sweep + evidence-infra H2 (Increment 2), and the scale-aggregate H1 + refinement check H3 (Increment 3) each get their own plan. The dose-response figure (spec §8 Figure 1) requires the cue-rate sweep and therefore lands in Increment 2; Increment 1 ships the spatial-activation and condition-overlay figures.

---

### Task 1: `ReceptorSubcell` process

**Files:**
- Create: `cpm/subcellular/receptor.py`
- Test: `tests/test_receptor_subcell.py`

**Interfaces:**
- Produces: `ReceptorSubcell(config: dict, core)` — a `process_bigraph.Process`. `inputs()={"ligand":"float"}`, `outputs()={"fate":"overwrite[integer]"}`. Method `occupancy(ligand: float) -> float` (fractional, in [0,1]). `update(state, interval) -> {"fate": int}`. Config keys: `kd, hill, conc_scale, activate_occupancy, naive_type, activated_type`.

- [ ] **Step 1: Write the failing test** — `tests/test_receptor_subcell.py`

```python
import process_bigraph as pb
from cpm.subcellular.receptor import ReceptorSubcell

CFG = {"kd": 10.0, "hill": 1.0, "conc_scale": 1.0, "activate_occupancy": 0.5,
       "naive_type": 2, "activated_type": 3}


def _proc(**over):
    cfg = {**CFG, **over}
    return ReceptorSubcell(cfg, core=pb.allocate_core())


def test_occupancy_half_at_kd():
    assert abs(_proc(kd=10.0).occupancy(10.0) - 0.5) < 1e-9


def test_occupancy_monotonic_and_bounded():
    p = _proc(kd=10.0)
    vals = [p.occupancy(c) for c in (0.0, 1.0, 5.0, 10.0, 50.0, 500.0)]
    assert vals[0] == 0.0
    assert all(0.0 <= v <= 1.0 for v in vals)
    assert all(b >= a for a, b in zip(vals, vals[1:]))


def test_fate_activates_above_threshold():
    p = _proc(kd=10.0, activate_occupancy=0.5)
    assert p.update({"ligand": 20.0}, 1.0)["fate"] == 3   # theta ~0.667 -> activated
    assert p.update({"ligand": 2.0}, 1.0)["fate"] == 2    # theta ~0.167 -> naive
    assert p.update({"ligand": 0.0}, 1.0)["fate"] == 2


def test_deterministic():
    p = _proc()
    assert p.update({"ligand": 12.3}, 1.0) == p.update({"ligand": 12.3}, 1.0)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_receptor_subcell.py -q`
Expected: FAIL — `ModuleNotFoundError: No module named 'cpm.subcellular.receptor'`.

- [ ] **Step 3: Write minimal implementation** — `cpm/subcellular/receptor.py`

```python
from process_bigraph import Process


class ReceptorSubcell(Process):
    """Minimal per-cell ligand-receptor occupancy model. Reads the local
    chemokine concentration, computes fractional receptor occupancy via a Hill
    relation, and emits a responder cell-type (fate): activated when occupancy
    crosses threshold, else naive. Deterministic.

        theta = c^h / (Kd^h + c^h),   c = ligand * conc_scale

    With hill=1, conc_scale=1 and activate_occupancy=0.5 the activation threshold
    sits exactly at c = Kd, so the cited Kd is directly visible in a dose-response.
    """

    config_schema = {
        "kd": {"_type": "float", "_default": 1.0},
        "hill": {"_type": "float", "_default": 1.0},
        "conc_scale": {"_type": "float", "_default": 1.0},
        "activate_occupancy": {"_type": "float", "_default": 0.5},
        "naive_type": {"_type": "integer", "_default": 2},
        "activated_type": {"_type": "integer", "_default": 3},
    }

    def initialize(self, config):
        self.kd = float(config["kd"])
        self.hill = float(config["hill"])
        self.conc_scale = float(config["conc_scale"])
        self.activate_occupancy = float(config["activate_occupancy"])
        self.naive_type = int(config["naive_type"])
        self.activated_type = int(config["activated_type"])

    def inputs(self):
        return {"ligand": "float"}

    def outputs(self):
        return {"fate": "overwrite[integer]"}

    def occupancy(self, ligand):
        c = max(0.0, float(ligand)) * self.conc_scale
        if c <= 0.0:
            return 0.0
        ch = c ** self.hill
        return ch / (self.kd ** self.hill + ch)

    def update(self, state, interval):
        theta = self.occupancy(float((state or {}).get("ligand", 0.0)))
        fate = self.activated_type if theta >= self.activate_occupancy else self.naive_type
        return {"fate": fate}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_receptor_subcell.py -q`
Expected: PASS (4 tests).

- [ ] **Step 5: Register the class so `local:` addresses resolve** — append to `cpm/core.py`'s registration (follow the pattern used for `BooleanSubcell`/`SBMLSubcell`; grep `register` in `cpm/core.py` and mirror it). Then verify:

Run: `pytest tests/test_core_registration.py -q`
Expected: PASS (existing registration test still green; ReceptorSubcell importable via the registry).

- [ ] **Step 6: Commit**

```bash
git add cpm/subcellular/receptor.py tests/test_receptor_subcell.py cpm/core.py
git commit -m "feat(cpm): ReceptorSubcell — Hill occupancy -> responder fate"
```

---

### Task 2: H4 down-scale coupling helper

**Files:**
- Create: `cpm/coupling.py`
- Test: `tests/test_coupling_helper.py`

**Interfaces:**
- Produces: `receptor_coupling(cell_ids: list[int], *, receptor_config: dict, receptor_addr: str = "local:!cpm.subcellular.receptor.ReceptorSubcell") -> dict`. Returns a composite-document fragment: one `ReceptorSubcell` node per cell id (keyed `receptor_<cid>`), each with `inputs={"ligand": ["field_at_cell", str(cid)]}` and `outputs={"fate": ["fates", str(cid)]}`, plus a `fates` store pre-initialized to `{str(cid): receptor_config["naive_type"]}` for every cell id (the H4 pre-init fallback). Consumed by Task 3.

- [ ] **Step 1: Write the failing test** — `tests/test_coupling_helper.py`

```python
from cpm.coupling import receptor_coupling

RCFG = {"kd": 10.0, "hill": 1.0, "conc_scale": 1.0, "activate_occupancy": 0.5,
        "naive_type": 2, "activated_type": 3}


def test_builds_one_subcell_per_cell_with_wiring():
    frag = receptor_coupling([4, 7], receptor_config=RCFG)
    assert frag["receptor_4"]["inputs"]["ligand"] == ["field_at_cell", "4"]
    assert frag["receptor_4"]["outputs"]["fate"] == ["fates", "4"]
    assert frag["receptor_7"]["inputs"]["ligand"] == ["field_at_cell", "7"]
    assert frag["receptor_4"]["address"].endswith("ReceptorSubcell")
    assert frag["receptor_4"]["config"]["kd"] == 10.0


def test_fates_store_preinitialized_for_every_cell():
    frag = receptor_coupling([4, 7], receptor_config=RCFG)
    assert frag["fates"] == {"4": 2, "7": 2}   # H4 pre-init fallback (naive_type)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_coupling_helper.py -q`
Expected: FAIL — `ModuleNotFoundError: No module named 'cpm.coupling'`.

- [ ] **Step 3: Write minimal implementation** — `cpm/coupling.py`

```python
"""Down-scale coupling helper: wire per-cell subcellular processes to the CPM
`fates` port. Encapsulates the map-key contract (a `fates` write to an absent
key is dropped on apply, so the store is pre-initialized per cell). When the
bigraph-schema keyed-map affordance (H5) lands, the pre-init here is removed and
the `fates` port declared with the affordance instead.
"""
from __future__ import annotations

RECEPTOR_ADDR = "local:!cpm.subcellular.receptor.ReceptorSubcell"


def receptor_coupling(cell_ids, *, receptor_config, receptor_addr=RECEPTOR_ADDR):
    frag = {}
    naive = int(receptor_config["naive_type"])
    fates = {}
    for cid in cell_ids:
        key = str(int(cid))
        frag[f"receptor_{cid}"] = {
            "_type": "process",
            "address": receptor_addr,
            "config": dict(receptor_config),
            "inputs": {"ligand": ["field_at_cell", key]},
            "outputs": {"fate": ["fates", key]},
        }
        fates[key] = naive  # H4 pre-init fallback
    frag["fates"] = fates
    return frag
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_coupling_helper.py -q`
Expected: PASS (2 tests).

- [ ] **Step 5: Commit**

```bash
git add cpm/coupling.py tests/test_coupling_helper.py
git commit -m "feat(cpm): receptor_coupling helper (fates wiring + pre-init)"
```

---

### Task 3: `recruitment_receptor` composite

**Files:**
- Create: `pbg_cpm_studies/composites/chemotaxis_receptor.py`
- Test: `tests/test_recruitment_receptor.py`

**Interfaces:**
- Consumes: `receptor_coupling` (Task 2); `pbg_cpm_studies.composites.chemotaxis.build_spec` and its module constants (`SOURCE_TYPE=1`, `RESPONDER_TYPE=2`, `CUE_RATE`, `CHEMO_LAMBDA`, `_cell`, the responder seed layout).
- Produces: `@composite_generator name="recruitment_receptor"` factory `recruitment_receptor(core=None, cue_rate=…, chemo_lambda=…, kd=…, seed=…, blocked=False)` returning a composite document. Helper `build_receptor_spec(*, cue_rate, chemo_lambda, blocked, seed)` returning the CPM `spec` with the naive(2)/activated(3) responder sub-types and matched contact-J. Module constants `ACTIVATED_TYPE=3`, `NAIVE_TYPE=2`, `KD_DEFAULT`.

- [ ] **Step 1: Write the failing tests** — `tests/test_recruitment_receptor.py`

```python
import process_bigraph as pb
from pbg_cpm_studies.composites import chemotaxis_receptor as CR


def test_contact_j_invariant_activated_matches_naive():
    """Activated responder (3) shares every contact-J with naive (2): only
    chemotaxis differs, so adhesion is not confounded."""
    spec = CR.build_receptor_spec(cue_rate=10.0, chemo_lambda=14.0, blocked=False, seed=17)
    j = {(min(r["a"], r["b"]), max(r["a"], r["b"])): r["j"] for r in spec["contact"]}
    for other in (0, 1, 2, 3):
        lo2, hi2 = min(2, other), max(2, other)
        lo3, hi3 = min(3, other), max(3, other)
        assert j[(lo2, hi2)] == j[(lo3, hi3)], f"J mismatch for type {other}"


def test_only_activated_type_chemotaxes():
    spec = CR.build_receptor_spec(cue_rate=10.0, chemo_lambda=14.0, blocked=False, seed=17)
    chemo = {c["type"]: c["lambda"] for c in spec["fields"][0]["chemotaxis"]}
    assert chemo.get(CR.NAIVE_TYPE, 0.0) == 0.0
    assert chemo.get(CR.ACTIVATED_TYPE) == 14.0


def test_blocked_zeroes_activated_lambda():
    spec = CR.build_receptor_spec(cue_rate=10.0, chemo_lambda=14.0, blocked=True, seed=17)
    chemo = {c["type"]: c["lambda"] for c in spec["fields"][0]["chemotaxis"]}
    assert chemo.get(CR.ACTIVATED_TYPE, 0.0) == 0.0


def test_document_preinitializes_fates_for_all_responders():
    doc = CR.recruitment_receptor(cue_rate=10.0, chemo_lambda=14.0)
    n_responders = sum(1 for c in doc["cpm"]["config"]["spec"]["cells"]
                       if c["type"] == CR.NAIVE_TYPE)
    assert len(doc["fates"]) == n_responders
    assert all(v == CR.NAIVE_TYPE for v in doc["fates"].values())


def test_smoke_run_baseline_recruits_more_than_blocked():
    from pbg_cpm_studies.chemotaxis import metrics as M
    core = pb.allocate_core()

    def final_index(blocked):
        doc = CR.recruitment_receptor(core=core, blocked=blocked, seed=17)
        comp = pb.Composite({"state": doc}, core=core)
        comp.run(40)
        world = comp.state["cpm"]["instance"].world
        return M.recruitment_index(world)

    assert final_index(blocked=False) > final_index(blocked=True)
```

*Note on the smoke test's world access:* mirror how `pbg_cpm_studies/chemotaxis/run.py` reaches the live `CPMProcess.world` after a run; if the access path differs, copy run.py's exact idiom rather than the sketch above. Keep the assertion (baseline recruits strictly more than blocked over a short run).

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_recruitment_receptor.py -q`
Expected: FAIL — module `chemotaxis_receptor` missing.

- [ ] **Step 3: Write the composite** — `pbg_cpm_studies/composites/chemotaxis_receptor.py`

Build on `chemotaxis.build_spec`: start from its cells + contact + fields, then (a) relabel each responder that should be able to activate — responders remain seeded as `NAIVE_TYPE=2`, and add contact rows for `ACTIVATED_TYPE=3` **copied from the type-2 rows** (pair-by-pair) so J is identical; (b) set the field `chemotaxis` list to `[{type: 3, lambda: chemo_lambda}]` (naive type 2 omitted ⇒ λ=0), or `[]`/λ=0 when `blocked`; (c) build the per-responder `ReceptorSubcell` wiring + pre-init via `receptor_coupling(responder_ids, receptor_config=…)`, merging its fragment into the document alongside the `cpm` node. Responder ids are the CPM cell ids of the type-2 cells (ids assigned in cell order; source is id 1, responders 2..N — confirm against `chemotaxis.build_spec`'s cell order and `CPMProcess` id assignment). Decorate with `@composite_generator(name="recruitment_receptor", parameters={…cue_rate, chemo_lambda, kd, blocked…}, visualizations=…)` mirroring `chemotaxis.recruitment`.

*(Full code assembled from the two references — `chemotaxis.build_spec`/`composite_document` and `receptor_coupling`. Keep `build_receptor_spec` and `recruitment_receptor` as the two public entry points named in the Interfaces block.)*

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_recruitment_receptor.py -q`
Expected: PASS (5 tests). If the smoke test is flaky at 40 steps, raise steps until baseline clearly separates from blocked, then keep that value.

- [ ] **Step 5: Commit**

```bash
git add pbg_cpm_studies/composites/chemotaxis_receptor.py tests/test_recruitment_receptor.py
git commit -m "feat(studies): recruitment_receptor composite (emergent response via receptor occupancy)"
```

---

### Task 4: Cited chemokine Kd (evidence translator) + run harness

**Files:**
- Modify: `workspace/references/papers.bib`
- Create: `pbg_cpm_studies/chemotaxis/run_receptor.py`
- Test: `tests/test_run_receptor.py`

**Interfaces:**
- Produces: `run_receptor(*, blocked: bool, seeds: list[int], steps: int, kd: float) -> dict` writing a summary JSON to `workspace/chemotaxis_data/results/receptor_<condition>.json` with keys `mcs: [..]`, `recruitment_index_mean: [..]`, `recruitment_index_ci: [[lo,hi],..]`, `final_mean`, `final_ci`, `activation_by_cell` (per-responder activation fraction over time for the spatial figure), `kd`, `seeds`. Mirrors `pbg_cpm_studies/chemotaxis/run.py`'s series construction.

- [ ] **Step 1: Add the cited chemokine–receptor entry** to `workspace/references/papers.bib`. Use a well-characterized chemokine–receptor affinity — default **fMLP–FPR1** (formyl-peptide receptor), a nanomolar Kd. Add a real BibTeX entry with `doi`, and record the numeric Kd (in the paper's units) in the entry's `note`. Confirm the Kd value against the cited paper before using it downstream (open question #2 — pick and verify one concrete affinity; do not invent a number).

- [ ] **Step 2: Write the failing test** — `tests/test_run_receptor.py`

```python
from pbg_cpm_studies.chemotaxis.run_receptor import run_receptor


def test_summary_shape_and_recruitment_ordering(tmp_path):
    base = run_receptor(blocked=False, seeds=[17, 29, 43], steps=60, kd=10.0)
    blk = run_receptor(blocked=True, seeds=[17, 29, 43], steps=60, kd=10.0)
    for s in (base, blk):
        assert len(s["mcs"]) == len(s["recruitment_index_mean"]) == len(s["recruitment_index_ci"])
        assert 0.0 <= s["final_mean"] <= 1.0
    assert base["final_mean"] > blk["final_mean"]
```

- [ ] **Step 3: Run test to verify it fails**

Run: `pytest tests/test_run_receptor.py -q`
Expected: FAIL — module missing.

- [ ] **Step 4: Implement `run_receptor.py`** by adapting `pbg_cpm_studies/chemotaxis/run.py`: for each seed, build `recruitment_receptor(blocked=…, seed=…, kd=…)`, run to `steps`, sample `metrics.recruitment_index_from_coms` at intervals, aggregate mean + 95% CI across seeds (normal-approx or seed spread), and record per-responder activation (fate==ACTIVATED_TYPE) per frame. Write the JSON summary.

- [ ] **Step 5: Run test to verify it passes** — `pytest tests/test_run_receptor.py -q`. Expected: PASS.

- [ ] **Step 6: Generate the real summaries** for ≥5 seeds and commit the data:

```bash
python -c "from pbg_cpm_studies.chemotaxis.run_receptor import run_receptor; \
run_receptor(blocked=False, seeds=[17,29,43,61,89], steps=500, kd=<KD>); \
run_receptor(blocked=True,  seeds=[17,29,43,61,89], steps=500, kd=<KD>)"
git add workspace/references/papers.bib pbg_cpm_studies/chemotaxis/run_receptor.py \
        tests/test_run_receptor.py workspace/chemotaxis_data/results/receptor_*.json
git commit -m "feat(studies): cited chemokine Kd + multi-seed receptor run harness"
```
(Use the confirmed Kd value from Step 1 in place of `<KD>`.)

---

### Task 5: Receptor studies (baseline + blocked) with calibration anchor

**Files:**
- Create: `workspace/studies/recruitment-receptor-baseline/study.yaml`
- Create: `workspace/studies/recruitment-receptor-blocked/study.yaml`
- Modify: `workspace/investigations/chemotactic-recruitment/investigation.yaml` (add the two studies + acceptance criteria)
- Test: `tests/test_receptor_studies.py`

**Interfaces:**
- Consumes: the composite `pbg_cpm_studies.composites.chemotaxis_receptor.recruitment_receptor`; the summaries from Task 4; the existing `recruitment_index` measure kind.
- Produces: two schema-v3 `study.yaml` files mirroring `workspace/studies/recruitment-baseline/study.yaml`, each with a `behavior_tests[].pass_if` and a `calibration_anchor{literature_target: <Kd>, cites: [<bibkey>], resolution: model}` on the receptor `kd` parameter; `model_settings[]` recording `kd/hill/conc_scale/activate_occupancy` with `cites`/provenance.

- [ ] **Step 1: Write the failing test** — `tests/test_receptor_studies.py`

```python
from pathlib import Path
import yaml

BASE = Path("workspace/studies/recruitment-receptor-baseline/study.yaml")


def test_baseline_study_wired_and_anchored():
    doc = yaml.safe_load(BASE.read_text())
    assert doc["baseline"][0]["composite"].endswith("chemotaxis_receptor.recruitment_receptor")
    test = doc["behavior_tests"][0]
    anchor = test.get("calibration_anchor") or {}
    assert anchor.get("cites"), "acceptance band must cite a source"
    assert anchor.get("literature_target") is not None
```

- [ ] **Step 2: Run test to verify it fails** — `pytest tests/test_receptor_studies.py -q`. Expected: FAIL (file missing).

- [ ] **Step 3: Author the two study.yaml files** mirroring `recruitment-baseline/study.yaml` (schema_version 3, `investigation: chemotactic-recruitment`). Baseline: `baseline[0].composite = pbg_cpm_studies.composites.chemotaxis_receptor.recruitment_receptor`, `params: {cue_rate: 10.0, chemo_lambda: 14.0, kd: <KD>}`; behavior_test `measure: {kind: recruitment_index, condition: recruitment-receptor-baseline, stat: final}`, `pass_if: {op: gt, value: 0.5}`, plus `calibration_anchor: {literature_target: <Kd>, cites: [<bibkey>], resolution: model}`. Blocked: `params: {cue_rate: 10.0, chemo_lambda: 14.0, kd: <KD>, blocked: true}`, `pass_if: {op: lt, value: 0.1}`. Fill observed values + CIs from the Task 4 summaries. Add both to `investigation.yaml` `studies:` and `acceptance_criteria:`, and add the receptor realization's `at_a_glance` rows.

- [ ] **Step 4: Run test to verify it passes** — `pytest tests/test_receptor_studies.py -q`. Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add workspace/studies/recruitment-receptor-baseline workspace/studies/recruitment-receptor-blocked \
        workspace/investigations/chemotactic-recruitment/investigation.yaml tests/test_receptor_studies.py
git commit -m "feat(studies): receptor baseline + blocked studies with cited Kd anchor"
```

---

### Task 6: Figures — spatial activation map + condition overlay

**Files:**
- Create: `pbg_cpm_studies/visualizations/receptor_studies.py`
- Test: `tests/test_receptor_viz.py`

**Interfaces:**
- Consumes: the Task 4 summaries (baked data, per the `chemotaxis_studies.py` convention).
- Produces: two `@as_visualization` functions — `ReceptorRecruitment` (recruitment index vs MCS with CI ribbons, baseline vs receptor-blocked vs no-cue) and `ReceptorActivationMap` (per-responder activation over space/time: the activation front). Registered `local:ReceptorRecruitment`, `local:ReceptorActivationMap`.

- [ ] **Step 1: Write the failing test** — `tests/test_receptor_viz.py`

```python
from pbg_cpm_studies.visualizations import receptor_studies as V


def test_visualizations_render_html():
    for fn in (V.ReceptorRecruitment, V.ReceptorActivationMap):
        html = fn()
        assert isinstance(html, str) and "plotly" in html.lower()
```

- [ ] **Step 2: Run test to verify it fails** — `pytest tests/test_receptor_viz.py -q`. Expected: FAIL (module missing).

- [ ] **Step 3: Implement the two viz functions** in `pbg_cpm_studies/visualizations/receptor_studies.py`, following `pbg_cpm_studies/visualizations/chemotaxis_studies.py` exactly (imports `from viva_superpowers.visualization import as_visualization`, baked series, Plotly trace dicts, dark theme). `ReceptorRecruitment`: three conditions with CI ribbons (`fill="tonexty"`) from `receptor_*.json`. `ReceptorActivationMap`: a heatmap / scatter of per-responder activation state vs distance-to-source over MCS (the activation front), from `activation_by_cell`.

- [ ] **Step 4: Run test to verify it passes** — `pytest tests/test_receptor_viz.py -q`. Expected: PASS.

- [ ] **Step 5: Wire the figures into the studies** — add `visualizations: [{name, address: local:ReceptorRecruitment}, {name, address: local:ReceptorActivationMap}]` to both receptor study.yaml files.

- [ ] **Step 6: Commit**

```bash
git add pbg_cpm_studies/visualizations/receptor_studies.py tests/test_receptor_viz.py \
        workspace/studies/recruitment-receptor-baseline/study.yaml \
        workspace/studies/recruitment-receptor-blocked/study.yaml
git commit -m "feat(viz): receptor recruitment overlay + spatial activation map"
```

---

### Task 7: Full-suite green + investigation sanity

- [ ] **Step 1:** Run the whole suite: `pytest -q`. Expected: all green (new + existing).
- [ ] **Step 2:** Confirm `ReceptorSubcell`, `receptor_coupling`, and `recruitment_receptor` are import-clean and the composite resolves through the registry (re-run `pytest tests/test_core_registration.py tests/test_recruitment_receptor.py -q`).
- [ ] **Step 3:** Sanity-read the two new study reports render (no missing-figure/observable errors) if a local dashboard is available; otherwise verify the study YAML loads and the `recruitment_index` measure resolves against the run summaries.
- [ ] **Step 4: Commit** any fixups, then stop — Increment 1 complete. Next: H5 (bigraph-schema keyed-map affordance) and Increment 2 (dose-response + evidence-infra H2) get their own plans.

---

## Self-review notes

- **Spec coverage (Increment 1 rows):** ReceptorSubcell (§2.1 → T1); H4 coupling helper (§4 H4 → T2); composite + contact-J invariant + fates pre-init (§2.2–2.3 → T3); evidence translator / cited Kd (§3.1 → T4/T5); baseline+blocked studies + multi-seed + CIs (§5 → T4/T5); figures (§8 Fig 2 + overlay → T6). Deferred rows (H1/H2/H3/H5, dose-response Fig 1) are explicitly out of Increment 1.
- **Placeholders:** `<KD>`/`<bibkey>` are intentional implementation-time literature values (open question #2), gated by an explicit "confirm against the paper, do not invent" step (T4-S1). The Task 3 composite body references two exact source modules rather than inlining assembled code; every public name it must expose is pinned in its Interfaces block.
- **Type consistency:** `NAIVE_TYPE=2`/`ACTIVATED_TYPE=3`, `recruitment_index[_from_coms]`, `receptor_coupling(...)->frag`, and the `fates`/`field_at_cell` string-keyed maps are used identically across tasks.
