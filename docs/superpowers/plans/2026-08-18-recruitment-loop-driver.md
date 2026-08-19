# Recruitment Loop Driver (Phase 2) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the adaptive-recruitment mechanism ladder (`static_lambda` → `hill_occupancy` → `adaptive_receptor`) as real `process_bigraph` code over the existing CPM engine, plus a deterministic model-building loop driver that climbs it from a locked contract and captures a real, emergent `trajectory.json`.

**Architecture:** One new per-cell Process `AdaptiveReceptorSubcell` subsumes the ladder — it extends the existing `ReceptorSubcell` with a persistent per-cell adaptation store `m` (dm/dt = ε(θ−m)) and a `background` offset added to the sensed ligand. `epsilon=0` reduces it to pure Hill occupancy (the middle rung AND the discriminating knockout); `epsilon>0` gives fold-change detection. Background-cue conditions are a subcell config (`background` added to the ligand), NOT an engine change — no Rust. The loop driver mirrors `viva-casebook/scripts/build_loop_demo.py`, using `viva_superpowers.{loop_state, test_audit, test_contract}`.

**Tech Stack:** Python, `process_bigraph` (Composite, Process, allocate_core), the compiled `cpm` engine (`cpm.subcellular.receptor.ReceptorSubcell`, `cpm.coupling.receptor_coupling`, `cpm.load_world`), `viva_superpowers` loop framework, pytest.

**Spec:** `docs/superpowers/specs/2026-08-18-recruitment-model-building-design.md` (Part B, §4 — B1/B2/B3).

## Global Constraints

- Worktree root: the repo checkout root; run Python via `.venv/bin/python`.
- **No new Rust.** The adaptive mechanism is pure Python over the existing `cpm_core` bindings. Background cue = a `background` float added to the ligand inside the subcell, NOT a new field API.
- Reuse, don't re-duplicate: import `cpm.subcellular.receptor.ReceptorSubcell` and call its `.occupancy()` for the Hill formula rather than copy-pasting `c**h/(kd**h+c**h)` a fourth time.
- Per-cell **map-write pre-init rule** (from `cpm/coupling.py`): a map-write to an absent key is dropped on apply. Any new per-cell map port (e.g. `adaptation`) MUST be pre-initialized in the coupling fragment, exactly as `receptor_coupling` pre-inits `fates`.
- Composites run through `pb.Composite({"state": doc}, core=core)` + `comp.run(n)`; the live CPM world is at `comp.state["cpm"]["instance"].world`. Per-cell subcells only fire on the Composite path (NOT raw `load_world`).
- Keep loop-driver runs SMALL (reduced lattice / fewer seeds / fewer steps) so `trajectory.json` regenerates in minutes and CI stays fast. The committed trajectory is generated offline, not in CI.
- Emitted metric keys (Phase 1 confirmed): `recruitment_index`, `mean_distance`, `activation_by_cell`. The distance readout's emitted key is `mean_distance`.

## Known risk (read before Task 4)

Whether the CPM adaptive mechanism actually exhibits the clean ladder (`adaptive` recruits at high background where `hill` saturates and fails) is an EMPIRICAL question the parameters must satisfy. Task 4 is the validation gate: it tunes `epsilon`, the `background` levels, `activate_occupancy`, and the recruitment `radius` until the ladder separates, and records the working values. If a clean separation is not achievable after tuning, STOP and report — the contract bands (Task 6) and the spec's ladder claim depend on it. Do not fake the separation.

---

## File structure

New:
- `cpm/subcellular/adaptive_receptor.py` — `AdaptiveReceptorSubcell` (Process)
- `cpm/coupling.py` — ADD `adaptation_coupling(...)` alongside `receptor_coupling` (same file)
- `pbg_cpm_studies/composites/chemotaxis_adaptive.py` — `build_adaptive_spec`, `composite_document`, `@composite_generator recruitment_adaptive`, the 5 background conditions
- `pbg_cpm_studies/model_building/__init__.py`
- `pbg_cpm_studies/model_building/mechanisms.py` — the LIBRARY (static/hill/adaptive → composite+config) + `simulate_condition`
- `pbg_cpm_studies/model_building/navigate.py` — deterministic NAVIGATE policy
- `workspace/studies/recruitment-adaptive/study.yaml` — the locked contract
- `scripts/build_recruitment_loop.py` — the loop driver → `trajectory.json`
- Tests: `tests/test_adaptive_receptor.py`, `tests/test_adaptation_coupling.py`, `tests/test_chemotaxis_adaptive.py`, `tests/test_recruitment_ladder.py` (the validation), `tests/test_recruitment_loop.py`

Modified:
- `workspace/references/papers.bib` — add `barkai1997`, `tu2008` (adaptation refs)

---

### Task 1: `AdaptiveReceptorSubcell` Process

**Files:**
- Create: `cpm/subcellular/adaptive_receptor.py`
- Test: `tests/test_adaptive_receptor.py`

**Interfaces:**
- Produces: `AdaptiveReceptorSubcell(Process)` with `config_schema` keys `kd, hill, conc_scale, activate_occupancy, epsilon, background, naive_type, activated_type`; `inputs() -> {"ligand": "float", "m_prev": "float"}`; `outputs() -> {"fate": "overwrite[integer]", "m": "overwrite[float]"}`; a `.signal(ligand, m_prev) -> (theta, m_new, response)` helper. Reuses `ReceptorSubcell.occupancy` for θ.

- [ ] **Step 1: Write the failing unit tests**

```python
# tests/test_adaptive_receptor.py
import process_bigraph as pb
from cpm.subcellular.adaptive_receptor import AdaptiveReceptorSubcell

def _proc(**over):
    cfg = {"kd": 2.9, "hill": 2.0, "conc_scale": 0.02, "activate_occupancy": 0.5,
           "epsilon": 0.1, "background": 0.0, "naive_type": 2, "activated_type": 3}
    cfg.update(over)
    return AdaptiveReceptorSubcell(cfg, core=pb.allocate_core())

def test_epsilon_zero_is_pure_hill():
    # with epsilon=0, m stays 0 so response == theta (the hill_occupancy rung / knockout)
    p = _proc(epsilon=0.0)
    theta, m_new, response = p.signal(ligand=200.0, m_prev=0.0)
    assert m_new == 0.0
    assert abs(response - theta) < 1e-9

def test_background_adds_to_ligand():
    # a uniform background raises sensed concentration -> higher occupancy
    p = _proc(epsilon=0.0, background=100.0)
    theta_bg, _, _ = p.signal(ligand=100.0, m_prev=0.0)
    theta_nobg, _, _ = _proc(epsilon=0.0, background=0.0).signal(ligand=100.0, m_prev=0.0)
    assert theta_bg > theta_nobg

def test_m_relaxes_toward_theta():
    # repeated updates at constant ligand drive m toward theta (adaptation)
    p = _proc(epsilon=0.3)
    theta, _, _ = p.signal(ligand=300.0, m_prev=0.0)
    m = 0.0
    for _ in range(50):
        _, m, _ = p.signal(ligand=300.0, m_prev=m)
    assert abs(m - theta) < 0.05          # adapted setpoint tracks occupancy
    _, _, response = p.signal(ligand=300.0, m_prev=m)
    assert response < 0.05                 # steady-state response adapts away

def test_update_emits_fate_and_m():
    p = _proc(epsilon=0.1)
    out = p.update({"ligand": 300.0, "m_prev": 0.0}, 1.0)
    assert set(out) == {"fate", "m"}
    assert out["fate"] in (2, 3)
```

- [ ] **Step 2: Run to verify they fail**

Run: `.venv/bin/python -m pytest tests/test_adaptive_receptor.py -q`
Expected: FAIL (module `cpm.subcellular.adaptive_receptor` does not exist).

- [ ] **Step 3: Implement `AdaptiveReceptorSubcell`**

```python
# cpm/subcellular/adaptive_receptor.py
from process_bigraph import Process
from cpm.subcellular.receptor import ReceptorSubcell

class AdaptiveReceptorSubcell(Process):
    """Receptor occupancy with slow adaptation (fold-change detection).

    theta = Hill occupancy of (ligand + background); m is a slow per-cell
    setpoint, dm/dt = epsilon*(theta - m); the response that drives chemotaxis
    is r = max(0, theta - m). epsilon=0 reduces this to pure Hill occupancy
    (response == theta) — the hill_occupancy rung and the adaptation knockout.
    A cell activates (naive->activated) when r >= activate_occupancy.
    """
    config_schema = {
        "kd": {"_type": "float", "_default": 2.9},
        "hill": {"_type": "float", "_default": 2.0},
        "conc_scale": {"_type": "float", "_default": 0.02},
        "activate_occupancy": {"_type": "float", "_default": 0.5},
        "epsilon": {"_type": "float", "_default": 0.1},
        "background": {"_type": "float", "_default": 0.0},
        "naive_type": {"_type": "integer", "_default": 2},
        "activated_type": {"_type": "integer", "_default": 3},
    }

    def initialize(self, config):
        # reuse the canonical Hill formula (no fourth copy of c^h/(kd^h+c^h))
        self._hill = ReceptorSubcell(
            {"kd": config["kd"], "hill": config["hill"], "conc_scale": config["conc_scale"],
             "activate_occupancy": config["activate_occupancy"],
             "naive_type": config["naive_type"], "activated_type": config["activated_type"]},
            core=self.core)
        self.epsilon = float(config["epsilon"]); self.background = float(config["background"])
        self.activate_occupancy = float(config["activate_occupancy"])
        self.naive_type = int(config["naive_type"]); self.activated_type = int(config["activated_type"])

    def inputs(self):  return {"ligand": "float", "m_prev": "float"}
    def outputs(self): return {"fate": "overwrite[integer]", "m": "overwrite[float]"}

    def signal(self, ligand, m_prev):
        theta = self._hill.occupancy(float(ligand) + self.background)
        m_new = float(m_prev) + self.epsilon * (theta - float(m_prev))
        response = max(0.0, theta - float(m_prev))
        return theta, m_new, response

    def update(self, state, interval):
        s = state or {}
        _, m_new, response = self.signal(s.get("ligand", 0.0), s.get("m_prev", 0.0))
        fate = self.activated_type if response >= self.activate_occupancy else self.naive_type
        return {"fate": fate, "m": m_new}
```

- [ ] **Step 4: Run to verify pass**

Run: `.venv/bin/python -m pytest tests/test_adaptive_receptor.py -q`
Expected: PASS (4 tests). If `test_m_relaxes_toward_theta` is flaky on the tolerance, adjust epsilon/iteration count in the TEST only (not the class).

- [ ] **Step 5: Commit**

```bash
git add cpm/subcellular/adaptive_receptor.py tests/test_adaptive_receptor.py
git commit -m "feat: AdaptiveReceptorSubcell — occupancy + slow adaptation (fold-change detection)"
```

---

### Task 2: `adaptation_coupling` wiring helper

**Files:**
- Modify: `cpm/coupling.py` (add function; do not change `receptor_coupling`)
- Test: `tests/test_adaptation_coupling.py`

**Interfaces:**
- Consumes: `AdaptiveReceptorSubcell` (Task 1).
- Produces: `adaptation_coupling(cell_ids, *, receptor_config, addr=ADAPTIVE_ADDR) -> dict` — per cell wires `ligand <- ["field_at_cell", key]`, `m_prev <- ["adaptation", key]`, `fate -> ["fates", key]`, `m -> ["adaptation", key]`; pre-inits BOTH `fates` (to naive_type) and `adaptation` (to 0.0) maps. `ADAPTIVE_ADDR = "local:!cpm.subcellular.adaptive_receptor.AdaptiveReceptorSubcell"`.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_adaptation_coupling.py
from cpm.coupling import adaptation_coupling, ADAPTIVE_ADDR

def test_fragment_shape_and_preinit():
    cfg = {"kd": 2.9, "hill": 2.0, "conc_scale": 0.02, "activate_occupancy": 0.5,
           "epsilon": 0.1, "background": 0.0, "naive_type": 2, "activated_type": 3}
    frag = adaptation_coupling([2, 3], receptor_config=cfg)
    node = frag["receptor_2"]
    assert node["address"] == ADAPTIVE_ADDR
    assert node["inputs"]["ligand"] == ["field_at_cell", "2"]
    assert node["inputs"]["m_prev"] == ["adaptation", "2"]
    assert node["outputs"]["fate"] == ["fates", "2"]
    assert node["outputs"]["m"] == ["adaptation", "2"]
    # BOTH per-cell maps pre-initialized (map-write-to-absent-key is dropped)
    assert frag["fates"] == {"2": 2, "3": 2}
    assert frag["adaptation"] == {"2": 0.0, "3": 0.0}
```

- [ ] **Step 2: Run to verify fail**

Run: `.venv/bin/python -m pytest tests/test_adaptation_coupling.py -q`
Expected: FAIL (ImportError: adaptation_coupling).

- [ ] **Step 3: Implement (append to `cpm/coupling.py`)**

```python
ADAPTIVE_ADDR = "local:!cpm.subcellular.adaptive_receptor.AdaptiveReceptorSubcell"

def adaptation_coupling(cell_ids, *, receptor_config, addr=ADAPTIVE_ADDR):
    """Wire one AdaptiveReceptorSubcell per responder cell. Mirrors
    receptor_coupling but adds a persistent per-cell adaptation store `m`
    (map port `adaptation`), pre-initialized to 0.0 (map-write-to-absent-key
    is dropped on apply, so both `fates` and `adaptation` must be pre-inited)."""
    frag = {}
    naive = int(receptor_config["naive_type"])
    fates, adaptation = {}, {}
    for cid in cell_ids:
        key = str(int(cid))
        frag[f"receptor_{cid}"] = {
            "_type": "process", "address": addr, "config": dict(receptor_config),
            "inputs": {"ligand": ["field_at_cell", key], "m_prev": ["adaptation", key]},
            "outputs": {"fate": ["fates", key], "m": ["adaptation", key]},
        }
        fates[key] = naive
        adaptation[key] = 0.0
    frag["fates"] = fates
    frag["adaptation"] = adaptation
    return frag
```

- [ ] **Step 4: Run to verify pass**

Run: `.venv/bin/python -m pytest tests/test_adaptation_coupling.py -q`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add cpm/coupling.py tests/test_adaptation_coupling.py
git commit -m "feat: adaptation_coupling — per-cell adaptive receptor wiring with pre-inited m store"
```

---

### Task 3: `chemotaxis_adaptive` composite + background conditions

**Files:**
- Create: `pbg_cpm_studies/composites/chemotaxis_adaptive.py`
- Test: `tests/test_chemotaxis_adaptive.py`

**Interfaces:**
- Consumes: `chemotaxis.build_spec` (via `import ... as CT`), `chemotaxis_receptor` patterns, `adaptation_coupling` (Task 2).
- Produces: `build_adaptive_spec(*, cue_rate, chemo_lambda, seed)` (reuses `CT.build_spec`, appends activated-type contact rows, sets chemotaxis on the activated type); `composite_document(*, cue_rate, chemo_lambda, kd, epsilon, background, seed, blocked)`; `@composite_generator(name="recruitment_adaptive", ...)`; `CONDITIONS: dict[str, dict]` mapping `low_bg/mid_bg/high_bg/high_bg_blocked/high_bg_knockout` → the config overrides (background level; `blocked`→chemotaxis off; `high_bg_knockout`→`epsilon=0`).

- [ ] **Step 1: Write the failing composite smoke test**

```python
# tests/test_chemotaxis_adaptive.py
import process_bigraph as pb
from pbg_cpm_studies.composites import chemotaxis_adaptive as CA
from pbg_cpm_studies.chemotaxis import metrics as M

def _final_index(condition, seed=17, steps=40):
    core = pb.allocate_core()
    cfg = CA.CONDITIONS[condition]
    doc = CA.recruitment_adaptive(core=core, seed=seed, **cfg)
    comp = pb.Composite({"state": doc}, core=core)
    comp.run(steps)
    world = comp.state["cpm"]["instance"].world
    return M.recruitment_index(world, responder_type=CA.ACTIVATED_TYPE)

def test_conditions_exist():
    assert set(CA.CONDITIONS) == {"low_bg", "mid_bg", "high_bg", "high_bg_blocked", "high_bg_knockout"}

def test_low_bg_recruits_and_blocked_does_not():
    assert _final_index("low_bg") > _final_index("high_bg_blocked")
```

- [ ] **Step 2: Run to verify fail**

Run: `.venv/bin/python -m pytest tests/test_chemotaxis_adaptive.py -q`
Expected: FAIL (module missing).

- [ ] **Step 3: Implement `chemotaxis_adaptive.py`**

Mirror `chemotaxis_receptor.py` closely. Key points: reuse `CT.build_spec`; `ACTIVATED_TYPE=3`, `NAIVE_TYPE=2`; copy the activated-type contact rows from naive; set `spec["fields"][0]["chemotaxis"] = [] if blocked else [{"type": ACTIVATED_TYPE, "lambda": chemo_lambda}]`; build `receptor_config = {kd, hill, conc_scale, activate_occupancy, epsilon, background, naive_type, activated_type}`; assemble the `cpm` node (identical to `chemotaxis_receptor.composite_document`); then `doc.update(adaptation_coupling(responder_ids, receptor_config=receptor_config))`. Then:

```python
CONDITIONS = {
    "low_bg":           {"background": 0.0,   "epsilon": 0.1, "blocked": False},
    "mid_bg":           {"background": 80.0,  "epsilon": 0.1, "blocked": False},
    "high_bg":          {"background": 200.0, "epsilon": 0.1, "blocked": False},
    "high_bg_blocked":  {"background": 200.0, "epsilon": 0.1, "blocked": True},   # receptor_gating
    "high_bg_knockout": {"background": 200.0, "epsilon": 0.0, "blocked": False},  # adaptation off
}

@composite_generator(name="recruitment_adaptive", default_n_steps=40, visualizations=_VIZ,
    description="Adaptive receptor recruitment (fold-change detection).",
    parameters={
        "cue_rate": {"type": "float", "default": CUE_RATE},
        "chemo_lambda": {"type": "float", "default": CHEMO_LAMBDA},
        "kd": {"type": "float", "default": 2.9},
        "epsilon": {"type": "float", "default": 0.1},
        "background": {"type": "float", "default": 0.0},
        "blocked": {"type": "boolean", "default": False},
    })
def recruitment_adaptive(core=None, cue_rate=CUE_RATE, chemo_lambda=CHEMO_LAMBDA,
                         kd=2.9, epsilon=0.1, background=0.0, seed=SEED, blocked=False):
    return composite_document(cue_rate=cue_rate, chemo_lambda=chemo_lambda, kd=kd,
                              epsilon=epsilon, background=background, seed=seed, blocked=blocked)
```

(The `background` values above are PLACEHOLDER starting points — Task 4 tunes them.)

- [ ] **Step 4: Run to verify pass + generator registers**

Run:
```bash
.venv/bin/python -m pytest tests/test_chemotaxis_adaptive.py -q
.venv/bin/python -c "import pbg_cpm_studies.composites.chemotaxis_adaptive; from viva_superpowers.composite_generator import discover_generators; print('recruitment_adaptive' in str(sorted(discover_generators(extra_packages=['pbg_cpm_studies.composites']))))"
```
Expected: tests PASS; generator discovery prints `True`.

- [ ] **Step 5: Commit**

```bash
git add pbg_cpm_studies/composites/chemotaxis_adaptive.py tests/test_chemotaxis_adaptive.py
git commit -m "feat: chemotaxis_adaptive composite + 5 background conditions"
```

---

### Task 4: Empirical ladder validation + parameter calibration (RISK GATE)

**Files:**
- Create: `tests/test_recruitment_ladder.py`
- Modify (calibration only): `pbg_cpm_studies/composites/chemotaxis_adaptive.py` (`CONDITIONS` background levels, default `epsilon`, `activate_occupancy`)

**Interfaces:**
- Produces: a proven parameter set where the ladder separates, encoded in `CONDITIONS` + defaults. This is the empirical foundation the contract (Task 6) cites.

- [ ] **Step 1: Write the ladder discrimination test (multi-seed, tolerant)**

```python
# tests/test_recruitment_ladder.py — the emergent ladder must separate
import process_bigraph as pb
from statistics import mean
from pbg_cpm_studies.composites import chemotaxis_adaptive as CA
from pbg_cpm_studies.chemotaxis import metrics as M

SEEDS = [17, 29, 43]; STEPS = 40

def _idx(condition, **over):
    vals = []
    for s in SEEDS:
        core = pb.allocate_core()
        cfg = dict(CA.CONDITIONS[condition]); cfg.update(over)
        comp = pb.Composite({"state": CA.recruitment_adaptive(core=core, seed=s, **cfg)}, core=core)
        comp.run(STEPS)
        vals.append(M.recruitment_index(comp.state["cpm"]["instance"].world,
                                        responder_type=CA.ACTIVATED_TYPE))
    return mean(vals)

def test_adaptive_beats_hill_at_high_background():
    # hill == adaptation knockout (epsilon=0) saturates and fails; adaptive recovers
    adaptive_high = _idx("high_bg")                 # epsilon>0
    hill_high     = _idx("high_bg_knockout")        # epsilon=0  (the middle rung)
    assert adaptive_high > 0.4                       # adaptive recruits at high bg
    assert adaptive_high - hill_high > 0.2           # and clearly beats hill there

def test_receptor_gating_abolishes():
    assert _idx("high_bg_blocked") < 0.2             # blocking the response kills recruitment

def test_low_background_recruits_for_both():
    assert _idx("low_bg") > 0.4                       # sanity: both rungs work at low bg
```

- [ ] **Step 2: Run it — EXPECT it may fail first, then calibrate**

Run: `.venv/bin/python -m pytest tests/test_recruitment_ladder.py -q`
- If it PASSES with the Task-3 defaults, record the values and skip to Step 4.
- If it FAILS, calibrate (Step 3). This is expected — the starting `background`/`epsilon` are guesses.

- [ ] **Step 3: Calibrate the parameters (empirical loop)**

Sweep, by editing `CONDITIONS`/defaults and re-running the test:
- `background` for `high_bg*`: high enough that `high_bg_knockout` (epsilon=0, pure Hill) SATURATES (θ→1, gradient contrast lost → recruitment drops below ~0.3), but the source gradient still rises above it locally. Start 200; try 150/300/500.
- `epsilon`: fast enough to adapt within `STEPS` ticks, slow enough to leave a transient gradient response. Try 0.05/0.1/0.2/0.4.
- `activate_occupancy`: the response threshold r=θ−m must be crossable at high bg for adaptive but not for hill. Try 0.2/0.3/0.5.
- Optionally raise `STEPS` (e.g. 60) if adaptive needs more ticks to climb.

Record each sweep line with `log`-style prints. Land on values where all three tests pass across the 3 seeds. Write the final chosen values into `CONDITIONS` + the generator defaults, with a comment citing this calibration.

**If no parameter set separates the ladder after a reasonable sweep (~12 configs):** STOP. Write a short findings note to `docs/superpowers/notes/2026-08-18-ladder-calibration.md` describing what was tried and the closest behavior, and report BLOCKED to the coordinator — the spec's ladder claim needs revisiting (e.g. a different adaptation form, or reframing hill's failure mode). Do not weaken the test thresholds to force a pass.

- [ ] **Step 4: Commit the validated parameters**

```bash
git add tests/test_recruitment_ladder.py pbg_cpm_studies/composites/chemotaxis_adaptive.py
git commit -m "test: empirical recruitment ladder separation + calibrated adaptive params"
```

---

### Task 5: Mechanism library + NAVIGATE policy

**Files:**
- Create: `pbg_cpm_studies/model_building/__init__.py`, `pbg_cpm_studies/model_building/mechanisms.py`, `pbg_cpm_studies/model_building/navigate.py`
- Test: extend `tests/test_recruitment_ladder.py` or a new `tests/test_navigate.py`

**Interfaces:**
- Consumes: `chemotaxis` (static), `chemotaxis_adaptive` (hill via epsilon=0, adaptive via epsilon>0), calibrated `CONDITIONS` (Task 4).
- Produces:
  - `mechanisms.LIBRARY: dict[str, dict]` — `{"static_lambda": {...}, "hill_occupancy": {...}, "adaptive_receptor": {...}}`, each `{"cite": <bib_key>, "test": <which contract test it satisfies>, "build": <callable(core, condition, seed)->doc>}`. `static_lambda`→`chemotaxis.recruitment` (no receptor gating). `hill_occupancy`→`chemotaxis_adaptive` with `epsilon=0`. `adaptive_receptor`→`chemotaxis_adaptive` with the calibrated epsilon.
  - `mechanisms.simulate_condition(mechanism, condition, *, seeds, steps) -> float` — mean recruitment_index (mirrors Task-4 `_idx`, mechanism-aware).
  - `navigate.next_mechanism(active, graded_axes) -> str | None` — deterministic policy: given the currently-active mechanism and the graded contract axes, return the next mechanism to install (the one the most-negative HARD failing test pulls on), or `None` when all hard tests pass. Ordering: a failing `receptor_gating` → install `hill_occupancy`; a failing `recruits_high` → install `adaptive_receptor`.

- [ ] **Step 1: Write failing tests for the policy**

```python
# tests/test_navigate.py
from pbg_cpm_studies.model_building import mechanisms, navigate

def test_library_has_three_rungs():
    assert set(mechanisms.LIBRARY) == {"static_lambda", "hill_occupancy", "adaptive_receptor"}

def test_navigate_climbs_static_to_hill_to_adaptive():
    # receptor_gating failing (static can't gate) -> install hill
    axes = [{"id": "receptor_gating", "severity": "hard", "verdict": "mismatch", "margin": -0.5},
            {"id": "recruits_high",   "severity": "hard", "verdict": "mismatch", "margin": -0.3}]
    assert navigate.next_mechanism("static_lambda", axes) == "hill_occupancy"
    # now receptor_gating passes, recruits_high still fails -> install adaptive
    axes2 = [{"id": "receptor_gating", "severity": "hard", "verdict": "within_tol", "margin": 0.1},
             {"id": "recruits_high",   "severity": "hard", "verdict": "mismatch", "margin": -0.3}]
    assert navigate.next_mechanism("hill_occupancy", axes2) == "adaptive_receptor"
    # all hard pass -> done
    axes3 = [{"id": "recruits_high", "severity": "hard", "verdict": "within_tol", "margin": 0.1}]
    assert navigate.next_mechanism("adaptive_receptor", axes3) is None
```

- [ ] **Step 2-4: Verify fail → implement `mechanisms.py` + `navigate.py` → verify pass**

`navigate.next_mechanism`: filter axes to hard mismatches; if none → `None`; pick the most-negative-margin failing axis; map its `id` → the mechanism that fixes it via a static table `{"receptor_gating": "hill_occupancy", "recruits_high": "adaptive_receptor", "recruits_low": "static_lambda"}`; if the mapped mechanism == active (already installed but still failing), fall through to the next rung. Run `.venv/bin/python -m pytest tests/test_navigate.py -q` → PASS.

- [ ] **Step 5: Commit**

```bash
git add pbg_cpm_studies/model_building tests/test_navigate.py
git commit -m "feat: mechanism library (static/hill/adaptive) + deterministic NAVIGATE policy"
```

---

### Task 6: The `recruitment-adaptive` contract study

**Files:**
- Create: `workspace/studies/recruitment-adaptive/study.yaml`
- Modify: `workspace/references/papers.bib` (add `barkai1997`, `tu2008`)

**Interfaces:**
- Consumes: calibrated bands (Task 4), the composite id `chemotaxis_adaptive.recruitment_adaptive`.
- Produces: a schema-3 study that participates in the loop — `question`, `purpose.mechanism`, `requires`, `sourcing`, `catalog`, `behavior_tests` (the 4-test contract from the spec §4/B2, with the CALIBRATED band values from Task 4), `conditions.baseline` (composite + calibrated params), `controls` (adaptation_knockout), `readouts`.

- [ ] **Step 1: Add bib entries** for `barkai1997` (Barkai & Leibler, robustness/perfect adaptation) and `tu2008` (Tu, Shimizu, Berg — fold-change detection). Commit.
- [ ] **Step 2: Author `study.yaml`** using the spec §4/B2 block, substituting the Task-4 calibrated band low/high values into `recruits_low`/`recruits_high` (`in_range`) and the `receptor_gating`/`adaptation_knockout_collapses` thresholds. Set `conditions.baseline.params` to the calibrated `{kd, epsilon}`.
- [ ] **Step 3: Verify** the study passes structural lint and the audit gate:
```bash
.venv/bin/python -m viva_superpowers.study_audit --workspace . 2>/dev/null | grep recruitment-adaptive
.venv/bin/python -m viva_superpowers.study_audit --workspace . --gate; echo "gate=$?"
```
Expected: `recruitment-adaptive [pass]` (composite resolves; canonical model schema; DAG intact); `gate=0`. (Note: this study is a fresh member — leave `status` at an honest in-progress value like `draft`/`in-progress` until the loop runs it in Task 7; it is NOT a completed result yet.)
- [ ] **Step 4: Commit** the study.

---

### Task 7: The loop driver + emergent trajectory

**Files:**
- Create: `scripts/build_recruitment_loop.py`
- Test: `tests/test_recruitment_loop.py`
- Generates (committed): `workspace/investigations/recruitment-model-building/trajectory.json`, `.pbg/loop/recruitment-adaptive.json`

**Interfaces:**
- Consumes: `mechanisms.LIBRARY`/`simulate_condition` (Task 5), `navigate.next_mechanism` (Task 5), the contract (Task 6), `viva_superpowers.{loop_state as ls, test_audit, test_contract as tc}`.
- Produces: `main()` writing a `model_build_trajectory/v2` `trajectory.json`; helper `run_loop(library, *, tag, withhold=None)` returning `(outcome, iterations, final_mechanism)`.

- [ ] **Step 1: Write the driver test (small, deterministic)**

```python
# tests/test_recruitment_loop.py
import importlib.util, pathlib
spec = importlib.util.spec_from_file_location("brl", pathlib.Path("scripts/build_recruitment_loop.py"))
brl = importlib.util.module_from_spec(spec); spec.loader.exec_module(brl)

def test_loop_reaches_adaptive():
    outcome, iters, final = brl.run_loop(brl.LIBRARY, tag="test")
    assert outcome == "DONE"
    assert final == "adaptive_receptor"
    assert len(iters) >= 2                      # climbed at least static->hill->adaptive

def test_giveup_companion_terminates_open():
    # withhold adaptive_receptor -> loop cannot pass recruits_high -> honest GIVE_UP
    outcome, iters, final = brl.run_loop(brl.LIBRARY, tag="giveup", withhold="adaptive_receptor")
    assert outcome == "GIVE_UP"
    assert final != "adaptive_receptor"

def test_audit_two_rounds():
    # the one-sided draft is flagged; the two-sided locked spec passes
    import viva_superpowers.test_audit as ta
    assert ta.audit_gate(ta.build_audit_report(brl.draft_spec())) in ("warn", "fail")
    assert ta.audit_gate(ta.build_audit_report(brl.locked_spec())) == "pass"
```

- [ ] **Step 2: Run to verify fail.** Expected: FAIL (script missing).

- [ ] **Step 3: Implement `scripts/build_recruitment_loop.py`** mirroring `viva-casebook/scripts/build_loop_demo.py`'s structure:
  - `draft_spec()` — the contract with a ONE-SIDED `recruits_high` (`op: ">="`) → audit flags discrimination drift.
  - `locked_spec()` — the two-sided banded contract (matches Task 6's study.yaml).
  - `grade(mechanism, *, use_locked=True)` — run each condition via `mechanisms.simulate_condition`, build `report_card_verdict/v2` axes via `tc.check(...)` against the bands (recruits_low, receptor_gating, recruits_high, adaptation_knockout).
  - `run_loop(library, *, tag, withhold=None)` — start at `static_lambda`; each iteration grade → `navigate.next_mechanism` → install (or DONE); `ls.record_iteration` each step; DONE when all hard axes within_tol; GIVE_UP when no further mechanism available or budget spent. `withhold` removes a rung from the library (the honest GIVE_UP companion).
  - `main()`: `rep1=audit(draft)`, `rep2=audit(locked)`; `ls.create` → `advance("AUDIT", ...)` → `advance("SELECT", sourcing={"decision":"compose","modules":["chemotaxis_receptor"]})` → `lock_tests` → `run_loop` (real) → GIVE_UP companion run → `perturbation` (re-grade final under a ±10% epsilon / different seed) → `ls.validate` → `ls.save` + write `trajectory.json` (schema `model_build_trajectory/v2`: `contract, draft, final_composite, select, audit, control, lock, tests, iterations, result, giveup_companion, timeseries`). Keep seeds/steps small.

- [ ] **Step 4: Run the driver, then the tests**

```bash
.venv/bin/python scripts/build_recruitment_loop.py    # generates trajectory.json + .pbg/loop/*.json
.venv/bin/python -m pytest tests/test_recruitment_loop.py -q
```
Expected: driver prints the emergent climb and writes the artifacts; tests PASS. If `run_loop` does not reach DONE at `adaptive_receptor`, the mechanism params (Task 4) or the bands (Task 6) need reconciling — fix there, not by weakening the loop.

- [ ] **Step 5: Full-suite + audit + commit**

```bash
.venv/bin/python -m pytest tests/ -q
.venv/bin/python -m viva_superpowers.study_audit --workspace . --gate; echo "gate=$?"
git add scripts/build_recruitment_loop.py tests/test_recruitment_loop.py \
  workspace/investigations/recruitment-model-building/trajectory.json .pbg/loop/recruitment-adaptive.json
git commit -m "feat: recruitment model-building loop driver + emergent trajectory"
```
(Confirm `.pbg/loop/*.json` is not gitignored; if it is, commit only the trajectory under `workspace/`.)

---

## Self-review notes

- **Spec coverage (Part B):** B1 mechanism library → Tasks 1,2,3,5; the emergent climb → Task 5 (navigate) + Task 7 (run_loop). B2 contract → Task 6, with the audit two-rounds + prereg-lock in Task 7. B3 driver + trajectory + GIVE_UP + perturbation → Task 7. The empirical ladder (the spec's load-bearing claim) → Task 4, gated explicitly.
- **Deferred to Phase 3:** the flagship `recruitment-pipeline-report.html` renderer (B4) and the `recruitment-model-building/investigation.yaml` (B5). Task 7 writes the `trajectory.json` those consume, and creates the investigation DIRECTORY (for the trajectory), but the investigation.yaml + report are Phase 3.
- **Risk:** Task 4 is the gate. If the ladder doesn't separate empirically, Task 4 reports BLOCKED rather than faking it — the honest-termination principle applies to the plan itself.
- **Type consistency:** `AdaptiveReceptorSubcell.signal` returns `(theta, m_new, response)` (Task 1) consumed by its own `update` and reused by the calibration test; `adaptation_coupling` (Task 2) wires the `m_prev`/`m` ports Task 1 declares; `CONDITIONS` keys (Task 3) are the condition names the contract's `measure.condition` fields (Task 6) and `simulate_condition` (Task 5) reference; `navigate.next_mechanism(active, axes)` (Task 5) is called by `run_loop` (Task 7).
