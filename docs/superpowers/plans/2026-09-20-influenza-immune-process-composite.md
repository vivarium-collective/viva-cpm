# Influenza Immune-Layer Process + full_model Composite — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Promote the influenza capstone science from `run_full_model`'s monolithic Python loop into a process-bigraph `Composite` whose immune layer is a swappable, agent-based `ImmuneProcess` with immune cells living off the CPM lattice (removing the space cap that holds spatial immune counts ~2–4× below the fig3b bands).

**Architecture:** Four processes wired through shared stores and executed by `process_bigraph.Composite`: `EpitheliumProcess` (extends `CPMProcess`; lattice + 4 fields + epithelial cell-fates), `ImmuneProcess` (off-lattice agents: chemotax/proximity-kill/secrete/ODE-driven-recruit, uncapped), `SystemicODEProcess` (wraps `GlobalODE`), and a rewritten `full_model` document driven by a new `run_full_model_composite`. All four diffusive fields stay inside the Rust `World`; immune agents read/write them through new point sample/deposit ports.

**Tech Stack:** Rust (`crates/cpm-core`, `crates/cpm-py`, PyO3/maturin), Python 3.12, `process_bigraph` (`Process`, `Composite`), numpy, pytest.

**Spec:** `docs/superpowers/specs/2026-09-20-influenza-immune-process-composite-design.md`

## Global Constraints

- **No tuning of source constants to fit figure bands.** Only the spatial *realization* changes (off-lattice agents, proximity killing). Recruitment/kill/secretion/ODE constants keep their `params.yaml` / `resolve_constants` values. Source dossier (`docs/cc3d-reference/sego2022-*`) is authority.
- **The old plain-function path stays** (`run_full_model`, `killing.py`, `recruitment.py`) until the Composite passes the Fig-3B parity gate (Task 3.4). Do not delete it before then.
- **Fast test suite stays green** at every task boundary. Heavy paper-scale runs go on the Mac mini, not in pytest.
- **`kill_radius` and agent step sizes are derived from CPM geometry** (cell footprint = `sqrt(cell_sites)` = 5 sites), never fitted to outcomes. Document the mapping.
- **Worktree discipline:** work in a dedicated worktree off current `origin/main`; Rust rebuild after any `crates/` change via `VIRTUAL_ENV="$PWD/.venv" uvx maturin develop --release`.
- **No AI attribution** in commits or PRs.
- Field site indexing is flat: `site = x + y*nx + z*nx*ny` (see `cpm_process.py:54`).

---

## File Structure

- `crates/cpm-core/src/field.rs` — add `Field::value_at` / `Field::add_source_at` (core).
- `crates/cpm-core/src/lib.rs` (World) + `crates/cpm-py/src/lib.rs` — expose `field_value_at` / `field_add_source_at` on the Python `World`.
- `pbg_cpm_studies/influenza/epithelium_process.py` (NEW) — `EpitheliumProcess(CPMProcess)`: all-field per-cell readout, `field_at_point` output, `field_deposit` input, and the epithelial cell-fate transitions moved out of `run.py`.
- `pbg_cpm_studies/influenza/immune_process.py` (NEW) — `ImmuneProcess`: agent state, chemotaxis, proximity killing, secretion deposition, ODE-driven recruitment.
- `pbg_cpm_studies/influenza/ode_process.py` (NEW) — `SystemicODEProcess`: wraps `price_ode.GlobalODE`.
- `pbg_cpm_studies/composites/influenza.py` — rewrite `full_model_composite_document` (`:452`) to wire the four processes; add `EPITHELIUM_ADDR`/`IMMUNE_ADDR`/`ODE_ADDR`.
- `pbg_cpm_studies/influenza/run.py` — add `run_full_model_composite(...)` driving `process_bigraph.Composite`; leave `run_full_model` intact.
- `tests/test_epithelium_process.py`, `tests/test_immune_process.py`, `tests/test_ode_process.py`, `tests/test_full_model_composite.py` (NEW).

---

## Phase 1 — Epithelium substrate: Rust point I/O + EpitheliumProcess

### Task 1.1: Rust — field point read/deposit

**Files:**
- Modify: `crates/cpm-core/src/field.rs` (add methods to `Field` and the field-owning struct near `field_conc` at `:176`)
- Modify: `crates/cpm-py/src/lib.rs` (expose on `World` near `:177`)
- Test: `tests/test_cpm_field_point.py` (Create)

**Interfaces:**
- Consumes: existing `World.add_field(name, d, decay) -> usize`, `World.field_conc(idx) -> list[f32]`, `World.dims() -> (nx,ny,nz)`.
- Produces: `World.field_value_at(field_idx: int, x: int, y: int, z: int) -> float` (concentration at that site) and `World.field_add_source_at(field_idx: int, x: int, y: int, z: int, amount: float)` (adds `amount` to that site's concentration; clamps site to domain).

- [ ] **Step 1: Write the failing test**

```python
# tests/test_cpm_field_point.py
import cpm
from cpm.cpm_core import World

def test_field_value_and_deposit_at_point():
    w = World(10, 10, 1, neighbor_order=2, temperature=10.0, seed=1)
    idx = w.add_field("c", 0.0, 0.0)          # no diffusion, no decay
    w.finalize(1)
    assert w.field_value_at(idx, 3, 4, 0) == 0.0
    w.field_add_source_at(idx, 3, 4, 0, 2.5)
    assert abs(w.field_value_at(idx, 3, 4, 0) - 2.5) < 1e-9
    # out-of-range coordinates are clamped, not a panic
    w.field_add_source_at(idx, 99, 99, 0, 1.0)
    assert w.field_value_at(idx, 9, 9, 0) >= 1.0
```

(Confirm the exact `World(...)` constructor kwargs against `crates/cpm-py/src/lib.rs`; match the signature `build.world_from_spec` uses if these differ.)

- [ ] **Step 2: Run test to verify it fails** — `pytest tests/test_cpm_field_point.py -v` → FAIL (`field_value_at` not defined).

- [ ] **Step 3: Implement core methods** in `field.rs`. Add to `Field`:

```rust
pub fn value_at(&self, site: usize) -> f64 { self.conc[site] as f64 }
pub fn add_source_at(&mut self, site: usize, amount: f64) {
    self.conc[site] += amount as f32;
}
```

(Use the same `conc` storage `field_conc`/`diffuse_step` read/write, `field.rs:44,49,176`.) On the field-owning struct add, next to `field_conc` (`:176`):

```rust
pub fn field_value_at(&self, field_idx: usize, site: usize) -> f64 {
    self.fields[field_idx].value_at(site)
}
pub fn field_add_source_at(&mut self, field_idx: usize, site: usize, amount: f64) {
    self.fields[field_idx].add_source_at(site, amount);
}
```

- [ ] **Step 4: Expose on the Python `World`** in `crates/cpm-py/src/lib.rs` (near `:177`), clamping to the domain and computing the flat site:

```rust
fn field_value_at(&self, field_idx: usize, x: usize, y: usize, z: usize) -> f64 {
    let (nx, ny, nz) = self.inner.dims();
    let (x, y, z) = (x.min(nx-1), y.min(ny-1), z.min(nz-1));
    self.inner.field_value_at(field_idx, x + y*nx + z*nx*ny)
}
fn field_add_source_at(&mut self, field_idx: usize, x: usize, y: usize, z: usize, amount: f64) {
    let (nx, ny, nz) = self.inner.dims();
    let (x, y, z) = (x.min(nx-1), y.min(ny-1), z.min(nz-1));
    self.inner.field_add_source_at(field_idx, x + y*nx + z*nx*ny, amount);
}
```

(Match the actual struct/field names — `self.inner` vs direct — against the existing `field_conc` wrapper.)

- [ ] **Step 5: Rebuild + run test** — `VIRTUAL_ENV="$PWD/.venv" uvx maturin develop --release` then `pytest tests/test_cpm_field_point.py -v` → PASS.

- [ ] **Step 6: Commit**

```bash
git add crates/cpm-core/src/field.rs crates/cpm-py/src/lib.rs tests/test_cpm_field_point.py
git commit -m "cpm: field point read/deposit (field_value_at, field_add_source_at)"
```

### Task 1.2: EpitheliumProcess — extended CPM ports (all-field readout, point sample, point deposit)

**Files:**
- Create: `pbg_cpm_studies/influenza/epithelium_process.py`
- Test: `tests/test_epithelium_process.py` (Create)

**Interfaces:**
- Consumes: `CPMProcess` (`cpm/processes/cpm_process.py`), `World.field_value_at/field_add_source_at` (Task 1.1), `World.field_mean_at_cell(idx, cid)`, `World.cell_coms()`.
- Produces: `EpitheliumProcess(CPMProcess)` with `outputs()` adding `field_at_cell_all: overwrite[map[list]]` (per cell id → `[mean_field0, …, mean_fieldN-1]`) and keeping `positions`/`types`/`volumes`; `inputs()` adding `field_deposit: list` (each item `[field_idx, x, y, z, amount]`). This task does NOT yet move epithelial fates (Task 1.3).

- [ ] **Step 1: Write the failing test**

```python
# tests/test_epithelium_process.py
import numpy as np
from pbg_cpm_studies.influenza.epithelium_process import EpitheliumProcess
from pbg_cpm_studies.influenza import build, immune

def _spec():
    return immune.build_cytotoxic_scenario_spec(
        epithelial_cells_per_side=4, n_infected=1, n_macrophages=2,
        n_nk=2, n_cd8=2, margin_sites=10, separation_sites=8, seed=1)

def test_epithelium_process_exposes_all_fields_and_deposits():
    spec = _spec(); spec["fields"] = None  # process adds fields per config
    p = EpitheliumProcess({"spec": spec, "mcs_per_update": 1, "n_fields": 4})
    out = p.update({"fates": {}, "field_deposit": []}, 1.0)
    # per-cell readout is a list of 4 field means per cell
    any_cell = next(iter(out["field_at_cell_all"]))
    assert len(out["field_at_cell_all"][any_cell]) == 4
    # depositing chemokine (field idx 2) at a point then reading it back rises
    p.update({"fates": {}, "field_deposit": [[2, 20, 20, 0, 5.0]]}, 1.0)
    assert p.world.field_value_at(2, 20, 20, 0) > 0.0
```

- [ ] **Step 2: Run test to verify it fails** — `pytest tests/test_epithelium_process.py -v` → FAIL (module missing).

- [ ] **Step 3: Implement `EpitheliumProcess`** subclassing `CPMProcess`, overriding `outputs()`, `inputs()`, and `update()` to (a) apply `field_deposit` before `world.step`, (b) emit all-field per-cell means:

```python
from process_bigraph import Process
from cpm.processes.cpm_process import CPMProcess

class EpitheliumProcess(CPMProcess):
    def inputs(self):
        return {"fates": "map[integer]", "field_deposit": "list"}

    def outputs(self):
        out = dict(super().outputs())
        out["field_at_cell_all"] = "overwrite[map[list]]"
        # raw flat fields + dims so ImmuneProcess can sample gradients in numpy
        # (process-bigraph has no within-update request/response, so the whole
        # field crosses the store boundary once per update).
        out["chemo_field"] = "overwrite[list]"
        out["virus_field"] = "overwrite[list]"
        out["dims"] = "overwrite[list]"
        return out

    def update(self, state, interval):
        state = state or {}
        for item in state.get("field_deposit") or []:
            fidx, x, y, z, amt = item
            self.world.field_add_source_at(int(fidx), int(x), int(y), int(z), float(amt))
        base = super().update(state, interval)   # applies fates, steps, base outputs
        types = base["types"]; n = len(types)
        base["field_at_cell_all"] = {
            str(cid): [self.world.field_mean_at_cell(f, cid) for f in range(self.n_fields)]
            for cid in range(1, n)
        }
        # field indices (per fields.py order): 0=virus 1=ifn 2=chemokine 3=il10
        base["chemo_field"] = list(self.world.field_conc(2)) if self.n_fields > 2 else []
        base["virus_field"] = list(self.world.field_conc(0)) if self.n_fields > 0 else []
        base["dims"] = list(self.world.dims())
        return base
```

Also extend Task 1.2's test to assert these outputs exist:

```python
    assert len(out["dims"]) == 3
    assert len(out["virus_field"]) == out["dims"][0]*out["dims"][1]*out["dims"][2]
```

- [ ] **Step 4: Run test to verify it passes** — `pytest tests/test_epithelium_process.py -v` → PASS.

- [ ] **Step 5: Commit**

```bash
git add pbg_cpm_studies/influenza/epithelium_process.py tests/test_epithelium_process.py
git commit -m "influenza: EpitheliumProcess extends CPMProcess with all-field readout + point deposit"
```

### Task 1.3: EpitheliumProcess — epithelial cell-fate transitions moved in

**Files:**
- Modify: `pbg_cpm_studies/influenza/epithelium_process.py`
- Reference (verbatim mechanism source): `pbg_cpm_studies/influenza/run.py` `_one_mcs` steps 2–8.5 (`run.py:1975`–`2097`) and the helpers `transitions.*`, `resistance.cell_resistance`, `signaling.*`, `allee.*`, `fields.il10_hill_constants`, `price_ode.hill`.
- Test: `tests/test_epithelium_process.py` (extend)

**Interfaces:**
- Consumes: everything Task 1.2 produces + `state["ode_state"]` (a dict with at least `"X"` for ROS) and `state["immune_kills"]` (list of infected cell ids the immune layer killed this record).
- Produces: `EpitheliumProcess.update` now runs the ordered epithelial fate pipeline internally each MCS: infection H→I, per-cell resistance, secretion scales (virus + macrophage-independent uninfected IL-10 gate), infected apoptosis I→D, Allee death/recovery, ROS death via `ode_state["X"]`. Immune kills are merged as an additional fate source at defined priority (applied after infection, before apoptosis). Adds config `mcs_per_step` and per-mechanism `enable` set. Emits `counts` (H,I,D) in outputs.

- [ ] **Step 1: Write the failing test** — infection depletes living epithelium and dead accumulates, mirroring `tests/test_influenza_full_model.py::test_full_model_infection_depletes_uninfected`:

```python
def test_epithelium_process_fates_deplete_living_epithelium():
    from pbg_cpm_studies.influenza.run import load_params
    from pbg_cpm_studies.influenza import price_ode, types
    spec = _spec(); spec["fields"] = None
    p = EpitheliumProcess({"spec": spec, "mcs_per_update": 1, "n_fields": 4,
                           "enable": ["infection","ifn","death","ros","allee"]})
    # feed a large ROS X so ROS death fires
    living0 = None
    for _ in range(25):
        out = p.update({"fates": {}, "field_deposit": [],
                        "ode_state": {"X": 50.0}, "immune_kills": []}, 1.0)
        H, I, D = out["counts"]["H"], out["counts"]["I"], out["counts"]["D"]
        living0 = living0 if living0 is not None else H + I
    assert D > 0 and (H + I) < living0
```

- [ ] **Step 2: Run test to verify it fails** — `pytest tests/test_epithelium_process.py -k fates -v` → FAIL.

- [ ] **Step 3: Implement the fate pipeline.** Move the body of `_one_mcs` steps 2–8.5 (`run.py:1975`–`2097`) into a private `_epithelial_fates(self, ode_state, immune_kills)` called inside `update` after `field_deposit` and before/around `world.step`. Preserve the exact source order and helper calls; read constants once in `initialize` via `load_params()`/`resolve_constants` (num_epithelial = number of epithelial cells in the spec). Apply `immune_kills` (I→D, add to a `dead_from_infected` set) right after the infection step. Do NOT change any rate. Emit `counts` from `world.cell_types()`.

(This is a mechanical move of already-tested logic; keep the RNG seeding streams identical — `infect_rng`, `death_rng`, `allee_rng`, `ros_rng` — so parity Task 1.4 can match.)

- [ ] **Step 4: Run test to verify it passes** — `pytest tests/test_epithelium_process.py -v` → PASS.

- [ ] **Step 5: Commit**

```bash
git add pbg_cpm_studies/influenza/epithelium_process.py tests/test_epithelium_process.py
git commit -m "influenza: move epithelial cell-fate transitions into EpitheliumProcess (ordered, source-faithful)"
```

### Task 1.4: Phase-1 parity gate — EpitheliumProcess vs run_full_model (immune off)

**Files:**
- Test: `tests/test_epithelium_process.py` (extend)

**Interfaces:**
- Consumes: `run.run_full_model(..., enable=<epithelial-only subset>)` and a bare `EpitheliumProcess` stepped the same number of MCS with `immune_kills=[]` and the same `ode_state["X"]` schedule.

- [ ] **Step 1: Write the parity test** — with immune subsystems disabled and a fixed injected `X` schedule (or ODE disabled → X=0), the epithelial series (H,I,D) from `EpitheliumProcess` driven mcs-for-mcs must match `run_full_model(enable=epithelial-only)` within a small tolerance for the same seed:

```python
def test_epithelium_process_parity_with_run_full_model_epithelial_only():
    import numpy as np
    from pbg_cpm_studies.influenza import run
    enable = ["infection","ifn","death","allee"]   # ros/ode/immune off for a clean parity
    r = run.run_full_model(cells_per_side=8, steps=20, seed=3,
                           init_infection_frac=0.1, enable=enable)
    # build the same scene + step EpitheliumProcess identically; compare dead[-1]
    # within tolerance (see helper _drive_epithelium_process below)
    dead_ref = r["counts"]["dead"][-1]
    dead_proc = _drive_epithelium_process(cells_per_side=8, steps=20, seed=3,
                                          init_infection_frac=0.1, enable=enable)
    assert abs(dead_proc - dead_ref) <= max(2, int(0.1 * max(dead_ref, 1)))
```

(Write `_drive_epithelium_process` in the test module: build the scene via the same `immune.build_cytotoxic_scenario_spec` + scatter that `run_full_model` uses, construct `EpitheliumProcess`, step `mcs_per_step` MCS per record, return final dead count. The tolerance encodes the documented parallel-within-MCS approximation; if `enable` has no cross-process interaction the match should be near-exact.)

- [ ] **Step 2: Run** — `pytest tests/test_epithelium_process.py -k parity -v`. If it fails beyond tolerance, reconcile RNG stream ordering (Task 1.3) until within tolerance; record the achieved delta in a comment.

- [ ] **Step 3: Run the full fast suite** — `pytest tests/ -q -k "influenza or recruit or cpm"` → all pass.

- [ ] **Step 4: Commit**

```bash
git add tests/test_epithelium_process.py
git commit -m "influenza: Phase-1 parity gate — EpitheliumProcess matches run_full_model epithelial-only"
```

---

## Phase 2 — Agent-based ImmuneProcess

### Task 2.1: ImmuneProcess — agent state + chemotaxis

**Files:**
- Create: `pbg_cpm_studies/influenza/immune_process.py`
- Test: `tests/test_immune_process.py` (Create)

**Interfaces:**
- Consumes: `state["field_at_point"]` is NOT used (agents sample via a callback is impossible across the store boundary); instead the Composite wires `field_samples` — see below. To keep the process pure, `ImmuneProcess` receives a **sampled chemokine gradient per agent** is too coupled; simpler: the process holds agent positions and emits `sample_requests`, but process-bigraph has no request/response within one update. **Resolution:** `ImmuneProcess` reads the whole chemokine field once per update via `state["chemo_field"]` (the flat `list[f32]` from `World.field_conc(2)`) plus `state["dims"]`, and samples/gradients itself in numpy. `EpitheliumProcess` therefore also outputs `chemo_field: overwrite[list]` and `dims: overwrite[list]`.
- Produces: `ImmuneProcess(Process)` with `config_schema` {`chemotaxis`: params, `seed`}, agent store `immune_agents` (list of `{id,type,x,y}`), and `update` that moves each agent one bounded step up the local chemokine gradient (macrophages up virus — add `virus_field` similarly). `inputs()`: `{"chemo_field":"list","virus_field":"list","dims":"list","immune_agents":"list"}`; `outputs()`: `{"immune_agents":"overwrite[list]"}`.

- [ ] **Step 1: Write the failing test** — an agent in a monotone chemokine gradient moves up-gradient:

```python
# tests/test_immune_process.py
import numpy as np
from pbg_cpm_studies.influenza.immune_process import ImmuneProcess
from pbg_cpm_studies.influenza import types

def _linear_field(nx, ny):
    # increasing in +x
    return [float(x) for _ in range(ny) for x in range(nx)]  # site=x+y*nx

def test_nk_agent_moves_up_chemokine_gradient():
    nx, ny = 30, 10
    p = ImmuneProcess({"seed": 1})
    agents = [{"id": 1, "type": int(types.K), "x": 10.0, "y": 5.0}]
    out = p.update({"chemo_field": _linear_field(nx, ny),
                    "virus_field": [0.0]*(nx*ny),
                    "dims": [nx, ny, 1], "immune_agents": agents}, 1.0)
    assert out["immune_agents"][0]["x"] > 10.0   # moved toward higher chemokine
```

- [ ] **Step 2: Run to verify it fails** — `pytest tests/test_immune_process.py -k gradient -v` → FAIL.

- [ ] **Step 3: Implement agent state + gradient step.** In `update`, reshape `chemo_field` to `(ny,nx)`, compute the local gradient at each agent's nearest site (central difference, clamped at edges), and move the agent `step_len * sign(grad)` with a small random component (`self.rng`). `step_len` from config (default 1.0 site). Macrophages (`types.M`) use `virus_field`; K/E use `chemo_field`. Keep agents within `[0,nx-1]×[0,ny-1]`.

- [ ] **Step 4: Run to verify it passes** — PASS.

- [ ] **Step 5: Commit**

```bash
git add pbg_cpm_studies/influenza/immune_process.py tests/test_immune_process.py
git commit -m "influenza: ImmuneProcess agent state + chemotaxis (field-gradient step)"
```

### Task 2.2: ImmuneProcess — proximity killing + secretion deposition

**Files:**
- Modify: `pbg_cpm_studies/influenza/immune_process.py`
- Reference: `pbg_cpm_studies/influenza/killing.py` (`contact_kill_rate`, `nearby_kill_rate`), `signaling.macrophage_secretion_scale`, `fields.il10_hill_constants`.
- Test: `tests/test_immune_process.py` (extend)

**Interfaces:**
- Consumes: `state["epithelial_positions"]` (map cell id → `[x,y]`, from `EpitheliumProcess.positions` filtered to infected), `state["infected_ids"]` (list), `state["sig_1"]` (float from the ODE record), the fields.
- Produces: `outputs()` gains `immune_kills: overwrite[list]` (infected cell ids killed this update) and `field_deposit: overwrite[list]` (chemokine idx 2 + IL-10 idx 3 point sources at macrophage agent positions). Killing uses the source rate at proximity: for each K/E agent within `kill_radius` (= `sqrt(cell_sites)` ≈ 5) of an infected cell, apply per-MCS `1 - exp(-rate)` with `rate` from `killing.contact_kill_rate` using a unit contact-area proxy (document the mapping; do not fit). Secretion: each macrophage agent deposits `b_c_per_site * scale` chemokine and `b_l_per_site * scale` IL-10, where `scale = macrophage_secretion_scale(L_at_agent, sig_1, g_1, g_2, d_2)` — the SAME function `run.py` uses (`fields.il10_hill_constants()`), constants unchanged.

- [ ] **Step 1: Write the failing tests** — (a) an NK agent adjacent to an infected cell can kill it (non-zero kill probability over many draws); (b) a macrophage agent with positive `sig_1` emits a positive chemokine `field_deposit`.

```python
def test_nk_agent_kills_adjacent_infected():
    p = ImmuneProcess({"seed": 2, "kill_radius": 5.0})
    agents = [{"id": 1, "type": int(types.K), "x": 10.0, "y": 5.0}]
    killed_any = False
    for _ in range(200):
        out = p.update({"chemo_field":[0.0]*300,"virus_field":[0.0]*300,"dims":[30,10,1],
                        "immune_agents": agents,
                        "epithelial_positions": {"7":[11.0,5.0]}, "infected_ids":[7],
                        "sig_1": 0.0}, 1.0)
        if 7 in out["immune_kills"]:
            killed_any = True; break
    assert killed_any

def test_macrophage_deposits_chemokine_when_sig1_positive():
    p = ImmuneProcess({"seed": 3})
    agents = [{"id": 1, "type": int(types.M), "x": 10.0, "y": 5.0}]
    out = p.update({"chemo_field":[0.0]*300,"virus_field":[0.0]*300,"dims":[30,10,1],
                    "immune_agents": agents, "epithelial_positions": {}, "infected_ids": [],
                    "sig_1": 0.02}, 1.0)
    deposits = [d for d in out["field_deposit"] if d[0] == 2]  # chemokine field idx
    assert deposits and deposits[0][4] > 0.0
```

- [ ] **Step 2: Run to verify they fail** → FAIL.

- [ ] **Step 3: Implement killing + secretion** per the interface above; read `g_ik`/`g_ie`/`tot_ec`/`cell_volume` and the IL-10 Hill constants from `load_params()`/`il10_hill_constants()` in `initialize`. Document the contact-area proxy choice in the docstring.

- [ ] **Step 4: Run to verify they pass** → PASS.

- [ ] **Step 5: Commit**

```bash
git add pbg_cpm_studies/influenza/immune_process.py tests/test_immune_process.py
git commit -m "influenza: ImmuneProcess proximity killing + macrophage chemokine/IL-10 deposition (source rates, proxy geometry documented)"
```

### Task 2.3: ImmuneProcess — ODE-driven recruitment (uncapped spawn/remove)

**Files:**
- Modify: `pbg_cpm_studies/influenza/immune_process.py`
- Reference: `pbg_cpm_studies/influenza/recruitment.py` (`macrophage_inflow`, `nk_inflow`, `cd8_inflow`, `*_outflow`, `ul_rate_to_prob`, `poisson_inflow_count`).
- Test: `tests/test_immune_process.py` (extend)

**Interfaces:**
- Consumes: `state["recruit_drivers"]` = `{"macro_inflow":…, "nk_inflow":…, "cd8_inflow":…, "macro_outflow":…, "nk_outflow":…, "cd8_outflow":…}` (rates from the ODE record), and `state["margin_box"]` = `[x0,y0,x1,y1]` (where to spawn). Recruitment runs when `state.get("apply_recruitment")` is truthy (the Composite sets it once per record so per-MCS cadence is applied by looping — see Task 3.3).
- Produces: `update` spawns `poisson_inflow_count(local_ratio*inflow)` new agents per type at random margin positions (new incrementing ids) and removes each existing agent of that type with prob `ul_rate_to_prob(local_ratio*outflow)`. **No reserve-pool cap.** `local_ratio` from `params["coupling"]["recruitment"]["local_ratios"]`. With agents uncapped, default `local_ratio → 1.0` for all types (drop the nearby surrogate) — record this in the docstring as the design decision from the spec.

- [ ] **Step 1: Write the failing test** — a large macrophage inflow driver grows the agent population well beyond any fixed pool:

```python
def test_recruitment_grows_agents_uncapped():
    p = ImmuneProcess({"seed": 4})
    agents = [{"id": 1, "type": int(types.M), "x": 5.0, "y": 5.0}]
    drivers = {"macro_inflow": 5.0, "nk_inflow": 0.0, "cd8_inflow": 0.0,
               "macro_outflow": 0.0, "nk_outflow": 0.0, "cd8_outflow": 0.0}
    for _ in range(50):
        out = p.update({"chemo_field":[0.0]*400,"virus_field":[0.0]*400,"dims":[20,20,1],
                        "immune_agents": agents, "epithelial_positions": {}, "infected_ids": [],
                        "sig_1": 0.0, "recruit_drivers": drivers, "margin_box": [1,1,19,19],
                        "apply_recruitment": True}, 1.0)
        agents = out["immune_agents"]
    n_macro = sum(1 for a in agents if a["type"] == int(types.M))
    assert n_macro > 100   # far beyond the old pool_per_type cap
```

- [ ] **Step 2: Run to verify it fails** → FAIL.

- [ ] **Step 3: Implement recruitment** per interface; reuse `recruitment.poisson_inflow_count`/`ul_rate_to_prob` verbatim (do not re-derive). Assign new agent ids from a monotone counter held on the process.

- [ ] **Step 4: Run to verify it passes** → PASS.

- [ ] **Step 5: Commit**

```bash
git add pbg_cpm_studies/influenza/immune_process.py tests/test_immune_process.py
git commit -m "influenza: ImmuneProcess ODE-driven recruitment (uncapped agent spawn/remove, local_ratio=1)"
```

---

## Phase 3 — SystemicODEProcess, Composite, and end-to-end parity

### Task 3.1: SystemicODEProcess

**Files:**
- Create: `pbg_cpm_studies/influenza/ode_process.py`
- Reference: `pbg_cpm_studies/influenza/price_ode.py` (`GlobalODE.__init__(consts, *, num_epithelial)`, `GlobalODE.step(state, inputs, dt_seconds) -> dict`, `INTEGRATED_STATES`, `resolve_constants`), and `run.py` `_ode_couple` (`run.py:2098`–`2128`) for the input assembly + the `recruit_drivers` derivation.
- Test: `tests/test_ode_process.py` (Create)

**Interfaces:**
- Consumes: `state` with spatial aggregates `H,I,M,K,E,DH` (ints), field integrals `V,F,C,L` (floats), and `B_ei,G_ki` (floats). Config: `consts` (from `resolve_constants`), `num_epithelial`, `dt_seconds`.
- Produces: `outputs()` = `{"ode_state":"overwrite[map[float]]","recruit_drivers":"overwrite[map[float]]","sig_1":"overwrite[float]"}`. `ode_state` carries the 10 integrated species; `recruit_drivers` carries the six rates (via `recruitment.*` at the freshly-integrated `C`,`P`); `sig_1 = a_11*T + a_12*D`.

- [ ] **Step 1: Write the failing test** — stepping with an APC-producing input advances `P` and yields positive macrophage inflow:

```python
# tests/test_ode_process.py
from pbg_cpm_studies.influenza.ode_process import SystemicODEProcess
from pbg_cpm_studies.influenza import price_ode
from pbg_cpm_studies.influenza.run import load_params

def test_ode_process_advances_and_emits_drivers():
    p = load_params()
    consts = price_ode.resolve_constants(p["price_ode"], num_epithelial=1225)
    proc = SystemicODEProcess({"consts": consts, "num_epithelial": 1225, "dt_seconds": 420.0})
    out = proc.update({"H":1000,"I":200,"M":10,"K":5,"E":2,"DH":0,
                       "V":50.0,"F":0.0,"C":0.5,"L":0.0,"B_ei":0.0,"G_ki":0.0}, 1.0)
    assert set(out["ode_state"]).issuperset(set(price_ode.INTEGRATED_STATES))
    assert "macro_inflow" in out["recruit_drivers"]
    assert out["recruit_drivers"]["macro_inflow"] > 0.0
```

- [ ] **Step 2: Run to verify it fails** → FAIL.

- [ ] **Step 3: Implement** by lifting the input-assembly + `GlobalODE.step` + `recruitment.*` rate computation from `_ode_couple` (`run.py:2098`–`2128`); hold `state` across updates on the process (`self.state`, initialized from the ODE ICs in `initialize`).

- [ ] **Step 4: Run to verify it passes** → PASS.

- [ ] **Step 5: Commit**

```bash
git add pbg_cpm_studies/influenza/ode_process.py tests/test_ode_process.py
git commit -m "influenza: SystemicODEProcess wraps GlobalODE (ode_state + recruit_drivers + sig_1)"
```

### Task 3.2: full_model Composite document

**Files:**
- Modify: `pbg_cpm_studies/composites/influenza.py` (`full_model_composite_document`, `:452`; add address constants near `:44`)
- Test: `tests/test_full_model_composite.py` (Create)

**Interfaces:**
- Consumes: `EpitheliumProcess`, `ImmuneProcess`, `SystemicODEProcess`.
- Produces: `full_model_composite_document(...)` returning a document with nodes `epithelium`, `immune`, `ode` and stores `fates`, `field_deposit`, `immune_agents`, `chemo_field`, `virus_field`, `dims`, `positions`, `types`, `epithelial_positions`, `infected_ids`, `ode_state`, `recruit_drivers`, `sig_1`, `immune_kills`, `margin_box`. Wire per the spec's data-flow diagram. `EPITHELIUM_ADDR = "local:!pbg_cpm_studies.influenza.epithelium_process.EpitheliumProcess"` etc.

- [ ] **Step 1: Write the failing test** — the document builds a `process_bigraph.Composite` and steps once without error, and all three process nodes are present:

```python
# tests/test_full_model_composite.py
from pbg_cpm_studies.composites.influenza import full_model_composite_document
from pbg_cpm_studies.core import build_core
from process_bigraph import Composite

def test_full_model_composite_builds_and_steps():
    doc = full_model_composite_document(cells_per_side=6, seed=1, init_infection_frac=0.1)
    assert set(doc).issuperset({"epithelium", "immune", "ode"})
    core = build_core()
    comp = Composite({"state": doc}, core=core)
    comp.run(1)   # one interval, no exception
```

- [ ] **Step 2: Run to verify it fails** → FAIL.

- [ ] **Step 3: Implement** the document. Reuse `_cpm_store`'s scene-building for the epithelium node config but point its address at `EPITHELIUM_ADDR`; build the scene via the same `immune.build_cytotoxic_scenario_spec` + scatter that `run_full_model` uses; seed `immune_agents` with the initial M/K/E agents; set `margin_box` from the scene dims. Confirm `build_core()` resolves `local:!…` addresses (see `pbg_cpm_studies/core.py`).

- [ ] **Step 4: Run to verify it passes** → PASS.

- [ ] **Step 5: Commit**

```bash
git add pbg_cpm_studies/composites/influenza.py tests/test_full_model_composite.py
git commit -m "influenza: full_model composite document wires Epithelium+Immune+ODE processes"
```

### Task 3.3: run_full_model_composite driver

**Files:**
- Modify: `pbg_cpm_studies/influenza/run.py` (add `run_full_model_composite`; do not touch `run_full_model`)
- Test: `tests/test_full_model_composite.py` (extend)

**Interfaces:**
- Consumes: `full_model_composite_document` (Task 3.2), `process_bigraph.Composite`.
- Produces: `run_full_model_composite(*, cells_per_side, steps, seed, init_infection_frac=None, init_viral_load=None, mcs_per_step=7, s_per_mcs=60.0) -> dict` returning the SAME result-dict shape as `run_full_model` (`counts`, `fields`, `ode`, `t_days`, `params`). Drives the Composite `mcs_per_step` MCS per record (Epithelium+Immune each MCS; ODE once per record — set `apply_recruitment` so Immune applies drivers each MCS), records observables from the stores.

- [ ] **Step 1: Write the failing test** — the driver returns the expected keys and a declining uninfected series under infection:

```python
def test_run_full_model_composite_shape_and_decline():
    from pbg_cpm_studies.influenza import run
    r = run.run_full_model_composite(cells_per_side=8, steps=12, seed=1, init_infection_frac=0.1)
    for k in ("counts","fields","ode","t_days","params"):
        assert k in r
    u = r["counts"]["uninfected"]
    assert u[-1] <= u[0]
```

- [ ] **Step 2: Run to verify it fails** → FAIL.

- [ ] **Step 3: Implement** the driver: build the document, construct `Composite`, loop `steps` records × `mcs_per_step` MCS, reading `types`/fields/`ode_state` from the composite state to fill the result dict (immune counts from `immune_agents` by type). Match the `t_days`/`params` conventions of `run_full_model` (`run.py:2147`).

- [ ] **Step 4: Run to verify it passes** → PASS.

- [ ] **Step 5: Commit**

```bash
git add pbg_cpm_studies/influenza/run.py tests/test_full_model_composite.py
git commit -m "influenza: run_full_model_composite drives the pb Composite (same result shape)"
```

### Task 3.4: End-to-end Fig-3B parity gate + immune-count validation

**Files:**
- Test: `tests/test_full_model_composite.py` (extend)
- Doc: `docs/cc3d-reference/chemokine-recruitment-scale-analysis.md` (append the composite result)

**Interfaces:**
- Consumes: `run.run_full_model_composite` and `run.run_full_model` (reference).

- [ ] **Step 1: Write the parity + validation tests.** (a) At reduced scale with immune subsystems weak, epithelial (H,I,D) from the composite matches `run_full_model` within tolerance (documented parallel-within-MCS bound). (b) The composite's macrophage count is NOT capped at the old ceiling — with recruitment on, it exceeds the reserve-pool number at a small scale.

```python
def test_composite_epithelial_parity_reduced_scale():
    from pbg_cpm_studies.influenza import run
    ref = run.run_full_model(cells_per_side=8, steps=15, seed=2, init_infection_frac=0.1)
    got = run.run_full_model_composite(cells_per_side=8, steps=15, seed=2, init_infection_frac=0.1)
    d_ref, d_got = ref["counts"]["dead"][-1], got["counts"]["dead"][-1]
    assert abs(d_got - d_ref) <= max(3, int(0.15 * max(d_ref, 1)))

def test_composite_immune_not_pool_capped():
    from pbg_cpm_studies.influenza import run
    got = run.run_full_model_composite(cells_per_side=12, steps=40, seed=1, init_infection_frac=0.05)
    assert max(got["counts"]["macrophage"]) >= 12   # grows with the tissue, not a fixed pool
```

- [ ] **Step 2: Run** — `pytest tests/test_full_model_composite.py -v`. Reconcile until parity holds within the stated tolerance; record the achieved delta in the test and the doc.

- [ ] **Step 3: Run the full fast suite** — `pytest tests/ -q -k "influenza or recruit or cpm or immune or epithelium or ode"` → all pass.

- [ ] **Step 4: Paper-scale immune check (mini, not pytest).** On the mini worktree, run `run_full_model_composite(cells_per_side=35, steps=720, seed=0, init_infection_frac=0.05)` and confirm macrophage/NK/CD8 rise materially above the ~176/196/188 ceilings from the pre-composite run; append the trajectory to `docs/cc3d-reference/chemokine-recruitment-scale-analysis.md`.

- [ ] **Step 5: Commit**

```bash
git add tests/test_full_model_composite.py docs/cc3d-reference/chemokine-recruitment-scale-analysis.md
git commit -m "influenza: Phase-3 parity gate + paper-scale immune-count validation (composite)"
```

---

## Self-Review notes (author)

- **Spec coverage:** EpitheliumProcess (Tasks 1.1–1.3), ImmuneProcess agents/kill/secrete/recruit (2.1–2.3), SystemicODEProcess (3.1), Composite + driver (3.2–3.3), parity gate + immune validation (1.4, 3.4), field-ownership-in-World (Task 1.1/1.2 point I/O), documented parallel-within-MCS ordering (1.4, 3.4 tolerances), no-tuning + old-path-retained (Global Constraints). Covered.
- **Cross-process field sampling:** resolved by having `EpitheliumProcess` also output `chemo_field`/`virus_field`/`dims` (flat `list[f32]`) so `ImmuneProcess` samples in numpy (no within-update request/response). Add these two outputs in Task 1.2's `outputs()` (extend the test there) — noted here so the Task-3.2 wiring has them.
- **kill_radius / step_len / secretion-per-site** are geometry-derived constants (Global Constraints), not fitted.
- **Open follow-up (not blocking):** once parity holds, migrate `repro_fig3b/5/7` to an opt-in `engine="composite"` and eventually retire the hand loop — a separate plan.
