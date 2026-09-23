# Influenza (Sego 2022) — Increment 8: Global ODE compartment + cellularization coupling

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add the Price-2015 global (non-spatial) ODE compartment as a hybrid coupled to the spatial CPM model, and use it to resolve the three stubs carried since Increments 5–7 — `sig_1` (TNF-driven signaling), immune-cell recruitment, and the nearby-population killing term — producing a `global-coupling` study that shows the integrated systemic species (TNF, ROS, antibody, APC, IL-12, IFN-γ, CD4⁺, B, neutrophils) evolve consistently with the spatial infection.

**Architecture:** The CC3D source runs a **hybrid**, not a standalone full ODE (dossier discrepancy #9). A libRoadRunner ODE integrates only the 10 *systemic* species `{NB, N, T, X, A, B, P, W, G, O}` (+ nearby surrogates `M_nb,K_nb,E_nb`); the spatial/field species `{H, I, M, E, K, L, C, F, V, DH}` are overwritten from CPM/field state every MCS. We reproduce this hybrid in Python: a `GlobalODEProcess` (scipy) integrated once per MCS (Δt = 1 min), fed by spatial aggregates, feeding back the three stubs. Additive to Increments 0–7; no Rust change (all primitives already exist).

**Tech Stack:** Python, scipy `solve_ivp` (LSODA/BDF, matching RoadRunner's stiff CVODE), numpy; existing `viva_cpm_studies/influenza/*` package; process-bigraph composite; matplotlib(Agg) viz.

**Spec:** `docs/superpowers/specs/2026-09-18-influenza-sego2022-reproduction-design.md` (Increment 8 row).

**Source authority (parameter dossier, in-repo):** `docs/cc3d-reference/sego2022-global-ode.md` — the verbatim ODE system (§2.3), the full ~90-constant table with scaling tokens (§2.4), initial conditions (§2.5), η/θ scaling (§3), and the bidirectional coupling map (§4). **Every numeric value and formula in this plan is sourced there; implementers read the dossier section cited, never guess.**

## Global Constraints

- **CC3D ViralInfectionVTM source is authoritative** over the paper for parameters/forms. The dossier is the transcription; cite `dossier §N` (which itself cites `file:symbol`).
- **Hybrid, not full ODE (discrepancy #9).** Integrate ONLY `{NB, N, T, X, A, B, P, W, G, O}` + nearby surrogates. The spatial species `{H, I, M, E, K, L, C, F, V, DH}` are *inputs* set from spatial state each step — never integrated by our process.
- **Use the SPATIAL model string's scalings (discrepancy #10).** Where the reference full-ODE string differs, use `immune_model_string` (spatial): `g_2 = base·s_v·s_l`, `d_2 = base·s_l`, `a_v = base/s_l`, `g_hv,g_vh,g_fi` carry `/s_l`, `b_h = base·s_t` (no `/s_v`), `R := F/(a_rf/s_l·s_v + F)`. The reference scalings appear only in an optional reference run.
- **η = `num_epithelial/250000`**, **θ = `1/250000/cell_volume = 1.6e-7`**, `s_t = s_to_mcs/86400 = 6.9444e-4` (Δt=1 min/MCS). η is derived from the domain, NOT a literal. Studies run at reduced scale (0.3 mm / ~900–1225 cells) per the ladder.
- **sig_1 is a saturating IL-10-inhibited secretion, not a linear source (discrepancy #11).** `Sigma1 = a_11·T + a_12·D` with `T` from the ODE and `D = tot_cell − H − I` from spatial counts; it drives the *existing* `signaling.macrophage_secretion_scale` Michaelis form. Do not change that functional form — only replace the stubbed constant input with the dynamic value.
- **Asymmetric recruitment (discrepancy #12).** CD8⁺ inflow is APC(`P`)-driven with NO homeostatic baseline; macrophage/NK are chemokine(`C`)-driven with a `mu·b` baseline. NK/CD8 outflow carry extra resistance-weighted `G_ki`/`B_ei` terms. Do NOT symmetrize.
- **Nearby-killing multiplies `cell_resist` DIRECTLY (discrepancy #7), not `(1−resist)`** — same as the Increment-7 local term. Preserve verbatim; flag in code + report.
- **Honesty:** the `global-coupling` study is `verdict: documented`, `biological_validation: PENDING` (quantitative Fig 3B reproduction is Increment 9). Report what the hybrid demonstrates (systemic species evolve, stubs resolved, loop closes) and what it does NOT (no quantitative Fig-3B match yet; any calibration gaps flagged for Incr 9). Keep ALL caveats visible.
- **No AI attribution** in commits/PRs. Work ONLY in worktree `/Users/eranagmon/code/viva-cpm--influenza-incr8`; NEVER touch canonical `/Users/eranagmon/code/viva-cpm`. Do NOT edit `viva_cpm_studies/visualizations/` or `viva_cpm_studies/composites/__init__.py` (peer-owned; composites auto-register via `@composite_generator`).
- **Per-worktree venv:** tests run with `/Users/eranagmon/code/viva-cpm--influenza-incr8/.venv/bin/python`. No Rust rebuild expected.

---

## File Structure

- `docs/cc3d-reference/sego2022-global-ode.md` — **already written** (the dossier); commit it in Task 8.0.
- `viva_cpm_studies/influenza/params.yaml` — ADD `price_ode:` + `coupling:` blocks (Task 8.0).
- `viva_cpm_studies/influenza/price_ode.py` — NEW. Pure ODE RHS for the 10 integrated states + nearby surrogates, a scaled-constant loader, ICs, and a `GlobalODE` integrator wrapper (Task 8.1).
- `viva_cpm_studies/influenza/recruitment.py` — NEW. Pure Hill inflow/outflow rate functions per immune type (Task 8.4).
- `viva_cpm_studies/influenza/killing.py` — MODIFY. Add `nearby_kill_rate` (well-mixed term) beside the Incr-7 `contact_kill_rate` (Task 8.5).
- `viva_cpm_studies/influenza/run.py` — MODIFY. Add `run_global_coupling(...)` driver wiring spatial→ODE→spatial each MCS (Tasks 8.2–8.5).
- `viva_cpm_studies/influenza/viz.py` — MODIFY. Add `global_coupling_figure` (Task 8.6).
- `viva_cpm_studies/composites/influenza.py` — MODIFY. Add a `global_coupling` live-demo composite (auto-registers) (Task 8.6).
- `tests/test_influenza_price_ode.py` — NEW (Tasks 8.1–8.3).
- `tests/test_influenza_recruitment.py` — NEW (Task 8.4).
- `tests/test_influenza_killing.py` — MODIFY (Task 8.5, add nearby-term tests; keep fast).
- `tests/test_influenza_viz.py` — MODIFY (Task 8.6).
- `tests/test_influenza_params.py` — MODIFY (Task 8.0).
- `workspace/studies/global-coupling/study.yaml` — NEW hand-authored study (Task 8.6).
- `workspace/investigations/influenza-sego2022/investigation.yaml` — MODIFY, add member (Task 8.6).

---

### Task 8.0: params `price_ode:` + `coupling:` blocks + dossier + Price-2015 reference

**Files:**
- Modify: `viva_cpm_studies/influenza/params.yaml`
- Modify: `tests/test_influenza_params.py`
- Add (already on disk): `docs/cc3d-reference/sego2022-global-ode.md`
- Reference: add Price et al. 2015 to the bib if the workspace uses one (check `workspace/*.bib` / expert-doc bib; if none applies, record the citation in the `price_ode.source` comment and skip).

**Interfaces:**
- Produces: `params.load_params()["price_ode"]` (dict with `integrated_states`, `spatialized_states`, `hill_exponents`, `integrator`, `scaling`, and a `constants` sub-map of `{name: {base, scale}}`), and `params.load_params()["coupling"]` (dict with `spatial_to_ode`, `sig_1`, `recruitment`, `nearby_killing`, `resistance`, `other_ode_to_spatial`). Consumed by Tasks 8.1–8.5.

- [ ] **Step 1: Failing test** — add to `tests/test_influenza_params.py`:

```python
def test_params_carry_price_ode_and_coupling():
    p = params.load_params()
    po = p["price_ode"]
    assert po["integrated_states"] == ["NB", "N", "T", "X", "A", "B", "P", "W", "G", "O"]
    assert set(po["spatialized_states"]) == {"H","I","M","E","K","L","C","F","V","DH"}
    assert po["hill_exponents"] == {"h_m":3,"h_x":2,"h_k":2,"h_e":3,"h_o":2}
    assert po["integrator"]["seconds_per_mcs"] == 60
    assert po["scaling"]["ode_epithelial_population"] == 250000
    # a handful of load-bearing constants carried as {base, scale}
    c = po["constants"]
    assert c["a_11"]["base"] == 0.00061091213762049 and c["a_11"]["scale"] is None
    assert c["g_ik"]["base"] == 0.0000308183563969158 and c["g_ik"]["scale"] == "s_t"
    cp = p["coupling"]
    assert cp["sig_1"]["definition"].startswith("Sigma1 = a_11*T + a_12*D")
    assert cp["recruitment"]["cd8"]["field"] == "apc_P"          # discrepancy #12
    assert "1-resist" not in cp["nearby_killing"]["note"].lower().replace(" ", "")  # resist DIRECT (#7)
```

- [ ] **Step 2: Run, verify it fails** (KeyError on `price_ode`).

- [ ] **Step 3: Add the blocks to `params.yaml`.** Transcribe from dossier §2.4 (full constant table), §2.5 (ICs), §3 (scaling), §4 (coupling), §7 (the proposed YAML). Requirements:
  - `price_ode.integrated_states`, `.spatialized_states`, `.integrated_nearby_surrogates`, `.hill_exponents` exactly as §7.
  - `price_ode.integrator`: `solver: roadrunner_cvode` (we approximate with scipy LSODA — note in a comment), `step_size: 1`, `steps_per_mcs: 1`, `seconds_per_mcs: 60`.
  - `price_ode.scaling`: `ode_epithelial_population: 250000`, `eta_formula: "num_epithelial / 250000"`, `eta_at_175: 0.0049`, `eta_at_500_1mm: 0.04`, `theta: 1.6e-7`, `s_t_formula: "s_to_mcs/86400"`, `cell_volume: 25`, `loc_to_global: "s_v/s_l"`.
  - `price_ode.constants`: EVERY constant in dossier §2.4 as `NAME: {base: <value>, scale: <token-or-null>}` where `scale` ∈ {`null`, `s_t`, `s_v`, `s_l`, `s_t*s_v`, `s_t/s_v`, `s_t/s_l`, `s_v*s_l`, `s_t*s_l/s_v`, `s_t/s_v/s_v`}. Split each `base * factor…` from §2.4 into `base` and the `scale` token string. Templated values (`b_m,b_k,b_p,V0`) → `scale: templated` with a `formula` note (from §2.5/§7).
  - `price_ode.homeostatic_pops_ode`: the five values in §2.5 (macro/nk/apc/neutro_blood/neutro), with a comment `# multiply by eta`.
  - `price_ode.initial_conditions`: the §2.5 IC expressions as strings (they reference other constants; the loader in 8.1 evaluates them).
  - `coupling`: the whole §7 `coupling:` block (spatial_to_ode, sig_1, recruitment with the asymmetric cd8, nearby_killing with the resist-DIRECT note, resistance, other_ode_to_spatial). Keep the `# discrepancy #N` comments.
  - Every value gets a trailing `# file:symbol:line` comment (copy from the dossier). **Beware YAML `#` in unquoted scalars** — quote any scalar containing `#` or `:` (the recurring parse bug in this package).

- [ ] **Step 4: Run tests, verify pass.** Also run `scripts/lint-workspace.py` — must stay OK.

- [ ] **Step 5: Commit** — `git add viva_cpm_studies/influenza/params.yaml tests/test_influenza_params.py docs/cc3d-reference/sego2022-global-ode.md && git commit -m "feat(influenza): cite Price-2015 global ODE + coupling params (Incr 8)"`

---

### Task 8.1: `GlobalODE` integrator core (10 systemic species)

**Files:**
- Create: `viva_cpm_studies/influenza/price_ode.py`
- Create: `tests/test_influenza_price_ode.py`

**Interfaces:**
- Consumes: `params.load_params()["price_ode"]` (Task 8.0).
- Produces:
  - `resolve_constants(price_ode: dict, *, num_epithelial: int, s_to_mcs: float = 60.0, cell_volume: float = 25.0) -> dict[str, float]` — resolves every `{base, scale}` into a concrete float using `s_t = s_to_mcs/86400`, `s_v = num_epithelial/250000`, `s_l = 1/250000/cell_volume`; resolves templated (`b_m,b_k,b_p`) from `homeostatic_pops_ode × s_v`; leaves `V0` to the caller.
  - `initial_state(consts: dict, *, num_epithelial: int, v0: float, resist0: float = 0.0) -> dict[str, float]` — the §2.5 ICs for the 10 integrated states + nearby surrogates (evaluate the IC expressions; `H,I,...` spatial vars are inputs, not part of the integrated state vector).
  - `class GlobalODE`: `__init__(self, consts, *, num_epithelial)`; `step(self, state: dict, inputs: dict, dt_seconds: float = 60.0) -> dict` — integrate the 10 integrated states forward `dt_seconds` (scipy `solve_ivp`, method `"LSODA"`), holding `inputs` (the spatialized species `H,I,M,K,E,DH,V,F,C,L,B_ei,G_ki` and derived `D,DI,Sigma1,Sigma2,R`) fixed over the step; return the new integrated-state dict. Only `{NB,N,T,X,A,B,P,W,G,O,M_nb,K_nb,E_nb}` are advanced.
  - `rhs(state_vec, consts, inputs) -> np.ndarray` — the d/dt for the integrated states, transcribed verbatim from dossier §2.3 (only the `[integrated]` lines: NB, N, T, X, A, B, P, W, G, O).

- [ ] **Step 1: Failing tests** — `tests/test_influenza_price_ode.py`:

```python
import math
import numpy as np
from viva_cpm_studies.influenza import params, price_ode

P = params.load_params()

def test_resolve_constants_scaling():
    c = price_ode.resolve_constants(P["price_ode"], num_epithelial=900)
    s_v = 900/250000
    assert c["a_11"] == 0.00061091213762049                      # scale None
    assert math.isclose(c["a_21"], 1.45266909729275 * s_v)       # scale s_v
    s_t = 60/86400
    assert math.isclose(c["g_ik"], 0.0000308183563969158 * s_t)  # scale s_t

def test_initial_state_healthy_homeostasis():
    c = price_ode.resolve_constants(P["price_ode"], num_epithelial=900)
    s0 = price_ode.initial_state(c, num_epithelial=900, v0=0.0)
    for k in ("NB","N","T","X","A","B","P","W","G","O","M_nb","K_nb","E_nb"):
        assert k in s0
    assert s0["T"] == 0.0 and s0["X"] == 0.0          # §2.5 ICs
    assert s0["A"] > 0 and s0["B"] > 0 and s0["P"] > 0  # homeostatic positives

def test_step_quiescent_is_near_stationary():
    # With no infection (I=0, V=0, no immune load) the healthy homeostasis
    # should barely move over one 1-min step.
    c = price_ode.resolve_constants(P["price_ode"], num_epithelial=900)
    s0 = price_ode.initial_state(c, num_epithelial=900, v0=0.0)
    ode = price_ode.GlobalODE(c, num_epithelial=900)
    inputs = dict(H=900, I=0, M=s0.get("M_nb",0), K=0, E=0, DH=0,
                  V=0.0, F=0.0, C=0.0, L=0.0, B_ei=0.0, G_ki=0.0)
    s1 = ode.step(s0, inputs, dt_seconds=60.0)
    for k in ("A","B","P","W","G","O"):
        assert math.isclose(s1[k], s0[k], rel_tol=1e-3, abs_tol=1e-6), k

def test_step_tnf_rises_with_infection():
    # TNF (T) is produced by b_t*M*Sigma2/(...); with macrophages present and
    # virus load driving Sigma2, T must increase from 0 over a step.
    c = price_ode.resolve_constants(P["price_ode"], num_epithelial=900)
    s0 = price_ode.initial_state(c, num_epithelial=900, v0=0.0)
    ode = price_ode.GlobalODE(c, num_epithelial=900)
    inputs = dict(H=800, I=100, M=50, K=0, E=0, DH=0,
                  V=1e4, F=0.0, C=1e3, L=0.0, B_ei=0.0, G_ki=0.0)
    s1 = ode.step(s0, inputs, dt_seconds=60.0)
    assert s1["T"] > s0["T"]        # TNF rises
    assert s1["X"] >= 0 and s1["A"] >= 0   # non-negative
```

- [ ] **Step 2: Run, verify failure** (module missing).

- [ ] **Step 3: Implement `price_ode.py`.**
  - `resolve_constants`: map each `scale` token to the product of `{s_t,s_v,s_l}` powers; `None`→1.0; `templated`→resolve `b_m,b_k,b_p` = `homeostatic_pops_ode[...] * s_v`.
  - Derived algebraic inputs each step (dossier §2.2): `D = num_epithelial − H − I`; `Sigma1 = a_11*T + a_12*D`; `Sigma2 = Sigma1 + a_21*V/(a_22+V)`; `DI = num_epithelial − H − I − DH`; `R` is supplied by the caller (from the spatial resistance) or computed `F/(a_rf/s_l*s_v + F)` when only `F` is given — expose a helper `derived_inputs(state, inputs, consts, num_epithelial)`.
  - `rhs`: transcribe the 10 `[integrated]` equations from §2.3 **exactly**. Nearby surrogates `M_nb,K_nb,E_nb` decay/accumulate — hold them constant in this process (they are maintained by recruitment in 8.4); their d/dt = 0 here.
  - `GlobalODE.step`: build the RHS closure with fixed `inputs`+`derived`, call `scipy.integrate.solve_ivp(fun, (0, dt), y0, method="LSODA", rtol=1e-6, atol=1e-9)`, take the final column, clamp tiny negatives to 0, return dict. On solver failure raise a clear error (the source stops the sim on ODE exception).
  - Hill helper `hill(v,a,h) = 0.0 if v<=0 else 1/(1+(a/v)**h)` (dossier appendix).

- [ ] **Step 4: Run tests, verify pass.**
- [ ] **Step 5: Commit** — `feat(influenza): Price-2015 global ODE integrator core (Incr 8 Task 8.1)`

---

### Task 8.2: spatial→ODE coupling + `run_global_coupling` driver

**Files:**
- Modify: `viva_cpm_studies/influenza/run.py`
- Modify: `tests/test_influenza_price_ode.py`

**Interfaces:**
- Consumes: `price_ode.GlobalODE` (8.1); the existing Increment-6/7 scenario builders in `immune.py`/`run.py` (macrophages + fields + optional NK/CD8).
- Produces: `run.run_global_coupling(*, num_epithelial=None, side=..., steps=..., seed=..., with_immune=True, **kw) -> dict` returning at least `{"mcs": [...], "ode": {sp: [...] for sp in integrated_states}, "spatial": {"H":[...],"I":[...],"M":[...],"K":[...],"E":[...],"DH":[...],"V":[...],"F":[...],"C":[...],"L":[...]}, "sigma1": [...]}`. Each MCS: (1) read spatial aggregates → build `inputs`; (2) `ode.step`; (3) record. (Back-couplings sig_1/recruitment/nearby-killing are wired in 8.3–8.5; this task establishes the loop and the spatial→ODE direction.)

- [ ] **Step 1: Failing test** — append to `tests/test_influenza_price_ode.py`:

```python
def test_run_global_coupling_smoke_and_shapes():
    from viva_cpm_studies.influenza import run
    r = run.run_global_coupling(side=30, steps=8, seed=1, with_immune=True)
    n = len(r["mcs"])
    assert n == 8
    for sp in ("T","X","A","P","G"):
        assert len(r["ode"][sp]) == n
    for sp in ("H","I","V","C","L"):
        assert len(r["spatial"][sp]) == n
    assert len(r["sigma1"]) == n
    # spatial aggregates are actually fed in: infected count is a non-trivial series
    assert all(v >= 0 for v in r["spatial"]["I"])
    # ODE advanced (TNF not stuck at exactly 0 once infection present) — allow 0 if no infection seeded
    assert all(t >= 0 for t in r["ode"]["T"])
```

- [ ] **Step 2: Run, verify failure.**
- [ ] **Step 3: Implement `run_global_coupling`.** Reuse the Increment-6/7 scenario (a small infected patch + macrophages secreting chemokine/IL-10; optionally NK/CD8). Compute `num_epithelial` from the sheet (or the `side` arg). Build `GlobalODE` with resolved constants. Each MCS:
  - advance the CPM world one step (existing `world.run(1)` / step call used by sibling drivers — match `run_macrophage_signaling`);
  - read spatial aggregates via the existing helpers: counts by type (`world.cells_by_type` or the count util the drivers already use), field integrals `sum(world.field_conc(fi))/dim.z` for V/F/C/L (match the per-pixel/z convention used in Increments 2/6), dead-uninfected count for `DH`, `B_ei = g_ei_sum`, `G_ki = g_ki_sum` from Σ resist over infected cells (use `params.coupling` constants);
  - `inputs = {...}`; `state = ode.step(state, inputs)`; append `state[sp]` and `spatial[...]` and `Sigma1`.
  - Keep the default fast (`side` small, `steps` small); document that paper-scale is Increment 9.

- [ ] **Step 4: Run tests, verify pass.**
- [ ] **Step 5: Commit** — `feat(influenza): spatial->ODE coupling + run_global_coupling driver (Incr 8 Task 8.2)`

---

### Task 8.3: resolve the `sig_1` stub with dynamic ODE TNF

**Files:**
- Modify: `viva_cpm_studies/influenza/run.py` (in `run_global_coupling` / the secretion wiring)
- Modify: `viva_cpm_studies/influenza/signaling.py` (docstring only — the stub note becomes "resolved in Incr 8")
- Modify: `tests/test_influenza_price_ode.py`

**Interfaces:**
- Consumes: `run_global_coupling` (8.2) now has live `Sigma1 = a_11*T + a_12*D` each MCS.
- Produces: in the driver, macrophage chemokine/IL-10 secretion scale is computed with the **dynamic** `sig_1` (from the ODE `T` + spatial dead count `D`) passed into the existing `signaling.macrophage_secretion_scale(il10_local, sig_1, g_1, g_2, d_2)` — replacing `params.il10.sig_1_stub`. Record `sigma1` and confirm it becomes dynamic (rises as infection/TNF build). Do NOT alter the Michaelis functional form (constraint / discrepancy #11).

- [ ] **Step 1: Failing test** — append:

```python
def test_sig1_is_dynamic_not_stubbed():
    from viva_cpm_studies.influenza import run, params
    P = params.load_params()
    stub = P["il10"]["sig_1_stub"]
    r = run.run_global_coupling(side=30, steps=20, seed=3, with_immune=True, seed_infection_frac=0.05)
    s1 = r["sigma1"]
    # sig_1 is no longer the constant stub: it varies across the run
    assert max(s1) != min(s1), "sig_1 must be dynamic (ODE TNF + spatial dead count)"
    assert any(abs(v - stub) > 1e-9 for v in s1)
    # and it should rise from ~baseline as TNF (T) and dead count build with infection
    assert s1[-1] >= s1[0]
```

- [ ] **Step 2: Run, verify failure** (sigma1 constant/absent).
- [ ] **Step 3: Implement.** In `run_global_coupling`, compute `sig_1 = a_11*T + a_12*D` from the current ODE `T` and spatial `D = num_epithelial - H - I` each MCS, and use it when the driver sets the macrophage secretion scale (via `world.set_cell_secretion_scale` for chemokine + IL-10, exactly as `run_macrophage_signaling` did, but with the dynamic `sig_1` instead of the stub constant). Record `r["sigma1"]`. Update `signaling.py`'s module docstring: the `sig_1 STUB` paragraph now reads that Increment 8 resolves it via the ODE (keep the historical note, mark RESOLVED, cite dossier §4b(i)). Leave `macrophage_secretion_scale` code unchanged.

- [ ] **Step 4: Run tests, verify pass** (plus the 8.1/8.2 tests still green).
- [ ] **Step 5: Commit** — `feat(influenza): resolve sig_1 stub with dynamic ODE TNF (Incr 8 Task 8.3)`

---

### Task 8.4: immune-cell recruitment (Hill inflows, asymmetric)

**Files:**
- Create: `viva_cpm_studies/influenza/recruitment.py`
- Create: `tests/test_influenza_recruitment.py`
- Modify: `viva_cpm_studies/influenza/run.py` (wire recruitment into `run_global_coupling`)

**Interfaces:**
- Consumes: `params.coupling.recruitment` (8.0); the ODE fields `C` (chemokine), `P` (APC) and outflow terms `G_ki`,`B_ei` (8.2).
- Produces:
  - `inflow_rate(cell_type, driver_field, consts, *, b_type, mu_type) -> float` and per-type wrappers `macrophage_inflow(C, consts)`, `nk_inflow(C, consts)`, `cd8_inflow(P, consts)`; `macrophage_outflow(consts)`, `nk_outflow(G_ki, consts)`, `cd8_outflow(B_ei, consts)` — the exact dossier §4b(ii) forms with `hill(v,a,h)=1/(1+(a/v)**h)`. CD8 inflow has NO homeostatic baseline (discrepancy #12).
  - In `run_global_coupling`: each MCS, convert inflow rates to expected new-cell counts (`ul_rate_to_prob` / Poisson draw, seeded RNG), seed new immune cells via the existing `immune.py` seeding helper (interior/near-lesion per the source `local_ratio`), and remove cells per outflow. Maintain the nearby surrogates `M_nb,K_nb,E_nb` from the `1 - local_ratio` fraction.

- [ ] **Step 1: Failing tests** — `tests/test_influenza_recruitment.py`:

```python
import math
from viva_cpm_studies.influenza import params, price_ode, recruitment
P = params.load_params()
C = price_ode.resolve_constants(P["price_ode"], num_epithelial=900)

def test_macrophage_inflow_monotone_in_chemokine():
    lo = recruitment.macrophage_inflow(10.0, C)
    hi = recruitment.macrophage_inflow(1e5, C)
    assert hi > lo                                  # more chemokine -> more inflow
    assert recruitment.macrophage_inflow(0.0, C) == C["mu_m"]*C["b_m"]  # hill(0)=0 -> baseline only

def test_cd8_inflow_apc_driven_no_baseline():
    # CD8 inflow is P-driven and has NO homeostatic baseline (discrepancy #12):
    assert recruitment.cd8_inflow(0.0, C) == 0.0     # hill(0)=0, no +mu*b term
    assert recruitment.cd8_inflow(1e5, C) > 0.0

def test_nk_inflow_has_baseline():
    assert recruitment.nk_inflow(0.0, C) == C["mu_k"]*C["b_k"]

def test_outflow_terms_carry_resistance_load():
    # NK/CD8 outflow rise with the resistance-weighted infected load (G_ki/B_ei)
    assert recruitment.nk_outflow(0.5, C) > recruitment.nk_outflow(0.0, C)
    assert recruitment.cd8_outflow(0.5, C) > recruitment.cd8_outflow(0.0, C)
    assert math.isclose(recruitment.macrophage_outflow(C), C["mu_m"])
```

- [ ] **Step 2: Run, verify failure.**
- [ ] **Step 3: Implement `recruitment.py`** (pure functions, dossier §4b(ii) forms). Then wire into `run_global_coupling`: seeded RNG; per MCS compute inflow counts, seed via `immune.py`'s existing new-cell helper (reuse the Increment-5/7 seeding; place per `local_ratios` {macro 1.0, nk 0.75, cd8 0.75}); apply outflow removals; update `M_nb/K_nb/E_nb`. Add an INTEGRATION assertion to `test_influenza_recruitment.py` (keep FAST — small side, ~15 steps): with chemokine/APC present, macrophage (and NK) counts do NOT collapse to zero and respond upward vs a no-signal control. Deterministic/seeded.

- [ ] **Step 4: Run tests, verify pass.**
- [ ] **Step 5: Commit** — `feat(influenza): ODE-driven immune recruitment (Hill inflows) (Incr 8 Task 8.4)`

---

### Task 8.5: nearby-population contact-killing term

**Files:**
- Modify: `viva_cpm_studies/influenza/killing.py`
- Modify: `viva_cpm_studies/influenza/run.py` (apply nearby term in `run_global_coupling`)
- Modify: `tests/test_influenza_killing.py` (add nearby-term tests; keep FAST)

**Interfaces:**
- Consumes: `params.coupling.nearby_killing` (8.0); `eta`, the nearby surrogates `K_nb`,`E_nb` (8.4).
- Produces: `nearby_kill_rate(g_i: float, eta: float, num_nearby: float, cell_resist: float) -> float` = `(g_i/eta) * num_nearby * cell_resist` (≥0; 0 when `num_nearby<=0` or `eta<=0`). **Multiplies `cell_resist` DIRECTLY (discrepancy #7).** In `run_global_coupling`, each infected cell's total death probability from cytotoxic cells combines the Incr-7 LOCAL contact term (`contact_kill_rate` via `cell_contact_area_by_type`) AND this well-mixed nearby term, via `ul_rate_to_prob(local + nearby) = 1 - exp(-(local+nearby))`.

- [ ] **Step 1: Failing tests** — add to `tests/test_influenza_killing.py`:

```python
from viva_cpm_studies.influenza import killing

def test_nearby_kill_rate_form():
    r = killing.nearby_kill_rate(g_i=1e-3, eta=0.0049, num_nearby=4.0, cell_resist=0.5)
    assert math.isclose(r, (1e-3/0.0049)*4.0*0.5)
    assert killing.nearby_kill_rate(1e-3, 0.0049, 0.0, 0.5) == 0.0   # no nearby -> 0
    assert killing.nearby_kill_rate(1e-3, 0.0, 4.0, 0.5) == 0.0      # eta<=0 -> 0

def test_nearby_kill_rate_resist_direct_discrepancy7():
    # DIRECT in resist (NOT 1-resist): rate rises with resist
    assert killing.nearby_kill_rate(1e-3,0.0049,4.0,0.9) > killing.nearby_kill_rate(1e-3,0.0049,4.0,0.1)
    assert killing.nearby_kill_rate(1e-3,0.0049,4.0,0.0) == 0.0
```

- [ ] **Step 2: Run, verify failure.**
- [ ] **Step 3: Implement `nearby_kill_rate`** in `killing.py` (flag discrepancy #7 in the docstring, mirroring `contact_kill_rate`). Wire into `run_global_coupling`: for each infected cell, `total_rate = contact_kill_rate(...) + nearby_kill_rate(g_ik, eta, K_nb, resist) [NK] + …[CD8]`; draw death with `1 - exp(-total_rate)`. Add a FAST integration assertion: with a nearby NK/CD8 surrogate population present, infected clearance is stronger than with `K_nb=E_nb=0` (same seed) — the nearby term does work. Do NOT add a multi-minute test.

- [ ] **Step 4: Run tests, verify pass.**
- [ ] **Step 5: Commit** — `feat(influenza): nearby-population contact-killing term (Incr 8 Task 8.5)`

---

### Task 8.6: `global-coupling` study + viz + composite + membership

**Files:**
- Modify: `viva_cpm_studies/influenza/viz.py` (add `global_coupling_figure`)
- Modify: `viva_cpm_studies/composites/influenza.py` (add `global_coupling` composite — auto-registers; do NOT edit `__init__.py`)
- Create: `workspace/studies/global-coupling/study.yaml`
- Modify: `workspace/investigations/influenza-sego2022/investigation.yaml` (add member)
- Modify: `tests/test_influenza_viz.py` (add a fast non-vacuous test)

**Interfaces:**
- Consumes: `run.run_global_coupling` result keys (8.2–8.5). Model the study on `workspace/studies/cytotoxic-killing/study.yaml` (schema_version 3).

- [ ] **Step 1: Failing viz test** — add to `tests/test_influenza_viz.py`: build a small `run_global_coupling(side=30, steps=8, seed=0)` result, assert `global_coupling_figure(result)` returns a `matplotlib.figure.Figure` with ≥2 axes (one systemic-species panel: T/X/A/P over MCS; one spatial panel: I/M/K/E + sigma1). Keep FAST.

- [ ] **Step 2: Run, verify failure.**
- [ ] **Step 3: Implement** `global_coupling_figure` (Agg; systemic-species trajectories + spatial aggregates + dynamic sig_1). Add a `global_coupling` live-demo composite in `composites/influenza.py` mirroring the `cytotoxic_response` precedent (reuse the field-spec helpers). Hand-author `workspace/studies/global-coupling/study.yaml`:
  - `schema_version: 3`, `investigation: influenza-sego2022`, phase `Simulate`, `status: in-progress`, `verdict: documented`, `biological_validation: PENDING`.
  - Report HONESTLY: the hybrid Price-2015 global ODE (10 systemic species) integrates once per MCS coupled to the spatial model; **the three stubs are RESOLVED** — sig_1 dynamic (ODE TNF + spatial dead count, saturating IL-10-inhibited form, discrepancy #11), recruitment engaged (asymmetric Hill inflows, discrepancy #12), nearby-killing term added (discrepancy #7 resist-direct). Systemic species (TNF/ROS/antibody/APC/IL-12/IFN-γ/CD4⁺/B/neutrophils) evolve consistently with infection.
  - **KEY CAVEATS (central, honest):** this is NOT the quantitative Fig-3B reproduction (Increment 9) — validation is qualitative consistency at reduced scale; discrepancies #9 (hybrid), #10 (spatial-string scalings), #11, #12 recorded; scipy-LSODA approximates RoadRunner-CVODE; any species with an unphysical trajectory or a calibration gap is flagged for Incr 9. Do NOT claim quantitative match. Target: Fig 3B global panels (qualitative).
  - Add the study to `investigation.yaml` members.

- [ ] **Step 4: Run** `tests/test_influenza_viz.py` + `scripts/lint-workspace.py` (OK; the new study may be `dashboard ✗` PENDING — expected). Verify.
- [ ] **Step 5: Commit** — `feat(influenza): global-coupling study + viz + composite (Incr 8 Task 8.6)`

---

## Self-Review notes (author)

- **Spec coverage:** Increment 8 row = "Full Price-2015 global ODE, η/θ scaling, bidirectional coupling; study `global-coupling`; Target Fig 3B global panels" → Tasks 8.1 (ODE), 8.0/8.1 (scaling), 8.2–8.5 (coupling + 3 stubs), 8.6 (study/target). Covered, with the source-mandated **hybrid** correction (discrepancy #9) made explicit.
- **Discrepancy #9 is load-bearing:** every task treats the 10 integrated states as the ODE and the spatial species as inputs. A reviewer must reject any task that integrates `{H,I,M,E,K,L,C,F,V,DH}`.
- **Honesty:** reduced-scale qualitative consistency only; PENDING quantitative repro (Incr 9). All four new discrepancies (#9–#12) surface in params/code/study.
- **Type consistency:** `resolve_constants` → `GlobalODE`/`recruitment`/`killing` all consume the same resolved-float `consts` dict; `run_global_coupling` result keys are fixed in 8.2 and reused by 8.6 viz.
