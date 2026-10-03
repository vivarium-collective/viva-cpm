# Merks ECM Reciprocity 2D — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Reproduce, in 2D with the real published static-adhesion TST-MD engine, the three mechanical-reciprocity results (Figs 2–4) of Keijzer & Merks 2026 (3D, unreleased), as a viva-cpm investigation consuming a new Dockerized `viva-tstmd` bridge.

**Architecture:** A sibling `viva-tstmd/` repo wraps the real TST-MD engine (TST C++ + HOOMD + MUSCLE3), built once in a Linux CPU-only Docker image and driven from Python via a file handshake (the `viva-chaste` pattern). The viva-cpm investigation `merks-ecm-reciprocity-2d` consumes that bridge: three studies sweep crosslinking / fiber stiffness and measure equilibrium area, percolation, fiber reorientation, and densification, rendered as baked, interactive visualizations.

**Tech Stack:** process-bigraph, Docker (Colima on arm64), TST-MD (Zenodo 7906973), HOOMD-blue, MUSCLE3, Python 3.12 (uv/.venv), pytest, plotly, `viva_superpowers.visualization`.

**Spec:** `docs/superpowers/specs/2026-10-03-merks-ecm-reciprocity-2d-design.md`

## Global Constraints

- Python `>=3.12`; `process-bigraph>=1.8.5`, `bigraph-schema>=1.4.4` (match viva-cpm pyproject).
- Real engine only — never a mock/stub/reimplementation of the dynamics (viva-expert default). Faked dynamics are an explicit, separate opt-in the user has NOT given.
- Fidelity is **directional/2D**, not quantitative/3D: acceptance bands encode monotonicity / threshold-existence / biphasic-shape, never absolute 3D numbers.
- Engine to bridge = **static-adhesion** TST-MD (ref [38]); NOT the dynamic-FA model (ref [39]).
- Bridge package name `pbg_tstmd` (matches `pbg_chaste`, `pbg_compucell3d`).
- No AI attribution in commits/PRs. Investigation branch is a living integration branch (draft PR, `investigation:` title prefix), not a merge target. Never auto-merge.
- viva-cpm commits only in worktree `~/code/viva-cpm--merks-ecm-reciprocity-2d`; run tests with that worktree on `PYTHONPATH`.
- CPU-tractable sim sizes (smaller lattice / fewer seeds than the paper's HPC sweeps); document the chosen sizes per study.

## Review Focus

- **Docker daemon down / image missing** at run time → bridge must raise a clear, actionable error (preflight guard), not hang or produce empty output. (Task 1.2)
- **A single coupled run fails / times out mid-sweep** → the sweep records the failure for that parameter point and continues, rather than aborting the whole study or silently dropping the point. (Task 1.3)
- **Degenerate ECM network** (zero crosslinks → disconnected graph) fed to the percolation metric → giant-component fraction is well-defined (largest component / n), not a divide-by-zero or crash. (Task 2.1)
- **Reorientation metric on an empty/near-cell fiber set** (distance bin with no fibers) → `q(r)/q0(r)` returns NaN/None for that bin, not a spurious 0 or exception. (Task 2.2)
- **Densification when `V_far` is empty or zero-area** → `ρ_dens` guards the ratio (returns NaN, flagged), not inf. (Task 2.3)

---

## Phase 0 — M0: engine up in Docker (GO/NO-GO gate) — via `/viva-expert`

> Not TDD-able in advance: this is a real-tool build + I/O discovery. Execute with the **viva-superpowers:viva-expert** skill (heavy mode, default = bridge the REAL tool), creating the sibling repo `~/code/viva-tstmd/`. The deliverables below are the gate; everything downstream consumes the facts discovered here.

### Task 0.1: Scaffold `viva-tstmd` + build the real engine in Docker

**Files:**
- Create: `~/code/viva-tstmd/Dockerfile`, repo scaffold (pyproject, pbg_tstmd/, tests/, README, NEXT_STEPS) per viva-expert heavy mode.

- [ ] Invoke `/viva-expert` to wrap TST-MD (static-adhesion) as a process-bigraph Process, sibling repo `viva-tstmd`.
- [ ] Author `Dockerfile`: Linux CPU-only build from `rmerks/Tissue-Simulation-Toolkit` branch `TST2.0` + `make with_adhesions` (no GPU/MPI), bundling HOOMD + MUSCLE3; source cross-checked against Zenodo 7906973.
- [ ] `colima start` if daemon down; `docker build` the image. Expect a long HOOMD compile.
- [ ] **GATE:** image builds successfully. If unbuildable after reasonable effort → STOP, report, consider fallback (native build on mini / scope cut). Do not proceed to studies.

### Task 0.2: One reference coupled run + characterize I/O

**Deliverable (facts the rest of the plan needs):** a committed `docs/engine-io.md` in `viva-tstmd` recording: the YMMSL settings keys for crosslink number, fiber stiffness `k`, cell target volume, adhesion count, MCS, output interval; the on-disk layout + Python structure of a pickle state dump (cell voxel/spin field, fiber bead coordinates, bond/crosslink lists); and how to extract cell area per output step.

- [ ] Run `muscle_manager --start-all ymmsl/adhesions.ymmsl ...` for a short contractile-cell run inside the container.
- [ ] Load a resulting pickle dump in Python; record its structure in `docs/engine-io.md`.
- [ ] Extract cell-area-over-time from the dumps; confirm the contractile cell contracts. **GATE passes** when area-over-time is extracted from a real run.

---

## Phase 1 — M1: the Python bridge (`viva-tstmd/pbg_tstmd/`)

### Task 1.1: `analysis.py` metrics — see Phase 2 (pure, synthetic-data TDD; build these first, they gate nothing on Docker)

### Task 1.2: `runtime.py` — `TstmdSession` + preflight guards

**Files:**
- Create: `~/code/viva-tstmd/pbg_tstmd/runtime.py`, `~/code/viva-tstmd/tests/test_runtime.py`
**Interfaces:**
- Produces: `docker_available() -> bool`, `image_present(tag:str) -> bool`, `class TstmdSession(params: dict)` with `.run() -> Path` (run dir of pickle dumps), `.close()`, context-manager support.
- Consumes: YMMSL keys + pickle layout from `docs/engine-io.md` (Task 0.2).

- [ ] **Step 1 (test, no Docker):** `test_run_without_docker_raises` — monkeypatch `docker_available` → False; `TstmdSession(...).run()` raises `RuntimeError` naming "docker"/"colima".
- [ ] **Step 2:** run → FAIL.
- [ ] **Step 3:** implement `docker_available()` (shell `docker info`), `image_present()`, and `TstmdSession` whose `run()` preflights then `docker run`s, writes the YMMSL settings file from `params`, invokes `muscle_manager`, awaits the run dir. (I/O exactly per `docs/engine-io.md`.)
- [ ] **Step 4:** run → PASS.
- [ ] **Step 5 (Docker-gated test):** `test_short_run_contracts` (mark `@pytest.mark.docker`, skip if `not docker_available()`) — a short run returns a run dir whose extracted area time-series ends below its start.
- [ ] **Step 6:** run with Docker up → PASS.
- [ ] **Step 7:** commit.

### Task 1.3: `processes.py` — `TstmdEcmProcess` + resilient sweep helper

**Files:**
- Create: `~/code/viva-tstmd/pbg_tstmd/processes.py`, `~/code/viva-tstmd/tests/test_processes.py`
**Interfaces:**
- Consumes: `TstmdSession` (1.2), `analysis.py` (2.x).
- Produces: `class TstmdEcmProcess(Process)` — `config_schema` keys `{lattice:int, n_cross:int, fiber_stiffness:float, target_volume:int, seed:int, mcs:int, output_interval:int}`; `inputs()/outputs()` expose `{time, cell_area, cell_voxels, fiber_beads, crosslinks}`; `update(state, interval)`. Plus `run_sweep(param:str, values:list, base:dict) -> list[dict]` that catches a failed point, records `{value, error}`, and continues.

- [ ] **Step 1 (test):** `test_process_registers_and_schema` — `core=allocate_core(); core.register_process('TstmdEcmProcess', TstmdEcmProcess)`; assert `config_schema`/`outputs()` contain the keys above. (No Docker.)
- [ ] **Step 2:** run → FAIL.
- [ ] **Step 3:** implement `TstmdEcmProcess` delegating to `TstmdSession` + `analysis.py`.
- [ ] **Step 4:** run → PASS.
- [ ] **Step 5 (test, no Docker):** `test_run_sweep_continues_on_failure` — inject a runner stub that raises on one value; assert the returned list has an `error` entry for it and results for the others.
- [ ] **Step 6:** implement `run_sweep`; run → PASS.
- [ ] **Step 7 (Docker-gated):** `test_process_update_contracts` — one `update({}, interval=…)` returns decreasing `cell_area`.
- [ ] **Step 8:** commit.

---

## Phase 2 — `analysis.py` metrics (pure functions, full TDD, no Docker)

**Files:** Create `~/code/viva-tstmd/pbg_tstmd/analysis.py`, `~/code/viva-tstmd/tests/test_analysis.py`.
**Interfaces — Produces:**
- `giant_component_fraction(nodes:list, edges:list[tuple]) -> float`
- `reorientation_q(fiber_centers, fiber_dirs, cell_centroid, lo=0.9) -> float` and `q_ratio_vs_distance(..., q0_by_bin, bins) -> dict[bin,float]`
- `densification_factor(bead_coords, cell_mask, close_px, far_margin_px) -> float`

### Task 2.1: giant component / percolation (Fig 2B)
- [ ] **Step 1:** `test_giant_component` — two triangles joined by one edge over 6 nodes → fraction `1.0`; fully disconnected 5 nodes → `1/5`; **zero edges** → `1/n` (no crash). 
- [ ] **Step 2:** run → FAIL.
- [ ] **Step 3:** implement via union-find; guard empty-edge case.
- [ ] **Step 4:** run → PASS. **Step 5:** commit.

### Task 2.2: reorientation q/q0 (Fig 4A–C, eq. 7)
- [ ] **Step 1:** `test_reorientation` — all fibers pointing at the centroid (γ≈1) → `q≈1.0`; uniformly random directions → `q≈0.1` (within tolerance); **a distance bin with no fibers** → `q_ratio_vs_distance` yields `None`/NaN for that bin, not 0.
- [ ] **Step 2:** FAIL. **Step 3:** implement `γ_i = |c×v|/|c·v|`-style alignment per eq. 7, `q=Prob(0.9≤γ≤1)`, ratio helper with empty-bin guard. **Step 4:** PASS. **Step 5:** commit.

### Task 2.3: densification factor (Fig 4D–E, eqs. 8–9)
- [ ] **Step 1:** `test_densification` — beads packed within `close_px` of the cell, sparse far away → `ρ_dens > 1`; uniform field → `ρ_dens ≈ 1`; **empty far region** → NaN (flagged), not inf.
- [ ] **Step 2:** FAIL. **Step 3:** implement subnode interpolation (eq. 8) + ρ_close/ρ_far ratio (eq. 9) in 2D (areas); guard zero far-area. **Step 4:** PASS. **Step 5:** commit.

---

## Phase 3 — Studies (viva-cpm worktree). Each study is independent once M1 is green.

For each study: create `workspace/studies/<slug>/study.yaml` (schema v3), `targets/<fig>.json` (directional oracle), a pytest module under `tests/` whose nodeid is the study's `behavior_tests[].evaluated_by`, bake the real sweep output into `viva_cpm_studies/visualizations/_tstmd_data.py`, and add it to `members[]` of the investigation. Study sizes (lattice, seeds, sweep points) chosen CPU-tractable and recorded in `study.yaml:description`.

### Task 3.1: `tstmd-contraction-vs-crosslinking` (Fig 2)
- [ ] Run the real `n_cross` sweep via `run_sweep`; record equilibrium area + giant-component per point into `_tstmd_data.py`.
- [ ] Write `targets/fig2.json`: assertions `equilibrium_area` monotonic-nondecreasing in `n_cross`; `percolation` shows a threshold (a step up); `n_cross=0` area ≈ full contraction.
- [ ] Write `tests/test_tstmd_fig2.py` evaluating the baked data against `targets/fig2.json`; run → PASS on real data (or record divergence honestly).
- [ ] Write `study.yaml` (report/verdict/findings wired to the test nodeid). Commit.

### Task 3.2: `tstmd-fiber-stiffness-drag` (Fig 3)
- [ ] Real `k` sweep, with and without crosslinks; record equilibrium area + area-vs-time curves.
- [ ] `targets/fig3.json`: equilibrium area monotonic in `k`; non-crosslinked high-`k` retains nonzero area (drag); intermediate `k` reaches ≈0 but slower (time-to-threshold ordering).
- [ ] `tests/test_tstmd_fig3.py`; `study.yaml`. Commit.

### Task 3.3: `tstmd-ecm-remodeling` (Fig 4)
- [ ] Real sweep over `n_cross`; record `q/q0` vs time, `q(r)/q0(r)` vs distance, `ρ_dens` vs time and vs `n_cross`.
- [ ] `targets/fig4.json`: `q/q0 > 1` and increasing with `n_cross`; `ρ_dens` vs `n_cross` is **biphasic** (argmax at an interior crosslink value).
- [ ] `tests/test_tstmd_fig4.py`; `study.yaml`. Commit.

---

## Phase 4 — Visualizations (viva-cpm worktree)

**Files:** Create `viva_cpm_studies/visualizations/tstmd_studies.py` + `_tstmd_data.py` (baked real output).

### Task 4.1: the figure functions
- [ ] One `@as_visualization` function per study (`TstmdContractionCrosslinking`, `TstmdFiberStiffnessDrag`, `TstmdEcmRemodeling`) returning `{"html": ...}` from baked data: interactive plotly (sweep curves + threshold + biphasic + q/q0), with cell+fiber-network snapshots as insets. Demo kwargs provided.
- [ ] Verify auto-discovery (render each via the workspace render path); wire `study.yaml:visualizations[].address = local:<Name>` and `embed_visualizations[].url`.
- [ ] Commit.

---

## Phase 5 — Resources + investigation + report

### Task 5.1: Register the 3 papers in Resources
- [ ] Copy PDFs → `workspace/references/papers/`: `keijzer-merks-2026-3d-cpm-ecm.pdf` (from `~/Downloads/2609.02375v1.pdf`), `keijzer-2025-focal-adhesions.pdf` (from `~/Downloads/fcell-12-1462277.pdf`); fetch bioRxiv PDF → `tsingos-2023-cpm-beadspring.pdf`.
- [ ] Add `@article{}` entries to `workspace/references/papers.bib` (keys `keijzer2026cpmecm3d`, `keijzer2025focaladhesions`, `tsingos2023beadspring`) with DOIs.
- [ ] Register each in `workspace.yaml:references_pdfs` as `{bib_key, path, sha256}` (compute sha256).
- [ ] Add `workspace/references/notes/<key>.md` one-paragraph summaries. Commit.

### Task 5.2: Investigation YAML + report
- [ ] Create `workspace/investigations/merks-ecm-reciprocity-2d/investigation.yaml` (schema v2): question, hypothesis, `executive` (verdict + fidelity caveat), `scientific_argument` (key_figures → the 3 studies), `members` = the 3 study slugs, `acceptance_criteria`, `at_a_glance`.
- [ ] Run `/viva-report` (reviewer-readiness audit → lint → render); fix flags; confirm `report.html` renders the 3 figures.
- [ ] Open draft investigation PR (`investigation:` prefix). Commit.

---

## Self-Review

- **Spec coverage:** §1 framing → Phase 3 studies + investigation.yaml (5.2). §3 bridge → Phases 0–2. §4 studies/bands → Phase 3. §5 viz → Phase 4. §6 Resources → 5.1. §7 milestones → phase order (M0 gate = Phase 0). §8 out-of-scope → honored (no 3D, no dynamic FA, no Fig 1). §9 process → Global Constraints. No gaps.
- **Placeholders:** engine-dependent bridge code (runtime/process bodies) intentionally defers exact I/O to the `docs/engine-io.md` deliverable from Task 0.2 — this is a real discovery dependency, not a hand-wave; all pure-function tasks (Phase 2) and YAML/oracle tasks (Phase 3/5) are concrete.
- **Type consistency:** `giant_component_fraction`, `reorientation_q`/`q_ratio_vs_distance`, `densification_factor`, `TstmdSession.run`, `TstmdEcmProcess`, `run_sweep` referenced consistently across Phases 1–3.
- **Review Focus:** 5 failure modes above each pinned to a task's tests (1.2, 1.3, 2.1, 2.2, 2.3).
