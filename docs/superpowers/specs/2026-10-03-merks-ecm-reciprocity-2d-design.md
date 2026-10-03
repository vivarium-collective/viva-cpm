# Design: `merks-ecm-reciprocity-2d` investigation + `viva-tstmd` bridge

**Date:** 2026-10-03
**Slug:** `merks-ecm-reciprocity-2d`
**Status:** approved design → spec (pre-implementation)
**Author:** Eran Agmon (+ Claude)

## 1. Purpose & research framing

Reproduce, in **2D with the real published engine**, the three core
mechanical-reciprocity results of the (unreleased, 3D) paper:

> Keijzer, K.A.E. & Merks, R.M.H. (2026). *A 3D hybrid cellular Potts model
> with a discrete deformable fiber network: modeling cell contraction and
> extracellular matrix remodeling.* arXiv:2609.02375v1 (in revision, BMMB).

The 3D paper has **no public code** and its results come from GPU/HPC
stochastic simulations. Its model is a 3D extension of a **published, open 2D
engine** — ref [38], Tsingos et al. 2023, *Biophys. J.* — whose code is
archived on Zenodo (record 7906973: `TST-MD` + `PyTST`). The 3D paper uses
**static adhesion sites**; it names *dynamic* focal adhesions (ref [39],
Keijzer et al. 2025, *Front. Cell Dev. Biol.*) as a limitation / future step.
Therefore the faithful engine to bridge is the **static-adhesion TST-MD**.

**Research question:** Does the published 2D hybrid CPM–bead-spring ECM
engine — the real ancestor of Merks' 3D model — reproduce the same
mechanical-reciprocity phenomena the 3D paper reports?

**Fidelity contract (stated up front in the report):** 2D, real-engine,
*directional* reproduction. Because the sims are stochastic and lower
dimensional than the 3D target, acceptance is **trend/shape-level**
(monotonicity, threshold existence, biphasic peak), **not** quantitative or
pixel identity with the 3D figures. The user chose "adapt the real 2D
upstream code" over a (GPU/HPC-only, unreachable) quantitative 3D redo.

### Figures reproduced (2D analogs)

- **Fig 2** — ECM crosslinking limits cell contraction: equilibrium cell
  size vs crosslinking, plus network percolation (giant component);
  high- vs low-crosslinking snapshots.
- **Fig 3** — fiber stiffness introduces an effective drag on contraction:
  equilibrium size vs fiber stiffness `k`; size-vs-time curves with and
  without crosslinkers.
- **Fig 4** — contractile cells remodel the ECM: fiber reorientation
  `q/q0` over time and vs distance `r`; biphasic densification factor
  `ρ_dens` vs time and vs crosslink count.

Out: Fig 1 (schematic — we may build our *own* model-overview viz instead).

## 2. Architecture — two deliverables

1. **`viva-tstmd/`** — a **new sibling repo** (`~/code/viva-tstmd/`): the
   reusable, Dockerized real-engine bridge, built via `/viva-expert` heavy
   mode, installed into the viva-cpm workspace catalog.
2. **viva-cpm investigation** `merks-ecm-reciprocity-2d` — consumes the
   bridge; holds the study YAMLs, baked real-run visualizations, the report,
   and the 3 papers in Resources. Branch:
   `investigation/merks-ecm-reciprocity-2d` (worktree
   `~/code/viva-cpm--merks-ecm-reciprocity-2d`).

Rationale: matches AGENTS.md "feature-first" (reusable infra = a module;
the investigation *consumes* it) and the sibling-repo precedent for real
external engines (`viva-chaste`, `viva-compucell3d`). Reusable beyond this
paper; isolates the heavy native toolchain from the studies package.

## 3. `viva-tstmd/` bridge repo

Mirrors `viva-chaste` (real C++ engine via Docker + file handshake).

### 3.1 Build — Docker, Linux, CPU-only

- `Dockerfile`: Linux CPU-only build of **static-adhesion TST-MD** from the
  Zenodo deposit / `rmerks/Tissue-Simulation-Toolkit` branch `TST2.0`,
  `make with_adhesions` (no GPU, no MPI — or OpenMPI for >2 cores), bundling
  HOOMD-blue + MUSCLE3. Built once; runs identically on laptop + mini.
- Rationale: native build of HOOMD+MUSCLE3 on Apple Silicon is the project's
  main risk; Docker (already used by `viva-chaste`) sidesteps it. Docker +
  Colima are installed on this machine (daemon started with `colima start`).
- Fallback (only if Docker proves unworkable): native build on the mini per
  the authors' macOS docs.

### 3.2 Python bridge (`pbg_tstmd/`)

- `runtime.py` — `TstmdSession`: `docker run` the image; write a **YMMSL
  settings** file from parameters (crosslink number, fiber stiffness `k`,
  cell target volume, adhesion count, MCS, `state_output_interval`);
  launch `muscle_manager --start-all ymmsl/adhesions.ymmsl ...`; await the
  **pickle state dumps** in the shared workdir; parse; teardown. Preflight
  guards: `docker_available()`, `image_present()`.
- `processes.py` — `TstmdEcmProcess(Process)`: `config_schema` (lattice size,
  `n_cross`, `fiber_stiffness`, `target_volume`, `seed`, `mcs`,
  `output_interval`); `inputs()/outputs()` (cell_area, cell voxel set, fiber
  bead positions, crosslink bond list, time); `update(state, interval)`.
- `analysis.py` — the paper's metrics, computed from state dumps:
  - **Percolation / giant component** (Fig 2B): fibers as nodes, bonds +
    crosslinks as edges; `n_connected / n`.
  - **Reorientation** `γ_i`, `q = Prob(0.9 ≤ γ ≤ 1)`, `q/q0`, and
    `q(r)/q0(r)` vs distance (Fig 4A–C; paper eq. 7).
  - **Densification** `ρ_dens = ρ_close / ρ_far` via interpolated subnodes
    (Fig 4D–E; paper eqs. 8–9), 2D areas in place of 3D volumes.
  - **Cell area** time series + equilibrium value (Figs 2, 3).

### 3.3 Tests, showcase, packaging

- `tests/test_analysis.py` — metrics unit-tested on **synthetic** networks
  (deterministic, no Docker): known giant-component fraction, a network with
  all fibers radial → `q/q0 > 1`, a densified patch → `ρ_dens > 1`.
- `tests/test_processes.py` — **Docker-gated** end-to-end (skips if daemon
  down): one short coupled run returns a monotone-ish contracting area.
- `composites/`, `visualizations.py`, `demo/demo_report.py` + `demo/report.html`,
  `pixi.toml` (native toolchain), `README.md`, `NEXT_STEPS.md`, `pyproject.toml`.
- Package name: `pbg_tstmd` (consistent with existing C++-bridge siblings
  `pbg_chaste`, `pbg_compucell3d`).

## 4. Investigation + three studies

`workspace/investigations/merks-ecm-reciprocity-2d/investigation.yaml`
(schema v2): question, hypothesis, `executive` verdict block,
`scientific_argument`, `members[]`, `acceptance_criteria[]`, `at_a_glance[]`.

Three member studies (schema v3), flat under `workspace/studies/`:

| Study slug | Fig | Sweep | Acceptance band (directional) |
|---|---|---|---|
| `tstmd-contraction-vs-crosslinking` | 2 | `n_cross` | eq. area ↑ with crosslinking; sharp percolation threshold exists; non-crosslinked → ~full contraction |
| `tstmd-fiber-stiffness-drag` | 3 | `k` (±crosslinks) | eq. area ↑ with `k`; non-crosslinked high-`k` retains area (drag); intermediate `k` → ~0 but slower |
| `tstmd-ecm-remodeling` | 4 | `n_cross`, time, `r` | `q/q0 > 1` and ↑ with crosslinking; densification **biphasic** in crosslinking |

Each study: `behavior_tests[]` with `measure`/`pass_if`, `evaluated_by` =
pytest nodeid; acceptance bands as **external JSON oracles**
(`workspace/studies/<slug>/targets/*.json`) encoding the directional
expectations above (monotonic-increasing flags, "threshold exists",
"argmax at interior crosslink value"), **not** absolute 3D numbers.

## 5. Visualizations

`@as_visualization` functions in `viva_cpm_studies/visualizations/`
(e.g. `tstmd_studies.py`), with **real TST-MD 2D output baked** into a
`_tstmd_data.py` helper (workspace convention → renders identically live and
in the published read-only dashboard). Aesthetic target: publication-grade,
**interactive (plotly)** — sweep curves, the percolation threshold, the
biphasic densification plot, `q/q0` vs time & distance — with cell+fiber
network snapshots as insets. Studies reference them via
`visualizations[].address = local:<Name>` and embed standalone HTML via
`embed_visualizations[].url -> workspace/studies/<slug>/viz/*.html`.

## 6. Resources (3 papers)

PDFs → `workspace/references/papers/`; `@article{}` → `papers.bib`;
registered in `workspace.yaml:references_pdfs` (`{bib_key, path, sha256}`);
`notes/<key>.md` summaries.

- `keijzer-merks-2026-3d-cpm-ecm.pdf` (target) — from `~/Downloads/2609.02375v1.pdf`.
  bib key `keijzer2026cpmecm3d`.
- `keijzer-2025-focal-adhesions.pdf` (context, ref [39]) — from
  `~/Downloads/fcell-12-1462277.pdf`. bib key `keijzer2025focaladhesions`.
- `tsingos-2023-cpm-beadspring.pdf` (**the engine**, ref [38]) — fetch the
  open bioRxiv PDF (2022.06.10.495667); cite the published Biophys. J. DOI
  `10.1016/j.bpj.2023.05.013`. bib key `tsingos2023beadspring`.

## 7. Milestones & the one hard gate

- **M0 — GO/NO-GO:** Docker image builds; one coupled run emits a pickle
  state dump; `analysis.py` extracts cell-area-over-time. If this cannot be
  made to work after reasonable effort → **stop and report** (fallbacks:
  native build on the mini; reduce scope). No faked/mocked dynamics.
- **M1** `viva-tstmd` bridge complete (Process + analysis + tests green;
  Docker-gated test passes locally).
- **M2** Fig 2 study. **M3** Fig 3 study. **M4** Fig 4 study — each baking
  real runs + wiring `behavior_tests` to acceptance bands.
- **M5** visualizations + Resources + `investigation.yaml`; `/viva-report`
  render; investigation PR (draft, `investigation:` prefix, NOT a merge
  target per AGENTS.md).

## 8. Out of scope

3D; dynamic focal adhesions ([39]); Fig 1 schematic; quantitative/pixel
identity; paper-scale HPC sweeps (we use smaller CPU-tractable
lattices / seed counts, documented per study).

## 9. Process / conventions

- All viva-cpm commits in the dedicated worktree
  `~/code/viva-cpm--merks-ecm-reciprocity-2d`; the shared
  `~/code/viva-cpm` checkout stays read-only. `viva-tstmd/` is a fresh repo.
- No AI attribution in commits/PRs. Never auto-merge. Investigation branch
  is a living integration branch, not a merge target.
- Editable-install caveat: run tests with the worktree on `PYTHONPATH` (or
  `-e` install the worktree) so they hit this tree, not the canonical one.
