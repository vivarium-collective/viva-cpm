# viva-cpm

<!-- BEGIN dashboard -->
> ## 📊 [**Live dashboard →**](https://vivarium-collective.github.io/viva-cpm/dashboard/)
> Browse every investigation & study interactively, or read the [published investigation reports](https://vivarium-collective.github.io/viva-cpm/). Auto-published from `main` on every merge.
<!-- END dashboard -->

A [process-bigraph](https://github.com/vivarium-collective/process-bigraph) **Cellular Potts Model** framework — a fast Rust CPM engine (2D/3D, thousands of cells) with a Python layer for pluggable subcellular models, structural constraints, schema-driven world construction, and analysis metrics. A modern, composable remake of CompuCell3D built to do better in 3D, and a **research workspace** where CPM models are wrapped as typed processes, composed, run as studies, and graded against acceptance-criteria tests.

## ▶ Live viewer

**[Explore the demos interactively in your browser →](https://vivarium-collective.github.io/viva-cpm/dashboard/)**
A living 3D colonic crypt (stem cells dividing at the base, differentiating, and sloughing at the mouth), 3D cell sorting, chemotaxis, growth & division, the structural-integrity constraints, and real tissue initialized from Human Reference Atlas / MIBI-TOF imaging — each rotatable, scrubbable, and cell-inspectable.

## Install

The importable engine is `cpm` (a compiled Rust extension) and the research package is `pbg_cpm_studies`. Install from the repo:

```bash
# with uv (recommended)
uv pip install "viva-cpm @ git+https://github.com/vivarium-collective/viva-cpm.git"

# extras: [sbml] SBML/ODE subcellular models · [ftu] Human Reference Atlas FTU→CPM · [all] everything
uv pip install "viva-cpm[all] @ git+https://github.com/vivarium-collective/viva-cpm.git"
```

From source (editable, requires a Rust toolchain + [maturin](https://www.maturin.rs)):

```bash
python -m venv .venv && source .venv/bin/activate
pip install maturin
maturin develop -m crates/cpm-py/Cargo.toml   # builds cpm.cpm_core
pytest                                          # Python suite; `cargo test` for the Rust core
```

## Use the engine from another project

```python
from cpm import load_world, cpm_core

spec = {
    "potts": {"dims": [50, 50, 1], "boundary": "periodic",
              "neighbor_order": 2, "temperature": 12.0, "seed": 0},
    "cells": [
        {"type": 1, "target_volume": 25, "lambda_volume": 1.0,
         "target_surface": 0, "lambda_surface": 0, "seed_block": [5, 5, 0, 13, 13, 1]},
    ],
    "contact": [{"a": 0, "b": 1, "j": 12.0}],
}
world = load_world(spec)
world.step(100)                     # run 100 Monte-Carlo sweeps
print(world.cell_volumes())
```

The engine itself is `cpm.cpm_core` (a compiled Rust extension). `load_world` builds a
world from a plain dict spec (cells or a seeded label array, contact energies, diffusion
fields, connectivity, basement membrane). Chemotaxis can operate on the raw field or, via
`set_chemotaxis_occupancy`, in receptor-**occupancy** space — the substrate for fold-change
detection (see the recruitment investigation below).

## Process-bigraph composites

Cells are wired as process-bigraph processes via import-path addresses, so any
process-bigraph `Composite` can embed them:

- `local:!cpm.processes.cpm_process.CPMProcess` — the CPM step as a process
- `local:!cpm.subcellular.sbml.SBMLSubcell` — a per-cell SBML/ODE model (needs `[sbml]`)
- `local:!cpm.subcellular.boolean.BooleanSubcell` — a per-cell Boolean fate network
- `local:!cpm.subcellular.adaptive_receptor.AdaptiveReceptorSubcell` — a per-cell receptor with slow adaptation

See `cpm/composites/crypt.py` for a full crypt-differentiation composite (CPM + SBML
stemness ODE + Boolean fate switch), run with the process-bigraph `Composite` engine.

## Research workspace: investigations & studies

`workspace/` is a process-bigraph research workspace: composites in `pbg_cpm_studies/`,
studies under `workspace/studies/`, grouped into **investigations**. Each study carries a
model, readouts, simulation runs, and acceptance-criteria **behavior tests** that grade a
run into a signed pass/fail verdict. Two investigations ship today (browse them live on the
[dashboard](https://vivarium-collective.github.io/viva-cpm/dashboard/)):

- **glazier-graner-1993** — an 11-study reproduction of the classic Glazier–Graner
  differential-adhesion results (annealing, global equilibration, checkerboard, cell
  sorting, engulfment, position reversal, partial sorting, dispersal, vacancy nucleation).
- **chemotactic-recruitment** — a secreted cue recruits responder cells, realized at three
  levels: a phenomenological chemotaxis-λ (baseline + inhibited + adversarial controls), a
  Kd-calibrated receptor-occupancy model (receptor-baseline + blocked), and an **adaptive**
  fold-change-detection refinement built by an agentic model-building loop. That loop —
  author a contract of tests → audit → feasibility spike → lock → build/run/evaluate →
  navigate — climbs an emergent mechanism ladder (`static → hill_occupancy → adaptive`) in
  which occupancy-space chemotaxis makes the fixed-`kd` rung collapse at high background and
  adaptation rescues it; the run is captured as a `model_build_trajectory`. The
  calibration tooling (`pbg_cpm_studies/model_building/calibrate.py`) is a sensitivity
  screen + common-random-numbers + refine, not a hand grid.

The loop, contract, audit, and grading machinery live in
[viva-superpowers](https://github.com/vivarium-collective/viva-superpowers); this repo is one
of its research workspaces.

## Structural constraints

- **Connectivity** (E1): forbids copy attempts that would fragment a cell or pinch off
  interior medium (gaps). `spec["connectivity"] = {"types": [1, 2], "medium": true}`.
- **Basement membrane** (E3a): a basal anchor energy keeping epithelial cells in a thin
  band hugging a fixed membrane surface. `spec["membrane"] = {"anchors": [...], "k": ...,
  "band": ..., "types": [...]}`.

## Layout

```
crates/            Rust workspace: cpm-core (engine) + cpm-py (pyo3 bindings → cpm.cpm_core)
cpm/               Python framework: schema, processes, subcellular, composites, metrics, ftu
pbg_cpm_studies/   research package: composites, model-building mechanisms + calibrate, visualizations
workspace/         the research workspace: studies/, investigations/, references/, reports/
demos/             runnable demos (each validates + exports a viewer model)
viewer/            browser 2D/3D viewer for the exported models
docs/              specs & implementation plans
tests/             Rust (cargo test) + Python (pytest) suites
```

## License

MIT
