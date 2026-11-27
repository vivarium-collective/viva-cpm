# viva-cpm

[![release](https://github.com/vivarium-collective/viva-cpm/actions/workflows/release.yml/badge.svg)](https://github.com/vivarium-collective/viva-cpm/actions/workflows/release.yml)
[![PyPI](https://img.shields.io/pypi/v/viva-cpm)](https://pypi.org/project/viva-cpm/)
[![license](https://img.shields.io/github/license/vivarium-collective/viva-cpm)](./LICENSE)

A [process-bigraph](https://github.com/vivarium-collective/process-bigraph) **Cellular Potts Model** engine — a fast Rust CPM core (2D/3D, thousands of cells) with a Python layer for pluggable subcellular models, structural constraints, schema-driven world construction, and analysis metrics. A modern, composable remake of CompuCell3D built to do better in 3D, with CPM cells wrapped as typed process-bigraph processes.

> **Studies & investigations live in a separate repo.** The research workspace built on
> this engine — the differential-adhesion (Glazier–Graner 1993), chemotactic-recruitment,
> influenza (Sego 2022) and endosymbiosis investigations — now lives in
> [**viva-cpm-studies**](https://github.com/vivarium-collective/viva-cpm-studies), which
> pins this engine by version. This repo is the engine only.

## Install

```bash
# from PyPI (prebuilt wheels: CPython 3.12 / 3.13, manylinux x86_64/aarch64, macOS arm64)
uv pip install viva-cpm

# extras: [sbml] SBML/ODE subcellular models · [ftu] Human Reference Atlas FTU→CPM · [all] everything
uv pip install "viva-cpm[all]"
```

From source (editable, requires a Rust toolchain + [maturin](https://www.maturin.rs)):

```bash
python -m venv .venv && source .venv/bin/activate
pip install maturin
maturin develop -m crates/cpm-py/Cargo.toml   # builds viva_cpm.cpm_core
pytest                                          # Python suite; `cargo test` for the Rust core
```

## Use the engine

```python
from viva_cpm import load_world, cpm_core

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

The engine itself is `viva_cpm.cpm_core` (a compiled Rust extension). `load_world` builds a
world from a plain dict spec (cells or a seeded label array, contact energies, diffusion
fields, connectivity, basement membrane). Chemotaxis can operate on the raw field or, via
`set_chemotaxis_occupancy`, in receptor-**occupancy** space.

## Process-bigraph composites

Cells are wired as process-bigraph processes via import-path addresses, so any
process-bigraph `Composite` can embed them:

- `local:!viva_cpm.processes.cpm_process.CPMProcess` — the CPM step as a process
- `local:!viva_cpm.subcellular.sbml.SBMLSubcell` — a per-cell SBML/ODE model (needs `[sbml]`)
- `local:!viva_cpm.subcellular.boolean.BooleanSubcell` — a per-cell Boolean fate network
- `local:!viva_cpm.subcellular.adaptive_receptor.AdaptiveReceptorSubcell` — a per-cell receptor with slow adaptation

See `viva_cpm/composites/crypt.py` for a full crypt-differentiation composite (CPM + SBML
stemness ODE + Boolean fate switch), run with the process-bigraph `Composite` engine.

## Structural constraints

- **Connectivity** (E1): forbids copy attempts that would fragment a cell or pinch off
  interior medium (gaps). `spec["connectivity"] = {"types": [1, 2], "medium": true}`.
- **Basement membrane** (E3a): a basal anchor energy keeping epithelial cells in a thin
  band hugging a fixed membrane surface. `spec["membrane"] = {"anchors": [...], "k": ...,
  "band": ..., "types": [...]}`.

## Layout

```
crates/   Rust workspace: cpm-core (engine) + cpm-py (pyo3 bindings → viva_cpm.cpm_core)
viva_cpm/ Python layer: schema, processes, subcellular, composites, metrics, ftu
demos/    runnable demos (each validates + exports a viewer model)
viewer/   browser 2D/3D viewer for the exported models
bench/    Rust benchmarks
tests/    Rust (cargo test) + Python (pytest) suites
```

## License

MIT
