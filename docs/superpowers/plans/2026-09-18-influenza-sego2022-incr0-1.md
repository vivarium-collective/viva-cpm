# Influenza (Sego 2022) reproduction — Increments 0 & 1 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Establish the cited parameter authority + reproduction targets (Increment 0) and a validated, performant epithelial-sheet CPM substrate at paper scale (Increment 1) for the native reproduction of Sego et al. 2022.

**Architecture:** A new `viva_cpm_studies/influenza/` package holds the canonical parameters (transcribed from the paper's CompuCell3D source), the digitized figure targets, and the epithelial-sheet builder that emits a `load_world` spec for the existing Rust `cpm_core.World`. No engine (Rust) changes in these two increments — the sheet uses primitives that already exist (`add_cell`, `set_contact`, volume constraint, `step_parallel`).

**Tech Stack:** Python 3.12, `cpm.cpm_core` (Rust/PyO3 CPM engine), `process_bigraph`, PyYAML, pytest, numpy. viva-superpowers dashboard skills (`/viva-investigation`, `/viva-study`, `/viva-cite-bands`) for the investigation shell.

**Spec:** `docs/superpowers/specs/2026-09-18-influenza-sego2022-reproduction-design.md`

## Global Constraints

- Reproduction target: quantitative match of Figs 3B / 5 / 7 at paper scale (1 mm² = 500×500-site 2 µm lattice, ~10k epithelial cells, 50 replicas). "250 k" in the paper is the ODE organism population for η/θ scaling, NOT the spatial cell count.
- Oracle: the paper's supplementary CompuCell3D source is the authoritative parameter/spec; published figures are the numeric acceptance targets.
- Fidelity convention (from `docs/cc3d-reference/demo-parameters.md`): where a per-version CC3D constant cannot be retrieved verbatim, use the CC3D-typical value, **mark it**, and make the *validated behavior* (not the exact constant) the fidelity criterion.
- Lattice: discretization length 2 µm; epithelial cells are uniform 5×5 sites (10 µm diameter). Volume constraint v_c = 100 µm² = 25 sites/cell; λ_volume = 9; intrinsic random motility (CPM temperature) H* = 10; neighbor order 2 (second-order Manhattan); BCs Neumann (no-flux) + periodic; Δt = 1 min/step.
- No AI attribution in commits or PRs (per repo policy). Commit messages end with no `Co-Authored-By` trailer.
- Work stays on branch `investigation/influenza-sego2022` in worktree `~/code/viva-cpm--influenza-sego2022`. One PR per increment.
- New cell-type integer codes (fixed across the whole investigation): `MEDIUM=0, H=1 (uninfected), I=2 (infected), D=3 (dead), M=4 (macrophage), K=5 (NK), E=6 (CD8+ T)`.

---

# Increment 0 — Spec authority & reproduction targets

No engine code. Deliverables: a documented CC3D-source acquisition note, a human-readable parameter transcription, a machine-readable cited `params.yaml` + loader, digitized figure targets, and the investigation shell with a `parameter-provenance` study.

### Task 0.1: Acquire & record the CompuCell3D source

**Files:**
- Create: `docs/cc3d-reference/sego2022-source-notes.md`

**Interfaces:**
- Produces: a documented location/availability of the authoritative CC3D model source, referenced by Tasks 0.2–0.3.

- [ ] **Step 1: Locate the supplementary source.** The paper (Appendix A) states the CC3D model is in the Supplementary Materials, package "Source Code". Try, in order: (a) the article's ScienceDirect supplementary at the DOI `10.1016/j.jtbi.2021.110918`; (b) the CompuCell3D cellularization models by T.J. Sego on GitHub (`github.com/tjsego`) and nanoHUB; (c) the Price et al. 2015 ODE model for the global compartment. Use WebFetch/WebSearch.

- [ ] **Step 2: Record findings.** Write `sego2022-source-notes.md` with: the exact URL(s) found, what the package contains (CC3DML + Python steppables? Antimony? parameter files?), and whether the exact per-parameter constants are retrievable. If the source is NOT retrievable, state that explicitly and note that parameters will be transcribed from the paper's Tables 1–4 and equations, marked CC3D-typical where the paper is silent (per the fidelity convention).

- [ ] **Step 3: Commit.**
```bash
git add docs/cc3d-reference/sego2022-source-notes.md
git commit -m "docs(influenza): record CC3D source acquisition for Sego 2022"
```

### Task 0.2: Transcribe parameters to the cc3d-reference doc

**Files:**
- Create: `docs/cc3d-reference/sego2022-parameters.md`

**Interfaces:**
- Produces: the human-readable, sourced transcription that Task 0.3's `params.yaml` mirrors.

- [ ] **Step 1: Write the transcription.** Following the structure of `docs/cc3d-reference/demo-parameters.md`, transcribe, each with its source (paper table/eq number, or CC3D source file, or "CC3D-typical, marked"):
  - CPM: lattice 2 µm, cell 5×5 sites, v_c=100 µm²=25 sites, λ_volume=9, temperature H*=10, neighbor order 2, BCs, Δt=1 min (Table 4; Table 3).
  - Adhesion J matrix (Table 3): uninfected–immune 20, infected–immune 10, dead–immune 20, homotypic-immune 25, heterotypic-immune 10. **Flag as gap:** epithelial–epithelial and epithelial–medium J are not in Table 3 → extract from the CC3D source (Task 0.1) or mark CC3D-typical.
  - Diffusion coefficients & lengths (Table 3): virus 0.0119 µm²/s (5 cell diam), chemokines 1.04 (10 cell diam), type-I IFN 0.520 (2 cell diam), IL-10 0.327 (2 cell diam). Record the cell-diameter unit conversion.
  - Chemotaxis λ (Table 3): macrophage–virus 5000, NK–chemokines 5000, CD8⁺–chemokines 10000. Record the paper's saturating/COM chemotaxis functional form (`λc·c/(1+c_CM)`) vs the engine's current linear form.
  - Field PDE source/decay terms (Table 1) for v, f, c, l.
  - Type-transition rates (Table 2): infection H→I, uninfected death, infected death (incl. contact killing γ terms), recovery D→H (Allee a_H).
  - Resistance ρ = 1 − f̄/(θ·a_rf + f̄).
  - Cellularization scaling: θ = 4×10⁻⁸ (from 250k ODE cells + v_c), η = 0.0049 (0.3 mm) / 0.04 (1.0 mm); seeding fraction 1%; local fractions M 100%, NK 75%, CD8⁺ 75% (Table 4).
  - Scenarios: initial viral load (1/10/100/1000/10000), initial infection fraction (0.001/0.005/0.01/0.05); patches 0.3 mm & 1.0 mm.

- [ ] **Step 2: Commit.**
```bash
git add docs/cc3d-reference/sego2022-parameters.md
git commit -m "docs(influenza): transcribe Sego 2022 parameters (Tables 1-4)"
```

### Task 0.3: Machine-readable parameters + loader + type constants

**Files:**
- Create: `viva_cpm_studies/influenza/__init__.py`
- Create: `viva_cpm_studies/influenza/types.py`
- Create: `viva_cpm_studies/influenza/params.yaml`
- Create: `viva_cpm_studies/influenza/params.py`
- Test: `tests/test_influenza_params.py`

**Interfaces:**
- Produces:
  - `viva_cpm_studies.influenza.types`: `MEDIUM=0, H=1, I=2, D=3, M=4, K=5, E=6` (ints).
  - `viva_cpm_studies.influenza.params.load_params() -> dict` — parsed `params.yaml`.
  - `viva_cpm_studies.influenza.params.PARAMS` — module-level cached dict.

- [ ] **Step 1: Write the failing test.**
```python
# tests/test_influenza_params.py
from viva_cpm_studies.influenza import params, types


def test_type_codes_are_distinct_and_ordered():
    codes = [types.MEDIUM, types.H, types.I, types.D, types.M, types.K, types.E]
    assert codes == [0, 1, 2, 3, 4, 5, 6]
    assert len(set(codes)) == 7


def test_params_load_core_cpm_values():
    p = params.load_params()
    cpm = p["cpm"]
    assert cpm["lattice_um"] == 2
    assert cpm["cell_sites"] == 25            # 5x5
    assert cpm["volume_constraint_um2"] == 100
    assert cpm["lambda_volume"] == 9
    assert cpm["temperature"] == 10
    assert cpm["neighbor_order"] == 2
    assert cpm["dt_min"] == 1


def test_params_carry_provenance_for_every_top_section():
    p = params.load_params()
    # every parameter section must cite a source string (paper table/eq or CC3D source)
    for section in ("cpm", "adhesion", "diffusion", "chemotaxis", "scaling"):
        assert p[section].get("source"), f"{section} missing provenance"


def test_diffusion_coefficients_match_table3():
    d = params.load_params()["diffusion"]
    assert d["virus_um2_s"] == 0.0119
    assert d["chemokine_um2_s"] == 1.04
    assert d["ifn_type1_um2_s"] == 0.520
    assert d["il10_um2_s"] == 0.327
```

- [ ] **Step 2: Run test to verify it fails.**
Run: `PYTHONPATH=~/code/viva-cpm--influenza-sego2022 pytest tests/test_influenza_params.py -v`
Expected: FAIL (module `viva_cpm_studies.influenza` not found).

- [ ] **Step 3: Create `types.py`.**
```python
# viva_cpm_studies/influenza/types.py
"""Cell-type integer codes for the influenza-sego2022 investigation.

Fixed across every increment so composites, transitions and viz agree.
"""
MEDIUM = 0
H = 1  # uninfected epithelial
I = 2  # infected epithelial
D = 3  # dead epithelial
M = 4  # macrophage
K = 5  # NK cell
E = 6  # CD8+ T cell
```

- [ ] **Step 4: Create `params.yaml`** (mirror the transcription in Task 0.2; every section carries a `source:` string). Minimal shape the test requires:
```yaml
# viva_cpm_studies/influenza/params.yaml
# Canonical parameters for the Sego et al. 2022 influenza reproduction.
# Mirrors docs/cc3d-reference/sego2022-parameters.md. Every section cites a source.
cpm:
  source: "Sego 2022 Tables 3-4"
  lattice_um: 2
  cell_sites: 25          # 5x5 sites; 10 um diameter
  volume_constraint_um2: 100
  lambda_volume: 9
  temperature: 10         # intrinsic random motility H*
  neighbor_order: 2
  dt_min: 1
adhesion:                 # contact energies J (Table 3)
  source: "Sego 2022 Table 3"
  uninfected_immune: 20
  infected_immune: 10
  dead_immune: 20
  homotypic_immune: 25
  heterotypic_immune: 10
  epithelial_epithelial: null   # TODO(Task 0.1): from CC3D source or mark CC3D-typical
  epithelial_medium: null       # TODO(Task 0.1): from CC3D source or mark CC3D-typical
diffusion:
  source: "Sego 2022 Table 3 (de Jong et al. 2006)"
  virus_um2_s: 0.0119
  chemokine_um2_s: 1.04
  ifn_type1_um2_s: 0.520
  il10_um2_s: 0.327
chemotaxis:
  source: "Sego 2022 Table 3"
  functional_form: "lambda_c * c / (1 + c_CM)  # saturating/COM; engine currently linear"
  macrophage_virus: 5000
  nk_chemokine: 5000
  cd8_chemokine: 10000
scaling:
  source: "Sego 2022 p.6, Table 4"
  ode_epithelial_population: 250000
  theta: 4.0e-8
  eta_0p3mm: 0.0049
  eta_1p0mm: 0.04
  seeding_fraction: 0.01
  local_fraction_macrophage: 1.00
  local_fraction_nk: 0.75
  local_fraction_cd8: 0.75
scenarios:
  source: "Sego 2022 Section 3"
  initial_viral_load: [1, 10, 100, 1000, 10000]
  initial_infection_fraction: [0.001, 0.005, 0.01, 0.05]
  patch_mm: [0.3, 1.0]
```
Note: the `null` adhesion entries are the one legitimate TODO carried forward — they are a data gap owned by Task 0.1, not a plan placeholder. If Task 0.1 retrieved them, fill them in and drop the TODO comment.

- [ ] **Step 5: Create `params.py` loader.**
```python
# viva_cpm_studies/influenza/params.py
"""Load the canonical Sego-2022 parameter set from params.yaml."""
from __future__ import annotations
from functools import lru_cache
from pathlib import Path
import yaml

_PARAMS_PATH = Path(__file__).with_name("params.yaml")


@lru_cache(maxsize=1)
def load_params() -> dict:
    return yaml.safe_load(_PARAMS_PATH.read_text())


PARAMS = load_params()
```

- [ ] **Step 6: Create `__init__.py`.**
```python
# viva_cpm_studies/influenza/__init__.py
"""Native viva-cpm reproduction of Sego et al. 2022 influenza model."""
```

- [ ] **Step 7: Run test to verify it passes.**
Run: `PYTHONPATH=~/code/viva-cpm--influenza-sego2022 pytest tests/test_influenza_params.py -v`
Expected: PASS (4 tests).

- [ ] **Step 8: Commit.**
```bash
git add viva_cpm_studies/influenza/ tests/test_influenza_params.py
git commit -m "feat(influenza): cited canonical parameter set + type codes"
```

### Task 0.4: Digitize figure targets + loader

**Files:**
- Create: `viva_cpm_studies/influenza/targets/fig3b.json`
- Create: `viva_cpm_studies/influenza/targets/fig5.json`
- Create: `viva_cpm_studies/influenza/targets/fig7.json`
- Create: `viva_cpm_studies/influenza/targets.py`
- Test: `tests/test_influenza_targets.py`

**Interfaces:**
- Consumes: `viva_cpm_studies/influenza/params.py` (scenario lists).
- Produces: `viva_cpm_studies.influenza.targets.load_target(name) -> dict` with, per observable, a list of `{t_days, value, lo, hi}` band points.

- [ ] **Step 1: Write the failing test.**
```python
# tests/test_influenza_targets.py
import pytest
from viva_cpm_studies.influenza import targets


@pytest.mark.parametrize("name", ["fig3b", "fig5", "fig7"])
def test_target_loads_and_has_bands(name):
    t = targets.load_target(name)
    assert t["figure"] == name
    assert t["observables"], "at least one observable"
    for obs, series in t["observables"].items():
        assert series, f"{obs} has points"
        for pt in series:
            assert pt["lo"] <= pt["value"] <= pt["hi"], f"{obs} band brackets value"
            assert pt["t_days"] >= 0


def test_fig3b_covers_key_observables():
    obs = targets.load_target("fig3b")["observables"]
    for key in ("uninfected_cells", "infected_cells", "extracellular_virus"):
        assert key in obs
```

- [ ] **Step 2: Run test to verify it fails.**
Run: `PYTHONPATH=~/code/viva-cpm--influenza-sego2022 pytest tests/test_influenza_targets.py -v`
Expected: FAIL (module `targets` not found).

- [ ] **Step 3: Digitize the targets.** From the paper PDF (`workspace/references/papers/sego-2022-influenza.pdf`), read representative points off each panel and record them as bands. Use the ODE reference curve (black line) as `value` and set `lo/hi` to the spatial replica ensemble spread shown (colored band); where only qualitative, use a generous ±0.5 log band and mark `"soft": true`. Fig 3B panels: uninfected/infected/dead cells, extracellular virus, APCs, macrophages, NK, CD8⁺, chemokines, type-I IFN, antibodies, IL-10 over 0–3.5 days. Fig 5: uninfected/infected/virus/antibodies across viral loads 1–10000 with the lethal-threshold note. Fig 7: same observables across infection fractions 0.001–0.05. JSON shape:
```json
{
  "figure": "fig3b",
  "scenario": {"initial_infection_fraction": 0.05, "patch_mm": 0.3, "replicas": 50},
  "observables": {
    "uninfected_cells": [
      {"t_days": 0.0, "value": 1225, "lo": 1100, "hi": 1300},
      {"t_days": 1.0, "value": 800,  "lo": 500,  "hi": 1000}
    ]
  }
}
```
(Populate all panels; the two points above are illustrative of the shape, not the full series.)

- [ ] **Step 4: Create `targets.py` loader.**
```python
# viva_cpm_studies/influenza/targets.py
"""Load digitized figure targets (acceptance bands) for the reproduction studies."""
from __future__ import annotations
import json
from pathlib import Path

_DIR = Path(__file__).with_name("targets")


def load_target(name: str) -> dict:
    return json.loads((_DIR / f"{name}.json").read_text())
```

- [ ] **Step 5: Run test to verify it passes.**
Run: `PYTHONPATH=~/code/viva-cpm--influenza-sego2022 pytest tests/test_influenza_targets.py -v`
Expected: PASS.

- [ ] **Step 6: Commit.**
```bash
git add viva_cpm_studies/influenza/targets* tests/test_influenza_targets.py
git commit -m "feat(influenza): digitized Fig 3B/5/7 reproduction targets + loader"
```

### Task 0.5: Scaffold the investigation shell + parameter-provenance study

**Files:**
- Create (via skills): `workspace/investigations/influenza-sego2022/investigation.yaml` and a `parameter-provenance` study.

**Interfaces:**
- Consumes: the reference `sego2022influenza`, `params.yaml`, and the targets.
- Produces: the dashboard investigation the later increments' studies attach to.

- [ ] **Step 1: Ensure the dashboard server is up** (`/viva-workbench status`; `/viva-workbench start` if not) — precondition for the study/investigation skills.

- [ ] **Step 2: Create the investigation** with `/viva-investigation new influenza-sego2022`, setting: title "Influenza infection & immune response (Sego 2022 reproduction)"; question = native reproduction of the CC3D cellularized influenza model to quantitative figure match; hypothesis + lead per the spec §1. Set the overview to reference the increment ladder (spec §5).

- [ ] **Step 3: Create the `parameter-provenance` study** with `/viva-study` in the Design phase: it documents the cited `params.yaml`, the CC3D-source note, and the digitized targets as the reproduction's parameter authority — no simulation. Link the `sego2022influenza` reference.

- [ ] **Step 4: Attach reference provenance** with `/viva-cite-bands`, linking `sego2022influenza` to the target acceptance bands.

- [ ] **Step 5: Commit** whatever YAML the skills wrote.
```bash
git add workspace/investigations/influenza-sego2022 workspace.yaml
git commit -m "feat(influenza): scaffold investigation + parameter-provenance study"
```

---

# Increment 1 — Epithelial sheet at scale

Build and validate the confluent epithelial CPM substrate (no fields, no immune cells yet), and establish the MCS-throughput budget that makes the capstone feasible. Uses only existing engine primitives.

### Task 1.1: Epithelial-sheet spec builder

**Files:**
- Create: `viva_cpm_studies/influenza/sheet.py`
- Test: `tests/test_influenza_sheet.py`

**Interfaces:**
- Consumes: `params.load_params()`, `types` (MEDIUM, H).
- Produces: `viva_cpm_studies.influenza.sheet.build_sheet_spec(patch_mm: float, seed: int = 17) -> dict` returning a `load_world` spec dict `{"potts": {...}, "cells": [...], "contact": [...]}` — same shape consumed by `cpm.processes.cpm_process` composites (see `viva_cpm_studies/composites/gg1993.py:build_spec`). All epithelial cells are type `H`, uniform 5×5-site blocks tiling the domain, `target_volume=25`, `lambda_volume=9`.

- [ ] **Step 1: Write the failing test.**
```python
# tests/test_influenza_sheet.py
from viva_cpm_studies.influenza import sheet, types


def test_sheet_lattice_and_cell_count_0p3mm():
    spec = sheet.build_sheet_spec(0.3)
    # 0.3 mm / 2 um = 150 sites per side; 5x5 cells -> 30 per side -> 900 cells
    assert spec["potts"]["dims"] == [150, 150, 1]
    assert len(spec["cells"]) == 900
    assert all(c["type"] == types.H for c in spec["cells"])
    assert all(c["target_volume"] == 25.0 for c in spec["cells"])
    assert all(c["lambda_volume"] == 9.0 for c in spec["cells"])


def test_sheet_lattice_and_cell_count_1mm():
    spec = sheet.build_sheet_spec(1.0)
    # 1.0 mm / 2 um = 500 sites per side; 5x5 cells -> 100 per side -> 10000 cells
    assert spec["potts"]["dims"] == [500, 500, 1]
    assert len(spec["cells"]) == 10000


def test_sheet_potts_uses_paper_cpm_constants():
    p = sheet.build_sheet_spec(0.3)["potts"]
    assert p["temperature"] == 10
    assert p["neighbor_order"] == 2
    assert p["boundary"] in ("noflux", "periodic")
```

- [ ] **Step 2: Run test to verify it fails.**
Run: `PYTHONPATH=~/code/viva-cpm--influenza-sego2022 pytest tests/test_influenza_sheet.py -v`
Expected: FAIL (module `sheet` not found).

- [ ] **Step 3: Implement `sheet.py`.**
```python
# viva_cpm_studies/influenza/sheet.py
"""Confluent epithelial sheet as a load_world spec for the Rust CPM engine.

A uniform grid of 5x5-site (10 um) cells of type H (uninfected), tiling a
square patch. No fields or immune cells yet (Increment 1). Adhesion for the
epithelial-only sheet keeps cells confluent; the immune J matrix (Table 3)
arrives in Increments 6-7.
"""
from __future__ import annotations

from .params import load_params
from . import types

CELL_SIDE_SITES = 5  # 5x5 = 25 sites = 10 um cell


def _sites_per_side(patch_mm: float, lattice_um: int) -> int:
    return int(round(patch_mm * 1000.0 / lattice_um))


def build_sheet_spec(patch_mm: float, seed: int = 17) -> dict:
    p = load_params()["cpm"]
    n = _sites_per_side(patch_mm, p["lattice_um"])          # sites per side
    cells_per_side = n // CELL_SIDE_SITES
    cells = []
    for gy in range(cells_per_side):
        for gx in range(cells_per_side):
            x0, y0 = gx * CELL_SIDE_SITES, gy * CELL_SIDE_SITES
            cells.append({
                "type": types.H,
                "target_volume": float(p["cell_sites"]),    # 25
                "lambda_volume": float(p["lambda_volume"]),  # 9
                "target_surface": 0.0, "lambda_surface": 0.0,
                "seed_block": [x0, y0, 0,
                               x0 + CELL_SIDE_SITES, y0 + CELL_SIDE_SITES, 1],
            })
    # Epithelial-only adhesion: cohesive H-H, higher H-medium so the sheet
    # stays confluent. Values marked CC3D-typical pending Task 0.1 extraction.
    contact = [
        {"a": types.MEDIUM, "b": types.MEDIUM, "j": 0.0},
        {"a": types.MEDIUM, "b": types.H, "j": 16.0},   # CC3D-typical, marked
        {"a": types.H, "b": types.H, "j": 4.0},         # CC3D-typical, marked
    ]
    return {
        "potts": {"dims": [n, n, 1], "boundary": "noflux",
                  "neighbor_order": int(p["neighbor_order"]),
                  "temperature": float(p["temperature"]), "seed": seed},
        "cells": cells,
        "contact": contact,
    }
```

- [ ] **Step 4: Run test to verify it passes.**
Run: `PYTHONPATH=~/code/viva-cpm--influenza-sego2022 pytest tests/test_influenza_sheet.py -v`
Expected: PASS (3 tests).

- [ ] **Step 5: Commit.**
```bash
git add viva_cpm_studies/influenza/sheet.py tests/test_influenza_sheet.py
git commit -m "feat(influenza): epithelial-sheet load_world spec builder"
```

### Task 1.2: Instantiate the world and validate geometry

**Files:**
- Modify: `tests/test_influenza_sheet.py` (add an engine-level test)
- Create: `viva_cpm_studies/influenza/build.py`

**Interfaces:**
- Consumes: `sheet.build_sheet_spec`, `cpm.cpm_core.World`.
- Produces: `viva_cpm_studies.influenza.build.world_from_spec(spec) -> cpm_core.World` (finalized, ready to `step`).

- [ ] **Step 1: Write the failing test** (append to `tests/test_influenza_sheet.py`).
```python
def test_world_builds_and_holds_volume_after_relaxation():
    from viva_cpm_studies.influenza import build, sheet
    spec = sheet.build_sheet_spec(0.3)
    w = build.world_from_spec(spec)
    assert w.n_cells() == 900
    w.step(50)  # short relaxation
    vols = w.cell_volumes()
    mean_v = sum(vols) / len(vols)
    assert 18 <= mean_v <= 32          # ~25 sites, confluent, no collapse
    assert min(vols) > 0               # no cell vanished
```

- [ ] **Step 2: Run test to verify it fails.**
Run: `PYTHONPATH=~/code/viva-cpm--influenza-sego2022 pytest tests/test_influenza_sheet.py::test_world_builds_and_holds_volume_after_relaxation -v`
Expected: FAIL (`build` module not found).

- [ ] **Step 3: Implement `build.py`.** Follow the low-level API used in `tests/test_bindings.py` (`World(dims, boundary, neighbor_order, temperature)`, `add_cell(type, target_volume, lambda_volume, target_surface, lambda_surface)`, `set_contact(a, b, j)`, `seed_block(id, x0,y0,z0,x1,y1,z1)`, `finalize(seed)`).
```python
# viva_cpm_studies/influenza/build.py
"""Instantiate a finalized cpm_core.World from a load_world sheet spec."""
from __future__ import annotations
from cpm import cpm_core


def world_from_spec(spec: dict):
    po = spec["potts"]
    w = cpm_core.World(tuple(po["dims"]), po["boundary"],
                       int(po["neighbor_order"]), float(po["temperature"]))
    for c in spec["contact"]:
        w.set_contact(int(c["a"]), int(c["b"]), float(c["j"]))
    for cell in spec["cells"]:
        cid = w.add_cell(int(cell["type"]), float(cell["target_volume"]),
                         float(cell["lambda_volume"]),
                         float(cell["target_surface"]), float(cell["lambda_surface"]))
        x0, y0, z0, x1, y1, z1 = cell["seed_block"]
        w.seed_block(cid, x0, y0, z0, x1, y1, z1)
    w.finalize(int(po.get("seed", 0)))
    return w
```

- [ ] **Step 4: Run test to verify it passes.**
Run: `PYTHONPATH=~/code/viva-cpm--influenza-sego2022 pytest tests/test_influenza_sheet.py -v`
Expected: PASS (4 tests). If mean volume falls outside the band, adjust the epithelial J values (Step-3 note in Task 1.1) and re-run — this is the confluence calibration, and the *behavior* (cells hold ~25 sites, none vanish) is the fidelity criterion, not the exact J.

- [ ] **Step 5: Commit.**
```bash
git add viva_cpm_studies/influenza/build.py tests/test_influenza_sheet.py
git commit -m "feat(influenza): build finalized CPM world from sheet spec + geometry test"
```

### Task 1.3: Throughput benchmark at 1 mm² (perf budget)

**Files:**
- Create: `tests/test_influenza_perf.py`

**Interfaces:**
- Consumes: `build.world_from_spec`, `sheet.build_sheet_spec`, `World.step_parallel`.
- Produces: the recorded MCS/sec throughput at 1 mm² — the acceptance gate the spec (§7) requires before the capstone.

- [ ] **Step 1: Write the benchmark test.**
```python
# tests/test_influenza_perf.py
import time
import pytest
from viva_cpm_studies.influenza import build, sheet


@pytest.mark.perf
def test_mcs_throughput_1mm():
    w = build.world_from_spec(sheet.build_sheet_spec(1.0))  # 500x500, 10k cells
    w.step_parallel(5, 64)          # warm up
    t0 = time.perf_counter()
    steps = 50
    w.step_parallel(steps, 64)
    rate = steps / (time.perf_counter() - t0)
    print(f"\n1mm^2 throughput: {rate:.2f} MCS/s")
    # Budget floor: the 2-week capstone is ~20160 steps. At >=2 MCS/s one
    # replica finishes in <~2.8 h; assert a conservative floor so a regression
    # is caught. Tune the floor to the measured baseline on first run.
    assert rate >= 2.0, f"throughput {rate:.2f} MCS/s below budget floor"
```

- [ ] **Step 2: Register the `perf` marker** in `pyproject.toml` under `[tool.pytest.ini_options]` (add `markers = ["perf: performance budget tests"]` if not present) so `pytest` does not warn.

- [ ] **Step 3: Run the benchmark.**
Run: `PYTHONPATH=~/code/viva-cpm--influenza-sego2022 pytest tests/test_influenza_perf.py -v -s -m perf`
Expected: PASS, and it prints the measured MCS/s. If the measured rate is far above 2.0, raise the floor to ~50% of measured to make it a real regression gate; if below 2.0, record the number in the study and note the Rust hot-path optimization (spec §7) as a follow-up before the capstone.

- [ ] **Step 4: Commit.**
```bash
git add tests/test_influenza_perf.py pyproject.toml
git commit -m "test(influenza): 1mm^2 MCS throughput budget gate"
```

### Task 1.4: Composite factory + `epithelial-sheet-baseline` study

**Files:**
- Create: `viva_cpm_studies/composites/influenza.py`
- Test: `tests/test_influenza_composite.py`
- Create (via skills): the `epithelial-sheet-baseline` study in the investigation.

**Interfaces:**
- Consumes: `sheet.build_sheet_spec`, the `CPMProcess` address `local:!cpm.processes.cpm_process.CPMProcess` (see `viva_cpm_studies/composites/gg1993.py`).
- Produces: `viva_cpm_studies.composites.influenza.epithelial_sheet_baseline() -> dict` — a process-bigraph composite document the dashboard can run; `build_spec(patch_mm)` for a modest live-demo aggregate.

- [ ] **Step 1: Write the failing test.**
```python
# tests/test_influenza_composite.py
from viva_cpm_studies.composites import influenza as inf


def test_composite_document_runs_and_emits():
    doc = inf.epithelial_sheet_baseline()
    assert isinstance(doc, dict)
    # smoke: the demo composite builds a modest sheet, not the full 1mm^2
    spec = inf.build_spec(patch_mm=0.1)   # 50x50 sites, 100 cells
    assert spec["potts"]["dims"] == [50, 50, 1]
    assert len(spec["cells"]) == 100
```

- [ ] **Step 2: Run test to verify it fails.**
Run: `PYTHONPATH=~/code/viva-cpm--influenza-sego2022 pytest tests/test_influenza_composite.py -v`
Expected: FAIL (module not found).

- [ ] **Step 3: Implement the composite factory**, mirroring `viva_cpm_studies/composites/gg1993.py` (same `CPM_ADDR`, `composite_generator` usage, `build_spec` for a modest demo). Reuse `sheet.build_sheet_spec` for the geometry; wrap it in a `CPMProcess` config (`mcs_per_update`, `n_fields=0`, `secretory_types=[]`). Keep the live-demo patch small (0.1 mm) so the dashboard "run baseline" is fast; the full-scale run comes from a driver later.
```python
# viva_cpm_studies/composites/influenza.py
"""process-bigraph composite factories for the influenza-sego2022 studies.

Increment 1: the confluent epithelial sheet baseline (CPM only, no fields).
"""
from __future__ import annotations
from process_bigraph.composite_generator import composite_generator
from ..influenza import sheet

CPM_ADDR = "local:!cpm.processes.cpm_process.CPMProcess"


def build_spec(patch_mm: float = 0.1):
    return sheet.build_sheet_spec(patch_mm)


def epithelial_sheet_baseline(patch_mm: float = 0.1) -> dict:
    spec = build_spec(patch_mm)
    # Follow gg1993.composite_document: embed CPMProcess over the world spec.
    # (Copy the exact wrapping gg1993.py uses for load_world + CPMProcess.)
    return composite_generator({
        "cpm": {
            "_type": "process", "address": CPM_ADDR,
            "config": {"mcs_per_update": 10, "n_fields": 0, "secretory_types": []},
            "inputs": {}, "outputs": {},
        },
        "world": spec,
    })
```
Note: match the exact composite wrapping in `gg1993.composite_document` (read it first); the dict above is the shape, adjust keys to whatever `gg1993.py` uses so the dashboard resolves it identically.

- [ ] **Step 4: Run test to verify it passes.**
Run: `PYTHONPATH=~/code/viva-cpm--influenza-sego2022 pytest tests/test_influenza_composite.py -v`
Expected: PASS.

- [ ] **Step 5: Create the `epithelial-sheet-baseline` study** with `/viva-study` (Build/Simulate phase): baseline composite = `viva_cpm_studies.composites.influenza.epithelial_sheet_baseline`; readouts = mean cell volume, cell count, MCS/s (from Task 1.3); acceptance = cells hold ~25 sites & confluent, throughput >= recorded floor. Link to the investigation.

- [ ] **Step 6: Add a visualization** with `/viva-viz`: a snapshot of the sheet (cell types colored) plus the cell-volume distribution histogram, so the increment ships a figure (spec §5).

- [ ] **Step 7: Commit.**
```bash
git add viva_cpm_studies/composites/influenza.py tests/test_influenza_composite.py workspace/
git commit -m "feat(influenza): epithelial-sheet-baseline composite + study + viz"
```

### Task 1.5: Full test run + increment PR

- [ ] **Step 1: Run the whole influenza test subset.**
Run: `PYTHONPATH=~/code/viva-cpm--influenza-sego2022 pytest tests/test_influenza_*.py -v`
Expected: all PASS (params, targets, sheet, composite; perf prints its rate).

- [ ] **Step 2: Run the repo lint** (`scripts/lint-workspace.py`) to confirm the investigation/study YAML is well-formed.

- [ ] **Step 3: Verify branch hygiene** (per repo policy): `git log --oneline origin/main..HEAD` shows only your commits; `git branch --show-current` is `investigation/influenza-sego2022`.

- [ ] **Step 4: Push and open the Increment-1 PR** (title "Influenza (Sego 2022) Increment 1: epithelial sheet at scale"), summarizing Increments 0+1 deliverables and the measured throughput. No AI attribution.

## Self-Review notes

- **Spec coverage:** Increment 0 (spec §5 row 0) → Tasks 0.1–0.5; Increment 1 (row 1) → Tasks 1.1–1.5. Later increments (2–9) are deferred to their own plans per spec §8. The cellularization scaling constants (θ, η) are transcribed in 0.3 but not yet *used* — correct; they're consumed in Increment 8.
- **The `null` adhesion values** in `params.yaml` are a data gap owned by Task 0.1, not a plan placeholder; the sheet builder uses marked CC3D-typical epithelial J and validates by behavior.
- **Type consistency:** `build_sheet_spec(patch_mm, seed)`, `world_from_spec(spec)`, `load_params()`, `load_target(name)`, `epithelial_sheet_baseline()`, `build_spec(patch_mm)` are used identically across tasks; cell-type codes come only from `influenza.types`.
- **Engine changes:** none in Increments 0–1 (uses existing primitives) — the first Rust work (virus field source-gating, then log-chemotaxis/haptotaxis/ρ) starts in Increment 2.
