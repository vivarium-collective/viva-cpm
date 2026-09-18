"""process-bigraph composite factory for the influenza-sego2022 studies.

Increment 1: the confluent epithelial-sheet baseline (CPM only, no fields,
no immune cells). The full multiscale composite (virus field, IFN,
chemokines, macrophages/NK/CD8+ T, global ODE coupling) arrives across
later increments.

Factory returns a composite *document* wrapping the real Rust CPM engine
(``cpm.processes.cpm_process.CPMProcess``) over the epithelial-sheet
geometry built by ``pbg_cpm_studies.influenza.sheet.build_sheet_spec``. The
live-demo patch (0.1 mm, 100 cells) is deliberately small so the dashboard's
"run baseline" is fast; the full-scale 1 mm^2 throughput characterization is
Task 1.3 (``tests/test_influenza_perf.py``), not this live-demo composite.

Referenced from the study's ``baseline[].composite`` as
``pbg_cpm_studies.composites.influenza.epithelial_sheet_baseline``.
"""
from __future__ import annotations

import numpy as np
from process_bigraph.composite_generator import composite_generator

from ..influenza import fields, sheet
from ..influenza import types as inf_types
from ..influenza.params import load_params

CPM_ADDR = "local:!cpm.processes.cpm_process.CPMProcess"
INFECTION_ADDR = "local:!pbg_cpm_studies.influenza.infection_process.InfectionProcess"

# modest live-demo aggregate (full-scale 1mm^2 throughput is characterized
# separately by tests/test_influenza_perf.py, not run live from the dashboard)
DEMO_PATCH_MM = 0.1
DEMO_SEED = 17
DEMO_INIT_INFECTED_FRAC = 0.05


def build_spec(patch_mm: float = DEMO_PATCH_MM) -> dict:
    """Return a load_world spec dict for a modest live demo of the sheet."""
    return sheet.build_sheet_spec(patch_mm)


def composite_document(patch_mm: float = DEMO_PATCH_MM) -> dict:
    """A process-bigraph composite document embedding the CPM engine over the
    confluent epithelial-sheet geometry."""
    spec = build_spec(patch_mm)
    return {
        "cpm": {
            "_type": "process",
            "address": CPM_ADDR,
            "config": {"spec": spec, "mcs_per_update": 10, "n_fields": 0,
                       "secretory_types": []},
            "inputs": {"fates": ["fates"]},
            "outputs": {
                "volumes": ["volumes"],
                "types": ["types"],
                "positions": ["positions"],
                "field_at_cell": ["field_at_cell"],
                "neighbor_secretory": ["neighbor_secretory"],
            },
        },
        "fates": {},
    }


@composite_generator(
    name="epithelial_sheet_baseline", default_n_steps=10,
    description=(
        "Influenza-sego2022 Increment 1: confluent epithelial sheet baseline "
        "(single CPMProcess, no fields, no immune cells). Validates the "
        "substrate geometry — not a biology reproduction."
    ),
    parameters={
        "patch_mm": {"type": "float", "default": DEMO_PATCH_MM,
                     "description": "square patch side length in mm (live demo; small by default)"},
    },
)
def epithelial_sheet_baseline(core=None, patch_mm: float = DEMO_PATCH_MM) -> dict:
    return composite_document(patch_mm)


# ---------------------------------------------------------------------------
# Increment 2 (Task 2.3): virus field + infection transition wired over the
# Increment-1 sheet. Dashboard/live-demo wrapper only -- the quantitative
# lesion-spread measurements (and the integration test) come from
# ``pbg_cpm_studies.influenza.run.run_virus_infection``, which drives the same
# mechanism directly against ``cpm_core.World`` without going through
# process-bigraph. Not a Fig-3B/5/7 reproduction (Increment 9).
# ---------------------------------------------------------------------------

def build_virus_infection_spec(patch_mm: float = DEMO_PATCH_MM, seed: int = DEMO_SEED,
                               init_infected_frac: float = DEMO_INIT_INFECTED_FRAC) -> dict:
    """A ``load_world`` spec for the sheet + virus field, with
    ``init_infected_frac`` of the (seeded-RNG-chosen) H cells pre-set to
    ``types.I`` so the live demo starts from a lesion, not a blank sheet."""
    spec = sheet.build_sheet_spec(patch_mm, seed=seed)
    spec["fields"] = [fields.virus_field_spec_entry()]

    n_cells = len(spec["cells"])
    n_init = int(round(n_cells * init_infected_frac))
    if n_init > 0:
        rng = np.random.default_rng(seed)
        infected_ids = rng.choice(np.arange(1, n_cells + 1), size=n_init, replace=False)
        for cid in infected_ids:
            spec["cells"][int(cid) - 1]["type"] = inf_types.I
    return spec


def virus_infection_composite_document(patch_mm: float = DEMO_PATCH_MM, seed: int = DEMO_SEED,
                                       init_infected_frac: float = DEMO_INIT_INFECTED_FRAC) -> dict:
    """CPMProcess (sheet + virus field) + InfectionProcess (2.2's stochastic
    H -> I transition), wired through the shared ``fates`` map -- matching
    the ``cpm.coupling`` down-scale-coupling convention used by
    ``chemotaxis_receptor.recruitment_receptor``, but with ONE population-
    level infection process (a single shared RNG for the whole sheet) rather
    than one subcell process per cell, since the transition is a population-
    wide stochastic draw, not an independent per-cell deterministic rule."""
    spec = build_virus_infection_spec(patch_mm, seed=seed, init_infected_frac=init_infected_frac)
    g_hv = float(load_params()["virus"]["infection_g_hv"])
    return {
        "cpm": {
            "_type": "process",
            "address": CPM_ADDR,
            "config": {"spec": spec, "mcs_per_update": 10, "n_fields": 1,
                       "secretory_types": [inf_types.I]},
            "inputs": {"fates": ["fates"]},
            "outputs": {
                "volumes": ["volumes"],
                "types": ["types"],
                "positions": ["positions"],
                "field_at_cell": ["field_at_cell"],
                "neighbor_secretory": ["neighbor_secretory"],
            },
        },
        "infection": {
            "_type": "process",
            "address": INFECTION_ADDR,
            "config": {"g_hv": g_hv, "seed": seed},
            "inputs": {"types": ["types"], "field_at_cell": ["field_at_cell"]},
            "outputs": {"fates": ["fates"]},
        },
        "fates": {},
    }


@composite_generator(
    name="virus_infection", default_n_steps=20,
    description=(
        "Influenza-sego2022 Increment 2: virus field (secretion/diffusion/decay "
        "from InfectedReleasing cells) + stochastic H -> I infection transition, "
        "wired over the Increment-1 confluent epithelial sheet. Live-demo/"
        "dashboard wrapper -- the lesion-spread measurements come from "
        "run_virus_infection, not this composite."
    ),
    parameters={
        "patch_mm": {"type": "float", "default": DEMO_PATCH_MM,
                     "description": "square patch side length in mm (live demo; small by default)"},
        "seed": {"type": "int", "default": DEMO_SEED,
                 "description": "RNG seed for both the sheet build and the infection transition"},
        "init_infected_frac": {"type": "float", "default": DEMO_INIT_INFECTED_FRAC,
                               "description": "fraction of H cells seeded as I (InfectedReleasing) at t=0"},
    },
)
def virus_infection(core=None, patch_mm: float = DEMO_PATCH_MM, seed: int = DEMO_SEED,
                    init_infected_frac: float = DEMO_INIT_INFECTED_FRAC) -> dict:
    return virus_infection_composite_document(patch_mm, seed=seed,
                                              init_infected_frac=init_infected_frac)
