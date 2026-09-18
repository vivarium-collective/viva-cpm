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

from process_bigraph.composite_generator import composite_generator

from ..influenza import sheet

CPM_ADDR = "local:!cpm.processes.cpm_process.CPMProcess"

# modest live-demo aggregate (full-scale 1mm^2 throughput is characterized
# separately by tests/test_influenza_perf.py, not run live from the dashboard)
DEMO_PATCH_MM = 0.1


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
