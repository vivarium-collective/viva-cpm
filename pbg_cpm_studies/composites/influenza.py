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

from ..influenza import fields, immune, sheet
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


# ---------------------------------------------------------------------------
# Increment 5 (Task 5.2, updated for Task 5.1's fix round 1 -- commit
# 4092d7c): macrophage-response live-demo composite. Wraps the Task-5.1
# non-confluent 2D-approximation scenario (``immune.
# build_macrophage_scenario_spec``) + virus field, with macrophage
# (``inf_types.M``) chemotaxis wired DECLARATIVELY via the field's
# ``chemotaxis`` list (the same ``cpm.schema.load_world`` convention
# ``pbg_cpm_studies.composites.chemotaxis.build_spec`` uses), rather than the
# ``World.set_chemotaxis``/``immune.set_macrophage_chemotaxis`` post-hoc call
# ``run.run_macrophage_response`` makes on a live ``cpm_core.World``.
#
# Geometry (fix round 1): BOTH the epithelial/infection patch and the
# macrophage cluster sit in the domain's INTERIOR (each >= margin_sites from
# every noflux wall, separated by separation_sites of open Medium) -- not a
# domain corner. An earlier version pinned the macrophage cluster in a corner
# diagonally opposite the infection, which biased the lambda=0 control's
# thermal/wall drift to net-point roughly toward the infection (a confound);
# interior placement removes that bias, matching ``run.
# run_macrophage_response``'s current scenario (see immune.py's docstring /
# task-5.1-report.md's fix note).
#
# Dashboard/live-demo wrapper ONLY -- the measured localization numbers cited
# by the macrophage-response study (chemotaxis-on closes the mean distance to
# the infection in 5/5 tried seeds, mean change -13.39 sites, vs a near-
# isotropic lambda=0 control, mean change -0.85 sites) come directly from
# ``run.run_macrophage_response`` across multiple seeds, which ALSO warms the
# virus field up via ``World.advance_fields`` before stepping and records
# ``mean_distance_to_infection``/``macrophage_com`` each update -- a plain
# CPMProcess has neither of those readouts, so this composite reproduces
# neither the warmup nor the localization measurement, only the scenario
# geometry + chemotaxis wiring. Not a Fig-2B/3A reproduction (Increment 9).
# ---------------------------------------------------------------------------

DEMO_MACROPHAGE_CELLS_PER_SIDE = 4
DEMO_N_INFECTED = 1
DEMO_N_MACROPHAGES = 6
DEMO_MARGIN_SITES = 30
DEMO_SEPARATION_SITES = 25
DEMO_CHEMOTAXIS_V_MACRO = 5000.0  # params.yaml macrophage.chemotaxis_v_macro; pass 0.0 for the control


def build_macrophage_response_spec(*, epithelial_cells_per_side: int = DEMO_MACROPHAGE_CELLS_PER_SIDE,
                                   n_infected: int = DEMO_N_INFECTED,
                                   n_macrophages: int = DEMO_N_MACROPHAGES,
                                   margin_sites: int = DEMO_MARGIN_SITES,
                                   separation_sites: int = DEMO_SEPARATION_SITES,
                                   seed: int = DEMO_SEED,
                                   chemotaxis_v_macro: float = DEMO_CHEMOTAXIS_V_MACRO) -> dict:
    """A ``load_world`` spec for the Task-5.1 macrophage scenario + virus
    field -- both the epithelial/infection patch and the macrophage cluster
    placed in the domain's INTERIOR (unbiased lambda=0-control geometry, see
    module note above) -- with macrophage (``inf_types.M``) chemotaxis on
    the virus field wired via the field's ``chemotaxis`` list.
    ``chemotaxis_v_macro=0.0`` reproduces the lambda=0 control (no directed
    chemotaxis)."""
    spec = immune.build_macrophage_scenario_spec(
        epithelial_cells_per_side=epithelial_cells_per_side, n_infected=n_infected,
        n_macrophages=n_macrophages, margin_sites=margin_sites,
        separation_sites=separation_sites, seed=seed)
    field_entry = fields.virus_field_spec_entry()
    field_entry["chemotaxis"] = [{"type": inf_types.M, "lambda": float(chemotaxis_v_macro)}]
    spec["fields"] = [field_entry]
    return spec


def macrophage_response_composite_document(*, epithelial_cells_per_side=DEMO_MACROPHAGE_CELLS_PER_SIDE,
                                            n_infected=DEMO_N_INFECTED, n_macrophages=DEMO_N_MACROPHAGES,
                                            margin_sites=DEMO_MARGIN_SITES,
                                            separation_sites=DEMO_SEPARATION_SITES, seed=DEMO_SEED,
                                            chemotaxis_v_macro=DEMO_CHEMOTAXIS_V_MACRO) -> dict:
    spec = build_macrophage_response_spec(
        epithelial_cells_per_side=epithelial_cells_per_side, n_infected=n_infected,
        n_macrophages=n_macrophages, margin_sites=margin_sites,
        separation_sites=separation_sites, seed=seed,
        chemotaxis_v_macro=chemotaxis_v_macro)
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
        "fates": {},
    }


@composite_generator(
    name="macrophage_response", default_n_steps=40,
    description=(
        "Influenza-sego2022 Increment 5: macrophage (type M) chemotaxis up the "
        "virus field toward the infection, over the Task-5.1 non-confluent "
        "2D-approximation scenario (both clusters placed in the domain interior, "
        "an unbiased lambda=0-control geometry). Live-demo/dashboard wrapper -- "
        "the measured localization numbers (chemotaxis-on closes the distance in "
        "5/5 tried seeds, mean -13.39 sites, vs a near-isotropic lambda=0 control, "
        "mean -0.85 sites) come from run.run_macrophage_response across multiple "
        "seeds, not this composite. chemotaxis_v_macro=0.0 reproduces the "
        "lambda=0 control."
    ),
    parameters={
        "chemotaxis_v_macro": {"type": "float", "default": DEMO_CHEMOTAXIS_V_MACRO,
                               "description": "macrophage chemotaxis strength on the virus field (0.0 = lambda=0 control)"},
        "seed": {"type": "int", "default": DEMO_SEED,
                 "description": "RNG seed for the scenario build"},
    },
)
def macrophage_response(core=None, chemotaxis_v_macro: float = DEMO_CHEMOTAXIS_V_MACRO,
                        seed: int = DEMO_SEED) -> dict:
    return macrophage_response_composite_document(seed=seed, chemotaxis_v_macro=chemotaxis_v_macro)


# ---------------------------------------------------------------------------
# Increment 7 (Task 7.3): NK/CD8 cytotoxic-response live-demo composite.
# Wraps the Task-7.1 three-cluster scenario (``immune.
# build_cytotoxic_scenario_spec`` -- epithelial/infection patch | macrophage
# cluster | NK+CD8 cluster, every cluster interior-placed) + virus field
# (macrophage chemotaxis, unchanged from Increment 5) + chemokine field
# (``fields.chemokine_field_spec_entry``, Task 7.3-new -- macrophage-secreted,
# Increment 6) with NK (``inf_types.K``)/CD8+ (``inf_types.E``) chemotaxis
# wired DECLARATIVELY via the chemokine field's ``chemotaxis`` list, the same
# ``cpm.schema.load_world`` convention ``macrophage_response`` above uses.
#
# Dashboard/live-demo wrapper ONLY, same limitation as ``macrophage_response``
# above (and now compounded): a plain declarative CPMProcess has no field
# warmup, no per-cell IL-10-Hill secretion regulation, no NK/CD8 CONTACT-
# KILLING step (``killing.contact_kill_rate`` is a custom per-update Python
# calculation over ``world.cell_contact_area_by_type``, not expressible in
# this declarative spec format), and none of ``run.run_cytotoxic_response``'s
# distance/n_infected readouts. This composite reproduces ONLY the scenario
# geometry + chemotaxis wiring (localization half); it does NOT reproduce or
# claim the killing mechanism, and does NOT reproduce the measured numbers
# cited by the cytotoxic-killing study (those come directly from
# ``run.run_cytotoxic_response``, a raw ``cpm_core.World`` loop, across
# multiple seeds/conditions -- see that study's ``model_change.notes``).
# ``chemotaxis_v_nk``/``chemotaxis_v_cd8`` default to Task 7.1's chosen
# ``run.NK_CD8_CHEMOTAXIS_ENGINE_SCALE=100x`` scale (params.yaml's literal
# 5000/10000 would be statistically undetectable at this field's
# concentration scale, see that constant's docstring); pass 0.0 for the
# lambda=0 localization control.
# ---------------------------------------------------------------------------

DEMO_N_NK = 6
DEMO_N_CD8 = 6
DEMO_CYTOTOXIC_MARGIN_SITES = 20
DEMO_CYTOTOXIC_SEPARATION_SITES = 20
DEMO_CHEMOTAXIS_V_NK = 500000.0    # params.yaml nk.chemotaxis_v_nk (5000) x 100 (run.NK_CD8_CHEMOTAXIS_ENGINE_SCALE)
DEMO_CHEMOTAXIS_V_CD8 = 1000000.0  # params.yaml cd8.chemotaxis_v_cd8 (10000) x 100 (same scale, preserves 2:1)


def build_cytotoxic_response_spec(*, epithelial_cells_per_side: int = DEMO_MACROPHAGE_CELLS_PER_SIDE,
                                  n_infected: int = DEMO_N_INFECTED,
                                  n_macrophages: int = DEMO_N_MACROPHAGES,
                                  n_nk: int = DEMO_N_NK, n_cd8: int = DEMO_N_CD8,
                                  margin_sites: int = DEMO_CYTOTOXIC_MARGIN_SITES,
                                  separation_sites: int = DEMO_CYTOTOXIC_SEPARATION_SITES,
                                  seed: int = DEMO_SEED,
                                  chemotaxis_v_macro: float = DEMO_CHEMOTAXIS_V_MACRO,
                                  chemotaxis_v_nk: float = DEMO_CHEMOTAXIS_V_NK,
                                  chemotaxis_v_cd8: float = DEMO_CHEMOTAXIS_V_CD8) -> dict:
    """A ``load_world`` spec for the Task-7.1 three-cluster cytotoxic
    scenario + virus field (macrophage chemotaxis) + chemokine field (NK/CD8
    chemotaxis) -- every cluster placed in the domain's INTERIOR (see
    ``immune.build_cytotoxic_scenario_spec``'s docstring). ``chemotaxis_v_nk``/
    ``chemotaxis_v_cd8=0.0`` reproduces the localization lambda=0 control (no
    directed NK/CD8 chemotaxis); NO killing step exists in this declarative
    composite (see module note above)."""
    spec = immune.build_cytotoxic_scenario_spec(
        epithelial_cells_per_side=epithelial_cells_per_side, n_infected=n_infected,
        n_macrophages=n_macrophages, n_nk=n_nk, n_cd8=n_cd8,
        margin_sites=margin_sites, separation_sites=separation_sites, seed=seed)
    virus_entry = fields.virus_field_spec_entry()
    virus_entry["chemotaxis"] = [{"type": inf_types.M, "lambda": float(chemotaxis_v_macro)}]
    chemo_entry = fields.chemokine_field_spec_entry()
    chemo_entry["chemotaxis"] = [
        {"type": inf_types.K, "lambda": float(chemotaxis_v_nk)},
        {"type": inf_types.E, "lambda": float(chemotaxis_v_cd8)},
    ]
    spec["fields"] = [virus_entry, chemo_entry]
    return spec


def cytotoxic_response_composite_document(*, epithelial_cells_per_side=DEMO_MACROPHAGE_CELLS_PER_SIDE,
                                          n_infected=DEMO_N_INFECTED, n_macrophages=DEMO_N_MACROPHAGES,
                                          n_nk=DEMO_N_NK, n_cd8=DEMO_N_CD8,
                                          margin_sites=DEMO_CYTOTOXIC_MARGIN_SITES,
                                          separation_sites=DEMO_CYTOTOXIC_SEPARATION_SITES, seed=DEMO_SEED,
                                          chemotaxis_v_macro=DEMO_CHEMOTAXIS_V_MACRO,
                                          chemotaxis_v_nk=DEMO_CHEMOTAXIS_V_NK,
                                          chemotaxis_v_cd8=DEMO_CHEMOTAXIS_V_CD8) -> dict:
    spec = build_cytotoxic_response_spec(
        epithelial_cells_per_side=epithelial_cells_per_side, n_infected=n_infected,
        n_macrophages=n_macrophages, n_nk=n_nk, n_cd8=n_cd8,
        margin_sites=margin_sites, separation_sites=separation_sites, seed=seed,
        chemotaxis_v_macro=chemotaxis_v_macro, chemotaxis_v_nk=chemotaxis_v_nk,
        chemotaxis_v_cd8=chemotaxis_v_cd8)
    return {
        "cpm": {
            "_type": "process",
            "address": CPM_ADDR,
            "config": {"spec": spec, "mcs_per_update": 10, "n_fields": 2,
                       "secretory_types": [inf_types.I, inf_types.M]},
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
    name="cytotoxic_response", default_n_steps=60,
    description=(
        "Influenza-sego2022 Increment 7: NK (type K) + CD8+ (type E) chemotaxis "
        "up the macrophage-released chemokine field, over the Task-7.1 "
        "three-cluster non-confluent scenario (epithelial/infection patch | "
        "macrophage cluster | NK+CD8 cluster, every cluster interior-placed). "
        "Live-demo/dashboard wrapper -- LOCALIZATION geometry/chemotaxis only, "
        "NO contact-killing step (not expressible in this declarative composite "
        "format; see killing.contact_kill_rate / run.run_cytotoxic_response for "
        "the actual killing mechanism and its measured numbers). "
        "chemotaxis_v_nk=chemotaxis_v_cd8=0.0 reproduces the localization "
        "lambda=0 control."
    ),
    parameters={
        "chemotaxis_v_macro": {"type": "float", "default": DEMO_CHEMOTAXIS_V_MACRO,
                               "description": "macrophage chemotaxis strength on the virus field"},
        "chemotaxis_v_nk": {"type": "float", "default": DEMO_CHEMOTAXIS_V_NK,
                            "description": "NK chemotaxis strength on the chemokine field (0.0 = lambda=0 control)"},
        "chemotaxis_v_cd8": {"type": "float", "default": DEMO_CHEMOTAXIS_V_CD8,
                             "description": "CD8+ chemotaxis strength on the chemokine field (0.0 = lambda=0 control)"},
        "seed": {"type": "int", "default": DEMO_SEED,
                 "description": "RNG seed for the scenario build"},
    },
)
def cytotoxic_response(core=None, chemotaxis_v_macro: float = DEMO_CHEMOTAXIS_V_MACRO,
                       chemotaxis_v_nk: float = DEMO_CHEMOTAXIS_V_NK,
                       chemotaxis_v_cd8: float = DEMO_CHEMOTAXIS_V_CD8,
                       seed: int = DEMO_SEED) -> dict:
    return cytotoxic_response_composite_document(
        seed=seed, chemotaxis_v_macro=chemotaxis_v_macro,
        chemotaxis_v_nk=chemotaxis_v_nk, chemotaxis_v_cd8=chemotaxis_v_cd8)


# ---------------------------------------------------------------------------
# Increment 8 (Task 8.6): global-coupling live-demo composite. Same
# limitation as ``cytotoxic_response`` above, now compounded further: a
# plain declarative ``CPMProcess`` has no per-MCS `price_ode.GlobalODE`
# integration step (Task 8.1), no spatial->ODE aggregate push / ODE->spatial
# dynamic sig_1 feedback (Tasks 8.2/8.3), no ODE-driven recruitment
# inflow/outflow (Task 8.4), and no NK/CD8 LOCAL+NEARBY cytotoxic killing
# (Task 8.5) -- none of these are expressible in this declarative spec
# format (custom per-MCS Python over ``price_ode.GlobalODE``/`recruitment`/
# `killing`, not a CPMProcess config). This composite reproduces ONLY the
# Task-7.1 three-cluster scenario geometry + the virus/chemokine field
# stack + NK/CD8/macrophage chemotaxis wiring (identical to
# ``cytotoxic_response`` above); it does NOT reproduce or claim the global
# ODE coupling, dynamic sig_1, recruitment, or nearby-killing mechanisms --
# those come directly from ``run.run_global_coupling`` (Tasks 8.1-8.5, a raw
# ``cpm_core.World`` + ``price_ode.GlobalODE`` loop), not from this
# composite. See that study's ``model_change.notes`` and
# ``pbg_cpm_studies/influenza/viz.py::global_coupling_figure``, which render
# ``run.run_global_coupling``'s actual output, not this composite's.
# ---------------------------------------------------------------------------


def global_coupling_composite_document(*, epithelial_cells_per_side=DEMO_MACROPHAGE_CELLS_PER_SIDE,
                                        n_infected=DEMO_N_INFECTED, n_macrophages=DEMO_N_MACROPHAGES,
                                        n_nk=DEMO_N_NK, n_cd8=DEMO_N_CD8,
                                        margin_sites=DEMO_CYTOTOXIC_MARGIN_SITES,
                                        separation_sites=DEMO_CYTOTOXIC_SEPARATION_SITES, seed=DEMO_SEED,
                                        chemotaxis_v_macro=DEMO_CHEMOTAXIS_V_MACRO,
                                        chemotaxis_v_nk=DEMO_CHEMOTAXIS_V_NK,
                                        chemotaxis_v_cd8=DEMO_CHEMOTAXIS_V_CD8) -> dict:
    """Reuses ``cytotoxic_response_composite_document`` verbatim -- see the
    module note above for why the global-ODE coupling/recruitment/killing
    mechanisms (Tasks 8.1-8.5) are not expressible in this declarative
    composite format."""
    return cytotoxic_response_composite_document(
        epithelial_cells_per_side=epithelial_cells_per_side, n_infected=n_infected,
        n_macrophages=n_macrophages, n_nk=n_nk, n_cd8=n_cd8,
        margin_sites=margin_sites, separation_sites=separation_sites, seed=seed,
        chemotaxis_v_macro=chemotaxis_v_macro, chemotaxis_v_nk=chemotaxis_v_nk,
        chemotaxis_v_cd8=chemotaxis_v_cd8)


@composite_generator(
    name="global_coupling", default_n_steps=60,
    description=(
        "Influenza-sego2022 Increment 8: live-demo/dashboard wrapper for the "
        "Task-7.1 three-cluster scenario (epithelial/infection patch | "
        "macrophage cluster | NK+CD8 cluster) + virus/chemokine field stack "
        "+ macrophage/NK/CD8 chemotaxis wiring -- SCENARIO GEOMETRY ONLY. "
        "Does NOT implement the Task 8.1-8.5 hybrid Price-2015 global ODE "
        "(10 systemic species), the spatial<->ODE bidirectional coupling, "
        "dynamic sig_1 secretion feedback, ODE-driven recruitment, or the "
        "NK/CD8 nearby-population killing term -- none are expressible in "
        "this declarative composite format (custom per-MCS Python over "
        "price_ode.GlobalODE/recruitment/killing). See "
        "run.run_global_coupling and viz.global_coupling_figure for the "
        "actual coupled mechanism and its measured (calibration-pending) "
        "numbers."
    ),
    parameters={
        "chemotaxis_v_macro": {"type": "float", "default": DEMO_CHEMOTAXIS_V_MACRO,
                               "description": "macrophage chemotaxis strength on the virus field"},
        "chemotaxis_v_nk": {"type": "float", "default": DEMO_CHEMOTAXIS_V_NK,
                            "description": "NK chemotaxis strength on the chemokine field (0.0 = lambda=0 control)"},
        "chemotaxis_v_cd8": {"type": "float", "default": DEMO_CHEMOTAXIS_V_CD8,
                             "description": "CD8+ chemotaxis strength on the chemokine field (0.0 = lambda=0 control)"},
        "seed": {"type": "int", "default": DEMO_SEED,
                 "description": "RNG seed for the scenario build"},
    },
)
def global_coupling(core=None, chemotaxis_v_macro: float = DEMO_CHEMOTAXIS_V_MACRO,
                    chemotaxis_v_nk: float = DEMO_CHEMOTAXIS_V_NK,
                    chemotaxis_v_cd8: float = DEMO_CHEMOTAXIS_V_CD8,
                    seed: int = DEMO_SEED) -> dict:
    return global_coupling_composite_document(
        seed=seed, chemotaxis_v_macro=chemotaxis_v_macro,
        chemotaxis_v_nk=chemotaxis_v_nk, chemotaxis_v_cd8=chemotaxis_v_cd8)
